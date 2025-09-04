#!/usr/bin/env python3

from Bio import SeqIO
import numpy as np
import pandas as pd
from itertools import combinations
import re
from collections import defaultdict
import os
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

class PromoterSimilarityAnalyzer:
    def __init__(self):
        """Complete promoter similarity analysis with TFBS patterns and normalization"""
        
        # Plant-specific TFBS patterns from PlantCARE, PLACE, JASPAR databases
        self.plant_tfbs_patterns = {
            # Core promoter elements
            'TATA_box': [r'TATAAA', r'TATAWAW', r'TATAWAR'],
            'CAAT_box': [r'CAAT', r'CCAAT', r'CAAAT'],
            'Initiator': [r'[CT][CT]A[ATGC][AT][CT][CT]'],
            
            # Major plant TF families
            'MYB_binding': [r'[CT]AACNG', r'YAACKG', r'CAACAG'],
            'bZIP_binding': [r'[TG]GACGT[CA]', r'ACGT', r'CACGTG'],
            'WRKY_binding': [r'[CT]TGAC[CT]', r'TTGACY'],
            'AP2_ERF': [r'GCCGCC', r'AGCCGCC'],
            'bHLH_binding': [r'CANNTG', r'CACGTG', r'CATGTG'],
            'NAC_binding': [r'CACG', r'CATGTG'],
            'GATA_binding': [r'[AT]GATA[AG]', r'HGATAR'],
            'DOF_binding': [r'AAAGG', r'AAAG'],
            'TCP_binding': [r'GGNCCCAC'],
            'HSF_binding': [r'GAANNTTC', r'GAANNTTCNNGAANNTTC'],
            
            # Hormone response elements
            'ABA_response': [r'ACGTGG[CT]', r'ACGTGTC'],
            'Auxin_response': [r'TGTCTC', r'GAGACA'],
            'GA_response': [r'TAACAAA', r'TATCCAC'],
            'Cytokinin_response': [r'TATTAG', r'AGATCCT'],
            'Ethylene_response': [r'ATTTCAAA', r'TAAGAGCCGCC'],
            
            # Stress response elements
            'Drought_response': [r'TACCGACAT', r'ACGTGG[CT]'],
            'Cold_response': [r'TGGCCGAC', r'CCGAC'],
            'Heat_response': [r'GAANNTTC'],
            'Light_response': [r'ACGT', r'CAATCA', r'ATCTTA'],
            
            # Tissue-specific elements
            'Root_specific': [r'ATATT', r'TTATT'],
            'Leaf_specific': [r'CAATCA', r'ATCTTA'],
            'Seed_specific': [r'ACAAAA', r'CATGCA'],
            
            # Other regulatory elements
            'CpG_island': [r'CG'],
            'GC_box': [r'GGGCGG', r'CCGCCC'],
            'E_box': [r'CANNTG'],
            'CCAAT_box': [r'CCAAT'],
            'Silencer': [r'AATTT'],
            'Enhancer': [r'TGACGT[CA]']
        }
    
    def scan_tfbs_patterns(self, sequence):
        """Scan sequence for known TFBS patterns"""
        sequence = str(sequence).upper()
        tfbs_matches = defaultdict(list)
        
        for tfbs_name, patterns in self.plant_tfbs_patterns.items():
            for pattern in patterns:
                # Convert IUPAC nucleotide codes to regex
                regex_pattern = pattern
                iupac_codes = {
                    'W': '[AT]', 'S': '[GC]', 'M': '[AC]', 'K': '[GT]',
                    'R': '[AG]', 'Y': '[CT]', 'B': '[CGT]', 'D': '[AGT]',
                    'H': '[ACT]', 'V': '[ACG]', 'N': '[ATGC]'
                }
                
                for code, replacement in iupac_codes.items():
                    regex_pattern = regex_pattern.replace(code, replacement)
                
                # Find all matches
                matches = list(re.finditer(regex_pattern, sequence))
                for match in matches:
                    tfbs_matches[tfbs_name].append({
                        'start': match.start(),
                        'end': match.end(),
                        'sequence': match.group()
                    })
        
        return dict(tfbs_matches)
    
    def analyze_promoter_regions(self, sequences):
        """Analyze proximal and distal promoter regions separately"""
        regional_analysis = {}
        
        for i, seq in enumerate(sequences):
            seq_str = str(seq.seq).upper()
            length = len(seq_str)
            
            # Split promoter into regions (assuming TSS at end)
            if length >= 800:
                # Standard split: distal (-800 to -200), proximal (-200 to TSS)
                distal_region = seq_str[:600]    # First 600bp = distal
                proximal_region = seq_str[600:]  # Last 400bp = proximal
            elif length >= 400:
                # Shorter promoters: split in half
                split_point = length // 2
                distal_region = seq_str[:split_point]
                proximal_region = seq_str[split_point:]
            else:
                # Very short: use entire sequence for both
                distal_region = seq_str
                proximal_region = seq_str
            
            # Scan for TFBS in each region
            distal_tfbs = self.scan_tfbs_patterns(distal_region)
            proximal_tfbs = self.scan_tfbs_patterns(proximal_region)
            
            regional_analysis[i] = {
                'seq_id': seq.id,
                'total_length': length,
                'distal_region': {
                    'length': len(distal_region),
                    'tfbs': distal_tfbs,
                    'tfbs_count': sum(len(matches) for matches in distal_tfbs.values())
                },
                'proximal_region': {
                    'length': len(proximal_region),
                    'tfbs': proximal_tfbs,
                    'tfbs_count': sum(len(matches) for matches in proximal_tfbs.values())
                }
            }
        
        return regional_analysis
    
    def calculate_tfbs_similarity(self, tfbs1, tfbs2):
        """Calculate similarity between two TFBS profiles"""
        
        # Get all TFBS types present in either promoter
        all_types = set(tfbs1.keys()) | set(tfbs2.keys())
        if not all_types:
            return 0.0
        
        # Jaccard similarity: shared TFBS types / total TFBS types
        shared_types = set(tfbs1.keys()) & set(tfbs2.keys())
        jaccard_similarity = len(shared_types) / len(all_types)
        
        # Correlation of TFBS densities
        counts1 = [len(tfbs1.get(tfbs_type, [])) for tfbs_type in all_types]
        counts2 = [len(tfbs2.get(tfbs_type, [])) for tfbs_type in all_types]
        
        if sum(counts1) > 0 and sum(counts2) > 0:
            # Normalize counts to get densities
            total1, total2 = sum(counts1), sum(counts2)
            norm_counts1 = [c/total1 for c in counts1]
            norm_counts2 = [c/total2 for c in counts2]
            
            # Calculate Pearson correlation
            if np.std(norm_counts1) > 0 and np.std(norm_counts2) > 0:
                correlation = np.corrcoef(norm_counts1, norm_counts2)[0, 1]
                correlation = max(0, correlation)  # Only positive correlations
            else:
                correlation = 0
        else:
            correlation = 0
        
        # Combined similarity: 60% shared types + 40% density correlation
        return 0.6 * jaccard_similarity + 0.4 * correlation
    
    def analyze_group_similarity(self, fasta_file, group_name):
        """Main analysis function with normalization"""
        try:
            # Read sequences
            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            
            if len(sequences) < 2:
                return {
                    'group_name': group_name,
                    'num_promoters': len(sequences),
                    'num_comparisons': 0,
                    'similarity_score': 0.0,
                    'normalized_similarity': 0.0,
                    'similarity_category': 'Insufficient sequences',
                    'status': 'error',
                    'error_msg': f'Need >=2 sequences, found {len(sequences)}'
                }
            
            # Analyze promoter regions
            regional_data = self.analyze_promoter_regions(sequences)
            
            # Calculate all pairwise similarities
            similarities = []
            for i, j in combinations(range(len(sequences)), 2):
                seq1_data = regional_data[i]
                seq2_data = regional_data[j]
                
                # Calculate regional similarities
                proximal_sim = self.calculate_tfbs_similarity(
                    seq1_data['proximal_region']['tfbs'],
                    seq2_data['proximal_region']['tfbs']
                )
                
                distal_sim = self.calculate_tfbs_similarity(
                    seq1_data['distal_region']['tfbs'],
                    seq2_data['distal_region']['tfbs']
                )
                
                # Weighted overall similarity (proximal region more important)
                overall_sim = 0.7 * proximal_sim + 0.3 * distal_sim
                
                similarities.append({
                    'proximal_similarity': proximal_sim,
                    'distal_similarity': distal_sim,
                    'overall_similarity': overall_sim
                })
            
            if not similarities:
                return {
                    'group_name': group_name,
                    'num_promoters': len(sequences),
                    'similarity_score': 0.0,
                    'normalized_similarity': 0.0,
                    'similarity_category': 'No similarities calculated',
                    'status': 'error',
                    'error_msg': 'Failed to calculate pairwise similarities'
                }
            
            # Calculate statistics
            df = pd.DataFrame(similarities)
            raw_similarity = df['overall_similarity'].mean() * 100
            similarity_std = df['overall_similarity'].std() * 100
            
            # NORMALIZATION based on group characteristics
            num_genes = len(sequences)
            num_comparisons = len(similarities)
            
            # 1. Size adjustment (optimal size ~5 genes)
            optimal_size = 5
            size_penalty = 1.0 - abs(num_genes - optimal_size) / max(num_genes, optimal_size) * 0.2
            size_adjusted = raw_similarity * max(0.6, size_penalty)
            
            # 2. Statistical robustness (more comparisons = more reliable)
            if num_comparisons >= 10:
                comparison_weight = 1.0
            elif num_comparisons >= 6:
                comparison_weight = 0.95
            elif num_comparisons >= 3:
                comparison_weight = 0.85
            else:
                comparison_weight = 0.7
            
            comparison_weighted = raw_similarity * comparison_weight
            
            # 3. Variance penalty (high variance = less consistent)
            if similarity_std > 20:
                std_penalty = 0.9
            elif similarity_std > 15:
                std_penalty = 0.95
            else:
                std_penalty = 1.0
            
            std_adjusted = raw_similarity * std_penalty
            
            # Combined normalized score
            normalized_similarity = (
                size_adjusted * 0.3 +           # 30% size adjustment
                comparison_weighted * 0.5 +     # 50% comparison reliability  
                std_adjusted * 0.2              # 20% variance penalty
            )
            
            # Get TFBS statistics
            all_tfbs_types = set()
            for data in regional_data.values():
                all_tfbs_types.update(data['proximal_region']['tfbs'].keys())
                all_tfbs_types.update(data['distal_region']['tfbs'].keys())
            
            # Categorize similarity (using normalized score)
            if normalized_similarity >= 70:
                category = "Very High Similarity"
            elif normalized_similarity >= 50:
                category = "High Similarity"
            elif normalized_similarity >= 30:
                category = "Moderate Similarity"
            elif normalized_similarity >= 15:
                category = "Low Similarity"
            else:
                category = "Very Low Similarity"
            
            return {
                'group_name': group_name,
                'num_promoters': num_genes,
                'num_comparisons': num_comparisons,
                'similarity_score': round(raw_similarity, 2),
                'normalized_similarity': round(normalized_similarity, 2),
                'similarity_std': round(similarity_std, 2),
                'size_adjustment_factor': round(size_penalty, 3),
                'comparison_weight': round(comparison_weight, 3),
                'std_adjustment_factor': round(std_penalty, 3),
                'mean_proximal_similarity': round(df['proximal_similarity'].mean(), 3),
                'mean_distal_similarity': round(df['distal_similarity'].mean(), 3),
                'num_tfbs_types_found': len(all_tfbs_types),
                'similarity_category': category,
                'status': 'success'
            }
            
        except Exception as e:
            return {
                'group_name': group_name,
                'num_promoters': 0,
                'similarity_score': 0.0,
                'normalized_similarity': 0.0,
                'similarity_category': 'Analysis Error',
                'status': 'error',
                'error_msg': str(e)
            }

def process_single_file(args):
    """Process single file for concurrent execution"""
    fasta_file, group_name = args
    
    # Basic file validation
    if not os.path.exists(fasta_file):
        return {
            'group_name': group_name,
            'similarity_score': 0.0,
            'normalized_similarity': 0.0,
            'status': 'error',
            'error_msg': 'File not found'
        }
    
    if os.path.getsize(fasta_file) == 0:
        return {
            'group_name': group_name,
            'similarity_score': 0.0,
            'normalized_similarity': 0.0,
            'status': 'error',
            'error_msg': 'Empty file'
        }
    
    # Run analysis
    analyzer = PromoterSimilarityAnalyzer()
    return analyzer.analyze_group_similarity(fasta_file, group_name)

def main():
    """Main function with concurrent processing and group analysis"""
    
    # Configuration
    promoter_directory = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_promotor_analysis/"
    
    # Find all promoter files
    fasta_files = glob.glob(os.path.join(promoter_directory, "*_promoters.fasta"))
    
    if not fasta_files:
        print("No promoter files found!")
        print(f"Searched in: {promoter_directory}")
        print("Pattern: *_promoters.fasta")
        
        # Try alternative patterns
        alternatives = ["*.fasta", "*.fa", "*promoter*.fasta"]
        for pattern in alternatives:
            alt_files = glob.glob(os.path.join(promoter_directory, pattern))
            if alt_files:
                print(f"\nFound {len(alt_files)} files with pattern {pattern}")
                response = input(f"Process these files? (y/n): ")
                if response.lower().startswith('y'):
                    fasta_files = alt_files
                break
        
        if not fasta_files:
            return
    
    print(f"Found {len(fasta_files)} promoter files to analyze")
    
    # Extract group names and classify by dataset
    group_names = [os.path.basename(f).replace('_promoters.fasta', '') for f in fasta_files]
    
    kegg_files = [(f, g) for f, g in zip(fasta_files, group_names) if g.startswith('MGC_CANDIDATE')]
    mibig_files = [(f, g) for f, g in zip(fasta_files, group_names) if g.startswith('BGC')]
    random_files = [(f, g) for f, g in zip(fasta_files, group_names) if g.startswith('RANDOM')]
    
    print(f"\nDataset breakdown:")
    print(f"  KEGG (MGC_CANDIDATE): {len(kegg_files)} files")
    print(f"  MiBIG (BGC): {len(mibig_files)} files")
    print(f"  Random (RANDOM): {len(random_files)} files")
    
    # Concurrent processing setup
    MAX_WORKERS = min(len(fasta_files), 6)
    args_list = list(zip(fasta_files, group_names))
    
    print(f"\nStarting concurrent analysis with {MAX_WORKERS} workers...")
    
    results = []
    completed = 0
    start_time = time.time()
    
    # Progress tracking by group
    group_counts = {'KEGG': 0, 'MiBIG': 0, 'Random': 0}
    
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_args = {executor.submit(process_single_file, args): args for args in args_list}
        
        for future in as_completed(future_to_args):
            args = future_to_args[future]
            group_name = args[1]
            
            try:
                result = future.result()
                results.append(result)
                completed += 1
                
                # Classify group type
                if group_name.startswith('MGC_CANDIDATE'):
                    group_counts['KEGG'] += 1
                    group_type = "KEGG"
                elif group_name.startswith('BGC'):
                    group_counts['MiBIG'] += 1
                    group_type = "MiBIG"
                elif group_name.startswith('RANDOM'):
                    group_counts['Random'] += 1
                    group_type = "Random"
                else:
                    group_type = "Other"
                
                # Progress calculation
                elapsed = time.time() - start_time
                rate = completed / elapsed if elapsed > 0 else 0
                eta = (len(fasta_files) - completed) / rate if rate > 0 else 0
                
                # Display progress
                if result.get('status') == 'success':
                    raw_score = result.get('similarity_score', 0)
                    norm_score = result.get('normalized_similarity', 0)
                    num_seqs = result.get('num_promoters', 0)
                    print(f"✅ [{completed:4d}/{len(fasta_files)}] {group_type:6s} {group_name}: Raw:{raw_score:.1f} Norm:{norm_score:.1f} ({num_seqs} seqs, {rate:.1f}/sec, ETA: {eta/60:.0f}m)")
                else:
                    error_msg = result.get('error_msg', 'Unknown error')[:40]
                    print(f"❌ [{completed:4d}/{len(fasta_files)}] {group_type:6s} {group_name}: {error_msg}")
                
            except Exception as e:
                print(f"❌ [{completed:4d}/{len(fasta_files)}] {group_name}: Exception - {str(e)[:40]}")
                completed += 1
    
    total_time = time.time() - start_time
    
    # Process and save results
    if results:
        df = pd.DataFrame(results)
        
        # Add dataset group classification
        def classify_group(name):
            if name.startswith('MGC_CANDIDATE'):
                return 'KEGG'
            elif name.startswith('BGC'):
                return 'MiBIG'
            elif name.startswith('RANDOM'):
                return 'Random'
            else:
                return 'Other'
        
        df['dataset_group'] = df['group_name'].apply(classify_group)
        successful = df[df['status'] == 'success']
        
        print(f"\n{'='*80}")
        print("ANALYSIS COMPLETE!")
        print(f"{'='*80}")
        print(f"Total time: {total_time/60:.1f} minutes")
        print(f"Processing rate: {len(results)/(total_time/60):.1f} files/minute")
        print(f"Successful analyses: {len(successful)}/{len(results)} ({len(successful)/len(results)*100:.1f}%)")
        
        if len(successful) > 0:
            # Group-wise statistics
            print(f"\nRESULTS BY DATASET:")
            print(f"{'Group':<8} {'Count':<6} {'Success':<8} {'Raw Mean':<10} {'Norm Mean':<10} {'Raw Std':<8} {'Norm Std':<8}")
            print("-" * 70)
            
            for group in ['KEGG', 'MiBIG', 'Random']:
                group_data = df[df['dataset_group'] == group]
                group_success = group_data[group_data['status'] == 'success']
                
                if len(group_success) > 0:
                    raw_mean = group_success['similarity_score'].mean()
                    norm_mean = group_success['normalized_similarity'].mean()
                    raw_std = group_success['similarity_score'].std()
                    norm_std = group_success['normalized_similarity'].std()
                    
                    print(f"{group:<8} {len(group_data):<6} {len(group_success):<8} {raw_mean:<10.1f} {norm_mean:<10.1f} {raw_std:<8.1f} {norm_std:<8.1f}")
                else:
                    print(f"{group:<8} {len(group_data):<6} {'0':<8} {'N/A':<10} {'N/A':<10} {'N/A':<8} {'N/A':<8}")
            
            # Statistical comparisons
            print(f"\nSTATISTICAL COMPARISONS:")
            
            kegg_raw = successful[successful['dataset_group'] == 'KEGG']['similarity_score']
            kegg_norm = successful[successful['dataset_group'] == 'KEGG']['normalized_similarity']
            mibig_raw = successful[successful['dataset_group'] == 'MiBIG']['similarity_score']
            mibig_norm = successful[successful['dataset_group'] == 'MiBIG']['normalized_similarity']
            random_raw = successful[successful['dataset_group'] == 'Random']['similarity_score']
            random_norm = successful[successful['dataset_group'] == 'Random']['normalized_similarity']
            
            if len(kegg_raw) > 0 and len(mibig_raw) > 0:
                print(f"  KEGG vs MiBIG:")
                print(f"    Raw scores   - KEGG: {kegg_raw.mean():.1f}±{kegg_raw.std():.1f}, MiBIG: {mibig_raw.mean():.1f}±{mibig_raw.std():.1f}, Diff: {kegg_raw.mean() - mibig_raw.mean():.1f}")
                print(f"    Normalized   - KEGG: {kegg_norm.mean():.1f}±{kegg_norm.std():.1f}, MiBIG: {mibig_norm.mean():.1f}±{mibig_norm.std():.1f}, Diff: {kegg_norm.mean() - mibig_norm.mean():.1f}")
            
            if len(kegg_raw) > 0 and len(random_raw) > 0:
                print(f"  KEGG vs Random:")
                print(f"    Raw scores   - KEGG: {kegg_raw.mean():.1f}±{kegg_raw.std():.1f}, Random: {random_raw.mean():.1f}±{random_raw.std():.1f}, Diff: {kegg_raw.mean() - random_raw.mean():.1f}")
                print(f"    Normalized   - KEGG: {kegg_norm.mean():.1f}±{kegg_norm.std():.1f}, Random: {random_norm.mean():.1f}±{random_norm.std():.1f}, Diff: {kegg_norm.mean() - random_norm.mean():.1f}")
            
            # Show normalization effects
            print(f"\nNORMALIZATION EFFECTS:")
            for group in ['KEGG', 'MiBIG', 'Random']:
                group_success = successful[successful['dataset_group'] == group]
                if len(group_success) > 0:
                    raw_mean = group_success['similarity_score'].mean()
                    norm_mean = group_success['normalized_similarity'].mean()
                    effect = norm_mean - raw_mean
                    print(f"  {group:6s}: {raw_mean:.1f} -> {norm_mean:.1f} (change: {effect:+.1f})")
            
            # Top performers by normalized score
            print(f"\nTOP 5 PERFORMERS BY NORMALIZED SCORE:")
            
            for group in ['KEGG', 'MiBIG', 'Random']:
                group_success = successful[successful['dataset_group'] == group]
                if len(group_success) > 0:
                    top_5 = group_success.nlargest(5, 'normalized_similarity')
                    print(f"\n  {group} Dataset:")
                    for i, (_, row) in enumerate(top_5.iterrows(), 1):
                        norm = row['normalized_similarity']
                        raw = row['similarity_score']
                        seqs = row['num_promoters']
                        print(f"    {i}. {row['group_name']}: {norm:.1f} (raw:{raw:.1f}, {seqs} seqs)")
        
        # Save results
        output_columns = [
            'group_name', 'dataset_group', 'num_promoters', 'num_comparisons',
            'similarity_score', 'normalized_similarity', 'similarity_category',
            'size_adjustment_factor', 'comparison_weight', 'std_adjustment_factor',
            'similarity_std', 'mean_proximal_similarity', 'mean_distal_similarity',
            'num_tfbs_types_found', 'status', 'error_msg'
        ]
        
        clean_df = df[[col for col in output_columns if col in df.columns]]
        clean_df.to_csv('promoter_similarity_results.csv', index=False)
        
        # Save separate files by group
        for group in ['KEGG', 'MiBIG', 'Random']:
            group_df = clean_df[clean_df['dataset_group'] == group]
            if len(group_df) > 0:
                filename = f'{group.lower()}_promoter_results.csv'
                group_df.to_csv(filename, index=False)
                print(f"Saved {group} results to: {filename}")
        
        print(f"Complete results saved to: promoter_similarity_results.csv")

if __name__ == "__main__":
    main()