#!/usr/bin/env python3
"""
E2P2 Feature Extractor - Extract EC number features from E2P2 .default.pf output files
"""

import pandas as pd
from pathlib import Path
import re
from collections import defaultdict
import argparse
import sys

class E2P2FeatureExtractor:
    def __init__(self, base_dir=None):
        if base_dir:
            self.base_dir = Path(base_dir)
        else:
            self.base_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test")
        
        self.e2p2_results_dir = self.base_dir / "e2p2_results"
        self.output_dir = self.base_dir / "e2p2_feature_extraction_results"
        self.output_dir.mkdir(exist_ok=True)
    
    def parse_e2p2_default_pf(self, e2p2_file):
        """Parse E2P2 results from .default.pf file and extract EC numbers"""
        results = []
        current_entry = {}
        
        try:
            with open(e2p2_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('ID\t'):
                        if current_entry and 'protein_id' in current_entry:
                            results.append(current_entry)
                        current_entry = {'protein_id': line[3:].strip()}
                    elif line.startswith('EC\t'):
                        # Handle multiple EC numbers per protein
                        if 'ec_numbers' not in current_entry:
                            current_entry['ec_numbers'] = []
                        current_entry['ec_numbers'].append(line[3:].strip())
                    elif line.startswith('NAME\t'):
                        current_entry['protein_name'] = line[5:].strip()
                    elif line.startswith('PRODUCT-TYPE\t'):
                        current_entry['product_type'] = line[13:].strip()
                    elif line == '//':
                        if current_entry and 'protein_id' in current_entry:
                            results.append(current_entry)
                        current_entry = {}
            
            # Add last entry if exists
            if current_entry and 'protein_id' in current_entry:
                results.append(current_entry)
                
        except Exception as e:
            print(f"Error parsing E2P2 file {e2p2_file}: {e}")
        
        return results
    
    def extract_ec_features(self, e2p2_results):
        """Extract EC number features from E2P2 results"""
        all_ec_numbers = []
        
        for result in e2p2_results:
            ec_numbers = result.get('ec_numbers', [])
            for ec_number in ec_numbers:
                if ec_number and ec_number != 'UNKNOWN':
                    # Clean EC number (remove any extra text)
                    ec_clean = re.sub(r'[^\d\.]', '', ec_number)
                    if re.match(r'^\d+\.\d+\.\d+\.\d+$', ec_clean):
                        all_ec_numbers.append(ec_clean)
        
        # Extract enzyme classes (first level: X)
        enzyme_classes = set()
        # Extract enzyme subclasses (first two levels: X.Y)
        enzyme_subclasses = set()
        # Extract enzyme families (first three levels: X.Y.Z)
        enzyme_families = set()
        # Extract enzyme subfamilies (all levels: X.Y.Z.W)
        enzyme_subfamilies = set()
        
        for ec in all_ec_numbers:
            parts = ec.split('.')
            if len(parts) >= 1:
                enzyme_classes.add(parts[0])
            if len(parts) >= 2:
                enzyme_subclasses.add(f"{parts[0]}.{parts[1]}")
            if len(parts) >= 3:
                enzyme_families.add(f"{parts[0]}.{parts[1]}.{parts[2]}")
            if len(parts) >= 4:
                enzyme_subfamilies.add(f"{parts[0]}.{parts[1]}.{parts[2]}.{parts[3]}")
        
        return {
            'num_distinct_enzyme_classes': len(enzyme_classes),
            'num_distinct_enzyme_subclasses': len(enzyme_subclasses),
            'num_distinct_enzyme_families': len(enzyme_families),
            'num_distinct_enzyme_subfamilies': len(enzyme_subfamilies),
            'enzyme_classes': sorted(list(enzyme_classes)),
            'enzyme_subclasses': sorted(list(enzyme_subclasses)),
            'enzyme_families': sorted(list(enzyme_families)),
            'enzyme_subfamilies': sorted(list(enzyme_subfamilies)),
            'total_ec_numbers': len(all_ec_numbers),
            'all_ec_numbers': all_ec_numbers
        }
    
    def determine_classification_label(self, cluster_name):
        """Determine classification label based on cluster name"""
        cluster_lower = cluster_name.lower()
        
        # Check for BGC or MGC_CANDIDATE patterns
        if any(pattern in cluster_lower for pattern in ['bgc', 'mgc_candidate']):
            return 1
        # Check for RANDOM patterns
        elif any(pattern in cluster_lower for pattern in ['random', 'random_mgc']):
            return 0
        else:
            return 1  # Default to BGC/MGC_CANDIDATE
    
    def process_single_cluster(self, cluster_dir):
        """Process a single cluster directory and extract features"""
        cluster_name = cluster_dir.name
        
        # Look for .default.pf file
        default_pf_file = cluster_dir / f"{cluster_name}.MaxWeightAbsoluteThreshold.default.pf"
        
        if not default_pf_file.exists():
            print(f"  ⚠ No .default.pf file found for {cluster_name}")
            return None
        
        try:
            # Parse E2P2 results
            e2p2_results = self.parse_e2p2_default_pf(default_pf_file)
            ec_features = self.extract_ec_features(e2p2_results)
            classification_label = self.determine_classification_label(cluster_name)
            
            # Create result entry
            result = {
                'cluster_name': cluster_name,
                'classification_label': classification_label,
                'num_distinct_enzyme_classes': ec_features['num_distinct_enzyme_classes'],
                'num_distinct_enzyme_subclasses': ec_features['num_distinct_enzyme_subclasses'],
                'num_distinct_enzyme_families': ec_features['num_distinct_enzyme_families'],
                'num_distinct_enzyme_subfamilies': ec_features['num_distinct_enzyme_subfamilies'],
                'total_ec_numbers': ec_features['total_ec_numbers'],
                'enzyme_classes': ';'.join(ec_features['enzyme_classes']),
                'enzyme_subclasses': ';'.join(ec_features['enzyme_subclasses']),
                'enzyme_families': ';'.join(ec_features['enzyme_families']),
                'enzyme_subfamilies': ';'.join(ec_features['enzyme_subfamilies']),
                'all_ec_numbers': ';'.join(ec_features['all_ec_numbers']),
                'e2p2_file': str(default_pf_file)
            }
            
            print(f"  ✓ {cluster_name}: {ec_features['num_distinct_enzyme_classes']} enzyme classes, {ec_features['num_distinct_enzyme_subfamilies']} subfamilies")
            return result
            
        except Exception as e:
            print(f"  ✗ Error processing {cluster_name}: {e}")
            return None
    
    def process_all_clusters(self):
        """Process all cluster directories and extract features"""
        if not self.e2p2_results_dir.exists():
            print(f"E2P2 results directory not found: {self.e2p2_results_dir}")
            return None
        
        # Find all cluster directories
        cluster_dirs = [d for d in self.e2p2_results_dir.iterdir() if d.is_dir()]
        
        if not cluster_dirs:
            print(f"No cluster directories found in {self.e2p2_results_dir}")
            return None
        
        print(f"Processing {len(cluster_dirs)} clusters...")
        
        results = []
        successful = 0
        failed = 0
        
        for cluster_dir in cluster_dirs:
            result = self.process_single_cluster(cluster_dir)
            if result:
                results.append(result)
                successful += 1
            else:
                failed += 1
        
        # Create DataFrame and save
        if results:
            results_df = pd.DataFrame(results)
            
            # Save detailed results
            detailed_file = self.output_dir / "e2p2_feature_extraction_detailed.csv"
            results_df.to_csv(detailed_file, index=False)
            
            # Create summary
            self.create_summary(results_df)
            
            print(f"\nFeature extraction completed:")
            print(f"  - Successful: {successful}")
            print(f"  - Failed: {failed}")
            print(f"  - Detailed results saved to: {detailed_file}")
            
            return results_df
        else:
            print("No results to process")
            return None
    
    def create_summary(self, results_df):
        """Create summary statistics"""
        # Calculate summary statistics
        total_clusters = len(results_df)
        avg_enzyme_classes = results_df['num_distinct_enzyme_classes'].mean()
        avg_enzyme_subclasses = results_df['num_distinct_enzyme_subclasses'].mean()
        avg_enzyme_families = results_df['num_distinct_enzyme_families'].mean()
        avg_enzyme_subfamilies = results_df['num_distinct_enzyme_subfamilies'].mean()
        avg_total_ec = results_df['total_ec_numbers'].mean()
        
        # Count classification labels
        bgc_mgc_count = len(results_df[results_df['classification_label'] == 1])
        random_count = len(results_df[results_df['classification_label'] == 0])
        
        # Create summary
        summary = {
            'total_clusters_analyzed': total_clusters,
            'avg_distinct_enzyme_classes': round(avg_enzyme_classes, 2),
            'avg_distinct_enzyme_subclasses': round(avg_enzyme_subclasses, 2),
            'avg_distinct_enzyme_families': round(avg_enzyme_families, 2),
            'avg_distinct_enzyme_subfamilies': round(avg_enzyme_subfamilies, 2),
            'avg_total_ec_numbers': round(avg_total_ec, 2),
            'bgc_mgc_candidate_clusters': bgc_mgc_count,
            'random_clusters': random_count,
            'classification_ratio': round(bgc_mgc_count / total_clusters, 3) if total_clusters > 0 else 0
        }
        
        # Save summary
        summary_file = self.output_dir / "e2p2_feature_extraction_summary.csv"
        summary_df = pd.DataFrame([summary])
        summary_df.to_csv(summary_file, index=False)
        
        print(f"  ✓ Summary saved to: {summary_file}")
        print(f"    - Average distinct enzyme classes: {avg_enzyme_classes:.2f}")
        print(f"    - Average distinct enzyme subclasses: {avg_enzyme_subclasses:.2f}")
        print(f"    - Average distinct enzyme families: {avg_enzyme_families:.2f}")
        print(f"    - Average distinct enzyme subfamilies: {avg_enzyme_subfamilies:.2f}")
        print(f"    - Average total EC number: {avg_total_ec:.2f}")
        print(f"    - BGC/MGC_CANDIDATE clusters: {bgc_mgc_count}")
        print(f"    - RANDOM clusters: {random_count}")

def main():
    parser = argparse.ArgumentParser(description='Extract EC features from E2P2 .default.pf files')
    parser.add_argument('--base_dir', type=str, help='Base directory containing e2p2_results')
    parser.add_argument('--output_dir', type=str, help='Output directory for results')
    
    args = parser.parse_args()
    
    # Initialize extractor
    extractor = E2P2FeatureExtractor(base_dir=args.base_dir)
    
    if args.output_dir:
        extractor.output_dir = Path(args.output_dir)
        extractor.output_dir.mkdir(exist_ok=True)
    
    print("E2P2 Feature Extractor")
    print("=" * 40)
    print(f"Base directory: {extractor.base_dir}")
    print(f"E2P2 results directory: {extractor.e2p2_results_dir}")
    print(f"Output directory: {extractor.output_dir}")
    print()
    
    # Process all clusters
    results_df = extractor.process_all_clusters()
    
    if results_df is not None:
        print(f"\nKey features extracted:")
        print(f"  1. Number of distinct enzyme classes (X): {results_df['num_distinct_enzyme_classes'].mean():.2f} average")
        print(f"  2. Number of distinct enzyme subclasses (X.Y): {results_df['num_distinct_enzyme_subclasses'].mean():.2f} average")
        print(f"  3. Number of distinct enzyme families (X.Y.Z): {results_df['num_distinct_enzyme_families'].mean():.2f} average")
        print(f"  4. Number of distinct enzyme subfamilies (X.Y.Z.W): {results_df['num_distinct_enzyme_subfamilies'].mean():.2f} average")
        print(f"  5. Total EC numbers: {results_df['total_ec_numbers'].mean():.2f} average")
        print(f"  6. Classification labels: 1=BGC/MGC_CANDIDATE, 0=RANDOM")
        print(f"     - BGC/MGC_CANDIDATE: {len(results_df[results_df['classification_label'] == 1])}")
        print(f"     - RANDOM: {len(results_df[results_df['classification_label'] == 0])}")

if __name__ == "__main__":
    main()
