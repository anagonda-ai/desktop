#!/usr/bin/env python3
"""
Protein Domain Detector - ONLY Domain Detection
No classification, just domain detection using Pfam and InterProScan
"""

import subprocess
import json
from pathlib import Path
from typing import Dict, List
import logging
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
import time

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class DomainHit:
    """Represents a protein domain hit"""
    protein_id: str
    domain_id: str
    domain_name: str
    domain_description: str
    start: int
    end: int
    e_value: float
    score: float
    database: str
    confidence: str = "low"

class DomainDetector:
    """
    Protein domain detection system - ONLY domain detection
    No classification, just finds domains using Pfam and InterProScan
    """
    
    def __init__(self, 
                 hmmer_path: str = "hmmscan",
                 interproscan_path: str = "interproscan.sh",
                 pfam_db_path: str = None,
                 temp_dir: str = None):
        """
        Initialize the domain detector
        
        Args:
            hmmer_path: Path to hmmscan executable
            interproscan_path: Path to InterProScan executable
            pfam_db_path: Path to Pfam HMM database
            temp_dir: Temporary directory for intermediate files
        """
        self.hmmer_path = hmmer_path
        self.interproscan_path = interproscan_path
        self.pfam_db_path = pfam_db_path or "/groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/plantismash/antismash/generic_modules/fullhmmer/Pfam-A.hmm"
        self.temp_dir = temp_dir or "/tmp/domain_detection"
        
        # Create temp directory
        Path(self.temp_dir).mkdir(parents=True, exist_ok=True)
    
    def check_dependencies(self) -> Dict[str, bool]:
        """
        Check if required tools are available
        
        Returns:
            Dictionary with tool availability status
        """
        dependencies = {}
        
        # Check HMMER
        try:
            result = subprocess.run([self.hmmer_path, "-h"], 
                                  capture_output=True, text=True, timeout=10)
            dependencies['hmmer'] = result.returncode == 0
        except:
            dependencies['hmmer'] = False
        
        # Check InterProScan
        try:
            result = subprocess.run([self.interproscan_path, "-version"], 
                                  capture_output=True, text=True, timeout=10)
            dependencies['interproscan'] = result.returncode == 0
        except:
            dependencies['interproscan'] = False
        
        # Check Pfam database
        dependencies['pfam_db'] = Path(self.pfam_db_path).exists()
        
        return dependencies
    
    def run_pfam_analysis(self, fasta_file: str) -> List[DomainHit]:
        """
        Run Pfam analysis using HMMER
        
        Args:
            fasta_file: Input FASTA file path
            
        Returns:
            List of DomainHit objects
        """
        logger.info("Running Pfam analysis...")
        
        output_file = Path(self.temp_dir) / f"{Path(fasta_file).stem}_pfam.txt"
        
        try:
            # Run hmmscan with Pfam database
            cmd = [
                self.hmmer_path,
                "--domtblout", str(output_file),
                "--cut_ga",  # Use gathering thresholds
                "--cpu", "4",
                self.pfam_db_path,
                fasta_file
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            
            if result.returncode != 0:
                logger.error(f"HMMER failed: {result.stderr}")
                return []
            
            # Parse results
            return self._parse_hmmer_output(output_file)
            
        except subprocess.TimeoutExpired:
            logger.error("HMMER scan timed out")
            return []
        except Exception as e:
            logger.error(f"Error running Pfam analysis: {e}")
            return []
    
    def run_interproscan_analysis(self, fasta_file: str) -> List[DomainHit]:
        """
        Run InterProScan analysis
        
        Args:
            fasta_file: Input FASTA file path
            
        Returns:
            List of DomainHit objects
        """
        logger.info("Running InterProScan analysis...")
        
        output_file = Path(self.temp_dir) / f"{Path(fasta_file).stem}_interpro.tsv"
        
        try:
            # Run InterProScan with limited databases to avoid database issues
            cmd = [
                self.interproscan_path,
                "-i", fasta_file,
                "-o", str(output_file),
                "-f", "tsv",
                "--cpu", "4",
                "--applications", "PfamA",  # Only use PfamA to avoid database issues
                "--goterms",
                "--pathways"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
            
            if result.returncode != 0:
                logger.error(f"InterProScan failed: {result.stderr}")
                return []
            
            # Check if output file exists and has content
            if not output_file.exists() or output_file.stat().st_size == 0:
                logger.info("InterProScan completed but found no additional domains")
                return []
            
            # Parse results
            return self._parse_interproscan_output(output_file)
            
        except subprocess.TimeoutExpired:
            logger.error("InterProScan timed out")
            return []
        except Exception as e:
            logger.error(f"Error running InterProScan: {e}")
            return []
    
    def _parse_hmmer_output(self, output_file: Path) -> List[DomainHit]:
        """
        Parse HMMER output file
        
        Args:
            output_file: Path to HMMER output file
            
        Returns:
            List of DomainHit objects
        """
        domain_hits = []
        
        try:
            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split()
                    if len(parts) >= 22:
                        domain_hit = DomainHit(
                            protein_id=parts[3],
                            domain_id=parts[0],
                            domain_name=parts[0],
                            domain_description=parts[22] if len(parts) > 22 else '',
                            start=int(parts[19]),
                            end=int(parts[20]),
                            e_value=float(parts[6]),
                            score=float(parts[7]),
                            database="Pfam",
                            confidence=self._calculate_confidence(float(parts[6]), float(parts[7]))
                        )
                        domain_hits.append(domain_hit)
            
            return domain_hits
            
        except Exception as e:
            logger.error(f"Error parsing HMMER output: {e}")
            return []
    
    def _parse_interproscan_output(self, output_file: Path) -> List[DomainHit]:
        """
        Parse InterProScan output file
        
        Args:
            output_file: Path to InterProScan output file
            
        Returns:
            List of DomainHit objects
        """
        domain_hits = []
        
        try:
            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) >= 12:
                        domain_hit = DomainHit(
                            protein_id=parts[0],
                            domain_id=parts[4],
                            domain_name=parts[4],
                            domain_description=parts[5],
                            start=int(parts[6]),
                            end=int(parts[7]),
                            e_value=float(parts[8]) if parts[8] != '-' else 1.0,
                            score=float(parts[8]) if parts[8] != '-' else 0.0,
                            database=parts[3],
                            confidence=self._calculate_confidence(
                                float(parts[8]) if parts[8] != '-' else 1.0,
                                float(parts[8]) if parts[8] != '-' else 0.0
                            )
                        )
                        domain_hits.append(domain_hit)
            
            return domain_hits
            
        except Exception as e:
            logger.error(f"Error parsing InterProScan output: {e}")
            return []
    
    def _calculate_confidence(self, e_value: float, score: float) -> str:
        """
        Calculate confidence level based on E-value and score
        
        Args:
            e_value: E-value from domain search
            score: Score from domain search
            
        Returns:
            Confidence level string
        """
        if e_value < 1e-10 and score > 50:
            return "high"
        elif e_value < 1e-5 and score > 30:
            return "medium"
        else:
            return "low"
    
    def detect_protein_domains(self, fasta_file: str) -> Dict:
        """
        Detect protein domains using Pfam and InterProScan
        
        Args:
            fasta_file: Input FASTA file path
            
        Returns:
            Dictionary with domain detection results
        """
        logger.info(f"Starting domain detection for {fasta_file}")
        
        # Check dependencies
        deps = self.check_dependencies()
        
        results = {
            'fasta_file': fasta_file,
            'domain_hits': [],
            'all_domains': [],
            'confidence_score': 0.0,
            'analysis_methods': []
        }
        
        # Run Pfam analysis
        if deps.get('hmmer', False) and deps.get('pfam_db', False):
            pfam_hits = self.run_pfam_analysis(fasta_file)
            results['domain_hits'].extend(pfam_hits)
            results['analysis_methods'].append('Pfam')
            logger.info(f"Pfam analysis: {len(pfam_hits)} domain hits")
        
        # Run InterProScan analysis (now enabled with fixed databases)
        if deps.get('interproscan', False):
            try:
                interpro_hits = self.run_interproscan_analysis(fasta_file)
                if interpro_hits:  # Only add if successful
                    results['domain_hits'].extend(interpro_hits)
                    results['analysis_methods'].append('InterProScan')
                    logger.info(f"InterProScan analysis: {len(interpro_hits)} domain hits")
                else:
                    logger.warning("InterProScan failed - using Pfam only")
            except Exception as e:
                logger.warning(f"InterProScan failed: {e} - using Pfam only")
        
        # Extract domain names
        if results['domain_hits']:
            domain_names = [hit.domain_name for hit in results['domain_hits']]
            results['all_domains'] = domain_names
            
            # Calculate overall confidence
            results['confidence_score'] = self._calculate_overall_confidence(results['domain_hits'])
        
        logger.info(f"Domain detection completed: {len(results['domain_hits'])} domain hits found")
        return results
    
    def _calculate_overall_confidence(self, domain_hits: List[DomainHit]) -> float:
        """
        Calculate overall confidence score for the analysis
        
        Args:
            domain_hits: List of domain hits
            
        Returns:
            Confidence score (0.0 to 1.0)
        """
        if not domain_hits:
            return 0.0
        
        # Weight by confidence level
        confidence_weights = {'high': 1.0, 'medium': 0.7, 'low': 0.3}
        weighted_scores = []
        
        for hit in domain_hits:
            weight = confidence_weights.get(hit.confidence, 0.3)
            # Convert E-value to score (lower E-value = higher score)
            e_score = max(0, 1 - (hit.e_value / 1e-10) if hit.e_value > 0 else 1)
            weighted_scores.append(weight * e_score)
        
        return sum(weighted_scores) / len(weighted_scores) if weighted_scores else 0.0

def process_single_fasta(detector, fasta_file, output_dir, file_index, total_files):
    """Process a single FASTA file and return results"""
    try:
        print(f"[{file_index}/{total_files}] Processing: {fasta_file.name}")
        
        # Run domain detection
        results = detector.detect_protein_domains(fasta_file)
        
        # Save results to JSON in output directory
        output_file = output_dir / f"{fasta_file.stem}_domains.json"
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        # Return results for summary
        result_info = {
            'fasta_file': str(fasta_file),
            'fasta_name': fasta_file.name,
            'fasta_stem': fasta_file.stem,
            'domain_hits': len(results['domain_hits']),
            'domains': ', '.join(results.get('all_domains', [])),
            'confidence_score': results['confidence_score'],
            'analysis_methods': ', '.join(results['analysis_methods']),
            'output_file': str(output_file),
            'success': True
        }
        
        print(f"  ✓ [{file_index}/{total_files}] {fasta_file.name}: {len(results['domain_hits'])} domains, confidence: {results['confidence_score']:.3f}")
        return result_info
        
    except Exception as e:
        print(f"  ✗ [{file_index}/{total_files}] Error processing {fasta_file.name}: {e}")
        return {
            'fasta_file': str(fasta_file),
            'fasta_name': fasta_file.name,
            'fasta_stem': fasta_file.stem,
            'domain_hits': 0,
            'domains': '',
            'confidence_score': 0.0,
            'analysis_methods': '',
            'output_file': '',
            'success': False,
            'error': str(e)
        }

def process_all_fasta_directories():
    """Process all FASTA files in both directories"""
    detector = DomainDetector()
    
    # Check dependencies first
    deps = detector.check_dependencies()
    print("Dependency check:")
    for tool, available in deps.items():
        print(f"  {tool}: {'✓' if available else '✗'}")
    
    # Set up output directory (same pattern as e2p2)
    base_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test")
    output_dir = base_dir / "domain_detection_results"
    output_dir.mkdir(exist_ok=True)
    
    fasta_directories = [
        base_dir / "mgc_candidates_fasta_files",
        base_dir / "random_mgc_candidates_fasta_files"
    ]
    
    all_results = []
    all_fasta_files = []
    
    # Collect all FASTA files from both directories
    for fasta_dir in fasta_directories:
        if not fasta_dir.exists():
            print(f"Directory not found: {fasta_dir}")
            continue
        
        print(f"\n=== Scanning FASTA directory: {fasta_dir.name} ===")
        
        # Collect FASTA files
        fasta_files = list(fasta_dir.glob("*.fasta"))
        
        if not fasta_files:
            print(f"No FASTA files found in {fasta_dir.name}")
            continue
        
        print(f"Found {len(fasta_files)} FASTA files")
        all_fasta_files.extend(fasta_files)
    
    if not all_fasta_files:
        print("No FASTA files found in any directory!")
        return
    
    print(f"\n=== Processing {len(all_fasta_files)} FASTA files with concurrency ===")
    
    # Process files with concurrency (use CPU count as max workers)
    max_workers = min(multiprocessing.cpu_count(), len(all_fasta_files))
    print(f"Using {max_workers} parallel workers (CPU count: {multiprocessing.cpu_count()})")
    
    start_time = time.time()
    successful_jobs = 0
    failed_jobs = 0
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_fasta = {
            executor.submit(process_single_fasta, detector, fasta_file, output_dir, i+1, len(all_fasta_files)): fasta_file 
            for i, fasta_file in enumerate(all_fasta_files)
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_fasta):
            fasta_file = future_to_fasta[future]
            try:
                result = future.result()
                all_results.append(result)
                
                if result['success']:
                    successful_jobs += 1
                else:
                    failed_jobs += 1
                    
            except Exception as e:
                print(f"  ✗ Unexpected error processing {fasta_file.name}: {e}")
                failed_jobs += 1
                all_results.append({
                    'fasta_file': str(fasta_file),
                    'fasta_name': fasta_file.name,
                    'fasta_stem': fasta_file.stem,
                    'domain_hits': 0,
                    'domains': '',
                    'confidence_score': 0.0,
                    'analysis_methods': '',
                    'output_file': '',
                    'success': False,
                    'error': str(e)
                })
    
    elapsed_time = time.time() - start_time
    print(f"\n=== Processing Complete ===")
    print(f"✓ Successful: {successful_jobs}")
    print(f"✗ Failed: {failed_jobs}")
    print(f"⏱️  Total time: {elapsed_time:.1f} seconds")
    print(f"⚡ Average time per file: {elapsed_time/len(all_fasta_files):.1f} seconds")
    
    # Save summary CSV
    if all_results:
        import pandas as pd
        summary_df = pd.DataFrame(all_results)
        summary_file = output_dir / "domain_detection_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        print(f"\n✓ Summary saved: {summary_file}")
        print(f"✓ Total files processed: {len(all_results)}")
        print(f"✓ Results directory: {output_dir}")

def main():
    """Process all FASTA files in both directories"""
    print("Domain Detection Batch Processor")
    print("=" * 40)
    
    process_all_fasta_directories()
    
    print("\nAll domain detection complete!")

if __name__ == "__main__":
    main()
