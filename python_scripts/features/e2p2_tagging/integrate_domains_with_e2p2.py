#!/usr/bin/env python3
"""
Simple Integration Script
Adds protein domain detection to your existing E2P2 results
"""

import sys
import pandas as pd
import json
from pathlib import Path
import logging

# Import domain detector
import importlib.util
spec = importlib.util.spec_from_file_location("domain_detector", "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/protein_domain_detection/domain_detector_only.py")
domain_detector = importlib.util.module_from_spec(spec)
spec.loader.exec_module(domain_detector)
DomainDetector = domain_detector.DomainDetector

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def parse_e2p2_results(e2p2_file: str) -> pd.DataFrame:
    """Parse your existing E2P2 results"""
    results = []
    
    with open(e2p2_file, 'r') as f:
        current_protein = None
        current_ec = None
        
        for line in f:
            line = line.strip()
            
            if line.startswith('>'):
                if current_protein and current_ec:
                    results.append({
                        'protein_id': current_protein,
                        'ec_number': current_ec
                    })
                current_protein = line[1:].split()[0]
                current_ec = None
            elif line.startswith('EC:'):
                current_ec = line.split(':')[1].strip()
        
        if current_protein and current_ec:
            results.append({
                'protein_id': current_protein,
                'ec_number': current_ec
            })
    
    return pd.DataFrame(results)

def classify_enzyme_with_ec_and_domains(domains: list, ec_number: str) -> str:
    """
    Classify enzyme as signature or tailoring using BOTH EC numbers AND domains
    
    Args:
        domains: List of detected domain names
        ec_number: EC number from E2P2
        
    Returns:
        Classification: 'signature', 'tailoring', or 'unknown'
    """
    # EC number patterns for classification
    signature_ec_patterns = [
        '2.3.1.',  # Acyltransferases (PKS)
        '6.1.1.',  # Aminoacyl-tRNA synthetases (NRPS)
        '4.2.1.',  # Carbon-carbon lyases (terpene synthases)
        '2.3.3.',  # Acyltransferases (fatty acid synthases)
        '4.1.1.',  # Carboxy-lyases
    ]
    
    tailoring_ec_patterns = [
        '2.1.1.',  # Methyltransferases
        '1.14.14.', # Cytochrome P450 monooxygenases
        '2.4.1.',  # Glycosyltransferases
        '1.14.13.', # Monooxygenases
        '2.7.1.',  # Phosphotransferases
    ]
    
    # Domain classification
    signature_domains = {
        'PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_ER', 'PKS_KR', 'PKS_ACP',
        'NRPS_A', 'NRPS_C', 'NRPS_PCP', 'NRPS_E', 'NRPS_TE',
        'Terpene_synth_C', 'Terpene_synth_N', 'FAS', 'KS', 'AT', 'DH', 'ER', 'KR', 'ACP'
    }
    
    tailoring_domains = {
        'Methyltransf_1', 'Methyltransf_2', 'Methyltransf_3',
        'P450', 'FAD_binding_3', 'FAD_binding_4', 'Cytochrome_P450',
        'Glycos_transf_1', 'Glycos_transf_2', 'GT1', 'GT2',
        'Acetyltransf_1', 'Acetyltransf_2', 'GNAT'
    }
    
    domain_set = set(domains)
    signature_domain_matches = domain_set.intersection(signature_domains)
    tailoring_domain_matches = domain_set.intersection(tailoring_domains)
    
    # Check EC number first (higher priority)
    if ec_number:
        for pattern in signature_ec_patterns:
            if ec_number.startswith(pattern):
                return 'signature'
        for pattern in tailoring_ec_patterns:
            if ec_number.startswith(pattern):
                return 'tailoring'
    
    # If EC number doesn't match, use domains
    if signature_domain_matches and not tailoring_domain_matches:
        return 'signature'
    elif tailoring_domain_matches and not signature_domain_matches:
        return 'tailoring'
    elif signature_domain_matches and tailoring_domain_matches:
        # Both present - use EC number if available, otherwise default to signature
        return 'signature'
    else:
        return 'unknown'

def classify_by_ec_only(ec_number: str) -> str:
    """Classify using only EC number"""
    signature_ec_patterns = ['2.3.1.', '6.1.1.', '4.2.1.', '2.3.3.', '4.1.1.']
    tailoring_ec_patterns = ['2.1.1.', '1.14.14.', '2.4.1.', '1.14.13.', '2.7.1.']
    
    if ec_number:
        for pattern in signature_ec_patterns:
            if ec_number.startswith(pattern):
                return 'signature'
        for pattern in tailoring_ec_patterns:
            if ec_number.startswith(pattern):
                return 'tailoring'
    return 'unknown'

def classify_by_domains_only(domains: list) -> str:
    """Classify using only domains"""
    signature_domains = {
        'PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_ER', 'PKS_KR', 'PKS_ACP',
        'NRPS_A', 'NRPS_C', 'NRPS_PCP', 'NRPS_E', 'NRPS_TE',
        'Terpene_synth_C', 'Terpene_synth_N', 'FAS', 'KS', 'AT', 'DH', 'ER', 'KR', 'ACP'
    }
    tailoring_domains = {
        'Methyltransf_1', 'Methyltransf_2', 'Methyltransf_3',
        'P450', 'FAD_binding_3', 'FAD_binding_4', 'Cytochrome_P450',
        'Glycos_transf_1', 'Glycos_transf_2', 'GT1', 'GT2',
        'Acetyltransf_1', 'Acetyltransf_2', 'GNAT'
    }
    
    domain_set = set(domains)
    signature_matches = domain_set.intersection(signature_domains)
    tailoring_matches = domain_set.intersection(tailoring_domains)
    
    if signature_matches and not tailoring_matches:
        return 'signature'
    elif tailoring_matches and not signature_matches:
        return 'tailoring'
    elif signature_matches and tailoring_matches:
        return 'mixed'
    else:
        return 'unknown'

def calculate_confidence_score(ec_number: str, domains: list, signature_domains: list, tailoring_domains: list) -> float:
    """Calculate confidence score based on EC + domains agreement"""
    ec_class = classify_by_ec_only(ec_number)
    domain_class = classify_by_domains_only(domains)
    
    # High confidence if both agree
    if ec_class == domain_class and ec_class != 'unknown':
        return 0.9
    # Medium confidence if one is unknown
    elif ec_class == 'unknown' or domain_class == 'unknown':
        return 0.6
    # Low confidence if they disagree
    elif ec_class != domain_class:
        return 0.3
    else:
        return 0.0

def add_domain_analysis_to_e2p2(e2p2_file: str, fasta_file: str, output_file: str):
    """
    Add domain analysis to your existing E2P2 results
    
    Args:
        e2p2_file: Your existing E2P2 results file
        fasta_file: Original FASTA file
        output_file: Output file with E2P2 + domain results
    """
    
    logger.info(f"Adding domain analysis to E2P2 results...")
    
    # Initialize domain detector
    detector = DomainDetector()
    
    # Check if tools are available
    deps = detector.check_dependencies()
    logger.info("Available tools:")
    for tool, available in deps.items():
        logger.info(f"  {tool}: {'✓' if available else '✗'}")
    
    # Parse your E2P2 results
    e2p2_df = parse_e2p2_results(e2p2_file)
    logger.info(f"Found {len(e2p2_df)} proteins with EC numbers")
    
    # Run domain analysis on the FASTA file
    if any(deps.values()):
        logger.info("Running domain analysis...")
        domain_results = detector.detect_protein_domains(fasta_file)
        
        # Create integrated results using E2P2 + domain detection
        integrated_results = []
        
        for _, row in e2p2_df.iterrows():
            protein_id = row['protein_id']
            ec_number = row['ec_number']
            
            # Find domains for this protein from domain detection
            protein_domains = []
            signature_domains = []
            tailoring_domains = []
            
            # Find domains for this protein from domain detection
            for domain_hit in domain_results.get('domain_hits', []):
                if domain_hit.get('protein_id') == protein_id:
                    domain_name = domain_hit.get('domain_name', '')
                    protein_domains.append(domain_name)
            
            # Categorize domains as signature or tailoring (THIS IS WHERE CLASSIFICATION HAPPENS)
            signature_domains = []
            tailoring_domains = []
            
            for domain_name in protein_domains:
                if domain_name in {
                    'PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_ER', 'PKS_KR', 'PKS_ACP',
                    'NRPS_A', 'NRPS_C', 'NRPS_PCP', 'NRPS_E', 'NRPS_TE',
                    'Terpene_synth_C', 'Terpene_synth_N', 'FAS', 'KS', 'AT', 'DH', 'ER', 'KR', 'ACP'
                }:
                    signature_domains.append(domain_name)
                elif domain_name in {
                    'Methyltransf_1', 'Methyltransf_2', 'Methyltransf_3',
                    'P450', 'FAD_binding_3', 'FAD_binding_4', 'Cytochrome_P450',
                    'Glycos_transf_1', 'Glycos_transf_2', 'GT1', 'GT2',
                    'Acetyltransf_1', 'Acetyltransf_2', 'GNAT'
                }:
                    tailoring_domains.append(domain_name)
            
            # Classify using E2P2 EC numbers + domain detection
            enzyme_type = classify_enzyme_with_ec_and_domains(protein_domains, ec_number)
            
            # Calculate confidence based on both EC and domains
            confidence = calculate_confidence_score(ec_number, protein_domains, signature_domains, tailoring_domains)
            
            integrated_results.append({
                'protein_id': protein_id,
                'ec_number': ec_number,
                'domains': ','.join(protein_domains),
                'enzyme_classification': enzyme_type,
                'signature_domains': ','.join(signature_domains),
                'tailoring_domains': ','.join(tailoring_domains),
                'confidence_score': confidence,
                'ec_based_classification': classify_by_ec_only(ec_number),
                'domain_based_classification': classify_by_domains_only(protein_domains)
            })
        
        # Save integrated results
        integrated_df = pd.DataFrame(integrated_results)
        integrated_df.to_csv(output_file, index=False)
        
        logger.info(f"Integrated results saved to {output_file}")
        
        # Print summary
        signature_count = len(integrated_df[integrated_df['enzyme_classification'] == 'signature'])
        tailoring_count = len(integrated_df[integrated_df['enzyme_classification'] == 'tailoring'])
        
        logger.info(f"Classification summary:")
        logger.info(f"  Signature enzymes: {signature_count}")
        logger.info(f"  Tailoring enzymes: {tailoring_count}")
        
    else:
        logger.warning("No domain detection tools available. Saving E2P2 results only.")
        e2p2_df.to_csv(output_file, index=False)

def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Add domain analysis to E2P2 results')
    parser.add_argument('--e2p2_file', required=True, help='Your E2P2 results file')
    parser.add_argument('--fasta_file', required=True, help='Original FASTA file')
    parser.add_argument('--output', required=True, help='Output file with integrated results')
    
    args = parser.parse_args()
    
    add_domain_analysis_to_e2p2(args.e2p2_file, args.fasta_file, args.output)

if __name__ == "__main__":
    main()
