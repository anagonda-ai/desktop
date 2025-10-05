#!/usr/bin/env python3
"""
Single Cluster Integration - Process one cluster with dynamic classification
"""

import sys
import pandas as pd
import json
from pathlib import Path
import logging
from dynamic_classifier import DynamicEnzymeClassifier, ClassificationResult

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s:%(message)s')
logger = logging.getLogger(__name__)

def is_valid_protein_id(protein_id: str) -> bool:
    """Check if protein ID is valid (not metadata)"""
    invalid_patterns = ['EC', 'METACYC', 'NAME', 'PRODUCT-TYPE']
    return not any(pattern in protein_id for pattern in invalid_patterns)

def parse_e2p2_results(e2p2_file):
    """Parse E2P2 results from .default.pf file"""
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
                    current_entry['ec_number'] = line[3:].strip()
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
        logger.error(f"Error parsing E2P2 file {e2p2_file}: {e}")
    
    return results

def load_domain_results(domain_file):
    """Load domain analysis results from JSON file"""
    try:
        with open(domain_file, 'r') as f:
            data = json.load(f)
        return data.get('domain_hits', [])
    except Exception as e:
        logger.error(f"Error loading domain file {domain_file}: {e}")
        return []

def load_fasta_proteins(fasta_file):
    """Extract protein IDs from FASTA file"""
    proteins = []
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    protein_id = line[1:].strip().split()[0]  # Get first part of header
                    proteins.append(protein_id)
    except Exception as e:
        logger.error(f"Error loading FASTA file {fasta_file}: {e}")
    return proteins

def classify_enzyme_dynamically(domains: list, ec_number: str) -> str:
    """Classify enzyme using dynamic web-based classification"""
    try:
        classifier = DynamicEnzymeClassifier()
        result = classifier.classify_enzyme(ec_number, domains)
        logger.info(f"Dynamic classification: EC {ec_number}, Domains {domains} -> {result.result.value} (confidence: {result.confidence})")
        return result.result.value
    except Exception as e:
        logger.warning(f"Dynamic classification failed: {e}. Using fallback classification.")
        # Simple fallback classification
        if ec_number.startswith('2.1.1.') or 'methyl' in str(domains).lower():
            return 'tailoring'
        elif ec_number.startswith('2.3.1.') or 'synth' in str(domains).lower():
            return 'signature'
        else:
            return 'unknown'

def calculate_confidence_score(ec_number, protein_domains, signature_domains, tailoring_domains):
    """Calculate confidence score for classification"""
    score = 0.5  # Base score
    
    # EC number confidence
    if ec_number and ec_number != 'UNKNOWN':
        score += 0.2
    
    # Domain confidence
    if signature_domains and tailoring_domains:
        score += 0.3  # Mixed evidence
    elif signature_domains or tailoring_domains:
        score += 0.2  # Single type of evidence
    
    return min(score, 1.0)

def integrate_cluster_data(cluster_name, e2p2_file, domain_file, fasta_file, output_file):
    """Integrate E2P2 and domain data for a single cluster"""
    logger.info(f"Processing cluster: {cluster_name}")
    
    # Load data
    e2p2_results = parse_e2p2_results(e2p2_file)
    domain_hits = load_domain_results(domain_file)
    fasta_proteins = load_fasta_proteins(fasta_file)
    
    logger.info(f"Found {len(e2p2_results)} E2P2 results")
    logger.info(f"Found {len(domain_hits)} domain hits")
    logger.info(f"Found {len(fasta_proteins)} proteins in FASTA")
    
    # Create domain lookup
    domain_lookup = {}
    for hit in domain_hits:
        protein_id = hit.get('protein_id')
        if protein_id not in domain_lookup:
            domain_lookup[protein_id] = []
        domain_lookup[protein_id].append(hit.get('domain_name', ''))
    
    # Integrate data
    integrated_results = []
    
    for e2p2_result in e2p2_results:
        protein_id = e2p2_result.get('protein_id', '')
        ec_number = e2p2_result.get('ec_number', 'UNKNOWN')
        
        # Skip invalid protein IDs
        if not is_valid_protein_id(protein_id):
            continue
        
        # Get domains for this protein
        protein_domains = domain_lookup.get(protein_id, [])
        
        # Classify domains
        signature_domains = []
        tailoring_domains = []
        
        for domain in protein_domains:
            if any(pattern in domain.lower() for pattern in ['pks', 'nrps', 'synth', 'transferase', 'ligase']):
                signature_domains.append(domain)
            elif any(pattern in domain.lower() for pattern in ['methyl', 'hydroxyl', 'glycosyl', 'acetyl', 'kinase', 'phosphatase']):
                tailoring_domains.append(domain)
        
        # Classify enzyme
        enzyme_type = classify_enzyme_dynamically(protein_domains, ec_number)
        confidence = calculate_confidence_score(ec_number, protein_domains, signature_domains, tailoring_domains)
        
        integrated_results.append({
            'cluster_id': cluster_name,
            'protein_id': protein_id,
            'protein_name': e2p2_result.get('protein_name', ''),
            'ec_number': ec_number,
            'product_type': e2p2_result.get('product_type', ''),
            'domains': ';'.join(protein_domains) if protein_domains else '',
            'signature_domains': ';'.join(signature_domains) if signature_domains else '',
            'tailoring_domains': ';'.join(tailoring_domains) if tailoring_domains else '',
            'enzyme_classification': enzyme_type,
            'confidence_score': confidence,
            'analysis_methods': 'dynamic_classifier'
        })
    
    # Create DataFrame and save
    if integrated_results:
        df = pd.DataFrame(integrated_results)
        df = df.drop_duplicates(subset=['cluster_id', 'protein_id'], keep='first')
        df.to_csv(output_file, index=False)
        logger.info(f"✓ Saved {len(df)} integrated results to {output_file}")
        
        # Print summary
        classification_counts = df['enzyme_classification'].value_counts()
        logger.info(f"Classification summary: {dict(classification_counts)}")
        return True
    else:
        logger.warning(f"No valid results found for cluster {cluster_name}")
        return False

def main():
    """Main function - process a single cluster"""
    if len(sys.argv) != 6:
        print("Usage: python integrate_single_cluster.py <cluster_name> <e2p2_file> <domain_file> <fasta_file> <output_file>")
        sys.exit(1)
    
    cluster_name = sys.argv[1]
    e2p2_file = Path(sys.argv[2])
    domain_file = Path(sys.argv[3])
    fasta_file = Path(sys.argv[4])
    output_file = Path(sys.argv[5])
    
    print(f"🧬 Processing cluster: {cluster_name}")
    print(f"E2P2 file: {e2p2_file}")
    print(f"Domain file: {domain_file}")
    print(f"FASTA file: {fasta_file}")
    print(f"Output file: {output_file}")
    print("=" * 50)
    
    success = integrate_cluster_data(cluster_name, e2p2_file, domain_file, fasta_file, output_file)
    
    if success:
        print("✓ Integration completed successfully")
        sys.exit(0)
    else:
        print("✗ Integration failed")
        sys.exit(1)

if __name__ == "__main__":
    main()