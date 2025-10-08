#!/usr/bin/env python3
"""
Database-Driven Enzyme Classification via Reference Proteins
===========================================================

Classification strategy:
1. Use curated reference proteins (polyketide synthases, methyltransferases, etc.)
2. Extract their GO term profiles
3. Compare query enzymes to reference profiles statistically
4. NO keyword matching - pure profile similarity

Author: Scientific Analysis System
"""

import requests
import json
import time
import logging
from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass
from enum import Enum
from collections import Counter
import math

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ClassificationResult(Enum):
    SIGNATURE = "signature"
    TAILORING = "tailoring"
    MIXED = "mixed"
    UNKNOWN = "unknown"

@dataclass
class Evidence:
    source: str
    classification: str
    confidence: float
    details: str

@dataclass
class ClassificationInfo:
    result: ClassificationResult
    confidence: float
    evidence: List[Evidence]
    reasoning: str

class DatabaseDrivenClassifier:
    """
    Classifier using reference protein GO profiles
    """
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({'User-Agent': 'DatabaseDrivenClassifier/1.0'})
        
        self.uniprot_base = "https://rest.uniprot.org"
        self.quickgo_base = "https://www.ebi.ac.uk/QuickGO/services"
        
        # Reference EC numbers with known classifications (from scientific literature)
        self.biosynthetic_references = {
            '2.3.1.85',   # Fatty acid synthase
            '2.3.1.41',   # Malonate CoA-transferase  
            '6.2.1.3',    # Long-chain fatty acid-CoA ligase
            '4.2.1.17',   # Enoyl-CoA hydratase
            '2.3.1.16',   # Acetyl-CoA acetyltransferase
        }
        
        self.modification_references = {
            '2.1.1.37',   # DNA methyltransferase
            '2.1.1.43',   # Histone-lysine N-methyltransferase
            '1.14.13.39', # Cytochrome P450 hydroxylase
            '2.4.1.17',   # Glucosyltransferase
            '2.3.1.48',   # Acetyltransferase
        }
        
        # Build reference GO profiles
        logger.info("Building reference GO profiles from curated proteins...")
        self.biosynthetic_go_profile = self._build_reference_profile(self.biosynthetic_references)
        self.modification_go_profile = self._build_reference_profile(self.modification_references)
        logger.info(f"Biosynthetic profile: {len(self.biosynthetic_go_profile)} GO terms")
        logger.info(f"Modification profile: {len(self.modification_go_profile)} GO terms")
    
    def _build_reference_profile(self, ec_numbers: Set[str]) -> Counter:
        """Build GO term frequency profile from reference EC numbers"""
        go_profile = Counter()
        
        for ec in ec_numbers:
            try:
                # Get proteins for this EC
                url = f"{self.uniprot_base}/uniprotkb/search"
                params = {'query': f'ec:{ec}', 'format': 'json', 'size': 20}
                
                response = self.session.get(url, params=params, timeout=10)
                if response.status_code != 200:
                    continue
                
                data = response.json()
                results = data.get('results', [])
                
                # Extract GO terms
                for protein in results:
                    go_refs = protein.get('uniProtKBCrossReferences', [])
                    for xref in go_refs:
                        if xref.get('database') == 'GO':
                            go_id = xref.get('id', '')
                            if go_id:
                                go_profile[go_id] += 1
                
                time.sleep(0.5)  # Rate limiting
                
            except Exception as e:
                logger.warning(f"Failed to process reference EC {ec}: {e}")
        
        return go_profile
    
    def classify_enzyme(self, ec_number: str, domains: List[str]) -> ClassificationInfo:
        """Classify enzyme by comparing to reference profiles"""
        logger.info(f"Classifying: EC={ec_number}, Domains={domains}")
        
        evidence = []
        
        # 1. GO profile similarity
        go_evidence = self._classify_by_go_similarity(ec_number)
        if go_evidence:
            evidence.append(go_evidence)
        
        # 2. Domain profile similarity
        domain_evidence = self._classify_by_domain_similarity(domains)
        evidence.extend(domain_evidence)
        
        # 3. EC number structural similarity
        ec_evidence = self._classify_by_ec_structure(ec_number)
        if ec_evidence:
            evidence.append(ec_evidence)
        
        # Statistical consensus
        return self._calculate_consensus(evidence)
    
    def _classify_by_go_similarity(self, ec_number: str) -> Optional[Evidence]:
        """Classify by GO term profile similarity to references"""
        try:
            # Get GO profile for query EC
            url = f"{self.uniprot_base}/uniprotkb/search"
            params = {'query': f'ec:{ec_number}', 'format': 'json', 'size': 30}
            
            response = self.session.get(url, params=params, timeout=15)
            if response.status_code != 200:
                return None
            
            data = response.json()
            results = data.get('results', [])
            
            if not results:
                return None
            
            # Build query GO profile
            query_profile = Counter()
            for protein in results:
                go_refs = protein.get('uniProtKBCrossReferences', [])
                for xref in go_refs:
                    if xref.get('database') == 'GO':
                        go_id = xref.get('id', '')
                        if go_id:
                            query_profile[go_id] += 1
            
            if not query_profile:
                return None
            
            # Calculate similarity to each reference profile
            bio_similarity = self._cosine_similarity(query_profile, self.biosynthetic_go_profile)
            mod_similarity = self._cosine_similarity(query_profile, self.modification_go_profile)
            
            logger.info(f"GO similarity - Biosynthetic: {bio_similarity:.3f}, Modification: {mod_similarity:.3f}")
            
            # Require minimum similarity threshold
            if max(bio_similarity, mod_similarity) < 0.1:
                return None
            
            if bio_similarity > mod_similarity * 1.2:  # 20% margin
                classification = "signature"
                confidence = min(0.9, bio_similarity * 1.5)
            elif mod_similarity > bio_similarity * 1.2:
                classification = "tailoring"
                confidence = min(0.9, mod_similarity * 1.5)
            else:
                classification = "mixed"
                confidence = 0.5
            
            return Evidence(
                source="go_profile",
                classification=classification,
                confidence=confidence,
                details=f"GO similarity: bio={bio_similarity:.2f}, mod={mod_similarity:.2f}"
            )
        
        except Exception as e:
            logger.warning(f"GO similarity analysis failed: {e}")
            return None
    
    def _cosine_similarity(self, profile1: Counter, profile2: Counter) -> float:
        """Calculate cosine similarity between two GO profiles"""
        # Get common terms
        common = set(profile1.keys()) & set(profile2.keys())
        
        if not common:
            return 0.0
        
        # Calculate dot product
        dot_product = sum(profile1[term] * profile2[term] for term in common)
        
        # Calculate magnitudes
        mag1 = math.sqrt(sum(count**2 for count in profile1.values()))
        mag2 = math.sqrt(sum(count**2 for count in profile2.values()))
        
        if mag1 == 0 or mag2 == 0:
            return 0.0
        
        return dot_product / (mag1 * mag2)
    
    def _classify_by_domain_similarity(self, domains: List[str]) -> List[Evidence]:
        """Classify by comparing domain annotations to reference proteins"""
        evidence = []
        
        for domain in domains:
            try:
                # Search for proteins with this domain
                url = f"{self.uniprot_base}/uniprotkb/search"
                params = {'query': f'family:{domain}', 'format': 'json', 'size': 50}
                
                response = self.session.get(url, params=params, timeout=10)
                if response.status_code != 200:
                    continue
                
                data = response.json()
                results = data.get('results', [])
                
                if not results:
                    continue
                
                # Build domain GO profile
                domain_profile = Counter()
                for protein in results:
                    go_refs = protein.get('uniProtKBCrossReferences', [])
                    for xref in go_refs:
                        if xref.get('database') == 'GO':
                            go_id = xref.get('id', '')
                            if go_id:
                                domain_profile[go_id] += 1
                
                if not domain_profile:
                    continue
                
                # Compare to references
                bio_sim = self._cosine_similarity(domain_profile, self.biosynthetic_go_profile)
                mod_sim = self._cosine_similarity(domain_profile, self.modification_go_profile)
                
                if max(bio_sim, mod_sim) < 0.05:
                    continue
                
                if bio_sim > mod_sim * 1.2:
                    classification = "signature"
                    confidence = min(0.8, bio_sim * 1.5)
                elif mod_sim > bio_sim * 1.2:
                    classification = "tailoring"
                    confidence = min(0.8, mod_sim * 1.5)
                else:
                    classification = "mixed"
                    confidence = 0.4
                
                evidence.append(Evidence(
                    source=f"domain_{domain}",
                    classification=classification,
                    confidence=confidence,
                    details=f"Domain profile: bio={bio_sim:.2f}, mod={mod_sim:.2f}"
                ))
                
            except Exception as e:
                logger.debug(f"Domain analysis failed for {domain}: {e}")
        
        return evidence
    
    def _classify_by_ec_structure(self, ec_number: str) -> Optional[Evidence]:
        """Classify by EC number pattern similarity to references"""
        try:
            parts = ec_number.split('.')
            if len(parts) != 4:
                return None
            
            main_class = parts[0]
            subclass = f"{parts[0]}.{parts[1]}"
            
            # Count how many reference ECs share main class/subclass
            bio_main_matches = sum(1 for ec in self.biosynthetic_references if ec.startswith(f"{main_class}."))
            mod_main_matches = sum(1 for ec in self.modification_references if ec.startswith(f"{main_class}."))
            
            bio_sub_matches = sum(1 for ec in self.biosynthetic_references if ec.startswith(f"{subclass}."))
            mod_sub_matches = sum(1 for ec in self.modification_references if ec.startswith(f"{subclass}."))
            
            # Weight subclass matches more heavily
            bio_score = bio_main_matches * 0.5 + bio_sub_matches * 2.0
            mod_score = mod_main_matches * 0.5 + mod_sub_matches * 2.0
            
            total = bio_score + mod_score
            if total < 0.5:
                return None
            
            bio_ratio = bio_score / total
            mod_ratio = mod_score / total
            
            if bio_ratio > 0.65:
                classification = "signature"
                confidence = bio_ratio * 0.6
            elif mod_ratio > 0.65:
                classification = "tailoring"
                confidence = mod_ratio * 0.6
            else:
                classification = "mixed"
                confidence = 0.3
            
            return Evidence(
                source="ec_structure",
                classification=classification,
                confidence=confidence,
                details=f"EC pattern: bio={bio_score:.1f}, mod={mod_score:.1f}"
            )
        
        except Exception as e:
            logger.debug(f"EC structure analysis failed: {e}")
            return None
    
    def _calculate_consensus(self, evidence: List[Evidence]) -> ClassificationInfo:
        """Calculate weighted consensus from evidence"""
        if not evidence:
            return ClassificationInfo(
                result=ClassificationResult.UNKNOWN,
                confidence=0.0,
                evidence=[],
                reasoning="No evidence from any source"
            )
        
        # Weight evidence by confidence
        sig_score = sum(e.confidence for e in evidence if e.classification == "signature")
        tail_score = sum(e.confidence for e in evidence if e.classification == "tailoring")
        mix_score = sum(e.confidence for e in evidence if e.classification == "mixed")
        
        total = sig_score + tail_score + mix_score
        
        if total == 0:
            return ClassificationInfo(
                result=ClassificationResult.UNKNOWN,
                confidence=0.0,
                evidence=evidence,
                reasoning="All evidence inconclusive"
            )
        
        sig_ratio = sig_score / total
        tail_ratio = tail_score / total
        mix_ratio = mix_score / total
        
        # Decision logic
        if sig_ratio > 0.6:
            result = ClassificationResult.SIGNATURE
            confidence = sig_ratio
            reasoning = f"Strong signature similarity ({sig_ratio:.1%} weighted)"
        elif tail_ratio > 0.6:
            result = ClassificationResult.TAILORING
            confidence = tail_ratio
            reasoning = f"Strong tailoring similarity ({tail_ratio:.1%} weighted)"
        elif sig_ratio > tail_ratio and sig_ratio > 0.4:
            result = ClassificationResult.SIGNATURE
            confidence = sig_ratio * 0.85
            reasoning = f"Moderate signature similarity ({sig_ratio:.1%} vs {tail_ratio:.1%})"
        elif tail_ratio > sig_ratio and tail_ratio > 0.4:
            result = ClassificationResult.TAILORING
            confidence = tail_ratio * 0.85
            reasoning = f"Moderate tailoring similarity ({tail_ratio:.1%} vs {sig_ratio:.1%})"
        else:
            result = ClassificationResult.MIXED
            confidence = max(sig_ratio, tail_ratio, mix_ratio)
            reasoning = f"Mixed: sig={sig_ratio:.1%}, tail={tail_ratio:.1%}, mix={mix_ratio:.1%}"
        
        return ClassificationInfo(
            result=result,
            confidence=confidence,
            evidence=evidence,
            reasoning=reasoning
        )

def test_classifier():
    """Test the classifier"""
    classifier = DatabaseDrivenClassifier()
    
    test_cases = [
        ("2.3.1.199", ["Acyl_transf_3"]),
        ("2.1.1.104", ["Methyltransf_3"]),
        ("3.1.1.31", ["Abhydrolase_1"]),
    ]
    
    for ec, domains in test_cases:
        print(f"\n{'='*60}")
        print(f"Testing: EC {ec}, Domains: {domains}")
        print('='*60)
        
        result = classifier.classify_enzyme(ec, domains)
        
        print(f"\nResult: {result.result.value.upper()}")
        print(f"Confidence: {result.confidence:.2%}")
        print(f"Reasoning: {result.reasoning}")
        print(f"\nEvidence ({len(result.evidence)} sources):")
        for i, ev in enumerate(result.evidence, 1):
            print(f"  {i}. {ev.source}: {ev.classification} ({ev.confidence:.2%})")
            print(f"     {ev.details}")

if __name__ == "__main__":
    test_classifier()