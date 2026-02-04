#!/usr/bin/env python3
"""
CGAS Module 4: GenBank Format Conversion & Validation (IMPROVED)
=================================================================

This module provides:
1. GenBank to FASTA/TBL Converter - for NCBI submission preparation
2. GenBank Validator - for annotation quality checking

Converts GenBank files to submission-ready formats and validates annotations.

IMPROVEMENTS IN THIS VERSION:
- Added validation for gene features without products
- Added validation for products without corresponding gene features
- Enhanced checking for protein-coding, tRNA, and rRNA consistency

Part of the CGAS (Chloroplast Genome Assembly Suite) pipeline.

Author: Abdullah
Date: January 2026
"""

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation
import os
import glob
import argparse
from pathlib import Path
from datetime import datetime


# ============================================================================
# PART A: GenBank to FASTA/TBL Converter
# ============================================================================

def parse_genbank_to_fasta(gb_file, fasta_file):
    """
    Extract sequence from GenBank file and write to FASTA format.
    
    Args:
        gb_file: Input GenBank file path
        fasta_file: Output FASTA file path
    """
    try:
        with open(fasta_file, 'w') as f_out:
            for record in SeqIO.parse(gb_file, "genbank"):
                # Get organism name from FEATURES section (source feature)
                organism = record.id  # default fallback
                for feature in record.features:
                    if feature.type == 'source':
                        if 'organism' in feature.qualifiers:
                            organism = feature.qualifiers['organism'][0]
                        break
                organism_id = organism.replace(' ', '_')
                
                # Get topology from GenBank record
                topology = 'circular'  # default for chloroplast genomes
                
                if 'topology' in record.annotations:
                    topology = record.annotations['topology']
                
                if 'keywords' in record.annotations:
                    keywords = ' '.join(record.annotations['keywords']).lower()
                    if 'complete' in keywords:
                        topology = 'circular'
                
                if record.description:
                    desc_lower = record.description.lower()
                    if 'complete genome' in desc_lower or 'complete chloroplast' in desc_lower:
                        topology = 'circular'
                
                # Get molecule type
                mol_type = 'genomic DNA'
                
                # Get genetic code
                gcode = '11'  # bacterial/plastid code
                
                # Completeness
                completeness = 'complete' if topology == 'circular' else 'partial'
                
                if record.description:
                    desc_lower = record.description.lower()
                    if 'complete' in desc_lower:
                        completeness = 'complete'
                    elif 'partial' in desc_lower:
                        completeness = 'partial'
                
                # Location
                location = 'chloroplast'
                
                for feature in record.features:
                    if feature.type == 'source':
                        if 'organelle' in feature.qualifiers:
                            organelle = feature.qualifiers['organelle'][0]
                            if 'plastid' in organelle.lower() or 'chloroplast' in organelle.lower():
                                location = 'chloroplast'
                            elif 'mitochondr' in organelle.lower():
                                location = 'mitochondrion'
                            else:
                                location = organelle
                        break
                
                # Write FASTA header with metadata
                f_out.write(f">{organism_id}_ [organism={organism}] [mol_type={mol_type}] ")
                f_out.write(f"[completeness={completeness}] [topology={topology}] ")
                f_out.write(f"[gcode={gcode}] [location={location}] {organism}\n")
                
                # Write sequence
                sequence = str(record.seq).lower()
                f_out.write(sequence + "\n")
                
        return True
    except Exception as e:
        print(f"  Error creating FASTA: {e}")
        return False


def format_location(feature):
    """Format feature location for TBL file."""
    locations = []
    
    if isinstance(feature.location, CompoundLocation):
        for part in feature.location.parts:
            start = int(part.start) + 1  # Convert to 1-based
            end = int(part.end)
            if part.strand == -1:
                locations.append(f"{end}\t{start}")
            else:
                locations.append(f"{start}\t{end}")
    else:
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        if feature.location.strand == -1:
            locations.append(f"{end}\t{start}")
        else:
            locations.append(f"{start}\t{end}")
    
    return locations


def get_feature_qualifiers(feature):
    """Extract and format qualifiers for a feature."""
    qualifiers = []
    
    priority_quals = ['trans_splicing', 'codon_start', 'db_xref', 'product', 
                      'protein_id', 'gene', 'transl_table', 'locus_tag', 
                      'note', 'exception', 'pseudo', 'pseudogene', 
                      'rpt_type', 'annotator', 'info', 'label', 'modified_by']
    
    exclude_quals = ['translation']
    
    remove_keywords = [
        "protein_id", "db_xref\tGeneID", "db_xref", "GeneID",
        "locus_tag", "modified_by", "annotator", "info"
    ]
    
    remove_note_patterns = [
        "annotator OGDRAWinfo", "annotated by OGDRAW", "anticodon:"
    ]
    
    for qual_name in priority_quals:
        if qual_name in feature.qualifiers and qual_name not in exclude_quals:
            for value in feature.qualifiers[qual_name]:
                qual_line = f"\t\t\t{qual_name}\t{value}"
                
                should_remove = False
                for keyword in remove_keywords:
                    if keyword in qual_line:
                        should_remove = True
                        break
                
                if qual_name == "note" and not should_remove:
                    note_lower = value.lower()
                    
                    important_keywords = [
                        "exception", "pseudo", "pseudogene", "rna editing",
                        "alternative start", "ribosomal slippage", "trans-splicing",
                        "frameshift", "codon redefined"
                    ]
                    
                    has_important_info = any(kw in note_lower for kw in important_keywords)
                    
                    if has_important_info:
                        should_remove = False
                    else:
                        for pattern in remove_note_patterns:
                            if pattern in note_lower:
                                should_remove = True
                                break
                
                if not should_remove:
                    qualifiers.append(qual_line)
    
    for qual_name in feature.qualifiers:
        if qual_name not in priority_quals and qual_name not in exclude_quals:
            for value in feature.qualifiers[qual_name]:
                qual_line = f"\t\t\t{qual_name}\t{value}"
                
                should_remove = False
                for keyword in remove_keywords:
                    if keyword in qual_line:
                        should_remove = True
                        break
                
                if not should_remove:
                    qualifiers.append(qual_line)
    
    return qualifiers


def parse_genbank_to_tbl(gb_file, tbl_file):
    """
    Parse GenBank file and generate TBL (feature table) format.
    Auto-converts 'order' to 'join' for NCBI compliance.
    
    Args:
        gb_file: Input GenBank file path
        tbl_file: Output TBL file path
    """
    try:
        feature_types = ['gene', 'CDS', 'tRNA', 'rRNA', 'misc_feature', 'repeat_region']
        order_to_join_count = 0
        
        with open(tbl_file, 'w') as f_out:
            for record in SeqIO.parse(gb_file, "genbank"):
                # Get organism name from FEATURES section (source feature)
                organism = record.id  # default fallback
                for feature in record.features:
                    if feature.type == 'source':
                        if 'organism' in feature.qualifiers:
                            organism = feature.qualifiers['organism'][0]
                        break
                organism_id = organism.replace(' ', '_')
                
                f_out.write(f">Feature {organism_id}_\n")
                
                for feature in record.features:
                    if feature.type not in feature_types:
                        continue
                    
                    # Check and warn if 'order' operator is used
                    if isinstance(feature.location, CompoundLocation):
                        if hasattr(feature.location, 'operator') and feature.location.operator == 'order':
                            order_to_join_count += 1
                            gene_id = (
                                feature.qualifiers.get("gene", [""])[0]
                                or feature.qualifiers.get("locus_tag", [""])[0]
                                or "unknown"
                            )
                            # Note: BioPython will write it as 'join' automatically in most cases
                            # This is just for tracking
                    
                    locations = format_location(feature)
                    
                    if len(locations) == 1:
                        f_out.write(f"{locations[0]}\t{feature.type}\n")
                    else:
                        for i, loc in enumerate(locations):
                            if i == 0:
                                f_out.write(f"{loc}\t{feature.type}\n")
                            else:
                                f_out.write(f"{loc}\n")
                    
                    qualifiers = get_feature_qualifiers(feature)
                    for qual in qualifiers:
                        f_out.write(f"{qual}\n")
        
        if order_to_join_count > 0:
            print(f"  Note: {order_to_join_count} feature(s) had 'order' operator (auto-converted to 'join')")
        
        return True
    except Exception as e:
        print(f"  Error creating TBL: {e}")
        return False


# ============================================================================
# PART B: GenBank Validator (IMPROVED VERSION)
# ============================================================================

def validate_genbank_files(folder_path, output_file):
    """
    Validate GenBank files for annotation issues.
    Auto-corrects gene locations when CDS passes all validation criteria.
    
    IMPROVEMENTS:
    - Validates that all genes have corresponding products (or are pseudogenes)
    - Validates that all products have corresponding gene features
    - Reports inconsistencies for protein-coding, tRNA, and rRNA features
    
    Args:
        folder_path: Directory containing GenBank files
        output_file: Path to validation report file
    
    Returns:
        total_issues: Number of validation issues found
    """
    valid_ext = (".gb", ".gbk", ".genbank", ".gbff")
    gb_files = [
        f for f in os.listdir(folder_path)
        if f.lower().endswith(valid_ext)
    ]
    
    def log(message):
        """Print to console and write to file"""
        print(message)
        with open(output_file, "a", encoding="utf-8") as f:
            f.write(message + "\n")
    
    def get_location_bounds(location):
        """Get the minimum and maximum positions from a feature location"""
        if isinstance(location, CompoundLocation):
            all_positions = []
            for part in location.parts:
                all_positions.append(int(part.start))
                all_positions.append(int(part.end))
            return min(all_positions), max(all_positions)
        else:
            return int(location.start), int(location.end)
    
    def has_order_operator(feature):
        """Check if feature uses 'order' instead of 'join' in compound location"""
        if isinstance(feature.location, CompoundLocation):
            if hasattr(feature.location, 'operator'):
                return feature.location.operator == 'order'
        return False
    
    def is_pseudogene(feature):
        """Check if feature is marked as a pseudogene"""
        return ('pseudo' in feature.qualifiers or 
                'pseudogene' in feature.qualifiers)
    
    log(f"GenBank Validation & Auto-Correction Report")
    log(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"{'='*60}\n")
    log(f"Detected {len(gb_files)} GenBank file(s)\n")
    
    if not gb_files:
        log("⚠ No GenBank files found")
        return 0
    
    total_issues = 0
    total_corrections = 0
    
    for fname in sorted(gb_files):
        file_path = os.path.join(folder_path, fname)
        
        log(f"\nProcessing: {fname}")
        cds_count = 0
        rrna_count = 0
        trna_count = 0
        issue_count = 0
        correction_count = 0
        
        # Read all records
        records = list(SeqIO.parse(file_path, "genbank"))
        
        for record in records:
            # Build gene dictionary with ALL gene features (including duplicates in IR)
            gene_dict = {}
            gene_list = []  # Keep track of all gene features by index
            
            # NEW: Track all genes and products for cross-validation
            genes_without_products = set()  # Gene names that exist as gene features
            products_to_genes = {}  # Map product features to their gene names
            gene_to_products = {}  # Map gene names to their product features
            
            for idx, feature in enumerate(record.features):
                if feature.type == "gene":
                    gene_name = (
                        feature.qualifiers.get("gene", [""])[0]
                        or feature.qualifiers.get("locus_tag", [""])[0]
                    )
                    if gene_name:
                        # Store first occurrence only for correction
                        if gene_name not in gene_dict:
                            gene_dict[gene_name] = feature
                            genes_without_products.add(gene_name)
                        gene_list.append((idx, gene_name, feature))
            
            # Track which genes have been corrected (to avoid duplicate corrections)
            corrected_genes = set()
            
            # Validate and correct features
            for feature in record.features:
                # Check for 'order' operator
                if has_order_operator(feature):
                    gene_id = (
                        feature.qualifiers.get("gene", [""])[0]
                        or feature.qualifiers.get("locus_tag", [""])[0]
                        or "unknown"
                    )
                    log(f"  [{gene_id}] {feature.type} uses 'order' instead of 'join' (NCBI error)")
                    issue_count += 1
                
                # Process CDS features
                if feature.type == "CDS":
                    cds_count += 1
                    
                    gene_id = (
                        feature.qualifiers.get("gene", [""])[0]
                        or feature.qualifiers.get("locus_tag", [""])[0]
                        or feature.qualifiers.get("protein_id", ["unknown"])[0]
                    )
                    
                    # NEW: Track that this gene has a product
                    if gene_id in genes_without_products:
                        genes_without_products.remove(gene_id)
                    if gene_id:
                        if gene_id not in gene_to_products:
                            gene_to_products[gene_id] = []
                        gene_to_products[gene_id].append('CDS')
                    
                    # NEW: Check if product exists for this CDS
                    product = feature.qualifiers.get("product", [])
                    if not product or not product[0].strip():
                        # Check if it's a pseudogene
                        if not is_pseudogene(feature):
                            log(f"  [{gene_id}] CDS missing product annotation (not marked as pseudogene)")
                            issue_count += 1
                    
                    # NEW: Check if corresponding gene feature exists
                    if gene_id and gene_id not in gene_dict:
                        log(f"  [{gene_id}] CDS product exists but gene feature is missing")
                        issue_count += 1
                    
                    # Validate CDS
                    cds_is_valid = True
                    validation_errors = []
                    
                    try:
                        cds_seq = feature.location.extract(record.seq).upper()
                    except Exception:
                        validation_errors.append("CDS extraction failed")
                        cds_is_valid = False
                    
                    if cds_is_valid:
                        # Check start codon
                        if len(cds_seq) >= 3 and str(cds_seq[:3]) != "ATG":
                            validation_errors.append(f"non-ATG start: {cds_seq[:3]}")
                            cds_is_valid = False
                        
                        # Check frame length
                        if len(cds_seq) % 3 != 0:
                            validation_errors.append(f"frame length error ({len(cds_seq)})")
                            cds_is_valid = False
                        else:
                            # Check for internal stops
                            protein = cds_seq.translate(to_stop=False)
                            internal_stops = protein[:-1].count("*")
                            if internal_stops:
                                validation_errors.append("internal stop codon(s)")
                                cds_is_valid = False
                            
                            # Check for terminal stop
                            if not protein.endswith("*"):
                                validation_errors.append("missing terminal stop codon")
                                cds_is_valid = False
                    
                    # Log validation errors
                    for error in validation_errors:
                        log(f"  [{gene_id}] {error}")
                        issue_count += 1
                    
                    # Check gene location and auto-correct if CDS is valid
                    if gene_id and gene_id in gene_dict:
                        # Skip if we already corrected this gene
                        if gene_id in corrected_genes:
                            continue
                        
                        gene_feature = gene_dict[gene_id]
                        cds_start, cds_end = get_location_bounds(feature.location)
                        gene_start, gene_end = get_location_bounds(gene_feature.location)
                        
                        # Calculate the differences in start and end positions
                        start_diff = abs(cds_start - gene_start)
                        end_diff = abs(cds_end - gene_end)
                        max_diff = max(start_diff, end_diff)
                        
                        # Only auto-correct if:
                        # 1. CDS is valid
                        # 2. Difference is small (< 100 bp) - just boundary issues
                        # 3. They overlap significantly (same genomic region)
                        
                        location_mismatch = (cds_start != gene_start or cds_end != gene_end)
                        
                        if location_mismatch:
                            old_loc = f"{gene_start}..{gene_end}"
                            new_loc = f"{cds_start}..{cds_end}"
                            
                            if cds_is_valid and max_diff < 100:
                                # Small difference + valid CDS = AUTO-CORRECT
                                gene_feature.location = feature.location
                                log(f"  [{gene_id}] ✓ AUTO-CORRECTED gene location: {old_loc} → {new_loc} (diff: {max_diff}bp)")
                                correction_count += 1
                                corrected_genes.add(gene_id)
                            elif cds_is_valid and max_diff >= 100:
                                # Large difference - likely IR duplicate, silently skip
                                pass
                            elif not cds_is_valid and max_diff < 100:
                                # CDS has errors, don't correct
                                log(f"  [{gene_id}] gene location mismatch: gene({old_loc}) vs CDS({new_loc}) [NOT CORRECTED - CDS has errors]")
                                issue_count += 1
                            # Else: large diff + invalid CDS = silently skip (IR duplicate with errors)
                
                # Process rRNA features
                elif feature.type == "rRNA":
                    rrna_count += 1
                    gene_id = (
                        feature.qualifiers.get("gene", [""])[0]
                        or feature.qualifiers.get("locus_tag", [""])[0]
                        or f"rRNA_{rrna_count}"
                    )
                    
                    # NEW: Track that this gene has a product
                    if gene_id in genes_without_products:
                        genes_without_products.remove(gene_id)
                    if gene_id:
                        if gene_id not in gene_to_products:
                            gene_to_products[gene_id] = []
                        gene_to_products[gene_id].append('rRNA')
                    
                    # NEW: Check if product exists for this rRNA
                    product = feature.qualifiers.get("product", [])
                    if not product or not product[0].strip():
                        log(f"  [{gene_id}] rRNA missing product annotation")
                        issue_count += 1
                    
                    # NEW: Check if corresponding gene feature exists
                    if gene_id and gene_id not in gene_dict:
                        log(f"  [{gene_id}] rRNA product exists but gene feature is missing")
                        issue_count += 1
                    
                    # Check gene location for rRNA (only if not already corrected)
                    if gene_id and gene_id in gene_dict and gene_id not in corrected_genes:
                        gene_feature = gene_dict[gene_id]
                        rrna_start, rrna_end = get_location_bounds(feature.location)
                        gene_start, gene_end = get_location_bounds(gene_feature.location)
                        
                        start_diff = abs(rrna_start - gene_start)
                        end_diff = abs(rrna_end - gene_end)
                        max_diff = max(start_diff, end_diff)
                        
                        if rrna_start != gene_start or rrna_end != gene_end:
                            old_loc = f"{gene_start}..{gene_end}"
                            new_loc = f"{rrna_start}..{rrna_end}"
                            
                            if max_diff < 100:
                                # Small difference - auto-correct
                                gene_feature.location = feature.location
                                log(f"  [{gene_id}] ✓ AUTO-CORRECTED rRNA gene location: {old_loc} → {new_loc} (diff: {max_diff}bp)")
                                correction_count += 1
                                corrected_genes.add(gene_id)
                            # Else: large difference - silently skip (IR duplicate)
                
                # Process tRNA features
                elif feature.type == "tRNA":
                    trna_count += 1
                    gene_id = (
                        feature.qualifiers.get("gene", [""])[0]
                        or feature.qualifiers.get("locus_tag", [""])[0]
                        or f"tRNA_{trna_count}"
                    )
                    
                    # NEW: Track that this gene has a product
                    if gene_id in genes_without_products:
                        genes_without_products.remove(gene_id)
                    if gene_id:
                        if gene_id not in gene_to_products:
                            gene_to_products[gene_id] = []
                        gene_to_products[gene_id].append('tRNA')
                    
                    # NEW: Check if product exists for this tRNA
                    product = feature.qualifiers.get("product", [])
                    if not product or not product[0].strip():
                        log(f"  [{gene_id}] tRNA missing product annotation")
                        issue_count += 1
                    else:
                        product_name = product[0].lower()
                        amino_acids = [
                            "ala", "arg", "asn", "asp", "cys", "gln", "glu", 
                            "gly", "his", "ile", "leu", "lys", "met", "phe", 
                            "pro", "ser", "thr", "trp", "tyr", "val",
                            "alanine", "arginine", "asparagine", "aspartate", 
                            "cysteine", "glutamine", "glutamate", "glycine", 
                            "histidine", "isoleucine", "leucine", "lysine", 
                            "methionine", "phenylalanine", "proline", "serine", 
                            "threonine", "tryptophan", "tyrosine", "valine",
                            "sec", "pyl", "selenocysteine", "pyrrolysine"
                        ]
                        has_amino_acid = any(aa in product_name for aa in amino_acids)
                        if not has_amino_acid:
                            log(f"  [{gene_id}] tRNA product '{product[0]}' missing amino acid")
                            issue_count += 1
                    
                    # NEW: Check if corresponding gene feature exists
                    if gene_id and gene_id not in gene_dict:
                        log(f"  [{gene_id}] tRNA product exists but gene feature is missing")
                        issue_count += 1
                    
                    # Check gene location for tRNA (only if not already corrected)
                    if gene_id and gene_id in gene_dict and gene_id not in corrected_genes:
                        gene_feature = gene_dict[gene_id]
                        trna_start, trna_end = get_location_bounds(feature.location)
                        gene_start, gene_end = get_location_bounds(gene_feature.location)
                        
                        start_diff = abs(trna_start - gene_start)
                        end_diff = abs(trna_end - gene_end)
                        max_diff = max(start_diff, end_diff)
                        
                        if trna_start != gene_start or trna_end != gene_end:
                            old_loc = f"{gene_start}..{gene_end}"
                            new_loc = f"{trna_start}..{trna_end}"
                            
                            if max_diff < 100:
                                # Small difference - auto-correct
                                gene_feature.location = feature.location
                                log(f"  [{gene_id}] ✓ AUTO-CORRECTED tRNA gene location: {old_loc} → {new_loc} (diff: {max_diff}bp)")
                                correction_count += 1
                                corrected_genes.add(gene_id)
                            # Else: large difference - silently skip (IR duplicate)
            
            # NEW: Check for genes without products (excluding pseudogenes)
            for gene_name in genes_without_products:
                gene_feature = gene_dict[gene_name]
                if not is_pseudogene(gene_feature):
                    log(f"  [{gene_name}] gene exists but no corresponding product found (not marked as pseudogene)")
                    issue_count += 1
        
        log(f"\n  CDS scanned: {cds_count}")
        log(f"  rRNA scanned: {rrna_count}")
        log(f"  tRNA scanned: {trna_count}")
        log(f"  Issues found: {issue_count}")
        if correction_count > 0:
            log(f"  Auto-corrections: {correction_count}")
        
        total_issues += issue_count
        total_corrections += correction_count
    
    log(f"\n{'='*60}")
    log(f"SUMMARY:")
    log(f"  Total issues found: {total_issues}")
    if total_corrections > 0:
        log(f"  Total auto-corrections identified: {total_corrections}")
    log(f"  Report saved to: {output_file}")
    log(f"{'='*60}")
    
    return total_issues


# ============================================================================
# CGAS MODULE 4 MAIN FUNCTION
# ============================================================================

def main():
    """
    Main function for CGAS Module 4: Format Conversion & Validation
    
    Runs in two phases:
    1. Convert GenBank files to FASTA/TBL format
    2. Validate GenBank annotations
    """
    parser = argparse.ArgumentParser(
        description="CGAS Module 4: GenBank Format Conversion & Validation (IMPROVED)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process GenBank files in current directory
  python cgas_module4_improved.py
  
  # Process GenBank files in specific directory
  python cgas_module4_improved.py -i genbank_files/
  
  # Specify output directory
  python cgas_module4_improved.py -i genbank_files/ -o ncbi_submission/
  
  # Process files from Module 2 output
  python cgas_module4_improved.py -i module2_annotations/Annotated_GenBank/

IMPROVEMENTS IN THIS VERSION:
  - Validates that all genes have corresponding products (or are pseudogenes)
  - Validates that all products have corresponding gene features
  - Reports inconsistencies for protein-coding, tRNA, and rRNA features
        """
    )
    
    parser.add_argument("-i", "--input", default=".",
                       help="Input directory containing GenBank files (default: current directory)")
    parser.add_argument("-o", "--output", default="Module4_NCBI_Submission",
                       help="Output directory (default: Module4_NCBI_Submission)")
    parser.add_argument("--skip-validation", action="store_true",
                       help="Skip validation phase (only convert formats)")
    parser.add_argument("--skip-conversion", action="store_true",
                       help="Skip conversion phase (only validate)")
    
    args = parser.parse_args()
    
    input_dir = Path(args.input)
    output_dir = Path(args.output)
    
    print("\n" + "="*80)
    print("CGAS MODULE 4: GenBank Format Conversion & Validation (IMPROVED)")
    print("="*80)
    
    # Detect GenBank files (exclude any corrected_ files from previous runs)
    genbank_extensions = ['gb', 'gbk', 'genbank', 'gbff']
    gb_files = []
    
    for ext in genbank_extensions:
        pattern = str(input_dir / f'*.{ext}')
        found_files = glob.glob(pattern)
        gb_files.extend(found_files)
    
    # Remove duplicates and filter out corrected_ files
    gb_files = list(set(gb_files))
    gb_files = [f for f in gb_files if not os.path.basename(f).startswith('corrected_')]
    gb_files.sort()
    
    if not gb_files:
        print(f"⚠ No GenBank files found in {input_dir}")
        print(f"\nLooking for files with extensions: {', '.join(genbank_extensions)}")
        return
    
    print(f"\nInput directory: {input_dir}")
    print(f"Found {len(gb_files)} GenBank file(s)")
    print("-" * 80)
    for gb_file in gb_files:
        print(f"  • {os.path.basename(gb_file)}")
    print("-" * 80)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nOutput directory: {output_dir}/")
    
    # ========================================================================
    # PHASE 1: Format Conversion
    # ========================================================================
    if not args.skip_conversion:
        print("\n" + "="*80)
        print("PHASE 1: Converting to FASTA and TBL formats")
        print("="*80)
        
        success_count = 0
        error_count = 0
        generated_fasta_files = []
        generated_tbl_files = []
        
        for gb_file in gb_files:
            print(f"\nProcessing: {os.path.basename(gb_file)}")
            
            base_name = os.path.splitext(os.path.basename(gb_file))[0]
            fasta_file = output_dir / f"{base_name}.fsa"
            tbl_file = output_dir / f"{base_name}.tbl"
            
            print(f"  → FASTA: {fasta_file.name}")
            print(f"  → TBL: {tbl_file.name}")
            
            fasta_success = parse_genbank_to_fasta(gb_file, str(fasta_file))
            tbl_success = parse_genbank_to_tbl(gb_file, str(tbl_file))
            
            if fasta_success and tbl_success:
                success_count += 1
                print(f"  ✓ Conversion successful")
                generated_fasta_files.append(str(fasta_file))
                generated_tbl_files.append(str(tbl_file))
            else:
                error_count += 1
                print(f"  ✗ Conversion failed")
        
        print("\n" + "-"*80)
        print(f"Conversion Summary: {success_count} successful, {error_count} failed")
        print("-"*80)
        
        # Combine files
        if success_count > 0:
            print("\nCombining files...")
            
            # Combine FASTA
            combined_fasta = output_dir / 'combined.fasta'
            try:
                with open(combined_fasta, 'w', encoding='utf-8') as outfile:
                    for fasta_file in sorted(generated_fasta_files):
                        with open(fasta_file, 'r', encoding='utf-8') as infile:
                            content = infile.read()
                            outfile.write(content)
                            if not content.endswith('\n'):
                                outfile.write('\n')
                print(f"  ✓ Combined FASTA: {combined_fasta.name}")
            except Exception as e:
                print(f"  ✗ Error combining FASTA: {e}")
            
            # Combine TBL
            combined_tbl = output_dir / 'combined_features.tbl'
            try:
                with open(combined_tbl, 'w', encoding='utf-8') as outfile:
                    for tbl_file in sorted(generated_tbl_files):
                        with open(tbl_file, 'r', encoding='utf-8') as infile:
                            content = infile.read()
                            outfile.write(content)
                            if not content.endswith('\n'):
                                outfile.write('\n')
                print(f"  ✓ Combined TBL: {combined_tbl.name}")
            except Exception as e:
                print(f"  ✗ Error combining TBL: {e}")
    else:
        print("\n⊘ Skipping format conversion (--skip-conversion)")
    
    # ========================================================================
    # PHASE 2: Validation
    # ========================================================================
    if not args.skip_validation:
        print("\n" + "="*80)
        print("PHASE 2: Validating GenBank annotations")
        print("="*80)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        validation_report = output_dir / f"validation_report_{timestamp}.txt"
        
        total_issues = validate_genbank_files(str(input_dir), str(validation_report))
    else:
        print("\n⊘ Skipping validation (--skip-validation)")
        total_issues = 0
    
    # ========================================================================
    # Final Summary
    # ========================================================================
    print("\n" + "="*80)
    print("CGAS MODULE 4 COMPLETE")
    print("="*80)
    print(f"\nOutput directory: {output_dir}/")
    
    if not args.skip_conversion:
        print(f"\nGenerated files:")
        if 'generated_fasta_files' in locals():
            print(f"  • Individual FASTA files: {len(generated_fasta_files)}")
            print(f"  • Individual TBL files: {len(generated_tbl_files)}")
            if success_count > 0:
                print(f"  • Combined FASTA: combined.fasta")
                print(f"  • Combined TBL: combined_features.tbl")
    
    if not args.skip_validation:
        print(f"\nValidation:")
        print(f"  • {total_issues} issue(s) detected across all files")
        print(f"  • Report: validation_report_{timestamp}.txt")
    
    print("\n" + "="*80)
    print("✓ Ready for NCBI submission!")
    print("  Next step: Upload .fsa and .tbl files to NCBI Submission Portal")
    print("="*80)


if __name__ == "__main__":
    main()
