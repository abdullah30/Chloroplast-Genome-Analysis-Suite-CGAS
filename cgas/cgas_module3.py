#!/usr/bin/env python3
"""
CGAS Module 3: Plastome Gene Comparison and Normalization
==========================================================

This module compares target GenBank annotations against a reference and:
1. Normalizes gene names to match reference nomenclature
2. Detects missing or extra genes
3. Validates product descriptions
4. Checks for intron presence/absence
5. Compares CDS lengths
6. Generates comprehensive comparison reports

Part of the CGAS (Chloroplast Genome Assembly Suite) pipeline.

Author: Abdullah
Date: January 2026
"""

import os
import re
import argparse
from copy import deepcopy
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation


# ===============================
# Paths
# ===============================

OUTPUT_DIR = "module_3"
REVISED_DIR = os.path.join(OUTPUT_DIR, "revise_annotations")
EXCEL_OUT = os.path.join(OUTPUT_DIR, "comparison_results.xlsx")


# ===============================
# Utilities
# ===============================

def ensure_dirs():
    os.makedirs(REVISED_DIR, exist_ok=True)


def is_genbank(fname):
    return fname.endswith((".gb", ".gbk"))


def clean_name(name):
    return re.sub(r"\s+", "", name.lower())


def cds_length(feature):
    if isinstance(feature.location, CompoundLocation):
        return sum(len(p) for p in feature.location.parts)
    return len(feature.location)


def has_intron(feature):
    return isinstance(feature.location, CompoundLocation)


# ===============================
# Amino acid name mapping
# ===============================

AA_THREE_TO_ONE = {
    'ala': 'A', 'arg': 'R', 'asn': 'N', 'asp': 'D',
    'cys': 'C', 'gln': 'Q', 'glu': 'E', 'gly': 'G',
    'his': 'H', 'ile': 'I', 'leu': 'L', 'lys': 'K',
    'met': 'M', 'phe': 'F', 'pro': 'P', 'ser': 'S',
    'thr': 'T', 'trp': 'W', 'tyr': 'Y', 'val': 'V',
    'sec': 'U', 'pyl': 'O', 'fmet': 'M'
}


# ===============================
# Reference parsing
# ===============================

def parse_reference(ref_gb):
    record = SeqIO.read(ref_gb, "genbank")

    ref = {
        "genes": {},      # gene_name -> feature
        "trna": {},       # gene_name -> feature
        "rrna": {},       # gene_name -> feature
        "intron_genes": set(),
        "cds_lengths": {}
    }

    for f in record.features:
        if "gene" not in f.qualifiers:
            continue

        gene = f.qualifiers["gene"][0]

        if f.type == "CDS":
            ref["genes"][gene] = f
            ref["cds_lengths"][gene] = cds_length(f)
            if has_intron(f) and gene != "rps12":
                ref["intron_genes"].add(gene)

        elif f.type == "tRNA":
            ref["trna"][gene] = f

        elif f.type == "rRNA":
            ref["rrna"][gene] = f

    return ref


# ===============================
# Matching logic
# ===============================

def match_gene(raw, reference_genes):
    """Match a gene name to reference and return the reference gene name."""
    raw = raw.strip()
    
    # Remove parentheses and their content
    raw = re.sub(r"\([^)]*\)", "", raw)
    
    # Remove non-alphanumeric (except - and _)
    raw = re.sub(r"[^A-Za-z0-9_-]", "", raw)
    
    # Remove trailing/leading separators
    raw = raw.strip("-_")
    
    raw_clean = clean_name(raw)
    
    # Try exact match (case-insensitive)
    for ref in reference_genes:
        if clean_name(ref) == raw_clean:
            return ref
    
    # Try without trailing digit (clpP1 -> clpP)
    raw_no_digit = raw_clean.rstrip("0123456789")
    for ref in reference_genes:
        if clean_name(ref) == raw_no_digit:
            return ref
    
    # Try prefix match for genes like rps16 vs rps16-exon1
    for ref in reference_genes:
        ref_clean = clean_name(ref)
        if raw_clean.startswith(ref_clean) or ref_clean.startswith(raw_clean):
            return ref
    
    return None


def match_trna(raw, reference_trnas):
    """Match a tRNA name to reference and return the reference gene name."""
    raw = raw.strip()
    
    # Remove parentheses and their content first
    raw = re.sub(r"\([^)]*\)", "", raw)
    
    # Normalize separators
    raw = raw.replace("_", "-").replace(" ", "")
    
    # Convert to lowercase for processing
    raw_lower = raw.lower()
    
    # Try to extract amino acid letter
    aa_letter = None
    anticodon = None
    
    # Pattern 1: trn[A](-)(anticodon)
    m = re.search(r'trn[-]?([a-z])(?:-([acgtu]{3}))?', raw_lower)
    if m:
        aa_letter = m.group(1).upper()
        if m.group(2):
            anticodon = m.group(2).upper()
    
    # Pattern 2: tRNA with full amino acid name
    if not aa_letter:
        m = re.search(r'trna[-\s]?([a-z]+)', raw_lower)
        if m:
            aa_name = m.group(1)
            # Try three-letter code
            if aa_name in AA_THREE_TO_ONE:
                aa_letter = AA_THREE_TO_ONE[aa_name]
            # Try full name
            else:
                for three, one in AA_THREE_TO_ONE.items():
                    if aa_name.startswith(three):
                        aa_letter = one
                        break
    
    if not aa_letter:
        return None
    
    # Match to reference based on amino acid letter
    candidates = [ref for ref in reference_trnas if ref.lower().startswith(f'trn{aa_letter.lower()}')]
    
    if not candidates:
        return None
    
    # If we have anticodon info, try to match it
    if anticodon:
        for ref in candidates:
            if anticodon.lower() in ref.lower():
                return ref
    
    # If only one candidate, use it
    if len(candidates) == 1:
        return candidates[0]
    
    # Multiple candidates - return first as best guess
    return candidates[0]


def match_rrna(raw, reference_rrnas):
    """Match an rRNA name to reference and return the reference gene name."""
    raw = raw.strip()
    
    # Remove parentheses and content
    raw = re.sub(r"\([^)]*\)", "", raw)
    
    # Normalize separators and spaces
    raw = re.sub(r"[-_\s]+", "", raw)
    
    raw_clean = clean_name(raw)
    
    # Try exact match
    for ref in reference_rrnas:
        if clean_name(ref) == raw_clean:
            return ref
    
    # Try partial match (rrn16 matches rrn16S, etc.)
    for ref in reference_rrnas:
        ref_clean = clean_name(ref)
        if raw_clean in ref_clean or ref_clean in raw_clean:
            # Check if it's a meaningful match
            if len(raw_clean) >= 4:  # Avoid matching just "rrn"
                return ref
    
    return None


# ===============================
# Target processing
# ===============================

def process_target(gb, reference):
    record = SeqIO.read(gb, "genbank")
    revised = deepcopy(record)

    observed = {
        "genes": {},
        "cds": {},
        "trna": {},
        "rrna": {},
        "issues": []
    }
    
    changes_made = 0

    for f in revised.features:
        if "gene" not in f.qualifiers:
            continue

        raw = f.qualifiers["gene"][0]
        matched_ref_name = None

        if f.type == "gene":
            # Standalone gene feature - try to match it to any reference gene
            # First try CDS
            matched_ref_name = match_gene(raw, reference["genes"].keys())
            if not matched_ref_name:
                # Try tRNA
                matched_ref_name = match_trna(raw, reference["trna"].keys())
            if not matched_ref_name:
                # Try rRNA
                matched_ref_name = match_rrna(raw, reference["rrna"].keys())
            
            if matched_ref_name:
                if raw != matched_ref_name:
                    print(f"  gene: '{raw}' → '{matched_ref_name}'")
                    changes_made += 1
                f.qualifiers["gene"] = [matched_ref_name]
            else:
                print(f"  gene: '{raw}' - NO MATCH FOUND")
                observed["issues"].append({
                    "type": "gene",
                    "original": raw,
                    "issue": f"Unmatched gene feature: '{raw}'"
                })

        elif f.type == "CDS":
            matched_ref_name = match_gene(raw, reference["genes"].keys())
            if matched_ref_name:
                if raw != matched_ref_name:
                    print(f"  CDS: '{raw}' → '{matched_ref_name}'")
                    changes_made += 1
                # Update gene name to reference
                f.qualifiers["gene"] = [matched_ref_name]
                # Copy reference product if target is missing it
                if "product" not in f.qualifiers and "product" in reference["genes"][matched_ref_name].qualifiers:
                    f.qualifiers["product"] = reference["genes"][matched_ref_name].qualifiers["product"]
                    print(f"    Added product: {f.qualifiers['product'][0]}")
                observed["genes"][matched_ref_name] = f
                observed["cds"][matched_ref_name] = f
            else:
                observed["issues"].append({
                    "type": "CDS",
                    "original": raw,
                    "issue": f"Unmatched gene name: '{raw}'"
                })
                print(f"  CDS: '{raw}' - NO MATCH FOUND")

        elif f.type == "tRNA":
            matched_ref_name = match_trna(raw, reference["trna"].keys())
            if matched_ref_name:
                if raw != matched_ref_name:
                    print(f"  tRNA: '{raw}' → '{matched_ref_name}'")
                    changes_made += 1
                # Update gene name to reference
                f.qualifiers["gene"] = [matched_ref_name]
                # Copy reference product if target is missing it
                if "product" not in f.qualifiers and "product" in reference["trna"][matched_ref_name].qualifiers:
                    f.qualifiers["product"] = reference["trna"][matched_ref_name].qualifiers["product"]
                    print(f"    Added product: {f.qualifiers['product'][0]}")
                observed["trna"][matched_ref_name] = f
            else:
                observed["issues"].append({
                    "type": "tRNA",
                    "original": raw,
                    "issue": f"Could not match tRNA: '{raw}'"
                })
                print(f"  tRNA: '{raw}' - NO MATCH FOUND")

        elif f.type == "rRNA":
            matched_ref_name = match_rrna(raw, reference["rrna"].keys())
            if matched_ref_name:
                if raw != matched_ref_name:
                    print(f"  rRNA: '{raw}' → '{matched_ref_name}'")
                    changes_made += 1
                # Update gene name to reference
                f.qualifiers["gene"] = [matched_ref_name]
                # Copy reference product if target is missing it
                if "product" not in f.qualifiers and "product" in reference["rrna"][matched_ref_name].qualifiers:
                    f.qualifiers["product"] = reference["rrna"][matched_ref_name].qualifiers["product"]
                    print(f"    Added product: {f.qualifiers['product'][0]}")
                observed["rrna"][matched_ref_name] = f
            else:
                observed["issues"].append({
                    "type": "rRNA",
                    "original": raw,
                    "issue": f"Unmatched rRNA: '{raw}'"
                })
                print(f"  rRNA: '{raw}' - NO MATCH FOUND")
    
    print(f"  Total changes: {changes_made}")
    return revised, observed


# ===============================
# Comparison
# ===============================

def compare(reference, observed, target_name):
    rows = []

    # Compare all CDS genes
    for gene in sorted(reference["genes"].keys()):
        row = {
            "Target": target_name,
            "Gene": gene,
            "Type": "CDS",
            "Present": gene in observed["genes"],
            "Reference_Intron": gene in reference["intron_genes"],
            "Target_Intron": False,
            "Ref_CDS_Length": reference["cds_lengths"].get(gene),
            "Target_CDS_Length": None,
            "Ref_Product": reference["genes"][gene].qualifiers.get("product", [""])[0],
            "Target_Product": None,
            "Issue": "",
            "Verification": ""
        }

        if gene not in observed["genes"]:
            row["Issue"] = "Missing gene"
            row["Verification"] = "May reflect annotation or assembly gap"
        else:
            cds = observed["cds"].get(gene)
            if cds:
                row["Target_CDS_Length"] = cds_length(cds)
                row["Target_Product"] = cds.qualifiers.get("product", [""])[0]
                
                if has_intron(cds) and gene != "rps12":
                    row["Target_Intron"] = True
                    
                if row["Reference_Intron"] and not row["Target_Intron"]:
                    row["Issue"] = "Intron loss"
                    row["Verification"] = "May reflect annotation incompleteness"
                    
                if not row["Target_Product"]:
                    if row["Issue"]:
                        row["Issue"] += "; Missing CDS product"
                    else:
                        row["Issue"] = "Missing CDS product"
                    row["Verification"] = "Annotation incomplete"
            else:
                row["Issue"] = "CDS missing"
                row["Verification"] = "Annotation error likely"

        rows.append(row)

    # Compare all tRNA genes
    for gene in sorted(reference["trna"].keys()):
        row = {
            "Target": target_name,
            "Gene": gene,
            "Type": "tRNA",
            "Present": gene in observed["trna"],
            "Reference_Intron": False,
            "Target_Intron": False,
            "Ref_CDS_Length": None,
            "Target_CDS_Length": None,
            "Ref_Product": reference["trna"][gene].qualifiers.get("product", [""])[0],
            "Target_Product": None,
            "Issue": "",
            "Verification": ""
        }

        if gene not in observed["trna"]:
            row["Issue"] = "Missing tRNA"
            row["Verification"] = "May reflect annotation or assembly gap"
        else:
            trna = observed["trna"][gene]
            row["Target_Product"] = trna.qualifiers.get("product", [""])[0]
            
            if not row["Target_Product"]:
                row["Issue"] = "Missing tRNA product"
                row["Verification"] = "Annotation incomplete"

        rows.append(row)

    # Compare all rRNA genes
    for gene in sorted(reference["rrna"].keys()):
        row = {
            "Target": target_name,
            "Gene": gene,
            "Type": "rRNA",
            "Present": gene in observed["rrna"],
            "Reference_Intron": False,
            "Target_Intron": False,
            "Ref_CDS_Length": None,
            "Target_CDS_Length": None,
            "Ref_Product": reference["rrna"][gene].qualifiers.get("product", [""])[0],
            "Target_Product": None,
            "Issue": "",
            "Verification": ""
        }

        if gene not in observed["rrna"]:
            row["Issue"] = "Missing rRNA"
            row["Verification"] = "May reflect annotation or assembly gap"
        else:
            rrna = observed["rrna"][gene]
            row["Target_Product"] = rrna.qualifiers.get("product", [""])[0]
            
            if not row["Target_Product"]:
                row["Issue"] = "Missing rRNA product"
                row["Verification"] = "Annotation incomplete"

        rows.append(row)

    return rows


# ===============================
# Main
# ===============================

def main():
    parser = argparse.ArgumentParser(
        description="CGAS Module 3: Plastome Gene Comparison and Normalization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare current directory GenBank files to reference
  python cgas_module3.py -r reference.gb
  
  # Compare specific directory to reference
  python cgas_module3.py -r reference.gb -t module2_annotations/Annotated_GenBank/
  
  # Use a close reference for best matching
  python cgas_module3.py -r Hibiscus_reference.gb -t my_genbanks/

Output:
  module_3/
  ├── revise_annotations/          # Normalized GenBank files
  │   ├── sample1.gb
  │   └── sample2.gb
  └── comparison_results.xlsx      # Comparison report (2 sheets)
      ├── Comparison               # Gene-by-gene comparison
      └── Normalization_Issues     # Name matching problems
        """
    )
    parser.add_argument(
        "-r", "--reference", required=True,
        help="Path to reference GenBank file (use closely related species)"
    )
    parser.add_argument(
        "-t", "--targets", default=None,
        help="Directory containing target GenBank files (default: current directory)"
    )

    args = parser.parse_args()
    ensure_dirs()

    reference = parse_reference(args.reference)
    ref_abs = os.path.abspath(args.reference)

    if args.targets:
        target_files = [
            os.path.join(args.targets, f)
            for f in os.listdir(args.targets)
            if is_genbank(f)
        ]
    else:
        target_files = [
            os.path.abspath(f)
            for f in os.listdir(".")
            if is_genbank(f)
        ]

    target_files = [f for f in target_files if os.path.abspath(f) != ref_abs]

    all_rows = []
    all_issues = []

    print(f"\nProcessing {len(target_files)} target file(s)...\n")

    for gb in target_files:
        print(f"Processing: {os.path.basename(gb)}")
        revised, observed = process_target(gb, reference)
        out_gb = os.path.join(REVISED_DIR, os.path.basename(gb))
        SeqIO.write(revised, out_gb, "genbank")
        print(f"  Saved to: {out_gb}\n")

        rows = compare(reference, observed, os.path.basename(gb))
        all_rows.extend(rows)
        
        # Collect normalization issues
        for issue_info in observed["issues"]:
            all_issues.append({
                "Target": os.path.basename(gb),
                "Feature_Type": issue_info["type"],
                "Original_Name": issue_info["original"],
                "Issue": issue_info["issue"]
            })

    df = pd.DataFrame(all_rows)
    
    # Write to Excel with multiple sheets if there are issues
    with pd.ExcelWriter(EXCEL_OUT, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Comparison', index=False)
        
        if all_issues:
            df_issues = pd.DataFrame(all_issues)
            df_issues.to_excel(writer, sheet_name='Normalization_Issues', index=False)

    print(f"\n{'='*80}")
    print(f"CGAS MODULE 3 COMPLETE")
    print(f"{'='*80}")
    print(f"Processed {len(target_files)} target file(s)")
    print(f"\nOutputs:")
    print(f"  • Normalized GenBank files: {REVISED_DIR}/")
    print(f"  • Comparison report: {EXCEL_OUT}")
    
    if all_issues:
        print(f"\n⚠ Found {len(all_issues)} normalization issues")
        print(f"  See 'Normalization_Issues' sheet in Excel report")
    else:
        print(f"\n✓ All gene names matched successfully!")
    
    print(f"{'='*80}\n")


if __name__ == "__main__":
    main()
