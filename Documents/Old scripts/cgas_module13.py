#!/usr/bin/env python3
"""
CGAS Module 13: Comprehensive Nucleotide Diversity Analysis (CORRECTED)
========================================================================

CORRECTIONS IN THIS VERSION (per module13_correction_protocol):
1. Genes are now grouped by gene name ONLY (not by name+position).
   This collapses IRa/IRb duplicate copies so each gene is counted
   once per species, eliminating IR-duplication bias in π values.
2. Reverse-strand sequences (strand = -1) are reverse-complemented
   before alignment so all orientations are comparable.
3. Pseudogenes are detected via the 'pseudo'/'pseudogene' qualifier.
4. Pseudogenes are EXCLUDED from gene diversity calculations but still
   retained in the feature list so intergenic boundaries are preserved.
5. Intergenic spacer names are formed by sorting the two flanking gene
   names alphabetically (e.g., both "trnR-trnN" and "trnN-trnR" become
   the canonical "trnN-trnR"), deduplicating IR spacer pairs.
6. Intron names are keyed by gene name rather than coordinates, so IRa
   and IRb copies of an intron collapse to a single entry.
7. Genomic coordinates are used only as secondary information; biological
   identity (gene name, spacer pair, intron gene) is the primary key.

PREVIOUS IMPROVEMENTS (v1.0):
- Nucleotide diversity values are ROUNDED to 1 decimal place (0.0, 0.1, 0.2, 0.3, etc.)
- Region markers (LSC/SSC/IR) are DISABLED as they were not working with clarity
- All other functionality remains EXACTLY the same

This script performs comprehensive nucleotide diversity (π) analysis for chloroplast
genomes by extracting genes, introns, and intergenic spacers from GenBank files,
aligning them with MAFFT, and calculating nucleotide diversity for each region.

Key Features:
1. Extracts all genes, introns, and intergenic spacers from GenBank files
2. Handles gene name synonyms and tRNA anticodon variations
3. Aligns sequences using MAFFT
4. Calculates nucleotide diversity (π) for each region
5. Generates comprehensive text reports
6. Creates publication-quality R visualizations with:
   - Panel A: Gene diversity plot
   - Panel B: Non-coding regions (introns + intergenic spacers) plot

Author: Abdullah
Version: 1.0.1 (Module 13 - IR-corrected)
Date: March 2026

Dependencies:
    Python: biopython, numpy
    External: MAFFT (alignment tool)
    R: ggplot2, readr, cowplot (for visualization)

Usage:
    python cgas_module13.py
    python cgas_module13.py -i genbank_files/
    python cgas_module13.py file1.gb file2.gb
    python cgas_module13.py -i data/ -o results/
    python cgas_module13.py --no-figures  # Skip R visualization

Output:
    Module13_Diversity_Analysis/
        - genes/ (extracted gene sequences)
        - introns/ (extracted intron sequences)
        - intergenic/ (extracted intergenic spacer sequences)
        - alignments/ (MAFFT alignments)
        - nucleotide_diversity_results.txt (comprehensive report)
        - nucleotide_diversity_by_position.txt (ordered by genomic position)
        - nucleotide_diversity_all_regions_by_position.txt (all regions combined)
        - diversity_data_for_r.txt (data formatted for R visualization)
        - Figures/ (if R available):
            - nucleotide_diversity_plot.pdf and .png

Notes:
    - Requires MAFFT to be installed and accessible in PATH
    - Handles gene name synonyms automatically
    - tRNA genes with anticodon mismatches are flagged in the report
    - All gene names in figures are properly italicized
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import numpy as np
import re

def to_roman(num):
    """Convert integer to Roman numeral (for intron numbering)."""
    val = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1]
    syms = ['M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I']
    roman_num = ''
    for i in range(len(val)):
        count = int(num / val[i])
        if count:
            roman_num += syms[i] * count
            num -= val[i] * count
    return roman_num

def normalize_gene_name(gene_name):
    """
    Normalize gene names to handle different naming conventions:
    - tRNA: trnN-GUU, trnN_GUU, trnN_guu, trnN(GUU) -> trnK-UUU (uppercase)
    - rRNA: normalize case
    - Gene synonyms: ycf3/pafI/paf1, ycf4/pafII, psbN/psb1, psbZ/lhbA, clpP/clpP1/clp1, etc.
    """
    if not gene_name:
        return gene_name
    
    # Gene synonym mapping - maps alternative names to the primary/standard name
    gene_synonyms = {
        # ycf genes
        'pafi': 'ycf3', 'paf1': 'ycf3', 'pafI': 'ycf3',
        'pafii': 'ycf4', 'paf2': 'ycf4', 'pafII': 'ycf4',
        
        # psb genes
        'psb1': 'psbN', 'psbTn': 'psbN',
        'lhba': 'psbZ', 'lhbA': 'psbZ',
        
        # clp genes
        'clp1': 'clpP', 'clpp1': 'clpP', 'clpP1': 'clpP',
        
        # rpo genes
        'rpoa': 'rpoA', 'rpob': 'rpoB', 'rpoc1': 'rpoC1', 'rpoc2': 'rpoC2',
        
        # rpl genes (ribosomal proteins large subunit)
        'rpl2': 'rpl2', 'rpl14': 'rpl14', 'rpl16': 'rpl16', 'rpl20': 'rpl20',
        'rpl22': 'rpl22', 'rpl23': 'rpl23', 'rpl32': 'rpl32', 'rpl33': 'rpl33',
        'rpl36': 'rpl36',
        
        # rps genes (ribosomal proteins small subunit)
        'rps2': 'rps2', 'rps3': 'rps3', 'rps4': 'rps4', 'rps7': 'rps7',
        'rps8': 'rps8', 'rps11': 'rps11', 'rps12': 'rps12', 'rps14': 'rps14',
        'rps15': 'rps15', 'rps16': 'rps16', 'rps18': 'rps18', 'rps19': 'rps19',
        
        # inf gene
        'infa': 'infA',
        
        # mat gene
        'matk': 'matK',
        
        # cem gene
        'cema': 'cemA',
        
        # acc gene
        'accd': 'accD',
        
        # ndh genes
        'ndha': 'ndhA', 'ndhb': 'ndhB', 'ndhc': 'ndhC', 'ndhd': 'ndhD',
        'ndhe': 'ndhE', 'ndhf': 'ndhF', 'ndhg': 'ndhG', 'ndhh': 'ndhH',
        'ndhi': 'ndhI', 'ndhj': 'ndhJ', 'ndhk': 'ndhK',
        
        # pet genes
        'peta': 'petA', 'petb': 'petB', 'petd': 'petD', 'petg': 'petG',
        'petl': 'petL', 'petn': 'petN',
        
        # psb genes (photosystem II)
        'psba': 'psbA', 'psbb': 'psbB', 'psbc': 'psbC', 'psbd': 'psbD',
        'psbe': 'psbE', 'psbf': 'psbF', 'psbh': 'psbH', 'psbi': 'psbI',
        'psbj': 'psbJ', 'psbk': 'psbK', 'psbl': 'psbL', 'psbm': 'psbM',
        'psbt': 'psbT', 'psbz': 'psbZ',
        
        # psa genes (photosystem I)
        'psaa': 'psaA', 'psab': 'psaB', 'psac': 'psaC', 'psai': 'psaI',
        'psaj': 'psaJ', 'psam': 'psaM',
        
        # atp genes
        'atpa': 'atpA', 'atpb': 'atpB', 'atpe': 'atpE', 'atpf': 'atpF',
        'atph': 'atpH', 'atpi': 'atpI',
        
        # rbc genes
        'rbcl': 'rbcL',
        
        # ccs gene
        'ccsa': 'ccsA',
        
        # ycf genes
        'ycf1': 'ycf1', 'ycf2': 'ycf2', 'ycf3': 'ycf3', 'ycf4': 'ycf4',
    }
    
    # Handle tRNA naming conventions for ALL tRNA genes
    if gene_name.lower().startswith('trn'):
        # Extract the amino acid and anticodon
        # Patterns: trnN-GUU, trnN_GUU, trnN(GUU), trnN_guu, trnN-guu, etc.
        match = re.match(r'(trn[A-Z])[_\-\(]?([A-Z]{3})\)?', gene_name, re.IGNORECASE)
        if match:
            amino_acid = match.group(1)  # Keep original case from match
            # Ensure proper capitalization: trnK not trnk
            amino_acid = 'trn' + amino_acid[3].upper()
            anticodon = match.group(2).upper()   # e.g., 'GUU'
            return f"{amino_acid}-{anticodon}"
        else:
            # Just trnN without anticodon (e.g., trnH, trnN)
            match = re.match(r'(trn[A-Z])', gene_name, re.IGNORECASE)
            if match:
                amino_acid = match.group(1)
                return 'trn' + amino_acid[3].upper()  # trnH not trnh
    
    # Handle rRNA naming conventions
    if gene_name.lower().startswith('rrn'):
        # Normalize to lowercase
        return gene_name.lower()
    
    # Check if gene has a synonym and use the standard name
    gene_lower = gene_name.lower()
    if gene_lower in gene_synonyms:
        return gene_synonyms[gene_lower]
    
    # For other genes, return as is (keep original case)
    return gene_name

def get_organism_label(record):
    """
    Return a clean species label from a GenBank record.
    Uses the organism annotation (e.g. 'Gossypium arboreum') with spaces
    replaced by underscores.  Falls back to the accession if not present.
    This label is used as the FASTA sequence ID so alignments are labelled
    by species, not by accession number.
    """
    organism = record.annotations.get("organism", "")
    if not organism:
        for feature in record.features:
            if feature.type == "source" and "organism" in feature.qualifiers:
                organism = feature.qualifiers["organism"][0]
                break
    if organism:
        # Replace spaces/special chars; keep it MAFFT-safe
        label = re.sub(r"[^A-Za-z0-9_]", "_", organism).strip("_")
        return label
    return record.id


def extract_features_from_genbank(genbank_files, output_dir):
    """
    Extract genes, introns, and intergenic spacers from multiple GenBank files.

    Core guarantee
    --------------
    Every alignment file will contain EXACTLY ONE sequence per species for
    each region.  If a region is absent or is a pseudogene in a species it
    simply will not appear in that species' contribution — the region will
    then have fewer than N sequences and will be noted in the report.

    How IR duplication is prevented
    --------------------------------
    - Genes   : grouped by normalised gene name; the first (lowest-position)
                copy is extracted; a (species_label, gene_name) guard blocks
                any second copy.
    - Introns : extracted from the same first copy; same (species_label,
                intron_name) guard.
    - Spacers : spacer name = sorted([gene1, gene2]) so IRa and IRb produce
                the same key; a (species_label, spacer_name) guard keeps only
                the first occurrence.

    Orientation normalisation
    -------------------------
    All sequences are stored in the FORWARD (5'→3') orientation relative to
    the gene / spacer as it appears in the LSC / forward-strand context:
    - Genes   : Bio.SeqFeature.extract() already reverse-complements minus-
                strand features automatically.
    - Introns : the raw genomic slice is reverse-complemented when the host
                gene is on the minus strand.
    - Spacers : the raw genomic slice is reverse-complemented when BOTH
                flanking genes are on the minus strand (i.e. the spacer is
                in the IRb context).

    Returns
    -------
    genes_dict, intron_dict, intergenic_dict,
    ambiguous_genes, gene_without_anticodon, synonym_mapping
    """
    genes_dict       = defaultdict(list)
    intron_dict      = defaultdict(list)
    intergenic_dict  = defaultdict(list)
    ambiguous_genes        = set()
    gene_without_anticodon = set()
    synonym_mapping        = defaultdict(set)

    # Per-species deduplication guards:
    # {(species_label, region_name)} — one entry means "already stored".
    seen_genes   = set()
    seen_introns = set()
    seen_spacers = set()

    for gb_file in genbank_files:
        print(f"Processing {gb_file}...")

        for record in SeqIO.parse(gb_file, "genbank"):
            sequence      = record.seq
            # Use organism name as the sequence label (not accession)
            species_label = get_organism_label(record)
            print(f"  Species: {species_label}")

            # ----------------------------------------------------------------
            # Collect all coding/RNA features, sorted by start position.
            # Keep "gene"-type features out — we only want CDS/tRNA/rRNA so
            # we don't double-count.  Pseudogenes stay so they anchor spacers.
            # ----------------------------------------------------------------
            all_features = [
                f for f in record.features
                if f.type in ["CDS", "tRNA", "rRNA"]
            ]
            all_features.sort(key=lambda f: int(f.location.start))

            # ----------------------------------------------------------------
            # Build features_by_gene: normalised_name → [feature, ...]
            # Grouping by name (not position) already merges IRa+IRb under
            # one key per species.
            # ----------------------------------------------------------------
            features_by_gene = defaultdict(list)

            for feature in all_features:
                raw_name = (
                    feature.qualifiers.get("gene",    [None])[0] or
                    feature.qualifiers.get("product", [None])[0] or
                    feature.qualifiers.get("locus_tag",[None])[0]
                )
                if raw_name is None:
                    raw_name = f"unknown_{feature.location.start}"
                    ambiguous_genes.add(raw_name)

                # Flag tRNA without anticodon info
                if raw_name.lower().startswith("trn"):
                    has_ac = ('-' in raw_name or '_' in raw_name or
                              '(' in raw_name or
                              'anticodon' in feature.qualifiers)
                    if not has_ac:
                        gene_without_anticodon.add(
                            f"{raw_name} (in {species_label} at "
                            f"{feature.location.start}-{feature.location.end})"
                        )

                norm_name = normalize_gene_name(raw_name)
                if raw_name != norm_name:
                    synonym_mapping[norm_name].add(raw_name)

                features_by_gene[norm_name].append(feature)

            # ----------------------------------------------------------------
            # Process each gene — use the FIRST (lowest-position) copy only.
            # ----------------------------------------------------------------
            for gene_name, gfeatures in features_by_gene.items():
                gfeatures.sort(key=lambda f: int(f.location.start))
                first = gfeatures[0]   # canonical copy (LSC / lowest position)
                strand = first.location.strand

                # Pseudogene check
                is_pseudo = (
                    "pseudo"     in first.qualifiers or
                    "pseudogene" in first.qualifiers
                )
                if is_pseudo:
                    print(f"    [pseudo] {gene_name} — excluded from diversity")

                # ----------------------------------------------------------
                # GENE sequence
                # extract() handles reverse-complement for minus-strand genes
                # automatically, so no manual RC needed here.
                # ----------------------------------------------------------
                gene_key = (species_label, gene_name)
                if not is_pseudo and gene_key not in seen_genes:
                    seen_genes.add(gene_key)
                    gene_seq = first.location.extract(sequence)
                    genes_dict[gene_name].append(
                        SeqRecord(gene_seq,
                                  id=species_label,
                                  description="")
                    )

                # ----------------------------------------------------------
                # INTRON sequences — derived from the same first copy
                # ----------------------------------------------------------
                exon_parts = (
                    list(first.location.parts)
                    if hasattr(first.location, "parts")
                    else [first.location]
                )
                # Sort exons by genomic coordinate (always ascending)
                exon_parts.sort(key=lambda x: int(x.start))

                if len(exon_parts) < 2:
                    continue   # no introns in this gene

                gene_span = int(exon_parts[-1].end) - int(exon_parts[0].start)

                # Skip trans-spliced genes (except rps12 handled below)
                if gene_name.lower() != "rps12" and gene_span > 10000:
                    print(f"    Skipping trans-spliced {gene_name} "
                          f"(span {gene_span} bp)")
                    continue

                # Collect valid introns
                valid_introns = []
                for i in range(len(exon_parts) - 1):
                    i_start = int(exon_parts[i].end)
                    i_end   = int(exon_parts[i + 1].start)
                    if i_end <= i_start:
                        continue
                    i_len = i_end - i_start
                    if i_len > 10000:
                        print(f"    Skipping trans-spliced intron in "
                              f"{gene_name} ({i_len} bp)")
                        continue
                    # Raw genomic slice, then orient to match gene strand
                    i_seq = sequence[i_start:i_end]
                    if strand == -1:
                        i_seq = i_seq.reverse_complement()
                    valid_introns.append((i_start, i_end, i_seq, i_len))

                num_introns = len(valid_introns)
                for idx, (i_start, i_end, i_seq, i_len) in \
                        enumerate(valid_introns, 1):
                    intron_name = (
                        f"{gene_name} intron"
                        if num_introns == 1
                        else f"{gene_name} intron {to_roman(idx)}"
                    )
                    intron_key = (species_label, intron_name)
                    if intron_key not in seen_introns:
                        seen_introns.add(intron_key)
                        intron_dict[intron_name].append(
                            SeqRecord(i_seq,
                                      id=species_label,
                                      description="")
                        )
                        print(f"    Extracted {intron_name}: "
                              f"{i_start+1}..{i_end} ({i_len} bp)")
                    else:
                        print(f"    Skipping duplicate {intron_name} "
                              f"for {species_label} (IR copy)")

            # ----------------------------------------------------------------
            # INTERGENIC SPACERS
            # Iterate consecutive feature pairs; normalise spacer name by
            # sorting the two flanking gene names so IRa and IRb yield the
            # same key.  Reverse-complement if both flanking genes are on
            # the minus strand (= IRb context).
            # ----------------------------------------------------------------
            for i in range(len(all_features) - 1):
                cur  = all_features[i]
                nxt  = all_features[i + 1]

                g1_raw = (cur.qualifiers.get("gene",    [None])[0] or
                          cur.qualifiers.get("product", [None])[0])
                g2_raw = (nxt.qualifiers.get("gene",    [None])[0] or
                          nxt.qualifiers.get("product", [None])[0])
                if not g1_raw or not g2_raw:
                    continue

                g1 = normalize_gene_name(g1_raw)
                g2 = normalize_gene_name(g2_raw)

                # Canonical spacer name: alphabetical order of the two genes
                spacer_name = "-".join(sorted([g1, g2]))

                sp_start = int(cur.location.end)
                sp_end   = int(nxt.location.start)

                if sp_end <= sp_start:
                    continue          # overlapping features — no spacer
                if (sp_end - sp_start) <= 10:
                    continue          # too short to be meaningful

                spacer_key = (species_label, spacer_name)
                if spacer_key in seen_spacers:
                    print(f"    Skipping duplicate spacer {spacer_name} "
                          f"for {species_label} (IR copy)")
                    continue

                sp_seq = sequence[sp_start:sp_end]

                # Orientation: reverse-complement when BOTH flanking genes
                # are on the minus strand (IRb context).
                if cur.location.strand == -1 and nxt.location.strand == -1:
                    sp_seq = sp_seq.reverse_complement()

                seen_spacers.add(spacer_key)
                intergenic_dict[spacer_name].append(
                    SeqRecord(sp_seq,
                              id=species_label,
                              description="")
                )

    return (genes_dict, intron_dict, intergenic_dict,
            ambiguous_genes, gene_without_anticodon, synonym_mapping)

def get_feature_order(genbank_file):
    """
    Build a unified ordered list of ALL regions in the order they appear in
    the reference genome:

        gene1  ->  gene1_intron(s)  ->  gene1-gene2_spacer  ->  gene2  ->  ...

    This ensures that IGS and intron results follow exactly the same
    left-to-right genomic pattern as the gene track.

    Returns
    -------
    gene_positions        : {gene_name:   genomic_start}
    intron_positions      : {intron_name: genomic_start}
    intergenic_positions  : {sorted_key:  genomic_start}
    canonical_spacer_names: {sorted_key -> directional_name}
                            Maps sorted extraction key (e.g. "trnN-trnR")
                            to the directional reference name ("trnR-trnN").
    unified_order         : [(region_type, display_name, genomic_pos), ...]
                            Full interleaved order used to sort Panel B.
    """
    gene_positions         = {}
    intron_positions       = {}
    intergenic_positions   = {}
    canonical_spacer_names = {}
    unified_order          = []

    for record in SeqIO.parse(genbank_file, "genbank"):
        all_features = [f for f in record.features
                        if f.type in ["CDS", "tRNA", "rRNA"]]
        all_features.sort(key=lambda x: int(x.location.start))

        # Deduplicate by gene name — keep only first (lowest-position) copy
        seen = {}
        deduped = []
        for feature in all_features:
            raw = (feature.qualifiers.get("gene",    [None])[0] or
                   feature.qualifiers.get("product", [None])[0])
            if not raw:
                continue
            norm = normalize_gene_name(raw)
            if norm not in seen:
                seen[norm] = True
                gene_positions[norm] = int(feature.location.start)
                deduped.append((norm, feature))

        # Build unified order: gene -> its introns -> spacer to next gene -> ...
        for idx, (gene_name, feature) in enumerate(deduped):
            gene_pos = int(feature.location.start)
            unified_order.append(("Gene", gene_name, gene_pos))

            # Introns of this gene
            exon_parts = (list(feature.location.parts)
                          if hasattr(feature.location, "parts")
                          else [feature.location])
            exon_parts.sort(key=lambda x: int(x.start))

            if len(exon_parts) >= 2:
                valid_introns = []
                for i in range(len(exon_parts) - 1):
                    i_start = int(exon_parts[i].end)
                    i_end   = int(exon_parts[i + 1].start)
                    i_len   = i_end - i_start
                    if 0 < i_len <= 10000:
                        valid_introns.append(i_start)

                num = len(valid_introns)
                for i2, i_start in enumerate(valid_introns, 1):
                    iname = (f"{gene_name} intron"
                             if num == 1
                             else f"{gene_name} intron {to_roman(i2)}")
                    if iname not in intron_positions:
                        intron_positions[iname] = i_start
                    unified_order.append(("Intron", iname, i_start))

            # Spacer between this gene and the next
            if idx + 1 < len(deduped):
                next_name, next_feat = deduped[idx + 1]
                sp_start = int(feature.location.end)
                sp_end   = int(next_feat.location.start)

                if sp_end > sp_start and (sp_end - sp_start) > 10:
                    directional = f"{gene_name}-{next_name}"
                    sorted_key  = "-".join(sorted([gene_name, next_name]))

                    if sorted_key not in intergenic_positions:
                        intergenic_positions[sorted_key] = sp_start
                    if sorted_key not in canonical_spacer_names:
                        canonical_spacer_names[sorted_key] = directional

                    unified_order.append(("Intergenic", sorted_key, sp_start))

        break  # Only the first record defines the reference order

    return (gene_positions, intron_positions, intergenic_positions,
            canonical_spacer_names, unified_order)

def write_fasta_files(features_dict, output_dir, prefix):
    """Write sequences to FASTA files for each feature."""
    fasta_files = []
    
    for feature_name, seq_records in features_dict.items():
        if len(seq_records) < 2:
            print(f"Skipping {feature_name}: only {len(seq_records)} sequence(s)")
            continue
        
        # Sanitize filename
        safe_name = feature_name.replace("/", "_").replace(" ", "_")
        fasta_file = os.path.join(output_dir, f"{prefix}_{safe_name}.fasta")
        
        SeqIO.write(seq_records, fasta_file, "fasta")
        fasta_files.append((feature_name, fasta_file))
    
    return fasta_files

def align_with_mafft(fasta_file, output_file):
    """Align sequences using MAFFT."""
    try:
        cmd = ["mafft", "--auto", "--quiet", fasta_file]
        with open(output_file, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)
        return True
    except subprocess.CalledProcessError:
        print(f"Error aligning {fasta_file}")
        return False
    except FileNotFoundError:
        print("MAFFT not found. Please install MAFFT.")
        return False

def calculate_nucleotide_diversity(alignment_file):
    """
    Calculate nucleotide diversity (π) from aligned sequences.
    π = average number of nucleotide differences per site between any two sequences.
    """
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
    except:
        return None, None, None
    
    n_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()
    
    if n_sequences < 2:
        return 0.0, 0, alignment_length
    
    # Calculate pairwise differences
    total_differences = 0
    total_comparisons = 0
    valid_sites = 0
    
    # Check each site
    for i in range(alignment_length):
        column = alignment[:, i]
        
        # Skip sites with gaps
        if '-' in column or 'N' in column.upper():
            continue
        
        valid_sites += 1
        
        # Count pairwise differences at this site
        for j in range(n_sequences):
            for k in range(j + 1, n_sequences):
                if column[j].upper() != column[k].upper():
                    total_differences += 1
                total_comparisons += 1
    
    # Calculate diversity
    if total_comparisons == 0 or valid_sites == 0:
        return 0.0, 0, alignment_length
    
    # π = average pairwise differences per site
    pi = total_differences / total_comparisons
    
    return pi, valid_sites, alignment_length

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='CGAS Module 13: Nucleotide Diversity Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s                          # Use current directory
  %(prog)s genbank_files/           # Process specific directory
  %(prog)s file1.gb file2.gb        # Process specific files
  %(prog)s -o results/              # Custom output directory
  %(prog)s --figures-only           # Regenerate figures from existing data
  %(prog)s --figures-only -o Module13_Diversity_Analysis/

For more information, see the CGAS documentation.
        '''
    )
    
    parser.add_argument('input_paths',
                       nargs='*',
                       default=['.'],
                       help='GenBank files or directories (default: current directory)')
    
    parser.add_argument('-o', '--output',
                       type=str,
                       default='Module13_Diversity_Analysis',
                       help='Output directory (default: Module13_Diversity_Analysis)')
    
    parser.add_argument('--figures-only',
                       action='store_true',
                       help='Skip analysis and regenerate figures from existing data files')
    
    args = parser.parse_args()
    
    output_folder = args.output
    
    # ========================================================================
    # FIGURES-ONLY MODE: Regenerate plots from existing data
    # ========================================================================
    if args.figures_only:
        print("\n" + "="*80)
        print("CGAS MODULE 13: REGENERATING FIGURES")
        print("="*80)
        print(f"\nOutput directory: {output_folder}")
        
        # Check if output directory exists
        if not os.path.exists(output_folder):
            print(f"\n✗ Error: Output directory '{output_folder}' not found")
            print("  Run the full analysis first, or specify correct output directory with -o")
            sys.exit(1)
        
        # Check for required data files
        genes_file = os.path.join(output_folder, "genes_for_R_plot.txt")
        noncoding_file = os.path.join(output_folder, "noncoding_for_R_plot.txt")
        
        if not os.path.exists(genes_file):
            print(f"\n✗ Error: Required file not found: {genes_file}")
            print("  Run the full analysis first")
            sys.exit(1)
        
        if not os.path.exists(noncoding_file):
            print(f"\n✗ Error: Required file not found: {noncoding_file}")
            print("  Run the full analysis first")
            sys.exit(1)
        
        print(f"\n✓ Found data files:")
        print(f"  - {os.path.basename(genes_file)}")
        print(f"  - {os.path.basename(noncoding_file)}")
        
        # Generate R visualizations
        generate_r_visualizations(genes_file, noncoding_file, output_folder)
        
        print("\n" + "="*80)
        print("FIGURE REGENERATION COMPLETE")
        print("="*80)
        print(f"\nFigures saved in: {output_folder}/Figures/")
        print("="*80 + "\n")
        
        return
    
    # ========================================================================
    # NORMAL MODE: Full analysis
    # ========================================================================
    
    # Configuration
    input_paths = args.input_paths if args.input_paths else ['.']
    
    genbank_files = []
    
    # Process input paths - handle both files and directories
    import glob
    for path in input_paths:
        if not os.path.exists(path):
            print(f"Error: Path {path} not found")
            sys.exit(1)
        
        if os.path.isdir(path):
            # If it's a directory, find all GenBank files
            print(f"Searching for GenBank files in directory: {path}")
            for ext in ['*.gb', '*.gbk', '*.genbank', '*.gbff']:
                found_files = glob.glob(os.path.join(path, ext))
                genbank_files.extend(found_files)
        else:
            # It's a file
            genbank_files.append(path)
    
    if not genbank_files:
        print("Error: No GenBank files found")
        print("Supported extensions: .gb, .gbk, .genbank, .gbff")
        print("\nUsage:")
        print("  Run without arguments to use current directory:")
        print("    python module9_diversity_analysis.py")
        print("  Or provide explicit paths:")
        print("    python module9_diversity_analysis.py <genbank_file1> <genbank_file2> ...")
        print("    python module9_diversity_analysis.py <directory>")
        print("    python module9_diversity_analysis.py *.gb")
        sys.exit(1)
    
    print(f"Processing {len(genbank_files)} GenBank files...")
    
    # Create output directories
    os.makedirs(output_folder, exist_ok=True)
    
    genes_dir = os.path.join(output_folder, "genes")
    intron_dir = os.path.join(output_folder, "introns")
    intergenic_dir = os.path.join(output_folder, "intergenic")
    alignments_dir = os.path.join(output_folder, "alignments")
    
    os.makedirs(genes_dir, exist_ok=True)
    os.makedirs(intron_dir, exist_ok=True)
    os.makedirs(intergenic_dir, exist_ok=True)
    os.makedirs(alignments_dir, exist_ok=True)
    
    # Extract features
    print("\nExtracting genes, introns, and intergenic regions...")
    genes_dict, intron_dict, intergenic_dict, ambiguous_genes, genes_without_anticodon, synonym_mapping = extract_features_from_genbank(genbank_files, output_folder)
    
    print(f"\nFound {len(genes_dict)} unique genes")
    print(f"Found {len(intron_dict)} unique introns")
    print(f"Found {len(intergenic_dict)} intergenic regions")
    
    # Report gene synonyms that were merged
    if synonym_mapping:
        print(f"\n✓ Merged {len(synonym_mapping)} genes with alternative names:")
        for standard_name, original_names in sorted(synonym_mapping.items()):
            if len(original_names) > 1:
                names_str = ", ".join(sorted(original_names))
                print(f"  {standard_name} ← [{names_str}]")
    
    # Report problematic genes
    if ambiguous_genes:
        print(f"\n⚠ Warning: {len(ambiguous_genes)} genes without clear identification:")
        for gene in sorted(ambiguous_genes):
            print(f"  - {gene}")
    
    if genes_without_anticodon:
        print(f"\n⚠ Warning: {len(genes_without_anticodon)} tRNA genes without anticodon information:")
        for gene in sorted(genes_without_anticodon):
            print(f"  - {gene}")
        print("  Note: These tRNAs may not align properly with tRNAs from other species")
        print("        that have anticodon information. Consider manual inspection.")
    
    # Write FASTA files
    print("\nWriting FASTA files...")
    gene_fasta_files = write_fasta_files(genes_dict, genes_dir, "gene")
    intron_fasta_files = write_fasta_files(intron_dict, intron_dir, "intron")
    intergenic_fasta_files = write_fasta_files(intergenic_dict, intergenic_dir, "spacer")
    
    # Align and calculate diversity
    results = []
    
    print("\nAligning genes and calculating diversity...")
    for gene_name, fasta_file in gene_fasta_files:
        alignment_file = os.path.join(alignments_dir, f"gene_{os.path.basename(fasta_file)}")
        n_seqs = len(genes_dict.get(gene_name, []))
        if align_with_mafft(fasta_file, alignment_file):
            pi, valid_sites, total_length = calculate_nucleotide_diversity(alignment_file)
            if pi is not None:
                results.append(("Gene", gene_name, pi, valid_sites, total_length, n_seqs))
                print(f"  {gene_name}: π = {pi:.6f}  (n={n_seqs})")
    
    print("\nAligning introns and calculating diversity...")
    for intron_name, fasta_file in intron_fasta_files:
        alignment_file = os.path.join(alignments_dir, f"intron_{os.path.basename(fasta_file)}")
        n_seqs = len(intron_dict.get(intron_name, []))
        if align_with_mafft(fasta_file, alignment_file):
            pi, valid_sites, total_length = calculate_nucleotide_diversity(alignment_file)
            if pi is not None:
                results.append(("Intron", intron_name, pi, valid_sites, total_length, n_seqs))
                print(f"  {intron_name}: π = {pi:.6f}  (n={n_seqs})")
    
    print("\nAligning intergenic regions and calculating diversity...")
    for spacer_name, fasta_file in intergenic_fasta_files:
        alignment_file = os.path.join(alignments_dir, f"spacer_{os.path.basename(fasta_file)}")
        n_seqs = len(intergenic_dict.get(spacer_name, []))
        
        if align_with_mafft(fasta_file, alignment_file):
            pi, valid_sites, total_length = calculate_nucleotide_diversity(alignment_file)
            if pi is not None:
                results.append(("Intergenic", spacer_name, pi, valid_sites, total_length, n_seqs))
                print(f"  {spacer_name}: π = {pi:.6f}  (n={n_seqs})")
    # ======================================================================
    # Build ordered results using the reference genome unified order
    # ======================================================================
    print(f"\nGenerating genomic position-ordered results...")
    gene_positions, intron_positions, intergenic_positions, \
        canonical_spacer_names, unified_order = get_feature_order(genbank_files[0])

    n_total_species = len(genbank_files)

    # Build lookup: (region_type, sorted_key) -> (pi, valid_sites, length, n_seqs)
    results_lookup = {}
    for tup in results:
        rtype, rname, pi, vs, ln, ns = tup
        results_lookup[(rtype, rname)] = (pi, vs, ln, ns)

    # Build ordered list following unified_order (gene -> intron(s) -> igs -> gene -> ...)
    ordered_all = []
    for rtype, rname, rpos in unified_order:
        if rtype == "Intergenic":
            display = canonical_spacer_names.get(rname, rname)
        else:
            display = rname
        key = (rtype, rname)
        if key in results_lookup:
            pi, vs, ln, ns = results_lookup[key]
            ordered_all.append((rpos, rtype, display, rname, pi, vs, ln, ns))

    # ======================================================================
    # Main results text file  (genomic order, with n_sequences column)
    # ======================================================================
    results_file = os.path.join(output_folder, "nucleotide_diversity_results.txt")
    with open(results_file, "w") as f:
        W = 98
        f.write("Nucleotide Diversity Analysis Results\n")
        f.write("=" * W + "\n\n")
        f.write(f"Input files   : {n_total_species} GenBank files\n")
        f.write(f"Total regions : {len(results)}\n\n")

        if synonym_mapping:
            f.write("GENE SYNONYMS MERGED:\n")
            f.write("-" * W + "\n")
            for std, origs in sorted(synonym_mapping.items()):
                if len(origs) > 1:
                    f.write(f"  {std} <- [{', '.join(sorted(origs))}]\n")
            f.write("\n")

        if ambiguous_genes or genes_without_anticodon:
            f.write("WARNINGS:\n")
            f.write("-" * W + "\n")
            if ambiguous_genes:
                f.write(f"\nGenes without clear identification ({len(ambiguous_genes)}):\n")
                for g in sorted(ambiguous_genes):
                    f.write(f"  - {g}\n")
            if genes_without_anticodon:
                f.write(f"\ntRNA genes without anticodon ({len(genes_without_anticodon)}):\n")
                for g in sorted(genes_without_anticodon):
                    f.write(f"  - {g}\n")
                f.write("\nNote: These tRNAs may not align properly across species.\n")
            f.write("\n")

        f.write("-" * W + "\n")
        f.write(f"{'Type':<12} {'Region':<38} {'n_seq':<7} {'pi':<12} {'Valid_Sites':<13} {'Aln_Len':<10}\n")
        f.write("-" * W + "\n")
        for rpos, rtype, display, rname, pi, vs, ln, ns in ordered_all:
            flag = f"  <- missing {n_total_species - ns}" if ns < n_total_species else ""
            f.write(f"{rtype:<12} {display:<38} {ns:<7} {pi:<12.6f} {vs:<13} {ln:<10}{flag}\n")
        f.write("-" * W + "\n\n")

        # Summary statistics
        gene_pi = [x[4] for x in ordered_all if x[1] == "Gene"]
        nc_pi   = [x[4] for x in ordered_all if x[1] != "Gene"]
        intr_pi = [x[4] for x in ordered_all if x[1] == "Intron"]
        igs_pi  = [x[4] for x in ordered_all if x[1] == "Intergenic"]

        if gene_pi:
            f.write("Gene Statistics:\n")
            f.write(f"  Regions  : {len(gene_pi)}\n")
            f.write(f"  Mean pi  : {np.mean(gene_pi):.6f}\n")
            f.write(f"  Median pi: {np.median(gene_pi):.6f}\n")
            f.write(f"  Min pi   : {np.min(gene_pi):.6f}\n")
            f.write(f"  Max pi   : {np.max(gene_pi):.6f}\n\n")
        if intr_pi:
            f.write("Intron Statistics:\n")
            f.write(f"  Regions  : {len(intr_pi)}\n")
            f.write(f"  Mean pi  : {np.mean(intr_pi):.6f}\n")
            f.write(f"  Median pi: {np.median(intr_pi):.6f}\n")
            f.write(f"  Min pi   : {np.min(intr_pi):.6f}\n")
            f.write(f"  Max pi   : {np.max(intr_pi):.6f}\n\n")
        if igs_pi:
            f.write("Intergenic Spacer Statistics:\n")
            f.write(f"  Regions  : {len(igs_pi)}\n")
            f.write(f"  Mean pi  : {np.mean(igs_pi):.6f}\n")
            f.write(f"  Median pi: {np.median(igs_pi):.6f}\n")
            f.write(f"  Min pi   : {np.min(igs_pi):.6f}\n")
            f.write(f"  Max pi   : {np.max(igs_pi):.6f}\n")

    print(f"\n{'='*80}")
    print(f"Analysis complete!")
    print(f"Results written to: {results_file}")
    print(f"FASTA files in: {genes_dir}, {intron_dir}, and {intergenic_dir}")
    print(f"Alignments in: {alignments_dir}")
    if ambiguous_genes or genes_without_anticodon:
        print(f"\n⚠ Check the results file for warnings")

    # ======================================================================
    # Ordered text file  (coding / non-coding separate sections)
    # ======================================================================
    ordered_results_file = os.path.join(output_folder, "nucleotide_diversity_by_position.txt")
    W = 98
    with open(ordered_results_file, "w") as f:
        f.write("Nucleotide Diversity Results Organised by Genomic Position\n")
        f.write("=" * W + "\n\n")
        f.write("PROTEIN-CODING AND RNA GENES\n")
        f.write("-" * W + "\n")
        f.write(f"{'Pos':<12} {'Gene':<38} {'n_seq':<7} {'pi':<12} {'Valid_Sites':<13} {'Aln_Len':<10}\n")
        f.write("-" * W + "\n")
        for rpos, rtype, display, rname, pi, vs, ln, ns in ordered_all:
            if rtype == "Gene":
                flag = f"  <- missing {n_total_species - ns}" if ns < n_total_species else ""
                f.write(f"{rpos:<12} {display:<38} {ns:<7} {pi:<12.6f} {vs:<13} {ln:<10}{flag}\n")
        f.write("\n\nNON-CODING REGIONS (introns and IGS in genome order)\n")
        f.write("-" * W + "\n")
        f.write(f"{'Pos':<12} {'Type':<12} {'Region':<38} {'n_seq':<7} {'pi':<12} {'Valid_Sites':<13} {'Aln_Len':<10}\n")
        f.write("-" * W + "\n")
        for rpos, rtype, display, rname, pi, vs, ln, ns in ordered_all:
            if rtype != "Gene":
                flag = f"  <- missing {n_total_species - ns}" if ns < n_total_species else ""
                f.write(f"{rpos:<12} {rtype:<12} {display:<38} {ns:<7} {pi:<12.6f} {vs:<13} {ln:<10}{flag}\n")

    # ======================================================================
    # Combined all-regions file  (full interleaved genomic order)
    # ======================================================================
    combined_ordered_file = os.path.join(output_folder, "nucleotide_diversity_all_regions_by_position.txt")
    with open(combined_ordered_file, "w") as f:
        f.write("All Regions in Genomic Order  (gene → intron → IGS → gene → ...)\n")
        f.write("=" * W + "\n\n")
        f.write(f"{'Pos':<12} {'Type':<12} {'Region':<38} {'n_seq':<7} {'pi':<12} {'Valid_Sites':<13} {'Aln_Len':<10}\n")
        f.write("-" * W + "\n")
        for rpos, rtype, display, rname, pi, vs, ln, ns in ordered_all:
            flag = f"  <- missing {n_total_species - ns}" if ns < n_total_species else ""
            f.write(f"{rpos:<12} {rtype:<12} {display:<38} {ns:<7} {pi:<12.6f} {vs:<13} {ln:<10}{flag}\n")
        f.write("-" * W + "\n")
        f.write(f"\nTotal regions: {len(ordered_all)}\n")

    print(f"✓ Ordered results (coding/non-coding separate): {ordered_results_file}")
    print(f"✓ Ordered results (all combined)              : {combined_ordered_file}")

    # ======================================================================
    # R-compatible output files  (Panel A: genes | Panel B: non-coding)
    # Both follow the unified genomic order.
    # IGS uses directional name (trnH-psbA not sorted psbA-trnH).
    # ======================================================================
    print(f"\nGenerating R-compatible data files for plotting...")

    r_genes_file = os.path.join(output_folder, "genes_for_R_plot.txt")
    with open(r_genes_file, "w") as f:
        f.write("Region\tValue\n")
        for rpos, rtype, display, rname, pi, vs, ln, ns in ordered_all:
            if rtype == "Gene":
                f.write(f"{display}\t{pi:.5f}\n")

    r_noncoding_file = os.path.join(output_folder, "noncoding_for_R_plot.txt")
    with open(r_noncoding_file, "w") as f:
        f.write("Region\tValue\n")
        for rpos, rtype, display, rname, pi, vs, ln, ns in ordered_all:
            if rtype != "Gene":
                f.write(f"{display}\t{pi:.5f}\n")

    print(f"✓ R-compatible genes data     : {r_genes_file}")
    print(f"✓ R-compatible non-coding data: {r_noncoding_file}")

    generate_r_visualizations(r_genes_file, r_noncoding_file, output_folder)

    print(f"{'='*80}")


# ============================================================================
# R VISUALIZATION FUNCTIONS
# ============================================================================

def check_r_packages() -> tuple[bool, list[str]]:
    """Check if required R packages are installed."""
    required_packages = ['ggplot2', 'readr', 'cowplot', 'scales']
    
    try:
        r_code = f'''
packages <- c({", ".join([f'"{pkg}"' for pkg in required_packages])})
missing <- packages[!packages %in% installed.packages()[,"Package"]]
cat(paste(missing, collapse=","))
'''
        result = subprocess.run(['R', '--vanilla', '--slave', '-e', r_code],
                              capture_output=True,
                              text=True,
                              timeout=10)
        
        missing = [pkg.strip() for pkg in result.stdout.strip().split(',') if pkg.strip()]
        return (len(missing) == 0, missing)
    except Exception:
        return (False, required_packages)


def install_r_packages(packages: list[str]) -> bool:
    """Attempt to install missing R packages."""
    print(f"\n  → Installing R packages: {', '.join(packages)}")
    try:
        r_code = f'''
packages <- c({", ".join([f'"{pkg}"' for pkg in packages])})
install.packages(packages, repos="https://cloud.r-project.org", quiet=TRUE)
'''
        result = subprocess.run(['R', '--vanilla', '--slave', '-e', r_code],
                              capture_output=True,
                              text=True,
                              timeout=300)
        return result.returncode == 0
    except Exception as e:
        print(f"  ✗ Installation failed: {e}")
        return False


def generate_r_visualizations(genes_file, noncoding_file, output_folder):
    """
    Generate R-based nucleotide diversity visualizations.
    
    Parameters:
    -----------
    genes_file : str
        Path to genes data file for R
    noncoding_file : str
        Path to non-coding data file for R
    output_folder : str
        Output directory
    """
    import subprocess
    
    print(f"\n{'='*80}")
    print("GENERATING R VISUALIZATIONS")
    print(f"{'='*80}")
    
    # Check R availability
    try:
        result = subprocess.run(['Rscript', '--version'],
                              capture_output=True,
                              text=True,
                              timeout=10)
        if result.returncode != 0:
            print("  ⚠ R not found. Skipping visualization.")
            print("  Install R to enable figure generation.")
            return
    except (subprocess.TimeoutExpired, FileNotFoundError):
        print("  ⚠ R not found. Skipping visualization.")
        print("  Install R from: https://www.r-project.org/")
        return
    
    print("\n  ✓ R is available")
    
    # Check R packages
    packages_ok, missing = check_r_packages()
    if not packages_ok:
        print(f"\n  ⚠ Missing R packages: {', '.join(missing)}")
        print("  → Attempting automatic installation...")
        if install_r_packages(missing):
            print("  ✓ R packages installed successfully")
        else:
            print("  ✗ Could not install R packages automatically")
            print("  → Please install manually in R:")
            print(f"     install.packages(c({', '.join([repr(p) for p in missing])}))")
            return
    else:
        print("  ✓ All required R packages are installed")
    
    # Create figures directory
    figures_dir = os.path.join(output_folder, "Figures")
    os.makedirs(figures_dir, exist_ok=True)
    
    # Create R script
    r_script_content = create_diversity_r_script(
        genes_file,
        noncoding_file
    )
    
    r_script_file = os.path.join(output_folder, "generate_diversity_plot.R")
    with open(r_script_file, 'w') as f:
        f.write(r_script_content)
    
    print(f"  ✓ R script created: {os.path.basename(r_script_file)}")
    
    # Execute R script
    print(f"\n  Executing R script...")
    try:
        # Use basename since cwd is set to output_folder
        result = subprocess.run(['Rscript', os.path.basename(r_script_file)],
                              capture_output=True,
                              text=True,
                              timeout=180,
                              cwd=output_folder)
        
        if result.stdout:
            print(result.stdout)
        
        if result.returncode == 0:
            print(f"\n  ✓ R visualization completed successfully")
            print(f"\n  Generated figure in Figures/:")
            print(f"    - nucleotide_diversity_plot.pdf")
            print(f"    - nucleotide_diversity_plot.png")
            
            # Verify files
            import glob
            pdf_files = glob.glob(os.path.join(figures_dir, "*.pdf"))
            png_files = glob.glob(os.path.join(figures_dir, "*.png"))
            if pdf_files or png_files:
                print(f"\n  Files created:")
                for f in sorted(pdf_files + png_files):
                    size = os.path.getsize(f)
                    print(f"    ✓ {os.path.basename(f)} ({size:,} bytes)")
        else:
            print(f"  ✗ R script execution failed (return code: {result.returncode})")
            if result.stderr:
                print(f"\n  Error output:")
                print(result.stderr[:500])
    
    except subprocess.TimeoutExpired:
        print(f"  ✗ R script timed out (>180s)")
    except Exception as e:
        print(f"  ✗ Error executing R script: {e}")


def create_diversity_r_script(genes_file, noncoding_file):
    """
    Create R script for nucleotide diversity visualization.
    
    Parameters:
    -----------
    genes_file : str
        Path to genes data file
    noncoding_file : str
        Path to non-coding regions data file
    """
    
    # Base R script
    script = f'''
# CGAS Module 13: Nucleotide Diversity Visualization
# Generated automatically

suppressPackageStartupMessages({{
  library(ggplot2)
  library(readr)
  library(cowplot)
  library(scales)
}})

cat("\\n========================================\\n")
cat("CGAS MODULE 13: DIVERSITY VISUALIZATION\\n")
cat("========================================\\n\\n")

# Read data
cat("Reading nucleotide diversity data...\\n")
df_A <- read_delim("{os.path.basename(genes_file)}", delim = "\\t", col_names = TRUE, show_col_types = FALSE)
df_B <- read_delim("{os.path.basename(noncoding_file)}", delim = "\\t", col_names = TRUE, show_col_types = FALSE)

cat(paste("  ✓ Loaded genes:", nrow(df_A), "regions\\n"))
cat(paste("  ✓ Loaded non-coding:", nrow(df_B), "regions\\n\\n"))

# Smart y-axis interval function
# Determines appropriate interval based on max value
# Each panel calculates independently since genes and non-coding have different ranges
get_y_axis_breaks <- function(max_val) {{
  if (max_val <= 0.04) {{
    interval <- 0.005
  }} else if (max_val <= 0.08) {{
    interval <- 0.01
  }} else if (max_val <= 0.12) {{
    interval <- 0.02
  }} else if (max_val <= 0.3) {{
    interval <- 0.05
  }} else if (max_val <= 0.6) {{
    interval <- 0.1
  }} else if (max_val <= 1.2) {{
    interval <- 0.2
  }} else if (max_val <= 2.0) {{
    interval <- 0.4
  }} else {{
    interval <- 0.5
  }}
  max_break <- ceiling(max_val / interval) * interval
  return(list(breaks = seq(0, max_break, by = interval), interval = interval))
}}

# Calculate optimal y-axis breaks for both panels
y_axis_A <- get_y_axis_breaks(max(df_A$Value, na.rm = TRUE))
y_axis_B <- get_y_axis_breaks(max(df_B$Value, na.rm = TRUE))

cat(sprintf("Y-axis intervals:\\n"))
cat(sprintf("  Panel A: max=%.4f, using %.3f intervals\\n", 
            max(df_A$Value, na.rm = TRUE), y_axis_A$interval))
cat(sprintf("  Panel B: max=%.4f, using %.3f intervals\\n\\n", 
            max(df_B$Value, na.rm = TRUE), y_axis_B$interval))

# Create figures directory
if (!dir.exists("Figures")) {{
  dir.create("Figures")
}}

# ============================================================================
# PANEL A: GENES  (blue line + points, your original fine scaling)
# ============================================================================

cat("Generating Panel A: Genes...\\n")

plot_A <- ggplot(df_A, aes(x = factor(Region, levels = Region), y = Value, group = 1)) +
  geom_line(color = "#2c7fb8", size = 0.5) +
  geom_point(size = 1.5, color = "#2c7fb8", fill = "white",
             shape = 21, stroke = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "A: Genes", x = NULL,
       y = expression(paste("Nucleotide diversity (", pi, ")"))) +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      size = 9, face = "italic",
      margin = margin(t = 2, b = 20)
    ),
    axis.text.y  = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 10, margin = margin(r = 5)),
    axis.line    = element_line(colour = "black", linewidth = 0.4),
    axis.ticks   = element_line(colour = "black", linewidth = 0.4),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor   = element_blank(),
    plot.title = element_text(
      size = 11, face = "bold", hjust = 0,
      margin = margin(b = 10)
    ),
    plot.margin = margin(t = 10, r = 10, b = 50, l = 10)
  ) +
  scale_y_continuous(
    breaks = y_axis_A$breaks,
    limits = c(0, max(df_A$Value, na.rm = TRUE) * 1.1),
    labels = scales::label_number(accuracy = y_axis_A$interval),
    expand = expansion(mult = c(0.05, 0.05))
  )

# ============================================================================
# PANEL B: INTERGENIC SPACERS AND INTRONS
# ============================================================================

cat("Generating Panel B: Non-coding regions...\\n")

plot_B <- ggplot(df_B, aes(x = factor(Region, levels = Region), y = Value, group = 1)) +
  geom_line(color = "#2c7fb8", size = 0.5) +
  geom_point(size = 1.5, color = "#2c7fb8", fill = "white", shape = 21, stroke = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "B: Intergenic spacers and introns", x = NULL,
       y = expression(paste("Nucleotide diversity (", pi, ")"))) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                               size = 9, face = "italic",
                               margin = margin(t = 2, b = 20)),
    axis.text.y  = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 10, margin = margin(r = 5)),
    axis.line    = element_line(colour = "black", linewidth = 0.4),
    axis.ticks   = element_line(colour = "black", linewidth = 0.4),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor   = element_blank(),
    plot.title = element_text(size = 11, face = "bold", hjust = 0, 
                              margin = margin(b = 10)),
    plot.margin = margin(t = 10, r = 10, b = 50, l = 10)
  ) +
  scale_y_continuous(breaks = y_axis_B$breaks,
                     limits = c(0, max(df_B$Value, na.rm = TRUE) * 1.1),
                     labels = scales::number_format(accuracy = y_axis_B$interval, decimal.mark = "."),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(clip = "off")

# ============================================================================
# COMBINE PANELS AND SAVE
# ============================================================================

cat("Combining panels and saving...\\n")

final_plot <- plot_grid(plot_A, plot_B, 
                        ncol = 1, 
                        rel_heights = c(1, 1),
                        align = "v",
                        axis  = "lr")

# Save as PDF
ggsave("Figures/nucleotide_diversity_plot.pdf",
       final_plot,
       width = 20,
       height = 10,
       units = "in",
       device = cairo_pdf,
       dpi = 600)

# Save as PNG
ggsave("Figures/nucleotide_diversity_plot.png",
       final_plot,
       width = 20,
       height = 10,
       units = "in",
       dpi = 600,
       bg = "white")

cat("\\n========================================\\n")
cat("Visualization completed successfully!\\n")
cat("========================================\\n")
cat("\\nOutput files:\\n")
cat("  - Figures/nucleotide_diversity_plot.pdf\\n")
cat("  - Figures/nucleotide_diversity_plot.png\\n")
'''
    
    return script


if __name__ == "__main__":
    main()