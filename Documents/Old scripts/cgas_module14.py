#!/usr/bin/env python3
"""
CGAS Module 14: Phylogenetic Matrix Builder and Tree Inference
===============================================================

Extracts genomic regions from GenBank files, aligns them, creates concatenated
matrices, and optionally performs phylogenetic inference with IQ-TREE.

Matrix Types Generated:
1. Protein-coding genes (CDS) only
2. Introns only
3. Intergenic spacer regions only
4. Complete matrix (all regions) - RECOMMENDED for phylogeny

Alignment Methods:
- Default: MAFFT (fast, accurate for all region types)
- --macse: MACSE for CDS only (codon-aware translation alignment)

Phylogenetic Inference (NEW):
- Automatic IQ-TREE integration with --iqtree flag
- Model selection: ModelFinder Plus (MFP)
- Support values: 1000 UFBoot + 1000 SH-aLRT
- Partitioned analysis with automatically generated partition files
- Multi-threading support (auto-detects CPU cores)

Features:
- Uses organism names from GenBank files (not accessions)
- Region-specific extraction: --genes-only, --introns-only, --intergenic-only, --complete-only
- Optional outgroup specification: -og "Species name"
- CDS validation and gene name normalization
- Partition files in 4 formats: RAxML, NEXUS, IQ-TREE, PHYLIP
- Compatible with all IQ-TREE versions (iqtree, iqtree2, iqtree-omp)

Part of CGAS (Chloroplast Genome Analysis Suite)
Author: Abdullah
Version: 2.0 - Module 14 with IQ-TREE integration
Date: January 2026

Dependencies:
    Python: biopython
    Required: MAFFT
    Optional: MACSE (for codon-aware alignment), IQ-TREE (for phylogeny)

Usage Examples:
    # Basic matrix generation
    python cgas_module14.py
    
    # With phylogenetic inference
    python cgas_module14.py --iqtree
    
    # Use MACSE + phylogeny
    python cgas_module14.py --macse --iqtree
    
    # Genes only with custom output
    python cgas_module14.py --genes-only --iqtree -o results/
    
    # Complete matrix only
    python cgas_module14.py --complete-only
    
    # Complete matrix with phylogeny
    python cgas_module14.py --complete-only --iqtree
    
    # Specify outgroup
    python cgas_module14.py --iqtree -og "Gossypium arboreum"
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict, OrderedDict
import re
import datetime


# ============================================================================
# CONSTANTS
# ============================================================================

MODULE_NAME = "Module14_Phylogeny"
MODULE_VERSION = "2.0"


def get_name_from_filename(gb_file):
    """
    Extract a unique identifier from the GenBank filename.
    
    This function uses the filename (without extension) as the identifier
    to handle cases where multiple accessions of the same species exist.
    
    Examples:
        Abutilon_theophrasti_MN641467.gb -> Abutilon_theophrasti_MN641467
        Abutilon_theophrasti_NC_053702.gb -> Abutilon_theophrasti_NC_053702
        species_name.gb -> species_name
    
    Args:
        gb_file: Path to GenBank file
    
    Returns:
        String with the filename (without path and extension)
    """
    import os
    # Get filename without directory path
    filename = os.path.basename(gb_file)
    # Remove extension (.gb, .gbk, .genbank, etc.)
    name = os.path.splitext(filename)[0]
    # Clean up: replace spaces and special characters
    name = name.replace(' ', '_').replace(',', '').replace(';', '').replace(':', '').replace('(', '').replace(')', '')
    return name


def get_organism_name_from_genbank(gb_file):
    """
    Extract the organism name from a GenBank file.
    
    NOTE: This function is kept for backward compatibility but is NO LONGER
    used for sequence naming. Use get_name_from_filename() instead to handle
    multiple accessions of the same species.
    
    Args:
        gb_file: Path to GenBank file
    
    Returns:
        String with the organism name (spaces replaced with underscores) or None if file cannot be read
    """
    try:
        for record in SeqIO.parse(gb_file, "genbank"):
            # Extract organism name from annotations
            organism = record.annotations.get('organism', 'unknown')
            # If no organism found in annotations, try source feature
            if organism == 'unknown':
                for feature in record.features:
                    if feature.type == "source":
                        if "organism" in feature.qualifiers:
                            organism = feature.qualifiers["organism"][0]
                            break
            # Use accession if organism still not found
            if organism == 'unknown':
                organism = record.id
            # Clean up organism name (replace spaces with underscores, remove special chars)
            organism = organism.replace(' ', '_').replace(',', '').replace(';', '').replace(':', '')
            return organism
    except Exception as e:
        print(f"Warning: Could not read organism name from {gb_file}: {e}")
        return None


def parse_outgroup_files(outgroup_args, genbank_files):
    """
    Parse outgroup file arguments and get their identifiers.
    
    Uses filename-based naming to handle multiple accessions of the same species.
    
    Args:
        outgroup_args: List of outgroup file names (can be full paths or just filenames)
        genbank_files: List of all GenBank file paths
    
    Returns:
        List of identifiers for outgroup species (based on filenames)
    """
    if not outgroup_args:
        return []
    
    outgroup_names = []
    
    # Create a mapping of filenames to full paths
    filename_to_path = {}
    for gb_file in genbank_files:
        filename = os.path.basename(str(gb_file))
        filename_to_path[filename] = gb_file
    
    for outgroup_arg in outgroup_args:
        # Check if it's a full path
        if os.path.exists(outgroup_arg):
            identifier = get_name_from_filename(outgroup_arg)
            if identifier:
                outgroup_names.append(identifier)
                print(f"  Outgroup: {identifier} (from {os.path.basename(outgroup_arg)})")
            else:
                print(f"  Warning: Could not get identifier from {outgroup_arg}")
        # Check if it's just a filename
        elif outgroup_arg in filename_to_path:
            identifier = get_name_from_filename(filename_to_path[outgroup_arg])
            if identifier:
                outgroup_names.append(identifier)
                print(f"  Outgroup: {identifier} (from {outgroup_arg})")
            else:
                print(f"  Warning: Could not get identifier from {outgroup_arg}")
        else:
            print(f"  Warning: Outgroup file not found: {outgroup_arg}")
    
    return outgroup_names


def normalize_gene_name(gene_name, log_replacements=None):
    """
    Normalize gene names to handle different naming conventions.
    Handles tRNA naming (anticodon variations), gene synonyms, and case normalization.
    
    Args:
        gene_name: Original gene name
        log_replacements: Optional dict to log synonym replacements
    
    Returns:
        Normalized gene name
    """
    if not gene_name:
        return gene_name
    
    original_name = gene_name
    
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
        # inf gene
        'infa': 'infA',
        # mat gene
        'matk': 'matK',
        # cem gene
        'cema': 'cemA',
        # acc gene
        'accd': 'accD',
    }
    
    # Handle tRNA naming conventions
    if gene_name.lower().startswith('trn'):
        match = re.match(r'(trn[A-Z])[_\-\(]?([A-Z]{3})\)?', gene_name, re.IGNORECASE)
        if match:
            amino_acid = 'trn' + match.group(1)[3].upper()
            anticodon = match.group(2).upper()
            normalized = f"{amino_acid}-{anticodon}"
            if log_replacements is not None and normalized != original_name:
                log_replacements[original_name] = normalized
            return normalized
        else:
            match = re.match(r'(trn[A-Z])', gene_name, re.IGNORECASE)
            if match:
                amino_acid = match.group(1)
                normalized = 'trn' + amino_acid[3].upper()
                if log_replacements is not None and normalized != original_name:
                    log_replacements[original_name] = normalized
                return normalized
    
    # Handle rRNA naming conventions
    if gene_name.lower().startswith('rrn'):
        normalized = gene_name.lower()
        if log_replacements is not None and normalized != original_name:
            log_replacements[original_name] = normalized
        return normalized
    
    # Check for gene synonyms
    gene_lower = gene_name.lower()
    if gene_lower in gene_synonyms:
        normalized = gene_synonyms[gene_lower]
        if log_replacements is not None:
            log_replacements[original_name] = normalized
        return normalized
    
    return gene_name


def reorder_alignment_with_outgroups(alignment, outgroup_names):
    """
    Reorder alignment to place outgroup sequences at the top.
    
    Args:
        alignment: MultipleSeqAlignment object
        outgroup_names: List of organism names to place at top
    
    Returns:
        Reordered MultipleSeqAlignment object
    """
    if not alignment or not outgroup_names:
        return alignment
    
    # Separate outgroup and ingroup sequences
    outgroup_seqs = []
    ingroup_seqs = []
    
    for record in alignment:
        if record.id in outgroup_names:
            outgroup_seqs.append(record)
        else:
            ingroup_seqs.append(record)
    
    # Sort outgroups according to the order specified
    sorted_outgroups = []
    for og_name in outgroup_names:
        for record in outgroup_seqs:
            if record.id == og_name:
                sorted_outgroups.append(record)
                break
    
    # Combine: outgroups first, then ingroups
    reordered_records = sorted_outgroups + ingroup_seqs
    
    # Create new alignment
    if reordered_records:
        return MultipleSeqAlignment(reordered_records)
    else:
        return alignment


def extract_features_from_genbank(genbank_files):
    """
    Extract genes, introns, and intergenic regions from GenBank files.
    
    Returns:
        protein_coding_genes: {gene_name: [SeqRecord1, SeqRecord2, ...]}
        introns: {intron_name: [SeqRecord1, SeqRecord2, ...]}
        intergenic: {spacer_name: [SeqRecord1, SeqRecord2, ...]}
        synonym_replacements: {original_name: normalized_name}
    """
    protein_coding_genes = defaultdict(list)
    introns = defaultdict(list)
    intergenic = defaultdict(list)
    synonym_replacements = {}
    
    for gb_file in genbank_files:
        print(f"Processing {os.path.basename(gb_file)}...")
        
        # Use filename-based identifier instead of organism name
        # This handles multiple accessions of the same species
        organism_id = get_name_from_filename(gb_file)
        
        for record in SeqIO.parse(gb_file, "genbank"):
            sequence = record.seq
            
            # Extract all relevant features
            all_features = [f for f in record.features if f.type in ["gene", "CDS", "tRNA", "rRNA"]]
            
            # Group features by gene
            features_by_gene = defaultdict(list)
            
            for feature in all_features:
                if feature.type not in ["CDS", "tRNA", "rRNA"]:
                    continue
                
                # Get gene name
                gene_name = None
                if "gene" in feature.qualifiers:
                    gene_name = feature.qualifiers["gene"][0]
                elif "product" in feature.qualifiers:
                    gene_name = feature.qualifiers["product"][0]
                elif "locus_tag" in feature.qualifiers:
                    gene_name = feature.qualifiers["locus_tag"][0]
                else:
                    gene_name = f"unknown_{feature.location.start}"
                
                # Normalize gene name and log any replacements
                gene_name = normalize_gene_name(gene_name, synonym_replacements)
                
                # Group by gene
                features_by_gene[gene_name].append(feature)
            
            # Process each gene
            for gene_name, features in features_by_gene.items():
                # Sort features by position
                features.sort(key=lambda f: f.location.start)
                
                # Check if gene has introns (multiple CDS features or join location)
                has_introns = False
                exon_locations = []
                
                for feature in features:
                    if feature.type == "CDS":
                        if hasattr(feature.location, 'parts') and len(feature.location.parts) > 1:
                            has_introns = True
                            exon_locations = list(feature.location.parts)
                            break
                
                if not has_introns and len([f for f in features if f.type == "CDS"]) > 1:
                    has_introns = True
                    exon_locations = [f.location for f in features if f.type == "CDS"]
                
                # Extract gene sequence
                gene_seq = None
                is_protein_coding = False
                
                for feature in features:
                    if feature.type == "CDS":
                        is_protein_coding = True
                        # CRITICAL FIX: Use feature.extract() which handles joined locations correctly
                        # This properly concatenates all exons for multi-exon CDS
                        try:
                            gene_seq = feature.extract(sequence)
                        except Exception as e:
                            print(f"    Warning: Failed to extract CDS for {gene_name}: {e}")
                            gene_seq = None
                        break
                    elif feature.type in ["tRNA", "rRNA"]:
                        try:
                            gene_seq = feature.extract(sequence)
                        except Exception as e:
                            print(f"    Warning: Failed to extract {feature.type} for {gene_name}: {e}")
                            gene_seq = None
                        break
                
                if gene_seq:
                    # Validate CDS sequences
                    if is_protein_coding:
                        # Check if length is multiple of 3 (valid reading frame)
                        if len(gene_seq) % 3 != 0:
                            print(f"    Warning: {gene_name} length {len(gene_seq)} not multiple of 3 (frameshift?)")
                    
                    # FASTA header format for Geneious compatibility:
                    # Use id field only, set description to empty string
                    # This ensures Geneious displays the full name correctly
                    seq_record = SeqRecord(
                        gene_seq,
                        id=organism_id,
                        description=""
                    )
                    
                    # Add to appropriate dictionaries
                    if is_protein_coding:
                        protein_coding_genes[gene_name].append(seq_record)
                
                # Extract introns if present
                if has_introns and exon_locations:
                    for i in range(len(exon_locations) - 1):
                        exon1 = exon_locations[i]
                        exon2 = exon_locations[i + 1]
                        
                        # Get intron boundaries
                        if exon1.strand == 1:
                            intron_start = exon1.end
                            intron_end = exon2.start
                        else:
                            intron_start = exon2.end
                            intron_end = exon1.start
                        
                        if intron_start < intron_end:
                            intron_seq = sequence[intron_start:intron_end]
                            
                            if exon1.strand == -1:
                                intron_seq = intron_seq.reverse_complement()
                            
                            intron_name = f"{gene_name}_intron{i+1}"
                            
                            intron_record = SeqRecord(
                                intron_seq,
                                id=organism_id,
                                description=""
                            )
                            
                            introns[intron_name].append(intron_record)
            
            # Extract intergenic regions
            # NOTE: Intergenic spacers are extracted in genomic forward orientation
            # regardless of flanking gene strand. This is a biological assumption.
            # For strand-specific analysis, filter spacers manually.
            
            # Get all gene boundaries
            gene_boundaries = []
            gene_features_list = []
            for feature in all_features:
                if feature.type in ["CDS", "tRNA", "rRNA"]:
                    gene_boundaries.append((feature.location.start, feature.location.end))
                    gene_features_list.append(feature)
            
            # Sort boundaries
            gene_boundaries.sort()
            
            # Extract intergenic spacers
            for i in range(len(gene_boundaries) - 1):
                end_gene1 = gene_boundaries[i][1]
                start_gene2 = gene_boundaries[i + 1][0]
                
                # Check if there's a gap (and not overlapping genes)
                if start_gene2 > end_gene1:
                    spacer_length = start_gene2 - end_gene1
                    
                    # Skip very short spacers (likely annotation artifacts)
                    if spacer_length < 10:
                        continue
                    
                    spacer_seq = sequence[end_gene1:start_gene2]
                    
                    # Create spacer name based on flanking genes
                    # Find the genes
                    gene1_name = None
                    gene2_name = None
                    
                    for feature in all_features:
                        if feature.type in ["CDS", "tRNA", "rRNA"]:
                            if feature.location.end == end_gene1:
                                if "gene" in feature.qualifiers:
                                    gene1_name = normalize_gene_name(feature.qualifiers["gene"][0])
                                elif "product" in feature.qualifiers:
                                    gene1_name = normalize_gene_name(feature.qualifiers["product"][0])
                            
                            if feature.location.start == start_gene2:
                                if "gene" in feature.qualifiers:
                                    gene2_name = normalize_gene_name(feature.qualifiers["gene"][0])
                                elif "product" in feature.qualifiers:
                                    gene2_name = normalize_gene_name(feature.qualifiers["product"][0])
                    
                    if gene1_name and gene2_name:
                        spacer_name = f"{gene1_name}_{gene2_name}"
                        
                        spacer_record = SeqRecord(
                            spacer_seq,
                            id=organism_id,
                            description=""
                        )
                        
                        intergenic[spacer_name].append(spacer_record)
    
    return protein_coding_genes, introns, intergenic, synonym_replacements


def check_mafft():
    """Check if MAFFT is available."""
    try:
        subprocess.run(["mafft", "--version"], 
                      stdout=subprocess.PIPE, 
                      stderr=subprocess.PIPE, 
                      check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def check_macse():
    """Check if MACSE is available."""
    try:
        # Try to run MACSE
        result = subprocess.run(["java", "-jar", "macse.jar", "-help"], 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE,
                              timeout=5)
        return True
    except:
        # Try alternative location
        try:
            result = subprocess.run(["macse", "-help"], 
                                  stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE,
                                  timeout=5)
            return True
        except:
            return False


def find_macse_path():
    """Find the path to MACSE jar file or executable."""
    # Check common locations
    common_paths = [
        "macse.jar",
        "macse_v2.jar",
        "/usr/local/bin/macse.jar",
        os.path.expanduser("~/bin/macse.jar"),
        os.path.expanduser("~/tools/macse.jar"),
    ]
    
    for path in common_paths:
        if os.path.exists(path):
            return path
    
    # Check if macse is in PATH as executable
    try:
        subprocess.run(["macse", "-help"], 
                      stdout=subprocess.PIPE, 
                      stderr=subprocess.PIPE,
                      timeout=5)
        return "macse"
    except:
        pass
    
    return None


def align_sequences_macse(seq_records, output_fasta, output_alignment, macse_path, outgroup_names=None):
    """
    Align protein-coding sequences using MACSE (codon-aware alignment).
    This creates a separate MACSE alignment file in addition to the MAFFT alignment.
    
    Args:
        seq_records: List of SeqRecord objects
        output_fasta: Path to unaligned FASTA (already created by MAFFT)
        output_alignment: Path to write MACSE aligned FASTA
        macse_path: Path to MACSE jar or executable
        outgroup_names: List of organism names to place at top of alignment
    
    Returns:
        MultipleSeqAlignment object or None if alignment fails
    """
    # Skip if less than 2 sequences
    if len(seq_records) < 2:
        return None
    
    # CRITICAL: Write unaligned sequences FIRST
    SeqIO.write(seq_records, output_fasta, "fasta")
    
    # Prepare MACSE output files
    macse_nt_output = output_alignment.replace(".fasta", "_NT.fasta")
    macse_aa_output = output_alignment.replace(".fasta", "_AA.fasta")
    
    try:
        # Run MACSE
        if macse_path.endswith('.jar'):
            cmd = [
                "java", "-jar", macse_path,
                "-prog", "alignSequences",
                "-seq", output_fasta,
                "-out_NT", macse_nt_output,
                "-out_AA", macse_aa_output
            ]
        else:
            cmd = [
                macse_path,
                "-prog", "alignSequences",
                "-seq", output_fasta,
                "-out_NT", macse_nt_output,
                "-out_AA", macse_aa_output
            ]
        
        subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            timeout=300  # 5 minute timeout
        )
        
        # Read MACSE nucleotide alignment
        alignment = AlignIO.read(macse_nt_output, "fasta")
        
        # MACSE uses ! for frameshifts and * for stops - convert to gaps for compatibility
        cleaned_records = []
        for record in alignment:
            cleaned_seq = str(record.seq).replace('!', '-').replace('*', '-')
            cleaned_record = SeqRecord(
                Seq(cleaned_seq),
                id=record.id,
                description=""  # Empty description for Geneious compatibility
            )
            cleaned_records.append(cleaned_record)
        
        alignment = MultipleSeqAlignment(cleaned_records)
        
        # Reorder alignment if outgroups specified
        if outgroup_names:
            alignment = reorder_alignment_with_outgroups(alignment, outgroup_names)
        
        # Write cleaned alignment
        AlignIO.write(alignment, output_alignment, "fasta")
        
        # Clean up intermediate files
        if os.path.exists(macse_aa_output):
            os.remove(macse_aa_output)
        if os.path.exists(macse_nt_output):
            os.remove(macse_nt_output)
        
        return alignment
    
    except subprocess.CalledProcessError as e:
        # MACSE command failed - show the error
        error_msg = e.stderr.decode() if e.stderr else str(e)
        print(f" (MACSE error: {error_msg[:100]})", end="")
        return None
    except subprocess.TimeoutExpired as e:
        print(f" (MACSE timeout)", end="")
        return None
    except Exception as e:
        print(f" (Error: {str(e)[:100]})", end="")
        return None


def align_sequences(seq_records, output_fasta, output_alignment, outgroup_names=None):
    """
    Align sequences using MAFFT.
    
    Args:
        seq_records: List of SeqRecord objects
        output_fasta: Path to write unaligned FASTA
        output_alignment: Path to write aligned FASTA
        outgroup_names: List of organism names to place at top of alignment
    
    Returns:
        MultipleSeqAlignment object or None if alignment fails
    """
    # Skip if less than 2 sequences
    if len(seq_records) < 2:
        print(f"  Skipping (only {len(seq_records)} sequence)")
        return None
    
    # Write unaligned sequences
    SeqIO.write(seq_records, output_fasta, "fasta")
    
    # Run MAFFT
    try:
        with open(output_alignment, "w") as out_file:
            subprocess.run(
                ["mafft", "--auto", "--quiet", output_fasta],
                stdout=out_file,
                stderr=subprocess.PIPE,
                check=True
            )
        
        # Read alignment
        alignment = AlignIO.read(output_alignment, "fasta")
        
        # Reorder alignment if outgroups specified
        if outgroup_names:
            alignment = reorder_alignment_with_outgroups(alignment, outgroup_names)
            # Write reordered alignment
            AlignIO.write(alignment, output_alignment, "fasta")
        
        return alignment
    
    except subprocess.CalledProcessError as e:
        print(f"  MAFFT alignment failed: {e}")
        return None


def concatenate_alignments(alignment_dict, region_order, output_file, matrix_type="all", outgroup_names=None):
    """
    Concatenate multiple alignments into a single matrix.
    
    Args:
        alignment_dict: {region_name: alignment_object}
        region_order: List of region names in desired order
        output_file: Path to output concatenated alignment
        matrix_type: Type of matrix being created (for reporting)
        outgroup_names: List of organism names to place at top of matrix
    
    Returns:
        Tuple of (alignment length, number of regions, partition info)
    """
    if not alignment_dict:
        print(f"  No alignments to concatenate for {matrix_type}")
        return 0, 0, []
    
    # Get all sequence IDs (organism names)
    all_ids = set()
    for alignment in alignment_dict.values():
        for record in alignment:
            all_ids.add(record.id)
    
    all_ids = sorted(list(all_ids))
    
    # Reorder to place outgroups first
    if outgroup_names:
        outgroup_present = [og_name for og_name in outgroup_names if og_name in all_ids]
        ingroup_ids = [seq_id for seq_id in all_ids if seq_id not in outgroup_names]
        all_ids = outgroup_present + sorted(ingroup_ids)
    
    # Create concatenated sequences
    concatenated_seqs = {seq_id: [] for seq_id in all_ids}
    partition_info = []
    current_position = 1
    
    # Process regions in order
    for region_name in region_order:
        if region_name not in alignment_dict:
            continue
        
        alignment = alignment_dict[region_name]
        alignment_length = alignment.get_alignment_length()
        
        # Record partition information
        end_position = current_position + alignment_length - 1
        partition_info.append((region_name, current_position, end_position))
        
        # Add sequences for each organism
        for seq_id in all_ids:
            # Find sequence in alignment
            seq_found = False
            for record in alignment:
                if record.id == seq_id:
                    concatenated_seqs[seq_id].append(str(record.seq))
                    seq_found = True
                    break
            
            # If sequence not found, add gaps
            if not seq_found:
                concatenated_seqs[seq_id].append("-" * alignment_length)
        
        current_position = end_position + 1
    
    # Write concatenated alignment
    total_length = current_position - 1
    
    with open(output_file, "w") as f:
        for seq_id in all_ids:
            concatenated_seq = "".join(concatenated_seqs[seq_id])
            f.write(f">{seq_id}\n")
            f.write(f"{concatenated_seq}\n")
    
    return total_length, len(partition_info), partition_info


def write_partition_file(partition_info, output_file, format_type="raxml"):
    """
    Write partition file for phylogenetic analysis in multiple formats.
    
    Args:
        partition_info: List of tuples (region_name, start, end)
        output_file: Path to output partition file
        format_type: 'raxml', 'nexus', 'iqtree', or 'phylip'
    """
    with open(output_file, "w") as f:
        if format_type == "raxml":
            # RAxML format
            for region_name, start, end in partition_info:
                f.write(f"DNA, {region_name} = {start}-{end}\n")
        
        elif format_type == "nexus":
            # NEXUS format (for MrBayes, PAUP*)
            f.write("#nexus\n")
            f.write("begin sets;\n")
            for region_name, start, end in partition_info:
                f.write(f"  charset {region_name} = {start}-{end};\n")
            f.write("end;\n")
        
        elif format_type == "iqtree":
            # IQ-TREE format (NEXUS-like but slightly different)
            f.write("#nexus\n")
            f.write("begin sets;\n")
            for region_name, start, end in partition_info:
                # Replace problematic characters in region names for IQ-TREE
                clean_name = region_name.replace('-', '_').replace('.', '_')
                f.write(f"  charset {clean_name} = {start}-{end};\n")
            
            # Add partition scheme with cleaned names
            f.write("  charpartition mine = ")
            partition_defs = [f"GTR+G:{region_name.replace('-', '_').replace('.', '_')}" 
                             for region_name, _, _ in partition_info]
            f.write(", ".join(partition_defs))
            f.write(";\n")
            f.write("end;\n")
        
        elif format_type == "phylip":
            # Simple PHYLIP-style format
            for region_name, start, end in partition_info:
                f.write(f"{region_name}\t{start}\t{end}\n")


# ============================================================================
# IQ-TREE PHYLOGENETIC INFERENCE
# ============================================================================

def check_iqtree_availability():
    """
    Check if IQ-TREE is installed and get version info.
    Tries multiple possible executable names to support different installations.
    
    Returns:
        tuple: (is_available, version_string, executable_name)
    """
    # Try different IQ-TREE executable names
    possible_names = ['iqtree2', 'iqtree', 'iqtree3', 'iqtree-omp', 'iqtree2-omp']
    
    for exe_name in possible_names:
        try:
            result = subprocess.run([exe_name, '--version'],
                                  capture_output=True,
                                  text=True,
                                  timeout=10)
            if result.returncode == 0:
                # Extract version from first line
                version_line = result.stdout.split('\n')[0] if result.stdout else "IQ-TREE"
                return (True, version_line, exe_name)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    
    return (False, None, None)


def run_iqtree_analysis(alignment_file, output_dir, partition_file=None, 
                        matrix_name="phylogeny", threads='AUTO'):
    """
    Run IQ-TREE phylogenetic inference with recommended settings.
    
    Parameters:
    -----------
    alignment_file : str
        Path to alignment file (FASTA format)
    output_dir : str
        Directory for output files
    partition_file : str, optional
        Path to partition file for partitioned analysis
    matrix_name : str
        Name for output prefix
    threads : str or int
        Number of threads ('AUTO' for automatic, or specific number)
    
    Returns:
    --------
    tuple: (success, tree_file, log_file)
    """
    # Check IQ-TREE availability
    is_available, version, exe_name = check_iqtree_availability()
    
    if not is_available:
        print("\n" + "="*80)
        print("IQ-TREE NOT FOUND")
        print("="*80)
        print("\n  ⚠ IQ-TREE is not installed or not in PATH")
        print("\n  Installation options:")
        print("    • Download from: http://www.iqtree.org/")
        print("    • Via conda: conda install -c bioconda iqtree")
        print("    • Via apt: sudo apt install iqtree")
        print("\n  Skipping phylogenetic inference...")
        print("="*80 + "\n")
        return (False, None, None)
    
    print("\n" + "="*80)
    print("IQ-TREE PHYLOGENETIC INFERENCE")
    print("="*80)
    print(f"\n  ✓ Found: {version}")
    print(f"  ✓ Executable: {exe_name}")
    
    # Set up output prefix
    output_prefix = os.path.join(output_dir, matrix_name)
    
    # Build IQ-TREE command with recommended settings
    cmd = [
        exe_name,
        '-s', alignment_file,           # Input alignment
        '--prefix', output_prefix,       # Output prefix
        '-m', 'MFP',                    # ModelFinder Plus (automatic model selection)
        '-bb', '1000',                  # Ultrafast bootstrap (1000 replicates)
        '-alrt', '1000',                # SH-aLRT test (1000 replicates)
        '-bnni',                        # Optimize UFBoot trees by NNI
        '-nt', str(threads)             # Number of threads
    ]
    
    # Add partition file if provided
    analysis_type = "Partitioned" if partition_file else "Concatenated"
    if partition_file and os.path.exists(partition_file):
        cmd.extend(['-p', partition_file])
        print(f"\n  → Analysis type: {analysis_type}")
        print(f"  → Partition file: {os.path.basename(partition_file)}")
        
        # Validate partition file
        try:
            with open(partition_file, 'r') as f:
                partition_content = f.read()
                if not partition_content.strip():
                    print(f"\n  ⚠ Warning: Partition file is empty, running without partitions")
                    cmd.remove('-p')
                    cmd.remove(partition_file)
                    partition_file = None
        except Exception as e:
            print(f"\n  ⚠ Warning: Could not read partition file: {e}")
            print(f"  Running without partitions...")
            if '-p' in cmd:
                cmd.remove('-p')
            if partition_file in cmd:
                cmd.remove(partition_file)
            partition_file = None
    else:
        print(f"\n  → Analysis type: {analysis_type} (no partitions)")
    
    # Display analysis settings
    print(f"\n  Analysis Settings:")
    print(f"    • Model selection: MFP (ModelFinder Plus)")
    print(f"    • Bootstrap support: 1000 ultrafast bootstrap replicates")
    print(f"    • Branch support: 1000 SH-aLRT replicates")
    print(f"    • BNNI optimization: Enabled")
    print(f"    • Threads: {threads}")
    
    # Run IQ-TREE
    print(f"\n  ⏳ Running IQ-TREE analysis...")
    print(f"     (This may take minutes to hours depending on dataset size)")
    print(f"     Output: {os.path.basename(output_prefix)}.*")
    
    try:
        result = subprocess.run(cmd,
                              capture_output=True,
                              text=True,
                              timeout=7200)  # 2 hour timeout
        
        # Check if analysis completed successfully
        tree_file = f"{output_prefix}.treefile"
        log_file = f"{output_prefix}.log"
        iqtree_file = f"{output_prefix}.iqtree"
        
        if result.returncode == 0 and os.path.exists(tree_file):
            print(f"\n  ✓ IQ-TREE analysis completed successfully!")
            return (True, tree_file, log_file)
        else:
            print(f"\n  ✗ IQ-TREE analysis failed (exit code: {result.returncode})")
            if result.stderr and len(result.stderr) > 0:
                print(f"\n  Error details:")
                # Show more of the error for partition issues
                error_lines = result.stderr.split('\n')
                for line in error_lines:
                    if 'ERROR:' in line or 'CharSet' in line or 'not found' in line:
                        print(f"  {line}")
                        if 'CharSet' in line and 'not found' in line:
                            print(f"\n  💡 Tip: This error usually means there's a mismatch between")
                            print(f"     the partition file and the alignment. Check that all")
                            print(f"     charset names in the partition file match the regions")
                            print(f"     in your alignment.")
            return (False, None, log_file if os.path.exists(log_file) else None)
            
    except subprocess.TimeoutExpired:
        print(f"\n  ✗ IQ-TREE exceeded time limit (2 hours)")
        print(f"  → Consider using a smaller dataset or fewer partitions")
        return (False, None, None)
    except Exception as e:
        print(f"\n  ✗ Error running IQ-TREE: {e}")
        return (False, None, None)


def summarize_iqtree_results(log_file, tree_file, matrix_name):
    """
    Extract and display key results from IQ-TREE log file.
    
    Parameters:
    -----------
    log_file : str
        Path to IQ-TREE .log file
    tree_file : str
        Path to IQ-TREE .treefile
    matrix_name : str
        Name of the matrix analyzed
    """
    if not os.path.exists(log_file):
        return
    
    print(f"\n{'='*80}")
    print(f"IQ-TREE RESULTS: {matrix_name}")
    print(f"{'='*80}")
    
    try:
        with open(log_file, 'r') as f:
            log_content = f.read()
        
        # Extract best model
        if 'Best-fit model' in log_content:
            for line in log_content.split('\n'):
                if 'Best-fit model' in line and ':' in line:
                    print(f"\n  {line.strip()}")
                    break
        
        # Extract log-likelihood
        for line in log_content.split('\n'):
            if 'BEST SCORE FOUND' in line or 'Final LogLikelihood' in line:
                print(f"  {line.strip()}")
                break
        
        # Extract tree stats
        for line in log_content.split('\n'):
            if 'Total tree length' in line or 'Sum of internal branch' in line:
                print(f"  {line.strip()}")
        
        # List output files
        output_dir = os.path.dirname(tree_file)
        prefix = os.path.basename(tree_file).replace('.treefile', '')
        
        print(f"\n  Output Files:")
        print(f"    • {prefix}.treefile      (ML tree with support values)")
        print(f"    • {prefix}.log           (Analysis log)")
        print(f"    • {prefix}.iqtree        (Full IQ-TREE report)")
        print(f"    • {prefix}.model.gz      (Best-fit model parameters)")
        print(f"    • {prefix}.mldist        (ML distances)")
        
        if os.path.exists(f"{output_dir}/{prefix}.contree"):
            print(f"    • {prefix}.contree       (Consensus tree)")
        
    except Exception as e:
        print(f"  ⚠ Could not parse log file: {e}")
    
    print(f"{'='*80}\n")


def main():
    """Main execution function."""
    
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='CGAS Module 14: Phylogenetic Matrix Builder with IQ-TREE Integration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic - uses current directory, processes all regions with MAFFT
  python cgas_module14.py
  
  # With IQ-TREE phylogenetic inference
  python cgas_module14.py --iqtree
  
  # Use MACSE for codon-aware CDS alignment + IQ-TREE
  python cgas_module14.py --macse --iqtree
  
  # Process only protein-coding genes with phylogeny
  python cgas_module14.py --genes-only --iqtree
  
  # Process only introns with phylogeny
  python cgas_module14.py --introns-only --iqtree
  
  # Complete matrix only
  python cgas_module14.py --complete-only
  
  # Complete matrix with phylogeny
  python cgas_module14.py --complete-only --iqtree
  
  # Specify outgroup
  python cgas_module14.py --iqtree -og outgroup.gb
  
  # Custom output directory
  python cgas_module14.py --iqtree -o my_phylogeny/
        """
    )
    
    parser.add_argument('input_dir',
                       nargs='?',
                       default='.',
                       help='Directory containing GenBank files (.gb or .gbk). Default: current directory')
    
    parser.add_argument('-o', '--output', 
                       dest='output_dir',
                       default='Module14_Phylogeny',
                       help='Output directory (default: Module14_Phylogeny)')
    
    parser.add_argument('-og', '--outgroup',
                       dest='outgroups',
                       nargs='+',
                       default=None,
                       help='Outgroup GenBank file(s). Can specify filename or full path. Multiple outgroups will be ordered as specified.')
    
    parser.add_argument('--macse',
                       action='store_true',
                       default=False,
                       help='Use MACSE for protein-coding genes ONLY (skips introns and intergenic). Default: MAFFT for all.')
    
    parser.add_argument('--genes-only',
                       action='store_true',
                       default=False,
                       help='Process only protein-coding genes, skip introns and intergenic regions.')
    
    parser.add_argument('--introns-only',
                       action='store_true',
                       default=False,
                       help='Process only introns, skip protein-coding genes and intergenic regions.')
    
    parser.add_argument('--intergenic-only',
                       action='store_true',
                       default=False,
                       help='Process only intergenic spacers, skip protein-coding genes and introns.')
    
    parser.add_argument('--complete-only',
                       action='store_true',
                       default=False,
                       help='Generate only the complete concatenated matrix (all regions), skip individual alignments.')
    
    parser.add_argument('--iqtree',
                       action='store_true',
                       default=False,
                       help='Run IQ-TREE phylogenetic inference on generated matrices (recommended for phylogeny)')
    
    parser.add_argument('--threads',
                       type=str,
                       default='AUTO',
                       help='Number of threads for IQ-TREE (default: AUTO for automatic detection)')
    
    # Parse arguments
    args = parser.parse_args()
    
    input_dir = args.input_dir
    output_dir = args.output_dir
    outgroup_files = args.outgroups
    use_macse_flag = args.macse
    genes_only = args.genes_only
    introns_only = args.introns_only
    intergenic_only = args.intergenic_only
    complete_only = args.complete_only
    run_iqtree = args.iqtree
    threads = args.threads
    
    print("=" * 80)
    print("CGAS Module 14: Phylogenetic Matrix Builder")
    print("=" * 80)
    
    # Show input directory
    if input_dir == '.':
        print(f"\nInput directory: Current directory")
    else:
        print(f"\nInput directory: {input_dir}")
    
    # Check if input directory exists
    if not os.path.isdir(input_dir):
        print(f"\nError: Input directory '{input_dir}' not found!")
        sys.exit(1)
    
    # Find GenBank files
    genbank_files = list(Path(input_dir).glob("*.gb")) + list(Path(input_dir).glob("*.gbk"))
    
    if not genbank_files:
        print(f"\nError: No GenBank files (.gb or .gbk) found in '{input_dir}'")
        if input_dir == '.':
            print("Tip: Specify a directory with GenBank files:")
            print("  python module10_phylogenetic_matrix_builder.py genbank_directory/")
        sys.exit(1)
    
    print(f"\nFound {len(genbank_files)} GenBank files")
    
    # Parse outgroup files
    outgroup_names = []
    if outgroup_files:
        print(f"\nProcessing outgroup specification:")
        outgroup_names = parse_outgroup_files(outgroup_files, genbank_files)
        
        if outgroup_names:
            print(f"\n✓ {len(outgroup_names)} outgroup(s) will be placed at the top of alignments")
        else:
            print("\n⚠ Warning: No valid outgroup sequences found")
            print("  Continuing without outgroup ordering")
    
    # Determine which regions to process
    if use_macse_flag:
        print("\n--macse flag specified: Will process ONLY protein-coding genes with MACSE")
        process_genes = True
        process_introns = False
        process_intergenic = False
        cds_only = True
    elif genes_only:
        print("\n--genes-only flag specified: Processing ONLY protein-coding genes")
        process_genes = True
        process_introns = False
        process_intergenic = False
        cds_only = True
    elif introns_only:
        print("\n--introns-only flag specified: Processing ONLY introns")
        process_genes = False
        process_introns = True
        process_intergenic = False
        cds_only = False
    elif intergenic_only:
        print("\n--intergenic-only flag specified: Processing ONLY intergenic spacers")
        process_genes = False
        process_introns = False
        process_intergenic = True
        cds_only = False
    elif complete_only:
        print("\n--complete-only flag specified: Processing all regions but generating ONLY the complete concatenated matrix")
        process_genes = True
        process_introns = True
        process_intergenic = True
        cds_only = False
    else:
        print("\nNo region flags specified: Processing ALL regions with MAFFT")
        process_genes = True
        process_introns = True
        process_intergenic = True
        cds_only = False
    
    # Check for aligners
    if not check_mafft():
        print("\nError: MAFFT not found!")
        print("Please install MAFFT: https://mafft.cbrc.jp/alignment/software/")
        sys.exit(1)
    
    # Determine if MACSE should be used
    use_macse = False
    macse_path = None
    
    if use_macse_flag:
        if check_macse():
            macse_path = find_macse_path()
            if macse_path:
                use_macse = True
                print(f"✓ MACSE found: {macse_path}")
                print(f"  Will use MACSE for protein-coding genes (codon-aware alignment)")
            else:
                print(f"\n❌ ERROR: --macse specified but MACSE path not found!")
                sys.exit(1)
        else:
            print(f"\n❌ ERROR: --macse specified but MACSE not available!")
            print(f"  Install MACSE: https://bioweb.supagro.inra.fr/macse/")
            sys.exit(1)
    else:
        print("  Using MAFFT for alignment")
    
    # Set output directory
    print(f"\nOutput directory: {output_dir}")
    
    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    
    # Always create subdirectories inside output_dir
    alignments_dir = os.path.join(output_dir, "01_individual_alignments")
    matrices_dir = os.path.join(output_dir, "02_concatenated_matrices")
    
    genes_dir = os.path.join(alignments_dir, "genes")
    introns_dir = os.path.join(alignments_dir, "introns")
    intergenic_dir = os.path.join(alignments_dir, "intergenic")
    
    os.makedirs(genes_dir, exist_ok=True)
    os.makedirs(introns_dir, exist_ok=True)
    os.makedirs(intergenic_dir, exist_ok=True)
    os.makedirs(matrices_dir, exist_ok=True)
    
    # Extract features
    print("\n" + "=" * 80)
    print("STEP 1: Extracting features from GenBank files")
    print("=" * 80)
    
    protein_coding_genes, introns, intergenic, synonym_replacements = extract_features_from_genbank(genbank_files)
    
    # Log synonym replacements if any occurred
    if synonym_replacements:
        print(f"\n⚠ Gene name normalization performed ({len(synonym_replacements)} replacements):")
        print(f"   Check SUMMARY_REPORT.txt for details")
    
    print(f"\nExtracted:")
    if process_genes:
        print(f"  - {len(protein_coding_genes)} protein-coding genes")
    if process_introns:
        print(f"  - {len(introns)} introns")
    if process_intergenic:
        print(f"  - {len(intergenic)} intergenic spacers")
    
    # Align sequences
    print("\n" + "=" * 80)
    print("STEP 2: Aligning sequences")
    print("=" * 80)
    
    if outgroup_names:
        print(f"(Outgroups will be placed at top of alignments)")
    
    # Store alignments for concatenation
    gene_alignments = {}
    intron_alignments = {}
    intergenic_alignments = {}
    
    # Track region order (for concatenation)
    gene_order = []
    intron_order = []
    intergenic_order = []
    
    # Always save individual alignments for user access
    # Align protein-coding genes
    if process_genes:
        if use_macse:
            print(f"\nAligning {len(protein_coding_genes)} protein-coding genes with MACSE (codon-aware)...")
        else:
            print(f"\nAligning {len(protein_coding_genes)} protein-coding genes with MAFFT...")
        
        for i, (gene_name, seq_records) in enumerate(sorted(protein_coding_genes.items()), 1):
            print(f"  [{i}/{len(protein_coding_genes)}] {gene_name}: {len(seq_records)} sequences", end=" ")
            
            output_fasta = os.path.join(genes_dir, f"{gene_name}.fasta")
            output_alignment = os.path.join(genes_dir, f"{gene_name}_aligned.fasta")
            
            if use_macse:
                alignment = align_sequences_macse(seq_records, output_fasta, output_alignment, macse_path, outgroup_names)
            else:
                alignment = align_sequences(seq_records, output_fasta, output_alignment, outgroup_names)
            
            if alignment:
                gene_alignments[gene_name] = alignment
                gene_order.append(gene_name)
                print("✓")
            else:
                print("✗ (alignment failed)")
        
        print(f"\n✓ Completed aligning {len(gene_alignments)} protein-coding genes")
    
    # Align introns
    if process_introns:
        print(f"\nAligning {len(introns)} introns with MAFFT...")
        
        for i, (intron_name, seq_records) in enumerate(sorted(introns.items()), 1):
            print(f"  [{i}/{len(introns)}] {intron_name}: {len(seq_records)} sequences", end=" ")
            
            output_fasta = os.path.join(introns_dir, f"{intron_name}.fasta")
            output_alignment = os.path.join(introns_dir, f"{intron_name}_aligned.fasta")
            
            alignment = align_sequences(seq_records, output_fasta, output_alignment, outgroup_names)
            
            if alignment:
                intron_alignments[intron_name] = alignment
                intron_order.append(intron_name)
                print("✓")
            else:
                print("✗ (alignment failed)")
        
        print(f"\n✓ Completed aligning {len(intron_alignments)} introns")
    
    # Align intergenic spacers
    if process_intergenic:
        print(f"\nAligning {len(intergenic)} intergenic spacers with MAFFT...")
        
        for i, (spacer_name, seq_records) in enumerate(sorted(intergenic.items()), 1):
            print(f"  [{i}/{len(intergenic)}] {spacer_name}: {len(seq_records)} sequences", end=" ")
            
            output_fasta = os.path.join(intergenic_dir, f"{spacer_name}.fasta")
            output_alignment = os.path.join(intergenic_dir, f"{spacer_name}_aligned.fasta")
            
            alignment = align_sequences(seq_records, output_fasta, output_alignment, outgroup_names)
            
            if alignment:
                intergenic_alignments[spacer_name] = alignment
                intergenic_order.append(spacer_name)
                print("✓")
            else:
                print("✗ (alignment failed)")
        
        print(f"\n✓ Completed aligning {len(intergenic_alignments)} intergenic spacers")
    
    # Generate concatenated matrices
    print("\n" + "=" * 80)
    print("STEP 3: Generating concatenated matrices")
    print("=" * 80)
    
    # Create partition files directory
    partitions_dir = os.path.join(matrices_dir, "partitions")
    os.makedirs(partitions_dir, exist_ok=True)
    
    # Generate individual matrices (skip if --complete-only)
    if not complete_only:
        print("\n" + "=" * 80)
        print("STEP 3b: Generating individual matrices")
        print("=" * 80)
        
        # Genes-only matrix
        if process_genes and gene_alignments:
            print("\nGenerating genes-only matrix...")
            genes_matrix = os.path.join(matrices_dir, "genes_only.fasta")
            length, regions, partition_info = concatenate_alignments(
                gene_alignments, gene_order, genes_matrix, "genes_only", outgroup_names
            )
            
            if length > 0:
                print(f"  ✓ Created genes-only matrix: {length} bp, {regions} regions")
                
                # Write partition files
                for fmt in ["raxml", "nexus", "iqtree", "phylip"]:
                    partition_file = os.path.join(partitions_dir, f"genes_only.{fmt}")
                    write_partition_file(partition_info, partition_file, fmt)
                
                # Run IQ-TREE if requested
                if run_iqtree:
                    partition_file = os.path.join(partitions_dir, "genes_only.iqtree")
                    success, tree_file, log_file = run_iqtree_analysis(
                        genes_matrix, matrices_dir, partition_file, "genes_only", threads
                    )
                    if success:
                        summarize_iqtree_results(log_file, tree_file, "genes-only matrix")
        
        # Introns-only matrix
        if process_introns and intron_alignments:
            print("\nGenerating introns-only matrix...")
            introns_matrix = os.path.join(matrices_dir, "introns_only.fasta")
            length, regions, partition_info = concatenate_alignments(
                intron_alignments, intron_order, introns_matrix, "introns_only", outgroup_names
            )
            
            if length > 0:
                print(f"  ✓ Created introns-only matrix: {length} bp, {regions} regions")
                
                # Write partition files
                for fmt in ["raxml", "nexus", "iqtree", "phylip"]:
                    partition_file = os.path.join(partitions_dir, f"introns_only.{fmt}")
                    write_partition_file(partition_info, partition_file, fmt)
                
                # Run IQ-TREE if requested
                if run_iqtree:
                    partition_file = os.path.join(partitions_dir, "introns_only.iqtree")
                    success, tree_file, log_file = run_iqtree_analysis(
                        introns_matrix, matrices_dir, partition_file, "introns_only", threads
                    )
                    if success:
                        summarize_iqtree_results(log_file, tree_file, "introns-only matrix")
        
        # Intergenic-only matrix
        if process_intergenic and intergenic_alignments:
            print("\nGenerating intergenic-only matrix...")
            intergenic_matrix = os.path.join(matrices_dir, "intergenic_only.fasta")
            length, regions, partition_info = concatenate_alignments(
                intergenic_alignments, intergenic_order, intergenic_matrix, "intergenic_only", outgroup_names
            )
            
            if length > 0:
                print(f"  ✓ Created intergenic-only matrix: {length} bp, {regions} regions")
                
                # Write partition files
                for fmt in ["raxml", "nexus", "iqtree", "phylip"]:
                    partition_file = os.path.join(partitions_dir, f"intergenic_only.{fmt}")
                    write_partition_file(partition_info, partition_file, fmt)
                
                # Run IQ-TREE if requested
                if run_iqtree:
                    partition_file = os.path.join(partitions_dir, "intergenic_only.iqtree")
                    success, tree_file, log_file = run_iqtree_analysis(
                        intergenic_matrix, matrices_dir, partition_file, "intergenic_only", threads
                    )
                    if success:
                        summarize_iqtree_results(log_file, tree_file, "intergenic-only matrix")
    
    # Complete matrix (generated when --complete-only is specified OR when no specific region flag is given)
    # This ensures that --genes-only, --introns-only, or --intergenic-only produce ONLY their respective matrices
    should_generate_complete = complete_only or not (genes_only or introns_only or intergenic_only or use_macse_flag)
    
    if should_generate_complete:
        print("\nGenerating complete concatenated matrix (all regions)...")
        all_alignments = {}
        all_order = []
        
        if gene_alignments:
            all_alignments.update(gene_alignments)
            all_order.extend(gene_order)
        
        if intron_alignments:
            all_alignments.update(intron_alignments)
            all_order.extend(intron_order)
        
        if intergenic_alignments:
            all_alignments.update(intergenic_alignments)
            all_order.extend(intergenic_order)
        
        if all_alignments:
            complete_matrix = os.path.join(matrices_dir, "complete_matrix.fasta")
            length, regions, partition_info = concatenate_alignments(
                all_alignments, all_order, complete_matrix, "complete", outgroup_names
            )
            
            if length > 0:
                print(f"  ✓ Created complete matrix: {length} bp, {regions} regions")
                
                # Write partition files
                for fmt in ["raxml", "nexus", "iqtree", "phylip"]:
                    partition_file = os.path.join(partitions_dir, f"complete_matrix.{fmt}")
                    write_partition_file(partition_info, partition_file, fmt)
                
                # Run IQ-TREE if requested
                if run_iqtree:
                    partition_file = os.path.join(partitions_dir, "complete_matrix.iqtree")
                    success, tree_file, log_file = run_iqtree_analysis(
                        complete_matrix, matrices_dir, partition_file, "complete_matrix", threads
                    )
                    if success:
                        summarize_iqtree_results(log_file, tree_file, "complete matrix")
        else:
            print("  ✗ No alignments found to concatenate")
    
    # Generate summary report
    print("\n" + "=" * 80)
    print("STEP 4: Generating summary report")
    print("=" * 80)
    
    summary_file = os.path.join(output_dir, "SUMMARY_REPORT.txt")
    
    with open(summary_file, "w") as f:
        f.write(f"CGAS Module 14: Phylogenetic Matrix Builder\n")
        f.write(f"{'='*80}\n\n")
        f.write(f"Analysis Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input Directory: {input_dir}\n")
        f.write(f"Output Directory: {output_dir}\n\n")
        
        f.write(f"Parameters:\n")
        f.write(f"  - Process Genes: {process_genes}\n")
        f.write(f"  - Process Introns: {process_introns}\n")
        f.write(f"  - Process Intergenic: {process_intergenic}\n")
        f.write(f"  - Use MACSE: {use_macse}\n")
        f.write(f"  - Run IQ-TREE: {run_iqtree}\n")
        f.write(f"  - Complete Only: {complete_only}\n")
        f.write(f"  - Threads: {threads}\n\n")
        
        f.write(f"Input Files ({len(genbank_files)}):\n")
        for gb_file in genbank_files:
            f.write(f"  - {os.path.basename(gb_file)}\n")
        
        if outgroup_names:
            f.write(f"\nOutgroups ({len(outgroup_names)}):\n")
            for og_name in outgroup_names:
                f.write(f"  - {og_name}\n")
        
        f.write(f"\nExtracted Regions:\n")
        if process_genes:
            f.write(f"  - Protein-coding genes: {len(protein_coding_genes)}\n")
        if process_introns:
            f.write(f"  - Introns: {len(introns)}\n")
        if process_intergenic:
            f.write(f"  - Intergenic spacers: {len(intergenic)}\n")
        
        f.write(f"\nAlignments Generated:\n")
        if process_genes:
            f.write(f"  - Protein-coding genes: {len(gene_alignments)}\n")
        if process_introns:
            f.write(f"  - Introns: {len(intron_alignments)}\n")
        if process_intergenic:
            f.write(f"  - Intergenic spacers: {len(intergenic_alignments)}\n")
        
        f.write(f"\nOutput Files:\n")
        f.write(f"  - Individual alignments: {alignments_dir}/\n")
        f.write(f"  - Concatenated matrices: {matrices_dir}/\n")
        f.write(f"  - Partition files: {partitions_dir}/\n")
        f.write(f"  - Summary report: {summary_file}\n")
        
        if synonym_replacements:
            f.write(f"\nGene Name Normalizations:\n")
            for original, normalized in sorted(synonym_replacements.items()):
                f.write(f"  - {original} → {normalized}\n")
    
    print(f"  ✓ Summary report written to: {summary_file}")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    
    if complete_only:
        print("\n✓ Complete matrix generated successfully")
        print("  Individual alignments are available in:", alignments_dir)
    else:
        print("\n✓ All matrices generated successfully")
    
    if run_iqtree:
        print("\n✓ Phylogenetic trees generated successfully")
        print("  Check output directories for tree files and reports")
    
    print("\nOutput directory:", output_dir)
    print("=" * 80)


if __name__ == "__main__":
    # Import datetime for summary report
    import datetime
    main()