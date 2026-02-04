#!/usr/bin/env python3
"""
CGAS Module 11: Gene and tRNA Intron Extraction and Analysis
=============================================================

This script extracts and analyzes intron positions and lengths from both gene 
(CDS/gene features) and tRNA features in GenBank format files. Results are saved 
in Excel format with separate sheets for gene introns and tRNA introns, featuring
publication-quality formatting.

Key Features:
1. Extracts intron data from CDS/gene features
2. Extracts intron data from tRNA features
3. Validates intron lengths (1-15000 bp)
4. Uses Roman numerals for multiple introns
5. Groups output by species
6. Generates publication-ready Excel files with:
   - Proper formatting (headers, italics for genes, species grouping)
   - Separate sheets for genes and tRNAs
   - Abbreviation footnotes
   - Auto-adjusted column widths

Author: Abdullah
Version: 2.0 (Module 11 - CGAS Integration)
Date: January 2026

Dependencies:
    Python: biopython, pandas, openpyxl

Usage:
    python cgas_module11.py
    python cgas_module11.py -i genbank_files/
    python cgas_module11.py -i data/ -o results/

Output:
    Module11_Intron_Analysis/
        - intron_data_YYYYMMDD_HHMMSS.xlsx
          * Sheet 1: Gene Introns (CDS/gene features)
          * Sheet 2: tRNA Introns (tRNA features)
          * Sheet 3: Gene Intron Size Comparison (if multiple species)
          * Sheet 4: tRNA Intron Size Comparison (if multiple species)

Notes:
    - Processes all GenBank files (.gb, .gbff, .genbank, .gbk) in input directory
    - Intron length must be 1-15000 bp
    - Handles multiple introns per gene with Roman numeral notation
    - Groups results by species for clarity
    - IR (Inverted Repeat) regions: Genes appearing twice (e.g., trnA-UGC, trnI-GAU)
      with identical intron sizes are automatically deduplicated and counted only
      once per species in comparative analysis
"""

import os
import sys
import argparse
import logging
from datetime import datetime
from pathlib import Path
from typing import List, Tuple, Dict

try:
    from Bio import SeqIO
    from Bio.SeqFeature import CompoundLocation
except ImportError as e:
    print(f"Error: BioPython not installed: {e}")
    print("Please install required packages using:")
    print("pip install biopython")
    sys.exit(1)

try:
    import pandas as pd
    import numpy as np
except ImportError as e:
    print(f"Error: pandas/numpy not installed: {e}")
    print("Please install required packages using:")
    print("pip install pandas numpy")
    sys.exit(1)

try:
    import openpyxl
    from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
except ImportError as e:
    print(f"Error: openpyxl not installed: {e}")
    print("Please install required packages using:")
    print("pip install openpyxl")
    sys.exit(1)


# ============================================================================
# CONSTANTS
# ============================================================================

OUTPUT_FOLDER = "Module11_Intron_Analysis"
MAX_INTRON_LENGTH = 15000
VALID_EXTENSIONS = ('.gb', '.gbff', '.genbank', '.gbk')


# ============================================================================
# LOGGING CONFIGURATION
# ============================================================================

def setup_logging():
    """Configure logging with proper encoding for all platforms."""
    # Handle Windows console encoding
    if sys.platform == 'win32':
        try:
            if hasattr(sys.stdout, 'reconfigure'):
                sys.stdout.reconfigure(encoding='utf-8')
        except (AttributeError, OSError):
            pass
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


logger = setup_logging()


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def to_roman(num: int) -> str:
    """
    Convert integer to Roman numeral for intron numbering.
    
    Parameters:
    -----------
    num : int
        Integer to convert (1-3999)
        
    Returns:
    --------
    str
        Roman numeral representation
        
    Examples:
    ---------
    >>> to_roman(1)
    'I'
    >>> to_roman(4)
    'IV'
    >>> to_roman(10)
    'X'
    """
    val = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1]
    syms = ['M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I']
    roman_num = ''
    
    for i in range(len(val)):
        count = int(num / val[i])
        if count:
            roman_num += syms[i] * count
            num -= val[i] * count
    
    return roman_num


def find_genbank_files(directory: str) -> List[str]:
    """
    Find all GenBank files in the specified directory.
    
    Parameters:
    -----------
    directory : str
        Path to directory to search
        
    Returns:
    --------
    List[str]
        List of full paths to GenBank files
    """
    genbank_files = []
    
    for file in os.listdir(directory):
        if file.endswith(VALID_EXTENSIONS):
            genbank_files.append(os.path.join(directory, file))
    
    return sorted(genbank_files)


# ============================================================================
# INTRON EXTRACTION FUNCTIONS
# ============================================================================

def extract_gene_introns_from_genbank(gb_file: str) -> List[List]:
    """
    Extract gene/CDS intron information from a GenBank file.
    
    This function identifies compound locations in CDS or gene features,
    calculates intron positions and lengths, and returns formatted data
    suitable for Excel output.
    
    Parameters:
    -----------
    gb_file : str
        Path to GenBank format file
        
    Returns:
    --------
    List[List]
        List of gene records where each record contains:
        [filename, species, gene_name, intron_label_1, start_1, end_1, 
         length_1, intron_label_2, start_2, end_2, length_2, ...]
        
    Notes:
    ------
    - Only processes introns with length 1-15000 bp
    - Uses Roman numerals for multiple introns (I, II, III, etc.)
    - Single introns are labeled simply as "Intron"
    """
    gene_introns = []
    filename = os.path.splitext(os.path.basename(gb_file))[0]
    
    try:
        for record in SeqIO.parse(gb_file, "genbank"):
            species = record.annotations.get('organism', 'Unknown species')
            
            for feature in record.features:
                if feature.type in ['CDS', 'gene']:
                    gene = feature.qualifiers.get('gene', ['Unknown'])[0]
                    
                    # Check if feature has compound location (multiple exons)
                    if isinstance(feature.location, CompoundLocation):
                        exons = sorted(feature.location.parts, key=lambda x: x.start)
                        
                        # Calculate valid introns
                        valid_introns = []
                        for i in range(len(exons) - 1):
                            intron_start = int(exons[i].end) + 1
                            intron_end = int(exons[i + 1].start)
                            intron_length = intron_end - intron_start + 1
                            
                            # Validate intron length
                            if 0 < intron_length <= MAX_INTRON_LENGTH:
                                valid_introns.append((intron_start, intron_end, intron_length))
                        
                        # Format intron data with labels
                        if valid_introns:
                            intron_info = []
                            num_introns = len(valid_introns)
                            
                            for idx, (intron_start, intron_end, intron_length) in enumerate(valid_introns, 1):
                                # Use Roman numerals for multiple introns
                                if num_introns == 1:
                                    intron_label = "Intron"
                                else:
                                    intron_label = f"Intron {to_roman(idx)}"
                                
                                intron_info.extend([
                                    intron_label,
                                    intron_start,
                                    intron_end,
                                    intron_length
                                ])
                            
                            # Add complete record
                            gene_data = [filename, species, gene] + intron_info
                            gene_introns.append(gene_data)
    
    except Exception as e:
        logger.warning(f"  Warning: Error processing {os.path.basename(gb_file)}: {e}")
    
    return gene_introns


def extract_tRNA_introns_from_genbank(gb_file: str) -> List[List]:
    """
    Extract tRNA intron information from a GenBank file.
    
    This function identifies compound locations in tRNA features,
    calculates intron positions and lengths, and returns formatted data
    suitable for Excel output.
    
    Parameters:
    -----------
    gb_file : str
        Path to GenBank format file
        
    Returns:
    --------
    List[List]
        List of tRNA records where each record contains:
        [filename, species, trna_gene, intron_label_1, start_1, end_1, 
         length_1, intron_label_2, start_2, end_2, length_2, ...]
        
    Notes:
    ------
    - Only processes introns with length 1-15000 bp
    - Uses Roman numerals for multiple introns (I, II, III, etc.)
    - Single introns are labeled simply as "Intron"
    - Attempts to get gene name from 'gene' qualifier, falls back to 'product'
    """
    trna_introns = []
    filename = os.path.splitext(os.path.basename(gb_file))[0]
    
    try:
        for record in SeqIO.parse(gb_file, "genbank"):
            species = record.annotations.get('organism', 'Unknown species')
            
            for feature in record.features:
                if feature.type == 'tRNA':
                    # Get gene name (try 'gene' first, then 'product')
                    gene = feature.qualifiers.get('gene', 
                           feature.qualifiers.get('product', ['Unknown']))[0]
                    
                    # Check if feature has compound location (multiple exons)
                    if isinstance(feature.location, CompoundLocation):
                        exons = sorted(feature.location.parts, key=lambda x: x.start)
                        
                        # Calculate valid introns
                        valid_introns = []
                        for i in range(len(exons) - 1):
                            intron_start = exons[i].end + 1
                            intron_end = exons[i + 1].start - 1
                            intron_length = int(intron_end - intron_start + 1)
                            
                            # Validate intron length
                            if intron_length > 0 and intron_length <= MAX_INTRON_LENGTH:
                                valid_introns.append((int(intron_start), int(intron_end), intron_length))
                        
                        # Format intron data with labels
                        if valid_introns:
                            intron_info = []
                            num_introns = len(valid_introns)
                            
                            for idx, (intron_start, intron_end, intron_length) in enumerate(valid_introns, 1):
                                # Use Roman numerals for multiple introns
                                if num_introns == 1:
                                    intron_label = "Intron"
                                else:
                                    intron_label = f"Intron {to_roman(idx)}"
                                
                                intron_info.extend([
                                    intron_label,
                                    intron_start,
                                    intron_end,
                                    intron_length
                                ])
                            
                            # Add complete record
                            gene_data = [filename, species, gene] + intron_info
                            trna_introns.append(gene_data)
    
    except Exception as e:
        logger.warning(f"  Warning: Error processing {os.path.basename(gb_file)}: {e}")
    
    return trna_introns


# ============================================================================
# COMPARATIVE ANALYSIS FUNCTIONS
# ============================================================================

def analyze_intron_size_ranges(data: List[List], is_trna: bool = False) -> Dict:
    """
    Analyze intron size ranges across multiple species for each gene.
    
    This function groups intron data by gene name and intron position (e.g., Intron I, II),
    then calculates statistics across all species for each gene-intron combination.
    
    IMPORTANT: Handles duplicate genes from Inverted Repeat (IR) regions by only 
    counting unique species-gene-size combinations once. Genes in IRs appear twice 
    with identical intron sizes but should only be counted once per species.
    
    Parameters:
    -----------
    data : List[List]
        List of intron records [filename, species, gene, intron_label, start, end, length, ...]
    is_trna : bool
        Whether this is tRNA data (affects naming in output)
        
    Returns:
    --------
    Dict
        Dictionary mapping (gene, intron_label) to statistics:
        {
            ('petB', 'Intron'): {
                'min_size': 507,
                'max_size': 526,
                'difference': 19,
                'mean_size': 516.5,
                'species_count': 2,
                'min_species': 'Species A',
                'max_species': 'Species B',
                'all_sizes': [507, 526]
            }
        }
    """
    # Dictionary to store intron data: {(gene, intron_label): set((species, length))}
    # Using set to automatically handle duplicates from IR regions
    intron_data = {}
    
    for record in data:
        filename, species, gene = record[0], record[1], record[2]
        intron_info = record[3:]
        
        # Process each intron in the record (groups of 4: label, start, end, length)
        for i in range(0, len(intron_info), 4):
            if i + 3 < len(intron_info):
                intron_label = intron_info[i]
                intron_length = intron_info[i + 3]
                
                key = (gene, intron_label)
                if key not in intron_data:
                    intron_data[key] = set()
                
                # Add as tuple to set - duplicates (same species + length) will be automatically removed
                # This handles IR regions where the same gene appears twice with identical intron sizes
                intron_data[key].add((species, intron_length))
    
    # Calculate statistics for each gene-intron combination
    results = {}
    for (gene, intron_label), size_set in intron_data.items():
        if len(size_set) > 0:  # Only process if we have data
            # Convert set back to list for processing
            size_list = list(size_set)
            lengths = [length for _, length in size_list]
            species_names = [species for species, _ in size_list]
            
            min_size = min(lengths)
            max_size = max(lengths)
            mean_size = sum(lengths) / len(lengths)
            
            # Find species with min and max sizes
            min_idx = lengths.index(min_size)
            max_idx = lengths.index(max_size)
            min_species = species_names[min_idx]
            max_species = species_names[max_idx]
            
            # Count unique species (not total entries)
            unique_species = set(species_names)
            
            results[(gene, intron_label)] = {
                'min_size': min_size,
                'max_size': max_size,
                'difference': max_size - min_size,
                'mean_size': round(mean_size, 2),
                'species_count': len(unique_species),  # Now counts unique species only
                'min_species': min_species,
                'max_species': max_species,
                'all_sizes': lengths
            }
    
    return results


def create_comparison_sheet(workbook, sheet_name: str, comparison_data: Dict, is_trna: bool = False):
    """
    Create a comparison sheet showing intron size ranges across species.
    
    Parameters:
    -----------
    workbook : openpyxl.Workbook
        Excel workbook object
    sheet_name : str
        Name for the new sheet
    comparison_data : Dict
        Dictionary from analyze_intron_size_ranges()
    is_trna : bool
        Whether this is tRNA data
    """
    ws = workbook.create_sheet(sheet_name)
    
    # Define styles
    header_font = Font(bold=True, size=11, color='FFFFFF')
    header_fill = PatternFill(start_color='4472C4', end_color='4472C4', fill_type='solid')
    header_alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
    
    gene_font = Font(italic=True, size=10)  # Gene names: italic
    species_font = Font(italic=True, size=10)  # Species names: italic
    regular_font = Font(size=10)
    center_alignment = Alignment(horizontal='center', vertical='center')
    left_alignment = Alignment(horizontal='left', vertical='center')
    
    # Thin border
    thin_border = Border(
        left=Side(style='thin'),
        right=Side(style='thin'),
        top=Side(style='thin'),
        bottom=Side(style='thin')
    )
    
    # Create headers
    headers = [
        'Gene' if not is_trna else 'tRNA Gene',
        'Intron',
        'Minimum Size (bp)',
        'Species (Min)',
        'Maximum Size (bp)',
        'Species (Max)',
        'Difference (bp)',
        'Mean Size (bp)',
        'No. of Species'
    ]
    
    # Write headers
    for col, header in enumerate(headers, 1):
        cell = ws.cell(1, col, header)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = header_alignment
        cell.border = thin_border
    
    # Sort data by gene name, then intron label
    sorted_keys = sorted(comparison_data.keys(), key=lambda x: (x[0], x[1]))
    
    # Write data rows
    row = 2
    for (gene, intron_label) in sorted_keys:
        stats = comparison_data[(gene, intron_label)]
        
        # Only show comparison if there are multiple species or variation
        if stats['species_count'] >= 1:
            # Gene name (italic)
            cell = ws.cell(row, 1, gene)
            cell.font = gene_font
            cell.alignment = left_alignment
            cell.border = thin_border
            
            # Intron label
            cell = ws.cell(row, 2, intron_label)
            cell.font = regular_font
            cell.alignment = center_alignment
            cell.border = thin_border
            
            # Minimum size
            cell = ws.cell(row, 3, stats['min_size'])
            cell.font = regular_font
            cell.alignment = center_alignment
            cell.border = thin_border
            
            # Species with minimum (italic)
            cell = ws.cell(row, 4, stats['min_species'])
            cell.font = species_font
            cell.alignment = center_alignment
            cell.border = thin_border
            
            # Maximum size
            cell = ws.cell(row, 5, stats['max_size'])
            cell.font = regular_font
            cell.alignment = center_alignment
            cell.border = thin_border
            
            # Species with maximum (italic)
            cell = ws.cell(row, 6, stats['max_species'])
            cell.font = species_font
            cell.alignment = center_alignment
            cell.border = thin_border
            
            # Difference
            cell = ws.cell(row, 7, stats['difference'])
            cell.font = regular_font
            cell.alignment = center_alignment
            cell.border = thin_border
            
            # Mean size
            cell = ws.cell(row, 8, stats['mean_size'])
            cell.font = regular_font
            cell.alignment = center_alignment
            cell.border = thin_border
            
            # Species count
            cell = ws.cell(row, 9, stats['species_count'])
            cell.font = regular_font
            cell.alignment = center_alignment
            cell.border = thin_border
            
            row += 1
    
    # Auto-adjust column widths
    for column in ws.columns:
        max_length = 0
        column_letter = column[0].column_letter
        for cell in column:
            if cell.value:
                max_length = max(max_length, len(str(cell.value)))
        # Set width with some padding
        adjusted_width = min(max_length + 3, 50)
        ws.column_dimensions[column_letter].width = adjusted_width
    
    # Add explanatory note
    note_row = row + 2
    note_cell = ws.cell(note_row, 1, "Note:")
    note_cell.font = Font(bold=True, size=10)
    
    note_row += 1
    note_text = (
        "This sheet compares intron sizes across multiple species for each gene. "
        "The difference column shows variation in intron length among species. "
        "Genes in Inverted Repeat (IR) regions are counted only once per species, "
        "even though they appear twice in the genome."
    )
    note_cell = ws.cell(note_row, 1, note_text)
    note_cell.font = Font(size=9)
    ws.merge_cells(start_row=note_row, start_column=1, end_row=note_row, end_column=9)
    
    logger.info(f"  ✓ Comparison sheet created: {sheet_name}")


# ============================================================================
# EXCEL FORMATTING AND OUTPUT
# ============================================================================

def add_abbreviation_footnotes(worksheet, start_row: int) -> int:
    """
    Add abbreviation footnotes to Excel worksheet for publication quality.
    
    Parameters:
    -----------
    worksheet : openpyxl.worksheet.worksheet.Worksheet
        Worksheet object to add footnotes to
    start_row : int
        Row number to start adding footnotes
        
    Returns:
    --------
    int
        Next available row number after footnotes
    """
    footnotes = [
        "Abbreviations:",
        "tRNA: Transfer RNA",
        "CDS: Coding Sequence",
        "bp: Base pairs"
    ]
    
    row = start_row
    for footnote in footnotes:
        cell = worksheet.cell(row, 1, footnote)
        if footnote == "Abbreviations:":
            cell.font = Font(bold=True, size=10)
        else:
            cell.font = Font(size=9)
        cell.alignment = Alignment(horizontal='left', vertical='top')
        row += 1
    
    return row


def create_combined_excel_output(gene_data: List[List], 
                                 trna_data: List[List], 
                                 output_excel: str) -> None:
    """
    Create publication-quality Excel file with gene and tRNA intron data.
    
    This function creates a multi-sheet Excel workbook with professional
    formatting including:
    - Separate sheets for gene and tRNA introns
    - Species grouping with merged cells and highlighting
    - Proper headers with bold formatting
    - Italic gene names
    - Auto-adjusted column widths
    - Abbreviation footnotes
    - Comparison sheets showing intron size ranges across species (if multiple species)
    
    Parameters:
    -----------
    gene_data : List[List]
        List of gene intron records
    trna_data : List[List]
        List of tRNA intron records
    output_excel : str
        Path to output Excel file
        
    Notes:
    ------
    - Records are grouped by species
    - Species names are displayed in merged cells with gray background
    - Gene names are italicized
    - Numerical data is center-aligned
    - Columns are auto-sized for readability
    - If multiple species are present, additional comparison sheets are created
      showing min/max/mean intron sizes for each gene across all species
    """
    wb = openpyxl.Workbook()
    
    # Remove default sheet
    if 'Sheet' in wb.sheetnames:
        wb.remove(wb['Sheet'])
    
    # Define styles
    header_font = Font(bold=True, size=11)
    header_alignment = Alignment(horizontal='center', vertical='center')
    species_font = Font(bold=True, italic=True, size=11)  # Species names: bold + italic
    species_alignment = Alignment(horizontal='center', vertical='center')
    species_fill = PatternFill(start_color='D9D9D9', end_color='D9D9D9', fill_type='solid')
    italic_font = Font(italic=True, size=10)  # Gene names: italic only
    regular_font = Font(size=10)
    left_alignment = Alignment(horizontal='left', vertical='center')
    center_alignment = Alignment(horizontal='center', vertical='center')
    
    # ========================================================================
    # Create Gene Introns sheet
    # ========================================================================
    if gene_data:
        ws = wb.create_sheet('Gene Introns')
        
        # Determine maximum number of introns across all genes
        max_introns = max((len(row) - 3) // 4 for row in gene_data)
        
        # Create headers (without File column for cleaner output)
        headers = ["Gene"]
        for i in range(1, max_introns + 1):
            if max_introns == 1:
                headers += ["Intron", "Start", "End", "Length (bp)"]
            else:
                roman = to_roman(i)
                headers += [f"Intron {roman}", f"Start", f"End", f"Length (bp)"]
        
        # Write headers
        for col, header in enumerate(headers, 1):
            cell = ws.cell(1, col, header)
            cell.font = header_font
            cell.alignment = header_alignment
        
        row = 2
        
        # Group data by species
        data_by_species = {}
        for record in gene_data:
            filename, species, gene = record[0], record[1], record[2]
            intron_data = record[3:]
            if species not in data_by_species:
                data_by_species[species] = []
            # Store without filename for cleaner output
            data_by_species[species].append([gene] + intron_data)
        
        # Write data grouped by species
        for species in sorted(data_by_species.keys()):
            species_records = data_by_species[species]
            
            # Add species row (merged across all columns)
            ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=len(headers))
            species_cell = ws.cell(row, 1, species)
            species_cell.font = species_font
            species_cell.alignment = species_alignment
            species_cell.fill = species_fill
            row += 1
            
            # Add data rows for this species
            for record in species_records:
                # Pad record to match column count
                while len(record) < len(headers):
                    record.append("")
                
                # Write gene name (italic)
                ws.cell(row, 1, record[0]).font = italic_font
                ws.cell(row, 1, record[0]).alignment = left_alignment
                
                # Write intron data (centered)
                for col_idx, value in enumerate(record[1:], 2):
                    cell = ws.cell(row, col_idx, value)
                    cell.font = regular_font
                    cell.alignment = center_alignment
                
                row += 1
        
        # Auto-adjust column widths
        for column in ws.columns:
            max_length = 0
            column_letter = column[0].column_letter
            for cell in column:
                if cell.value:
                    max_length = max(max_length, len(str(cell.value)))
            ws.column_dimensions[column_letter].width = min(max_length + 3, 50)
        
        # Add abbreviation footnotes
        add_abbreviation_footnotes(ws, row + 2)
    
    # ========================================================================
    # Create tRNA Introns sheet
    # ========================================================================
    if trna_data:
        ws = wb.create_sheet('tRNA Introns')
        
        # Determine maximum number of introns across all tRNAs
        max_introns = max((len(row) - 3) // 4 for row in trna_data)
        
        # Create headers (without File column for cleaner output)
        headers = ["tRNA Gene"]
        for i in range(1, max_introns + 1):
            if max_introns == 1:
                headers += ["Intron", "Start", "End", "Length (bp)"]
            else:
                roman = to_roman(i)
                headers += [f"Intron {roman}", f"Start", f"End", f"Length (bp)"]
        
        # Write headers
        for col, header in enumerate(headers, 1):
            cell = ws.cell(1, col, header)
            cell.font = header_font
            cell.alignment = header_alignment
        
        row = 2
        
        # Group data by species
        data_by_species = {}
        for record in trna_data:
            filename, species, gene = record[0], record[1], record[2]
            intron_data = record[3:]
            if species not in data_by_species:
                data_by_species[species] = []
            # Store without filename for cleaner output
            data_by_species[species].append([gene] + intron_data)
        
        # Write data grouped by species
        for species in sorted(data_by_species.keys()):
            species_records = data_by_species[species]
            
            # Add species row (merged across all columns)
            ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=len(headers))
            species_cell = ws.cell(row, 1, species)
            species_cell.font = species_font
            species_cell.alignment = species_alignment
            species_cell.fill = species_fill
            row += 1
            
            # Add data rows for this species
            for record in species_records:
                # Pad record to match column count
                while len(record) < len(headers):
                    record.append("")
                
                # Write gene name (italic)
                ws.cell(row, 1, record[0]).font = italic_font
                ws.cell(row, 1, record[0]).alignment = left_alignment
                
                # Write intron data (centered)
                for col_idx, value in enumerate(record[1:], 2):
                    cell = ws.cell(row, col_idx, value)
                    cell.font = regular_font
                    cell.alignment = center_alignment
                
                row += 1
        
        # Auto-adjust column widths
        for column in ws.columns:
            max_length = 0
            column_letter = column[0].column_letter
            for cell in column:
                if cell.value:
                    max_length = max(max_length, len(str(cell.value)))
            ws.column_dimensions[column_letter].width = min(max_length + 3, 50)
        
        # Add abbreviation footnotes
        add_abbreviation_footnotes(ws, row + 2)
    
    # ========================================================================
    # Create Comparison Sheets (if multiple species present)
    # ========================================================================
    
    # Check if we have multiple species in the data
    def get_unique_species(data: List[List]) -> set:
        """Extract unique species from data."""
        species_set = set()
        for record in data:
            if len(record) >= 2:
                species_set.add(record[1])  # Species is at index 1
        return species_set
    
    # Gene intron comparison
    if gene_data:
        unique_species = get_unique_species(gene_data)
        if len(unique_species) > 1:
            logger.info(f"  Creating gene intron comparison sheet ({len(unique_species)} species)...")
            gene_comparison = analyze_intron_size_ranges(gene_data, is_trna=False)
            if gene_comparison:
                create_comparison_sheet(wb, 'Gene Intron Size Comparison', gene_comparison, is_trna=False)
    
    # tRNA intron comparison
    if trna_data:
        unique_species = get_unique_species(trna_data)
        if len(unique_species) > 1:
            logger.info(f"  Creating tRNA intron comparison sheet ({len(unique_species)} species)...")
            trna_comparison = analyze_intron_size_ranges(trna_data, is_trna=True)
            if trna_comparison:
                create_comparison_sheet(wb, 'tRNA Intron Size Comparison', trna_comparison, is_trna=True)
    
    # Save workbook
    wb.save(output_excel)
    logger.info(f"  ✓ Excel file created: {os.path.basename(output_excel)}")


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def process_all_genbank_files(input_dir: str, output_dir: str) -> Tuple[int, int]:
    """
    Process all GenBank files in the specified directory.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing GenBank files
    output_dir : str
        Directory for output files
        
    Returns:
    --------
    Tuple[int, int]
        Tuple of (number of gene introns, number of tRNA introns)
        
    Raises:
    -------
    FileNotFoundError
        If no GenBank files found
    ValueError
        If no introns found in any file
    """
    all_gene_introns = []
    all_trna_introns = []
    
    # Find GenBank files
    gb_files = find_genbank_files(input_dir)
    
    if not gb_files:
        raise FileNotFoundError(
            f"No GenBank files found in {input_dir}\n"
            f"Expected extensions: {', '.join(VALID_EXTENSIONS)}"
        )
    
    logger.info(f"Found {len(gb_files)} GenBank file(s):")
    for gb_file in gb_files:
        logger.info(f"  - {os.path.basename(gb_file)}")
    
    logger.info(f"\n{'='*70}")
    logger.info("PROCESSING GENBANK FILES")
    logger.info(f"{'='*70}")
    
    # Process each file
    for idx, gb_file in enumerate(gb_files, 1):
        logger.info(f"[{idx}/{len(gb_files)}] Processing: {os.path.basename(gb_file)}")
        
        try:
            gene_rows = extract_gene_introns_from_genbank(gb_file)
            trna_rows = extract_tRNA_introns_from_genbank(gb_file)
            
            all_gene_introns.extend(gene_rows)
            all_trna_introns.extend(trna_rows)
            
            logger.info(f"  ✓ Genes with introns: {len(gene_rows)}")
            logger.info(f"  ✓ tRNAs with introns: {len(trna_rows)}")
        
        except Exception as e:
            logger.warning(f"  ⚠ Error processing {os.path.basename(gb_file)}: {e}")
            continue
    
    # Validate results
    if not all_gene_introns and not all_trna_introns:
        raise ValueError("No introns found in any GenBank file")
    
    # Generate output file
    logger.info(f"\n{'='*70}")
    logger.info("GENERATING OUTPUT")
    logger.info(f"{'='*70}")
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"intron_data_{timestamp}.xlsx")
    
    create_combined_excel_output(all_gene_introns, all_trna_introns, output_file)
    
    return len(all_gene_introns), len(all_trna_introns)


def main():
    """Main pipeline execution function with command-line arguments."""
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='CGAS Module 11: Gene and tRNA Intron Extraction',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s                          # Process current directory
  %(prog)s -i genbank_files/        # Process specific directory
  %(prog)s -i data/ -o results/     # Custom input and output

Output:
  Creates timestamped Excel file with two sheets:
    - Gene Introns: CDS/gene feature introns
    - tRNA Introns: tRNA feature introns

For more information, see the CGAS documentation.
        '''
    )
    
    parser.add_argument('-i', '--input',
                       type=str,
                       default='.',
                       help='Input directory containing GenBank files (default: current directory)')
    
    parser.add_argument('-o', '--output',
                       type=str,
                       default=None,
                       help='Output directory (default: Module11_Intron_Analysis in input directory)')
    
    args = parser.parse_args()
    
    # Set up directories
    input_dir = os.path.abspath(args.input)
    
    if args.output:
        output_dir = os.path.abspath(args.output)
    else:
        output_dir = os.path.join(input_dir, OUTPUT_FOLDER)
    
    # Print header
    print(f"\n{'='*70}")
    print("CGAS MODULE 11: INTRON EXTRACTION AND ANALYSIS")
    print(f"{'='*70}")
    print(f"\nInput directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Max intron length: {MAX_INTRON_LENGTH} bp")
    
    # Create output folder
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput folder: {os.path.basename(output_dir)}/")
    print()
    
    try:
        # Process GenBank files
        num_gene_introns, num_trna_introns = process_all_genbank_files(input_dir, output_dir)
        
        # Print summary
        print(f"\n{'='*70}")
        print("ANALYSIS COMPLETE")
        print(f"{'='*70}")
        print(f"✓ Gene introns found: {num_gene_introns}")
        print(f"✓ tRNA introns found: {num_trna_introns}")
        print(f"✓ Output saved in: {output_dir}/")
        print(f"{'='*70}\n")
    
    except FileNotFoundError as e:
        logger.error(f"\n❌ ERROR: {e}")
        logger.error("Please ensure GenBank files are in the input directory.")
        sys.exit(1)
    
    except ValueError as e:
        logger.error(f"\n❌ ERROR: {e}")
        sys.exit(1)
    
    except Exception as e:
        logger.error(f"\n❌ FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


# ============================================================================
# SCRIPT ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠ Analysis interrupted by user.")
        sys.exit(0)
    except Exception as e:
        print(f"\n\n❌ FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
