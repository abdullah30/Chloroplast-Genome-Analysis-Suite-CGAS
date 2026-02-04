# CGAS Module 12: Comprehensive Chloroplast SSR Analysis
## Complete Documentation and User Guide

---

## 📋 Table of Contents
1. [Introduction](#introduction)
2. [Quick Start - Just Run It!](#quick-start---just-run-it)
3. [Command-Line Usage with Examples](#command-line-usage-with-examples)
4. [Jupyter Notebook Usage](#jupyter-notebook-usage)
5. [Input Requirements](#input-requirements)
6. [Output Structure](#output-structure)
7. [Detailed Feature Explanation](#detailed-feature-explanation)
8. [SSR Detection Algorithm](#ssr-detection-algorithm)
9. [Genomic Location Classification](#genomic-location-classification)
10. [Gene Notation System](#gene-notation-system)
11. [R Visualization](#r-visualization)
12. [Troubleshooting](#troubleshooting)
13. [Examples](#examples)
14. [FAQ](#faq)
15. [Technical Specifications](#technical-specifications)
16. [References](#references)

---

## Introduction

**CGAS Module 12** is a comprehensive tool for detecting and analyzing Simple Sequence Repeats (SSRs) in chloroplast genomes from GenBank files. SSRs (also known as microsatellites) are short, tandemly repeated DNA sequences that are valuable molecular markers for evolutionary studies, population genetics, and comparative genomics.

This module performs complete SSR analysis with:
- **SSR detection** using customizable thresholds
- **Multi-level classification** by genomic region, motif type, and location
- **Precise annotation** within genes, introns, tRNAs, rRNAs, and intergenic spacers
- **Publication-ready outputs** with professional formatting and visualizations

### Scientific Applications:
- **Evolutionary Studies**: Compare SSR patterns across species
- **Population Genetics**: Identify polymorphic SSR markers
- **Genome Annotation**: Annotate repeat regions in chloroplast genomes
- **Comparative Genomics**: Study SSR conservation and evolution
- **Marker Development**: Identify potential SSR markers for PCR primers

### Part of CGAS (Chloroplast Genome Analysis Suite)
CGAS is a comprehensive suite of tools designed for chloroplast genome analysis, providing seamless integration across different analysis modules.

---

## Quick Start - Just Run It!

### ⚡ One Command Does It All!

```bash
# Option 1: Using cgas command
cgas --module 12

# Option 2: Using cgas-ssr shortcut
cgas-ssr

# Option 3: Using python directly
python cgas_module12.py
```

**What happens when you run this:**
1. ✅ Looks for GenBank files in your **current directory**
2. ✅ Creates `Module12_SSR_Analysis` folder automatically
3. ✅ Processes all `.gb`, `.gbff`, `.genbank`, `.gbk` files
4. ✅ Uses default thresholds: Mono=10, Di=5, Tri=4, Tetra=3, Penta=3, Hexa=3
5. ✅ Generates publication-ready Excel files with italic species names
6. ✅ Creates visualizations if R is installed

### 📁 Real-World Folder Examples

#### Example 1: Basic Usage (Easiest!)
```bash
# Your folder structure:
abdullah/
├── research_project/
│   ├── arabidopsis.gb
│   ├── oryza.gb
│   └── zea.gb

# Navigate to your data folder
cd /home/abdullah/research_project/

# Run the script (no arguments needed!)
python cgas_module12.py

# Output created automatically:
# Module12_SSR_Analysis/
# ├── all_ssrs_detailed.xlsx
# ├── ssr_summary.xlsx
# ├── ssr_by_species.xlsx
# └── Figures/
```

#### Example 2: Specify Input Folder
```bash
# Your folder structure:
/home/abdullah/data_analysis/
├── chloroplast_genomes/
│   ├── species1.gb
│   ├── species2.gb
│   └── species3.gb

# Run from ANYWHERE with -i flag
python cgas_module12.py -i /home/abdullah/data_analysis/chloroplast_genomes/

# Output created in:
# /home/abdullah/data_analysis/chloroplast_genomes/Module12_SSR_Analysis/
```

#### Example 3: Custom Input AND Output Folders
```bash
# You control both input and output locations
python cgas_module12.py -i abdullah/data_analysis -o SSRs_results/

# What this does:
# 1. Reads from: abdullah/data_analysis/*.gb
# 2. Saves to: SSRs_results/ (instead of Module12_SSR_Analysis)
```

#### Example 4: Custom SSR Thresholds
```bash
# Want more sensitive detection? Change thresholds!
python cgas_module12.py -t 8,4,3,3,3,3

# Threshold format: MONO,DI,TRI,TETRA,PENTA,HEXA
# Default: 10,5,4,3,3,3 ← More strict
# Custom: 8,4,3,3,3,3 ← More sensitive
```

**Threshold Comparison Table:**
| Motif Type | Default | Custom (More Sensitive) | What It Means |
|------------|---------|-------------------------|---------------|
| **Mono** | 10 repeats | 8 repeats | Detects shorter A/T repeats |
| **Di** | 5 repeats | 4 repeats | Finds shorter AT/TA repeats |
| **Tri** | 4 repeats | 3 repeats | Catches shorter trinucleotides |
| **Tetra** | 3 repeats | 3 repeats | Same for tetra |
| **Penta** | 3 repeats | 3 repeats | Same for penta |
| **Hexa** | 3 repeats | 3 repeats | Same for hexa |

### 🎯 Quick Start Cheat Sheet

```bash
# ====== cgas command ======
cgas --module 12                                         # Process current directory
cgas --module 12 -i your_data_folder/                    # Specify input folder
cgas --module 12 -i abdullah/data_analysis -o SSRs_results/  # Custom I/O
cgas --module 12 -t 8,4,3,3,3,3                          # Custom SSR thresholds
cgas --module 12 -i data/ -o out/ -t 12,6,5,4,4,4        # Combine everything

# ====== cgas-ssr shortcut ======
cgas-ssr                                                 # Process current directory
cgas-ssr -i your_data_folder/                            # Specify input folder
cgas-ssr -i chloroplast_data/ -o final_results/          # Custom I/O
cgas-ssr -t 8,4,3,3,3,3                                  # Custom SSR thresholds

# ====== python command ======
python cgas_module12.py                                  # Process current directory
python cgas_module12.py -i your_data_folder/             # Specify input folder
python cgas_module12.py -i abdullah/data_analysis -o SSRs_results/  # Custom I/O
python cgas_module12.py -t 8,4,3,3,3,3                   # Custom SSR thresholds
python cgas_module12.py --help                           # Get help
```

### 📊 What You Get (Output Files)
After running ANY of the above commands, you'll get:

```
Output_Folder/  # Either Module12_SSR_Analysis or your custom name
├── 📊 all_ssrs_detailed.xlsx      # Complete SSR catalog with all details
├── 📈 ssr_summary.xlsx            # Statistics table per species
├── 📁 ssr_by_species.xlsx         # Individual sheets for each species
├── 📄 ssr_summary_for_r.csv       # Data for R visualization
├── 📜 generate_ssr_figures.R      # R script (if R installed)
└── 🖼️  Figures/                   # Publication-quality figures (PDF+PNG)
    ├── SSR_Genomic_Regions.pdf
    ├── SSR_Motif_Types.pdf
    └── SSR_Genomic_Locations.pdf
```

---

---

## Command-Line Usage with Examples

### Complete Command Reference

```bash
# ====================================================================
# USING cgas COMMAND
# ====================================================================

# Process current directory
cgas --module 12

# Specify input directory
cgas --module 12 -i /home/user/chloroplast_genomes/

# Custom input and output directories
cgas --module 12 -i abdullah/data_analysis -o SSRs_results/

# Custom SSR thresholds
cgas --module 12 -t 8,4,3,3,3,3
```

```bash
# ====================================================================
# USING cgas-ssr SHORTCUT
# ====================================================================

# Process current directory
cgas-ssr

# With specific input and output
cgas-ssr -i abdullah/data_analysis -o SSRs_results/

# With custom thresholds
cgas-ssr -t 8,4,3,3,3,3
```

```bash
# ====================================================================
# USING python COMMAND
# ====================================================================

# ====================================================================
# BASIC USAGE
# ====================================================================

# 1. Run in current directory (simplest!)
python cgas_module12.py

# 2. Specify input folder
python cgas_module12.py -i genbank_files/

# 3. Specify BOTH input and output folders
python cgas_module12.py -i abdullah/data_analysis -o SSRs_results/

# 4. Custom SSR thresholds
python cgas_module12.py -t 8,4,3,3,3,3

# 5. Get help
python cgas_module12.py --help


# ====================================================================
# REAL-WORLD EXAMPLES WITH PATHS
# ====================================================================

# Example 1: Process data in "chloroplast_genomes" folder
python cgas_module12.py -i /home/user/chloroplast_genomes/

# Example 2: Save results to custom "analysis_results" folder
python cgas_module12.py -i data/ -o analysis_results/

# Example 3: Run from different directory
cd /home/abdullah/scripts/
python cgas_module12.py -i ../data/genbank/ -o ../results/ssrs/

# Example 4: Sensitive detection for detailed analysis
python cgas_module12.py -i monocots/ -o monocot_ssr_analysis/ -t 8,4,3,3,3,3

# Example 5: Strict detection for reliable markers
python cgas_module12.py -i dicots/ -o dicot_ssr_markers/ -t 12,6,5,4,4,4

# Example 6: Windows users
python cgas_module12.py -i "C:\Users\abdullah\data" -o "C:\Results\ssrs"
```

---

## Input Requirements

### Supported File Formats
- `.gb` (GenBank)
- `.gbff` (GenBank Flat File)
- `.genbank`
- `.gbk` (GenBank)

### GenBank File Requirements
Standard GenBank files from NCBI with:
- Sequence data
- Feature annotations (genes, CDS, tRNA, rRNA)
- Organism information

### File Organization
```
your_data_folder/
├── Arabidopsis_thaliana.gb      # Species name recommended in filename
├── Oryza_sativa.gb
├── Zea_mays.gb
└── outgroup_species.gb
```

**Note**: Species names are extracted from the `/organism=` qualifier in GenBank files.

---

## Output Structure

### Directory Structure
```
Module12_SSR_Analysis/           # Default output folder
├── all_ssrs_detailed.xlsx       # Complete SSR catalog
├── ssr_summary.xlsx              # Statistical summary per species
├── ssr_by_species.xlsx           # Individual species sheets
├── ssr_summary_for_r.csv         # Data for R visualization
├── generate_ssr_figures.R        # R script (editable)
└── Figures/                      # Publication-quality figures
    ├── SSR_Genomic_Regions.pdf   # SSR distribution by LSC/SSC/IR
    ├── SSR_Genomic_Regions.png   # (600 DPI)
    ├── SSR_Motif_Types.pdf       # SSR types (mono, di, tri, etc.)
    ├── SSR_Motif_Types.png       # (600 DPI)
    ├── SSR_Genomic_Locations.pdf # SSR locations (gene, intron, etc.)
    └── SSR_Genomic_Locations.png # (600 DPI)
```

### Excel File Details

#### 1. all_ssrs_detailed.xlsx
Complete catalog with columns:
- **Species** (*italicized*)
- **SSR_ID**
- **Motif**
- **Motif_Type**
- **Repeat_Count**
- **SSR_Length**
- **Start_Position**
- **End_Position**
- **Genomic_Region** (LSC, SSC, IRa, IRb)
- **Location_Type** (Gene, Intron, tRNA, rRNA, Intergenic)
- **Gene_Name** (*italicized* if applicable)

#### 2. ssr_summary.xlsx
Statistical overview:
- **Species** (*italicized*)
- **Total_SSRs**
- **Mono**, **Di**, **Tri**, **Tetra**, **Penta**, **Hexa** counts
- **LSC**, **SSC**, **IRa**, **IRb** counts
- **Gene**, **Intron**, **tRNA**, **rRNA**, **Intergenic** counts

#### 3. ssr_by_species.xlsx
Individual sheets for each species with detailed SSR lists.

---

## Detailed Feature Explanation

### 1. SSR Detection Algorithm
The module detects SSRs using a sliding window approach:

1. **Motif Identification**: Detects 1-6 bp repeating units
2. **Repeat Counting**: Counts consecutive repeats of each motif
3. **Threshold Filtering**: Keeps SSRs meeting minimum repeat thresholds
4. **Overlap Prevention**: Ensures no overlapping SSR calls
5. **Validation**: Checks for valid DNA sequences (ATCG only)

### 2. Genomic Region Classification
SSRs are classified into chloroplast structural regions:

- **LSC** (Large Single Copy): Main region with most genes
- **SSC** (Small Single Copy): Smaller region with specific genes
- **IRa/IRb** (Inverted Repeats): Duplicated regions flanking SSC

Classification based on feature coordinates from GenBank annotations.

### 3. Location Type Classification
SSRs are annotated by their genomic context:

- **Gene**: Within protein-coding genes (CDS features)
- **Intron**: Within introns of split genes
- **tRNA**: Within tRNA genes
- **rRNA**: Within ribosomal RNA genes (16S, 23S, 5S, 4.5S)
- **Intergenic**: Between annotated features

### 4. Gene Notation System
- Gene names are *italicized* following biological nomenclature
- Species names are *italicized* (e.g., *Arabidopsis thaliana*)
- Follows standard publication formatting

---

## SSR Detection Algorithm

### Motif Types and Thresholds

**Default Thresholds:**
| Motif Type | Size | Min Repeats | Example | Total Length |
|------------|------|-------------|---------|--------------|
| Mono | 1 bp | 10× | (A)₁₀ | ≥10 bp |
| Di | 2 bp | 5× | (AT)₅ | ≥10 bp |
| Tri | 3 bp | 4× | (ATG)₄ | ≥12 bp |
| Tetra | 4 bp | 3× | (ATGC)₃ | ≥12 bp |
| Penta | 5 bp | 3× | (ATGCA)₃ | ≥15 bp |
| Hexa | 6 bp | 3× | (ATGCAT)₃ | ≥18 bp |

### Detection Steps

1. **Sequence Scanning**
   ```python
   for position in genome:
       for motif_size in [1, 2, 3, 4, 5, 6]:
           extract_motif()
           count_consecutive_repeats()
           if repeats >= threshold:
               record_SSR()
   ```

2. **Overlap Prevention**
   - Longer SSRs take precedence
   - No overlapping SSR calls
   - Position-based conflict resolution

3. **Validation**
   - Valid DNA (ATCG only)
   - Minimum length requirements
   - Consistency checks

---

## Genomic Location Classification

### Region Detection Method

**Using GenBank Annotations:**
```
FEATURES             Location/Qualifiers
     source          1..154478
                     /organism="Arabidopsis thaliana"
     gene            complement(1..1512)
                     /gene="trnH-GUG"
                     /region="LSC"
```

### Classification Hierarchy

1. **Primary**: Check if SSR overlaps with annotated features
2. **Gene**: Matches CDS or gene feature
3. **Intron**: Within introns of split genes (join locations)
4. **tRNA/rRNA**: Specific RNA gene features
5. **Intergenic**: Between features (default)

### Inverted Repeat Handling

IRa and IRb are mirror images:
- Automatically classified based on coordinates
- Handles both standard and unusual IR configurations
- Accounts for IR junction variations

---

## Gene Notation System

### Formatting Rules

**Species Names:**
- *Italicized* in Excel output
- Full binomial nomenclature
- Example: *Arabidopsis thaliana*

**Gene Names:**
- *Italicized* when applicable
- Standard nomenclature (*atpF*, *rbcL*, *psbA*)
- Empty for intergenic SSRs

**Excel Formatting:**
- Applied using openpyxl
- Preserved when opening in Microsoft Excel
- May not display in Google Sheets

---

## R Visualization

### Automatic Figure Generation

If R is installed with required packages, three figures are automatically generated:

#### 1. SSR Distribution by Genomic Regions
- Bar plot showing SSR counts in LSC, SSC, IRa, IRb
- Grouped by species
- Shows structural variation

#### 2. SSR Motif Type Distribution
- Bar plot of mono-, di-, tri-, tetra-, penta-, hexa-nucleotide repeats
- Comparison across species
- Identifies dominant SSR types

#### 3. SSR Genomic Location Distribution
- Bar plot showing SSRs in genes, introns, tRNAs, rRNAs, intergenic
- Functional context analysis
- Comparative distribution

### Figure Specifications
- **Format**: PDF + PNG (600 DPI)
- **Size**: Optimized for publications
- **Colors**: Professional color schemes
- **Fonts**: Publication-ready typography

### Manual R Script Execution
```bash
# If you want to customize figures
cd Module12_SSR_Analysis/
Rscript generate_ssr_figures.R

# Edit the R script first if needed
nano generate_ssr_figures.R
```

---

## Troubleshooting

### Common Issues and Solutions

#### 1. No GenBank Files Found
```bash
❌ ERROR: No GenBank files found in directory
```
**Solution:**
```bash
# Check files are present
ls *.gb

# Verify file extensions
ls -la | grep -E "\.gb$|\.gbff$|\.genbank$|\.gbk$"

# Try with explicit path
python cgas_module12.py -i /full/path/to/genbank/files/
```

#### 2. R Not Found (Warning Only)
```bash
⚠ R not found. Skipping figure generation.
```
**Solution:**
- This is just a warning - Excel files are still created
- Install R if you want figures:
  ```bash
  # Ubuntu/Debian
  sudo apt-get install r-base
  
  # macOS
  brew install r
  
  # Then install R packages
  Rscript -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'reshape2'))"
  ```

#### 3. Wrong Threshold Format
```bash
❌ ERROR: Invalid threshold format
```
**Solution:**
```bash
# Use exactly 6 numbers, comma-separated (no spaces)
python cgas_module12.py -t 10,5,4,3,3,3  # ✓ CORRECT
python cgas_module12.py -t 10 5 4 3 3 3  # ✗ WRONG
python cgas_module12.py -t "10, 5, 4"    # ✗ WRONG
```

#### 4. Permission Denied
```bash
PermissionError: [Errno 13] Permission denied
```
**Solution:**
```bash
# Check folder permissions
ls -la output_folder/

# Use different output folder with write permissions
python cgas_module12.py -o ~/my_results/

# Or fix permissions (Linux/Mac)
chmod 755 Module12_SSR_Analysis/
```

#### 5. ModuleNotFoundError
```bash
ModuleNotFoundError: No module named 'openpyxl'
```
**Solution:**
```bash
# Install missing package
pip install openpyxl

# Or install all dependencies
pip install biopython pandas numpy openpyxl
```

### Quick Fixes Table

| Problem | Symptom | Quick Fix |
|---------|---------|-----------|
| No files found | "No GenBank files found" | Check file extensions (.gb, .gbff, etc.) |
| R missing | "R not found" warning | Install R or ignore (Excel still works) |
| Wrong thresholds | "Invalid format" error | Use `-t 10,5,4,3,3,3` format |
| Permission issues | "Permission denied" | Use `-o` for writable folder |
| Memory error | Script crashes | Process fewer files at once |
| Excel formatting | Species not italic | Open in Microsoft Excel (not Google Sheets) |

### Get Help
```bash
# Show all options
python cgas_module12.py --help

# Output:
# usage: cgas_module12.py [-h] [-i INPUT] [-o OUTPUT] [-t THRESHOLDS]
# 
# CGAS Module 12: Comprehensive Chloroplast SSR Analysis
# 
# options:
#   -h, --help            show this help message and exit
#   -i INPUT, --input INPUT
#                         Input directory containing GenBank files (default: current directory)
#   -o OUTPUT, --output OUTPUT
#                         Output directory (default: Module12_SSR_Analysis in input directory)
#   -t THRESHOLDS, --thresholds THRESHOLDS
#                         SSR detection thresholds (format: mono,di,tri,tetra,penta,hexa)
#                         Default: 10,5,4,3,3,3
```

---

## Examples

### Example 1: Complete Workflow
```bash
# 1. Prepare your data
mkdir -p /home/abdullah/chloroplast_data/
cp *.gb /home/abdullah/chloroplast_data/

# 2. Run analysis (simplest way)
cd /home/abdullah/chloroplast_data/
python cgas_module12.py

# 3. Check results
ls Module12_SSR_Analysis/
open Module12_SSR_Analysis/ssr_summary.xlsx
```

### Example 2: Sensitive Analysis for Detailed Study
```bash
# Want to find ALL possible SSRs?
python cgas_module12.py -t 8,4,3,3,3,3

# This finds:
# - Shorter mononucleotide repeats (8× instead of 10×)
# - Shorter dinucleotide repeats (4× instead of 5×)
# - Shorter trinucleotide repeats (3× instead of 4×)
# Perfect for detailed genome analysis
```

### Example 3: Strict Analysis for Marker Development
```bash
# Looking for reliable PCR markers?
python cgas_module12.py -t 12,6,5,4,4,4

# This finds:
# - Longer repeats (more stable)
# - Better for primer design
# - Fewer false positives
# Ideal for population genetics studies
```

### Example 4: Batch Processing Multiple Datasets
```bash
#!/bin/bash
# Process multiple folders

for dataset in dataset1 dataset2 dataset3; do
    echo "Processing $dataset..."
    python cgas_module12.py -i $dataset/ -o ${dataset}_ssr_results/
    echo "Done with $dataset"
done

echo "All datasets processed!"
```

### Example 5: Jupyter Notebook Analysis
```python
# In Jupyter notebook

# Cell 1: Run analysis
%run cgas_module12.py -i ./genomes/ -o ./analysis/

# Cell 2: Load results
import pandas as pd
summary = pd.read_excel('./analysis/ssr_summary.xlsx')
detailed = pd.read_excel('./analysis/all_ssrs_detailed.xlsx')

# Cell 3: Display with formatting
display(summary.style.background_gradient())

# Cell 4: Custom visualization
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.bar(summary['Species'], summary['Total_SSRs'])
plt.xticks(rotation=45, ha='right', fontstyle='italic')
plt.xlabel('Species', fontstyle='italic')
plt.ylabel('Total SSRs')
plt.title('Total SSRs per Species')
plt.tight_layout()
plt.show()
```

---

## FAQ

### Q1: What's the simplest way to run this?
**A:** Just `python cgas_module12.py` in your data folder. It automatically creates everything.

### Q2: What if I don't have R installed?
**A:** No problem! All Excel files are created. Only the PDF/PNG figures are skipped.

### Q3: Can I change where results are saved?
**A:** Yes! Use `-o` flag: `python cgas_module12.py -o my_results/`

### Q4: What thresholds should I use?
**A:**
- **Default (10,5,4,3,3,3)**: Good for most chloroplast studies and publications
- **Sensitive (8,4,3,3,3,3)**: Find more SSRs, good for detailed genomic analysis
- **Strict (12,6,5,4,4,4)**: Fewer but more reliable SSRs, good for marker development

### Q5: How long does it take?
**A:** Typically seconds to minutes per file, depending on genome size (~150kb = 5-10 seconds).

### Q6: Can I run this in Jupyter?
**A:** Yes! Use `%run cgas_module12.py` in a Jupyter cell.

### Q7: What file formats are supported?
**A:** .gb, .gbff, .genbank, .gbk (standard GenBank formats from NCBI)

### Q8: Why aren't species names italic in my Excel file?
**A:** Italics work in Microsoft Excel but may not display in Google Sheets or LibreOffice. Open with Excel 2010+.

### Q9: Can I analyze multiple species together?
**A:** Yes! Put all GenBank files in one folder and run the script once. Results are automatically organized.

### Q10: How do I cite this tool?
**A:** 
```
Abdullah. (2026). CGAS Module 12: Comprehensive Chloroplast SSR Analysis.
Version 1.0.1. Chloroplast Genome Analysis Suite (CGAS).
```

### Q11: What Python version do I need?
**A:** Python 3.7 or higher. Tested on Python 3.7-3.12.

### Q12: Can I customize the R figures?
**A:** Yes! Edit the `generate_ssr_figures.R` file in the output folder and re-run it with `Rscript generate_ssr_figures.R`.

---

## Technical Specifications

### Performance
- **File size**: Handles typical chloroplast genomes (~150-200kb)
- **Speed**: ~5-15 seconds per file (depending on genome size and SSR density)
- **Memory**: Minimal (processes files sequentially, ~100-200MB RAM)
- **Output**: Professional Excel files with rich formatting

### Compatibility
- **Python**: 3.9+ (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating Systems**: Windows 10/11, macOS 10.15+, Linux (Ubuntu 18.04+)
- **Excel**: Compatible with Microsoft Excel 2010+ (.xlsx format)
- **R**: 4.0+ (optional, for figures)

### Limits
- **Max files**: No practical limit (processes sequentially)
- **Max genome size**: Tested up to 500kb (works with larger genomes)
- **SSR length**: Theoretically unlimited (tested up to 100+ bp repeats)
- **Excel rows**: Up to 1,048,576 (Excel format limit)
- **Species**: No limit (tested with 100+ species)

### Quality Features
- ✅ No overlapping SSRs (automatic conflict resolution)
- ✅ Valid motif checking (ATCG nucleotides only)
- ✅ Trans-splicing detection in split genes
- ✅ IR region handling (IRa/IRb classification)
- ✅ Professional Excel formatting (italics, colors, borders)
- ✅ Comprehensive error handling and recovery
- ✅ Publication-ready outputs

---

## References

### GenBank Format
- [NCBI GenBank Format Guide](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)
- [INSDC Feature Table Documentation](https://www.insdc.org/)

### SSR Analysis Background
- **Chloroplast SSRs**: Valuable markers for plant phylogenetics and population genetics
- **Microsatellites**: Tandem repeats used extensively in genetic marker development
- **SSR Evolution**: Studies on repeat instability and evolutionary dynamics
- **Marker Development**: Guidelines for selecting SSRs for molecular markers

### Related SSR Tools
- **MISA** (MicroSatellite identification tool): Perl-based SSR detector
- **CPStool** A tool detect SSRs from the genbank file. 

### Python Libraries
- **BioPython**: [Official Documentation](https://biopython.org/)
- **pandas**: [User Guide](https://pandas.pydata.org/docs/)
- **openpyxl**: [Documentation](https://openpyxl.readthedocs.io/)
- **NumPy**: [Reference](https://numpy.org/doc/)

### R Packages (Optional)
- **ggplot2**: [Documentation](https://ggplot2.tidyverse.org/)
- **dplyr**: [Documentation](https://dplyr.tidyverse.org/)
- **tidyr**: [Documentation](https://tidyr.tidyverse.org/)

### CGAS Suite
- **CGAS**: Chloroplast Genome Analysis Suite - integrated tools for comprehensive chloroplast genomics

### Citation
If using CGAS Module 12 in publications, please cite:
```
Abdullah, Yan, R., & Tian, X. (2025). CGAS (Chloroplast Genome Analysis Suite): 
An Automated Python Pipeline for Comprehensive Comparative Chloroplast Genomics. 
bioRxiv. https://doi.org/10.64898/2025.12.21.695765
```

---

## Support & Contact

### Getting Help
```bash
# 1. First check built-in help
python cgas_module12.py --help

# 2. Review this documentation
# 3. Check example outputs in Examples section
# 4. Verify dependencies are installed
pip list | grep -E "biopython|pandas|openpyxl|numpy"
```

### Common Issues Solved Here
- ✅ File not found? Use `-i` flag with correct path
- ✅ No output? Check file extensions (.gb, .gbff, .genbank, .gbk)
- ✅ R warnings? Install R or ignore (Excel files still work)
- ✅ Formatting issues? Open in Microsoft Excel (not Google Sheets)
- ✅ Wrong thresholds? Use format: `-t 10,5,4,3,3,3`

### Quick Reference Card
```bash
# ⚡ QUICK START ⚡
cgas --module 12                                 # cgas command
cgas-ssr                                         # shortcut command
python cgas_module12.py                          # python command
python cgas_module12.py -i data/                 # Specify input
python cgas_module12.py -o results/              # Specify output
python cgas_module12.py -t 8,4,3,3,3,3           # Sensitive detection

# 🔬 JUPYTER NOTEBOOK 🔬
%run cgas_module12.py                            # %run works for Module 12
!cgas --module 12                                # ! also works
!cgas-ssr -t 8,4,3,3,3,3                         # Custom thresholds in notebook

# 📊 OUTPUT 📊
# Module12_SSR_Analysis/
# ├── all_ssrs_detailed.xlsx      # Complete catalog
# ├── ssr_summary.xlsx             # Statistics
# ├── ssr_by_species.xlsx          # Per-species sheets
# └── Figures/                     # PDF + PNG (if R installed)

# 🎯 TIPS 🎯
# - GenBank files: .gb, .gbff, .genbank, .gbk
# - Install R for automatic figure generation
# - Adjust thresholds based on your research needs
# - Species/gene names are italicized (Excel 2010+)
# - Default thresholds: 10,5,4,3,3,3 (good for most studies)
```

---

**Last Updated**: January 2026  
**Version**: 1.0.1  
**Author**: Abdullah  
**Part of**: CGAS (Chloroplast Genome Analysis Suite)

---

*Happy SSR Analysis! 🧬✨*
