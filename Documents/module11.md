# CGAS Module 11: Gene and tRNA Intron Extraction and Analysis
## Complete Documentation and User Guide

---

## 📋 Table of Contents
1. [Introduction](#introduction)
2. [Quick Start - Run It Now!](#quick-start---run-it-now)
3. [Command-Line Usage with Examples](#command-line-usage-with-examples)
4. [Jupyter Notebook Usage](#jupyter-notebook-usage)
5. [Input Requirements](#input-requirements)
6. [Output Structure](#output-structure)
7. [Troubleshooting](#troubleshooting)
8. [Examples](#examples)
9. [FAQ](#faq)
10. [Technical Specifications](#technical-specifications)
11. [References](#references)

---

## Quick Start - Run It Now!

### ⚡ One Command Does Everything!

```bash
# Option 1: Using cgas command
cgas --module 11

# Option 2: Using cgas-intron shortcut
cgas-intron

# Option 3: Using python directly
python cgas_module11.py
```

**What happens when you run this simple command:**
1. ✅ Searches for GenBank files in your **current directory**
2. ✅ Automatically creates `Module11_Intron_Analysis` folder
3. ✅ Processes all `.gb`, `.gbff`, `.genbank`, `.gbk` files
4. ✅ Extracts introns from both protein-coding genes AND tRNA genes
5. ✅ Groups results by species with professional formatting
6. ✅ Creates comparison sheets if multiple species present
7. ✅ Saves timestamped Excel file with all results

### 📁 Real-World Folder Examples

#### Example 1: Basic Usage (Easiest!)
```bash
# Your folder structure:
abdullah/
├── chloroplast_project/
│   ├── arabidopsis.gb
│   ├── oryza.gb
│   └── maize.gb

# Navigate to your data folder
cd /home/abdullah/chloroplast_project/

# Run the script (no arguments needed!)
python cgas_module11.py

# Output created automatically:
# Module11_Intron_Analysis/
# └── intron_data_20240115_143022.xlsx
#     ├── Gene Introns sheet
#     ├── tRNA Introns sheet
#     ├── Gene Intron Size Comparison (if multiple species)
#     └── tRNA Intron Size Comparison (if multiple species)
```

#### Example 2: Specify Input Folder
```bash
# Your folder structure:
/home/abdullah/data_analysis/
├── genbank_files/
│   ├── species1.gb
│   ├── species2.gb
│   └── species3.gb

# Run from ANYWHERE with -i flag
python cgas_module11.py -i /home/abdullah/data_analysis/genbank_files/

# Output created in:
# /home/abdullah/data_analysis/genbank_files/Module11_Intron_Analysis/
```

#### Example 3: Custom Input AND Output Folders
```bash
# You control both input and output locations
python cgas_module11.py -i abdullah/data_analysis -o intron_results/

# What this does:
# 1. Reads from: abdullah/data_analysis/*.gb
# 2. Saves to: intron_results/ (instead of default folder)
```

### 🎯 Quick Start Cheat Sheet

```bash
# ====== cgas command ======
cgas --module 11                                 # Process current directory
cgas --module 11 -i your_data_folder/            # Specify input folder
cgas --module 11 -i abdullah/data_analysis -o intron_results/  # Custom I/O

# ====== cgas-intron shortcut ======
cgas-intron                                      # Process current directory
cgas-intron -i your_data_folder/                 # Specify input folder
cgas-intron -i chloroplast_data/ -o final_results/  # Custom I/O

# ====== python command ======
python cgas_module11.py                          # Process current directory
python cgas_module11.py -i your_data_folder/     # Specify input folder
python cgas_module11.py -i abdullah/data_analysis -o intron_results/  # Custom I/O
python cgas_module11.py --help                   # Get help
```

### 📊 What You Get (Output Files)
After running ANY of the above commands, you'll get:

```
Module11_Intron_Analysis/  # Or your custom output folder
└── 📊 intron_data_YYYYMMDD_HHMMSS.xlsx
    ├── 📄 Gene Introns               # All gene/CDS introns
    ├── 📄 tRNA Introns               # All tRNA introns
    ├── 📈 Gene Intron Size Comparison # Stats if multiple species
    └── 📈 tRNA Intron Size Comparison # Stats if multiple species
```

---

## Introduction

**CGAS Module 11** is a specialized tool for extracting and analyzing intron positions and lengths from chloroplast genomes. It processes GenBank format files to identify introns in both **protein-coding genes (CDS/gene features)** and **tRNA genes**, producing publication-ready Excel reports with comprehensive comparative analysis.

### Key Applications:
- **Evolutionary Studies**: Compare intron size variation across species
- **Genome Annotation**: Verify intron-exon boundaries
- **Comparative Genomics**: Analyze intron conservation patterns
- **Publication**: Generate formatted tables for scientific papers

### Why Use This Tool?
1. **Dual Analysis**: Extracts from both genes AND tRNAs in one run
2. **Smart Deduplication**: Handles duplicate genes in Inverted Repeat regions
3. **Roman Numerals**: Professional notation for multiple introns (Intron I, II, III)
4. **Publication-Ready**: Excel files with italic species/gene names
5. **Comparative Statistics**: Shows min/max/mean intron sizes across species

---

## Jupyter Notebook Usage

> **Note:** Module 11 is a pure Python script, so both `%run` and the `!` shell operator work. Use `%run` for direct script execution, or `!` for the `cgas` / `cgas-intron` commands. Run the cell from the directory containing your GenBank files, or pass `-i` to specify the input path.

```python
# Using %run (executes the script directly)
%run cgas_module11.py
%run cgas_module11.py -i genbank_files/
%run cgas_module11.py -i data/ -o results/

# Using ! operator with cgas command
!cgas --module 11
!cgas --module 11 -i genbank_files/
!cgas --module 11 -i data/ -o results/

# Using ! operator with cgas-intron shortcut
!cgas-intron
!cgas-intron -i genbank_files/
!cgas-intron -i data/ -o results/

# Using ! operator with python
!python cgas_module11.py
!python cgas_module11.py -i genbank_files/ -o results/
```

### Advanced: Display Results in Notebook
```python
# After running, display the Excel output
import pandas as pd
from IPython.display import display

# Load the Excel file (use your actual timestamp)
excel_file = './Module11_Intron_Analysis/intron_data_20240115_143022.xlsx'
xls = pd.ExcelFile(excel_file)

# Display each sheet
for sheet_name in xls.sheet_names:
    print(f"\n=== {sheet_name} ===")
    df = pd.read_excel(xls, sheet_name=sheet_name)
    display(HTML(df.head(10).to_html(classes='table table-striped')))
```

---

## Command-Line Usage with Examples

### Complete Command Reference

```bash
# ====================================================================
# USING cgas COMMAND
# ====================================================================

# Process current directory
cgas --module 11

# Specify input directory
cgas --module 11 -i /home/user/chloroplast_genomes/

# Custom input and output directories
cgas --module 11 -i abdullah/data_analysis -o intron_results/
```

```bash
# ====================================================================
# USING cgas-intron SHORTCUT
# ====================================================================

# Process current directory
cgas-intron

# With specific input and output
cgas-intron -i abdullah/data_analysis -o intron_results/
```

```bash
# ====================================================================
# USING python COMMAND
# ====================================================================

# ====================================================================
# BASIC USAGE
# ====================================================================

# 1. Run in current directory (simplest!)
python cgas_module11.py

# 2. Specify input folder
python cgas_module11.py -i genbank_files/

# 3. Specify BOTH input and output folders
python cgas_module11.py -i abdullah/data_analysis -o intron_results/

# 4. Get help
python cgas_module11.py --help


# ====================================================================
# REAL-WORLD EXAMPLES WITH PATHS
# ====================================================================

# Example 1: Process data in "chloroplast_genomes" folder
python cgas_module11.py -i /home/user/chloroplast_genomes/

# Example 2: Save results to custom "analysis_results" folder
python cgas_module11.py -i data/ -o analysis_results/

# Example 3: Run from different directory
cd /home/abdullah/scripts/
python cgas_module11.py -i ../data/genbank/ -o ../results/introns/

# Example 4: Process specific set of species
python cgas_module11.py -i monocots/ -o monocot_intron_analysis/
python cgas_module11.py -i dicots/ -o dicot_intron_analysis/

# Example 5: Windows users
python cgas_module11.py -i "C:\Users\abdullah\chloroplast_data" -o "C:\Results\introns"
```

---

## Input Requirements

### Supported File Formats
- `.gb` (GenBank)
- `.gbff` (GenBank Flat File)
- `.genbank`
- `.gbk` (GenBank)

### GenBank Feature Requirements
For intron detection, GenBank files must contain:

1. **Gene Features** with introns:
```
gene            join(2001..2500,2701..3200)
                /gene="atpF"
```

2. **CDS Features** with introns:
```
CDS             join(2001..2500,2701..3200)
                /gene="atpF"
                /product="ATP synthase subunit I"
```

3. **tRNA Features** with introns:
```
tRNA            join(5001..5037,5078..5149)
                /gene="trnV-UAC"
                /product="tRNA-Val"
```

### File Organization
```
your_data_folder/
├── Arabidopsis_thaliana.gb      # Species name in filename recommended
├── Oryza_sativa.gb
├── Zea_mays.gb
└── outgroup_species.gb
```

**Note**: Species names are extracted from the `/organism=` qualifier in GenBank files, not filenames.

---

## Output Structure

### Directory Structure
```
Module11_Intron_Analysis/
└── intron_data_20260128_154523.xlsx    # Timestamped filename

# Custom output folder example:
my_custom_results/
└── intron_data_20260128_154523.xlsx
```

### Excel File Structure

#### For Single Species (2 sheets):
1. **Gene Introns** - All protein-coding gene introns
2. **tRNA Introns** - All tRNA introns

#### For Multiple Species (4 sheets):
1. **Gene Introns** - All protein-coding gene introns
2. **tRNA Introns** - All tRNA introns
3. **Gene Intron Size Comparison** - Statistical comparison across species
4. **tRNA Intron Size Comparison** - Statistical comparison across species

### Sheet Column Definitions

#### Gene Introns & tRNA Introns Sheets:
| Column | Description | Example |
|--------|-------------|---------|
| **Species** | Organism name (*italicized*) | *Arabidopsis thaliana* |
| **Gene** | Gene name (*italicized*) | *atpF*, *trnV-UAC* |
| **Intron** | Intron identifier | Intron I, Intron II |
| **Intron Size (bp)** | Length in base pairs | 842 |
| **Intron Position** | Start-End coordinates | 2501-3342 |

#### Comparison Sheets:
| Column | Description | Example |
|--------|-------------|---------|
| **Gene** | Gene name (*italicized*) | *atpF* |
| **Min (bp)** | Minimum size across species | 735 |
| **Max (bp)** | Maximum size across species | 945 |
| **Mean (bp)** | Average size | 842.3 |
| **Species Count** | Number of species with this gene | 15 |

---

## Troubleshooting

### Common Issues and Solutions

#### Problem 1: "No GenBank files found"
**Cause**: Script can't find `.gb`, `.gbff`, `.genbank`, or `.gbk` files

**Solution**:
```bash
# Check files are present
ls *.gb

# Verify file extensions
ls -la

# Try with explicit path
python cgas_module11.py -i /full/path/to/genbank/files/
```

#### Problem 2: "No introns found in any files"
**Cause**: GenBank features don't contain `join()` locations

**Solution**:
- Open GenBank file in text editor
- Search for `join(` 
- If not found, genome may have no introns
- Example valid entry:
```
gene            join(2001..2500,2701..3200)
                /gene="atpF"
```

#### Problem 3: "ModuleNotFoundError: No module named 'openpyxl'"
**Cause**: Missing Python package

**Solution**:
```bash
pip install openpyxl
# Or install all dependencies
pip install biopython pandas numpy openpyxl
```

#### Problem 4: Excel file won't open
**Cause**: File may be locked by another process

**Solution**:
- Close Excel completely
- Try reopening the file
- Check file permissions
- Re-run the script

#### Problem 5: Permission denied when creating output folder
**Cause**: Insufficient permissions

**Solution**:
```bash
# Use custom output folder with write permissions
python cgas_module11.py -o ~/intron_results/

# Or change permissions (Linux/Mac)
chmod 755 Module11_Intron_Analysis/
```

#### Problem 6: Very large memory usage
**Cause**: Processing many large GenBank files

**Solution**:
- Process files in batches
- Split large GenBank files
- Use `-i` to process one folder at a time

### Quick Fixes Table

| Problem | Symptom | Quick Fix |
|---------|---------|-----------|
| No files found | "No GenBank files found" | Check file extensions (.gb, .gbff, etc.) |
| No introns found | "No introns found" | Verify `join()` in GenBank features |
| Missing packages | "ModuleNotFoundError" | `pip install biopython pandas openpyxl` |
| Permission issues | "Permission denied" | Use `-o` for custom writable folder |
| Excel issues | File won't open | Close Excel, check file isn't locked |

### Get Detailed Help
```bash
# Show all options and examples
python cgas_module11.py --help

# Output:
# usage: cgas_module11.py [-h] [-i INPUT] [-o OUTPUT]
# 
# CGAS Module 11: Gene and tRNA Intron Extraction
# 
# options:
#   -h, --help            show this help message and exit
#   -i INPUT, --input INPUT
#                         Input directory containing GenBank files (default: current directory)
#   -o OUTPUT, --output OUTPUT
#                         Output directory (default: Module11_Intron_Analysis in input directory)
# 
# Examples:
#   python cgas_module11.py                          # Process current directory
#   python cgas_module11.py -i genbank_files/        # Process specific directory
#   python cgas_module11.py -i data/ -o results/     # Custom input and output
```

---

## Examples

### Example 1: Complete Workflow
```bash
# 1. Prepare your data
mkdir -p /home/abdullah/chloroplast_introns/
cp *.gb /home/abdullah/chloroplast_introns/

# 2. Run analysis (simplest way)
cd /home/abdullah/chloroplast_introns/
python cgas_module11.py

# 3. Check results
ls Module11_Intron_Analysis/
# intron_data_20260128_143022.xlsx

# 4. Open in Excel
# Shows formatted sheets with all introns
```

### Example 2: Analyzing Specific Dataset
```bash
# Analyze only monocot species
python cgas_module11.py -i /data/monocots/ -o monocot_introns/

# Analyze only dicot species
python cgas_module11.py -i /data/dicots/ -o dicot_introns/

# Compare results from different groups
```

### Example 3: Jupyter Notebook Analysis
```python
# In Jupyter notebook

# Cell 1: Run analysis
%run cgas_module11.py -i ./genomes/ -o ./intron_analysis/

# Cell 2: Load and display results
import pandas as pd
import glob

# Find the latest output file
output_files = glob.glob('./intron_analysis/intron_data_*.xlsx')
latest_file = sorted(output_files)[-1] if output_files else None

if latest_file:
    # Read Excel
    xls = pd.ExcelFile(latest_file)
    
    # Display each sheet
    for sheet in xls.sheet_names:
        print(f"\n{'='*50}")
        print(f"Sheet: {sheet}")
        print('='*50)
        df = pd.read_excel(xls, sheet_name=sheet)
        display(df.head())
```

### Example 4: Batch Processing
```bash
#!/bin/bash
# Process multiple datasets

datasets=("set1" "set2" "set3" "set4")

for dataset in "${datasets[@]}"; do
    echo "Processing $dataset..."
    python cgas_module11.py -i "data/$dataset/" -o "results/${dataset}_introns/"
    echo "✓ Completed $dataset"
done

echo "All datasets processed!"
```

### Example 5: Checking Results
```bash
# After running, check what you got
ls -la Module11_Intron_Analysis/

# Check file size
du -h Module11_Intron_Analysis/*.xlsx

# Count sheets in Excel (requires pandas)
python -c "
import pandas as pd
import glob
files = glob.glob('Module11_Intron_Analysis/intron_data_*.xlsx')
if files:
    xls = pd.ExcelFile(files[0])
    print(f'Sheets: {xls.sheet_names}')
    print(f'Number of sheets: {len(xls.sheet_names)}')
"
```

---

## FAQ

### Q1: What's the simplest way to run this?
**A:** Just `python cgas_module11.py` in your data folder. It automatically creates everything.

### Q2: What file formats are supported?
**A:** .gb, .gbff, .genbank, .gbk (standard GenBank formats from NCBI)

### Q3: Can I change where results are saved?
**A:** Yes! Use `-o` flag: `python cgas_module11.py -o my_results/`

### Q4: What if I only have one species?
**A:** You get 2 sheets: Gene Introns and tRNA Introns (no comparison sheets).

### Q5: What if I have multiple species?
**A:** You get 4 sheets: two data sheets + two comparison sheets with statistics.

### Q6: How are duplicate genes in IR regions handled?
**A:** Automatically deduplicated. Each gene-intron combination counted once per species.

### Q7: Can I run this in Jupyter?
**A:** Yes! Use `%run cgas_module11.py` in a Jupyter cell.

### Q8: Why Roman numerals for multiple introns?
**A:** Standard biological notation (e.g., "*Intron I*", "*Intron II*") used in publications.

### Q9: What's the intron length limit?
**A:** 1-15,000 bp. This covers all known chloroplast introns.

### Q10: Can I analyze nuclear genes?
**A:** The script works for any GenBank file, but it's optimized for chloroplast genomes. Nuclear introns may have different characteristics and could exceed the 15kb limit.

### Q11: Why are gene and species names italicized?
**A:** Following biological nomenclature standards for publication-ready output.

### Q12: What Python version do I need?
**A:** Python 3.7 or higher. Tested on Python 3.7-3.12.

---

## Technical Specifications

### Performance
- **File size**: Handles typical chloroplast genomes (~150-200kb)
- **Speed**: ~5 seconds per file (depending on complexity)
- **Memory**: Minimal (processes files sequentially, ~50-100MB RAM)
- **Output**: Professional Excel files with rich formatting

### Compatibility
- **Python**: 3.9+ (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating Systems**: Windows 10/11, macOS 10.15+, Linux (Ubuntu 18.04+)
- **Excel**: Compatible with Microsoft Excel 2010+ (.xlsx format)
- **GenBank**: NCBI GenBank flat file format (standard)

### Limits
- **Max intron length**: 15,000 bp
- **Min intron length**: 1 bp
- **Excel rows**: Up to 1,048,576 (Excel limit)
- **Files processed**: No practical limit (sequential processing)
- **Species**: No limit (tested with 100+ species)

### Quality Features
- ✅ Validates intron lengths (1-15,000 bp range check)
- ✅ Handles IR region duplicates automatically
- ✅ Professional Excel formatting (italics, colors, borders)
- ✅ Roman numeral notation for multiple introns
- ✅ Comparative statistics for multiple species
- ✅ Comprehensive error handling and recovery
- ✅ Timestamped outputs to prevent overwriting

---

## References

### GenBank Format
- [NCBI GenBank Format Guide](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)
- [INSDC Feature Table Documentation](https://www.insdc.org/submitting-standards/feature-table/)

### Python Libraries
- **BioPython**: [Official Documentation](https://biopython.org/)
- **pandas**: [User Guide](https://pandas.pydata.org/docs/)
- **openpyxl**: [Documentation](https://openpyxl.readthedocs.io/)
- **NumPy**: [Reference](https://numpy.org/doc/)

### Related Tools
- **CGAS Suite**: Complete chloroplast genome analysis suite
- **DOGMA**: Dual Organellar Genome Annotator
- **GeSeq**: Genome Sequencing annotation tool
- **Geneious**: Commercial genome analysis software
- **ApE**: Free plasmid editor with GenBank support

### Citation
If using CGAS Module 11 in publications, please cite:
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
python cgas_module11.py --help

# 2. Review this documentation
# 3. Check example outputs in Examples section
# 4. Verify dependencies are installed
pip list | grep -E "biopython|pandas|openpyxl"
```

### Quick Reference Card
```bash
# ⚡ QUICK START ⚡
cgas --module 11                                 # cgas command
cgas-intron                                      # shortcut command
python cgas_module11.py                          # python command
python cgas_module11.py -i data/                 # Specify input
python cgas_module11.py -o results/              # Specify output
python cgas_module11.py -i data/ -o out/         # Both custom

# 🔬 JUPYTER NOTEBOOK 🔬
%run cgas_module11.py                            # %run works for Module 11
!cgas --module 11                                # ! also works
!cgas-intron -i data/                            # Shortcut in notebook

# 📊 OUTPUT 📊
# Module11_Intron_Analysis/
# └── intron_data_YYYYMMDD_HHMMSS.xlsx
#     ├── Gene Introns
#     ├── tRNA Introns
#     ├── Gene Intron Size Comparison (if multiple species)
#     └── tRNA Intron Size Comparison (if multiple species)

# 🎯 TIPS 🎯
# - GenBank files need .gb, .gbff, .genbank, or .gbk extension
# - Check for `join()` in GenBank features for intron detection
# - Use -o for custom output folder location
# - Multiple species automatically generates comparison sheets
# - Gene and species names are italicized (publication standard)
```

---

**Last Updated**: January 2026  
**Version**: 1.0.1  
**Author**: Abdullah  
**Part of**: CGAS (Chloroplast Genome Analysis Suite)

---

*Happy Intron Analysis! 🧬✨*
