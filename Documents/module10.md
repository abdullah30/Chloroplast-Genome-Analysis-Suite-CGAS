# CGAS Module 10: SNP/Substitution Analysis with Graphical Visualization
## Complete Documentation and User Guide

---

## 📋 Table of Contents
1. [Introduction](#introduction)
2. [Quick Start - Just Run It!](#quick-start---just-run-it)
3. [Command-Line Usage with Examples](#command-line-usage-with-examples)
4. [Jupyter Notebook Usage](#jupyter-notebook-usage)
5. [Input Requirements](#input-requirements)
6. [Output Structure](#output-structure)
7. [Substitution Types](#substitution-types)
8. [Transition vs Transversion](#transition-vs-transversion)
9. [R Visualization](#r-visualization)
10. [Troubleshooting](#troubleshooting)
11. [Examples](#examples)
12. [FAQ](#faq)
13. [Technical Specifications](#technical-specifications)
14. [References](#references)

---

## Introduction

**CGAS Module 10** is a comprehensive tool for analyzing nucleotide substitutions (SNPs) in aligned chloroplast genome sequences. This module identifies all types of nucleotide changes, calculates transition/transversion ratios, and generates publication-quality visualizations.

This module performs complete SNP analysis with:
- **Pairwise sequence comparison**: Analyzes aligned FASTA files
- **12 substitution types**: Identifies all possible nucleotide changes (A→G, G→A, etc.)
- **Transition/Transversion analysis**: Calculates Ts/Tv ratios
- **Position tracking**: Records exact positions of each substitution
- **Excel reports**: Individual and merged analysis spreadsheets
- **R visualizations**: Publication-quality graphs (PDF + PNG, 600 DPI)

### Key Features:
- **Automatic detection**: Finds all FASTA files in directory
- **Comprehensive analysis**: All 12 nucleotide substitution types
- **Ts/Tv ratio**: Critical metric for evolutionary studies
- **Position tracking**: Exact coordinates of each substitution
- **Batch processing**: Analyzes multiple alignments at once
- **Merged analysis**: Comparative statistics across all samples
- **Professional visualizations**: Grouped bar charts and ratio plots

### Scientific Applications:
- **Molecular Evolution**: Study nucleotide substitution patterns
- **Phylogenetic Analysis**: Compare substitution rates across species
- **Population Genetics**: Identify polymorphic sites
- **Selection Pressure**: Detect synonymous vs. non-synonymous changes
- **Comparative Genomics**: Analyze divergence between species
- **Mutation Rate Studies**: Quantify substitution frequencies

### Part of CGAS (Chloroplast Genome Analysis Suite)
CGAS is a comprehensive suite of tools designed for chloroplast genome analysis, providing seamless integration across different analysis modules.

---

## Quick Start - Just Run It!

### ⚡ One Command Does It All!

```bash
# Option 1: Using cgas command
cgas --module 10

# Option 2: Using cgas-snp shortcut
cgas-snp

# Option 3: Using python directly
python cgas_module10.py
```

**What happens when you run this:**
1. ✅ Finds all FASTA alignment files in current directory
2. ✅ Analyzes each pairwise alignment for substitutions
3. ✅ Creates individual Excel reports per alignment
4. ✅ Generates merged comparative analysis
5. ✅ Produces CSV file for R visualization
6. ✅ Creates publication-quality graphs (if R installed)

### 📁 Real-World Folder Examples

#### Example 1: Basic Usage (Easiest!)
```bash
# Your folder structure:
/home/abdullah/snp_analysis/
├── aligned_sequences/
│   ├── rbcL_aligned.fasta
│   ├── matK_aligned.fasta
│   └── ndhF_aligned.fasta

# Navigate to your folder
cd /home/abdullah/snp_analysis/aligned_sequences/

# Run the module (no arguments needed!)
python cgas_module10.py

# Output created automatically:
# Module10_SNP_Analysis/
# ├── rbcL_aligned_Substitutions.xlsx
# ├── matK_aligned_Substitutions.xlsx
# ├── ndhF_aligned_Substitutions.xlsx
# ├── Complete_Substitution_Analysis.xlsx
# ├── Complete_Substitution_Analysis.csv
# └── Figures/
#     ├── Substitution_Types.pdf
#     ├── Substitution_Types.png
#     ├── Ts_Tv_Ratio.pdf
#     └── Ts_Tv_Ratio.png
```

#### Example 2: Specify Input Folder
```bash
# Process alignments from specific directory
python cgas_module10.py -i alignments/

# Output created in:
# alignments/Module10_SNP_Analysis/
```

#### Example 3: Custom Input and Output
```bash
# Full control over directories
python cgas_module10.py -i aligned_genes/ -o SNP_Results/

# Input from: aligned_genes/
# Output to: SNP_Results/
```

#### Example 4: Skip Visualizations
```bash
# Generate only Excel files (no R graphs)
python cgas_module10.py --no-figures

# Useful if R is not installed
```

### 🎯 Quick Start Cheat Sheet

```bash
# ====== cgas command ======
cgas --module 10                                 # Process current directory
cgas --module 10 -i alignments/                  # Specify input folder
cgas --module 10 -i data/ -o results/            # Custom input and output
cgas --module 10 --no-figures                    # Skip R visualizations

# ====== cgas-snp shortcut ======
cgas-snp                                         # Process current directory
cgas-snp -i alignments/                          # Specify input folder
cgas-snp -i data/ -o results/                    # Custom input and output
cgas-snp --no-figures                            # Skip R visualizations

# ====== python command ======
python cgas_module10.py                          # Process current directory
python cgas_module10.py -i alignments/           # Specify input folder
python cgas_module10.py -i data/ -o results/     # Custom input and output
python cgas_module10.py --no-figures             # Skip R visualizations
python cgas_module10.py --help                   # Get help
```

### 📊 What You Get (Output Files)

```
Module10_SNP_Analysis/                    # Created automatically
├── 📊 Individual Reports (one per alignment)
│   ├── rbcL_aligned_Substitutions.xlsx    # Complete substitution analysis
│   ├── matK_aligned_Substitutions.xlsx
│   └── ndhF_aligned_Substitutions.xlsx
│
├── 📈 Merged Analysis
│   ├── Complete_Substitution_Analysis.xlsx  # All alignments combined
│   └── Complete_Substitution_Analysis.csv   # For R visualization
│
└── 🖼️  Figures/ (if R installed)
    ├── Substitution_Types.pdf             # Publication-quality (vector)
    ├── Substitution_Types.png             # High-resolution (600 DPI)
    ├── Ts_Tv_Ratio.pdf                    # Ts/Tv comparison
    └── Ts_Tv_Ratio.png                    # High-resolution (600 DPI)
```
---

## Command-Line Usage with Examples

### Complete Command Reference

```bash
# ====================================================================
# USING cgas COMMAND
# ====================================================================

# Process current directory
cgas --module 10

# Specify input directory
cgas --module 10 -i /home/abdullah/alignments/

# Custom input and output directories
cgas --module 10 -i aligned_genes/ -o SNP_Results/

# Skip R visualization
cgas --module 10 --no-figures
```

```bash
# ====================================================================
# USING cgas-snp SHORTCUT
# ====================================================================

# Process current directory
cgas-snp

# With specific input and output
cgas-snp -i aligned_genes/ -o SNP_Results/

# Skip R visualization
cgas-snp --no-figures
```

```bash
# ====================================================================
# USING python COMMAND
# ====================================================================

# ====================================================================
# BASIC USAGE
# ====================================================================

# 1. Process current directory (simplest!)
python cgas_module10.py

# 2. Specify input folder
python cgas_module10.py -i alignments/

# 3. Custom input and output folders
python cgas_module10.py -i aligned_genes/ -o SNP_Results/

# 4. Skip R visualization
python cgas_module10.py --no-figures

# 5. Get help
python cgas_module10.py --help


# ====================================================================
# REAL-WORLD EXAMPLES
# ====================================================================

# Example 1: Process gene alignments
python cgas_module10.py -i gene_alignments/

# Example 2: Save to specific output folder
python cgas_module10.py -i data/ -o ../results/snps/

# Example 3: No R visualization (faster)
python cgas_module10.py -i alignments/ --no-figures

# Example 4: Windows users
python cgas_module10.py -i "C:\Users\abdullah\alignments"
```

### Parameter Details

| Parameter | Short | Default | Description |
|-----------|-------|---------|-------------|
| `--input` | `-i` | `.` (current) | Input directory with FASTA alignments |
| `--output` | `-o` | `Module10_SNP_Analysis` | Output directory for results |
| `--no-figures` | - | False | Skip R visualization generation |

---

## Jupyter Notebook Usage

> **Note:** Module 10 is a pure Python script, so both `%run` and the `!` shell operator work. Use `%run` for direct script execution, or `!` for the `cgas` / `cgas-snp` commands. Run the cell from the directory containing your FASTA alignment files, or pass `-i` to specify the input path.

```python
# Using %run (executes the script directly)
%run cgas_module10.py
%run cgas_module10.py -i alignments/
%run cgas_module10.py -i data/ -o results/
%run cgas_module10.py --no-figures

# Using ! operator with cgas command
!cgas --module 10
!cgas --module 10 -i alignments/
!cgas --module 10 --no-figures

# Using ! operator with cgas-snp shortcut
!cgas-snp
!cgas-snp -i aligned_genes/ -o SNP_Results/
!cgas-snp --no-figures

# Using ! operator with python
!python cgas_module10.py
!python cgas_module10.py -i aligned_genes/ -o SNP_Results/
!python cgas_module10.py --no-figures
```

---

## Input Requirements

### Supported File Formats
FASTA alignment files:
- `.fasta`
- `.fa`
- `.fna`

### Critical Requirements

**1. Exactly 2 Sequences Per File**
```fasta
>Reference_Species
ATGCATGCATGCATGC
>Query_Species
ATGGATGCTTGCATGC
```
- First sequence: Reference
- Second sequence: Query
- ❌ Files with 1, 3, or more sequences will be skipped

**2. Sequences Must Be Pre-Aligned**
- Same length required
- Gaps (-) allowed and will be ignored
- Alignment must be done beforehand (use MAFFT, MUSCLE, etc.)

**3. Standard Nucleotides**
- Valid bases: A, T, G, C
- Gaps: - (ignored in analysis)
- Ambiguous bases (N, R, Y, etc.) are ignored

### File Organization

```bash
alignments/
├── rbcL_aligned.fasta          # Gene 1 alignment
├── matK_aligned.fasta          # Gene 2 alignment
├── ndhF_aligned.fasta          # Gene 3 alignment
└── trnH-psbA_aligned.fasta     # Intergenic spacer
```

### Example Alignment Format

```fasta
>Arabidopsis_thaliana
ATGGCTAGCTAGCTAGCTAGCTAGC-TAGCTAGC
>Oryza_sativa
ATGGCTAGCTAGGTAGCTAGGTAGC-TAGGTAGC
```

**What gets analyzed:**
- Position 13: G→G (identical, not counted)
- Position 14: C→G (substitution: C→G)
- Position 20: C→G (substitution: C→G)
- Position 24: C→G (substitution: C→G)
- Gap at position 26: Ignored

---

## Output Structure

### Directory Organization

```
Module10_SNP_Analysis/                         # Main output folder
│
├── Individual Analysis Files (one per alignment)
│   ├── rbcL_aligned_Substitutions.xlsx        # rbcL gene analysis
│   │   └── Contains:
│   │       - Summary sheet (Ts/Tv ratio, total substitutions)
│   │       - Detailed sheet (all 12 substitution types)
│   │       - Position lists (exact coordinates)
│   │
│   ├── matK_aligned_Substitutions.xlsx        # matK gene analysis
│   └── ndhF_aligned_Substitutions.xlsx        # ndhF gene analysis
│
├── Merged Comparative Analysis
│   ├── Complete_Substitution_Analysis.xlsx    # All samples combined
│   │   └── Contains:
│   │       - Comparison table (all alignments)
│   │       - Ts/Tv ratios for each
│   │       - Total substitution counts
│   │
│   └── Complete_Substitution_Analysis.csv     # R-friendly format
│
└── Figures/ (if R available)
    ├── Substitution_Types.pdf                 # Vector format (scalable)
    ├── Substitution_Types.png                 # Raster format (600 DPI)
    ├── Ts_Tv_Ratio.pdf                        # Vector format
    └── Ts_Tv_Ratio.png                        # Raster format (600 DPI)
```

### Key Output Files Explained

#### 1. Individual Substitution Files (.xlsx)

**Sheet 1: Summary**
```
File Name: rbcL_aligned.fasta
Reference: Arabidopsis_thaliana
Query: Oryza_sativa

Total Positions Compared: 1,428
Identical Positions: 1,380
Substitutions: 48

Transitions (Ts): 32
Transversions (Tv): 16
Ts/Tv Ratio: 2.00
```

**Sheet 2: Substitution Details**
| Substitution Type | Count | Percentage | Positions |
|-------------------|-------|------------|-----------|
| A→G (Transition) | 8 | 16.7% | 45, 123, 456, ... |
| G→A (Transition) | 10 | 20.8% | 78, 234, 567, ... |
| C→T (Transition) | 7 | 14.6% | 90, 345, ... |
| ... | ... | ... | ... |

#### 2. Complete Substitution Analysis (.xlsx)

**Merged comparison table:**
| Sample | A→G | G→A | C→T | T→C | ... | Total | Ts | Tv | Ts/Tv |
|--------|-----|-----|-----|-----|-----|-------|----|----|-------|
| rbcL | 8 | 10 | 7 | 7 | ... | 48 | 32 | 16 | 2.00 |
| matK | 5 | 8 | 6 | 5 | ... | 35 | 24 | 11 | 2.18 |
| ndhF | 12 | 15 | 9 | 10 | ... | 62 | 46 | 16 | 2.88 |

#### 3. Visualizations (Figures/)

**Substitution_Types graph:**
- Grouped bar chart showing all 12 substitution types
- Complementary pairs grouped together (e.g., A→G with G→A)
- Different colors for each comparison
- Publication-ready formatting

**Ts_Tv_Ratio graph:**
- Bar chart comparing Ts/Tv ratios across samples
- Horizontal reference line at theoretical ratio (2.0)
- Clear labeling for publication
---

## Substitution Types

### All 12 Nucleotide Substitution Types

**Transitions (4 types):**
1. **A → G** (purine to purine)
2. **G → A** (purine to purine)
3. **C → T** (pyrimidine to pyrimidine)
4. **T → C** (pyrimidine to pyrimidine)

**Transversions (8 types):**
5. **A → C** (purine to pyrimidine)
6. **C → A** (pyrimidine to purine)
7. **A → T** (purine to pyrimidine)
8. **T → A** (pyrimidine to purine)
9. **G → C** (purine to pyrimidine)
10. **C → G** (pyrimidine to purine)
11. **G → T** (purine to pyrimidine)
12. **T → G** (pyrimidine to purine)

### Chemical Classification

**Purines (larger molecules):**
- Adenine (A)
- Guanine (G)

**Pyrimidines (smaller molecules):**
- Cytosine (C)
- Thymine (T)

**Transition:** Same class (purine↔purine or pyrimidine↔pyrimidine)
**Transversion:** Different class (purine↔pyrimidine)

---

## Transition vs Transversion

### What Are They?

**Transitions (Ts):**
- Substitution within the same chemical class
- A ↔ G (both purines)
- C ↔ T (both pyrimidines)
- **More common** in nature (~2x more frequent)

**Transversions (Tv):**
- Substitution between different classes
- Purine ↔ Pyrimidine
- A/G ↔ C/T
- **Less common** but more diverse (8 types)

### Ts/Tv Ratio

**What it means:**
```
Ts/Tv Ratio = Transition Count / Transversion Count

Example:
Transitions: 32
Transversions: 16
Ts/Tv Ratio = 32/16 = 2.00
```

**Typical values:**
- **~2.0**: Expected for random mutations
- **>2.0**: Excess transitions (common in chloroplasts)
- **<2.0**: Excess transversions (may indicate selection or sequencing errors)

**Biological significance:**
- **High Ts/Tv (~2-3)**: Neutral evolution
- **Low Ts/Tv (<1.5)**: May indicate strong selection or errors
- **Very high (>4)**: Possible CpG hypermutation

### Why Transitions Are More Common

1. **Chemical similarity**: Same ring structure
2. **DNA replication errors**: More likely
3. **Spontaneous deamination**: C→T is common
4. **Lower energy barrier**: Easier transition

---

## R Visualization

### Automatic Figure Generation

If R is installed with required packages, two publication-quality figures are automatically generated.

#### Figure 1: Substitution Types Graph

**Features:**
- **Grouped bar chart**: All 12 substitution types
- **Complementary grouping**: A→G paired with G→A, etc.
- **Color-coded samples**: Different color per alignment
- **Professional styling**: Publication-ready
- **Clear labels**: Substitution types in A→G format

**Purpose:** Compare substitution patterns across genes/samples

#### Figure 2: Ts/Tv Ratio Graph

**Features:**
- **Bar chart**: One bar per sample
- **Reference line**: Horizontal line at 2.0 (theoretical ratio)
- **Clear comparison**: Easy to spot deviations
- **Professional styling**: Publication-ready

**Purpose:** Compare evolutionary rates across genes/samples

### Figure Specifications

**Format:**
- PDF: Vector format (scalable, perfect for publications)
- PNG: Raster format (600 DPI, high resolution)

**Size:**
- Optimized for publication (8" × 6")
- Professional fonts and colors
- Clear axis labels and legends

### Manual R Script Execution

If you want to customize figures:
```bash
# R script is generated in output folder
cd Module10_SNP_Analysis/

# Edit the R script if desired
nano generate_snp_figures.R

# Run manually
Rscript generate_snp_figures.R
```

---

## Troubleshooting

### Common Issues and Solutions

#### 1. openpyxl Not Installed
```bash
❌ ERROR: Required package not installed: No module named 'openpyxl'
```
**Solution:**
```bash
pip install openpyxl

# Verify
python -c "import openpyxl"
```

#### 2. No FASTA Files Found
```bash
❌ ERROR: No FASTA files found!
```
**Solution:**
```bash
# Check file extensions
ls *.fasta
ls *.fa
ls *.fna

# Make sure files are in the correct directory
python cgas_module10.py -i /full/path/to/alignments/
```

#### 3. File Has Wrong Number of Sequences
```bash
⚠️ Skipping file.fasta: Expected 2 sequences, found 3
```
**Solution:**
```
Each FASTA file must contain EXACTLY 2 sequences:
- First sequence: Reference
- Second sequence: Query

If you have multiple sequences, create separate pairwise alignment files.
```

#### 4. Sequences Not Aligned (Different Lengths)
```bash
⚠️ Skipping file.fasta: Sequences have different lengths
```
**Solution:**
```bash
# Align sequences first using MAFFT, MUSCLE, or other alignment tools

# Example with MAFFT:
mafft --auto input.fasta > aligned.fasta

# Then run Module 10:
python cgas_module10.py -i .
```

#### 5. R Not Found (Warning Only)
```bash
⚠️ R not found. Skipping figure generation.
```
**Solution:**
- This is just a warning - Excel files are still created
- Install R if you want figures:
  ```bash
  # Ubuntu/Debian
  sudo apt-get install r-base
  
  # macOS
  brew install r
  
  # Install R packages
  Rscript -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'reshape2'))"
  ```

#### 6. R Packages Missing
```bash
⚠️ R package 'ggplot2' not found
```
**Solution:**
```bash
# Install missing packages
Rscript -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'reshape2'), repos='https://cloud.r-project.org')"
```

#### 7. No Substitutions Found
```bash
⚠️ No substitutions found in alignment.fasta
```
**Solution:**
```
This means sequences are identical.
- Check if you're comparing the same species
- Verify alignment contains different sequences
```

### Quick Fixes Table

| Problem | Symptom | Quick Fix |
|---------|---------|-----------|
| openpyxl missing | Import error | `pip install openpyxl` |
| No FASTA files | "No files found" | Check file extensions (.fasta, .fa) |
| Wrong seq count | "Expected 2, found X" | Ensure exactly 2 sequences per file |
| Not aligned | "Different lengths" | Align sequences with MAFFT first |
| R missing | R warning | Install R or use --no-figures |
| R packages missing | Package error | Install with Rscript |
| No substitutions | "No substitutions" | Sequences are identical |

---

## Examples

### Example 1: Complete Workflow
```bash
# 1. Prepare aligned sequences
mkdir -p /home/abdullah/snp_analysis/
cd /home/abdullah/snp_analysis/

# 2. Align your sequences (if not already aligned)
for gene in rbcL matK ndhF; do
    mafft --auto ${gene}.fasta > ${gene}_aligned.fasta
done

# 3. Run SNP analysis
python cgas_module10.py

# 4. Check results
ls Module10_SNP_Analysis/
open Module10_SNP_Analysis/Complete_Substitution_Analysis.xlsx
open Module10_SNP_Analysis/Figures/Substitution_Types.pdf
```

### Example 2: Process Specific Folder
```bash
# Alignments in separate folder
python cgas_module10.py -i aligned_genes/ -o snp_results/

# Check what was created
ls snp_results/
```

### Example 3: Skip Visualizations (Faster)
```bash
# Generate only Excel files (no R)
python cgas_module10.py --no-figures

# Useful for quick analysis or when R is not available
```

### Example 4: Batch Processing with Script
```bash
#!/bin/bash
# Process multiple datasets

for dataset in monocots dicots ferns; do
    echo "Analyzing $dataset..."
    python cgas_module10.py \
        -i ${dataset}_alignments/ \
        -o ${dataset}_snps/
    echo "Completed $dataset"
done

echo "All datasets analyzed!"
```

### Example 5: From MAFFT to SNP Analysis
```bash
# Complete pipeline from raw sequences to SNP analysis

# Step 1: Align with MAFFT
mafft --auto rbcL.fasta > rbcL_aligned.fasta

# Step 2: Run SNP analysis
python cgas_module10.py

# Step 3: View results
open Module10_SNP_Analysis/rbcL_aligned_Substitutions.xlsx
```

---

## FAQ

### Q1: What's the simplest way to run this?
**A:** Just `python cgas_module10.py` in a folder with aligned FASTA files.

### Q2: Do my sequences need to be aligned?
**A:** Yes! Sequences must be pre-aligned and of equal length. Use MAFFT, MUSCLE, or similar tools first.

### Q3: How many sequences should be in each FASTA file?
**A:** Exactly 2 sequences - one reference and one query. Files with more or fewer will be skipped.

### Q4: What if I don't have R installed?
**A:** No problem! Excel reports are still generated. Use `--no-figures` to suppress R warnings.

### Q5: What does Ts/Tv ratio mean?
**A:** Transition/Transversion ratio. Typical values are ~2.0. Higher values suggest neutral evolution, lower values may indicate selection or errors.

### Q6: Can I analyze protein sequences?
**A:** No, this module is specifically for nucleotide sequences (DNA). It expects A, T, G, C bases.

### Q7: How are gaps handled?
**A:** Gaps (-) are completely ignored. Only positions with A, T, G, or C in both sequences are analyzed.

### Q8: Can I process multiple gene alignments at once?
**A:** Yes! Put all FASTA files in one folder and run the module once. It will process all files and create a merged analysis.

### Q9: What file formats are supported?
**A:** FASTA files with extensions: .fasta, .fa, or .fna

### Q10: How do I cite this tool?
**A:** See Citation section below.

### Q11: Why are some substitutions more common than others?
**A:** Transitions (A↔G, C↔T) are chemically easier and more common than transversions (purine↔pyrimidine).

### Q12: Can I use this for mitochondrial genomes?
**A:** Yes! It works with any aligned DNA sequences, including mitochondrial, nuclear, or viral genomes.

---

## Technical Specifications

### Performance
- **Processing speed**: ~5 seconds per alignment
- **Memory usage**: <100 MB RAM
- **Disk space**: Minimal (<10 MB per analysis)

### Compatibility
- **Python**: 3.9+ (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating Systems**: Windows, macOS, Linux
- **Excel**: Generates .xlsx files (Excel 2010+)
- **R**: 4.0+ (optional, for visualization)

### Input Limits
- **Max alignment length**: No practical limit (tested with 200kb+)
- **Max files**: No practical limit (tested with 100+ files)
- **Sequences per file**: Must be exactly 2

### Quality Features
- ✅ All 12 substitution types tracked
- ✅ Position-by-position analysis
- ✅ Transition/transversion classification
- ✅ Exact position recording
- ✅ Gap handling (ignored)
- ✅ Batch processing
- ✅ Merged comparative analysis
- ✅ Professional visualizations

---

## References
### Python Libraries
- **openpyxl**: [Official Documentation](https://openpyxl.readthedocs.io/)
- Excel file creation and manipulation

### R Packages
- **ggplot2**: [Documentation](https://ggplot2.tidyverse.org/) - Data visualization
- **dplyr**: [Documentation](https://dplyr.tidyverse.org/) - Data manipulation
- **tidyr**: [Documentation](https://tidyr.tidyverse.org/) - Data tidying
- **reshape2**: Data transformation

### Sequence Alignment Tools
- **MAFFT**: [Katoh & Standley (2013)](https://doi.org/10.1093/molbev/mst010)
- **MUSCLE**: Edgar (2004). *Nucleic Acids Research*, 32(5), 1792-1797.

### CGAS Suite
- **CGAS**: Chloroplast Genome Analysis Suite - integrated tools for comprehensive chloroplast genomics

### Citation
If using CGAS Module 10 in publications, please cite:
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
python cgas_module10.py --help

# 2. Verify openpyxl is installed
python -c "import openpyxl"

# 3. Check FASTA files
ls *.fasta *.fa *.fna

# 4. Verify alignment format
head -20 your_alignment.fasta
```

### Common Issues Solved Here
- ✅ openpyxl not installed? Run `pip install openpyxl`
- ✅ No files found? Check file extensions (.fasta, .fa, .fna)
- ✅ Wrong sequence count? Each file needs exactly 2 sequences
- ✅ Not aligned? Use MAFFT or MUSCLE first
- ✅ R missing? Install R or use --no-figures
- ✅ No substitutions? Sequences may be identical

### Quick Reference Card
```bash
# ⚡ QUICK START ⚡
cgas --module 10                                 # cgas command
cgas-snp                                         # shortcut command
python cgas_module10.py                          # python command
python cgas_module10.py -i alignments/           # Specify input
python cgas_module10.py -i data/ -o results/     # Custom output
python cgas_module10.py --no-figures             # Skip R graphs

# 🔬 JUPYTER NOTEBOOK 🔬
%run cgas_module10.py                            # %run works for Module 10
!cgas --module 10                                # ! also works
!cgas-snp --no-figures                           # Skip R in notebook

# 📊 OUTPUT 📊
# Module10_SNP_Analysis/
# ├── [gene]_Substitutions.xlsx       # Individual reports
# ├── Complete_Substitution_Analysis.xlsx  # Merged analysis
# ├── Complete_Substitution_Analysis.csv   # R data
# └── Figures/                        # PDF + PNG graphs

# 🎯 TIPS 🎯
# - Sequences must be PRE-ALIGNED
# - Each FASTA needs EXACTLY 2 sequences
# - Gaps (-) are automatically ignored
# - Typical Ts/Tv ratio: ~2.0
# - Higher Ts/Tv: more transitions (common in chloroplasts)
# - Complementary substitutions grouped in graphs
```

---

**Last Updated**: January 2026  
**Version**: 1.0.1  
**Author**: Abdullah  
**Part of**: CGAS (Chloroplast Genome Analysis Suite)

---

*Happy SNP Analysis! 🧬✨*
