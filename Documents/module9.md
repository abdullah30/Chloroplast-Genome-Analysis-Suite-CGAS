# CGAS Module 9: Comprehensive Amino Acid Analysis with R Visualization
## Complete Documentation and User Guide

---

## 📋 Table of Contents
1. [Introduction](#introduction)
2. [Quick Start - Just Run It!](#quick-start---just-run-it)
3. [Command-Line Usage with Examples](#command-line-usage-with-examples)
4. [Jupyter Notebook Usage](#jupyter-notebook-usage)
5. [Input Requirements](#input-requirements)
6. [Output Structure](#output-structure)
7. [Amino Acid Properties](#amino-acid-properties)
8. [R Visualization](#r-visualization)
9. [Troubleshooting](#troubleshooting)
10. [Examples](#examples)
11. [FAQ](#faq)
12. [Technical Specifications](#technical-specifications)
13. [References](#references)

---

## Introduction

**CGAS Module 9** is a comprehensive tool for analyzing amino acid composition in chloroplast genomes. This module extracts all coding sequences (CDS), translates them to amino acids, calculates composition percentages, and generates publication-quality visualizations for comparative analysis.

This module performs complete amino acid analysis with:
- **CDS extraction**: Identifies all protein-coding sequences from GenBank files
- **Translation**: Converts DNA to amino acid sequences using standard genetic code
- **Composition analysis**: Calculates percentage of all 20 amino acids
- **Individual reports**: Excel files for each genome
- **Merged analysis**: Comparative statistics across all genomes
- **Heatmap visualization**: Publication-quality heatmap showing composition patterns
- **Bar plot comparison**: Direct comparison of representative species

### Key Features:
- **Automatic CDS detection**: Finds all coding sequences in GenBank files
- **20 amino acids**: Complete composition analysis (A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V)
- **Batch processing**: Analyzes multiple genomes at once
- **Excel reports**: Individual and merged spreadsheets
- **Professional visualizations**: Heatmap and bar plots (PDF + PNG, 600 DPI)
- **Full names**: Both single-letter codes and full amino acid names

### Scientific Applications:
- **Comparative Genomics**: Compare amino acid usage across species
- **Protein Evolution**: Study amino acid composition patterns
- **Codon Usage Bias**: Understand translation efficiency preferences
- **Phylogenetic Analysis**: Identify amino acid composition signatures
- **Metabolic Studies**: Analyze biosynthetic costs of proteins
- **Environmental Adaptation**: Detect composition shifts in different habitats

### Part of CGAS (Chloroplast Genome Analysis Suite)
CGAS is a comprehensive suite of tools designed for chloroplast genome analysis, providing seamless integration across different analysis modules.

---

## Quick Start - Just Run It!

### ⚡ One Command Does It All!

```bash
# Option 1: Using cgas command
cgas --module 9

# Option 2: Using cgas-amino shortcut
cgas-amino

# Option 3: Using python directly
python cgas_module9.py
```

**What happens when you run this:**
1. ✅ Finds all GenBank files in current directory
2. ✅ Extracts all CDS features from each genome
3. ✅ Translates DNA sequences to amino acids
4. ✅ Calculates amino acid composition percentages
5. ✅ Creates individual Excel reports per genome
6. ✅ Generates merged comparative analysis
7. ✅ Produces publication-quality visualizations (if R installed)

### 📁 Real-World Folder Examples

#### Example 1: Basic Usage (Easiest!)
```bash
# Your folder structure:
/home/abdullah/amino_acid_analysis/
├── genomes/
│   ├── Arabidopsis_thaliana.gb
│   ├── Oryza_sativa.gb
│   ├── Zea_mays.gb
│   └── Nicotiana_tabacum.gb

# Navigate to your folder
cd /home/abdullah/amino_acid_analysis/genomes/

# Run the module (no arguments needed!)
python cgas_module9.py

# Output created automatically:
# Module9_Amino_Acid_Analysis/
# ├── Arabidopsis_thaliana_AminoAcid.xlsx
# ├── Oryza_sativa_AminoAcid.xlsx
# ├── Zea_mays_AminoAcid.xlsx
# ├── Nicotiana_tabacum_AminoAcid.xlsx
# ├── Complete_Amino_Acid_Analysis.xlsx
# ├── Complete_Amino_Acid_Analysis.csv
# └── Figures/
#     ├── AA_Composition_Heatmap.pdf
#     ├── AA_Composition_Heatmap.png
#     ├── Amino_Acid_Comparison.pdf
#     └── Amino_Acid_Comparison.png
```

#### Example 2: Specify Input Folder
```bash
# Process genomes from specific directory
python cgas_module9.py -i chloroplast_genomes/

# Output created in:
# chloroplast_genomes/Module9_Amino_Acid_Analysis/
```

#### Example 3: Custom Input and Output
```bash
# Full control over directories
python cgas_module9.py -i genbank_files/ -o AA_Results/

# Input from: genbank_files/
# Output to: AA_Results/
```

#### Example 4: Skip Visualizations
```bash
# Generate only Excel files (no R graphs)
python cgas_module9.py --no-figures

# Useful if R is not installed
```

### 🎯 Quick Start Cheat Sheet

```bash
# ====== cgas command ======
cgas --module 9                                  # Process current directory
cgas --module 9 -i genomes/                      # Specify input folder
cgas --module 9 -i data/ -o results/             # Custom input and output
cgas --module 9 --no-figures                     # Skip R visualizations

# ====== cgas-amino shortcut ======
cgas-amino                                       # Process current directory
cgas-amino -i genomes/                           # Specify input folder
cgas-amino -i data/ -o results/                  # Custom input and output
cgas-amino --no-figures                          # Skip R visualizations

# ====== python command ======
python cgas_module9.py                           # Process current directory
python cgas_module9.py -i genomes/               # Specify input folder
python cgas_module9.py -i data/ -o results/      # Custom input and output
python cgas_module9.py --no-figures              # Skip R visualizations
python cgas_module9.py --help                    # Get help
```

### 📊 What You Get (Output Files)

```
Module9_Amino_Acid_Analysis/              # Created automatically
├── 📊 Individual Reports (one per genome)
│   ├── Arabidopsis_thaliana_AminoAcid.xlsx
│   ├── Oryza_sativa_AminoAcid.xlsx
│   └── Zea_mays_AminoAcid.xlsx
│
├── 📈 Merged Analysis
│   ├── Complete_Amino_Acid_Analysis.xlsx  # All genomes combined
│   └── Complete_Amino_Acid_Analysis.csv   # For R visualization
│
└── 🖼️  Figures/ (if R installed)
    ├── AA_Composition_Heatmap.pdf         # Heatmap (vector)
    ├── AA_Composition_Heatmap.png         # Heatmap (600 DPI)
    ├── Amino_Acid_Comparison.pdf          # Bar plot (vector)
    └── Amino_Acid_Comparison.png          # Bar plot (600 DPI)
```
---

## Command-Line Usage with Examples

### Complete Command Reference

```bash
# ====================================================================
# USING cgas COMMAND
# ====================================================================

# Process current directory
cgas --module 9

# Specify input directory
cgas --module 9 -i /home/abdullah/genomes/

# Custom input and output directories
cgas --module 9 -i genbank_files/ -o AA_Results/

# Skip R visualization
cgas --module 9 --no-figures
```

```bash
# ====================================================================
# USING cgas-amino SHORTCUT
# ====================================================================

# Process current directory
cgas-amino

# With specific input and output
cgas-amino -i genbank_files/ -o AA_Results/

# Skip R visualization
cgas-amino --no-figures
```

```bash
# ====================================================================
# USING python COMMAND
# ====================================================================

# ====================================================================
# BASIC USAGE
# ====================================================================

# 1. Process current directory (simplest!)
python cgas_module9.py

# 2. Specify input folder
python cgas_module9.py -i genomes/

# 3. Custom input and output folders
python cgas_module9.py -i genbank_files/ -o AA_Results/

# 4. Skip R visualization
python cgas_module9.py --no-figures

# 5. Get help
python cgas_module9.py --help


# ====================================================================
# REAL-WORLD EXAMPLES
# ====================================================================

# Example 1: Analyze chloroplast genomes
python cgas_module9.py -i chloroplast_genomes/

# Example 2: Save to specific output folder
python cgas_module9.py -i data/ -o ../results/amino_acids/

# Example 3: No R visualization (faster)
python cgas_module9.py -i genomes/ --no-figures

# Example 4: Windows users
python cgas_module9.py -i "C:\Users\abdullah\genomes"
```

### Parameter Details

| Parameter | Short | Default | Description |
|-----------|-------|---------|-------------|
| `--input` | `-i` | `.` (current) | Input directory with GenBank files |
| `--output` | `-o` | `Module9_Amino_Acid_Analysis` | Output directory for results |
| `--no-figures` | - | False | Skip R visualization generation |

---

## Jupyter Notebook Usage

> **Note:** Module 9 is a pure Python script, so both `%run` and the `!` shell operator work. Use `%run` for direct script execution, or `!` for the `cgas` / `cgas-amino` commands. Run the cell from the directory containing your GenBank files, or pass `-i` to specify the input path.

```python
# Using %run (executes the script directly)
%run cgas_module9.py
%run cgas_module9.py -i genomes/
%run cgas_module9.py -i data/ -o results/
%run cgas_module9.py --no-figures

# Using ! operator with cgas command
!cgas --module 9
!cgas --module 9 -i genomes/
!cgas --module 9 --no-figures

# Using ! operator with cgas-amino shortcut
!cgas-amino
!cgas-amino -i genbank_files/ -o AA_Results/
!cgas-amino --no-figures

# Using ! operator with python
!python cgas_module9.py
!python cgas_module9.py -i genbank_files/ -o AA_Results/
!python cgas_module9.py --no-figures
```

---

## Input Requirements

### Supported File Formats
GenBank annotation files:
- `.gb` (GenBank)
- `.gbf` (GenBank Flat)
- `.gbk` (GenBank)
- `.genbank`

### File Requirements

**Must contain:**
- Complete GenBank annotation
- CDS (Coding Sequence) features with translations
- Valid DNA sequences

**Example GenBank structure:**
```
LOCUS       NC_000932             154478 bp    DNA     circular PLN
DEFINITION  Arabidopsis thaliana chloroplast, complete genome.
...
FEATURES             Location/Qualifiers
     source          1..154478
                     /organism="Arabidopsis thaliana"
     gene            1..1485
                     /gene="psbA"
     CDS             1..1485
                     /gene="psbA"
                     /product="photosystem II protein D1"
                     /translation="MTAILER..."
```

### File Organization

```bash
genomes/
├── Arabidopsis_thaliana.gb       # Annotated genome 1
├── Oryza_sativa.gb               # Annotated genome 2
├── Zea_mays.gb                   # Annotated genome 3
└── Nicotiana_tabacum.gb          # Annotated genome 4
```

### What Gets Analyzed

**CDS features:**
- All features with `type="CDS"`
- Extracted and translated to amino acids
- Translation uses standard genetic code
- Stops at first stop codon

**Excluded:**
- tRNA, rRNA genes (not translated)
- Introns (removed during translation)
- Intergenic regions
- Hypothetical proteins (included if annotated as CDS)

---

## Output Structure

### Directory Organization

```
Module9_Amino_Acid_Analysis/                    # Main output folder
│
├── Individual Analysis Files (one per genome)
│   ├── Arabidopsis_thaliana_AminoAcid.xlsx    # Arabidopsis analysis
│   │   └── Contains:
│   │       - Amino acid composition (20 AA)
│   │       - Percentages for each amino acid
│   │       - Full names and single-letter codes
│   │
│   ├── Oryza_sativa_AminoAcid.xlsx            # Rice analysis
│   └── Zea_mays_AminoAcid.xlsx                # Maize analysis
│
├── Merged Comparative Analysis
│   ├── Complete_Amino_Acid_Analysis.xlsx      # All genomes combined
│   │   └── Contains:
│   │       - Comparison table (all genomes)
│   │       - Percentages for each AA across genomes
│   │       - Easy comparison of composition
│   │
│   └── Complete_Amino_Acid_Analysis.csv       # R-friendly format
│
└── Figures/ (if R available)
    ├── AA_Composition_Heatmap.pdf             # Vector format (scalable)
    ├── AA_Composition_Heatmap.png             # Raster format (600 DPI)
    ├── Amino_Acid_Comparison.pdf              # Vector format
    └── Amino_Acid_Comparison.png              # Raster format (600 DPI)
```

### Key Output Files Explained

#### 1. Individual Amino Acid Files (.xlsx)

**Table structure:**
| AA_Code | Amino_Acid | Percentage |
|---------|------------|------------|
| A | Alanine | 8.45% |
| R | Arginine | 5.23% |
| N | Asparagine | 4.67% |
| D | Aspartate | 5.12% |
| C | Cysteine | 1.23% |
| Q | Glutamine | 3.89% |
| E | Glutamate | 6.78% |
| G | Glycine | 7.34% |
| H | Histidine | 2.11% |
| I | Isoleucine | 5.67% |
| L | Leucine | 9.23% |
| K | Lysine | 6.45% |
| M | Methionine | 2.34% |
| F | Phenylalanine | 4.56% |
| P | Proline | 4.78% |
| S | Serine | 6.89% |
| T | Threonine | 5.34% |
| W | Tryptophan | 1.45% |
| Y | Tyrosine | 3.67% |
| V | Valine | 6.12% |

**Total:** 100.00% (all 20 standard amino acids)

#### 2. Complete Amino Acid Analysis (.xlsx)

**Merged comparison table:**
| AA | Amino_Acid | Arabidopsis | Oryza | Zea | Nicotiana | Average |
|----|------------|-------------|-------|-----|-----------|---------|
| A | Alanine | 8.45% | 8.62% | 8.33% | 8.51% | 8.48% |
| R | Arginine | 5.23% | 5.31% | 5.18% | 5.27% | 5.25% |
| ... | ... | ... | ... | ... | ... | ... |

**Uses:**
- Compare composition across species
- Identify conserved vs. variable amino acids
- Detect compositional biases
- Publication-ready table

#### 3. Visualizations (Figures/)

**Heatmap (AA_Composition_Heatmap):**
- Rows: All 20 amino acids
- Columns: All analyzed genomes
- Color intensity: Percentage (higher = darker)
- Purpose: Visual comparison of composition patterns

**Bar Plot (Amino_Acid_Comparison):**
- 3 representative species side-by-side
- All 20 amino acids shown
- Different colors per species
- Purpose: Direct quantitative comparison
---

## Amino Acid Properties

### 20 Standard Amino Acids

#### Nonpolar (Hydrophobic)
1. **A - Alanine**: Small, hydrophobic
2. **V - Valine**: Branched-chain, hydrophobic
3. **L - Leucine**: Branched-chain, hydrophobic (most common)
4. **I - Isoleucine**: Branched-chain, hydrophobic
5. **M - Methionine**: Contains sulfur, hydrophobic
6. **F - Phenylalanine**: Aromatic, hydrophobic
7. **W - Tryptophan**: Aromatic, hydrophobic (rarest)
8. **P - Proline**: Unique cyclic structure

#### Polar (Hydrophilic)
9. **G - Glycine**: Smallest, flexible (no side chain)
10. **S - Serine**: Small polar
11. **T - Threonine**: Polar, hydroxyl group
12. **C - Cysteine**: Contains sulfur, forms disulfide bonds
13. **Y - Tyrosine**: Aromatic, polar
14. **N - Asparagine**: Amide group
15. **Q - Glutamine**: Amide group

#### Charged (Ionizable)
16. **D - Aspartate**: Negatively charged (acidic)
17. **E - Glutamate**: Negatively charged (acidic)
18. **K - Lysine**: Positively charged (basic)
19. **R - Arginine**: Positively charged (basic)
20. **H - Histidine**: Positively charged (basic, pKa ~6)

### Typical Chloroplast Protein Composition

**Most abundant:**
- **Leucine (L)**: ~9-10% (most common)
- **Alanine (A)**: ~8-9%
- **Glycine (G)**: ~7-8%
- **Serine (S)**: ~6-7%

**Least abundant:**
- **Tryptophan (W)**: ~1-2% (rarest)
- **Cysteine (C)**: ~1-2%
- **Methionine (M)**: ~2-3%

**Biological significance:**
- **High leucine**: Important for protein structure
- **High glycine**: Flexibility in proteins
- **Low cysteine**: Fewer disulfide bonds in chloroplasts
- **Low tryptophan**: Biosynthetically expensive

---

## R Visualization

### Automatic Figure Generation

If R is installed with required packages, two publication-quality figures are automatically generated.

#### Figure 1: Amino Acid Composition Heatmap

**Features:**
- **Heatmap layout**: Rows = amino acids, Columns = genomes
- **Color gradient**: Low (light) to High (dark) percentage
- **Clustered**: Related genomes and amino acids grouped
- **Professional styling**: Publication-ready
- **Color scheme**: Blue gradient (pheatmap default)

**Purpose:** 
- Identify compositional patterns across genomes
- Detect conserved amino acid usage
- Visualize clustering of related species

#### Figure 2: Amino Acid Comparison Bar Plot

**Features:**
- **Side-by-side bars**: Compare 3 representative species
- **All 20 amino acids**: Complete profile
- **Color-coded**: Different color per species
- **Clear labels**: Full amino acid names on x-axis
- **Professional styling**: Publication-ready

**Purpose:** 
- Direct quantitative comparison
- Highlight differences between species
- Easy interpretation for readers

### Figure Specifications

**Format:**
- PDF: Vector format (scalable, perfect for publications)
- PNG: Raster format (600 DPI, high resolution)

**Size:**
- Heatmap: Optimized for number of genomes
- Bar plot: Standard 8" × 6"
- Professional fonts and colors

### Manual R Script Execution

If you want to customize figures:
```bash
# R script is generated in output folder
cd Module9_Amino_Acid_Analysis/

# Edit the R script if desired
nano generate_aa_figures.R

# Run manually
Rscript generate_aa_figures.R
```

---

## Troubleshooting

### Common Issues and Solutions

#### 1. BioPython Not Installed
```bash
❌ ERROR: Required package not installed: No module named 'Bio'
```
**Solution:**
```bash
pip install biopython

# Or with --break-system-packages
pip install biopython --break-system-packages

# Verify
python -c "from Bio import SeqIO"
```

#### 2. pandas or openpyxl Not Installed
```bash
❌ ERROR: No module named 'pandas' or 'openpyxl'
```
**Solution:**
```bash
pip install pandas openpyxl

# Verify
python -c "import pandas; import openpyxl"
```

#### 3. No GenBank Files Found
```bash
❌ ERROR: No GenBank files found!
```
**Solution:**
```bash
# Check file extensions
ls *.gb
ls *.gbk
ls *.gbf

# Make sure files are in the correct directory
python cgas_module9.py -i /full/path/to/genomes/
```

#### 4. No CDS Features Found
```bash
⚠️ Warning: No CDS features found in genome.gb
```
**Solution:**
```
GenBank file lacks CDS annotations.
- Check if file is properly annotated
- Make sure it's not just raw sequence
- Use annotated genomes from NCBI
- Or annotate with Module 2 first
```

#### 5. Translation Errors
```bash
⚠️ Warning: Error translating CDS: Sequence length not a multiple of three
```
**Solution:**
```
CDS sequence is incomplete or corrupted.
- This is just a warning - other CDS still processed
- Check GenBank file quality
- May indicate annotation issues
- Usually safe to ignore if only a few CDS affected
```

#### 6. R Not Found (Warning Only)
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
  Rscript -e "install.packages(c('ggplot2', 'pheatmap', 'RColorBrewer', 'reshape2', 'dplyr', 'tidyr'))"
  ```

#### 7. R Packages Missing
```bash
⚠️ R package 'pheatmap' not found
```
**Solution:**
```bash
# Install missing packages
Rscript -e "install.packages(c('ggplot2', 'pheatmap', 'RColorBrewer', 'reshape2', 'dplyr', 'tidyr'), repos='https://cloud.r-project.org')"
```

### Quick Fixes Table

| Problem | Symptom | Quick Fix |
|---------|---------|-----------|
| BioPython missing | Import error | `pip install biopython` |
| pandas missing | Import error | `pip install pandas openpyxl` |
| No GenBank files | "No files found" | Check file extensions (.gb, .gbk) |
| No CDS features | Warning message | Use annotated GenBank files |
| Translation error | Warning message | Usually safe to ignore |
| R missing | R warning | Install R or use --no-figures |
| R packages missing | Package error | Install with Rscript |

---

## Examples

### Example 1: Complete Workflow
```bash
# 1. Prepare annotated genomes
mkdir -p /home/abdullah/aa_analysis/
cd /home/abdullah/aa_analysis/

# 2. Copy GenBank files (from NCBI or Module 2)
cp /path/to/genomes/*.gb genomes/

# 3. Run amino acid analysis
python cgas_module9.py -i genomes/

# 4. Check results
ls Module9_Amino_Acid_Analysis/
open Module9_Amino_Acid_Analysis/Complete_Amino_Acid_Analysis.xlsx
open Module9_Amino_Acid_Analysis/Figures/AA_Composition_Heatmap.pdf
```

### Example 2: Process Specific Folder
```bash
# Genomes in separate folder
python cgas_module9.py -i chloroplast_genomes/ -o aa_results/

# Check what was created
ls aa_results/
```

### Example 3: Skip Visualizations (Faster)
```bash
# Generate only Excel files (no R)
python cgas_module9.py --no-figures

# Useful for quick analysis or when R is not available
```

### Example 4: Batch Processing with Script
```bash
#!/bin/bash
# Process multiple datasets

for group in monocots dicots ferns; do
    echo "Analyzing $group..."
    python cgas_module9.py \
        -i ${group}_genomes/ \
        -o ${group}_aa_analysis/
    echo "Completed $group"
done

echo "All groups analyzed!"
```

### Example 5: From NCBI to Analysis
```bash
# Complete pipeline

# Step 1: Download GenBank files from NCBI
# (manually or using scripts)

# Step 2: Run amino acid analysis
python cgas_module9.py -i genbank_files/

# Step 3: View results
open Module9_Amino_Acid_Analysis/Complete_Amino_Acid_Analysis.xlsx
```

---

## FAQ

### Q1: What's the simplest way to run this?
**A:** Just `python cgas_module9.py` in a folder with GenBank files.

### Q2: Do I need annotated genomes?
**A:** Yes! GenBank files must have CDS features annotated. Use NCBI genomes or annotate with Module 2.

### Q3: What if I don't have R installed?
**A:** No problem! Excel reports are still generated. Use `--no-figures` to suppress R warnings.

### Q4: How many amino acids are analyzed?
**A:** All 20 standard amino acids (A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V).

### Q5: Can I analyze nuclear genes?
**A:** Yes! This module works with any CDS-annotated GenBank file, including nuclear or mitochondrial genomes.

### Q6: Why are some amino acids more common?
**A:** Leucine, alanine, and glycine are typically most abundant due to codon usage bias and protein structural requirements.

### Q7: Can I process hundreds of genomes?
**A:** Yes! The module handles batch processing. However, the heatmap may become crowded with >50 genomes.

### Q8: What if translation fails for some CDS?
**A:** Warnings are shown but analysis continues. Only affects that specific CDS; others are still processed.

### Q9: How do I cite this tool?
**A:** See Citation section below.

### Q10: Can I compare chloroplast vs. nuclear amino acid composition?
**A:** Yes! Just run the module separately on chloroplast and nuclear GenBank files, then compare the Excel outputs.

### Q11: Why is my heatmap very large?
**A:** If you have many genomes (>30), the heatmap automatically adjusts size. You can manually edit the R script to customize dimensions.

---

## Technical Specifications

### Performance
- **Processing speed**: ~5 seconds per genome
- **Memory usage**: <200 MB RAM
- **Disk space**: Minimal (<10 MB per analysis)

### Compatibility
- **Python**: 3.9+ (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating Systems**: Windows, macOS, Linux
- **Excel**: Generates .xlsx files (Excel 2010+)
- **R**: 4.0+ (optional, for visualization)

### Input Limits
- **Max genomes**: No practical limit (tested with 100+ genomes)
- **Max CDS per genome**: No practical limit (tested with 200+ CDS)
- **Genome size**: Works with any size (tested with 200kb+ genomes)

### Quality Features
- ✅ All 20 amino acids tracked
- ✅ Automatic CDS extraction
- ✅ Standard genetic code translation
- ✅ Stop codon handling
- ✅ Batch processing
- ✅ Merged comparative analysis
- ✅ Professional visualizations
- ✅ Full amino acid names

---

## References
### Python Libraries
- **BioPython**: [Official Documentation](https://biopython.org/)
  - Cock, P. J., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423.
- **pandas**: [User Guide](https://pandas.pydata.org/docs/)
- **openpyxl**: [Documentation](https://openpyxl.readthedocs.io/)

### R Packages
- **ggplot2**: [Documentation](https://ggplot2.tidyverse.org/)
- **pheatmap**: [CRAN](https://cran.r-project.org/package=pheatmap) - Pretty heatmaps
- **RColorBrewer**: Color palettes for visualization

### CGAS Suite
- **CGAS**: Chloroplast Genome Analysis Suite - integrated tools for comprehensive chloroplast genomics

### Citation
If using CGAS Module 9 in publications, please cite:
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
python cgas_module9.py --help

# 2. Verify all packages are installed
python -c "from Bio import SeqIO; import pandas; import openpyxl"

# 3. Check GenBank files
ls *.gb *.gbk *.gbf

# 4. Verify GenBank annotation
grep "CDS" your_genome.gb | head
```

### Common Issues Solved Here
- ✅ BioPython not installed? Run `pip install biopython`
- ✅ No files found? Check file extensions (.gb, .gbk, .gbf)
- ✅ No CDS features? Use annotated GenBank files
- ✅ Translation errors? Usually safe to ignore warnings
- ✅ R missing? Install R or use --no-figures
- ✅ Want custom colors? Edit the generated R script

### Quick Reference Card
```bash
# ⚡ QUICK START ⚡
cgas --module 9                                  # cgas command
cgas-amino                                       # shortcut command
python cgas_module9.py                           # python command
python cgas_module9.py -i genomes/               # Specify input
python cgas_module9.py -i data/ -o results/      # Custom output
python cgas_module9.py --no-figures              # Skip R graphs

# 🔬 JUPYTER NOTEBOOK 🔬
%run cgas_module9.py                             # %run works for Module 9
!cgas --module 9                                 # ! also works
!cgas-amino --no-figures                         # Skip R in notebook

# 📊 OUTPUT 📊
# Module9_Amino_Acid_Analysis/
# ├── [genome]_AminoAcid.xlsx           # Individual reports
# ├── Complete_Amino_Acid_Analysis.xlsx # Merged analysis
# ├── Complete_Amino_Acid_Analysis.csv  # R data
# └── Figures/                          # Heatmap + bar plot

# 🎯 TIPS 🎯
# - Requires annotated GenBank files with CDS features
# - Analyzes all 20 standard amino acids
# - Leucine typically most common (~9-10%)
# - Tryptophan typically rarest (~1-2%)
# - Heatmap clusters related genomes/amino acids
# - Bar plot shows 3 representative species
# - Stop codons excluded from analysis
```

---

**Last Updated**: January 2026  
**Version**: 1.0.1  
**Author**: Abdullah  
**Part of**: CGAS (Chloroplast Genome Analysis Suite)

---

*Happy Amino Acid Analysis! 🧬✨*
