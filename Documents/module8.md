

# CGAS Module 8: Codon Usage Analysis (RSCU) with R Visualization
## Complete Documentation and User Guide

---

## 📋 Table of Contents
1. [Introduction](#introduction)
2. [Quick Start - Just Run It!](#quick-start---just-run-it)
3. [Command-Line Usage with Examples](#command-line-usage-with-examples)
4. [Jupyter Notebook Usage](#jupyter-notebook-usage)
5. [Input Requirements](#input-requirements)
6. [Output Structure](#output-structure)
7. [RSCU Calculation Method](#rscu-calculation-method)
8. [Amino Acid Grouping](#amino-acid-grouping)
9. [R Visualization](#r-visualization)
10. [Troubleshooting](#troubleshooting)
11. [Examples](#examples)
12. [FAQ](#faq)
13. [Technical Specifications](#technical-specifications)
14. [References](#references)

---

## Introduction

**CGAS Module 8** is a comprehensive tool for analyzing codon usage patterns in chloroplast genomes by calculating Relative Synonymous Codon Usage (RSCU) values. This module extracts coding sequences from GenBank files, calculates codon usage statistics, and generates publication-quality visualizations.

This module performs complete codon usage analysis with:
- **CDS extraction**: Automatically extracts coding sequences from GenBank files
- **RSCU calculation**: Computes Relative Synonymous Codon Usage values for all codons
- **Codon counting**: Tallies occurrences of all 61 sense codons
- **Comparative analysis**: Merges data across multiple genomes
- **Excel reports**: Individual and merged analysis spreadsheets
- **R visualizations**: Publication-quality heatmaps and comparison plots (PDF + PNG, 600 DPI)

### Key Features:
- **Automatic CDS extraction**: Finds all coding sequences in GenBank files
- **Comprehensive analysis**: All 61 sense codons across 20 amino acids
- **RSCU values**: Standard metric for codon bias analysis
- **Batch processing**: Analyzes multiple genomes at once
- **Merged analysis**: Comparative statistics across all samples
- **Professional visualizations**: Heatmaps and stacked bar charts

### Scientific Applications:
- **Molecular Evolution**: Study codon usage bias patterns
- **Phylogenetic Analysis**: Compare codon preferences across species
- **Gene Expression**: Infer expression levels from codon usage
- **Genome Evolution**: Analyze mutation and selection pressures
- **Comparative Genomics**: Examine codon usage divergence
- **Biotechnology**: Optimize gene sequences for expression

### Part of CGAS (Chloroplast Genome Analysis Suite)
CGAS is a comprehensive suite of tools designed for chloroplast genome analysis, providing seamless integration across different analysis modules.

---

## Quick Start - Just Run It!

### ⚡ One Command Does It All!

```bash
# Option 1: Using cgas command
cgas --module 8

# Option 2: Using cgas-codon shortcut
cgas-codon

# Option 3: Using python directly
python cgas_module8.py
```

**What happens when you run this:**
1. ✅ Finds all GenBank files in current directory
2. ✅ Extracts coding sequences from each file
3. ✅ Calculates RSCU values for all codons
4. ✅ Creates individual Excel reports per genome
5. ✅ Generates merged comparative analysis
6. ✅ Produces CSV file for R visualization
7. ✅ Creates publication-quality graphs (if R installed)

### 📁 Real-World Folder Examples

#### Example 1: Basic Usage (Easiest!)
```bash
# Your folder structure:
/home/abdullah/codon_analysis/
├── chloroplast_genomes/
│   ├── Arabidopsis_thaliana.gb
│   ├── Oryza_sativa.gb
│   └── Zea_mays.gb

# Navigate to your folder
cd /home/abdullah/codon_analysis/chloroplast_genomes/

# Run the module (no arguments needed!)
python cgas_module8.py

# Output created automatically:
# Module8_Codon_Usage_Analysis/
# ├── Arabidopsis_thaliana_RSCU.xlsx
# ├── Oryza_sativa_RSCU.xlsx
# ├── Zea_mays_RSCU.xlsx
# ├── Complete_Codon_Usage_Analysis.xlsx
# ├── Complete_Codon_Usage_Analysis.csv
# └── Figures/
#     ├── RSCU_Heatmap.pdf
#     ├── RSCU_Heatmap.png
#     ├── Codon_Usage_Comparison.pdf
#     └── Codon_Usage_Comparison.png
```

#### Example 2: Specify Input Folder
```bash
# Process genomes from specific directory
python cgas_module8.py -i genbank_files/

# Output created in:
# genbank_files/Module8_Codon_Usage_Analysis/
```

#### Example 3: Custom Input and Output
```bash
# Full control over directories
python cgas_module8.py -i chloroplast_genomes/ -o Codon_Results/

# Input from: chloroplast_genomes/
# Output to: Codon_Results/
```

#### Example 4: Skip Visualizations
```bash
# Generate only Excel files (no R graphs)
python cgas_module8.py --no-figures

# Useful if R is not installed
```

### 🎯 Quick Start Cheat Sheet

```bash
# ====== cgas command ======
cgas --module 8                                  # Process current directory
cgas --module 8 -i genbank_files/                # Specify input folder
cgas --module 8 -i data/ -o results/             # Custom input and output
cgas --module 8 --no-figures                     # Skip R visualizations

# ====== cgas-codon shortcut ======
cgas-codon                                       # Process current directory
cgas-codon -i genbank_files/                     # Specify input folder
cgas-codon -i data/ -o results/                  # Custom input and output
cgas-codon --no-figures                          # Skip R visualizations

# ====== python command ======
python cgas_module8.py                           # Process current directory
python cgas_module8.py -i genbank_files/         # Specify input folder
python cgas_module8.py -i data/ -o results/      # Custom input and output
python cgas_module8.py --no-figures              # Skip R visualizations
python cgas_module8.py --help                    # Get help
```

### 📊 What You Get (Output Files)

```
Module8_Codon_Usage_Analysis/               # Created automatically
├── 📊 Individual Reports (one per genome)
│   ├── Arabidopsis_thaliana_RSCU.xlsx      # Complete codon usage analysis
│   ├── Oryza_sativa_RSCU.xlsx
│   └── Zea_mays_RSCU.xlsx
│
├── 📈 Merged Analysis
│   ├── Complete_Codon_Usage_Analysis.xlsx  # All genomes combined
│   └── Complete_Codon_Usage_Analysis.csv   # For R visualization
│
└── 🖼️  Figures/ (if R installed)
    ├── RSCU_Heatmap.pdf                    # Publication-quality (vector)
    ├── RSCU_Heatmap.png                    # High-resolution (600 DPI)
    ├── Codon_Usage_Comparison.pdf          # Representative species
    └── Codon_Usage_Comparison.png          # High-resolution (600 DPI)
```
---

## Command-Line Usage with Examples

### Complete Command Reference

```bash
# ====================================================================
# USING cgas COMMAND
# ====================================================================

# Process current directory
cgas --module 8

# Specify input directory
cgas --module 8 -i /home/abdullah/chloroplast_genomes/

# Custom input and output directories
cgas --module 8 -i chloroplast_genomes/ -o Codon_Results/

# Skip R visualization
cgas --module 8 --no-figures
```

```bash
# ====================================================================
# USING cgas-codon SHORTCUT
# ====================================================================

# Process current directory
cgas-codon

# With specific input and output
cgas-codon -i chloroplast_genomes/ -o Codon_Results/

# Skip R visualization
cgas-codon --no-figures
```

```bash
# ====================================================================
# USING python COMMAND
# ====================================================================

# ====================================================================
# BASIC USAGE
# ====================================================================

# 1. Process current directory (simplest!)
python cgas_module8.py

# 2. Specify input folder
python cgas_module8.py -i genbank_files/

# 3. Custom input and output folders
python cgas_module8.py -i chloroplast_genomes/ -o Codon_Results/

# 4. Skip R visualization
python cgas_module8.py --no-figures

# 5. Get help
python cgas_module8.py --help


# ====================================================================
# REAL-WORLD EXAMPLES
# ====================================================================

# Example 1: Process chloroplast genomes
python cgas_module8.py -i chloroplast_genomes/

# Example 2: Save to specific output folder
python cgas_module8.py -i data/ -o ../results/codon_usage/

# Example 3: No R visualization (faster)
python cgas_module8.py -i genbank_files/ --no-figures

# Example 4: Windows users
python cgas_module8.py -i "C:\Users\abdullah\genbank_files"
```

### Parameter Details

| Parameter | Short | Default | Description |
|-----------|-------|---------|-------------|
| `--input` | `-i` | `.` (current) | Input directory with GenBank files |
| `--output` | `-o` | `Module8_Codon_Usage_Analysis` | Output directory for results |
| `--no-figures` | - | False | Skip R visualization generation |

---

## Jupyter Notebook Usage

> **Note:** Module 8 is a pure Python script, so both `%run` and the `!` shell operator work. Use `%run` for direct script execution, or `!` for the `cgas` / `cgas-codon` commands. Run the cell from the directory containing your GenBank files, or pass `-i` to specify the input path.

```python
# Using %run (executes the script directly)
%run cgas_module8.py
%run cgas_module8.py -i chloroplast_genomes/
%run cgas_module8.py -i data/ -o results/
%run cgas_module8.py --no-figures

# Using ! operator with cgas command
!cgas --module 8
!cgas --module 8 -i chloroplast_genomes/
!cgas --module 8 --no-figures

# Using ! operator with cgas-codon shortcut
!cgas-codon
!cgas-codon -i chloroplast_genomes/ -o Codon_Results/
!cgas-codon --no-figures

# Using ! operator with python
!python cgas_module8.py
!python cgas_module8.py -i chloroplast_genomes/ -o Codon_Results/
!python cgas_module8.py --no-figures
```

---

## Input Requirements

### Supported File Formats
GenBank files:
- `.gb`
- `.gbf`
- `.gbk`
- `.genbank`

### Critical Requirements

**1. Complete GenBank Files**
```gb
LOCUS       NC_000932               154478 bp DNA     circular PLN 15-JUL-2022
DEFINITION  Arabidopsis thaliana chloroplast, complete genome.
ACCESSION   NC_000932
VERSION     NC_000932.1
...
FEATURES             Location/Qualifiers
     source          1..154478
                     /organism="Arabidopsis thaliana"
                     /organelle="plastid"
                     /mol_type="genomic DNA"
     CDS             864..2690
                     /gene="rbcL"
                     /product="ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"
                     /translation="MASSSNSS..."
ORIGIN
        1 atggcgacga cgttcgtcgt cgtttgtcga tctcgtctga cttcagcctg atcggtagca
       61 atggtaagtt ggtggtgtga aagcagctga cgggagcatt cggatgtaga tttggagaaa
...
//
```

**2. CDS Features Required**
- Files must contain CDS (coding sequence) features
- Each CDS should have proper start/end coordinates
- Translation annotation is optional but helpful

**3. Standard Nucleotides**
- Valid bases: A, T, G, C
- Ambiguous bases (N, R, Y, etc.) are counted but may generate warnings

### File Organization

```bash
genbank_files/
├── Arabidopsis_thaliana.gb          # Chloroplast genome 1
├── Oryza_sativa.gb                  # Chloroplast genome 2
├── Zea_mays.gb                      # Chloroplast genome 3
└── Nicotiana_tabacum.gb             # Chloroplast genome 4
```

### Example GenBank Format

```gb
FEATURES             Location/Qualifiers
     source          1..154478
                     /organism="Arabidopsis thaliana"
                     /organelle="plastid"
                     /mol_type="genomic DNA"
     CDS             864..2690
                     /gene="rbcL"
                     /product="ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"
     CDS             2770..3624
                     /gene="atpB"
                     /product="ATP synthase subunit beta"
```

**What gets analyzed:**
- All CDS features are extracted
- Coding sequences are concatenated
- Codons are counted across all genes
- RSCU values are calculated for each codon

---

## Output Structure

### Directory Organization

```
Module8_Codon_Usage_Analysis/                # Main output folder
│
├── Individual Analysis Files (one per genome)
│   ├── Arabidopsis_thaliana_RSCU.xlsx       # Arabidopsis genome analysis
│   │   └── Contains:
│   │       - RSCU values for all 61 codons
│   │       - Codon counts for each amino acid
│   │       - Amino acid groupings
│   │
│   ├── Oryza_sativa_RSCU.xlsx               # Rice genome analysis
│   └── Zea_mays_RSCU.xlsx                   # Maize genome analysis
│
├── Merged Comparative Analysis
│   ├── Complete_Codon_Usage_Analysis.xlsx   # All samples combined
│   │   └── Contains:
│   │       - RSCU comparison table (all genomes)
│   │       - Amino acid groupings
│   │       - Codon-level comparisons
│   │
│   └── Complete_Codon_Usage_Analysis.csv    # R-friendly format
│
└── Figures/ (if R available)
    ├── RSCU_Heatmap.pdf                     # Vector format (scalable)
    ├── RSCU_Heatmap.png                     # Raster format (600 DPI)
    ├── Codon_Usage_Comparison.pdf           # Representative species
    └── Codon_Usage_Comparison.png           # Raster format (600 DPI)
```

### Key Output Files Explained

#### 1. Individual RSCU Files (.xlsx)

**RSCU Analysis sheet:**
```
AA    Amino_Acid    Codon    Count    RSCU
Ala   Alanine       GCT      45       0.8235
                    GCC      78       1.4265
                    GCA      23       0.4206
                    GCG      67       1.2294
Arg   Arginine      CGT      12       0.5714
                    CGC      34       1.6190
                    ...
```

**What the columns mean:**
- **AA**: 3-letter amino acid code
- **Amino_Acid**: Full amino acid name
- **Codon**: Specific codon sequence
- **Count**: Number of occurrences in all CDS
- **RSCU**: Relative Synonymous Codon Usage value

#### 2. Complete Codon Usage Analysis (.xlsx)

**Merged comparison table:**
```
Amino_Acid    Codon    Arabidopsis_thaliana    Oryza_sativa    Zea_mays
Alanine       GCT      0.8235                  0.9123          0.7456
              GCC      1.4265                  1.3456          1.5678
              GCA      0.4206                  0.3891          0.4567
              GCG      1.2294                  1.3530          1.2299
Arginine      CGT      0.5714                  0.6234          0.5123
              ...
```

#### 3. Visualizations (Figures/)

**RSCU_Heatmap graph:**
- Color-coded heatmap of RSCU values
- Rows: Codons grouped by amino acid
- Columns: Different genomes/species
- Color scale: Blue (low RSCU) to Red (high RSCU)
- Publication-ready formatting

**Codon_Usage_Comparison graph:**
- Stacked bar charts for 3 representative species
- Each amino acid shown as separate facet
- Codon usage proportions visualized
- Clear labeling for publication

---
## RSCU Calculation Method

### What is RSCU?

**Relative Synonymous Codon Usage (RSCU)** is a widely used metric in codon bias analysis:

```
RSCU = (Observed frequency of codon) / (Expected frequency if all synonymous codons used equally)
```

**Interpretation:**
- **RSCU = 1.0**: No bias (codon used exactly as expected)
- **RSCU > 1.0**: Positive bias (codon preferred)
- **RSCU < 1.0**: Negative bias (codon avoided)

### Example Calculation

For Alanine (4 synonymous codons: GCT, GCC, GCA, GCG):

```
Observed counts:
GCT: 45 codons
GCC: 78 codons
GCA: 23 codons
GCG: 67 codons
Total: 213 codons

Expected frequency if equal usage: 1/4 = 0.25 for each codon

RSCU calculations:
GCT: (45/213) / 0.25 = 0.8235
GCC: (78/213) / 0.25 = 1.4265
GCA: (23/213) / 0.25 = 0.4206
GCG: (67/213) / 0.25 = 1.2294
```

### Biological Significance

**Why RSCU matters:**
- **Gene expression**: Codon bias correlates with expression levels
- **Evolutionary patterns**: Different species show distinct preferences
- **Genome composition**: Reflects mutational pressures and selection
- **Biotechnology**: Optimizing gene sequences for heterologous expression

**Typical patterns in chloroplast genomes:**
- **Highly expressed genes**: Strong codon bias
- **AT-rich genomes**: Preference for A/T-ending codons
- **Conserved patterns**: Similar bias across related species

---

## Amino Acid Grouping

### Standard Genetic Code

**Synonymous Codon Groups:**

| Amino Acid | 3-Letter Code | Synonymous Codons | Number |
|------------|---------------|-------------------|--------|
| Alanine    | Ala           | GCT, GCC, GCA, GCG | 4 |
| Arginine   | Arg           | CGT, CGC, CGA, CGG, AGA, AGG | 6 |
| Asparagine | Asn           | AAT, AAC | 2 |
| Aspartate  | Asp           | GAT, GAC | 2 |
| Cysteine   | Cys           | TGT, TGC | 2 |
| Glutamine  | Gln           | CAA, CAG | 2 |
| Glutamate  | Glu           | GAA, GAG | 2 |
| Glycine    | Gly           | GGT, GGC, GGA, GGG | 4 |
| Histidine  | His           | CAT, CAC | 2 |
| Isoleucine | Ile           | ATT, ATC, ATA | 3 |
| Leucine    | Leu           | TTA, TTG, CTT, CTC, CTA, CTG | 6 |
| Lysine     | Lys           | AAA, AAG | 2 |
| Methionine | Met           | ATG | 1 |
| Phenylalanine | Phe        | TTT, TTC | 2 |
| Proline    | Pro           | CCT, CCC, CCA, CCG | 4 |
| Serine     | Ser           | TCT, TCC, TCA, TCG, AGT, AGC | 6 |
| Threonine  | Thr           | ACT, ACC, ACA, ACG | 4 |
| Tryptophan | Trp           | TGG | 1 |
| Tyrosine   | Tyr           | TAT, TAC | 2 |
| Valine     | Val           | GTT, GTC, GTA, GTG | 4 |
| Stop       | -             | TAA, TAG, TGA | 3 |

### Chemical Properties

**Amino Acid Categories:**
- **Nonpolar, aliphatic**: Ala, Val, Leu, Ile, Met, Pro, Gly
- **Aromatic**: Phe, Tyr, Trp
- **Polar, uncharged**: Ser, Thr, Cys, Asn, Gln
- **Positively charged**: Lys, Arg, His
- **Negatively charged**: Asp, Glu

**Codon Usage Patterns:**
- **GC-rich codons**: Often preferred in GC-rich genomes
- **AT-rich codons**: Often preferred in AT-rich genomes
- **Ending nucleotides**: Third position often shows strongest bias

---

## R Visualization

### Automatic Figure Generation

If R is installed with required packages, two publication-quality figures are automatically generated.

#### Figure 1: RSCU Heatmap

**Features:**
- **Color-coded heatmap**: Visual representation of RSCU values
- **Amino acid grouping**: Codons grouped by amino acid
- **Genome clustering**: Species clustered by similarity
- **Color scale**: Blue (low RSCU) to Red (high RSCU)
- **Professional styling**: Publication-ready

**Purpose:** Compare codon usage patterns across species

#### Figure 2: Codon Usage Comparison

**Features:**
- **Stacked bar charts**: 3 representative species
- **Amino acid facets**: Each amino acid shown separately
- **Proportional visualization**: Relative usage of synonymous codons
- **Professional styling**: Publication-ready

**Purpose:** Show detailed codon preferences for key species

### Figure Specifications

**Format:**
- PDF: Vector format (scalable, perfect for publications)
- PNG: Raster format (600 DPI, high resolution)

**Size:**
- Heatmap: Optimized for number of genomes
- Comparison: 28" × 12" for detailed visualization
- Professional fonts and colors
- Clear axis labels and legends

### Manual R Script Execution

If you want to customize figures:
```bash
# R script is generated in output folder
cd Module8_Codon_Usage_Analysis/

# Edit the R script if desired
nano generate_plots.R

# Run manually
Rscript generate_plots.R
```

---

## Troubleshooting

### Common Issues and Solutions

#### 1. Missing Python Packages
```bash
❌ ERROR: Required package not installed: No module named 'Bio'
```
**Solution:**
```bash
pip install biopython pandas openpyxl

# Verify
python -c "import Bio, pandas, openpyxl"
```

#### 2. No GenBank Files Found
```bash
❌ ERROR: No GenBank files found!
```
**Solution:**
```bash
# Check file extensions
ls *.gb
ls *.gbf
ls *.gbk
ls *.genbank

# Make sure files are in the correct directory
python cgas_module8.py -i /full/path/to/genbank_files/
```

#### 3. No CDS Features Found
```bash
⚠ WARNING: No coding sequences found!
```
**Solution:**
```
Check that your GenBank files contain CDS features:
- Open file in a text editor
- Look for "CDS" in the FEATURES section
- Verify proper start/end coordinates
```

#### 4. Unusual Codons Detected
```bash
⚠ WARNING: Found 5 unusual codon(s): NTG, CGA, ...
```
**Solution:**
```
This indicates non-standard nucleotides in your sequences:
- Check sequencing quality
- Verify genome assembly
- These codons are counted but may affect RSCU accuracy
```

#### 5. R Not Found (Warning Only)
```bash
⚠ R is not installed or not in PATH
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
  Rscript -e "install.packages(c('ggplot2', 'pheatmap', 'RColorBrewer', 'reshape2', 'grid', 'dplyr', 'tidyr', 'zoo'))"
  ```

#### 6. R Packages Missing
```bash
⚠ Missing R packages: ggplot2, pheatmap
```
**Solution:**
```bash
# Install missing packages
Rscript -e "install.packages(c('ggplot2', 'pheatmap', 'RColorBrewer', 'reshape2', 'grid', 'dplyr', 'tidyr', 'zoo'), repos='https://cloud.r-project.org')"
```

### Quick Fixes Table

| Problem | Symptom | Quick Fix |
|---------|---------|-----------|
| Missing packages | Import error | `pip install biopython pandas openpyxl` |
| No GenBank files | "No files found" | Check file extensions (.gb, .gbf, .gbk) |
| No CDS features | "No coding sequences" | Verify GenBank files have CDS features |
| Unusual codons | "Unusual codon(s)" | Check sequence quality |
| R missing | R warning | Install R or use --no-figures |
| R packages missing | Package error | Install with Rscript |

---

## Examples

### Example 1: Complete Workflow
```bash
# 1. Prepare GenBank files
mkdir -p /home/abdullah/codon_analysis/
cd /home/abdullah/codon_analysis/

# 2. Run codon usage analysis
python cgas_module8.py

# 3. Check results
ls Module8_Codon_Usage_Analysis/
open Module8_Codon_Usage_Analysis/Complete_Codon_Usage_Analysis.xlsx
open Module8_Codon_Usage_Analysis/Figures/RSCU_Heatmap.pdf
```

### Example 2: Process Specific Folder
```bash
# GenBank files in separate folder
python cgas_module8.py -i chloroplast_genomes/ -o codon_results/

# Check what was created
ls codon_results/
```

### Example 3: Skip Visualizations (Faster)
```bash
# Generate only Excel files (no R)
python cgas_module8.py --no-figures

# Useful for quick analysis or when R is not available
```

### Example 4: Batch Processing with Script
```bash
#!/bin/bash
# Process multiple datasets

for dataset in monocots dicots ferns; do
    echo "Analyzing $dataset..."
    python cgas_module8.py \
        -i ${dataset}_genomes/ \
        -o ${dataset}_codon_usage/
    echo "Completed $dataset"
done

echo "All datasets analyzed!"
```

### Example 5: From Raw Sequences to Codon Analysis
```bash
# Complete pipeline from raw sequences to codon usage analysis

# Step 1: Annotate genomes (if not already done)
# Use annotation tools like GeSeq, DOGMA, etc.

# Step 2: Run codon usage analysis
python cgas_module8.py

# Step 3: View results
open Module8_Codon_Usage_Analysis/Complete_Codon_Usage_Analysis.xlsx
```

---

## FAQ

### Q1: What's the simplest way to run this?
**A:** Just `python cgas_module8.py` in a folder with GenBank files.

### Q2: Do my GenBank files need to be annotated?
**A:** Yes! Files must contain CDS (coding sequence) features with proper coordinates.

### Q3: What if my GenBank files have unusual codons?
**A:** These are counted but flagged with warnings. Check sequence quality for best results.

### Q4: What if I don't have R installed?
**A:** No problem! Excel reports are still generated. Use `--no-figures` to suppress R warnings. However, If you use simple command still the results are accurate. 

### Q5: What does RSCU value mean?
**A:** Relative Synonymous Codon Usage. RSCU = 1.0 means no bias, >1.0 means preferred codon, <1.0 means avoided codon.

### Q6: Can I analyze mitochondrial genomes?
**A:** Yes! It works with any GenBank files containing CDS features, including mitochondrial or nuclear genomes.

### Q7: How are stop codons handled?
**A:** Stop codons (TAA, TAG, TGA) are included in the analysis but not used for RSCU calculations.

### Q8: Can I process multiple genomes at once?
**A:** Yes! Put all GenBank files in one folder and run the module once. It will process all files and create a merged analysis.

### Q9: What file formats are supported?
**A:** GenBank files with extensions: .gb, .gbf, .gbk, .genbank

### Q10: How do I cite this tool?
**A:** See Citation section below.

### Q11: Why are some codons preferred over others?
**A:** Codon preferences can be due to tRNA abundance, mutational biases, or selection for translational efficiency.

### Q12: How are representative species selected for comparison?
**A:** The script automatically selects three species that represent the diversity of the dataset based on clustering analysis. Because plant lineages generally exhibit high similarity in codon usage patterns, representation by three species is sufficient for visualization in the main article. Data for all analyzed genomes can be provided as supplementary material. In addition, the heatmap presents codon usage information for all species included in the study. 

---

## Technical Specifications

### Performance
- **Processing speed**: ~5 seconds per genome
- **Memory usage**: <500 MB RAM
- **Disk space**: Minimal (<50 MB per analysis)

### Compatibility
- **Python**: 3.9+ (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating Systems**: Windows, macOS, Linux
- **Excel**: Generates .xlsx files (Excel 2010+)
- **R**: 4.0+ (optional, for visualization)

### Input Limits
- **Max genome size**: No practical limit (tested with 200kb+ genomes)
- **Max files**: No practical limit (tested with 50+ files)
- **CDS features**: No limit (all CDS features are processed)

### Quality Features
- ✅ All 61 sense codons analyzed
- ✅ RSCU values calculated for each codon
- ✅ Amino acid grouping maintained
- ✅ Unusual codon detection and reporting
- ✅ Batch processing
- ✅ Merged comparative analysis
- ✅ Professional visualizations

---

## References

### Python Libraries
- **Biopython**: [Official Documentation](https://biopython.org/) - GenBank file parsing
- **pandas**: [Documentation](https://pandas.pydata.org/) - Data manipulation
- **openpyxl**: [Documentation](https://openpyxl.readthedocs.io/) - Excel file creation

### R Packages
- **ggplot2**: [Documentation](https://ggplot2.tidyverse.org/) - Data visualization
- **pheatmap**: [Documentation](https://cran.r-project.org/web/packages/pheatmap/) - Heatmap generation
- **RColorBrewer**: [Documentation](https://cran.r-project.org/web/packages/RColorBrewer/) - Color palettes
- **reshape2**: [Documentation](https://cran.r-project.org/web/packages/reshape2/) - Data transformation

### GenBank Format
- **NCBI**: [GenBank Format Documentation](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/) - File format specifications

### CGAS Suite
- **CGAS**: Chloroplast Genome Analysis Suite - integrated tools for comprehensive chloroplast genomics

### Citation
If using CGAS Module 8 in publications, please cite:
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
python cgas_module8.py --help

# 2. Verify packages are installed
python -c "import Bio, pandas, openpyxl"

# 3. Check GenBank files
ls *.gb *.gbf *.gbk *.genbank

# 4. Verify GenBank format
head -50 your_genbank.gb
```

### Common Issues Solved Here
- ✅ Missing packages? Run `pip install biopython pandas openpyxl`
- ✅ No files found? Check file extensions (.gb, .gbf, .gbk, .genbank)
- ✅ No CDS features? Verify GenBank files have CDS annotations
- ✅ Unusual codons? Check sequence quality
- ✅ R missing? Install R or use --no-figures
- ✅ R packages missing? Install with Rscript

### Quick Reference Card
```bash
# ⚡ QUICK START ⚡
cgas --module 8                                  # cgas command
cgas-codon                                       # shortcut command
python cgas_module8.py                           # python command
python cgas_module8.py -i genbank_files/         # Specify input
python cgas_module8.py -i data/ -o results/      # Custom output
python cgas_module8.py --no-figures              # Skip R graphs

# 🔬 JUPYTER NOTEBOOK 🔬
%run cgas_module8.py                             # %run works for Module 8
!cgas --module 8                                 # ! also works
!cgas-codon --no-figures                         # Skip R in notebook

# 📊 OUTPUT 📊
# Module8_Codon_Usage_Analysis/
# ├── [species]_RSCU.xlsx              # Individual reports
# ├── Complete_Codon_Usage_Analysis.xlsx  # Merged analysis
# ├── Complete_Codon_Usage_Analysis.csv   # R data
# └── Figures/                          # PDF + PNG graphs

# 🎯 TIPS 🎯
# - GenBank files must contain CDS features
# - RSCU = 1.0 means no codon bias
# - RSCU > 1.0 means preferred codon
# - RSCU < 1.0 means avoided codon
# - Chloroplast genomes often show AT bias
```

---

**Last Updated**: January 2026  
**Version**: 1.0.1  
**Author**: Abdullah  
**Part of**: CGAS (Chloroplast Genome Analysis Suite)

---

*Happy Codon Usage Analysis! 🧬✨*