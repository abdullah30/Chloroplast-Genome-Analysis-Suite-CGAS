# CGAS Module 13: Comprehensive Nucleotide Diversity Analysis
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
8. [Nucleotide Diversity Algorithm](#nucleotide-diversity-algorithm)
9. [Region Extraction System](#region-extraction-system)
10. [Gene Notation System](#gene-notation-system)
11. [R Visualization](#r-visualization)
12. [Troubleshooting](#troubleshooting)
13. [Examples](#examples)
14. [FAQ](#faq)
15. [Technical Specifications](#technical-specifications)
16. [References](#references)

---

## Introduction

**CGAS Module 13** is a comprehensive tool for calculating nucleotide diversity (π) in chloroplast genomes from GenBank files. Nucleotide diversity measures the average number of nucleotide differences per site between any two sequences, providing crucial insights into genetic variation and evolutionary patterns.

This module performs complete nucleotide diversity analysis with:
- **Three region types**: Genes, introns, and intergenic spacers
- **Automatic extraction** of coding and non-coding regions
- **Multiple sequence alignment** using MAFFT
- **Publication-ready outputs** with professional formatting and visualizations

### Key Features:
- **Rounded π values** to 1 decimal place (0.0, 0.1, 0.2, etc.) for clarity
- **Enhanced gene synonym handling** for better cross-species comparisons
- **Figures-only mode** to regenerate plots without rerunning analysis
- **Automatic R package installation** if packages are missing
- **LSC/SSC/IR region annotation** in visualizations (when junction positions are consistent)

### Scientific Applications:
- **Evolutionary Studies**: Measure genetic diversity across species
- **Comparative Genomics**: Identify conserved vs. variable regions
- **Population Genetics**: Assess genetic variation within populations
- **Marker Development**: Find polymorphic regions for primer design
- **Functional Genomics**: Compare diversity in coding vs. non-coding regions

### Part of CGAS (Chloroplast Genome Analysis Suite)
CGAS is a comprehensive suite of tools designed for chloroplast genome analysis, providing seamless integration across different analysis modules.

---

## Quick Start - Just Run It!

### ⚡ One Command Does It All!

```bash
# Option 1: Using cgas command
cgas --module 13

# Option 2: Using cgas-diversity shortcut
cgas-diversity

# Option 3: Using python directly
python cgas_module13.py
```

**What happens when you run this:**
1. ✅ Looks for GenBank files in your **current directory**
2. ✅ Creates `Module13_Diversity_Analysis` folder automatically
3. ✅ Processes all `.gb`, `.gbff`, `.genbank`, `.gbk` files
4. ✅ Extracts genes, introns, and intergenic spacers
5. ✅ Aligns sequences using MAFFT
6. ✅ Calculates nucleotide diversity (π)
7. ✅ Generates comprehensive text reports
8. ✅ Creates visualizations if R is installed (auto-installs packages if needed)

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
python cgas_module13.py

# Output created automatically:
# Module13_Diversity_Analysis/
# ├── nucleotide_diversity_results.txt
# ├── nucleotide_diversity_by_position.txt
# ├── genes_for_R_plot.txt
# ├── noncoding_for_R_plot.txt
# ├── genes/ (extracted gene sequences)
# ├── introns/ (extracted intron sequences)
# ├── intergenic/ (extracted intergenic spacers)
# ├── alignments/ (MAFFT alignments)
# └── Figures/ (if R installed)
#     ├── nucleotide_diversity_plot.pdf
#     └── nucleotide_diversity_plot.png (600 DPI)
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
python cgas_module13.py /home/abdullah/data_analysis/chloroplast_genomes/

# Output created in:
# /home/abdullah/data_analysis/chloroplast_genomes/Module13_Diversity_Analysis/
```

#### Example 3: Figures-Only Mode (Regenerate from existing data)
```bash
# Already ran analysis? Just regenerate figures!
python cgas_module13.py --figures-only

# Or specify output folder:
python cgas_module13.py --figures-only -o Module13_Diversity_Analysis/

# This skips extraction, alignment, and diversity calculation
# Only regenerates R plots from existing data files
```

### 🎯 Quick Start Cheat Sheet

```bash
# ====== cgas command ======
cgas --module 13                                 # Process current directory
cgas --module 13 your_data_folder/               # Specify input folder
cgas --module 13 data/ -o results/               # Custom output folder
cgas --module 13 --figures-only                  # Regenerate plots only

# ====== cgas-diversity shortcut ======
cgas-diversity                                   # Process current directory
cgas-diversity your_data_folder/                 # Specify input folder
cgas-diversity data/ -o results/                 # Custom output folder
cgas-diversity --figures-only                    # Regenerate plots only

# ====== python command ======
python cgas_module13.py                          # Process current directory
python cgas_module13.py your_data_folder/        # Specify input folder
python cgas_module13.py data/ -o results/        # Custom output folder
python cgas_module13.py --figures-only           # Regenerate plots only
python cgas_module13.py --help                   # Get help
```

### 📊 What You Get (Output Files)
After running ANY of the above commands, you'll get:

```
Module13_Diversity_Analysis/  # Created automatically
├── 📊 nucleotide_diversity_results.txt             # Complete report
├── 📈 nucleotide_diversity_by_position.txt         # Ordered by position
├── 📁 nucleotide_diversity_all_regions_by_position.txt  # All combined
├── 📄 genes_for_R_plot.txt                         # Data for R (genes)
├── 📄 noncoding_for_R_plot.txt                     # Data for R (non-coding)
├── 🧬 genes/                                       # Extracted gene sequences
├── 🧬 introns/                                     # Extracted intron sequences
├── 🧬 intergenic/                                  # Extracted intergenic spacers
├── 🧬 alignments/                                  # MAFFT alignments
├── 📜 generate_diversity_plot.R                    # R script (editable)
└── 🖼️  Figures/                                   # Publication-quality figures
    ├── nucleotide_diversity_plot.pdf
    └── nucleotide_diversity_plot.png (600 DPI)
```
---

## Jupyter Notebook Usage

> **Note:** Module 13 is a pure Python script, so both `%run` and the `!` shell operator work. Use `%run` for direct script execution, or `!` for the `cgas` / `cgas-diversity` commands. Run the cell from the directory containing your GenBank files, or pass the input folder as a positional argument.

```python
# Using %run (executes the script directly)
%run cgas_module13.py
%run cgas_module13.py genbank_files/
%run cgas_module13.py data/ -o results/
%run cgas_module13.py --figures-only

# Using ! operator with cgas command
!cgas --module 13
!cgas --module 13 genbank_files/
!cgas --module 13 data/ -o results/
!cgas --module 13 --figures-only

# Using ! operator with cgas-diversity shortcut
!cgas-diversity
!cgas-diversity genbank_files/
!cgas-diversity --figures-only

# Using ! operator with python
!python cgas_module13.py
!python cgas_module13.py genbank_files/ -o results/
```

### Advanced: Display Results in Notebook
```python
# After running analysis, display results
from IPython.display import display, Markdown, Image
import pandas as pd

# Display text results
with open('./Module13_Diversity_Analysis/nucleotide_diversity_results.txt', 'r') as f:
    display(Markdown(f"```\n{f.read()}\n```"))

# Load and display genes data
genes_df = pd.read_csv('./Module13_Diversity_Analysis/genes_for_R_plot.txt', sep='\t')
display(genes_df.head(10))

# Display the figure
display(Image('./Module13_Diversity_Analysis/Figures/nucleotide_diversity_plot.png'))
```

---

## Command-Line Usage with Examples

### Complete Command Reference

```bash
# ====================================================================
# USING cgas COMMAND
# ====================================================================

# Process current directory
cgas --module 13

# Specify input directory
cgas --module 13 /home/user/chloroplast_data/

# Custom input and output directories
cgas --module 13 abdullah/data_analysis -o Diversity_Results/

# Figures-only mode
cgas --module 13 --figures-only
```

```bash
# ====================================================================
# USING cgas-diversity SHORTCUT
# ====================================================================

# Process current directory
cgas-diversity

# Specify input directory
cgas-diversity genbank_files/

# With custom output
cgas-diversity data/ -o diversity_results/

# Figures-only mode
cgas-diversity --figures-only
```

```bash
# ====================================================================
# USING python COMMAND
# ====================================================================

# ====================================================================
# BASIC USAGE
# ====================================================================

# 1. Run in current directory (simplest!)
python cgas_module13.py

# 2. Specify input folder
python cgas_module13.py genbank_files/

# 3. Specify BOTH input and output folders
python cgas_module13.py abdullah/data_analysis -o Diversity_Results/

# 4. Figures-only mode (regenerate from existing data)
python cgas_module13.py --figures-only

# 5. Figures-only with custom output folder
python cgas_module13.py --figures-only -o Module13_Diversity_Analysis/

# 6. Get help
python cgas_module13.py --help


# ====================================================================
# REAL-WORLD EXAMPLES WITH PATHS
# ====================================================================

# Example 1: Process data in "chloroplast_data" folder
python cgas_module13.py /home/user/chloroplast_data/

# Example 2: Save results to custom "diversity_results" folder
python cgas_module13.py data/ -o diversity_results/

# Example 3: Regenerate figures from existing analysis
python cgas_module13.py --figures-only -o Module13_Diversity_Analysis/

# Example 4: Windows users
python cgas_module13.py "C:\Users\abdullah\data" -o "C:\Results\diversity"
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
- Complete sequence data
- Feature annotations (genes, CDS, tRNA, rRNA, introns)
- Organism information

**Minimum requirements:**
- At least 2 GenBank files (need multiple samples for diversity calculation)
- Feature annotations with gene names
- Valid DNA sequences (ATCG)

### File Organization
```
your_data_folder/
├── Arabidopsis_thaliana.gb      # Recommended naming
├── Oryza_sativa.gb
├── Zea_mays.gb
└── outgroup_species.gb
```

**Note**: Species names are extracted from the `/organism=` qualifier in GenBank files.

---

## Output Structure

### Directory Structure
```
Module13_Diversity_Analysis/           # Default output folder
├── nucleotide_diversity_results.txt           # Main report
├── nucleotide_diversity_by_position.txt       # Position-ordered
├── nucleotide_diversity_all_regions_by_position.txt  # All combined
├── genes_for_R_plot.txt                       # R data (genes)
├── noncoding_for_R_plot.txt                   # R data (non-coding)
├── genes/                                     # Individual gene FASTA
│   ├── atpF.fasta
│   ├── rbcL.fasta
│   └── ...
├── introns/                                   # Individual intron FASTA
│   ├── atpF_intron.fasta
│   └── ...
├── intergenic/                                # Intergenic spacer FASTA
│   ├── trnH-psbA.fasta
│   └── ...
├── alignments/                                # MAFFT alignments
│   ├── atpF_aligned.fasta
│   └── ...
├── generate_diversity_plot.R                  # R script (editable)
└── Figures/                                   # Publication-quality
    ├── nucleotide_diversity_plot.pdf
    └── nucleotide_diversity_plot.png (600 DPI)
```

### Key Output Files Explained

#### 1. nucleotide_diversity_results.txt
Main report with three sections:
- **Genes**: π values for all protein-coding genes
- **Introns**: π values for all introns
- **Intergenic Spacers**: π values for intergenic regions

Format:
```
Gene/Region Name    Nucleotide Diversity (π)
atpF                0.3
rbcL                0.1
trnH-psbA           2.5
```

#### 2. nucleotide_diversity_by_position.txt
Genes and non-coding regions ordered by their genomic position in the first species.

#### 3. genes_for_R_plot.txt & noncoding_for_R_plot.txt
Tab-separated files for R visualization:
```
Region_Name	Diversity	Type
atpF	0.3	Gene
rbcL	0.1	Gene
atpF_intron	0.5	Intron
trnH-psbA	2.5	Intergenic
```

JSON file containing LSC/SSC/IR boundary markers for visualization (if junction positions are consistent across species).

---

## Detailed Feature Explanation

### 1. Nucleotide Diversity (π) Calculation

**Definition:** π (pi) is the average number of nucleotide differences per site between any two DNA sequences in a sample.

**Formula:**
```
π = Σ(pairwise_differences) / (number_of_comparisons × alignment_length)
```

**Interpretation:**
- π = 0.0: No variation (all sequences identical)
- π = 0.1: Low variation (1 difference per 10 sites)
- π = 1.0: Moderate variation (1 difference per site)
- π = 2.0+: High variation (common in intergenic regions)

### 2. Three Region Types

**Genes (Protein-Coding):**
- Extracted from CDS/gene features
- Generally lower diversity (selection pressure)
- Examples: *atpF*, *rbcL*, *psbA*

**Introns:**
- Extracted from split genes (join locations)
- Moderate diversity (less constrained)
- Examples: *atpF* intron, *ycf3* intron I

**Intergenic Spacers:**
- Regions between annotated features
- Typically highest diversity (least constrained)
- Examples: *trnH*-*psbA*, *trnS*-*trnG*

### 3. Gene Synonym Handling

Module 13 automatically merges gene synonyms:
- *psaA* = *psa1*
- *tRNA-Gly(UCC)* = *trnG-UCC*
- *rps16* = *rps16* (various naming conventions)

This ensures cross-species comparison even when different labs use different gene names.

### 4. Gene Notation System
- Gene names in text: *atpF*, *rbcL*, *psbA* (italicized)
- Species names in text: *Arabidopsis thaliana* (italicized)
- Follows standard biological nomenclature

---

## Nucleotide Diversity Algorithm

### Step-by-Step Process

**Step 1: Region Extraction**
```
For each GenBank file:
  Extract genes → genes/
  Extract introns → introns/
  Extract intergenic spacers → intergenic/
```

**Step 2: Multiple Sequence Alignment**
```
For each region (if present in ≥2 species):
  Collect sequences from all species
  Run MAFFT alignment
  Save to alignments/
```

**Step 3: Diversity Calculation**
```
For each aligned region:
  For each pair of sequences (i, j):
    Count differences
  Calculate average differences per site
  π = total_differences / (n_comparisons × alignment_length)
  Round to 1 decimal place
```

### Handling Missing Data

- **Gaps (-)**: Ignored in pairwise comparisons
- **Ambiguous bases (N)**: Treated as differences
- **Missing regions**: Excluded from analysis for that region

### Minimum Requirements

- At least 2 sequences per region (need pairs for comparison)
- Alignment length ≥ 50 bp (statistical reliability)
- Valid nucleotides (ATCG, gaps allowed)

---

## Region Extraction System

### Gene Extraction

**From CDS/gene features:**
```
gene            1001..1500
                /gene="atpF"
```

**Handling split genes (introns):**
```
gene            join(2001..2500,2701..3200)
                /gene="atpF"
```
- Exon 1: 2001-2500
- Intron: 2501-2700 (extracted separately)
- Exon 2: 2701-3200

### Intron Extraction

**Detected from join() locations:**
- For `join(A..B,C..D)`: intron = B+1 to C-1
- Handles multiple introns per gene
- Names: *atpF_intron*, *ycf3_intron_I*, *ycf3_intron_II*

### Intergenic Spacer Extraction

**Between adjacent features:**
```
Feature 1: trnH (positions 1..72)
Feature 2: psbA (positions 300..1500)
Intergenic spacer: 73..299 (named "trnH-psbA")
```

---

## Gene Notation System

### Formatting Rules (in text)

**Species Names:**
- *Italicized* in documentation
- Full binomial nomenclature
- Example: *Arabidopsis thaliana*, *Oryza sativa*

**Gene Names:**
- *Italicized* when mentioned in text
- Standard chloroplast gene nomenclature
- Examples: *atpF*, *rbcL*, *psbA*, *trnH*

**File Names (no italics):**
- `Arabidopsis_thaliana.gb` (filenames)
- `atpF.fasta` (data files)
- `atpF_aligned.fasta` (alignments)

### Naming Conventions

**Genes:**
- Format: *geneName* (e.g., *atpF*, *rbcL*)
- Lowercase with capital first letter for some (e.g., *rbcL*)

**tRNAs:**
- Format: *trnX-XXX* (e.g., *trnH-GUG*, *trnL-UAA*)
- Three-letter amino acid + anticodon

**Intergenic:**
- Format: *gene1*-*gene2* (e.g., *trnH*-*psbA*)
- Order follows genomic position

---

## R Visualization

### Automatic Figure Generation

If R is installed with required packages (auto-installed if missing), a publication-quality figure is automatically generated.

#### Nucleotide Diversity Plot
- **Panel A**: Gene diversity (coding regions)
- **Panel B**: Non-coding diversity (introns + intergenic spacers)
- **X-axis**: Genomic position (if consistent) or region order
- **Y-axis**: Nucleotide diversity (π)
- **Annotations**: LSC/SSC/IR regions (if junctions consistent)

### Figure Specifications
- **Format**: PDF + PNG (600 DPI)
- **Size**: Optimized for publications (10" × 8")
- **Colors**: Professional color schemes
- **Fonts**: Arial/Helvetica, publication-ready

### LSC/SSC/IR Region Markers

**Displayed when:**
- IR junction positions are consistent across ≥70% of species
- Marker genes detected at boundaries (e.g., *psbA* at LSC start)

**Not displayed when:**
- Junction positions vary significantly
- IR regions absent or unusual structure
- Warning file created: `Figures/IR_junction_variability.txt`

### Manual R Script Execution
```bash
# If you want to customize figures
cd Module13_Diversity_Analysis/
Rscript generate_diversity_plot.R

# Edit the R script first if needed
nano generate_diversity_plot.R
```

### Figures-Only Mode

**Purpose:** Regenerate plots without rerunning the entire analysis (useful for tweaking visualizations).

```bash
# After initial analysis is complete:
python cgas_module13.py --figures-only

# Or specify custom output folder:
python cgas_module13.py --figures-only -o Module13_Diversity_Analysis/
```

**Requirements:**
- Existing `genes_for_R_plot.txt`
- Existing `noncoding_for_R_plot.txt`

---

## Troubleshooting

### Common Issues and Solutions

#### 1. MAFFT Not Found (Critical Error)
```bash
❌ ERROR: MAFFT not found. Please install MAFFT.
```
**Solution:**
```bash
# Ubuntu/Debian
sudo apt-get install mafft

# macOS
brew install mafft

# Conda (all platforms)
conda install -c bioconda mafft

# Verify
mafft --version
```

#### 2. No GenBank Files Found
```bash
❌ ERROR: No GenBank files found in directory
```
**Solution:**
```bash
# Check files exist
ls *.gb

# Verify file extensions
ls -la | grep -E "\.gb$|\.gbff$|\.genbank$|\.gbk$"

# Try with explicit path
python cgas_module13.py /full/path/to/genbank/files/
```

#### 3. R Not Found (Warning Only)
```bash
⚠ R not found. Skipping figure generation.
```
**Solution:**
- This is just a warning - analysis still completes
- Text reports are generated
- Install R if you want figures:
  ```bash
  # Ubuntu/Debian
  sudo apt-get install r-base
  
  # macOS
  brew install r
  ```

#### 4. R Packages Missing
```bash
⚠ Missing R packages: readr, cowplot
→ Attempting automatic installation...
```
**Solution:**
- Script attempts auto-installation
- If fails, install manually:
  ```bash
  Rscript -e "install.packages(c('ggplot2', 'readr', 'cowplot'), repos='https://cloud.r-project.org')"
  ```

#### 5. No Introns Found
```bash
⚠ Warning: No introns found in any species
```
**Solution:**
- This is normal for some chloroplast genomes
- Not all genes have introns
- Analysis continues with genes and intergenic regions

#### 6. Figures-Only Mode Fails
```bash
✗ Error: Required file not found: genes_for_R_plot.txt
  Run the full analysis first
```
**Solution:**
```bash
# Run full analysis first
python cgas_module13.py

# Then use figures-only mode
python cgas_module13.py --figures-only
```

### Quick Fixes Table

| Problem | Symptom | Quick Fix |
|---------|---------|-----------|
| MAFFT missing | "MAFFT not found" | Install MAFFT (see above) |
| No files found | "No GenBank files found" | Check file extensions |
| R missing | "R not found" warning | Install R or ignore |
| R packages missing | Package error | Auto-installs or manual install |
| No diversity | All π = 0.0 | Check if species are identical |
| Alignment fails | MAFFT error | Check sequence quality |

### Get Help
```bash
# Show all options
python cgas_module13.py --help

# Output:
# usage: cgas_module13.py [-h] [-o OUTPUT] [--figures-only] [input_paths ...]
# 
# CGAS Module 13: Nucleotide Diversity Analysis
# 
# positional arguments:
#   input_paths           GenBank files or directories (default: current directory)
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -o OUTPUT, --output OUTPUT
#                         Output directory (default: Module13_Diversity_Analysis)
#   --figures-only        Regenerate figures from existing data (skip analysis)
```

---

## Examples

### Example 1: Complete Workflow
```bash
# 1. Prepare your data
mkdir -p /home/abdullah/chloroplast_diversity/
cp *.gb /home/abdullah/chloroplast_diversity/

# 2. Run analysis (simplest way)
cd /home/abdullah/chloroplast_diversity/
python cgas_module13.py

# 3. Check results
ls Module13_Diversity_Analysis/
less Module13_Diversity_Analysis/nucleotide_diversity_results.txt

# 4. View figure
open Module13_Diversity_Analysis/Figures/nucleotide_diversity_plot.pdf
```

### Example 2: Batch Processing Multiple Datasets
```bash
#!/bin/bash
# Process multiple datasets

for dataset in monocots dicots ferns; do
    echo "Processing $dataset..."
    python cgas_module13.py ${dataset}_genomes/ -o ${dataset}_diversity/
    echo "Completed $dataset"
done

echo "All datasets processed!"
```

### Example 3: Jupyter Notebook Analysis
```python
# In Jupyter notebook

# Cell 1: Run analysis
%run cgas_module13.py -i ./genomes/

# Cell 2: Load and display results
import pandas as pd
from IPython.display import Image, display

genes = pd.read_csv('./Module13_Diversity_Analysis/genes_for_R_plot.txt', sep='\t')
noncoding = pd.read_csv('./Module13_Diversity_Analysis/noncoding_for_R_plot.txt', sep='\t')

print("Top 5 most diverse genes:")
display(genes.nlargest(5, 'Diversity'))

print("\nTop 5 most diverse non-coding regions:")
display(noncoding.nlargest(5, 'Diversity'))

# Cell 3: Display figure
display(Image('./Module13_Diversity_Analysis/Figures/nucleotide_diversity_plot.png'))
```

### Example 4: Regenerate Figures with Custom Edits
```bash
# 1. Run initial analysis
python cgas_module13.py data/

# 2. Edit R script to customize plot
cd Module13_Diversity_Analysis/
nano generate_diversity_plot.R
# (Make your changes)

# 3. Regenerate figures
Rscript generate_diversity_plot.R

# Or use figures-only mode:
cd ..
python cgas_module13.py --figures-only
```

---

## FAQ

### Q1: What's the simplest way to run this?
**A:** Just `python cgas_module13.py` in your data folder. Make sure MAFFT is installed first.

### Q2: Do I need R installed?
**A:** No, R is optional. All text reports are generated without R. R is only needed for the publication-quality figure.

### Q3: What if MAFFT is not installed?
**A:** MAFFT is **required**. Install it with `sudo apt-get install mafft` (Linux) or `brew install mafft` (Mac).

### Q4: How many GenBank files do I need?
**A:** Minimum 2 files. Nucleotide diversity requires comparing multiple sequences.

### Q5: What does π = 0.0 mean?
**A:** All sequences are identical for that region (no variation).

### Q6: What does π = 2.0 mean?
**A:** High variation - on average, 2 differences per nucleotide site when comparing sequences.

### Q7: Can I run this in Jupyter?
**A:** Yes! Use `%run cgas_module13.py` in a Jupyter cell.

### Q8: What is "figures-only" mode?
**A:** It regenerates plots without rerunning the entire analysis. Useful for tweaking visualizations.

### Q9: Why aren't LSC/SSC/IR regions shown in my plot?
**A:** Junction positions vary across your species. Check `Figures/IR_junction_variability.txt` for details.

### Q10: How do I cite this tool?
**A:** 
```
Abdullah. (2026). CGAS Module 13: Comprehensive Nucleotide Diversity Analysis.
Version 1.0.1. Chloroplast Genome Analysis Suite (CGAS).
```

### Q11: What Python version do I need?
**A:** Python 3.7 or higher. Tested on Python 3.7-3.12.

### Q12: Can I analyze nuclear genes?
**A:** The tool works for any GenBank file, but it's optimized for chloroplast genomes.

---

## Technical Specifications

### Performance
- **File size**: Handles typical chloroplast genomes (~150-200kb)
- **Speed**: ~30-120 seconds per species (depending on gene count)
- **Memory**: Moderate (MAFFT alignments, ~500MB-2GB RAM)
- **Output**: Comprehensive text reports + publication figures

### Compatibility
- **Python**: 3.9+ (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating Systems**: Windows 10/11, macOS 10.15+, Linux (Ubuntu 18.04+)
- **MAFFT**: v7.0+ (REQUIRED)
- **R**: 4.0+ (optional, for figures)

### Limits
- **Min species**: 2 (need for pairwise comparison)
- **Max species**: No practical limit (tested with 100+ species)
- **Alignment length**: No limit (MAFFT handles long sequences)
- **Regions**: Processes all genes/introns/intergenic found

### Quality Features
- ✅ Gene synonym normalization (cross-species consistency)
- ✅ Automatic MAFFT alignment
- ✅ Rounded π values (1 decimal place)
- ✅ LSC/SSC/IR annotation (when consistent)
- ✅ Figures-only mode (rapid plot regeneration)
- ✅ Automatic R package installation
- ✅ Comprehensive error handling

---

## References

### GenBank Format
- [NCBI GenBank Format Guide](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)
- [INSDC Feature Table Documentation](https://www.insdc.org/)

### Nucleotide Diversity
- **Nei & Li (1979)**: Mathematical model for studying genetic variation in terms of restriction endonucleases.
- **MAFFT**: Katoh & Standley (2013) Multiple alignment program
- **Chloroplast genomes**: Typical structure and gene content

### Related Tools
- **MAFFT**: [Multiple sequence alignment](https://mafft.cbrc.jp/)
- **Biopython**: [Python bioinformatics tools](https://biopython.org/)
- **ggplot2**: [R data visualization](https://ggplot2.tidyverse.org/)
- **DnaSP**: Desktop software for population genetics analysis
- **MEGA**: Molecular evolutionary genetics analysis
- **CPStool**: Determine nucleotide diversity

### Python Libraries
- **BioPython**: [Official Documentation](https://biopython.org/)
- **pandas**: [User Guide](https://pandas.pydata.org/docs/)
- **NumPy**: [Reference](https://numpy.org/doc/)

### R Packages (Optional)
- **ggplot2**: [Documentation](https://ggplot2.tidyverse.org/)
- **readr**: [Documentation](https://readr.tidyverse.org/)
- **cowplot**: [Documentation](https://wilkelab.org/cowplot/)

### CGAS Suite
- **CGAS**: Chloroplast Genome Analysis Suite - integrated tools for comprehensive chloroplast genomics

### Citation
If using CGAS Module 13 in publications, please cite:

**CGAS Suite:**
```
Abdullah, Yan, R., & Tian, X. (2025). CGAS (Chloroplast Genome Analysis Suite): 
An Automated Python Pipeline for Comprehensive Comparative Chloroplast Genomics. 
bioRxiv. https://doi.org/10.64898/2025.12.21.695765
```

**MAFFT (used by Module 13):**
```
Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software 
version 7: improvements in performance and usability. Molecular Biology and 
Evolution, 30(4), 772-780.
```

---

## Support & Contact

### Getting Help
```bash
# 1. First check built-in help
python cgas_module13.py --help

# 2. Review this documentation
# 3. Check example outputs in Examples section
# 4. Verify dependencies are installed
pip list | grep -E "biopython|pandas|numpy"
mafft --version
```

### Common Issues Solved Here
- ✅ MAFFT not found? Install MAFFT first (REQUIRED)
- ✅ No files found? Use `-i` flag with correct path
- ✅ R warnings? Install R or ignore (text reports still generated)
- ✅ Alignment issues? Check GenBank file quality and annotations
- ✅ No diversity? Check if all sequences are identical
- ✅ Figures-only fails? Run full analysis first

### Quick Reference Card
```bash
# ⚡ QUICK START ⚡
cgas --module 13                                 # cgas command
cgas-diversity                                   # shortcut command
python cgas_module13.py                          # python command
python cgas_module13.py data/                    # Specify input
python cgas_module13.py -o results/              # Specify output
python cgas_module13.py --figures-only           # Regenerate plots only

# 🔬 JUPYTER NOTEBOOK 🔬
%run cgas_module13.py                            # %run works for Module 13
!cgas --module 13                                # ! also works
!cgas-diversity --figures-only                   # Regenerate in notebook

# 📊 OUTPUT 📊
# Module13_Diversity_Analysis/
# ├── nucleotide_diversity_results.txt     # Main report
# ├── nucleotide_diversity_by_position.txt  # Position-ordered
# ├── genes_for_R_plot.txt, noncoding_for_R_plot.txt  # R data
# ├── genes/, introns/, intergenic/        # Extracted sequences
# ├── alignments/                          # MAFFT alignments
# └── Figures/                             # PDF + PNG

# 🎯 TIPS 🎯
# - Install MAFFT first (REQUIRED for alignment)
# - Need at least 2 GenBank files
# - R is optional (auto-installs packages)
# - Use --figures-only to regenerate plots
# - π values rounded to 1 decimal place
# - Gene names follow biological nomenclature (italic in text)
```

---

**Last Updated**: January 2026  
**Version**: 1.0.1  
**Author**: Abdullah  
**Part of**: CGAS (Chloroplast Genome Analysis Suite)

---

*Happy Diversity Analysis! 🧬✨*
