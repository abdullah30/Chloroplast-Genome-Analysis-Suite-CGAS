# CGAS Module 14: Phylogenetic Matrix Builder with IQ-TREE Integration
## Complete Documentation and User Guide

---

## 📋 Table of Contents
1. [Introduction](#introduction)
2. [Quick Start - Just Run It!](#quick-start---just-run-it)
3. [Installation Guide](#installation-guide)
4. [Command-Line Usage with Examples](#command-line-usage-with-examples)
5. [Jupyter Notebook Usage](#jupyter-notebook-usage)
6. [Input Requirements](#input-requirements)
7. [Output Structure](#output-structure)
8. [Detailed Feature Explanation](#detailed-feature-explanation)
9. [Alignment Methods](#alignment-methods)
10. [IQ-TREE Phylogenetic Inference](#iq-tree-phylogenetic-inference)
11. [Gene Notation System](#gene-notation-system)
12. [Region Types](#region-types)
13. [Troubleshooting](#troubleshooting)
14. [Examples](#examples)
15. [FAQ](#faq)
16. [Technical Specifications](#technical-specifications)
17. [References](#references)

---

## Introduction

**CGAS Module 14** is a comprehensive tool for building phylogenetic matrices and inferring phylogenetic trees from chloroplast genomes. This module extracts genomic regions, aligns them, creates concatenated matrices, and performs phylogenetic inference using IQ-TREE - all in one automated workflow.

This module performs complete phylogenetic analysis with:
- **Four matrix types**: Genes-only, introns-only, intergenic-only, and complete concatenated
- **Two alignment methods**: MAFFT (default) and MACSE (codon-aware for CDS)
- **Automatic phylogenetic inference** with IQ-TREE
- **Partitioned analysis** with automatically generated partition files
- **Professional output organization** with all files properly categorized

### Key Features:
- **Automatic region extraction**: *atpA*, *rbcL*, *matK*, introns, intergenic spacers
- **Smart gene name normalization**: Handles tRNA naming variants and gene synonyms
- **Multiple sequence alignment**: MAFFT for all regions, MACSE for codon-aware CDS alignment
- **Concatenated matrices**: Builds supermatrices from aligned regions
- **IQ-TREE integration**: Automatic phylogenetic inference with support values
- **Outgroup specification**: Place outgroups at top of trees and alignments
- **Multi-threading**: Utilizes all CPU cores for faster analysis

### Scientific Applications:
- **Phylogenetic reconstruction**: Build species trees from chloroplast genomes
- **Comparative genomics**: Study evolutionary relationships
- **Marker development**: Identify phylogenetically informative regions
- **Evolutionary studies**: Analyze patterns of molecular evolution
- **Systematics**: Support taxonomic classifications

### Part of CGAS (Chloroplast Genome Analysis Suite)
CGAS is a comprehensive suite of tools designed for chloroplast genome analysis, providing seamless integration across different analysis modules.

---

## Quick Start - Just Run It!

### ⚡ One Command Does It All!

```bash
# Option 1: Using cgas command
cgas --module 14

# Option 2: Using cgas-phylogeny shortcut
cgas-phylogeny

# Option 3: Using python directly
python cgas_module14.py
```

**What happens when you run this:**
1. ✅ Looks for GenBank files in your **current directory**
2. ✅ Creates `Module14_Phylogeny` folder automatically
3. ✅ Processes all `.gb`, `.gbk` files
4. ✅ Extracts genes, introns, and intergenic spacers
5. ✅ Aligns sequences using MAFFT
6. ✅ Generates 4 concatenated matrices
7. ✅ Creates partition files in 4 formats
8. ✅ Organizes outputs in professional directory structure

### 📁 Real-World Folder Examples

#### Example 1: Basic Usage (Easiest!)
```bash
# Your folder structure:
abdullah/
├── research_project/
│   ├── Arabidopsis_thaliana.gb
│   ├── Oryza_sativa.gb
│   └── Zea_mays.gb

# Navigate to your data folder
cd /home/abdullah/research_project/

# Run the script (no arguments needed!)
python cgas_module14.py

# Output created automatically:
# Module14_Phylogeny/
# ├── 01_individual_alignments/
# ├── 02_concatenated_matrices/
# ├── partitions/
# └── SUMMARY_REPORT.txt
```

#### Example 2: With Phylogenetic Inference
```bash
# Generate matrices AND infer phylogenetic trees!
python cgas_module14.py --iqtree

# This will:
# 1. Generate all matrices (genes-only, introns-only, etc.)
# 2. Run IQ-TREE on each matrix
# 3. Generate trees with support values (UFBoot/SH-aLRT)
# 4. Create comprehensive output reports
```

#### Example 3: Complete Matrix Only with Phylogeny (RECOMMENDED)
```bash
# Skip individual matrices, just build complete matrix and run phylogeny
python cgas_module14.py --complete-only --iqtree

# This is the RECOMMENDED approach for phylogenetic analysis:
# - Builds complete concatenated matrix (all regions)
# - Runs partitioned IQ-TREE analysis
# - Generates maximum likelihood tree with support values
# - Fastest workflow for phylogeny
```

#### Example 4: With Outgroup Specification
```bash
# Specify outgroup for proper tree rooting
python cgas_module14.py --iqtree -og Gossypium_arboreum.gb

# Multiple outgroups also supported:
python cgas_module14.py --iqtree -og outgroup1.gb outgroup2.gb outgroup3.gb
```

### 🎯 Quick Start Cheat Sheet

```bash
# ====== cgas command ======
cgas --module 14                                           # Generate matrices only
cgas --module 14 --iqtree                                  # With phylogenetic inference
cgas --module 14 --complete-only --iqtree                  # Complete matrix + phylogeny (RECOMMENDED)
cgas --module 14 --genes-only --iqtree                     # Genes only + phylogeny
cgas --module 14 --iqtree -og Nicotiana_tabacum.gb         # With outgroup
cgas --module 14 --macse --iqtree                          # MACSE codon-aware alignment
cgas --module 14 --iqtree -o my_phylogeny_results/         # Custom output
cgas --module 14 --iqtree --threads 8                      # Custom thread count

# ====== cgas-phylogeny shortcut ======
cgas-phylogeny                                             # Generate matrices only
cgas-phylogeny --iqtree                                    # With phylogenetic inference
cgas-phylogeny --complete-only --iqtree                    # Complete matrix + phylogeny (RECOMMENDED)
cgas-phylogeny --genes-only --iqtree -og outgroup.gb       # With outgroup
cgas-phylogeny --macse --iqtree                            # MACSE alignment

# ====== python command ======
python cgas_module14.py                                    # Generate matrices only
python cgas_module14.py --iqtree                           # With phylogenetic inference
python cgas_module14.py --complete-only --iqtree           # Complete matrix + phylogeny (RECOMMENDED)
python cgas_module14.py --genes-only --iqtree              # Genes only + phylogeny
python cgas_module14.py --iqtree -og Nicotiana_tabacum.gb  # With outgroup
python cgas_module14.py --macse --iqtree                   # MACSE codon-aware alignment
python cgas_module14.py --iqtree -o my_phylogeny_results/  # Custom output
python cgas_module14.py --iqtree --threads 8               # Custom thread count
python cgas_module14.py --help                             # Get help
```

### 📊 What You Get (Output Files)
After running ANY of the above commands, you'll get:

```
Module14_Phylogeny/  # Created automatically
├── 📄 SUMMARY_REPORT.txt                    # Comprehensive analysis summary
├── 📁 01_individual_alignments/            # Individual region alignments
│   ├── genes/                              # Gene alignments
│   │   ├── atpA_aligned.fasta
│   │   ├── rbcL_aligned.fasta
│   │   └── ...
│   ├── introns/                            # Intron alignments
│   └── intergenic/                         # Intergenic spacer alignments
├── 📁 02_concatenated_matrices/            # Concatenated supermatrices
│   ├── matrix_protein_coding.fasta         # Genes-only matrix
│   ├── matrix_introns.fasta                # Introns-only matrix
│   ├── matrix_intergenic.fasta             # Intergenic-only matrix
│   ├── matrix_complete.fasta               # Complete matrix (RECOMMENDED)
│   ├── partitions_protein_coding.txt       # RAxML format
│   ├── partitions_protein_coding.nex       # NEXUS format
│   ├── partitions_protein_coding_iqtree.nex # IQ-TREE format
│   ├── partitions_protein_coding.phylip    # PHYLIP format
│   └── ... (partition files for each matrix)
└── 🌳 If --iqtree used:                    # Phylogenetic trees
    ├── phylogeny_protein_coding.treefile   # ML tree with supports
    ├── phylogeny_protein_coding.log        # Analysis log
    ├── phylogeny_protein_coding.iqtree     # Full IQ-TREE report
    ├── phylogeny_complete.treefile         # Complete matrix tree
    └── ...
```

---

## Installation Guide

> **Note:** If you have already installed the full CGAS tool (using `environment.yml` or `environment-minimal.yml`), all dependencies for Module 14 are already included — no additional installation except MACSE alinger is needed. The steps below are only for running Module 14 as a standalone script.

### Prerequisites
- **Python 3.9 or higher** (tested on Python 3.9–3.12)
- **pip** (Python package manager, usually included with Python)
- **MAFFT** (multiple sequence alignment tool) - **REQUIRED**
- **Optional**: IQ-TREE 2.0+ (for phylogenetic inference)
- **Optional**: MACSE (for codon-aware CDS alignment)

### Step-by-Step Installation

#### 1. Install Python Dependencies (One Command)
```bash
pip install biopython pandas numpy
```

**Recommended versions:**
```bash
# Or with specific versions
pip install biopython>=1.79 pandas>=1.3.0 numpy>=1.21.0
```

#### 2. Install MAFFT (Required - Critical!)
```bash
# Ubuntu/Debian
sudo apt-get install mafft

# macOS
brew install mafft

# Windows: Download from https://mafft.cbrc.jp/alignment/software/
# Or use conda (all platforms): conda install -c bioconda mafft

# Verify installation
mafft --version
# Should show: MAFFT v7.x or higher
```

**Why MAFFT is required:** Module 14 aligns sequences for concatenated matrix building. Without MAFFT, the module cannot create alignments.

#### 3. Install IQ-TREE (Optional, but recommended for phylogeny)
```bash
# Ubuntu/Debian
sudo apt-get install iqtree

# macOS
brew install iqtree

# Windows: Download from http://www.iqtree.org/
# Or use conda: conda install -c bioconda iqtree

# Verify installation
iqtree --version
# Should show: IQ-TREE multicore version 2.x
```

**IQ-TREE versions supported:** iqtree, iqtree2, iqtree-omp, iqtree2-omp (auto-detected)

#### 4. Install MACSE (Optional, for codon-aware alignment)
```bash
# Download from: https://bioweb.supagro.inra.fr/macse/
# Place macse.jar in your working directory or in PATH

# Verify installation
java -jar macse.jar -help
```

**Note:** MACSE only works with CDS (protein-coding genes). When using `--macse`, introns and intergenic spacers are skipped.

#### 5. Download the Script
```bash
# Download cgas_module14.py to your working directory
# Or copy from your source
```

#### 6. Verify Installation
```bash
# Test Python dependencies
python cgas_module14.py --help

# Test MAFFT
mafft --version

# Test IQ-TREE (if installed)
iqtree --version 2>/dev/null || echo "IQ-TREE not installed (optional)"
```

### Dependency Details
| Package | Version | Purpose | Required |
|---------|---------|---------|----------|
| **biopython** | ≥1.79 | GenBank parsing, sequence handling | Yes |
| **pandas** | ≥1.3.0 | Data manipulation | Yes |
| **numpy** | ≥1.21.0 | Numerical operations | Yes |
| **MAFFT** | v7.0+ | Multiple sequence alignment | **Yes** |
| **IQ-TREE** | v2.0+ | Phylogenetic inference | Optional |
| **MACSE** | v2.0+ | Codon-aware CDS alignment | Optional |

---

## Jupyter Notebook Usage

> **Note:** Module 14 is a pure Python script, so both `%run` and the `!` shell operator work. Use `%run` for direct script execution, or `!` for the `cgas` / `cgas-phylogeny` commands. Run the cell from the directory containing your GenBank files, or pass the input folder with `-i`.

```python
# Using %run (executes the script directly)
%run cgas_module14.py
%run cgas_module14.py --iqtree
%run cgas_module14.py --complete-only --iqtree
%run cgas_module14.py --iqtree -og Gossypium_arboreum.gb
%run cgas_module14.py -i data/ -o results/ --iqtree

# Using ! operator with cgas command
!cgas --module 14
!cgas --module 14 --iqtree
!cgas --module 14 --complete-only --iqtree
!cgas --module 14 --iqtree -og outgroup.gb
!cgas --module 14 --macse --iqtree

# Using ! operator with cgas-phylogeny shortcut
!cgas-phylogeny
!cgas-phylogeny --iqtree
!cgas-phylogeny --complete-only --iqtree
!cgas-phylogeny --genes-only --iqtree -og outgroup.gb

# Using ! operator with python
!python cgas_module14.py
!python cgas_module14.py --iqtree -i data/ -o results/
```

### Advanced: Display Results in Notebook
```python
# After running analysis, display results
from IPython.display import display, Markdown
from Bio import Phylo
import matplotlib.pyplot as plt
import os

# Display summary report
with open('./Module14_Phylogeny/SUMMARY_REPORT.txt', 'r') as f:
    display(Markdown(f"```\n{f.read()}\n```"))

# Display tree (if generated)
if os.path.exists('./Module14_Phylogeny/phylogeny_complete.treefile'):
    tree = Phylo.read('./Module14_Phylogeny/phylogeny_complete.treefile', 'newick')
    fig, ax = plt.subplots(figsize=(12, 8))
    Phylo.draw(tree, axes=ax, do_show=False)
    plt.title('Phylogenetic Tree - Complete Matrix', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()
```

---

## Command-Line Usage with Examples

### Complete Command Reference

```bash
# ====================================================================
# USING cgas COMMAND
# ====================================================================

# Generate matrices only
cgas --module 14

# Generate matrices + phylogenetic inference
cgas --module 14 --iqtree

# Complete matrix only + phylogeny (RECOMMENDED)
cgas --module 14 --complete-only --iqtree

# Genes only + phylogeny
cgas --module 14 --genes-only --iqtree

# With outgroup specification
cgas --module 14 --iqtree -og outgroup.gb

# MACSE codon-aware alignment
cgas --module 14 --macse --iqtree

# Custom output folder
cgas --module 14 --iqtree -o my_results/
```

```bash
# ====================================================================
# USING cgas-phylogeny SHORTCUT
# ====================================================================

# Generate matrices only
cgas-phylogeny

# Generate matrices + phylogenetic inference
cgas-phylogeny --iqtree

# Complete matrix only + phylogeny (RECOMMENDED)
cgas-phylogeny --complete-only --iqtree

# Genes only + phylogeny with outgroup
cgas-phylogeny --genes-only --iqtree -og outgroup.gb

# MACSE codon-aware alignment
cgas-phylogeny --macse --iqtree
```

```bash
# ====================================================================
# USING python COMMAND
# ====================================================================

# ====================================================================
# BASIC USAGE
# ====================================================================

# 1. Generate matrices only (no phylogeny)
python cgas_module14.py

# 2. Generate matrices for complete genome, introns, IGS, and genes + phylogenetic inference
python cgas_module14.py --iqtree

# 3. Complete matrix only + phylogeny (RECOMMENDED for population genetic or closely related species of a genus)
python cgas_module14.py --complete-only --iqtree

# 4. Genes only + phylogeny
python cgas_module14.py --genes-only --iqtree #for only matrice
python cgas_module14.py --macse --genes-only --iqtree
python cgas_module14.py --macse --genes-only --iqtree -og outgroup.gb #with outgroup
# 5. Introns only
python cgas_module14.py --introns-only #for only matrice
python cgas_module14.py --introns-only --iqtree


# 6. Intergenic only
python cgas_module14.py --intergenic-only #for only matrice
python cgas_module14.py --intergenic-only --iqtree
# 7. With outgroup specification
python cgas_module14.py --iqtree -og outgroup.gb

# 8. Custom output folder
python cgas_module14.py --iqtree -o my_results/

# 9. Custom threads
python cgas_module14.py --iqtree --threads 16

# 10. Get help
python cgas_module14.py --help


# ====================================================================
# REAL-WORLD EXAMPLES
# ====================================================================

# Example 1: Complete phylogenetic workflow (RECOMMENDED)
python cgas_module14.py --complete-only --iqtree -og Gossypium_arboreum.gb

# Example 2: Traditional genes-only approach
python cgas_module14.py --genes-only --iqtree

# Example 3: Use MACSE for codon-aware CDS alignment
python cgas_module14.py --genes-only --macse --iqtree

# Example 4: Process specific directory
python cgas_module14.py -i /home/user/chloroplast_data/ --iqtree

# Example 5: Multiple outgroups
python cgas_module14.py --iqtree -og outgroup1.gb outgroup2.gb outgroup3.gb

# Example 6: Windows users
python cgas_module14.py --iqtree -i "C:\Users\abdullah\data"
```

---

## Input Requirements

### Supported File Formats
- `.gb` (GenBank)
- `.gbk` (GenBank)
- `.gbff` (GenBank Flat File)
- `.genbank`

### GenBank File Requirements
Standard GenBank files from NCBI with:
- Complete sequence data
- Feature annotations (genes, CDS, tRNA, rRNA, introns)
- Organism information

**Minimum requirements:**
- At least 2 GenBank files (need multiple samples for phylogeny)
- Feature annotations with gene names
- Valid DNA sequences (ATCG)

### File Organization
```
your_data_folder/
├── Arabidopsis_thaliana.gb      # Recommended naming
├── Oryza_sativa.gb
├── Zea_mays.gb
├── outgroup_species.gb           # Can specify as outgroup
└── ...
```

**Note**: Organism names are extracted from the `/organism=` qualifier in GenBank files and used for sequence naming in alignments and trees.

---

## Output Structure

### Directory Structure (Default: Module14_Phylogeny)
```
Module14_Phylogeny/
├── SUMMARY_REPORT.txt                          # Main analysis summary
├── 01_individual_alignments/                   # Individual region alignments
│   ├── genes/                                  # Gene alignments
│   │   ├── atpA_aligned.fasta
│   │   ├── atpF_aligned.fasta
│   │   ├── rbcL_aligned.fasta
│   │   └── ... (all extracted genes)
│   ├── introns/                                # Intron alignments
│   │   ├── atpF_intron_aligned.fasta
│   │   ├── ycf3_intron_I_aligned.fasta
│   │   └── ...
│   └── intergenic/                             # Intergenic spacer alignments
│       ├── trnH-psbA_aligned.fasta
│       ├── trnS-trnG_aligned.fasta
│       └── ...
├── 02_concatenated_matrices/                   # Concatenated matrices
│   ├── matrix_protein_coding.fasta             # Genes only (CDS)
│   ├── matrix_introns.fasta                    # Introns only
│   ├── matrix_intergenic.fasta                 # Intergenic only
│   ├── matrix_complete.fasta                   # ALL REGIONS (RECOMMENDED)
│   ├── partitions_protein_coding.txt           # RAxML partition format
│   ├── partitions_protein_coding.nex           # NEXUS partition format
│   ├── partitions_protein_coding_iqtree.nex    # IQ-TREE partition format
│   ├── partitions_protein_coding.phylip        # PHYLIP partition format
│   └── ... (partition files for each matrix)
└── IQ-TREE output (if --iqtree used):
    ├── phylogeny_protein_coding.treefile       # ML tree (Newick)
    ├── phylogeny_protein_coding.log            # Analysis log
    ├── phylogeny_protein_coding.iqtree         # Full report
    ├── phylogeny_protein_coding.model.gz       # Model parameters
    ├── phylogeny_protein_coding.mldist         # ML distances
    ├── phylogeny_complete.treefile             # Complete matrix tree
    └── ...
```

### Key Output Files Explained

#### 1. SUMMARY_REPORT.txt
Comprehensive analysis summary containing:
- List of processed species
- Extracted genes and regions
- Matrix statistics (taxa, alignment length)
- Partition information
- Outgroup ordering (if specified)
- Command used

#### 2. Individual Alignments (01_individual_alignments/)
- Separate FASTA file for each extracted region
- Aligned with MAFFT or MACSE
- Named by gene/region (*atpF*, *rbcL*, *trnH-psbA*, etc.)

#### 3. Concatenated Matrices (02_concatenated_matrices/)
Four matrix types:
- **matrix_protein_coding.fasta**: Only protein-coding genes
- **matrix_introns.fasta**: Only introns
- **matrix_intergenic.fasta**: Only intergenic spacers
- **matrix_complete.fasta**: All regions combined (RECOMMENDED)

#### 4. Partition Files
Four formats for each matrix:
- **.txt**: RAxML format
- **.nex**: NEXUS format (MrBayes, PAUP*)
- **_iqtree.nex**: IQ-TREE format
- **.phylip**: PHYLIP format

#### 5. Phylogenetic Trees (if --iqtree used)
- **.treefile**: Maximum likelihood tree with support values
- **.log**: Analysis log with model selection details
- **.iqtree**: Complete IQ-TREE report
- **.model.gz**: Best-fit model parameters
- **.mldist**: ML distance matrix

---

## Detailed Feature Explanation

### 1. Automatic Region Extraction

**Genes (Protein-Coding):**
- Extracted from CDS/gene features
- Examples: *atpA*, *atpF*, *rbcL*, *matK*, *ndhF*, *psbA*
- Handles split genes (with introns)

**Introns:**
- Detected from join() locations in GenBank
- Named: *atpF_intron*, *ycf3_intron_I*, *ycf3_intron_II*
- Multiple introns per gene supported

**Intergenic Spacers:**
- Regions between adjacent features
- Named: *trnH-psbA*, *trnS-trnG*, *atpF-atpH*
- Ordered by genomic position

### 2. Gene Name Normalization

Module 14 handles naming variants:
- tRNA genes: *trnH-GUG*, *tRNA-His(GUG)*, *tRNAHis* → normalized
- Gene synonyms: *psaA*/*psa1*, *rps16* variants
- Case sensitivity handled

### 3. Outgroup Placement

When outgroups specified with `-og`:
- Outgroups placed at TOP of all alignments
- Outgroups placed FIRST in all matrices
- Order preserved as specified
- Helps with tree rooting

### 4. Sequence Naming

**In alignments and matrices:**
- Uses organism names from GenBank: *Arabidopsis thaliana*
- NOT accession numbers
- Handles spaces and special characters
- Maintains consistency across all files

---

## Alignment Methods

### MAFFT (Default)

**When used:**
- Default for all region types
- Fast and accurate
- Handles all sequence types

**Command:**
```bash
# MAFFT used automatically for all regions
python cgas_module14.py
```

**Advantages:**
- Fast (seconds per gene)
- Handles large datasets
- Works with genes, introns, intergenic

### MACSE (Optional - Codon-Aware)

**When used:**
- Only with `--macse` flag
- Only processes CDS (protein-coding genes)
- Skips introns and intergenic regions

**Command:**
```bash
# Use MACSE for CDS alignment
python cgas_module14.py --macse
```

**Advantages:**
- Codon-aware alignment
- Preserves reading frames
- Better for protein-coding genes
- Recommended for genes-only phylogeny

**Limitations:**
- Only works with CDS
- Requires Java
- Slower than MAFFT
- Cannot process introns/intergenic

**When to use MACSE:**
- Genes-only phylogeny
- Codon-aware analysis needed
- Protein phylogeny
- Selection analysis

**When to use MAFFT:**
- Complete matrix (default, recommended)
- Including introns/intergenic
- Faster processing needed
- Standard phylogenetic analysis

---

## IQ-TREE Phylogenetic Inference

### Automatic Phylogenetic Analysis

**Enable with `--iqtree` flag:**
```bash
python cgas_module14.py --iqtree
```

### IQ-TREE Settings Used

**Model Selection:**
- `-m MFP`: ModelFinder Plus (automatic best-fit model selection)
- Tests all standard models
- Selects best by BIC/AIC

**Support Values:**
- `-bb 1000`: 1000 ultrafast bootstrap (UFBoot) replicates
- `-alrt 1000`: 1000 SH-aLRT test replicates
- Both displayed on tree (UFBoot/SH-aLRT format)

**Optimization:**
- `-bnni`: Optimize UFBoot trees by nearest neighbor interchange
- Reduces impact of severe model violations

**Threading:**
- `-nt AUTO`: Automatic thread detection (default)
- Or specify: `--threads 8`

**Partitioning:**
- `-p partition_file`: Partitioned analysis
- Each gene/region gets own model
- More accurate than concatenated single-model

### IQ-TREE Output Files

For each matrix analyzed:

**phylogeny_[matrix].treefile:**
- Maximum likelihood tree
- Newick format
- Support values: UFBoot/SH-aLRT
- Can open in FigTree, iTOL, etc.

**phylogeny_[matrix].log:**
- Analysis log
- Model selection results
- Log-likelihood scores
- Running time

**phylogeny_[matrix].iqtree:**
- Complete IQ-TREE report
- Full model details
- Tree statistics
- Bootstrap convergence

**phylogeny_[matrix].model.gz:**
- Best-fit model parameters
- Compressed format

**phylogeny_[matrix].mldist:**
- ML distance matrix
- Pairwise distances

### Support Value Interpretation

**UFBoot (Ultrafast Bootstrap):**
- ≥95%: Strong support
- 80-95%: Moderate support
- <80%: Weak support

**SH-aLRT (Shimodaira-Hasegawa approximate likelihood ratio test):**
- ≥80%: Strong support
- 70-80%: Moderate support
- <70%: Weak support

**On tree displayed as:** `UFBoot/SH-aLRT`
- Example: `100/100` = maximum support
- Example: `85/92` = moderate/strong support

---

## Gene Notation System

### Formatting Rules

**In Text:**
- Species names: *Arabidopsis thaliana*, *Oryza sativa* (italic)
- Gene names: *atpF*, *rbcL*, *psbA*, *trnH* (italic)
- Standard biological nomenclature

**In Filenames:**
- `Arabidopsis_thaliana.gb` (no italics - filename)
- `atpF_aligned.fasta` (no italics - file)
- `matrix_complete.fasta` (no italics - data file)

**In Trees:**
- Organism names from GenBank
- *Arabidopsis thaliana* (as labeled)
- Spaces and special characters handled

---

## Region Types

### 1. Protein-Coding Genes (CDS)

**Characteristics:**
- Translate to proteins
- Generally conserved (selection pressure)
- Most chloroplast genes

**Examples:**
- *atpA*, *atpF*, *atpH* (ATP synthase)
- *rbcL* (RuBisCO large subunit)
- *matK* (maturase K)
- *ndhF*, *ndhH* (NADH dehydrogenase)
- *psbA*, *psbD* (photosystem II)

**Uses:**
- Traditional phylogenetics
- Protein phylogeny
- Functional studies

### 2. Introns

**Characteristics:**
- Removed during RNA processing
- Less conserved than exons
- Not all genes have introns

**Examples:**
- *atpF_intron*
- *ycf3_intron_I*, *ycf3_intron_II*
- *rps12_intron*

**Uses:**
- Resolve recent divergences
- Population genetics
- Lower-level phylogeny

### 3. Intergenic Spacers

**Characteristics:**
- Between genes
- Least conserved (neutral)
- Fastest evolution

**Examples:**
- *trnH-psbA* (most variable)
- *trnS-trnG*
- *atpF-atpH*

**Uses:**
- Species delimitation
- Population genetics
- DNA barcoding
- Recent divergences

### 4. Complete Matrix (RECOMMENDED)

**Includes:**
- All genes
- All introns
- All intergenic spacers

**Advantages:**
- Maximum phylogenetic signal
- Multiple evolutionary rates
- Partitioned analysis
- Best resolution

**When to use:**
- Modern phylogenomics
- All-evidence approach
- Conflicting signals
- Deep + shallow phylogeny

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

#### 2. IQ-TREE Not Found (Warning)
```bash
⚠ IQ-TREE not found. Skipping phylogenetic inference.
```
**Solution:**
- This is a warning - matrices still generated
- Install IQ-TREE if you want phylogeny:
  ```bash
  # Ubuntu/Debian
  sudo apt-get install iqtree
  
  # macOS
  brew install iqtree
  
  # Conda
  conda install -c bioconda iqtree
  ```

#### 3. No GenBank Files Found
```bash
❌ ERROR: No GenBank files found in directory
```
**Solution:**
```bash
# Check files exist
ls *.gb

# Verify file extensions
ls -la | grep -E "\.gb$|\.gbk$"

# Try with explicit path
python cgas_module14.py -i /full/path/to/files/
```

#### 4. MACSE Not Working
```bash
⚠ MACSE not found or failed
```
**Solution:**
```bash
# Download MACSE
wget https://bioweb.supagro.inra.fr/macse/releases/macse_v2.06.jar

# Rename to macse.jar
mv macse_v2.06.jar macse.jar

# Test
java -jar macse.jar -help
```

#### 5. Outgroup Not Found
```bash
❌ ERROR: Outgroup file not found: outgroup.gb
```
**Solution:**
```bash
# Check file exists
ls -la outgroup.gb

# Use full path if needed
python cgas_module14.py --iqtree -og /full/path/to/outgroup.gb

# Or use filename if in same directory
python cgas_module14.py --iqtree -og Gossypium_arboreum.gb
```

#### 6. IQ-TREE Timeout
```bash
✗ IQ-TREE exceeded time limit (2 hours)
```
**Solution:**
- Dataset too large
- Try with fewer partitions
- Use complete-only for faster analysis
- Increase timeout in code if needed

### Quick Fixes Table

| Problem | Symptom | Quick Fix |
|---------|---------|-----------|
| MAFFT missing | "MAFFT not found" | Install MAFFT (required) |
| IQ-TREE missing | "IQ-TREE not found" | Install IQ-TREE or skip --iqtree |
| No files found | "No GenBank files" | Check .gb/.gbk extensions |
| Outgroup issue | "Outgroup not found" | Use full path to file |
| MACSE fails | MACSE warning | Check Java, download MACSE |
| Memory error | Script crashes | Process fewer files |

### Get Help
```bash
# Show all options
python cgas_module14.py --help

# Output:
# usage: cgas_module14.py [-h] [-i INPUT] [-o OUTPUT] [--iqtree] 
#                         [--genes-only] [--introns-only] [--intergenic-only]
#                         [--complete-only] [--macse] [-og OUTGROUPS [OUTGROUPS ...]]
#                         [--threads THREADS]
# 
# CGAS Module 14: Phylogenetic Matrix Builder with IQ-TREE Integration
```

---

## Examples

### Example 1: Complete Workflow (RECOMMENDED)
```bash
# 1. Prepare data
mkdir -p /home/abdullah/phylogeny/
cp *.gb /home/abdullah/phylogeny/

# 2. Run complete analysis with phylogeny
cd /home/abdullah/phylogeny/
python cgas_module14.py --complete-only --iqtree -og outgroup.gb

# 3. Check results
ls Module14_Phylogeny/02_concatenated_matrices/
cat Module14_Phylogeny/SUMMARY_REPORT.txt

# 4. View tree
# Open phylogeny_complete.treefile in FigTree
```

### Example 2: Traditional Genes-Only Approach
```bash
# Use only protein-coding genes with MACSE
python cgas_module14.py --genes-only --macse --iqtree

# Good for:
# - Traditional chloroplast phylogeny
# - Protein-based analysis
# - Codon-aware alignment
```

### Example 3: Batch Processing Multiple Datasets
```bash
#!/bin/bash
# Process multiple taxonomic groups

for group in monocots dicots ferns; do
    echo "Processing $group..."
    python cgas_module14.py \
        -i ${group}_genomes/ \
        -o ${group}_phylogeny/ \
        --complete-only \
        --iqtree \
        -og ${group}_outgroup.gb
    echo "Completed $group"
done

echo "All phylogenies completed!"
```

### Example 4: Jupyter Notebook Analysis
```python
# In Jupyter notebook

# Cell 1: Run analysis
%run cgas_module14.py --complete-only --iqtree

# Cell 2: Display results
from IPython.display import display, Markdown
with open('./Module14_Phylogeny/SUMMARY_REPORT.txt', 'r') as f:
    display(Markdown(f"```\n{f.read()}\n```"))

# Cell 3: Visualize tree
from Bio import Phylo
import matplotlib.pyplot as plt

tree = Phylo.read('./Module14_Phylogeny/phylogeny_complete.treefile', 'newick')
fig, ax = plt.subplots(figsize=(14, 10))
Phylo.draw(tree, axes=ax)
plt.title('Phylogenetic Tree - Complete Matrix', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.show()
```

---

## FAQ

### Q1: What's the simplest way to run this?
**A:** `python cgas_module14.py` generates matrices. Add `--iqtree` for phylogeny.

### Q2: Which matrix should I use for phylogeny?
**A:** Complete matrix (`--complete-only --iqtree`) is RECOMMENDED. Includes all information.

### Q3: Do I need IQ-TREE installed?
**A:** No, IQ-TREE is optional. Without it, you get matrices that you can analyze separately.

### Q4: What's the difference between MAFFT and MACSE?
**A:** MAFFT is fast and works for all regions. MACSE is codon-aware but only for CDS (genes).

### Q5: How do I specify an outgroup?
**A:** Use `-og outgroup.gb` flag. Multiple outgroups: `-og og1.gb og2.gb og3.gb`

### Q6: How many GenBank files do I need?
**A:** Minimum 2 files for alignment and phylogeny. More is better (typically 10-100+).

### Q7: Can I run this in Jupyter?
**A:** Yes! Use `%run cgas_module14.py --iqtree` in a Jupyter cell.

### Q8: What are support values on the tree?
**A:** Format is UFBoot/SH-aLRT. Values ≥95/80 indicate strong support.

### Q9: How long does analysis take?
**A:** 
- Matrix generation: Minutes (depends on genome count)
- IQ-TREE: Minutes to hours (depends on matrix size and partitions). The module 14 just facilitate the phylogeny work. The phylogenetic inference depend on the work of iqtree. So, please follow iqtree documentation for this. 

### Q10: How do I cite this tool?
**A:** See Citation section below.

### Q11: What Python version do I need?
**A:** Python 3.7 or higher. Tested on Python 3.7-3.12.

### Q12: Can I use nuclear genes?
**A:** Module is optimized for chloroplast genomes, but will work with any GenBank file.

---

## Technical Specifications

### Performance
- **File size**: Handles typical chloroplast genomes (~150-200kb)
- **Alignment speed**: ~1-5 seconds per gene with MAFFT
- **Matrix generation**: ~1-10 minutes (depending on genome count)
- **IQ-TREE analysis**: ~10 minutes to 2 hours (depends on matrix size)
- **Memory**: Moderate (~1-4GB RAM for typical datasets)

### Compatibility
- **Python**: 3.9+ (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating Systems**: Windows 10/11, macOS 10.15+, Linux (Ubuntu 18.04+)
- **MAFFT**: v7.0+ (REQUIRED)
- **IQ-TREE**: v2.0+ (optional, auto-detects iqtree/iqtree2 variants)
- **MACSE**: v2.0+ (optional, requires Java)

### Limits
- **Min species**: 2 (need for alignment)
- **Max species**: No practical limit (tested with 200+ species)
- **Max alignment length**: No limit (tested with 200kb+ alignments)
- **IQ-TREE timeout**: 2 hours (can be modified in code)

### Quality Features
- ✅ Gene name normalization (cross-species consistency)
- ✅ Automatic MAFFT/MACSE alignment
- ✅ Four partition file formats
- ✅ Outgroup ordering support
- ✅ Multi-threading (IQ-TREE)
- ✅ Comprehensive error handling
- ✅ Professional output organization

---

## References

### GenBank Format
- [NCBI GenBank Format Guide](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)
- [INSDC Feature Table Documentation](https://www.insdc.org/)

### Alignment Software
- **MAFFT**: [Katoh & Standley (2013)](https://doi.org/10.1093/molbev/mst010) - Multiple alignment program
- **MACSE**: [Ranwez et al. (2011)](https://doi.org/10.1371/journal.pone.0022594) - Codon-aware alignment

### Phylogenetic Software
- **IQ-TREE**: [Nguyen et al. (2015)](https://doi.org/10.1093/molbev/msu300) - Fast ML phylogenetic inference
- **IQ-TREE 2**: [Minh et al. (2020)](https://doi.org/10.1093/molbev/msaa015) - Version 2 improvements
- **ModelFinder**: [Kalyaanamoorthy et al. (2017)](https://doi.org/10.1038/nmeth.4285) - Model selection

### Chloroplast Phylogenomics
- **Complete plastid genomes**: Revolutionized plant phylogeny
- **Partitioned analysis**: Improved phylogenetic inference
- **Concatenated matrices**: Standard in modern phylogenomics

### Python Libraries
- **BioPython**: [Official Documentation](https://biopython.org/)
- **pandas**: [User Guide](https://pandas.pydata.org/docs/)
- **NumPy**: [Reference](https://numpy.org/doc/)

### Related Tools
- **RAxML**: Maximum likelihood phylogeny software
- **MrBayes**: Bayesian phylogenetic inference
- **PAUP**: Phylogenetic analysis using parsimony
- **ASTRAL**: Species tree estimation

### CGAS Suite
- **CGAS**: Chloroplast Genome Analysis Suite - integrated tools for comprehensive chloroplast genomics

### Citation
If using CGAS Module 14 in publications, please cite:
```
Abdullah, Yan, R., & Tian, X. (2025). CGAS (Chloroplast Genome Analysis Suite): 
An Automated Python Pipeline for Comprehensive Comparative Chloroplast Genomics. 
bioRxiv. https://doi.org/10.64898/2025.12.21.695765
```

**Also cite the tools used:**

**MAFFT (always used):**
```
Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software 
version 7: improvements in performance and usability. Molecular Biology and 
Evolution, 30(4), 772-780.
```

**IQ-TREE (if --iqtree used):**
```
Nguyen, L.T., Schmidt, H.A., von Haeseler, A., & Minh, B.Q. (2015). IQ-TREE: 
A fast and effective stochastic algorithm for estimating maximum likelihood 
phylogenies. Molecular Biology and Evolution, 32(1), 268-274.

Minh, B.Q., Schmidt, H.A., Chernomor, O., Schrempf, D., Woodhams, M.D., 
von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient 
methods for phylogenetic inference in the genomic era. Molecular Biology and 
Evolution, 37(5), 1530-1534.
```

**MACSE (if --macse used):**
```
Ranwez, V., Harispe, S., Delsuc, F., & Douzery, E.J. (2011). MACSE: Multiple 
Alignment of Coding SEquences accounting for frameshifts and stop codons. 
PLoS ONE, 6(9), e22594.
```

---

## Support & Contact

### Getting Help
```bash
# 1. First check built-in help
python cgas_module14.py --help

# 2. Review this documentation
# 3. Check example outputs in Module14_Phylogeny/
# 4. Verify dependencies are installed
pip list | grep -E "biopython|pandas|numpy"
mafft --version
iqtree --version 2>/dev/null || echo "IQ-TREE not installed"
```

### Common Issues Solved Here
- ✅ MAFFT not found? Install MAFFT first (REQUIRED)
- ✅ IQ-TREE missing? Install IQ-TREE or omit --iqtree flag
- ✅ No files found? Check file extensions (.gb, .gbk)
- ✅ Outgroup issues? Use full path or check filename
- ✅ MACSE fails? Check Java installation and MACSE download
- ✅ Memory issues? Process fewer genomes or use complete-only

### Quick Reference Card
```bash
# ⚡ QUICK START ⚡
cgas --module 14                                      # cgas command
cgas-phylogeny                                        # shortcut command
python cgas_module14.py                               # python command
python cgas_module14.py --iqtree                      # Matrices + phylogeny
python cgas_module14.py --complete-only --iqtree      # Complete matrix (RECOMMENDED)

# 🔬 JUPYTER NOTEBOOK 🔬
%run cgas_module14.py                                 # %run works for Module 14
!cgas --module 14 --iqtree                            # ! also works
!cgas-phylogeny --complete-only --iqtree -og out.gb   # Complete workflow

# 🌳 PHYLOGENY OPTIONS 🌳
# Complete matrix (all regions) - RECOMMENDED
python cgas_module14.py --complete-only --iqtree

# Genes only (traditional)
python cgas_module14.py --genes-only --iqtree

# Genes with MACSE (codon-aware)
python cgas_module14.py --genes-only --macse --iqtree

# With outgroup
python cgas_module14.py --complete-only --iqtree -og outgroup.gb

# 📊 OUTPUT 📊
# Module14_Phylogeny/
# ├── 01_individual_alignments/
# ├── 02_concatenated_matrices/
# │   ├── matrix_complete.fasta          (RECOMMENDED)
# │   └── partitions_complete_iqtree.nex
# ├── phylogeny_complete.treefile        (ML tree)
# └── SUMMARY_REPORT.txt

# 🎯 TIPS 🎯
# - Install MAFFT first (REQUIRED)
# - Install IQ-TREE for phylogeny (recommended)
# - Use --complete-only for best results
# - Specify outgroup for proper rooting
# - Complete matrix recommended over genes-only
# - Gene names follow biological nomenclature (italic in text)
```

---

**Last Updated**: January 2026  
**Version**: 1.0.1  
**Author**: Abdullah  
**Part of**: CGAS (Chloroplast Genome Analysis Suite)

---

*Happy Phylogenetic Analysis! 🌳✨*
