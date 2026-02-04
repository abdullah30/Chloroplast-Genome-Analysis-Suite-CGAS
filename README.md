# CGAS - Chloroplast Genome Analysis Suite

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Version](https://img.shields.io/badge/version-1.0.1-green.svg)](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS/releases)
[![Platform](https://img.shields.io/badge/platform-linux%20%7C%20macos%20%7C%20windows-lightgrey.svg)](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS)

> A comprehensive command-line toolkit for chloroplast genome assembly, annotation, and comparative analysis.

**CGAS analyzes 50+ chloroplast genomes in minutes, generating publication-ready results through 14 integrated modules.**

---

## Table of Contents

- [🌟 Features](#-features)
- [📋 Modules Overview](#-modules-overview)
- [🚀 Quick Start](#-quick-start)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
  - [External Dependencies](#external-dependencies)
- [💻 Usage Guide](#-usage-guide)
  - [Command Line Interface](#command-line-interface)
  - [Jupyter Notebook Usage](#jupyter-notebook-usage)
- [📖 Detailed Module Commands](#-detailed-module-commands)
  - [Module 1: Chloroplast Genome Assembly](#module-1-chloroplast-genome-assembly)
  - [Module 2: Plastome Annotation](#module-2-plastome-annotation)
  - [Module 14: Phylogenetic Matrix Builder](#module-14-phylogenetic-matrix-builder)
- [🐛 Troubleshooting](#-troubleshooting)
- [🤝 Contributing](#-contributing)
- [📄 Citation](#-citation)
- [📜 License](#-license)
- [🙏 Acknowledgments](#-acknowledgments)

---

## 🌟 Features

✨ **Complete Workflow** - From raw sequencing reads to publication-ready phylogenetic trees

🔧 **14 Specialized Modules** - Each module handles a specific aspect of chloroplast genome analysis

🚀 **Automated Pipeline** - Streamlined analysis with minimal user intervention

📊 **Publication-Ready Outputs** - Generate Excel tables, Word documents, and high-quality figures

🐍 **Python-Based** - Works seamlessly with Jupyter notebooks

💻 **Cross-Platform** - Compatible with Linux, macOS, and Windows

🎯 **Outgroup Support** - Automatic handling of outgroup species for phylogenetic analysis

📈 **Visualization Integration** - R-based visualizations for modules where researchers provide figures during publication

---

## 📋 Modules Overview

### Preparation Modules (1-4) - Run Individually

| Module | Name | Input | Output | Description |
|--------|------|-------|--------|-------------|
| **1** | Genome Assembly & QC in Batch | FASTQ (raw reads) | Assembled genomes + QC reports + clean reads + coverage analysis | Quality analysis of raw data (comprehensive QC reports ) → adapter trimming → quality filtering → clean reads → chloroplast genome assembly → coverage depth analysis → Detail summary|
| **2** | Annotation | FASTA | GenBank files | Annotates plastomes using PGA |
| **3** | Gene Comparison | GenBank | Normalized GenBank | Standardizes gene names across genomes |
| **4** | Format Conversion | GenBank | FASTA/TBL for NCBI | Converts to NCBI submission format |

 Main Analysis Modules (5-14) - Run with `cgas 5,6,7,8,9,10,11,12,13` ##Run phylogeny separately Ifyou do not to include outgroup species in the comparative analysis. 

| Module | Name | Purpose |
|--------|------|---------|
| **5** | Gene Comparative Analysis | Compare gene content across species |
| **6** | Gene Content Tables | Publication-ready tables (Word) |
| **7** | Genome Structure Analysis | LSC/SSC/IR regions and functional genes, GC content |
| **8** | Codon Usage (RSCU) | Codon bias analysis with visualization |
| **9** | Amino Acid Analysis | Amino acid composition patterns with visulization |
| **10** | SNP Analysis | Substitution detection with quality visulization|
| **11** | Intron Analysis | Gene and tRNA introns data and complete comparison|
| **12** | SSR Analysis | Detect microsatellites, assign to genomic/functional regions, compare motifs, and generate high-quality visualization|
| **13** | Nucleotide Diversity | Assessment of nucleotide diversity (π) with high-quality visualizations |
| **14** | Phylogenetic Analysis | Build matrices and reconstruct trees using IQ-TREE |

---

## 🚀 Quick Start

### Prerequisites

- **Python 3.9 or higher** (recommended: Python 3.9+)
- Required bioinformatics tools (see [External Dependencies](#external-dependencies))
- R (optional, for visualization modules 8, 9, 12, and 13)
- Java 8+ (optional, MACSE optional in module 14 for codon-aware alignment)

### Installation

#### Option 1: Complete Conda Environment (Recommended)

Includes all dependencies - ideal for servers and high-performance systems.

```bash
# Clone the repository
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas
cd cgas

# Create conda environment with all dependencies
conda env create -f environment.yml
conda activate cgas

# Install CGAS
pip install -e .

# Verify installation
cgas --list
```

**Installs:** Python 3.9, all dependencies, fastp, GetOrganelle, BLAST+, MAFFT, MACSE, IQ-TREE, R with 12 packages

#### Option 2: Minimal Conda Environment

Lightweight setup ideal for standard laptops and desktops. Works for modules 2-14.

```bash
# Clone the repository
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas
cd cgas

# Create minimal conda environment
conda env create -f environment-minimal.yml
conda activate cgas-minimal

# Install CGAS
pip install -e .

# Verify installation
cgas --list
```

#### Option 3: Manual Installation (Without Conda)

```bash
# Clone the repository
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas
cd cgas

# Create virtual environment (recommended)
python3.9 -m venv cgas_env
source cgas_env/bin/activate  # On Windows: cgas_env\Scripts\activate

# Install Python dependencies
pip install -e .
pip install -r requirements.txt

# Verify installation
cgas --list
```

**Note:** External bioinformatics tools must be installed separately (see below).

### External Dependencies

Required tools for specific modules:

#### Ubuntu/Debian
```bash
# Basic tools: For detailed installation instructions, please refer to each tool’s official website.
sudo apt-get update
sudo apt-get install python3-dev python3-pip r-base mafft iqtree ncbi-blast+

# Additional tools for Assembly (Module 1)
sudo apt-get install fastp samtools bwa

# Optional tools
sudo apt-get install fasttree
```

#### macOS (using Homebrew)
```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install tools
brew install python@3.9 r mafft iqtree blast fastp samtools bwa fasttree

# Java (for MACSE)
brew install openjdk
```

#### Additional Tools (manual installation)

1. **GetOrganelle** (for Module 1):
   ```bash
   pip install get_organelle
   ```

2. **PGA** (for Module 2):
   ```bash
   git clone https://github.com/quxiaojian/PGA.git
   # Follow installation instructions in PGA repository
   ```

3. **MACSE** (optional, for Module 14):
   - Download from: https://bioweb.supagro.inra.fr/macse/
   - Requires Java 8 or higher

**Full installation guide:** [INSTALL.md](INSTALL.md)

---

## 💻 Usage Guide

### Command Line Interface

After installation, you can use the `cgas` command:

```bash
# Show help
cgas --help

# List all available modules
cgas --list

# PHASE 1: Preparation Modules (run individually)
cgas --module 1 -i raw_reads/ -o results/
# --module 1       : Run Module 1 (raw read processing and genome assembly)
# -i raw_reads/    : Input directory containing paired-end or single-end raw reads
# -o results/      : Output directory where all results will be saved
# Note: Keep the output directory outside of the input directory. 
# CGAS efficiently processes multiple samples from raw reads to assembled genomes 
# and quality-controlled repeats in a single command, avoiding repeated runs.
#Results produced in the following structure
# Output created automatically:
# results/
# ├── 01_QC/                           # Fastp reports
# ├── 02_clean_reads/                  # Cleaned reads
# ├── 03_assemblies/                   # GetOrganelle outputs
# ├── 04_mapping/                      # BWA mapping results
# ├── 05_cp_reads/                     # Chloroplast reads
# ├── 06_reports/                      # Summary reports
# └── 07_assembled_genomes/            # Complete assemblies only

# Default: skip samples that are already processed
cgas module 1 -i raw_reads/ -o results/ --skip-existing
#For detail look to command line content

# 1. Annotate genomes using reference genomes
cgas --module 2 -i 7_assembled_genomes/ -r reference_genomes/  --pga /home/yanrushan/tools/PGA/PGA.pl #you should use your actualy path
cgas --module 2 -i 7_assembled_genomes/ -r reference_genomes/ #If you follow the installation document for pga then the path is not needed. Look at documentation. 
# --module 2               : Run Module 2 (genome annotation)
# -i 7_assembled_genomes/  : Input directory containing assembled genomes from Module 1
# -r reference_genomes/    : Directory containing reference genomes for annotation

# 2. Annotate genomes with organism name mapping
cgas --module 2 -i 7_assembled_genomes/ -r reference_genomes/ --organism-file organisms.txt  --pga /home/yanrushan/tools/PGA/PGA.pl #you should use your actualy path
# --organism-file organisms.txt : Optional file mapping sample IDs to organism names
# Useful for proper labeling in downstream analyses
#This is also done in external tools such as Geneious
# - The script manually curates the LSC and SSC regions based on IR regions identified by PGA.
# - Annotations of the reference genome are expected to be of very high quality.
# - For more details on the impact of reference selection and curation steps, 
#   see the PGA GitHub documentation.

# Run Module 3 with a reference genome
cgas --module 3 -r reference.gb #The file name should not include space for example Hibiscus syriacus.fasta is not acceptable. This should Hibscus_syriacus.fasta
# --module 3        : Run Module 3 (gene annotation comparison and normalization)
# -r reference.gb   : Reference genome file used for comparison
#                     All genomes to be analyzed should be placed in a single folder.
#                     One genome should be chosen as the reference.
#                     Reference filename should use underscores instead of spaces, e.g.,
#                     Abutilon_grandifolium.gb, NOT Abutilon grandifolium.gb.
# Module Description:
# CGAS Module 3 standardizes chloroplast genome gene annotations across multiple genomes.
# - Standardizes gene names across different annotation conventions
# - Detects missing or extra genes
# - Validates product descriptions
# - Shows intron presence/absence and coding region lengths
# - Provides statistics about genes
# - Generates comprehensive comparison reports
# This ensures consistent, high-quality annotation across all genomes in the dataset.
# Can also be applied to genomes downloaded from NCBI to carefully curate differences in gene naming and annotations.

#### Module 4: Annotation Validation and NCBI Submission Preparation

# Run Module 4 to validate annotations and generate submission files

cgas --module 4

# --module 4        : Validate genome annotations and prepare files required for NCBI submission
#                     Run this in the directory containing all finalized GenBank files.
#                     The module detects annotation errors, predicts issues, and checks overall quality.
#                     Processing is fast (1–2 minutes for ~20 genomes).

# Notes:
# - Apply this separately to genomes intended for NCBI submission, 
#   as it generates combined files for batch submission.
# - For genomes where you only want to validate annotations (not submit to NCBI),
#   you can still run this module for quality checks.
# - This module streamlines submission, replacing the need to upload single species manually.
 

### 🎯 PHASE 2: Main Analysis

cd all_genomes/
# Place all genomes in a single directory for comparative analysis.
# Do not include the outgroup genome in this folder if you do not want it included in the analysis.

# Run specific modules
cgas --module 5
cgas --module 8
#same for 5-13; 

cgas --module 14 --macse --genes-only --iqtree  -og outgroup_species.gb #To run with macse aligner
cgas --module 14 --genes-only --iqtree  -og outgroup_species.gb #you can specify one or more outgroup
#                                                               Reference filename should use underscores instead of spaces
#                                                               e.g., Abutilon_grandifolium.gb, NOT Abutilon grandifolium.gb.
# You can specify multiple outgroups; for details, see [Module 14 command line]



# Run multiple modules
cgas --modules 5,6,7
cgas --modules 8,9,12,13
```

### Jupyter Notebook Usage

For Jupyter notebook users, CGAS can be used in several ways:

#### Method 1: Using %run magic command for individual modules

```python
# In a Jupyter notebook cell

# For Modules 3-14 (works directly); Install all dependencies
%pip install biopython>=1.79 pandas>=1.3.0 openpyxl>=3.1 numpy>=1.24 python-docx>=0.8.11

# Running Modules 3 to 14
# After installing dependencies, you can now run Modules 3 to 14.
# Ensure that MAFFT and IQ-TREE are installed and available in your system PATH.
# Alternatively, Windows users can run CGAS via WSL2 with the minimal conda environment (`cgas-minimal`) as described above.
# A tutorial video is available specifically for WSL2 users; please watch it for step-by-step guidance.
%run cgas_module3.py -r reference.gb
%run cgas_module5.py
%run cgas_module8.py #follow same for all. 
%run cgas_module14.py --genes-only --macse --iqtree

# Running Modules 1 and 2

# Modules 1 and 2 require special handling and use specialized software.
# Linux or macOS systems are recommended for best compatibility.
# WSL2 may work on Windows, but compatibility is not guaranteed—GetOrganelle in particular may not function reliably.

```
```bash
# Direct module commands
cgas-assembly          # Module 1: Genome Assembly
cgas-annotate          # Module 2: Plastome Annotation
cgas-compare           # Module 3: Gene Comparison
cgas-convert           # Module 4: Format Conversion
cgas-gene-compare      # Module 5: Gene Comparative Analysis
cgas-gene-table        # Module 6: Gene Content Tables
cgas-genome-compare    # Module 7: Genome Structure Analysis
cgas-codon             # Module 8: Codon Usage (RSCU)
cgas-amino             # Module 9: Amino Acid Analysis
cgas-snp               # Module 10: SNP Analysis
cgas-intron            # Module 11: Intron Analysis
cgas-ssr               # Module 12: SSR Analysis
cgas-diversity         # Module 13: Nucleotide Diversity
cgas-phylogeny         # Module 14: Phylogenetic Analysis
```

---

## 📖 Detailed Module Commands

### Module 1: Chloroplast Genome Assembly

**Assemble chloroplast genomes from raw FASTQ sequencing reads.**

#### Command Line

```bash

### 🎯 Module 1: Raw Reads Processing & Genome Assembly

# 1. Basic assembly
cgas --module 1 -i /home/abdullah/raw_reads/ -o /home/abdullah/results

# 2. Use more threads for faster processing
cgas --module 1 -i raw_reads/ -o results/ -t 16

# 3. Skip already processed samples when adding new data
cgas --module 1 -i raw_reads/ -o results/ --skip-existing

# 4. Custom k-mer values for difficult samples
cgas --module 1 -i raw_reads/ -o results/ -k "21,33,55,77,99"

# 5. More GetOrganelle rounds for challenging assemblies
cgas --module 1 -i raw_reads/ -o results/ -R 20

# 6. Force mapping of incomplete assemblies
cgas --module 1 -i raw_reads/ -o results/ --force-mapping
# 7. Use nohup for large datasets when running on a server
cgas --module 1 -i large_dataset/ -o assemblies/ --use-nohup

# 8. Direct assembly command (shortcut)
cgas-assembly -i raw_reads/ -o results/

# With additional options
cgas-assembly \
  -i raw_reads/ \
  -o assemblies/ \
  -t 16 \
  -F embplant_pt \
  -k 21,45,65,85,105 \
  -R 15 \
  --trim-poly-g
  #For more details on Module 1 parameters and usage, please refer to the Module 1 README.md in the Documents folder.  
  #You can also run Module 1 with additional options by executing the Python script within the CGAS conda environment.

```
---

### Module 2: Plastome Annotation

**Annotate assembled plastomes using PGA (Plastid Genome Annotator).**

#### Command Line

```bash
# Using cgas command
cgas --module 2 -i  7_assembled_genomes/ -r reference_genomes/  --pga /home/yanrushan/tools/PGA/PGA.pl

# Using direct command
cgas-annotate -i  7_assembled_genomes/ -r reference_genomes/  --pga /home/yanrushan/tools/PGA/PGA.pl

# With custom reference and organism names
cgas-annotate \
  -i assemblies/ \
  -o annotations/ \
  -r custom_reference.gb \
  --organism-file organisms.txt \
  -t 16  --pga /home/yanrushan/tools/PGA/PGA.pl
```

#### Options

```
Required:
  -i, --input           Input directory with FASTA files
  -o, --output          Output directory for annotations

Optional:
  -r, --reference       Custom reference GenBank file
  -t, --threads         Number of CPU threads (default: 4)
  --organism-file       Text file with organism names (one per line)
```

---

### Module 3: Plastome Gene Comparison

**Normalize gene names across multiple plastomes for consistency.**

#### Command Line

```bash
# Using cgas command
cgas --module 3 -r /home/abdullah/reference/Abutilon_grandifolium.gb #all the genomes in the folder will be consider as target

# Using direct command
cgas-compare -r reference/Abutilon_grandifolium.gb

# With specific target directory
cgas-compare -r reference.gb
```

#### Jupyter Notebook

```python
# Module 3 works with %run
%run cgas_module3.py -r /home/abdullah/reference/Abutilon_grandifolium.gb

# Or with cgas command
!cgas --module 3 -r reference/Abutilon_grandifolium.gb
```

#### Options

```
Required:
  -r, --reference       Reference GenBank file for gene names

Optional:
  -t, --targets         Directory with target GenBank files (default: current dir) 
  --output              Output directory (default: Module3_GeneComparison/)
```


### Module 14: Phylogenetic Matrix Builder

**Build concatenated gene matrices and phylogenetic trees with MACSE and IQ-TREE.**

#### Command Line

```bash
# Using cgas command with MACSE (codon-aware alignment, recommended)
cgas --module 14 --macse --genes-only --iqtree -og outgroup.gb

# Using direct command
cgas-phylogeny --macse --genes-only --iqtree -og reference/Abutilon_grandifolium.gb

# Complete genome phylogeny
cgas-phylogeny --complete-genome --iqtree

# Custom gene selection
cgas-phylogeny --genes matK,rbcL,atpB,ndhF --macse --iqtree -og outgroup.gb

# MAFFT alignment instead of MACSE
cgas-phylogeny --genes-only --iqtree -og outgroup.gb
```

#### Jupyter Notebook

```python
# Module 14 works with %run
%run cgas_module14.py --genes-only --macse --iqtree -og /home/abdullah/reference/outgroup.gb

# Or with cgas command
!cgas --module 14 --macse --genes-only --iqtree -og outgroup.gb
```

#### Options

```
Optional:
  -i, --input           Input directory with GenBank files
  -o, --output          Output directory (default: Module14_Phylogeny/)
  -og, --outgroup       Outgroup GenBank file(s)
  --genes-only          Use only protein-coding genes
  --complete-genome     Use complete genome sequences
  --genes               Comma-separated list of specific genes
  --macse               Use MACSE for codon-aware alignment (recommended)
  --iqtree              Run IQ-TREE for phylogenetic inference
  --complete-only       Only use complete matrix (no missing data)
  --format              Output format (fasta, phylip, nexus)
```

#### Output Files

```

---

## 🐛 Troubleshooting

For troubleshooting specific to individual modules, please refer to the README file in each module's directory. Each module has detailed documentation covering:

- Common errors and solutions
- Parameter optimization tips
- Input/output format requirements
- Performance tuning guidelines

**General Support:**
- Check module-specific README files for detailed troubleshooting
- Visit [GitHub Issues](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS/issues) for reported problems
- Create a new issue with your error message and system details

---

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

---

## 📄 Citation

If you use CGAS in your research, please cite:

**Preferred Citation (bioRxiv):**
```
Abdullah, Rushan Yan, Xiaoxuan Tian (2025). CGAS (Chloroplast Genome Analysis Suite): 
An Automated Python Pipeline for Comprehensive Comparative Chloroplast Genomics. 
bioRxiv. doi: https://doi.org/10.64898/2025.12.21.695765
```

**Alternative Citation Formats:**

*APA Style:*
```
Abdullah, Yan, R., & Tian, X. (2025). CGAS (Chloroplast Genome Analysis Suite): 
An Automated Python Pipeline for Comprehensive Comparative Chloroplast Genomics. 
bioRxiv. https://doi.org/10.64898/2025.12.21.695765
```

*MLA Style:*
```
Abdullah, Rushan Yan, and Xiaoxuan Tian. "CGAS (Chloroplast Genome Analysis Suite): 
An Automated Python Pipeline for Comprehensive Comparative Chloroplast Genomics." 
bioRxiv (2025). https://doi.org/10.64898/2025.12.21.695765
```

*Chicago Style:*
```
Abdullah, Rushan Yan, and Xiaoxuan Tian. 2025. "CGAS (Chloroplast Genome Analysis Suite): 
An Automated Python Pipeline for Comprehensive Comparative Chloroplast Genomics." 
bioRxiv. https://doi.org/10.64898/2025.12.21.695765
```

**BibTeX:**
```bibtex
@article{abdullah2025cgas,
  title={CGAS (Chloroplast Genome Analysis Suite): An Automated Python Pipeline for Comprehensive Comparative Chloroplast Genomics},
  author={Abdullah and Yan, Rushan and Tian, Xiaoxuan},
  journal={bioRxiv},
  year={2025},
  doi={10.64898/2025.12.21.695765},
  url={https://doi.org/10.64898/2025.12.21.695765}
}
```

**Dependencies to cite:**

- **GetOrganelle**: Jin et al. (2020) *Genome Biology*
- **PGA**: Qu et al. (2019) *Bioinformatics*
- **MAFFT**: Katoh & Standley (2013) *Molecular Biology and Evolution*
- **MACSE**: Ranwez et al. (2011) *PLoS ONE*
- **IQ-TREE**: Nguyen et al. (2015) *Molecular Biology and Evolution*
- **Biopython**: Cock et al. (2009) *Bioinformatics*


---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

CGAS integrates and builds upon several excellent tools:

- **GetOrganelle** - Chloroplast assembly from NGS data
- **PGA** - Plastid genome annotation
- **MAFFT** - Multiple sequence alignment
- **MACSE** - Codon-aware alignment for coding sequences
- **IQ-TREE** - Maximum likelihood phylogenetic inference
- **Biopython** - Python tools for biological computation
- **R** - Statistical computing and visualization

Special thanks to all developers and the bioinformatics community.

---

## 📧 Contact & Support

- **Issues:** [GitHub Issues](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS/issues)
- **Repository:** [GitHub Repository](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS)
- **Author:** Abdullah

---

## 📌 Version History

### v1.0.1 (2026)
- ✅ Initial public release
- ✅ 14 integrated modules
- ✅ Complete conda support
- ✅ MACSE integration for phylogeny
- ✅ Comprehensive R visualization
- ✅ Outgroup handling mechanism
- ✅ Python 3.9+ support
- ✅ Jupyter notebook compatibility

---

<p align="center">
  <strong>Made with ❤️ for chloroplast genomics research</strong>
</p>

<p align="center">
  <a href="#-quick-start">Quick Start</a> •
  <a href="#-modules-overview">Modules</a> •
  <a href="#-usage-guide">Usage</a> •
  <a href="#-complete-workflow-example">Workflow</a> •
  <a href="#-citation">Citation</a>
</p>
