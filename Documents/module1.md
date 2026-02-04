d

# CGAS Module 1: Chloroplast Genome Assembly and Quality Control Pipeline
## Complete Documentation and User Guide

---

## 📋 Table of Contents
1. [Introduction](#introduction)
2. [Quick Start - Just Run It!](#quick-start---just-run-it)
3. [Installation Guide](#installation-guide)
4. [Command-Line Usage with Examples](#command-line-usage-with-examples)
6. [Input Requirements](#input-requirements)
7. [Output Structure](#output-structure)
8. [Detailed Feature Explanation](#detailed-feature-explanation)
9. [Assembly Quality Assessment](#assembly-quality-assessment)
10. [Read Mapping and Coverage Analysis](#read-mapping-and-coverage-analysis)
11. [Troubleshooting](#troubleshooting)
12. [Examples](#examples)
13. [FAQ](#faq)
14. [Technical Specifications](#technical-specifications)
15. [References](#references)

---

## Introduction

**CGAS Module 1** is a comprehensive pipeline designed for chloroplast genome assembly and quality control. This module processes raw sequencing reads through quality control, assembly with GetOrganelle, validation, read mapping, and coverage analysis to produce high-quality chloroplast genome assemblies ready for downstream analysis.

This module performs automated chloroplast genome assembly with:
- **Fastp integration**: Advanced read quality control and trimming
- **GetOrganelle assembly**: Specialized chloroplast genome assembly
- **Assembly validation**: Completeness assessment and circularity detection
- **Read mapping**: BWA-based mapping to assembled genomes
- **Coverage analysis**: Comprehensive coverage statistics
- **Quality reporting**: Detailed assembly and QC metrics

### Key Features:
- **Complete pipeline**: From raw reads to validated assemblies
- **Quality control**: Advanced read filtering with Fastp
- **Specialized assembly**: GetOrganelle optimized for organellar genomes
- **Automatic validation**: Detects complete, circular, and flip-flop assemblies
- **Coverage analysis**: BWA mapping and depth calculation
- **Comprehensive reporting**: Statistics for each step of the pipeline
- **Batch processing**: Handles multiple samples simultaneously

### Scientific Applications:
- **De novo assembly**: Assemble chloroplast genomes from raw reads
- **Quality control**: Assess assembly completeness and quality
- **Coverage analysis**: Evaluate sequencing depth and uniformity
- **Structure detection**: Identify LSC, SSC, and IR regions
- **Comparative genomics**: Prepare standardized assemblies for comparison
- **Large-scale projects**: Batch processing of multiple samples

### Part of CGAS (Chloroplast Genome Analysis Suite)
CGAS is a comprehensive suite of tools designed for chloroplast genome analysis, providing seamless integration across different analysis modules.

---

## Quick Start - Just Run It!

### ⚡ One Command Does It All!

```bash
# Option 1: Using cgas command
cgas --module 1 -i raw_reads/ -o output/

# Option 2: Using cgas-assembly shortcut
cgas-assembly -i raw_reads/ -o output/

# Option 3: Using python directly
python cgas_module1.py -i raw_reads/ -o output/
```

**What happens when you run this:**
1. ✅ Finds all read pairs in input directory
2. ✅ Runs Fastp quality control on all samples
3. ✅ Assembles chloroplast genomes with GetOrganelle
4. ✅ Validates assemblies and detects completeness
5. ✅ Maps reads back to assemblies with BWA
6. ✅ Calculates coverage statistics
7. ✅ Generates comprehensive reports
8. ✅ Copies complete assemblies to organized folder

### 📁 Real-World Folder Examples

#### Example 1: Basic Usage (Easiest!)
```bash
# Your folder structure:
/home/abdullah/assembly_project/
├── raw_reads/
│   ├── SRR8666784_1.fastq.gz
│   ├── SRR8666784_2.fastq.gz
│   ├── SRR8666785_1.fastq.gz
│   ├── SRR8666785_2.fastq.gz
│   └── SRR8666786_1.fastq.gz
│   └── SRR8666786_2.fastq.gz

# Navigate to project folder
cd /home/abdullah/assembly_project/

# Run assembly pipeline
python cgas_module1.py -i raw_reads/ -o results/

# Output created automatically:
# results/
# ├── 01_QC/                           # Fastp reports
# ├── 02_clean_reads/                  # Cleaned reads
# ├── 03_assemblies/                   # GetOrganelle outputs
# ├── 04_mapping/                      # BWA mapping results
# ├── 05_cp_reads/                     # Chloroplast reads
# ├── 06_reports/                      # Summary reports
# └── 07_assembled_genomes/            # Complete assemblies only
```

#### Example 2: Custom Parameters
```bash
# Use custom parameters for difficult samples
python cgas_module1.py -i raw_reads/ -o results/ -t 16 -k "21,33,55,77,99,121" -R 20
cgas --module 1 i raw_reads/ -o results/ -t 16 -k "21,33,55,77,99,121" -R 20
```

#### Example 3: Force Rerun
```bash
# Force processing of all samples (skip existing check)
python cgas_module1.py -i raw_reads/ -o results/ --no-skip-existing
#skipping the already completed assembly and quality checking
python cgas_module1.py -i raw_reads/ -o output/ --skip-existing
```

### 🎯 Quick Start Cheat Sheet

```bash
# ====== cgas command ======
cgas --module 1 -i raw_reads/ -o output/                    # Basic
cgas --module 1 -i raw_reads/ -o output/ -t 16              # More threads
cgas --module 1 -i raw_reads/ -o output/ -k "21,33,55,77,99"  # Custom k-mers
cgas --module 1 -i raw_reads/ -o output/ -R 20              # More rounds
cgas --module 1 -i raw_reads/ -o output/ --force-mapping    # Force mapping
cgas --module 1 -i raw_reads/ -o output/ --skip-existing    # Skip done samples

# ====== cgas-assembly shortcut ======
cgas-assembly -i raw_reads/ -o output/                      # Basic
cgas-assembly -i raw_reads/ -o output/ -t 16 -R 15 --trim-poly-g

# ====== python command ======
python cgas_module1.py -i raw_reads/ -o output/             # Basic
python cgas_module1.py -i raw_reads/ -o output/ -t 16       # More threads
python cgas_module1.py -i raw_reads/ -o output/ -k "21,33,55,77,99"  # Custom k-mers
python cgas_module1.py -i raw_reads/ -o output/ -R 20       # More rounds
python cgas_module1.py -i raw_reads/ -o output/ --force-mapping  # Force mapping
python cgas_module1.py --help                               # Get help
```

### 📊 What You Get (Output Files)

```
results/                                    # Created automatically
├── 📁 Quality Control Reports
│   ├── 01_QC/                             # Fastp QC reports
│   │   ├── SRR8666784/
│   │   │   ├── SRR8666784_fastp.html      # HTML QC report
│   │   │   └── SRR8666784_fastp.json      # JSON QC metrics
│   │   ├── SRR8666785/
│   │   └── SRR8666786/
│
├── 📁02_clean_reads/                    # Cleaned FASTQ files
│   │   ├── SRR8666784_clean_R1.fq.gz      # Cleaned forward reads
│   │   ├── SRR8666784_clean_R2.fq.gz      # Cleaned reverse reads
│   │   ├── SRR8666785_clean_R1.fq.gz
│   │   ├── SRR8666785_clean_R2.fq.gz
│   │   ├── SRR8666786_clean_R1.fq.gz
│   │   └── SRR8666786_clean_R2.fq.gz
│
├── 📁03_assemblies/                     # GetOrganelle outputs
│   │   ├── SRR8666784_assembly/
│   │   │   ├── SRR8666784.path_sequence.fasta  # Assembly FASTA
│   │   │   ├── get_org.log.txt            # GetOrganelle log
│   │   │   └── ...                        # Other GetOrganelle files
│   │   ├── SRR8666785_assembly/
│   │   └── SRR8666786_assembly/
│
├── 📁 04_mapping/                        # BWA mapping results
│   │   ├── SRR8666784/
│   │   │   ├── SRR8666784.sorted.bam      # Sorted BAM file
│   │   │   └── SRR8666784.sorted.bam.bai  # BAM index
│   │   ├── SRR8666785/
│   │   └── SRR8666786/
│
├── 📁 05_cp_reads/                       # Extracted chloroplast reads
│   │   ├── SRR8666784/
│   │   │   ├── SRR8666784_cp_R1.fq.gz     # Mapped forward reads
│   │   │   └── SRR8666784_cp_R2.fq.gz     # Mapped reverse reads
│   │   ├── SRR8666785/
│   │   └── SRR8666786/
│
├── 📁06_reports/                        # Summary reports
│   │   ├── 00_SUMMARY_ALL_SAMPLES.tsv     # Tab-delimited summary
│   │   ├── 00_SUMMARY_ALL_SAMPLES.xlsx    # Excel summary
│   │   ├── SRR8666784_report.txt          # Individual sample report
│   │   ├── SRR8666785_report.txt
│   │   ├── SRR8666786_report.txt
│   │   └── 00_PROBLEMATIC_ASSEMBLIES.txt  # Problematic assemblies list
│
└── 📁07_assembled_genomes/              # Complete assemblies only
        ├── SRR8666784_1.fasta             # Complete assembly 1
        ├── SRR8666785_1.fasta             # Complete assembly 2
        └── SRR8666786_1.fasta             # Complete assembly 3
```

---

## Installation Guide

> **Note:** If you have already installed the full CGAS tool (using `environment.yml`), all dependencies for Module 1 are already included — no additional installation is needed. The steps below are only for running Module 1 as a standalone script. However, If you want to run just module 1 then the conda enviroment or your system should has the following tools. Any change to link is possible. Please consider visiting webiste of each tools for complete installations guide. 
#These are external tools, please site them along CGAS to give full credit to orignal researchers. 

### Prerequisites
- **Python 3.9 or higher** (tested on Python 3.9–3.12)
- **Fastp** for read quality control
- **GetOrganelle** for chloroplast genome assembly
- **BWA** for read mapping
- **Samtools** for BAM file processing

### Step-by-Step Installation

#### 1. Install Python Dependencies
```bash
# Install required packages
pip install biopython pandas openpyxl

# Verify installation
python -c "import Bio, pandas, openpyxl; print('All packages installed successfully')"
```

**Recommended versions:**
```bash
pip install biopython>=1.79 pandas>=1.3.0 openpyxl>=3.0.9
```

#### 2. Install Fastp

**Ubuntu/Debian:**
```bash
# Download fastp
wget http://github.com/OpenGene/fastp/releases/download/v0.23.2/fastp
chmod a+x ./fastp
sudo mv fastp /usr/local/bin/

# Verify installation
fastp --version
```

**macOS:**
```bash
# Install with Homebrew
brew install fastp

# Verify installation
fastp --version
```

**From source (all platforms):**
```bash
# Clone and compile
git clone https://github.com/OpenGene/fastp.git
cd fastp
make
sudo make install

# Verify installation
fastp --version
```

#### 3. Install GetOrganelle

```bash
# Install GetOrganelle
pip install getorganelle

# Or install with conda
conda install -c bioconda getorganelle

# Verify installation
get_organelle_from_reads.py --help
```

**Alternative installation:**
```bash
# Clone from GitHub
git clone https://github.com/Kinggerm/GetOrganelle.git
cd GetOrganelle
pip install -e .
```

#### 4. Install BWA

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install bwa
```

**macOS:**
```bash
brew install bwa
```

**CentOS/RHEL:**
```bash
sudo yum install bwa
```

**From source (all platforms):**
```bash
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download
tar -xvjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
sudo cp bwa /usr/local/bin/
```

#### 5. Install Samtools

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install samtools
```

**macOS:**
```bash
brew install samtools
```

**CentOS/RHEL:**
```bash
sudo yum install samtools
```

**From source (all platforms):**
```bash
wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2
tar -xvjf samtools-1.15.tar.bz2
cd samtools-1.15
./configure
make
sudo make install
```

#### 7. Verify Complete Installation
```bash
# Test Python
python --version                    # Should be 3.9+

# Test packages
python -c "import Bio, pandas, openpyxl"  # Should not give errors

# Test tools
fastp --version                     # Should show fastp version
get_organelle_from_reads.py --help  # Should show GetOrganelle help
bwa                                 # Should show BWA help
samtools --version                  # Should show samtools version
```

### Dependency Details

| Package | Version | Purpose | Required |
|---------|---------|---------|----------|
| **Python** | ≥3.9 | Script execution | Yes |
| **Biopython** | ≥1.79 | FASTA file parsing | Yes |
| **pandas** | ≥1.3.0 | Data manipulation | Yes |
| **openpyxl** | ≥3.0.9 | Excel file generation | Yes |
| **Fastp** | ≥0.20.0 | Read quality control | Yes |
| **GetOrganelle** | ≥1.7.0 | Chloroplast assembly | Yes |
| **BWA** | ≥0.7.17 | Read mapping | Yes |
| **Samtools** | ≥1.15 | BAM processing | Yes |

---

## Command-Line Usage with Examples

### Complete Command Reference

```bash
# ====================================================================
# USING cgas COMMAND
# ====================================================================

# Basic assembly
cgas --module 1 -i /home/abdullah/raw_reads/ -o /home/abdullah/results/

# More threads for faster processing
cgas --module 1 -i raw_reads/ -o results/ -t 16

# Skip already processed samples when adding new data
cgas --module 1 -i raw_reads/ -o results/ --skip-existing

# Custom k-mer values for difficult samples
cgas --module 1 -i raw_reads/ -o results/ -k "21,33,55,77,99"

# More GetOrganelle rounds for challenging assemblies
cgas --module 1 -i raw_reads/ -o results/ -R 20

# Force mapping of incomplete assemblies
cgas --module 1 -i raw_reads/ -o results/ --force-mapping

# Use nohup for large datasets on a server
cgas --module 1 -i large_dataset/ -o assemblies/ --use-nohup
```

```bash
# ====================================================================
# USING cgas-assembly SHORTCUT
# ====================================================================

# Basic assembly
cgas-assembly -i raw_reads/ -o results/

# With all options
cgas-assembly \
  -i raw_reads/ \
  -o assemblies/ \
  -t 16 \
  -F embplant_pt \
  -k 21,45,65,85,105 \
  -R 15 \
  --trim-poly-g
```

```bash
# ====================================================================
# USING python COMMAND
# ====================================================================

# ====================================================================
# BASIC USAGE
# ====================================================================

# 1. Basic assembly (required parameters)
python cgas_module1.py -i raw_reads/ -o output/

# 2. Custom number of threads
python cgas_module1.py -i raw_reads/ -o output/ -t 16

# 3. Custom k-mer values
python cgas_module1.py -i raw_reads/ -o output/ -k "21,33,55,77,99"

# 4. More GetOrganelle rounds
python cgas_module1.py -i raw_reads/ -o output/ -R 20

# 5. Force mapping of incomplete assemblies
python cgas_module1.py -i raw_reads/ -o output/ --force-mapping

# 6. Get help
python cgas_module1.py --help


# ====================================================================
# REAL-WORLD EXAMPLES
# ====================================================================

# Example 1: High-performance assembly
python cgas_module1.py -i hibiscus_reads/ -o hibiscus_assemblies/ -t 24 -k "21,33,55,77,99,121" -R 25

# Example 2: Low-quality data processing
python cgas_module1.py -i difficult_samples/ -o results/ --trim-poly-g --cut-mean-quality 15 -R 30

# Example 3: Custom tool paths
python cgas_module1.py -i reads/ -o assemblies/ --fastp /custom/path/fastp --getorganelle /custom/path/get_organelle_from_reads.py


# ====================================================================
# ADVANCED PARAMETERS
# ====================================================================

# Example 1: NovaSeq data with poly-G trimming
python cgas_module1.py -i novaseq_reads/ -o assemblies/ --trim-poly-g

# Example 2: Background processing with nohup
python cgas_module1.py -i large_dataset/ -o assemblies/ --use-nohup

# Example 3: Skip existing samples
python cgas_module1.py -i reads/ -o assemblies/ --skip-existing

# Example 4: Force rerun all samples
python cgas_module1.py -i reads/ -o assemblies/ --no-skip-existing

# Example 5: Different genome type
python cgas_module1.py -i reads/ -o assemblies/ --genome-type embplant_mt
```

### Parameter Details

| Parameter | Short | Required | Default | Description |
|-----------|-------|----------|---------|-------------|
| `--input` | `-i` | Yes | - | Directory with raw FASTQ files |
| `--output` | `-o` | Yes | - | Output directory |
| `--threads` | `-t` | No | 8 | Number of threads for parallel processing |
| `--fastp` | - | No | `fastp` | Path to fastp executable |
| `--getorganelle` | - | No | `get_organelle_from_reads.py` | Path to GetOrganelle script |
| `--bwa` | - | No | `bwa` | Path to BWA executable |
| `--samtools` | - | No | `samtools` | Path to samtools executable |
| `--genome-type` | - | No | `embplant_pt` | Genome type for GetOrganelle |
| `--k-values` | `-k` | No | `21,45,65,85,105` | K-mer values for assembly |
| `--rounds` | `-R` | No | 15 | Number of GetOrganelle rounds |
| `--use-nohup` | - | No | False | Use nohup for long-running assemblies |
| `--log-level` | - | No | INFO | Logging level |
| `--force-mapping` | - | No | False | Force mapping even if assembly is incomplete |
| `--trim-poly-g` | - | No | False | Enable poly-G tail trimming |
| `--cut-mean-quality` | - | No | 20 | Mean quality requirement for sliding window |
| `--skip-existing` | - | No | True | Skip samples that have already been processed |

---

---

## Input Requirements

### Supported File Formats

**Input FASTQ files:**
- `.fastq`
- `.fastq.gz`
- `.fq`
- `.fq.gz`

**Read naming conventions supported:**
- `sample_1.fastq.gz` / `sample_2.fastq.gz`
- `sample_R1.fastq.gz` / `sample_R2.fastq.gz`
- `sample_1.clean.fq.gz` / `sample_2.clean.fq.gz`

### Critical Requirements

**1. Raw Sequencing Reads**
```bash
# Example of paired-end reads
sample_1.fastq.gz
sample_2.fastq.gz

# Or with different naming
sample_R1.fastq.gz
sample_R2.fastq.gz
```

**Requirements:**
- Illumina sequencing data recommended
- Paired-end reads preferred (150bp or longer)
- Minimum 30x chloroplast genome coverage recommended
- Raw or previously QC'd reads (both supported)

**2. Single-end Reads (Supported)**
```bash
# Single-end reads #If place in the main directory of reads; this will still proceed along pairend reads without any issue
sample.fastq.gz
```

**Note:** Single-end assemblies may require more rounds and higher coverage for successful assembly.

### File Organization

```bash
project_folder/
└── raw_reads/                      # Input FASTQ files
    ├── SRR8666784_1.fastq.gz
    ├── SRR8666784_2.fastq.gz
    ├── SRR8666785_1.fastq.gz
    ├── SRR8666785_2.fastq.gz
    ├── SRR8666786_1.fastq.gz
    └── SRR8666786_2.fastq.gz
```

### Read Quality Requirements
#Please check getorganelle manual for details
**Recommended quality metrics:**
- **Q30 rate**: >85% (before QC)
- **Read length**: 150bp or longer
- **Insert size**: 300-500bp for paired-end
- **Coverage**: 30x chloroplast genome or higher
- **Contamination**: Minimal nuclear/mitochondrial DNA

**Challenging samples may benefit from:**
- Increasing GetOrganelle rounds (`-R 20` or higher)
- Adjusting k-mer values (`-k "21,33,55,77,99,121"`)
- Enabling poly-G trimming (`--trim-poly-g`)
- Lowering quality cutoff (`--cut-mean-quality 15`)

---

## Output Structure

### Directory Organization

```
results/                                    # Created automatically
├── 📁 Quality Control Reports
│   ├── 01_QC/                             # Fastp QC reports
│   │   ├── SRR8666784/
│   │   │   ├── SRR8666784_fastp.html      # Interactive HTML report
│   │   │   └── SRR8666784_fastp.json      # Machine-readable metrics
│   │   │       └── Content:
│   │   │           {
│   │   │               "summary": {
│   │   │                   "before_filtering": {
│   │   │                       "total_reads": 5000000,
│   │   │                       "total_bases": 750000000,
│   │   │                       "q20_rate": 0.95,
│   │   │                       "q30_rate": 0.90
│   │   │                   },
│   │   │                   "after_filtering": {
│   │   │                       "total_reads": 4800000,
│   │   │                       "total_bases": 720000000,
│   │   │                       "q20_rate": 0.97,
│   │   │                       "q30_rate": 0.93
│   │   │                   },
│   │   │                   "duplication": {"rate": 0.15}
│   │   │               }
│   │   │           }
│   │   ├── SRR8666785/
│   │   └── SRR8666786/
│
├── 📁 Clean Reads
│   ├── 02_clean_reads/                    # Cleaned FASTQ files
│   │   ├── SRR8666784_clean_R1.fq.gz      # Cleaned forward reads
│   │   ├── SRR8666784_clean_R2.fq.gz      # Cleaned reverse reads
│   │   ├── SRR8666785_clean_R1.fq.gz
│   │   ├── SRR8666785_clean_R2.fq.gz
│   │   ├── SRR8666786_clean_R1.fq.gz
│   │   └── SRR8666786_clean_R2.fq.gz
│
├── 📁 Assemblies
│   ├── 03_assemblies/                     # GetOrganelle outputs
│   │   ├── SRR8666784_assembly/
│   │   │   ├── SRR8666784.path_sequence.fasta  # Assembly FASTA
│   │   │   │   └── Content:
│   │   │   │       >SRR8666784_1 circular=true
│   │   │   │       ATGGCGACGACGTTCGTCGTCGTTTGTCGATCTCGTCTGACTTCAGCCTGATCGGTAGCA
│   │   │   │       ...
│   │   │   ├── get_org.log.txt            # GetOrganelle log
│   │   │   │   └── Content:
│   │   │   │       Total:LSC:SSC:Repeat(bp) = 160373:89521:20300:25276
│   │   │   │       ...
│   │   │   └── ...                        # Other GetOrganelle files
│   │   ├── SRR8666785_assembly/
│   │   └── SRR8666786_assembly/
│
├── 📁 Mapping Results
│   ├── 04_mapping/                        # BWA mapping results
│   │   ├── SRR8666784/
│   │   │   ├── SRR8666784.sorted.bam      # Sorted BAM file
│   │   │   └── SRR8666784.sorted.bam.bai  # BAM index
│   │   ├── SRR8666785/
│   │   └── SRR8666786/
│
├── 📁 Chloroplast Reads
│   ├── 05_cp_reads/                       # Extracted chloroplast reads
│   │   ├── SRR8666784/
│   │   │   ├── SRR8666784_cp_R1.fq.gz     # Mapped forward reads
│   │   │   └── SRR8666784_cp_R2.fq.gz     # Mapped reverse reads
│   │   ├── SRR8666785/
│   │   └── SRR8666786/
│
├── 📁 Reports
│   ├── 06_reports/                        # Summary reports
│   │   ├── 00_SUMMARY_ALL_SAMPLES.tsv     # Tab-delimited summary
│   │   │   └── Content:
│   │   │       Sample	Reads_Before_QC	Reads_After_QC	Bases_After_QC_Gb	Q20_Rate_%	Q30_Rate_%	Num_Contigs	Assembly_Status	Circular	Mapped_Reads	Mapping_Rate_%	Mean_Coverage	Median_Coverage
│   │   │       SRR8666784_1	5000000	4800000	0.72	97.00	93.00	1	COMPLETE	Yes	4500000	93.75	850.5	842.0
│   │   │       SRR8666785_1	4500000	4300000	0.65	96.50	92.00	2	COMPLETE	Yes	4000000	93.02	780.2	775.5
│   │   │       SRR8666786_1	5200000	5100000	0.77	97.50	94.00	1	COMPLETE	Yes	4900000	96.08	920.8	915.0
│   │   ├── 00_SUMMARY_ALL_SAMPLES.xlsx    # Excel summary with formatting
│   │   ├── SRR8666784_report.txt          # Individual sample report
│   │   │   └── Content:
│   │   │       ================================================================================
│   │   │       CGAS Module 1 Report: SRR8666784
│   │   │       ================================================================================
│   │   │       
│   │   │       1. QUALITY CONTROL (Fastp)
│   │   │       --------------------------------------------------------------------------------
│   │   │       Total reads before QC:  5,000,000
│   │   │       Total reads after QC:   4,800,000
│   │   │       Total bases before QC:  0.75 Gb
│   │   │       Total bases after QC:   0.72 Gb
│   │   │       Q20 rate before QC:     95.00%
│   │   │       Q20 rate after QC:      97.00%
│   │   │       Q30 rate before QC:     90.00%
│   │   │       Q30 rate after QC:      93.00%
│   │   │       GC content:             37.50%
│   │   │       Duplication rate:       15.00%
│   │   │       
│   │   │       2. CHLOROPLAST GENOME ASSEMBLY (GetOrganelle)
│   │   │       --------------------------------------------------------------------------------
│   │   │       Number of contigs:      1
│   │   │       Assembly status:        COMPLETE
│   │   │       Circular:               Yes
│   │   │       Total length:           160,373 bp
│   │   │       LSC length:             89,521 bp
│   │   │       SSC length:             20,300 bp
│   │   │       IR length:              25,276 bp
│   │   │       
│   │   │       3. READ MAPPING AND COVERAGE
│   │   │       --------------------------------------------------------------------------------
│   │   │       Total reads:            4,800,000
│   │   │       Mapped reads:           4,500,000
│   │   │       Mapping rate:           93.75%
│   │   │       Mean coverage:          850.5x
│   │   │       Median coverage:        842.0x
│   │   │       CP reads saved:         4,500,000
│   │   │       
│   │   │       4. RECOMMENDATIONS
│   │   │       --------------------------------------------------------------------------------
│   │   │       ✓ Assembly is COMPLETE and ready for downstream CGAS analyses.
│   │   │       ✓ Coverage is excellent (>800x)
│   │   │       
│   │   │       ================================================================================
│   │   ├── SRR8666785_report.txt
│   │   ├── SRR8666786_report.txt
│   │   └── 00_PROBLEMATIC_ASSEMBLIES.txt  # Problematic assemblies list
│   │       └── Content:
│   │           ================================================================================
│   │           PROBLEMATIC ASSEMBLIES REPORT
│   │           ================================================================================
│   │           
│   │           Total problematic assemblies: 0
│   │           
│   │           These samples were NOT copied to 07_assembled_genomes for the following reasons:
│   │           
│   │           ISSUE TYPES:
│   │             1. >2 contigs: More than 2 path_sequence.fasta files
│   │                (Expected: 1 circular or 2 for SSC flip-flop)
│   │           
│   │             2. Fragmented assembly: Multiple sequences (headers) in single FASTA file
│   │                (Expected: 1 header per path_sequence.fasta for complete genome)
│   │           
│   │           --------------------------------------------------------------------------------
│   │           
│   │           No problematic assemblies found.
│
└── 📁 Complete Assemblies
    └── 07_assembled_genomes/              # Complete assemblies only
        ├── SRR8666784_1.fasta             # Complete assembly 1
        │   └── Content:
        │       >SRR8666784_1
        │       ATGGCGACGACGTTCGTCGTCGTTTGTCGATCTCGTCTGACTTCAGCCTGATCGGTAGCA
        │       ...
        ├── SRR8666785_1.fasta             # Complete assembly 1 of 2
        ├── SRR8666785_2.fasta             # Complete assembly 2 of 2 (SSC flip-flop)
        └── SRR8666786_1.fasta             # Complete assembly 1
```

### Key Output Files Explained

#### 1. Quality Control Reports (01_QC)

**Fastp HTML report:**
- Interactive visualization of QC metrics
- Before/after quality plots
- Adapter content analysis
- Duplication levels
- K-mer content

**Fastp JSON report:**
- Machine-readable metrics
- Used for downstream analysis
- Contains detailed statistics

#### 2. Clean Reads (02_clean_reads)

**Quality-filtered reads:**
- Adapter sequences removed
- Low-quality bases trimmed
- Poly-G tails trimmed (if enabled)
- Ready for assembly

#### 3. Assemblies (03_assemblies)

**GetOrganelle output:**
- Assembly FASTA files
- Assembly logs with LSC/SSC/IR information
- Intermediate files for debugging

#### 4. Mapping Results (04_mapping)

**BWA alignment files:**
- Sorted BAM files
- BAM index files
- Used for coverage analysis

#### 5. Chloroplast Reads (05_cp_reads)

**Extracted chloroplast reads:**
- Reads that mapped to chloroplast genome
- Useful for downstream analyses
- Can be used for re-assembly if needed

#### 6. Reports (06_reports)

**Summary reports:**
- Tab-delimited and Excel formats
- Individual sample reports
- Problematic assemblies list

#### 7. Complete Assemblies (07_assembled_genomes)

**Validated complete assemblies:**
- Only complete assemblies are copied here
- Clean naming convention
- Ready for downstream analysis (Module 2)

---

## Detailed Feature Explanation

### 1. Read Quality Control with Fastp

**Fastp integration:**
```python
# Fastp command construction
cmd = [
    self.fastp_path,
    "-i", str(read_pair.forward_read),     # Input forward reads
    "-o", str(clean_r1),                   # Output forward reads
    "-j", str(json_report),                # JSON report
    "-h", str(html_report),                # HTML report
    "-w", str(self.threads),               # Number of threads
    "--cut_mean_quality", str(self.cut_mean_quality),  # Quality cutoff
]

# Add poly-G trimming for NovaSeq data
if self.trim_poly_g:
    cmd.append("--trim_poly_g")

# Add paired-end specific options
if read_pair.is_paired:
    cmd.append("--detect_adapter_for_pe")
    cmd.extend(["-I", str(read_pair.reverse_read), "-O", str(clean_r2)])
```

**Fastp parameters:**
- **Adapter detection**: Automatic detection for paired-end reads
- **Quality trimming**: Sliding window based quality cutting
- **Poly-G trimming**: For NovaSeq two-color chemistry
- **Length filtering**: Automatically filters very short reads
- **Duplication analysis**: Identifies duplicate reads

### 2. Read Pair Detection

**Automatic read pair detection:**
```python
# Supported read file patterns
FORWARD_PATTERNS = [
    r'(.+?)_1\.fastq$',
    r'(.+?)_1\.fastq\.gz$',
    r'(.+?)_1\.fq$',
    r'(.+?)_1\.fq\.gz$',
    r'(.+?)_1\.clean\.fq\.gz$',
    r'(.+?)_R1\.fastq$',
    r'(.+?)_R1\.fastq\.gz$',
    r'(.+?)_R1\.fq$',
    r'(.+?)_R1\.fq\.gz$',
]

def _find_reverse_read(self, forward_path: Path) -> Optional[Path]:
    """Find corresponding reverse read for a forward read"""
    forward_name = forward_path.name
    
    # Define replacement patterns
    replacements = [
        ('_1.fastq', '_2.fastq'),
        ('_1.fastq.gz', '_2.fastq.gz'),
        ('_1.fq', '_2.fq'),
        ('_1.fq.gz', '_2.fq.gz'),
        ('_1.clean.fq.gz', '_2.clean.fq.gz'),
        ('_R1.fastq', '_R2.fastq'),
        ('_R1.fastq.gz', '_R2.fastq.gz'),
        ('_R1.fq', '_R2.fq'),
        ('_R1.fq.gz', '_R2.fq.gz'),
    ]
```

**Read pair detection logic:**
1. Scan input directory for forward read patterns
2. For each forward read, find corresponding reverse read
3. Handle various naming conventions
4. Support single-end reads

### 3. GetOrganelle Assembly

**GetOrganelle integration:**
```python
# Build GetOrganelle command
cmd = [
    self.getorganelle_path,
    "-o", str(assembly_out),              # Output directory
    "-F", self.genome_type,                # Genome type
    "-R", str(self.rounds),                # Number of rounds
    "-k", self.k_values,                   # K-mer values
    "-t", str(self.threads),               # Number of threads
    "--overwrite"                          # Overwrite existing
]

# Add read files based on paired or single-end
if is_paired:
    cmd.extend(["-1", str(clean_r1), "-2", str(clean_r2)])
else:
    cmd.extend(["-u", str(clean_r1)])
```

**GetOrganelle parameters:**
- **Genome type**: embplant_pt (embryophyte plant chloroplast)
- **K-mer values**: Default "21,45,65,85,105"
- **Rounds**: Default 15 iterations
- **Working mode**: Paired-end or single-end

### 4. Assembly Validation

**Completeness assessment:**
```python
def _parse_getorganelle_output(self, assembly_dir: Path, sample_name: str) -> AssemblyStats:
    """Parse GetOrganelle output to assess assembly completeness"""
    
    # Look for assembly files - prioritize path_sequence.fasta
    fasta_files = list(assembly_dir.glob("*.path_sequence.fasta"))
    
    # Parse FASTA file
    contigs = self._parse_fasta(assembly_file)
    num_contigs = len(contigs)
    
    # Check for circularity in FASTA headers
    is_circular = any('circular' in header.lower() for header, _ in contigs)
    
    # Parse GetOrganelle log for LSC/SSC/IR information
    log_file = assembly_dir / "get_org.log.txt"
    if log_file.exists():
        lsc_length, ssc_length, ir_length = self._parse_getorganelle_log(log_file)
    
    # Enhanced SSC flip-flop detection
    has_flipflop = False
    if num_contigs == 2 and lsc_length and ssc_length and ir_length:
        # Check if contig sizes are consistent with flip-flop
        expected_size = lsc_length + ssc_length + (2 * ir_length)
        contig_sizes = sorted([len(seq) for _, seq in contigs])
        
        # Both contigs should be approximately the same size
        size_diff = abs(contig_sizes[0] - contig_sizes[1])
        avg_size = sum(contig_sizes) / 2
        
        if size_diff < 100 and abs(avg_size - expected_size) < 500:
            has_flipflop = True
    
    # Determine completeness
    is_complete = False
    if has_flipflop:
        is_complete = True
    elif num_contigs == 1 and is_circular:
        is_complete = True
    elif num_contigs == 1 and lsc_length and ssc_length and ir_length:
        # Single contig with structure info
        expected_size = lsc_length + ssc_length + (2 * ir_length)
        if abs(total_length - expected_size) < 500:
            is_complete = True
    
    return AssemblyStats(...)
```

**Assembly validation logic:**
1. Count number of contigs
2. Check for circularity in headers
3. Parse LSC/SSC/IR information from log
4. Detect SSC flip-flop (2 contigs)
5. Determine assembly completeness

---

## Assembly Quality Assessment

### Completeness Criteria

**Complete assembly indicators:**
- **Single circular contig**: 1 contig marked as circular
- **SSC flip-flop**: 2 contigs with appropriate sizes
- **Structure match**: Size matches expected LSC+SSC+2*IR

**Incomplete assembly indicators:**
- **Multiple contigs**: >2 contigs
- **Non-circular**: Single contig not marked circular
- **Size mismatch**: Doesn't match expected structure

### Structure Detection

**LSC/SSC/IR parsing:**
```python
def _parse_getorganelle_log(self, log_file: Path) -> Tuple[Optional[int], Optional[int], Optional[int]]:
    """Parse GetOrganelle log to extract LSC, SSC, IR lengths"""
    lsc_length = None
    ssc_length = None
    ir_length = None
    
    with open(log_file, 'r') as f:
        content = f.read()
        
        # Parse from "Detecting large repeats" line
        # Format: Total:LSC:SSC:Repeat(bp) = 160373:89521:20300:25276
        repeat_match = re.search(r'Total:LSC:SSC:Repeat\(bp\)\s*=\s*(\d+):(\d+):(\d+):(\d+)', content)
        if repeat_match:
            total_length = int(repeat_match.group(1))
            lsc_length = int(repeat_match.group(2))
            ssc_length = int(repeat_match.group(3))
            ir_length = int(repeat_match.group(4))
    
    return lsc_length, ssc_length, ir_length
```

**Structure detection logic:**
1. Parse GetOrganelle log for structure information
2. Extract LSC, SSC, and IR lengths
3. Validate assembly size against expected structure
4. Identify SSC flip-flop assemblies

### Problematic Assemblies

**Problematic assembly detection:**
```python
def _copy_assembled_genome(self, sample_name: str, assembly_dir: Path, stats: AssemblyStats):
    """Copy complete assembled genomes to 07_assembled_genomes with clean naming"""
    
    # Find path_sequence.fasta files
    fasta_files = sorted(assembly_dir.glob("*.path_sequence.fasta"))
    
    num_contigs = len(fasta_files)
    
    # Check if assembly has too many contigs
    if num_contigs > 2:
        self.logger.warning(f"⚠ {sample_name} has {num_contigs} contigs (expected ≤2 for complete assembly)")
        self.problematic_assemblies[sample_name] = f"{num_contigs} contigs (expected 1-2)"
        return
    
    # Check if file has multiple headers (fragmented assembly)
    header_count = 0
    with open(fasta_file, 'r') as f_in:
        for line in f_in:
            if line.startswith('>'):
                header_count += 1
    
    # Check if file has multiple headers (fragmented assembly)
    if header_count > 1:
        self.logger.warning(f"⚠ {sample_name}: {fasta_file.name} has {header_count} headers (fragmented assembly)")
        self.problematic_assemblies[sample_name] = f"Fragmented assembly ({header_count} sequences in {fasta_file.name})"
        return
```

**Problematic assembly types:**
1. **Too many contigs**: >2 path_sequence.fasta files
2. **Fragmented assembly**: Multiple sequences in single FASTA file
3. **Incomplete structure**: Missing LSC/SSC/IR components

---

## Read Mapping and Coverage Analysis

### BWA Mapping

**BWA integration:**
```python
# Index reference
subprocess.run([self.bwa_path, "index", str(assembly_fasta)], check=True)

# Map reads
with open(sam_file, 'w') as sam_out:
    if is_paired:
        subprocess.run(
            [self.bwa_path, "mem", "-t", str(self.threads), 
             str(assembly_fasta), str(clean_r1), str(clean_r2)],
            stdout=sam_out, stderr=subprocess.PIPE, text=True, check=True
        )
    else:
        subprocess.run(
            [self.bwa_path, "mem", "-t", str(self.threads), 
             str(assembly_fasta), str(clean_r1)],
            stdout=sam_out, stderr=subprocess.PIPE, text=True, check=True
        )

# Convert to BAM and sort
subprocess.run([self.samtools_path, "view", "-b", str(sam_file), "-o", str(bam_file)], check=True)
subprocess.run([self.samtools_path, "sort", str(bam_file), "-o", str(sorted_bam)], check=True)
subprocess.run([self.samtools_path, "index", str(sorted_bam)], check=True)
```

**Mapping process:**
1. Index reference assembly with BWA
2. Map clean reads to assembly with BWA mem
3. Convert SAM to BAM format
4. Sort and index BAM file
5. Extract mapped reads for downstream analysis

### Coverage Calculation

**Coverage analysis:**
```python
def _calculate_coverage(self, sample_name: str, bam_file: Path, assembly_stats: AssemblyStats) -> MappingStats:
    """Calculate coverage statistics from BAM file"""
    
    # Get mapping statistics
    flagstat_output = subprocess.run(
        [self.samtools_path, "flagstat", str(bam_file)],
        capture_output=True, text=True, check=True
    ).stdout
    
    # Parse flagstat output
    total_reads = 0
    mapped_reads = 0
    for line in flagstat_output.split('\n'):
        if 'in total' in line:
            total_reads = int(line.split()[0])
        elif 'mapped (' in line and 'primary' not in line:
            mapped_reads = int(line.split()[0])
    
    mapping_rate = (mapped_reads / total_reads * 100) if total_reads > 0 else 0.0
    
    # Calculate coverage depth
    depth_output = subprocess.run(
        [self.samtools_path, "depth", str(bam_file)],
        capture_output=True, text=True, check=True
    ).stdout
    
    depths = []
    for line in depth_output.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 3:
                depths.append(int(parts[2]))
    
    mean_coverage = sum(depths) / len(depths) if depths else 0.0
    median_coverage = sorted(depths)[len(depths)//2] if depths else 0.0
    
    return MappingStats(...)
```

**Coverage metrics:**
- **Total reads**: Number of reads processed
- **Mapped reads**: Number of reads that mapped to assembly
- **Mapping rate**: Percentage of reads that mapped
- **Mean coverage**: Average depth across assembly
- **Median coverage**: Median depth across assembly

### Chloroplast Read Extraction

**CP read extraction:**
```python
def _extract_cp_reads(self, sample_name: str, bam_file: Path, clean_r1: Path, clean_r2: Path, is_paired: bool):
    """Extract reads that mapped to chloroplast genome"""
    
    # Get mapped read names
    mapped_bam = cp_reads_sample_dir / "mapped.bam"
    subprocess.run([self.samtools_path, "view", "-b", "-F", "4", 
                   str(bam_file), "-o", str(mapped_bam)], check=True)
    
    # Convert to FASTQ
    cp_r1 = cp_reads_sample_dir / f"{sample_name}_cp_R1.fq.gz"
    
    if is_paired:
        cp_r2 = cp_reads_sample_dir / f"{sample_name}_cp_R2.fq.gz"
        subprocess.run([self.samtools_path, "fastq", "-1", str(cp_r1), 
                      "-2", str(cp_r2), "-0", "/dev/null", "-s", "/dev/null",
                      str(mapped_bam)], check=True)
    else:
        # For single-end: use stdout redirection
        with open(cp_r1, "wb") as out_fq:
            subprocess.run([self.samtools_path, "fastq", str(mapped_bam)],
                         stdout=out_fq, check=True)
```

**Read extraction process:**
1. Filter BAM for mapped reads only
2. Convert mapped reads to FASTQ format
3. Separate forward and reverse reads for paired-end
4. Save chloroplast reads for downstream analysis

---

## Troubleshooting

### Common Issues and Solutions

#### 1. Missing Dependencies
```bash
❌ ERROR: fastp not found. Please install fastp
```
**Solution:**
```bash
# Download and install fastp
wget http://github.com/OpenGene/fastp/releases/download/v0.23.2/fastp
chmod a+x ./fastp
sudo mv fastp /usr/local/bin/

# Verify installation
fastp --version
```

#### 2. GetOrganelle Not Found
```bash
❌ ERROR: get_organelle_from_reads.py not found or not working
```
**Solution:**
```bash
# Install GetOrganelle
pip install getorganelle

# Or install with conda
conda install -c bioconda getorganelle

# Verify installation
get_organelle_from_reads.py --help

# Or specify custom path
python cgas_module1.py --getorganelle /custom/path/get_organelle_from_reads.py -i reads/ -o output/
```

#### 3. No Read Pairs Found
```bash
❌ No read pairs found in /path/to/directory
```
**Solution:**
```bash
# Check file extensions
ls *.fastq.gz *.fq.gz

# Check file naming
# Forward reads should have _1 or _R1 in name
# Reverse reads should have _2 or _R2 in name

# Supported patterns:
# sample_1.fastq.gz / sample_2.fastq.gz
# sample_R1.fastq.gz / sample_R2.fastq.gz
# sample_1.clean.fq.gz / sample_2.clean.fq.gz
```

#### 4. Assembly Fails
```bash
⚠ WARNING: GetOrganelle failed for sample1
```
**Solution:**
```
Common GetOrganelle issues:
- Low chloroplast coverage (<30x)
- Poor read quality
- High contamination
- Distant reference genome

Try these solutions:
- Increase GetOrganelle rounds (-R 20 or higher)
- Adjust k-mer values (-k "21,33,55,77,99,121")
- Lower quality cutoff (--cut-mean-quality 15)
- Enable poly-G trimming (--trim-poly-g)
```

#### 5. Incomplete Assembly
```bash
⚠ WARNING: Assembly for sample1 is INCOMPLETE
```
**Solution:**
```
Incomplete assembly indicators:
- Multiple contigs (>2)
- Non-circular assembly
- Size mismatch

Try these solutions:
- Increase GetOrganelle rounds (-R 25)
- Adjust k-mer values
- Check read quality
- Verify sequencing depth
```

#### 6. Low Mapping Rate
```bash
⚠ WARNING: Low mapping rate: 45%
```
**Solution:**
```
Low mapping rate causes:
- Poor assembly quality
- High contamination
- Incorrect genome type

Try these solutions:
- Check assembly completeness
- Verify correct genome type
- Review Fastp QC reports
- Consider re-assembly with different parameters
```

### Quick Fixes Table

| Problem | Symptom | Quick Fix |
|---------|---------|-----------|
| Missing fastp | "fastp not found" | Install fastp |
| GetOrganelle missing | "get_organelle not found" | Install GetOrganelle |
| No read pairs | "No read pairs found" | Check file naming |
| Assembly fails | "GetOrganelle failed" | Increase rounds, adjust k-mers |
| Incomplete assembly | "INCOMPLETE assembly" | Increase rounds, check quality |
| Low mapping rate | "Low mapping rate" | Check assembly quality |
| BWA missing | "bwa not found" | Install BWA |
| Samtools missing | "samtools not found" | Install Samtools |

---

## Examples

### Example 1: Complete Assembly Workflow
```bash
# 1. Prepare project structure
mkdir -p /home/abdullah/hibiscus_assembly/
cd /home/abdullah/hibiscus_assembly/

# 2. Organize files
mkdir -p raw_reads
# Copy FASTQ files to raw_reads/

# 3. Run assembly pipeline
python cgas_module1.py -i raw_reads/ -o results/ -t 16

# 4. Check results
ls results/07_assembled_genomes/
open results/06_reports/00_SUMMARY_ALL_SAMPLES.xlsx
```

### Example 2: Difficult Sample Processing
```bash
# For challenging samples with low quality or coverage
python cgas_module1.py -i difficult_reads/ -o results/ \
    -t 24 \
    -k "21,33,55,77,99,121" \
    -R 25 \
    --trim-poly-g \
    --cut-mean-quality 15 \
    --force-mapping
```

### Example 3: Single-end Data Processing
```bash
# For single-end reads
python cgas_module1.py -i single_end_reads/ -o single_end_results/ \
    -R 30 \
    -k "21,45,65,85,105,125"
```

### Example 4: Background Processing
```bash
# For large datasets or long-running assemblies
python cgas_module1.py -i large_dataset/ -o results/ \
    -t 32 \
    --use-nohup

# Check progress
tail -f results/03_assemblies/sample_getorganelle.log
```

### Example 5: Custom Tool Paths
```bash
# If tools are installed in non-standard locations
python cgas_module1.py -i reads/ -o output/ \
    --fastp /custom/bin/fastp \
    --getorganelle /custom/bin/get_organelle_from_reads.py \
    --bwa /custom/bin/bwa \
    --samtools /custom/bin/samtools
```

---

## FAQ

### Q1: What's the most important factor for successful assembly?
**A:** Read quality and chloroplast coverage. Clean reads with >100x chloroplast coverage give the best results. Please check getorganelle manual for detail. 

### Q2: How many GetOrganelle rounds should I use?
**A:** 15 rounds is default. Use 20-25 for difficult samples, 30+ for very challenging data. Please check getorganelle manual for detail. 

### Q3: What k-mer values work best?
**A:** Default "21,45,65,85,105" works for most samples. For difficult samples, try "21,33,55,77,99,121". Please check getorganelle manual for detail. 

### Q4: My assembly has multiple contigs. Is this bad?
**A:** 1 or 2 contig is ideal. 2 contigs can indicate SSC flip-flop (normal). >2 contigs suggests incomplete assembly or may be some indels or SNPs difference among sample. But If circular, you can resolve it. Please check getorganelle manual for detail. 

### Q5: What coverage depth is recommended?
**A:** Minimum 100x chloroplast coverage. High coverage is ideal for difficult samples. Sometime species biology also effect assembly. 

### Q6: How do I handle NovaSeq data?
**A:** Use `--trim-poly-g` to remove poly-G tails from NovaSeq two-color chemistry.

### Q7: Can I use this for mitochondrial genomes?
**A:** Yes, use `--genome-type embplant_mt` for plant mitochondrial genomes. Please check getorganelle manual for detail. However, module 2 is completley related to chloroplast. 

### Q8: What if my assembly is incomplete?
**A:** Try increasing GetOrganelle rounds, adjusting k-mer values, or checking read quality and depth. Please check getorganelle manual for detail. 

### Q9: How do I process single-end reads?
**A:** The pipeline automatically detects single-end reads. Increase rounds and k-mer values for best results. You just keep them in working directory along pairend reads. 

### Q10: What mapping rate is considered good?
**A:** The acceptable mapping rate depends on the biology of the plant under analysis. There is no fixed coverage depth that guarantees success. Based on our experience, a coverage depth of up to ~100× generally allows high-quality assembly. However, even with very high coverage (e.g., 500×), assembly issues can still occur when using our Python script or GetOrganelle, depending on genome complexity, repeats, or sequence quality. 

---

## Technical Specifications

### Performance
- **Processing speed**: ~30-60 minutes per sample (depends on data size and parameters). Here, the script just combine all analysis. Please see details on specific website of implemeted tools. A specific claim for us (CGAS developer), may not be possible. 
- **Memory usage**: 4-16 GB RAM (depends on data size and threads)
- **Disk space**: ~2-5x input data size (including intermediate files)

### Compatibility
- **Python**: 3.9+ (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating Systems**: Windows, macOS, Linux
- **Fastp**: 0.20.0 or higher
- **GetOrganelle**: 1.7.0 or higher
- **BWA**: 0.7.17 or higher
- **Samtools**: 1.15 or higher

### Input Limits
- **Max read length**: No practical limit
- **Max file size**: No practical limit
- **Max samples**: No practical limit (tested with 100+ samples)

### Quality Features
- ✅ Fastp integration for advanced QC
- ✅ GetOrganelle for specialized organellar assembly
- ✅ Automatic assembly validation
- ✅ BWA-based read mapping
- ✅ Comprehensive coverage analysis
- ✅ Batch processing
- ✅ Detailed reporting
- ✅ Error recovery

---

## References

### Fastp
- **Chen, S., et al. (2018).** Fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884-i890.
- **Fastp GitHub**: [https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp) - Official Fastp repository

### GetOrganelle
- **Jin, J., et al. (2020).** GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organellar genomes. *Genome Biology*, 21: 241.
- **GetOrganelle GitHub**: [https://github.com/Kinggerm/GetOrganelle](https://github.com/Kinggerm/GetOrganelle) - Official GetOrganelle repository

### BWA
- **Li, H., & Durbin, R. (2009).** Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics*, 25(14), 1754-1760.
- **BWA GitHub**: [https://github.com/lh3/bwa](https://github.com/lh3/bwa) - Official BWA repository

### Samtools
- **Li, H., et al. (2009).** The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078-2079.
- **Samtools GitHub**: [https://github.com/samtools/samtools](https://github.com/samtools/samtools) - Official Samtools repository

### Python Libraries
- **Biopython**: [Official Documentation](https://biopython.org/) - FASTA file parsing
- **pandas**: [Documentation](https://pandas.pydata.org/) - Data manipulation
- **openpyxl**: [Documentation](https://openpyxl.readthedocs.io/) - Excel file creation

### CGAS Suite
- **CGAS**: Chloroplast Genome Analysis Suite - integrated tools for comprehensive chloroplast genomics

### Citation
If using CGAS Module 1 in publications, please cite:
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
python cgas_module1.py --help

# 2. Verify dependencies
fastp --version
get_organelle_from_reads.py --help
bwa
samtools --version

# 3. Check files
ls raw_reads/*.fastq.gz

# 4. Test with small dataset
python cgas_module1.py -i test_reads/ -o test_output/
```

### Common Issues Solved Here
- ✅ Missing fastp? Install fastp
- ✅ GetOrganelle not found? Install GetOrganelle
- ✅ No read pairs? Check file naming
- ✅ Assembly fails? Increase rounds, adjust k-mers
- ✅ Incomplete assembly? Increase rounds, check quality
- ✅ Low mapping rate? Check assembly quality
- ✅ BWA missing? Install BWA
- ✅ Samtools missing? Install Samtools

### Quick Reference Card
```bash
# ⚡ QUICK START ⚡
cgas --module 1 -i raw_reads/ -o output/                    # cgas command
cgas-assembly -i raw_reads/ -o output/                      # shortcut command
python cgas_module1.py -i raw_reads/ -o output/             # python command
python cgas_module1.py -i raw_reads/ -o output/ -t 16       # More threads
python cgas_module1.py -i raw_reads/ -o output/ -R 20       # More rounds
python cgas_module1.py -i raw_reads/ -o output/ --force-mapping  # Force mapping

# 🔬 JUPYTER NOTEBOOK 🔬
# %run does NOT work for Module 1 — use ! operator instead
!cgas --module 1 -i raw_reads/ -o output/
!python cgas_module1.py -i raw_reads/ -o output/ -t 8

# 📊 OUTPUT 📊
# output/
# ├── 01_QC/                     # Fastp reports
# ├── 02_clean_reads/            # Cleaned reads
# ├── 03_assemblies/             # GetOrganelle outputs
# ├── 04_mapping/                # BWA mapping results
# ├── 05_cp_reads/               # Chloroplast reads
# ├── 06_reports/                # Summary reports
# └── 07_assembled_genomes/      # Complete assemblies only

# 🎯 TIPS 🎯
# - Use high-quality reads (>Q30)
# - Minimum 30x chloroplast coverage
# - 1 contig = ideal assembly
# - 2 contigs = SSC flip-flop (normal)
# - >2 contigs = incomplete assembly
# - Good mapping rate: >80%
```

---

**Last Updated**: January 2026  
**Version**: 1.0.1  
**Author**: Abdullah  
**Part of**: CGAS (Chloroplast Genome Analysis Suite)

---

*Happy Chloroplast Assembly! 🧬✨*