# CGAS Module 2: Simplified Plastome Annotation using PGA
## Complete Documentation and User Guide

---

## 📋 Table of Contents
1. [Introduction](#introduction)
2. [Quick Start - Just Run It!](#quick-start---just-run-it)
3. [Installation Guide](#installation-guide)
4. [Command-Line Usage with Examples](#command-line-usage-with-examples)
5. [Input Requirements](#input-requirements)
6. [Output Structure](#output-structure)
7. [Structural Feature Addition](#structural-feature-addition)
8. [Organism Name Handling](#organism-name-handling)
9. [Troubleshooting](#troubleshooting)
10. [Examples](#examples)
11. [FAQ](#faq)
12. [Technical Specifications](#technical-specifications)
13. [References](#references)

---

## Introduction

**CGAS Module 2** is a streamlined tool for annotating assembled chloroplast genomes using PGA (**Plastid Genome Annotator**). This module takes clean, assembled genome sequences from Module 1 and performs comprehensive gene annotation, identifies structural features (including large single copy [LSC] and small single copy [SSC] regions), and GenBank files with standardized organism names. Incorporation of organism names and accurate annotation of LSC and SSC regions are handled automatically as part of the **CGAS** workflow.

This module performs automated annotation with:
- **PGA integration**: Uses PGA for annotation
- **Batch processing**: Annotates multiple genomes simultaneously
- **Structural annotation**: Adds LSC, SSC, and IR region features
- **Organism standardization**: Handles species name mapping
- **Quality reporting**: Generates comprehensive annotation statistics
- **Publication-ready output**: GenBank files with complete annotations. However, require some manual curations using module 3 and module 4 of CGAS. 

### Key Features:
- **Simplified workflow**: Direct annotation from assembled genomes
- **PGA-powered**: Uses proven annotation tool
- **Automatic structural features**: Adds LSC/SSC/IR annotations
- **Organism name mapping**: Supports custom species names
- **Comprehensive reporting**: Statistics and warnings for each genome
- **Excel summaries**: Detailed annotation reports in Excel format

### Scientific Applications:
- **Genome annotation**: Annotate newly assembled chloroplast genomes
- **Comparative genomics**: Prepare standardized annotations for comparison
- **Database submission**: Generate GenBank files for deposition
- **Publication preparation**: Create publication-ready annotations
- **Quality control**: Assess annotation completeness and quality
- **Large-scale projects**: Batch annotation of multiple genomes

### Part of CGAS (Chloroplast Genome Analysis Suite)
CGAS is a comprehensive suite of tools designed for chloroplast genome analysis, providing seamless integration across different analysis modules.

---

## Quick Start - Just Run It!

### ⚡ One Command Does It All!

```bash
# Option 1: Using cgas command
cgas --module 2 -i assembled_genomes/ -r reference_genomes/ --pga /path/to/PGA.pl # PGA need perl therefore not install in conda environment
cgas --module 2 -i assembled_genomes/ -r reference_genomes/ --organism-file organisms.txt --pga /path/to/PGA.pl

# Option 2: Using cgas-annotate shortcut
cgas-annotate -i assembled_genomes/ -r reference_genomes/ --pga /path/to/PGA.pl

#Please look at README_PGA.md in Documents about installation of pga and setting up path for the CGAS. After path setting just run:
cgas --module 2 -i assembled_genomes/ -r reference_genomes/ --organism-file organisms.txt
cgas --module 2 -i assembled_genomes/ -r reference_genomes/
# Option 3: Using python directly
python cgas_module2.py -i assembled_genomes/ -r reference_genomes/ --pga /path/to/PGA.pl
python cgas_module2.py -i assembled_genomes/ -r reference_genomes/ --organism-file organisms.txt --pga /path/to/PGA.pl

## Create the Alias
```bash
echo 'alias pga="perl /home/abdullah/PGA/PGA.pl"' >> ~/.bashrc
source ~/.bashrc
```

> Replace `abdullah` with your username if different. Find the actual path but did not install in any other directory; so your path should follow this pattern: /home/abdullah/PGA/PGA.pl"

---
```

**What happens when you run this:**
1. ✅ Finds all FASTA files in input directory
2. ✅ Runs PGA annotation for each genome
3. ✅ Updates organism names from mapping file
4. ✅ Adds structural features (LSC, SSC, IR)
5. ✅ Creates individual GenBank files
6. ✅ Generates comprehensive Excel reports
7. ✅ Produces annotation statistics

### 📁 Real-World Folder Examples

#### Example 1: Basic Usage (Easiest!)
```bash
# Your folder structure:
/home/abdullah/annotation_project/
├── assembled_genomes/
│   ├── SRR8666784_1.fasta
│   ├── SRR8666785_1.fasta
│   └── SRR8666786_1.fasta
├── reference_genomes/
│   └── Nicotiana_tabacum.gb
└── organisms.txt

# Navigate to project folder
cd /home/abdullah/annotation_project/

# Run annotation
python cgas_module2.py -i assembled_genomes/ -r reference_genomes/ --organism-file organisms.txt #Any command can be used based on your choice

# Output created automatically:
# module2_annotations/
# ├── Annotated_GenBank/
# │   ├── Hibiscus_coatesii_SRR8666784_1.gb
# │   ├── Hibiscus_goldsworthii_SRR8666785_1.gb
# │   └── Hibiscus_schizopetalus_SRR8666786_1.gb
# ├── Annotation_Logs/
# ├── Reports/
# │   ├── 00_ANNOTATION_SUMMARY.tsv
# │   └── 00_ANNOTATION_SUMMARY.xlsx
```

#### Example 2: Quick Annotation (No organism mapping)
```bash
# Simple annotation without species names
python cgas_module2.py -i assembled_genomes/ -r reference_genomes/

# The fasta header will place instead of  Organism If organims file is not provided. 
```

#### Example 3: Custom Output Directory
```bash
# Full control over directories
python cgas_module2.py -i data/ -r refs/ -o my_annotations/ --organism-file species.txt

# Input from: data/ and refs/
# Output to: my_annotations/
```

### 🎯 Quick Start Cheat Sheet

```bash
# ====== cgas command ======
cgas --module 2 -i assembled_genomes/ -r refs/ --pga /path/to/PGA.pl                                       # Basic
cgas --module 2 -i genomes/ -r refs/ --organism-file species.txt --pga /path/to/PGA.pl                     # With organism names
cgas --module 2 -i data/ -r refs/ -o results/ --pga /path/to/PGA.pl                                        # Custom output
cgas --module 2 -i genomes/ -r refs/ --force-rerun --pga /path/to/PGA.pl                                   # Force re-annotate

# ====== cgas-annotate shortcut ======
cgas-annotate -i assembled_genomes/ -r refs/ --pga /path/to/PGA.pl                                         # Basic
cgas-annotate -i genomes/ -r refs/ --organism-file species.txt --pga /path/to/PGA.pl                       # With organism names

# ====== python command ======
python cgas_module2.py -i assembled_genomes/ -r reference_genomes/                                         # Basic
python cgas_module2.py -i genomes/ -r refs/ --organism-file species.txt                                    # With organism names
python cgas_module2.py -i data/ -r refs/ -o results/                                                       # Custom output
python cgas_module2.py -i genomes/ -r refs/ --force-rerun                                                  # Force re-annotate
python cgas_module2.py --help                                                                              # Get help
```

### 📊 What You Get (Output Files)

```
module2_annotations/                           # Created automatically
├── 📁 Annotated GenBank Files
│   ├── Annotated_GenBank/                     # Individual GenBank files
│   │   ├── Species1_Sample1.gb              # Annotated genome 1
│   │   ├── Species2_Sample2.gb              # Annotated genome 2
│   │   └── Species3_Sample3.gb              # Annotated genome 3
│
├── 📁 Annotation Logs
│   ├── Annotation_Logs/                      # PGA warning logs
│   │   ├── Species1_Sample1_warning.log     # PGA warnings
│   │   ├── Species2_Sample2_warning.log     # PGA warnings
│   │   └── Species3_Sample3_warning.log     # PGA warnings
│
└── 📊 Reports
    ├── Reports/                              # Summary reports
    │   ├── 00_ANNOTATION_SUMMARY.tsv         # Tab-delimited summary
    │   └── 00_ANNOTATION_SUMMARY.xlsx         # Excel summary
```

---

## Installation Guide

> **Note:** If you have already installed any version of CGAS tool (using `environment.yml`or `environment-minimal.yml`). Then just install pearl (If not installed as wsl include by default in most cases) and PGA.

### Prerequisites
- **Python 3.9 or higher** (tested on Python 3.9–3.12)
- **BLAST+ 2.8.1 or higher** (required by PGA)
- **Perl 5 or higher** (required by PGA)
- **PGA (Plastid Genome Annotator)**

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

#### 2. Install BLAST+ (Required by PGA)

**Ubuntu/Debian:**
```bash
#All methods are provided to the best of our knowledge and rely on external tools. Users are encouraged to consult the original tool providers for detailed configuration options and updates.
sudo apt-get update 
sudo apt-get install ncbi-blast+
```

**macOS:**
```bash
brew install blast
```

**CentOS/RHEL:**
```bash
sudo yum install ncbi-blast+
```

**From source (all platforms):**
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
tar -zxvf ncbi-blast-*.tar.gz
cd ncbi-blast-*/
./configure
make
sudo make install
```

#### 3. Install Perl (Required by PGA)

**Ubuntu/Debian:**
```bash
sudo apt-get install perl
```

**macOS:**
```bash
# Perl usually comes pre-installed
perl --version
```

**Windows:**
- Download from https://www.perl.org/get.html
- Or use WSL (Windows Subsystem for Linux)

#### 4. Install PGA (Plastid Genome Annotator)

```bash
# Download PGA
wget https://github.com/quayc/pga/raw/master/PGA.pl

##Create the Alias
```bash
echo 'alias pga="perl /home/abdullah/PGA/PGA.pl"' >> ~/.bashrc
source ~/.bashrc
```
# Replace `abdullah` with your username and install in home directory instead of making any specific folder. 
# Test PGA
perl PGA.pl -h 

#Make executable
chmod +x PGA.pl
#If not work then 
# Test PGA
perl PGA.pl -h
```

**Alternative installation:**
```bash
# Clone from GitHub
git clone https://github.com/quayc/pga.git
cd pga
chmod +x PGA.pl
```

#### 5. Usage of CGAS Module 2

# Make it executable (optional)
chmod +x cgas_module2.py

# Verify
python cgas_module2.py --help
```

#### 6. Verify Complete Installation
```bash
# Test Python
python --version                    # Should be 3.9+

# Test packages
python -c "import Bio, pandas, openpyxl"  # Should not give errors

# Test BLAST
blastn -version                       # Should show BLAST version

# Test Perl
perl --version                        # Should show Perl version

# Test PGA
perl PGA.pl -h                        # Should show PGA help
```

### Dependency Details

| Package | Version | Purpose | Required |
|---------|---------|---------|----------|
| **Python** | ≥3.9 | Script execution | Yes |
| **Biopython** | ≥1.79 | GenBank file parsing | Yes |
| **pandas** | ≥1.3.0 | Data manipulation | Yes |
| **openpyxl** | ≥3.0.9 | Excel file generation | Yes |
| **BLAST+** | ≥2.8.1 | Sequence similarity search | Yes (PGA) |
| **Perl** | ≥5 | PGA execution | Yes (PGA) |
| **PGA** | Latest | Plastid genome annotation | Yes |

---

## Command-Line Usage with Examples

### Complete Command Reference

```bash
# ====================================================================
# USING cgas COMMAND
# ====================================================================

# Basic annotation
cgas --module 2 -i /home/abdullah/7_assembled_genomes/ -r /home/abdullah/reference_genomes/ --pga /path/to/PGA.pl

# With organism name mapping
cgas --module 2 -i 7_assembled_genomes/ -r reference_genomes/ --organism-file organisms.txt --pga /path/to/PGA.pl

# Custom output directory
cgas --module 2 -i genomes/ -r refs/ -o annotations/ --pga /path/to/PGA.pl

# Force re-annotation of all files
cgas --module 2 -i genomes/ -r refs/ --force-rerun --pga /path/to/PGA.pl
```

```bash
# ====================================================================
# USING cgas-annotate SHORTCUT
# ====================================================================

# Basic annotation
cgas-annotate -i 7_assembled_genomes/ -r reference_genomes/ --pga /path/to/PGA.pl

# With all options
cgas-annotate \
  -i assembled_genomes/ \
  -o annotations/ \
  -r reference_genomes/ \
  --organism-file organisms.txt \
  --pga /path/to/PGA.pl
```

```bash
# ====================================================================
# USING python COMMAND
# ====================================================================

# ====================================================================
# BASIC USAGE
# ====================================================================

# 1. Basic annotation (required parameters)
python cgas_module2.py -i assembled_genomes/ -r reference_genomes/

# 2. With organism name mapping
python cgas_module2.py -i genomes/ -r refs/ --organism-file species.txt

# 3. Custom output directory
python cgas_module2.py -i data/ -r refs/ -o annotations/

# 4. Force re-annotation of all files
python cgas_module2.py -i genomes/ -r refs/ --force-rerun

# 5. Get help
python cgas_module2.py --help


# ====================================================================
# REAL-WORLD EXAMPLES
# ====================================================================

# Example 1: Annotate Hibiscus species
python cgas_module2.py -i hibiscus_assemblies/ -r Malva.gb --organism-file hibiscus_species.txt

# Example 2: Batch annotation of multiple genera
python cgas_module2.py -i all_assemblies/ -r reference_set/ -o all_annotations/

# Example 3: Quick test with single genome
python cgas_module2.py -i test_genome.fasta -r reference.gb 

# Example 4: Custom PGA path
python cgas_module2.py -i genomes/ -r refs/ --pga /custom/path/PGA.pl
```

### Parameter Details

| Parameter | Short | Required | Default | Description |
|-----------|-------|----------|---------|-------------|
| `--input` | `-i` | Yes | - | Directory with assembled FASTA files |
| `--reference` | `-r` | Yes | - | Directory with reference GenBank files |
| `--output` | `-o` | No | `module2_annotations` | Output directory |
| `--organism-file` | - | Yes | - | specify as mentioned in command line |
| `--pga` | - | No | `PGA.pl` | Path to PGA.pl script |
| `--min-ir` | - | No | 1000 | Minimum IR length (bp) |
| `-p` | - | No | 40 | Minimum percent identity |
| `-q` | - | No | 0.5,2 | Query coverage range |
| `-f` | - | No | circular | Genome form (circular/linear) |
| `--force-rerun` | - | No | False | Force re-annotation of all samples |
| `--log-level` | - | No | INFO | Logging level |

---

---

## Input Requirements

### Supported File Formats

**Input FASTA files:**
- `.fasta`
- `.fa`
- `.fna`

**Reference GenBank files:**
- `.gb`
- `.gbk`
- `.genbank`

### Critical Requirements

**1. Assembled FASTA Files**
```fasta
>SRR8666784_1
ATGGCGACGACGTTCGTCGTCGTTTGTCGATCTCGTCTGACTTCAGCCTGATCGGTAGCA
ATGGTAAGTTGGTGGTGTGAAAGCAGCTGACGGGAGCATTCGGATGTAGATTTGGAGAAA
...
```

**Requirements:**
- Complete or near-complete chloroplast genome assembly
- Clean sequence (no N's in critical regions)
- Standard nucleotides (A, T, G, C)
- Reasonable length (120-200 kb typical for chloroplasts)

**2. Reference GenBank Files**
```gb
LOCUS       NC_007877               155939 bp DNA     circular PLN 15-JUL-2022
DEFINITION  Nicotiana tabacum chloroplast, complete genome.
ACCESSION   NC_007877
VERSION     NC_007877.1
...
FEATURES             Location/Qualifiers
     source          1..155939
                     /organism="Nicotiana tabacum"
                     /organelle="plastid"
     gene            864..2690
                     /gene="rbcL"
     CDS             864..2690
                     /gene="rbcL"
                     /product="ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"
ORIGIN
        1 atggcgacga cgttcgtcgt cgtttgtcga tctcgtctga cttcagcctg atcggtagca
...
//
```

**Requirements:**
- Well-annotated chloroplast genome
- Complete gene set
- Standard gene nomenclature
- Preferably from closely related species

**3. Organism Mapping File (Optional)**

**Format (txt):**
```txt
accession	organism
SRR8666784	Hibiscus coatesii
SRR8666785	Hibiscus goldsworthii
SRR8666786	Hibiscus schizopetalus
```

**Requirements:**
- Tab-separated values
- Two columns: accession and organism
- No header line (or header will be skipped)
- Accession should match FASTA header prefix

### File Organization

```bash
project_folder/
├── assembled_genomes/          # Input FASTA files
│   ├── SRR8666784_1.fasta
│   ├── SRR8666785_1.fasta
│   └── SRR8666786_1.fasta
├── reference_genomes/          # Reference GenBank files
│   ├── Nicotiana_tabacum.gb
│   └── Arabidopsis_thaliana.gb
└── organisms.txt               # Optional organism mapping
```

### Reference Selection Guidelines

**Best practices:**
- Choose closely related species (same family/genus)
- Prefer well-annotated reference genomes
- Consider genome size similarity
- Use reference with complete gene set
- Avoid references with annotation issues

**Reference quality indicators:**
- Complete gene annotation
- Standard gene names
- No pseudogene issues
- Recent annotation version
- Published in peer-reviewed journal

---

## Output Structure

### Directory Organization

```
module2_annotations/                           # Created automatically
├── 📁 Annotated GenBank Files
│   ├── Annotated_GenBank/                     # Individual GenBank files
│   │   ├── Hibiscus_coatesii_SRR8666784_1.gb   # Annotated genome 1
│   │   │   └── FEATURES:
│   │   │       source          1..160000
│   │   │                       /organism="Hibiscus coatesii"
│   │   │                       /organelle="plastid"
│   │   │       gene            864..2690
│   │   │                       /gene="rbcL"
│   │   │       CDS             864..2690
│   │   │                       /gene="rbcL"
│   │   │                       /product="ribulose-1,5-bisphosphate..."
│   │   │       misc_feature    1..88000
│   │   │                       /note="large single copy (LSC)"
│   │   │       misc_feature    88001..95000
│   │   │                       /note="inverted repeat B (IRB)"
│   │   │       misc_feature    95001..155000
│   │   │                       /note="small single copy (SSC)"
│   │   │       misc_feature    155001..160000
│   │   │                       /note="inverted repeat A (IRA)"
│   │   ├── Hibiscus_goldsworthii_SRR8666785_1.gb
│   │   └── Hibiscus_schizopetalus_SRR8666786_1.gb
│
├── 📁 Annotation Logs
│   ├── Annotation_Logs/                      # PGA warning logs
│   │   ├── Hibiscus_coatesii_SRR8666784_1_warning.log
│   │   │   └── Content:
│   │   │       Total number of genes in the reference plastome(s): 85
│   │   │       Total number of genes annotated in the target plastome: 83
│   │   │       Warning: trnT-GGT not annotated
│   │   │       Warning: ycf15 not annotated
│   │   ├── Hibiscus_goldsworthii_SRR8666785_1_warning.log
│   │   └── Hibiscus_schizopetalus_SRR8666786_1_warning.log
│
└── 📊 Reports
    ├── Reports/                              # Summary reports
    │   ├── 00_ANNOTATION_SUMMARY.tsv         # Tab-delimited summary
    │   │   └── Content:
    │   │       Sample	Status	Genes_Reference	Genes_Annotated	Annotation_Rate_%	Warnings	Unannotated_Genes	GenBank_File
    │   │       SRR8666784_1	SUCCESS	85	83	97.65	2	trnT-GGT, ycf15	Hibiscus_coatesii_SRR8666784_1.gb
    │   │       SRR8666785_1	SUCCESS	85	84	98.82	1	trnT-UGU	Hibiscus_goldsworthii_SRR8666785_1.gb
    │   │       SRR8666786_1	SUCCESS	85	85	100.00	0		Hibiscus_schizopetalus_SRR8666786_1.gb
    │   └── 00_ANNOTATION_SUMMARY.xlsx         # Excel summary with formatting
```

### Key Output Files Explained

#### 1. Annotated GenBank Files (.gb)

**Complete annotation with structural features:**
```gb
LOCUS       Sample1_1               160000 bp DNA     circular PLN 15-JUL-2022
DEFINITION  Sample1_1 chloroplast, complete genome.
ACCESSION   Sample1_1
VERSION     Sample1_1.1
KEYWORDS    RefSeq.
SOURCE      Sample1_1
  ORGANISM  Hibiscus coatesii
            Eukaryota; Viridiplantae; Streptophyta; Embryophyta;
            Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons;
            core eudicots; rosids; Malvales; Malvaceae; Hibiscus.
FEATURES             Location/Qualifiers
     source          1..160000
                     /organism="Hibiscus coatesii"
                     /organelle="plastid"
                     /mol_type="genomic DNA"
     misc_feature    1..88000
                     /note="large single copy (LSC)"
     misc_feature    88001..95000
                     /note="inverted repeat B (IRB)"
     misc_feature    95001..155000
                     /note="small single copy (SSC)"
     misc_feature    155001..160000
                     /note="inverted repeat A (IRA)"
     gene            864..2690
                     /gene="rbcL"
     CDS             864..2690
                     /gene="rbcL"
                     /product="ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"
                     /protein_id="Sample1_1_001"
                     /translation="MASSSNSS..."
ORIGIN
        1 atggcgacga cgttcgtcgt cgtttgtcga tctcgtctga cttcagcctg atcggtagca
...
//
```

**Key features:**
- Updated organism name
- Added structural region annotations (LSC, SSC, IRa, IRb)
- Complete gene annotation from PGA
- Standard GenBank format

#### 2. Annotation Logs (.log)

**PGA execution logs:**
```
Total number of genes in the reference plastome(s): 85
Total number of genes annotated in the target plastome: 83
Warning: trnT-GGT not annotated
Warning: ycf15 not annotated
```

**Information provided:**
- Reference gene count
- Annotated gene count
- Annotation success rate
- List of unannotated genes
- PGA warnings and errors

---

## Structural Feature Addition

### IR Region Detection

**PGA IR detection:**
- PGA automatically detects IR regions
- Outputs as repeat_region features
- May use various naming conventions

**Standardization process:**
```python
# Update IR annotations to standard wording
ir_features_sorted[0].qualifiers["note"] = ["inverted repeat B (IRB)"]
ir_features_sorted[1].qualifiers["note"] = ["inverted repeat A (IRA)"]
```

**IR naming convention:**
- **IRB**: First IR (typically 3' end)
- **IRA**: Second IR (typically 5' end)
- Based on genomic position, not orientation

### LSC/SSC Inference

**Gap analysis:**
```python
# Calculate regions between IRs
region1_start = 1
region1_end = irb_start - 1
region2_start = irb_end + 1
region2_end = ira_start - 1
region3_start = ira_end + 1
region3_end = genome_length

# Determine which is LSC (larger) and SSC (smaller)
if region1_length > region2_length:
    lsc_start, lsc_end = region1_start, region1_end
    ssc_start, ssc_end = region2_start, region2_end
else:
    lsc_start, lsc_end = region2_start, region2_end
    ssc_start, ssc_end = region1_start, region1_end
```

**LSC/SSC characteristics:**
- **LSC**: Large Single Copy (80-90 kb typical)
- **SSC**: Small Single Copy (15-30 kb typical)
- **IR**: Inverted Repeats (20-30 kb each)

### Feature Integration

**Insertion order:**
```python
# Insert LSC and SSC features after source feature
source_idx = 0
for i, feat in enumerate(record.features):
    if feat.type == "source":
        source_idx = i + 1
        break

# Add features in order
record.features.insert(source_idx, lsc_feature)
record.features.insert(source_idx + 1, ssc_feature)
```

**Feature format:**
```gb
misc_feature    1..88000
                /note="large single copy (LSC)"
misc_feature    88001..95000
                /note="inverted repeat B (IRB)"
misc_feature    95001..155000
                /note="small single copy (SSC)"
misc_feature    155001..160000
                /note="inverted repeat A (IRA)"
```

---

## Organism Name Handling

### Organism Mapping File

**Format specification:**
```tsv
# Tab-separated values
# No header line (automatically skipped if present)
accession	organism
SRR8666784	Hibiscus coatesii
SRR8666785	Hibiscus goldsworthii
SRR8666786	Hibiscus schizopetalus
```


#### 2. PGA Not Found
```bash
❌ ERROR: PGA not found or not working at: PGA.pl
```
**Solution:**
```bash
# Download PGA
wget https://github.com/quayc/pga/raw/master/PGA.pl
chmod +x PGA.pl

# Test PGA
perl PGA.pl -h

# Or specify custom path
python cgas_module2.py --pga /path/to/PGA.pl -i genomes/ -r refs/
```

#### 3. No FASTA Files Found
```bash
❌ No FASTA files found in /path/to/directory
```
**Solution:**
```bash
# Check file extensions
ls *.fasta *.fa *.fna

# Check directory
ls -la /path/to/assembled_genomes/

# Specify correct directory
python cgas_module2.py -i /full/path/to/genomes/ -r refs/
```

#### 4. No Reference GenBank Files
```bash
❌ No GenBank reference files found in: /path/to/reference
```
**Solution:**
```bash
# Check reference directory
ls -la /path/to/reference_genomes/

# Check file extensions
ls *.gb *.gbk *.genbank

# Add reference files
cp /path/to/reference.gb /path/to/reference_genomes/
```

#### 5. PGA Annotation Fails
```bash
⚠ WARNING: PGA failed for sample1
```
**Solution:**
```
Common PGA issues:
- Poor assembly quality
- Incomplete genome
- Distant reference species
- Unusual genome structure

Check the warning log for details and consider:
- Using closer reference species
- Improving assembly quality
- Manual annotation review
```

#### 6. Low Annotation Rate
```bash
⚠ WARNING: Low annotation rate: 65%
```
**Solution:**
```
This indicates problems:
- Check assembly completeness
- Verify reference suitability
- Review FASTA file quality
- Consider different reference species
```

### Quick Fixes Table

| Problem | Symptom | Quick Fix |
|---------|---------|-----------|
| Missing BLAST | "BLAST not found" | Install BLAST+ |
| PGA not found | "PGA not working" | Download PGA.pl |
| No FASTA files | "No FASTA files" | Check input directory |
| No reference | "No reference files" | Add reference GenBank |
| PGA fails | "PGA failed" | Check assembly quality |
| Low rate | "Low annotation rate" | Use closer reference |

---

## Examples

### Example 1: Complete Annotation Workflow
```bash
# 1. Prepare project structure
mkdir -p /home/abdullah/hibiscus_annotation/
cd /home/abdullah/hibiscus_annotation/

# 2. Organize files
mkdir -p assembled_genomes reference_genomes
# Copy FASTA files to assembled_genomes/
# Copy reference GenBank to reference_genomes/

# 3. Create organism mapping
cat > organisms.txt << EOF
SRR8666784	Hibiscus coatesii
SRR8666785	Hibiscus goldsworthii
SRR8666786	Hibiscus schizopetalus
EOF

# 4. Run annotation
python cgas_module2.py -i assembled_genomes/ -r reference_genomes/ --organism-file organisms.txt

# 5. Check results
ls module2_annotations/
open module2_annotations/Reports/00_ANNOTATION_SUMMARY.xlsx
```

### Example 2: Quick Test Annotation
```bash
# Quick test with single genome
python cgas_module2.py -i test_genome.fasta -r reference.gb

# Check if annotation succeeded
ls module2_annotations/Annotated_GenBank/
```

### Example 3: Batch Annotation of Multiple Genera
```bash
# Annotate multiple genera with different references
python cgas_module2.py -i malvales_assemblies/ -r malvales_reference.gb -o malvales_annotations/
python cgas_module2.py -i rosaceae_assemblies/ -r rosaceae_reference.gb -o rosaceae_annotations/
```

### Example 4: Custom PGA Parameters
```bash
# Use custom PGA parameters for difficult genomes
python cgas_module2.py -i genomes/ -r refs/ -p 30 -q 0.3,3.0 --min-ir 800
```

### Example 5: Re-annotation After Updates
```bash
# Force re-annotation of all files
python cgas_module2.py -i genomes/ -r refs/ --force-rerun
```

---

## FAQ

### Q1: What's the most important factor for successful annotation?
**A:** Assembly quality. Clean, complete assemblies give the best results.

### Q2: How do I choose a reference genome?
**A:** Use a closely related species with complete, well-annotated chloroplast genome.

### Q3: What annotation rate is considered good?
**A:** >90% is good, >95% is excellent. <80% may need review. Please follow module 3 and 4 for finilizing the annotations before manual visulization.

### Q4: Can I use mitochondrial genomes as reference?
**A:** PGA is designed for chloroplast genomes. Use chloroplast references.

### Q5: What if PGA misses some genes?
**A:** Check the warning log. Common issues: assembly gaps, distant reference, unusual genome structure.

### Q6: Can I customize PGA parameters?
**A:** Yes, use parameters like `-p` (identity), `-q` (coverage), `--min-ir` (IR length).

### Q7: How are structural features added?
**A:** Automatically detected from PGA IR annotations, LSC/SSC inferred from gaps.

### Q8: What if I have multiple reference genomes?
**A:** Place all in reference directory. PGA will use all available references.

### Q9: Can I annotate non-circular genomes?
**A:** Yes, use `-f linear` parameter.

### Q10: How do I handle very large genomes?
**A:** PGA can handle larger genomes, but may need adjusted parameters.

---

## Technical Specifications

### Performance
- **Processing speed**: Annotation requires approximately ~1 minute per genome, depending on genome size, reference selection, and available computational resources. Detailed performance characteristics are provided on the PGA website. The CGAS-specific integrations, including annotation of large single copy (LSC) and small single copy (SSC) regions and incorporation of species names from a text file, require only a few additional seconds per genome.
- **Memory usage**: <500 MB RAM per annotation
- **Disk space**: ~10 MB per annotated genome

### Compatibility
- **Python**: 3.9+ (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating Systems**: Windows, macOS, Linux
- **PGA**: Compatible with latest version
- **BLAST+**: 2.8.1 or higher required

### Input Limits
- **Max genome size**: No practical limit (tested with 300kb+ genomes)
- **Max genomes**: No practical limit (tested with 100+ genomes)
- **Reference genomes**: No limit (all references used)

### Quality Features
- ✅ PGA integration for proven annotation
- ✅ Automatic structural feature addition
- ✅ Organism name standardization
- ✅ Comprehensive reporting
- ✅ Batch processing
- ✅ Quality metrics
- ✅ Error recovery

---

## References

### PGA (Plastid Genome Annotator)
- **Qu, C. J., et al. (2019).** PGA: a software package for rapid and accurate annotation of plastid genomes. *Plant Methods*, 15(1), 1-12.
- **PGA GitHub**: [https://github.com/quayc/pga](https://github.com/quayc/pga) - Official PGA repository

### BLAST+
- **NCBI**: [BLAST+ User Manual](https://www.ncbi.nlm.nih.gov/books/NBK279684/) - BLAST+ documentation

### Python Libraries
- **Biopython**: [Official Documentation](https://biopython.org/) - GenBank file parsing
- **pandas**: [Documentation](https://pandas.pydata.org/) - Data manipulation
- **openpyxl**: [Documentation](https://openpyxl.readthedocs.io/) - Excel file creation

### CGAS Suite
- **CGAS**: Chloroplast Genome Analysis Suite - integrated tools for comprehensive chloroplast genomics

### Citation
If using CGAS Module 2 in publications, please cite:
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
python cgas_module2.py --help

# 2. Verify dependencies
blastn -version
perl --version
python -c "import Bio, pandas, openpyxl"

# 3. Check files
ls assembled_genomes/*.fasta
ls reference_genomes/*.gb

# 4. Test PGA
perl PGA.pl -h
```

### Common Issues Solved Here
- ✅ Missing BLAST? Install BLAST+
- ✅ PGA not found? Download PGA.pl
- ✅ No FASTA files? Check input directory
- ✅ No reference? Add reference GenBank files
- ✅ PGA fails? Check assembly quality
- ✅ Low rate? Use closer reference species

### Quick Reference Card
```bash
# ⚡ QUICK START ⚡
cgas --module 2 -i genomes/ -r refs/ --pga /path/to/PGA.pl                                  # cgas command
cgas-annotate -i genomes/ -r refs/ --pga /path/to/PGA.pl                                    # shortcut command
python cgas_module2.py -i genomes/ -r refs/                                                  # python command
python cgas_module2.py -i genomes/ -r refs/ --organism-file species.txt                      # With names
python cgas_module2.py -i genomes/ -r refs/ --force-rerun                                    # Re-annotate

# 🔬 JUPYTER NOTEBOOK 🔬
%run cgas_module2.py -i genomes/ -r refs/                                                    # %run works for Module 2
!cgas --module 2 -i genomes/ -r refs/ --pga /path/to/PGA.pl                                 # ! also works

# 📊 OUTPUT 📊
# module2_annotations/
# ├── Annotated_GenBank/          # Annotated genomes
# ├── Annotation_Logs/           # PGA logs
# └── Reports/                   # Summary reports
#     ├── 00_ANNOTATION_SUMMARY.tsv
#     └── 00_ANNOTATION_SUMMARY.xlsx

# 🎯 TIPS 🎯
# - Use clean, complete assemblies
# - Choose closely related reference species
# - Good annotation rate: >90%
# - Check warning logs for issues
# - Structural features added automatically
# - Organism names updated from mapping file
```

---

**Last Updated**: January 2026  
**Version**: 1.0.1  
**Author**: Abdullah  
**Part of**: CGAS (Chloroplast Genome Analysis Suite)

---

*Happy Plastome Annotation! 🧬✨*