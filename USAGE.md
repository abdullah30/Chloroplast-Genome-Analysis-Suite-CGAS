# CGAS Command-Line Usage Guide

Comprehensive guide for using the Chloroplast Genome Analysis Suite (CGAS).

## Table of Contents
- [Installation](#installation)
- [Basic Usage](#basic-usage)
- [Output Structure](#output-structure)
- [Individual Module Commands](#individual-module-commands)
- [Module-Specific Options](#module-specific-options)
- [Advanced Usage](#advanced-usage)
- [Practical Examples](#practical-examples)
- [Troubleshooting](#troubleshooting)

---

## Installation

### Method 1: Install from GitHub

```bash
pip install git+https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git
```

### Method 2: Clone and Install

```bash
cd ~
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS
cd Chloroplast-Genome-Analysis-Suite-CGAS
pip install -r requirements.txt
pip install -e .
```

---

## Basic Usage

### Run All Modules

```bash
# Navigate to directory with your GenBank/FASTA files
cd /path/to/your/data

# Run all available modules
cgas
```

### Run Specific Module

```bash
# Run only Module 1 (Gene Count)
cgas --module 1

# Run only Module 4 (Codon Usage)
cgas --module 4

# Run Module 8 (SSR Analysis)
cgas --module 8
```

### List Available Modules

```bash
# Show all modules with descriptions
cgas --list
```

---

## Output Structure

CGAS organizes outputs into dedicated folders for each module:

```
your_working_directory/
├── [Input Files]
│   ├── species1.gb
│   ├── species2.gb
│   └── species3.gb
│
├── [Module Outputs]
│   ├── Module1_Gene_Count_Analysis/
│   │   ├── Chloroplast_Gene_Analysis_20251226_143022.xlsx
│   │   └── Gene_Normalization_Log.xlsx
│   │
│   ├── Module2_Gene_Table/
│   │   └── Gene_Table_20251226_143045.docx
│   │
│   ├── Module3_Comparative_Analysis/
│   │   └── Comparative_Genome_Analysis.xlsx
│   │
│   ├── Module4_Codon_Usage/
│   │   └── Codon_Usage_RSCU_20251226_143115.xlsx
│   │
│   ├── Module5_Amino_Acid/
│   │   └── Amino_Acid_Composition_20251226_143145.xlsx
│   │
│   ├── Module6_SNP_Analysis/
│   │   └── SNP_Analysis_20251226_143215.xlsx
│   │
│   ├── Module7_Intron_Analysis/
│   │   └── Intron_Analysis_20251226_143245.xlsx
│   │
│   ├── Module8_SSR_Analysis_20251226_143315/
│   │   ├── SSR_Summary.xlsx
│   │   ├── SSR_Detailed.xlsx
│   │   └── SSR_Statistics.xlsx
│   │
│   └── Module9_Nucleotide_Diversity/
│       ├── Nucleotide_Diversity_Summary.xlsx
│       └── alignment_files/
```

**Benefits:**
- ✅ Clean separation of inputs and outputs
- ✅ Easy to locate specific analysis results
- ✅ No clutter in working directory
- ✅ Timestamped outputs prevent overwriting

---

## Individual Module Commands

Each module can be run independently:

### Module 1: Gene Count and Summary
```bash
cgas --module 1

# Output location:
# Module1_Gene_Count_Analysis/
#   ├── Chloroplast_Gene_Analysis_TIMESTAMP.xlsx
#   └── Gene_Normalization_Log.xlsx
```

**Output includes:**
- Complete gene inventory
- Gene counts per genome
- IR-duplicated genes
- Pseudogene detection
- Gene name normalization log

### Module 2: Gene Table Generation
```bash
cgas --module 2

# Output location:
# Module2_Gene_Table/
#   └── Gene_Table_TIMESTAMP.docx
```

**Output includes:**
- Publication-ready Word document
- Formatted gene tables
- Gene locations and annotations

### Module 3: Comparative Genome Analysis
```bash
cgas --module 3

# Output location:
# Module3_Comparative_Analysis/
#   └── Comparative_Genome_Analysis.xlsx
```

**Output includes:**
- Genome length comparisons (LSC, SSC, IR)
- GC content analysis
- Region-specific statistics
- Professional formatting with headers

### Module 4: Codon Usage Analysis
```bash
cgas --module 4

# Output location:
# Module4_Codon_Usage/
#   └── Codon_Usage_RSCU_TIMESTAMP.xlsx
```

**Output includes:**
- RSCU values for all codons
- Codon frequency tables
- Statistical summaries

### Module 5: Amino Acid Composition
```bash
cgas --module 5

# Output location:
# Module5_Amino_Acid/
#   └── Amino_Acid_Composition_TIMESTAMP.xlsx
```

**Output includes:**
- Amino acid frequencies
- Compositional analysis
- Cross-genome comparisons

### Module 6: SNP/Substitution Analysis
```bash
cgas --module 6

# Requires FASTA alignment files

# Output location:
# Module6_SNP_Analysis/
#   └── SNP_Analysis_TIMESTAMP.xlsx
```

**Output includes:**
- SNP/substitution matrices
- Variant positions
- Frequency analysis

### Module 7: Intron Extraction
```bash
cgas --module 7

# Output location:
# Module7_Intron_Analysis/
#   └── Intron_Analysis_TIMESTAMP.xlsx
```

**Output includes:**
- Intron/exon structures
- Intron sequences
- Gene architecture details

### Module 8: SSR Analysis
```bash
# Default thresholds (10,5,4,3,3,3)
cgas --module 8

# Custom thresholds
cgas --module 8 -t 12,6,5,4,4,4

# Stricter thresholds
cgas --module 8 -t 15,7,6,5,4,4

# Output location:
# Module8_SSR_Analysis_TIMESTAMP/
#   ├── SSR_Summary.xlsx
#   ├── SSR_Detailed.xlsx
#   └── SSR_Statistics.xlsx
```

**Threshold format:** `-t mono,di,tri,tetra,penta,hexa`
- Mononucleotide: minimum repeat count (default: 10)
- Dinucleotide: minimum repeat count (default: 5)
- Trinucleotide: minimum repeat count (default: 4)
- Tetranucleotide: minimum repeat count (default: 3)
- Pentanucleotide: minimum repeat count (default: 3)
- Hexanucleotide: minimum repeat count (default: 3)

**Output includes:**
- SSR summary by type and region
- Detailed SSR list with locations
- Distribution statistics
- Region-specific analysis (LSC/SSC/IR)

### Module 9: Nucleotide Diversity
```bash
cgas --module 9

# Requires:
# - MAFFT installed
# - 2+ GenBank files

# Output location:
# Module9_Nucleotide_Diversity/
#   ├── Nucleotide_Diversity_Summary.xlsx
#   └── alignment_files/
```

**Output includes:**
- Pi (nucleotide diversity)
- Theta (Watterson's estimator)
- Gene-by-gene diversity
- Alignment files for verification

---

## Module-Specific Options

### SSR Analysis (Module 8)

```bash
# Standard usage
cgas --module 8

# Custom thresholds
cgas --module 8 -t 10,5,4,3,3,3

# Very strict (fewer, longer SSRs)
cgas --module 8 -t 20,10,8,6,5,5

# More permissive (more, shorter SSRs)
cgas --module 8 -t 8,4,3,3,3,3
```

**Choosing thresholds:**
- **Standard (10,5,4,3,3,3):** Balanced for most studies
- **Strict (15,7,6,5,4,4):** Focus on longer, more stable SSRs
- **Permissive (8,4,3,3,3,3):** Capture more potential SSRs

---

## Advanced Usage

### Working with Jupyter Notebook

```python
# In Jupyter Notebook - use %run command

# Run complete analysis
%run "unified_analyzer.py"

# Run individual modules
%run "module1_gene_count.py"
%run "module3_comparative_analysis.py"
%run "module8_ssr_analysis.py" -t 10,5,4,3,3,3
```

⚠️ **Important:** Always use `%run`, not `import`

### Batch Processing Multiple Datasets

```bash
#!/bin/bash
# Process multiple genome sets

for dataset in dataset1 dataset2 dataset3; do
    echo "Processing $dataset..."
    cd $dataset
    cgas
    cd ..
done
```

### Selective Module Execution

```bash
# Quick gene count only
cd my_genomes/
cgas --module 1

# Codon usage for comparative study
cd comparative_study/
cgas --module 4

# SSR mining with custom settings
cd ssr_project/
cgas --module 8 -t 12,6,5,4,4,4
```

---

## Practical Examples

### Example 1: Complete Genome Analysis

```bash
# Setup
mkdir ~/chloroplast_project
cd ~/chloroplast_project

# Add your GenBank files
cp /path/to/genomes/*.gb .

# Verify files
ls *.gb

# Run complete analysis
cgas

# Check results
ls -d Module*_*/
```

**Output folders created:**
```
Module1_Gene_Count_Analysis/
Module2_Gene_Table/
Module3_Comparative_Analysis/
Module4_Codon_Usage/
Module5_Amino_Acid/
Module6_SNP_Analysis/        (if FASTA files present)
Module7_Intron_Analysis/
Module8_SSR_Analysis_TIMESTAMP/
Module9_Nucleotide_Diversity/ (if MAFFT installed and 2+ files)
```

### Example 2: SSR Mining Project

```bash
# Navigate to genome folder
cd ~/ssr_study

# Place GenBank files
ls *.gb
# species1.gb  species2.gb  species3.gb

# Run SSR analysis with standard thresholds
cgas --module 8

# Results saved to:
# Module8_SSR_Analysis_20251226_143315/
#   ├── SSR_Summary.xlsx          # Overview by type and region
#   ├── SSR_Detailed.xlsx         # Full list with positions
#   └── SSR_Statistics.xlsx       # Distribution stats
```

### Example 3: Comparative Study

```bash
# Project structure
mkdir ~/comparative_analysis
cd ~/comparative_analysis

# Add genomes from different species
cp ~/data/family1/*.gb .
cp ~/data/family2/*.gb .

# Run comparative modules
cgas --module 1  # Gene counts
cgas --module 3  # Genome structure comparison
cgas --module 4  # Codon usage patterns

# Results organized by module:
ls -d Module*_*/
```

### Example 4: Gene Structure Analysis

```bash
# Focus on gene architecture
cd ~/gene_structure_project

# Run relevant modules
cgas --module 1  # Gene inventory
cgas --module 2  # Gene tables
cgas --module 7  # Intron analysis

# Outputs:
# Module1_Gene_Count_Analysis/
#   └── Complete gene inventory with IR duplications
# Module2_Gene_Table/
#   └── Publication-ready Word document
# Module7_Intron_Analysis/
#   └── Intron/exon structures
```

### Example 5: Custom SSR Thresholds

```bash
# Different thresholds for different studies

# Study 1: Microsatellite markers (longer repeats)
cgas --module 8 -t 15,7,6,5,4,4

# Study 2: Population genetics (standard)
cgas --module 8 -t 10,5,4,3,3,3

# Study 3: Comprehensive survey (permissive)
cgas --module 8 -t 8,4,3,3,3,3
```

---

## Working with Different File Types

### GenBank Files (Modules 1-5, 7-9)

```bash
# Accepted extensions
ls *.gb *.gbk *.genbank *.gbff

# Verify format
head -n 5 species1.gb
# Should show: LOCUS, DEFINITION, ACCESSION...

# Run analysis
cgas
```

### FASTA Files (Module 6)

```bash
# Accepted extensions
ls *.fasta *.fa *.fna *.fas

# Must be aligned sequences
head -n 10 alignment.fasta

# Run SNP analysis
cgas --module 6
```

---

## Pipeline Integration

### Shell Script

```bash
#!/bin/bash
# analyze_chloroplasts.sh

DATA_DIR="./genomes"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="analysis_${TIMESTAMP}.log"

echo "Starting CGAS analysis at $(date)" | tee $LOG_FILE

# Navigate to data directory
cd "$DATA_DIR"

# Count input files
GB_COUNT=$(ls -1 *.gb 2>/dev/null | wc -l)
echo "Found $GB_COUNT GenBank files" | tee -a $LOG_FILE

# Run analysis
cgas 2>&1 | tee -a $LOG_FILE

echo "Analysis complete at $(date)" | tee -a $LOG_FILE
echo "Results saved in Module*_*/ folders" | tee -a $LOG_FILE
```

Run it:
```bash
chmod +x analyze_chloroplasts.sh
./analyze_chloroplasts.sh
```

### Makefile

```makefile
.PHONY: all count table compare codon ssr clean

all:
	@echo "Running complete CGAS analysis..."
	cgas

count:
	@echo "Running gene count analysis..."
	cgas --module 1

table:
	@echo "Generating gene tables..."
	cgas --module 2

compare:
	@echo "Running comparative analysis..."
	cgas --module 3

codon:
	@echo "Analyzing codon usage..."
	cgas --module 4

ssr:
	@echo "Running SSR analysis..."
	cgas --module 8 -t 10,5,4,3,3,3

clean:
	@echo "Cleaning output directories..."
	rm -rf Module*_*/

help:
	@echo "CGAS Makefile Commands:"
	@echo "  make all     - Run all modules"
	@echo "  make count   - Module 1: Gene count"
	@echo "  make table   - Module 2: Gene tables"
	@echo "  make compare - Module 3: Comparative analysis"
	@echo "  make codon   - Module 4: Codon usage"
	@echo "  make ssr     - Module 8: SSR analysis"
	@echo "  make clean   - Remove output directories"
```

Run it:
```bash
make all
make count
make ssr
```

---

## Logging and Monitoring

### Save Output to Log File

```bash
# Save all output with timestamp
cgas 2>&1 | tee "cgas_analysis_$(date +%Y%m%d_%H%M%S).log"

# Save only errors
cgas 2> errors.log

# Append to existing log
cgas 2>&1 | tee -a analysis.log
```

### Background Execution

```bash
# Run in background
nohup cgas > cgas.log 2>&1 &

# Check progress
tail -f cgas.log

# Find process
ps aux | grep cgas

# Kill if needed
pkill -f cgas
```

---

## Troubleshooting

### No GenBank Files Found

```bash
# Check current directory
pwd

# List GenBank files
ls *.gb *.gbk *.genbank 2>/dev/null

# If files are elsewhere
cd /path/to/genbank/files
cgas
```

### Module 9 MAFFT Error

```bash
# Check MAFFT installation
mafft --version

# Install MAFFT
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install mafft

# macOS
brew install mafft

# Windows
# Download from: https://mafft.cbrc.jp/alignment/software/
```

### Permission Denied

```bash
# Check permissions
ls -la

# Fix permissions
chmod +w .
chmod +x cgas

# Check output folder permissions
ls -la Module*_*/
```

### Import Errors

```bash
# Check installation
pip list | grep -E "biopython|pandas|openpyxl"

# Reinstall dependencies
pip install -r requirements.txt --upgrade

# Reinstall CGAS
cd ~/Chloroplast-Genome-Analysis-Suite-CGAS
pip install -e . --force-reinstall
```

### Output Folders Not Created

```bash
# Check write permissions in current directory
touch test_file && rm test_file

# Check disk space
df -h .

# Run with verbose output
cgas --module 1 -v
```

---

## Tips and Best Practices

### 1. Organize Your Data

```bash
# Good directory structure
project/
├── raw_data/
│   └── *.gb
├── analysis/
│   ├── Module1_Gene_Count_Analysis/
│   ├── Module2_Gene_Table/
│   └── ...
└── scripts/
```

### 2. Use Version Control

```bash
# Track your analysis
git init
git add *.gb
git commit -m "Initial genomes"

# After analysis
git add Module*_*/
git commit -m "Completed CGAS analysis $(date +%Y%m%d)"
```

### 3. Backup Results

```bash
# Create archive of results
timestamp=$(date +%Y%m%d_%H%M%S)
tar -czf cgas_results_${timestamp}.tar.gz Module*_*/

# Copy to backup location
cp cgas_results_${timestamp}.tar.gz ~/backups/
```

### 4. Document Your Workflow

```bash
# Create analysis log
cat > ANALYSIS_LOG.md << EOF
# CGAS Analysis Log

**Date:** $(date)
**Genomes:** $(ls *.gb | wc -l) files
**Modules Run:** 1-9

## Parameters:
- SSR thresholds: 10,5,4,3,3,3

## Results:
- See Module*_*/ folders

EOF
```

---

## Getting Help

```bash
# Main help
cgas --help

# List modules
cgas --list

# Check version
cgas --version

# Module-specific help
cgas --module 8 --help
```

---

## Update Package

```bash
# Update to latest version
cd ~/Chloroplast-Genome-Analysis-Suite-CGAS
git pull
pip install -e . --upgrade

# Verify update
cgas --version
```

---

For more information:
- [README.md](README.md) - Overview and features
- [INSTALL.md](INSTALL.md) - Installation instructions
- [GitHub Repository](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS)
- [Report Issues](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS/issues)
