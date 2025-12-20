# Command-Line Usage Guide

## Installation

First, install the package:

```bash
pip install git+https://github.com/yourusername/chloroplast-analyzer.git
```

## Basic Usage

### Run All Modules

```bash
# Navigate to directory with your GenBank/FASTA files
cd /path/to/your/data

# Run all available modules
chloroplast-analyze
```

### Run Specific Module

```bash
# Run only Module 1 (Gene Count)
chloroplast-analyze --module 1

# Run only Module 4 (Codon Usage)
chloroplast-analyze --module 4
```

### Run Multiple Modules

```bash
# Run modules 1, 4, and 8
chloroplast-analyze --modules 1,4,8

# Run modules 1 through 5
chloroplast-analyze --modules 1,2,3,4,5
```

## Individual Module Commands

Each module has its own command for convenience:

### Module 1: Gene Count and Summary
```bash
chloroplast-gene-count
# or
cp-count
```

### Module 2: Gene Table Generation
```bash
chloroplast-gene-table
# or
cp-table
```

### Module 3: Comparative Genome Analysis
```bash
chloroplast-compare
# or
cp-compare
```

### Module 4: Codon Usage Analysis
```bash
chloroplast-codon
# or
cp-codon
```

### Module 5: Amino Acid Composition
```bash
chloroplast-aminoacid
# or
cp-aa
```

### Module 6: SNP/Substitution Analysis
```bash
chloroplast-snp
# or
cp-snp
```

### Module 7: Intron Extraction
```bash
chloroplast-intron
# or
cp-intron
```

### Module 8: SSR Analysis
```bash
# Default thresholds
chloroplast-ssr

# Custom thresholds
chloroplast-ssr --mono 15 --di 7 --tri 5 --tetra 4

# or short version
cp-ssr --mono 15 --di 7
```

### Module 9: Nucleotide Diversity
```bash
chloroplast-diversity
# or
cp-diversity

# Note: Requires MAFFT and 2+ GenBank files
```

## Advanced Options

### Specify Input Directory

```bash
# Analyze files in a specific directory
chloroplast-analyze --input /path/to/data

# Short version
chloroplast-analyze -i /path/to/data
```

### Specify Output Directory

```bash
# Save results to specific directory
chloroplast-analyze --output /path/to/results

# Short version
chloroplast-analyze -o /path/to/results
```

### Combined Options

```bash
# Run specific modules with custom directories
chloroplast-analyze --modules 1,4,8 --input ./data --output ./results

# Short version
chloroplast-analyze -m 4 -i ./genomes -o ./outputs
```

### List Available Modules

```bash
# Show all modules with descriptions
chloroplast-analyze --list
```

### Show Version

```bash
chloroplast-analyze --version
```

### Verbose Output

```bash
chloroplast-analyze --verbose
# or
chloroplast-analyze -v
```

## Practical Examples

### Example 1: Quick Analysis

```bash
# You have GenBank files in current directory
ls *.gb
# species1.gb  species2.gb  species3.gb

# Run all analyses
chloroplast-analyze
```

### Example 2: Gene Count Only

```bash
# Quick gene count
cd my_genomes/
chloroplast-gene-count
```

### Example 3: Codon Usage Analysis

```bash
# Analyze codon usage for chloroplast genomes
cd chloroplast_data/
chloroplast-codon
```

### Example 4: SSR Analysis with Custom Thresholds

```bash
# Run SSR with stricter thresholds
cd genomes/
chloroplast-ssr --mono 20 --di 10 --tri 6 --tetra 5 --penta 4 --hexa 4
```

### Example 5: Multiple Genomes Comparison

```bash
# Compare multiple species
cd comparative_study/
chloroplast-compare
```

### Example 6: Complete Pipeline

```bash
# Step 1: Organize data
mkdir -p ~/analysis/chloroplast_study
cd ~/analysis/chloroplast_study

# Step 2: Copy GenBank files
cp ~/data/genomes/*.gb .

# Step 3: Run specific modules
chloroplast-analyze --modules 1,3,4,8

# Step 4: Check results
ls *.xlsx
```

### Example 7: Batch Processing

```bash
# Process multiple datasets
for dataset in dataset1 dataset2 dataset3; do
    echo "Processing $dataset..."
    cd $dataset
    chloroplast-analyze
    cd ..
done
```

### Example 8: Analysis with Different Input Types

```bash
# Analyze GenBank files
cd genbank_files/
chloroplast-gene-count

# Analyze FASTA files (for SNP analysis)
cd ../fasta_files/
chloroplast-snp
```

## Module-Specific Options

### SSR Analysis Options

```bash
chloroplast-ssr [OPTIONS]

Options:
  -i, --input DIR       Input directory (default: current)
  -o, --output DIR      Output directory
  --mono INT           Mononucleotide threshold (default: 10)
  --di INT             Dinucleotide threshold (default: 5)
  --tri INT            Trinucleotide threshold (default: 4)
  --tetra INT          Tetranucleotide threshold (default: 3)
  --penta INT          Pentanucleotide threshold (default: 3)
  --hexa INT           Hexanucleotide threshold (default: 3)

Examples:
  chloroplast-ssr
  chloroplast-ssr --mono 15 --di 7
  chloroplast-ssr -i ./data -o ./results --tri 6
```

## Working with Different File Types

### GenBank Files

Most modules work with GenBank files:

```bash
# Ensure files have correct extensions
ls *.gb *.gbk *.genbank *.gbff

# Run analysis
chloroplast-analyze
```

### FASTA Files

Module 6 (SNP Analysis) uses FASTA files:

```bash
# Check FASTA files
ls *.fasta *.fa *.fna *.fas

# Run SNP analysis
chloroplast-snp
```

## Output Files

### Default Output Locations

```bash
# Excel files in current directory
Chloroplast_Gene_Count_YYYYMMDD_HHMMSS.xlsx
Codon_Usage_RSCU_YYYYMMDD_HHMMSS.xlsx

# Module-specific folders
Module8_SSR_Analysis_YYYYMMDD_HHMMSS/
Module9_Nucleotide_Diversity/
```

### Custom Output Location

```bash
# Specify output directory
chloroplast-analyze --output ./results

# Results will be saved to ./results/
```

## Pipeline Integration

### Shell Script

```bash
#!/bin/bash
# analyze_chloroplasts.sh

DATA_DIR="./genomes"
OUTPUT_DIR="./results_$(date +%Y%m%d)"

mkdir -p "$OUTPUT_DIR"

echo "Starting analysis..."
chloroplast-analyze --input "$DATA_DIR" --output "$OUTPUT_DIR"

echo "Analysis complete!"
echo "Results saved to: $OUTPUT_DIR"
```

Run it:
```bash
chmod +x analyze_chloroplasts.sh
./analyze_chloroplasts.sh
```

### Makefile

```makefile
.PHONY: all count codon ssr clean

all:
	chloroplast-analyze

count:
	chloroplast-gene-count

codon:
	chloroplast-codon

ssr:
	chloroplast-ssr --mono 15 --di 7

clean:
	rm -f *.xlsx *.docx
	rm -rf Module8_SSR_Analysis_* Module9_Nucleotide_Diversity
```

Run it:
```bash
make all
make count
make ssr
```

## Logging and Monitoring

### Save Output to Log File

```bash
# Save all output
chloroplast-analyze 2>&1 | tee analysis.log

# Save only errors
chloroplast-analyze 2> errors.log

# Save with timestamp
chloroplast-analyze 2>&1 | tee "analysis_$(date +%Y%m%d_%H%M%S).log"
```

### Background Execution

```bash
# Run in background
nohup chloroplast-analyze > analysis.log 2>&1 &

# Check progress
tail -f analysis.log

# Find process
ps aux | grep chloroplast
```

## Tips and Tricks

### 1. Quick Module Check

```bash
# List modules before running
chloroplast-analyze --list
```

### 2. Test with Single Module

```bash
# Test setup with quick module
chloroplast-gene-count
```

### 3. Organize Results

```bash
# Create timestamped result folder
RESULT_DIR="results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULT_DIR"
chloroplast-analyze --output "$RESULT_DIR"
```

### 4. Process Multiple Directories

```bash
# Process all subdirectories
for dir in */; do
    echo "Processing $dir"
    (cd "$dir" && chloroplast-analyze)
done
```

### 5. Conditional Module Execution

```bash
# Only run if GenBank files exist
if ls *.gb 1> /dev/null 2>&1; then
    chloroplast-analyze
else
    echo "No GenBank files found"
fi
```

## Getting Help

```bash
# Main help
chloroplast-analyze --help

# Module help
chloroplast-ssr --help

# List all commands
compgen -c | grep chloroplast

# Or
compgen -c | grep cp-
```

## Troubleshooting Commands

```bash
# Check installation
which chloroplast-analyze
chloroplast-analyze --version

# Test import
python -c "import chloroplast_analyzer; print(chloroplast_analyzer.__version__)"

# Check files
chloroplast-analyze --input . --list

# Verify MAFFT (for Module 9)
mafft --version
```

## Update Package

```bash
# Update to latest version
pip install --upgrade git+https://github.com/yourusername/chloroplast-analyzer.git

# Or reinstall
pip uninstall chloroplast-analyzer
pip install git+https://github.com/yourusername/chloroplast-analyzer.git
```

---

For more information, see:
- [README.md](README.md) - Overview and features
- [INSTALL.md](INSTALL.md) - Installation instructions
- GitHub: https://github.com/yourusername/chloroplast-analyzer
