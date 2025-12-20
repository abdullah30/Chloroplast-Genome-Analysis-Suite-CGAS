# Quick Start Guide

## 🚀 Install in 30 Seconds

```bash
pip install git+https://github.com/yourusername/chloroplast-analyzer.git
```

## ✅ Verify Installation

```bash
chloroplast-analyze --version
chloroplast-analyze --list
```

## 📊 Run Analysis

### All Modules
```bash
cd /path/to/your/genbank/files
chloroplast-analyze
```

### Single Module
```bash
chloroplast-analyze --module 1
```

### Multiple Modules
```bash
chloroplast-analyze --modules 1,4,8
```

## 🎯 Quick Commands

| Command | Module | What it does |
|---------|--------|--------------|
| `chloroplast-gene-count` | 1 | Count genes, tRNAs, rRNAs |
| `chloroplast-gene-table` | 2 | Generate gene tables (Word) |
| `chloroplast-compare` | 3 | Compare genomes |
| `chloroplast-codon` | 4 | Codon usage (RSCU) |
| `chloroplast-aminoacid` | 5 | Amino acid composition |
| `chloroplast-snp` | 6 | SNP analysis |
| `chloroplast-intron` | 7 | Extract introns |
| `chloroplast-ssr` | 8 | Find SSRs |
| `chloroplast-diversity` | 9 | Nucleotide diversity |

## 💡 Examples

### Example 1: Gene Count
```bash
cd my_genomes/
chloroplast-gene-count
# Output: Chloroplast_Gene_Count_20241220.xlsx
```

### Example 2: Codon Usage
```bash
cd my_genomes/
chloroplast-codon
# Output: Codon_Usage_RSCU_20241220.xlsx
```

### Example 3: SSR with Custom Settings
```bash
cd my_genomes/
chloroplast-ssr --mono 15 --di 7 --tri 5
# Output: Module8_SSR_Analysis_20241220/
```

### Example 4: Multiple Modules
```bash
cd my_genomes/
chloroplast-analyze --modules 1,3,4,8
```

## 📁 File Requirements

### GenBank Files (.gb, .gbk, .genbank)
- For modules 1, 2, 3, 4, 5, 7, 8, 9
- Must have CDS, tRNA, rRNA annotations

### FASTA Files (.fasta, .fa, .fna)
- For module 6 (SNP analysis)

## 🔧 Options

```bash
# Specify input directory
chloroplast-analyze --input ./data

# Specify output directory  
chloroplast-analyze --output ./results

# Run specific module
chloroplast-analyze --module 4

# Run multiple modules
chloroplast-analyze --modules 1,4,8

# List all modules
chloroplast-analyze --list

# Show version
chloroplast-analyze --version

# Help
chloroplast-analyze --help
```

## 📋 System Requirements

- Python 3.7+
- BioPython, Pandas, OpenPyXL, NumPy, python-docx
- MAFFT (only for Module 9)

### Install MAFFT
```bash
# Ubuntu/Debian
sudo apt-get install mafft

# macOS
brew install mafft
```

## 🆘 Troubleshooting

### Command not found?
```bash
export PATH="$HOME/.local/bin:$PATH"
```

### Import errors?
```bash
pip install biopython pandas openpyxl numpy python-docx
```

### No GenBank files found?
```bash
# Check file extensions
ls *.gb *.gbk
# Files must end in .gb, .gbk, .genbank, or .gbff
```

## 📚 More Information

- **Full Documentation**: [README.md](README.md)
- **Installation Guide**: [INSTALL.md](INSTALL.md)
- **Command Reference**: [USAGE.md](USAGE.md)
- **GitHub**: https://github.com/yourusername/chloroplast-analyzer

---

**Need help?** Open an issue on GitHub or check the full documentation!
