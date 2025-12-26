# Quick Start Guide

## 🚀 Install in 30 Seconds

```bash
pip install git+https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git
```

## ✅ Verify Installation

```bash
cgas --version
cgas --list
```

## 📊 Run Analysis

### All Modules
```bash
cd /path/to/your/genbank/files
cgas
```

### Single Module
```bash
cgas --module 1
```

### Multiple Modules
```bash
cgas --modules 1,4,8
```

## 🎯 Quick Commands

| Command | Module | What it does |
|---------|--------|--------------|
| `cgas-count` | 1 | Count genes, tRNAs, rRNAs |
| `cgas-table` | 2 | Generate gene tables (Word) |
| `cgas-compare` | 3 | Compare genomes |
| `cgas-codon` | 4 | Codon usage (RSCU) |
| `cgas-aa` | 5 | Amino acid composition |
| `cgas-snp` | 6 | SNP analysis |
| `cgas-intron` | 7 | Extract introns |
| `cgas-ssr` | 8 | Find SSRs |
| `cgas-diversity` | 9 | Nucleotide diversity |

## 💡 Examples

### Example 1: Gene Count
```bash
cd my_genomes/
cgas-count
# Output: Chloroplast_Gene_Count_20241220.xlsx
```

### Example 2: Codon Usage
```bash
cd my_genomes/
cgas-codon
# Output: Codon_Usage_RSCU_20241220.xlsx
```

### Example 3: SSR with Custom Settings
```bash
cd my_genomes/
cgas-ssr --mono 15 --di 7 --tri 5
# Output: Module8_SSR_Analysis_20241220/
```

### Example 4: Multiple Modules
```bash
cd my_genomes/
cgas --modules 1,3,4,8
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
cgas --input ./data

# Specify output directory  
cgas --output ./results

# Run specific module
cgas --module 4

# Run multiple modules
cgas --modules 1,4,8

# List all modules
cgas --list

# Show version
cgas --version

# Help
cgas --help
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
- **GitHub**: https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS

---

**Need help?** Open an issue on GitHub or check the full documentation!
