cd ~/chloroplast-analyzer

cat > README.md << 'EOF'
# Chloroplast Genome Analyzer

A comprehensive toolkit for analyzing chloroplast genomes with 9 specialized modules.

## Installation
```bash
pip install -e .
```

## Usage
```bash
# Run all modules
chloroplast-analyze

# Run specific module
chloroplast-analyze --module 1

# List all modules
chloroplast-analyze --list
```

## Modules

1. Gene Count and Summary
2. Gene Table Generation
3. Comparative Genome Analysis
4. Codon Usage Analysis (RSCU)
5. Amino Acid Composition
6. SNP/Substitution Analysis
7. Intron Extraction
8. SSR Analysis
9. Nucleotide Diversity

## Requirements

- Python 3.7+
- biopython
- pandas
- openpyxl
- numpy
- python-docx

## Author

Abdullah - 2024
EOF

pip install -e .