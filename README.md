# Chloroplast Genome Analysis Suite (CGAS)

**CGAS: A comprehensive toolkit for analyzing chloroplast genomes with 9 specialized modules.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

## Features

✨ **Flexible Usage** - Unified script, module imports, individual scripts, Jupyter notebooks  
📊 **Publication-Ready** - Formatted Excel and Word outputs  
🧬 **Smart Handling** - Correctly handles trans-spliced genes (rps12)  
🖥️ **Cross-Platform** - Windows, Linux, macOS  
📚 **Well-Documented** - Comprehensive guides included  
🗂️ **Organized Outputs** - Each module creates its own output folder

## Installation

```bash
cd ~
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS
cd Chloroplast-Genome-Analysis-Suite-CGAS
pip install -r requirements.txt
pip install -e .
```

## Quick Start

```bash
# Run all modules
cgas

# Run specific module
cgas --module 1
cgas --module 2  # Similarly change the number 3 to 9

# List all modules
cgas --list

# For detailed help
cgas --help
```

## Modules

1. **Gene Count and Summary** - Comprehensive gene inventory with IR detection
2. **Gene Table Generation** - Publication-ready Word document tables
3. **Comparative Genome Analysis** - Multi-genome comparison with GC content
4. **Codon Usage Analysis (RSCU)** - Relative synonymous codon usage
5. **Amino Acid Composition** - Amino acid frequency analysis
6. **SNP/Substitution Analysis** - Variant detection from FASTA alignments
7. **Intron Extraction** - Intron/exon structure analysis
8. **SSR Analysis** - Simple sequence repeat detection and classification
9. **Nucleotide Diversity** - Pi and Theta calculations (requires MAFFT)

## Output Structure

Each module creates its own organized output directory:

```
your_working_directory/
├── input_file1.gb
├── input_file2.gb
├── input_file3.gb
├── Module1_Gene_Count_Analysis/
│   ├── Chloroplast_Gene_Analysis_20251226_143022.xlsx
│   └── Gene_Normalization_Log.xlsx
├── Module2_Gene_Table/
│   └── Gene_Table_20251226_143045.docx
├── Module3_Comparative_Analysis/
│   └── Comparative_Genome_Analysis.xlsx
├── Module4_Codon_Usage/
│   └── Codon_Usage_RSCU_20251226_143115.xlsx
├── Module5_Amino_Acid/
│   └── Amino_Acid_Composition_20251226_143145.xlsx
├── Module6_SNP_Analysis/
│   └── SNP_Analysis_20251226_143215.xlsx
├── Module7_Intron_Analysis/
│   └── Intron_Analysis_20251226_143245.xlsx
├── Module8_SSR_Analysis_20251226_143315/
│   ├── SSR_Summary.xlsx
│   ├── SSR_Detailed.xlsx
│   └── SSR_Statistics.xlsx
└── Module9_Nucleotide_Diversity/
    ├── Nucleotide_Diversity_Summary.xlsx
    └── alignment_files/
```

**Benefits of organized outputs:**
- Clean separation of inputs and outputs
- Easy to find results for each analysis
- No clutter in your working directory
- Version-controlled output with timestamps

## Usage Examples

### Command Line Interface

```bash
# Navigate to your data directory
cd /path/to/your/genbank/files

# Run all available modules
cgas

# Run only Module 1 (Gene Count)
cgas --module 1

# Run only Module 8 (SSR Analysis)
cgas --module 8

# Run with custom SSR thresholds
cgas --module 8 -t 12,6,5,4,4,4
```

### Jupyter Notebook

All CGAS modules can be executed in Jupyter Notebook using the `%run` magic command.

> ⚠️ **Important:** Always use `%run` for notebook execution.  
> Direct importing is not supported due to command-line argument parsing.

```python
# Recommended: Run unified analyzer
%run "unified_analyzer.py"

# Or run individual modules
%run "module1_gene_count.py"
%run "module2_gene_table.py"
%run "module3_comparative_analysis.py"
%run "module4_codon_usage.py"
%run "module5_amino_acid.py"
%run "module6_snp_analysis.py"
%run "module7_intron_extraction.py"
%run "module8_ssr_analysis.py" -t 10,5,4,3,3,3
%run "module9_diversity_analysis.py"
```

### Module 8: SSR Analysis Options

```bash
# Default thresholds
cgas --module 8

# Custom thresholds format: -t mono,di,tri,tetra,penta,hexa
cgas --module 8 -t 10,5,4,3,3,3

# Example: Stricter thresholds
cgas --module 8 -t 15,7,6,5,4,4
```

**SSR Threshold Meaning:**
- Mononucleotide repeats ≥ 10 (default)
- Dinucleotide repeats ≥ 5 (default)
- Trinucleotide repeats ≥ 4 (default)
- Tetranucleotide repeats ≥ 3 (default)
- Pentanucleotide repeats ≥ 3 (default)
- Hexanucleotide repeats ≥ 3 (default)

## Requirements

### Python Packages
```
Python 3.7+
biopython
pandas
openpyxl
numpy
python-docx
```

### External Tools (Optional)
- **MAFFT** - Required only for Module 9 (Nucleotide Diversity)
  ```bash
  # Ubuntu/Debian
  sudo apt-get install mafft
  
  # macOS
  brew install mafft
  
  # Windows
  # Download from https://mafft.cbrc.jp/alignment/software/
  ```

## Input File Requirements

### GenBank Files (Modules 1-5, 7-9)
- **Extensions:** `.gb`, `.gbk`, `.genbank`, `.gbff`
- **Required annotations:** Gene features with proper qualifiers
- **Format:** Standard GenBank format from NCBI

### FASTA Files (Module 6)
- **Extensions:** `.fasta`, `.fa`, `.fna`, `.fas`
- **Format:** Multi-FASTA alignment
- **Requirement:** Pre-aligned sequences

## Typical Workflow

```bash
# 1. Organize your data
mkdir ~/chloroplast_analysis
cd ~/chloroplast_analysis

# 2. Copy or download GenBank files
cp /path/to/*.gb .
# or download from NCBI

# 3. Run complete analysis
cgas

# 4. Results are organized in separate folders
ls -d Module*_*/
```

## Performance

On a standard workstation (Intel Core i5, 16GB RAM):
- **50 genomes:** ~5 minutes for complete analysis
- **10 genomes:** ~1 minute for complete analysis
- **Single genome:** ~10-15 seconds per module

## Output Files

### Excel Workbooks (.xlsx)
- **Module 1:** Gene count summary with normalization log
- **Module 3:** Comparative genome table with GC content
- **Module 4:** RSCU values for all codons
- **Module 5:** Amino acid frequencies
- **Module 6:** SNP/substitution matrices
- **Module 7:** Intron structure data
- **Module 8:** SSR summary, detailed list, and statistics
- **Module 9:** Nucleotide diversity (Pi, Theta)

### Word Documents (.docx)
- **Module 2:** Publication-ready gene tables

### Features of Output Files
- ✅ Professional formatting
- ✅ Color-coded headers
- ✅ Frozen header rows
- ✅ Auto-sized columns
- ✅ Comprehensive footnotes
- ✅ Timestamped filenames

## Troubleshooting

### No GenBank files found
```bash
# Check file extensions
ls *.gb *.gbk *.genbank

# Ensure files are in current directory
pwd
```

### Module 9 fails
```bash
# Check MAFFT installation
mafft --version

# Install if missing
sudo apt-get install mafft  # Linux
brew install mafft          # macOS
```

### Permission errors
```bash
# Ensure write permissions
chmod +w .

# Check folder permissions
ls -la
```

### Import errors
```bash
# Reinstall dependencies
pip install -r requirements.txt --upgrade
```

## Citation

If you use CGAS in your research, please cite:

**APA Style:**
```
Abdullah, Yan, R., & Tian, X. (2025). CGAS (Chloroplast Genome Analysis Suite): 
    An automated Python pipeline for comprehensive comparative chloroplast genomics. 
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

**Plain Text:**
```
Abdullah, Rushan Yan, Xiaoxuan Tian (2025). CGAS (Chloroplast Genome Analysis Suite): 
An Automated Python Pipeline for Comprehensive Comparative Chloroplast Genomics. 
bioRxiv. doi: https://doi.org/10.64898/2025.12.21.695765
```

## Contributing

Contributions are welcome! Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Abdullah**  
📧 [Contact via GitHub Issues](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS/issues)  
🌐 [GitHub Profile](https://github.com/abdullah30)

## Acknowledgments

- BioPython community for excellent genomic data handling
- NCBI for providing comprehensive GenBank annotations
- Scientific community for chloroplast genome research

## Version History

- **1.0.0** (December 2025) - Initial release
  - 9 comprehensive analysis modules
  - Organized output folder structure
  - Publication-ready outputs
  - Cross-platform compatibility

## Support

- 📖 [Documentation](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS)
- 🐛 [Report Issues](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS/issues)
- 💬 [Discussions](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS/discussions)

---

**Happy Analyzing! 🧬🔬📊**
