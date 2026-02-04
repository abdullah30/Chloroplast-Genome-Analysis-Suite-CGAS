# CGAS Quick Start Guide
## Get Started in 5 Minutes! 🚀

---

## 🎯 What is CGAS?

**CGAS (Chloroplast Genome Analysis Suite)** is a comprehensive Python pipeline for comparative chloroplast genomics. It provides **14 integrated modules** covering everything from raw read assembly to phylogenetic tree construction.

---

## ⚡ Installation (Choose One Method)

### Option 1: Complete Installation (Recommended)
```bash
# Clone repository
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git
cd Chloroplast-Genome-Analysis-Suite-CGAS

# Create conda environment with all dependencies
conda env create -f environment.yml
conda activate cgas

# Install CGAS
pip install -e .
```

### Option 2: Minimal Installation (not include assembly and quality check of reads and assembly)
```bash
# For modules 2-14 only (no assembly and quality analysis)
conda env create -f environment-minimal.yml
conda activate cgas
pip install -e .
```

### Option 3: Manual Installation
```bash
# Install Python dependencies
pip install biopython pandas numpy openpyxl python-docx

# Install external tools (Ubuntu/Debian)
sudo apt-get install mafft

# Optional: For phylogeny and visualization
sudo apt-get install iqtree r-base
Rscript -e "install.packages(c('ggplot2', 'pheatmap', 'RColorBrewer', 'reshape2'))"

# Install CGAS
pip install git+https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git
```

---

## 📋 The 14 Modules at a Glance

### 🔧 **PHASE 1: Preparation Modules (1-4)**
Run these individually, as needed:

| Module | Command | Purpose |
|--------|---------|---------|
| **1** | `cgas-assembly` | Assemble genomes from raw FASTQ reads |
| **2** | `cgas-annotate` | Annotate assembled plastomes |
| **3** | `cgas-compare` | Compare genes across species |
| **4** | `cgas-convert` | Convert formats & prepare NCBI submissions |

### 🧬 **PHASE 2: Analysis Modules (5-14)**
Run these together on annotated genomes:

| Module | Command | Analysis |
|--------|---------|----------|
| **5** | `cgas-gene-compare` | Gene content comparison |
| **6** | `cgas-gene-table` | Publication-ready gene tables |
| **7** | `cgas-genome-compare` | Genome structure & GC content |
| **8** | `cgas-codon` | Codon usage (RSCU) + R plots |
| **9** | `cgas-amino` | Amino acid composition + R plots |
| **10** | `cgas-snp` | SNP/substitution analysis |
| **11** | `cgas-intron` | Intron characterization |
| **12** | `cgas-ssr` | SSR detection + R plots |
| **13** | `cgas-diversity` | Nucleotide diversity + R plots |
| **14** | `cgas-phylogeny` | Phylogenetic analysis + trees |

---

## 🚀 Quick Start Examples

### Example 1: Complete Analysis (Most Common)

**You have:** Annotated GenBank files (`.gb`)

```bash
# Put all your .gb files in a folder
cd my_genomes/

# Run all analysis modules at once
cgas --module 5-14

# Or run modules individually:
cgas-gene-compare          # Module 5
cgas-gene-table            # Module 6
cgas-genome-compare        # Module 7
cgas-codon                 # Module 8 (with R plots)
cgas-amino                 # Module 9 (with R plots)
cgas-snp                   # Module 10
cgas-intron                # Module 11
cgas-ssr                   # Module 12 (with R plots)
cgas-diversity             # Module 13 (with R plots)
cgas-phylogeny --iqtree    # Module 14 (with phylogenetic tree)
```

**Output:** Each module creates its own folder with publication-ready results!

### Example 2: Start from Raw Reads

**You have:** Raw sequencing data (`.fastq`)

```bash
# Step 1: Assemble genomes (Module 1)
cgas-assembly -i raw_reads/ -o assembled_genomes/

# Step 2: Annotate genomes (Module 2)
cgas-annotate -i assembled_genomes/ -r reference_genome.gb

# Step 3: Run all analyses (Modules 5-14)
cd annotated_genomes/
cgas --module 5-14
```

### Example 3: Phylogenetic Analysis Only

**You have:** GenBank files, need a phylogenetic tree

```bash
# Quick phylogeny with complete matrix (RECOMMENDED)
cgas-phylogeny --complete-only --iqtree -og outgroup.gb

# Or genes-only approach
cgas-phylogeny --genes-only --iqtree -og outgroup.gb

# With codon-aware alignment (MACSE)
cgas-phylogeny --genes-only --macse --iqtree -og outgroup.gb
```

**Output:** `Module14_Phylogeny/phylogeny_complete.treefile` + beautiful tree visualizations!

### Example 4: Just Nucleotide Diversity

**You have:** GenBank files, need π values and plots

```bash
cgas-diversity

# Or regenerate plots only (after running once)
cgas-diversity --figures-only
```

**Output:** 
- Text reports with π values
- Beautiful R-generated heatmaps and plots (PDF + PNG, 600 DPI)

### Example 5: Jupyter Notebook

```python
# In a Jupyter notebook cell:

# Run individual modules
%run cgas_module8.py      # Codon usage analysis
%run cgas_module13.py     # Nucleotide diversity

# Or use the ! operator
!cgas-codon
!cgas-diversity
!cgas-phylogeny --complete-only --iqtree
```

---

## 📁 Typical Workflow

### Complete Pipeline (Raw Reads → Publication)

```bash
# Start with raw reads in fastq_files/
mkdir analysis/
cd analysis/

# PHASE 1: Preparation
cgas-assembly -i ../fastq_files/ -o 1_assembled/
cgas-annotate -i 1_assembled/ -r reference.gb -o 2_annotated/

# PHASE 2: Analysis
cd 2_annotated/
cgas --module 5-14

# Results are in:
# - Module5_Gene_Comparison/
# - Module6_Gene_Table/
# - Module7_Genome_Structure/
# - Module8_Codon_Usage_Analysis/    (with R plots!)
# - Module9_Amino_Acid_Analysis/     (with R plots!)
# - Module10_SNP_Analysis/
# - Module11_Intron_Analysis/
# - Module12_SSR_Analysis/           (with R plots!)
# - Module13_Diversity_Analysis/     (with R plots!)
# - Module14_Phylogeny/              (with phylogenetic tree!)
```

---

## 💡 Pro Tips

### 1. **Check if R Visualization is Working**
```bash
# Test R installation
Rscript --version

# Test R packages
Rscript -e "library(ggplot2); library(pheatmap)"

# If packages are missing:
Rscript -e "install.packages(c('ggplot2', 'pheatmap', 'RColorBrewer', 'reshape2'), repos='https://cloud.r-project.org')"
```

### 2. **Run Modules Selectively**
```bash
# Just the visualization modules
cgas-codon          # Module 8
cgas-amino          # Module 9
cgas-ssr            # Module 12
cgas-diversity      # Module 13
cgas-phylogeny --iqtree  # Module 14
```

### 3. **Speed Up Analysis**
```bash
# Use multiple threads for phylogeny
cgas-phylogeny --iqtree --threads 16

# Skip figures if R isn't installed
cgas-codon --no-figures
cgas-diversity --no-figures
```

### 4. **Regenerate Plots Without Rerunning Analysis**
```bash
# For modules with --figures-only support
cgas-diversity --figures-only
```

---

## 🎨 What You Get (Output Highlights)

### Excel Files 📊
- **Gene tables** (Word + Excel formats)
- **Codon usage data** (RSCU values)
- **Amino acid frequencies**
- **SSR catalogs** (detailed classification)
- **Diversity statistics** (π values)

### Word Documents 📄
- **Publication-ready gene tables** (formatted, italicized gene names)

### R Visualizations 📈
Automatically generated when R is installed:
- **Codon usage heatmaps** (Module 8)
- **Amino acid composition plots** (Module 9)
- **SSR distribution charts** (Module 12)
- **Nucleotide diversity plots** (Module 13)
- **Phylogenetic trees** (Module 14)

All figures in **PDF** (vector) + **PNG** (600 DPI) formats!

---

## ❓ Common Questions

### Q: Do I need all 14 modules?
**A:** No! Use what you need:
- Just analyzing existing genomes? Start with modules 5-14
- Just want a phylogenetic tree? Use module 14 only
- Just need codon usage? Use module 8 only

### Q: What if I don't have R?
**A:** Most analyses still work! Excel/Word outputs are always generated. Only visualizations require R. You can:
- Skip R plots: `cgas-codon --no-figures`
- Or install R later and regenerate: `cgas-diversity --figures-only`

### Q: How long does it take?
**A:** 
- 10 genomes: ~10 minutes (modules 5-13)
- 50 genomes: ~30 minutes (complete analysis)
- Most time: nucleotide diversity + phylogeny (require alignment)

### Q: Can I run on Windows/Mac/Linux?
**A:** Yes! CGAS is fully cross-platform.

---

## 🆘 Getting Help

### Built-in Help
```bash
cgas --help                    # General help
cgas-phylogeny --help          # Module-specific help
cgas --module 14 --help        # Alternative syntax
```

### Documentation
- **GitHub:** https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS
- **README:** Check the main README.md for detailed documentation
- **Module READMEs:** Each module has its own documentation

### Test Installation
```bash
# Verify CGAS is installed
cgas --version

# Check dependencies
python -c "import Bio, pandas, numpy, openpyxl; print('Core dependencies OK')"
mafft --version
iqtree --version
Rscript --version
```

---

## 🎯 Next Steps

1. **Install CGAS** (choose method above)
2. **Prepare your data** (GenBank files in a folder)
3. **Run the modules** you need
4. **Check the output folders** for results
5. **Use the figures and tables** in your manuscript!

---

## 📚 Citation

If you use CGAS in your research, please cite:

```
Abdullah, Yan, R., & Tian, X. (2025). CGAS (Chloroplast Genome Analysis Suite): 
An Automated Python Pipeline for Comprehensive Comparative Chloroplast Genomics. 
bioRxiv. https://doi.org/10.64898/2025.12.21.695765
```

---

## 🌟 That's It!

You're ready to analyze chloroplast genomes! Start with the examples above, and explore the full documentation as needed.

**Happy analyzing! 🧬✨**
