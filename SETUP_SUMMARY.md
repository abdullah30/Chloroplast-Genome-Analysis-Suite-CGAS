# CGAS Setup Summary

> **Quick reference guide for installing and configuring CGAS**

This document provides a condensed overview of CGAS setup for quick reference. For detailed instructions, see [INSTALL.md](INSTALL.md).

---

## Table of Contents

- [System Requirements](#system-requirements)
- [Quick Installation](#quick-installation)
- [Installation Methods Comparison](#installation-methods-comparison)
- [External Tools by Module](#external-tools-by-module)
- [Platform-Specific Quick Setup](#platform-specific-quick-setup)
- [Verification Checklist](#verification-checklist)
- [Common Issues Quick Fixes](#common-issues-quick-fixes)
- [Environment Activation](#environment-activation)
- [First Steps After Installation](#first-steps-after-installation)

---

## System Requirements

### Minimum Setup
- **Python**: 3.9 or higher
- **RAM**: 4GB (8GB+ recommended)
- **Disk**: 5GB for software + data space
- **OS**: Linux, macOS 10.14+, or Windows 10+

### For Complete Functionality
- **Python**: 3.9 or higher
- **RAM**: 16GB+
- **CPU**: 8+ cores
- **Disk**: 200GB+

---

## Quick Installation

### Option 1: Complete Setup (Recommended for Servers)

**Everything included - Python, tools, R, packages**

```bash
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas
cd cgas
conda env create -f environment.yml
conda activate cgas
pip install -e .
cgas --list
```

**Time**: ~15-30 minutes  
**Modules**: All (1-14)

---

### Option 2: Minimal Setup (Recommended for low specification Laptops)
#I used the completed version with on Corei5 i13 inspiron 15 3530 having 32 GB RAM. So, the minimal version will require less resources. 
**Lightweight - core functionality only wihout**

```bash
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas
cd cgas
conda env create -f environment-minimal.yml
conda activate cgas-minimal
pip install -e .
cgas --list
```

**Time**: ~5-10 minutes  
**Modules**: 3-14 (add tools for 1-2 as needed)

---

### Option 3: Manual Setup (Maximum Control)

**Custom installation without conda**

```bash
# Install Python 3.9+ first, then:
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas
cd cgas
python3.9 -m venv cgas_env
source cgas_env/bin/activate  # Windows: cgas_env\Scripts\activate
pip install -e .
pip install -r requirements.txt
cgas --list
```

**Time**: ~5 minutes + external tools  
**Modules**: Python parts work; install external tools separately

---

## Installation Methods Comparison

| Feature | Complete Conda | Minimal Conda | Manual |
|---------|----------------|---------------|--------|
| **Setup Time** | 15-30 min | 5-10 min | 5 min + tools |
| **Disk Space** | ~5GB | ~2GB | ~1GB + tools |
| **All Modules** | ✓ | Partial | Partial |
| **Python Deps** | ✓ | ✓ | ✓ |
| **External Tools** | ✓ | MAFFT, IQ-TREE | Manual install |
| **R + Packages** | ✓ | ✗ | Manual install |
| **GetOrganelle** | ✓ | ✗ | Manual install |
| **Best For** | Servers, complete workflow | Laptops, modules 3-14 | Custom setups |

---

## External Tools by Module

### Required Tools Matrix

| Module | Tool | Install Priority | Installation Method |
|--------|------|------------------|---------------------|
| **1** | fastp | High | `apt-get`, `brew`, or `pip install fastp` |
| **1** | GetOrganelle | High | `pip install getorganelle` |
| **1** | BWA | High | `apt-get`, `brew`, or conda |
| **1** | SAMtools | High | `apt-get`, `brew`, or conda |
| **2** | BLAST+ | High | `apt-get`, `brew`, or conda |
| **2** | PGA | High | Manual (see below) |
| **3-7** | None | - | Python only |
| **8, 9** | R + packages | Medium | See R Setup |
| **10-11** | None | - | Python only |
| **12, 13** | R + packages | Medium | See R Setup |
| **14** | MAFFT | High | `apt-get`, `brew`, or conda |
| **14** | IQ-TREE | High | `apt-get`, `brew`, or conda |
| **14** | MACSE | Low (optional) | Manual + Java |

### Quick Install Commands

**Ubuntu/Debian:**
```bash
# Core tools (modules 3-14)
sudo apt-get install mafft ncbi-blast+ r-base iqtree

# Assembly tools (module 1)
sudo apt-get install fastp samtools bwa
pip install getorganelle
get_organelle_config.py --add embplant_pt
```

**macOS:**
```bash
# Core tools
brew install mafft blast iqtree r

# Assembly tools
brew install fastp samtools bwa
pip install getorganelle
get_organelle_config.py --add embplant_pt
```

**Windows (Chocolatey):**
```powershell
choco install mafft r.project
# Additional tools may require manual installation
```

---

## Platform-Specific Quick Setup

### Ubuntu/Debian One-Liner

```bash
sudo apt-get update && \
sudo apt-get install -y python3.9 python3.9-venv git mafft ncbi-blast+ r-base iqtree fastp samtools bwa && \
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas && \
cd cgas && \
conda env create -f environment.yml && \
conda activate cgas && \
pip install -e . && \
echo "CGAS installed successfully!"
```

### macOS (Homebrew) One-Liner

```bash
brew install python@3.9 git mafft blast iqtree r fastp samtools bwa && \
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas && \
cd cgas && \
python3.9 -m venv cgas_env && \
source cgas_env/bin/activate && \
pip install -e . && \
echo "CGAS installed successfully!"
```

---

## R Setup (for Visualization Modules)

### Quick R Package Installation

```bash
# Install all required R packages at once
R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'readr', 'writexl', 'seqinr', 'patchwork', 'RColorBrewer', 'ggrepel', 'scales', 'cowplot', 'gridExtra'), repos='https://cran.rstudio.com/')"
```

### Module-Specific R Packages

| Modules | Required Packages |
|---------|-------------------|
| **8** | ggplot2, seqinr, dplyr, tidyr, RColorBrewer, patchwork |
| **9** | ggplot2, dplyr, tidyr, scales |
| **12** | ggplot2, dplyr, tidyr, RColorBrewer, patchwork, ggrepel, scales |
| **13** | ggplot2, dplyr, tidyr, cowplot, gridExtra |

### Install by Module

```bash
# Module 8 only
R -e "install.packages(c('ggplot2', 'seqinr', 'dplyr', 'tidyr', 'RColorBrewer', 'patchwork'), repos='https://cran.rstudio.com/')"

# Module 9 only
R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'scales'), repos='https://cran.rstudio.com/')"

# Module 12 only
R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'RColorBrewer', 'patchwork', 'ggrepel', 'scales'), repos='https://cran.rstudio.com/')"

# Module 13 only
R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'cowplot', 'gridExtra'), repos='https://cran.rstudio.com/')"
```

---

## PGA Installation (Module 2)

### Quick Install

```bash
# Clone PGA
git clone https://github.com/quxiaojian/PGA.git ~/tools/PGA
cd ~/tools/PGA

# Install Perl dependencies (Ubuntu/Debian)
sudo apt-get install perl libdbi-perl libdbd-mysql-perl libxml-simple-perl bioperl

# Make executable
chmod +x PGA.pl

# Test installation
perl PGA.pl

# Use with CGAS
cgas --module 2 -i genomes/ -r reference/ --pga ~/tools/PGA/PGA.pl
```

---

## MACSE Installation (Module 14 - Optional)

### Quick Install

```bash
# Install Java (if needed)
sudo apt-get install default-jre  # Ubuntu
brew install openjdk              # macOS

# Download MACSE
mkdir -p ~/tools && cd ~/tools
wget https://bioweb.supagro.inra.fr/macse/releases/macse_v2.07.jar

# Create executable wrapper
cat > ~/tools/macse << 'EOF'
#!/bin/bash
java -Xmx4G -jar ~/tools/macse_v2.07.jar "$@"
EOF

chmod +x ~/tools/macse

# Add to PATH
echo 'export PATH="$HOME/tools:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Test
macse -help
```

---

## Verification Checklist

### ✓ Core Installation

```bash
# Check CGAS
cgas --version           # Should show v1.0.1
cgas --list              # Should list modules 1-14

# Check Python packages
python -c "import Bio; import pandas; import openpyxl; print('Python packages OK')"
```

### ✓ External Tools

```bash
# Quick verification
which python && echo "✓ Python"
which mafft && echo "✓ MAFFT"
which blastn && echo "✓ BLAST+"
which iqtree && echo "✓ IQ-TREE" || which iqtree2 && echo "✓ IQ-TREE2"
which R && echo "✓ R"

# For Module 1
which fastp && echo "✓ fastp"
which get_organelle_from_reads.py && echo "✓ GetOrganelle"
which bwa && echo "✓ BWA"
which samtools && echo "✓ SAMtools"

# For Module 14 (optional)
which java && echo "✓ Java"
which macse && echo "✓ MACSE" || echo "⚠ MACSE not found (optional)"
```

### ✓ R Packages

```bash
R -e "
packages <- c('ggplot2', 'dplyr', 'tidyr', 'seqinr')
for (pkg in packages) {
  if (require(pkg, character.only=TRUE, quietly=TRUE)) {
    cat(sprintf('✓ %s\n', pkg))
  } else {
    cat(sprintf('✗ %s MISSING\n', pkg))
  }
}
"
```

### ✓ Test Run

```bash
# Create test directory
mkdir -p ~/cgas_test && cd ~/cgas_test

# Test Module 5 (requires GenBank files)
# Copy some GenBank files here, then:
cgas --module 5

# Check output
ls -R Module5_Gene_Comparative_Analysis/
```

---

## Common Issues Quick Fixes

### Issue: "cgas: command not found"

```bash
# Activate environment
conda activate cgas           # If using conda
source cgas_env/bin/activate  # If using venv

# Or reinstall
cd /path/to/cgas
pip install -e .
```

### Issue: "Module X requires tool Y"

```bash
# Install missing tool
sudo apt-get install <tool>   # Ubuntu
brew install <tool>            # macOS
conda install -c bioconda <tool>  # Any platform with conda
```

### Issue: "Python package not found"

```bash
# Reinstall requirements
pip install -r requirements.txt

# Or specific package
pip install biopython pandas openpyxl
```

### Issue: "R package missing"

```bash
# Install missing package
R -e "install.packages('PACKAGE_NAME', repos='https://cran.rstudio.com/')"
```

### Issue: "Permission denied"

```bash
# Fix ownership
sudo chown -R $USER:$USER ~/.local

# Or use --user flag
pip install --user -e .
```

### Issue: GetOrganelle database not found

```bash
# Download database
get_organelle_config.py --add embplant_pt

# Verify
get_organelle_config.py --list
```

---

## Environment Activation

### Every Time You Use CGAS

**If using conda:**
```bash
conda activate cgas
# or
conda activate cgas-minimal
```

**If using virtual environment:**
```bash
source cgas_env/bin/activate      # Linux/macOS
cgas_env\Scripts\activate         # Windows
```

**Verify activation:**
```bash
which cgas    # Should show path in your environment
cgas --list   # Should work without errors
```

**Deactivation:**
```bash
conda deactivate          # For conda
deactivate                # For venv
```

---

## First Steps After Installation

### 1. Verify Installation

```bash
cgas --version
cgas --list
```

### 2. Prepare Your Data

```bash
# Create project structure
mkdir -p ~/my_cgas_project/{raw_reads,genomes,reference,results}
cd ~/my_cgas_project
```

### 3. Check Available Tools

```bash
# Run verification script
cat > check_setup.sh << 'EOF'
#!/bin/bash
echo "=== CGAS Setup Verification ==="
echo ""
echo "CGAS:"
cgas --version 2>/dev/null && echo "✓" || echo "✗ NOT FOUND"
echo ""
echo "External Tools:"
which mafft >/dev/null 2>&1 && echo "✓ MAFFT" || echo "✗ MAFFT"
which blastn >/dev/null 2>&1 && echo "✓ BLAST+" || echo "✗ BLAST+"
which iqtree >/dev/null 2>&1 && echo "✓ IQ-TREE" || which iqtree2 >/dev/null 2>&1 && echo "✓ IQ-TREE2" || echo "✗ IQ-TREE"
which R >/dev/null 2>&1 && echo "✓ R" || echo "✗ R"
which fastp >/dev/null 2>&1 && echo "✓ fastp" || echo "⚠ fastp (Module 1)"
which get_organelle_from_reads.py >/dev/null 2>&1 && echo "✓ GetOrganelle" || echo "⚠ GetOrganelle (Module 1)"
echo ""
echo "=== Setup complete! ==="
EOF

chmod +x check_setup.sh
./check_setup.sh
```

### 4. Read Documentation

```bash
# View help
cgas --help

# Check module-specific help
cgas --module 1 --help
cgas --module 5 --help
```

### 5. Run Your First Analysis

**If you have GenBank files:**
```bash
cd genomes_directory/
cgas --module 5
```

**If you have raw reads:**
```bash
cgas --module 1 -i raw_reads/ -o results/
```

---

## Quick Reference Card

### Essential Commands

```bash
# Help and info
cgas --help                              # General help
cgas --list                              # List all modules
cgas --module X --help                   # Module-specific help
cgas --version                           # Show version

# Running modules
cgas --module 5                          # Single module
cgas --modules 5,6,7                     # Multiple modules
cgas --module 1 -i input/ -o output/     # With options

# Environment
conda activate cgas                      # Activate (conda)
source cgas_env/bin/activate            # Activate (venv)
conda deactivate / deactivate           # Deactivate
```

### File Requirements

| Module | Input Format | Output Format |
|--------|-------------|---------------|
| **1** | FASTQ (.fastq.gz) | FASTA + reports |
| **2** | FASTA (.fasta, .fa) | GenBank (.gb) |
| **3-14** | GenBank (.gb) | Various (Excel, PDF, etc.) |

### Common Paths to Remember

```bash
# Installation
~/cgas/                                  # CGAS installation directory
~/cgas_env/                              # Virtual environment (if used)

# External tools (examples)
~/tools/PGA/PGA.pl                       # PGA location
~/tools/macse                            # MACSE wrapper
~/.GetOrganelle/                         # GetOrganelle database

# Environment files
~/.bashrc or ~/.zshrc                    # Shell configuration
```

---

## Module Availability Matrix

| Module | Python Only | Needs R | Needs External Tools | Works in Jupyter |
|--------|-------------|---------|---------------------|------------------|
| **1** | ✗ | ✗ | fastp, GetOrganelle, BWA, SAMtools | ⚠ (Linux/macOS) |
| **2** | ✗ | ✗ | BLAST+, PGA | ⚠ (Linux/macOS) |
| **3** | ✓ | ✗ | ✗ | ✓ |
| **4** | ✓ | ✗ | ✗ | ✓ |
| **5** | ✓ | ✗ | ✗ | ✓ |
| **6** | ✓ | ✗ | ✗ | ✓ |
| **7** | ✓ | ✗ | ✗ | ✓ |
| **8** | ✗ | ✓ | ✗ | ✓ |
| **9** | ✗ | ✓ | ✗ | ✓ |
| **10** | ✓ | ✗ | ✗ | ✓ |
| **11** | ✓ | ✗ | ✗ | ✓ |
| **12** | ✗ | ✓ | ✗ | ✓ |
| **13** | ✗ | ✓ | ✗ | ✓ |
| **14** | ✗ | ✗ | MAFFT, IQ-TREE, MACSE (opt) | ✓ |

**Legend:**
- ✓ = Fully supported
- ✗ = Not required / Not applicable
- ⚠ = Limited support / Use with caution

---

## Next Steps

After completing setup:

1. ✅ **Verify** all checks pass
2. ✅ **Read** [USAGE.md](USAGE.md) for detailed usage instructions
3. ✅ **Prepare** your input data
4. ✅ **Start** with appropriate module for your data type
5. ✅ **Follow** the complete workflow in USAGE.md

---

## Quick Support

- **Detailed Installation**: [INSTALL.md](INSTALL.md)
- **Usage Instructions**: [USAGE.md](USAGE.md)
- **GitHub Issues**: [Report Problems](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS/issues)
- **Module READMEs**: Check individual module directories for specific guidance

---

**Installation Time Estimates:**
- Complete Conda: 15-30 minutes
- Minimal Conda: 5-10 minutes  
- Manual: 5 minutes + time for external tools
- R packages: 5-10 minutes
- Total (everything): 30-45 minutes

**Disk Space Requirements:**
- CGAS + Python packages: ~500MB
- Complete Conda environment: ~5GB
- Minimal Conda environment: ~2GB
- External tools: ~500MB-1GB
- Working space: Varies by project

---

*Last updated: January 2026*

**Ready to start? → See [USAGE.md](USAGE.md) for next steps!**
