# Installation Guide for CGAS v1.0.1

> **Chloroplast Genome Analysis Suite (CGAS)** - Complete installation instructions for all platforms

This guide provides step-by-step installation instructions for CGAS on Linux, macOS, and Windows. Choose the installation method that best fits your needs.

---

## Table of Contents

- [System Requirements](#system-requirements)
- [Quick Start - Recommended Methods](#quick-start---recommended-methods)
- [Installation Options](#installation-options)
  - [Option 1: Complete Conda Environment (Recommended)](#option-1-complete-conda-environment-recommended)
  - [Option 2: Minimal Conda Environment](#option-2-minimal-conda-environment)
  - [Option 3: Manual Installation Without Conda](#option-3-manual-installation-without-conda)
- [External Dependencies Installation](#external-dependencies-installation)
- [Platform-Specific Instructions](#platform-specific-instructions)
- [R Package Installation](#r-package-installation)
- [Optional Tools Installation](#optional-tools-installation)
- [Verification and Testing](#verification-and-testing)
- [Troubleshooting](#troubleshooting)
- [Docker Installation](#docker-installation)
- [Uninstallation](#uninstallation)

---

## System Requirements

### Minimum Requirements
- **Python**: 3.9 or higher
- **RAM**: 4GB minimum (8GB+ recommended for large datasets)
- **Disk Space**: 5GB for CGAS + dependencies, plus additional space for your data
- **Operating System**: 
  - Linux (Ubuntu 18.04+, CentOS 7+, or similar)
  - macOS (10.14 Mojave or later)
  - Windows 10/11

### Recommended Requirements
- **Python**: 3.9 or higher
- **RAM**: 16GB or more (for processing multiple genomes)
- **CPU**: 4+ cores for parallel processing
- **Disk Space**: 20GB+ for analyses and results
- **Internet**: Required for initial installation and downloading reference databases

### Module-Specific Requirements

| Modules | Tools Required | Purpose |
|---------|---------------|---------|
| **1** | fastp, GetOrganelle, BWA, SAMtools | Assembly and QC |
| **2** | BLAST+, PGA | Annotation |
| **3-7** | Python only | Comparative analysis |
| **8, 9** | R with packages | Codon usage and amino acid analysis |
| **10-11** | Python only | SNP and intron analysis |
| **12, 13** | R with packages | SSR and diversity analysis |
| **14** | MAFFT, IQ-TREE, MACSE (optional) | Phylogenetic analysis |

---

## Quick Start - Recommended Methods

### For Servers and High-Performance Systems
**→ Use [Option 1: Complete Conda Environment](#option-1-complete-conda-environment-recommended)**
- Installs all dependencies automatically
- Best for running all 14 modules
- Includes R with visualization packages

### For Standard Laptops and Desktops
**→ Use [Option 2: Minimal Conda Environment](#option-2-minimal-conda-environment)**
- Lightweight installation
- Works for modules 3-14
- Install additional tools as needed

### For Custom Installations
**→ Use [Option 3: Manual Installation](#option-3-manual-installation-without-conda)**
- Full control over dependencies
- Requires manual installation of external tools

---

## Installation Options

### Option 1: Complete Conda Environment (Recommended)

This method installs **all dependencies automatically**, including Python, bioinformatics tools, and R packages. Ideal for servers and complete workflows.

**What gets installed:**
- Python 3.9
- All Python dependencies (Biopython, pandas, openpyxl, etc.)
- Bioinformatics tools: fastp, GetOrganelle, BLAST+, MAFFT, IQ-TREE, HyPhy
- Optional tools: MACSE, FastTree
- R 4.0+ with 12 required packages

#### Step-by-Step Installation

```bash
# 1. Clone the repository
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas
cd cgas

# 2. Create conda environment with all dependencies
conda env create -f environment.yml

# 3. Activate the environment
conda activate cgas

# 4. Install CGAS package
pip install -e .

# 5. Verify installation
cgas --list
```

#### What This Installs

The `environment.yml` includes:
- **Python packages**: biopython, pandas, numpy, openpyxl, python-docx, matplotlib, seaborn
- **Bioinformatics tools**: fastp, blast, mafft, iqtree, hyphy
- **Assembly tools**: bwa, samtools, getorganelle
- **R and packages**: ggplot2, dplyr, tidyr, readr, writexl, seqinr, patchwork, RColorBrewer, ggrepel, scales, cowplot, gridExtra

#### Verification

```bash
# Check CGAS version
cgas --version

# List all modules
cgas --list

# Verify external tools
which mafft
which blastn
which iqtree
which R
```

---

### Option 2: Minimal Conda Environment

Lightweight setup ideal for modules 3-14. This option is perfect for standard laptops and desktops where you don't need assembly (Module 1) or annotation (Module 2).

**What gets installed:**
- Python 3.9
- Core Python dependencies
- MAFFT and IQ-TREE for phylogenetic analysis
- Basic R installation (you'll install R packages separately if needed)

#### Step-by-Step Installation

```bash
# 1. Clone the repository
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas
cd cgas

# 2. Create minimal conda environment
conda env create -f environment-minimal.yml

# 3. Activate the environment
conda activate cgas-minimal

# 4. Install CGAS package
pip install -e .

# 5. Verify installation
cgas --list
```

#### Adding Additional Tools Later

If you need tools for specific modules:

```bash
# For Module 1 (Assembly)
conda install -c bioconda fastp getorganelle bwa samtools

# For Module 2 (Annotation)
conda install -c bioconda blast

# For visualization modules (8, 9, 12, 13)
# See R Package Installation section below
```

---

### Option 3: Manual Installation Without Conda

Full control over your installation. Choose this if you prefer managing dependencies yourself or need specific versions of tools.

#### Step 1: Install Python 3.9+

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install python3.9 python3.9-dev python3.9-venv python3-pip git
```

**macOS:**
```bash
# Install Homebrew if needed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Python
brew install python@3.9 git
```

**Windows:**
1. Download Python 3.9+ from [python.org](https://www.python.org/downloads/windows/)
2. Run installer with "Add Python to PATH" checked
3. Install Git from [git-scm.com](https://git-scm.com/download/win)

#### Step 2: Create Virtual Environment

**Linux/macOS:**
```bash
# Clone repository
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas
cd cgas

# Create virtual environment
python3.9 -m venv cgas_env
source cgas_env/bin/activate

# Install CGAS and Python dependencies
pip install -e .
pip install -r requirements.txt
```

**Windows:**
```cmd
# Clone repository
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas
cd cgas

# Create virtual environment
python -m venv cgas_env
cgas_env\Scripts\activate

# Install CGAS and Python dependencies
pip install -e .
pip install -r requirements.txt
```

#### Step 3: Verify CGAS Installation

```bash
# Check installation
cgas --list
cgas --help

# You should see all 14 modules listed
```

**Note:** External bioinformatics tools must be installed separately (see [External Dependencies Installation](#external-dependencies-installation)).

---

## External Dependencies Installation

Install these tools based on which CGAS modules you plan to use.

### Quick Reference Table

| Tool | Modules | Installation Priority | Required For |
|------|---------|----------------------|--------------|
| **fastp** | 1 | High (if using Module 1) | Read quality control |
| **GetOrganelle** | 1 | High (if using Module 1) | Genome assembly |
| **BWA** | 1 | High (if using Module 1) | Read mapping |
| **SAMtools** | 1 | High (if using Module 1) | BAM processing |
| **BLAST+** | 2 | High (if using Module 2) | Gene annotation |
| **PGA** | 2 | High (if using Module 2) | Plastome annotation |
| **MAFFT** | 13, 14 | Medium | Multiple alignment |
| **IQ-TREE** | 14 | Medium | Phylogenetic trees |
| **R** | 8, 9, 12, 13 | Medium | Visualization |
| **MACSE** | 14 | Low (optional) | Codon-aware alignment |
| **Java 8+** | 14 | Low (for MACSE) | MACSE runtime |
| **FastTree** | 14 | Low (optional) | Alternative phylogeny |

### Platform-Specific Installation

#### Ubuntu/Debian (20.04+)

**Core Tools (all modules except 1-2):**
```bash
sudo apt-get update
sudo apt-get install -y \
    mafft \
    ncbi-blast+ \
    r-base \
    iqtree
```

**Assembly Tools (Module 1):**
```bash
# Install system tools
sudo apt-get install -y \
    fastp \
    samtools \
    bwa

# Install GetOrganelle via pip
pip install getorganelle

# Download GetOrganelle database
get_organelle_config.py --add embplant_pt
```

**Annotation Tools (Module 2):**
```bash
# BLAST+ (if not already installed)
sudo apt-get install -y ncbi-blast+

# PGA - requires manual installation
# See PGA Installation section below
```

#### macOS (Homebrew)

**Core Tools:**
```bash
brew install mafft blast iqtree r
```

**Assembly Tools (Module 1):**
```bash
brew install fastp samtools bwa

# GetOrganelle via pip
pip install getorganelle
get_organelle_config.py --add embplant_pt
```

**Java (for MACSE):**
```bash
brew install openjdk
```

#### Windows (Chocolatey)

**Install Chocolatey first:**
```powershell
Set-ExecutionPolicy Bypass -Scope Process -Force
[System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072
iex ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1'))
```

**Install Tools:**
```powershell
# Core tools
choco install -y mafft r.project

# Note: Some tools may require manual installation on Windows
# See manual installation instructions below
```

---

## R Package Installation

CGAS modules 8, 9, 12, and 13 require R with specific packages for visualization.

### Installing R

**Ubuntu/Debian:**
```bash
sudo apt-get install r-base r-base-dev
```

**macOS:**
```bash
brew install r
```

**Windows:**
Download and install from [CRAN](https://cran.r-project.org/bin/windows/base/)

### Installing Required R Packages

Run this in R or from command line:

```bash
# Install all required packages
R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'readr', 'writexl', 'seqinr', 'patchwork', 'RColorBrewer', 'ggrepel', 'scales', 'cowplot', 'gridExtra'), repos='https://cran.rstudio.com/')"
```

Or install packages individually as needed:

```bash
# For Module 8 (Codon Usage Analysis)
R -e "install.packages(c('ggplot2', 'seqinr', 'dplyr', 'tidyr', 'RColorBrewer', 'patchwork'), repos='https://cran.rstudio.com/')"

# For Module 9 (Amino Acid Analysis)
R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'scales'), repos='https://cran.rstudio.com/')"

# For Module 12 (SSR Analysis)
R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'RColorBrewer', 'patchwork', 'ggrepel', 'scales'), repos='https://cran.rstudio.com/')"

# For Module 13 (Nucleotide Diversity)
R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'cowplot', 'gridExtra'), repos='https://cran.rstudio.com/')"
```

### Verify R Installation

```bash
R --version
R -e "library(ggplot2); library(dplyr); library(tidyr)"
```

---

## Optional Tools Installation

### PGA (Plastome Genome Annotator) - Module 2

PGA requires manual installation from GitHub.

```bash
# Clone PGA repository
git clone https://github.com/quxiaojian/PGA.git ~/tools/PGA
cd ~/tools/PGA

# Install Perl dependencies (Ubuntu/Debian)
sudo apt-get install -y \
    perl \
    libdbi-perl \
    libdbd-mysql-perl \
    libxml-simple-perl \
    bioperl

# Make PGA executable
chmod +x PGA.pl

# Test PGA
perl PGA.pl

# Note the path for use with CGAS
# Example: --pga /home/username/tools/PGA/PGA.pl
```

**macOS:**
```bash
# Install Perl modules via CPAN
cpan install DBI DBD::mysql XML::Simple Bio::Perl
```

### MACSE (Codon-Aware Alignment) - Module 14

MACSE enables codon-aware alignment for phylogenetic analysis. Requires Java 8+.

```bash
# 1. Install Java if needed
# Ubuntu/Debian:
sudo apt-get install default-jre

# macOS:
brew install openjdk

# Windows: Download from java.com

# 2. Download MACSE
mkdir -p ~/tools
cd ~/tools
wget https://bioweb.supagro.inra.fr/macse/releases/macse_v2.07.jar

# 3. Create executable script (Linux/macOS)
cat > ~/tools/macse << 'EOF'
#!/bin/bash
java -jar ~/tools/macse_v2.07.jar "$@"
EOF

chmod +x ~/tools/macse

# 4. Add to PATH
echo 'export PATH="$HOME/tools:$PATH"' >> ~/.bashrc
source ~/.bashrc

# 5. Test MACSE
macse -help
```

### IQ-TREE 2 (Latest Version) - Module 14

If the system package is outdated, install the latest version:

```bash
# Download latest release
cd ~/tools
wget https://github.com/iqtree/iqtree2/releases/download/v2.3.0/iqtree-2.3.0-Linux.tar.gz
tar -xzf iqtree-2.3.0-Linux.tar.gz

# Add to PATH
echo 'export PATH="$HOME/tools/iqtree-2.3.0-Linux/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Verify
iqtree2 --version
```

### FastTree (Alternative Phylogenetic Tool) - Module 14

```bash
# Ubuntu/Debian
sudo apt-get install fasttree

# macOS
brew install fasttree

# From source (Linux/macOS)
wget http://www.microbesonline.org/fasttree/FastTree
chmod +x FastTree
sudo mv FastTree /usr/local/bin/
```

---

## Verification and Testing

### Step 1: Verify CGAS Installation

```bash
# Check CGAS version
cgas --version
# Expected: CGAS v1.0.1

# List all modules
cgas --list
# Expected: Shows modules 1-14 with descriptions

# View help
cgas --help
```

### Step 2: Check Python Dependencies

```bash
# Create verification script
cat > verify_python.py << 'EOF'
import sys
print("Checking Python dependencies...")

dependencies = {
    'Bio': 'Biopython',
    'pandas': 'pandas',
    'numpy': 'numpy',
    'openpyxl': 'openpyxl',
    'docx': 'python-docx',
    'matplotlib': 'matplotlib',
    'seaborn': 'seaborn',
}

missing = []
for module, name in dependencies.items():
    try:
        __import__(module)
        print(f"✓ {name}")
    except ImportError:
        print(f"✗ {name} - MISSING")
        missing.append(name)

if missing:
    print(f"\nMissing packages: {', '.join(missing)}")
    print("Install with: pip install " + " ".join(missing))
    sys.exit(1)
else:
    print("\nAll Python dependencies installed!")
EOF

python verify_python.py
```

### Step 3: Check External Tools

```bash
# Create dependency checker
cat > check_tools.sh << 'EOF'
#!/bin/bash
echo "Checking external bioinformatics tools..."
echo "=========================================="

check_tool() {
    if command -v $1 &> /dev/null; then
        echo "✓ $1 installed"
        return 0
    else
        echo "✗ $1 NOT FOUND"
        return 1
    fi
}

# Core tools
check_tool python
check_tool mafft
check_tool blastn
check_tool iqtree || check_tool iqtree2
check_tool R

# Assembly tools (Module 1)
echo -e "\nModule 1 (Assembly) tools:"
check_tool fastp
check_tool get_organelle_from_reads.py
check_tool bwa
check_tool samtools

# Optional tools
echo -e "\nOptional tools:"
check_tool java
check_tool fasttree

echo "=========================================="
EOF

chmod +x check_tools.sh
./check_tools.sh
```

### Step 4: Run Test Analysis

Create test data and run a simple module:

```bash
# Create test directory
mkdir -p ~/cgas_test
cd ~/cgas_test

# If you have test GenBank files, copy them here
# Or download example chloroplast genomes from NCBI

# Test Module 5 (Gene Comparative Analysis)
# This assumes you have GenBank files in current directory
cgas --module 5

# Check output
ls -R Module5_Gene_Comparative_Analysis/
```

### Step 5: Verify R Packages (if using modules 8, 9, 12, 13)

```bash
# Test R package installation
R -e "
packages <- c('ggplot2', 'dplyr', 'tidyr', 'readr', 'writexl', 
              'seqinr', 'patchwork', 'RColorBrewer', 'ggrepel', 
              'scales', 'cowplot', 'gridExtra')

for (pkg in packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf('✓ %s installed\n', pkg))
  } else {
    cat(sprintf('✗ %s MISSING\n', pkg))
  }
}
"
```

---

## Troubleshooting

### Common Installation Issues

#### Issue 1: Python Version Conflicts

**Problem:** Multiple Python versions causing conflicts

**Solution:**
```bash
# Check all Python versions
python --version
python3 --version
python3.9 --version

# Use specific Python version for venv
python3.9 -m venv cgas_env

# Or specify Python in conda
conda create -n cgas python=3.9
```

#### Issue 2: Permission Denied Errors

**Problem:** Cannot write to system directories

**Solution:**
```bash
# Linux/macOS - Install in user directory
pip install --user -e .

# Or fix ownership
sudo chown -R $USER:$USER ~/.local

# Windows - Run PowerShell as Administrator
# Right-click PowerShell → "Run as Administrator"
```

#### Issue 3: Conda Environment Creation Fails

**Problem:** `environment.yml` fails to create environment

**Solution:**
```bash
# Update conda
conda update conda

# Try with specific channel priority
conda env create -f environment.yml --channel-priority flexible

# Or create minimal environment first
conda create -n cgas python=3.9
conda activate cgas
conda install -c bioconda -c conda-forge mafft blast iqtree
pip install -e .
```

#### Issue 4: "Module Not Found" After Installation

**Problem:** CGAS commands not recognized

**Solution:**
```bash
# Ensure CGAS installed in editable mode
cd /path/to/cgas
pip install -e .

# Check if in PATH
which cgas

# Verify installation
pip show chloroplast-genome-analysis-suite

# If using virtual environment, ensure it's activated
source cgas_env/bin/activate  # Linux/macOS
cgas_env\Scripts\activate      # Windows
```

#### Issue 5: External Tool Not in PATH

**Problem:** Tools installed but not found by CGAS

**Solution:**
```bash
# Find tool location
which mafft
which blastn

# If not found, add to PATH temporarily
export PATH=$PATH:/path/to/tool/bin

# Add permanently (Linux/macOS - add to ~/.bashrc or ~/.zshrc)
echo 'export PATH=$PATH:/path/to/tool/bin' >> ~/.bashrc
source ~/.bashrc

# Windows - Add to System Environment Variables
# System Properties → Environment Variables → Path → Edit → New
```

#### Issue 6: R Package Installation Fails

**Problem:** Cannot install R packages

**Solution:**
```bash
# Ubuntu/Debian - Install R development tools
sudo apt-get install r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev

# Install packages with dependencies
R -e "install.packages('ggplot2', dependencies=TRUE, repos='https://cran.rstudio.com/')"

# If specific package fails, install system dependencies
# Example for XML packages:
sudo apt-get install libxml2-dev

# macOS - May need Xcode Command Line Tools
xcode-select --install
```

#### Issue 7: GetOrganelle Database Not Found

**Problem:** Module 1 fails with database errors

**Solution:**
```bash
# Download and configure GetOrganelle database
get_organelle_config.py --add embplant_pt

# Check database location
get_organelle_config.py --list

# If needed, manually download
get_organelle_config.py --add embplant_pt,embplant_mt
```

#### Issue 8: MACSE Java Errors

**Problem:** MACSE fails to run or Java errors

**Solution:**
```bash
# Check Java version (need Java 8 or higher)
java -version

# Test MACSE manually
java -jar /path/to/macse_v2.07.jar -help

# Increase Java memory if needed
java -Xmx4G -jar /path/to/macse_v2.07.jar -help

# Update shell script with memory allocation
cat > macse << 'EOF'
#!/bin/bash
java -Xmx4G -jar ~/tools/macse_v2.07.jar "$@"
EOF
```

### Getting Help

If you continue to experience issues:

1. **Check existing issues**: [GitHub Issues](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS/issues)

2. **Create a new issue** with the following information:
   - Operating system and version
   - Python version (`python --version`)
   - CGAS version (`cgas --version`)
   - Installation method used (Conda/Manual/Docker)
   - Complete error message and traceback
   - Output of `cgas --list`
   - Output of dependency checker scripts

3. **Module-specific issues**: Check the README file in each module's directory for module-specific troubleshooting

---

## Docker Installation

Docker provides a consistent environment across all platforms. Ideal for reproducibility and avoiding dependency conflicts.

### Prerequisites

Install Docker:
- **Linux**: Follow [official Docker docs](https://docs.docker.com/engine/install/)
- **macOS**: Install [Docker Desktop](https://docs.docker.com/desktop/install/mac-install/)
- **Windows**: Install [Docker Desktop](https://docs.docker.com/desktop/install/windows-install/)

### Building CGAS Docker Image

```bash
# Clone repository
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas
cd cgas

# Build Docker image
docker build -t cgas:v1.0.1 .

# Or use provided Dockerfile
```

**Sample Dockerfile:**

```dockerfile
FROM ubuntu:22.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3.9 \
    python3.9-dev \
    python3.9-venv \
    python3-pip \
    git \
    mafft \
    ncbi-blast+ \
    r-base \
    r-base-dev \
    samtools \
    bwa \
    fasttree \
    wget \
    curl \
    default-jre \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir \
    biopython \
    pandas \
    numpy \
    openpyxl \
    python-docx \
    matplotlib \
    seaborn \
    getorganelle

# Install IQ-TREE
RUN wget https://github.com/iqtree/iqtree2/releases/download/v2.3.0/iqtree-2.3.0-Linux.tar.gz && \
    tar -xzf iqtree-2.3.0-Linux.tar.gz && \
    mv iqtree-2.3.0-Linux/bin/iqtree2 /usr/local/bin/ && \
    rm -rf iqtree-2.3.0-Linux*

# Install R packages
RUN R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'readr', 'writexl', 'seqinr', 'patchwork', 'RColorBrewer', 'ggrepel', 'scales', 'cowplot', 'gridExtra'), repos='https://cran.rstudio.com/')"

# Clone and install CGAS
WORKDIR /opt
RUN git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git cgas && \
    cd cgas && \
    pip3 install -e .

# Setup GetOrganelle database
RUN get_organelle_config.py --add embplant_pt

# Set working directory
WORKDIR /data
ENV PATH="/opt/cgas:$PATH"

# Default command
ENTRYPOINT ["cgas"]
CMD ["--help"]
```

### Using Docker Image

```bash
# Run CGAS help
docker run --rm cgas:v1.0.1 --help

# Run with mounted data directory
docker run --rm -v $(pwd):/data cgas:v1.0.1 --module 5

# Interactive session
docker run --rm -it -v $(pwd):/data cgas:v1.0.1 bash

# Inside container:
cgas --list
cgas --module 5
```

### Docker Compose (Optional)

Create `docker-compose.yml`:

```yaml
version: '3.8'
services:
  cgas:
    build: .
    image: cgas:v1.0.1
    volumes:
      - ./data:/data
      - ./results:/results
    working_dir: /data
    command: ["--help"]
```

Run with:
```bash
docker-compose run cgas --module 5
```

---

## Uninstallation

### Remove CGAS

**If installed via pip:**
```bash
pip uninstall chloroplast-genome-analysis-suite
```

**If installed with conda:**
```bash
# Remove entire conda environment
conda env remove -n cgas

# Or just remove from current environment
conda remove chloroplast-genome-analysis-suite
```

**If installed manually:**
```bash
# Remove source directory
rm -rf ~/cgas

# Remove virtual environment
rm -rf ~/cgas_env
```

### Remove External Dependencies

#### Ubuntu/Debian
```bash
sudo apt-get remove --purge \
    mafft \
    ncbi-blast+ \
    r-base \
    r-base-dev \
    samtools \
    bwa \
    fasttree \
    iqtree

sudo apt-get autoremove
sudo apt-get clean
```

#### macOS
```bash
brew uninstall mafft blast r samtools bwa fasttree iqtree
brew cleanup
```

#### Windows
Use "Add or Remove Programs" in Windows Settings to remove installed applications.

### Remove Configuration and Cache Files

```bash
# Remove configuration files (if any)
rm -f ~/.cgas_config.yaml
rm -rf ~/.cgas_cache

# Remove GetOrganelle database
rm -rf ~/.GetOrganelle

# Remove R packages (optional)
R -e "remove.packages(c('ggplot2', 'dplyr', 'tidyr', 'seqinr', 'patchwork', 'RColorBrewer', 'ggrepel', 'scales', 'cowplot', 'gridExtra'))"
```

### Complete Cleanup

```bash
# Remove all CGAS-related files and dependencies
rm -rf ~/cgas ~/cgas_env ~/.cgas* ~/.GetOrganelle
pip uninstall chloroplast-genome-analysis-suite
conda env remove -n cgas
```

---

## Post-Installation Steps

### 1. Activate Your Environment

**If using conda:**
```bash
conda activate cgas
```

**If using virtual environment:**
```bash
source cgas_env/bin/activate  # Linux/macOS
cgas_env\Scripts\activate      # Windows
```

### 2. Prepare Your Data

Organize your input files:
```bash
mkdir -p ~/my_project/{raw_reads,genomes,reference,results}
```

### 3. Read the Documentation

- Main README: `README.md`
- Module-specific guides: Check each module directory
- Usage examples: `cgas --help` and `cgas --list`

### 4. Start Using CGAS

```bash
# List all modules
cgas --list

# Get help for specific module
cgas --module 1 --help

# Run your first analysis
cgas --module 5 -i genomes/ -o results/
```

---

## Next Steps

After successful installation:

1. ✅ **Verify installation** with `cgas --list`
2. 📖 **Read the README** for module descriptions
3. 🧪 **Run a test analysis** with sample data
4. 📊 **Process your data** starting with Module 1 or Module 3
5. 🎯 **Explore all 14 modules** for complete chloroplast genome analysis

For detailed usage instructions, see the main [README.md](README.md) file and module-specific documentation.

---

## Support and Contact

- **Issues**: [GitHub Issues](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS/issues)
- **Repository**: [GitHub](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS)
- **Documentation**: Check module README files for specific guidance

---

**Made with ❤️ for chloroplast genomics research**

*Last updated: January 2026*
