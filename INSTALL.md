# Installation Guide

## Quick Install (Recommended)

### From GitHub (Direct Install)

```bash
# Install directly from GitHub
pip install git+https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git

# Or with specific version/branch
pip install git+https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git@v1.0.0
pip install git+https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git@main
```

### From PyPI (when published)

```bash
pip install Chloroplast-Genome-Analysis-Suite-CGAS
```

## Manual Installation

### 1. Clone the Repository

```bash
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git
cd Chloroplast-Genome-Analysis-Suite-CGAS
```

### 2. Install in Development Mode

```bash
# Install in editable mode (recommended for development)
pip install -e .

# Or install normally
pip install .
```

### 3. Install with Optional Dependencies

```bash
# Install with development tools
pip install -e ".[dev]"
```

## System Requirements

- **Python**: 3.7 or higher
- **Operating System**: Linux, macOS, or Windows
- **MAFFT**: Required for Module 9 (Nucleotide Diversity)

## Installing MAFFT

### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install mafft
```

### macOS
```bash
brew install mafft
```

### CentOS/RHEL
```bash
sudo yum install mafft
```

### Windows
Download from: https://mafft.cbrc.jp/alignment/software/

### Verify Installation
```bash
mafft --version
```

## Virtual Environment (Recommended)

### Using venv

```bash
# Create virtual environment
python3 -m venv chloroplast_env

# Activate (Linux/macOS)
source chloroplast_env/bin/activate

# Activate (Windows)
chloroplast_env\Scripts\activate

# Install package
pip install git+https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git

# When done
deactivate
```

### Using Conda

```bash
# Create conda environment
conda create -n chloroplast python=3.9

# Activate
conda activate chloroplast

# Install dependencies
conda install -c conda-forge biopython pandas openpyxl numpy python-docx

# Install package
pip install git+https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git

# When done
conda deactivate
```

## Verification

After installation, verify the package is working:

```bash
# Check version
cgas --version

# List available modules
cgas --list

# Show help
cgas --help
```

## Troubleshooting

### "command not found: cgas"

The installation path is not in your PATH. Try:

```bash
# Find where pip installed scripts
pip show Chloroplast-Genome-Analysis-Suite-CGAS

# Add to PATH (add to ~/.bashrc or ~/.zshrc)
export PATH="$HOME/.local/bin:$PATH"

# Or use python -m
python -m chloroplast_analyzer.cli --help
```

### "No module named 'Bio'"

BioPython not installed:

```bash
pip install biopython
```

### Permission Denied

Use `--user` flag:

```bash
pip install --user git+https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git
```

### Upgrade Existing Installation

```bash
pip install --upgrade git+https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git
```

### Uninstall

```bash
pip uninstall Chloroplast-Genome-Analysis-Suite-CGAS
```

## For Developers

### Clone for Development

```bash
git clone https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS.git
cd Chloroplast-Genome-Analysis-Suite-CGAS
pip install -e ".[dev]"
```

### Run Tests

```bash
pytest tests/
```

### Code Formatting

```bash
black chloroplast_analyzer/
flake8 chloroplast_analyzer/
```

## Next Steps

After installation, see [README.md](README.md) for usage instructions.
