#!/bin/bash
# setup_github_repo.sh
# Quick setup script to organize chloroplast-analyzer for GitHub

set -e  # Exit on error

echo "=========================================="
echo "Chloroplast Analyzer - GitHub Setup"
echo "=========================================="
echo ""

# Configuration
REPO_NAME="chloroplast-analyzer"
PACKAGE_NAME="chloroplast_analyzer"

# Check if we're in the right place
if [ ! -f "chloroplast_unified_analyzer.py" ]; then
    echo "Error: chloroplast_unified_analyzer.py not found in current directory"
    echo "Please run this script from the directory containing your module files"
    exit 1
fi

echo "Creating repository structure..."

# Create main repo directory
mkdir -p "$REPO_NAME"
cd "$REPO_NAME"

# Create package directory
mkdir -p "$PACKAGE_NAME"
mkdir -p tests
mkdir -p examples/sample_data
mkdir -p docs/modules
mkdir -p scripts

echo "✓ Directory structure created"

# Copy module files
echo ""
echo "Copying module files..."

# Check for module files in parent directory
cd ..
MODULE_FILES=(
    "module1_gene_count.py"
    "module2_gene_table.py"
    "module3_comparative_analysis.py"
    "module4_codon_usage.py"
    "module5_amino_acid.py"
    "module6_snp_analysis.py"
    "module7_intron_extraction.py"
    "module8_ssr_analysis.py"
    "module9_diversity_analysis.py"
)

cd "$REPO_NAME"

for module in "${MODULE_FILES[@]}"; do
    if [ -f "../$module" ]; then
        cp "../$module" "$PACKAGE_NAME/"
        echo "✓ Copied $module"
    else
        echo "⚠ Warning: $module not found (will need to be added manually)"
    fi
done

# Copy and rename unified analyzer
if [ -f "../chloroplast_unified_analyzer.py" ]; then
    cp "../chloroplast_unified_analyzer.py" "$PACKAGE_NAME/unified_analyzer.py"
    echo "✓ Copied unified_analyzer.py"
fi

# Create setup files (these should be created by your previous commands)
echo ""
echo "Creating setup files..."

# If setup files exist in parent, copy them
SETUP_FILES=(
    "setup.py"
    "pyproject.toml"
    "requirements.txt"
    "LICENSE"
    "MANIFEST.in"
    ".gitignore"
    "README.md"
    "INSTALL.md"
    "USAGE.md"
)

for file in "${SETUP_FILES[@]}"; do
    if [ -f "../$file" ]; then
        cp "../$file" .
        echo "✓ Copied $file"
    else
        echo "⚠ Warning: $file not found"
    fi
done

# Copy package __init__ and cli
if [ -f "../__init__.py" ]; then
    cp "../__init__.py" "$PACKAGE_NAME/"
    echo "✓ Copied __init__.py"
fi

if [ -f "../cli.py" ]; then
    cp "../cli.py" "$PACKAGE_NAME/"
    echo "✓ Copied cli.py"
fi

# Initialize git repository
echo ""
echo "Initializing git repository..."
git init
echo "✓ Git repository initialized"

# Create initial commit
echo ""
echo "Creating initial commit..."
git add .
git commit -m "Initial commit: Chloroplast Analyzer v2.1.0

- Added all 9 analysis modules
- Setup package structure for pip installation
- Added CLI commands for each module
- Created comprehensive documentation"

echo "✓ Initial commit created"

# Print summary
echo ""
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "Repository created at: $(pwd)"
echo ""
echo "Next steps:"
echo ""
echo "1. Review and update files:"
echo "   - Edit README.md with your GitHub username"
echo "   - Update setup.py email address"
echo "   - Check that all module files are present"
echo ""
echo "2. Create GitHub repository:"
echo "   - Go to https://github.com/new"
echo "   - Name: $REPO_NAME"
echo "   - Don't initialize with README (we already have one)"
echo ""
echo "3. Push to GitHub:"
echo "   git remote add origin https://github.com/YOUR_USERNAME/$REPO_NAME.git"
echo "   git branch -M main"
echo "   git push -u origin main"
echo ""
echo "4. Test installation:"
echo "   pip install git+https://github.com/YOUR_USERNAME/$REPO_NAME.git"
echo "   chloroplast-analyze --version"
echo ""
echo "5. Optional - Create release tag:"
echo "   git tag -a v2.1.0 -m 'Release version 2.1.0'"
echo "   git push origin v2.1.0"
echo ""
echo "=========================================="
echo ""
echo "📚 Documentation files:"
ls -1 *.md
echo ""
echo "📦 Package contents:"
ls -1 "$PACKAGE_NAME/"
echo ""
echo "Happy analyzing! 🧬"
