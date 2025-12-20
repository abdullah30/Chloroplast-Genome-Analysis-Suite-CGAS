# Contributing to Chloroplast Analyzer

Thank you for your interest in contributing to the Chloroplast Genome Analyzer! This document provides guidelines and instructions for contributing.

## 📋 Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)

## Code of Conduct

### Our Pledge

We are committed to providing a welcoming and inspiring community for all. Please be respectful and constructive in all interactions.

### Expected Behavior

- Be respectful and inclusive
- Provide constructive feedback
- Focus on what is best for the community
- Show empathy towards other community members

## Getting Started

### Prerequisites

- Python 3.7 or higher
- Git
- Familiarity with BioPython
- Basic understanding of genomics concepts

### Setting Up Development Environment

1. **Fork the repository**
   ```bash
   # Click "Fork" on GitHub
   # Clone your fork
   git clone https://github.com/YOUR_USERNAME/chloroplast-analyzer.git
   cd chloroplast-analyzer
   ```

2. **Create virtual environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # Linux/macOS
   venv\Scripts\activate     # Windows
   ```

3. **Install in development mode**
   ```bash
   pip install -e ".[dev]"
   ```

4. **Install pre-commit hooks** (optional but recommended)
   ```bash
   pip install pre-commit
   pre-commit install
   ```

## How to Contribute

### Types of Contributions

We welcome various types of contributions:

#### 🐛 Bug Reports
- Use GitHub Issues
- Include clear title and description
- Provide steps to reproduce
- Include error messages and screenshots
- Specify your environment (OS, Python version)

#### 💡 Feature Requests
- Use GitHub Issues with "enhancement" label
- Clearly describe the feature
- Explain the use case
- Consider implementation details

#### 📝 Documentation
- Fix typos or clarify existing docs
- Add examples and tutorials
- Improve API documentation
- Translate documentation

#### 🔧 Code Contributions
- Bug fixes
- New analysis modules
- Performance improvements
- Code refactoring

#### 🧪 Test Contributions
- Add unit tests
- Add integration tests
- Improve test coverage

## Development Setup

### Branch Naming Convention

```
feature/add-new-module        # New features
bugfix/fix-codon-calculation  # Bug fixes
docs/update-readme            # Documentation
refactor/improve-cli          # Code refactoring
test/add-module4-tests        # Tests
```

### Creating a Branch

```bash
# Update main branch
git checkout main
git pull origin main

# Create feature branch
git checkout -b feature/your-feature-name
```

## Coding Standards

### Python Style Guide

We follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) with some modifications:

- **Line Length**: 100 characters (not 79)
- **Indentation**: 4 spaces
- **Quotes**: Double quotes for strings
- **Imports**: Grouped and sorted

### Code Formatting

Use **Black** for code formatting:

```bash
black chloroplast_analyzer/
```

### Linting

Use **flake8** for linting:

```bash
flake8 chloroplast_analyzer/
```

### Type Hints

Use type hints where appropriate:

```python
def analyze_genome(file_path: str, threshold: int = 10) -> dict:
    """
    Analyze genome file.
    
    Args:
        file_path: Path to GenBank file
        threshold: Minimum threshold value
        
    Returns:
        Dictionary with analysis results
    """
    pass
```

### Docstrings

Use Google-style docstrings:

```python
def process_sequence(sequence, quality_threshold=0.8):
    """Process DNA sequence with quality filtering.
    
    This function processes a DNA sequence and filters based on
    quality scores.
    
    Args:
        sequence (str): DNA sequence string
        quality_threshold (float): Minimum quality score (0-1)
        
    Returns:
        str: Filtered sequence
        
    Raises:
        ValueError: If quality_threshold is not between 0 and 1
        
    Example:
        >>> process_sequence("ATCG", 0.9)
        "ATCG"
    """
    pass
```

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_module1.py

# Run with coverage
pytest --cov=chloroplast_analyzer tests/

# Run with verbose output
pytest -v
```

### Writing Tests

Place tests in `tests/` directory:

```python
# tests/test_module1.py
import pytest
from chloroplast_analyzer import module1_gene_count

def test_gene_count():
    """Test gene counting functionality."""
    result = module1_gene_count.count_genes("test_data/sample.gb")
    assert result["total_genes"] > 0
    assert "cds_count" in result

def test_invalid_file():
    """Test handling of invalid file."""
    with pytest.raises(FileNotFoundError):
        module1_gene_count.count_genes("nonexistent.gb")
```

### Test Coverage

Aim for >80% test coverage:

```bash
pytest --cov=chloroplast_analyzer --cov-report=html
# Open htmlcov/index.html to view coverage report
```

## Documentation

### Code Documentation

- Document all public functions and classes
- Include examples in docstrings
- Update README.md for new features
- Add entries to CHANGELOG.md

### Documentation Files

- **README.md** - Overview and quick start
- **INSTALL.md** - Installation instructions
- **USAGE.md** - Detailed usage guide
- **CHANGELOG.md** - Version history
- **docs/** - Additional documentation

### Building Documentation

If we add Sphinx documentation:

```bash
cd docs
make html
```

## Submitting Changes

### Commit Messages

Write clear, descriptive commit messages:

```
feat: Add nucleotide diversity sliding window analysis

- Implement sliding window algorithm
- Add window size parameter
- Update documentation
- Add unit tests

Closes #123
```

**Commit Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `test`: Adding tests
- `refactor`: Code refactoring
- `style`: Formatting changes
- `perf`: Performance improvements
- `chore`: Maintenance tasks

### Pull Request Process

1. **Create Pull Request**
   ```bash
   # Push your branch
   git push origin feature/your-feature-name
   
   # Create PR on GitHub
   ```

2. **PR Description Template**
   ```markdown
   ## Description
   Brief description of changes
   
   ## Type of Change
   - [ ] Bug fix
   - [ ] New feature
   - [ ] Documentation update
   - [ ] Performance improvement
   
   ## Testing
   - [ ] Tests added/updated
   - [ ] All tests passing
   - [ ] Documentation updated
   
   ## Checklist
   - [ ] Code follows style guidelines
   - [ ] Self-review completed
   - [ ] Comments added for complex code
   - [ ] Documentation updated
   - [ ] No new warnings generated
   
   ## Related Issues
   Closes #123
   ```

3. **Review Process**
   - Maintainers will review your PR
   - Address any requested changes
   - Keep PR focused and reasonably sized
   - Be responsive to feedback

4. **After Approval**
   - PR will be merged by maintainers
   - Branch will be deleted
   - Your contribution will be in next release!

## Adding New Modules

### Module Template

```python
#!/usr/bin/env python3
"""
Module X: [Module Name]

Description of what this module does.

Author: Your Name
Date: YYYY-MM-DD
"""

from pathlib import Path
from Bio import SeqIO
import pandas as pd

def main():
    """Main entry point for Module X."""
    print("="*80)
    print("MODULE X: [Module Name]")
    print("="*80)
    
    # Your code here
    
    print("\n✓ Module X completed")

if __name__ == "__main__":
    main()
```

### Module Guidelines

1. **Input**: Use standard GenBank or FASTA files
2. **Output**: Generate timestamped files
3. **Error Handling**: Catch and report errors gracefully
4. **Documentation**: Include docstrings and comments
5. **Testing**: Add unit tests
6. **CLI Integration**: Add entry point in `cli.py`

### Adding Module to System

1. Create `moduleX_name.py` in `chloroplast_analyzer/`
2. Add entry to `unified_analyzer.py` MODULES dict
3. Add CLI command in `cli.py`
4. Update `setup.py` entry_points
5. Add tests in `tests/test_moduleX.py`
6. Update documentation

## Development Workflow

### Typical Workflow

```bash
# 1. Create feature branch
git checkout -b feature/new-analysis

# 2. Make changes
# Edit files...

# 3. Test changes
pytest

# 4. Format code
black chloroplast_analyzer/
flake8 chloroplast_analyzer/

# 5. Commit changes
git add .
git commit -m "feat: Add new analysis feature"

# 6. Push to GitHub
git push origin feature/new-analysis

# 7. Create Pull Request on GitHub
```

### Keeping Your Fork Updated

```bash
# Add upstream remote (once)
git remote add upstream https://github.com/yourusername/chloroplast-analyzer.git

# Update your fork
git checkout main
git fetch upstream
git merge upstream/main
git push origin main
```

## Release Process

### Version Numbering

We use [Semantic Versioning](https://semver.org/):

- **MAJOR.MINOR.PATCH** (e.g., 2.1.0)
- Increment MAJOR for breaking changes
- Increment MINOR for new features
- Increment PATCH for bug fixes

### Creating a Release

1. Update version in:
   - `setup.py`
   - `pyproject.toml`
   - `chloroplast_analyzer/__init__.py`
   - `CHANGELOG.md`

2. Create release commit:
   ```bash
   git add .
   git commit -m "chore: Bump version to 2.2.0"
   ```

3. Create tag:
   ```bash
   git tag -a v2.2.0 -m "Release version 2.2.0"
   git push origin v2.2.0
   ```

4. Create GitHub Release with changelog

## Getting Help

### Communication Channels

- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: General questions and discussions
- **Email**: your.email@example.com

### Resources

- [BioPython Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
- [Python Packaging Guide](https://packaging.python.org/)
- [Git Documentation](https://git-scm.com/doc)

## Recognition

Contributors will be:
- Listed in CONTRIBUTORS.md
- Mentioned in release notes
- Given credit in academic citations (if applicable)

## Questions?

Don't hesitate to ask! We're here to help:
- Open a GitHub Issue
- Start a Discussion
- Email the maintainers

Thank you for contributing! 🧬✨
