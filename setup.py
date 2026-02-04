#!/usr/bin/env python3
from setuptools import setup, find_packages
from pathlib import Path

# ---------------------------------------------------------------------
# Basic package metadata
# ---------------------------------------------------------------------
PACKAGE_NAME = "cgas"
VERSION = "1.0.1"
DESCRIPTION = "CGAS - Chloroplast Genome Analysis Suite (Modules 1-14)"

# Long description from README
this_dir = Path(__file__).parent
long_description = (this_dir / "README.md").read_text(encoding="utf-8")

# ---------------------------------------------------------------------
# Setup configuration
# ---------------------------------------------------------------------
setup(
    name=PACKAGE_NAME,
    version=VERSION,
    author="Your Name",
    author_email="your.email@example.com",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/CGAS",
    license="MIT",

    # -----------------------------------------------------------------
    # Package discovery
    # -----------------------------------------------------------------
    packages=find_packages(),
    include_package_data=True,

    # -----------------------------------------------------------------
    # Python requirements
    # -----------------------------------------------------------------
    python_requires=">=3.8",

    # -----------------------------------------------------------------
    # Runtime Python dependencies
    # -----------------------------------------------------------------
    install_requires=[
        "biopython>=1.80",
        "pandas>=1.5",
        "numpy>=1.21",
        "openpyxl",
        "python-docx",
    ],

    # -----------------------------------------------------------------
    # Console entry points
    # -----------------------------------------------------------------
    entry_points={
        "console_scripts": [
            # Main unified CLI
            "cgas=cgas.cli:main",

            # Individual module commands
            "cgas-assembly=cgas.cli:assembly",
            "cgas-annotate=cgas.cli:annotate",
            "cgas-compare=cgas.cli:compare",
            "cgas-convert=cgas.cli:convert",
            "cgas-gene-compare=cgas.cli:gene_compare",
            "cgas-gene-table=cgas.cli:gene_table",
            "cgas-genome-compare=cgas.cli:genome_compare",
            "cgas-codon=cgas.cli:codon",
            "cgas-amino=cgas.cli:amino",
            "cgas-snp=cgas.cli:snp",
            "cgas-intron=cgas.cli:intron",
            "cgas-ssr=cgas.cli:ssr",
            "cgas-diversity=cgas.cli:diversity",
            "cgas-phylogeny=cgas.cli:phylogeny",
        ]
    },

    # -----------------------------------------------------------------
    # Classifiers
    # -----------------------------------------------------------------
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],

    # -----------------------------------------------------------------
    # Keywords
    # -----------------------------------------------------------------
    keywords=[
        "chloroplast",
        "plastome",
        "phylogenomics",
        "comparative genomics",
        "bioinformatics",
        "plant genomics",
    ],
)