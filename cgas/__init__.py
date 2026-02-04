"""
CGAS - Chloroplast Genome Analysis Suite
==========================================

A comprehensive toolkit for analyzing chloroplast genomes with 14 specialized modules.

Modules:
    1. Chloroplast Genome Assembly
    2. Plastome Annotation
    3. Plastome Gene Comparison
    4. GenBank Format Conversion
    5. Gene Comparative Analysis
    6. Gene Content Tables
    7. Comparative Genome Analysis
    8. Codon Usage Analysis (RSCU)
    9. Amino Acid Analysis
    10. SNP/Substitution Analysis
    11. Gene and tRNA Intron Analysis
    12. Comprehensive SSR Analysis
    13. Nucleotide Diversity Analysis
    14. Phylogenetic Matrix Builder

Author: Abdullah
Version: 1.0.1
"""

__version__ = '1.0.1'
__author__ = 'Abdullah'
__email__ = 'your.email@example.com'  # Replace with your actual email
__description__ = 'A comprehensive toolkit for analyzing chloroplast genomes'
__url__ = 'https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS'

# Import the unified analyzer
try:
    from . import unified_analyzer
except ImportError:
    unified_analyzer = None

# Make CLI accessible
try:
    from . import cli
except ImportError:
    cli = None

# Define what gets imported with `from cgas import *`
__all__ = [
    '__version__',
    '__author__',
    '__email__',
    '__description__',
    '__url__',
]

# Only include modules in __all__ if they imported successfully
if unified_analyzer is not None:
    __all__.append('unified_analyzer')

if cli is not None:
    __all__.append('cli')