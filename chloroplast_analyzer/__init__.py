"""
Chloroplast Genome Unified Analysis Suite
==========================================

A comprehensive toolkit for analyzing chloroplast genomes.

Modules:
    1. Gene Count and Summary
    2. Gene Table Generation
    3. Comparative Genome Analysis
    4. Codon Usage Analysis (RSCU)
    5. Amino Acid Composition
    6. SNP/Substitution Analysis
    7. Intron Extraction
    8. SSR Analysis
    9. Nucleotide Diversity

Author: Abdullah
Version: 1.0.0
"""

__version__ = '1.0.0'
__author__ = 'Abdullah'

# Import the unified analyzer
try:
    from . import unified_analyzer
except ImportError:
    pass

# Make CLI accessible
try:
    from . import cli
except ImportError:
    pass

__all__ = [
    'unified_analyzer',
    'cli',
    '__version__',
    '__author__',
]
