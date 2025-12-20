#!/usr/bin/env python3
"""
Setup configuration for Chloroplast Genome Analysis Suite (CGAS)
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

setup(
    name='chloroplast-analyzer',
    version='1.0.0',
    author='Abdullah',
    author_email='your.email@example.com',
    description='Chloroplast Genome Analysis Suite (CGAS) - A comprehensive toolkit for analyzing chloroplast genomes',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/chloroplast-analyzer',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    python_requires='>=3.7',
    install_requires=[
        'biopython>=1.79',
        'pandas>=1.3.0',
        'openpyxl>=3.0.9',
        'numpy>=1.21.0',
        'python-docx>=0.8.11',
    ],
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.12',
            'black>=21.0',
            'flake8>=3.9',
        ],
    },
    entry_points={
        'console_scripts': [
            # Main unified command - CGAS
            'cgas=chloroplast_analyzer.cli:main',
            'cgas-analyze=chloroplast_analyzer.cli:main',
            
            # Individual module commands
            'cgas-count=chloroplast_analyzer.cli:module1',
            'cgas-table=chloroplast_analyzer.cli:module2',
            'cgas-compare=chloroplast_analyzer.cli:module3',
            'cgas-codon=chloroplast_analyzer.cli:module4',
            'cgas-aa=chloroplast_analyzer.cli:module5',
            'cgas-snp=chloroplast_analyzer.cli:module6',
            'cgas-intron=chloroplast_analyzer.cli:module7',
            'cgas-ssr=chloroplast_analyzer.cli:module8',
            'cgas-diversity=chloroplast_analyzer.cli:module9',
        ],
    },
    include_package_data=True,
    package_data={
        'chloroplast_analyzer': ['data/*'],
    },
    keywords='chloroplast genome bioinformatics genomics analysis plastid CGAS',
    project_urls={
        'Bug Reports': 'https://github.com/yourusername/chloroplast-analyzer/issues',
        'Source': 'https://github.com/yourusername/chloroplast-analyzer',
        'Documentation': 'https://github.com/yourusername/chloroplast-analyzer/blob/main/README.md',
    },
)
