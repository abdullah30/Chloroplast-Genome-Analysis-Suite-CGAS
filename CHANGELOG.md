# Changelog

All notable changes to CGAS will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.1] - 2026-01-30

### Added
- **NEW: Modules 1-4 - Complete Preparation Pipeline**
  - **Module 1: Chloroplast Genome Assembly**
    - Raw read quality control with Fastp
    - De novo plastome assembly using GetOrganelle
    - Automated quality assessment and filtering
    - Support for paired-end and single end FASTQ files
    - Batch processed of multiple genome assemblies at once
  - **Module 2: Plastome Annotation**
    - Genome annotation using PGA (Plastid Genome Annotator)
    - BLAST-based gene identification and validation
    - Comprehensive feature annotation (genes, tRNAs, rRNAs)
  - **Module 3: Plastome Gene Comparison**
    - Interactive genome comparison and visualization
    - Gene normalization across samples
    - Annotation error detection to avoid inconsistency in the genomes downloaded during study from NCBI
    - Manual curation support
  - **Module 4: GenBank Format Conversion**
    - Convert GenBank to FASTA format
    - Generate NCBI-compliant TBL and fasta files
    - Batch processing of multiple
    - Format validation for submission readiness to avoid annotation errors

- **NEW: Module 14 - Phylogenetic Matrix Builder**
  - Automated phylogenetic tree construction
  - MAFFT-based sequence alignment
  - MACSE integration for codon-aware CDS alignment
  - IQ-TREE support for maximum likelihood trees
  - Outgroup handling and tree rooting
  - Multiple alignment matrix outputs (genes, CDS, intergenic)

- **Enhanced R Visualizations**
  - **Module 8 (Codon Usage)**: Publication-ready RSCU heatmaps and barplots
  - **Module 9 (Amino Acid)**: Amino acid composition heatmaps and comparative plots
  - **Module 10 (SNP Analysis)**: SNP distribution plots and substitution matrices
  - **Module 12 (SSR Analysis)**: SSRs distribution in genomic regions, functional regions, and various type of motifs
  - **Module 13 (Nucleotide Diversity)**: Pi and Theta diversity plots with sliding windows

- **Unified CLI interface** with subcommands for all 14 modules
- **Complete Conda support** with environment.yml
- **MACSE integration** for codon-aware CDS alignment
- **Comprehensive R visualization** with 12 packages
  - Core: ggplot2, reshape2, pheatmap, RColorBrewer
  - Data manipulation: dplyr, tidyr, zoo, readr
  - Enhancement: cowplot, scales

- **Outgroup handling mechanism** for phylogenetic analysis
  - Auto-excludes outgroups from comparative analysis (Modules 5-13)
  - Includes outgroups in phylogenetic trees (Module 14)
  - Supports up to 4 outgroup files

- **Module-specific options** exposed via CLI
  - Module 1: Full assembly parameter control (--read1, --read2, --threads, --kmers)
  - Module 2: --organism-file option for custom reference
  - Module 14: --macse, --genes-only, --complete-only, --iqtree options

- **Interactive installation script** (install.sh)
- **Conda setup script** (setup_conda.sh)
- **Comprehensive documentation**
  - README.md - Main documentation
  - INSTALL.md - Installation guide
  - QUICKREF.md - Command reference
  - CONDA_INSTALL.md - Detailed conda guide
  - CONDA_QUICKSTART.md - Quick conda setup

- **Python 3.9+ optimization** for low-RAM systems
- **requirements.txt** for pip installation

### Changed
- **Major workflow restructure**: Two-phase approach
  - **Phase 1 - Preparation (Modules 1-4)**: Run individually with manual curation
    1. Assembly → 2. Annotation → 3. Comparison → 4. Conversion
  - **Phase 2 - Analysis (Modules 5-14)**: Run together with `cgas all`
    
- **`cgas all` command** now runs modules 5-14 only (main analysis pipeline)
- **Module numbering** expanded from 9 to 14 modules
  - Modules 1-4: Preparation workflow (NEW)
  - Modules 5-13: Core analysis (renumbered)
  - Module 14: Phylogenetics (NEW)

- **Old numbering → New numbering:**
  - Gene Count (1) → (5)
  - Gene Table (2) → (6)
  - Genome Compare (3) → (7)
  - Codon Usage (4) → (8)
  - Amino Acid (5) → (9)
  - SNP Analysis (6) → (10)
  - Intron Analysis (7) → (11)
  - SSR Analysis (8) → (12)
  - Nucleotide Diversity (9) → (13)
  - Phylogeny (NEW) → (14)

- **Version bumped** to 1.1.0 for major feature release

### Enhanced
- **Module 8 (Codon Usage)**
  - Added comprehensive RSCU heatmap generation
  - Interactive barplots for codon preference
  - Statistical summaries with R integration
  
- **Module 9 (Amino Acid Composition)**
  - Publication-ready amino acid heatmaps
  - Comparative composition analysis
  - Enhanced statistical outputs

- **Module 10 (SNP/Substitution Analysis)**
  - Visual SNP distribution plots
  - Substitution type matrices
  - Transition/transversion ratio calculations

  - **Module 12 (Simple Sequence Repeats)**
  - Genenrate figure for SSRs motifs
  - Generate figure for SSRs distibution in genomic regions
  - Generate figure for SSRs distibution in coding and non-coding regions

- **Module 13 (Nucleotide Diversity)**
  - Sliding window diversity plots
  - Pi and Theta comparative visualizations
  - Region-specific diversity analysis

### Fixed
- Module import handling for better compatibility
- Argument passing to individual modules
- File path handling across platforms
- R visualization rendering on different systems
- Memory optimization for large genome datasets

### Documentation
- Complete rewrite of all documentation files
- Detail separate readme file is provided for each module
- Added conda installation guides
- Enhanced usage examples with all 14 modules
- Platform-specific instructions (Linux, macOS, Windows/WSL2)
- Module-by-module workflow guide

---

## [1.0.1] - 2025-01-30

### Initial Public Release
- 9 core analysis modules (Modules 1-9 in old numbering)
- Basic CLI interface
- Initial conda support
- Basic R visualizations files

---

## Release Notes

### v1.0.1 Highlights

This is a major feature release of CGAS with:
- ✅ 14 fully integrated modules (4 preparation + 10 analysis)
- ✅ Complete end-to-end workflow from raw reads to phylogeny
- ✅ Enhanced R visualizations for all analytical modules
- ✅ Phylogenetic tree construction with multiple regions, including complete genome, introns, genes, and protein-coding genes
- ✅ Assembly and annotation pipeline integration
- ✅ Professional publication-ready outputs
- ✅ Comprehensive conda environment
- ✅ Cross-platform support

### Upgrade Guide

**From v1.0.0 to v1.0.1:**

1. **Module renumbering**: If you have scripts using old module numbers (1-9), update them:
   - Old Modules 1-9 → New Modules 5-13
   - Add new Modules 1-4 (preparation)
   - Add new Module 14 (phylogeny)

2. **Update conda environment**:
   ```bash
   conda env update -f environment.yml
   ```

3. **New dependencies**: Modules 1-2 require:
   - Fastp (for Module 1)
   - GetOrganelle (for Module 1)
   - BWA (for module 1)
   - Samtool (for module 1)
   - PGA (for Module 2)

4. **Workflow changes**:
   - Use `cgas assembly`, `cgas annotate`, `cgas compare`, `cgas convert` for preparation
   - Continue using `cgas all` for main analysis (now runs Modules 5-14)

### Breaking Changes
- Module numbers changed (1-9 → 5-13)
- `cgas all` now runs Modules 5-14 instead of all modules
- New preparation workflow required before analysis

### Known Limitations
- PGA (Plastid Genome Annotator) must be installed separately for Module 2 and required pearl. Alternatively, you can use online platform such as Geseq etc.
- GetOrganelle requires separate installation for Module 1
- Some platforms may require manual R package installation. However, require R files for figures are generated in which only change to file path can provide you complete visulization results using RStudio. 
- For best compatibility, Windows users are strongly recommended to use WSL2. In WSL2, installing a fully functional, up-to-date version of Getorganelle is difficult and often fails due to dependency conflicts.  
- Module 14 requires at least 4 GenBank files

### Future Updates
Future updates will be guided by community needs and may include additional alignment algorithms, enhanced visualization features, and support for expanded output formats.

---

