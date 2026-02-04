# CGAS Usage Guide

> **Chloroplast Genome Analysis Suite (CGAS)** - Complete usage instructions for all 14 modules

This guide provides detailed usage instructions for CGAS modules, from raw read processing to phylogenetic analysis.

---

## Table of Contents

- [Getting Started](#getting-started)
- [General Usage](#general-usage)
- [Complete Workflow](#complete-workflow)
- [Phase 1: Preparation Modules (1-4)](#phase-1-preparation-modules-1-4)
  - [Module 1: Genome Assembly & QC](#module-1-genome-assembly--qc)
  - [Module 2: Plastome Annotation](#module-2-plastome-annotation)
  - [Module 3: Gene Comparison & Normalization](#module-3-gene-comparison--normalization)
  - [Module 4: Format Conversion & NCBI Submission](#module-4-format-conversion--ncbi-submission)
- [Phase 2: Main Analysis Modules (5-14)](#phase-2-main-analysis-modules-5-14)
  - [Module 5: Gene Comparative Analysis](#module-5-gene-comparative-analysis)
  - [Module 6: Gene Content Tables](#module-6-gene-content-tables)
  - [Module 7: Genome Structure Analysis](#module-7-genome-structure-analysis)
  - [Module 8: Codon Usage (RSCU)](#module-8-codon-usage-rscu)
  - [Module 9: Amino Acid Analysis](#module-9-amino-acid-analysis)
  - [Module 10: SNP Analysis](#module-10-snp-analysis)
  - [Module 11: Intron Analysis](#module-11-intron-analysis)
  - [Module 12: SSR Analysis](#module-12-ssr-analysis)
  - [Module 13: Nucleotide Diversity](#module-13-nucleotide-diversity)
  - [Module 14: Phylogenetic Analysis](#module-14-phylogenetic-analysis)
- [Jupyter Notebook Usage](#jupyter-notebook-usage)
- [Direct Module Commands](#direct-module-commands)
- [Best Practices](#best-practices)
- [Tips and Tricks](#tips-and-tricks)
- [Common Workflows](#common-workflows)
- [Troubleshooting](#troubleshooting)

---

## Getting Started

### Prerequisites

Before using CGAS, ensure you have:

1. **CGAS installed** (see [INSTALL.md](INSTALL.md))
2. **Environment activated**:
   ```bash
   conda activate cgas  # If using conda
   # OR
   source cgas_env/bin/activate  # If using virtual environment
   ```
3. **Input data prepared** in appropriate directories

### Quick Help

```bash
# Show CGAS help
cgas --help

# List all available modules
cgas --list

# Get help for specific module
cgas --module 1 --help
cgas --module 5 --help
```

---

## General Usage

### Basic Command Structure

```bash
cgas --module <MODULE_NUMBER> [OPTIONS]
```

### Common Options

```bash
-i, --input       Input directory or file
-o, --output      Output directory (default: ModuleX_Name/)
-t, --threads     Number of CPU threads (default: 4)
-h, --help        Show help message
```

### Running Multiple Modules

```bash
# Run modules sequentially
cgas --modules 5,6,7

# Run all comparative analysis modules
cgas --modules 5,6,7,8,9,10,11,12,13
```

---

## Complete Workflow

CGAS follows a two-phase workflow:

### **Phase 1: Preparation** (Modules 1-4)
Run these modules **individually** in sequence to prepare your data.

```
Raw Reads → Assembly (1) → Annotation (2) → Normalization (3) → Format Conversion (4)
```

### **Phase 2: Main Analysis** (Modules 5-14)
Run these modules on **finalized GenBank files** for comparative analysis.

```
GenBank Files → Comparative Analysis (5-13) → Phylogeny (14)
```

---

## Phase 1: Preparation Modules (1-4)

These modules prepare your data for analysis. **Run them individually** in the correct order.

---

### Module 1: Genome Assembly & QC

**Purpose**: Process raw sequencing reads and assemble chloroplast genomes with comprehensive quality control.

**Input**: FASTQ files (paired-end or single-end raw reads)  
**Output**: Assembled genomes, QC reports, coverage analysis

#### Basic Usage

```bash
# Simple assembly
cgas --module 1 -i raw_reads/ -o results/
```

#### Input Directory Structure

Your raw reads directory should look like this:

```
raw_reads/
├── Sample1_R1.fastq.gz
├── Sample1_R2.fastq.gz
├── Sample2_R1.fastq.gz
├── Sample2_R2.fastq.gz
└── ...
```

**File naming conventions**:
- Paired-end: `*_R1.fastq.gz` and `*_R2.fastq.gz` (or `*_1.fastq.gz` / `*_2.fastq.gz`)
- Single-end: `*.fastq.gz`

#### Output Structure

```
results/
├── 01_QC/                    # Fastp quality control reports
├── 02_clean_reads/           # Quality-filtered reads
├── 03_assemblies/            # GetOrganelle assembly outputs
├── 04_mapping/               # BWA mapping results
├── 05_cp_reads/              # Extracted chloroplast reads
├── 06_reports/               # Summary statistics
└── 07_assembled_genomes/     # Final assembled genomes (use for Module 2)
```

#### Advanced Options

```bash
# Use more threads for faster processing
cgas --module 1 -i raw_reads/ -o results/ -t 16

# Skip already processed samples (useful when adding new data)
cgas --module 1 -i raw_reads/ -o results/ --skip-existing

# Custom k-mer values for difficult assemblies
cgas --module 1 -i raw_reads/ -o results/ -k "21,33,55,77,99"

# More GetOrganelle rounds for challenging samples
cgas --module 1 -i raw_reads/ -o results/ -R 20

# Force mapping even for incomplete assemblies
cgas --module 1 -i raw_reads/ -o results/ --force-mapping

# Run in background with nohup (for large datasets on servers)
cgas --module 1 -i large_dataset/ -o assemblies/ --use-nohup

# Complete example with all options
cgas --module 1 \
  -i raw_reads/ \
  -o results/ \
  -t 16 \
  -F embplant_pt \
  -k "21,45,65,85,105" \
  -R 15 \
  --trim-poly-g \
  --skip-existing
```

#### Direct Command

```bash
# Alternative direct command
cgas-assembly -i raw_reads/ -o results/ -t 16
```

#### Key Options Explained

| Option | Description | Default |
|--------|-------------|---------|
| `-t, --threads` | Number of CPU threads | 4 |
| `-k, --kmer` | K-mer values for assembly | "21,45,65,85,105" |
| `-R, --rounds` | GetOrganelle assembly rounds | 15 |
| `-F, --organelle` | Organelle genome type | embplant_pt |
| `--skip-existing` | Skip already processed samples | False |
| `--force-mapping` | Map reads even if assembly incomplete | False |
| `--trim-poly-g` | Remove poly-G tails (for NextSeq) | False |
| `--use-nohup` | Run in background with nohup | False |

#### Tips for Module 1

- **Keep output directory outside input directory** to avoid confusion
- **Use `--skip-existing`** when adding new samples to an existing project
- **Increase `-R` (rounds)** for difficult samples or low coverage
- **Use `--trim-poly-g`** for NextSeq or NovaSeq data with poly-G tails
- **Check `06_reports/`** for assembly quality summaries
- **Only complete assemblies** are copied to `07_assembled_genomes/`

---

### Module 2: Plastome Annotation

**Purpose**: Annotate assembled plastomes using PGA (Plastid Genome Annotator) with reference-based annotation.

**Input**: FASTA files from Module 1 (assembled genomes)  
**Output**: GenBank files with complete annotations

#### Basic Usage

```bash
# Using reference genomes directory
cgas --module 2 \
  -i 07_assembled_genomes/ \
  -r reference_genomes/ \
  --pga /path/to/PGA/PGA.pl
```

#### Input Requirements

1. **Assembled genomes**: From Module 1's `07_assembled_genomes/` directory
2. **Reference genomes**: Directory containing high-quality reference GenBank files
3. **PGA path**: Path to PGA.pl script

#### Output Structure

```
Module2_Annotations/
├── Sample1/
│   ├── Sample1.fasta         # Input genome
│   ├── Sample1.gb            # Annotated GenBank file
│   ├── Sample1_curation.txt  # Curation log
│   └── PGA_output/           # PGA raw output
├── Sample2/
│   └── ...
└── summary_report.txt        # Overall annotation summary
```

#### Advanced Options

```bash
# With organism name mapping
cgas --module 2 \
  -i 07_assembled_genomes/ \
  -r reference_genomes/ \
  --organism-file organisms.txt \
  --pga /path/to/PGA/PGA.pl

# Using specific output directory
cgas --module 2 \
  -i assemblies/ \
  -o my_annotations/ \
  -r custom_reference.gb \
  -t 16 \
  --pga /path/to/PGA/PGA.pl
```

#### Organism File Format

Create a text file (`organisms.txt`) mapping sample IDs to organism names:

```
Sample1	Species name 1
Sample2	Species name 2
Sample3	Species name 3
```

Format: `SampleID[TAB]Organism Name`

#### Direct Command

```bash
cgas-annotate -i assemblies/ -r reference/ --pga /path/to/PGA/PGA.pl
```

#### Important Notes

- **Reference quality is critical**: Use high-quality, well-annotated reference genomes
- **Manual curation**: The script automatically curates LSC/SSC regions based on IR regions identified by PGA
- **Alternative tools**: You can also annotate using external tools like Geneious or GeSeq, then use Module 3 for normalization
- **For NCBI submission**: Annotations should be carefully reviewed before submission

#### Key Options

| Option | Description | Required |
|--------|-------------|----------|
| `-i, --input` | Input directory with FASTA files | Yes |
| `-r, --reference` | Reference GenBank file(s) or directory | Yes |
| `--pga` | Path to PGA.pl script | Yes |
| `-o, --output` | Output directory | No (default: Module2_Annotations/) |
| `-t, --threads` | Number of threads | No (default: 4) |
| `--organism-file` | File mapping samples to organism names | No |

---

### Module 3: Gene Comparison & Normalization

**Purpose**: Standardize gene names and annotations across multiple chloroplast genomes for consistency.

**Input**: GenBank files (from Module 2 or downloaded from NCBI)  
**Output**: Normalized GenBank files with standardized gene names

#### Basic Usage

```bash
# All genomes should be in one directory
# Choose one as reference (use underscores, not spaces in filename)
cgas --module 3 -r reference_genome.gb
```

#### What Module 3 Does

- ✅ Standardizes gene names across different annotation conventions
- ✅ Detects missing or extra genes
- ✅ Validates product descriptions
- ✅ Shows intron presence/absence
- ✅ Compares coding region lengths
- ✅ Generates comprehensive comparison reports

#### File Naming Convention

**Important**: Reference filename must use underscores, not spaces:

```bash
# ✓ Correct
cgas --module 3 -r Abutilon_grandifolium.gb

# ✗ Wrong
cgas --module 3 -r Abutilon grandifolium.gb
```

#### Directory Structure

```
working_directory/
├── reference_genome.gb       # Your reference genome
├── genome1.gb               # Target genome 1
├── genome2.gb               # Target genome 2
├── genome3.gb               # Target genome 3
└── ...
```

Run Module 3 from this directory:

```bash
cd working_directory/
cgas --module 3 -r reference_genome.gb
```

#### Output Structure

```
Module3_GeneComparison/
├── normalized_genbank/          # Standardized GenBank files
│   ├── reference_genome.gb
│   ├── genome1.gb
│   ├── genome2.gb
│   └── ...
├── comparison_report.txt        # Detailed comparison
├── gene_presence_matrix.xlsx    # Gene presence/absence
├── missing_genes_report.xlsx    # Missing gene details
└── statistics_summary.txt       # Overall statistics
```

#### Direct Command

```bash
cgas-compare -r reference_genome.gb
```

#### Advanced Usage

```bash
# Specify target directory (if not in current directory)
cgas --module 3 -r reference.gb -t target_genomes/

# Custom output directory
cgas --module 3 -r reference.gb --output custom_output/
```

#### Use Cases

1. **After Module 2**: Normalize newly annotated genomes
2. **NCBI genomes**: Standardize downloaded genomes with inconsistent annotations
3. **Mixed sources**: Harmonize genomes from different annotation pipelines
4. **Quality control**: Identify annotation errors or missing genes

#### Tips for Module 3

- **Choose a high-quality reference**: Your reference should have complete, accurate annotations
- **Review comparison reports**: Check `comparison_report.txt` for discrepancies
- **Missing genes**: Check `missing_genes_report.xlsx` to identify genuinely missing genes vs. annotation errors
- **Use normalized output**: Always use GenBank files from `normalized_genbank/` for downstream modules

---

### Module 4: Format Conversion & NCBI Submission

**Purpose**: Validate genome annotations and prepare files for NCBI GenBank submission.

**Input**: Finalized GenBank files  
**Output**: FASTA, TBL files, and submission-ready formats

#### Basic Usage

```bash
# Run in directory containing finalized GenBank files
cd finalized_genomes/
cgas --module 4
```

#### What Module 4 Does

- ✅ Validates genome annotations
- ✅ Detects annotation errors
- ✅ Predicts potential submission issues
- ✅ Checks overall quality
- ✅ Generates FASTA files for NCBI
- ✅ Creates TBL files for submissions
- ✅ Combines files for batch submission

#### Output Structure

```
Module4_Format_Conversion/
├── fasta/                       # FASTA files for each genome
│   ├── genome1.fasta
│   ├── genome2.fasta
│   └── ...
├── tbl/                         # TBL annotation files
│   ├── genome1.tbl
│   ├── genome2.tbl
│   └── ...
├── combined/                    # Combined files for batch submission
│   ├── all_genomes.fsa
│   └── all_genomes.tbl
├── validation_report.txt        # Quality and error report
└── ncbi_submission_guide.txt   # Submission instructions
```

#### Direct Command

```bash
cgas-convert -i genomes/ -o ncbi_submission/
```

#### Advanced Options

```bash
# Specify input directory
cgas --module 4 -i finalized_genomes/ -o submission_files/

# With custom options
cgas-convert \
  -i genomes/ \
  -o ncbi_ready/ \
  --validate-only  # Only validate, don't convert
```

#### NCBI Submission Workflow

1. **Run Module 4** in your finalized genomes directory
2. **Review validation report** - fix any errors
3. **Use combined files** from `combined/` directory for batch submission
4. **Follow NCBI guidelines** in `ncbi_submission_guide.txt`

#### Important Notes

- **Apply separately**: For NCBI submission genomes, run Module 4 separately from comparative analysis genomes
- **Quality check**: Even if not submitting to NCBI, run this for quality validation
- **Fast processing**: Takes 1-2 minutes for ~20 genomes
- **Streamlines submission**: Replaces manual single-species uploads

#### Tips for Module 4

- **Fix errors first**: Address all issues in `validation_report.txt` before submission
- **Separate projects**: Keep NCBI submission files separate from comparative analysis
- **Batch submission**: Use combined files for easier submission of multiple genomes
- **Double-check**: Verify file formats match NCBI requirements

---

## Phase 2: Main Analysis Modules (5-14)

These modules perform comparative analysis on your finalized, normalized GenBank files. Most can be **run together** or individually.

### Setup for Phase 2

```bash
# 1. Create analysis directory
mkdir all_genomes/
cd all_genomes/

# 2. Copy normalized GenBank files here
cp /path/to/Module3_GeneComparison/normalized_genbank/*.gb .

# 3. Run analysis modules
# Important: Do NOT include outgroup genome in this directory 
# if you don't want it in comparative analysis
```

### Running Multiple Modules

```bash
# Run all comparative modules (5-13) at once
cgas --modules 5,6,7,8,9,10,11,12,13

# Run specific subsets
cgas --modules 5,6,7          # Gene and structure analysis
cgas --modules 8,9            # Codon and amino acid
cgas --modules 10,11,12,13    # SNP, intron, SSR, diversity

# Run phylogeny separately (Module 14)
cgas --module 14 --macse --genes-only --iqtree -og ../outgroup.gb
```

---

### Module 5: Gene Comparative Analysis

**Purpose**: Compare gene content across all species in your dataset.

**Input**: GenBank files  
**Output**: Comparative gene analysis tables

#### Basic Usage

```bash
cd all_genomes/
cgas --module 5
```

#### Output Structure

```
Module5_Gene_Comparative_Analysis/
├── gene_comparison_matrix.xlsx     # Gene presence/absence matrix
├── gene_statistics.xlsx            # Gene counts and statistics
├── unique_genes_per_species.xlsx   # Species-specific genes
├── core_genes.txt                  # Genes present in all species
├── variable_genes.txt              # Genes absent in some species
└── summary_report.txt              # Overall summary
```

#### Direct Command

```bash
cgas-gene-compare -i genomes/ -o gene_analysis/
```

#### What You Get

- **Gene presence/absence matrix**: Shows which genes are present in each species
- **Core genes**: Genes conserved across all species
- **Variable genes**: Genes with presence/absence variation
- **Unique genes**: Species-specific genes
- **Statistics**: Comprehensive gene content statistics

---

### Module 6: Gene Content Tables

**Purpose**: Generate publication-ready gene content tables in Word format.

**Input**: GenBank files  
**Output**: Formatted Word documents

#### Basic Usage

```bash
cd all_genomes/
cgas --module 6
```

#### Output Structure

```
Module6_Gene_Content_Tables/
├── Table1_Gene_List.docx           # Complete gene list with details
├── Table2_Gene_Categories.docx     # Genes by functional category
├── Table3_tRNA_rRNA.docx          # tRNA and rRNA genes
└── summary_statistics.xlsx         # Supporting statistics
```

#### Direct Command

```bash
cgas-gene-table -i genomes/ -o gene_tables/
```

#### What You Get

- **Publication-ready tables**: Formatted Word documents ready for manuscripts
- **Gene categories**: Genes organized by function (photosynthesis, translation, etc.)
- **RNA genes**: Separate tables for tRNA and rRNA genes
- **Consistent formatting**: Professional tables with proper headers and structure

---

### Module 7: Genome Structure Analysis

**Purpose**: Analyze and compare genome structure including LSC, SSC, IR regions, and GC content.

**Input**: GenBank files  
**Output**: Structural comparison tables and statistics

#### Basic Usage

```bash
cd all_genomes/
cgas --module 7
```

#### Output Structure

```
Module7_Genome_Structure_Analysis/
├── genome_structure_comparison.xlsx    # LSC/SSC/IR sizes
├── gc_content_analysis.xlsx           # GC% by region and gene
├── functional_gene_distribution.xlsx  # Genes by region and function
├── region_statistics.xlsx             # Detailed region stats
└── structure_summary.txt              # Overall summary
```

#### Direct Command

```bash
cgas-genome-compare -i genomes/ -o structure_analysis/
```

#### What You Get

- **Region sizes**: LSC, SSC, IR lengths for all species
- **GC content**: Overall and region-specific GC percentages
- **Functional distribution**: Gene distribution by region and function
- **Comparative statistics**: Structural variation across species

---

### Module 8: Codon Usage (RSCU)

**Purpose**: Analyze codon usage bias with RSCU (Relative Synonymous Codon Usage) calculations and visualization.

**Input**: GenBank files  
**Output**: RSCU tables and high-quality visualizations

**Requirements**: R with packages (ggplot2, seqinr, dplyr, tidyr, RColorBrewer, patchwork)

#### Basic Usage

```bash
cd all_genomes/
cgas --module 8
```

#### Output Structure

```
Module8_Codon_Usage/
├── rscu_values.xlsx               # RSCU values for all species
├── codon_frequency.xlsx           # Raw codon frequencies
├── rscu_heatmap.pdf              # Heatmap visualization
├── codon_preference.pdf          # Codon preference plot
└── statistics_summary.txt         # Summary statistics
```

#### Direct Command

```bash
cgas-codon -i genomes/ -o codon_analysis/
```

#### What You Get

- **RSCU values**: Relative synonymous codon usage for each species
- **Codon frequencies**: Absolute codon counts
- **High-quality plots**: Publication-ready PDF visualizations
- **Bias analysis**: Identification of codon usage preferences

#### Tips for Module 8

- **Requires R**: Ensure R and required packages are installed
- **Check plots**: Review PDFs for codon usage patterns
- **Compare across species**: Look for conserved vs. variable codon preferences

---

### Module 9: Amino Acid Analysis

**Purpose**: Analyze amino acid composition patterns across species with visualization.

**Input**: GenBank files  
**Output**: Amino acid composition tables and plots

**Requirements**: R with packages (ggplot2, dplyr, tidyr, scales)

#### Basic Usage

```bash
cd all_genomes/
cgas --module 9
```

#### Output Structure

```
Module9_Amino_Acid_Analysis/
├── amino_acid_composition.xlsx    # AA composition for all species
├── aa_frequency.xlsx              # Raw AA frequencies
├── aa_composition_plot.pdf        # Composition visualization
├── aa_comparison_heatmap.pdf      # Species comparison heatmap
└── statistics_summary.txt          # Summary statistics
```

#### Direct Command

```bash
cgas-amino -i genomes/ -o aa_analysis/
```

#### What You Get

- **AA composition**: Percentage of each amino acid
- **Frequency tables**: Raw counts of amino acids
- **Visualizations**: High-quality comparative plots
- **Pattern analysis**: Identification of compositional biases

---

### Module 10: SNP Analysis

**Purpose**: Detect and analyze single nucleotide polymorphisms (SNPs) across aligned genomes.

**Input**: GenBank files  
**Output**: SNP tables, substitution matrices, and visualizations

#### Basic Usage

```bash
cd all_genomes/
cgas --module 10
```

#### Output Structure

```
Module10_SNP_Analysis/
├── snp_summary.xlsx               # Overall SNP statistics
├── substitution_matrix.xlsx       # Substitution types and counts
├── snp_distribution.xlsx          # SNPs by gene and region
├── snp_density_plot.pdf          # SNP density visualization
├── substitution_heatmap.pdf      # Substitution pattern heatmap
└── variable_sites.txt             # List of all variable positions
```

#### Direct Command

```bash
cgas-snp -i genomes/ -o snp_analysis/
```

#### What You Get

- **SNP detection**: Identification of all variable positions
- **Substitution patterns**: Types and frequencies of substitutions
- **Gene-level SNPs**: SNPs mapped to specific genes
- **Quality plots**: Visualization of SNP distribution and patterns

---

### Module 11: Intron Analysis

**Purpose**: Analyze intron presence, absence, and structure in genes and tRNAs.

**Input**: GenBank files  
**Output**: Comprehensive intron comparison tables

#### Basic Usage

```bash
cd all_genomes/
cgas --module 11
```

#### Output Structure

```
Module11_Intron_Analysis/
├── gene_introns_comparison.xlsx   # Introns in protein-coding genes
├── trna_introns_comparison.xlsx   # Introns in tRNA genes
├── intron_statistics.xlsx         # Overall intron stats
├── intron_presence_matrix.xlsx    # Intron presence/absence
└── summary_report.txt              # Summary of findings
```

#### Direct Command

```bash
cgas-intron -i genomes/ -o intron_analysis/
```

#### What You Get

- **Intron detection**: Complete intron annotation for all genes
- **Gene-level analysis**: Intron counts and positions per gene
- **tRNA introns**: Separate analysis for tRNA introns
- **Comparative tables**: Intron presence/absence across species

---

### Module 12: SSR Analysis

**Purpose**: Detect microsatellites (SSRs), assign to genomic/functional regions, compare motifs, and generate visualizations.

**Input**: GenBank files  
**Output**: SSR tables, motif comparisons, and high-quality plots

**Requirements**: R with packages (ggplot2, dplyr, tidyr, RColorBrewer, patchwork, ggrepel, scales)

#### Basic Usage

```bash
cd all_genomes/
cgas --module 12
```

#### Output Structure

```
Module12_SSR_Analysis/
├── ssr_summary.xlsx               # Overall SSR statistics
├── ssr_by_motif.xlsx             # SSRs grouped by motif type
├── ssr_by_region.xlsx            # SSRs by genomic region
├── ssr_distribution.xlsx          # SSR density and distribution
├── ssr_comparison_plot.pdf        # Comparative visualization
├── ssr_motif_heatmap.pdf         # Motif frequency heatmap
└── ssr_regions_barplot.pdf        # SSRs by region and function
```

#### Direct Command

```bash
cgas-ssr -i genomes/ -o ssr_analysis/
```

#### What You Get

- **SSR detection**: All microsatellites identified and characterized
- **Motif analysis**: SSRs grouped by repeat motif
- **Regional assignment**: SSRs mapped to LSC/SSC/IR and functional regions
- **High-quality plots**: Publication-ready visualizations

#### Tips for Module 12

- **Check motif types**: Common plastid SSRs include A/T mononucleotides
- **Regional patterns**: Compare SSR density across genomic regions
- **Functional relevance**: Note SSRs in coding vs. non-coding regions

---

### Module 13: Nucleotide Diversity

**Purpose**: Calculate nucleotide diversity (π) across the genome with high-quality visualizations.

**Input**: GenBank files  
**Output**: Diversity statistics, sliding window analysis, and plots

**Requirements**: R with packages (ggplot2, dplyr, tidyr, cowplot, gridExtra)

#### Basic Usage

```bash
cd all_genomes/
cgas --module 13
```

#### Output Structure

```
Module13_Nucleotide_Diversity/
├── nucleotide_diversity.xlsx      # π values by gene and region
├── sliding_window_pi.xlsx         # Sliding window analysis
├── diversity_by_region.xlsx       # π for LSC/SSC/IR
├── diversity_plot.pdf             # Genome-wide π visualization
├── sliding_window_plot.pdf        # Sliding window plot
└── hotspot_regions.txt            # High diversity regions
```

#### Direct Command

```bash
cgas-diversity -i genomes/ -o diversity_analysis/
```

#### What You Get

- **Nucleotide diversity (π)**: Calculated for genes, regions, and genome-wide
- **Sliding window analysis**: Fine-scale diversity patterns
- **Hotspot identification**: Regions of high genetic diversity
- **Visualizations**: High-quality plots showing diversity patterns

#### Tips for Module 13

- **Look for hotspots**: High π regions may be under positive selection
- **Compare regions**: IR typically shows lower diversity than LSC/SSC
- **Gene-level analysis**: Identify rapidly evolving genes

---

### Module 14: Phylogenetic Analysis

**Purpose**: Build gene matrices and reconstruct phylogenetic trees using IQ-TREE with optional MACSE codon-aware alignment.

**Input**: GenBank files + outgroup  
**Output**: Aligned matrices, phylogenetic trees, and support values

**Requirements**: MAFFT, IQ-TREE, MACSE (optional), Java (for MACSE)

#### Basic Usage

```bash
# Using MACSE (codon-aware alignment, recommended for coding genes)
cgas --module 14 --macse --genes-only --iqtree -og ../outgroup.gb

# Using MAFFT (standard alignment)
cgas --module 14 --genes-only --iqtree -og ../outgroup.gb
```

#### Alignment Options

**Option 1: Protein-coding genes only (Recommended)**
```bash
cgas --module 14 --macse --genes-only --iqtree -og outgroup.gb
```

**Option 2: Complete genome sequences**
```bash
cgas --module 14 --complete-genome --iqtree -og outgroup.gb
```

**Option 3: Specific genes**
```bash
cgas --module 14 --genes matK,rbcL,atpB,ndhF --macse --iqtree -og outgroup.gb
```

**Option 4: Without outgroup (comparative analysis only)**
```bash
cgas --module 14 --genes-only --macse --iqtree
```

#### Output Structure

```
Module14_Phylogeny/
├── alignments/
│   ├── gene1_aligned.fasta
│   ├── gene2_aligned.fasta
│   └── ...
├── concatenated_matrix.fasta      # Complete concatenated alignment
├── concatenated_matrix.phy        # PHYLIP format
├── partition_file.txt             # Gene partition information
├── tree.treefile                  # Best ML tree (from IQ-TREE)
├── tree.contree                   # Consensus tree
├── tree.log                       # IQ-TREE log file
├── bootstrap_support.txt          # Bootstrap values
└── phylogeny_summary.txt          # Summary of analysis
```

#### Direct Command

```bash
cgas-phylogeny --macse --genes-only --iqtree -og outgroup.gb
```

#### Advanced Options

```bash
# Multiple outgroups
cgas --module 14 --macse --genes-only --iqtree \
  -og outgroup1.gb,outgroup2.gb

# Complete matrix only (no missing data)
cgas --module 14 --genes-only --iqtree --complete-only

# Custom output format
cgas --module 14 --genes-only --iqtree --format nexus

# With custom parameters
cgas-phylogeny \
  --genes-only \
  --macse \
  --iqtree \
  -og outgroup.gb \
  --complete-only \
  --format phylip \
  -o custom_phylogeny/
```

#### Key Options

| Option | Description |
|--------|-------------|
| `--genes-only` | Use only protein-coding genes |
| `--complete-genome` | Use complete genome sequences |
| `--genes GENE1,GENE2` | Use specific genes (comma-separated) |
| `--macse` | Use MACSE for codon-aware alignment (recommended) |
| `--iqtree` | Run IQ-TREE for phylogenetic inference |
| `-og, --outgroup` | Outgroup GenBank file(s) |
| `--complete-only` | Only use complete matrix (no missing data) |
| `--format` | Output format (fasta, phylip, nexus) |

#### Important Notes

- **Outgroup file naming**: Use underscores, not spaces (e.g., `Abutilon_grandifolium.gb`)
- **Outgroup location**: Keep outgroup files **separate** from main analysis directory
- **MACSE recommended**: For protein-coding genes, MACSE provides better codon-aware alignment
- **Missing data**: By default, includes all genes; use `--complete-only` for complete matrix
- **Multiple outgroups**: Separate with commas: `-og outgroup1.gb,outgroup2.gb`

#### Tips for Module 14

- **Choose appropriate genes**: matK, rbcL, atpB, ndhF are common phylogenetic markers
- **Use MACSE for coding genes**: Maintains codon structure in alignment
- **Check alignment quality**: Review alignments in `alignments/` directory
- **Bootstrap support**: Check `bootstrap_support.txt` for node confidence
- **Tree visualization**: Use FigTree, iTOL, or R packages to visualize trees

---

## Jupyter Notebook Usage

CGAS can be used in Jupyter notebooks for interactive analysis.

### Setup

```python
# Install dependencies (if not using conda)
%pip install biopython>=1.79 pandas>=1.3.0 openpyxl>=3.1 numpy>=1.24 python-docx>=0.8.11
```

### Running Modules

#### Method 1: Using %run Magic Command

```python
# Module 3: Gene comparison
%run cgas_module3.py -r reference.gb

# Module 5: Gene comparative analysis
%run cgas_module5.py

# Module 8: Codon usage
%run cgas_module8.py

# Module 14: Phylogeny
%run cgas_module14.py --genes-only --macse --iqtree -og outgroup.gb
```

#### Method 2: Using ! Shell Commands

```python
# List modules
!cgas --list

# Run specific module
!cgas --module 5

# Run with options
!cgas --module 14 --macse --genes-only --iqtree -og ../outgroup.gb
```

#### Method 3: Programmatic Access

```python
import subprocess

# Run CGAS module
result = subprocess.run(
    ['cgas', '--module', '5'],
    capture_output=True,
    text=True
)

print(result.stdout)
```

### Loading and Analyzing Results

```python
import pandas as pd
from Bio import SeqIO

# Load CGAS output
gene_comparison = pd.read_excel('Module5_Gene_Comparative_Analysis/gene_comparison_matrix.xlsx')
print(gene_comparison.head())

# Load GenBank files
genomes = list(SeqIO.parse('genome.gb', 'genbank'))
print(f"Loaded {len(genomes)} genome(s)")

# Custom analysis on CGAS results
# ... your analysis code ...
```

### Important Notes for Jupyter

- **Modules 1 and 2**: Require specialized software; Linux or macOS recommended
- **WSL2 on Windows**: May work but not guaranteed (especially GetOrganelle)
- **Modules 3-14**: Work reliably in Jupyter with proper dependencies
- **External tools**: Ensure MAFFT, IQ-TREE, R are in system PATH

---

## Direct Module Commands

Each CGAS module has a direct command for convenient access:

```bash
cgas-assembly          # Module 1: Genome Assembly
cgas-annotate          # Module 2: Plastome Annotation
cgas-compare           # Module 3: Gene Comparison
cgas-convert           # Module 4: Format Conversion
cgas-gene-compare      # Module 5: Gene Comparative Analysis
cgas-gene-table        # Module 6: Gene Content Tables
cgas-genome-compare    # Module 7: Genome Structure Analysis
cgas-codon             # Module 8: Codon Usage (RSCU)
cgas-amino             # Module 9: Amino Acid Analysis
cgas-snp               # Module 10: SNP Analysis
cgas-intron            # Module 11: Intron Analysis
cgas-ssr               # Module 12: SSR Analysis
cgas-diversity         # Module 13: Nucleotide Diversity
cgas-phylogeny         # Module 14: Phylogenetic Analysis
```

### Examples

```bash
# Module 5
cgas-gene-compare -i genomes/ -o gene_analysis/

# Module 8
cgas-codon -i genomes/ -o codon_analysis/

# Module 14
cgas-phylogeny --macse --genes-only --iqtree -og outgroup.gb
```

---

## Best Practices

### 1. Directory Organization

```
project/
├── raw_reads/              # Original FASTQ files
├── assemblies/             # Module 1 output
├── annotations/            # Module 2 output
├── normalized/             # Module 3 output
├── analysis/               # Modules 5-14
│   ├── all_genomes/        # GenBank files for analysis
│   ├── Module5_*/
│   ├── Module6_*/
│   └── ...
└── references/
    ├── outgroup.gb
    └── reference_genomes/
```

### 2. Workflow Sequence

**Complete workflow**:
```bash
# Phase 1: Preparation
cgas --module 1 -i raw_reads/ -o assemblies/          # Assembly
cgas --module 2 -i assemblies/ -r references/          # Annotation
cgas --module 3 -r reference.gb                        # Normalization
cgas --module 4                                        # Validation (optional)

# Phase 2: Analysis
cd analysis/all_genomes/
cgas --modules 5,6,7,8,9,10,11,12,13                  # Comparative analysis
cgas --module 14 --macse --genes-only --iqtree -og outgroup.gb  # Phylogeny
```

### 3. Quality Control

- **After Module 1**: Check `06_reports/` for assembly quality
- **After Module 2**: Verify annotations in GenBank files
- **After Module 3**: Review comparison reports for discrepancies
- **Before Module 14**: Ensure all GenBank files are properly formatted

### 4. File Naming

- **Use underscores**, not spaces: `Species_name.gb` not `Species name.gb`
- **Be consistent**: Same naming convention throughout
- **Avoid special characters**: Stick to letters, numbers, underscores

### 5. Resource Management

- **Use appropriate threads**: `-t` option (don't exceed available cores)
- **Monitor disk space**: Assemblies and alignments can be large
- **Use `--skip-existing`**: For Module 1 when adding new samples
- **Background processing**: Use `--use-nohup` for long-running jobs

---

## Tips and Tricks

### Parallel Processing

```bash
# Run multiple independent modules in background
cgas --module 5 &
cgas --module 6 &
cgas --module 7 &
wait  # Wait for all to complete
```

### Batch Processing

```bash
# Process multiple reference genomes
for ref in reference_genomes/*.gb; do
    cgas --module 3 -r "$ref"
done
```

### Combining Outputs

```bash
# Merge results from multiple runs
cat Module5_*/summary_report.txt > combined_summary.txt
```

### Quick Quality Checks

```bash
# Count genes per genome
for gb in *.gb; do
    echo "$gb: $(grep -c "     CDS     " $gb) CDSs"
done

# Check genome sizes
for gb in *.gb; do
    python -c "from Bio import SeqIO; print('$gb:', len(list(SeqIO.parse('$gb', 'genbank'))[0]))"
done
```

### Handling Large Datasets

```bash
# Use screen or tmux for long-running jobs
screen -S cgas_assembly
cgas --module 1 -i large_dataset/ -o results/ --use-nohup
# Ctrl+A, D to detach
# screen -r cgas_assembly to reattach
```

---

## Common Workflows

### Workflow 1: Complete Analysis from Scratch

```bash
# 1. Assembly
cgas --module 1 -i raw_reads/ -o assemblies/ -t 16

# 2. Annotation
cgas --module 2 -i assemblies/07_assembled_genomes/ \
    -r reference_genomes/ --pga /path/to/PGA/PGA.pl

# 3. Normalization
cd Module2_Annotations/
cgas --module 3 -r reference_species.gb

# 4. Analysis
mkdir ../analysis/all_genomes/
cp Module3_GeneComparison/normalized_genbank/*.gb ../analysis/all_genomes/
cd ../analysis/all_genomes/

# 5. Comparative analysis
cgas --modules 5,6,7,8,9,10,11,12,13

# 6. Phylogeny
cgas --module 14 --macse --genes-only --iqtree \
    -og ../../references/outgroup.gb
```

### Workflow 2: Analysis of NCBI Genomes

```bash
# 1. Download genomes from NCBI (use browser or EDirect)
mkdir ncbi_genomes/
cd ncbi_genomes/
# ... download GenBank files ...

# 2. Normalize annotations
cgas --module 3 -r best_annotated_genome.gb

# 3. Analysis
mkdir ../analysis/
cp Module3_GeneComparison/normalized_genbank/*.gb ../analysis/
cd ../analysis/

# 4. Run all analyses
cgas --modules 5,6,7,8,9,10,11,12,13
cgas --module 14 --macse --genes-only --iqtree -og ../outgroup.gb
```

### Workflow 3: Quick Phylogenetic Analysis

```bash
# 1. Place all GenBank files in one directory
mkdir phylogeny_analysis/
cp normalized_genomes/*.gb phylogeny_analysis/
cd phylogeny_analysis/

# 2. Run phylogeny
cgas --module 14 --macse --genes-only --iqtree -og ../outgroup.gb

# 3. View tree
# Use FigTree, iTOL, or:
cat Module14_Phylogeny/tree.treefile
```

### Workflow 4: Targeted Gene Analysis

```bash
# 1. Extract specific genes for analysis
cgas --module 14 --genes matK,rbcL,atpB,ndhF \
    --macse --iqtree -og outgroup.gb

# 2. Run codon usage on same genes
cgas --module 8

# 3. Check SNPs in these genes
cgas --module 10
```

---

## Troubleshooting

### Module 1 Issues

**Problem**: GetOrganelle fails to assemble  
**Solution**:
```bash
# Try more rounds
cgas --module 1 -i reads/ -o results/ -R 20

# Try different k-mers
cgas --module 1 -i reads/ -o results/ -k "21,33,45,55,77,99"

# Check input quality
fastp -i sample_R1.fastq.gz -o clean_R1.fastq.gz
```

**Problem**: Low coverage  
**Solution**: Check if you're using the correct organelle type (`-F embplant_pt`)

### Module 2 Issues

**Problem**: PGA not found  
**Solution**:
```bash
# Verify PGA path
ls -la /path/to/PGA/PGA.pl
which perl

# Test PGA
perl /path/to/PGA/PGA.pl
```

**Problem**: Poor annotation quality  
**Solution**: Use a higher-quality reference genome; check reference annotations manually

### Module 3 Issues

**Problem**: Many missing genes  
**Solution**: Check if genomes are from the same organelle type; verify reference quality

**Problem**: Gene name conflicts  
**Solution**: Module 3 handles this automatically; check `comparison_report.txt`

### Modules 5-13 Issues

**Problem**: No GenBank files found  
**Solution**: Ensure you're in the correct directory with `.gb` files

**Problem**: R visualization fails  
**Solution**:
```bash
# Install R packages
R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr'), repos='https://cran.rstudio.com/')"
```

### Module 14 Issues

**Problem**: MACSE fails  
**Solution**:
```bash
# Check Java
java -version  # Need Java 8+

# Test MACSE
java -jar /path/to/macse.jar -help
```

**Problem**: IQ-TREE not found  
**Solution**:
```bash
# Check installation
which iqtree
which iqtree2

# Add to PATH if needed
export PATH=$PATH:/path/to/iqtree/bin
```

**Problem**: Alignment quality issues  
**Solution**: Try MAFFT instead of MACSE, or manually review alignments

### General Issues

**Problem**: Out of memory  
**Solution**: Reduce number of threads (`-t`) or process fewer genomes at once

**Problem**: Permission denied  
**Solution**:
```bash
chmod +x cgas
chmod 755 output_directory/
```

**Problem**: Module output not found  
**Solution**: Check current directory; modules create output in current location by default

---

## Getting More Help

### Documentation

- **Main README**: [README.md](README.md)
- **Installation Guide**: [INSTALL.md](INSTALL.md)
- **Module-specific READMEs**: Check each module's directory

### Support Channels

- **GitHub Issues**: [Report bugs or ask questions](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS/issues)
- **Module Documentation**: Each module folder contains detailed README files

### Useful Commands

```bash
# View module help
cgas --module X --help

# Check CGAS version
cgas --version

# List all modules with descriptions
cgas --list

# View command history
history | grep cgas
```

---

## Next Steps

After reading this guide:

1. ✅ **Start with Phase 1** if you have raw reads
2. ✅ **Or start with Phase 2** if you have GenBank files
3. ✅ **Follow the complete workflow** for best results
4. ✅ **Check module-specific READMEs** for advanced options
5. ✅ **Join the community** on GitHub for support

---

**Happy analyzing! 🧬**

*For questions or issues, visit the [GitHub repository](https://github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS)*

*Last updated: January 2026*
