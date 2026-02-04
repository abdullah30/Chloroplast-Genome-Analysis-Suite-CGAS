# CGAS Module 5: Chloroplast Genome Gene Comparative Analysis
## Complete Documentation and User Guide

---

## 📋 Table of Contents
1. [Introduction](#introduction)
2. [Quick Start - Just Run It!](#quick-start---just-run-it)
4. [Command-Line Usage with Examples](#command-line-usage-with-examples)
5. [Jupyter Notebook Usage](#jupyter-notebook-usage)
6. [Input Requirements](#input-requirements)
7. [Output Structure](#output-structure)
8. [Detailed Feature Explanation](#detailed-feature-explanation)
9. [Gene Name Normalization](#gene-name-normalization)
10. [Pseudogene Detection](#pseudogene-detection)
11. [IR Duplication Analysis](#ir-duplication-analysis)
12. [Troubleshooting](#troubleshooting)
13. [Examples](#examples)
14. [FAQ](#faq)
15. [Technical Specifications](#technical-specifications)
16. [References](#references)

---

## Introduction

**CGAS Module 5** is a comprehensive tool for comparative analysis of chloroplast genome gene content across multiple species. This module identifies and counts genes, classifies them as functional or pseudogenes, detects IR-mediated duplications, and normalizes gene names across different annotation conventions to produce standardized comparative reports.

This module performs advanced gene analysis with:
- **Gene identification**: Extracts all genes from GenBank files
- **Name normalization**: Standardizes gene names across different conventions
- **Functional classification**: Distinguishes functional genes from pseudogenes
- **Duplication detection**: Identifies IR-mediated gene duplications
- **Comparative statistics**: Generates summary tables across all genomes
- **Normalization tracking**: Logs all name changes for transparency

### Key Features:
- **Comprehensive gene detection**: Finds all gene types (CDS, tRNA, rRNA)
- **Intelligent normalization**: Handles various gene naming conventions
- **Pseudogene identification**: Multiple criteria for pseudogene detection
- **IR duplication analysis**: Detects and characterizes gene duplications
- **Batch processing**: Analyzes multiple genomes simultaneously
- **Excel reporting**: Publication-ready spreadsheets with detailed statistics

### Scientific Applications:
- **Comparative genomics**: Compare gene content across species
- **Phylogenetic analysis**: Use gene presence/absence as characters
- **Evolutionary studies**: Track gene loss, gain, and pseudogenization
- **Taxonomic research**: Identify lineage-specific gene patterns
- **Genome annotation**: Standardize gene names across datasets
- **Publication preparation**: Generate comprehensive gene content tables

### Part of CGAS (Chloroplast Genome Analysis Suite)
CGAS is a comprehensive suite of tools designed for chloroplast genome analysis, providing seamless integration across different analysis modules.

---

## Quick Start - Just Run It!

### ⚡ One Command Does It All!

```bash
# Option 1: Using cgas command
cgas --module 5

# Option 2: Using cgas-gene-compare shortcut
cgas-gene-compare

# Option 3: Using python directly
python cgas_module5.py
```

**What happens when you run this:**
1. ✅ Finds all GenBank files in current directory
2. ✅ Extracts and normalizes all gene names
3. ✅ Classifies genes as functional or pseudogenes
4. ✅ Detects IR-mediated duplications
5. ✅ Generates comparative Excel reports
6. ✅ Creates normalization log file

### 📁 Real-World Folder Examples

#### Example 1: Basic Usage (Easiest!)
```bash
# Your folder structure:
/home/abdullah/gene_comparison/
├── chloroplast_genomes/
│   ├── Arabidopsis_thaliana.gb
│   ├── Oryza_sativa.gb
│   └── Zea_mays.gb

# Navigate to your folder
cd /home/abdullah/gene_comparison/chloroplast_genomes/

# Run the module (no arguments needed!)
python cgas_module5.py

# Output created automatically:
# Module5_Gene_Comparative_Analysis/
# ├── Chloroplast_Gene_Analysis_20260115_143022.xlsx
# └── Gene_Normalization_Log.xlsx
```

#### Example 2: Specify Input Folder
```bash
# Process genomes from specific directory
python cgas_module5.py -i genbank_files/

# Output created in:
# genbank_files/Module5_Gene_Comparative_Analysis/
```

#### Example 3: Custom Input and Output
```bash
# Full control over directories
python cgas_module5.py -i chloroplast_genomes/ -o Gene_Results/

# Input from: chloroplast_genomes/
# Output to: Gene_Results/
```

### 🎯 Quick Start Cheat Sheet

```bash
# ====== cgas command ======
cgas --module 5                                  # Process current directory
cgas --module 5 -i genbank_files/                # Specify input folder
cgas --module 5 -i data/ -o results/             # Custom input and output

# ====== cgas-gene-compare shortcut ======
cgas-gene-compare                                # Process current directory
cgas-gene-compare -i genbank_files/              # Specify input folder
cgas-gene-compare -i data/ -o results/           # Custom input and output

# ====== python command ======
python cgas_module5.py                           # Process current directory
python cgas_module5.py -i genbank_files/         # Specify input folder
python cgas_module5.py -i data/ -o results/      # Custom input and output
python cgas_module5.py --help                    # Get help
```

### 📊 What You Get (Output Files)

```
Module5_Gene_Comparative_Analysis/                 # Created automatically
├── 📊 Comparative Analysis
│   ├── Chloroplast_Gene_Analysis_YYYYMMDD_HHMMSS.xlsx  # Main analysis
│   └── Gene_Normalization_Log.xlsx                    # Normalization tracking
```

---
## Command-Line Usage with Examples

### Complete Command Reference

```bash
# ====================================================================
# USING cgas COMMAND
# ====================================================================

# Process current directory
cgas --module 5

# Specify input directory
cgas --module 5 -i /home/abdullah/chloroplast_genomes/

# Custom input and output directories
cgas --module 5 -i chloroplast_genomes/ -o Gene_Results/
```

```bash
# ====================================================================
# USING cgas-gene-compare SHORTCUT
# ====================================================================

# Process current directory
cgas-gene-compare

# With specific input and output
cgas-gene-compare -i chloroplast_genomes/ -o Gene_Results/
```

```bash
# ====================================================================
# USING python COMMAND
# ====================================================================

# ====================================================================
# BASIC USAGE
# ====================================================================

# 1. Process current directory (simplest!)
python cgas_module5.py

# 2. Specify input folder
python cgas_module5.py -i genbank_files/

# 3. Custom input and output folders
python cgas_module5.py -i chloroplast_genomes/ -o Gene_Results/

# 4. Get help
python cgas_module5.py --help


# ====================================================================
# REAL-WORLD EXAMPLES
# ====================================================================

# Example 1: Process chloroplast genomes
python cgas_module5.py -i chloroplast_genomes/

# Example 2: Save to specific output folder
python cgas_module5.py -i data/ -o ../results/gene_comparison/

# Example 3: Windows users
python cgas_module5.py -i "C:\Users\abdullah\genbank_files"
```

### Parameter Details

| Parameter | Short | Default | Description |
|-----------|-------|---------|-------------|
| `--input` | `-i` | `.` (current) | Input directory with GenBank files |
| `--output` | `-o` | `Module5_Gene_Comparative_Analysis` | Output directory for results |

---

## Jupyter Notebook Usage

> **Note:** Module 5 is a pure Python script, so both `%run` and the `!` shell operator work. Use `%run` for direct script execution, or `!` for the `cgas` / `cgas-gene-compare` commands. Run the cell from the directory containing your GenBank files, or pass `-i` to specify the input path.

```python
# Using %run (executes the script directly)
%run cgas_module5.py
%run cgas_module5.py -i chloroplast_genomes/
%run cgas_module5.py -i data/ -o results/

# Using ! operator with cgas command
!cgas --module 5
!cgas --module 5 -i chloroplast_genomes/

# Using ! operator with cgas-gene-compare shortcut
!cgas-gene-compare
!cgas-gene-compare -i chloroplast_genomes/ -o Gene_Results/

# Using ! operator with python
!python cgas_module5.py
!python cgas_module5.py -i chloroplast_genomes/ -o Gene_Results/
```

---

## Input Requirements

### Supported File Formats
GenBank files:
- `.gb`
- `.gbf`
- `.gbk`
- `.genbank`

### Critical Requirements

**1. Complete GenBank Files**
```gb
LOCUS       NC_000932               154478 bp DNA     circular PLN 15-JUL-2022
DEFINITION  Arabidopsis thaliana chloroplast, complete genome.
ACCESSION   NC_000932
VERSION     NC_000932.1
...
FEATURES             Location/Qualifiers
     source          1..154478
                     /organism="Arabidopsis thaliana"
                     /organelle="plastid"
                     /mol_type="genomic DNA"
     gene            864..2690
                     /gene="rbcL"
     CDS             864..2690
                     /gene="rbcL"
                     /product="ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"
     gene            5000..5200
                     /gene="trnI-GAU"
     tRNA            5000..5200
                     /gene="trnI-GAU"
                     /product="tRNA-Ile (GAU)"
     gene            6000..6200
                     /gene="ycf1"
                     /pseudo
     CDS             6000..6200
                     /gene="ycf1"
                     /pseudo
ORIGIN
        1 atggcgacga cgttcgtcgt cgtttgtcga tctcgtctga cttcagcctg atcggtagca
...
//
```

**2. Gene Features Required**
- Files must contain `gene` features
- Each gene should have proper `/gene` qualifier
- Optional: `CDS`, `tRNA`, `rRNA` features for classification

**3. Various Naming Conventions**
- Module handles different annotation styles automatically
- Gene names are normalized to standard conventions

### File Organization

```bash
genbank_files/
├── Arabidopsis_thaliana.gb          # Chloroplast genome 1
├── Oryza_sativa.gb                  # Chloroplast genome 2
├── Zea_mays.gb                      # Chloroplast genome 3
└── Nicotiana_tabacum.gb             # Chloroplast genome 4
```

### Example Gene Name Variations

**Different annotation styles handled:**
```gb
# Standard names
gene            /gene="rbcL"
gene            /gene="matK"
gene            /gene="trnI-GAU"

# Alternative names
gene            /gene="pafI"          # Will normalize to ycf3
gene            /gene="lhbA"          # Will normalize to psbZ
gene            /gene="TRNI-GAU"       # Will normalize to trnI-GAU
gene            /gene="rrn16S"         # Will normalize to rrn16
```

**What gets analyzed:**
- All gene features are extracted
- Names are normalized to standard conventions
- Genes classified as functional or pseudogenes
- Duplications identified and characterized

---

## Output Structure

### Directory Organization

```
Module5_Gene_Comparative_Analysis/                 # Main output folder
├── 📊 Comparative Analysis
│   ├── Chloroplast_Gene_Analysis_YYYYMMDD_HHMMSS.xlsx  # Main analysis
│   └── Gene_Normalization_Log.xlsx                    # Normalization tracking
```

### Key Output Files Explained

#### 1. Chloroplast_Gene_Analysis_YYYYMMDD_HHMMSS.xlsx

**Multiple worksheets:**

**Sheet 1: Gene_Table**
```
┌─────────────────────┬─────────────┬───────┬─────────────────────┬─────────────────────┐
│ Species            │ Gene Name   │ Copies│ Functional/Pseudo │ Duplicate Status    │
├─────────────────────┼─────────────┼───────┼─────────────────────┼─────────────────────┤
│ Arabidopsis_thaliana│ rbcL        │ 1     │ functional         │ single copy         │
│ Arabidopsis_thaliana│ trnI-GAU    │ 2     │ functional         │ 2 functional copies │
│ Arabidopsis_thaliana│ ycf1        │ 1     │ pseudogene         │ single copy         │
│ Arabidopsis_thaliana│ rps12       │ 2     │ functional         │ trans-spliced       │
└─────────────────────┴─────────────┴───────┴─────────────────────┴─────────────────────┘
```

**Sheet 2: Summary_Statistics**
```
┌─────────────────────┬───────────────────┬──────────┬─────────────┬─────────────────────┐
│ Genome             │ Total Unique Genes│ Functional│ Pseudogenes │ IR-duplicated Genes │
├─────────────────────┼───────────────────┼──────────┼─────────────┼─────────────────────┤
│ Arabidopsis_thaliana│ 79               │ 77       │ 2           │ 18 [includes ycf1] │
│ Oryza_sativa       │ 78               │ 76       │ 2           │ 17                   │
│ Zea_mays           │ 80               │ 79       │ 1           │ 19                   │
└─────────────────────┴───────────────────┴──────────┴─────────────┴─────────────────────┘
```

**Sheet 3: Genome_Specific_Genes**
```
┌─────────────────────┬─────────────┬───────┬─────────────────────┬─────────────────────┐
│ Species            │ Gene Name   │ Copies│ Functional/Pseudo │ Duplicate Status    │
├─────────────────────┼─────────────┼───────┼─────────────────────┼─────────────────────┤
│ Arabidopsis_thaliana│ accD        │ 1     │ functional         │ single copy         │
│ Oryza_sativa       │ ndhF        │ 1     │ functional         │ single copy         │
└─────────────────────┴─────────────┴───────┴─────────────────────┴─────────────────────┘
```

**Sheet 4: Missing_Genes_Report**
```
┌─────────────┬─────────────┬──────────┬─────────────────────────────┬─────────────┐
│ Gene Name   │ Present_In  │ Missing_From│ Present_In_Genomes         │ Status      │
├─────────────┼─────────────┼──────────┼─────────────────────────────┼─────────────┤
│ accD        │ 1           │ 2        │ Arabidopsis_thaliana       │ Unique      │
│ ndhF        │ 2           │ 1        │ Arabidopsis, Oryza        │ Variable    │
│ rbcL        │ 3           │ 0        │ All genomes                │ Universal  │
└─────────────┴─────────────┴──────────┴─────────────────────────────┴─────────────┘
```

#### 2. Gene_Normalization_Log.xlsx

**Tracks all name changes:**
```
┌─────────────────────┬─────────────┬───────────────────┬──────────────────┬─────────┐
│ Species            │ Original    │ Normalized        │ Normalization   │ Warning │
│                     │ Name        │ Name              │ Type            │         │
├─────────────────────┼─────────────┼───────────────────┼──────────────────┼─────────┤
│ Arabidopsis_thaliana│ pafI        │ ycf3              │ gene synonym    │         │
│ Arabidopsis_thaliana│ TRNI-GAU    │ trnI-GAU          │ tRNA format     │         │
│ Arabidopsis_thaliana│ rrn16S      │ rrn16             │ rRNA case       │         │
└─────────────────────┴─────────────┴───────────────────┴──────────────────┴─────────┘
```

---

## Detailed Feature Explanation

### 1. Gene Normalization

- tRNA names: `trnX-YYY` format (uppercase anticodon)
- rRNA names: `rrnXX` format (lowercase, no trailing 'S')
- Protein genes: Apply synonym mapping
- Track all changes for transparency

### 2. Pseudogene Detection
**Detection criteria:**
- Explicit `/pseudo` qualifier
- "pseudo" or "pseudogene" in product/note
- Gene feature without corresponding CDS
- Internal stop codons

### 3. IR Duplication Analysis

**Characterization:**
- Identifies genes in IRa, IRb, or both
- Distinguishes true duplications from single copies
- Special handling for trans-spliced genes (rps12)

### 4. Comparative Statistics

**Statistics generated:**
- Total unique functional genes
- Number of functional vs pseudogenes
- IR-duplicated genes with detailed notes
- Genome-specific genes
- Nearly universal genes

---

## Gene Name Normalization

### Why Normalization is Needed

Different annotation pipelines use different naming conventions:
- **tRNA**: `trnI-GAU`, `TRNI-GAU`, `trnI_GAU`, `trnI(GAU)`
- **rRNA**: `rrn16`, `rrn16S`, `RRN16`
- **Protein genes**: `pafI` vs `ycf3`, `lhbA` vs `psbZ`

**Examples:**
- `trnI-GAU` → `trnI-GAU` (already correct)
- `TRNI-GAU` → `trnI-GAU`
- `trnI_GAU` → `trnI-GAU`
- `trnI(GAU)` → `trnI-GAU`
- `trnI-GUU` → `trnI-GUU` (DNA to RNA)

**2. rRNA Genes**
**Examples:**
- `rrn16` → `rrn16`
- `rrn16S` → `rrn16`
- `RRN16` → `rrn16`
- `rrn4.5s` → `rrn4.5`

**3. Protein-Coding Genes**
```python
GENE_SYNONYMS
    "pafI": "ycf3",
    "paf1": "ycf3",
    "pafII": "ycf4",
    "ycf10": "cemA",
    "lhbA": "psbZ",
    "pbf1": "psbN",
    "clpP1": "clpP",
    "infA": "infA",
    "matK": "matK",
    "accD": "accD"
}
```

**Examples:**
- `pafI` → `ycf3`
- `lhbA` → `psbZ`
- `ycf10` → `cemA`

### Normalization Tracking

All changes are logged in `Gene_Normalization_Log.xlsx`:
- Original name
- Normalized name
- Type of normalization
- Any warnings

This ensures transparency and reproducibility.

---

## Pseudogene Detection

### Comprehensive Detection Strategy

The module uses multiple criteria to identify pseudogenes:

**1. Explicit Annotation**
```gb
gene            /gene="ycf15"
                /pseudo
CDS             /gene="ycf15"
                /pseudo
```

**2. Implicit Markers**
```gb
gene            /gene="infA"
CDS             /gene="infA"
                /product="pseudogene"
```

**3. Missing CDS**
```gb
gene            /gene="accD"
# No corresponding CDS feature → pseudogene
```

**4. Internal Stop Codons**
```gb
CDS             /gene="ycf1"
                /transl_except="(pos:1234,aa:TERM)"
```

### Detection Logic

```python
def classify_gene_locus(features, gene_name):
    has_gene_feature = False
    has_cds = False
    has_functional_cds = False
    explicit_pseudo = False
    has_internal_stops = False
    
    for feature in features:
        # Track feature types
        if feature.type == "gene":
            has_gene_feature = True
        
        if feature.type == "CDS":
            has_cds = True
            # Check for functionality indicators
            
        # Check for explicit pseudo markers
        if "pseudo" in feature.qualifiers:
            explicit_pseudo = True
        
        # Check product/note for pseudo keywords
        for field in ["product", "note"]:
            if field in feature.qualifiers:
                text = " ".join(feature.qualifiers[field]).lower()
                if any(keyword in text for keyword in ["pseudo", "pseudogene", "non-functional"]):
                    explicit_pseudo = True
        
        # Check for internal stops
        if "transl_except" in feature.qualifiers:
            text = " ".join(feature.qualifiers["transl_except"]).lower()
            if "stop" in text or "ter" in text:
                has_internal_stops = True
    
    # Classification decision
    if explicit_pseudo or has_internal_stops:
        return "pseudogene"
    
    # For protein-coding genes: gene without CDS = pseudogene
    if not gene_name.lower().startswith(('trn', 'rrn')):
        if has_gene_feature and not has_cds:
            return "pseudogene"
    
    # Default to functional
    return "functional"
```

### Special Cases

**rps12 (trans-spliced):**
- Has two separate CDS features
- Considered functional despite split nature
- Not counted as duplicated

**tRNA/rRNA pseudogenes:**
- Must have explicit pseudo markers
- Missing product annotation suggests pseudogene

**Partial genes:**
- Marked as pseudogenes if truncated
- Check for "partial" in product/note

---

## IR Duplication Analysis

### Inverted Repeat Detection

The module first identifies IR regions:
```python
def detect_inverted_repeat_regions(record):
    ira_ranges = []
    irb_ranges = []
    
    for feature in record.features:
        if feature.type.lower() in ("repeat_region", "misc_feature"):
            annotation_text = get_annotation_text(feature)
            
            if "inverted repeat" in annotation_text or "ira" in annotation_text:
                span = (int(feature.location.start), int(feature.location.end))
                if "ira" in annotation_text:
                    ira_ranges.append(span)
                elif "irb" in annotation_text:
                    irb_ranges.append(span)
    
    # Merge overlapping ranges
    return _merge_genomic_ranges(ira_ranges), _merge_genomic_ranges(irb_ranges)
```

### Gene Location Analysis

For each gene, determine its location:
```python
for gene_name, features in gene_features.items():
    for feature in features:
        span = (int(feature.location.start), int(feature.location.end))
        
        in_ira = any(check_span_overlap(span, ir) for ir in ira_ranges)
        in_irb = any(check_span_overlap(span, ir) for ir in irb_ranges)
        
        if in_ira and in_irb:
            region = "IR(both)"  # Spans both IRs (rare)
        elif in_ira:
            region = "IRa"
        elif in_irb:
            region = "IRb"
        else:
            region = "single_copy"
```

### Duplication Characterization

**True duplications:**
- Gene present in both IRa and IRb
- Typically identical or nearly identical sequences
- Marked as "IRa+IRb"

**Single copies:**
- Gene only in IRa or only in IRb
- Gene in single-copy regions
- Marked accordingly

**Special cases:**
- Genes spanning both IR regions
- Genes with partial duplications
- Genes with pseudogene copies

### Duplicate Status Reporting

The module provides detailed duplication status:
- `single copy` - Not duplicated
- `2 functional copies (IRa+IRb)` - True duplication
- `1 functional + 1 pseudogene copies (IRa+IRb)` - Mixed
- `2 pseudogene copies (IRa+IRb)` - Both pseudogenes
- `trans-spliced` - Special case (rps12)

---

## Troubleshooting

### Common Issues and Solutions

#### 1. Missing Python Packages
```bash
❌ ERROR: Required package not installed: No module named 'Bio'
```
**Solution:**
```bash
pip install biopython pandas openpyxl

# Verify
python -c "import Bio, pandas, openpyxl"
```

#### 2. No GenBank Files Found
```bash
❌ ERROR: No GenBank files found!
```
**Solution:**
```bash
# Check file extensions
ls *.gb
ls *.gbf
ls *.gbk
ls *.genbank

# Make sure files are in the correct directory
python cgas_module5.py -i /full/path/to/genbank_files/
```

#### 3. No Gene Features Found
```bash
⚠ Warning: No genes found in genome.gb
```
**Solution:**
```
Check that your GenBank files contain gene features:
- Open file in a text editor
- Look for "gene" in the FEATURES section
- Verify proper /gene qualifiers
```

#### 4. Unusual Gene Names
```bash
⚠ Warning: Non-standard tRNA format: trnX
```
**Solution:**
```
This indicates unusual gene naming:
- Check the normalization log for details
- Verify annotation quality
- Some manual curation may be needed
```

#### 5. Many Pseudogenes Detected
```bash
⚠ Warning: High number of pseudogenes detected
```
**Solution:**
```
This may indicate:
- Poor annotation quality
- Ancient pseudogenization
- Overly strict detection criteria
- Review the normalization log for details
```

### Quick Fixes Table

| Problem | Symptom | Quick Fix |
|---------|---------|-----------|
| Missing packages | Import error | `pip install biopython pandas openpyxl` |
| No GenBank files | "No files found" | Check file extensions (.gb, .gbf, .gbk) |
| No gene features | "No genes found" | Verify GenBank files have gene features |
| Unusual names | Warning messages | Check normalization log |
| Many pseudogenes | High count | Review annotation quality |

---

## Examples

### Example 1: Complete Workflow
```bash
# 1. Prepare GenBank files
mkdir -p /home/abdullah/gene_comparison/
cd /home/abdullah/gene_comparison/

# 2. Run gene comparative analysis
python cgas_module5.py

# 3. Check results
ls Module5_Gene_Comparative_Analysis/
open Module5_Gene_Comparative_Analysis/Chloroplast_Gene_Analysis_*.xlsx
```

### Example 2: Process Specific Folder
```bash
# GenBank files in separate folder
python cgas_module5.py -i chloroplast_genomes/ -o gene_results/

# Check what was created
ls gene_results/
```

### Example 3: Analyzing Normalization
```bash
# After running analysis:
# Check the normalization log to see:
# - Which names were changed
# - Why they were changed
# - Any unusual patterns
open Module5_Gene_Comparative_Analysis/Gene_Normalization_Log.xlsx
```

### Example 4: Interpreting Results
```bash
# Key things to look for in the output:
# 1. Genes with different functional status across species
# 2. Species-specific genes
# 3. Differences in IR duplication patterns
# 4. Unusual pseudogene patterns
```

---

## FAQ

### Q1: What's the simplest way to run this?
**A:** Just `python cgas_module5.py` in a folder with GenBank files.

### Q2: Do my GenBank files need gene features?
**A:** Yes! Files must contain `gene` features with proper `/gene` qualifiers.

### Q3: What does gene normalization do?
**A:** It standardizes gene names across different annotation conventions (e.g., pafI → ycf3).

### Q4: How are pseudogenes detected?
**A:** Using multiple criteria: explicit markers, missing CDS, internal stops, etc.

### Q5: Can I customize the gene synonyms?
**A:** Yes, you can modify the `GENE_SYNONYMS` dictionary in the script.

### Q6: What if I disagree with a pseudogene classification?
**A:** Check the normalization log for the detection criteria used.

### Q7: Can I analyze mitochondrial genomes?
**A:** Yes! It works with any GenBank files with gene features, though optimized for chloroplasts.

### Q8: What file formats are supported?
**A:** GenBank files with extensions: .gb, .gbf, .gbk, .genbank

### Q9: How do I cite this tool?
**A:** See Citation section below.

### Q10: Why is rps12 special?
**A:** rps12 is trans-spliced with parts in different locations, not a true duplication. The annotator also mostly do inaccurate annotation. So, careful evaluation of annotation needed. 
### Q11: Does the genome contain inverted repeats regoin annotatoins? 
**A:** Yes, this is important. If you performed annotations using module 2 of CGAS then the genome already contain annotations for LSC, SSC, and IRs (IRa and IRb). However, If you downloaded from NCBI, then please check the annotatoins of these regions carefully. 

---

## Technical Specifications

### Performance
- **Processing speed**: ~5-10 seconds per genome
- **Memory usage**: <200 MB RAM
- **Disk space**: Minimal (<10 MB per analysis)

### Compatibility
- **Python**: 3.9+ (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating Systems**: Windows, macOS, Linux
- **Excel**: Generates .xlsx files (Excel 2010+)

### Input Limits
- **Max genome size**: No practical limit (tested with 300kb+ genomes)
- **Max files**: No practical limit (tested with 50+ files)
- **Gene features**: No limit (all gene features are processed)

### Quality Features
- ✅ Comprehensive gene detection
- ✅ Intelligent name normalization
- ✅ Multi-criteria pseudogene detection
- ✅ IR duplication analysis
- ✅ Comparative statistics
- ✅ Normalization tracking
- ✅ Publication-ready output

---

## References
### Python Libraries
- **Biopython**: [Official Documentation](https://biopython.org/) - GenBank file parsing
- **pandas**: [Documentation](https://pandas.pydata.org/) - Data manipulation
- **openpyxl**: [Documentation](https://openpyxl.readthedocs.io/) - Excel file creation

### GenBank Format
- **NCBI**: [GenBank Format Documentation](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/) - File format specifications

### CGAS Suite
- **CGAS**: Chloroplast Genome Analysis Suite - integrated tools for comprehensive chloroplast genomics

### Citation
If using CGAS Module 5 in publications, please cite:
```
Abdullah, Yan, R., & Tian, X. (2025). CGAS (Chloroplast Genome Analysis Suite): 
An Automated Python Pipeline for Comprehensive Comparative Chloroplast Genomics. 
bioRxiv. https://doi.org/10.64898/2025.12.21.695765
```

---

## Support & Contact

### Getting Help
```bash
# 1. First check built-in help
python cgas_module5.py --help

# 2. Verify packages are installed
python -c "import Bio, pandas, openpyxl"

# 3. Check GenBank files
ls *.gb *.gbf *.gbk *.genbank

# 4. Verify GenBank format
head -50 your_genbank.gb
```

### Common Issues Solved Here
- ✅ Missing packages? Run `pip install biopython pandas openpyxl`
- ✅ No files found? Check file extensions (.gb, .gbf, .gbk, .genbank)
- ✅ No gene features? Verify GenBank files have gene annotations
- ✅ Unusual names? Check normalization log
- ✅ Many pseudogenes? Review annotation quality

### Quick Reference Card
```bash
# ⚡ QUICK START ⚡
cgas --module 5                                  # cgas command
cgas-gene-compare                                # shortcut command
python cgas_module5.py                           # python command
python cgas_module5.py -i genbank_files/         # Specify input
python cgas_module5.py -i data/ -o results/      # Custom output

# 🔬 JUPYTER NOTEBOOK 🔬
%run cgas_module5.py                             # %run works for Module 5
!cgas --module 5                                 # ! also works

# 📊 OUTPUT 📊
# Module5_Gene_Comparative_Analysis/
# ├── Chloroplast_Gene_Analysis_YYYYMMDD_HHMMSS.xlsx  # Main analysis
# └── Gene_Normalization_Log.xlsx                    # Normalization log

# 🎯 TIPS 🎯
# - GenBank files must have gene features
# - Gene names are automatically normalized
# - Pseudogenes detected by multiple criteria
# - IR duplications identified and characterized
# - Check normalization log for all name changes
```

---

**Last Updated**: January 2026  
**Version**: 1.0.1  
**Author**: Abdullah  
**Part of**: CGAS (Chloroplast Genome Analysis Suite)

---

*Happy Comparative Gene Analysis! 🧬✨*