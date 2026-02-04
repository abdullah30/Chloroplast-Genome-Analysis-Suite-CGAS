# CGAS Module 6: Publication-Quality Gene Content Tables
## Complete Documentation and User Guide

---

## 📋 Table of Contents
1. [Introduction](#introduction)
2. [Quick Start - Just Run It!](#quick-start---just-run-it)
3. [Command-Line Usage with Examples](#command-line-usage-with-examples)
4. [Jupyter Notebook Usage](#jupyter-notebook-usage)
5. [Input Requirements](#input-requirements)
6. [Output Structure](#output-structure)
7. [Gene Categorization System](#gene-categorization-system)
8. [Special Features Handling](#special-features-handling)
9. [Word Document Formatting](#word-document-formatting)
10. [Troubleshooting](#troubleshooting)
11. [Examples](#examples)
12. [FAQ](#faq)
13. [Technical Specifications](#technical-specifications)
14. [References](#references)

---

## Introduction

**CGAS Module 6** is a specialized tool for generating publication-quality gene content tables from chloroplast genome annotations. This module categorizes genes by function, detects special features like introns and pseudogenes, and creates professionally formatted Microsoft Word documents suitable for direct inclusion in scientific publications.

This module performs comprehensive gene content analysis with:
- **Gene categorization**: Organizes genes by functional groups
- **Feature detection**: Identifies introns, duplications, and pseudogenes
- **Duplicate handling**: Marks genes in inverted repeat regions
- **Publication formatting**: Creates Word documents with proper formatting
- **Combined output**: Merges all species tables into one document
- **Special markers**: Uses superscripts and symbols for special features

### Key Features:
- **Automatic categorization**: Groups genes into Self-replication, Photosynthesis, and Other
- **Intron detection**: Identifies genes with 1 or 2 introns
- **Pseudogene identification**: Detects and marks pseudogenes with Ψ symbol
- **Duplicate marking**: Uses superscript 'a' for IR-duplicated genes
- **Professional formatting**: Italicized gene names, proper table styling
- **Batch processing**: Handles multiple genomes simultaneously

### Scientific Applications:
- **Genome annotation**: Standardized gene content reporting
- **Comparative genomics**: Compare gene content across species
- **Phylogenetic studies**: Use gene content as phylogenetic characters
- **Publication preparation**: Generate tables directly for manuscripts
- **Taxonomic studies**: Characterize gene families in different lineages
- **Evolutionary analysis**: Track gene loss/gain patterns

### Part of CGAS (Chloroplast Genome Analysis Suite)
CGAS is a comprehensive suite of tools designed for chloroplast genome analysis, providing seamless integration across different analysis modules.

---

## Quick Start - Just Run It!

### ⚡ One Command Does It All!

```bash
# Option 1: Using cgas command
cgas --module 6

# Option 2: Using cgas-gene-table shortcut
cgas-gene-table

# Option 3: Using python directly
python cgas_module6.py
```

**What happens when you run this:**
1. ✅ Finds all GenBank files in current directory
2. ✅ Extracts and categorizes all genes
3. ✅ Detects introns, duplications, and pseudogenes
4. ✅ Creates individual Word documents per species
5. ✅ Generates combined master document
6. ✅ Applies publication-quality formatting

### 📁 Real-World Folder Examples

#### Example 1: Basic Usage (Easiest!)
```bash
# Your folder structure:
/home/abdullah/gene_content/
├── chloroplast_genomes/
│   ├── Arabidopsis_thaliana.gb
│   ├── Oryza_sativa.gb
│   └── Zea_mays.gb

# Navigate to your folder
cd /home/abdullah/gene_content/chloroplast_genomes/

# Run the module (no arguments needed!)
python cgas_module6.py

# Output created automatically:
# Module6_Gene_Content_Tables/
# ├── Table_Arabidopsis_thaliana.docx
# ├── Table_Oryza_sativa.docx
# ├── Table_Zea_mays.docx
# └── Complete_Gene_Content_Tables.docx
```

#### Example 2: Specify Input Folder
```bash
# Process genomes from specific directory
python cgas_module6.py -i genbank_files/

# Output created in:
# genbank_files/Module6_Gene_Content_Tables/
```

#### Example 3: Custom Input and Output
```bash
# Full control over directories
python cgas_module6.py -i chloroplast_genomes/ -o Gene_Tables/

# Input from: chloroplast_genomes/
# Output to: Gene_Tables/
```

### 🎯 Quick Start Cheat Sheet

```bash
# ====== cgas command ======
cgas --module 6                                  # Process current directory
cgas --module 6 -i genbank_files/                # Specify input folder
cgas --module 6 -i data/ -o results/             # Custom input and output

# ====== cgas-gene-table shortcut ======
cgas-gene-table                                  # Process current directory
cgas-gene-table -i genbank_files/                # Specify input folder
cgas-gene-table -i data/ -o results/             # Custom input and output

# ====== python command ======
python cgas_module6.py                           # Process current directory
python cgas_module6.py -i genbank_files/         # Specify input folder
python cgas_module6.py -i data/ -o results/      # Custom input and output
python cgas_module6.py --help                    # Get help
```

### 📊 What You Get (Output Files)

```
Module6_Gene_Content_Tables/                 # Created automatically
├── 📄 Individual Tables (one per species)
│   ├── Table_Arabidopsis_thaliana.docx     # Gene content table
│   ├── Table_Oryza_sativa.docx             # Gene content table
│   └── Table_Zea_mays.docx                 # Gene content table
│
└── 📑 Combined Document
    └── Complete_Gene_Content_Tables.docx   # All species combined
```

---

## Command-Line Usage with Examples

### Complete Command Reference

```bash
# ====================================================================
# USING cgas COMMAND
# ====================================================================

# Process current directory
cgas --module 6

# Specify input directory
cgas --module 6 -i /home/abdullah/chloroplast_genomes/

# Custom input and output directories
cgas --module 6 -i chloroplast_genomes/ -o Gene_Tables/
```

```bash
# ====================================================================
# USING cgas-gene-table SHORTCUT
# ====================================================================

# Process current directory
cgas-gene-table

# With specific input and output
cgas-gene-table -i chloroplast_genomes/ -o Gene_Tables/
```

```bash
# ====================================================================
# USING python COMMAND
# ====================================================================

# ====================================================================
# BASIC USAGE
# ====================================================================

# 1. Process current directory (simplest!)
python cgas_module6.py

# 2. Specify input folder
python cgas_module6.py -i genbank_files/

# 3. Custom input and output folders
python cgas_module6.py -i chloroplast_genomes/ -o Gene_Tables/

# 4. Get help
python cgas_module6.py --help


# ====================================================================
# REAL-WORLD EXAMPLES
# ====================================================================

# Example 1: Process chloroplast genomes
python cgas_module6.py -i chloroplast_genomes/

# Example 2: Save to specific output folder
python cgas_module6.py -i data/ -o ../results/gene_tables/

# Example 3: Windows users
python cgas_module6.py -i "C:\Users\abdullah\genbank_files"
```

### Parameter Details

| Parameter | Short | Default | Description |
|-----------|-------|---------|-------------|
| `--input` | `-i` | `.` (current) | Input directory with GenBank files |
| `--output` | `-o` | `Module6_Gene_Content_Tables` | Output directory for results |

---

## Jupyter Notebook Usage

> **Note:** Module 6 is a pure Python script, so both `%run` and the `!` shell operator work. Use `%run` for direct script execution, or `!` for the `cgas` / `cgas-gene-table` commands. Run the cell from the directory containing your GenBank files, or pass `-i` to specify the input path.

```python
# Using %run (executes the script directly)
%run cgas_module6.py
%run cgas_module6.py -i chloroplast_genomes/
%run cgas_module6.py -i data/ -o results/

# Using ! operator with cgas command
!cgas --module 6
!cgas --module 6 -i chloroplast_genomes/

# Using ! operator with cgas-gene-table shortcut
!cgas-gene-table
!cgas-gene-table -i chloroplast_genomes/ -o Gene_Tables/

# Using ! operator with python
!python cgas_module6.py
!python cgas_module6.py -i chloroplast_genomes/ -o Gene_Tables/
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
     tRNA            100..200
                     /gene="trnI"
                     /product="tRNA-Ile (CAU)"
     rRNA            300..1500
                     /gene="rrn16"
                     /product="16S ribosomal RNA"
ORIGIN
        1 atggcgacga cgttcgtcgt cgtttgtcga tctcgtctga cttcagcctg atcggtagca
...
//
```

**2. Gene Features Required**
- Files must contain `gene` features
- Each gene should have proper `/gene` qualifier
- Optional: `CDS`, `tRNA`, `rRNA` features for classification

**3. Standard Gene Names**
- Use standard chloroplast gene nomenclature
- Common gene names: rbcL, matK, atpB, psbA, trnL, rrn16, etc.
- Gene names should be consistent across files

### File Organization

```bash
genbank_files/
├── Arabidopsis_thaliana.gb          # Chloroplast genome 1
├── Oryza_sativa.gb                  # Chloroplast genome 2
├── Zea_mays.gb                      # Chloroplast genome 3
└── Nicotiana_tabacum.gb             # Chloroplast genome 4
```

### Example Gene Features

```gb
     gene            864..2690
                     /gene="rbcL"
     CDS             864..2690
                     /gene="rbcL"
                     /product="ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"
                     /translation="MASSSNSS..."

     gene            5000..5150
                     /gene="trnL"
                     /note="contains one intron"
     tRNA            5000..5150
                     /gene="trnL"
                     /product="tRNA-Leu (UAA)"

     gene            6000..6200
                     /gene="ycf1"
                     /pseudo
     CDS             6000..6200
                     /gene="ycf1"
                     /product="hypothetical protein"
                     /pseudo
```

**What gets analyzed:**
- All gene features are extracted and categorized
- Introns detected from multi-part CDS/tRNA features
- Duplications identified from gene counts
- Pseudogenes detected from `/pseudo` qualifier or absence of CDS

---

## Output Structure

### Directory Organization

```
Module6_Gene_Content_Tables/                 # Main output folder
├── 📄 Individual Tables (one per species)
│   ├── Table_Arabidopsis_thaliana.docx     # Individual table
│   │   └── Contains:
│   │       - Gene content table with categories
│   │       - Italicized gene names
│   │       - Special markers for introns/duplicates/pseudogenes
│   │       - Species-specific title
│   │
│   ├── Table_Oryza_sativa.docx             # Individual table
│   └── Table_Zea_mays.docx                 # Individual table
│
└── 📑 Combined Document
    └── Complete_Gene_Content_Tables.docx   # All species in one document
        └── Contains:
            - All individual tables
            - Sequential numbering (Table S1, S2, etc.)
            - Page breaks between species
```

### Key Output Files Explained

#### 1. Individual Tables (.docx)

**Table structure:**
```
Table S1. Gene content of the chloroplast genome of Arabidopsis thaliana

┌─────────────────────┬──────────────────────────┬─────────────────────────┬───────┐
│ Category for genes  │ Group of genes           │ Name of genes           │ Total │
├─────────────────────┼──────────────────────────┼─────────────────────────┼───────┤
│ Self-replication    │ Large subunit of ribosome│ rpl2, rpl16, rpl20      │ 3     │
│                     │ Small subunit of ribosome│ rps2, rps3, rps4        │ 3     │
│                     │ rRNA genes               │ rrn4.5, rrn5, rrn16     │ 4     │
│                     │ tRNA genes               │ trnA-UGC, trnC-GCA...   │ 30    │
│ Photosynthesis      │ Photosystem I            │ psaA, psaB, psaC        │ 3     │
│                     │ Photosystem II           │ psbA, psbB, psbC        │ 12    │
│ Other genes         │ Conserved ORFs           │ ycf1, ycf2, ycf4        │ 3     │
│                     │ Pseudogenes              │ ycf15Ψ                  │ 1     │
│                     │                          │                         │       │
│                     │                          │ Total number of genes   │ 79    │
└─────────────────────┴──────────────────────────┴─────────────────────────┴───────┘

Note: * and ** indicate genes containing one and two introns, respectively; 
a and Ψ indicate duplicated genes in inverted repeat (IR) regions and 
pseudogenes, respectively. The rps12 gene is a trans-spliced gene and is 
not marked as duplicated despite appearing in multiple locations.
```

**Special markers:**
- `*` = Gene with one intron
- `**` = Gene with two introns
- `^a` = Gene duplicated in IR regions (superscript)
- `Ψ` = Pseudogene (psi symbol)

#### 2. Complete_Gene_Content_Tables.docx

**Combined document features:**
- All species tables in one document
- Sequential table numbering (Table S1, S2, S3...)
- Page breaks between species
- Consistent formatting across all tables
- Species names italicized in titles

---


---

## Gene Categorization System

### Main Categories

**1. Self-replication**
- **Large subunit of ribosome**: `rpl2`, `rpl16`, `rpl20`, etc.
- **Small subunit of ribosome**: `rps2`, `rps3`, `rps4`, etc.
- **DNA dependent RNA polymerase**: `rpoA`, `rpoB`, `rpoC1`, `rpoC2`
- **rRNA genes**: `rrn4.5`, `rrn5`, `rrn16`, `rrn23`
- **tRNA genes**: `trnA-UGC`, `trnC-GCA`, `trnD-GUC`, etc.

**2. Photosynthesis**
- **Photosystem I**: `psaA`, `psaB`, `psaC`, `psaI`, `psaJ`
- **Photosystem II**: `psbA`, `psbB`, `psbC`, `psbD`, `psbE`, etc.
- **NADPH dehydrogenase**: `ndhA`, `ndhB`, `ndhC`, etc.
- **Cytochrome b/f complex**: `petA`, `petB`, `petD`, `petG`, `petL`, `petN`
- **Subunits of ATP synthase**: `atpA`, `atpB`, `atpE`, `atpF`, `atpH`, `atpI`
- **Large subunit of Rubisco**: `rbcL`
- **Photosynthesis assembly genes**: `ycf3` (pafI), `ycf4` (pafII)

**3. Other genes**
- **Protease**: `clpP`
- **Maturase**: `matK`
- **Envelop membrane protein**: `cemA`
- **Subunit of Acetyl-CoA-carboxylase**: `accD`
- **C-type cytochrome synthesis gene**: `ccsA`
- **Translation initiation factor**: `infA`
- **Conserved open reading frames**: `ycf1`, `ycf2`, `ycf15`, etc.
- **Pseudogenes**: All identified pseudogenes

**4. Excluding genes**
- Genes that don't fit into any category
- Unusual or lineage-specific genes
- Genes with non-standard names

### Gene Name Mapping

**Alternative names handled:**
```python
name_map = {
    "pafI": "ycf3 (pafI)",
    "pafII": "ycf4 (pafII)", 
    "lhbA": "psbZ (lhbA)",
    "pbf1": "psbN (pbf1)",
    "ycf3": "ycf3 (pafI)",
    "ycf4": "ycf4 (pafII)",
    "psbZ": "psbZ (lhbA)",
    "psbN": "psbN (pbf1)"
}
```

This ensures consistent naming across different annotation styles.

---

## Special Features Handling

### 1. Intron-Containing Genes

**Detection method:**
- Multi-part CDS or tRNA features indicate introns
- Number of parts - 1 = number of introns
- Common intron-containing genes: `ycf3`, `clpP`, `rps12`, `trnI`, `trnK`, etc.

**Marking system:**
- `*` = One intron (2 exons)
- `**` = Two introns (3 exons)
- Example: `ycf3**`, `clpP**`, `trnK*`

### 2. Duplicated Genes (IR Regions)

**Detection method:**
- Genes with count > 1 are considered duplicated
- Typically located in inverted repeat regions
- Special case: `rps12` (trans-spliced, not marked as duplicate)

**Marking system:**
- Superscript `a` after gene name
- Example: `rpl2^a`, `rrn16^a`, `ycf2^a`

### 3. Pseudogenes

**Detection criteria:**
- Explicit `/pseudo` qualifier
- "pseudo" or "pseudogene" in product/note
- Gene feature without corresponding CDS (protein-coding genes only)

**Marking system:**
- Psi symbol (Ψ) after gene name
- Example: `ycf15Ψ`, `infAΨ`

### 4. Combined Markers

**Multiple features can be combined:**
- `ycf3**^a` = Two introns + duplicated
- `trnI*^aΨ` = One intron + duplicated + pseudogene
- Order of markers: introns first, then duplication, then pseudogene

### 5. Special Cases

**rps12 gene:**
- Trans-spliced gene (parts in different locations)
- Never marked as duplicate despite appearing multiple times
- Important exception to duplication rule

**ycf3 and ycf4:**
- Placed in "Photosynthesis assembly genes" group
- Alternative names (pafI, pafII) handled correctly

---

## Word Document Formatting

### Professional Table Styling

**Table structure:**
- 4 columns: Category, Group, Genes, Total
- Merged cells for categories
- Centered alignment for most columns
- Left alignment for gene names

**Text formatting:**
- **Gene names**: Italicized
- **Special markers**: Superscripts for duplicates
- **Headers**: Bold
- **Species names**: Italicized in titles

### Document Layout

**Individual tables:**
- Species-specific title: "Table S1. Gene content of..."
- Compact layout optimized for publication
- Footnote explaining markers

**Combined document:**
- Sequential numbering: Table S1, S2, S3...
- Page breaks between species
- Consistent formatting throughout

### Advanced Formatting Features

**Superscript handling:**
```python
def add_superscript(run, text):
    run.font.superscript = True
    run.text = text
```

**Italic gene names:**
```python
run = p.add_run(base_name)
run.italic = True
```

**Cell merging:**
- Category cells merged across multiple rows
- Vertical centering for merged cells
- Professional appearance

### Publication-Ready Output

**Features for direct manuscript inclusion:**
- Standard scientific table format
- Clear legend/footnote
- Consistent gene nomenclature
- Proper special character handling

**Compatibility:**
- Microsoft Word 2010+
- LibreOffice (with some formatting differences)
- Google Docs (limited formatting support)

---

## Troubleshooting

### Common Issues and Solutions

#### 1. Missing Python Packages
```bash
❌ ERROR: Required package not installed: No module named 'docx'
```
**Solution:**
```bash
pip install biopython pandas python-docx

# Verify
python -c "import Bio, pandas, docx"
```

#### 2. No GenBank Files Found
```bash
❌ No GenBank files found!
```
**Solution:**
```bash
# Check file extensions
ls *.gb
ls *.gbf
ls *.gbk
ls *.genbank

# Make sure files are in the correct directory
python cgas_module6.py -i /full/path/to/genbank_files/
```

#### 3. No Gene Features Found
```bash
⚠ Warning: Empty gene content table for species.gb
```
**Solution:**
```
Check that your GenBank files contain gene features:
- Open file in a text editor
- Look for "gene" in the FEATURES section
- Verify proper /gene qualifiers
```

#### 4. Word Document Issues
```bash
⚠ Warning: Error creating Word document
```
**Solution:**
```
This may happen if:
- python-docx version is too old
- File permissions are restrictive
- Disk space is limited

Try updating: pip install --upgrade python-docx
```

#### 5. Special Character Display
```bash
⚠ Warning: Psi symbol (Ψ) not displaying correctly
```
**Solution:**
```
This is a font/encoding issue:
- Open in Microsoft Word for best results
- Ensure Unicode font is used
- Some PDF viewers may not render correctly
```

### Quick Fixes Table

| Problem | Symptom | Quick Fix |
|---------|---------|-----------|
| Missing packages | Import error | `pip install biopython pandas python-docx` |
| No GenBank files | "No files found" | Check file extensions (.gb, .gbf, .gbk) |
| No gene features | Empty table | Verify GenBank files have gene features |
| Word doc error | Document creation fails | Update python-docx |
| Display issues | Special characters wrong | Use Microsoft Word |

---

## Examples

### Example 1: Complete Workflow
```bash
# 1. Prepare GenBank files
mkdir -p /home/abdullah/gene_content/
cd /home/abdullah/gene_content/

# 2. Run gene content analysis
python cgas_module6.py

# 3. Check results
ls Module6_Gene_Content_Tables/
open Module6_Gene_Content_Tables/Complete_Gene_Content_Tables.docx
```

### Example 2: Process Specific Folder
```bash
# GenBank files in separate folder
python cgas_module6.py -i chloroplast_genomes/ -o gene_tables/

# Check what was created
ls gene_tables/
```

### Example 3: Batch Processing with Script
```bash
#!/bin/bash
# Process multiple datasets

for dataset in monocots dicots ferns; do
    echo "Analyzing $dataset..."
    python cgas_module6.py \
        -i ${dataset}_genomes/ \
        -o ${dataset}_gene_tables/
    echo "Completed $dataset"
done

echo "All datasets analyzed!"
```

### Example 4: Interpreting Special Markers
```bash
# After running analysis, you might see:
# - ycf3** = ycf3 with two introns
# - rpl2^a = rpl2 duplicated in IR regions
# - ycf15Ψ = ycf15 is a pseudogene
# - rps12 = trans-spliced (not marked as duplicate)
```

---

## FAQ

### Q1: What's the simplest way to run this?
**A:** Just `python cgas_module6.py` in a folder with GenBank files.

### Q2: Do my GenBank files need gene features?
**A:** Yes! Files must contain `gene` features with proper `/gene` qualifiers.

### Q3: How are introns detected?
**A:** From multi-part CDS or tRNA features. Number of parts - 1 = number of introns.

### Q4: What if a gene is duplicated but not in IRs?
**A:** It will still be marked with superscript 'a'. The module doesn't verify IR location.

### Q5: Can I customize the gene categories?
**A:** Yes, you can modify the categories dictionary in the script.

### Q6: Why is rps12 not marked as duplicate?
**A:** rps12 is a trans-spliced gene with parts in different locations, not a true duplication.

### Q7: Can I analyze mitochondrial genomes?
**A:** Yes! It works with any GenBank files with gene features, though categories are optimized for chloroplasts.

### Q8: What file formats are supported?
**A:** GenBank files with extensions: .gb, .gbf, .gbk, .genbank

### Q9: How do I cite this tool?
**A:** See Citation section below.

### Q10: Can I edit the Word documents?
**A:** Yes! The documents are fully editable in Microsoft Word.

---

## Technical Specifications

### Performance
- **Processing speed**: ~5-10 seconds per genome
- **Memory usage**: <100 MB RAM
- **Disk space**: Minimal (<1 MB per document)

### Compatibility
- **Python**: 3.9+ (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating Systems**: Windows, macOS, Linux
- **Word**: Microsoft Word 2010+ recommended

### Input Limits
- **Max genome size**: No practical limit (tested with 300kb+ genomes)
- **Max files**: No practical limit (tested with 50+ files)
- **Gene features**: No limit (all gene features are processed)

### Quality Features
- ✅ Automatic gene categorization
- ✅ Intron detection and marking
- ✅ Pseudogene identification
- ✅ Duplicate gene marking
- ✅ Publication-quality formatting
- ✅ Combined document generation
- ✅ Special character handling

---

## References

### Python Libraries
- **Biopython**: [Official Documentation](https://biopython.org/) - GenBank file parsing
- **pandas**: [Documentation](https://pandas.pydata.org/) - Data manipulation
- **python-docx**: [Documentation](https://python-docx.readthedocs.io/) - Word document creation

### CGAS Suite
- **CGAS**: Chloroplast Genome Analysis Suite - integrated tools for comprehensive chloroplast genomics

### Citation
If using CGAS Module 6 in publications, please cite:
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
python cgas_module6.py --help

# 2. Verify packages are installed
python -c "import Bio, pandas, docx"

# 3. Check GenBank files
ls *.gb *.gbf *.gbk *.genbank

# 4. Verify GenBank format
head -50 your_genbank.gb
```

### Common Issues Solved Here
- ✅ Missing packages? Run `pip install biopython pandas python-docx`
- ✅ No files found? Check file extensions (.gb, .gbf, .gbk, .genbank)
- ✅ No gene features? Verify GenBank files have gene annotations
- ✅ Word doc issues? Update python-docx package
- ✅ Display problems? Use Microsoft Word for best results

### Quick Reference Card
```bash
# ⚡ QUICK START ⚡
cgas --module 6                                  # cgas command
cgas-gene-table                                  # shortcut command
python cgas_module6.py                           # python command
python cgas_module6.py -i genbank_files/         # Specify input
python cgas_module6.py -i data/ -o results/      # Custom output

# 🔬 JUPYTER NOTEBOOK 🔬
%run cgas_module6.py                             # %run works for Module 6
!cgas --module 6                                 # ! also works

# 📊 OUTPUT 📊
# Module6_Gene_Content_Tables/
# ├── Table_[species].docx              # Individual tables
# └── Complete_Gene_Content_Tables.docx # Combined document

# 🎯 TIPS 🎯
# - GenBank files must have gene features
# - * = one intron, ** = two introns
# - ^a = duplicated in IR regions
# - Ψ = pseudogene
# - rps12 is trans-spliced (not marked as duplicate)
```

---

**Last Updated**: January 2026  
**Version**: 1.0.1  
**Author**: Abdullah  
**Part of**: CGAS (Chloroplast Genome Analysis Suite)

---

*Happy Gene Content Table Generation! 🧬✨*