# PGA (Plastome Genome Annotator) — Installation Guide on WSL with Anaconda

## Prerequisites
- Windows Subsystem for Linux (WSL)
- Anaconda installed
- A conda environment named `cgas`

---

## Step 1: Activate Your Conda Environment
```bash
conda activate cgas
```

---

## Step 2: Install BLAST+
```bash
conda install -c bioconda blast -y #No need If you are working in CGAS conda enviroment
```

---

## Step 3: Install Perl (if not already installed)
```bash
sudo apt update && sudo apt install -y perl 
```

---

## Step 4: Clone PGA
```bash
git clone https://github.com/quxiaojian/PGA.git
```

---

## Step 5: Create the Alias
```bash
echo 'alias pga="perl /home/abdullah/PGA/PGA.pl"' >> ~/.bashrc
source ~/.bashrc
```

> Replace `abdullah` with your username if different.

---

## Step 6: Test PGA
```bash
conda activate cgas
pga -h
```

If the help message prints, PGA is ready to use.

---

## Running the CGAS Module 2 Script

### ✨ Auto-Detection Feature

**The script now automatically detects your PGA installation!** If you followed the installation steps above, you can simply run:

```bash
python cgas_module2.py -i <input_directory> -r <reference_directory>
```

The script will automatically search for PGA in common locations:
- `~/PGA/PGA.pl`
- System PATH (using the `pga` alias)
- `/usr/local/bin/pga`
- And other standard locations

### Basic Usage

**Simple command (recommended - auto-detects PGA):**
```bash
python cgas_module2.py -i genomes/ -r refs/
```

**With organism names:**
```bash
python cgas_module2.py -i genomes/ -r refs/ --organism-file species.txt
```

**With explicit PGA path (optional, if auto-detection fails):**
```bash
python cgas_module2.py -i genomes/ -r refs/ --pga "/home/abdullah/PGA/PGA.pl"
```

**With output directory:**
```bash
python cgas_module2.py -i genomes/ -r refs/ -o my_annotations/
```

### Script Arguments

| Flag | Description | Required |
|------|-------------|----------|
| `-i` / `--input` | Input directory with FASTA files | Yes |
| `-r` / `--reference` | Directory with reference GenBank files | Yes |
| `-o` / `--output` | Output directory (default: module2_annotations) | No |
| `--organism-file` | TSV file mapping accession to organism name | No |
| `--pga` | Path to PGA (optional - auto-detects if not provided) | No |
| `--force-rerun` | Force re-annotation of all samples | No |
| `--log-level` | Logging level (DEBUG/INFO/WARNING/ERROR) | No |

---

## Organism File Format

The `--organism-file` should be a TSV (tab-separated) file with this format:

```
accession	organism
27	Hibiscus rosa-sinensis
SRR123	Arabidopsis thaliana
```

**For SSC flip-flop variants:** If you have files like `27_1.fasta` and `27_2.fasta` from the same species, you only need one entry:

```
accession	organism
27	Hibiscus rosa-sinensis
```

The script will intelligently apply "Hibiscus rosa-sinensis" to both `27_1` and `27_2`.

---

## Finding Your PGA Path (Manual Method)

If you need to provide the PGA path manually (using `--pga`), here's how to find it:

### Method 1: Check where you cloned PGA
```bash
cd ~/PGA
pwd
```

This will show something like: `/home/abdullah/PGA`

Your full PGA path will be: `/home/abdullah/PGA/PGA.pl`

### Method 2: Use the `which` command (if alias is set up)
```bash
which pga
```

This might show: `alias pga='perl /home/abdullah/PGA/PGA.pl'`

The path you need is: `/home/abdullah/PGA/PGA.pl`

### Method 3: Find command
```bash
find ~ -name "PGA.pl" 2>/dev/null
```

This will search your home directory and show the exact path.

### Method 4: Check your current username

If you cloned PGA in your home directory:
```bash
echo "$HOME/PGA/PGA.pl"
```

**Copy the full path shown** — you can use it with the `--pga` argument if needed.

---

## Example Workflows

### Single FASTA file
```bash
python cgas_module2.py -i targets/sample.fasta -r refs/
```

### Directory of FASTA files
```bash
python cgas_module2.py -i targets/ -r refs/ -o annotations/
```

### With organism names and custom output
```bash
python cgas_module2.py -i targets/ -r refs/ -o my_output/ --organism-file species.txt
```

### Force re-annotation (skip existing files)
```bash
python cgas_module2.py -i targets/ -r refs/ --force-rerun
```

### Debug mode (detailed logging)
```bash
python cgas_module2.py -i targets/ -r refs/ --log-level DEBUG
```

---

## PGA Direct Usage (Optional)

You can also use PGA directly without the wrapper script:
```bash
pga -r <reference_folder> -t <target_folder>
```

| Flag | Description |
|------|-------------|
| `-r` | Input directory with GenBank-formatted reference plastome(s) |
| `-t` | Input directory with FASTA-formatted target plastome(s) |
| `-i` | Minimum IR length (default: 1000) |
| `-p` | Minimum TBLASTN percent identity (default: 40) |
| `-q` | Query coverage thresholds (default: 0.5, 2) |
| `-o` | Output directory name (default: gb) |
| `-f` | Plastome form — circular or linear (default: circular) |
| `-l` | Log file name (default: warning) |

---

## Notes
- Always activate `cgas` before running PGA or the script.
- Replace `anaconda3` with `miniconda3` if you use Miniconda.
- The script auto-detects PGA in most cases - manual path only needed if auto-detection fails
- **Always use the full absolute path** if providing `--pga` argument manually

---

## Troubleshooting

### Error: PGA not found / Could not find PGA installation
**Solution:**

The script searches multiple locations automatically. If it still can't find PGA:

1. Verify your PGA installation:
```bash
ls ~/PGA/PGA.pl
```

2. If PGA exists, provide the path manually:
```bash
python cgas_module2.py -i targets/ -r refs/ --pga "$HOME/PGA/PGA.pl"
```

3. If PGA doesn't exist, reinstall following Steps 4-6 above.

### Error: Permission denied
**Solution:**
```bash
chmod +x cgas_module2.py
```

### Error: Command not found
**Solution:**
- Make sure you've activated the conda environment:
```bash
conda activate cgas
```

### How to verify your PGA installation is correct
```bash
# Test if PGA works:
pga -h
# OR
perl ~/PGA/PGA.pl -h

# If it shows help message, PGA is installed correctly!
```

### Script finds PGA but annotation fails
**Solution:**
1. Check that reference directory contains `.gb` or `.gbk` files:
```bash
ls refs/*.gb
```

2. Check that input directory contains `.fasta` or `.fa` files:
```bash
ls targets/*.fasta
```

3. Run with debug logging to see details:
```bash
python cgas_module2.py -i targets/ -r refs/ --log-level DEBUG
```

---

## Output Structure

After successful annotation, you'll find:

```
module2_annotations/
├── Annotated_GenBank/          # Annotated .gb files
│   ├── Hibiscus_rosa-sinensis_27_1.gb
│   └── Hibiscus_rosa-sinensis_27_2.gb
├── Annotation_Logs/            # PGA warning logs
│   ├── Hibiscus_rosa-sinensis_27_1_warning.log
│   └── Hibiscus_rosa-sinensis_27_2_warning.log
└── Reports/                    # Summary reports
    └── 00_ANNOTATION_SUMMARY.tsv
```
