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

## Finding Your PGA Path

Before running the CGAS Module 2 script, you need to know the full path to your PGA installation.

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

**Copy the full path shown** — you'll need this for the `-pga` argument.

---

## Running the CGAS Module 2 Script

### Basic Usage
```bash
python cgas_module2.py -i <input_fasta> -o <output_directory> -pga <path_to_pga>
```

### Example Commands

**On WSL/Linux:**
```bash
python cgas_module2.py -i input.fasta -o output_results -pga "/home/abdullah/PGA/PGA.pl"
```

**If your username is different:**
```bash
# First, find your path:
echo "$HOME/PGA/PGA.pl"

# Then use that path:
python cgas_module2.py -i input.fasta -o output_results -pga "/home/YOUR_USERNAME/PGA/PGA.pl"
```

**On Windows (if running directly with Windows Python):**
```bash
python cgas_module2.py -i input.fasta -o output_results -pga "C:\Path\To\PGA_windows.exe"
```

### Script Arguments

| Flag | Description | Required |
|------|-------------|----------|
| `-i` / `--input` | Input FASTA file | Yes |
| `-o` / `--output` | Output directory | Yes |
| `-pga` / `--pga_path` | Full path to PGA executable or Perl script | Yes |

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
- The `-pga` path should point to either:
  - The Perl script: `/home/username/PGA/PGA.pl`
  - Or the executable if using a compiled version
- **Always use the full absolute path** for the `-pga` argument

---

## Troubleshooting

### Error: PGA path not found
**Solution:**
1. Verify your PGA path exists:
```bash
   ls /home/abdullah/PGA/PGA.pl
```
2. If file not found, use the find command:
```bash
   find ~ -name "PGA.pl" 2>/dev/null
```
3. Use the absolute path shown, not relative paths

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

### How to verify your path is correct
```bash
# Test if the path works:
perl /home/abdullah/PGA/PGA.pl -h

# If it shows help message, your path is correct!
```