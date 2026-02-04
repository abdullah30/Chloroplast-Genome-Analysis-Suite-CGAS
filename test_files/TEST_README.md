## Test Data Usage by Module

The following datasets are provided for testing and validation of different CGAS modules:

### 1. GenBank files (`gb_files`)
Can be used to test **Modules 4–14**, excluding Module 10.

### 2. FASTA alignments (`fasta_files`)
Used specifically for **Module 10**, which requires aligned FASTA sequences.

### 3. GenBank files for annotation validation (`annotations`)
Can be used to test genome annotations generated via **PGA**.

### 4. GenBank files for normalization testing (`normalization`)
Can be used to verify the **normalization step in Module 3**.  
**Note:** Please use `Erigeron_acris.gb` as reference.

### 5. SRA datasets
The following SRA datasets can be used as input data for **Module 1**:
- `SRR8666554`
- `SRR8666754`
- `SRR8591790`

#### Download Links

**SRR8666554:**
```bash
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR866/004/SRR8666554/SRR8666554_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR866/004/SRR8666554/SRR8666554_2.fastq.gz
```

**SRR8666754:**
```bash
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR866/004/SRR8666754/SRR8666754_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR866/004/SRR8666754/SRR8666754_2.fastq.gz
```

**SRR8591790:**
```bash
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR859/000/SRR8591790/SRR8591790_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR859/000/SRR8591790/SRR8591790_2.fastq.gz
```

---

**Quick Download Script:**
```bash
# Download all SRA datasets at once
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR866/004/SRR8666554/SRR8666554_{1,2}.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR866/004/SRR8666754/SRR8666754_{1,2}.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR859/000/SRR8591790/SRR8591790_{1,2}.fastq.gz
```
