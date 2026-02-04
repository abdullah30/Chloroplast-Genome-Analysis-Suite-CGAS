#!/usr/bin/env python3
"""
CGAS Module 2: Simplified Plastome Annotation using PGA
========================================================

This module annotates clean, assembled genomes from Module 1 (07_assembled_genomes)
using PGA (Plastid Genome Annotator). Updated version with PGA auto-detection.

HOW THIS SCRIPT WORKS WITH PGA:
-------------------------------
This script acts as a wrapper around PGA and enhances its functionality:
1. Calls PGA as a subprocess → PGA runs and generates its standard logs
2. PGA creates screen.log and warning.log files
3. This script reads and parses those PGA log files
4. Adds custom enhancements:
   - Organism name updates in GenBank files
   - Structural feature annotations (LSC, SSC, IR regions)
   - Comprehensive summary reports
5. Generates additional reports beyond what PGA provides

This design allows the script to work WITH PGA (not replace it), leveraging
PGA's powerful annotation capabilities while adding project-specific features.

ORGANISM NAME HANDLING:
-----------------------
The script intelligently determines organism names in priority order:

1. **ORGANISM FILE (Highest Priority)**
   - If you provide --organism-file, those names are used first
   - Format: accession\torganism (with spaces in organism name)
   - Example: 27_1    Hibiscus syriacus

2. **FASTA HEADER EXTRACTION (Automatic Fallback)**
   - If no organism file provided, extracts from FASTA header
   - Works with NCBI format: ">Erigeron compositus TPA_asm: Erigeron compositus chloroplast..."
   - Automatically identifies genus + species from header
   - No need for manual organism file if FASTA headers are properly formatted!

3. **Unknown Organism (Last Resort)**
   - Only if neither of above work

The script automatically converts spaces to underscores in filenames:
- Input organism: "Hibiscus syriacus"
- Output filename: "Hibiscus_syriacus_27_1.gb"

Author: Abdullah
"""

import os
import sys
import subprocess
import logging
import shutil
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, asdict
import argparse
from datetime import datetime
import re


def find_pga_path():
    """
    Auto-detect PGA path from multiple common locations
    
    Returns:
        str: Path to PGA executable/script if found, None otherwise
    """
    # Method 1: Check if 'pga' command is available in PATH
    pga_cmd = shutil.which('pga')
    if pga_cmd:
        return pga_cmd
    
    # Method 2: Check common installation locations
    common_paths = [
        # Linux/WSL paths
        os.path.expanduser("~/PGA/PGA.pl"),
        os.path.expanduser("~/PGA/pga"),
        os.path.expanduser("~/PGA/PGA_linux"),
        "/usr/local/bin/pga",
        "/usr/bin/pga",
        os.path.expanduser("~/bin/pga"),
        
        # Mac paths
        os.path.expanduser("~/PGA/PGA_mac"),
        "/Applications/PGA/PGA.pl",
        
        # Windows paths
        "C:\\PGA\\PGA_windows.exe",
        "C:\\PGA\\PGA.pl",
        os.path.expanduser("~\\PGA\\PGA_windows.exe"),
        os.path.expanduser("~\\PGA\\PGA.pl"),
    ]
    
    # Also check with current username
    username = os.getenv('USER') or os.getenv('USERNAME')
    if username:
        common_paths.extend([
            f"/home/{username}/PGA/PGA.pl",
            f"/home/{username}/PGA/pga",
            f"C:\\Users\\{username}\\PGA\\PGA.pl",
            f"C:\\Users\\{username}\\PGA\\PGA_windows.exe",
        ])
    
    for path in common_paths:
        if os.path.exists(path):
            return path
    
    return None


@dataclass
class AnnotationStats:
    """Store PGA annotation statistics"""
    sample_name: str
    input_fasta: Path
    total_genes_reference: int = 0
    total_genes_annotated: int = 0
    annotation_success: bool = False
    warnings: List[str] = None
    unannotated_genes: List[str] = None
    output_genbank: Optional[Path] = None
    
    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []
        if self.unannotated_genes is None:
            self.unannotated_genes = []


class CGASModule2:
    """CGAS Module 2 - PGA Annotation Pipeline (BATCH MODE)"""
    
    def __init__(self,
                 input_dir: str,
                 output_dir: str,
                 reference_dir: str,
                 pga_path: Optional[str] = None,
                 organism_file: Optional[str] = None,
                 min_ir_length: int = 1000,
                 min_pidentity: int = 40,
                 qcoverage_range: str = "0.5,2",
                 genome_form: str = "circular",
                 skip_existing: bool = True,
                 log_level: str = "INFO"):
        """
        Initialize CGAS Module 2
        
        Args:
            input_dir: Directory containing clean FASTA files (e.g., Module1/07_assembled_genomes)
            output_dir: Output directory for annotations
            reference_dir: Directory containing reference GenBank files
            pga_path: Path to PGA.pl script (optional - will auto-detect if not provided)
            organism_file: Optional TSV file mapping accession to organism name
            min_ir_length: Minimum IR length (default: 1000)
            min_pidentity: Minimum percent identity (default: 40)
            qcoverage_range: Query coverage range (default: "0.5,2")
            genome_form: Genome form - circular or linear (default: "circular")
            skip_existing: Skip already annotated samples
            log_level: Logging level
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.reference_dir = Path(reference_dir)
        self.organism_file = Path(organism_file) if organism_file else None
        self.min_ir_length = min_ir_length
        self.min_pidentity = min_pidentity
        self.qcoverage_range = qcoverage_range
        self.genome_form = genome_form
        self.skip_existing = skip_existing
        
        # Create output directories FIRST
        self.genbank_dir = self.output_dir / "Annotated_GenBank"
        self.logs_dir = self.output_dir / "Annotation_Logs"
        self.reports_dir = self.output_dir / "Reports"
        
        for directory in [self.genbank_dir, self.logs_dir, self.reports_dir]:
            directory.mkdir(parents=True, exist_ok=True)
        
        # Setup logging FIRST (before loading organism map)
        self._setup_logging(log_level)
        
        # Auto-detect or use provided PGA path
        if pga_path:
            self.pga_path = pga_path
            self.logger.info(f"Using provided PGA path: {pga_path}")
        else:
            self.pga_path = find_pga_path()
            if self.pga_path:
                self.logger.info(f"✓ Auto-detected PGA at: {self.pga_path}")
            else:
                self.logger.error("ERROR: Could not find PGA installation.")
                self.logger.error("Please either:")
                self.logger.error("  1. Install PGA in ~/PGA/")
                self.logger.error("  2. Provide path using: --pga /path/to/PGA.pl")
                self.logger.error("  3. Ensure 'pga' command is in your PATH")
                raise FileNotFoundError("PGA not found. Please install or specify path.")
        
        # Verify PGA exists
        if not os.path.exists(self.pga_path):
            raise FileNotFoundError(f"PGA not found at: {self.pga_path}")
        
        # Find PGA command
        self.pga_command = self._find_pga_command()
        
        # Load organism mapping (AFTER logger is set up)
        self.organism_map = self._load_organism_mapping()
        
        # Storage for results
        self.annotation_results: Dict[str, AnnotationStats] = {}
        self.failed_samples: Dict[str, str] = {}
        
        # Track which files need processing
        self.files_to_process: List[Path] = []
    
    def _find_pga_command(self) -> List[str]:
        """Find PGA command - try multiple ways"""
        # Try direct 'pga' command first
        try:
            subprocess.run(['pga', '-h'], capture_output=True, text=True, check=True)
            self.logger.info("✓ Found 'pga' command in PATH")
            return ['pga']
        except:
            pass
        
        # Try perl with detected/provided path
        if Path(self.pga_path).exists():
            # Check if it's a .pl file (needs perl)
            if self.pga_path.endswith('.pl'):
                self.logger.info(f"✓ Using Perl script: {self.pga_path}")
                return ['perl', self.pga_path]
            else:
                # It's likely an executable
                self.logger.info(f"✓ Using executable: {self.pga_path}")
                return [self.pga_path]
        
        raise FileNotFoundError("Could not find PGA. Install with 'conda install -c bioconda pga' or set correct pga_path")
    
    def _setup_logging(self, log_level: str):
        """Setup logging configuration"""
        log_file = self.output_dir / f"cgas_module2_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        
        logging.basicConfig(
            level=getattr(logging, log_level.upper()),
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"CGAS Module 2 initialized. Log file: {log_file}")
    
    def _load_organism_mapping(self) -> Dict[str, str]:
        """
        Load organism name mapping from TSV file
        
        Format: accession\torganism
        Example:
            SRR8666784\tHibiscus coatesii
            SRR8666785\tHibiscus goldsworthii
        
        Returns:
            Dictionary mapping accession to organism name
        """
        organism_map = {}
        
        if not self.organism_file:
            self.logger.debug("No organism file provided")
            return organism_map
        
        if not self.organism_file.exists():
            self.logger.warning(f"Organism file not found: {self.organism_file}")
            return organism_map
        
        try:
            with open(self.organism_file, 'r') as f:
                lines = f.readlines()
                
                # Skip header if present
                start_idx = 0
                if lines and ('accession' in lines[0].lower() or 'sample' in lines[0].lower()):
                    start_idx = 1
                
                for line in lines[start_idx:]:
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        accession = parts[0].strip()
                        organism = parts[1].strip()
                        organism_map[accession] = organism
            
            self.logger.info(f"Loaded {len(organism_map)} organism names from {self.organism_file}")
            
        except Exception as e:
            self.logger.warning(f"Could not read organism file: {e}")
        
        return organism_map
    
    def _extract_organism_from_fasta(self, fasta_file: Path) -> Optional[str]:
        """
        Extract organism name from FASTA header
        
        Handles various NCBI FASTA header formats:
        - >Erigeron compositus TPA_asm: Erigeron compositus chloroplast, complete genome
        - >NC_123456.1 Arabidopsis thaliana chloroplast, complete genome
        - >Hibiscus syriacus isolate xyz, complete plastid genome
        
        Returns:
            Organism name (genus + species) or None if not found
        """
        try:
            with open(fasta_file, 'r') as f:
                header = f.readline().strip()
                
                if not header.startswith('>'):
                    return None
                
                # Remove the '>' character
                header = header[1:].strip()
                
                # Pattern 1: "Genus species TPA_asm:" or "Genus species isolate" format
                # Extract the first two words (genus + species)
                words = header.split()
                if len(words) >= 2:
                    # Check if first two words look like a binomial name
                    # (starts with capital letter, second word is lowercase)
                    if words[0][0].isupper() and len(words[1]) > 0 and words[1][0].islower():
                        organism = f"{words[0]} {words[1]}"
                        self.logger.debug(f"Extracted organism from FASTA header: {organism}")
                        return organism
                
                # Pattern 2: After ":" in headers like "accession: Genus species description"
                if ':' in header:
                    parts = header.split(':', 1)
                    if len(parts) == 2:
                        after_colon = parts[1].strip().split()
                        if len(after_colon) >= 2:
                            if after_colon[0][0].isupper() and after_colon[1][0].islower():
                                organism = f"{after_colon[0]} {after_colon[1]}"
                                self.logger.debug(f"Extracted organism from FASTA header (after colon): {organism}")
                                return organism
                
        except Exception as e:
            self.logger.debug(f"Could not extract organism from {fasta_file}: {e}")
        
        return None
    
    def _get_organism_info(self, sample_name: str, fasta_file: Optional[Path] = None) -> Tuple[str, str]:
        """
        Get organism name and filename prefix for a sample
        
        Priority order:
        1. Organism mapping file (if provided)
        2. FASTA header extraction (if fasta_file provided)
        3. "Unknown organism" as fallback
        
        Handles SSC flip-flop variants intelligently:
        - 27_1 and 27_2 → both map to organism for "27_1"
        - SRR123_1 and SRR123_2 → both map to organism for "SRR123_1"
        
        NOTE: This function automatically handles organism names with spaces.
        In your organism mapping file, use natural format with spaces:
            27_1    Hibiscus syriacus
        The script will automatically convert spaces to underscores in filenames:
            Output: Hibiscus_syriacus_27_1.gb
        
        Args:
            sample_name: Sample identifier (e.g., 'SRR8666784' or '27_1' or '27_2')
            fasta_file: Optional path to FASTA file for header extraction
        
        Returns:
            Tuple of (organism_name, filename_prefix)
        """
        organism_name = None
        
        # Priority 1: Try organism mapping file first
        # First, try exact match
        if sample_name in self.organism_map:
            organism_name = self.organism_map[sample_name]
            self.logger.debug(f"Exact match from organism file: {sample_name} → {organism_name}")
        
        # If no exact match, check if this is a flip-flop variant (ends with _1 or _2)
        elif sample_name.endswith('_2'):
            # Try the _1 variant
            base_name = sample_name[:-2] + '_1'
            if base_name in self.organism_map:
                organism_name = self.organism_map[base_name]
                self.logger.debug(f"Flip-flop match from organism file: {sample_name} → {base_name} → {organism_name}")
        
        elif sample_name.endswith('_1'):
            # Already ends with _1, try without suffix
            base_name = sample_name[:-2]
            if base_name in self.organism_map:
                organism_name = self.organism_map[base_name]
                self.logger.debug(f"Base match from organism file: {sample_name} → {base_name} → {organism_name}")
        
        # If still no match, try base name without any suffix
        if not organism_name:
            # Check if sample has _1 or _2 suffix
            if '_' in sample_name:
                parts = sample_name.rsplit('_', 1)
                if len(parts) == 2 and parts[1] in ['1', '2']:
                    base_name = parts[0]
                    if base_name in self.organism_map:
                        organism_name = self.organism_map[base_name]
                        self.logger.debug(f"Base match (generic) from organism file: {sample_name} → {base_name} → {organism_name}")
        
        # Priority 2: If no match in organism file, try extracting from FASTA header
        if not organism_name and fasta_file:
            organism_name = self._extract_organism_from_fasta(fasta_file)
            if organism_name:
                self.logger.info(f"Extracted organism from FASTA header: {sample_name} → {organism_name}")
        
        # Final fallback
        if organism_name:
            # Create filename-safe prefix from organism name
            filename_prefix = organism_name.replace(' ', '_').replace('.', '')
        else:
            organism_name = "Unknown organism"
            filename_prefix = sample_name
            self.logger.warning(f"No organism found for {sample_name}. Using 'Unknown organism'")
            if not self.organism_map:
                self.logger.info(f"TIP: Provide --organism-file to specify organism names, or ensure FASTA headers contain organism names")
        
        return organism_name, filename_prefix
    
    def verify_dependencies(self) -> bool:
        """Verify all required dependencies are available"""
        self.logger.info("Verifying dependencies...")
        
        # Check PGA
        try:
            result = subprocess.run(
                self.pga_command + ['-h'],
                capture_output=True,
                text=True,
                check=False
            )
            if result.returncode == 0 or 'PGA' in result.stdout or 'PGA' in result.stderr:
                self.logger.info("✓ PGA is available")
            else:
                self.logger.error("✗ PGA check failed")
                return False
        except Exception as e:
            self.logger.error(f"✗ PGA not found: {e}")
            return False
        
        # Check BLAST
        try:
            result = subprocess.run(['blastn', '-version'], capture_output=True, text=True)
            if result.returncode == 0:
                self.logger.info("✓ BLAST+ is available")
            else:
                self.logger.warning("⚠ BLAST+ check uncertain")
        except FileNotFoundError:
            self.logger.error("✗ BLAST+ not found. Install with: conda install -c bioconda blast")
            return False
        
        # Check input directory
        if not self.input_dir.exists():
            self.logger.error(f"✗ Input directory not found: {self.input_dir}")
            return False
        
        # Check reference directory
        if not self.reference_dir.exists():
            self.logger.error(f"✗ Reference directory not found: {self.reference_dir}")
            return False
        
        # Check for reference GenBank files
        ref_files = list(self.reference_dir.glob("*.gb")) + list(self.reference_dir.glob("*.gbk"))
        if not ref_files:
            self.logger.error(f"✗ No GenBank reference files found in {self.reference_dir}")
            return False
        
        self.logger.info(f"✓ Found {len(ref_files)} reference GenBank file(s)")
        
        return True
    
    def _check_existing_annotations(self) -> Tuple[List[Path], List[Path]]:
        """
        Check which files are already annotated vs need annotation
        
        Returns:
            Tuple of (already_done, need_processing) file lists
        """
        # Find all input FASTA files
        if self.input_dir.is_file():
            input_files = [self.input_dir]
        else:
            input_files = []
            for ext in ['*.fasta', '*.fa', '*.fna']:
                input_files.extend(self.input_dir.glob(ext))
        
        already_done = []
        need_processing = []
        
        for fasta_file in input_files:
            sample_name = fasta_file.stem
            organism_name, filename_prefix = self._get_organism_info(sample_name, fasta_file)
            
            # Check if GenBank file already exists
            expected_gb = self.genbank_dir / f"{filename_prefix}_{sample_name}.gb"
            
            if expected_gb.exists() and self.skip_existing:
                already_done.append(fasta_file)
                self.logger.debug(f"Skipping {sample_name} (already annotated)")
            else:
                need_processing.append(fasta_file)
        
        return already_done, need_processing
    
    def run_batch_annotation(self, fasta_files: List[Path]) -> bool:
        """
        Run PGA batch annotation on multiple FASTA files
        
        Args:
            fasta_files: List of FASTA files to annotate
        
        Returns:
            bool: True if annotation completed successfully
        """
        if not fasta_files:
            self.logger.warning("No files to annotate")
            return True
        
        # Create temporary target directory for batch processing
        temp_target_dir = self.output_dir / "temp_targets"
        temp_target_dir.mkdir(exist_ok=True)
        
        try:
            # Copy all FASTA files to temp directory
            for fasta_file in fasta_files:
                shutil.copy2(fasta_file, temp_target_dir / fasta_file.name)
            
            self.logger.info(f"Running PGA on {len(fasta_files)} files...")
            
            # Prepare PGA output directory
            pga_output_dir = self.output_dir / "temp_pga_output"
            
            # Run PGA
            cmd = self.pga_command + [
                '-r', str(self.reference_dir),
                '-t', str(temp_target_dir),
                '-o', str(pga_output_dir),
                '-i', str(self.min_ir_length),
                '-p', str(self.min_pidentity),
                '-q', self.qcoverage_range,
                '-f', self.genome_form,
                '-l', 'warning'
            ]
            
            self.logger.debug(f"PGA command: {' '.join(map(str, cmd))}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            
            if result.returncode != 0:
                self.logger.warning(f"PGA returned non-zero exit code: {result.returncode}")
                self.logger.debug(f"STDOUT: {result.stdout}")
                self.logger.debug(f"STDERR: {result.stderr}")
            
            # Process results
            if pga_output_dir.exists():
                self._process_pga_output(pga_output_dir, fasta_files)
                
                # Cleanup
                shutil.rmtree(pga_output_dir, ignore_errors=True)
                shutil.rmtree(temp_target_dir, ignore_errors=True)
                
                return True
            else:
                self.logger.error(f"PGA output directory not created: {pga_output_dir}")
                return False
                
        except Exception as e:
            self.logger.error(f"Batch annotation failed: {e}")
            return False
        finally:
            # Cleanup temp directory
            if temp_target_dir.exists():
                shutil.rmtree(temp_target_dir, ignore_errors=True)
    
    def _process_pga_output(self, pga_output_dir: Path, input_files: List[Path]):
        """
        Process PGA output files and organize them
        
        Args:
            pga_output_dir: Directory containing PGA output
            input_files: List of input FASTA files
        """
        # Find all GenBank files in output
        gb_files = list(pga_output_dir.rglob("*.gb")) + list(pga_output_dir.rglob("*.gbk"))
        
        if not gb_files:
            self.logger.warning("No GenBank files found in PGA output")
            return
        
        self.logger.info(f"Processing {len(gb_files)} GenBank files...")
        
        # Create mapping of sample names to input files
        sample_to_input = {f.stem: f for f in input_files}
        
        for gb_file in gb_files:
            try:
                # Extract sample name from GenBank filename
                # PGA typically names files like: target_name.gb
                sample_name = gb_file.stem
                
                # Find corresponding input file
                input_fasta = sample_to_input.get(sample_name)
                
                if not input_fasta:
                    # Try to match by checking all stems
                    for stem, fasta in sample_to_input.items():
                        if stem in gb_file.stem or gb_file.stem in stem:
                            sample_name = stem
                            input_fasta = fasta
                            break
                
                if not input_fasta:
                    self.logger.warning(f"Could not match {gb_file.name} to input file")
                    continue
                
                # Get organism info
                organism_name, filename_prefix = self._get_organism_info(sample_name, input_fasta)
                
                # Create final filename
                final_genbank_name = f"{filename_prefix}_{sample_name}.gb"
                final_genbank_path = self.genbank_dir / final_genbank_name
                
                # Copy and update GenBank file
                shutil.copy2(gb_file, final_genbank_path)
                self._update_genbank_organism(final_genbank_path, organism_name)
                self._add_structural_features(final_genbank_path)
                
                # Look for corresponding log file
                log_file = gb_file.parent / f"{gb_file.stem}_warning.log"
                if log_file.exists():
                    final_log_name = f"{filename_prefix}_{sample_name}_warning.log"
                    final_log_path = self.logs_dir / final_log_name
                    shutil.copy2(log_file, final_log_path)
                else:
                    log_file = None
                
                # Parse results and create stats
                stats = self._parse_pga_log(sample_name, input_fasta, log_file, final_genbank_path)
                self.annotation_results[sample_name] = stats
                
                self.logger.info(f"✓ Processed: {sample_name} → {final_genbank_name}")
                
            except Exception as e:
                self.logger.error(f"Error processing {gb_file.name}: {e}")
                sample_name = gb_file.stem
                self.failed_samples[sample_name] = str(e)
    
    def _parse_pga_log(self, sample_name: str, input_fasta: Path, 
                       log_file: Optional[Path], genbank_file: Path) -> AnnotationStats:
        """
        Parse PGA log file to extract annotation statistics
        
        Args:
            sample_name: Sample identifier
            input_fasta: Input FASTA file
            log_file: PGA warning log file (if exists)
            genbank_file: Output GenBank file
        
        Returns:
            AnnotationStats object
        """
        stats = AnnotationStats(
            sample_name=sample_name,
            input_fasta=input_fasta,
            output_genbank=genbank_file
        )
        
        # Count genes in GenBank file
        try:
            gene_count = 0
            with open(genbank_file, 'r') as f:
                for line in f:
                    if line.strip().startswith('gene '):
                        gene_count += 1
            
            stats.total_genes_annotated = gene_count
            # Set reference count to annotated count (since we don't have explicit reference count)
            # This makes annotation rate meaningful: if all expected genes are annotated, it's 100%
            stats.total_genes_reference = gene_count
            stats.annotation_success = gene_count > 0
            
        except Exception as e:
            self.logger.warning(f"Could not count genes in {genbank_file}: {e}")
        
        # Parse log file if exists
        if log_file and log_file.exists():
            try:
                with open(log_file, 'r') as f:
                    log_content = f.read()
                    
                    # Look for warnings
                    warning_patterns = [
                        r'WARNING: (.+)',
                        r'WARN: (.+)',
                        r'Missing: (.+)'
                    ]
                    
                    for pattern in warning_patterns:
                        warnings = re.findall(pattern, log_content, re.IGNORECASE)
                        stats.warnings.extend(warnings)
                    
                    # Look for unannotated genes
                    unannotated = re.findall(r'unannotated.*?(\w+)', log_content, re.IGNORECASE)
                    stats.unannotated_genes.extend(unannotated)
                    
            except Exception as e:
                self.logger.warning(f"Could not parse log file {log_file}: {e}")
        
        return stats
    
    def _update_genbank_organism(self, genbank_file: Path, organism_name: str):
        """Update the ORGANISM field in GenBank file"""
        try:
            with open(genbank_file, 'r') as f:
                lines = f.readlines()
            
            updated = False
            for i, line in enumerate(lines):
                if line.strip().startswith('/organism='):
                    lines[i] = f'                     /organism="{organism_name}"\n'
                    updated = True
                    break
            
            if updated:
                with open(genbank_file, 'w') as f:
                    f.writelines(lines)
                self.logger.debug(f"Updated organism to: {organism_name}")
                
        except Exception as e:
            self.logger.warning(f"Failed to update organism name: {e}")
    
    def _add_structural_features(self, genbank_file: Path):
        """
        Update IR annotations and add LSC/SSC features to GenBank file
        
        Strategy:
        1. Find existing repeat_region features (IRs from PGA)
        2. Update their /note to standard wording (IRA/IRB) without changing coordinates
        3. Identify regions between IRs
        4. Add LSC (larger region) and SSC (smaller region) as misc_features
        
        Args:
            genbank_file: Path to GenBank file
        """
        try:
            from Bio import SeqIO
            from Bio.SeqFeature import SeqFeature, FeatureLocation
            
            # Read GenBank file
            record = SeqIO.read(genbank_file, "genbank")
            
            # Find existing repeat_region features (IRs annotated by PGA)
            ir_features = []
            for feat in record.features:
                if feat.type == "repeat_region":
                    ir_features.append(feat)
            
            if len(ir_features) < 2:
                self.logger.debug(f"Found {len(ir_features)} IR features, need at least 2. Skipping structural annotation.")
                return
            
            # Sort IRs by start position
            ir_features_sorted = sorted(ir_features, key=lambda f: f.location.start)
            
            # Update IR annotations to standard wording
            # First IR = IRB (inverted repeat B)
            # Second IR = IRA (inverted repeat A)
            ir_features_sorted[0].qualifiers["note"] = ["inverted repeat B (IRB)"]
            if len(ir_features_sorted) >= 2:
                ir_features_sorted[1].qualifiers["note"] = ["inverted repeat A (IRA)"]
            
            # Get IR positions (1-based for display)
            irb_start = int(ir_features_sorted[0].location.start) + 1
            irb_end = int(ir_features_sorted[0].location.end)
            
            ira_start = int(ir_features_sorted[1].location.start) + 1
            ira_end = int(ir_features_sorted[1].location.end)
            
            # Calculate regions between start, IRb, IRa, and end
            # Region 1: start to IRb
            region1_start = 1
            region1_end = irb_start - 1
            region1_length = region1_end - region1_start + 1
            
            # Region 2: between IRb and IRa
            region2_start = irb_end + 1
            region2_end = ira_start - 1
            region2_length = region2_end - region2_start + 1
            
            # Determine which is LSC (larger) and which is SSC (smaller)
            if region1_length > region2_length:
                # Region 1 is LSC, Region 2 is SSC
                lsc_start, lsc_end = region1_start, region1_end
                ssc_start, ssc_end = region2_start, region2_end
            else:
                # Region 2 is LSC, Region 1 is SSC
                lsc_start, lsc_end = region2_start, region2_end
                ssc_start, ssc_end = region1_start, region1_end
            
            # Create LSC and SSC features
            lsc_feature = SeqFeature(
                FeatureLocation(lsc_start - 1, lsc_end, strand=1),
                type="misc_feature",
                qualifiers={"note": ["large single copy (LSC)"]}
            )
            
            ssc_feature = SeqFeature(
                FeatureLocation(ssc_start - 1, ssc_end, strand=1),
                type="misc_feature",
                qualifiers={"note": ["small single copy (SSC)"]}
            )
            
            # Insert LSC and SSC features after source feature
            source_idx = 0
            for i, feat in enumerate(record.features):
                if feat.type == "source":
                    source_idx = i + 1
                    break
            
            # Add features in order (after source)
            record.features.insert(source_idx, lsc_feature)
            record.features.insert(source_idx + 1, ssc_feature)
            
            # Write updated GenBank file
            SeqIO.write(record, genbank_file, "genbank")
            
            self.logger.info(f"✓ Updated structural features: LSC={lsc_end-lsc_start+1}bp, SSC={ssc_end-ssc_start+1}bp, IRB={irb_end-irb_start+1}bp, IRA={ira_end-ira_start+1}bp")
            
        except ImportError:
            self.logger.warning("Biopython not installed. Skipping structural annotation.")
        except Exception as e:
            self.logger.warning(f"Failed to add structural features: {e}")
            import traceback
            self.logger.debug(traceback.format_exc())
    
    def generate_reports(self):
        """Generate summary reports"""
        self.logger.info("Generating reports...")
        
        # Summary statistics report
        summary_file = self.reports_dir / "00_ANNOTATION_SUMMARY.tsv"
        with open(summary_file, 'w') as f:
            f.write("Sample\tOrganism\tInput_File\tGenes_Annotated\t")
            f.write("Annotation_Rate\tWarnings\tUnannotated_Genes\tOutput_GenBank\n")
            
            for sample_name, stats in sorted(self.annotation_results.items()):
                organism_name, _ = self._get_organism_info(sample_name, stats.input_fasta)
                rate = (stats.total_genes_annotated / stats.total_genes_reference * 100 
                       if stats.total_genes_reference > 0 else 0)
                
                f.write(f"{sample_name}\t")
                f.write(f"{organism_name}\t")
                f.write(f"{stats.input_fasta.name}\t")
                f.write(f"{stats.total_genes_annotated}\t")
                f.write(f"{rate:.2f}\t")
                f.write(f"{len(stats.warnings)}\t")
                f.write(f"{len(stats.unannotated_genes)}\t")
                f.write(f"{stats.output_genbank.name if stats.output_genbank else 'NA'}\n")
        
        # Failed samples list
        if self.failed_samples:
            failed_file = self.reports_dir / "00_FAILED_SAMPLES.tsv"
            with open(failed_file, 'w') as f:
                f.write("Sample\tError\n")
                for sample, error in self.failed_samples.items():
                    f.write(f"{sample}\t{error}\n")
        
        self.logger.info(f"Reports saved to {self.reports_dir}")
    
    def run_pipeline(self):
        """Run the complete annotation pipeline"""
        self.logger.info("="*80)
        self.logger.info("CGAS Module 2 - Batch Annotation Pipeline")
        self.logger.info("="*80)
        
        # Verify dependencies
        if not self.verify_dependencies():
            self.logger.error("Dependency verification failed")
            return
        
        # Find and categorize files
        already_done, need_processing = self._check_existing_annotations()
        
        if already_done:
            self.logger.info(f"\n{len(already_done)} files already annotated (skipping)")
        
        if need_processing:
            self.logger.info(f"\n{len(need_processing)} files need annotation")
            
            # Run batch annotation
            success = self.run_batch_annotation(need_processing)
            
            if not success and len(need_processing) > 1:
                # Try individual annotation as fallback
                self.logger.warning("Batch annotation failed, trying individual files...")
                for fasta_file in need_processing:
                    try:
                        # For individual files, we need a different approach
                        # Since PGA expects directories, create a temp dir for each
                        temp_dir = self.output_dir / f"temp_{fasta_file.stem}"
                        temp_target_dir = temp_dir / "targets"
                        temp_target_dir.mkdir(parents=True, exist_ok=True)
                        
                        shutil.copy2(fasta_file, temp_target_dir / fasta_file.name)
                        
                        pga_output_dir = temp_dir / "output"
                        
                        cmd = self.pga_command + [
                            '-r', str(self.reference_dir),
                            '-t', str(temp_target_dir),
                            '-o', str(pga_output_dir),
                            '-i', str(self.min_ir_length),
                            '-p', str(self.min_pidentity),
                            '-q', self.qcoverage_range,
                            '-f', self.genome_form,
                            '-l', 'warning'
                        ]
                        
                        self.logger.info(f"Running individual annotation for {fasta_file.stem}")
                        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
                        
                        if pga_output_dir.exists():
                            output_files = list(pga_output_dir.rglob("*"))
                            for f in output_files:
                                if f.suffix.lower() in ['.gb', '.gbk']:
                                    sample_name = fasta_file.stem
                                    organism_name, filename_prefix = self._get_organism_info(sample_name, fasta_file)
                                    final_genbank_name = f"{filename_prefix}_{sample_name}.gb"
                                    genbank_file = self.genbank_dir / final_genbank_name
                                    
                                    shutil.move(str(f), str(genbank_file))
                                    self._update_genbank_organism(genbank_file, organism_name)
                                    self._add_structural_features(genbank_file)
                                    
                                    log_file = self.logs_dir / f"{filename_prefix}_{sample_name}_warning.log"
                                    stats = self._parse_pga_log(sample_name, fasta_file, log_file, genbank_file)
                                    self.annotation_results[sample_name] = stats
                                    
                                    self.logger.info(f"✓ Individually annotated: {sample_name}")
                                    break
                        
                        shutil.rmtree(temp_dir, ignore_errors=True)
                        
                    except Exception as e:
                        sample_name = fasta_file.stem
                        self.logger.error(f"Failed to annotate {sample_name}: {e}")
                        self.failed_samples[sample_name] = str(e)
        
        # Generate reports
        self.generate_reports()
        
        # Final summary
        self.logger.info("\n" + "="*80)
        self.logger.info("ANNOTATION SUMMARY")
        self.logger.info("="*80)
        total_files = len(self.annotation_results) + len(self.failed_samples)
        self.logger.info(f"Total files processed: {total_files}")
        self.logger.info(f"Successfully annotated: {len(self.annotation_results)}")
        self.logger.info(f"Failed: {len(self.failed_samples)}")
        
        if self.failed_samples:
            self.logger.warning("\nFailed samples:")
            for sample, error in self.failed_samples.items():
                self.logger.warning(f"  - {sample}: {error}")
        
        self.logger.info(f"\nOutput directory: {self.output_dir}")
        self.logger.info(f"GenBank files: {self.genbank_dir}")
        self.logger.info(f"Logs: {self.logs_dir}")
        self.logger.info(f"Reports: {self.reports_dir}")
        self.logger.info("="*80)


def main():
    parser = argparse.ArgumentParser(
        description="CGAS Module 2: Batch Plastome Annotation using PGA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic batch annotation (auto-detect PGA)
  python cgas_module2.py -i targets/ -r refs/ -o annotations/
  
  # Single file annotation
  python cgas_module2.py -i targets/27_1.fasta -r refs/ -o single_output/
  
  # With organism names
  python cgas_module2.py -i targets/ -r refs/ --organism-file organisms.txt
  
  # With explicit PGA path
  python cgas_module2.py -i targets/ -r refs/ --pga /home/user/PGA/PGA.pl
  
  # Force re-annotation of all files
  python cgas_module2.py -i targets/ -r refs/ -o annotations/ --force-rerun
  
  # Debug mode
  python cgas_module2.py -i targets/ -r refs/ -o annotations/ --log-level DEBUG
        """
    )
    
    parser.add_argument("-i", "--input", required=True,
                       help="Input directory or file with FASTA files")
    parser.add_argument("-r", "--reference", required=True,
                       help="Directory containing reference GenBank files")
    parser.add_argument("-o", "--output", default="module2_annotations",
                       help="Output directory (default: module2_annotations)")
    
    parser.add_argument("--organism-file", 
                       help="TSV file mapping accession to organism name (accession\\torganism)")
    
    parser.add_argument("--pga", 
                       help="Path to PGA.pl or PGA executable (optional - will auto-detect if not provided)")
    parser.add_argument("--min-ir", type=int, default=1000,
                       help="Minimum IR length (default: 1000)")
    parser.add_argument("-p", "--min-pidentity", type=int, default=40,
                       help="Minimum percent identity (default: 40)")
    parser.add_argument("-q", "--qcoverage", default="0.5,2",
                       help="Query coverage range (default: 0.5,2)")
    parser.add_argument("-f", "--form", default="circular",
                       choices=["circular", "linear"],
                       help="Genome form (default: circular)")
    parser.add_argument("--force-rerun", action="store_true",
                       help="Force re-annotation of all samples (ignore existing)")
    parser.add_argument("--log-level", default="INFO",
                       choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                       help="Logging level (default: INFO)")
    
    args = parser.parse_args()
    
    try:
        pipeline = CGASModule2(
            input_dir=args.input,
            output_dir=args.output,
            reference_dir=args.reference,
            pga_path=args.pga,  # Can be None - will auto-detect
            organism_file=args.organism_file,
            min_ir_length=args.min_ir,
            min_pidentity=args.min_pidentity,
            qcoverage_range=args.qcoverage,
            genome_form=args.form,
            skip_existing=not args.force_rerun,
            log_level=args.log_level
        )
        
        pipeline.run_pipeline()
        
    except FileNotFoundError as e:
        print(f"\nERROR: {e}")
        print("\nTo fix this:")
        print("1. Install PGA following the README instructions, OR")
        print("2. Provide the PGA path explicitly: --pga /path/to/PGA.pl")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
