#!/usr/bin/env python3
"""
CGAS - Chloroplast Genome Analysis Suite (Modules 1-14)
Command-line Interface
"""

import os
import sys
import argparse
import importlib.util
from pathlib import Path

def print_banner():
    """Print CGAS banner."""
    print("=" * 80)
    print("     CGAS - CHLOROPLAST GENOME ANALYSIS SUITE v1.0.1")
    print("=" * 80)
    print()

def run_module(module_num, module_args=None):
    """Run a specific module with optional arguments."""
    print(f"\n" + "═" * 60)
    print(f"Running Module {module_num}")
    print("═" * 60)
    
    try:
        # Import the module
        module_name = f"cgas_module{module_num}"
        
        # Try to import from cgas package first
        try:
            module = __import__(f"cgas.{module_name}", fromlist=[''])
        except ImportError:
            # Fallback to direct file import
            cli_dir = Path(__file__).resolve().parent
            module_path = cli_dir / f"{module_name}.py"
            
            if not module_path.exists():
                print(f"Error: Module file not found: {module_path}")
                return False
            
            spec = importlib.util.spec_from_file_location(module_name, module_path)
            module = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = module
            spec.loader.exec_module(module)
        
        # Save original sys.argv
        original_argv = sys.argv.copy()
        
        # Set sys.argv for the module
        if module_args:
            sys.argv = [f"cgas_module{module_num}.py"] + module_args
        else:
            sys.argv = [f"cgas_module{module_num}.py"]
        
        # Run the module
        module.main()
        
        # Restore original sys.argv
        sys.argv = original_argv
        
        return True
        
    except Exception as e:
        print(f"Error running Module {module_num}: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main entry point - handles command-line arguments."""
    parser = argparse.ArgumentParser(
        description='CGAS - Chloroplast Genome Analysis Suite v1.0.1 (Modules 1-14)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  cgas --module 1          # Run Module 1
  cgas --module 2          # Run Module 2
  cgas --modules 5,6,7     # Run Modules 5, 6, and 7
  cgas --all               # Run all modules (1-14)
  cgas --list              # List available modules
  cgas --help              # Show this help message
  
  # With module-specific arguments:
  cgas --module 1 -i raw_reads/ -o results/
  cgas --module 3 -r Abutilon_grandifolium.gb
  cgas --module 14 --genes-only -og Abutilon_grandifolium.gb
  cgas --modules 5,6 -i input_dir/ -o output_dir/
  
Individual module commands:
  cgas-assembly            # Run Module 1
  cgas-annotate            # Run Module 2
  cgas-compare             # Run Module 3
  cgas-convert             # Run Module 4
  cgas-gene-compare        # Run Module 5
  cgas-gene-table          # Run Module 6
  cgas-genome-compare      # Run Module 7
  cgas-codon               # Run Module 8
  cgas-amino               # Run Module 9
  cgas-snp                 # Run Module 10
  cgas-intron              # Run Module 11
  cgas-ssr                 # Run Module 12
  cgas-diversity           # Run Module 13
  cgas-phylogeny           # Run Module 14
  
Usage:
  1. Navigate to directory with your data files
  2. Run: cgas --module <number> [module arguments]
  3. Results are saved in Module<number>_ folders
        """
    )
    
    parser.add_argument(
        '-m', '--module',
        type=int,
        choices=range(1, 15),
        help='Run a specific module (1-14)'
    )
    
    parser.add_argument(
        '--modules',
        type=str,
        help='Run specific modules (comma-separated, e.g., 1,3,5)'
    )
    
    parser.add_argument(
        '--all',
        action='store_true',
        help='Run all modules 1-14 (complete workflow)'
    )
    
    parser.add_argument(
        '--list',
        action='store_true',
        help='List all available modules and exit'
    )
    
    # Parse known args to allow module-specific arguments to pass through
    args, module_args = parser.parse_known_args()
    
    # Print banner
    print_banner()
    
    # Handle arguments
    if args.list:
        print("Available Modules (1-14):")
        modules = [
            "1.  Chloroplast Genome Assembly (from FASTQ reads)",
            "2.  Plastome Annotation (from FASTA assemblies)", 
            "3.  Plastome Gene Comparison (normalize gene names)",
            "4.  GenBank Format Conversion (for NCBI submission)",
            "5.  Gene Comparative Analysis",
            "6.  Gene Content Tables",
            "7.  Comparative Genome Analysis",
            "8.  Codon Usage Analysis (RSCU)",
            "9.  Amino Acid Analysis",
            "10. SNP/Substitution Analysis",
            "11. Gene and tRNA Intron Analysis",
            "12. Comprehensive SSR Analysis",
            "13. Nucleotide Diversity Analysis",
            "14. Phylogenetic Matrix Builder"
        ]
        for module in modules:
            print(f"  {module}")
        return
    
    elif args.module:
        # Run single module with optional arguments
        run_module(args.module, module_args if module_args else None)
        return
    
    elif args.modules:
        # Run multiple modules with optional arguments
        try:
            module_nums = [int(m.strip()) for m in args.modules.split(',')]
            for num in module_nums:
                if 1 <= num <= 14:
                    run_module(num, module_args if module_args else None)
                else:
                    print(f"Warning: Invalid module number {num}")
        except ValueError:
            print("Error: Invalid module numbers format")
            return
    
    elif args.all:
        # Run all modules
        confirm = input("Run ALL modules (1-14)? This may take a while. (y/N): ").strip().lower()
        if confirm == 'y':
            for module_num in range(1, 15):
                run_module(module_num, module_args if module_args else None)
        return
    
    else:
        # No arguments: show help
        parser.print_help()

# Individual module commands for setup.py entry points
def assembly(): 
    print_banner()
    # Get arguments after the command name
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(1, module_args)

def annotate(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(2, module_args)

def compare(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(3, module_args)

def convert(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(4, module_args)

def gene_compare(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(5, module_args)

def gene_table(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(6, module_args)

def genome_compare(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(7, module_args)

def codon(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(8, module_args)

def amino(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(9, module_args)

def snp(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(10, module_args)

def intron(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(11, module_args)

def ssr(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(12, module_args)

def diversity(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(13, module_args)

def phylogeny(): 
    print_banner()
    module_args = sys.argv[1:] if len(sys.argv) > 1 else None
    run_module(14, module_args)

if __name__ == '__main__':
    main()
