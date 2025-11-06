#!/usr/bin/env python3
"""
OGD: Orthologous Group Delineation Pipeline
"""

import sys
import argparse
import warnings
import logging
from pathlib import Path
import traceback
from ete4 import Tree

# Local application import for the main pipeline logic
from ogd.ogd_core import run_ogd_pipeline

# --- Global Configuration ---

# Set a high recursion limit for deep phylogenetic tree traversal
sys.setrecursionlimit(10000)

# Suppress common runtime warnings from scientific libraries
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Configure basic logging for user feedback
# The format ensures all messages are clear and timestamped.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# --- Argument Parsing ---

def get_args() -> argparse.Namespace:
    """Defines and parses command-line arguments for the OGD pipeline."""
    
    parser = argparse.ArgumentParser(
        description="OGD: A tool for delineating Orthologous Groups from gene phylogenies.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-v', '--version', action='version', version='OGD 1.0')

    # --- Required Arguments ---
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--tree', required=True, type=Path,
                        help="Path to the input gene tree file (e.g., Newick).")
    required.add_argument('--output_path', dest='out_path', required=True, type=Path,
                        help="Path to the output directory.")

    # --- Taxonomy Options ---
    taxonomy = parser.add_argument_group('Taxonomy options')
    taxonomy.add_argument('--taxonomy_type', choices=['NCBI', 'GTDB'], default='NCBI',
                        help="Taxonomy database to use (default: NCBI).")
    taxonomy.add_argument('--user_taxonomy', type=Path, default=None,
                        help="Path to a custom local taxonomy database file.")
    taxonomy.add_argument('--user_taxonomy_counter', type=Path, default=None,
                        help="Path to a pre-computed taxonomy counter JSON file.")
    taxonomy.add_argument('--reftree', type=Path, default=None,
                        help="Path to a user-provided species reference tree.")
    taxonomy.add_argument('--sp_delimitator', dest='sp_delim', default='.',
                        help="Delimiter for species ID in leaf names (default: '.').")

    # --- Delineation Parameters ---
    params = parser.add_argument_group('Delineation algorithm parameters')
    params.add_argument('--sp_ovlap_all', default=0.1, dest='so_all', type=float,
                        help='Default species overlap threshold (default: 0.1).')
    params.add_argument('--sp_ovlap_euk', dest='so_euk', type=float,
                        help="SO threshold for Eukaryote-only nodes.")
    params.add_argument('--sp_ovlap_bact', dest='so_bact', type=float,
                        help="SO threshold for Bacteria-only nodes.")
    params.add_argument('--sp_ovlap_arq', dest='so_arq', type=float,
                        help="SO threshold for Archaea-only nodes.")
    params.add_argument('--lineage_threshold', dest='lineage_thr', default=0.05, type=float,
                        help="Threshold for taxonomic outlier detection (default: 0.05).")
    params.add_argument('--best_taxa_threshold', dest='best_tax_thr', default=0.9, type=float,
                        help="Threshold for LCA annotation confidence (default: 0.9).")
    params.add_argument('--species_losses_perct', dest='sp_loss_perc', default=0.7, type=float,
                        help="Threshold for species loss percentage (default: 0.7).")
    params.add_argument('--no_inherit_outliers',  action='store_true',
                        help="Flag to prevent outliers from being inherited by parent nodes.")
    
    # --- Tree Processing Options ---
    tree_opts = parser.add_argument_group('Tree processing options')
    tree_opts.add_argument('--rooting', choices=['Midpoint', 'MinVar'], default='Midpoint',
                         help="Rooting method for the tree (default: Midpoint).")

    # --- Optional Modules ---
    modules = parser.add_argument_group('Optional modules')
    modules.add_argument('--raw_alg', dest='alg', type=Path,
                         help="Path to alignment file (required for recovery and emapper).")
    modules.add_argument('--run_recovery', choices=["run-align", "skip-align"],
                         help="Enable the sequence recovery module.")
    modules.add_argument('--run_emapper', action='store_true',
                         help="Run emapper for functional annotation.")
    modules.add_argument('--emapper_dmnd_db', dest='emapper_dmnd', type=Path,
                         help="Path to emapper's DIAMOND database.")
    modules.add_argument('--emapper_pfam_db', dest='emapper_pfam', type=Path,
                         help="Path to emapper's PFAM database.")
    modules.add_argument('--path2emapper_main', type=Path,
                         help="Annotate with pre-computed emapper results (main table).")
    modules.add_argument('--path2emapper_pfams', type=Path,
                         help="Annotate with pre-computed emapper results (PFAM table).")
    modules.add_argument('--skip_get_pairs', action='store_true',
                         help="Skip the generation of the ortholog pairs table.")
    modules.add_argument('--open_visualization', action='store_true',
                         help="Automatically open the ETE smart view GUI after analysis.")
    modules.add_argument('--only_visualization', action='store_true',
                         help="Automatically open the ETE smart view GUI after analysis.")

    return parser.parse_args()

def _validate_args(args: argparse.Namespace):
    """
    Validates that all required files and dependent arguments are logical.
    Raises FileNotFoundError or ValueError if a check fails.
    """
    logging.info("Validating input arguments...")
    if not args.tree.is_file():
        raise FileNotFoundError(f"Input tree file not found: {args.tree}")

    # --- Check dependencies for optional modules ---
    alignment_needed = args.run_recovery or args.run_emapper
    
    if alignment_needed and not args.alg:
        raise ValueError("--run_recovery or --run_emapper requires --raw_alg to be set.")
    
    if args.alg and not args.alg.is_file():
        raise FileNotFoundError(f"Alignment file not found: {args.alg}")

    # --- Check emapper dependencies ---
    if args.run_emapper:
        if not args.emapper_dmnd or not args.emapper_pfam:
            raise ValueError("--run_emapper requires --emapper_dmnd_db and --emapper_pfam_db.")
        if not args.emapper_dmnd.is_file():
            raise FileNotFoundError(f"Emapper diamond DB not found: {args.emapper_dmnd}")
        if not args.emapper_pfam.is_file():
            raise FileNotFoundError(f"Emapper Pfam DB not found: {args.emapper_pfam}")

    # --- Check pre-computed emapper dependencies ---
    if args.path2emapper_main or args.path2emapper_pfams:
        if not (args.path2emapper_main and args.path2emapper_pfams):
            raise ValueError("Both --path2emapper_main and --path2emapper_pfams must be provided together.")
        if not args.path2emapper_main.is_file():
            raise FileNotFoundError(f"Emapper main table not found: {args.path2emapper_main}")
        if not args.path2emapper_pfams.is_file():
            raise FileNotFoundError(f"Emapper Pfam table not found: {args.path2emapper_pfams}")

    logging.info("Arguments validated successfully.")

# --- Main Execution ---

def main():
    """Parses arguments, validates them, and executes the OGD pipeline."""
    args = get_args()
    
    try:
        # --- 1. Argument Validation ---
        _validate_args(args)

        # --- 2. Output Directory Creation ---
        args.out_path.mkdir(parents=True, exist_ok=True)
        
        if args.only_visualization:
            from ogd.run_ete4_smartview import run_smartview
            processed_tree = Tree(open(args.tree))
            alignment_path = str(args.alg) if args.alg else None
            run_smartview(processed_tree, alignment_path)
            sys.exit(1)
            

        # --- 3. Run Pipeline ---
        # Unified logging for pipeline status
        logging.info("\n🚀 Starting OG Delineation pipeline...")
        run_ogd_pipeline(args)
        logging.info(f"\n✅ Pipeline finished successfully. Results are in: {args.out_path}\n")

    except (FileNotFoundError, ValueError) as e:
        # Errors related to configuration or missing files
        logging.error(f"❌ CONFIGURATION ERROR: {e}")
        sys.exit(1)
    except OSError as e:
        # Errors related to file system operations (e.g., creating directory)
        logging.error(f"❌ FILE SYSTEM ERROR: {e}")
        sys.exit(1)
    except Exception as e:
        # Catch-all for any other unexpected error during the pipeline
        logging.error(f"❌ An unexpected error occurred: {e}")
        # Provide full traceback in debug mode (if user sets logging level to DEBUG)
        logging.debug(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()