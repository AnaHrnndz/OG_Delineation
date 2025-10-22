#!/usr/bin/env python3
"""
OGD: Orthologous Group Delineation Pipeline

This script serves as the main entry point for the OGD pipeline.
It parses command-line arguments and calls the core pipeline function.
"""

import sys
import argparse
import warnings
from pathlib import Path
import logging

# Local application import for the main pipeline logic
from ogd.ogd_core import run_ogd_pipeline

# --- Global Configuration ---

# Set a high recursion limit for deep phylogenetic tree traversal.
sys.setrecursionlimit(10000)

# Suppress common runtime warnings from scientific libraries.
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Configure basic logging for user feedback
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
    params.add_argument('--no_inherit_outliers', dest='inherit_outliers', action='store_false',
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

    return parser.parse_args()


# --- Main Execution ---

def main():
    """Parses arguments and executes the OGD pipeline."""
    args = get_args()
    
    try:
        # --- Input Validation ---
        if not args.tree.is_file():
            raise FileNotFoundError(f"Input tree file not found: {args.tree}")

        # --- Output Directory Creation ---
        args.out_path.mkdir(parents=True, exist_ok=True)
        
        print("\n🚀 Starting OG Delineation pipeline...")
        run_ogd_pipeline(args)
        print(f"✅ Pipeline finished successfully. Results are in: {args.out_path}\n")

    except FileNotFoundError as e:
        print(f"❌ ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        # This will catch other errors during the pipeline execution
        print(f"❌ An unexpected error occurred: {e}", file=sys.stderr)
        # For debugging, you might want to uncomment the next line to see the full traceback
        # import traceback; traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

