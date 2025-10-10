#!/usr/bin/env python3

import sys
import os
import argparse
import warnings
import pathlib

# Import the core logic function
from ogd.ogd_core import run_ogd 

# Set recursion limit early for safety in tree traversal
sys.setrecursionlimit(10000)

warnings.filterwarnings("ignore", category=RuntimeWarning)


# --- Argument Parsing ---

def get_args():
    """Defines and parses command-line arguments for the OGD pipeline."""
    
    parser = argparse.ArgumentParser(
        description="OG_Delineation: Orthologous Group Delineation (OGD) Pipeline. Identifies OGs based on duplication events."
    )
    
    # 1. Pipeline Control Arguments (Keep Required first)
    parser.add_argument('--tree', dest='tree', required=True,
                        help="Input tree file path.")
    parser.add_argument('--output_path', dest='out_path', required=True,
                        help="Output path for results.")
    
    # 2. Thresholds (Clean up the help text for better formatting)
    thresholds = parser.add_argument_group('Species Overlap Thresholds')
    thresholds.add_argument('--sp_ovlap_all', default=0.1, dest='so_all', type=float,
                            help='Species overlap threshold used for all nodes.')
    thresholds.add_argument('--sp_ovlap_euk', dest='so_euk', type=float,
                            help="SO threshold for Eukaryotic-only internal nodes.")
    thresholds.add_argument('--sp_ovlap_bact', dest='so_bact', type=float,
                            help="SO threshold for Bacterial-only internal nodes.")
    thresholds.add_argument('--sp_ovlap_arq', dest='so_arq', type=float,
                            help="SO threshold for Archaeal-only internal nodes.")
                            
    # 3. Core Parameters
    core_params = parser.add_argument_group('Core OGD Parameters')
    core_params.add_argument('--lineage_threshold', dest='lineage_thr', default=0.05, type=float)
    core_params.add_argument('--best_taxa_threshold', dest='best_tax_thr', default=0.9, type=float)
    core_params.add_argument('--inherit_outliers', dest='inherit_out', choices=['Yes', 'No'], default='Yes', type=str)
    core_params.add_argument('--species_losses_perct', dest='sp_loss_perc', default=0.7, type=float)
    core_params.add_argument('--rooting', choices=['Midpoint', 'MinVar', ''], help="Method used for tree rooting.")
    
    # 4. Taxonomy & References
    tax_ref = parser.add_argument_group('Taxonomy and References')
    tax_ref.add_argument('--taxonomy_type', choices=['NCBI', 'GTDB'], default='NCBI', help="Taxonomy database to use.")
    tax_ref.add_argument('--user_taxonomy', default=None, help="Path to user-defined taxonomy file.")
    tax_ref.add_argument('--user_taxonomy_counter', default=None)
    tax_ref.add_argument('--reftree', default=None, help="Path to reference species tree.")
    tax_ref.add_argument('--sp_delimitator', dest='sp_delim', default='.', help="Character used to separate gene name from species name.")
    
    # 5. Optional Pipelines
    opt_pipes = parser.add_argument_group('Optional Modules')
    opt_pipes.add_argument('--raw_alg', dest='alg', help="Input alignment. Needed for emapper and recovery modules.")
    opt_pipes.add_argument('--run_emapper', action='store_true', help="Run the eggNOG-mapper annotation pipeline.")
    opt_pipes.add_argument('--emapper_dmnd_db', dest='emapper_dmnd', help="Path to emapper Diamond DB.")
    opt_pipes.add_argument('--emapper_pfam_db', dest='emapper_pfam', help="Path to emapper Pfam DB.")
    opt_pipes.add_argument('--run_treeprofiler_emapper_annotation', dest='path2emapper_main')
    opt_pipes.add_argument('--run_treeprofiler_emapper_pfams', dest='path2emapper_pfams')
    opt_pipes.add_argument('--run_recovery', dest='run_recovery', choices=["run-align", "skip-align"], help="Run the sequence recovery pipeline.")
    opt_pipes.add_argument('--skip_get_pairs', action='store_true', help="Skip the step of generating all ortholog pairs.")
    opt_pipes.add_argument('--open_visualization', action='store_true', help="Open the visualization tool after completion.")
    
    # 6. Utility
    parser.add_argument('-v', '--version', action='store_true', help="Show version and exit.")


    args = parser.parse_args()
    
    # NOTE: You should implement the version check here
    # if args.version:
    #    print("OGD Version 1.0")
    #    sys.exit(0)

    return args


# --- Main Execution ---

def main():
    """
    Parses arguments and initiates the OGD pipeline execution.
    """
    args = get_args()

    run_ogd(args)


if __name__ == "__main__":
    main()

