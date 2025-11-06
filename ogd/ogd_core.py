"""
Core pipeline logic for Orthologous Group Delineation (OGD).
"""

import argparse
import logging
from pathlib import Path

# --- Local Application Imports ---

import ogd.utils as utils
from ogd.timer import Timer
from ogd.handle_input import load_all_input_files
from ogd.tree_setup import run_setup
from ogd.outliers_scores import run_outliers_and_scores
from ogd..s import run_get_main_dups
from ogd.orthologs_groups import get_all_ogs
import ogd.pairs as pairs
from ogd.emapper_annotate import annotate_with_emapper
from ogd.recovery import recover_sequences
import ogd.prepare_outputs as prepare_outputs 
from ogd.run_ete4_smartview import run_smartview

# ---  Main Pipeline Function ---

def run_ogd_pipeline(args: argparse.Namespace):
    """
    Executes the entire OG delineation pipeline .
    """
    total_timer = Timer('Total_Time')
    total_timer.start()

    logging.info("Pipeline started. Preparing run...")
    clean_tree_name = utils.remove_file_extension(args.tree.name)
    temp_dir = utils.create_temp_dir(args.out_path)

    # --- Step 1: Load All Input Data ---
    logging.info("\n--- Step 1: Loading Input Data ---")
    phylo_tree, tax_db, species_tree, level2species = load_all_input_files(args)
    
    # --- Step 2: Tree Preprocessing ---
    logging.info("\n--- Step 2: Preprocessing Tree (Setup) ---")
    processed_tree, species_set, total_seqs, num_total_species = run_setup(
        phylo_tree, tax_db, temp_dir, args
    )
    
    # --- Step 3: Outlier Detection and Scoring ---
    logging.info("\n--- Step 3: Detecting Outliers and Scoring Nodes ---")
    processed_tree, _, total_outliers = run_outliers_and_scores(
        processed_tree, tax_db, num_total_species, level2species, args
    )
    
    # --- Step 4: Identify High-Quality Duplications ---
    logging.info("\n--- Step 4: Identifying High-Quality Duplications ---")
    processed_tree, taxids_defining_ogs = run_get_main_dups(processed_tree, tax_db, args)
    
    
    # --- Step 5: Delineate Orthologous Groups ---
    logging.info("\n--- Step 5: Delineating Orthologous Groups (OGs) ---")
    processed_tree, ogs_info, seqs_in_ogs = get_all_ogs(processed_tree, tax_db)

    # --- Step 6: (Optional) Ortholog Pairs ---
    if not args.skip_get_pairs:
        logging.info("\n--- Step 6 (Optional): Extracting Ortholog Pairs ---")
        clean_pairs, strict_pairs = pairs.get_all_pairs(processed_tree, total_outliers)
        pairs.write_pairs_table(clean_pairs, strict_pairs, args.out_path, clean_tree_name)
    else:
        logging.info("\n--- Step 6 (Optional): Skipped Ortholog Pair Extraction ---")

    # --- Step 7: (Optional) Sequence Recovery ---
    recovered_seqs = set()
    seqs_to_recover = total_seqs - seqs_in_ogs
    if args.run_recovery and args.alg and seqs_to_recover and seqs_in_ogs:
        logging.info(f"\n--- Step 7 (Optional): Running Sequence Recovery for {len(seqs_to_recover)} unassigned sequences ---")
        recovered_seqs, ogs_info = recover_sequences(
            processed_tree, str(args.alg), ogs_info, seqs_to_recover, args.tree.name, args.out_path, args.run_recovery
        )
        seqs_in_ogs.update(recovered_seqs)
    elif args.run_recovery:
        logging.info("\n--- Step 7 (Optional): Skipped Sequence Recovery (no sequences to recover or no OGs found) ---")
    else:
        logging.info("\n--- Step 7 (Optional): Skipped Sequence Recovery ---")

    # --- Step 8: (Optional) Emapper Annotation ---
    if args.run_emapper:
        logging.info("\n--- Step 8 (Optional): Running eggNOG-Mapper Annotation ---")
        processed_tree = annotate_with_emapper(processed_tree, str(args.alg), temp_dir, args.emapper_dmnd, args.emapper_pfam)
    else:
        logging.info("\n--- Step 8 (Optional): Skipped eggNOG-Mapper Annotation ---")

    # --- Step 9: Finalization & Output ---
    logging.info("\n--- Step 9: Finalizing and Writing Output Files ---")
    prepare_outputs.finalize_and_write_outputs(
        tree=processed_tree,
        ogs_info=ogs_info,
        total_seqs=total_seqs,
        seqs_in_ogs=seqs_in_ogs,
        recovered_seqs=recovered_seqs,
        clean_tree_name=clean_tree_name,
        output_path=args.out_path,
        args=args
    )

    total_timer.stop()
    logging.info(f"Execution Timers: {Timer.timers}")

    # --- Step 10: (Optional) Visualization ---
    if args.open_visualization:
        logging.info("\n--- Step 10 (Optional): Opening ETE Smart View Visualization ---")
        # Robustly pass alignment path only if it exists
        alignment_path = str(args.alg) if args.alg else None
        run_smartview(processed_tree, alignment_path)
    
    