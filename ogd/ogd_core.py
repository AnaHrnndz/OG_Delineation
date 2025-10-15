# -*- coding: utf-8 -*-
"""
ogd_core.py - Main orchestrator for the OGD pipeline.
"""

from ogd.timer import Timer

# --> Import the specific utility functions and new handler classes
import ogd.new_utils as utils
from ogd.new_utils import NCBITaxonomyHandler, GTDBTaxonomyHandler

# --> Import other pipeline modules
import ogd.handle_input as handle_input
from ogd.tree_setup import run_setup
from ogd.new_outliers import run_outliers_and_scores
from ogd.select_duplications import run_get_main_dups
from ogd.orthologs_groups import get_all_ogs
import ogd.pairs as pairs
from ogd.emapper_annotate import annotate_with_emapper, annot_treeprofiler
from ogd.recovery import recover_sequences
import ogd.final_annotations as final_annotations
import ogd.prepare_outputs as prepare_outputs
from ogd.run_ete4_smartview import run_smartview

def run_ogd(args):
    """
    Run all steps of the OGD analysis pipeline.
    """
    _t = Timer('Total_time')
    _t.start()
    
    # ----------------------------------------------------
    # 1. LOAD INPUTS AND PREPARE ENVIRONMENT
    # ----------------------------------------------------
    print("--- 1. LOADING INPUTS AND DATABASES ---")
    
    # --> Load the raw database object first
    t, raw_taxonomy_db, reftree, level2sp_mem = handle_input.load_all_input_files(args)
    
    # Initialize the correct taxonomy handler based on user arguments.
    # This handler will be passed to all downstream functions.
    if args.taxonomy_type == 'NCBI':
        tax_handler = NCBITaxonomyHandler(raw_taxonomy_db)
    elif args.taxonomy_type == 'GTDB':
        tax_handler = GTDBTaxonomyHandler(raw_taxonomy_db)
    else:
        # --> Fail early if the taxonomy type is unsupported.
        raise ValueError(f"Unsupported taxonomy type: '{args.taxonomy_type}'")

    # --> Use pathlib-aware utils functions. `args.tree` is a Path object.
    clean_name_tree = utils.remove_file_extension(args.tree.name)
    tmpdir = utils.create_tmp(args.out_path)
    
    
    # ----------------------------------------------------
    # 2. TREE SETUP (Pre-analysis)
    # ----------------------------------------------------
    print("\n--- 2. PREPARING TREE AND ANNOTATIONS ---")
    
    t, sp_set, total_mems_in_tree, num_total_sp = run_setup(
        t, tax_handler, args, tmpdir
    )
    
    # ----------------------------------------------------
    # 3. OUTLIERS AND DUPLICATION SCORES
    # ----------------------------------------------------
    print("--- 3. CALCULATING SCORES AND OUTLIERS ---")
    
    t, CONTENT, total_outliers = run_outliers_and_scores(
        t, tax_handler, num_total_sp, level2sp_mem, args
    )
    
    # ----------------------------------------------------
    # 4. DETECT HQ-DUPLICATIONS 
    # ----------------------------------------------------
    print("\n--- 4. DELINEATING ORTHOLOG GROUPS (OGs) ---")
   
    t, taxid_dups_og = run_get_main_dups(t, tax_handler, total_mems_in_tree, args)

    # ----------------------------------------------------
    # 5. GET OGS
    # ----------------------------------------------------
    
    t, ogs_info, seqs_in_ogs = get_all_ogs(t, tax_handler)

    # ----------------------------------------------------
    # 6. OPTIONAL: ORTHOLOG PAIRS
    # ----------------------------------------------------
    if not args.skip_get_pairs:
        print("--- 5. GENERATING ORTHOLOG PAIRS ---")
        clean_pairs, strict_pairs = pairs.get_all_pairs(t, total_outliers)
        pairs.write_pairs_table(clean_pairs, strict_pairs, args.out_path, clean_name_tree)

    # ----------------------------------------------------
    # 7. OPTIONAL: SEQUENCE RECOVERY
    # ----------------------------------------------------
    if args.run_recovery:
        print("--- 6. RUNNING SEQUENCE RECOVERY PIPELINE ---")
        recover_seqs, ogs_info = recover_sequences(
            t, ogs_info, total_mems_in_tree, seqs_in_ogs, clean_name_tree, args
        )
        seqs_in_ogs.update(recover_seqs)
    else:
        recover_seqs = set() 

    # ----------------------------------------------------
    # 8. OPTIONAL: EMAPPER ANNOTATION
    # ----------------------------------------------------
    if args.run_emapper or args.path2emapper_main:
        print("--- 7. RUNNING EMAPPER ANNOTATION ---")
        if args.run_emapper:
            t = annotate_with_emapper(t, tmpdir, args)
        if args.path2emapper_main:
            t = annot_treeprofiler(t, tmpdir, args)

    # ----------------------------------------------------
    # 9. & 10. FINAL ANNOTATIONS
    # ----------------------------------------------------
    print("--- 8. FINALIZING TREE ANNOTATIONS ---")
    
    final_annotations.annotate_root(
        t, ogs_info, clean_name_tree, total_mems_in_tree, sp_set, 
        seqs_in_ogs, recover_seqs, tax_handler, args
    )
    t = final_annotations.flag_seqs_out_og(t, seqs_in_ogs, total_mems_in_tree)

    # ----------------------------------------------------
    # 11. WRITE OUTPUTS AND CLEANUP
    # ----------------------------------------------------
    print("--- 9. WRITING FINAL OUTPUTS AND CLEANUP ---")

    prepare_outputs.write_ogs_info(ogs_info, clean_name_tree, args.out_path)
    seq2ogs = prepare_outputs.get_seq2og(ogs_info)
    prepare_outputs.write_seq2ogs(seq2ogs, args.out_path, clean_name_tree)

    post_tree_path = prepare_outputs.clean_and_write_tree(t, clean_name_tree, args.out_path)

    if args.open_visualization:
        run_smartview(t, args.alg)

    # TODO: add a cleanup function for the temporary directory
    # import shutil
    # shutil.rmtree(tmpdir)

    _t.stop()
    print("\n" + str(Timer.timers))
    print(f"Pipeline complete. All results are located in: {args.out_path.resolve()}")