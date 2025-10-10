import os

from ogd.timer import Timer
import ogd.utils as utils

import ogd.handle_input as handle_input
from ogd.tree_setup  import run_setup
from ogd.outliers_scores import run_outliers_and_scores
from ogd.select_duplications import  run_get_main_dups
from ogd.orthologs_groups import get_all_ogs
import  ogd.pairs as pairs
from ogd.emapper_annotate import annotate_with_emapper, annot_treeprofiler
from ogd.recovery import recover_sequences
import ogd.final_annotations as final_annotations
import ogd.prepare_outputs as prepare_outputs
from ogd.run_ete4_smartview import run_smartview





def run_ogd(args):

    """
    MAIN FUNCTION, run all steps of the analysis:
        1. Load files needed for the analysis:
            Tree, taxonomy, reference tree (species tree), taxonomy counter
        2. Tree setup (Pre-analysis):  : Run some basic analysis:
            Resolve polytomies, rooting, ncbi annotation, save original species and original leaf names
        3. Outliers and Dups score:
            Detect long branches, taxonomical outliers, calculate species overlap, score1, score2 and inpalalogs_rate
        4. Detect HQ-Duplications:
            Select high quality duplication that create Monophyletic Orthologs Groups (OGs)
        5. Get OGs
        6. Optionally skip get all orthologs pairs
        7. Optionally modify OGs by recovering sequences (see RECOVERY PIPELINE)
        8. Optionally add annotations from emapper  (see EMAPPER ANNOTATION)
            8.1 Run emapper and annotate with treeprofiler
            8.2 Only annotate with treeprofiler
        9. Annotate root
        10. Flag seqs out OGs
        11. Write output files
    """

    _t = Timer('Total_time')
    _t.start()

    # 1. Load files and DBs
    mssg = f"""
    0. Load info"""
    print(mssg)

    tree = args.tree
    path_out = args.out_path

    clean_name_tree = utils.remove_file_extension(os.path.basename(tree))

    t = handle_input.load_tree_local(tree = tree, taxonomy = args.taxonomy_type, sp_delimitator = args.sp_delim)
    taxonomy_db = handle_input.load_taxonomy(taxonomy = args.taxonomy_type, user_taxonomy= args.user_taxonomy)
    reftree = handle_input.load_reftree(rtree = args.reftree, t = t, taxonomy_db = taxonomy_db)
    level2sp_mem = handle_input.load_taxonomy_counter(reftree=reftree, user_taxonomy_counter = args.user_taxonomy_counter)
    tmpdir = utils.create_tmp(path_out)



    # 2. Tree setup (Pre-analysis):  resolve polytomies, rooting, ncbi annotation, etc
    t_nw , sp_set, total_mems_in_tree, num_total_sp = run_setup(t, clean_name_tree, taxonomy_db, path_out, tmpdir, args)
   

    # 3. Outliers and Dups score functions
    t, CONTENT, total_outliers = run_outliers_and_scores(t_nw, taxonomy_db, num_total_sp, level2sp_mem, args)
    
    
    # 4. Detect HQ-Duplications
    t,  taxid_dups_og  = run_get_main_dups(t, taxonomy_db, total_mems_in_tree, args)


    # 5. Get OGs
    t, ogs_info, seqs_in_ogs = get_all_ogs(t, taxonomy_db)


    # 6. Optionally: Skip get all orthologs pairs
    if not args.skip_get_pairs:
        clean_pairs, strict_pairs = pairs.get_all_pairs(t, total_outliers)
        pairs.write_pairs_table(clean_pairs, strict_pairs, path_out, clean_name_tree)


    # 7. Optionally: Recovering sequences
    recover_seqs = set()
    seqs2recover = (total_mems_in_tree.difference(seqs_in_ogs))
    if args.run_recovery:
        if args.alg and len(seqs2recover)>0 and len(seqs_in_ogs) > 0:
            recover_seqs, ogs_info = recover_sequences(t, args.alg, ogs_info, seqs2recover, clean_name_tree, path_out, args.run_recovery)
            seqs_in_ogs.update(recover_seqs)
            

    # 8. Optionally: Run eggnog-mapper and annotate with treeprofiler
        # 8.1 Run emapper 
    if args.run_emapper:
        t = annotate_with_emapper(t, args.alg, tmpdir, args.emapper_dmnd, args.emapper_pfam)
        
        # 8.2 Only annotate with treeprofiler
    if args.path2emapper_main:
        main_table = args.path2emapper_main
        pfam_table = args.path2emapper_pfams
        t = annot_treeprofiler(t, args.alg, main_table, pfam_table, tmpdir) 
    

    # 9. Annotate root 
    final_annotations.annotate_root(ogs_info, t, clean_name_tree, total_mems_in_tree, sp_set, seqs_in_ogs, recover_seqs, taxonomy_db, args)

    
    # 10. Flag seqs out OGs
    t = final_annotations.flag_seqs_out_og(t, seqs_in_ogs, total_mems_in_tree)

    # 11. Write output files
    mssg4 = f"""
    5. Writing output files
    """
    print(mssg4)

    prepare_outputs.write_ogs_info(ogs_info, clean_name_tree, path_out)

    seq2ogs = prepare_outputs.get_seq2og(ogs_info)
    prepare_outputs.write_seq2ogs(seq2ogs, path_out,  clean_name_tree)

    t, all_props = utils.run_clean_properties(t)

    post_tree_path = utils.run_write_post_tree(t, clean_name_tree, path_out, all_props)

    if args.open_visualization:
        run_smartview(t, args.alg)


    _t.stop()
    print(Timer.timers)