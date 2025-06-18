#!/usr/bin/env python3

from ete4 import  PhyloTree, Tree
from ete4 import NCBITaxa, GTDBTaxa
from collections import defaultdict
import sys
import os
import json
import argparse
import warnings
import re
import csv
import pathlib
import tempfile
import sys
sys.setrecursionlimit(10000)

from ogd.emapper_annotate import annotate_with_emapper, annot_treeprofiler
from ogd.recovery import recover_sequences
import  ogd.pairs as pairs
import ogd.utils as utils
from ogd.tree_setup  import run_setup
from ogd.outliers_scores import run_outliers_and_scores
from ogd.select_duplications import  run_get_main_dups
from ogd.orthologs_groups import get_all_ogs
import ogd.prepare_outputs as prepare_outputs
from ogd.timer import Timer


import ete4
from ete4.smartview import Layout, explorer

from ogd.import_layouts import all_layouts

#from messy_ogs import get_messy_groups, annotate_messy_og

cwd =  str(pathlib.Path(__file__).parent.resolve())

bin_path =  cwd+'/bin'


sys.path.append(bin_path)

warnings.filterwarnings("ignore", category=RuntimeWarning)


## 1. Load initial info   ##

def create_tmp(path_out):

    """
    Create a temporary directory.

    Tries to create the temporary directory inside the results folder (`path_out`).
    If it fails (e.g., due to permissions or when running in a container),
    falls back to creating it in /tmp.

    Args:
        path_out (str): Path to the desired output directory.

    Returns:
        str: Full path to the created temporary directory, ending with a slash.
    
    """

    try:
        results_abs_path = os.path.abspath(path_out)
        tmpdir = tempfile.mkdtemp(dir=results_abs_path)

    except Exception as e:
        print(f"Warning: Failed to create tmp dir in '{path_out}', using /tmp instead. Error: {e}")
        tmpdir = tempfile.mkdtemp(prefix="ogd_", dir="/tmp")
    
    return tmpdir + '/'
    
    

    


def load_tree_local(tree=None, taxonomy = None, sp_delimitator = None):

    """
        Load Tree from file
    """
    
    mssg1 = f"""        -Load tree: {os.path.basename(tree)}"""
    print(mssg1)
  
    
    #t = PhyloTree(open(tree), parser = 0)
    t = PhyloTree(open(tree))
    t.resolve_polytomy()
   
    t.set_species_naming_function(lambda node: node.name.split(sp_delimitator)[0])

    return t


def load_taxonomy(taxonomy=None, user_taxonomy=None):

    """
        Load taxonomy from local server or from user file
        Local server:
            -Eggnog 5
            -Eggnog 6
            -GTDB v207
    """

    if taxonomy == 'NCBI':
        if user_taxonomy != None:
            taxonomy_db = NCBITaxa(user_taxonomy, memory = True)
        else:
            taxonomy_db = NCBITaxa(memory = True)

    elif taxonomy == 'GTDB':
        if user_taxonomy != None:
            taxonomy_db = GTDBTaxa(user_taxonomy, memory = True)
        else:
            taxonomy_db = GTDBTaxa()

    mssg2 = f"""        -Load taxonomy: {taxonomy_db}"""
    print(mssg2)
    

    return taxonomy_db


def load_reftree(rtree=None, t=None, taxonomy_db=None):

    """
        Get reference tree (species tree) from input tree or user provide it
    """
    
    if rtree != None:
        mssg3 = f"""        -Load reftree: from user"""
        print(mssg3)
        
        reftree = PhyloTree(open(rtree), parser = 0)
    else:
        mssg3 = f"""        -Load reftree: from gene tree"""
        print(mssg3)
        
        reftree = get_reftree(t, taxonomy_db)


    taxonomy_db.annotate_tree(reftree,  taxid_attr="name")
    
    return reftree


def get_reftree(t, taxonomy_db):

    """
        Create reference tree (species tree) if user do not provide it
    """


    taxid_list = t.get_species()

    reftree = taxonomy_db.get_topology(taxid_list)
    
    return reftree


def load_taxonomy_counter(reftree=None, user_taxonomy_counter=None):

    """
        Get number of species per each taxonomic level
        Get it from reference tree or user provide it
    """
    
    if user_taxonomy_counter:
        mssg5 = f"""        -Load taxonomy counter: from user"""
        print(mssg5)

        if isinstance(user_taxonomy_counter, dict):
            level2sp_mem = user_taxonomy_counter
        else:
            with open(user_taxonomy_counter) as levels:
                level2sp_mem = json.load(levels)
    else:
        mssg5 = f"""        -Load taxonomy counter: from gene tree"""
        print(mssg5)

        level2sp_mem = get_taxonomy_counter(reftree)

    return level2sp_mem


def get_taxonomy_counter(reftree, taxonomy = None):

    """
        Create Taxonomy counter if user de not provide it
    """

    level2sp_mem = defaultdict(set)
    for l in reftree:  
        
        if l.name.isdigit():
            lin = l.props.get('lineage')
        else:
            lin = l.props.get('named_lineage')
        
        for tax in lin:
            level2sp_mem[str(tax)].add(l.name)
    
    return level2sp_mem


##  8.Annotate root ##
def annotate_root(ogs_info, t, name_tree, total_mems_in_tree, sp_set, seqs_in_ogs, recover_seqs, taxonomy_db, args):

    """
    Add properties in root
    Add info about the root in ogs_info dict
    Need for web
    """

    # Save parameters as property

    parametrs2save = ["lineage_thr@"+str(args.lineage_thr), "best_tax_thr@"+str(args.best_tax_thr), 
                    "sp_loss_perc@"+str(args.sp_loss_perc), "inherit_out@"+str(args.inherit_out), "sp_ovlap_all@"+str(args.so_all)]
    
    if args.so_bact != None:
        parametrs2save.append("sp_ovlap_bact@"+str(args.so_bact))
    if args.so_euk != None:
        parametrs2save.append("sp_ovlap_euk@"+str(args.so_euk))
    if args.so_arq != None:
        parametrs2save.append("sp_ovlap_arq@"+str(args.so_arq))
    
    parameters_str = '|'.join(parametrs2save)

    t.add_prop('parameters', parameters_str)

    # Save general results as properties, needed for web

    seqs_out_ogs = total_mems_in_tree.difference(seqs_in_ogs)
    
    results2save = ["Tree_name@"+name_tree, "Total_seqs@"+str(len(total_mems_in_tree)), "Total_species@"+str(len(sp_set)),
                    "Seqs_in_OGs@"+str(len(seqs_in_ogs)), "Recovery_seqs@"+str(len(recover_seqs)), "Seqs_out_OGs@"+str(len(seqs_out_ogs)), 
                    "Num_OGs@"+str(len(ogs_info))]
    
    general_results_str = '|'.join(results2save)
    t.add_prop('general_result', general_results_str)


    # Save taxlevel OGs info, needed for web

    lca_all_dups = set()
    taxlev2ogs = defaultdict(set)
    taxlev2mems = defaultdict(set)

    for dup_node in t.search_nodes(node_create_og='True'):
        lca_dup = dup_node.props.get('lca_node')
        mems_dup = set(dup_node.props.get('leaves_in'))
        lca_all_dups.add(lca_dup)
        taxlev2ogs[lca_dup].add(dup_node.name)
        taxlev2mems[lca_dup].update(mems_dup)

    taxlev_list = []
    for taxlev, og in taxlev2ogs.items():
        
        #if isinstance(taxlev, str):
        taxlev = int(taxlev)

        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            sci_name = taxonomy_db.get_taxid_translator([taxlev])[taxlev]
        
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
            sci_name = taxlev
          
        sci_name_taxid = sci_name+'_'+str(taxlev)
        ogs_str = '_'.join(list(og))
        num_mems = len(taxlev2mems[taxlev])

        tax_str = '|'.join([sci_name_taxid, ogs_str,str(num_mems)])

        taxlev_list.append(tax_str)

    taxlev_str = '@'.join(taxlev_list)
    t.add_prop('taxlev2ogs', taxlev_str)

    t.add_prop("OGD_annot", True)


##  9. Flag seqs out OGs ##
def flag_seqs_out_og(t, seqs_in_ogs, total_mems_in_tree):

    """
        Flags seqs that do no belong to any OG
        Could be for taxonomic outlier or for branch lenght
    """
    seqs_out_og = total_mems_in_tree.difference(seqs_in_ogs)
    for leaf in t:
        if leaf.name in seqs_out_og:
            leaf.add_prop('seq_out_og', "true")

    t, props = utils.run_clean_properties(t)

    return t




###############################################
###############################################
###############################################



def run_app(tree, abs_path, name_tree, path_out, args):

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

    # 1. Load files and DBs
    mssg = f"""
    0. Load info"""
    print(mssg)


    clean_name_tree = utils.remove_problematic_characters(name_tree)

    t = load_tree_local(tree = tree, taxonomy = args.taxonomy_type, sp_delimitator = args.sp_delim)
    taxonomy_db = load_taxonomy(taxonomy = args.taxonomy_type, user_taxonomy= args.user_taxonomy)
    reftree = load_reftree(rtree = args.reftree, t = t, taxonomy_db = taxonomy_db)
    level2sp_mem = load_taxonomy_counter(reftree=reftree, user_taxonomy_counter = args.user_taxonomy_counter)
    tmpdir = create_tmp(path_out)

  

    # 2. Tree setup (Pre-analysis):  resolve polytomies, rooting, ncbi annotation, etc
    t_nw , sp_set, total_mems_in_tree, num_total_sp = run_setup(t, name_tree, taxonomy_db, path_out, tmpdir, args)
   

    # 3. Outliers and Dups score functions
    t, CONTENT, total_outliers = run_outliers_and_scores(t_nw, taxonomy_db, num_total_sp, level2sp_mem, args)
    
    
    # 4. Detect HQ-Duplications
    t,  taxid_dups_og  = run_get_main_dups(t, taxonomy_db, total_mems_in_tree, args)


    # 5. Get OGs
    t, ogs_info, seqs_in_ogs = get_all_ogs(t, taxonomy_db)


    # 6. Optionally skip get all orthologs pairs
    if not args.skip_get_pairs:
        clean_pairs, strict_pairs = pairs.get_all_pairs(t, total_outliers)
        pairs.write_pairs_table(clean_pairs, strict_pairs, path_out, clean_name_tree)


    # 7. Optionally modify Basal-OGs by recovering sequences
    recover_seqs = set()
    seqs2recover = (total_mems_in_tree.difference(seqs_in_ogs))
    if args.run_recovery:
        if args.alg and len(seqs2recover)>0 and len(seqs_in_ogs) > 0:
            recover_seqs, ogs_info = recover_sequences(t, args.alg, ogs_info, seqs2recover, name_tree, path_out, args.run_recovery)
            seqs_in_ogs.update(recover_seqs)
            

    # 8. Optionally add annotations from emapper
        # 8.1 Run emapper and annotate with treeprofiler
    if args.run_emapper:
        t = annotate_with_emapper(t, args.alg, tmpdir, args.emapper_data)
        
        # 8.2 Only annotate with treeprofiler
    if args.path2emapper_main:
        main_table = args.path2emapper_main
        pfam_table = args.path2emapper_pfams
        t = annot_treeprofiler(t, args.alg, main_table, pfam_table, tmpdir) 
    

    # 9. Annotate root 
    annotate_root(ogs_info, t, name_tree, total_mems_in_tree, sp_set, seqs_in_ogs, recover_seqs, taxonomy_db, args)

    
    # 10. Flag seqs out OGs
    t = flag_seqs_out_og(t, seqs_in_ogs, total_mems_in_tree)

    # 11. Write output files
    mssg4 = f"""
    5. Writing output files
    """
    print(mssg4)

    prepare_outputs.write_ogs_info(ogs_info, clean_name_tree, path_out)

    seq2ogs = prepare_outputs.get_seq2og(ogs_info)
    prepare_outputs.write_seq2ogs(seq2ogs, path_out,  clean_name_tree)

    t, all_props = utils.run_clean_properties(t)

    utils.run_write_post_tree(t, clean_name_tree, path_out, all_props)

    if args.open_visualization:
        props_popup = ['node_is_og', 'dist', 'species_losses', 'node_create_og', 'lca_node_name', 'len_leaves_in', 
        'taxid', 'sci_name', 'lineage', 'lca_dup', 'inparalogs_rate', 'ch1_name', 'dup_node_name', 'is_root', 
        'total_leaves', 'dups_up', 'ogs_up', 'common_name', 'dups_down', 'so_score_dup','ogs_down', 'score1', 
        'rank', 'dup_lineage', 'lca_node', 'len_leaves_out', 'species_losses_percentage', 'name', 'ch2_name', 
        'score2', 'sp_out', 'so_score', 'leaves_out','dup_score', 'overlap', 'evoltype_2', 'mOG', 'len_sp_in', 
        'best_tax', 'node_is_mog', 'recover_seqs', 'recover_in', 'Preferred_name', 'Preferred_name_counter',
        'eggNOG_OGs_counter', 'eggNOG_OGs', 'long_branch_outlier']
        t.explore(name = clean_name_tree, layouts = all_layouts, show_leaf_name = False , include_props = props_popup, keep_server=True , host = 'localhost', port = 5000)
       
        






def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-v', '--version', action='store_true',
                        help="show version and exit.")

    parser.add_argument('--tree', dest = 'tree', required = True,
                    help="Input tree" )

    parser.add_argument('--raw_alg', dest = 'alg', 
                    help="Input alignment. Needed for run emapper and recovery modules")

    parser.add_argument('--output_path', dest = 'out_path', required= True,
                    help="Output path")
    
    parser.add_argument('--sp_ovlap_all', default = 0.01, dest = 'so_all', type = float,
                    help='Species overlap treshold used for all nodes. ' )

    parser.add_argument('--sp_ovlap_euk', dest = 'so_euk', type = float, 
                    help="Species overlap treshold used for internal nodes that contain only eukariotas. "
                    "For all other internal nodes, the sp_ovlap_all threshold is applied.")

    parser.add_argument('--sp_ovlap_bact', dest = 'so_bact' ,  type = float, 
                        help="Species overlap treshold used for internal nodes that contain only bacterias. "
                        "For all other internal nodes, the sp_ovlap_all threshold is applied.")

    parser.add_argument('--sp_ovlap_arq', dest = 'so_arq', type = float,
                        help="Species overlap treshold used for internal nodes that contain only archeas. "
                        "For all other internal nodes, the sp_ovlap_all threshold is applied.")

    parser.add_argument('--lineage_threshold', dest = 'lineage_thr' , default = 0.05, type = float)
    parser.add_argument('--best_taxa_threshold', dest = 'best_tax_thr' , default = 0.9, type = float)
    parser.add_argument('--inherit_outliers', dest = 'inherit_out', choices = ['Yes', 'No'], default = 'Yes', type = str)
    parser.add_argument('--species_losses_perct', dest = 'sp_loss_perc' , default = 0.3, type = float)
    parser.add_argument('--rooting', choices = ['Midpoint', 'MinVar', ''])
    parser.add_argument('--taxonomy_type', choices=['NCBI', 'GTDB'], default='NCBI')
    parser.add_argument('--user_taxonomy', default= None)
    parser.add_argument('--user_taxonomy_counter', default=None)
    parser.add_argument('--reftree', default=None)
    parser.add_argument('--run_emapper', action='store_true')
    parser.add_argument('--emapper_datapath', dest = 'emapper_data')
    parser.add_argument('--run_treeprofiler_emapper_annotation', dest='path2emapper_main')
    parser.add_argument('--run_treeprofiler_emapper_pfams', dest='path2emapper_pfams')
    parser.add_argument('--run_recovery', dest = 'run_recovery',  choices= ["run-align", "skip-align"])
    parser.add_argument('--skip_get_pairs', action='store_true')
    parser.add_argument('--sp_delimitator', dest = 'sp_delim', default='-')
    parser.add_argument('--open_visualization', action='store_true')

    return parser.parse_args()


def main():

    """
        Parse argument and call run_app
    """

    args = get_args()

    init_tree = args.tree
    name_tree = os.path.basename(init_tree)
    abs_path = os.path.abspath(init_tree)
    
    _t = Timer('Total_time')
    _t.start()

    
    run_app(init_tree, abs_path, name_tree, args.out_path, args)

    _t.stop()
    print(Timer.timers)



if __name__ == "__main__":
    main()


