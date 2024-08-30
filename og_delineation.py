#!/usr/bin/env python3

from ete4 import  PhyloTree, SeqGroup, Tree
from ete4 import NCBITaxa, GTDBTaxa
from collections import Counter, OrderedDict, defaultdict
import sys
import os
import json
import argparse
import warnings
import re
import glob
import subprocess
import csv
import pathlib

cwd =  str(pathlib.Path(__file__).parent.resolve())
sys.path.append(cwd+'/ogd')
from emapper_annotate import annotate_with_emapper, annot_treeprofiler
from recovery import recover_sequences
import  pairs
import utils
from tree_setup  import run_setup
from outliers_scores import run_outliers_and_scores
from select_duplications import  run_get_main_dups
from orthologs_groups import get_ogs
from timer import Timer
from messy_ogs import get_messy_groups, annotate_messy_og

cwd =  str(pathlib.Path(__file__).parent.resolve())

bin_path =  cwd+'/bin'
data_path = cwd+'/data'
sys.path.append(bin_path)

warnings.filterwarnings("ignore", category=RuntimeWarning)


## 1. Load initial info   ##

def create_tmp(path_out):
    path = os.path.join(path_out, 'tmp_dir') 
  
# Create the directory 
# 'Nikhil' 
    try: 
        os.makedirs(path, exist_ok = True) 
        print("Directory '%s' created successfully" % path) 
    except OSError as error: 
        print("Directory '%s' can not be created" % path)

    return path 

def load_tree_local(tree=None, taxonomy = None, sp_delimitator = None):

    """
        Load Tree from file
    """

    print(' -Load tree:', os.path.basename(tree))
    
    t = PhyloTree(open(tree), parser = 0)

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

    print(' -Load taxonomy:  ')

    return taxonomy_db


def load_reftree(rtree=None, t=None, taxonomy_db=None):

    """
        Get reference tree (species tree) from input tree or user provide it
    """

    print(' -Load reftree:   ')
    if rtree:
            reftree = PhyloTree(open(rtree))
    else:
        reftree = get_reftree(t, taxonomy_db)


    taxonomy_db.annotate_tree(reftree,  taxid_attr="name")
    
    return reftree


def get_reftree(t, taxonomy_db):

    """
        Create reference tree (species tree) if user do not provide it
    """

    print('\t**Create referente species tree from gene tree')

    taxid_list = t.get_species()

    reftree = taxonomy_db.get_topology(taxid_list)

    
    return reftree


def load_taxonomy_counter(reftree=None, user_taxonomy_counter=None):

    """
        Get number of species per each taxonomic level
        Get it from reference tree or user provide it
    """

    print(' -Load taxonomy counter:  ')
    if user_taxonomy_counter:
        if isinstance(user_taxonomy_counter, dict):
            level2sp_mem = user_taxonomy_counter
        else:
            with open(user_taxonomy_counter) as levels:
                level2sp_mem = json.load(levels)
    else:
        level2sp_mem = get_taxonomy_counter(reftree)

    
    return level2sp_mem


def get_taxonomy_counter(reftree, taxonomy = None):

    """
        Create Taxonomy counter if user de not provide it
    """

    print('\t**Create taxonomy counter  from gene tree')

    level2sp_mem = defaultdict(set)
    for l in reftree:  
        
        if l.name.isdigit():
            lin = l.props.get('lineage')
        else:
            lin = l.props.get('named_lineage')
        
        for tax in lin:
            level2sp_mem[str(tax)].add(l.name)
    
    return level2sp_mem



##  6. Add info about Core-OGs up and down  ##

def add_nodes_up_down(t):

    """
        Traverse the tree, for each node,  detect ogs that are up or below
    """

    for node in t.traverse('preorder'):
        if node.props.get('node_is_og'):
            node_name = node.name

            #Detect OGs below node
            ogs_down = set()
            dups_down = list()
            for n in node.search_nodes(node_is_og="True"):
                if node_name != n.name:
                    ogs_down.add(n.name)
                    dups_down.append(n.up.name)

            if len(ogs_down) > 0:
                node.add_prop('ogs_down',ogs_down)
                node.add_prop('dups_down', dups_down)
            else:
                node.add_prop('ogs_down','-')
                node.add_prop('dups_down', '-')

            #Detect OGs up
            ogs_up = set()
            dups_up = list()
            ogs_up, dups_up = utils.check_nodes_up(node)

            if len(ogs_up) > 0:
                node.add_prop('ogs_up', ogs_up)
                node.add_prop('dups_up', dups_up)
            else:
                node.add_prop('ogs_up', '-')
                node.add_prop('dups_up', '-')

    ogs_down = set()
    dups_down = list()
    for n in t.search_nodes(node_is_og="True"):
        if n.name != t.name:
            ogs_down.add(n.name)
            dups_down.append(n.up.name)

    if len(ogs_down) > 0:
        t.add_prop('ogs_down',ogs_down)
        t.add_prop('dups_down', dups_down)
    else:
        t.add_prop('ogs_down','-')
        t.add_prop('dups_down', '-')


    return t


# 7. Annotate base-OGs
def annot_ogs(t, base_ogs, taxonomy_db):

    ##OG_name   TaxoLevel   AssocNode  lensp_in_OG OG_up   OG_down num_OG_mems    members

    annot_base_ogs = defaultdict(dict)
    total_mems_in_ogs = set()
    for taxa, ogs in base_ogs.items():
        for og_name, og_info in ogs.items():
            
            if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
                sci_name_taxa =  taxonomy_db.get_taxid_translator([taxa])[int(taxa)]
        
            elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
                sci_name_taxa = taxa
            
            
            ogs_down = '|'.join(list((t[og_info[0]].props.get('ogs_down','-'))))
            ogs_up = '|'.join(list((t[og_info[0]].props.get('ogs_up', '-'))))
            recover_seqs = list(t[og_info[0]].props.get('recovery_seqs', list()))



            mems = og_info[1].split('|')
            sp_in_og = set()
            for l in mems:
                sp_in_og.add(l.split('.')[0])

            if not recover_seqs:
                len_recover_seqs = '0'
                recover_seqs.append('-')
            else:
                len_recover_seqs = str(len(recover_seqs))


            annot_base_ogs[og_name]['TaxoLevel'] = taxa
            annot_base_ogs[og_name]['SciName_TaxoLevel'] = sci_name_taxa.replace(' ', '_')
            annot_base_ogs[og_name]['AssocNode'] = og_info[0]
            annot_base_ogs[og_name]['NumSP'] = len(sp_in_og)
            annot_base_ogs[og_name]['OG_down'] = ogs_down
            annot_base_ogs[og_name]['OG_up'] = ogs_up
            annot_base_ogs[og_name]['NumMems'] = len(mems)
            annot_base_ogs[og_name]['Mems'] = '|'.join(mems)
            annot_base_ogs[og_name]['RecoverySeqs'] = '|'.join(recover_seqs)
            annot_base_ogs[og_name]['NumRecoverySeqs'] = len_recover_seqs

            total_mems_in_ogs.update(set(mems))

    return annot_base_ogs, total_mems_in_ogs


##  8.Annotate root ##
def annotate_root(base_ogs_annot, t, name_tree, total_mems_in_tree, sp_set, total_mems_in_ogs, recovery_seqs, taxonomy_db, args):

    """
    Add properties in root
    Add info about the root in ogs_info dict
    """

    # Save parameters as property

    parameters_str = '|'.join(["outliers_node@"+str(args.outliers_node), "outliers_reftree@"+str(args.outliers_reftree), "sp_loss_perc@"+str(args.sp_loss_perc),
                                "so_cell_org@"+str(args.so_cell_org), "so_arq@"+str(args.so_arq), "so_bact@"+str(args.so_bact), "so_euk@"+str(args.so_euk), "inherit_out@"+str(args.inherit_out)])

    t.add_prop('parameters', parameters_str)

    # Save general results as properties, needed for web

    seqs_out_ogs = total_mems_in_tree.difference(total_mems_in_ogs)

    general_results_str = '|'.join(["Tree_name@"+name_tree, "Total_seqs@"+str(len(total_mems_in_tree)), "Total_species@"+str(len(sp_set)),
                                    "Seqs_in_OGs@"+str(len(total_mems_in_ogs)), "Recovery_seqs@"+str(len(recovery_seqs)), "Seqs_out_OGs@"+str(len(seqs_out_ogs)), "Num_OGs@"+str(len(base_ogs_annot)) ])

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
        
        if isinstance(taxlev, str):
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
def flag_seqs_out_og(t, total_mems_in_ogs, total_mems_in_tree):

    """
        Flags seqs that do no belong to any OG
    """
    seqs_out_og = total_mems_in_tree.difference(total_mems_in_ogs)
    for leaf in t:
        if leaf.name in seqs_out_og:
            leaf.add_prop('seq_out_og', "true")

    t, props = utils.run_clean_properties(t)

    return t



#####   FUNCIONTS TO PREPARE OUTPUTS FILES (newick, etc)    ####




def get_seq2og_v2(base_ogs_annot, messy_ogs):

    seq2ogs = defaultdict(set)
    for og, info in base_ogs_annot.items():
        for s in info['Mems'].split('|'):
            seq2ogs[s].add(og+'@'+info['TaxoLevel'])
        for s in info['RecoverySeqs'].split('|'):
            if s != '-':
                seq2ogs[s].add(og+'@'+info['TaxoLevel'])

    for mog, info in messy_ogs.items():
        for s in info['Mems']:
            seq2ogs[s].add(mog+'@'+info['TaxoLevel'])
        for s in info['RecoverySeqs']:
            if s != '-':
                seq2ogs[s].add(og+'@'+info['TaxoLevel'])

    return seq2ogs


def write_seq2ogs(seq2ogs, path, name_tree):

    """
        Write a table with seq2ogs info
    """

    name_fam = name_tree.split('.')[0]
    seq2ogs_out = open(path+'/'+name_fam+'.seq2ogs.tsv', 'w')
    for seq, ogs in seq2ogs.items():
        seq2ogs_out.write(seq+'\t'+'|'.join(list(ogs))+'\n')

    seq2ogs_out.close()


def write_ogs_info(base_ogs_annot, messy_ogs, name_tree, path):

    name_fam = name_tree.split('.',1)[0]
    name_out =  path+'/'+name_fam+'.ogs_info.tsv'


    with open(name_out, "w",  newline='') as myfile:
        w = csv.writer(myfile, delimiter='\t')
        w.writerow(('#OG_name', 'TaxoLevel', 'SciName_TaxoLevel', 'AssocNode',  'NumSP', 'OG_down', 'OG_up', 'NumSeqs', 'NumRecoverySeqs','Seqs', 'RecoverySeqs'))

        for og_name in base_ogs_annot.keys():

            w.writerow((
                og_name,    #1
                base_ogs_annot[og_name]['TaxoLevel'], #2
                base_ogs_annot[og_name]['SciName_TaxoLevel'], #3
                base_ogs_annot[og_name]['AssocNode'],   #4
                base_ogs_annot[og_name]['NumSP'],   #5
                base_ogs_annot[og_name]['OG_down'], #6
                base_ogs_annot[og_name]['OG_up'] ,  #7
                base_ogs_annot[og_name]['NumMems'], #8
                base_ogs_annot[og_name]['NumRecoverySeqs'], #9
                base_ogs_annot[og_name]['Mems'], #10
                base_ogs_annot[og_name]['RecoverySeqs'] #11
            ))

        for mog_name in messy_ogs.keys():
            w.writerow((
                mog_name, #1
                messy_ogs[mog_name]['TaxoLevel'],  #2
                messy_ogs[mog_name]['SciName_TaxoLevel'],  #3
                messy_ogs[mog_name]['AssocNode'], #4
                messy_ogs[mog_name]['NumSP'],  #5
                '|'.join(list(messy_ogs[mog_name]['OG_down'])), #6
                '|'.join(list(messy_ogs[mog_name]['OG_up'])), #7
                messy_ogs[mog_name]['NumMems'], #8
                str(messy_ogs[mog_name]['NumRecoverySeqs']) , #9
                '|'.join(list(messy_ogs[mog_name]['Mems'])), #10
                '|'.join((messy_ogs[mog_name]['RecoverySeqs']))#11
            ))





###############################################
###############################################
###############################################



def run_app(tree, abs_path, name_tree, reftree, user_counter, user_taxo, taxonomy_type, path_out, args):

    """
    MAIN FUNCTION, run all steps of the analysis:
        1. Load files needed for the analysis:
                Tree, taxonomy, reference tree (species tree), taxonomy counter
        2. Tree setup (Pre-analysis):  : Run some basic analysis:
                Resolve polytomies, rooting, ncbi annotation, save original species and original leaf names
        3. Outliers and Dups score:
                Detect long branches, taxonomical outliers, calculate species overlap, score1, score2 and inpalalogs_rate
        4. Detect HQ-Duplications:
                Select high quality duplication that create Basal Orthologs Groups (OGs)
        5. Get OGs for all taxonomical levels
        6. Add info about Basal-OGs up and down:
                For each Basal-OG detect upper and below Basal-OGs (Basal-OGs have nested structure)
        7. Annotate Basal-OGs with taxonomy, etc
        8. Write a table with the results for the Basal-OGs
        9. Optionally modify Basal-OGs by recovering sequences (see RECOVERY PIPELINE)
        10. Optionally add annotations from emapper (see EMAPPER ANNOTATION)
        11. Optionally  Get all orthologs pairs
        12. Annotate root:
            Root is a special node that has to be annotated diffent
        13. Flag seqs out OGs
        14. Write output files:
            annot_tree
            seq2ogs

    """

    # 1. Load files and DBs
    print('\n0. Load info')
    t = load_tree_local(tree = tree, taxonomy = taxonomy_type, sp_delimitator = args.sp_delim)
    taxonomy_db = load_taxonomy(taxonomy = taxonomy_type, user_taxonomy= user_taxo)
    reftree = load_reftree(rtree = reftree, t = t, taxonomy_db = taxonomy_db)
    level2sp_mem = load_taxonomy_counter(reftree=reftree, user_taxonomy_counter = user_counter)
    tmp_path = create_tmp(path_out)
    


    # 2. Tree setup (Pre-analysis):  resolve polytomies, rooting, ncbi annotation, etc
    t_nw , sp_set, total_mems_in_tree, NUM_TOTAL_SP, user_props = run_setup(t, name_tree, taxonomy_db, args.rooting, path_out, abs_path, args.sp_delim)

    # 3. Outliers and Dups score functions
    t, CONTENT = run_outliers_and_scores(t_nw, taxonomy_db, NUM_TOTAL_SP, level2sp_mem, args)

    # 4. Detect HQ-Duplications
    #t, total_mems_in_ogs, taxid_dups_og  = run_get_main_dups(t, taxonomy_db, total_mems_in_tree, args)
    t,  taxid_dups_og  = run_get_main_dups(t, taxonomy_db, total_mems_in_tree, args)


    # 5. Get Basal-OGs for all taxonomical levels
    base_ogs = get_ogs(t, level2sp_mem,taxonomy_db)

    # 5. Add info about OGs up and down
    t = add_nodes_up_down(t)

    # 6. Annotate Basal-OGs
    #base_ogs_annot = annot_ogs(t, base_ogs, taxonomy_db)
    base_ogs_annot, total_mems_in_ogs = annot_ogs(t, base_ogs, taxonomy_db)

    # 7. Get messy OGs
    t, messy_ogs, seqs_in_messy_ogs = get_messy_groups(t, taxonomy_db)
    total_mems_in_ogs.update(seqs_in_messy_ogs)


    # 8. Optionally modify Basal-OGs by recovering sequences
    recovery_seqs = set()
    best_match = defaultdict()
    if args.run_recovery:
        if args.alg and len(total_mems_in_tree.difference(total_mems_in_ogs)) > 0 and len(total_mems_in_ogs) > 0:
            recovery_seqs, best_match = recover_sequences(tree, t, args.alg, total_mems_in_tree, total_mems_in_ogs, name_tree, path_out, base_ogs_annot, args.mode)
            #base_ogs_annot = annot_ogs(t, base_ogs, taxonomy_db)
            base_ogs_annot, total_mems_in_ogs = annot_ogs(t, base_ogs, taxonomy_db)


    else:
        recovery_seqs = set()
        best_match = defaultdict()


    # 9. Optionally add annotations from emapper
        # 9.1 Run emapper and annotate with treeprofiler
    if args.run_emapper:
        t = annotate_with_emapper(t, args.alg, tmp_path)
        
        # 9.2 Only annotate with treeprofiler
    if args.path2emapper_main:
        main_table = args.path2emapper_main
        pfam_table = args.path2emapper_pfams
        t = annot_treeprofiler(t, args.alg, main_table, pfam_table, tmp_path)
        

    # 10. Optionally skip get all orthologs pairs
    if not args.skip_get_pairs:
        clean_pairs = pairs.get_all_pairs(CONTENT, total_mems_in_ogs)
        pairs.write_pairs_table(clean_pairs, path_out, name_tree)


    # 11. Annotate root & messy_ogs
    annotate_root(base_ogs_annot, t, name_tree, total_mems_in_tree, sp_set, total_mems_in_ogs, recovery_seqs, taxonomy_db, args)
    messy_ogs_annot = annotate_messy_og(t, messy_ogs)

    # 12. Flag seqs out OGs
    total_mems_in_ogs.update(recovery_seqs)
    t = flag_seqs_out_og(t, total_mems_in_ogs, total_mems_in_tree)

    # 13. Write output files
    print('\n4. Writing outputfiles\n')

    write_ogs_info(base_ogs_annot, messy_ogs_annot, name_tree, path_out)

    seq2ogs = get_seq2og_v2(base_ogs_annot, messy_ogs_annot)
    write_seq2ogs(seq2ogs, path_out,  name_tree)

    t, all_props = utils.run_clean_properties(t)

    utils.run_write_post_tree(t, name_tree, path_out, all_props)


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', dest = 'tree', required = True)
    parser.add_argument('--raw_alg', dest = 'alg')
    parser.add_argument('--output_path', dest = 'out_path', required= True)
    parser.add_argument('--mode', dest='mode', choices= ["fast", "regular"], default = "regular", type = str)
    parser.add_argument('--species_overlap_cell_org', dest = 'so_cell_org' , default = 0.1, type = float)
    parser.add_argument('--species_overlap_euk', dest = 'so_euk' , default = 0.1, type = float)
    parser.add_argument('--species_overlap_bact', dest = 'so_bact' , default = 0.1, type = float)
    parser.add_argument('--species_overlap_arq', dest = 'so_arq' , default = 0.1, type = float)
    parser.add_argument('--outliers_in_node', dest = 'outliers_node' , default = 0.05, type = float)
    parser.add_argument('--outliers_in_reftree', dest = 'outliers_reftree' , default = 0.01, type = float)
    parser.add_argument('--inherit_outliers', dest = 'inherit_out', choices = ['Yes', 'No'], default = 'Yes', type = str)
    parser.add_argument('--species_losses_perct', dest = 'sp_loss_perc' , default = 0.7, type = float)
    parser.add_argument('--rooting', choices = ['Midpoint', 'MinVar', ''])
    parser.add_argument('--taxonomy_type', choices=['NCBI', 'GTDB'], default='NCBI')
    parser.add_argument('--user_taxonomy')
    parser.add_argument('--user_taxonomy_counter')
    parser.add_argument('--reftree')
    parser.add_argument('--run_emapper', action='store_true')
    parser.add_argument('--run_treeprofiler_emapper_annotation', dest='path2emapper_main')
    parser.add_argument('--run_treeprofiler_emapper_pfams', dest='path2emapper_pfams')
    parser.add_argument('--run_recovery', action='store_true')
    parser.add_argument('--skip_get_pairs', action='store_true')
    parser.add_argument('--sp_delimitator', dest = 'sp_delim', default='-')



    return parser.parse_args()


def main():

    """
        Parse argument and call run_app
    """

    args = get_args()

    init_tree = args.tree
    name_tree = os.path.basename(init_tree)
    abs_path = os.path.abspath(init_tree)
    so_cell_org = args.so_cell_org
    so_euk = args.so_euk
    so_bact = args.so_bact
    so_arq = args.so_arq
    outliers_node = args.outliers_node
    outliers_reftree = args.outliers_reftree
    sp_loss_perc = args.sp_loss_perc
    taxonomy_type = args.taxonomy_type
    rooting = args.rooting
    path_out = args.out_path
    inherit_out = args.inherit_out
    alg = args.alg
    mode = args.mode

    if args.out_path:
        path_out = args.out_path

    if args.user_taxonomy:
        user_taxo = args.user_taxonomy
    else:
        user_taxo = None

    if args.reftree:
        rtree = args.reftree
    else:
        rtree = None

    if args.user_taxonomy:
        user_counter = args.user_taxonomy_counter
    else:
        user_counter = None


    _t = Timer('Total_time')
    _t.start()

    run_app(init_tree, abs_path, name_tree, rtree, user_counter, user_taxo, taxonomy_type, path_out, args)

    _t.stop()
    print(Timer.timers)



if __name__ == "__main__":
    main()


