import re
import csv
import pathlib
import json
from collections import defaultdict




"""
Prepares and writes all output files for the OGD pipeline.

This module handles the final annotations on the tree, the generation of
OG tables, sequence-to-OG mappings, and the final annotated Newick tree.
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, Set

from ete4 import PhyloTree

import ogd.utils as utils

# --- Main Orchestrator for Outputs ---

def finalize_and_write_outputs(
    tree: PhyloTree,
    ogs_info: Dict,
    total_seqs: Set[str],
    seqs_in_ogs: Set[str],
    recovered_seqs: Set[str],
    clean_tree_name: str,
    output_path: Path,
    args: argparse.Namespace
):
    """
    Orchestrates all finalization and output writing steps.

    Args:
        tree: The final, processed PhyloTree object.
        ogs_info: Dictionary containing information about each OG.
        total_seqs: Set of all sequence names in the original tree.
        seqs_in_ogs: Set of sequence names assigned to an OG.
        recovered_seqs: Set of sequence names recovered by the recovery step.
        clean_tree_name: The base name for output files.
        output_path: The directory to write results to.
        args: The command-line arguments namespace.
    """
    logging.info("--- Step 9: Finalizing Annotations and Writing Outputs ---")

    # 1. Add final summary annotations to the tree root
    _finalize_tree_annotations(tree, ogs_info, total_seqs, seqs_in_ogs, recovered_seqs, args)

    # 2. Flag all sequences that were not assigned to any OG
    _flag_unassigned_sequences(tree, seqs_in_ogs, total_seqs)

    # 3. Write OG information and sequence mappings to tables
    write_ogs_info(ogs_info, clean_tree_name, output_path)
    seq2ogs = get_seq2og(ogs_info)
    write_seq2ogs(seq2ogs, output_path, clean_tree_name)

    # 4. Sanitize and write the final annotated tree
    sanitized_tree, all_props = utils.sanitize_tree_properties(tree)
    utils.write_annotated_tree(sanitized_tree, output_path, clean_tree_name, all_props)


# --- Helper Functions (moved from ogd_core.py) ---

def _finalize_tree_annotations(
    tree: PhyloTree, ogs_info: dict, total_seqs: set,
    seqs_in_ogs: set, recovered_seqs: set, args: argparse.Namespace
):
    """Adds summary information as properties to the tree's root node."""
    params = [f"lineage_thr@{args.lineage_thr}", f"best_tax_thr@{args.best_tax_thr}",
              f"sp_loss_perc@{args.sp_loss_perc}", f"inherit_out@{args.inherit_outliers}",
              f"sp_ovlap_all@{args.so_all}"]
    if args.so_bact is not None: params.append(f"sp_ovlap_bact@{args.so_bact}")
    if args.so_euk is not None: params.append(f"sp_ovlap_euk@{args.so_euk}")
    if args.so_arq is not None: params.append(f"sp_ovlap_arq@{args.so_arq}")
    tree.add_prop('parameters', '|'.join(params))

    seqs_out_ogs = total_seqs - seqs_in_ogs
    results = [f"Tree_name@{args.tree.name}", f"Total_seqs@{len(total_seqs)}",
               f"Total_species@{len(tree.get_species())}", f"Seqs_in_OGs@{len(seqs_in_ogs)}",
               f"Recovered_seqs@{len(recovered_seqs)}", f"Seqs_out_OGs@{len(seqs_out_ogs)}",
               f"Num_OGs@{len(ogs_info)}"]
    tree.add_prop('general_result', '|'.join(results))
    tree.add_prop("OGD_annot", True)

def _flag_unassigned_sequences(tree: PhyloTree, seqs_in_ogs: set, total_seqs: set):
    """Flags leaf nodes that were not assigned to any OG."""
    unassigned_seqs = total_seqs - seqs_in_ogs
    for leaf in tree:
        if leaf.name in unassigned_seqs:
            leaf.add_prop('seq_out_og', "true")
    logging.info(f"Flagged {len(unassigned_seqs)} sequences as not belonging to any OG.")




#####   FUNCIONTS TO PREPARE OUTPUTS FILES (newick, etc)    ####
def get_seq2og(ogs_info):

    seq2ogs = defaultdict(set)
    for og, info in ogs_info.items():
        for s in info['Mems']:
            seq2ogs[s].add(og+'@'+str(info['TaxoLevel']))
        for s in info['RecoverySeqs']:
            if s != '-':
                seq2ogs[s].add(og+'@'+str(info['TaxoLevel']))

    return seq2ogs



def write_seq2ogs(seq2ogs, path, clean_name_tree):

    """
        Write a table with seq2ogs info
    """

    
    #clean_name_tree = clean_name_tree.split('.')[0]
    seq2og_list =[]
    seq2ogs_out = open(path+'/'+clean_name_tree+'.seq2ogs.tsv', 'w')
    for seq, ogs in seq2ogs.items():
        seq2og_list.append({seq:list(ogs)})
        seq2ogs_out.write(seq+'\t'+','.join(list(ogs))+'\n')

    seq2ogs_out.close()

    seq2ogs_json = (path+'/'+clean_name_tree+'.seq2ogs.jsonl')
    with open(seq2ogs_json, "w") as file:
        for seq in seq2og_list:
            file.write(json.dumps(seq) + "\n")
    





def write_ogs_info(ogs_info, clean_name_tree, path):

    
    name_out =  path+'/'+clean_name_tree+'.ogs_info.tsv'


    with open(name_out, "w",  newline='') as myfile:
        w = csv.writer(myfile, delimiter='\t')
        w.writerow(('#OG_name', 'Lca_Dup','TaxoLevel', 'SciName_TaxoLevel', 'AssocNode',  'NumSP', 'OG_down', 'OG_up', 
        'NumSeqs', 'NumRecoverySeqs',  'Species_Outliers', 'Num_SP_Outliers', 'Inparalogs_Rate', 'SP_overlap_dup',
        'Seqs', 'RecoverySeqs'))

        for og_name in ogs_info.keys():
            
            og_name_extend = clean_name_tree+'@'+og_name
            ogs_down = ','.join(list(ogs_info[og_name]['OG_down']))
            ogs_up = ','.join(list(ogs_info[og_name]['OG_up']))
            sp_outliers = ','.join(ogs_info[og_name]['Species_Outliers'])
            seqs = ','.join(ogs_info[og_name]['Mems'])
            if len(ogs_info[og_name]['RecoverySeqs']) == 0:
                rec_seqs = '-'
            else:
                rec_seqs = ','.join(ogs_info[og_name]['RecoverySeqs'])
            
            w.writerow((
                og_name_extend,    #1
                ogs_info[og_name]['Lca_Dup'],
                ogs_info[og_name]['TaxoLevel'], #2
                ogs_info[og_name]['SciName_TaxoLevel'], #3
                ogs_info[og_name]['AssocNode'],   #4
                ogs_info[og_name]['NumSP'],   #5
                ogs_down, #6
                ogs_up ,  #7
                ogs_info[og_name]['NumMems'], #8
                ogs_info[og_name]['NumRecoverySeqs'], #9
                sp_outliers, #11
                ogs_info[og_name]['Num_SP_Outliers'], #12
                ogs_info[og_name]['Inparalogs_Rate'], #13
                ogs_info[og_name]['SP_overlap_dup'], #14
                seqs, #16
                rec_seqs #16
            ))
