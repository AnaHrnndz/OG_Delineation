"""
Prepares and writes all output files for the OGD pipeline.

This module handles the final annotations on the tree, the generation of
OG tables, sequence-to-OG mappings, and the final annotated Newick tree.
"""

import argparse
import csv
import json
import logging
import re
from collections import defaultdict
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
        output_path: The directory (Path object) to write results to.
        args: The command-line arguments namespace.
    """

    
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


# --- Helper Functions (Final Annotations) ---

def _finalize_tree_annotations(
    tree: PhyloTree, ogs_info: dict, total_seqs: set,
    seqs_in_ogs: set, recovered_seqs: set, args: argparse.Namespace
):
    
    
    """Adds summary information as properties to the tree's root node."""
    params = [f"lineage_thr@{args.lineage_thr}", f"best_tax_thr@{args.best_tax_thr}",
              f"sp_loss_perc@{args.sp_loss_perc}", f"inherit_out@{args.no_inherit_outliers}",
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


# --- Worker Functions (File Writing) ---

def get_seq2og(ogs_info: Dict) -> Dict[str, Set[str]]:
    """
    Generates a dictionary mapping each sequence to a set of OG names.
    """
    seq2ogs = defaultdict(set)
    for og, info in ogs_info.items():
        og_tag = f"{og}@{info['TaxoLevel']}"
        
        for s in info['Mems']:
            seq2ogs[s].add(og_tag)
        
        for s in info['RecoverySeqs']:
            if s and s != '-': # Check for empty or placeholder strings
                seq2ogs[s].add(og_tag)

    return seq2ogs


def write_seq2ogs(seq2ogs: Dict[str, Set[str]], output_path: Path, clean_tree_name: str):
    """
    Writes the sequence-to-OG mapping to both .tsv and .jsonl files.
    """
    tsv_path = output_path / f"{clean_tree_name}.seq2ogs.tsv"
    jsonl_path = output_path / f"{clean_tree_name}.seq2ogs.jsonl"

    try:
        # Open both files at once for efficient, single-pass writing
        with open(tsv_path, 'w') as tsv_file, open(jsonl_path, 'w') as json_file:
            # Add a header to the TSV file for clarity
            tsv_file.write("Sequence\tOGs_at_TaxaLevel\n")
            
            for seq, ogs in seq2ogs.items():
                ogs_list = sorted(list(ogs)) # Sort for consistent output
                
                # Write TSV line
                tsv_file.write(f"{seq}\t{','.join(ogs_list)}\n")
                
                # Write JSONL line
                json_line = json.dumps({seq: ogs_list})
                json_file.write(json_line + "\n")
                
    except IOError as e:
        logging.error(f"Error writing seq2ogs files: {e}")


def write_ogs_info(ogs_info: Dict, clean_tree_name: str, output_path: Path):
    """
    Writes the main OG information table to a .tsv file.
    """
    
    def _join_helper(data):
        """Safely joins list/set elements or returns '-'."""
        if isinstance(data, (list, set)) and data:
            # Ensure all elements are strings before joining
            return ','.join(map(str, data))
        if isinstance(data, str) and data == '-':
            return '-'
        return '-' # Default for empty lists/sets

    tsv_path = output_path / f"{clean_tree_name}.ogs_info.tsv"

    try:
        with open(tsv_path, "w", newline='') as myfile:
            w = csv.writer(myfile, delimiter='\t')
            # Write the header
            w.writerow((
                '#OG_name', 'Lca_Dup', 'TaxoLevel', 'SciName_TaxoLevel', 'AssocNode', 
                'NumSP', 'OG_down', 'OG_up', 'NumSeqs', 'NumRecoverySeqs', 
                'Species_Outliers', 'Num_SP_Outliers', 'Inparalogs_Rate', 
                'SP_overlap_dup', 'Seqs', 'RecoverySeqs'
            ))

            # Write the data rows
            for og_name, info in ogs_info.items():
                
                og_name_extend = f"{clean_tree_name}@{og_name}"
                
                # Use the robust helper to format list-like data
                ogs_down = _join_helper(info['OG_down'])
                ogs_up = _join_helper(info['OG_up'])
                sp_outliers = _join_helper(info['Species_Outliers'])
                seqs = _join_helper(info['Mems'])
                rec_seqs = _join_helper(info['RecoverySeqs'])
                
                w.writerow((
                    og_name_extend,
                    info['Lca_Dup'],
                    info['TaxoLevel'],
                    info['SciName_TaxoLevel'],
                    info['AssocNode'],
                    info['NumSP'],
                    ogs_down,
                    ogs_up,
                    info['NumMems'],
                    info['NumRecoverySeqs'],
                    sp_outliers,
                    info['Num_SP_Outliers'],
                    info['Inparalogs_Rate'],
                    info['SP_overlap_dup'],
                    seqs,
                    rec_seqs
                ))
    except IOError as e:
        logging.error(f"Error writing OG info file: {e}")
    except KeyError as e:
        logging.error(f"Missing expected key in ogs_info dictionary: {e}. Check orthologs_groups.py")