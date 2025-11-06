import logging
from collections import defaultdict
from typing import Any, Dict, Set, Tuple
import re

from ete4 import PhyloTree

import ogd.utils as utils
from ogd.utils import TaxonomyDB


# --- Main Orchestrator Function ---

def get_all_ogs(tree: PhyloTree, tax_db: TaxonomyDB) -> Tuple[PhyloTree, Dict, Set[str]]:
    """
    Finds and annotates all OGs from the processed tree.
    
    This function orchestrates the collection of:
    1. Monophyletic OGs (from 'monophyletic_og=True' nodes)
    2. Paraphyletic OGs (members not captured by monophyletic OGs)
    3. Root-level OGs (ancestral groups and root-as-OG)
    
    It then calculates the hierarchy and returns the final OG dictionary.

    Args:
        tree: The fully annotated PhyloTree.
        tax_db: The initialized taxonomy database.

    Returns:
        A tuple containing:
        - The tree (potentially with new 'og_name' annotations).
        - The comprehensive dictionary of all OGs (ogs_info).
        - A set of all sequence names included in any OG.
    """
    
    ogs_info = {}

    # 1. Get Monophyletic OGs
    logging.info("Identifying monophyletic OGs...")
    _add_monophyletic_ogs(tree, tax_db, ogs_info)
    
    # 2. Find and Annotate Paraphyletic OGs
    # This process is split in two: find candidates, then annotate them.
    logging.info("Identifying paraphyletic OGs...")
    pog_candidates = _find_paraphyletic_ogs(tree)
    _annot_paraphyletic_ogs(tree, tax_db, ogs_info, pog_candidates)
    
    # 3. Process OGs related to the root
    logging.info("Processing root-level OGs...")
    _process_root_ogs(tree, tax_db, ogs_info)
    
    # 4. Calculate hierarchical relationships (up/down)
    logging.info("Calculating OG hierarchy...")
    _calculate_hierarchy(ogs_info, tax_db)
    
    # 5. Get all unique sequences in OGs
    seqs_in_ogs = set()
    for og, info in ogs_info.items():
        seqs_in_ogs.update(info['Mems'])
        info['NumRecoverySeqs'] = 0
        info['RecoverySeqs'] = []
    
    logging.info(f"Defined {len(ogs_info)} total OGs, containing {len(seqs_in_ogs)} sequences.")

    
    return tree, ogs_info, seqs_in_ogs

# --- 1. Monophyletic OGs ---

def _add_monophyletic_ogs(tree: PhyloTree, tax_db: TaxonomyDB, ogs_info: Dict):
    """
    Finds all nodes marked as 'monophyletic_og' and adds them to ogs_info.
    """
    
    for node in tree.search_nodes(monophyletic_og=True):
        if node.is_root:
            continue  # Root node is handled separately by _process_root_ogs

        og_name = f"{node.props.get('lca_dup')}|{node.name}"
        node.add_prop('og_name', og_name)
        
        lca_og = str(node.props.get('lca_node'))
        sci_name_lca = utils.get_scientific_name(tax_db, lca_og)
        
        values = [
            og_name,
            lca_og,
            sci_name_lca,
            node.name,
            node.props.get('sp_in', set()),
            node.props.get('leaves_in', set()),
            node.props.get('lca_dup'),
            node.props.get('sp_out', set()),
            node.props.get('inparalogs_rate'),
            node.props.get('so_score_dup')
        ]
        
        _add_entry_to_ogs_info(ogs_info, values)



# --- 2. Paraphyletic OGs ---

def _find_paraphyletic_ogs(tree: PhyloTree) -> Dict[str, Any]:
    """
    Finds nodes and member lists that constitute paraphyletic OGs.
    
    A POG is formed by leaves under a duplication child that is *not* a
    monophyletic OG, excluding any leaves that belong to *other* nested OGs.
    """
    
    paraphyletic_ogs = {}
    
    for node in tree.traverse():
        # We only check children of duplication nodes
        if node.props.get('evoltype_2') != 'D':
            continue

        lca_target = node.props.get('lca_node')
            
        for child in node.children:
            # If the child is already a monophyletic OG, skip it
            if child.props.get('monophyletic_og') is True:
                continue

            all_mems = child.props.get('leaves_in', set())
            mems_to_remove = set()
            
            # Find all nested OGs to subtract their members
            # OGs defined from a duplication at the *same* tax level
            nested_ogs = child.search_nodes(monophyletic_og=True, lca_dup=lca_target)
            # Other duplication nodes at the *same* tax level
            nested_dups = child.search_nodes(evoltype_2='D', lca_node=lca_target)

            for n in nested_ogs:
                mems_to_remove.update(n.props.get('leaves_in', set()))
            for n in nested_dups:
                mems_to_remove.update(n.props.get('leaves_in', set()))

            
            # The remaining members form the paraphyletic group
            pog_members = set(all_mems).difference(mems_to_remove)
            
            if len(pog_members) > 1:
                og_name = f"{child.props.get('lca_node')}|{child.name}!" # '!' denotes POG
                child.add_prop('paraphyletic_og', True)
                child.add_prop('og_name', og_name)
                paraphyletic_ogs[og_name] = [child, list(pog_members)]
                
    return paraphyletic_ogs

def _annot_paraphyletic_ogs(
    tree: PhyloTree, tax_db: TaxonomyDB, ogs_info: Dict, pog_candidates: Dict):
    """
    Validates members of POG candidates and adds them to the ogs_info dict.
    """
    
    for pog_name, (pog_node, member_names) in pog_candidates.items():
        
        lca_pog = pog_node.props.get('lca_node')
        sci_name_lca = utils.get_scientific_name(tax_db, lca_pog)

        retained_species = set()
        retained_leaves = []
        outlier_species = set() # Store species IDs of outliers
        
        # Get the set of leaves that were considered outliers *at this node*
        node_outlier_leaves = pog_node.props.get('leaves_out', set())

        for leaf_name in member_names:
            
            leaf_node = tree[leaf_name]
            taxid = leaf_node.props.get('taxid')

            # Check if this leaf was an outlier *for this POG's node*
            if leaf_name not in node_outlier_leaves:
                retained_leaves.append(leaf_name)
                retained_species.add(taxid)
                
                # Annotate the leaf with its POG
                old_pog_annot = leaf_node.props.get('pOG', '')
                new_pog_annot = f"{old_pog_annot}@{pog_name}"
                leaf_node.add_prop('pOG', new_pog_annot)
            else:
                outlier_species.add(taxid)
            
        values = [
            pog_name,
            lca_pog,
            sci_name_lca,
            pog_node.name,
            retained_species,
            retained_leaves,
            '-',  # lca_dup
            list(outlier_species),
            pog_node.props.get('inparalogs_rate'),
            pog_node.props.get('so_score')
        ]
                  
        _add_entry_to_ogs_info(ogs_info, values)



# --- 3. Root OGs ---

def _process_root_ogs(tree: PhyloTree, tax_db: TaxonomyDB, ogs_info: Dict):
    """
    Handles two cases for the root:
    1. If the root's LCA is not LUCA, create OGs for its ancestors.
    2. If the root itself is marked as a monophyletic OG, add it.
    """
    
    lca_tree = tree.props.get('lca_node')
    
    # Case 1: Create ancestral OGs above the root's LCA
    # (e.g., if tree is 'Bacteria', create OGs for 'cellular organisms')
    if lca_tree not in ['1', '131567', 'r_root', 'Unk']:
        lineage = set( map(str, utils.get_lineage(tax_db, lca_tree) ))



        if '1' in lineage: lineage.remove('1') # Remove 'root'
        
        # We want ancestors *not* including the root's own LCA
        try:
            lineage.remove(lca_tree) 
        except ValueError:
            pass 
            
        all_leaves = tree.props.get('leaves_in', [])
        
        for taxid in lineage:
            og_name = f"{taxid}|{tree.name}*" # '*' denotes ancestral OG
            sci_name = utils.get_scientific_name(tax_db, taxid)
            values = [
                og_name, taxid, sci_name, tree.name,
                tree.props.get('sp_in', set()),
                all_leaves,
                '-', # lca_dup
                tree.props.get('sp_out', set()),
                tree.props.get('inparalogs_rate'),
                tree.props.get('so_score')
            ]
            _add_entry_to_ogs_info(ogs_info, values)

    # Case 2: The root itself is a monophyletic OG
    if tree.props.get('monophyletic_og') is True:
        sci_name_lca = utils.get_scientific_name(tax_db, lca_tree)
        og_name = f"{lca_tree}|{tree.name}"
        values = [
            og_name, lca_tree, sci_name_lca, tree.name,
            tree.props.get('sp_in', set()),
            tree.props.get('leaves_in', set()),
            '-', # lca_dup
            tree.props.get('sp_out', set()),
            tree.props.get('inparalogs_rate'),
            tree.props.get('so_score')
        ]
        _add_entry_to_ogs_info(ogs_info, values)



# --- 4. Hierarchy Calculation ---

def _calculate_hierarchy(ogs_info: Dict, tax_db: TaxonomyDB):
    """Calculates and adds 'OG_up' and 'OG_down' relationships."""
    
    # Cache lineages to avoid redundant DB calls
    lineage_cache = {}
    og_names = list(ogs_info.keys())
    
    for og_name in og_names:
        info = ogs_info[og_name]
        ogs_up = set()
        ogs_down = set()
        
        target_mems = set(info['Mems'])
        target_depth = _get_depth_from_og_name(og_name)
        target_taxid = _get_taxid_from_og_name(og_name)
        if target_taxid not in lineage_cache:
            lineage_cache[target_taxid] = set(utils.get_lineage(tax_db, target_taxid))
        
        for other_og in og_names:
            if og_name == other_og:
                continue
            
            other_info = ogs_info[other_og]
            
            # Check for shared members
            if not target_mems.isdisjoint(other_info['Mems']):
                other_depth = _get_depth_from_og_name(other_og)
                
                if target_depth < other_depth:
                    ogs_down.add(other_og)
                elif target_depth > other_depth:
                    ogs_up.add(other_og)
                else: # Depths are equal, must compare taxonomy
                    other_taxid = _get_taxid_from_og_name(other_og)
                    if other_taxid == target_taxid:
                        continue
                    
                    if other_taxid not in lineage_cache:
                        lineage_cache[other_taxid] = set(utils.get_lineage(tax_db, other_taxid))
                    
                    if target_taxid in lineage_cache[other_taxid]:
                        ogs_down.add(other_og) 
                    elif other_taxid in lineage_cache[target_taxid]:
                        ogs_up.add(other_og) 
        
        info['OG_down'] = list(ogs_down) if ogs_down else '-'
        info['OG_up'] = list(ogs_up) if ogs_up else '-'


def _get_depth_from_og_name(og_name: str) -> int:
    """Helper to safely parse depth from an OG name."""
    
    try:
        # Name is like 'LCA|NODE_NAME' where NODE_NAME is 'A-1' or 'A-1!'
        node_name = og_name.split('|', 1)[1]
        depth_str = node_name.split('-', 1)[1]    
        return int(re.sub(r'[\*!]', '', depth_str)) # Remove markers and convert
    
    except Exception:
        logging.warning(f"Could not parse depth from OG name: {og_name}")
        return -1

def _get_taxid_from_og_name(og_name: str) -> int:
    """Helper to safely parse taxid from an OG name."""
    try:
        return int(og_name.split('|', 1)[0])
    except Exception:
        return -1


# --- 5. Utility ---

def _add_entry_to_ogs_info(ogs_info: Dict, values: list):
    """Utility function to create a structured dictionary entry for an OG."""
    
    (og_name, lca, sci_name_lca, assoc_node_name, sp_in_og, 
     list_seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap) = values

    if not sp_outliers:
        num_sp_out = 0
        sp_outliers_list = ['-']
    else:
        num_sp_out = len(sp_outliers)
        sp_outliers_list = list(sp_outliers)
    
    
    

    ogs_info[og_name] = {
        'TaxoLevel': str(lca),
        'SciName_TaxoLevel': sci_name_lca.replace(' ', '_'),
        'AssocNode': assoc_node_name,
        'NumSP': len(sp_in_og),
        'NumMems': len(list_seqs_in),
        'Mems': list(list_seqs_in),
        'Lca_Dup': str(lca_dup),
        'Species_Outliers': sp_outliers_list,
        'Num_SP_Outliers': num_sp_out,
        'Inparalogs_Rate': float(inparalogs_rate),
        'SP_overlap_dup': species_overlap
    }

    