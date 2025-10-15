# -*- coding: utf-8 -*-
"""
outliers_scores.py - Step 3 of the pipeline.
Detects outliers, calculates duplication and evolutionary scores,
and annotates the tree with these properties.
"""

import sys
from collections import defaultdict, Counter
from typing import Any, Tuple, Set, List, Dict

import numpy as np
from ete4 import PhyloTree
import ogd.new_utils as utils
from ogd.new_utils import BaseTaxonomyHandler, run_clean_properties

# --- Constants ---
LONG_BRANCH_MULTIPLIER = 50.0

# --- Main Function ---

def run_outliers_and_scores(
    t: PhyloTree, 
    tax_handler: BaseTaxonomyHandler, 
    num_total_sp: int, 
    level2sp_mem: Dict[str, Set[str]], 
    args: Any
) -> Tuple[PhyloTree, Dict, Set[str]]:
    """
    Main function for detecting outliers and calculating scores.

    Args:
        t (PhyloTree): The processed PhyloTree object from the setup step.
        tax_handler (BaseTaxonomyHandler): The initialized taxonomy handler.
        num_total_sp (int): Total number of unique species in the tree.
        level2sp_mem (Dict[str, Set[str]]): Maps taxonomic levels to species sets.
        args (Any): The command-line arguments object.

    Returns:
        A tuple containing:
        - The fully annotated PhyloTree object.
        - The cached content of the tree.
        - A set of all outlier leaf names.
    """
    print(f"""        - Outlier thresholds:
            - Best Taxa: {args.best_tax_thr}
            - Lineage: {args.lineage_thr}
        - Species loss threshold: {args.sp_loss_perc}""")

    CONTENT = t.get_cached_content()
    t.add_prop('is_root', True) # Use boolean, not string
    
    # Initialize root properties
    t.add_prop('sp_in', {l.props.get('taxid') for l in t})

    # 1. Detect long branches, which are considered outliers.
    total_outliers = detect_long_branches(t)

    # 2. Traverse the tree to process each internal node.
    for n in t.traverse("preorder"):
        
        _clean_and_initialize_node(n, tax_handler)

        if n.is_leaf:
            continue

        # --- OUTLIER DETECTION AND LEAF ASSIGNMENT ---
        leaves_out_set, leaves_in_set, sp_in_taxids_list = _process_node_outliers(
            n, CONTENT, total_outliers, tax_handler, args, level2sp_mem
        )
        total_outliers.update(leaves_out_set)
        print(len(total_outliers))

        # --- RE-CALCULATE TAXONOMY AND SPECIES OVERLAP (SO) ---
        sp_in_set = set(sp_in_taxids_list)
        update_taxonomical_props(n, sp_in_set, tax_handler)

        ch1, ch2 = n.children
        
        # Determine species/leaves in each child among the IN-GROUP members
        leaves_ch1_in = leaves_in_set.intersection(CONTENT[ch1])
        leaves_ch2_in = leaves_in_set.intersection(CONTENT[ch2])
        
        sp_ch1_in = {leaf.props.get('taxid') for leaf in leaves_ch1_in}
        sp_ch2_in = {leaf.props.get('taxid') for leaf in leaves_ch2_in}

        overlaped_species = sp_ch1_in & sp_ch2_in
        so_score = len(overlaped_species) / len(sp_in_set) if sp_in_set else 0.0

        # --- SAVE PRIMARY PROPERTIES ---
        _save_node_properties(n, locals())

        # --- SCORE CALCULATION AND EVOLTYPE ASSIGNMENT ---
        _calculate_scores_and_evoltype(n, sp_in_taxids_list, num_total_sp, level2sp_mem, tax_handler, args)

    # 3. Final cleaning of property strings.
    t, _ = run_clean_properties(t)

    return t, CONTENT, total_outliers

# --- Helper and Sub-functions ---

def detect_long_branches(t: PhyloTree) -> Set[str]:
    """
    Detects leaves with branches longer than a dynamic threshold.
    A branch is "long" if it's > N times the mean branch length of the tree.
    """
    all_distances = [n.dist for n in t.traverse() if n.dist is not None and n.dist > 0]
    if not all_distances:
        return set()

    mean_length = np.mean(all_distances)
    threshold = mean_length * LONG_BRANCH_MULTIPLIER
    
    long_leaves = set()
    for leaf in t:
        if leaf.dist is not None and leaf.dist > threshold:
            long_leaves.add(leaf.name)
            leaf.add_prop('long_branch_outlier', True)
            
    return long_leaves

def _clean_and_initialize_node(n: PhyloTree, tax_handler: BaseTaxonomyHandler) -> None:
    """Cleans unneeded properties and initializes default properties for a node."""
    # Use handler type to decide which properties to clean
    if tax_handler.__class__.__name__ == 'NCBITaxonomyHandler':
        n.del_prop('named_lineage')
    n.del_prop('_speciesFunction')

    if n.is_leaf:
        n.add_prop('lca_node', n.props.get('taxid'))
    else:
        n.add_prop('old_lca_name', n.props.get('sci_name'))
        n.add_prop('leaves_out', [])
        n.add_prop('sp_out', [])


def _process_node_outliers(
    n: PhyloTree, CONTENT: Dict, long_leaves: Set[str], tax_handler: BaseTaxonomyHandler, 
    args: Any, level2sp_mem: Dict[str, Set[str]]
) -> Tuple[Set[str], Set[str], List[str]]:
    """Detects and separates leaves into outliers and in-group members for a single node."""
    
    sp_out_current = set()
   
    
    if args.inherit_out: # Check boolean directly
        sp_out_current.update(add_upper_outliers(n, CONTENT))
    

    newly_detected_outliers = outliers_detection(
        n, args.lineage_thr, args.best_tax_thr, CONTENT, 
        level2sp_mem, sp_out_current, tax_handler
    )
    sp_out_current.update(newly_detected_outliers)

    leaves_out = set()
    leaves_in = set()
    sp_in_taxids = []

    for leaf in CONTENT[n]:
        is_long_branch = leaf.name in long_leaves
        is_tax_outlier = leaf.props.get('taxid') in sp_out_current
        
        if is_long_branch or is_tax_outlier:
            leaves_out.add(leaf.name)
            if is_tax_outlier:
                leaf.add_prop('taxo_outlier', True)
        else:
            leaves_in.add(leaf.name)
            sp_in_taxids.append(leaf.props.get('taxid'))
            
    n.add_prop('leaves_out', list(leaves_out))
    n.add_prop('sp_out', list(sp_out_current))
    
    return leaves_out, leaves_in, sp_in_taxids



def add_upper_outliers(n: PhyloTree, CONTENT: Dict) -> Set[str]:
    """Inherits species outliers from the parent node if they are present in the current node."""
    if not n.up or 'sp_out' not in n.up.props:
        return set()
    
    sp_out_up = set(n.up.props['sp_out'])
    node_species = {leaf.props.get('taxid') for leaf in CONTENT[n]}
    
    return sp_out_up.intersection(node_species)

def outliers_detection(
    n: PhyloTree, 
    lineage_thr: float, 
    best_tax_thr: float, 
    CONTENT: Dict, 
    level2sp_mem: Dict, 
    sp_out_up: Set[str], 
    tax_handler: BaseTaxonomyHandler
) -> Set[str]:
    """
    Detects taxonomic outliers for a given node.

    This process involves three main steps:
    1. Find the "best" taxonomic level that represents the majority of species in the node.
    2. Identify species in the node that do not belong to this best-fit lineage (candidates for removal).
    3. For each candidate, check if its own lineage is poorly represented in the node
       compared to its representation in the entire tree. If so, it's an outlier.

    Returns:
        A set of species taxids identified as outliers for this node.
    """
    
    # 1. Create a taxonomy counter for the non-outlier species currently in the node.
    sp_per_level_in_node = defaultdict(set)
    sp_in_node = set()
    for leaf in CONTENT[n]:
        taxid = leaf.props.get('taxid')
        if taxid not in sp_out_up:
            sp_in_node.add(taxid)
            # Use the handler to update species counts for each taxonomic level in the leaf's lineage
            tax_handler.update_species_per_level(leaf, sp_per_level_in_node)
            
    if not sp_in_node:
        return set()

    # 2. Determine the best taxonomic level for the node.
    # The "best" level is the most specific (deepest) one that contains at least `best_tax_thr` % of the species.
    ptax = defaultdict(dict)
    for tax, species_set in sp_per_level_in_node.items():
        # The depth is simply the length of the lineage.
        depth = len(tax_handler.get_lineage(tax))
        perc_in_node = len(species_set) / len(sp_in_node)
        ptax[depth][tax] = perc_in_node

    best_tax = ''
    # Iterate from the deepest (most specific) level upwards.
    for depth in sorted(ptax.keys(), reverse=True):
        
        for tax, perc in ptax[depth].items():
            if perc >= best_tax_thr:
                best_tax = str(tax)
                break
        if best_tax:
            break
            
    if not best_tax:
        # If no level meets the threshold, no outliers can be determined this way.
        return set()
    
   
    n.add_prop('best_tax', best_tax)

    # 3. Identify outliers among the candidate species.
    # Candidates are species in the node that are NOT in the best-fit lineage.
    species_to_keep = sp_per_level_in_node.get(str(best_tax), set())
    
    candidates_to_remove = sp_in_node - species_to_keep
    species_to_remove = set()
    best_tax_lineage = set(tax_handler.get_lineage(best_tax))
    
    for sp in candidates_to_remove:
        # For each candidate, find the parts of its lineage that are more specific than the node's best_tax.
        total_lineage = set(tax_handler.get_lineage(sp))
        lineage_to_check = total_lineage - best_tax_lineage
        
        for level in lineage_to_check:
            # How many species of this 'level' are in the node vs. the whole tree?
            species_of_level_in_node = sp_per_level_in_node.get(str(level), set())
            species_of_level_in_tree = level2sp_mem.get(str(level), set())
            
            if not species_of_level_in_tree:
                continue

            # Calculate the rareness of this lineage within the current node.
            lineage_rareness = len(species_of_level_in_node) / len(species_of_level_in_tree)

            # If the lineage is very rare in this node, the species is an outlier.
            if lineage_rareness < lineage_thr:
                species_to_remove.add(sp)
                break # Move to the next candidate species
    
    return species_to_remove



def _save_node_properties(node: PhyloTree, local_vars: Dict) -> None:
    """Helper to save all calculated properties to a node."""
    # This pattern makes the main loop cleaner by centralizing property assignment
    props_to_add = {
        'so_score': round(local_vars.get('so_score', 0.0), 4),
        'overlaped_species': list(local_vars.get('overlaped_species')),
        'sp_in': local_vars.get('sp_in_set'),
        'len_sp_in': len(local_vars.get('sp_in_set')),
        'leaves_in': local_vars.get('leaves_in_set'),
        'len_leaves_in': len(local_vars.get('leaves_in_set')),
        'len_leaves_out': len(local_vars.get('leaves_out_set')),
        'total_leaves': len(node),
        'ch1_name': local_vars.get('ch1').name,
        'ch2_name': local_vars.get('ch2').name,
        'leaves_ch1': {leaf.name for leaf in local_vars.get('leaves_ch1_in')},
        'leaves_ch2': {leaf.name for leaf in local_vars.get('leaves_ch2_in')},
        'sp_in_ch1': local_vars.get('sp_ch1_in'),
        'sp_in_ch2': local_vars.get('sp_ch2_in')
    }
    for key, value in props_to_add.items():
        node.add_prop(key, value)








def get_lca_from_species(sp_list: Set[str], tax_handler: BaseTaxonomyHandler) -> str:
    """Finds the Last Common Ancestor (LCA) for a given set of species."""
    if not sp_list:
        return 'Unk'

    lineages = Counter()
    for sp in sp_list:
        lineages.update(tax_handler.get_lineage(sp))
    
    # The LCA is the last taxon present in all lineages
    lca_candidates = [taxon for taxon, count in lineages.items() if count == len(sp_list)]
    return lca_candidates[-1] if lca_candidates else 'Unk'


def update_taxonomical_props(n: PhyloTree, sp_in_set: Set[str], tax_handler: BaseTaxonomyHandler) -> None:
    """Updates taxonomic properties (LCA, rank, lineage) of a node based on its in-group species."""
    if not sp_in_set:
        lca_node, lca_name, rank, lineage = 'Empty', 'Empty', 'Empty', ['Empty']
    else:
        lca_node = get_lca_from_species(sp_in_set, tax_handler)
        if lca_node == 'r_root' or lca_node == 'Unk':
            lca_name, rank, lineage = lca_node, lca_node, [lca_node]
        else:
            lca_name = tax_handler.get_lca_node_name(lca_node)
            rank = tax_handler.get_rank(lca_node)
            lineage = tax_handler.get_lineage(lca_node)

    n.add_prop('lineage', lineage)
    n.add_prop('taxid', lca_node)
    n.add_prop('lca_node', lca_node)
    n.add_prop('lca_node_name', lca_name.replace(":", ""))
    n.add_prop('rank', rank)
    n.props['sci_name'] = n.props.get('lca_node_name')



def _calculate_scores_and_evoltype(
    n: PhyloTree, sp_in_taxids_list: List[str], num_total_sp: int, 
    level2sp_mem: Dict, tax_handler: BaseTaxonomyHandler, args: Any
) -> None:
    """Calculates all scores and determines the evolutionary type for a node."""
    
    sp_in_set = set(sp_in_taxids_list)
    len_sp_in = len(sp_in_set)
    len_leaves_in = n.props.get('len_leaves_in', 0)

    # Basic scores
    score1 = len_sp_in / num_total_sp if num_total_sp > 0 else 0.0
    score2 = len_sp_in / len_leaves_in if len_leaves_in > 0 else 0.0
    
    # Inparalogs
    dups_per_sp = Counter(sp_in_taxids_list)
    inparalogs_rate = np.median(list(dups_per_sp.values())) if dups_per_sp else 0.0

    # Duplication Score
    sp1 = n.props.get('sp_in_ch1', set())
    sp2 = n.props.get('sp_in_ch2', set())
    shared_len = len(sp1 & sp2)
    min_len = min(len(s) for s in (sp1, sp2) if s) if (sp1 or sp2) else 0
    dup_score = shared_len / min_len if min_len > 0 else 0.0

    # Species Losses
    lca_node = str(n.props.get('lca_node'))
    expected_sp = level2sp_mem.get(lca_node, set())
    lost_sp = expected_sp - sp_in_set
    perc_lost = len(lost_sp) / len(expected_sp) if expected_sp else 1.0

    # Save scores
    n.add_prop('score1', score1)
    n.add_prop('score2', score2)
    n.add_prop('inparalogs_rate', inparalogs_rate)
    n.add_prop('dup_score', dup_score)
    n.add_prop('species_losses', len(lost_sp))
    n.add_prop('species_losses_percentage', perc_lost)

    # Determine evoltype
    lin_lca = n.props.get('lineage')
    so_threshold_to_use = utils.determine_so_threshold(lin_lca, args, tax_handler)
    
    if n.props.get('so_score', 0.0) >= so_threshold_to_use:
        n.add_prop('evoltype_2', 'FD' if perc_lost > args.sp_loss_perc else 'D')
    else:
        n.add_prop('evoltype_2', 'S')

    # Note: best_lin_lost function logic would also be refactored to use the tax_handler
    # but is omitted here for brevity as it's a large, complex function.
        
