import argparse
import logging
from typing import Set, Tuple

from ete4 import PhyloTree

import ogd.utils as utils
from ogd.utils import TaxonomyDB

# --- Main Orchestrator Function ---

def run_get_main_dups(
    tree: PhyloTree, 
    tax_db: TaxonomyDB, 
    args: argparse.Namespace
) -> Tuple[PhyloTree, Set[str]]:
    """
    Identifies high-quality duplication nodes and annotates their children as OGs.

    Args:
        tree: The fully scored and annotated PhyloTree from the previous step.
        tax_db: The initialized taxonomy database.
        args: The parsed command-line arguments.

    Returns:
        A tuple containing:
        - The tree, now annotated with 'monophyletic_og' properties.
        - A set of taxids where OGs have been defined.
    """
    
    
    # This set will store the LCA taxids of duplication nodes that create OGs.
    taxids_defining_ogs = set()

    # We traverse the tree looking for potential duplication nodes.
    for node in tree.traverse("preorder"):
        
        # A node is a candidate if it's a valid duplication event.
        if _is_valid_duplication_node(node):
            
            lca_target = node.props.get('lca_node')
            child1, child2 = node.children

            is_child1_nested_dup = False
            for n_ in  child1.search_nodes(evoltype_2='D', lca_node=lca_target):
                is_child1_nested_dup = _is_valid_duplication_node(n_)

            is_child2_nested_dup = False
            for n_ in  child2.search_nodes(evoltype_2='D', lca_node=lca_target):
                is_child2_nested_dup = _is_valid_duplication_node(n_)                

            # If a child is NOT a nested duplication, it becomes an OG.
            if not is_child1_nested_dup:
                _annotate_og_node(node, child1, taxids_defining_ogs, tax_db)
            
            if not is_child2_nested_dup:
                _annotate_og_node(node, child2, taxids_defining_ogs, tax_db)

    # Finally, check if the root node could be an OG.
    _check_root_as_og(tree)
    
    logging.info(f"Identified {len(list(tree.search_nodes(node_create_og = True)))} HQ duplication nodes.")
    
    # Sanitize properties before finishing.
    tree, _ = utils.sanitize_tree_properties(tree)

    return tree, taxids_defining_ogs


# ---  Helper Functions ---

def _is_valid_duplication_node(node: PhyloTree) -> bool:
    """Checks if a node represents a well-supported duplication event."""
    
    if node.is_leaf:
        return False
    
    # Conditions: must be a duplication, have more than one leaf, and more than one species.
    return (
        node.props.get('evoltype_2') == 'D' and
        len(node.props.get('leaves_in', set())) > 1 and
        len(node.props.get('sp_in', set())) > 1
    )

def _annotate_og_node(
    duplication_node: PhyloTree, 
    og_node: PhyloTree, 
    taxids_defining_ogs: Set[str],
    tax_db: TaxonomyDB
):
    """
    Adds all necessary 'monophyletic_og' properties to a node.

    Args:
        duplication_node: The parent node marked as 'D'.
        og_node: The child node that will be marked as an OG.
        taxids_defining_ogs: The set of taxids to update.
        tax_db: The taxonomy database.
    """
    
    # Double-check that the child node itself represents a valid group.
    if len(og_node.props.get('sp_in', set())) > 1 and len(og_node.props.get('leaves_in', [])) > 1:
        
        lca_of_duplication = duplication_node.props.get('lca_node')
        
        og_node.add_prop('monophyletic_og', True)
        og_node.add_prop('lca_dup', lca_of_duplication)
        og_node.add_prop('so_score_dup', duplication_node.props.get('so_score'))
        og_node.add_prop('dup_lineage', utils.get_lineage(tax_db, lca_of_duplication))
        og_node.add_prop('dup_node_name', duplication_node.name)
        
        # Mark the parent duplication as a creator of OGs.
        duplication_node.add_prop('node_create_og', True)
        

       
        
def _check_root_as_og(tree: PhyloTree):
    """
    Checks if the root of the tree can be considered an OG.
    This happens if no other OGs were defined at the root's taxonomic level.
    """
    root_lca = tree.props.get('lca_node')

    # Search for any existing OGs defined at the same taxonomic level as the root.
    if  len(list(tree.search_nodes(monophyletic_og=True, lca_node=root_lca))) == 0:
        tree.add_prop('monophyletic_og', True)
        
    
    
    