import ogd.utils as utils
from ete4 import PhyloTree # t is assumed to be a PhyloTree
from typing import Any, Tuple, Set, List, Dict

## 4. Detect Duplications and Core-OGs & PGs ##

def run_get_main_dups(t, taxonomy_db, total_mems_in_tree, args):
    """
    Selects high-quality duplication nodes (HQ-Dups) to be marked as Monophyletic OGs.

    HQ-Dups criteria:
        1. Must be a duplication node (evoltype_2='D').
        2. Must have >1 leaf and >1 species (non-trivial group).
        3. Logic: If child nodes do NOT contain another HQ-Dup at the same LCA level, 
           the children nodes are marked as Monophyletic OGs.
    """

    taxid_dups_og = set()
    
    
    so_euk = args.so_euk if args.so_euk is not None else 'None'
    so_bact = args.so_bact if args.so_bact is not None else 'None'
    so_arq = args.so_arq if args.so_arq is not None else 'None'

    mssg = f"""
       -Species overlap threshold:
            General: {args.so_all}
            Euk: {so_euk}
            Bact: {so_bact}
            Arq: {so_arq}"""
    print(mssg)
    
    # Optimized traversal to find HQ duplications
    for node in t.traverse("preorder"):

        # The root is handled at the end; leaf nodes are skipped
        if node.is_leaf:
            continue
        
        # Criterion 1 & 2: Must be a Duplication and non-trivial
        is_hq_candidate = (
            node.props.get('evoltype_2') == 'D' and 
            len(node.props.get('leaves_in', [])) > 1 and 
            len(node.props.get('sp_in', [])) > 1
        )

        if is_hq_candidate:
            
            lca_target = node.props.get('lca_node')
            
            # --- Criterion 3: Check for HQ duplications below the current node ---
            
            # 3.1. Search for D duplications under the current node
            ch1, ch2 = node.children[0], node.children[1]

            # Search for D-nodes under CH1 and CH2 with the same LCA_TARGET
            dups_under_ch1 = list(ch1.search_nodes(evoltype_2='D', lca_node=lca_target))
            dups_under_ch2 = list(ch2.search_nodes(evoltype_2='D', lca_node=lca_target))
            
            # 3.2. Validate if the duplications found under the children meet the HQ requirement
            
            # Determine if an HQ-Dup exists under CH1 with the same LCA
            save_dups_ch1 = _count_hq_dups_under_child(dups_under_ch1)
            
            # Determine if an HQ-Dup exists under CH2 with the same LCA
            save_dups_ch2 = _count_hq_dups_under_child(dups_under_ch2)

            
            # 3.3. Annotate Monophyletic OGs
            
            # If there are no HQ-Dups under CH1, CH1 is a monophyletic OG
            if save_dups_ch1 == 0:
                _annotate_og_child(taxid_dups_og, node, ch1, taxonomy_db)

            # If there are no HQ-Dups under CH2, CH2 is a monophyletic OG
            if save_dups_ch2 == 0:
                _annotate_og_child(taxid_dups_og, node, ch2, taxonomy_db)


            # Annotate that this parent node created OGs
            if save_dups_ch1 == 0 or save_dups_ch2 == 0:
                node.add_prop('node_create_og', 'True')

    # 3. Decide if the root is a monophyletic OG
    _annotate_root_og(t)

    # 4. Final cleanup
    t, props = utils.run_clean_properties(t)

    return t, taxid_dups_og




def _count_hq_dups_under_child(dups_under_child):
    """
    Counts how many duplication nodes in the list meet the HQ criteria
    (>1 leaf and >1 species).
    """
    count = 0
    for n_ in dups_under_child:
        if len(n_.props.get('leaves_in', [])) > 1 and len(n_.props.get('sp_in', [])) > 1:
            count += 1
    return count


def _annotate_og_child(taxid_dups_og, parent_node, target_node, taxonomy_db):
    """
    Adds properties to a child node that has been identified as a Monophyletic OG.
    """
    
    # The child must also be non-trivial
    sp_ch = target_node.props.get('sp_in', [])
    og_ch_mems = target_node.props.get('leaves_in', [])
    
    if len(sp_ch) > 1 and len(og_ch_mems) > 1:
        target_node.add_prop('monophyletic_og', 'True')
        target_node.add_prop('lca_dup', parent_node.props.get('lca_node'))
        target_node.add_prop('so_score_dup', parent_node.props.get('so_score'))

        #  Get and save the lineage only once
        dup_lca = parent_node.props.get('lca_node')
        dup_lin = utils.get_lineage(taxonomy_db, dup_lca)
        target_node.add_prop('dup_lineage', dup_lin)
        
        target_node.add_prop('dup_node_name', parent_node.props.get('name'))

        taxid_dups_og.add(dup_lca)
    

def _annotate_root_og(t):
    """
    Decides if the root node should be marked as Monophyletic OG.
    """
    lca_root = t.props.get('lca_node')
    
    # Use the 'lca_node' property directly in the search.
    # If there are no monophyletic OGs sharing the root's LCA, the root is an OG.
    if not list(t.search_nodes(monophyletic_og='True', lca_node=lca_root)):
        t.add_prop('monophyletic_og', 'True')