

from collections import Counter, defaultdict
from typing import Set, Dict, Any, Tuple, List
import logging
import numpy as np
from ete4 import PhyloTree

import ogd.utils as utils

# --- Main  Function ---

def run_outliers_and_scores(
    t: PhyloTree, 
    taxonomy_db: Any, 
    num_total_sp: int, 
    level2sp_mem: Dict, 
    args: Any
) -> Tuple[PhyloTree, Dict, Set]:
    """
    Main function for Step 3: Detects outliers, calculates scores, and annotates the tree.
    """
    logging.info(f"Parameters: BestTaxaThr={args.best_tax_thr}, LineageThr={args.lineage_thr}, SpLossPerc={args.sp_loss_perc}.")
    logging.info(f"Species Overlap Thresholds: General: {args.so_all}, Eukaryotes: {args.so_euk or 'N/A'}, Bacteria: {args.so_bact or 'N/A'}, Archaea: {args.so_arq or 'N/A'}  ")
    
    content = t.get_cached_content()
    _initialize_root_node(t, content)
    
    long_leaves = detect_long_leaves_branches(t)
    total_leaf_outliers = set(long_leaves)

    for node in t.traverse("preorder"):
        
        # Clean up properties to ensure a clean run
        if not utils.is_gtdb(taxonomy_db):
            node.del_prop('named_lineage')
        node.del_prop('_speciesFunction')

        if node.is_leaf:
            node.add_prop('lca_node', node.props.get('taxid'))
        else:
            parent_outliers = set()
            if node.up:
                parent_outliers = set( node.up.props.get('sp_out') )
               
            # The processing of internal nodes 
            leaves_out_node = _process_internal_node(
                node, parent_outliers, taxonomy_db, num_total_sp, level2sp_mem, content, long_leaves, args
            )
            total_leaf_outliers.update(leaves_out_node)
    
    logging.info(f"Identified {len(long_leaves)} leaves as long-branch outliers.")
    logging.info(f"Identified {len(total_leaf_outliers)} leaves as taxonomical outliers.")
    
    t, props = utils.sanitize_tree_properties(t)
    
    return t, content, total_leaf_outliers




def _initialize_root_node(t: PhyloTree, content: Dict):
    """Adds the initial properties to the tree's root node."""
    t.add_prop('is_root', 'True')
    sp_in_root = [str(leaf.props.get('taxid')) for leaf in content[t]]
    t.add_prop('sp_in', sp_in_root)
    

def detect_long_leaves_branches(t: PhyloTree) -> Set[str]:
    """Detects leaves with long branches (50x the tree's average)."""
    all_dists = [n.dist for n in t.traverse() if n.dist is not None]
    if not all_dists:
        return set()
    
    mean_length = np.mean(all_dists)
    threshold = mean_length * 50
    
    long_leaves = set()
    for leaf in t:
        if leaf.dist > threshold:
            long_leaves.add(leaf.name)
            leaf.add_prop('long_branch_outlier', 'True')
    
    return long_leaves


def _process_internal_node(
    n: PhyloTree, parent_outliers: Set[int], taxonomy_db: Any, num_total_sp: int, level2sp_mem: Dict, 
    content: Dict, long_leaves: Set[str], args: Any
) -> Set[str]:
    """
    Processes a single internal node: detects outliers, classifies leaves, 
    calculates scores, and determines the evolutionary event type.
    """

    n.add_prop('old_lca_name', n.props.get('sci_name'))

    # 1. Outlier detection
    sp_out = _detect_node_outliers(n, content, parent_outliers, level2sp_mem, taxonomy_db, args)
    
    # 2. Leaf classification 
    node_classification = _classify_leaves(n, content, sp_out, long_leaves)

  
    # 3. Overlap calculation and taxonomy update
    so_score = _calculate_species_overlap(n, node_classification['sp_in_ch1'], node_classification['sp_in_ch2'], node_classification['sp_in'])
    
    update_taxonomical_props(n, node_classification['sp_in'], taxonomy_db)
    
    # 4. Score calculation and species loss
    _calculate_and_set_scores(n, node_classification, sp_out, num_total_sp, level2sp_mem, taxonomy_db)
    
    # 5. Determination of the evolutionary event type (D, FD, S)
    _determine_evolutionary_event(n, taxonomy_db, so_score, args)

    # 6. Save remaining properties
    _set_final_node_properties(n, node_classification, sp_out, so_score)

    return node_classification['leaves_out']




# --- Sub-functions of _process_internal_node ---
# 1. Outlier detection
def _detect_node_outliers(n: PhyloTree, content: Dict, parent_outliers: Set[int], level2sp_mem: Dict, taxonomy_db: Any, args: Any) -> Set[str]:
    """Detects and returns a set of outlier species for the node."""
    
    sp_out = set()
    
    #1.1 Inherit taxonomic outliers from parent node
    if not args.no_inherit_outliers:
        parent_out_in_node = add_upper_outliers(n, content, parent_outliers)
        sp_out.update(parent_out_in_node)
    
    #1.2 Detect node-specific taxonomic outliers  
    sp_out.update(outliers_detection(n, args.lineage_thr, args.best_tax_thr, content, level2sp_mem, sp_out, taxonomy_db))
    
    return sp_out

def add_upper_outliers(n: PhyloTree, content: Dict, parent_outliers: Set[int]) -> Set[int]:
    """Inherits outliers from the parent if its species are in the current node ."""
    
    node_species = set(map(str, [leaf.props.get('taxid') for leaf in content[n]]))
    
    return parent_outliers.intersection(node_species)

def outliers_detection(n, lineage_thr, best_tax_thr, content, level2sp_mem, sp_out_up, taxonomy_db):

    """
        Detect outliers for each taxid from lineages
        1. Create taxonomy counter for the node -> sp_per_level_in_node
        2. Select the best taxonomic level for the node. Best level its the most recent level that represent >90% sp in the node
        3. For the sp that do not belong to the best tax level (candidates2remove), check if are outliers or not: 
            for each sp in candidates2remove, get the lineage, for each taxlev in the lineage calculate:
                sp_in_tree
                lin_in_tree
    """

    sp_in_node = set()

    """
    Create a taxonomy counter for the node. 
    For each lineage that is present in the node, count how many species there are in the node.
    """
    sp_per_level_in_node = defaultdict(set)
    
    for l in content[n]:
        if str(l.props.get('taxid')) not in sp_out_up:
            sp_in_node.add(str(l.props.get('taxid')))
            utils.update_sp_per_level_in_node(sp_per_level_in_node, taxonomy_db, l)

    
    """
    Select best taxonomic level for the node. 
    Best level its the most recent level that grouped >90% sp in the node
    To know that is the most recel, we need to know the "depth" of the taxa, 
    Ex. cell org depth's == 1, Euk, bact and arq == 2, Metazoa == 4
    """
    best_tax = str()
    ptax = defaultdict(dict)
    
    for tax, sp_list in sp_per_level_in_node.items():
        
        n_sp_list = len(sp_list)
        depth_tax = len(utils.get_lineage(taxonomy_db, tax))

        perc_tax_in_node = n_sp_list/len(sp_in_node)
        ptax[depth_tax][tax] = perc_tax_in_node
        

    """
    Once I have the depth of all the lineages, sort the list in reverse to traverse 
    from the most recent to the oldest (the oldest is always cell_org). 
    The first level with more than 90% of species will be set as the 'best taxon'.
    """
    
    depths_list = (list(ptax.keys()))
    depths_list.sort(reverse=True)

    best_tax = str()
    best_depth = int()
    for depth in depths_list:
        for tax, perc in ptax[depth].items():
            if perc >=best_tax_thr:
                best_tax = (tax)
                best_depth = depth
                break
        if best_tax != '':
            break

    
    if best_depth in ptax.keys():
        n.add_prop('best_tax', str(best_tax)+'_'+str(ptax[best_depth][best_tax]))
    else:
        n.add_prop('best_tax', 'NO LEAVES')


    best_tax_lineage =  utils.get_lineage(taxonomy_db, best_tax)
    
    # For the sp that do not belong to the best_tax level, check if are outliers or not
    sp_keep = sp_per_level_in_node[best_tax]
    candidates2remove = sp_in_node.difference(sp_keep)
    sp2remove = set()

    for sp in candidates2remove:
        
        total_lin = utils.get_lineage(taxonomy_db, sp)
        
        lin2use = set(total_lin).difference(set(best_tax_lineage))
        
        for l in lin2use:

            """
            For a given lineage, we evaluate its representation in a node compared to the entire tree.
            For example, if the tree contains 100 bacterial species, we check how many of these bacteria are present in the node. 
            If a lineage in a node includes less than 5% of its species from the entire tree, 
            all species from that lineage are considered outliers.
            """
            
            lin_rareness  = len(sp_per_level_in_node[str(l)]) / len(level2sp_mem[(str(l))])

            if lin_rareness < lineage_thr:
                sp2remove.add(sp)
                break
    
    return sp2remove


# 2. Leaf classification 
def _classify_leaves(n: PhyloTree, content: Dict, sp_out: Set[str], long_leaves: Set[str]) -> Dict:
    """Classifies the node's leaves into 'in' and 'out' and groups them by child."""
    
    classification = {
        'leaves_in': set(), 'leaves_in_names': set(), 'leaves_out': set(), 'sp_in':set(),
        'sp_in_ch1': set(), 'leaves_in_ch1': set(),
        'sp_in_ch2': set(), 'leaves_in_ch2': set()
    }
    
    ch1_leaf_names = set(n.children[0].leaf_names())
    
    n.add_prop('sp_in_ch1',set())
    n.add_prop('sp_in_ch2',set())
    n.add_prop('sp_in',set())
    
    for leaf in content[n]:
        
        taxid_str = str(leaf.props.get('taxid'))
        is_outlier = (leaf.name in long_leaves) or (taxid_str in sp_out)
        
        if is_outlier:
            classification['leaves_out'].add(leaf.name)
            
            if taxid_str in sp_out:
                leaf.add_prop('taxo_outlier', 'true')
        else:
            classification['leaves_in'].add(leaf)
            classification['leaves_in_names'].add(leaf.name)
            classification['sp_in'].add(taxid_str)
            n.props.get('sp_in').add(taxid_str)    
            if leaf.name in ch1_leaf_names:
                classification['sp_in_ch1'].add(taxid_str)
                classification['leaves_in_ch1'].add(leaf.name)
                n.props.get('sp_in_ch1').add(taxid_str)
            else:
                classification['sp_in_ch2'].add(taxid_str)
                classification['leaves_in_ch2'].add(leaf.name)
                n.props.get('sp_in_ch2').add(taxid_str)
    
    return classification


# 3. Overlap calculation and taxonomy update
def _calculate_species_overlap(n: PhyloTree, sp_in_ch1: Set[str], sp_in_ch2: Set[str], sp_in_set: Set[str]) -> float:
    """Calculates and returns the species overlap score."""
    
    overlapped_species = sp_in_ch1.intersection(sp_in_ch2)
    n.add_prop('overlaped_species', list(overlapped_species))
    
    return round(len(overlapped_species) / len(sp_in_set), 4) if sp_in_set else 0.0

def update_taxonomical_props(n, sp_in, taxonomy_db):

    """
        Update taxonomical information after detect outliers
    """

    if len(sp_in) > 0:
        
        lca_node = utils.get_lca_node(sp_in, taxonomy_db)
        
        if lca_node == 'r_root':
            lin_lca = ['r_root']
            rank = ['r_root']

        else:

            rank = utils.get_rank(taxonomy_db, lca_node)
            lin_lca = utils.get_lineage(taxonomy_db, lca_node)
            #lca_node_name = utils.get_lca_node_name(taxonomy_db, lca_node)
            lca_node_name = utils.get_scientific_name(taxonomy_db, lca_node)

        n.add_prop('lineage', lin_lca)
        n.add_prop('taxid', lca_node)
        n.add_prop('lca_node', lca_node)
        n.add_prop('lca_node_name', lca_node_name.replace(":", ""))
        n.add_prop('rank', rank)
        n.props['sci_name'] = n.props.get('lca_node_name')
    
    elif len(sp_in) == 0:
    
        n.add_prop('lineage', ['Empty'])
        n.add_prop('taxid', 'Empty')
        n.add_prop('lca_node', 'Empty')
        n.add_prop('lca_node_name', 'Empty')
        n.add_prop('rank', 'Empty')
        n.props['sci_name'] = n.props.get('Empty')


# 4. Scores calculation and species loss
def _calculate_and_set_scores(n: PhyloTree, node_classification: Dict, sp_out: Set[str], num_total_sp: int, level2sp_mem: Dict, taxonomy_db: Any):
    """Calculates and adds all scores (inparalogs, dup_score, s1, s2, etc.) to the node."""
    
    n.add_prop('inparalogs_rate', get_inparalogs_rate(node_classification['leaves_in']))
    
    n.add_prop('score1', get_score1(len(node_classification['sp_in']), num_total_sp))
    
    n.add_prop('score2', get_score2(len(node_classification['sp_in']), len(node_classification['leaves_in']))) 
    
    n.add_prop('dup_score', get_dup_score(n))
    
    sp_lost(n, sp_out, level2sp_mem)
   
    if n.up and n.up.props.get('sp_in'):
        
        lost_lineage = best_lin_lost(
            expected_sp=set(n.up.props['sp_in']),
            found_sp=node_classification['sp_in'],
            taxonomy_db=taxonomy_db
        )
       
        n.add_prop('lost_from_uppernode', '@'.join(map(str, lost_lineage)))


def get_inparalogs_rate(leaves_in):
    """Calculate the median number of inparalogs inparalogs = seqs that belong to the same species"""

    list_sp_in = [l.props.get('taxid') for l in leaves_in]
    
    dups_per_sp = Counter(list_sp_in)
    
    inparalogs_rate = np.median(list(dups_per_sp.values()))
    
    return inparalogs_rate

def get_score1(len_sp_in, num_total_sp):
   
    score1 = float(len_sp_in/num_total_sp)
    return score1

def get_score2(len_sp_in, nseqs):

    score2 = float(len_sp_in/nseqs) if nseqs else 0.0
    return score2

def get_dup_score(n):

    """
       Duplication Score: number of shared species in both children nodes divided by the smallest number of species
       Similiar to species overlap, but species overlap uses the sum of all the species found in children nodes
       Duplication score use the number of species of the children node with less species.
    """

    sp1 = n.props.get('sp_in_ch1')
    sp2 = n.props.get('sp_in_ch2')

    
    if len(sp1) == 0 and len(sp2) == 0:
        return 0

    a = np.array([len(sp1), len(sp2)])
    minval = int(np.min(a[np.nonzero(a)]))

    dup_score = float(len(sp1 & sp2) / minval)
    
    return dup_score

def sp_lost(n, sp_out, level2sp_mem):
    """
    Calculate the number and percentage of species that have been lost in the node for the lca.
    If the lca of the node is bacteria and there are 10 species of bacteria, but in the tree there are in total 20 species of bacteria
    10 species will have been lost, corresponding to 50%.
    """
    
    lca_node = str(n.props.get('lca_node'))
    sp_in = set(n.props.get('sp_in'))
    
    if lca_node == 'Empty':
        diff_sp =  sp_out
        perc_diff_sp = 1.0

    else:
        
        expected_sp = level2sp_mem[lca_node]

        if len(sp_in) == 0 or len(expected_sp) == 0:
            diff_sp = expected_sp.difference(sp_in)
            perc_diff_sp = 1.0
        else:
            diff_sp = expected_sp.difference(sp_in)
            perc_diff_sp = len(diff_sp) / len(expected_sp)

    n.add_prop('species_losses', len(diff_sp))
    n.add_prop('species_losses_percentage', perc_diff_sp)


def best_lin_lost(expected_sp: Set[str], found_sp: Set[str], taxonomy_db: Any) -> list:
    """Identifies the most recent and significant lineage loss event."""
    
    lost_species = expected_sp.difference(found_sp)
    if not lost_species:
        return []

    sp_per_taxon = defaultdict(set)
    lineage_cache = {}

    for sp in expected_sp:
        lineage = utils.get_lineage(taxonomy_db, sp)
        lineage_cache[sp] = lineage
        for taxon in lineage:
            sp_per_taxon[taxon].add(sp)

    taxa_with_losses = {taxon for sp in lost_species for taxon in lineage_cache.get(sp, [])}
    
    loss_candidates = []
    for taxon in taxa_with_losses:
        expected_in_taxon = sp_per_taxon[taxon]
        lost_in_taxon = lost_species.intersection(expected_in_taxon)
        percentage_lost = len(lost_in_taxon) / len(expected_in_taxon)

        if percentage_lost >= 0.80:
            # Get depth from any species in the taxon group to avoid extra DB calls
            depth = len(lineage_cache[next(iter(expected_in_taxon))])
            if depth > 2:
                loss_candidates.append((depth, taxon, percentage_lost))

    if not loss_candidates:
        return []

    loss_candidates.sort() # Sorts by depth (the first element of the tuple)
    best_loss = loss_candidates[0]
    return [best_loss[1], best_loss[2]] # Return [taxon, percentage]


# 5. Determination of the evolutionary event type
def _determine_evolutionary_event(n: PhyloTree, taxonomy_db: Any, so_score: float, args: Any):
    """Determines and adds the 'evoltype_2' property (D, FD, S) to the node."""
    so_threshold = utils.determine_so_threshold(taxonomy_db, n.props.get('lineage', []), args)
    
    if so_score >= so_threshold:
        
        is_false_dup = float(n.props.get('species_losses_percentage', 0)) > args.sp_loss_perc
        n.add_prop('evoltype_2', 'FD' if is_false_dup else 'D')
        
    else:
        n.add_prop('evoltype_2', 'S')


# 6. Save remaining properties
def _set_final_node_properties(n: PhyloTree, classification: Dict, sp_out: Set[str], so_score: float):
    """Adds all final calculated properties to the node."""
    n.add_prop('sp_in', list(set(classification['sp_in'])))
    n.add_prop('len_sp_in', len(set(classification['sp_in'])))
    n.add_prop('sp_in_ch1', list(classification['sp_in_ch1']))
    n.add_prop('sp_in_ch2', list(classification['sp_in_ch2']))
    n.add_prop('ch1_name', n.children[0].props.get('name'))
    n.add_prop('ch2_name', n.children[1].props.get('name'))
    n.add_prop('leaves_in', list(classification['leaves_in_names']))
    n.add_prop('total_leaves', len(n))
    n.add_prop('len_leaves_in', len(classification['leaves_in']))
    n.add_prop('len_leaves_out', len(classification['leaves_out']))
    n.add_prop('leaves_ch1', list(classification['leaves_in_ch1']))
    n.add_prop('leaves_ch2', list(classification['leaves_in_ch2']))
    n.add_prop('leaves_out', list(classification['leaves_out']))
    n.add_prop('so_score', so_score)
    n.add_prop('sp_out', list(sp_out))

