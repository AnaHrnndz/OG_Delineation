from ete4 import PhyloTree
import numpy as np
from collections import defaultdict, Counter, OrderedDict
from typing import Any, Tuple, Set, List, Dict
import ogd.utils as utils
import sys




## --- Main  Function ---

def run_outliers_and_scores(t_nw: str, taxonomy_db: Any, num_total_sp: int, level2sp_mem: Dict[str, Set[str]], args: Any) -> Tuple[PhyloTree, Dict, Set[str]]:
    """
    Main function for detecting outliers, calculating scores, and assigning properties.
    """

    print(f"""
        - Outliers thresholds:
            Best Taxa Threshold: {args.best_tax_thr}
            Lineage Threshold: {args.lineage_thr} 
        - Species losses percentage threshold: {args.sp_loss_perc}
    """)

    t = PhyloTree(t_nw, parser = 0)
    CONTENT = t.get_cached_content()
   
    t.add_prop('is_root', str('True'))
    
    # Initialize root properties
    t.add_prop('sp_in', [l.props.get('taxid') for l in t])

    # 1. Detect long branches
    total_outliers = set()
    long_leaves = detect_long_branches(t)
    total_outliers.update(long_leaves)
    

    # 2. Traverse and process nodes
    for n in t.traverse("preorder"):
        
        _clean_and_initialize_node(n, taxonomy_db)

        if n.is_leaf:
            # total_outliers.update([n.name] if n.props.get('long_branch_outlier') == 'True' else set())
            continue

        # --- OUTLIERS DETECTION AND LEAF ASSIGNMENT ---
        sp_out_set, leaves_out_set, sp_in_taxids_list, leaves_in_set = _process_node_outliers(
            n, CONTENT, long_leaves, taxonomy_db, args, level2sp_mem
        )
        
        total_outliers.update(leaves_out_set) # Update master list of outliers

        # --- RE-CALCULATE TAXONOMY AND SPECIES OVERLAP (SO) ---
        sp_in_set = set(sp_in_taxids_list)
        update_taxonomical_props(n, sp_in_set, taxonomy_db)

        ch1, ch2 = n.children[0], n.children[1]
        
        # Determine species/leaves in each child among the IN-GROUP members
        leaves_ch1_in = leaves_in_set.intersection(set(ch1.leaf_names()))
        leaves_ch2_in = leaves_in_set.intersection(set(ch2.leaf_names()))
        
        sp_ch1_in = set([n[l].props.get('taxid') for l in leaves_ch1_in])
        sp_ch2_in = set([n[l].props.get('taxid') for l in leaves_ch2_in])

        # Species Overlap (SO) Score Calculation
        overlaped_spces = sp_ch1_in & sp_ch2_in
        if sp_in_set:
            so_score = round(float(len(overlaped_spces) / len(sp_in_set)), 4)
        else:
            so_score = 0.0

        # --- SAVE PRIMARY PROPERTIES ---
        n.add_prop('so_score', so_score)
        n.add_prop('overlaped_species', list(overlaped_spces))
        n.add_prop('sp_in', sp_in_set)
        n.add_prop('len_sp_in', len(sp_in_set))
        n.add_prop('leaves_in', leaves_in_set)
        n.add_prop('len_leaves_in', len(leaves_in_set))
        n.add_prop('len_leaves_out', len(leaves_out_set))
        n.add_prop('total_leaves', len(n))
        
        n.add_prop('ch1_name', ch1.props.get('name'))
        n.add_prop('ch2_name', ch2.props.get('name'))
        n.add_prop('leaves_ch1', list(leaves_ch1_in))
        n.add_prop('leaves_ch2', list(leaves_ch2_in))
        n.add_prop('sp_in_ch1', sp_ch1_in)
        n.add_prop('sp_in_ch2', sp_ch2_in)


        # --- SCORE CALCULATION AND EVOLTYPE ASSIGNMENT ---
        _calculate_scores(n, sp_in_taxids_list, num_total_sp, level2sp_mem)
        _determine_and_assign_evoltype(n, taxonomy_db, args)

    # 3. Final cleaning
    t, props = utils.run_clean_properties(t)

    return t, CONTENT, total_outliers






def detect_long_branches(t):

    """
        Detect leaves with long branches
        Long branches leaves are when a leave's branch is 50 times bigger than the mean of the branches for the whole tree
    """
    # 1. Collect all distances 
    total_dist = [n.dist for n in t.traverse() if n.dist is not None]
    
    if not total_dist:
        return set()

    # 2. Calculate mean length and threshold
    mean_length = np.mean(total_dist)
    threshold = mean_length * 50
    long_leaves = set()

    # 3. Identify and tag long branches
    for l in t.leaves():
        if l.dist is not None and l.dist > threshold:
            long_leaves.add(l.name)
            l.add_prop('long_branch_outlier', 'True')

    return long_leaves



def _clean_and_initialize_node(n, taxonomy_db):
    """Clean unneeded properties and initialize default properties."""
    
    # Clean properties added by ETE4 during initial annotation/setup
    if "ncbi_taxonomy" in str(taxonomy_db):
        n.del_prop('named_lineage')
    n.del_prop('_speciesFunction')

    # Initialize properties
    if n.is_leaf:
        n.add_prop('lca_node', n.props.get('taxid'))
    else:
        # Save original LCA name before outlier-based re-calculation
        n.add_prop('old_lca_name', n.props.get('sci_name'))
        n.add_prop('leaves_out', list())
        n.add_prop('sp_out', list())



def _process_node_outliers(n, CONTENT, long_leaves, taxonomy_db, args, level2sp_mem):
    """Detects and separates leaves into outliers and in-group members."""
    
    sp_out_current = set()
    leaves_out = set()
    leaves_in = set()
    sp_in_taxids = [] 

    # 1. Inherit species outliers from upper nodes
    if args.inherit_out == 'Yes':
        sp_out_current.update(add_upper_outliers(n, CONTENT))

    # 2. Detect taxonomic outliers
    newly_detected_outliers = outliers_detection(n, args.lineage_thr, args.best_tax_thr, CONTENT, level2sp_mem, sp_out_current, taxonomy_db)
    sp_out_current.update(newly_detected_outliers)

    # 3. Classify all leaves based on long branches or taxo outliers
    all_leaves = CONTENT[n]
    for leaf in all_leaves:
        leaf_name = leaf.props.get('name')
        taxid = str(leaf.props.get('taxid'))

        is_long_branch = leaf_name in long_leaves
        is_tax_outlier = taxid in sp_out_current

        if is_long_branch or is_tax_outlier:
            leaves_out.add(leaf_name)
            if is_tax_outlier:
                leaf.add_prop('taxo_outlier', 'true')
        else:
            sp_in_taxids.append(taxid)
            leaves_in.add(leaf_name)
    
    # 4. Save final outlier set to the node
    n.add_prop('leaves_out', list(leaves_out))
    n.add_prop('sp_out', list(sp_out_current))
    
    return sp_out_current, leaves_out, sp_in_taxids, leaves_in

def add_upper_outliers(n, CONTENT):

    """
        Get outliers species in upper nodes
    """

    sp_out_inherit = set()
    if n.up and n.up.props.get('sp_out'):
        sp_out_up = n.up.props.get('sp_out')
        for sp in sp_out_up:
            if (sp) in [leaf.props.get('taxid') for leaf in CONTENT[n]]:
                sp_out_inherit.add(sp)
    
    return sp_out_inherit

def outliers_detection(n, lineage_thr, best_tax_thr, CONTENT, level2sp_mem, sp_out_up, taxonomy_db):

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
    
    for l in CONTENT[n]:
        if l.props.get('taxid') not in sp_out_up:
            sp_in_node.add(l.props.get('taxid'))
            utils.update_sp_per_level_in_node(sp_per_level_in_node, taxonomy_db, l)
            
    """
    Select best taxonomic level for the node. 
    Best level its the most recent level that grouped >90% sp in the node
    To know that is the most recel, we need to know the "depth" of the taxa, 
    Ex. cell org depth's == 1, Euk, bact and arq == 2, Metazoa == 4
    """
    best_tax = str()
    ptax = defaultdict(dict)
    
    for tax, num in sp_per_level_in_node.items():
        
        depth_tax = len(utils.get_lineage(taxonomy_db, tax))

        perc_tax_in_node = len(num)/len(sp_in_node)
        ptax[depth_tax][tax] = perc_tax_in_node


    """
    Once that I have the depth of all the lineages, sort the list in reverse 
    to traverse from more recent to the oldest (oldest always will be cell_org)
    First level with more that 90% of species will be the best_tax
    """
    depths_list = (list(ptax.keys()))
    depths_list.sort(reverse=True)

    best_tax = str()
    best_depth = int()
    for depth in depths_list:
        for tax, perc in ptax[depth].items():
            if perc >=best_tax_thr:
                best_tax = str(tax)
                best_depth = depth
                break
        if best_tax != '':
            break

    if best_depth in ptax.keys():
        n.add_prop('best_tax', best_tax+'_'+str(ptax[best_depth][best_tax]))
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

def update_taxonomical_props(n, sp_in, taxonomy_db):

    """
        Update taxonomical information after detect outliers
    """

    if len(sp_in) > 0:
        lca_node = get_lca_node(sp_in, taxonomy_db)
        
        if lca_node == 'r_root':
            lin_lca = ['r_root']
            rank = ['r_root']

        else:

            rank = utils.get_rank(taxonomy_db, lca_node)
            lin_lca = utils.get_lineage(taxonomy_db, lca_node)
            lca_node_name = utils.get_lca_node_name(taxonomy_db, lca_node)

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

def get_lca_node(sp_list, taxonomy_db):

    """
        Get Last Common Ancestor (lca) from a set of species
    """

    lineages = utils.OrderedCounter()
    nspecies = 0
    for sp in sp_list:
        nspecies += 1

        lineages.update(utils.get_lineage(taxonomy_db, sp))

    lca = [l for l, count in lineages.items() if count == nspecies]
    
    if not lca:
        lca = ['Unk']

    
    return lca[-1]



def _calculate_scores(n, sp_in_taxids, num_total_sp, level2sp_mem):
    """Calculates Score1, Score2, DupScore, Inparalogs, and Species Losses."""
    
    sp_in_set = set(sp_in_taxids)
    len_sp_in = len(sp_in_set)
    len_leaves_in = len(n.props.get('leaves_in'))

    # Score 1 and 2
    n.add_prop('score1', get_score1(len_sp_in, num_total_sp))
    n.add_prop('score2', get_score2(len_sp_in, len_leaves_in))

    # Inparalogs Rate
    n.add_prop('inparalogs_rate', str(get_inparalogs_rate(sp_in_taxids)))

    # Dup Score
    n.add_prop('dup_score', get_dup_score(n))

    # Species Losses (uses node's new LCA)
    sp_lost(n, level2sp_mem)

def get_score1(len_sp_in, num_total_sp):
   
    score1 = float(len_sp_in/num_total_sp)
    return score1

def get_score2(len_sp_in, nseqs):

    score2 = float(len_sp_in/nseqs) if nseqs else 0.0
    return score2

def get_inparalogs_rate(sp_in):

    """
        Calculate the median number of inparalogs
        inparalogs = seqs that belong to the same species
    """

    dups_per_sp = Counter(sp_in)
    
    inparalogs_rate = np.median(list(dups_per_sp.values()))
    
    return inparalogs_rate

def get_dup_score(n):

    """
       Get Duplication Score (Dup score)
       Duplication Score: number of shared species in both children nodes divided by the smallest number of species
       Similiar to species overlap, but species overlap uses the sum of all the species found in children nodes
       Duplication score use the number of species of the children node with less species.
    """

    sp1 = n.props.get('sp_in_ch1')
    sp2 = n.props.get('sp_in_ch2')
    
    len1 = len(sp1)
    len2 = len(sp2)

    # Intersection length
    shared_len = len(sp1 & sp2)
    
    if len1 == 0 and len2 == 0:
        return 0.0
    
    # Find the minimum non-zero length
    min_len = min(l for l in [len1, len2] if l > 0)
    
    # Calculate score
    dup_score = shared_len / min_len
    
    return float(dup_score)

def sp_lost(n, level2sp_mem):

    """
    Calculate the number and percentage of species that have been lost in the node for the lca.
    If the lca of the node is bacteria and there are 10 species of bacteria, but in the tree there are in total 20 species of bacteria
    10 species will have been lost, corresponding to 50%.
    """

    lca_node = str(n.props.get('lca_node'))
    sp_in = set(n.props.get('sp_in'))

    expected_sp = set(level2sp_mem[lca_node])
    
    if len(sp_in) == 0 or len(expected_sp) == 0:
        diff_sp = expected_sp.difference(sp_in)
        perc_diff_sp = 1.0
    else:
        diff_sp = expected_sp.difference(sp_in)
        perc_diff_sp = len(diff_sp) / len(expected_sp)

    n.add_prop('species_losses', len(diff_sp))
    n.add_prop('species_losses_percentage', perc_diff_sp)



def _determine_and_assign_evoltype(n, taxonomy_db, args):
    """Determines the evolutionary type (D, FD, or S) based on SO score."""
    
    so_score = n.props.get('so_score')
    
    lin_lca = n.props.get('lineage')
    so_2_use = utils.determine_so_threshold(taxonomy_db, lin_lca, args)
    
    if so_score >= so_2_use:
        sp_loss_perc = n.props.get('species_losses_percentage')
        if float(sp_loss_perc) > args.sp_loss_perc:
            n.add_prop('evoltype_2', 'FD') # False Duplication
        else:
            n.add_prop('evoltype_2', 'D') # Duplication
    else:
        n.add_prop('evoltype_2', 'S') # Speciation

    # Calculate Lineage Lost from upper node
    if n.up and n.up.props.get('sp_in'):
        sp_in_pnode = n.up.props.get('sp_in')
        lin_lost_from_pnode = '@'.join(map(str, best_lin_lost(expected_sp=sp_in_pnode, found_sp=n.props.get('sp_in'), taxonomy_db=taxonomy_db)))
        n.add_prop('lost_from_uppernode', lin_lost_from_pnode)

def best_lin_lost(expected_sp, found_sp, taxonomy_db):

    diff_sp = expected_sp.difference(found_sp)

    sp_lost_per_level = defaultdict(set)
    for sp in diff_sp:

        lin2use = utils.get_lineage(taxonomy_db, sp)
        
        for tax in lin2use:
            sp_lost_per_level[tax].add(sp)

    sp_exp_per_level = defaultdict(set)
    for sp in expected_sp:

        lin2use = utils.get_lineage(taxonomy_db, sp)
        
        for tax in lin2use:
            sp_exp_per_level[tax].add(sp)

    ptax = defaultdict(dict)
    for tax, num in sp_lost_per_level.items():
        depth_tax = len(utils.get_lineage(taxonomy_db, sp))
        
        perc_tax_in_node = len(num)/len(sp_exp_per_level[tax]) 
        ptax[depth_tax][tax] = perc_tax_in_node

    depths_list = (list(ptax.keys()))
    depths_list.sort()

    best_loss = []

    for depth in depths_list:
        if depth not in [1, 2]:
            for tax, perc in ptax[depth].items():

                if  perc >=0.80:
                    best_loss = [tax, perc]
                    break
    
    return best_loss

