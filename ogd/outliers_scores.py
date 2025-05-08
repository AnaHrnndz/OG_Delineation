from ete4 import PhyloTree
import numpy as np
from collections import defaultdict, Counter, OrderedDict
import ogd.utils as utils



class OrderedCounter(Counter, OrderedDict):
     'Counter that remembers the order elements are first seen'
     def __repr__(self):
         return '%s(%r)' % (self.__class__.__name__,
                            OrderedDict(self))
     def __reduce__(self):
         return self.__class__, (OrderedDict(self),)



##  3. Outliers and Scores
def run_outliers_and_scores(t_nw, taxonomy_db, num_total_sp, level2sp_mem, args):

    mssg = f"""
    2. Detect outlierts and get scores:
        -Outliers thresholds:
            Best Taxa Threshold: {args.best_tax_thr}
            Lineage Threshold: {args.lineage_thr} 
        -Species losses percentage threshold: {args.sp_loss_perc}"""
    print(mssg)

    t = PhyloTree(t_nw, parser = 0)

    CONTENT = t.get_cached_content()
   
    t.add_prop('is_root', str('True'))

    sp_in_root = [l.props.get('taxid') for l in t]
    t.add_prop('sp_in', sp_in_root)
    

    # Detect long branches
    long_leaves = detect_long_branches(t)

    for n in t.traverse("preorder"):

        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            n.del_prop('named_lineage')
        
        n.del_prop('_speciesFunction')

        if n.is_leaf:
            n.add_prop('lca_node', n.props.get('taxid'))

        else:
            n.add_prop('old_lca_name',(n.props.get('sci_name')))

            #DETECT OUTLIERS
            sp_out = set()
            leaves_out = set()
            sp_in = list()
            leaves_in = set()
            
            #Add outliers from upper nodes
            if args.inherit_out == 'Yes':
                sp_out_inherit = add_upper_outliers(n, CONTENT)
                sp_out.update(sp_out_inherit)

            #   Detect outliers
            ch1 = n.children[0]
            ch2 = n.children[1]

            # Call outlierts detection and update the set with the species taxid that are outliers
            sp_out.update(outliers_detection(n, args.lineage_thr, args.best_tax_thr, CONTENT, level2sp_mem, sp_out, taxonomy_db))

            #  Save info for children_node_1 and children_node_2
            ch1_name = n.children[0].props.get('name')
            ch1_leaf_names = list(ch1.leaf_names())
            sp_ch1 = set()
            leaves_ch1 = set()
            
            ch2_name = n.children[1].props.get('name')
            ch2_leaf_names = list(ch2.leaf_names())
            sp_ch2 = set()
            leaves_ch2 = set()

            all_leafs = CONTENT[n]
            for l in all_leafs:
                if l.name in long_leaves:
                    leaves_out.add(l.props.get('name'))

                # Leaves that are outliers have the prop called "taxo_outlier"
                elif  str(l.props.get('taxid')) in sp_out:
                    leaves_out.add(l.props.get('name'))
                    l.add_prop('taxo_outlier', 'true')

                else:
                    sp_in.append(str(l.props.get('taxid')))
                    leaves_in.add(l.props.get('name'))
                    
                    if l.name in ch1_leaf_names:
                        sp_ch1.add(str(l.props.get('taxid')))
                        leaves_ch1.add(l.props.get('name'))
                    elif l.name in ch2_leaf_names:
                        sp_ch2.add(str(l.props.get('taxid')))
                        leaves_ch2.add(l.props.get('name'))

           
            #   Re-calculate species overlap after detect outliers
            overlaped_spces = set(sp_ch1 & sp_ch2)

            if len(overlaped_spces)>0:
                so_score = round(float(len(overlaped_spces) / len(set(sp_in))),4)
                n.add_prop('overlaped_species', list(overlaped_spces))
            else:
                so_score = 0.0


            #   Update last common ancestor, rank and lineage for the node after detect outliers
            update_taxonomical_props(n, set(sp_in), taxonomy_db)

            
            #   Save properties
            n.add_prop('sp_in', set(sp_in))
            n.add_prop('len_sp_in', len(set(sp_in)))
            n.add_prop('sp_in_ch1', sp_ch1)
            n.add_prop('sp_in_ch2', sp_ch2)
            n.add_prop('ch1_name', ch1_name)
            n.add_prop('ch2_name', ch2_name)
            n.add_prop('leaves_in', leaves_in)
            n.add_prop('total_leaves', len(n))
            n.add_prop('len_leaves_in', len(leaves_in))
            n.add_prop('len_leaves_out', len(leaves_out))
            n.add_prop('leaves_ch1',list(leaves_ch1))
            n.add_prop('leaves_ch2',list(leaves_ch2))
            n.add_prop('leaves_out', list(leaves_out))
            n.add_prop('so_score', so_score)
            n.add_prop('sp_out', list(sp_out))



            #   Load_node_scores add properties: score1, score2,  inpalalogs_rate, dup_score and species losses
            inparalogs_rate = get_inparalogs_rate(sp_in)
            n.add_prop('inparalogs_rate', str(inparalogs_rate))

            score1 = get_score1(len(set(sp_in)), num_total_sp)
            n.add_prop('score1', score1)

            score2 = get_score2(len(set(sp_in)), len(leaves_in))
            n.add_prop('score2', score2)

            dup_score = get_dup_score(n)
            n.add_prop('dup_score', dup_score)

            sp_lost(n, level2sp_mem)


            if n.up:
                sp_in_pnode = n.up.props.get('sp_in')
                lin_lost_from_pnode = '@'.join(map(str, best_lin_lost(expected_sp=sp_in_pnode, found_sp=set(sp_in), taxonomy_db=taxonomy_db)))
                n.add_prop('lost_from_uppernode', lin_lost_from_pnode)


            #   Based on species overlap threshold define duplication, false duplication and speciation nodes
            lin_lca = n.props.get('lineage')
            lca_node = n.props.get('lca_node')

            so_2_use = utils.get_so2use(taxonomy_db, lin_lca, lca_node,args)
        
            if so_score >= so_2_use:
                
                if float(n.props.get('species_losses_percentage')) > args.sp_loss_perc:
                    n.add_prop('evoltype_2', 'FD')
                else:
                    n.add_prop('evoltype_2', 'D')
            else:
                n.add_prop('evoltype_2', 'S')

    t, props = utils.run_clean_properties(t)

    return t, CONTENT



def detect_long_branches(t):

    """
        Detect leaves with long branches
        Long branches leaves are when a leave's branch is 50 times bigger than the mean of the branches for the whole tree
        Remove entire nodes do not perform well when basal branches are too large (ex Phou fam)
    """

    length_leaves = list()
    total_dist = list()
    for n in t.traverse():
        # Skip root dist (Root dist is None)
        if n.dist != None:
            if n.is_leaf:
                length_leaves.append(n.dist)
                total_dist.append(n.dist)
            else:
                total_dist.append(n.dist)

    
    mean_length = np.mean(total_dist)
    long_leaves = set()
    # Only check dist for leaves, not for internal nodes
    for l in t:
        if l.dist > mean_length*50:
            long_leaves.add(l.name)
            l.add_prop('long_branch_outlier', 'True')

    return long_leaves



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

    lineages = OrderedCounter()
    nspecies = 0
    for sp in sp_list:
        nspecies += 1

        lineages.update(utils.get_lineage(taxonomy_db, sp))

    lca = [l for l, count in lineages.items() if count == nspecies]
    
    if not lca:
        lca = ['Unk']

    
    return lca[-1]


def get_inparalogs_rate(sp_in):

    """
        Calculate the median number of inparalogs
        inparalogs = seqs that belong to the same species
    """

    dups_per_sp = Counter(sp_in)
    
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
       Get Duplication Score (Dup score)
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

    