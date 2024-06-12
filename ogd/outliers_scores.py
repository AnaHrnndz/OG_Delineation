from ete4 import PhyloTree
import numpy as np
from collections import defaultdict, Counter, OrderedDict
import utils



class OrderedCounter(Counter, OrderedDict):
     'Counter that remembers the order elements are first seen'
     def __repr__(self):
         return '%s(%r)' % (self.__class__.__name__,
                            OrderedDict(self))
     def __reduce__(self):
         return self.__class__, (OrderedDict(self),)



##  3. Outliers and Scores
def run_outliers_and_scores(t_nw, taxonomy_db, NUM_TOTAL_SP, level2sp_mem, args):

    print('\n'+'2. Detect outlierts and get scores')
    print(' -Outliers thresholds:')
    print('\t'+'Node: '+ str(args.outliers_node))
    print('\t'+'Reftree: '+ str(args.outliers_reftree))
    print(' -Species losses percentage threshold: ' + str(args.sp_loss_perc))

    t = PhyloTree(t_nw)

    CONTENT = t.get_cached_content()
   
    t.add_prop('is_root', str('True'))

    sp_in_root = [l.props.get('taxid') for l in t]
    t.add_prop('_sp_in', sp_in_root)

    # Detect long branches
    long_leaves = detect_long_branches(t)
    
    # if args.user_taxonomy_counter:
        # root_outliers = root_taxono_outliers(t, NUM_TOTAL_SP, level2sp_mem)
    # else:
    #    root_outliers = set()    

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
            sp_in = set()
            leaves_in = set()
            leaves_in_nodes = set()

            #Add outliers from upper nodes
            if args.inherit_out == 'Yes':
                #sp_out.update(root_outliers)
                sp_out_inherit = add_upper_outliers(n, CONTENT)
                sp_out.update(sp_out_inherit)

            #   Detect outliers
            ch1 = n.children[0]
            ch2 = n.children[1]

            sp_out.update(outliers_detection(n, args.outliers_node, args.outliers_reftree, CONTENT, level2sp_mem, sp_out, taxonomy_db))

            #   Save info for children_node_1
            ch1_name = n.children[0].props.get('name')
            ch1_leaf_names = list(ch1.leaf_names())
            sp_ch1 = set()
            leaves_ch1 = set()

            #   Save info for children_node_2
            ch2_name = n.children[1].props.get('name')
            ch2_leaf_names = list(ch2.leaf_names())
            sp_ch2 = set()
            leaves_ch2 = set()

            all_leafs = CONTENT[n]
            for l in all_leafs:
                if l.name in long_leaves:
                    leaves_out.add(l.props.get('name'))

                elif  str(l.props.get('taxid')) in sp_out:
                    leaves_out.add(l.props.get('name'))
                    l.add_prop('taxo_outlier', 'true')

                else:
                    sp_in.add(str(l.props.get('taxid')))
                    leaves_in.add(l.props.get('name'))
                    leaves_in_nodes.add(l)

                    if l.name in ch1_leaf_names:
                        sp_ch1.add(str(l.props.get('taxid')))
                        leaves_ch1.add(l.props.get('name'))
                    elif l.name in ch2_leaf_names:
                        sp_ch2.add(str(l.props.get('taxid')))
                        leaves_ch2.add(l.props.get('name'))


            #   Re-calculate species overlap after detect outliers
            overlaped_spces = set(sp_ch1 & sp_ch2)

            if len(overlaped_spces)>0:
                so_score = float(len(overlaped_spces) / len(sp_in))
                n.add_prop('_overlap', list(overlaped_spces))
            else:
                so_score = 0.0

            #   Update last common ancestor, rank and lineage for the node after detect outliers
            update_taxonomical_props(n, sp_in, taxonomy_db)

            #   Save properties
            n.add_prop('_sp_in', sp_in)
            n.add_prop('len_sp_in', len(sp_in))
            n.add_prop('_sp_in_ch1', sp_ch1)
            n.add_prop('_sp_in_ch2', sp_ch2)
            n.add_prop('_ch1_name', ch1_name)
            n.add_prop('_ch2_name', ch2_name)
            n.add_prop('_leaves_in', leaves_in)
            n.add_prop('_leaves_in_nodes', leaves_in_nodes)
            n.add_prop('total_leaves', len(n))
            n.add_prop('len_leaves_in', len(leaves_in))
            n.add_prop('len_leaves_out', len(leaves_out))
            n.add_prop('_leaves_ch1',list(leaves_ch1))
            n.add_prop('_leaves_ch2',list(leaves_ch2))
            n.add_prop('_leaves_out', list(leaves_out))
            n.add_prop('so_score', so_score)
            if len(sp_out) == 0:
                n.add_prop('sp_out', ['None'])
            else:
                n.add_prop('sp_out', list(sp_out))


            #   Load_node_scores add properties: score1, score2 and inpalalogs_rate

            inparalogs_rate = get_inparalogs_rate(list(sp_in))
            n.add_prop('inparalogs_rate', str(inparalogs_rate))

            score1 = get_score1(len(sp_in), NUM_TOTAL_SP)
            n.add_prop('score1', score1)

            score2 = get_score2(len(sp_in), len(leaves_in))
            n.add_prop('score2', score2)

            dup_score = get_dup_score(n)
            n.add_prop('dup_score', dup_score)

            #call_species_lost(n, reftree)

            sp_lost_v2(n, level2sp_mem)

            #call_lineage_lost(n, reftree, taxonomy_db)

            #lin_lost_from_root = best_lin_lost(expected_sp=sp_in_root, found_sp=sp_in, taxonomy_db=taxonomy_db)
            #n.add_prop('lost_from_root', lin_lost_from_root)
            if n.up:
                sp_in_pnode = n.up.props.get('_sp_in')
                lin_lost_from_pnode = best_lin_lost(expected_sp=sp_in_pnode, found_sp=sp_in, taxonomy_db=taxonomy_db)
                n.add_prop('lost_from_uppernode', lin_lost_from_pnode)


            #   Based on species overlap threshold define duplication, false duplication and speciation nodes
            lin_lca = n.props.get('lineage')
            lca_node = n.props.get('lca_node')

            so_2_use = utils.get_so2use(taxonomy_db, lin_lca, args)
        
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
        if n.is_leaf:
            length_leaves.append(n.dist)
        else:
            total_dist.append(n.dist)

    mean_length = np.mean(total_dist)

    long_leaves = set()
    for l in t:
        if l.dist > mean_length*50:
            long_leaves.add(l.name)
            l.add_prop('long_branch_outlier', 'True')


    return long_leaves

# def root_taxono_outliers(t, NUM_TOTAL_SP, level2sp_mem):
    # root_outliers = set()
    # for taxid in ['2', '2759', '2157']:

        # if taxid in level2sp_mem.keys():
            # if (len(level2sp_mem[taxid])/ NUM_TOTAL_SP) < 0.1:
                # root_outliers.update(level2sp_mem[taxid])

    # return root_outliers

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

def outliers_detection(n, outliers_node, outliers_reftree, CONTENT, level2sp_mem, sp_out_up, taxonomy_db):

    """
        Detect outliers for each taxid from lineages
        1. Create taxonomy counter for the node -> sp_per_level
        2. Select the best taxonomic level for the node. Best level its the most recent level that represent >90% sp in the node
        3. For the sp that do not belong to the best tax level (candidates2remove), check if are outliers or not: 
            for each sp in candidates2remove, get the lineage, for each taxlev in the lineage calculate:
                sp_in_tree
                lin_in_tree
    """


    sp_in_node = set()

    # Create taxonomy counter for the node. For each lineage that its present in the node, count how many species there are in the node
    sp_per_level = defaultdict(set)

    for l in CONTENT[n]:
        if l.props.get('taxid') not in sp_out_up:
            sp_in_node.add(l.props.get('taxid'))

            if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
                for tax in l.props.get('lineage').split('|'):
                    sp_per_level[tax].add(str(l.props.get('taxid')))
            elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy': 
                for tax in l.props.get('named_lineage').split('|'):
                    sp_per_level[tax].add(str(l.props.get('taxid')))


    # Select best taxonomic level for the node. 
    # Best level its the most recent level that represent >90% sp in the node
    best_tax = str()
    ptax = defaultdict(dict)
    
    for tax, num in sp_per_level.items():
        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            depth_tax = len(taxonomy_db.get_lineage(tax))
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy': 
            depth_tax = len(taxonomy_db.get_name_lineage([tax])[0][tax])

        perc_tax_in_node = len(num)/len(sp_in_node)
        ptax[depth_tax][tax] = perc_tax_in_node

    depths_list = (list(ptax.keys()))
    depths_list.sort(reverse=True)

    best_tax = str()
    best_depth = int()
    for depth in depths_list:
        for tax, perc in ptax[depth].items():
            if perc >=0.90:
                best_tax = str(tax)
                best_depth = depth
                break
        if best_tax != '':
            break

    if best_depth in ptax.keys():
        n.add_prop('best_tax', best_tax+'_'+str(ptax[best_depth][best_tax]))
    else:
        n.add_prop('best_tax', 'NO LEAVES')


    # For the sp that do not belong to the best tax level, check if are outliers or not
    sp_keep = sp_per_level[best_tax]
    candidates2remove = sp_in_node.difference(sp_keep)
    sp2remove = set()

    for sp in candidates2remove:
        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            lin2use =taxonomy_db.get_lineage(sp)
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
            lin2use = taxonomy_db.get_name_lineage([sp])[0][sp]
        
        for l in lin2use:

            # this rule prevent for remove strange lineages with few representatives in the tree
            # for each taxlev of the species, check how specific is that taxlev in the tree,
            # if taxlev has few representatives in the tree (close to 1), then sp_in_tree ~1
            sp_in_tree = 1/len(level2sp_mem[(str(l))])

            # For the lineage, how well represented is in the node in comparison with the whole tree
            # if lineage has almost all species lin_in_tree~1, if not lin_in_tree~0
            lin_in_tree = len(sp_per_level[str(l)]) / len(level2sp_mem[(str(l))])

            if sp_in_tree < 0.05 and lin_in_tree <0.05:
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
            if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
                rank = utils.clean_string(taxonomy_db.get_rank([lca_node])[lca_node])
                lin_lca = taxonomy_db.get_lineage(lca_node)
                lca_node_name = taxonomy_db.get_taxid_translator([lca_node])[lca_node]
            elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
                rank = utils.get_gtdb_rank(lca_node)
                lin_lca = taxonomy_db.get_name_lineage([lca_node])[0][lca_node]
                lca_node_name = lca_node

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
        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            lineages.update(taxonomy_db.get_lineage(sp))
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
            lineages.update(taxonomy_db.get_name_lineage([sp])[0][sp])

    lca = [l for l, count in lineages.items() if count == nspecies]
    
    if not lca:
        lca = ['Unk']

    # if lca[-1] == 'root':
        # lca[-1] = 'r_root'

    
    return lca[-1]

def get_inparalogs_rate(leaf_targets):

    """
        Calculate the median number of inparalogs
        inparalogs = seqs that belong to the same species
    """

    dups_per_sp = Counter(leaf_targets)
    inparalogs_rate = np.median(list(dups_per_sp.values()))

    return inparalogs_rate

def get_score1(sp_in, NUM_TOTAL_SP):

    score1 = float(sp_in/NUM_TOTAL_SP)
    return score1

def get_score2(sp_in, nseqs):

    score2 = float(sp_in/nseqs) if nseqs else 0.0
    return score2

def get_dup_score(n):

    """
       Get Duplication Score (Dup score)
       Duplication Score: number of shared species in both children nodes divided by the smallest number of species
       Similiar to species overlap, but species overlap uses the sum of all the species found in children nodes
       Suplication score use the number of species of the children node with less species.
    """

    sp1 = n.props.get('_sp_in_ch1')
    sp2 = n.props.get('_sp_in_ch2')

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
        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            lin2use =taxonomy_db.get_lineage(sp)
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
            lin2use = taxonomy_db.get_name_lineage([sp])[0][sp]
        
        for tax in lin2use:
            sp_lost_per_level[tax].add(sp)

    sp_exp_per_level = defaultdict(set)
    for sp in expected_sp:
        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            lin2use =taxonomy_db.get_lineage(sp)
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
            lin2use = taxonomy_db.get_name_lineage([sp])[0][sp]
        
        for tax in lin2use:
            
            sp_exp_per_level[tax].add(sp)

    ptax = defaultdict(dict)
    for tax, num in sp_lost_per_level.items():
        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            depth_tax = len(taxonomy_db.get_lineage(sp))
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
            depth_tax = len(taxonomy_db.get_name_lineage([sp])[0][sp])
        
        perc_tax_in_node = len(num)/len(sp_exp_per_level[tax]) 
        ptax[depth_tax][tax] = perc_tax_in_node

    depths_list = (list(ptax.keys()))
    depths_list.sort()

    best_loss = defaultdict(dict)

    for depth in depths_list:
        if depth not in [1, 2]:
            for tax, perc in ptax[depth].items():

                if  perc >=0.80:
                    best_loss[depth][tax] = perc

    return best_loss

def sp_lost_v2(n, level2sp_mem):

    lca_node = str(n.props.get('lca_node'))
    sp_in = set(n.props.get('_sp_in'))

    expected_sp = set(level2sp_mem[lca_node])
    

    if len(sp_in) == 0 or len(expected_sp) == 0:
        diff_sp = expected_sp.difference(sp_in)
        perc_diff_sp = 1.0
    else:
        diff_sp = expected_sp.difference(sp_in)
        perc_diff_sp = len(diff_sp) / len(expected_sp)

   

    n.add_prop('species_losses', len(diff_sp))
    n.add_prop('species_losses_percentage', perc_diff_sp)

    

