from collections import defaultdict

# 5. Get OGs for all taxonomical levels
def get_ogs(t, level2sp_mem, taxonomy_db):

    my_descendant = get_my_descendant(level2sp_mem, taxonomy_db)
    set_trees  = set(t.search_nodes(node_is_og='True'))

    #set_trees.add(t)
    ogs = defaultdict(dict)
    count = 0

    for subtree in set_trees:

        lca_subtree = str(subtree.props.get('lca_dup'))
        if lca_subtree == 'None':
            lca_subtree = str(subtree.props.get('lca_node'))

        lin_lca_subtree = my_descendant[lca_subtree]
        taxa2remove = set()


        '''
            Si el lca_subtree es Opistok y luego hay otros OG a nivel xejem Vertebrata
            Solo quieres visitar el linaje desde Opist hasta Vertebrara
        '''
        for ogs_in_subtree in subtree.search_nodes(node_is_og='True'):
            if ogs_in_subtree.name != subtree.name:
                lca_dup = str(ogs_in_subtree.props.get('lca_dup'))
                taxa2remove.update(my_descendant[lca_dup])


        lin2check = set(lin_lca_subtree).difference(taxa2remove)
        if lca_subtree in lin2check:
            lin2check.remove(lca_subtree)


        name = 'OG_'+str(count)
        count+=1

        mems = get_members(subtree, lca_subtree)
        if len(mems) >=2:
            ogs[str(lca_subtree)][name] = (subtree.name, mems)

        for taxa in lin2check:

            if str(taxa) in level2sp_mem.keys():
                if len(list(subtree.search_nodes(lca_node=str(taxa), node_create_og='True'))) == 0:
                    name = 'OG_'+str(count)
                    count+=1
                    mems = get_members(subtree, taxa)
                    if len(mems.split('|')) >=2:
                        ogs[str(taxa)][name] = (subtree.name, (mems))

    base_ogs = from_lca_tree2cellorg(t, ogs, taxonomy_db, count)

    return base_ogs

def get_my_descendant(level2sp_mem, taxonomy_db):

    my_descendant = defaultdict(set)

    for taxa in level2sp_mem.keys():
        total_ancest = taxonomy_db.get_lineage(taxa)

        for ancest in total_ancest:
            if str(ancest) in level2sp_mem.keys():
                my_descendant[str(ancest)].add(taxa)

    return my_descendant

def get_members(node, taxa):

    all_leafs = node.props.get('_leaves_in_nodes')

    mems_set = set()
    for l in all_leafs:
        if str(taxa) in l.props.get('lineage').split('|'):
            mems_set.add(l.name)

    mems = '|'.join(list(mems_set))

    return mems

def from_lca_tree2cellorg(t, base_ogs, taxonomy_db, count):

    lca_tree = t.props.get('lca_node')
    if '131567' not in base_ogs:
        taxa2add = taxonomy_db.get_lineage(lca_tree)
        mems = get_members(t, str(lca_tree))
        for taxa in taxa2add:
            if str(taxa) not in base_ogs.keys():
                count +=1
                name='OG_'+str(count)
                base_ogs[str(taxa)][name] = (t.name, mems)

    return base_ogs


