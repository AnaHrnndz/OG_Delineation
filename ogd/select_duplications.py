import utils



##  4. Detect Duplications and Core-OGs & PGs ##

def run_get_main_dups(t, taxonomy_db, total_mems_in_tree, args):
    """
        Select high-quality duplication nodes:
            - Species overlap higher than min threshold
            - More than one leave and more than one specie
            - No more duplication nodes at the same taxonomic level below the node
    """

    # taxid_dups_og, set with the lca of the nodes that are OGs
    taxid_dups_og = set()
    

    if args.so_euk == None:
        so_euk = 'None'
    else:
        so_euk = args.so_euk

    if args.so_bact == None:
        so_bact = 'None'
    else:
        so_bact = args.so_bact
    
    if args.so_arq == None:
        so_arq = 'None'
    else:
        so_arq = args.so_arq

    mssg = f"""
    3. Select high quality duplication nodes
       -Species overlap threshold:
            General: {args.so_all}
            Euk: {so_euk}
            Bact: {so_bact}
            Arq: {so_arq}"""
    print(mssg)
    

    # Traverse tree to find the nodes that are "good" duplications and generate OGs.
    for node in t.traverse("preorder"):

        #if node.is_root:
        #    pass

        if not node.is_leaf and node.props.get('evoltype_2') == 'D'  \
        and len(node.props.get('leaves_in')) >1 and len(node.props.get('sp_in')) > 1 :

            dups_under_node = []
            for n in node.search_nodes(evoltype_2='D'):

                if n.name!= node.name:
                    dups_under_node.append(n)

            # There are more dups under the node
            if len(dups_under_node) > 0:

                lca_target = node.props.get('lca_node')

                #Save Dups under child1 and child2 that have the same lca_node
                dups_under_ch1 = list(node.children[0].search_nodes(evoltype_2='D', lca_node=lca_target))
                dups_under_ch2 = list(node.children[1].search_nodes(evoltype_2='D', lca_node=lca_target))

                save_dups_ch1 = 0 
                save_dups_ch2 = 0 

                # Check that dups under child1 and child2 (that have the same lca) fit all requirements : species overlap min requirement,
                # more than 1 leaves and more than 1 species

                for n_ in  dups_under_ch1:
                    if len(n_.props.get('leaves_in')) >1 and len(n_.props.get('sp_in'))> 1 :
                        save_dups_ch1 += 1

                for n_ in  dups_under_ch2:
                    if  len(n_.props.get('leaves_in'))>1 and len(n_.props.get('sp_in'))> 1 :
                        save_dups_ch2 += 1

                if save_dups_ch1 == 0:
                    annotate_dups_ch(taxid_dups_og, node, 'ch1', taxonomy_db)

                if save_dups_ch2 == 0 :
                    annotate_dups_ch(taxid_dups_og ,node, 'ch2', taxonomy_db)


            elif len(dups_under_node) == 0:
                annotate_dups_ch(taxid_dups_og, node, 'ch1', taxonomy_db)
                annotate_dups_ch(taxid_dups_og, node, 'ch2', taxonomy_db)


    #  Now decide if root is OG or not,
    #  if there are not dups at the same level that root's lca
    #  then root is a monophyletic_og 
    if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
        lca_root = int(t.props.get('lca_node'))
    elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
        lca_root = (t.props.get('lca_node'))

    #if len(list(t.search_nodes(evoltype_2='D', lca_node=lca_root))) == 0:
    if len(list(t.search_nodes(monophyletic_og='True', lca_node=lca_root))) == 0:
        t.add_prop('monophyletic_og', 'True')

    
    t, props = utils.run_clean_properties(t)

    return t, taxid_dups_og




def annotate_dups_ch(taxid_dups_og, node, ch_node, taxonomy_db):

    """
    Add props and save info about the children node that is OG
    node is the duplication node
    target node is the node that is OG
    node = parent
    target node = child
    """

    if ch_node == 'ch1':
        og_name_ch = node.props.get('ch1_name')
        og_ch_mems = node.props.get('leaves_ch1')
        sp_ch = node.props.get('sp_in_ch1')
        target_node = next(node.search_nodes(name=og_name_ch))

    elif ch_node == 'ch2':
        og_name_ch = node.props.get('ch2_name')
        og_ch_mems = node.props.get('leaves_ch2')
        sp_ch = node.props.get('sp_in_ch2')
        target_node = next(node.search_nodes(name=og_name_ch))

    if  len(sp_ch) > 1 and len(og_ch_mems) > 1:
        target_node.add_prop('monophyletic_og', 'True')
        target_node.add_prop('lca_dup', node.props.get('lca_node'))
        target_node.add_prop('so_score_dup', node.props.get('so_score'))

        
        #if node.props.get('lca_node') != 'r_root':
        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            target_node.add_prop('dup_lineage', list(taxonomy_db.get_lineage(node.props.get('lca_node'))))
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
            target_node.add_prop('dup_lineage', list(taxonomy_db.get_name_lineage([node.props.get('lca_node')])[0][node.props.get('lca_node')]))
        # else:
            # target_node.add_prop('dup_lineage', ['r_root'])

        target_node.add_prop('dup_node_name', node.props.get('name'))

        taxid_dups_og.add(node.props.get('lca_node'))
        node.add_prop('node_create_og', 'True')
        




    # """
        # Find all OGs created under root node
    # """
    # ogs_down = set()
    # dups_down = list()
    # for n in t.search_nodes(monophyletic_og="True"):
        # if n.name != t.name:
            # ogs_down.add(n.name)
            # dups_down.append(n.up.name)

    # ogs_down_value = ogs_down if len(ogs_down) > 0 else '-'
    # t.add_prop('ogs_down', ogs_down_value)
    
    # dups_down_value = dups_down if len(dups_down) > 0 else '-'
    # t.add_prop('dups_down', dups_down_value)

    

