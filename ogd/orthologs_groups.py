from collections import defaultdict
import utils




def get_all_ogs(t, taxonomy_db, clean_name_tree):

    #   5.1. Get Monophyletic OGs
    monophyletic_ogs_annot, seqs_in_mono_ogs, count = get_monophyletic_ogs(t, taxonomy_db, clean_name_tree)
 
    #   5.2 Get Paraphyletic OGs
    t,  paraphyletic_ogs_annot, seqs_in_para_ogs, c = get_paraphyletic_ogs(t, taxonomy_db, clean_name_tree)

    #   5.3 Chech OGs in root
    annot_ogs_in_root, seqs_in_root = root_ogs(t, taxonomy_db, count, c, clean_name_tree)
    

    #   5.3 Join OGs dictionaries
    ogs_info = defaultdict()
    ogs_info.update(monophyletic_ogs_annot)
    ogs_info.update(paraphyletic_ogs_annot)
    ogs_info.update(annot_ogs_in_root)


    #   5.4 Get hierarchical structure of OGs (mono and paraphyletic)
    #       and add info to the tree and to the ogs_info dict
    t, ogs_info = hierarchy_ogs(t, ogs_info)


    #   5.5 Join all seqs in OGs (mono and para phyletic)
    seqs_in_ogs = set()
    seqs_in_ogs.update(seqs_in_mono_ogs)
    seqs_in_ogs.update(seqs_in_para_ogs)
    seqs_in_ogs.update(seqs_in_root)


    for og, info in ogs_info.items():
        ogs_info[og]['NumRecoverySeqs'] = '0'
        ogs_info[og]['RecoverySeqs'] = list()

    
    return t, ogs_info, seqs_in_ogs






#  Functions for monophyletic OGs
def get_monophyletic_ogs(t,  taxonomy_db, clean_name_tree):

    """
        Return a dict with the info for all monophyletic OGs 
    """

    set_trees  = set(t.search_nodes(monophyletic_og='True'))
    

    monophyletic_ogs = defaultdict(dict)
    count = 0

    for subtree in set_trees:
        
        if subtree.is_root:
            pass
        else:
            all_leaves_subtree = subtree.props.get('leaves_in')

            if subtree.props.get('lca_dup') != None:
                lca_subtree = str(subtree.props.get('lca_dup'))

            elif  subtree.props.get('lca_dup') == None:
                lca_subtree = str(subtree.props.get('lca_node'))

            name = clean_name_tree+'@mono_OG_'+str(count)
            count+=1

            if len(all_leaves_subtree) >=1:
                subtree.add_prop('mog_name', name)
                monophyletic_ogs[name] = [subtree.name, str(lca_subtree), (list(all_leaves_subtree))]          


    # Annot all monophyletic OGs
    monophyletic_ogs_annot, seqs_in_mono_ogs = annot_monophyletic_ogs(t, monophyletic_ogs, taxonomy_db) 

    return monophyletic_ogs_annot, seqs_in_mono_ogs, count 


def annot_monophyletic_ogs(t, monophyletic_ogs, taxonomy_db):


    """
        Create dict to save all the info for each monpphyletic OGs
        monophyletic_ogs = og_name: assoc_node : X ; mems: str with all mems 
        this info will be written in *.ogd_info.tsv result table with this header:
            ##OG_name   TaxoLevel   AssocNode  lensp_in_OG  OG_up   OG_down  num_OG_mems    members
    """
    
    annot_mono_ogs = defaultdict(dict)
    seqs_in_mono_ogs = set()
    
    for og_name, og_info in monophyletic_ogs.items():

        subtree_name = og_info[0]
        lca_subtree = og_info[1]
        list_seqs_in = og_info[2]

        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            sci_name_taxa =  taxonomy_db.get_taxid_translator([lca_subtree])[int(lca_subtree)]
    
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
            sci_name_taxa = lca_subtree
        
        
        sp_in_og = t[subtree_name].props.get('sp_in')
        lca_dup = t[subtree_name].props.get('lca_dup')
        sp_outliers = t[subtree_name].props.get('sp_out')
        len_sp_outliers = len(sp_outliers)
        if len_sp_outliers == 0:
                sp_outliers = '-'
        inparalogs_rate = t[subtree_name].props.get('inparalogs_rate')
        species_overlap = t[subtree_name].props.get('so_score_dup')
            
        annot_mono_ogs[og_name]['TaxoLevel'] = lca_subtree
        annot_mono_ogs[og_name]['SciName_TaxoLevel'] = sci_name_taxa.replace(' ', '_')
        annot_mono_ogs[og_name]['AssocNode'] = subtree_name
        annot_mono_ogs[og_name]['NumSP'] = len(sp_in_og)
        annot_mono_ogs[og_name]['NumMems'] = len(list_seqs_in)
        annot_mono_ogs[og_name]['Mems'] = list_seqs_in
        annot_mono_ogs[og_name]['Lca_Dup'] = lca_dup
        annot_mono_ogs[og_name]['Species_Outliers'] = sp_outliers
        annot_mono_ogs[og_name]['Num_SP_Outliers'] = len_sp_outliers
        annot_mono_ogs[og_name]['Inparalogs_Rate'] = inparalogs_rate
        annot_mono_ogs[og_name]['SP_overlap_dup'] = species_overlap

        seqs_in_mono_ogs.update(set(list_seqs_in))

    return annot_mono_ogs, seqs_in_mono_ogs



#  Functions for paraphyletic OGs
def get_paraphyletic_ogs(t, taxonomy_db, clean_name_tree):

    """
    Recuperar las seqs que se hayan podido "escapar" entre OGs
    Ignoramos los nodos que son OGs,
    para los nodos que no son OGs, hago un set con todas las hojas (excepto outliers)
    busco por debajo de este nodo si hay OGs, si los hay hago otro set con las seqs que si estan en OG
    hago la diferencia de los sets y asi saco los paraphyletic_ogs
    """

    paraphyletic_ogs = defaultdict(dict)
    c=0
    
    for n in t.traverse():

        lca_target = n.props.get('lca_node')

        if not n.is_root:

            if n.props.get('evoltype_2')=='D':

                if n.children[0].props.get('monophyletic_og') == 'True' and n.children[1].props.get('monophyletic_og') == 'True':
                    pass

                else:

                    for ch in n.children:
                        if ch.props.get('monophyletic_og') == 'True':
                            pass
                        else:
                            all_mems = ch.props.get('leaves_in', set())
                            lca_ch = ch.props.get('lca_node')
                            mem2remove = set()
                            nodes2ckech = set(ch.search_nodes(monophyletic_og='True', lca_dup=lca_target))
                            nodes2ckech.update(ch.search_nodes(evoltype_2='D', lca_node=lca_target))
                            for n_ in nodes2ckech:
                                mem2remove.update(n_.props.get('leaves_in'))

                            diff = all_mems.difference(mem2remove)
                            if  len(diff)>1:

                                name = clean_name_tree+'@para_OG_'+str(c)
                                c+=1
                                ch.add_prop('paraphyletic_og', 'True')
                                ch.add_prop('pog_name', name)

                                paraphyletic_ogs[name] = [ch.name, str(lca_ch), (list(diff))]

    # Annot all the paraphyletic OGs
    paraphyletic_ogs_annot, seqs_in_para_ogs = annot_paraphyletic_ogs(t, paraphyletic_ogs, taxonomy_db) 

    return t,  paraphyletic_ogs_annot, seqs_in_para_ogs, c


def annot_paraphyletic_ogs(t, paraphyletic_ogs, taxonomy_db):

    annot_para_ogs = defaultdict(dict)
    seqs_in_para_ogs = set()

    for pog_name, info in paraphyletic_ogs.items():
        assoc_node = info[0]
        lca_pog = int(info[1])
        mems2check = info[2]

        pog_node = t[assoc_node]

        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            sci_name_taxa= taxonomy_db.get_taxid_translator([lca_pog])[lca_pog]
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
            sci_name_taxa = lca_pog

        sp_in_og = set()
        ogs_down = set()
        ogs_up = set()
        seqs_in = set()
    
        for lname in mems2check:
            if lname not in pog_node.props.get('leaves_out'):
                
                seqs_in.add(lname)
                seqs_in_para_ogs.add(lname)
                sp_in_og.add(lname.split('.')[0])
                
                # Add to the leaves a prop with all the paraphyletic OGs
                old_pOG = pog_node[lname].props.get('pOG', str())
                new_pOG = old_pOG+'@'+pog_name
                pog_node[lname].add_prop('pOG', new_pOG)
                
        annot_para_ogs[pog_name]['TaxoLevel'] = lca_pog
        annot_para_ogs[pog_name]['SciName_TaxoLevel'] = sci_name_taxa.replace(' ', '_')
        annot_para_ogs[pog_name]['AssocNode'] = assoc_node
        annot_para_ogs[pog_name]['NumSP'] = len(sp_in_og)
        annot_para_ogs[pog_name]['NumMems'] = len(seqs_in)
        annot_para_ogs[pog_name]['Mems'] = list(seqs_in)
        annot_para_ogs[pog_name]['Lca_Dup'] = '-'
        annot_para_ogs[pog_name]['Species_Outliers'] = '-'
        annot_para_ogs[pog_name]['Num_SP_Outliers'] = '0'
        annot_para_ogs[pog_name]['Inparalogs_Rate'] = '-'
        annot_para_ogs[pog_name]['SP_overlap_dup'] = '-'

    return annot_para_ogs, seqs_in_para_ogs


#   OGs in root
def root_ogs(t, taxonomy_db, count, c, clean_name_tree):

    """
    If root's lca is different than LUCA, 
    Recovery one OGs for each taxlev from root's lca to LUCA 
    """
    
    lca_tree = t.props.get('lca_node')

    annot_ogs_in_root = defaultdict(dict)
    seqs_in_root = set()

    # First chech if root is monophyletic og
    if 'monophyletic_og' in t.props.keys():
        
        # Return one OG per level from root's lca to LUCA 
        taxa2add = list()
        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            taxa2add = taxonomy_db.get_lineage(lca_tree)
            
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy': 
            if lca_tree != 'r_root':
                taxa2add = taxonomy_db.get_name_lineage([lca_tree])[0][lca_tree]
                
        
        list_seqs_in = t.props.get('leaves_in')
        seqs_in_root.update(list_seqs_in)
        for taxa in taxa2add:
            if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
                sci_name_taxa =  taxonomy_db.get_taxid_translator([taxa])[int(taxa)]
            elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy': 
                sci_name_taxa = taxa
            

            count +=1
            og_name = clean_name_tree+'@mono_OG_'+str(count)

            if taxa != lca_tree:
                og_name = og_name+'*'
            
            sp_in_og = t.props.get('sp_in')
            lca_dup = '-'
            sp_outliers = t.props.get('sp_out', '-')
            len_sp_outliers = len(sp_outliers)
            if len_sp_outliers == 0:
                sp_outliers = '-'
            inparalogs_rate = t.props.get('inparalogs_rate')
            species_overlap = t.props.get('so_score')
            
            annot_ogs_in_root[og_name]['TaxoLevel'] = taxa
            annot_ogs_in_root[og_name]['SciName_TaxoLevel'] = sci_name_taxa.replace(' ', '_')
            annot_ogs_in_root[og_name]['AssocNode'] = t.name
            annot_ogs_in_root[og_name]['NumSP'] = len(sp_in_og)
            annot_ogs_in_root[og_name]['NumMems'] = len(list_seqs_in)
            annot_ogs_in_root[og_name]['Mems'] = list_seqs_in
            annot_ogs_in_root[og_name]['Lca_Dup'] = lca_dup
            annot_ogs_in_root[og_name]['Species_Outliers'] = sp_outliers
            annot_ogs_in_root[og_name]['Num_SP_Outliers'] = len_sp_outliers
            annot_ogs_in_root[og_name]['Inparalogs_Rate'] = inparalogs_rate
            annot_ogs_in_root[og_name]['SP_overlap_dup'] = species_overlap

    
    else:
        para_ogs_root = dict()
        # Check if root child's are monophyletic OG, if both childs are OGs, then there are not any paraphyletic OG at root level 
        if t.children[0].props.get('monophyletic_og') == 'True' and t.children[1].props.get('monophyletic_og') == 'True':
            pass
        else:
            # If root is Duplication, then search in each child for paraphyletic OGs independently
            if t.props.get('evoltype_2') == 'D':
                for ch in t.children:
                    if ch.props.get('monophyletic_og') == 'True':
                        pass
                    else:
                        lca_ch = ch.props.get('lca_node')
                        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
                            sci_name_taxa =  taxonomy_db.get_taxid_translator([lca_ch])[int(lca_ch)]
                        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy': 
                            sci_name_taxa = lca_ch
                        

                        all_mems = ch.props.get('leaves_in', set())
                        mem2remove = set()
                        nodes2ckech = set(ch.search_nodes(monophyletic_og='True', lca_dup=lca_ch))
                        nodes2ckech.update(ch.search_nodes(evoltype_2='D', lca_node=lca_ch))
                        for n_ in nodes2ckech:
                            mem2remove.update(n_.props.get('leaves_in'))
                        diff = all_mems.difference(mem2remove)
                        if  len(diff)>1:
                                name = clean_name_tree+'@para_OG_'+str(c)
                                c+=1
                                ch.add_prop('paraphyletic_og', 'True')
                                ch.add_prop('pog_name', name)

                                para_ogs_root[name] = [ch.name, str(lca_ch), (list(diff))]
            
            # If root is no Duplication, then search for one paraphyletic OGs 
            elif t.props.get('evoltype_2') != 'D':
                    all_mems = t.props.get('leaves_in', set())
                    mem2remove = set()
                    nodes2ckech = set(t.search_nodes(monophyletic_og='True', lca_dup=lca_tree))
                    nodes2ckech.update(t.search_nodes(evoltype_2='D', lca_node=lca_tree))
                    for n_ in nodes2ckech:
                        mem2remove.update(n_.props.get('leaves_in'))
                    diff = all_mems.difference(mem2remove)
                    if  len(diff)>1:
                            name = clean_name_tree+'@para_OG_'+str(c)
                            c+=1
                            t.add_prop('paraphyletic_og', 'True')
                            t.add_prop('pog_name', name)
                            
                            para_ogs_root[name] = [t.name, str(lca_tree), (list(diff))]

        annot_ogs_in_root, seqs_in_root = annot_paraphyletic_ogs(t, para_ogs_root, taxonomy_db) 
    
    return annot_ogs_in_root, seqs_in_root



def hierarchy_ogs(t, ogs_info):

    """
       Traverse the tree finding monophyletic and paraphyeltic OGs,
       add the info to the tree and to the ogs_info dict 
       return the tree and the ogs_info with the info
    """

    # First annote the hierarchy in the tree
    for node in t.traverse('preorder'):
        '''
        For root node, search all OGs 
        '''
        if node.is_root:
            ogs_down = set()
            dups_down = list()
            lca_node = node.props.get('lca_node')
            # Search for monophyletic OGs
            for n in node.search_nodes(monophyletic_og="True"):
                if n.name != t.name:
                    #ogs_down.add(n.name)
                    ogs_down.add(n.props.get('mog_name'))
                    dups_down.append(n.up.name)
            
            # Search for paraphyletic OGs
            for n in node.search_nodes(paraphyletic_og="True"):
                if n.props.get('lca_node') != lca_node and n.name != t.name:
                #if n.name != t.name and n:
                    ogs_down.add(n.props.get('pog_name'))
                

            ogs_down_value = ogs_down if len(ogs_down) > 0 else ['-']
            node.add_prop('ogs_down', ogs_down_value)
    
            dups_down_value = dups_down if len(dups_down) > 0 else ['-']
            node.add_prop('dups_down', dups_down_value)
            
        else:
            '''
            Only for internal root that are OGs, search for OGs up and below
            '''
            if node.props.get('monophyletic_og') or node.props.get('paraphyletic_og'):
                node_name = node.name

                #Detect OGs below node
                ogs_down = set()
                dups_down = list()
                # Search for monophyletic OGs 
                for n in node.search_nodes(monophyletic_og="True"):
                    if node_name != n.name:
                        #ogs_down.add(n.name)
                        ogs_down.add(n.props.get('mog_name'))
                        dups_down.append(n.up.name)
                # Search for paraphyletic OGs 
                for n in node.search_nodes(paraphyletic_og="True"):
                    if n.props.get('lca_node') != lca_node and n.name != t.name and n:
                    #if node_name != n.name:
                        ogs_down.add(n.props.get('pog_name'))
                        

                ogs_down_value = ogs_down if len(ogs_down) > 0 else ['-']
                node.add_prop('ogs_down', ogs_down_value)

                dups_down_value = dups_down if len(dups_down) > 0 else ['-']
                node.add_prop('dups_down', dups_down_value)


                #Detect OGs up
                # Search for monophyletic OGs 
                ogs_up = set()
                dups_up = list()
                ogs_up, dups_up = utils.check_nodes_up(node)
                
               
                # Search for paraphyletic OGs 
                para_ogs_up = check_pogs_up(node)
                ogs_up.update(para_ogs_up)

               
                ogs_up_value = ogs_up if len(ogs_up) > 0 else ['-']
                node.add_prop('ogs_up', ogs_up_value)
             
                dups_up_value = dups_up if len(dups_up) > 0 else ['-']
                node.add_prop('dups_up', dups_up_value)


    # Then annote the hierarchy in the ogs_info dict
    for og_name, info in ogs_info.items():
        
        assoc_node = info['AssocNode']
        ogs_down = t[assoc_node].props.get('ogs_down', ['-'])
        ogs_up = t[assoc_node].props.get('ogs_up', ['-'])

        info['OG_down'] = ogs_down
        info['OG_up'] = ogs_up


    return t, ogs_info


def check_pogs_up(node):

    """
        Find OGs in upper nodes
    """

    pogs_up = set()
    lca_included = set()
    lca_node = node.props.get('lca_node')
    
    while node.up:

        if node.up.props.get('paraphyletic_og') and node.up.props.get('lca_node') != lca_node:
            if node.up.props.get('is_root'):
                if node.up.props.get('lca_node') not in lca_included:
                    pogs_up.add(node.up.props.get('pog_name'))
                    lca_included.add(node.up.props.get('lca_node'))
            else:
                pogs_up.add(node.up.props.get('pog_name'))
                lca_included.add(node.up.props.get('lca_node'))
                
        node = node.up

    return pogs_up



