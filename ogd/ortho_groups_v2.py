from collections import defaultdict
import utils




def get_all_ogs(t, taxonomy_db, clean_name_tree):

    #   5.1. Get Monophyletic OGs
    ogs_info = defaultdict()
    seqs_in_mono_ogs, count = get_monophyletic_ogs(t, taxonomy_db, clean_name_tree, ogs_info)
 
    #   5.2 Get Paraphyletic OGs
    t, seqs_in_para_ogs, c = get_paraphyletic_ogs(t, taxonomy_db, clean_name_tree, ogs_info)

    #   5.3 Chech OGs in root
    seqs_in_root = root_ogs(t, taxonomy_db, count, c, clean_name_tree, ogs_info)
    

    #   5.3 Join OGs dictionaries
    
    # ogs_info.update(monophyletic_ogs_annot)
    # ogs_info.update(paraphyletic_ogs_annot)
    # # ogs_info.update(annot_ogs_in_root)


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



def add_entry_to_ogs_info(ogs_info, values):
    """
    Añade una nueva entrada al diccionario ogs_info con la estructura dada.
    """
    og_name, lca, sci_name_taxa, assoc_node, sp_in_og, seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap = values

    if len(sp_outliers) == 0:
        len_sp_out = '0'
    else:
        len_sp_out = len(sp_outliers)
    
    ogs_info[og_name] = {
        'TaxoLevel': lca,
        'SciName_TaxoLevel': sci_name_taxa.replace(' ', '_'),
        'AssocNode': assoc_node,
        'NumSP': len(sp_in_og),
        'NumMems': len(seqs_in),
        'Mems': list(seqs_in),
        'Lca_Dup': lca_dup,
        'Species_Outliers': sp_outliers,
        'Num_SP_Outliers': len_sp_out,
        'Inparalogs_Rate': inparalogs_rate,
        'SP_overlap_dup': species_overlap
    } 



#  Functions for monophyletic OGs
def get_monophyletic_ogs(t,  taxonomy_db, clean_name_tree, ogs_info):

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
    seqs_in_mono_ogs = annot_monophyletic_ogs(t, monophyletic_ogs, taxonomy_db, ogs_info) 

    return  seqs_in_mono_ogs, count 


def annot_monophyletic_ogs(t, monophyletic_ogs, taxonomy_db, ogs_info):


    """
        Create dict to save all the info for each monpphyletic OGs
        monophyletic_ogs = og_name: assoc_node : X ; mems: str with all mems 
        this info will be written in *.ogd_info.tsv result table with this header:
            ##OG_name   TaxoLevel   AssocNode  lensp_in_OG  OG_up   OG_down  num_OG_mems    members
    """
    
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
            
        values = [og_name, lca_subtree, sci_name_taxa, subtree_name, sp_in_og, list_seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap]
        add_entry_to_ogs_info(ogs_info, values)


        seqs_in_mono_ogs.update(set(list_seqs_in))

    return seqs_in_mono_ogs




#  Functions for paraphyletic OGs
def get_paraphyletic_ogs(t, taxonomy_db, clean_name_tree, ogs_info) :

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
    seqs_in_para_ogs = annot_paraphyletic_ogs(t, paraphyletic_ogs, taxonomy_db, ogs_info) 

    return t, seqs_in_para_ogs, c


def annot_paraphyletic_ogs(t, paraphyletic_ogs, taxonomy_db, ogs_info):

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

        lca_dup = '-'
        sp_outliers = list()
        inparalogs_rate = '-'
        species_overlap = '-'
        values = [pog_name, lca_pog, sci_name_taxa, assoc_node, sp_in_og, seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap]
        add_entry_to_ogs_info(ogs_info, values)
                
    return seqs_in_para_ogs


#   OGs in root
def root_ogs(t, taxonomy_db, count, c, clean_name_tree, ogs_info):

    """
    If root's lca is different than LUCA, 
    Recovery one OGs for each taxlev from root's lca to LUCA 
    """
    
    lca_tree = t.props.get('lca_node')

    fromluca2lcatree = list()
    if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
        fromluca2lcatree = taxonomy_db.get_lineage(lca_tree)
            
    elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy': 
        if lca_tree != 'r_root':
            fromluca2lcatree = taxonomy_db.get_name_lineage([lca_tree])[0][lca_tree]

    # Remove taxid 1 == root
    fromluca2lcatree.remove(1)
    annot_ogs_in_root = defaultdict(dict)
    seqs_in_root = set()


    #Si el lca es cell org, 
    if t.props.get("lca_node") == 131567:
        if t.props.get('monophyletic_og') == 'True':
           
            list_seqs_in = t.props.get('leaves_in')
            seqs_in_root.update(list_seqs_in)
            for taxa in fromluca2lcatree:
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


                values = [og_name, taxa, sci_name_taxa, t.name, sp_in_og, list_seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap]
                add_entry_to_ogs_info(ogs_info, values)

                
                
        #if both root child's are monophyletic OG, then there are not paraphyletic OG at cell org level 
        if t.children[0].props.get('monophyletic_og') == 'True' and t.children[1].props.get('monophyletic_og') == 'True':
            pass
        else:
            para_ogs_root = dict()
    
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
                if  len(diff)>1 and len(nodes2ckech)>0:
                        name = clean_name_tree+'@para_OG_'+str(c)
                        c+=1
                        t.add_prop('paraphyletic_og', 'True')
                        t.add_prop('pog_name', name)
                        
                        para_ogs_root[name] = [t.name, str(lca_tree), (list(diff))]
                

            seqs_paraog_in_root = annot_paraphyletic_ogs(t, para_ogs_root, taxonomy_db, ogs_info) 
            
    
    
    else:
        list_seqs_in = t.props.get('leaves_in')
        seqs_in_root.update(list_seqs_in)


        for taxa in fromluca2lcatree:
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


            values = [og_name, taxa, sci_name_taxa, t.name, sp_in_og, list_seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap]
            add_entry_to_ogs_info(ogs_info, values)
                
    return  seqs_in_root



 
                          

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



