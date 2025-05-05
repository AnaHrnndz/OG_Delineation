from collections import defaultdict
import ogd.utils as utils



def get_all_ogs(t, taxonomy_db):

    #   5.1. Get Monophyletic OGs
    ogs_info = defaultdict()

    add_monophyletic_ogs(t, taxonomy_db, ogs_info)
    
    #   5.2 Get Paraphyletic OGs
    add_paraphyletic_ogs(t, taxonomy_db, ogs_info)
    
    #   5.3 Chech OGs in root
    root_ogs(t, taxonomy_db, ogs_info)
    
    #   5.4 Get hierarchical structure of OGs and add info to the tree and to the ogs_info dict
    hierarchy(ogs_info, taxonomy_db)
    

    #   5.5 Get all seqs in OGs 
    seqs_in_ogs = set()
    for og, info in ogs_info.items():
        seqs_in_ogs.update(set(info['Mems']))
        ogs_info[og]['NumRecoverySeqs'] = '0'
        ogs_info[og]['RecoverySeqs'] = list()
    
    return t, ogs_info, seqs_in_ogs




#  Functions for monophyletic OGs
def add_monophyletic_ogs(t,  taxonomy_db, ogs_info):

    """
        Busca todos los nodos internos marcados como monophyletic_og para anotarlos
    """

    set_trees  = set(t.search_nodes(monophyletic_og='True'))
    
    monophyletic_ogs = defaultdict()
    
    for subtree in set_trees:
        
        if subtree.is_root:
            pass
        else: 
            og_name = str(subtree.props.get('lca_dup'))+'|'+subtree.name

            subtree.add_prop('og_name', og_name)
            monophyletic_ogs[og_name] = subtree.name          

    # Annot nodes that are monophyletic OGs
    annot_monophyletic_ogs(t, monophyletic_ogs, taxonomy_db, ogs_info) 

    

def annot_monophyletic_ogs(t, monophyletic_ogs, taxonomy_db, ogs_info):


    """
        Create dict to save all the info for each monpphyletic OGs
        monophyletic_ogs = og_name: assoc_node : X ; mems: str with all mems 
        this info will be written in *.ogd_info.tsv result table with this header:
            ##OG_name   TaxoLevel   AssocNode  lensp_in_OG  OG_up   OG_down  num_OG_mems    members
    """
    
    for og_name, node_name in monophyletic_ogs.items():

        subtree_name = node_name

        lca_subtree = str(t[subtree_name].props.get('lca_node'))
        sci_name_lca = utils.get_sci_name(taxonomy_db, lca_subtree)
        
        list_seqs_in = t[subtree_name].props.get('leaves_in')
        sp_in_og = t[subtree_name].props.get('sp_in')
        lca_dup = t[subtree_name].props.get('lca_dup')
        sp_outliers = t[subtree_name].props.get('sp_out')
        inparalogs_rate = t[subtree_name].props.get('inparalogs_rate')
        species_overlap = t[subtree_name].props.get('so_score_dup')
            
        values = [og_name, lca_subtree, sci_name_lca, subtree_name, sp_in_og, list_seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap]
                  
        add_entry_to_ogs_info(ogs_info, values)

 


#  Functions for paraphyletic OGs
def add_paraphyletic_ogs(t, taxonomy_db, ogs_info) :

    """
    Recuperar las seqs que se hayan podido "escapar" entre OGs
    Ignoramos los nodos que son OGs,
    para los nodos que no son OGs, hago un set con todas las hojas (excepto outliers)
    busco por debajo de este nodo si hay OGs, si los hay hago otro set con las seqs que si estan en OG
    hago la diferencia de los sets y asi saco los paraphyletic_ogs
    """

    paraphyletic_ogs = defaultdict(dict)
    
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
                            
                            mem2remove = set()
                            nodes2ckech = set(ch.search_nodes(monophyletic_og='True', lca_dup=lca_target))
                            nodes2ckech.update(ch.search_nodes(evoltype_2='D', lca_node=lca_target))
                            for n_ in nodes2ckech:
                                mem2remove.update(n_.props.get('leaves_in'))

                            # Remove seqs that belong to some OGs
                            diff = all_mems.difference(mem2remove)
                            
                            if  len(diff)>1:
                                og_name = str(ch.props.get('lca_node'))+'|'+ch.name+'!'
                                ch.add_prop('paraphyletic_og', 'True')
                                ch.add_prop('og_name', og_name)
                                paraphyletic_ogs[og_name] = [ch.name, list(diff)]

    # Annot all the paraphyletic OGs
    annot_paraphyletic_ogs(t, paraphyletic_ogs, taxonomy_db, ogs_info) 



def annot_paraphyletic_ogs(t, paraphyletic_ogs, taxonomy_db, ogs_info):

    for pog_name, info in paraphyletic_ogs.items():
        
        assoc_node_name = info[0]
        mems2check = info[1]
       
        pog_node = t[assoc_node_name]
        
        lca_pog = pog_node.props.get('lca_node')

        sci_name_lca = utils.get_sci_name(taxonomy_db, lca_pog)

        sp_in_og = set()
        list_seqs_in = list()
        sp_outliers = list()
    
        for lname in mems2check:
            taxid = t[lname].props.get('taxid')

            # Check that seq is not an taxo_outlier or long branch
            if lname not in pog_node.props.get('leaves_out'):
                
                list_seqs_in.append(lname)
                sp_in_og.add(taxid)
                
                # Add to the leaves a prop with all the paraphyletic OGs
                old_pOG = pog_node[lname].props.get('pOG', str())
                new_pOG = old_pOG+'@'+pog_name
                pog_node[lname].add_prop('pOG', new_pOG)
            
            elif lname not in pog_node.props.get('leaves_out'):
                sp_outliers.add(taxid)
            

        lca_dup = '-'
        inparalogs_rate = pog_node.props.get('inparalogs_rate')
        species_overlap = pog_node.props.get('so_score')
        
        values = [pog_name, lca_pog, sci_name_lca, assoc_node_name, sp_in_og, list_seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap]
                  
        add_entry_to_ogs_info(ogs_info, values)
                
    return 


#   OGs in root
def root_ogs(t, taxonomy_db, ogs_info):

    """
    If root's lca is different than LUCA, 
    Recovery one OGs for each taxlev from root's lca to LUCA 
    """

    lca_tree = t.props.get('lca_node')

    fromluca2lcatree = utils.get_lineage(taxonomy_db, lca_tree)
    
    # Remove taxid 1 == root
    fromluca2lcatree.remove(1)
    fromluca2lcatree.reverse()
   
    
    # If lca's root is not luca, then retrieve one OG per taxa level from lca's root to luca
    if lca_tree not in  [131567, 'r_root']:
        list_seqs_in = list(t.props.get('leaves_in'))

        # Avoid repeat lca's root
        fromluca2lcatree.remove(lca_tree)
       
        for taxa in fromluca2lcatree:
            
            sci_name_taxa = utils.get_sci_name(taxonomy_db, taxa)
            
            og_name =str(taxa)+'|'+t.name+'*'
            
            sp_in_og = list(t.props.get('sp_in'))
            lca_dup = '-'
            sp_outliers = t.props.get('sp_out')
            inparalogs_rate = t.props.get('inparalogs_rate')
            species_overlap = t.props.get('so_score')
            
            values = [og_name, taxa, sci_name_taxa, t.name, sp_in_og, list_seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap]
            add_entry_to_ogs_info(ogs_info, values)
        
        

        
    
    if t.props.get('monophyletic_og') == 'True':
       
        list_seqs_in = t.props.get('leaves_in')
        lca_root = t.props.get("lca_node")
        sci_name_lca = utils.get_sci_name(taxonomy_db, lca_root)
       
        og_name = str(lca_root)+'|'+t.name
            
        sp_in_og = list(t.props.get('sp_in'))
        lca_dup = '-'
        sp_outliers = t.props.get('sp_out')
        inparalogs_rate = t.props.get('inparalogs_rate')
        species_overlap = t.props.get('so_score')
        
        values = [og_name, lca_root, sci_name_lca, t.name, sp_in_og, list_seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap]
        add_entry_to_ogs_info(ogs_info, values)

        
    

    else:
       
        if t.children[0].props.get('monophyletic_og') == 'True' and t.children[1].props.get('monophyletic_og') == 'True':
            pass
       
        else:

            
            for ch in t.children:
                if ch.props.get('monophyletic_og') == 'True':
                    pass
                else:
                    lca_ch = ch.props.get('lca_node')
                    
                    sci_name_taxa = utils.get_sci_name(taxonomy_db, lca_ch)
                    
                    all_mems = ch.props.get('leaves_in', set())
                    mem2remove = set()
                    nodes2ckech = set(ch.search_nodes(monophyletic_og='True', lca_dup=lca_ch))
                    nodes2ckech.update(ch.search_nodes(evoltype_2='D', lca_node=lca_ch))
                    for n_ in nodes2ckech:
                        mem2remove.update(n_.props.get('leaves_in'))
                    
                    list_seqs_in = all_mems.difference(mem2remove)
                    
                    if  len(list_seqs_in)>1:
                       
                        og_name = str(lca_ch)+'|'+ch.name+'!'
                
                        ch.add_prop('paraphyletic_og', 'True')
                        ch.add_prop('og_name', og_name)
                        lca_dup = '-'
                        sp_outliers = list()
                        inparalogs_rate = ch.props.get('inparalogs_rate')
                        species_overlap = ch.props.get('so_score')
                        sp_in_og = set()
                        for l in list_seqs_in:
                            taxid = t[l].props.get('taxid')
                            sp_in_og.add(taxid)
                            old_pOG = t[l].props.get('pOG', str())
                            new_pOG = old_pOG+'@'+og_name
                            t[l].add_prop('pOG', new_pOG)
                        
                        values = [og_name, lca_ch, sci_name_taxa, ch.name, list(sp_in_og), list_seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap]
                        add_entry_to_ogs_info(ogs_info, values)

                       

    
    




def add_entry_to_ogs_info(ogs_info, values):

    """
    Añade una nueva entrada al diccionario ogs_info con la estructura dada.
    """
    og_name, lca, sci_name_lca, assoc_node_name, sp_in_og, list_seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap = values

    if len(sp_outliers) == 0:
        len_sp_out = '0'
        sp_outliers = ['-']
    else:
        len_sp_out = len(sp_outliers)
    
    ogs_info[og_name] = {
        'TaxoLevel': lca,
        'SciName_TaxoLevel': sci_name_lca.replace(' ', '_'),
        'AssocNode': assoc_node_name,
        'NumSP': len(sp_in_og),
        'NumMems': len(list_seqs_in),
        'Mems': list(list_seqs_in),
        'Lca_Dup': lca_dup,
        'Species_Outliers': sp_outliers,
        'Num_SP_Outliers': len_sp_out,
        'Inparalogs_Rate': inparalogs_rate,
        'SP_overlap_dup': species_overlap
    } 
                          


def hierarchy(ogs_info, taxonomy_db):

    for og_name, info in ogs_info.items():
        ogs_up = set()
        ogs_down = set()
        target_mems = set(info['Mems'])
        og_level = int(og_name.split('|')[1].split('-')[1].replace('*', '').replace('!', ''))
        tax_level = int(og_name.split('|')[0])
        for o, i in ogs_info.items():
            if len(target_mems&set(i['Mems'])) > 0:
                
                o_lev = int(o.split('|')[1].split('-')[1].replace('*', '').replace('!', ''))
                t_lev = int(o.split('|')[0])

                if og_level < o_lev:
                    ogs_down.add(o)
                elif og_level > o_lev:
                    ogs_up.add(o)
                else:
                    if tax_level == t_lev:
                        pass
                    else:
                        o_lin = utils.get_lineage(taxonomy_db, t_lev)
                        
                        if tax_level not in o_lin:
                            ogs_up.add(o)
                        elif tax_level in o_lin:
                            ogs_down.add(o)

        if len(ogs_down) == 0:
            ogs_down = '-'

        if len(ogs_up) == 0:
            ogs_up = '-'
                        
        info['OG_down'] = ogs_down
        info['OG_up'] = ogs_up

       
            




