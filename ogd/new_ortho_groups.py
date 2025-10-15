from collections import defaultdict
import ogd.utils as utils
from ete4 import PhyloTree
from typing import Any, Tuple, Set, List, Dict, Optional

# --- Main  Function ---

def get_all_ogs(t, taxonomy_db, og_indexes) :
    """
    Delineates all OGs (Monophyletic, Paraphyletic, and Root-based) and determines their hierarchy.
    Utilizes pre-calculated hq-dups for faster lookups.
    """
    ogs_info = defaultdict(dict)

    # 5.1. Get Monophyletic OGs (Uses annotations from Step 4)
    add_monophyletic_ogs(t, taxonomy_db, ogs_info)
    
    # 5.2 Get Paraphyletic OGs (Optimized using og_indexes)
    add_paraphyletic_ogs(t, taxonomy_db, ogs_info, og_indexes)
    
    # 5.3 Check OGs in root (Optimized using og_indexes)
    root_ogs(t, taxonomy_db, ogs_info, og_indexes)
    
    # 5.4 Get hierarchical structure of OGs (Optimized)
    hierarchy(ogs_info, taxonomy_db)
    

    # 5.5 Get all seqs in OGs and clean up
    seqs_in_ogs = set()
    final_ogs_info = dict(ogs_info) # Convert to standard dict for return
    
    for og, info in final_ogs_info.items():
        # Ensure Mems is iterable
        mems = info.get('Mems', [])
        seqs_in_ogs.update(set(mems))
        
        info['NumRecoverySeqs'] = '0'
        info['RecoverySeqs'] = []
    
    return t, final_ogs_info, seqs_in_ogs


# --- Functions for Monophyletic OGs ---

def add_monophyletic_ogs(t, taxonomy_db, ogs_info):
    """
    Finds and annotates all internal nodes marked as 'monophyletic_og' based on previous steps.
    """

    # Search is done once across the whole tree
    set_trees = set(t.search_nodes(monophyletic_og='True'))
    
    monophyletic_ogs = {} # og_name: node object
    
    for subtree in set_trees:
        
        if subtree.is_root:
            # Root OGs are handled in root_ogs for consistency.
            continue
        
        lca_dup = subtree.props.get('lca_dup')
        if lca_dup:
            og_name = f"{lca_dup}|{subtree.name}"
            subtree.add_prop('og_name', og_name)
            monophyletic_ogs[og_name] = subtree
        # Note: In the optimized version, we pass the node object, not the name string.

    # Annotate nodes that are monophyletic OGs
    _annot_monophyletic_ogs(monophyletic_ogs, taxonomy_db, ogs_info) 

    

def _annot_monophyletic_ogs(monophyletic_ogs, taxonomy_db, ogs_info):
    """
    Creates entries in the ogs_info dictionary for Monophyletic OGs.
    """
    
    for og_name, subtree_node in monophyletic_ogs.items():

        # Access node properties directly from the object
        lca_subtree = str(subtree_node.props.get('lca_node'))
        sci_name_lca = utils.get_sci_name(taxonomy_db, lca_subtree)
        
        values = [
            og_name, 
            lca_subtree, 
            sci_name_lca, 
            subtree_node.name, 
            subtree_node.props.get('sp_in', set()),
            subtree_node.props.get('leaves_in', []),
            subtree_node.props.get('lca_dup', '-'),
            subtree_node.props.get('sp_out', []),
            subtree_node.props.get('inparalogs_rate', '-'),
            subtree_node.props.get('so_score_dup', 0.0)
        ]
                  
        add_entry_to_ogs_info(ogs_info, values)

 


# --- Functions for Paraphyletic OGs ---

def add_paraphyletic_ogs(t, taxonomy_db, ogs_info, og_indexes) :

    """
    Recovers sequences that 'escaped' between OGs. 
    We identify nodes that are duplication nodes (D) but whose children are not fully OGs.
    Optimized to use og_indexes instead of repeated subtree searches.
    """

    paraphyletic_ogs = {} 
    
    for n in t.traverse():
        if n.is_leaf or n.is_root:
            continue

        lca_target = n.props.get('lca_node')

        if n.props.get('evoltype_2')=='D':

            # Skip if both children are already Monophyletic OGs
            if (n.children[0].props.get('monophyletic_og') == 'True' and 
                n.children[1].props.get('monophyletic_og') == 'True'):
                continue

            for ch in n.children:
                if ch.props.get('monophyletic_og') == 'True':
                    continue
                
                all_mems = set(ch.props.get('leaves_in', [])) # Sequences in the current child
                mem2remove = set()               
                
                # Nodos Monofiléticos bajo 'ch'
                monophyletic_nodes = og_indexes['monophyletic_og'].get('True', set())
                
                # Nodos de Duplicación (D) cuyo LCA coincide
                dup_nodes =  og_indexes.get('evoltype_2', {}).get('D', set())

                nodes_to_check = set()
                # Find relevant nodes under 'ch'
                for mn in monophyletic_nodes:
                    if mn in ch: 
                        if mn.props.get('lca_dup') == lca_target: 
                            nodes_to_check.add(mn)
                            
                for dn in dup_nodes:
                    if dn in ch: 
                        if dn.props.get('lca_node') == lca_target: 
                            nodes_to_check.add(dn)

                # Collect members to remove
                for n_ in nodes_to_check:
                    mem2remove.update(n_.props.get('leaves_in', []))

                # Remove seqs that belong to descendant OGs/Dups
                diff = all_mems.difference(mem2remove)
                
                if  len(diff)>1:
                    og_name = f"{ch.props.get('lca_node')}|{ch.name}!"
                    ch.add_prop('paraphyletic_og', 'True')
                    ch.add_prop('og_name', og_name)
                    # Store node object and members set
                    paraphyletic_ogs[og_name] = [ch, list(diff)] 

    # Annotate all the paraphyletic OGs
    _annot_paraphyletic_ogs(t, paraphyletic_ogs, taxonomy_db, ogs_info) 



def _annot_paraphyletic_ogs(t, paraphyletic_ogs, taxonomy_db, ogs_info):

    for pog_name, info in paraphyletic_ogs.items():
        
        pog_node = info[0]      # PhyloTree node object
        mems2check = set(info[1]) # Set of leaf names
       
        lca_pog = pog_node.props.get('lca_node')
        sci_name_lca = utils.get_sci_name(taxonomy_db, lca_pog)

        sp_in_og = set()
        list_seqs_in = []
        sp_outliers = []
        
        leaves_out = pog_node.props.get('leaves_out', set())
    
        for lname in mems2check:
            leaf = t[lname]

            # Check that seq is not a taxonomic outlier or long branch
            if lname not in leaves_out:
                taxid = leaf.props.get('taxid')

                list_seqs_in.append(lname)
                sp_in_og.add(taxid)
                
                # Add to the leaves a prop with all the paraphyletic OGs
                old_pOG = leaf.props.get('pOG', '')
                new_pOG = f"{old_pOG}@{pog_name}" if old_pOG else pog_name
                leaf.add_prop('pOG', new_pOG)
            else:
                # If it is an outlier, record its species ID
                sp_outliers.append(leaf.props.get('taxid')) 
            

        values = [
            pog_name, 
            lca_pog, 
            sci_name_lca, 
            pog_node.name, 
            sp_in_og, 
            list_seqs_in, 
            '-', 
            sp_outliers, 
            pog_node.props.get('inparalogs_rate', '-'), 
            pog_node.props.get('so_score', 0.0)
        ]
                  
        add_entry_to_ogs_info(ogs_info, values)
                
    return 


# --- OGs in root ---
# ⚠️ Requires the og_indexes pre-calculated in the main function (Step 4)

def root_ogs(t, taxonomy_db, ogs_info, og_indexes):

    """
    Handles special OGs related to the root:
    1. Recovery of OGs for taxonomic levels between the tree's LCA and LUCA.
    2. Paraphyletic OGs in the root's children (if they aren't monophyletic OGs).
    3. Monophyletic OG at the root (if marked).
    """

    lca_tree = t.props.get('lca_node')

    # 1. Recovery OGs from LUCA to LCA of the tree
    fromluca2lcatree = utils.get_lineage(taxonomy_db, lca_tree)
    
    lineage_copy = list(fromluca2lcatree)
    if 1 in lineage_copy: lineage_copy.remove(1) # Remove taxid 1 == root
    lineage_copy.reverse()
   
    
    if lca_tree not in [131567, 'r_root']:
        list_seqs_in = t.props.get('leaves_in', [])
        sp_in_og = set(t.props.get('sp_in', []))

        if lca_tree in lineage_copy: lineage_copy.remove(lca_tree) # Avoid repeating lca's root
       
        for taxa in lineage_copy:
            
            sci_name_taxa = utils.get_sci_name(taxonomy_db, taxa)
            og_name = f"{taxa}|{t.name}*"
            
            values = [og_name, taxa, sci_name_taxa, t.name, sp_in_og, list_seqs_in, '-', t.props.get('sp_out', []), t.props.get('inparalogs_rate', '-'), t.props.get('so_score', 0.0)]
            add_entry_to_ogs_info(ogs_info, values)
        
    
    # 2. Monophyletic OG at the root
    if t.props.get('monophyletic_og') == 'True':
       
        lca_root = t.props.get("lca_node")
        sci_name_lca = utils.get_sci_name(taxonomy_db, lca_root)
       
        og_name = f"{lca_root}|{t.name}"
            
        values = [og_name, lca_root, sci_name_lca, t.name, set(t.props.get('sp_in', [])), t.props.get('leaves_in', []), '-', t.props.get('sp_out', []), t.props.get('inparalogs_rate', '-'), t.props.get('so_score', 0.0)]
        add_entry_to_ogs_info(ogs_info, values)
    

    # 3. Paraphyletic OGs in the root's children (if neither child is a Monophyletic OG)
    else:
       
        if not (t.children[0].props.get('monophyletic_og') == 'True' and t.children[1].props.get('monophyletic_og') == 'True'):

            # Use the optimized logic similar to add_paraphyletic_ogs
            _add_root_paraphyletic_ogs(t, taxonomy_db, ogs_info, og_indexes)

def _add_root_paraphyletic_ogs(t, taxonomy_db, ogs_info, og_indexes):
    """Auxiliary: Adds paraphyletic OGs in the root's children."""
    
    lca_target = t.props.get('lca_node')
    
    for ch in t.children:
        if ch.props.get('monophyletic_og') == 'True':
            continue
        
        lca_ch = ch.props.get('lca_node')
        sci_name_taxa = utils.get_sci_name(taxonomy_db, lca_ch)
        all_mems = set(ch.props.get('leaves_in', []))
        mem2remove = set()
        
        # --- OPTIMIZATION: Use Indexes ---
        monophyletic_nodes = og_indexes['monophyletic_og'].get('True', set())
        dup_nodes = og_indexes['evoltype_d'].get(lca_target, set())

        nodes_to_check = set()
        for mn in monophyletic_nodes:
            if ch in mn.get_ancestors(): 
                if mn.props.get('lca_dup') == lca_ch: 
                    nodes_to_check.add(mn)
                    
        for dn in dup_nodes:
            if dn in ch:
                if dn.props.get('lca_node') == lca_ch: 
                    nodes_to_check.add(dn)
                    
        for n_ in nodes_to_check:
            mem2remove.update(n_.props.get('leaves_in', []))
        
        list_seqs_in = all_mems.difference(mem2remove)
        
        if len(list_seqs_in) > 1:
           
            og_name = f"{lca_ch}|{ch.name}!"
            ch.add_prop('paraphyletic_og', 'True')
            ch.add_prop('og_name', og_name)
            
            sp_in_og = set()
            sp_outliers = list()
            
            for lname in list_seqs_in:
                leaf = t.get_leaf(lname)
                taxid = leaf.props.get('taxid')
                sp_in_og.add(taxid)
                old_pOG = leaf.props.get('pOG', '')
                new_pOG = f"{old_pOG}@{og_name}" if old_pOG else og_name
                leaf.add_prop('pOG', new_pOG)
            
            values = [og_name, lca_ch, sci_name_taxa, list(ch.name), list(sp_in_og), list(list_seqs_in), '-', sp_outliers, ch.props.get('inparalogs_rate', '-'), ch.props.get('so_score', 0.0)]
            add_entry_to_ogs_info(ogs_info, values)


# --- General Utilities ---

def add_entry_to_ogs_info(ogs_info, values):

    """
    Adds a new entry to the ogs_info dictionary with the defined structure.
    """
    og_name, lca, sci_name_lca, assoc_node_name, sp_in_og, list_seqs_in, lca_dup, sp_outliers, inparalogs_rate, species_overlap = values

    len_sp_out = len(sp_outliers)
    if len_sp_out == 0:
        sp_outliers = ['-'] # Use a marker in the list
    
    ogs_info[og_name] = {
        'TaxoLevel': lca,
        'SciName_TaxoLevel': sci_name_lca.replace(' ', '_'),
        'AssocNode': assoc_node_name,
        'NumSP': len(sp_in_og),
        'NumMems': len(list_seqs_in),
        'Mems': list_seqs_in, # It is already a list or a set converted to list
        'Lca_Dup': lca_dup,
        'Species_Outliers': sp_outliers,
        'Num_SP_Outliers': len_sp_out,
        'Inparalogs_Rate': inparalogs_rate,
        'SP_overlap_dup': species_overlap
    } 
                          

def hierarchy(ogs_info, taxonomy_db):
    """
    Determines the vertical (parent/child) relationship between OGs based on overlap and taxonomy.
    Optimized: Parses OG names once and uses set operations for speed.
    """

    parsed_ogs = {}
    
    # 1. Parse OG names and members once
    for og_name, info in ogs_info.items():
        try:
            parts = og_name.split('|')
            tax_level = parts[0]
            og_level_raw = parts[1].replace('*', '').replace('!', '')
            # Assuming format: NAME-X, where X is the level number
            og_level = int(og_level_raw.split('-')[-1]) if '-' in og_level_raw else 0 
            
            parsed_ogs[og_name] = {
                'tax_level': tax_level, 
                'og_level': og_level, 
                'mems': set(info['Mems'])
            }
        except Exception as e:
            # Handle potential KeyError or ValueError from split/int conversion
            print(f"Warning: Skipping OG {og_name} due to unexpected name format: {e}")
            continue
    
    # 2. Compare every OG pair
    for og_name, target in parsed_ogs.items():
        ogs_up = set()
        ogs_down = set()
        
        for o_name, other in parsed_ogs.items():
            if og_name == o_name:
                continue

            # Check for members overlap
            if target['mems'] & other['mems']:
                
                # Check 1: Hierarchy by OG Level (based on node name position)
                if target['og_level'] < other['og_level']:
                    # Target node is "older" (less derived) -> Other is DOWNSTREAM
                    ogs_down.add(o_name)
                elif target['og_level'] > other['og_level']:
                    # Target node is "younger" (more derived) -> Other is UPSTREAM
                    ogs_up.add(o_name)
                
                # Check 2: Hierarchy by Taxonomic Level (for same OG Level)
                else: 
                    # If og_levels are the same, compare taxids
                    tax_level = target['tax_level']
                    t_lev = other['tax_level']

                    if tax_level != t_lev:
                        
                        # Determine lineage relationship
                        o_lin = utils.get_lineage(taxonomy_db, t_lev)
                        
                        if tax_level not in o_lin:
                            # Target's tax level is NOT ancestor of Other's tax level -> Target is 'up'
                            ogs_up.add(o_name)
                        elif tax_level in o_lin:
                            # Target's tax level IS ancestor of Other's tax level -> Target is 'down'
                            ogs_down.add(o_name)

        # 3. Assign results back to the original ogs_info dict
        ogs_info[og_name]['OG_down'] = ogs_down if ogs_down else '-'
        ogs_info[og_name]['OG_up'] = ogs_up if ogs_up else '-'