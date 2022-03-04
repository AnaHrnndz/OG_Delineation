from ete4 import  PhyloTree, SeqGroup, Tree
from ete4 import NCBITaxa, GTDBTaxa
from collections import Counter, OrderedDict, defaultdict
import sys
import os
import numpy as np
import pandas as pd
import json
import string
import random
import time
import argparse
import warnings
import pickle
from django.utils.crypto import get_random_string

warnings.filterwarnings("ignore", category=RuntimeWarning) 


class OrderedCounter(Counter, OrderedDict):
     'Counter that remembers the order elements are first seen'
     def __repr__(self):
         return '%s(%r)' % (self.__class__.__name__,
                            OrderedDict(self))
     def __reduce__(self):
         return self.__class__, (OrderedDict(self),)


class TimerError(Exception):

    """A custom exception used to report errors in use of Timer class"""

class Timer:
    timers = dict()

    def __init__(
        self,
        name=None,
        text="Elapsed time: {:0.4f} seconds",
        logger=None,
    ):
        self._start_time = None
        self.name = name
        self.text = text
        self.logger = logger

        # Add new named timers to dictionary of timers
        if name:
            self.timers.setdefault(name, 0)

    def start(self):

        """Start a new timer"""

        if self._start_time is not None:

            raise TimerError(f"Timer is running. Use .stop() to stop it")


        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        if self.logger:
            self.logger(self.text.format(elapsed_time))
        if self.name:
            self.timers[self.name] += elapsed_time

        return elapsed_time

def start_timers():

    timers = defaultdict()

    total_time = Timer("Total")
    _t0 = Timer("Start")
    _t1 = Timer("Outliers")
    _t2 = Timer("Iter_dups")
    _t3 = Timer("OGs_file")
    _t4 = Timer("Fastas")
    _t5 = Timer("Post_tree")
    _t6 = Timer('species_out')
    _t7 = Timer("lineage losses")
    _t8 = Timer("species losses")
    _t9 = Timer("test")

    

def clean_string(string):
    return string.replace("'","").replace("[", "(").replace("]", ")").replace("=","").replace(".","").replace("-"," ").replace(":","")


def id_generator():
    id_ = get_random_string(length=12)
    return id_


def get_reftree(t,taxonomy_db):
    print('GETTING REFERENTE TREE FROM USER TREE')
    taxid_list = list()
    for leaf in t:
        
        taxid_list.append(int(leaf.name.split('.')[0]))

    #Extract the smallest tree that connects all your query taxid
    reftree = taxonomy_db.get_topology(taxid_list)
    
    
    return reftree


def get_taxonomy_counter(reftree):
    print('GETTING TAXONOMY COUNTER FROM USER TREE')
    level2sp_mem = defaultdict(set)
    for l in reftree:
        lin = l.props.get('lineage')
        for tax in lin:
            level2sp_mem[tax].add(l.name)


    level2sp_num = defaultdict()
    for taxid, mems in level2sp_mem.items():
        level2sp_num[taxid] = len(mems)

    
    return level2sp_num
        

def get_lca_node(sp_list, taxonomy_db):
    lineages = OrderedCounter()
    nspecies = 0
    for sp in sp_list:
        nspecies += 1
        lineages.update(taxonomy_db.get_lineage(sp))

    lca = [l for l, count in lineages.items() if count == nspecies]
    if not lca:
        lca = ['Unk']
    return lca[-1]


def load_node_scores(n, SPTOTAL):
    leaf_targets = []
    for l in n.props.get('_leaves_in'):
        leaf_targets.append(l.split('.')[0])
    nspcs = len(set(leaf_targets))   
    
    dups_per_sp = Counter(leaf_targets)
    inparalogs_rate = np.median(list(dups_per_sp.values()))
    n.add_prop('inparalogs_rate', str(inparalogs_rate))    
    
    nseqs = len(n.props.get('_leaves_in'))    
    n.add_prop('score1', float(nspcs/SPTOTAL)) 
    n.add_prop('score2', float(nspcs/nseqs) if nseqs else 0.0)


#FUNCTIONS TO CALCULATE NODE LOSSES

def get_dup_score(n): 

    sp1 = n.props.get('_sp_in_ch1')
    sp2 = n.props.get('_sp_in_ch2')

    if len(sp1) == 0 and len(sp2) == 0:
        return 0

    a = np.array([len(sp1), len(sp2)]) 
    minval = int(np.min(a[np.nonzero(a)]))

    dup_score = float(len(sp1 & sp2) / minval)
    
    return dup_score


def count_lineage_losses(expected_sp, found_sp, reftree):
    def is_leaf_2(_n):
        if not _n.children: 
            return True
        elif len(found_sp & set(_n.get_leaf_names())) == 0:
            return True
        else: 
            return False

    if len(expected_sp) == 1: 
        return 0, 0 

    root = reftree.get_common_ancestor(expected_sp)

    losses = 0
    lin_losses = set()

    for leaf in root.iter_leaves(is_leaf_fn=is_leaf_2):
        losses +=1
        
        # Encontrar los linajes perdidos ralentiza mucho el script
        #anc = get_lca_node(leaf.get_leaf_names())
        #lin_losses.add(anc)

    return int(losses), lin_losses


def losses(node, reftree):

    sp1 = node.props.get('_sp_in_ch1')
    sp2 = node.props.get('_sp_in_ch2')
    
    if float(node.props.get('dup_score')) > 0.0:   
        losses1, lin_losses1 = count_lineage_losses(expected_sp=sp1|sp2, found_sp=sp1, reftree=reftree)
        losses2, lin_losses2 = count_lineage_losses(expected_sp=sp1|sp2, found_sp=sp2, reftree=reftree)
        node.add_prop('num_lineage_losses', [losses1, losses2])
        node.add_prop('_taxid_lineage_losses', [lin_losses1, lin_losses2])
        


def count_species_losses(expected_sp, found_sp, reftree):
    if len(expected_sp) == 1: 
        return 0, 0

    root = reftree.get_common_ancestor(expected_sp)
   
    losses = len(root) - len(found_sp)
    per_loss = losses / len(root)
    
    return int(losses), float(per_loss)


def percentage_losses(node, reftree):
  
    sp1 = node.props.get('_sp_in_ch1')
    sp2 = node.props.get('_sp_in_ch2')
    
    #Calculate number of losses -> number of expected species - number of found species
    #Calculate percentage of losses -> losses / number of expected species
    
    if float(node.props.get('dup_score')) > 0.0:    
        losses1, per_loss1 = count_species_losses(expected_sp=sp1|sp2, found_sp=sp1, reftree=reftree)
        losses2, per_loss2 = count_species_losses(expected_sp=sp1|sp2, found_sp=sp2, reftree=reftree) 
        node.add_prop('species_losses', [losses1, losses2])
        node.add_prop('species_losses_percentage', [per_loss1, per_loss2])

    else:
        node.add_prop('species_losses', [0, 0])
        node.add_prop('species_losses_percentage', [0.0, 0.0])



#FUNCTIONS TO DETECT OUTLIERS
def outliers_detection(n, outliers_node, outliers_reftree, CONTENT, taxonomy_counter):
    
    #count_lin : Count for each taxid how many seqs (Bilateria is level 6 -> {6:{33213:x}})
    count_lin = defaultdict(int)

    #count_lin_mem: for each index level, for each taxid at that index level, save leaves names
    count_lin_mem = defaultdict(lambda: defaultdict(list))

    #Save sp in node,  it doesnt matter that for the same sp there are several seqs
    sp_in_node = set()
    #For each taxonomic level in the leaves's lineage, save sp under that level 
    sp_per_level = defaultdict(set)

    for l in CONTENT[n]:
        sp_in_node.add(l.props.get('taxid'))
        for index, tax in enumerate(l.props.get('lineage')):
            count_lin[tax] += 1
            sp_per_level[tax].add(str(l.props.get('taxid')))
                
    sp2remove = set()
        
    # Detect outliers for each taxid from lineages
    for tax, num in count_lin.items():
        #num = How many sp there are at that taxonomic level in the node 
        #sp_in_node = all sp in the node 
        #levels_eggnog = number of sp at that taxonomic level in reftree (ej sp number of bilateria in eggnog reftree )
        
        #Not in use:
            #global_linages = numero de sp para ese nivel taxonomico en todo el arbol
            #SPTOTAL = numero de sp en total del arbo
        
        #% lineage in the node
        per_Node = num / len(sp_in_node)
       
        #% linege in tree
        #per_Tree = num / global_linages[level][tax]
        #per_Tree_total = num /SPTOTAL
        
        #% lineage in reftree
        per_Egg = num / taxonomy_counter.get(str(tax), 2)
        #per_Egg_total = levels_eggnog[str(tax)] / len(reftree)
    
        # If % of sp in the node at that taxonomic level is below 1%  and % of sp in reftree is below 5%, remove that species, that means:
            # 1. The taxonomic level has to be rare in the node (a few porifera species inside a node with all bilateria)
            # 2. Also there has to be few representatives from reftree at that taxonomic level, those leaves will be remove     
            #   (in the node there are 2 porifera species but in reftree there are 2000 porifera species)
            # 3. Taxonomic level that are rare in reftree, will be preserved
        if per_Egg < outliers_reftree and per_Node < outliers_node:
            sp2remove.update(sp_per_level[tax])
                    
    return sp2remove

        
#Add propeties to the root node
def annotate_root(ogs_info, node, taxonomy_db, total_mems_in_tree, outliers_node, outliers_reftree, sp_loss_perc, so_arq, so_bact, so_euk):
    name = node.props.get('name')
    node.add_prop('node_is_og', 'True')
    node.add_prop('so_score_dup', node.props.get('so_score'))
    node.add_prop('dup_lineage', list(taxonomy_db.get_lineage(node.props.get('lca_node'))))
    node.add_prop('_mems_og', total_mems_in_tree)
    node.add_prop('_dup_node_name', node.props.get('name'))
    
    
    node.add_prop("outliers_node", outliers_node)
    node.add_prop("outliers_reftree", outliers_reftree)
    node.add_prop("sp_loss_perc", sp_loss_perc)
    node.add_prop("so_arq", so_arq)
    node.add_prop("so_bact", so_bact)
    node.add_prop("so_euk", so_euk)
    
   


#Add properties to the  intertal nodes
def annotate_dups_ch(ogs_info, total_mems_in_ogs,taxid_dups_og, node, ch_node, taxonomy_db):

    if ch_node == 'ch1':
        og_name_ch = node.props.get('_ch1_name')
        og_ch_mems = node.props.get('_leaves_ch1')
        sp_ch = node.props.get('_sp_in_ch1')
        target_node = node.search_nodes(name=og_name_ch)[0]
    elif ch_node == 'ch2':
        og_name_ch = node.props.get('_ch2_name')
        og_ch_mems = node.props.get('_leaves_ch2')
        sp_ch = node.props.get('_sp_in_ch2')
        target_node = node.search_nodes(name=og_name_ch)[0]

    
    if len(sp_ch) > 2 and len(og_ch_mems) > 2:
        target_node.add_prop('node_is_og', 'True')
        target_node.add_prop('lca_dup', node.props.get('lca_node'))
        target_node.add_prop('so_score_dup', node.props.get('so_score'))
        target_node.add_prop('dup_lineage', list(taxonomy_db.get_lineage(node.props.get('lca_node'))))
        target_node.add_prop('_mems_og', list((og_ch_mems)))
        target_node.add_prop('_dup_node_name', node.props.get('name'))

        taxid_dups_og.add(node.props.get('lca_node'))
        node.add_prop('node_create_og', 'True')

        sp_in_og = target_node.props.get('_sp_in')
        sp_in_anc_og = node.props.get('_sp_in')
        diff_sp = list(set(sp_in_anc_og).difference(set(sp_in_og)))

        if len(diff_sp) >=1:
            target_node.add_prop('_miss_sp', diff_sp)
            
        total_mems_in_ogs.update(set(og_ch_mems))

        #ogs_info["name_og"] = target_node.props.get("name")
        ogs_info[og_name_ch]["name_dup"] = node.props.get("name")
        ogs_info[og_name_ch]["lca_dup"] = node.props.get("lca_node")
        ogs_info[og_name_ch]["lca_dup_lineage"] = node.props.get("lineage")
        ogs_info[og_name_ch]["tax_scope_og"] = target_node.props.get("lca_node")
        ogs_info[og_name_ch]["ogs_mems"] = target_node.props.get("_mems_og")



def expand_ogs(t, taxid_dups_og, total_mems_in_ogs,taxonomy_db):
    count_ogs = defaultdict(dict)
   
    
    for taxid in taxid_dups_og:
        count_ogs[taxid]['ogs_names'] = set()
        count_ogs[taxid]['mems'] = set()
        count_ogs[taxid]['ogs_down_included'] = set()
        count_ogs[taxid]["sci_name"] = taxonomy_db.get_taxid_translator([int(taxid)])[int(taxid)]
        
        for node in t.traverse('preorder'):
            if node.props.get('node_is_og') and taxid == node.props.get("lca_dup"):
                
                count_ogs[taxid]['ogs_names'].add(node.name)
                count_ogs[taxid]['mems'].update(node.props.get('_mems_og'))
                
                if node.props.get('ogs_down'):
                    count_ogs[taxid]['ogs_down_included'].update(node.props.get('ogs_down'))
                
            elif node.props.get('node_create_og') and taxid != node.props.get("lca_node") and taxid in node.props.get("lineage"):
                
                candidates_nodes = node.search_nodes(node_is_og="True")
                save_candidates = list()
                for candidate in candidates_nodes:
                    if taxid in candidate.props.get('lineage') and candidate.props.get('name') not in count_ogs[taxid]['ogs_down_included']:
                        save_candidates.append(candidate.props.get('name'))
                
                if len(save_candidates) >0:
                    count_ogs[taxid]['ogs_names'].add(node.props.get('name'))
                    count_ogs[taxid]['mems'].update(node.props.get('_leaves_in'))
                    count_ogs[taxid]['ogs_down_included'].update(node.props.get('ogs_down'))
                    total_mems_in_ogs.update(node.props.get('_leaves_in'))
        
    return(count_ogs)

def expand_ogs_annotated_tree(t, taxid_dups_og, total_mems_in_ogs,taxonomy_db):
    count_ogs = defaultdict(dict)
    
    for taxid in taxid_dups_og:
        count_ogs[taxid]['ogs_names'] = set()
        count_ogs[taxid]['mems'] = set()
        count_ogs[taxid]['ogs_down_included'] = set()
        count_ogs[taxid]["sci_name"] = taxonomy_db.get_taxid_translator([int(taxid)])[int(taxid)]
        
        for node in t.traverse('preorder'):
            if node.props.get('node_is_og') and taxid == node.props.get("lca_dup"):
                
                count_ogs[taxid]['ogs_names'].add(node.name)
                count_ogs[taxid]['mems'].update(set(node.props.get('_mems_og').split('|')))
                
                if node.props.get('ogs_down'):
                    count_ogs[taxid]['ogs_down_included'].update(set(node.props.get('ogs_down').split('|')))
                
            elif node.props.get('node_create_og') and taxid != node.props.get("lca_node") and taxid in node.props.get("lineage"):
                
                candidates_nodes = node.search_nodes(node_is_og="True")
                save_candidates = list()
                for candidate in candidates_nodes:
                    if taxid in candidate.props.get('lineage') and candidate.props.get('name') not in count_ogs[taxid]['ogs_down_included']:
                        save_candidates.append(candidate.props.get('name'))
                
                if len(save_candidates) >0:
                    count_ogs[taxid]['ogs_names'].add(node.props.get('name'))
                    count_ogs[taxid]['mems'].update(set(node.props.get('_leaves_in').split('|')))
                    count_ogs[taxid]['ogs_down_included'].update(set(node.props.get('ogs_down').split('|')))
                    total_mems_in_ogs.update(set(node.props.get('_leaves_in').split('|')))
                
             
    return(count_ogs)






def add_ogs_up_down(t, ogs_info):
    for node in t.traverse('preorder'):
        
        ogs_down = set()
        dups_down = list()
        for n in node.search_nodes(node_is_og="True"):
            if node.name != n.name:
                ogs_down.add(n.name)
                dups_down.append(n.up.name)
        if len(ogs_down) > 0:
            if node.props.get('node_is_og') and not node.props.get('is_root'):
                ogs_info[node.name]['ogs_down'] = ','.join(list(ogs_down))
                ogs_info[node.name]['dups_down'] = dups_down
            node.add_prop('ogs_down',ogs_down)
            node.add_prop('dups_down', dups_down)
        else:
            if node.props.get('node_is_og') and not node.props.get('is_root'):
                ogs_info[node.name]['ogs_down'] = '-'
                ogs_info[node.name]['dups_down'] = '-'
        ogs_up = set()
        dups_up = list()
        while node.up:
            if node.up.props.get('node_is_og'):
                if not node.up.props.get('is_root'):
                    ogs_up.add(node.up.props.get('name'))
            node = node.up
            dups_up.append(node.up)
        if len(ogs_up) > 0:
            if node.props.get('node_is_og') and not node.props.get('is_root'):
                ogs_info[node.name]['ogs_up'] = ','.join(list(ogs_up))
                ogs_info[node.name]['dups_up'] = dups_up
            node.add_prop('ogs_up', ogs_up)
            node.add_prop('dups_up', dups_up)
        else:
            if node.props.get('node_is_og') and not node.props.get('is_root'):
                ogs_info[node.name]['ogs_up'] = '-'
                ogs_info[node.name]['dups_up'] = '-'
    return t, ogs_info 



def run_load_tree(tree=None):
    print('-Load tree: ')
    if os.path.basename(tree).split('.')[-1] == 'pickle': 
        with open(args.tree, "rb") as handle:
            print('\t'+'FORMAT PICKLE')
            t = pickle.load(handle)

    elif os.path.basename(tree).split('.')[-1] == 'nw': 
        print('\t'+'FORMAT NEWICK')
        t = PhyloTree(args.tree)

    else:
        print('\t'+'FORMAT string')
        t = PhyloTree(tree, format = 1)

    return t


def run_load_taxonomy(taxonomy=None, user_taxonomy=None):
    print('-Load taxonomy:  ')
   
    if taxonomy == 'NCBI':
        if user_taxonomy != None:
            print('\t'+os.path.dirname(os.path.abspath(user_taxonomy)))
            taxonomy_db = NCBITaxa(user_taxonomy, memory = True)
        else:
            taxonomy_db = NCBITaxa(memory = True)
    
    elif taxonomy == 'GTDB':
        if user_taxonomy != None:
            print('\t'+os.path.dirname(os.path.abspath(user_taxonomy)))
            taxonomy_db = GTDBTaxa(user_taxonomy, memory = True)
        else:
            taxonomy_db = GTDBTaxa()

    
    return taxonomy_db


def run_load_reftree(rtree=None, t=None, taxonomy_db=None):

    print('-Load reftree:   ')
    if rtree:
            reftree = PhyloTree(rtree)
    else:
        reftree = get_reftree(t, taxonomy_db)
    print('\t'+'Len reftree: ', len(reftree))
    
    taxonomy_db.annotate_tree(reftree,  taxid_attr="name")
    

    return reftree

def run_load_taxonomy_counter(reftree=None, user_taxonomy_counter=None):
    print('-Load taxonomy counter:  ')
    if user_taxonomy_counter:
        if isinstance(user_taxonomy_counter, dict):
            taxonomy_counter = user_taxonomy_counter
        else:        
            with open(user_taxonomy_counter) as levels:
                taxonomy_counter = json.load(levels)
    else:
        taxonomy_counter = get_taxonomy_counter(reftree)
    return taxonomy_counter



def run_preanalysis(t, name_tree, taxonomy_db, midpoint):
    
    
    #_t0.start()
    
    print('\n'+'START PROCESSING TREE: '+name_tree)
    
    t.resolve_polytomy()

    if  midpoint == "Yes":
        print('-Run midpoint rooting')
        midpoint = t.get_midpoint_outgroup()
        t.set_outgroup(midpoint)


    # Parsing function used to extract species name from a node’s name.
    t.set_species_naming_function(lambda x: x.split('.')[0])
    
    # Adding additional information to any internal a leaf node (sci_name, taxid, named_lineage, lineage, rank)
    taxonomy_db.annotate_tree(t,  taxid_attr="species")

    # Total species in tree and total members in tree
    sp_set = set()
    total_mems_in_tree = set()

    for n in t.traverse("preorder"):
        if not n.is_leaf():
            #Create an ID for each internal node
            name = id_generator()
            n.add_prop('name', name)

        else:
            total_mems_in_tree.add(n.name)
            sp_set.add(n.props.get('taxid'))


   
    SPTOTAL = len(sp_set)
    print('-Len tree: ', len(t))
    print('-Total species in tree:', SPTOTAL)


    # Let's cache the list of leaves under each internal node. 
    # This is a global variable.

    CONTENT = t.get_cached_content()
    #_t0.stop()

    return t, sp_set, total_mems_in_tree, SPTOTAL, CONTENT    


def run_preanalysis_annot_tree(t, name_tree):
    sp_set = set()
    total_mems_in_tree = set()
    for n in t.traverse("preorder"):
        if n.is_leaf():
            total_mems_in_tree.add(n.name)
            sp_set.add(n.props.get('taxid'))

    SPTOTAL = len(sp_set)
    CONTENT = t.get_cached_content()
    return t, sp_set, total_mems_in_tree, SPTOTAL, CONTENT    


def get_analysis_parameters(t):
    root_node = t.get_tree_root()

    parameters = {
        'outliers_node': root_node.props.get("outliers_node"),
        'outliers_reftree': root_node.props.get("outliers_reftree"),
        'sp_loss_perc': root_node.props.get("sp_loss_perc"), 
        'so_arq': root_node.props.get("so_arq"), 
        'so_bact': root_node.props.get("so_bact"), 
        'so_euk': root_node.props.get("so_euk")
    }

    return parameters


def get_og_info(t):
    ogs_info = defaultdict(dict)
    taxid_dups_og = set()
    total_mems_in_ogs = set()

    for node in t.traverse():
        if node.props.get('node_is_og') and not node.is_root():
            
            
            ogs_info[node.name]["name_dup"] = node.props.get("_dup_node_name")
            ogs_info[node.name]["lca_dup"] = node.props.get("lca_dup")
            ogs_info[node.name]["lca_dup_lineage"] = node.props.get("dup_lineage")
            ogs_info[node.name]["tax_scope_og"] = node.props.get("lca_node")
            ogs_info[node.name]["ogs_mems"] = set(node.props.get('_mems_og').split('|'))
            ogs_info[node.name]["ogs_down"] = node.props.get('ogs_down')
            ogs_info[node.name]["dups_down"] = node.props.get('dups_down')
            ogs_info[node.name]["ogs_up"] = node.props.get('ogs_up')
            ogs_info[node.name]["dups_up"] = node.props.get('dups_up')
            
            total_mems_in_ogs.update(set(node.props.get('_mems_og').split('|')))

            if node.props.get("lca_dup"):
                taxid_dups_og.add(node.props.get("lca_dup"))

    return ogs_info, taxid_dups_og, total_mems_in_ogs



def run_outliers_dup_score(t, outliers_node, outliers_reftree, sp_loss_perc, so_arq, so_bact, so_euk, CONTENT, taxonomy_counter, taxonomy_db, SPTOTAL, reftree):

    
    #_t1.start()

    print('\n'+'DETECT OUTLIERS AND DUPS SCORE FUNCTION')
    print('\t'+'--Outliers threshodls:')
    print('\t'+'\t'+'--Node: '+ str(outliers_node))
    print('\t'+'\t'+'--Reftree: '+ str(outliers_reftree))
    print('\t'+'--Species losses percentage threshold: ' + str(sp_loss_perc))

    CONTENT = t.get_cached_content()
    #root_name = t.name
    t.add_prop('is_root', str('True'))


    for n in t.traverse("preorder"):

        
        n.del_prop('named_lineage')
        #n.del_prop('lineage')
        n.del_prop('_speciesFunction')
        n.del_prop('rank')

        sci_name = clean_string(n.props.get('sci_name'))
        n.props["sci_name"] = sci_name
        # n.add_prop('_sci_name', sci_name)
        # n.del_prop('sci_name')


        if n.is_leaf():
            n.add_prop('lca_node', n.props.get('taxid'))

        else:
            
            n.del_prop('common_name')
            #DETECT OUTLIERS
            sp_out = set()
            leaves_out = set()
            sp_in = set()
            leaves_in = set()

            #Add outliers from upper nodes
            sp_out_up = set()
            if n.up and n.up.props.get('sp_out'):
                sp_out_up = n.up.props.get('sp_out')
                for sp in sp_out_up:
                    if int(sp) in [leaf.props.get('taxid') for leaf in CONTENT[n]]:
                        sp_out.add(sp)


            
            
            #_t3.start()
            #Detect outliers in each of children nodes
            ch1 = n.children[0]
            sp_out.update(outliers_detection(ch1, outliers_node, outliers_reftree, CONTENT, taxonomy_counter))
            ch2 = n.children[1]
            sp_out.update(outliers_detection(ch2, outliers_node, outliers_reftree, CONTENT, taxonomy_counter))
           # _t3.stop()

            
            if len(sp_out) > 0:
                all_leafs = CONTENT[n]
                for l in all_leafs:
                    if str(l.props.get('taxid')) in sp_out:
                        leaves_out.add(l.props.get('name'))   
                    else:
                        sp_in.add(str(l.props.get('taxid')))
                        leaves_in.add(l.props.get('name'))
            else:
                for l in  CONTENT[n]:
                    sp_in.add(str(l.props.get('taxid')))
                    leaves_in.add(l.props.get('name'))


            all_spcs = set()

            #Save info for children_node_1
            ch1_name = n.children[0].props.get('name')
            sp_ch1 = set()
            leaves_ch1 = set()

            for l in ch1:
                if str(l.props.get('taxid')) not in sp_out:
                    all_spcs.add(str(l.props.get('taxid')))
                    sp_ch1.add(str(l.props.get('taxid')))
                    leaves_ch1.add(l.props.get('name'))

            #Save info for children_node_2
            ch2_name = n.children[1].props.get('name')
            sp_ch2 = set()
            leaves_ch2 = set()

            for l in ch2:
                if str(l.props.get('taxid')) not in sp_out:
                    all_spcs.add(str(l.props.get('taxid')))
                    sp_ch2.add(str(l.props.get('taxid')))
                    leaves_ch2.add(l.props.get('name')) 
            
            #Re-calculate species overlap after detect outliers
            overlaped_spces = set(sp_ch1 & sp_ch2)
            if len(overlaped_spces)>0:
                so_score = float(len(overlaped_spces) / len(all_spcs))
                n.add_prop('_overlap', list(overlaped_spces))
            else:
                so_score = 0.0
            
            #_t7.start()
            #Calculate common ancestor and rank for the species included in that node
            if len(sp_in) > 0:
                lca_node = get_lca_node(sp_in, taxonomy_db)
                #lca_node = taxonomy_db.get_topology(list(sp_in)).props.get('taxid')
                rank = clean_string(taxonomy_db.get_rank([lca_node])[lca_node])
                lin_lca = taxonomy_db.get_lineage(lca_node)
                n.add_prop('lineage', lin_lca)
                n.add_prop('taxid', lca_node)
            #_t7.stop()

            #SAVE PROPERTIES
            n.add_prop('lca_node', lca_node)
            lca_node_name = taxonomy_db.get_taxid_translator([lca_node])[lca_node]
            n.add_prop('lca_node_name', lca_node_name.replace(":", ""))
            n.add_prop('_rank', rank)
            n.add_prop('_sp_in', sp_in)
            n.add_prop('_sp_in_ch1', sp_ch1)
            n.add_prop('_sp_in_ch2', sp_ch2)
            n.add_prop('_ch1_name', ch1_name)
            n.add_prop('_ch2_name', ch2_name)
            n.add_prop('_leaves_in', leaves_in)
            n.add_prop('_leaves_ch1',list(leaves_ch1))
            n.add_prop('_leaves_ch2',list(leaves_ch2))
            n.add_prop('sp_out', list(sp_out))
            n.add_prop('_leaves_out', list(leaves_out))
            n.add_prop('so_score', so_score)

            
            # Load_node_scores add properties: score1, score2 and inpalalogs_rate
            #_t2.start()
            load_node_scores(n, SPTOTAL)
            dup_score = get_dup_score(n)        
            n.add_prop('dup_score', dup_score)
            #_t2.stop()


          
            #_t4.start()
            percentage_losses(n, reftree)
            #_t4.stop()


          
            #_t5.start()
            losses(n, reftree)
            #_t5.stop()


            if 2 in lin_lca:
                so_2_use = so_bact
            elif 2759 in lin_lca:
                so_2_use = so_euk
            elif 2157 in lin_lca:
                so_2_use = so_arq
            else:
                so_2_use = 0.2

            if so_score > so_2_use:
                if float(n.props.get('species_losses_percentage')[0]) > sp_loss_perc and float(n.props.get('species_losses_percentage')[1]) > sp_loss_perc:
                    n.add_prop('evoltype_2', 'FD')
                else:                
                    n.add_prop('evoltype_2', 'D')

            else:
                n.add_prop('evoltype_2', 'S')

    #_t1.stop()

    return t

def run_dups_and_ogs(t, outliers_node, outliers_reftree, sp_loss_perc, so_arq, so_bact, so_euk, taxonomy_db, total_mems_in_tree):

   
    #_t6.start()
    #ogs_info, save info of nodes that are OGs
    ogs_info = defaultdict(dict)

    #taxid_dups_og, set with the lca of the nodes that are OGs
    taxid_dups_og = set()
    total_mems_in_ogs = set()

    #Traverse tree to find the nodes that are "good" duplications and generate OGs.
    print('\n'+'ITER DUPS')
    print('\t'+'--species overlap used:')
    print('\t'+'\t'+'--Euk: '+ str(so_euk))
    print('\t'+'\t'+'--Bact :' + str(so_bact))
    print('\t'+'\t'+'--Arq :' + str(so_arq))


    for node in t.traverse("preorder"):
        node.del_prop('_speciesFunction')
        lin2check = node.props.get('lineage')
        
        if isinstance(lin2check, str):
            lin2check = lin2check.split('|')

        if 2 in lin2check:
            so_2_use = so_bact
        elif 2759 in lin2check:
            so_2_use = so_euk
        elif 2157 in lin2check:
            so_2_use = so_arq
        else:
            so_2_use = 0.2

        node.add_prop('min_sp_overlap', so_2_use)

        if node.is_root():
            annotate_root(ogs_info, node, taxonomy_db, total_mems_in_tree, outliers_node, outliers_reftree, sp_loss_perc, so_arq, so_bact, so_euk)

        if not node.is_leaf() and node.props.get('evoltype_2') == 'D' and so_2_use < float(node.props.get('so_score')) \
        and len(node.props.get('_leaves_in')) >2 and len(node.props.get('_sp_in')) > 2 :

            dups_under_node = []
            for n in node.search_nodes(evoltype_2='D'):
                if n.name!= node.name:
                    dups_under_node.append(n)

            #There are more dups under the node
            if len(dups_under_node) > 0:

                lca_target = node.props.get('lca_node')

                #Save Dups under child 1 and child2 that have the same lca_node
                dups_under_ch1 = node.children[0].search_nodes(evoltype_2='D', lca_node=lca_target)
                dups_under_ch2 = node.children[1].search_nodes(evoltype_2='D', lca_node=lca_target)

                save_dups_ch1 = defaultdict()
                save_dups_ch2 = defaultdict()

                # Check that dups under child1 and child2 (that have the same lca) fit out requirements : species overlap min requirement,
                # more than 2 leaves and more than 2 species
                if len(dups_under_ch1)> 0:

                    for n_ in  dups_under_ch1:
                        if float(n_.props.get('so_score'))> so_2_use  \
                        and len(n_.props.get('_leaves_ch1')) >2 and len(n_.props.get('_sp_in_ch1'))> 2 :
                            root2node = node.get_distance(node, n_, topology_only=True)
                            save_dups_ch1[n_.name] = root2node

                if len(dups_under_ch2)> 0:
                    for n_ in  dups_under_ch2:
                        if float(n_.props.get('so_score'))> so_2_use \
                        and len(n_.props.get('_leaves_ch2'))>2 and len(n_.props.get('_sp_in_ch2'))> 2 :
                            root2node = node.get_distance(node, n_, topology_only=True)
                            save_dups_ch2[n_.name] = root2node


                #If dups under child1 do not achieve our requirement, then child1 is OG
                if len(save_dups_ch1) == 0:  
                    annotate_dups_ch(ogs_info, total_mems_in_ogs, taxid_dups_og,node, 'ch1', taxonomy_db)
                    


                #If dups under child2 do not achieve our requirement, then child2 is OG
                if len(save_dups_ch2) == 0 :  
                    annotate_dups_ch(ogs_info, total_mems_in_ogs, taxid_dups_og,node, 'ch2', taxonomy_db)
                   



            # No more dups under node with the same lca_node
            elif len(dups_under_node) == 0:
                annotate_dups_ch(ogs_info, total_mems_in_ogs, taxid_dups_og,node, 'ch1', taxonomy_db)
                annotate_dups_ch(ogs_info, total_mems_in_ogs, taxid_dups_og,node, 'ch2', taxonomy_db)
                

    #_t6.stop()

    return t, total_mems_in_ogs, ogs_info, taxid_dups_og


def run_clean_properties(t):
    # clean strings
    all_props = set()
    for n in t.traverse():
        for string in ("sci_name", "lca_node_name", "common_name"):
            prop = n.props.get(string)
            if prop:
                n.props[string] = clean_string(prop)
        all_props.update(set(n.props.keys()))

    return t , all_props


def run_write_post_tree(t, name_tree, path_out, all_props):
    print('-Post tree')

    name_fam = name_tree.split('.')[0]
    if os.path.basename(args.tree).split('.')[-1] == 'pickle': 
        with open(path_out+'/post_'+name_fam+'.pickle', "wb") as handle:
            pickle.dump(t, handle)

    post_tree = path_out+'/'+'post_'+name_fam+'.nw'
    t.write(format=1, outfile=post_tree, properties=all_props, format_root_node = True)

    return

def run_write_fastas(fasta, name_tree, path_out, total_mems_in_ogs, total_mems_in_tree):
    print('-Fastas')

    fasta = SeqGroup(args.fasta)

    diff = (total_mems_in_tree.difference(total_mems_in_ogs))

    name_fam = name_tree.split('.')[0]
    not_og_fasta = path_out+'/'+name_fam+'.notOG.faa'

    with open(not_og_fasta, 'w') as  f_out:
        for name_seq in diff:
            aa = fasta.get_seq(name_seq)
            f_out.write('>'+name_seq+'\n'+aa+'\n')

    for name_og, info in ogs_info.items():
        if not 'is_root' in ogs_info[name_og].keys():
            lca = str(info['lca_dup'])
            with open(path_out+'/'+name_og+'_'+lca+'.faa', 'w') as f_out:

                for m in info["mems"]:
                    aa = fasta.get_seq(m)
                    f_out.write('>'+m+'\n'+aa+'\n')

    return

    
def get_newick(t, all_props):
    print(len(t))
    t = t.write(format=1, properties = all_props, format_root_node = True)
    return t


def run_app(tree, name_tree, outliers_node, outliers_reftree, sp_loss_perc, so_arq, so_bact, so_euk, reftree, user_counter, user_taxo, taxonomy_type, midpoint, path_out):


    #Load files and DBs:
    t = run_load_tree(tree = tree)
    taxonomy_db = run_load_taxonomy(taxonomy = taxonomy_type, user_taxonomy= user_taxo)
    reftree = run_load_reftree(rtree= reftree, t = t, taxonomy_db = taxonomy_db)
    taxonomy_counter = run_load_taxonomy_counter(reftree=reftree, user_taxonomy_counter = user_counter)


    #Pre-analysis: rooting, annotate tree, etc
    t , sp_set, total_mems_in_tree, SPTOTAL, CONTENT = run_preanalysis(t, name_tree, taxonomy_db, midpoint)

    #Outliers and Dups score functions
    t =  run_outliers_dup_score(t, outliers_node, outliers_reftree, sp_loss_perc, so_arq, so_bact, so_euk, CONTENT, taxonomy_counter, taxonomy_db, SPTOTAL, reftree)

    #Detect duplications and OGs
    t, total_mems_in_ogs, ogs_info = run_dups_and_ogs(t, outliers_node, outliers_reftree, sp_loss_perc, so_arq, so_bact, so_euk, taxonomy_db, total_mems_in_tree)
                                           
    #Clean properties
    t, all_props = run_clean_properties(t)

    #Write post_tree
    run_write_post_tree(t, name_tree, path_out, all_props)

    #Write OGs raw fasta files
    if args.fasta:
        fasta = args.fasta
        run_write_fastas(fasta, name_tree, path_out, total_mems_in_ogs, total_mems_in_tree)


    return(t)




if __name__ == "__main__":


    print('\n'+'READ TREE, SEQ FILE, REF TREE')


    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', dest = 'tree', required = True)
    parser.add_argument('--raw_fasta', dest = 'fasta')
    parser.add_argument('--output_path', dest = 'out_path', required= True)
    parser.add_argument('--species_overlap_euk', dest = 'so_euk' , default = 0.2, type = float)
    parser.add_argument('--species_overlap_bact', dest = 'so_bact' , default = 0.2, type = float)
    parser.add_argument('--species_overlap_arq', dest = 'so_arq' , default = 0.2, type = float)
    parser.add_argument('--outliers_in_node', dest = 'outliers_node' , default = 0.01, type = float)
    parser.add_argument('--outliers_in_reftree', dest = 'outliers_reftree' , default = 0.05, type = float)
    parser.add_argument('--species_losses_perct', dest = 'sp_loss_perc' , default = 0.9, type = float)
    parser.add_argument('--run_midpoint', dest = 'midpoint', choices = ['Yes', 'No'], default = 'no')
    parser.add_argument('--taxonomy_type', dest = 'taxonomy_type', choices = ['NCBI', 'GTDB'], default = 'NCBI')
    parser.add_argument('--user_taxonomy', dest = 'user_taxonomy')
    parser.add_argument('--user_taxonomy_counter', dest = 'user_taxonomy_counter')
    parser.add_argument('--reftree', dest = 'reftree')

    
    args = parser.parse_args()

    init_tree = args.tree
    name_tree = os.path.basename(init_tree)
    
    so_euk = args.so_euk
    so_bact = args.so_bact
    so_arq = args.so_arq
    outliers_node = args.outliers_node
    outliers_reftree = args.outliers_reftree
    sp_loss_perc = args.sp_loss_perc
    taxonomy_type = args.taxonomy_type
    midpoint = args.midpoint
    path_out = args.out_path

    if args.out_path:
        path_out = args.out_path
    
    if args.user_taxonomy:
        user_taxo = args.user_taxonomy
    else:
        user_taxo = None
    
    if args.reftree:
        rtree = args.reftree
    else:
        rtree = None

    if args.user_taxonomy:
        user_counter = args.user_taxonomy_counter
    else:
        user_counter = None

    
    _t = Timer('Total_time')
    _t0 = Timer('Preanalysis')
    _t1 = Timer('Outliers')
    _t2 = Timer('dup score')
    _t3 = Timer('outlier detection')
    _t4 = Timer('Perct losses')
    _t5 = Timer('Losses')
    _t6 = Timer('Iter dups')
    _t7 = Timer('test')

    _t.start()

    run_app(init_tree, name_tree, outliers_node, outliers_reftree, sp_loss_perc, so_arq, so_bact, so_euk, rtree, user_counter, user_taxo, taxonomy_type, midpoint, path_out)

    _t.stop()
    print(Timer.timers)  
