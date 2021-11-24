from ete4 import  PhyloTree, SeqGroup
from ete4 import NCBITaxa
from ete4.treeview import random_color
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
warnings.filterwarnings("ignore", category=RuntimeWarning) 

start_time = time.time()

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



def clean_string(string):
    return string.replace("'","").replace("[", "(").replace("]", ")").replace("=","").replace(".","").replace("-"," ").replace(":","")



#I have to change this function
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))



# def get_lca(n):
    # lineages = OrderedCounter()
    # nleaves = 0
    # #for n in n.get_leaves():        
    # for n in n.props.get('leaves_in'):
        # nleaves += 1
        # lineages.update(n.lineage)

    # lca = [l for l, count in lineages.items() if count == nleaves]
    # if not lca: 
        # lca = ['Unk']
    # return lca[-1]



def load_node_scores(n):
    leaf_targets = []
    for l in n.props.get('_leaves_in'):
        leaf_targets.append(l.split('.')[0])
    nspcs = len(set(leaf_targets))   
    
    dups_per_sp = Counter(leaf_targets)
    n.add_prop('inparalogs_rate', np.median(list(dups_per_sp.values())))    
    
    nseqs = len(n.props.get('_leaves_in'))    
    n.add_prop('score1', (nspcs/SPTOTAL)) 
    n.add_prop('score2', nspcs/nseqs if nseqs else 0.0)



#FUNCTIONS TO CALCULATE NODE LOSSES

# def count_losses(expected_sp, found_sp):
    # def is_leaf_2(_n):       
        # # Returns true if node is an actual leaf, or a branch without target species
        # if not _n.children: 
            # return True
        # elif not (refn2sp[node] & found_sp):
            # return True
        # else: 
            # return False

    # if len(expected_sp) == 1: 
        # return 0
    # root = reftree.get_common_ancestor(expected_sp)
    
    # losses = 0
    # n2losses = {}
    # for post, node in reftree.iter_prepostorder():
        # if post: 
            # n2losses[node] = sum([n2losses[ch] for ch in node.children])
        # else: 
            # if not node.children:
                # if node.name not in found_sp:
                    # n2losses[node] = 1
                # else: 
                    # n2losses[node] = 0
    # miss_sp = set()
    # for node in root.traverse(is_leaf_fn=is_leaf_2): 
        # if not (refn2sp[node] & found_sp) and (refn2sp[node] & expected_sp):
            # losses += 1
    # return losses


def get_dup_score(n): 

    sp1 = n.props.get('_sp_in_ch1')
    sp2 = n.props.get('_sp_in_ch2')

    if len(sp1) == 0 and len(sp2) == 0:
        return 0, None, None

    a = np.array([len(sp1), len(sp2)]) 
    minval = np.min(a[np.nonzero(a)])

    dup_score = len(sp1 & sp2) / minval
    print(n.name, dup_score)
    
    #Idea -> add_prop instead of return¿?
   
    return dup_score#, sp1, sp2



def count_lineage_losses(expected_sp, found_sp):
    def is_leaf_2(_n):
        if not _n.children: 
            return True
        elif len(found_sp & set(_n.get_leaf_names())) == 0:
            return True
        else: 
            return False

    if len(expected_sp) == 1: 
        return 0, None, 

    root = reftree.get_common_ancestor(expected_sp)
    root.annotate_ncbi_taxa(taxid_attr="name")
    
    losses = 0
    lin_losses = []
    for node in root.traverse(is_leaf_fn=is_leaf_2):
        if len(found_sp & set(node.get_leaf_names())) == 0:
            lin_losses.append(node.props.get('sci_name'))
            lin_losses.append(node.props.get('taxid'))
            losses += 1
    
    return losses, lin_losses

def losses(node):

    #dup_score, sp1, sp2 = get_dup_score(node)
    sp1 = n.props.get('_sp_in_ch1')
    sp2 = n.props.get('_sp_in_ch2')
    
    if float(node.props.get('dup_score')) > 0.0:   
        losses1, lin_losses1 = count_lineage_losses(expected_sp=sp1|sp2, found_sp=sp1)
        losses2, lin_losses2= count_lineage_losses(expected_sp=sp1|sp2, found_sp=sp2)
        # Add properties instead of return¿?
        #return dup_score, losses1, losses2
        node.add_prop('num_lineage_losses', [losses1, losses2])
        node.add_prop('_name_lineage_losses', [lin_losses1, lin_losses2])
        #return losses1, losses2
    # else:
        # return None, None, None


def count_species_losses(expected_sp, found_sp):
    def is_leaf_2(_n):
        if not _n.children: 
            return True
        elif len(found_sp & set(_n.get_leaf_names())) == 0:
            return True
        else: 
            return False

    if len(expected_sp) == 1: 
        return 0, None, 

    root = reftree.get_common_ancestor(expected_sp)
    root.annotate_ncbi_taxa(taxid_attr="name")
   
    losses = len(root) - len(found_sp)
    per_loss = losses / len(root)
    
    return losses, per_loss



def percentage_losses(node):
  
    #dup_score, sp1, sp2 = get_dup_score(node)
    sp1 = n.props.get('_sp_in_ch1')
    sp2 = n.props.get('_sp_in_ch2')
    
    
    #Calculate number of losses -> number of expected species - number of found species
    #Calculate percentage of losses -> losses / number of expected species
    if float(node.props.get('dup_score')) > 0.0:    
        losses1, per_loss1 = count_species_losses(expected_sp=sp1|sp2, found_sp=sp1)
        losses2, per_loss2 = count_species_losses(expected_sp=sp1|sp2, found_sp=sp2) 
        # Add properties instead of return¿?
        return dup_score, losses1, losses2, per_loss1, per_loss2
    else:
        return None, None, None, None, None


#FUNCTIONS TO DETECT OUTLIERS
def outliers_detection(n):
   
    #count_lin : Count for each index level, for each taxid at that index level, how many seqs (Bilateria is level 6 -> {6:{33213:x}})
    count_lin = defaultdict(lambda: defaultdict(int))

    #count_lin_mem: for each index level, for each taxid at that index level, save leaves names
    count_lin_mem = defaultdict(lambda: defaultdict(list))

    #Save sp in node,  it doesnt matter that for the same sp there are several seqs
    sp_in_node = set()
    #For each taxonomic level in the leaves's lineage, save sp under that level 
    sp_per_level = defaultdict(set)

    #leaves = CONTENT[n]
    for l in CONTENT[n]:

        #if l.props.get('taxid') not in sp_in_node: 
        sp_in_node.add(l.props.get('taxid'))
        for index, tax in enumerate(l.props.get('lineage')):
            count_lin[index][tax] += 1
            count_lin_mem[index][tax].append(l.props.get('name'))
            sp_per_level[tax].add(str(l.props.get('taxid')))
                
    sp2remove = set()

    #For each taxid from lineages
    for level, linages in count_lin.items():
            
        #Detect outliers only for index levels with more than one taxid, (level 0, all seqs belong to root, but level 1 seqs belong to euk and bact)
        if len(linages.keys())>1:
                
            for tax, num in linages.items():
                #num = How many sp there are at that taxonomic level in the node 
                #sp_in_node = all sp in the node 
                #levels_eggnog = number of sp at that taxonomic level in reftree (ej sp number of bilateria in eggnog reftree )
                
                #Not in use:
                    #global_linages = numero de sp para ese nivel taxonomico en todo el arbol
                    #SPTOTAL = numero de sp en total del arbol


                #% lineage in the node
                per_Node = num / len(sp_in_node)
                node_lin[n.name][tax] = num

                #% linege in tree
                #per_Tree = num / global_linages[level][tax]
                #per_Tree_total = num /SPTOTAL
                
                #% lineage in reftree
                per_Egg = num / levels_eggnog.get(str(tax), 2)
                #per_Egg_total = levels_eggnog[str(tax)] / len(reftree)
                

                # If % of sp in the node at that taxonomic level is below 1%  and % of sp in reftree is below 5%, remove that species, that means:
                    # 1. The taxonomic level has to be rare in the node (a few porifera species inside a node with all bilateria)
                    # 2. Also there has to be few representatives from reftree at that taxonomic level, those leaves will be remove     
                    #   (in the node there are 2 porifera species but in reftree there are 2000 porifera species)
                    # 3. Taxonomic level that are rare in reftree, will be preserved 

                if per_Egg <0.05 and per_Node < 0.01:
                    sp2remove.update(sp_per_level[tax])
                    
    return sp2remove

#Add propeties to the root node
def annotate_root(og_dict, node):
    name = node.props.get('name')
    og_dict[name]['mems'] = node.props.get('_all_mems')
    og_dict[name]['lca'] = node.props.get('lca_node')
    og_dict[name]['dup_lineage'] = list(ncbi.get_lineage(node.props.get('lca_node')))
    og_dict[name]['so_score'] = node.props.get('so_score')
    og_dict[name]['is_root'] = 'True'
    og_dict[name]['dup_node_name'] = name 

#Add properties to the  intertal nodes
def annotate_dups_ch(og_dict, total_mems_in_ogs, node, ch_node):

    if ch_node == 'ch1':
        og_name_ch = node.props.get('ch1_name')
        og_ch_mems = node.props.get('_leaves_ch1')
        sp_ch = node.props.get('_sp_in_ch1')
        target_node = node.search_nodes(name=og_name_ch)[0]
    elif ch_node == 'ch2':
        og_name_ch = node.props.get('ch2_name')
        og_ch_mems = node.props.get('_leaves_ch2')
        sp_ch = node.props.get('_sp_in_ch2')
        target_node = node.search_nodes(name=og_name_ch)[0]

    
    if len(sp_ch) > 2 and len(og_ch_mems) > 2:
        og_dict[og_name_ch]['mems'] = list((og_ch_mems))
        og_dict[og_name_ch]['tax_scope_og'] = target_node.props.get('lca_node')
        og_dict[og_name_ch]['lca_dup'] = node.props.get('lca_node')
        og_dict[og_name_ch]['dup_lineage'] = list(ncbi.get_lineage(node.props.get('lca_node')))
        og_dict[og_name_ch]['so_score_dup'] = node.props.get('so_score')
        og_dict[og_name_ch]['dup_node_name'] = node.props.get('name')
        
        target_node.add_prop('node_is_og', 'True')
        target_node.add_prop('lca_dup', node.props.get('lca_node'))
        target_node.add_prop('so_score_dup', node.props.get('so_score'))

        sp_in_og = target_node.props.get('_sp_in')
        sp_in_anc_og = node.props.get('_sp_in')
        diff_sp = list(set(sp_in_anc_og).difference(set(sp_in_og)))

        if len(diff_sp) >=1:
            
            miss_lin = defaultdict(list)
            for sp in diff_sp:
                lin = ncbi.get_lineage(sp)
                for l in lin:
                    if not l in ncbi.get_lineage(node.props.get('lca_node')):
                        l_name = ncbi.get_taxid_translator([l])[l]
                        miss_lin[l].append(sp)
                        break

            miss_lin_perc = defaultdict()
            for lin, sp_miss in miss_lin.items():
                total_node_lin = node_lin[node.props.get('name')][lin]
                perc_loss = len(sp_miss) / int(total_node_lin)
                l_name = ncbi.get_taxid_translator([lin])[lin]
                miss_lin_perc[l_name] = perc_loss
            
            miss_lineages = list(miss_lin_perc.keys())

            target_node.add_prop('miss_lineage', miss_lineages)
            target_node.add_prop('miss_sp', diff_sp)
           
        else:
            target_node.add_prop('miss_lineage', '-')
            
        total_mems_in_ogs.update(set(og_ch_mems))



print('READ TREE, SEQ FILE, REF TREE')


parser = argparse.ArgumentParser()
parser.add_argument('--tree', dest='tree', required=True)
parser.add_argument('--raw_fasta', dest='fasta', required=True)
parser.add_argument('--output_path', dest='out_path', required=True)
parser.add_argument('--midpoint', dest='midpoint', required=True, choices=['yes', 'no'])
parser.add_argument('--taxonomy', dest='taxonomy', default='/home/plaza/projects/eggnog6/pfamA_families/eggnog_experiments/data/levels2numSp.json')
parser.add_argument('--reftree', dest='reftree', default='/home/plaza/projects/eggnog6/pfamA_families/eggnog_experiments/data/totalNCBITree.nw')
parser.add_argument('--ncbitaxa', dest='ncbitaxa')
args = parser.parse_args()

t = PhyloTree(args.tree)
name_tree = os.path.basename(args.tree)
fasta = SeqGroup(args.fasta)
path_out = args.out_path
taxonomy = args.taxonomy
rtree = args.reftree

#timers
total_time = Timer("Total")
_t0 = Timer("Start")
_t1 = Timer("Outliers")
_t2 = Timer("Iter_dups")
_t3 = Timer("OGs_file")
_t4 = Timer("Fastas")
_t5 = Timer("Post_tree")

total_time.start()
_t0.start()

if args.ncbitaxa:
    print(str(args.ncbitaxa))
    ncbi = NCBITaxa(args.ncbitaxa)
else:
    ncbi = NCBITaxa()


print('-Load taxonomy')
levels_eggnog = {}
with open(taxonomy) as levels:
    levels_eggnog = json.load(levels)


print('-Load reftree')
reftree = PhyloTree(rtree)
#refn2sp = reftree.get_cached_content(store_attr='name', container_type=set)
print('--Len reftree: ', len(reftree))


print('START PROCESSING TREE')
t.resolve_polytomy()
if  args.midpoint  and args.midpoint == "yes":
    print('-Midpoint rooting')
    midpoint = t.get_midpoint_outgroup()
    t = t.set_outgroup(midpoint)

# Parsing function used to extract species name from a node’s name.
t.set_species_naming_function(lambda x: x.split('.')[0])
# Adding additional information to any internal a leaf node (sci_name, taxid, named_lineage, lineage, rank)
ncbi.annotate_tree(t,  taxid_attr="species")

# ref_leaves_names = []
# for l in reftree.get_leaves():
    # ref_leaves_names.append(l.name)


# Count how many species are in each taxonomic level
#counter_taxono = defaultdict(int)
#global_linages = defaultdict(lambda: defaultdict(int))
#taxa = {}

# Total species in tree and total members in tree
sp_set = set()
total_mems_in_tree = set()

# Check that all species has lineage in ncbi taxonomy
miss_lineages = []

# Detect all evolutionary events (Speciation or Duplicaton node) in the tree
t.get_descendant_evol_events()

for n in t.traverse("preorder"):
    if not n.is_leaf():
        #Create an ID for each internal node
        name = id_generator()
        n.add_prop('name', name)
    else:
        total_mems_in_tree.add(n.name)
        sp_set.add(n.props.get('taxid'))

        if len(n.props.get('lineage')) == 0:
            miss_lineages.append(n.props.get('taxid'))
        
        # if n.props.get('taxid') not in sp_set:
            # sp_set.add(n.props.get('taxid'))
            # for index, term in enumerate(n.props.get('lineage')):
                # global_linages[index][term] += 1
                # counter_taxono[term] +=1

print('-Missing Lineages: ', len(miss_lineages))
SPTOTAL = len(sp_set)
print('-Total species in tree:', SPTOTAL)

# Let's cache the list of leaves under each internal node. 
# This is a global variable.
CONTENT = t.get_cached_content()

# For each node, for each taxid in lineage, save number of species under that taxid,
# This dict is fill during outliers detection and its used in annotate_dups_ch() to annotate missing lineages
node_lin = defaultdict(dict)

_t0.stop()        

print('DETECT OUTLIERS AND DUPS SCORE FUNCTION')
_t1.start()
root_name = t.name
t.add_prop('is_root', 'True')
t.add_prop('_all_mems', clean_string(str(list(total_mems_in_tree))))

for n in t.traverse("preorder"):

    n.del_prop('named_lineage')
    sci_name = clean_string(n.props.get('sci_name'))
    n.add_prop('sci_name', sci_name)
    
    if n.is_leaf():
        n.add_prop('lca_node', n.props.get('taxid'))
        
    else:
        
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

        #Only detect outliers if node is a duplication node (For saving time)
        if n.props.get('evoltype') == 'D':
            sp_out.update(outliers_detection(n))
            
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
        else:
            for l in  CONTENT[n]:
                sp_in.add(str(l.props.get('taxid')))
                leaves_in.add(l.props.get('name'))
            
        # sp_out_sciName = []
        # if sp_out:
            # for sp in sp_out:
                # sciName = ncbi.get_taxid_translator([int(sp)])[int(sp)]
                # sp_out_sciName.append(sciName.replace("'","").replace("[", "(").replace("]", ")").replace("=","").replace(".","").replace("-"," ").replace(":"," "))
        
        # n.add_prop('sp_out_name', sp_out_sciName)

        #Save info for children_node_1
        all_spcs = set()
        ch1 = n.children[0]
        ch1_name = n.children[0].props.get('name')
        sp_ch1 = set()
        leaves_ch1 = set()

        for l in ch1:
            if str(l.props.get('taxid')) not in sp_out:
                all_spcs.add(str(l.props.get('taxid')))
                sp_ch1.add(str(l.props.get('taxid')))
                leaves_ch1.add(l.props.get('name'))

        #Save info for children_node_2
        ch2 = n.children[1]
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
        else:
            so_score = 0.0

        
        #Calculate common ancestor and rank for the species included in that node
        #######SLOW STEP#########
        if len(sp_in) > 0:
            lca_node = ncbi.get_topology(list(sp_in)).props.get('taxid')
            rank = ncbi.get_rank([lca_node])[lca_node]
            lin_lca = ncbi.get_lineage(lca_node)
        #########################
        
        
        #If there is upper node, save upper node ancestor
        # if n.up:
            # lca_node_up = n.up.props.get('lca_node')
            # lin_lca_node_up = ncbi.get_lineage(lca_node_up)
        # else:
            # lin_lca_node_up = []
        
        
        #If node is a speciation node, outliers are not detected, so ancestor node can be more ancient that previous ancestor node
        #If ancestor of current node is more ancient that ancestor of upper node, change ancestor of current node to ancestor of upper node
        # if len(lin_lca) < len(lin_lca_node_up):
            # n.add_prop('lca_node', lca_node_up)
            # lca_node_name = ncbi.get_taxid_translator([lca_node_up])[lca_node_up]
            
        # else:
        n.add_prop('lca_node', lca_node)
        lca_node_name = ncbi.get_taxid_translator([lca_node])[lca_node]
        

        #SAVE PROPERTIES
        n.add_prop('lca_node_name', lca_node_name.replace(":", ""))
        n.add_prop('rank', rank)
        n.add_prop('_sp_in', sp_in)
        n.add_prop('_sp_in_ch1', sp_ch1)
        n.add_prop('_sp_in_ch2', sp_ch2)
        n.add_prop('ch1_name', ch1_name)
        n.add_prop('ch2_name', ch2_name)
        n.add_prop('_leaves_in', leaves_in)
        n.add_prop('_leaves_ch1',list(leaves_ch1))
        n.add_prop('_leaves_ch2',list(leaves_ch2))
        n.add_prop('sp_out', list(sp_out))
        n.add_prop('_leaves_out', list(leaves_out))
        n.add_prop('so_score', so_score)
        n.add_prop('lineage', clean_string(str(n.props.get('lineage'))))

        if so_score > 0:
            n.add_prop('evoltype_2', 'D')
        else:
            n.add_prop('evoltype_2', 'S')
        

        #load_node_scores add properties: score1, score2 and inpalalogs_rate
        load_node_scores(n)
        dup_score = get_dup_score(n)
        
        n.add_prop('dup_score', dup_score)

       
        ######SLOW STEP#######
        dup_score, losses1, losses2, per_loss1, per_loss2 = percentage_losses(n)

        if per_loss1 != None and per_loss1 > 0.7 and per_loss2 != None  and per_loss2 >0.7:
            n.add_prop('evoltype_2', 'FD')

        n.add_prop('species_loss',  [dup_score, losses1, losses2, per_loss1, per_loss2])
        
        ######SLOW STEP#######
        #losses() add propeties 'num_lineage_losses' and '_name_lineage_losses' 
        losses(n)
        
        
        # if dup_score != None:
            # n.add_prop('lineage_loss', list(loss1, loss2))
        # else:
            # n.add_prop('lineage_loss', [])
        
        
_t1.stop()

_t2.start()
#og_dict, save info of nodes that are OGs
og_dict = defaultdict(dict)

#tax_dups_og, set with the ancestors of the nodes that are OGs
tax_dups_og = set()
total_mems_in_ogs = set()

#Traverse tree to find the nodes that are "good" duplications and generate OGs.
print('ITER DUPS')
for node in t.traverse("preorder"):

    if node.is_root():
        annotate_root(og_dict, node)
                        
    if not node.is_leaf() and node.props.get('evoltype_2') == 'D' and 0.2 < float(node.props.get('so_score')) \
    and len(node.props.get('_leaves_in')) >2 and len(node.props.get('_sp_in')) > 2 :
    #and float( len(properties[node.name]['sp_in']) / counter_taxono[(node.lca_node)] ) > 0.6:
        
        dups_under_node = []
        for n in node.search_nodes(evoltype_2='D'):
            if n.name!= node.name:
                dups_under_node.append(n)
        
        #Hay mas dups bajo ese nodo
        if len(dups_under_node) > 0:
                
            lca_target = node.props.get('lca_node')

            dups_under_ch1 = node.children[0].search_nodes(evoltype_2='D', lca_node=lca_target)
            dups_under_ch2 = node.children[1].search_nodes(evoltype_2='D', lca_node=lca_target)
            
            save_dups_ch1 = defaultdict()
            save_dups_ch2 = defaultdict()
            
            if len(dups_under_ch1)> 0:
                #dups por debajo del nodo principal anotadas al mismo nivel
                for n_ in  dups_under_ch1:
                    if float(n_.props.get('so_score'))> 0.2  \
                    and len(n_.props.get('_leaves_ch1')) >2 and len(n_.props.get('_sp_in_ch1'))> 2 :
                    #and float( len(properties[n_.name]['sp_in']) / counter_taxono[(lca_target)] ) > 0.6:

                        root2node = node.get_distance(node, n_, topology_only=True)
                        save_dups_ch1[n_.name] = root2node
            
            if len(dups_under_ch2)> 0:
                for n_ in  dups_under_ch2:
                    if float(n_.props.get('so_score'))> 0.2  \
                    and len(n_.props.get('_leaves_ch2'))>2 and len(n_.props.get('_sp_in_ch2'))> 2 :
                    #and float( len(properties[n_.name]['sp_in']) / counter_taxono[(lca_target)] ) > 0.6:

                        root2node = node.get_distance(node, n_, topology_only=True)
                        save_dups_ch2[n_.name] = root2node

            #si las dups que hay por debajo del hijo1 no cumplen los requisitos, entonces hijo1 es OG
            if len(save_dups_ch1) == 0:  
                annotate_dups_ch(og_dict, total_mems_in_ogs, node, 'ch1')
                
            #si las dups que hay por debajo del hijo2 no cumplen los requisitos, entonces hijo2 es OG    
            if len(save_dups_ch2) == 0 :  
                annotate_dups_ch(og_dict, total_mems_in_ogs, node, 'ch2')
                
        #No hay mas dups   
        elif len(dups_under_node) == 0:
            annotate_dups_ch(og_dict, total_mems_in_ogs, node, 'ch1')
            annotate_dups_ch(og_dict, total_mems_in_ogs, node, 'ch2')

_t2.stop()


###########################################################################################################3

_t3.start()
print('ogs tree')
for name_og, info in og_dict.items():
    mems_og = set(info['mems'])
    og_dict[name_og]['anc_og'] = []
    node_intersection = []
    max_dist = t.get_distance(t,name_og, topology_only=True)
    og_dict[name_og]['dist'] = max_dist

    #save all og with shared members
    for n, i in og_dict.items():
        og2compare = set(i['mems'])
        if mems_og.intersection(og2compare):
            node_intersection.append(n)

    #get distance from node with shared memers to root
    save_dups = defaultdict()
    for node in node_intersection:
        root2node = t.get_distance(t,node, topology_only=True)
        if root2node < max_dist:
            save_dups[node] = root2node
    sort_dups = {k: v for k, v in sorted(save_dups.items(), key=lambda item: item[1] ,reverse = True)}
    prev_og = defaultdict(dict)

    if len(sort_dups) > 0:
        
        for og_anc, dist in sort_dups.items():
            og_anc_node = t.search_nodes(name=og_anc)[0]
            lca = og_anc_node.props.get('lca_node')
            prev_og[og_anc]['dist'] = dist
            prev_og[og_anc]['lca'] = lca

        og_dict[name_og]['anc_og'] = prev_og
        
    else:
        root_node = t.get_tree_root()
        root_lca = root_node.props.get('lca_node')
        prev_og[root_node.name]['dist'] = 0.0
        prev_og[root_node.name]['lca'] = root_lca
        og_dict[root_node.name]['anc_og'] = prev_og
        
og_mems_out = path_out+'/'+name_tree+'_mems.json'
with open(og_mems_out, 'w') as f:
    json.dump(og_dict, f, indent=2)

_t3.stop()

        
######################################################################################

_t4.start()
print('#Writing fastas')
diff = (total_mems_in_tree.difference(total_mems_in_ogs))

name_fam = name_tree.split('.')[0]
not_og_fasta = path_out+'/'+name_fam+'.notOG.faa'

with open(not_og_fasta, 'w') as  f_out:
    for name_seq in diff:
        aa = fasta.get_seq(name_seq)
        f_out.write('>'+name_seq+'\n'+aa+'\n')

for name_og, info in og_dict.items():
    if not 'is_root' in og_dict[name_og].keys():
        lca = str(info['lca_dup'])
        with open(path_out+'/'+name_og+'_'+lca+'.faa', 'w') as f_out:

            for m in info["mems"]:
                aa = fasta.get_seq(m)
                f_out.write('>'+m+'\n'+aa+'\n')

_t4.stop()

        
############################################################################

_t5.start()
# clean strings
all_props = set()
for n in t.traverse():
    for string in ("sci_name", "lca_node_name"):
        prop = n.props.get(string)
        if prop:
            n.props[string] = clean_string(prop)
    all_props.update(set(n.props.keys()))


print('#Writing post tree')
post_tree = path_out+'/'+'post_'+name_fam+'.nw'
t.write(format=1, outfile=post_tree, properties=all_props)


_t5.stop()  
total_time.stop()
print(Timer.timers)     


        