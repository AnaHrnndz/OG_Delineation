
from ete4 import  PhyloTree, SeqGroup, Tree
from ete4 import NCBITaxa, GTDBTaxa
from collections import Counter, OrderedDict, defaultdict
import sys
import os
import numpy as np
#import pandas as pd
import json
import string
import random
import time
import argparse
import warnings
import pickle
import re
import glob
import subprocess
from django.utils.crypto import get_random_string
import csv
import tarfile
import os
import shutil
import pathlib


cwd =  str(pathlib.Path(__file__).parent.resolve())
print(cwd)
bin_path =  cwd+'/bin'
data_path = cwd+'/data'
sys.path.append(bin_path)

warnings.filterwarnings("ignore", category=RuntimeWarning)

chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'

OGD_props = list()
with open(data_path+'/OGD_props.txt') as f:
    for line in f:
        OGD_props.append(line.strip())


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

    """
        Remove problematic characters for newick format
    """
    #return string.replace("'","").replace("[", "(").replace("]", ")").replace("=","").replace(".","").replace("-"," ").replace(":","")
    clean_string = re.sub(r"'|\[|\]|\=|\.|\-|\:", "", string)
    return clean_string

def get_newick(t, all_props):

    """
        Return tree in newick format with annotations
    """

    t = t.write(props=all_props, format_root_node=True)
    return t

def check_nodes_up(node):

    """
        Find OGs in upper nodes
    """

    ogs_up = set()
    dups_up = list()
    while node.up:
        if node.up.props.get('node_is_og'):
            if not node.up.props.get('is_root'):
                ogs_up.add(node.up.props.get('name'))
                dups_up.append(node.up.up.props.get('name'))
        node = node.up

    return ogs_up, dups_up


####    FUNCTIONS FOR OG DELINEATION WEB    ####

def clean_folder(dir, aln_path):

    """
        Remove files in tmp dir created for web analysis
    """

    filelist = glob.glob(os.path.join(dir, "*"))
    for f in filelist:

        if f != aln_path:
            os.remove(f)


def run_preanalysis_annot_tree(t, name_tree):

    """
        Obtain basic info from the tree
        TODO: return less things
    """

    sp_set = set()
    total_mems_in_tree = set()
    taxid_dups = set()
    names_ogs = set()


    for n in t.traverse("preorder"):

        if n.is_leaf:

            total_mems_in_tree.add(n.name)
            sp_set.add(n.props.get('taxid'))
        else:
            if n.props.get('node_is_og'):
                names_ogs.add(n.name)
                taxid_dups.add(n.props.get('lca_dup'))

    SPTOTAL = len(sp_set)
    CONTENT = t.get_cached_content()
    t.dist = 0.01

    return t, sp_set, total_mems_in_tree, SPTOTAL, CONTENT, taxid_dups, names_ogs


def get_analysis_parameters(t):

    """
        Get the parameters of OGD analysis from the root
    """

    root_node = t.get_tree_root()

    parameters = {
        'outliers_node': root_node.props.get("outliers_node"),
        'outliers_reftree': root_node.props.get("outliers_reftree"),
        'sp_loss_perc': root_node.props.get("sp_loss_perc"),
        'so_cell_org': root_node.props.get("so_cell_org"),
        'so_arq': root_node.props.get("so_arq"),
        'so_bact': root_node.props.get("so_bact"),
        'so_euk': root_node.props.get("so_euk"),
        'inherit_outliers': root_node.props.get('inherit_out')
    }

    return parameters


def get_taxlevel2ogs(t, taxid_dups_og, total_mems_in_ogs, taxonomy_db):

    """
        Table needed in web app
        For each taxlevel that create ogs, keep:
            ogs names, all mems in that ogs, ogs below that taxlevel, sci_name
    """

    taxlev2ogs = defaultdict(dict)


    for taxid in taxid_dups_og:
        taxlev2ogs[taxid]['ogs_names'] = set()
        taxlev2ogs[taxid]['mems'] = set()
        taxlev2ogs[taxid]['ogs_down_included'] = set()
        taxlev2ogs[taxid]["sci_name"] = taxonomy_db.get_taxid_translator([int(taxid)])[int(taxid)]

        for node in t.traverse('preorder'):
            #Keep ogs at that taxlevel
            if node.props.get('node_is_og') and taxid == node.props.get("lca_dup"):

                taxlev2ogs[taxid]['ogs_names'].add(node.name)
                taxlev2ogs[taxid]['mems'].update(set(node.props.get('_mems_og').split('|')))

                if node.props.get('_ogs_down'):
                    taxlev2ogs[taxid]['ogs_down_included'].update(node.props.get('_ogs_down'))

            #Keep ogs below taxlevel
            elif node.props.get('node_create_og') and taxid != node.props.get("lca_node") and taxid in node.props.get("lineage"):

                candidates_nodes = node.search_nodes(node_is_og="True")
                save_candidates = list()
                for candidate in candidates_nodes:
                    if isinstance(candidate.props.get('lineage'), list) :
                        lin_cadidate = candidate.props.get('lineage')
                    else:
                        lin_cadidate = candidate.props.get('lineage').split('|')

                    if taxid in  lin_cadidate  and candidate.props.get('name') not in taxlev2ogs[taxid]['ogs_down_included']:
                        save_candidates.append(candidate.props.get('name'))

                if len(save_candidates) >0:
                    taxlev2ogs[taxid]['ogs_names'].add(node.props.get('name'))
                    taxlev2ogs[taxid]['mems'].update(node.props.get('_leaves_in'))
                    taxlev2ogs[taxid]['ogs_down_included'].update(node.props.get('_ogs_down'))



    return(taxlev2ogs)


def get_og_info(t):

    """
        Get info about OGs for to show in results web app
    """

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
            ogs_info[node.name]["ogs_down"] = node.props.get('_ogs_down')
            ogs_info[node.name]["dups_down"] = node.props.get('_dups_down')
            ogs_info[node.name]["ogs_up"] = node.props.get('_ogs_up')
            ogs_info[node.name]["dups_up"] = node.props.get('_dups_up')
            ogs_info[node.name]["total_leaves"] = node.props.get('total_leaves')
            ogs_info[node.name]["num_sp_OGs"] = node.props.get('len_sp_in')
            ogs_info[node.name]["num_sp_out"] = node.props.get('sp_out','0')
            if node.props.get('recover_seqs'):
                ogs_info[node.name]['recovery_mems'] = set(node.props.get('recover_seqs').split('|'))
                total_mems_in_ogs.update(set(node.props.get('recover_seqs').split('|')))
            else:
                ogs_info[node.name]['recovery_mems'] = set()

            if node.props.get('_ogs_down'):
                for n_name_down in node.props.get('_ogs_down').split('|'):
                    n_down  = t.search_nodes(name = n_name_down)[0]
                    ogs_info[node.name]["ogs_mems"].update(set(n_down.props.get('_mems_og').split('|')))


            total_mems_in_ogs.update(set(node.props.get('_mems_og').split('|')))

            if node.props.get("lca_dup"):
                taxid_dups_og.add(node.props.get("lca_dup"))


    return ogs_info, taxid_dups_og, total_mems_in_ogs


def prune_tree(t, total_mems_in_ogs):

    """
        Prune tree to keep only mems that belong to at least one OG
        We need this version of the tree to visualization
    """

    t.prune(list(total_mems_in_ogs))
    return t

def clean_tree(t):
    '''
    If tree was annotated with OGD, delete OGD props to run new OGD analysis
    '''

    for node in t.traverse():
        for prop in OGD_props:
                node.del_prop(prop)

    return t

#  def load_node_scores(n, SPTOTAL):

    # """
        # Calculate:
            # inparalogs_rate
            # score 1
            # score 2

        # TODO: split function
            #   change names of score1 and score2
    # """

    # leaf_targets = []
    # for l in n.props.get('_leaves_in'):
        # leaf_targets.append(l.split('.')[0])
    # nspcs = len(set(leaf_targets))

    # dups_per_sp = Counter(leaf_targets)

    # inparalogs_rate = np.median(list(dups_per_sp.values()))
    # n.add_prop('inparalogs_rate', str(inparalogs_rate))

    # nseqs = len(n.props.get('_leaves_in'))
    # n.add_prop('score1', float(nspcs/SPTOTAL))
    # n.add_prop('score2', float(nspcs/nseqs) if nseqs else 0.0)


# def taxlev2ogs_annotated_tree(t, taxid_dups_og, taxonomy_db):

    # """

        # For each taxonomic level get number of OGs and members included
        # We need that info for the table in the web app
        # TODO: Delete??? I do not use this function in web app
    # """

    # taxlev2ogs = defaultdict(dict)
    # total_mems_in_ogs = set()

    # for taxid in taxid_dups_og:
        # taxlev2ogs[taxid]['ogs_names'] = set()
        # taxlev2ogs[taxid]['mems'] = set()
        # taxlev2ogs[taxid]['ogs_down_included'] = set()
        # taxlev2ogs[taxid]["sci_name"] = taxonomy_db.get_taxid_translator([int(taxid)])[int(taxid)]

        # for node in t.traverse('preorder'):
            # if node.props.get('node_is_og') and taxid == node.props.get("lca_dup"):

                # taxlev2ogs[taxid]['ogs_names'].add(node.name)
                # taxlev2ogs[taxid]['mems'].update(set(node.props.get('_mems_og').split('|')))
                # if node.props.get('recover_seqs'):
                    # taxlev2ogs[taxid]['mems'].update(node.props.get('recover_seqs').split('|'))


                # if node.props.get('_ogs_down'):
                    # taxlev2ogs[taxid]['ogs_down_included'].update(set(node.props.get('_ogs_down').split('|')))


                    # for n_name_down in node.props.get('_ogs_down').split('|'):
                        # n_down  = t.search_nodes(name = n_name_down)[0]
                        # taxlev2ogs[taxid]['mems'].update(set(n_down.props.get('_mems_og').split('|')))
                        # if n_down.props.get('recover_seqs'):
                            # taxlev2ogs[taxid]['mems'].update(n_down.props.get('recover_seqs').split('|'))

        # total_mems_in_ogs.update(taxlev2ogs[taxid]['mems'])

    # return(taxlev2ogs, total_mems_in_ogs)





####   MAIN FUNCTIONS  ####

## 1. Load initial info   ##

def load_tree_local(tree=None):

    """
        Load Tree from file
    """

    print('-Load tree:', tree)

    # if os.path.basename(tree).split('.')[-1] == 'pickle':
        # with open(args.tree, "rb") as handle:
            # print('\t'+'FORMAT: PICKLE')
            # t = pickle.load(handle)

    # elif os.path.basename(tree).split('.')[-1] != 'pickle':
    t = PhyloTree(open(tree))

    # else:
        # print('\t'+'WRONG TREE FORMAT')

    return t


def load_tree_web(tree=None):

    """
        Load tree from newick format
    """

    print('-Load tree: ')

    t = PhyloTree((tree))

    return t


def load_taxonomy(taxonomy=None, user_taxonomy=None):

    """
        Load taxonomy from local server or from user file
        Local server:
            -Eggnog 5
            -Eggnog 6
            -GTDB????
    """

    if taxonomy == 'NCBI':
        if user_taxonomy != None:
            taxonomy_db = NCBITaxa(user_taxonomy, memory = True)
        else:
            taxonomy_db = NCBITaxa(memory = True)

    #TODO: try gtdb taxonomy
    elif taxonomy == 'GTDB':
        if user_taxonomy != None:
            taxonomy_db = GTDBTaxa(user_taxonomy, memory = True)
        else:
            taxonomy_db = GTDBTaxa()

    print('-Load taxonomy:  ')


    return taxonomy_db


def load_reftree(rtree=None, t=None, taxonomy_db=None):

    """
        Get reference tree (species tree) from input tree or user provide it
    """

    print('-Load reftree:   ')
    if rtree:
            reftree = PhyloTree(rtree)
    else:
        reftree = get_reftree(t, taxonomy_db)
    print('\t'+'Len reftree: ', len(reftree))

    taxonomy_db.annotate_tree(reftree,  taxid_attr="name")

    return reftree

def get_reftree(t,taxonomy_db):

    """
        Create reference tree (species tree) if user do not provide it
    """

    #print(' GETTING REFERENCE SPECIES TREE FROM INPUT GENE TREE')
    print('**Create referente species tree from gene tree')

    taxid_list = list()
    clase = str(taxonomy_db.__class__)
    if clase.split('.')[1] == 'ncbi_taxonomy':
        for leaf in t:
            taxid_list.append(int(leaf.name.split('.')[0]))

    elif clase.split('.')[1] == 'gtdb_taxonomy':
        for leaf in t:
            taxid_list.append(str(leaf.name))

    #Extract the smallest tree that connects all your query taxid
    reftree = taxonomy_db.get_topology(taxid_list)

    return reftree


def load_taxonomy_counter(reftree=None, user_taxonomy_counter=None):

    """
        Get number of species per each taxonomic level
        Get it from reference tree or user provide it
    """

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

def get_taxonomy_counter(reftree):

    """
        Create Taxonomy counter if user de not provide it
    """

    print('GETTING TAXONOMY COUNTER FROM INPUT GENE TREE')
    level2sp_mem = defaultdict(set)
    for l in reftree:
        lin = l.props.get('lineage')
        for tax in lin:
            level2sp_mem[str(tax)].add(l.name)

    level2sp_num = defaultdict()
    for taxid, mems in level2sp_mem.items():

        level2sp_num[taxid] = len(mems)

    return level2sp_num


## 2. Preanalysis  ##

def run_preanalysis(t, name_tree, taxonomy_db, rooting):

    print('\n'+'START PROCESSING TREE: '+name_tree)

    t = run_rooting(t, rooting)

    t = add_taxomical_annotation(t, taxonomy_db)

    # Total species in tree and total members in tree
    sp_set = set()
    total_mems_in_tree = set(t.leaf_names())

    for i, n in enumerate(t.traverse()):
        if not n.is_leaf:
            #Create an ID for each internal node
            n.name = '%s-%d' % (make_name(i), get_depth(n))

        else:
            sp_set.add(n.props.get('taxid'))


    SPTOTAL = len(sp_set)
    print('-Len tree: ', len(t))
    print('-Total species in tree:', SPTOTAL)

    t, props = run_clean_properties(t)
    tree_nw = get_newick(t, props)

    return tree_nw, sp_set, total_mems_in_tree, SPTOTAL, props


def run_rooting(t, rooting):

    if  rooting == "Midpoint":
        print('-Run midpoint rooting')
        root_mid = t.get_midpoint_outgroup()
        t.set_outgroup(root_mid)

    elif rooting == "MinVar":
        print('-Run MinVar rooting')

    t.dist = 0.01

    return t

def add_taxomical_annotation(t, taxonomy_db):

    # Parsing function used to extract species name from a node’s name.
    t.set_species_naming_function(lambda x: x.split('.')[0])

    # Adding additional information to any internal a leaf node (sci_name, taxid, named_lineage, lineage, rank)
    taxonomy_db.annotate_tree(t,  taxid_attr="species")

    return t

def make_name(i):

    """
        Create names for internal nodes
    """

    name = ''
    while i >= 0:
        name = chars[i % len(chars)] + name
        i = i // len(chars) - 1
    return name

def get_depth(node):

    """
        Get depth of internal nodes
        Depth = number nodes from root node to target node
    """
    depth = 0
    while node is not None:
        depth += 1
        node = node.up
    return depth


##  3. Outliers and Scores

def run_outliers_dup_score(t_nw, outliers_node, outliers_reftree, sp_loss_perc, so_cell_org, so_arq, so_bact, so_euk, taxonomy_counter, taxonomy_db, SPTOTAL, reftree, inherit_out):

    print('\n'+'DETECT OUTLIERS AND DUPS SCORE FUNCTION')
    print('\t'+'--Outliers threshodls:')
    print('\t'+'\t'+'--Node: '+ str(outliers_node))
    print('\t'+'\t'+'--Reftree: '+ str(outliers_reftree))
    print('\t'+'--Species losses percentage threshold: ' + str(sp_loss_perc))

    t = PhyloTree(t_nw)

    CONTENT = t.get_cached_content()

    t.add_prop('is_root', str('True'))

    long_leaves = detect_long_branches(t)


    for n in t.traverse("preorder"):

        n.del_prop('named_lineage')
        n.del_prop('_speciesFunction')

        #sci_name = clean_string(n.props.get('sci_name'))
        #n.props["sci_name"] = sci_name

        n.props["sci_name"] = n.props.get('sci_name')

        if n.is_leaf:
            n.add_prop('lca_node', n.props.get('taxid'))

        else:

            n.del_prop('common_name')

            #DETECT OUTLIERS
            sp_out = set()
            leaves_out = set()
            sp_in = set()
            sp_in_list = list()
            leaves_in = set()

            #Add outliers from upper nodes
            if inherit_out == 'Yes':
                sp_out_inherit = add_upper_outliers(n, CONTENT)
                sp_out.update(sp_out_inherit)
                # if n.up and n.up.props.get('sp_out'):
                    # sp_out_up = n.up.props.get('sp_out')
                    # for sp in sp_out_up:
                        # if (sp) in [leaf.props.get('taxid') for leaf in CONTENT[n]]:
                            # sp_out.add(sp)
            else:
                sp_out_inherit = set()


            # n_up_lin: Depth of parent node taxonomic level, needed for the outliers detection
            if n.up:
                n_up_lin = len(n.up.props.get('lineage'))
            else:
                n_up_lin = len(n.props.get('lineage').split('|'))


            #   Detect outliers in each of children nodes
            ch1 = n.children[0]
            sp_out.update(outliers_detection(ch1, outliers_node, outliers_reftree, CONTENT, taxonomy_counter, sp_out_inherit, n_up_lin,taxonomy_db))
            ch2 = n.children[1]
            sp_out.update(outliers_detection(ch2, outliers_node, outliers_reftree, CONTENT, taxonomy_counter, sp_out_inherit, n_up_lin,taxonomy_db))


            #Save leaves and species that are going to be in and out of species overlap
            # if len(sp_out) > 0:
                # #all_leafs = CONTENT[n]
                # for l in CONTENT[n]:
                    # if str(l.props.get('taxid')) in sp_out or l.name in long_leaves:
                        # leaves_out.add(l.props.get('name'))
                    # else:
                        # sp_in.add(str(l.props.get('taxid')))
                        # leaves_in.add(l.props.get('name'))

            # else:
                # for l in  CONTENT[n]:
                    # if l.name not in long_leaves :
                        # sp_in.add(str(l.props.get('taxid')))
                        # sp_in_list.append(str(l.props.get('taxid')))
                        # leaves_in.add(l.props.get('name'))



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

                else:
                    sp_in.add(str(l.props.get('taxid')))
                    sp_in_list.append(str(l.props.get('taxid')))
                    leaves_in.add(l.props.get('name'))

                    if l.name in ch1_leaf_names:
                        sp_ch1.add(str(l.props.get('taxid')))
                        leaves_ch1.add(l.props.get('name'))
                    elif l.name in ch2_leaf_names:
                        sp_ch2.add(str(l.props.get('taxid')))
                        leaves_ch2.add(l.props.get('name'))


            # for l in ch1:
                # if str(l.props.get('taxid')) not in sp_out and l.name not in long_leaves:
                    # sp_ch1.add(str(l.props.get('taxid')))
                    # leaves_ch1.add(l.props.get('name'))

            # for l in ch2:
                # if str(l.props.get('taxid')) not in sp_out and l.name not in long_leaves:
                    # sp_ch2.add(str(l.props.get('taxid')))
                    # leaves_ch2.add(l.props.get('name'))



            #   Re-calculate species overlap after detect outliers
            overlaped_spces = set(sp_ch1 & sp_ch2)

            if len(overlaped_spces)>0:
               # so_score = float(len(overlaped_spces) / len(all_spcs))
                so_score = float(len(overlaped_spces) / len(sp_in))
                n.add_prop('_overlap', list(overlaped_spces))
            else:
                so_score = 0.0


            #   Update last common ancestor, rank and lineage for the node after detect outliers
            update_taxonomical_props(n, sp_in, taxonomy_db)

            # if len(sp_in) > 0:
                # lca_node = get_lca_node(sp_in, taxonomy_db)
                # rank = clean_string(taxonomy_db.get_rank([lca_node])[lca_node])
                # lin_lca = taxonomy_db.get_lineage(lca_node)
                # n.add_prop('lineage', lin_lca)
                # n.add_prop('taxid', lca_node)
            # elif len(sp_in) == 0:
                # lca_node = 1
                # rank = clean_string(taxonomy_db.get_rank([lca_node])[lca_node])
                # lin_lca = taxonomy_db.get_lineage(lca_node)
                # n.add_prop('lineage', lin_lca)
                # n.add_prop('taxid', lca_node)


            #   Save properties
            # n.add_prop('lca_node', lca_node)
            # lca_node_name = taxonomy_db.get_taxid_translator([lca_node])[lca_node]

            # n.add_prop('lca_node_name', lca_node_name.replace(":", ""))

            # n.props['sci_name'] = n.props.get('lca_node_name')
            #n.add_prop('rank', rank)

            n.add_prop('_sp_in', sp_in)
            n.add_prop('len_sp_in', len(sp_in))
            n.add_prop('_sp_in_ch1', sp_ch1)
            n.add_prop('_sp_in_ch2', sp_ch2)
            n.add_prop('_ch1_name', ch1_name)
            n.add_prop('_ch2_name', ch2_name)
            n.add_prop('_leaves_in', leaves_in)
            n.add_prop('total_leaves', len(n))
            n.add_prop('len_leaves_in', len(leaves_in))
            n.add_prop('_leaves_ch1',list(leaves_ch1))
            n.add_prop('_leaves_ch2',list(leaves_ch2))
            if len(sp_out) == 0:
                n.add_prop('sp_out', ['None'])
            else:
                n.add_prop('sp_out', list(sp_out))
            n.add_prop('_leaves_out', list(leaves_out))
            n.add_prop('so_score', so_score)


            #   Load_node_scores add properties: score1, score2 and inpalalogs_rate

            #load_node_scores(n, SPTOTAL)

            inparalogs_rate = get_inparalogs_rate(sp_in_list)
            n.add_prop('inparalogs_rate', str(inparalogs_rate))

            score1 = get_score1(len(sp_in), SPTOTAL)
            n.add_prop('score1', score1)

            score2 = get_score2(len(sp_in), len(leaves_in))
            n.add_prop('score2', score2)

            dup_score = get_dup_score(n)
            n.add_prop('dup_score', dup_score)

            call_species_lost(n, reftree)

            call_lineage_lost(n, reftree)


            #   Based on species overlap threshold define duplication, false duplication and speciation nodes
            lin_lca = n.props.get('lineage')
            lca_node = n.props.get('lca_node')
            if 2 in lin_lca:
                so_2_use = so_bact
            elif 2759 in lin_lca:
                so_2_use = so_euk
            elif 2157 in lin_lca:
                so_2_use = so_arq
            elif lca_node == 131567:
                so_2_use = so_cell_org
            else:
                so_2_use = 0.2

            if so_score >= so_2_use:
                if float(n.props.get('species_losses_percentage')[0]) > sp_loss_perc and float(n.props.get('species_losses_percentage')[1]) > sp_loss_perc:
                    n.add_prop('evoltype_2', 'FD')
                else:
                    n.add_prop('evoltype_2', 'D')
            else:
                n.add_prop('evoltype_2', 'S')


    t, props = run_clean_properties(t)

    return t

def detect_long_branches(t):

    """
        Detect leaves with long lenght branchs
        Long lenght branchs is when a branch is 10 times bigger than the mean for the whole tree
    """

    length_leaves = list()
    for l in t:
        length_leaves.append(l.dist)

    mean_length = np.mean(length_leaves)

    long_leaves = set()
    for l in t:
        if l.dist > mean_length*10:
            long_leaves.add(l.name)

    return long_leaves

def add_upper_outliers(n, CONTENT):
    sp_out_inherit = set()
    if n.up and n.up.props.get('sp_out'):
        sp_out_up = n.up.props.get('sp_out')
        for sp in sp_out_up:
            if (sp) in [leaf.props.get('taxid') for leaf in CONTENT[n]]:
                sp_out_inherit.add(sp)
    return sp_out_inherit

def outliers_detection(n, outliers_node, outliers_reftree, CONTENT_, taxonomy_counter, sp_out_up, n_up_lin, taxonomy_db):


    #Save sp in node,  it doesnt matter that for the same sp there are several seqs
    sp_in_node = set()

    #For each taxonomic level in the leaves's lineage, save sp under that level
    sp_per_level = defaultdict(set)

    for l in CONTENT_[n]:
        if l.props.get('taxid') not in sp_out_up:
            sp_in_node.add(l.props.get('taxid'))

            for tax in l.props.get('lineage').split('|'):
                sp_per_level[tax].add(str(l.props.get('taxid')))



    # Detect outliers for each taxid from lineages

    # best_rep = defaultdict()
    # best_rep_num = defaultdict()
    # best_rep_node = defaultdict()
     #best_rep_val = 0.0

    tax_out = list()
    sp2remove = set()

    for tax, num in sp_per_level.items():

        #Only visit taxononomic levels two levels below parent node taxonomic level
        if len(taxonomy_db.get_lineage(tax)) <= int(n_up_lin)+2:

            #Percentage of lineage in the node
            per_node = len(num) / len(sp_in_node)

            #Percentage of lineage in reference species tree
            per_tree = len(num) / taxonomy_counter[str(tax)]

            # If % of sp in the node at that taxonomic level is below 1%  and % of sp in reftree is below 5%, remove that species, that means:
                # 1. The taxonomic level has to be rare in the node (a few porifera species inside a node with all bilateria)
                # 2. Also there has to be few representatives from reftree at that taxonomic level, those leaves will be remove
                #   (in the node there are 2 porifera species but in reftree there are 2000 porifera species)
                # 3. Taxonomic level that are rare in reftree, will be preserved

            if outliers_node == 0.0:
                if per_tree < outliers_reftree:# and per_Node < outliers_node:
                    sp2remove.update(sp_per_level[tax])
                    tax_out.append(tax)
            elif outliers_reftree == 0.0:
                if per_node < outliers_node:# and per_Node < outliers_node:
                    sp2remove.update(sp_per_level[tax])
                    tax_out.append(tax)

            elif outliers_reftree != 0.0 and outliers_node != 0.0:
                if per_tree < outliers_reftree and per_node < outliers_node:
                    sp2remove.update(sp_per_level[tax])
                    tax_out.append(tax)


    n.add_prop('outliers_tax', tax_out)
    return sp2remove

def update_taxonomical_props(n, sp_in, taxonomy_db):

    if len(sp_in) > 0:
        lca_node = get_lca_node(sp_in, taxonomy_db)
        rank = clean_string(taxonomy_db.get_rank([lca_node])[lca_node])
        lin_lca = taxonomy_db.get_lineage(lca_node)
        n.add_prop('lineage', lin_lca)
        n.add_prop('taxid', lca_node)
    elif len(sp_in) == 0:
        lca_node = 1
        rank = clean_string(taxonomy_db.get_rank([lca_node])[lca_node])
        lin_lca = taxonomy_db.get_lineage(lca_node)
        n.add_prop('lineage', lin_lca)
        n.add_prop('taxid', lca_node)

    n.add_prop('lca_node', lca_node)
    lca_node_name = taxonomy_db.get_taxid_translator([lca_node])[lca_node]
    n.add_prop('lca_node_name', lca_node_name.replace(":", ""))
    n.props['sci_name'] = n.props.get('lca_node_name')
    n.add_prop('rank', rank)

def get_lca_node(sp_list, taxonomy_db):

    """
        Get Last Common Ancestor (lca) from a set of species
    """

    lineages = OrderedCounter()
    nspecies = 0
    for sp in sp_list:
        nspecies += 1
        lineages.update(taxonomy_db.get_lineage(sp))

    lca = [l for l, count in lineages.items() if count == nspecies]
    if not lca:
        lca = ['Unk']

    return lca[-1]

def get_inparalogs_rate(leaf_targets):

    dups_per_sp = Counter(leaf_targets)
    inparalogs_rate = np.median(list(dups_per_sp.values()))

    return inparalogs_rate

def get_score1(sp_in, SPTOTAL):

    score1 = float(sp_in/SPTOTAL)
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

def call_lineage_lost (node, reftree):
    """
        Call count_lineage_losses() to count how many lineages are missing in each children node
        TODO: CHANGE NAME FUNCTION
            if sentence, change dup_score maybe better with species overlap??? or maybe just delte the if and calculate count_lienage_losses for all nodes (duplication and speciacion nodes)
    """

    sp1 = node.props.get('_sp_in_ch1')
    sp2 = node.props.get('_sp_in_ch2')

    if float(node.props.get('dup_score')) > 0.0:
        losses1, lin_losses1 = count_lineage_lost(expected_sp=sp1|sp2, found_sp=sp1, reftree=reftree)
        losses2, lin_losses2 = count_lineage_lost(expected_sp=sp1|sp2, found_sp=sp2, reftree=reftree)
        node.add_prop('num_lineage_losses', [losses1, losses2])
        node.add_prop('_taxid_lineage_losses', [lin_losses1, lin_losses2])

def count_lineage_lost(expected_sp, found_sp, reftree):

    """
        Find how many lineages are missing
        Lineage: all species under internal node
        If all species under internal node are lost, then losses +=1
    """

    def is_leaf_2(_n):
        if not _n.children:
            return True
        elif len(found_sp & set(_n.leaf_names())) == 0:
            return True
        else:
            return False


    losses = 0
    lca_lin_lost = set()

    if len(expected_sp) == 1:
        return losses, lca_lin_lost

    elif len(expected_sp) > 1:

        expected_reftree = reftree.common_ancestor(expected_sp)

        for leaf in expected_reftree.leaves(is_leaf_fn=is_leaf_2):
            losses +=1

            # Find missinf lineages are very slow
            #anc = get_lca_node(leaf.get_leaf_names())
            #lca_lin_lost.add(anc)

        return losses, lca_lin_lost

def call_species_lost(node, reftree):
    """
        Call count_species_losses() for each children node
        TODO: if sentence?? delete?? change dup_score for species overlap???
    """

    sp1 = node.props.get('_sp_in_ch1')
    sp2 = node.props.get('_sp_in_ch2')

    if float(node.props.get('dup_score')) > 0.0:
        losses1, per_loss1 = count_species_lost(expected_sp=sp1|sp2, found_sp=sp1, reftree=reftree)
        losses2, per_loss2 = count_species_lost(expected_sp=sp1|sp2, found_sp=sp2, reftree=reftree)
        node.add_prop('species_losses', [losses1, losses2])
        node.add_prop('species_losses_percentage', [per_loss1, per_loss2])

    else:
        node.add_prop('species_losses', [0, 0])
        node.add_prop('species_losses_percentage', [0.0, 0.0])

def count_species_lost(expected_sp, found_sp, reftree):

    """
        Calculate number of species_losses -> number of expected species - number of found species
        Calculate percentage of species_losses_percentage -> losses / number of expected species
    """

    if len(expected_sp) == 1:
        return 0, 0

    root = reftree.common_ancestor(expected_sp)

    losses = len(root) - len(found_sp)
    per_loss = losses / len(root)

    return int(losses), float(per_loss)


##  4. Detect Duplications and Core-OGs ##

def run_dups_and_ogs(t, outliers_node, outliers_reftree, sp_loss_perc, so_cell_org, so_arq, so_bact, so_euk, taxonomy_db, total_mems_in_tree):

    """
        Select high-quality duplication nodes
    """

    #ogs_info, save info of nodes that are OGs
    ogs_info = defaultdict(dict)

    #taxid_dups_og, set with the lca of the nodes that are OGs
    taxid_dups_og = set()
    total_mems_in_ogs = set()

    #Traverse tree to find the nodes that are "good" duplications and generate OGs.
    print('\n'+'ITER DUPS')
    print('\t'+'--species overlap used:')
    print('\t'+'\t'+'--Cell Org: '+ str(so_cell_org))
    print('\t'+'\t'+'--Euk: '+ str(so_euk))
    print('\t'+'\t'+'--Bact :' + str(so_bact))
    print('\t'+'\t'+'--Arq :' + str(so_arq))


    for node in t.traverse("preorder"):
        #node.del_prop('_speciesFunction')
        lin2check = node.props.get('lineage')

        lca_node = node.props.get('lca_node')

        if isinstance(lin2check, str):
            lin2check = lin2check.split('|')

        if 2 in lin2check:
            so_2_use = so_bact
        elif 2759 in lin2check:
            so_2_use = so_euk
        elif 2157 in lin2check:
            so_2_use = so_arq
        elif lca_node == 131567:
            so_2_use = so_cell_org
        else:
            so_2_use = 0.2

        node.add_prop('min_sp_overlap', so_2_use)



        #Requisites for high- quality dups:
            #species overlap higher than min threshold
            #more than one leave and more than one specie
            #no more dups at the same taxonomic level below the node

        if not node.is_leaf and node.props.get('evoltype_2') == 'D' and so_2_use <= float(node.props.get('so_score')) \
        and len(node.props.get('_leaves_in')) >1 and len(node.props.get('_sp_in')) > 1 :


            dups_under_node = []
            for n in node.search_nodes(evoltype_2='D'):

                if n.name!= node.name:
                    dups_under_node.append(n)

            #There are more dups under the node
            if len(dups_under_node) > 0:

                lca_target = node.props.get('lca_node')

                #Save Dups under child 1 and child2 that have the same lca_node
                dups_under_ch1 = list(node.children[0].search_nodes(evoltype_2='D', lca_node=lca_target))
                dups_under_ch2 = list(node.children[1].search_nodes(evoltype_2='D', lca_node=lca_target))

                save_dups_ch1 = defaultdict()
                save_dups_ch2 = defaultdict()

                # Check that dups under child1 and child2 (that have the same lca) fit out requirements : species overlap min requirement,
                # more than 1 leaves and more than 1 species
                #if len(dups_under_ch1)> 0:
                for n_ in  dups_under_ch1:
                    if float(n_.props.get('so_score'))>= so_2_use \
                    and len(n_.props.get('_leaves_in')) >1 and len(n_.props.get('_sp_in'))> 1 :
                        root2node = node.get_distance(node, n_, topological=True)
                        save_dups_ch1[n_.name] = root2node


                #if len(dups_under_ch2)> 0:
                for n_ in  dups_under_ch2:
                    if float(n_.props.get('so_score'))>= so_2_use  \
                    and len(n_.props.get('_leaves_in'))>1 and len(n_.props.get('_sp_in'))> 1 :
                        root2node = node.get_distance(node, n_, topological=True)
                        save_dups_ch2[n_.name] = root2node

                #If dups under child1 do not achieve our requirement, then child1 is OG
                if len(save_dups_ch1) == 0:
                    annotate_dups_ch(ogs_info, total_mems_in_ogs, taxid_dups_og,node, 'ch1', taxonomy_db)

                #If dups under child2 do not achieve our requirement, then child2 is OG
                if len(save_dups_ch2) == 0 :
                    annotate_dups_ch(ogs_info, total_mems_in_ogs, taxid_dups_og,node, 'ch2', taxonomy_db)


            # No more dups under node with the same lca_node, both ch nodes are OGs
            elif len(dups_under_node) == 0:
                annotate_dups_ch(ogs_info, total_mems_in_ogs, taxid_dups_og,node, 'ch1', taxonomy_db)
                annotate_dups_ch(ogs_info, total_mems_in_ogs, taxid_dups_og,node, 'ch2', taxonomy_db)

    t, props = run_clean_properties(t)

    return t, total_mems_in_ogs, ogs_info, taxid_dups_og

def annotate_dups_ch(ogs_info, total_mems_in_ogs, taxid_dups_og, node, ch_node, taxonomy_db):

    """
    Add props and save info about the children node that is OG
    node is the duplication node
    target node is the node that is OG
    """

    if ch_node == 'ch1':
        og_name_ch = node.props.get('_ch1_name')
        og_ch_mems = node.props.get('_leaves_ch1')
        sp_ch = node.props.get('_sp_in_ch1')

        # print(len(node.search_nodes(name=og_name_ch)))
        # for n in node.search_nodes(name=og_name_ch):
            # print(n)
            # target_

        target_node = next(node.search_nodes(name=og_name_ch)) #[0]
        #support_target_node = float(target_node.props.get('support'))

    elif ch_node == 'ch2':
        og_name_ch = node.props.get('_ch2_name')
        og_ch_mems = node.props.get('_leaves_ch2')
        sp_ch = node.props.get('_sp_in_ch2')
        target_node = next(node.search_nodes(name=og_name_ch)) #[0]

    if  len(sp_ch) > 1 and len(og_ch_mems) > 1: #support_target_node > 50.0 and
        target_node.add_prop('node_is_og', 'True')
        target_node.add_prop('lca_dup', node.props.get('lca_node'))
        target_node.add_prop('so_score_dup', node.props.get('so_score'))
        target_node.add_prop('dup_lineage', list(taxonomy_db.get_lineage(node.props.get('lca_node'))))
        target_node.add_prop('_mems_og', '|'.join(list((og_ch_mems))))
        target_node.add_prop('_dup_node_name', node.props.get('name'))

        taxid_dups_og.add(node.props.get('lca_node'))
        node.add_prop('node_create_og', 'True')
        sp_in_og = target_node.props.get('_sp_in')
        sp_in_anc_og = node.props.get('_sp_in')

        total_mems_in_ogs.update(set(og_ch_mems))
        ogs_info[og_name_ch]["name_dup"] = node.props.get("name")
        ogs_info[og_name_ch]["lca_dup"] = node.props.get("lca_node")
        ogs_info[og_name_ch]["lca_dup_lineage"] = node.props.get("lineage")
        ogs_info[og_name_ch]["tax_scope_og"] = target_node.props.get("lca_node")
        ogs_info[og_name_ch]["ogs_mems"] = target_node.props.get("_mems_og")
        ogs_info[og_name_ch]["total_leaves"] =target_node.props.get('total_leaves')
        ogs_info[og_name_ch]["num_sp_OGs"] = target_node.props.get('len_sp_in')
        ogs_info[og_name_ch]["num_sp_out"] = target_node.props.get('sp_out','0')
        ogs_info[og_name_ch]['recovery_mems'] = set()


##  5. Add info about Core-OGs up and down  ##

def add_ogs_up_down(t, ogs_info):

    """
        Traverse the tree, for each node,  detect ogs that are up or below
    """

    for node in t.traverse('preorder'):

        node_name = node.name

        #Detect OGs below node
        ogs_down = set()
        dups_down = list()
        for n in node.search_nodes(node_is_og="True"):
            if node_name != n.name:
                ogs_down.add(n.name)
                dups_down.append(n.up.name)


        if node_name != n.name:
            ogs_down.add(n.name)
            dups_down.append(n.up.name)

        if len(ogs_down) > 0:
            if node.props.get('node_is_og'):
                ogs_info[node_name]['ogs_down'] = ','.join(list(ogs_down))
                ogs_info[node_name]['dups_down'] = dups_down
            node.add_prop('_ogs_down',ogs_down)
            node.add_prop('_dups_down', dups_down)
        else:
            if node.props.get('node_is_og'):
                ogs_info[node_name]['ogs_down'] = '-'
                ogs_info[node_name]['dups_down'] = '-'


        #Detect OGs up
        ogs_up = set()
        dups_up = list()
        ogs_up, dups_up = check_nodes_up(node)

        if len(ogs_up) > 0:
            if node.props.get('node_is_og'):
                ogs_info[node_name]['ogs_up'] = ','.join(list(ogs_up))
                ogs_info[node_name]['dups_up'] = dups_up
            node.add_prop('_ogs_up', ogs_up)
            node.add_prop('_dups_up', dups_up)
        else:
            if node.props.get('node_is_og'):
                ogs_info[node_name]['ogs_up'] = '-'
                ogs_info[node_name]['dups_up'] = '-'

    #t, props = run_clean_properties(t)

    return t, ogs_info


##  8.Annotate root ##

def annotate_root(ogs_info, t, total_mems_in_tree, sp_set, total_mems_in_ogs, recovery_seqs, outliers_node, outliers_reftree, sp_loss_perc, so_cell_org, so_arq, so_bact, so_euk, inherit_out, taxonomy_db):

    """
    Add properties in root
    Add info about the root in ogs_info dict
    """

    # node.add_prop('so_score_dup', node.props.get('so_score'))
    # node.add_prop('dup_lineage', list(taxonomy_db.get_lineage(node.props.get('lca_node'))))
    # node.add_prop('_mems_og', total_mems_in_tree)
    # node.add_prop('_dup_node_name', node.props.get('name'))



    # Save parameters as property

    parameters_str = '|'.join(["outliers_node@"+str(outliers_node), "outliers_reftree@"+str(outliers_reftree), "sp_loss_perc@"+str(sp_loss_perc),
                                "so_cell_org@"+str(so_cell_org), "so_arq@"+str(so_arq), "so_bact@"+str(so_bact), "so_euk@"+str(so_euk), "inherit_out@"+str(inherit_out)])

    t.add_prop('parameters', parameters_str)
    # name = node.props.get('name')
    # node.add_prop("outliers_node", outliers_node)
    # node.add_prop("outliers_reftree", outliers_reftree)
    # node.add_prop("sp_loss_perc", sp_loss_perc)
    # node.add_prop("so_cell_org", so_cell_org)
    # node.add_prop("so_arq", so_arq)
    # node.add_prop("so_bact", so_bact)
    # node.add_prop("so_euk", so_euk)
    # node.add_prop("inherit_out", inherit_out)


    # Save general results as properties, needed for web

    seqs_out_ogs = total_mems_in_tree.difference(total_mems_in_ogs)

    general_results_str = '|'.join(["Tree_name@"+t.props.get('name'), "Total_seqs@"+str(len(total_mems_in_tree)), "Total_species@"+str(len(sp_set)),
                                    "Seqs_in_OGs@"+str(len(total_mems_in_ogs)), "Recovery_seqs@"+str(len(recovery_seqs)), "Seqs_out_OGs@"+str(len(seqs_out_ogs)), "Num_OGs@"+str(len(ogs_info)) ])

    t.add_prop('general_result', general_results_str)


    # Save taxlevel 2 OGs info, needed for web

    lca_all_dups = set()
    taxlev2ogs = defaultdict(set)
    taxlev2mems = defaultdict(set)
    for og, info in ogs_info.items():
        lca_all_dups.add(info['lca_dup'])
        taxlev2ogs[info['lca_dup']].add(og)
        taxlev2mems[info['lca_dup']].update(set(info["ogs_mems"]))

    taxlev_list = []
    for taxlev, og in taxlev2ogs.items():
        sci_name = taxonomy_db.get_taxid_translator([taxlev])[taxlev]
        sci_name_taxid = sci_name+'_'+str(taxlev)
        ogs_str = '_'.join(list(og))
        num_mems = len(taxlev2mems[taxlev])

        tax_str = '|'.join([sci_name_taxid, ogs_str,str(num_mems)])

        taxlev_list.append(tax_str)

    taxlev_str = '@'.join(taxlev_list)
    t.add_prop('taxlev2ogs', taxlev_str)


    # If no OGs in tree, or the taxonomic level of OGs are under taxonomic level of the root
    # Root is  OG

    name = t.props.get("name")

    if len(ogs_info.keys()) == 0 or t.props.get('lca_node') not in lca_all_dups:
        t.add_prop('root_is_og', 'True')
        ogs_info[name]["name_dup"] = t.props.get("name")
        ogs_info[name]["lca_dup"] = t.props.get("lca_node")
        ogs_info[name]["lca_dup_lineage"] = t.props.get("lineage")
        ogs_info[name]["tax_scope_og"] = t.props.get("lca_node")
        ogs_info[name]["ogs_mems"] = t.props.get("_mems_og")
        ogs_info[name]["total_leaves"] = t.props.get('total_leaves')
        ogs_info[name]["num_sp_OGs"] = t.props.get('len_sp_in')
        ogs_info[name]["num_sp_out"] = t.props.get('sp_out','0')
        ogs_info[name]['recovery_mems'] = set()
        #total_mems_in_ogs.update(set(node.props.get("_mems_og")))
        ogs_info[name]['ogs_up'] = '-'
        ogs_info[name]['dups_up'] = '-'
        ogs_info[name]['ogs_down'] = '-'
        ogs_info[name]['dups_down'] = '-'


    # elif  node.props.get('lca_node') not in lca_all_dups:
        # node.add_prop('root_is_og', 'True')
        # ogs_info[name]["name_dup"] = node.props.get("name")
        # ogs_info[name]["lca_dup"] = node.props.get("lca_node")
        # ogs_info[name]["lca_dup_lineage"] = node.props.get("lineage")
        # ogs_info[name]["tax_scope_og"] = node.props.get("lca_node")
        # ogs_info[name]["ogs_mems"] = node.props.get("_mems_og")
        # ogs_info[name]["total_leaves"] = node.props.get('total_leaves')
        # ogs_info[name]["num_sp_OGs"] = node.props.get('len_sp_in')
        # ogs_info[name]["num_sp_out"] = node.props.get('sp_out','0')
        # ogs_info[name]['recovery_mems'] = set()
        # #total_mems_in_ogs.update(set(node.props.get("_mems_og")))
        # ogs_info[name]['ogs_up'] = '-'
        # ogs_info[name]['dups_up'] = '-'
        # ogs_info[name]['ogs_down'] = '-'
        # ogs_info[name]['dups_down'] = '-'


    ogs_info[t.props.get('name')]['is_root'] = 'True'

    t.add_prop("OGD_annot", True)


##  9. Flag seqs out OGs ##

def flag_seqs_out_og(t, total_mems_in_ogs, total_mems_in_tree):

    """
        Flags seqs that do no belong to any OG
    """
    seqs_out_og = total_mems_in_tree.difference(total_mems_in_ogs)
    for leaf in t:
        if leaf.name in seqs_out_og:
            leaf.add_prop('seq_out_og', "true")

    t, props = run_clean_properties(t)

    return t



####    FUNCTIONS FOR RECOVERING SEQUENCES PIPELINE ####

def recover_sequences(tree, t, alg, total_mems_in_tree, total_mems_in_ogs, name_tree, path_out, ogs_info, mode):

    name_fam = os.path.basename(tree).split('.')[0]
    pathout = path_out+'/'+name_fam+'_aln_hmm'
    if not os.path.exists(pathout):
        os.mkdir(pathout)
    fasta = alg

    # 1. Write a fasta file
    run_write_fastas(fasta, name_tree, pathout, ogs_info, total_mems_in_ogs, total_mems_in_tree, mode)

    if mode == "regular":  # regular mode: Re-aling sequences
        run_alignments(pathout)
    else:
        pass  # we just go with what we had

    # 2. Create HMM file for each Core-OGs fasta file and Build HMM DB
    run_create_hmm_og(pathout)

    # 3. Run Hmmscan to assign seqs out Core-OGs to Core-OGs
    tblfile = run_hmmscan(pathout)

    # 4. For each seq out Core-OGs, select best Core-OGs
    best_match = get_best_match(tblfile)

    # 5. Create a dict with all the OGs that have recover some seq
    og_info_recovery, recovery_seqs = expand_hmm(best_match, ogs_info)#, total_recovery_seqs)

    # 6. Update og_info with all the seqs taht have been recovered
    og_info_updated = update_og_info(ogs_info, og_info_recovery)

    total_mems_in_ogs.update(recovery_seqs)

    # 7. Update tree with all the seqs taht have been recovered
    t = update_tree(t, og_info_recovery)

    # 8. write_ogs_info for recovery
    write_ogs_info(og_info_updated, 'recovery', path_out)

    # 9. write_best_match
    write_best_match(best_match, path_out)

    output_filename = path_out+'/'+name_fam+'.tar.gz'

    # 10. Tar the folder with all the fasta file, HMMs, etc
    make_tarfile(output_filename, pathout)
    shutil.rmtree(pathout)

    return recovery_seqs, best_match


def write_outog_seqs(fasta, diff, not_og_fasta):

    """
        Recovery pipeline
        Write fasta file with seqs that do not belong to any OG
    """

    with open(not_og_fasta, 'w') as  f_out:
        for name_seq in diff:
            aa = fasta.get_seq(name_seq)
            clean_aa = aa.replace('-','')
            f_out.write('>'+name_seq+'\n'+clean_aa+'\n')


def write_og_seqs_regular_mode(fasta, ogs_info, path_out):

    """
        Write one fasta file per each OG
        In regular mode, sequences will be realing, so gaps are removed
    """

    for name_og, info in ogs_info.items():
        if not 'is_root' in ogs_info[name_og].keys():
            lca = str(info['lca_dup'])
            with open(path_out+'/'+name_og+'_'+lca+'.raw_fasta.faa', 'w') as f_out:
                if isinstance(info["ogs_mems"], str):
                    list_mems = (info["ogs_mems"]).split('|')
                else:
                    list_mems = (info["ogs_mems"])
                if len(list_mems) >0:
                    for m in list_mems:
                        aa = fasta.get_seq(m)
                        clean_aa = aa.replace('-','')
                        f_out.write('>'+m+'\n'+clean_aa+'\n')


def write_og_seqs_fast_mode(fasta, ogs_info, path_out):

    """
        Write one fasta file per each OG
        In fast mode, sequences wont be realing, so gaps are keept
    """

    for name_og, info in ogs_info.items():
        if not 'is_root' in ogs_info[name_og].keys():
            lca = str(info['lca_dup'])
            with open(path_out+'/'+name_og+'_'+lca+'.aln', 'w') as f_out:
                list_mems = info["ogs_mems"].split('|')
                if len(list_mems) >0:
                    for m in list_mems:
                        aa = fasta.get_seq(m)
                        f_out.write('>'+m+'\n'+aa+'\n')


def run_write_fastas(fasta, name_tree, path_out, ogs_info, total_mems_in_ogs, total_mems_in_tree, mode):

    """
        Function calls  write_outog_seqs() AND write_og_seqs_regular_mode() OR write_og_seqs_fast_mode()
    """

    print('-Fastas')

    fasta = SeqGroup(fasta)

    diff = (total_mems_in_tree.difference(total_mems_in_ogs))

    name_fam = name_tree.split('.')[0]
    not_og_fasta = path_out+'/'+name_fam+'.notOG.faa'

    write_outog_seqs(fasta, diff, not_og_fasta)

    if mode == "regular":
        write_og_seqs_regular_mode(fasta, ogs_info, path_out)
    elif mode == "fast":
        write_og_seqs_fast_mode(fasta, ogs_info, path_out)

    return


def run_mafft(raw_fasta, path_out, name_og):

    """
        Run mafft software for OGs < 100 seqs
    """

    aln_fasta = path_out+'/'+name_og+'.aln'
    subprocess.run(("mafft --auto --anysymbol %s > %s" %(raw_fasta, aln_fasta)), shell = True)


def run_famsa(raw_fasta,path_out, name_og):

    """
        Run famsa software for OGs > 100 seqs
    """

    aln_fasta = path_out+'/'+name_og+'.aln'
    subprocess.run(("famsa %s %s" %(raw_fasta, aln_fasta)), shell = True)


def run_alignments(path_out):

    """
        Function call run_mafft() or run_famsa() depend of the size of OG
    """

    total_aln = glob.glob(path_out+'/*.raw_fasta.faa')
    for aln in total_aln:
        aln_seq = SeqGroup(aln)
        name_og = os.path.basename(aln).split('.')[0]

        if len(aln_seq) <= 100:
            run_mafft(aln,path_out, name_og)
        else:
            run_famsa(aln, path_out, name_og)


def run_create_hmm_og(path):

    """
        Build HHMM file for each OG-aln fasta
        Create HMM DB with all the HMM files

        TODO: split in 2¿?
    """

    aln_list = glob.glob(path+'/*.aln')
    for aln in aln_list:
        out_hmm = aln+'.hmm'
        subprocess.run(("hmmbuild %s %s" %(out_hmm, aln)), shell = True)

    hmm_list = glob.glob(path+'/*.hmm')
    hmm_db = path+'/hmm_db.hmm'
    for hmm in hmm_list:
        subprocess.run(("cat %s >>%s" %(hmm, hmm_db)), shell = True)


    subprocess.run(("hmmpress %s" %(hmm_db)), shell = True)


def run_hmmscan(path):

    """
        Run Hmmscan
    """

    hmm_db = path+'/hmm_db.hmm'
    outfile = path+'/result_hmmscan.tsv'
    tblfile = path+'/tblout.tsv'
    domtblfile = path+'/domfile.tsv'
    seqs_not_og = glob.glob(path+'/*.notOG.faa')[0]
    print(("hmmscan -o %s %s %s" %(outfile, hmm_db, seqs_not_og)))
    subprocess.run(("hmmscan -o %s --tblout %s --domtblout %s  %s %s" %(outfile, tblfile, domtblfile, hmm_db, seqs_not_og)), shell = True)
    return tblfile


def get_best_match(tblout):

    """
        For each seq,  get best hit from Hmmscan
    """

    best_match = defaultdict(dict)

    with open(tblout) as f:
        for line in f:
            if not line.startswith('#'):
                line = re.sub(' +','\t',line)
                info = line.split('\t')
                name_og = info[0]
                name_seq = info[2]
                score = float(info[5])

                k_ = str()
                score_ = float()
                if name_seq in best_match.keys():
                    for k, val in best_match[name_seq].items():
                        k_ = k
                        score_ = val

                    if score > score_:
                        best_match[name_seq].pop(k_)
                        best_match[name_seq][name_og] = score
                else:
                    best_match[name_seq][name_og] = score


    return(best_match)


def expand_hmm(best_match, ogs_info):#, total_recovery_seqs):

    """
        Create a dict with all the OGs that have recover some seq

        TODO: change name
    """

    total_recovery_seqs = set()
    og_info_recovery = defaultdict(set)

    #total_recovery_seqs = set()

    for seq_name, best_og in best_match.items():
        for k in best_og.keys():
            best_og_name = k.split('_')[0]
        og_info_recovery[best_og_name].add(seq_name)
        total_recovery_seqs.add(seq_name)


    return og_info_recovery, total_recovery_seqs


def update_taxlevel2ogs(glob_taxlev2ogs, og_info_recovery, glob_og_info) :

    """
        update taxlevel2og after recovered step
        Needed for web
    """

    for tax, info in glob_taxlev2ogs.items():
        for og in info['ogs_names']:
            if og in og_info_recovery.keys():
                glob_taxlev2ogs[tax]['mems'].update(set(og_info_recovery[og]))

                ogs_up = glob_og_info[og]['ogs_up'].split(',')
                for og_up in ogs_up:

                    if og_up != '-':
                        lca_tax = glob_og_info[og_up]['lca_dup']
                        glob_taxlev2ogs[lca_tax]['mems'].update(set(glob_og_info[og_up]['ogs_mems']))

    return(glob_taxlev2ogs)


def update_og_info(og_info, og_info_recovery):

    """
        Add recovered seqs to respective OGs
    """

    for og, info in og_info.items():
        if 'is_root' not in info.keys():
            if isinstance(info['ogs_mems'], str):
                str2set = set(info['ogs_mems'].split('|') )
            else:
                str2set = set(info['ogs_mems'])
            info['ogs_mems'] = str2set


    for og, recover_seqs in og_info_recovery.items():

        if 'ogs_mems' in og_info[og].keys():
            og_info[og]['ogs_mems'].update(recover_seqs)
            og_info[og]['recovery_mems'].update(recover_seqs)

        else:
            og_info[og]['ogs_mems'] = recover_seqs
            og_info[og]['recovery_mems'] = recover_seqs

        if 'ogs_up' in og_info[og] and og_info[og]['ogs_up'] != '-':
            for og_up in (og_info[og]['ogs_up']).split(','):

                og_info[og_up]['ogs_mems'].update(recover_seqs)
                og_info[og_up]['recovery_mems'].update(recover_seqs)



    return og_info


def update_tree(t, og_info_recovery):

    """
        Update tree properties after recovery pipeline
    """

    for node2update in t.traverse('preorder'):

        if node2update.name in og_info_recovery.keys():

            recover_seq = og_info_recovery[node2update.name]
            if isinstance(node2update.props.get('_mems_og'), str):
                set2update = (set(node2update.props.get('_mems_og').split('|')))
            else:
                set2update = (set(node2update.props.get('_mems_og')))

            for s in recover_seq:
                set2update.add(s)

            node2update.add_prop('_mems_og', '|'.join(list(set2update)))

            node2update.add_prop('recover_seqs', recover_seq )


            if '_ogs_up' in node2update.props.keys():
                set_ogs_up =  ((node2update.props.get('_ogs_up')))
                for l in set_ogs_up:
                    node_up2update = next(t.search_nodes(name=l)) # list(t.search_nodes(name=l)[0]
                    if isinstance(node_up2update.props.get('_mems_og'), str):
                        node_up_seqs = set(node_up2update.props.get('_mems_og').split('|'))
                    else:
                        node_up_seqs = set(node_up2update.props.get('_mems_og'))

                    node_up_seqs.update(recover_seq)

                    node_up2update.add_prop('_mems_og', '|'.join(list(node_up_seqs)))


    return t


def make_tarfile(output_filename, source_dir):

    """
        Create tar.gz file
    """

    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))



####    EMAPPER APP     ####

def annotate_with_emapper(t, alg, path_out):
    """"
        Add to the leaves of tree t the information coming from emapper
        and return the annotated tree.
    """
    # In the main program, we are actually only interested in annotations
    # the pfam domains, but we have to take it all from emapper.
    path2pfam_table, path2main_table = run_emapper(alg, path_out)

    t = annot_tree_pfam_table(t, path2pfam_table, alg)

    t = annot_tree_main_table(t, path2main_table)

    # TODO: Check if we really do not need this (because it
    # is done in the main function already?):
    #   t, all_props = run_clean_properties(t)
    #   run_write_post_tree(t, name_tree, path_out, all_props)

    return t


def run_emapper(alg_fasta,path_out):

    """
        Run eggnog-mapper
    """

    fasta = SeqGroup(alg_fasta)
    path2raw = path_out+'/total_raw_fasta.faa'
    raw_fasta = open(path2raw, 'w')
    for num, (name, aa, _) in enumerate(fasta):
        clean_aa = aa.replace('-','')
        raw_fasta.write('>'+name+'\n'+clean_aa+'\n')
    raw_fasta.close()

    subprocess.run(("python /home/plaza/soft/eggnog-mapper/emapper.py  --hmm_maxhits 0 --cut_ga --hmm_maxseqlen 60000 --pfam_realign denovo --clean_overlaps clans --dbtype hmmdb --qtype seq \
        -i %s -o %s --output_dir %s" %(path2raw, 'result_emapper', path_out)), shell = True)

    path2pfam_table = path_out+'/result_emapper.emapper.pfam'
    path2main_table = path_out+'/result_emapper.emapper.annotations'
    return path2pfam_table, path2main_table


def annot_tree_pfam_table(post_tree, pfam_table, alg_fasta):

    """
        Annotate tree with pfam table from eggnog-mapper
    """

    fasta = SeqGroup(alg_fasta)
    raw2alg = defaultdict(dict)
    for num, (name, seq, _) in enumerate(fasta):

        p_raw = 1
        for p_alg, (a) in enumerate(seq, 1):
            if a != '-':
                raw2alg[name][p_raw] = p_alg
                p_raw +=1

    seq2doms = defaultdict(list)
    with open(pfam_table) as f_in:
        for line in f_in:
            if not line.startswith('#'):
                info = line.strip().split('\t')
                seq_name = info[0]
                dom_name = info[1]
                dom_start = int(info[7])
                dom_end = int(info[8])

                trans_dom_start = raw2alg[seq_name][dom_start]
                trans_dom_end = raw2alg[seq_name][dom_end]

                dom_info_string = '@'.join([dom_name, str(trans_dom_start), str(trans_dom_end)])
                seq2doms[seq_name].append(dom_info_string)



    for l in post_tree:
        if l.name in seq2doms.keys():
            domains = seq2doms[l.name]
            domains_string = '|'.join(domains)
            l.add_prop('dom_arq', domains_string)

    for n in post_tree.traverse():
        if not n.is_leaf:
            random_seq_name = random.choice(n.leaf_names())
            random_node = post_tree.search_nodes(name=random_seq_name)[0]
            random_node_domains = random_node.props.get('dom_arq', 'none@none@none')
            n.add_prop('dom_arq', random_node_domains)

    return post_tree


def annot_tree_main_table(post_tree, main_table):

    """
        Annotate tree with main table from eggnog-mapper
    """

    seq2info = defaultdict(dict)
    with open(main_table) as fin:
        for line in fin:
            if not line.startswith('#'):
                info = line.strip().split('\t')
                seq_name = info[0]
                eggnog_ogs = info[4]
                for og in eggnog_ogs.split(','):
                    level = og.split('|')[0].split('@')[1]
                    if level in ['2759', '2', '2157'] :
                        basal_og = og.split('|')[0].split('@')[0]
                pref_name = info[8]
                kegg_pathway = info[12]

                seq2info[seq_name]['basal_og'] = basal_og
                seq2info[seq_name]['pref_name'] = pref_name
                seq2info[seq_name]['kegg_path'] = kegg_pathway


    for l in post_tree:
        if l.name in seq2info.keys():
            info_dict = seq2info[l.name]
            for i, val in info_dict.items():
                l.add_prop(i, val)

    for n in post_tree.traverse():
        if not n.is_leaf:
            random_seq_name = random.choice(n.leaf_names())
            random_node = post_tree.search_nodes(name=random_seq_name)[0]
            random_node_basal_og = random_node.props.get('basal_og', 'none@none@none')
            random_node_pref_name = random_node.props.get('pref_name', 'none@none@none')
            random_node_kegg_path = random_node.props.get('kegg_path', 'none@none@none')

            n.add_prop('basal_og', random_node_basal_og)
            n.add_prop('pref_name', random_node_pref_name)
            n.add_prop('kegg_path', random_node_kegg_path)

    return post_tree



#####   FUNCIONTS TO PREPARE OUTPUTS FILES (newick, etc)    ####

def run_clean_properties(t):

    """
        Clean problematic characters from some properties
        Call clean_string()
    """

    # clean strings
    all_props = set()

    #if tree came from web server is  str format,
    if isinstance(t, str):
        t = PhyloTree(t)

    for n in t.traverse():
        for string in ("sci_name", "lca_node_name", "common_name"):
            prop = n.props.get(string)
            if prop:
                n.props[string] = clean_string(prop)

        all_props.update(set(n.props.keys()))


    return t, all_props


def run_write_post_tree(t, name_tree, path_out, all_props):

    """
        Write newick file after the analysis

        TODO: return nothing????
    """

    print('-Post tree')

    name_fam = name_tree.split('.')[0]
    # if os.path.basename(args.tree).split('.')[-1] == 'pickle':
        # with open(path_out+'/post_'+name_fam+'.pickle', "wb") as handle:
            # pickle.dump(t, handle)

    post_tree = path_out+'/'+'post_'+name_fam+'.nw'
    t.write(outfile=post_tree, props=all_props, format_root_node = True)

    return


def get_seq2og(t, best_match=dict()):

    """
        For each seq get all the OG in which are included

        TODO: simplificar???
    """

    seq2ogs = defaultdict(dict)
    nodes_total_leaves = t.leaves()

    for l in nodes_total_leaves:
        ogs_up = list()
        if l.name in best_match.keys():

            for og_recovery, score in best_match[l.name].items():
                og_recovery_name = og_recovery.split('_',1)[0]
                n_recovery = next(t.search_nodes(name=og_recovery_name))#[0]
                ogs_up, dups_up = check_nodes_up(n_recovery)
                ogs_up.add(og_recovery_name)

        else:
            ogs_up, dups_up = check_nodes_up(l)


        for n_up_name in ogs_up:
            n_up = next(t.search_nodes(name=n_up_name))#[0]
            taxlev = n_up.props.get('lca_dup')
            seq2ogs[l.name][taxlev] = n_up_name

        if t.props.get('node_is_og') == 'True':
            lca_root = t.props.get('lca_node')
            root_name = t.props.get('name')
            seq2ogs[l.name][lca_root] = root_name

    if len(seq2ogs) == 0:
        for l in nodes_total_leaves:
            lca_root = t.props.get('lca_node')
            root_name = t.props.get('name')
            seq2ogs[l.name][lca_root] = root_name


    return seq2ogs


def write_seq2ogs(seq2ogs, path, name_tree):

    """
        Write a table with seq2ogs info
    """

    name_fam = name_tree.split('.')[0]
    seq2ogs_out = open(path+'/'+name_fam+'_seq2ogs.tsv', 'w')
    for seq, ogs in seq2ogs.items():
        ogs_out = list()
        for taxlev, og_name in ogs.items():
            ogs_out.append(og_name+'|'+str(taxlev))

        seq2ogs_out.write(seq+'\t'+'@'.join(ogs_out)+'\n')

    seq2ogs_out.close()


def write_ogs_info(ogs_info, pipeline, name_tree,path):

    """
        Write a table with all the info for each OG
    """

    #name_fam = os.path.basename(args.tree).split('.')[0]
    name_fam = name_tree
    if pipeline == 'recovery':
        name_out =  path+'/'+name_fam+'.recovery_ogs_info.tsv'

    elif pipeline == 'original':
        name_out = path+'/'+name_fam+'.original_ogs_info.tsv'

    with open(name_out, "w",  newline='') as myfile:
        w = csv.writer(myfile, delimiter='\t')
        w.writerow(('#OG_name', 'Tax_scope_OG','Dup_name','Dup_lca',  'num_OG_mems', 'Recovery_seqs','OG_up', 'OG_down', 'total_leaves', 'sp_in_OG', 'sp_out_OG', 'members'))


        for og_name, info in ogs_info.items():
            if not 'is_root' in info.keys():

                if isinstance(info['ogs_mems'], str):
                    members_str = info['ogs_mems']
                    num_og_mems = len(members_str.split('|'))

                elif isinstance(info['ogs_mems'], set):
                    members_str = '|'.join(list(info['ogs_mems']))
                    num_og_mems = len(info['ogs_mems'])

                w.writerow((
                    og_name,
                    info['tax_scope_og'],
                    info['name_dup'],
                    info['lca_dup'],
                    num_og_mems,
                    len(info['recovery_mems']),
                    info['ogs_up'],
                    info['ogs_down'],
                    info['total_leaves'],
                    info['num_sp_OGs'],
                    '|'.join(info['num_sp_out']),
                    members_str

                ))
        # syield data.getvalue()
        # sdata.seek(0)
        # sdata.truncate(0)

    # for og_name, info in ogs_info.items():
        # print (og_name, info)


def write_best_match(best_match, path):

    """
        Write a table with the best match result from hmmscan for each seq
    """

    with open(path+'/best_match.tsv', 'w') as fout:
        for seq, best_og in best_match.items():
            for best_og_name in best_og.keys():
                fout.write('\t'.join([seq, best_og_name])+'\n')






###############################################
###############################################
###############################################


def run_app(tree, alg ,name_tree, outliers_node, outliers_reftree, sp_loss_perc, so_cell_org,so_arq, so_bact, so_euk, reftree, user_counter, user_taxo, taxonomy_type, rooting, path_out, inherit_out, mode, run_emapper):

    """
    MAIN FUNCTION, run all steps of the analysis:
        1. Load files needed for the analysis: Tree, taxonomy, reference tree (species tree), taxonomy counter
        2. Pre-analysis: Run some basic analysis:
            rooting, taxonomical annotation, save original species and original members
        3. Outliers and Dups score: Detect long branches, taxonomical outliers, calculate species overlap, score1, score2 and inpalalogs_rate
        4. Detect Duplications and Core-OGs: Select high quality duplication nodes and Core Orthologs Groups (OGs)
        5. Add info about Core-OGs up and down: For each Core-OG detect upper and below Core-OGs (Core-OGs have nested structure)
        6. Write a reference table with the current results for the Core-OGs (for debugging)
        7. Optionally modify the Core-OGs by recovering sequences (see RECOVERY PIPELINE)
        8. Annotate root: Root is a special node that has to be annotated diffent
        9. Flag seqs out OGs: If some seqs still out at least one Core-OGs, add a flag for visualization
        10. Optionally add annotations from emapper (see EMAPPER ANNOTATION)
        11. Write output files:
            seq2ogs
            tree

    RECOVERY PIPELINE
        - Recovery Pipeline: If user provide a alignment file AND some sequences (members) are not included in at least one Core-OGs
                1. Write a fasta file with all raw seqs left out of the Core-OGs and for each Core-OGs write a fasta file
                    1.1 Regular mode: Re-aling sequences
                    1.2 Fast mode: use aligned sequnces
                2. Create HMM file for each Core-OGs fasta file and Build HMM DB
                3. Run Hmmscan to assign seqs out Core-OGs to Core-OGs
                4. For each seq out Core-OGs, select best Core-OGs
                5. Expand_hmm: Create a dict with all the OGs that have recover some seq           TODO: change name expand_hmm
                6. update_og_info: Update og_info with all the seqs taht have been recovered
                7. update_tree: Update tree with all the seqs taht have been recovered
                8. write_ogs_info for recovery
                9. write_best_match
                10. Tar the folder with all the fasta file, HMMs, etc

    EMAPPER ANNOTATION
        - Eggnog-Mapper: annotate tree sequences with eggnog-mapper, user needs to provide alignment file

        TODO
        -   Load OGD_props
        -   Functions with multiple arguments:
                - run_app
                - outliers_detection
                - annotate_root
                - annotate_dups_ch
                - run_preanalysis
                - run_outliers_dup_score
                - run_dups_and_ogs
                - get_taxlevel2ogs (web app)
                - run_write_fastas
    """

    # 1. Load files and DBs
    t = load_tree_local(tree = tree)
    taxonomy_db = load_taxonomy(taxonomy = taxonomy_type, user_taxonomy= user_taxo)
    reftree = load_reftree(rtree = reftree, t = t, taxonomy_db = taxonomy_db)
    taxonomy_counter = load_taxonomy_counter(reftree=reftree, user_taxonomy_counter = user_counter)

    # 2. Pre-analysis: rooting, annotate tree, etc
    t_nw , sp_set, total_mems_in_tree, SPTOTAL, user_props = run_preanalysis(t, name_tree, taxonomy_db, rooting)

    # 3. Outliers and Dups score functions
    t = run_outliers_dup_score(t_nw, outliers_node, outliers_reftree, sp_loss_perc, so_cell_org,so_arq, so_bact, so_euk, taxonomy_counter, taxonomy_db, SPTOTAL, reftree, inherit_out)

    # 4. Detect duplications and OGs
    t, total_mems_in_ogs, ogs_info, taxid_dups_og  = run_dups_and_ogs(t, outliers_node, outliers_reftree, sp_loss_perc, so_cell_org,so_arq, so_bact, so_euk, taxonomy_db, total_mems_in_tree)

    # 5. Add info about Core-OGs up and down
    t, ogs_info = add_ogs_up_down(t, ogs_info)

    # 6. Write a table with the results for the Core-OGs
    write_ogs_info(ogs_info, 'original', name_tree, path_out)

    # 7. Optionally modify the Core-OGs by recovering sequences
    if alg and len(total_mems_in_tree.difference(total_mems_in_ogs)) > 0 and len(total_mems_in_ogs) > 0:
        recovery_seqs, best_match = recover_sequences(tree, t, alg, total_mems_in_tree, total_mems_in_ogs, name_tree, path_out, ogs_info, mode)
    else:
        recovery_seqs = set()
        best_match = defaultdict()

    # 8. Annotate root
    annotate_root(ogs_info, t, total_mems_in_tree, sp_set, total_mems_in_ogs, recovery_seqs, outliers_node, outliers_reftree, sp_loss_perc, so_cell_org, so_arq, so_bact, so_euk, inherit_out, taxonomy_db)

    # 9. Flag seqs out OGs
    t = flag_seqs_out_og(t, total_mems_in_ogs, total_mems_in_tree)

    # 10. Optionally add annotations from emapper
    if run_emapper:
        t = annotate_with_emapper(t, alg, path_out)

    # 11. Write output files
    seq2ogs = get_seq2og(t, best_match)
    write_seq2ogs(seq2ogs, path_out,  name_tree)

    t, all_props = run_clean_properties(t)
    run_write_post_tree(t, name_tree, path_out, all_props)


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', dest = 'tree', required = True)
    parser.add_argument('--raw_alg', dest = 'alg')
    parser.add_argument('--output_path', dest = 'out_path', required= True)
    parser.add_argument('--mode', dest='mode', choices= ["fast", "regular"], default = "regular", type = str)
    parser.add_argument('--species_overlap_cell_org', dest = 'so_cell_org' , default = 0.2, type = float)
    parser.add_argument('--species_overlap_euk', dest = 'so_euk' , default = 0.2, type = float)
    parser.add_argument('--species_overlap_bact', dest = 'so_bact' , default = 0.2, type = float)
    parser.add_argument('--species_overlap_arq', dest = 'so_arq' , default = 0.2, type = float)
    parser.add_argument('--outliers_in_node', dest = 'outliers_node' , default = 0.01, type = float)
    parser.add_argument('--outliers_in_reftree', dest = 'outliers_reftree' , default = 0.05, type = float)
    parser.add_argument('--inherit_outliers', dest = 'inherit_out', choices = ['Yes', 'No'], default = 'Yes', type = str)
    parser.add_argument('--species_losses_perct', dest = 'sp_loss_perc' , default = 0.8, type = float)
    parser.add_argument('--rooting', choices = ['Midpoint', 'MinVar'])
    parser.add_argument('--taxonomy_type', choices=['NCBI', 'GTDB'], default='NCBI')
    parser.add_argument('--user_taxonomy')
    parser.add_argument('--user_taxonomy_counter')
    parser.add_argument('--reftree')
    parser.add_argument('--run_emapper', action='store_true')

    return parser.parse_args()


def main():

    """
        Parse argument and call run_app
    """

    print('\n'+'READ TREE, SEQ FILE, REF TREE')

    args = get_args()

    print(args)

    init_tree = args.tree
    name_tree = os.path.basename(init_tree)
    so_cell_org = args.so_cell_org
    so_euk = args.so_euk
    so_bact = args.so_bact
    so_arq = args.so_arq
    outliers_node = args.outliers_node
    outliers_reftree = args.outliers_reftree
    sp_loss_perc = args.sp_loss_perc
    taxonomy_type = args.taxonomy_type
    rooting = args.rooting
    path_out = args.out_path
    inherit_out = args.inherit_out
    alg = args.alg
    mode = args.mode

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

    run_app(init_tree, alg, name_tree, outliers_node, outliers_reftree, sp_loss_perc, so_cell_org,so_arq, so_bact, so_euk, rtree, user_counter, user_taxo, taxonomy_type, rooting, path_out, inherit_out, mode, args.run_emapper)

    _t.stop()
    print(Timer.timers)



if __name__ == "__main__":
    main()


