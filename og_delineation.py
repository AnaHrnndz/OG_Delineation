#!/usr/bin/env python3

from ete4 import  PhyloTree, SeqGroup, Tree
from ete4 import NCBITaxa, GTDBTaxa
from collections import Counter, OrderedDict, defaultdict
import sys
import os
import numpy as np
import json
#import string
import argparse
import warnings
import re
import glob
import subprocess
import csv
import pathlib


sys.path.append('/data/projects/og_delineation/ogd')
from emapper_annotate import annotate_with_emapper
from recovery import recover_sequences
import  pairs
from timer import Timer

cwd =  str(pathlib.Path(__file__).parent.resolve())

bin_path =  cwd+'/bin'
data_path = cwd+'/data'
sys.path.append(bin_path)

warnings.filterwarnings("ignore", category=RuntimeWarning)

chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'

# OGD_props = list()
# with open(data_path+'/OGD_props.txt') as f:
    # for line in f:
        # OGD_props.append(line.strip())


class OrderedCounter(Counter, OrderedDict):
     'Counter that remembers the order elements are first seen'
     def __repr__(self):
         return '%s(%r)' % (self.__class__.__name__,
                            OrderedDict(self))
     def __reduce__(self):
         return self.__class__, (OrderedDict(self),)


####    FUNCTIONS FOR OG DELINEATION WEB    ####

def clean_folder(dir, aln_path):

    """
        Remove files in tmp dir created for web analysis
    """

    filelist = glob.glob(os.path.join(dir, "*"))
    for f in filelist:

        if f != aln_path:
            os.remove(f)

def load_annotated_tree(t, name_tree):

    sp_set = set(t.species)
    total_mems_in_tree = set(t.leaf_names())
    NUM_TOTAL_SP = len(sp_set)


    return t, sp_set, total_mems_in_tree, NUM_TOTAL_SP

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

def prune_tree(t, total_mems_in_ogs):

    """
        Prune tree to keep only mems that belong to at least one OG
        We need this version of the tree to visualization
    """

    t.prune(list(total_mems_in_ogs))
    return t




####   MAIN FUNCTIONS  ####

def parse_taxid(node):

    #TODO: add argument for split gen name
    return node.name.split('.')[0]

def clean_string(string):

    """
        Remove problematic characters for newick format
    """

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
            else:
                ogs_up.add(node.up.props.get('name'))
        node = node.up

    return ogs_up, dups_up



## 1. Load initial info   ##

def load_tree_local(tree=None):

    """
        Load Tree from file
    """

    print(' -Load tree:', os.path.basename(tree))

    t = PhyloTree(open(tree), parser = 0)

    t.set_species_naming_function(parse_taxid)

    return t


def load_tree_web(tree=None):

    """
        Load tree from newick format
    """

    print('-Load tree: ')

    t = PhyloTree((tree))

    t.set_species_naming_function(parse_taxid)

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

    print(' -Load taxonomy:  ')

    return taxonomy_db


def load_reftree(rtree=None, t=None, taxonomy_db=None):

    """
        Get reference tree (species tree) from input tree or user provide it
    """

    print(' -Load reftree:   ')
    if rtree:
            reftree = PhyloTree(open(rtree))
    else:
        reftree = get_reftree(t, taxonomy_db)


    taxonomy_db.annotate_tree(reftree,  taxid_attr="name")

    return reftree

def get_reftree(t, taxonomy_db):

    """
        Create reference tree (species tree) if user do not provide it
    """

    print('\t**Create referente species tree from gene tree')

    taxid_list = t.get_species()

    reftree = taxonomy_db.get_topology(taxid_list)

    return reftree


def load_taxonomy_counter(reftree=None, user_taxonomy_counter=None):

    """
        Get number of species per each taxonomic level
        Get it from reference tree or user provide it
    """

    print(' -Load taxonomy counter:  ')
    if user_taxonomy_counter:
        if isinstance(user_taxonomy_counter, dict):
            level2sp_mem = user_taxonomy_counter
        else:
            with open(user_taxonomy_counter) as levels:
                level2sp_mem = json.load(levels)
    else:
        level2sp_mem = get_taxonomy_counter(reftree)

    return level2sp_mem

def get_taxonomy_counter(reftree):

    """
        Create Taxonomy counter if user de not provide it
    """

    print('\t**Create taxonomy counter  from gene tree')

    level2sp_mem = defaultdict(set)
    for l in reftree:
        lin = l.props.get('lineage')
        for tax in lin:
            level2sp_mem[str(tax)].add(l.name)

    return level2sp_mem


## 2. Preanalysis  ##

def run_preanalysis(t, name_tree, taxonomy_db, rooting, path_out, abs_path):


    """
        Preanalysis include several steps:
            rooting
            add taxonomy annotation
            Get original seqs and species in the tree
            Name internal nodes
    """

    print('\n'+'1. Pre-analysis:')

    t = run_rooting(t, rooting, path_out, abs_path)

    t = add_taxomical_annotation(t, taxonomy_db)

    # Total members(leafs name) in tree
    total_mems_in_tree = set(t.leaf_names())

    for i, n in enumerate(t.traverse()):
        if not n.is_leaf:
            #Create an ID for each internal node
            n.name = '%s-%d' % (make_name(i), get_depth(n))


    set_sp_total = t.get_species()
    NUM_TOTAL_SP = len(set_sp_total)

    print(' -Len tree: ', len(t))
    print(' -Total species in tree:', NUM_TOTAL_SP)

    t, props = run_clean_properties(t)
    tree_nw = get_newick(t, props)

    return tree_nw, set_sp_total, total_mems_in_tree, NUM_TOTAL_SP, props

def run_rooting(t, rooting,  path_out, abs_path):

    """
        Tree rooting.
        TODO: Add other methods to root the tree as MinVar
    """

    if  rooting == "Midpoint":
        print(' -Run midpoint rooting')
        root_mid = t.get_midpoint_outgroup()
        t.set_outgroup(root_mid)

    elif rooting == "MinVar":
        print(' -Run MinVar rooting')
        t = run_minvar(t,path_out, abs_path)

    else:
        print(' -No rooting')

    #¿Change?
    t.dist = 0.01

    return t

def run_minvar(t, path_out, abs_path):

    path2tree = abs_path
    path2tmptree = path_out+'/tmp_tree.nw'

    subprocess.run(("python /data/soft/FastRoot/MinVar-Rooting/FastRoot.py -i %s -o %s" \
        %(path2tree, path2tmptree)), shell = True)

    t = load_tree_local(tree = path2tmptree)

    return t


def add_taxomical_annotation(t, taxonomy_db):

    """
        Add taxonomical annotation to nodes
        Parsing function used to extract species name from a node’s name.
    """

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
def run_outliers_and_scores(t_nw, taxonomy_db, NUM_TOTAL_SP, level2sp_mem, args):

    print('\n'+'2. Detect outlierts and get scores')
    print(' -Outliers thresholds:')
    print('\tNode: '+ str(args.outliers_node))
    print('\tReftree: '+ str(args.outliers_reftree))
    print(' -Species losses percentage threshold: ' + str(args.sp_loss_perc))

    t = PhyloTree(t_nw)

    CONTENT = t.get_cached_content()

    t.add_prop('is_root', str('True'))

    long_leaves = detect_long_branches(t)
    if args.user_taxonomy_counter:
        root_outliers = root_taxono_outliers(t, NUM_TOTAL_SP, level2sp_mem)
    else:
        root_outliers = set()

    sp_in_root = t.species


    t.add_prop('_sp_in', sp_in_root)

    for n in t.traverse("preorder"):

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
                sp_out.update(root_outliers)
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

            if 2 in lin_lca:
                so_2_use = args.so_bact
            elif 2759 in lin_lca:
                so_2_use = args.so_euk
            elif 2157 in lin_lca:
                so_2_use = args.so_arq
            elif lca_node == 131567:
                so_2_use = args.so_cell_org
            else:
                so_2_use = 0.1

            if so_score >= so_2_use:
                if float(n.props.get('species_losses_percentage')) > args.sp_loss_perc:
                    n.add_prop('evoltype_2', 'FD')
                else:
                    n.add_prop('evoltype_2', 'D')
            else:
                n.add_prop('evoltype_2', 'S')

    t, props = run_clean_properties(t)

    return t, CONTENT

def detect_long_branches(t):

    """
        Detect leaves with long branches
        Long branches leaves are when a leave's branch is 10 times bigger than the mean of the branches for the whole tree
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

def root_taxono_outliers(t, NUM_TOTAL_SP, level2sp_mem):
    root_outliers = set()

    for taxid in ['2', '2759', '2157']:

        if taxid in level2sp_mem.keys():
            if (len(level2sp_mem[taxid])/ NUM_TOTAL_SP) < 0.1:
                root_outliers.update(level2sp_mem[taxid])

    return root_outliers

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

        - outliers_node: default 0.01
        - outliers_reftree: default 0.05

        The algorithm traverse the taxonomy counter created for the node.
        Only get outliers two taxonomical levels (depht) below node lca
            -> If node is metazoa, outlierts will be evaluated in Eumetazoa, Porifera... and Bilateria, Cnidaria etc
            -> Outliers wont be evaluated for level such as Mammalia, Gnathostomata, etc
                    Example: cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata

            Calculate:
                - per_node: Percentage of lineage in the node. How important is the lineage in comparison with all species that are present in the node
                - per_tree: Percentage of lineage in reference species tree. How well represented is the lineage in the node in comparison with the same lineage in the reftree

            If a lineage in the node is below 1% and below 5% in the reftree, remove that species, that means:
                1. The the species in that lineage are rare in the node (below 1%)
                    Few porifera species inside a node with all bilateria species
                2. There are few species representatives from reftree at that taxonomic level (below 5%)
                    There are 2 porifera species but in reftree there are 2000 porifera species -> 2/200=0.01
                    Taxonomic level that are rare in reftree (just 1 or 2 species of one lineage), will be preserved.
                        If in the reftree you have two cnidaria, you always keeps cnidaria species -> 1/2=0.5

    """


    sp_in_node = set()

    # Create taxonomy counter for the node. For each lineage that its present in the node, count how many species there are in the node
    sp_per_level = defaultdict(set)

    for l in CONTENT[n]:
        if l.props.get('taxid') not in sp_out_up:
            sp_in_node.add(l.props.get('taxid'))

            for tax in l.props.get('lineage').split('|'):
                sp_per_level[tax].add(str(l.props.get('taxid')))


    best_tax = str()
    ptax = defaultdict(dict)

    for tax, num in sp_per_level.items():
        depth_tax = len(taxonomy_db.get_lineage(tax))
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

    sp_keep = sp_per_level[best_tax]
    candidates2remove = sp_in_node.difference(sp_keep)
    sp2remove = set()

    for sp in candidates2remove:
        for l in taxonomy_db.get_lineage(sp):

            # Representacion de esa especie para cada taxid de su linaje en el arbol, si es una sp que tiene pocos representantes en el arbol sp_in_tree~1.
            sp_in_tree = 1/len(level2sp_mem[(str(l))])

            # Representacion del linaje en el nodo, en comparacion con el linaje en el arbol. Si ese linaje contiene casi todos los representantes del linaje en el arbol, lin_in_tree~1
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
        for tax in taxonomy_db.get_lineage(sp):
            sp_lost_per_level[tax].add(sp)


    sp_exp_per_level = defaultdict(set)
    for sp in expected_sp:
        for tax in taxonomy_db.get_lineage(sp):
            sp_exp_per_level[tax].add(sp)

    ptax = defaultdict(dict)
    for tax, num in sp_lost_per_level.items():
        depth_tax = len(taxonomy_db.get_lineage(tax))
        perc_tax_in_node = len(num)/len(sp_exp_per_level[tax]) # ????????
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

    diff_sp = expected_sp.difference(sp_in)
    perc_diff_sp = len(diff_sp) / len(expected_sp)

    n.add_prop('species_losses', len(diff_sp))
    n.add_prop('species_losses_percentage', perc_diff_sp)




##  4. Detect Duplications and Core-OGs & PGs ##

def run_get_main_dups(t, taxonomy_db, total_mems_in_tree, args):
    """
        Select high-quality duplication nodes:
            - Species overlap higher than min threshold
            - More than one leave and more than one specie
            - No more duplication nodes at the same taxonomic level below the node
    """



    #taxid_dups_og, set with the lca of the nodes that are OGs
    taxid_dups_og = set()
    total_mems_in_ogs = set()

    #Traverse tree to find the nodes that are "good" duplications and generate OGs.
    print('\n'+'3. Select high quality duplication nodes')
    print(' -Species overlap threshold:')
    print('\t'+'Cell Org: '+ str(args.so_cell_org))
    print('\t'+'Euk: '+ str(args.so_euk))
    print('\t'+'Bact :' + str(args.so_bact))
    print('\t'+'Arq :' + str(args.so_arq))


    for node in t.traverse("preorder"):


        if not node.is_leaf and node.props.get('evoltype_2') == 'D'  \
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

                save_dups_ch1 = 0 #defaultdict(int)
                save_dups_ch2 = 0 #defaultdict(int)


                # Check that dups under child1 and child2 (that have the same lca) fit out requirements : species overlap min requirement,
                # more than 1 leaves and more than 1 species
                #if len(dups_under_ch1)> 0:
                for n_ in  dups_under_ch1:
                    #if float(n_.props.get('so_score'))>= so_2_use \
                    if len(n_.props.get('_leaves_in')) >1 and len(n_.props.get('_sp_in'))> 1 :
                        save_dups_ch1 += 1


                #if len(dups_under_ch2)> 0:
                for n_ in  dups_under_ch2:
                   # if float(n_.props.get('so_score'))>= so_2_use  \
                    if  len(n_.props.get('_leaves_in'))>1 and len(n_.props.get('_sp_in'))> 1 :
                        save_dups_ch2 += 1

                #If dups under child1 do not achieve our requirement, then child1 is OG
                if save_dups_ch1 == 0:
                    annotate_dups_ch(total_mems_in_ogs, taxid_dups_og,node, 'ch1', taxonomy_db)

                #If dups under child2 do not achieve our requirement, then child2 is OG
                if save_dups_ch2 == 0 :
                    annotate_dups_ch(total_mems_in_ogs, taxid_dups_og,node, 'ch2', taxonomy_db)


            # No more dups under node with the same lca_node, both ch nodes are OGs
            elif len(dups_under_node) == 0:
                annotate_dups_ch(total_mems_in_ogs, taxid_dups_og,node, 'ch1', taxonomy_db)
                annotate_dups_ch(total_mems_in_ogs, taxid_dups_og,node, 'ch2', taxonomy_db)


    lca_root = int(t.props.get('lca_node'))
    if len(list(t.search_nodes(evoltype_2='D', lca_node=lca_root))) == 0:
        t.add_prop('node_is_og', 'True')
        total_mems_in_ogs.update(set(t.props.get('_leaves_in')))

    t, props = run_clean_properties(t)

    return t, total_mems_in_ogs, taxid_dups_og

def annotate_dups_ch(total_mems_in_ogs, taxid_dups_og, node, ch_node, taxonomy_db):

    """
    Add props and save info about the children node that is OG
    node is the duplication node
    target node is the node that is OG
    """

    if ch_node == 'ch1':
        og_name_ch = node.props.get('_ch1_name')
        og_ch_mems = node.props.get('_leaves_ch1')
        sp_ch = node.props.get('_sp_in_ch1')
        target_node = next(node.search_nodes(name=og_name_ch))

    elif ch_node == 'ch2':
        og_name_ch = node.props.get('_ch2_name')
        og_ch_mems = node.props.get('_leaves_ch2')
        sp_ch = node.props.get('_sp_in_ch2')
        target_node = next(node.search_nodes(name=og_name_ch))

    if  len(sp_ch) > 1 and len(og_ch_mems) > 1:
        target_node.add_prop('node_is_og', 'True')
        target_node.add_prop('lca_dup', node.props.get('lca_node'))
        target_node.add_prop('so_score_dup', node.props.get('so_score'))
        target_node.add_prop('dup_lineage', list(taxonomy_db.get_lineage(node.props.get('lca_node'))))
        #target_node.add_prop('_mems_og', '|'.join(list((og_ch_mems))))
        target_node.add_prop('_dup_node_name', node.props.get('name'))

        taxid_dups_og.add(node.props.get('lca_node'))
        node.add_prop('node_create_og', 'True')
        sp_in_og = target_node.props.get('_sp_in')
        sp_in_anc_og = node.props.get('_sp_in')

        total_mems_in_ogs.update(set(og_ch_mems))

def run_ogs_v3(t, level2sp_mem, taxonomy_db):

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



##  5. Add info about Core-OGs up and down  ##

def add_nodes_up_down(t):

    """
        Traverse the tree, for each node,  detect ogs that are up or below
    """

    for node in t.traverse('preorder'):
        if node.props.get('node_is_og'):
            node_name = node.name

            #Detect OGs below node
            ogs_down = set()
            dups_down = list()
            for n in node.search_nodes(node_is_og="True"):
                if node_name != n.name:
                    ogs_down.add(n.name)
                    dups_down.append(n.up.name)

            if len(ogs_down) > 0:
                node.add_prop('ogs_down',ogs_down)
                node.add_prop('dups_down', dups_down)
            else:
                node.add_prop('ogs_down','-')
                node.add_prop('dups_down', '-')

            #Detect OGs up
            ogs_up = set()
            dups_up = list()
            ogs_up, dups_up = check_nodes_up(node)

            if len(ogs_up) > 0:
                node.add_prop('ogs_up', ogs_up)
                node.add_prop('dups_up', dups_up)
            else:
                node.add_prop('ogs_up', '-')
                node.add_prop('dups_up', '-')

    ogs_down = set()
    dups_down = list()
    for n in t.search_nodes(node_is_og="True"):
        if n.name != t.name:
            ogs_down.add(n.name)
            dups_down.append(n.up.name)

    if len(ogs_down) > 0:
        t.add_prop('ogs_down',ogs_down)
        t.add_prop('dups_down', dups_down)
    else:
        t.add_prop('ogs_down','-')
        t.add_prop('dups_down', '-')


    return t


def annot_ogs(t, base_ogs, taxonomy_db):

    ##OG_name   TaxoLevel   AssocNode  len_sp_in_OG OG_up   OG_down num_OG_mems    members

    annot_base_ogs = defaultdict(dict)
    for taxa, ogs in base_ogs.items():
        for og_name, og_info in ogs.items():
            sci_name_taxa =  taxonomy_db.get_taxid_translator([taxa])[int(taxa)]
            ogs_down = '|'.join(list((t[og_info[0]].props.get('ogs_down','-'))))
            ogs_up = '|'.join(list((t[og_info[0]].props.get('ogs_up', '-'))))
            recover_seqs = t[og_info[0]].props.get('recovery_seqs','-')
            mems = og_info[1].split('|')
            sp_in_og = set()
            for l in mems:
                sp_in_og.add(l.split('.')[0])


            annot_base_ogs[og_name]['TaxoLevel'] = taxa
            annot_base_ogs[og_name]['SciName_TaxoLevel'] = sci_name_taxa.replace(' ', '_')
            annot_base_ogs[og_name]['AssocNode'] = og_info[0]
            annot_base_ogs[og_name]['NumSP'] = len(sp_in_og)
            annot_base_ogs[og_name]['OG_down'] = ogs_down
            annot_base_ogs[og_name]['OG_up'] = ogs_up
            annot_base_ogs[og_name]['NumMems'] = len(mems)
            annot_base_ogs[og_name]['Mems'] = '|'.join(mems)
            annot_base_ogs[og_name]['RecoverySeqs'] = '|'.join(recover_seqs)
            annot_base_ogs[og_name]['NumRecoverySeqs'] = str(len(recover_seqs))

    return annot_base_ogs



##  8.Annotate root ##
def annotate_root(base_ogs_annot, t, name_tree, total_mems_in_tree, sp_set, total_mems_in_ogs, recovery_seqs, taxonomy_db, args):

    """
    Add properties in root
    Add info about the root in ogs_info dict
    """

    # Save parameters as property

    parameters_str = '|'.join(["outliers_node@"+str(args.outliers_node), "outliers_reftree@"+str(args.outliers_reftree), "sp_loss_perc@"+str(args.sp_loss_perc),
                                "so_cell_org@"+str(args.so_cell_org), "so_arq@"+str(args.so_arq), "so_bact@"+str(args.so_bact), "so_euk@"+str(args.so_euk), "inherit_out@"+str(args.inherit_out)])

    t.add_prop('parameters', parameters_str)

    # Save general results as properties, needed for web

    seqs_out_ogs = total_mems_in_tree.difference(total_mems_in_ogs)

    general_results_str = '|'.join(["Tree_name@"+name_tree, "Total_seqs@"+str(len(total_mems_in_tree)), "Total_species@"+str(len(sp_set)),
                                    "Seqs_in_OGs@"+str(len(total_mems_in_ogs)), "Recovery_seqs@"+str(len(recovery_seqs)), "Seqs_out_OGs@"+str(len(seqs_out_ogs)), "Num_OGs@"+str(len(base_ogs_annot)) ])

    t.add_prop('general_result', general_results_str)


    # Save taxlevel OGs info, needed for web

    lca_all_dups = set()
    taxlev2ogs = defaultdict(set)
    taxlev2mems = defaultdict(set)

    for dup_node in t.search_nodes(node_create_og='True'):
        lca_dup = dup_node.props.get('lca_node')
        mems_dup = set(dup_node.props.get('_leaves_in'))
        lca_all_dups.add(lca_dup)
        taxlev2ogs[lca_dup].add(dup_node.name)
        taxlev2mems[lca_dup].update(mems_dup)

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
    name_fam = name_tree.split('.',1)[0]

    post_tree = path_out+'/'+name_fam+'.tree_annot.nw'
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
    seq2ogs_out = open(path+'/'+name_fam+'.seq2ogs.tsv', 'w')
    for seq, ogs in seq2ogs.items():
        ogs_out = list()
        for taxlev, og_name in ogs.items():
            ogs_out.append(og_name+'|'+str(taxlev))

        seq2ogs_out.write(seq+'\t'+'@'.join(ogs_out)+'\n')

    seq2ogs_out.close()


def write_ogs_info(base_ogs, pipeline, name_tree, path):

    """
        Write a table with all the info for each OG
    """


    name_fam = name_tree
    if pipeline == 'recovery':
        name_out =  path+'/'+name_fam+'.recovery_ogs_info.tsv'

    elif pipeline == 'base':
        name_out = path+'/'+name_fam+'.base_ogs_info.tsv'

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


def write_ogs_info_v2(base_ogs_annot, name_tree, path):

    name_fam = name_tree.split('.',1)[0]
    name_out =  path+'/'+name_fam+'.ogs_info.tsv'

    # if pipeline == 'recovery':
        # name_out =  path+'/'+name_fam+'.recovery_ogs_info.tsv'

    # elif pipeline == 'base':
        # name_out = path+'/'+name_fam+'.base_ogs_info.tsv'

    with open(name_out, "w",  newline='') as myfile:
        w = csv.writer(myfile, delimiter='\t')
        w.writerow(('#OG_name', 'TaxoLevel', 'SciName_TaxoLevel', 'AssocNode',  'NumSP', 'OG_down', 'OG_up', 'NumSeqs', 'NumRecoverySeqs','Seqs', 'RecoverySeqs'))

        for og_name in base_ogs_annot.keys():

            w.writerow((
                og_name,
                base_ogs_annot[og_name]['TaxoLevel'],
                base_ogs_annot[og_name]['SciName_TaxoLevel'],
                base_ogs_annot[og_name]['AssocNode'],
                base_ogs_annot[og_name]['NumSP'],
                base_ogs_annot[og_name]['OG_down'],
                base_ogs_annot[og_name]['OG_up'] ,
                base_ogs_annot[og_name]['NumMems'],
                base_ogs_annot[og_name]['NumRecoverySeqs'],
                base_ogs_annot[og_name]['Mems'],
                base_ogs_annot[og_name]['RecoverySeqs']
            ))





###############################################
###############################################
###############################################



def run_app(tree, abs_path, name_tree, reftree, user_counter, user_taxo, taxonomy_type, path_out, args):
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
                5. Expand_hmm: Create a dict with all the OGs that have recover some seq
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
                - run_outliers_and_scores
                - run_get_main_dups
                - get_taxlevel2ogs (web app)
                - run_write_fastas
    """

    # 1. Load files and DBs
    print('\n0. Load info')
    t = load_tree_local(tree = tree)
    taxonomy_db = load_taxonomy(taxonomy = taxonomy_type, user_taxonomy= user_taxo)
    reftree = load_reftree(rtree = reftree, t = t, taxonomy_db = taxonomy_db)
    level2sp_mem = load_taxonomy_counter(reftree=reftree, user_taxonomy_counter = user_counter)

    # 2. Pre-analysis: rooting, annotate tree, etc
    t_nw , sp_set, total_mems_in_tree, NUM_TOTAL_SP, user_props = run_preanalysis(t, name_tree, taxonomy_db, args.rooting, path_out, abs_path)

    # 3. Outliers and Dups score functions
    t, CONTENT = run_outliers_and_scores(t_nw, taxonomy_db, NUM_TOTAL_SP, level2sp_mem, args)

    # 4. Detect duplications and OGs
    t, total_mems_in_ogs, taxid_dups_og  = run_get_main_dups(t, taxonomy_db, total_mems_in_tree, args)

    # 5. Retrieve OGs
    base_ogs = run_ogs_v3(t, level2sp_mem,taxonomy_db)

    # 6. Add info about Nodes that split OGs up and down
    t = add_nodes_up_down(t)

    # 7. Annotate base-OGs
    base_ogs_annot = annot_ogs(t, base_ogs, taxonomy_db)

    # 8. Write a table with the results for the Core-OGs
    write_ogs_info_v2(base_ogs_annot, name_tree, path_out)


    # 9. Optionally modify the Core-OGs by recovering sequences
    recovery_seqs = set()
    best_match = defaultdict()
    if args.run_recovery:
        if args.alg and len(total_mems_in_tree.difference(total_mems_in_ogs)) > 0 and len(total_mems_in_ogs) > 0:
            recovery_seqs, best_match = recover_sequences(tree, t, args.alg, total_mems_in_tree, total_mems_in_ogs, name_tree, path_out, base_ogs_annot, args.mode)
            recovery_ogs_annot = annot_ogs(t, base_ogs, taxonomy_db)
            write_ogs_info_v2(base_ogs_annot, name_tree, path_out)

    else:
        recovery_seqs = set()
        best_match = defaultdict()

    # 10. Annotate root
    annotate_root(base_ogs_annot, t, name_tree, total_mems_in_tree, sp_set, total_mems_in_ogs, recovery_seqs, taxonomy_db, args)

    # 11. Flag seqs out OGs
    t = flag_seqs_out_og(t, total_mems_in_ogs, total_mems_in_tree)

    # 12. Optionally add annotations from emapper
    if args.run_emapper:
        t = annotate_with_emapper(t, args.alg, path_out)

    if args.get_pairs:
        clean_pairs = pairs_module.get_all_pairs(CONTENT, total_mems_in_ogs)
        pairs_module.write_pairs_table(clean_pairs, path_out, name_tree)

    # 13. Write output files
    print('\n4. Writing outputfiles\n')
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
    parser.add_argument('--species_overlap_cell_org', dest = 'so_cell_org' , default = 0.1, type = float)
    parser.add_argument('--species_overlap_euk', dest = 'so_euk' , default = 0.1, type = float)
    parser.add_argument('--species_overlap_bact', dest = 'so_bact' , default = 0.1, type = float)
    parser.add_argument('--species_overlap_arq', dest = 'so_arq' , default = 0.1, type = float)
    parser.add_argument('--outliers_in_node', dest = 'outliers_node' , default = 0.05, type = float)
    parser.add_argument('--outliers_in_reftree', dest = 'outliers_reftree' , default = 0.01, type = float)
    parser.add_argument('--inherit_outliers', dest = 'inherit_out', choices = ['Yes', 'No'], default = 'Yes', type = str)
    parser.add_argument('--species_losses_perct', dest = 'sp_loss_perc' , default = 0.7, type = float)
    parser.add_argument('--rooting', choices = ['Midpoint', 'MinVar', ''])
    parser.add_argument('--taxonomy_type', choices=['NCBI', 'GTDB'], default='NCBI')
    parser.add_argument('--user_taxonomy')
    parser.add_argument('--user_taxonomy_counter')
    parser.add_argument('--reftree')
    parser.add_argument('--run_emapper', action='store_true')
    parser.add_argument('--run_recovery', action='store_true')
    parser.add_argument('--get_pairs', action='store_true')

    return parser.parse_args()


def main():

    """
        Parse argument and call run_app
    """

    args = get_args()

    init_tree = args.tree
    name_tree = os.path.basename(init_tree)
    abs_path = os.path.abspath(init_tree)
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
    _t.start()

    run_app(init_tree, abs_path, name_tree, rtree, user_counter, user_taxo, taxonomy_type, path_out, args)

    _t.stop()
    print(Timer.timers)



if __name__ == "__main__":
    main()


