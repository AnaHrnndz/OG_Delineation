import subprocess
from ete4 import  PhyloTree
import utils
import re

## 2. Preanalysis - Tree setup  ##



chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'

def run_setup(t, name_tree, taxonomy_db, rooting, path_out, abs_path):


    """
        Preanalysis include several steps:
            resolve polytomies: need for sp overlap
            rooting
            add taxonomy annotation
            Get original seqs and species in the tree
            Name internal nodes
    """

    print('\n'+'1. Pre-analysis:')

    t.resolve_polytomy()

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

    t, props = utils.run_clean_properties(t)
    tree_nw = utils.get_newick(t, props)

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

    t.dist = 0.01

    return t

def run_minvar(t, path_out, abs_path):

    path2tree = abs_path
    path2tmptree = path_out+'/tmp_tree.nw'

    # subprocess.run(("python /data/soft/FastRoot/MinVar-Rooting/FastRoot.py -i %s -o %s" \
        # %(path2tree, path2tmptree)), shell = True)

    subprocess.run(("FastRoot.py -i %s -o %s" \
        %(path2tree, path2tmptree)), shell = True)

    #t = load_tree_local(tree = path2tmptree)

    t_minvar = PhyloTree(open(path2tmptree), parser = 0)
    t_minvar.set_species_naming_function(utils.parse_taxid)


    return t_minvar

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
