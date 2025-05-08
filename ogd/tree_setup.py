import subprocess
from ete4 import  PhyloTree, GTDBTaxa
import ogd.utils as utils
import re

## 2. Preanalysis - Tree setup  ##

chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'

def run_setup(t, name_tree, taxonomy_db, path_out, tmp_path, args):


    """
        Preanalysis include several steps:
            resolve polytomies: need for sp overlap
            rooting: Midpoint o MinVar
            add taxonomy annotation
            Get original seqs and species in the tree
            Name internal nodes
    """

    rooting = args.rooting

    sp_delimitator = args.sp_delim

    t.resolve_polytomy()

    t = check_branch_legth(t)

    t = run_rooting(t, rooting, tmp_path, sp_delimitator)

    t = add_taxomical_annotation(t, taxonomy_db)

    # Total members(leafs name) in tree
    total_mems_in_tree = set(t.leaf_names())

    #Create an ID for each internal node
    for i, n in enumerate(t.traverse()):
        if not n.is_leaf:   
            n.name = '%s-%d' % (make_name(i), get_depth(n))

    set_sp_total = t.get_species()
    num_total_sp = len(set_sp_total)

    mssg = f"""
    1. Pre-analysis
        -Rooting: {rooting}
        -Len tree: {len(t)}
        -Total species in tree: {num_total_sp}"""
    print(mssg)
    

    # Newick format need for web
    t, props = utils.run_clean_properties(t)
    tree_nw = utils.get_newick(t, props)

    return tree_nw, set_sp_total, total_mems_in_tree, num_total_sp



def check_branch_legth(t):

    """
    Preguntar Jordi
    """
    lenghts = list()
    for n in t.traverse():
        if n.dist != None:
            lenghts.append(n.dist)
        

    if (all(v == 0.0 for v in lenghts)) == True:
        for n in t.traverse():
            if n.dist != None:
                n.dist = 1.0
            else:
                print(len(n))

    return t


def run_rooting(t, rooting, tmp_path, sp_delimitator):

    """
        Tree rooting.
        Midpoint or MinVar methods availables
    """

    if  rooting == "Midpoint":
        root_mid = t.get_midpoint_outgroup()
        try:
            t.set_outgroup(root_mid)
        except:
            print('Error in Midpoint')

    elif rooting == "MinVar":
        t = run_minvar(t, tmp_path, sp_delimitator)

    else:
        print('No rooting')

    return t


def run_minvar(t, tmp_path, sp_delimitator):

    """
        With MinVar rooting, you need to :
        1. Write the tree, polytomies had been resolved in previous step
        2. Run MinVar
        3. Open it again with PhyloTree
    
    """
   
    
    input_tree_minvar = utils.write_tree_for_minvar_rootin(t, tmp_path)
    path2tree = input_tree_minvar
    path2tmptree = tmp_path+'output_minvar_tree.nw'
    stdout_file = open(tmp_path+'minvar.stdout', 'w')
    stderr_file = open(tmp_path+'minvar.stderr', 'w')

    subprocess.run(("python3  $(which FastRoot.py) -i %s -o %s" \
        %(path2tree, path2tmptree)), shell = True, stdout = stdout_file, stderr= stderr_file)
    
    stdout_file.close()
    stderr_file.close()

    t_minvar = PhyloTree(open(path2tmptree), parser = 0)

    #check that minvar run and its not oppening an old tree from previous run

    t_minvar.set_species_naming_function(lambda node: node.name.split(sp_delimitator)[0])

    return t_minvar


def add_taxomical_annotation(t, taxonomy_db):

    """
        Add taxonomical annotation to nodes (sci_name, taxid, named_lineage, lineage, rank)
        Parsing function used to extract species name from a node’s name.
    """

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
