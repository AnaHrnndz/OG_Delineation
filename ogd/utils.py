import re

def parse_taxid(node):

    return node.name.split('.')[0]

def parse_taxid_gtdb(node):

    #TODO: add argument for split gen name
    return node.name.split('@')[0]


def get_gtdb_rank(gtdb_code):
    init_rank = gtdb_code.split('__')[0]
    if init_rank == 'r_root':
        rank = 'r_root'
    elif init_rank == 'd':
        rank = 'Domain'
    elif init_rank == 'p':
        rank = 'Phylum'
    elif init_rank == 'c':
        rank = 'Class'
    elif init_rank == 'o':
        rank = 'Order'
    elif init_rank == 'f':
        rank = 'Family'
    elif init_rank == 'g':
        rank = 'Genus'
    elif init_rank == 's':
        rank = 'Species'
    else:
        rank = 'Unk'

    return rank


def get_so2use(taxonomy_db, lin_lca, args):
    print(lin_lca)
    if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
        if 2 in lin_lca:
            so_2_use = args.so_bact
        elif 2759 in lin_lca:
            so_2_use = args.so_euk
        elif 2157 in lin_lca:
            so_2_use = args.so_arq
        elif 'Empty' in lin_lca:
            so_2_use = 0.0
        

    elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
        
        if 'd__Bacteria' in lin_lca:
            so_2_use = args.so_bact
        elif 'd__Archaea' in lin_lca:
            so_2_use = args.so_arq
        elif 'root' in lin_lca:
            so_2_use = args.so_cell_org
        elif 'Empty' in lin_lca:
            so_2_use = 0.0

    else:
        so_2_use = 0.1

    return so_2_use


def get_newick(t, all_props):

    """
        Return tree in newick format with annotations
    """

    t = t.write(props=all_props, format_root_node=True)
    return t

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


def clean_string(string):


    """
        Remove problematic characters for newick format
    """

    clean_string = re.sub(r"'|\[|\]|\=|\.|\-|\:", "", string)
    return clean_string


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