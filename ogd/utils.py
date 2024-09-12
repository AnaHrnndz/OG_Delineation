import re
from ete4 import PhyloTree

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


def get_so2use(taxonomy_db, lin_lca, lca_node,args):
    
    if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
        
        if args.so_euk != None:
            if 2759 in lin_lca:
                so_2_use = args.so_euk
            else:
                so_2_use = args.so_all

        elif args.so_bact != None:
            if 2 in lin_lca:
                so_2_use = args.so_bact
            else:
                so_2_use = args.so_all

        elif args.so_arq != None:
            if 2157 in lin_lca:
                so_2_use = args.so_arq
            else:
                so_2_use = args.so_all

        else:
            so_2_use = args.so_all

        # elif args.so_specf == None and  args.so_lin == None:
            # so_2_use = args.so_all

        # elif args.so_lin != None and args.so_specf == None:
            # lev = int(args.so_lin[0])
            # so = args.so_lin[1]
            # if lev in lin_lca:
                # so_2_use = so
            # else:
                # so_2_use = args.so_all

        # elif args.so_specf != None and args.so_lin == None:
            # lev = int(args.so_specf[0])
            # so = args.so_specf[1]
            # if lev == lca_node :
                # so_2_use = so
            # else:
                # so_2_use = args.so_all

        # elif args.so_specf != None and args.so_lin != None:
            # lev_only = int(args.so_specf[0])
            # lev_from = int(args.so_lin[0])
           
            # if lev_only == lca_node :
                # so = args.so_specf[1]
                # so_2_use = so
            # elif lev_from in lin_lca:
                # so = args.so_lin[1]
                # so_2_use = so
            # else:
                # so_2_use = args.so_all

        
        
        #if  in lin_lca:
        #    so_2_use = args.so_bact
        #elif 2759 in lin_lca:
        #    so_2_use = args.so_euk
        #elif 2157 in lin_lca:
        #    so_2_use = args.so_arq
        #elif 131567 in lin_lca:
        #    so_2_use = args.so_cell_org
        #elif 'Empty' in lin_lca:
        #    so_2_use = 0.0
        #
        #else:
        #    so_2_use = args.so_all
        

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
            so_2_use = args.so_all



    
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



def run_write_post_tree(t, name_tree, path_out, all_props):

    """
        Write newick file after the analysis
        TODO: return nothing????
    """
    name_fam = name_tree.split('.',1)[0]

    post_tree = path_out+'/'+name_fam+'.tree_annot.nw'
    t.write(outfile=post_tree, props=all_props, format_root_node = True)

    return