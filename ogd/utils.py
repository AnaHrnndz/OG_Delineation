import re
from ete4 import PhyloTree
from collections import defaultdict

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




def remove_problematic_characters(name):
    clean_name = name.replace('-','_').replace('.faa', '').replace('.nw', '').replace('.fa', '').replace('.fasta','')

    return clean_name



def run_write_post_tree(t, clean_name_tree, path_out, all_props):

    """
        Write newick file after the analysis
        
    """
    

    post_tree_path = path_out+'/'+clean_name_tree+'.tree_annot.nw'
    t.write(outfile=post_tree_path, props=all_props, format_root_node = True)

    return post_tree_path


def write_tree_for_minvar_rootin(t, tmpdir):

    input_tree_minvar = tmpdir+'input_tree_minvar.nw'
    t.write(outfile=input_tree_minvar, format_root_node = True)

    return input_tree_minvar



def get_sci_name(taxonomy_db, taxa):

    if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
        sci_name_taxa =  taxonomy_db.get_taxid_translator([taxa])[int(taxa)]
    elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy': 
        sci_name_taxa = taxa

    return sci_name_taxa


def get_lineage(taxonomy_db, taxid):
    if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
        lin = taxonomy_db.get_lineage(taxid)
            
    elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy': 
        if lca_tree != 'r_root':
            lin =  taxonomy_db.get_name_lineage([taxid])[0][taxid]

    return lin


def get_rank(taxonomy_db, taxid):
    if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
        rank = clean_string(taxonomy_db.get_rank([taxid])[taxid])

    elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
        rank = get_gtdb_rank(taxid)

    return rank


def get_lca_node_name(taxonomy_db, taxid):
    if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
        lca_node_name = taxonomy_db.get_taxid_translator([taxid])[taxid]

    elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
        lca_node_name = lca_node


    return lca_node_name


def update_sp_per_level_in_node(sp_per_level_in_node, taxonomy_db, l):

    
    if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
        for tax in l.props.get('lineage').split('|'):
            sp_per_level_in_node[tax].add(str(l.props.get('taxid')))
    elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy': 
        for tax in l.props.get('named_lineage').split('|'):
            sp_per_level_in_node[tax].add(str(l.props.get('taxid')))

    