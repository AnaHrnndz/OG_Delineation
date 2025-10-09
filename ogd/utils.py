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


def determine_so_threshold(taxonomy_db, lin_lca, args):
    """
    Determines the Species Overlap (SO) threshold to use based on the 
    taxonomy type (NCBI or GTDB) and the ancestral lineage (lin_lca) of the LCA.
    
    Args:
        taxonomy_db: The taxonomy database object (e.g., NCBI Taxonomy object).
        lin_lca: A list of tax IDs (int) or tax names (str) in the LCA lineage.
        args: Object containing the SO thresholds (so_euk, so_bact, so_arq, so_all, etc.).
        
    Returns:
        The Species Overlap threshold (float) to use.
    """
    
    # 1. Initialize with the default threshold
    so_to_use = args.so_all

    # Use the string representation to determine the taxonomy type
    taxonomy_name = str(taxonomy_db).lower()

    # 2. Define search maps for tax IDs/names to SO arguments
    
    # Map for NCBI Taxonomy (Tax IDs)
    NCBI_MAP = {
        2759: args.so_euk,   # Eukaryota
        2:    args.so_bact,   # Bacteria
        2157: args.so_arq    # Archaea
    }

    # Map for GTDB Taxonomy (Domain names)
    GTDB_MAP = {
        'd__Bacteria': args.so_bact,
        'd__Archaea': args.so_arq,
        'root': args.so_all 
    }

    # 3. Selection logic based on taxonomy database
    
    if 'ncbi_taxonomy' in taxonomy_name:
        
        # Iterate over the NCBI map to find the most specific available threshold
        for tax_id, so_threshold in NCBI_MAP.items():
            
            # Condition: User provided the threshold AND the tax_id is in the lineage
            if so_threshold is not None and tax_id in lin_lca:
                so_to_use = so_threshold
                break 

    elif 'gtdb_taxonomy' in taxonomy_name:
        
        # Iterate over the GTDB map to find the matching threshold
        for tax_name, so_threshold in GTDB_MAP.items():
            if tax_name in lin_lca:
                so_to_use = so_threshold
                break
                
        # Handle the special case 'Empty' for GTDB
        if so_to_use == args.so_all and 'Empty' in lin_lca:
            so_to_use = 0.0

    return so_to_use


def get_newick(t, all_props):

    """
        Return tree in newick format with annotations
    """

    t = t.write(props=all_props, format_root_node=True)
    return t


def run_clean_properties(tree):

    # 1. Input Handling: Ensure the input is a PhyloTree object
    if isinstance(tree, str):
        
        t = PhyloTree(tree)
    else:
        t = tree 

    
    PROPERTIES_TO_CLEAN = ["sci_name", "lca_node_name", "common_name"]
    
    all_prop_keys= set()

    # 2. Traverse and Clean Properties
    for node in t.traverse():
        # Update the set of all property keys present in the tree
        all_prop_keys.update(node.props.keys())
        
        # Iterate only through the required properties to clean
        for prop_name in PROPERTIES_TO_CLEAN:
            # Check if the property exists and is not None
            prop_value = node.props.get(prop_name)
            
            if prop_value:
                # Apply the cleaning function and update the property value
                node.props[prop_name] = clean_string(prop_value)

    return t, all_prop_keys


def clean_string(string):
    
    clean_string = re.sub(r"'|\[|\]|\=|\-|\:", "", string)
    return clean_string.replace(' ', '_') 


def remove_file_extension(name):
    
    extensions = ['.faa', '.nw', '.fa', '.fasta']
    
    removed = True
    while removed:
        removed = False
        for ext in extensions:
           
            if name.endswith(ext):
               
                name = name.removesuffix(ext)
                removed = True
               
                
    return name[:len(name)]



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

    