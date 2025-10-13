

from ete4 import PhyloTree, NCBITaxa, GTDBTaxa
from collections import defaultdict
import os
import json
import sys


def load_all_input_files(args):
    """
    Consolidates loading of the gene tree, taxonomy database, reference species tree, 
    and the species counter for the OGD pipeline.
    
    Returns:
        A tuple: (gene_tree, taxonomy_db, ref_species_tree, level2sp_mem)
    """
    print("    - Loading inputs and databases...")
    
    # 1. Load Taxonomy Database
    taxonomy_db = load_taxonomy(taxonomy=args.taxonomy_type, user_taxonomy=args.user_taxonomy)
    
    # 2. Load Gene Tree
    t = load_gene_tree(tree_path=args.tree, sp_delimitator=args.sp_delim)
    
    # 3. Load/Generate Reference Species Tree
    reftree = load_reftree(rtree_path=args.reftree, t=t, taxonomy_db=taxonomy_db)
    
    # 4. Load/Generate Taxonomy Counter
    level2sp_mem = load_taxonomy_counter(reftree=reftree, user_taxonomy_counter=args.user_taxonomy_counter)

    return t, taxonomy_db, reftree, level2sp_mem




def load_gene_tree(tree_path, sp_delimitator ):
    """
    Loads Gene Tree from file, resolves polytomies, and sets the species naming function.
    """
    sys.stdout.write(f"        - Load gene tree: {os.path.basename(tree_path)}\n")
    
    
    t = PhyloTree(tree_path)
    
    # Resolving polytomies is critical step.
    t.resolve_polytomy() 
   
    # Set the species naming function 
    t.set_species_naming_function(lambda node: node.name.split(sp_delimitator)[0])

    return t


def load_taxonomy(taxonomy, user_taxonomy):
    """
    Loads taxonomy database (NCBI or GTDB) from local server or user file.
    """
    sys.stdout.write(f"        - Load taxonomy: {taxonomy}\n")

    if taxonomy == 'NCBI':
        # NCBITaxa automatically handles whether to use a custom file or default DB
        taxonomy_db = NCBITaxa(user_taxonomy, memory=True)
    elif taxonomy == 'GTDB':
        
        taxonomy_db = GTDBTaxa(user_taxonomy)
    else:
        raise ValueError(f"Unsupported taxonomy type: {taxonomy}")

    return taxonomy_db


def load_reftree(rtree_path, t, taxonomy_db):
    """
    Gets the reference species tree either from a user-provided file or generated 
    from the gene tree's species set.
    """
    
    if rtree_path:
        sys.stdout.write("        - Load reftree: from user provided file\n")
        
        reftree = PhyloTree(rtree_path, parser=0) 
    else:
        sys.stdout.write("        - Load reftree: generated from gene tree\n")
        reftree = get_reftree(t, taxonomy_db)

    # Annotate the species tree with tax IDs/names from the database
    taxonomy_db.annotate_tree(reftree, taxid_attr="name")
    
    return reftree


def get_reftree(t, taxonomy_db):
    """
    Creates reference tree (species tree) from the species present in the gene tree.
    """
    
    taxid_list = t.get_species() 
    reftree = taxonomy_db.get_topology(taxid_list)
    
    return reftree


def load_taxonomy_counter(reftree, user_taxonomy_counter):
    """
    Gets the taxonomy counter (number of species per level) from a user file or 
    generates it from the reference tree.
    """
    
    if user_taxonomy_counter:
        sys.stdout.write("        - Load taxonomy counter: from user file\n")

        if isinstance(user_taxonomy_counter, dict):
            
            level2sp_mem = user_taxonomy_counter
        else:
            try:
                with open(user_taxonomy_counter, 'r') as levels:
                    level2sp_mem = json.load(levels)
            except Exception as e:
                sys.stderr.write(f"Error loading taxonomy counter JSON: {e}\n")
                sys.exit(1)
    else:
        sys.stdout.write("        - Load taxonomy counter: generated from reference tree\n")
        level2sp_mem = get_taxonomy_counter(reftree)

    return level2sp_mem


def get_taxonomy_counter(reftree):
    """
    Creates Taxonomy counter (taxonomic level mapped to set of species) from the reference tree.
    """

    level2sp_mem = defaultdict(set)
    
    for leaf in reftree:
        
        # Determine which lineage attribute to use based on the leaf name type
        if leaf.name.isdigit():
            # If name is a taxid, use 'lineage' property
            lin = leaf.props.get('lineage') 
        else:
            # If name is a species name (string), use 'named_lineage'
            lin = leaf.props.get('named_lineage') 
        
        # Ensure lineage exists before iterating
        if lin:
            for tax in lin:
                # Store taxid/name as string key, and the species name as a member
                level2sp_mem[str(tax)].add(leaf.name) 
    
    return level2sp_mem