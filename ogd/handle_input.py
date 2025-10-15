# -*- coding: utf-8 -*-
"""
handle_input.py - Loads and preprocesses all input files for the OGD pipeline.
"""

import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Tuple, Dict, Any

from ete4 import PhyloTree, NCBITaxa, GTDBTaxa

TaxonomyDB = Any # Represents NCBITaxa or GTDBTaxa objects

def load_all_input_files(args: Any) -> Tuple[PhyloTree, TaxonomyDB, PhyloTree, Dict[str, set]]:
    """
    Consolidates loading of all primary input files for the pipeline.

    Args:
        args: The command-line arguments object from argparse.

    Returns:
        A tuple containing the loaded gene tree, taxonomy database object,
        reference species tree, and the taxonomy counter dictionary.
    """
    print("    - Loading inputs and databases...")
    
    # 1. Load Taxonomy Database
    taxonomy_db = load_taxonomy(
        taxonomy_type=args.taxonomy_type, 
        user_taxonomy_path=args.user_taxonomy
    )
    
    # 2. Load Gene Tree
    gene_tree = load_gene_tree(
        tree_path=args.tree, 
        sp_delimitator=args.sp_delim
    )
    
    # 3. Load or Generate Reference Species Tree
    ref_species_tree = load_reftree(
        rtree_path=args.reftree, 
        gene_tree=gene_tree, 
        taxonomy_db=taxonomy_db
    )
    
    # 4. Load or Generate Taxonomy Counter
    level2sp_mem = load_taxonomy_counter(
        ref_species_tree=ref_species_tree, 
        user_counter_path=args.user_taxonomy_counter
    )

    return gene_tree, taxonomy_db, ref_species_tree, level2sp_mem


def load_gene_tree(tree_path: Path, sp_delimitator: str) -> PhyloTree:
    """
    Loads a gene tree, resolves polytomies, and sets the species naming function.
    
    Args:
        tree_path (Path): Path to the newick tree file.
        sp_delimitator (str): Character separating species ID from gene ID in leaf names.

    Returns:
        An ete4 PhyloTree object.
    """
    print(f"        - Loading gene tree: {tree_path.name}")
    try:
        t = PhyloTree(str(tree_path))
    except Exception as e:
        print(f"ERROR: Could not parse the gene tree file at '{tree_path}'.\nDetails: {e}", file=sys.stderr)
        sys.exit(1)
        
    # Resolving polytomies is a critical preprocessing step.
    t.resolve_polytomy()
    
    # Set the function to extract species names from leaves.
    t.set_species_naming_function(lambda node: node.name.split(sp_delimitator)[0])

    return t


def load_taxonomy(taxonomy_type: str, user_taxonomy_path: Path = None) -> TaxonomyDB:
    """
    Loads the specified taxonomy database (NCBI or GTDB).

    Args:
        taxonomy_type (str): The type of taxonomy to load ('NCBI' or 'GTDB').
        user_taxonomy_path (Path, optional): Path to a custom taxonomy file. Defaults to None.

    Returns:
        An initialized NCBITaxa or GTDBTaxa object.
    """
    print(f"        - Loading taxonomy: {taxonomy_type}")
    
    db_path = str(user_taxonomy_path) if user_taxonomy_path else None

    if taxonomy_type == 'NCBI':
        # NCBITaxa can handle a None path, using its default database.
        return NCBITaxa(dbfile=db_path, memory=True)
    elif taxonomy_type == 'GTDB':
        return GTDBTaxa(dbfile=db_path)
    else:
        raise ValueError(f"Unsupported taxonomy type: {taxonomy_type}")


def load_reftree(rtree_path: Path, gene_tree: PhyloTree, taxonomy_db: TaxonomyDB) -> PhyloTree:
    """
    Loads a reference species tree from a file or generates it from the gene tree.

    Args:
        rtree_path (Path): Path to a user-provided species tree file.
        gene_tree (PhyloTree): The loaded gene tree.
        taxonomy_db (TaxonomyDB): The initialized taxonomy database.

    Returns:
        The annotated reference species tree.
    """
    if rtree_path:
        print("        - Loading reference tree: from user-provided file")
        ref_tree = PhyloTree(newick=rtree_path)
    else:
        print("        - Loading reference tree: generating from gene tree species")
        species_list = gene_tree.get_species()
        ref_tree = taxonomy_db.get_topology(species_list)

    # Annotate the tree with lineage and other info from the database.
    taxonomy_db.annotate_tree(ref_tree, taxid_attr="name")
    return ref_tree


def load_taxonomy_counter(ref_species_tree: PhyloTree, user_counter_path: Path = None) -> Dict[str, set]:
    """
    Loads a taxonomy counter from a file or generates it from the reference tree.

    Args:
        ref_species_tree (PhyloTree): The annotated reference species tree.
        user_counter_path (Path, optional): Path to a user-provided JSON counter file.

    Returns:
        A dictionary mapping each taxonomic level to a set of species.
    """
    if user_counter_path:
        print("        - Loading taxonomy counter: from user file")
        try:
            with user_counter_path.open('r') as f:
                # Convert list values back to sets after loading from JSON
                json_data = json.load(f)
                return {level: set(species) for level, species in json_data.items()}
        except Exception as e:
            print(f"ERROR: Could not load taxonomy counter from '{user_counter_path}'.\nDetails: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("        - Loading taxonomy counter: generating from reference tree")
        return get_taxonomy_counter(ref_species_tree)


def get_taxonomy_counter(ref_species_tree: PhyloTree) -> Dict[str, set]:
    """
    Generates a taxonomy counter from an annotated reference species tree.

    The counter maps each taxonomic level (e.g., a taxid or a rank name) to the
    set of species belonging to it.

    Args:
        ref_species_tree (PhyloTree): An ete4 tree already annotated by a taxonomy DB.

    Returns:
        A dictionary mapping levels to sets of species names.
    """
    level2sp_mem = defaultdict(set)
    for leaf in ref_species_tree:
        # Check for both possible lineage attributes to be robust
        lineage = leaf.props.get('lineage') or leaf.props.get('named_lineage')
        
        if lineage:
            for tax_level in lineage:
                level2sp_mem[str(tax_level)].add(leaf.name)
                
    return level2sp_mem