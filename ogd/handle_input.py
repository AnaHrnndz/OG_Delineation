import argparse
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Tuple, Union

from ete4 import GTDBTaxa, NCBITaxa, PhyloTree
from ogd.utils import TaxonomyDB

# --- Main Orchestrator Function ---

def load_all_input_files(args: argparse.Namespace) -> Tuple[PhyloTree, TaxonomyDB, PhyloTree, Dict[str, set]]:
    """
    Load of all primary input files and databases for the pipeline.

    Args:
        args: The command-line arguments object from argparse.

    Returns:
        A tuple containing:
        - The loaded gene tree.
        - The initialized taxonomy database object.
        - The loaded or generated reference species tree.
        - The taxonomy counter dictionary.
    """
    
    # 1. Load Taxonomy Database
    taxonomy_db = _load_taxonomy(
        taxonomy_type=args.taxonomy_type, 
        user_taxonomy_path=args.user_taxonomy
    )
    
    # 2. Load Gene Tree
    gene_tree = _load_gene_tree(
        tree_path=args.tree, 
        species_delimiter=args.sp_delim
    )
    
    # 3. Load or Generate Reference Species Tree
    ref_species_tree = _load_reftree(
        rtree_path=args.reftree, 
        gene_tree=gene_tree, 
        taxonomy_db=taxonomy_db
    )
    
    # 4. Load or Generate Taxonomy Counter
    level2species = _load_taxonomy_counter(
        ref_species_tree=ref_species_tree, 
        user_counter_path=args.user_taxonomy_counter
    )

    logging.info("All input data loaded successfully.")
    return gene_tree, taxonomy_db, ref_species_tree, level2species

# ---  Helper Functions ---

def _load_gene_tree(tree_path: Path, species_delimiter: str) -> PhyloTree:
    """
    Loads a gene tree, resolves polytomies, and sets the species naming function.
    """
    logging.info(f"Loading gene tree from: {tree_path.name}")
    try:
        t = PhyloTree(open(tree_path))
    except Exception as e:
        logging.error(f"Failed to parse the gene tree file at '{tree_path}'.")
        logging.error(f"Error details: {e}")
        logging.error("Please ensure the file is a valid Newick format.")
        sys.exit(1)

    #  Resolve polytomies (essential for many downstream analyses)
    t.resolve_polytomy()
    
    # This function sets how species are identified from leaf names
    t.set_species_naming_function(lambda node: node.name.split(species_delimiter)[0])

    return t

def _load_taxonomy(taxonomy_type: str, user_taxonomy_path: Path = None) -> TaxonomyDB:
    """
    Loads the specified taxonomy database (NCBI or GTDB).
    """
    logging.info(f"Initializing {taxonomy_type} taxonomy database...")
    
    db_path = str(user_taxonomy_path) if user_taxonomy_path else None

    try:
        if taxonomy_type == 'NCBI':
            return NCBITaxa(dbfile=db_path, memory=True)
        elif taxonomy_type == 'GTDB':
            return GTDBTaxa(dbfile=db_path)
        else:
            raise ValueError(f"Unsupported taxonomy type: {taxonomy_type}")
    except Exception as e:
        logging.error(f"Failed to load taxonomy database from: {db_path or 'default path'}")
        logging.error(f"Error details: {e}")
        sys.exit(1)

def _load_reftree(rtree_path: Path, gene_tree: PhyloTree, taxonomy_db: TaxonomyDB) -> PhyloTree:
    """
    Loads a reference species tree from a file or generates it from the gene tree.
    """
    if rtree_path and rtree_path.is_file():
        logging.info(f"Loading reference species tree from user file: {rtree_path.name}")
        try:
            ref_tree = PhyloTree(open(rtree_path))
        except Exception as e:
            logging.error(f"Failed to parse the reference tree file at '{rtree_path}'.")
            logging.error(f"Error details: {e}")
            sys.exit(1)
    else:
        if rtree_path:
             logging.warning(f"User-provided reftree not found at '{rtree_path}'.")
        logging.info("Generating reference species tree from gene tree content.")
        species_list = gene_tree.get_species()
        ref_tree = taxonomy_db.get_topology(species_list)

    # Annotate the tree with lineage and other info from the database.
    taxonomy_db.annotate_tree(ref_tree, taxid_attr="name")
    return ref_tree

def _load_taxonomy_counter(ref_species_tree: PhyloTree, user_counter_path: Path = None) -> Dict[str, set]:
    """
    Loads a taxonomy counter from a file or generates it from the reference tree.
    """
    if user_counter_path and user_counter_path.is_file():
        logging.info(f"Loading taxonomy counter from user file: {user_counter_path.name}")
        try:
            with user_counter_path.open('r') as f:
                json_data = json.load(f)
                return {level: set(species) for level, species in json_data.items()}
        except (json.JSONDecodeError, IOError) as e:
            logging.error(f"Could not load or parse taxonomy counter from '{user_counter_path}'.")
            logging.error(f"Error details: {e}")
            sys.exit(1)
    else:
        logging.info("Generating taxonomy counter from reference tree.")
        return _get_taxonomy_counter(ref_species_tree)


def _get_taxonomy_counter(ref_species_tree: PhyloTree) -> Dict[str, set]:
    """
    Generates a taxonomy counter from an annotated reference species tree.
    Keys (internal taxonomic IDs) and values (species taxid) are stored as strings.
    """
    level2species = defaultdict(set)
    for leaf in ref_species_tree:
        lineage = leaf.props.get('lineage') or leaf.props.get('named_lineage')
        if lineage:
            for tax_level in lineage:
                level2species[str(tax_level)].add(leaf.name)
    
    return dict(level2species)