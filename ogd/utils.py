"""
Utility functions for the Orthologous Group Delineation (OGD) pipeline.

This module provides a collection of helper functions for tasks such as
file system operations, string manipulation, tree processing, and interacting
with taxonomy databases in a generic way.
"""

import argparse
import json
import logging
import re
import tempfile
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path
from typing import Set, Union

from ete4 import GTDBTaxa, NCBITaxa, PhyloTree

# ---  Aliases ---
TaxonomyDB = Union[NCBITaxa, GTDBTaxa]

# --- Constants ---
# Used for generating unique names for internal nodes
_NODE_NAME_CHARS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
# Node properties to sanitize before writing to Newick format
_PROPERTIES_TO_SANITIZE = ["sci_name", "lca_node_name", "common_name"]


# --- Classes ---

class OrderedCounter(Counter, OrderedDict):
    """A Counter that remembers the order elements are first seen."""
    def __repr__(self):
        return f'{self.__class__.__name__}({OrderedDict(self)!r})'

    def __reduce__(self):
        return self.__class__, (OrderedDict(self),)


# --- File System and String Manipulation ---

def create_temp_dir(output_path: Path) -> Path:
    """
    Creates a temporary directory for intermediate files.
    Tries to create the directory within the specified output path. 
    If it fails, back to the system's default temporary directory.
    """
    try:
        tmpdir = tempfile.mkdtemp(dir=output_path)
        logging.info(f"Temporary directory created: {tmpdir}")
    except OSError as e:
        tmpdir = tempfile.mkdtemp(prefix="ogd_")
        logging.warning(f"Could not create temp dir in '{output_path}'. Using system default: {tmpdir}. Reason: {e}")
    return Path(tmpdir)

def sanitize_string_for_newick(text: str) -> str:
    """Removes characters that can break the Newick format and replaces spaces."""
    if not isinstance(text, str):
        return str(text)
    # Remove problematic characters: 
    sanitized = re.sub(r"['\[\]=\-:]", "", text)
    # Replace spaces with underscores
    return sanitized.replace(' ', '_')

def remove_file_extension(filename: str) -> str:
    """
    Removes known extensions from a filename string, handling multiple extensions.
    """
    p = Path(filename)
    
    while p.suffix in {'.faa', '.nw', '.fa', '.fasta'}:
        p = p.with_suffix('')
    return p.name


# --- Taxonomy Abstraction Layer ---

def is_gtdb(tax_db: TaxonomyDB) -> bool:
    """Checks if the provided taxonomy database is GTDB."""
    return isinstance(tax_db, GTDBTaxa)



def get_rank(tax_db: TaxonomyDB, taxid: Union[int, str]) -> str:
    """Gets the rank for a given taxid, abstracting the db type."""
    if is_gtdb(tax_db):
        rank_code = str(taxid).split('__')[0]
        rank_map = {
            'd': 'Domain', 'p': 'Phylum', 'c': 'Class', 'o': 'Order',
            'f': 'Family', 'g': 'Genus', 's': 'Species', 'r_root': 'root'
        }
        return rank_map.get(rank_code, 'Unknown')
    try:
        raw_rank = tax_db.get_rank([int(taxid)])[int(taxid)]
        return sanitize_string_for_newick(raw_rank)
    except (ValueError, KeyError):
        return "Unknown"

def get_lineage(tax_db: TaxonomyDB, taxid: Union[int, str]) -> list:
    """Gets the lineage for a given taxid, abstracting the db type."""
    # Return empty list for invalid/empty taxids
    if not taxid or taxid == 'Unk' or taxid == 'Empty':
        return []
        
    if is_gtdb(tax_db):
        try:
            # GTDB uses string names as taxids
            return tax_db.get_name_lineage([str(taxid)])[0].get(str(taxid), [])
        except (KeyError, IndexError):
            return [] # Handle if GTDB name is not found
    else:
        try:
            # NCBI uses integer taxids
            return tax_db.get_lineage(int(taxid))
        except (ValueError, KeyError):
            return [] # Handle if taxid is not an int or not found

def get_lca_node(species_set: Set[str], tax_db: TaxonomyDB) -> str:
    """Finds the Last Common Ancestor (LCA) for a given set of species taxids."""
    if not species_set:
        return 'Unk'

    lineage_counts = OrderedCounter()
    for species_id in species_set:
        
        lineage = get_lineage(tax_db, species_id) 
        if lineage:
            lineage_counts.update(lineage)

    if not lineage_counts:
        return 'Unk'

    num_species = len(species_set)
    lca_candidates = [taxon for taxon, count in lineage_counts.items() if count == num_species]
    
    return str(lca_candidates[-1]) if lca_candidates else 'Unk'

# def get_lca_node_name(tax_db: TaxonomyDB, taxid: Union[int, str]) -> str:
    # """
    # Gets the scientific name for a given taxid.
    # NOTE: This is now an alias for get_scientific_name to fix a bug
    # and remove redundant logic.
    # """
    # return get_scientific_name(tax_db, taxid)


def get_scientific_name(tax_db: TaxonomyDB, taxid: Union[int, str]) -> str:
    """Gets the scientific name for a given taxid, abstracting the db type."""
    if is_gtdb(tax_db):
        return str(taxid)
    try:
        # NCBI taxids are integers
        return tax_db.get_taxid_translator([int(taxid)])[int(taxid)]
    except (ValueError, KeyError):
        return "Unknown" # Handle if taxid is invalid

def update_sp_per_level_in_node(sp_per_level_in_node: defaultdict, tax_db: TaxonomyDB, leaf_node: PhyloTree):
    """Updates a species counter dict based on a leaf node's lineage."""
    
    if is_gtdb(tax_db):
        lineage = leaf_node.props.get('named_lineage', [])
    else:
        lineage = leaf_node.props.get('lineage', [])
    
    taxid = leaf_node.props.get('taxid')
    if not taxid or not lineage:
        return # Skip leaves with no data
        
    for tax in lineage:
        # Ensure keys and values are strings for consistency
        sp_per_level_in_node[str(tax)].add(str(taxid))


# --- Tree Processing ---

def generate_internal_node_name(counter: int) -> str:
    """
    Creates a unique, short name for an internal node based on a counter.
    Example: 0->A, 1->B, 51->Z, 52->BA
    """
    if counter < 0: return ""
    name = ''
    i = counter
    while True:
        name = _NODE_NAME_CHARS[i % len(_NODE_NAME_CHARS)] + name
        i //= len(_NODE_NAME_CHARS)
        if i == 0: break
        i -= 1 # Decrement for base conversion
    return name

def get_node_depth(node: PhyloTree) -> int:
    """
    Get depth of nodes (number of nodes from root to target).
    Root node has a depth of 1.
    """
    depth = 0
    current = node
    while current is not None:
        depth += 1
        current = current.up
    return depth

def sanitize_tree_properties(tree: Union[PhyloTree, str]) -> tuple[PhyloTree, Set[str]]:
    """
    Cleans specific node properties in a tree to ensure Newick compatibility.
    """
    phylo_tree = tree if isinstance(tree, PhyloTree) else PhyloTree(tree)
    
    all_prop_keys = set()
    for node in phylo_tree.traverse():
        all_prop_keys.update(node.props.keys())
        for prop_name in _PROPERTIES_TO_SANITIZE:
            # Use .get() to safely check for property existence
            if prop_value := node.props.get(prop_name):
                node.props[prop_name] = sanitize_string_for_newick(prop_value)

    return phylo_tree, all_prop_keys

def write_annotated_tree(tree: PhyloTree, output_dir: Path, filename_base: str, props: Set[str]) -> Path:
    """
    Writes the final annotated tree to a Newick file.
    """
    output_path = output_dir / f"{filename_base}.tree_annot.nw"
    try:
        # Convert props set to list for ete4 compatibility
        tree.write(outfile=str(output_path), props=list(props), format_root_node=True)
        logging.info(f"Annotated tree written to: {output_path}")
    except IOError as e:
        logging.error(f"Failed to write annotated tree to {output_path}: {e}")
    return output_path


# --- Algorithm-specific Helpers ---

def determine_so_threshold(tax_db: TaxonomyDB, lineage: list, args: argparse.Namespace) -> float:
    """
    Determines the Species Overlap (SO) threshold based on the node's lineage.
    """
    lineage_set = set(lineage)
    
    if is_gtdb(tax_db):
        if args.so_bact is not None and 'd__Bacteria' in lineage_set: return args.so_bact
        if args.so_arq is not None and 'd__Archaea' in lineage_set: return args.so_arq
    else: # NCBI
        # Use integer taxids for NCBI-specific checks
        int_lineage_set = {int(t) for t in lineage_set if str(t).isdigit()}
        if args.so_euk is not None and 2759 in int_lineage_set: return args.so_euk
        if args.so_bact is not None and 2 in int_lineage_set: return args.so_bact
        if args.so_arq is not None and 2157 in int_lineage_set: return args.so_arq
            
    return args.so_all