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
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Set, Union

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
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Set, Union

from ete4 import GTDBTaxa, NCBITaxa, PhyloTree

# --- Type Aliases ---
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

    Tries to create the directory within the specified output path. If it fails
    (e.g., due to permissions), it falls back to the system's default temporary
    directory.

    Args:
        output_path: The Path object for the main output directory.

    Returns:
        The Path object for the created temporary directory.
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
    # Remove problematic characters: ', [, ], =, -, :
    sanitized = re.sub(r"['\[\]=\-:]", "", text)
    # Replace spaces with underscores
    return sanitized.replace(' ', '_')

def remove_file_extension(filename: str) -> str:
    """
    Removes known extensions from a filename string, handling multiple extensions.
    """
    p = Path(filename)
    # Sequentially remove suffixes until no known extension is found
    while p.suffix in {'.gz', '.bz2', '.zip', '.faa', '.nw', '.fa', '.fasta'}:
        p = p.with_suffix('')
    return p.name


# --- Taxonomy Abstraction ---

def is_gtdb(tax_db: TaxonomyDB) -> bool:
    """Checks if the provided taxonomy database is GTDB."""
    return isinstance(tax_db, GTDBTaxa)

def get_scientific_name(tax_db: TaxonomyDB, taxid: Union[int, str]) -> str:
    """Gets the scientific name for a given taxid, abstracting the db type."""
    if is_gtdb(tax_db):
        return str(taxid)
    return tax_db.get_taxid_translator([int(taxid)])[int(taxid)]

def get_rank(tax_db: TaxonomyDB, taxid: Union[int, str]) -> str:
    """Gets the rank for a given taxid, abstracting the db type."""
    if is_gtdb(tax_db):
        rank_code = str(taxid).split('__')[0]
        rank_map = {
            'd': 'Domain', 'p': 'Phylum', 'c': 'Class', 'o': 'Order',
            'f': 'Family', 'g': 'Genus', 's': 'Species', 'r_root': 'root'
        }
        return rank_map.get(rank_code, 'Unknown')
    raw_rank = tax_db.get_rank([int(taxid)])[int(taxid)]
    return sanitize_string_for_newick(raw_rank)

def get_lca_node(species_set: Set[str], tax_db: TaxonomyDB) -> str:
    """
    Finds the Last Common Ancestor (LCA) for a given set of species taxids.
    """
    if not species_set:
        return 'Unk'

    # Count occurrences of each taxon in all lineages
    lineage_counts = OrderedCounter()
    for species_id in species_set:
        lineage = get_lineage(tax_db, species_id) 
        if lineage:
            lineage_counts.update(lineage)

    # The LCA is the most recent (last) taxon present in ALL lineages
    num_species = len(species_set)
    lca_candidates = [taxon for taxon, count in lineage_counts.items() if count == num_species]
    
    return lca_candidates[-1] if lca_candidates else 'Unk'



def get_lineage(taxonomy_db, taxid):
    if 'ncbi_taxonomy' in str(taxonomy_db):
        lin = taxonomy_db.get_lineage(taxid)
        
            
    elif 'gtdb_taxonomy' in str(taxonomy_db) : 
        #if lca_tree != 'r_root':
        lin =  taxonomy_db.get_name_lineage([taxid])[0][taxid]

    
    return lin

def update_sp_per_level_in_node(sp_per_level_in_node, taxonomy_db, l):

    
    if 'ncbi_taxonomy' in str(taxonomy_db):
        
        for tax in l.props.get('lineage'):
            sp_per_level_in_node[tax].add(l.props.get('taxid'))

    elif 'gtdb_taxonomy' in str(taxonomy_db): 
        for tax in l.props.get('named_lineage'):
            sp_per_level_in_node[tax].add((l.props.get('taxid')))

    
def get_lca_node_name(taxonomy_db, taxid):
    if 'ncbi_taxonomy' in str(taxonomy_db):
        lca_node_name = taxonomy_db.get_taxid_translator([taxid])[taxid]

    elif 'gtdb_taxonomy' in str(taxonomy_db): 
        lca_node_name = lca_node


    return lca_node_name

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
        i -= 1
    return name


def get_node_depth(node):

    """
        Get depth of internal nodes
        Depth = number nodes from root node to target node
    """
    depth = 0
    while node is not None:
        depth += 1
        node = node.up
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
            if prop_value := node.props.get(prop_name):
                node.props[prop_name] = sanitize_string_for_newick(prop_value)

    return phylo_tree, all_prop_keys

def write_annotated_tree(tree: PhyloTree, output_dir: Path, filename_base: str, props: Set[str]) -> Path:
    """
    Writes the final annotated tree to a Newick file.
    """
    output_path = output_dir / f"{filename_base}.tree_annot.nw"
    # Convert props set to list for ete4 compatibility
    tree.write(outfile=str(output_path), props=list(props), format_root_node=True)
    logging.info(f"Annotated tree written to: {output_path}")
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
        if args.so_euk is not None and 2759 in lineage_set: return args.so_euk
        if args.so_bact is not None and 2 in lineage_set: return args.so_bact
        if args.so_arq is not None and 2157 in lineage_set: return args.so_arq
            
    return args.so_all