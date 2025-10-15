
"""
utils.py - A collection of helper functions for the OGD pipeline.
"""

import re
import tempfile
import os
from pathlib import Path
from collections import defaultdict, Counter, OrderedDict
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Set, Union, Tuple

from ete4 import PhyloTree

# --- Constants ---

# For make_name function
VALID_NODE_CHARS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'

# For run_clean_properties function
PROPERTIES_TO_CLEAN = ["sci_name", "lca_node_name", "common_name"]

# For remove_file_extension function
FILE_EXTENSIONS = ('.faa', '.nw', '.fa', '.fasta')

# For clean_string function - pre-compiled for efficiency
INVALID_CHARS_PATTERN = re.compile(r"['\[\]=\-:]")

# For get_gtdb_rank function
GTDB_RANK_MAP = {
    'r_root': 'r_root', 'd': 'Domain', 'p': 'Phylum',
    'c': 'Class', 'o': 'Order', 'f': 'Family',
    'g': 'Genus', 's': 'Species'
}


# --- Taxonomy Handling Abstraction ---
# Note: For even better organization, these classes could be moved to their own
# file, for example, `taxonomy_handlers.py`.

class BaseTaxonomyHandler(ABC):
    """Abstract base class for taxonomy database handlers."""
    def __init__(self, db_object: Any):
        self.db = db_object

    @abstractmethod
    def get_lineage(self, taxid: Union[int, str]) -> List[Union[int, str]]:
        """Returns the lineage for a given taxid."""
        pass

    @abstractmethod
    def get_rank(self, taxid: Union[int, str]) -> str:
        """Returns the rank for a given taxid."""
        pass

    @abstractmethod
    def get_scientific_name(self, taxid: Union[int, str]) -> str:
        """Returns the scientific name for a given taxid."""
        pass

    def get_lca_node_name(self, taxid: Union[int, str]) -> str:
        """By default, the scientific name is used for the LCA node name."""
        return self.get_scientific_name(taxid)
        
    @abstractmethod
    def update_species_per_level(self, node: PhyloTree, sp_per_level: Dict[str, set]) -> None:
        """Updates a dictionary mapping taxonomic levels to species."""
        pass

class NCBITaxonomyHandler(BaseTaxonomyHandler):
    """Handler for NCBI Taxonomy database."""
    def get_lineage(self, taxid: int) -> List[int]:
        return self.db.get_lineage(taxid)
    
    def get_rank(self, taxid: int) -> str:
        rank_map = self.db.get_rank([taxid])
        return clean_string(rank_map.get(taxid, 'Unk'))

    def get_scientific_name(self, taxid: int) -> str:
        name_map = self.db.get_taxid_translator([taxid])
        return name_map.get(int(taxid), 'Unknown')
        
    def update_species_per_level(self, node: PhyloTree, sp_per_level: Dict[str, set]) -> None:
        lineage_str = node.props.get('lineage', '')
        taxid_str = str(node.props.get('taxid'))
       
        for tax in lineage_str:
            sp_per_level[str(tax)].add(taxid_str)

class GTDBTaxonomyHandler(BaseTaxonomyHandler):
    """Handler for GTDB Taxonomy database."""
    def get_lineage(self, taxid: str) -> List[str]:
        if taxid == 'r_root':
            return []
        # Assumes get_name_lineage returns a list with one dictionary
        return self.db.get_name_lineage([taxid])[0].get(taxid, [])
    
    def get_rank(self, taxid: str) -> str:
        return get_gtdb_rank(taxid)

    def get_scientific_name(self, taxid: str) -> str:
        # For GTDB, the taxid (name) is the scientific name
        return taxid
        
    def update_species_per_level(self, node: PhyloTree, sp_per_level: Dict[str, set]) -> None:
        lineage_str = node.props.get('named_lineage', '')
        taxid_str = str(node.props.get('taxid'))
        for tax in lineage_str.split('|'):
            sp_per_level[tax].add(taxid_str)

# --- Helper Classes ---

class OrderedCounter(Counter, OrderedDict):
    """Counter that remembers the order elements are first seen."""
    def __repr__(self):
        return f'{self.__class__.__name__}({OrderedDict(self)})'
    def __reduce__(self):
        return self.__class__, (OrderedDict(self),)


# --- File and Directory Management ---

def create_tmp(path_out: Path) -> Path:
    """
    Creates a temporary directory.

    Tries to create it inside the results folder (`path_out`), falling back
    to the system's default temporary location on failure.

    Args:
        path_out (Path): Path to the desired output directory.

    Returns:
        Path: A Path object for the created temporary directory.
    """
    try:
        results_abs_path = path_out.resolve()
        tmpdir = tempfile.mkdtemp(dir=results_abs_path)
    except Exception as e:
        print(f"Warning: Failed to create tmp dir in '{path_out}', using system default. Error: {e}")
        tmpdir = tempfile.mkdtemp(prefix="ogd_")
    
    return Path(tmpdir)

def run_write_post_tree(t: PhyloTree, clean_name_tree: str, path_out: Path, all_props: Set[str]) -> Path:
    """Writes the annotated newick tree after analysis."""
    post_tree_path = path_out / f'{clean_name_tree}.tree_annot.nw'
    t.write(outfile=str(post_tree_path), props=list(all_props), format_root_node=True)
    return post_tree_path

def write_tree_for_minvar_rooting(t: PhyloTree, tmpdir: Path) -> Path:
    """Writes a clean newick tree required for MinVar rooting."""
    input_tree_minvar = tmpdir / 'input_tree_minvar.nw'
    t.write(outfile=str(input_tree_minvar), format_root_node=True)
    return input_tree_minvar


# --- String and Name Parsing ---

def parse_taxid(node: PhyloTree) -> str:
    """Parses taxid assuming a 'taxid.genename' format."""
    return node.name.split('.')[0]

def parse_taxid_gtdb(node: PhyloTree) -> str:
    """Parses taxid assuming a 'taxid@genename' format."""
    # TODO: add argument for split gen name
    return node.name.split('@')[0]

def clean_string(input_string: str) -> str:
    """Removes special characters and replaces spaces with underscores."""
    clean = INVALID_CHARS_PATTERN.sub("", input_string)
    return clean.replace(' ', '_')

def remove_file_extension(name: str) -> str:
    """Removes common sequence or tree file extensions from a string."""
    stem = name
    # Loop to remove multiple extensions like .fasta.faa
    while stem.endswith(FILE_EXTENSIONS):
        for ext in FILE_EXTENSIONS:
            if stem.endswith(ext):
                stem = stem[:-len(ext)]
                break
    return stem
    
def get_gtdb_rank(gtdb_code: str) -> str:
    """Gets the rank name from a GTDB code (e.g., 'p__Proteobacteria')."""
    prefix = gtdb_code.split('__')[0]
    return GTDB_RANK_MAP.get(prefix, 'Unk')


# --- Tree Processing ---

def get_newick(t: PhyloTree, all_props: Set[str]) -> str:
    """Returns tree in newick format with specified annotations."""
    return t.write(props=list(all_props), format_root_node=True)

def run_clean_properties(tree: Union[str, PhyloTree]) -> Tuple[PhyloTree, Set[str]]:
    """
    Cleans specified node property strings in a tree.

    Iterates through all nodes and cleans special characters from a defined
    list of properties. Also collects all unique property keys found in the tree.

    Args:
        tree (Union[str, PhyloTree]): An ete4 PhyloTree object or a newick string.

    Returns:
        Tuple[PhyloTree, Set[str]]: The modified tree and a set of all property keys.
    """
    t = PhyloTree(tree) if isinstance(tree, str) else tree
    
    all_prop_keys = set()
    for node in t.traverse():
        all_prop_keys.update(node.props.keys())
        
        for prop_name in PROPERTIES_TO_CLEAN:
            prop_value = node.props.get(prop_name)
            if prop_value and isinstance(prop_value, str):
                node.props[prop_name] = clean_string(prop_value)

    return t, all_prop_keys

def make_name(i: int) -> str:
    """Creates a short, unique name for an internal node based on an integer."""
    if i < 0: return ""
    
    name = ''
    num_chars = len(VALID_NODE_CHARS)
    while True:
        name = VALID_NODE_CHARS[i % num_chars] + name
        i = i // num_chars
        if i == 0: break
        i -= 1 # Decrement to handle the base conversion correctly
    return name

def get_depth(node: PhyloTree) -> int:
    """
    Calculates the depth of a node from the root.
    Depth is the number of nodes from the root to the target node.
    """
    depth = 0
    current = node
    while current.up:
        depth += 1
        current = current.up
    return depth

# --- Logic and Thresholding ---

def determine_so_threshold(lin_lca: List[Union[int, str]], args: Any, handler: BaseTaxonomyHandler) -> float:
    """
    Determines the Species Overlap (SO) threshold based on LCA lineage.

    Args:
        lin_lca: Lineage of the Last Common Ancestor.
        args: Command-line arguments object with so_* thresholds.
        handler: The taxonomy handler being used (NCBI or GTDB).

    Returns:
        The appropriate Species Overlap threshold (float).
    """
    so_to_use = args.so_all

    if isinstance(handler, NCBITaxonomyHandler):
        # Map tax IDs to their corresponding SO threshold argument
        NCBI_MAP = {
            2759: args.so_euk,   # Eukaryota
            2:    args.so_bact,  # Bacteria
            2157: args.so_arq    # Archaea
        }
        for tax_id, so_threshold in NCBI_MAP.items():
            if so_threshold is not None and tax_id in lin_lca:
                so_to_use = so_threshold
                break

    elif isinstance(handler, GTDBTaxonomyHandler):
        # Map domain names to their SO threshold argument
        GTDB_MAP = {
            'd__Bacteria': args.so_bact,
            'd__Archaea': args.so_arq,
            'root': args.so_all
        }
        for tax_name, so_threshold in GTDB_MAP.items():
            if so_threshold is not None and tax_name in lin_lca:
                so_to_use = so_threshold
                break
        
        # Handle the special case 'Empty' for GTDB
        if so_to_use == args.so_all and 'Empty' in lin_lca:
            so_to_use = 0.0
            
    return so_to_use