"""
Performs the initial setup and preprocessing of the phylogenetic tree.

This module is responsible for critical initial steps such as rooting the tree,
annotating it with taxonomic information, and preparing the nodes for
subsequent analysis.
"""

import argparse
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Set, Tuple

from ete4 import PhyloTree

# Local application imports
from ogd.utils import (TaxonomyDB, generate_internal_node_name, get_node_depth,
                       sanitize_tree_properties)



# --- Main Setup Function ---

def run_setup(
    tree: PhyloTree,
    tax_db: TaxonomyDB,
    tmpdir: Path,
    args: argparse.Namespace
) -> Tuple[PhyloTree, Set[str], Set[str], int]:
    """
    Executes the entire tree preprocessing pipeline.

    This includes:
    1. Resolving polytomies.
    2. Rooting the tree.
    3. Annotating nodes with taxonomy.
    4. Naming all internal nodes uniquely.

    Args:
        tree: The initial PhyloTree object loaded from the input file.
        tax_db: The initialized taxonomy database.
        tmpdir: Path to the temporary directory.
        args: The parsed command-line arguments.

    Returns:
        A tuple containing:
        - The fully processed and annotated PhyloTree object.
        - A set of all species names found in the tree.
        - A set of all leaf names (members) in the tree.
        - The total count of unique species.
    """
    # 1. Resolve polytomies (essential for many downstream analyses)
    tree.resolve_polytomy()
    
    # 2. Apply the chosen rooting method
    tree = _apply_rooting(tree, args.rooting, tmpdir, args.sp_delim)

    # 3. Annotate the entire tree with taxonomic information
    tax_db.annotate_tree(tree, taxid_attr="species")
    
    # 4. Get total leaf (member) and species counts
    total_members_in_tree = set(tree.leaf_names())
    species_set = tree.get_species()
    total_species_count = len(species_set)

    # 5. Assign a unique, readable name to each internal node
    logging.info("Assigning names to internal nodes...")
    for i, node in enumerate(tree.traverse()):
        if not node.is_leaf:
            # Name combines a generated ID and the node's depth
            node.name = f"{generate_internal_node_name(i)}-{get_node_depth(node)}"

    logging.info(f"Tree setup complete. Rooting: {args.rooting}, "
                 f"Leaves: {len(total_members_in_tree)}, Species: {total_species_count}")
    
    # Return the processed PhyloTree OBJECT, not a string representation.
    return tree, species_set, total_members_in_tree, total_species_count


def _run_minvar(tree: PhyloTree, tmpdir: Path, species_delimiter: str) -> PhyloTree:
    """
    Roots a tree using the external FastRoot.py (MinVar) tool.

    This involves writing the tree to a temporary file, running the external
    script, and then reading the resulting rooted tree back into memory.

    Args:
        tree: The PhyloTree object to be rooted.
        tmpdir: The temporary directory for intermediate files.
        species_delimiter: Delimiter to correctly parse species from leaf names
                           in the re-loaded tree.

    Returns:
        A new, rooted PhyloTree object.
    """
    # 1. Find the FastRoot.py executable
    fastroot_path = shutil.which('FastRoot.py')
    if not fastroot_path:
        logging.error("FastRoot.py not found in system PATH. Cannot perform MinVar rooting.")
        raise FileNotFoundError("FastRoot.py is required for MinVar rooting.")

    # 2. Define temporary file paths
    input_tree_path = tmpdir / "minvar_input.nw"
    output_tree_path = tmpdir / "minvar_output.nw"
    stdout_log = tmpdir / "minvar.stdout"
    stderr_log = tmpdir / "minvar.stderr"

    # 3. Write the unrooted tree to a file
    tree.write(outfile=str(input_tree_path), format_root_node=True)

    # 4. Build and run the command safely
    command = ["python3", fastroot_path, "-i", str(input_tree_path), "-o", str(output_tree_path)]
    
    logging.info(f"Running MinVar rooting with command: {' '.join(command)}")
    try:
        with open(stdout_log, 'w') as f_out, open(stderr_log, 'w') as f_err:
            subprocess.run(
                command, 
                shell=False,  # shell=False is safer
                check=True,   # Raises an exception if the command fails
                stdout=f_out, 
                stderr=f_err
            )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logging.error(f"MinVar rooting failed. Check logs in {tmpdir}. Error: {e}")
        # Return the original tree as a fallback
        return tree

    # 5. Load the newly rooted tree
    rooted_tree = PhyloTree(str(output_tree_path))
    rooted_tree.set_species_naming_function(lambda node: node.name.split(species_delimiter)[0])
    return rooted_tree


def _apply_rooting(tree: PhyloTree, method: str, tmpdir: Path, species_delimiter: str) -> PhyloTree:
    """
    Roots a tree using one of the available methods.

    Args:
        tree: The PhyloTree object to root.
        method: The rooting method ('Midpoint' or 'MinVar').
        tmpdir: Temporary directory, required for MinVar.
        species_delimiter: Delimiter for species names, required for MinVar.

    Returns:
        The rooted PhyloTree object.
    """
    logging.info(f"Applying '{method}' rooting...")
    if method == "Midpoint":
        try:
            midpoint = tree.get_midpoint_outgroup()
            if midpoint:
                tree.set_outgroup(midpoint)
            else:
                logging.warning("Could not find a midpoint for rooting. Tree remains unrooted.")
        except Exception as e:
            logging.error(f"Midpoint rooting failed with an exception: {e}")
    
    elif method == "MinVar":
        tree = _run_minvar(tree, tmpdir, species_delimiter)
    
    else:
        logging.warning(f"No valid rooting method specified ('{method}'). Tree will not be rooted.")
        
    return tree