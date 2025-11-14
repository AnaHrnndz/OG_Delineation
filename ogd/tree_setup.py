import argparse
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Set, Tuple

from ete4 import PhyloTree

# Local application imports

from ogd.utils import TaxonomyDB, generate_internal_node_name, get_node_depth, sanitize_tree_properties

# --- Main Setup Function ---

def run_setup(
    tree: PhyloTree,
    tax_db: TaxonomyDB,
    tmpdir: Path,
    args: argparse.Namespace
) -> Tuple[PhyloTree, Set[str], Set[str], int]:
    """
    Executes the entire tree preprocessing pipeline.

    Includes:  rooting the tree, add taxonomic annotation, and naming internal nodes.
    """
    logging.info("Starting tree preprocessing...")
    
    
    # 1. Apply the chosen rooting method
    logging.info(f"Applying '{args.rooting}' rooting method...")
    tree = _apply_rooting(tree, args.rooting, tmpdir, args.sp_delim)

    # 2. Annotate the entire tree with taxonomic information
    logging.info("Annotating tree nodes with taxonomy...")
    tax_db.annotate_tree(tree, taxid_attr="species")
    
    # 3. Get total leaf (member) and species counts
    total_members_in_tree = set(tree.leaf_names())
    species_set = tree.get_species() # Retrieves species based on the naming function
    total_species_count = len(species_set)

    # 4. Assign a unique, readable name to each internal node
    logging.info("Assigning unique names to internal nodes...")
    node_counter = 0
    for node in tree.traverse():
        if not node.is_leaf:
            # Name combines a generated ID and the node's depth
            node.name = f"{generate_internal_node_name(node_counter)}-{get_node_depth(node)}"
            node_counter += 1 # Ensure unique IDs

    logging.info(f"Tree setup complete. Rooting applied: {args.rooting}. "
                 f"Total leaves: {len(total_members_in_tree)}, Total species: {total_species_count}")
    
    # Return the processed PhyloTree  
    return tree, species_set, total_members_in_tree, total_species_count


# ---  Helper Functions for Rooting ---

def _apply_rooting(tree: PhyloTree, method: str, tmpdir: Path, species_delimiter: str) -> PhyloTree:
    """
    Roots a tree using either Midpoint or MinVar methods.
    """
    if method == "Midpoint":
        try:
            midpoint = tree.get_midpoint_outgroup()
            if midpoint:
                tree.set_outgroup(midpoint)
                logging.info("Midpoint rooting applied successfully.")
            else:
                logging.warning("Could not determine a midpoint for rooting. Tree remains unrooted.")
        except Exception as e:
            # Catch potential errors during midpoint calculation or setting outgroup
            logging.error(f"Midpoint rooting failed with an exception: {e}")
            logging.warning("Tree remains unrooted.")
    
    elif method == "MinVar":
        # _run_minvar handles its own logging and error reporting
        tree = _run_minvar(tree, tmpdir, species_delimiter)
    
    elif method is None or method == '':
        logging.info("No rooting method specified. Tree will not be rooted.")
    else:
        # Should not happen if argparse choices are set correctly, but good practice
        logging.warning(f"Unknown rooting method specified: '{method}'. Tree will not be rooted.")
        
    return tree


def _run_minvar(tree: PhyloTree, tmpdir: Path, species_delimiter: str) -> PhyloTree:
    """
    Roots a tree using the external FastRoot.py (MinVar) tool.
    """
    # 1. Find the FastRoot.py executable
    fastroot_path = shutil.which('FastRoot.py')
    if not fastroot_path:
        logging.error("FastRoot.py not found in system PATH. Cannot perform MinVar rooting.")
        raise FileNotFoundError("FastRoot.py is required for MinVar rooting.")

    # 2. Define temporary file paths using pathlib
    input_tree_path = tmpdir / "minvar_input.nw"
    output_tree_path = tmpdir / "minvar_output.nw"
    stdout_log = tmpdir / "minvar.stdout"
    stderr_log = tmpdir / "minvar.stderr"

    # 3. Write the unrooted tree to a file
    try:
        tree.write(outfile=str(input_tree_path), format_root_node=True)
    except IOError as e:
        logging.error(f"Failed to write temporary tree for MinVar: {e}")
        return tree # Return original tree if writing fails

    # 4. Build and run the command safely
    command = ["python3", fastroot_path, "-i", str(input_tree_path), "-o", str(output_tree_path)]
    
    logging.info(f"Running MinVar rooting via FastRoot.py...")
    try:
        with open(stdout_log, 'w') as f_out, open(stderr_log, 'w') as f_err:
            subprocess.run(
                command, 
                shell=False,  
                check=True,   
                stdout=f_out, 
                stderr=f_err
            )
    except FileNotFoundError:
        # Specific error if python3 or FastRoot.py itself isn't found after shutil.which check
        logging.error(f"Command execution failed: 'python3' or '{fastroot_path}' not found.")
        return tree # Fallback to original tree
    except subprocess.CalledProcessError as e:
        logging.error(f"MinVar rooting failed. FastRoot.py exited with error code {e.returncode}.")
        logging.error(f"Check logs for details: {stdout_log} and {stderr_log}")
        return tree # Fallback to original tree
    except Exception as e: # Catch any other unexpected errors during subprocess run
        logging.error(f"An unexpected error occurred during MinVar rooting: {e}")
        return tree

    # 5. Load the newly rooted tree using 
    logging.info("MinVar rooting successful. Loading rooted tree...")
    try:
        #rooted_tree = PhyloTree(open(output_tree_path))
        with open(output_tree_path, 'r') as f:
            rooted_tree = PhyloTree(f)
        # Re-apply the species naming function as it's a new tree object
        rooted_tree.set_species_naming_function(lambda node: node.name.split(species_delimiter)[0])
        return rooted_tree
    except Exception as e:
        logging.error(f"Failed to load the MinVar rooted tree from {output_tree_path}. Error: {e}")
        return tree # Fallback to original tree if loading fails


