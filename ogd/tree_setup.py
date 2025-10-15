# -*- coding: utf-8 -*-
"""
tree_setup.py - Pre-analysis steps for the gene tree.
Includes rooting, taxonomic annotation, and internal node naming.
"""

import subprocess
import sys
import shutil
from pathlib import Path
from typing import Tuple, Set, Any

from ete4 import PhyloTree
from ogd.new_utils import BaseTaxonomyHandler # --> Import the base handler for type hinting

# --> Import specific utils functions for clarity
from ogd.new_utils import make_name, get_depth, run_clean_properties, write_tree_for_minvar_rooting

def run_setup(
    t: PhyloTree, 
    tax_handler: BaseTaxonomyHandler, 
    args: Any, 
    tmpdir: Path
) -> Tuple[PhyloTree, Set[str], Set[str], int]:
    """
    Run all pre-analysis steps on the gene tree.

    Args:
        t (PhyloTree): The input gene tree object.
        tax_handler (BaseTaxonomyHandler): The initialized taxonomy handler.
        args (Any): The command-line arguments object.
        tmpdir (Path): The path to the temporary directory.

    Returns:
        A tuple containing:
        - The processed and annotated PhyloTree object.
        - A set of all species in the tree.
        - A set of all leaf names (members) in the tree.
        - The total number of unique species.
    """
    
    # 1. Root the tree based on the user's chosen method.
    t = run_rooting(t, args.rooting, tmpdir, args.sp_delim)

    # 2. Add taxonomic annotations to all nodes using the handler.
    t = add_taxomical_annotation(t, tax_handler)

    # 3. Get original sequence and species sets *after* processing.
    total_mems_in_tree = set(t.leaf_names())
    set_sp_total = t.get_species()
    num_total_sp = len(set_sp_total)

    # 4. Name all internal nodes for unique identification.
    for i, n in enumerate(t.traverse()):
        if not n.is_leaf:
            # Use utility functions to create a consistent, informative name.
            n.name = f"{make_name(i)}-{get_depth(n)}"
    
    # 5. Clean special characters from annotations.
    t, _ = run_clean_properties(t)

    print(f"""        - Total leaves: {len(t)}
        - Total species: {num_total_sp}""")
    
    # --> Return the PhyloTree object directly instead of converting to newick
    return t, set_sp_total, total_mems_in_tree, num_total_sp


def run_rooting(t: PhyloTree, rooting_method: str, tmpdir: Path, sp_delimitator: str) -> PhyloTree:
    """Applies the specified rooting method to the tree."""
    
    if rooting_method == "Midpoint":
        try:
            root_mid = t.get_midpoint_outgroup()
            if root_mid:
                t.set_outgroup(root_mid)
                print("        - Rooting: Midpoint successful.")
            else:
                print("        - WARNING: Midpoint rooting failed (could not find midpoint). No rooting applied.", file=sys.stderr)
        except Exception as e:
            print(f"        - WARNING: An error occurred during Midpoint rooting: {e}. No rooting applied.", file=sys.stderr)

    elif rooting_method == "MinVar":
        print("        - Rooting: Running MinVar...")
        t = run_minvar(t, tmpdir, sp_delimitator)
        
    else:
        print("        - No rooting applied.")
        
    return t


def run_minvar(t: PhyloTree, tmpdir: Path, sp_delimitator: str) -> PhyloTree:
    """Executes the external FastRoot.py script for MinVar rooting."""
    
    # 1. Locate the FastRoot.py executable in the system's PATH.
    minvar_exec = shutil.which("FastRoot.py")
    if not minvar_exec:
        print("        - ERROR: FastRoot.py executable not found in PATH. Skipping MinVar rooting.", file=sys.stderr)
        return t

    # 2. Prepare input and output file paths using pathlib.
    input_tree_path = write_tree_for_minvar_rooting(t, tmpdir)
    output_tree_path = tmpdir / 'output_minvar_tree.nw'
    stdout_log = tmpdir / 'minvar.stdout'
    stderr_log = tmpdir / 'minvar.stderr'

    # 3. Construct and run the command.
    command = [
        sys.executable, str(minvar_exec),
        "-i", str(input_tree_path),
        "-o", str(output_tree_path)
    ]
    
    try:
        with stdout_log.open('w') as stdout_f, stderr_log.open('w') as stderr_f:
            subprocess.run(
                command,
                check=True,       # Raise an exception if the command fails
                capture_output=False, # We are redirecting manually
                stdout=stdout_f,
                stderr=stderr_f
            )
        
        # 4. Load the newly rooted tree.
        t_minvar = PhyloTree(str(output_tree_path))
        t_minvar.set_species_naming_function(lambda node: node.name.split(sp_delimitator)[0])
        print("        - MinVar rooting successful.")
        return t_minvar

    except FileNotFoundError:
        print(f"        - ERROR: Could not find the MinVar output tree at '{output_tree_path}'.", file=sys.stderr)
        return t # Return the original tree
    except subprocess.CalledProcessError:
        print(f"        - FATAL ERROR: MinVar rooting failed. Check logs in '{tmpdir}'.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"        - FATAL ERROR: An unexpected error occurred during MinVar processing: {e}", file=sys.stderr)
        sys.exit(1)


def add_taxomical_annotation(t: PhyloTree, tax_handler: BaseTaxonomyHandler) -> PhyloTree:
    """
    Adds taxonomic annotations to tree nodes using the provided handler.
    
    Args:
        t (PhyloTree): The tree to be annotated.
        tax_handler (BaseTaxonomyHandler): The taxonomy handler (NCBI or GTDB).

    Returns:
        The annotated PhyloTree object.
    """
    print("        - Adding taxonomical annotations...")
    
    # --> The handler's internal db object is used for the annotation.
    # This keeps the core logic agnostic to the database type.
    tax_handler.db.annotate_tree(t, taxid_attr="species")
    
    return t