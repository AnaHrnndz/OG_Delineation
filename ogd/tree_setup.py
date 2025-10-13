import subprocess
from ete4 import PhyloTree
from typing import Any, Tuple, Set, Optional, Dict
import ogd.utils as utils
import os
import sys
import shutil 



def run_setup(t, taxonomy_db, reftree, clean_name_tree, args, tmpdir):
    """
    Preanalysis steps: rooting, taxonomy annotation, set up node properties, 
    and name internal nodes.
    """
    
    rooting = args.rooting
    sp_delimitator = args.sp_delim
    

    # 1. Rooting
    t = run_rooting(t, rooting, tmpdir, sp_delimitator)

    # 2. Add Taxonomy Annotation (Do this first, as rooting/scores might need it)
    t = add_taxomical_annotation(t, taxonomy_db)

    # 3. Get original seqs and species set
    total_mems_in_tree = set(t.leaf_names())
    set_sp_total = t.get_species()
    num_total_sp = len(set_sp_total)

    # 4. Name internal nodes (using imported utils functions)
    for i, n in enumerate(t.traverse()):
        if not n.is_leaf:   
            # Assuming make_name and get_depth were moved to utils
            n.name = f"{utils.make_name(i)}-{utils.get_depth(n)}" 
    

    # 5. Logging and Output

    print(f"""        - Total leaves: {len(t)}
        - Total species: {num_total_sp}
    """)
    
    # 6. Clean properties and get Newick for later use
    t, props = utils.run_clean_properties(t)

    # Assuming get_newick now takes the tree and returns the string
    tree_nw = utils.get_newick(t, props) 

    return tree_nw, set_sp_total, total_mems_in_tree, num_total_sp


def run_rooting(t, rooting, tmpdir, sp_delimitator):
    """
    Applies tree rooting using Midpoint, MinVar, or no rooting.
    """

    if rooting == "Midpoint":
        try:
            root_mid = t.get_midpoint_outgroup()
            t.set_outgroup(root_mid)
            sys.stdout.write("        - Rooting: Midpoint successful.\n")
        except Exception as e:
            sys.stderr.write(f"        - WARNING: Error during Midpoint rooting: {e}\n")
            # Continue without rooting if it fails
            sys.stdout.write("        - No rooting applied.\n")

    elif rooting == "MinVar":
        sys.stdout.write("        - Rooting: Running MinVar...\n")
        t = run_minvar(t, tmpdir, sp_delimitator)
        
    else:
        sys.stdout.write("        - No rooting applied.\n")
    
    return t


def run_minvar(t, tmpdir, sp_delimitator):
    """
    Executes external FastRoot.py for MinVar rooting.
    """
    
    # 1. Locate FastRoot.py executable
    minvar_exec = shutil.which("FastRoot.py")
    if not minvar_exec:
        sys.stderr.write("        - ERROR: FastRoot.py executable not found in PATH. Skipping MinVar rooting.\n")
        return t

    # 2. Write the tree for MinVar input
    input_tree_minvar = utils.write_tree_for_minvar_rootin(t, tmpdir)
    path2tmptree = os.path.join(tmpdir, 'output_minvar_tree.nw')
    stdout_file = os.path.join(tmpdir, 'minvar.stdout')
    stderr_file = os.path.join(tmpdir, 'minvar.stderr')

    # 3. Run MinVar using the safer list-form for subprocess.run
    try:
        with open(stdout_file, 'w') as stdout_f, open(stderr_file, 'w') as stderr_f:
            subprocess.run(
                [sys.executable, minvar_exec, "-i", input_tree_minvar, "-o", path2tmptree],
                check=True,  # Raise error if MinVar fails
                stdout=stdout_f,
                stderr=stderr_f,
                text=True
            )
        
        # 4. Open the rooted tree
        t_minvar = PhyloTree(open(path2tmptree), parser=0)
        t_minvar.set_species_naming_function(lambda node: node.name.split(sp_delimitator)[0])
        ("        - MinVar rooting successful.\n")
        return t_minvar

    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"        - FATAL ERROR: MinVar rooting failed (Non-zero exit code). Check logs in {tmpdir}. {e}\n")
        # STOP EXECUTION
        sys.exit(1)

    except Exception as e:
        sys.stderr.write(f"        - FATAL ERROR: Failed to process MinVar output or load tree: {e}\n")
        # STOP EXECUTION
        sys.exit(1)


def add_taxomical_annotation(t, taxonomy_db):
    """
    Adds taxonomical annotation to nodes (sci_name, taxid, named_lineage, rank, etc.).
    """
    sys.stdout.write("        - Adding taxonomical annotations...\n")
    
    taxonomy_db.annotate_tree(t, taxid_attr="species") 
    
    return t

