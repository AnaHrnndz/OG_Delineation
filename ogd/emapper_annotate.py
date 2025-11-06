"""
Step 8 (Optional): Functional Annotation with eggNOG-Mapper.

This module provides functions to:
1. Run emapper.py and hmm_mapper.py on the input sequences.
2. Annotate the phylogenetic tree with the functional results.
3. Provide an alternative annotation method using 'treeprofiler'.
"""

import glob
import logging
import os
import random
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from ete4 import PhyloTree, SeqGroup, Tree

import ogd.utils as utils

# --- Main Annotation Orchestrator ---

def annotate_with_emapper(
    tree: PhyloTree, 
    alignment_path: str, 
    tmpdir: Path, 
    emapper_dmnd_db: str, 
    emapper_pfam_db: str
) -> PhyloTree:
    """
    Runs the full emapper annotation pipeline and adds results to the tree.

    Args:
        tree: The PhyloTree object to annotate.
        alignment_path: Path to the input sequence alignment file.
        tmpdir: Path to the temporary directory.
        emapper_dmnd_db: Path to the emapper diamond database (.dmnd file).
        emapper_pfam_db: Path to the emapper PFAM database (e.g., Pfam-A.hmm).

    Returns:
        The annotated PhyloTree object.
    """
    logging.info("\n--- Step 8.1: Running eggNOG-Mapper Annotation ---")
    
    # 1. Convert alignment to ungapped fasta
    raw_fasta_path = _alignment_to_raw_fasta(alignment_path, tmpdir)

    # 2. Run emapper.py (for main annotations)
    main_table_path = _run_emapper(raw_fasta_path, tmpdir, emapper_dmnd_db)

    # 3. Run hmm_mapper.py (for Pfam domains)
    pfam_table_path = _run_hmm_mapper(raw_fasta_path, tmpdir, emapper_pfam_db)

    # 4. Annotate the tree with the results
    logging.info("Annotating tree with emapper main results...")
    tree = _annot_tree_main_table(tree, main_table_path)
    logging.info("Annotating tree with Pfam domain results...")
    tree = _annot_tree_pfam_table(tree, pfam_table_path, alignment_path)
    
    return tree

# --- Emapper Execution Functions ---

def _alignment_to_raw_fasta(alignment_path: str, tmpdir: Path) -> Path:
    """
    Converts a gapped alignment file to an ungapped (raw) fasta file.
    """
    fasta = SeqGroup(str(alignment_path))
    raw_fasta_path = tmpdir / 'total_raw_fasta.faa'
    
    try:
        with open(raw_fasta_path, 'w') as raw_fasta:
            for name, seq, _ in fasta:
                clean_seq = seq.replace('-', '')
                raw_fasta.write(f">{name}\n{clean_seq}\n")
    except IOError as e:
        logging.error(f"Failed to write raw fasta file: {e}")
        raise
        
    return raw_fasta_path

def _run_emapper(raw_fasta_path: Path, tmpdir: Path, emapper_dmnd_db: str) -> Path:
    """
    Runs the main emapper.py script.
    """
    logging.info("Running emapper.py...")
    # Get the data directory (parent of the .dmnd file)
    data_dir = Path(emapper_dmnd_db).parent.resolve()
    output_file_base = "result_emapper"
    
    command = [
        "emapper.py",
        "--sensmode", "fast",
        "--cpu", "8",
        "--data_dir", str(data_dir),
        "--temp_dir", str(tmpdir),
        "-i", str(raw_fasta_path),
        "-o", output_file_base,
        "--output_dir", str(tmpdir)
    ]

    try:
        subprocess.run(command, shell=False, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        logging.error("emapper.py not found. Is it in your system's PATH?")
        raise
    except subprocess.CalledProcessError as e:
        logging.error(f"emapper.py failed with exit code {e.returncode}:")
        logging.error(f"STDOUT: {e.stdout}")
        logging.error(f"STDERR: {e.stderr}")
        raise

    output_table = tmpdir / f"{output_file_base}.emapper.annotations"
    if not output_table.exists():
        raise FileNotFoundError(f"emapper.py ran but output file not found: {output_table}")
        
    return output_table

def _run_hmm_mapper(raw_fasta_path: Path, tmpdir: Path, emapper_pfam_db: str) -> Path:
    """
    Runs the hmm_mapper.py script for Pfam annotations.
    """
    logging.info("Running hmm_mapper.py for Pfam domains...")
    output_file_base = "result_emapper"
    
    # Robustly find the data directory and db name
    db_path = Path(emapper_pfam_db)
    data_dir = db_path.parent.parent
    db_name = db_path.name

    command = [
        "hmm_mapper.py",
        "--cut_ga", "--clean_overlaps", "clans", "--usemem",
        "--num_servers", "1", "--num_workers", "4", "--cpu", "4",
        "--dbtype", "hmmdb",
        "--data_dir", str(data_dir),
        "-d", str(db_path),
        "--hmm_maxhits", "0", "--hmm_maxseqlen", "60000",
        "--qtype", "seq",
        "-i", str(raw_fasta_path),
        "-o", output_file_base,
        "--output_dir", str(tmpdir)
    ]
    
    try:
        subprocess.run(command, shell=False, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        logging.error("hmm_mapper.py not found. Is it in your system's PATH?")
        raise
    except subprocess.CalledProcessError as e:
        logging.error(f"hmm_mapper.py failed with exit code {e.returncode}:")
        logging.error(f"STDERR: {e.stderr}")
        raise

    output_table = tmpdir / f"{output_file_base}.emapper.hmm_hits"
    if not output_table.exists():
        raise FileNotFoundError(f"hmm_mapper.py ran but output file not found: {output_table}")
        
    return output_table

# --- Annotation Functions (Internal Logic) ---

def _annot_tree_pfam_table(tree: PhyloTree, pfam_table_path: Path, alignment_path: str) -> PhyloTree:
    """
    Annotates the tree with Pfam domains, translating coordinates to the alignment.
    """
    
    # 1. Build a map from raw (ungapped) sequence position to alignment (gapped) position
    fasta = SeqGroup(str(alignment_path))
    raw2alg_map = defaultdict(dict)
    for name, seq, _ in fasta:
        p_raw = 1
        for p_alg, char in enumerate(seq, 1):
            if char != '-':
                raw2alg_map[name][p_raw] = p_alg
                p_raw += 1

    # 2. Parse domain hits and store them by sequence
    seq2doms = defaultdict(list)
    tree_leaves = set(tree.leaf_names())
    with open(pfam_table_path) as f_in:
        for line in f_in:
            if line.startswith('#'):
                continue
            
            info = line.strip().split('\t')
            seq_name = info[0]
            
            # Ensure this sequence is actually in our tree
            if seq_name not in tree_leaves:
                continue
                
            dom_name = info[1]
            try:
                dom_start_raw = int(info[7])
                dom_end_raw = int(info[8])

                # Translate coordinates from raw to alignment space
                trans_dom_start = raw2alg_map[seq_name].get(dom_start_raw, dom_start_raw)
                trans_dom_end = raw2alg_map[seq_name].get(dom_end_raw, dom_end_raw)
                
                dom_info_string = f"{dom_name}@{trans_dom_start}@{trans_dom_end}"
                seq2doms[seq_name].append(dom_info_string)
            except (ValueError, KeyError, IndexError) as e:
                logging.warning(f"Could not parse Pfam hit for {seq_name}: {line.strip()}. Error: {e}")

    # 3. Annotate leaf nodes
    for leaf in tree:
        if leaf.name in seq2doms:
            domains_string = '|'.join(seq2doms[leaf.name])
            leaf.add_prop('dom_arq', domains_string)

    # 4. Annotate internal nodes (using user's random sampling logic)
    for node in tree.traverse():
        if not node.is_leaf:
            random_seq_name = random.choice(list(node.leaf_names()))
            random_leaf_node = tree[random_seq_name] # More efficient than search_nodes
            random_node_domains = random_leaf_node.props.get('dom_arq', 'none@none@none')
            node.add_prop('dom_arq', random_node_domains)

    return tree


def _annot_tree_main_table(tree: PhyloTree, main_table_path: Path) -> PhyloTree:
    """
    Annotates the tree with main emapper results (OGs, KOs, pathways).
    """

    def calculate_best_terms(seq_ids, annotation_dicts, annot_type):
        """Inner function to find the most common annotation among descendants."""
        term_counter = Counter()
        for seq_id in seq_ids:
            annotations_for_seq = annotation_dicts[annot_type].get(seq_id, [None])
            term_counter.update(annotations_for_seq)
    
        if None in term_counter:
             del term_counter[None]
        
        if term_counter:
            term, count = term_counter.most_common(1)[0]
            percentage = round(((count / len(seq_ids)) * 100), 3)
            return f"{term}|{percentage}"
        else:
            return None

    # 1. Parse the emapper annotation table into a dictionary
    seq2info = defaultdict(lambda: defaultdict(list))
    with open(main_table_path) as fin:
        for line in fin:
            if line.startswith('#'):
                continue

            info = line.strip().split('\t')
            seq_name = info[0]
            
            basal_og = None # Fix for potential UnboundLocalError
            eggnog_ogs = info[4]
            for og in eggnog_ogs.split(','):
                try:
                    level = og.split('|')[0].split('@')[1]
                    if level in ['2759', '2', '2157']:
                        basal_og = og.split('|')[0].split('@')[0]
                except IndexError:
                    continue # Skip malformed OG strings
            
            # Store annotations
            seq2info['basal_og'][seq_name] = [basal_og]
            seq2info['pref_name'][seq_name] = [info[8] or None] # Use None if string is empty
            seq2info['kegg_path'][seq_name] = [p for p in info[12].split(',') if p] or [None]
            seq2info['kegg_ko'][seq_name] = [k for k in info[11].split(',') if k] or [None]

    # 2. Annotate the tree
    for node in tree.traverse():
        leaf_names = list(node.leaf_names())
        
        # Add direct annotation to leaf nodes
        if node.is_leaf:
            for annot_type in seq2info.keys():
                lannot = seq2info[annot_type].get(node.name, [None])
                node.add_prop(annot_type, lannot)
        
        # Add consensus annotation to *all* nodes (overwrites leaf-specific lists)
        kko_top_term = calculate_best_terms(leaf_names, seq2info, "kegg_ko")
        kpath_top_term = calculate_best_terms(leaf_names, seq2info, "kegg_path")
        pname_top_term = calculate_best_terms(leaf_names, seq2info, "pref_name")
        basal_og_top_term = calculate_best_terms(leaf_names, seq2info, "basal_og")

        if kko_top_term: node.add_prop('kegg_ko', kko_top_term)
        if kpath_top_term: node.add_prop('kegg_path', kpath_top_term)
        if pname_top_term: node.add_prop('pref_name', pname_top_term)
        if basal_og_top_term: node.add_prop('basal_og', basal_og_top_term)
        
    return tree

# --- Alternative Annotation Method ---

def annot_treeprofiler(
    tree: PhyloTree, 
    alignment_path: str, 
    main_table_path: str, 
    pfam_table_path: str, 
    tmpdir: Path
) -> Tree:
    """
    Annotates a tree using the external 'treeprofiler' tool.
    
    This function is separate from the main emapper pipeline above and is
    called directly by ogd_core.py if pre-computed emapper files are provided.
    """
    logging.info("--- Step 8.2: Annotating with TreeProfiler ---")
    
    # 1. Get all property keys and write a temporary tree
    tree, all_props = utils.sanitize_tree_properties(tree)
    tmp_tree_name = "tmp_treeprofiler_in"
    tmp_tree_path = tmpdir / f"{tmp_tree_name}.tree_annot.nw"
    
    # Use the refactored utility function
    utils.write_annotated_tree(tree, tmpdir, tmp_tree_name, all_props)
    
    # 2. Define expected output path
    output_dir = tmpdir / "treeprofiler_out"
    output_dir.mkdir(exist_ok=True)
    expected_output_tree = output_dir / f"{tmp_tree_path.name}_annotated.nw"

    # 3. Run treeprofiler
    command = [
        "treeprofiler", "annotate",
        "--counter-stat", "none",
        "-t", str(tmp_tree_path),
        "-o", str(output_dir),
        "--alignment", str(alignment_path),
        "--emapper-pfam", str(pfam_table_path),
        "--emapper-annotations", str(main_table_path)
    ]
    
    try:
        subprocess.run(command, shell=False, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        logging.error("treeprofiler not found. Is it in your system's PATH?")
        raise
    except subprocess.CalledProcessError as e:
        logging.error(f"treeprofiler failed with exit code {e.returncode}:")
        logging.error(f"STDERR: {e.stderr}")
        raise
        
    # 4. Load the newly annotated tree
    if not expected_output_tree.exists():
        raise FileNotFoundError(f"treeprofiler ran but output file not found: {expected_output_tree}")
    
    # Using 'Tree' as in the original script, not 'PhyloTree'
    annotated_tree = Tree(str(expected_output_tree))
    
    return annotated_tree