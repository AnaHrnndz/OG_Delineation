"""
Step 6: Extract Orthologous Pairs.

This module traverses the tree to find all speciation nodes and extracts all
pairs of orthologous sequences derived from them. It generates two files:
- 'total_pairs.tsv': All pairs found.
- 'strict_pairs.tsv': A subset of pairs excluding any sequences that
  were flagged as outliers in previous steps.
"""

import logging
from collections import defaultdict
from pathlib import Path
from typing import Set, Tuple

from ete4 import PhyloTree

# --- Main Function ---

def get_all_pairs(tree: PhyloTree, total_outliers: Set[str]) -> Tuple[Set[Tuple], Set[Tuple]]:
    """
    Extracts all orthologous pairs from speciation nodes in the tree.

    An orthologous pair is defined as two sequences from different species
    that descend from a common speciation node.

    Args:
        tree: The fully annotated PhyloTree.
        total_outliers: A set of all leaf names (sequences) flagged as outliers.

    Returns:
        A tuple containing two sets:
        - total_pairs: All orthologous pairs found.
        - strict_pairs: Pairs where neither sequence is an outlier.
    """
    total_pairs = set()
    strict_pairs = set()

    for node in tree.traverse():
        # Only speciation nodes ('S') generate orthologs
        if node.props.get('evoltype_2') != 'S':
            continue

        # Get the retained leaves from each child branch
        # These properties were set in the outliers_scores step
        source_seqs = node.props.get('leaves_ch1', [])
        ortho_seqs = node.props.get('leaves_ch2', [])

        if not source_seqs or not ortho_seqs:
            continue

        # 1. Pre-calculate the taxid for each sequence
        # We use node[seq_name] to get the leaf node object from its name
        source_taxids = {s: tree[s].props.get('taxid') for s in source_seqs}
        ortho_taxids = {o: tree[o].props.get('taxid') for o in ortho_seqs}

        # 2. Pre-group sequences by taxid to determine orthology type
        source_tax2mems = defaultdict(list)
        for seq, taxid in source_taxids.items():
            source_tax2mems[taxid].append(seq)

        orthologs_tax2mems = defaultdict(list)
        for seq, taxid in ortho_taxids.items():
            orthologs_tax2mems[taxid].append(seq)

        # 3. Generate all pairwise combinations
        for l_source in source_seqs:
            source_tax = source_taxids[l_source]
            
            # Determine the first part of the orthology type (e.g., "one-to-")
            otype_prefix = "one-to-" if len(source_tax2mems[source_tax]) == 1 else "many-to-"

            for l_ortho in ortho_seqs:
                ortho_tax = ortho_taxids[l_ortho]

                # Orthologs must be from different species
                if source_tax == ortho_tax:
                    continue

                # Determine the full orthology type (e.g., "one-to-many")
                otype_suffix = "one" if len(orthologs_tax2mems[ortho_tax]) == 1 else "many"
                orthology_type = otype_prefix + otype_suffix
                
                # Store the pair information
                pair_info = (l_source, l_ortho, orthology_type, node.name)
                total_pairs.add(pair_info)

                # If neither sequence is an outlier, add to 'strict_pairs'
                if l_source not in total_outliers and l_ortho not in total_outliers:
                    strict_pairs.add(pair_info)

    
    logging.info(f"Total pairs found: {len(total_pairs)}")
    logging.info(f"Strict pairs (excluding outliers): {len(strict_pairs)}")
    
    return total_pairs, strict_pairs


def write_pairs_table(
    total_pairs: Set[Tuple], 
    strict_pairs: Set[Tuple], 
    output_path: Path, 
    tree_name_base: str
):
    """
    Writes the total and strict pairs sets to two separate TSV files.

    Args:
        total_pairs: The set of all orthologous pairs.
        strict_pairs: The set of strict orthologous pairs.
        output_path: The directory (Path object) to write files to.
        tree_name_base: The base name for the output files.
    """
    # --- Write total pairs file ---
    total_pairs_path = output_path / f"{tree_name_base}.pairs.tsv"
    with open(total_pairs_path, 'w') as fout:
        # Add a header for clarity
        fout.write("Seq1\tSeq2\tOrthologyType\tSpeciationNode\n")
        for pair in sorted(list(total_pairs)): # Sort for consistent output
            # .join works directly on tuples, no need for list()
            fout.write('\t'.join(pair) + '\n')

    # --- Write strict pairs file ---
    # !! Bug Fix: Corrected typo from '.stric_pairs.tsv' to '.strict_pairs.tsv'
    strict_pairs_path = output_path / f"{tree_name_base}.strict_pairs.tsv"
    with open(strict_pairs_path, 'w') as fout:
        # Add a header for clarity
        fout.write("Seq1\tSeq2\tOrthologyType\tSpeciationNode\n")
        for pair in sorted(list(strict_pairs)): # Sort for consistent output
            fout.write('\t'.join(pair) + '\n')