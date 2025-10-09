#!/usr/bin/env python3

from ete4.smartview import Layout, explorer
from ogd.import_layouts import all_layouts
import ogd.emapper_layouts as el
import sys
from ete4 import PhyloTree


def run_smartview(tree, alignment):

    props_popup = ['node_is_og', 'dist', 'species_losses', 'node_create_og', 'lca_node_name', 'len_leaves_in', 
    'taxid', 'sci_name', 'lineage', 'lca_dup', 'inparalogs_rate', 'ch1_name', 'dup_node_name', 'is_root', 
    'total_leaves', 'dups_up', 'ogs_up', 'common_name', 'dups_down', 'so_score_dup','ogs_down', 'score1', 
    'rank', 'dup_lineage', 'lca_node', 'len_leaves_out', 'species_losses_percentage', 'name', 'ch2_name', 
    'score2', 'sp_out', 'so_score', 'leaves_out','dup_score', 'overlap', 'evoltype_2', 'mOG', 'len_sp_in', 
    'best_tax', 'node_is_mog', 'recover_seqs', 'recover_in', 'Preferred_name', 'Preferred_name_counter',
    'eggNOG_OGs_counter', 'eggNOG_OGs', 'long_branch_outlier']

    
    tree.link_to_alignment(alignment=alignment, alg_format='fasta') 
    for l in tree:
        len_alg = len(l.props['sequence'])
        break

    
    pfam_layout = Layout('PFAM', draw_node=el.draw_pfam_domains(tree, len_alg=len_alg))
    all_layouts.append(pfam_layout)

   
    tree.explore(layouts=all_layouts, show_leaf_name=False, include_props=props_popup, keep_server=True, host='138.4.138.141', port=5000)


if __name__ == "__main__":
    
    
    tpath = sys.argv[1]
    apath = sys.argv[2]
    
    
    try:
        t = PhyloTree(tpath)
        
        
    except Exception as e:
        print(f"Error al cargar los archivos: {e}")
        sys.exit(1)

    
    run_smartview(t, apath)