import random
import subprocess
from collections import defaultdict, Counter
from ete4 import SeqGroup, PhyloTree, Tree
import glob
import ogd.utils as utils
import os
import pathlib


"""
EMAPPER ANNOTATION
        - Eggnog-Mapper: annotate tree sequences with eggnog-mapper, user needs to provide alignment file
"""


def annotate_with_emapper(t, alg, tmpdir, emapper_dmnd, emapper_pfam):
    
    """"
        Add to the leaves of tree t the information coming from emapper
        and return the annotated tree.
    """
   
    path2raw = alg2rawfasta(alg, tmpdir)

    path2main_table = run_emapper(path2raw, tmpdir, emapper_dmnd)
    path2pfam_table = run_hmm_mapper(path2raw, tmpdir, emapper_pfam)

    t = annot_tree_main_table(t, path2main_table)
    #t = annot_tree_pfam_table(t, path2pfam_table, alg)

    
    #t = annot_treeprofiler(t, alg, path2main_table, path2pfam_table, tmpdir)

    

    return t

def alg2rawfasta(alg, tmpdir):
    
    """
    Remove gaps '-' from alignment to get raw fasta 
    """

    fasta = SeqGroup(alg)
    path2raw = tmpdir+'total_raw_fasta.faa'
    raw_fasta = open(path2raw, 'w')
    for num, (name, aa, _) in enumerate(fasta):
        clean_aa = aa.replace('-','')
        raw_fasta.write('>'+name+'\n'+clean_aa+'\n')
    raw_fasta.close()

    return path2raw


def run_emapper(path2raw, tmpdir, emapper_dmnd):

    """
        Run eggnog-mapper:

        emapper-2.1.12-f8b9fa5 / Expected eggNOG DB version: 5.0.2 / Installed eggNOG DB version: unknown /
        Diamond version found: diamond version 2.0.11 / MMseqs2 version found: 113e3212c137d026e297c7540e1fcd039f6812b1 /
        Compatible novel families DB version: 1.0.1

    """

    data_path = pathlib.Path(emapper_dmnd).parent.resolve()
    #subprocess.run(("python /data/soft/eggnog-mapper_2.1.12/eggnog-mapper/emapper.py --sensmode fast 
    subprocess.run(("emapper.py --sensmode fast --cpu 8 \
        --data_dir %s --temp_dir %s -i %s -o %s --output_dir %s" %(data_path, tmpdir, path2raw, 'result_emapper', tmpdir)), shell = True)

    path2main_table = tmpdir+'result_emapper.emapper.annotations'
    return  path2main_table



def run_hmm_mapper(path2raw, tmpdir, emapper_pfam):

    """
        Pfam annotation with hmm_mapper from eggnog-mapper scripts
    """
    
    data_path = str(pathlib.Path(emapper_pfam).parent.resolve()).replace('pfam', '')
    
    #subprocess.run(("python /data/soft/eggnog-mapper_2.1.12/eggnog-mapper/hmm_mapper.py \
    subprocess.run(("hmm_mapper.py \
        --cut_ga --clean_overlaps clans --usemem --num_servers 1 --num_workers 4 --cpu 4 \
        --dbtype hmmdb  --data_dir %s -d %s \
        --hmm_maxhits 0 --hmm_maxseqlen 60000 \
        --qtype seq -i %s -o %s --output_dir %s" %(data_path, emapper_pfam, path2raw, 'result_emapper', tmpdir)), shell = True)

    path2pfam_table = tmpdir+'result_emapper.emapper.hmm_hits'

    return path2pfam_table



def annot_treeprofiler(t, aln, path2main_table, path2pfam_table, tmpdir ):


    t, all_props = utils.run_clean_properties(t)
    tmp_nw = 'tmp_tree'
    utils.run_write_post_tree(t, tmp_nw, tmpdir, all_props)
    path2tmp_nw = tmpdir+'tmp_tree.tree_annot.nw'
   
    subprocess.run(("treeprofiler annotate --counter-stat none \
                    -t %s  -o %s --alignment %s --emapper-pfam %s --emapper-annotations %s" %(path2tmp_nw, tmpdir , aln, path2pfam_table, path2main_table)), shell = True)

    #Open again the tree:
    path2tree_treprofiler = glob.glob(tmpdir+'*_annotated.nw')[0]
    
    t = Tree(open(path2tree_treprofiler), parser = 0)
    
    for n in t.traverse():
        if n.props.get('node_create_og', 'false') == 'True':
            print(n.props['lca_node'])
    return t
    




def annot_tree_pfam_table(post_tree, pfam_table, alg_fasta):

    """
        Deprecated
        Annotate tree with pfam table from eggnog-mapper
    """

    fasta = SeqGroup(alg_fasta)
    raw2alg = defaultdict(dict)
    for num, (name, seq, _) in enumerate(fasta):

        p_raw = 1
        for p_alg, (a) in enumerate(seq, 1):
            if a != '-':
                raw2alg[name][p_raw] = p_alg
                p_raw +=1

    seq2doms = defaultdict(list)
    with open(pfam_table) as f_in:
        for line in f_in:
            if not line.startswith('#'):
                info = line.strip().split('\t')
                if info [0] in post_tree.leaf_names():
                    seq_name = info[0]
                    dom_name = info[1]
                elif info[1] in post_tree.leaf_names():
                    seq_name = info[1]
                    dom_name = info[0]

                dom_start = int(info[7])
                dom_end = int(info[8])

                trans_dom_start = raw2alg[seq_name][dom_start]
                trans_dom_end = raw2alg[seq_name][dom_end]

                dom_info_string = '@'.join([dom_name, str(trans_dom_start), str(trans_dom_end)])
                seq2doms[seq_name].append(dom_info_string)



    for l in post_tree:
        if l.name in seq2doms.keys():
            domains = seq2doms[l.name]
            domains_string = '|'.join(domains)
            l.add_prop('dom_arq', domains_string)

    for n in post_tree.traverse():
        if not n.is_leaf:
            random_seq_name = random.choice(list(n.leaf_names()))
            random_node = next(post_tree.search_nodes(name=random_seq_name))
            random_node_domains = random_node.props.get('dom_arq', 'none@none@none')
            n.add_prop('dom_arq', random_node_domains)

    return post_tree


def annot_tree_main_table(t, main_table):

    """
        Deprecated
        Annotate tree with main table from eggnog-mapper
    """

    def calculate_best_terms(seq_ids, annotation_dicts, annot_type):
        
        
        term_counter = Counter()
        for seq_id in seq_ids:
            
            annotations_for_seq = annotation_dicts[annot_type].get(seq_id, [])
            
            term_counter.update(annotations_for_seq)

    
        if None in term_counter: # Check if None is actually a key
             del term_counter[None]
        
        if term_counter: # Check if term_counter is not empty
            term, count = term_counter.most_common(1)[0]
            
            percentage = round(((count / len(seq_ids)) * 100), 3)
            best_terms = '|'.join([term, str(percentage)])
       

        else:
            best_terms = None

       
        return best_terms



    seq2info = defaultdict(dict)
    with open(main_table) as fin:
        for line in fin:
            if line.startswith('#'):
                continue

            info = line.strip().split('\t')
            seq_name = info[0]
            eggnog_ogs = info[4]
            for og in eggnog_ogs.split(','):
                
                level = og.split('|')[0].split('@')[1]
                if level in ['2759', '2', '2157'] :
                    basal_og = og.split('|')[0].split('@')[0]
                    
            pref_name = info[8]
            kegg_pathway = info[12]
            kegg_ko = info[11]
            seq2info['basal_og'][seq_name] = [basal_og]
            seq2info['pref_name'][seq_name] = [pref_name]
            seq2info['kegg_path'][seq_name] = kegg_pathway.split(',')
            seq2info['kegg_ko'][seq_name] = kegg_ko.split(',')


    for n in t.traverse():
        
        if n.is_leaf:
            for annot_type in seq2info.keys():
                lannot = seq2info[annot_type].get(n.name)
                if lannot is not None:
                    n.add_prop(annot_type, lannot)

            
        
        leaves = list(n.leaf_names())

        kko_top_term = calculate_best_terms(leaves, seq2info, "kegg_ko")
        kpath_top_term = calculate_best_terms(leaves, seq2info, "kegg_path")
        pname_top_term = calculate_best_terms(leaves, seq2info, "pref_name")
        basal_og_top_term = calculate_best_terms(leaves, seq2info, "basal_og")

        if kko_top_term: 
            n.add_prop('kegg_ko', kko_top_term)

        if kpath_top_term: 
            n.add_prop('kegg_path', kpath_top_term)

        if pname_top_term: 
            n.add_prop('pref_name', pname_top_term)

        if basal_og_top_term: 
            n.add_prop('basal_og', basal_og_top_term)
 
        
        
    return t

