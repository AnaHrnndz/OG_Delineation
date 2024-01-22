import random
import subprocess
from collections import defaultdict
from ete4 import SeqGroup

"""
EMAPPER ANNOTATION
        - Eggnog-Mapper: annotate tree sequences with eggnog-mapper, user needs to provide alignment file
"""


def annotate_with_emapper(t, alg, path_out):
    """"
        Add to the leaves of tree t the information coming from emapper
        and return the annotated tree.
    """
    # In the main program, we are actually only interested in annotations
    # the pfam domains, but we have to take it all from emapper.
    path2raw = alg2rawfasta(alg, path_out)

    path2main_table = run_emapper(path2raw, path_out)
    path2pfam_table = run_hmm_mapper(path2raw, path_out)

    t = annot_tree_pfam_table(t, path2pfam_table, alg)

    t = annot_tree_main_table(t, path2main_table)

    # TODO: Check if we really do not need this (because it
    # is done in the main function already?):
    #   t, all_props = run_clean_properties(t)
    #   run_write_post_tree(t, name_tree, path_out, all_props)

    return t

def alg2rawfasta(alg, path_out):

    fasta = SeqGroup(alg)
    path2raw = path_out+'/total_raw_fasta.faa'
    raw_fasta = open(path2raw, 'w')
    for num, (name, aa, _) in enumerate(fasta):
        clean_aa = aa.replace('-','')
        raw_fasta.write('>'+name+'\n'+clean_aa+'\n')
    raw_fasta.close()

    return path2raw


def run_emapper(path2raw, path_out):

    """
        Run eggnog-mapper:

        emapper-2.1.12-f8b9fa5 / Expected eggNOG DB version: 5.0.2 / Installed eggNOG DB version: unknown /
        Diamond version found: diamond version 2.0.11 / MMseqs2 version found: 113e3212c137d026e297c7540e1fcd039f6812b1 /
        Compatible novel families DB version: 1.0.1

    """

    subprocess.run(("python /data/soft/eggnog-mapper_2.1.12/eggnog-mapper/emapper.py --sensmode fast  \
        --data_dir /data/soft/eggnog-mapper_2.1.12/eggnog-mapper/data  \
        -i %s -o %s --output_dir %s" %(path2raw, 'result_emapper', path_out)), shell = True)


    path2main_table = path_out+'/result_emapper.emapper.annotations'
    return  path2main_table



def run_hmm_mapper(path2raw, path_out):

    """
        Pfam annotation with hmm_mapper from eggnog-mapper scripts
    """

    subprocess.run(("python /data/soft/eggnog-mapper_2.1.12/eggnog-mapper/hmm_mapper.py \
        --cut_ga --clean_overlaps clans --usemem --num_servers 1 --num_workers 8 --cpu 8 \
        --dbtype hmmdb  -d /data/soft/eggnog-mapper_2.1.9/data/pfam/Pfam-A.hmm \
        --hmm_maxhits 0 --hmm_maxseqlen 60000 \
        --qtype seq -i %s -o %s --output_dir %s" %(path2raw, 'result_emapper', path_out)), shell = True)

    path2pfam_table = path_out+'/result_emapper.emapper.hmm_hits'

    return path2pfam_table



def annot_tree_pfam_table(post_tree, pfam_table, alg_fasta):

    """
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


def annot_tree_main_table(post_tree, main_table):

    """
        Annotate tree with main table from eggnog-mapper
    """

    seq2info = defaultdict(dict)
    with open(main_table) as fin:
        for line in fin:
            if not line.startswith('#'):
                info = line.strip().split('\t')
                seq_name = info[0]
                eggnog_ogs = info[4]
                for og in eggnog_ogs.split(','):
                    level = og.split('|')[0].split('@')[1]
                    if level in ['2759', '2', '2157'] :
                        basal_og = og.split('|')[0].split('@')[0]
                pref_name = info[8]
                kegg_pathway = info[12]

                seq2info[seq_name]['basal_og'] = basal_og
                seq2info[seq_name]['pref_name'] = pref_name
                seq2info[seq_name]['kegg_path'] = kegg_pathway


    for l in post_tree:
        if l.name in seq2info.keys():
            info_dict = seq2info[l.name]
            for i, val in info_dict.items():
                l.add_prop(i, val)

    for n in post_tree.traverse():
        if not n.is_leaf:
            random_seq_name = random.choice(list(n.leaf_names()))
            random_node = next(post_tree.search_nodes(name=random_seq_name))
            random_node_basal_og = random_node.props.get('basal_og', 'none@none@none')
            random_node_pref_name = random_node.props.get('pref_name', 'none@none@none')
            random_node_kegg_path = random_node.props.get('kegg_path', 'none@none@none')

            n.add_prop('basal_og', random_node_basal_og)
            n.add_prop('pref_name', random_node_pref_name)
            n.add_prop('kegg_path', random_node_kegg_path)

    return post_tree

