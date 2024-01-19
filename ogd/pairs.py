#!/usr/bin/env python3

import itertools



####    GET ALL ORTHOLOGS PAIRS ####


def get_all_pairs(CONTENT, total_mems_in_ogs):

    'Recovery seqs will not be included'

    def removeDuplicates(lst):
        return [t for t in (set(tuple(i) for i in lst))]


    total_pairs = set()
    for n in CONTENT:
        if n.props.get('evoltype_2') ==   'S':
            leaves0 = n[0].props.get('_leaves_in')
            leaves1 = n[1].props.get('_leaves_in')
            if leaves0 != None and leaves1 != None:
                total_pairs.update(itertools.product(leaves0, leaves1))


    clean_pairs = removeDuplicates(total_pairs)
    print('TOTAL PAIRS: ', len(clean_pairs))
    return clean_pairs

def write_pairs_table(clean_pairs, path_out, name_tree):

    name_fam = name_tree.split('.',1)[0]
    pairs_results = path_out+'/'+name_fam+'.pairs.tsv'

    with open(pairs_results, 'w') as fout:
        for pair in clean_pairs:
            fout.write('\t'.join(list(pair))+'\n')