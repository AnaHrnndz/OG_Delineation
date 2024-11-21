import itertools
from collections import defaultdict

####    GET ALL ORTHOLOGS PAIRS ####


def get_all_pairs(CONTENT):

    'Recovery seqs will not be included'

    total_pairs = set()
    total_outliers = set()
    strict_pairs = set()

    for n in CONTENT:

        total_outliers.update(n.props.get('leaves_out', set()))

        if n.props.get('evoltype_2') == 'S':
            
            source_seqs = n.props.get('leaves_ch1')
            ortho_seqs = n.props.get('leaves_ch2')
            
            source = defaultdict(list)
            for s_seq in source_seqs:
                taxid =  n[s_seq].props.get('taxid')
                source[taxid].append(s_seq)

            orthologs = defaultdict(list)
            for o_seq in ortho_seqs:
                taxid =  n[o_seq].props.get('taxid')
                orthologs[taxid].append(o_seq)
            
            for l_source in source_seqs:
                source_tax = n[l_source].props.get('taxid')
                
                if len(source[source_tax]) == 1:
                    _otype = "one-to-"
                else:
                    _otype = "many-to-"
                
                for l_ortho in ortho_seqs:
                    ortho_tax = n[l_ortho].props.get('taxid')
                    if len(orthologs[ortho_tax]) == 1:
                        otype = _otype + "one"
                    else:
                        otype = _otype + "many"
                    
                    if source_tax != ortho_tax:
                        tupla = (l_source, l_ortho, otype, n.name)
                        total_pairs.add(tupla)


    # remove from total_pairs, seqs that are outliers

    for p in total_pairs:
        if len(total_outliers.intersection(p)) == 0:
            strict_pairs.add(p)

    mssg = f"""
    4. Get Pairs
        -Total pairs: {len(total_pairs)} """
    print(mssg)
   
    
    return total_pairs, strict_pairs


def write_pairs_table(clean_pairs, strict_pairs,path_out, name_tree):

    name_fam = name_tree.split('.',1)[0]
    pairs_results = path_out+'/'+name_fam+'.pairs.tsv'

    with open(pairs_results, 'w') as fout:
        for pair in clean_pairs:
            fout.write('\t'.join(list(pair))+'\n')

    strict_pairs_results = path_out+'/'+name_fam+'.stric_pairs.tsv'

    with open(strict_pairs_results, 'w') as fout:
        for pair in strict_pairs:
            fout.write('\t'.join(list(pair))+'\n')