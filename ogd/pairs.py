import itertools
from collections import defaultdict

####    GET ALL ORTHOLOGS PAIRS ####


def get_all_pairs(t, total_outliers):

    'Recovery seqs will not be included'

    total_pairs = set()
    strict_pairs = set()

    for n in t.traverse():
        
        if n.props.get('evoltype_2') != 'S':
            continue

        source_seqs = n.props.get('leaves_ch1') or []
        ortho_seqs = n.props.get('leaves_ch2') or []

        if not source_seqs or not ortho_seqs:
            continue

        # Precalcular taxid de cada secuencia solo una vez
        source_taxids = {s: n[s].props.get('taxid') for s in source_seqs}
        ortho_taxids = {o: n[o].props.get('taxid') for o in ortho_seqs}

        # Agrupar por taxid
        source_tax2mems = defaultdict(list)
        for s, taxid in source_taxids.items():
            source_tax2mems[taxid].append(s)

        orthologs_tax2mems = defaultdict(list)
        for o, taxid in ortho_taxids.items():
            orthologs_tax2mems[taxid].append(o)

        # Recorrer combinaciones
        for l_source in source_seqs:
            source_tax = source_taxids[l_source]
            _otype = "one-to-" if len(source_tax2mems[source_tax]) == 1 else "many-to-"

            for l_ortho in ortho_seqs:
                ortho_tax = ortho_taxids[l_ortho]

                if source_tax == ortho_tax:
                    continue  # mismo taxón, no interesa

                otype = _otype + ("one" if len(orthologs_tax2mems[ortho_tax]) == 1 else "many")
                pair_info = (l_source, l_ortho, otype, n.name)
                total_pairs.add(pair_info)

                if l_source not in total_outliers and l_ortho not in total_outliers:
                    strict_pairs.add(pair_info)

    

    mssg = f"""
    4. Get Pairs
        -Total pairs: {len(total_pairs)}
        -Strict pairs: {len(strict_pairs)} """
    print(mssg)
   
    
    return total_pairs, strict_pairs


def write_pairs_table(clean_pairs, strict_pairs,path_out, clean_name_tree):

    #name_fam = name_tree.split('.',1)[0]
    pairs_results = path_out+'/'+clean_name_tree+'.pairs.tsv'

    with open(pairs_results, 'w') as fout:
        for pair in clean_pairs:
            fout.write('\t'.join(list(pair))+'\n')

    strict_pairs_results = path_out+'/'+clean_name_tree+'.stric_pairs.tsv'

    with open(strict_pairs_results, 'w') as fout:
        for pair in strict_pairs:
            fout.write('\t'.join(list(pair))+'\n')