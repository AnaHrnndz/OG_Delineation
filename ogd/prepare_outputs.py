import re
import csv
import pathlib
from collections import defaultdict



#####   FUNCIONTS TO PREPARE OUTPUTS FILES (newick, etc)    ####
def get_seq2og(ogs_info):

    seq2ogs = defaultdict(set)
    for og, info in ogs_info.items():
        for s in info['Mems']:
            seq2ogs[s].add(og+'@'+str(info['TaxoLevel']))
        for s in info['RecoverySeqs']:
            if s != '-':
                seq2ogs[s].add(og+'@'+str(info['TaxoLevel']))

    return seq2ogs



def write_seq2ogs(seq2ogs, path, name_tree):

    """
        Write a table with seq2ogs info
    """

    name_fam = name_tree.split('.')[0]
    seq2ogs_out = open(path+'/'+name_fam+'.seq2ogs.tsv', 'w')
    for seq, ogs in seq2ogs.items():
        seq2ogs_out.write(seq+'\t'+'|'.join(list(ogs))+'\n')

    seq2ogs_out.close()



def write_ogs_info(ogs_info, name_tree, path):

    name_fam = name_tree.split('.',1)[0]
    name_out =  path+'/'+name_fam+'.ogs_info.tsv'

    with open(name_out, "w",  newline='') as myfile:
        w = csv.writer(myfile, delimiter='\t')
        w.writerow(('#OG_name', 'TaxoLevel', 'SciName_TaxoLevel', 'AssocNode',  'NumSP', 'OG_down', 'OG_up', 
        'NumSeqs', 'NumRecoverySeqs', 'Lca_Dup', 'Species_Outliers', 'Num_SP_Outliers', 'Inparalogs_Rate', 'SP_overlap_dup',
        'Seqs', 'RecoverySeqs'))

        for og_name in ogs_info.keys():

            ogs_down = '|'.join(list(ogs_info[og_name]['OG_down']))
            ogs_up = '|'.join(list(ogs_info[og_name]['OG_up']))
            seqs = ','.join(ogs_info[og_name]['Mems'])
            if len(ogs_info[og_name]['RecoverySeqs']) == 0:
                rec_seqs = '-'
            else:
                rec_seqs = ','.join(ogs_info[og_name]['RecoverySeqs'])
            
            w.writerow((
                og_name,    #1
                ogs_info[og_name]['TaxoLevel'], #2
                ogs_info[og_name]['SciName_TaxoLevel'], #3
                ogs_info[og_name]['AssocNode'],   #4
                ogs_info[og_name]['NumSP'],   #5
                ogs_down, #6
                ogs_up ,  #7
                ogs_info[og_name]['NumMems'], #8
                ogs_info[og_name]['NumRecoverySeqs'], #9
                ogs_info[og_name]['Lca_Dup'],
                ogs_info[og_name]['Species_Outliers'],
                ogs_info[og_name]['Num_SP_Outliers'],
                ogs_info[og_name]['Inparalogs_Rate'],
                ogs_info[og_name]['SP_overlap_dup'],
                seqs, #10
                rec_seqs #11
            ))
