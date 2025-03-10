import re
import csv
import pathlib
import json
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
    seq2og_list =[]
    seq2ogs_out = open(path+'/'+name_fam+'.seq2ogs.tsv', 'w')
    for seq, ogs in seq2ogs.items():
        seq2og_list.append({seq:list(ogs)})
        seq2ogs_out.write(seq+'\t'+','.join(list(ogs))+'\n')

    seq2ogs_out.close()

    seq2ogs_json = (path+'/'+name_fam+'.seq2ogs.jsonl')
    with open(seq2ogs_json, "w") as file:
        for seq in seq2og_list:
            file.write(json.dumps(seq) + "\n")
    





def write_ogs_info(ogs_info, clean_name_tree, path):

    
    name_out =  path+'/'+clean_name_tree+'.ogs_info.tsv'

    with open(name_out, "w",  newline='') as myfile:
        w = csv.writer(myfile, delimiter='\t')
        w.writerow(('#OG_name', 'Lca_Dup','TaxoLevel', 'SciName_TaxoLevel', 'AssocNode',  'NumSP', 'OG_down', 'OG_up', 
        'NumSeqs', 'NumRecoverySeqs',  'Species_Outliers', 'Num_SP_Outliers', 'Inparalogs_Rate', 'SP_overlap_dup',
        'Seqs', 'RecoverySeqs'))

        for og_name in ogs_info.keys():
            
            og_name_extend = clean_name_tree+'@'+og_name
            ogs_down = ','.join(list(ogs_info[og_name]['OG_down']))
            ogs_up = ','.join(list(ogs_info[og_name]['OG_up']))
            sp_outliers = ','.join(ogs_info[og_name]['Species_Outliers'])
            seqs = ','.join(ogs_info[og_name]['Mems'])
            if len(ogs_info[og_name]['RecoverySeqs']) == 0:
                rec_seqs = '-'
            else:
                rec_seqs = ','.join(ogs_info[og_name]['RecoverySeqs'])
            
            w.writerow((
                og_name_extend,    #1
                ogs_info[og_name]['Lca_Dup'],
                ogs_info[og_name]['TaxoLevel'], #2
                ogs_info[og_name]['SciName_TaxoLevel'], #3
                ogs_info[og_name]['AssocNode'],   #4
                ogs_info[og_name]['NumSP'],   #5
                ogs_down, #6
                ogs_up ,  #7
                ogs_info[og_name]['NumMems'], #8
                ogs_info[og_name]['NumRecoverySeqs'], #9
                sp_outliers, #11
                ogs_info[og_name]['Num_SP_Outliers'], #12
                ogs_info[og_name]['Inparalogs_Rate'], #13
                ogs_info[og_name]['SP_overlap_dup'], #14
                seqs, #16
                rec_seqs #16
            ))
