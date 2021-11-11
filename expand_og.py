from ete3 import SeqGroup
import sys
import os
from collections import defaultdict


expand_table = sys.argv[1]
not_og_fasta = SeqGroup(sys.argv[2])
path2wd = os.path.dirname(os.path.abspath(expand_table))
ogs = defaultdict(dict)

with open(expand_table) as f:
    for line in f:
        info = line.split('\t')
        
        seq2add_name = info[0]
        seq2add_fasta = not_og_fasta.get_seq(seq2add_name)
        og_name = info[1].strip()

        if og_name in ogs.keys():
            ogs[og_name].set_seq(seq2add_name,seq2add_fasta)
            
        else:
            path2og = path2wd+'/'+og_name
            og2add = SeqGroup(path2og)
            ogs[og_name] = og2add

        
        
for k,val in ogs.items():
    path2write = path2wd+'/'+'final_'+k
    val.write(format='fasta', outfile=path2write)
    
    
        


