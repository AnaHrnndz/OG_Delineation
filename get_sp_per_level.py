from ete3 import PhyloTree, NCBITaxa
from collections import defaultdict
import json
import sys


t = PhyloTree(sys.argv[1])

if sys.argv[2]:
    print('\t'+os.path.dirname(sys.argv[2)))
    ncbi = NCBITaxa(sys.argv[2)
else:
    ncbi = NCBITaxa()


ncbi.annotate_tree(t)

print(len(t))


level2sp_mem = defaultdict(set)
level2sp_num = defaultdict()

for l in t:
    lin = l.lineage
    lin_2 = ncbi.get_lineage(l.name)
    
    for tax in lin:
        level2sp_mem[tax].add(l.name)
    
for taxid, mems in level2sp_mem.items():
    level2sp_num[taxid] = len(mems)


with open('levels2sp.json', 'w') as fp:
    json.dump(level2sp_num, fp)

