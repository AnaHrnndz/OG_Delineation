from ete4 import  PhyloTree, SeqGroup, Tree
from ete4 import NCBITaxa
import sys


taxid_ori_list = list()
with open(sys.argv[1]) as f_:
    for line in f_:
        taxid_ori_list.append(line.strip())

if sys.argv[2]:
    print('\t'+os.path.dirname(sys.argv[2)))
    ncbi = NCBITaxa(sys.argv[2)
else:
    ncbi = NCBITaxa()


print ('len taxid ncbi', len(taxid_ori_list))
tree = ncbi.get_topology(taxid_ori_list)
print('len tree ', len(tree))

tree.write(format = 1 , outfile = 'NCBITree.nw')

out_taxid = open('total_taxid.txt', 'w')
for leaf in tree:
    out_taxid.write(leaf.name+'\n')
out_taxid.close()
