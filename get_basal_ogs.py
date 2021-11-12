import json
import sys
from collections import defaultdict
from ete3 import PhyloTree, NCBITaxa
import os



ncbi = NCBITaxa('/data/jhc/cold/eggnog6/build/00_level_clades/ref_trees_damian_taxonomy/etetoolkit/taxa.sqlite')



with open(sys.argv[1]) as f:
    data = json.load(f)


target_level = sys.argv[2]

tree = PhyloTree(sys.argv[3], format = 1)
tree.set_species_naming_function(lambda x: x.split('.')[0])
ncbi.annotate_tree(tree,  taxid_attr="species")

   
save_dups = {}
og_level = {}
for og, val in data.items():
    if data[og]["dist"] != 0.0:

        if data[og]["lca"] == int(target_level):
             og_level[og] = [data[og]["dist"], len(data[og]["mems"])]
    
        elif int(target_level) in data[og]["dup_lineage"]:
            save_dups[og] =  data[og]["dist"]
        
sort_dups = {k: v for k, v in sorted(save_dups.items(), key=lambda item: item[1])}

for og, dist in sort_dups.items():
    if not set(og_level.keys()).intersection(set(data[og]["anc_og"].keys())):
        og_level[og] = [data[og]["dist"], len(data[og]["mems"])]


f_out_name = os.path.basename(sys.argv[1]).split('.')[0]+'_OGs_'+str(target_level)+'.json'
with open(f_out_name, 'w') as out:
    json.dump(og_level, out)

print(len(og_level))
mems_in_ogs = set()
total_mems = set()

for og in og_level.keys():
    mems_in_ogs.update(data[og]["mems"])


for n in tree.traverse():
    if n.is_leaf():
        lin = ncbi.get_lineage(n.taxid)
        if int(target_level) in lin:
            total_mems.add(n.name)

diff = total_mems.difference(mems_in_ogs)
print(len(mems_in_ogs), len(total_mems), len(diff))




    


   
