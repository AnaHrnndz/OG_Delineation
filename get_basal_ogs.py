import json
import sys
from collections import defaultdict
from ete4 import PhyloTree, NCBITaxa
import os



ncbi = NCBITaxa('data/taxonomy/e5.taxa.sqlite')



with open(sys.argv[1]) as f:
    og_dict = json.load(f)

target_level = sys.argv[2]


t = PhyloTree(sys.argv[3], format = 1)
t.set_species_naming_function(lambda x: x.split('.')[0])
ncbi.annotate_tree(t,  taxid_attr="species")

for name_og, info in og_dict.items():
    mems_og = set(info['mems'])
    og_dict[name_og]['anc_og'] = []
    node_intersection = []
    max_dist = t.get_distance(t,name_og, topology_only=True)
    og_dict[name_og]['dist'] = max_dist

    #save all og with shared members
    for n, i in og_dict.items():
        og2compare = set(i['mems'])
        if mems_og.intersection(og2compare):
            node_intersection.append(n)

    #get distance from node with shared memers to root
    save_dups = defaultdict()
    for node in node_intersection:
        root2node = t.get_distance(t,node, topology_only=True)
        if root2node < max_dist:
            save_dups[node] = root2node
    sort_dups = {k: v for k, v in sorted(save_dups.items(), key=lambda item: item[1] ,reverse = True)}
    prev_og = defaultdict(dict)

    if len(sort_dups) > 0:
        
        for og_anc, dist in sort_dups.items():
            og_anc_node = t.search_nodes(name=og_anc)[0]
            lca = og_anc_node.props.get('lca_node')
            prev_og[og_anc]['dist'] = dist
            prev_og[og_anc]['lca'] = lca

        og_dict[name_og]['anc_og'] = prev_og
        
    else:
        root_node = t.get_tree_root()
        root_lca = root_node.props.get('lca_node')
        prev_og[root_node.name]['dist'] = 0.0
        prev_og[root_node.name]['lca'] = root_lca
        og_dict[root_node.name]['anc_og'] = prev_og




save_dups = {}
og_level = {}
for og, val in og_dict.items():
    if og_dict[og]["dist"] != 0.0:
        if og_dict[og]["lca_dup"] == int(target_level):
            og_level[og] = [og_dict[og]["dist"], len(og_dict[og]["mems"])]

        elif int(target_level) in og_dict[og]["dup_lineage"]:
            save_dups[og] =  og_dict[og]["dist"]
        
sort_dups = {k: v for k, v in sorted(save_dups.items(), key=lambda item: item[1])}

print(og_level)
for og, dist in sort_dups.items():
    if not set(og_level.keys()).intersection(set(og_dict[og]["anc_og"].keys())):
        og_level[og] = [og_dict[og]["dist"], len(og_dict[og]["mems"])]


f_out_name = os.path.basename(sys.argv[1]).split('.')[0]+'_OGs_'+str(target_level)+'.json'
with open(f_out_name, 'w') as out:
    json.dump(og_level, out)

print(len(og_level))
mems_in_ogs = set()
total_mems = set()

for og in og_level.keys():
    mems_in_ogs.update(og_dict[og]["mems"])


for n in t.traverse():
    if n.is_leaf():
        lin = ncbi.get_lineage(n.props.get('taxid'))
        if int(target_level) in lin:
            total_mems.add(n.name)

diff = total_mems.difference(mems_in_ogs)
print(len(mems_in_ogs), len(total_mems), len(diff))




    


   
