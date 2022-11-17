import sys 
from ete4 import Tree

t = Tree(sys.argv[1], format = 1)

for n in t.traverse():
    if n.props.get('node_is_og') == 'True'  and n.props.get('lca_dup') == '33208':

        leaves = n.get_leaves()
        refog = set()
        for l in leaves:
            rog = (l.props.get('RefOG'))
            if rog != '-':
                refog.add(rog)

        if len(refog) >1:
            print(n.name, len(leaves), refog)
 