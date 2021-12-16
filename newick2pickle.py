from pickle import load, dump
import sys
from ete4 import PhyloTree
import os



t = PhyloTree(sys.argv[1]) 
name = os.path.basename(sys.argv[1])
name_no_extension = name.replace('.faa.aln.nw', '')
name_out = name_no_extension+'.pickle'

with open(name_out, "wb") as handle:
    dump(t, handle)

# with open("tree.pickle", "rb") as handle:
    # t = load(handle)