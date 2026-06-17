import os
import sys
import random
import itertools
import json
from tempfile import NamedTemporaryFile, TemporaryDirectory

import unittest
from pathlib import Path
from ete4 import PhyloTree, NCBITaxa
from . import datasets  as ds


from ogd.handle_input  import _load_gene_tree, _load_reftree, _load_taxonomy_counter
from ogd.outliers_scores import detect_long_leaves_branches

DATABASE_PATH = '/data/databases/ETE_taxonomy/EggNOG6/e6.taxa.sqlite'