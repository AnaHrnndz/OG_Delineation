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
from ogd.tree_setup import _apply_rooting

DATABASE_PATH = '/data/databases/ETE_taxonomy/EggNOG6/e6.taxa.sqlite'

class TestSetUP(unittest.TestCase):
    def test_load_tree_local(self):
        sp_delim = '.'
        with NamedTemporaryFile() as f_tree:  
            f_tree.write(ds.nw1.encode('utf8'))
            f_tree.flush()
            tpath = Path(f_tree.name)
            t = _load_gene_tree(tpath, sp_delim)

        self.assertEqual(ds.nw2, t.write())

    # def test_load_reftree(self):
        # sp_delim = '.'
        # ncbi = NCBITaxa(dbfile=DATABASE_PATH, memory=True)
        # with NamedTemporaryFile() as f_tree:  
            # f_tree.write(ds.nw1.encode('utf8'))
            # f_tree.flush()
            # tpath = Path(f_tree.name)
            # gene_tree = _load_gene_tree(tpath, sp_delim)
            # rtree = _load_reftree(rtree_path=None, gene_tree=gene_tree, taxonomy_db=ncbi)            
            
        # self.assertEqual(ds.reftree, rtree.write())



    def test_counter_sp(self):
        sp_delim = '.'
        ncbi = NCBITaxa(dbfile=DATABASE_PATH, memory=True)
        with NamedTemporaryFile() as f_tree:  
            f_tree.write(ds.nw1.encode('utf8'))
            f_tree.flush()
            tpath = Path(f_tree.name)
            gene_tree = _load_gene_tree(tpath, sp_delim)
            rtree = _load_reftree(rtree_path=None, gene_tree=gene_tree, taxonomy_db=ncbi)
            taxocounter = _load_taxonomy_counter(rtree)

            expected = ds.counter

            self.assertDictEqual(expected, taxocounter)


    def test_midpoint(self):
        sp_delim = '.'
        with NamedTemporaryFile() as f_tree,  TemporaryDirectory() as temp_dir: 
            f_tree.write(ds.ur_p53.encode('utf8'))
            f_tree.flush()
            
            tpath = Path(f_tree.name)
    
            t = _load_gene_tree(tpath, sp_delim)
            mp_nw = _apply_rooting(t, 'Midpoint', temp_dir, sp_delim)

        self.assertEqual(ds.mp_p53, mp_nw.write())


    def test_minvar(self):
        sp_delim = '.'
        with NamedTemporaryFile() as f_tree,  TemporaryDirectory() as temp_dir: 
            f_tree.write(ds.ur_p53.encode('utf8'))
            f_tree.flush()
            
            tpath = Path(f_tree.name)
           
            t = _load_gene_tree(tpath, sp_delim)
            mv_nw = _apply_rooting(t, 'MinVar', Path(temp_dir), sp_delim)

        self.assertEqual(ds.mv_p53, mv_nw.write())



if __name__ == '__main__':
    unittest.main()






# def test_long_branches():
    # filename = "P53.faa.nw"
    # t = og_delineation.load_tree_local(filename)

    # long_leaves = og_delineation.detect_long_branches(t)

    # assert len(long_leaves) == 16


# def test_annot_main_table():
    # filename = "post_P53.nw"

    # t = og_delineation.load_tree_local(filename)
    # main_table = "result_emapper.emapper.annotations"

    # og_delineation.annot_tree_main_table(t, main_table)

    # print

# def test_dummy_run():
    # assert os.system("""../og_delineation.py --tree egg6_ogd_vs_egg6_oficial/original_trees/D6J8K.nw --user_taxonomy /data/databases/ETE_taxonomy/EggNOG6/e6.taxa.sqlite \
                    # --output_path ./tests/ --rooting MinVar --raw_alg egg6_ogd_vs_egg6_oficial/original_trees/D6J8K.faa.aln  --get_pairs --mode fast""") == 0
