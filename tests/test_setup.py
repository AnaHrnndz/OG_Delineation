import os
import sys
import random
import itertools
import json
from tempfile import NamedTemporaryFile
import unittest
from pathlib import Path
from ete4 import PhyloTree
from . import datasets  as ds


from ogd.handle_input  import _load_gene_tree


class TestSetUP(unittest.TestCase):
    def test_load_tree_local(self):
        sp_delim = '.'
        with NamedTemporaryFile() as f_tree:  # test reading from file
            f_tree.write(ds.nw_p53.encode('utf8'))
            f_tree.flush()
            tpath = Path(f_tree.name)
            t = _load_gene_tree(tpath, sp_delim)

       
        self.assertEqual(ds.p_nwpost, t.write())



if __name__ == '__main__':
    unittest.main()



# def test_load_taxonomy():
    # ncbi_file = "e6.taxa.sqlite"
    # ncbi = og_delineation.load_taxonomy('NCBI', ncbi_file)

    # assert 'Eukaryota' == ncbi.get_taxid_translator(['2759'])[2759]


# def test_load_reftree():
    # filename = "P53.faa.nw"
    # t = og_delineation.load_tree_local(filename)

    # ncbi_file = "e6.taxa.sqlite"
    # ncbi = og_delineation.load_taxonomy('NCBI', ncbi_file)

    # reftree = og_delineation.load_reftree(t=t, taxonomy_db=ncbi)

    # assert len(reftree) == 441


# def test_preanalysis():
    # filename = "P53.faa.nw"
    # t = og_delineation.load_tree_local(filename)

    # ncbi_file = "e6.taxa.sqlite"
    # ncbi = og_delineation.load_taxonomy('NCBI', ncbi_file)

    # rooting = "Midpoint"

    # tree_nw, sp_set, total_mems_in_tree, SPTOTAL, props = og_delineation.run_preanalysis(t, filename, ncbi, rooting)

    # assert len(total_mems_in_tree) == 1092



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
