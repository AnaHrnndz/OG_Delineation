import sys
sys.path.append("..")

import og_delineation


def test_load_tree_local():
    filename = "P53.faa.nw"
    t = og_delineation.load_tree_local(filename)

    assert len(t) == 1092


def test_load_taxonomy():
    ncbi_file = "e6.taxa.sqlite"
    ncbi = og_delineation.load_taxonomy('NCBI', ncbi_file)

    assert 'Eukaryota' == ncbi.get_taxid_translator(['2759'])[2759]


def test_load_reftree():
    filename = "P53.faa.nw"
    t = og_delineation.load_tree_local(filename)

    ncbi_file = "e6.taxa.sqlite"
    ncbi = og_delineation.load_taxonomy('NCBI', ncbi_file)

    reftree = og_delineation.load_reftree(t=t, taxonomy_db=ncbi)

    assert len(reftree) == 441