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


def test_preanalysis():
    filename = "P53.faa.nw"
    t = og_delineation.load_tree_local(filename)

    ncbi_file = "e6.taxa.sqlite"
    ncbi = og_delineation.load_taxonomy('NCBI', ncbi_file)

    rooting = "Midpoint"

    tree_nw, sp_set, total_mems_in_tree, SPTOTAL, props = og_delineation.run_preanalysis(t, filename, ncbi, rooting)

    assert len(total_mems_in_tree) == 1092



def test_long_branches():
    filename = "P53.faa.nw"
    t = og_delineation.load_tree_local(filename)

    long_leaves = og_delineation.detect_long_branches(t)

    assert len(long_leaves) == 16


def test_annot_main_table():
    filename = "post_P53.nw"

    t = og_delineation.load_tree_local(filename)
    main_table = "result_emapper.emapper.annotations"

    og_delineation.annot_tree_main_table(t, main_table)

    print


