"""
Tests for bug fixes identified during static analysis.

Each test targets a specific bug fix and verifies:
  - OS-2: sp_lost exits with sys.exit(1) when lca_node is absent from level2sp_mem
  - OS-3: outliers_detection skips lineage taxons missing from level2sp_mem (no crash)
  - OG-1: _process_root_ogs does not crash when lca_tree is not in its lineage set
"""

import sys
import unittest
from unittest.mock import MagicMock
from pathlib import Path

from ete4 import PhyloTree, NCBITaxa

from ogd.outliers_scores import sp_lost, outliers_detection
from ogd.orthologs_groups import _process_root_ogs

DATABASE_PATH = '/data/databases/ETE_taxonomy/EggNOG6/e6.taxa.sqlite'


class TestBugFixes(unittest.TestCase):

    # -------------------------------------------------------------------------
    # OS-2: sp_lost — KeyError when lca_node absent from level2sp_mem
    # -------------------------------------------------------------------------

    def test_sp_lost_exits_on_missing_lca(self):
        """
        OS-2: sp_lost should call sys.exit(1) when lca_node is absent from
        level2sp_mem. The lca_node should always be present (it is generated
        from the tree's own taxids), so its absence indicates bad input data.
        """
        tree = PhyloTree('(562.a:0.1,287.b:0.1);')
        tree.add_prop('lca_node', '2')       # Bacteria — intentionally absent from level2sp_mem
        tree.add_prop('sp_in', ['562', '287'])

        level2sp_mem = {}  # Empty — '2' missing
        sp_out = set()

        with self.assertRaises(SystemExit) as ctx:
            sp_lost(tree, sp_out, level2sp_mem)

        self.assertEqual(ctx.exception.code, 1)

    def test_sp_lost_empty_lca_no_exit(self):
        """
        OS-2 (complementary): sp_lost should NOT exit when lca_node is 'Empty'
        — this is the expected case for nodes where all species are outliers.
        """
        tree = PhyloTree('(562.a:0.1,287.b:0.1);')
        tree.add_prop('lca_node', 'Empty')
        tree.add_prop('sp_in', [])

        sp_out = {'562', '287'}
        level2sp_mem = {}  # Empty — irrelevant for 'Empty' lca_node

        sp_lost(tree, sp_out, level2sp_mem)  # Should not raise

        self.assertEqual(tree.props.get('species_losses_percentage'), 1.0)

    # -------------------------------------------------------------------------
    # OS-3: outliers_detection — KeyError when lineage taxon absent from counter
    # -------------------------------------------------------------------------

    def test_outliers_detection_missing_lineage_no_crash(self):
        """
        OS-3: outliers_detection should not crash when a lineage taxon of a
        candidate species is absent from level2sp_mem. This happens when the
        reference tree does not cover all taxa in the taxonomy DB.

        Setup:
          - 9 bacterial species + 1 eukaryote (Homo sapiens, 9606)
          - Bacteria cover 9/10 = 90% → best_tax = Bacteria
          - Human is a candidate to remove (not in best_tax)
          - Human's lineage contains Eukaryota (2759), Metazoa (33208), etc.
          - None of these eukaryote lineage taxons are in level2sp_mem
          → old code: KeyError on level2sp_mem[str(l)]
          → new code: continue, returns a set without crashing
        """
        ncbi = NCBITaxa(dbfile=DATABASE_PATH, memory=True)

        # 9 bacteria (different species) + 1 eukaryote
        nw = (
            '(((562.b1:0.1,287.b2:0.1):0.1,(1423.b3:0.1,1280.b4:0.1):0.1):0.1,'
            '((1718.b5:0.1,210.b6:0.1):0.1,(194.b7:0.1,1773.b8:0.1):0.1,'
            '(1313.b9:0.1,9606.euk:0.1):0.1):0.1);'
        )
        tree = PhyloTree(nw)
        tree.set_species_naming_function(lambda n: n.name.split('.')[0])
        tree.resolve_polytomy()
        ncbi.annotate_tree(tree, taxid_attr='species')
        content = tree.get_cached_content()

        # level2sp_mem only has the top bacterial taxons.
        # All eukaryote lineage taxons (2759, 33208, 7711, ...) are intentionally
        # absent — exactly the scenario that triggered the KeyError bug.
        bacteria_sps = {'562', '287', '1423', '1280', '1718', '210', '194', '1773', '1313'}
        all_sps = bacteria_sps | {'9606'}
        level2sp_mem = {
            '1':      all_sps,        # root of life
            '131567': all_sps,        # cellular organisms
            '2':      bacteria_sps,   # Bacteria — covers 9/10 = 90% of species
        }

        # Should not raise KeyError for any absent eukaryote lineage taxon
        result = outliers_detection(tree, 0.05, 0.9, content, level2sp_mem, set(), ncbi)

        self.assertIsInstance(result, set)

    # -------------------------------------------------------------------------
    # OG-1: _process_root_ogs — KeyError from set.remove() not caught by ValueError
    # -------------------------------------------------------------------------

    def test_process_root_ogs_empty_lca_no_crash(self):
        """
        OG-1: _process_root_ogs should not crash when lca_tree is 'Empty'.

        'Empty' is not in the exclusion list ['1','131567','r_root','Unk'], so
        the code enters the ancestral-OG block and tries to remove lca_tree from
        the lineage set. Since get_lineage('Empty') returns [], the set is empty
        and set.remove('Empty') raised KeyError — not caught by except ValueError.

        With the fix (discard), no exception is raised.
        """
        tree = PhyloTree('(562.a:0.1,287.b:0.1);')
        tree.name = 'A-1'
        tree.add_prop('lca_node', 'Empty')
        tree.add_prop('leaves_in', ['562.a', '287.b'])
        tree.add_prop('sp_in', ['562', '287'])
        tree.add_prop('sp_out', [])
        tree.add_prop('inparalogs_rate', 1.0)
        tree.add_prop('so_score', 0.0)

        ogs_info = {}
        mock_tax_db = MagicMock()

        # Old code: lineage.remove('Empty') → KeyError (not caught by except ValueError)
        # New code: lineage.discard('Empty') → no exception
        _process_root_ogs(tree, mock_tax_db, ogs_info)

        # get_lineage('Empty') returns [] → lineage is empty → no ancestral OGs created
        # monophyletic_og not set → root OG block also skipped
        self.assertEqual(ogs_info, {})


if __name__ == '__main__':
    unittest.main()
