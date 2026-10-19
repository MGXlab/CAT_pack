import unittest
from decimal import Decimal

import numpy as np

from CAT_pack7.classification import ClassificationEngine
from CAT_pack7.io.parsers import AlignmentInput
from CAT_pack7.results import ClassificationStatus
from CAT_pack7.tax import ORF, Taxonomy


class BasicClassificationTests(unittest.TestCase):
    def test_missing_markers_do_not_hide_second_orf_or_third_taxon(self):
        taxonomy = Taxonomy(
            parents=np.array([0, 0, 1, 1], dtype=np.int32),
            taxid_index={"1": 0, "2": 1, "10": 2, "20": 3},
            taxids=["1", "2", "10", "20"],
        )
        alignment = AlignmentInput(
            orf_names=["a_1", "b_1"], orf_index={"a_1": 0, "b_1": 1},
            orf_offsets=np.array([0, 1, 2], dtype=np.int64),
            hit_ids=np.array([0, 1], dtype=np.int32),
            bitscores=np.array([1000, 1000], dtype=np.int64),
            hit_names=["ref_a", "ref_b"],
        )
        engine = ClassificationEngine(taxonomy, np.array([2, 3], dtype=np.int32), Decimal("0.5"))
        result = engine.classify("contig", {"a": ["a_1"], "b": ["b_1"]}, alignment)
        self.assertEqual(result.orf_lca.tolist(), [2, 3])
        self.assertEqual(result.n_assigned, 2)

    def setUp(self):
        self.taxonomy = Taxonomy(
            parents=np.array([0, 0], dtype=np.int32),
            taxid_index={"1": 0, "2": 1},
            taxids=["1", "2"],
        )

    def test_orf_without_hits_leaves_its_contig_unassigned(self):
        alignment = AlignmentInput(
            orf_names=[],
            orf_index={},
            orf_offsets=np.array([0], dtype=np.int64),
            hit_ids=np.array([], dtype=np.int32),
            bitscores=np.array([], dtype=np.int64),
            hit_names=[],
        )
        engine = ClassificationEngine(
            taxonomy=self.taxonomy,
            fastaid2taxid=np.array([], dtype=np.int32),
            fraction=Decimal("0.3"),
        )

        result = engine.classify("contig", {"contig": ["contig_1"]}, alignment)

        self.assertEqual(result.get_status(0), ClassificationStatus.NO_HITS)
        self.assertEqual(result.orf_lca.tolist(), [ORF.NO_TAXID])
        self.assertEqual(result.orf_n_hits.tolist(), [0])
        self.assertEqual(list(result.get_assignments(0)), [])

    def test_single_hit_assigns_its_reference_taxid(self):
        alignment = AlignmentInput(
            orf_names=["contig_1"],
            orf_index={"contig_1": 0},
            orf_offsets=np.array([0, 1], dtype=np.int64),
            hit_ids=np.array([0], dtype=np.int32),
            bitscores=np.array([1000], dtype=np.int64),  # A score of 100 in tenths.
            hit_names=["protein_a"],
        )
        engine = ClassificationEngine(
            taxonomy=self.taxonomy,
            fastaid2taxid=np.array([1], dtype=np.int32),  # protein_a maps to taxid 2.
            fraction=Decimal("0.3"),
        )

        result = engine.classify("contig", {"contig": ["contig_1"]}, alignment)

        self.assertEqual(result.get_status(0), ClassificationStatus.ASSIGNED)
        self.assertEqual(result.orf_n_hits.tolist(), [1])
        self.assertEqual(result.orf_top_bitscore.tolist(), [1000])
        self.assertEqual(list(result.get_assignments(0)), [0])
        assigned_taxon = result.weighted.assignment_taxon[0]
        self.assertEqual(result.taxonomy.taxids[assigned_taxon], "2")


if __name__ == "__main__":
    unittest.main()
