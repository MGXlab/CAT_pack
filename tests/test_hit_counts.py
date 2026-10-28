import io
import unittest
from decimal import Decimal

import numpy as np

from CAT_pack.classification import ClassificationEngine
from CAT_pack.io.parsers import AlignmentInput
from CAT_pack.io.writers import ClassificationWriter
from CAT_pack.results import ClassificationStatus
from CAT_pack.tax import ORF, Taxonomy
from CAT_pack.utils.typing import EntityType


class SequenceHitCountTests(unittest.TestCase):
    def setUp(self):
        self.taxonomy = Taxonomy(
            parents=np.array([0, 0, 1, 1], dtype=np.int32),
            taxid_index={"1": 0, "2": 1, "10": 2, "20": 3},
            taxids=["1", "2", "10", "20"],
        )
        self.alignment = AlignmentInput(
            orf_names=["mixed", "unknown", "root"],
            orf_index={"mixed": 0, "unknown": 1, "root": 2},
            orf_offsets=np.array([0, 5, 7, 8], dtype=np.int64),
            hit_ids=np.array([0, 1, 1, 2, 3, 0, 2, 4], dtype=np.int32),
            bitscores=np.array(
                [2000, 1980, 1960, 1940, 1920, 1000, 990, 900], dtype=np.int64,
            ),
            hit_names=["unmapped_reference", "ref10", "missing_taxid", "ref20", "root_ref"],
        )
        self.hit_taxids = self.taxonomy.taxid_of_hits(
            self.alignment.hit_names,
            {"ref10": "10", "ref20": "20", "missing_taxid": "999", "root_ref": "1"},
        )
        self.entities = {
            "mixed": ["mixed"], "unknown": ["unknown"],
            "missing": ["absent"], "root": ["root"],
        }

    def classify(self, entity_type=EntityType.SEQUENCE):
        return ClassificationEngine(
            self.taxonomy, self.hit_taxids, Decimal("0.5"),
        ).classify(entity_type, self.entities, self.alignment)

    def render(self, result, orf_out=None):
        output = io.StringIO()
        writer = ClassificationWriter(
            output, orf_out, entity_type=result.entity_type,
            fraction=Decimal("0.5"), range_=Decimal("10"),
            branches={"1", "2"}, no_stars=True,
        )
        for entity in range(len(result)):
            writer.write(result, entity)
        return output.getvalue()

    def test_sequence_reason_counts_only_hits_with_resolved_taxids(self):
        result = self.classify()
        output = self.render(result)
        self.assertIn("mixed\ttaxid assigned\tbased on 3 hits\t", output)
        self.assertIn("root\ttaxid assigned\tbased on 1 hits\t", output)
        counts = dict(zip(result.orf_names, result.n_w_tax.tolist()))
        self.assertEqual(counts, {"absent": 0, "mixed": 3, "root": 1, "unknown": 0})
        mixed = result.orf_names.index("mixed")
        self.assertEqual(result.orf_n_hits[mixed], 5)
        self.assertEqual(result.orf_top_bitscore[mixed], 2000)
        self.assertEqual(result.orf_lca[mixed], self.taxonomy.taxid_index["2"])
        self.assertEqual(result.weighted.total[result.entities.index("mixed")], 2000)

    def test_unknown_taxids_and_missing_hits_remain_distinct(self):
        result = self.classify()
        self.assertEqual(
            result.get_status(result.entities.index("unknown")),
            ClassificationStatus.NO_TAXIDS,
        )
        self.assertEqual(
            result.get_status(result.entities.index("missing")),
            ClassificationStatus.NO_HITS,
        )
        for name, total in (("unknown", 2), ("absent", 0)):
            orf = result.orf_names.index(name)
            self.assertEqual(result.orf_n_hits[orf], total)
            self.assertEqual(result.n_w_tax[orf], 0)

    def test_empty_alignment_and_empty_hit_groups(self):
        for offsets in ([0], [0, 0, 0]):
            with self.subTest(offsets=offsets):
                lcas, counts = self.taxonomy.find_LCA_for_ORF(
                    np.array(offsets, dtype=np.int64),
                    np.array([], dtype=np.int32),
                )
                self.assertEqual(counts.tolist(), [0] * (len(offsets) - 1))
                self.assertEqual(lcas.tolist(), [ORF.NO_TAXID] * (len(offsets) - 1))

    def test_parallel_batches_preserve_counts_and_classifications(self):
        size = 40_000
        names = [f"q{i}" for i in range(size)]
        alignment = AlignmentInput(
            orf_names=names, orf_index={name: i for i, name in enumerate(names)},
            orf_offsets=np.arange(size + 1, dtype=np.int64),
            hit_ids=np.arange(size, dtype=np.int32) % 2,
            bitscores=np.full(size, 1000, dtype=np.int64),
            hit_names=["known", "unknown"],
        )
        hit_taxids = np.array([2, ORF.NO_TAXID], dtype=np.int32)
        entities = {name: [name] for name in names}
        serial = ClassificationEngine(
            self.taxonomy, hit_taxids, Decimal("0.5"), threads=1,
        ).classify(EntityType.SEQUENCE, entities, alignment)
        parallel = ClassificationEngine(
            self.taxonomy, hit_taxids, Decimal("0.5"), threads=2,
        ).classify(EntityType.SEQUENCE, entities, alignment)
        expected = [int(name[1:]) % 2 == 0 for name in parallel.orf_names]
        np.testing.assert_array_equal(parallel.n_w_tax, expected)
        for attribute in (
            "n_orfs_with_taxid", "orf_n_hits", "orf_lca", "status", "orf_top_bitscore",
        ):
            np.testing.assert_array_equal(
                getattr(serial, attribute), getattr(parallel, attribute),
            )
        self.assertEqual(parallel.n_assigned, size // 2)

    def test_cat_and_bat_orf_reports_keep_total_hit_counts(self):
        for entity_type in (EntityType.CONTIG, EntityType.BIN):
            with self.subTest(entity_type=entity_type):
                result = self.classify(entity_type)
                orfs = io.StringIO()
                output = self.render(result, orfs)
                self.assertIn("mixed\ttaxid assigned\tbased on 1/1 ORFs\t", output)
                prefix = "mixed\tmixed" if entity_type == EntityType.BIN else "mixed"
                self.assertIn(f"{prefix}\t5\t1;2\t200\n", orfs.getvalue())


if __name__ == "__main__":
    unittest.main()
