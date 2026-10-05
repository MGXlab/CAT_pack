import unittest
from decimal import Decimal

from CAT_pack7.classification import ClassificationEngine, ORFStatus


class ClassificationTests(unittest.TestCase):
    def setUp(self):
        self.engine = ClassificationEngine(
            taxid2parent={"1": "1", "2": "1"},
            fastaid2taxid={"protein_a": "2"},
            fraction=Decimal("0.3"),
        )

    def test_orf_without_hits(self):
        result = self.engine.classify_orf("contig_1", [])
        self.assertEqual(result.status, ORFStatus.NO_HIT)

    def test_orf_with_one_hit(self):
        result = self.engine.classify_orf("contig_1", [("protein_a", Decimal(100))])
        self.assertEqual(result.status, ORFStatus.ASSIGNED)
        self.assertEqual(result.taxid, "2")


if __name__ == "__main__":
    unittest.main()
