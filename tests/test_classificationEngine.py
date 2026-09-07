from decimal import Decimal
from pathlib import Path
import sys
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "CAT_pack"))

import classification


TAXID_TO_PARENT = {
    "1": "1",
    "2": "1",
    "10": "2",
    "11": "10",
    "12": "10",
    "13": "11",
    "14": "11",
}

FASTAID_TO_TAXID = {"species_a": "13", "species_b": "14"}


class ClassificationEngineTests(unittest.TestCase):
    def setUp(self):
        self.engine = classification.ClassificationEngine(
            taxid2parent=TAXID_TO_PARENT,
            fastaid2taxid=FASTAID_TO_TAXID,
            fraction=Decimal("0.3"),
        )

    def test_classify_orf_without_hits_reports_no_hit(self):
        result = self.engine.classify_orf("contig_1", None)

        self.assertEqual(result.status, classification.ORFStatus.NO_HIT)
        self.assertEqual(result.orf_id, "contig_1")

if __name__ == "__main__":
    unittest.main()
