from decimal import Decimal
from pathlib import Path
import sys
import tempfile
import unittest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "CAT_pack"))

import shared


def alignment_row(query, hit, bitscore):
    return "\t".join([query, hit] + ["0"] * 9 + [str(bitscore)])

class SharedFileTests(unittest.TestCase):
    def test_import_orfs_groups_headers_by_their_contig_prefix(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            protein_file = Path(temporary_directory) / "proteins.faa"
            protein_file.write_text(
                ">contig_CAT_1 # metadata\nCAT_PEPTIDE\n>contig_CAT_2\nCAT_BAT_PEPTIDE\n>contig_RAT_1\nRAT_PEPTIDE\n",
                encoding="utf-8",
            )

            contig2orfs = shared.import_ORFs(str(protein_file), log_file=None, quiet=True)

        self.assertEqual(
            contig2orfs,
            {"contig_CAT": ["contig_CAT_1", "contig_CAT_2"], "contig_RAT": ["contig_RAT_1"]},
        )

    def test_parse_tabular_alignment_keeps_only_hits_inside_the_bitscore_range(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            alignment_file = Path(temporary_directory) / "alignment.tsv"
            alignment_file.write_text(
                "\n".join(
                    [
                        alignment_row("orf_a", "hit_1", 100),
                        alignment_row("orf_a", "hit_2", 95),
                        alignment_row("orf_a", "hit_3", 94),
                        alignment_row("orf_a", "hit_4", 93),
                        alignment_row("orf_b", "hit_5", 50),
                    ]
                ) + "\n",
                encoding="utf-8",
            )

            orf2hits, all_hits = shared.parse_tabular_alignment(
                str(alignment_file), Decimal("0.95"), log_file=None, quiet=True
            )

        self.assertEqual(set(orf2hits), {"orf_a", "orf_b"})
        self.assertEqual({hit for hit, _ in orf2hits["orf_a"]}, {"hit_1", "hit_2"})
        self.assertEqual(all_hits, {"hit_1", "hit_2", "hit_5"})

if __name__ == "__main__":
    unittest.main()
