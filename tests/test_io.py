from decimal import Decimal
from pathlib import Path
import os
import sys
import tempfile
import unittest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "CAT_pack"))

import shared
import bins


def alignment_row(query, hit, bitscore):
    return "\t".join([query, hit] + ["0"] * 9 + [str(bitscore)])


class BatFileTests(unittest.TestCase):
    def test_import_bins_reads_fasta_headers_and_ignores_non_bins(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            bin_directory = Path(temporary_directory)
            # I would suggest translating one of these sequences into amino acids
            # (Reading frame one) ;)
            (bin_directory / "bin_a.fna").write_text(
                ">contig_b description\nGATGCGCGCACCCATGTGGCGGATGAACGC\n"
                ">contig_a\nGATGCGCGCACCCATGTGGCGGATGAACGC\n", encoding="utf-8"
            )
            (bin_directory / "notes.txt").write_text("not a bin\n", encoding="utf-8")
            (bin_directory / ".hidden.fna").write_text(">hidden\nACGT\n", encoding="utf-8")

            bin2contigs, contig_names = bins.import_bins(
                str(bin_directory), ".fna", log_file=None, quiet=True
            )

        self.assertEqual(set(bin2contigs), {"bin_a.fna"})
        self.assertEqual(set(bin2contigs["bin_a.fna"]), {"contig_a", "contig_b"})
        self.assertEqual(contig_names, {"contig_a", "contig_b"})

    def test_import_bins_rejects_a_contig_in_two_bins(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            bin_directory = Path(temporary_directory)
            (bin_directory / "bin_a.fna").write_text(">shared_contig\nACGT\n", encoding="utf-8")
            (bin_directory / "bin_b.fna").write_text(">shared_contig\nACGT\n", encoding="utf-8")

            with self.assertRaises(SystemExit):
                bins.import_bins(str(bin_directory), ".fna", log_file=None, quiet=True)

    def test_extract_bin_orfs_sorts_contigs_and_skips_contigs_without_orfs(self):
        orfs = bins.extract_bin_orfs(
            ["contig_b", "contig_missing", "contig_a"],
            {"contig_a": ["contig_a_1"], "contig_b": ["contig_b_1", "contig_b_2"]},
        )

        self.assertEqual(orfs, ["contig_a_1", "contig_b_1", "contig_b_2"])

    def test_make_concatenated_fasta_preserves_sequences_and_normalizes_headers(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            bin_directory = Path(temporary_directory)
            (bin_directory / "bin_b.fna").write_text(">contig_b description\nTBB\n", encoding="utf-8")
            (bin_directory / "bin_a.fna").write_text(">contig_a another description\nMGX\n", encoding="utf-8")
            output_file = bin_directory / "concatenated.fna"

            bins.make_concatenated_fasta(
                str(output_file),
                {"bin_b.fna": ["contig_b"], "bin_a.fna": ["contig_a"]},
                str(bin_directory) + os.sep,
                log_file=None,
                quiet=True,
            )

            contents = output_file.read_text(encoding="utf-8")

        self.assertEqual(contents, ">contig_a\nTBB\n>contig_b\nMGX\n")


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
