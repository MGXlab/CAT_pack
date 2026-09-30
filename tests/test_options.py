import unittest
from pathlib import Path

from CAT_pack7.options import BatOptions


class OptionsTests(unittest.TestCase):
    def test_log_uses_output_prefix(self):
        options = BatOptions(bins=Path("bins"), database=Path("db"),
                             output_prefix=Path("sample"))
        self.assertEqual(options.log_path, Path("sample.log"))

    def test_custom_log_file(self):
        options = BatOptions(bins=Path("bins"), database=Path("db"),
                             log_file=Path("my_run.log"))
        self.assertEqual(options.log_path, Path("my_run.log"))


if __name__ == "__main__":
    unittest.main()
