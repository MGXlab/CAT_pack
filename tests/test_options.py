import unittest
from pathlib import Path

from CAT_pack.config.options import BatOptions


class LogPathTests(unittest.TestCase):
    def test_log_path_is_derived_from_output_prefix(self):
        options = BatOptions(
            bins=Path("bins"),
            database=Path("db"),
            output_prefix=Path("sample"),
        )

        self.assertEqual(options.log_path, Path("sample.log"))

    def test_explicit_log_file_overrides_output_prefix(self):
        options = BatOptions(
            bins=Path("bins"),
            database=Path("db"),
            output_prefix=Path("sample"),
            log_file=Path("my_run.log"),
        )

        self.assertEqual(options.log_path, Path("my_run.log"))


if __name__ == "__main__":
    unittest.main()
