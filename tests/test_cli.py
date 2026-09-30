import unittest

from typer.testing import CliRunner

from CAT_pack7.cli import app


class CliTests(unittest.TestCase):
    def test_bat_help(self):
        result = CliRunner().invoke(app, ["bat", "--help"])
        self.assertEqual(result.exit_code, 0)
        self.assertIn("Run Bin Annotation Tool (BAT).", result.output)

    def test_bat_requires_inputs(self):
        result = CliRunner().invoke(app, ["bat"])
        self.assertEqual(result.exit_code, 2)
        self.assertIn("Missing option", result.output)


if __name__ == "__main__":
    unittest.main()
