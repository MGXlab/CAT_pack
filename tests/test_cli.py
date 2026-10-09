import unittest
from unittest.mock import patch

from typer.testing import CliRunner

from CAT_pack7.cli import app


class CliTests(unittest.TestCase):
    def test_shared_flags_in_all_command_help(self):
        for command in ("prepare", "cat", "bat", "bins"):
            with self.subTest(command=command):
                result = CliRunner().invoke(app, [command, "--help"], color=False)
                self.assertEqual(result.exit_code, 0, result.output)
                for flag in ("--threads", "--verbose", "--quiet", "--debug"):
                    self.assertIn(flag, result.output)

    def test_annotation_commands_receive_execution_flags(self):
        for command, input_flag in (("cat", "-c"), ("bat", "-b"), ("bins", "-b")):
            with self.subTest(command=command), patch("CAT_pack7.cli.run_annotation_cli") as run:
                result = CliRunner().invoke(app, [
                    command, input_flag, "input.fna", "-d", "database",
                    "-n", "8", "--quiet", "--verbose", "--debug",
                    "--diamond-mode", "sensitive", "--block-size", "4",
                    "--index-chunks", "2", "--no-self-hits", "--path-to-diamond", "diamond.exe",
                    "--sensitivity", "6.5", "--split-memory-limit", "1G",
                    "--path-to-mmseqs", "mmseqs.exe",
                ])
                self.assertEqual(result.exit_code, 0, result.output)
                run.assert_called_once()
                options = run.call_args.args[0]
                self.assertEqual(options.threads, 8)
                self.assertTrue(options.quiet)
                self.assertTrue(options.verbose)
                self.assertTrue(options.debug)
                self.assertEqual(options.diamond.mode, "sensitive")
                self.assertEqual(options.diamond.block_size, 4)
                self.assertEqual(options.diamond.index_chunks, 2)
                self.assertTrue(options.diamond.no_self_hits)
                self.assertEqual(str(options.diamond.path_to_diamond), "diamond.exe")
                self.assertEqual(options.mmseqs.sensitivity, 6.5)
                self.assertEqual(options.mmseqs.split_memory_limit, "1G")
                self.assertEqual(str(options.mmseqs.executable), "mmseqs.exe")

    def test_aligner_group_defaults(self):
        from CAT_pack7.config.options import DiamondOptions, MMseqsOptions

        for command, input_flag in (("cat", "-c"), ("bat", "-b")):
            with self.subTest(command=command), patch("CAT_pack7.cli.run_annotation_cli") as run:
                result = CliRunner().invoke(app, [command, input_flag, "input.fna", "-d", "database"])
                self.assertEqual(result.exit_code, 0, result.output)
                options = run.call_args.args[0]
                self.assertEqual(options.diamond, DiamondOptions())
                self.assertEqual(options.mmseqs, MMseqsOptions())

    def test_group_validation_rejects_invalid_aligner_values(self):
        for flag, value in (("--index-chunks", "0"), ("--sensitivity", "8")):
            with self.subTest(flag=flag), patch("CAT_pack7.cli.run_annotation_cli") as run:
                result = CliRunner().invoke(app, ["cat", "-c", "input.fna", "-d", "database", flag, value])
                self.assertEqual(result.exit_code, 2, result.output)
                run.assert_not_called()

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
