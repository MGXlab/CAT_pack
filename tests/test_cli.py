import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from rich.text import Text
from typer.testing import CliRunner

from CAT_pack.cli import app
from CAT_pack.config.options import DiamondOptions, MMseqsOptions
from CAT_pack.progress import ProgressTracker
from CAT_pack.utils.logging import Status

ANNOTATION_COMMANDS = (("cat", "-c"), ("bat", "-b"), ("bins", "-b"))


class CliTests(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()

    def test_command_help_lists_execution_flags(self):
        for command in ("prepare", "cat", "bat", "bins"):
            for force_color in (False, True):
                with self.subTest(command=command, force_color=force_color), patch.multiple(
                    "typer.rich_utils", FORCE_TERMINAL=force_color,
                    COLOR_SYSTEM="standard" if force_color else None, MAX_WIDTH=80,
                ):
                    result = self.runner.invoke(app, [command, "--help"], color=False)

                    self.assertEqual(result.exit_code, 0, result.output)
                    self.assertEqual("\x1b[" in result.output, force_color)
                    # Rich can force ANSI styling in CI despite color=False.
                    help_text = Text.from_ansi(result.output).plain
                    for flag in ("--threads", "--verbose", "--quiet", "--debug"):
                        self.assertIn(flag, help_text)

    def test_annotation_commands_forward_execution_flags(self):
        for command, input_flag in ANNOTATION_COMMANDS:
            with self.subTest(command=command), patch("CAT_pack.cli.run_annotation_cli") as run_annotation:
                arguments = [
                    command, input_flag, "input.fna",
                    "-d", "database",
                    "-n", "8",
                    "--quiet", "--verbose", "--debug",
                ]

                result = self.runner.invoke(app, arguments)

                self.assertEqual(result.exit_code, 0, result.output)
                run_annotation.assert_called_once()
                options = run_annotation.call_args.args[0]
                self.assertEqual(options.threads, 8)
                self.assertTrue(options.quiet)
                self.assertTrue(options.verbose)
                self.assertTrue(options.debug)

    def test_annotation_commands_forward_diamond_options(self):
        for command, input_flag in ANNOTATION_COMMANDS:
            with self.subTest(command=command), patch("CAT_pack.cli.run_annotation_cli") as run_annotation:
                arguments = [
                    command, input_flag, "input.fna",
                    "-d", "database",
                    "--diamond-mode", "sensitive",
                    "--block-size", "4",
                    "--index-chunks", "2",
                    "--no-self-hits",
                    "--path-to-diamond", "diamond.exe",
                ]

                result = self.runner.invoke(app, arguments)

                self.assertEqual(result.exit_code, 0, result.output)
                run_annotation.assert_called_once()
                options = run_annotation.call_args.args[0]
                self.assertEqual(options.diamond.mode, "sensitive")
                self.assertEqual(options.diamond.block_size, 4)
                self.assertEqual(options.diamond.index_chunks, 2)
                self.assertTrue(options.diamond.no_self_hits)
                self.assertEqual(options.diamond.path_to_diamond, Path("diamond.exe"))

    def test_annotation_commands_forward_mmseqs_options(self):
        for command, input_flag in ANNOTATION_COMMANDS:
            with self.subTest(command=command), patch("CAT_pack.cli.run_annotation_cli") as run_annotation:
                arguments = [
                    command, input_flag, "input.fna",
                    "-d", "database",
                    "--aligner", "mmseqs2",
                    "--sensitivity", "6.5",
                    "--split-memory-limit", "1G",
                    "--path-to-mmseqs", "mmseqs.exe",
                ]

                result = self.runner.invoke(app, arguments)

                self.assertEqual(result.exit_code, 0, result.output)
                run_annotation.assert_called_once()
                options = run_annotation.call_args.args[0]
                self.assertEqual(options.aligner, "mmseqs2")
                self.assertEqual(options.mmseqs.sensitivity, 6.5)
                self.assertEqual(options.mmseqs.split_memory_limit, "1G")
                self.assertEqual(options.mmseqs.executable, Path("mmseqs.exe"))

    def test_annotation_commands_use_default_aligner_options(self):
        for command, input_flag in ANNOTATION_COMMANDS:
            with self.subTest(command=command), patch("CAT_pack.cli.run_annotation_cli") as run_annotation:
                arguments = [command, input_flag, "input.fna", "-d", "database"]

                result = self.runner.invoke(app, arguments)

                self.assertEqual(result.exit_code, 0, result.output)
                run_annotation.assert_called_once()
                options = run_annotation.call_args.args[0]
                self.assertEqual(options.diamond, DiamondOptions())
                self.assertEqual(options.mmseqs, MMseqsOptions())

    def test_zero_index_chunks_is_rejected_before_annotation(self):
        arguments = ["cat", "-c", "input.fna", "-d", "database", "--index-chunks", "0"]

        with patch("CAT_pack.cli.run_annotation_cli") as run_annotation:
            result = self.runner.invoke(app, arguments)

        self.assertEqual(result.exit_code, 2, result.output)
        run_annotation.assert_not_called()

    def test_mmseqs_sensitivity_above_the_limit_is_rejected(self):
        arguments = ["cat", "-c", "input.fna", "-d", "database", "--sensitivity", "8"]

        with patch("CAT_pack.cli.run_annotation_cli") as run_annotation:
            result = self.runner.invoke(app, arguments)

        self.assertEqual(result.exit_code, 2, result.output)
        run_annotation.assert_not_called()

    def test_bat_help_describes_the_command(self):
        result = self.runner.invoke(app, ["bat", "--help"])

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertIn("Run Bin Annotation Tool (BAT).", result.output)

    def test_bat_requires_inputs(self):
        result = self.runner.invoke(app, ["bat"])

        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("Missing option", result.output)

    def test_download_keeps_file_callbacks_with_step_tracker(self):
        def download(args, report, report_file, queue):
            self.assertIsInstance(report.__self__, ProgressTracker)
            path = Path("example.faa.gz")
            queue([path])
            report("Downloading", Status.RUNNING)
            report_file(path, 5, 10)
            report_file(path, 10, 10)
            report("Downloading", Status.COMPLETE, 1, 1)
            return SimpleNamespace(db_fasta=path, names=path, nodes=path, acc2tax=path)

        for quiet in (False, True):
            with self.subTest(quiet=quiet), patch("CAT_pack.cli.run_download", side_effect=download):
                args = ["download", "--db", "nr", "--out-dir", "example-db"]
                if quiet:
                    args.append("--quiet")
                result = self.runner.invoke(app, args)
                self.assertEqual(
                    result.exit_code, 0, f"{result.output}\n{result.exception!r}",
                )
                if quiet:
                    self.assertEqual(result.output, "")
                else:
                    self.assertIn("Download completed", result.output)


if __name__ == "__main__":
    unittest.main()
