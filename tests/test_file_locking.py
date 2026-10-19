import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from textwrap import dedent

from CAT_pack7.utils.locking import lock_outputs


class OutputLockingTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory(dir=Path(__file__).parent)
        self.addCleanup(self.temp_dir.cleanup)
        self.directory = Path(self.temp_dir.name)
        self.output_prefix = self.directory / "sample.CAT"

    def test_another_process_can_lock_only_after_release(self):
        child_code = dedent("""
            import sys
            from pathlib import Path

            from CAT_pack7.utils.locking import lock_outputs
            from CAT_pack7.utils.errors import InputError

            try:
                with lock_outputs([Path(sys.argv[1])]):
                    pass
            except InputError as error:
                print(error, file=sys.stderr)
                sys.exit(2)
        """)
        command = [sys.executable, "-B", "-c", child_code, str(self.output_prefix)]
        repository = Path(__file__).resolve().parents[1]

        with lock_outputs([self.output_prefix]):
            blocked = subprocess.run(
                command, cwd=repository, capture_output=True, text=True, timeout=10,
            )
            self.assertEqual(blocked.returncode, 2, blocked.stderr)
            self.assertIn("Another run is using", blocked.stderr)

        released = subprocess.run(
            command, cwd=repository, capture_output=True, text=True, timeout=10,
        )
        self.assertEqual(released.returncode, 0, released.stderr)

    def test_exception_releases_the_lock(self):
        with self.assertRaisesRegex(RuntimeError, "run failed"):
            with lock_outputs([self.output_prefix]):
                raise RuntimeError("run failed")

        with lock_outputs([self.output_prefix]):
            pass

    def test_different_outputs_can_run_together(self):
        other_prefix = self.directory / "other.CAT"

        with lock_outputs([self.output_prefix]):
            with lock_outputs([other_prefix]):
                pass


if __name__ == "__main__":
    unittest.main()
