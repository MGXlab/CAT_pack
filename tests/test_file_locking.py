import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from CAT_pack7.config.options import CatOptions
from CAT_pack7.pipeline import run_annotation
from CAT_pack7.utils.errors import InputError
from CAT_pack7.utils.locking import lock_outputs


class LockingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.output = self.root / "output"

    def test_if_another_process_is_blocked(self):
        code =  """
            import sys
            from pathlib import Path
            from CAT_pack7.utils.locking import lock_outputs
            from CAT_pack7.utils.errors import InputError
            
            try:
                with lock_outputs([Path(sys.argv[1])]):
                    pass
            except InputError:
                sys.exit(1)
            """
        command = [sys.executable, "-c", code, str(self.output)]
        with lock_outputs([self.output]):
            result = subprocess.run(command, capture_output=True, timeout=10)
            self.assertEqual(result.returncode, 1, result.stderr)
        result = subprocess.run(command, capture_output=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_exception_releases_the_lock(self):
        with self.assertRaises(RuntimeError):
            with lock_outputs([self.output]):
                raise RuntimeError("run failed")
        with lock_outputs([self.output]):
            pass

    def test_different_outputs_can_run_together(self):
        with lock_outputs([self.output]):
            with lock_outputs([self.root / "other"]):
                pass

    def test_conflicting_annotation_does_not_write_its_log(self):
        options = CatOptions(contigs=self.root / "contigs.fna", database=self.root / "db",
                             output_prefix=self.output)
        with lock_outputs([self.output]):
            with self.assertRaises(InputError):
                run_annotation(options, lambda *event: None)
        self.assertFalse(options.log_path.exists())


if __name__ == "__main__":
    unittest.main()
