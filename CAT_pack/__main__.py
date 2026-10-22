import sys
from pathlib import Path

root = str(Path(__file__).resolve().parent.parent)
if root not in sys.path:
    sys.path.insert(0, root)

from CAT_pack7.cli import app

if __name__ == "__main__":
    app()
