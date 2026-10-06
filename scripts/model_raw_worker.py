"""Compatibility entry point; the installable implementation lives in one place."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qcforever_model_workers.worker import main

if __name__ == "__main__":
    main()
