"""Run finite checks and scalar Lean proofs; not a full-paper verifier."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
for script in ("houghton_labels.py", "gadgetcheck.py", "tagcheck.py", "accuracycheck.py", "scalarcheck.py"):
    subprocess.run([sys.executable, script], cwd=ROOT, check=True, timeout=120)
subprocess.run(["lean", "-M", "512", "-j", "1", "-T", "100000", "Resources.lean"],
               cwd=ROOT / "lean", check=True, timeout=120)
print("All specified finite and scalar checks passed; see README.md for exclusions.")
