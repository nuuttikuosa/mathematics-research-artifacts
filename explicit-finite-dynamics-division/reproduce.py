"""Run the five checks in a fresh directory and compare their complete reports."""

import json
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parent
CHECKS = (
    ("explicit_tables.py", [], "explicit_tables.txt"),
    ("cycle_reference.py", [], "exact_check.json"),
    ("three_plane_cycles.py", [], "three_plane_cycles.json"),
    ("check_checksum.py", ["--output", "check_checksum.json"], "check_checksum.json"),
    ("check_routed.py", ["--blind-decoder", "three_plane_cycles.py",
                         "--output", "check_routed.json"], "check_routed.json"),
)


def main():
    with TemporaryDirectory() as directory:
        work = Path(directory)
        for script, _, _ in CHECKS:
            shutil.copyfile(ROOT / script, work / script)
        for script, args, report in CHECKS:
            # Ignore Python environment flags (in particular PYTHONOPTIMIZE)
            # so the scripts' mathematical assertions are always enabled.
            run = subprocess.run([sys.executable, "-E", "-B", script, *args],
                                 cwd=work, stdout=subprocess.PIPE, text=True, check=True)
            expected = (ROOT / "expected" / report).read_text(encoding="utf-8")
            if report.endswith(".json"):
                actual = (work / report).read_text(encoding="utf-8")
                matches = json.loads(actual) == json.loads(expected)
            else:
                matches = run.stdout.splitlines() == expected.splitlines()
            if not matches:
                raise SystemExit(f"{script}: output differs from expected/{report}")
            print(f"PASS  {script} matches expected/{report}", flush=True)


if __name__ == "__main__":
    main()
