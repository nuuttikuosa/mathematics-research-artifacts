"""Run all published checks, or one suite, and verify the reference outputs."""
import argparse
from difflib import unified_diff
from pathlib import Path
import os
import subprocess
import sys
from time import monotonic


ROOT = Path(__file__).resolve().parent
SUITES = {
    "anchors": ("local_checks/check_anchor_cuts.py", "expected/anchors.txt"),
    "local": ("local_checks/compare_local.py", "expected/local.txt"),
    "specialization": ("local_checks/check_specialization.py", "expected/specialization.txt"),
    "global": ("global_checks/run_global.py", "global_checks/expected.txt"),
}


def run(name):
    script, reference = SUITES[name]
    print(f"Running {name} ...", flush=True)
    started = monotonic()
    result = subprocess.run(
        [sys.executable, "-B", str(ROOT / script)], cwd=ROOT, capture_output=True,
        env={**os.environ, "PYTHONUTF8": "1", "PYTHONOPTIMIZE": ""},
        text=True, encoding="utf-8",
    )
    if result.returncode:
        sys.stderr.write(result.stdout + result.stderr)
        raise SystemExit(result.returncode)
    print(result.stdout, end="")
    actual = result.stdout.strip().splitlines()
    expected = (ROOT / reference).read_text(encoding="utf-8-sig").strip().splitlines()
    if name == "local" and actual and actual[0].startswith("Python "):
        # The interpreter/library banner is informative; only the results are data.
        actual, expected = actual[1:], expected[1:]
    if actual != expected:
        diff = unified_diff(expected, actual, fromfile=reference, tofile=script, lineterm="")
        raise SystemExit("\n".join(diff))
    print(f"Passed {name} ({monotonic() - started:.1f} s).", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=SUITES)
    selected = parser.parse_args().suite
    for name in ([selected] if selected else SUITES):
        run(name)


if __name__ == "__main__":
    main()
