"""Run the finite checks and compare each JSON result with its reference output."""

import difflib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
CHECKS = ("direct_blocks", "pairs", "seeded_pairs", "weighted_tables", "routing", "addresses")


def formatted(value):
    return json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True).splitlines()


def main():
    for name in CHECKS:
        script = ROOT / "checks" / f"check_{name}.py"
        reference = ROOT / "expected" / f"{name}.json"
        run = subprocess.run([sys.executable, "-I", "-B", str(script)],
                             cwd=ROOT, capture_output=True, encoding="utf-8")
        if run.returncode:
            sys.stderr.write(run.stderr or run.stdout)
            print(f"FAIL {name}: exit code {run.returncode}", file=sys.stderr)
            return 1
        actual = json.loads(run.stdout)
        expected = json.loads(reference.read_text(encoding="utf-8"))
        if actual != expected:
            diff = difflib.unified_diff(formatted(expected), formatted(actual),
                                        fromfile=reference.name, tofile="computed")
            print("\n".join(diff), file=sys.stderr)
            print(f"FAIL {name}: different JSON values", file=sys.stderr)
            return 1
        print(f"PASS {name}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
