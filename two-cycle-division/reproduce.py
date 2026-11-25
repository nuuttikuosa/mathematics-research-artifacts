"""Run the finite checks and compare their output with the reference files."""

import difflib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
REFERENCES = (
    "direct_blocks.json", "pairs.json", "seeded_pairs.json",
    "weighted_tables.json", "routing.json", "addresses.json",
    "shared_length.txt", "exact_covers.txt", "gcd_scaling.txt", "mixed_axes.txt",
    "consecutive_divisors.txt", "normal_form_successor_tables.txt",
    "collision_inverse.txt", "split_quota_dp.txt", "transient_targets.txt", "rotations_and_size.txt",
    "second_quotient.txt", "lifts_and_cost.txt", "crt_blocks.txt",
    "factoring_queries.txt", "final_crosschecks.json",
)


def output_lines(text, extension):
    if extension == ".json":
        text = json.dumps(json.loads(text), sort_keys=True, indent=2, ensure_ascii=True)
    return text.splitlines()


def main():
    for filename in REFERENCES:
        reference = ROOT / "expected" / filename
        name = reference.stem
        script = ROOT / "checks" / f"check_{name}.py"
        run = subprocess.run([sys.executable, "-I", "-B", str(script)],
                             cwd=ROOT, capture_output=True, encoding="utf-8")
        if run.returncode:
            sys.stderr.write(run.stderr or run.stdout)
            print(f"FAIL {name}: exit code {run.returncode}", file=sys.stderr)
            return 1
        actual = output_lines(run.stdout, reference.suffix)
        expected = output_lines(reference.read_text(encoding="utf-8"), reference.suffix)
        if actual != expected:
            diff = difflib.unified_diff(expected, actual, fromfile=filename,
                                        tofile="computed", lineterm="")
            print("\n".join(diff), file=sys.stderr)
            print(f"FAIL {name}: output differs from {filename}", file=sys.stderr)
            return 1
        print(f"PASS {name}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
