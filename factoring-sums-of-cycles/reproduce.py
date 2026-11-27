"""Run the finite checks in isolation and compare every deterministic result.

Requires Python 3.10+ and the standard library. No network access is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
SUITES = {
    "algebra": ("checks/algebra/check_algebra.py", "results.json", {"seconds"}),
    "sat": ("checks/sat/independent_routing_check.py", "results.json", {"seconds"}),
    "parameters": ("checks/parameters/parameter_check.py", "parameter_check.json", set()),
    "network": ("checks/network/network_geometry_check.py", "network_geometry_results.json", {"source_sha256"}),
    "blocks": ("checks/blocks/blocks_geometry_check.py", "blocks_geometry_results.json", set()),
    "mass": ("checks/mass/independent_mass_check.py", "independent_mass_results.json", set()),
    "associates": ("checks/associates/independent_associate_check.py", "independent_associate_results.json", {"source_sha256"}),
}


def semantic_result(result, metadata_fields, executed_source):
    """Validate run metadata before omitting it from cross-machine comparison."""
    if not isinstance(result, dict):
        raise ValueError("The checker did not return a JSON object")
    result = dict(result)
    for field in metadata_fields:
        if field not in result:
            raise ValueError(f"Missing expected run metadata: {field}")
        value = result.pop(field)
        if field == "seconds":
            if (isinstance(value, bool) or not isinstance(value, (int, float))
                    or not math.isfinite(value) or value < 0):
                raise ValueError("Invalid elapsed-time value")
        elif field == "source_sha256":
            expected_hash = hashlib.sha256(executed_source.read_bytes()).hexdigest()
            if value != expected_hash:
                raise ValueError("The reported source hash does not match the executed file")
    return result


def first_difference(expected, actual, path="$"):
    """Return a useful location without suppressing any result field."""
    if type(expected) is not type(actual):
        return f"{path}: expected type {type(expected).__name__}, got {type(actual).__name__}"
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            missing = sorted(expected.keys() - actual.keys())
            extra = sorted(actual.keys() - expected.keys())
            return f"{path}: missing keys {missing}; extra keys {extra}"
        for key in expected:
            difference = first_difference(expected[key], actual[key], f"{path}.{key}")
            if difference:
                return difference
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            return f"{path}: expected {len(expected)} items, got {len(actual)}"
        for index, (left, right) in enumerate(zip(expected, actual)):
            difference = first_difference(left, right, f"{path}[{index}]")
            if difference:
                return difference
    elif expected != actual:
        return f"{path}: expected {expected!r}, got {actual!r}"
    return None


def run_suite(name, output_dir=None, timeout=600):
    relative_source, output_name, metadata_fields = SUITES[name]
    source = ROOT / relative_source
    expected = json.loads((ROOT / "expected" / f"{name}.json").read_text(encoding="utf-8"))
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix=f"cycle-factor-{name}-") as directory:
        working = Path(directory)
        copied_source = working / source.name
        shutil.copyfile(source, copied_source)
        completed = subprocess.run(
            [sys.executable, "-I", "-B", str(copied_source)],
            cwd=working, capture_output=True, text=True, encoding="utf-8",
            timeout=timeout,
        )
        if completed.returncode:
            raise RuntimeError(
                f"{name}: checker exited with {completed.returncode}\n"
                f"{completed.stdout}{completed.stderr}"
            )
        if completed.stderr:
            raise RuntimeError(f"{name}: unexpected stderr\n{completed.stderr}")
        output = working / output_name
        actual = json.loads(output.read_text(encoding="utf-8"))
        semantic = semantic_result(actual, metadata_fields, copied_source)
        difference = first_difference(expected, semantic)
        if difference:
            raise RuntimeError(f"{name}: reference mismatch at {difference}")
        if output_dir is not None:
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / f"{name}.json").write_bytes(output.read_bytes())
            (output_dir / f"{name}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    elapsed = time.perf_counter() - started
    print(f"PASS {name}: all deterministic result fields match ({elapsed:.2f} s)", flush=True)
    return elapsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", action="append", choices=SUITES,
                        help="Run only this suite; repeat the option to select several")
    parser.add_argument("--output-dir", type=Path,
                        help="Optionally retain raw JSON results and stdout in this directory")
    parser.add_argument("--timeout", type=float, default=600,
                        help="Maximum seconds per suite (default: 600)")
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be a finite positive number")
    names = list(dict.fromkeys(args.suite or SUITES))
    try:
        durations = [run_suite(name, args.output_dir, args.timeout) for name in names]
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 1
    print(f"Verified {len(names)} suite(s) in {sum(durations):.2f} s. Finite checks do not replace the proofs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
