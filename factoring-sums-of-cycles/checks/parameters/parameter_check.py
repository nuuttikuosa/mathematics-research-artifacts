"""Reproduce finite exact-arithmetic parameter checks; not an asymptotic proof.

Run with Python 3.10+ and its standard library:
    python -B parameter_check.py
The output is written beside this script, independently of the working directory.
The argument R in parameters() is the paper's formula-size parameter Lambda.
"""
from pathlib import Path
import json


def ceil_log(base, target):
    exponent, power = 0, 1
    while power < target:
        exponent += 1
        power *= base
    return exponent


def parameters(R):
    lam = 3
    while pow(lam, lam) < R:
        lam += 1
    log2lam = (lam - 1).bit_length()
    u = max(3, (lam + log2lam - 1) // log2lam)
    g = max(1, ceil_log(u, R))
    h = 3 * (2 * g - 1)
    k = ceil_log(2 * u - 1, h)
    D = g + 2 * k + 3
    assert lam == 3 or pow(lam - 1, lam - 1) < R
    assert pow(lam, lam) >= R
    assert pow(u, g) >= R and pow(u, g) < R * u
    assert pow(2 * u - 1, k) >= h
    assert k >= 1 and D >= 2 and 3 <= u <= lam
    prime_count = (g + 2 * k + 2) * u + 4
    assert prime_count <= D * (u + 1)
    return dict(R_decimal_digits=len(str(R)), lam=lam, u=u, g=g, h=h,
                k=k, D=D, prime_count=prime_count)


def main():
    inputs = set(range(3, 10001))
    for s in range(3, 101):
        inputs.update([s ** s - 1, s ** s, s ** s + 1])
    for e in list(range(1, 101)) + [200, 500, 1000]:
        inputs.update([10 ** e - 1, 10 ** e, 10 ** e + 1])
    rows = [parameters(R) for R in sorted(inputs)]
    data = dict(
        cases=len(inputs),
        all_exact_arithmetic_checks_passed=True,
        max_k_observed=max(r['k'] for r in rows),
        largest_example=rows[-1],
        note=('Finite arithmetic boundary checks only; asymptotic and '
              'bit-complexity claims are proved in the paper.'),
    )
    out = Path(__file__).resolve().with_name('parameter_check.json')
    out.write_text(json.dumps(data, indent=2) + '\n', encoding='utf8')
    print(out.read_text(encoding='utf8'))


if __name__ == '__main__':
    main()
