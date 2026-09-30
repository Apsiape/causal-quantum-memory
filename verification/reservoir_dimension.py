"""Exact sufficient dimension selection for the rolling reservoir.

Inputs B >= the analytic buffer b and h >= every export size must be supplied.
This checks dimension/padding inequalities, not the existence of good mixers.
Only integer and rational arithmetic is used for acceptance.
"""
if not __debug__:
    raise SystemExit("This check uses assert statements; run it without python -O.")
import argparse
from fractions import Fraction
import json


def margins(m, r, n, h, B, epsilon, q):
    port = 1 << h
    physical = 1 << q
    D = m * port * ((1 << (q-h)) // m)
    LD = 3 + D.bit_length() + 2*(n+1)*h
    Lc = 3 + (n+2)*r.bit_length()
    AD = 3*(1 + 8*port**2*n*n*LD)
    delta = epsilon / 8
    left = D*D*delta.numerator**2
    right = 32*m*m*r*r*n*n*(1 << B)*AD*Lc*delta.denominator**2
    padding_ok = 4*m*port*epsilon.denominator <= physical*epsilon.numerator
    return D, left, right, padding_ok


def select_dimension(m, r, n, h, B, epsilon, max_qubits=100000):
    epsilon = Fraction(epsilon)
    if min(m, r, n) < 1 or min(h, B) < 0 or not 0 < epsilon <= 1:
        raise ValueError("Positive m,r,n; nonnegative h,B; 0 < epsilon <= 1 required")
    q = h
    # q >= h + ceil(log2(4*m/epsilon)), without a logarithm evaluation.
    while (1 << (q-h))*epsilon.numerator < 4*m*epsilon.denominator:
        q += 1
        if q > max_qubits:
            raise ValueError("Explicit qubit cap exceeded")
    while q <= max_qubits:
        D, left, right, pad = margins(m, r, n, h, B, epsilon, q)
        if pad and left >= right:
            return q, D
        q += 1
    raise ValueError("Explicit qubit cap exceeded")


def self_test():
    cases = [(1, 1, 1, 1, 0, "1"), (3, 4, 2, 2, 3, "1/10"),
             (5, 16, 100, 3, 30, "1/100"), (15, 30, 1000, 4, 100, "1/1000")]
    for m, r, n, h, B, epsilon in cases:
        ep = Fraction(epsilon)
        q, D = select_dimension(m, r, n, h, B, ep)
        _, left, right, pad = margins(m, r, n, h, B, ep, q)
        assert left >= right and pad and D % (m*(1 << h)) == 0
        # The immediately preceding candidate must fail at least one test.
        _, prev_l, prev_r, prev_pad = margins(m, r, n, h, B, ep, q-1)
        assert prev_l < prev_r or not prev_pad
        print(f"Exact reservoir selector: m={m}, r={r}, n={n}, B={B}, epsilon={ep}, q={q}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    for key in ("m", "r", "n", "export-qubits", "buffer-bits"):
        parser.add_argument("--"+key, type=int)
    parser.add_argument("--epsilon")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        if any(x is None for x in (args.m, args.r, args.n, args.export_qubits, args.buffer_bits, args.epsilon)):
            parser.error("Supply all dimension inputs, or use --self-test")
        q, D = select_dimension(args.m, args.r, args.n, args.export_qubits, args.buffer_bits, args.epsilon)
        print(json.dumps({"qubits": q, "active_dimension": D, "physical_dimension": 1 << q}))
