"""Sparse polynomials from PY-2 section 3.4: a polynomial is a dict {power: coefficient}."""
from __future__ import annotations

from collections import defaultdict

Poly = dict[int, float]


def poly_clean(p: Poly) -> Poly:
    """Drop zero coefficients so that equal polynomials are equal dicts."""
    return {k: c for k, c in p.items() if c != 0}


def poly_add(p: Poly, q: Poly) -> Poly:
    """p + q: coefficients add power by power."""
    r = dict(p)
    for k, c in q.items():
        r[k] = r.get(k, 0) + c
    return poly_clean(r)


def poly_mul(p: Poly, q: Poly) -> Poly:
    """p * q: c_i x^i times d_j x^j contributes c_i d_j to the coefficient of x^(i+j)."""
    r: defaultdict[int, float] = defaultdict(int)
    for i, a in p.items():
        for j, b in q.items():
            r[i + j] += a * b
    return poly_clean(r)


def poly_eval(p: Poly, x: float) -> float:
    """Sparse Horner: walk powers downward, multiplying by x**gap between consecutive terms."""
    if not p:
        return 0
    result, prev = 0, max(p)
    for k in sorted(p, reverse=True):
        result = result * x ** (prev - k) + p[k]
        prev = k
    return result * x**prev


def poly_deriv(p: Poly) -> Poly:
    """d/dx: c x^k -> k c x^(k-1); constants vanish."""
    return {k - 1: k * c for k, c in p.items() if k > 0}
