"""Scalar root finders from PY-1: bisection, Newton and secant.

Each returns ``(root, iterations)``.  They take the function ``f`` as an ordinary
argument, which is the whole point of PY-3 section 3.2.
"""
from __future__ import annotations

from collections.abc import Callable

Func = Callable[[float], float]


def bisect(f: Func, a: float, b: float, tol: float = 1e-12, max_iter: int = 200) -> tuple[float, int]:
    """Root of f in [a, b] by halving the bracket (iteration count: PY-1 eq. 3.9); f(a) f(b) < 0 is required."""
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError("f(a) and f(b) must have opposite signs")
    for k in range(1, max_iter + 1):
        m = a + (b - a) / 2
        fm = f(m)
        if fm == 0 or (b - a) / 2 <= tol:
            return m, k
        if fa * fm < 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    raise RuntimeError("bisection did not reach the tolerance")


def newton(f: Func, df: Func, x0: float, tol: float = 1e-12, max_iter: int = 50) -> tuple[float, int]:
    """Root of f near x0 by x <- x - f(x)/f'(x) (PY-1 eq. 3.10)."""
    x = x0
    for k in range(1, max_iter + 1):
        step = f(x) / df(x)
        x -= step
        if abs(step) <= tol:
            return x, k
    raise RuntimeError("Newton's method did not converge")


def secant(f: Func, x0: float, x1: float, tol: float = 1e-12, max_iter: int = 50) -> tuple[float, int]:
    """Root of f from two starting points; Newton with f' replaced by a difference quotient."""
    f0, f1 = f(x0), f(x1)
    for k in range(1, max_iter + 1):
        if f1 == f0:
            raise ZeroDivisionError("secant line is horizontal")
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(x2 - x1) <= tol:
            return x2, k
        x0, f0, x1, f1 = x1, f1, x2, f(x2)
    raise RuntimeError("secant method did not converge")
