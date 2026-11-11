"""Composite quadrature rules, written test-first (PY-5, notes section 3.10).

Public API
----------
trapezoid(f, a, b, n)
    Composite trapezoid rule on ``n`` panels (error ~ h**2).
simpson(f, a, b, n)
    Composite Simpson rule on an even number ``n`` of panels (error ~ h**4).
adaptive_simpson(f, a, b, tol)
    Recursive Simpson rule that splits an interval until a Richardson error
    estimate is below its share of ``tol``.
SolverError, ConvergenceError
    The exception hierarchy of this package: one base class, so a caller can
    write ``except SolverError`` without also catching ``ValueError``.

Every public function has type hints and a NumPy-style docstring, and
``python -m doctest integrators.py`` runs the examples in them.
"""

import logging
import operator
from collections.abc import Callable

__all__ = [
    "ConvergenceError",
    "Integrand",
    "SolverError",
    "adaptive_simpson",
    "simpson",
    "trapezoid",
]

log = logging.getLogger(__name__)

#: A function of one real variable, e.g. ``math.sin``.
Integrand = Callable[[float], float]

#: Teaching switch for PY-5 section 3.10. When ``True``, ``simpson`` drops the
#: last odd-indexed sample (the off-by-one the notes plant). Leave it ``False``.
PLANT_OFF_BY_ONE: bool = False


class SolverError(Exception):
    """Base class for every numerical failure this package raises on purpose."""


class ConvergenceError(SolverError):
    """An iterative method stopped without meeting its tolerance.

    Parameters
    ----------
    message : str
        What went wrong, in a sentence.
    iterations : int
        Iterations (or recursion depth) used when the method gave up.
    residual : float
        The error estimate or residual at that point.
    """

    def __init__(self, message: str, iterations: int, residual: float) -> None:
        super().__init__(message, iterations, residual)
        self.message = message
        self.iterations = iterations
        self.residual = residual

    def __str__(self) -> str:
        return f"{self.message} (iterations={self.iterations}, residual={self.residual:.3g})"


def _panels(n: int, *, even: bool = False) -> int:
    """Validate a panel count and return it as a plain ``int``."""
    n = operator.index(n)  # TypeError for 2.5 or "4"; accepts NumPy integers
    if n < 1:
        raise ValueError(f"need at least 1 panel, got n={n}")
    if even and n % 2:
        raise ValueError(f"Simpson's rule needs an even number of panels, got n={n}")
    return n


def trapezoid(f: Integrand, a: float, b: float, n: int) -> float:
    r"""Composite trapezoid rule for :math:`\int_a^b f(x)\,dx`.

    Parameters
    ----------
    f : callable
        Integrand, called as ``f(x)`` with a float.
    a, b : float
        Limits of integration.
    n : int
        Number of panels, at least 1.

    Returns
    -------
    float
        :math:`h(\tfrac12 y_0 + y_1 + \dots + y_{n-1} + \tfrac12 y_n)`
        with :math:`h = (b-a)/n`.

    Raises
    ------
    TypeError
        If ``n`` is not an integer.
    ValueError
        If ``n < 1``.

    Examples
    --------
    >>> trapezoid(lambda x: x, 0.0, 1.0, 4)
    0.5
    """
    n = _panels(n)
    h = (b - a) / n
    ys = [f(a + i * h) for i in range(n + 1)]
    return h * (0.5 * ys[0] + sum(ys[1:-1]) + 0.5 * ys[-1])


def simpson(f: Integrand, a: float, b: float, n: int) -> float:
    r"""Composite Simpson rule for :math:`\int_a^b f(x)\,dx`.

    Parameters
    ----------
    f : callable
        Integrand, called as ``f(x)`` with a float.
    a, b : float
        Limits of integration.
    n : int
        Number of panels; must be even and at least 2.

    Returns
    -------
    float
        :math:`\tfrac h3(y_0 + 4y_1 + 2y_2 + \dots + 4y_{n-1} + y_n)`
        with :math:`h = (b-a)/n`. Exact for polynomials of degree at most 3.

    Raises
    ------
    TypeError
        If ``n`` is not an integer.
    ValueError
        If ``n`` is odd or ``n < 1``.

    Examples
    --------
    >>> simpson(lambda x: x**3, 0.0, 1.0, 2)
    0.25
    """
    n = _panels(n, even=True)
    h = (b - a) / n
    ys = [f(a + i * h) for i in range(n + 1)]
    stop = -2 if PLANT_OFF_BY_ONE else -1
    odd = sum(ys[1:stop:2])  # y_1, y_3, ..., y_{n-1}: weight 4
    even = sum(ys[2:-1:2])  # y_2, y_4, ..., y_{n-2}: weight 2
    return h / 3 * (ys[0] + 4 * odd + 2 * even + ys[-1])


def adaptive_simpson(
    f: Integrand, a: float, b: float, tol: float = 1e-10, max_depth: int = 40
) -> float:
    r"""Adaptive Simpson quadrature of :math:`\int_a^b f(x)\,dx`.

    On each interval compare one Simpson panel :math:`S_1` with two half-size
    panels :math:`S_2`. If :math:`|S_2 - S_1| \le 15\,\mathrm{tol}` accept
    :math:`S_2 + (S_2 - S_1)/15` (Richardson extrapolation); otherwise split
    the interval and give each half ``tol / 2``.

    Parameters
    ----------
    f : callable
        Integrand, called as ``f(x)`` with a float.
    a, b : float
        Limits of integration.
    tol : float, optional
        Target absolute error for the whole interval; must be positive.
    max_depth : int, optional
        Deepest allowed bisection level.

    Returns
    -------
    float
        The integral estimate.

    Raises
    ------
    ValueError
        If ``tol`` is not positive.
    ConvergenceError
        If an interval still fails its test at ``max_depth`` (or the
        integrand returns NaN, which fails every comparison).

    Examples
    --------
    >>> import math
    >>> round(adaptive_simpson(math.sin, 0.0, math.pi, tol=1e-9), 8)
    2.0
    """
    if not tol > 0:
        raise ValueError(f"tol must be positive, got {tol!r}")
    stats = {"evals": 3, "depth": 0}

    def panel(lo: float, hi: float, f_lo: float, f_mid: float, f_hi: float) -> float:
        return (hi - lo) / 6 * (f_lo + 4 * f_mid + f_hi)

    def recurse(
        lo: float,
        hi: float,
        f_lo: float,
        f_mid: float,
        f_hi: float,
        whole: float,
        eps: float,
        depth: int,
    ) -> float:
        mid = 0.5 * (lo + hi)
        f_left, f_right = f(0.5 * (lo + mid)), f(0.5 * (mid + hi))
        stats["evals"] += 2
        stats["depth"] = max(stats["depth"], depth)
        left = panel(lo, mid, f_lo, f_left, f_mid)
        right = panel(mid, hi, f_mid, f_right, f_hi)
        delta = left + right - whole
        if abs(delta) <= 15 * eps:
            return left + right + delta / 15
        if depth >= max_depth:
            raise ConvergenceError("adaptive_simpson reached max_depth", depth, abs(delta) / 15)
        log.debug(
            "depth %d: [%.5f, %.5f] |delta|=%.2e > %.2e, split", depth, lo, hi, abs(delta), 15 * eps
        )
        return recurse(lo, mid, f_lo, f_left, f_mid, left, eps / 2, depth + 1) + recurse(
            mid, hi, f_mid, f_right, f_hi, right, eps / 2, depth + 1
        )

    f_a, f_m, f_b = f(a), f(0.5 * (a + b)), f(b)
    result = recurse(a, b, f_a, f_m, f_b, panel(a, b, f_a, f_m, f_b), tol, 0)
    log.info(
        "adaptive_simpson on [%g, %g]: tol=%g, %d evaluations, deepest level %d",
        a,
        b,
        tol,
        stats["evals"],
        stats["depth"],
    )
    return result
