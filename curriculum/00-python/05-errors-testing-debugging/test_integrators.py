"""Tests for integrators.py. Written before the code (PY-5, notes section 3.10).

Run from this folder with ``python -m pytest -q``.
"""

import math
from collections.abc import Callable

import numpy as np
import pytest

import integrators as ig

Fn = Callable[[float], float]
Rule = Callable[..., float]


def observed_order(
    rule: Rule, f: Fn, a: float, b: float, exact: float, ns: tuple[int, ...]
) -> float:
    """Slope of log(error) against log(h): the order of the rule, measured."""
    hs = [(b - a) / n for n in ns]
    errs = [abs(rule(f, a, b, n) - exact) for n in ns]
    return float(np.polyfit(np.log(hs), np.log(errs), 1)[0])


# --- trapezoid: tests written first -------------------------------------------------------
LINES = [
    pytest.param(lambda x: 3.0, 0.0, 2.0, 6.0, id="constant"),
    pytest.param(lambda x: 2 * x + 1, 0.0, 2.0, 6.0, id="line"),
    pytest.param(lambda x: x, -1.0, 3.0, 4.0, id="x-on-[-1,3]"),
]


@pytest.mark.parametrize("n", [1, 2, 7])
@pytest.mark.parametrize("f, a, b, exact", LINES)
def test_trapezoid_exact_on_lines(f: Fn, a: float, b: float, exact: float, n: int) -> None:
    assert ig.trapezoid(f, a, b, n) == pytest.approx(exact, rel=1e-12, abs=1e-12)


def test_trapezoid_error_order() -> None:
    order = observed_order(ig.trapezoid, math.sin, 0.0, math.pi, 2.0, (16, 32, 64, 128))
    assert order == pytest.approx(2.0, abs=0.1)


@pytest.mark.parametrize(
    "n, error", [(0, ValueError), (-3, ValueError), (2.5, TypeError), ("4", TypeError)]
)
def test_trapezoid_rejects_bad_n(n: object, error: type[Exception]) -> None:
    with pytest.raises(error):
        ig.trapezoid(math.sin, 0.0, 1.0, n)  # type: ignore[arg-type]  # a wrong type is the point


# --- Simpson: tests written second ----------------------------------------------------------
CUBICS = [
    pytest.param(lambda x: 1.0, 0.0, 1.0, 1.0, id="1"),
    pytest.param(lambda x: x, 0.0, 2.0, 2.0, id="x"),
    pytest.param(lambda x: x**2, 0.0, 3.0, 9.0, id="x^2"),
    pytest.param(lambda x: x**3 - 2 * x, -1.0, 2.0, 0.75, id="x^3-2x"),
]


@pytest.mark.parametrize("n", [2, 4, 10])
@pytest.mark.parametrize("f, a, b, exact", CUBICS)
def test_simpson_exact_on_cubics(f: Fn, a: float, b: float, exact: float, n: int) -> None:
    assert ig.simpson(f, a, b, n) == pytest.approx(exact, rel=1e-12, abs=1e-12)


def test_simpson_is_richardson_of_two_trapezoids() -> None:
    # S_{2m} = (4 T_{2m} - T_m) / 3 on the same nodes (notes, section 3.10).
    for m in (2, 5, 12):
        s = ig.simpson(math.exp, 0.0, 1.0, 2 * m)
        t = (4 * ig.trapezoid(math.exp, 0.0, 1.0, 2 * m) - ig.trapezoid(math.exp, 0.0, 1.0, m)) / 3
        assert s == pytest.approx(t, rel=1e-13)


def test_simpson_error_order() -> None:
    order = observed_order(ig.simpson, math.sin, 0.0, math.pi, 2.0, (8, 16, 32, 64))
    assert order == pytest.approx(4.0, abs=0.2)


CASES = [
    pytest.param(math.sin, 0.0, math.pi, 2.0, id="sin"),
    pytest.param(math.exp, 0.0, 1.0, math.e - 1.0, id="exp"),
    pytest.param(lambda x: 1 / (1 + x * x), 0.0, 1.0, math.pi / 4, id="arctan"),
    pytest.param(lambda x: x * math.exp(-x * x), 0.0, 1.0, (1 - math.exp(-1.0)) / 2, id="gauss-x"),
    pytest.param(math.cos, 0.0, math.pi / 2, 1.0, id="cos"),
    pytest.param(math.log, 1.0, 2.0, 2 * math.log(2.0) - 1.0, id="log"),
]


@pytest.mark.parametrize("f, a, b, exact", CASES)
def test_simpson_known_integrals(f: Fn, a: float, b: float, exact: float) -> None:
    # Bound (b-a) h^4 M4 / 180 at n = 64 is about 5e-8 relative at worst (sin); rel=1e-6 is 20x.
    assert ig.simpson(f, a, b, 64) == pytest.approx(exact, rel=1e-6)


@pytest.mark.parametrize("n", [1, 3, 7])
def test_simpson_rejects_odd_n(n: int) -> None:
    with pytest.raises(ValueError, match="even"):
        ig.simpson(math.sin, 0.0, 1.0, n)


@pytest.mark.parametrize("n, error", [(0, ValueError), (-2, ValueError), (2.0, TypeError)])
def test_simpson_rejects_bad_n(n: object, error: type[Exception]) -> None:
    with pytest.raises(error):
        ig.simpson(math.sin, 0.0, 1.0, n)  # type: ignore[arg-type]  # a wrong type is the point


# --- regression: added after the planted bug was found ------------------------------------
@pytest.mark.parametrize("k", range(5))
def test_simpson_weights_are_1_4_2_4_1(k: int) -> None:
    """Pin the weight of every node: f is 1 at node k, 0 elsewhere (n = 4 on [0, 1])."""
    nodes = [i / 4 for i in range(5)]  # dyadic, so the grid points are exact floats
    weights = [1, 4, 2, 4, 1]

    def indicator(x: float) -> float:
        return 1.0 if x == nodes[k] else 0.0

    assert ig.simpson(indicator, 0.0, 1.0, 4) == pytest.approx(weights[k] * 0.25 / 3, rel=1e-14)


# --- invariants checked on random inputs (a hand-rolled property test) ----------------------
@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(12345)


@pytest.mark.parametrize("rule", [ig.trapezoid, ig.simpson], ids=["trapezoid", "simpson"])
def test_linearity_reversal_and_additivity(rule: Rule, rng: np.random.Generator) -> None:
    f, g = math.sin, lambda x: x * x + 1.0
    for _ in range(25):
        al, be = (float(v) for v in rng.normal(size=2))
        a, c = (float(v) for v in np.sort(rng.uniform(-2, 2, size=2)))
        b = 0.5 * (a + c)

        def mix(x: float, al: float = al, be: float = be) -> float:
            return al * f(x) + be * g(x)

        # linearity
        expected = al * rule(f, a, c, 8) + be * rule(g, a, c, 8)
        assert rule(mix, a, c, 8) == pytest.approx(expected, rel=1e-10, abs=1e-12)
        # swapping the limits flips the sign
        assert rule(f, c, a, 8) == pytest.approx(-rule(f, a, c, 8), rel=1e-12, abs=1e-14)
        # two half-intervals with 4 panels each use the same nodes as one interval with 8
        halves = rule(f, a, b, 4) + rule(f, b, c, 4)
        assert halves == pytest.approx(rule(f, a, c, 8), rel=1e-10, abs=1e-12)


# --- adaptive Simpson and the exception hierarchy -------------------------------------------
def counting(f: Fn) -> tuple[Fn, list[int]]:
    calls = [0]

    def wrapped(x: float) -> float:
        calls[0] += 1
        return f(x)

    return wrapped, calls


def test_adaptive_matches_scipy_quad_on_a_peak() -> None:
    from scipy.integrate import quad

    def peak(x: float) -> float:
        return math.exp(-100.0 * (x - 0.5) ** 2)

    reference, _ = quad(peak, 0.0, 1.0, epsabs=1e-13, epsrel=1e-13)
    assert abs(ig.adaptive_simpson(peak, 0.0, 1.0, tol=1e-10) - reference) < 1e-8


def test_adaptive_accepts_a_cubic_without_splitting() -> None:
    f, calls = counting(lambda x: x**3 - x)
    assert ig.adaptive_simpson(f, 0.0, 2.0) == pytest.approx(2.0, rel=1e-12)
    assert calls[0] == 5  # 3 evaluations for the whole interval, 2 for the half-panels


def test_adaptive_gives_up_with_a_convergence_error() -> None:
    with pytest.raises(ig.ConvergenceError) as info:
        ig.adaptive_simpson(math.sin, 0.0, 50.0, tol=1e-14, max_depth=3)
    err = info.value
    assert isinstance(err, ig.SolverError) and not isinstance(err, ValueError)
    assert err.iterations == 3 and err.residual > 1e-14


@pytest.mark.parametrize("tol", [0.0, -1e-6, float("nan")])
def test_adaptive_rejects_bad_tol(tol: float) -> None:
    with pytest.raises(ValueError, match="tol"):
        ig.adaptive_simpson(math.sin, 0.0, 1.0, tol=tol)
