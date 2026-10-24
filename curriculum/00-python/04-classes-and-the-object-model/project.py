# %% [markdown]
# # PY-4 Project: `Vec`, `Rational`, `Dual`, `Matrix`, and why automatic differentiation exists
#
# Builds four classes in pure Python (notes §3.2-3.8), tests each against a reference
# (LA-1 identities, `fractions.Fraction`, SymPy and finite differences, `numpy.linalg.solve`),
# and measures the one result that justifies dual numbers: a difference quotient has an error
# floor near 1e-8 whatever the step, a `Dual` does not (figure `error_vs_h.png`). Then it
# covers the dataclass, ABC and protocol examples and asserts every number quoted in §4.
# NumPy, SymPy and Matplotlib appear only in the checks. Run headless:
# `MPLBACKEND=Agg python3 project.py`. Importing this file re-runs the checks (a few seconds).

# %%
import functools
import math
import numbers
import operator
import time
import tracemalloc
from abc import ABC, abstractmethod
from collections import UserDict
from collections.abc import Callable, Iterable
from dataclasses import FrozenInstanceError, dataclass, field, replace
from fractions import Fraction
from math import gcd
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from numpy.polynomial import Polynomial as NPPoly

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.perf_counter()

# %% [markdown]
# ## Part 0: The classes (notes §3.2-§3.8)
# `Vec` (immutable, hashable), `Rational` (normalized, hash-compatible with int), `Dual` with
# `dmath` (§3.7), and the mutable `Matrix` (§3.8). Read the docstrings; the notes derive them.

# %%
class Vec:
    """An immutable vector in R^n; operators follow LA-1 (§3.1, §3.5)."""

    def __init__(self, *components):
        self._c = tuple(float(c) for c in components)

    def __repr__(self):
        return f"Vec({', '.join(map(repr, self._c))})"

    def __str__(self):
        return "(" + ", ".join(f"{c:g}" for c in self._c) + ")"

    def __eq__(self, other):
        if not isinstance(other, Vec):
            return NotImplemented
        return self._c == other._c

    def __hash__(self):
        return hash(self._c)

    def __add__(self, other):
        if not isinstance(other, Vec):
            return NotImplemented
        return Vec(*(a + b for a, b in zip(self._c, other._c, strict=True)))

    def __sub__(self, other):
        return self + (-other) if isinstance(other, Vec) else NotImplemented

    def __neg__(self):
        return Vec(*(-a for a in self._c))

    def __mul__(self, k):
        if not isinstance(k, numbers.Real):
            return NotImplemented
        return Vec(*(k * a for a in self._c))

    __rmul__ = __mul__

    def __truediv__(self, k):
        if not isinstance(k, numbers.Real):
            return NotImplemented
        return Vec(*(a / k for a in self._c))

    def __abs__(self):
        return math.hypot(*self._c)

    def dot(self, other):
        return sum(a * b for a, b in zip(self._c, other._c, strict=True))

    def __matmul__(self, other):
        return self.dot(other) if isinstance(other, Vec) else NotImplemented

    def cross(self, other):
        if len(self._c) != 3 or len(other._c) != 3:
            raise ValueError("cross product needs two 3-vectors")
        (a1, a2, a3), (b1, b2, b3) = self._c, other._c
        return Vec(a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1)

    def __iter__(self):
        return iter(self._c)

    def __len__(self):
        return len(self._c)

    def __bool__(self):
        return any(self._c)

    x = property(lambda self: self._c[0])
    y = property(lambda self: self._c[1])
    z = property(lambda self: self._c[2])

    @classmethod
    def from_polar(cls, r, theta):
        return cls(r * math.cos(theta), r * math.sin(theta))


def _coerced(method):
    """Accept a Rational or an int (promoted to Rational); decline anything else."""
    @functools.wraps(method)
    def wrapper(self, other):
        if isinstance(other, numbers.Integral):
            other = Rational(other)
        elif not isinstance(other, Rational):
            return NotImplemented
        return method(self, other)
    return wrapper


@functools.total_ordering
class Rational:
    """n/d in lowest terms with d > 0; equal to ints, so it hashes like them."""

    def __init__(self, n, d=1):
        n, d = operator.index(n), operator.index(d)
        if d == 0:
            raise ZeroDivisionError("Rational with zero denominator")
        self._n, self._d = Rational._lowest_terms(n, d)

    @staticmethod
    def _lowest_terms(n, d):
        g = gcd(n, d) if d > 0 else -gcd(n, d)
        return n // g, d // g

    numerator = property(lambda self: self._n)
    denominator = property(lambda self: self._d)

    def __repr__(self):
        return f"Rational({self._n}, {self._d})"

    def __str__(self):
        return str(self._n) if self._d == 1 else f"{self._n}/{self._d}"

    @_coerced
    def __eq__(self, other):
        return (self._n, self._d) == (other._n, other._d)

    def __hash__(self):
        return hash(self._n) if self._d == 1 else hash((self._n, self._d))

    @_coerced
    def __lt__(self, other):
        return self._n * other._d < other._n * self._d

    def __bool__(self):
        return self._n != 0

    def __neg__(self):
        return Rational(-self._n, self._d)

    @_coerced
    def __add__(self, other):
        return Rational(self._n * other._d + other._n * self._d, self._d * other._d)

    __radd__ = __add__

    @_coerced
    def __sub__(self, other):
        return self + (-other)

    @_coerced
    def __rsub__(self, other):
        return other + (-self)

    @_coerced
    def __mul__(self, other):
        return Rational(self._n * other._n, self._d * other._d)

    __rmul__ = __mul__

    @_coerced
    def __truediv__(self, other):
        return Rational(self._n * other._d, self._d * other._n)

    @_coerced
    def __rtruediv__(self, other):
        return other / self

    def __pow__(self, k):
        if not isinstance(k, int):
            return NotImplemented
        return Rational(self._n ** k, self._d ** k) if k >= 0 else Rational(self._d ** -k, self._n ** -k)


def _dual_operand(method):
    """Promote a real number to Dual(number, 0); decline anything else."""
    @functools.wraps(method)
    def wrapper(self, other):
        if isinstance(other, numbers.Real):
            other = Dual(other, 0.0)
        elif not isinstance(other, Dual):
            return NotImplemented
        return method(self, other)
    return wrapper


class Dual:
    """a + b*eps with eps**2 = 0: a value a and its derivative b (notes §3.7)."""
    __slots__ = ("a", "b")

    def __init__(self, a, b=0.0):
        self.a, self.b = a, b

    def __repr__(self):
        return f"Dual({self.a!r}, {self.b!r})"

    def __eq__(self, other):
        if not isinstance(other, Dual):
            return NotImplemented
        return (self.a, self.b) == (other.a, other.b)

    def __hash__(self):
        return hash((self.a, self.b))

    def __neg__(self):
        return Dual(-self.a, -self.b)

    @_dual_operand
    def __add__(self, o):
        return Dual(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    @_dual_operand
    def __sub__(self, o):
        return self + (-o)

    @_dual_operand
    def __rsub__(self, o):
        return o + (-self)

    @_dual_operand
    def __mul__(self, o):                                         # (3.1)
        return Dual(self.a * o.a, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    @_dual_operand
    def __truediv__(self, o):                                     # (3.2)
        return Dual(self.a / o.a, (self.b * o.a - self.a * o.b) / (o.a * o.a))

    @_dual_operand
    def __rtruediv__(self, o):
        return o / self

    def __pow__(self, n):                                         # integer powers
        if not isinstance(n, int):
            return NotImplemented
        if n == 0:
            return Dual(1.0, 0.0)
        return Dual(self.a ** n, n * self.a ** (n - 1) * self.b)


def _elementary(f, df):
    """Dual-aware f: f(a + b eps) = f(a) + df(a) * b eps  (3.3)."""
    def wrapped(x):
        if isinstance(x, Dual):
            return Dual(wrapped(x.a), df(x.a) * x.b)
        return f(x)
    return wrapped


dmath = SimpleNamespace()                       # stands in for a module dmath.py
dmath.sin = _elementary(math.sin, lambda a: dmath.cos(a))
dmath.cos = _elementary(math.cos, lambda a: -dmath.sin(a))
dmath.exp = _elementary(math.exp, lambda a: dmath.exp(a))
dmath.log = _elementary(math.log, lambda a: 1 / a)
dmath.sqrt = _elementary(math.sqrt, lambda a: 0.5 / dmath.sqrt(a))
dmath.tanh = _elementary(math.tanh, lambda a: 1 - dmath.tanh(a) ** 2)
dmath.atan = _elementary(math.atan, lambda a: 1 / (1 + a * a))


def derivative(f, x):
    """Return (f(x), f'(x)) from one evaluation of f at Dual(x, 1)."""
    y = f(Dual(x, 1.0))
    return (y.a, y.b) if isinstance(y, Dual) else (y, 0.0)


def jacobian(f, x):
    """Rows of the Jacobian of f: R^n -> R^m, one forward pass per input (LA-14 preview)."""
    cols = []
    for j in range(len(x)):
        out = f([Dual(xi, 1.0 if i == j else 0.0) for i, xi in enumerate(x)])
        cols.append([y.b if isinstance(y, Dual) else 0.0 for y in out])
    return [list(row) for row in zip(*cols)]


def newton_ad(f, x0, tol=1e-14, max_iter=50):
    """Newton's method where f'(x) is the Dual part of f(Dual(x, 1))."""
    x = x0
    for k in range(1, max_iter + 1):
        y = f(Dual(x, 1.0))
        step = y.a / y.b
        x -= step
        if abs(step) <= tol:
            return x, k
    raise RuntimeError("Newton did not converge")


class Matrix:
    """Dense matrix, list of float rows; mutable, hence unhashable (notes §3.8)."""

    def __init__(self, rows):
        self._rows = [[float(x) for x in row] for row in rows]
        if len({len(r) for r in self._rows}) > 1:
            raise ValueError("rows have different lengths")

    shape = property(lambda self: (len(self._rows), len(self._rows[0])))
    T = property(lambda self: Matrix(zip(*self._rows)))

    def __repr__(self):
        return f"Matrix({self._rows!r})"

    def __eq__(self, other):
        if not isinstance(other, Matrix):
            return NotImplemented
        return self._rows == other._rows

    __hash__ = None                                  # mutable, so unhashable

    def __getitem__(self, key):
        if isinstance(key, tuple):                   # A[i, j]
            i, j = key
            return self._rows[i][j]
        return Vec(*self._rows[key])                 # A[i] is row i

    def __setitem__(self, key, value):
        i, j = key
        self._rows[i][j] = float(value)

    def __iadd__(self, other):
        if not isinstance(other, Matrix):
            return NotImplemented
        if self.shape != other.shape:
            raise ValueError(f"shapes {self.shape} and {other.shape} differ")
        for r, s in zip(self._rows, other._rows):
            r[:] = [a + b for a, b in zip(r, s)]
        return self                                  # forgetting this binds the name to None

    def __matmul__(self, other):
        if isinstance(other, Vec):
            return Vec(*(sum(a * b for a, b in zip(row, other)) for row in self._rows))
        if not isinstance(other, Matrix):
            return NotImplemented
        if self.shape[1] != other.shape[0]:
            raise ValueError(f"shapes {self.shape} and {other.shape} do not align")
        cols = list(zip(*other._rows))
        return Matrix([[sum(a * b for a, b in zip(row, col)) for col in cols] for row in self._rows])

    def lu(self):
        """Return (perm, L, U) with P A = L U; row i of P A is row perm[i] of A (LA-2 Thm 3.13)."""
        n = self.shape[0]
        if self.shape[1] != n:
            raise ValueError("LU needs a square matrix")
        U = [row[:] for row in self._rows]
        L = [[0.0] * n for _ in range(n)]
        perm = list(range(n))
        for k in range(n):
            p = max(range(k, n), key=lambda i: abs(U[i][k]))        # partial pivoting
            if U[p][k] == 0.0:
                raise ValueError("matrix is singular")
            if p != k:
                U[k], U[p] = U[p], U[k]
                L[k], L[p] = L[p], L[k]                              # swap the multipliers found so far
                perm[k], perm[p] = perm[p], perm[k]
            for i in range(k + 1, n):
                L[i][k] = U[i][k] / U[k][k]                          # the multiplier
                for j in range(k, n):
                    U[i][j] -= L[i][k] * U[k][j]
        for i in range(n):
            L[i][i] = 1.0
        return perm, Matrix(L), Matrix(U)

    def solve(self, b):
        """Solve A x = b: L y = P b forward (LA-2 (3.4)), then U x = y backward (LA-2 (3.3))."""
        perm, L, U = self.lu()
        n, b = len(perm), list(b)
        y = [0.0] * n
        for i in range(n):
            y[i] = float(b[perm[i]]) - sum(L._rows[i][j] * y[j] for j in range(i))
        x = [0.0] * n
        for i in reversed(range(n)):
            x[i] = (y[i] - sum(U._rows[i][j] * x[j] for j in range(i + 1, n))) / U._rows[i][i]
        return Vec(*x)


@dataclass(frozen=True, order=True)
class Measurement:
    value: float
    unit: str
    time: float


@dataclass(frozen=True)
class Polynomial:
    """c[0] + c[1] x + ...; ascending coefficients, trailing zeros removed."""
    coeffs: tuple

    def __post_init__(self):
        c = [float(a) for a in self.coeffs]
        while len(c) > 1 and c[-1] == 0.0:
            c.pop()
        object.__setattr__(self, "coeffs", tuple(c))

    def __call__(self, x):
        acc = 0.0
        for c in reversed(self.coeffs):
            acc = acc * x + c
        return acc

    def deriv(self):
        return Polynomial(tuple(k * c for k, c in enumerate(self.coeffs))[1:] or (0.0,))


print(f"Part 0: classes defined ({time.perf_counter() - T0:.2f} s)")

# %% [markdown]
# ## Part 1: `Vec` against LA-1 and the data model (notes §3.1-§3.5, Examples 4.1-4.3)
# Every method returns a new `Vec` or `NotImplemented`. The random checks are LA-1's theorems:
# Cauchy-Schwarz (Theorem 3.11), the parallelogram law, Lagrange's identity, `u x v` orthogonal.

# %%
def raises(exc, fn):
    try:
        fn()
    except exc:
        return True
    return False


v, w = Vec(1, 2, 3), Vec(4, 5, 6)
assert repr(v) == "Vec(1.0, 2.0, 3.0)" and str(v) == "(1, 2, 3)" and eval(repr(v)) == v
assert hash(v) == hash(Vec(1.0, 2.0, 3.0)) and len({v, Vec(1, 2, 3), w}) == 2
assert (v + w, v - w, -v, v * 2, 2 * v) == (Vec(5, 7, 9), Vec(-3, -3, -3), Vec(-1, -2, -3), Vec(2, 4, 6), Vec(2, 4, 6))
assert 2 * v + w / 2 == Vec(4, 6.5, 9) and abs(Vec(3, 4)) == 5.0 and v @ w == 32.0 and v.cross(w) == Vec(-3, 6, -3)
assert v.cross(w) @ v == 0.0 and v.cross(w) @ w == 0.0
x, y, z = v
assert (x, y, z) == (v.x, v.y, v.z) == (1.0, 2.0, 3.0) and len(v) == 3 and not Vec(0, 0) and max(Vec(3, 9, 4)) == 9.0
p = Vec.from_polar(2, math.pi / 3)
assert str(p) == "(1, 1.73205)" and math.isclose(abs(p), 2.0)
assert raises(TypeError, lambda: v * w) and raises(TypeError, lambda: v + 1) and raises(ValueError, lambda: v + Vec(1, 2))
assert v.__eq__([1, 2, 3]) is NotImplemented and v != [1, 2, 3] and raises(TypeError, lambda: sorted([v, w]))
v.__add__ = lambda other: "hijacked"            # special-method lookup is on the type, not the instance
assert v + w == Vec(5, 7, 9)
del v.__add__

U, W = rng.normal(size=(2000, 3)), rng.normal(size=(2000, 3))
for a, b in zip(U, W):
    u, t = Vec(*a), Vec(*b)
    c = u.cross(t)
    assert abs(u @ t) <= abs(u) * abs(t) + 1e-12                                            # Cauchy-Schwarz
    assert math.isclose(abs(u + t) ** 2 + abs(u - t) ** 2, 2 * abs(u) ** 2 + 2 * abs(t) ** 2, rel_tol=1e-12)
    assert math.isclose(abs(c) ** 2, abs(u) ** 2 * abs(t) ** 2 - (u @ t) ** 2, rel_tol=1e-9, abs_tol=1e-12)
    assert abs(c @ u) < 1e-9 and abs(c @ t) < 1e-9
    assert np.allclose(list(c), np.cross(a, b)) and math.isclose(u @ t, float(a @ b), abs_tol=1e-12)
    assert math.isclose(abs(u), np.linalg.norm(a))

LOG = []                                          # Example 4.3: which method runs, in which order


def logged(name, result=NotImplemented):
    def method(self, other):
        LOG.append(name)
        return result
    return method


Left = type("Left", (), {"__add__": logged("Left.__add__"), "__radd__": logged("Left.__radd__", "L")})
Right = type("Right", (), {"__radd__": logged("Right.__radd__", "R")})
SubLeft = type("SubLeft", (Left,), {"__radd__": logged("SubLeft.__radd__", "S")})
assert Left() + Right() == "R" and LOG == ["Left.__add__", "Right.__radd__"]
LOG.clear()
assert 1 + Right() == "R" and LOG == ["Right.__radd__"]
LOG.clear()
assert Left() + SubLeft() == "S" and LOG == ["SubLeft.__radd__"]
LOG.clear()
assert raises(TypeError, lambda: Left() + Left()) and LOG == ["Left.__add__"]


class IndexedVec(Vec):                            # NumPy and the sequence protocol (§3.10)
    def __getitem__(self, i): return self._c[i]


class QuietVec(Vec):
    __array_ufunc__ = None


arr = np.array([1.0, 2.0, 3.0])
r = arr * Vec(1, 2, 3)
assert np.asarray(v).shape == () and np.asarray(IndexedVec(1, 2, 3)).shape == (3,)
assert r.dtype == object and isinstance(r[0], Vec) and raises(TypeError, lambda: arr * QuietVec(1, 2, 3))
assert isinstance(v, Iterable) and not isinstance(3, Iterable)
print(f"Part 1: Vec ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 2: `Rational` against `fractions.Fraction` (notes §3.6)

# %%
a = Rational(6, -8)
assert (str(a), repr(a), str(Rational(0, 5))) == ("-3/4", "Rational(-3, 4)", "0") and Rational(2, 4) == Rational(1, 2)
assert sorted([Rational(1, 2), Rational(1, 3), Rational(-1, 2)]) == [Rational(-1, 2), Rational(1, 3), Rational(1, 2)]
assert Rational(3) == 3 and hash(Rational(3)) == hash(3) and {3: "x"}[Rational(3)] == "x" and Rational(1, 2) != 0.5
assert Fraction(1, 2) == 0.5 and hash(Fraction(1, 2)) == hash(0.5)
assert [str(t) for t in (Rational(1, 3) + Rational(1, 6), 1 - Rational(1, 3), 2 * Rational(3, 4), 3 / Rational(1, 2),
                         Rational(2, 3) ** -2)] == ["1/2", "2/3", "3/2", "6", "9/4"]
assert raises(TypeError, lambda: Rational(0.5)) and raises(ZeroDivisionError, lambda: Rational(1, 0))
assert not Rational(0) and Rational(1, 2).numerator == 1 and Rational(1, 2).denominator == 2


def same(r, f):
    return (r.numerator, r.denominator) == (f.numerator, f.denominator)


for p, q, r_, s in rng.integers(-50, 51, size=(1000, 4)):
    if q == 0 or s == 0:
        continue
    A, B, FA, FB = Rational(p, q), Rational(r_, s), Fraction(int(p), int(q)), Fraction(int(r_), int(s))
    assert same(A + B, FA + FB) and same(A - B, FA - FB) and same(A * B, FA * FB)
    assert (B == 0 or same(A / B, FA / FB)) and (A < B) == (FA < FB) and (A == B) == (FA == FB)
    assert (A <= B) == (FA <= FB) and (A > B) == (FA > FB) and (A >= B) == (FA >= FB)
    assert same(A + int(p), FA + int(p)) and same(int(p) - A, int(p) - FA)
print(f"Part 2: Rational ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 3: `Dual` (notes §3.7, Examples 4.4 and 4.7)
# (a) The algebra (3.1)-(3.2) against 2x2 matrices (eps = [[0,1],[0,0]]); (b) ten functions, each
# against SymPy `lambdify` and a central difference, with errors measured against SymPy evaluated
# to 40 digits at the exact double `x0`; (c) the V-shaped error curve of a difference quotient;
# (d) a Jacobian by one pass per input; (e) Newton's method with no hand-written derivative.

# %%
d3 = Dual(3.0, 1.0)
assert d3 * d3 * d3 == Dual(27.0, 27.0) and d3 ** 3 == Dual(27.0, 27.0) and 2 * d3 - 5 == Dual(1.0, 2.0)
assert (d3 + 1) / (d3 - 1) == Dual(2.0, -0.5) and 1 / d3 == Dual(1 / 3, -1 / 9) and d3 ** 0 == Dual(1.0, 0.0)
assert Dual(1, 2) != 1 and Dual(1, 0).__eq__(1) is NotImplemented and raises(TypeError, lambda: Dual(1, 1) + "a")
assert raises(TypeError, lambda: Dual(1, 1) ** 0.5) and hash(Dual(1.0, 2.0)) == hash(Dual(1.0, 2.0))
sx = dmath.sin(Dual(0.5, 1.0))
assert math.isclose(sx.a, math.sin(0.5)) and math.isclose(sx.b, math.cos(0.5))
sq = dmath.sin(Dual(1.0, 1.0) ** 2)                                       # chain rule: d/dx sin(x^2) = 2x cos(x^2)
assert math.isclose(sq.a, math.sin(1.0)) and math.isclose(sq.b, 2 * math.cos(1.0))
assert not hasattr(Dual(1.0), "__dict__") and raises(AttributeError, lambda: setattr(Dual(1.0), "c", 3))


def as_matrix(d):                                  # a + b eps  <->  aI + bE
    return np.array([[d.a, d.b], [0.0, d.a]])


for _ in range(100):
    p, q = Dual(*rng.normal(size=2)), Dual(*rng.normal(size=2))
    assert np.allclose(as_matrix(p) @ as_matrix(q), as_matrix(p * q))
    assert np.allclose(as_matrix(p) @ np.linalg.inv(as_matrix(q)), as_matrix(p / q))

FUNCS = {
    "x^3 - 2x + 1": lambda m, x: x ** 3 - 2 * x + 1,
    "sin x cos x": lambda m, x: m.sin(x) * m.cos(x),
    "exp(-x^2)": lambda m, x: m.exp(-x * x),
    "log(1 + x^2)": lambda m, x: m.log(1 + x * x),
    "sqrt(1 + x^2)": lambda m, x: m.sqrt(1 + x * x),
    "x sin x / (1 + x^2)": lambda m, x: x * m.sin(x) / (1 + x ** 2),
    "tanh x": lambda m, x: m.tanh(x),
    "atan(x^2)": lambda m, x: m.atan(x ** 2),
    "exp(x) cos x": lambda m, x: m.exp(x) * m.cos(x),
    "sigmoid x": lambda m, x: 1 / (1 + m.exp(-x)),
}
X, PTS = sp.symbols("x"), [0.3, 0.7, 1.3, 1.9, 2.5]
errs = {}
for name, F in FUNCS.items():
    expr = F(sp, X)
    dexpr = sp.diff(expr, X)
    f_num, df_num = sp.lambdify(X, expr, "math"), sp.lambdify(X, dexpr, "math")
    e_dual = e_fd = 0.0
    for x0 in PTS:
        val, der = derivative(lambda t: F(dmath, t), x0)
        h = 1e-5
        fd = (F(math, x0 + h) - F(math, x0 - h)) / (2 * h)
        exact = float(dexpr.subs(X, sp.Rational(x0)).evalf(40))               # exact double x0, 40 digits
        assert abs(val - f_num(x0)) <= 1e-13 * max(1.0, abs(val)) and abs(der - df_num(x0)) <= 1e-12 * max(1.0, abs(der))
        e_dual = max(e_dual, abs(der - exact) / max(1.0, abs(exact)))
        e_fd = max(e_fd, abs(fd - exact) / max(1.0, abs(exact)))
    assert e_dual < 1e-14 and e_fd < 1e-7
    errs[name] = (e_dual, e_fd)
    print(f"  {name:22s} Dual {e_dual:8.1e}   central difference (h=1e-5) {e_fd:8.1e}")

# Example 4.4 and §1: f(x) = x sin x / (1 + x^2) at x = 1.3
f = lambda t: t * math.sin(t) / (1 + t ** 2)
x0 = 1.3
u, vv = Dual(1.3, 1.0) * dmath.sin(Dual(1.3, 1.0)), 1 + Dual(1.3, 1.0) ** 2
yy = u / vv
f_exact = sp.lambdify(X, FUNCS["x sin x / (1 + x^2)"](sp, X), "math")(1.3)
df_exact = float(sp.diff(FUNCS["x sin x / (1 + x^2)"](sp, X), X).subs(X, sp.Rational(1.3)).evalf(40))
assert abs(yy.a - f_exact) < 1e-14 and abs(yy.b - df_exact) < 1e-14
assert math.isclose(u.a, 1.3 * math.sin(1.3)) and math.isclose(u.b, math.sin(1.3) + 1.3 * math.cos(1.3))
assert math.isclose(vv.a, 2.69) and math.isclose(vv.b, 2.6)
assert math.isclose((u.b * vv.a - u.a * vv.b) / vv.a ** 2, yy.b, rel_tol=1e-14)         # (3.2) by hand
assert abs((1.31131 * 2.69 - 1.25263 * 2.6) / 2.69 ** 2 - 0.037394) < 1e-5
assert abs(df_exact - 0.0373942117663269) < 1e-16

hs = np.logspace(-1, -14, 120)
fwd = np.array([abs((f(x0 + h) - f(x0)) / h - df_exact) for h in hs])
cen = np.array([abs((f(x0 + h) - f(x0 - h)) / (2 * h) - df_exact) for h in hs])
dual_err = abs(yy.b - df_exact)
for k, e_f, e_c in zip((2, 5, 8, 11), (3.41e-3, 3.43e-6, 9.29e-9, 3.65e-6), (1.17e-5, 5.19e-12, 6.51e-9, 1.90e-6)):
    h = 10.0 ** -k                                  # the §1 table (orders of magnitude; last digits are platform noise)
    assert 0.1 < abs((f(x0 + h) - f(x0)) / h - df_exact) / e_f < 10
    assert 0.01 < abs((f(x0 + h) - f(x0 - h)) / (2 * h) - df_exact) / e_c < 100
print(f"  forward floor {fwd.min():.1e} at h={hs[fwd.argmin()]:.0e}; "
      f"central floor {cen.min():.1e} at h={hs[cen.argmin()]:.0e}; Dual {dual_err:.1e}")
assert 1e-10 < fwd.min() < 1e-7 and 1e-9 < hs[fwd.argmin()] < 1e-6      # floor near sqrt(eps) ~ 1.5e-8
assert cen.min() < 1e-9 and 1e-7 < hs[cen.argmin()] < 1e-4              # floor near eps^(2/3), at h ~ eps^(1/3)
assert min(fwd[0], fwd[-1]) > 1e3 * fwd.min() and min(cen[0], cen[-1]) > 1e3 * cen.min()
assert dual_err < 1e-15 and fwd.min() > 1e6 * dual_err                  # the V: truncation falls, rounding rises

fig, ax = plt.subplots(figsize=(6.4, 4))
ax.loglog(hs, fwd, label="forward difference")
ax.loglog(hs, cen, label="central difference")
ax.axhline(max(dual_err, 1e-18), color="k", ls="--", label="Dual (no step size)")
ax.set(xlabel="step size h", ylabel="|error| in f'(1.3)", title="Why automatic differentiation exists")
ax.legend()
fig.tight_layout()
fig.savefig(FIG / "error_vs_h.png", dpi=120)
plt.close(fig)

polar = lambda q: [q[0] * dmath.cos(q[1]), q[0] * dmath.sin(q[1])]
J = jacobian(polar, [2.0, 0.5])
assert np.allclose(J, [[math.cos(0.5), -2 * math.sin(0.5)], [math.sin(0.5), 2 * math.cos(0.5)]])
assert math.isclose(J[0][0] * J[1][1] - J[0][1] * J[1][0], 2.0)
g = lambda q: [q[0] * q[1] + dmath.sin(q[2]), q[0] ** 2 / (1 + q[1] * q[1]), 3.0]     # R^3 -> R^3, one constant
sx_, sy_, sz_ = sp.symbols("a b c")
Jsym = sp.Matrix([sx_ * sy_ + sp.sin(sz_), sx_ ** 2 / (1 + sy_ ** 2), 3]).jacobian([sx_, sy_, sz_])
assert np.allclose(jacobian(g, [0.4, 1.1, 2.0]), np.array(Jsym.subs({sx_: 0.4, sy_: 1.1, sz_: 2.0}), dtype=float))

root, k1 = newton_ad(lambda t: dmath.cos(t) - t, 1.0)
e_, M_ = 0.6, 1.0
E, k2 = newton_ad(lambda t: t - e_ * dmath.sin(t) - M_, M_)
assert abs(root - 0.7390851332151607) < 1e-15 and k1 == 5 and k2 == 6
assert abs(E - e_ * math.sin(E) - M_) < 1e-15 and abs(E - 1.5997485482275295) < 1e-14
print(f"Part 3: Dual ok ({time.perf_counter() - T0:.1f} s); figure saved")

# %% [markdown]
# ## Part 4: `Matrix` LU against NumPy (notes §3.8, Example 4.5)
# `solve` is LA-2's elimination with partial pivoting. Checks: Example 4.5 to the digit, 200 random
# 6x6 systems, the singular case, in-place `+=`, and the cost of pure Python (measured, not quoted).

# %%
A = Matrix([[2, 1, 1], [4, -6, 0], [-2, 7, 2]])
b = [5, -2, 9]
perm, L, U = A.lu()
assert perm == [1, 0, 2] and L == Matrix([[1, 0, 0], [0.5, 1, 0], [-0.5, 1, 1]])
assert U == Matrix([[4, -6, 0], [0, 4, 1], [0, 0, 1]])
assert Matrix([A._rows[i] for i in perm]) == L @ U
y = L.solve([b[i] for i in perm])
assert y == Vec(-2, 6, 2) and U.solve(y) == Vec(1, 1, 2) and A.solve(b) == Vec(1, 1, 2)
assert np.allclose(np.linalg.solve(np.array(A._rows), np.array(b)), [1, 1, 2]) and A @ A.solve(b) - Vec(*b) == Vec(0, 0, 0)

worst = 0.0
for _ in range(200):
    M_, rhs = rng.normal(size=(6, 6)), rng.normal(size=6)
    perm, L, U = Matrix(M_).lu()
    assert np.allclose(np.array(L._rows) @ np.array(U._rows), M_[perm])                  # P A = L U
    worst = max(worst, float(np.max(np.abs(np.array(list(Matrix(M_).solve(rhs))) - np.linalg.solve(M_, rhs)))))
assert worst < 1e-8
assert raises(ValueError, lambda: Matrix([[1, 2], [2, 4]]).solve([1, 2])) and raises(ValueError, lambda: Matrix([[1, 2], [3]]))
Q = Matrix([[1, 2], [3, 4]])
alias = Q
Q += Matrix([[10, 0], [0, 10]])
assert alias is Q and Q == Matrix([[11, 2], [3, 14]]) and Q[1, 0] == 3.0 and Q[1] == Vec(3, 14) and Q.T[0, 1] == 3.0
Q[0, 0] = 5
assert alias[0, 0] == 5.0 and raises(TypeError, lambda: {Q: 1}) and Matrix.__hash__ is None
assert raises(ValueError, lambda: Q.__iadd__(Matrix([[1]])))
va = vb = Vec(1, 2)
va += Vec(1, 1)                                      # Vec has no __iadd__: a new object, vb unchanged
assert vb == Vec(1, 2) and va == Vec(2, 3) and va is not vb

sizes, t_py, t_np = (30, 60, 120), [], []
for n in sizes:
    P_, R_ = rng.normal(size=(n, n)), rng.normal(size=(n, n))
    t0 = time.perf_counter()
    C = Matrix(P_) @ Matrix(R_)
    t_py.append(time.perf_counter() - t0)
    t0 = time.perf_counter()
    for _ in range(20):
        C_np = P_ @ R_
    t_np.append((time.perf_counter() - t0) / 20)
    assert np.allclose(np.array(C._rows), C_np)
slope = np.polyfit(np.log(sizes), np.log(t_py), 1)[0]
print(f"  pure-Python matmul, n = {sizes}: {[f'{t * 1e3:.1f} ms' for t in t_py]}, log-log slope {slope:.2f}; "
      f"NumPy is {t_py[-1] / t_np[-1]:.0f}x faster at n = {sizes[-1]}")
assert 2.0 < slope < 3.8 and t_py[-1] / t_np[-1] > 20
print(f"Part 4: Matrix ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 5: dataclasses, inheritance, protocols, `__slots__` (notes §3.9-§3.11, Examples 4.6)

# %%
m = Measurement(21.5, "degC", 3.0)
assert str(m) == "Measurement(value=21.5, unit='degC', time=3.0)" and m == Measurement(21.5, "degC", 3.0)
assert hash(m) == hash(Measurement(21.5, "degC", 3.0)) and replace(m, value=22.0).value == 22.0
assert raises(FrozenInstanceError, lambda: setattr(m, "value", 0.0))
readings = [Measurement(21.5, "degC", 3.0), Measurement(19.0, "degC", 1.0),
            Measurement(21.5, "degC", 0.5), Measurement(20.2, "degC", 2.0)]
assert [r.time for r in sorted(readings)] == [1.0, 2.0, 0.5, 3.0] and max(readings).time == 3.0
assert [r.time for r in sorted(readings, key=operator.attrgetter("time"))] == [0.5, 1.0, 2.0, 3.0]
try:
    @dataclass
    class BadRun:
        name: str
        readings: list = []
except ValueError as err:
    assert str(err) == "mutable default <class 'list'> for field readings is not allowed: use default_factory"


@dataclass
class Run:
    name: str
    readings: list[Measurement] = field(default_factory=list)


r1, r2 = Run("a"), Run("b")
r1.readings.append(m)
assert r2.readings == [] and r1 == Run("a", [m]) and Run.__hash__ is None

pol = Polynomial((1, -3, 0, 2, 0))
assert pol.coeffs == (1.0, -3.0, 0.0, 2.0) and pol(2.0) == 11.0
assert pol.deriv() == Polynomial((-3, 0, 6)) and pol.deriv()(2.0) == 21.0
assert pol(Dual(2.0, 1.0)) == Dual(11.0, 21.0) and isinstance(pol, Callable)
for _ in range(50):
    c, x_ = rng.normal(size=int(rng.integers(1, 7))), float(rng.normal())
    assert math.isclose(Polynomial(tuple(c))(x_), NPPoly(c)(x_), rel_tol=1e-9, abs_tol=1e-12)
    assert math.isclose(Polynomial(tuple(c)).deriv()(x_), NPPoly(c).deriv()(x_), rel_tol=1e-9, abs_tol=1e-12)
    assert math.isclose(Polynomial(tuple(c))(Dual(x_, 1.0)).b, NPPoly(c).deriv()(x_), rel_tol=1e-9, abs_tol=1e-12)


class Integrator(ABC):
    def __init__(self, n):
        self.n = n

    def integrate(self, f, a, b):
        h = (b - a) / self.n
        return self._combine([f(a + i * h) for i in range(self.n + 1)], h)

    @abstractmethod
    def _combine(self, ys, h): ...


class Trapezoid(Integrator):
    def _combine(self, ys, h):
        return h * (sum(ys) - 0.5 * (ys[0] + ys[-1]))


class Simpson(Integrator):
    def _combine(self, ys, h):
        return h / 3 * (ys[0] + ys[-1] + 4 * sum(ys[1:-1:2]) + 2 * sum(ys[2:-1:2]))


trapezoid, simpson = Trapezoid(2)._combine, Simpson(2)._combine         # the rules as plain functions


def integrate(f, a, b, n, rule):
    h = (b - a) / n
    return rule([f(a + i * h) for i in range(n + 1)], h)


assert raises(TypeError, lambda: Integrator(8))
for cls, rule, ratio, first in ((Trapezoid, trapezoid, 4.0, 2.58e-2), (Simpson, simpson, 16.0, 2.69e-4)):
    e = [abs(cls(n).integrate(math.sin, 0.0, math.pi) - 2.0) for n in (8, 16, 32)]
    assert all(cls(n).integrate(math.sin, 0.0, math.pi) == integrate(math.sin, 0.0, math.pi, n, rule) for n in (8, 16, 32))
    assert abs(e[0] / first - 1) < 0.01 and all(abs(e[i] / e[i + 1] / ratio - 1) < 0.05 for i in (0, 1))


class CountingDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.writes = 0

    def __setitem__(self, key, value):
        self.writes += 1
        super().__setitem__(key, value)


class CountingUserDict(UserDict):
    writes = 0

    def __setitem__(self, key, value):
        self.writes += 1
        super().__setitem__(key, value)


for cls, expected in ((CountingDict, 1), (CountingUserDict, 4)):
    d = cls()
    d["a"] = 1
    d.update(b=2, c=3)
    d.setdefault("e", 5)
    assert len(d) == 4 and d.writes == expected


class Base:                                       # a diamond: Both(Left2, Right2) over Base
    def path(self):
        return ["Base"]


class Left2(Base):
    def path(self):
        return ["Left"] + super().path()


class Right2(Base):
    def path(self):
        return ["Right"] + super().path()


class Both(Left2, Right2):
    def path(self):
        return ["Both"] + super().path()


assert Both().path() == ["Both", "Left", "Right", "Base"]
assert [c.__name__ for c in Both.__mro__] == ["Both", "Left2", "Right2", "Base", "object"]


class Plain:
    def __init__(self, x, y):
        self.x, self.y = x, y


class Slotted:
    __slots__ = ("x", "y")

    def __init__(self, x, y):
        self.x, self.y = x, y


def bytes_per_instance(cls, n=100_000):
    tracemalloc.start()
    objs = [cls(1.0, 2.0) for _ in range(n)]
    current, _ = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return current / len(objs)


bp, bs = bytes_per_instance(Plain), bytes_per_instance(Slotted)
print(f"  bytes per instance: plain {bp:.0f}, slotted {bs:.0f}")
assert bs < 0.8 * bp and hasattr(Plain(1, 2), "__dict__") and not hasattr(Slotted(1, 2), "__dict__")
print(f"Part 5: dataclasses, ABCs, protocols ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Your turn
# Reference solutions are in `solutions.md` (B1-B4); try first.
#
# 1. **Slices.** Give `Vec` a `__getitem__`: a float for an `int`, a `Vec` for a `slice` (`v[1:]`).
#    Then check `np.asarray(v).shape` becomes `(3,)` (notes §3.10).
# 2. **`Complex`.** Write `Complex(re, im)` with `+ - * /`, `abs`, `conjugate`, `==`, `hash`, `repr`, and the
#    reflected operators needed for `2 * z`. Compare with built-in `complex` on 1000 random pairs.
# 3. **Real and dual exponents.** Let `Dual ** p` accept a float `p` (base `a > 0`), a `Dual` `p` (`x ** x`),
#    and add `__rpow__` (`2 ** x`). Check against SymPy.
# 4. **Second derivatives.** Differentiate twice by nesting: seed `Dual(Dual(x, 1), Dual(1, 0))`.
#    Which entry of the result holds `f''(x)`? Check `x * exp(x)` against SymPy.
# 5. **Exercise B7.** Add an eleventh function to `FUNCS` in Part 3 and draw its error-versus-`h` curve.

# %%
print(f"All checks passed in {time.perf_counter() - T0:.1f} s.")
