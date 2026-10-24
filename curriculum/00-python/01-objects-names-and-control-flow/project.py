# %% [markdown]
# # PY-1 Project: A Floating-Point Lab
#
# Rebuilds the numerical tools of PY-1 and checks every number quoted in the notes:
# the floor rule (3.3), the grid of doubles (3.4), `0.1 + 0.2` as a tie (Example 4.3),
# cancellation in the quadratic formula (3.5), Fast2Sum and Kahan summation (3.6),
# bisection, Newton and Horner (3.12), and a formatted report (Example 4.4).
# Standard library plus NumPy for three comparisons; no figures in this module.
# Style: PEP 8 with a 100-character line limit. Run: `MPLBACKEND=Agg python3 project.py`.
# Build first, then compare: write your own epsilon loop, kahan_sum, bisect and horner.

# %%
import copy
import math
import random
import struct
import sys
from decimal import Decimal, localcontext
from fractions import Fraction

import numpy as np

rng = random.Random(0)


def raises(exc, func, *args):
    """Return True if func(*args) raises exc."""
    try:
        func(*args)
    except exc:
        return True
    return False


def assign0(obj):
    """Try obj[0] = 9 (tuples and strings refuse)."""
    obj[0] = 9


# %% [markdown]
# ## Part 1: Integers and the floor rule (notes 3.3, eq. 3.1)

# %%
for a in range(-30, 31):
    for b in (-7, -3, -2, -1, 1, 2, 3, 7):
        q, r = a // b, a % b
        assert a == b * q + r                       # eq. (3.1)
        assert r == 0 or (r > 0) == (b > 0)         # r has the sign of the divisor
        assert abs(r) < abs(b)
assert (-7 // 2, -7 % 2, 7 // -2, 7 % -2) == (-4, 1, -4, -1)
assert divmod(-7, 2) == (-4, 1) and math.fmod(-7, 2) == -1.0
assert (12 & 10, 12 | 10, 12 ^ 10, ~12, -7 >> 1, 40 & 39, 32 & 31) == (8, 14, 6, -13, -4, 32, 0)
assert 10**20 + 1 == 100000000000000000001 and float(10**20) + 1 == 1e20

# %% [markdown]
# ## Part 2: The grid of doubles (notes 3.4, eqs. 3.2-3.4)


# %%
def machine_epsilon():
    """Gap between 1.0 and the next double, found by halving (sequence 2)."""
    eps = 1.0
    while 1.0 + eps / 2 != 1.0:
        eps /= 2
    return eps


EPS = machine_epsilon()
U = EPS / 2                                         # unit roundoff, 2**-53
assert EPS == sys.float_info.epsilon == 2.0**-52 == math.ulp(1.0) == np.finfo(float).eps
assert 1 + 2**-52 != 1 and 1 + 2**-53 == 1          # the tie goes to 1.0 (even)
assert (1 + 2**-53) + 2**-53 == 1 and 1 + (2**-53 + 2**-53) != 1
assert float(2**53 + 1) == float(2**53) and 2.0**53 + 1 == 2.0**53
assert math.nextafter(1.0, 2.0) == 1.0000000000000002 and math.ulp(0.0) == 5e-324
assert math.nextafter(1.0, 0.0) == 0.9999999999999999

print(f"{'x':>10} {'ulp(x)':>12} {'ulp(x)/x':>12}")
for x in (1e-10, 1e-3, 1.0, 1e3, 1e8, 1e16):
    spacing = math.ulp(x)
    print(f"{x:>10.1e} {spacing:>12.3e} {spacing / x:>12.3e}")
    assert spacing == 2.0 ** (math.frexp(x)[1] - 1 - 52)    # eq. (3.3)
    assert U < spacing / x <= EPS                            # relative spacing band

worst = 0.0                                          # eq. (3.4) against exact arithmetic
for _ in range(2000):
    q = Fraction(rng.randint(1, 10**9), rng.randint(1, 10**9)) * Fraction(10) ** rng.randint(-20, 20)
    delta = abs(Fraction(float(q)) - q) / q          # float(Fraction) is correctly rounded
    assert delta <= Fraction(1, 2**53)
    worst = max(worst, float(delta))
print(f"largest relative rounding error, 2000 random rationals: {worst:.3e} (u = {U:.3e})")


def fields(x):
    """Return (sign, exponent field E, fraction field F) of the double x."""
    (n,) = struct.unpack(">Q", struct.pack(">d", x))
    return n >> 63, (n >> 52) & 0x7FF, n & ((1 << 52) - 1)


def decode(x):
    """Rebuild a normal double exactly from its fields, eq. (3.2)."""
    s, e, f = fields(x)
    assert 1 <= e <= 2046
    return Fraction(-1) ** s * (1 + Fraction(f, 2**52)) * Fraction(2) ** (e - 1023)


assert all(decode(x) == Fraction(x) for x in (1.0, 0.5, -2.5, 0.1, math.pi, 1e300, 1e-300))
assert fields(1.0) == (0, 1023, 0) and fields(0.5)[1] == 1022 and fields(-2.5) == (1, 1024, 2**50)
assert fields(0.1)[1] == 1019 and f"{fields(0.1)[2]:052b}" == "1001" * 12 + "1010"
assert fields(5e-324) == (0, 0, 1) and fields(float("inf")) == (0, 2047, 0)
assert fields(float("nan"))[1:] == (2047, 2**51)

# %% [markdown]
# ## Part 3: Why 0.1 + 0.2 != 0.3 (Example 4.3), comparison, inf and nan

# %%
assert str(Decimal(0.1)) == "0.1000000000000000055511151231257827021181583404541015625"
assert str(Decimal(0.2)) == "0.200000000000000011102230246251565404236316680908203125"
assert str(Decimal(0.3)) == "0.299999999999999988897769753748434595763683319091796875"
assert str(Decimal(0.1 + 0.2)) == "0.3000000000000000444089209850062616169452667236328125"
assert (0.1).hex() == "0x1.999999999999ap-4" and (0.2).hex() == "0x1.999999999999ap-3"
assert (0.3).hex() == "0x1.3333333333333p-2" and (0.1 + 0.2).hex() == "0x1.3333333333334p-2"
assert (0.1).as_integer_ratio() == (3602879701896397, 2**55)
lo, hi = Fraction(0.3), Fraction(0.1 + 0.2)
assert hi - lo == Fraction(math.ulp(0.3))            # computed sum is one ulp above 0.3
assert Fraction(0.1) + Fraction(0.2) == (lo + hi) / 2    # ... because the exact sum is a tie
assert 0.1 + 0.2 != 0.3 and math.isclose(0.1 + 0.2, 0.3) and math.isclose(1e10 + 1, 1e10)
assert not math.isclose(1e-20, 0.0) and math.isclose(1e-20, 0.0, abs_tol=1e-12)
assert round(2.675, 2) == 2.67 and Decimal(2.675) < Decimal("2.675")

inf, nan = float("inf"), float("nan")
assert 1e308 * 10 == inf and 5e-324 / 2 == 0.0 and math.isnan(inf - inf)
assert nan != nan and not nan == nan and not nan < 1 and math.isnan(nan)
assert nan in [nan] and float("nan") not in [float("nan")]
assert raises(OverflowError, lambda: 2.0**10000) and raises(OverflowError, math.exp, 1000)
assert raises(ZeroDivisionError, lambda: 1 / 0) and raises(ValueError, math.sqrt, -1)

# %% [markdown]
# ## Part 4: Catastrophic cancellation, the quadratic formula (notes 3.5)


# %%
def roots_naive(a, b, c):
    """Roots of a x^2 + b x + c by the textbook formula."""
    d = math.sqrt(b * b - 4 * a * c)
    return (-b - d) / (2 * a), (-b + d) / (2 * a)


def roots_stable(a, b, c):
    """Same roots without cancellation, eq. (3.7): q = -(b + sign(b) d) / 2."""
    d = math.sqrt(b * b - 4 * a * c)
    q = -(b + math.copysign(d, b)) / 2
    return q / a, c / q


def roots_exact(a, b, c):
    """Reference roots from 60-digit decimal arithmetic."""
    with localcontext() as ctx:
        ctx.prec = 60
        a, b, c = Decimal(a), Decimal(b), Decimal(c)
        d = (b * b - 4 * a * c).sqrt()
        return (-b - d) / (2 * a), (-b + d) / (2 * a)


def rel_err(x, ref):
    """Relative error of the float x against the Decimal ref."""
    with localcontext() as ctx:
        ctx.prec = 60
        return float(abs((Decimal(x) - ref) / ref))


print(f"{'b':>8} {'naive small':>12} {'stable small':>13} {'naive big':>10} {'stable big':>11}")
for k in range(2, 11):
    b = 10.0**k
    big_ref, small_ref = roots_exact(1.0, b, 1.0)
    (nb, ns), (sb, ss) = roots_naive(1.0, b, 1.0), roots_stable(1.0, b, 1.0)
    e_ns, e_ss = rel_err(ns, small_ref), rel_err(ss, small_ref)
    e_nb, e_sb = rel_err(nb, big_ref), rel_err(sb, big_ref)
    print(f"{b:>8.0e} {e_ns:>12.2e} {e_ss:>13.2e} {e_nb:>10.2e} {e_sb:>11.2e}")
    assert max(e_ss, e_nb, e_sb) < 1e-15            # stable: full accuracy; big root always fine
    if k == 8:
        assert e_ns > 1e-8 and math.isclose(e_ns, 0.255, abs_tol=0.001)
        assert (ns, ss, nb, sb) == (-7.450580596923828e-09, -1e-08, -1e8, -1e8)
        d_true = (Decimal(b) ** 2 - 4).sqrt()       # bound (3.6) with x = b, y = d
        kappa = float((Decimal(b) + d_true) / (Decimal(b) - d_true))
        print(f"b=1e8: amplification factor (3.6) = {kappa:.2e}; kappa*u = {kappa * U:.3f}")
        assert e_ns <= kappa * U + U
assert math.sqrt(1e16 - 4) == 99999999.99999999 and math.ulp(1e8) == 2.0**-26

# %% [markdown]
# ## Part 5: Summation error, Fast2Sum, Kahan (notes 3.6)


# %%
def naive_sum(xs):
    """Left-to-right accumulation."""
    total = 0.0
    for x in xs:
        total += x
    return total


def kahan_sum(xs):
    """Compensated summation; c carries minus the rounding error of the last addition."""
    total = c = 0.0
    for x in xs:
        y = x - c
        t = total + y
        c = (t - total) - y
        total = t
    return total


def fast_two_sum(s, y):
    """Return (t, e) with t = fl(s + y) and s + y == t + e exactly, if |s| >= |y|; eq. (3.8)."""
    t = s + y
    return t, y - (t - s)


for _ in range(20000):
    p, q = (rng.uniform(-1, 1) * 10.0 ** rng.randint(-8, 8) for _ in range(2))
    s, y = (p, q) if abs(p) >= abs(q) else (q, p)
    t, e = fast_two_sum(s, y)
    assert Fraction(s) + Fraction(y) == Fraction(t) + Fraction(e)

N = 10**7
xs = [0.1] * N
ref = math.fsum(xs)
naive, kahan = naive_sum(xs), kahan_sum(xs)
err_naive, err_kahan = abs(naive - ref), abs(kahan - ref)
err_np = abs(float(np.sum(np.full(N, 0.1))) - ref)   # NumPy sums pairwise
print(f"fsum {ref!r}; naive {naive!r}; errors: naive {err_naive:.3e}, kahan {err_kahan:.3e}, "
      f"numpy {err_np:.3e}")
assert ref == 1000000.0 and math.ulp(1e6) == 2.0**-33
assert math.isclose(err_naive, 1.610e-4, rel_tol=1e-2) and 1e-5 < err_naive < 1e-3
assert err_kahan < 1e-12                             # measured: exactly 0.0
assert err_np < err_naive / 100
if sys.version_info >= (3, 12):
    assert abs(sum(xs) - ref) < 1e-9                 # built-in sum is compensated from 3.12
else:
    assert sum(xs) == naive

# %% [markdown]
# ## Part 6: Root finders (notes 3.12, eqs. 3.9-3.11)


# %%
def bisect(f, a, b, tol=1e-12, max_iter=200):
    """Root of f in [a, b] by bisection; returns (midpoint, iterations). Count: eq. (3.9)."""
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


def newton(f, df, x0, tol=1e-12, max_iter=50):
    """Root of f near x0 by x <- x - f(x)/f'(x), eq. (3.10); returns (x, iterations)."""
    x = x0
    for k in range(1, max_iter + 1):
        step = f(x) / df(x)
        x -= step
        if abs(step) <= tol:
            return x, k
    raise RuntimeError("Newton's method did not converge")


# name, f, f', f'', bracket, Newton start, true root (nearest double)
CASES = [
    ("cos x - x", lambda x: math.cos(x) - x, lambda x: -math.sin(x) - 1, lambda x: -math.cos(x),
     (0.0, 1.0), 1.0, 0.7390851332151607),
    ("x^2 - 2", lambda x: x * x - 2, lambda x: 2 * x, lambda x: 2.0, (1.0, 2.0), 1.5, math.sqrt(2)),
    ("x^3 - 2x - 5", lambda x: x**3 - 2 * x - 5, lambda x: 3 * x * x - 2, lambda x: 6 * x,
     (2.0, 3.0), 2.0, 2.0945514815423265),
]
TOL = 1e-10
print(f"{'function':>14} {'bisect k':>9} {'newton k':>9}   Newton errors e_n (then e_{{n+1}}/e_n^2 -> C)")
for name, f, df, d2f, (a, b), x0, root in CASES:
    r_b, k_b = bisect(f, a, b, tol=TOL)
    assert k_b == math.ceil(math.log2((b - a) / TOL)) and abs(r_b - root) <= TOL   # eq. (3.9)
    r_n, k_n = newton(f, df, x0, tol=1e-14)
    assert abs(r_n - root) < 1e-14
    errs, x = [], x0
    for _ in range(5):
        errs.append(abs(x - root))
        x -= f(x) / df(x)
    C = abs(d2f(root) / (2 * df(root)))                           # constant of eq. (3.11)
    ratios = [errs[i + 1] / errs[i] ** 2 for i in range(4) if 1e-6 < errs[i] < 0.05]
    assert ratios and all(abs(r - C) < 0.25 * C for r in ratios)
    print(f"{name:>14} {k_b:>9} {k_n:>9}   " + " ".join(f"{e:.1e}" for e in errs) + f"   C = {C:.3f}")

f_cos, df_cos = CASES[0][1], CASES[0][2]
assert bisect(f_cos, 0.0, 1.0, tol=1e-6) == (0.7390851974487305, 20)
calls = []
bisect(lambda x: calls.append(x) or f_cos(x), 0.0, 1.0, tol=1e-3)   # Example 4.5
mids = calls[2:]                                     # calls[0], calls[1] are f(a), f(b)
assert len(mids) == 10 == math.ceil(math.log2(1 / 1e-3)) and mids[:3] == [0.5, 0.75, 0.625]
assert f"{mids[-1]:.6f}" == "0.739258" and f"{f_cos(mids[1]):.6f}" == "-0.018311"
x, errors = 1.0, []                                  # the Newton demo of notes 3.12
for _ in range(5):
    errors.append(abs(x - 0.7390851332151607))
    x -= f_cos(x) / df_cos(x)
assert [f"{e:.2e}" for e in errors] == ["2.61e-01", "1.13e-02", "2.78e-05", "1.70e-10", "0.00e+00"]
assert raises(RuntimeError, bisect, lambda x: x * x - 2, 1.0, 2.0, 1e-20)   # tol below the grid
assert raises(ValueError, bisect, lambda x: x * x + 1, -1.0, 1.0)
assert raises(ZeroDivisionError, newton, lambda x: x * x - 2, lambda x: 2 * x, 0.0)
true_root = Decimal("0.73908513321516064165531208767387")
assert abs(Decimal(0.7390851332151607) - true_root) < Decimal(math.ulp(0.74)) / 2   # nearest double

# %% [markdown]
# ## Part 7: Horner's rule and operation counts (notes 3.12, eq. 3.12)


# %%
def poly_naive(coeffs, x):
    """coeffs = [a_n, ..., a_0], term by term; returns (value, multiplications, additions)."""
    n = len(coeffs) - 1
    total, mults, adds = coeffs[-1], 0, 0
    for k in range(1, n + 1):
        power = x
        for _ in range(k - 1):
            power *= x
            mults += 1
        total += coeffs[n - k] * power
        mults += 1
        adds += 1
    return total, mults, adds


def horner(coeffs, x):
    """Nested form (3.12); returns (value, multiplications, additions)."""
    result, mults, adds = coeffs[0], 0, 0
    for a in coeffs[1:]:
        result = result * x + a
        mults += 1
        adds += 1
    return result, mults, adds


assert horner([2, -6, 2, -1], 3)[0] == poly_naive([2, -6, 2, -1], 3)[0] == 5
for n in (3, 10, 100):
    coeffs = [1.0 / (k + 1) for k in range(n + 1)]
    (v1, m1, a1), (v2, m2, a2) = poly_naive(coeffs, 1.01), horner(coeffs, 1.01)
    assert (m1, m2, a1, a2) == (n * (n + 1) // 2, n, n, n) and math.isclose(v1, v2, rel_tol=1e-12)
    print(f"n={n:>3}: naive {m1:>4} multiplications, Horner {m2:>3}")
coeffs = [1.0 / (k + 1) for k in range(21)]
assert math.isclose(float(np.polyval(coeffs, 1.01)), horner(coeffs, 1.01)[0], rel_tol=1e-13)

# %% [markdown]
# ## Part 8: Formatting, the report generator (Example 4.4)


# %%
def format_table(rows):
    """Return the table of Example 4.4 as a list of lines."""
    lines = [f"{'quantity':<20}{'fixed ,.4f':>38}{'sci .4e':>13}{'general .6g':>14}", "-" * 85]
    for name, value in rows:
        lines.append(f"{name:<20}{value:>38,.4f}{value:>13.4e}{value:>14.6g}")
    return lines


ROWS = [("pi", math.pi), ("speed of light, m/s", 299_792_458), ("Avogadro, 1/mol", 6.02214076e23),
        ("Planck, J s", 6.62607015e-34), ("machine epsilon", sys.float_info.epsilon)]
table = format_table(ROWS)
print("\n".join(table))
assert [len(line) for line in table] == [85] * 7
assert [line[20:].split() for line in table[2:]] == [
    ["3.1416", "3.1416e+00", "3.14159"],
    ["299,792,458.0000", "2.9979e+08", "2.99792e+08"],
    ["602,214,075,999,999,987,023,872.0000", "6.0221e+23", "6.02214e+23"],
    ["0.0000", "6.6261e-34", "6.62607e-34"],
    ["0.0000", "2.2204e-16", "2.22045e-16"],
]
x = 1234.56789
assert f"{x:.2f}|{x:10.2f}|{x:<10.2f}|{x:^10.2f}|{x:010.2f}|{x:+.1f}" == \
    "1234.57|   1234.57|1234.57   | 1234.57  |0001234.57|+1234.6"
assert f"{x:,.1f}|{x:.3e}|{x:.3g}|{0.256:.1%}|{255:#x}|{5:08b}" == \
    "1,234.6|1.235e+03|1.23e+03|25.6%|0xff|00000101"
assert f"{x=:.1f}|{'ab':*^8}" == "x=1234.6|***ab***"

# %% [markdown]
# ## Part 9: The language examples (notes 3.2, 3.7-3.11 and Examples 4.1, 4.2, 4.6)

# %%
QUIZ = [                                             # sequence 1: (setup, expression, value)
    ("a = [1, 2]; b = a; b.append(3)", "a", [1, 2, 3]),
    ("a = [1, 2]; b = a; b = b + [3]", "a", [1, 2]),
    ("a = [1, 2]; b = a; b += [3]", "a", [1, 2, 3]),
    ("x = 5; y = x; y += 1", "(x, y)", (5, 6)),
    ("grid = [[0] * 2] * 2; grid[0][0] = 9", "grid", [[9, 0], [9, 0]]),
    ("def grow(i): i.append(1)\ndef regrow(i): i = i + [1]\nlst = []; grow(lst); regrow(lst)", "lst", [1]),
    ("t = (1, [2, 3]); t[1].append(4)", "t", (1, [2, 3, 4])),
    ("def collect(x, bucket=[]):\n    bucket.append(x)\n    return bucket\ncollect(1)",
     "(collect(2), collect.__defaults__)", ([1, 2], ([1, 2],))),
    ("counter = 0\ndef bump_global():\n    global counter\n    counter += 1\n"
     "bump_global(); bump_global()", "counter", 2),
    ("def make_counter():\n    count = 0\n    def step():\n        nonlocal count\n"
     "        count += 1\n        return count\n    return step\ns = make_counter()", "(s(), s(), s())", (1, 2, 3)),
    ("total = 0\nfor k in range(1, 10):\n    if k % 2 == 0:\n        continue\n    total += k * k", "total", 165),
    ("n, steps = 27, 0\nwhile n != 1:\n    n = 3 * n + 1 if n % 2 else n // 2\n    steps += 1", "steps", 111),
    ("primes = []\nfor n in range(2, 30):\n    for d in range(2, n):\n        if n % d == 0:\n"
     "            break\n    else:\n        primes.append(n)", "primes", [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]),
]
for setup, expr, value in QUIZ:
    ns = {}
    exec(setup, ns)
    assert eval(expr, ns) == value, setup
assert raises(TypeError, assign0, (1, [2])) and raises(TypeError, assign0, "abc")
x = 10


def bump():
    x = x + 1                                        # noqa: F823  (UnboundLocalError is the point)
    return x


assert raises(UnboundLocalError, bump) and x == 10
assert (0 or 10, 5 or 10, 0 and 10, 5 and 10) == (10, 5, 0, 10)
assert not any([0, 0.0, "", [], {}, None]) and all([0.1, " ", [0], "False", nan])
z = 3 + 4j
assert (z.real, z.imag, abs(z), z.conjugate(), z * 1j) == (3.0, 4.0, 5.0, 3 - 4j, -4 + 3j)
assert __import__("cmath").exp(1j * math.pi) == -1 + 1.2246467991473532e-16j
assert abs(((-8) ** (1 / 3)).imag - 1.7320508075688772) < 1e-15
s = "Hello, Python"
assert (len(s), s[0], s[-1], s[7:], s[:5], s[::-1]) == (13, "H", "n", "Python", "Hello", "nohtyP ,olleH")
assert s.upper() == "HELLO, PYTHON" and s == "Hello, Python"
assert list(range(8, 3, -2)) == [8, 6, 4] and len(range(3, 8)) == 5

a = [1, 2]                                           # Example 4.2
b = a
c = list(a)
assert a is b and a is not c and a == c and a[0] is c[0]
a.append(3)
assert (a, b, c) == ([1, 2, 3], [1, 2, 3], [1, 2])
m = [[1, 2], [3]]
shallow, deep = list(m), copy.deepcopy(m)
m[0].append(99)
assert shallow == [[1, 2, 99], [3]] and deep == [[1, 2], [3]]


def parse(line):                                     # Example 4.1: the traceback program
    values = []
    for token in line.split(","):
        values.append(float(token))
    return values


def summarize(lines):
    return [sum(parse(line)) / len(parse(line)) for line in lines]


assert summarize(["1.0,2.0,3.0", "4.5,5.5"]) == [2.0, 5.0]
try:
    summarize("1.0,2.0,3.0")                         # a str, so it iterates over characters
except ValueError as err:
    assert str(err) == "could not convert string to float: '.'"
else:
    raise AssertionError("expected ValueError")


def record(value, log=[]):                           # Example 4.6
    log.append(value)
    return log


def record_fixed(value, log=None):
    if log is None:
        log = []
    log.append(value)
    return log


run1, run2 = record(1.5), record(2.5)
assert run1 is run2 and run1 == [1.5, 2.5] and record.__defaults__ == ([1.5, 2.5],)
run3, run4 = record_fixed(1.5), record_fixed(2.5)
mine = []
assert run3 is not run4 and record_fixed.__defaults__ == (None,)
assert record_fixed(7.0, mine) is mine and mine == [7.0]

# %% [markdown]
# ## Your turn
#
# 1. **Secant method (exercise B3).** After `newton` in Part 6, add `secant(f, x0, x1, tol, max_iter)`
#    that replaces f'(x_n) by the slope through the last two points. Compare iteration counts
#    with Newton and bisection on the three CASES, and estimate its order of convergence.
# 2. **Series with a stopping rule (exercise B5).** Sum 1/k^2 in a `while` loop until a term is
#    below 1e-8. Compare the error with pi^2/6. When does "stop when the next term is tiny"
#    stop too early, and for which series is it safe?
# 3. **Largest n with `float(10**n) == 10**n`.** Predict it from 5**n < 2**53, then run it.
# 4. **Extend the lab (exercise B7).** Add a cell that measures how the error of `naive_sum` on
#    0.1 grows with n = 10^3 ... 10^7. Solutions to all four are in solutions.md.

# %%
if __name__ == "__main__":
    print("All checks passed.")
