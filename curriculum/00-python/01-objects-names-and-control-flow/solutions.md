# PY-1 Solutions

Every code block below was run, and its output is pasted from that run. Blocks that end in asserts passed them. `roots.py` is the file from notes §3.12. Code lines may run to 100 columns here (PEP 8 allows a team limit; the notes use 88).

## A. Warm-up

### A1. Aliasing

```python
# file: a1.py
x = [1, 2, 3]
y = x
y[0] = 99
z = x[:]
z[1] = 77
print(x, y, z)


def add_one(lst):
    lst.append(1)
    return lst


a = [0]
b = add_one(a)
print(a is b, a)

s = "abc"
t = s
t += "d"
print(s, t)

d = {"k": [1]}
e = dict(d)
e["k"].append(2)
e["j"] = 0
print(d, e)
```

```text
[99, 2, 3] [99, 2, 3] [99, 77, 3]
True [0, 1]
abc abcd
{'k': [1, 2]} {'k': [1, 2], 'j': 0}
```

Line 1: `y = x` ties a second name to the list (R1), so writing through `y` shows in `x`; `x[:]` builds a new list, so `z[1] = 77` does not. Line 2: R3 ties `lst` to the caller's list, `append` mutates it (R2), and `return lst` hands back the same object, so `a is b`. Line 3: strings are immutable, so `t += "d"` builds a new string and rebinds `t`; `s` is untouched. Line 4: `dict(d)` is a shallow copy: `e` has its own keys (adding `"j"` does not touch `d`), but both dictionaries hold the same list under `"k"`, so the `append` shows in both.

### A2. Floor rule

```python
# file: a2.py
print(-17 // 5, -17 % 5, 17 // -5, 17 % -5, divmod(-17, 5))
print(2 ** -2, 7 // 2.0, -7.5 % 2, -3 ** 2, (-3) ** 2)


def hour_after(h, delta):
    """Clock hour in range(24) that is delta hours after hour h."""
    return (h + delta) % 24


assert hour_after(22, 5) == 3
assert hour_after(2, -5) == 21
assert hour_after(0, -48) == 0
assert all(0 <= hour_after(h, d) < 24 for h in range(24) for d in range(-100, 100))
```

```text
-4 3 -4 -3 (-4, 3)
0.25 3.0 0.5 -9 9
```

By the floor rule (3.1), `-17 // 5 = floor(-3.4) = -4` and the remainder is `-17 - 5(-4) = 3`; `17 // -5 = floor(-3.4) = -4` with remainder `17 - (-5)(-4) = -3`. `2 ** -2` is the float `0.25`; `7 // 2.0` is the float `3.0`; `-7.5 % 2` is `-7.5 - 2(-4) = 0.5`. Precedence: `-3 ** 2` is `-(3 ** 2)`. The one-liner works because for a positive divisor `a % 24` always lies in `range(24)`, even for negative `a`; a truncating remainder would give negative hours.

### A3. The grid

```python
# file: a3.py
print(0.1 + 0.1 + 0.1 == 0.3)
print(1e16 + 1 == 1e16)
print((1e16 + 1) + 1 == 1e16, 1e16 + (1 + 1) == 1e16)
print(0.5 + 0.25 == 0.75)
print(float("inf") == float("inf") + 1)
print(2**53 + 1 == float(2**53 + 1))
```

```text
False
True
True False
True
True
False
```

1. `False`: `0.1 + 0.1` is exactly `0.2`, but `0.2 + 0.1` ties between two doubles and rounds to the even one, `0.30000000000000004` (Example 4.3).
2. `True`: near $10^{16}$ the spacing is 2, so `1e16 + 1` is a tie, and round-half-even gives `1e16`.
3. `True, False`: adding `1` twice, one at a time, rounds back to `1e16` both times; adding `2` first lands on the next double. Floating-point addition is not associative.
4. `True`: $0.5$, $0.25$ and $0.75$ are sums of powers of two, exactly representable, so no rounding occurs.
5. `True`: `inf + 1` is `inf`.
6. `False`: the left side is the exact integer $2^{53} + 1$; the right side is that integer converted to a double, which rounds to $2^{53}$ (a tie, even). Comparing an `int` with a `float` in Python is exact, so they differ.

### A4. Scope

```python
# file: a4.py
total = 0


def add(values):
    for v in values:
        total += v
    return total


print(add([1, 2, 3]))
```

```text
Traceback (most recent call last):
  File "/tmp/py1/a4.py", line 10, in <module>
    print(add([1, 2, 3]))
          ~~~^^^^^^^^^^^
  File "/tmp/py1/a4.py", line 6, in add
    total += v
    ^^^^^
UnboundLocalError: cannot access local variable 'total' where it is not associated with a value
```

`total += v` assigns to `total`, so `total` is local to `add` for the whole body (§3.11); the first execution reads a local with no value: `UnboundLocalError`. Two fixes:

```python
# file: a4_fixed.py
total = 0


def add_global(values):
    global total              # rebinds the module-level name
    for v in values:
        total += v
    return total


def add_local(values):
    result = 0                # a local accumulator; nothing outside changes
    for v in values:
        result += v
    return result


assert add_local([1, 2, 3]) == 6
assert add_global([1, 2, 3]) == 6 and add_global([1, 2, 3]) == 12   # remembers state!
names = ["a", "b"]


def rename():
    names.append("c")         # mutation: no declaration needed


rename()
print(names, total)
```

```text
['a', 'b', 'c'] 12
```

Ship `add_local` (or the built-in `sum(values)`): `add_global` returns 12 on its second call, because it depends on hidden state, which is exactly why `global` is avoided. `rename` needs no `global` because `append` mutates the list that `names` points at (R2); no name is rebound.

### A5. The bits of a double

```python
# file: a5.py
import struct


def fields(x):
    """Return (sign, E, F) for the IEEE 754 double x, eq. (3.2)."""
    n = int.from_bytes(struct.pack(">d", x), "big")
    return n >> 63, (n >> 52) & 0x7FF, n & ((1 << 52) - 1)


for x in (1.0, 0.5, 8.0, 0.1, -0.0):
    s, e, f = fields(x)
    print(f"{x!r:>5}: sign={s} E={e:>4} F={f:013x}")
assert fields(1.0)[1] == 1023 and fields(0.5)[1] == 1022 and fields(8.0)[1] == 1026
assert fields(-0.0) == (1, 0, 0)
```

```text
  1.0: sign=0 E=1023 F=0000000000000
  0.5: sign=0 E=1022 F=0000000000000
  8.0: sign=0 E=1026 F=0000000000000
  0.1: sign=0 E=1019 F=999999999999a
 -0.0: sign=1 E=   0 F=0000000000000
```

The exponent field of $8.0 = 1.0 \times 2^3$ is $1023 + 3 = 1026$, so `E - 1023` is the power of two: $3$. `0.1` has $E = 1019$, $2^{-4}$, and a fraction (printed in hex) of repeating `9`s ending in `a`. `-0.0` is a real value with only the sign bit set; it compares equal to `0.0`.

### A6. Truthiness and formatting

```python
# file: a6.py
print(bool([]), bool([0]), bool(0.0), bool("0"), bool(None))
print(0 or "default", "" or 0 or [], [1] and "yes", None and 1 / 0)
print(f"{3.14159:8.3f}|{42:<6d}|{42:>6d}|{0.00012345:.2e}|{1234567:,}|{0.5:.0%}")


def connect(timeout=None):
    if timeout is None:       # only a *missing* argument gets the default
        timeout = 30
    return timeout


assert connect() == 30 and connect(timeout=0) == 0
```

```text
False True False True False
default [] yes None
   3.142|42    |    42|1.23e-04|1,234,567|50%
```

Only empty containers, zeros and `None` are falsy: `[0]` and `"0"` are non-empty, so truthy. `or` returns the first truthy operand or the last operand (`[]`); `and` returns the first falsy operand or the last, and short-circuits, so `1 / 0` never runs. The format line: width 8 with 3 decimals, left and right alignment, scientific with 2 decimals, thousands separators, and `.0%` multiplies by 100 and rounds. The bug in `timeout = timeout or 30`: `connect(timeout=0)` (a legitimate "do not wait") silently becomes `30`.

## B. Core

### B1. Aliasing audit

```python
# file: b1.py
import copy


def describe_sharing(a, b):
    """Describe how two lists relate: identity, equality, and shared elements."""
    info = {"same_object": a is b, "equal": a == b}
    if len(a) == len(b):
        info["shared"] = [i for i in range(len(a)) if a[i] is b[i]]
    return info


rows_bad = [[0] * 3] * 3
rows_good = [[0] * 3 for _ in range(3)]
print("distinct rows:", len({id(r) for r in rows_bad}), "vs", len({id(r) for r in rows_good}))
assert len({id(r) for r in rows_bad}) == 1 and len({id(r) for r in rows_good}) == 3

m = [[1, 2], [3, 4]]
for name, c in (("m itself", m), ("list(m)", list(m)), ("m[:]", m[:]), ("deepcopy", copy.deepcopy(m))):
    print(f"{name:>9}:", describe_sharing(m, c))
assert describe_sharing(m, m)["same_object"]
assert describe_sharing(m, list(m)) == {"same_object": False, "equal": True, "shared": [0, 1]}
assert describe_sharing(m, m[:]) == {"same_object": False, "equal": True, "shared": [0, 1]}
assert describe_sharing(m, copy.deepcopy(m)) == {"same_object": False, "equal": True, "shared": []}
```

```text
distinct rows: 1 vs 3
 m itself: {'same_object': True, 'equal': True, 'shared': [0, 1]}
  list(m): {'same_object': False, 'equal': True, 'shared': [0, 1]}
     m[:]: {'same_object': False, 'equal': True, 'shared': [0, 1]}
 deepcopy: {'same_object': False, 'equal': True, 'shared': []}
```

`[[0] * 3] * 3` has one inner row three times (§3.2 snippet 5). `list(m)` and `m[:]` are shallow: a new outer list, the same inner lists, so `shared` is `[0, 1]`. `deepcopy` shares nothing. **Rule:** a shallow copy is safe when the elements are immutable (numbers, strings, tuples of those) or when you will only replace elements, never mutate them; use `deepcopy` when you will mutate nested mutable elements and the original must not see it.

### B2. `isclose` from scratch

```python
# file: b2.py
import math
import random


def my_isclose(a, b, rel_tol=1e-9, abs_tol=0.0):
    """True if |a - b| <= max(rel_tol * max(|a|, |b|), abs_tol), eq. (3.5)."""
    if a == b:
        return True               # includes inf == inf
    if math.isinf(a) or math.isinf(b):
        return False              # one infinite, or opposite infinities
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


inf, nan = float("inf"), float("nan")
assert my_isclose(inf, inf) and not my_isclose(inf, 1e308) and not my_isclose(inf, -inf)
assert not my_isclose(nan, nan) and not my_isclose(nan, 1.0)
assert my_isclose(0.0, -0.0) and not my_isclose(1e-20, 0.0) and my_isclose(1e-20, 0.0, abs_tol=1e-12)

rng = random.Random(0)
mismatches = 0
for i in range(100_000):
    a = rng.uniform(-1, 1) * 10.0 ** rng.randint(-12, 12)
    if i % 2:
        b = a * (1 + rng.uniform(-3, 3) * 1e-9)       # relative gap within 3x of rel_tol
    else:
        b = rng.uniform(-1, 1) * 10.0 ** rng.randint(-12, 12)
    mismatches += my_isclose(a, b) != math.isclose(a, b)
    assert my_isclose(a, b) == my_isclose(b, a)       # symmetric
print("mismatches with math.isclose:", mismatches)
assert mismatches == 0

a, b = 1.0, 1.0 + 1e-12
print(my_isclose(a, b), my_isclose(a, b, rel_tol=0))
assert my_isclose(a, b) and not my_isclose(a, b, rel_tol=0)
```

```text
mismatches with math.isclose: 0
True False
```

The pair `1.0` and `1.0 + 1e-12` is close under the default relative tolerance and not under `rel_tol=0`, where only `abs_tol` (default `0`) remains. For voltages around $10^{-6}$ V, a relative tolerance alone is fine away from zero, but a reading of exactly `0.0` V can only be compared with an `abs_tol` from the instrument: its noise floor or resolution, say `1e-9` V (a thousandth of the typical value), never a number picked to make a test pass.

### B3. The secant method

```python
# file: b3.py
import math

from roots import bisect, newton


def secant(f, x0, x1, tol=1e-12, max_iter=50):
    """Root of f from two starting points; returns (x, iterations)."""
    f0, f1 = f(x0), f(x1)
    for k in range(1, max_iter + 1):
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(x2 - x1) <= tol:
            return x2, k
        x0, f0, x1, f1 = x1, f1, x2, f(x2)
    raise RuntimeError("secant method did not converge")


cases = [
    ("cos x - x", lambda x: math.cos(x) - x, lambda x: -math.sin(x) - 1, (0.0, 1.0), 1.0),
    ("x^3 - 2x - 5", lambda x: x**3 - 2 * x - 5, lambda x: 3 * x * x - 2, (2.0, 3.0), 2.0),
]
print(f"{'function':>14} {'bisect':>7} {'newton':>7} {'secant':>7}")
for name, f, df, (a, b), x0 in cases:
    kb = bisect(f, a, b, tol=1e-12)[1]
    kn = newton(f, df, x0, tol=1e-12)[1]
    rs, ks = secant(f, a, b, tol=1e-12)
    print(f"{name:>14} {kb:>7} {kn:>7} {ks:>7}")
    assert abs(f(rs)) < 1e-12

# order of convergence on x^2 - 2 from (1, 2)
f, root = (lambda x: x * x - 2), math.sqrt(2)
xs = [1.0, 2.0]
for _ in range(6):
    x0, x1 = xs[-2], xs[-1]
    xs.append(x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0)))
errs = [abs(x - root) for x in xs]
print("errors:", " ".join(f"{e:.1e}" for e in errs))
ps = [math.log(errs[i + 1] / errs[i]) / math.log(errs[i] / errs[i - 1])
      for i in range(3, len(errs) - 1) if errs[i + 1] > 1e-12]     # skip the transient and the floor
print("order estimates:", " ".join(f"{p:.2f}" for p in ps), f"(mean {sum(ps) / len(ps):.2f})")
assert 1.4 < sum(ps) / len(ps) < 1.9, ps
```

```text
      function  bisect  newton  secant
     cos x - x      40       5       6
  x^3 - 2x - 5      40       5       7
errors: 4.1e-01 5.9e-01 8.1e-02 1.4e-02 4.2e-04 2.1e-06 3.2e-10 2.2e-16
order estimates: 2.02 1.50 1.67 (mean 1.73)
```

Secant needs more iterations than Newton and far fewer than bisection. Its order is $(1+\sqrt5)/2 \approx 1.618$ (against Newton's 2); single-step estimates scatter between 1.5 and 2.0, and their mean above agrees with the theory. The right comparison is cost per *evaluation of f*: Newton needs `f` and `df` each step, secant one new `f`, so secant's effective rate per evaluation, $1.618$, beats Newton's $\sqrt2 \approx 1.414$ when `df` costs as much as `f`. If `f` takes a minute and no derivative formula exists, use secant (or `brentq`, PY-8, which adds a bracket).

### B4. Remainders in real use

```python
# file: b4.py
import random
from datetime import timedelta


def to_clock(seconds):
    """Split an integer number of seconds into (days, hours, minutes, secs), all but days >= 0."""
    minutes, secs = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    return days, hours, minutes, secs


rng = random.Random(0)
for _ in range(10_000):
    s = rng.randint(-10**6, 10**6)
    days, hours, minutes, secs = to_clock(s)
    td = timedelta(seconds=s)
    assert days == td.days and hours * 3600 + minutes * 60 + secs == td.seconds
    assert 0 <= hours < 24 and 0 <= minutes < 60 and 0 <= secs < 60
print(to_clock(-1), to_clock(90_061), to_clock(-90_061))
assert to_clock(-1) == (-1, 23, 59, 59)
```

```text
(-1, 23, 59, 59) (1, 1, 1, 1) (-2, 22, 58, 59)
```

Each `divmod` is the floor rule (3.1) with a positive divisor, so each remainder lies in `range(divisor)` whatever the sign of the input; the sign is carried by the quotient and ends up in `days`. A truncating remainder (`math.fmod(-1, 60)` is `-1.0`) would have returned `-1` seconds, and a negative `hours` and `minutes` would follow, which is not a clock reading.

### B5. A series with a stopping rule

```python
# file: b5.py
import math

# (a) e = sum 1/k!: stop when the next term is below 1e-12
total, term, k = 0.0, 1.0, 0
while term >= 1e-12:
    total += term
    k += 1
    term /= k
print(f"(a) e: {k} terms, error {abs(total - math.e):.2e}")
assert abs(total - math.e) < 2e-12

# (b) sum 1/k^2: stop when the next term is below 1e-8
total, k = 0.0, 1
while 1 / k**2 >= 1e-8:
    total += 1 / k**2
    k += 1
err_b = abs(total - math.pi**2 / 6)
print(f"(b) pi^2/6: {k - 1} terms, error {err_b:.2e}")
assert 0.5e-4 < err_b < 2e-4

# (c) stop when the tail bound 1/N, valid after N terms, is below 1e-6
total, n = 0.0, 0
while 1 / max(n, 1) >= 1e-6:
    n += 1
    total += 1 / n**2
err_c = abs(total - math.pi**2 / 6)
print(f"(c) tail-bound rule: {n} terms, error {err_c:.2e}")
assert err_c < 1e-6

# (d) e^-20 from its own series loses everything; the reciprocal of e^20 does not
def exp_series(x):
    total, term, k = 0.0, 1.0, 0
    while abs(term) >= 1e-17 * abs(total) or k < 1:
        total += term
        k += 1
        term *= x / k
    return total


direct = exp_series(-20.0)
recip = 1 / exp_series(20.0)
exact = math.exp(-20.0)
print(f"(d) series at -20: {direct:.6e}; 1/series at +20: {recip:.6e}; math.exp: {exact:.6e}")
print(f"    relative errors: {abs(direct - exact) / exact:.1e} and {abs(recip - exact) / exact:.1e}")
assert abs(direct - exact) / exact > 1e-3 and abs(recip - exact) / exact < 1e-14
```

```text
(a) e: 15 terms, error 8.15e-13
(b) pi^2/6: 10000 terms, error 1.00e-04
(c) tail-bound rule: 1000001 terms, error 1.00e-06
(d) series at -20: 6.147562e-09; 1/series at +20: 2.061154e-09; math.exp: 2.061154e-09
    relative errors: 2.0e+00 and 4.0e-16
```

(a) works because the tail after a term is about the next term times $1/(k+1)$, so "term below `1e-12`" means "error about `1e-12`". (b) fails because the terms of $1/k^2$ shrink slowly: after $N$ terms the tail is about $1/N$, which is $N$ times the last term. Stopping when the term is $10^{-8}$ ($N = 10^4$) leaves an error near $10^{-4}$. A term-size rule certifies the error only when the tail is bounded by a small multiple of the next term (geometric-type decay). (c) uses the bound $\sum_{k>N}1/k^2 < 1/N$ to stop on a quantity that really bounds the error. (d) is §3.5 in disguise: at $x = -20$ the terms reach $4\times10^7$ in magnitude and alternate in sign, each carrying an absolute rounding error near $4\times10^7 \times 10^{-16} \approx 4\times10^{-9}$, more than the answer ($2\times10^{-9}$). Summing the positive series for $e^{20}$ and taking the reciprocal avoids the subtraction.

### B6. A report with computed column widths

```python
# file: b6.py
data = [("copper", 8.96, 0.02), ("iron", 7.874, 0.004), ("aluminium", 2.70, 0.01), ("lead", 11.34, 0.05)]


def report(rows):
    """Return an aligned density table as a list of lines."""
    w = max(len(name) for name, _, _ in rows)
    lines = [f"{'material':<{w}} {'density':>8} {'± unc':>8}", "-" * (w + 18)]
    for name, rho, unc in rows:
        lines.append(f"{name:<{w}} {rho:>8.3f} {unc:>8.3f}")
    mean = sum(rho for _, rho, _ in rows) / len(rows)
    lines += ["-" * (w + 18), f"{'mean':<{w}} {mean:>8.3f}"]
    return lines


text = "\n".join(report(data))
print(text)
assert text == """\
material   density    ± unc
---------------------------
copper       8.960    0.020
iron         7.874    0.004
aluminium    2.700    0.010
lead        11.340    0.050
---------------------------
mean         7.719"""
```

```text
material   density    ± unc
---------------------------
copper       8.960    0.020
iron         7.874    0.004
aluminium    2.700    0.010
lead        11.340    0.050
---------------------------
mean         7.719
```

The nested field `{w}` inside the spec makes the width a variable, so the name column fits the longest name. The mean of the decimal inputs is exactly $7.7185$, a tie at three decimals; the doubles involved make the computed value `7.718500000000001` (the reverse order and `math.fsum` agree), a hair above the tie, so it prints as `7.719`. Decimal values that sit on a formatting tie behave like `round(2.675, 2)`: the stored double decides.

### B7. Extending `project.py`

Add this cell after Part 5 (it uses `naive_sum`, `math` and `sys`, all defined above; the block below is self-contained so it can be run alone):

```python
# file: b7.py
import math


def naive_sum(xs):
    total = 0.0
    for x in xs:
        total += x
    return total


rows = []
print(f"{'n':>10} {'abs error':>12} {'bound n*u*sum|x|':>17}")
for k in range(3, 8):
    n = 10**k
    xs = [0.1] * n
    err = abs(naive_sum(xs) - math.fsum(xs))
    bound = n * 2.0**-53 * math.fsum(xs)
    rows.append((math.log10(n), math.log10(err)))
    print(f"{n:>10} {err:>12.3e} {bound:>17.3e}")
    assert err <= bound

mx = sum(x for x, _ in rows) / len(rows)
my = sum(y for _, y in rows) / len(rows)
slope = sum((x - mx) * (y - my) for x, y in rows) / sum((x - mx) ** 2 for x, _ in rows)
print(f"least-squares slope of log10(error) against log10(n): {slope:.3f}")
assert 1.8 < slope < 2.2
```

```text
         n    abs error  bound n*u*sum|x|
      1000    1.407e-12         1.110e-11
     10000    1.588e-10         1.110e-09
    100000    1.885e-08         1.110e-07
   1000000    1.333e-06         1.110e-05
  10000000    1.610e-04         1.110e-03
least-squares slope of log10(error) against log10(n): 2.004
```

The slope is about 2: the error grows like $n^2$, not $n$. That matches the naive bound $n\,u\sum|x_i|$ from §3.6: here $\sum|x_i| = 0.1n$ also grows with $n$, so the bound is proportional to $n^2$, and the measured error stays within a factor of about seven of it, because the same $0.1$ is rounded the same way each time. Each addition errs by up to half the spacing at the size of the running total, and the total is proportional to the number of terms added so far.

## C. Challenge

### C1. Walking the grid

```python
# file: c1.py
import math
import random
import struct


def to_bits(x):
    return int.from_bytes(struct.pack(">d", x), "big")


def from_bits(n):
    return struct.unpack(">d", n.to_bytes(8, "big"))[0]


def next_up(x):
    """Next double above x, by integer arithmetic on the bit pattern."""
    if x != x or x == math.inf:
        return x
    if x == 0.0:
        return from_bits(1)           # +0.0 and -0.0 both go to the smallest subnormal
    n = to_bits(x)
    return from_bits(n + 1) if n >> 63 == 0 else from_bits(n - 1)


rng = random.Random(0)
checked = 0
for _ in range(100_000):
    n = rng.getrandbits(64)
    x = from_bits(n)
    if math.isnan(x) or math.isinf(x):
        continue
    assert to_bits(next_up(x)) == to_bits(math.nextafter(x, math.inf)), x
    checked += 1
for x in (0.0, -0.0, 5e-324, -5e-324, 1.0, -1.0, 1.7976931348623157e308, math.inf, -math.inf):
    assert to_bits(next_up(x)) == to_bits(math.nextafter(x, math.inf)), x
assert math.isnan(next_up(math.nan))
print("checked", checked, "random patterns and 9 special values")
```

```text
checked 99952 random patterns and 9 special values
```

Why it works: the 64-bit pattern of a positive double reads sign (0), then $E$, then $F$, with $E$ in the more significant bits. For positive doubles, a larger pattern means a larger number, because a larger $E$ outweighs any $F$ and, within one $E$, a larger $F$ is a larger significand. Adding 1 to the integer increments $F$; when $F$ is all ones the carry propagates into $E$ and $F$ becomes 0, which is exactly the move from $1.111\ldots_2 \times 2^e$ to $1.0 \times 2^{e+1}$. The same carry takes the largest double to `inf`, and from the subnormals ($E = 0$) to the smallest normal. Negative numbers have the sign bit set and the same magnitude ordering, so "up" means a smaller magnitude: subtract 1; `-5e-324` goes to `-0.0`. Zeros are handled separately because both have magnitude pattern 0.

### C2. Design: `x or default`

```python
# file: c2.py
def configure_or(timeout=None, retries=None, name=None):
    """The `or` idiom: replaces every falsy value, including legitimate ones."""
    return timeout or 30, retries or 3, name or "anon"


def configure(timeout=30, retries=3, name="anon"):
    """Defaults live in the signature; they are immutable, so this is safe (R4)."""
    return timeout, retries, name


def configure_none_ok(timeout=None, retries=None, name=None):
    """Same, for callers that may pass None explicitly (say, from parsed JSON)."""
    if timeout is None:
        timeout = 30
    if retries is None:
        retries = 3
    if name is None:
        name = "anon"
    return timeout, retries, name


assert configure_or(timeout=0, retries=0, name="") == (30, 3, "anon")     # explicit values lost
assert configure(timeout=0, retries=0, name="") == (0, 0, "")
assert configure() == configure_none_ok() == (30, 3, "anon")
assert configure_none_ok(timeout=0, retries=None, name="") == (0, 3, "")
print(configure_or(timeout=0, retries=0, name=""), configure(timeout=0, retries=0, name=""))
```

```text
(30, 3, 'anon') (0, 0, '')
```

(a) `timeout=0` is a legitimate "do not wait"; `retries=0` is "never retry"; `name=""` may be a deliberate empty label (an empty list passed for a collection is also valid input). The idiom turns each into the default, silently, and the call site gives no hint. (b) For immutable defaults put them in the signature: `configure(timeout=30, ...)`. If callers may pass `None` explicitly, test `is None`. For mutable defaults use the `None` sentinel (Example 4.6).

(c) *Rule for review:* accept `x or default` only when every falsy value of `x` is meaningless as an input, so replacing it loses nothing, and say so in the docstring. A reasonable case is `label = user_label or "untitled"`, where an empty label and a missing label should both display as "untitled". Reject it whenever `0`, `0.0`, `""`, `[]` or `False` could be a real request: counts, timeouts, offsets, thresholds, flags. The cost of the longer `is None` test is two lines; the cost of the bug is a wrong value that no exception reports.

### C3. Does summation order matter?

```python
# file: c3.py
import math
from itertools import accumulate

N = 10**6


def naive_sum(xs):
    total = 0.0
    for x in xs:
        total += x
    return total


def kahan_sum(xs):
    total = c = 0.0
    for x in xs:
        y = x - c
        t = total + y
        c = (t - total) - y
        total = t
    return total


def pairwise_sum(xs, lo=0, hi=None):
    """Split in halves; add blocks of at most 8 terms directly."""
    if hi is None:
        hi = len(xs)
    if hi - lo <= 8:
        s = 0.0
        for i in range(lo, hi):
            s += xs[i]
        return s
    mid = (lo + hi) // 2
    return pairwise_sum(xs, lo, mid) + pairwise_sum(xs, mid, hi)


series = {
    "1/k^2": [1.0 / (k * k) for k in range(1, N + 1)],
    "alternating 1/k": [(-1) ** (k + 1) / k for k in range(1, N + 1)],
}
results = {}
for label, terms in series.items():
    ref = math.fsum(terms)
    methods = {
        "forward": naive_sum(terms),
        "backward": naive_sum(reversed(terms)),
        "sorted(terms)": naive_sum(sorted(terms)),
        "sorted by |x|": naive_sum(sorted(terms, key=abs)),
        "kahan": kahan_sum(terms),
        "pairwise": pairwise_sum(terms),
    }
    results[label] = {name: abs(v - ref) for name, v in methods.items()}
    print(f"{label}: fsum = {ref!r}")
    for name, err in results[label].items():
        print(f"  {name:>14}: error {err:.2e}")

low = min(accumulate(sorted(series["alternating 1/k"])))     # lowest running total, sorted by value
print(f"alternating, sorted(terms): lowest running total {low:.3f}, grid spacing there {math.ulp(low):.2e}")
assert -6.9 < low < -6.8
r = results["1/k^2"]
assert r["kahan"] < 1e-15 and r["forward"] > 1e-14
assert max(r["backward"], r["sorted(terms)"], r["pairwise"]) <= 4.5e-16
a = results["alternating 1/k"]
assert a["sorted(terms)"] > 1e-14 and a["sorted(terms)"] > a["forward"] > 1e-14
assert a["kahan"] < 1e-15 and a["sorted by |x|"] < 1e-15
```

```text
1/k^2: fsum = 1.6449330668487265
         forward: error 4.35e-14
        backward: error 2.22e-16
   sorted(terms): error 2.22e-16
   sorted by |x|: error 2.22e-16
           kahan: error 0.00e+00
        pairwise: error 2.22e-16
alternating 1/k: fsum = 0.6931466805601953
         forward: error 5.72e-14
        backward: error 0.00e+00
   sorted(terms): error 1.57e-13
   sorted by |x|: error 0.00e+00
           kahan: error 0.00e+00
        pairwise: error 0.00e+00
alternating, sorted(terms): lowest running total -6.850, grid spacing there 8.88e-16
```

The error of each addition is about half the grid spacing at the size of the *running total*. For $1/k^2$ the total climbs to $1.64$ quickly and then barely changes, so a forward loop adds $10^6$ small terms to a total whose grid spacing ($2.2\times10^{-16}$) is coarse compared with them: error $4\times10^{-14}$. Backward, ascending and pairwise orders add small numbers to small numbers; the total is large only at the very end, so the error is about one final rounding, $2.2\times10^{-16}$. Kahan carries each rounding error forward and is exact here (as `fsum`).

For the alternating series the ranking changes. `sorted(terms)` sorts by *value*, so all the negative terms come first: the running total travels to about $-6.85$ (half the harmonic sum $H_{500000}$) before the positive terms bring it back to $0.69$, and the error is set by the spacing at $6.85$ ($8.9\times10^{-16}$), giving $1.6\times10^{-13}$, worse than forward. Sorting by magnitude (`key=abs`), Kahan and pairwise are exact or nearly so. The lesson: the bound is governed by the size of the *partial sums*, not of the final answer, so keep partial sums small, or use `math.fsum`.

### Project, "Your turn" item 3: the largest exact power of ten

```python
# file: pow10.py
exact = [n for n in range(300) if float(10**n) == 10**n]   # int == float compares exactly
print(max(exact), exact == list(range(23)))
assert max(exact) == 22 and 5**22 < 2**53 < 5**23
```

```text
22 True
```

$10^n = 2^n \cdot 5^n$. The factor $2^n$ only shifts the exponent, so `float(10**n)` is exact precisely when the odd part $5^n$ fits in the 53 significant bits of a double: $5^n < 2^{53} \approx 9.007\times10^{15}$. Since $5^{22} \approx 2.38\times10^{15}$ fits and $5^{23} \approx 1.19\times10^{16}$ does not (and $5^n$ is odd, so it cannot be shortened by trailing zeros), $n = 22$ is the largest. The scan up to 299 confirms that no larger $n$ is exact.
