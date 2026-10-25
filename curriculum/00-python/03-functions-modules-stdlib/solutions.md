# PY-3 Solutions

Every code block below was run, and the outputs shown are the real ones. Blocks with `assert`s passed them. Timing-based solutions assert only orderings and wide margins, so they pass on any ordinary machine but print numbers that differ from yours. Shell sessions ran in the scratch directory `/tmp/py3` (the module folder for `numtools`), so tracebacks show that path.

## A. Warm-up

**A1.** All three lambdas in the first list share one variable `i`, and the loop leaves it at 2, so each adds 2. In the second list `i=i` is a default argument, evaluated once per iteration when the `lambda` expression runs.

```python
fs = []
for i in range(3):
    fs.append(lambda x: x + i)
print([f(10) for f in fs], "i =", i)
gs = [lambda x, i=i: x + i for i in range(3)]
print([g(10) for g in gs])
assert [f(10) for f in fs] == [12, 12, 12] and [g(10) for g in gs] == [10, 11, 12]
```
```text
[12, 12, 12] i = 2
[10, 11, 12]
```

**A2.** `wrapper` replaced `greet`, so the name and the docstring are `wrapper`'s. `functools.wraps(func)` copies them back (and sets `__wrapped__`).

```python
import functools

def shout(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs).upper()
    return wrapper

@shout
def greet(name):
    """Say hello."""
    return f"hello, {name}"

print(greet("ada"), greet.__name__, greet.__doc__)
assert (greet.__name__, greet.__doc__, greet.__wrapped__("x")) == ("greet", "Say hello.", "hello, x")
```
```text
HELLO, ADA greet Say hello.
```

Without `wraps` the last two values would be `wrapper` and `None`.

**A3.** `a` is positional-only, `b` is normal, `c` is keyword-only, and extra keywords go to `kw`. In the last call `a=1` cannot bind `a`, so with `**kw` present it lands in `kw` and `a` is reported missing.

```pycon
>>> def f(a, /, b, *, c=3, **kw):
...     return a, b, c, kw
...
>>> f(1, 2)
(1, 2, 3, {})
>>> f(1, b=2, c=5, z=0)
(1, 2, 5, {'z': 0})
>>> f(*(1, 2), **{"c": 4})
(1, 2, 4, {})
>>> f(1, 2, 3)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: f() takes 2 positional arguments but 3 were given
>>> f(a=1, b=2)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: f() missing 1 required positional argument: 'a'
```

**A4.** `tri(5)` calls `tri(5), tri(4), ..., tri(0)`: 6 calls (a linear chain, $n + 1$ calls). For `g`, let $G(n)$ be the number of calls: $G(0) = 1$ and $G(n) = 2G(n-1) + 1$. Then $G(n) + 1 = 2(G(n-1) + 1)$, so $G(n) + 1 = 2^n(G(0) + 1) = 2^{n+1}$ and $G(n) = 2^{n+1} - 1$: $G(3) = 15$.

```python
count = 0

def tri(n):
    global count
    count += 1
    return 0 if n == 0 else n + tri(n - 1)

def g(n):
    global count
    count += 1
    return 1 if n == 0 else g(n - 1) + g(n - 1)

for func, n in ((tri, 5), (g, 3), (g, 10)):
    count = 0
    func(n)
    print(func.__name__, n, count)
count = 0; tri(5); assert count == 6
count = 0; g(3); assert count == 15
count = 0; g(10); assert count == 2**11 - 1
```
```text
tri 5 6
g 3 15
g 10 2047
```

**A5.** `import util` runs `util.py` once (printing with `__name__ == "util"`); the other two imports find it in `sys.modules`. Run as a script, the same file has `__name__ == "__main__"`, so the guard fires.

```python
# file: util.py
print("util loaded", __name__)

def f():
    return 1

if __name__ == "__main__":
    print("util main")
```

```python
# file: main.py
import util
from util import f
import util as u

print(f(), u is util)
```

```console
$ python3 main.py
util loaded util
1 True
$ python3 util.py
util loaded __main__
util main
```

**A6.** The returned object is a function. `n` and `f` live in cells of the closure, which stay alive as long as the returned function does.

```python
def apply_n(f, n):
    def run(x):
        for _ in range(n):
            x = f(x)
        return x
    return run

double5 = apply_n(lambda x: 2 * x, 5)
print(double5(1), type(double5).__name__, double5.__code__.co_freevars, [c.cell_contents for c in double5.__closure__][1])
assert double5(1) == 32 and apply_n(abs, 0)(-3) == -3 and apply_n(abs, 3)(-3) == 3
```
```text
32 function ('f', 'n') 5
```

## B. Core

**B1.** `counted` stores its state on the wrapper (a function has a `__dict__`). The decorated name must be the one the recursion calls, which it is when `@counted` sits on the `def`.

```python
import functools

def counted(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        wrapper.calls += 1
        return func(*args, **kwargs)
    wrapper.calls = 0
    wrapper.reset = lambda: setattr(wrapper, "calls", 0)
    return wrapper

@counted
def fib(n):
    return n if n < 2 else fib(n - 1) + fib(n - 2)

def fib_iter(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

for n in range(21):
    fib.reset()
    fib(n)
    assert fib.calls == 2 * fib_iter(n + 1) - 1                  # C(n) = 2 F(n+1) - 1
print("C(20) =", fib.calls, "| name kept:", fib.__name__)

def trace(func):
    depth = 0
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        nonlocal depth
        shown = ", ".join([repr(a) for a in args] + [f"{k}={v!r}" for k, v in kwargs.items()])
        print(f"{'  ' * depth}{func.__name__}({shown})")
        depth += 1
        try:
            result = func(*args, **kwargs)
        finally:
            depth -= 1
        print(f"{'  ' * depth}-> {result!r}")
        return result
    return wrapper

@trace
def tri(n):
    return 0 if n == 0 else n + tri(n - 1)

@trace
def gcd(a, b):
    return a if b == 0 else gcd(b, a % b)

assert tri(3) == 6
assert gcd(48, 18) == 6
```
```text
C(20) = 21891 | name kept: fib
tri(3)
  tri(2)
    tri(1)
      tri(0)
      -> 0
    -> 1
  -> 3
-> 6
gcd(48, 18)
  gcd(18, 12)
    gcd(12, 6)
      gcd(6, 0)
      -> 6
    -> 6
  -> 6
-> 6
```

`C(20) = 21891` agrees with $2F_{21} - 1 = 21891$ (`fib_iter(21)` is 10946). The trace of `gcd` shows Euclid's algorithm as a chain of calls whose arguments shrink: the last non-zero remainder is the answer.

**B2.** `OrderedDict` keeps insertion order and offers `move_to_end(key)` and `popitem(last=False)`, both $O(1)$: that is "mark as most recent" and "evict the least recent" without scanning. A plain dict can move a key to the end in $O(1)$ (`d[k] = d.pop(k)`), but its oldest entry is `next(iter(d))`, which walks over the slots left by earlier deletions and slows down as the cache churns; `OrderedDict` keeps a linked list, so both operations stay $O(1)$.

```python
import functools
import random
from collections import OrderedDict

def lru(maxsize):
    def decorator(func):
        cache = OrderedDict()

        @functools.wraps(func)
        def wrapper(*args):
            if args in cache:
                cache.move_to_end(args)
                wrapper.hits += 1
                return cache[args]
            wrapper.misses += 1
            value = cache[args] = func(*args)
            if len(cache) > maxsize:
                cache.popitem(last=False)
            return value

        def cache_clear():
            cache.clear()
            wrapper.hits = wrapper.misses = 0

        wrapper.hits = wrapper.misses = 0
        wrapper.cache_clear = cache_clear
        return wrapper
    return decorator

rng = random.Random(0)
requests = [rng.randint(0, 99) for _ in range(10_000)]
for m in (5, 20, 50):
    mine = lru(m)(lambda k: k * k)
    ref = functools.lru_cache(maxsize=m)(lambda k: k * k)
    for k in requests:
        assert mine(k) == ref(k) == k * k
    info = ref.cache_info()
    print(f"maxsize {m:>2}: mine ({mine.hits}, {mine.misses})   functools ({info.hits}, {info.misses})")
    assert (mine.hits, mine.misses) == (info.hits, info.misses)
```
```text
maxsize  5: mine (500, 9500)   functools (500, 9500)
maxsize 20: mine (1976, 8024)   functools (1976, 8024)
maxsize 50: mine (4992, 5008)   functools (4992, 5008)
```

**B3.** The factory has three levels; the check goes in the innermost. Written without parentheses, `@bounded` calls `bounded(func)` with one argument. If `hi` is a required parameter, Python itself reports a missing `hi` before our code runs, which does not say what to do, and a `callable(lo)` check in the body is never reached. Giving `hi` a `None` default lets the body run, so we can detect the single callable argument and raise a message that names the fix.

```python
import functools
import math

from numtools import bisect

def bounded(lo, hi=None):
    if hi is None or callable(lo):
        raise TypeError("bounded needs two bounds: write @bounded(lo, hi), not @bounded")
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            value = func(*args, **kwargs)
            if not lo <= value <= hi:
                raise ValueError(f"{func.__name__} returned {value!r}, outside [{lo}, {hi}]")
            return value
        return wrapper
    return decorator

@bounded(0.0, 1.0)
def fixed_point(tol=1e-12):
    return bisect(lambda x: math.cos(x) - x, 0.0, 1.0, tol=tol)[0]

@bounded(0.0, 0.5)
def too_big():
    return fixed_point()

print(round(fixed_point(tol=1e-6), 6))
for thunk in (too_big, lambda: bounded(lambda: 1)):
    try:
        thunk()
    except (ValueError, TypeError) as err:
        print(type(err).__name__ + ":", err)
assert fixed_point.__name__ == "fixed_point" and 0 <= fixed_point() <= 1
```
```text
0.739085
ValueError: too_big returned 0.7390851332156672, outside [0.0, 0.5]
TypeError: bounded needs two bounds: write @bounded(lo, hi), not @bounded
```

**B4.** (a) With the default limit of 1000, all three functions fail at the same depth on this build, 997 (the script itself and the measuring loop use the remaining frames; do not memorize 997, since versions differ by a few). The limit counts *Python frames*: the number of locals changes the memory per frame, not the count, and the C-implemented `lru_cache` wrapper adds no Python frame of its own. The cache is cleared before every trial; otherwise each call would find its predecessor cached and never recurse deeply (a cache does not make recursion shallower, it only avoids repeated work). (b) A plain recursion needs $n$ frames. A loop needs none. Splitting the range in half at each step needs only $\log_2 n$ frames, because the two halves are processed one after the other, not stacked.

```python
import functools
import sys

def max_depth(func, reset=lambda: None):
    """Largest n for which func(n) does not raise RecursionError (reset() empties any cache first)."""
    n = 1
    while True:
        reset()
        try:
            func(n)
        except RecursionError:
            return n - 1
        n += 1

def down(n):
    return 0 if n == 0 else down(n - 1)

def down3(n):
    a, b, c = n, n + 1, n + 2
    return 0 if n == 0 else down3(n - 1) + (a - a) + (b - b) + (c - c)

@functools.lru_cache(maxsize=None)
def down_cached(n):
    return 0 if n == 0 else down_cached(n - 1)

print("limit", sys.getrecursionlimit())
depths = {}
for func, reset in ((down, lambda: None), (down3, lambda: None), (down_cached, down_cached.cache_clear)):
    depths[func.__name__] = max_depth(func, reset)
print(depths)

def sum_to(n):
    return 0 if n == 0 else n + sum_to(n - 1)

def sum_loop(n):
    total = 0
    for k in range(1, n + 1):
        total += k
    return total

def sum_halves(lo, hi, depth=1):
    """Sum of lo..hi inclusive, and the deepest recursion level reached."""
    if lo == hi:
        return lo, depth
    mid = (lo + hi) // 2
    left, d1 = sum_halves(lo, mid, depth + 1)
    right, d2 = sum_halves(mid + 1, hi, depth + 1)
    return left + right, max(d1, d2)

try:
    sum_to(10_000)
except RecursionError as err:
    print("sum_to(10000):", type(err).__name__)
total, d = sum_halves(1, 10**6)
print(total, "deepest level", d)
assert total == sum_loop(10**6) == 10**6 * (10**6 + 1) // 2 and d <= 21
assert depths["down"] == depths["down3"] and all(v < 1000 for v in depths.values())
```
```text
limit 1000
{'down': 997, 'down3': 997, 'down_cached': 997}
sum_to(10000): RecursionError
500000500000 deepest level 21
```

**B5.** The package lives in `geom/`. `stats.py` imports `dot` with a relative import, so it only works when imported as part of the package.

```python
# file: geom/vectors.py
import math

def dot(u, v):
    return sum(a * b for a, b in zip(u, v, strict=True))

def norm(u):
    return math.sqrt(dot(u, u))
```

```python
# file: geom/stats.py
from .vectors import dot

def mean(xs):
    return sum(xs) / len(xs)

def variance(xs):
    centered = [x - mean(xs) for x in xs]
    return dot(centered, centered) / (len(xs) - 1)
```

```python
# file: geom/__init__.py
"""geom: tiny vector and statistics helpers."""
from .stats import mean, variance
from .vectors import dot, norm

__all__ = ["dot", "norm", "mean", "variance"]
```

```python
# file: geom/__main__.py
from . import mean, norm, variance

print("norm of (3, 4):", norm([3, 4]))
print("mean and variance of 1..5:", mean([1, 2, 3, 4, 5]), variance([1, 2, 3, 4, 5]))
```

```console
$ python3 -c "import geom; print(geom.__all__)"
['dot', 'norm', 'mean', 'variance']
$ python3 -c "from geom import norm; print(norm([3, 4]))"
5.0
$ python3 -m geom
norm of (3, 4): 5.0
mean and variance of 1..5: 3.0 2.5
$ python3 geom/stats.py
Traceback (most recent call last):
  File "/tmp/py3/geom/stats.py", line 2, in <module>
    from .vectors import dot
ImportError: attempted relative import with no known parent package
$ cd geom && python3 -c "import geom"
Traceback (most recent call last):
  File "<string>", line 1, in <module>
    import geom
ModuleNotFoundError: No module named 'geom'
```

Reading the results with `sys.path[0]`: in the first three commands Python starts in the parent directory, so `sys.path[0]` is the current directory (for `-c` and `-m`), which contains the folder `geom`: found. `python3 geom/stats.py` makes `geom/` itself `sys.path[0]` and runs `stats.py` as a top-level script, so `from .vectors import dot` has no package to be relative to. And inside `geom/`, `sys.path[0]` is `geom/`, which contains `vectors.py` but not a folder called `geom`, so `import geom` fails. The working directory that makes everything work is the parent of the package.

**B6.** The parser has two subcommands that share their options through a helper; `main(argv)` returns the status, so tests need no subprocess.

```python
# file: summarize.py
import argparse
import csv
import statistics
import sys
from pathlib import Path


def read_column(path, column, delimiter, skip_header):
    values, skipped = [], 0
    with Path(path).open(newline="", encoding="utf-8") as f:
        rows = csv.reader(f, delimiter=delimiter)
        if skip_header:
            next(rows, None)
        for row in rows:
            try:
                values.append(float(row[column]))
            except (ValueError, IndexError):
                skipped += 1
    return values, skipped


def show_mean(values):
    print(f"mean = {statistics.mean(values):.6g}")


def show_minmax(values):
    print(f"min = {min(values):.6g}, max = {max(values):.6g}")


def build_parser():
    parser = argparse.ArgumentParser(prog="summarize", description="Summarize one column of a delimited file.")
    sub = parser.add_subparsers(dest="command", required=True)
    for name, func, help_text in (("mean", show_mean, "arithmetic mean"), ("minmax", show_minmax, "smallest and largest")):
        p = sub.add_parser(name, help=help_text)
        p.add_argument("file")
        p.add_argument("--column", type=int, default=0, help="zero-based column (default: %(default)s)")
        p.add_argument("--delimiter", default=",", help="field separator (default: %(default)r)")
        p.add_argument("--skip-header", action="store_true", help="ignore the first row")
        p.set_defaults(func=func)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        values, skipped = read_column(args.file, args.column, args.delimiter, args.skip_header)
    except OSError as err:
        print(f"error: {err.strerror}: {args.file!r}", file=sys.stderr)
        return 2
    if skipped:
        print(f"warning: skipped {skipped} non-numeric cell(s)", file=sys.stderr)
    if not values:
        print("error: no numeric values found", file=sys.stderr)
        return 2
    args.func(values)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

```python
import contextlib
import io
from pathlib import Path

import summarize

Path("data.csv").write_text("name,value\na,1.5\nb,2.5\nc,n/a\nd,4.0\n")

def run(*argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        status = summarize.main(list(argv))
    return status, out.getvalue().strip(), err.getvalue().strip()

print(run("mean", "data.csv", "--column", "1", "--skip-header"))
print(run("minmax", "data.csv", "--column", "1", "--skip-header"))
print(run("mean", "nofile.csv"))
assert run("mean", "data.csv", "--column", "1", "--skip-header") == (0, "mean = 2.66667", "warning: skipped 1 non-numeric cell(s)")
assert run("minmax", "data.csv", "--column", "1", "--skip-header")[:2] == (0, "min = 1.5, max = 4")
assert run("mean", "nofile.csv")[0] == 2 and run("mean", "data.csv", "--column", "7")[0] == 2
```
```text
(0, 'mean = 2.66667', 'warning: skipped 1 non-numeric cell(s)')
(0, 'min = 1.5, max = 4', 'warning: skipped 1 non-numeric cell(s)')
(2, '', "error: No such file or directory: 'nofile.csv'")
```

Skipping the header, the numeric cells of column 1 are $1.5, 2.5, 4.0$, mean $8/3 = 2.66667$; the `n/a` cell is skipped and counted. Without `--skip-header`, the word `value` would be skipped as non-numeric, with the same effect and a warning for two cells. A bad column index is an `IndexError` per row, handled the same way until no values remain.

**B7.** Add the following to Part 4 of `project.py`, after the `t_naive`, `t_memo` lists. Here it runs as a standalone script that imports the package. The cache is cleared before every timed run so that every run is cold.

```python
import functools
import time

from numtools import memoize

def naive(n):
    return n if n < 2 else naive(n - 1) + naive(n - 2)

@memoize
def memo(n):
    return n if n < 2 else memo(n - 1) + memo(n - 2)

@functools.lru_cache(maxsize=None)
def cached(n):
    return n if n < 2 else cached(n - 1) + cached(n - 2)

def fib_iter(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

def best(func, n, repeat, before=lambda: None):
    times = []
    for _ in range(repeat):
        before()
        start = time.perf_counter()
        func(n)
        times.append(time.perf_counter() - start)
    return min(times)

n = 27
t = {"naive": best(naive, n, 3),
     "memoize": best(memo, n, 50, memo.cache.clear),
     "lru_cache": best(cached, n, 50, cached.cache_clear),
     "iterative": best(fib_iter, n, 50)}
for name, value in t.items():
    print(f"{name:>10}: {value * 1e6:>12.1f} us")
print(f"memoize / lru_cache = {t['memoize'] / t['lru_cache']:.1f}   lru_cache / iterative = {t['lru_cache'] / t['iterative']:.1f}")
assert t["iterative"] < t["lru_cache"] < t["memoize"] < t["naive"]
assert memo(n) == cached(n) == fib_iter(n) == 196418
```
```text
     naive:      32788.4 us
   memoize:         22.3 us
 lru_cache:          6.1 us
 iterative:          0.9 us
memoize / lru_cache = 3.6   lru_cache / iterative = 6.5
```

`lru_cache`'s wrapper is written in C: it builds the key and does the lookup without creating a Python frame for the wrapper, while `memoize` pays for a Python call, an `*args` tuple, a `**kwargs` dict, a `sorted` call on its items, a tuple key, and two dictionary lookups per call. The loop does none of that: it makes no function calls, 27 additions on small integers. Memoization removes the exponential; choosing the loop removes the constant.

## C. Challenge

**C1.** With `func=None` as the first parameter and everything else keyword-only, the decorator sees one of two calls: `counted(f)` (bare use: `func` is the function) or `counted(label="x")` (with options: `func` is `None`, and we return the real decorator). A positional argument that is not callable, such as `@counted("fib")`, is neither: without a check, `counted("fib")` would try to wrap the string and fail at `f.__name__` with `AttributeError: 'str' object has no attribute '__name__'` (`functools.wraps` itself tolerates missing attributes), a message that does not name the real mistake. We test `callable(func)` first. A function passed positionally, `@counted(fib)`, cannot be told from bare use: it wraps `fib` and returns the wrapper, which is then applied to the decorated function and silently misbehaves; the signature cannot prevent it, which is why the options are keyword-only.

```python
import functools

def counted(func=None, *, label=None):
    if func is not None and not callable(func):
        raise TypeError(f"counted() got {func!r}; pass options by keyword: @counted(label=...)")

    def decorate(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            wrapper.calls += 1
            return f(*args, **kwargs)
        wrapper.calls, wrapper.label = 0, label or f.__name__
        return wrapper

    return decorate if func is None else decorate(func)

@counted
def double(x):
    return 2 * x

@counted(label="triple!")
def triple(x):
    return 3 * x

class Scale:
    def __init__(self, k):
        self.k = k

    @counted
    def apply(self, x):
        return self.k * x

a, b = Scale(2), Scale(5)
assert (double(2), double(3), triple(1)) == (4, 6, 3) and (a.apply(1), b.apply(1), a.apply(2)) == (2, 5, 4)
print(double.calls, double.label, triple.calls, triple.label, Scale.apply.calls)
try:
    counted("fib")
except TypeError as err:
    print("TypeError:", err)
assert (double.calls, triple.calls, Scale.apply.calls, Scale.apply.__name__) == (2, 1, 3, "apply")
```
```text
2 double 1 triple! 3
TypeError: counted() got 'fib'; pass options by keyword: @counted(label=...)
```

The method decorated this way counts calls across *all* instances (the wrapper lives on the class), and `self` needs no special handling because it is the first element of `*args`.

**C2.** Once the cache is full it holds the $m$ most recently used distinct keys. If requests are independent and uniform over $N$ keys, the next request is one of those $m$ with probability exactly $m/N$, whatever the order of past requests. For a Zipf-like workload, popular keys are requested so often that they are almost always among the most recent: the hit rate is far above $m/N$.

```python
import functools
import random

N, REQUESTS = 1000, 100_000
workloads = {
    "uniform": random.Random(0).choices(range(N), k=REQUESTS),
    "zipf": random.Random(0).choices(range(N), weights=[1 / k for k in range(1, N + 1)], k=REQUESTS),
}
print(f"{'m':>4} {'m/N':>6} {'uniform':>8} {'zipf':>8}")
for m in (50, 100, 250, 500):
    rates = {}
    for name, keys in workloads.items():
        cached = functools.lru_cache(maxsize=m)(lambda k: k)
        for k in keys:
            cached(k)
        info = cached.cache_info()
        rates[name] = info.hits / (info.hits + info.misses)
    print(f"{m:>4} {m / N:>6.3f} {rates['uniform']:>8.3f} {rates['zipf']:>8.3f}")
    assert abs(rates["uniform"] - m / N) < 0.02 and rates["zipf"] > rates["uniform"] + 0.2
```
```text
   m    m/N  uniform     zipf
  50  0.050    0.050    0.463
 100  0.100    0.100    0.576
 250  0.250    0.248    0.735
 500  0.500    0.497    0.864
```

The uniform column matches $m/N$ to within sampling error (the first $m$ requests are forced misses, a bias of about $m/\text{requests}$, at most 0.5% here). The Zipf column is far larger: popular keys keep returning before they are evicted. The gap is widest for mid-sized caches and closes only when the cache holds nearly everything. Real workloads are skewed, which is why a small cache is worth having.

**C3.** The recursive generator fixes the first element and delegates the rest with `yield from`; the iterative version keeps one iterator per level of the recursion in an explicit stack, which is exactly what the Python call stack was holding.

```python
import itertools
import time

def perms(items):
    items = list(items)
    if len(items) <= 1:
        yield tuple(items)
        return
    for i, x in enumerate(items):
        yield from ((x,) + rest for rest in perms(items[:i] + items[i + 1:]))

def perms_iter(items):
    items, n = tuple(items), len(items)
    if n == 0:
        yield ()
        return
    used, prefix, stack = [False] * n, [], [iter(range(n))]
    while stack:
        for i in stack[-1]:
            if used[i]:
                continue
            used[i] = True
            prefix.append(i)
            if len(prefix) == n:
                yield tuple(items[j] for j in prefix)
                prefix.pop()
                used[i] = False
                continue
            stack.append(iter(range(n)))
            break
        else:
            stack.pop()
            if prefix:
                used[prefix.pop()] = False

expected = list(itertools.permutations(range(6)))
assert list(perms(range(6))) == expected == list(perms_iter(range(6))) and len(expected) == 720
assert list(perms("")) == list(perms_iter("")) == [()] and list(perms_iter("a")) == [("a",)]

def best_time(func, repeat=3):
    times = []
    for _ in range(repeat):
        start = time.perf_counter()
        for _ in func(range(8)):
            pass
        times.append(time.perf_counter() - start)
    return min(times)

t = {"recursive": best_time(perms), "explicit stack": best_time(perms_iter), "itertools": best_time(itertools.permutations)}
for name, value in t.items():
    print(f"{name:>14}: {value * 1e3:8.2f} ms for 40,320 permutations")
assert t["itertools"] * 5 < min(t["recursive"], t["explicit stack"])
```
```text
     recursive:    71.22 ms for 40,320 permutations
explicit stack:    63.43 ms for 40,320 permutations
     itertools:     1.88 ms for 40,320 permutations
```

`itertools.permutations` is written in C and wins by a large factor (more than 30 here); use it. The two Python versions run at about the same speed. The recursive one is much shorter and mirrors the definition (choose the first element, permute the rest); its depth is only $n = 8$ generators, so there is no depth problem to solve and the explicit stack buys nothing. Convert recursion to an explicit stack when the depth is proportional to the input size (the 50,000-deep nested list in `project.py`'s "Your turn"), when you need to pause and resume or limit memory, or when the recursion limit is the actual obstacle; do not convert for speed before measuring.
