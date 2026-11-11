# %% [markdown]
# # PY-3 Project: `numtools`, functions as values, and recursion with measurements
#
# Exercises the package in `numtools/` (PY-1 root finders, PY-2 sparse polynomial, the decorators
# `timer`, `memoize`, `retry`, `compose`/`pipeline`, and a CLI) and checks every number quoted in
# the notes: closures and late binding (3.4, Example 4.2), the dispatch table (4.1), `timer` on the
# root finders (4.3), the call count 2 F(n+1) - 1 and its time ratio (3.7), `lru_cache` (4.4), the
# cofactor determinant against NumPy (4.7), binary search with an executable invariant, merge sort,
# regular-expression parsing (4.5), the import system and the command line (3.8, 4.6).
# Two figures go to figures/. NumPy and Matplotlib only.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Build first, then compare: write your own
# timer, memoize and retry before reading numtools/decorators.py.

# %%
import ast
import cmath
import contextlib
import functools
import importlib
import inspect
import io
import math
import operator
import os
import pprint
import random
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import textwrap
import time
from bisect import bisect_left
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
# `python3 project.py` already puts this folder first on sys.path. Run cell by cell from another
# working directory and it is not there, so we insert it explicitly (notes 3.8).
sys.path.insert(0, str(HERE))

import numtools  # noqa: E402
from numtools import bisect, compose, memoize, newton, pipeline, retry, secant, timer  # noqa: E402
from numtools.cli import main as cli_main  # noqa: E402

PHI = (1 + 5**0.5) / 2

def raises(exc, func, *args, **kwargs):
    """Return the exception message if func(*args, **kwargs) raises exc, else None."""
    try:
        func(*args, **kwargs)
    except exc as err:
        return str(err)
    return None

def run_py(*args, cwd=HERE, check=True):
    """Run a fresh interpreter with normal sys.path[0] rules and no PYTHONPATH."""
    env = dict(os.environ, COLUMNS="80", PYTHONDONTWRITEBYTECODE="1")
    env.pop("PYTHONPATH", None)
    proc = subprocess.run([sys.executable, *args], cwd=cwd, env=env, capture_output=True, text=True, timeout=120)
    assert not check or proc.returncode == 0, proc.stderr
    return proc

def best(func, arg, repeat):
    """Smallest wall time of func(arg) over `repeat` runs (seconds)."""
    times = []
    for _ in range(repeat):
        start = time.perf_counter()
        func(arg)
        times.append(time.perf_counter() - start)
    return min(times)

def save_line_plot(name, title, xs, series, ylabel="seconds"):
    """Save a semilog line plot of {label: (ys, style)} to figures/<name>."""
    fig, ax = plt.subplots(figsize=(6, 4))
    for label, (ys, style) in series.items():
        ax.semilogy(xs, ys, style, label=label, lw=1 if style == "k--" else 1.5)
    ax.set(xlabel="n", ylabel=ylabel, title=title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / name, dpi=120)
    plt.close(fig)

# %% [markdown]
# ## Part 1: The import system, packages and the command line (notes 3.8, 3.9, Example 4.6)

# %%
out = run_py("-c", "import numtools, sys; print(numtools.__all__); print('numtools' in sys.modules)").stdout.splitlines()
assert ast.literal_eval(out[0]) == numtools.__all__ and out[1] == "True" and len(numtools.__all__) == 12
assert numtools.bisect.__module__ == "numtools.roots" and run_py("-m", "numtools", "--version").stdout == "numtools 0.1\n"
help_text = run_py("-m", "numtools", "--help").stdout                      # Example 4.6
assert help_text.startswith("usage: python -m numtools [-h] [--version] {bisect,newton} ...\n")
assert "bisection on a bracket [a, b]" in help_text and "derivative by central differences" in help_text

# (arguments, exit status, expected first lines of stdout, expected text in stderr)
CLI_CASES = [
    (("bisect", "cos(x)-x", "0", "1", "-v"), 0, ["root = 0.739085133216", "iterations = 40"], ""),
    (("newton", "x*x-2", "1.5", "--tol", "1e-14", "-v"), 0, ["root = 1.41421356237", "iterations = 5"], ""),
    (("bisect", "x*x+1", "-1", "1"), 2, [], "error: f(a) and f(b) must have opposite signs"),
    (("bisect", "cos(x)-x", "0"), 2, [], "the following arguments are required: b"),
]
for args, status, lines, err in CLI_CASES:
    proc = run_py("-m", "numtools", *args, check=False)
    assert proc.returncode == status and proc.stdout.splitlines()[:len(lines)] == lines and err in proc.stderr
assert 40 == math.ceil(math.log2(1 / 1e-12))
proc = run_py("numtools/cli.py", check=False)                              # relative import with no package
assert proc.returncode != 0 and "attempted relative import with no known parent package" in proc.stderr

buffer = io.StringIO()
with contextlib.redirect_stdout(buffer):                                    # main(argv) is testable in-process
    assert cli_main(["bisect", "x*x-2", "1", "2", "--tol", "1e-6"]) == 0
assert buffer.getvalue().strip() == "root = 1.41421413422"
with contextlib.redirect_stderr(io.StringIO()):
    assert raises(SystemExit, cli_main, ["bisect", "x", "one", "2"]) == "2"  # argparse's own error: exit 2
    assert cli_main(["newton", "log(x)", "-1"]) == 2                          # math domain error caught in main

FILES = {"chatty.py": "import builtins\nbuiltins.CHATTY_RUNS = getattr(builtins, 'CHATTY_RUNS', 0) + 1\nVALUE = 42\n",
         "tool.py": "def double(x):\n    return 2 * x\n\nif __name__ == '__main__':\n    print('running as a script:', double(21))\n",
         "noisy.py": "import random\nrandom.seed(0)\n",
         "proj/tools/helpers.py": "def hello():\n    return 'hello from tools.helpers'\n",
         "proj/scripts/run.py": "from tools.helpers import hello\nprint(hello())\n"}
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    for name, text in FILES.items():
        (tmp / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp / name).write_text(text)
    import builtins

    sys.path.insert(0, str(tmp))
    try:
        import chatty
        import chatty as c
        from chatty import VALUE
        assert builtins.CHATTY_RUNS == 1 and c is chatty and VALUE == 42 and sys.modules["chatty"] is chatty
        importlib.reload(chatty)
        assert builtins.CHATTY_RUNS == 2                                     # only an explicit reload runs it again
    finally:
        sys.path.remove(str(tmp))
        sys.modules.pop("chatty", None)
        del builtins.CHATTY_RUNS
    proj = tmp / "proj"
    for args, cwd, expected in [(("tool.py",), tmp, "running as a script: 42\n"),
                                (("-c", "import tool; print(tool.double(4), tool.__name__)"), tmp, "8 tool\n"),
                                (("-m", "tool"), tmp, "running as a script: 42\n"),
                                (("-m", "scripts.run"), proj, "hello from tools.helpers\n"),
                                (("-c", "import sys; print(sys.path[0] == '')"), proj, "True\n"),
                                (("-c", "import random; random.seed(123); import noisy; "
                                  "print(random.random() == random.Random(0).random())"), tmp, "True\n")]:
        assert run_py(*args, cwd=cwd).stdout == expected
    proc = run_py("scripts/run.py", cwd=proj, check=False)                  # sys.path[0] is scripts/, not proj/
    assert proc.returncode != 0 and "ModuleNotFoundError: No module named 'tools'" in proc.stderr
    proc = run_py("-c", "from os import *; open('notes.txt')", cwd=tmp, check=False)
    assert "open() missing required argument 'flags'" in proc.stderr
    if sys.platform != "win32":                                              # section 3.10: a venv is a private prefix
        run_py("-m", "venv", "--without-pip", ".venv", cwd=tmp)
        inside = subprocess.run([str(tmp / ".venv" / "bin" / "python"), "-c", "import sys; print(sys.prefix != sys.base_prefix)"],
                                capture_output=True, text=True)
        assert inside.stdout == "True\n"
print("Part 1: package, import system and CLI checks passed.")

# %% [markdown]
# ## Part 2: Functions as values (notes 1, 3.1 to 3.4, Examples 4.1 and 4.2)

# %%
CASES = [("cos x - x", lambda x: math.cos(x) - x, lambda x: -math.sin(x) - 1, 0.0, 1.0),
         ("x^2 - 2", lambda x: x * x - 2, lambda x: 2 * x, 1.0, 2.0),
         ("x^3 - 2x - 5", lambda x: x**3 - 2 * x - 5, lambda x: 3 * x * x - 2, 2.0, 3.0)]
expected = {"bisect": [34, 34, 34], "secant": [6, 7, 7]}
for k, (cname, f, df, a, b) in enumerate(CASES):                                   # the table of section 1
    for sname, solve in {"bisect": bisect, "secant": secant}.items():
        root, iters = solve(f, a, b, tol=1e-10)
        assert abs(f(root)) < 1e-8 and iters == expected[sname][k]

def area(r, pi=3.14159):
    """Area of a circle of radius r."""
    return pi * r * r

area.units = "m^2"
assert (type(area).__name__, area.__defaults__, area.__code__.co_varnames, area.__dict__) == (
    "function", (3.14159,), ("r", "pi"), {"units": "m^2"})
assert area(2) == 12.56636 and callable(area) and not callable(3)

def one():
    return 1

alias = one

def one():  # noqa: F811
    return 2

assert (alias(), one()) == (1, 2)
nums = [1, 2, 3, 4]
squares = map(lambda x: x * x, nums)
assert (list(squares), list(squares)) == ([1, 4, 9, 16], []) and functools.reduce(lambda acc, x: acc * 10 + x, nums) == 1234
words = ["pear", "Fig", "apple"]
assert [sorted(words), sorted(words, key=str.lower), sorted(words, key=len)] == [
    ["Fig", "apple", "pear"], ["apple", "Fig", "pear"], ["Fig", "pear", "apple"]]
assert functools.reduce(lambda acc, c: acc * 3 + c, [2, -6, 2, -1]) == 5          # Horner at x = 3 (PY-1)
to_kelvin = pipeline(str.strip, float, lambda c: c + 273.15)
assert math.isclose(compose(math.sin, math.sqrt)(2.0), 0.9877659459927356) and abs(to_kelvin("  21.5 ") - 294.65) < 1e-12
assert compose(len, str.split)("a b c d") == 4 and pipeline()(7) == 7 and compose(max, lambda *a: a)(3, 9, 4) == 9

def f_args(a, b=2, /, c=3, *args, d, e=5, **kw):
    return a, b, c, args, d, e, kw

assert f_args(1, d=4) == (1, 2, 3, (), 4, 5, {}) and f_args(1, 2, 3, 4, 5, d=6, z=7) == (1, 2, 3, (4, 5), 6, 5, {"z": 7})
assert f_args(*[1, 2], **{"d": 0}) == (1, 2, 3, (), 0, 5, {}) and str(inspect.signature(divmod)) == "(x, y, /)"
assert raises(TypeError, f_args, d=2, a=1) == "f_args() missing 1 required positional argument: 'a'"
assert raises(TypeError, f_args, 1) == "f_args() missing 1 required keyword-only argument: 'd'"

def make_adder(n):
    def add(x):
        return x + n
    return add

add5, add7 = make_adder(5), make_adder(7)
assert (add5(10), add7(10), add5.__code__.co_freevars, add5.__closure__[0].cell_contents) == (15, 17, ("n",), 5)
assert add5.__closure__[0] is not add7.__closure__[0] and make_adder.__closure__ is None
fs = [lambda: i for i in range(3)]                                                 # the late-binding trap
assert [f() for f in fs] == [2, 2, 2] and fs[0].__closure__[0] is fs[2].__closure__[0]

def power(x, k):
    return x**k

def make_power(k):
    return lambda x: x**k

variants = {"late": [lambda x: x**k for k in (1, 2, 3)], "default": [lambda x, k=k: x**k for k in (1, 2, 3)],
            "partial": [functools.partial(power, k=k) for k in (1, 2, 3)], "factory": [make_power(k) for k in (1, 2, 3)]}
assert [[f(2) for f in fs] for fs in variants.values()] == [[8, 8, 8], [2, 4, 8], [2, 4, 8], [2, 4, 8]]
assert variants["partial"][2].func is power and variants["partial"][2].keywords == {"k": 3}

OPS = {"+": operator.add, "-": operator.sub, "*": operator.mul, "/": operator.truediv, "**": operator.pow}

def op(name):
    def register(func):
        OPS[name] = func
        return func
    return register

@op("hypot")
def hypot(a, b):
    return math.hypot(a, b)

def calc(line):
    left, name, right = line.split()
    try:
        func = OPS[name]
    except KeyError:
        raise ValueError(f"unknown operator {name!r}; choose from {sorted(OPS)}") from None
    return func(float(left), float(right))

assert [calc(s) for s in ("3 + 4", "2 ** 10", "3 hypot 4")] == [7.0, 1024.0, 5.0] and hypot(6, 8) == 10.0
assert raises(ValueError, calc, "1 ? 2") == "unknown operator '?'; choose from ['*', '**', '+', '-', '/', 'hypot']"

def make_counter(start=0):
    count = start
    def step(inc=1):
        nonlocal count
        count += inc
        return count
    return step

def broken_counter():
    count = 0
    def step():
        count += 1  # noqa: F823
        return count
    return step

ca, cb = make_counter(), make_counter(100)
assert [ca(), ca(), ca(10), cb(), ca()] == [1, 2, 12, 101, 13]
cell = ca.__closure__[0]
assert (ca.__code__.co_freevars, cell.cell_contents, cb.__closure__[0].cell_contents) == (("count",), 13, 101)
cell.cell_contents = 0                                                              # cells are writable
assert ca() == 1 and "cannot access local variable 'count'" in raises(UnboundLocalError, broken_counter())

rng = random.Random(0)                                                              # sparse polynomial vs NumPy (PY-2)
dense = lambda d: [d.get(k, 0) for k in range(max(d) + 1)]  # noqa: E731
for _ in range(200):
    p = {k: c for k, c in zip(rng.sample(range(8), 4), rng.sample(range(-5, 6), 4)) if c}
    q = {k: c for k, c in zip(rng.sample(range(8), 3), rng.sample(range(-5, 6), 3)) if c}
    if p and q:
        assert dense(numtools.poly_mul(p, q)) == [int(v) for v in np.polynomial.polynomial.polymul(dense(p), dense(q))]
        assert numtools.poly_eval(p, 2) == np.polynomial.polynomial.polyval(2, dense(p))
assert numtools.poly_deriv({3: 2, 1: 5, 0: 7}) == {2: 6, 0: 5} and numtools.poly_add({1: 2}, {1: -2}) == {}
print("Part 2: functions-as-values checks passed.")

# %% [markdown]
# ## Part 3: Decorators and caching (notes 3.5, 3.6, Example 4.3)

# %%
@timer
def slow_sum(n):
    """Sum of squares below n."""
    return sum(i * i for i in range(n))

assert [slow_sum(1000) for _ in range(3)] == [332833500] * 3
assert (slow_sum.__name__, slow_sum.__doc__, len(slow_sum.timings), slow_sum.__wrapped__.__name__) == (
    "slow_sum", "Sum of squares below n.", 3, "slow_sum")

@timer
def explode():
    raise KeyError("x")

assert raises(KeyError, explode) is not None and len(explode.timings) == 1          # try/finally records failures

def bare_timer(func):                                                                # no functools.wraps
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper

assert (bare_timer(slow_sum.__wrapped__).__name__, bare_timer(slow_sum.__wrapped__).__doc__) == ("wrapper", None)
calls = {"plain": 0, "memo": 0}

def fib_plain(n):
    calls["plain"] += 1
    return n if n < 2 else fib_plain(n - 1) + fib_plain(n - 2)

@memoize
def fib_memo(n):
    calls["memo"] += 1
    return n if n < 2 else fib_memo(n - 1) + fib_memo(n - 2)

assert (fib_plain(25), fib_memo(25), calls, len(fib_memo.cache)) == (75025, 75025, {"plain": 242785, "memo": 26}, 26)
assert fib_memo.__name__ == "fib_memo" and raises(TypeError, memoize(sum), [1, 2]) == "unhashable type: 'list'"
total = functools.lru_cache(maxsize=None)(sum)
assert total((1, 2, 3)) == total((1, 2, 3)) == 6 and str(total.cache_info()) == "CacheInfo(hits=1, misses=1, maxsize=None, currsize=1)"
assert raises(TypeError, total, [1, 2, 3]) == "unhashable type: 'list'"

@functools.cache
def zeros(n):
    return [0] * n

zeros(3).append(9)
assert zeros(3) == [0, 0, 0, 9]                                                      # the aliasing trap
pair, pair_typed = functools.lru_cache(maxsize=None)(lambda x, y: (x, y)), functools.lru_cache(typed=True)(lambda x, y: (x, y))
for args in ((1, 2), (1.0, 2.0)):
    pair(*args), pair_typed(*args)
assert pair.cache_info().currsize == 1 and pair_typed.cache_info().currsize == 2     # equal keys share an entry

attempts = []

def flaky_body():
    attempts.append(1)
    if len(attempts) < 3:
        raise ConnectionError(f"attempt {len(attempts)} failed")
    return "ok"

for n_try, ok in ((3, True), (2, False)):                                            # Section 3.5: flaky() succeeds on attempt 3
    attempts.clear()
    flaky = retry(n_try, exceptions=(ConnectionError,))(flaky_body)
    assert (flaky() == "ok") if ok else raises(ConnectionError, flaky) == "attempt 2 failed"
    assert len(attempts) == n_try
assert raises(ValueError, retry, 0) == "n must be at least 1"
assert raises(KeyError, retry(3, exceptions=(ConnectionError,))(explode.__wrapped__)) is not None   # not retried

# Example 4.3: timer applied to the PY-1 solvers, 200 runs each, median time
timed = {"bisect": timer(bisect), "newton": timer(newton), "secant": timer(secant)}
call = {"bisect": lambda s, f, df, a, b: s(f, a, b, tol=1e-10), "newton": lambda s, f, df, a, b: s(f, df, b, tol=1e-10),
        "secant": lambda s, f, df, a, b: s(f, a, b, tol=1e-10)}
iteration_counts = {"bisect": [34, 34, 34], "newton": [5, 5, 6], "secant": [6, 7, 7]}
print(f"{'function':>13} {'solver':>7} {'iterations':>10} {'median us':>10}")
for k, (cname, f, df, a, b) in enumerate(CASES):
    medians = {}
    for sname, solver in timed.items():
        solver.timings.clear()
        for _ in range(200):
            root, iterations = call[sname](solver, f, df, a, b)
        medians[sname] = statistics.median(solver.timings)
        assert iterations == iteration_counts[sname][k] and len(solver.timings) == 200
        print(f"{cname:>13} {sname:>7} {iterations:>10} {medians[sname] * 1e6:>10.1f}")
    assert medians["newton"] < medians["bisect"] and medians["secant"] < medians["bisect"]
assert math.ceil(math.log2(1e10)) == 34
print("Part 3: decorator and caching checks passed.")

# %% [markdown]
# ## Part 4: Recursion (notes 3.7, Examples 4.4 and 4.7)

# %%
def fib_calls(n):
    return 1 if n < 2 else fib_calls(n - 1) + fib_calls(n - 2) + 1

def fib_iter(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

assert all(fib_calls(n) == 2 * fib_iter(n + 1) - 1 for n in range(26))               # C(n) = 2 F(n+1) - 1
assert [fib_calls(n) for n in (20, 25, 30)] == [21891, 242785, 2692537] and fib_iter(21) == 10946
assert abs(fib_calls(25) / fib_calls(20) - PHI**5) < 0.05 and round(PHI**5, 2) == 11.09

def make_cached_fib():
    body_runs = []

    @functools.cache
    def fibc(n):
        body_runs.append(n)
        return n if n < 2 else fibc(n - 1) + fibc(n - 2)

    return fibc, body_runs

for n in (4, 10, 30):                                                                # n + 1 bodies, 2n - 1 calls
    fibc, runs = make_cached_fib()
    value, info = fibc(n), fibc.cache_info()
    assert value == fib_iter(n) and len(runs) == n + 1 and info.misses == n + 1 and info.hits + info.misses == 2 * n - 1

@functools.lru_cache(maxsize=None)
def fib(n):
    return n if n < 2 else fib(n - 1) + fib(n - 2)

assert (fib(30), str(fib.cache_info())) == (832040, "CacheInfo(hits=28, misses=31, maxsize=None, currsize=31)")
assert fib(100) == 354224848179261915075 and str(fib.cache_info()) == "CacheInfo(hits=99, misses=101, maxsize=None, currsize=101)"
fib.cache_clear()
assert raises(RecursionError, fib, 2000) is not None                                   # cold start: depth 2000
fib.cache_clear()
for k in range(2001):
    fib(k)
assert len(str(fib(2000))) == 418 and str(fib.cache_info()) == "CacheInfo(hits=3999, misses=2001, maxsize=None, currsize=2001)"

def depth(n):
    return 0 if n == 0 else 1 + depth(n - 1)

assert sys.getrecursionlimit() == 1000 and depth(900) == 900
assert raises(RecursionError, depth, 5000) == "maximum recursion depth exceeded"

def naive(n):
    return n if n < 2 else naive(n - 1) + naive(n - 2)

@memoize
def memo(n):
    return n if n < 2 else memo(n - 1) + memo(n - 2)

def cold_memo(n):
    memo.cache.clear()
    return memo(n)

ns = list(range(10, 28))
t_naive, t_memo = [best(naive, n, 3) for n in ns], [best(cold_memo, n, 30) for n in ns]
ratio = t_naive[ns.index(25)] / t_naive[ns.index(20)]
print(f"naive time(25)/time(20) = {ratio:.1f}   (phi^5 = {PHI**5:.2f})")
assert PHI**5 / 2 < ratio < PHI**5 * 2 and t_memo[-1] < t_naive[-1] / 100 and cold_memo(20) == naive(20)
save_line_plot("fibonacci_timing.png", "fib(n): the call tree against the cache", ns, {
    "naive (exponential)": (t_naive, "o-"), "memoized (linear)": (t_memo, "s-"),
    "$\\varphi^n$ reference": ([t_naive[0] * PHI ** (n - ns[0]) for n in ns], "k--")}, "seconds (best of repeats)")

def lower_bound(a, x, check=False):
    lo, hi, steps = 0, len(a), 0
    while lo < hi:
        if check:
            assert all(v < x for v in a[:lo]) and all(v >= x for v in a[hi:])             # the invariant
        steps += 1
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo, steps

def lower_bound_rec(a, x, lo=0, hi=None):
    hi = len(a) if hi is None else hi
    if lo == hi:
        return lo
    mid = (lo + hi) // 2
    if a[mid] < x:
        return lower_bound_rec(a, x, mid + 1, hi)
    return lower_bound_rec(a, x, lo, mid)

a = [1, 3, 3, 5, 8, 8, 8, 13]
assert [lower_bound(a, x) for x in (0, 3, 4, 14)] == [(0, 4), (1, 4), (3, 3), (8, 3)] and math.floor(math.log2(len(a))) + 1 == 4
rng = random.Random(0)
for trial in range(3000):
    arr = sorted(rng.randint(0, 50) for _ in range(rng.randint(1, 1000)))
    x = rng.randint(-5, 55)
    i, steps = lower_bound(arr, x, check=trial < 300)
    assert i == lower_bound_rec(arr, x) == bisect_left(arr, x) and steps <= math.floor(math.log2(len(arr))) + 1
assert lower_bound([], 1) == (0, 0)

def merge_sort(xs, key=lambda v: v, counter=None):
    if len(xs) <= 1:
        return list(xs)
    mid = len(xs) // 2
    left, right = merge_sort(xs[:mid], key, counter), merge_sort(xs[mid:], key, counter)
    out, i, j = [], 0, 0
    while i < len(left) and j < len(right):
        if counter is not None:
            counter[0] += 1
        if key(left[i]) <= key(right[j]):                                                  # <= keeps the sort stable
            out.append(left[i])
            i += 1
        else:
            out.append(right[j])
            j += 1
    return out + left[i:] + right[j:]

rng = random.Random(1)
data = [rng.randint(0, 999) for _ in range(10_000)]
pairs = [(rng.randint(0, 20), idx) for idx in range(2000)]
count = [0]
merge_sort([rng.random() for _ in range(1024)], counter=count)
assert merge_sort([5, 2, 9, 1, 5, 6, 0]) == [0, 1, 2, 5, 5, 6, 9] and merge_sort([]) == [] and merge_sort(data) == sorted(data)
assert merge_sort(pairs, key=lambda p: p[0]) == sorted(pairs, key=lambda p: p[0]) and count[0] <= 1024 * math.ceil(math.log2(1024))
t_ms, t_sorted = best(merge_sort, data, 3), best(sorted, data, 3)
print(f"merge sort {t_ms * 1e3:.1f} ms vs sorted {t_sorted * 1e3:.2f} ms for 10^4 items; {count[0]} comparisons at n=1024")
assert t_sorted < t_ms

mults = [0]

def det_cofactor(A):
    n = len(A)
    if n == 1:
        return A[0][0]
    total = 0
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in A[1:]]
        mults[0] += 1
        total += (-1) ** j * A[0][j] * det_cofactor(minor)
    return total

def T(n):
    return sum(math.factorial(n) // math.factorial(k) for k in range(1, n))               # LA-6 section 3.5

assert det_cofactor([[0, 2, 1], [1, 1, 3], [2, -1, 1]]) == 7 and [T(n) for n in range(2, 8)] == [2, 9, 40, 205, 1236, 8659]
assert all(T(n) == n * T(n - 1) + n for n in range(2, 12)) and 4.1e18 < T(20) < 4.3e18
assert abs(T(20) / (math.factorial(20) * (math.e - 1)) - 1) < 1e-12
rng = np.random.default_rng(0)
det_ns, t_cof, t_np = list(range(2, 9)), [], []
for n in det_ns:
    M = rng.integers(-3, 4, size=(n, n))
    mults[0] = 0
    ref = np.linalg.det(M.astype(float))
    assert abs(det_cofactor(M.tolist()) - ref) < 1e-8 * max(1.0, abs(ref)) and mults[0] == T(n)
    t_cof.append(best(det_cofactor, M.tolist(), 2))
    t_np.append(best(np.linalg.det, M.astype(float), 20))
print("cofactor seconds, n=2..8:", " ".join(f"{t:.1e}" for t in t_cof))
assert t_cof[-1] / t_cof[-3] > 15 and t_cof[-1] > 100 * t_np[-1]                           # factorial growth (about 7 * 8)
save_line_plot("determinant_growth.png", "Determinant: recursion against elimination", det_ns, {
    "cofactor recursion (pure Python)": (t_cof, "o-"), "numpy.linalg.det (LU)": (t_np, "s-"),
    "$\\propto n!\\sum 1/k!$": ([T(n) * t_cof[0] / T(2) for n in det_ns], "k--")})

def is_even(n):
    return True if n == 0 else is_odd(n - 1)

def is_odd(n):
    return False if n == 0 else is_even(n - 1)

def count_to(n):
    i = 0
    while i < n:
        i += 1
        yield i

def flatten(x):
    for item in x:
        if isinstance(item, list):
            yield from flatten(item)
        else:
            yield item

gen = count_to(3)
assert (is_even(10), is_odd(7), is_even(7)) == (True, True, False) and (next(gen), next(gen)) == (1, 2)
assert gen.gi_frame.f_locals == {"n": 3, "i": 2} and list(flatten([1, [2, [3, [4, 5]], 6], [], 7])) == [1, 2, 3, 4, 5, 6, 7]
print("Part 4: recursion checks passed.")

# %% [markdown]
# ## Part 5: Text and the standard library (notes 3.9, Example 4.5)

# %%
LINE = re.compile(r"""
    (?P<ts>\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})   # ISO timestamp
    \s+
    sensor_(?P<id>\d+)                            # sensor number
    \s+
    (?P<value>-?\d+(?:\.\d+)?)                    # signed decimal
    """, re.VERBOSE)

def parse_log(lines):
    """Return (records, rejected line indexes); a record is (datetime, sensor id, value)."""
    records, rejected = [], []
    for idx, text in enumerate(lines):
        m = LINE.fullmatch(text)
        if m is None:
            rejected.append(idx)
        else:
            records.append((datetime.fromisoformat(m["ts"]), int(m["id"]), float(m["value"])))
    return records, rejected

line = "2024-03-01T12:00:00 sensor_7 23.5"
assert LINE.fullmatch(line).groupdict() == {"ts": "2024-03-01T12:00:00", "id": "7", "value": "23.5"}
assert parse_log([line])[0] == [(datetime(2024, 3, 1, 12, 0), 7, 23.5)]
assert re.sub(r"sensor_(\d+)", r"S\1", line) == "2024-03-01T12:00:00 S7 23.5"
small = ["2024-03-01T12:00:00 sensor_7 23.5", "2024-03-01T12:00:10 sensor_7 23.9", "2024-03-01T12:00:10 sensor_12 -4.25",
         "garbage line", "2024-03-01T12:00:20 sensor_7 ERR", "2024-03-01T12:00:30 sensor_12 -4.0"]
recs, bad = parse_log(small)
vals = defaultdict(list)
for _, sid, value in recs:
    vals[sid].append(value)
assert [i + 1 for i in bad] == [4, 5] and {s: statistics.mean(v) for s, v in vals.items()} == {7: 23.7, 12: -4.125}

rng = random.Random(0)                                                                    # 1000 generated lines
t0, truth, log, bad_lines = datetime(2024, 3, 1, 12, 0, 0), [], [], []
for idx in range(1000):
    if idx % 13 == 5:
        bad_lines.append(idx)
        log.append(["garbage line", "2024-03-01T12:00:20 sensor_7 ERR", "sensor_3 4.5", ""][(idx // 13) % 4])
        continue
    ts, sid, value = t0 + timedelta(seconds=10 * idx), rng.randint(1, 12), round(rng.gauss(20, 5), 2)
    truth.append((ts, sid, value))
    log.append(f"{ts.isoformat()} sensor_{sid} {value:.2f}")
parsed, rejected = parse_log(log)
assert rejected == bad_lines and len(parsed) == 1000 - len(bad_lines) and parsed == truth
means = lambda rows: {s: statistics.mean(v for _, sid, v in rows if sid == s) for s in range(1, 13)}  # noqa: E731
assert means(parsed) == means(truth) and re.match(r"\d+", "apples 42") is None
assert re.search(r"\d+", "apples 42").group() == "42" and re.findall(r"\d+", "a1 b22 c333") == ["1", "22", "333"]
assert re.match(r"(\w+)=(\d+)", "alpha=7").groups() == ("alpha", "7") and (len("\n"), len(r"\n")) == (1, 2)

assert (math.hypot(3, 4), math.comb(5, 2), math.isclose(0.1 + 0.2, 0.3), cmath.sqrt(-4)) == (5.0, 10, True, 2j)
assert Fraction(1, 3) + Fraction(1, 6) == Fraction(1, 2) and Decimal("0.1") + Decimal("0.2") == Decimal("0.3") != 0.1 + 0.2
data = [2, 4, 4, 4, 5, 5, 7, 9]
assert (statistics.mean(data), statistics.median(data), statistics.pstdev(data), round(statistics.stdev(data), 4)) == (5, 4.5, 2.0, 2.1381)
r1, r2 = random.Random(42), random.Random(42)
assert [r1.randint(1, 6) for _ in range(5)] == [r2.randint(1, 6) for _ in range(5)]
d = datetime.fromisoformat("2024-03-01T12:00:00")
assert d + timedelta(days=30) == datetime(2024, 3, 31, 12, 0) and (datetime(2024, 3, 1) - datetime(2024, 2, 1)).days == 29
assert d.replace(tzinfo=timezone.utc).timestamp() == 1709294400.0 and os.environ.get("NO_SUCH_VARIABLE", "unset") == "unset"
assert textwrap.fill("Python is a language that lets you work quickly and integrate systems.", width=30) == \
    "Python is a language that lets\nyou work quickly and integrate\nsystems."
assert pprint.pformat({"gamma": "a string", "alpha": list(range(10)), "beta": {"x": 1, "y": [1, 2, 3]}}, width=50) == \
    "{'alpha': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],\n 'beta': {'x': 1, 'y': [1, 2, 3]},\n 'gamma': 'a string'}"
print("Part 5: text and standard-library checks passed.")

# %% [markdown]
# ## Your turn
#
# 1. **A call counter.** Write `counted(func)` that adds a `.calls` attribute (with `functools.wraps`) and
#    use it to count the calls of naive `fib(10)`. Predict the count from 2 F(n+1) - 1 first.
# 2. **Hit rate against `maxsize`.** Run `functools.lru_cache(maxsize=m)` on the recursion of `fib(30)` for
#    several `m` and tabulate `hits / (hits + misses)`. Predict the shape before you run it.
# 3. **Explicit stack.** Rewrite `sum_nested` (recursive sum of a nested list) with a list as a stack and
#    run it on a list nested 50,000 deep, where the recursion fails.
# 4. **A new subcommand.** Add `secant EXPR X0 X1` to `numtools/cli.py` (a patched copy is tested below).
#
# Reference solutions follow, each asserted. Try first.

# %%
def counted(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        wrapper.calls += 1
        return func(*args, **kwargs)
    wrapper.calls = 0
    return wrapper

@counted
def fib_counted(n):
    return n if n < 2 else fib_counted(n - 1) + fib_counted(n - 2)

assert fib_counted(10) == 55 and fib_counted.calls == 2 * 89 - 1 == 177 and fib_counted.__name__ == "fib_counted"

measured = {}
for size in (1, 2, 3, 4, 8, None):
    @functools.lru_cache(maxsize=size)
    def fibm(n):
        return n if n < 2 else fibm(n - 1) + fibm(n - 2)
    fibm(30)
    info = fibm.cache_info()
    measured[size] = (info.hits, info.misses)
    print(f"maxsize={size!s:>4}: hits {info.hits:>4}, misses {info.misses:>7}, hit rate {info.hits / (info.hits + info.misses):.3f}")
# maxsize=1 never hits and repeats all 2,692,537 calls; maxsize=2 thrashes; three slots already give the
# unbounded cache's 28 hits and 31 misses, because fib(n) only ever needs its two most recent results.
assert measured[1] == (0, 2692537) and measured[2] == (8656, 41641)
assert all(measured[s] == (28, 31) for s in (3, 4, 8, None))

def sum_nested(x):
    return sum(sum_nested(item) if isinstance(item, list) else item for item in x)

def sum_nested_iter(x):
    total, stack = 0, [x]
    while stack:
        for item in stack.pop():
            if isinstance(item, list):
                stack.append(item)
            else:
                total += item
    return total

deep = [1]
for _ in range(50_000):
    deep = [1, deep]
assert raises(RecursionError, sum_nested, deep) is not None and sum_nested_iter(deep) == 50_001
assert sum_nested([1, [2, [3, [4]]], 5]) == sum_nested_iter([1, [2, [3, [4]]], 5]) == 15

PATCHES = [("from .roots import bisect, newton\n", "from .roots import bisect, newton, secant\n"),
           ("def build_parser()", "def run_secant(args: argparse.Namespace) -> tuple[float, int]:\n"
            "    return secant(make_function(args.expr), args.x0, args.x1, tol=args.tol)\n\n\ndef build_parser()"),
           ('metavar="{bisect,newton}"', 'metavar="{bisect,newton,secant}"'),
           ("    return parser\n", '    ps = sub.add_parser("secant", help="secant method from two starting points")\n'
            '    common(ps)\n    ps.add_argument("x0", type=float, help="first starting point")\n'
            '    ps.add_argument("x1", type=float, help="second starting point")\n    ps.set_defaults(func=run_secant)\n'
            "    return parser\n")]
with tempfile.TemporaryDirectory() as tmp:                                                  # a patched copy of the package
    shutil.copytree(HERE / "numtools", Path(tmp) / "numtools", ignore=shutil.ignore_patterns("__pycache__"))
    cli_path = Path(tmp) / "numtools" / "cli.py"
    src = cli_path.read_text()
    for old, new in PATCHES:
        assert src.count(old) == 1, old
        src = src.replace(old, new)
    cli_path.write_text(src)
    proc = run_py("-m", "numtools", "secant", "x*x-2", "1", "2", "-v", cwd=tmp)
    assert proc.stdout.splitlines()[:2] == ["root = 1.41421356237", "iterations = 7"]
    assert "{bisect,newton,secant}" in run_py("-m", "numtools", "--help", cwd=tmp).stdout
print("Your-turn reference solutions passed.")

# %%
if __name__ == "__main__":
    print("All checks passed.")
