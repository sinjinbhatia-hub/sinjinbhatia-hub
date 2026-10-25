# PY-5: Errors, Testing, Debugging, and Type Hints — Writing Code You Can Trust

**Time:** ~8 hours · **Prerequisites:** PY-4 · **Unlocks:** PY-6 onward (every later project is written test-first); strongly recommended before LA-2's project

> **Conventions in this module.** Python 3.11+. Every output below was pasted by a script from a real run (Python 3.13, pytest 9.1, mypy 2.4, ruff). Three liberties: tracebacks omit the `~~^^^` markers (`PYTHONNODEBUGRANGES=1`) and show file names, not paths; pytest's `platform`/`rootdir`/`plugins` header lines are dropped; in pdb transcripts the commands are re-inserted after each `(Pdb)` prompt (stdin was piped) and pdb's own frames are marked omitted. Failure diffs are verbatim. `$` lines are shell commands. The code under test is `integrators.py`, its tests `test_integrators.py`, both in this folder (`pip install pytest mypy ruff`; `project.py` runs without the last two). A **public** function has no leading underscore. A **panel** is a subinterval of width $h = (b-a)/n$; Simpson's rule needs an even number $n$ of them.

## 0. Where this sits

PY-1 to PY-4 taught you to predict what Python does. This module teaches your code to say so when it cannot do what you asked, and gives you the means to check that it does what you think. You need PY-4's classes (exceptions are classes, PY-4 §3.10; `Vec` of PY-4 §3.3 is annotated in §3.11, and its `Evaluable` protocol is checked for real) and PY-3's decorators and generators (`@contextmanager`). It makes possible every project from PY-6 on (each ships a `test_*.py`), LA-2's LU project, CM-1's conservation-law tests, and ML-2's gradient checks.

## 1. The problem that forced this idea

Here is the Simpson rule of PY-4 §3.10 with one character wrong, integrating $\sin$ over $[0, \pi]$ (exact value 2).

```python
import math


def simpson(f, a, b, n):
    """Composite Simpson rule, n even. (Your PY-4 exercise, with one character wrong.)"""
    h = (b - a) / n
    ys = [f(a + i * h) for i in range(n + 1)]
    return h / 3 * (ys[0] + 4 * sum(ys[1:-2:2]) + 2 * sum(ys[2:-1:2]) + ys[-1])


for n in (12, 24, 48, 96):
    value = simpson(math.sin, 0.0, math.pi, n)
    print(f"n = {n:3d}   integral = {value:.6f}   error = {value - 2.0:+.2e}")
```
```text
n =  12   integral = 1.909708   error = -9.03e-02
n =  24   integral = 1.977222   error = -2.28e-02
n =  48   integral = 1.994293   error = -5.71e-03
n =  96   integral = 1.998572   error = -1.43e-03
```

Nothing was raised or printed. Every number is plausible: at $n = 48$ the answer is 1.994293, 0.29% low, and the sequence converges as $n$ grows. The one clue is that the error shrinks by about 4 per doubling of $n$ where Simpson's rule gives 16 (PY-4 §3.10's table). Unless you look for it, the 0.3% travels into everything downstream.

**A numerical routine that returns a plausible wrong number is worse than one that crashes.** A crash costs minutes; a wrong number costs a week, after the simulation or paper rests on it. Python cannot help: `ys[1:-2:2]` is a legal slice. Four tools close the gap:

- **Exceptions** (§3.1–3.7): the code announces its own failure instead of returning a guess.
- **Tests** (§3.10): small programs that pin down what the code must do and fail loudly otherwise.
- **A debugger and logging** (§3.8–3.9): ways to see what the code is doing.
- **Type hints** (§3.11): a checker that catches wrong-type mistakes before anything runs.

The forcing problem: make the off-by-one above impossible to miss. By §3.10 a test fails the moment it is planted and names the missing sample.

## 2. Intuition

**An exception is a signal that unwinds the call stack until something catches it.** Each call pushes a frame; `raise` stops the current line and pops frames, running `finally` blocks on the way, until a frame has a matching `except`.

```
call stack (newest at top)        what the exception does as it climbs
┌───────────────────────────┐
│ read_sample("x3")         │  float("x3") raises ValueError. No try here: frame discarded.  ①
├───────────────────────────┤
│ mean_of([...])            │  has try/finally: the finally block runs, frame discarded.     ②
├───────────────────────────┤
│ report([...])             │  no try: frame discarded.                                      ③
├───────────────────────────┤
│ <module>                  │  try: ... except ValueError: ← matches. Unwinding stops.       ④
└───────────────────────────┘
```

**A test is a small program whose only job is to fail loudly.** **A debugger pauses the program** so you can inspect the frames a traceback only names. **A log** is what the program writes about itself. **A type hint is a comment the tools can check:** Python ignores it; a checker reads it before the program runs.

**Guess before §3.10.** Simpson's rule on $[-h, h]$ is $\frac h3[f(-h) + 4f(0) + f(h)]$. Try $f = 1, x, x^2, x^3$ against the exact integrals. If the rule is exact on all four, what test would have caught the bug of §1 without knowing anything about $\sin$?

## 3. Formal development

### 3.1 Exceptions: what they are and how they climb

1. An **exception** is an instance of a class derived from `BaseException`. `raise obj` (or `raise Cls`, which calls `Cls()`) stops the current line. Python looks for a matching handler in the current frame; failing that it discards the frame (running its `finally` block) and looks in the caller, and so on. Falling off the bottom prints the traceback and exits with status 1 (Example 4.1).
2. `except C as err` matches when `isinstance(exc, C)`, so subclasses match; `C` may be a tuple. Clauses are tried top to bottom and the first match wins: write the subclass first. After the clause, the name `err` is deleted.

▶ **Try it first.** Predict the parent chain of `KeyError`, and whether a bare `except:` catches `KeyboardInterrupt`.

```python
print(sorted(cls.__name__ for cls in BaseException.__subclasses__()))
for cls in (ZeroDivisionError, KeyError, FileNotFoundError, StopIteration, KeyboardInterrupt):
    print(" -> ".join(c.__name__ for c in cls.__mro__[:-1]))
```
```text
['BaseExceptionGroup', 'Exception', 'GeneratorExit', 'KeyboardInterrupt', 'SystemExit']
ZeroDivisionError -> ArithmeticError -> Exception -> BaseException
KeyError -> LookupError -> Exception -> BaseException
FileNotFoundError -> OSError -> Exception -> BaseException
StopIteration -> Exception -> BaseException
KeyboardInterrupt -> BaseException
```

■ **What just happened.** `KeyError` and `IndexError` share `LookupError`; `FileNotFoundError` is an `OSError`. `KeyboardInterrupt`, `SystemExit` and `GeneratorExit` sit beside `Exception`, not under it, because they mean "stop", not "failed". So `except Exception` leaves Ctrl-C working; a bare `except:` is `except BaseException` and does not.

### 3.2 `try` / `except` / `else` / `finally`: control flow

▶ **Try it first.** Write the list returned for `x = 4`, `0`, `5`; then say what `x = "a"` does and when `finally` prints.

```python
def attempt(x):
    log = []
    try:
        log.append("try")
        y = 10 / x
        if y == 2:
            return log + ["return inside try"]
    except ZeroDivisionError:
        log.append("except")
    else:
        log.append("else")
    finally:
        print(f"  (finally runs for x = {x!r})")
    return log


for x in (4, 0, 5):
    print(f"x = {x}: {attempt(x)}")
try:
    attempt("a")
except TypeError as err:
    print("x = 'a': TypeError escaped:", err)
```
```text
  (finally runs for x = 4)
x = 4: ['try', 'else']
  (finally runs for x = 0)
x = 0: ['try', 'except']
  (finally runs for x = 5)
x = 5: ['try', 'return inside try']
  (finally runs for x = 'a')
x = 'a': TypeError escaped: unsupported operand type(s) for /: 'int' and 'str'
```

1. The `try` body runs top to bottom; at the first exception the rest is skipped.
2. A matching `except` runs, then control continues *after the whole statement*. No match: the exception keeps climbing (`x = "a"`).
3. `else` runs only if the body finished normally (no exception, no `return`/`break`/`continue`), and exceptions raised in it are *not* caught by this statement's `except` clauses. So put into `try` only the line that can fail; the handler then cannot mask a bug elsewhere.
4. `finally` runs on every way out: normal, after a handled exception, while an unhandled one climbs, and on `return`/`break`/`continue`. The return value is computed *before* it runs (`x = 5`).
5. A `return` in `finally` replaces whatever was leaving, even a climbing exception: `try: raise ValueError` / `finally: return "x"` returns `"x"` and the error vanishes (`project.py` asserts it). Never write one.

### 3.3 Catching the right thing

▶ **Try it first.** One function returns a plausible wrong answer. Which, and what is the right mean?

```python
def mean_reading(lines):
    values = []
    for line in lines:
        try:
            values.append(float(line))
        except Exception:
            pass
    return sum(values) / len(values)


def mean_reading_strict(lines):
    values = []
    for number, line in enumerate(lines, start=1):
        try:
            values.append(float(line))
        except ValueError as err:
            raise ValueError(f"line {number}: cannot read {line!r} as a number") from err
    return sum(values) / len(values)


lines = ["1.0", "2.0", "3,5", "4.0"]  # one European decimal comma
print(mean_reading(lines))  # what the mean of 1, 2, 3.5, 4 should be: 2.625
try:
    mean_reading_strict(lines)
except ValueError as err:
    print("ValueError:", err)
```
```text
2.3333333333333335
ValueError: line 3: cannot read '3,5' as a number
```

`except Exception: pass` turned a crash into 2.3333 where the mean of 1, 2, 3.5, 4 is 2.625: the decimal comma in `3,5` was dropped without a word, the module's forcing problem in a new costume. The repair names what it expects (`ValueError`), says which line failed, and keeps the original with `from err` (§3.4).

- Catch the narrowest class naming a failure you can handle. A handler that can do nothing useful should not exist; let the exception climb. A bare `except:` is never right.
- `except Exception` has two legitimate homes: the edge of a loop or worker that must survive one bad item, and the `main()` of a script that prints one clear line and exits nonzero. Log the traceback there (`log.exception`, §3.9).

### 3.4 Raising, chaining, and your own exceptions

1. Put the *value* in the message: `raise ValueError(f"n must be even, got {n}")`. A bare `raise` inside `except` re-raises the current exception.
2. `raise New(...) from err` stores `err` as `New.__cause__`; the traceback says "the direct cause of the following exception" (Example 4.5). Raising inside an `except` *without* `from` still links the old one (`__context__`: "during handling of the above exception, another exception occurred"), which usually means the handler has a bug. `from None` hides the context.
3. **Custom exceptions** subclass `Exception`. Give a package **one base class**, so callers can write `except SolverError`, and a subclass per failure a caller may treat differently. Use built-ins when they fit: a bad argument is a `ValueError`; a `ConvergenceError` is *not*, because the arguments were fine and the algorithm failed. Carry the data a handler needs as attributes, and pass every constructor argument to `super().__init__` so `copy` and `pickle` (which rebuild from `args`) work.

▶ **Build it yourself (2): exceptions for the PY-1 root finders.** PY-1's Newton raised `RuntimeError("did not converge")`: no data, indistinguishable from any other `RuntimeError`. Write the hierarchy, then make Newton turn its `ZeroDivisionError` into a `ConvergenceError` with `raise ... from`. Predict which of the three problems succeed, which raise, and what `__cause__` is.

```python
class SolverError(Exception):
    ...


class ConvergenceError(SolverError):
    def __init__(self, message: str, iterations: int, residual: float) -> None:
        super().__init__(message, iterations, residual)
        self.message = message
        self.iterations = iterations
        self.residual = residual

    def __str__(self) -> str:
        return f"{self.message} (iterations={self.iterations}, residual={self.residual:.3g})"
```

```python
from integrators import ConvergenceError


def newton(f, df, x0, tol=1e-12, max_iter=50):
    """Newton's method for f(x) = 0; raises ConvergenceError (never RuntimeError)."""
    x = x0
    for k in range(1, max_iter + 1):
        try:
            step = f(x) / df(x)
        except ZeroDivisionError as err:
            raise ConvergenceError("derivative is zero", k, abs(f(x))) from err
        x -= step
        if abs(step) < tol:
            return x
    raise ConvergenceError("no convergence", max_iter, abs(f(x)))
```
```python
from integrators import SolverError
from rootfind import newton

problems = {
    "x^2 - 2, from 1": (lambda x: x * x - 2, lambda x: 2 * x, 1.0),
    "x^2 - 2, from 0": (lambda x: x * x - 2, lambda x: 2 * x, 0.0),
    "x^3 - 2x + 2, from 0": (lambda x: x**3 - 2 * x + 2, lambda x: 3 * x * x - 2, 0.0),
}
for name, (f, df, x0) in problems.items():
    try:
        print(f"{name:21s} root = {newton(f, df, x0):.10f}")
    except SolverError as err:
        print(f"{name:21s} {type(err).__name__}: {err}; cause: {err.__cause__!r}")
```
```text
x^2 - 2, from 1       root = 1.4142135624
x^2 - 2, from 0       ConvergenceError: derivative is zero (iterations=1, residual=2); cause: ZeroDivisionError('float division by zero')
x^3 - 2x + 2, from 0  ConvergenceError: no convergence (iterations=50, residual=2); cause: None
```

■ **What just happened.** From $x_0 = 0$ the derivative is zero and `ZeroDivisionError` became a `ConvergenceError` whose `__cause__` keeps the original. The cubic from $x_0 = 0$ cycles $0 \to 1 \to 0 \to \cdots$ for 50 iterations, with no cause because nothing raised. One `except SolverError` catches both and reads `iterations` and `residual` as numbers.

### 3.5 EAFP, LBYL, and exceptions that are not errors

Python prefers **EAFP** ("easier to ask forgiveness than permission": try it, handle failure) to **LBYL** ("look before you leap": test first). Neither always wins.

▶ **Try it first.** Which lines print `None`, which raise, and what does `LBYL still failed` mean?

```python
from pathlib import Path

import integrators as ig

RULES = {"trapezoid": ig.trapezoid, "simpson": ig.simpson}

# 1. Dictionary. A missing key that is NORMAL: ask with .get. A missing key that is a MISTAKE: try/except.
counts: dict[str, int] = {}
for word in "to be or not to be".split():
    counts[word] = counts.get(word, 0) + 1
print(counts)


def pick_rule(name):
    try:
        return RULES[name]
    except KeyError:
        raise ValueError(f"unknown rule {name!r}; choose from {sorted(RULES)}") from None


print(pick_rule("simpson").__name__)
try:
    pick_rule("simpsons")
except ValueError as err:
    print("ValueError:", err)

# 2. File. Looking first leaves a gap in which the world can change.
path = Path("config.txt")
path.write_text("tol = 1e-8\n")


def read_lbyl(p):
    if p.exists():  # look ...
        p.unlink()  # (meanwhile another process deletes the file)
        return p.read_text()  # ... and leap
    return None


def read_eafp(p):
    try:
        return p.read_text()
    except FileNotFoundError:
        return None


try:
    read_lbyl(path)
except FileNotFoundError as err:
    print("LBYL still failed:", type(err).__name__)
print("EAFP on the missing file:", read_eafp(path))
```
```text
{'to': 2, 'be': 2, 'or': 1, 'not': 1}
simpson
ValueError: unknown rule 'simpsons'; choose from ['simpson', 'trapezoid']
LBYL still failed: FileNotFoundError
EAFP on the missing file: None
```

1. **A missing key that is normal**: `.get`, `defaultdict` or `Counter`. **A missing key that is a mistake**: `try`/`except KeyError`, re-raising a domain error that lists the valid choices (`pick_rule`; `from None` because the `KeyError` adds nothing).
2. **When the world can change between check and use** (files, sockets, threads), LBYL is wrong, not merely slow: `exists()` then `read_text()` leaves a gap. Do the operation and catch the *specific* exception (`FileNotFoundError`, not `OSError`, not `Exception`).
3. **When the check is cheap and you cannot recover locally**, check first: `integrators._panels` tests `n % 2` before any work.

**Exceptions that are not errors.** A `for` loop is `iter`, `next`, and a catch of `StopIteration` (PY-2 §3.4). Build one:

```python
def my_for(iterable, body):
    """A for loop made of iter, next and StopIteration."""
    it = iter(iterable)
    while True:
        try:
            item = next(it)
        except StopIteration:
            break
        body(item)


my_for("abc", lambda ch: print("got", ch))


def gen():
    yield 1
    raise StopIteration  # a bug: a generator must return instead


try:
    list(gen())
except RuntimeError as err:
    print("RuntimeError:", err)
```
```text
got a
got b
got c
RuntimeError: generator raised StopIteration
```

A stray `StopIteration` inside a generator could silently end it early, so Python converts it to `RuntimeError` (PEP 479), as the last line shows; end a generator with `return`. For a default on an empty iterator, `next(it, default)` avoids the exception.

### 3.6 `assert` and `warnings`

▶ **Try it first.** Two functions reject a negative number, one with `assert`, one with `raise`. Predict their output normally and under `python3 -O`.

```python
def sqrt_assert(x):
    assert x >= 0, "x must be non-negative"  # WRONG tool: this is input validation
    return x**0.5


def sqrt_raise(x):
    if x < 0:
        raise ValueError(f"x must be non-negative, got {x}")
    return x**0.5


for fn in (sqrt_assert, sqrt_raise):
    try:
        print(fn.__name__, "->", fn(-4.0))
    except (AssertionError, ValueError) as err:
        print(fn.__name__, "->", type(err).__name__, err)
```
```console
$ python3 assert_demo.py
sqrt_assert -> AssertionError x must be non-negative
sqrt_raise -> ValueError x must be non-negative, got -4.0
$ python3 -O assert_demo.py
sqrt_assert -> (1.2246467991473532e-16+2j)
sqrt_raise -> ValueError x must be non-negative, got -4.0
```

`-O` removes every `assert`, so the validator stopped validating and `sqrt(-4)` returned a *complex number* with a rounding-error real part: a plausible wrong answer once more.

- `assert cond, msg` states an **internal invariant**, true if *your own code* is right (`assert len(ys) == n + 1`). If it fires, the bug is yours.
- **Input validation is `if ...: raise ValueError(...)`**: the caller erred, and the check must survive `-O`.
- `assert (cond, "msg")` is always true (a non-empty tuple); Python warns at compile time.

**Warnings** say "this still works, but..." without stopping. `warnings.warn(msg, Category, stacklevel=2)` blames the caller's line. `DeprecationWarning` is for developers (shown by default only when triggered from `__main__`; pytest displays it), `FutureWarning` for end users, `UserWarning` is the default, `RuntimeWarning` flags suspicious numerics. A filter can promote any warning to an exception (`warnings.simplefilter("error")`, `pytest -W error`).

NumPy's warnings are the ones you will meet most: `np.log(0.0)` warns (`RuntimeWarning: divide by zero encountered in log`) and returns `-inf`, and inside `with np.errstate(divide="raise"):` it raises `FloatingPointError`. A `DeprecationWarning` raised with `stacklevel=2` names the *caller's* line. `project.py` Part 4 asserts all three behaviors (PY-6 returns to NumPy's).

### 3.7 Context managers: cleanup that cannot be forgotten

`with EXPR as v: BODY` is approximately (the language reference, "The with statement", has the exact form):

```text
mgr = EXPR;  v = type(mgr).__enter__(mgr)
try:     BODY
except BaseException as e:
    if not type(mgr).__exit__(mgr, type(e), e, e.__traceback__): raise
else:    type(mgr).__exit__(mgr, None, None, None)
```

`__exit__` runs however the body ends; a true return value *suppresses* the exception (`contextlib.suppress`). `@contextmanager` turns a generator into a manager: code before `yield` is `__enter__`, the yielded value is `v`, code after is `__exit__`. An exception in the body is thrown *into the generator at the `yield`*, so cleanup needs `try/finally`.

▶ **Build it yourself (3): a timer twice, and a temporary directory.** Predict what the last two lines print.

```python
import shutil
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace


class Timer:
    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.elapsed = time.perf_counter() - self.start
        return False  # do not suppress an exception


@contextmanager
def timer():
    box = SimpleNamespace(elapsed=None)
    start = time.perf_counter()
    try:
        yield box
    finally:  # runs even if the body raised
        box.elapsed = time.perf_counter() - start


@contextmanager
def temp_directory():
    path = Path(tempfile.mkdtemp(prefix="py5_"))
    try:
        yield path
    finally:
        shutil.rmtree(path)


with Timer() as t1:
    time.sleep(0.02)
with timer() as t2:
    time.sleep(0.02)
print("both slept at least 20 ms:", t1.elapsed >= 0.02, t2.elapsed >= 0.02)

try:
    with temp_directory() as d:
        (d / "scratch.txt").write_text("x")
        print("inside the block: exists =", d.exists())
        raise RuntimeError("body failed")
except RuntimeError:
    print("after the failure: exists =", d.exists())
```
```text
both slept at least 20 ms: True True
inside the block: exists = True
after the failure: exists = False
```

■ **What just happened.** Class and generator do the same job; the generator is shorter, the class easier to extend. The body raised and the directory was still removed; without `try/finally` around the `yield` it would not have been. (Library versions: `tempfile.TemporaryDirectory`, `contextlib.chdir`.)

### 3.8 Debugging

**Method.** Reproduce first, with the smallest input that still fails (halve the input or the code and repeat: bisection; `git bisect` does it over history, PY-10). State one hypothesis, look at the value that confirms or kills it, fix, then pin the bug with a test. **Read the traceback bottom-up:** the last line is *what*, the frame above it *where it was noticed*, the rest *how we got there*. The cause is often earlier than the line blamed.

**`print` debugging, done right:** name, `repr`, one fact per line.

```python
tol = "1e-8 "  # as read from a config file
n = 8
err = 0.1 + 0.2 - 0.3

print("tol", tol, "n", n)  # which is which? is tol a string?
print(f"{tol=} {n=} {err=}")  # name and repr in one go
print(f"{err=:.3e} {len(tol)=}")  # format specs and expressions work too
```
```text
tol 1e-8  n 8
tol='1e-8 ' n=8 err=5.551115123125783e-17
err=5.551e-17 len(tol)=5
```

`tol='1e-8 '` exposes a *string with a trailing space* where a number belongs. Delete debugging prints before committing, or route them through `logging`. Explaining the problem aloud to a rubber duck works because it forces you to state your assumptions. **"It works in the notebook"**: cells run out of order and deleted cells leave their variables alive, so the script fails for lack of stale state. Restart and run all cells before you trust a notebook.

**pdb.** `breakpoint()` (PEP 553) calls `pdb.set_trace()`; `PYTHONBREAKPOINT=0` disables every call. `python3 -m pdb script.py` starts paused on line 1 and, if the program raises, stops in the failing frame (**post-mortem**); in a REPL `pdb.pm()` does that for the last exception; `pytest --pdb` opens it at each failing test. An empty line repeats the last command.

| Command | Does |
|---|---|
| `l` / `ll` | list source around the current line / the whole function |
| `n` / `s` | run the current line, stepping *over* calls / *into* a call |
| `c` / `r` | continue to the next breakpoint / until this function returns |
| `b 16, k == 1` | breakpoint at line 16 when `k == 1`; `b func` breaks at a function |
| `p expr` / `pp expr` / `a` | print / pretty-print / this function's arguments |
| `w` / `u` / `d` | stack (newest last) / up to the caller / back down |
| `q` | quit |

▶ **Build it yourself (4): a pdb session on a planted bug.** `next_step` picks a new step size and has a branch that returns nothing. Before reading the sessions, say what `h` must be at the failing line and how you could find who set it.

```python
def estimate_error(h):
    return 3.0 * h**3  # local error model c * h^3


def next_step(h, err, tol):
    if err > tol:
        return h / 2
    elif err < tol / 10:
        return h * 2


def run(h, tol, steps):
    t = 0.0
    for k in range(steps):
        err = estimate_error(h)
        h = next_step(h, err, tol)
        t += h
    return t


print(run(h=0.1, tol=1e-3, steps=6))
```
```text
Traceback (most recent call last):
  File "adapt.py", line 21, in <module>
    print(run(h=0.1, tol=1e-3, steps=6))
  File "adapt.py", line 17, in run
    t += h
TypeError: unsupported operand type(s) for +=: 'float' and 'NoneType'
```

One session does both jobs: post-mortem first, then a restart with a conditional breakpoint.

```console
$ python3 -m pdb adapt.py
> adapt.py(1)<module>()
-> def estimate_error(h):
(Pdb) c
[the program's traceback, as above]
Uncaught exception. Entering post mortem debugging
Running 'cont' or 'step' will restart the program
> adapt.py(17)run()
-> t += h
(Pdb) p h
None
(Pdb) p (err, tol, tol / 10)
(0.0003750000000000001, 0.001, 0.0001)
(Pdb) restart
Restarting adapt.py with arguments:

> adapt.py(1)<module>()
-> def estimate_error(h):
(Pdb) b 16, k == 1
Breakpoint 1 at adapt.py:16
(Pdb) c
> adapt.py(16)run()
-> h = next_step(h, err, tol)
(Pdb) p h, err
(0.05, 0.0003750000000000001)
(Pdb) s
--Call--
> adapt.py(5)next_step()
-> def next_step(h, err, tol):
(Pdb) n
> adapt.py(6)next_step()
-> if err > tol:
(Pdb) n
> adapt.py(8)next_step()
-> elif err < tol / 10:
(Pdb) n
--Return--
> adapt.py(8)next_step()->None
-> elif err < tol / 10:
(Pdb) w
  bdb.py(680)run()
-> exec(cmd, globals, locals)
  <string>(1)<module>()
  adapt.py(21)<module>()
-> print(run(h=0.1, tol=1e-3, steps=6))
  adapt.py(16)run()
-> h = next_step(h, err, tol)
> adapt.py(8)next_step()->None
-> elif err < tol / 10:
(Pdb) u
> adapt.py(16)run()
-> h = next_step(h, err, tol)
(Pdb) q
```

`c` ran to the exception and pdb stopped in `run`, line 17 (post-mortem). `p h` prints `None`: symptom confirmed, cause elsewhere. `p (err, tol, tol / 10)` shows $3.75\times10^{-4}$ between $10^{-4}$ and $10^{-3}$, the one case `next_step` never handles. `restart` reruns the program; `b 16, k == 1` stops at line 16 on the second pass only, and `p h, err` shows the arguments about to be passed. `s` enters `next_step` (`--Call--`); the first `n` reaches the `if`, the second runs it (false) and reaches the `elif`, the third runs that (false), falls off the end, and prints `--Return--` with `next_step()->None`. In `w` the first two entries (`bdb.py`, `<string>`) are pdb's own; yours are `<module>`, `run`, `next_step`, and `u` moves up to the caller.

■ **What just happened.** An `if`/`elif` with no `else` returns `None` on the missing branch, and the `TypeError` appeared a line later. The fix is `else: return h`. Five commands (`n`, `s`, `c`, `p`, `w`/`u`) replace twenty prints; the habit that matters is asking where a value *came from*.

### 3.9 Logging

`print` writes to stdout among your results, with no severity, filter or destination choice. `logging` separates **what happened** (the call) from **whether and where it appears** (configuration). Levels: `DEBUG` 10, `INFO` 20, `WARNING` 30 (the default threshold), `ERROR` 40, `CRITICAL` 50.

▶ **Build it yourself (5): logging in a simulation loop.** Explicit Euler for $x'' = -x$ should conserve $E = \frac12(x^2 + v^2)$ but does not. Log every step at `DEBUG`, a summary every 4 steps at `INFO`, and a `WARNING` once $E$ has grown by 10%. Predict the lines at each level.

```python
import logging
import sys

log = logging.getLogger("sim")


def simulate(h, steps):
    """Explicit Euler for x'' = -x. The exact energy 0.5(x^2 + v^2) is constant; Euler's grows."""
    x, v = 1.0, 0.0
    e0 = 0.5 * (x * x + v * v)
    warned = False
    log.info("start: h=%g steps=%d E0=%.4f", h, steps, e0)
    for k in range(1, steps + 1):
        x, v = x + h * v, v - h * x
        energy = 0.5 * (x * x + v * v)
        log.debug("step %d  x=%+.4f  v=%+.4f  E=%.4f", k, x, v, energy)
        if k % 4 == 0:
            log.info("step %d  E=%.4f  (drift %+.1f%%)", k, energy, 100 * (energy / e0 - 1))
        if energy > 1.1 * e0 and not warned:
            log.warning("energy grew by more than 10%% at step %d", k)
            warned = True
    return x, v


for level in (logging.INFO, logging.DEBUG):
    print(f"--- logging.basicConfig(level=logging.{logging.getLevelName(level)}) ---")
    logging.basicConfig(level=level, format="%(levelname)-7s %(name)s: %(message)s", stream=sys.stdout, force=True)
    simulate(h=0.2, steps=4)
```
```text
--- logging.basicConfig(level=logging.INFO) ---
INFO    sim: start: h=0.2 steps=4 E0=0.5000
WARNING sim: energy grew by more than 10% at step 3
INFO    sim: step 4  E=0.5849  (drift +17.0%)
--- logging.basicConfig(level=logging.DEBUG) ---
INFO    sim: start: h=0.2 steps=4 E0=0.5000
DEBUG   sim: step 1  x=+1.0000  v=-0.2000  E=0.5200
DEBUG   sim: step 2  x=+0.9600  v=-0.4000  E=0.5408
DEBUG   sim: step 3  x=+0.8800  v=-0.5920  E=0.5624
WARNING sim: energy grew by more than 10% at step 3
DEBUG   sim: step 4  x=+0.7616  v=-0.7680  E=0.5849
INFO    sim: step 4  E=0.5849  (drift +17.0%)
```

■ **What just happened.** Only the configuration changed between runs. The warning's cause is a one-line derivation: one Euler step is $x' = x + hv$, $v' = v - hx$, so

$$x'^2 + v'^2 = (x^2 + 2hxv + h^2v^2) + (v^2 - 2hxv + h^2x^2) = (1 + h^2)(x^2 + v^2). \tag{3.1}$$

The cross terms cancel, so $E_k = (1+h^2)^k E_0$ exactly. For $h = 0.2$ the factor is $1.04$: $1.04^3 = 1.125$ (the warning, $E = 0.5624$) and $1.04^4 - 1 = 17.0\%$ (the `INFO` drift). `project.py` asserts these.

**Handlers in one paragraph.** A *logger* (named, hierarchical by dots: `sim.euler` is a child of `sim`) creates records, which go to its handlers and then, by `propagate`, to its ancestors' handlers. A *handler* has its own level and *formatter*. So one logger can feed a console that shows warnings and a file that keeps everything: give the parent logger `sim` a `StreamHandler` with `setLevel(logging.WARNING)` and a `FileHandler` with `setLevel(logging.DEBUG)`; a record logged on `sim.euler` reaches both. `project.py` Part 6 does exactly this and writes `figures/run.log`.

**Rules.** A library does `log = logging.getLogger(__name__)` and never configures handlers (`integrators.py`; Example 4.6). The *application* calls `logging.basicConfig(...)` once, at its entry point (`force=True` replaces existing handlers). Pass arguments, not formatted strings: `log.debug("x=%s", x)` formats only if emitted. In an `except`, `log.exception("...")` logs at `ERROR` with the traceback. Do not log and re-raise; the catcher reports it.

### 3.10 Testing with pytest

**What pytest does.** It *discovers* tests by name, runs each as a function, and rewrites your `assert` statements so a failure prints the values on both sides (Example 4.3). Files: `test_*.py` or `*_test.py`; functions whose names *start with* `test`; methods of classes named `Test*` (no `__init__`). Exit status: 0 all passed, 1 failures, 2 interrupted, 3 internal error, 4 usage error, 5 nothing collected.

▶ **Try it first.** `test_basic.py` holds `test_a`, `check_b`, `testing_helper`, `TestC.test_c` and `Helper.test_d`. Which are collected?

```console
$ python3 -m pytest --collect-only -q
test_basic.py::test_a
test_basic.py::testing_helper
test_basic.py::TestC::test_c

3 tests collected in 0.00s
$ python3 -m pytest -q; echo "exit status $?"
no tests ran in 0.00s
exit status 5
```

`check_b` (wrong prefix) and `Helper` (wrong class prefix) are not; `testing_helper` *is*, since the prefix is `test`, not `test_`. A folder whose only file is `check_extra.py` yields "no tests ran", status 5. "Pytest doesn't see my tests" is almost always one of these rules. `python -m pytest` also puts the current directory on `sys.path`, so `import integrators` works here.

**Build it yourself (1): test-driven development of the trapezoid and Simpson rules.** **Red** (write a test, watch it fail for the right reason), **green** (least code that passes), **refactor** (clean up, tests as a net). A test that never failed has not shown it can.

*Step 0: decide what to assert, by deriving it.*

1. **Simpson is exact on cubics.** On $[-h, h]$ the rule is $\frac h3[f(-h) + 4f(0) + f(h)]$. For $f = 1$: $\frac h3\cdot 6 = 2h$, the integral. For $f = x^2$: $\frac h3\cdot 2h^2 = \frac23h^3$, the integral. For $f = x, x^3$ both sides are $0$ (odd). By linearity it is exact for every polynomial of degree $\le 3$, on any panel pair (shift and scale).
2. **Its error is $O(h^4)$.** For $f = x^4$ the rule gives $\frac h3\cdot 2h^4 = \frac23h^5$, the integral is $\frac25h^5$: error $\frac4{15}h^5$ per double panel. A smooth $f$ is a cubic plus $\frac{f^{(4)}}{24}x^4 + \cdots$, so the double-panel error is $\frac{f^{(4)}h^5}{90}$. With $n/2 = (b-a)/2h$ double panels,
$$S - I \approx \frac{(b-a)\,h^4}{180}\,f^{(4)}(\xi). \tag{3.2}$$
The same computation for the trapezoid rule (error $\frac{h^3}{6}\cdot\frac{f''}{2}$ per panel) gives $T - I \approx \frac{(b-a)h^2}{12}f''$. Doubling $n$ divides the trapezoid error by 4 and Simpson's by 16: *orders* 2 and 4. (This is a Taylor argument; the rigorous bound uses the Peano kernel, in any numerical analysis text.)
3. **Simpson is a Richardson combination of two trapezoids.** With $E$ the end samples halved, $O$ the odd-index and $P$ the even interior samples, $T_{2m} = h(E + O + P)$ and $T_m = 2h(E + P)$ on the same nodes. Then $4T_{2m} - T_m = h(2E + 4O + 2P)$, so
$$S_{2m} = \frac{4\,T_{2m} - T_m}{3} = \frac h3\big(y_0 + 4y_1 + 2y_2 + \cdots + 4y_{n-1} + y_n\big). \tag{3.3}$$
A second route to Simpson's weights, so a test can compare two implementations.

The trapezoid tests, written *before* `integrators.py` exists. The Simpson tests in the file follow the same pattern: four cubics (fact 1), equation (3.3) as a cross-check, the order (fact 2).

```python
def observed_order(
    rule: Rule, f: Fn, a: float, b: float, exact: float, ns: tuple[int, ...]
) -> float:
    """Slope of log(error) against log(h): the order of the rule, measured."""
    hs = [(b - a) / n for n in ns]
    errs = [abs(rule(f, a, b, n) - exact) for n in ns]
    return float(np.polyfit(np.log(hs), np.log(errs), 1)[0])
```
```python
# --- trapezoid: tests written first -------------------------------------------------------
LINES = [
    pytest.param(lambda x: 3.0, 0.0, 2.0, 6.0, id="constant"),
    pytest.param(lambda x: 2 * x + 1, 0.0, 2.0, 6.0, id="line"),
    pytest.param(lambda x: x, -1.0, 3.0, 4.0, id="x-on-[-1,3]"),
]


@pytest.mark.parametrize("n", [1, 2, 7])
@pytest.mark.parametrize("f, a, b, exact", LINES)
def test_trapezoid_exact_on_lines(f, a, b, exact, n):
    assert ig.trapezoid(f, a, b, n) == pytest.approx(exact, rel=1e-12, abs=1e-12)


def test_trapezoid_error_order():
    order = observed_order(ig.trapezoid, math.sin, 0.0, math.pi, 2.0, (16, 32, 64, 128))
    assert order == pytest.approx(2.0, abs=0.1)


@pytest.mark.parametrize(
    "n, error", [(0, ValueError), (-3, ValueError), (2.5, TypeError), ("4", TypeError)]
)
def test_trapezoid_rejects_bad_n(n, error):
    with pytest.raises(error):
        ig.trapezoid(math.sin, 0.0, 1.0, n)
```

*Red.* No module exists, so collection fails:

```console
$ python3 -m pytest -q --tb=line | grep -E "^E|error"
E   ModuleNotFoundError: No module named 'integrators'
ERROR test_integrators.py
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Interrupted: 1 error during collection !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1 error in 0.15s
```

*Green, naively.* The textbook formula:

```python
def trapezoid(f, a, b, n):
    h = (b - a) / n
    ys = [f(a + i * h) for i in range(n + 1)]
    return h * (0.5 * ys[0] + sum(ys[1:-1]) + 0.5 * ys[-1])
```
```console
$ python3 -m pytest -q --tb=no -rf
..........FF..                                                                               [100%]
===================================== short test summary info ======================================
FAILED test_integrators.py::test_trapezoid_rejects_bad_n[0-ValueError] - ZeroDivisionError: float...
FAILED test_integrators.py::test_trapezoid_rejects_bad_n[-3-ValueError] - IndexError: list index ...
2 failed, 12 passed in 0.08s
```

Twelve of 14 pass; the failures are informative: `n = 0` gave `ZeroDivisionError` and `n = -3` an `IndexError` from deep inside, where the contract says `ValueError` (`2.5` and `"4"` raised `TypeError` by luck). *Refactor*: validate in one helper that Simpson shares.

```python
def _panels(n: int, *, even: bool = False) -> int:
    n = operator.index(n)  # TypeError for 2.5 or "4"; accepts NumPy integers
    if n < 1:
        raise ValueError(f"need at least 1 panel, got n={n}")
    if even and n % 2:
        raise ValueError(f"Simpson's rule needs an even number of panels, got n={n}")
    return n


def trapezoid(f: Integrand, a: float, b: float, n: int) -> float:
    n = _panels(n)
    h = (b - a) / n
    ys = [f(a + i * h) for i in range(n + 1)]
    return h * (0.5 * ys[0] + sum(ys[1:-1]) + 0.5 * ys[-1])
```
```console
$ python3 -m pytest -q
..............                                                                               [100%]
14 passed in 0.08s
```

*Simpson, red then green.* With the Simpson tests added, 26 fail (no `simpson`); then:

```console
$ python3 -m pytest -q --tb=no | tail -1
26 failed, 14 passed in 0.15s
```
```python
def simpson(f: Integrand, a: float, b: float, n: int) -> float:
    n = _panels(n, even=True)
    h = (b - a) / n
    ys = [f(a + i * h) for i in range(n + 1)]
    stop = -2 if PLANT_OFF_BY_ONE else -1
    odd = sum(ys[1:stop:2])  # y_1, y_3, ..., y_{n-1}: weight 4
    even = sum(ys[2:-1:2])  # y_2, y_4, ..., y_{n-2}: weight 2
    return h / 3 * (ys[0] + 4 * odd + 2 * even + ys[-1])
```
```console
$ python3 -m pytest -q
........................................                                                     [100%]
40 passed in 0.10s
```

Two lines carry the story. `ys[1:stop:2]` collects the odd-index samples $y_1, y_3, \dots, y_{n-1}$ (weight 4); `stop = -1` ends the slice before $y_n$. `ys[2:-1:2]` collects the even interior samples (weight 2). **The planted bug**: `PLANT_OFF_BY_ONE = True` makes `stop = -2`, which drops the last odd sample $y_{n-1}$: the `ys[1:-2:2]` of §1. (The switch exists only for this module.)

*Plant it.* `pytest -p plant_bug` loads a module before the tests; this one sets the flag, so no source is edited. Below: the failures grouped by test function, then the test that names the missing sample. (`project.py` Part 3 runs this in-process and asserts exactly this set.)

```python
import integrators

integrators.PLANT_OFF_BY_ONE = True  # simpson now drops y_{n-1}
```
```console
$ python3 -m pytest -q --tb=no -rf -p plant_bug test_integrators.py | grep FAILED | sed 's/.*:://; s/\[.*//; s/ - .*//' | sort | uniq -c
      1 test_linearity_reversal_and_additivity
      1 test_simpson_error_order
     12 test_simpson_exact_on_cubics
      1 test_simpson_is_richardson_of_two_trapezoids
      6 test_simpson_known_integrals
      1 test_simpson_weights_are_1_4_2_4_1
$ python3 -m pytest -q --tb=short -p plant_bug "test_integrators.py::test_simpson_weights_are_1_4_2_4_1[3]"
F                                                                                            [100%]
============================================= FAILURES =============================================
______________________________ test_simpson_weights_are_1_4_2_4_1[3] _______________________________
test_integrators.py:120: in test_simpson_weights_are_1_4_2_4_1
    assert ig.simpson(indicator, 0.0, 1.0, 4) == pytest.approx(weights[k] * 0.25 / 3, rel=1e-14)
E   assert 0.0 == 0.3333333333333333 ± 1.0e-12
E
E     comparison failed
E     Obtained: 0.0
E     Expected: 0.3333333333333333 ± 1.0e-12
===================================== short test summary info ======================================
FAILED test_integrators.py::test_simpson_weights_are_1_4_2_4_1[3] - assert 0.0 == 0.3333333333333...
1 failed in 0.11s
```

■ **What just happened.** The bugged rule fails 22 of 53 test cases, in six test functions, and *no* trapezoid test. The plainest: the constant 1 on $[0, 1]$ with $n = 2$ integrates to $\frac13$, because the weight-4 sample was dropped. The most informative is `test_simpson_weights_are_1_4_2_4_1[3]`, the **regression test** (a test per bug fixed) written *after* the bug was understood: it feeds a function that is 1 at node $k$ and 0 elsewhere, so the answer *is* node $k$'s weight, and for $k = 3$ (index $n-1$) it returns $0.0$ instead of $\frac{4h}{3}$. And `test_linearity_reversal_and_additivity` found the bug with no known integral: two half-intervals of 4 panels must equal one of 8, but the bug drops a different sample in each. The tolerance trap: at $n = 100$ the buggy rule returns 1.998684, so `abs(value - 2) < 1e-2` *passes* it, while its relative error, $6.6\times10^{-4}$, fails a test at $10^{-7}$ (the correct rule's error there is about $5\times10^{-9}$).

**Tolerances are claims.** `pytest.approx(expected, rel=..., abs=...)` passes if the difference is within the *larger* tolerance (given only `rel`, an absolute floor of $10^{-12}$ still applies; given only `abs`, there is no relative part). For arrays, `np.testing.assert_allclose(actual, desired, rtol=..., atol=...)` passes if $|a - d| \le \text{atol} + \text{rtol}\,|d|$. A tolerance states how accurate your *algorithm* is: derive it from the method's error formula (Example 4.3) and the problem's conditioning (LA-13), and a failure then means a bug, not bad luck.

**`pytest.raises`** tests the error path, half the contract: `with pytest.raises(ValueError, match="even"):` passes only if that exception is raised with a message matching the regular expression; `as info` exposes `info.value` (see `test_adaptive_gives_up_with_a_convergence_error`). **`parametrize`** runs one body over a table, each row its own test (Example 4.4). **Fixtures** give tests set-up and tear-down by name: code before `yield` is set-up, after it tear-down. Built-ins: `tmp_path` (a fresh directory), `monkeypatch` (undoable changes), `capsys` (captured output), `caplog` (captured log records).

```python
@pytest.fixture
def table(tmp_path: Path) -> Iterator[Path]:
    path = tmp_path / "errors.csv"
    save_errors(path, [(4, 0.1), (8, 0.025)])
    print("\n  [setup]    wrote", path.name)
    yield path  # the test runs here, receiving `path`
    print("\n  [teardown] test finished")


def test_round_trip(table):
    assert load_errors(table) == [(4, 0.1), (8, 0.025)]
```
```console
$ python3 -m pytest -q -s
  [setup]    wrote errors.csv
.
  [teardown] test finished

1 passed in 0.01s
```

**Organization:** tests beside the code for a small project, a `tests/` directory for a package, shared fixtures in `conftest.py`. **Subsets:** `-k expr` selects by name (Example 4.4), `-x` stops at the first failure, `-q`/`-v` set verbosity, `--lf` re-runs last failures, `file.py::test_name[id]` runs one case.

**What to test.** *Behavior you promise*: boundaries (`n = 1`, `n = 2`), known closed forms ($\int_0^\pi\sin = 2$), **invariants** (linearity, swapping the limits flips the sign, splitting the interval, equation (3.3)), every **error path**, and a regression test per fixed bug. *Not implementation details* such as a private helper's name: those tests break when you improve the code and prove nothing.

**Property-based testing** is the idea behind the invariant test: draw many random inputs (25 from a seeded generator here) and check a property that must hold for all. The library `hypothesis` does this at scale and adds *shrinking* (it searches for the smallest failing input). Not required here; Your turn has you write the idea by hand.

### 3.11 Type hints

A **type hint** annotates a parameter or return value: `def simpson(f: Integrand, a: float, b: float, n: int) -> float`. **Python does not check them.** They are metadata (`__annotations__`) for readers and for tools such as mypy and pyright, which read your source and report inconsistent uses without running it. After `def double(x: int) -> int: return x * 2`, the call `double("ab")` returns `'abab'`, as `project.py` Part 7 asserts.

Vocabulary (3.11 spellings): `int`, `float`, `str`; `list[float]`, `dict[str, int]`, `tuple[int, ...]` (any length), `tuple[float, float]` (fixed); `X | None` (older: `Optional[X]`); `Callable[[float], float]` (from `collections.abc`); read-only inputs `Sequence[float]` and `Iterable[float]`, which accept lists, tuples and generators where `list[float]` accepts only lists. An `int` is accepted where a `float` is expected. A `TypeVar` means "the same type in and out": `def first(xs: Sequence[T]) -> T` with `T = TypeVar("T")`. `typing.Protocol` means "anything with these methods", no inheritance.

▶ **Build it yourself (6): hints on `Vec`, and one mypy error.** Annotate the PY-4 `Vec` (`from __future__ import annotations` lets a method name its own class; 3.14 makes that the default). Then write a script with two mistakes and predict which line mypy reports and which line Python reports.

```python
from __future__ import annotations

import math


class Vec:
    """An immutable vector in R^n: the PY-4 class, with type hints."""

    def __init__(self, *components: float) -> None:
        self._c: tuple[float, ...] = tuple(float(c) for c in components)

    def __repr__(self) -> str:
        return f"Vec({', '.join(map(repr, self._c))})"

    def __add__(self, other: Vec) -> Vec:
        if not isinstance(other, Vec):
            return NotImplemented
        return Vec(*(a + b for a, b in zip(self._c, other._c, strict=True)))

    def __abs__(self) -> float:
        return math.hypot(*self._c)

    def dot(self, other: Vec) -> float:
        return sum(a * b for a, b in zip(self._c, other._c, strict=True))
```
```python
from vec_typed import Vec

v = Vec(1, 2, 3)
w = v + Vec(4, 5, 6)  # fine
n: float = abs(w)  # fine
shifted = v + 1  # a mistake
d: Vec = v.dot(w)  # another mistake: dot returns a float
print(w, n)
print(shifted)
```
```console
$ python3 -m mypy use_vec.py
use_vec.py:6: error: Unsupported operand types for + ("Vec" and "int")  [operator]
use_vec.py:7: error: Incompatible types in assignment (expression has type "float", variable has type "Vec")  [assignment]
Found 2 errors in 1 file (checked 1 source file)
$ python3 use_vec.py
Traceback (most recent call last):
  File "use_vec.py", line 6, in <module>
    shifted = v + 1  # a mistake
TypeError: unsupported operand type(s) for +: 'Vec' and 'int'
```

Read the first error left to right: file and line; `error:`; `Unsupported operand types for + ("Vec" and "int")` (mypy looked up `Vec.__add__`, saw `other: Vec`, and `1` is an `int`); `[operator]`, the error code, which `# type: ignore[operator]` silences on one line (sparingly). The second error: a `float` assigned to a variable declared `Vec`. mypy reported *both* without running anything; Python reported only the first (a `TypeError`, via the `NotImplemented` of PY-4 §3.3) and would never reach the second.

**`Protocol` made real.** PY-4's `Evaluable` accepted `math.sin`; a checker is stricter than a run:

```python
import math
from typing import Protocol


class Evaluable(Protocol):  # PY-4 section 3.10, now checked
    def __call__(self, x: float, /) -> float: ...


def midpoint(f: Evaluable, a: float, b: float) -> float:
    return (b - a) * f(0.5 * (a + b))


print(midpoint(math.sin, 0.0, math.pi))
print(midpoint(math.atan2, 0.0, 1.0))  # atan2 needs two arguments
```
```console
$ python3 -m mypy proto.py
proto.py:14: error: Argument 1 to "midpoint" has incompatible type "Callable[[SupportsFloat | SupportsIndex, SupportsFloat | SupportsIndex], float]"; expected "Evaluable"  [arg-type]
proto.py:14: note: "Evaluable.__call__" has type "Callable[[float], float]"
Found 1 error in 1 file (checked 1 source file)
```

(My first draft omitted the `/`; mypy then rejected `math.sin`, whose argument the standard library declares *positional-only*, while the protocol promised `f(x=...)` works too. No runtime test would show that.) The error left is `math.atan2`, which needs two arguments, found before the call ran. `integrators.py` passes `mypy --strict` (`project.py` checks it when mypy is installed). Hints are **documentation that gets checked**, not enforcement.

### 3.12 Docstrings, doctests, formatting and linting

Use **NumPy-style docstrings**: summary, then `Parameters`, `Returns`, `Raises`, `Examples` (see `integrators.py`). **`doctest`** runs the `>>>` lines and compares output, so examples cannot rot (`python3 -m doctest integrators.py` is silent on success; `project.py` asserts it). `ruff format` (Black's style) rewrites code into one layout and `ruff check` flags unused imports, undefined names and bug-prone patterns; both files of this module pass them with `--line-length 100`. **Consistency beats preference:** with one layout a diff shows changes of meaning, not whitespace.

## 4. Worked examples

### Example 4.1: An exception climbing through three frames

The diagram of §2 as a program. `read_sample` fails; `mean_of` has a `finally`; `report` has nothing; the module level catches.

```python
import traceback


def read_sample(text):  # innermost frame
    return float(text)


def mean_of(texts):  # middle frame
    try:
        values = [read_sample(t) for t in texts]
    finally:
        print("  mean_of: finally block runs as the exception passes through")
    return sum(values) / len(values)


def report(texts):  # outer frame
    print(f"mean = {mean_of(texts):.3f}")


try:
    report(["1.5", "2.5", "x3"])
except ValueError as err:
    print("caught at the top:", err)
    for frame in traceback.extract_tb(err.__traceback__):
        print(f"  frame {frame.name:12s} line {frame.lineno}: {frame.line}")
```
```text
  mean_of: finally block runs as the exception passes through
caught at the top: could not convert string to float: 'x3'
  frame <module>     line 21: report(["1.5", "2.5", "x3"])
  frame report       line 17: print(f"mean = {mean_of(texts):.3f}")
  frame mean_of      line 10: values = [read_sample(t) for t in texts]
  frame read_sample  line 5: return float(text)
```

The `finally` message precedes "caught at the top", as in §2: cleanup runs as the stack unwinds, then the handler. `traceback.extract_tb` lists the four frames the exception crossed, outermost first, which is the order of a printed traceback (so you read it from the bottom, as in Example 4.5 and §3.8).

### Example 4.2: EAFP versus LBYL, a dictionary and a file

§3.5's program is the example; read its output as decisions. `counts.get` for a normal miss; `pick_rule("simpsons")` turns a typo's `KeyError` into a `ValueError` listing the choices; `read_lbyl` checked `exists()`, the file vanished in the gap, and `FileNotFoundError` came anyway, while `read_eafp` caught that one exception and returned `None`. `project.py` Part 4 times a dictionary hit and a miss under `try/except`: a miss costs several times a hit, so with many misses on a hot path prefer `.get`.

### Example 4.3: A failing pytest run with `approx`, then passing

Two tests every numerical programmer writes first. The expected values are right; the *tolerances* are wrong.

```console
$ python3 -m pytest -q --tb=short
FF                                                                                           [100%]
============================================= FAILURES =============================================
_______________________________________ test_plain_equality ________________________________________
test_tolerance.py:10: in test_plain_equality
    assert total == 0.3
E   assert 0.30000000000000004 == 0.3
_________________________________________ test_simpson_n8 __________________________________________
test_tolerance.py:14: in test_simpson_n8
    assert ig.simpson(math.sin, 0.0, math.pi, 8) == pytest.approx(2.0, rel=1e-6)
E   assert 2.0002691699483877 == 2.0 ± 2.0e-06
E
E     comparison failed
E     Obtained: 2.0002691699483877
E     Expected: 2.0 ± 2.0e-06
===================================== short test summary info ======================================
FAILED test_tolerance.py::test_plain_equality - assert 0.30000000000000004 == 0.3
FAILED test_tolerance.py::test_simpson_n8 - assert 2.0002691699483877 == 2.0 ± 2.0e-06
2 failed in 0.02s
$ python3 -m pytest -q
..                                                                                           [100%]
2 passed in 0.01s
```

(1) Plain `==` failed: the trapezoid rule on the constant $0.1$ with 7 panels sums to `0.30000000000000004`, and pytest printed both sides unasked. (2) `approx(2.0, rel=1e-6)` was a tolerance nobody derived. By (3.2) with $\max|f^{(4)}| = 1$ for $\sin$, Simpson at $n = 8$ is within $\frac{\pi(\pi/8)^4}{180} = 4.15\times10^{-4}$ absolute, $2.1\times10^{-4}$ relative; the measured error, $2.69\times10^{-4}$ ($1.3\times10^{-4}$ relative), sits inside, so `rel=2.2e-4` is a defensible claim and `1e-6` never was. After the two edits (the script's `sed`), both pass. `assert_allclose` prints similar diagnostics for arrays.

### Example 4.4: `parametrize` over six (function, exact integral) pairs

`CASES` in `test_integrators.py`: $\sin$ on $[0, \pi]$ (2); $e^x$ on $[0, 1]$ ($e - 1$); $1/(1+x^2)$ on $[0, 1]$ ($\pi/4$); $xe^{-x^2}$ on $[0, 1]$ ($(1-e^{-1})/2$); $\cos$ on $[0, \pi/2]$ (1); $\ln x$ on $[1, 2]$ ($2\ln 2 - 1$).

```console
$ python3 -m pytest -v -k known_integrals test_integrators.py
======================================= test session starts ========================================
collecting ... collected 53 items / 47 deselected / 6 selected

test_integrators.py::test_simpson_known_integrals[sin] PASSED                                [ 16%]
test_integrators.py::test_simpson_known_integrals[exp] PASSED                                [ 33%]
test_integrators.py::test_simpson_known_integrals[arctan] PASSED                             [ 50%]
test_integrators.py::test_simpson_known_integrals[gauss-x] PASSED                            [ 66%]
test_integrators.py::test_simpson_known_integrals[cos] PASSED                                [ 83%]
test_integrators.py::test_simpson_known_integrals[log] PASSED                                [100%]

================================= 6 passed, 47 deselected in 0.12s =================================
```

Each row is its own test with its own id; a failure names the row. `rel=1e-6` came from (3.2) at $n = 64$: the worst relative bound of the six is $5.1\times10^{-8}$ ($\sin$), so $10^{-6}$ is 20 times the worst case. `project.py` computes all six bounds with SymPy.

### Example 4.5: A chained traceback from `raise ... from`

`chain_demo.py` calls Newton from $x_0 = 0$ with nobody catching, so you see the cause chain.

```python
from rootfind import newton

newton(lambda x: x * x - 2, lambda x: 2 * x, 0.0)
```
```text
Traceback (most recent call last):
  File "rootfind.py", line 9, in newton
    step = f(x) / df(x)
ZeroDivisionError: float division by zero

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "chain_demo.py", line 3, in <module>
    newton(lambda x: x * x - 2, lambda x: 2 * x, 0.0)
  File "rootfind.py", line 11, in newton
    raise ConvergenceError("derivative is zero", k, abs(f(x))) from err
integrators.ConvergenceError: derivative is zero (iterations=1, residual=2)
```

Read a chain from the bottom: the final exception is what the caller sees (`ConvergenceError: derivative is zero (iterations=1, residual=2)`); the block above it, separated by "the direct cause of the following exception", is what broke inside. Without `from`, the separator reads "During handling of the above exception, another exception occurred", which says the handler itself failed: a different bug.

### Example 4.6: A logging transcript, INFO versus DEBUG

`adaptive_simpson` logs through `logging.getLogger("integrators")` and configures nothing; the application chooses what to see. `INFO` gives one line per call, `DEBUG` every split. Here $\int_0^\pi\sin = 2$ with a coarse tolerance, to keep the output short.

```python
import logging
import math
import sys

import integrators

for level in (logging.INFO, logging.DEBUG):
    print(f"--- logging.basicConfig(level=logging.{logging.getLevelName(level)}) ---")
    logging.basicConfig(level=level, format="%(levelname)-7s %(name)s: %(message)s", stream=sys.stdout, force=True)
    value = integrators.adaptive_simpson(math.sin, 0.0, math.pi, tol=1e-4)
    print(f"value = {value:.8f}")
```
```text
--- logging.basicConfig(level=logging.INFO) ---
INFO    integrators: adaptive_simpson on [0, 3.14159]: tol=0.0001, 17 evaluations, deepest level 2
value = 1.99999975
--- logging.basicConfig(level=logging.DEBUG) ---
DEBUG   integrators: depth 0: [0.00000, 3.14159] |delta|=8.98e-02 > 1.50e-03, split
DEBUG   integrators: depth 1: [0.00000, 1.57080] |delta|=2.15e-03 > 7.50e-04, split
DEBUG   integrators: depth 1: [1.57080, 3.14159] |delta|=2.15e-03 > 7.50e-04, split
INFO    integrators: adaptive_simpson on [0, 3.14159]: tol=0.0001, 17 evaluations, deepest level 2
value = 1.99999975
```

`INFO` tells you the cost: 17 evaluations, deepest level 2. `DEBUG` tells you where: $[0, \pi]$ failed its test at depth 0 ($|\delta| = 0.0898$ against $15\,\text{tol} = 1.5\times10^{-3}$), each half failed at depth 1 (against $15\,\text{tol}/2$), and the four quarters passed. Same library, two views, no edits. `project.py` writes a `DEBUG` run to `figures/run.log` and asserts line counts by level.

## 5. Why it matters and where it connects

- **PY-6 onward.** NumPy code is tested with `np.testing.assert_allclose` against a slow pure-Python version, as PY-4's `Matrix` was against `numpy.linalg.solve`; `np.errstate` (§3.6) makes NumPy fail loudly. **PY-10** builds packaging and continuous integration on §3.10's layout and exit codes.
- **LA-2, LA-13.** The LU project is checked by comparing $LU$ with $A$ and the solve with `numpy.linalg.solve`, `rtol` tied to the condition number: a tolerance is a claim about conditioning as much as about your code.
- **CM-1.** Energy conservation is an integrator's test oracle; the $(1+h^2)^k$ growth of (3.1) is Euler's failure on the oscillator (CM-1 §3.11 measures it and fixes it with symplectic Euler; CM-8 §3.11 explains why that works).
- **ML-2, LA-14.** A *gradient check* (an analytic gradient against PY-4 §3.7's dual-number derivative, with `approx`) tests every optimizer you write.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are predict-the-output, bug-fix and short-function drills; do all. B1–B7 are the module; B1 (test-driven midpoint rule), B2 (an exception hierarchy), B6 (type hints and invariance) and B7 (extend `project.py` with a mutation-testing table) are essential. C1–C3: a property-test runner with shrinking, the cost of exceptions measured, and hardening the integrators.

## 7. Coding project

`project.py` imports `integrators`, runs `test_integrators.py`, and asserts the module's claims: (1) the suite, through `pytest.main`, exits 0 (a fallback runner if pytest is missing); (2) convergence orders as log-log slopes, $2.0\pm0.1$ and $4.0\pm0.2$ (`figures/convergence.png`); (3) the planted bug is caught by exactly the six test functions of §3.10; (4) the exception hierarchy, `raise ... from`, EAFP timings; (5) context managers; (6) logging to `figures/run.log`, and the Euler numbers of §3.9; (7) type hints on every public function, mypy `--strict` if installed, the two errors of §3.11; (8) adaptive Simpson on $e^{-100(x-0.5)^2}$ against `scipy.integrate.quad` to $10^{-8}$, and the coarse-`tol` failure of §8.

**Build it yourself:** before opening the file, write the tests for a *new* rule (B1, midpoint) first and watch them fail. Then Your turn: a regression test for a bug you have hit; a property-based test for `sorted` with random lists; make the integrators raise on NaN; put a `breakpoint()` in the adaptive recursion and inspect the depth.

## 8. Where you'll struggle

- **Tests feel slow until the first one saves you.** Dozens of cases for a six-line function feel like overkill until a one-character slice error shows up as a red line instead of a wrong paper. §3.10 is designed to be that day.
- **`except Exception: pass`** turns a crash into a wrong answer (§3.3). If you must catch broadly, log the traceback. A `try` of more than a few lines is a smell: the handler cannot know which line failed. And **`assert` is not validation**: `python -O` removes it.
- **People avoid pdb for years and lose hours to print statements.** Learn five commands (`n`, `s`, `c`, `p`, `w`/`u`). The first session feels clumsy because you must decide what to look at; that is also why it is faster.
- **Float equality in tests needs `approx`, and the tolerance is a claim about your algorithm, not a fudge.** Loosening it until the test passes deletes the test (the lazy `1e-2` of §3.10 accepted the buggy rule). Work out the error bound; pick a tolerance a few times above it.
- **Type hints do nothing at runtime**, and students expect them to. The benefit comes from running mypy or pyright. Also `list[int]` is not accepted where `list[float]` is declared (lists are mutable, so the types are invariant); declare read-only inputs as `Sequence[float]` (B6).
- **"pytest can't find my tests"** is almost always a filename or function-name rule (§3.10); run `pytest --collect-only -q` first. **A test that never failed** may not be able to: break the code on purpose and watch it go red (B7).
- **An adaptive routine's `tol` is an estimate, not a guarantee.** `adaptive_simpson(peak, 0, 1, tol=1e-3)` on $e^{-100(x-0.5)^2}$ returns 0.152564 for a true value of 0.177245: the first samples sit on the flanks of a narrow peak, the two Simpson estimates agree, and the recursion stops. Check such routines against an independent reference (`project.py` Part 8; C3).

## 9. Check yourself

1. In which order do the `try`, `except`, `else` and `finally` bodies run for (a) no exception, (b) a matched exception, (c) an unmatched one? What happens to a climbing exception if `finally` executes `return`?
2. Why is `except Exception` usually wrong and a bare `except:` never right? Where is `except Exception` acceptable, and what must you do inside?
3. What does `raise B from A` record? What does the traceback say if you raise inside an `except` without `from`?
4. Give a case where EAFP is plainly right and one where LBYL is. Why is `if path.exists(): open(path)` wrong?
5. Why must `assert` not validate arguments?
6. In `@contextmanager`, why must cleanup sit in `try/finally` around the `yield`? What does a true return from `__exit__` do?
7. Name the pdb command that runs a line without entering calls; that enters a call; that runs until the function returns; that moves to the caller's frame.
8. List the five logging levels. Why does a library use `getLogger(__name__)` and never `basicConfig`?
9. Which names does pytest collect (file, function, class)? What is the exit status when it collects nothing?
10. Why is `assert integral == 0.3` a bad test, and why is `rel=1e-6` a "claim"? How would you choose a tolerance for Simpson's rule?
11. State the three facts of §3.10 Step 0, and the invariant that found the off-by-one with no known integral.
12. What does a hint do at run time, and what does mypy do with it? Why `Sequence[float]` for a read-only argument?

## 10. Pair with

- **Okken, *Python Testing with pytest*, 2nd ed., ch. 1–5.** The book-length version of §3.10 (test functions, fixtures, built-in fixtures, parametrization); better on fixture scope and on organizing a growing suite.
- **The pytest documentation ("Get started", "How to use fixtures", "How to parametrize fixtures and test functions").** The authority on every flag.
- **Python Tutorial §8, "Errors and Exceptions"; the language reference ("The try statement", "The with statement").** The exact statements behind §3.2 and §3.7.
- **Slatkin, *Effective Python*, the robustness, testing and typing items the roadmap lists** (numbering differs between editions). Short numbered rules; a good review checklist.
- **The `logging` HOWTO and the `pdb` documentation.** The best short treatment of loggers, handlers and propagation; every pdb command.
- **The mypy cheat sheet and the `typing` documentation.** Annotation spellings by example.
- **Software Carpentry, "Programming with Python" (errors, defensive programming, debugging).** A gentler pass over the same ideas.
- **Hunt and Thomas, *The Pragmatic Programmer*, the debugging chapter (reference only).** The "reproduce, minimize, don't assume" attitude of §3.8.
