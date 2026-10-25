# PY-5 Solutions

Every code block below was run, and every output was pasted by a script from that run (Python 3.13, pytest 9.1, mypy 2.4); tracebacks and pytest headers are shortened as in the notes. Files named in a solution sit next to `integrators.py` and `test_integrators.py`. Where a solution reports timings (C2, C3) the numbers belong to one machine and one run; the asserts check only their direction.

## A. Warm-up

**A1.**

```python
def f(items):
    out = []
    for it in items:
        try:
            out.append(10 // it)
        except ZeroDivisionError:
            out.append("zero")
            continue
        except TypeError:
            break
        else:
            out.append("ok")
        finally:
            out.append("f")
    return out


result = f([5, 0, "a", 2])
print(result)
assert result == [2, "ok", "f", "zero", "f", "f"]
```
```text
[2, 'ok', 'f', 'zero', 'f', 'f']
```

`5`: `10 // 5` is `2`, the `else` adds `"ok"`, and `finally` adds `"f"`. `0`: `ZeroDivisionError` is caught, `"zero"` is appended, and `continue` leaves the loop body, but `finally` still runs on the way out (rule 4 of §3.2), so another `"f"`. `"a"`: `10 // "a"` is a `TypeError`, the second `except` runs `break`, and `finally` runs again before the loop is left: the third `"f"`. The loop is over, so `2` is never processed. `else` was skipped for `0` because the `try` body did not finish normally.

**A2.**

```python
def load_tol(text: str | None, default: float = 1e-8) -> float:
    """Parse a tolerance; '' means 'use the default'; garbage is an error; None is a TypeError."""
    if text is None:
        raise TypeError("load_tol needs a string, got None")
    if text.strip() == "":
        return default
    try:
        return float(text)
    except ValueError as err:
        raise ValueError(f"cannot read a tolerance from {text!r}") from err


assert load_tol("1e-6") == 1e-6 and load_tol("") == 1e-8 and load_tol("  ") == 1e-8
for bad, exc in (("abc", ValueError), (None, TypeError)):
    try:
        load_tol(bad)
    except exc as err:
        print(f"{bad!r:6} -> {type(err).__name__}: {err}")
    else:
        raise AssertionError(bad)
```
```text
'abc'  -> ValueError: cannot read a tolerance from 'abc'
None   -> TypeError: load_tol needs a string, got None
```

The bare `except:` also caught `KeyboardInterrupt` and `SystemExit`: pressing Ctrl-C inside a loop that calls `load_tol` would have been swallowed and replaced by `1e-8`, and `load_tol(None)` returned the default instead of reporting the caller's bug. The repair separates the three cases: empty means "use the default" (decided by an `if`, not by an exception), unparsable text raises a `ValueError` that quotes the text and keeps the original with `from err`, and `None` raises a `TypeError` explicitly.

**A3.** The first assertion, `len(xs) > 0`, validates the *caller's* input: it must be `if ...: raise ValueError`. The second, `len(out) == len(xs)`, can only fail if this function's own loop is wrong: that is an invariant, and `assert` is right.

```python
def cumulative_mean(xs: list[float]) -> list[float]:
    if len(xs) == 0:  # input validation: the caller's mistake, so raise (survives -O)
        raise ValueError("cumulative_mean needs at least one value")
    total, out = 0.0, []
    for k, x in enumerate(xs, start=1):
        total += x
        out.append(total / k)
    assert len(out) == len(xs)  # internal invariant: only a bug in this function can break it
    return out


print(cumulative_mean([1.0, 2.0, 6.0]))
try:
    cumulative_mean([])
except ValueError as err:
    print("ValueError:", err)
```
```console
$ python3 a3.py
[1.0, 1.5, 3.0]
ValueError: cumulative_mean needs at least one value
$ python3 a3_O.py 2>&1 | tail -1
AssertionError: need data
$ python3 -O a3_O.py
[]
```

(`a3_O.py` is the exercise's original function followed by `print(cumulative_mean([]))`.) Under `python -O` both assertions vanish, so `cumulative_mean([])` returns `[]` instead of failing: no error, no data, a silent wrong answer for whoever divides by `len` next.

**A4.**

```python
def clamp(x: float, lo: float, hi: float) -> float:
    """Limit x to the interval [lo, hi]."""
    if lo > hi:
        raise ValueError(f"empty interval: lo={lo} > hi={hi}")
    return min(max(x, lo), hi)
```
```python
import pytest

from clamp import clamp


@pytest.mark.parametrize(
    "x, expected",
    [(0.5, 0.5), (-1.0, 0.0), (2.0, 1.0), (0.0, 0.0), (1.0, 1.0), (1e-300, 1e-300)],
    ids=["inside", "below", "above", "on-lo", "on-hi", "tiny"],
)
def test_clamp_to_unit_interval(x, expected):
    assert clamp(x, 0.0, 1.0) == expected


def test_clamp_rejects_empty_interval():
    with pytest.raises(ValueError, match="empty interval"):
        clamp(0.5, 1.0, 0.0)


def test_degenerate_interval_is_allowed():
    assert clamp(7.0, 2.0, 2.0) == 2.0
```
```console
$ python3 -m pytest -v
======================================= test session starts ========================================
collecting ... collected 8 items

test_clamp.py::test_clamp_to_unit_interval[inside] PASSED                                    [ 12%]
test_clamp.py::test_clamp_to_unit_interval[below] PASSED                                     [ 25%]
test_clamp.py::test_clamp_to_unit_interval[above] PASSED                                     [ 37%]
test_clamp.py::test_clamp_to_unit_interval[on-lo] PASSED                                     [ 50%]
test_clamp.py::test_clamp_to_unit_interval[on-hi] PASSED                                     [ 62%]
test_clamp.py::test_clamp_to_unit_interval[tiny] PASSED                                      [ 75%]
test_clamp.py::test_clamp_rejects_empty_interval PASSED                                      [ 87%]
test_clamp.py::test_degenerate_interval_is_allowed PASSED                                    [100%]

======================================== 8 passed in 0.03s =========================================
```

`min(max(x, lo), hi)` is correct when `lo <= hi`; the check comes first so an empty interval cannot return a number that looks plausible. The `parametrize` ids make a failure name its row; `match` is a regular expression searched in the message. The degenerate interval `lo == hi` is a boundary case worth its own test.

**A5.**

```console
$ python3 -m pytest --collect-only -q
b_test.py::test_u
test_a.py::test_x
test_a.py::TestZ::test_w

3 tests collected in 0.01s
$ python3 -m pytest -q --tb=no bench.py; echo "exit status $?"
F                                                                                            [100%]
===================================== short test summary info ======================================
FAILED bench.py::test_t - assert False
1 failed in 0.01s
exit status 1
```

Collected: `test_x`, `TestZ::test_w` (class name starts with `Test`), and `b_test.py::test_u` (the `*_test.py` file pattern). Not collected: `check_y` (the name does not start with `test`), `Helper.test_v` (the class name does not start with `Test`). `bench.py` does not match either file pattern, so bare `pytest` would never open it, but a file *named on the command line* is collected whatever its name: `test_t` ran, failed, and the exit status is 1.

**A6.**

```console
$ python3 -m mypy a6.py
a6.py:10: error: Argument 1 to "scale" has incompatible type "str"; expected "Sequence[float]"  [arg-type]
a6.py:11: error: Argument 2 to "scale" has incompatible type "str"; expected "float"  [arg-type]
Found 2 errors in 1 file (checked 1 source file)
```

`scale([1, 2, 3], 2)` is accepted: a `list[int]` is a `Sequence[float]` (`Sequence` is read-only, hence covariant, and `int` is accepted where `float` is expected). `scale((1.0, 2.0), 2.5)` is accepted. `scale("abc", 2)` is rejected: a `str` is a `Sequence[str]`, not of floats. `scale([1, 2], "2")` is rejected on the second argument. Python would have run both rejected calls without a word: `scale("abc", 2)` returns `['aa', 'bb', 'cc']` and `scale([1, 2], "2")` returns `['2', '22']` (sequence repetition), plausible-looking wrong answers that only the checker catches.

## B. Core

**B1.** *Derivation.* On one panel centred at $x_m$ with width $h$, write $f(x_m + t) = f_m + f'_m t + \tfrac12 f''_m t^2 + \cdots$ for $t \in [-h/2, h/2]$. The exact integral is $hf_m + \tfrac12 f''_m\cdot\tfrac{2}{3}(h/2)^3 = hf_m + \tfrac{1}{24}f''_m h^3$, the midpoint rule gives $hf_m$, so the midpoint error is $-\tfrac1{24}f''h^3$ per panel. The trapezoid rule on the same panel gives $\tfrac h2(f(x_m - h/2) + f(x_m + h/2)) = hf_m + \tfrac18f''_m h^3$, an error of $\tfrac18 - \tfrac1{24} = +\tfrac1{12}f''h^3$ (this is the $T - I \approx \frac{(b-a)h^2}{12}f''$ of §3.10). The ratio of errors is $\frac{1/12}{-1/24} = -2$: the trapezoid error is $-2$ times the midpoint error.

```python
import math

import numpy as np
import pytest

import integrators as ig
from midpoint_rule import midpoint

LINES = [
    pytest.param(lambda x: 3.0, 0.0, 2.0, 6.0, id="constant"),
    pytest.param(lambda x: 2 * x + 1, 0.0, 2.0, 6.0, id="line"),
    pytest.param(lambda x: x, -1.0, 3.0, 4.0, id="x-on-[-1,3]"),
]


@pytest.mark.parametrize("n", [1, 2, 7])
@pytest.mark.parametrize("f, a, b, exact", LINES)
def test_exact_on_lines(f, a, b, exact, n):
    assert midpoint(f, a, b, n) == pytest.approx(exact, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize(
    "f, a, b, exact",
    [(math.sin, 0.0, math.pi, 2.0), (math.exp, 0.0, 1.0, math.e - 1.0)],
    ids=["sin", "exp"],
)
def test_observed_order_is_two(f, a, b, exact):
    ns = (16, 32, 64, 128)
    hs = [(b - a) / n for n in ns]
    errs = [abs(midpoint(f, a, b, n) - exact) for n in ns]
    assert np.polyfit(np.log(hs), np.log(errs), 1)[0] == pytest.approx(2.0, abs=0.1)


@pytest.mark.parametrize("n", [32, 64, 128])
def test_error_is_minus_half_the_trapezoid_error(n):
    e_mid = midpoint(math.sin, 0.0, math.pi, n) - 2.0
    e_trap = ig.trapezoid(math.sin, 0.0, math.pi, n) - 2.0
    assert e_trap / e_mid == pytest.approx(-2.0, rel=1e-2)


@pytest.mark.parametrize("n, error", [(0, ValueError), (-1, ValueError), (1.5, TypeError)])
def test_bad_n(n, error):
    with pytest.raises(error):
        midpoint(math.sin, 0.0, 1.0, n)
```
```console
$ python3 -m pytest -q --tb=line | grep -E "^E"
E   ModuleNotFoundError: No module named 'midpoint_rule'
ERROR test_midpoint.py
$ python3 -m pytest -q
.................                                                                            [100%]
17 passed in 0.12s
$ python3 -m pytest -q --tb=no -rf -p plant_mid | grep FAILED | sed 's/.*:://; s/ - .*//; s/\[[^e].*//' | sort | uniq -c
      3 test_error_is_minus_half_the_trapezoid_error
      9 test_exact_on_lines
      1 test_observed_order_is_two[exp]
```

The first run is the red step (no module); then `midpoint_rule.py`:

```python
import operator

from integrators import Integrand


def midpoint(f: Integrand, a: float, b: float, n: int) -> float:
    """Composite midpoint rule: h * sum of f at the panel centres."""
    n = operator.index(n)
    if n < 1:
        raise ValueError(f"need at least 1 panel, got n={n}")
    h = (b - a) / n
    return h * sum(f(a + (i + 0.5) * h) for i in range(n))
```

All 17 pass. The planted bug is a loop that stops one panel early (`range(n - 1)`), loaded without editing any source by `plant_mid.py`, a `-p` plugin that replaces `midpoint_rule.midpoint`:

```python
import midpoint_rule

_original = midpoint_rule.midpoint


def buggy_midpoint(f, a, b, n):
    _original(f, a, b, n)  # keep the real validation
    h = (b - a) / n
    return h * sum(f(a + (i + 0.5) * h) for i in range(n - 1))  # off by one: last panel dropped


midpoint_rule.midpoint = buggy_midpoint
```

It is caught by the nine line tests, the three ratio tests, and `test_observed_order_is_two[exp]`, but **not** by `test_observed_order_is_two[sin]`: the integrand $\sin$ vanishes at $\pi$, so the dropped panel contributes only $h\sin(h/2)\approx h^2/2$, which still looks like order 2. That is why the order test has a second integrand ($e^x$ on $[0, 1]$, which does not vanish at an end). *A test integrand that is zero at the endpoints is blind to errors that live at the endpoints.*

**B2.**

```python
from collections.abc import Callable

from integrators import ConvergenceError, SolverError


class BracketError(SolverError, ValueError):
    """f(a) and f(b) have the same sign, so [a, b] need not contain a root."""

    def __init__(self, a: float, b: float, fa: float, fb: float) -> None:
        super().__init__(a, b, fa, fb)  # every constructor argument, so pickle and copy work
        self.a, self.b, self.fa, self.fb = a, b, fa, fb

    def __str__(self) -> str:
        return f"no sign change on [{self.a}, {self.b}]: f(a)={self.fa:.3g}, f(b)={self.fb:.3g}"


def _evaluate(f: Callable[[float], float], x: float) -> float:
    try:
        return f(x)
    except (ZeroDivisionError, OverflowError) as err:
        raise SolverError(f"f raised {type(err).__name__} at x={x!r}") from err


def bisect(f: Callable[[float], float], a: float, b: float, tol: float = 1e-12, max_iter: int = 100) -> float:
    """A root of f in [a, b] by bisection, assuming f(a) and f(b) differ in sign."""
    fa, fb = _evaluate(f, a), _evaluate(f, b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if (fa < 0) == (fb < 0):
        raise BracketError(a, b, fa, fb)
    for _ in range(max_iter):
        m = 0.5 * (a + b)
        fm = _evaluate(f, m)
        if fm == 0.0 or 0.5 * (b - a) < tol:
            return m
        if (fa < 0) != (fm < 0):
            b, fb = m, fm
        else:
            a, fa = m, fm
    raise ConvergenceError("bisection ran out of iterations", max_iter, b - a)
```
```python
import math
import pickle

import pytest

from integrators import ConvergenceError, SolverError
from roots import BracketError, bisect


def test_finds_sqrt_two():
    assert bisect(lambda x: x * x - 2, 1.0, 2.0) == pytest.approx(math.sqrt(2), abs=1e-11)


def test_same_sign_raises_bracket_error_catchable_two_ways():
    with pytest.raises(BracketError) as info:
        bisect(lambda x: x * x + 1, 0.0, 1.0)
    err = info.value
    assert isinstance(err, SolverError) and isinstance(err, ValueError)
    assert (err.a, err.b, err.fa, err.fb) == (0.0, 1.0, 1.0, 2.0)
    assert "no sign change" in str(err)
    assert pickle.loads(pickle.dumps(err)).fb == 2.0
    with pytest.raises(ValueError):  # the same exception, caught as the built-in
        bisect(lambda x: x * x + 1, 0.0, 1.0)


def test_out_of_iterations_is_a_convergence_error_not_a_value_error():
    with pytest.raises(ConvergenceError) as info:
        bisect(lambda x: x - 0.3, 0.0, 1.0, tol=1e-300, max_iter=5)
    err = info.value
    assert (err.iterations, err.residual) == (5, 1 / 32)
    assert not isinstance(err, ValueError) and err.__cause__ is None


def test_a_pole_is_reported_with_its_cause():
    with pytest.raises(SolverError) as info:  # 1/x changes sign at 0, but it is a pole, not a root
        bisect(lambda x: 1 / x, -1.0, 1.0)
    assert isinstance(info.value.__cause__, ZeroDivisionError)


def test_overflow_at_an_endpoint_is_converted():
    with pytest.raises(SolverError) as info:
        bisect(lambda x: math.exp(x) - 2, 0.0, 1000.0)
    assert isinstance(info.value.__cause__, OverflowError)
    assert "x=1000.0" in str(info.value)


def test_root_at_an_endpoint():
    assert bisect(lambda x: x, 0.0, 1.0) == 0.0
```
```console
$ python3 -m pytest -q
......                                                                                       [100%]
6 passed in 0.01s
```

`class BracketError(SolverError, ValueError)` is catchable both ways, and `test_same_sign_raises...` proves it by catching one exception as `SolverError` and then as `ValueError`. `ConvergenceError` is deliberately not a `ValueError` (the arguments were fine). `1/x` on $[-1, 1]$ changes sign at a *pole*, and the first midpoint is exactly 0, so `f` raises `ZeroDivisionError`; the caller sees a `SolverError` whose `__cause__` is that `ZeroDivisionError`. `math.exp(1000)` raises `OverflowError` at the endpoint and is converted the same way. Every constructor argument goes to `super().__init__`, so the pickle round-trip in the test works.

**B3.**

```python
import os
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from types import TracebackType


class cd:
    """Class-based: change the working directory for a with block, then restore it."""

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = Path(path)

    def __enter__(self) -> Path:
        self.old = Path.cwd()
        os.chdir(self.path)  # if this fails, __exit__ is never called and nothing changed
        return self.path

    def __exit__(
        self, exc_type: type[BaseException] | None, exc: BaseException | None, tb: TracebackType | None
    ) -> None:
        os.chdir(self.old)  # returns None, so a body exception propagates


@contextmanager
def cd_gen(path: str | os.PathLike[str]) -> Iterator[Path]:
    """Generator-based equivalent."""
    old = Path.cwd()
    os.chdir(path)
    try:
        yield Path(path)
    finally:  # runs when the body raises, too
        os.chdir(old)
```
```python
import contextlib
from pathlib import Path

import pytest

from cdtools import cd, cd_gen

MANAGERS = [cd, cd_gen, contextlib.chdir]


@pytest.fixture(autouse=True)
def start_in_tmp(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # every test starts in a known directory and cannot leak a chdir
    (tmp_path / "a" / "b").mkdir(parents=True)
    return tmp_path


@pytest.mark.parametrize("manager", MANAGERS, ids=lambda m: m.__name__)
def test_changes_and_restores(manager, start_in_tmp):
    with manager("a"):
        assert Path.cwd() == start_in_tmp / "a"
    assert Path.cwd() == start_in_tmp


@pytest.mark.parametrize("manager", MANAGERS, ids=lambda m: m.__name__)
def test_restores_when_the_body_raises(manager, start_in_tmp):
    with pytest.raises(RuntimeError, match="boom"):
        with manager("a"):
            raise RuntimeError("boom")
    assert Path.cwd() == start_in_tmp


@pytest.mark.parametrize("manager", MANAGERS, ids=lambda m: m.__name__)
def test_nesting(manager, start_in_tmp):
    with manager("a"):
        with manager("b"):
            assert Path.cwd() == start_in_tmp / "a" / "b"
        assert Path.cwd() == start_in_tmp / "a"
    assert Path.cwd() == start_in_tmp


@pytest.mark.parametrize("manager", MANAGERS, ids=lambda m: m.__name__)
def test_missing_directory_changes_nothing(manager, start_in_tmp):
    with pytest.raises(FileNotFoundError):
        with manager("nope"):
            pytest.fail("the body must not run")
    assert Path.cwd() == start_in_tmp
```
```console
$ python3 -m pytest -q
............                                                                                 [100%]
12 passed in 0.03s
```

The `autouse` fixture uses `monkeypatch.chdir` so each test starts in a known directory and, being undone by pytest, cannot leak a changed working directory into the next test. The class version relies on `__exit__` returning `None` (falsy), so a body exception propagates; if `os.chdir` fails in `__enter__`, `__exit__` is never called and nothing has changed (the last test). The generator version puts `os.chdir(path)` before the `try`, with the same effect. All three managers, including `contextlib.chdir`, pass the same tests.

**B4.**

```python
import logging
from collections.abc import Callable

from integrators import ConvergenceError

log = logging.getLogger(__name__)


def newton(
    f: Callable[[float], float],
    df: Callable[[float], float],
    x0: float,
    tol: float = 1e-12,
    max_iter: int = 50,
) -> float:
    """Newton's method for f(x) = 0, with a module logger."""
    x = x0
    for k in range(1, max_iter + 1):
        fx, dfx = f(x), df(x)
        log.debug("iteration %d: x=%.12g f=%.3e f'=%.3e", k, x, fx, dfx)
        if abs(dfx) < 1e-3:
            log.warning("small derivative |f'|=%.1e at x=%.6g (iteration %d)", abs(dfx), x, k)
        try:
            step = fx / dfx
        except ZeroDivisionError as err:
            raise ConvergenceError("derivative is zero", k, abs(fx)) from err
        x -= step
        if abs(step) < tol:
            log.info("converged to %.12g in %d iterations", x, k)
            return x
    raise ConvergenceError("no convergence", max_iter, abs(f(x)))
```
```python
import logging

import pytest

from integrators import ConvergenceError
from newton_log import newton

SQ = (lambda x: x * x - 2, lambda x: 2 * x)


def test_debug_per_iteration_then_one_info(caplog):
    calls = []
    f = lambda x: (calls.append(x), SQ[0](x))[1]  # counts iterations: f is called once per iteration
    with caplog.at_level(logging.DEBUG, logger="newton_log"):
        root = newton(f, SQ[1], 1.0)
    levels = [r.levelname for r in caplog.records]
    assert levels == ["DEBUG"] * len(calls) + ["INFO"]
    assert "converged to 1.41421356237" in caplog.records[-1].getMessage()
    assert root == pytest.approx(2**0.5, abs=1e-12)


def test_small_derivative_warns_once(caplog):
    with caplog.at_level(logging.DEBUG, logger="newton_log"):
        newton(*SQ, 1e-4)
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "iteration 1" in warnings[0].getMessage()


def test_quiet_by_default(caplog):
    newton(*SQ, 1.0)  # DEBUG and INFO are below the default WARNING threshold
    assert caplog.records == []


def test_failure_is_an_exception_not_an_error_record(caplog):
    with caplog.at_level(logging.DEBUG, logger="newton_log"):
        with pytest.raises(ConvergenceError):
            newton(lambda x: x**3 - 2 * x + 2, lambda x: 3 * x * x - 2, 0.0)
    assert all(r.levelno < logging.ERROR for r in caplog.records)
```
```console
$ python3 -m pytest -q
....                                                                                         [100%]
4 passed in 0.01s
```

`caplog.at_level(logging.DEBUG, logger="newton_log")` lowers that logger's threshold for the block. The count of `DEBUG` records equals the count of calls to `f` because `f` is evaluated once per iteration. Without `at_level`, the default `WARNING` threshold filters `DEBUG` and `INFO`, so a healthy run leaves `caplog.records` empty, which is the property that makes the logging free for library users. The failing run raises `ConvergenceError` and logs nothing at `ERROR`: the exception is the report.

**B5.** The session sets a conditional breakpoint on the `append` line for the second pass, steps to the end of the loop, and compares what the function produced with what it should have:

```console
$ python3 -m pdb ma.py
> ma.py(1)<module>()
-> def moving_average(xs, w):
(Pdb) b 4, i == 3
Breakpoint 1 at ma.py:4
(Pdb) c
> ma.py(4)moving_average()
-> out.append(sum(xs[i - w:i]) / w)
(Pdb) p i, xs[i - w:i]
(3, [2, 3])
(Pdb) n
> ma.py(3)moving_average()
-> for i in range(w, len(xs)):
(Pdb) n
> ma.py(5)moving_average()
-> return out
(Pdb) p len(out), len(xs) - w + 1
(2, 3)
(Pdb) p list(range(w, len(xs)))
[2, 3]
(Pdb) r
--Return--
> ma.py(5)moving_average()->[1.5, 2.5]
-> return out
(Pdb) q
```

`p i, xs[i - w:i]` prints `(3, [2, 3])`: on the second (and last) pass the window is `xs[1:3] = [2, 3]`, whose mean 2.5 is the last value returned. `p len(out), len(xs) - w + 1` prints `(2, 3)`: three windows should exist and two were produced. `p list(range(w, len(xs)))` prints `[2, 3]`: the loop never visits `i = 4`, which is the index that closes the last window `xs[2:4]` (mean 3.5). The slice `xs[i - w:i]` ends *before* `i`, so the loop index must run to `len(xs)` inclusive: `range(w, len(xs) + 1)`. `r` then shows the return value `[1.5, 2.5]` for the record. Fix and regression test (the first parametrized row is the bug):

```python
def moving_average(xs: list[float], w: int) -> list[float]:
    """Means of every window of w consecutive values: len(xs) - w + 1 of them."""
    if w < 1:
        raise ValueError(f"window must be at least 1, got {w}")
    return [sum(xs[i - w : i]) / w for i in range(w, len(xs) + 1)]  # was range(w, len(xs))
```
```python
import pytest

from ma_fixed import moving_average


@pytest.mark.parametrize(
    "xs, w, expected",
    [
        ([1, 2, 3, 4], 2, [1.5, 2.5, 3.5]),  # the regression: the last window was missing
        ([1, 2, 3, 4], 4, [2.5]),  # the window is the whole list
        ([1, 2, 3, 4], 1, [1, 2, 3, 4]),
        ([1, 2, 3], 5, []),  # window longer than the data
    ],
)
def test_window_count_and_values(xs, w, expected):
    assert moving_average(xs, w) == pytest.approx(expected)
    assert len(moving_average(xs, w)) == max(len(xs) - w + 1, 0)


def test_rejects_empty_window():
    with pytest.raises(ValueError, match="at least 1"):
        moving_average([1, 2, 3], 0)
```
```console
$ python3 -m pytest -q
.....                                                                                        [100%]
5 passed in 0.01s
```

**B6.**

```python
from collections.abc import Callable, Sequence


def mean(xs: Sequence[float]) -> float:
    return sum(xs) / len(xs)


def normalize(xs: Sequence[float]) -> list[float]:
    total = sum(xs)
    return [x / total for x in xs]


def first_positive(xs: Sequence[float]) -> float | None:
    for x in xs:
        if x > 0:
            return x
    return None


def apply_all(x: float, funcs: Sequence[Callable[[float], float]]) -> float:
    for fn in funcs:
        x = fn(x)
    return x
```
```python
from typed import first_positive, mean


def mean_list(xs: list[float]) -> float:  # the tempting annotation
    return sum(xs) / len(xs)


ints: list[int] = [1, 2, 3]
print(mean(ints))  # fine: Sequence is read-only, and int is accepted where float is expected
print(mean_list(ints))  # rejected: list is invariant

p = first_positive([-1.0, 2.0])
print(p + 1)  # rejected: p may be None
if p is not None:  # the narrowing that fixes it
    print(p + 1)
```
```console
$ python3 -m mypy --strict typed.py
Success: no issues found in 1 source file
$ python3 -m mypy --strict misuse.py
misuse.py:10: error: Argument 1 to "mean_list" has incompatible type "list[int]"; expected "list[float]"  [arg-type]
misuse.py:10: note: "list" is invariant -- see https://mypy.readthedocs.io/en/stable/common_issues.html#variance
misuse.py:10: note: Consider using "Sequence" instead, which is covariant
misuse.py:13: error: Unsupported operand types for + ("None" and "int")  [operator]
misuse.py:13: note: Left operand is of type "float | None"
Found 2 errors in 1 file (checked 1 source file)
$ python3 -m mypy --strict fixed.py
Success: no issues found in 1 source file
$ python3 fixed.py
2.0
3.0
```

(a) `list[int]` is not a `list[float]`: if it were, a function holding the list as `list[float]` could append `2.5` to it and the caller's `list[int]` would hold a float. Lists are mutable, so the type is *invariant*. A read-only parameter should say so with `Sequence[float]`, which is covariant. mypy's note even suggests this. (b) `first_positive` returns `float | None`, and `None + 1` is an error; the `if p is not None:` test *narrows* the type inside the block, which is the fix shown in `fixed.py`.

**B7.** Run from the module folder, or append the body (from `ORIGINAL = ...` on) to `project.py` as a new part before "Your turn"; the code uses only names `project.py` already imports, and was tested both ways.

```python
"""Mutation testing by hand: run the suite against deliberately broken copies of simpson."""
import contextlib
import io
import math
import sys
from collections.abc import Callable

import pytest

import integrators as ig

ORIGINAL = ig.simpson
Rule = Callable[[Callable[[float], float], float, float, int], float]


def samples(f: Callable[[float], float], a: float, b: float, n: int, validate: bool = True):
    if validate:
        ORIGINAL(f, a, b, n)  # the real function, run only to raise its validation errors
    h = (b - a) / n
    return h, [f(a + i * h) for i in range(n + 1)]


def swapped_weights(f, a, b, n):
    h, ys = samples(f, a, b, n)
    return h / 3 * (ys[0] + 2 * sum(ys[1:-1:2]) + 4 * sum(ys[2:-1:2]) + ys[-1])


def wrong_step(f, a, b, n):
    ORIGINAL(f, a, b, n)
    h = (b - a) / (n + 1)
    ys = [f(a + i * h) for i in range(n + 1)]
    return h / 3 * (ys[0] + 4 * sum(ys[1:-1:2]) + 2 * sum(ys[2:-1:2]) + ys[-1])


def first_sample_dropped(f, a, b, n):
    h, ys = samples(f, a, b, n)
    return h / 3 * (4 * sum(ys[1:-1:2]) + 2 * sum(ys[2:-1:2]) + ys[-1])


def no_even_check(f, a, b, n):
    h, ys = samples(f, a, b, n, validate=False)
    return h / 3 * (ys[0] + 4 * sum(ys[1:-1:2]) + 2 * sum(ys[2:-1:2]) + ys[-1])


def stop_minus_two(f, a, b, n):  # ys[2:-2:2] instead of ys[2:-1:2]
    h, ys = samples(f, a, b, n)
    return h / 3 * (ys[0] + 4 * sum(ys[1:-1:2]) + 2 * sum(ys[2:-2:2]) + ys[-1])


MUTANTS: dict[str, Rule] = {
    "swapped weights 2/4": swapped_weights,
    "h = (b-a)/(n+1)": wrong_step,
    "first sample dropped": first_sample_dropped,
    "even-n check removed": no_even_check,
    "ys[2:-2:2] for ys[2:-1:2]": stop_minus_two,
}


class Collect:
    def __init__(self) -> None:
        self.failed: set[str] = set()

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call" and report.failed:
            self.failed.add(report.nodeid.split("::")[1].split("[")[0])


def kills(rule: Rule) -> set[str]:
    sys.modules.pop("test_integrators", None)
    collector = Collect()
    ig.simpson = rule
    try:
        with contextlib.redirect_stdout(io.StringIO()):  # keep pytest's own report out of our table
            pytest.main(["-q", "--tb=no", "-rN", "-p", "no:cacheprovider", "test_integrators.py"], plugins=[collector])
    finally:
        ig.simpson = ORIGINAL
    return collector.failed


print(f"{'mutant':28s} killed by")
survivors = []
for name, rule in MUTANTS.items():
    failed = kills(rule)
    print(f"{name:28s} {', '.join(sorted(t.removeprefix('test_') for t in failed)) or 'NOBODY (survives)'}")
    if not failed:
        survivors.append(name)

# The only survivor is an equivalent mutant: for even n the two slices select the same samples.
assert survivors == ["ys[2:-2:2] for ys[2:-1:2]"]
for n in (2, 4, 10, 48):
    assert stop_minus_two(math.sin, 0.0, math.pi, n) == ORIGINAL(math.sin, 0.0, math.pi, n)
print("survivor is equivalent to the original for every even n, so no test can kill it")
```
```console
$ python3 mutants.py
mutant                       killed by
swapped weights 2/4          linearity_reversal_and_additivity, simpson_error_order, simpson_exact_on_cubics, simpson_is_richardson_of_two_trapezoids, simpson_known_integrals, simpson_weights_are_1_4_2_4_1
h = (b-a)/(n+1)              linearity_reversal_and_additivity, simpson_error_order, simpson_exact_on_cubics, simpson_is_richardson_of_two_trapezoids, simpson_known_integrals, simpson_weights_are_1_4_2_4_1
first sample dropped         linearity_reversal_and_additivity, simpson_exact_on_cubics, simpson_is_richardson_of_two_trapezoids, simpson_known_integrals, simpson_weights_are_1_4_2_4_1
even-n check removed         simpson_rejects_bad_n, simpson_rejects_odd_n
ys[2:-2:2] for ys[2:-1:2]    NOBODY (survives)
survivor is equivalent to the original for every even n, so no test can kill it
```

Four mutants are killed, each by the tests you would expect: the weight and step mutants by the exactness, order and known-integral tests; the dropped sample by exactness and the weights test; the removed `n` check only by the two `rejects` tests (the numerical tests all use even `n`). The fifth is the **equivalent mutant** that was added to show the limit of the method: `ys[2:-2:2]` selects the same samples as `ys[2:-1:2]` whenever `n` is even (the last even interior index is $n-2$ either way), so *no test can kill it* and it is not a gap in the suite. A survivor needs a look, not automatic panic: write the test if the mutant changes behavior, and recognize it as equivalent if it does not. (Compare the real bug of the notes, `ys[1:-2:2]`: there the odd samples *do* differ, because the last odd index $n-1$ is dropped.)

## C. Challenge

**C1.**

```python
from collections import Counter
from collections.abc import Callable, Iterator

import numpy as np

IntList = list[int]


def gen_list(rng: np.random.Generator) -> IntList:
    return [int(v) for v in rng.integers(-20, 21, size=int(rng.integers(0, 12)))]


def candidates(xs: IntList) -> Iterator[IntList]:
    """Strictly smaller neighbours: delete one element, or move one number toward 0."""
    for i in range(len(xs)):
        yield xs[:i] + xs[i + 1 :]
    for i, v in enumerate(xs):
        for w in (0, v // 2, v - (1 if v > 0 else -1)):
            if abs(w) < abs(v):
                yield xs[:i] + [w] + xs[i + 1 :]


def forall(gen: Callable[[np.random.Generator], IntList], prop: Callable[[IntList], bool],
           n: int = 200, seed: int = 0) -> dict[str, object] | None:
    """Check prop on n random inputs; on the first failure, shrink it to a local minimum."""
    rng = np.random.default_rng(seed)
    for trial in range(n):
        x = gen(rng)
        if prop(x):
            continue
        original, progress = x, True
        while progress:  # greedy: take the first smaller candidate that still fails
            progress = False
            for cand in candidates(x):
                if not prop(cand):
                    x, progress = cand, True
                    break
        return {"trial": trial, "original": original, "shrunk": x}
    return None


def sorter_properties(sorter: Callable[[IntList], IntList]) -> Callable[[IntList], bool]:
    def prop(xs: IntList) -> bool:
        out = sorter(list(xs))
        ordered = all(a <= b for a, b in zip(out, out[1:]))
        return ordered and Counter(out) == Counter(xs) and sorter(list(out)) == out
    return prop


def bad_insertion_sort(xs: IntList) -> IntList:
    xs = list(xs)
    for i in range(1, len(xs) - 1):  # bug: the last element is never inserted
        j = i
        while j > 0 and xs[j - 1] > xs[j]:
            xs[j - 1], xs[j] = xs[j], xs[j - 1]
            j -= 1
    return xs


print("sorted:            ", forall(gen_list, sorter_properties(sorted)))
result = forall(gen_list, sorter_properties(bad_insertion_sort))
print("bad_insertion_sort:", result)
assert forall(gen_list, sorter_properties(sorted)) is None
assert result is not None and result["shrunk"] == [1, 0]
assert len(result["original"]) > 2  # the raw counterexample was bigger
```
```console
$ python3 forall.py
sorted:             None
bad_insertion_sort: {'trial': 0, 'original': [6, 0, -9, -8, -19, -17, -20, -13, 13, 6], 'shrunk': [1, 0]}
```

`sorted` passes 200 random lists (`None` means no counterexample). The buggy sort fails on the very first draw, a ten-element list; greedy shrinking, which deletes an element or moves a number toward 0 whenever the property still fails, reduces it to `[1, 0]`, the smallest unsorted list over small integers. That is also the diagnosis: with two elements the loop `range(1, len(xs) - 1)` is empty, so the last element is never inserted. This is what `hypothesis` automates (with smarter strategies and a database of past failures). Shrinking is *greedy*, so it finds a local minimum, not necessarily the global one.

**C2.**

```python
import platform
import time
from collections.abc import Callable

import numpy as np

N, SIZE = 200_000, 1000
table = {i: i for i in range(SIZE)}


def keys_with_miss_rate(rate: float) -> list[int]:
    rng = np.random.default_rng(0)
    hit = rng.integers(0, SIZE, N)
    miss = -1 - rng.integers(0, SIZE, N)  # negative keys are never in the table
    return np.where(rng.random(N) < rate, miss, hit).tolist()


def with_try(keys: list[int]) -> int:
    total = 0
    for k in keys:
        try:
            total += table[k]
        except KeyError:
            pass
    return total


def with_in(keys: list[int]) -> int:
    total = 0
    for k in keys:
        if k in table:
            total += table[k]
    return total


def with_get(keys: list[int]) -> int:
    total = 0
    for k in keys:
        v = table.get(k)
        if v is not None:
            total += v
    return total


def micros(fn: Callable[[list[int]], int], keys: list[int]) -> float:
    best = float("inf")
    for _ in range(5):
        t0 = time.perf_counter()
        fn(keys)
        best = min(best, time.perf_counter() - t0)
    return best / len(keys) * 1e6


METHODS = {"try/except": with_try, "if k in d": with_in, "d.get(k)": with_get}
grid = [k / 10 for k in range(10)]  # miss rates 0%, 10%, ..., 90%
timing = {name: [micros(fn, keys_with_miss_rate(rate)) for rate in grid] for name, fn in METHODS.items()}

print(f"Python {platform.python_version()}, {platform.machine()}; microseconds per lookup")
print(f"{'miss rate':>9} " + " ".join(f"{name:>11}" for name in METHODS))
for rate in (0.0, 0.1, 0.5, 0.9):
    print(f"{rate:9.0%} " + " ".join(f"{timing[name][int(rate * 10)]:11.3f}" for name in METHODS))

# Straight-line fits t(rate) = a + b * rate smooth out timing noise; try/except stops being the fastest
# where its line crosses the lower of the other two.
fits = {name: np.polyfit(grid, ts, 1) for name, ts in timing.items()}  # (slope, intercept) per method
slope = {name: fit[0] for name, fit in fits.items()}
intercept = {name: fit[1] for name, fit in fits.items()}
crossings = [(intercept[o] - intercept["try/except"]) / (slope["try/except"] - slope[o])
             for o in ("if k in d", "d.get(k)") if slope["try/except"] > slope[o]]
crossover = max(0.0, min(crossings))
print(f"try/except stops being the fastest at a miss rate of about {crossover:.0%} (straight-line fits)")

at = lambda name, rate: timing[name][int(rate * 10)]  # noqa: E731
assert at("try/except", 0.9) == max(at(n, 0.9) for n in METHODS)  # at 90% misses it is the slowest
assert at("try/except", 0.9) > 1.3 * at("try/except", 0.1)  # its cost grows with the miss rate
assert all(at("try/except", 0.9) > 1.5 * at(n, 0.9) for n in ("if k in d", "d.get(k)"))
assert 0.0 <= crossover <= 0.5
```
```console
$ python3 exc_cost.py
Python 3.13.16, x86_64; microseconds per lookup
miss rate  try/except   if k in d    d.get(k)
       0%       0.055       0.085       0.070
      10%       0.067       0.080       0.066
      50%       0.113       0.060       0.059
      90%       0.128       0.033       0.040
try/except stops being the fastest at a miss rate of about 8% (straight-line fits)
```

In this run `try/except` is the fastest while misses are rare (a hit costs one lookup, with no function call and no second `in` test), its cost climbs steadily with the miss rate, and the last line gives the miss rate at which it stops being the fastest. `if k in d` and `d.get(k)` hardly depend on the miss rate (misses are even a little cheaper, since they skip the addition). *What it says:* raising and catching an exception costs several lookups' worth of time, so when misses are common, ask with `in` or `.get`. *What it does not say:* that exceptions are slow in general. The cost matters only in a hot loop with a high failure rate; a 2-microsecond handler on a rare error is invisible. EAFP is chosen for correctness (no gap between check and use, §3.5) and clarity before speed, and the right question is the one in §3.5: is a miss normal (use `.get`) or a mistake (use `try`)? Measure before optimizing: the crossover moves with the Python version and the machine, and the numbers above come from one run.

**C3.**

```python
import math
from collections.abc import Callable

import integrators as ig

Fn = Callable[[float], float]


class NonFiniteError(ig.SolverError, ValueError):
    """A limit, a sample, or a result was NaN or infinite."""


def _finite(value: float, what: str) -> None:
    if not math.isfinite(value):
        raise NonFiniteError(f"{what} is not finite: {value!r}")


def checked(f: Fn) -> Fn:
    """Wrap f so that the first non-finite sample raises, naming the x where it happened."""

    def wrapper(x: float) -> float:
        y = f(x)
        if not math.isfinite(y):
            raise NonFiniteError(f"f({x!r}) is not finite: {y!r}")
        return y

    return wrapper


def trapezoid(f: Fn, a: float, b: float, n: int) -> float:
    """O(1) overhead: a non-finite sample makes the sum non-finite, so test only the result."""
    _finite(a, "a"), _finite(b, "b")
    result = ig.trapezoid(f, a, b, n)
    _finite(result, "the integral (a sample was NaN/inf, or the sum overflowed)")
    return result


def simpson(f: Fn, a: float, b: float, n: int) -> float:
    _finite(a, "a"), _finite(b, "b")
    result = ig.simpson(f, a, b, n)
    _finite(result, "the integral (a sample was NaN/inf, or the sum overflowed)")
    return result


def adaptive_simpson(f: Fn, a: float, b: float, tol: float = 1e-10, max_depth: int = 40, min_depth: int = 0) -> float:
    """adaptive_simpson that fails fast on non-finite samples and bisects min_depth times first."""
    _finite(a, "a"), _finite(b, "b")
    if min_depth < 0:
        raise ValueError(f"min_depth must be >= 0, got {min_depth}")
    parts = 2**min_depth
    edges = [a + (b - a) * i / parts for i in range(parts + 1)]
    g = checked(f)
    return sum(
        ig.adaptive_simpson(g, lo, hi, tol / parts, max_depth) for lo, hi in zip(edges, edges[1:])
    )
```
```python
import math

import pytest

import harden
import integrators as ig

RULES = [harden.trapezoid, harden.simpson]


@pytest.mark.parametrize("rule", RULES, ids=["trapezoid", "simpson"])
@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_non_finite_sample_raises(rule, bad):
    def f(x: float) -> float:
        return bad if x == 0.5 else x

    with pytest.raises(harden.NonFiniteError, match="not finite"):
        rule(f, 0.0, 1.0, 4)  # x = 0.5 is the middle node


@pytest.mark.parametrize("rule", RULES, ids=["trapezoid", "simpson"])
@pytest.mark.parametrize("a, b", [(math.nan, 1.0), (0.0, math.inf)])
def test_non_finite_limits_raise(rule, a, b):
    with pytest.raises(ValueError):  # NonFiniteError is also a ValueError
        rule(math.sin, a, b, 4)


@pytest.mark.parametrize("rule, base", [(harden.trapezoid, ig.trapezoid), (harden.simpson, ig.simpson)])
def test_finite_results_are_bit_identical(rule, base):
    assert rule(math.exp, 0.0, 1.0, 64) == base(math.exp, 0.0, 1.0, 64)


def test_adaptive_names_the_bad_point_instead_of_running_out_of_depth():
    with pytest.raises(harden.NonFiniteError, match=r"f\(0\.5\)"):
        harden.adaptive_simpson(lambda x: math.nan if x == 0.5 else x, 0.0, 1.0)


def test_error_hierarchy():
    assert issubclass(harden.NonFiniteError, ig.SolverError) and issubclass(harden.NonFiniteError, ValueError)
    assert not issubclass(harden.NonFiniteError, ig.ConvergenceError)
```
```console
$ python3 -m pytest -q
..............                                                                               [100%]
14 passed in 0.02s
$ python3 measure.py
trapezoid checked / plain time = 0.985
simpson   checked / plain time = 1.002
min_depth  |error|      requested tol = 1e-3
        0  2.47e-02
        1  2.47e-02
        2  2.68e-05
        3  1.42e-05
        4  8.08e-14
```

(a) A NaN or an infinity in any sample makes the final sum non-finite, so `trapezoid` and `simpson` need only an $O(1)$ test on the result (plus the limits): the measured cost ratio is about 1.00. The price is a worse message (the error cannot say *where* the bad sample was); the adaptive routine evaluates few points, so it wraps `f` with `checked`, which names the offending $x$ and replaces the misleading `ConvergenceError` of §8. `NonFiniteError(SolverError, ValueError)` is catchable as both, like `BracketError` in B2. The tests are written first, check all three non-finite values, and include a bit-for-bit comparison that the finite path is unchanged.

(b) The failure at `tol=1e-3` is a failure of the *first estimate*, not of the recursion: the algorithm compares one Simpson estimate on the whole interval with two on its halves, and when all five nodes sit on the flat flanks of the peak or at an end of it, the two agree and the recursion stops. With `min_depth=1` the peak sits at the shared endpoint $x = 0.5$ of both halves, and each half's three nodes still miss its shape: the error is unchanged. With `min_depth=2`, the intervals $[0.25, 0.5]$ and $[0.5, 0.75]$ have their midpoints at $0.375$ and $0.625$, on the peak's flank, the two estimates now disagree, and the recursion refines; the error drops by three orders of magnitude, and at `min_depth=4` to rounding level. The moral is the one of §8: `tol` is an estimate, so give the routine a first look fine enough to see the features.
