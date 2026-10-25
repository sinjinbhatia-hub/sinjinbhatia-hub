# %% [markdown]
# # PY-5 Project: trust, but verify (integrators, tests, exceptions, logging, types)
#
# Runs the test suite for `integrators.py` (trapezoid, Simpson, adaptive Simpson), then checks every
# claim the notes make: the planted off-by-one is caught by exactly six test functions (§3.10); the
# convergence orders are 2 and 4; the exception hierarchy, `raise ... from`, context managers,
# logging, type hints and `assert` under `-O` behave as described; adaptive Simpson matches
# `scipy.integrate.quad`, and its `tol` is only an estimate (§8). Numbers quoted in the notes are
# asserted here. Run headless: `MPLBACKEND=Agg python3 project.py` (about 10 s). pytest is needed
# for Parts 1 and 3 (a tiny fallback runs without it); mypy and ruff are used if installed.

# %%
import doctest
import importlib.util
import inspect
import logging
import math
import pickle
import shutil
import subprocess
import sys
import tempfile
import time
import timeit
import typing
import warnings
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import integrators as ig  # noqa: E402

FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
T0 = time.perf_counter()
A, B, EXACT = 0.0, math.pi, 2.0
SIX = {"test_simpson_exact_on_cubics", "test_simpson_is_richardson_of_two_trapezoids",
       "test_simpson_error_order", "test_simpson_known_integrals", "test_simpson_weights_are_1_4_2_4_1",
       "test_linearity_reversal_and_additivity"}


E = typing.TypeVar("E", bound=BaseException)


def raises(exc: type[E], fn: Callable[..., object], *args: object) -> E:
    """Call fn(*args), require it to raise exc (a hand-made `pytest.raises`), return the exception."""
    try:
        fn(*args)
    except exc as err:
        return err
    raise AssertionError(f"{fn.__name__}{args} did not raise {exc.__name__}")


# %% [markdown]
# ## Part 1: run the test suite (notes §3.10)
# `pytest.main` in-process, with a small plugin that records which test ids fail. Without pytest a
# handful of direct checks stands in (the real suite needs pytest).

# %%
HAVE_PYTEST = importlib.util.find_spec("pytest") is not None


def run_suite() -> tuple[int, list[str], list[str]]:
    """Run test_integrators.py in-process; return (exit code, failed node ids, passed node ids)."""
    import pytest

    class Collect:
        failed: list[str] = []
        passed: list[str] = []

        def pytest_runtest_logreport(self, report: "pytest.TestReport") -> None:
            if report.when == "call":
                (Collect.failed if report.failed else Collect.passed).append(report.nodeid)

    Collect.failed, Collect.passed = [], []
    sys.modules.pop("test_integrators", None)  # import it afresh on every run
    args = ["-q", "--tb=no", "-rN", "-p", "no:cacheprovider", str(HERE / "test_integrators.py")]
    code = pytest.main(args, plugins=[Collect()])
    return int(code), Collect.failed, Collect.passed


if HAVE_PYTEST:
    code, failed, passed_ids = run_suite()
    N_TESTS = len(passed_ids)
    assert code == 0 and not failed, f"suite failed: {failed}"
    print(f"Part 1: suite passed ({N_TESTS} test cases)")
    group = lambda *prefixes: sum(i.split("::")[1].startswith(prefixes) for i in passed_ids)  # noqa: E731
    assert N_TESTS == 53 and group("test_trapezoid") == 14  # the 14 cases of §3.10's first green step
    assert group("test_simpson_exact", "test_simpson_is_richardson", "test_simpson_error", "test_simpson_known",
                 "test_simpson_rejects") == 26  # the 26 Simpson cases that were red before simpson existed
else:  # minimal stand-in: the core contracts, called directly
    N_TESTS = 53
    assert abs(ig.simpson(lambda t: t**3, 0.0, 2.0, 2) - 4.0) < 1e-12
    raises(ValueError, ig.simpson, math.sin, 0.0, 1.0, 3)
    print("Part 1: pytest not installed; ran fallback checks only")

# %% [markdown]
# ## Part 2: convergence orders and error bounds (notes §3.10 Step 0, equations (3.2), (3.3))
# Slope of log(error) against log(h) for the integral of sin over [0, pi]: 2 for the trapezoid rule,
# 4 for Simpson. The constants of (3.2) are derived with SymPy; the bounds give the tolerances.

# %%
NS = [8, 16, 32, 64, 128, 256]
H = np.array([(B - A) / n for n in NS])
err_trap = np.array([abs(ig.trapezoid(math.sin, A, B, n) - EXACT) for n in NS])
err_simp = np.array([abs(ig.simpson(math.sin, A, B, n) - EXACT) for n in NS])
slope_trap, slope_simp = (float(np.polyfit(np.log(H), np.log(e), 1)[0]) for e in (err_trap, err_simp))
print(f"Part 2: slopes trapezoid {slope_trap:.3f}, Simpson {slope_simp:.3f}")
assert abs(slope_trap - 2.0) <= 0.1 and abs(slope_simp - 4.0) <= 0.2
assert 3.9 < err_trap[0] / err_trap[1] < 4.1 and 15.5 < err_simp[0] / err_simp[1] < 16.5

x, h = sp.symbols("x h", positive=True)  # x^4 on [-h, h]: Simpson error (4/15) h^5; x^2 on [0, h]: trapezoid h^3/6
assert sp.simplify(h / 3 * 2 * h**4 - sp.integrate(x**4, (x, -h, h)) - sp.Rational(4, 15) * h**5) == 0
assert sp.simplify(h * h**2 / 2 - sp.integrate(x**2, (x, 0, h)) - h**3 / 6) == 0
bound_n8 = math.pi * (math.pi / 8) ** 4 / 180  # (b-a) h^4 max|f''''| / 180 with max|sin''''| = 1
err_n8 = ig.simpson(math.sin, A, B, 8) - EXACT
assert abs(bound_n8 - 4.15e-4) < 5e-6 and abs(err_n8 - 2.69e-4) < 5e-6 and err_n8 < bound_n8
assert abs(err_n8 / 2 - 1.3e-4) < 1e-5 and abs(bound_n8 / 2 - 2.1e-4) < 1e-5  # the relative figures of Example 4.3
assert bound_n8 / 2 < 2.2e-4 and f"{ig.simpson(math.sin, A, B, 8):.16f}" == "2.0002691699483877"
assert ig.trapezoid(lambda t: 0.1, 0.0, 3.0, 7) == 0.30000000000000004 != 0.3
for m in (2, 5, 12):  # equation (3.3): S_2m = (4 T_2m - T_m) / 3
    rich = (4 * ig.trapezoid(math.exp, 0.0, 1.0, 2 * m) - ig.trapezoid(math.exp, 0.0, 1.0, m)) / 3
    assert abs(ig.simpson(math.exp, 0.0, 1.0, 2 * m) - rich) < 1e-13

u = sp.symbols("u")  # Example 4.4: relative error bound (b-a) h^4 M4 / 180 / |I| at n = 64, six cases
cases = {"sin": (sp.sin(u), 0, sp.pi), "exp": (sp.exp(u), 0, 1), "arctan": (1 / (1 + u**2), 0, 1),
         "gauss-x": (u * sp.exp(-(u**2)), 0, 1), "cos": (sp.cos(u), 0, sp.pi / 2), "log": (sp.log(u), 1, 2)}
rel_bounds = {}
for name, (f_expr, lo, hi) in cases.items():
    f4 = sp.lambdify(u, sp.diff(f_expr, u, 4), "numpy")
    m4 = float(np.max(np.abs(f4(np.linspace(float(lo), float(hi), 20001)))))
    width = float(hi) - float(lo)
    rel_bounds[name] = width * (width / 64) ** 4 * m4 / 180 / abs(float(sp.integrate(f_expr, (u, lo, hi))))
worst = max(rel_bounds, key=lambda k: rel_bounds[k])
assert worst == "sin" and abs(rel_bounds["sin"] - 5.1e-8) < 1e-9 and 19 < 1e-6 / rel_bounds["sin"] < 21

# %% [markdown]
# ## Part 3: the planted off-by-one (notes §1, §3.10)
# `PLANT_OFF_BY_ONE = True` drops the sample y_{n-1}. The whole suite runs again; exactly six test
# functions must fail, no trapezoid test among them, and the weights test only for node 3.

# %%
def simpson_buggy_inline(f: Callable[[float], float], a: float, b: float, n: int) -> float:
    """The §1 version: `ys[1:-2:2]` instead of `ys[1:-1:2]`."""
    step = (b - a) / n
    ys = [f(a + i * step) for i in range(n + 1)]
    return step / 3 * (ys[0] + 4 * sum(ys[1:-2:2]) + 2 * sum(ys[2:-1:2]) + ys[-1])


for n_, shown, err_ in ((12, "1.909708", -9.03e-2), (24, "1.977222", -2.28e-2), (48, "1.994293", -5.71e-3),
                        (96, "1.998572", -1.43e-3)):
    v_inline = simpson_buggy_inline(math.sin, A, B, n_)
    assert f"{v_inline:.6f}" == shown and abs((v_inline - 2.0) - err_) < 5e-5

ig.PLANT_OFF_BY_ONE = True
try:
    assert all(simpson_buggy_inline(math.sin, A, B, n_) == ig.simpson(math.sin, A, B, n_) for n_ in (12, 48))
    bug_err = np.array([abs(ig.simpson(math.sin, A, B, n) - EXACT) for n in NS])
    rel_48 = abs(ig.simpson(math.sin, A, B, 48) - 2.0) / 2.0
    value_100 = ig.simpson(math.sin, A, B, 100)
    code_bug, failed_bug, passed_bug = run_suite() if HAVE_PYTEST else (1, [], [])
finally:
    ig.PLANT_OFF_BY_ONE = False

assert 0.00285 < rel_48 < 0.00286  # "0.29% low" at n = 48
assert 3.9 < bug_err[0] / bug_err[1] < 4.1  # the error falls like h^2 (ratio 4), not h^4 (ratio 16)
assert abs(value_100 - 1.998684) < 5e-7 and abs(value_100 - 2.0) < 1e-2  # a lazy tolerance accepts it ...
assert abs(value_100 - 2.0) / 2.0 > 1e-7  # ... an honest one does not
assert abs((ig.simpson(math.sin, A, B, 100) - 2.0) / 2.0 - 5.4e-9) < 1e-9  # the correct rule at n = 100
if HAVE_PYTEST:
    names = {nid.split("::")[1].split("[")[0] for nid in failed_bug}
    print(f"Part 3: planted bug fails {len(failed_bug)} of {N_TESTS} cases, in {len(names)} test functions")
    assert code_bug != 0 and len(failed_bug) == 22 and len(passed_bug) == N_TESTS - 22 == 31
    assert names == SIX and not any(nm.startswith(("test_trapezoid", "test_adaptive", "test_simpson_rejects"))
                                    for nm in names)
    assert [n for n in failed_bug if "weights" in n] == ["test_integrators.py::test_simpson_weights_are_1_4_2_4_1[3]"]
    assert sum("exact_on_cubics" in n for n in failed_bug) == 12 and sum("known_integrals" in n for n in failed_bug) == 6
assert ig.PLANT_OFF_BY_ONE is False

fig, ax = plt.subplots(figsize=(5.5, 4))
for errs, style, label in ((err_trap, "o-", f"trapezoid, slope {slope_trap:.2f}"),
                           (err_simp, "s-", f"Simpson, slope {slope_simp:.2f}"),
                           (bug_err, "^--", "Simpson with the planted bug (order 2)")):
    ax.loglog(H, errs, style, label=label)
ax.set(xlabel="h", ylabel="|error|", title="PY-5: measured convergence orders")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "convergence.png", dpi=130)
plt.close(fig)

# %% [markdown]
# ## Part 4: exceptions, control flow, EAFP (notes §3.1-§3.6, Examples 4.2 and 4.5)

# %%
def newton(f: Callable[[float], float], df: Callable[[float], float], x0: float,
           tol: float = 1e-12, max_iter: int = 50) -> float:
    """Newton's method for f(x) = 0; raises ConvergenceError, never RuntimeError."""
    xk = x0
    for k in range(1, max_iter + 1):
        try:
            step = f(xk) / df(xk)
        except ZeroDivisionError as err:
            raise ig.ConvergenceError("derivative is zero", k, abs(f(xk))) from err
        xk -= step
        if abs(step) < tol:
            return xk
    raise ig.ConvergenceError("no convergence", max_iter, abs(f(xk)))


sq, dsq = (lambda t: t * t - 2), (lambda t: 2 * t)
assert f"{newton(sq, dsq, 1.0):.10f}" == "1.4142135624"
err = raises(ig.SolverError, newton, sq, dsq, 0.0)
assert type(err) is ig.ConvergenceError and (err.iterations, err.residual) == (1, 2.0)
assert isinstance(err.__cause__, ZeroDivisionError) and str(err) == "derivative is zero (iterations=1, residual=2)"
assert pickle.loads(pickle.dumps(err)).iterations == 1  # all constructor args went to super().__init__
err = raises(ig.ConvergenceError, newton, lambda t: t**3 - 2 * t + 2, lambda t: 3 * t * t - 2, 0.0)
assert (err.iterations, err.residual, err.__cause__) == (50, 2.0, None)  # the cycle 0 -> 1 -> 0 -> ...
assert issubclass(ig.ConvergenceError, ig.SolverError) and not issubclass(ig.ConvergenceError, ValueError)
assert "even" in str(raises(ValueError, ig.simpson, math.sin, 0.0, 1.0, 3))


def naive_trapezoid(f: Callable[[float], float], a: float, b: float, n: int) -> float:
    step = (b - a) / n
    ys = [f(a + i * step) for i in range(n + 1)]
    return step * (0.5 * ys[0] + sum(ys[1:-1]) + 0.5 * ys[-1])


for bad_n, naive_exc in ((0, ZeroDivisionError), (-3, IndexError), (2.5, TypeError), ("4", TypeError)):
    raises(naive_exc, naive_trapezoid, math.sin, 0.0, 1.0, bad_n)  # §3.10 "Green, naively"
    raises(ValueError if isinstance(bad_n, int) else TypeError, ig.trapezoid, math.sin, 0.0, 1.0, bad_n)


def mean_reading(lines: list[str]) -> float:
    values = []
    for line in lines:
        try:
            values.append(float(line))
        except Exception:  # the anti-pattern of §3.3, kept to show the wrong answer
            pass
    return sum(values) / len(values)


assert mean_reading(["1.0", "2.0", "3,5", "4.0"]) == 2.3333333333333335 and abs(np.mean([1, 2, 3.5, 4]) - 2.625) < 1e-15


def swallow() -> str:
    try:
        raise ValueError("this error is about to vanish")
    finally:
        return "x"  # noqa: B012


def attempt(xv: object) -> list[str]:
    trail: list[str] = []
    try:
        trail.append("try")
        if 10 / xv == 2:  # type: ignore[operator]
            return trail + ["return inside try"]
    except ZeroDivisionError:
        trail.append("except")
    else:
        trail.append("else")
    finally:
        trail.append("finally")  # runs after the return value was built, so it is not in it
    return trail


with warnings.catch_warnings():
    warnings.simplefilter("ignore")  # newer Pythons warn about return-in-finally at compile time
    assert swallow() == "x"
assert attempt(4) == ["try", "else", "finally"] and attempt(0) == ["try", "except", "finally"]
assert attempt(5) == ["try", "return inside try"]
raises(TypeError, attempt, "a")


def next_step(h_: float, err_: float, tol_: float) -> float | None:
    if err_ > tol_:
        return h_ / 2
    elif err_ < tol_ / 10:
        return h_ * 2
    return None  # the missing branch of §3.8: falling off the end returns None, silently


def run_steps(h_: float, tol_: float, steps: int) -> float:
    t = 0.0
    for _ in range(steps):
        h_ = next_step(h_, 3.0 * h_**3, tol_)  # type: ignore[assignment]
        t += h_  # type: ignore[operator]
    return t


err_second = 3.0 * 0.05**3  # the second pass of §3.8's pdb session: h = 0.05
assert abs(err_second - 3.75e-4) < 1e-9 and 1e-4 < err_second < 1e-3 and next_step(0.05, err_second, 1e-3) is None
assert "NoneType" in str(raises(TypeError, run_steps, 0.1, 1e-3, 6))


def gen() -> Iterator[int]:
    yield 1
    raise StopIteration


assert str(raises(RuntimeError, list, gen())) == "generator raised StopIteration"  # PEP 479

snippet = "def f(x):\n    assert x >= 0\n    return x**0.5\nprint(f(-4.0))\n"  # assert is removed by -O
plain = subprocess.run([sys.executable, "-c", snippet], capture_output=True, text=True)
opt = subprocess.run([sys.executable, "-O", "-c", snippet], capture_output=True, text=True)
assert plain.returncode != 0 and "AssertionError" in plain.stderr
assert opt.returncode == 0 and opt.stdout.strip() == "(1.2246467991473532e-16+2j)"

with warnings.catch_warnings(record=True) as caught:  # NumPy warns and returns -inf; errstate raises
    warnings.simplefilter("always")
    assert np.log(0.0) == -np.inf
assert issubclass(caught[0].category, RuntimeWarning)
with np.errstate(divide="raise"):
    raises(FloatingPointError, np.log, 0.0)


def old_rule(v: float) -> float:
    warnings.warn("old_rule is deprecated", DeprecationWarning, stacklevel=2)
    return 2 * v


with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    old_rule(3)
assert Path(caught[0].filename).name == "project.py"  # stacklevel=2 blames the caller

with tempfile.TemporaryDirectory() as tmp:  # EAFP versus LBYL (Example 4.2)
    cfg = Path(tmp) / "config.txt"
    cfg.write_text("tol = 1e-8\n")

    def read_lbyl(p: Path) -> str | None:
        if p.exists():
            p.unlink()  # another process deletes the file in the gap
            return p.read_text()
        return None

    def read_eafp(p: Path) -> str | None:
        try:
            return p.read_text()
        except FileNotFoundError:
            return None

    raises(FileNotFoundError, read_lbyl, cfg)
    assert read_eafp(cfg) is None

table = {"a": 1}
stmt = "try:\n    table[{!r}]\nexcept KeyError:\n    pass"
t_hit, t_miss = (min(timeit.repeat(stmt.format(k), globals={"table": table}, number=200_000, repeat=5)) for k in "az")
print(f"Part 4: a dictionary miss under try/except costs {t_miss / t_hit:.1f}x a hit")
assert t_miss > 2 * t_hit

# %% [markdown]
# ## Part 5: context managers (notes §3.7)

# %%
class Timer:
    def __enter__(self) -> "Timer":
        self.start = time.perf_counter()
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        self.elapsed = time.perf_counter() - self.start  # returns None: an exception is not suppressed


@contextmanager
def timer() -> Iterator[SimpleNamespace]:
    box = SimpleNamespace(elapsed=None)
    start = time.perf_counter()
    try:
        yield box
    finally:  # runs even if the body raised
        box.elapsed = time.perf_counter() - start


@contextmanager
def temp_directory() -> Iterator[Path]:
    path = Path(tempfile.mkdtemp(prefix="py5_"))
    try:
        yield path
    finally:
        shutil.rmtree(path)


class Suppress:
    def __init__(self, *excs: type[BaseException]) -> None:
        self.excs = excs

    def __enter__(self) -> None:
        return None

    def __exit__(self, exc_type: type[BaseException] | None, exc: object, tb: object) -> bool:
        return exc_type is not None and issubclass(exc_type, self.excs)  # a true return suppresses


with Timer() as t1:
    time.sleep(0.02)
with timer() as t2:
    time.sleep(0.02)
assert t1.elapsed >= 0.02 and t2.elapsed >= 0.02
kept = None
try:
    with temp_directory() as d:
        (d / "scratch.txt").write_text("x")
        assert d.exists()
        kept = d
        raise RuntimeError("body failed")
except RuntimeError:
    assert kept is not None and not kept.exists()
empty: dict[str, int] = {}
with Suppress(KeyError):
    empty["missing"]
print("Part 5: Timer, timer(), temp_directory() and a suppressing __exit__ behave as in §3.7")

# %% [markdown]
# ## Part 6: logging (notes §3.9, Example 4.6) and the Euler energy numbers
# Equation (3.1): one Euler step multiplies x^2 + v^2 by 1 + h^2. Handlers: a console that shows
# WARNING and up, and a file (`figures/run.log`) that keeps everything.

# %%
xs_, vs_, hs_ = sp.symbols("x v h")
assert sp.expand((xs_ + hs_ * vs_) ** 2 + (vs_ - hs_ * xs_) ** 2 - (1 + hs_**2) * (xs_**2 + vs_**2)) == 0
sim_log = logging.getLogger("sim")


def simulate(step: float, steps: int) -> float:
    """Explicit Euler for x'' = -x, logging as in notes §3.9. Returns the final energy."""
    xk, vk = 1.0, 0.0
    e0 = energy = 0.5 * (xk * xk + vk * vk)
    warned = False
    sim_log.info("start: h=%g steps=%d E0=%.4f", step, steps, e0)
    for k in range(1, steps + 1):
        xk, vk = xk + step * vk, vk - step * xk
        energy = 0.5 * (xk * xk + vk * vk)
        assert abs(energy - e0 * (1 + step**2) ** k) < 1e-12  # equation (3.1), iterated
        sim_log.debug("step %d  x=%+.4f  v=%+.4f  E=%.4f", k, xk, vk, energy)
        if k % 4 == 0:
            sim_log.info("step %d  E=%.4f  (drift %+.1f%%)", k, energy, 100 * (energy / e0 - 1))
        if energy > 1.1 * e0 and not warned:
            sim_log.warning("energy grew by more than 10%% at step %d", k)
            warned = True
    return energy


class ListHandler(logging.Handler):
    def __init__(self, level: int = logging.NOTSET) -> None:
        super().__init__(level)
        self.lines: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.lines.append(f"{record.levelname} {record.getMessage()}")


def capture(logger: logging.Logger, level: int, fn: Callable[[], object], *extra: logging.Handler) -> list[str]:
    """Run fn with a list handler on `logger` at `level`; return the lines it collected."""
    handler = ListHandler()
    for hd in (handler, *extra):
        logger.addHandler(hd)
    logger.setLevel(level)
    logger.propagate = False
    try:
        fn()
    finally:
        for hd in (handler, *extra):
            logger.removeHandler(hd)
        logger.setLevel(logging.NOTSET)
        logger.propagate = True
    return handler.lines


info_lines = capture(sim_log, logging.INFO, lambda: simulate(0.2, 4))
debug_lines = capture(sim_log, logging.DEBUG, lambda: simulate(0.2, 4))
assert len(info_lines) == 3 and len(debug_lines) == 7
assert info_lines[1].startswith("WARNING") and info_lines[1].endswith("step 3")
assert "E=0.5624" in debug_lines[3] and abs(1.04**3 - 1.125) < 1e-3
assert info_lines[2] == "INFO step 4  E=0.5849  (drift +17.0%)"

file_handler = logging.FileHandler(FIG / "run.log", mode="w")  # console shows WARNING+, the file keeps all
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(logging.Formatter("%(levelname)-7s %(name)s: %(message)s"))
console = ListHandler(logging.WARNING)


def run_both() -> None:
    simulate(0.2, 8)
    logging.getLogger("sim.euler").debug("a child logger's record propagates up to sim's handlers")


capture(sim_log, logging.DEBUG, run_both, file_handler, console)
file_handler.close()
log_text = (FIG / "run.log").read_text().splitlines()
n_debug, n_info, n_warn = (sum(ln.startswith(lv) for ln in log_text) for lv in ("DEBUG", "INFO", "WARNING"))
print(f"Part 6: run.log has {n_debug} DEBUG, {n_info} INFO, {n_warn} WARNING lines; console saw {len(console.lines)}")
assert (n_debug, n_info, n_warn) == (9, 3, 1) and len(console.lines) == 1  # 8 steps + 1 child; start, 4, 8
assert any("sim.euler" in ln for ln in log_text) and "36.9%" in " ".join(log_text)

sin_lines = capture(logging.getLogger("integrators"), logging.DEBUG,  # Example 4.6
                    lambda: ig.adaptive_simpson(math.sin, 0.0, math.pi, tol=1e-4))
assert f"{ig.adaptive_simpson(math.sin, 0.0, math.pi, tol=1e-4):.8f}" == "1.99999975"
assert sum(ln.startswith("DEBUG") for ln in sin_lines) == 3 and "|delta|=8.98e-02" in sin_lines[0]
assert "17 evaluations, deepest level 2" in sin_lines[-1]

# %% [markdown]
# ## Part 7: type hints, doctest, lint (notes §3.11-§3.12)

# %%
def double(v: int) -> int:
    return v * 2  # hints are not enforced at run time: double("ab") works


assert double("ab") == "abab"  # type: ignore[arg-type]
for fn in (ig.trapezoid, ig.simpson, ig.adaptive_simpson, ig.ConvergenceError.__init__):
    sig = inspect.signature(fn)
    assert all(p.annotation is not inspect.Parameter.empty for nm, p in sig.parameters.items() if nm != "self"), fn
    assert sig.return_annotation is not inspect.Signature.empty, fn
    typing.get_type_hints(fn)
results = doctest.testmod(ig)
assert results.failed == 0 and results.attempted == 4

VEC_SRC = '''from __future__ import annotations

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
'''
USE_SRC = """from vec_typed import Vec

v = Vec(1, 2, 3)
w = v + Vec(4, 5, 6)  # fine
n: float = abs(w)  # fine
shifted = v + 1  # a mistake
d: Vec = v.dot(w)  # another mistake: dot returns a float
print(w, n)
print(shifted)
"""
demo = FIG / "typed_demo"
demo.mkdir(exist_ok=True)
for fname, src in (("vec_typed.py", VEC_SRC), ("use_vec.py", USE_SRC)):
    (demo / fname).write_text(src)
run_use = subprocess.run([sys.executable, "use_vec.py"], cwd=demo, capture_output=True, text=True)
assert run_use.returncode != 0 and "unsupported operand type(s) for +: 'Vec' and 'int'" in run_use.stderr

if importlib.util.find_spec("mypy") is not None:
    from mypy import api

    cache = str(FIG / ".mypy_cache")
    out, _, status = api.run(["--strict", "--cache-dir", cache, str(HERE / "integrators.py")])
    assert status == 0 and out.startswith("Success"), out
    out, _, status = api.run(["--cache-dir", cache, str(demo / "use_vec.py")])
    use_errs = [ln.split("use_vec.py:")[-1] for ln in out.splitlines() if "use_vec.py:" in ln]
    assert status == 1 and len(use_errs) == 2, out
    assert use_errs[0].startswith('6: error: Unsupported operand types for + ("Vec" and "int")  [operator]')
    assert use_errs[1].startswith('7: error: Incompatible types in assignment (expression has type "float", '
                                  'variable has type "Vec")  [assignment]')
    print("Part 7: mypy --strict is clean on integrators.py; both Vec errors reproduced")
else:
    print("Part 7: mypy not installed; skipped the static checks (pip install mypy)")
if importlib.util.find_spec("ruff") is not None:
    for cmd in (["format", "--check"], ["check", "--select", "E,F,B,W"]):
        res = subprocess.run([sys.executable, "-m", "ruff", *cmd, "--line-length", "100", "integrators.py",
                              "test_integrators.py"], cwd=HERE, capture_output=True, text=True)
        assert res.returncode == 0, res.stdout + res.stderr
    print("        ruff format and ruff check pass on integrators.py and test_integrators.py")

# %% [markdown]
# ## Part 8: adaptive Simpson against SciPy, and the limits of `tol` (notes §8)

# %%
def peak(t: float) -> float:
    return math.exp(-100.0 * (t - 0.5) ** 2)


ref, _ = quad(peak, 0.0, 1.0, epsabs=1e-13, epsrel=1e-13)
assert abs(ref - math.sqrt(math.pi) / 10 * math.erf(5.0)) < 1e-14 and f"{ref:.6f}" == "0.177245"
visited: list[float] = []


def tracked(t: float) -> float:
    visited.append(t)
    return peak(t)


value = ig.adaptive_simpson(tracked, 0.0, 1.0, tol=1e-10)
assert abs(value - ref) < 1e-8, abs(value - ref)
n_uniform = next(n for n in range(2, 4000, 2) if abs(ig.simpson(peak, 0.0, 1.0, n) - ref) < 1e-8)
print(f"Part 8: adaptive error {abs(value - ref):.1e} with {len(visited)} evaluations; "
      f"uniform Simpson needs n = {n_uniform} for 1e-8")
coarse = ig.adaptive_simpson(peak, 0.0, 1.0, tol=1e-3)  # §8: the estimate is fooled by the narrow peak
assert f"{coarse:.6f}" == "0.152564" and abs(coarse - ref) > 0.02
assert abs(ig.adaptive_simpson(peak, 0.0, 1.0, tol=1e-4) - ref) < 1e-5
assert abs(ig.adaptive_simpson(lambda t: t**3 - t, 0.0, 2.0) - 2.0) < 1e-12
t_nan = time.perf_counter()  # NaN fails every comparison, so depth runs out: a misleading message (Your turn 3)
err = raises(ig.ConvergenceError, ig.adaptive_simpson, lambda t: float("nan"), 0.0, 1.0)
assert err.iterations == 40 and math.isnan(err.residual) and time.perf_counter() - t_nan < 1.0
err = raises(ig.ConvergenceError, ig.adaptive_simpson, math.sin, 0.0, 50.0, 1e-14, 3)
assert err.iterations == 3 and err.residual > 1e-14 and isinstance(err, ig.SolverError)

tols = [10.0**-k for k in range(2, 11)]
adaptive_errs = [abs(ig.adaptive_simpson(peak, 0.0, 1.0, tol=tl) - ref) for tl in tols]
assert adaptive_errs[1] > 10 * tols[1]  # tol = 1e-3: the actual error is far above the request
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.5))
grid_x = np.linspace(0, 1, 400)
ax1.plot(grid_x, [peak(t) for t in grid_x], label="integrand")
ax1.plot(visited, np.full(len(visited), -0.04), "|", color="C3", label=f"{len(visited)} adaptive samples")
ax1.set_title("adaptive_simpson puts samples where it must")
ax2.loglog(tols, adaptive_errs, "o-", label="actual error")
ax2.loglog(tols, tols, "k--", label="requested tol")
ax2.invert_xaxis()
ax2.set(xlabel="tol", ylabel="|error|", title="tol is an estimate (fooled at 1e-2, 1e-3)")
for axis in (ax1, ax2):
    axis.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "adaptive.png", dpi=130)
plt.close(fig)

elapsed = time.perf_counter() - T0
print(f"runtime {elapsed:.1f} s")
assert elapsed < 90

# %% [markdown]
# ## Your turn
# 1. Add a test to `test_integrators.py` that would have caught a bug you have actually hit (a
#    wrong sign, an off-by-one in a loop bound, a mutated default). Plant that bug, watch the test
#    fail, remove the bug, watch it pass.
# 2. Write a property-based test for `sorted` with random lists (seeded `np.random.default_rng`):
#    the output is ordered, has the same length, is a permutation of the input (compare
#    `collections.Counter`), and sorting twice changes nothing. Then break a copy of `sorted` on
#    purpose and find the smallest list that exposes the bug by hand (what `hypothesis` automates).
# 3. Make the integrators raise a clear `ValueError` on NaN input (check `a`, `b`, and the sample
#    values). Part 8 shows that `adaptive_simpson` now raises a `ConvergenceError` for a NaN
#    integrand after 40 levels, with a misleading message, and `simpson` quietly returns NaN.
# 4. Put a `breakpoint()` inside `recurse` in `adaptive_simpson`, run
#    `python3 -c "import math, integrators as ig; ig.adaptive_simpson(math.sin, 0, math.pi, 1e-4)"`
#    and inspect `depth`, `lo`, `hi`, `delta` and `eps` at each stop (`c` continues; set
#    `PYTHONBREAKPOINT=0` to disable). Which interval splits last?

# %%
if __name__ == "__main__":
    print("All checks passed.")
