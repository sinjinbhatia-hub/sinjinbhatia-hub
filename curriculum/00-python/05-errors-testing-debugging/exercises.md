# PY-5 Exercises

Conventions as in the notes: Python 3.11+, `integrators.py` and `test_integrators.py` from this folder (put your files beside them), plain `assert` in tests, run pytest as `python -m pytest -q`. Write your prediction down before you run anything. "Write a test" means: write it first, watch it fail for the right reason, then make it pass.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Control flow.) Predict the printed list, then run it. Explain where each `"f"` comes from and why `2` never appears.

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


print(f([5, 0, "a", 2]))
```

**A2.** (Catch the right thing.) Rewrite `load_tol` so that `""` gives the default `1e-8`, `"1e-6"` gives `1e-6`, any other unparsable text raises a `ValueError` that quotes the text, and `None` raises a `TypeError` (let it propagate). Say what the bare `except:` did to Ctrl-C.

```python
def load_tol(text):
    try:
        return float(text)
    except:
        return 1e-8
```

**A3.** (`assert` or `raise`?) In the function below, one `assert` states an internal invariant and one validates input. Say which is which, rewrite the function, and say what `python -O` does to the original on `cumulative_mean([])`.

```python
def cumulative_mean(xs):
    assert len(xs) > 0, "need data"
    total, out = 0.0, []
    for k, x in enumerate(xs, start=1):
        total += x
        out.append(total / k)
    assert len(out) == len(xs)
    return out
```

**A4.** (Test a small function.) Write `clamp(x, lo, hi)` (limit `x` to `[lo, hi]`; `ValueError` if `lo > hi`) in 5–10 lines, and a `test_clamp.py` with one `parametrize` over at least five cases (inside, below, above, on each boundary) and one `pytest.raises(..., match=...)`.

**A5.** (Discovery.) A folder holds `test_a.py` with `def test_x`, `def check_y`, `class TestZ` (with `def test_w`) and `class Helper` (with `def test_v`), and `b_test.py` with `def test_u`, and `bench.py` with `def test_t`. Predict the output of `pytest --collect-only -q` and the exit status of `pytest -q bench.py` when `test_t` contains `assert False`.

**A6.** (Type hints.) Given `def scale(v: Sequence[float], k: float) -> list[float]`, predict which of these four calls mypy rejects, and why: `scale([1, 2, 3], 2)`, `scale((1.0, 2.0), 2.5)`, `scale("abc", 2)`, `scale([1, 2], "2")`. Check by running mypy.

## B. Core (do all; these are the module)

**B1.** (Test-driven development.) Add a composite **midpoint rule** `midpoint(f, a, b, n)`, $h\sum_{i=0}^{n-1} f(a + (i+\tfrac12)h)$, in the style of `integrators.py` (same validation). Write the tests *first*: exact on lines, observed order $2.0\pm0.1$ on $\sin$, `ValueError` for `n < 1`, and a relation to the trapezoid rule that you derive: for a smooth $f$ the midpoint error is about $-\tfrac12$ of the trapezoid error with the same $n$ (show the one-panel Taylor computation in your write-up). Watch the tests fail, implement, watch them pass, then plant an off-by-one in your rule and report which tests catch it.

**B2.** (Exception hierarchy.) Write `bisect(f, a, b, tol, max_iter)` that raises `BracketError` when $f(a)$ and $f(b)$ have the same sign, and `ConvergenceError` (from `integrators`) when it runs out of iterations. `BracketError` must be catchable both as `SolverError` and as `ValueError` (multiple inheritance). If `f` itself raises `ZeroDivisionError` or `OverflowError`, convert it to a `SolverError` with `raise ... from`. Test every path with `pytest.raises`, checking `__cause__` and the exception's attributes.

**B3.** (Context managers.) Write `cd(path)`, which changes the working directory for a `with` block and restores it afterwards even if the block raises. Implement it twice (class and `@contextmanager`), test both with `tmp_path` (including a body that raises, and nesting), and compare with `contextlib.chdir`.

**B4.** (Logging.) Give the `newton` of §3.4 a module logger. Log each iteration at `DEBUG`, the result at `INFO`, a `WARNING` when $|f'(x)| < 10^{-3}$ at an iterate, and nothing at `ERROR` (the exception is the report). Test with pytest's `caplog` fixture: the levels, the number of records, and that nothing is logged when the logger is not enabled for `DEBUG`.

**B5.** (pdb.) `moving_average([1, 2, 3, 4], 2)` should return `[1.5, 2.5, 3.5]`; the function below returns two values. Without reading the loop twice, use pdb (a conditional breakpoint, `n`, `p`, `w`, `r`) to find the bug. Write down the exact commands you typed and what each showed, fix the function, and add a regression test.

```python
def moving_average(xs, w):
    out = []
    for i in range(w, len(xs)):
        out.append(sum(xs[i - w:i]) / w)
    return out
```

**B6.** (Type hints.) Annotate `mean`, `normalize`, `first_positive` (returns the first positive element or `None`) and `apply_all` (applies a list of `Callable[[float], float]` in turn) so that `mypy --strict` passes. Then: (a) show that passing a `list[int]` to a parameter declared `list[float]` is rejected and explain why (invariance), and fix it with `Sequence[float]`; (b) show the error mypy gives for `first_positive(xs) + 1`, and write the narrowing that fixes it.

**B7.** (Extend `project.py`: mutation testing.) Write four mutants of `simpson` as functions (weights 2 and 4 swapped, `h` computed as `(b - a) / (n + 1)`, the first sample dropped, an even `n` check removed), run the whole suite against each by assigning to `integrators.simpson` inside `try/finally`, and print a table of mutant against the failing test functions. Every mutant must be killed by at least one test; if one survives, write the test that kills it.

## C. Challenge (do at least one)

**C1.** (A property-test runner.) Write `forall(gen, prop, n=200, seed=0)` that draws `n` inputs with a seeded NumPy generator, and on the first failing input **shrinks** it (delete elements, then move numbers toward 0) while the property still fails, and returns the smallest counterexample found. Use it on `sorted` (properties: ordered, same multiset, idempotent) and on a buggy insertion sort whose outer loop is `range(1, len(xs) - 1)`. Report the counterexample before and after shrinking.

**C2.** (Measure the cost of exceptions.) For a dictionary lookup in a loop of 200,000 keys with a miss rate of 0%, 10%, 50% and 90%, time (a) `try/except KeyError`, (b) `if k in d: d[k]`, (c) `d.get(k)`. Tabulate microseconds per lookup, find the miss rate at which (a) stops being the fastest, and write half a page on what this does and does not say about EAFP. State your machine and Python version.

**C3.** (Harden the integrators.) (a) Add `NonFiniteError(SolverError, ValueError)`, raised when `a`, `b`, or any sample is NaN or infinite, for all three rules, without slowing the finite path measurably (measure). Tests first. (b) Fix the coarse-`tol` failure of §8 by adding a `min_depth` argument to `adaptive_simpson` that bisects the interval $2^{\text{min\_depth}}$ times before testing; tabulate the error of $\int_0^1 e^{-100(x-0.5)^2}\,dx$ at `tol=1e-3` for `min_depth` $=0,\dots,4$ against `scipy.integrate.quad`, and explain why the failure disappears.
