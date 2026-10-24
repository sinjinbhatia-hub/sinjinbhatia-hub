# PY-3 Exercises

Conventions as in the notes: Python 3.11+, standard library only except where an exercise says otherwise. Write your prediction down before you run anything. "Write a function" means: write it, then test it with `assert`s you chose, including an edge case. Several exercises create files; use a scratch directory (`tempfile.TemporaryDirectory()` or a folder you delete afterwards).

## A. Warm-up (do all; 5 to 10 minutes each)

**A1.** (Late binding.) Predict both printed lists, then say what `i` is after the first loop.

```python
fs = []
for i in range(3):
    fs.append(lambda x: x + i)
print([f(10) for f in fs])
gs = [lambda x, i=i: x + i for i in range(3)]
print([g(10) for g in gs])
```

**A2.** (A decorator that forgot something.) Predict the output. Then repair `shout` so that the last two values become `greet` and `Say hello.`, and say which standard-library decorator does the repair.

```python
def shout(func):
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs).upper()
    return wrapper

@shout
def greet(name):
    """Say hello."""
    return f"hello, {name}"

print(greet("ada"), greet.__name__, greet.__doc__)
```

**A3.** (Argument kinds.) For `def f(a, /, b, *, c=3, **kw): return a, b, c, kw`, predict the result or the error of `f(1, 2)`, `f(1, b=2, c=5, z=0)`, `f(*(1, 2), **{"c": 4})`, `f(1, 2, 3)` and `f(a=1, b=2)`.

**A4.** (Counting calls.) Without running anything, how many calls (including the first) do `tri(5)` and `g(3)` make?

```python
def tri(n):
    return 0 if n == 0 else n + tri(n - 1)

def g(n):
    return 1 if n == 0 else g(n - 1) + g(n - 1)
```

Find the formula for $g(n)$'s calls, then confirm your counts with a global counter.

**A5.** (Imports.) File `util.py` contains `print("util loaded", __name__)`, `def f(): return 1`, and `if __name__ == "__main__": print("util main")`. File `main.py` contains `import util`, `from util import f`, `import util as u`, and `print(f(), u is util)`. Predict the output of `python3 main.py`, and of `python3 util.py`.

**A6.** (Write a function.) Write `apply_n(f, n)` that returns the function applying `f` $n$ times (`apply_n(f, 0)` is the identity). Check `apply_n(lambda x: 2 * x, 5)(1) == 32`. What is the type of the thing it returns, and where does `n` live after `apply_n` has returned?

## B. Core (do all; these are the module)

**B1.** (`counted` and `trace`.) Write `counted(func)`, which adds `.calls` and `.reset()` to the wrapped function and keeps its name and docstring. Use it to confirm $C(n) = 2F_{n+1} - 1$ for $n = 0, \dots, 20$ on naive Fibonacci. Then write a `trace` decorator that works for any number of arguments, prints `name(arg1, arg2)` indented by depth, and prints the return value on exit. Run it on `tri(3)` and on `gcd(48, 18)` from Euclid's algorithm.

**B2.** (LRU by hand.) Write `lru(maxsize)`, a decorator factory that caches results in a `collections.OrderedDict`, moves a key to the end on every hit, and evicts the oldest entry when the cache is full. Give the wrapper `hits`, `misses` and `cache_clear()`. On a seeded workload of $10^4$ calls with arguments from `random.Random(0).randint(0, 99)`, check that `(hits, misses)` equals `functools.lru_cache(maxsize=m).cache_info()` for $m = 5, 20, 50$. Why is `OrderedDict` the right container?

**B3.** (A decorator with arguments.) Write `bounded(lo, hi)`: the decorated function's return value must lie in $[lo, hi]$, otherwise raise `ValueError` whose message names the function and the value. Use `functools.wraps`, forward keyword arguments, and test it on a function that returns a root found by `numtools.bisect` (bracket $[0, 1]$) and on one that violates the bound. What happens if someone writes `@bounded` without parentheses, and can you make that error message clearer?

**B4.** (The recursion limit.) (a) Measure the largest `n` for which each of `down(n)` (one argument), `down3(n)` (three local variables) and a `functools.lru_cache`-wrapped version succeeds. Which differ, and why? (b) Write `sum_to(n)` as a plain recursion, show that it fails for $n = 10^4$, then write two fixes: a loop, and a divide-and-conquer recursion on halves that works for $n = 10^6$. State the depth of each.

**B5.** (A package from scratch.) In a temporary directory, create a package `geom` with modules `vectors.py` (`dot`, `norm`) and `stats.py` (`mean`, `variance`), an `__init__.py` that re-exports the four functions and defines `__all__`, and a `__main__.py` that prints a small demonstration. Verify from a subprocess started in the parent directory that `import geom`, `from geom import norm` and `python -m geom` work, and that running `python geom/vectors.py` directly fails when `vectors.py` uses a relative import. Which working directory makes each command succeed?

**B6.** (An `argparse` tool.) Write `summarize.py` with `main(argv=None)` returning an exit status. Subcommands: `mean FILE --column N` and `minmax FILE --column N`, both with `--delimiter` (default `,`) and `--skip-header`. Non-numeric cells are skipped with a warning on stderr and counted; a missing file is an error message on stderr with status 2. Test `main([...])` on a temporary CSV, capturing output with `contextlib.redirect_stdout`.

**B7.** (Extend `project.py`.) In Part 4, add `functools.lru_cache` and the iterative `fib_iter` to the Fibonacci timing study. At $n = 27$ report the ratios (your `memoize` time) : (`lru_cache` time) and (cold `lru_cache`) : (iterative), and assert the order iterative < `lru_cache` < `memoize` < naive. Why is the C-implemented `lru_cache` faster than your dict-based `memoize`, and why is the loop faster than both?

## C. Challenge (do at least one)

**C1.** (Optional arguments.) Write `counted` so that both `@counted` and `@counted(label="fib")` work, and so that it also works on methods (`self` is just the first positional argument). The trick is that the first call receives either the function or keyword arguments. Test both forms and a method. What goes wrong with `@counted(fib)` (a function passed positionally), and how do you make that an explicit error?

**C2.** (A cache-hit-rate study.) Simulate $10^5$ requests with keys drawn (a) uniformly from $N = 1000$ items and (b) from a Zipf-like distribution (weights $1/k$ for key $k = 1..N$, using `random.Random(0).choices`) through `functools.lru_cache(maxsize=m)` for $m = 50, 100, 250, 500$. Argue why the uniform hit rate should be close to $m/N$ once the cache is full, check it within 0.02, and explain the gap for the Zipf workload.

**C3.** (Recursion to iteration, with measurements.) Write a recursive generator `perms(items)` using `yield from`, and an iterative version with an explicit stack of iterators that yields the same permutations in the same order. Check both against `itertools.permutations` for `range(6)`. Measure time for `range(8)` and the maximum depth of the recursive version. Which is fastest, which is shortest, and what does that say about when to convert?
