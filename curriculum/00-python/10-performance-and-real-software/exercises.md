# PY-10 Exercises

Conventions as in the notes: work from the module folder with `sys.path.insert(0, "nbody/src")` (or `pip install -e nbody` in a venv); time with `time.perf_counter` or `timeit` and report the minimum of repeats; put every pool under `if __name__ == "__main__":`. **Write each prediction down before you run anything**, and turn every "check" into an `assert`. Your timings will differ from the solutions'; your ratios should be of the same order.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (A timing that lies.) A colleague times NumPy's sort of a million doubles like this:

```python
import time
import numpy as np

a = np.random.default_rng(0).random(1_000_000)
t0 = time.time(); a.sort(); first = time.time() - t0
later = []
for _ in range(5):
    t0 = time.time(); a.sort(); later.append(time.time() - t0)
print(f"first {first * 1e3:.1f} ms, mean of the next five {sum(later) / 5 * 1e3:.1f} ms")
```

Predict the two numbers relative to each other. Name three things wrong with this measurement, and rewrite it with `timeit.repeat` so that every repeat sorts the same unsorted data.

**A2.** (Read a profile.) Running `analyse()` below under `cProfile` in this container printed the listing that follows.

```python
def load(n, rng):
    return rng.normal(size=n).tolist()

def smooth(xs, k):
    out = []
    for i in range(len(xs) - k + 1):
        out.append(sum(xs[i:i + k]) / k)
    return out

def fit(ys):
    t = np.arange(len(ys))
    return np.polyfit(t, ys, 3)

def analyse(n=200_000, k=50):
    rng = np.random.default_rng(1)
    return fit(smooth(load(n, rng), k))
```
```text
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.005    0.005    0.533    0.533 a2.py:23(analyse)
        1    0.192    0.192    0.420    0.420 a2.py:11(smooth)
   199951    0.180    0.000    0.180    0.000 {built-in method builtins.sum}
        1    0.000    0.000    0.098    0.098 a2.py:18(fit)
        1    0.041    0.041    0.097    0.097 _polynomial_impl.py:460(polyfit)
   199951    0.048    0.000    0.048    0.000 {method 'append' of 'list' objects}
```

(a) Which function do you optimize first, and what is the largest speedup of `analyse` you could get by making it infinitely fast? The same for `fit`. (b) Where inside `smooth` does its time go? (c) Rewrite `smooth` with NumPy (a cumulative sum or `np.convolve`), assert it equals the original to $10^{-12}$, and measure the new total.

**A3.** (The aliasing bug.) Predict the three printed values, then explain them:

```python
from nbody.forces import ForceWorkspace

work = ForceWorkspace(3)
pos, mass = np.array([[0.0, 0, 0], [1, 0, 0], [0, 2, 0]]), np.ones(3)
history = []
for step in range(3):
    pos = pos + 0.1 * step
    history.append(work(pos, mass))
print(history[0] is history[1], np.shares_memory(history[0], history[2]), np.allclose(history[0], history[2]))
```

Then predict what goes wrong if you pass `accel=ForceWorkspace(2)` to `nbody.integrate.rk4` for a Kepler orbit, check it against `accel=accel_numpy`, and fix it without touching `rk4`.

**A4.** (Preallocate.) Build the array $y_i = \sqrt{i}$, $i = 0, \dots, n-1$, four ways: `np.append` in a loop, a list then `np.array`, `np.empty(n)` filled in a loop, and `np.sqrt(np.arange(n))`. Predict how each time changes when $n$ goes from $10^4$ to $4\times10^4$, then measure.

**A5.** (A missing guard.) Predict what this script prints when run with the `fork` start method (Linux, Python ≤ 3.13) and with `spawn` (macOS, Windows), then fix it. Why must `work` be a module-level function and not a `lambda`?

```python
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

import numpy as np


def work(seed: int) -> float:
    return float(np.random.default_rng(seed).random(1000).mean())


print("starting in", mp.current_process().name)
with ProcessPoolExecutor(max_workers=2) as pool:
    print([round(x, 4) for x in pool.map(work, range(4))])
```

## B. Core (do all; these are the module)

**B1.** (Profile, then fix, a slow program.) The radial distribution of $N$ random points in a unit cube: the histogram of all pair distances below $r_{\max} = 0.5$ in 50 bins, divided by the expected counts for uniform points, $\frac{N(N-1)}{2}\cdot 4\pi r^2\,\Delta r$ (no boundary correction). A first version:

```python
def pair_distances(points):
    out = []
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            out.append(math.dist(points[i], points[j]))
    return out

def histogram(values, edges):
    counts = [0] * (len(edges) - 1)
    for v in values:
        for k in range(len(edges) - 1):
            if edges[k] <= v < edges[k + 1]:
                counts[k] += 1
                break
    return counts
```

Before profiling, guess which function dominates at $N = 600$. Profile, fix the top item, profile again, fix the next, and report a table of total time after each fix. Every version must produce the same counts.

**B2.** (Memory traffic of the potential energy.) Compute $U = -\sum_{i<j} m_i m_j/\sqrt{r_{ij}^2+\varepsilon^2}$ three ways: (a) the full `(N, N, 3)` broadcast and a sum over the upper triangle; (b) coordinate-first `(3, N, N)` differences and `np.triu_indices`; (c) a loop over $i$ with a vectorized sum over $j > i$, which never holds more than $O(N)$ memory. For $N = 500, 1000, 2000$, measure the time and the peak memory (`tracemalloc`), assert the three agree to $10^{-12}$, and explain which you would put in `nbody.integrate.energy` and why.

**B3.** (Threads, processes, NumPy.) Run eight `simulate(seed, n=200, steps=200)` jobs serially, on a `ThreadPoolExecutor(4)` and on a `ProcessPoolExecutor(4)`; then the same three ways for eight calls of the pure-Python `accel_loop` at $N = 200$. Predict all six times before running. Explain each result in terms of the GIL.

**B4.** (A `--profile` flag, with a test.) In a copy of the package, add `--profile PATH` to `nbody/cli.py`: when given, the run happens under `cProfile` and the stats are written to `PATH`. Add a test in `tests/test_cli.py` that calls `main([... "--profile", str(tmp_path / "run.prof")])`, asserts the exit code is 0 and the file exists, and that `pstats.Stats` can read it and finds `leapfrog` among its functions. Run the whole suite.

**B5.** (float32 for a long orbit.) Integrate the $e = 0.5$ Kepler orbit of `project.py` Part 3 for 50 orbits with leapfrog at 1000 steps per orbit, once in float64 and once in float32. Report the relative energy error and the error in the angle of the relative position at the end (against the float64 run with 4000 steps per orbit as reference). Which error does float32 damage more, and why?

**B6.** (A git history for a bug fix.) In a temporary copy of the package: commit it; make a commit that silently flips the sign of the force in `accel_numpy`; make two more harmless commits. Then let `git bisect run python -m pytest -q` find the bad commit, and fix it with a commit whose message says what broke and why the tests caught it. Paste the transcript.

## C. Challenge (do at least one)

**C1.** (Is a cutoff good enough for gravity?) For $N = 4000$ bodies of `uniform_ball`, compute the accelerations with `accel_short_range(method="tree")` for cutoffs $r_c = 0.1, 0.2, 0.4, 0.8, 1.6$ and compare with the full `accel_numpy`: median relative error per body and time. Which $r_c$ gives 1% accuracy, and what does it cost against the workspace? Explain the result physically, and what Barnes–Hut does instead.

**C2.** (A model for `chunksize`.) Time `pool.map(math.sqrt, xs, chunksize=c)` over 10,000 floats with 4 workers for $c = 1, 4, 16, 64, 256, 1024, 2500$. Fit $T(c) = T_0 + a\,n/c$ by least squares on the measurements and interpret $T_0$ and $a$. Where does the model fail, and why?

**C3.** (Numba in parallel.) Write `accel_numba_parallel` with `@njit(parallel=True)` and `numba.prange` over $i$. Check it against `accel_loop`, then time it at $N = 200$ and $N = 2000$ with `numba.set_num_threads(t)` for $t = 1, 2, 4$, against the serial `accel_numba` and the workspace. Where does parallelism pay, and why not everywhere?
