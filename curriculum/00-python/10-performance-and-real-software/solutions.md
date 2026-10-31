# PY-10 Solutions

Every block below was run, in order, in one session from the module folder (Python 3.13.16, NumPy 2.5.3, numba 0.68.1, 4 CPUs in this shared container), and its output pasted from that run; its asserts passed. Blocks marked as scripts were saved to a file and run from the module folder with `PYTHONPATH=nbody/src python3 file.py`, because process pools need the `__main__` guard. **Timings are the writer's measurements on this container**; yours will differ, and the ratios are what to compare. The preamble:

```python
import math
import os
import subprocess
import sys
import tempfile
import time
import timeit
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path("nbody/src").resolve()))
from nbody.forces import ForceWorkspace, accel_loop, accel_numba, accel_numpy, accel_short_range
from nbody.integrate import energy, leapfrog, rk4, simulate, uniform_ball


def bench(candidates: dict, number: int = 1, repeat: int = 5) -> dict:
    """Seconds per call: minimum over `repeat` interleaved rounds of `number` calls."""
    best = {name: math.inf for name in candidates}
    for _ in range(repeat):
        for name, f in candidates.items():
            t0 = time.perf_counter()
            for _ in range(number):
                f()
            best[name] = min(best[name], (time.perf_counter() - t0) / number)
    return best


def two_body(a: float = 1.0, e: float = 0.5, m1: float = 0.8, m2: float = 0.2, dtype=np.float64):
    """Two bodies at pericentre of a Kepler orbit (semi-major axis a, eccentricity e), CM at rest."""
    big_m = m1 + m2
    rp, vp = a * (1 - e), math.sqrt(big_m * (1 + e) / (a * (1 - e)))
    pos = np.array([[-m2 / big_m * rp, 0, 0], [m1 / big_m * rp, 0, 0]], dtype=dtype)
    vel = np.array([[0, -m2 / big_m * vp, 0], [0, m1 / big_m * vp, 0]], dtype=dtype)
    return pos, vel, np.array([m1, m2], dtype=dtype)


print(sys.version.split()[0], np.__version__)
```
```text
3.13.16 2.5.3
```

## A. Warm-up

**A1.**

```python
a = np.random.default_rng(0).random(1_000_000)
t0 = time.time(); a.sort(); first = time.time() - t0
later = []
for _ in range(5):
    t0 = time.time(); a.sort(); later.append(time.time() - t0)
print(f"as written: first {first * 1e3:.1f} ms, mean of the next five {sum(later) / 5 * 1e3:.1f} ms")

raw = np.random.default_rng(0).random(1_000_000)
fair = timeit.repeat("b.sort()", setup="b = raw.copy()", globals={"raw": raw}, number=1, repeat=7)
print(f"fair: min {min(fair) * 1e3:.1f} ms, max {max(fair) * 1e3:.1f} ms over 7 sorts of fresh unsorted copies")
```
```text
as written: first 8.7 ms, mean of the next five 7.3 ms
fair: min 7.3 ms, max 9.0 ms over 7 sorts of fresh unsorted copies
```

Here the two numbers came out close (8.7 and 7.3 ms), which is the danger: the measurement looks plausible while the last five sorts time a different problem, because after the first call the array *is sorted*. Whether sorted input is faster or slower depends on the algorithm (classic quicksort degrades on it; NumPy's sort here barely notices), so you cannot reason your way to the right number; you have to measure the right thing. The three faults: (1) the timed operation mutates its input, so repeats do not repeat the same work; (2) `time.time` is the calendar clock, adjustable and coarser than `perf_counter`; (3) a mean folds noise into the estimate, where the minimum is the best estimate of cost. `timeit.repeat` with a `setup` statement runs the setup before every repeat, so with `number=1` each timing sorts a fresh unsorted copy, and the copy is not timed. The fair minimum (7.3 ms) is what to report.

**A2.**

(a) `smooth`: cumtime 0.420 s of 0.533 s, 79% of the run. Making it infinitely fast leaves $0.533 - 0.420 = 0.113$ s, a speedup of at most $0.533/0.113 = 4.7\times$. Making `fit` infinitely fast leaves $0.533 - 0.098 = 0.435$ s: at most $1.23\times$. This bound, $1/(1-f)$ for a part taking fraction $f$, is Amdahl's law, and it is why you optimize the top of the cumtime list.

(b) `smooth` has tottime 0.192 s (its own loop, slicing `xs[i:i + k]` into a new 50-element list each time, the division) and calls `sum` 199,951 times (0.180 s) and `list.append` 199,951 times (0.048 s); the three add up to its cumtime. The slices and the 200,000 calls into `sum` are the cost, not the additions.

(c) A moving average is a difference of cumulative sums: $\bar y_i = (S_{i+k} - S_i)/k$ with $S_j = \sum_{l<j} x_l$.

```python
import cProfile
import pstats


def load(n: int, rng: np.random.Generator) -> list[float]:
    return rng.normal(size=n).tolist()


def smooth(xs: list[float], k: int) -> list[float]:
    out = []
    for i in range(len(xs) - k + 1):
        out.append(sum(xs[i:i + k]) / k)
    return out


def smooth_np(xs, k: int) -> np.ndarray:
    """Moving average of width k: (S[i+k] - S[i]) / k with S the cumulative sum, S[0] = 0."""
    s = np.concatenate(([0.0], np.cumsum(xs)))
    return (s[k:] - s[:-k]) / k


def fit(ys) -> np.ndarray:
    return np.polyfit(np.arange(len(ys)), ys, 3)


xs = load(200_000, np.random.default_rng(1))
assert np.allclose(smooth_np(xs, 50), smooth(xs, 50), rtol=0, atol=1e-12)
t = bench({"load": lambda: load(200_000, np.random.default_rng(1)), "smooth": lambda: smooth(xs, 50),
           "smooth_np": lambda: smooth_np(xs, 50), "fit": lambda: fit(smooth_np(xs, 50))}, repeat=3)
old_total, new_total = t["load"] + t["smooth"] + t["fit"], t["load"] + t["smooth_np"] + t["fit"]
print("unprofiled: " + ", ".join(f"{k} {v * 1e3:.1f} ms" for k, v in t.items()))
print(f"analyse: {old_total * 1e3:.0f} ms -> {new_total * 1e3:.0f} ms ({old_total / new_total:.1f}x); "
      f"smooth's unprofiled share was {t['smooth'] / old_total:.0%}")
```
```text
unprofiled: load 7.3 ms, smooth 118.7 ms, smooth_np 7.5 ms, fit 30.5 ms
analyse: 156 ms -> 45 ms (3.5x); smooth's unprofiled share was 76%
```

Without the profiler, `smooth` was 76% of the run (cProfile's per-call overhead had inflated it to 79%), so Amdahl's bound for removing it was $1/(1 - 0.76) = 4.2\times$; replacing it with a 7.5 ms cumulative sum gave $3.5\times$ (156 ms → 45 ms). `fit` is now two thirds of what remains, and a cubic fit to 200,000 points is already one least-squares solve: this is where you stop. The cumulative-sum trick changes the summation order, which is why the check is to $10^{-12}$ and not `==`.

**A3.**

```python
work = ForceWorkspace(3)
pos3, mass3 = np.array([[0.0, 0, 0], [1, 0, 0], [0, 2, 0]]), np.ones(3)
history = []
for step in range(3):
    pos3 = pos3 + 0.1 * step
    history.append(work(pos3, mass3))
print(history[0] is history[1], np.shares_memory(history[0], history[2]), np.allclose(history[0], history[2]))

p2, v2, m2 = two_body()
period = 2 * np.pi
e0 = energy(p2, v2, m2)
_, v_ok, e_ok = rk4(p2, v2, m2, period / 200, 200, accel=accel_numpy, record_every=200)
_, v_bad, e_bad = rk4(p2, v2, m2, period / 200, 200, accel=ForceWorkspace(2), record_every=200)
w2 = ForceWorkspace(2)
_, _, e_fix = rk4(p2, v2, m2, period / 200, 200, accel=lambda p, m, eps: w2(p, m, eps).copy(), record_every=200)
for name, e in (("accel_numpy", e_ok), ("ForceWorkspace", e_bad), ("workspace + copy", e_fix)):
    print(f"RK4, one orbit, {name:<17}: |dE/E| = {abs(e[-1] / e0 - 1):.1e}")
assert abs(e_fix[-1] - e_ok[-1]) < 1e-14
```
```text
False True True
RK4, one orbit, accel_numpy      : |dE/E| = 9.3e-07
RK4, one orbit, ForceWorkspace   : |dE/E| = 8.5e+00
RK4, one orbit, workspace + copy : |dE/E| = 9.3e-07
```

`False True True`. Each call returns `self.acc.T`, a *new view object* (so `is` is `False`) of the *same buffer* (so `shares_memory` is `True`), and since every entry of `history` is a view of that one buffer, all three hold the last forces (`allclose` is `True`). In `rk4` the same thing happens inside one step: `k1v` to `k4v` all end up as views of the buffer holding the last evaluation, so the velocity update becomes $\mathbf v + \Delta t\,\mathbf k_4$, and after one orbit the energy is off by 850% instead of $10^{-6}$, with no error message. The positions are still updated correctly, which keeps the bug quiet. The fix wraps the workspace in a function that copies, and agrees with `accel_numpy` to $10^{-14}$; `rk4` itself is untouched.

**A4.**

```python
def by_append(n: int) -> np.ndarray:
    y = np.array([])
    for i in range(n):
        y = np.append(y, math.sqrt(i))
    return y


def by_list(n: int) -> np.ndarray:
    return np.array([math.sqrt(i) for i in range(n)])


def by_fill(n: int) -> np.ndarray:
    y = np.empty(n)
    for i in range(n):
        y[i] = math.sqrt(i)
    return y


def by_vector(n: int) -> np.ndarray:
    return np.sqrt(np.arange(n))


rows = {}
for n in (10_000, 40_000):
    rows[n] = bench({f.__name__: (lambda f=f: f(n)) for f in (by_append, by_list, by_fill, by_vector)}, repeat=3)
    assert np.array_equal(by_append(1000), by_vector(1000)) and np.array_equal(by_fill(n), by_list(n))
for name in rows[10_000]:
    print(f"{name:<10} n=1e4 {rows[10_000][name] * 1e3:8.2f} ms   n=4e4 {rows[40_000][name] * 1e3:8.2f} ms   "
          f"ratio {rows[40_000][name] / rows[10_000][name]:5.1f}")
```
```text
by_append  n=1e4    29.53 ms   n=4e4   249.49 ms   ratio   8.4
by_list    n=1e4     0.79 ms   n=4e4     3.17 ms   ratio   4.0
by_fill    n=1e4     0.80 ms   n=4e4     3.10 ms   ratio   3.9
by_vector  n=1e4     0.06 ms   n=4e4     0.12 ms   ratio   2.0
```

`np.append` allocates a new array and copies everything each time, so building $n$ elements costs $\sum_i i = O(n^2)$ in copying, on top of a fixed cost per call. Quadrupling $n$ multiplied its time by 8.4 here, on the way to 16 as the copying term takes over, while the list and the fill scale by 4 (linear), and the vectorized line by 2 (at $10^4$ elements it is still mostly fixed overhead). Filling a preallocated array element by element is no faster than the list, because each `y[i] = ...` is still one Python-level operation; only the vectorized line removes the loop, and it is 13 to 26 times faster than either.

**A5.**

```python
script = '''import multiprocessing as mp
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np


def work(seed: int) -> float:
    return float(np.random.default_rng(seed).random(1000).mean())


print("starting in", mp.current_process().name)
with ProcessPoolExecutor(max_workers=2, mp_context=mp.get_context(sys.argv[1])) as pool:
    print([round(x, 4) for x in pool.map(work, range(4))])
'''
fixed = script.replace('print("starting in"', 'if __name__ == "__main__":\n    print("starting in"').replace(
    "\nwith ProcessPoolExecutor", "\n    with ProcessPoolExecutor").replace("\n    print([round", "\n        print([round")
with tempfile.TemporaryDirectory() as tmp:
    for name, text in (("unguarded", script), ("guarded", fixed)):
        path = Path(tmp) / f"{name}.py"
        path.write_text(text)
        for method in ("fork", "spawn"):
            out = subprocess.run([sys.executable, str(path), method], capture_output=True, text=True, timeout=120)
            lines = (out.stdout + out.stderr).splitlines()
            shown = [ln for ln in lines if ln.startswith(("starting", "[", "concurrent"))]
            print(f"{name}, {method}: exit {out.returncode}; " + " | ".join(shown))
print(fixed.split("\n\n\n")[-1])
```
```text
unguarded, fork: exit 0; starting in MainProcess | [0.5169, 0.5028, 0.5022, 0.4944]
unguarded, spawn: exit 1; starting in MainProcess | starting in SpawnProcess-1 | starting in SpawnProcess-2 | concurrent.futures.process.BrokenProcessPool: A process in the process pool was terminated abruptly while the future was running or pending.
guarded, fork: exit 0; starting in MainProcess | [0.5169, 0.5028, 0.5022, 0.4944]
guarded, spawn: exit 0; starting in MainProcess | [0.5169, 0.5028, 0.5022, 0.4944]
if __name__ == "__main__":
    print("starting in", mp.current_process().name)
    with ProcessPoolExecutor(max_workers=2, mp_context=mp.get_context(sys.argv[1])) as pool:
        print([round(x, 4) for x in pool.map(work, range(4))])
```

With `fork` the unguarded script works: workers are copies of the parent made *after* the top level ran, so nothing is re-executed. With `spawn` each worker imports the script as `__mp_main__` and re-runs its top level, printing "starting in SpawnProcess-1" and "-2", and tries to start a pool of its own; Python refuses with the `RuntimeError` of notes §3.6, and the parent gets a `BrokenProcessPool`. The guarded version prints the same four means under both. `work` must be a module-level function because the pool sends it to the worker by *pickling a reference* (module name plus function name) that the worker re-imports; a `lambda` has no importable name and fails to pickle.

## B. Core

**B1.** The guess most people make is `pair_distances`, the $O(N^2)$ double loop. The profile:

```python
import random


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


def histogram_direct(values, rmax: float, nbins: int):
    """Fix 1: the bin of v is int(v / dr); no search over bins."""
    counts = [0] * nbins
    dr = rmax / nbins
    for v in values:
        if v < rmax:
            counts[int(v / dr)] += 1
    return counts


def rdf(points, rmax=0.5, nbins=50, version=0):
    n = len(points)
    dr = rmax / nbins
    if version == 0:
        counts = histogram(pair_distances(points), [k * dr for k in range(nbins + 1)])
    elif version == 1:
        counts = histogram_direct(pair_distances(points), rmax, nbins)
    else:                                   # fix 2: distances and histogram in NumPy
        from scipy.spatial.distance import pdist
        counts = np.histogram(pdist(np.asarray(points)), bins=nbins, range=(0, rmax))[0]
    r = (np.arange(nbins) + 0.5) * dr
    return np.asarray(counts) / (n * (n - 1) / 2 * 4 * np.pi * r**2 * dr)


random.seed(3)
points = [[random.random() for _ in range(3)] for _ in range(600)]
prof = cProfile.Profile()
prof.enable()
g0 = rdf(points, version=0)
prof.disable()
pstats.Stats(prof).strip_dirs().sort_stats("cumulative").print_stats(5)
t = {v: bench({"r": lambda v=v: rdf(points, version=v)}, repeat=2)["r"] for v in (0, 1, 2)}
for v in (1, 2):
    assert np.allclose(rdf(points, version=v), g0, rtol=0, atol=1e-12)
print(f"original {t[0]:.3f} s | fix 1 (direct bin index) {t[1]:.3f} s ({t[0] / t[1]:.1f}x) | "
      f"fix 2 (pdist + np.histogram) {t[2] * 1e3:.1f} ms ({t[0] / t[2]:.0f}x)")
print("g(r) in the first five bins:", np.round(g0[:5], 3))
```
```text
         539709 function calls in 0.758 seconds

   Ordered by: cumulative time
   List reduced from 9 to 5 due to restriction <5>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.002    0.002    0.758    0.758 <block 7>:32(rdf)
        1    0.614    0.614    0.636    0.636 <block 7>:12(histogram)
        1    0.063    0.063    0.120    0.120 <block 7>:4(pair_distances)
   179700    0.035    0.000    0.035    0.000 {built-in method math.dist}
   179700    0.022    0.000    0.022    0.000 {method 'append' of 'list' objects}


original 0.583 s | fix 1 (direct bin index) 0.040 s (14.7x) | fix 2 (pdist + np.histogram) 7.9 ms (73x)
g(r) in the first five bins: [0.    0.984 0.921 0.976 0.962]
```

(The functions were defined in the solution session, hence the `<block 7>` file names.) The profile says the histogram (0.636 s of 0.758 s), not the distances (0.120 s): for each of the 179,700 distances, `histogram` walks the bin list until it finds the right bin, so its cost is $O(N^2\times\text{bins})$. Fix 1 replaces the search with arithmetic (the bin is `int(v / dr)`) and gives 15×; the profile then points at `pair_distances`, and fix 2 moves both remaining loops into compiled code (`pdist` returns the condensed distance vector, `np.histogram` bins it) for 73× in total. All three versions give the same $g(r)$. The first bin is empty because its shell, $r < 0.01$, expects only about 0.6 pairs; $g(r)$ sits a little below 1 because the cube's walls remove neighbours, the boundary correction the problem told you to ignore.

**B2.**

```python
import tracemalloc


def u_broadcast(pos, mass, eps):
    d = pos[np.newaxis, :, :] - pos[:, np.newaxis, :]
    r = np.sqrt((d * d).sum(axis=-1) + eps**2)
    return -np.sum(np.triu(np.outer(mass, mass) / r, k=1))


def u_triu(pos, mass, eps):
    x = np.ascontiguousarray(pos.T)
    i, j = np.triu_indices(len(mass), k=1)
    r2 = sum((x[k, j] - x[k, i]) ** 2 for k in range(3))
    return -np.sum(mass[i] * mass[j] / np.sqrt(r2 + eps**2))


def u_rows(pos, mass, eps):
    total = 0.0
    for i in range(len(mass) - 1):
        d = pos[i + 1:] - pos[i]
        total -= mass[i] * np.sum(mass[i + 1:] / np.sqrt(np.einsum("jk,jk->j", d, d) + eps**2))
    return total


def peak_mb(f) -> float:
    tracemalloc.start()
    f()
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    return peak / 1e6


for n in (500, 1000, 2000):
    p, _, m = uniform_ball(n, np.random.default_rng(n))
    vals = [f(p, m, 0.05) for f in (u_broadcast, u_triu, u_rows)]
    assert max(abs(v / vals[0] - 1) for v in vals) < 1e-12
    t = bench({f.__name__: (lambda f=f: f(p, m, 0.05)) for f in (u_broadcast, u_triu, u_rows)}, repeat=3)
    print(f"N = {n}: " + "; ".join(f"{f.__name__} {t[f.__name__] * 1e3:6.1f} ms, peak {peak_mb(lambda f=f: f(p, m, 0.05)):6.1f} MB"
                                   for f in (u_broadcast, u_triu, u_rows)))
```
```text
N = 500: u_broadcast    9.8 ms, peak   14.0 MB; u_triu    4.2 ms, peak    6.0 MB; u_rows    6.1 ms, peak    0.0 MB
N = 1000: u_broadcast   49.4 ms, peak   56.0 MB; u_triu   20.0 ms, peak   24.0 MB; u_rows   16.0 ms, peak    0.1 MB
N = 2000: u_broadcast  204.0 ms, peak  224.0 MB; u_triu   89.7 ms, peak   96.0 MB; u_rows   48.8 ms, peak    0.1 MB
```

The broadcast holds several $N\times N$ and $N\times N\times 3$ arrays at once: its peak memory grows as $N^2$ (14, 56, 224 MB). The `triu_indices` version computes only the $N(N-1)/2$ pairs it needs but stores index arrays and per-pair temporaries of that length, so it is still $O(N^2)$ in memory, at less than half the broadcast's peak. The row loop runs $N$ Python iterations, each a vectorized pass over at most $N$ elements: its peak is $O(N)$ (0.1 MB), and it is the *fastest* from $N = 1000$ on, because $N$ Python iterations are cheap next to $N^2/2$ pair terms and its arrays stay in cache (PY-6 §3.10, "loop over the short axis"). For `nbody.integrate.energy`, called rarely and required not to fail at large $N$, the row loop is the right choice; the package uses the broadcast, and this exercise is the evidence for changing it.

**B3.** Script (`b3.py`):

<!-- file: b3.py -->
```python
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import numpy as np

from nbody.forces import accel_loop
from nbody.integrate import simulate, uniform_ball


def numpy_job(seed: int) -> float:
    return simulate(seed, n=200, steps=200)["rel_energy_change"]


def python_job(seed: int) -> float:
    pos, _, mass = uniform_ball(200, np.random.default_rng(seed))
    return accel_loop(pos.tolist(), mass.tolist(), 0.05)[0][0]


def timed(executor_cls, job) -> tuple[float, list[float]]:
    t0 = time.perf_counter()
    if executor_cls is None:
        out = [job(s) for s in range(8)]
    else:
        with executor_cls(max_workers=4) as ex:
            out = list(ex.map(job, range(8)))
    return time.perf_counter() - t0, out


if __name__ == "__main__":
    for job in (numpy_job, python_job):
        t_s, ref = timed(None, job)
        t_t, out_t = timed(ThreadPoolExecutor, job)
        t_p, out_p = timed(ProcessPoolExecutor, job)
        assert out_t == ref and out_p == ref
        print(f"{job.__name__:<11} serial {t_s:.2f} s | 4 threads {t_t:.2f} s ({t_s / t_t:.1f}x) | "
              f"4 processes {t_p:.2f} s ({t_s / t_p:.1f}x)")
```
```text
numpy_job   serial 0.94 s | 4 threads 0.58 s (1.6x) | 4 processes 0.30 s (3.1x)
python_job  serial 0.10 s | 4 threads 0.12 s (0.8x) | 4 processes 0.06 s (1.7x)
```

For the pure-Python force loop, threads lose (0.8×: they take turns holding the GIL) and processes give 1.7×, less than 4 because each 12-ms job is short next to the pool's start-up. For the NumPy simulation, threads help partly (1.6×): inside each force evaluation the ufunc loops over 40,000-element arrays release the GIL, so several threads compute at once, but every step also runs Python bytecode (the leapfrog updates, the calls), during which the threads queue for the lock. Processes have no lock to share and give 3.1×. How much threads help NumPy code depends on how much of its time is inside large array operations; measure it, as here.

**B4.**

```python
patch_cli = [
    ('    p.add_argument("-v", "--verbose"',
     '    p.add_argument("--profile", type=Path, help="run under cProfile and write the stats to this file")\n'
     '    p.add_argument("-v", "--verbose"'),
    ('''    results = run_ensemble(seeds, workers=cfg.workers, n=cfg.n, steps=cfg.steps, dt=cfg.dt,
                           eps=cfg.eps, dtype=cfg.dtype)''',
     '''    run = lambda: run_ensemble(seeds, workers=cfg.workers, n=cfg.n, steps=cfg.steps, dt=cfg.dt,
                               eps=cfg.eps, dtype=cfg.dtype)
    if args.profile:
        import cProfile
        prof = cProfile.Profile()
        results = prof.runcall(run)
        prof.dump_stats(args.profile)
        log.info("profile written to %s", args.profile)
    else:
        results = run()'''),
]
test_cli = '''import pstats

from nbody.cli import main


def test_profile_flag_writes_readable_stats(tmp_path):
    out = tmp_path / "run.prof"
    assert main(["--n", "20", "--steps", "10", "--profile", str(out)]) == 0
    assert out.exists()
    names = {func[2] for func in pstats.Stats(str(out)).stats}
    assert "leapfrog" in names
'''
with tempfile.TemporaryDirectory() as tmp:
    pkg = Path(tmp) / "nbody"
    subprocess.run(["cp", "-r", "nbody", str(pkg)], check=True)
    cli_file = pkg / "src" / "nbody" / "cli.py"
    text = cli_file.read_text()
    for old, new in patch_cli:
        assert text.count(old) == 1
        text = text.replace(old, new)
    cli_file.write_text(text)
    (pkg / "tests" / "test_cli.py").write_text(test_cli)
    out = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"], cwd=pkg,
                         capture_output=True, text=True)
    print(out.stdout.strip().splitlines()[-1])
    run = subprocess.run([sys.executable, "-m", "nbody", "--n", "200", "--steps", "100", "--profile", "run.prof"],
                         cwd=pkg, capture_output=True, text=True, env={**os.environ, "PYTHONPATH": str(pkg / "src")})
    print(run.stderr.strip().splitlines()[-1])
    stats = pstats.Stats(str(pkg / "run.prof"))
    stats.strip_dirs().sort_stats("tottime").print_stats(3)
```
```text
15 passed in 1.12s
INFO nbody: 1 run(s), N = 200, 100 steps, 1 worker(s): 0.06 s
Sat Oct 31 17:12:38 2026    /tmp/tmpq0xktp65/nbody/run.prof

         1154 function calls (1152 primitive calls) in 0.056 seconds

   Ordered by: internal time
   List reduced from 90 to 3 due to restriction <3>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
      101    0.041    0.000    0.047    0.000 forces.py:83(__call__)
      113    0.007    0.000    0.007    0.000 {method 'reduce' of 'numpy.ufunc' objects}
        2    0.004    0.002    0.006    0.003 integrate.py:37(energy)
```

`Profile.runcall(f)` profiles exactly one call, which keeps argument parsing and logging out of the profile; `dump_stats` writes the binary file that `python -m pstats run.prof` (or `snakeviz`) reads. The test uses pytest's `tmp_path` fixture (PY-5 §3.10) so it never writes into the repository, and it checks that the file is *readable* and *contains the integrator*, not merely that it exists. The profile of the real run says what Build 3 predicts: 0.041 s of 0.056 s is in the workspace's `__call__` (the ufuncs it calls are C functions, counted in its own time), and the energy diagnostic is next.

**B5.**

```python
def angle_and_energy(dtype, per_orbit: int, orbits: int = 50):
    p, v, m = two_body(dtype=dtype)
    e0 = energy(p, v, m)
    p, v, _ = leapfrog(p, v, m, 2 * np.pi / per_orbit, orbits * per_orbit)
    rel = (p[1] - p[0]).astype(np.float64)
    return math.atan2(rel[1], rel[0]), abs(energy(p, v, m) / e0 - 1)


ref_angle, _ = angle_and_energy(np.float64, 4000)
angles = {}
for dtype in (np.float64, np.float32):
    angles[dtype], de = angle_and_energy(dtype, 1000)
    print(f"{np.dtype(dtype).name}: |dE/E| = {de:.1e}, angle error after 50 orbits = {abs(angles[dtype] - ref_angle):.2e} rad")
print(f"float32 minus float64, same step: {abs(angles[np.float32] - angles[np.float64]):.2e} rad")
```
```text
float64: |dE/E| = 3.7e-06, angle error after 50 orbits = 1.65e-01 rad
float32: |dE/E| = 3.6e-05, angle error after 50 orbits = 2.15e-01 rad
float32 minus float64, same step: 4.94e-02 rad
```

In float64 at 1000 steps per orbit, the energy error is leapfrog's bounded oscillation ($3.7\times10^{-6}$), and the angle error after 50 orbits, 0.165 rad, is the method's own $O(\Delta t^2)$ period error accumulating linearly in time. float32 multiplies the energy error by 10 ($3.6\times10^{-5}$), and adds 0.049 rad of phase: rounding positions and velocities at relative size $10^{-7}$ every step perturbs the energy, hence the semi-major axis and the period (Kepler's third law), and a period error becomes a phase error that grows every orbit. The energy is damaged more in relative terms; the phase error is the one that grows without bound. In a long integration, errors that change a *rate* accumulate, and float32 is the wrong place to save memory.

**B6.**

```python
def git(*args: str, cwd: Path) -> str:
    out = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True)
    return (out.stdout + out.stderr).strip()


with tempfile.TemporaryDirectory() as tmp:
    pkg = Path(tmp) / "nbody"
    subprocess.run(["cp", "-r", "nbody", str(pkg)], check=True)
    for d in pkg.rglob("__pycache__"):
        subprocess.run(["rm", "-rf", str(d)], check=True)
    forces = pkg / "src" / "nbody" / "forces.py"
    readme = pkg / "README.md"
    steps = [("init", "-q", "-b", "main"), ("config", "user.name", "You"), ("config", "user.email", "you@example.com"),
             ("add", "."), ("commit", "-q", "-m", "Import the nbody package")]
    for s in steps:
        git(*s, cwd=pkg)
    forces.write_text(forces.read_text().replace("d = x[:, np.newaxis, :] - x[:, :, np.newaxis]",
                                                 "d = x[:, :, np.newaxis] - x[:, np.newaxis, :]"))
    git("commit", "-q", "-am", "Tidy the broadcast in accel_numpy", cwd=pkg)
    readme.write_text(readme.read_text() + "\nSee notes PY-10 for the derivation of every kernel.\n")
    git("commit", "-q", "-am", "README: point to the notes", cwd=pkg)
    readme.write_text(readme.read_text() + "Units: G = 1, total mass 1.\n")
    git("commit", "-q", "-am", "README: state the units", cwd=pkg)
    print("$ git log --oneline\n" + git("log", "--oneline", cwd=pkg))
    print("$ git bisect start HEAD HEAD~3\n" + git("bisect", "start", "HEAD", "HEAD~3", cwd=pkg))
    out = git("bisect", "run", sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", cwd=pkg)
    print("$ git bisect run python -m pytest -q -x\n" + "\n".join(
        ln for ln in out.splitlines() if "is the first bad commit" in ln or ln.startswith(("Bisecting", "    Tidy"))))
    git("bisect", "reset", cwd=pkg)
    forces.write_text(forces.read_text().replace("d = x[:, :, np.newaxis] - x[:, np.newaxis, :]",
                                                 "d = x[:, np.newaxis, :] - x[:, :, np.newaxis]"))
    msg = ("Fix sign of pairwise differences in accel_numpy\n\n"
           "d[k, i, j] must be x_j - x_i so that gravity attracts. The 'tidy' commit swapped the operands,\n"
           "which turned attraction into repulsion; test_numpy_versions_match_loop and\n"
           "test_two_bodies_inverse_square caught it, and git bisect located the commit.")
    git("commit", "-q", "-am", msg, cwd=pkg)
    print("$ git log --oneline -1\n" + git("log", "--oneline", "-1", cwd=pkg))
    print(subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"], cwd=pkg,
                         capture_output=True, text=True).stdout.strip().splitlines()[-1])
```
```text
$ git log --oneline
7533852 README: state the units
768daac README: point to the notes
08c777b Tidy the broadcast in accel_numpy
b3a6c27 Import the nbody package
$ git bisect start HEAD HEAD~3
Bisecting: 0 revisions left to test after this (roughly 1 step)
[768daac57d201aacde00ec7c43894757396ed3f4] README: point to the notes
$ git bisect run python -m pytest -q -x
Bisecting: 0 revisions left to test after this (roughly 0 steps)
08c777b2e1c51036318b02214a40de409d42eca4 is the first bad commit
    Tidy the broadcast in accel_numpy
$ git log --oneline -1
9ba0082 Fix sign of pairwise differences in accel_numpy
14 passed in 1.07s
```

`git bisect start BAD GOOD` marks the endpoints; `git bisect run CMD` checks out midpoints, runs the command, and treats exit code 0 as good and 1 to 127 (except 125, "skip") as bad, halving the range each time; with three candidate commits it needed two test runs. Bisection needs a test that fails on the bug, which is why the swap was "silent" only to a reader: the test comparing `accel_numpy` with the loop failed at once. The fix's message says what broke, why, and what caught it; its first line is what `git log --oneline` shows a year later.

## C. Challenge

**C1.**

```python
n = 4000
p4, _, m4 = uniform_ball(n, np.random.default_rng(4))
full = accel_numpy(p4, m4, 0.05)
work4 = ForceWorkspace(n)
t_full = bench({"w": lambda: work4(p4, m4, 0.05)}, repeat=3)["w"]
print(f"N = {n}: workspace (all pairs) {t_full * 1e3:.0f} ms")
for rc in (0.1, 0.2, 0.4, 0.8, 1.6):
    approx = accel_short_range(p4, m4, rc, eps=0.05, method="tree")
    err = np.median(np.linalg.norm(approx - full, axis=1) / np.linalg.norm(full, axis=1))
    t_rc = bench({"t": lambda rc=rc: accel_short_range(p4, m4, rc, eps=0.05)}, repeat=3)["t"]
    print(f"rc = {rc:3.1f}: median relative error {err:8.2e}, time {t_rc * 1e3:7.1f} ms ({t_rc / t_full:5.2f} of the workspace)")
```
```text
N = 4000: workspace (all pairs) 375 ms
rc = 0.1: median relative error 9.95e-01, time     3.4 ms ( 0.01 of the workspace)
rc = 0.2: median relative error 9.69e-01, time     8.9 ms ( 0.02 of the workspace)
rc = 0.4: median relative error 8.78e-01, time    44.6 ms ( 0.12 of the workspace)
rc = 0.8: median relative error 5.19e-01, time   464.9 ms ( 1.24 of the workspace)
rc = 1.6: median relative error 2.18e-02, time  1491.2 ms ( 3.97 of the workspace)
```

A cutoff is a poor approximation to gravity. Inside a uniform ball, the mass within any sphere around a body that lies entirely inside the ball pulls equally in all directions, so its net force is (up to graininess) zero: the whole force comes from the large-scale asymmetry of the distribution, the far field. A cutoff removes exactly the part that matters, which is why the error is still 88% at $r_c = 0.4$. Even $r_c = 1.6$, most of the ball's diameter, leaves a 2.2% error, and by then the tree finds most of the eight million pairs and costs 4 times the workspace; 1% needs a cutoff covering the whole ball, i.e. all pairs. Cutoffs work for forces that decay faster than the number of neighbours grows (Lennard-Jones, screened Coulomb). Barnes–Hut keeps the far field instead of dropping it: it replaces each distant cell by its total mass at its centre of mass, with an opening angle $\theta$ controlling the error.

**C2.** Script (`c2.py`):

<!-- file: c2.py -->
```python
import math
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

if __name__ == "__main__":
    xs = [float(i) for i in range(10_000)]
    chunks = np.array([1, 4, 16, 64, 256, 1024, 2500])
    times = []
    for c in chunks:
        best = math.inf
        for _ in range(3):
            t0 = time.perf_counter()
            with ProcessPoolExecutor(max_workers=4) as pool:
                out = list(pool.map(math.sqrt, xs, chunksize=int(c)))
            best = min(best, time.perf_counter() - t0)
        assert out[-1] == math.sqrt(xs[-1])
        times.append(best)
    times = np.array(times)
    design = np.column_stack([np.ones(len(chunks)), len(xs) / chunks])
    (t0_fit, a_fit), *_ = np.linalg.lstsq(design, times, rcond=None)
    for c, t, f in zip(chunks, times, design @ [t0_fit, a_fit]):
        print(f"chunksize {c:5d}: {t * 1e3:8.1f} ms   model {f * 1e3:8.1f} ms")
    print(f"T0 = {t0_fit * 1e3:.1f} ms (pool start-up and shutdown), a = {a_fit * 1e6:.0f} us per chunk")
```
```text
chunksize     1:   2019.9 ms   model   2014.9 ms
chunksize     4:    492.0 ms   model    514.4 ms
chunksize    16:    147.3 ms   model    139.3 ms
chunksize    64:     50.6 ms   model     45.5 ms
chunksize   256:     23.3 ms   model     22.0 ms
chunksize  1024:     16.8 ms   model     16.2 ms
chunksize  2500:     17.3 ms   model     15.0 ms
T0 = 14.2 ms (pool start-up and shutdown), a = 200 us per chunk
```

The fit is good over three decades of `chunksize`: a fixed cost $T_0 = 14$ ms of starting and stopping four worker processes, plus $a = 0.2$ ms per chunk (pickle the chunk, send it through a pipe, wake a worker, unpickle, send the results back), next to which the square roots are negligible. With $n/c$ chunks the per-chunk term falls as $1/c$ until it is buried under $T_0$. The model fails at the large-$c$ end, where it predicts 15.0 ms and the measurement is 17.3 ms: with `chunksize = 2500` there are four chunks for four workers, so the last chunk to arrive sets the time, and with real work per item, large chunks would leave workers idle while the slowest finishes. A common choice, about `n / (4 * workers)`, sits between the two failure modes.

**C3.**

```python
import numba


@numba.njit(parallel=True)
def accel_numba_parallel(pos, mass, eps=0.0):
    """Eq. (3.1): each i is independent, so prange splits the outer loop across threads."""
    n = pos.shape[0]
    acc = np.zeros((n, 3))
    for i in numba.prange(n):
        ax = 0.0
        ay = 0.0
        az = 0.0
        for j in range(n):
            if j != i:
                dx = pos[j, 0] - pos[i, 0]
                dy = pos[j, 1] - pos[i, 1]
                dz = pos[j, 2] - pos[i, 2]
                r2 = dx * dx + dy * dy + dz * dz + eps * eps
                w = mass[j] / (r2 * math.sqrt(r2))
                ax += w * dx
                ay += w * dy
                az += w * dz
        acc[i, 0] = ax
        acc[i, 1] = ay
        acc[i, 2] = az
    return acc


for n in (200, 2000):
    p, _, m = uniform_ball(n, np.random.default_rng(n))
    ref = np.array(accel_loop(p.tolist(), m.tolist(), 0.05)) if n == 200 else accel_numba(p, m, 0.05)
    assert np.max(np.abs(accel_numba_parallel(p, m, 0.05) - ref)) < 1e-12 * np.max(np.abs(ref))
    wk = ForceWorkspace(n)
    reps = 20 if n == 200 else 2
    row = bench({"serial numba": lambda: accel_numba(p, m, 0.05), "workspace": lambda: wk(p, m, 0.05)}, number=reps)
    for threads in (1, 2, 4):
        numba.set_num_threads(threads)
        row[f"prange x{threads}"] = bench({"x": lambda: accel_numba_parallel(p, m, 0.05)}, number=reps)["x"]
    print(f"N = {n}: " + ", ".join(f"{k} {v * 1e3:.2f} ms" for k, v in row.items()))
```
```text
N = 200: serial numba 0.16 ms, workspace 0.46 ms, prange x1 0.16 ms, prange x2 0.12 ms, prange x4 0.06 ms
N = 2000: serial numba 15.67 ms, workspace 103.36 ms, prange x1 15.69 ms, prange x2 8.07 ms, prange x4 4.93 ms
```

At $N = 2000$ each outer iteration is 2000 inner iterations of real work, so the threads stay busy: 4 threads gave 3.2× over one, and with them the compiled loop is 21 times faster than the workspace (whose `(3, N, N)` arrays, 96 MB each, no longer fit in cache). At $N = 200$ a whole evaluation takes 0.16 ms, so waking and joining the threads is a visible fraction of it: 4 threads gave 2.7×. The serial `accel_numba` and `prange x1` are the same loop, with the same time. Parallelize the loops that take long enough, and measure where that line is.
