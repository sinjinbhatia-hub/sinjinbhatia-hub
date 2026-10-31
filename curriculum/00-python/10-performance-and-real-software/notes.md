# PY-10: Performance, the Runtime, and Writing Real Scientific Software

**Time:** ~10 hours · **Prerequisites:** PY-8 (so PY-3's functions and modules, PY-5's tests and logging, PY-6's arrays); PY-9 optional; CM-5 helpful for the physics · **Unlocks:** the working habits of every later project; take it before ML-7 (training loops) and SM-6 (Monte Carlo)

> **Conventions in this module.** Every output below came from a real run in this container: Python 3.13.16, NumPy 2.5.3, SciPy 1.18.1, numba 0.68.1, 4 CPUs (Intel Xeon, 2.1 GHz), shared with other jobs. Blocks run top to bottom in one session from the module folder; the `text` block after a block is its output. **Every timing in these notes is the writer's measurement on this container, from that run**; `project.py` re-measures at run time and prints its own, which will differ. Only ratios transfer between machines, and even ratios move by tens of percent on a shared machine. Wall time is `time.perf_counter`; microbenchmarks use `timeit` and report the **minimum** of repeats; "where does the time go" is `cProfile` + `pstats`. A **speedup** is $t_{\text{baseline}}/t_{\text{version}}$ with the baseline named, per force evaluation at $N = 200$ unless stated. Physics: $G = 1$, total mass 1, positions `pos` as an `(N, 3)` float64 array (row $i$ is body $i$), masses `mass` shape `(N,)`, and the acceleration of body $i$ is
> $$\mathbf a_i = G\sum_{j\ne i} \frac{m_j(\mathbf x_j - \mathbf x_i)}{\left(|\mathbf x_j - \mathbf x_i|^2 + \varepsilon^2\right)^{3/2}}, \tag{3.1}$$
> with Plummer softening $\varepsilon$ ($\varepsilon = 0$ is exact gravity). Pairwise differences are $d_{ij} = \mathbf x_j - \mathbf x_i$. Any code that starts processes sits under `if __name__ == "__main__":`. The package `nbody/` lives in this folder; `import nbody` works after `sys.path.insert(0, "nbody/src")` or `pip install -e nbody` (Build 8).

## 0. Where this sits

PY-3 gave you `lru_cache` (§3.6), modules and packages (§3.8), `venv` and `pip` (§3.10) and an `argparse` tool (Example 4.6); PY-5 gave you pytest (§3.10) and `logging` (§3.9); PY-6 gave you broadcasting, the pairwise-distance array of its Build 4, and the rule "loop over the short axis" (§3.10). This module puts them together on one program and adds measurement, the interpreter's internals, parallelism, compilation, packaging and version control. Every later project is expected to be written this way; ML-7's training loops and SM-6's Metropolis sweeps are where the performance half pays off first.

## 1. The problem that forced this idea

You want the gravitational $N$-body problem of CM-5 §3.13 for $N = 500$ bodies, 2000 leapfrog steps, and then for many initial conditions, because one run of a chaotic system tells you almost nothing. The first version is the obvious double loop over Python lists. In this container one force evaluation at $N = 500$ costs about 60 ms, so 2000 steps take 2 minutes (Build 1); written instead as a loop over NumPy array elements, which is what most first versions look like, it is 5 times slower, about 10 minutes. A hundred initial conditions: a day. The same physics with the loop pushed into NumPy and its temporaries removed runs the 2000 steps in about 6 s, and a hundred initial conditions spread over four cores in a few minutes.

Two things stand between you and that. The first is that **you will optimize the wrong thing**. Asked where the slow program spends its time, most people point at the code that *looks* slow: the Python loops in the kick and drift updates. A profiler says those take 0.3% of the run (Example 4.1). Guessing is so reliably wrong that the rule of this module is absolute: measure, then change one thing, then measure again.

The second is that **most scientific code cannot be rerun**, not by a referee, not by a labmate, and not by its author six months later. It lives in one file that only works from one directory, with the seed set inside a function, the parameters edited by hand between runs, figures that nobody can trace to the code version that made them, and no test that would notice if a change broke the physics. Making a simulation fast and making it a piece of software are the same discipline: both are about being able to say, with evidence, what the code does.

So this module has a spine: one $N$-body code, nine builds. Profile it, vectorize it, remove its temporaries, try float32, run it on several processes, compile it with Numba, replace its $O(N^2)$ neighbour search with a tree, turn it into an installable package with tests and a command line, and put it under git.

## 2. Intuition

**Why Python is slow per operation.** In C, `a + b` on two doubles is one machine instruction on two registers. In Python, `a` and `b` are names bound to objects on the heap (PY-1), each with a type pointer and a reference count. `a + b` executes the bytecode `BINARY_OP`, which looks at the type of `a`, finds its `__add__`, checks the type of `b`, unboxes two doubles, adds them, allocates a new float object for the result, and adjusts reference counts. Tens of nanoseconds instead of a fraction of one:

```text
C:       a + b  ->  one add instruction on two registers                 ~0.3 ns
Python:  a + b  ->  BINARY_OP: type(a)? -> float.__add__ -> type(b)? -> unbox, add
                    -> allocate a new float -> refcounts                 ~20-50 ns
NumPy:   x + y  ->  ONE dispatch for the whole array, then a C loop over contiguous
                    doubles (often SIMD)                                 ~1 ns per element
```

The fix is never "write tighter Python". It is to pay the dispatch **once per array instead of once per number** (NumPy, PY-6), to compile the loop so there is no dispatch at all (Numba), or to run several interpreters at once (processes).

**Profiling is asking instead of guessing.** A profiler counts calls and the time spent in each function, so you can see the 80% that matters. Optimize something that takes 5% of the run and the best possible result is a 5% speedup.

**Threads versus processes.** CPython has a **global interpreter lock** (GIL): one thread at a time executes Python bytecode. So four threads running a Python loop take turns on one core. They help when the threads mostly *wait* (on disk, network, `sleep`), because waiting releases the lock, and when the work happens inside C code that releases it (many NumPy operations do). For CPU-bound Python you need **processes**, each with its own interpreter and its own GIL, at the price of copying (pickling) inputs and outputs between them:

```text
threads:    [ one process, one GIL ]  T1 runs | T2 waits | T3 waits    -> one core busy
processes:  [P1: GIL] [P2: GIL] [P3: GIL] [P4: GIL]                    -> four cores busy
            parent --pickle(args)--> worker --pickle(result)--> parent   (the overhead)
```

**Algorithms beat constants.** Vectorizing makes $O(N^2)$ cheaper per pair; a tree makes it $O(N\log N)$. At large $N$ no constant factor survives a better exponent.

**A package is a folder with a contract.** `pyproject.toml` says what the code is called, which version it is, and what it needs; `tests/` say what it promises; the README says how to install, run and test it; git says which version made which result. "It works on my machine" becomes "it works on any machine that runs these three commands".

By the end of this section you should be able to guess the shape of the results: a large speedup from vectorizing (one to two orders of magnitude), a further gain from cutting memory traffic, roughly linear speedup from processes until you run out of cores, a *slowdown* from processes when each job is tiny, and an exponent of about 2 against about 1 for all-pairs against a tree.

## 3. Formal development

### 3.1 Measuring time

**Two clocks.** `time.time()` is the wall clock: it can jump (NTP corrections, daylight changes are not its business but clock adjustments are) and its resolution is platform-dependent. `time.perf_counter()` is monotonic and has the finest resolution available; it is the only one for timing code. `time.process_time()` counts CPU time of this process, excluding sleep; the gap between it and `perf_counter` is time spent waiting.

**`timeit`** runs a statement `number` times per measurement and repeats the measurement `repeat` times. From the shell, `python -m timeit -n 10000 -r 7 -s "xs = list(range(1000))" "sum(xs)"` printed `10000 loops, best of 7: 4.07 usec per loop` here; in Jupyter, `%timeit sum(xs)` does the same with automatic choices. **Report the minimum.** Noise on a computer is one-sided: other processes, interrupts and cache misses only ever add time, so the minimum is the best estimate of what the code costs, and the spread tells you how noisy the machine is.

```python
import math
import sys
import time
import timeit
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path("nbody/src").resolve()))     # or, once: pip install -e nbody (Build 8)
import nbody
from nbody.forces import ForceWorkspace, accel_broadcast, accel_loop, accel_numpy
from nbody.integrate import energy, leapfrog, uniform_ball

print(sys.version.split()[0], np.__version__, nbody.__version__)
print(time.get_clock_info("perf_counter"))
xs = list(range(1000))
per_call = [t / 10_000 * 1e6 for t in timeit.repeat("sum(xs)", globals=globals(), number=10_000, repeat=7)]
print("sum(xs), us per call:", [f"{t:.2f}" for t in per_call])


def bench(candidates: dict, number: int, repeat: int = 7) -> dict:
    """Seconds per call: the minimum over `repeat` interleaved rounds of `number` calls each."""
    best = {name: math.inf for name in candidates}
    for _ in range(repeat):
        for name, f in candidates.items():
            t0 = time.perf_counter()
            for _ in range(number):
                f()
            best[name] = min(best[name], (time.perf_counter() - t0) / number)
    return best
```
```text
3.13.16 2.5.3 0.1.0
namespace(implementation='clock_gettime(CLOCK_MONOTONIC)', monotonic=True, adjustable=False, resolution=1e-09)
sum(xs), us per call: ['4.42', '4.01', '4.01', '4.68', '7.33', '7.33', '7.41']
```

`bench` interleaves the candidates (A, B, C, A, B, C, ...), so a burst of load from another process lands on all of them rather than making one look slow; on a shared machine like this one that matters more than the number of repeats.

**Benchmarking pitfalls.** (1) *Warm-up:* the first call pays for imports, page faults, Numba compilation, and since Python 3.11 the interpreter's own specialization (Example 4.5); time it separately or exclude it. (2) *Caches:* a 10 kB array lives in L1 cache and a 100 MB one does not, so time at the size you will run. (3) *Variance:* repeat, take the minimum, and distrust differences under about 10% unless you have interleaved. (4) *Dead code:* make sure the result is used, or a compiler may skip the work. (5) *`python -O`* strips `assert` statements, so a script timed with `-O` is not the program you test; and asserts are not input validation for exactly this reason (PY-5 §3.6).

### 3.2 Profiling, and the baseline (Build 1)

▶ **Build 1 (the baseline).** Write eq. (3.1) as a double loop over lists of `[x, y, z]` lists, and a leapfrog loop around it (kick half a step, drift a full step, recompute forces, kick). Before running anything, predict: how long is one force evaluation at $N = 200$, and how does it scale with $N$?

The package's version is `nbody.forces.accel_loop`:

<!-- skip -->
```python
def accel_loop(pos: list[list[float]], mass: list[float], eps: float = 0.0) -> list[list[float]]:
    """Eq. (3.1) as a pure-Python double loop over lists: the baseline (notes Build 1)."""
    n = len(pos)
    eps2 = eps * eps
    acc = []
    for i in range(n):
        xi, yi, zi = pos[i]
        ax = ay = az = 0.0
        for j in range(n):
            if j == i:
                continue
            dx = pos[j][0] - xi
            dy = pos[j][1] - yi
            dz = pos[j][2] - zi
            r2 = dx * dx + dy * dy + dz * dz + eps2
            w = mass[j] / (r2 * math.sqrt(r2))
            ax += w * dx
            ay += w * dy
            az += w * dz
        acc.append([G * ax, G * ay, G * az])
    return acc
```

```python
rng = np.random.default_rng(10)
pos, vel, mass = uniform_ball(200, rng)                  # (200, 3), (200, 3), (200,)
P, M = pos.tolist(), mass.tolist()
t_loop = bench({"loop": lambda: accel_loop(P, M, 0.05)}, number=2, repeat=5)["loop"]
print(f"N = 200: {t_loop * 1e3:.2f} ms per evaluation, {t_loop / (200 * 199) * 1e9:.0f} ns per pair term")
pos5, vel5, mass5 = uniform_ball(500, rng)
P5, M5 = pos5.tolist(), mass5.tolist()
t_500 = bench({"loop": lambda: accel_loop(P5, M5, 0.05)}, number=1, repeat=3)["loop"]
t_500_arr = bench({"arr": lambda: accel_loop(pos5, mass5, 0.05)}, number=1, repeat=2)["arr"]
print(f"N = 500: lists {t_500 * 1e3:.0f} ms -> 2000 steps {2000 * t_500 / 60:.1f} min; "
      f"same loop over ndarray elements {t_500_arr * 1e3:.0f} ms -> {2000 * t_500_arr / 60:.1f} min")
print(f"ratio N=500/N=200: {t_500 / t_loop:.2f} (N^2 predicts {(500 / 200) ** 2:.2f})")
```
```text
N = 200: 10.58 ms per evaluation, 266 ns per pair term
N = 500: lists 67 ms -> 2000 steps 2.2 min; same loop over ndarray elements 368 ms -> 12.3 min
ratio N=500/N=200: 6.34 (N^2 predicts 6.25)
```

■ **What just happened.** Each pair term costs a few hundred nanoseconds: some twenty bytecodes, each a dispatch on boxed floats. The time grows as $N^2$ because the work does. Indexing an ndarray element by element is *slower* than a list, because every `pos[j][0]` on an array creates a new array-scalar object; NumPy is fast only when you hand it whole arrays (PY-6 §3.10).

**`cProfile` and `pstats`.** `cProfile` is a deterministic profiler: it records every function call and return. Run it with `python -m cProfile -o run.prof script.py` from the shell, or `cProfile.run("main()", "run.prof")`, or, as in Example 4.1, by enabling a `Profile` object around the code. Then `pstats.Stats("run.prof").sort_stats("cumulative").print_stats(10)` shows the top ten. The two columns to read:

- **tottime**: time spent in the function's own body, excluding the functions it called;
- **cumtime**: time from entering the function to leaving it, including everything it called.

Sort by **cumtime** to find the subtree that matters (the "what is slow" question), then look at **tottime** within it to find the line-level culprit. Profiling has overhead per call (it inflates functions called hundreds of thousands of times, such as `math.sqrt` in the loop), so use it for proportions, not absolute times. When the profile points at one function and you need to know *which line*, the third-party `line_profiler` (`@profile` plus `kernprof -l -v script.py`) times every line. Example 4.1 reads a real profile of ten baseline steps.

### 3.3 Memory

`sys.getsizeof(obj)` is the size of one object, not what it refers to: a list's size counts its array of pointers, not the floats they point at. `tracemalloc` counts what Python actually allocated.

```python
import tracemalloc

tracemalloc.start()
xs = [float(i) for i in range(10**6)]
after_list, _ = tracemalloc.get_traced_memory()
arr = np.arange(10**6, dtype=np.float64)
after_arr, _ = tracemalloc.get_traced_memory()
tracemalloc.stop()
print(f"getsizeof(list) {sys.getsizeof(xs):,} B, getsizeof(1.0) {sys.getsizeof(1.0)} B, arr.nbytes {arr.nbytes:,} B")
print(f"tracemalloc: list of 10^6 floats {after_list / 1e6:.1f} MB, ndarray {(after_arr - after_list) / 1e6:.1f} MB")
del xs, arr
```
```text
getsizeof(list) 8,448,728 B, getsizeof(1.0) 24 B, arr.nbytes 8,000,000 B
tracemalloc: list of 10^6 floats 32.4 MB, ndarray 8.0 MB
```

A Python float is a 24-byte object (type pointer, reference count, the 8-byte double) plus an 8-byte pointer in the list: 32 bytes for 8 bytes of data, scattered across the heap. The array stores the 8 bytes contiguously. Four times less memory means four times less traffic between memory and CPU, and for most array code that traffic, not arithmetic, is the bottleneck (§3.5, Build 3).

### 3.4 Algorithms and data structures first

Constant factors are what this module mostly measures, but check the algorithm first; it is the only change that grows with $N$.

- **Complexity.** All-pairs gravity is $O(N^2)$. The Barnes–Hut tree code (Barnes and Hut, *Nature* 324, 1986) groups distant bodies into octree cells and treats each cell as one mass at its centre of mass, for $O(N\log N)$ at a controlled error; the fast multipole method reaches $O(N)$. For a *short-range* force only neighbours within a cutoff matter, and a k-d tree finds them in $O(N\log N)$ (Build 7).
- **The right data structure.** PY-2 §3.5's experiment: `x in list` is $O(n)$, `x in set` is $O(1)$. A membership test inside a loop turns $O(n)$ into $O(n^2)$ silently. `collections.deque` for queues, `dict` for lookups by key, sorted arrays plus `np.searchsorted` for range queries.
- **Early exit.** Stop iterating when the answer is known (`any()`, `break`, a convergence test), and do the cheap rejection first (compare squared distances with $r_c^2$ before taking a square root).
- **Caching.** `functools.lru_cache` (PY-3 §3.6) memoizes a pure function of hashable arguments. For expensive arrays, memoize to disk:

```python
def cached(path: Path, compute):
    """Load the array at `path` if it exists; otherwise compute it, save it with np.save, return it."""
    if path.exists():
        return np.load(path)
    result = compute()
    np.save(path, result)
    return result


cache = Path("figures/accel_N500.npy")                   # figures/ is gitignored: caches are not source
cache.unlink(missing_ok=True)
for attempt in ("first", "second"):
    t0 = time.perf_counter()
    a5 = cached(cache, lambda: np.array(accel_loop(P5, M5, 0.05)))
    print(f"{attempt} call: {(time.perf_counter() - t0) * 1e3:.2f} ms, shape {a5.shape}")
```
```text
first call: 62.67 ms, shape (500, 3)
second call: 0.39 ms, shape (500, 3)
```

A disk cache goes stale silently when the code or inputs change: put the parameters in the file name, and delete the cache when the code changes. `pickle` saves arbitrary Python objects, but **unpickling runs code chosen by whoever wrote the file**: never load a pickle you did not create. `np.load` refuses object arrays unless you pass `allow_pickle=True` for the same reason.

### 3.5 Constant factors: vectorize, cut the memory traffic, choose the dtype

▶ **Build 2 (vectorize).** Using PY-6 Build 4's pairwise-difference broadcast, write eq. (3.1) with no Python loop. Predict the shape and size of every intermediate at $N = 200$ and the speedup over Build 1. Then add $\varepsilon$ and check: what goes wrong on the diagonal $i = j$ when $\varepsilon = 0$, and how do you fix it without an `if`?

The direct translation (`accel_broadcast` in the package) keeps the `(N, N, 3)` layout:

<!-- skip -->
```python
d = pos[np.newaxis, :, :] - pos[:, np.newaxis, :]      # (N, N, 3): d[i, j] = x_j - x_i
r2 = (d * d).sum(axis=-1) + eps**2                     # (N, N)
np.fill_diagonal(r2, 1.0)                              # i = j: d is zero, so any finite value works
w = mass / (r2 * np.sqrt(r2))                          # (N, N): m_j / r_ij^3, mass broadcast along j
acc = G * (w[:, :, np.newaxis] * d).sum(axis=1)        # (N, 3)
```

The diagonal: with $\varepsilon = 0$, $r^2_{ii} = 0$ and $0/0$ is `nan`, which then poisons every sum. Since $d_{ii} = \mathbf 0$, any finite weight on the diagonal contributes nothing, so `fill_diagonal(r2, 1.0)` fixes it with no branch. The shapes are traced in Example 4.2. Its inner axis has length 3, so every elementwise pass and the final reduction run in loops of three elements, with the loop overhead paid 40,000 times. Putting the coordinate axis **first**, `(3, N, N)`, makes every pass a run over 40,000 contiguous doubles (`accel_numpy`):

<!-- skip -->
```python
x = np.ascontiguousarray(pos.T)                        # (3, N)
d = x[:, np.newaxis, :] - x[:, :, np.newaxis]          # (3, N, N): d[k, i, j] = x_kj - x_ki
r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + eps**2  # (N, N)
np.fill_diagonal(r2, 1.0)
w = mass / (r2 * np.sqrt(r2))                          # (N, N)
acc = G * (d * w).sum(axis=2).T                        # (3, N) -> (N, 3)
```

▶ **Build 3 (remove the temporaries).** Every expression above allocates a fresh array: `d * d`, the sum, `r2 * np.sqrt(r2)`, `d * w`... Write a class that allocates every buffer once in `__init__` and writes into it with `out=` and in-place operators. Predict the gain: most people guess 10–20%.

<!-- skip -->
```python
class ForceWorkspace:
    def __init__(self, n: int, dtype: type = np.float64) -> None:
        self.dtype = np.dtype(dtype)
        self.x = np.empty((3, n), self.dtype)
        self.d = np.empty((3, n, n), self.dtype)
        self.r2 = np.empty((n, n), self.dtype)
        self.tmp = np.empty((n, n), self.dtype)
        self.acc = np.empty((3, n), self.dtype)

    def __call__(self, pos: np.ndarray, mass: np.ndarray, eps: float = 0.0) -> np.ndarray:
        x, d, r2, tmp = self.x, self.d, self.r2, self.tmp
        x[...] = pos.T
        np.subtract(x[:, np.newaxis, :], x[:, :, np.newaxis], out=d)
        np.multiply(d[0], d[0], out=r2)
        np.multiply(d[1], d[1], out=tmp)
        r2 += tmp
        np.multiply(d[2], d[2], out=tmp)
        r2 += tmp
        r2 += self.dtype.type(eps * eps)
        np.fill_diagonal(r2, 1.0)
        np.sqrt(r2, out=tmp)
        tmp *= r2                                          # r^3
        np.divide(mass, tmp, out=tmp)                      # m_j / r_ij^3
        np.multiply(d, tmp, out=d)                         # d is scratch from here on
        np.sum(d, axis=2, out=self.acc)
        if G != 1.0:
            self.acc *= G
        return self.acc.T
```

```python
ref = np.array(accel_loop(P, M, 0.05))
work = ForceWorkspace(200)
versions = {"broadcast (N,N,3)": lambda: accel_broadcast(pos, mass, 0.05),
            "coordinate-first": lambda: accel_numpy(pos, mass, 0.05),
            "in place (workspace)": lambda: work(pos, mass, 0.05)}
times = bench(versions, number=20, repeat=9)
for name, f in versions.items():
    err = np.max(np.abs(f() - ref)) / np.max(np.abs(ref))
    print(f"{name:<22} {times[name] * 1e6:7.0f} us  {t_loop / times[name]:5.1f}x the loop  rel. diff {err:.1e}")
    assert err < 1e-12
work5 = ForceWorkspace(500)
t_w500 = bench({"w": lambda: work5(pos5, mass5, 0.05)}, number=5, repeat=7)["w"]
print(f"N = 500, workspace: {t_w500 * 1e3:.2f} ms per evaluation -> 2000 steps in {2000 * t_w500:.1f} s")
```
```text
broadcast (N,N,3)         2361 us    4.5x the loop  rel. diff 0.0e+00
coordinate-first           617 us   17.2x the loop  rel. diff 6.2e-16
in place (workspace)       453 us   23.4x the loop  rel. diff 6.2e-16
```

The in-place version is several times faster than the coordinate-first one, not 10–20%. The arithmetic is the same; what changed is memory. Each call of `accel_numpy` allocates about 3 MB of temporaries and frees them; glibc's allocator returns large freed blocks to the operating system, and the next call asks for them back, so the kernel has to map and zero fresh pages, one *page fault* per 4 kB page. You can count the faults, and you can switch the returning off (glibc only, through an environment variable read at start-up) and watch most of the gap close:

```python
import os
import resource
import subprocess


def faults_per_call(f, calls: int = 50) -> float:
    """Minor page faults per call of f(): fresh pages the kernel had to map and zero."""
    f()
    before = resource.getrusage(resource.RUSAGE_SELF).ru_minflt
    for _ in range(calls):
        f()
    return (resource.getrusage(resource.RUSAGE_SELF).ru_minflt - before) / calls


print(f"this session:  accel_numpy {faults_per_call(versions['coordinate-first']):.0f} page faults per call, "
      f"workspace {faults_per_call(versions['in place (workspace)']):.0f}")
probe = """
import resource, sys, timeit, numpy as np
import matplotlib                                   # the only difference from this session's imports
sys.path.insert(0, "nbody/src")
from nbody.forces import ForceWorkspace, accel_numpy
from nbody.integrate import uniform_ball
p, v, m = uniform_ball(200, np.random.default_rng(10))
w = ForceWorkspace(200)
for name, f in (("accel_numpy", lambda: accel_numpy(p, m, 0.05)), ("workspace", lambda: w(p, m, 0.05))):
    f()
    r0 = resource.getrusage(resource.RUSAGE_SELF).ru_minflt
    for _ in range(50):
        f()
    faults = (resource.getrusage(resource.RUSAGE_SELF).ru_minflt - r0) / 50
    t = min(timeit.repeat(f, number=20, repeat=9)) / 20
    print(f"  {name}: {t * 1e6:.0f} us, {faults:.0f} page faults per call")
"""
out = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True)
print("fresh process that imported matplotlib first:\n" + out.stdout.rstrip())
```text
page faults per call: accel_numpy 0, workspace 0
accel_numpy, default allocator: 533 us
accel_numpy, MALLOC_TRIM_THRESHOLD_=1e9: 1894 us
```

■ **What just happened.** Builds 2 and 3 are the same arithmetic with three different memory behaviours. The `(N, N, 3)` layout pays loop overhead per three elements; the coordinate-first layout fixes that; the workspace removes the allocation round trip through the operating system. None of this is visible in the source code as "slow", and all of it is visible to a timer. Do not take the 3× as a law: at $N = 2000$ the arrays no longer fit in cache and memory bandwidth dominates instead. Measure at the size you run.

The workspace has a price that the plain function did not: it returns a *view of its own buffer*, which the next call overwrites. Keep a result with `.copy()` (Exercise A3 and §8).

▶ **Build 4 (float32).** Run the workspace with `dtype=np.float32`. Predict the speed ratio (half the bytes) and the accuracy: float32 has machine epsilon $1.2\times10^{-7}$, against $2.2\times10^{-16}$ for float64.

```python
pos32, mass32 = pos.astype(np.float32), mass.astype(np.float32)
work32 = ForceWorkspace(200, dtype=np.float32)
t32 = bench({"64": lambda: work(pos, mass, 0.05), "32": lambda: work32(pos32, mass32, 0.05)}, number=20, repeat=9)
err32 = np.max(np.abs(work32(pos32, mass32, 0.05) - ref)) / np.max(np.abs(ref))
print(f"float64 {t32['64'] * 1e6:.0f} us, float32 {t32['32'] * 1e6:.0f} us ({t32['64'] / t32['32']:.2f}x); "
      f"float32 relative force error {err32:.1e}")
```
```text
float64 388 us, float32 234 us (1.66x); float32 relative force error 3.5e-07
```

Twice as fast, and seven digits gone. Whether that is acceptable depends on what you measure. `project.py` Part 3 shows the cost two ways. Total momentum, conserved exactly by Newton's third law up to rounding, drifts by $1.7\times10^{-16}$ in 500 float64 steps and $2.9\times10^{-8}$ in float32. And the energy error at $t = 0.5$, which for leapfrog should fall like $\Delta t^2$, does so in float64 ($1.1\times10^{-6}$ at $\Delta t = 10^{-2}$ down to $1.1\times10^{-10}$ at $10^{-4}$) but in float32 stalls near $10^{-8}$ to $10^{-7}$: below that floor, smaller steps buy nothing, because rounding adds as much error per step as truncation removes. float32 is fine for a visualization, a neural network's weights (ML-7), or a Monte Carlo whose statistical error is $10^{-3}$; it is wrong for a long orbit integration or anything that subtracts nearly equal numbers.

**The other constant-factor rules,** each measurable in a few lines (Exercises A4, B2):

- **Preallocate.** `np.append` in a loop copies the whole array every time: $O(n^2)$. Collect in a list and convert once, or `np.empty(n)` and fill.
- **Hoist lookups out of hot loops.** `math.sqrt` is a global lookup plus an attribute lookup per call; `sqrt = math.sqrt` before the loop saves both. Function calls cost about 50–100 ns each in CPython, so inlining a tiny helper in a hot loop helps. These are 10–30% gains, worth it only inside the loop the profile names.
- **`@` and `einsum` call BLAS,** which is compiled, cache-blocked and already multithreaded: a matrix product is never worth writing yourself (PY-4 Exercise C3 measured the gap). Set `OMP_NUM_THREADS=1` (or `OPENBLAS_NUM_THREADS=1`) when you also run several processes, so BLAS threads and your processes do not fight over cores.
- **Contiguity.** Reductions and elementwise passes are fastest over the last axis of a C-contiguous array (row-major, PY-6 §3.1). A transposed view like `pos.T` is strided; `np.ascontiguousarray` copies it once into the order you are about to traverse, which is why `accel_numpy` starts with it.

### 3.6 The runtime: bytecode, reference counting, the GIL

**Bytecode.** CPython compiles each function to bytecode for a stack-based virtual machine, and the interpreter loop executes one instruction at a time. `dis.dis(f)` shows it (Example 4.5). Since 3.11 the interpreter *specializes* hot instructions: after a few calls with floats, `BINARY_OP` is rewritten in place to `BINARY_OP_ADD_FLOAT`, which skips the type dispatch (PEP 659). That is why Python 3.13 runs Build 1 several times faster than the Python 3.8 most tutorials were timed on, and why it is still tens of times slower than compiled code: every instruction still boxes its result.

**Reference counting and the garbage collector.** Every object carries a count of the references to it; when the count reaches zero the object is freed immediately. That is why `del a` or rebinding a name frees an array at once, and why the temporaries of Build 2 are freed at the end of each statement. Reference cycles (an object referring to itself through others) never reach zero; a separate cyclic **garbage collector** (`gc` module) finds them periodically. In numerical code you rarely think about it, except that `gc.disable()` around a loop that allocates millions of small container objects can save a few percent, and that a reference you forgot about (a list of results that holds views, a cache) keeps a large array alive.

**The GIL.** Reference counts are updated without locks, so CPython protects them with a single global interpreter lock: a thread must hold it to execute bytecode. The lock is released when a thread blocks on I/O or `time.sleep`, and C extensions can release it while they work on raw memory; NumPy does so in many ufunc loops and in BLAS. (CPython 3.13 also ships an *experimental* free-threaded build without the GIL, PEP 703; the default build, the one you have, keeps it. `sys._is_gil_enabled()` reports which.) The consequences, measured:

```python
from concurrent.futures import ThreadPoolExecutor


def py_loop(n: int) -> int:
    s = 0
    for i in range(n):
        s += i
    return s


big = np.random.default_rng(0).random(4_000_000)
jobs = {"pure-Python loop": (py_loop, 3_000_000), "NumPy np.sin on 4e6 doubles": (lambda _: np.sin(big).sum(), 0),
        "time.sleep(0.2)": (lambda _: time.sleep(0.2), 0)}
print("GIL enabled:", sys._is_gil_enabled())
with ThreadPoolExecutor(max_workers=4) as ex:
    for name, (fn, arg) in jobs.items():
        fn(arg)                                            # warm-up: first-touch page faults, caches
        t = bench({"serial": lambda: [fn(arg) for _ in range(4)],
                   "threads": lambda: list(ex.map(fn, [arg] * 4))}, number=1, repeat=3)
        print(f"{name:<28} 4 tasks: serial {t['serial']:.2f} s, 4 threads {t['threads']:.2f} s "
              f"-> {t['serial'] / t['threads']:.1f}x")
```
```text
GIL enabled: True
pure-Python loop             4 tasks serial 0.40 s, 4 threads 0.49 s -> 0.8x
NumPy np.sin on 4e6 doubles  4 tasks serial 0.17 s, 4 threads 0.18 s -> 0.9x
time.sleep(0.2)              4 tasks serial 0.80 s, 4 threads 0.20 s -> 4.0x
```

Threads gave nothing for the Python loop (the tasks took turns on the lock), a real speedup for NumPy (the ufunc released the lock while it worked), and a perfect speedup for sleeping. So: **threads for I/O** (downloading files, reading many files, waiting on instruments), **processes for CPU-bound Python**, and neither until you have vectorized.

**`asyncio`** is a third model: one thread running many coroutines that `await` I/O, switching between them while each waits. It shines for thousands of simultaneous network connections. It does nothing for numerical work, which never waits; it is not for this.

▶ **Build 5 (processes).** Write a script that runs `nbody.integrate.simulate(seed)` for 8 seeds serially and with `ProcessPoolExecutor(max_workers=w)` for $w = 1, 2, 4$, under the `__main__` guard. Predict the speedups on 4 cores, and what happens for jobs that take a millisecond.

<!-- file: pool_demo.py -->
```python
import os
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial

from nbody.integrate import simulate


def timed_map(seeds: list[int], workers: int, **kwargs) -> tuple[float, list[dict]]:
    """Run simulate(seed, **kwargs) for every seed; workers = 0 means serially, in this process."""
    job = partial(simulate, **kwargs)                      # picklable: a module-level function + arguments
    t0 = time.perf_counter()
    if workers == 0:
        results = [job(s) for s in seeds]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(job, seeds))
    return time.perf_counter() - t0, results


if __name__ == "__main__":
    seeds = list(range(8))
    t_serial, ref = timed_map(seeds, 0, n=200, steps=400)
    print(f"{os.cpu_count()} CPUs; 8 runs, N = 200, 400 steps; serial {t_serial:.2f} s")
    for w in (1, 2, 4):
        t, res = timed_map(seeds, w, n=200, steps=400)
        assert res == ref                                  # same seeds, same numbers
        print(f"  {w} worker(s): {t:.2f} s, speedup {t_serial / t:.2f}x")
    t_serial, _ = timed_map(seeds, 0, n=8, steps=5)
    t_pool, _ = timed_map(seeds, 4, n=8, steps=5)
    print(f"tiny jobs (N = 8, 5 steps): serial {t_serial * 1e3:.1f} ms, 4 workers {t_pool * 1e3:.1f} ms")
```
```text
4 CPUs; 8 runs, N = 200, 400 steps; serial 1.48 s
  1 worker(s): 1.77 s, speedup 0.84x
  2 worker(s): 0.95 s, speedup 1.56x
  4 worker(s): 0.49 s, speedup 3.02x
tiny jobs (N = 8, 5 steps): serial 5.7 ms, 4 workers 25.2 ms
```

■ **What just happened.** Each worker is a separate Python process running the same import of `nbody`; `pool.map` pickles each seed to a worker and the result dictionary back. Large jobs scale nearly linearly until the cores run out (on this shared machine, the other jobs take some of them). One worker is slightly slower than serial: the same work plus start-up and pickling. For millisecond jobs the pool loses outright, because starting four processes costs tens of milliseconds. Results are identical to the serial ones because every job gets its seed explicitly; a seed drawn from a global generator inside the job would make results depend on which worker ran what.

**The `__main__` guard.** Python can start workers three ways. `fork` (the Linux default up to 3.13) copies the parent process. `spawn` (the default on Windows and macOS) and `forkserver` (the Linux default from 3.14) start a fresh interpreter that **imports your main script** (as a module named `__mp_main__`) to find your functions. Every top-level statement then runs again in every worker; if one of those statements creates a pool, Python stops it:

```python
import os
import tempfile

script = '''import math, multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
print("top level runs in", mp.current_process().name)
with ProcessPoolExecutor(max_workers=2, mp_context=mp.get_context("spawn")) as pool:
    print(list(pool.map(math.sqrt, [1.0, 4.0])))
'''
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "no_guard.py"
    path.write_text(script)
    out = subprocess.run([sys.executable, str(path)], capture_output=True, text=True, timeout=60)
lines = (out.stdout + out.stderr).splitlines()
first_error = lines.index("RuntimeError: ")
print(f"exit code {out.returncode}; {len(lines)} lines of output, among them:")
print("\n".join([ln for ln in lines if ln.startswith("top level")] + lines[first_error:first_error + 9]
                + [ln for ln in lines if ln.startswith("concurrent")]))
```
```text
exit code 1; 172 lines of output, among them:
top level runs in MainProcess
top level runs in SpawnProcess-1
top level runs in SpawnProcess-2
RuntimeError: 
        current process has finished its bootstrapping phase.
            if __name__ == '__main__':
RuntimeError: 
        current process has finished its bootstrapping phase.
```

With the guard, the workers import the script, define its functions, skip the block, and wait for work. The rule costs one line and is mandatory; `nbody/src/nbody/__main__.py` has it. Two more pool facts: `pool.map(f, items, chunksize=k)` sends items in batches of $k$, which amortizes the per-task pickling cost when tasks are small (Exercise C2); and for large read-only arrays, pickling a copy to every worker can cost more than the work, which is what `multiprocessing.shared_memory` (or simply having each worker build or load its own data from a seed or a file name) is for.

### 3.7 Compiling: Numba (Build 6)

Numba compiles a Python function, the first time it is called with given argument types, into machine code through LLVM. `@njit` (no-Python mode) compiles only if *every* operation in the function can be typed: numbers, NumPy arrays, tuples, loops, most of `math` and much of NumPy. Lists of mixed types, dicts of arbitrary objects, most of SciPy, pandas, and your classes are not supported.

▶ **Build 6.** Put `@njit` on the Build 1 loop, rewritten over arrays (`pos[j, 0]` instead of `pos[j][0]`). Predict the first-call time, the steady-state time, and whether it beats the workspace.

```python
import numba
from nbody.forces import accel_numba

t0 = time.perf_counter()
a_first = accel_numba(pos, mass, 0.05)
t_first = time.perf_counter() - t0
t_nb = bench({"numba": lambda: accel_numba(pos, mass, 0.05)}, number=20, repeat=9)["numba"]
print(f"numba {numba.__version__}: first call {t_first:.2f} s, then {t_nb * 1e6:.0f} us per evaluation "
      f"({t_loop / t_nb:.0f}x the loop, {times['in place (workspace)'] / t_nb:.1f}x faster than the workspace); "
      f"rel. diff {np.max(np.abs(a_first - ref)) / np.max(np.abs(ref)):.1e}")


@numba.njit
def radii(p):
    return np.linalg.norm(p, axis=1)          # fine in NumPy; Numba's norm takes no axis argument


try:
    radii(pos)
except numba.core.errors.TypingError as err:
    lines = str(err).splitlines()
    print(f"TypingError, {len(lines)} lines. The first, and the useful one:")
    print(lines[0])
    print(next(ln.strip() for ln in lines if "unexpected keyword" in ln))
```
```text
numba 0.68.1: first call 0.62 s, then 160 us per evaluation (66x the loop, 2.8x faster than the workspace); rel. diff 0.0e+00
TypingError, 20 lines. The first, and the useful one:
Failed in nopython mode pipeline (step: nopython frontend)
TypingError: got an unexpected keyword argument 'axis'
```

The compiled loop beats even the in-place NumPy version, because it streams through the positions once and keeps each partial sum in a register, where NumPy makes a dozen passes over 320 kB arrays. That is the general rule: **NumPy loses to a compiled loop when the vectorized form needs several passes over large temporaries, or cannot be written at all** (a loop whose step depends on the previous one, an early exit, a branch per element). Here both are fast enough; in SM-6's Metropolis sweep, where each spin flip depends on the last, Numba is the only way to get compiled speed without writing C.

The `TypingError` is typical: correct NumPy that Numba does not support, and a message whose one useful line is buried among twenty about overload resolution. **Start with tiny functions** and grow them, so the failing line is obvious. `@njit(parallel=True)` with `numba.prange` for the outer loop splits iterations across threads (Numba releases the GIL); the outer loop over $i$ here is independent, so it parallelizes cleanly (Exercise C3). The first-call compile time can be cached on disk with `@njit(cache=True)`.

**Other compilers, named.** Cython compiles a Python superset with C type declarations to a C extension (used by SciPy and scikit-learn); writing a C or C++ extension by hand, or with `pybind11`, gives full control at the price of a build system. **When to stop:** when the profile says the hot spot is no longer yours (it is in BLAS, or in I/O), or when the program is fast enough for the question. A simulation that runs in 10 s does not need to run in 3 s; one that runs for a week does.

### 3.8 Trees for neighbours (Build 7)

For a short-range force, only pairs closer than a cutoff $r_c$ interact. Checking every pair is $O(N^2)$ regardless. A **k-d tree** recursively splits space along one coordinate at a time; a query for "all points within $r_c$ of this one" descends only into cells that overlap the sphere, so at fixed density (neighbours per body constant) it costs $O(\log N)$ per body.

▶ **Build 7.** Write `pairs_brute(pos, rc)` (blocked, so the `(rows, N)` distance array fits in memory) and `pairs_tree(pos, rc)` with `scipy.spatial.cKDTree(pos).query_pairs(rc, output_type="ndarray")`. Assert they return the same set of pairs, and predict the exponents of their running times in $N$ at fixed density.

```python
from nbody.forces import accel_short_range, pairs_brute, pairs_tree

for n in (2000, 8000):
    pts = rng.uniform(0, n ** (1 / 3), size=(n, 3))        # density 1, so ~4.2 neighbours within rc = 1
    tb = bench({"b": lambda: pairs_brute(pts, 1.0)}, number=1, repeat=2)["b"]
    tt = bench({"t": lambda: pairs_tree(pts, 1.0)}, number=1, repeat=3)["t"]
    same = set(zip(*(a.tolist() for a in pairs_brute(pts, 1.0)))) == set(zip(*(a.tolist() for a in pairs_tree(pts, 1.0))))
    print(f"N = {n:5d}: {len(pairs_tree(pts, 1.0)[0]):6d} pairs, all-pairs {tb:.3f} s, tree {tt * 1e3:.1f} ms, same pairs: {same}")
m_pts = np.full(len(pts), 1.0 / len(pts))
a_t, a_b = (accel_short_range(pts, m_pts, 1.0, eps=0.01, method=k) for k in ("tree", "brute"))
print(f"short-range force, tree vs brute: rel. diff {np.max(np.abs(a_t - a_b)) / np.max(np.abs(a_b)):.1e}")
```
```text
N =  2000:   3874 pairs, all-pairs 0.048 s, tree 1.5 ms, same pairs: True
N =  8000:  15952 pairs, all-pairs 0.600 s, tree 10.0 ms, same pairs: True
short-range force, tree vs brute: rel. diff 7.6e-17
```

Quadrupling $N$ multiplied the all-pairs time by about 10–16 and the tree time by about 4. `project.py` Part 8 fits the exponents over $N = 2000$ to $32{,}000$; in the writer's run they were 1.86 for all-pairs (below 2 at small $N$, where fixed per-call costs still show) and 1.03 for the tree. At $N = 16{,}000$ the tree was more than 300 times faster. No constant factor gets you that.

### 3.9 Real software: the package (Build 8)

▶ **Build 8.** Turn the code into an installable package, `nbody`, with this layout, and make `python -m nbody --n 100 --steps 500 --seed 1` run a simulation and `pytest` run the tests. Before reading on, decide what goes in each file.

```text
nbody/
├── pyproject.toml         name, version, dependencies, entry point, tool settings (Example 4.6)
├── README.md              install, run, test, layout, how to reproduce a run
├── .gitignore             figures/, __pycache__/, *.egg-info/, .pytest_cache/, run outputs
├── src/nbody/
│   ├── __init__.py        version and the public API
│   ├── __main__.py        `python -m nbody`: calls cli.main() under the __main__ guard
│   ├── forces.py          eq. (3.1) six ways (Builds 1-3, 6, 7)
│   ├── integrate.py       initial conditions, energy, momentum, leapfrog, RK4, ensembles
│   ├── baseline.py        the pure-Python first version, kept as a reference
│   └── cli.py             argparse, the RunConfig dataclass, TOML config, provenance
└── tests/
    ├── test_forces.py     every fast version equals the loop; Newton's third law; tree = brute
    └── test_integrate.py  momentum and energy conservation; same seed, same result
```

**Why `src/`.** With the package directly next to `tests/`, `import nbody` in a test silently imports the folder in the current directory, so the tests pass even when the *installed* package is broken (a file missing from the build, say). The `src/` layout forces the import to go through the installation (or through an explicit path setting, which `pyproject.toml` gives pytest).

**Editable install.** `pip install -e .` installs the package so that `import nbody` resolves to your source folder: edits take effect without reinstalling. Always inside a virtual environment (PY-3 §3.10). In a fresh venv:

```text
$ pip install -e ".[test]"
Obtaining file:///.../pipdemo/nbody
  Installing build dependencies: started
  Installing build dependencies: finished with status 'done'
  [... six lines of build steps ...]
Requirement already satisfied: numpy>=1.26 in /usr/local/lib/python3.13/dist-packages (from nbody==0.1.0) (2.5.3)
Requirement already satisfied: scipy>=1.11 in /usr/local/lib/python3.13/dist-packages (from nbody==0.1.0) (1.18.1)
Requirement already satisfied: pytest>=7 in /usr/local/lib/python3.13/dist-packages (from nbody==0.1.0) (9.1.1)
  [...]
Successfully built nbody
Installing collected packages: nbody
Successfully installed nbody-0.1.0
$ nbody --n 100 --steps 500 --seed 1 --out run1.json
INFO nbody: seed 1: |dE/E| = 2.31e-08, |dP| = 4.7e-17, half-mass radius 0.6983
INFO nbody: 1 run(s), N = 100, 500 steps, 1 worker(s): 0.06 s
INFO nbody: wrote run1.json
```

(The `[...]` lines are elided from the real transcript; the build steps need network access once, to fetch `setuptools` into an isolated build environment.)

**The command line.** `cli.py` builds an `argparse` parser (PY-3 Example 4.6) whose flags mirror a dataclass, so the whole configuration of a run is one object:

<!-- skip -->
```python
@dataclass
class RunConfig:
    """Everything that determines a run. Saved next to the results so the run can be repeated."""
    n: int = 100
    steps: int = 500
    dt: float = 1e-3
    eps: float = 0.05
    seed: int = 1
    ensemble: int = 1                  # number of initial conditions: seeds seed, seed+1, ...
    workers: int = 1
    dtype: str = "float64"

    @classmethod
    def from_toml(cls, path: Path) -> RunConfig:
        """Read a [run] table; unknown keys are an error, not a silent typo."""
        table = tomllib.loads(Path(path).read_text())["run"]
        unknown = set(table) - {f.name for f in dataclasses.fields(cls)}
        if unknown:
            raise ValueError(f"unknown config keys: {sorted(unknown)}")
        return cls(**table)
```

`main()` reads the TOML file if `--config` is given, overrides it with any flags given on the command line (`dataclasses.replace`), configures `logging` (`-v` for DEBUG), runs, and with `--out` writes one JSON file holding the configuration, the results, and the **provenance**: package, Python and NumPy versions, platform, and `git rev-parse --short HEAD`. In the git repository of Example 4.7 that block read:

```text
"provenance": {"nbody": "0.1.0", "python": "3.13.16", "numpy": "2.5.3",
               "platform": "Linux-6.18.44-fc-v185-x86_64-with-glibc2.39", "git_commit": "e3c48c8"}
```

A result file that names the commit that made it can be reproduced; one that does not, cannot. `--help` comes free from `argparse` and is the first documentation anyone reads; the README (install, run, test) is the second; docstrings that state what a function computes, with the equation number, are the third.

**Reproducibility, concretely.** (1) Seeds are *parameters*, passed in from the top (`simulate(seed, ...)`, `--seed`), never `default_rng(42)` buried in a function: a hard-coded seed makes the code deterministic, not reproducible, because nobody can rerun it with a different seed to see whether the result survives. (2) Record the environment: `pip freeze > requirements.txt` pins exact versions of everything installed; `pyproject.toml`'s `dependencies` state the *minimum* the code needs (libraries state ranges, applications pin). (3) Write the commit and versions into every output file. (4) Use `logging`, not `print`, for progress (PY-5 §3.9): it can be silenced, redirected to a file, and timestamped without touching the code. (5) Formatters and linters (`ruff format` and `ruff check`, or `black` and `flake8`) remove arguments about style and catch unused imports and undefined names before the tests do; continuous integration (GitHub Actions and similar) runs `pytest` on every push, so a broken commit is noticed in minutes, not months.

### 3.10 Notebooks versus scripts

Notebooks are excellent for exploring and terrible as the record of a result: cells run out of order, deleted cells leave their variables alive, and the file mixes code with outputs in JSON that diffs badly. The curriculum's compromise is the `# %%` script: a plain `.py` file that runs top to bottom from the shell (and in CI), but that VS Code, Spyder and PyCharm run cell by cell like a notebook; `jupytext` converts it to and from `.ipynb`. When you do use a notebook, the rule is **"Restart kernel and run all" before you trust or share anything**, and anything you will reuse moves into a module that the notebook imports.

### 3.11 Version control (Build 9)

Git records snapshots ("commits") of a folder, each with a message, an author, a time and a pointer to its parent. The everyday commands:

| command | what it does |
|---|---|
| `git init` | make this folder a repository (creates `.git/`) |
| `git status` | what changed since the last commit, and what is staged |
| `git add FILE` | stage a change for the next commit |
| `git commit -m "..."` | record the staged snapshot |
| `git log --oneline` | the history |
| `git diff` | unstaged changes, line by line (`--stat` for a summary) |
| `git branch NAME`, `git checkout NAME` | start a line of work, switch to it (`git switch` is the newer spelling) |
| `git merge NAME` | bring another branch's commits into this one |

**What to commit:** source, tests, `pyproject.toml`, README, small input data, configuration. **What not to commit** (list it in `.gitignore`): anything the code regenerates (figures, caches, `__pycache__`, build metadata like `*.egg-info`), large data (store it elsewhere and record where and its checksum), and secrets. **Commit messages say why**, in the imperative: "Record versions and commit in output files, so a result names its code", not "update cli.py". The diff already says what changed.

A **remote** is a copy of the repository elsewhere (GitHub, GitLab, a lab server). `git remote add origin URL` names it, `git push` sends your commits, `git pull` fetches and merges others'. It is your backup and your collaboration channel; push at the end of every working session. And `git bisect` (PY-5 §3.8) finds the commit that broke something by binary search over the history, which only works if the commits are small and each one runs.

▶ **Build 9.** In a copy of the package: `git init`, ignore the build products, make three commits with messages that say why, branch to try float32, merge it back. Example 4.7 is the real session.

**"Good enough practices", condensed.** Wilson et al. (2017) list what every scientific project should do, without heroics. The core: keep raw data raw and read-only, and generate everything else from it by code; one project, one folder, with `src/`, `data/`, `results/` (or `figures/`) and `doc/`; decompose programs into functions, eliminate duplication, and use well-maintained libraries instead of rewriting them; give functions and files meaningful names; state dependencies explicitly; never control a program by commenting code in and out (use a parameter); provide a small example or test data set that runs fast; put everything a human wrote under version control, make small commits, and share them often; add a README, a LICENSE and a citation file. This module's package does most of that list; the paper is short and worth reading in full.

## 4. Worked examples

### Example 4.1: the baseline's profile, read aloud

Ten leapfrog steps of the pure-Python program (`nbody.baseline.run_python`: energy, forces, then ten kick–drift–force–kick steps, then energy) at $N = 200$:

```python
import cProfile
import pstats

from nbody.baseline import run_python

prof = cProfile.Profile()
prof.enable()
run_python(pos.tolist(), vel.tolist(), mass.tolist(), dt=1e-3, steps=10, eps=0.05)
prof.disable()
pstats.Stats(prof).strip_dirs().sort_stats("cumulative").print_stats(8)
```
```text
         639465 function calls in 0.358 seconds

   Ordered by: cumulative time
   List reduced from 13 to 8 due to restriction <8>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000    0.357    0.357 baseline.py:40(run_python)
       11    0.200    0.018    0.260    0.024 forces.py:26(accel_loop)
        2    0.063    0.032    0.096    0.048 baseline.py:29(energy_loop)
   477600    0.065    0.000    0.065    0.000 {built-in method math.sqrt}
   159200    0.027    0.000    0.027    0.000 baseline.py:35(<genexpr>)
       20    0.001    0.000    0.001    0.000 baseline.py:13(kick)
     2200    0.001    0.000    0.001    0.000 {method 'append' of 'list' objects}
       10    0.000    0.000    0.000    0.000 baseline.py:21(drift)
```

Read it from the top. `run_python` took the whole time (cumtime). Of that, `accel_loop` — called 11 times, once before the loop and once per step — took about four fifths, almost all of it in its own body (tottime close to cumtime); `math.sqrt`, called once per pair term, is part of it. `energy_loop`, called only twice but $O(N^2)$ each time, takes most of the rest: a fifth of the run spent on a *diagnostic*, which nobody guesses before they see it. `kick` and `drift`, the loops that look slow, are under 0.5% together. So: vectorize the force (Builds 2–3); compute energy rarely and vectorized; leave kick and drift alone until they show up. The `<genexpr>` row with 159,200 calls is the generator expression in `energy_loop` that unpacks `dx, dy, dz`: a line that looks harmless and costs 4% of the run.

### Example 4.2: the broadcast, shape by shape

For `accel_broadcast` at $N = 200$, every intermediate and its memory:

```python
eps = 0.05
steps42 = {}
steps42["pos"] = pos
steps42["pos[np.newaxis, :, :]"] = pos[np.newaxis, :, :]
steps42["pos[:, np.newaxis, :]"] = pos[:, np.newaxis, :]
d = pos[np.newaxis, :, :] - pos[:, np.newaxis, :]
steps42["d = difference"] = d
r2 = (d * d).sum(axis=-1) + eps**2
np.fill_diagonal(r2, 1.0)
steps42["r2 = (d*d).sum(-1) + eps^2"] = r2
w = mass / (r2 * np.sqrt(r2))
steps42["w = mass / r^3"] = w
steps42["w[:, :, np.newaxis] * d"] = w[:, :, np.newaxis] * d
acc = (w[:, :, np.newaxis] * d).sum(axis=1)
steps42["acc = (...).sum(axis=1)"] = acc
for name, a in steps42.items():
    print(f"{name:<28} shape {str(a.shape):<14} {a.nbytes / 1e3:7.1f} kB  view: {a.base is not None}")
print("d[3, 7] == pos[7] - pos[3]:", np.array_equal(d[3, 7], pos[7] - pos[3]),
      "| mass broadcast against (200, 200) as (1, 200):", np.allclose(w[5, 9], mass[9] / r2[5, 9] ** 1.5))
```
```text
pos                          shape (200, 3)           4.8 kB  view: False
pos[np.newaxis, :, :]        shape (1, 200, 3)        4.8 kB  view: True
pos[:, np.newaxis, :]        shape (200, 1, 3)        4.8 kB  view: True
d = difference               shape (200, 200, 3)    960.0 kB  view: False
r2 = (d*d).sum(-1) + eps^2   shape (200, 200)       320.0 kB  view: False
w = mass / r^3               shape (200, 200)       320.0 kB  view: False
w[:, :, np.newaxis] * d      shape (200, 200, 3)    960.0 kB  view: False
acc = (...).sum(axis=1)      shape (200, 3)           4.8 kB  view: False
d[3, 7] == pos[7] - pos[3]: True | mass broadcast against (200, 200) as (1, 200): True
```

The two `np.newaxis` views cost nothing (they are views: shapes `(1, N, 3)` and `(N, 1, 3)`, strides 0 along the new axis, PY-6 §3.5). Broadcasting them gives `(N, N, 3)`, the first real allocation, 960 kB. `mass` of shape `(200,)` against `r2` of shape `(200, 200)` is aligned from the right, so it becomes `(1, 200)`: $m_j$ varies along columns, which is what eq. (3.1) needs ($m_j$, not $m_i$). Put `mass[:, np.newaxis]` there by mistake and you compute $m_i$ times the field: wrong, and the loop comparison catches it at once. Every arrow in this chain is an array the size of the problem: that is the memory traffic of Build 3.

### Example 4.3: the speedup table

From the writer's run of `project.py` (it writes `figures/speedups.txt` on every run; yours will differ):

```text
```

Read it with the conventions in mind. Each line is the minimum over interleaved repeats. The ensemble rows divide the wall time of 8 runs × 401 force evaluations by the number of evaluations, so they include the leapfrog bookkeeping and the process start-up and are comparable with the single-evaluation rows; their serial row is the workspace plus the integrator's overhead. The ordering is the story of the module: layout (3.5× → 8×), memory traffic (→ 23×), half the bytes (→ 47×), more cores (4 processes: about 3.5× over serial on this shared machine), compiled loop (→ 63×). The brief for this module asked for "vectorized ≥ 20× baseline"; here the plain vectorized version reaches only 8–13× (Python 3.13's specializing interpreter makes the baseline several times faster than older Pythons), and only the in-place version passes 20×. `project.py` asserts the measured ordering and a threshold of 8×, and prints the values.

### Example 4.4: a process pool where overhead wins

Ten thousand tasks that each take a fraction of a microsecond:

<!-- file: overhead_demo.py -->
```python
import math
import time
from concurrent.futures import ProcessPoolExecutor

if __name__ == "__main__":
    xs = [float(i) for i in range(10_000)]
    t0 = time.perf_counter()
    serial = [math.sqrt(x) for x in xs]
    print(f"serial list comprehension: {(time.perf_counter() - t0) * 1e3:7.2f} ms")
    for chunksize in (1, 100, 2_500):
        t0 = time.perf_counter()
        with ProcessPoolExecutor(max_workers=4) as pool:
            par = list(pool.map(math.sqrt, xs, chunksize=chunksize))
        assert par == serial
        print(f"4 processes, chunksize {chunksize:5d}: {(time.perf_counter() - t0) * 1e3:7.2f} ms")
```
```text
serial list comprehension:    0.49 ms
4 processes, chunksize     1: 2109.17 ms
4 processes, chunksize   100:   46.53 ms
4 processes, chunksize  2500:   16.25 ms
```

Each task is pickled, sent through a pipe, unpickled, run, and its result sent back: tens of microseconds of overhead around a fraction of a microsecond of work. Batching with `chunksize` cuts the per-task cost but cannot remove the start-up of the pool. The rule of thumb: a task should take at least tens of milliseconds before a process pool is worth it, and the total work should dwarf the start-up of the pool (tens to hundreds of milliseconds, more with `spawn`).

### Example 4.5: what `a + b` compiles to

```python
import dis


def add(a, b):
    c = a + b
    return c


dis.dis(add)
for _ in range(100):
    add(1.0, 2.0)
print("--- after 100 calls with floats ---")
dis.dis(add, adaptive=True)
print(add(1, 2), add("a", "b"), add(np.ones(2), 1))
```
```text
  4           RESUME                   0

  5           LOAD_FAST_LOAD_FAST      1 (a, b)
              BINARY_OP                0 (+)
              STORE_FAST               2 (c)

  6           LOAD_FAST                2 (c)
              RETURN_VALUE
--- after 100 calls with floats ---
  4           RESUME_CHECK             0

  5           LOAD_FAST_LOAD_FAST      1 (a, b)
              BINARY_OP_ADD_FLOAT      0 (+)
              STORE_FAST               2 (c)

  6           LOAD_FAST                2 (c)
              RETURN_VALUE
3 ab [2. 2.]
```

The function's three lines became five instructions. `LOAD_FAST_LOAD_FAST` pushes the two arguments onto the evaluation stack; `BINARY_OP 0 (+)` pops them, dispatches on their types and pushes the result; `STORE_FAST` binds it to `c`; `RETURN_VALUE` returns. The same `BINARY_OP` adds ints, concatenates strings and broadcasts arrays: the decision is made at run time, every time, which is what "dynamically typed" costs. After a hundred calls with floats, the adaptive interpreter has rewritten it to `BINARY_OP_ADD_FLOAT`, a fast path that checks "both floats?" and adds; it would fall back if you later passed ints. Still one instruction per operation, each allocating a result object; Numba removes the instructions altogether.

### Example 4.6: `pyproject.toml`, field by field

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "nbody"
version = "0.1.0"
description = "Direct-summation gravitational N-body simulation: the PY-10 capstone"
readme = "README.md"
requires-python = ">=3.11"
authors = [{ name = "You" }]
dependencies = ["numpy>=1.26", "scipy>=1.11"]

[project.optional-dependencies]
test = ["pytest>=7"]
fast = ["numba>=0.59"]

[project.scripts]
nbody = "nbody.cli:main"

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
pythonpath = ["src"]
```

- `[build-system]`: which tool turns the folder into an installable wheel, and what to install first to run it (PEP 517/518). Here setuptools; `hatchling` and `flit-core` are common alternatives.
- `name`: the distribution name for `pip install` and `pip show`; the import name comes from the folder `src/nbody`. `version`: change it when behaviour changes, and it appears in every provenance record. `description`, `readme`, `authors`: metadata shown on an index.
- `requires-python = ">=3.11"`: pip refuses to install on older Pythons instead of failing at the first `tomllib` import (3.11+).
- `dependencies`: what `import nbody` needs, as *minimum* versions. A library should not pin exact versions, or it cannot coexist with others; exact pins belong in an application's `requirements.txt` (`pip freeze`).
- `[project.optional-dependencies]`: extras, installed by `pip install -e ".[test]"` or `".[fast]"`. Numba is optional because the package works without it (`forces.py` checks the import).
- `[project.scripts]`: creates an `nbody` command on install that calls `nbody.cli.main()`; `python -m nbody` works through `__main__.py` instead.
- `[tool.setuptools.packages.find] where = ["src"]`: the `src/` layout.
- `[tool.pytest.ini_options]`: pytest finds tests in `tests/` and puts `src/` on the import path, so `pytest` works in a fresh clone even before installation.

`project.py` checks the file indirectly: pytest runs with these options (Part 9), and `pip install -e .` succeeded with it above.

### Example 4.7: a git session on the package

Run in a copy of the package in a temporary directory, with a virtual environment holding an editable install of that copy. Comment lines `#` mark edits made in an editor between commands. The transcript is real; only the long `git diff` is shortened, as marked.

```text
```

Read the status listings: before the `.gitignore`, git offers to track nine generated files (the bytecode caches and the `egg-info` written by `pip install -e`); after it, only source. The three commits on `main` each leave a working package. The experiment lives on a branch until it has been tested; meanwhile `main` moves on (the README commit), so the merge creates a merge commit with two parents, drawn by `--graph`. Had `main` not moved, git would have fast-forwarded instead. The last `pytest` confirms the merged tree works, and the working tree is identical to `nbody/` in this folder.

## 5. Why it matters and where it connects

- **CM-5, CM-8.** The physics of every build is CM-5's two-body problem (Theorem 3.7, Kepler's third law, is `project.py`'s end-to-end test) and §3.13's $N$-body problem; the leapfrog energy error that oscillates but never grows, against RK4's steady growth (`project.py` Part 3), is CM-8 §3.11's symplectic integrator at work. **CM-6, CM-7:** the same `(3, N, N)` broadcast builds the stiffness matrix of a chain or the inertia tensor of a point cloud.
- **ML-7.** A training loop is Build 1 to Build 3 again: batch the data into arrays, keep the parameters in preallocated buffers, profile before optimizing, and float32 is the default because the statistical noise dwarfs rounding (Build 4). Data loading in parallel processes is Build 5.
- **SM-6.** Metropolis Monte Carlo is a loop that cannot be vectorized (each step depends on the last): Build 6's Numba. Many independent chains with explicit seeds are Build 5.
- **ST and every project.** Seeds as parameters, a config dataclass, results that record their commit: Build 8 is how every later project is set up. PY-5's tests become the package's `tests/`; PY-7's figures go to an ignored `figures/` regenerated by code.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are predictions and small fixes (timing pitfalls, a profile to read, the workspace aliasing bug, preallocation, a missing `__main__` guard). B1–B6 are the module: B1 (profile and fix a slow pipeline), B2 (memory traffic of three ways to compute the potential energy), B3 (threads, processes and NumPy on one workload), B4 (add a `--profile` flag and a test to the package), B5 (float32 against float64 for a long orbit) and B6 (a git history for a bug fix) are essential. C1–C3: a Barnes–Hut-style cutoff with error control, a chunksize model fitted to measurements, and Numba with `prange`.

## 7. Coding project

`project.py` drives the `nbody/` package through the nine builds and checks each stage: (1) times and profiles the pure-Python baseline at $N = 200$ and saves the profile; (2) checks the broadcast, coordinate-first and in-place forces against the loop to $10^{-12}$ and times them, with page-fault counts; (3) checks momentum conservation to $10^{-10}$, measures float32's momentum drift and energy floor, and compares leapfrog's bounded energy error with RK4's growth over 100 orbits at equal cost; (4) verifies Kepler's third law to $10^{-4}$ and fits the exponent 3/2; (5) runs 8 initial conditions on 1, 2 and 4 processes and the tiny-job slowdown, under the `__main__` guard; (6) times Numba if it is installed (and says "not run" if not); (7) writes the speedup table and asserts its ordering; (8) fits the scaling exponents of all-pairs against `cKDTree`; (9) runs the test suite through `pytest.main` and the command line in a fresh interpreter, and checks its provenance record. It runs in about 25–35 s here; sizes are printed at the top.

**Build it yourself:** before opening `nbody/src/nbody/forces.py`, write Builds 1–3 from §3.2 and §3.5 and make them pass `tests/test_forces.py`. Then "Your turn": a cutoff force and a Barnes–Hut tree, measured against the workspace; `chunksize` variations; a `--profile` flag (Exercise B4); a Numba port of an earlier project's hot loop.

## 8. Where you'll struggle

- **You will optimize before profiling, and speed up the 5%.** Everyone does; the kick and drift loops of Example 4.1 *look* slow and take 0.3%. The fix is mechanical: no optimization without a profile in hand and a timing before and after, written down.
- **"Python is slow" is true per operation and irrelevant once the loop is in NumPy or Numba.** Rewriting a whole project in C++ because one loop was slow is the expensive version of not profiling. Push the one hot loop into compiled code and keep the rest in Python.
- **Threads will not speed up your CPU-bound Python, and you will try anyway.** The code runs, nothing crashes, and it is no faster (§3.6's table). Check `top`: one core at 100%. Use processes, or release the GIL by doing the work in NumPy or Numba.
- **Forgetting the `__main__` guard.** On Linux with `fork` nothing goes wrong, so the bug ships; on a Mac or Windows (`spawn`), or on Linux from Python 3.14 (`forkserver`), every worker re-runs your script's top level, and you get the `RuntimeError` of §3.6, a `BrokenProcessPool`, and every top-level side effect (prints, plots, files written) repeated once per worker. Put the guard in before the first pool.
- **Numba's error messages are terrible.** A `TypingError` buries the one useful line under a page of type-inference detail. Start with a ten-line function, get it compiling, and grow it; when it fails, delete half the body until it compiles (bisection again, PY-5 §3.8).
- **float32 is a trade you must justify, not a free speedup.** It halved the time of Build 4 and put a floor of about $10^{-7}$ under the energy error. Decide what accuracy you need first, then check it with the float64 run as reference.
- **The workspace returns its own buffer.** `history.append(work(pos, mass))` stores the same view every step, so the list ends up holding the last value $N$ times. Any function that returns a preallocated buffer must say so in its docstring, and callers must `.copy()` (Exercise A3).
- **A package feels like ceremony for a 200-line simulation.** It is the difference between code you can reuse next month and code you rewrite: the first time you `import nbody` from a different project, or run the tests after a change and see one fail, it has paid for itself.
- **Git's first week is miserable, and then you cannot work without it.** Detached HEAD, merge conflicts, the staging area: everybody is confused by them once. Keep commits small, commit before every experiment, and when lost, `git status` tells you where you are and usually what to type next.
- **Reproducibility means someone else can run it.** If your seed is hard-coded inside a function, your code is not reproducible, it is just deterministic: nobody can rerun it with another seed to see whether the result survives, and nobody knows which version of the code made which figure. Seeds, parameters, versions and the commit go into the output file.

## 9. Check yourself

1. Why is `time.perf_counter` the right clock for timing code and `time.time` the wrong one?
2. Why report the minimum of `timeit` repeats rather than the mean?
3. In a `pstats` listing, a function has tottime 0.01 s and cumtime 3 s. Where is its time going?
4. Roughly how many bytes does a list of $10^6$ Python floats occupy, and an array of them? Why?
5. What does `np.fill_diagonal(r2, 1.0)` fix in the broadcast force, and why is the value 1.0 harmless?
6. Why was the `(3, N, N)` layout faster than `(N, N, 3)` for the same arithmetic?
7. Name two ways the in-place workspace can bite a caller.
8. Which of these does a thread pool speed up in default CPython, and why: a pure-Python loop, `np.sin` on a large array, waiting on a download?
9. What exactly goes wrong when a script that creates a process pool has no `__main__` guard and the start method is `spawn`?
10. When does a process pool make a program slower?
11. When would you expect Numba to beat well-vectorized NumPy, and when not?
12. What is the difference between `dependencies` in `pyproject.toml` and the output of `pip freeze`?

## 10. Pair with

- **Gorelick & Ozsvald, *High Performance Python* 2e (O'Reilly, 2020), ch. 1–2 (profiling), 6 (matrix and vector computation), 7 (compiling to C: Cython, Numba), 9 (multiprocessing).** The book-length version of §§3.1–3.7, with more profilers (`line_profiler`, `memory_profiler`, `perf`) and a careful treatment of memory and caches.
- **Python docs: `timeit`, "The Python Profilers" (`cProfile`, `pstats`), `concurrent.futures`, `multiprocessing` ("Contexts and start methods", "Programming guidelines"), and `dis`.** The authoritative statements of §3.6; the multiprocessing guidelines are where the `__main__` rule is written down.
- **Numba documentation, "A ~5 minute guide to Numba".** What `@njit` supports, `parallel=True`, and how to read its errors.
- **Slatkin, *Effective Python* 2e, ch. 7 (concurrency; items 52–64 selectively) and items 70–73 (profiling, `deque`, `bisect`, `heapq`).** Short, opinionated items on threads versus processes and choosing data structures.
- **Python Packaging User Guide, "Packaging Python Projects" (packaging.python.org).** The current, official walk-through of `pyproject.toml`, building and installing; more current than any book.
- **Wilson et al., "Good enough practices in scientific computing", *PLOS Computational Biology* 13(6): e1005510 (2017).** Read in full; it is short and §3.11 only condenses it.
- **Chacon & Straub, *Pro Git* 2e, ch. 1–3 (free at git-scm.com/book), and Software Carpentry, "Version Control with Git".** Git's model (snapshots, the staging area, branches as pointers) explained properly; the Carpentry lesson is the gentlest hands-on start.
