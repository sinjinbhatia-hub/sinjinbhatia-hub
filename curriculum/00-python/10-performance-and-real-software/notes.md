# PY-10: Performance, the Runtime, and Writing Real Scientific Software

**Time:** ~10 hours · **Prerequisites:** PY-8 (so PY-3's functions and modules, PY-5's tests and logging, PY-6's arrays); PY-9 optional; CM-5 helpful for the physics · **Unlocks:** the working habits of every later project; take it before ML-7 (training loops) and SM-6 (Monte Carlo)

> **Conventions in this module.** Every output below came from a real run in this container (Python 3.13.16, NumPy 2.5.3, SciPy 1.18.1, numba 0.68.1, 4 CPUs, shared with other jobs); blocks run top to bottom in one session from the module folder, and the `text` block after a block is its output. **Every timing here is the writer's measurement on this container**; `project.py` re-measures and prints its own, which will differ. Only ratios transfer, and on a shared machine even they move by tens of percent. Wall time: `time.perf_counter`; microbenchmarks: `timeit`, reporting the **minimum** of repeats; where time goes: `cProfile` + `pstats`. A **speedup** is $t_{\text{baseline}}/t_{\text{version}}$, baseline named, per force evaluation at $N = 200$ unless stated. Physics: $G = 1$, total mass 1, positions `pos` as an `(N, 3)` float64 array (row $i$ is body $i$), masses `mass` shape `(N,)`, and the acceleration of body $i$ is
> $$\mathbf a_i = G\sum_{j\ne i} \frac{m_j(\mathbf x_j - \mathbf x_i)}{\left(|\mathbf x_j - \mathbf x_i|^2 + \varepsilon^2\right)^{3/2}}, \tag{3.1}$$
> with Plummer softening $\varepsilon$ ($\varepsilon = 0$ is exact gravity). Pairwise differences are $d_{ij} = \mathbf x_j - \mathbf x_i$. Any code that starts processes sits under `if __name__ == "__main__":`. The package `nbody/` lives in this folder; `import nbody` works after `sys.path.insert(0, "nbody/src")` or `pip install -e nbody` (Build 8).

## 0. Where this sits

PY-3 gave you `lru_cache` (§3.6), modules and packages (§3.8), `venv` and `pip` (§3.10) and an `argparse` tool (Example 4.6); PY-5 gave you pytest (§3.10) and `logging` (§3.9); PY-6 gave you broadcasting, the pairwise-distance array of its Build 4, and the rule "loop over the short axis" (§3.10). This module puts them together on one program and adds measurement, the interpreter's internals, parallelism, compilation, packaging and version control. Every later project is expected to be written this way; ML-7's training loops and SM-6's Metropolis sweeps are where the performance half pays off first.

## 1. The problem that forced this idea

You want the gravitational $N$-body problem of CM-5 §3.13 for $N = 500$ bodies and 2000 leapfrog steps, and then for many initial conditions, because one run of a chaotic system tells you almost nothing. The first version is the obvious double loop over Python lists. In the run below (Build 1) one force evaluation at $N = 500$ took 89 ms, so 2000 steps take 3 minutes; written as a loop over NumPy array elements, which is what most first versions look like, it took 315 ms, so 10 minutes. A hundred initial conditions: most of a day. The same physics with the loop pushed into NumPy and its temporaries removed runs the 2000 steps in 7 s (Build 3), and a hundred initial conditions on four cores in about three minutes.

Two things stand in the way. First, **you will optimize the wrong thing.** Asked where the slow program spends its time, most people point at the code that *looks* slow: the Python loops of the kick and drift updates. A profiler says they take 0.3% of the run (Example 4.1). So the rule here is absolute: measure, change one thing, measure again.

Second, **most scientific code cannot be rerun**, not by a referee, not by a labmate, not by its author six months later. It lives in one file that works from one directory, with the seed set inside a function, parameters edited by hand between runs, figures nobody can trace to the code that made them, and no test that would notice a change breaking the physics. Making a simulation fast and making it software are one discipline: being able to say, with evidence, what the code does.

So the module has a spine: one $N$-body code, nine builds. Profile it, vectorize it, remove its temporaries, try float32, run it on several processes, compile it with Numba, replace its $O(N^2)$ neighbour search with a tree, make it an installable package with tests and a command line, and put it under git.

## 2. Intuition

**Why Python is slow per operation.** In C, `a + b` on two doubles is one machine instruction on two registers. In Python, `a` and `b` are names bound to heap objects (PY-1), each with a type pointer and a reference count, and `a + b` executes the bytecode `BINARY_OP`: look up the type of `a`, find its `__add__`, check the type of `b`, unbox two doubles, add, allocate a new float for the result, adjust reference counts. Build 1 measures the consequence: about 300 ns per pair term of eq. (3.1), some twenty bytecodes.

```text
C:       a + b  ->  one add instruction on two registers
Python:  a + b  ->  BINARY_OP: type(a)? -> float.__add__ -> type(b)? -> unbox, add
                    -> allocate a new float -> refcounts                 (per number)
NumPy:   x + y  ->  ONE dispatch for the whole array, then a C loop over
                    contiguous doubles                                    (per array)
```

The fix is never "tighter Python". It is to pay the dispatch **once per array instead of once per number** (NumPy, PY-6), to compile the loop so there is no dispatch (Numba), or to run several interpreters at once (processes).

**Profiling is asking instead of guessing.** A profiler counts calls and time per function. Optimize something that takes 5% of the run and the best possible result is a 5% speedup.

**Threads versus processes.** CPython has a **global interpreter lock** (GIL): one thread at a time executes Python bytecode, so four threads running a Python loop take turns on one core. Threads help when they mostly *wait* (disk, network, `sleep`), and when the work runs in C code that releases the lock (many NumPy operations). CPU-bound Python needs **processes**, each with its own interpreter and GIL, at the price of copying (pickling) inputs and outputs between them:

```text
threads:    [ one process, one GIL ]  T1 runs | T2 waits | T3 waits    -> one core busy
processes:  [P1: GIL] [P2: GIL] [P3: GIL] [P4: GIL]                    -> four cores busy
            parent --pickle(args)--> worker --pickle(result)--> parent   (the overhead)
```

**Algorithms beat constants.** Vectorizing makes $O(N^2)$ cheaper per pair; a tree makes it $O(N\log N)$. At large $N$ no constant factor survives a better exponent.

**A package is a folder with a contract.** `pyproject.toml` says what the code is called, its version and what it needs; `tests/` say what it promises; the README says how to install, run and test it; git says which version made which result. "It works on my machine" becomes "it works on any machine that runs these three commands".

You should now be able to guess the results: one to two orders of magnitude from vectorizing, more from cutting memory traffic, near-linear speedup from processes until the cores run out, a *slowdown* from processes when jobs are tiny, and running-time exponents of about 2 for all-pairs against about 1 for a tree.

## 3. Formal development

### 3.1 Measuring time

**Two clocks.** `time.time()` is the calendar clock: it can be adjusted (NTP corrections) and so can jump. `time.perf_counter()` is monotonic with the finest resolution available: the clock for timing code. `time.process_time()` counts this process's CPU time; its gap to `perf_counter` is time spent waiting.

**`timeit`** runs a statement `number` times per measurement, `repeat` measurements. From the shell, `python -m timeit -n 10000 -r 7 -s "xs = list(range(1000))" "sum(xs)"` printed `10000 loops, best of 7: 4.07 usec per loop` here; in Jupyter, `%timeit` does the same. **Report the minimum:** noise only ever adds time (other processes, interrupts, cache misses), so the minimum best estimates the code's cost.

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
sum(xs), us per call: ['3.99', '3.95', '3.99', '4.59', '7.62', '7.45', '4.01']
```

The seven repeats of the same statement ranged from 3.95 to 7.62 µs: that spread is the machine, not the code. `bench` interleaves its candidates (A, B, C, A, B, C, ...), so a burst of load from another process lands on all of them instead of making one look slow.

**Benchmarking pitfalls.** (1) *Warm-up:* the first call pays for imports, page faults, Numba compilation and the interpreter's own specialization (Example 4.5); time it separately. (2) *Size:* a 10 kB array lives in cache and a 100 MB one does not; time at the size you will run. (3) *Variance:* take the minimum, and distrust differences under about 10% unless interleaved. (4) *Context:* the same function can cost twice as much in another process (Build 3). (5) *`python -O`* strips `assert` statements, so a script timed with `-O` is not the program you test (PY-5 §3.6).

### 3.2 Profiling, and the baseline (Build 1)

▶ **Build 1 (the baseline).** Write eq. (3.1) as a double loop over lists of `[x, y, z]` lists, and a leapfrog loop around it (kick half a step, drift, recompute forces, kick). Predict the time of one force evaluation at $N = 200$ and its scaling with $N$. (Reference: `accel_loop` in `nbody/src/nbody/forces.py`, and `nbody/src/nbody/baseline.py`.)

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
N = 200: 13.49 ms per evaluation, 339 ns per pair term
N = 500: lists 89 ms -> 2000 steps 3.0 min; same loop over ndarray elements 315 ms -> 10.5 min
ratio N=500/N=200: 6.60 (N^2 predicts 6.25)
```

■ **What just happened.** Each pair term costs about 300 ns: some twenty bytecodes, each a dispatch on boxed floats. Going from 200 to 500 bodies multiplied the time by 6.6, close to $N^2$'s 6.25. Indexing an ndarray element by element was 3.5 times *slower* than lists, because every `pos[j][0]` on an array builds a new array-scalar object; NumPy is fast only when you hand it whole arrays (PY-6 §3.10).

**`cProfile` and `pstats`.** `cProfile` records every function call and return. Run it as `python -m cProfile -o run.prof script.py`, or `cProfile.run("main()", "run.prof")`, or by enabling a `Profile` object around the code (Example 4.1); then `pstats.Stats("run.prof").sort_stats("cumulative").print_stats(10)` shows the top ten. Two columns matter:

- **tottime**: time spent in the function's own body, excluding the functions it called;
- **cumtime**: time from entering the function to leaving it, including everything it called.

Sort by **cumtime** to find the subtree that matters, then read **tottime** inside it to find the culprit. The profiler adds overhead per call (it inflates `math.sqrt`, called half a million times), so trust its proportions, not its absolute times. To find the *line*, the third-party `line_profiler` (`@profile`, then `kernprof -l -v script.py`) times every line.

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

A Python float is a 24-byte object (type pointer, reference count, the double) plus an 8-byte pointer in the list: 32 bytes per 8 bytes of data, scattered over the heap. The array stores the 8 bytes contiguously. For most array code, memory traffic rather than arithmetic is the bottleneck (Build 3).

### 3.4 Algorithms and data structures first

Check the algorithm first; it is the only change whose payoff grows with $N$.

- **Complexity.** All-pairs gravity is $O(N^2)$. The Barnes–Hut tree code (Barnes and Hut, *Nature* 324, 1986) groups distant bodies into octree cells and treats each cell as one mass at its centre of mass, for $O(N\log N)$ at a controlled error; the fast multipole method reaches $O(N)$. For a *short-range* force only neighbours within a cutoff matter, and a k-d tree finds them in $O(N\log N)$ (Build 7).
- **The right data structure.** PY-2 §3.5: `x in list` is $O(n)$, `x in set` is $O(1)$, so a membership test inside a loop silently makes it $O(n^2)$. `deque` for queues, `dict` for lookups, sorted arrays with `np.searchsorted` for range queries.
- **Early exit.** Stop when the answer is known (`any()`, `break`, a convergence test); do the cheap rejection first (compare $r^2$ with $r_c^2$ before any square root).
- **Caching.** `functools.lru_cache` (PY-3 §3.6) memoizes a pure function of hashable arguments. For an expensive array, memoize to disk: if `path.exists()`, return `np.load(path)`; otherwise compute it, `np.save(path, result)` and return it (keep such caches in an ignored folder like `figures/`).

A disk cache goes stale silently when code or inputs change: put the parameters in the file name, and delete caches when the code changes. `pickle` saves arbitrary objects, but **unpickling runs code chosen by whoever wrote the file**: never load a pickle you did not make. That is why `np.load` refuses object arrays unless you pass `allow_pickle=True`.

### 3.5 Constant factors: vectorize, cut the memory traffic, choose the dtype

▶ **Build 2 (vectorize).** With PY-6 Build 4's pairwise-difference broadcast, write eq. (3.1) with no Python loop. Predict every intermediate's shape and size at $N = 200$, and the speedup. What goes wrong on the diagonal $i = j$ when $\varepsilon = 0$, and how do you fix it without an `if`?

The direct translation, `accel_broadcast` in the package, keeps PY-6's `(N, N, 3)` layout; Example 4.2 traces it shape by shape. The diagonal: with $\varepsilon = 0$, $r^2_{ii} = 0$ and $0/0$ is `nan`, which poisons every sum. Since $d_{ii} = \mathbf 0$, any finite weight there contributes nothing, so `np.fill_diagonal(r2, 1.0)` fixes it with no branch. That layout's inner axis has length 3, so every pass runs in loops of three elements and pays the loop overhead 40,000 times. Putting the coordinate axis **first**, `(3, N, N)`, makes every pass a run over 40,000 contiguous doubles (`accel_numpy`):

```python
eps = 0.05
x = np.ascontiguousarray(pos.T)                        # (3, N)
d = x[:, np.newaxis, :] - x[:, :, np.newaxis]          # (3, N, N): d[k, i, j] = x_kj - x_ki
r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + eps**2  # (N, N)
np.fill_diagonal(r2, 1.0)
w = mass / (r2 * np.sqrt(r2))                          # (N, N)
acc = (d * w).sum(axis=2).T                            # (3, N) -> (N, 3); G = 1
print(d.shape, acc.shape, np.array_equal(acc, accel_numpy(pos, mass, eps)))
```
```text
(3, 200, 200) (200, 3) True
```

▶ **Build 3 (remove the temporaries).** Every expression above allocates a fresh array. Write a class `ForceWorkspace(n)` that allocates every buffer once in `__init__` and in `__call__` writes into them with `out=` and in-place operators (`np.subtract(a, b, out=d)`, `r2 += tmp`, `np.sum(d, axis=2, out=acc)`). Predict the gain. (Reference: `ForceWorkspace` in `forces.py`.)

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
broadcast (N,N,3)         2244 us    6.0x the loop  rel. diff 0.0e+00
coordinate-first           546 us   24.7x the loop  rel. diff 6.2e-16
in place (workspace)       472 us   28.6x the loop  rel. diff 6.2e-16
N = 500, workspace: 3.56 ms per evaluation -> 2000 steps in 7.1 s
```

Layout bought a factor of 4 over the `(N, N, 3)` broadcast; the workspace then gained 16%, about what people guess. But a temporary can cost far more than its arithmetic. Memory fresh from the operating system arrives unmapped: the first write to each 4 kB page traps into the kernel, which maps and zeroes it (a *minor page fault*). glibc's allocator takes large blocks directly from the OS (`mmap`) and returns them on `free`, with a size threshold that moves with the process's allocation history. Count the faults, and force the worst case with glibc's `MALLOC_MMAP_THRESHOLD_` (Linux only):

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
for label, extra in (("fresh process, default allocator", {}),
                     ("fresh process, MALLOC_MMAP_THRESHOLD_=131072", {"MALLOC_MMAP_THRESHOLD_": "131072"})):
    out = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True, env={**os.environ, **extra})
    print(label + ":\n" + out.stdout.rstrip())
```
```text
this session:  accel_numpy 0 page faults per call, workspace 0
fresh process, default allocator:
  accel_numpy: 670 us, 2 page faults per call
  workspace: 442 us, 0 page faults per call
fresh process, MALLOC_MMAP_THRESHOLD_=131072:
  accel_numpy: 1870 us, 1181 page faults per call
  workspace: 425 us, 0 page faults per call
```

■ **What just happened.** With every temporary fresh from the OS, `accel_numpy` paid 1181 page faults (4.8 MB) per call and ran 4.4 times slower than the workspace, which did not change. `project.py`'s process sits in between: 593 faults per call, `accel_numpy` 2.7 times slower than the workspace (in the writer's bisection, an earlier `import matplotlib` was enough to switch this on). Builds 2 and 3 are one arithmetic with three memory behaviours; the workspace's virtue is that its cost does not depend on allocator history you cannot see. Time code in the process it will run in. The price: the workspace returns a *view of its own buffer*, overwritten by the next call; `.copy()` what you keep (Exercise A3).

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
float64 444 us, float32 199 us (2.23x); float32 relative force error 3.5e-07
```

2.2 times as fast, and seven digits gone. `project.py` Part 3 prices it. Total momentum, conserved by Newton's third law up to rounding, drifts by $1.7\times10^{-16}$ in 500 float64 steps and $2.9\times10^{-8}$ in float32. The energy error at $t = 0.5$ falls like $\Delta t^2$ in float64 ($1.1\times10^{-6}$ at $\Delta t = 10^{-2}$, $1.1\times10^{-10}$ at $10^{-4}$) but in float32 stalls between $10^{-8}$ and $10^{-7}$, where rounding adds as much error per step as truncation removes. Fine for a picture, a network's weights (ML-7) or a Monte Carlo with statistical error $10^{-3}$; wrong for a long orbit.

**Other constant-factor rules** (Exercises A4, B2 measure some):

- **Preallocate.** `np.append` in a loop copies the array every time: $O(n^2)$. Collect in a list and convert once, or `np.empty(n)` and fill.
- **Hoist lookups out of hot loops.** `sqrt = math.sqrt` before the loop saves a global and an attribute lookup per call; inlining a tiny helper saves a call. Small gains, only in the loop the profile names.
- **Use BLAS.** `@` (and `np.einsum(..., optimize=True)` when it maps to a matrix product) calls compiled, cache-blocked, multithreaded BLAS (PY-4 Exercise C3). Set `OMP_NUM_THREADS=1` when you also run several processes.
- **Contiguity.** Passes are fastest along the last axis of a C-contiguous array (PY-6 §3.1); `np.ascontiguousarray(pos.T)` copies the strided view once into traversal order, as `accel_numpy` does.

### 3.6 The runtime: bytecode, reference counting, the GIL

**Bytecode.** CPython compiles each function to bytecode for a stack machine and executes it one instruction at a time; `dis.dis(f)` shows it (Example 4.5). Since 3.11 hot instructions are *specialized* in place (PEP 659): after a few float calls, `BINARY_OP` becomes `BINARY_OP_ADD_FLOAT`, skipping the type dispatch ("What's New in Python 3.11" reports 1.25× on average over 3.10). Each instruction still allocates its result, so compiled code stays far ahead.

**Reference counting and the garbage collector.** Every object counts the references to it and is freed the moment the count reaches zero: `del a`, or rebinding the name, frees an array at once, and Build 2's temporaries die at the end of each statement. Reference cycles never reach zero; the cyclic **garbage collector** (`gc`) finds them periodically. In numerical code the practical point is the reverse: a reference you forgot (a list of views, a cache) keeps a large array alive.

**The GIL.** Reference counts are updated without locks, so CPython protects them with one global interpreter lock: a thread must hold it to execute bytecode. A thread releases it while blocked on I/O or `time.sleep`, and C extensions may release it while they work on raw memory, as NumPy does in many ufunc loops and in BLAS. (CPython 3.13 also ships an *experimental* free-threaded build, PEP 703; the default build keeps the GIL, and `sys._is_gil_enabled()` says which you have.) Measured:

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
pure-Python loop             4 tasks: serial 0.38 s, 4 threads 0.47 s -> 0.8x
NumPy np.sin on 4e6 doubles  4 tasks: serial 0.17 s, 4 threads 0.04 s -> 4.1x
time.sleep(0.2)              4 tasks: serial 0.80 s, 4 threads 0.20 s -> 4.0x
```

Threads made the Python loop *slower* (0.8×: the tasks took turns on the lock, plus switching costs), sped NumPy up 4.1× (the ufunc released the lock while it worked), and sleeping 4.0×. So: **threads for I/O** (downloads, many files, instruments), **processes for CPU-bound Python**, and neither before you have vectorized. **`asyncio`**, one thread running many coroutines that `await` I/O, is built for thousands of simultaneous network connections; numerical work never waits, so it is not for this.

▶ **Build 5 (processes).** Write a script that runs `nbody.integrate.simulate(seed)` for 8 seeds serially and on `ProcessPoolExecutor(max_workers=w)`, $w = 1, 2, 4$, under the `__main__` guard. Predict the speedups on 4 cores, and for millisecond jobs.

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
4 CPUs; 8 runs, N = 200, 400 steps; serial 1.56 s
  1 worker(s): 1.34 s, speedup 1.17x
  2 worker(s): 0.74 s, speedup 2.10x
  4 worker(s): 0.43 s, speedup 3.61x
tiny jobs (N = 8, 5 steps): serial 4.0 ms, 4 workers 15.9 ms
```

■ **What just happened.** Each worker is a separate process with its own import of `nbody`; `pool.map` pickles each seed to a worker and each result dictionary back. Four workers gave 3.6× here (`project.py`'s run: 2.9×; the other jobs on this machine take cores). One worker should be slightly *slower* than serial (same work plus start-up and pickling); it came out 1.17× faster, which tells you the size of the noise here (`project.py`: 1.03×). For 4-ms jobs the pool lost, 15.9 ms against 4.0 ms: starting four processes costs more than the work. Results equal the serial ones because each job receives its seed explicitly; a seed drawn from a global generator inside the job would make results depend on which worker ran what.

**The `__main__` guard.** Workers start one of three ways. `fork` (the Linux default up to 3.13) copies the parent. `spawn` (Windows, macOS) and `forkserver` (the Linux default from 3.14) start a fresh interpreter that **imports your main script** (as `__mp_main__`) to find your functions, so every top-level statement runs again in every worker; if one creates a pool, Python stops it:

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
top level runs in SpawnProcess-2
top level runs in SpawnProcess-1
RuntimeError: 
        An attempt has been made to start a new process before the
        current process has finished its bootstrapping phase.

        This probably means that you are not using fork to start your
        child processes and you have forgotten to use the proper idiom
        in the main module:

            if __name__ == '__main__':
concurrent.futures.process.BrokenProcessPool: A process in the process pool was terminated abruptly while the future was running or pending.
```

With the guard, workers import the script, define its functions, skip the block and wait for work. One line, mandatory; `nbody/src/nbody/__main__.py` has it. Two more facts: `pool.map(f, items, chunksize=k)` sends items in batches of $k$, amortizing the per-task cost (Example 4.4, Exercise C2); and pickling a large read-only array to every worker can cost more than the work, which is what `multiprocessing.shared_memory` is for (or let each worker load its own data from a file name or a seed).

### 3.7 Compiling: Numba (Build 6)

Numba compiles a function to machine code (through LLVM) the first time it is called with given argument types. `@njit` compiles only if *every* operation can be typed: numbers, NumPy arrays, tuples, loops, most of `math` and much of NumPy. Mixed-type lists, arbitrary dicts, SciPy, pandas and your own classes are out.

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
numba 0.68.1: first call 0.54 s, then 153 us per evaluation (88x the loop, 3.1x faster than the workspace); rel. diff 0.0e+00
TypingError, 20 lines. The first, and the useful one:
Failed in nopython mode pipeline (step: nopython frontend)
TypingError: got an unexpected keyword argument 'axis'
```

Half a second to compile, then 3.1 times faster than the workspace. The brief for this module predicted "often similar for this problem"; here it was not, because the compiled loop reads the positions once and keeps each partial sum in a register, where NumPy makes a dozen passes over 320 kB arrays. The general rule: **a compiled loop beats NumPy when the vectorized form needs many passes over large temporaries, or cannot be written at all** (a step that depends on the previous one, an early exit, a branch per element). In SM-6's Metropolis sweep each spin flip depends on the last, and Numba is the way to compiled speed without writing C.

The `TypingError` is typical: correct NumPy that Numba does not support, and one useful line buried among twenty. **Start with tiny functions** and grow them. `@njit(parallel=True)` with `numba.prange` on the outer loop splits its iterations over Numba's own threads; here the loop over $i$ is independent, so that works (Exercise C3). `@njit(cache=True)` saves the compiled code to disk.

**Other compilers, named.** Cython compiles a typed Python superset to a C extension (SciPy and scikit-learn use it); a hand-written C/C++ extension (or `pybind11`) gives full control at the price of a build system. **When to stop:** when the profile says the hot spot is no longer yours (BLAS, I/O), or the program is fast enough for the question. A 10-s simulation does not need to run in 3 s; a week-long one does.

### 3.8 Trees for neighbours (Build 7)

For a short-range force only pairs closer than a cutoff $r_c$ interact, but checking every pair is still $O(N^2)$. A **k-d tree** splits space recursively, one coordinate at a time; a query for "all points within $r_c$" descends only into cells that overlap the sphere, so at fixed density it costs $O(\log N)$ per body.

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
N =  2000:   3874 pairs, all-pairs 0.050 s, tree 2.2 ms, same pairs: True
N =  8000:  15952 pairs, all-pairs 0.568 s, tree 8.6 ms, same pairs: True
short-range force, tree vs brute: rel. diff 7.6e-17
```

Quadrupling $N$ multiplied the all-pairs time by 11 and the tree time by 4. `project.py` Part 8 fits the exponents over $N = 2000$ to $32{,}000$; in the writer's run they were 1.92 for all-pairs (a little under 2: fixed per-call costs still show at small $N$) and 1.10 for the tree, and at $N = 16{,}000$ the tree was 340 times faster. No constant factor gets you that.

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

**Why `src/`.** With the package next to `tests/`, a test's `import nbody` silently picks up the folder in the current directory, and tests pass even when the *installed* package is broken. The `src/` layout forces imports through the installation (or an explicit path, which `pyproject.toml` gives pytest).

**Editable install.** `pip install -e .` installs the package so that `import nbody` resolves to your source folder: edits take effect without reinstalling. Always inside a virtual environment (PY-3 §3.10):

```text
$ pip install -e ".[test]"
Obtaining file:///.../pipdemo/nbody
  [... 8 lines of build steps; "Requirement already satisfied" for numpy, scipy, pytest and its dependencies ...]
Successfully built nbody
Installing collected packages: nbody
Successfully installed nbody-0.1.0
$ nbody --n 100 --steps 500 --seed 1 --out run1.json
INFO nbody: seed 1: |dE/E| = 2.31e-08, |dP| = 4.7e-17, half-mass radius 0.6983
INFO nbody: 1 run(s), N = 100, 500 steps, 1 worker(s): 0.06 s
INFO nbody: wrote run1.json
```

(Lines in `[...]` are elided from the real transcript. The build step fetches `setuptools`, so it needs network access once.)

**The command line.** `cli.py` holds a `RunConfig` dataclass with one field per parameter (`n`, `steps`, `dt`, `eps`, `seed`, `ensemble`, `workers`, `dtype`), a classmethod `from_toml(path)` that reads a `[run]` table with `tomllib` and rejects unknown keys (a typo must be an error), and an `argparse` parser (PY-3 Example 4.6) with one flag per field.
`main()` reads the TOML file if `--config` is given, overrides it with the flags given (`dataclasses.replace`), configures `logging` (`-v` for DEBUG), runs, and with `--out` writes one JSON file with the configuration, the results and the **provenance**: package, Python and NumPy versions, platform, and `git rev-parse --short HEAD`. In Example 4.7's repository:

```text
"provenance": {"nbody": "0.1.0", "python": "3.13.16", "numpy": "2.5.3",
               "platform": "Linux-6.18.44-fc-v185-x86_64-with-glibc2.39", "git_commit": "e3c48c8"}
```

A result file that names its commit can be reproduced. Documentation comes in three layers: `--help` (free from `argparse`), the README (install, run, test), and docstrings that state what a function computes, with the equation number.

**Reproducibility, concretely.** (1) Seeds are *parameters* passed in from the top (`simulate(seed, ...)`, `--seed`), never `default_rng(42)` buried in a function. (2) Record the environment: `pip freeze > requirements.txt` pins every installed version; `pyproject.toml`'s `dependencies` state *minimum* versions (libraries state ranges, applications pin). (3) Write commit and versions into every output file. (4) `logging`, not `print`, for progress (PY-5 §3.9): it can be silenced, redirected and timestamped without touching the code. (5) A formatter and linter (`ruff format`, `ruff check`) end style arguments and catch unused imports and undefined names; continuous integration (GitHub Actions and similar) runs `pytest` on every push, so a broken commit is noticed in minutes.

### 3.10 Notebooks versus scripts

Notebooks are good for exploring and bad as the record of a result: cells run out of order, deleted cells leave variables alive, and the JSON file diffs badly. The curriculum's compromise is the `# %%` script: a plain `.py` file that runs top to bottom from the shell and in CI, and that VS Code, Spyder or PyCharm run cell by cell; `jupytext` converts it to and from `.ipynb`. In a notebook, **"Restart kernel and run all" before you trust or share anything**, and move anything you reuse into a module the notebook imports.

### 3.11 Version control (Build 9)

Git records snapshots ("commits") of a folder, each with a message, an author, a time and a pointer to its parent. The everyday commands:

`git init` (make the folder a repository), `git status` (what changed, what is staged), `git add FILE` (stage a change), `git commit -m "..."` (record the staged snapshot), `git log --oneline` (history), `git diff` (unstaged changes; `--stat` summarizes), `git branch NAME` and `git checkout NAME` (start a line of work and switch to it; `git switch` is the newer spelling), `git merge NAME` (bring another branch's commits into this one).

**Commit** source, tests, `pyproject.toml`, README, configuration and small inputs. **Ignore** (`.gitignore`) whatever the code regenerates (figures, caches, `__pycache__`, `*.egg-info`), large data (record where it lives and its checksum instead) and secrets. **Commit messages say why**, in the imperative: "Record versions and commit in output files, so a result names its code", not "update cli.py"; the diff already says what.

A **remote** is a copy elsewhere (GitHub, GitLab, a lab server): `git remote add origin URL`, then `git push` to send commits and `git pull` to fetch and merge others'. It is your backup; push at the end of every session. `git bisect` (PY-5 §3.8) finds the commit that broke something by binary search, which works only if commits are small and each one runs.

▶ **Build 9.** In a copy of the package: `git init`, ignore the build products, make three commits with messages that say why, branch to try float32, merge it back. Example 4.7 is the real session.

**"Good enough practices", condensed** (Wilson et al. 2017). Keep raw data raw and generate everything else from it by code; one project, one folder, with `src/`, `data/`, `results/` and `doc/`; decompose programs into functions, eliminate duplication, and use well-maintained libraries rather than rewriting them; give functions and files meaningful names; state dependencies explicitly; control behaviour with parameters, never by commenting code in and out; provide a small example data set that runs fast; put everything a human wrote under version control, in small commits shared often; add a README, a LICENSE and a citation file.

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
         639465 function calls in 0.303 seconds

   Ordered by: cumulative time
   List reduced from 13 to 8 due to restriction <8>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000    0.303    0.303 baseline.py:40(run_python)
       11    0.176    0.016    0.225    0.020 forces.py:26(accel_loop)
        2    0.049    0.025    0.076    0.038 baseline.py:29(energy_loop)
   477600    0.054    0.000    0.054    0.000 {built-in method math.sqrt}
   159200    0.022    0.000    0.022    0.000 baseline.py:35(<genexpr>)
       20    0.001    0.000    0.001    0.000 baseline.py:13(kick)
     2200    0.001    0.000    0.001    0.000 {method 'append' of 'list' objects}
       10    0.000    0.000    0.000    0.000 baseline.py:21(drift)
```

Read it from the top. `run_python` is the whole run, 0.303 s (cumtime). `accel_loop`, called 11 times (once before the loop, once per step), has cumtime 0.225 s: 74% of the run, most of it in its own body (tottime 0.176 s) and the rest in `math.sqrt`, called once per pair term. `energy_loop`, called only twice but $O(N^2)$ each time, takes 0.076 s: a quarter of the run spent on a *diagnostic*, which nobody guesses before seeing it. Inside it, the `<genexpr>` row (159,200 calls) is the generator that unpacks `dx, dy, dz`, a harmless-looking line costing 7% of the run. `kick` and `drift`, the loops that look slow, take 0.001 s together. So: vectorize the force (Builds 2–3), compute the energy rarely and vectorized, and leave kick and drift alone.

### Example 4.2: the broadcast, shape by shape

For `accel_broadcast` at $N = 200$, every intermediate and its memory:

```python
d = pos[np.newaxis, :, :] - pos[:, np.newaxis, :]
r2 = (d * d).sum(axis=-1) + 0.05**2
np.fill_diagonal(r2, 1.0)
w = mass / (r2 * np.sqrt(r2))
wd = w[:, :, np.newaxis] * d
acc = wd.sum(axis=1)
for name, a in [("pos", pos), ("pos[np.newaxis, :, :]", pos[np.newaxis, :, :]),
                ("pos[:, np.newaxis, :]", pos[:, np.newaxis, :]), ("d", d), ("r2", r2), ("w", w),
                ("w[:, :, np.newaxis] * d", wd), ("acc = (...).sum(axis=1)", acc)]:
    print(f"{name:<24} shape {str(a.shape):<14} {a.nbytes / 1e3:6.1f} kB  view: {a.base is not None}")
print("d[3, 7] == pos[7] - pos[3]:", np.array_equal(d[3, 7], pos[7] - pos[3]),
      "| w[5, 9] == m_9 / r_59^3:", np.isclose(w[5, 9], mass[9] / r2[5, 9] ** 1.5))
```
```text
pos                      shape (200, 3)          4.8 kB  view: False
pos[np.newaxis, :, :]    shape (1, 200, 3)       4.8 kB  view: True
pos[:, np.newaxis, :]    shape (200, 1, 3)       4.8 kB  view: True
d                        shape (200, 200, 3)   960.0 kB  view: False
r2                       shape (200, 200)      320.0 kB  view: False
w                        shape (200, 200)      320.0 kB  view: False
w[:, :, np.newaxis] * d  shape (200, 200, 3)   960.0 kB  view: False
acc = (...).sum(axis=1)  shape (200, 3)          4.8 kB  view: False
d[3, 7] == pos[7] - pos[3]: True | w[5, 9] == m_9 / r_59^3: True
```

The two `np.newaxis` operands are views of shapes `(1, N, 3)` and `(N, 1, 3)` and cost nothing (PY-6 §3.5). Broadcasting them allocates `(N, N, 3)`: 960 kB. `mass`, shape `(200,)`, against `r2`, shape `(200, 200)`, is aligned from the right and acts as `(1, 200)`, so $m_j$ varies along columns, as eq. (3.1) needs. Write `mass[:, np.newaxis]` by mistake and you get $m_i$ instead: wrong, and the comparison with the loop catches it at once. Four of these arrays are the size of the problem; that is Build 3's memory traffic.

### Example 4.3: the speedup table

From the writer's run of `project.py` (it writes `figures/speedups.txt` on every run; yours will differ):

```text
N = 200, this machine (4 CPUs), Python 3.13.16, NumPy 2.5.3
version                                           time per eval   speedup
baseline: pure-Python loop over lists                 9547.6 us      1.0x
vectorized, (N, N, 3) broadcast                       2891.2 us      3.3x
vectorized, coordinate axis first                     1189.6 us      8.0x
vectorized + in place (workspace)                      440.2 us     21.7x
in place, float32                                      204.6 us     46.7x
ensemble of 8, serial (per evaluation)                 478.3 us     20.0x
ensemble of 8, 1 process(es) (per evaluation)          463.1 us     20.6x
ensemble of 8, 2 process(es) (per evaluation)          250.8 us     38.1x
ensemble of 8, 4 process(es) (per evaluation)          167.0 us     57.2x
Numba @njit loop (after compiling)                     153.3 us     62.3x
```

Each row is a minimum over interleaved repeats. The ensemble rows divide the wall time of 8 runs by their 8 × 401 force evaluations, so they include the integrator's bookkeeping and the pool's start-up; the serial ensemble row is the workspace plus that overhead. Read down the column: layout (3.3× → 8.0×), no temporaries (21.7×), half the bytes (46.7×), four processes (57.2×, which is 2.9× over the serial ensemble on this shared machine), a compiled loop (62.3×). The brief for this module asked for "vectorized ≥ 20×"; in `project.py`'s process the coordinate-first version reached 8.0×, because it paid 593 page faults per call there (Build 3), while in the notes session, without the faults, it reached 24.7×. Only the in-place version passes 20× in both. `project.py` asserts the ordering baseline ≤ vectorized ≤ in place and a threshold of 5× for the vectorized version, and prints the measured values.

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
serial list comprehension:    0.53 ms
4 processes, chunksize     1: 2228.92 ms
4 processes, chunksize   100:   41.15 ms
4 processes, chunksize  2500:   16.82 ms
```

With `chunksize=1` each task is pickled, sent through a pipe, unpickled, run and sent back: about 0.2 ms of overhead (2.2 s over 10,000 tasks) around a fraction of a microsecond of work, 4000 times slower than the list comprehension. Batching cuts the per-task cost (41 ms, then 17 ms) but cannot remove the pool's start-up. Rule of thumb: each task should take at least tens of milliseconds, and the total should dwarf the pool's start-up.

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

Three source lines, six instructions. `RESUME` marks the function's start (where tracing and optimization checks happen); `LOAD_FAST_LOAD_FAST` pushes both arguments onto the evaluation stack; `BINARY_OP 0 (+)` pops them, dispatches on their types and pushes the result; `STORE_FAST` binds it to `c`; `LOAD_FAST` and `RETURN_VALUE` return it. The same `BINARY_OP` adds ints, concatenates strings and broadcasts arrays: the decision is made at run time, every time. After a hundred calls with floats the adaptive interpreter has rewritten it to `BINARY_OP_ADD_FLOAT`, a fast path that checks "both floats?" and adds, falling back if the types change. Still one instruction per operation, each allocating its result; Numba removes the instructions altogether.

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

- `[build-system]`: the tool that builds the installable wheel, and what to install before running it (PEP 517/518). Here setuptools; `hatchling` and `flit-core` are alternatives.
- `name`: the distribution name for `pip install` and `pip show`; the *import* name comes from the folder `src/nbody`. `version`: bump it when behaviour changes; it appears in every provenance record. `description`, `readme`, `authors`: metadata.
- `requires-python = ">=3.11"`: pip refuses older Pythons up front instead of failing at `import tomllib` (new in 3.11).
- `dependencies`: what `import nbody` needs, as *minimum* versions. A library that pins exact versions cannot coexist with others; exact pins belong in an application's `requirements.txt`.
- `[project.optional-dependencies]`: extras, installed by `pip install -e ".[test]"` or `".[fast]"`. Numba is optional because `forces.py` works without it.
- `[project.scripts]`: installs an `nbody` command that calls `nbody.cli.main()`; `python -m nbody` goes through `__main__.py` instead.
- `[tool.setuptools.packages.find] where = ["src"]`: the `src/` layout.
- `[tool.pytest.ini_options]`: tests live in `tests/`, and `src/` goes on the import path, so `pytest` works in a fresh clone before installation.

The file is checked twice: `pip install -e .` succeeded with it (Build 8), and `project.py` Part 9 runs pytest under these options.

### Example 4.7: a git session on the package

Run in a copy of the package in a temporary directory, with a virtual environment holding an editable install of that copy. Comment lines `#` mark edits made in an editor between commands. The transcript is real; three long outputs are shortened, as marked.

```text
$ git init -q -b main
$ git config user.name 'You'
$ git config user.email you@example.com
$ pip install -q -e .
$ python -m pytest -q
..............                                                           [100%]
14 passed in 1.06s
$ git status --short --untracked-files=all
?? README.md
?? pyproject.toml
?? src/nbody.egg-info/PKG-INFO
[... 5 more files in src/nbody.egg-info/ ...]
?? src/nbody/__init__.py
?? src/nbody/__pycache__/__init__.cpython-313.pyc
?? src/nbody/__pycache__/forces.cpython-313.pyc
?? src/nbody/__pycache__/integrate.cpython-313.pyc
?? src/nbody/forces.py
?? src/nbody/integrate.py
?? tests/__pycache__/test_forces.cpython-313-pytest-9.1.1.pyc
?? tests/__pycache__/test_integrate.cpython-313-pytest-9.1.1.pyc
?? tests/test_forces.py
?? tests/test_integrate.py
# __pycache__ and nbody.egg-info are build products: ignore them, and figures and run outputs
$ cat .gitignore
[... the six lines of nbody/.gitignore: figures/, __pycache__/, *.egg-info/, .pytest_cache/, *.prof, run*.json ...]
$ git status --short
?? .gitignore
?? README.md
?? pyproject.toml
?? src/
?? tests/
$ git add .
$ git commit -q -m "Add softened-gravity kernels, leapfrog/RK4 and their tests"
# write src/nbody/cli.py and src/nbody/__main__.py
$ git add src/nbody/cli.py src/nbody/__main__.py
$ git commit -q -m "Add a command line so runs are repeatable from the shell"
# edit cli.py: record versions and the git commit in every --out file
$ git diff --stat
 src/nbody/cli.py | 23 +++++++++++++++++++++--
 1 file changed, 21 insertions(+), 2 deletions(-)
$ git commit -q -am "Record versions and commit in output files, so a result names its code"
$ git log --oneline
4ee7cff Record versions and commit in output files, so a result names its code
fc2b9c5 Add a command line so runs are repeatable from the shell
129bf4b Add softened-gravity kernels, leapfrog/RK4 and their tests
$ git branch float32-experiment
$ git checkout -q float32-experiment
# edit cli.py: add --dtype
$ git diff
[... 24 lines: three hunks adding the dtype field, the --dtype flag, and dtype=cfg.dtype in the run_ensemble call ...]
$ python -m nbody --n 200 --steps 300 --dtype float32
INFO nbody: seed 1: |dE/E| = 1.15e-07, |dP| = 3.4e-08, half-mass radius 0.7495
INFO nbody: 1 run(s), N = 200, 300 steps, 1 worker(s): 0.07 s
$ git commit -q -am "Allow float32 runs (--dtype); default stays float64, see energy floor"
$ git checkout -q main
# meanwhile on main: document how to reproduce a run
$ git commit -q -am "README: how to reproduce a run from its output file"
$ git merge float32-experiment -m "Merge float32-experiment: keep --dtype as an option"
Merge made by the 'ort' strategy.
 src/nbody/cli.py | 4 +++-
 1 file changed, 3 insertions(+), 1 deletion(-)
$ git log --oneline --graph --all
*   e3c48c8 Merge float32-experiment: keep --dtype as an option
|\  
| * 131d6e3 Allow float32 runs (--dtype); default stays float64, see energy floor
* | 87204fc README: how to reproduce a run from its output file
|/  
* 4ee7cff Record versions and commit in output files, so a result names its code
* fc2b9c5 Add a command line so runs are repeatable from the shell
* 129bf4b Add softened-gravity kernels, leapfrog/RK4 and their tests
$ python -m pytest -q
..............                                                           [100%]
14 passed in 0.91s
```

Before the `.gitignore`, git offers to track eleven generated files (bytecode caches, and the `egg-info` that `pip install -e` wrote); after it, only source. Each commit on `main` leaves a working package. The experiment stays on a branch until tested; meanwhile `main` moved (the README commit), so the merge made a commit with two parents, drawn by `--graph` (had `main` not moved, git would have fast-forwarded). The merged tree passes its tests and is identical to `nbody/` in this folder.

## 5. Why it matters and where it connects

- **CM-5, CM-8.** The physics here is CM-5's two-body problem (its Theorem 3.7, Kepler's third law, is `project.py`'s end-to-end test) and §3.13's $N$-body problem; leapfrog's energy error that oscillates but never grows, against RK4's steady growth (`project.py` Part 3), is CM-8 §3.11's symplectic integrator at work.
- **ML-7.** A training loop is Builds 1–3 again: batch into arrays, preallocate, profile first; float32 is the default there because statistical noise dwarfs rounding (Build 4). Parallel data loading is Build 5.
- **SM-6.** A Metropolis sweep cannot be vectorized (each step depends on the last): Build 6's Numba. Independent chains with explicit seeds: Build 5.
- **Every later project.** Seeds as parameters, a config dataclass, results that record their commit, tests in `tests/`, figures regenerated into an ignored `figures/`: Build 8 is how they are all set up.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are predictions and small fixes (a timing that lies, a profile to read, the workspace aliasing bug, preallocation, a missing `__main__` guard). B1–B6 are the module; B1 (profile and fix a slow program), B2 (memory traffic of the potential energy), B4 (a `--profile` flag with a test) and B6 (a git history for a bug fix) are essential. C1–C3: a cutoff force with error control, a chunksize model, Numba with `prange`.

## 7. Coding project

`project.py` drives `nbody/` through the builds and checks each stage: (1) times and profiles the baseline at $N = 200$; (2) checks the broadcast, coordinate-first and in-place forces against the loop to $10^{-12}$, times them and counts page faults; (3) checks momentum to $10^{-10}$, measures float32's momentum drift and energy floor, and compares leapfrog with RK4 over 100 orbits at equal cost; (4) verifies Kepler's third law to $10^{-4}$ and the exponent 3/2; (5) runs 8 initial conditions on 1, 2 and 4 processes, and the tiny-job slowdown; (6) times Numba, or says "not run"; (7) writes and checks the speedup table; (8) fits all-pairs against `cKDTree` exponents; (9) runs pytest and the command line, and checks the provenance record. About 25 s here.

**Build it yourself:** write Builds 1–3 before opening `nbody/src/nbody/forces.py`, and make them pass `tests/test_forces.py`. Then "Your turn": a cutoff force and a Barnes–Hut tree against the workspace; `chunksize` variations; a `--profile` flag; a Numba port of an earlier project's hot loop.

## 8. Where you'll struggle

- **You will optimize before profiling, and speed up the 5%.** Everyone does: kick and drift *look* slow and took 0.3% (Example 4.1). No optimization without a profile in hand and a timing before and after, written down.
- **"Python is slow" is true per operation and irrelevant once the loop is in NumPy or Numba.** Rewriting a project in C++ for one slow loop is the expensive version of not profiling.
- **Threads will not speed up your CPU-bound Python, and you will try anyway.** It runs, nothing crashes, it is no faster (§3.6), and `top` shows one core at 100%. Use processes, or do the work in NumPy or Numba.
- **Forgetting the `__main__` guard.** On Linux with `fork` nothing goes wrong, so the bug ships; under `spawn` (macOS, Windows) or `forkserver` (Linux from 3.14) every worker re-runs your top level, and you get §3.6's `RuntimeError`, a `BrokenProcessPool`, and every top-level side effect once per worker.
- **Numba's error messages are terrible.** Start with a ten-line function, get it compiling, grow it; when it fails, delete half the body until it compiles (bisection, PY-5 §3.8).
- **float32 is a trade you must justify, not a free speedup.** It halved Build 4's time and put a floor near $10^{-7}$ under the energy error. Decide the accuracy you need, then check against a float64 run.
- **Preallocated buffers alias.** `history.append(work(pos, mass))` stores the same view every step, so the list holds the last value $N$ times (Exercise A3). Document it, and `.copy()`.
- **A package feels like ceremony for a 200-line simulation.** It is the difference between code you reuse and code you rewrite; the first `import nbody` from another project, or the first failing test after a change, repays it.
- **Git's first week is miserable, and then you cannot work without it.** Detached HEAD, conflicts, the staging area confuse everyone once. Commit small, commit before every experiment, and when lost, `git status` says where you are.
- **Reproducibility means someone else can run it.** A seed hard-coded inside a function makes code deterministic, not reproducible: nobody can rerun it with another seed to see whether the result survives, or tell which code made which figure. Seeds, parameters, versions and the commit go into the output file.

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

- **Gorelick & Ozsvald, *High Performance Python* 2e (O'Reilly, 2020), ch. 1–2 (profiling), 6 (matrix and vector computation), 7 (compiling to C), 9 (multiprocessing).** The book-length version of §§3.1–3.7, with more profilers and a careful treatment of memory and caches.
- **Python docs: `timeit`, "The Python Profilers", `concurrent.futures`, `multiprocessing` ("Contexts and start methods", "Programming guidelines"), `dis`.** The authoritative semantics; the multiprocessing guidelines are where the `__main__` rule is written down.
- **Numba, "A ~5 minute guide to Numba".** What `@njit` supports, `parallel=True`, and reading its errors.
- **Slatkin, *Effective Python* 2e, ch. 7 (items 52–64, selectively) and items 70–73.** Short, opinionated items on threads versus processes, profiling, and data structures.
- **Python Packaging User Guide, "Packaging Python Projects" (packaging.python.org).** The official, current walk-through of `pyproject.toml`.
- **Wilson et al., "Good enough practices in scientific computing", *PLOS Computational Biology* 13(6): e1005510 (2017).** Read in full; §3.11 only condenses it.
- **Chacon & Straub, *Pro Git* 2e, ch. 1–3 (free at git-scm.com/book); Software Carpentry, "Version Control with Git".** Git's model (snapshots, staging, branches as pointers) explained properly, and the gentlest hands-on start.
