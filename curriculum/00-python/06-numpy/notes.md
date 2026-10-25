# PY-6: NumPy — Arrays, Memory, Broadcasting, Vectorization, and Views vs Copies

**Time:** ~10 hours · **Prerequisites:** PY-5 (LA-1 and LA-2 help; take this by LA-5 at the latest) · **Unlocks:** PY-7, PY-8, PY-10, LA-13, and every physics and ML project

> **Conventions in this module.** Python 3.11+; every output below came from a real run under Python 3.13 and NumPy 2.5.3 on 64-bit Linux. `import numpy as np`; random numbers come from `rng = np.random.default_rng(seed)`, never the legacy `np.random.seed`. Shapes are tuples `(n, m)`; **axis 0 runs down the rows**; strides are in bytes. A reduction **loses the axis you name** (§3.6). Code blocks run top to bottom in one session, so names from earlier blocks stay defined; the `text` block after a block is its real output. Timings come from one run and will differ on your machine; the ratios will be of the same order. Solve $A\mathbf{x}=\mathbf{b}$ with `np.linalg.solve`, never `inv` (§3.8; LA-2, LA-13). The default integer here is `int64`.

## 0. Where this sits

PY-1 gave you names and objects, PY-2 the list and its copying slices, PY-4 the data model (`__getitem__` with a tuple of slices, `__matmul__`), PY-5 asserts. NumPy is where they meet: an array is an object whose slices do *not* copy, whose operators loop in C, and whose bugs you find by asserting shapes. PY-7 plots arrays, PY-8 (SciPy) is built on them, PY-9's columns are arrays, PY-10 profiles the loops you could not remove; LA-2 onward does its numerical work here, LA-12's contractions are `einsum` (§3.8), and LA-13 builds on §3.9.

## 1. The problem that forced this idea

Sum the squares of a million numbers with a loop over a list, then with an array:

```python
import sys
import time

import numpy as np

Shape = tuple[int, ...]            # a shape or a strides tuple, for type hints
n = 1_000_000
xs = [float(i) for i in range(n)]
arr = np.arange(n, dtype=np.float64)
print(f"list: {sys.getsizeof(xs[1])}-byte float objects plus 8-byte pointers; array: {arr.itemsize} bytes per element, one block of {arr.nbytes:,} bytes")


def best(f, *args, repeat: int = 3) -> float:
    """Smallest wall time of `repeat` calls of f(*args), in seconds."""
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        f(*args)
        times.append(time.perf_counter() - t0)
    return min(times)


def sumsq_loop(v: list[float]) -> float:
    s = 0.0
    for t in v:
        s += t * t
    return s


assert np.isclose(sumsq_loop(xs), np.dot(arr, arr))
t_loop, t_np = best(sumsq_loop, xs), best(lambda v: np.dot(v, v), arr)
print(f"loop over the list: {t_loop * 1e3:7.1f} ms")
print(f"np.dot on the array: {t_np * 1e3:7.3f} ms  ({t_loop / t_np:,.0f}x faster)")
```
```text
list: 24-byte float objects plus 8-byte pointers; array: 8 bytes per element, one block of 8,000,000 bytes
loop over the list:    23.9 ms
np.dot on the array:   0.102 ms  (235x faster)
```

A list of floats is a list of *pointers* to separate 24-byte objects scattered through memory; each trip round the loop looks one up, checks its type, multiplies, and allocates a new float. An array is one contiguous block of raw 8-byte doubles with a recipe for reading it, and `np.dot` is a compiled loop over it that creates no Python object per element. The price: *you must stop writing the loop and express the computation as operations on whole arrays.*

The forcing problem is the random-walk ensemble from F5: $10^4$ walkers, $10^3$ steps each, $10^7$ steps in all, to measure how the mean-square displacement grows. The loop version first, measured on a small piece:

```python
import random


def walks_loop(n_walkers: int, n_steps: int, seed: int = 0) -> list[list[int]]:
    r = random.Random(seed)
    paths = []
    for _ in range(n_walkers):
        x, path = 0, []
        for _ in range(n_steps):
            x += r.choice((-1, 1))
            path.append(x)
        paths.append(path)
    return paths


def walks_array(n_walkers: int, n_steps: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    steps = rng.choice([-1, 1], size=(n_walkers, n_steps))
    return steps.cumsum(axis=1)


def walks_array_fast(n_walkers: int, n_steps: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    steps = 2 * rng.integers(0, 2, size=(n_walkers, n_steps), dtype=np.int16) - 1
    return steps.cumsum(axis=1, dtype=np.int16)


t_small = best(walks_loop, 200, 1000, repeat=2)
per_step = t_small / (200 * 1000)
print(f"loop:  200 x 1000 steps took {t_small:.3f} s = {per_step * 1e6:.2f} us per step")
print(f"       at that rate 10^7 steps would take about {per_step * 1e7:.1f} s (a projection)")
t_big = best(walks_array, 10_000, 1_000, repeat=2)
print(f"array: 10^4 x 10^3 steps took {t_big:.3f} s  (measured, all of them)")
t_fast = best(walks_array_fast, 10_000, 1_000, repeat=2)
print(f"array, cheaper dtype and generator call: {t_fast:.3f} s")
```
```text
loop:  200 x 1000 steps took 0.067 s = 0.34 us per step
       at that rate 10^7 steps would take about 3.4 s (a projection)
array: 10^4 x 10^3 steps took 0.116 s  (measured, all of them)
array, cheaper dtype and generator call: 0.064 s
```

The array version has no loop: one call draws every step, one `cumsum` along axis 1 turns steps into positions. The second version only narrows the dtype (16 bits hold any position up to $10^3$) and draws integers directly. Writing the array version, and trusting it, is what this module teaches. It rests on three ideas (§2): an array is *memory plus a recipe*, the recipe can change without touching the memory (views, broadcasting), and loops move inside whole-array operations.

## 2. Intuition

**An array is a block of memory plus a recipe.** Take `np.arange(6).reshape(2, 3)` of 8-byte integers: six numbers in one run of memory, and a recipe: shape `(2, 3)`; to find element `[i, j]`, jump `24*i + 8*j` bytes from the start. Those jump sizes are the **strides**.

```text
memory (48 bytes):   | 0 | 1 | 2 | 3 | 4 | 5 |
byte offset            0   8  16  24  32  40

a    shape (2, 3)  strides (24, 8)     a[i, j]  at 24*i + 8*j
a.T  shape (3, 2)  strides (8, 24)     a.T[j, i] at 8*j + 24*i    (the same bytes)
a[:, 1]  shape (2,)  strides (24,)  starts at byte 8   -> elements 1 and 4
```

The transpose did not move a byte; it swapped two numbers in the recipe. A slice did not move a byte either; it changed the start, the shape and the strides. So **slices are views: new recipes over the same memory.** Write through the view and the original changes. Python lists work the opposite way (`lst[1:3]` builds a new list), the most common source of NumPy bugs in people who know lists first.

**Broadcasting is a stride of zero.** To add a row `[10, 20, 30]` to each row of a $4\times3$ array, NumPy does not build four copies of the row; it reads the row with strides `(0, 8)`, so advancing the row index does not move the offset. "Stretching" a length-1 axis means giving it stride 0.

**Reductions delete axes.** `a.sum(axis=0)` for shape `(2, 3)` adds *down* the rows and returns shape `(3,)`: the axis you name disappears.

**Vectorization is moving the loop.** A `for` loop that does the same thing to each element can usually be written as an operation on whole arrays, and then the loop runs in C. The skill is the translation: find the loop's index variables and ask which axes they would be.

You can already guess the rules: slicing is cheap and aliases; selecting by a list of indices or a mask cannot be described by one stride per axis, so it must copy; shapes align from the right; reduced axes vanish. §3 states each exactly and has you build the machinery before using the library's.


## 3. Formal development

Each subsection states the semantics, then has you predict or build before it shows the library result. Eight **Build** sequences are marked; do them in order.

### 3.1 The ndarray: a buffer, a dtype, a shape, strides

An `ndarray` is four things: a **buffer** of bytes, a **dtype** saying how to read an element from those bytes, a **shape** $(n_0,\dots,n_{d-1})$, and **strides** $(s_0,\dots,s_{d-1})$, the number of bytes to jump when index $k$ increases by one. Element $(i_0,\dots,i_{d-1})$ is read from byte

$$\text{offset}(i_0,\dots,i_{d-1}) \;=\; \sum_{k=0}^{d-1} i_k\, s_k \tag{3.1}$$

counted from where this array starts in the buffer. A **C-contiguous** array (row-major: the last index changes fastest in memory) has

$$s_{d-1} = \text{itemsize},\qquad s_k = s_{k+1}\, n_{k+1}. \tag{3.2}$$

(**Fortran order**, first index fastest, is `order="F"`, LAPACK's convention; we use C order.)

```python
a = np.arange(6).reshape(2, 3)
print(a.shape, a.ndim, a.size, a.dtype, a.itemsize, a.nbytes, a.strides)   # in that order
print("C_CONTIGUOUS:", a.flags["C_CONTIGUOUS"], "| OWNDATA:", a.flags["OWNDATA"])
print("a.base:", a.base)
```
```text
(2, 3) 2 6 int64 8 48 (24, 8)
C_CONTIGUOUS: True | OWNDATA: False
a.base: [0 1 2 3 4 5]
```

`a` does not own its data: `reshape` returned a view of the 1-D array `arange` made, which is `a.base`. Views are the subject of §3.4.

**▶ Try it first.** For this `a`, predict `a.T.strides`, `a[:, 1].strides`, `a[:, ::2].strides` and `a[::-1].strides`, using (3.1) and the picture in §2. Then run.

```python
for label, v in [("a.T", a.T), ("a[:, 1]", a[:, 1]), ("a[:, ::2]", a[:, ::2]), ("a[::-1]", a[::-1])]:
    print(f"{label:10s} shape {v.shape}  strides {v.strides}")
```
```text
a.T        shape (3, 2)  strides (8, 24)
a[:, 1]    shape (2,)  strides (24,)
a[:, ::2]  shape (2, 2)  strides (24, 16)
a[::-1]    shape (2, 3)  strides (-24, 8)
```

**■ What just happened.** Every second column doubles the column stride (8 to 16). Reversing the rows *negates* the row stride: the view starts at the last row and walks backwards through the same bytes. Nothing was copied.

**Build 1 (addressing).** Before reading on, write `c_strides(shape, itemsize)` implementing (3.2) and `read_element(raw, dtype, strides, idx, start)` implementing (3.1) on the raw bytes of the buffer. Test both on every element of `a`, `a.T`, `a[::2, 1:]` and `a[::-1]`. Then compare:

```python
def c_strides(shape: Shape, itemsize: int) -> Shape:
    """Strides in bytes of a C-contiguous array, equation (3.2)."""
    strides, step = [], itemsize
    for n_k in reversed(shape):
        strides.append(step)
        step *= n_k
    return tuple(reversed(strides))


def read_element(raw: bytes, dtype: np.dtype, strides: Shape, idx: Shape, start: int = 0) -> np.generic:
    """Element idx of the array stored in bytes `raw`, equation (3.1)."""
    offset = start + sum(i * s for i, s in zip(idx, strides))
    return np.frombuffer(raw, dtype=dtype, count=1, offset=offset)[0]


raw = a.base.tobytes()                       # the 48 bytes of the buffer
assert c_strides(a.shape, a.itemsize) == a.strides == (24, 8)
for view in (a, a.T, a[::2, 1:], a[::-1]):
    start = view.ctypes.data - a.ctypes.data     # where this view begins in the buffer
    for idx in np.ndindex(view.shape):
        assert read_element(raw, a.dtype, view.strides, idx, start) == view[idx]
    print(f"shape {view.shape}, strides {view.strides}, starts at byte {start}: all elements match")
```
```text
shape (2, 3), strides (24, 8), starts at byte 0: all elements match
shape (3, 2), strides (8, 24), starts at byte 0: all elements match
shape (1, 2), strides (48, 8), starts at byte 8: all elements match
shape (2, 3), strides (-24, 8), starts at byte 24: all elements match
```

Your function reproduced NumPy's indexing on four views of one buffer: (3.1) *is* how indexing works.

### 3.2 dtypes: what a byte means, and what silently goes wrong

The dtypes you will meet: `float64` (about 15.9 decimal digits; the default), `float32` (about 7.2 digits; images, GPUs), `int64`, `int8`/`uint8` ($-128..127$ and $0..255$; pixels, and our overflow demonstration), `bool` (masks, §3.4) and `complex128` (F3's FFTs, quantum states).

**▶ Try it first.** With `x = np.array([100, 100], dtype=np.int8)`, predict `x + x`, `x.sum()` and `np.int8(100) + np.int8(100)`. The largest `int8` is 127.

```python
x = np.array([100, 100], dtype=np.int8)
print(x + x)
print(x.sum(), x.sum().dtype)
print(np.int8(100) + np.int8(100))
try:
    x + 200
except OverflowError as err:
    print("OverflowError:", err)
print(np.finfo(np.float64).eps, np.finfo(np.float32).eps, np.iinfo(np.int8).max)
```
```text
[-56 -56]
200 int64
RuntimeWarning: overflow encountered in scalar add
-56
OverflowError: Python integer 200 out of bounds for int8
2.220446049250313e-16 1.1920929e-07 127
```

**■ What just happened.** Integer arrays do arithmetic modulo $2^{8}$: $200 - 256 = -56$, with no error. Scalar arithmetic at least warns; array arithmetic is silent, and silence is the danger (a `uint8` image brightened by `+ 100` wraps to dark pixels). `sum` accumulates in `int64`, so it returns 200 while `x + x` wrapped. In NumPy 2 a Python integer is "weak": it adopts the array's dtype instead of promoting it, so `x + 200` raises rather than wrapping. Floats do not wrap; they lose digits (§3.9; `eps` above is their relative spacing at 1).

The promotion rules in one block:

```python
f32 = np.array([0.1, 0.2], dtype=np.float32)
print((f32 + 1.0).dtype, (f32 + np.float64(1.0)).dtype, (np.arange(3) + 0.5).dtype)
print(np.result_type(np.int8, np.uint8), np.result_type(np.int64, np.float32))
print(np.array([1.9, -1.9]).astype(int), np.array([300, -1]).astype(np.uint8))
```
```text
float32 float64 float64
int16 float64
[ 1 -1] [ 44 255]
```

Two arrays (or NumPy scalars) promote to the smallest dtype holding both (`int8` with `uint8` gives `int16`; `int64` with `float32`, `float64`). A Python scalar never promotes the array: `float32 + 1.0` stays `float32`. `astype` copies, truncates floats toward zero, and wraps integers that do not fit.

### 3.3 Creating arrays

```python
print(np.zeros((2, 3), dtype=int).tolist(), np.full(3, 7.5), np.eye(3)[1], np.diag([1, 2, 3])[1])
print(np.linspace(0, 1, 5), np.linspace(0, 1, 5, endpoint=False))
r = np.arange(1, 1.3, 0.1)
print(len(r), r[-1], (1.3 - 1) / 0.1)
x_ax, y_ax = np.arange(3), np.arange(2)
X, Y = np.meshgrid(x_ax, y_ax)
Xi, Yi = np.meshgrid(x_ax, y_ax, indexing="ij")
gy, gx = np.ogrid[0:2, 0:3]
print(X.shape, Xi.shape, gy.shape, gx.shape)
rng = np.random.default_rng(42)
print(rng.random(2), rng.normal(10, 2, size=2), rng.integers(0, 6, size=5), rng.permutation(5))
print(rng.choice(["a", "b", "c"], size=4), rng.choice(10, size=3, replace=False))
```
```text
[[0, 0, 0], [0, 0, 0]] [7.5 7.5 7.5] [0. 1. 0.] [0 2 0]
[0.   0.25 0.5  0.75 1.  ] [0.  0.2 0.4 0.6 0.8]
4 1.3000000000000003 3.0000000000000004
(2, 3) (3, 2) (2, 1) (1, 3)
[0.77395605 0.43887844] [11.50090239 11.88112943] [1 0 3 5 4] [1 4 2 3 0]
['b' 'b' 'b' 'a'] [7 6 8]
```

- `np.empty(shape)` allocates **without initializing**: the contents are whatever bytes were in that memory, different on different runs (so the block above does not print them). Use it only when you will overwrite every element.
- `np.arange(start, stop, step)` has length $\lceil (\text{stop}-\text{start})/\text{step}\rceil$. With a float step, $0.3/0.1 = 3.0000000000000004$ rounds *up* and the "excluded" endpoint appears (`r` has 4 elements ending at 1.3000000000000003). For floats use `np.linspace(start, stop, num)` (endpoint included; `endpoint=False` for periodic grids such as FFT sample points, F3), or an integer `arange` times the step.
- `np.meshgrid(x, y)` defaults to `indexing="xy"`, shape `(len(y), len(x))` (plotting); `indexing="ij"` gives `(len(x), len(y))` (matrices). `np.ogrid` returns the *open* grid, shapes `(2, 1)` and `(1, 3)`, which broadcast (§3.5) to the full grid without storing it.
- `rng.integers(low, high)` excludes `high`; `rng.choice(a, size, replace=True, p=None)` samples; `rng.permutation(n)` shuffles `0..n-1`. Pass `rng` into functions; old code's `np.random.seed` and `np.random.rand` share a hidden global state.


### 3.4 Indexing: views, copies, and masks

Four rules decide everything; call them (3.3).

1. **Basic indexing** (integers, slices, `None`, `...`) returns a **view**: a new shape, strides and start over the same buffer. An integer removes its axis, a slice keeps it, `None` inserts a length-1 axis, `...` means "all the axes I did not mention".
2. **Advanced indexing** (an integer list or array, or a boolean array) returns a **copy**: an arbitrary list of positions such as `[2, 0, 2]` is not "start plus a constant stride per axis", so no recipe over the old buffer describes it.
3. **Assignment** `a[key] = value` writes into `a` for every kind of key, with `value` broadcast (§3.5).
4. `.T`, `transpose`, `swapaxes`, `np.broadcast_to` and `sliding_window_view` are views; `reshape` and `ravel` are views *when the strides allow* and copies otherwise; `flatten`, `copy`, `astype`, arithmetic and `np.sort` copy.

```python
a = np.arange(12).reshape(3, 4)
print(a[1], a[1, 2], a[-1, -1], a[:, 1], a[1:, ::2], sep="\n")
print(a[:, None, 0].shape, a[..., 0].shape, a[None].shape)
```
```text
[4 5 6 7]
6
11
[1 5 9]
[[ 4  6]
 [ 8 10]]
(3, 1) (3,) (1, 3, 4)
```

**Build 2 (views versus copies).** **▶ Try it first.** For this `a`, write down which of the eight results below share memory with `a`, before running:

```python
ops = {
    "a[1]": lambda a: a[1],
    "a[:, 1]": lambda a: a[:, 1],
    "a[::2, ::2]": lambda a: a[::2, ::2],
    "a[[0, 2]]": lambda a: a[[0, 2]],
    "a[a % 2 == 0]": lambda a: a[a % 2 == 0],
    "a.T": lambda a: a.T,
    "a.reshape(4, 3)": lambda a: a.reshape(4, 3),
    "a.T.reshape(12)": lambda a: a.T.reshape(12),
}
for name, f in ops.items():
    b = f(a)
    print(f"{name:16s} shares memory: {np.shares_memory(a, b)!s:5}  owns its data: {b.flags['OWNDATA']}")
```
```text
a[1]             shares memory: True   owns its data: False
a[:, 1]          shares memory: True   owns its data: False
a[::2, ::2]      shares memory: True   owns its data: False
a[[0, 2]]        shares memory: False  owns its data: True
a[a % 2 == 0]    shares memory: False  owns its data: True
a.T              shares memory: True   owns its data: False
a.reshape(4, 3)  shares memory: True   owns its data: False
a.T.reshape(12)  shares memory: False  owns its data: False
```

**■ What just happened.** Five views, three copies. `a.T.reshape(12)` copies because the transposed layout (strides `(8, 32)`) cannot be read as one run with a single stride, while `a.reshape(4, 3)` reads the same contiguous bytes with new strides. `np.shares_memory` is exact; `np.may_share_memory` only compares memory bounds. Do not use `OWNDATA` or `b.base` as the test: the last row is a copy that does not own its data (`reshape` copied into a temporary and returned a view of *that*).

Aliasing, observed:

```python
v = a[:, 1]
v[:] = -1
print(a)                      # column 1 of a changed: v is a window onto a
a = np.arange(12).reshape(3, 4)
f = a[[0, 2]]
f[:] = -1
print(a[0])                   # unchanged: f is a copy
```
```text
[[ 0 -1  2  3]
 [ 4 -1  6  7]
 [ 8 -1 10 11]]
[0 1 2 3]
```

**Masks and positions.**

```python
mask = a % 3 == 0
print(a[mask], mask.sum(), np.count_nonzero(mask))
print(np.where(mask, a, -1))                          # three arguments: choose elementwise
print(np.nonzero(mask))                               # one argument: row and column indices of True
b = a.copy()
b[b > 8] = 0                                          # assignment through a mask
print(b[2], a[[0, 2], [1, 3]], a[np.ix_([0, 2], [1, 3])], sep="\n")
c = np.zeros(3, dtype=int)
c[[0, 0, 1]] += 1
d = np.zeros(3, dtype=int)
np.add.at(d, [0, 0, 1], 1)
print(c, d)
```
```text
[0 3 6 9] 4 4
[[ 0 -1 -1  3]
 [-1 -1  6 -1]
 [-1  9 -1 -1]]
(array([0, 0, 1, 2]), array([0, 3, 2, 1]))
[8 0 0 0]
[ 1 11]
[[ 1  3]
 [ 9 11]]
[1 1 0] [2 1 0]
```

`True` counts as 1, so `mask.sum()` counts. `a[[0, 2], [1, 3]]` pairs the lists into two points, `(0,1)` and `(2,3)`; `np.ix_` builds the 2×2 block of rows `{0,2}` × columns `{1,3}`. In `c[[0, 0, 1]] += 1` all reads happen first, so the repeated index increments once; `np.add.at` accumulates.

**When you need `.copy()`.** (i) You will modify the slice (or the original) and the other must not change. (ii) The slice is small but its base is huge: a view keeps the whole base alive. (iii) You return an array and the caller must not be able to change your internal state. An unintended copy only costs time; the aliasing bug costs correctness and surfaces far from its cause (Example 4.2).

### 3.5 Broadcasting

**Broadcasting rule (3.4).** To combine arrays of shapes $s$ and $t$:

1. Right-align the shapes; pad the shorter on the left with 1s.
2. Compare axis by axis: the lengths must be equal, or one of them must be 1. Otherwise `ValueError`.
3. The result's length on each axis is the larger one (a length-1 axis is stretched; if the other length is 0, the result is 0).

Stretching is virtual: a length-1 axis is read with **stride 0**, so no data is copied. As an index formula, element $(i,j,k)$ of `a + b` is `a[i, j, k] + b[i, j, k]`, with index 0 used wherever that array's axis has length 1 or is missing.

**▶ Try it first.** By hand, give the broadcast shape (or "error") for each pair: `(3,4)+(4,)`, `(3,4)+(3,)`, `(3,1)+(4,)`, `(2,1,3)+(4,1)`, `(8,1,6,1)+(7,1,5)`.

**Build 3 (broadcasting by hand).** Write `bshape(s, t)` from rules 1 to 3 and `my_add(a, b)` that uses `np.broadcast_to` and checks against `a + b`. Compare:

```python
def bshape(s: Shape, t: Shape) -> Shape:
    """Broadcast shape of shapes s and t by rules 1-3 of (3.4); ValueError if incompatible."""
    n = max(len(s), len(t))
    s, t = (1,) * (n - len(s)) + tuple(s), (1,) * (n - len(t)) + tuple(t)
    out = []
    for p, q in zip(s, t):
        if p != q and p != 1 and q != 1:
            raise ValueError(f"shapes {s} and {t} do not broadcast")
        out.append(q if p == 1 else p)
    return tuple(out)


def my_add(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    shape = bshape(a.shape, b.shape)
    return np.broadcast_to(a, shape) + np.broadcast_to(b, shape)


pairs = [((3, 4), (4,)), ((3, 4), (3,)), ((3, 1), (4,)), ((2, 1, 3), (4, 1)), ((8, 1, 6, 1), (7, 1, 5))]
for s, t in pairs:
    try:
        mine = bshape(s, t)
    except ValueError:
        mine = "error"
    try:
        real = np.broadcast_shapes(s, t)
    except ValueError:
        real = "error"
    print(f"{s!s:13} + {t!s:10} -> {mine!s:13} NumPy: {real!s:13} {'agree' if mine == real else 'DISAGREE'}")
    if mine != "error":
        p, q = np.ones(s), np.arange(float(np.prod(t))).reshape(t)
        assert np.array_equal(my_add(p, q), p + q)
try:
    np.ones((3, 4)) + np.ones(3)
except ValueError as err:
    print(err)
row = np.arange(3)
B = np.broadcast_to(row, (4, 3))
print(B.strides, B.flags["WRITEABLE"], np.shares_memory(B, row))
```
```text
(3, 4)        + (4,)       -> (3, 4)        NumPy: (3, 4)        agree
(3, 4)        + (3,)       -> error         NumPy: error         agree
(3, 1)        + (4,)       -> (3, 4)        NumPy: (3, 4)        agree
(2, 1, 3)     + (4, 1)     -> (2, 4, 3)     NumPy: (2, 4, 3)     agree
(8, 1, 6, 1)  + (7, 1, 5)  -> (8, 7, 6, 5)  NumPy: (8, 7, 6, 5)  agree
operands could not be broadcast together with shapes (3,4) (3,) 
(0, 8) False True
```

**■ What just happened.** Your five answers and NumPy's agree. The error message *is* the diagnosis: it prints both shapes, and the mismatch sits on a trailing axis (4 against 3). `B` has strides `(0, 8)`: four "copies" of one row cost nothing, and the result is read-only because a write would change all four at once.

The patterns you will use daily:

```python
x = np.array([1.0, 2.0, 4.0])
print((x[:, None] * x[None, :]).shape, (x[:, None] - x[None, :])[0])      # outer product; all differences x_0 - x_j
X = np.array([[1.0, 10.0], [2.0, 20.0], [6.0, 60.0]])
print(X - X.mean(axis=0))                       # centre each column: (3,2) - (2,)
print(X - X.mean(axis=1, keepdims=True))        # centre each row: (3,2) - (3,1)
```
```text
(3, 3) [ 0. -1. -3.]
[[ -2. -20.]
 [ -1. -10.]
 [  3.  30.]]
[[ -4.5   4.5]
 [ -9.    9. ]
 [-27.   27. ]]
```

`X.mean(axis=0)` has shape `(2,)`, which right-aligns against the columns. `X.mean(axis=1)` has shape `(3,)` and would align against the 2 columns and fail (the error above); `keepdims=True` keeps the length-1 axis. Rule of thumb: **to stretch an operand along axis $k$, it needs a length-1 axis at position $k$.**

**Build 4 (pairwise distances three ways).** For $n$ points $\mathbf{x}_i\in\mathbb{R}^d$ we want $D_{ij}=\|\mathbf{x}_i-\mathbf{x}_j\|$. **▶ Try it first.** Write the double loop and the broadcast form `X[:, None, :] - X[None, :, :]` (shape $(n,n,d)$; square, sum over the last axis, root), and predict the order of the three timings and the size of the $(n,n,d)$ array for $n=300$, $d=3$. Then compare with `scipy.spatial.distance.cdist`:

```python
from scipy.spatial.distance import cdist


def dist_loops(X: np.ndarray) -> np.ndarray:
    n = len(X)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = np.sqrt(((X[i] - X[j]) ** 2).sum())
    return D


def dist_broadcast(X: np.ndarray) -> np.ndarray:
    diff = X[:, None, :] - X[None, :, :]          # (n, n, d): the memory cost is n*n*d floats
    return np.sqrt((diff ** 2).sum(axis=-1))


X = np.random.default_rng(5).standard_normal((300, 3))
D1, D2, D3 = dist_loops(X), dist_broadcast(X), cdist(X, X)
assert np.allclose(D1, D2) and np.allclose(D2, D3)
for name, f in [("double loop", dist_loops), ("broadcasting", dist_broadcast), ("cdist", lambda Z: cdist(Z, Z))]:
    print(f"{name:13s} {best(f, X) * 1e3:9.2f} ms")
print(f"the (n, n, d) intermediate is {300 * 300 * 3 * 8 / 1e6:.1f} MB; at n = 10,000 it would be {10_000 ** 2 * 3 * 8 / 1e9:.1f} GB")
```
```text
double loop      195.74 ms
broadcasting       3.07 ms
cdist              0.19 ms
the (n, n, d) intermediate is 2.2 MB; at n = 10,000 it would be 2.4 GB
```

Three answers, equal to rounding. The loop pays Python overhead $n^2$ times. Broadcasting pays $O(n^2 d)$ memory instead; `cdist` stores nothing extra. When the intermediate no longer fits, chunk the work or call a purpose-built routine (Exercise B3: the Gram-matrix alternative, and where it loses digits).


### 3.6 Ufuncs and reductions: the "which axis" rule

A **ufunc** (universal function) applies a scalar function elementwise, with broadcasting, in a compiled loop: `+`, `*`, `np.sqrt`, `np.exp`, `np.maximum`, comparisons. Two keywords and four methods are worth knowing.

```python
x = np.array([1.0, 4.0, 9.0, -1.0])
out = np.zeros_like(x)
np.sqrt(x, out=out, where=x >= 0)                 # write into `out`; skip entries where the mask is False
print(out)
v = np.array([1, 2, 3, 4])
print(np.add.reduce(v), np.add.accumulate(v), np.multiply.outer(v, v)[1], np.maximum.accumulate([3, 1, 4, 1, 5]))
```
```text
[1. 2. 3. 0.]
10 [ 1  3  6 10] [2 4 6 8] [3 3 4 4 5]
```

`out=` reuses memory instead of allocating a temporary; `where=` computes only where the mask is true and leaves `out` untouched elsewhere (hence `zeros_like`; without `out`, skipped entries are uninitialized). `reduce` collapses an axis with the operation (`np.add.reduce` is `sum`), `accumulate` keeps the running result (`cumsum`), `outer` applies it to all pairs, `at` is the unbuffered in-place form of §3.4.

**Reductions** (`sum`, `prod`, `mean`, `std`, `var`, `min`, `max`, `argmin`, `argmax`, `any`, `all`) take `axis`. The rule that removes the confusion is (3.5):

> **The result loses the axis you name.** `f(a, axis=k)` returns shape `a.shape` with entry $k$ deleted (entry replaced by 1 if `keepdims=True`). `axis=None` (the default) collapses everything to a scalar; a tuple of axes deletes each. `cumsum`, `cumprod`, `sort` and `argsort` *keep* the axis: they work along it without removing it.

"Down the rows" and "along the columns" are the phrasings people mix up; delete-the-axis is unambiguous. For shape `(n, m)`: `axis=0` gives `(m,)`, one number per column; `axis=1` gives `(n,)`, one per row.

**Build 5 (the axis drill).** **▶ Try it first.** For `t = np.arange(24).reshape(2, 3, 4)`, write the shape of each result before running.

```python
t = np.arange(24).reshape(2, 3, 4)
drill = ["t.sum(axis=0)", "t.sum(axis=1)", "t.sum(axis=2)", "t.sum(axis=(0, 2))", "t.sum(axis=1, keepdims=True)",
         "t.mean(axis=-1)", "t.argmax(axis=0)", "t.cumsum(axis=1)", "(t > 10).any(axis=2)", "t.std(axis=0, ddof=1)"]
for expr in drill:
    print(f"{expr:30s} -> {np.shape(eval(expr))}")
print(t.sum(), t.max(), t.sum(axis=None) == t.sum())
```
```text
t.sum(axis=0)                  -> (3, 4)
t.sum(axis=1)                  -> (2, 4)
t.sum(axis=2)                  -> (2, 3)
t.sum(axis=(0, 2))             -> (3,)
t.sum(axis=1, keepdims=True)   -> (2, 1, 4)
t.mean(axis=-1)                -> (2, 3)
t.argmax(axis=0)               -> (3, 4)
t.cumsum(axis=1)               -> (2, 3, 4)
(t > 10).any(axis=2)           -> (2, 3)
t.std(axis=0, ddof=1)          -> (3, 4)
276 23 True
```

**■ What just happened.** Delete the named axis from `(2, 3, 4)` each time; `keepdims` leaves a 1; `cumsum` deletes nothing; `argmax(axis=0)` returns *positions* along the deleted axis; `axis=-1` counts from the end. When a result has the wrong shape, print `.shape` of every intermediate.

**`ddof`.** `np.var` and `np.std` divide by $N - \text{ddof}$, and `ddof` defaults to **0** (the population formula). Python's `statistics.stdev` divides by $N-1$. For a sample you want the unbiased variance, `ddof=1` (ST-1 derives why):

```python
import statistics

data = [1.0, 2.0, 3.0, 4.0]
print(np.std(data), np.std(data, ddof=1), statistics.stdev(data), statistics.pstdev(data))
```
```text
1.118033988749895 1.2909944487358056 1.2909944487358056 1.118033988749895
```

### 3.7 Sorting, searching, reshaping

```python
v = np.array([30, 10, 20, 10])
print(np.sort(v), np.argsort(v, kind="stable"), np.unique(v, return_counts=True))
s = np.sort(v)
print(np.searchsorted(s, [5, 10, 15, 40]), np.searchsorted(s, 10, side="right"))
print(np.partition(v, 1)[:2])
a = np.arange(12)
m = a.reshape(3, -1)
print(m.shape, np.shares_memory(a, m), m.T.shape)
print(np.concatenate([m, m], axis=0).shape, np.hstack([m, m]).shape, np.stack([m, m]).shape,
      [p.shape for p in np.split(m, 2, axis=1)])
t = np.arange(24).reshape(2, 3, 4)
print(t.transpose(2, 0, 1).shape, t.swapaxes(0, 2).shape, m.ravel().flags["OWNDATA"], np.shares_memory(m, m.T.ravel()))
```
```text
[10 10 20 30] [1 3 2 0] (array([10, 20, 30]), array([2, 1, 1]))
[0 0 2 4] 2
[10 10]
(3, 4) True (4, 3)
(6, 4) (3, 8) (2, 3, 4) [(3, 2), (3, 2)]
(4, 2, 3) (4, 3, 2) False False
```

- `np.sort(v)` returns a sorted copy, `v.sort()` sorts in place; both sort along the last axis. `np.argsort` returns the indices that would sort (Example 4.4); `kind="stable"` keeps ties in input order.
- `np.searchsorted(s, q)` binary-searches where each $q$ would be inserted to keep sorted `s` sorted (`side="left"` before equal entries, `"right"` after). It silently returns garbage if `s` is not sorted.
- `np.unique` returns sorted distinct values (`return_counts`, `return_inverse`). `np.partition(v, k)` puts the $k$ smallest values first in $O(N)$ time without sorting the rest.
- `reshape` takes one `-1` ("infer this length"). `concatenate` joins along an *existing* axis, `stack` creates a *new* one, `hstack`/`vstack` are shorthand for axis 1/0, `split` returns a list of views. `transpose(2, 0, 1)` lists the old axes in their new order. `ravel` is a view when it can be (`m.ravel()` does not own its data) and a copy when it must (`m.T.ravel()` shares nothing with `m`, the last two printed values); `flatten` always copies.


### 3.8 Linear algebra and `einsum`

`@` is matrix multiplication; `*` is elementwise and never matrix multiplication (the most common first-month bug, silent when shapes happen to agree). `@` on two 1-D arrays is the dot product; on 2-D arrays the matrix product; with 3-D or higher it treats leading axes as a batch and broadcasts them. `np.dot` agrees for 1-D and 2-D and does something different (no batching) above that; use `@`.

**Build 6 (matrix multiply, three ways).** Write $C_{ik}=\sum_j A_{ij}B_{jk}$ as three nested loops on lists of lists, then call `np.dot`, then `@`, and time all three for $n = 50, 100, 200$. **▶ Try it first:** guess the loops-to-`@` ratio at $n=200$ (10, 100, 1000, or more?).

```python
def matmul_loops(A: list[list[float]], B: list[list[float]]) -> list[list[float]]:
    """C[i][k] = sum_j A[i][j] * B[j][k] on lists of lists."""
    n, m, p = len(A), len(B), len(B[0])
    C = [[0.0] * p for _ in range(n)]
    for i in range(n):
        for k in range(p):
            s = 0.0
            for j in range(m):
                s += A[i][j] * B[j][k]
            C[i][k] = s
    return C


rng = np.random.default_rng(1)
print(f"{'n':>5} {'loops (s)':>10} {'np.dot (s)':>11} {'@ (s)':>10} {'loops / @':>10}")
for n in (50, 100, 200):
    A, B = rng.standard_normal((n, n)), rng.standard_normal((n, n))
    Al, Bl = A.tolist(), B.tolist()
    t0 = time.perf_counter()
    C = matmul_loops(Al, Bl)
    t_loops = time.perf_counter() - t0
    assert np.allclose(C, A @ B) and np.allclose(np.dot(A, B), A @ B)
    t_dot, t_at = best(np.dot, A, B, repeat=5), best(lambda: A @ B, repeat=5)
    print(f"{n:5d} {t_loops:10.4f} {t_dot:11.6f} {t_at:10.6f} {t_loops / t_at:10,.0f}")
```
```text
    n  loops (s)  np.dot (s)      @ (s)  loops / @
   50     0.0066    0.000006   0.000006      1,035
  100     0.0519    0.000039   0.000038      1,361
  200     0.4408    0.000135   0.000114      3,884
```

The flop count is $2n^3$ either way. The gap grows with $n$ because the library routine (BLAS) blocks the computation to stay in cache, uses SIMD instructions and may use several cores, while the loop pays an interpreter step per multiply. If a formula is a sum over a repeated index, it is a `@` or an `einsum`, not a loop.

| Task | Call | Destination |
|---|---|---|
| solve $A\mathbf{x}=\mathbf{b}$ | `np.linalg.solve(A, b)` | LA-2 (elimination, LU) |
| least squares $\min\|M\mathbf{c}-\mathbf{y}\|$ | `np.linalg.lstsq(M, y, rcond=None)` | LA-5 (projections) |
| orthogonal factorization | `np.linalg.qr(M)` | LA-5 (Gram–Schmidt) |
| determinant | `np.linalg.det(A)` | LA-6 |
| eigenpairs | `np.linalg.eig` (general), `np.linalg.eigh` (symmetric/Hermitian: real eigenvalues ascending, orthonormal vectors) | LA-7, LA-8 |
| SVD | `np.linalg.svd(M, full_matrices=False)` | LA-9 |
| norms | `np.linalg.norm(x, ord=None, axis=None)`: 2-norm for vectors, Frobenius for matrices by default | LA-1, LA-9 |
| inverse | `np.linalg.inv(A)`: almost never (below) | LA-2, LA-13 |

`project.py` checks each call against its defining property (for example `S @ V == V * w` for `eigh`). Note `V * w`: a `(4,4)` array times a `(4,)` array scales *column* $i$ by $w_i$ (`w` aligns with the last axis), which is $V\,\mathrm{diag}(w)$ without building the diagonal.

**Why not `inv`.** Computing $A^{-1}$ and multiplying costs several times one solve (measured below) and is less accurate. We build matrices with a chosen condition number $\kappa$ (the ratio of largest to smallest singular value, LA-9) and compare:

```python
def make_matrix(n: int, cond: float, rng: np.random.Generator) -> np.ndarray:
    U, _ = np.linalg.qr(rng.standard_normal((n, n)))
    V, _ = np.linalg.qr(rng.standard_normal((n, n)))
    return (U * np.logspace(0, -np.log10(cond), n)) @ V.T      # singular values from 1 down to 1/cond


rng = np.random.default_rng(0)
for cond in (1e2, 1e8, 1e12):
    A = make_matrix(200, cond, rng)
    x_true = rng.standard_normal(200)
    b = A @ x_true
    x_s, x_i = np.linalg.solve(A, b), np.linalg.inv(A) @ b
    res = lambda z: np.linalg.norm(A @ z - b) / np.linalg.norm(b)
    err = lambda z: np.linalg.norm(z - x_true) / np.linalg.norm(x_true)
    print(f"kappa {cond:6.0e}: residual solve {res(x_s):.0e} inv {res(x_i):.0e} | error solve {err(x_s):.0e} inv {err(x_i):.0e}")
A, b = rng.standard_normal((1000, 1000)), rng.standard_normal(1000)
print(f"n = 1000: solve {best(np.linalg.solve, A, b) * 1e3:.0f} ms, inv then @ {best(lambda: np.linalg.inv(A) @ b) * 1e3:.0f} ms")
```
```text
kappa  1e+02: residual solve 2e-15 inv 1e-14 | error solve 1e-14 inv 5e-14
kappa  1e+08: residual solve 7e-16 inv 2e-09 | error solve 1e-09 inv 1e-08
kappa  1e+12: residual solve 5e-16 inv 1e-05 | error solve 6e-06 inv 3e-04
n = 1000: solve 13 ms, inv then @ 41 ms
```

`solve` keeps the *residual* $\|A\hat x-b\|/\|b\|$ at machine precision whatever $\kappa$ is (LU with partial pivoting is backward stable, LA-13), so its error is only what the conditioning allows, about $\kappa\,\epsilon_{\text{mach}}$. `inv(A) @ b` loses the small residual and more digits, and costs about 3 times as much here. Use `inv` only when you need the inverse itself.

**`einsum`.** Einstein summation writes a tensor expression as index letters. Rules (3.6): (i) a letter appearing in several operands means "multiply along it"; (ii) a letter absent from the output is **summed over**; (iii) a letter repeated within one operand takes the diagonal; (iv) the output indices appear in the order written after `->`. Five examples against their named counterparts:

```python
rng = np.random.default_rng(7)
M3, A, B = rng.standard_normal((3, 3)), rng.standard_normal((3, 4)), rng.standard_normal((4, 2))
S, T = rng.standard_normal((5, 3, 4)), rng.standard_normal((5, 4, 2))
u, w, Q, W = rng.standard_normal(3), rng.standard_normal(4), rng.standard_normal((3, 4, 5)), rng.standard_normal((4, 5))
cases = [
    ("trace          ii->",         np.einsum("ii->", M3),                 np.trace(M3)),
    ("matmul         ij,jk->ik",    np.einsum("ij,jk->ik", A, B),          A @ B),
    ("batch matmul   bij,bjk->bik", np.einsum("bij,bjk->bik", S, T),       S @ T),
    ("outer          i,j->ij",      np.einsum("i,j->ij", u, w),            np.multiply.outer(u, w)),
    ("contraction    ijk,jk->i",    np.einsum("ijk,jk->i", Q, W),          (Q * W).sum(axis=(1, 2))),
]
for label, got, want in cases:
    assert got.shape == np.shape(want) and np.allclose(got, want)
    print(f"{label:28s} result shape {np.shape(got)}")
```
```text
trace          ii->          result shape ()
matmul         ij,jk->ik     result shape (3, 2)
batch matmul   bij,bjk->bik  result shape (5, 3, 2)
outer          i,j->ij       result shape (3, 4)
contraction    ijk,jk->i     result shape (3,)
```

The last case contracts two indices at once, $\sum_{j,k}Q_{ijk}W_{jk}$, which has no one-word NumPy name; this is the operation LA-12's tensors are made of, and ML-9's attention scores are `einsum("bqd,bkd->bqk", ...)`. For matrix products prefer `@` (it calls BLAS); `optimize=True` picks a cheaper pairing order for three or more operands.

### 3.9 Floating point in arrays

Arrays are float64 unless you say otherwise, so PY-1's floating-point section holds elementwise. Two NumPy-specific points.

**Compare with tolerance.** `np.isclose(a, b)` tests $|a-b|\le \text{atol} + \text{rtol}\,|b|$ with defaults `rtol=1e-5, atol=1e-8`; `np.allclose` is `isclose(...).all()`. Note the asymmetry (the tolerance scales with `b`, the second argument) and the absolute floor:

```python
print(0.1 + 0.2 == 0.3, np.isclose(0.1 + 0.2, 0.3), np.isclose(1e-10, 2e-10), np.isclose(1e-10, 2e-10, atol=0))
```
```text
False True True False
```

Numbers differing by a factor of 2 are "close" because both lie below `atol`; pass `atol=0` for small quantities.

**Summation order matters.** Left-to-right addition of $N$ floats has a worst-case error bound growing like $N\epsilon$; a balanced tree, like $\epsilon\log_2 N$ (Higham, *Accuracy and Stability of Numerical Algorithms*, ch. 4; LA-13). `np.sum` is pairwise; `np.cumsum` and a Python loop are sequential; `math.fsum` is exact to rounding; Python's `sum` of floats has been compensated since 3.12. In float32 the difference shows:

```python
import math

x = np.random.default_rng(11).random(10**6).astype(np.float32)
exact = math.fsum(x.astype(np.float64))                 # the true sum of the stored float32 values, rounded once
rel = lambda s: abs(float(s) - exact) / exact
print(f"np.sum (pairwise, float32):      {rel(np.sum(x)):.1e}")
print(f"np.cumsum[-1] (sequential, f32): {rel(np.cumsum(x)[-1]):.1e}")
print(f"np.sum(..., dtype=float64):      {rel(np.sum(x, dtype=np.float64)):.1e}")
print(f"Python sum of the same values:   {rel(sum(x.tolist())):.1e}")
```
```text
np.sum (pairwise, float32):      4.2e-09
np.cumsum[-1] (sequential, f32): 8.5e-06
np.sum(..., dtype=float64):      1.2e-16
Python sum of the same values:   0.0e+00
```

Pairwise summation is about 2000 times more accurate than sequential accumulation here, at no cost. Accumulating in float64 (`dtype=`) is the other lever. Vectorizing a loop changes the summation order, so results differ in the last digits: assert with `isclose`, not `==`.


### 3.10 Vectorization in practice, and when a loop is fine

**Build 7 (the random-walk ensemble).** §1 gave the loop version; the array version draws all steps and takes `cumsum` along axis 1 (position is the running sum of steps along time). The physics check is F5's. With $x_t=\sum_{s=1}^{t}\xi_s$ and independent steps $\xi_s=\pm1$ of mean 0 and variance 1,

$$\langle x_t^2\rangle=\sum_{s=1}^{t}\sum_{r=1}^{t}\langle\xi_s\xi_r\rangle=\sum_{s=1}^{t}\langle\xi_s^2\rangle=t, \tag{3.7}$$

because the cross terms $\langle\xi_s\xi_r\rangle=\langle\xi_s\rangle\langle\xi_r\rangle=0$ for $s\ne r$. The mean-square displacement (MSD) grows linearly, so a log-log fit has slope 1 (**▶ Try it first:** predict MSD at $t=1000$ and the fitted slope). Over $10^4$ walkers each MSD value has relative standard error about $\sqrt{2/10^4}\approx1.4\%$ (for large $t$, $x_t$ is nearly normal with variance $t$, and $\mathrm{Var}(g^2)=2\sigma^4$ for $g\sim N(0,\sigma^2)$, so $\mathrm{Var}(x_t^2)\approx2t^2$):

```python
paths = walks_array_fast(10_000, 1_000, seed=2)            # shape (walkers, time)
t = np.arange(1, 1001)
msd = (paths.astype(np.float64) ** 2).mean(axis=0)         # average over walkers: axis 0 disappears
slope, intercept = np.polyfit(np.log(t), np.log(msd), 1)
print(paths.shape, msd.shape, f"MSD(t=1) = {msd[0]:.3f}, MSD(t=1000) = {msd[-1]:.1f}")
print(f"log-log slope = {slope:.4f}, intercept = {intercept:.4f} (theory: 1 and 0)")
```
```text
(10000, 1000) (1000,) MSD(t=1) = 1.000, MSD(t=1000) = 992.8
log-log slope = 0.9965, intercept = 0.0239 (theory: 1 and 0)
```

`paths` is (walker, time), so averaging over walkers is `axis=0` and leaves a vector indexed by time.

**Three more rules:**

```python
def relu_scalar(v: float) -> float:
    return v if v > 0 else 0.0


x = np.random.default_rng(0).standard_normal(10**6)
print(f"np.vectorize:       {best(np.vectorize(relu_scalar), x, repeat=2) * 1e3:8.1f} ms")
print(f"list comprehension: {best(lambda: [relu_scalar(v) for v in x.tolist()], repeat=2) * 1e3:8.1f} ms")
print(f"np.maximum(x, 0):   {best(np.maximum, x, 0.0) * 1e3:8.3f} ms")

X = np.random.default_rng(1).random((100_000, 10))


def normalize_loop(X: np.ndarray) -> np.ndarray:
    out = np.empty_like(X)
    for i in range(len(X)):
        out[i] = X[i] / X[i].sum()
    return out


assert np.allclose(normalize_loop(X), X / X.sum(axis=1, keepdims=True))
print(f"row-normalise 10^5 rows: loop {best(normalize_loop, X, repeat=1) * 1e3:.0f} ms, broadcast {best(lambda: X / X.sum(axis=1, keepdims=True)) * 1e3:.1f} ms")
```
```text
np.vectorize:          127.9 ms
list comprehension:    107.1 ms
np.maximum(x, 0):      1.757 ms
row-normalise 10^5 rows: loop 206 ms, broadcast 5.1 ms
```

1. **`np.vectorize` is a convenience, not a speedup:** it calls the Python function once per element. The fast form of a branch is `np.where(x > 0, x, 0.0)` or a ufunc like `np.maximum`.
2. **Preallocate, or collect in a list and convert once;** `np.append` copies the whole array on every call, so $n$ appends move $O(n^2)$ bytes (`project.py` times it).
3. **Do not loop over rows or elements** to apply the same arithmetic; a loop body indexing `X[i]` is the signal. Reduce along an axis and broadcast back.

**When a loop is fine.** Loop over the *short* axis and vectorize over the long one. (i) Step $t+1$ depends on step $t$ (time integration, the Mandelbrot iteration of Example 4.6). (ii) The loop runs a few dozen times and each pass is a large array operation, so interpreter overhead vanishes. (iii) The vectorized form needs an intermediate too big for memory (Build 4's $(n,n,d)$ array): loop over blocks of rows. (iv) The data are small, the code runs once, and clarity wins. A hot loop that cannot be vectorized is PY-10's (profiling, Numba).

### 3.11 Windows, strides, and convolution

A moving average over $w$ samples is a $(N-w+1,\,w)$ array whose row $i$ is `x[i:i+w]`. Row $i+1$ starts one element later, so the rows overlap in memory: a strided view with shape `(N-w+1, w)` and strides `(s, s)`.

```python
from numpy.lib.stride_tricks import as_strided, sliding_window_view

x = np.arange(8.0)
Wd = sliding_window_view(x, 3)
print(Wd.shape, Wd.strides, np.shares_memory(Wd, x), Wd.flags["WRITEABLE"])
print(Wd.mean(axis=-1), np.convolve(x, np.ones(3) / 3, mode="valid"))
assert np.array_equal(Wd, as_strided(x, shape=(6, 3), strides=(8, 8)))
```
```text
(6, 3) (8, 8) True False
[1. 2. 3. 4. 5. 6.] [1. 2. 3. 4. 5. 6.]
```

`sliding_window_view` checks its arguments and returns a read-only view. `as_strided` builds the same thing unchecked: a wrong `shape` or `strides` reads whatever bytes follow the array. Prefer `sliding_window_view`, and never write into overlapping views.

**Build 8 (2D convolution).** For an image $I$ and an $h\times w$ kernel $K$, the "valid" convolution is

$$(I*K)_{ij}=\sum_{p=0}^{h-1}\sum_{q=0}^{w-1} I_{i+p,\;j+q}\;K_{h-1-p,\;w-1-q}, \tag{3.8}$$

with output shape $(H-h+1,\;W-w+1)$. The kernel is *flipped*; without the flip it is cross-correlation, which deep-learning libraries call "convolution" (ML-8). **▶ Try it first:** what is the output shape for a $200\times200$ image and a $5\times5$ kernel? Write the two-loop version, then `sliding_window_view` plus `einsum`, then compare with SciPy:

```python
from scipy.signal import convolve2d


def conv2d_loops(I: np.ndarray, K: np.ndarray) -> np.ndarray:
    h, w = K.shape
    Kf = K[::-1, ::-1]
    out = np.empty((I.shape[0] - h + 1, I.shape[1] - w + 1))
    for i in range(out.shape[0]):
        for j in range(out.shape[1]):
            out[i, j] = (I[i:i + h, j:j + w] * Kf).sum()
    return out


def conv2d_windows(I: np.ndarray, K: np.ndarray) -> np.ndarray:
    win = sliding_window_view(I, K.shape)                  # shape (H-h+1, W-w+1, h, w), no copy
    return np.einsum("ijpq,pq->ij", win, K[::-1, ::-1])


rng = np.random.default_rng(2)
I, K = rng.standard_normal((200, 200)), rng.standard_normal((5, 5))
ref = convolve2d(I, K, mode="valid")
for name, f in [("two loops", conv2d_loops), ("windows + einsum", conv2d_windows), ("scipy convolve2d", lambda I, K: convolve2d(I, K, mode="valid"))]:
    assert np.allclose(f(I, K), ref)
    print(f"{name:18s} {best(f, I, K, repeat=2) * 1e3:8.2f} ms   output shape {f(I, K).shape}")
```
```text
two loops             93.92 ms   output shape (196, 196)
windows + einsum       2.70 ms   output shape (196, 196)
scipy convolve2d       5.26 ms   output shape (196, 196)
```

All three agree to rounding; the windows version does the loop's arithmetic with both loops moved into `einsum`. (For large kernels SciPy switches to FFT methods, F3.)

### 3.12 Saving and loading

`np.save("a.npy", arr)` and `np.load` store dtype, shape and raw bytes: exact and fast. `np.savez("bundle.npz", arr=arr, ids=ids)` zips several named arrays (`np.load` of it returns a lazy mapping). `np.savetxt(path, arr, delimiter=",")` and `np.loadtxt` handle text for humans and other programs; the default `%.18e` format round-trips `float64` exactly, `fmt="%.3f"` does not. The project asserts both round trips.


## 4. Worked examples

### Example 4.1: the layout and strides of a 3×4 array and its transpose

By (3.2), a C-contiguous `(3, 4)` array of 8-byte floats has strides $(4\cdot8,\,8)=(32,8)$; its transpose has shape `(4, 3)` and strides $(8,32)$.

```text
bytes   0   8  16  24 | 32  40  48  56 | 64  72  80  88
value   0   1   2   3 |  4   5   6   7 |  8   9  10  11
a[i, j]    at 32*i + 8*j          a.T[j, i]  at 8*j + 32*i
```

```python
a = np.arange(12, dtype=np.float64).reshape(3, 4)
print(a.shape, a.strides, "|", a.T.shape, a.T.strides)
print("a   C/F-contiguous:", a.flags["C_CONTIGUOUS"], a.flags["F_CONTIGUOUS"], "| a.T:", a.T.flags["C_CONTIGUOUS"], a.T.flags["F_CONTIGUOUS"])
raw = a.tobytes()
for i, j in [(0, 0), (1, 2), (2, 3)]:
    off = i * a.strides[0] + j * a.strides[1]
    assert off == j * a.T.strides[0] + i * a.T.strides[1]
    assert np.frombuffer(raw, dtype=np.float64, count=1, offset=off)[0] == a[i, j] == a.T[j, i]
    print(f"a[{i},{j}] = a.T[{j},{i}] = {a[i, j]:4.1f} at byte {off}")
print(np.shares_memory(a, a.T), a.T.reshape(12).strides, np.shares_memory(a, a.T.reshape(12)))
```
```text
(3, 4) (32, 8) | (4, 3) (8, 32)
a   C/F-contiguous: True False | a.T: False True
a[0,0] = a.T[0,0] =  0.0 at byte 0
a[1,2] = a.T[2,1] =  6.0 at byte 48
a[2,3] = a.T[3,2] = 11.0 at byte 88
True (8,) False
```

`a[2, 3]` and `a.T[3, 2]` are both byte $2\cdot32+3\cdot8=88$, the last element: the transpose reads the same bytes in another order, so `a.T` is free. The cost comes later: walking along its rows jumps 32 bytes at a time, and `reshape(12)` must copy (strides `(8,)`).

### Example 4.2: a view mutation and its copy fix

```python
def center_columns(X: np.ndarray) -> np.ndarray:
    """Intended: return X with each column's mean subtracted."""
    X -= X.mean(axis=0)               # in place: this edits the caller's array
    return X


data = np.array([[1.0, 2.0], [3.0, 6.0], [5.0, 10.0]])
train = data[:2]                      # a view of rows 0 and 1
centered = center_columns(train)
print(data)
print(np.shares_memory(data, centered), centered is train)


data = np.array([[1.0, 2.0], [3.0, 6.0], [5.0, 10.0]])
print((data[:2] - data[:2].mean(axis=0))[0], data[0])      # builds a new array instead
```
```text
[[-1. -2.]
 [ 1.  2.]
 [ 5. 10.]]
True True
[-1. -2.] [1. 2.]
```

`train` is a window onto rows 0 and 1 of `data`, and `X -= ...` modifies it in place, so the caller's first two rows were overwritten and the third was not. Nothing raised an error; the data set is half centered, and the wrong number shows up later in someone else's model fit. The fix is a new array (`X - X.mean(axis=0)`) or `X = X.copy()` first. The rule: a function either documents that it modifies its argument or never uses `-=`, `X[...] =` or `.sort()` on it.

### Example 4.3: broadcasting shape algebra for `(3, 1, 5) + (4, 1)`

```text
a:       3  1  5
b:   (-) 4  1       pad with a 1 on the left:   1  4  1
result:  3  4  5    (rule 3: the larger length on each axis)
```

Axis 2: 5 against 1 gives 5. Axis 1: 1 against 4 gives 4. Axis 0: 3 against the padded 1 gives 3. Element $(i,j,k)$ is `a[i, 0, k] + b[j, 0]`.

```python
a = np.arange(15).reshape(3, 1, 5)
b = np.arange(4).reshape(4, 1) * 100
r = a + b
print(np.broadcast_shapes(a.shape, b.shape), r.shape, r[2, 3, 4], a[2, 0, 4] + b[3, 0])
A, B = np.broadcast_arrays(a, b)
print(A.strides, B.strides)
try:
    a + np.arange(4)
except ValueError as err:
    print(err)
print((a + np.arange(4)[:, None]).shape)
```
```text
(3, 4, 5) (3, 4, 5) 314 314
(40, 0, 8) (0, 8, 0)
operands could not be broadcast together with shapes (3,1,5) (4,) 
(3, 4, 5)
```

The stretched axes carry stride 0: `a` along axis 1, `b` along axes 0 and 2. The variant `(3,1,5) + (4,)` aligns 4 against 5 and fails, as the message says; `[:, None]` turns `(4,)` into `(4, 1)`.

### Example 4.4: `argsort` to sort one array by another

```python
names = np.array(["Ada", "Ben", "Cy", "Di"])
scores = np.array([72, 95, 88, 95])
order = np.argsort(-scores, kind="stable")           # positions, best first; ties keep input order
print(order, names[order], scores[order])
ranks = np.empty_like(order)
ranks[order] = np.arange(len(order))                  # the inverse permutation: each person's rank
print(ranks)

K = np.array([[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 2.0]])      # three equal masses, springs (CM-6)
w2, V = np.linalg.eigh(K)                             # omega^2 in ascending order
idx = np.argsort(w2)[::-1]                            # fastest mode first
print(np.sqrt(w2[idx]).round(4))
assert np.allclose(K @ V[:, idx], V[:, idx] * w2[idx])
```
```text
[1 3 2 0] ['Ben' 'Di' 'Cy' 'Ada'] [95 95 88 72]
[3 0 2 1]
[1.8478 1.4142 0.7654]
```

`np.argsort(-scores)` returns *which positions* would sort, not the sorted values. Indexing any same-length array with `order` applies that permutation, so one `argsort` reorders names, scores and anything else in step; `ranks[order] = arange` inverts it. In the physics half, `V[:, idx]` reorders eigenvector *columns* by the same permutation, and the assert confirms $K\mathbf v_i=\omega_i^2\mathbf v_i$ still holds. The frequencies $\sqrt{2+\sqrt2},\sqrt2,\sqrt{2-\sqrt2}$ are the chain's three normal modes.

### Example 4.5: `einsum("ij,jk->ik")` and `einsum("ii")`

Let $A=\begin{pmatrix}1&2\\3&4\end{pmatrix}$, $B=\begin{pmatrix}0&1\\1&0\end{pmatrix}$. The string `"ij,jk->ik"` says $C_{ik}=\sum_jA_{ij}B_{jk}$: $C_{00}=1\cdot0+2\cdot1=2$, $C_{01}=1\cdot1+2\cdot0=1$, $C_{10}=3\cdot0+4\cdot1=4$, $C_{11}=3\cdot1+4\cdot0=3$. `"ii"` repeats a letter within one operand and has no output: take the diagonal and sum it, $1+4=5$.

```python
A = np.array([[1, 2], [3, 4]])
B = np.array([[0, 1], [1, 0]])
print(np.einsum("ij,jk->ik", A, B).tolist(), np.einsum("ii", A), np.einsum("ii->i", A), np.einsum("ij->ji", A).tolist())
```
```text
[[2, 1], [4, 3]] 5 [1 4] [[1, 3], [2, 4]]
```

The inertia tensor of CM-7, $I_{ab}=\sum_n m_n(\|\mathbf r_n\|^2\delta_{ab}-r_{na}r_{nb})$, is `einsum("n,na,na->", m, r, r) * eye(3) - einsum("n,na,nb->ab", m, r, r)`; `project.py` checks it against a loop and against $\mathrm{diag}(8,1,9)$ for $m=(1,2)$ at $(1,0,0),(0,2,0)$.

### Example 4.6: Mandelbrot by boolean masks, row by row

$c\in\mathbb C$ is in the Mandelbrot set if $z_0=0$, $z_{k+1}=z_k^2+c$ stays bounded; once $|z_k|>2$ it diverges. Each step depends on the last, so loop over iterations $k$ (§3.10, case i) and handle all pixels together, with a boolean array `alive` for pixels that have not escaped. On four points the whole state fits on screen:

```python
c = np.array([0, 1, 1 + 1j, -2])
z = np.zeros_like(c)
alive = np.ones(c.shape, dtype=bool)
count = np.zeros(c.shape, dtype=int)
print("c     =", c)
for k in range(4):
    z[alive] = z[alive] ** 2 + c[alive]        # read through the mask (a copy), square, write back through it
    alive &= np.abs(z) <= 2                    # an escaped pixel becomes False and stays False
    count += alive                             # True counts as 1: iterations survived so far
    print(f"k={k}  z={np.round(z, 2)}  alive={alive.astype(int)}  count={count}")
```
```text
c     = [ 0.+0.j  1.+0.j  1.+1.j -2.+0.j]
k=0  z=[ 0.+0.j  1.+0.j  1.+1.j -2.+0.j]  alive=[1 1 1 1]  count=[1 1 1 1]
k=1  z=[0.+0.j 2.+0.j 1.+3.j 2.+0.j]  alive=[1 1 0 1]  count=[2 2 1 2]
k=2  z=[0.+0.j 5.+0.j 1.+3.j 2.+0.j]  alive=[1 0 0 1]  count=[3 2 1 3]
k=3  z=[0.+0.j 5.+0.j 1.+3.j 2.+0.j]  alive=[1 0 0 1]  count=[4 2 1 4]
```

By hand: $c=0$ stays at 0. $c=1$ gives $z=1,2,5$: it survives $k=0,1$ ($|2|\le2$) and escapes at $k=2$, count 2. $c=1+i$ gives $z_1=1+i$, $z_2=(1+i)^2+(1+i)=1+3i$, $|z_2|=\sqrt{10}>2$, count 1. $c=-2$ gives $z=-2,2,2,\dots$ and never escapes. The output matches. The right-hand `z[alive]` is a mask *read* (a copy), the left-hand `z[alive] = ...` a mask *write* (modifies `z`); updating only live pixels keeps diverged ones from overflowing. On an image:

```python
def mandelbrot_counts(c: np.ndarray, max_iter: int) -> tuple[np.ndarray, list[int]]:
    """Iterations survived (|z| <= 2) by z -> z^2 + c from z = 0, capped at max_iter; also live pixels per step."""
    z = np.zeros_like(c)
    alive = np.ones(c.shape, dtype=bool)
    count = np.zeros(c.shape, dtype=int)
    n_alive = []
    for _ in range(max_iter):
        z[alive] = z[alive] ** 2 + c[alive]
        alive &= np.abs(z) <= 2
        count += alive
        n_alive.append(int(alive.sum()))
    return count, n_alive


xs, ys = np.linspace(-2.2, 0.8, 600), np.linspace(-1.2, 1.2, 400)
C = xs[None, :] + 1j * ys[:, None]               # (400, 600): rows are Im c, columns are Re c
counts, n_alive = mandelbrot_counts(C, 80)
print(C.shape, "live pixels after k = 1, 3, 10, 80:", [n_alive[k - 1] for k in (1, 3, 10, 80)])
print("counts at c = 0, 1, 1+1j:", mandelbrot_counts(np.array([0, 1, 1 + 1j]), 80)[0])
```
```text
(400, 600) live pixels after k = 1, 3, 10, 80: [213584, 137632, 67222, 51626]
counts at c = 0, 1, 1+1j: [80  2  1]
```

An unmasked version (`z = z*z + c; count += abs(z) <= 2`, warnings silenced) gives identical counts; `project.py` asserts that and times both, and the masked form is not meaningfully faster here. Masks buy warning-free correctness, and a real saving only when most pixels die early.

## 5. Why it matters and where it connects

- **PY-7** plots what you build here (`imshow` of `counts`, `loglog` of `msd`); **PY-8:** SciPy takes and returns arrays (`cdist`, `convolve2d`); **PY-9:** a pandas column is an array plus an index; **PY-10:** profiling finds the loop you could not vectorize and Numba compiles the §3.10 loops, with these timings as baseline.
- **LA-2, LA-5, LA-7, LA-9:** `solve`, `lstsq`/`qr`, `eigh`, `svd` are the library forms of those modules' algorithms; each project checks its own against them. **LA-13:** §3.8's `solve` versus `inv` and §3.9's summation orders are the first two experiments of its stability study. **LA-12:** tensor contractions are `einsum` strings.
- **F3:** DFT sample points are `linspace(..., endpoint=False)`. **F5:** Build 7 is the random-walk MSD. **CM-1, CM-5:** $N$ particles are an `(N, 3)` array and the pairwise-force sum is Build 4's `(N, N, 3)` broadcast (Exercise C3); **CM-6:** Example 4.4's normal modes; **CM-7:** the inertia tensor is one `einsum`.
- **ML-1 onward:** a batch is a leading axis, and most deep-learning bugs are shape bugs ("90% tensor-shape bookkeeping", ROADMAP §7). **ML-8:** Build 8's convolution without the flip.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are predict-the-output and bug-fix drills (strides, views, broadcasting, overflow, axis shapes, `arange`); do all of them on paper first. B1–B7 are the module; B1 (strides), B2 (z-scores by broadcasting), B3 (Gram-matrix distances), B4 (stable softmax), B5 (moving average by `cumsum`) and B7 (extend `project.py`) are essential. C1–C3: reimplement `einsum`, find the vectorization break-even size with measurements and a written conclusion, and an $N$-body acceleration with memory-bounded chunking.

## 7. Coding project

`project.py` is the vectorization lab. In order it builds and asserts: (1) Example 4.1's strides, a ten-row views-versus-copies table with `np.shares_memory`, `int8` overflow, Example 4.2; (2) matmul with loops, `np.dot` and `@` for $n=50,100,200$, timed and plotted; (3) broadcasting by hand on the five pairs, Example 4.3, the axis drill; (4) pairwise distances three ways; (5) the $10^4\times10^3$ random-walk ensemble, MSD slope $1\pm0.05$, plot saved; (6) Mandelbrot by masks, counts at known points, masked equals unmasked, figure saved; (7) convolution four ways and a strides moving average against `np.convolve`; (8) `einsum` and the inertia tensor, `argsort`, summation orders, `.npy` and text round trips.

**Build it yourself:** before opening the file, write the triple-loop matmul, `bshape` and the masked Mandelbrot loop from the descriptions, and compare. Then "Your turn": vectorize a provided loopy nearest-centre function; a covariance matrix by broadcasting against `np.cov`; one Game of Life step with `sliding_window_view`; and find the `axis` bug in a provided function (it passes silently on square input).

## 8. Where you'll struggle

- **Views are the bug.** You slice, modify, and the original changes; the symptom appears three functions away (Example 4.2). When a result changes that you did not touch, suspect a view and test with `np.shares_memory`. When you hand a slice to code that may modify it, decide whether you wanted `.copy()`.
- **`axis=` will be wrong half the time for a month.** Use the rule *the result loses the axis you name*, and check `result.shape` against your prediction every time. Redo Build 5 until the ten answers are instant.
- **Broadcasting errors read as noise.** "operands could not be broadcast together with shapes (3,4) (3,)" is the whole diagnosis: right-align the shapes and find the trailing axes that disagree. The worse case is silent: `(n,) - (n, 1)` broadcasts to `(n, n)` with no error. Cure: `keepdims=True`, `[:, None]`, printed shapes.
- **Small traps.** `np.empty` is not zeros, though it often looks like it on your first test. Integer overflow in arrays is silent. `np.arange` with a float step is off by one unpredictably: use `linspace`. `np.vectorize` fools everyone once. `*` is not `@`, and square arrays will not warn you. `np.std` divides by $N$.
- **The hard part is not the API, it is thinking in arrays.** The fifty functions take a week's drill. Seeing a loop's axes (which index is the walker, which is time, which is the pixel) takes months, and the only cure is volume: whenever you write `for i in range(len(...))` over an array, write the vectorized version beside it, even if you keep the loop.

## 9. Check yourself

1. State equation (3.1). For `float64` `a` of shape `(3, 4)`, what are the shape and strides of `a.T`?
2. Which of `a[2:5]`, `a[[2, 3, 4]]`, `a[a > 0]`, `a.reshape(-1)`, `a.T.reshape(-1)` are views for a C-contiguous `a`? Why must advanced indexing copy?
3. What is `np.array([100], dtype=np.int8) * 2`, and why is it silent? Why does `np.arange(1, 1.3, 0.1)` have 4 elements?
4. State the broadcasting rule. Give the results of `(5,1,3)+(4,1)`, `(3,)+(3,1)` and `(2,3)+(3,2)`.
5. Why is `np.vectorize(f)(x)` not fast? What is the fast form of `x if x > 0 else 0`?
6. What goes wrong with `x += 1` inside a function whose argument is a slice of the caller's array?
7. Why prefer `np.linalg.solve(A, b)` to `np.linalg.inv(A) @ b`?
8. Why does `np.sum` of a float32 array beat a sequential sum? What does `np.isclose(1e-10, 2e-10)` return, and why?
9. Write `einsum` strings for the trace, the batch product `S @ T`, and $\sum_{jk}Q_{ijk}W_{jk}$.
10. When is a Python loop over an array the right choice?

## 10. Pair with

- **NumPy user guide: "The absolute basics for beginners", "Indexing on ndarrays", "Broadcasting", "Copies and views" (numpy.org/doc).** The authoritative statements of §3.1, §3.4, §3.5; the indexing page covers advanced-indexing combinations these notes skip. Read "Copies and views" after Build 2.
- **VanderPlas, *Python Data Science Handbook* 2e, ch. 2 (free online).** Gentler and more example-driven; extra practice on fancy indexing and sorting.
- **Harris et al., "Array programming with NumPy," *Nature* 585 (2020).** The broadcasting figure and the ecosystem; reference only.
- **SciPy Lecture Notes §1.4; Johansson, *Numerical Python* 3e, ch. 2.** A compact tour with more on strides; `linalg` and `einsum` in one chapter.
- **Oliphant, *Guide to NumPy*, ch. 2–3 (reference only).** The array object and strides, by NumPy's author.
