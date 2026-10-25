# PY-6 Solutions

Every code block was run, and the ones with `assert`s passed them; the outputs shown are the real ones (timings come from one machine and will differ on yours). The blocks share one session, top to bottom: start your own file with the imports below. Solutions that need helpers from `project.py` re-define them here.

```python
import math
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.special
from numpy.lib.stride_tricks import as_strided, sliding_window_view
from scipy.spatial.distance import cdist

Shape = tuple[int, ...]
FIG = Path("figures")
FIG.mkdir(exist_ok=True)


def best(f, *args, repeat: int = 3) -> float:
    """Smallest wall time (s) of `repeat` calls of f(*args)."""
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        f(*args)
        times.append(time.perf_counter() - t0)
    return min(times)
```

## A. Warm-up

**A1.** `float32` has itemsize 4, so `a` has strides $(3\cdot4\cdot4,\,4\cdot4,\,4)=(48,16,4)$. `transpose(2, 0, 1)` lists the old axes in the new order, so its shape is `(4, 2, 3)` and its strides are the old strides in that order, `(4, 48, 16)`: not C-contiguous. `a[:, ::2, 1]` drops axis 2 (integer index), doubles the axis-1 stride: shape `(2, 2)`, strides `(48, 32)`. `a[..., ::-1]` negates the last stride: `(48, 16, -4)`. All four share memory with `a`.

```python
a = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
want = {"a": ((2, 3, 4), (48, 16, 4), True), "a.transpose(2, 0, 1)": ((4, 2, 3), (4, 48, 16), False),
        "a[:, ::2, 1]": ((2, 2), (48, 32), False), "a[..., ::-1]": ((2, 3, 4), (48, 16, -4), False)}
for label, v in [("a", a), ("a.transpose(2, 0, 1)", a.transpose(2, 0, 1)), ("a[:, ::2, 1]", a[:, ::2, 1]), ("a[..., ::-1]", a[..., ::-1])]:
    print(f"{label:20s} shape {v.shape!s:10} strides {v.strides!s:13} C-contiguous {v.flags['C_CONTIGUOUS']!s:5} shares {np.shares_memory(a, v)}")
    assert (v.shape, v.strides, v.flags["C_CONTIGUOUS"]) == want[label] and np.shares_memory(a, v)
```
```text
a                    shape (2, 3, 4)  strides (48, 16, 4)   C-contiguous True  shares True
a.transpose(2, 0, 1) shape (4, 2, 3)  strides (4, 48, 16)   C-contiguous False shares True
a[:, ::2, 1]         shape (2, 2)     strides (48, 32)      C-contiguous False shares True
a[..., ::-1]         shape (2, 3, 4)  strides (48, 16, -4)  C-contiguous False shares True
```

**A2.** `y = x` binds a second name to the *same* array (PY-1), so `y[y < 0] = 0` writes into the caller's `data`; the function returns the very array it was given. `row = M[1]` is a basic-indexing view, and `row *= 10` is in place, so row 1 of `M` changes.

```python
def zero_negatives(x: np.ndarray) -> np.ndarray:
    y = x
    y[y < 0] = 0
    return y


data = np.array([3, -1, 4, -1, 5])
out = zero_negatives(data)
print(data, out is data)
M = np.arange(6).reshape(2, 3)
row = M[1]
row *= 10
print(M)
```
```text
[3 0 4 0 5] True
[[ 0  1  2]
 [30 40 50]]
```

Fixes: copy first, or never write into the argument.

```python
def zero_negatives_copy(x: np.ndarray) -> np.ndarray:
    y = x.copy()
    y[y < 0] = 0
    return y


def zero_negatives_pure(x: np.ndarray) -> np.ndarray:
    return np.where(x < 0, 0, x)          # or np.maximum(x, 0): builds a new array


data = np.array([3, -1, 4, -1, 5])
for f in (zero_negatives_copy, zero_negatives_pure):
    assert f(data).tolist() == [3, 0, 4, 0, 5] and data.tolist() == [3, -1, 4, -1, 5] and f(data) is not data
M = np.arange(6).reshape(2, 3)
row = M[1].copy()                          # or: row = M[1] * 10
row *= 10
assert M.tolist() == [[0, 1, 2], [3, 4, 5]]
```

**A3.** Right-align and compare.

```python
cases = [((4, 3), (4,)), ((4, 3), (4, 1)), ((5, 1), (1, 6)), ((2, 3, 4), (3, 1)), ((2, 3, 4), (3,)), ((6,), (6, 1))]
for s, t in cases:
    try:
        print(f"{s!s:10} + {t!s:7} -> {np.broadcast_shapes(s, t)}")
    except ValueError:
        print(f"{s!s:10} + {t!s:7} -> error")
v = np.arange(4.0)
assert (np.ones((4, 3)) + v[:, None]).shape == (4, 3)                  # (4,) -> (4, 1): one value per row
assert (np.ones((2, 3, 4)) + np.arange(3.0)[:, None]).shape == (2, 3, 4)
```
```text
(4, 3)     + (4,)    -> error
(4, 3)     + (4, 1)  -> (4, 3)
(5, 1)     + (1, 6)  -> (5, 6)
(2, 3, 4)  + (3, 1)  -> (2, 3, 4)
(2, 3, 4)  + (3,)    -> error
(6,)       + (6, 1)  -> (6, 6)
```

`(4,3)+(4,)` aligns 4 against 3: error; make `v` a column, `v[:, None]` (shape `(4, 1)`, which is the second line). `(2,3,4)+(3,1)` works: `(3,1)` becomes `(1,3,1)`. `(2,3,4)+(3,)` aligns 3 against 4: error; fix with `[:, None]` as in the last assert. `(6,)+(6,1)` is *not* an error: it is a `(6, 6)` outer sum, the silent trap of §8; if you meant elementwise, make both `(6,)` or both `(6, 1)`.

**A4.**

```python
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    r1 = np.uint8(250) + np.uint8(10)
print("scalar:", r1, "| warning:", caught[0].message)
r2 = np.array([250], dtype=np.uint8) + 10
print("array + 10:", r2)
try:
    np.array([250], dtype=np.uint8) + 300
except OverflowError as err:
    print("array + 300:", type(err).__name__, err)
w8 = np.arange(5, dtype=np.int8) * 100
print(w8, w8.sum(), w8.sum().dtype)
print(np.float32(16777216) + np.float32(1))
assert r1 == 4 and r2.tolist() == [4] and w8.tolist() == [0, 100, -56, 44, -112] and w8.sum() == -24
assert int(np.uint8(250)) + 10 == 260 and (np.array([250], dtype=np.uint8).astype(np.int16) + 10).tolist() == [260]
assert (np.arange(5, dtype=np.int64) * 100).tolist() == [0, 100, 200, 300, 400]
```
```text
scalar: 4 | warning: overflow encountered in scalar add
array + 10: [4]
array + 300: OverflowError Python integer 300 out of bounds for uint8
[   0  100  -56   44 -112] -24 int64
1.6777216e+07
```

The scalar sum wraps with a warning ($260-256=4$); the array sum wraps silently; adding the Python integer 300 raises in NumPy 2 because 300 does not fit `uint8`. The `int8` product wraps ($200\to-56$, $300\to44$, $400\to-112$), and `.sum()` then adds the *wrapped* values in `int64`: $0+100-56+44-112=-24$. A `float32` has a 24-bit significand, so $2^{24}+1=16777217$ rounds back to $16777216$. Correct versions: convert first (`astype(np.int16)`, or `int(...)`, or make the original `int64`).

**A5.**

```python
rng = np.random.default_rng(0)
X = rng.standard_normal((6, 4))
mean_feat = X.mean(axis=0)                           # (4,)    axis 0 (samples) disappears
norm_sample = np.linalg.norm(X, axis=1)              # (6,)
argmax_feat = X.argmax(axis=1)                       # (6,)
X_unit = X / norm_sample[:, None]                    # (6, 4)  the (6,) norms need a column shape
top_var = X.var(axis=0, ddof=1).argmax()             # ()      one feature index
running = X.cumsum(axis=0)                           # (6, 4)  cumsum keeps the axis
print(mean_feat.shape, norm_sample.shape, argmax_feat.shape, X_unit.shape, np.shape(top_var), running.shape)
for i in range(6):
    assert np.isclose(norm_sample[i], math.sqrt(sum(v * v for v in X[i]))) and argmax_feat[i] == max(range(4), key=lambda j: X[i, j])
assert np.allclose(np.linalg.norm(X_unit, axis=1), 1) and np.allclose(running[-1], X.sum(axis=0))
assert all(np.isclose(mean_feat[j], sum(X[:, j]) / 6) for j in range(4))
assert top_var == max(range(4), key=lambda j: sum((X[:, j] - X[:, j].mean()) ** 2))
```
```text
(4,) (6,) (6,) (6, 4) () (6, 4)
```

**A6.** $\lceil(1-0)/0.1\rceil=\lceil10.0\rceil=10$, but $(0.4-0.1)/0.1=3.0000000000000004$ and $\lceil\cdot\rceil=4$, so the "excluded" endpoint 0.4 appears.

```python
print(len(np.arange(0, 1, 0.1)), len(np.arange(0.1, 0.4, 0.1)), np.arange(0.1, 0.4, 0.1)[-1], (0.4 - 0.1) / 0.1)


def grid(a: float, b: float, n: int) -> np.ndarray:
    """n equally spaced points from a to b, both included."""
    return np.linspace(a, b, n)


def periodic_grid(a: float, b: float, n: int) -> np.ndarray:
    """n equally spaced points from a, excluding b (one period of a periodic function)."""
    return np.linspace(a, b, n, endpoint=False)


g, p = grid(0.1, 0.4, 4), periodic_grid(0.0, 2 * np.pi, 8)
assert len(g) == 4 and np.isclose(g[-1], 0.4) and np.allclose(np.diff(g), 0.1)
assert len(p) == 8 and p[-1] < 2 * np.pi and np.allclose(np.diff(p), 2 * np.pi / 8)
```
```text
10 4 0.4 3.0000000000000004
```


## B. Core

**B1.** Strides are cumulative products of the shape, starting from the fast axis (last for C, first for F). Contiguity: walking from the fastest axis, each stride must equal itemsize times the product of the lengths already passed; an axis of length 1 is never stepped along, so its stride is irrelevant.

```python
def strides_for(shape: Shape, itemsize: int, order: str = "C") -> Shape:
    """Strides in bytes: (3.2) for order 'C'; the mirror image for 'F'."""
    n = len(shape)
    strides, step = [0] * n, itemsize
    for k in (range(n - 1, -1, -1) if order == "C" else range(n)):
        strides[k] = step
        step *= shape[k]
    return tuple(strides)


def is_c_contiguous(shape: Shape, strides: Shape, itemsize: int) -> bool:
    expected = itemsize
    for n_k, s_k in zip(reversed(shape), reversed(strides)):
        if n_k == 1:
            continue
        if s_k != expected:
            return False
        expected *= n_k
    return True


base = np.zeros((3, 4))
odd = as_strided(np.zeros(12), shape=(3, 1, 4), strides=(32, 999, 8))        # a bogus stride on a length-1 axis
arrays = {"C (3,4)": base, "F (3,4)": np.zeros((3, 4), order="F"), "transpose": base.T, "a[::2]": base[::2],
          "a[:, 1]": base[:, 1], "a[:1]": base[:1], "length-1 axis": odd, "3-D C": np.zeros((2, 3, 4), dtype=np.float32)}
for label, z in arrays.items():
    mine = is_c_contiguous(z.shape, z.strides, z.itemsize)
    print(f"{label:14s} strides {z.strides!s:14} mine {mine!s:5} numpy {z.flags['C_CONTIGUOUS']}")
    assert mine == z.flags["C_CONTIGUOUS"]
for shape in [(3, 4), (2, 3, 4), (5,), (2, 1, 3)]:
    assert strides_for(shape, 8, "C") == np.zeros(shape).strides
    assert strides_for(shape, 8, "F") == np.zeros(shape, order="F").strides
```
```text
C (3,4)        strides (32, 8)        mine True  numpy True
F (3,4)        strides (8, 24)        mine False numpy False
transpose      strides (8, 32)        mine False numpy False
a[::2]         strides (64, 8)        mine False numpy False
a[:, 1]        strides (32,)          mine False numpy False
a[:1]          strides (32, 8)        mine True  numpy True
length-1 axis  strides (32, 999, 8)   mine True  numpy True
3-D C          strides (48, 16, 4)    mine True  numpy True
```

NumPy's own flag ignores length-1 axes exactly as the function does, which is why `odd` counts as contiguous.

**B2.** Subtract the mean and divide by the standard deviation, both with `keepdims=True` so they broadcast along `axis`. The trap is the constant column: its computed standard deviation is not exactly 0 because the mean of copies of 0.1 is not exactly 0.1, so dividing by that tiny nonzero std produces $\pm1$ garbage instead of 0 (and `X - mean` is a tiny nonzero too). Test against a tolerance, not equality, and set those entries to 0 explicitly.

```python
def zscore(X: np.ndarray, ddof: int = 1, axis: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Z = (X - mean) / std along `axis`; constant slices (std ~ 0 relative to scale) give zeros. Returns (Z, mean, std)."""
    mean = X.mean(axis=axis, keepdims=True)
    std = X.std(axis=axis, ddof=ddof, keepdims=True)
    constant = std <= 1e-12 * np.maximum(1.0, np.abs(mean))
    Z = np.where(constant, 0.0, (X - mean) / np.where(constant, 1.0, std))
    return Z, mean, np.where(constant, 0.0, std)


rng = np.random.default_rng(1)
X = rng.normal(5.0, 3.0, size=(50, 4))
X[:, 2] = 0.1                                                    # a constant column
Z, mean, std = zscore(X)
for j in range(4):                                               # loop reference
    col = X[:, j]
    s = math.sqrt(sum((v - col.mean()) ** 2 for v in col) / 49)
    ref = np.zeros(50) if j == 2 else (col - col.mean()) / s
    assert np.allclose(Z[:, j], ref)
print("naive std of the constant column:", X[:, 2].std(ddof=1))
assert Z.shape == (50, 4) and mean.shape == (1, 4) and (Z[:, 2] == 0).all()
Zr, mr, sr = zscore(X[:, [0, 1, 3]], axis=1)                     # per-row standardisation
assert Zr.shape == (50, 3) and mr.shape == (50, 1) and np.allclose(Zr.mean(axis=1), 0) and np.allclose(Zr.std(axis=1, ddof=1), 1)
assert np.allclose(Z[:, [0, 1, 3]].mean(axis=0), 0) and np.allclose(Z[:, [0, 1, 3]].std(axis=0, ddof=1), 1)
```
```text
naive std of the constant column: 2.803736533363884e-17
```

For `axis=1` the means have shape `(50, 1)`: `keepdims=True` is what lets `X - mean` align without the `[:, None]` bookkeeping.

**B3.** (a) and (c) agree with `cdist`; (b) fails because with coordinates near $10^8$ each $\|\mathbf x_i\|^2\approx3\times10^{16}$, where the spacing of doubles is about 4, while the true $D_{ij}^2$ is of order 1 to 10: the difference of three huge numbers keeps no correct digits. Centering first shrinks the norms to order 3.

```python
def dist_gram(X: np.ndarray, center: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """D_ij from D2 = |x_i|^2 + |x_j|^2 - 2 x_i.x_j with one matmul. Returns (D, D2 before clipping)."""
    if center:
        X = X - X.mean(axis=0)
    sq = (X * X).sum(axis=1)
    D2 = sq[:, None] + sq[None, :] - 2 * (X @ X.T)
    return np.sqrt(np.maximum(D2, 0.0)), D2


def dist_broadcast(X: np.ndarray) -> np.ndarray:
    diff = X[:, None, :] - X[None, :, :]
    return np.sqrt((diff ** 2).sum(axis=-1))


rng = np.random.default_rng(2)
X = rng.standard_normal((300, 3))
D, _ = dist_gram(X)
assert np.allclose(D, cdist(X, X), atol=1e-7)                           # (a)
Xs = X + 1e8
D_bad, D2_bad = dist_gram(Xs)
print(f"(b) offset 1e8: max |error| = {np.abs(D_bad - cdist(Xs, Xs)).max():.2f}, min D2 = {D2_bad.min():.1f}")
D_fix, _ = dist_gram(Xs, center=True)
print(f"(c) centred:    max |error| = {np.abs(D_fix - cdist(Xs, Xs)).max():.1e}")
assert np.abs(D_bad - cdist(Xs, Xs)).max() > 0.5 and D2_bad.min() < 0 and np.allclose(D_fix, cdist(Xs, Xs), atol=1e-6)

X = np.random.default_rng(3).standard_normal((2000, 3))
t_gram, t_bc, t_cd = best(dist_gram, X), best(dist_broadcast, X), best(lambda Z: cdist(Z, Z), X)
assert np.allclose(dist_gram(X)[0], cdist(X, X), atol=1e-6) and np.allclose(dist_broadcast(X), cdist(X, X))
print(f"(d) n=2000: gram {t_gram * 1e3:.0f} ms, broadcast {t_bc * 1e3:.0f} ms, cdist {t_cd * 1e3:.0f} ms")
print(f"    temporaries: gram D2 {2000 * 2000 * 8 / 1e6:.0f} MB; broadcast diff {2000 * 2000 * 3 * 8 / 1e6:.0f} MB (plus diff**2: another {2000 * 2000 * 3 * 8 / 1e6:.0f} MB); cdist none")
```
```text
(b) offset 1e8: max |error| = 3.96, min D2 = -16.0
(c) centred:    max |error| = 6.5e-15
(d) n=2000: gram 50 ms, broadcast 156 ms, cdist 9 ms
    temporaries: gram D2 32 MB; broadcast diff 96 MB (plus diff**2: another 96 MB); cdist none
```

Conclusion: the Gram form does its work in one BLAS call and stores only an $n\times n$ array, which pays off when $d$ is large (hundreds of features; at $d=3$ it was slower than `cdist` here), but it needs centering and clipping and loses digits for nearby points; the broadcast form is exact and simplest but stores $O(n^2d)$; `cdist` is exact, fast, and stores only the result, so it is the default unless you need something it lacks.

**B4.** The naive form overflows: `np.exp(1000)` is `inf` and `inf / inf` is `nan`. Subtracting the row maximum changes nothing mathematically (the factor $e^{-m}$ cancels) and makes the largest exponent 0.

```python
def softmax(Z: np.ndarray, axis: int = -1) -> np.ndarray:
    """softmax(z)_i = exp(z_i - m) / sum_j exp(z_j - m), m = max z along axis."""
    E = np.exp(Z - Z.max(axis=axis, keepdims=True))
    return E / E.sum(axis=axis, keepdims=True)


def log_softmax(Z: np.ndarray, axis: int = -1) -> np.ndarray:
    """z_i - log sum_j exp(z_j), computed as z_i - m - log sum_j exp(z_j - m)."""
    m = Z.max(axis=axis, keepdims=True)
    return Z - m - np.log(np.exp(Z - m).sum(axis=axis, keepdims=True))


z = np.array([1000.0, 1001.0, 1002.0])
with np.errstate(all="ignore"):
    print("naive:", np.exp(z) / np.exp(z).sum())
print("stable:", softmax(z))
rng = np.random.default_rng(4)
Z = rng.standard_normal((2, 5, 7)) * 5
S = softmax(Z)
assert S.shape == (2, 5, 7) and np.allclose(S.sum(axis=-1), 1) and (S >= 0).all()
assert np.allclose(softmax(Z + 1000.0), S) and np.allclose(softmax(Z + rng.standard_normal((2, 5, 1)) * 100), S)
assert np.allclose(S, scipy.special.softmax(Z, axis=-1)) and np.allclose(log_softmax(Z), scipy.special.log_softmax(Z, axis=-1))
assert np.allclose(log_softmax(z), np.log(softmax(z))) and np.allclose(softmax(Z, axis=0).sum(axis=0), 1)
```
```text
naive: [nan nan nan]
stable: [0.09003057 0.24472847 0.66524096]
```

**B5.**

```python
def ma_windows(x: np.ndarray, w: int) -> np.ndarray:
    return sliding_window_view(x, w).mean(axis=-1)


def ma_convolve(x: np.ndarray, w: int) -> np.ndarray:
    return np.convolve(x, np.ones(w) / w, mode="valid")


def ma_cumsum(x: np.ndarray, w: int) -> np.ndarray:
    """(S[i+w] - S[i]) / w with S the running sum, S[0] = 0, kept in x's dtype."""
    S = np.concatenate([np.zeros(1, dtype=x.dtype), np.cumsum(x)])
    return (S[w:] - S[:-w]) / w


rng = np.random.default_rng(5)
x = rng.standard_normal(10_000)
ref = ma_windows(x, 50)
assert np.allclose(ma_convolve(x, 50), ref) and np.allclose(ma_cumsum(x, 50), ref)

big = (1e6 + rng.standard_normal(10 ** 6)).astype(np.float32)           # data already stored in float32
truth = ma_windows(big.astype(np.float64), 1000)
for f in (ma_windows, ma_convolve, ma_cumsum):
    print(f"{f.__name__:12s} max |error| vs float64: {np.abs(f(big, 1000) - truth).max():10.3g}   {best(f, big, 1000, repeat=1) * 1e3:8.1f} ms")
```
```text
ma_windows   max |error| vs float64:      0.107      245.6 ms
ma_convolve  max |error| vs float64:   4.66e-10      136.0 ms
ma_cumsum    max |error| vs float64:    1.7e+04        4.2 ms
```

`ma_cumsum` loses almost everything: its running sum reaches $10^{12}$ in `float32`, where adjacent representable numbers are $2^{16}=65536$ apart, and each average is a difference of two such numbers. `ma_windows` is off by about 0.1, under two spacings of `float32` near $10^6$ (0.0625), which is as good as the stored data allows. `ma_convolve` is accurate to $10^{-10}$ because `np.ones(w) / w` is `float64` and promotes the arithmetic. The windows and convolution methods do about $w$ multiply-adds per output, the cumulative sum about two operations, which is why `ma_cumsum` is the fast one (3.9 ms against 233 ms here) and why it is only safe in `float64`, or with centred data.

**B6.**

```python
x = np.linspace(0, 1, 40)
rng = np.random.default_rng(1)
y = np.sin(2 * np.pi * x) + 0.05 * rng.standard_normal(40)
print(f"{'deg':>4} {'cond(V)':>10} {'lstsq vs polyfit':>17} {'normal eq. vs lstsq':>20}")
for deg in (3, 8, 15):
    V = x[:, None] ** np.arange(deg + 1)                              # V[i, j] = x_i ** j
    c_ls = np.linalg.lstsq(V, y, rcond=None)[0]
    c_pf = np.polyfit(x, y, deg)[::-1]                                # polyfit lists the highest power first
    try:
        c_ne = np.linalg.solve(V.T @ V, V.T @ y)
        err_ne = np.linalg.norm(c_ne - c_ls) / np.linalg.norm(c_ls)
    except np.linalg.LinAlgError:
        err_ne = float("inf")
    err_pf = np.linalg.norm(c_ls - c_pf) / np.linalg.norm(c_ls)
    print(f"{deg:4d} {np.linalg.cond(V):10.1e} {err_pf:17.1e} {err_ne:20.1e}")
    assert V.shape == (40, deg + 1) and err_pf < 1e-6 * max(1.0, np.linalg.cond(V))
```
```text
 deg    cond(V)  lstsq vs polyfit  normal eq. vs lstsq
   3    1.2e+02           1.7e-15              9.2e-14
   8    6.2e+05           4.2e-12              9.8e-07
  15    1.5e+11           7.0e-07              1.0e+00
```

Squaring the matrix squares the condition number ($\kappa(V^\top V)=\kappa(V)^2$, LA-9), so the normal equations lose twice as many digits as `lstsq` (which uses an orthogonal factorization, LA-5) and break down once $\kappa(V)^2$ approaches $1/\epsilon_{\text{mach}}\approx4.5\times10^{15}$: the degree-8 error is already visible and degree 15 is useless. Monomials are also a bad basis on $[0,1]$ (the columns become nearly parallel), which is why $\kappa(V)$ is large to begin with; orthogonal polynomials (`numpy.polynomial.Chebyshev.fit`) avoid it.

**B7.** (a) The four "Your turn" functions, checked here with the same data as the project. (b) The 2-D walk.

```python
def loopy_assign(X: np.ndarray, C: np.ndarray) -> np.ndarray:
    out = []
    for x in X:
        best_k, best_d = 0, float("inf")
        for k, c in enumerate(C):
            d = sum((x[j] - c[j]) ** 2 for j in range(len(x)))
            if d < best_d:
                best_k, best_d = k, d
        out.append(best_k)
    return np.array(out)


def assign_vec(X: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Index of the nearest row of C for each row of X: argmin over the k axis of an (n, k, d) difference."""
    return (((X[:, None, :] - C[None, :, :]) ** 2).sum(axis=-1)).argmin(axis=1)


def cov_broadcast(X: np.ndarray) -> np.ndarray:
    """Sample covariance (1/(n-1)) Xc^T Xc of the centred data."""
    Xc = X - X.mean(axis=0)
    return Xc.T @ Xc / (len(X) - 1)


def life_step(grid: np.ndarray) -> np.ndarray:
    """One Game of Life step, dead borders: count the 8 neighbours with a 3x3 window sum on a padded copy."""
    nbrs = sliding_window_view(np.pad(grid, 1), (3, 3)).sum(axis=(-2, -1)) - grid
    return ((nbrs == 3) | ((grid == 1) & (nbrs == 2))).astype(grid.dtype)


def row_normalize(X: np.ndarray) -> np.ndarray:
    return X / X.sum(axis=1, keepdims=True)


Xa, Ca = np.random.default_rng(20).standard_normal((60, 3)), np.random.default_rng(21).standard_normal((5, 3))
assert np.array_equal(assign_vec(Xa, Ca), loopy_assign(Xa, Ca))
Xcov = np.random.default_rng(22).standard_normal((50, 4))
assert np.allclose(cov_broadcast(Xcov), np.cov(Xcov, rowvar=False))
blinker = np.zeros((5, 5), dtype=int)
blinker[2, 1:4] = 1
once = life_step(blinker)
assert once[1:4, 2].tolist() == [1, 1, 1] and once.sum() == 3 and np.array_equal(life_step(once), blinker)
glider = np.zeros((8, 8), dtype=int)
glider[0, 1] = glider[1, 2] = glider[2, 0] = glider[2, 1] = glider[2, 2] = 1
g = glider
for _ in range(4):
    g = life_step(g)
assert np.array_equal(g[1:, 1:], glider[:-1, :-1]) and g.sum() == 5          # four steps move a glider by (1, 1)
Xsq = np.arange(1.0, 10.0).reshape(3, 3)
assert all(np.allclose(row_normalize(X).sum(axis=1), 1) for X in (Xsq, np.ones((3, 4))))
```

The bug in `row_normalize_buggy`: `X.sum(axis=1)` has shape `(n,)`, and `X / (n,)` aligns that vector with the *last* axis, the columns; so column $j$ is divided by the sum of row $j$. For a square matrix the shapes agree and nothing is raised; for `(3, 4)` the `3` meets the `4` and NumPy refuses. `keepdims=True` gives shape `(n, 1)`, which stretches along the columns as intended.

```python
def walks_2d(n_walkers: int, n_steps: int, seed: int) -> np.ndarray:
    """Positions of 2-D lattice walkers, shape (walkers, time, 2): table lookup, then cumsum along time."""
    table = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]], dtype=np.int16)
    choices = np.random.default_rng(seed).integers(0, 4, size=(n_walkers, n_steps))
    return table[choices].cumsum(axis=1, dtype=np.int16)               # (W, T) integers index a (4, 2) table -> (W, T, 2)


pos = walks_2d(5000, 500, seed=7)
t = np.arange(1, 501)
msd = (pos.astype(np.float64) ** 2).sum(axis=2).mean(axis=0)           # |r|^2 over the coordinate axis, then average over walkers
slope = np.polyfit(np.log(t), np.log(msd), 1)[0]
se = math.sqrt(1 - 1 / 500) / math.sqrt(5000)                          # relative standard error of the MSD at t=500
print(pos.shape, f"MSD(500) = {msd[-1]:.1f}, slope = {slope:.4f}, 5 standard errors = {5 * se:.3f}")
assert pos.shape == (5000, 500, 2) and abs(slope - 1) < 0.05 and abs(msd[-1] / 500 - 1) < 5 * se
fig, ax = plt.subplots(figsize=(5, 5))
for p in pos[:10]:
    ax.plot(p[:, 0], p[:, 1], lw=0.7)
ax.set(aspect="equal", xlabel="x", ylabel="y", title="10 of the 5000 two-dimensional walkers")
fig.savefig(FIG / "walk2d.png", dpi=110)
plt.close(fig)
```
```text
(5000, 500, 2) MSD(500) = 505.1, slope = 1.0057, 5 standard errors = 0.071
```

Each step has $\|\text{step}\|^2=1$ and the steps are independent with mean 0, so as in (3.7) $\langle\|\mathbf r_t\|^2\rangle=t$. For this lattice walk $\mathrm{Var}(\|\mathbf r_t\|^2)=t^2-t$, so the relative standard error of an average over $W$ walkers is $\sqrt{1-1/t}/\sqrt W$, which is the `se` used in the assert.


## C. Challenge

**C1.** Version 1 puts every letter on its own axis: an operand lacking a letter gets a length-1 axis there, the product broadcasts to the full index space, and the letters absent from the output are summed. That is the definition of Einstein summation, and its cost is the size of the full index space.

```python
def my_einsum(spec: str, A: np.ndarray, B: np.ndarray, method: str = "broadcast") -> np.ndarray:
    """Two-operand einsum with an explicit output and no repeated letter inside an operand."""
    ins, out = spec.replace(" ", "").split("->")
    sa, sb = ins.split(",")
    letters = sorted(set(sa) | set(sb))                        # the full index space, alphabetical

    def expand(T: np.ndarray, s: str) -> np.ndarray:                                          # T's axes -> alphabetical, with a 1 for each missing letter
        T = np.transpose(T, np.argsort([letters.index(c) for c in s]))
        return T.reshape([T.shape[sorted(s, key=letters.index).index(c)] if c in s else 1 for c in letters])

    if method == "broadcast":
        full = expand(A, sa) * expand(B, sb)                   # shape over all letters
        res = full.sum(axis=tuple(i for i, c in enumerate(letters) if c not in out))
        kept = [c for c in letters if c in out]
        return np.transpose(res, [kept.index(c) for c in out])
    # method == "tensordot": sum away letters private to one operand, then contract the shared letters
    priv_a = [c for c in sa if c not in sb and c not in out]
    priv_b = [c for c in sb if c not in sa and c not in out]
    A2, sa2 = A.sum(axis=tuple(sa.index(c) for c in priv_a)) if priv_a else A, "".join(c for c in sa if c not in priv_a)
    B2, sb2 = B.sum(axis=tuple(sb.index(c) for c in priv_b)) if priv_b else B, "".join(c for c in sb if c not in priv_b)
    shared = [c for c in sa2 if c in sb2]
    if any(c in out for c in shared):
        raise NotImplementedError("batch letters need the broadcast method")
    res = np.tensordot(A2, B2, axes=([sa2.index(c) for c in shared], [sb2.index(c) for c in shared]))
    order = "".join(c for c in sa2 if c not in shared) + "".join(c for c in sb2 if c not in shared)
    return np.transpose(res, [order.index(c) for c in out])


rng = np.random.default_rng(7)
A, B = rng.standard_normal((3, 4)), rng.standard_normal((4, 2))
S, T = rng.standard_normal((5, 3, 4)), rng.standard_normal((5, 4, 2))
u, w, Q, W = rng.standard_normal(3), rng.standard_normal(4), rng.standard_normal((3, 4, 5)), rng.standard_normal((4, 5))
tests = [("ij,jk->ik", A, B), ("bij,bjk->bik", S, T), ("i,j->ij", u, w), ("ijk,jk->i", Q, W), ("ij,jk->ki", A, B), ("ij,kl->ik", A, rng.standard_normal((2, 6)))]
for spec, X, Y in tests:
    want = np.einsum(spec, X, Y)
    assert np.allclose(my_einsum(spec, X, Y), want)
    try:
        assert np.allclose(my_einsum(spec, X, Y, method="tensordot"), want)
        kind = "both"
    except NotImplementedError:
        kind = "broadcast only"
    print(f"{spec:14s} ok ({kind})")
```
```text
ij,jk->ik      ok (both)
bij,bjk->bik   ok (broadcast only)
i,j->ij        ok (both)
ijk,jk->i      ok (both)
ij,jk->ki      ok (both)
ij,kl->ik      ok (both)
```

Cost, for `ij,jk->ik` at $n=300$:

```python
n = 300
P, R = rng.standard_normal((n, n)), rng.standard_normal((n, n))
print(f"Version 1 intermediate: {n ** 3 * 8 / 1e6:.0f} MB")
rows = [("my_einsum, broadcast", lambda: my_einsum("ij,jk->ik", P, R)), ("my_einsum, tensordot", lambda: my_einsum("ij,jk->ik", P, R, method="tensordot")),
        ("np.einsum", lambda: np.einsum("ij,jk->ik", P, R)), ("np.einsum, optimize=True", lambda: np.einsum("ij,jk->ik", P, R, optimize=True)), ("P @ R", lambda: P @ R)]
for label, f in rows:
    print(f"{label:26s} {best(f, repeat=2) * 1e3:9.2f} ms")
assert np.allclose(my_einsum("ij,jk->ik", P, R), P @ R)
```
```text
Version 1 intermediate: 216 MB
my_einsum, broadcast           97.07 ms
my_einsum, tensordot            7.98 ms
np.einsum                      11.83 ms
np.einsum, optimize=True        4.90 ms
P @ R                           7.94 ms
```

`@` wins because it reaches BLAS, which blocks for cache and uses vector instructions; `tensordot` reduces to the same `@` after reshaping, so it matches it; plain `np.einsum` loops in C but without BLAS; the broadcast version also writes an $n^3$-element intermediate (216 MB) and then sums it, so it is the slowest and the only one that can run out of memory. `einsum` is still the right tool when no named function exists (multi-index contractions, batched products with permuted axes, traces and diagonals), and `optimize=True` lets it call `tensordot` for you.

**C2.** The loop's cost per element is constant; NumPy has a fixed cost per call (three calls here: `sin`, `**2`, `sum`) plus a small cost per element.

```python
def loop_sum(xs: list[float]) -> float:
    s = 0.0
    for v in xs:
        s += math.sin(v) ** 2
    return s


def np_sum(x: np.ndarray) -> float:
    return (np.sin(x) ** 2).sum()


def np_dot(x: np.ndarray) -> float:
    s = np.sin(x)
    return np.dot(s, s)


def timed(f, arg, budget: float = 0.02) -> float:
    """Seconds per call, averaging enough calls to fill `budget` seconds; best of 3 rounds."""
    t0 = time.perf_counter()
    f(arg)
    one = max(time.perf_counter() - t0, 1e-8)
    reps = max(1, int(budget / one))
    best_t = float("inf")
    for _ in range(3):
        t0 = time.perf_counter()
        for _ in range(reps):
            f(arg)
        best_t = min(best_t, (time.perf_counter() - t0) / reps)
    return best_t


Ns = np.unique(np.logspace(0, 6, 25).astype(int))
rng = np.random.default_rng(0)
t_loop, t_sum, t_dot = [], [], []
for N in Ns:
    x = rng.random(N)
    xl = x.tolist()
    assert np.isclose(loop_sum(xl), np_sum(x)) and np.isclose(np_dot(x), np_sum(x))
    t_loop.append(timed(loop_sum, xl))
    t_sum.append(timed(np_sum, x))
    t_dot.append(timed(np_dot, x))
t_loop, t_sum, t_dot = map(np.array, (t_loop, t_sum, t_dot))
cross = Ns[np.argmax(t_sum < t_loop)]
print(f"NumPy fixed cost (N=1): {t_sum[0] * 1e6:.1f} us per call; loop at N=1: {t_loop[0] * 1e6:.2f} us")
print(f"asymptotic per element: loop {t_loop[-1] / Ns[-1] * 1e9:.0f} ns, numpy {t_sum[-1] / Ns[-1] * 1e9:.1f} ns")
print(f"crossover: NumPy is faster from N = {cross}")
fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(Ns, t_loop / Ns, "o-", ms=3, label="Python loop, math.sin")
ax.loglog(Ns, t_sum / Ns, "s-", ms=3, label="np.sin(x)**2 then .sum()")
ax.loglog(Ns, t_dot / Ns, "^-", ms=3, label="np.dot(sin x, sin x)")
ax.set(xlabel="N", ylabel="seconds per element", title="Per-call overhead dominates below the crossover")
ax.grid(alpha=0.3, which="both")
ax.legend()
fig.tight_layout()
fig.savefig(FIG / "break_even.png", dpi=110)
plt.close(fig)
assert cross < 200 and t_sum[0] > t_loop[0] and t_loop[-1] > 3 * t_sum[-1]
```
```text
NumPy fixed cost (N=1): 3.9 us per call; loop at N=1: 0.15 us
asymptotic per element: loop 85 ns, numpy 11.0 ns
crossover: NumPy is faster from N = 56
```

For the 3-vector question, measure the cost of one small operation both ways:

```python
a3, b3 = (1.0, 2.0, 3.0), (4.0, 5.0, 6.0)
A3, B3 = np.array(a3), np.array(b3)
t_py = best(lambda: a3[0] * b3[0] + a3[1] * b3[1] + a3[2] * b3[2], repeat=5)
t_npd = best(lambda: np.dot(A3, B3), repeat=5)
print(f"3-vector dot product: tuple arithmetic {t_py * 1e9:.0f} ns, np.dot {t_npd * 1e9:.0f} ns")
assert t_npd > 2 * t_py
```
```text
3-vector dot product: tuple arithmetic 191 ns, np.dot 892 ns
```

Conclusion: NumPy costs about 1.6 µs per call (three calls here) and the loop about 0.16 µs per element, so the crossover is N = 31 on this machine; beyond it NumPy wins by up to the ratio of the asymptotic per-element costs (83 ns against 11 ns, a factor of about 7, because `sin` dominates both). A loop that does a 3-vector operation per step for $10^6$ steps should keep each vector as plain floats (or a tuple) and use `math`: the NumPy dot product of two 3-vectors cost 966 ns against 190 ns for the written-out arithmetic. The way to make such a simulation fast is not to call NumPy on 3 numbers a million times but to put the long axis (particles, time, walkers) inside the array and loop only over what truly cannot be vectorized (steps of an integrator, §3.10 case i): CM-1's trajectory is a `(steps, 3)` array filled by a short loop or a single `solve_ivp` call (PY-8).

**C3.**

```python
def accel_loops(r: np.ndarray, m: np.ndarray, eps: float, G: float = 1.0) -> np.ndarray:
    N = len(r)
    a = np.zeros_like(r)
    for i in range(N):
        for j in range(N):
            if i != j:
                d = r[j] - r[i]
                a[i] += G * m[j] * d / (d @ d + eps ** 2) ** 1.5
    return a


def accel_full(r: np.ndarray, m: np.ndarray, eps: float, G: float = 1.0) -> np.ndarray:
    """a_i = G sum_j m_j d_ij / (|d_ij|^2 + eps^2)^(3/2), d_ij = r_j - r_i, via an (N, N, 3) array."""
    d = r[None, :, :] - r[:, None, :]
    w = (np.einsum("ijk,ijk->ij", d, d) + eps ** 2) ** -1.5
    np.fill_diagonal(w, 0.0)
    return G * np.einsum("ij,j,ijk->ik", w, m, d)


def accel_chunked(r: np.ndarray, m: np.ndarray, eps: float, G: float = 1.0, block: int = 256) -> np.ndarray:
    """Same, `block` rows of i at a time: the temporary is (block, N, 3), not (N, N, 3)."""
    a = np.empty_like(r)
    for lo in range(0, len(r), block):
        hi = min(lo + block, len(r))
        d = r[None, :, :] - r[lo:hi, None, :]
        w = (np.einsum("ijk,ijk->ij", d, d) + eps ** 2) ** -1.5
        w[np.arange(hi - lo), np.arange(lo, hi)] = 0.0
        a[lo:hi] = G * np.einsum("ij,j,ijk->ik", w, m, d)
    return a


rng = np.random.default_rng(8)
eps = 0.05
print(f"{'N':>5} {'loops (s)':>10} {'full (s)':>9} {'chunked (s)':>12} {'d array (MB)':>13} {'chunk d (MB)':>13} {'|sum m a|/sum |m a|':>20}")
for N in (200, 1000, 3000):
    r, m = rng.standard_normal((N, 3)), rng.random(N) + 0.5
    block = max(1, int(50e6 // (24 * N)))
    a_full, a_chunk = accel_full(r, m, eps), accel_chunked(r, m, eps, block=block)
    assert np.allclose(a_full, a_chunk, rtol=1e-12, atol=0)
    t_loop = None
    if N <= 1000:
        t0 = time.perf_counter()
        a_loop = accel_loops(r, m, eps)
        t_loop = time.perf_counter() - t0
        assert np.allclose(a_loop, a_full, rtol=1e-10, atol=1e-12)
    force_scale = np.abs(m[:, None] * a_full).sum()
    third = np.abs((m[:, None] * a_full).sum(axis=0)).max() / force_scale
    assert third < 1e-12
    t_loop_s = f"{t_loop:10.3f}" if t_loop is not None else f"{'skipped':>10}"
    print(f"{N:5d} {t_loop_s} {best(accel_full, r, m, eps, repeat=2):9.4f} {best(accel_chunked, r, m, eps, 1.0, block, repeat=2):12.4f} "
          f"{N * N * 3 * 8 / 1e6:13.1f} {min(block, N) * N * 3 * 8 / 1e6:13.1f} {third:20.1e}")
print(f"full (N, N, 3) array reaches 1 GB at N = {math.sqrt(1e9 / 24):.0f}")
```
```text
    N  loops (s)  full (s)  chunked (s)  d array (MB)  chunk d (MB)  |sum m a|/sum |m a|
  200      0.203    0.0014       0.0016           1.0           1.0              2.2e-17
 1000      5.010    0.0349       0.0342          24.0          24.0              4.0e-17
 3000    skipped    0.3454       0.3302         216.0          50.0              1.5e-17
full (N, N, 3) array reaches 1 GB at N = 6455
```

The loop result, the full broadcast and the chunked version agree to $10^{-12}$ (relative), and $\sum_i m_i\mathbf a_i$ vanishes to rounding relative to the size of the forces: the pairwise terms $m_im_j\mathbf d_{ij}w_{ij}$ are antisymmetric in $(i,j)$ (Newton's third law), a physics check on the code that needs no reference implementation. The full broadcast needs $24N^2$ bytes for `d` alone (plus $8N^2$ for `w`), which passes 1 GB near $N\approx6500$; the chunked version keeps the temporary under 50 MB for any $N$ and costs about the same time, because the work per pair is unchanged and each chunk is still a large array operation (§3.10, case iii).
