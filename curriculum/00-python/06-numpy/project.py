# %% [markdown]
# # PY-6 Project: the vectorization lab
#
# Builds, from raw bytes up, the machinery of notes §3 and checks every number quoted in §1-§4:
# strides and addressing (Build 1), views versus copies (Build 2), broadcasting by hand (Build 3),
# pairwise distances three ways (Build 4), the axis drill (Build 5), matrix multiply (Build 6),
# the random-walk ensemble with its mean-square displacement (Build 7) and 2D convolution (Build 8),
# plus the Mandelbrot set by boolean masks (Example 4.6). Timings are measured at run time and the
# asserts use wide margins (several times below what the notes quote). Figures go to `figures/`.
# Run headless: `MPLBACKEND=Agg python3 project.py` (about 15 s).

# %%
import math
import random
import sys
import tempfile
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.stride_tricks import as_strided, sliding_window_view
from scipy.signal import convolve2d
from scipy.spatial.distance import cdist

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.perf_counter()
NUMPY2 = int(np.__version__.split(".")[0]) >= 2


def best(f, *args, repeat=3):
    """Smallest wall time (s) of `repeat` calls of f(*args)."""
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        f(*args)
        times.append(time.perf_counter() - t0)
    return min(times)


# %% [markdown]
# ## Part 1: memory, strides, dtypes (notes §1, §3.1-§3.3, Example 4.1)
# Build 1: write equations (3.1) and (3.2) yourself and read raw bytes with them.

# %%
def c_strides(shape, itemsize):
    """Strides in bytes of a C-contiguous array: s_{d-1} = itemsize, s_k = s_{k+1} n_{k+1} (3.2)."""
    strides, step = [], itemsize
    for n_k in reversed(shape):
        strides.append(step)
        step *= n_k
    return tuple(reversed(strides))


def read_element(raw, dtype, strides, idx, start=0):
    """Element idx read from the bytes `raw`: byte offset = start + sum_k i_k s_k (3.1)."""
    offset = start + sum(i * s for i, s in zip(idx, strides))
    return np.frombuffer(raw, dtype=dtype, count=1, offset=offset)[0]


n = 1_000_000
xs = [float(i) for i in range(n)]
arr = np.arange(n, dtype=np.float64)
assert sys.getsizeof(xs[1]) == 24 and arr.itemsize == 8 and arr.nbytes == 8_000_000


def sumsq_loop(v):
    s = 0.0
    for t in v:
        s += t * t
    return s


assert np.isclose(sumsq_loop(xs), np.dot(arr, arr))
t_loop, t_np = best(sumsq_loop, xs), best(lambda v: np.dot(v, v), arr)
assert t_loop > 20 * t_np, (t_loop, t_np)
print(f"sum of squares, 10^6 floats: loop {t_loop * 1e3:.1f} ms, np.dot {t_np * 1e3:.3f} ms ({t_loop / t_np:.0f}x)")

a = np.arange(6).reshape(2, 3)
assert (a.shape, a.ndim, a.size, a.dtype, a.itemsize, a.nbytes, a.strides) == ((2, 3), 2, 6, np.dtype("int64"), 8, 48, (24, 8))
assert a.flags["C_CONTIGUOUS"] and not a.flags["OWNDATA"] and a.base is not None
assert (a.T.strides, a[:, 1].strides, a[:, ::2].strides, a[::-1].strides) == ((8, 24), (24,), (24, 16), (-24, 8))
raw = a.base.tobytes()
assert c_strides(a.shape, a.itemsize) == a.strides
for view in (a, a.T, a[::2, 1:], a[::-1]):
    start = view.ctypes.data - a.ctypes.data
    for idx in np.ndindex(view.shape):
        assert read_element(raw, a.dtype, view.strides, idx, start) == view[idx]
assert c_strides((2, 3, 4), 8) == (96, 32, 8) == np.zeros((2, 3, 4)).strides

# Example 4.1: the 3x4 array and its transpose
a34 = np.arange(12, dtype=np.float64).reshape(3, 4)
assert (a34.strides, a34.T.shape, a34.T.strides) == ((32, 8), (4, 3), (8, 32))
assert a34.flags["C_CONTIGUOUS"] and not a34.flags["F_CONTIGUOUS"]
assert a34.T.flags["F_CONTIGUOUS"] and not a34.T.flags["C_CONTIGUOUS"]
raw34 = a34.tobytes()
for i, j in [(0, 0), (1, 2), (2, 3)]:
    off = i * 32 + j * 8
    assert off == j * a34.T.strides[0] + i * a34.T.strides[1]
    assert np.frombuffer(raw34, dtype=np.float64, count=1, offset=off)[0] == a34[i, j] == a34.T[j, i]
assert (2 * 32 + 3 * 8) == 88
assert np.shares_memory(a34, a34.T)
flat = a34.T.reshape(12)
assert flat.strides == (8,) and not np.shares_memory(a34, flat)

# dtypes (notes §3.2)
x8 = np.array([100, 100], dtype=np.int8)
assert (x8 + x8).tolist() == [-56, -56] and x8.sum() == 200 and x8.sum().dtype == np.int64
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    assert np.int8(100) + np.int8(100) == -56
assert any("overflow" in str(w.message) for w in caught)
if NUMPY2:
    try:
        x8 + 200
        raise AssertionError("NumPy 2 should refuse a Python int that does not fit int8")
    except OverflowError:
        pass
assert np.finfo(np.float64).eps == 2.220446049250313e-16 and np.iinfo(np.int8).max == 127
f32 = np.array([0.1, 0.2], dtype=np.float32)
assert (f32 + 1.0).dtype == np.float32 and (f32 + np.float64(1.0)).dtype == np.float64
assert np.result_type(np.int8, np.uint8) == np.int16 and np.result_type(np.int64, np.float32) == np.float64
assert np.array([1.9, -1.9]).astype(int).tolist() == [1, -1]
assert np.array([300, -1]).astype(np.uint8).tolist() == [44, 255]
img = np.array([200, 250], dtype=np.uint8)
assert (img + 100).tolist() == [44, 94]           # a brightened uint8 image wraps to dark pixels

# creation (notes §3.3)
r = np.arange(1, 1.3, 0.1)
assert len(r) == 4 and r[-1] == 1.3000000000000003 and (1.3 - 1) / 0.1 == 3.0000000000000004
assert np.allclose(np.linspace(0, 1, 5, endpoint=False), [0, 0.2, 0.4, 0.6, 0.8])
xg, yg = np.arange(3), np.arange(2)
assert np.meshgrid(xg, yg)[0].shape == (2, 3) and np.meshgrid(xg, yg, indexing="ij")[0].shape == (3, 2)
gy, gx = np.ogrid[0:2, 0:3]
assert (gy.shape, gx.shape) == ((2, 1), (1, 3))
assert rng.integers(0, 6, size=1000).max() <= 5
print(f"Part 1 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 2: views versus copies (Build 2, notes §3.4, Example 4.2)
# Ten operations, predicted then asserted with `np.shares_memory`.

# %%
a = np.arange(12).reshape(3, 4)
VIEW, COPY = True, False
table = {
    "a[1]": (lambda a: a[1], VIEW),
    "a[:, 1]": (lambda a: a[:, 1], VIEW),
    "a[::2, ::2]": (lambda a: a[::2, ::2], VIEW),
    "a[[0, 2]]": (lambda a: a[[0, 2]], COPY),
    "a[a % 2 == 0]": (lambda a: a[a % 2 == 0], COPY),
    "a.T": (lambda a: a.T, VIEW),
    "a.reshape(4, 3)": (lambda a: a.reshape(4, 3), VIEW),
    "a.T.reshape(12)": (lambda a: a.T.reshape(12), COPY),
    "a.ravel()": (lambda a: a.ravel(), VIEW),
    "a.flatten()": (lambda a: a.flatten(), COPY),
}
for name, (f, expected) in table.items():
    assert np.shares_memory(a, f(a)) is expected, name
    print(f"  {name:16s} {'view' if expected else 'copy'}")
assert not np.shares_memory(a, a.T.ravel()) and np.shares_memory(a, a.T)
assert not a.T.reshape(12).flags["OWNDATA"]       # a copy that does not own its data: OWNDATA is not the test

b = a.copy()
v = b[:, 1]
v[:] = -1
assert (b[:, 1] == -1).all()                      # a view writes through
b = a.copy()
f = b[[0, 2]]
f[:] = -1
assert (b[0] == [0, 1, 2, 3]).all()               # a copy does not

mask = a % 3 == 0
assert a[mask].tolist() == [0, 3, 6, 9] and mask.sum() == np.count_nonzero(mask) == 4
assert np.where(mask, a, -1)[0].tolist() == [0, -1, -1, 3]
assert [i.tolist() for i in np.nonzero(mask)] == [[0, 0, 1, 2], [0, 3, 2, 1]]
b = a.copy()
b[b > 8] = 0
assert b[2].tolist() == [8, 0, 0, 0]
assert a[[0, 2], [1, 3]].tolist() == [1, 11] and a[np.ix_([0, 2], [1, 3])].tolist() == [[1, 3], [9, 11]]
c = np.zeros(3, dtype=int)
c[[0, 0, 1]] += 1
d = np.zeros(3, dtype=int)
np.add.at(d, [0, 0, 1], 1)
assert c.tolist() == [1, 1, 0] and d.tolist() == [2, 1, 0]
assert (a[:, None, 0].shape, a[..., 0].shape, a[None].shape) == ((3, 1), (3,), (1, 3, 4))


def center_columns(X):
    """Intended: X minus its column means. The bug: `-=` edits the caller's array."""
    X -= X.mean(axis=0)
    return X


def center_columns_safe(X):
    return X - X.mean(axis=0)


data = np.array([[1.0, 2.0], [3.0, 6.0], [5.0, 10.0]])
train = data[:2]
centered = center_columns(train)
assert data.tolist() == [[-1.0, -2.0], [1.0, 2.0], [5.0, 10.0]]       # Example 4.2: rows 0-1 overwritten
assert np.shares_memory(data, centered) and centered is train
data = np.array([[1.0, 2.0], [3.0, 6.0], [5.0, 10.0]])
assert center_columns_safe(data[:2])[0].tolist() == [-1.0, -2.0] and data[0].tolist() == [1.0, 2.0]
lst = [1, 2, 3]
assert lst[:2] is not lst and lst[:2] == [1, 2]                         # lists copy on slicing
print(f"Part 2 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 3: matrix multiply, three ways (Build 6, notes §3.8)
# $C_{ik}=\sum_j A_{ij}B_{jk}$ with three loops, `np.dot`, `@`. Times for n = 50, 100, 200 (and NumPy alone beyond).

# %%
def matmul_loops(A, B):
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


rng_mm = np.random.default_rng(1)
sizes_loop, t_loops, t_at = (50, 100, 200), [], []
for n in sizes_loop:
    A, B = rng_mm.standard_normal((n, n)), rng_mm.standard_normal((n, n))
    t0 = time.perf_counter()
    C = matmul_loops(A.tolist(), B.tolist())
    t_loops.append(time.perf_counter() - t0)
    assert np.allclose(C, A @ B) and np.allclose(np.dot(A, B), A @ B)
    t_at.append(best(lambda: A @ B, repeat=5))
    print(f"  n={n:4d}: loops {t_loops[-1]:.4f} s, @ {t_at[-1]:.6f} s, ratio {t_loops[-1] / t_at[-1]:,.0f}")
assert t_loops[-1] > 100 * t_at[-1]
sizes_np = (50, 100, 200, 400, 800, 1600)
t_np_big = [best(lambda A=rng_mm.standard_normal((n, n)): A @ A, repeat=2) for n in sizes_np]

fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(sizes_loop, t_loops, "o-", label="three Python loops")
ax.loglog(sizes_np, t_np_big, "s-", label="NumPy @ (BLAS)")
ax.loglog(sizes_loop, [t_loops[0] * (n / 50) ** 3 for n in sizes_loop], "k:", label=r"$\propto n^3$")
ax.set(xlabel="matrix size n", ylabel="seconds", title="Matrix multiply: the same 2n^3 flops")
ax.grid(alpha=0.3, which="both")
ax.legend()
fig.tight_layout()
fig.savefig(FIG / "matmul_timings.png", dpi=120)
plt.close(fig)
print(f"Part 3 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 4: broadcasting by hand, the axis drill, pairwise distances (Builds 3-5, notes §3.5-§3.6)

# %%
def bshape(s, t):
    """Broadcast shape of shapes s and t by rules 1-3 of (3.4); ValueError if incompatible."""
    n = max(len(s), len(t))
    s, t = (1,) * (n - len(s)) + tuple(s), (1,) * (n - len(t)) + tuple(t)
    out = []
    for p, q in zip(s, t):
        if p != q and p != 1 and q != 1:
            raise ValueError(f"shapes {s} and {t} do not broadcast")
        out.append(q if p == 1 else p)
    return tuple(out)


def my_add(a, b):
    """a + b via explicit broadcast_to: every element is a[i or 0, ...] + b[i or 0, ...]."""
    shape = bshape(a.shape, b.shape)
    return np.broadcast_to(a, shape) + np.broadcast_to(b, shape)


pairs = [((3, 4), (4,)), ((3, 4), (3,)), ((3, 1), (4,)), ((2, 1, 3), (4, 1)), ((8, 1, 6, 1), (7, 1, 5))]
answers = [(3, 4), None, (3, 4), (2, 4, 3), (8, 7, 6, 5)]
for (s, t), want in zip(pairs, answers):
    if want is None:
        for fn in (bshape, np.broadcast_shapes):
            try:
                fn(s, t)
                raise AssertionError("should not broadcast")
            except ValueError:
                pass
    else:
        assert bshape(s, t) == np.broadcast_shapes(s, t) == want
        p, q = np.ones(s), np.arange(float(np.prod(t))).reshape(t)
        assert np.array_equal(my_add(p, q), p + q)
try:
    np.ones((3, 4)) + np.ones(3)
except ValueError as err:
    assert "(3,4) (3,)" in str(err)
row = np.arange(3)
Bc = np.broadcast_to(row, (4, 3))
assert Bc.strides == (0, 8) and not Bc.flags["WRITEABLE"] and np.shares_memory(Bc, row)

# Example 4.3
a3 = np.arange(15).reshape(3, 1, 5)
b3 = np.arange(4).reshape(4, 1) * 100
r3 = a3 + b3
assert r3.shape == np.broadcast_shapes(a3.shape, b3.shape) == (3, 4, 5)
assert r3[2, 3, 4] == a3[2, 0, 4] + b3[3, 0] == 314
A3, B3 = np.broadcast_arrays(a3, b3)
assert (A3.strides, B3.strides) == ((40, 0, 8), (0, 8, 0))
try:
    a3 + np.arange(4)
    raise AssertionError("(3,1,5)+(4,) must fail")
except ValueError as err:
    assert "(3,1,5) (4,)" in str(err)
assert (a3 + np.arange(4)[:, None]).shape == (3, 4, 5)

# patterns
xv = np.array([1.0, 2.0, 4.0])
assert (xv[:, None] * xv[None, :]).shape == (3, 3) and (xv[:, None] - xv[None, :])[0].tolist() == [0.0, -1.0, -3.0]
Xc = np.array([[1.0, 10.0], [2.0, 20.0], [6.0, 60.0]])
assert np.allclose((Xc - Xc.mean(axis=0)).mean(axis=0), 0) and np.allclose((Xc - Xc.mean(axis=1, keepdims=True)).mean(axis=1), 0)
assert (Xc - Xc.mean(axis=1, keepdims=True))[0].tolist() == [-4.5, 4.5]
try:
    Xc - Xc.mean(axis=1)
    raise AssertionError("(3,2)-(3,) must fail")
except ValueError:
    pass
assert ((np.arange(3) - np.arange(3)[:, None]).shape) == (3, 3)      # the silent (n,) - (n,1) trap

# Build 5: the axis drill
t = np.arange(24).reshape(2, 3, 4)
drill = {
    "t.sum(axis=0)": (3, 4), "t.sum(axis=1)": (2, 4), "t.sum(axis=2)": (2, 3), "t.sum(axis=(0, 2))": (3,),
    "t.sum(axis=1, keepdims=True)": (2, 1, 4), "t.mean(axis=-1)": (2, 3), "t.argmax(axis=0)": (3, 4),
    "t.cumsum(axis=1)": (2, 3, 4), "(t > 10).any(axis=2)": (2, 3), "t.std(axis=0, ddof=1)": (3, 4),
}
for expr, shape in drill.items():
    assert np.shape(eval(expr)) == shape, expr
assert t.sum() == 276 and t.max() == 23 and t.sum(axis=None) == t.sum()
import statistics

data4 = [1.0, 2.0, 3.0, 4.0]
assert math.isclose(np.std(data4), 1.118033988749895) and math.isclose(np.std(data4, ddof=1), statistics.stdev(data4))
assert math.isclose(statistics.pstdev(data4), np.std(data4))

# Build 4: pairwise distances three ways
def dist_loops(X):
    n = len(X)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = np.sqrt(((X[i] - X[j]) ** 2).sum())
    return D


def dist_broadcast(X):
    """D_ij = ||x_i - x_j||_2 via an (n, n, d) difference array."""
    diff = X[:, None, :] - X[None, :, :]
    return np.sqrt((diff ** 2).sum(axis=-1))


Xp = np.random.default_rng(5).standard_normal((300, 3))
D1, D2, D3 = dist_loops(Xp), dist_broadcast(Xp), cdist(Xp, Xp)
assert np.allclose(D1, D2) and np.allclose(D2, D3) and np.allclose(D1, D1.T) and np.allclose(np.diag(D1), 0)
tl, tb, tc = best(dist_loops, Xp, repeat=1), best(dist_broadcast, Xp), best(lambda Z: cdist(Z, Z), Xp)
assert tl > 10 * tb and tb > tc
assert round(300 * 300 * 3 * 8 / 1e6, 1) == 2.2 and round(10_000 ** 2 * 3 * 8 / 1e9, 1) == 2.4
print(f"  distances: loops {tl * 1e3:.1f} ms, broadcast {tb * 1e3:.2f} ms, cdist {tc * 1e3:.2f} ms")
print(f"Part 4 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 5: ufuncs, sorting, linear algebra, einsum, floating point (notes §3.6-§3.9, Examples 4.4-4.5)

# %%
xu = np.array([1.0, 4.0, 9.0, -1.0])
out = np.zeros_like(xu)
np.sqrt(xu, out=out, where=xu >= 0)
assert out.tolist() == [1.0, 2.0, 3.0, 0.0]
vv = np.array([1, 2, 3, 4])
assert np.add.reduce(vv) == 10 and np.add.accumulate(vv).tolist() == [1, 3, 6, 10]
assert np.multiply.outer(vv, vv)[1].tolist() == [2, 4, 6, 8] and np.maximum.accumulate([3, 1, 4, 1, 5]).tolist() == [3, 3, 4, 4, 5]

# sorting and reshaping (3.7)
v = np.array([30, 10, 20, 10])
assert np.sort(v).tolist() == [10, 10, 20, 30] and np.argsort(v, kind="stable").tolist() == [1, 3, 2, 0]
uniq, counts = np.unique(v, return_counts=True)
assert uniq.tolist() == [10, 20, 30] and counts.tolist() == [2, 1, 1]
s = np.sort(v)
assert np.searchsorted(s, [5, 10, 15, 40]).tolist() == [0, 0, 2, 4] and np.searchsorted(s, 10, side="right") == 2
assert np.partition(v, 1)[:2].tolist() == [10, 10]
m = np.arange(12).reshape(3, -1)
assert m.shape == (3, 4) and np.shares_memory(np.arange(12), m) is False and m.T.shape == (4, 3)
assert np.concatenate([m, m], axis=0).shape == (6, 4) and np.hstack([m, m]).shape == (3, 8) and np.stack([m, m]).shape == (2, 3, 4)
assert [p.shape for p in np.split(m, 2, axis=1)] == [(3, 2), (3, 2)]
assert t.transpose(2, 0, 1).shape == (4, 2, 3) and t.swapaxes(0, 2).shape == (4, 3, 2)
assert not m.ravel().flags["OWNDATA"] and not np.shares_memory(m, m.T.ravel())

# Example 4.4: argsort sorts one array by another
names = np.array(["Ada", "Ben", "Cy", "Di"])
scores = np.array([72, 95, 88, 95])
order = np.argsort(-scores, kind="stable")
assert order.tolist() == [1, 3, 2, 0] and names[order].tolist() == ["Ben", "Di", "Cy", "Ada"] and scores[order].tolist() == [95, 95, 88, 72]
ranks = np.empty_like(order)
ranks[order] = np.arange(len(order))
assert ranks.tolist() == [3, 0, 2, 1]
K = np.array([[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 2.0]])
w2, V = np.linalg.eigh(K)
idx = np.argsort(w2)[::-1]
assert np.allclose(np.sqrt(w2[idx]), [math.sqrt(2 + math.sqrt(2)), math.sqrt(2), math.sqrt(2 - math.sqrt(2))])
assert np.sqrt(w2[idx]).round(4).tolist() == [1.8478, 1.4142, 0.7654]
assert np.allclose(K @ V[:, idx], V[:, idx] * w2[idx])

# linear algebra properties (3.8)
rng_la = np.random.default_rng(3)
A4 = rng_la.standard_normal((4, 4))
b4 = rng_la.standard_normal(4)
x4 = np.linalg.solve(A4, b4)
M, y = rng_la.standard_normal((6, 3)), rng_la.standard_normal(6)
coef, _, rank, _ = np.linalg.lstsq(M, y, rcond=None)
S = A4 + A4.T
w, Vs = np.linalg.eigh(S)
U, sv, Vt = np.linalg.svd(M, full_matrices=False)
Q, R = np.linalg.qr(M)
assert np.allclose(A4 @ x4, b4) and rank == 3 and np.allclose(M.T @ (M @ coef - y), 0)
assert np.allclose(S @ Vs, Vs * w) and np.allclose(Vs.T @ Vs, np.eye(4)) and np.all(np.diff(w) >= 0)
assert (U.shape, sv.shape, Vt.shape) == ((6, 3), (3,), (3, 3)) and np.allclose((U * sv) @ Vt, M) and np.allclose(Q @ R, M)
assert np.allclose(np.linalg.norm(M, axis=0), np.sqrt((M ** 2).sum(axis=0)))
assert np.isclose(np.linalg.det(A4), np.prod(np.linalg.eigvals(A4)).real)


def make_matrix(n, cond, rng):
    """Random n x n matrix with singular values logspace(0, -log10(cond)): condition number cond (LA-9)."""
    U, _ = np.linalg.qr(rng.standard_normal((n, n)))
    V, _ = np.linalg.qr(rng.standard_normal((n, n)))
    return (U * np.logspace(0, -np.log10(cond), n)) @ V.T


rng_inv = np.random.default_rng(0)
res_solve, res_inv = {}, {}
for cond in (1e2, 1e8, 1e12):
    Ac = make_matrix(200, cond, rng_inv)
    x_true = rng_inv.standard_normal(200)
    bc = Ac @ x_true
    x_s, x_i = np.linalg.solve(Ac, bc), np.linalg.inv(Ac) @ bc
    res_solve[cond] = np.linalg.norm(Ac @ x_s - bc) / np.linalg.norm(bc)
    res_inv[cond] = np.linalg.norm(Ac @ x_i - bc) / np.linalg.norm(bc)
    assert res_solve[cond] < 1e-13
assert res_inv[1e12] > 1e-8 and res_inv[1e12] > 1e3 * res_solve[1e12] and res_inv[1e8] > 10 * res_solve[1e8]
Ab, bb = np.random.default_rng(9).standard_normal((1000, 1000)), np.random.default_rng(10).standard_normal(1000)
t_solve, t_inv = best(np.linalg.solve, Ab, bb), best(lambda: np.linalg.inv(Ab) @ bb)
assert t_inv > 1.5 * t_solve
print(f"  residuals at kappa=1e12: solve {res_solve[1e12]:.1e}, inv {res_inv[1e12]:.1e}; n=1000 solve {t_solve * 1e3:.0f} ms, inv {t_inv * 1e3:.0f} ms")

# einsum (3.8) and Example 4.5
rng_e = np.random.default_rng(7)
M3, Ae, Be = rng_e.standard_normal((3, 3)), rng_e.standard_normal((3, 4)), rng_e.standard_normal((4, 2))
Se, Te = rng_e.standard_normal((5, 3, 4)), rng_e.standard_normal((5, 4, 2))
ue, we, Qe, We = rng_e.standard_normal(3), rng_e.standard_normal(4), rng_e.standard_normal((3, 4, 5)), rng_e.standard_normal((4, 5))
cases = [
    (np.einsum("ii->", M3), np.trace(M3), ()),
    (np.einsum("ij,jk->ik", Ae, Be), Ae @ Be, (3, 2)),
    (np.einsum("bij,bjk->bik", Se, Te), Se @ Te, (5, 3, 2)),
    (np.einsum("i,j->ij", ue, we), np.multiply.outer(ue, we), (3, 4)),
    (np.einsum("ijk,jk->i", Qe, We), (Qe * We).sum(axis=(1, 2)), (3,)),
]
for got, want, shape in cases:
    assert np.shape(got) == shape and np.allclose(got, want)
Ai, Bi = np.array([[1, 2], [3, 4]]), np.array([[0, 1], [1, 0]])
assert np.einsum("ij,jk->ik", Ai, Bi).tolist() == [[2, 1], [4, 3]] == (Ai @ Bi).tolist()
assert np.einsum("ii", Ai) == 5 and np.einsum("ii->i", Ai).tolist() == [1, 4] and np.einsum("ij->ji", Ai).tolist() == [[1, 3], [2, 4]]
mm = np.array([1.0, 2.0])
rr = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
I_e = np.einsum("n,na,na->", mm, rr, rr) * np.eye(3) - np.einsum("n,na,nb->ab", mm, rr, rr)
I_loop = sum(mn * (rn @ rn * np.eye(3) - np.outer(rn, rn)) for mn, rn in zip(mm, rr))
assert np.allclose(I_e, I_loop) and np.allclose(np.diag(I_e), [8, 1, 9]) and np.allclose(I_e, np.diag([8.0, 1.0, 9.0]))

# floating point (3.9)
assert (0.1 + 0.2 == 0.3) is False and np.isclose(0.1 + 0.2, 0.3)
assert np.isclose(1e-10, 2e-10) and not np.isclose(1e-10, 2e-10, atol=0)
xf = np.random.default_rng(11).random(10 ** 6).astype(np.float32)
exact = math.fsum(xf.astype(np.float64))
rel = lambda s: abs(float(s) - exact) / exact
e_pair, e_seq, e_64 = rel(np.sum(xf)), rel(np.cumsum(xf)[-1]), rel(np.sum(xf, dtype=np.float64))
assert e_pair < 1e-7 and e_seq > 100 * e_pair and e_64 < 1e-14
if sys.version_info >= (3, 12):
    assert rel(sum(xf.tolist())) < 1e-15           # Python's sum() is compensated from 3.12
print(f"  float32 sum rel. error: pairwise {e_pair:.1e}, sequential {e_seq:.1e}, float64 accumulator {e_64:.1e}")
print(f"Part 5 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 6: the random-walk ensemble (Build 7, notes §1 and §3.10, F5)
# $\langle x_t^2\rangle = t$ for $\pm1$ steps (3.7): a log-log fit of the MSD has slope 1.

# %%
def walks_loop(n_walkers, n_steps, seed=0):
    r = random.Random(seed)
    paths = []
    for _ in range(n_walkers):
        x, path = 0, []
        for _ in range(n_steps):
            x += r.choice((-1, 1))
            path.append(x)
        paths.append(path)
    return paths


def walks_array(n_walkers, n_steps, seed=0):
    """Positions of independent +-1 walkers, shape (walkers, time): one draw, one cumsum along axis 1."""
    gen = np.random.default_rng(seed)
    return gen.choice([-1, 1], size=(n_walkers, n_steps)).cumsum(axis=1)


def walks_array_fast(n_walkers, n_steps, seed=0):
    gen = np.random.default_rng(seed)
    steps = 2 * gen.integers(0, 2, size=(n_walkers, n_steps), dtype=np.int16) - 1
    return steps.cumsum(axis=1, dtype=np.int16)


t_small = best(walks_loop, 200, 1000, repeat=2)
per_step = t_small / (200 * 1000)
t_big, t_fast = best(walks_array, 10_000, 1_000, repeat=2), best(walks_array_fast, 10_000, 1_000, repeat=2)
assert t_big < per_step * 1e7 / 5 and t_fast < t_big
print(f"  loop {per_step * 1e6:.2f} us/step (projection for 1e7 steps: {per_step * 1e7:.1f} s); arrays {t_big:.3f} s, narrow dtype {t_fast:.3f} s")

paths_small = np.array(walks_loop(5, 50, seed=3))
assert paths_small.shape == (5, 50) and (np.abs(np.diff(paths_small, axis=1)) == 1).all()
paths = walks_array_fast(10_000, 1_000, seed=2)
assert paths.shape == (10_000, 1_000) and paths.dtype == np.int16 and (np.abs(paths[:, 0]) == 1).all()
tt = np.arange(1, 1001)
msd = (paths.astype(np.float64) ** 2).mean(axis=0)
slope, intercept = np.polyfit(np.log(tt), np.log(msd), 1)
assert msd.shape == (1000,) and msd[0] == 1.0
assert abs(slope - 1) < 0.05, slope
assert abs(msd[-1] / 1000 - 1) < 5 * math.sqrt(2 / 10_000)          # within 5 standard errors of t
print(f"  MSD(1)={msd[0]:.3f}, MSD(1000)={msd[-1]:.1f}, log-log slope {slope:.4f}, intercept {intercept:.4f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(np.arange(1, 1001), paths[:20].T, lw=0.6, alpha=0.8)
ax1.set(xlabel="step t", ylabel="position", title="20 of the 10,000 walkers")
ax2.loglog(tt, msd, ".", ms=3, label="measured MSD")
ax2.loglog(tt, tt, "k--", label=r"theory $\langle x_t^2\rangle = t$")
ax2.set(xlabel="step t", ylabel="mean-square displacement", title=f"log-log slope = {slope:.3f}")
ax2.grid(alpha=0.3, which="both")
ax2.legend()
fig.tight_layout()
fig.savefig(FIG / "random_walk_msd.png", dpi=120)
plt.close(fig)
print(f"Part 6 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 7: the Mandelbrot set by boolean masks (Example 4.6)

# %%
def mandelbrot_counts(c, max_iter):
    """Iterations survived (|z| <= 2) by z -> z^2 + c from z = 0, capped at max_iter; plus live pixels per step."""
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


def mandelbrot_unmasked(c, max_iter):
    """Same counts without masks: every pixel is updated every step (overflow warnings silenced)."""
    z = np.zeros_like(c)
    count = np.zeros(c.shape, dtype=int)
    with np.errstate(all="ignore"):
        for _ in range(max_iter):
            z = z * z + c
            count += np.abs(z) <= 2
    return count


c4 = np.array([0, 1, 1 + 1j, -2])
z4, alive4, count4 = np.zeros_like(c4), np.ones(4, dtype=bool), np.zeros(4, dtype=int)
trace = []
for k in range(4):
    z4[alive4] = z4[alive4] ** 2 + c4[alive4]
    alive4 &= np.abs(z4) <= 2
    count4 += alive4
    trace.append((alive4.astype(int).tolist(), count4.tolist()))
assert trace[0] == ([1, 1, 1, 1], [1, 1, 1, 1]) and trace[1] == ([1, 1, 0, 1], [2, 2, 1, 2])
assert trace[2] == ([1, 0, 0, 1], [3, 2, 1, 3]) and trace[3] == ([1, 0, 0, 1], [4, 2, 1, 4])
assert z4.tolist() == [0, 5, 1 + 3j, 2] and abs(1 + 3j) == math.sqrt(10)

xs_m, ys_m = np.linspace(-2.2, 0.8, 600), np.linspace(-1.2, 1.2, 400)
C = xs_m[None, :] + 1j * ys_m[:, None]
t_mask = best(mandelbrot_counts, C, 80, repeat=2)
counts, n_alive = mandelbrot_counts(C, 80)
t_unmask = best(mandelbrot_unmasked, C, 80, repeat=2)
assert np.array_equal(counts, mandelbrot_unmasked(C, 80)) and 0.2 < t_mask / t_unmask < 5
assert C.shape == counts.shape == (400, 600)
assert [n_alive[k - 1] for k in (1, 3, 10, 80)] == [213584, 137632, 67222, 51626]
assert mandelbrot_counts(np.array([0, 1, 1 + 1j]), 80)[0].tolist() == [80, 2, 1]        # the three known points
assert mandelbrot_counts(np.array([-2 + 0j, -1 + 0j, 0.25 + 0j]), 80)[0].tolist() == [80, 80, 80]     # in the set

fig, ax = plt.subplots(figsize=(7, 4.8))
im = ax.imshow(counts, extent=(-2.2, 0.8, -1.2, 1.2), origin="lower", cmap="magma")
ax.set(xlabel="Re c", ylabel="Im c", title="Mandelbrot escape counts (80 iterations)")
fig.colorbar(im, ax=ax, label="iterations survived")
fig.tight_layout()
fig.savefig(FIG / "mandelbrot.png", dpi=120)
plt.close(fig)
print(f"  masked {t_mask:.2f} s, unmasked {t_unmask:.2f} s")
print(f"Part 7 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 8: windows, strides and convolution (Build 8, notes §3.11)
# Four 2D convolutions (3.8): two loops, windows plus a product-and-sum, windows plus `einsum`, `scipy.signal.convolve2d`.

# %%
x8w = np.arange(8.0)
Wd = sliding_window_view(x8w, 3)
assert Wd.shape == (6, 3) and Wd.strides == (8, 8) and np.shares_memory(Wd, x8w) and not Wd.flags["WRITEABLE"]
assert np.allclose(Wd.mean(axis=-1), np.convolve(x8w, np.ones(3) / 3, mode="valid"))
assert np.array_equal(Wd, as_strided(x8w, shape=(6, 3), strides=(8, 8)))
sig = np.random.default_rng(6).standard_normal(1000)
for w in (3, 10, 51):
    ma_win = sliding_window_view(sig, w).mean(axis=-1)
    assert ma_win.shape == (1000 - w + 1,) and np.allclose(ma_win, np.convolve(sig, np.ones(w) / w, mode="valid"))


def conv2d_loops(I, K):
    h, w = K.shape
    Kf = K[::-1, ::-1]
    out = np.empty((I.shape[0] - h + 1, I.shape[1] - w + 1))
    for i in range(out.shape[0]):
        for j in range(out.shape[1]):
            out[i, j] = (I[i:i + h, j:j + w] * Kf).sum()
    return out


def conv2d_sum(I, K):
    return (sliding_window_view(I, K.shape) * K[::-1, ::-1]).sum(axis=(-2, -1))


def conv2d_windows(I, K):
    """(I*K)_ij = sum_pq I_{i+p,j+q} K_{h-1-p,w-1-q} (3.8) with no Python loop."""
    return np.einsum("ijpq,pq->ij", sliding_window_view(I, K.shape), K[::-1, ::-1])


rng_c = np.random.default_rng(2)
Im, Ker = rng_c.standard_normal((200, 200)), rng_c.standard_normal((5, 5))
ref = convolve2d(Im, Ker, mode="valid")
times_c = {}
for name, f in [("loops", conv2d_loops), ("window-sum", conv2d_sum), ("einsum", conv2d_windows), ("scipy", lambda I, K: convolve2d(I, K, mode="valid"))]:
    got = f(Im, Ker)
    assert got.shape == (196, 196) and np.allclose(got, ref)
    times_c[name] = best(f, Im, Ker, repeat=2)
assert times_c["loops"] > 10 * times_c["einsum"]
corr = convolve2d(Im, Ker[::-1, ::-1], mode="valid")                    # without the flip: cross-correlation
assert np.allclose(corr, np.einsum("ijpq,pq->ij", sliding_window_view(Im, (5, 5)), Ker))
print("  conv2d ms: " + ", ".join(f"{k} {v * 1e3:.2f}" for k, v in times_c.items()))
print(f"Part 8 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 9: performance rules and I/O (notes §3.10, §3.12)

# %%
def relu_scalar(v):
    return v if v > 0 else 0.0


xr = np.random.default_rng(0).standard_normal(10 ** 6)
t_vec = best(np.vectorize(relu_scalar), xr, repeat=2)
t_lc = best(lambda: [relu_scalar(v) for v in xr.tolist()], repeat=2)
t_max = best(np.maximum, xr, 0.0)
assert np.array_equal(np.vectorize(relu_scalar)(xr), np.maximum(xr, 0.0)) and np.array_equal(np.where(xr > 0, xr, 0.0), np.maximum(xr, 0.0))
assert t_vec > 0.3 * t_lc and t_vec > 10 * t_max                   # np.vectorize is a loop, not a speedup


def grow_append(n):
    a = np.array([])
    for i in range(n):
        a = np.append(a, i)
    return a


def preallocate(n):
    a = np.empty(n)
    for i in range(n):
        a[i] = i
    return a


assert np.array_equal(grow_append(20_000), preallocate(20_000))
assert best(grow_append, 20_000, repeat=2) > 5 * best(preallocate, 20_000, repeat=2)

Xn = np.random.default_rng(1).random((100_000, 10))


def normalize_loop(X):
    out = np.empty_like(X)
    for i in range(len(X)):
        out[i] = X[i] / X[i].sum()
    return out


assert np.allclose(normalize_loop(Xn), Xn / Xn.sum(axis=1, keepdims=True))
assert best(normalize_loop, Xn, repeat=1) > 10 * best(lambda: Xn / Xn.sum(axis=1, keepdims=True))

arr_io = np.random.default_rng(4).standard_normal((3, 2))
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    np.save(tmp / "a.npy", arr_io)
    np.savez(tmp / "bundle.npz", arr=arr_io, ids=np.arange(3))
    np.savetxt(tmp / "a.csv", arr_io, delimiter=",")
    np.savetxt(tmp / "short.csv", arr_io, delimiter=",", fmt="%.3f")
    with np.load(tmp / "bundle.npz") as z:
        assert z.files == ["arr", "ids"] and np.array_equal(z["arr"], arr_io) and z["ids"].tolist() == [0, 1, 2]
    assert np.array_equal(np.load(tmp / "a.npy"), arr_io)                                  # .npy: exact
    assert np.array_equal(np.loadtxt(tmp / "a.csv", delimiter=","), arr_io)                # default %.18e: exact
    short = np.loadtxt(tmp / "short.csv", delimiter=",")
    assert not np.array_equal(short, arr_io) and np.allclose(short, arr_io, atol=5e-4)     # %.3f: lossy
print(f"Part 9 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Your turn
# Reference solutions are in `solutions.md` (B7); try first. Set each `YOUR_...` name below to your
# own function and re-run; the cells check it against the provided reference.
#
# 1. **Vectorize a loop.** `loopy_assign(X, C)` returns, for each row of `X` (n, d), the index of the
#    nearest row of `C` (k, d). Write `assign_vec(X, C)` with no Python loop (hint: Build 4 plus `argmin`).
# 2. **Covariance by broadcasting.** Write `cov_broadcast(X)` for an (n, d) data matrix: centre with a
#    broadcast, then one `@`; divide by n - 1. Compare with `np.cov(X, rowvar=False)`.
# 3. **Game of Life.** Write `life_step(grid)` for a 0/1 array with dead borders using
#    `sliding_window_view` on a padded copy; check a blinker oscillates with period 2 and a glider moves.
# 4. **Find the bug.** `row_normalize_buggy` below passes silently on a square matrix. Explain, then fix.

# %%
def loopy_assign(X, C):
    """For each row x of X, the index of the nearest row of C (squared Euclidean distance)."""
    out = []
    for x in X:
        best_k, best_d = 0, float("inf")
        for k, c in enumerate(C):
            d = sum((x[j] - c[j]) ** 2 for j in range(len(x)))
            if d < best_d:
                best_k, best_d = k, d
        out.append(best_k)
    return np.array(out)


def row_normalize_buggy(X):
    """Intended: divide each row by its sum, so rows sum to 1."""
    return X / X.sum(axis=0)


YOUR_ASSIGN = None          # e.g. YOUR_ASSIGN = assign_vec
YOUR_COV = None             # e.g. YOUR_COV = cov_broadcast
YOUR_LIFE = None            # e.g. YOUR_LIFE = life_step
YOUR_ROW_NORMALIZE = None   # your fixed version of row_normalize_buggy

Xa, Ca = np.random.default_rng(20).standard_normal((60, 3)), np.random.default_rng(21).standard_normal((5, 3))
ref_assign = loopy_assign(Xa, Ca)
if YOUR_ASSIGN is not None:
    assert np.array_equal(YOUR_ASSIGN(Xa, Ca), ref_assign)
    print("Your turn 1 passed")
Xcov = np.random.default_rng(22).standard_normal((50, 4))
if YOUR_COV is not None:
    assert np.allclose(YOUR_COV(Xcov), np.cov(Xcov, rowvar=False))
    print("Your turn 2 passed")
blinker = np.zeros((5, 5), dtype=int)
blinker[2, 1:4] = 1
if YOUR_LIFE is not None:
    once = YOUR_LIFE(blinker)
    assert once[1:4, 2].tolist() == [1, 1, 1] and once.sum() == 3 and np.array_equal(YOUR_LIFE(once), blinker)
    print("Your turn 3 passed")
Xsq = np.arange(1.0, 10.0).reshape(3, 3)
bad = row_normalize_buggy(Xsq)
assert not np.allclose(bad.sum(axis=1), 1)                       # silently wrong on a square matrix
try:
    row_normalize_buggy(np.ones((3, 4)))
    raise AssertionError("expected a broadcast error")
except ValueError:
    pass                                                          # ...and a loud error on a non-square one
if YOUR_ROW_NORMALIZE is not None:
    assert np.allclose(YOUR_ROW_NORMALIZE(Xsq).sum(axis=1), 1) and np.allclose(YOUR_ROW_NORMALIZE(np.ones((3, 4))).sum(axis=1), 1)
    print("Your turn 4 passed")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.perf_counter() - T0:.1f} s.")
