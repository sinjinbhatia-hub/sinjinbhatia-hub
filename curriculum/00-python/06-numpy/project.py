# %% [markdown]
# # PY-6 Project: the vectorization lab
#
# Builds, from raw bytes up, the machinery of notes §3 and asserts every number quoted in §1-§4:
# addressing and strides (Build 1), views versus copies (2), broadcasting by hand (3), pairwise
# distances three ways (4), the axis drill (5), matrix multiply (6), the random-walk ensemble and its
# mean-square displacement (7), 2D convolution (8), and the Mandelbrot set by boolean masks (Example 4.6).
# Timings are measured at run time; asserts use margins several times wider than the notes' numbers.
# Figures go to `figures/`. Run headless: `MPLBACKEND=Agg python3 project.py` (about 15 s).

# %%
import math
import random
import statistics
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

Shape = tuple[int, ...]
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.perf_counter()
NUMPY2 = int(np.__version__.split(".")[0]) >= 2


def best(f, *args, repeat: int = 3) -> float:
    """Smallest wall time (s) of `repeat` calls of f(*args)."""
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        f(*args)
        times.append(time.perf_counter() - t0)
    return min(times)


# %% [markdown]
# ## Part 1: memory, strides, dtypes (notes §1, §3.1-§3.3, Example 4.1)
# Build 1: write equations (3.1) and (3.2) yourself, then read raw bytes with them.

# %%
def c_strides(shape: Shape, itemsize: int) -> Shape:
    """Strides in bytes of a C-contiguous array: s_{d-1} = itemsize, s_k = s_{k+1} n_{k+1} (3.2)."""
    strides, step = [], itemsize
    for n_k in reversed(shape):
        strides.append(step)
        step *= n_k
    return tuple(reversed(strides))


def read_element(raw: bytes, dtype: np.dtype, strides: Shape, idx: Shape, start: int = 0) -> np.generic:
    """Element idx read from bytes `raw`: byte offset = start + sum_k i_k s_k (3.1)."""
    offset = start + sum(i * s for i, s in zip(idx, strides))
    return np.frombuffer(raw, dtype=dtype, count=1, offset=offset)[0]


def sumsq_loop(v: list[float]) -> float:
    s = 0.0
    for t in v:
        s += t * t
    return s


xs, arr = [float(i) for i in range(10 ** 6)], np.arange(10 ** 6, dtype=np.float64)
assert (sys.getsizeof(xs[1]), arr.itemsize, arr.nbytes) == (24, 8, 8_000_000) and np.isclose(sumsq_loop(xs), np.dot(arr, arr))
t_loop, t_np = best(sumsq_loop, xs), best(lambda v: np.dot(v, v), arr)
assert t_loop > 20 * t_np
print(f"sum of squares, 10^6 floats: loop {t_loop * 1e3:.1f} ms, np.dot {t_np * 1e3:.3f} ms ({t_loop / t_np:.0f}x)")

a = np.arange(6).reshape(2, 3)
assert (a.shape, a.ndim, a.size, a.dtype, a.itemsize, a.nbytes, a.strides) == ((2, 3), 2, 6, np.dtype("int64"), 8, 48, (24, 8))
assert a.flags["C_CONTIGUOUS"] and not a.flags["OWNDATA"] and a.base is not None
assert (a.T.strides, a[:, 1].strides, a[:, ::2].strides, a[::-1].strides) == ((8, 24), (24,), (24, 16), (-24, 8))
assert c_strides(a.shape, 8) == a.strides and c_strides((2, 3, 4), 8) == np.zeros((2, 3, 4)).strides
raw = a.base.tobytes()
for view in (a, a.T, a[::2, 1:], a[::-1]):
    start = view.ctypes.data - a.ctypes.data
    assert all(read_element(raw, a.dtype, view.strides, i, start) == view[i] for i in np.ndindex(view.shape))

# Example 4.1: the 3x4 array and its transpose
a34 = np.arange(12, dtype=np.float64).reshape(3, 4)
assert (a34.strides, a34.T.shape, a34.T.strides) == ((32, 8), (4, 3), (8, 32))
flags = lambda z: (z.flags["C_CONTIGUOUS"], z.flags["F_CONTIGUOUS"])
assert flags(a34) == (True, False) and flags(a34.T) == (False, True)
for i, j in [(0, 0), (1, 2), (2, 3)]:
    off = i * 32 + j * 8
    assert off == j * a34.T.strides[0] + i * a34.T.strides[1]
    assert np.frombuffer(a34.tobytes(), dtype=np.float64, count=1, offset=off)[0] == a34[i, j] == a34.T[j, i]
assert 2 * 32 + 3 * 8 == 88 and np.shares_memory(a34, a34.T)
assert a34.T.reshape(12).strides == (8,) and not np.shares_memory(a34, a34.T.reshape(12))

# dtypes (3.2): arrays wrap silently, scalars warn, NumPy 2 refuses a Python int that does not fit
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
f32 = np.array([0.1, 0.2], dtype=np.float32)
assert (f32 + 1.0).dtype == np.float32 and (f32 + np.float64(1.0)).dtype == np.float64
assert (np.result_type(np.int8, np.uint8), np.result_type(np.int64, np.float32)) == (np.int16, np.float64)
assert np.array([1.9, -1.9]).astype(int).tolist() == [1, -1] and np.array([300, -1]).astype(np.uint8).tolist() == [44, 255]
# creation (3.3): the float-step arange pitfall
r = np.arange(1, 1.3, 0.1)
assert len(r) == 4 and r[-1] == 1.3000000000000003 and (1.3 - 1) / 0.1 == 3.0000000000000004
assert tuple(g.shape for g in np.ogrid[0:2, 0:3]) == ((2, 1), (1, 3))
print(f"Part 1 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 2: views versus copies (Build 2, notes §3.4, Example 4.2)
# Ten operations, predicted then asserted with `np.shares_memory`.

# %%
a = np.arange(12).reshape(3, 4)
table = {  # operation: does the result share memory with a?
    "a[1]": (lambda a: a[1], True), "a[:, 1]": (lambda a: a[:, 1], True), "a[::2, ::2]": (lambda a: a[::2, ::2], True),
    "a[[0, 2]]": (lambda a: a[[0, 2]], False), "a[a % 2 == 0]": (lambda a: a[a % 2 == 0], False),
    "a.T": (lambda a: a.T, True), "a.reshape(4, 3)": (lambda a: a.reshape(4, 3), True),
    "a.T.reshape(12)": (lambda a: a.T.reshape(12), False), "a.ravel()": (lambda a: a.ravel(), True),
    "a.flatten()": (lambda a: a.flatten(), False),
}
for name, (f, is_view) in table.items():
    assert np.shares_memory(a, f(a)) is is_view, name
    print(f"  {name:16s} {'view' if is_view else 'copy'}")
assert not np.shares_memory(a, a.T.ravel()) and not a.T.reshape(12).flags["OWNDATA"]    # OWNDATA is not the test

b = a.copy()
b[:, 1][:] = -1                                    # through a view: writes through
assert (b[:, 1] == -1).all()
b = a.copy()
b[[0, 2]][:] = -1                                  # through a copy: only the temporary changes
assert (b == a).all()
mask = a % 3 == 0
assert a[mask].tolist() == [0, 3, 6, 9] and mask.sum() == np.count_nonzero(mask) == 4
assert np.where(mask, a, -1)[0].tolist() == [0, -1, -1, 3] and [i.tolist() for i in np.nonzero(mask)] == [[0, 0, 1, 2], [0, 3, 2, 1]]
b = a.copy()
b[b > 8] = 0
assert b[2].tolist() == [8, 0, 0, 0] and a[[0, 2], [1, 3]].tolist() == [1, 11] and a[np.ix_([0, 2], [1, 3])].tolist() == [[1, 3], [9, 11]]
c, d = np.zeros(3, dtype=int), np.zeros(3, dtype=int)
c[[0, 0, 1]] += 1
np.add.at(d, [0, 0, 1], 1)
assert c.tolist() == [1, 1, 0] and d.tolist() == [2, 1, 0]
assert (a[:, None, 0].shape, a[..., 0].shape, a[None].shape) == ((3, 1), (3,), (1, 3, 4))


def center_columns(X: np.ndarray) -> np.ndarray:
    """Intended: X minus its column means. The bug: `-=` edits the caller's array."""
    X -= X.mean(axis=0)
    return X


data = np.array([[1.0, 2.0], [3.0, 6.0], [5.0, 10.0]])
train = data[:2]
centered = center_columns(train)
assert data.tolist() == [[-1.0, -2.0], [1.0, 2.0], [5.0, 10.0]]       # Example 4.2: rows 0-1 overwritten
assert np.shares_memory(data, centered) and centered is train
data = np.array([[1.0, 2.0], [3.0, 6.0], [5.0, 10.0]])
assert (data[:2] - data[:2].mean(axis=0))[0].tolist() == [-1.0, -2.0] and data[0].tolist() == [1.0, 2.0]
assert [1, 2, 3][:2] == [1, 2]                                        # lists copy on slicing, arrays alias
print(f"Part 2 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 3: matrix multiply, three ways (Build 6, notes §3.8)
# $C_{ik}=\sum_j A_{ij}B_{jk}$ with three loops, `np.dot`, `@`; n = 50, 100, 200, then NumPy alone up to 1600.

# %%
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
assert t_loops[-1] > 30 * t_at[-1]
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
    """a + b via explicit broadcast_to."""
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

a3, b3 = np.arange(15).reshape(3, 1, 5), np.arange(4).reshape(4, 1) * 100          # Example 4.3
r3 = a3 + b3
assert r3.shape == np.broadcast_shapes(a3.shape, b3.shape) == (3, 4, 5) and r3[2, 3, 4] == a3[2, 0, 4] + b3[3, 0] == 314
assert tuple(z.strides for z in np.broadcast_arrays(a3, b3)) == ((40, 0, 8), (0, 8, 0))
try:
    a3 + np.arange(4)
    raise AssertionError("(3,1,5)+(4,) must fail")
except ValueError as err:
    assert "(3,1,5) (4,)" in str(err)
assert (a3 + np.arange(4)[:, None]).shape == (3, 4, 5)

Xc = np.array([[1.0, 10.0], [2.0, 20.0], [6.0, 60.0]])                           # centring patterns
assert np.allclose((Xc - Xc.mean(axis=0)).mean(axis=0), 0) and (Xc - Xc.mean(axis=1, keepdims=True))[0].tolist() == [-4.5, 4.5]
assert (np.arange(3) - np.arange(3)[:, None]).shape == (3, 3)                    # the silent (n,) - (n,1) trap
try:
    Xc - Xc.mean(axis=1)
    raise AssertionError("(3,2)-(3,) must fail")
except ValueError:
    pass

t = np.arange(24).reshape(2, 3, 4)                                                # Build 5: the axis drill
drill = {"t.sum(axis=0)": (3, 4), "t.sum(axis=1)": (2, 4), "t.sum(axis=2)": (2, 3), "t.sum(axis=(0, 2))": (3,),
         "t.sum(axis=1, keepdims=True)": (2, 1, 4), "t.mean(axis=-1)": (2, 3), "t.argmax(axis=0)": (3, 4),
         "t.cumsum(axis=1)": (2, 3, 4), "(t > 10).any(axis=2)": (2, 3), "t.std(axis=0, ddof=1)": (3, 4)}
assert all(np.shape(eval(expr)) == shape for expr, shape in drill.items())
assert (t.sum(), t.max()) == (276, 23)
d4 = [1.0, 2.0, 3.0, 4.0]
assert math.isclose(np.std(d4), 1.118033988749895) and math.isclose(np.std(d4, ddof=1), statistics.stdev(d4))


def dist_loops(X: np.ndarray) -> np.ndarray:
    n = len(X)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = np.sqrt(((X[i] - X[j]) ** 2).sum())
    return D


def dist_broadcast(X: np.ndarray) -> np.ndarray:
    """D_ij = ||x_i - x_j||_2 via an (n, n, d) difference array."""
    diff = X[:, None, :] - X[None, :, :]
    return np.sqrt((diff ** 2).sum(axis=-1))


Xp = np.random.default_rng(5).standard_normal((300, 3))                           # Build 4
D1, D2, D3 = dist_loops(Xp), dist_broadcast(Xp), cdist(Xp, Xp)
assert np.allclose(D1, D2) and np.allclose(D2, D3) and np.allclose(D1, D1.T) and np.allclose(np.diag(D1), 0)
tl, tb, tc = best(dist_loops, Xp, repeat=1), best(dist_broadcast, Xp), best(lambda Z: cdist(Z, Z), Xp)
assert tl > 10 * tb and tb > tc
assert round(300 * 300 * 3 * 8 / 1e6, 1) == 2.2 and round(10_000 ** 2 * 3 * 8 / 1e9, 1) == 2.4
print(f"  distances: loops {tl * 1e3:.1f} ms, broadcast {tb * 1e3:.2f} ms, cdist {tc * 1e3:.2f} ms")
print(f"Part 4 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 5: sorting, linear algebra, einsum, floating point (notes §3.6-§3.9, Examples 4.4-4.5)

# %%
names, scores = np.array(["Ada", "Ben", "Cy", "Di"]), np.array([72, 95, 88, 95])          # Example 4.4
order = np.argsort(-scores, kind="stable")
ranks = np.empty_like(order)
ranks[order] = np.arange(len(order))
assert order.tolist() == [1, 3, 2, 0] and names[order].tolist() == ["Ben", "Di", "Cy", "Ada"] and ranks.tolist() == [3, 0, 2, 1]
K = np.array([[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 2.0]])
w2, V = np.linalg.eigh(K)
idx = np.argsort(w2)[::-1]
assert np.allclose(np.sqrt(w2[idx]), [math.sqrt(2 + math.sqrt(2)), math.sqrt(2), math.sqrt(2 - math.sqrt(2))])
assert np.sqrt(w2[idx]).round(4).tolist() == [1.8478, 1.4142, 0.7654] and np.allclose(K @ V[:, idx], V[:, idx] * w2[idx])

rng_la = np.random.default_rng(3)                                                         # linalg properties (3.8)
A4, b4 = rng_la.standard_normal((4, 4)), rng_la.standard_normal(4)
M, y = rng_la.standard_normal((6, 3)), rng_la.standard_normal(6)
coef, _, rank, _ = np.linalg.lstsq(M, y, rcond=None)
w, Vs = np.linalg.eigh(A4 + A4.T)
U, sv, Vt = np.linalg.svd(M, full_matrices=False)
Q, R = np.linalg.qr(M)
assert np.allclose(A4 @ np.linalg.solve(A4, b4), b4) and rank == 3 and np.allclose(M.T @ (M @ coef - y), 0)
assert np.allclose((A4 + A4.T) @ Vs, Vs * w) and np.allclose(Vs.T @ Vs, np.eye(4)) and np.all(np.diff(w) >= 0)
assert np.allclose((U * sv) @ Vt, M) and np.allclose(Q @ R, M) and np.isclose(np.linalg.det(A4), np.prod(np.linalg.eigvals(A4)).real)


def make_matrix(n: int, cond: float, rng: np.random.Generator) -> np.ndarray:
    """Random n x n matrix with singular values logspace(0, -log10(cond)): condition number cond (LA-9)."""
    U, _ = np.linalg.qr(rng.standard_normal((n, n)))
    V, _ = np.linalg.qr(rng.standard_normal((n, n)))
    return (U * np.logspace(0, -np.log10(cond), n)) @ V.T


rng_inv, res_solve, res_inv = np.random.default_rng(0), {}, {}
for cond in (1e2, 1e8, 1e12):
    Ac = make_matrix(200, cond, rng_inv)
    bc = Ac @ rng_inv.standard_normal(200)
    res_solve[cond] = np.linalg.norm(Ac @ np.linalg.solve(Ac, bc) - bc) / np.linalg.norm(bc)
    res_inv[cond] = np.linalg.norm(Ac @ (np.linalg.inv(Ac) @ bc) - bc) / np.linalg.norm(bc)
assert max(res_solve.values()) < 1e-13 and res_inv[1e12] > 1e-8 and res_inv[1e12] > 1e3 * res_solve[1e12] and res_inv[1e8] > 10 * res_solve[1e8]
Ab, bb = np.random.default_rng(9).standard_normal((1000, 1000)), np.random.default_rng(10).standard_normal(1000)
t_solve, t_inv = best(np.linalg.solve, Ab, bb), best(lambda: np.linalg.inv(Ab) @ bb)
assert t_inv > 1.5 * t_solve          # inv then @ costs several times one solve
print(f"  residuals at kappa=1e12: solve {res_solve[1e12]:.1e}, inv {res_inv[1e12]:.1e}; n=1000 solve {t_solve * 1e3:.0f} ms, inv {t_inv * 1e3:.0f} ms")

rng_e = np.random.default_rng(7)                                                          # einsum (3.8), Example 4.5
M3, Ae, Be = rng_e.standard_normal((3, 3)), rng_e.standard_normal((3, 4)), rng_e.standard_normal((4, 2))
Se, Te = rng_e.standard_normal((5, 3, 4)), rng_e.standard_normal((5, 4, 2))
ue, we, Qe, We = rng_e.standard_normal(3), rng_e.standard_normal(4), rng_e.standard_normal((3, 4, 5)), rng_e.standard_normal((4, 5))
for got, want, shape in [(np.einsum("ii->", M3), np.trace(M3), ()), (np.einsum("ij,jk->ik", Ae, Be), Ae @ Be, (3, 2)),
                         (np.einsum("bij,bjk->bik", Se, Te), Se @ Te, (5, 3, 2)), (np.einsum("i,j->ij", ue, we), np.multiply.outer(ue, we), (3, 4)),
                         (np.einsum("ijk,jk->i", Qe, We), (Qe * We).sum(axis=(1, 2)), (3,))]:
    assert np.shape(got) == shape and np.allclose(got, want)
Ai, Bi = np.array([[1, 2], [3, 4]]), np.array([[0, 1], [1, 0]])
assert np.einsum("ij,jk->ik", Ai, Bi).tolist() == [[2, 1], [4, 3]] == (Ai @ Bi).tolist()
assert np.einsum("ii", Ai) == 5 and np.einsum("ii->i", Ai).tolist() == [1, 4] and np.einsum("ij->ji", Ai).tolist() == [[1, 3], [2, 4]]
mm, rr = np.array([1.0, 2.0]), np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])             # inertia tensor of two point masses (CM-7)
I_e = np.einsum("n,na,na->", mm, rr, rr) * np.eye(3) - np.einsum("n,na,nb->ab", mm, rr, rr)
I_loop = sum(mn * (rn @ rn * np.eye(3) - np.outer(rn, rn)) for mn, rn in zip(mm, rr))
assert np.allclose(I_e, I_loop) and np.allclose(I_e, np.diag([8.0, 1.0, 9.0]))

assert (0.1 + 0.2 == 0.3) is False and np.isclose(0.1 + 0.2, 0.3)                          # floating point (3.9)
assert np.isclose(1e-10, 2e-10) and not np.isclose(1e-10, 2e-10, atol=0)
xf = np.random.default_rng(11).random(10 ** 6).astype(np.float32)
exact = math.fsum(xf.astype(np.float64))
rel = lambda s: abs(float(s) - exact) / exact
e_pair, e_seq, e_64 = rel(np.sum(xf)), rel(np.cumsum(xf)[-1]), rel(np.sum(xf, dtype=np.float64))
assert e_pair < 1e-7 and e_seq > 100 * e_pair and e_64 < 1e-14
if sys.version_info >= (3, 12):
    assert rel(sum(xf.tolist())) < 1e-15                                                  # Python's sum() is compensated from 3.12
print(f"  float32 sum, relative error: pairwise {e_pair:.1e}, sequential {e_seq:.1e}, float64 accumulator {e_64:.1e}")
print(f"Part 5 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 6: the random-walk ensemble (Build 7, notes §1 and §3.10, F5)
# $\langle x_t^2\rangle = t$ for $\pm1$ steps (3.7), so the MSD has log-log slope 1.

# %%
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
    """Positions of independent +-1 walkers, shape (walkers, time): one draw, one cumsum along axis 1."""
    return np.random.default_rng(seed).choice([-1, 1], size=(n_walkers, n_steps)).cumsum(axis=1)


def walks_array_fast(n_walkers: int, n_steps: int, seed: int = 0) -> np.ndarray:
    steps = 2 * np.random.default_rng(seed).integers(0, 2, size=(n_walkers, n_steps), dtype=np.int16) - 1
    return steps.cumsum(axis=1, dtype=np.int16)


per_step = best(walks_loop, 200, 1000, repeat=2) / (200 * 1000)
t_big, t_fast = best(walks_array, 10_000, 1_000, repeat=2), best(walks_array_fast, 10_000, 1_000, repeat=2)
assert t_big < per_step * 1e7 / 2 and t_fast < 1.5 * t_big
print(f"  loop {per_step * 1e6:.2f} us/step (projection for 1e7 steps: {per_step * 1e7:.1f} s); arrays {t_big:.3f} s, narrow dtype {t_fast:.3f} s")
small = np.array(walks_loop(5, 50, seed=3))
assert small.shape == (5, 50) and (np.abs(np.diff(small, axis=1)) == 1).all()
paths = walks_array_fast(10_000, 1_000, seed=2)
tt = np.arange(1, 1001)
msd = (paths.astype(np.float64) ** 2).mean(axis=0)
slope, intercept = np.polyfit(np.log(tt), np.log(msd), 1)
assert paths.shape == (10_000, 1_000) and paths.dtype == np.int16 and msd.shape == (1000,) and msd[0] == 1.0
assert abs(slope - 1) < 0.05 and abs(msd[-1] / 1000 - 1) < 5 * math.sqrt(2 / 10_000)
print(f"  MSD(1)={msd[0]:.3f}, MSD(1000)={msd[-1]:.1f}, log-log slope {slope:.4f}, intercept {intercept:.4f}")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(tt, paths[:20].T, lw=0.6, alpha=0.8)
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
def mandelbrot_counts(c: np.ndarray, max_iter: int) -> tuple[np.ndarray, list[int]]:
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


def mandelbrot_unmasked(c: np.ndarray, max_iter: int) -> np.ndarray:
    """Same counts with every pixel updated every step (overflow warnings silenced)."""
    z = np.zeros_like(c)
    count = np.zeros(c.shape, dtype=int)
    with np.errstate(all="ignore"):
        for _ in range(max_iter):
            z = z * z + c
            count += np.abs(z) <= 2
    return count


c4 = np.array([0, 1, 1 + 1j, -2])                                  # the four-point trace of Example 4.6
z4, alive4, count4, trace = np.zeros_like(c4), np.ones(4, dtype=bool), np.zeros(4, dtype=int), []
for k in range(4):
    z4[alive4] = z4[alive4] ** 2 + c4[alive4]
    alive4 &= np.abs(z4) <= 2
    count4 += alive4
    trace.append((alive4.astype(int).tolist(), count4.tolist()))
assert trace == [([1, 1, 1, 1], [1, 1, 1, 1]), ([1, 1, 0, 1], [2, 2, 1, 2]), ([1, 0, 0, 1], [3, 2, 1, 3]), ([1, 0, 0, 1], [4, 2, 1, 4])]
assert z4.tolist() == [0, 5, 1 + 3j, 2] and abs(1 + 3j) == math.sqrt(10)

xs_m, ys_m = np.linspace(-2.2, 0.8, 600), np.linspace(-1.2, 1.2, 400)
C = xs_m[None, :] + 1j * ys_m[:, None]
counts, n_alive = mandelbrot_counts(C, 80)
t_mask, t_unmask = best(mandelbrot_counts, C, 80, repeat=2), best(mandelbrot_unmasked, C, 80, repeat=2)
assert np.array_equal(counts, mandelbrot_unmasked(C, 80)) and 0.2 < t_mask / t_unmask < 5 and counts.shape == (400, 600)
assert [n_alive[k - 1] for k in (1, 3, 10, 80)] == [213584, 137632, 67222, 51626]
assert mandelbrot_counts(np.array([0, 1, 1 + 1j]), 80)[0].tolist() == [80, 2, 1]               # the three known points
assert mandelbrot_counts(np.array([-2 + 0j, -1 + 0j, 0.25 + 0j]), 80)[0].tolist() == [80, 80, 80]   # more points in the set
fig, ax = plt.subplots(figsize=(7, 4.8))
im = ax.imshow(counts, extent=(-2.2, 0.8, -1.2, 1.2), origin="lower", cmap="magma")
ax.set(xlabel="Re c", ylabel="Im c", title="Mandelbrot escape counts (80 iterations)")
fig.colorbar(im, ax=ax, label="iterations survived")
fig.tight_layout()
fig.savefig(FIG / "mandelbrot.png", dpi=120)
plt.close(fig)
print(f"  masked {t_mask:.2f} s, unmasked {t_unmask:.2f} s\nPart 7 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 8: windows, strides and convolution (Build 8, notes §3.11)
# Four 2D convolutions (3.8): two loops, windows plus product-and-sum, windows plus `einsum`, `convolve2d`.

# %%
x8w = np.arange(8.0)
Wd = sliding_window_view(x8w, 3)
assert Wd.shape == (6, 3) and Wd.strides == (8, 8) and np.shares_memory(Wd, x8w) and not Wd.flags["WRITEABLE"]
assert np.array_equal(Wd, as_strided(x8w, shape=(6, 3), strides=(8, 8)))
sig = np.random.default_rng(6).standard_normal(1000)
for w in (3, 10, 51):                                              # strides moving average against np.convolve
    assert np.allclose(sliding_window_view(sig, w).mean(axis=-1), np.convolve(sig, np.ones(w) / w, mode="valid"))


def conv2d_loops(I: np.ndarray, K: np.ndarray) -> np.ndarray:
    h, w = K.shape
    Kf = K[::-1, ::-1]
    out = np.empty((I.shape[0] - h + 1, I.shape[1] - w + 1))
    for i in range(out.shape[0]):
        for j in range(out.shape[1]):
            out[i, j] = (I[i:i + h, j:j + w] * Kf).sum()
    return out


def conv2d_sum(I: np.ndarray, K: np.ndarray) -> np.ndarray:
    return (sliding_window_view(I, K.shape) * K[::-1, ::-1]).sum(axis=(-2, -1))


def conv2d_windows(I: np.ndarray, K: np.ndarray) -> np.ndarray:
    """(I*K)_ij = sum_pq I_{i+p,j+q} K_{h-1-p,w-1-q} (3.8) with no Python loop."""
    return np.einsum("ijpq,pq->ij", sliding_window_view(I, K.shape), K[::-1, ::-1])


rng_c = np.random.default_rng(2)
Im, Ker = rng_c.standard_normal((200, 200)), rng_c.standard_normal((5, 5))
ref = convolve2d(Im, Ker, mode="valid")
times_c = {}
for name, f in [("loops", conv2d_loops), ("window-sum", conv2d_sum), ("einsum", conv2d_windows), ("scipy", lambda I, K: convolve2d(I, K, mode="valid"))]:
    assert f(Im, Ker).shape == (196, 196) and np.allclose(f(Im, Ker), ref)
    times_c[name] = best(f, Im, Ker, repeat=2)
assert times_c["loops"] > 10 * times_c["einsum"]
assert np.allclose(convolve2d(Im, Ker[::-1, ::-1], mode="valid"), np.einsum("ijpq,pq->ij", sliding_window_view(Im, (5, 5)), Ker))   # no flip = correlation
print("  conv2d ms: " + ", ".join(f"{k} {v * 1e3:.2f}" for k, v in times_c.items()))
print(f"Part 8 ok ({time.perf_counter() - T0:.1f} s)")

# %% [markdown]
# ## Part 9: performance rules and I/O (notes §3.10, §3.12)

# %%
def relu_scalar(v: float) -> float:
    return v if v > 0 else 0.0


def grow_append(n: int) -> np.ndarray:
    a = np.array([])
    for i in range(n):
        a = np.append(a, i)
    return a


def preallocate(n: int) -> np.ndarray:
    a = np.empty(n)
    for i in range(n):
        a[i] = i
    return a


def normalize_loop(X: np.ndarray) -> np.ndarray:
    out = np.empty_like(X)
    for i in range(len(X)):
        out[i] = X[i] / X[i].sum()
    return out


xr = np.random.default_rng(0).standard_normal(10 ** 6)
t_vec, t_lc, t_max = best(np.vectorize(relu_scalar), xr, repeat=2), best(lambda: [relu_scalar(v) for v in xr.tolist()], repeat=2), best(np.maximum, xr, 0.0)
assert np.array_equal(np.vectorize(relu_scalar)(xr), np.maximum(xr, 0.0)) and np.array_equal(np.where(xr > 0, xr, 0.0), np.maximum(xr, 0.0))
assert t_vec > 0.3 * t_lc and t_vec > 10 * t_max                   # np.vectorize is a loop, not a speedup
assert np.array_equal(grow_append(20_000), preallocate(20_000)) and best(grow_append, 20_000, repeat=2) > 5 * best(preallocate, 20_000, repeat=2)
Xn = np.random.default_rng(1).random((100_000, 10))
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
# Reference solutions are in `solutions.md` (B7); try first. Set each `YOUR_...` name below to your own
# function and re-run; the cell checks it against the provided reference.
#
# 1. **Vectorize a loop.** `loopy_assign(X, C)` returns, for each row of `X` (n, d), the index of the nearest
#    row of `C` (k, d). Write `assign_vec(X, C)` with no Python loop (hint: Build 4, then `argmin`).
# 2. **Covariance by broadcasting.** `cov_broadcast(X)` for an (n, d) matrix: centre with a broadcast, one `@`,
#    divide by n - 1. Compare with `np.cov(X, rowvar=False)`.
# 3. **Game of Life.** `life_step(grid)` for a 0/1 array with dead borders, using `sliding_window_view` on a
#    padded copy. A blinker must have period 2.
# 4. **Find the bug.** `row_normalize_buggy` passes silently on a square matrix. Explain, then fix.

# %%
def loopy_assign(X: np.ndarray, C: np.ndarray) -> np.ndarray:
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


def row_normalize_buggy(X: np.ndarray) -> np.ndarray:
    """Intended: divide each row by its sum, so every row sums to 1."""
    return X / X.sum(axis=1)


YOUR_ASSIGN = YOUR_COV = YOUR_LIFE = YOUR_ROW_NORMALIZE = None      # set each to your function, e.g. YOUR_ASSIGN = assign_vec

Xa, Ca = np.random.default_rng(20).standard_normal((60, 3)), np.random.default_rng(21).standard_normal((5, 3))
if YOUR_ASSIGN is not None:
    assert np.array_equal(YOUR_ASSIGN(Xa, Ca), loopy_assign(Xa, Ca))
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
assert not np.allclose(row_normalize_buggy(Xsq).sum(axis=1), 1)    # silently wrong on a square matrix...
try:
    row_normalize_buggy(np.ones((3, 4)))
    raise AssertionError("expected a broadcast error")
except ValueError:
    pass                                                            # ...and a loud error on a non-square one
if YOUR_ROW_NORMALIZE is not None:
    assert all(np.allclose(YOUR_ROW_NORMALIZE(X).sum(axis=1), 1) for X in (Xsq, np.ones((3, 4))))
    print("Your turn 4 passed")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.perf_counter() - T0:.1f} s.")
