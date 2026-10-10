# %% [markdown]
# # LA-14 Project: A Reverse-Mode Autodiff Engine, and Every Identity Checked
#
# This script builds the computational side of LA-14. Parts: (1) finite
# differences and the error-vs-h curve, (2) every gradient identity of §3.4
# checked on random inputs, with Examples 4.2 and 4.3, (3) a scalar `Value`
# autodiff class (micrograd-style), (4) a `Tensor` class with matmul,
# broadcasting and the VJP rules of §3.6–3.7, (5) Example 4.4: two-layer
# backprop by hand vs. the engine vs. finite differences, (6) the batch
# version trained on y = sin(x1) + x2^2, (7) Example 4.1: Newton vs gradient
# descent, (8) Example 4.5: Gauss–Newton on an exponential decay fit,
# (9) Gaussian MLE gradients and gradient ascent, (10) forward vs reverse cost.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import least_squares
import sympy as sp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
EPS = np.finfo(float).eps

# %% [markdown]
# ## Part 1: Finite differences (§3.8)
#
# Central difference: (f(x + h e_i) - f(x - h e_i)) / 2h = f_i + O(h^2).
# Truncation error ~ h^2 |f'''|/6, rounding error ~ eps |f| / h, so the best h
# is about eps^(1/3) ~ 6e-6 (central) or eps^(1/2) ~ 1.5e-8 (forward).
# Relative error metric (3.27): ||a - b|| / (||a|| + ||b||).


# %%
def num_grad(f, x, h=None):
    """Central-difference gradient of scalar f at array x (any shape), same shape as x."""
    x = np.asarray(x, float)
    h = EPS ** (1 / 3) if h is None else h
    g = np.zeros_like(x)
    it = np.nditer(x, flags=["multi_index"])
    for _ in it:
        idx = it.multi_index
        xp = x.copy(); xm = x.copy()
        step = h * max(1.0, abs(x[idx]))
        xp[idx] += step; xm[idx] -= step
        g[idx] = (f(xp) - f(xm)) / (2 * step)
    return g


def num_jac(F, x, h=None):
    """Central-difference Jacobian of vector F at vector x: J_ij = dF_i/dx_j (m x n)."""
    x = np.asarray(x, float)
    h = EPS ** (1 / 3) if h is None else h
    cols = []
    for j in range(x.size):
        e = np.zeros_like(x); e[j] = h * max(1.0, abs(x[j]))
        cols.append((np.asarray(F(x + e)) - np.asarray(F(x - e))) / (2 * e[j]))
    return np.column_stack(cols)


def rel_error(a, b):
    """||a - b|| / (||a|| + ||b||), eq. (3.27). Below 1e-7 is a pass for smooth f."""
    a = np.asarray(a, float).ravel(); b = np.asarray(b, float).ravel()
    return np.linalg.norm(a - b) / (np.linalg.norm(a) + np.linalg.norm(b) + 1e-300)


def grad_check(f, grad_f, x, tol=1e-6):
    """Compare analytic grad_f(x) with num_grad(f, x). Returns the relative error."""
    e = rel_error(grad_f(x), num_grad(f, x))
    assert e < tol, f"gradient check failed: rel error {e:.3e}"
    return e


# Error vs h for f(x) = e^x sin x at x = 1 (the §8 plot).
f1 = lambda x: np.exp(x) * np.sin(x)
x0 = 1.0
fp_exact = np.exp(x0) * (np.sin(x0) + np.cos(x0))
hs = np.logspace(-12, 0, 49)
err_c = np.abs((f1(x0 + hs) - f1(x0 - hs)) / (2 * hs) - fp_exact) / abs(fp_exact)
err_f = np.abs((f1(x0 + hs) - f1(x0)) / hs - fp_exact) / abs(fp_exact)
fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(hs, err_c, "o-", ms=3, label="central (f(x+h)-f(x-h))/2h")
ax.loglog(hs, err_f, "s-", ms=3, label="forward (f(x+h)-f(x))/h")
ax.loglog(hs, hs**2 / 6 * abs(np.exp(x0) * 2 * (np.cos(x0) - np.sin(x0))) / abs(fp_exact), "k--", lw=.8, label="h^2 |f'''|/6 (truncation)")
ax.loglog(hs, EPS * abs(f1(x0)) / hs / abs(fp_exact), "k:", lw=.8, label="eps |f|/h (rounding)")
ax.axvline(EPS ** (1 / 3), color="C0", alpha=.3); ax.axvline(EPS ** .5, color="C1", alpha=.3)
ax.set_xlabel("h"); ax.set_ylabel("relative error"); ax.set_title("Finite-difference error vs step size")
ax.legend(fontsize=7); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "fd_error_vs_h.png", dpi=120); plt.close(fig)
h_best_c = hs[err_c.argmin()]; h_best_f = hs[err_f.argmin()]
assert 1e-7 < h_best_c < 1e-4 and err_c.min() < 1e-9, (h_best_c, err_c.min())
assert 1e-10 < h_best_f < 1e-6, h_best_f
print(f"Part 1: best central h = {h_best_c:.1e} (eps^(1/3) = {EPS**(1/3):.1e}); best forward h = {h_best_f:.1e} (sqrt eps = {EPS**.5:.1e})")

# %% [markdown]
# ## Part 2: The identity catalogue (§3.4), Examples 4.2 and 4.3
#
# Every identity derived with differentials in the notes is checked here on
# random inputs by comparing to central differences. The Jacobian ones use
# `num_jac`. Nothing below is allowed to fail.

# %%
n, m = 5, 4
a = rng.standard_normal(n); x = rng.standard_normal(n); y = rng.standard_normal(m)
A = rng.standard_normal((n, n)); S = A + A.T                 # S symmetric
M = rng.standard_normal((m, n)); b = rng.standard_normal(m)   # for ||Mx - b||^2
B = rng.standard_normal((n, n)); P = rng.standard_normal((n, m))
Xmat = rng.standard_normal((n, n)) + 3 * np.eye(n)             # well away from singular

identities = {   # name: (f, analytic gradient, point)
    "(3.9)  grad a^T x = a":                 (lambda v: a @ v,                       lambda v: a,                      x),
    "(3.10) grad x^T A x = (A + A^T) x":     (lambda v: v @ A @ v,                   lambda v: (A + A.T) @ v,          x),
    "(3.10') grad x^T S x = 2 S x (S sym)":  (lambda v: v @ S @ v,                   lambda v: 2 * S @ v,              x),
    "(3.11) grad ||x||^2 = 2x":              (lambda v: v @ v,                       lambda v: 2 * v,                  x),
    "(3.12) grad ||Mx-b||^2 = 2M^T(Mx-b)":   (lambda v: np.sum((M @ v - b) ** 2),    lambda v: 2 * M.T @ (M @ v - b),  x),
    "(3.14) grad_A x^T A y = x y^T":         (lambda W: x @ W @ y,                   lambda W: np.outer(x, y),         P),
    "(3.15) grad_A tr(AB) = B^T":            (lambda W: np.trace(W @ B),             lambda W: B.T,                    Xmat),
    "(3.16) grad_A tr(A^T B A) = (B+B^T)A":  (lambda W: np.trace(W.T @ B @ W),       lambda W: (B + B.T) @ W,          Xmat),
    "(3.17) grad_A log det A = A^{-T}":      (lambda W: np.log(abs(np.linalg.det(W))), lambda W: np.linalg.inv(W).T,   Xmat),
    "(3.18) grad_A ||A||_F^2 = 2A":          (lambda W: np.sum(W ** 2),              lambda W: 2 * W,                  Xmat),
    "(3.19) grad a^T tanh(x) = diag(tanh'(x)) a": (lambda v: a @ np.tanh(v),         lambda v: (1 - np.tanh(v) ** 2) * a, x),
}
for name, (f, g, pt) in identities.items():
    e = grad_check(f, g, pt)
    print(f"  {name:48s} rel err {e:.1e}")

# Jacobians: d(Ax) = A dx; elementwise map has diagonal Jacobian (3.13), (3.19).
assert rel_error(num_jac(lambda v: M @ v, x), M) < 1e-7
assert rel_error(num_jac(np.tanh, x), np.diag(1 - np.tanh(x) ** 2)) < 1e-7


def softmax(z):
    """s_i = e^{z_i} / sum_k e^{z_k}, computed stably."""
    z = np.asarray(z, float); e = np.exp(z - z.max()); return e / e.sum()


def softmax_jacobian(z):
    """J = diag(s) - s s^T, eq. (3.20)."""
    s = softmax(z); return np.diag(s) - np.outer(s, s)


# Example 4.2: softmax Jacobian at s = (0.2, 0.3, 0.5), i.e. z = log s.
z_ex = np.log([0.2, 0.3, 0.5])
J_ex = softmax_jacobian(z_ex)
J_expected = np.array([[0.16, -0.06, -0.10], [-0.06, 0.21, -0.15], [-0.10, -0.15, 0.25]])
assert np.allclose(J_ex, J_expected) and np.allclose(J_ex.sum(1), 0) and np.allclose(J_ex, J_ex.T)
assert rel_error(num_jac(softmax, z_ex), J_ex) < 1e-7
assert np.linalg.matrix_rank(J_ex) == 2 and np.allclose(J_ex @ np.ones(3), 0)  # null vector is 1
print("Example 4.2: softmax Jacobian at s=(0.2,0.3,0.5):\n", np.round(J_ex, 4))
# Cross-entropy with softmax: grad_z (-y^T log s) = s - y, eq. (3.21).
for _ in range(10):
    z = rng.standard_normal(6); yk = rng.integers(6); y1 = np.eye(6)[yk]
    ce = lambda v: -np.log(softmax(v)[yk])
    grad_check(ce, lambda v: softmax(v) - y1, z)
e_ce = np.array([0.2, 0.3, 0.5]) - np.array([0, 0, 1.])
assert np.allclose(softmax(z_ex) - np.eye(3)[2], e_ce)
print("Example 4.2 (cont.): cross-entropy gradient s - y with y = e_3:", e_ce)

# Example 4.3: grad log det for 2x2, by SymPy on ad - bc, and numerically.
sa, sb, sc, sd = sp.symbols("a b c d", real=True)
Msym = sp.Matrix([[sa, sb], [sc, sd]])
logdet = sp.log(Msym.det())
grad_sym = sp.Matrix(2, 2, lambda i, j: sp.diff(logdet, Msym[i, j]))
assert sp.simplify(grad_sym - Msym.inv().T) == sp.zeros(2, 2)
A2 = np.array([[2., 1.], [1., 3.]])
g_logdet = np.linalg.inv(A2).T
assert np.allclose(g_logdet, [[0.6, -0.2], [-0.2, 0.4]])
assert rel_error(num_grad(lambda W: np.log(np.linalg.det(W)), A2), g_logdet) < 1e-7
print("Example 4.3: grad log det [[2,1],[1,3]] =", g_logdet.tolist(), " (SymPy: d/dA log(ad-bc) == A^{-T})")

# Differential rules (3.5)-(3.7) as first-order approximations with O(||E||^2) error.
E = 1e-4 * rng.standard_normal((n, n))
Ai = np.linalg.inv(Xmat)
assert np.allclose(np.linalg.inv(Xmat + E) - Ai, -Ai @ E @ Ai, atol=1e-6)      # d(A^{-1}) = -A^{-1} dA A^{-1}
assert np.isclose(np.linalg.det(Xmat + E) - np.linalg.det(Xmat),
                  np.linalg.det(Xmat) * np.trace(Ai @ E), rtol=1e-3)            # Jacobi: d det = det tr(A^{-1} dA)
# Hessian of x^T A x / 2 is (A + A^T)/2, eq. (3.22): Jacobian of the gradient.
H_num = num_jac(lambda v: 0.5 * (A + A.T) @ v, x)
assert np.allclose(H_num, (A + A.T) / 2)
# The vec trick (§3.10): F(X) = A X B has Jacobian B^T kron A w.r.t. column-major vec.
A3 = rng.standard_normal((3, 2)); B3 = rng.standard_normal((4, 2)); X3 = rng.standard_normal((2, 4))
vecF = lambda v: (A3 @ v.reshape(X3.shape, order="F") @ B3).ravel(order="F")
assert rel_error(num_jac(vecF, X3.ravel(order="F")), np.kron(B3.T, A3)) < 1e-7
print("Part 2: all identities, Jacobi's formula, d(A^-1), Hessian, vec trick verified.")

# %% [markdown]
# ## Part 3: `Value` — scalar reverse-mode autodiff
#
# Each `Value` remembers its children and a closure that pushes `out.grad`
# to them by the chain rule (one VJP per node). `backward()` topologically
# sorts the graph and runs the closures from the output backwards: this is
# reverse mode (§3.6), and the order is the "cheap order" of the thesis.
#
# **Build it yourself first:** write `__add__`, `__mul__`, `tanh` and
# `backward` before reading this.

# %%
class Value:
    """A scalar node in a computation graph with reverse-mode gradient."""

    def __init__(self, data, children=(), op=""):
        self.data = float(data); self.grad = 0.0
        self._children = tuple(children); self._backward = lambda: None; self._op = op

    def __repr__(self): return f"Value({self.data:.6g}, grad={self.grad:.6g})"

    @staticmethod
    def _w(v): return v if isinstance(v, Value) else Value(v)

    def __add__(self, o):
        o = Value._w(o); out = Value(self.data + o.data, (self, o), "+")
        def _b(): self.grad += out.grad; o.grad += out.grad
        out._backward = _b; return out

    def __mul__(self, o):
        o = Value._w(o); out = Value(self.data * o.data, (self, o), "*")
        def _b(): self.grad += o.data * out.grad; o.grad += self.data * out.grad
        out._backward = _b; return out

    def __pow__(self, k):
        out = Value(self.data ** k, (self,), f"**{k}")
        def _b(): self.grad += k * self.data ** (k - 1) * out.grad
        out._backward = _b; return out

    def exp(self):
        out = Value(np.exp(self.data), (self,), "exp")
        def _b(): self.grad += out.data * out.grad
        out._backward = _b; return out

    def log(self):
        out = Value(np.log(self.data), (self,), "log")
        def _b(): self.grad += out.grad / self.data
        out._backward = _b; return out

    def tanh(self):
        out = Value(np.tanh(self.data), (self,), "tanh")
        def _b(): self.grad += (1 - out.data ** 2) * out.grad
        out._backward = _b; return out

    def relu(self):
        out = Value(max(0.0, self.data), (self,), "relu")
        def _b(): self.grad += (self.data > 0) * out.grad
        out._backward = _b; return out

    __radd__ = __add__; __rmul__ = __mul__
    def __neg__(self): return self * -1
    def __sub__(self, o): return self + (-Value._w(o))
    def __rsub__(self, o): return Value._w(o) + (-self)
    def __truediv__(self, o): return self * Value._w(o) ** -1
    def __rtruediv__(self, o): return Value._w(o) * self ** -1

    def backward(self):
        """Reverse-mode sweep: seed dL/dL = 1, then one VJP per node in reverse topological order."""
        topo, seen = [], set()
        def build(v):
            if id(v) not in seen:
                seen.add(id(v))
                for c in v._children: build(c)
                topo.append(v)
        build(self)
        self.grad = 1.0
        for v in reversed(topo): v._backward()


def value_fn_check(expr, point):
    """Check d expr / d x_i from Value against central differences at `point`."""
    xs = [Value(p) for p in point]; out = expr(*xs); out.backward()
    g_num = num_grad(lambda p: expr(*[Value(q) for q in p]).data, np.array(point, float))
    return rel_error([v.grad for v in xs], g_num)


tests = [
    (lambda u, v: (u * v + u ** 3).tanh() * v.exp() - (u * u + 1).log() / v, [0.7, -1.3]),
    (lambda u, v, w: ((u - v) ** 2 + (v * w).relu() + (w / u).exp()).tanh(), [1.2, 0.4, 0.9]),
    (lambda u: u.exp().log() - u, [2.0]),  # identically zero: gradient should be 0
]
for expr, pt in tests:
    assert value_fn_check(expr, pt) < 1e-7
# A gradient that reuses a variable (the "+=" in _backward matters): f = x*x + x*y + x.
xv, yv = Value(3.0), Value(-2.0); fv = xv * xv + xv * yv + xv; fv.backward()
assert np.isclose(xv.grad, 2 * 3 - 2 + 1) and np.isclose(yv.grad, 3)
print("Part 3: Value autodiff agrees with finite differences; fan-out accumulates correctly.")

# %% [markdown]
# ## Part 4: `Tensor` — matrices, matmul, broadcasting
#
# The VJP of each op is the identification rule applied once (§3.6, eq. 3.24):
# for C = AB and dL = tr(G^T dC), grad_A = G B^T and grad_B = A^T G.
# Broadcasting a bias row b^T onto every row of a matrix means dZ contains
# db^T in every row, so the adjoint sums over the rows (eq. 3.26).

# %%
class Tensor:
    """An ndarray node with reverse-mode gradient. 2-D for matmul; broadcasting for + and *."""

    def __init__(self, data, children=(), op=""):
        self.data = np.array(data, dtype=float); self.grad = np.zeros_like(self.data)
        self._children = tuple(children); self._backward = lambda: None; self._op = op

    @property
    def shape(self): return self.data.shape

    @staticmethod
    def _w(v): return v if isinstance(v, Tensor) else Tensor(v)

    @staticmethod
    def _unbroadcast(g, shape):
        """Sum g over the axes broadcasting created, giving an array of `shape` (eq. 3.26)."""
        g = np.asarray(g, float)
        while g.ndim > len(shape): g = g.sum(axis=0)
        for ax, s in enumerate(shape):
            if s == 1 and g.shape[ax] != 1: g = g.sum(axis=ax, keepdims=True)
        return g.reshape(shape)

    def __add__(self, o):
        o = Tensor._w(o); out = Tensor(self.data + o.data, (self, o), "+")
        def _b():
            self.grad += Tensor._unbroadcast(out.grad, self.shape)
            o.grad += Tensor._unbroadcast(out.grad, o.shape)
        out._backward = _b; return out

    def __mul__(self, o):
        o = Tensor._w(o); out = Tensor(self.data * o.data, (self, o), "*")
        def _b():
            self.grad += Tensor._unbroadcast(out.grad * o.data, self.shape)
            o.grad += Tensor._unbroadcast(out.grad * self.data, o.shape)
        out._backward = _b; return out

    def __matmul__(self, o):
        assert self.data.ndim == 2 and o.data.ndim == 2, "matmul needs 2-D operands"
        out = Tensor(self.data @ o.data, (self, o), "@")
        def _b():
            self.grad += out.grad @ o.data.T      # grad_A = G B^T
            o.grad += self.data.T @ out.grad      # grad_B = A^T G
        out._backward = _b; return out

    def T(self):
        out = Tensor(self.data.T, (self,), "T")
        def _b(): self.grad += out.grad.T
        out._backward = _b; return out

    def sum(self):
        out = Tensor(self.data.sum(), (self,), "sum")
        def _b(): self.grad += out.grad * np.ones_like(self.data)
        out._backward = _b; return out

    def __pow__(self, k):
        out = Tensor(self.data ** k, (self,), f"**{k}")
        def _b(): self.grad += k * self.data ** (k - 1) * out.grad
        out._backward = _b; return out

    def _unary(self, fval, fder, op):
        out = Tensor(fval(self.data), (self,), op)
        def _b(): self.grad += fder(self.data, out.data) * out.grad   # Jacobian of an elementwise map is diagonal
        out._backward = _b; return out

    def tanh(self): return self._unary(np.tanh, lambda x, y: 1 - y ** 2, "tanh")
    def relu(self): return self._unary(lambda x: np.maximum(x, 0), lambda x, y: (x > 0).astype(float), "relu")
    def exp(self): return self._unary(np.exp, lambda x, y: y, "exp")
    def log(self): return self._unary(np.log, lambda x, y: 1 / x, "log")

    __radd__ = __add__; __rmul__ = __mul__
    def __neg__(self): return self * -1.0
    def __sub__(self, o): return self + (-Tensor._w(o))
    def __rsub__(self, o): return Tensor._w(o) + (-self)
    def __truediv__(self, k): return self * (1.0 / k)   # scalar divisor only

    def backward(self):
        topo, seen = [], set()
        def build(v):
            if id(v) not in seen:
                seen.add(id(v))
                for c in v._children: build(c)
                topo.append(v)
        build(self)
        self.grad = np.ones_like(self.data)
        for v in reversed(topo): v._backward()


def tensor_grads(loss_fn, arrays):
    """Wrap arrays as Tensors, run loss_fn, backward, return the gradients."""
    ts = [Tensor(a) for a in arrays]; L = loss_fn(*ts); L.backward()
    return float(L.data), [t.grad for t in ts]


# Random checks of every Tensor op against finite differences.
Xr, Wr, br = rng.standard_normal((6, 3)), rng.standard_normal((3, 4)), rng.standard_normal((1, 4))
loss_t = lambda X, W, bb: (((X @ W + bb).tanh() ** 2).exp().sum() * 0.1 + ((X * X).relu() + 1).log().sum() + (X @ W).T().sum())
Lval, (gX, gW, gb) = tensor_grads(loss_t, [Xr, Wr, br])
loss_np = lambda X, W, bb: float(loss_t(Tensor(X), Tensor(W), Tensor(bb)).data)
assert rel_error(gX, num_grad(lambda v: loss_np(v, Wr, br), Xr)) < 1e-6
assert rel_error(gW, num_grad(lambda v: loss_np(Xr, v, br), Wr)) < 1e-6
assert rel_error(gb, num_grad(lambda v: loss_np(Xr, Wr, v), br)) < 1e-6
print("Part 4: Tensor engine (matmul, broadcast bias, tanh, relu, exp, log, pow, sum, T) agrees with finite differences.")

# %% [markdown]
# ## Part 5: Example 4.4 — backprop by hand, three ways
#
# Network (§3.7): z1 = W1 x + b1, a1 = tanh(z1), z2 = W2 a1 + b2, L = ½||z2 - y||².
# Hand-derived (eqs. 3.28–3.33): δ2 = z2 - y, grad_W2 = δ2 a1^T, grad_b2 = δ2,
# δ1 = (W2^T δ2) ⊙ σ'(z1), grad_W1 = δ1 x^T, grad_b1 = δ1.  Shapes in comments.


# %%
def forward_single(W1, b1, W2, b2, x, act=np.tanh):
    z1 = W1 @ x + b1          # (h,1)
    a1 = act(z1)              # (h,1)
    z2 = W2 @ a1 + b2         # (m,1)
    return z1, a1, z2


def backprop_single(W1, b1, W2, b2, x, y):
    """Hand-derived gradients of L = ½||z2 - y||² for the two-layer tanh net, eqs. (3.28)-(3.33)."""
    z1, a1, z2 = forward_single(W1, b1, W2, b2, x)
    d2 = z2 - y                           # (m,1)   eq. (3.28)
    gW2 = d2 @ a1.T                       # (m,h)   eq. (3.29)
    gb2 = d2                              # (m,1)   eq. (3.30)
    d1 = (W2.T @ d2) * (1 - a1 ** 2)      # (h,1)   eq. (3.31): the transpose is the universal first bug
    gW1 = d1 @ x.T                        # (h,d)   eq. (3.32)
    gb1 = d1                              # (h,1)   eq. (3.33)
    L = 0.5 * float(np.sum((z2 - y) ** 2))
    return L, dict(W1=gW1, b1=gb1, W2=gW2, b2=gb2), dict(z1=z1, a1=a1, z2=z2, d2=d2, d1=d1)


# Example 4.4 numbers.
x_e = np.array([[1.], [2.]]); y_e = np.array([[1.]])
W1_e = np.array([[1., 1.], [2., -1.]]); b1_e = np.array([[-1.], [1.]])
W2_e = np.array([[2., -1.]]); b2_e = np.array([[1.]])
L_e, g_hand, mid = backprop_single(W1_e, b1_e, W2_e, b2_e, x_e, y_e)
print("Example 4.4 intermediates:")
for k, v in mid.items(): print(f"  {k} = {np.round(v.ravel(), 6)}")
print(f"  L = {L_e:.6f}")
for k, v in g_hand.items(): print(f"  grad_{k} =\n{np.round(v, 6)}")
# The numbers quoted in the notes (4 decimals):
assert np.allclose(mid["z1"].ravel(), [2, 1]) and np.allclose(mid["a1"].ravel(), [0.9640, 0.7616], atol=5e-5)
assert np.isclose(mid["z2"].item(), 2.1665, atol=5e-5) and np.isclose(L_e, 0.6803, atol=5e-5)
assert np.isclose(mid["d2"].item(), 1.1665, atol=5e-5)
assert np.allclose(g_hand["W2"], [[1.1245, 0.8884]], atol=5e-5) and np.isclose(g_hand["b2"].item(), 1.1665, atol=5e-5)
assert np.allclose((W2_e.T @ mid["d2"]).ravel(), [2.3329, -1.1665], atol=5e-5)
assert np.allclose((1 - mid["a1"] ** 2).ravel(), [0.0707, 0.4200], atol=5e-5)
assert np.allclose(mid["d1"].ravel(), [0.1648, -0.4899], atol=5e-5)
assert np.allclose(g_hand["W1"], [[0.1648, 0.3296], [-0.4899, -0.9798]], atol=5e-5)
assert np.allclose(g_hand["b1"].ravel(), [0.1648, -0.4899], atol=5e-5)

# (ii) Tensor engine on the same graph.
net_loss = lambda W1, b1, W2, b2: (((W2 @ (W1 @ Tensor(x_e) + b1).tanh() + b2) - Tensor(y_e)) ** 2).sum() * 0.5
L_t, g_t = tensor_grads(net_loss, [W1_e, b1_e, W2_e, b2_e])
assert np.isclose(L_t, L_e)
for gh, gt in zip([g_hand[k] for k in ("W1", "b1", "W2", "b2")], g_t):
    assert np.allclose(gh, gt, atol=1e-12), "hand-derived backprop disagrees with the Tensor engine"

# (iii) Value engine, one scalar node at a time (the network written without any matrices).
def net_values(W1, b1, W2, b2):
    xs = [Value(float(v)) for v in x_e.ravel()]
    z1 = [sum(W1[i][j] * xs[j] for j in range(2)) + b1[i] for i in range(2)]
    a1 = [z.tanh() for z in z1]
    z2 = sum(W2[0][j] * a1[j] for j in range(2)) + b2[0]
    return (z2 - float(y_e.item())) ** 2 * 0.5
W1v = [[Value(v) for v in row] for row in W1_e]; b1v = [Value(v) for v in b1_e.ravel()]
W2v = [[Value(v) for v in row] for row in W2_e]; b2v = [Value(v) for v in b2_e.ravel()]
Lv = net_values(W1v, b1v, W2v, b2v); Lv.backward()
assert np.isclose(Lv.data, L_e)
assert np.allclose([[v.grad for v in row] for row in W1v], g_hand["W1"])
assert np.allclose([v.grad for v in b1v], g_hand["b1"].ravel())
assert np.allclose([[v.grad for v in row] for row in W2v], g_hand["W2"])
assert np.allclose([v.grad for v in b2v], g_hand["b2"].ravel())

# (iv) Finite differences on every parameter.
def L_of(W1, b1, W2, b2): return backprop_single(W1, b1, W2, b2, x_e, y_e)[0]
assert rel_error(g_hand["W1"], num_grad(lambda v: L_of(v, b1_e, W2_e, b2_e), W1_e)) < 1e-7
assert rel_error(g_hand["b1"], num_grad(lambda v: L_of(W1_e, v, W2_e, b2_e), b1_e)) < 1e-7
assert rel_error(g_hand["W2"], num_grad(lambda v: L_of(W1_e, b1_e, v, b2_e), W2_e)) < 1e-7
assert rel_error(g_hand["b2"], num_grad(lambda v: L_of(W1_e, b1_e, W2_e, v), b2_e)) < 1e-7
print("Example 4.4: hand backprop == Tensor engine == Value engine == finite differences.")

# %% [markdown]
# ## Part 6: Batch backprop and training on y = sin(x1) + x2²
#
# Batch version (§3.7, eqs. 3.34–3.37) with X as n×d and L = (1/2n)||Z2 - Y||_F²:
# Z1 = X W1^T + 1 b1^T, A1 = σ(Z1), Z2 = A1 W2^T + 1 b2^T,
# Δ2 = (Z2 - Y)/n, grad_W2 = Δ2^T A1, grad_b2 = Δ2^T 1,
# Δ1 = (Δ2 W2) ⊙ σ'(Z1), grad_W1 = Δ1^T X, grad_b1 = Δ1^T 1.


# %%
def forward_batch(p, X):
    Z1 = X @ p["W1"].T + p["b1"].T        # (n,h)
    A1 = np.tanh(Z1)                      # (n,h)
    Z2 = A1 @ p["W2"].T + p["b2"].T       # (n,m)
    return Z1, A1, Z2


def backprop_batch(p, X, Y):
    """Hand-derived batch gradients of L = (1/2n)||Z2 - Y||_F², eqs. (3.34)-(3.37)."""
    n = X.shape[0]
    Z1, A1, Z2 = forward_batch(p, X)
    D2 = (Z2 - Y) / n                         # (n,m)
    gW2 = D2.T @ A1                           # (m,h)
    gb2 = D2.T @ np.ones((n, 1))              # (m,1): column sums
    D1 = (D2 @ p["W2"]) * (1 - A1 ** 2)       # (n,h)
    gW1 = D1.T @ X                            # (h,d)
    gb1 = D1.T @ np.ones((n, 1))              # (h,1)
    L = 0.5 * np.sum((Z2 - Y) ** 2) / n
    return L, dict(W1=gW1, b1=gb1, W2=gW2, b2=gb2)


def init_params(d, h, m, scale=0.5):
    return dict(W1=scale * rng.standard_normal((h, d)), b1=np.zeros((h, 1)),
                W2=scale * rng.standard_normal((m, h)), b2=np.zeros((m, 1)))


# Data: y = sin(x1) + x2^2 on [-2, 2]^2.
N = 256
Xd = rng.uniform(-2, 2, size=(N, 2)); Yd = (np.sin(Xd[:, 0]) + Xd[:, 1] ** 2)[:, None]
params = init_params(2, 32, 1)

# Check the batch formulas against the Tensor engine once, at initialization.
def batch_loss_t(W1, b1, W2, b2):
    Z2 = (Tensor(Xd) @ W1.T() + b1.T()).tanh() @ W2.T() + b2.T()
    return ((Z2 - Tensor(Yd)) ** 2).sum() * (0.5 / N)
L0, g0 = backprop_batch(params, Xd, Yd)
L0_t, g0_t = tensor_grads(batch_loss_t, [params[k] for k in ("W1", "b1", "W2", "b2")])
assert np.isclose(L0, L0_t)
for k, gt in zip(("W1", "b1", "W2", "b2"), g0_t): assert np.allclose(g0[k], gt, atol=1e-12), k
assert rel_error(g0["W1"], num_grad(lambda v: backprop_batch({**params, "W1": v}, Xd, Yd)[0], params["W1"])) < 1e-6

# Train with plain gradient descent on the hand-derived gradients.
lr, steps = 0.05, 4000
losses = []
for it in range(steps):
    L, g = backprop_batch(params, Xd, Yd)
    losses.append(L)
    for k in params: params[k] -= lr * g[k]
L_final = backprop_batch(params, Xd, Yd)[0]
assert L_final < 0.02 * losses[0] and L_final < 0.05, (losses[0], L_final)
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].semilogy(losses); axes[0].set_xlabel("GD step"); axes[0].set_ylabel("½ mean squared error"); axes[0].set_title("Training loss (hand-derived batch backprop)"); axes[0].grid(alpha=.3)
gx, gy = np.meshgrid(np.linspace(-2, 2, 60), np.linspace(-2, 2, 60))
Gpts = np.c_[gx.ravel(), gy.ravel()]
pred = forward_batch(params, Gpts)[2].reshape(gx.shape); truth = np.sin(gx) + gy ** 2
cs = axes[1].contourf(gx, gy, pred, 20); axes[1].contour(gx, gy, truth, 10, colors="w", linewidths=.6)
axes[1].set_title("Net output (filled) vs sin(x1)+x2² (white)"); fig.colorbar(cs, ax=axes[1])
fig.tight_layout(); fig.savefig(FIG / "training_sin_x2sq.png", dpi=120); plt.close(fig)
print(f"Part 6: batch formulas match the engine; training loss {losses[0]:.3f} -> {L_final:.4f} after {steps} GD steps.")

# %% [markdown]
# ## Part 7: Example 4.1 — Newton vs gradient descent on a quadratic
#
# f(x) = ½ x^T A x - b^T x with A = [[2,1],[1,3]], b = (1,1). grad f = Ax - b,
# H = A. Newton from 0: x1 = A^{-1} b = (0.4, 0.2), the exact minimizer, f* = -0.3.

# %%
Aq = np.array([[2., 1.], [1., 3.]]); bq = np.array([1., 1.])
fq = lambda v: 0.5 * v @ Aq @ v - bq @ v
gq = lambda v: Aq @ v - bq
grad_check(fq, gq, rng.standard_normal(2))
assert np.allclose(num_jac(gq, np.zeros(2)), Aq)                           # Hessian = A
x_star = np.linalg.solve(Aq, bq)
x_newton = np.zeros(2) - np.linalg.solve(Aq, gq(np.zeros(2)))               # one Newton step, eq. (3.23)
assert np.allclose(x_newton, [0.4, 0.2]) and np.allclose(x_newton, x_star) and np.isclose(fq(x_star), -0.3)
assert np.allclose(gq(np.zeros(2)), [-1, -1])
lam = np.linalg.eigvalsh(Aq)
assert np.allclose(lam, [(5 - np.sqrt(5)) / 2, (5 + np.sqrt(5)) / 2])


def gradient_descent(grad, x0, eta, tol=1e-8, x_true=None, maxit=10000):
    xs = [np.array(x0, float)]
    while len(xs) < maxit:
        xs.append(xs[-1] - eta * grad(xs[-1]))
        if np.linalg.norm(xs[-1] - x_true) < tol: break
    return np.array(xs)


path_gd = gradient_descent(gq, np.zeros(2), 1 / lam[1], x_true=x_star)
path_gd_opt = gradient_descent(gq, np.zeros(2), 2 / (lam[0] + lam[1]), x_true=x_star)
n_gd, n_gd_opt = len(path_gd) - 1, len(path_gd_opt) - 1
assert n_gd == 36 and n_gd_opt == 22, (n_gd, n_gd_opt)
# Linear convergence rate with the optimal step is (kappa-1)/(kappa+1) = (lmax-lmin)/(lmax+lmin).
rate = (lam[1] - lam[0]) / (lam[1] + lam[0])
errs = np.linalg.norm(path_gd_opt - x_star, axis=1)
assert np.allclose(errs[1:12] / errs[0:11], rate, atol=0.02)
fig, ax = plt.subplots(figsize=(5.5, 5))
g1, g2 = np.meshgrid(np.linspace(-0.3, 0.9, 200), np.linspace(-0.4, 0.7, 200))
Fq = 0.5 * (Aq[0, 0] * g1 ** 2 + 2 * Aq[0, 1] * g1 * g2 + Aq[1, 1] * g2 ** 2) - bq[0] * g1 - bq[1] * g2
ax.contour(g1, g2, Fq, levels=25, linewidths=.6)
ax.plot(path_gd[:, 0], path_gd[:, 1], "o-", ms=3, label=f"GD, step 1/λmax ({n_gd} its)")
ax.plot([0, x_newton[0]], [0, x_newton[1]], "r*-", ms=12, label="Newton (1 step)")
ax.set_aspect("equal"); ax.legend(fontsize=8); ax.set_title("Example 4.1: Newton vs gradient descent"); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "newton_vs_gd.png", dpi=120); plt.close(fig)
print(f"Example 4.1: Newton lands on {x_newton} in 1 step; GD needs {n_gd} (step 1/λmax) or {n_gd_opt} (optimal step) iterations to 1e-8; rate {rate:.4f}.")

# %% [markdown]
# ## Part 8: Example 4.5 — Gauss–Newton on y = a e^{-bt}
#
# r_i(a,b) = a e^{-b t_i} - y_i, J_r = [e^{-b t_i}, -a t_i e^{-b t_i}],
# grad f = J_r^T r, H ≈ J_r^T J_r, step h = -(J_r^T J_r)^{-1} J_r^T r (eq. 3.25).

# %%
t_d = np.array([0., 1., 2.]); y_d = np.array([2.0, 1.2, 0.8])
r_fn = lambda p: p[0] * np.exp(-p[1] * t_d) - y_d
J_fn = lambda p: np.column_stack([np.exp(-p[1] * t_d), -p[0] * t_d * np.exp(-p[1] * t_d)])
assert rel_error(num_jac(r_fn, np.array([1.5, 0.5])), J_fn(np.array([1.5, 0.5]))) < 1e-7
f_ls = lambda p: 0.5 * r_fn(p) @ r_fn(p)
grad_check(f_ls, lambda p: J_fn(p).T @ r_fn(p), np.array([1.5, 0.5]))   # grad f = J^T r


def gauss_newton(r, J, x0, tol=1e-12, maxit=50):
    """Undamped Gauss–Newton: solve (J^T J) h = -J^T r, x <- x + h. Returns the iterate trace."""
    xs = [np.array(x0, float)]
    for _ in range(maxit):
        Jk, rk = J(xs[-1]), r(xs[-1])
        h = -np.linalg.solve(Jk.T @ Jk, Jk.T @ rk)
        xs.append(xs[-1] + h)
        if np.linalg.norm(h) < tol: break
    return np.array(xs)


p0 = np.array([1.5, 0.5])
trace = gauss_newton(r_fn, J_fn, p0)
J0, r0 = J_fn(p0), r_fn(p0)
print("Example 4.5: r(x0) =", np.round(r0, 4), "\n  J_r(x0) =\n", np.round(J0, 4))
print("  J^T J =\n", np.round(J0.T @ J0, 4), "\n  J^T r =", np.round(J0.T @ r0, 4), "\n  step =", np.round(trace[1] - trace[0], 4))
assert np.allclose(r0, [-0.5, -0.2902, -0.2482], atol=5e-5)
assert np.allclose(J0, [[1, 0], [0.6065, -0.9098], [0.3679, -1.1036]], atol=5e-5)
assert np.allclose(J0.T @ J0, [[1.5032, -0.9578], [-0.9578, 2.0457]], atol=5e-5)
assert np.allclose(J0.T @ r0, [-0.7673, 0.5379], atol=5e-5)
assert np.allclose(trace[1], [1.9887, 0.4659], atol=5e-5)
assert np.allclose(trace[-1], [1.9884, 0.4743], atol=5e-5) and np.isclose(f_ls(trace[-1]), 1.2153e-3, rtol=1e-3)
sol = least_squares(r_fn, p0, jac=J_fn, xtol=1e-14, ftol=1e-14, gtol=1e-14)
assert np.allclose(sol.x, trace[-1], atol=1e-8)
steps_gn = np.linalg.norm(np.diff(trace, axis=0), axis=1)
n_gn = int(np.argmax(steps_gn < 1e-6)) + 1
# The dropped term sum_i r_i Hess(r_i) at the solution is small compared with J^T J.
a_s, b_s = trace[-1]
Hr = [np.array([[0, -ti * np.exp(-b_s * ti)], [-ti * np.exp(-b_s * ti), a_s * ti ** 2 * np.exp(-b_s * ti)]]) for ti in t_d]
dropped = sum(ri * Hi for ri, Hi in zip(r_fn(trace[-1]), Hr))
JtJ_s = J_fn(trace[-1]).T @ J_fn(trace[-1])
assert np.linalg.norm(dropped) < 0.02 * np.linalg.norm(JtJ_s)
assert np.allclose(num_jac(lambda p: J_fn(p).T @ r_fn(p), trace[-1]), JtJ_s + dropped, atol=1e-6)  # exact Hessian
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
tt = np.linspace(0, 2.5, 100)
axes[0].plot(t_d, y_d, "ko", label="data"); axes[0].plot(tt, a_s * np.exp(-b_s * tt), label=f"fit a={a_s:.4f}, b={b_s:.4f}")
axes[0].plot(tt, p0[0] * np.exp(-p0[1] * tt), "--", label="start (1.5, 0.5)"); axes[0].legend(); axes[0].set_xlabel("t"); axes[0].grid(alpha=.3)
axes[1].semilogy(steps_gn, "o-"); axes[1].set_xlabel("Gauss–Newton iteration"); axes[1].set_ylabel("||step||"); axes[1].set_title("Quadratic-looking convergence"); axes[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "gauss_newton_fit.png", dpi=120); plt.close(fig)
print(f"  converged to (a,b) = {np.round(trace[-1], 6)} in {n_gn} steps (|step|<1e-6); f* = {f_ls(trace[-1]):.6e}; scipy agrees.")
print(f"  dropped Hessian term at solution: {np.round(dropped, 4).tolist()} vs J^T J {np.round(JtJ_s, 4).tolist()}")

# %% [markdown]
# ## Part 9: Gaussian log-likelihood gradients and the MLE by gradient ascent (§3.9)
#
# l(μ, Σ) = -½ log det Σ - ½ (x-μ)^T Σ^{-1} (x-μ) + const.
# grad_μ l = Σ^{-1}(x - μ),  grad_Σ l = ½(Σ^{-1}(x-μ)(x-μ)^T Σ^{-1} - Σ^{-1})   (eqs. 3.38, 3.39).


# %%
def gauss_loglik(x, mu, Sig):
    u = x - mu; p = len(x)
    return -0.5 * p * np.log(2 * np.pi) - 0.5 * np.log(np.linalg.det(Sig)) - 0.5 * u @ np.linalg.solve(Sig, u)


def grad_mu(x, mu, Sig): return np.linalg.solve(Sig, x - mu)


def grad_Sig(x, mu, Sig):
    Si = np.linalg.inv(Sig); u = x - mu
    return 0.5 * (Si @ np.outer(u, u) @ Si - Si)


xg = rng.standard_normal(3); mug = rng.standard_normal(3)
Lg = rng.standard_normal((3, 3)); Sg = Lg @ Lg.T + np.eye(3)
grad_check(lambda mm: gauss_loglik(xg, mm, Sg), lambda mm: grad_mu(xg, mm, Sg), mug)
# grad_Sig is the unconstrained gradient (all p^2 entries free), so compare with FD over all entries.
grad_check(lambda SS: gauss_loglik(xg, mug, SS), lambda SS: grad_Sig(xg, mug, SS), Sg)

# Sample, then recover xbar and S/n by gradient ascent on the average log-likelihood.
mu_true = np.array([1., -2.]); Sig_true = np.array([[2., .6], [.6, 1.]])
Xs = rng.multivariate_normal(mu_true, Sig_true, size=400)
xbar = Xs.mean(0); S_hat = (Xs - xbar).T @ (Xs - xbar) / len(Xs)
mu_k, Sig_k, eta = np.zeros(2), np.eye(2), 0.5
for k in range(5000):
    U = Xs - mu_k; Si = np.linalg.inv(Sig_k)
    g_mu = Si @ U.mean(0)                              # average of grad_mu over the sample
    g_S = 0.5 * (Si @ (U.T @ U / len(Xs)) @ Si - Si)   # average of grad_Sig
    mu_k = mu_k + eta * g_mu; Sig_k = Sig_k + eta * g_S
    assert np.all(np.linalg.eigvalsh(Sig_k) > 0), "left the PD cone"
    if np.linalg.norm(g_mu) < 1e-10 and np.linalg.norm(g_S) < 1e-10: break
assert np.allclose(mu_k, xbar, atol=1e-8) and np.allclose(Sig_k, S_hat, atol=1e-8)
print(f"Part 9: Gaussian gradients verified; gradient ascent ({k} steps) recovers mu_hat = {np.round(mu_k, 4)}, Sigma_hat = {np.round(Sig_k, 4).tolist()}")

# %% [markdown]
# ## Part 10: Forward vs reverse mode cost (§3.6)
#
# For a chain of k Jacobians of size d×d ending in a 1×d row, multiplying
# from the left (reverse mode: vector–Jacobian products) costs ~k d^2 flops;
# forming the matrix products (or k JVPs for all d directions) costs ~k d^3.

# %%
d_, k_ = 400, 8
Js = [rng.standard_normal((d_, d_)) / np.sqrt(d_) for _ in range(k_)]
row = rng.standard_normal((1, d_))
t0 = time.perf_counter()
Pm = np.eye(d_)
for Jk in Js[::-1]: Pm = Pm @ Jk            # the full d×d product, right to left (never do this)
grad_full = (row @ Pm).ravel()
t_full = time.perf_counter() - t0
t0 = time.perf_counter()
v = row
for Jk in Js[::-1]: v = v @ Jk               # reverse mode: a row times a Jacobian, outermost first
t_rev = time.perf_counter() - t0
assert np.allclose(v.ravel(), grad_full)
print(f"Part 10: d={d_}, k={k_}: full Jacobian product {t_full*1e3:.1f} ms vs reverse-mode VJP chain {t_rev*1e3:.2f} ms "
      f"(wall-clock ratio ~{t_full/t_rev:.0f}; flop ratio is d = {d_}, BLAS matmul is more efficient per flop so the clock ratio is smaller).")

# Grid deformation figure for §2: f(x,y) = (x² - y², 2xy), Jacobian at (1, 0.5).
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
F2 = lambda P: np.c_[P[:, 0] ** 2 - P[:, 1] ** 2, 2 * P[:, 0] * P[:, 1]]
for c in np.linspace(0.2, 1.6, 8):
    s_ = np.linspace(0.2, 1.6, 100)
    for line in (np.c_[s_, np.full_like(s_, c)], np.c_[np.full_like(s_, c), s_]):
        axes[0].plot(line[:, 0], line[:, 1], color="0.7", lw=.7); img = F2(line); axes[1].plot(img[:, 0], img[:, 1], color="0.7", lw=.7)
p_ = np.array([1., .5]); Jp = np.array([[2 * p_[0], -2 * p_[1]], [2 * p_[1], 2 * p_[0]]]); r_ = 0.15
sq = p_ + r_ * np.array([[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, -1]])
axes[0].plot(sq[:, 0], sq[:, 1], "C3", lw=2); axes[1].plot(*F2(sq).T, "C3", lw=2, label="image of the square")
lin = F2(p_[None])[0] + (sq - p_) @ Jp.T
axes[1].plot(lin[:, 0], lin[:, 1], "C0--", lw=1.5, label="f(p) + J (q - p)")
axes[0].set_title("domain: grid and a small square at p=(1, 0.5)"); axes[1].set_title("image under f(x,y) = (x²-y², 2xy)"); axes[1].legend(fontsize=8)
for ax_ in axes: ax_.set_aspect("equal"); ax_.grid(alpha=.2)
fig.tight_layout(); fig.savefig(FIG / "jacobian_grid.png", dpi=120); plt.close(fig)
assert rel_error(num_jac(lambda q: F2(q[None])[0], p_), Jp) < 1e-7

# %% [markdown]
# ## Your turn
#
# 1. Add a ReLU option and a second hidden layer to `forward_batch` /
#    `backprop_batch` (so z1 → a1 → z2 → a2 → z3). Derive δ2 = (W3^T δ3) ⊙ σ'(z2)
#    first, then re-verify against the Tensor engine and `num_grad` exactly as
#    Part 6 does. `Tensor.relu` already exists.
# 2. Derive grad_A ||A X - B||_F^2 = 2 (A X - B) X^T with differentials
#    (exercise B3) and add it to the `identities` dict in Part 2.
# 3. Forward mode: write `jvp(F, x, v)` returning J_F(x) v by central
#    differences, then compute the full gradient of f: R^1000 → R,
#    f(x) = sum(tanh(W x)) with W 1000×1000, two ways: 1000 JVPs vs one
#    reverse pass through `Tensor`. Time both. Which is the ~d factor of §3.6?
# 4. Replace `tanh` with the sigmoid in Part 5 and redo Example 4.4 by hand;
#    only σ'(z1) changes. Check it with the engine.

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
