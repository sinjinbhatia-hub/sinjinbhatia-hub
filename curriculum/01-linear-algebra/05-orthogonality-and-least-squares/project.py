# %% [markdown]
# # LA-5 Project: Projections, Least Squares, Gram–Schmidt, and QR
#
# This script builds every object in LA-5 from scratch and checks each
# numerical claim in the notes. Parts: (1) the four subspaces meet at right
# angles, (2) projection onto a subspace (Example 4.1), (3) least squares by
# the normal equations and by calculus (Example 4.2), (4) polynomial fitting:
# a sensor calibration (Example 4.5) and noisy synthetic data with a residual
# plot, (5) classical vs modified Gram–Schmidt (Example 4.3) and the loss of
# orthogonality on Hilbert matrices, (6) QR and least squares via R x = Q^T b
# (Example 4.4), with an ill-conditioned fit that shows why QR beats the normal
# equations, (7) orthogonal matrices and Householder reflections, (8) orthogonal
# functions: the start of Fourier series.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import scipy.linalg as sla
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)

# %% [markdown]
# ## Part 1: The four subspaces, now with right angles
#
# Theorem 3.9, eq. (3.3): N(A) = C(A^T)^⊥ and N(A^T) = C(A)^⊥. We build a 6×4
# matrix of rank 2, get orthonormal bases of all four subspaces from SciPy, and
# check (i) the dimensions add up (Rank–Nullity, LA-4), (ii) every basis vector
# of one subspace is orthogonal to every basis vector of its complement, and
# (iii) the Fredholm alternative: Ax = b is solvable iff b ⟂ N(A^T).


# %%
def four_subspaces(A):
    """Orthonormal bases (as columns) of C(A), N(A^T), C(A^T), N(A), in that order."""
    A = np.asarray(A, float)
    return sla.orth(A), sla.null_space(A.T), sla.orth(A.T), sla.null_space(A)


m, n, r = 6, 4, 2
while True:  # random integer matrix of rank exactly r
    A_r = (rng.integers(-3, 4, size=(m, r)) @ rng.integers(-3, 4, size=(r, n))).astype(float)
    if np.linalg.matrix_rank(A_r) == r:
        break
col, leftnull, row, null = four_subspaces(A_r)
assert col.shape[1] == r and row.shape[1] == r
assert null.shape[1] == n - r and leftnull.shape[1] == m - r      # Rank–Nullity, in R^n and in R^m
assert np.allclose(row.T @ null, 0)                                # C(A^T) ⟂ N(A)
assert np.allclose(col.T @ leftnull, 0)                            # C(A)   ⟂ N(A^T)
b_in, b_out = A_r @ rng.standard_normal(n), rng.standard_normal(m)
assert np.allclose(leftnull.T @ b_in, 0) and not np.allclose(leftnull.T @ b_out, 0)   # Fredholm
# Theorem 3.7: every x splits uniquely as (row-space part) + (nullspace part), and A only sees the first.
x = rng.standard_normal(n)
x_row, x_null = row @ (row.T @ x), null @ (null.T @ x)
assert np.allclose(x_row + x_null, x) and np.isclose(x_row @ x_null, 0)
assert np.allclose(A_r @ x, A_r @ x_row)
print(f"Part 1: rank {r}; dim N(A) = {null.shape[1]}, dim N(A^T) = {leftnull.shape[1]}; "
      "all four right angles verified; Fredholm alternative verified.")

# %% [markdown]
# ## Part 2: Projection onto a subspace
#
# (3.4) line: p = (a·b / a·a) a.  (3.6) normal equations A^T A x̂ = A^T b.
# (3.8) P = A (A^T A)^{-1} A^T.  Theorems 3.12–3.14: P² = P, P^T = P, I − P
# projects onto the complement, and p is the closest point of C(A) to b.
#
# **Build it yourself first:** write `project_onto_subspace(A, b)` returning
# the projection p, the coefficients x̂, the error e = b − p, and the matrix P.
# Solve the normal equations with `np.linalg.solve`; do not form the inverse.


# %%
def project_onto_line(b, a):
    """Projection of b onto the line through a: p = (a·b / a·a) a, eq. (3.4). LA-1 Example 4.5."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    return (a @ b) / (a @ a) * a


def project_onto_subspace(A, b):
    """Orthogonal projection of b onto C(A) for A with independent columns.
    Solves the normal equations (3.6) for x̂, then p = A x̂.  Returns (p, x_hat, e, P) with
    e = b − p ⟂ C(A) and P = A (A^T A)^{-1} A^T the projection matrix (3.8)."""
    A = np.asarray(A, float); b = np.asarray(b, float)
    G = A.T @ A                        # Gram matrix; invertible iff columns independent (Theorem 3.11)
    x_hat = np.linalg.solve(G, A.T @ b)
    p = A @ x_hat
    P = A @ np.linalg.solve(G, A.T)    # (A^T A)^{-1} A^T computed by a solve, never an explicit inverse
    return p, x_hat, b - p, P


# LA-1 Example 4.5 recap: project (3, 1) onto the line through (1, 1).
assert np.allclose(project_onto_line([3, 1], [1, 1]), [2, 2])

# Example 4.1: b = (1, 2, 3) onto the plane spanned by (1, 0, 1) and (1, 1, 0).
A1 = np.array([[1, 1], [0, 1], [1, 0]], float)
b1 = np.array([1, 2, 3], float)
p1, xh1, e1, P1 = project_onto_subspace(A1, b1)
assert np.allclose(A1.T @ A1, [[2, 1], [1, 2]]) and np.allclose(A1.T @ b1, [4, 3])
assert np.allclose(xh1, [5 / 3, 2 / 3])
assert np.allclose(p1, [7 / 3, 2 / 3, 5 / 3])
assert np.allclose(e1, [-4 / 3, 4 / 3, 4 / 3])
assert np.allclose(A1.T @ e1, 0)                                   # error ⟂ both columns, eq. (3.5)
assert np.allclose(P1, np.array([[2, 1, 1], [1, 2, -1], [1, -1, 2]]) / 3)
assert np.allclose(P1 @ P1, P1) and np.allclose(P1.T, P1)          # Theorem 3.12
assert np.allclose(P1 @ b1, p1)
assert np.isclose(np.linalg.norm(b1) ** 2, np.linalg.norm(p1) ** 2 + np.linalg.norm(e1) ** 2)  # Pythagoras
I3 = np.eye(3)
n_normal = np.cross(A1[:, 0], A1[:, 1])                            # (-1, 1, 1), normal to the plane
assert np.allclose(I3 - P1, np.outer(n_normal, n_normal) / (n_normal @ n_normal))  # Theorem 3.13
assert np.allclose((I3 - P1) @ b1, e1)
assert np.allclose(P1 + (I3 - P1), I3) and np.allclose(P1 @ (I3 - P1), 0)
# Theorem 3.14: p is closer to b than any other point of the plane.
for _ in range(200):
    v = A1 @ rng.standard_normal(2)
    assert np.linalg.norm(b1 - v) >= np.linalg.norm(b1 - p1) - 1e-12
print("Example 4.1: x̂ =", xh1, " p =", p1, " e =", e1, " |e|² =", round(e1 @ e1, 12), "(= 16/3)")

# Figure: the plane, b, its shadow p, and the perpendicular error.
fig = plt.figure(figsize=(6, 5)); ax = fig.add_subplot(projection="3d")
s, t = np.meshgrid(np.linspace(-0.5, 2.5, 8), np.linspace(-0.5, 1.5, 8))
pts = s[..., None] * A1[:, 0] + t[..., None] * A1[:, 1]
ax.plot_surface(pts[..., 0], pts[..., 1], pts[..., 2], alpha=0.25, color="C0")
ax.quiver(0, 0, 0, *b1, color="k", arrow_length_ratio=0.08, label="b")
ax.quiver(0, 0, 0, *p1, color="C3", arrow_length_ratio=0.1, label="p = Pb")
ax.quiver(*p1, *e1, color="C2", arrow_length_ratio=0.1, label="e = b − p  (perpendicular to plane)")
ax.set_xlim(0, 3); ax.set_ylim(0, 3); ax.set_zlim(0, 3); ax.legend(fontsize=8)
ax.set_title("Example 4.1: projection onto a plane in R³")
fig.tight_layout(); fig.savefig(FIG / "projection_plane.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 3: Least squares by the normal equations, and by calculus
#
# (3.9) f(x) = ‖Ax − b‖².  Geometry (§3.6): the minimizer makes Ax̂ the
# projection of b, so A^T(b − Ax̂) = 0.  Calculus: ∇f = 2A^T(Ax − b), eq. (3.10),
# and the Hessian is 2A^TA.  Both roads give the normal equations; SymPy checks
# the gradient formula symbolically.


# %%
def least_squares_normal(A, b):
    """x̂ = argmin ‖Ax − b‖² via the normal equations (3.6), A^T A x̂ = A^T b.
    Uses a linear solve, never the explicit inverse written in (3.7)."""
    A = np.asarray(A, float); b = np.asarray(b, float)
    return np.linalg.solve(A.T @ A, A.T @ b)


def ls_objective(A, x, b):
    """f(x) = ‖Ax − b‖², eq. (3.9)."""
    return float(np.sum((A @ x - b) ** 2))


def ls_gradient(A, x, b):
    """∇f(x) = 2 A^T (Ax − b), eq. (3.10)."""
    return 2 * A.T @ (A @ x - b)


# Example 4.2: least-squares line through (0, 1), (1, 2), (2, 2), (3, 4).
t2 = np.array([0, 1, 2, 3], float); y2 = np.array([1, 2, 2, 4], float)
A2 = np.c_[np.ones(4), t2]                                          # design matrix, columns 1 and t
xh2 = least_squares_normal(A2, y2)
assert np.allclose(A2.T @ A2, [[4, 6], [6, 14]]) and np.allclose(A2.T @ y2, [9, 18])
assert np.allclose(xh2, [0.9, 0.9])
r2 = y2 - A2 @ xh2
assert np.allclose(r2, [0.1, 0.2, -0.7, 0.4]) and np.isclose(r2 @ r2, 0.7)
assert np.allclose(A2.T @ r2, 0)                                    # residual ⟂ both columns
assert np.allclose(xh2, np.linalg.lstsq(A2, y2, rcond=None)[0])     # against NumPy
assert np.allclose(ls_gradient(A2, xh2, y2), 0)                     # stationary point of f
for _ in range(100):                                                # and a minimum
    assert ls_objective(A2, xh2 + 0.3 * rng.standard_normal(2), y2) > ls_objective(A2, xh2, y2)
print("Example 4.2: C, D =", xh2, " residual =", r2, " |r|² =", round(r2 @ r2, 12))

# SymPy: gradient and Hessian of ‖Ax − b‖² for a symbolic 3×2 system, eq. (3.10).
As = sp.Matrix(3, 2, sp.symbols("a0:6")); xs = sp.Matrix(sp.symbols("x0:2")); bs = sp.Matrix(sp.symbols("b0:3"))
res = As * xs - bs
f_sym = (res.T * res)[0]
grad_sym = sp.Matrix([sp.diff(f_sym, xi) for xi in xs])
assert sp.simplify(grad_sym - 2 * As.T * res) == sp.zeros(2, 1)
assert sp.simplify(sp.hessian(f_sym, xs) - 2 * As.T * As) == sp.zeros(2, 2)
print("SymPy: ∇‖Ax−b‖² = 2Aᵀ(Ax−b) and Hessian = 2AᵀA verified symbolically.")

fig, ax = plt.subplots(figsize=(5.5, 4))
ax.plot(t2, y2, "ko", label="data")
tt = np.linspace(-0.3, 3.3, 50); ax.plot(tt, xh2[0] + xh2[1] * tt, "C0", label=f"fit y = {xh2[0]:.1f} + {xh2[1]:.1f} t")
for ti, yi, pi in zip(t2, y2, A2 @ xh2):
    ax.plot([ti, ti], [pi, yi], "C3", lw=2)
ax.plot([], [], "C3", lw=2, label="residuals")
ax.grid(alpha=.3); ax.legend(); ax.set_xlabel("t"); ax.set_ylabel("y"); ax.set_title("Example 4.2: least-squares line")
fig.tight_layout(); fig.savefig(FIG / "least_squares_line.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 4: Polynomial fits: a sensor calibration, then noisy data
#
# Example 4.5: a pressure transducer read at P = 0, 1, ..., 5 bar gives voltages
# V; fit V = a + bP + cP² with the design matrix of columns 1, P, P².
# Then: 40 noisy samples of a quadratic, fitted and checked with a residual plot
# (the residual is orthogonal to every column, so it has zero mean and zero
# correlation with P — if it shows a pattern, your model is missing a term).


# %%
def design_matrix(t, degree):
    """Design (Vandermonde) matrix with columns 1, t, ..., t^degree (§3.6)."""
    return np.vander(np.asarray(t, float), degree + 1, increasing=True)


P5 = np.arange(6, dtype=float)
V5 = np.array([0.45, 1.39, 2.22, 3.04, 3.95, 5.05])
A5 = design_matrix(P5, 2)
assert np.allclose(A5.T @ A5, [[6, 15, 55], [15, 55, 225], [55, 225, 979]])
assert np.allclose(A5.T @ V5, [16.10, 56.00, 227.08])
c5 = least_squares_normal(A5, V5)
assert np.allclose(c5, [0.5, 0.8, 0.02])
r5 = V5 - A5 @ c5
assert np.allclose(r5, [-0.05, 0.07, 0.04, -0.04, -0.07, 0.05]) and np.isclose(r5 @ r5, 0.018)
assert np.allclose(A5.T @ r5, 0)
print("Example 4.5: a, b, c =", c5, " residuals =", r5, " RMS =", round(np.sqrt(r5 @ r5 / 6), 4), "V")

# Noisy synthetic data: true model 0.5 + 2t − 1.5t², noise σ = 0.08.
t_syn = np.linspace(0, 2, 40)
y_true = 0.5 + 2 * t_syn - 1.5 * t_syn ** 2
y_syn = y_true + 0.08 * rng.standard_normal(t_syn.size)
A_syn = design_matrix(t_syn, 2)
c_syn = least_squares_normal(A_syn, y_syn)
r_syn = y_syn - A_syn @ c_syn
assert np.allclose(A_syn.T @ r_syn, 0, atol=1e-10)                 # residual ⟂ 1, t, t²: mean zero, uncorrelated with t
assert np.allclose(c_syn, np.polyfit(t_syn, y_syn, 2)[::-1])       # against NumPy's polyfit
assert np.linalg.norm(c_syn - [0.5, 2, -1.5]) < 0.15               # recovered the generating coefficients
# A model that is missing a term: fit a line to the same data and look at the residual pattern.
c_lin = least_squares_normal(A_syn[:, :2], y_syn)
r_lin = y_syn - A_syn[:, :2] @ c_lin
assert r_lin @ r_lin > 10 * (r_syn @ r_syn)                        # the pattern costs a lot of squared error
print(f"Noisy quadratic fit: coefficients {np.round(c_syn, 3)}; RMS residual {np.sqrt(r_syn @ r_syn / 40):.3f}"
      f" (quadratic) vs {np.sqrt(r_lin @ r_lin / 40):.3f} (line, wrong model).")

fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), sharex="col")
axes[0, 0].plot(P5, V5, "ko", label="readings"); pp = np.linspace(0, 5, 100)
axes[0, 0].plot(pp, design_matrix(pp, 2) @ c5, "C0", label=f"V = {c5[0]:.2f} + {c5[1]:.2f}P + {c5[2]:.2f}P²")
axes[0, 0].set_ylabel("V (volts)"); axes[0, 0].legend(fontsize=8); axes[0, 0].set_title("Example 4.5: sensor calibration")
axes[1, 0].stem(P5, r5); axes[1, 0].axhline(0, color="k", lw=.5); axes[1, 0].set_ylabel("residual (V)"); axes[1, 0].set_xlabel("P (bar)")
axes[0, 1].plot(t_syn, y_syn, "k.", label="noisy data"); axes[0, 1].plot(t_syn, A_syn @ c_syn, "C0", label="quadratic fit")
axes[0, 1].plot(t_syn, A_syn[:, :2] @ c_lin, "C3--", label="line (wrong model)"); axes[0, 1].legend(fontsize=8)
axes[0, 1].set_title("Noisy synthetic data")
axes[1, 1].plot(t_syn, r_syn, "C0.", label="quadratic residual"); axes[1, 1].plot(t_syn, r_lin, "C3x", label="line residual")
axes[1, 1].axhline(0, color="k", lw=.5); axes[1, 1].legend(fontsize=8); axes[1, 1].set_xlabel("t")
for ax in axes.ravel():
    ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "polyfit_residuals.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 5: Gram–Schmidt, classical and modified
#
# (3.14): v_j = a_j − Σ_{i<j} (q_i·a_j) q_i,  q_j = v_j / ‖v_j‖.  The only
# difference between the two versions is *which* vector the inner product
# q_i·(·) is taken with: the original a_j (classical) or the current, partly
# deflated v (modified). Identical in exact arithmetic (Exercise C4); very
# different in floating point.
#
# **Build it yourself first:** write `gram_schmidt_classical(A)` returning Q, R.


# %%
def gram_schmidt_classical(A):
    """Classical Gram–Schmidt, eq. (3.14): A = QR with Q orthonormal columns, R upper triangular,
    diag(R) > 0.  r_ij = q_i^T a_j uses the ORIGINAL column a_j: numerically unstable (§3.8, §8)."""
    A = np.asarray(A, float); m, n = A.shape
    Q = np.zeros((m, n)); R = np.zeros((n, n))
    for j in range(n):
        v = A[:, j].copy()
        for i in range(j):
            R[i, j] = Q[:, i] @ A[:, j]
            v -= R[i, j] * Q[:, i]
        R[j, j] = np.linalg.norm(v)
        Q[:, j] = v / R[j, j]
    return Q, R


def gram_schmidt_modified(A):
    """Modified Gram–Schmidt: same A = QR, but r_ij = q_i^T v uses the CURRENT, already-deflated v.
    Loses orthogonality like κ(A)·ε instead of κ(A)²·ε (Björck 1967; stated in §3.8, proved in LA-13)."""
    A = np.asarray(A, float); m, n = A.shape
    Q = np.zeros((m, n)); R = np.zeros((n, n))
    for j in range(n):
        v = A[:, j].copy()
        for i in range(j):
            R[i, j] = Q[:, i] @ v          # the only line that differs from the classical version
            v -= R[i, j] * Q[:, i]
        R[j, j] = np.linalg.norm(v)
        Q[:, j] = v / R[j, j]
    return Q, R


# Example 4.3: Gram–Schmidt on (1,1,1,1), (1,3,3,1), (1,1,1,4) in R^4.
A3 = np.array([[1, 1, 1], [1, 3, 1], [1, 3, 1], [1, 1, 4]], float)
Q3_expected = np.array([[.5, -.5, -1 / np.sqrt(2)], [.5, .5, 0], [.5, .5, 0], [.5, -.5, 1 / np.sqrt(2)]])
R3_expected = np.array([[2, 4, 3.5], [0, 2, -1.5], [0, 0, 3 / np.sqrt(2)]])
for gs in (gram_schmidt_classical, gram_schmidt_modified):
    Q3, R3 = gs(A3)
    assert np.allclose(Q3, Q3_expected) and np.allclose(R3, R3_expected)
    assert np.allclose(Q3.T @ Q3, np.eye(3)) and np.allclose(Q3 @ R3, A3)
print("Example 4.3: Q =\n", np.round(Q3, 6), "\nR =\n", np.round(R3, 6), "  (r33 = 3/√2 ≈ 2.1213)")

# Random well-conditioned matrices: both versions agree with each other and give A = QR, Q^TQ = I.
for _ in range(50):
    mm, nn = rng.integers(2, 8), rng.integers(1, 6); nn = min(mm, nn)
    A_rand = rng.standard_normal((mm, nn))
    Qc, Rc = gram_schmidt_classical(A_rand); Qm, Rm = gram_schmidt_modified(A_rand)
    assert np.allclose(Qc, Qm) and np.allclose(Rc, Rm)
    assert np.allclose(Qm.T @ Qm, np.eye(nn)) and np.allclose(Qm @ Rm, A_rand)
    assert np.allclose(Rm, np.triu(Rm)) and np.all(np.diag(Rm) > 0)
    # Theorem 3.16 / eq. (3.12): projection onto C(A) is Q Q^T, and it equals A (A^T A)^{-1} A^T.
    assert np.allclose(Qm @ Qm.T, project_onto_subspace(A_rand, np.zeros(mm))[3])

# Loss of orthogonality on Hilbert matrices H_ij = 1/(i + j − 1): cond(H_10) ≈ 1.6e13.
loss = lambda Q: np.linalg.norm(Q.T @ Q - np.eye(Q.shape[1]))
ns = np.arange(2, 13)
loss_c, loss_m, loss_h = [], [], []
for nh in ns:
    H = sla.hilbert(nh)
    loss_c.append(loss(gram_schmidt_classical(H)[0]))
    loss_m.append(loss(gram_schmidt_modified(H)[0]))
    loss_h.append(loss(np.linalg.qr(H)[0]))          # Householder QR (LA-13)
i10 = list(ns).index(10)
print(f"Hilbert 10×10, cond = {np.linalg.cond(sla.hilbert(10)):.2e}:  ‖QᵀQ − I‖  classical GS = {loss_c[i10]:.2e}, "
      f"modified GS = {loss_m[i10]:.2e}, Householder = {loss_h[i10]:.2e}")
assert loss_c[i10] > 1.0 and loss_m[i10] < 1e-2 and loss_h[i10] < 1e-13
assert all(c >= m_ for c, m_ in zip(loss_c[3:], loss_m[3:]))        # classical never beats modified once it matters

fig, ax = plt.subplots(figsize=(6, 4))
ax.semilogy(ns, loss_c, "o-", label="classical Gram–Schmidt"); ax.semilogy(ns, loss_m, "s-", label="modified Gram–Schmidt")
ax.semilogy(ns, loss_h, "^-", label="Householder (np.linalg.qr)")
ax.set_xlabel("n (Hilbert matrix size)"); ax.set_ylabel("‖QᵀQ − I‖"); ax.grid(alpha=.3, which="both"); ax.legend()
ax.set_title("Loss of orthogonality on Hilbert matrices")
fig.tight_layout(); fig.savefig(FIG / "gram_schmidt_loss.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 6: QR, and least squares via R x̂ = Qᵀb
#
# (3.15) A = QR.  (3.16) R x̂ = Qᵀb (derived from the normal equations in §3.10).
# `least_squares_qr` runs modified Gram–Schmidt on the augmented matrix [A | b]:
# its last column of R is Qᵀb and ‖e‖ computed *inside* the factorization, which
# is the trick that makes Gram–Schmidt least squares stable (Björck). The
# ill-conditioned fit at the end shows normal equations losing 11 digits where
# QR loses 5.


# %%
def qr_decompose(A):
    """Reduced QR factorization A = QR, eq. (3.15), by modified Gram–Schmidt.
    Sign convention: diag(R) > 0 (Theorem 3.19 makes it unique). NumPy's may differ by signs."""
    return gram_schmidt_modified(A)


def least_squares_qr(A, b):
    """x̂ = argmin ‖Ax − b‖² via QR: solve R x̂ = Qᵀb, eq. (3.16), by back substitution.
    Qᵀb and ‖b − Ax̂‖ are read off the QR of the augmented matrix [A | b] (§3.10)."""
    A = np.asarray(A, float); b = np.asarray(b, float); n = A.shape[1]
    _, R_aug = qr_decompose(np.c_[A, b])
    R, Qtb, resid_norm = R_aug[:n, :n], R_aug[:n, n], R_aug[n, n]
    return sla.solve_triangular(R, Qtb, lower=False), resid_norm


# Example 4.4(a): QR of the 3×2 from Example 4.1 and the same projection via R x̂ = Qᵀb.
Q4, R4 = qr_decompose(A1)
assert np.allclose(Q4, [[1 / np.sqrt(2), 1 / np.sqrt(6)], [0, 2 / np.sqrt(6)], [1 / np.sqrt(2), -1 / np.sqrt(6)]])
assert np.allclose(R4, [[np.sqrt(2), 1 / np.sqrt(2)], [0, np.sqrt(6) / 2]])
assert np.allclose(Q4.T @ b1, [2 * np.sqrt(2), 2 / np.sqrt(6)])
xh4, rn4 = least_squares_qr(A1, b1)
assert np.allclose(xh4, [5 / 3, 2 / 3]) and np.allclose(xh4, xh1)            # matches Example 4.1
assert np.allclose(Q4 @ (Q4.T @ b1), p1) and np.allclose(Q4 @ Q4.T, P1)     # QQᵀ is the same P
assert np.isclose(rn4, np.linalg.norm(e1)) and np.isclose(rn4, 4 / np.sqrt(3))
# Example 4.4(b): the 4×2 line fit of Example 4.2 via QR.
Q5, R5 = qr_decompose(A2)
assert np.allclose(R5, [[2, 3], [0, np.sqrt(5)]])
assert np.allclose(Q5[:, 1], np.array([-3, -1, 1, 3]) / (2 * np.sqrt(5)))
assert np.allclose(Q5.T @ y2, [4.5, 9 / (2 * np.sqrt(5))])
xh5, rn5 = least_squares_qr(A2, y2)
assert np.allclose(xh5, [0.9, 0.9]) and np.isclose(rn5 ** 2, 0.7)
print("Example 4.4: QR of the 3×2 gives x̂ =", xh4, "; QR of the 4×2 gives C, D =", xh5, " |r| =", round(rn5, 6))

# Against np.linalg.qr, up to the sign of each column of Q / row of R (§8).
for A_test in (A1, A2, A3, rng.standard_normal((7, 4))):
    Q_np, R_np = np.linalg.qr(A_test)
    D = np.diag(np.sign(np.diag(R_np)))              # D² = I, so (Q_np D)(D R_np) is also a QR
    Q_ours, R_ours = qr_decompose(A_test)
    assert np.allclose(Q_np @ D, Q_ours) and np.allclose(D @ R_np, R_ours)
    b_test = rng.standard_normal(A_test.shape[0])
    assert np.allclose(least_squares_qr(A_test, b_test)[0], np.linalg.lstsq(A_test, b_test, rcond=None)[0])
    assert np.allclose(least_squares_qr(A_test, b_test)[0], least_squares_normal(A_test, b_test))

# κ(AᵀA) = κ(A)² (stated in §3.10, proved in LA-13): check while both are computable.
for shift in (0.0, 10.0):
    A_k = design_matrix(shift + np.linspace(0, 5, 40), 2)
    assert np.isclose(np.linalg.cond(A_k.T @ A_k), np.linalg.cond(A_k) ** 2, rtol=1e-5)

# Ill-conditioned fit: a quadratic in t ∈ [100, 105] with columns 1, t, t² (κ(A) ≈ 6e7, κ(AᵀA) ≈ 3e15).
t_bad = 100 + np.linspace(0, 5, 40); A_bad = design_matrix(t_bad, 2)
x_true = np.array([1.0, -2.0, 0.5]); b_bad = A_bad @ x_true           # consistent system: exact answer known
rel = lambda x_: np.linalg.norm(x_ - x_true) / np.linalg.norm(x_true)
err_normal = rel(least_squares_normal(A_bad, b_bad))
Q_naive, R_naive = qr_decompose(A_bad)
err_naive = rel(sla.solve_triangular(R_naive, Q_naive.T @ b_bad))     # MGS, then Qᵀb formed separately
err_aug = rel(least_squares_qr(A_bad, b_bad)[0])                      # MGS on [A | b]
Q_hh, R_hh = np.linalg.qr(A_bad)
err_hh = rel(sla.solve_triangular(R_hh, Q_hh.T @ b_bad))              # Householder QR
err_lstsq = rel(np.linalg.lstsq(A_bad, b_bad, rcond=None)[0])
print(f"Ill-conditioned fit, κ(A) = {np.linalg.cond(A_bad):.1e}, κ(AᵀA) = {np.linalg.cond(A_bad.T @ A_bad):.1e}: "
      f"relative error  normal eqns {err_normal:.1e} | MGS then Qᵀb {err_naive:.1e} | MGS on [A|b] {err_aug:.1e} | "
      f"Householder {err_hh:.1e} | lstsq {err_lstsq:.1e}")
assert err_normal > 1e3 * err_aug and err_normal > 1e3 * err_hh
assert err_naive > 1e3 * err_aug      # naive MGS least squares is no better than the normal equations here (§8)

# %% [markdown]
# ## Part 7: Orthogonal matrices and Householder reflections
#
# Theorem 3.18: Q square with QᵀQ = I preserves lengths and dot products.
# Rotations, reflections, permutations. (3.17) H = I − 2vvᵀ/vᵀv is symmetric and
# orthogonal (Theorem 3.21), reflects across the hyperplane v^⊥, and the choice
# v = x − ‖x‖e₁ sends x to ‖x‖e₁ (the first step of Householder QR, LA-13).


# %%
def householder(v):
    """Householder reflection H = I − 2 v vᵀ / (vᵀ v), eq. (3.17). Symmetric, orthogonal, H v = −v, H u = u for u ⟂ v."""
    v = np.asarray(v, float).reshape(-1, 1)
    return np.eye(v.size) - 2 * (v @ v.T) / (v.T @ v)


theta = 0.7
Rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
Ref = np.array([[np.cos(theta), np.sin(theta)], [np.sin(theta), -np.cos(theta)]])
Perm = np.eye(4)[[2, 0, 3, 1]]
Q_big = np.linalg.qr(rng.standard_normal((5, 5)))[0]
for Q_ in (Rot, Ref, Perm, Q_big, Q_big.T, Q_big @ Q_big):          # transposes and products stay orthogonal (B7)
    k = Q_.shape[0]
    assert np.allclose(Q_.T @ Q_, np.eye(k)) and np.allclose(Q_ @ Q_.T, np.eye(k))
    for _ in range(20):
        u, w = rng.standard_normal(k), rng.standard_normal(k)
        assert np.isclose(np.linalg.norm(Q_ @ u), np.linalg.norm(u))          # lengths
        assert np.isclose((Q_ @ u) @ (Q_ @ w), u @ w)                           # dot products, hence angles
    # Parseval (Theorem 3.17): ‖b‖² = Σ (q_i·b)² for an orthonormal basis.
    assert np.isclose(np.sum((Q_.T @ u) ** 2), u @ u)
assert np.allclose(Ref @ Ref, np.eye(2)) and not np.allclose(Rot @ Rot, np.eye(2))   # a reflection is its own inverse

for k in (2, 3, 6):
    v = rng.standard_normal(k); H = householder(v)
    assert np.allclose(H.T, H) and np.allclose(H @ H, np.eye(k))     # symmetric and orthogonal (so H⁻¹ = H)
    assert np.allclose(H @ v, -v)
    u = rng.standard_normal(k); u -= (u @ v) / (v @ v) * v           # u ⟂ v
    assert np.allclose(H @ u, u)
    x = rng.standard_normal(k); e1 = np.eye(k)[0]
    v_x = x - np.linalg.norm(x) * e1
    assert np.allclose(householder(v_x) @ x, np.linalg.norm(x) * e1) # x is sent to the first axis
H34 = householder(np.array([3, 4.]) - 5 * np.array([1, 0.]))
assert np.allclose(H34 @ [3, 4], [5, 0])
print("Part 7: orthogonal matrices preserve lengths and angles; Householder H(3,4 → 5,0) =\n", np.round(H34, 4))

# Figure: reflect a shape across the mirror line v^⊥ in R².
v2 = np.array([1.0, 2.0]); H2 = householder(v2)
shape = np.array([[0.5, 0.5], [2.5, 0.5], [2.5, 1.0], [1.2, 1.0], [1.2, 1.5], [2.0, 1.5], [2.0, 2.0], [0.5, 2.0], [0.5, 0.5]]).T
refl = H2 @ shape
fig, ax = plt.subplots(figsize=(5.5, 5.5))
ax.plot(*shape, "C0-", lw=2, label="original"); ax.plot(*refl, "C1-", lw=2, label="H × original")
mirror = np.array([-v2[1], v2[0]]) / np.linalg.norm(v2)
ax.plot([-3.5 * mirror[0], 3.5 * mirror[0]], [-3.5 * mirror[1], 3.5 * mirror[1]], "k--", label="mirror: v^⊥")
ax.quiver(0, 0, *v2, angles="xy", scale_units="xy", scale=1, color="C3", label="v")
for p_, q_ in zip(shape.T[:-1], refl.T[:-1]):
    ax.plot([p_[0], q_[0]], [p_[1], q_[1]], color="gray", lw=.5)
ax.set_aspect("equal"); ax.grid(alpha=.3); ax.legend(fontsize=8); ax.set_xlim(-3.5, 3.5); ax.set_ylim(-3.5, 3.5)
ax.set_title("Householder reflection H = I − 2vvᵀ/vᵀv")
fig.tight_layout(); fig.savefig(FIG / "householder_2d.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 8: Orthogonal functions, the start of Fourier series (F2)
#
# With ⟨f, g⟩ = ∫_{−π}^{π} f g dx, the functions sin x, sin 2x, ... are
# orthogonal, eq. (3.18), each with ‖sin nx‖² = π. So the projection of f onto
# their span has coefficients ⟨f, sin nx⟩/π: the same q_i·b of (3.12), with the
# integral as the dot product. We check (3.18) symbolically and project a square
# wave.


# %%
xsym = sp.symbols("x", real=True)
for mm in range(1, 5):
    for nn in range(1, 5):
        val = sp.integrate(sp.sin(mm * xsym) * sp.sin(nn * xsym), (xsym, -sp.pi, sp.pi))
        assert sp.simplify(val - (sp.pi if mm == nn else 0)) == 0, (mm, nn, val)
# Square wave f(x) = sign(x) on (−π, π): coefficients b_n = ⟨f, sin nx⟩/π = 4/(nπ) for odd n, 0 for even n.
coef = {}
for nn in range(1, 8):
    b_n = 2 * sp.integrate(sp.sin(nn * xsym), (xsym, 0, sp.pi)) / sp.pi     # f is odd, so ∫_{-π}^{π} = 2∫_0^{π}
    coef[nn] = b_n
    assert sp.simplify(b_n - (4 / (nn * sp.pi) if nn % 2 else 0)) == 0
print("Part 8: ⟨sin mx, sin nx⟩ = π δ_mn verified for m, n ≤ 4; square-wave coefficients b_n =",
      {k: str(v) for k, v in coef.items()})

xx = np.linspace(-np.pi, np.pi, 800); sq = np.sign(xx)
fig, ax = plt.subplots(figsize=(6, 4)); ax.plot(xx, sq, "k", lw=1, label="square wave")
for N in (1, 3, 7):
    proj = sum(float(coef[nn]) * np.sin(nn * xx) for nn in range(1, N + 1))
    ax.plot(xx, proj, label=f"projection onto sin x … sin {N}x")
ax.grid(alpha=.3); ax.legend(fontsize=8); ax.set_title("Orthogonal projection onto a trigonometric basis (F2 preview)")
fig.tight_layout(); fig.savefig(FIG / "fourier_preview.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Your turn
#
# 1. **Weighted least squares.** Minimize Σ w_i r_i² with weights w_i > 0 (a
#    reading you trust more gets a bigger weight). Show that scaling row i of A
#    and entry i of b by √w_i turns this into an ordinary least-squares problem,
#    then implement `least_squares_weighted(A, b, w)` and check: all weights 1
#    reproduces `least_squares_qr`; a weight of 1e6 on one reading pulls the fit
#    through that point. (Exercise B8; a reference is below.)
# 2. **Damped oscillation envelope.** Generate y(t) = A e^{−γt} cos(ωt) + noise,
#    pick out the peaks, and fit ln|y_peak| = ln A − γ t_peak by least squares.
#    This is the log-decrement method from vibration testing (Exercise C3).
#    Compare the fitted γ with the true one; try noise levels 0.01 and 0.1.
# 3. **P and I − P.** For a random 6×3 A, build P and check P + (I − P) = I,
#    P(I − P) = 0, and that (I − P) is the projection onto N(Aᵀ) (compare with
#    `four_subspaces`). Then show P b and (I − P) b are orthogonal for random b.
# 4. **Centering.** Refit Example 4.2 with the column t replaced by t − mean(t).
#    AᵀA becomes diagonal. Why? Which theorem in the notes says the fitted line
#    cannot change? (Exercise B5.)
# 5. **Hilbert, again.** Re-run Part 5 with `gram_schmidt_classical` applied
#    twice ("CGS2": orthogonalize, then orthogonalize the result again). It
#    recovers orthogonality to machine precision. Why does the second pass work
#    when the first did not?


# %%
def least_squares_weighted(A, b, w):
    """argmin Σ w_i (Ax − b)_i² = argmin ‖W^{1/2}(Ax − b)‖²: scale row i by √w_i, then ordinary least squares.
    Exercise B8 reference."""
    s = np.sqrt(np.asarray(w, float))
    return least_squares_qr(s[:, None] * np.asarray(A, float), s * np.asarray(b, float))[0]


assert np.allclose(least_squares_weighted(A2, y2, np.ones(4)), xh2)
w_heavy = np.array([1, 1, 1e6, 1.0])
xw = least_squares_weighted(A2, y2, w_heavy)
assert abs((A2 @ xw)[2] - y2[2]) < 1e-5                                # the fit passes through point 3
assert not np.allclose(xw, xh2)

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
