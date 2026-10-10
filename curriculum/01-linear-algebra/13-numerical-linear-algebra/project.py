# %% [markdown]
# # LA-13 Project: Floating Point, Conditioning, Stability, and the Big Algorithms
#
# This script checks every number in the LA-13 notes and builds the module's
# algorithms from scratch before comparing them to NumPy/SciPy. Parts:
# (1) IEEE double and the rounding model fl(x∘y) = (x∘y)(1+δ); catastrophic
# cancellation in the quadratic formula in 8-digit arithmetic (Example 4.1);
# (2) conditioning: the 2×2 with κ = 10⁴ (Example 4.2), forward error ≈ κ ×
# backward error for LU across twelve decades of κ, growth factors;
# (3) Householder reflections (Example 4.3), the sign choice, Householder QR vs
# CGS vs MGS on Hilbert matrices, normal equations vs QR at κ(A) = 10⁸;
# (4) eigenvalues: Wilkinson's polynomial, the unshifted QR algorithm on a 2×2
# (Example 4.4), QR = simultaneous iteration, tridiagonalization, the practical
# QR algorithm with Wilkinson shifts and deflation, Golub–Kahan bidiagonalization;
# (5) iterative methods: the 2D Laplacian as a sparse matrix, Jacobi / Gauss–Seidel /
# CG, CG by hand (Example 4.5), the Chebyshev bound, Arnoldi + GMRES;
# (6) randomized SVD and the Halko–Martinsson–Tropp bounds.
# "Your turn" at the end: Rayleigh quotient iteration, preconditioned CG, Lanczos.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.linalg as sla
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import sympy as sp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T_START = time.time()

# A fixed categorical palette (same hue order in every figure).
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
plt.rcParams.update({"axes.prop_cycle": plt.cycler(color=PAL), "lines.linewidth": 1.6,
                     "axes.grid": True, "grid.alpha": 0.25, "legend.fontsize": 8,
                     "axes.titlesize": 10, "axes.labelsize": 9})

# Convention (notes §0 box): ε_mach = 2^-53 ≈ 1.11e-16 is the unit roundoff.
# NumPy's np.finfo(float).eps = 2^-52 is the *gap* from 1 to the next double = 2 ε_mach.
EPS = 2.0 ** -53

# %% [markdown]
# ## Part 1: IEEE double precision and the rounding model (§3.1)
#
# The 1 + ε == 1 search finds the gap between 1 and the next double, 2^-52.
# The largest relative rounding error is half that gap, ε_mach = 2^-53 (eq. 3.1).
# We check fl(x∘y) = (x∘y)(1+δ), |δ| ≤ ε_mach, against *exact* rational
# arithmetic (Python's `fractions`), for +, -, ×, ÷ on random doubles.


# %%
def machine_epsilon_search():
    """Halve ε until 1 + ε rounds back to 1; return the last ε that did not."""
    eps = 1.0
    while 1.0 + eps / 2 != 1.0:
        eps /= 2
    return eps


gap = machine_epsilon_search()
assert gap == 2.0 ** -52 == np.finfo(float).eps
assert 1.0 + 2.0 ** -53 == 1.0          # exactly half the gap: a tie, rounds to even (= 1)
assert 1.0 + 2.0 ** -53 * 1.000001 > 1.0  # anything past the tie rounds up
fi = np.finfo(np.float64)
print(f"gap 1→next double = {gap:.3e} = 2^-52;  ε_mach = 2^-53 = {EPS:.3e}")
print(f"significand bits = {fi.nmant + 1} (52 stored + hidden 1); range ≈ [{fi.tiny:.2e}, {fi.max:.2e}]")
assert fi.nmant == 52

# Rounding model check against exact rational arithmetic.
worst = 0.0
for _ in range(2000):
    x, y = rng.standard_normal(2) * 10.0 ** rng.integers(-5, 6, size=2)
    for op, fop in (("+", lambda a, b: a + b), ("-", lambda a, b: a - b),
                    ("*", lambda a, b: a * b), ("/", lambda a, b: a / b)):
        exact = fop(Fraction(x), Fraction(y))
        if exact == 0:
            continue
        computed = fop(x, y)
        delta = abs((Fraction(computed) - exact) / exact)
        worst = max(worst, float(delta))
        assert delta <= Fraction(EPS) * (1 + Fraction(1, 10**9)), (op, x, y, float(delta))
print(f"rounding model fl(x∘y) = (x∘y)(1+δ): worst |δ| over 8000 operations = {worst:.3e} ≤ ε_mach ✓")

# %% [markdown]
# ### Example 4.1: catastrophic cancellation in the quadratic formula, 8 significant digits
#
# x² + 9000x + 1 = 0. The small root is -c/(big root) ≈ -1.1111111e-4.
# In 8-digit arithmetic (Python's `decimal` with prec = 8 does exactly this),
# -b + sqrt(b² - 4ac) subtracts two numbers that agree to 8 digits and keeps 1.


# %%
def quadratic_8digit(a, b, c):
    """Both roots of ax²+bx+c by the textbook formula, every operation rounded to 8 significant digits.
    Also the stable small root c/(a·r_big) (eq. 3.5)."""
    getcontext().prec = 8
    a, b, c = Decimal(a), Decimal(b), Decimal(c)
    disc = +(+(b * b) - +(4 * a * c))          # unary + rounds to the current precision
    s = +disc.sqrt()
    r_small_naive = +(+(-b + s) / (2 * a))
    r_big = +(+(-b - s) / (2 * a))
    r_small_stable = +(c / (a * r_big))
    return disc, s, r_small_naive, r_big, r_small_stable


disc, s, r_naive, r_big, r_stable = quadratic_8digit("1", "9000", "1")
true_small = -2.0 / (9000.0 + np.sqrt(9000.0**2 - 4))   # the stable formula in double precision
print(f"disc = {disc}, sqrt = {s}, naive small root = {r_naive}, big root = {r_big}, stable small root = {r_stable}")
print(f"true small root = {true_small:.10e}")
assert str(disc) == "80999996" and str(s) == "8999.9998"
assert str(r_naive) == "-0.0001" and str(r_stable) == "-0.00011111111"
rel_naive = abs(float(r_naive) - true_small) / abs(true_small)
rel_stable = abs(float(r_stable) - true_small) / abs(true_small)
print(f"relative error: naive = {rel_naive:.3f} (one correct digit), stable = {rel_stable:.1e} (all eight: that is the 8-digit rounding itself)")
assert 0.09 < rel_naive < 0.11 and rel_stable < 5e-8   # 5e-8 = half a unit in the 8th digit
# Conditioning of the subtraction itself: κ_sub = (|x|+|y|)/|x-y|, eq. (3.4)
kappa_sub = (9000 + 8999.99977778) / abs(9000 - 8999.99977778)
print(f"κ of the subtraction 9000 - 8999.9998 = {kappa_sub:.1e}; times the 8-digit rounding level 5e-9 ≈ {kappa_sub*5e-9:.2f}")
assert 7e7 < kappa_sub < 9e7

# %% [markdown]
# ## Part 2: Conditioning and backward stability (§3.2–3.3)
#
# Example 4.2: A = Q diag(1, 1e-4) Qᵀ with Q the 45° rotation, so κ₂(A) = 10⁴
# exactly. The bound ‖δx‖/‖x‖ ≤ κ(A)‖δb‖/‖b‖ (Theorem 3.3) is *attained* when
# δb points along the small singular direction and loose by 10⁴ when it points
# along the large one.

# %%
A42 = np.array([[0.50005, 0.49995], [0.49995, 0.50005]])
b42 = np.array([1.0, 1.0])
kappa42 = np.linalg.cond(A42)
x42 = np.linalg.solve(A42, b42)
assert np.isclose(kappa42, 1e4) and np.allclose(x42, [1, 1])
db_bad = np.array([1e-4, -1e-4])     # along (1,-1): the σ_min direction
db_good = np.array([1e-4, 1e-4])     # along (1, 1): the σ_max direction
for name, db in (("worst-case δb ∥ (1,-1)", db_bad), ("benign δb ∥ (1,1)", db_good)):
    dx = np.linalg.solve(A42, b42 + db) - x42
    rel_in = np.linalg.norm(db) / np.linalg.norm(b42)
    rel_out = np.linalg.norm(dx) / np.linalg.norm(x42)
    print(f"Example 4.2 {name}: ‖δb‖/‖b‖ = {rel_in:.1e}, ‖δx‖/‖x‖ = {rel_out:.3e}, bound κ·‖δb‖/‖b‖ = {kappa42*rel_in:.3e}")
    assert rel_out <= kappa42 * rel_in * (1 + 1e-9)
x_bad = np.linalg.solve(A42, b42 + db_bad)
assert np.allclose(x_bad, [2, 0]) and np.isclose(np.linalg.norm(x_bad - x42) / np.linalg.norm(x42), 1.0)
print("   worst case: x = (1,1) becomes x̃ = (2,0): 100% relative error from a 0.01% perturbation.")

# %% [markdown]
# ### Forward error ≲ κ × backward error (Theorem 3.5) for LU with partial pivoting
#
# For κ from 1 to 10¹² we solve Ax = b with `np.linalg.solve` (LAPACK getrf/getrs).
# The backward error ‖b - Ax̂‖/(‖A‖‖x̂‖) stays at ε_mach regardless of κ; the
# forward error ‖x̂ - x‖/‖x‖ grows like κ ε_mach.


# %%
def matrix_with_condition(n, kappa, rng):
    """Random n×n matrix A = U diag(σ) Vᵀ with σ log-spaced in [1/κ, 1], so κ₂(A) = κ."""
    U, _ = np.linalg.qr(rng.standard_normal((n, n)))
    V, _ = np.linalg.qr(rng.standard_normal((n, n)))
    sig = np.logspace(0, -np.log10(kappa), n)
    return (U * sig) @ V.T


n = 60
kappas = np.logspace(0, 12, 13)
fwd, bwd = [], []
for kap in kappas:
    A = matrix_with_condition(n, kap, rng)
    x_true = rng.standard_normal(n)
    b = A @ x_true
    x_hat = np.linalg.solve(A, b)
    fwd.append(np.linalg.norm(x_hat - x_true) / np.linalg.norm(x_true))
    bwd.append(np.linalg.norm(b - A @ x_hat) / (np.linalg.norm(A, 2) * np.linalg.norm(x_hat)))
fwd, bwd = np.array(fwd), np.array(bwd)
assert np.all(bwd < 50 * n * EPS), "LU should be backward stable here"
assert np.all(fwd < 50 * n * kappas * EPS), "forward error should be ≲ κ ε"
assert fwd[-1] > 1e-6 and fwd[0] < 1e-13, "forward error must actually grow with κ"
print(f"LU: backward error in [{bwd.min():.1e}, {bwd.max():.1e}] for all κ; forward error {fwd[0]:.1e} at κ=1 → {fwd[-1]:.1e} at κ=1e12")

fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(kappas, fwd, "o-", label="forward error ‖x̂−x‖/‖x‖")
ax.loglog(kappas, bwd, "s-", label="backward error ‖b−Ax̂‖/(‖A‖‖x̂‖)")
ax.loglog(kappas, kappas * EPS, "--", color="gray", label="κ · ε_mach")
ax.set_xlabel("κ₂(A)"); ax.set_ylabel("relative error"); ax.legend()
ax.set_title("Solve Ax=b by LU (n=60): forward error ≈ κ × backward error")
fig.tight_layout(); fig.savefig(FIG / "forward_vs_kappa.png", dpi=120); plt.close(fig)


# %% [markdown]
# ### Growth factor of Gaussian elimination with partial pivoting (§3.3)
#
# ρ = max|u_ij| / max|a_ij|. Wilkinson's bound ρ ≤ 2^{n-1} is attained by the
# matrix with 1 on the diagonal, -1 below it, and a last column of 1s. Random
# matrices show growth of order n^{1/2}–n^{2/3}, nowhere near 2^{n-1}.

# %%
def lu_partial_pivoting(A):
    """PA = LU by Gaussian elimination with partial pivoting; returns P, L, U and the growth factor ρ (eq. 3.9)."""
    A = np.array(A, float); n = A.shape[0]
    U = A.copy(); L = np.eye(n); P = np.eye(n)
    for k in range(n - 1):
        p = k + np.argmax(np.abs(U[k:, k]))
        if p != k:
            U[[k, p], k:] = U[[p, k], k:]; L[[k, p], :k] = L[[p, k], :k]; P[[k, p]] = P[[p, k]]
        for i in range(k + 1, n):
            L[i, k] = U[i, k] / U[k, k]
            U[i, k:] -= L[i, k] * U[k, k:]
    rho = np.abs(U).max() / np.abs(A).max()
    return P, L, U, rho


def wilkinson_growth_matrix(n):
    W = np.eye(n) - np.tril(np.ones((n, n)), -1)
    W[:, -1] = 1.0
    return W


for nn in (10, 20):
    W = wilkinson_growth_matrix(nn)
    P, L, U, rho = lu_partial_pivoting(W)
    assert np.allclose(P @ W, L @ U) and np.isclose(rho, 2.0 ** (nn - 1))
    print(f"Wilkinson growth matrix n={nn}: ρ = {rho:.0f} = 2^{nn-1}")
A = rng.standard_normal((100, 100))
P, L, U, rho = lu_partial_pivoting(A)
assert np.allclose(P @ A, L @ U)
print(f"random Gaussian 100×100: ρ = {rho:.1f}  (vs worst case 2^99 ≈ {2.0**99:.1e})")
assert rho < 100

# %% [markdown]
# ## Part 3: Householder reflections and QR (§3.4–3.5)
#
# Example 4.3: x = (2, 1, 2)ᵀ, ‖x‖ = 3, v = x + sign(x₁)‖x‖e₁ = (5, 1, 2)ᵀ,
# vᵀv = 30, Hx = -3e₁. H = I - vvᵀ/15.


# %%
def householder_vector(x):
    """v such that H = I - 2vvᵀ/vᵀv reflects x onto -sign(x₁)‖x‖ e₁ (eq. 3.12, sign chosen to avoid cancellation)."""
    x = np.asarray(x, float)
    v = x.copy()
    v[0] += (1.0 if x[0] >= 0 else -1.0) * np.linalg.norm(x)
    return v


def householder_matrix(v):
    v = np.asarray(v, float)
    return np.eye(len(v)) - 2.0 * np.outer(v, v) / (v @ v)


x43 = np.array([2.0, 1.0, 2.0])
v43 = householder_vector(x43)
H43 = householder_matrix(v43)
assert np.allclose(v43, [5, 1, 2]) and np.isclose(v43 @ v43, 30)
assert np.allclose(H43 @ x43, [-3, 0, 0])
assert np.allclose(H43, np.array([[-10, -5, -10], [-5, 14, -2], [-10, -2, 11]]) / 15)
assert np.allclose(H43 @ H43, np.eye(3)) and np.allclose(H43, H43.T)   # involutory, symmetric → orthogonal
print("Example 4.3: v =", v43, " vᵀv =", v43 @ v43, " Hx =", H43 @ x43)
# The other sign also works in exact arithmetic: v' = x - ‖x‖e₁ = (-1, 1, 2), H'x = +3e₁.
v_other = x43 - 3 * np.eye(3)[0]
assert np.allclose(householder_matrix(v_other) @ x43, [3, 0, 0])

# Why the sign matters: x almost parallel to e₁. The "wrong" sign computes v₁ = x₁ - ‖x‖
# as a difference of nearly equal numbers; v points in the wrong direction and Hx misses.
x_near = np.array([1.0, 1e-8, 1e-8])
nx = np.linalg.norm(x_near)
v_bad = x_near.copy(); v_bad[0] -= nx           # cancellation: 1 - (1 + 1e-16) ≈ lost
v_good = householder_vector(x_near)              # 1 + (1 + 1e-16): no cancellation
err_bad = np.linalg.norm(householder_matrix(v_bad) @ x_near - nx * np.eye(3)[0]) / nx
err_good = np.linalg.norm(householder_matrix(v_good) @ x_near + nx * np.eye(3)[0]) / nx
print(f"sign choice on x = (1, 1e-8, 1e-8): ‖Hx − target‖/‖x‖ = {err_bad:.1e} (bad sign) vs {err_good:.1e} (good sign)")
assert err_bad > 1e-10 and err_good < 1e-15
# Both H's are exactly orthogonal; the bad one is a perfect reflection across the wrong plane.
assert np.allclose(householder_matrix(v_bad).T @ householder_matrix(v_bad), np.eye(3))


# %% [markdown]
# ### Householder QR, classical and modified Gram–Schmidt (Algorithm 3.1; §3.5)
#
# **Build it yourself first:** write `householder_qr(A)` returning Q (m×m) and R
# (m×n) by applying n reflections to A and then, in reverse order, to I.


# %%
def householder_qr(A):
    """Full QR factorization A = QR by Householder reflections (Algorithm 3.1, T&B Alg. 10.1).
    Q is m×m orthogonal (formed as Q₁Q₂…Qₙ applied to I), R is m×n upper triangular."""
    R = np.array(A, float); m, n = R.shape
    vs = []
    for k in range(min(n, m)):
        v = householder_vector(R[k:, k])
        if v @ v > 0:
            R[k:, k:] -= 2.0 * np.outer(v, v @ R[k:, k:]) / (v @ v)
        vs.append(v)
    Q = np.eye(m)
    for k in reversed(range(len(vs))):
        v = vs[k]
        if v @ v > 0:
            Q[k:, :] -= 2.0 * np.outer(v, v @ Q[k:, :]) / (v @ v)
    R = np.triu(R)
    return Q, R


def gram_schmidt_classical(A):
    """Reduced QR by classical Gram–Schmidt: q_j = (a_j - Σ_{i<j} (q_iᵀa_j) q_i)/‖·‖ (LA-5)."""
    A = np.array(A, float); m, n = A.shape
    Q = np.zeros((m, n)); R = np.zeros((n, n))
    for j in range(n):
        v = A[:, j].copy()
        for i in range(j):
            R[i, j] = Q[:, i] @ A[:, j]        # projections computed from the ORIGINAL column
            v -= R[i, j] * Q[:, i]
        R[j, j] = np.linalg.norm(v); Q[:, j] = v / R[j, j]
    return Q, R


def gram_schmidt_modified(A):
    """Reduced QR by modified Gram–Schmidt: project the *updated* vector, one q at a time (T&B Alg. 8.1)."""
    V = np.array(A, float); m, n = V.shape
    Q = np.zeros((m, n)); R = np.zeros((n, n))
    for i in range(n):
        R[i, i] = np.linalg.norm(V[:, i]); Q[:, i] = V[:, i] / R[i, i]
        for j in range(i + 1, n):
            R[i, j] = Q[:, i] @ V[:, j]
            V[:, j] -= R[i, j] * Q[:, i]
    return Q, R


# Correctness on a well-conditioned random matrix.
A = rng.standard_normal((7, 4))
Q, R = householder_qr(A)
assert np.allclose(Q @ R, A) and np.allclose(Q.T @ Q, np.eye(7)) and np.allclose(np.tril(R, -1), 0)
for f in (gram_schmidt_classical, gram_schmidt_modified):
    Qr, Rr = f(A)
    assert np.allclose(Qr @ Rr, A) and np.allclose(Qr.T @ Qr, np.eye(4))

# Orthogonality loss on Hilbert matrices n = 2..14.
ns = np.arange(2, 15)
loss = {"CGS": [], "MGS": [], "Householder": []}
conds = []
for nn in ns:
    Hn = sla.hilbert(nn); conds.append(np.linalg.cond(Hn))
    Qc, _ = gram_schmidt_classical(Hn); Qm, _ = gram_schmidt_modified(Hn); Qh, Rh = householder_qr(Hn)
    loss["CGS"].append(np.linalg.norm(Qc.T @ Qc - np.eye(nn)))
    loss["MGS"].append(np.linalg.norm(Qm.T @ Qm - np.eye(nn)))
    loss["Householder"].append(np.linalg.norm(Qh.T @ Qh - np.eye(nn)))
    assert np.linalg.norm(Qh @ Rh - Hn) / np.linalg.norm(Hn) < 1e-14   # backward stable: residual ~ ε
conds = np.array(conds)
for k in loss:
    loss[k] = np.array(loss[k])
assert loss["Householder"].max() < 1e-14
assert loss["CGS"][ns == 10][0] > 1 and loss["MGS"][ns == 10][0] < 1e-2 and loss["MGS"][ns == 10][0] > loss["Householder"][ns == 10][0] * 1e8
print(f"Hilbert n=10 (κ={conds[ns==10][0]:.1e}): ‖QᵀQ−I‖ = CGS {loss['CGS'][ns==10][0]:.1e}, MGS {loss['MGS'][ns==10][0]:.1e}, Householder {loss['Householder'][ns==10][0]:.1e}")

fig, ax = plt.subplots(figsize=(6.5, 4))
for name, mk in zip(loss, ("o", "s", "^")):
    ax.semilogy(ns, np.maximum(loss[name], 1e-17), mk + "-", label=name)
ax.semilogy(ns, np.minimum(conds * EPS, 10), "--", color="gray", label="κ(H) · ε_mach")
ax.set_xlabel("n (Hilbert matrix size)"); ax.set_ylabel("‖QᵀQ − I‖₂"); ax.legend()
ax.set_title("Loss of orthogonality: CGS ~ κ²ε, MGS ~ κε, Householder ~ ε")
fig.tight_layout(); fig.savefig(FIG / "orthogonality_loss.png", dpi=120); plt.close(fig)

# %% [markdown]
# ### Least squares: normal equations (κ²) vs QR (κ), Theorem 3.8
#
# A is 60×12 with singular values log-spaced from 1 to 1e-8, so κ(A) = 1e8 and
# κ(AᵀA) = 1e16 ≈ 1/ε_mach. The normal equations lose everything; QR keeps 8 digits.

# %%
m, n = 60, 12
U, _ = np.linalg.qr(rng.standard_normal((m, n)))
V, _ = np.linalg.qr(rng.standard_normal((n, n)))
sig = np.logspace(0, -8, n)
A_ls = (U * sig) @ V.T
assert np.isclose(np.linalg.cond(A_ls), 1e8, rtol=1e-6)
assert np.isclose(np.linalg.cond(A_ls.T @ A_ls), 1e16, rtol=1e-1)   # κ(AᵀA) = κ(A)², to rounding
x_true = rng.standard_normal(n)
# A tiny residual ⟂ range(A). (With a LARGE residual the *problem* itself has κ² sensitivity, T&B Thm 18.1,
# and no algorithm can beat it; keeping r small isolates the algorithmic gap.)
b_ls = A_ls @ x_true + 1e-10 * (np.eye(m) - U @ U.T) @ rng.standard_normal(m)
x_ne = np.linalg.solve(A_ls.T @ A_ls, A_ls.T @ b_ls)
Qf, Rf = householder_qr(A_ls)
x_qr = sla.solve_triangular(Rf[:n], (Qf.T @ b_ls)[:n])
x_ref = np.linalg.lstsq(A_ls, b_ls, rcond=None)[0]
err_ne = np.linalg.norm(x_ne - x_true) / np.linalg.norm(x_true)
err_qr = np.linalg.norm(x_qr - x_true) / np.linalg.norm(x_true)
err_ref = np.linalg.norm(x_ref - x_true) / np.linalg.norm(x_true)
print(f"least squares, κ(A)=1e8: error normal equations = {err_ne:.2e}, Householder QR = {err_qr:.2e}, lstsq (SVD) = {err_ref:.2e}")
assert err_ne > 1e-3 and err_qr < 1e-6 and err_ne / err_qr > 1e3

# %% [markdown]
# ## Part 4: Eigenvalue algorithms (§3.6–3.8)
#
# ### Wilkinson's polynomial: w(x) = Π_{k=1}^{20} (x − k)
# Perturb the x¹⁹ coefficient (−210) by −2⁻²³ ≈ −1.2e-7. The roots 10…19
# become five complex pairs; the pair near 16.73 has imaginary part ±2.81.
# We get the exact roots of the *perturbed* polynomial with SymPy at 30 digits,
# then show that np.roots in double precision gets similar numbers.

# %%
xs = sp.symbols("x")
wilk = sp.Poly(sp.prod([xs - k for k in range(1, 21)]), xs)
coeffs = wilk.all_coeffs()
assert coeffs[1] == -210 and coeffs[-1] == sp.factorial(20)
pert = list(coeffs); pert[1] = pert[1] - sp.Rational(1, 2**23)
roots_pert = [complex(r) for r in sp.Poly(pert, xs).nroots(n=30, maxsteps=200)]
roots_pert.sort(key=lambda z: (z.real, z.imag))
cplx = [z for z in roots_pert if abs(z.imag) > 1e-9]
print("perturbed Wilkinson roots (30-digit):")
for z in roots_pert:
    print(f"   {z.real:12.7f} {z.imag:+12.7f}i")
assert len(cplx) == 10
im_max = max(abs(z.imag) for z in cplx)
z_big = max(cplx, key=lambda z: abs(z.imag))
assert np.isclose(z_big.real, 16.7307375, atol=1e-6) and np.isclose(abs(z_big.imag), 2.8126249, atol=1e-6)
assert all(abs(z - round(z.real)) < 1e-6 for z in roots_pert[:4])    # roots 1..4 barely move
print(f"largest imaginary part = {im_max:.6f} at real part {z_big.real:.6f}; roots 1–4 unchanged to 1e-6")

c_float = np.array([float(v) for v in coeffs])
c_pert = c_float.copy(); c_pert[1] -= 2.0 ** -23
r_np_exact = np.roots(c_float)          # double-precision rounding of the coefficients already perturbs the roots
r_np_pert = np.roots(c_pert)
print(f"np.roots on the UNperturbed coefficients: max |root − nearest integer| = {np.abs(r_np_exact - np.round(r_np_exact.real)).max():.3f} (coefficients up to {c_float.max():.1e} > 2^53 are not even representable)")
assert np.abs(r_np_exact - np.round(r_np_exact.real)).max() > 1e-3
fig, ax = plt.subplots(figsize=(7, 3.6))
ax.plot(np.arange(1, 21), np.zeros(20), "o", ms=7, mfc="none", label="exact roots 1…20")
ax.plot([z.real for z in roots_pert], [z.imag for z in roots_pert], "x", ms=7, label="perturbed by 2⁻²³ in x¹⁹ (30-digit)")
ax.plot(r_np_pert.real, r_np_pert.imag, "+", ms=7, label="np.roots, double precision")
ax.axhline(0, color="gray", lw=0.6); ax.set_xlabel("Re"); ax.set_ylabel("Im"); ax.legend(loc="upper left")
ax.set_title("Wilkinson's polynomial: a 1e-7 coefficient change moves roots by ~3 in the complex plane")
fig.tight_layout(); fig.savefig(FIG / "wilkinson_roots.png", dpi=120); plt.close(fig)


# %% [markdown]
# ### Example 4.4: the unshifted QR algorithm on a 2×2 symmetric matrix
#
# A = [[3,1],[1,1]], eigenvalues 2 ± √2. A₁ = RQ = [[3.4, 0.2],[0.2, 0.6]],
# A₂ = [[99/29, 1/29],[1/29, 17/29]]. The off-diagonal shrinks by λ₂/λ₁ = 0.1716 per step.
# (QR with positive diagonal in R, as in the hand calculation.)


# %%
def qr_positive(A):
    """QR with R having a positive diagonal (fixes the sign ambiguity so hand and machine agree)."""
    Q, R = np.linalg.qr(A)
    s = np.sign(np.diag(R)); s[s == 0] = 1
    return Q * s, (R.T * s).T


def qr_algorithm_unshifted(A, steps):
    """A_{k+1} = R_k Q_k where A_k = Q_k R_k (eq. 3.17). Returns the iterates and the accumulated Q⁽ᵏ⁾ = Q₁…Q_k."""
    Ak = np.array(A, float); Qacc = np.eye(len(Ak)); hist = [Ak.copy()]
    for _ in range(steps):
        Q, R = qr_positive(Ak); Ak = R @ Q; Qacc = Qacc @ Q; hist.append(Ak.copy())
    return hist, Qacc


A44 = np.array([[3.0, 1.0], [1.0, 1.0]])
hist, Qacc = qr_algorithm_unshifted(A44, 4)
assert np.allclose(hist[1], [[3.4, 0.2], [0.2, 0.6]])
assert np.allclose(hist[2], np.array([[99, 1], [1, 17]]) / 29)
lam = np.sort(np.linalg.eigvalsh(A44))[::-1]
ratio = lam[1] / lam[0]
print(f"Example 4.4: A₁ = {hist[1].tolist()}, A₂ = {np.round(hist[2], 5).tolist()}; off-diagonals {[round(h[1,0],5) for h in hist]}")
print(f"   λ₂/λ₁ = {ratio:.4f}; successive off-diagonal ratios: {[round(hist[k+1][1,0]/hist[k][1,0], 4) for k in range(4)]}")
assert np.isclose(ratio, 3 - 2 * np.sqrt(2)) and np.isclose(hist[2][1, 0] / hist[1][1, 0], ratio, atol=2e-3)
# Theorem 3.11: A_{k+1} = (Q⁽ᵏ⁾)ᵀ A Q⁽ᵏ⁾ and A^k = Q⁽ᵏ⁾ R⁽ᵏ⁾, so Q⁽ᵏ⁾e₁ is normalized power iteration.
assert np.allclose(Qacc.T @ A44 @ Qacc, hist[4])
p = np.linalg.matrix_power(A44, 4) @ np.eye(2)[0]; p /= np.linalg.norm(p)
assert np.allclose(np.abs(Qacc[:, 0]), np.abs(p))
print("   Q⁽⁴⁾ e₁ equals normalized A⁴ e₁ (QR algorithm = simultaneous power iteration) ✓")

# Decay of the subdiagonal for a 6×6 symmetric matrix with known eigenvalues 1..6: entry (j+1,j) decays like (λ_{j+1}/λ_j)^k.
Q6, _ = np.linalg.qr(rng.standard_normal((6, 6)))
A6 = Q6 @ np.diag([6.0, 5.0, 4.0, 3.0, 2.0, 1.0]) @ Q6.T
hist6, _ = qr_algorithm_unshifted(A6, 60)
sub = np.array([[abs(h[j + 1, j]) for j in range(5)] for h in hist6])
fig, ax = plt.subplots(figsize=(6.5, 4))
for j in range(5):
    ax.semilogy(sub[:, j], label=f"|a_{j+2},{j+1}| ~ ({6-j-1}/{6-j})^k")
ax.set_xlabel("QR iteration k"); ax.set_ylabel("subdiagonal entry"); ax.legend()
ax.set_title("Unshifted QR on a 6×6 with eigenvalues 1..6: each subdiagonal decays at its own rate")
fig.tight_layout(); fig.savefig(FIG / "qr_subdiagonal_decay.png", dpi=120); plt.close(fig)
# Check the predicted geometric rates over iterations 20..60.
for j in range(5):
    rate = (sub[60, j] / sub[20, j]) ** (1 / 40)
    assert abs(rate - (5 - j) / (6 - j)) < 0.03, (j, rate)
print("   6×6 check: subdiagonal (j+1, j) decays at rate λ_{j+1}/λ_j for every j ✓")


# %% [markdown]
# ### The practical symmetric QR algorithm: tridiagonalize, Wilkinson shift, deflate (§3.8)
#
# **Build it yourself first:** `householder_tridiagonalize(A)` (n−2 reflections,
# applied from both sides) and `qr_algorithm_symmetric(A)`.


# %%
def householder_tridiagonalize(A):
    """Symmetric A = Q T Qᵀ with T tridiagonal, by n−2 two-sided Householder reflections (T&B Alg. 26.1)."""
    T = np.array(A, float); n = T.shape[0]; Q = np.eye(n)
    for k in range(n - 2):
        v = householder_vector(T[k + 1:, k])
        if v @ v == 0:
            continue
        v = v / np.linalg.norm(v)
        T[k + 1:, k:] -= 2.0 * np.outer(v, v @ T[k + 1:, k:])    # left
        T[:, k + 1:] -= 2.0 * np.outer(T[:, k + 1:] @ v, v)      # right
        Q[:, k + 1:] -= 2.0 * np.outer(Q[:, k + 1:] @ v, v)
    T = np.triu(np.tril(T, 1), -1)
    return 0.5 * (T + T.T), Q


def qr_step_tridiagonal(T, mu):
    """One explicit shifted QR step T − μI = QR, T' = RQ + μI on a symmetric tridiagonal block, by Givens
    rotations (eq. 3.19). Applied to full rows/columns for clarity: O(m²) per step here, O(m) with a banded code."""
    B = T - mu * np.eye(len(T)); G = []
    for k in range(len(B) - 1):
        a, b = B[k, k], B[k + 1, k]
        r = np.hypot(a, b)
        c, s = (1.0, 0.0) if r == 0 else (a / r, b / r)
        rk, rk1 = B[k, :].copy(), B[k + 1, :].copy()
        B[k, :] = c * rk + s * rk1; B[k + 1, :] = -s * rk + c * rk1
        G.append((c, s))
    for k, (c, s) in enumerate(G):                 # R Q = R G₀ᵀ G₁ᵀ …
        ck, ck1 = B[:, k].copy(), B[:, k + 1].copy()
        B[:, k] = c * ck + s * ck1; B[:, k + 1] = -s * ck + c * ck1
    B = np.triu(np.tril(B, 1), -1)
    return 0.5 * (B + B.T) + mu * np.eye(len(T))


def wilkinson_shift(T):
    """Eigenvalue of the trailing 2×2 block of T closer to its last entry (T&B eq. 29.10)."""
    a, b, c = T[-2, -2], T[-1, -2], T[-1, -1]
    delta = 0.5 * (a - c)
    sgn = 1.0 if delta >= 0 else -1.0
    return c - sgn * b * b / (abs(delta) + np.hypot(delta, b))


def qr_algorithm_symmetric(A, tol=1e-14):
    """Eigenvalues of symmetric A by Householder tridiagonalization + shifted QR with Wilkinson shifts and
    deflation (§3.8). Returns (sorted eigenvalues, number of QR steps)."""
    T, _ = householder_tridiagonalize(A)
    n = T.shape[0]; eigs = []; steps = 0; hi = n - 1
    while hi >= 0:
        lo = hi
        while lo > 0 and abs(T[lo, lo - 1]) > tol * (abs(T[lo, lo]) + abs(T[lo - 1, lo - 1])):
            lo -= 1
        if lo == hi:                                  # 1×1 block has converged: deflate
            eigs.append(T[hi, hi]); hi -= 1; continue
        blk = T[lo:hi + 1, lo:hi + 1]
        T[lo:hi + 1, lo:hi + 1] = qr_step_tridiagonal(blk, wilkinson_shift(blk))
        steps += 1
    return np.sort(np.array(eigs)), steps


A = rng.standard_normal((8, 8)); A = A + A.T
T, Qt = householder_tridiagonalize(A)
assert np.allclose(Qt @ T @ Qt.T, A) and np.allclose(Qt.T @ Qt, np.eye(8))
assert np.allclose(np.triu(T, 2), 0) and np.allclose(np.tril(T, -2), 0)
# One shifted step preserves eigenvalues (similarity) and keeps the tridiagonal form.
T1 = qr_step_tridiagonal(T.copy(), wilkinson_shift(T))
assert np.allclose(np.sort(np.linalg.eigvalsh(T1)), np.sort(np.linalg.eigvalsh(T)))
for nn in (10, 50, 150):
    A = rng.standard_normal((nn, nn)); A = A + A.T
    t0 = time.time(); ev, steps = qr_algorithm_symmetric(A); t_qr = time.time() - t0
    ev_ref = np.linalg.eigvalsh(A)
    err = np.abs(ev - ev_ref).max() / np.abs(ev_ref).max()
    assert err < 1e-12, err
    print(f"qr_algorithm_symmetric n={nn}: {steps} QR steps = {steps/nn:.2f} per eigenvalue, max rel error vs eigh {err:.1e}, {t_qr:.2f}s")
    assert steps < 4 * nn

# Golub–Kahan bidiagonalization: A = U B Vᵀ with B upper bidiagonal; σ(B) = σ(A). (SVD sketch, §3.8.)
def golub_kahan_bidiagonalize(A):
    """A = U B Vᵀ with B upper bidiagonal, by alternating left and right Householder reflections (T&B Alg. 31.1)."""
    B = np.array(A, float); m, n = B.shape; U = np.eye(m); V = np.eye(n)
    for k in range(n):
        v = householder_vector(B[k:, k]); v /= np.linalg.norm(v)
        B[k:, k:] -= 2.0 * np.outer(v, v @ B[k:, k:]); U[:, k:] -= 2.0 * np.outer(U[:, k:] @ v, v)
        if k < n - 2:
            w = householder_vector(B[k, k + 1:]); w /= np.linalg.norm(w)
            B[k:, k + 1:] -= 2.0 * np.outer(B[k:, k + 1:] @ w, w); V[:, k + 1:] -= 2.0 * np.outer(V[:, k + 1:] @ w, w)
    B = np.triu(np.tril(B, 1))
    return U, B, V


A = rng.standard_normal((9, 5))
Ub, Bb, Vb = golub_kahan_bidiagonalize(A)
assert np.allclose(Ub @ Bb @ Vb.T, A)
assert np.allclose(np.linalg.svd(Bb[:5], compute_uv=False), np.linalg.svd(A, compute_uv=False))
print("Golub–Kahan: A = U B Vᵀ with B bidiagonal and the same singular values ✓")

# %% [markdown]
# ## Part 5: Iterative methods on the 2D Poisson equation (§3.9–3.12)
#
# −u_xx − u_yy = f on the unit square, zero Dirichlet data, N×N interior grid,
# h = 1/(N+1). The 5-point stencil gives A = (I⊗L + L⊗I)/h² with L = tridiag(−1, 2, −1):
# an N²×N² SPD matrix with ≤ 5 nonzeros per row (CSR). Eigenvalues (eq. 3.24):
# λ_jk = (4/h²)(sin²(jπh/2) + sin²(kπh/2)), so κ = cot²(πh/2) ≈ 4/(π²h²).


# %%
def laplacian_2d(N):
    """Sparse 5-point Laplacian on an N×N interior grid of the unit square (eq. 3.23), and h."""
    h = 1.0 / (N + 1)
    L = sps.diags([-np.ones(N - 1), 2 * np.ones(N), -np.ones(N - 1)], [-1, 0, 1], format="csr")
    I = sps.identity(N, format="csr")
    return ((sps.kron(I, L) + sps.kron(L, I)) / h**2).tocsr(), h


def laplacian_condition(N):
    h = 1.0 / (N + 1)
    return (np.cos(np.pi * h / 2) / np.sin(np.pi * h / 2)) ** 2


A8, h8 = laplacian_2d(8)
jj, kk = np.meshgrid(np.arange(1, 9), np.arange(1, 9))
lam_formula = np.sort(((4 / h8**2) * (np.sin(jj * np.pi * h8 / 2) ** 2 + np.sin(kk * np.pi * h8 / 2) ** 2)).ravel())
assert np.allclose(lam_formula, np.linalg.eigvalsh(A8.toarray()))
assert np.isclose(lam_formula[-1] / lam_formula[0], laplacian_condition(8))
print(f"2D Laplacian N=8: eigenvalue formula (3.24) matches eigvalsh ✓; CSR arrays: data {A8.data.size}, indices {A8.indices.size}, indptr {A8.indptr.size}")


def conjugate_gradient(A, b, tol=1e-8, maxiter=10000, x0=None, callback=None):
    """CG for SPD A (Algorithm 3.3). Returns x, iterations, and relative residual history ‖r_k‖/‖b‖."""
    x = np.zeros(b.shape[0]) if x0 is None else x0.astype(float).copy()
    r = b - A @ x; p = r.copy(); rr = r @ r; nb = np.linalg.norm(b); hist = [np.sqrt(rr) / nb]
    for k in range(maxiter):
        Ap = A @ p
        alpha = rr / (p @ Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rr_new = r @ r
        hist.append(np.sqrt(rr_new) / nb)
        if callback is not None:
            callback(x)
        if hist[-1] < tol:
            return x, k + 1, np.array(hist)
        p = r + (rr_new / rr) * p
        rr = rr_new
    return x, maxiter, np.array(hist)


# Example 4.5: CG by hand on [[4,1],[1,3]] x = [1,2]: x₁ = (1/4, 1/2), x₂ = (1/11, 7/11), r₂ = 0.
A45 = np.array([[4.0, 1.0], [1.0, 3.0]]); b45 = np.array([1.0, 2.0])
iterates = []
x_cg, its, hist = conjugate_gradient(A45, b45, tol=1e-15, callback=lambda x: iterates.append(x.copy()))
assert its == 2 and np.allclose(iterates[0], [0.25, 0.5]) and np.allclose(iterates[1], [1 / 11, 7 / 11])
assert np.allclose(A45 @ x_cg, b45)
r0 = b45; p0 = r0; Ap0 = A45 @ p0; alpha0 = (r0 @ r0) / (p0 @ Ap0); r1 = r0 - alpha0 * Ap0
beta0 = (r1 @ r1) / (r0 @ r0); p1 = r1 + beta0 * p0
assert np.allclose(Ap0, [6, 7]) and alpha0 == 0.25 and np.allclose(r1, [-0.5, 0.25]) and beta0 == 0.0625
assert np.allclose(p1, [-0.4375, 0.375]) and np.isclose(r1 @ r0, 0) and np.isclose(p1 @ A45 @ p0, 0)
alpha1 = (r1 @ r1) / (p1 @ A45 @ p1)
assert np.isclose(alpha1, 4 / 11) and np.allclose(A45 @ p1, [-1.375, 0.6875])
print(f"Example 4.5: α₀={alpha0}, x₁={iterates[0]}, r₁={r1}, β₀={beta0}, p₁={p1}, α₁={alpha1:.6f}, x₂={iterates[1]} = (1/11, 7/11); r₁ ⟂ r₀ and p₁ A-conjugate to p₀ ✓")

# %% [markdown]
# ### Jacobi vs Gauss–Seidel vs CG on the N = 20 grid (§3.9, Theorem 3.13)
#
# Splitting A = M − N: Jacobi M = D, Gauss–Seidel M = D + L. Convergence iff
# ρ(M⁻¹N) < 1. For the 5-point Laplacian ρ_J = cos(πh) and ρ_GS = cos²(πh).

# %%
N = 20
A, h = laplacian_2d(N); n = N * N
b = rng.standard_normal(n); nb = np.linalg.norm(b)
x_star = spla.spsolve(A.tocsc(), b)
D = A.diagonal(); DL = sps.tril(A, 0, format="csr"); Uo = sps.triu(A, 1, format="csr")

def stationary(solve_M, maxiter=5000, tol=1e-8):
    """x_{k+1} = x_k + M⁻¹ r_k (eq. 3.21). solve_M(r) returns M⁻¹r."""
    x = np.zeros(n); hist = []
    for k in range(maxiter):
        r = b - A @ x; hist.append(np.linalg.norm(r) / nb)
        if hist[-1] < tol:
            break
        x = x + solve_M(r)
    return x, np.array(hist)

t0 = time.time(); x_j, hist_j = stationary(lambda r: r / D); t_j = time.time() - t0
t0 = time.time(); x_gs, hist_gs = stationary(lambda r: spla.spsolve_triangular(DL, r, lower=True)); t_gs = time.time() - t0
t0 = time.time(); x_c, it_c, hist_c = conjugate_gradient(A, b); t_c = time.time() - t0
for xx in (x_j, x_gs, x_c):
    assert np.linalg.norm(xx - x_star) / np.linalg.norm(x_star) < 1e-5
rho_J = np.cos(np.pi * h); rho_GS = rho_J**2
BJ = np.eye(n) - A.toarray() / D[:, None]
assert np.isclose(np.abs(np.linalg.eigvals(BJ)).max(), rho_J)
BGS = np.linalg.solve(DL.toarray(), -Uo.toarray())
assert np.isclose(np.abs(np.linalg.eigvals(BGS)).max(), rho_GS, atol=1e-10)
# Observed convergence rates (per iteration, late in the run) match ρ.
rate_j = (hist_j[-1] / hist_j[-201]) ** (1 / 200); rate_gs = (hist_gs[-1] / hist_gs[-101]) ** (1 / 100)
print(f"N=20 (n=400, κ={laplacian_condition(N):.0f}): Jacobi {len(hist_j)-1} its ({t_j:.2f}s), Gauss–Seidel {len(hist_gs)-1} its ({t_gs:.2f}s), CG {it_c} its ({t_c:.3f}s)")
print(f"   ρ_J = cos πh = {rho_J:.5f} (observed rate {rate_j:.5f}); ρ_GS = cos² πh = {rho_GS:.5f} (observed {rate_gs:.5f})")
assert abs(rate_j - rho_J) < 2e-3 and abs(rate_gs - rho_GS) < 4e-3
assert len(hist_j) > 2 * len(hist_gs) * 0.9 and it_c < len(hist_gs) / 5

fig, ax = plt.subplots(figsize=(6.5, 4))
ax.semilogy(hist_j, label=f"Jacobi ({len(hist_j)-1} its)")
ax.semilogy(hist_gs, label=f"Gauss–Seidel ({len(hist_gs)-1} its)")
ax.semilogy(hist_c, label=f"CG ({it_c} its)")
ax.set_xlabel("iteration"); ax.set_ylabel("‖r_k‖ / ‖b‖"); ax.legend()
ax.set_title(f"2D Poisson, 20×20 grid (κ ≈ {laplacian_condition(N):.0f}): stationary methods vs CG")
fig.tight_layout(); fig.savefig(FIG / "stationary_vs_cg.png", dpi=120); plt.close(fig)

# %% [markdown]
# ### CG on grids 20×20 … 200×200: residual curves, iterations ∝ √κ, timing vs dense solve
#
# Theorem 3.16: ‖e_k‖_A ≤ 2((√κ−1)/(√κ+1))^k ‖e_0‖_A, so iterations to a fixed
# tolerance grow like √κ ∝ N. We check the A-norm error against the bound directly.

# %%
grids = [20, 40, 60, 80, 100, 150, 200]
results = []
curves = {}
for N in grids:
    A, h = laplacian_2d(N); n = N * N
    b = rng.standard_normal(n)
    t0 = time.time(); x_c, it_c, hist_c = conjugate_gradient(A, b); t_c = time.time() - t0
    t_dense = np.nan; x_ref = None
    if N <= 60:
        t0 = time.time(); Ad = A.toarray(); x_ref = np.linalg.solve(Ad, b); t_dense = time.time() - t0
        assert np.linalg.norm(x_c - x_ref) / np.linalg.norm(x_ref) < 1e-6
    kap = laplacian_condition(N)
    results.append((N, n, kap, it_c, t_c, t_dense)); curves[N] = hist_c
    print(f"N={N:3d} n={n:6d} κ={kap:8.0f} √κ={np.sqrt(kap):6.1f}  CG its={it_c:4d} ({it_c/np.sqrt(kap):.2f}√κ)  CG {t_c:.3f}s  dense {'—' if np.isnan(t_dense) else f'{t_dense:.3f}s'}")
res = np.array(results)
ratios = res[:, 3] / np.sqrt(res[:, 2])
assert ratios.std() / ratios.mean() < 0.05, "iterations/√κ should be roughly constant"
print(f"iterations / √κ = {ratios.mean():.2f} ± {ratios.std():.2f} across all grids: CG cost scales like √κ ∝ N")

# Chebyshev bound check on N = 40 with a known solution (so we can measure the A-norm error).
A, h = laplacian_2d(40); x_true = rng.standard_normal(1600); b = A @ x_true
kap = laplacian_condition(40); rho = (np.sqrt(kap) - 1) / (np.sqrt(kap) + 1)
errs = []
e0 = x_true.copy(); e0A = np.sqrt(e0 @ (A @ e0))
conjugate_gradient(A, b, tol=1e-12, callback=lambda x: errs.append(np.sqrt((x_true - x) @ (A @ (x_true - x))) / e0A))
errs = np.array(errs); ks = np.arange(1, len(errs) + 1)
bound = 2 * rho**ks
assert np.all(errs <= bound * (1 + 1e-8)), "CG violated the Chebyshev bound"
print(f"Chebyshev bound (Thm 3.16) holds at every step on N=40: e.g. k=50: ‖e‖_A/‖e₀‖_A = {errs[49]:.2e} ≤ 2ρ^50 = {bound[49]:.2e}, ρ = (√κ−1)/(√κ+1) = {rho:.4f}")

fig, axs = plt.subplots(1, 3, figsize=(13, 4))
for N in grids:
    axs[0].semilogy(curves[N], label=f"N={N}")
axs[0].set_xlabel("CG iteration"); axs[0].set_ylabel("‖r_k‖/‖b‖"); axs[0].legend(); axs[0].set_title("Residual histories")
axs[1].plot(np.sqrt(res[:, 2]), res[:, 3], "o", label="measured")
axs[1].plot(np.sqrt(res[:, 2]), ratios.mean() * np.sqrt(res[:, 2]), "--", color="gray", label=f"{ratios.mean():.2f}·√κ")
axs[1].set_xlabel("√κ(A)"); axs[1].set_ylabel("iterations to ‖r‖/‖b‖ < 1e-8"); axs[1].legend(); axs[1].set_title("Iterations ∝ √κ")
axs[2].loglog(res[:, 1], res[:, 4], "o-", label="CG (sparse)")
mask = ~np.isnan(res[:, 5])
axs[2].loglog(res[mask, 1], res[mask, 5], "s-", label="dense LU solve")
axs[2].set_xlabel("n = N²"); axs[2].set_ylabel("seconds"); axs[2].legend(); axs[2].set_title("Wall time")
fig.tight_layout(); fig.savefig(FIG / "cg_laplacian.png", dpi=120); plt.close(fig)

# %% [markdown]
# ### Arnoldi and GMRES for a nonsymmetric matrix (§3.11, sketch)
#
# Convection–diffusion: −u'' + β u' on a 1D grid gives a nonsymmetric tridiagonal.
# Arnoldi builds AQ_k = Q_{k+1} H̃_k; GMRES minimizes ‖b − Ax‖ over the Krylov space
# by solving the small (k+1)×k least-squares problem with H̃_k.


# %%
def arnoldi(A, b, k):
    """k steps of Arnoldi (modified Gram–Schmidt): A Q_k = Q_{k+1} H̃_k with H̃_k (k+1)×k upper Hessenberg (eq. 3.27)."""
    n = b.shape[0]; Q = np.zeros((n, k + 1)); H = np.zeros((k + 1, k))
    Q[:, 0] = b / np.linalg.norm(b)
    for j in range(k):
        v = A @ Q[:, j]
        for i in range(j + 1):
            H[i, j] = Q[:, i] @ v; v = v - H[i, j] * Q[:, i]
        H[j + 1, j] = np.linalg.norm(v)
        if H[j + 1, j] < 1e-14:
            return Q[:, :j + 2], H[:j + 2, :j + 1]
        Q[:, j + 1] = v / H[j + 1, j]
    return Q, H


def gmres(A, b, k):
    """GMRES without restarts: x_k = Q_k y minimizing ‖‖b‖e₁ − H̃_k y‖ (eq. 3.28). Returns x_k and residual history."""
    Q, H = arnoldi(A, b, k); kk = H.shape[1]
    beta = np.linalg.norm(b); hist = []
    for j in range(1, kk + 1):
        y, *_ = np.linalg.lstsq(H[:j + 1, :j], beta * np.eye(j + 1)[0], rcond=None)
        hist.append(np.linalg.norm(beta * np.eye(j + 1)[0] - H[:j + 1, :j] @ y))
    x = Q[:, :kk] @ y
    return x, np.array(hist) / beta


m = 200; hh = 1.0 / (m + 1); beta_cd = 30.0
Acd = (sps.diags([-np.ones(m - 1), 2 * np.ones(m), -np.ones(m - 1)], [-1, 0, 1]) / hh**2
       + beta_cd * sps.diags([-np.ones(m - 1), np.ones(m - 1)], [-1, 1]) / (2 * hh)).tocsr()
assert abs(Acd - Acd.T).max() > 0            # genuinely nonsymmetric
bcd = np.ones(m)
Qa, Ha = arnoldi(Acd, bcd, 30)
assert np.allclose(Acd @ Qa[:, :30], Qa @ Ha) and np.allclose(Qa.T @ Qa, np.eye(31))
x_g, hist_g = gmres(Acd, bcd, 120)
x_cd = spla.spsolve(Acd.tocsc(), bcd)
print(f"GMRES on a 200×200 convection–diffusion matrix: relative residual {hist_g[-1]:.1e} after {len(hist_g)} steps; error vs direct solve {np.linalg.norm(x_g - x_cd)/np.linalg.norm(x_cd):.1e}")
assert hist_g[-1] < 1e-8 and np.all(np.diff(hist_g) <= 1e-12)   # GMRES residuals never increase

# %% [markdown]
# ## Part 6: Randomized SVD (§3.13)
#
# Range finder: Y = (AAᵀ)^q A Ω with Gaussian Ω (n × (k+p)), Q = qr(Y), B = QᵀA,
# SVD of the small B. Test matrix: 400×300 with σ_j = 1/j (slow decay, where
# oversampling p and power iterations q both matter). Error ‖A − QQᵀA‖₂ vs the
# optimum σ_{k+1}; the HMT bounds (Theorem 3.17) are checked on the mean over trials.


# %%
def randomized_svd(A, k, p=10, q=0, rng=rng):
    """Rank-k approximation A ≈ U Σ Vᵀ by the randomized range finder with oversampling p and q power iterations
    (Algorithm 3.4; HMT 2011 Algorithms 4.3/4.4 + 5.1). Returns U (m×k), s (k), Vt (k×n), and Q (m×(k+p))."""
    m, n = A.shape
    Omega = rng.standard_normal((n, k + p))
    Y = A @ Omega
    Q, _ = np.linalg.qr(Y)
    for _ in range(q):                                   # power iterations with re-orthonormalization
        Z, _ = np.linalg.qr(A.T @ Q)
        Q, _ = np.linalg.qr(A @ Z)
    B = Q.T @ A
    Ub, s, Vt = np.linalg.svd(B, full_matrices=False)
    return Q @ Ub[:, :k], s[:k], Vt[:k], Q


m, n, k = 400, 300, 10
Ur, _ = np.linalg.qr(rng.standard_normal((m, n))); Vr, _ = np.linalg.qr(rng.standard_normal((n, n)))
sig = 1.0 / np.arange(1, n + 1)
A_test = (Ur * sig) @ Vr.T
sig_exact = np.linalg.svd(A_test, compute_uv=False)
assert np.allclose(sig_exact, sig)
opt = sig[k]                                            # σ_{k+1}: Eckart–Young optimum (LA-9)
ps = [2, 5, 10, 20, 40]; qs = [0, 1, 2]; trials = 10
err = np.zeros((len(qs), len(ps))); errF = np.zeros_like(err)
tailF = np.sqrt(np.sum(sig[k:] ** 2))
for i, q in enumerate(qs):
    for j, p in enumerate(ps):
        e2 = []; ef = []
        for _ in range(trials):
            _, _, _, Q = randomized_svd(A_test, k, p, q)
            E = A_test - Q @ (Q.T @ A_test)
            e2.append(np.linalg.norm(E, 2)); ef.append(np.linalg.norm(E))
        err[i, j] = np.mean(e2); errF[i, j] = np.mean(ef)
        # HMT Theorem 10.5 (Frobenius) and 10.6 / Cor. 10.10 (spectral) bounds in expectation.
        boundF = np.sqrt(1 + k / (p - 1)) * tailF
        tail_q = np.sqrt(np.sum(sig[k:] ** (2 * (2 * q + 1))))
        bound2 = ((1 + np.sqrt(k / (p - 1))) * sig[k] ** (2 * q + 1) + np.e * np.sqrt(k + p) / p * tail_q) ** (1 / (2 * q + 1))
        if q == 0:
            assert errF[i, j] <= boundF * 1.05, (q, p, errF[i, j], boundF)
        assert err[i, j] <= bound2 * 1.05, (q, p, err[i, j], bound2)
        if q == 0:
            bound11 = (1 + 4 * np.sqrt(k + p) / (p - 1) * np.sqrt(min(m, n))) * sig[k]     # HMT Theorem 1.1
            assert err[i, j] <= bound11
print("randomized SVD, k=10, σ_j = 1/j: mean ‖A − QQᵀA‖₂ / σ_{k+1}  (rows q = 0,1,2; columns p = " + ", ".join(map(str, ps)) + ")")
for i, q in enumerate(qs):
    print(f"   q={q}: " + "  ".join(f"{v/opt:5.2f}" for v in err[i]))
assert np.all(np.diff(err, axis=1) < 0.05) and np.all(err[2] < err[0])    # more p and more q both help
assert err[2, -1] / opt < 1.15                                              # q=2, p=40: within 15% of optimal
bound11_p10 = (1 + 4 * np.sqrt(k + 10) / 9 * np.sqrt(300)) * sig[k]
print(f"   HMT Theorem 1.1 bound at p=10, q=0: {bound11_p10/opt:.1f} σ_{k+1} (loose but true: measured {err[0, 2]/opt:.2f} σ_{k+1})")
U_r, s_r, Vt_r, _ = randomized_svd(A_test, k, 10, 2)
assert np.max(np.abs(s_r - sig[:k]) / sig[:k]) < 0.02
print(f"   leading 10 singular values from randomized_svd (p=10, q=2): max relative error {np.max(np.abs(s_r - sig[:k]) / sig[:k]):.1e}")

fig, ax = plt.subplots(figsize=(6.5, 4))
for i, q in enumerate(qs):
    ax.plot(ps, err[i] / opt, "o-", label=f"q = {q} power iterations")
ax.axhline(1, color="gray", ls="--", lw=1, label="optimal: σ_{k+1}")
ax.set_xlabel("oversampling p"); ax.set_ylabel("‖A − QQᵀA‖₂ / σ_{k+1}"); ax.set_xscale("log"); ax.legend()
ax.set_title("Randomized range finder, k=10 on σ_j = 1/j (mean of 10 trials)")
fig.tight_layout(); fig.savefig(FIG / "randomized_svd.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Your turn
#
# 1. **Rayleigh quotient iteration.** Write `rayleigh_quotient_iteration(A, v0)`:
#    μ ← vᵀAv/vᵀv, solve (A − μI)w = v, v ← w/‖w‖. Print |μ_k − λ| at each step
#    and the number of correct digits; it should roughly *triple* per step (§3.7).
# 2. **Preconditioned CG.** Implement PCG with M = diag(A) (Algorithm 3.3′) and
#    run it on the variable-coefficient problem −∇·(a∇u) = f with a = 1 or 10³ in
#    a sub-square. Compare iteration counts with plain CG. Then explain why Jacobi
#    does *nothing* on the constant-coefficient Laplacian (what is D⁻¹A there?).
# 3. **Lanczos.** Implement 30 steps of Lanczos on the 100×100-grid Laplacian
#    (10,000 × 10,000) and compare the extreme Ritz values with λ_min, λ_max from
#    (3.24). Which end converges first, and why (think about eigenvalue spacing)?
# 4. Exercise B9 asks you to add a Hessenberg reduction for *nonsymmetric* A to
#    `householder_tridiagonalize` and run the shifted QR algorithm on it.
#
# Reference implementations follow (write yours before reading them).

# %%
def rayleigh_quotient_iteration(A, v0, steps=6):
    """RQI (Algorithm 3.2): shift = Rayleigh quotient, one inverse-iteration step per update. Returns the shifts."""
    v = v0 / np.linalg.norm(v0); mus = []
    for _ in range(steps):
        mu = v @ A @ v
        mus.append(mu)
        try:
            w = np.linalg.solve(A - mu * np.eye(len(v)), v)
        except np.linalg.LinAlgError:
            break                                        # exactly singular: we have converged
        v = w / np.linalg.norm(w)
    return np.array(mus), v


Qr6, _ = np.linalg.qr(rng.standard_normal((6, 6)))
A6 = Qr6 @ np.diag([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]) @ Qr6.T
v0 = Qr6[:, 3] + 0.3 * rng.standard_normal(6)       # a rough guess at the eigenvector for λ = 4
mus, _ = rayleigh_quotient_iteration(A6, v0, steps=5)
ev6 = np.linalg.eigvalsh(A6); lam_hit = ev6[np.argmin(np.abs(ev6 - mus[-1]))]
errs_rqi = np.abs(mus - lam_hit)
print("RQI |μ_k − λ|:", " ".join(f"{e:.1e}" for e in errs_rqi), " → correct digits:", [int(-np.log10(max(e, 1e-17))) for e in errs_rqi])
# Cubic convergence: e_{k+1} ≈ C e_k³, so the digits roughly triple each step until ε_mach.
converging = errs_rqi[(errs_rqi > 1e-15) & (errs_rqi < 1e-1)]
if len(converging) >= 3:
    for e_prev, e_next in zip(converging[:-1], converging[1:]):
        assert e_next < 5 * e_prev ** 2.5, (e_prev, e_next)          # at least much faster than quadratic
assert errs_rqi[-1] < 1e-12


def variable_coefficient_poisson(N, a_func):
    """Sparse 5-point discretization of −∇·(a∇u) on the unit square with Dirichlet data: the flux across each cell
    face uses a evaluated at the face midpoint (eq. 3.29). Returns A (N²×N², SPD) and h."""
    h = 1.0 / (N + 1); x = np.arange(1, N + 1) * h
    X, Y = np.meshgrid(x, x, indexing="ij")
    aE = a_func(X + h / 2, Y); aW = a_func(X - h / 2, Y); aN = a_func(X, Y + h / 2); aS = a_func(X, Y - h / 2)
    idx = np.arange(N * N).reshape(N, N)
    diag = (aE + aW + aN + aS).ravel() / h**2
    rows = [idx.ravel()]; cols = [idx.ravel()]; vals = [diag]
    rows += [idx[:-1, :].ravel(), idx[1:, :].ravel(), idx[:, :-1].ravel(), idx[:, 1:].ravel()]
    cols += [idx[1:, :].ravel(), idx[:-1, :].ravel(), idx[:, 1:].ravel(), idx[:, :-1].ravel()]
    vals += [-aE[:-1, :].ravel() / h**2, -aW[1:, :].ravel() / h**2, -aN[:, :-1].ravel() / h**2, -aS[:, 1:].ravel() / h**2]
    A = sps.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N * N, N * N))
    return A, h


def preconditioned_cg(A, b, Minv, tol=1e-8, maxiter=10000):
    """PCG (Algorithm 3.3′): CG on the preconditioned system without forming M^{-1/2}. Minv(r) returns M⁻¹r."""
    x = np.zeros(b.shape[0]); r = b.copy(); z = Minv(r); p = z.copy(); rz = r @ z; nb = np.linalg.norm(b)
    for k in range(maxiter):
        Ap = A @ p; alpha = rz / (p @ Ap)
        x = x + alpha * p; r = r - alpha * Ap
        if np.linalg.norm(r) / nb < tol:
            return x, k + 1
        z = Minv(r); rz_new = r @ z
        p = z + (rz_new / rz) * p; rz = rz_new
    return x, maxiter


a_jump = lambda X, Y: np.where((X > 0.25) & (X < 0.75) & (Y > 0.25) & (Y < 0.75), 1e3, 1.0)
Nv = 40
Av, hv = variable_coefficient_poisson(Nv, a_jump)
assert abs(Av - Av.T).max() < 1e-9
assert np.allclose(variable_coefficient_poisson(10, lambda X, Y: np.ones_like(X))[0].toarray(), laplacian_2d(10)[0].toarray())
bv = rng.standard_normal(Nv * Nv)
x_plain, it_plain, _ = conjugate_gradient(Av, bv)
Dv = Av.diagonal()
x_pcg, it_pcg = preconditioned_cg(Av, bv, lambda r: r / Dv)
x_dir = spla.spsolve(Av.tocsc(), bv)
assert np.linalg.norm(x_pcg - x_dir) / np.linalg.norm(x_dir) < 1e-5 and np.linalg.norm(x_plain - x_dir) / np.linalg.norm(x_dir) < 1e-5
kap_v = np.linalg.cond(Av.toarray()); kap_pv = np.linalg.cond((Av.toarray() / Dv[:, None]) / np.sqrt(Dv)[None, :] * np.sqrt(Dv)[:, None])
print(f"variable-coefficient Poisson (a = 1 or 1e3), N={Nv}: κ(A) = {kap_v:.1e}, κ(D^-1/2 A D^-1/2) = {kap_pv:.1e}; CG {it_plain} its, Jacobi-PCG {it_pcg} its")
assert it_pcg < 0.5 * it_plain
# On the constant-coefficient Laplacian D = (4/h²) I, so D⁻¹A = (h²/4) A has the same κ: Jacobi does nothing.
A, h = laplacian_2d(40); b = rng.standard_normal(1600)
_, it_a, _ = conjugate_gradient(A, b); _, it_b = preconditioned_cg(A, b, lambda r: r / A.diagonal())
assert it_a == it_b
print(f"   constant-coefficient Laplacian N=40: CG {it_a} its, Jacobi-PCG {it_b} its (identical: D is a multiple of I)")


def lanczos(A, b, k):
    """k steps of Lanczos (Algorithm 3.5, T&B Alg. 36.1): A Q_k = Q_k T_k + β_k q_{k+1} e_kᵀ, T_k symmetric tridiagonal.
    No reorthogonalization. Returns alphas (diag), betas (off-diag), Q."""
    n = b.shape[0]; Q = np.zeros((n, k + 1)); alphas = np.zeros(k); betas = np.zeros(k)
    Q[:, 0] = b / np.linalg.norm(b); beta = 0.0; q_prev = np.zeros(n)
    for j in range(k):
        v = A @ Q[:, j] - beta * q_prev
        alphas[j] = Q[:, j] @ v
        v = v - alphas[j] * Q[:, j]
        beta = np.linalg.norm(v); betas[j] = beta
        q_prev = Q[:, j]
        if beta < 1e-14:
            break
        Q[:, j + 1] = v / beta
    return alphas, betas[:-1], Q


N = 100
A, h = laplacian_2d(N)                                   # 10,000 × 10,000
lam_min = (8 / h**2) * np.sin(np.pi * h / 2) ** 2; lam_max = (8 / h**2) * np.cos(np.pi * h / 2) ** 2
b = rng.standard_normal(N * N)
al, be, QL = lanczos(A, b, 30)
Tk = np.diag(al) + np.diag(be, 1) + np.diag(be, -1)
assert np.allclose(QL[:, :30].T @ A @ QL[:, :30], Tk, atol=1e-8)       # T_k = Q_kᵀ A Q_k
ritz_hist = [np.linalg.eigvalsh(np.diag(al[:j]) + np.diag(be[:j - 1], 1) + np.diag(be[:j - 1], -1)) for j in range(1, 31)]
ritz = ritz_hist[-1]
print(f"Lanczos, 30 steps on the 10,000×10,000 Laplacian: Ritz extremes [{ritz[0]:.4f}, {ritz[-1]:.4f}] vs exact [{lam_min:.4f}, {lam_max:.4f}]")
print(f"   relative error: λ_max {abs(ritz[-1]-lam_max)/lam_max:.1e} (top of the spectrum is well separated),  λ_min {abs(ritz[0]-lam_min)/lam_min:.1e} (bottom is clustered: slower)")
assert abs(ritz[-1] - lam_max) / lam_max < 1e-4 and ritz[0] >= lam_min * (1 - 1e-10)   # Ritz values lie inside the spectrum
assert abs(ritz[0] - lam_min) / lam_min < 0.5
# Lanczos-estimated condition number vs exact (an inexpensive κ estimate for CG stopping/iteration prediction)
print(f"   κ estimate from Ritz values: {ritz[-1]/ritz[0]:.0f} vs exact {laplacian_condition(N):.0f}")

fig, ax = plt.subplots(figsize=(6.5, 4))
for j, rv in enumerate(ritz_hist, start=1):
    ax.plot(np.full(len(rv), j), rv, ".", color=PAL[0], ms=3)
ax.axhline(lam_min, color="gray", ls="--", lw=1, label="λ_min, λ_max exact"); ax.axhline(lam_max, color="gray", ls="--", lw=1)
ax.set_xlabel("Lanczos step k"); ax.set_ylabel("Ritz values of T_k"); ax.legend(loc="center right")
ax.set_title("Lanczos on the 10⁴×10⁴ Laplacian: Ritz values converge to the extremes first")
fig.tight_layout(); fig.savefig(FIG / "lanczos_ritz.png", dpi=120); plt.close(fig)

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T_START:.1f} s. Figures in {FIG}")
