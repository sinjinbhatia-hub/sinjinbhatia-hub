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
# "Your turn": Rayleigh quotient iteration, preconditioned CG, Lanczos. Part 7 checks
# every number quoted in solutions.md.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import os
# Single-threaded BLAS: CG below is single-threaded, so this keeps the CG-vs-dense timing fair, and it
# avoids multithreaded dot products on mid-sized vectors, which can be 100× slower on shared machines.
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import math
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
assert np.isclose(kappa_sub, 8.1e7, rtol=1e-2)                            # Example 4.1 quotes 8.1e7

# %% [markdown]
# ## Part 2: Conditioning and backward stability (§3.3–3.4)
#
# Example 4.2: A = Q diag(1, 1e-4) Qᵀ with Q the 45° rotation, so κ₂(A) = 10⁴
# exactly. The bound ‖δx‖/‖x‖ ≤ κ(A)‖δb‖/‖b‖ (Theorem 3.5, eq. 3.7) is *attained* when
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
# ### Forward error ≲ κ × backward error (Theorem 3.8, eq. 3.9) for LU with partial pivoting
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
# The numbers quoted in §3.3: backward error in [4.5e-17, 1.5e-15]; forward error 1.5e-15 at κ = 1, 5.5e-6 at κ = 1e12.
assert bwd.max() < 2e-15 and fwd[0] < 2e-15 and 4e-6 < fwd[-1] < 7e-6
print(f"LU: backward error in [{bwd.min():.1e}, {bwd.max():.1e}] for all κ; forward error {fwd[0]:.1e} at κ=1 → {fwd[-1]:.1e} at κ=1e12")

fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(kappas, fwd, "o-", label="forward error ‖x̂−x‖/‖x‖")
ax.loglog(kappas, bwd, "s-", label="backward error ‖b−Ax̂‖/(‖A‖‖x̂‖)")
ax.loglog(kappas, kappas * EPS, "--", color="gray", label="κ · ε_mach")
ax.set_xlabel("κ₂(A)"); ax.set_ylabel("relative error"); ax.legend()
ax.set_title("Solve Ax=b by LU (n=60): forward error ≈ κ × backward error")
fig.tight_layout(); fig.savefig(FIG / "forward_vs_kappa.png", dpi=120); plt.close(fig)


# %% [markdown]
# ### Growth factor of Gaussian elimination with partial pivoting (§3.4)
#
# ρ = max|u_ij| / max|a_ij|. Wilkinson's bound ρ ≤ 2^{n-1} is attained by the
# matrix with 1 on the diagonal, -1 below it, and a last column of 1s. Random
# matrices show growth of order n^{1/2}–n^{2/3}, nowhere near 2^{n-1}.

# %%
def lu_partial_pivoting(A):
    """PA = LU by Gaussian elimination with partial pivoting; returns P, L, U and the growth factor ρ (eq. 3.10)."""
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
assert 3 < rho < 6                     # §3.4 quotes ρ ≈ 4.2 for this matrix

# %% [markdown]
# ## Part 3: Householder reflections and QR (§3.5–3.6)
#
# Example 4.3: x = (2, 1, 2)ᵀ, ‖x‖ = 3, v = x + sign(x₁)‖x‖e₁ = (5, 1, 2)ᵀ,
# vᵀv = 30, Hx = -3e₁. H = I - vvᵀ/15.


# %%
def householder_vector(x):
    """v such that H = I - 2vvᵀ/vᵀv (eq. 3.11) reflects x onto -sign(x₁)‖x‖ e₁ (eq. 3.13: sign chosen to avoid cancellation)."""
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
assert v43 @ x43 == 15.0 and np.allclose(H43[:, 0], -x43 / 3)         # vᵀx = vᵀv/2; first column of H is −x/‖x‖
print("Example 4.3: v =", v43, " vᵀv =", v43 @ v43, " vᵀx =", v43 @ x43, " Hx =", H43 @ x43)
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
# fl(‖x‖) = 1 exactly, so the bad v is (0, 1e-8, 1e-8), H x = (1, −1e-8, −1e-8): residual √2·1e-8 (Example 4.3 quotes 1.4e-8).
assert np.isclose(err_bad, np.sqrt(2) * 1e-8, rtol=1e-6) and err_good < 1e-20
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
    loss["CGS"].append(np.linalg.norm(Qc.T @ Qc - np.eye(nn), 2))
    loss["MGS"].append(np.linalg.norm(Qm.T @ Qm - np.eye(nn), 2))
    loss["Householder"].append(np.linalg.norm(Qh.T @ Qh - np.eye(nn), 2))
    assert np.linalg.norm(Qh @ Rh - Hn) / np.linalg.norm(Hn) < 1e-14   # backward stable: residual ~ ε
conds = np.array(conds)
for k in loss:
    loss[k] = np.array(loss[k])
assert loss["Householder"].max() < 1e-14
assert loss["CGS"][ns == 10][0] > 1 and loss["MGS"][ns == 10][0] < 1e-2 and loss["MGS"][ns == 10][0] > loss["Householder"][ns == 10][0] * 1e8
assert loss["CGS"][ns == 7][0] > 0.5      # notes §2, §3.5: CGS has lost all orthogonality by n = 7
# The numbers quoted in §3.5: n = 7 CGS 0.58; n = 10 (κ = 1.6e13) CGS 3.0, MGS 1.8e-4, Householder 1.4e-15.
assert 0.55 < loss["CGS"][ns == 7][0] < 0.61 and np.isclose(conds[ns == 10][0], 1.6e13, rtol=0.02)
assert 2.9 < loss["CGS"][ns == 10][0] < 3.1 and 1.7e-4 < loss["MGS"][ns == 10][0] < 1.9e-4 and loss["Householder"].max() < 2e-15
print(f"Hilbert n=7: ‖QᵀQ−I‖₂ CGS = {loss['CGS'][ns==7][0]:.2f}; Householder max over n=2..14 = {loss['Householder'].max():.1e}")
print(f"Hilbert n=10 (κ={conds[ns==10][0]:.1e}): ‖QᵀQ−I‖₂ = CGS {loss['CGS'][ns==10][0]:.1e}, MGS {loss['MGS'][ns==10][0]:.1e}, Householder {loss['Householder'][ns==10][0]:.1e}")

fig, ax = plt.subplots(figsize=(6.5, 4))
for name, mk in zip(loss, ("o", "s", "^")):
    ax.semilogy(ns, np.maximum(loss[name], 1e-17), mk + "-", label=name)
ax.semilogy(ns, np.minimum(conds * EPS, 10), "--", color="gray", label="κ(H) · ε_mach")
ax.set_xlabel("n (Hilbert matrix size)"); ax.set_ylabel("‖QᵀQ − I‖₂"); ax.legend()
ax.set_title("Loss of orthogonality: CGS ~ κ²ε, MGS ~ κε, Householder ~ ε")
fig.tight_layout(); fig.savefig(FIG / "orthogonality_loss.png", dpi=120); plt.close(fig)

# %% [markdown]
# ### Least squares: normal equations (κ²) vs QR (κ), Theorem 3.12
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
assert 5e-3 < err_ne < 7e-3 and 6.5e-10 < err_qr < 8.5e-10            # §3.6 quotes 6e-3 and 7.5e-10

# %% [markdown]
# ## Part 4: Eigenvalue algorithms (§3.7–3.8)
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

# First-order root sensitivity (eq. 3.15) at r = 16 and r = 1 with respect to a₁₉ = −210.
dp16 = math.factorial(15) * math.factorial(4)            # |w'(16)| = Π_{k≠16}|16 − k| = 15!·4!
amp16 = 16**19 / dp16                                     # |r^19 / w'(r)|
cond16 = amp16 * 210 / 16; cond1 = 210 / math.factorial(19)
print(f"eq. (3.15): |w'(16)| = {dp16:.2e}, 16^19 = {16**19:.2e}, |r^19/w'(r)| = {amp16:.2e}, κ(root 16) = {cond16:.2e}, "
      f"κ(root 1) = {cond1:.2e}; linear prediction of the move = {amp16 * 2**-23:.0f}; relative coefficient change = {2**-23/210:.1e}")
assert np.isclose(dp16, 3.14e13, rtol=1e-2) and np.isclose(16**19, 7.56e22, rtol=1e-2) and np.isclose(amp16, 2.41e9, rtol=1e-2)
assert np.isclose(cond16, 3.16e10, rtol=1e-2) and np.isclose(cond1, 1.73e-15, rtol=1e-2)
assert round(amp16 * 2**-23) == 287 and np.isclose(2**-23 / 210, 5.7e-10, rtol=1e-2)

# Inverse iteration with a shift (eq. 3.16): eigenvalues 1, 2, 3 and μ = 2.9 (§3.7).
lam3 = np.array([1.0, 2.0, 3.0]); mu3 = 2.9
tr3 = 1 / (lam3 - mu3)
rate_inv = np.sort(np.abs(tr3))[-2] / np.abs(tr3).max()
steps_inv = math.ceil(np.log(1e-10) / np.log(rate_inv)); steps_pow = math.ceil(np.log(1e-10) / np.log(2 / 3))
assert np.allclose(tr3, [-0.526, -1.111, 10.0], atol=1e-3) and np.isclose(rate_inv, 1 / 9) and (steps_inv, steps_pow) == (11, 57)
rng3 = np.random.default_rng(3)                           # separate stream: keeps every later number in the notes fixed
Q3, _ = np.linalg.qr(rng3.standard_normal((3, 3))); A3m = Q3 @ np.diag(lam3) @ Q3.T
lu3 = sla.lu_factor(A3m - mu3 * np.eye(3)); v3 = rng3.standard_normal(3)
for _ in range(steps_inv):
    v3 = sla.lu_solve(lu3, v3); v3 /= np.linalg.norm(v3)
err3 = min(np.linalg.norm(v3 - Q3[:, 2]), np.linalg.norm(v3 + Q3[:, 2]))
print(f"inverse iteration, μ = 2.9: transformed eigenvalues {np.round(tr3, 3)}, rate {rate_inv:.3f}; "
      f"{steps_inv} steps vs {steps_pow} for power iteration to reach 1e-10; eigenvector error after {steps_inv} steps = {err3:.1e}")
assert err3 < 1e-8


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
    """A_{k+1} = R_k Q_k where A_k = Q_k R_k, A_0 = A (eq. 3.17). Returns the iterates A_0..A_steps and the accumulated
    Q⁽ˢᵗᵉᵖˢ⁾ = Q_0 Q_1 … Q_{steps−1} (Theorem 3.15)."""
    Ak = np.array(A, float); Qacc = np.eye(len(Ak)); hist = [Ak.copy()]
    for _ in range(steps):
        Q, R = qr_positive(Ak); Ak = R @ Q; Qacc = Qacc @ Q; hist.append(Ak.copy())
    return hist, Qacc


A44 = np.array([[3.0, 1.0], [1.0, 1.0]])
hist, Qacc = qr_algorithm_unshifted(A44, 4)
assert np.allclose(hist[1], [[3.4, 0.2], [0.2, 0.6]])
assert np.allclose(hist[2], np.array([[99, 1], [1, 17]]) / 29)
# The hand factors: step 1 Q = [[3,−1],[1,3]]/√10, R = [[10,4],[0,2]]/√10; step 2 Q = [[17,−1],[1,17]]/√290, R = [[58,4],[0,10]]/√290.
Qh1, Rh1 = qr_positive(A44); Qh2, Rh2 = qr_positive(hist[1])
assert np.allclose(Qh1 * np.sqrt(10), [[3, -1], [1, 3]]) and np.allclose(Rh1 * np.sqrt(10), [[10, 4], [0, 2]])
assert np.allclose(Qh2 * np.sqrt(290), [[17, -1], [1, 17]]) and np.allclose(Rh2 * np.sqrt(290), [[58, 4], [0, 10]])
lam = np.sort(np.linalg.eigvalsh(A44))[::-1]
ratio = lam[1] / lam[0]
print(f"Example 4.4: A₁ = {hist[1].tolist()}, A₂ = {np.round(hist[2], 5).tolist()}; off-diagonals {[round(float(h[1,0]),5) for h in hist]}")
print(f"   λ₂/λ₁ = {ratio:.4f}; successive off-diagonal ratios: {[round(float(hist[k+1][1,0]/hist[k][1,0]), 4) for k in range(4)]}")
assert np.isclose(ratio, 3 - 2 * np.sqrt(2)) and np.isclose(hist[2][1, 0] / hist[1][1, 0], ratio, atol=2e-3)
# Theorems 3.15–3.16 with Q⁽ᵏ⁾ = Q_0…Q_{k−1}: A_k = (Q⁽ᵏ⁾)ᵀ A Q⁽ᵏ⁾ and A^k = Q⁽ᵏ⁾ R⁽ᵏ⁾, so Q⁽ᵏ⁾e₁ is normalized power iteration.
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
    assert abs(rate - (5 - j) / (6 - j)) < 0.015, (j, rate)        # tight enough to separate 5/6 from 4/5
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
    rotations (eq. 3.18). Applied to full rows/columns for clarity: O(m²) per step here, O(m) with a banded code."""
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
    """Eigenvalue of the trailing 2×2 block of T closer to its last entry (eq. 3.19; T&B eq. 29.10)."""
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
    assert 1.9 <= steps / nn <= 2.2 and err < 1e-14                 # §3.8(c): about two steps per eigenvalue, error ~6e-15

# Example 4.4, continued: the diagonal of A₂ is off by O(off-diagonal²); the unshifted iteration needs 11 steps to
# push the off-diagonal below 1e-8; the Wilkinson shift (eq. 3.19) for this matrix is the exact eigenvalue 2 − √2.
assert np.isclose(99 / 29 - (2 + np.sqrt(2)), -4.2e-4, atol=5e-6)
hist15, _ = qr_algorithm_unshifted(A44, 15)
k_1e8 = next(k for k, hk in enumerate(hist15) if abs(hk[1, 0]) < 1e-8)
assert k_1e8 == 11 and np.isclose(wilkinson_shift(A44), 2 - np.sqrt(2))
print(f"Example 4.4: A₂[0,0] − (2+√2) = {99/29 - (2 + np.sqrt(2)):.1e}; unshifted QR needs {k_1e8} steps for |off-diag| < 1e-8; "
      f"Wilkinson shift = {wilkinson_shift(A44):.4f} = 2 − √2")

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
# an N²×N² SPD matrix with ≤ 5 nonzeros per row (CSR). Eigenvalues (eq. 3.22):
# λ_jk = (4/h²)(sin²(jπh/2) + sin²(kπh/2)), so κ = cot²(πh/2) ≈ 4/(π²h²).


# %%
def laplacian_2d(N):
    """Sparse 5-point Laplacian on an N×N interior grid of the unit square (eq. 3.21), and h."""
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
print(f"2D Laplacian N=8: eigenvalue formula (3.22) matches eigvalsh ✓; CSR arrays: data {A8.data.size}, indices {A8.indices.size}, indptr {A8.indptr.size}")


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
assert np.isclose(alpha1, 4 / 11) and np.allclose(A45 @ p1, [-1.375, 0.6875]) and np.isclose(p1 @ A45 @ p1, 55 / 64)
assert np.allclose(r1 - alpha1 * (A45 @ p1), 0)                         # r₂ = 0 exactly (Theorem 3.21 with n = 2)
print(f"Example 4.5: α₀={alpha0}, x₁={iterates[0]}, r₁={r1}, β₀={beta0}, p₁={p1}, α₁={alpha1:.6f}, x₂={iterates[1]} = (1/11, 7/11); r₁ ⟂ r₀ and p₁ A-conjugate to p₀ ✓")

# %% [markdown]
# ### Jacobi vs Gauss–Seidel vs CG on the N = 20 grid (§3.9, Theorem 3.19)
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
    """x_{k+1} = x_k + M⁻¹ r_k (eq. 3.23). solve_M(r) returns M⁻¹r."""
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
assert abs(rate_j - rho_J) < 5e-6 and abs(rate_gs - rho_GS) < 5e-6          # "match to five digits" (§3.9)
assert abs(len(hist_j) - 1 - 1416) <= 2 and abs(len(hist_gs) - 1 - 703) <= 2  # the counts quoted in §3.9 and Example 4.6
assert len(hist_j) > 2 * len(hist_gs) * 0.9 and it_c < len(hist_gs) / 5

fig, ax = plt.subplots(figsize=(6.5, 4))
ax.semilogy(hist_j, label=f"Jacobi ({len(hist_j)-1} its)")
ax.semilogy(hist_gs, label=f"Gauss–Seidel ({len(hist_gs)-1} its)")
ax.semilogy(hist_c, label=f"CG ({it_c} its)")
ax.set_xlabel("iteration"); ax.set_ylabel("‖r_k‖ / ‖b‖"); ax.legend()
ax.set_title(f"2D Poisson, 20×20 grid (κ ≈ {laplacian_condition(N):.0f}): stationary methods vs CG")
fig.tight_layout(); fig.savefig(FIG / "stationary_vs_cg.png", dpi=120); plt.close(fig)

# Example 4.6: what the theory predicts for this plate before running anything.
kap20 = laplacian_condition(20); rho_cg20 = (np.sqrt(kap20) - 1) / (np.sqrt(kap20) + 1)
pred = {"Jacobi": np.log(1e-8) / np.log(rho_J), "Gauss–Seidel": np.log(1e-8) / np.log(rho_GS),
        "CG (bound 3.30)": np.log(1e-8 / 2) / np.log(rho_cg20)}
meas = {"Jacobi": len(hist_j) - 1, "Gauss–Seidel": len(hist_gs) - 1, "CG (bound 3.30)": it_c}
print(f"Example 4.6: κ = {kap20:.1f}, digits available ≈ {16 - np.log10(kap20):.1f}; per-digit iterations Jacobi "
      f"{-1/np.log10(rho_J):.0f}, Gauss–Seidel {-1/np.log10(rho_GS):.0f}; ρ_CG = {rho_cg20:.4f}")
for key in pred:
    print(f"   {key:16s} predicted ≤ {pred[key]:6.1f} iterations for 1e-8, measured {meas[key]}")
    assert meas[key] <= pred[key]
assert np.isclose(kap20, 178.1, atol=0.05) and round(pred["Jacobi"]) == 1640 and round(pred["Gauss–Seidel"]) == 820
assert math.ceil(pred["CG (bound 3.30)"]) == 128
# (3.30) bounds the A-norm error. For the residual, ‖r_k‖/‖r_0‖ ≤ √κ ‖e_k‖_A/‖e_0‖_A costs an extra factor √κ.
k_res = math.ceil(np.log(1e-8 / (2 * np.sqrt(kap20))) / np.log(rho_cg20))
assert k_res == 145 and it_c == 63
print(f"   CG residual version of (3.30), extra factor √κ = {np.sqrt(kap20):.2f}: ≤ {k_res} iterations; measured {it_c}")
print(f"   Jacobi/CG iteration ratio = {meas['Jacobi'] / it_c:.1f}; h² = {(1 / 21) ** 2:.1e} (order of the discretization error)")
assert 22 < meas["Jacobi"] / it_c < 23

# %% [markdown]
# ### CG on grids 20×20 … 200×200: residual curves, iterations ∝ √κ, timing vs dense solve
#
# Theorem 3.24: ‖e_k‖_A ≤ 2((√κ−1)/(√κ+1))^k ‖e_0‖_A, so iterations to a fixed
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
assert ratios.std() / ratios.mean() < 0.02, "iterations/√κ should be roughly constant"   # §3.10 quotes 4.7√κ ± 1.5%
assert 4.6 < ratios.mean() < 4.85 and abs(res[-1, 3] - 600) <= 6
print(f"iterations / √κ = {ratios.mean():.2f} ± {ratios.std():.2f} across all grids: CG cost scales like √κ ∝ N")
t60 = res[res[:, 0] == 60, 5][0]
print(f"N=200: CG {res[-1, 4]:.2f}s; dense LU extrapolated from the N=60 time by n³: ≈ {t60 * (40000 / 3600)**3 / 60:.0f} min, "
      f"and the dense matrix alone needs {40000**2 * 8 / 1e9:.1f} GB")

# Chebyshev bound check on N = 40 with a known solution (so we can measure the A-norm error).
A, h = laplacian_2d(40); x_true = rng.standard_normal(1600); b = A @ x_true
kap = laplacian_condition(40); rho = (np.sqrt(kap) - 1) / (np.sqrt(kap) + 1)
errs = []
e0 = x_true.copy(); e0A = np.sqrt(e0 @ (A @ e0))
conjugate_gradient(A, b, tol=1e-12, callback=lambda x: errs.append(np.sqrt((x_true - x) @ (A @ (x_true - x))) / e0A))
errs = np.array(errs); ks = np.arange(1, len(errs) + 1)
bound = 2 * rho**ks
assert np.all(errs <= bound * (1 + 1e-8)), "CG violated the Chebyshev bound"
assert 350 < bound[49] / errs[49] < 480                                  # "loose by a factor of about 400 at k = 50"
print(f"Chebyshev bound (Thm 3.24) holds at every step on N=40: e.g. k=50: ‖e‖_A/‖e₀‖_A = {errs[49]:.2e} ≤ 2ρ^50 = {bound[49]:.2e}, ρ = (√κ−1)/(√κ+1) = {rho:.4f}")

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
    """k steps of Arnoldi (modified Gram–Schmidt): A Q_k = Q_{k+1} H̃_k with H̃_k (k+1)×k upper Hessenberg (eq. 3.31)."""
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
    """GMRES without restarts: x_k = Q_k y minimizing ‖‖b‖e₁ − H̃_k y‖ (eq. 3.32). Returns x_k and residual history."""
    Q, H = arnoldi(A, b, k); kk = H.shape[1]
    beta = np.linalg.norm(b); hist = []
    for j in range(1, kk + 1):
        y, *_ = np.linalg.lstsq(H[:j + 1, :j], beta * np.eye(j + 1)[0], rcond=None)
        hist.append(np.linalg.norm(beta * np.eye(j + 1)[0] - H[:j + 1, :j] @ y))
    x = Q[:, :kk] @ y
    return x, np.array(hist) / beta


# (a) Finite termination: on a 60×60 convection–diffusion matrix GMRES must reach the exact solution at step ≤ 60.
m = 60; hh = 1.0 / (m + 1); beta_cd = 30.0
Acd = (sps.diags([-np.ones(m - 1), 2 * np.ones(m), -np.ones(m - 1)], [-1, 0, 1]) / hh**2
       + beta_cd * sps.diags([-np.ones(m - 1), np.ones(m - 1)], [-1, 1]) / (2 * hh)).tocsr()
assert abs(Acd - Acd.T).max() > 0            # genuinely nonsymmetric
bcd = np.ones(m)
Qa, Ha = arnoldi(Acd, bcd, 30)
assert np.allclose(Acd @ Qa[:, :30], Qa @ Ha) and np.allclose(Qa.T @ Qa, np.eye(31))
x_g, hist_g = gmres(Acd, bcd, m)
x_cd = spla.spsolve(Acd.tocsc(), bcd)
assert np.all(np.diff(hist_g) <= 1e-10)          # residuals never increase (nested subspaces)
assert hist_g[-1] < 1e-10 and hist_g[-2] > 1e-6  # ... and this one needs all 60 steps: the spectrum is spread out
assert 3e-2 < hist_g[-2] < 6e-2 and hist_g[-1] < 5e-13                 # §3.11 quotes 4.4e-2 and 1.1e-13
print(f"GMRES on a 60×60 convection–diffusion matrix: residual {hist_g[-2]:.1e} at step 59, {hist_g[-1]:.1e} at step 60 (finite termination); error vs direct {np.linalg.norm(x_g - x_cd)/np.linalg.norm(x_cd):.1e}")
# (b) A clustered spectrum converges fast: A = 2I + B with ‖B‖₂ = 1. The polynomial p(z) = (1 − z/2)^k has p(0) = 1
# and |p(A)| = ‖(−B/2)^k‖ ≤ 2^-k, so GMRES (which minimizes over all such p, eq. 3.32) satisfies ‖r_k‖ ≤ 2^-k ‖b‖.
nB = 300
Bm = rng.standard_normal((nB, nB)); Bm /= np.linalg.norm(Bm, 2)
A2B = 2 * np.eye(nB) + Bm; b2B = rng.standard_normal(nB)
x_g2, hist_g2 = gmres(A2B, b2B, 40)
assert np.all(hist_g2 <= 2.0 ** -np.arange(1, 41) * (1 + 1e-9)) and hist_g2[-1] < 1e-10
print(f"GMRES on 2I + B (‖B‖=1, n=300): ‖r_k‖/‖b‖ at k=10,20,30: {hist_g2[9]:.1e}, {hist_g2[19]:.1e}, {hist_g2[29]:.1e}  ≤ 2^-k = {2.0**-10:.1e}, {2.0**-20:.1e}, {2.0**-30:.1e}")

# %% [markdown]
# ## Part 6: Randomized SVD (§3.13)
#
# Range finder: Y = (AAᵀ)^q A Ω with Gaussian Ω (n × (k+p)), Q = qr(Y), B = QᵀA,
# SVD of the small B. Test matrix: 400×300 with σ_j = 1/j (slow decay, where
# oversampling p and power iterations q both matter). Error ‖A − QQᵀA‖₂ vs the
# optimum σ_{k+1}; the HMT bounds (Theorem 3.25) are checked on the mean over trials.


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
table_313 = np.array([[2.12, 1.77, 1.42, 0.98, 0.59], [1.05, 0.87, 0.64, 0.45, 0.28], [0.94, 0.78, 0.58, 0.40, 0.25]])
assert np.allclose(err / opt, table_313, atol=0.01)                         # the table printed in §3.13
bound11_p10 = (1 + 4 * np.sqrt(k + 10) / 9 * np.sqrt(300)) * sig[k]
print(f"   HMT Theorem 1.1 bound at p=10, q=0: {bound11_p10/opt:.1f} σ_{k+1} (loose but true: measured {err[0, 2]/opt:.2f} σ_{k+1})")
U_r, s_r, Vt_r, _ = randomized_svd(A_test, k, 10, 2)
assert np.max(np.abs(s_r - sig[:k]) / sig[:k]) < 5e-4                     # §3.13 quotes 4.1e-4
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
#    and the number of correct digits; it should roughly *triple* per step (Theorem 3.14).
# 2. **Preconditioned CG.** Implement PCG with M = diag(A) (Algorithm 3.3′, §3.12) and
#    run it on the variable-coefficient problem −∇·(a∇u) = f with a = 1 or 10³ in
#    a sub-square. Compare iteration counts with plain CG. Then explain why Jacobi
#    does *nothing* on the constant-coefficient Laplacian (what is D⁻¹A there?).
# 3. **Lanczos.** Implement 30 steps of Lanczos on the 100×100-grid Laplacian
#    (10,000 × 10,000) and compare the extreme Ritz values with λ_min, λ_max from
#    (3.22). Which end converges first, and why (think about eigenvalue spacing)?
# 4. Exercise B8 asks you to turn `householder_tridiagonalize` into a Hessenberg
#    reduction for *nonsymmetric* A and run the QR algorithm on the result.
#
# Reference implementations follow (write yours before reading them).

# %%
def rayleigh_quotient_iteration(A, v0, steps=6):
    """RQI (Algorithm 3.2, Theorem 3.14): shift = Rayleigh quotient, one inverse-iteration step per update. Returns the shifts."""
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
digits_rqi = [int(-np.log10(max(e, 1e-17))) for e in errs_rqi]
assert digits_rqi[:3] == [0, 2, 8] and errs_rqi[3] < 2e-15                  # §3.7: 0, 2, 8, 15 digits


def variable_coefficient_poisson(N, a_func):
    """Sparse 5-point discretization of −∇·(a∇u) on the unit square with Dirichlet data: the flux across each cell
    face uses a evaluated at the face midpoint (§3.12). Returns A (N²×N², SPD) and h."""
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


a_smooth = lambda X, Y: np.exp(6 * (X + Y))          # conductivity varying by e^12 ≈ 1.6e5 across the square
a_jump = lambda X, Y: np.where((X > 0.25) & (X < 0.75) & (Y > 0.25) & (Y < 0.75), 1e3, 1.0)
Nv = 40
assert np.allclose(variable_coefficient_poisson(10, lambda X, Y: np.ones_like(X))[0].toarray(), laplacian_2d(10)[0].toarray())
for name, a_func in (("a = exp(6(x+y))", a_smooth), ("a = 1 or 1e3 (jump)", a_jump)):
    Av, hv = variable_coefficient_poisson(Nv, a_func)
    assert abs(Av - Av.T).max() < 1e-9 * abs(Av).max()
    bv = rng.standard_normal(Nv * Nv)
    x_plain, it_plain, _ = conjugate_gradient(Av, bv)
    Dv = Av.diagonal()
    x_pcg, it_pcg = preconditioned_cg(Av, bv, lambda r: r / Dv)
    x_dir = spla.spsolve(Av.tocsc(), bv)
    assert np.linalg.norm(x_pcg - x_dir) / np.linalg.norm(x_dir) < 1e-5 and np.linalg.norm(x_plain - x_dir) / np.linalg.norm(x_dir) < 1e-5
    Ad = Av.toarray(); kap_v = np.linalg.cond(Ad)
    kap_pv = np.linalg.cond(Ad / np.sqrt(Dv)[:, None] / np.sqrt(Dv)[None, :])      # κ(D^-1/2 A D^-1/2) = κ(D⁻¹A)
    print(f"variable-coefficient Poisson, {name}, N={Nv}: κ(A) = {kap_v:.1e} → κ(D⁻¹A) = {kap_pv:.1e} (√ ratio {np.sqrt(kap_v/kap_pv):.0f}); CG {it_plain} its → Jacobi-PCG {it_pcg} its (ratio {it_plain/it_pcg:.0f})")
    if a_func is a_smooth:
        assert kap_v / kap_pv > 1e3 and it_plain / it_pcg > 20
        # §3.12 quotes κ 1.3e6 → 3.6e2 and 4650 → 122 iterations.
        assert 1.25e6 < kap_v < 1.35e6 and 3.5e2 < kap_pv < 3.7e2 and abs(it_plain - 4650) <= 50 and abs(it_pcg - 122) <= 2
    else:   # κ barely moves but the preconditioned spectrum is clustered: CG responds to the whole spectrum, not only κ
        assert kap_v / kap_pv < 3 and it_plain / it_pcg > 5
        # §3.12 quotes κ 5.7e5 → 3.2e5 and 1504 → 155 iterations.
        assert 5.6e5 < kap_v < 5.8e5 and 3.1e5 < kap_pv < 3.3e5 and abs(it_plain - 1504) <= 15 and abs(it_pcg - 155) <= 2
        ev_p = np.linalg.eigvalsh(Ad / np.sqrt(Dv)[:, None] / np.sqrt(Dv)[None, :])
        print(f"   jump case: eigenvalues of D^-1/2 A D^-1/2: one outlier {ev_p[0]:.1e}, the other {len(ev_p)-1} in "
              f"[{ev_p[1]:.1e}, {ev_p[-1]:.2f}] (spread {ev_p[-1]/ev_p[1]:.0f})")
        assert ev_p[0] < 1e-5 and ev_p[1] > 5e-3 and 300 < ev_p[-1] / ev_p[1] < 400
# On the constant-coefficient Laplacian D = (4/h²) I, so D⁻¹A = (h²/4) A has the same κ: Jacobi does nothing.
A, h = laplacian_2d(40); b = rng.standard_normal(1600)
_, it_a, _ = conjugate_gradient(A, b); _, it_b = preconditioned_cg(A, b, lambda r: r / A.diagonal())
assert it_a == it_b and abs(it_a - 126) <= 1                # §3.12 quotes 126 either way
print(f"   constant-coefficient Laplacian N=40: CG {it_a} its, Jacobi-PCG {it_b} its (identical: D is a multiple of I)")


def lanczos(A, b, k):
    """k steps of Lanczos (eq. 3.33, T&B Alg. 36.1): A Q_k = Q_k T_k + β_k q_{k+1} e_kᵀ, T_k symmetric tridiagonal.
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
assert np.isclose(lam_min + lam_max, 8 / h**2)           # the spectrum is symmetric about 4/h²: both ends are equally clustered
b = rng.standard_normal(N * N)
K_LZ = 200
al, be, QL = lanczos(A, b, K_LZ)
T30 = np.diag(al[:30]) + np.diag(be[:29], 1) + np.diag(be[:29], -1)
assert np.allclose(QL[:, :30].T @ A @ QL[:, :30], T30, atol=1e-6)      # T_k = Q_kᵀ A Q_k (Lanczos = Arnoldi for symmetric A)
ritz_hist = [np.linalg.eigvalsh(np.diag(al[:j]) + np.diag(be[:j - 1], 1) + np.diag(be[:j - 1], -1)) for j in range(1, K_LZ + 1)]
for k in (30, 60, 100, 200):
    rv = ritz_hist[k - 1]
    print(f"Lanczos k={k:3d} on the 10,000×10,000 Laplacian: Ritz extremes [{rv[0]:9.3f}, {rv[-1]:10.2f}] vs exact [{lam_min:.3f}, {lam_max:.2f}];"
          f" abs errors {rv[0]-lam_min:8.2f} / {lam_max-rv[-1]:7.2f};  rel errors {(rv[0]-lam_min)/lam_min:.1e} / {(lam_max-rv[-1])/lam_max:.1e}")
    assert rv[0] >= lam_min * (1 - 1e-10) and rv[-1] <= lam_max * (1 + 1e-10)   # Ritz values lie inside the spectrum (B5)
r30 = ritz_hist[29]
# At 30 steps: both ends have converged to the same *absolute* accuracy (the spectrum is symmetric), a few hundred
# out of a width of 8/h² ≈ 81,600. That is 0.1–0.3% of λ_max but several hundred percent of λ_min ≈ 19.7.
assert (r30[0] - lam_min) < 400 and (lam_max - r30[-1]) < 400 and 0.25 < (r30[0] - lam_min) / (lam_max - r30[-1]) < 4
assert (lam_max - r30[-1]) / lam_max < 5e-3 and (r30[0] - lam_min) / lam_min > 1
r200 = ritz_hist[-1]
assert abs(r200[-1] - lam_max) / lam_max < 1e-5 and abs(r200[0] - lam_min) / lam_min < 3e-5   # §3.11: 3e-5 relative or better
# The numbers quoted in §3.11 and B5: Ritz extremes 141.9 and 81,425.9 at k = 30 (abs. errors 122 and 162); κ estimate 574 → 4134.
assert abs(r30[0] - 141.9) < 0.1 and abs(r30[-1] - 81425.9) < 0.1 and round(r30[-1] / r30[0]) == 574
assert abs(r200[-1] / r200[0] - laplacian_condition(N)) < 1 and round(laplacian_condition(N)) == 4134
print(f"   κ estimate from 30 Ritz values: {r30[-1]/r30[0]:.0f}; from 200: {r200[-1]/r200[0]:.0f}; exact {laplacian_condition(N):.0f}")
print("   (Lanczos finds the extremes first in ABSOLUTE terms; the relative gap (λ₂−λ₁)/(λₙ−λ₁) ≈ 3.6e-4 is what sets the speed, and it is the same at both ends.)")

fig, ax = plt.subplots(figsize=(6.5, 4))
for j, rv in enumerate(ritz_hist[:100], start=1):
    ax.semilogy(np.full(len(rv), j), rv, ".", color=PAL[0], ms=2.5)
ax.axhline(lam_min, color="gray", ls="--", lw=1, label="λ_min, λ_max exact"); ax.axhline(lam_max, color="gray", ls="--", lw=1)
ax.set_xlabel("Lanczos step k"); ax.set_ylabel("Ritz values of T_k (log scale)"); ax.legend(loc="center right")
ax.set_title("Lanczos on the 10⁴×10⁴ Laplacian: both ends converge at the same absolute rate")
fig.tight_layout(); fig.savefig(FIG / "lanczos_ritz.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: every number quoted in solutions.md
#
# A separate generator (seed 7) keeps Parts 1–6 unchanged. Exercise B8's reference
# `householder_hessenberg` is here too: write yours first.

# %%
rng7 = np.random.default_rng(7)

# A1: spacing of doubles, the first integer that is not a double, 0.1 + 0.2.
assert np.spacing(1024.0) == 2.0**-42 and np.spacing(2.0**-10) == 2.0**-62
assert int(float(2**53 + 1)) == 2**53 and int(float(2**53)) == 2**53       # 2^53 + 1 rounds (tie → even) to 2^53
assert 0.1 + 0.2 != 0.3 and (0.1 + 0.2) - 0.3 == 2.0**-54 == np.spacing(0.3)
print(f"A1: spacing in [2^10, 2^11) = {2.0**-42:.3e}, in [2^-10, 2^-9) = {2.0**-62:.3e}; (0.1+0.2) − 0.3 = {(0.1 + 0.2) - 0.3:.3e} = one ulp of 0.3")
exact_sum = Fraction(0.1) + Fraction(0.2)                  # the exact sum of the two stored doubles
assert exact_sum - Fraction(0.3) == Fraction(0.1 + 0.2) - exact_sum   # a tie: exactly halfway between two doubles
print(f"   stored 0.1 = {Decimal(0.1)}, 0.2 = {Decimal(0.2)}, 0.3 = {Decimal(0.3)}; fl(0.1+0.2) = {Decimal(0.1 + 0.2)}; "
      f"exact sum lies {float(exact_sum - Fraction(0.3)):.4e} above fl(0.3): a tie, rounded to even")

# A2: 1 − cos x at x = 1e-8.
xa = 1e-8
naive, stable = 1 - np.cos(xa), 2 * np.sin(xa / 2) ** 2
kappa_f = xa * np.sin(xa) / stable                       # κ(f; x) = |x f'(x) / f(x)|
kappa_s = (1 + 1) / stable                                # κ_sub of 1 − cos x, eq. (3.4): (|1| + |cos x|)/|1 − cos x| with the true difference
assert naive == 0.0 and np.isclose(stable, 5e-17, rtol=1e-8) and np.isclose(kappa_f, 2.0, rtol=1e-8) and np.isclose(kappa_s, 4e16, rtol=1e-6)
print(f"A2: 1 − cos(1e-8) naive = {naive}, 2 sin²(x/2) = {stable:.4e}; κ(f) = {kappa_f:.4f}; κ of the subtraction ≈ {kappa_s:.1e}")

# A3: κ of [[1, 1], [1, 1.0001]] in the ∞- and 2-norms.
A3x = np.array([[1.0, 1.0], [1.0, 1.0001]]); A3inv = np.linalg.inv(A3x)
k_inf = np.linalg.norm(A3x, np.inf) * np.linalg.norm(A3inv, np.inf); k_2 = np.linalg.cond(A3x)
ev3 = np.linalg.eigvalsh(A3x)
assert np.allclose(A3inv, [[10001, -10000], [-10000, 10000]], rtol=1e-8)
assert np.isclose(k_inf, 40004.0001, rtol=1e-8) and np.isclose(k_2, 40002.0, rtol=1e-5)
assert np.isclose(ev3[0], 4.99988e-5, rtol=1e-5) and np.isclose(ev3[1], 2.00005, rtol=1e-6)
print(f"A3: κ∞ = {k_inf:.4f}, eigenvalues {ev3[0]:.6e}, {ev3[1]:.6f}, κ₂ = {k_2:.1f}, digits left ≈ {16 - np.log10(k_2):.1f}")

# A4: Householder on x = (−3, 4, 0).
xa4 = np.array([-3.0, 4.0, 0.0]); va4 = householder_vector(xa4); Ha4 = householder_matrix(va4)
assert np.allclose(va4, [-8, 4, 0]) and np.isclose(va4 @ va4, 80) and np.allclose(Ha4 @ xa4, [5, 0, 0])
assert np.allclose(Ha4, [[-0.6, 0.8, 0], [0.8, 0.6, 0], [0, 0, 1]])
print("A4: v =", va4, " Hx =", np.round(Ha4 @ xa4, 12))

# A5: one Jacobi and one Gauss–Seidel step on Example 4.5's system; spectral radii.
xJ = b45 / np.diag(A45); xGS = np.linalg.solve(np.tril(A45), b45)
BJ5 = np.eye(2) - A45 / np.diag(A45)[:, None]; BGS5 = np.linalg.solve(np.tril(A45), -np.triu(A45, 1))
rJ5, rGS5 = np.abs(np.linalg.eigvals(BJ5)).max(), np.abs(np.linalg.eigvals(BGS5)).max()
assert np.allclose(xJ, [1 / 4, 2 / 3]) and np.allclose(xGS, [1 / 4, 7 / 12]) and np.allclose(BGS5, [[0, -1 / 4], [0, 1 / 12]])
assert np.isclose(rJ5, 1 / np.sqrt(12)) and np.isclose(rGS5, 1 / 12)
print(f"A5: Jacobi x₁ = {xJ}, GS x₁ = {xGS}; ρ_J = {rJ5:.4f} ({-np.log10(rJ5):.2f} digits/step), ρ_GS = {rGS5:.4f} ({-np.log10(rGS5):.2f} digits/step)")

# A6: power vs shifted inverse iteration for eigenvalues 10, 9, 1.
st_pow = math.ceil(np.log(1e-6) / np.log(0.9)); st_inv = math.ceil(np.log(1e-6) / np.log(0.2 / 0.8))
assert (st_pow, st_inv) == (132, 10)
st_94 = math.ceil(np.log(1e-6) / np.log(0.4 / 0.6))     # μ = 9.4: closest eigenvalue is 9, rate 0.4/0.6
assert st_94 == 35
print(f"A6: power iteration {st_pow} steps, inverse iteration with μ = 9.8: {st_inv} steps (to 1e-6); μ = 9.4 finds λ = 9 at rate 2/3 in {st_94} steps")

# B1: the residual is the backward error; Example 4.2's bad solution x̂ = (2, 0).
xh = np.array([2.0, 0.0]); rb1 = b42 - A42 @ xh
dA1 = np.outer(rb1, xh) / (xh @ xh)
bwd1 = np.linalg.norm(rb1) / (np.linalg.norm(A42, 2) * np.linalg.norm(xh))
fwd1 = np.linalg.norm(xh - x42) / np.linalg.norm(xh)
assert np.allclose((A42 + dA1) @ xh, b42) and np.isclose(np.linalg.norm(dA1, 2), np.linalg.norm(rb1) / np.linalg.norm(xh))
assert np.allclose(rb1, [-1e-4, 1e-4]) and np.isclose(bwd1, 7.0711e-5, rtol=1e-4) and np.isclose(fwd1, kappa42 * bwd1) and np.isclose(fwd1, 1 / np.sqrt(2))
print(f"B1: r = {rb1}, backward error = {bwd1:.4e}, κ × backward = {kappa42 * bwd1:.4f} = forward error ‖x̂ − x‖/‖x̂‖ = {fwd1:.4f}")

# B2(c): the nearest singular matrix to A is at relative distance 1/κ.
U42, s42, V42t = np.linalg.svd(A42)
A42s = A42 - s42[1] * np.outer(U42[:, 1], V42t[1])
assert np.isclose(s42[1], 1e-4) and np.linalg.svd(A42s, compute_uv=False)[1] < 1e-15 and np.allclose(A42s, 0.5)
print(f"B2: σ_min(A) = {s42[1]:.1e}; A − σ₂u₂v₂ᵀ has σ_min = {np.linalg.svd(A42s, compute_uv=False)[1]:.1e} (singular)")

# B3: the growth matrix, n = 4 by hand and n = 60 in LAPACK (np.linalg.solve).
P4, L4, U4, rho4 = lu_partial_pivoting(wilkinson_growth_matrix(4))
assert np.allclose(P4, np.eye(4)) and np.allclose(U4, [[1, 0, 0, 1], [0, 1, 0, 2], [0, 0, 1, 4], [0, 0, 0, 8]]) and rho4 == 8
assert np.allclose(L4, np.eye(4) - np.tril(np.ones((4, 4)), -1))
W60 = wilkinson_growth_matrix(60); xt60 = rng7.standard_normal(60); b60 = W60 @ xt60
x60 = np.linalg.solve(W60, b60); Q60, R60 = householder_qr(W60); x60q = sla.solve_triangular(R60, Q60.T @ b60)
e_lu = np.linalg.norm(x60 - xt60) / np.linalg.norm(xt60); e_qr = np.linalg.norm(x60q - xt60) / np.linalg.norm(xt60)
bw_lu = np.linalg.norm(b60 - W60 @ x60) / (np.linalg.norm(W60, 2) * np.linalg.norm(x60))
print(f"B3: n=60 growth matrix, κ = {np.linalg.cond(W60):.1f}: GEPP forward error {e_lu:.2f}, backward error {bw_lu:.1e}; Householder QR {e_qr:.1e}")
assert np.linalg.cond(W60) < 30 and e_lu > 0.05 and bw_lu > 1e-3 and e_qr < 1e-13
assert 26.7 < np.linalg.cond(W60) < 26.9 and 0.2 < e_lu < 0.32 and 1e-2 < bw_lu < 2e-2 and e_qr < 1e-14   # as quoted in B3(c)

# B4: Läuchli's matrix, δ = 1e-8.
dl = 1e-8; AL = np.array([[1.0, 1.0], [dl, 0.0], [0.0, dl]]); bL = np.array([2.0, dl, dl])
GL = AL.T @ AL
assert np.array_equal(GL, np.ones((2, 2)))                # 1 + δ² rounds to 1: the computed AᵀA is exactly singular
try:
    np.linalg.solve(GL, AL.T @ bL); ne_failed = False
except np.linalg.LinAlgError:
    ne_failed = True
QL4, RL4 = householder_qr(AL); xL = sla.solve_triangular(RL4[:2], (QL4.T @ bL)[:2])
assert ne_failed and np.allclose(xL, [1, 1], rtol=1e-14, atol=0) and np.isclose(np.linalg.cond(AL), np.sqrt(2 + dl**2) / dl)
print(f"B4: fl(AᵀA) = {GL.tolist()} (singular, solve raises LinAlgError); κ(A) = {np.linalg.cond(AL):.3e}; Householder QR x = {xL}")

# B5: Ritz values from the Lanczos run in Your-turn item 3 (T_k is a leading principal submatrix of T_{k+1}).
tmax = np.array([rv[-1] for rv in ritz_hist]); tmin = np.array([rv[0] for rv in ritz_hist])
assert np.all(np.diff(tmax) >= -1e-9 * lam_max) and np.all(np.diff(tmin) <= 1e-9 * lam_max)
print(f"B5: largest Ritz value nondecreasing and smallest nonincreasing over k = 1..{K_LZ} ✓; κ estimates {r30[-1]/r30[0]:.0f} (k=30) ≤ {laplacian_condition(N):.0f}")

# B6: a shifted QR step with an exact eigenvalue as shift deflates at once; the Rayleigh shift stalls on [[0,1],[1,0]].
mu6 = 2 - np.sqrt(2); Q6s, R6s = qr_positive(A44 - mu6 * np.eye(2)); A1s = R6s @ Q6s + mu6 * np.eye(2)
assert abs(R6s[1, 1]) < 1e-14 and np.isclose(R6s[0, 0], np.sqrt(4 + 2 * np.sqrt(2))) and np.isclose(R6s[0, 1], 2 * np.sqrt(2) / np.sqrt(4 + 2 * np.sqrt(2)))
assert np.allclose(A1s, np.diag([2 + np.sqrt(2), 2 - np.sqrt(2)]), atol=1e-14)
Asw = np.array([[0.0, 1.0], [1.0, 0.0]]); Qsw, Rsw = qr_positive(Asw)
assert np.allclose(Rsw, np.eye(2)) and np.allclose(Rsw @ Qsw, Asw) and wilkinson_shift(Asw) == -1.0
Qw, Rw = qr_positive(Asw + np.eye(2)); A1w = Rw @ Qw - np.eye(2)          # one step with the Wilkinson shift μ = −1
assert np.allclose(A1w, np.diag([1.0, -1.0])) and abs(Rw[1, 1]) < 1e-15
print(f"B6: R = {np.round(R6s, 4).tolist()}, A₁ = {np.round(A1s, 6).tolist()}; [[0,1],[1,0]] is a fixed point of unshifted QR; its Wilkinson shift is −1")

# B7: steepest descent on Example 4.5's system: two steps, the zigzag, and the Kantorovich factor.
xs7 = [np.zeros(2)]
for _ in range(2):
    r7 = b45 - A45 @ xs7[-1]; xs7.append(xs7[-1] + (r7 @ r7) / (r7 @ A45 @ r7) * r7)
x_ex = np.array([1 / 11, 7 / 11]); eA = [np.sqrt((x_ex - x) @ A45 @ (x_ex - x)) for x in xs7]
kap45 = np.linalg.cond(A45)
assert np.allclose(xs7[1], [1 / 4, 1 / 2]) and np.allclose(xs7[2], [1 / 12, 7 / 12]) and np.allclose(b45 - A45 @ xs7[2], b45 / 12)
assert np.isclose(eA[0] ** 2, 15 / 11) and np.isclose(eA[1] ** 2, 5 / 44) and np.isclose(eA[2] ** 2, 5 / 528)
assert np.isclose(eA[1] / eA[0], 1 / np.sqrt(12)) and np.isclose(eA[2] / eA[1], 1 / np.sqrt(12))
assert np.isclose(kap45, (7 + np.sqrt(5)) / (7 - np.sqrt(5))) and np.isclose((kap45 - 1) / (kap45 + 1), np.sqrt(5) / 7)
print(f"B7: steepest descent x₂ = {xs7[2]}, ‖e_k‖_A ratios {eA[1]/eA[0]:.4f}, {eA[2]/eA[1]:.4f} ≤ (κ−1)/(κ+1) = {(kap45-1)/(kap45+1):.4f} (κ = {kap45:.4f})")


# B8: Hessenberg reduction for nonsymmetric A, then the unshifted QR algorithm keeps the Hessenberg form (Theorem 3.17).
def householder_hessenberg(A):
    """A = Q H Qᵀ with H upper Hessenberg, by n−2 two-sided Householder reflections (§3.8(a); T&B Alg. 26.1)."""
    H = np.array(A, float); n = H.shape[0]; Q = np.eye(n)
    for k in range(n - 2):
        v = householder_vector(H[k + 1:, k])
        if v @ v == 0:
            continue
        v = v / np.linalg.norm(v)
        H[k + 1:, k:] -= 2.0 * np.outer(v, v @ H[k + 1:, k:])     # left: rows k+1..n-1
        H[:, k + 1:] -= 2.0 * np.outer(H[:, k + 1:] @ v, v)       # right: columns k+1..n-1
        Q[:, k + 1:] -= 2.0 * np.outer(Q[:, k + 1:] @ v, v)
    return np.triu(H, -1), Q


X8 = rng7.standard_normal((6, 6)); A8n = X8 @ np.diag([6.0, 5.0, 4.0, 3.0, 2.0, 1.0]) @ np.linalg.inv(X8)
H8, Q8 = householder_hessenberg(A8n)
rec8 = np.linalg.norm(Q8 @ H8 @ Q8.T - A8n) / np.linalg.norm(A8n)
hist8, _ = qr_algorithm_unshifted(H8, 150)
below = max(np.abs(np.tril(h, -2)).max() for h in hist8) / np.abs(A8n).max()
diag_err = np.abs(np.diag(hist8[-1]) - np.array([6.0, 5, 4, 3, 2, 1])).max()
assert rec8 < 1e-14 and np.allclose(Q8.T @ Q8, np.eye(6)) and below < 1e-12 and diag_err < 1e-12   # B8(b): 2.7e-13 < (5/6)^150
print(f"B8: ‖QHQᵀ − A‖/‖A‖ = {rec8:.1e}; max entry below the subdiagonal over 150 QR steps = {below:.1e} (relative); "
      f"diag(A₁₅₀) − (6,…,1) max = {diag_err:.1e}; slowest rate (5/6)^150 = {(5 / 6) ** 150:.1e}")

# C1: an SPD matrix on which Jacobi diverges and Gauss–Seidel converges.
aC = 0.8; AC = (1 - aC) * np.eye(3) + aC * np.ones((3, 3))
rJC = np.abs(np.linalg.eigvals(np.eye(3) - AC / np.diag(AC)[:, None])).max()
rGSC = np.abs(np.linalg.eigvals(np.linalg.solve(np.tril(AC), -np.triu(AC, 1)))).max()
assert np.allclose(np.linalg.eigvalsh(AC), [0.2, 0.2, 2.6]) and np.isclose(rJC, 1.6) and np.isclose(rGSC, 0.7155, atol=1e-4)
print(f"C1: a = 0.8: eig(A) = {np.round(np.linalg.eigvalsh(AC), 6)}, ρ(B_J) = {rJC:.4f} (diverges), ρ(B_GS) = {rGSC:.4f} (converges)")

# C2: exact Chebyshev factor vs the bound (3.30), κ = 178.1 (Example 4.6).
gam = (kap20 + 1) / (kap20 - 1)
k_cheb = math.ceil(np.arccosh(1e8) / np.arccosh(gam)); k_bound = math.ceil(np.log(1e-8 / 2) / np.log(rho_cg20))
assert np.isclose(np.arccosh(gam), np.log(1 / rho_cg20)) and k_cheb == k_bound == 128
print(f"C2: arccosh γ = {np.arccosh(gam):.5f} = ln((√κ+1)/(√κ−1)); 1/T_k(γ) ≤ 1e-8 needs k = {k_cheb}, the bound 2ρ^k needs {k_bound}")

# C3: power iterations without re-orthonormalization (σ_j = 1/j test matrix of Part 6, k = 10, p = 10).
c3 = {}
for q in (0, 4, 8, 12):
    plain, reorth = [], []
    for _ in range(5):
        Om = rng7.standard_normal((A_test.shape[1], 20))
        Y = A_test @ Om
        for _ in range(q):
            Y = A_test @ (A_test.T @ Y)                        # (AAᵀ)^q AΩ formed directly
        Qp, _ = np.linalg.qr(Y)
        plain.append(np.linalg.norm(A_test - Qp @ (Qp.T @ A_test), 2) / opt)
        Qo, _ = np.linalg.qr(A_test @ Om)
        for _ in range(q):
            Z, _ = np.linalg.qr(A_test.T @ Qo); Qo, _ = np.linalg.qr(A_test @ Z)
        reorth.append(np.linalg.norm(A_test - Qo @ (Qo.T @ A_test), 2) / opt)
    c3[q] = (np.mean(plain), np.mean(reorth))
print("C3: mean ‖A − QQᵀA‖₂/σ₁₁ (no re-orthonormalization, with): " + "; ".join(f"q={q}: {a:.2f}, {b:.2f}" for q, (a, b) in c3.items()))
assert c3[12][0] > 1.5 * c3[12][1] and c3[8][0] > c3[4][0] and abs(c3[4][0] - c3[4][1]) < 1e-3
print(f"   (σ₁/σ₂₀)^(2q+1) = 20^(2q+1): q=4 → {20.0**9:.0e}, q=8 → {20.0**17:.0e} vs 1/ε_mach = {1/EPS:.0e}")
j8, j12 = (1 / EPS) ** (1 / 17), (1 / EPS) ** (1 / 25)    # directions j with j^(2q+1) < 1/ε survive (σ_j = 1/j)
assert 8.6 < j8 < 8.7 and 4.3 < j12 < 4.4 and 1.0 < c3[8][0] < 1.25 and 2.0 < c3[12][0] < 2.3
print(f"   surviving directions (1/ε)^(1/(2q+1)): q=8 → {j8:.2f}, q=12 → {j12:.2f}; predicted error ≈ σ₉..σ₁₀ (ratio {11/10:.2f}–{11/9:.2f}) and σ₅ (ratio {11/5:.2f})")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T_START:.1f} s. Figures in {FIG}")
