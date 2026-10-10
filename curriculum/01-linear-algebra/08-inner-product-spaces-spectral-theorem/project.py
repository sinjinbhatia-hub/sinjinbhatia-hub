# %% [markdown]
# # LA-8 Project: Inner Products, the Spectral Theorem, and Positive Definiteness
#
# This script builds the computational side of LA-8 and checks every number in
# the notes. Parts: (1) inner products on vectors, functions, and matrices;
# (2) Legendre polynomials by Gram–Schmidt in function space; (3) Jacobi's
# eigenvalue algorithm for symmetric matrices; (4) the spectral decomposition
# and functions of a matrix; (5) a quadratic-form visualizer with principal
# axes; (6) the Rayleigh quotient on the unit circle; (7) three tests for
# positive definiteness; (8) Pauli matrices, unitaries, and the Bloch sphere;
# (9) the 2D stress tensor and Mohr's circle. "Your turn" at the end.
#
# Convention (notes box): complex inner products conjugate the FIRST slot,
# <u, v> = u^† v. Eigenvalues are ordered DESCENDING in the notes; np.linalg.eigh
# returns them ascending, so `eigh_desc` below flips them.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.linalg as sla
import scipy.special
import sympy as sp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)


def eigh_desc(A):
    """Eigenvalues (descending, as in the notes) and orthonormal eigenvectors of a Hermitian A."""
    w, V = np.linalg.eigh(A)
    return w[::-1], V[:, ::-1]


def inner(u, v):
    """Standard inner product <u, v> = u^† v, eq. (3.1): conjugate the FIRST slot."""
    return np.vdot(u, v)  # np.vdot conjugates its first argument


# %% [markdown]
# ## Part 1: Inner products on vectors, functions, and matrices
#
# (3.1) <u, v> = u^† v on C^n.  (3.2) <f, g> = ∫ conj(f) g dx.  (3.3) <A, B>_F = tr(A^† B).
# Theorem 3.4 (Cauchy–Schwarz) and Theorem 3.6 (coordinates are inner products; Parseval).


# %%
def frobenius(A, B):
    """<A, B>_F = tr(A^† B), eq. (3.3)."""
    return np.trace(A.conj().T @ B)


def inner_fn(f, g, a, b, n=64):
    """<f, g> = ∫_a^b conj(f) g dx, eq. (3.2), by n-point Gauss–Legendre quadrature."""
    xs, ws = np.polynomial.legendre.leggauss(n)
    x = 0.5 * (b - a) * xs + 0.5 * (b + a)
    return 0.5 * (b - a) * np.sum(ws * np.conj(f(x)) * g(x))


# Exercise A1: <A, B>_F two ways.
A1 = np.array([[1, 2], [0, -1.]]); B1 = np.array([[3, 0], [1, 1.]])
assert np.isclose(frobenius(A1, B1), 2.0) and np.isclose(np.sum(A1 * B1), 2.0)
assert np.isclose(np.linalg.norm(A1, "fro"), np.sqrt(6)) and np.isclose(np.linalg.norm(B1, "fro"), np.sqrt(11))
print(f"A1: <A,B>_F = {frobenius(A1, B1):.0f}, cos(angle) = {2 / np.sqrt(66):.4f}")

# Fourier's orthogonality (§1) and Exercise A2.
assert abs(inner_fn(np.sin, lambda x: np.sin(2 * x), -np.pi, np.pi)) < 1e-12
assert np.isclose(inner_fn(np.sin, np.sin, -np.pi, np.pi), np.pi)
ip_x_x2 = inner_fn(lambda x: x, lambda x: x**2, 0, 1)
assert np.isclose(ip_x_x2, 1 / 4)
cos_a2 = ip_x_x2 / np.sqrt(inner_fn(lambda x: x, lambda x: x, 0, 1) * inner_fn(lambda x: x**2, lambda x: x**2, 0, 1))
assert np.isclose(cos_a2, np.sqrt(15) / 4)
print(f"A2: <x, x^2> on [0,1] = {ip_x_x2:.4f}, angle = {np.degrees(np.arccos(cos_a2)):.2f} degrees")

# Cauchy–Schwarz on random complex vectors (Theorem 3.4), and conjugate-linearity in the first slot.
for _ in range(200):
    u = rng.standard_normal(5) + 1j * rng.standard_normal(5)
    v = rng.standard_normal(5) + 1j * rng.standard_normal(5)
    assert abs(inner(u, v)) <= np.linalg.norm(u) * np.linalg.norm(v) * (1 + 1e-12)
    c = 0.3 - 1.7j
    assert np.isclose(inner(c * u, v), np.conj(c) * inner(u, v))   # conjugate-linear first slot
    assert np.isclose(inner(u, c * v), c * inner(u, v))            # linear second slot

# Theorem 3.6: coordinates c_i = <q_i, v> and Parseval, with a random unitary basis.
Z = rng.standard_normal((6, 6)) + 1j * rng.standard_normal((6, 6))
U_rand, _ = np.linalg.qr(Z)
assert np.allclose(U_rand.conj().T @ U_rand, np.eye(6))
v = rng.standard_normal(6) + 1j * rng.standard_normal(6)
coeffs = np.array([inner(U_rand[:, i], v) for i in range(6)])        # (3.5): basis vector in FIRST slot
assert np.allclose(U_rand @ coeffs, v)                              # v = Σ <q_i, v> q_i
assert np.isclose(np.sum(np.abs(coeffs) ** 2), np.linalg.norm(v) ** 2)   # (3.6) Parseval
print("Part 1 checks passed: Frobenius, function inner products, Cauchy–Schwarz, Parseval.")

# %% [markdown]
# ## Part 2: Legendre polynomials by Gram–Schmidt in function space
#
# Gram–Schmidt (3.7) on 1, x, x², ... with <f, g> = ∫_{-1}^{1} f g dx. Example 4.1
# predicts u_2 = x² − 1/3 (norm² 8/45) and u_3 = x³ − 3x/5 (norm² 8/175), and after
# rescaling to P_n(1) = 1 the classical Legendre polynomials with ||P_n||² = 2/(2n+1).
#
# **Build it yourself first:** write `gram_schmidt_poly(k)` returning the first k
# orthogonal polynomials as `np.polynomial.Polynomial` objects.


# %%
Poly = np.polynomial.Polynomial
_xs, _ws = np.polynomial.legendre.leggauss(24)   # exact for polynomial integrands of degree <= 47


def ip_poly(p, q):
    """<p, q> = ∫_{-1}^{1} p(x) q(x) dx, exact for polynomials via Gauss–Legendre quadrature."""
    return float(np.sum(_ws * p(_xs) * q(_xs)))


def gram_schmidt_poly(k):
    """First k monic orthogonal polynomials on [-1, 1] from (3.7) applied to 1, x, ..., x^{k-1}."""
    basis = []
    for n in range(k):
        v = Poly([0] * n + [1])          # x^n
        u = v
        for b in basis:
            u = u - (ip_poly(b, v) / ip_poly(b, b)) * b
        basis.append(u)
    return basis


U = gram_schmidt_poly(6)
# Example 4.1 coefficients (Polynomial stores coefficients lowest degree first).
assert np.allclose(U[2].coef, [-1 / 3, 0, 1]), U[2].coef
assert np.allclose(U[3].coef, [0, -3 / 5, 0, 1]), U[3].coef
assert np.allclose([ip_poly(u, u) for u in U[:4]], [2, 2 / 3, 8 / 45, 8 / 175])
# Mutual orthogonality.
for i in range(6):
    for j in range(i):
        assert abs(ip_poly(U[i], U[j])) < 1e-12

# Rescale to P_n(1) = 1 and compare with scipy.special.legendre (poly1d, highest degree first).
P_legendre = [u / u(1.0) for u in U]
for n, Pn in enumerate(P_legendre):
    ref = scipy.special.legendre(n).coeffs[::-1]
    assert np.allclose(Pn.coef, ref, atol=1e-12), (n, Pn.coef, ref)
    assert np.isclose(ip_poly(Pn, Pn), 2 / (2 * n + 1))
assert np.allclose(P_legendre[2].coef, [-0.5, 0, 1.5])        # (3x² − 1)/2
assert np.allclose(P_legendre[3].coef, [0, -1.5, 0, 2.5])     # (5x³ − 3x)/2

# Exact symbolic check of Example 4.1 with SymPy.
x = sp.symbols("x")
sym_ip = lambda f, g: sp.integrate(f * g, (x, -1, 1))
u2 = x**2 - sym_ip(1, x**2) / sym_ip(1, 1) - sym_ip(x, x**2) / sym_ip(x, x) * x
assert sp.simplify(u2 - (x**2 - sp.Rational(1, 3))) == 0 and sym_ip(u2, u2) == sp.Rational(8, 45)
u3 = x**3 - sym_ip(1, x**3) / 2 - sym_ip(x, x**3) / sym_ip(x, x) * x - sym_ip(u2, x**3) / sym_ip(u2, u2) * u2
assert sp.simplify(u3 - (x**3 - sp.Rational(3, 5) * x)) == 0 and sym_ip(u3, u3) == sp.Rational(8, 175)

xx = np.linspace(-1, 1, 400)
fig, ax = plt.subplots(figsize=(6, 4))
for n, Pn in enumerate(P_legendre):
    ax.plot(xx, Pn(xx), label=f"$P_{n}$")
ax.set_title("Legendre polynomials from Gram–Schmidt on $1, x, x^2, \\dots$")
ax.grid(alpha=.3); ax.legend(ncol=3, fontsize=8); ax.set_xlabel("x")
fig.tight_layout(); fig.savefig(FIG / "legendre.png", dpi=120); plt.close(fig)
print("Part 2 checks passed: Gram–Schmidt reproduces P_0..P_5 and ||P_n||^2 = 2/(2n+1).")

# %% [markdown]
# ## Part 3: Jacobi's eigenvalue algorithm
#
# Theorem 3.11 says a symmetric A can be rotated to diagonal form. Jacobi (1846) does
# it one plane at a time. Let G be the identity except for the 2×2 block in rows and
# columns p < q, where it is the rotation [[c, −s], [s, c]], c = cos θ, s = sin θ.
# Then A' = Gᵀ A G is again symmetric with the same eigenvalues (orthogonal similarity).
# Only rows/columns p, q change, and a direct computation of the (p, q) entry gives
#
#     a'_pq = (c² − s²) a_pq + c s (a_qq − a_pp).
#
# Setting a'_pq = 0 and using c² − s² = cos 2θ, 2cs = sin 2θ:
#
#     tan 2θ = 2 a_pq / (a_pp − a_qq).
#
# For numerical stability pick the root with |θ| ≤ π/4: with φ = (a_pp − a_qq)/(2 a_pq) = cot 2θ,
# t = tan θ solves t² + 2φt − 1 = 0, and the small root is t = sign(φ)/(|φ| + √(1 + φ²)).
# Each rotation lowers the off-diagonal Frobenius norm by exactly 2 a_pq² (Exercise B7),
# so always choosing the largest |a_pq| forces off(A) → 0. Accumulating Q ← Q G gives Q
# with A = Q Λ Qᵀ.
#
# **Build it yourself first:** write `jacobi_rotate(A, p, q)` before reading on.


# %%
def jacobi_rotate(A, p, q):
    """One Jacobi rotation: returns (G^T A G, G) where G is the plane rotation in (p, q)
    with tan 2θ = 2 a_pq / (a_pp − a_qq), chosen so that (G^T A G)[p, q] = 0."""
    n = A.shape[0]
    if A[p, q] == 0:
        return A.copy(), np.eye(n)
    phi = (A[p, p] - A[q, q]) / (2 * A[p, q])
    t = np.sign(phi) / (abs(phi) + np.sqrt(1 + phi**2)) if phi != 0 else 1.0
    c = 1 / np.sqrt(1 + t**2); s = t * c
    G = np.eye(n)
    G[p, p] = c; G[q, q] = c; G[p, q] = -s; G[q, p] = s
    return G.T @ A @ G, G


def off_norm2(A):
    """off(A)^2 = sum of squared off-diagonal entries."""
    return np.sum(A**2) - np.sum(np.diag(A) ** 2)


def jacobi_eigh(A, tol=1e-12, max_sweeps=100):
    """Classical Jacobi: eigenvalues (descending) and orthonormal eigenvectors of symmetric A,
    by repeatedly zeroing the largest off-diagonal entry. Also returns the off-norm history."""
    A = np.array(A, float); n = A.shape[0]
    Q = np.eye(n)
    history = [off_norm2(A)]
    for _ in range(max_sweeps * n * n):
        off = np.abs(A - np.diag(np.diag(A)))
        p, q = np.unravel_index(np.argmax(off), off.shape)
        if off[p, q] < tol:
            break
        apq = A[p, q]
        A_new, G = jacobi_rotate(A, p, q)
        assert abs(A_new[p, q]) < 1e-10 * max(1, np.abs(A).max())   # the entry really is zeroed
        # Exercise B7(b),(d): off-norm drops by exactly 2 a_pq^2.
        assert np.isclose(off_norm2(A_new), off_norm2(A) - 2 * apq**2, atol=1e-9 * max(1, np.sum(A**2)))
        A, Q = A_new, Q @ G
        history.append(off_norm2(A))
    w = np.diag(A)
    order = np.argsort(w)[::-1]
    return w[order], Q[:, order], np.array(history)


# Example 4.2.
A42 = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2.]])
w_j, Q_j, hist = jacobi_eigh(A42)
w_ref, Q_ref = eigh_desc(A42)
assert np.allclose(w_j, [4, 1, 1]) and np.allclose(w_j, w_ref)
assert np.allclose(Q_j.T @ Q_j, np.eye(3))                        # Q orthogonal
assert np.allclose(A42 @ Q_j, Q_j @ np.diag(w_j))                 # A Q = Q Λ
assert np.allclose(Q_j @ np.diag(w_j) @ Q_j.T, A42)               # (3.9)
assert np.allclose(np.abs(Q_j[:, 0]), 1 / np.sqrt(3))             # q_1 = ±(1,1,1)/√3
assert np.allclose(Q_j[:, 1:].T @ np.ones(3), 0)                  # λ=1 eigenspace is x+y+z=0
print(f"Example 4.2 via Jacobi: eigenvalues {np.round(w_j, 12)}, {len(hist) - 1} rotations")

# Random symmetric matrices, compared with np.linalg.eigh.
for n in (2, 4, 6, 8):
    for _ in range(5):
        M = rng.standard_normal((n, n)); M = M + M.T
        w_j, Q_j, _ = jacobi_eigh(M)
        w_ref, _ = eigh_desc(M)
        assert np.allclose(w_j, w_ref, atol=1e-9)
        assert np.allclose(Q_j.T @ Q_j, np.eye(n)) and np.allclose(Q_j @ np.diag(w_j) @ Q_j.T, M)

M8 = rng.standard_normal((8, 8)); M8 = M8 + M8.T
_, _, hist8 = jacobi_eigh(M8)
fig, ax = plt.subplots(figsize=(6, 4))
ax.semilogy(hist8, ".-")
ax.set_xlabel("rotation"); ax.set_ylabel("off(A)$^2$")
ax.set_title("Classical Jacobi on a random symmetric 8×8: off-diagonal norm")
ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(FIG / "jacobi_convergence.png", dpi=120); plt.close(fig)
print("Part 3 checks passed: Jacobi agrees with eigh; off-norm drops by 2 a_pq^2 per rotation.")

# %% [markdown]
# ## Part 4: Spectral decomposition and functions of a matrix
#
# Theorem 3.13: A = Σ λ_i q_i q_iᵀ = Σ λ_i P_i, with P_i orthogonal projections summing to I,
# and f(A) = Σ f(λ_i) P_i, eq. (3.11). Example 4.2: A = I + J, A⁻¹ = I − J/4, √A = I + J/3.


# %%
def spectral_projectors(A):
    """Eigenvalues (descending) and the rank-one projectors P_i = q_i q_i^T of Theorem 3.13."""
    w, Q = eigh_desc(A)
    return w, [np.outer(Q[:, i], Q[:, i]) for i in range(len(w))]


def matrix_function(A, f):
    """f(A) = Σ f(λ_i) P_i, eq. (3.11), for symmetric A."""
    w, Ps = spectral_projectors(A)
    return sum(f(wi) * Pi for wi, Pi in zip(w, Ps))


w, Ps = spectral_projectors(A42)
J = np.ones((3, 3)); I3 = np.eye(3)
for i, Pi in enumerate(Ps):
    assert np.allclose(Pi @ Pi, Pi) and np.allclose(Pi, Pi.T)
    for j, Pj in enumerate(Ps):
        if i != j:
            assert np.allclose(Pi @ Pj, 0)
assert np.allclose(sum(Ps), I3)                                   # completeness
assert np.allclose(Ps[0], J / 3)                                  # P_1 = J/3
assert np.allclose(sum(wi * Pi for wi, Pi in zip(w, Ps)), A42)    # (3.10)
assert np.allclose(matrix_function(A42, lambda t: 1 / t), I3 - J / 4)
assert np.allclose(matrix_function(A42, np.sqrt) @ matrix_function(A42, np.sqrt), A42)
assert np.allclose(matrix_function(A42, np.sqrt), I3 + J / 3)
assert np.allclose(matrix_function(A42, np.exp), sla.expm(A42))
assert np.allclose(matrix_function(A42, lambda t: t**5), np.linalg.matrix_power(A42, 5))
print("Part 4 checks passed: A^{-1} = I - J/4, sqrt(A) = I + J/3, exp via (3.11) matches expm.")

# %% [markdown]
# ## Part 5: Quadratic-form visualizer with principal axes
#
# Theorem 3.14: xᵀAx = Σ λ_i y_i² in the eigenbasis; the level set xᵀAx = 1 has semi-axes
# 1/√λ_i along q_i when λ_i > 0. Example 4.3: 5x² + 4xy + 2y² = 1 has λ = 6, 1,
# axes along (2,1)/√5 and (1,−2)/√5, semi-axes 1/√6 ≈ 0.408 and 1.


# %%
def quadratic_form(A, X, Y):
    """xᵀ A x evaluated on a grid, with x = (X, Y)."""
    return A[0, 0] * X**2 + (A[0, 1] + A[1, 0]) * X * Y + A[1, 1] * Y**2


cases = {
    "positive definite (Example 4.3)": np.array([[5, 2], [2, 2.]]),
    "indefinite": np.array([[1, 2], [2, 1.]]),
    "positive semidefinite": np.array([[1, 1], [1, 1.]]),
}
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
g = np.linspace(-1.6, 1.6, 300); X, Y = np.meshgrid(g, g)
for ax, (name, A) in zip(axes, cases.items()):
    w, Q = eigh_desc(A)
    Zq = quadratic_form(A, X, Y)
    ax.contour(X, Y, Zq, levels=np.linspace(-3, 6, 19), cmap="coolwarm", linewidths=.8)
    ax.contour(X, Y, Zq, levels=[1.0], colors="k", linewidths=2)   # the level set xᵀAx = 1
    for i in range(2):
        q = Q[:, i]
        ax.plot([-2 * q[0], 2 * q[0]], [-2 * q[1], 2 * q[1]], "--", color="gray", lw=1)
        if w[i] > 1e-12:
            tip = q / np.sqrt(w[i])
            assert np.isclose(tip @ A @ tip, 1.0)      # the semi-axis tip lies on the level set
            ax.plot([-tip[0], tip[0]], [-tip[1], tip[1]], "o-", color=f"C{i}", lw=2,
                    label=f"$\\lambda_{i+1}={w[i]:.3g}$, semi-axis $1/\\sqrt{{\\lambda}}={1/np.sqrt(w[i]):.3g}$")
        else:
            ax.plot([-2 * q[0], 2 * q[0]], [-2 * q[1], 2 * q[1]], "-", color=f"C{i}", lw=2,
                    label=f"$\\lambda_{i+1}={w[i]:.2g}$: unbounded direction")
    ax.set_aspect("equal"); ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6)
    ax.set_title(name, fontsize=10); ax.legend(fontsize=7, loc="lower left"); ax.grid(alpha=.2)
fig.suptitle("Level sets of $x^T A x$ (black: $x^T A x = 1$) with eigenvector axes")
fig.tight_layout(); fig.savefig(FIG / "quadratic_forms.png", dpi=120); plt.close(fig)

# Example 4.3 numbers.
A43 = cases["positive definite (Example 4.3)"]
w43, Q43 = eigh_desc(A43)
assert np.allclose(w43, [6, 1])
assert np.allclose(np.abs(Q43[:, 0]), np.abs([2, 1]) / np.sqrt(5))
assert np.allclose(np.abs(Q43[:, 1]), np.abs([1, -2]) / np.sqrt(5))
ang = np.degrees(np.arctan2(Q43[1, 0], Q43[0, 0])) % 180
assert np.isclose(ang, np.degrees(np.arctan(0.5)))      # minor axis at 26.57 degrees
# The symmetric part is all that matters (Theorem 3.14a): A and [[5,4],[0,2]] give the same form.
A_nonsym = np.array([[5, 4], [0, 2.]])
pts = rng.standard_normal((100, 2))
assert np.allclose(np.einsum("ni,ij,nj->n", pts, A43, pts), np.einsum("ni,ij,nj->n", pts, A_nonsym, pts))
print(f"Example 4.3: eigenvalues {w43}, semi-axes {1/np.sqrt(w43)}, minor-axis angle {ang:.2f} degrees")

# %% [markdown]
# ## Part 6: The Rayleigh quotient on the unit circle
#
# Theorem 3.17: λ_min ≤ R(x) ≤ λ_max with equality exactly at eigenvectors. For a 2×2 symmetric
# A, parametrize the unit circle by θ and plot R(θ) = x(θ)ᵀ A x(θ).


# %%
def rayleigh(A, x):
    """R(x) = x^T A x / x^T x, eq. (3.13) (x^† for complex)."""
    return np.real(inner(x, A @ x) / inner(x, x))


theta = np.linspace(0, np.pi, 1441)
Xc = np.c_[np.cos(theta), np.sin(theta)]
R = np.einsum("ni,ij,nj->n", Xc, A43, Xc)
assert np.isclose(R.max(), 6, atol=1e-5) and np.isclose(R.min(), 1, atol=1e-5)
assert np.isclose(np.degrees(theta[R.argmax()]), np.degrees(np.arctan(0.5)), atol=0.2)   # argmax at q_1
assert np.isclose(np.degrees(theta[R.argmin()]), 180 - np.degrees(np.arctan(2)), atol=0.2)  # argmin at q_2
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(np.degrees(theta), R)
ax.axhline(6, ls="--", color="C1", label="$\\lambda_{max} = 6$"); ax.axhline(1, ls="--", color="C2", label="$\\lambda_{min} = 1$")
ax.axvline(np.degrees(np.arctan(0.5)), color="C1", alpha=.4); ax.axvline(180 - np.degrees(np.arctan(2)), color="C2", alpha=.4)
ax.set_xlabel("angle of unit vector x (degrees)"); ax.set_ylabel("R(x) = $x^T A x$")
ax.set_title("Rayleigh quotient of [[5,2],[2,2]] on the unit circle"); ax.legend(); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "rayleigh_circle.png", dpi=120); plt.close(fig)

# Exercise A6 and random checks of the bounds.
A6 = np.array([[3, 1], [1, 3.]])
assert np.allclose([rayleigh(A6, np.array(v, float)) for v in [(1, 0), (1, 1), (1, -1), (2, 1)]], [3, 4, 2, 3.8])
for _ in range(20):
    M = rng.standard_normal((5, 5)); M = M + M.T
    w, Q = eigh_desc(M)
    xs = rng.standard_normal((500, 5))
    Rs = np.einsum("ni,ij,nj->n", xs, M, xs) / np.einsum("ni,ni->n", xs, xs)
    assert np.all(Rs <= w[0] + 1e-10) and np.all(Rs >= w[-1] - 1e-10)
    assert np.isclose(rayleigh(M, Q[:, 0]), w[0]) and np.isclose(rayleigh(M, Q[:, -1]), w[-1])
    # Courant–Fischer (3.14), k = 2: on S = span(q_1, q_2), R >= λ_2 everywhere with min exactly λ_2.
    alpha = np.linspace(0, np.pi, 181)                       # includes alpha = pi/2, i.e. x = q_2
    S = Q[:, :2] @ np.vstack([np.cos(alpha), np.sin(alpha)])
    R_S = np.einsum("in,ij,jn->n", S, M, S) / np.einsum("in,in->n", S, S)
    assert np.all(R_S >= w[1] - 1e-10) and np.isclose(R_S.min(), w[1])
print("Part 6 checks passed: Rayleigh quotient extremes are the eigenvalues at the eigenvectors.")

# %% [markdown]
# ## Part 7: Three tests for positive definiteness
#
# Theorem 3.15: (i) eigenvalues > 0, (iii) Cholesky A = RᵀR exists, (iv) pivots > 0.
# Example 4.6: tridiag(−1, 2, −1) has eigenvalues 2 − √2, 2, 2 + √2 and pivots 2, 3/2, 4/3.


# %%
def pivots(A):
    """Pivots of Gaussian elimination WITHOUT row exchanges (Theorem 3.15 (iv)).
    Stops and returns what it has if a zero pivot appears."""
    U = np.array(A, float); n = U.shape[0]; piv = []
    for k in range(n):
        if abs(U[k, k]) < 1e-14:
            piv.append(0.0); break
        piv.append(U[k, k])
        for i in range(k + 1, n):
            U[i, k:] -= U[i, k] / U[k, k] * U[k, k:]
    return np.array(piv)


def is_pd_eig(A):
    return bool(np.all(np.linalg.eigvalsh(A) > 0))


def is_pd_chol(A):
    try:
        np.linalg.cholesky(A); return True
    except np.linalg.LinAlgError:
        return False


def is_pd_pivots(A):
    p = pivots(A)
    return bool(len(p) == A.shape[0] and np.all(p > 0))


A46 = np.array([[2, -1, 0], [-1, 2, -1], [0, -1, 2.]])
assert np.allclose(eigh_desc(A46)[0], [2 + np.sqrt(2), 2, 2 - np.sqrt(2)])
assert np.allclose(pivots(A46), [2, 1.5, 4 / 3])
assert np.isclose(np.prod(pivots(A46)), np.linalg.det(A46)) and np.isclose(np.linalg.det(A46), 4)
R46 = np.linalg.cholesky(A46).T            # numpy returns lower L with A = L L^T; R = L^T
R_hand = np.array([[np.sqrt(2), -1 / np.sqrt(2), 0], [0, np.sqrt(1.5), -np.sqrt(2 / 3)], [0, 0, np.sqrt(4 / 3)]])
assert np.allclose(R46, R_hand) and np.allclose(R46.T @ R46, A46)
assert np.allclose(np.diag(R46) ** 2, pivots(A46))   # diagonal of the Cholesky factor = sqrt(pivots)
# (ii) as a sum of squares: x1² + (x1−x2)² + (x2−x3)² + x3².
xs = rng.standard_normal((200, 3))
assert np.allclose(np.einsum("ni,ij,nj->n", xs, A46, xs),
                   xs[:, 0]**2 + (xs[:, 0] - xs[:, 1])**2 + (xs[:, 1] - xs[:, 2])**2 + xs[:, 2]**2)

n_pd = n_not = 0
for _ in range(60):
    n = int(rng.integers(2, 7))
    if rng.random() < 0.5:
        B = rng.standard_normal((n, n)); M = B.T @ B + 0.05 * np.eye(n)   # SPD by Theorem 3.15 (iii)
    else:
        M = rng.standard_normal((n, n)); M = M + M.T                       # usually indefinite
    verdicts = (is_pd_eig(M), is_pd_chol(M), is_pd_pivots(M))
    assert len(set(verdicts)) == 1, (M, verdicts)
    n_pd += verdicts[0]; n_not += not verdicts[0]
assert n_pd > 10 and n_not > 10
# Theorem 3.16: Gram matrices are PSD; PD iff independent columns.
G_dep = rng.standard_normal((3, 5)); G_dep = G_dep.T @ G_dep      # 5x5 of rank 3: PSD, not PD
assert np.all(np.linalg.eigvalsh(G_dep) > -1e-10) and not is_pd_eig(G_dep)
# Definiteness is about eigenvalues, not entries (§8).
assert not is_pd_eig(np.array([[1, 2], [2, 1.]])) and np.all(np.linalg.eigvalsh(np.array([[1, -1], [-1, 1.]])) > -1e-12)
print(f"Part 7 checks passed: eigenvalue, Cholesky, and pivot tests agree on {n_pd} SPD and {n_not} non-SPD matrices.")

# %% [markdown]
# ## Part 8: Pauli matrices, unitaries, normal matrices, and the Bloch sphere
#
# Proposition 3.22: σ_x, σ_y, σ_z are Hermitian AND unitary, so their eigenvalues are real and on
# the unit circle: ±1. For a unit vector n, n·σ has eigenvectors (cos θ/2, e^{iφ} sin θ/2) and
# (sin θ/2, −e^{iφ} cos θ/2). Theorem 3.20: unitary eigenvalues have |λ| = 1. Theorem 3.21:
# U Λ U^† is normal.


# %%
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2)
for s in (sx, sy, sz):
    assert np.allclose(s, s.conj().T)                     # Hermitian
    assert np.allclose(s.conj().T @ s, I2)                # unitary
    assert np.allclose(s @ s, I2)                         # squares to I
    assert np.allclose(sorted(np.linalg.eigvals(s).real), [-1, 1])
assert not np.allclose(sy, sy.T)                          # Hermitian but NOT symmetric (§8)

# σ_y eigenvectors and spectral decomposition (Example 4.4).
yp = np.array([1, 1j]) / np.sqrt(2); ym = np.array([1, -1j]) / np.sqrt(2)
assert np.allclose(sy @ yp, yp) and np.allclose(sy @ ym, -ym)
assert np.isclose(inner(yp, ym), 0)                       # orthogonal (conjugate the first slot!)
Pp = np.outer(yp, yp.conj()); Pm = np.outer(ym, ym.conj())   # |q><q| = q q^†
assert np.allclose(Pp, 0.5 * np.array([[1, -1j], [1j, 1]])) and np.allclose(Pp + Pm, I2)
assert np.allclose(Pp - Pm, sy)                           # σ_y = (+1) P_+ + (−1) P_−

# Exercise A5: symmetric complex matrices misbehave.
S_a5 = np.array([[0, 1j], [1j, 0]]); N_a5 = np.array([[1, 1j], [1j, -1]])
assert np.allclose(S_a5, S_a5.T) and not np.allclose(S_a5, S_a5.conj().T)
assert np.allclose(sorted(np.linalg.eigvals(S_a5).imag), [-1, 1])        # eigenvalues ±i
assert np.allclose(N_a5 @ N_a5, 0) and np.linalg.matrix_rank(N_a5) == 1  # nilpotent, one eigenvector


def n_dot_sigma(theta, phi):
    """n·σ for n = (sinθ cosφ, sinθ sinφ, cosθ), eq. (3.15)."""
    n = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
    return n, n[0] * sx + n[1] * sy + n[2] * sz


def bloch_eigvecs(theta, phi):
    """|n,+> = (cos θ/2, e^{iφ} sin θ/2),  |n,-> = (sin θ/2, −e^{iφ} cos θ/2)  (Proposition 3.22)."""
    plus = np.array([np.cos(theta / 2), np.exp(1j * phi) * np.sin(theta / 2)])
    minus = np.array([np.sin(theta / 2), -np.exp(1j * phi) * np.cos(theta / 2)])
    return plus, minus


for _ in range(50):
    th, ph = rng.uniform(0, np.pi), rng.uniform(0, 2 * np.pi)
    n, H = n_dot_sigma(th, ph)
    assert np.allclose(H, H.conj().T) and np.allclose(H @ H, I2)
    assert np.allclose(sorted(np.linalg.eigvalsh(H)), [-1, 1])
    vp, vm = bloch_eigvecs(th, ph)
    assert np.allclose(H @ vp, vp) and np.allclose(H @ vm, -vm)
    assert np.isclose(inner(vp, vm), 0) and np.isclose(inner(vp, vp), 1)
    # The Bloch vector of |n,+> is n itself: <ψ|σ_k|ψ> = n_k.
    assert np.allclose([np.real(inner(vp, s @ vp)) for s in (sx, sy, sz)], n)

# Example 4.4: θ = π/3, φ = π/2.
n44, H44 = n_dot_sigma(np.pi / 3, np.pi / 2)
assert np.allclose(n44, [0, np.sqrt(3) / 2, 0.5])
assert np.allclose(H44, [[0.5, -1j * np.sqrt(3) / 2], [1j * np.sqrt(3) / 2, -0.5]])
vp44, vm44 = bloch_eigvecs(np.pi / 3, np.pi / 2)
assert np.allclose(vp44, [np.sqrt(3) / 2, 0.5j]) and np.allclose(vm44, [0.5, -1j * np.sqrt(3) / 2])
ket0 = np.array([1, 0], complex)
p_plus, p_minus = abs(inner(vp44, ket0)) ** 2, abs(inner(vm44, ket0)) ** 2      # (3.16)
assert np.isclose(p_plus, 0.75) and np.isclose(p_minus, 0.25) and np.isclose(p_plus + p_minus, 1)
assert np.isclose(np.real(inner(ket0, H44 @ ket0)), (+1) * p_plus + (-1) * p_minus)   # <A> = Σ λ p
print(f"Example 4.4: Born probabilities for |0> measured along n(θ=π/3, φ=π/2): p+ = {p_plus:.3f}, p- = {p_minus:.3f}")

# Theorem 3.20: random unitaries have eigenvalues on the unit circle and preserve inner products.
for _ in range(20):
    Z = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))
    Uq, _ = np.linalg.qr(Z)
    assert np.allclose(np.abs(np.linalg.eigvals(Uq)), 1)
    u, v = rng.standard_normal(4) + 0j, rng.standard_normal(4) + 0j
    assert np.isclose(inner(Uq @ u, Uq @ v), inner(u, v))
    # Theorem 3.21: U Λ U^† with any complex diagonal Λ is normal.
    Lam = np.diag(rng.standard_normal(4) + 1j * rng.standard_normal(4))
    Nn = Uq @ Lam @ Uq.conj().T
    assert np.allclose(Nn @ Nn.conj().T, Nn.conj().T @ Nn)
# A diagonalizable but non-normal matrix: eigenvectors not orthogonal, cannot be U Λ U^†.
Nd = np.array([[1, 1], [0, 2.]])
assert not np.allclose(Nd @ Nd.T, Nd.T @ Nd)
# Exercise B6(c): exp(−iθ n·σ/2) = cos(θ/2) I − i sin(θ/2) n·σ.
th = 1.234; n, H = n_dot_sigma(0.7, 2.1)
assert np.allclose(sla.expm(-1j * th / 2 * H), np.cos(th / 2) * I2 - 1j * np.sin(th / 2) * H)
assert np.allclose(sla.expm(-1j * 2 * np.pi / 2 * H), -I2)   # θ = 2π gives −I

# Bloch sphere figure: a few n directions, each labelled by its |n,+> eigenvector.
fig = plt.figure(figsize=(6, 6)); ax = fig.add_subplot(111, projection="3d")
uu, vv = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
ax.plot_wireframe(np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv), color="lightgray", lw=.4)
points = {"+z: |0>": (0, 0), "-z: |1>": (np.pi, 0), "+x: (|0>+|1>)/√2": (np.pi / 2, 0),
          "+y: (|0>+i|1>)/√2": (np.pi / 2, np.pi / 2), "Example 4.4 n": (np.pi / 3, np.pi / 2)}
for k, (label, (th, ph)) in enumerate(points.items()):
    n, _ = n_dot_sigma(th, ph)
    ax.quiver(0, 0, 0, *n, color=f"C{k}", arrow_length_ratio=.12)
    ax.text(*(1.15 * n), label, fontsize=7, color=f"C{k}")
ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1); ax.set_box_aspect((1, 1, 1))
ax.set_title("Bloch sphere: n ↔ eigenvector $|n,+\\rangle$ of $n\\cdot\\sigma$")
fig.tight_layout(); fig.savefig(FIG / "bloch_sphere.png", dpi=120); plt.close(fig)
print("Part 8 checks passed: Pauli matrices, n·σ eigenvectors, unitaries, normality.")

# %% [markdown]
# ## Part 9: The 2D stress tensor and Mohr's circle (Example 4.5)
#
# σ = [[80, 30], [30, 20]] MPa. Principal stresses = eigenvalues 50 ± 30√2 ≈ 92.43, 7.57 MPa at
# θ_p = 22.5°. The normal stress on a plane with unit normal n is σ_n = nᵀσn, a Rayleigh quotient,
# so Theorem 3.17 bounds it by the principal stresses; the shear is τ_n = tᵀσn. As θ varies,
# (σ_n, τ_n) traces the circle centered at (50, 0) with radius 30√2: Mohr's circle.


# %%
sigma = np.array([[80, 30], [30, 20.]])
w_s, Q_s = eigh_desc(sigma)
m, d, r = 50.0, 30.0, 30 * np.sqrt(2)
assert np.allclose(w_s, [m + r, m - r])
assert np.isclose(w_s[0], 92.4264, atol=1e-4) and np.isclose(w_s[1], 7.5736, atol=1e-4)
theta_p = np.degrees(np.arctan2(Q_s[1, 0], Q_s[0, 0])) % 180
assert np.isclose(theta_p, 22.5) and np.isclose(np.degrees(0.5 * np.arctan2(2 * 30, 80 - 20)), 22.5)
assert np.isclose(np.tan(np.radians(22.5)), np.sqrt(2) - 1)
assert np.allclose(sigma @ Q_s[:, 0], w_s[0] * Q_s[:, 0])

th = np.linspace(0, np.pi, 721)
Nn = np.c_[np.cos(th), np.sin(th)]; Tt = np.c_[-np.sin(th), np.cos(th)]
sig_n = np.einsum("ni,ij,nj->n", Nn, sigma, Nn)
tau_n = np.einsum("ni,ij,nj->n", Tt, sigma, Nn)
assert np.allclose(sig_n, m + d * np.cos(2 * th) + 30 * np.sin(2 * th))
assert np.allclose(tau_n, -d * np.sin(2 * th) + 30 * np.cos(2 * th))
assert np.allclose((sig_n - m) ** 2 + tau_n**2, r**2)          # the circle
assert np.isclose(sig_n.max(), w_s[0]) and np.isclose(sig_n.min(), w_s[1])   # Rayleigh bounds
assert np.isclose(np.abs(tau_n).max(), r)                      # max shear = radius
assert np.isclose(np.degrees(th[sig_n.argmax()]), 22.5, atol=0.3)
assert np.isclose(np.degrees(th[np.abs(tau_n).argmax()]) % 90, 67.5, atol=0.3)   # 45° from principal
assert abs(tau_n[sig_n.argmax()]) < 0.1                        # no shear on a principal plane

fig, ax = plt.subplots(figsize=(6.5, 5))
ax.plot(sig_n, tau_n, label="$(\\sigma_n, \\tau_n)$ over all planes")
ax.plot([w_s[0], w_s[1]], [0, 0], "ko"); ax.plot(m, 0, "k+", ms=10)
ax.annotate(f"$\\sigma_1 = {w_s[0]:.2f}$", (w_s[0], 0), textcoords="offset points", xytext=(5, 8))
ax.annotate(f"$\\sigma_2 = {w_s[1]:.2f}$", (w_s[1], 0), textcoords="offset points", xytext=(-60, 8))
ax.plot([80], [30], "rs", label="x-face $(\\sigma_{xx}, \\tau_{xy}) = (80, 30)$")
ax.plot([20], [-30], "r^", label="y-face $(\\sigma_{yy}, -\\tau_{xy}) = (20, -30)$")
ax.plot([80, 20], [30, -30], "r--", lw=.8)
ax.plot([m, m], [r, -r], "g:", label=f"$\\tau_{{max}} = {r:.2f}$")
ax.set_aspect("equal"); ax.axhline(0, color="k", lw=.5); ax.grid(alpha=.3); ax.legend(fontsize=8, loc="upper left")
ax.set_xlabel("normal stress $\\sigma_n$ (MPa)"); ax.set_ylabel("shear stress $\\tau_n$ (MPa)")
ax.set_title("Mohr's circle = the spectral theorem for the stress tensor")
fig.tight_layout(); fig.savefig(FIG / "mohr_circle.png", dpi=120); plt.close(fig)
print(f"Example 4.5: principal stresses {w_s.round(2)} MPa at θ_p = {theta_p:.1f}°, max shear {r:.2f} MPa")

# %% [markdown]
# ## Your turn
#
# 1. **Symmetric ⟂ antisymmetric.** Using the Frobenius inner product (3.3), show numerically
#    that the symmetric part S = (A + Aᵀ)/2 and antisymmetric part K = (A − Aᵀ)/2 of a random A
#    are orthogonal, that ||A||_F² = ||S||_F² + ||K||_F² (Pythagoras), and that S is the symmetric
#    matrix nearest to A: ||A − S||_F ≤ ||A − M||_F for random symmetric M. (Exercise B2.)
# 2. **Hessian classification.** With SymPy, compute the Hessian of f(x, y) = x³ − 3xy² at a
#    general point, find its eigenvalues, locate the critical points, and classify them using
#    Theorem 3.19; then do h(x, y) = x³ − 3x + y². (Exercise B3.)
# 3. **AᵀA and AAᵀ share nonzero eigenvalues.** For a random 5×3 A, compare the eigenvalues of
#    AᵀA (3×3) and AAᵀ (5×5). (Exercise B5; the SVD in LA-9 starts here.)
# 4. (Extension) Implement *cyclic* Jacobi (sweep all (p, q) pairs in order instead of searching
#    for the largest) and compare the number of rotations with the classical version.
# 5. (Extension) Check Cauchy interlacing (Exercise C1) on random symmetric 6×6 matrices.
#
# Reference implementations of 1–3 follow so you can compare.


# %%
# 1. Symmetric and antisymmetric parts are Frobenius-orthogonal (B2).
A = rng.standard_normal((4, 4))
S, K = 0.5 * (A + A.T), 0.5 * (A - A.T)
assert np.isclose(frobenius(S, K), 0)
assert np.isclose(np.linalg.norm(A, "fro") ** 2, np.linalg.norm(S, "fro") ** 2 + np.linalg.norm(K, "fro") ** 2)
for _ in range(50):
    M = rng.standard_normal((4, 4)); M = M + M.T
    assert np.linalg.norm(A - S, "fro") <= np.linalg.norm(A - M, "fro") + 1e-12
# dim Sym + dim Antisym = n(n+1)/2 + n(n-1)/2 = n^2.
assert 4 * 5 // 2 + 4 * 3 // 2 == 16

# 2. Hessians with SymPy (B3, Theorem 3.19).
xs_, ys_ = sp.symbols("x y", real=True)
f = xs_**3 - 3 * xs_ * ys_**2
grad = [sp.diff(f, v) for v in (xs_, ys_)]
Hf = sp.hessian(f, (xs_, ys_))
assert Hf == sp.Matrix([[6 * xs_, -6 * ys_], [-6 * ys_, -6 * xs_]])
crit = sp.solve(grad, [xs_, ys_], dict=True)
assert crit == [{xs_: 0, ys_: 0}]                                  # the only critical point
assert Hf.subs({xs_: 0, ys_: 0}) == sp.zeros(2)                    # Hessian vanishes: test inconclusive
eigs_H = list(Hf.eigenvals().keys())
assert {sp.simplify(e**2) for e in eigs_H} == {36 * xs_**2 + 36 * ys_**2}   # ±6 sqrt(x²+y²): indefinite off the origin
assert Hf.subs({xs_: 1, ys_: 0}) == sp.diag(6, -6)
t = sp.symbols("t", real=True)
assert f.subs({xs_: t, ys_: 0}) == t**3                            # changes sign through 0: a saddle (monkey saddle)
h = xs_**3 - 3 * xs_ + ys_**2
crit_h = sp.solve([sp.diff(h, v) for v in (xs_, ys_)], [xs_, ys_], dict=True)
assert sorted(c[xs_] for c in crit_h) == [-1, 1] and all(c[ys_] == 0 for c in crit_h)
Hh = sp.hessian(h, (xs_, ys_))
assert Hh.subs({xs_: 1, ys_: 0}) == sp.diag(6, 2)                  # PD: local minimum
assert Hh.subs({xs_: -1, ys_: 0}) == sp.diag(-6, 2)                # indefinite: saddle
print("Your turn 2: f = x^3 - 3xy^2 has H = 0 at its only critical point (monkey saddle);",
      "h = x^3 - 3x + y^2 has a min at (1,0) and a saddle at (-1,0).")

# 3. A^T A and A A^T share nonzero eigenvalues (B5).
A53 = rng.standard_normal((5, 3))
w_small = np.sort(np.linalg.eigvalsh(A53.T @ A53))[::-1]
w_big = np.sort(np.linalg.eigvalsh(A53 @ A53.T))[::-1]
assert np.allclose(w_small, w_big[:3])                             # same three nonzero eigenvalues
assert np.allclose(w_big[3:], 0, atol=1e-10)                       # plus m − n = 2 extra zeros
assert np.all(w_small > 0)                                         # Theorem 3.16: PD for full column rank
assert np.allclose(np.sqrt(w_small), np.linalg.svd(A53, compute_uv=False))   # LA-9: singular values
print("Your turn 3: nonzero eigenvalues of A^T A and A A^T agree:", w_small.round(4))

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
