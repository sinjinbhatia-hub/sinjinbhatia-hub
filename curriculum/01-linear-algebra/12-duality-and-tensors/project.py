# %% [markdown]
# # LA-12 Project: Duality, Tensors, Kronecker Products, Index Notation
#
# This script makes every theorem of LA-12 numerical. Parts: (1) dual bases and
# why they are not the transpose, (2) change of basis for operators vs bilinear
# forms and Sylvester's law, (3) Kronecker identities and the vec trick,
# (4) Schmidt decomposition and entanglement entropy of two-qubit and bipartite
# states, (5) einsum: contraction, the inertia tensor, Levi-Civita, the metric,
# and the transformation law, (6) the stress tensor and Mohr's circle,
# (7) SO(3), Rodrigues, SU(2) and the 2-to-1 map.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.linalg import expm, solve_sylvester
from scipy.spatial.transform import Rotation

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)

# %% [markdown]
# ## Part 1: The dual basis is the rows of E^{-1}, not of E^T
#
# Eq. (3.3): if the basis vectors are the columns of E, the dual basis
# functionals (as row vectors) are the rows of E^{-1}. They coincide with the
# rows of E^T only for an orthonormal basis. Example 4.1.
#
# **Build it yourself first:** write `dual_basis(E)` before reading it.


# %%
def dual_basis(E):
    """Rows of the returned matrix are the dual basis functionals e^i to the
    basis given by the columns of E, eq. (3.3): F = E^{-1}, so F E = I."""
    E = np.asarray(E, float)
    return np.linalg.inv(E)


def check_dual(E, F):
    """Assert e^i(e_j) = delta^i_j, eq. (3.2)."""
    assert np.allclose(F @ E, np.eye(E.shape[0])), "not a dual basis"


# Example 4.1, orthogonal but not orthonormal basis {(1,1), (1,-1)}.
E1 = np.array([[1, 1], [1, -1.]])          # columns are the basis vectors
F1 = dual_basis(E1)
check_dual(E1, F1)
assert np.allclose(F1, [[0.5, 0.5], [0.5, -0.5]])
assert not np.allclose(F1, E1.T)           # right directions, wrong scale
assert np.allclose(F1, E1.T / 2)
v = np.array([3, 1.])
assert np.allclose(F1 @ v, [2, 1]) and np.allclose(E1 @ (F1 @ v), v)   # v = 2 e_1 + 1 e_2
print("Example 4.1: dual basis of {(1,1),(1,-1)} =", F1.tolist(), " (E^T would give", E1.T.tolist(), ")")

# Example 4.1, non-orthogonal basis {(1,0), (1,1)}: wrong directions entirely.
E2 = np.array([[1, 1], [0, 1.]])
F2 = dual_basis(E2)
check_dual(E2, F2)
assert np.allclose(F2, [[1, -1], [0, 1]])
assert np.isclose(F2[0] @ E2[:, 1], 0)     # e^1 is orthogonal to e_2, not parallel to e_1
print("Example 4.1: dual basis of {(1,0),(1,1)} =", F2.tolist())

# Random bases: dual basis == transpose iff orthonormal.
for _ in range(20):
    E = rng.standard_normal((4, 4))
    check_dual(E, dual_basis(E))
    Q, _ = np.linalg.qr(E)                 # an orthonormal basis
    assert np.allclose(dual_basis(Q), Q.T)
print("Part 1 passed: dual basis = rows of E^{-1}; equals E^T only for orthonormal bases.")

# %% [markdown]
# ## Part 2: Change of basis, two ways (Theorem 3.11), and Sylvester (Theorem 3.12)
#
# Under e'_j = sum_i P_ij e_i: operator A -> P^{-1} A P, bilinear form B -> P^T B P.
# Example 4.2: same numbers M, same P, two answers. Then: similarity preserves
# eigenvalues; congruence preserves only the signature (p, q, z).


# %%
def change_basis_operator(A, P):
    """Matrix of the operator A in the basis whose vectors are the columns of P: P^{-1} A P, eq. (3.5)."""
    return np.linalg.solve(P, A @ P)


def change_basis_form(B, P):
    """Matrix of the bilinear form B in the basis whose vectors are the columns of P: P^T B P, eq. (3.5)."""
    return P.T @ B @ P


def signature(B, tol=1e-9):
    """(p, q, z) = numbers of positive, negative, zero eigenvalues of a symmetric B (Theorem 3.12)."""
    w = np.linalg.eigvalsh((B + B.T) / 2)
    return int(np.sum(w > tol)), int(np.sum(w < -tol)), int(np.sum(np.abs(w) <= tol))


M = np.array([[2, 1], [1, 3.]])
P = np.array([[1, 1], [0, 1.]])
A_new = change_basis_operator(M, P)
B_new = change_basis_form(M, P)
assert np.allclose(A_new, [[1, -1], [1, 4]]), A_new
assert np.allclose(B_new, [[2, 3], [3, 7]]), B_new
assert np.isclose(np.trace(A_new), 5) and np.isclose(np.trace(B_new), 9)
assert np.allclose(np.sort(np.linalg.eigvals(A_new).real), np.sort(np.linalg.eigvalsh(M)))      # (5 ± sqrt5)/2
assert np.allclose(np.linalg.eigvalsh(B_new), [(9 - np.sqrt(61)) / 2, (9 + np.sqrt(61)) / 2])
assert np.isclose(np.linalg.det(B_new), np.linalg.det(P) ** 2 * np.linalg.det(M))
assert signature(M) == signature(B_new) == (2, 0, 0)
print("Example 4.2: P^-1 M P =", A_new.tolist(), "  P^T M P =", B_new.tolist())

# Random checks. Build symmetric forms of every signature in R^4.
for p_, q_ in [(4, 0), (3, 1), (2, 2), (2, 1), (1, 1), (0, 3)]:
    z_ = 4 - p_ - q_
    D = np.diag([1.0] * p_ + [-1.0] * q_ + [0.0] * z_)
    for _ in range(10):
        P = rng.standard_normal((4, 4))
        while abs(np.linalg.det(P)) < 0.1:
            P = rng.standard_normal((4, 4))
        B = change_basis_form(D, P)
        assert signature(B) == (p_, q_, z_), "Sylvester violated"
        # eigenvalues of the form are NOT invariant (unless P is orthogonal)
        assert not np.allclose(np.sort(np.linalg.eigvalsh(B)), np.sort(np.diag(D)))
        # ... but for an operator, eigenvalues ARE invariant under P^{-1} D P
        Aop = change_basis_operator(D, P)
        assert np.allclose(np.sort(np.linalg.eigvals(Aop).real), np.sort(np.diag(D)))
        assert np.isclose(np.trace(Aop), np.trace(D))
        # Orthogonal P: the two laws coincide.
        Q, _ = np.linalg.qr(P)
        assert np.allclose(change_basis_operator(D, Q), change_basis_form(D, Q))
# Exercise B4: q = x^2 + 4xy + 2xz + 2y^2 + z^2 has signature (2, 1, 0); the P from completing squares diagonalizes it.
B4 = np.array([[1, 2, 1], [2, 2, 0], [1, 0, 1.]])
r2_ = np.sqrt(2)
P4 = np.array([[1, -r2_, 1 / r2_], [0, 1 / r2_, -1 / r2_], [0, 0, 1 / r2_]])
assert np.allclose(change_basis_form(B4, P4), np.diag([1, -1, 1])) and signature(B4) == (2, 1, 0)
assert np.isclose(np.linalg.det(B4), -4)
print("Part 2 passed: congruence keeps only the signature; similarity keeps the spectrum; they agree for orthogonal P.")

# %% [markdown]
# ## Part 3: Kronecker products (Theorem 3.17) and the vec trick (Theorem 3.18)
#
# (A kron B)_{(i,k),(j,l)} = A_ij B_kl with row-major pair ordering = np.kron.
# Identities (a)-(f), Example 4.4, and vec(AXB) = (B^T kron A) vec(X) in
# column-major form, plus the row-major variant (A kron B^T).


# %%
def same_multiset(x, y, tol=1e-7):
    """True if the complex arrays x and y agree as multisets (greedy nearest matching, tolerance tol)."""
    x = list(np.asarray(x, complex)); y = list(np.asarray(y, complex))
    if len(x) != len(y):
        return False
    for a in x:
        k = int(np.argmin(np.abs(np.array(y) - a)))
        if abs(y[k] - a) > tol:
            return False
        y.pop(k)
    return True


def kron_by_hand(A, B):
    """Kronecker product from eq. (3.11): entry ((i,k),(j,l)) = A_ij B_kl, row-major pair ordering."""
    A = np.asarray(A); B = np.asarray(B)
    m, n = A.shape; p, q = B.shape
    K = np.zeros((m * p, n * q), dtype=np.result_type(A, B))
    for i in range(m):
        for j in range(n):
            K[i * p:(i + 1) * p, j * q:(j + 1) * q] = A[i, j] * B      # block (i,j) is A_ij B, eq. (3.12)
    return K


def ket(bits):
    """Computational basis state |b_0 b_1 ...> with qubit 0 the LEFTMOST Kronecker factor."""
    v = np.array([1.0])
    for b in bits:
        v = np.kron(v, np.eye(2)[int(b)])
    return v


def vec(X):
    """Column-major vectorization (stack columns), the convention of eq. (3.13)."""
    return np.asarray(X).flatten(order="F")


def unvec(x, shape):
    return np.asarray(x).reshape(shape, order="F")


def solve_sylvester_vec(A, B, C):
    """Solve AX + XB = C via (I kron A + B^T kron I) vec X = vec C, eq. (3.13). (Exercise B7 reference.)"""
    n, m = C.shape
    K = np.kron(np.eye(m), A) + np.kron(B.T, np.eye(n))
    return unvec(np.linalg.solve(K, vec(C)), (n, m))


# Example 4.4.
A = np.array([[2, 1], [0, 3.]]); B = np.array([[1, 2], [2, 1.]])
K = np.kron(A, B)
assert np.allclose(K, kron_by_hand(A, B))
assert np.allclose(K, [[2, 4, 1, 2], [4, 2, 2, 1], [0, 0, 3, 6], [0, 0, 6, 3]])
assert np.allclose(np.sort(np.linalg.eigvals(K).real), [-3, -2, 6, 9])
assert np.isclose(np.trace(K), 10) and np.isclose(np.linalg.det(K), 324)
assert np.allclose(K @ np.kron([1, 0], [1, 1]), 6 * np.kron([1, 0], [1, 1]))     # x kron y is an eigenvector
print("Example 4.4: eigenvalues of A kron B =", np.sort(np.linalg.eigvals(K).real), " tr =", np.trace(K))

fig, ax = plt.subplots(figsize=(4.2, 4))
ax.imshow(K, cmap="coolwarm", vmin=-6, vmax=6)
for (i, j), val in np.ndenumerate(K):
    ax.text(j, i, f"{val:g}", ha="center", va="center", fontsize=11)
for s in (1.5,):
    ax.axhline(s, color="k", lw=1.5); ax.axvline(s, color="k", lw=1.5)
ax.set_xticks(range(4)); ax.set_yticks(range(4))
ax.set_xticklabels(["(0,0)", "(0,1)", "(1,0)", "(1,1)"]); ax.set_yticklabels(["(0,0)", "(0,1)", "(1,0)", "(1,1)"])
ax.set_title("A kron B: block (i,j) is A_ij B, eq. (3.12)")
fig.tight_layout(); fig.savefig(FIG / "kron_blocks.png", dpi=120); plt.close(fig)

# Qubit ordering: |01> is index 1, |10> is index 2.
assert np.allclose(ket("01"), [0, 1, 0, 0]) and np.allclose(ket("10"), [0, 0, 1, 0])
sx = np.array([[0, 1], [1, 0.]]); sz = np.diag([1., -1.])
assert np.allclose(np.kron(sx, sz) @ ket("00"), ket("10"))          # sigma_x acts on qubit 0 (A4)
CNOT = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0.]])
assert np.allclose(CNOT @ ket("10"), ket("11")) and np.allclose(CNOT @ ket("01"), ket("01"))

# Theorem 3.17 on random matrices.
for _ in range(20):
    m, n, p, q, r, s = rng.integers(2, 5, size=6)
    A = rng.standard_normal((m, n)); C = rng.standard_normal((n, p))
    B = rng.standard_normal((q, r)); D = rng.standard_normal((r, s))
    assert np.allclose(np.kron(A, B) @ np.kron(C, D), np.kron(A @ C, B @ D))     # (a) mixed product
    assert np.allclose(np.kron(A, B).T, np.kron(A.T, B.T))                        # (b) transpose
    assert np.allclose(np.kron(A, B), kron_by_hand(A, B))                         # (3.11) by hand
    A = rng.standard_normal((3, 3)); B = rng.standard_normal((4, 4))
    assert np.allclose(np.linalg.inv(np.kron(A, B)), np.kron(np.linalg.inv(A), np.linalg.inv(B)))   # (c)
    assert np.isclose(np.trace(np.kron(A, B)), np.trace(A) * np.trace(B))                           # (d)
    lam = np.linalg.eigvals(A); mu = np.linalg.eigvals(B)
    assert same_multiset(np.linalg.eigvals(np.kron(A, B)), np.outer(lam, mu).ravel())             # (e)
    assert np.isclose(np.linalg.det(np.kron(A, B)), np.linalg.det(A) ** 4 * np.linalg.det(B) ** 3)  # (f)
    # vec trick, column-major (3.13) and row-major variant.
    A = rng.standard_normal((3, 4)); X = rng.standard_normal((4, 2)); B = rng.standard_normal((2, 5))
    assert np.allclose(vec(A @ X @ B), np.kron(B.T, A) @ vec(X))
    assert np.allclose((A @ X @ B).ravel(), np.kron(A, B.T) @ X.ravel())
    assert not np.allclose(vec(A @ X @ B), np.kron(B.T, A) @ X.ravel())   # mixing conventions fails

# Sylvester equation via the vec trick (B7).
A = rng.standard_normal((3, 3)); B = rng.standard_normal((2, 2)); C = rng.standard_normal((3, 2))
X = solve_sylvester_vec(A, B, C)
assert np.abs(A @ X + X @ B - C).max() < 1e-10
assert np.allclose(X, solve_sylvester(A, B, C))
Ksyl = np.kron(np.eye(2), A) + np.kron(B.T, np.eye(3))
assert same_multiset(np.linalg.eigvals(Ksyl), np.add.outer(np.linalg.eigvals(B), np.linalg.eigvals(A)).ravel())  # B7(c)
print("Part 3 passed: Kronecker identities (a)-(f), vec trick both orders, Sylvester equation solved via kron.")

# %% [markdown]
# ## Part 4: Schmidt decomposition (Theorem 3.16) and entanglement entropy
#
# psi in V kron W with dims (n, m): the coefficient matrix is psi.reshape(n, m)
# (row-major, qubit 0 = row index). Its SVD C = U S V^dagger gives
# psi = sum_k s_k a_k kron b_k with a_k = U[:, k], b_k = Vh[k, :], eq. (3.9).
#
# **Build it yourself first:** three lines.


# %%
def schmidt_decomposition(psi, dims):
    """Schmidt decomposition psi = sum_k s_k a_k (x) b_k, eq. (3.9), via SVD of the
    coefficient matrix (3.8). Returns (s, a, b): s[k] >= 0 descending, a[k] in C^n,
    b[k] in C^m orthonormal. Schmidt rank = len(s)."""
    n, m = dims
    C = np.asarray(psi).reshape(n, m)
    U, s, Vh = np.linalg.svd(C)
    r = int(np.sum(s > 1e-12))
    return s[:r], U[:, :r].T, Vh[:r, :]


def reconstruct(s, a, b):
    return sum(s_k * np.kron(a_k, b_k) for s_k, a_k, b_k in zip(s, a, b))


def entanglement_entropy(s):
    """S = -sum_k s_k^2 log s_k^2 (natural log) of the Schmidt coefficients; a preview of QC-7."""
    p = np.asarray(s) ** 2
    p = p[p > 1e-15]
    return float(-np.sum(p * np.log(p)))


def is_product(psi, dims):
    return len(schmidt_decomposition(psi, dims)[0]) == 1


plus = np.array([1, 1]) / np.sqrt(2); minus = np.array([1, -1]) / np.sqrt(2)
zero, one = np.eye(2)

# Product state |+>|0>: coefficient matrix [[1,0],[1,0]]/sqrt2, rank 1, entropy 0.
psi_prod = np.kron(plus, zero)
s, a, b = schmidt_decomposition(psi_prod, (2, 2))
assert len(s) == 1 and np.isclose(s[0], 1) and np.isclose(entanglement_entropy(s), 0)
assert np.allclose(psi_prod.reshape(2, 2), [[1, 0], [1, 0]] / np.sqrt(2))

# Four Bell states: Schmidt coefficients (1/sqrt2, 1/sqrt2), entropy log 2.
bells = {"Phi+": (ket("00") + ket("11")) / np.sqrt(2), "Phi-": (ket("00") - ket("11")) / np.sqrt(2),
         "Psi+": (ket("01") + ket("10")) / np.sqrt(2), "Psi-": (ket("01") - ket("10")) / np.sqrt(2)}
for name, psi in bells.items():
    s, a, b = schmidt_decomposition(psi, (2, 2))
    assert len(s) == 2 and np.allclose(s, 1 / np.sqrt(2))
    assert np.isclose(entanglement_entropy(s), np.log(2))
    assert np.allclose(reconstruct(s, a, b), psi)
    assert not is_product(psi, (2, 2))
assert np.allclose(bells["Phi+"].reshape(2, 2), np.eye(2) / np.sqrt(2))

# Example 4.3: (|00> + |01> + |10> - |11>)/2 is maximally entangled.
psi43 = (ket("00") + ket("01") + ket("10") - ket("11")) / 2
C43 = psi43.reshape(2, 2)
assert np.allclose(C43, 0.5 * np.array([[1, 1], [1, -1]]))
s, a, b = schmidt_decomposition(psi43, (2, 2))
assert np.allclose(s, [1 / np.sqrt(2), 1 / np.sqrt(2)]) and np.isclose(entanglement_entropy(s), np.log(2))
assert np.allclose(reconstruct(s, a, b), psi43)
# The hand decomposition (|0>|+> + |1>|->)/sqrt2 is one valid choice (degenerate singular values).
assert np.allclose((np.kron(zero, plus) + np.kron(one, minus)) / np.sqrt(2), psi43)
print("Example 4.3: Schmidt coefficients", np.round(s, 4), " entropy", round(entanglement_entropy(s), 4), "= ln 2")
print("            SVD's Schmidt vectors a_k =", np.round(a, 4).tolist(), " b_k =", np.round(b, 4).tolist())

# Exercise B3's state.
psiB3 = (2 * ket("00") + ket("01") + ket("10") + 2 * ket("11")) / np.sqrt(10)
s, a, b = schmidt_decomposition(psiB3, (2, 2))
assert np.allclose(s, [3 / np.sqrt(10), 1 / np.sqrt(10)])
assert np.isclose(entanglement_entropy(s), -(0.9 * np.log(0.9) + 0.1 * np.log(0.1)))

# Random two-qubit states (complex) and a random 3x4 bipartite state.
spectra = {}
for name, dims in [("random 2x2", (2, 2)), ("random 3x4", (3, 4))]:
    psi = rng.standard_normal(dims[0] * dims[1]) + 1j * rng.standard_normal(dims[0] * dims[1])
    psi /= np.linalg.norm(psi)
    s, a, b = schmidt_decomposition(psi, dims)
    assert len(s) == min(dims)                                   # generic: full Schmidt rank
    assert np.allclose(reconstruct(s, a, b), psi)
    assert np.isclose(np.sum(s ** 2), 1)                        # ||psi||^2 = sum s_k^2
    assert np.allclose(a.conj() @ a.T, np.eye(len(s))) and np.allclose(b.conj() @ b.T, np.eye(len(s)))
    assert 0 < entanglement_entropy(s) <= np.log(min(dims)) + 1e-12
    spectra[name] = s
    # Reduced density matrix of the first factor: C C^dagger has eigenvalues s_k^2.
    C = psi.reshape(dims)
    assert np.allclose(np.sort(np.linalg.eigvalsh(C @ C.conj().T))[::-1][:len(s)], s ** 2)
# Random product states have Schmidt rank 1.
for _ in range(5):
    u = rng.standard_normal(3) + 1j * rng.standard_normal(3); w = rng.standard_normal(4) + 1j * rng.standard_normal(4)
    assert is_product(np.kron(u, w), (3, 4))
# Simple tensors form an (n + m - 1)-dimensional surface, not a subspace (§3.11): the Jacobian of
# (v, w) -> v kron w is [I_n kron w | v kron I_m], of rank n + m - 1, with kernel spanned by (v, -w)
# because (c v) kron w = v kron (c w). Two qubits: 3 < 4; a 3x4 system: 6 < 12.
for n_, m_ in [(2, 2), (3, 4)]:
    v_ = rng.standard_normal(n_); w_ = rng.standard_normal(m_)
    J = np.hstack([np.kron(np.eye(n_), w_[:, None]), np.kron(v_[:, None], np.eye(m_))])
    assert J.shape == (n_ * m_, n_ + m_) and np.linalg.matrix_rank(J) == n_ + m_ - 1
    assert np.allclose(J @ np.concatenate([v_, -w_]), 0)
    h = 1e-6; dv = rng.standard_normal(n_); dw = rng.standard_normal(m_)      # J really is the derivative
    assert np.allclose((np.kron(v_ + h * dv, w_ + h * dw) - np.kron(v_ - h * dv, w_ - h * dw)) / (2 * h),
                       J @ np.concatenate([dv, dw]), atol=1e-6)

fig, ax = plt.subplots(figsize=(7, 3.6))
cases = [("|+>|0>", schmidt_decomposition(psi_prod, (2, 2))[0]), ("Bell", schmidt_decomposition(bells["Phi+"], (2, 2))[0]),
         ("Ex 4.3", schmidt_decomposition(psi43, (2, 2))[0]), ("B3", schmidt_decomposition(psiB3, (2, 2))[0]),
         ("random 2x2", spectra["random 2x2"]), ("random 3x4", spectra["random 3x4"])]
w = 0.8
for k, (name, s) in enumerate(cases):
    ax.bar(k + w * (np.arange(len(s)) - (len(s) - 1) / 2) / max(len(s), 1), s ** 2, width=w / max(len(s), 1) * 0.9,
           color="C0", edgecolor="k")
    ax.text(k, 1.02, f"S={entanglement_entropy(s):.2f}", ha="center", fontsize=8)
ax.set_xticks(range(len(cases))); ax.set_xticklabels([c[0] for c in cases])
ax.set_ylabel(r"Schmidt coefficients $\sigma_k^2$"); ax.set_ylim(0, 1.12)
ax.set_title("Schmidt spectra: one bar = product state; equal bars = maximally entangled")
fig.tight_layout(); fig.savefig(FIG / "schmidt_coefficients.png", dpi=120); plt.close(fig)
print("Part 4 passed: Schmidt rank 1 for products, (1/sqrt2, 1/sqrt2) and entropy ln 2 for Bell states, reconstruction exact.")

# %% [markdown]
# ## Part 5: einsum is index notation (§3.15–3.17)
#
# Contraction = repeated letter. The inertia tensor of Example 4.6 is
# I_ij = sum_a m_a (r_a^2 delta_ij - x_ai x_aj), i.e. 'a,a,ij->ij' minus 'a,ai,aj->ij'.
# The cross product is (a x b)_i = eps_ijk a_j b_k, eq. (3.23). The metric
# g_ij = <e_i, e_j> lowers indices, eq. (3.16). The transformation law (3.18):
# upper indices contract with P^{-1}, lower with P.

# %%
# Trace, matrix product, bilinear form.
A = rng.standard_normal((4, 4)); B = rng.standard_normal((4, 4)); u = rng.standard_normal(4); v = rng.standard_normal(4)
assert np.isclose(np.einsum("ii", A), np.trace(A))
assert np.allclose(np.einsum("ij,jk->ik", A, B), A @ B)
assert np.isclose(np.einsum("ij,i,j", B, u, v), u @ B @ v)                 # B(u, v), eq. (3.4)
assert np.isclose(np.einsum("ij,ij", A, B), np.trace(A.T @ B))             # Frobenius inner product = full contraction

# Levi-Civita symbol and the cross product (3.23); the eps-delta identity (B6).
eps = np.zeros((3, 3, 3))
for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    eps[i, j, k] = 1
for (i, j, k) in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
    eps[i, j, k] = -1
a3 = rng.standard_normal(3); b3 = rng.standard_normal(3); c3 = rng.standard_normal(3)
assert np.allclose(np.einsum("ijk,j,k->i", eps, a3, b3), np.cross(a3, b3))
d = np.eye(3)
assert np.allclose(np.einsum("ijk,imn->jkmn", eps, eps), np.einsum("jm,kn->jkmn", d, d) - np.einsum("jn,km->jkmn", d, d))
assert np.allclose(np.cross(a3, np.cross(b3, c3)), b3 * (a3 @ c3) - c3 * (a3 @ b3))
assert np.isclose(np.einsum("ijk,i,j,k", eps, a3, b3, c3), np.linalg.det(np.c_[a3, b3, c3]))   # volume form

# Example 4.6: inertia tensor of three point masses.
m = np.array([1., 2., 1.])
X = np.array([[1, 1, 0], [-1, 0, 1], [0, 1, 1.]])                 # rows are the positions r_a
r2 = np.einsum("ai,ai->a", X, X)
I = np.einsum("a,a,ij->ij", m, r2, np.eye(3)) - np.einsum("a,ai,aj->ij", m, X, X)   # eq. (3.22)
assert np.allclose(I, [[5, -1, 2], [-1, 6, -1], [2, -1, 5]])
moments, axes = np.linalg.eigh(I)
assert np.allclose(moments, [3, 5, 8])
for ax_, target in zip(axes.T, [[1, 0, -1], [1, 2, 1], [1, -1, 1]]):
    t = np.array(target, float) / np.linalg.norm(target)
    assert np.isclose(abs(ax_ @ t), 1), (ax_, t)                   # parallel up to sign
omega = rng.standard_normal(3)
assert np.isclose(omega @ I @ omega, np.sum(m * np.linalg.norm(np.cross(X, omega), axis=1) ** 2))   # 2T >= 0
print("Example 4.6: I =", I.tolist(), " principal moments", moments)

# Metric tensor in a skewed basis (§3.16): e_1 = (1,0), e_2 = (1,1).
E = np.array([[1, 1], [0, 1.]])                                   # columns e_1, e_2
g = E.T @ E                                                       # g_ij = <e_i, e_j>, eq. (3.15)
assert np.allclose(g, [[1, 1], [1, 2]]) and np.allclose(np.linalg.inv(g), [[2, -1], [-1, 1]])
vstd = np.array([2, 1.])
v_up = np.linalg.solve(E, vstd)                                   # contravariant components v^i
v_dn = g @ v_up                                                   # covariant components v_i = g_ij v^j, eq. (3.16)
assert np.allclose(v_up, [1, 1]) and np.allclose(v_dn, [2, 3])
assert np.allclose(v_dn, E.T @ vstd)                              # v_j = <e_j, v>
assert np.isclose(v_up @ v_dn, vstd @ vstd)                       # |v|^2 = v^i v_i = 5
assert np.allclose(np.linalg.inv(g) @ v_dn, v_up)                 # raising brings it back

# Transformation law (3.18) vs matrix formulas (3.5), for (1,1), (0,2), and (1,2) tensors.
for _ in range(10):
    n = 3
    P = rng.standard_normal((n, n)); Pinv = np.linalg.inv(P)
    T11 = rng.standard_normal((n, n)); B02 = rng.standard_normal((n, n)); phi = rng.standard_normal(n)
    assert np.allclose(np.einsum("ik,kl,lj->ij", Pinv, T11, P), Pinv @ T11 @ P)       # T'^i_j = P^-1 T P
    assert np.allclose(np.einsum("ki,kl,lj->ij", P, B02, P), P.T @ B02 @ P)           # B'_ij = P^T B P
    assert not np.allclose(Pinv @ T11 @ P, P.T @ T11 @ P)                            # the two laws differ
    vec_up = rng.standard_normal(n)
    assert np.allclose(np.einsum("ik,k->i", Pinv, vec_up), Pinv @ vec_up)             # v'^i = (P^-1 v)^i
    assert np.allclose(np.einsum("ki,k->i", P, phi), P.T @ phi)                       # phi'_i = (P^T phi)_i
    # (1,2) tensor T^a_{bc} = A^a_b phi_c (exercise B5).
    T12 = np.einsum("ab,c->abc", T11, phi)
    T12_law = np.einsum("ia,bj,ck,abc->ijk", Pinv, P, P, T12)
    assert np.allclose(T12_law, np.einsum("ij,k->ijk", Pinv @ T11 @ P, P.T @ phi))
    # Contraction commutes with the law: T'^i_i = T^i_i.
    assert np.isclose(np.einsum("ii", Pinv @ T11 @ P), np.einsum("ii", T11))
    # Orthogonal P: upper and lower indices transform identically.
    Q, _ = np.linalg.qr(P)
    assert np.allclose(np.einsum("ik,kl,lj->ij", Q.T, T11, Q), np.einsum("ki,kl,lj->ij", Q, T11, Q))
print("Part 5 passed: einsum contractions, Levi-Civita, metric lowering, transformation law for (1,1), (0,2), (1,2).")

# %% [markdown]
# ## Part 6: The stress tensor and Mohr's circle (§3.18, Example 4.5)
#
# sigma' = P^T sigma P for a rotation P by theta (eq. 3.19) must agree with the
# Mohr parametrization (3.21): sigma'_xx = c0 + r cos(2theta - 2theta_p),
# sigma'_xy = -r sin(2theta - 2theta_p). Principal stresses = eigenvalues.


# %%
def rot2(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


def mohr(sigma, theta):
    """(sigma'_x'x', sigma'_x'y') from Mohr's circle, eq. (3.21)."""
    sxx, syy, sxy = sigma[0, 0], sigma[1, 1], sigma[0, 1]
    c0 = (sxx + syy) / 2; r = np.hypot((sxx - syy) / 2, sxy)
    two_theta_p = np.arctan2(sxy, (sxx - syy) / 2)
    phi = 2 * theta - two_theta_p
    return c0 + r * np.cos(phi), -r * np.sin(phi), c0, r, two_theta_p / 2


sigma = np.array([[80, 30], [30, 20.]])                       # MPa
theta = np.radians(30)
P = rot2(theta)
sig_rot = P.T @ sigma @ P                                     # eq. (3.19)
sxx_m, sxy_m, c0, r, theta_p = mohr(sigma, theta)
assert np.allclose(sig_rot, [[90.9808, -10.9808], [-10.9808, 9.0192]], atol=1e-4), sig_rot
assert np.isclose(sig_rot[0, 0], sxx_m) and np.isclose(sig_rot[0, 1], sxy_m)
assert np.isclose(np.trace(sig_rot), 100)
assert np.isclose(c0, 50) and np.isclose(r, 30 * np.sqrt(2)) and np.isclose(np.degrees(theta_p), 22.5)
princ = np.linalg.eigvalsh(sigma)
assert np.allclose(princ, [c0 - r, c0 + r])                   # Mohr's circle IS the 2D spectral theorem
# Principal directions: at theta = theta_p the shear vanishes and sigma'_xx = c0 + r.
Pp = rot2(theta_p)
assert np.allclose(Pp.T @ sigma @ Pp, np.diag([c0 + r, c0 - r]), atol=1e-9)
# The Mohr formula agrees with P^T sigma P for every angle and every symmetric sigma.
for _ in range(50):
    S = rng.standard_normal((2, 2)); S = S + S.T; th = rng.uniform(0, np.pi)
    Sr = rot2(th).T @ S @ rot2(th)
    mx, mxy, *_ = mohr(S, th)
    assert np.isclose(Sr[0, 0], mx) and np.isclose(Sr[0, 1], mxy)
print(f"Example 4.5: sigma' at 30 deg = {np.round(sig_rot, 2).tolist()} MPa (matrix) ; Mohr gives ({sxx_m:.2f}, {sxy_m:.2f})")
print(f"             center {c0:.1f}, radius {r:.2f}, principal stresses {princ[::-1].round(2)}, theta_p = {np.degrees(theta_p):.1f} deg")

fig, ax = plt.subplots(figsize=(6.5, 5))
ang = np.linspace(0, 2 * np.pi, 400)
ax.plot(c0 + r * np.cos(ang), r * np.sin(ang), "k", lw=1.2)
ax.axhline(0, color="gray", lw=.8)
ax.plot([sigma[0, 0], sigma[1, 1]], [sigma[0, 1], -sigma[0, 1]], "C0o-", label="θ = 0: faces x (80, 30) and y (20, −30)")
ax.plot([sig_rot[0, 0], sig_rot[1, 1]], [sig_rot[0, 1], -sig_rot[0, 1]], "C1s-",
        label=f"θ = 30°: x' ({sig_rot[0,0]:.2f}, {sig_rot[0,1]:.2f}) and y'")
ax.plot(princ, [0, 0], "C3^", ms=9, label=f"principal: {princ[1]:.2f}, {princ[0]:.2f} MPa (θ_p = 22.5°)")
ax.plot(c0, 0, "k+", ms=10)
ax.annotate("2θ = 60° on the circle\n(clockwise with τ plotted up)", xy=(sig_rot[0, 0], sig_rot[0, 1]),
            xytext=(60, -45), fontsize=8, arrowprops=dict(arrowstyle="->", lw=.8))
ax.set_aspect("equal"); ax.grid(alpha=.3); ax.legend(fontsize=8, loc="upper left")
ax.set_xlabel("normal stress σ (MPa)"); ax.set_ylabel("shear stress τ (MPa)")
ax.set_title("Mohr's circle for σ = [[80, 30], [30, 20]] MPa, eq. (3.21)")
fig.tight_layout(); fig.savefig(FIG / "mohr_circle.png", dpi=120); plt.close(fig)
print("Part 6 passed: P^T sigma P agrees with Mohr's circle for all angles; principal stresses are the eigenvalues.")

# %% [markdown]
# ## Part 7: SO(3), Rodrigues, SU(2), and the 2-to-1 map (Appendix A)
#
# R(n, theta) = I + sin(theta) K + (1 - cos theta) K^2, K = [n]_x, eq. (A.4).
# Tangent space at I is antisymmetric (A.1); [L_i, L_j] = eps_ijk L_k (A.3);
# U = cos(theta/2) I - i sin(theta/2) n.sigma (A.7); U (v.sigma) U^dagger = (R v).sigma;
# theta = 2 pi gives U = -I but R = I.
#
# **Build it yourself first:** `rodrigues(n, theta)`.


# %%
def hat(n):
    """[n]_x, the antisymmetric matrix with [n]_x v = n x v; equals n . L with L from (A.2)."""
    n = np.asarray(n, float)
    return np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])


def rodrigues(n, theta):
    """Rotation by theta about the unit axis n, eq. (A.4)."""
    n = np.asarray(n, float); n = n / np.linalg.norm(n)
    K = hat(n)
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)


sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
pauli = np.array([sx, sy, sz])


def su2(n, theta):
    """U(n, theta) = exp(-i theta n.sigma / 2) = cos(theta/2) I - i sin(theta/2) n.sigma, eq. (A.7)."""
    n = np.asarray(n, float); n = n / np.linalg.norm(n)
    ns = np.einsum("j,jab->ab", n, pauli)
    return np.cos(theta / 2) * np.eye(2) - 1j * np.sin(theta / 2) * ns


# Example 4.7: 90 degrees about (1,1,1)/sqrt3.
n111 = np.ones(3) / np.sqrt(3)
R = rodrigues(n111, np.pi / 2)
third, s3 = 1 / 3, 1 / np.sqrt(3)
R_expected = np.array([[third, third - s3, third + s3], [third + s3, third, third - s3], [third - s3, third + s3, third]])
assert np.allclose(R, R_expected)
assert np.allclose(R @ n111, n111) and np.isclose(np.trace(R), 1) and np.isclose(np.linalg.det(R), 1)
assert np.allclose(R.T @ R, np.eye(3)) and np.allclose(np.linalg.matrix_power(R, 4), np.eye(3))
assert np.allclose(hat(n111) @ hat(n111), np.outer(n111, n111) - np.eye(3))          # K^2 = n n^T - I
assert np.allclose(rodrigues(n111, 2 * np.pi / 3), [[0, 0, 1], [1, 0, 0], [0, 1, 0]])  # 120 deg = cyclic permutation
assert np.allclose(R, Rotation.from_rotvec(np.pi / 2 * n111).as_matrix())
assert np.allclose(R, expm(np.pi / 2 * hat(n111)))                                   # R = exp(theta K)
print("Example 4.7: R(90 deg about (1,1,1)) =\n", np.round(R, 4))

fig = plt.figure(figsize=(5, 4.5)); ax = fig.add_subplot(projection="3d")
for k, (col, lab) in enumerate(zip("rgb", "xyz")):
    e = np.eye(3)[k]
    ax.quiver(0, 0, 0, *e, color=col, alpha=.35, arrow_length_ratio=.1)
    ax.quiver(0, 0, 0, *(R @ e), color=col, arrow_length_ratio=.1, label=f"R e_{lab}")
ax.quiver(0, 0, 0, *(1.3 * n111), color="k", arrow_length_ratio=.08, label="axis n")
ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1); ax.legend(fontsize=8)
ax.set_title("90° rotation about (1,1,1)/√3: faint = before, solid = after")
fig.tight_layout(); fig.savefig(FIG / "rotation_rodrigues.png", dpi=120); plt.close(fig)

# Random axis/angle: Rodrigues == scipy == expm; tangent vector at I is antisymmetric (A.1).
for _ in range(20):
    n = rng.standard_normal(3); n /= np.linalg.norm(n); th = rng.uniform(-np.pi, np.pi)
    Rn = rodrigues(n, th)
    assert np.allclose(Rn, Rotation.from_rotvec(th * n).as_matrix())
    assert np.allclose(Rn, expm(th * hat(n)))
    assert np.allclose(Rn.T @ Rn, np.eye(3)) and np.isclose(np.linalg.det(Rn), 1)
    h = 1e-6
    dR = (rodrigues(n, h) - rodrigues(n, -h)) / (2 * h)          # derivative of R(t) at t = 0
    assert np.allclose(dR, -dR.T, atol=1e-8) and np.allclose(dR, hat(n), atol=1e-8)
    # Composition about the same axis adds angles; about different axes it does not commute.
    assert np.allclose(rodrigues(n, th) @ rodrigues(n, 0.7), rodrigues(n, th + 0.7))

# Generators (A.2) and commutators (A.3): (L_i)_{jk} = -eps_ijk, [L_i, L_j] = eps_ijk L_k.
L = -eps
assert np.allclose(L[2], [[0, -1, 0], [1, 0, 0], [0, 0, 0]])
assert np.allclose(np.einsum("i,ijk->jk", n111, L), hat(n111))         # n.L = [n]_x
for i in range(3):
    for j in range(3):
        assert np.allclose(L[i] @ L[j] - L[j] @ L[i], np.einsum("k,kab->ab", eps[i, j], L))
# SU(2) generators -i sigma_j / 2 obey the same relations (A.6).
G = -0.5j * pauli
for i in range(3):
    for j in range(3):
        assert np.allclose(G[i] @ G[j] - G[j] @ G[i], np.einsum("k,kab->ab", eps[i, j], G))
assert np.allclose(pauli[0] @ pauli[1], 1j * pauli[2])                   # sigma_j sigma_k = delta I + i eps sigma_l

# SU(2) -> SO(3): U (v.sigma) U^dagger = (R v).sigma for random axis, angle, and v.
for _ in range(20):
    n = rng.standard_normal(3); n /= np.linalg.norm(n); th = rng.uniform(-2 * np.pi, 2 * np.pi)
    U = su2(n, th)
    assert np.allclose(U, expm(-0.5j * th * np.einsum("j,jab->ab", n, pauli)))
    assert np.allclose(U.conj().T @ U, np.eye(2)) and np.isclose(np.linalg.det(U), 1)
    Rn = rodrigues(n, th)
    vv = rng.standard_normal(3)
    lhs = U @ np.einsum("j,jab->ab", vv, pauli) @ U.conj().T
    rhs = np.einsum("j,jab->ab", Rn @ vv, pauli)
    assert np.allclose(lhs, rhs)
    assert np.allclose(su2(n, th + 2 * np.pi), -U)                      # the 2-to-1 property
    assert np.allclose(rodrigues(n, th + 2 * np.pi), Rn)
U2pi = su2(n111, 2 * np.pi)
assert np.allclose(U2pi, -np.eye(2)) and np.allclose(rodrigues(n111, 2 * np.pi), np.eye(3))
assert np.allclose(su2(n111, 4 * np.pi), np.eye(2))
# The e_z demonstration of A.6: U sigma_x U^dagger = cos(theta) sigma_x + sin(theta) sigma_y.
th = 0.9; Uz = su2([0, 0, 1], th)
assert np.allclose(Uz @ sx @ Uz.conj().T, np.cos(th) * sx + np.sin(th) * sy)
print("Part 7 passed: Rodrigues == scipy == expm; tangent space antisymmetric; [L_i,L_j] = eps_ijk L_k for both algebras;")
print("              U (v.sigma) U^dag = (R v).sigma; theta = 2 pi gives U = -I, R = I.")

# %% [markdown]
# ## Your turn
#
# 1. **Metric in a skewed basis of R^3.** Take e_1 = (1,0,0), e_2 = (1,1,0),
#    e_3 = (1,1,1). Compute g_ij = <e_i, e_j> (as `E.T @ E`), its inverse g^{ij},
#    and for v = (1, 2, 3) in standard coordinates find the contravariant
#    components v^i (solve E v^ = v) and the covariant components v_i = g_ij v^j.
#    Check v_i = <e_i, v> and v^i v_i = |v|^2 = 14. Then change to a second
#    skewed basis with a shear P and confirm that v^i picks up P^{-1} and v_i
#    picks up P^T (B1).
# 2. **Thin rectangular plate.** Sample N = 200_000 points uniformly on the
#    rectangle [-a/2, a/2] x [-b/2, b/2] in the xy-plane with a = 2, b = 1 and
#    total mass 1 (each point mass 1/N). Build I with the einsum of Part 5 and
#    compare its eigenvalues with the textbook values m b^2/12, m a^2/12,
#    m (a^2 + b^2)/12 (perpendicular-axis theorem: I_zz = I_xx + I_yy). Then
#    rotate the plate by 30 degrees with `rodrigues([0,0,1], pi/6)` applied to
#    every point, recompute I, and check it equals R I R^T (eq. 3.19). The
#    principal moments must not change.
# 3. **Sylvester equation.** Call `solve_sylvester_vec` on a 5x5 A, 4x4 B,
#    5x4 C. Then set B = -A (same size) and watch `np.linalg.solve` fail or
#    return garbage: explain with B7(c). Finally implement the row-major
#    version (A kron I + I kron B^T) acting on `C.ravel()` and check it agrees.
# 4. **Partial trace.** For a random two-qubit psi, compute the reduced density
#    matrix rho_0 = C C^dagger of qubit 0 and verify its eigenvalues are the
#    squared Schmidt coefficients and that -tr(rho_0 log rho_0) equals
#    `entanglement_entropy`. (QM-7 and QC-7 call this the von Neumann entropy of the reduced state.)

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
