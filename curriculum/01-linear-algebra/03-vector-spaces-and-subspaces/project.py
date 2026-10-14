# %% [markdown]
# # LA-3 Project: The Four Fundamental Subspaces from Scratch
#
# This script builds everything in LA-3 §3.3–§3.10 with no black boxes, then
# checks it against SciPy. Parts: (1) `rref_with_E` with partial pivoting on
# [A | I], returning R, the pivot columns, and E with EA = R (`rref` drops E); (2) bases for C(A), N(A),
# C(A^T), N(A^T) built on `rref`, with the dimensions r, n-r, r, m-r asserted on
# random matrices of known rank, the orthogonality N(A) ⟂ C(A^T),
# N(A^T) ⟂ C(A) and the factorization A = C R~ (exercise B6) checked numerically;
# (3) every number in worked examples 4.1–4.5 reproduced (4.5 symbolically with
# SymPy); (4) pictures of N(A) and
# C(A^T) in R^3; (5) an independence tester, `extend_to_basis`, and the rank
# inequality rank(AB) <= min(rank A, rank B).
#
# `scipy.linalg.null_space` and `scipy.linalg.orth` are used only as reference
# answers. They compute the same subspaces by the SVD (LA-9), which is the
# numerically right tool; RREF with a tolerance is the pedagogically right one.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.linalg as sla
import sympy as sp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)

# %% [markdown]
# ## Part 1: Reduced row echelon form, recording the row operations
#
# Gauss–Jordan on the augmented matrix [A | I_m]. When the left block reaches
# RREF R, the right block is the invertible E with EA = R (notes §3.3 and
# Theorem 3.31). `partial_pivoting=True` picks the largest available entry as
# the pivot (what real codes do; LA-13 explains why); `False` picks the first
# nonzero entry, which reproduces hand computations exactly. R is the same
# either way (the RREF is unique: the corollary after Theorem 3.28, and LA-2
# §3.9); E is not.
#
# **Build it yourself first:** write `rref` before reading this one.


# %%
def rref_with_E(A, tol=1e-10, partial_pivoting=True):
    """Return (R, pivots, E): R the RREF of A, pivots its pivot-column indices,
    E invertible with E @ A == R. Entries with |x| <= tol are treated as 0."""
    A = np.asarray(A, float)
    m, n = A.shape
    M = np.hstack([A, np.eye(m)])
    pivots, row = [], 0
    for col in range(n):
        if row >= m:
            break
        if partial_pivoting:
            p = row + int(np.argmax(np.abs(M[row:, col])))
        else:
            nz = np.flatnonzero(np.abs(M[row:, col]) > tol)
            p = row + int(nz[0]) if nz.size else row
        if abs(M[p, col]) <= tol:          # no pivot in this column: it is free
            M[row:, col] = 0.0
            continue
        M[[row, p]] = M[[p, row]]          # swap the pivot row up
        M[row] /= M[row, col]              # scale the pivot to 1
        for i in range(m):                 # eliminate above and below
            if i != row:
                M[i] -= M[i, col] * M[row]
        pivots.append(col)
        row += 1
    R, E = M[:, :n].copy(), M[:, n:].copy()
    R[np.abs(R) <= tol] = 0.0
    return R, pivots, E


def rref(A, tol=1e-10):
    """RREF of A and the list of pivot columns (LA-2). Returns (R, pivots)."""
    R, pivots, _ = rref_with_E(A, tol)
    return R, pivots


def is_rref(R, pivots, tol=1e-10):
    """Check the defining properties of RREF used in Theorems 3.8, 3.28, 3.29."""
    m, n = R.shape
    r = len(pivots)
    if np.any(np.abs(R[r:]) > tol):                       # rows below r are zero
        return False
    for i, p in enumerate(pivots):
        e = np.zeros(m); e[i] = 1.0
        if not np.allclose(R[:, p], e, atol=tol):         # pivot column p_i is e_i
            return False
        if np.any(np.abs(R[i, :p]) > tol):                # zeros left of the pivot
            return False
    return True


for _ in range(200):
    m, n = rng.integers(1, 7, size=2)
    r_true = int(rng.integers(0, min(m, n) + 1))
    A = rng.standard_normal((m, r_true)) @ rng.standard_normal((r_true, n)) if r_true else np.zeros((m, n))
    for pp in (True, False):
        R, piv, E = rref_with_E(A, partial_pivoting=pp)
        assert is_rref(R, piv), (A, R, piv)
        assert np.allclose(E @ A, R, atol=1e-8), "E A must equal R"
        assert len(piv) == r_true == np.linalg.matrix_rank(A), (len(piv), r_true)
        assert abs(np.linalg.det(E)) > 1e-12, "E must be invertible"
    # The RREF is unique (corollary after Theorem 3.28): both pivoting strategies give the same R.
    assert np.allclose(rref_with_E(A, partial_pivoting=True)[0], rref_with_E(A, partial_pivoting=False)[0], atol=1e-8)
print("Part 1: rref passes structure, EA = R, invertibility, and rank checks on 200 random matrices.")

# %% [markdown]
# ## Part 2: Bases of the four fundamental subspaces
#
# * `column_space_basis`: the pivot columns of A itself (Theorem 3.28).
# * `null_space_basis`: the special solutions, eq. (3.2) (Theorem 3.18).
# * `row_space_basis`: the nonzero rows of R (Theorem 3.29).
# * `left_null_space_basis`: the last m - r rows of E (Theorem 3.31).
#
# Random matrices of known rank r are built as products B C with B an m x r and
# C an r x n Gaussian matrix (rank(BC) <= r by (3.7), and = r with probability 1).


# %%
def column_space_basis(A, tol=1e-10):
    """Basis of C(A) as the columns of an m x r array: the pivot columns of A (Theorem 3.28)."""
    A = np.asarray(A, float)
    _, piv = rref(A, tol)
    return A[:, piv]


def null_space_basis(A, tol=1e-10):
    """Basis of N(A) as the columns of an n x (n-r) array: the special solutions s_j,
    (s_j)_{f_j} = 1, other free entries 0, (s_j)_{p_i} = -R[i, f_j]  (eq. 3.2)."""
    A = np.asarray(A, float)
    m, n = A.shape
    R, piv = rref(A, tol)
    free = [j for j in range(n) if j not in piv]
    S = np.zeros((n, len(free)))
    for k, f in enumerate(free):
        S[f, k] = 1.0
        for i, p in enumerate(piv):
            S[p, k] = -R[i, f]
    return S


def row_space_basis(A, tol=1e-10):
    """Basis of C(A^T) as the columns of an n x r array: the nonzero rows of R (Theorem 3.29)."""
    R, piv = rref(A, tol)
    return R[:len(piv), :].T.copy()


def left_null_space_basis(A, tol=1e-10, partial_pivoting=True):
    """Basis of N(A^T) as the columns of an m x (m-r) array: rows r+1..m of E, EA = R (Theorem 3.31)."""
    R, piv, E = rref_with_E(A, tol, partial_pivoting)
    return E[len(piv):, :].T.copy()


def same_span(B1, B2, tol=1e-8):
    """True iff the columns of B1 and B2 span the same subspace:
    rank(B1) == rank(B2) == rank([B1 B2]). (matrix_rank is SVD-based; LA-9.)"""
    B1 = np.asarray(B1, float); B2 = np.asarray(B2, float)
    if B1.shape[1] == 0 or B2.shape[1] == 0:
        return B1.shape[1] == B2.shape[1]
    r1 = np.linalg.matrix_rank(B1, tol); r2 = np.linalg.matrix_rank(B2, tol)
    return r1 == r2 == np.linalg.matrix_rank(np.hstack([B1, B2]), tol)


def four_subspaces(A, tol=1e-10):
    """Dictionary of bases for C(A), N(A), C(A^T), N(A^T)."""
    return {"C(A)": column_space_basis(A, tol), "N(A)": null_space_basis(A, tol),
            "C(A^T)": row_space_basis(A, tol), "N(A^T)": left_null_space_basis(A, tol)}


n_checked = 0
for _ in range(300):
    m, n = rng.integers(1, 8, size=2)
    r = int(rng.integers(0, min(m, n) + 1))
    A = rng.standard_normal((m, r)) @ rng.standard_normal((r, n)) if r else np.zeros((m, n))
    S = four_subspaces(A)
    # Theorem 3.30: dimensions r, n - r, r, m - r.
    assert S["C(A)"].shape == (m, r) and S["N(A)"].shape == (n, n - r)
    assert S["C(A^T)"].shape == (n, r) and S["N(A^T)"].shape == (m, m - r)
    # Each basis actually lies in its subspace and is independent.
    assert np.allclose(A @ S["N(A)"], 0, atol=1e-8)
    assert np.allclose(A.T @ S["N(A^T)"], 0, atol=1e-8)
    for B in S.values():
        assert np.linalg.matrix_rank(B, 1e-8) == B.shape[1] if B.shape[1] else True
    # Eq. (3.6): N(A) ⟂ C(A^T) and N(A^T) ⟂ C(A).
    assert np.allclose(S["N(A)"].T @ S["C(A^T)"], 0, atol=1e-8)
    assert np.allclose(S["N(A^T)"].T @ S["C(A)"], 0, atol=1e-8)
    # Compare with SciPy's SVD-based answers: same subspaces, different bases.
    assert same_span(S["N(A)"], sla.null_space(A))
    assert same_span(S["C(A)"], sla.orth(A))
    assert same_span(S["C(A^T)"], sla.orth(A.T))
    assert same_span(S["N(A^T)"], sla.null_space(A.T))
    # Two routes to N(A^T): rows of E, or special solutions of A^T.
    assert same_span(S["N(A^T)"], null_space_basis(A.T))
    # Rank–Nullity (3.5) read off the shapes.
    assert S["C(A)"].shape[1] + S["N(A)"].shape[1] == n
    # Exercise B6: A = C R~, pivot columns of A times nonzero rows of R.
    assert np.allclose(S["C(A)"] @ S["C(A^T)"].T, A, atol=1e-8)
    n_checked += 1
print(f"Part 2: dimensions r, n-r, r, m-r, orthogonality (3.6), A = C R~, and agreement with SciPy on {n_checked} random matrices.")

# %% [markdown]
# ## Part 3: The worked examples, number by number
#
# Example 4.2 (3 x 4 matrix), 4.3 (rank one), 4.1 (independence in R^4),
# 4.4 (two bases of P_2), 4.5 (y'' + y = 0, with SymPy).


# %%
# ---- Example 4.2: all four subspaces of a 3 x 4 matrix ----
A42 = np.array([[1, 2, 0, 1],
                [2, 4, 1, 4],
                [1, 2, 1, 3]], float)
R, piv, E = rref_with_E(A42, partial_pivoting=False)      # the hand computation's path
assert np.allclose(R, [[1, 2, 0, 1], [0, 0, 1, 2], [0, 0, 0, 0]]) and piv == [0, 2]
assert np.allclose(E, [[1, 0, 0], [-2, 1, 0], [1, -1, 1]]), E
assert np.allclose(E @ A42, R)
S42 = four_subspaces(A42)
assert np.allclose(S42["C(A)"], [[1, 0], [2, 1], [1, 1]])                       # columns 1, 3 of A
assert np.allclose(S42["N(A)"], np.array([[-2, 1, 0, 0], [-1, 0, -2, 1]]).T)   # special solutions
assert np.allclose(S42["C(A^T)"], np.array([[1, 2, 0, 1], [0, 0, 1, 2]]).T)    # nonzero rows of R
y = left_null_space_basis(A42, partial_pivoting=False)[:, 0]
assert np.allclose(y, [1, -1, 1]), y
assert same_span(S42["N(A^T)"], np.array([[1, -1, 1]]).T)                       # pivoting changes E, not the span
assert np.allclose(y @ A42, 0)
# Free columns are the combinations R says they are (Theorem 3.28 step 4).
assert np.allclose(A42[:, 1], 2 * A42[:, 0]) and np.allclose(A42[:, 3], 1 * A42[:, 0] + 2 * A42[:, 2])
# Exercise B6(b): A = C R~ with C the pivot columns of A and R~ the nonzero rows of R.
C42, Rt42 = A42[:, piv], R[:len(piv), :]
assert np.allclose(C42, [[1, 0], [2, 1], [1, 1]]) and np.allclose(Rt42, [[1, 2, 0, 1], [0, 0, 1, 2]])
assert np.allclose(C42 @ Rt42, A42)
# Complete solution (3.3) for b = (1, 3, 2); inconsistency for b' = (1, 3, 3).
b = np.array([1, 3, 2.]); Eb = E @ b
assert np.allclose(Eb, [1, 1, 0]), Eb
x_p = np.zeros(4); x_p[piv] = Eb[:len(piv)]
assert np.allclose(x_p, [1, 0, 1, 0]) and np.allclose(A42 @ x_p, b)
for _ in range(20):                                        # x_p + N(A) all solve Ax = b ...
    c = rng.standard_normal(2)
    assert np.allclose(A42 @ (x_p + S42["N(A)"] @ c), b)
x_other = np.linalg.lstsq(A42, b, rcond=None)[0]           # ... and any solver's answer is x_p + something in N(A)
assert np.allclose(A42 @ x_other, b)
assert same_span(S42["N(A)"], np.column_stack([S42["N(A)"], x_other - x_p]))
b_bad = np.array([1, 3, 3.])
assert abs((E @ b_bad)[2]) > 1e-12 and np.isclose(y @ b_bad, 1.0)             # not in C(A)
print("Example 4.2: R, E, all four bases, A = C R~, x_p = (1,0,1,0), and the inconsistent b' reproduced.")

# ---- Example 4.3: rank-one matrix u v^T ----
u = np.array([1, 2, 3.]); v = np.array([1, -1, 2.])
A43 = np.outer(u, v)
R, piv, E = rref_with_E(A43, partial_pivoting=False)
assert np.allclose(R, [[1, -1, 2], [0, 0, 0], [0, 0, 0]]) and piv == [0]
assert np.allclose(E, [[1, 0, 0], [-2, 1, 0], [-3, 0, 1]])
S43 = four_subspaces(A43)
assert same_span(S43["C(A)"], u[:, None]) and same_span(S43["C(A^T)"], v[:, None])
assert np.allclose(S43["N(A)"], np.array([[1, 1, 0], [-2, 0, 1]]).T)
LN = left_null_space_basis(A43, partial_pivoting=False)
assert np.allclose(LN, np.array([[-2, 1, 0], [-3, 0, 1]]).T)
assert np.allclose(v @ S43["N(A)"], 0) and np.allclose(u @ LN, 0)            # N(A) = plane ⟂ v, N(A^T) = plane ⟂ u
assert [B.shape[1] for B in S43.values()] == [1, 2, 1, 2]
print("Example 4.3: rank-one dimensions 1, 2, 1, 2; N(A) ⟂ v and N(A^T) ⟂ u.")

# %% [markdown]
# ### Independence tester (Theorem 3.12) and Example 4.1
#
# Put the vectors as columns; they are independent iff every column is a pivot
# column. If not, a special solution is an explicit relation.


# %%
def is_independent(vectors, tol=1e-10):
    """(independent?, relation): vectors are the rows of the input (or a list of 1-D arrays).
    If dependent, `relation` c satisfies sum_j c_j v_j = 0 with c != 0 (a special solution)."""
    vectors = [np.asarray(v, float) for v in vectors]
    if len(vectors) == 0:
        return True, None                                        # the empty list is independent (Definition 3.11)
    V = np.array(vectors).T                                      # columns are the vectors
    N = null_space_basis(V, tol)
    if N.shape[1] == 0:
        return True, None
    return False, N[:, 0]


v1, v2, v3 = [1, 2, 0, 1], [0, 1, 1, 1], [2, 5, 1, 3]
indep, rel = is_independent([v1, v2, v3])
assert not indep and np.allclose(rel, [-2, -1, 1]), rel
assert np.allclose(np.array([v1, v2, v3]).T @ rel, 0)
w = [2, 5, 1, 4]
assert is_independent([v1, v2, w])[0]
print("Example 4.1: {v1, v2, v3} dependent with relation", rel, "; {v1, v2, w} independent.")

# Theorem 3.13(a): any list containing the zero vector is dependent.
for _ in range(20):
    k = int(rng.integers(1, 5)); n = int(rng.integers(k, 7))
    L = [rng.standard_normal(n) for _ in range(k)] + [np.zeros(n)]
    rng.shuffle(L)
    assert not is_independent(L)[0]
# Theorem 3.14: any n + 1 vectors in R^n are dependent; n random ones are (almost surely) independent.
for n in range(1, 8):
    assert not is_independent(rng.standard_normal((n + 1, n)))[0]
    assert is_independent(rng.standard_normal((n, n)))[0]
print("Theorems 3.13(a) and 3.14 confirmed on random lists.")

# %% [markdown]
# ### Example 4.4: two bases of P_2, coordinates converted
#
# Coordinates in (1, x, x^2) turn polynomials into vectors in R^3 (Theorem 3.17).
# M has the (1, x, x^2)-coordinates of 1, 1+x, 1+x+x^2 as columns; [p]_B = M [p]_C.


# %%
M = np.array([[1, 1, 1], [0, 1, 1], [0, 0, 1]], float)
assert is_independent(M.T)[0], "1, 1+x, 1+x+x^2 are independent (three pivots)"
p_B = np.array([2, 3, -1.])                     # p = 2 + 3x - x^2
p_C = np.linalg.solve(M, p_B)
assert np.allclose(p_C, [-1, 4, -1]), p_C
Minv = np.linalg.inv(M)
assert np.allclose(Minv, [[1, -1, 0], [0, 1, -1], [0, 0, 1]])
assert np.allclose(Minv @ p_B, p_C) and np.allclose(M @ p_C, p_B)
x = sp.symbols("x")
lhs = -1 * 1 + 4 * (1 + x) - 1 * (1 + x + x**2)
assert sp.expand(lhs - (2 + 3 * x - x**2)) == 0
print("Example 4.4: [p]_C =", p_C, "and -1 + 4(1+x) - (1+x+x^2) == 2 + 3x - x^2 symbolically.")

# %% [markdown]
# ### Example 4.5: the solution space of y'' + y = 0, symbolically
#
# cos t and sin t solve the ODE; their Wronskian is 1 (so they are independent,
# a second proof besides evaluating at t = 0 and t = pi/2); the energy identity
# E' = 2 z'(z'' + z) holds for an *arbitrary* smooth z, which is what the
# spanning proof uses (so E is constant along every solution, not only along
# combinations of cos and sin); the coordinates of a cos t + b sin t are its
# initial conditions (a, b); and SymPy's general solution of the forced
# equation y'' + y = t is particular + homogeneous, i.e. eq. (3.3).


# %%
t = sp.symbols("t", real=True)
a, b_ = sp.symbols("a b", real=True)
L = lambda f: sp.diff(f, t, 2) + f                       # the linear operator L[y] = y'' + y
assert sp.simplify(L(sp.cos(t))) == 0 and sp.simplify(L(sp.sin(t))) == 0
W = sp.Matrix([[sp.cos(t), sp.sin(t)], [sp.diff(sp.cos(t), t), sp.diff(sp.sin(t), t)]]).det()
assert sp.simplify(W) == 1
# Energy identity for an arbitrary function z(t): E' = 2 z' (z'' + z) = 2 z' L[z].
zf = sp.Function("z")(t)
Eng = zf.diff(t)**2 + zf**2
assert sp.simplify(Eng.diff(t) - 2 * zf.diff(t) * L(zf)) == 0
# Coordinates of a solution are its initial conditions: (z(0), z'(0)) = (a, b).
z = a * sp.cos(t) + b_ * sp.sin(t)
assert sp.simplify(z.subs(t, 0) - a) == 0 and sp.simplify(sp.diff(z, t).subs(t, 0) - b_) == 0
# Forced oscillator: general solution = t + homogeneous (Theorem 3.9 for the operator L).
yf = sp.Function("y")
sol = sp.dsolve(sp.Eq(yf(t).diff(t, 2) + yf(t), t), yf(t)).rhs
assert sp.simplify(L(sol) - t) == 0
consts = sorted(sol.free_symbols - {t}, key=str)
assert len(consts) == 2 and sp.simplify(sol.subs({c: 0 for c in consts}) - t) == 0
assert sp.simplify(L(sol - t)) == 0                       # the rest is in N(L) = span{cos, sin}
print("Example 4.5: cos, sin solve the ODE, W = 1, energy identity E' = 2z'(z''+z); dsolve(y''+y=t) =", sol)

# %% [markdown]
# ## Part 4: Pictures of N(A) and C(A^T) in R^3
#
# For a 1 x 3 matrix, N(A) is the plane perpendicular to the single row and
# C(A^T) is the line along it. For a 2 x 3 matrix of rank 2, N(A) is a line and
# C(A^T) is the plane spanned by the rows. In both pictures the two subspaces
# are perpendicular (eq. 3.6) and their dimensions add to n = 3.


# %%
def draw_subspace(ax, B, color, label, extent=2.0):
    """Draw span of the columns of B (3 x k, k in {1, 2}) as a line segment or a plane patch."""
    if B.shape[1] == 1:
        d = B[:, 0] / np.linalg.norm(B[:, 0])
        pts = np.outer(np.array([-extent, extent]), d)
        ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=color, lw=2.5, label=label)
    elif B.shape[1] == 2:
        s = np.linspace(-extent, extent, 9)
        S1, S2 = np.meshgrid(s, s)
        Q, _ = np.linalg.qr(B)                              # orthonormal basis of the same plane, for a square patch
        P = S1[..., None] * Q[:, 0] + S2[..., None] * Q[:, 1]
        ax.plot_surface(P[..., 0], P[..., 1], P[..., 2], color=color, alpha=0.35, linewidth=0)
        ax.plot([], [], [], color=color, lw=6, alpha=0.5, label=label)


fig = plt.figure(figsize=(11, 5))
for k, (A, title) in enumerate([(np.array([[1, -1, 2.]]), "A is 1 x 3: N(A) plane, C(A^T) line"),
                                (np.array([[1, 0, 1.], [0, 1, 1.]]), "A is 2 x 3: N(A) line, C(A^T) plane")], 1):
    N, Rw = null_space_basis(A), row_space_basis(A)
    assert np.allclose(N.T @ Rw, 0)                          # (3.6)
    assert N.shape[1] + Rw.shape[1] == 3                     # (3.5)
    ax = fig.add_subplot(1, 2, k, projection="3d")
    draw_subspace(ax, N, "C0", f"N(A), dim {N.shape[1]}")
    draw_subspace(ax, Rw, "C3", f"C(A^T), dim {Rw.shape[1]}")
    ax.quiver(0, 0, 0, *(A.T / np.abs(A).max()).T.reshape(-1, 3).T, color="k", arrow_length_ratio=0.15)
    ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.set_zlim(-2, 2)
    ax.set_title(title, fontsize=10); ax.legend(fontsize=8, loc="upper left")
fig.suptitle("Null space and row space are perpendicular; dimensions add to n = 3")
fig.tight_layout()
fig.savefig(FIG / "null_vs_row_space.png", dpi=120)
plt.close(fig)
print("Part 4: figure saved to", FIG / "null_vs_row_space.png")

# %% [markdown]
# ## Part 5: Extending to a basis, and the rank of a product
#
# `extend_to_basis` is the constructive version of Theorem 3.25: walk through
# the standard basis and keep e_i iff it is not in the current span
# (Lemma 3.20 keeps the list independent). `in_span` is LA-1 Theorem 3.6 plus
# LA-2 Theorem 3.7(a): v is in span(S) iff [S | v] has no pivot in the last column.


# %%
def in_span(v, S, tol=1e-10):
    """Is v in the span of the columns of S? Pivot test on the augmented matrix [S | v]."""
    S = np.asarray(S, float).reshape(len(v), -1)
    _, piv = rref(np.column_stack([S, np.asarray(v, float)]), tol)
    return S.shape[1] not in piv


def extend_to_basis(vectors, n=None, tol=1e-10):
    """Extend an independent list (rows of `vectors`, or list of 1-D arrays) to a basis of R^n,
    adjoining standard basis vectors (Theorem 3.25). Returns an n x n array whose columns are the basis.
    `n` is needed only when the list is empty."""
    vectors = [np.asarray(v, float) for v in vectors]
    if n is None:
        n = len(vectors[0])
    U = np.array(vectors).T if vectors else np.zeros((n, 0))
    assert is_independent(U.T, tol)[0], "input list must be independent"
    for i in range(n):
        e = np.zeros(n); e[i] = 1.0
        if not in_span(e, U, tol):
            U = np.column_stack([U, e])
    assert U.shape == (n, n) and len(rref(U, tol)[1]) == n
    return U


B4 = extend_to_basis([[1, 2, 0, 1], [0, 1, 1, 1], [2, 5, 1, 4]])       # Example 4.1's independent triple
assert B4.shape == (4, 4)
print("Example 4.1 extended to a basis of R^4 by adjoining", B4[:, 3])
for _ in range(50):
    n = int(rng.integers(1, 7)); k = int(rng.integers(0, n + 1))
    U = rng.standard_normal((k, n))
    B = extend_to_basis(U, n=n)
    assert np.allclose(B[:, :k], U.T) and len(rref(B)[1]) == n           # original vectors kept, result is a basis
print("extend_to_basis: 50 random independent lists extended to bases.")


def rank(A, tol=1e-10):
    """rank = number of pivots = dim C(A) = dim C(A^T) (Theorems 3.28, 3.29)."""
    return len(rref(A, tol)[1])


# Theorem 3.32: rank(AB) <= min(rank A, rank B), on random products of known ranks.
strict = 0
for _ in range(200):
    m, p, n = rng.integers(1, 7, size=3)
    ra, rb = int(rng.integers(0, min(m, p) + 1)), int(rng.integers(0, min(p, n) + 1))
    A = rng.standard_normal((m, ra)) @ rng.standard_normal((ra, p)) if ra else np.zeros((m, p))
    B = rng.standard_normal((p, rb)) @ rng.standard_normal((rb, n)) if rb else np.zeros((p, n))
    rAB = rank(A @ B)
    assert rAB <= min(ra, rb) and rAB == np.linalg.matrix_rank(A @ B)
    # Step 1 and 2 of the proof: C(AB) ⊆ C(A) and N(B) ⊆ N(AB).
    assert all(in_span(col, column_space_basis(A)) for col in (A @ B).T) or ra == 0
    NB = null_space_basis(B)
    assert np.allclose(A @ B @ NB, 0, atol=1e-8)
    strict += rAB < min(ra, rb)
A_s = np.array([[1, 0], [0, 0.]]); B_s = np.array([[0, 0], [0, 1.]])
assert rank(A_s) == rank(B_s) == 1 and rank(A_s @ B_s) == 0
# Random products are generically equality cases. Strictness needs structure: A = u v^T, B = w z^T
# with w ⟂ v gives AB = u (v·w) z^T = 0, so rank(AB) = 0 < 1. (C(B) = span{w} sits inside N(A).)
for _ in range(50):
    u, v, z = rng.standard_normal((3, 5)); w = rng.standard_normal(5)
    w -= (w @ v) / (v @ v) * v                              # make w orthogonal to v (LA-1 Example 4.5)
    assert rank(np.outer(u, v) @ np.outer(w, z)) == 0 < 1 == rank(np.outer(u, v)) == rank(np.outer(w, z))
print(f"Theorem 3.32 holds on 200 random products ({strict} strict) and is strict on 50 constructed ones.")

# %% [markdown]
# ## Your turn
#
# 1. **Symmetric 2 x 2 matrices.** Flatten a 2 x 2 matrix to (a, b, c, d). The
#    symmetric ones are the null space of the 1 x 4 matrix [0, 1, -1, 0]
#    (the equation b - c = 0). Use `null_space_basis` to get a basis, read off
#    the dimension, and check it spans the same space as the vectorized
#    {E11, E12 + E21, E22}. Then do 3 x 3 (dimension 6) the same way.
# 2. **Polynomials of degree <= 3 vanishing at x = 1.** In coefficient
#    coordinates (a0, a1, a2, a3), the condition p(1) = 0 is the 1 x 4 matrix
#    [1, 1, 1, 1]. Show the subspace has dimension 3 and compare the special
#    solutions with the basis {x - 1, x^2 - 1, x^3 - 1}. Why is "p(1) = 1"
#    NOT a subspace? (Try to find it as a null space; you cannot.)
# 3. **Extend {(1,1,0), (0,1,1)} to a basis of R^3** with `extend_to_basis`.
#    Which standard basis vector gets adjoined, and why not e_2 or e_3? Then
#    adjoin a *random* vector instead and confirm you still get a basis.
# 4. (B8) Write `basis_of_sum(U, W)` and `basis_of_intersection(U, W)` for two
#    subspaces given by bases, and check the dimension formula (3.8) of
#    Theorem 3.35 on random pairs of subspaces of R^6.
#
# Reference answers for 1–3 follow so you can compare.

# %%
# 1. Symmetric 2 x 2: dimension 3.
Sym2 = null_space_basis(np.array([[0, 1, -1, 0.]]))
assert Sym2.shape == (4, 3)
E11, E12p21, E22 = np.array([1, 0, 0, 0.]), np.array([0, 1, 1, 0.]), np.array([0, 0, 0, 1.])
assert same_span(Sym2, np.column_stack([E11, E12p21, E22]))
# 3 x 3 symmetric: three equations a12 = a21, a13 = a31, a23 = a32 on 9 unknowns -> dimension 6.
C9 = np.zeros((3, 9)); C9[0, 1] = 1; C9[0, 3] = -1; C9[1, 2] = 1; C9[1, 6] = -1; C9[2, 5] = 1; C9[2, 7] = -1
assert null_space_basis(C9).shape == (9, 6)
# Skew 3 x 3 (A^T = -A: zero diagonal, a_ij = -a_ji) has dimension 3, and M_3x3 = Sym_3 ⊕ Skew_3
# (exercise B3(d), Theorem 3.34): the two bases together are 6 + 3 = 9 independent vectors in R^9.
K9 = np.zeros((6, 9)); K9[[0, 1, 2], [0, 4, 8]] = 1
K9[3, 1] = K9[3, 3] = 1; K9[4, 2] = K9[4, 6] = 1; K9[5, 5] = K9[5, 7] = 1
assert null_space_basis(K9).shape == (9, 3)
assert rank(np.hstack([null_space_basis(C9), null_space_basis(K9)])) == 9
# 2. Cubics vanishing at 1: dimension 3, same span as {x - 1, x^2 - 1, x^3 - 1}.
V1 = null_space_basis(np.array([[1, 1, 1, 1.]]))
assert V1.shape == (4, 3)
alt = np.array([[-1, 1, 0, 0], [-1, 0, 1, 0], [-1, 0, 0, 1.]]).T   # coefficients of x-1, x^2-1, x^3-1
assert same_span(V1, alt)
# 3. Extend {(1,1,0), (0,1,1)}: e_1 is adjoined (it is not in the span; check: a(1,1,0)+b(0,1,1)=(1,0,0) forces a=1, b=-1, b=0).
B3 = extend_to_basis([[1, 1, 0], [0, 1, 1]])
assert np.allclose(B3[:, 2], [1, 0, 0]) and len(rref(B3)[1]) == 3
assert not in_span([1, 0, 0], B3[:, :2]) and in_span([1, 2, 1], B3[:, :2])
print("Your turn reference: dim Sym(2x2) = 3, dim Sym(3x3) = 6, dim Skew(3x3) = 3, dim {p in P_3 : p(1) = 0} = 3, basis extended by", B3[:, 2])

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
