# %% [markdown]
# # LA-2 Project: Elimination, LU, Inverses, and Why Pivoting Matters
#
# This script builds the algorithms of LA-2 from scratch and checks every number
# in the notes. Parts: (1) elementary matrices and left-vs-right multiplication,
# (2) LU with partial pivoting, triangular substitution, and `solve`, checked
# against SciPy/NumPy and against Examples 4.1 and 4.4, (3) RREF with a
# tolerance and the three outcomes of Theorem 3.7, (4) Gauss–Jordan inverses and
# the case against computing them, (5) the 1e-20 pivoting catastrophe,
# (6) operation counts and timing, (7) checks of the exercise numbers.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
import numpy as np
import scipy.linalg as sla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)

# %% [markdown]
# ## Part 1: Elementary matrices, eq. (3.1)
#
# E_ij(l) = I - l e_i e_j^T subtracts l * row j from row i.  P_ij swaps rows i, j.
# D_i(c) scales row i by c.  Theorem 3.1: EA acts on rows, AE acts on columns.
# Indices here are 0-based (the notes are 1-based).


# %%
def unit(n, i, j):
    """The matrix e_i e_j^T: a single 1 in position (i, j)."""
    M = np.zeros((n, n)); M[i, j] = 1.0; return M


def E(n, i, j, l):
    """E_ij(l) = I - l e_i e_j^T, eq. (3.1): row i <- row i - l * row j."""
    return np.eye(n) - l * unit(n, i, j)


def Pswap(n, i, j):
    """P_ij from (3.1): the identity with rows i and j exchanged."""
    return np.eye(n) - unit(n, i, i) - unit(n, j, j) + unit(n, i, j) + unit(n, j, i)


def D(n, i, c):
    """D_i(c) from (3.1): the identity with c in position (i, i)."""
    return np.eye(n) + (c - 1) * unit(n, i, i)


A = rng.standard_normal((4, 4))
for (i, j) in [(1, 0), (3, 1), (0, 2)]:
    l = rng.standard_normal()
    EA = E(4, i, j, l) @ A
    expect = A.copy(); expect[i] -= l * A[j]                 # (R1) on rows
    assert np.allclose(EA, expect), "left multiplication must act on rows"
    AE = A @ E(4, i, j, l)
    expect = A.copy(); expect[:, j] -= l * A[:, i]           # right multiplication: column j <- column j - l column i
    assert np.allclose(AE, expect), "right multiplication must act on columns"
    # Inverses derived in §3.1.
    assert np.allclose(E(4, i, j, l) @ E(4, i, j, -l), np.eye(4))
    assert np.allclose(Pswap(4, i, j) @ Pswap(4, i, j), np.eye(4))
    assert np.allclose(D(4, i, 3.0) @ D(4, i, 1 / 3.0), np.eye(4))
    assert np.allclose(Pswap(4, i, j) @ A, A[[j if k == i else i if k == j else k for k in range(4)]])

# Example 4.1: the forward product E32 E31 E21 is NOT "minus L" (§3.5, step 4).
l21, l31, l32 = 2.0, -1.0, -3.0
Efwd = E(3, 2, 1, l32) @ E(3, 2, 0, l31) @ E(3, 1, 0, l21)
assert np.allclose(Efwd, [[1, 0, 0], [-2, 1, 0], [-5, 3, 1]]), Efwd
assert np.isclose(Efwd[2, 0], -l31 + l32 * l21)              # the surviving cross term
L_41 = np.linalg.inv(Efwd)
assert np.allclose(L_41, [[1, 0, 0], [2, 1, 0], [-1, -3, 1]])  # L = (E32 E31 E21)^(-1) holds +l_ij
print("Part 1: elementary matrices act on rows (left) / columns (right); inverses verified;")
print("        Example 4.1: (E32 E31 E21)[3,1] =", Efwd[2, 0], " but L[3,1] =", L_41[2, 0])

# %% [markdown]
# ## Part 2: PA = LU with partial pivoting, substitution, solve
#
# `lu_decompose` implements Algorithm 3.5 on a square matrix and returns P, L, U
# with P A = L U, eq. (3.6). The swap of the already-computed columns of L is
# Theorem 3.13. `forward_sub` is (3.4), `back_sub` is (3.3). Flops are counted
# as in §3.7: at step k (0-based), m = n-k-1 rows below, m divisions and 2m^2
# for the rank-one update of the trailing block.
#
# **Build it yourself first:** write these four functions before reading them.


# %%
def lu_decompose(A, pivot=True, tol=0.0, count_flops=False):
    """P A = L U (eq. 3.6) by elimination with partial pivoting (Algorithm 3.5 + Theorem 3.13).

    Returns P (permutation), L (unit lower, L[i,j] = multiplier l_ij), U (upper).
    Raises LinAlgError if a pivot has |u_kk| <= tol (singular to that tolerance).
    With count_flops=True also returns the flop count of §3.7.
    """
    U = np.array(A, dtype=float)
    n, n2 = U.shape
    assert n == n2, "square matrices only; use rref() for rectangular systems"
    L = np.eye(n)
    perm = np.arange(n)                    # perm[i] = which original row now sits in row i
    flops = 0
    for k in range(n):
        if pivot:
            p = k + int(np.argmax(np.abs(U[k:, k])))
            if p != k:
                U[[k, p], :] = U[[p, k], :]
                L[[k, p], :k] = L[[p, k], :k]   # Theorem 3.13: swap earlier multipliers too
                perm[[k, p]] = perm[[p, k]]
        if abs(U[k, k]) <= tol:
            raise np.linalg.LinAlgError(f"pivot {k} is {U[k, k]:.3g}: singular to tolerance {tol}")
        if k == n - 1:
            break
        l = U[k + 1:, k] / U[k, k]              # multipliers l_ik (eq. 3.5)
        L[k + 1:, k] = l
        U[k + 1:, k + 1:] -= np.outer(l, U[k, k + 1:])   # row i <- row i - l_ik row k, columns k+1..n
        U[k + 1:, k] = 0.0                      # known to be zero; not computed (§3.7)
        m = n - k - 1
        flops += m + 2 * m * m
    P = np.eye(n)[perm]                         # row i of P is e_perm[i]^T, so (P A)[i] = A[perm[i]]
    return (P, L, U, flops) if count_flops else (P, L, U)


def forward_sub(L, b):
    """Solve L c = b for unit lower triangular L, eq. (3.4)."""
    n = len(b); c = np.zeros(n)
    for i in range(n):
        c[i] = b[i] - L[i, :i] @ c[:i]
    return c


def back_sub(U, c):
    """Solve U x = c for upper triangular U with nonzero diagonal, eq. (3.3)."""
    n = len(c); x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (c[i] - U[i, i + 1:] @ x[i + 1:]) / U[i, i]
    return x


def solve(A, b, pivot=True, tol=0.0):
    """Solve A x = b: factor P A = L U, then L c = P b, then U x = c."""
    P, L, U = lu_decompose(A, pivot=pivot, tol=tol)
    return back_sub(U, forward_sub(L, P @ np.asarray(b, float)))


# Random checks against SciPy (which returns A = P_s L U, i.e. P_s = P^T) and NumPy.
for n in (2, 3, 5, 8, 20, 50):
    for _ in range(3):
        A = rng.standard_normal((n, n)); b = rng.standard_normal(n)
        P, L, U = lu_decompose(A)
        assert np.allclose(P @ A, L @ U), "PA = LU failed"
        assert np.allclose(P @ P.T, np.eye(n)) and np.allclose(P.T, np.linalg.inv(P))   # Theorem 3.23
        assert np.allclose(np.tril(L, -1) + np.eye(n), L) and np.allclose(np.triu(U), U)
        assert np.abs(L).max() <= 1 + 1e-12, "partial pivoting must give |l_ij| <= 1"
        Ps, Ls, Us = sla.lu(A)
        assert np.allclose(P, Ps.T) and np.allclose(L, Ls) and np.allclose(U, Us)      # Theorem 3.12 (uniqueness)
        x = solve(A, b)
        assert np.allclose(x, np.linalg.solve(A, b), atol=1e-9)
        # LA-6 preview: det A = (sign of P) * product of pivots.
        assert np.isclose(np.linalg.det(A), np.linalg.det(P) * np.prod(np.diag(U)))
print("Part 2: PA = LU matches scipy.linalg.lu and solve matches np.linalg.solve for n up to 50.")

# Example 4.1 to the digit.
A41 = np.array([[2, 1, 1], [4, 1, 0], [-2, 2, 1]], float); b41 = np.array([3, 6, 1.])
P, L, U = lu_decompose(A41, pivot=False)
assert np.allclose(L, [[1, 0, 0], [2, 1, 0], [-1, -3, 1]]) and np.allclose(U, [[2, 1, 1], [0, -1, -2], [0, 0, -4]])
c = forward_sub(L, b41); assert np.allclose(c, [3, 0, 4])
x = back_sub(U, c); assert np.allclose(x, [1, 2, -1]) and np.allclose(A41 @ x, b41)
assert np.allclose(np.linalg.inv(Efwd), L)          # L is the inverse of the forward product
print("Example 4.1: pivots", np.diag(U), " c =", c, " x =", x)

# Example 4.4: the circuit. Pivots are effective conductances (Kron reduction).
R1, R2, R3, R4 = 1.0, 2.0, 2.0, 1.0
g1, g2, g3, g4 = 1 / R1, 1 / R2, 1 / R3, 1 / R4
G = np.array([[g1 + g2, -g2, 0], [-g2, g2 + g3, -g3], [0, -g3, g3 + g4]])
i_src = np.array([3.0, 0, 0])
P, L, U = lu_decompose(G, pivot=False)
assert np.allclose(np.diag(U), [3 / 2, 5 / 6, 6 / 5])
assert np.allclose(L[1, 0], -1 / 3) and np.allclose(L[2, 1], -3 / 5)
assert np.isclose(U[1, 1], 1 / (R1 + R2) + g3), "pivot 2 = (R1+R2 in series) parallel with R3"
assert np.isclose(U[2, 2], 1 / (R1 + R2 + R3) + g4), "pivot 3 = (R1+R2+R3 in series) parallel with R4"
v = solve(G, i_src, pivot=False)
assert np.allclose(v, [2.5, 1.5, 0.5])
assert np.isclose(v[0] / R1 + (v[0] - v[1]) / R2, 3.0)   # KCL at node 1
assert np.isclose(np.linalg.det(G), 1.5)
print("Example 4.4: v =", v, "V; pivots =", np.diag(U), "= effective conductances")

# %% [markdown]
# ## Part 3: RREF with a tolerance, and the three outcomes (Theorem 3.7)
#
# `rref` returns the reduced row echelon form and the pivot columns, treating
# any candidate pivot below `tol` as zero. `classify` applies it to [A | b]:
# 'none' if the augmented column is a pivot column, 'unique' if every column of
# A is a pivot column, else 'infinite' with x_p and the special solutions s_f.


# %%
def rref(A, tol=1e-10):
    """Reduced row echelon form (Definition 3.4) and the list of pivot columns."""
    R = np.array(A, dtype=float); m, n = R.shape
    pivots = []; r = 0
    for j in range(n):
        if r == m:
            break
        p = r + int(np.argmax(np.abs(R[r:, j])))
        if abs(R[p, j]) < tol:
            R[r:, j] = 0.0                       # numerically zero column below the staircase: free
            continue
        R[[r, p]] = R[[p, r]]
        R[r] = R[r] / R[r, j]                    # scale the pivot to 1
        others = [i for i in range(m) if i != r]
        R[others] -= np.outer(R[others, j], R[r])   # clear above and below
        pivots.append(j); r += 1
    return R, pivots


def classify(A, b, tol=1e-10):
    """Theorem 3.7 for A x = b: kind in {'none', 'unique', 'infinite'}, pivot/free columns,
    a particular solution x_p (free variables 0) and the special solutions s_f."""
    A = np.asarray(A, float); b = np.asarray(b, float).reshape(-1); m, n = A.shape
    R, piv = rref(np.column_stack([A, b]), tol)
    if n in piv:
        return {"kind": "none", "pivot_cols": [j for j in piv if j < n], "free_cols": None}
    free = [j for j in range(n) if j not in piv]
    x_p = np.zeros(n)
    for k, j in enumerate(piv):
        x_p[j] = R[k, n]
    specials = []
    for f in free:
        s = np.zeros(n); s[f] = 1.0
        for k, j in enumerate(piv):
            s[j] = -R[k, f]
        specials.append(s)
    return {"kind": "unique" if not free else "infinite", "pivot_cols": piv, "free_cols": free,
            "x_p": x_p, "specials": specials}


def rank(A, tol=1e-10):
    """Provisional rank (Definition 3.20): number of pivots above tol."""
    return len(rref(A, tol)[1])


# Example 4.2: singular 3x3.
A42 = np.array([[1, 2, 3], [2, 4, 7], [3, 6, 10]], float)
R42, piv42 = rref(A42)
assert np.allclose(R42, [[1, 2, 0], [0, 0, 1], [0, 0, 0]]) and piv42 == [0, 2] and rank(A42) == 2
out = classify(A42, [1, 3, 4])
assert out["kind"] == "infinite" and out["free_cols"] == [1]
assert np.allclose(out["x_p"], [-2, 0, 1]) and np.allclose(out["specials"][0], [-2, 1, 0])
assert np.allclose(A42 @ out["x_p"], [1, 3, 4]) and np.allclose(A42 @ out["specials"][0], 0)
assert classify(A42, [1, 3, 5])["kind"] == "none"
try:
    lu_decompose(A42, tol=1e-12)            # no usable pivot in column 2 -> singular
    raise AssertionError("should have raised")
except np.linalg.LinAlgError as err:
    print("Example 4.2: lu_decompose raised:", err)
print("Example 4.2: RREF pivots", piv42, "| b=(1,3,4):", out["kind"], "x_p =", out["x_p"],
      "s =", out["specials"][0], "| b=(1,3,5):", classify(A42, [1, 3, 5])["kind"])

# LA-1's three 2x2 outcomes (LA-1 exercise B1), now classified mechanically.
assert classify([[1, 0], [0, 1]], [1, 2])["kind"] == "unique"
assert classify([[1, 1], [1, 1]], [1, 2])["kind"] == "none"
assert classify([[1, 1], [1, 1]], [1, 1])["kind"] == "infinite"

# Theorem 3.7(c) on random rectangular systems of known rank: A x_p = b, A s_f = 0, n - r specials.
for _ in range(20):
    m, n, r = rng.integers(2, 7), rng.integers(2, 7), 0
    r = rng.integers(1, min(m, n) + 1)
    A = rng.standard_normal((m, r)) @ rng.standard_normal((r, n))      # exact rank r
    b = A @ rng.standard_normal(n)                                     # consistent by construction
    out = classify(A, b, tol=1e-8)
    assert out["kind"] in ("unique", "infinite") and len(out["pivot_cols"]) == r == np.linalg.matrix_rank(A)
    assert np.allclose(A @ out["x_p"], b, atol=1e-8) and len(out["specials"]) == n - r
    for s in out["specials"]:
        assert np.allclose(A @ s, 0, atol=1e-8)
    if m > r:   # a right-hand side off the column space is inconsistent
        b_bad = b + 1e3 * rng.standard_normal(m)
        assert classify(A, b_bad, tol=1e-8)["kind"] == "none" or np.linalg.matrix_rank(np.column_stack([A, b_bad])) == r
print("Part 3: Theorem 3.7 verified on random systems; rank (pivot count) matches np.linalg.matrix_rank.")

# %% [markdown]
# ## Part 4: Gauss–Jordan inverses, Theorem 3.19, and the case against inv()
#
# [A | I] -> [I | A^{-1}], eq. (3.8). Then the three arguments of §3.7:
# cost (inv costs ~3 LU's and saves nothing per solve), structure (the inverse
# of a sparse matrix is dense), and accuracy (printed, not asserted: look).


# %%
def gauss_jordan_inverse(A, tol=1e-12):
    """A^{-1} by row-reducing [A | I] to [I | A^{-1}] (Corollary 3.17, eq. 3.8)."""
    A = np.asarray(A, float); n = A.shape[0]
    R, piv = rref(np.column_stack([A, np.eye(n)]), tol)
    if piv[:n] != list(range(n)):
        raise np.linalg.LinAlgError("fewer than n pivots: singular")
    return R[:, n:]


A43 = np.array([[2, 1, 0], [1, 2, 1], [0, 1, 2]], float)
Ainv = gauss_jordan_inverse(A43)
assert np.allclose(Ainv, np.array([[3, -2, 1], [-2, 4, -2], [1, -2, 3]]) / 4)
assert np.allclose(A43 @ Ainv, np.eye(3)) and np.allclose(Ainv, Ainv.T)
P, L, U = lu_decompose(A43, pivot=False)
assert np.allclose(np.diag(U), [2, 3 / 2, 4 / 3]) and np.isclose(np.prod(np.diag(U)), 4) and np.isclose(np.linalg.det(A43), 4)
assert np.allclose(L, [[1, 0, 0], [1 / 2, 1, 0], [0, 2 / 3, 1]])
Dm = np.diag(np.diag(U))
assert np.allclose(L @ Dm @ L.T, A43), "Theorem 3.22: symmetric A = L D L^T"
assert np.allclose(np.linalg.inv(Dm) @ U, L.T), "U' = D^{-1} U equals L^T"
print("Example 4.3: A^{-1} =\n", Ainv, "\n   pivots", np.diag(U), " product", np.prod(np.diag(U)), " (= det A)")

# Theorem 3.19 on random matrices.
for _ in range(10):
    A = rng.standard_normal((5, 5)); B = rng.standard_normal((5, 5))
    Ai, Bi = gauss_jordan_inverse(A), gauss_jordan_inverse(B)
    assert np.allclose(Ai, np.linalg.inv(A))
    assert np.allclose(gauss_jordan_inverse(A @ B), Bi @ Ai)             # (AB)^{-1} = B^{-1} A^{-1}
    assert not np.allclose(gauss_jordan_inverse(A @ B), Ai @ Bi)        # ... and the order matters
    assert np.allclose(gauss_jordan_inverse(A.T), Ai.T)                 # (A^T)^{-1} = (A^{-1})^T
    assert np.allclose(gauss_jordan_inverse(Ai), A)                     # (A^{-1})^{-1} = A
try:
    gauss_jordan_inverse(A42); raise AssertionError
except np.linalg.LinAlgError:
    pass
print("Part 4: Gauss–Jordan matches np.linalg.inv; Theorem 3.19 verified; singular A42 rejected.")

# The case against inv(): (1) structure. Tridiagonal A has 3n-2 nonzeros; its inverse has n^2.
n = 50
T = 2 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
P, L, U = lu_decompose(T, pivot=False)
nnz = lambda M: int(np.sum(np.abs(M) > 1e-14))
assert nnz(T) == 3 * n - 2 and nnz(L) + nnz(U) - n == 3 * n - 2     # L and U are bidiagonal: same storage as T
assert nnz(np.linalg.inv(T)) == n * n                               # the inverse is completely dense
print(f"Never invert (structure): tridiagonal {n}x{n} has {nnz(T)} nonzeros, L and U together {nnz(L)+nnz(U)-n}, inverse {nnz(np.linalg.inv(T))}.")

# (2) cost: time inv vs. a factor-and-solve (LAPACK for both, so the comparison is fair).
n = 1000
A = rng.standard_normal((n, n)); b = rng.standard_normal(n)
t0 = time.perf_counter(); Ai = np.linalg.inv(A); x_inv = Ai @ b; t_inv = time.perf_counter() - t0
t0 = time.perf_counter(); x_sol = np.linalg.solve(A, b); t_sol = time.perf_counter() - t0
print(f"Never invert (cost): n={n}: inv(A)@b took {t_inv:.3f}s, solve(A,b) took {t_sol:.3f}s (ratio {t_inv/t_sol:.1f}; §3.7 predicts about 3).")

# (3) accuracy: the relative residual of both routes on an ill-conditioned matrix. Printed honestly.
Q1, _ = np.linalg.qr(rng.standard_normal((200, 200))); Q2, _ = np.linalg.qr(rng.standard_normal((200, 200)))
A_ill = Q1 @ np.diag(np.logspace(0, -12, 200)) @ Q2.T; b = rng.standard_normal(200)
x_sol = np.linalg.solve(A_ill, b); x_inv = np.linalg.inv(A_ill) @ b
rr = lambda x: np.linalg.norm(A_ill @ x - b) / (np.linalg.norm(A_ill, 2) * np.linalg.norm(x))
print(f"Never invert (accuracy): kappa={np.linalg.cond(A_ill):.1e}; relative residual via solve {rr(x_sol):.1e}, via inv {rr(x_inv):.1e}.")
print("        (Backward stability guarantees the first is at rounding level; the second has no such guarantee, though it is often close.)")

# %% [markdown]
# ## Part 5: The 1e-20 example (§3.6, Example 4.5)
#
# A = [[eps, 1], [1, 1]], b = (1, 2). Exact solution x1 = 1/(1-eps), x2 = (1-2eps)/(1-eps).
# Without pivoting the multiplier is 1/eps and fl(1 - 1/eps) loses the 1.


# %%
def pivot_experiment(eps):
    A = np.array([[eps, 1.0], [1.0, 1.0]]); b = np.array([1.0, 2.0])
    x_true = np.array([1 / (1 - eps), (1 - 2 * eps) / (1 - eps)])
    x_np = solve(A, b, pivot=False); x_p = solve(A, b, pivot=True)
    return np.abs(x_np - x_true).max(), np.abs(x_p - x_true).max()


eps = 1e-20
A_eps = np.array([[eps, 1.0], [1.0, 1.0]])
P, L, U = lu_decompose(A_eps, pivot=False)
assert L[1, 0] == 1e20 and U[1, 1] == -1e20, (L, U)                  # fl(1 - 1e20) = -1e20
assert (L @ U)[1, 1] == 0.0, "the computed LU has a 0 where A has a 1"
x_bad = solve(A_eps, [1, 2], pivot=False); x_good = solve(A_eps, [1, 2], pivot=True)
assert np.allclose(x_bad, [0, 1]) and np.allclose(x_good, [1, 1])
print("Example 4.5: without pivoting LU =", (L @ U).tolist(), " x =", x_bad, "; with pivoting x =", x_good)

exps = np.arange(1, 21)
errs = np.array([pivot_experiment(10.0 ** -k) for k in exps])
assert errs[-1, 0] > 0.5 and errs[:, 1].max() < 1e-12
fig, ax = plt.subplots(figsize=(6.4, 4))
ax.semilogy(exps, np.maximum(errs[:, 0], 1e-18), "o-", lw=2, ms=6, color="C1", label="no pivoting")
ax.semilogy(exps, np.maximum(errs[:, 1], 1e-18), "s-", lw=2, ms=6, color="C0", label="partial pivoting")
ax.axvline(16, color="0.6", lw=1, ls="--"); ax.text(16.2, 1e-9, "machine\nprecision", fontsize=8, color="0.4")
ax.set_xlabel(r"$k$ in $\varepsilon = 10^{-k}$"); ax.set_ylabel(r"max error in $x$ (floored at $10^{-18}$)")
ax.set_title(r"Elimination on $[[\varepsilon, 1],[1, 1]]\,x = (1, 2)$"); ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "pivoting_error.png", dpi=120); plt.close(fig)
print("Part 5: pivoting_error.png saved; unpivoted error at eps=1e-20 is", errs[-1, 0])

# %% [markdown]
# ## Part 6: Operation counts (§3.7) and factor-once-solve-many timing
#
# Exact flop count of Algorithm 3.5 without pivoting: sum_{m=1}^{n-1} (m + 2m^2)
# = n(n-1)/2 + n(n-1)(2n-1)/3  ~  (2/3) n^3.  One triangular solve ~ n^2.


# %%
def predicted_flops(n):
    return n * (n - 1) // 2 + n * (n - 1) * (2 * n - 1) // 3


ns = np.array([3, 5, 10, 20, 50, 100, 200])
flops = []
for n in ns:
    A = rng.standard_normal((n, n)) + n * np.eye(n)
    *_, f = lu_decompose(A, pivot=False, count_flops=True)
    assert f == predicted_flops(n), (n, f, predicted_flops(n))
    flops.append(f)
flops = np.array(flops, float)
assert predicted_flops(3) == 13                         # exercise B8(a): Example 4.1 by hand
fig, ax = plt.subplots(figsize=(6.4, 4))
ax.semilogx(ns, flops / ns.astype(float) ** 3, "o-", lw=2, ms=6, color="C0", label="counted flops / $n^3$")
ax.axhline(2 / 3, color="C1", lw=2, ls="--", label="2/3 (leading term)")
ax.set_xlabel("$n$"); ax.set_ylabel("flops / $n^3$"); ax.set_title("LU cost approaches $\\frac{2}{3}n^3$ from below")
ax.grid(alpha=.3); ax.legend(); fig.tight_layout(); fig.savefig(FIG / "flop_count.png", dpi=120); plt.close(fig)
print("Part 6: flop counter matches the exact formula; flops/n^3 at n=200 is", round(flops[-1] / 200**3, 4))

# Timing: LU once + k substitutions vs. k solves from scratch (our own code for both).
k = 10
sizes = [50, 100, 200, 300, 400]
t_once, t_scratch = [], []
for n in sizes:
    A = rng.standard_normal((n, n)); B = rng.standard_normal((n, k))
    t0 = time.perf_counter()
    P, L, U = lu_decompose(A)
    X1 = np.column_stack([back_sub(U, forward_sub(L, P @ B[:, j])) for j in range(k)])
    t_once.append(time.perf_counter() - t0)
    t0 = time.perf_counter()
    X2 = np.column_stack([solve(A, B[:, j]) for j in range(k)])
    t_scratch.append(time.perf_counter() - t0)
    assert np.allclose(X1, X2) and np.allclose(A @ X1, B, atol=1e-8)
fig, ax = plt.subplots(figsize=(6.4, 4))
ax.loglog(sizes, t_scratch, "o-", lw=2, ms=6, color="C1", label=f"{k} solves from scratch")
ax.loglog(sizes, t_once, "s-", lw=2, ms=6, color="C0", label=f"LU once + {k} substitutions")
ax.set_xlabel("$n$"); ax.set_ylabel("seconds"); ax.set_title("Factor once, solve many")
ax.grid(alpha=.3, which="both"); ax.legend(); fig.tight_layout(); fig.savefig(FIG / "timing_lu_vs_solve.png", dpi=120); plt.close(fig)
print("Timing (n, once, scratch):", [(n, round(a, 3), round(b, 3)) for n, a, b in zip(sizes, t_once, t_scratch)])

# %% [markdown]
# ## Part 7: Exercise numbers (so solutions.md and this script agree)

# %%
# A1
P, L, U = lu_decompose([[3, 1], [6, 5]], pivot=False)
assert np.allclose(L, [[1, 0], [2, 1]]) and np.allclose(U, [[3, 1], [0, 3]]) and np.allclose(solve([[3, 1], [6, 5]], [4, 11]), [1, 1])
# A2
A = np.arange(1, 10).reshape(3, 3).astype(float)
assert np.allclose(E(3, 2, 0, 4) @ A, [[1, 2, 3], [4, 5, 6], [3, 0, -3]])
assert np.allclose(A @ E(3, 2, 0, 4), [[-11, 2, 3], [-20, 5, 6], [-29, 8, 9]])
# A3
R, piv = rref([[1, 2, 1, 3], [2, 4, 0, 2], [1, 2, -1, -1]])
assert np.allclose(R, [[1, 2, 0, 1], [0, 0, 1, 2], [0, 0, 0, 0]]) and piv == [0, 2]
# A4
x = solve(A41, [2, 7, -5], pivot=False); assert np.allclose(x, [7 / 4, 0, -3 / 2])
assert np.allclose(forward_sub(L_41, [2, 7, -5]), [2, 3, 6])
# A5
assert np.allclose(gauss_jordan_inverse([[2, 3], [1, 4]]), np.array([[4, -3], [-1, 2]]) / 5)
# A6
Pc = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], float)
assert np.allclose(Pc @ Pc, Pc.T) and np.allclose(Pc @ Pc @ Pc, np.eye(3))
# B2
A_b2 = np.array([[1, 2, 0], [3, 4, 1], [2, 5, -1]], float)
P, L, U = lu_decompose(A_b2, pivot=False)
assert np.allclose(L, [[1, 0, 0], [3, 1, 0], [2, -1 / 2, 1]]) and np.allclose(U, [[1, 2, 0], [0, -2, 1], [0, 0, -1 / 2]])
assert np.allclose(solve(A_b2, [-1, 1, -5], pivot=False), [1, -1, 2])
assert np.allclose(forward_sub(L, [-1, 1, -5]), [-1, 4, -1])
Efwd_b2 = E(3, 2, 1, -1 / 2) @ E(3, 2, 0, 2) @ E(3, 1, 0, 3)
assert np.allclose(Efwd_b2, [[1, 0, 0], [-3, 1, 0], [-7 / 2, 1 / 2, 1]])
# B3(b)
A_b3 = np.array([[0, 2, 1], [1, 1, 0], [2, 3, 1]], float)
P, L, U = lu_decompose(A_b3)
assert np.allclose(P, [[0, 0, 1], [1, 0, 0], [0, 1, 0]])
assert np.allclose(L, [[1, 0, 0], [0, 1, 0], [1 / 2, -1 / 4, 1]]) and np.allclose(U, [[2, 3, 1], [0, 2, 1], [0, 0, -1 / 4]])
# B4
A_b4 = np.array([[1, 1, 1], [1, 2, 3], [2, 3, 4]], float)
out = classify(A_b4, [1, 2, 3]); assert out["kind"] == "infinite"
assert np.allclose(out["x_p"], [0, 1, 0]) and np.allclose(out["specials"][0], [1, -2, 1]) and rank(A_b4) == 2
assert classify(A_b4, [1, 2, 4])["kind"] == "none"            # b3 != b1 + b2
# B6
A_b6 = np.array([[1, 0, 2], [2, 1, 3], [0, 1, 1]], float)
assert np.allclose(gauss_jordan_inverse(A_b6), [[-1, 1, -1], [-1, 1 / 2, 1 / 2], [1, -1 / 2, 1 / 2]])
P, L, U = lu_decompose(A_b6, pivot=False)
assert np.allclose(L, [[1, 0, 0], [2, 1, 0], [0, 1, 1]]) and np.allclose(np.diag(U), [1, 1, 2])
# B7(b): x^T A x as weighted squares for Example 4.3
x = rng.standard_normal(3)
q_direct = x @ A43 @ x
q_squares = 2 * (x[0] + x[1] / 2) ** 2 + 1.5 * (x[1] + 2 * x[2] / 3) ** 2 + (4 / 3) * x[2] ** 2
assert np.isclose(q_direct, q_squares) and np.isclose(q_direct, 2 * x[0]**2 + 2 * x[0] * x[1] + 2 * x[1]**2 + 2 * x[1] * x[2] + 2 * x[2]**2)
# B8(b): k at which from-scratch is 100x the factor-once cost for n = 500
n = 500; k_star = (100 * 2 * n**3 / 3) / (2 * n**3 / 3 - 200 * n**2)
assert np.isclose(k_star, 250)
# C1: Wilkinson's matrix, growth 2^(n-1) with no swaps
def wilkinson(n):
    W = np.eye(n) - np.tril(np.ones((n, n)), -1); W[:, -1] = 1.0; return W
P, L, U = lu_decompose(wilkinson(4))
assert np.allclose(P, np.eye(4)) and np.allclose(U[:, -1], [1, 2, 4, 8])
for n in (10, 30):
    P, L, U = lu_decompose(wilkinson(n)); assert np.allclose(P, np.eye(n)) and np.isclose(np.abs(U).max(), 2.0 ** (n - 1))
W60 = wilkinson(60); b60 = W60 @ np.ones(60)
x60 = sla.lu_solve(sla.lu_factor(W60), b60)
print(f"C1: Wilkinson n=60: max|U| = {np.abs(sla.lu(W60)[2]).max():.3g}, max|x-1| = {np.abs(x60-1).max():.3g}, "
      f"rel. residual = {np.linalg.norm(W60@x60-b60)/np.linalg.norm(b60):.2g}")
# C2: Schur complements of the circuit matrix
S1 = G[1:, 1:] - np.outer(G[1:, 0], G[0, 1:]) / G[0, 0]
assert np.isclose(S1[0, 0], 5 / 6) and np.allclose(S1, [[5 / 6, -1 / 2], [-1 / 2, 3 / 2]])
S2 = G[2:, 2:] - G[2:, :2] @ np.linalg.solve(G[:2, :2], G[:2, 2:])
assert np.isclose(S2[0, 0], 6 / 5)
# C4: Sherman–Morrison
A = rng.standard_normal((6, 6)); u = rng.standard_normal(6); v = rng.standard_normal(6)
Ai = np.linalg.inv(A); sigma = 1 + v @ Ai @ u
assert np.allclose(np.linalg.inv(A + np.outer(u, v)), Ai - np.outer(Ai @ u, v @ Ai) / sigma)
print("Part 7: all exercise numbers verified.")

# %% [markdown]
# ## Your turn
#
# 1. **LDL^T (exercise B9).** Implement `ldlt(A)` for symmetric A with no row
#    swaps: `L, d` with A = L diag(d) L^T, computing only the lower triangle
#    (Theorem 3.22). A reference implementation follows so you can compare.
# 2. **Rank by pivot count.** Build a 6x6 matrix of exact rank 2 as B @ C with B
#    6x2. Call `rank(A, tol)` with tol = 1e-10 and tol = 0. Explain the difference.
# 3. **A random singular matrix.** Feed the same A to `solve(A, b)` with
#    tol=0 and with tol=1e-10. Look at the size of x and of the residual A x - b
#    in the first case. What is the code "solving"?
# 4. Add a `pivot="complete"` option to `lu_decompose` that also swaps columns
#    (P A Q = L U). What has to happen to x afterwards?


# %%
def ldlt(A):
    """A = L diag(d) L^T for symmetric A without row swaps (Theorem 3.22, eq. 3.10).
    Only the lower triangle of A is read; about n^3/3 flops."""
    A = np.asarray(A, float); n = A.shape[0]
    assert np.allclose(A, A.T), "ldlt needs a symmetric matrix"
    L = np.eye(n); d = np.zeros(n)
    for j in range(n):
        d[j] = A[j, j] - np.sum(L[j, :j] ** 2 * d[:j])          # pivot j, (3.10) entry (j, j)
        if d[j] == 0:
            raise np.linalg.LinAlgError("zero pivot: LDL^T without pivoting fails")
        for i in range(j + 1, n):
            L[i, j] = (A[i, j] - np.sum(L[i, :j] * L[j, :j] * d[:j])) / d[j]   # (3.10) entry (i, j)
    return L, d


for n in (3, 6, 12):
    B = rng.standard_normal((n, n)); S = B @ B.T + np.eye(n)
    L, d = ldlt(S)
    assert np.allclose(L @ np.diag(d) @ L.T, S, atol=1e-10)
    P, L2, U2 = lu_decompose(S, pivot=False)
    assert np.allclose(L, L2) and np.allclose(d, np.diag(U2)), "d must be the pivots of elimination"
    assert np.all(d > 0)        # positive definite: all pivots positive (LA-8 preview, exercise B7)
L, d = ldlt(A43)
assert np.allclose(d, [2, 3 / 2, 4 / 3]) and np.isclose(L[1, 0], 1 / 2) and np.isclose(L[2, 1], 2 / 3)
print("Your turn 1: ldlt reproduces Example 4.3 (d =", d, ") and random SPD matrices.")

# Rank by pivot count needs a tolerance: a matrix of exact rank 2 has a third "pivot" of size ~1e-16.
B = rng.standard_normal((6, 2)); C = rng.standard_normal((2, 6)); A_sing = B @ C
print("Your turn 2: exact rank 2; pivot count with tol=1e-10:", rank(A_sing, 1e-10), "; with tol=0:", rank(A_sing, 0.0))
assert rank(A_sing, 1e-10) == 2 and np.linalg.matrix_rank(A_sing) == 2
P, L, U = lu_decompose(A_sing, pivot=True, tol=0.0)
print("             pivots of the singular matrix:", np.array2string(np.diag(U), precision=2))
try:
    solve(A_sing, np.ones(6), tol=1e-10); raise AssertionError
except np.linalg.LinAlgError as err:
    print("Your turn 3: with tol=1e-10 solve raises:", err)
x_garbage = solve(A_sing, np.ones(6), tol=0.0)
print(f"             with tol=0 solve returns |x| = {np.linalg.norm(x_garbage):.2e}, residual |Ax-b| = {np.linalg.norm(A_sing@x_garbage-np.ones(6)):.2e}")

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
