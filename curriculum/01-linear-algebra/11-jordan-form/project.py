# %% [markdown]
# # LA-11 Project: Jordan Form, Nilpotent Chains, Cayley–Hamilton, Functions of Matrices
#
# This script builds the machinery of LA-11 in exact arithmetic (SymPy) and
# then tests it against floating-point reference implementations (SciPy).
# Parts: (1) kernel dimensions and the block-counting formula (3.4),
# (2) the chain construction of Theorem 3.4, built by hand,
# (3) the ODE y'' - 2y' + y = 0 and where t e^t comes from (Example 4.1),
# (4) the 3x3 and 4x4 worked examples (4.2, 4.3), (5) Cayley–Hamilton
# (Example 4.4, random matrices, and the Schur mechanism of the proof),
# (6) the minimal polynomial two ways (Example 4.5) and why p, m, rank are not
# enough to pin down a Jordan form, (7) f(A) via the Jordan formula (3.7),
# exp(At) against scipy.linalg.expm, and the critically damped oscillator,
# (8) the discontinuity of the Jordan form in the entries.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from scipy.linalg import expm, schur
from scipy.integrate import solve_ivp
from pathlib import Path

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
x = sp.symbols("x")
t = sp.symbols("t", real=True)


def M(rows):
    """Exact SymPy matrix from a nested list of integers/rationals."""
    return sp.Matrix(rows)


def to_np(A):
    return np.array(A.tolist(), dtype=float)


# %% [markdown]
# ## Part 1: Kernel dimensions and the block-counting formula
#
# For an eigenvalue lam, write d_k = dim ker (A - lam I)^k (exact rank, so
# exact dimension). Theorem 3.9: the number of Jordan blocks for lam of size
# >= k is d_k - d_{k-1}, eq. (3.4). This is the "easy part" of the module and
# the one to learn cold.


# %%
def kernel_dims(A, lam):
    """d_k = dim ker (A - lam I)^k for k = 0, 1, ..., until the sequence stops
    growing (Lemma 3.2: once equal, equal forever). Exact arithmetic."""
    n = A.rows
    N = A - lam * sp.eye(n)
    dims = [0]
    Nk = sp.eye(n)
    while True:
        Nk = Nk * N
        dims.append(n - Nk.rank())
        if dims[-1] == dims[-2]:
            return dims[:-1]


def jordan_block_sizes(A, lam):
    """Sizes of the Jordan blocks of A for eigenvalue lam, from eq. (3.4):
    #blocks of size >= k equals d_k - d_{k-1}. Returned largest first."""
    d = kernel_dims(A, lam)
    at_least = [d[k] - d[k - 1] for k in range(1, len(d))]  # at_least[k-1] = #blocks of size >= k
    sizes = []
    for k in range(len(at_least), 0, -1):
        exactly = at_least[k - 1] - (at_least[k] if k < len(at_least) else 0)
        sizes += [k] * exactly
    assert sum(sizes) == d[-1], "block sizes must add up to dim G_lam"
    return sizes


def jordan_structure(A):
    """{eigenvalue: block sizes} computed only from ranks (no eigenvectors)."""
    return {lam: jordan_block_sizes(A, lam) for lam in A.eigenvals()}


def blocks_of_J(J):
    """Read the (eigenvalue, size) blocks off a Jordan matrix J, in order."""
    n = J.rows
    blocks, i = [], 0
    while i < n:
        s = 1
        while i + s < n and J[i + s - 1, i + s] != 0:
            s += 1
        blocks.append((J[i, i], s))
        i += s
    return blocks


def sympy_structure(A):
    """Same dictionary as jordan_structure, but read off sympy's Matrix.jordan_form()."""
    _, J = A.jordan_form()
    out = {}
    for lam, s in blocks_of_J(J):
        out.setdefault(lam, []).append(s)
    return {lam: sorted(v, reverse=True) for lam, v in out.items()}


# Single Jordan blocks: d_k = min(s, k), so exactly one block of size s.
for s in (1, 2, 3, 5):
    Js = sp.jordan_cell(sp.Rational(7, 2), s)
    assert kernel_dims(Js, sp.Rational(7, 2)) == [min(s, k) for k in range(s + 1)]
    assert jordan_block_sizes(Js, sp.Rational(7, 2)) == [s]

# A direct sum of blocks: the formula must recover the sizes we put in.
Jmix = sp.diag(sp.jordan_cell(2, 3), sp.jordan_cell(2, 1), sp.jordan_cell(2, 2), sp.jordan_cell(-1, 2))
assert jordan_block_sizes(Jmix, 2) == [3, 2, 1]
assert jordan_block_sizes(Jmix, -1) == [2]
assert kernel_dims(Jmix, 2) == [0, 3, 5, 6]        # Exercise A3's numbers
print("Part 1: block-counting formula (3.4) recovers the sizes of every test matrix.")

# %% [markdown]
# ## Part 2: The chain construction (Theorem 3.4), built by hand
#
# **Build it yourself first.** For a nilpotent N, work from the top level
# m (where ker N^m is everything) down to 1. At level k you already hold the
# vectors N*(level k+1 vectors); add new vectors from ker N^k until, together
# with ker N^{k-1}, they span ker N^k. Each new vector starts a chain of length
# k. Then push everything down one level by applying N. The bottoms of the
# chains land in ker N; Jordan basis = chains written bottom-first.
#
# `dim` lets the same routine run on N = A - lam I restricted to the
# generalized eigenspace G_lam (dimension = algebraic multiplicity), where
# N is nilpotent even though it is not nilpotent on all of C^n.


# %%
def nilpotent_chains(N, dim=None):
    """Jordan chains of the nilpotent operator N (on ker N^m, of dimension dim).
    Returns a list of chains, each [v, N v, N^2 v, ..., N^{s-1} v] (top first),
    with N^s v = 0. Proof/algorithm: Theorem 3.4 and §3.6 of the notes."""
    n = N.rows
    dim = n if dim is None else dim
    kers = [sp.zeros(n, 0)]           # kers[k] = basis matrix of ker N^k
    Nk, k = sp.eye(n), 0
    while kers[-1].cols < dim:
        Nk, k = Nk * N, k + 1
        assert k <= n, "N is not nilpotent on a space of this dimension"
        basis = Nk.nullspace()
        kers.append(sp.Matrix.hstack(*basis) if basis else sp.zeros(n, 0))
    m = k                              # index of nilpotency
    chains = []
    for k in range(m, 0, -1):
        current = [ch[-1] for ch in chains]              # already at level k
        base = sp.Matrix.hstack(kers[k - 1], *current)
        assert base.rank() == base.cols, "level-k vectors must be independent mod ker N^{k-1}"
        new = []
        for j in range(kers[k].cols):
            v = kers[k][:, j]
            cand = sp.Matrix.hstack(base, *new, v)
            if cand.rank() == cand.cols:
                new.append(v)
        # eq. (3.4): number of chains of length >= k is d_k - d_{k-1}
        assert len(current) + len(new) == kers[k].cols - kers[k - 1].cols
        chains += [[v] for v in new]
        if k > 1:
            for ch in chains:
                ch.append(N * ch[-1])
    for ch in chains:
        assert N * ch[-1] == sp.zeros(n, 1)
    return chains


def jordan_basis(A):
    """P and J with P^-1 A P = J, assembled from chains in each generalized
    eigenspace (Theorem 3.8). Also returns the (eigenvalue, size) block list."""
    n = A.rows
    cols, blocks = [], []
    eig = sorted(A.eigenvals().items(), key=lambda kv: (float(sp.re(kv[0])), float(sp.im(kv[0]))))
    for lam, mult in eig:
        for ch in nilpotent_chains(A - lam * sp.eye(n), dim=mult):
            cols += list(reversed(ch))                   # bottom first: superdiagonal ones
            blocks.append((lam, len(ch)))
    P = sp.Matrix.hstack(*cols)
    J = sp.diag(*[sp.jordan_cell(lam, s) for lam, s in blocks])
    assert P.rank() == n and P.inv() * A * P == J
    return P, J, blocks


def check_against_sympy(A):
    """Our rank-only structure, our chain construction, and sympy must agree."""
    ours = jordan_structure(A)
    theirs = sympy_structure(A)
    assert ours == theirs, (ours, theirs)
    P, J, blocks = jordan_basis(A)
    built = {}
    for lam, s in blocks:
        built.setdefault(lam, []).append(s)
    assert {k: sorted(v, reverse=True) for k, v in built.items()} == ours
    return P, J


# A nilpotent test: shift of size 3 plus shift of size 2, scrambled by a change of basis.
S = M([[1, 2, 0, 0, 1], [0, 1, 1, 0, 0], [1, 0, 1, 1, 0], [0, 0, 1, 1, 1], [1, 1, 0, 0, 1]])
assert S.det() != 0
Ntest = S * sp.diag(sp.jordan_cell(0, 3), sp.jordan_cell(0, 2)) * S.inv()
chains = nilpotent_chains(Ntest)
assert sorted(len(c) for c in chains) == [2, 3]
check_against_sympy(Ntest)
print("Part 2: chain construction recovers blocks [3, 2] for a scrambled nilpotent 5x5.")

# %% [markdown]
# ## Part 3: Example 4.1, the ODE y'' - 2y' + y = 0
#
# A = [[0, 1], [-1, 2]] has p_A = (x - 1)^2 and one eigenvector. Jordan form
# J_2(1); e^{At} = e^t [[1 - t, t], [-t, 1 + t]]; y(t) = e^t (y0 + t (y0' - y0)).


# %%
A_ode = M([[0, 1], [-1, 2]])
assert sp.factor(A_ode.charpoly(x).as_expr()) == (x - 1) ** 2
assert kernel_dims(A_ode, 1) == [0, 1, 2]
assert jordan_block_sizes(A_ode, 1) == [2]
P_notes = M([[1, 0], [1, 1]])                       # chain v2 = e2 -> v1 = (1, 1)
assert (A_ode - sp.eye(2)) * M([0, 1]) == M([1, 1])
assert P_notes.inv() * A_ode * P_notes == sp.jordan_cell(1, 2)
check_against_sympy(A_ode)

E_formula = sp.exp(t) * M([[1 - t, t], [-t, 1 + t]])
assert sp.simplify((A_ode * t).exp() - E_formula) == sp.zeros(2, 2)
assert sp.simplify(sp.diff(E_formula, t) - A_ode * E_formula) == sp.zeros(2, 2)   # X' = A X
assert E_formula.subs(t, 0) == sp.eye(2)

# Numerical: y(0) = 1, y'(0) = 0 gives y = e^t (1 - t). Compare the Jordan
# solution with scipy.linalg.expm and with a black-box ODE integrator.
A_np = to_np(A_ode)
ts = np.linspace(0, 3, 61)
y_exact = np.exp(ts) * (1 - ts)
y_expm = np.array([(expm(A_np * tv) @ [1.0, 0.0])[0] for tv in ts])
sol = solve_ivp(lambda _, z: A_np @ z, (0, 3), [1.0, 0.0], t_eval=ts, rtol=1e-10, atol=1e-12)
assert np.allclose(y_expm, y_exact, rtol=1e-10, atol=1e-10)
assert np.allclose(sol.y[0], y_exact, rtol=1e-6, atol=1e-7)

fig, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(ts, np.exp(ts), "--", label=r"$e^{t}$ (the eigenvector direction)")
ax.plot(ts, ts * np.exp(ts), "--", label=r"$t\,e^{t}$ (from the shift $N$)")
ax.plot(ts, y_exact, "k", lw=2, label=r"$y(t) = e^{t}(1-t)$,  $y(0)=1, y'(0)=0$")
ax.plot(ts[::6], sol.y[0][::6], "o", ms=4, label="solve_ivp")
ax.set_xlabel("t"); ax.set_ylabel("y"); ax.set_ylim(-45, 25); ax.grid(alpha=.3); ax.legend(fontsize=8)
ax.set_title("Example 4.1: repeated eigenvalue, one eigenvector, a $t\\,e^{t}$ term")
fig.tight_layout(); fig.savefig(FIG / "ode_repeated_root.png", dpi=120); plt.close(fig)
print("Part 3: e^{At} = e^t[[1-t, t],[-t, 1+t]] verified symbolically, against expm, and against solve_ivp.")

# %% [markdown]
# ## Part 4: Examples 4.2 and 4.3
#
# 4.2: A = [[1,1,0],[-1,2,1],[-1,0,3]], eigenvalue 2 with a_2 = 3, g_2 = 1, so one
# 3-block. Chain from e3: (1,1,1) <- (0,1,1) <- (0,0,1).
# 4.3: nilpotent 4x4 with d = (2, 3, 4): one 3-block and one 1-block.
# Chain from e3: (1,0,1,0) <- (0,1,-1,1) <- (0,0,1,0); the 1-block is (1,0,0,1).


# %%
A3 = M([[1, 1, 0], [-1, 2, 1], [-1, 0, 3]])
assert sp.factor(A3.charpoly(x).as_expr()) == (x - 2) ** 3
N3 = A3 - 2 * sp.eye(3)
assert N3 == M([[-1, 1, 0], [-1, 0, 1], [-1, 0, 1]])
assert N3 ** 2 == M([[0, -1, 1], [0, -1, 1], [0, -1, 1]]) and N3 ** 3 == sp.zeros(3, 3)
assert kernel_dims(A3, 2) == [0, 1, 2, 3] and jordan_block_sizes(A3, 2) == [3]
e3 = M([0, 0, 1])
v2, v1 = N3 * e3, N3 ** 2 * e3
assert v2 == M([0, 1, 1]) and v1 == M([1, 1, 1]) and N3 * v1 == sp.zeros(3, 1)
P3 = sp.Matrix.hstack(v1, v2, e3)
assert P3 == M([[1, 0, 0], [1, 1, 0], [1, 1, 1]])
assert P3.inv() * A3 * P3 == sp.jordan_cell(2, 3)
check_against_sympy(A3)
print("Example 4.2: one 3-block for lambda = 2; chain (1,1,1) <- (0,1,1) <- (0,0,1).")

A4 = M([[0, 1, 0, 0], [-1, 0, 1, 1], [1, 1, -1, -1], [-1, 0, 1, 1]])
assert A4 ** 2 == M([[-1, 0, 1, 1], [0, 0, 0, 0], [-1, 0, 1, 1], [0, 0, 0, 0]])
assert A4 ** 3 == sp.zeros(4, 4)
d4 = kernel_dims(A4, 0)
assert d4 == [0, 2, 3, 4], d4
assert [d4[k] - d4[k - 1] for k in (1, 2, 3)] == [2, 1, 1]   # blocks of size >= 1, 2, 3
assert jordan_block_sizes(A4, 0) == [3, 1]
e3 = M([0, 0, 1, 0])
u2, u1 = A4 * e3, A4 ** 2 * e3
assert u2 == M([0, 1, -1, 1]) and u1 == M([1, 0, 1, 0]) and A4 * u1 == sp.zeros(4, 1)
w = M([1, 0, 0, 1])
assert A4 * w == sp.zeros(4, 1)
P4 = sp.Matrix.hstack(u1, u2, e3, w)
assert P4.det() == 1
assert P4.inv() * A4 * P4 == sp.diag(sp.jordan_cell(0, 3), sp.jordan_cell(0, 1))
chains4 = nilpotent_chains(A4)
assert sorted(len(c) for c in chains4) == [1, 3]
check_against_sympy(A4)
print("Example 4.3: d = (2, 3, 4) => blocks [3, 1]; P^-1 A P = J_3(0) + J_1(0) with the notes' P.")

# Young-diagram picture of the kernel dimensions (the "staircase" of §2).
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
for ax, (name, Amat) in zip(axes, [("Example 4.3 (4x4)", A4), ("Example 4.2 (3x3), N = A - 2I", N3)]):
    d = kernel_dims(Amat, 0)
    ax.step(range(len(d)), d, where="post", lw=2)
    ax.plot(range(len(d)), d, "o")
    for k in range(1, len(d)):
        ax.annotate(f"+{d[k]-d[k-1]}", (k - .5, (d[k] + d[k-1]) / 2), ha="center", fontsize=9)
    ax.set_xticks(range(len(d))); ax.set_xlabel("k"); ax.set_ylabel(r"$d_k = \dim\ker N^k$")
    ax.set_title(name); ax.grid(alpha=.3)
fig.suptitle("Jumps $d_k - d_{k-1}$ = number of blocks of size $\\geq k$")
fig.tight_layout(); fig.savefig(FIG / "kernel_dimensions.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 5: Cayley–Hamilton
#
# Example 4.4: A = [[2,1,0],[1,1,0],[0,1,1]], p_A = x^3 - 4x^2 + 4x - 1,
# A^-1 = A^2 - 4A + 4I, A^5 = 33A^2 - 44A + 12I = [[89,55,0],[55,34,0],[33,22,1]].
# Then p_A(A) = 0 numerically for random matrices, and the mechanism of the
# Schur proof (Theorem 3.15): (T - l_1 I)...(T - l_k I) kills span(e_1..e_k).


# %%
A5 = M([[2, 1, 0], [1, 1, 0], [0, 1, 1]])
p5 = A5.charpoly(x).as_expr()
assert sp.expand(p5) == x ** 3 - 4 * x ** 2 + 4 * x - 1
assert A5 ** 3 - 4 * A5 ** 2 + 4 * A5 - sp.eye(3) == sp.zeros(3, 3)       # p_A(A) = 0
assert A5 ** 2 == M([[5, 3, 0], [3, 2, 0], [1, 2, 1]])
assert A5.inv() == A5 ** 2 - 4 * A5 + 4 * sp.eye(3) == M([[1, -1, 0], [-1, 2, 0], [1, -2, 1]])
_, rem = sp.div(sp.Poly(x ** 5, x), sp.Poly(p5, x))
assert rem.as_expr() == 33 * x ** 2 - 44 * x + 12
assert A5 ** 5 == 33 * A5 ** 2 - 44 * A5 + 12 * sp.eye(3) == M([[89, 55, 0], [55, 34, 0], [33, 22, 1]])
print("Example 4.4: A^-1 = A^2 - 4A + 4I and A^5 = 33A^2 - 44A + 12I verified exactly.")


def charpoly_coeffs(A):
    """Coefficients of p_A(x) = det(xI - A), highest degree first (via eigenvalues)."""
    return np.poly(A)


def polyval_matrix(coeffs, A):
    """p(A) by Horner's rule; coeffs highest degree first."""
    n = A.shape[0]
    out = np.zeros_like(A, dtype=complex)
    for c in coeffs:
        out = out @ A + c * np.eye(n)
    return out


for n in (2, 3, 5, 8):
    for _ in range(20):
        A = rng.standard_normal((n, n))
        c = charpoly_coeffs(A)
        res = polyval_matrix(c, A)
        scale = sum(abs(ck) * np.linalg.norm(A, 2) ** (n - k) for k, ck in enumerate(c))
        assert np.linalg.norm(res) / scale < 1e-12, "Cayley–Hamilton residual too large"

# The Schur mechanism: with T upper triangular (eigenvalues on the diagonal),
# M_k = prod_{j<=k} (T - l_j I) annihilates V_k = span(e_1, ..., e_k).
A = rng.standard_normal((6, 6))
T, Z = schur(A, output="complex")
lam = np.diag(T)
Mk = np.eye(6, dtype=complex)
for k in range(6):
    Mk = Mk @ (T - lam[k] * np.eye(6))
    assert np.linalg.norm(Mk[:, : k + 1]) < 1e-10 * np.linalg.norm(T) ** (k + 1), k
assert np.linalg.norm(Mk) < 1e-9 * np.linalg.norm(T) ** 6            # M_n = p_A(T) = 0
print("Part 5: p_A(A) = 0 on random matrices; the successive products in the Schur proof kill V_1, V_2, ...")

# %% [markdown]
# ## Part 6: The minimal polynomial, two ways
#
# Way 1 (Theorem 3.14): m_A = prod (x - lam)^{largest block for lam}.
# Way 2 (Definition 3.12): smallest d with I, A, ..., A^d linearly dependent;
# the (monic) dependence is m_A. They must agree. Example 4.5:
# J_2(1) + J_1(1) + J_1(3) has p = (x-1)^3 (x-3) and m = (x-1)^2 (x-3).


# %%
def minimal_polynomial_from_blocks(A):
    """m_A(x) = prod_lam (x - lam)^{s_max(lam)}, eq. (3.5)."""
    return sp.expand(sp.prod([(x - lam) ** max(sizes) for lam, sizes in jordan_structure(A).items()]))


def minimal_polynomial_by_dependence(A):
    """The monic polynomial of least degree with m(A) = 0, found as the first
    linear dependence among I, A, A^2, ... (vectorized to n^2-vectors)."""
    n = A.rows
    powers = [sp.eye(n)]
    for d in range(1, n + 1):
        powers.append(powers[-1] * A)
        K = sp.Matrix.hstack(*[Pk.reshape(n * n, 1) for Pk in powers])
        null = K.nullspace()
        if null:
            c = null[0]
            assert len(null) == 1 and c[-1] != 0
            c = c / c[-1]
            return sp.expand(sum(c[k] * x ** k for k in range(d + 1)))
    raise AssertionError("Cayley–Hamilton guarantees a dependence by degree n")


def poly_of_matrix(poly, A):
    """poly(A) exactly, for a sympy polynomial expression in x."""
    n = A.rows
    out = sp.zeros(n, n)
    for k, c in enumerate(reversed(sp.Poly(poly, x).all_coeffs())):
        out += c * A ** k
    return out


J45 = sp.diag(sp.jordan_cell(1, 2), sp.jordan_cell(1, 1), sp.jordan_cell(3, 1))
Pr = M([[1, 2, 0, 1], [0, 1, 1, 0], [1, 0, 1, 1], [0, 1, 0, 1]])
assert Pr.det() != 0
A45 = Pr * J45 * Pr.inv()                               # same Jordan form, dense entries
m45 = minimal_polynomial_from_blocks(A45)
assert m45 == sp.expand((x - 1) ** 2 * (x - 3))
assert minimal_polynomial_by_dependence(A45) == m45
assert sp.factor(A45.charpoly(x).as_expr()) == (x - 3) * (x - 1) ** 3
assert poly_of_matrix(m45, A45) == sp.zeros(4, 4)
assert poly_of_matrix((x - 1) * (x - 3), A45) != sp.zeros(4, 4)      # one factor too few
assert jordan_structure(A45) == {1: [2, 1], 3: [1]}
# m_A divides p_A (Theorem 3.13 + Cayley–Hamilton):
assert sp.rem(sp.Poly(A45.charpoly(x).as_expr(), x), sp.Poly(m45, x)).is_zero
# Numerically, in floating point, m(A) and p(A) are zero to rounding:
A45f = to_np(A45)
m_res = np.linalg.norm(polyval_matrix([1, -5, 7, -3], A45f))          # (x-1)^2 (x-3) = x^3 - 5x^2 + 7x - 3
assert sp.expand((x - 1) ** 2 * (x - 3)) == x ** 3 - 5 * x ** 2 + 7 * x - 3
assert m_res < 1e-10
print(f"Example 4.5: p = (x-3)(x-1)^3, m = (x-1)^2(x-3) by both methods; ||m(A)|| in floating point = {m_res:.1e}")

# Both methods agree on every matrix so far, and on a diagonalizable one.
for Amat in (A_ode, A3, A4, A5, Jmix):
    assert minimal_polynomial_from_blocks(Amat) == minimal_polynomial_by_dependence(Amat)
Adiag = M([[2, 0, 0], [0, 2, 0], [0, 0, 5]])
assert minimal_polynomial_by_dependence(Adiag) == sp.expand((x - 2) * (x - 5))   # distinct roots <=> diagonalizable

# p, m and rank do NOT determine the Jordan form; the full sequence d_k does.
J7a = sp.diag(sp.jordan_cell(0, 3), sp.jordan_cell(0, 3), sp.jordan_cell(0, 1))
J7b = sp.diag(sp.jordan_cell(0, 3), sp.jordan_cell(0, 2), sp.jordan_cell(0, 2))
assert J7a.charpoly(x) == J7b.charpoly(x) == sp.Poly(x ** 7, x)
assert minimal_polynomial_by_dependence(J7a) == minimal_polynomial_by_dependence(J7b) == x ** 3
assert J7a.rank() == J7b.rank() == 4
assert kernel_dims(J7a, 0) == [0, 3, 5, 7] and kernel_dims(J7b, 0) == [0, 3, 6, 7]
assert jordan_structure(J7a) != jordan_structure(J7b)
print("Part 6: two 7x7 matrices with equal p, m and rank but d_2 = 5 vs 6: not similar.")

# Symmetric matrices never have blocks (Proposition 3.11): the ring of three masses.
K_ring = M([[2, -1, -1], [-1, 2, -1], [-1, -1, 2]])
assert K_ring.eigenvals() == {0: 1, 3: 2}
assert jordan_structure(K_ring) == {0: [1], 3: [1, 1]}
assert kernel_dims(K_ring, 3) == [0, 2]                 # ker (K - 3I)^2 = ker (K - 3I)

# %% [markdown]
# ## Part 7: Functions of matrices via the Jordan form
#
# Theorem 3.16, eq. (3.7): f(J_s(lam)) = sum_{k<s} f^{(k)}(lam)/k! N^k. Then
# f(A) = P f(J) P^-1 (3.10). We build exp(tA) this way with sympy's own
# jordan_form and compare with scipy.linalg.expm, then solve the critically
# damped oscillator x'' + 2w x' + w^2 x = 0 both ways.


# %%
def f_jordan_block(f, lam, s):
    """f(J_s(lam)) from the terminating Taylor series (3.7); f is an expression in x."""
    out = sp.zeros(s, s)
    for k in range(s):
        coeff = sp.diff(f, x, k).subs(x, lam) / sp.factorial(k)
        for i in range(s - k):
            out[i, i + k] = coeff
    return out


def f_of_matrix(A, f):
    """f(A) = P f(J) P^-1 using sympy's Jordan form (exact, for small matrices)."""
    P, J = A.jordan_form()
    FJ = sp.diag(*[f_jordan_block(f, lam, s) for lam, s in blocks_of_J(J)])
    return P * FJ * P.inv()


def expm_jordan(A, tval):
    """e^{t A} via the Jordan form; tval may be a number or a symbol."""
    return f_of_matrix(A, sp.exp(x * tval))


# The two block formulas in the notes, (3.8) and (3.9):
lam_s = sp.symbols("lambda")
assert sp.simplify(f_jordan_block(sp.exp(x * t), lam_s, 2) - sp.exp(lam_s * t) * M([[1, t], [0, 1]])) == sp.zeros(2, 2)
assert sp.simplify(f_jordan_block(sp.exp(x * t), lam_s, 3)
                   - sp.exp(lam_s * t) * M([[1, t, t ** 2 / 2], [0, 1, t], [0, 0, 1]])) == sp.zeros(3, 3)
# A polynomial f is reproduced exactly (the series really terminates):
assert f_of_matrix(A3, x ** 3 - 2 * x + 1) == A3 ** 3 - 2 * A3 + sp.eye(3)
# e^{tA} by Jordan == the closed form of Example 4.1 == scipy:
assert sp.simplify(expm_jordan(A_ode, t) - E_formula) == sp.zeros(2, 2)
for Amat in (A_ode, A3, A4, A45):
    for tv in (0.3, 1.0, 2.5):
        ours = to_np(expm_jordan(Amat, sp.Rational(tv).limit_denominator(1000)).evalf(30))
        ref = expm(to_np(Amat) * tv)
        assert np.allclose(ours, ref, rtol=1e-9, atol=1e-9), (Amat, tv)
# The cos of a 2x2 block, from (3.7): cos(lam) I - sin(lam) N  (reference for "Your turn" 3)
assert sp.simplify(f_jordan_block(sp.cos(x), lam_s, 2) - M([[sp.cos(lam_s), -sp.sin(lam_s)], [0, sp.cos(lam_s)]])) == sp.zeros(2, 2)
print("Part 7: f(J) formula gives (3.8), (3.9); e^{tA} via Jordan matches scipy.linalg.expm on four matrices.")

# Critically damped oscillator: x'' + 2w x' + w^2 x = 0, A = [[0,1],[-w^2,-2w]].
w = sp.symbols("omega", positive=True)
A_cd = M([[0, 1], [-w ** 2, -2 * w]])
assert sp.factor(A_cd.charpoly(x).as_expr()) == (x + w) ** 2
assert (A_cd + w * sp.eye(2)).rank() == 1                                   # one eigenvector: a 2-block
E_cd = sp.exp(-w * t) * M([[1 + w * t, t], [-w ** 2 * t, 1 - w * t]])      # Example 4.6
assert sp.simplify((A_cd * t).exp() - E_cd) == sp.zeros(2, 2)
assert sp.simplify(sp.diff(E_cd, t) - A_cd * E_cd) == sp.zeros(2, 2)

omega = 2.0
A_cd_np = np.array([[0, 1], [-omega ** 2, -2 * omega]])
ts = np.linspace(0, 5, 251)
fig, ax = plt.subplots(figsize=(7, 4.2))
for x0, v0, style in ((0.0, 1.0, "k"), (1.0, 0.0, "C0"), (1.0, -4.0, "C3")):
    x_formula = np.exp(-omega * ts) * ((1 + omega * ts) * x0 + ts * v0)
    x_expm = np.array([(expm(A_cd_np * tv) @ [x0, v0])[0] for tv in ts])
    sol = solve_ivp(lambda _, z: A_cd_np @ z, (0, 5), [x0, v0], t_eval=ts, rtol=1e-10, atol=1e-12)
    assert np.allclose(x_formula, x_expm, atol=1e-10) and np.allclose(x_formula, sol.y[0], atol=1e-7)
    ax.plot(ts, x_formula, style, lw=2, label=rf"$x(0)={x0:g},\ \dot x(0)={v0:g}$")
ax.plot(ts, ts * np.exp(-omega * ts), "--", color="gray", label=r"$t\,e^{-\omega t}$")
ax.plot(ts, np.exp(-omega * ts), ":", color="gray", label=r"$e^{-\omega t}$")
ax.axhline(0, color="k", lw=.5)
ax.set_xlabel("t"); ax.set_ylabel("x(t)"); ax.grid(alpha=.3); ax.legend(fontsize=8)
ax.set_title(r"Critically damped: $x(t) = e^{-\omega t}[(1+\omega t)x_0 + t\,\dot x_0]$, $\omega = 2$")
fig.tight_layout(); fig.savefig(FIG / "critically_damped.png", dpi=120); plt.close(fig)
# The pure-shift solution t e^{-wt} peaks at t = 1/w with value 1/(e w):
tpk = ts[np.argmax(ts * np.exp(-omega * ts))]
assert abs(tpk - 1 / omega) < 0.03 and np.isclose((ts * np.exp(-omega * ts)).max(), 1 / (np.e * omega), rtol=1e-3)
print("Part 7: critically damped oscillator solved via the Jordan block, expm and solve_ivp agree.")

# %% [markdown]
# ## Part 8: The Jordan form is discontinuous in the entries
#
# A_eps = [[1, 1], [eps, 1]] has eigenvalues 1 +- sqrt(eps) and eigenvectors
# (1, +-sqrt(eps)): diagonalizable for every eps != 0, a 2-block at eps = 0.
# The angle between the eigenvectors is 2 arctan(sqrt(eps)) -> 0 and the
# condition number of the eigenvector matrix is exactly 1/sqrt(eps) -> infinity.


# %%
eps_vals = np.logspace(-14, -1, 27)
angles, conds = [], []
for e in eps_vals:
    vals, V = np.linalg.eig(np.array([[1.0, 1.0], [e, 1.0]]))
    v1, v2 = V[:, 0], V[:, 1]
    angles.append(np.arccos(np.clip(abs(v1 @ v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), 0, 1)))
    conds.append(np.linalg.cond(V))
    assert np.allclose(sorted(vals), [1 - np.sqrt(e), 1 + np.sqrt(e)], rtol=1e-6, atol=1e-12)
angles, conds = np.array(angles), np.array(conds)
assert np.allclose(angles, 2 * np.arctan(np.sqrt(eps_vals)), rtol=1e-4)
assert np.allclose(conds, 1 / np.sqrt(eps_vals), rtol=1e-4)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
ax1.loglog(eps_vals, angles, "o-", label="angle between eigenvectors")
ax1.loglog(eps_vals, 2 * np.sqrt(eps_vals), "--", label=r"$2\sqrt{\varepsilon}$")
ax1.set_xlabel(r"$\varepsilon$"); ax1.set_ylabel("radians"); ax1.legend(); ax1.grid(alpha=.3)
ax2.loglog(eps_vals, conds, "o-", label=r"cond of eigenvector matrix $V$")
ax2.loglog(eps_vals, 1 / np.sqrt(eps_vals), "--", label=r"$1/\sqrt{\varepsilon}$")
ax2.set_xlabel(r"$\varepsilon$"); ax2.legend(); ax2.grid(alpha=.3)
fig.suptitle(r"$A_\varepsilon = [[1, 1], [\varepsilon, 1]]$: diagonalizable for every $\varepsilon \neq 0$, a Jordan block at $0$")
fig.tight_layout(); fig.savefig(FIG / "jordan_discontinuity.png", dpi=120); plt.close(fig)

# Exact arithmetic sees the switch too: eps = 10^-16 is "diagonalizable" to sympy.
_, J_tiny = sp.Matrix([[1, 1], [sp.Rational(1, 10 ** 16), 1]]).jordan_form()
assert J_tiny.is_diagonal()
_, J_zero = sp.Matrix([[1, 1], [0, 1]]).jordan_form()
assert not J_zero.is_diagonal()
# And floating point cannot tell: eig on the exact block returns two nearly parallel eigenvectors.
_, V_block = np.linalg.eig(np.array([[1.0, 1.0], [0.0, 1.0]]))
cond_block = np.linalg.cond(V_block)
assert cond_block > 1e15
print(f"Part 8: angle ~ 2 sqrt(eps), cond(V) = 1/sqrt(eps) exactly; cond(V) for the exact block in floating point = {cond_block:.1e}")

# %% [markdown]
# ## Your turn
#
# 1. **Three masses on a ring.** `K_ring` above is the stiffness matrix of three
#    equal masses joined in a ring by equal springs; it has frequency sqrt(3) twice.
#    Compute its eigenvectors with `K_ring.eigenvects()`, confirm the 3-eigenspace
#    is two-dimensional, and run `kernel_dims(K_ring, 3)`: the sequence stops at
#    2 immediately, so there is no generalized eigenvector. Then prove in a
#    comment why ker (A - lam I)^2 = ker (A - lam I) for every real symmetric A
#    (hint: ||(A - lam I) v||^2 = v^T (A - lam I)^2 v).
# 2. **A nilpotent 3x3 of index 2.** Any rank-one N = u v^T with v^T u = 0 works,
#    e.g. u = (1, 1, 2), v = (1, 1, -1). Check N^2 = 0, run `kernel_dims(N, 0)`
#    and `nilpotent_chains(N)`, and predict the Jordan form before you look.
# 3. **cos(A) for a 2x2 block.** Using `f_jordan_block(sp.cos(x), lam, 2)` you get
#    cos(lam) I - sin(lam) N. Now take A = P J_2(lam) P^-1 for your own P and
#    check `f_of_matrix(A, sp.cos(x))` against the power series
#    sum_k (-1)^k A^{2k} / (2k)! truncated at k = 12, in floating point.
# 4. **(Exercise B7)** Write `minimal_polynomial_by_dependence` yourself before
#    reading the one above, and make it return the degree-d coefficient vector
#    as well. On which of the test matrices is deg m_A < deg p_A?
# 5. **Break it.** Feed `jordan_block_sizes` a float matrix such as
#    `sp.Matrix([[1.0, 1.0], [1e-15, 1.0]])`. Explain what you see in terms of §3.11.

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
