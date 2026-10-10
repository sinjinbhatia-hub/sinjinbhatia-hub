# %% [markdown]
# # F1 Project: Proof Checkers and the Geometry of Complex Numbers
#
# This script checks every numerical claim in the F1 notes and builds the
# pictures that make complex numbers concrete. Parts: (1) a SymPy-verified
# induction identity checker, (2) the quantifier game on a finite set,
# (3) roots of unity and their vanishing sum, (4) Euler's formula numerically
# and via partial sums of the series, (5) complex multiplication as a 2x2
# rotation-scaling matrix, (6) phasor addition, (7) a damped oscillator solved
# by complex exponential and compared with scipy, (8) injective/surjective
# tests for finite functions, (9) the remaining worked-example numbers.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sympy as sp
from scipy.integrate import solve_ivp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)

# %% [markdown]
# ## Part 1: Induction identity checker (Theorem 3.8, Example 4.2)
#
# An induction proof has pieces a machine can check: the base case P(1), and the
# inductive step, which for a sum identity says closed(m+1) - closed(m) == term(m+1).
# SymPy verifies both symbolically and also computes the sum itself; exact integer
# arithmetic brute-forces n <= 50. Exercise B3 asks you to add sum k^3.

# %%
n, k = sp.symbols("n k", integer=True, positive=True)


def check_induction_identity(name, term, closed, n_max=50):
    """Verify sum_{k=1}^{n} term(k) == closed(n), where term is an expression in k
    and closed an expression in n. Four checks: (i) base case n = 1; (ii) the
    inductive step closed(n+1) - closed(n) - term(n+1) == 0 as a symbolic identity;
    (iii) SymPy's own summation equals closed; (iv) exact brute force for n <= n_max."""
    base = sp.simplify(term.subs(k, 1) - closed.subs(n, 1))
    step = sp.simplify(closed.subs(n, n + 1) - closed - term.subs(k, n + 1))
    symbolic = sp.simplify(sp.summation(term, (k, 1, n)) - closed)
    assert base == 0, f"{name}: base case fails ({base})"
    assert step == 0, f"{name}: inductive step fails ({step})"
    assert symbolic == 0, f"{name}: closed form is wrong ({symbolic})"
    running = sp.Integer(0)
    for m in range(1, n_max + 1):
        running += term.subs(k, m)
        assert running == closed.subs(n, m), f"{name}: fails at n = {m}"
    print(f"  {name}: base ok, inductive step ok, symbolic sum ok, brute force n <= {n_max} ok")


identities = {
    "sum k = n(n+1)/2": (k, n * (n + 1) / 2),
    "sum k^2 = n(n+1)(2n+1)/6": (k**2, n * (n + 1) * (2 * n + 1) / 6),
}
print("Part 1: induction identities")
for name, (term, closed) in identities.items():
    check_induction_identity(name, term, closed)

# The checker must reject a false closed form (otherwise it is a yes-machine).
try:
    check_induction_identity("wrong formula", k**2, n * (n + 1) * (2 * n + 1) / 5)
    raise RuntimeError("checker accepted a false identity")
except AssertionError:
    print("  a deliberately wrong closed form is rejected: ok")

# %% [markdown]
# ## Part 2: The quantifier game (§3.4, Theorem 3.6, rule (3.1))
#
# On the finite set X = {0,...,5} with a relation R(x, y), "for all x there is y"
# is a nested all/any; "there is y for all x" is any/all. They differ. Theorem 3.6
# says exists-forall implies forall-exists, never the reverse; the negation rule
# says not(forall x exists y R) == exists x forall y not R. Both are checked on
# 200 random relations.

# %%
X = list(range(6))


def forall_x_exists_y(R):
    """∀x ∃y R(x, y): the opponent picks x, you answer with y."""
    return all(any(R(x, y) for y in X) for x in X)


def exists_y_forall_x(R):
    """∃y ∀x R(x, y): you commit to y first, then the opponent picks any x."""
    return any(all(R(x, y) for x in X) for y in X)


def exists_x_forall_y(R):
    """∃x ∀y R(x, y) (used for the negation rule)."""
    return any(all(R(x, y) for y in X) for x in X)


def forall_y_exists_x(R):
    """∀y ∃x R(x, y) (used for the negation rule)."""
    return all(any(R(x, y) for x in X) for y in X)


R_sum = lambda x, y: (x + y) % 6 == 0        # every x has an additive inverse mod 6 ...
assert forall_x_exists_y(R_sum) and not exists_y_forall_x(R_sum)  # ... but no universal one
R_gt = lambda x, y: y > x                    # x = 5 has no y above it
assert not forall_x_exists_y(R_gt) and not exists_y_forall_x(R_gt)
R_ge = lambda x, y: y >= x                   # y = 5 works for everyone
assert forall_x_exists_y(R_ge) and exists_y_forall_x(R_ge)
print("Part 2: quantifier game")
print(f"  R(x,y): x+y = 0 mod 6   forall-exists = {forall_x_exists_y(R_sum)}, exists-forall = {exists_y_forall_x(R_sum)}")
print(f"  R(x,y): y > x           forall-exists = {forall_x_exists_y(R_gt)}, exists-forall = {exists_y_forall_x(R_gt)}")
print(f"  R(x,y): y >= x          forall-exists = {forall_x_exists_y(R_ge)}, exists-forall = {exists_y_forall_x(R_ge)}")

n_differ = 0
for _ in range(200):
    T = rng.random((6, 6)) < 0.5
    R = lambda x, y, T=T: bool(T[x, y])
    notR = lambda x, y, T=T: not T[x, y]
    fe, ef = forall_x_exists_y(R), exists_y_forall_x(R)
    assert (not ef) or fe, "Theorem 3.6 violated: exists-forall without forall-exists"
    n_differ += fe != ef
    # Negation rule (3.1), applied twice: ¬∀x∃y R ≡ ∃x∀y ¬R  and  ¬∃y∀x R ≡ ∀y∃x ¬R.
    assert (not fe) == exists_x_forall_y(notR), "negation rule (3.1) violated"
    assert (not ef) == forall_y_exists_x(notR), "negation rule (3.1) violated"
print(f"  200 random relations: Theorem 3.6 and rule (3.1) hold on all; the two orders differ on {n_differ}.")

# %% [markdown]
# ## Part 3: Roots of unity (Theorems 3.17, 3.18; exercise B7(d))
#
# omega_k = exp(2 pi i k / n). Each satisfies z^n = 1, they are distinct, they are
# the powers of omega_1, and they sum to zero for n >= 2.
#
# **Build it yourself first:** write `roots_of_unity(n)` in one line.

# %%
def roots_of_unity(n):
    """The n-th roots of unity omega_k = exp(2 pi i k / n), k = 0..n-1, eq. (3.15)."""
    return np.exp(2j * np.pi * np.arange(n) / n)


fig, axes = plt.subplots(2, 3, figsize=(10, 6.6))
circle_t = np.linspace(0, 2 * np.pi, 300)
for ax, nn in zip(axes.flat, range(3, 9)):
    w = roots_of_unity(nn)
    assert np.allclose(w**nn, 1), "not all roots satisfy z^n = 1"
    assert abs(w.sum()) < 1e-12, "roots of unity should sum to zero (Theorem 3.18)"
    assert len(np.unique(np.round(w, 9))) == nn, "roots should be distinct"
    assert np.allclose(w, w[1] ** np.arange(nn)), "omega_k should equal omega^k"
    ax.plot(np.cos(circle_t), np.sin(circle_t), color="0.75", lw=1)
    closed = np.append(w, w[0])
    ax.plot(closed.real, closed.imag, "-", color="C0", lw=1.5)
    ax.plot(w.real, w.imag, "o", color="C0", ms=6)
    ax.plot(0, 0, "k+", ms=8)
    ax.set_title(f"n = {nn}:  |sum| = {abs(w.sum()):.0e}", fontsize=10)
    ax.set_aspect("equal"); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3); ax.grid(alpha=.3)
    ax.set_xlabel("Re"); ax.set_ylabel("Im")
fig.suptitle("n-th roots of unity: vertices of a regular n-gon, centroid at the origin")
fig.tight_layout()
fig.savefig(FIG / "roots_of_unity.png", dpi=120)
plt.close(fig)

# Exercise B7(d): sum_k omega^{km} = n if n | m, else 0 (DFT orthogonality).
w6 = roots_of_unity(6)
for m in range(-12, 13):
    s = np.sum(w6**m)
    expected = 6 if m % 6 == 0 else 0
    assert abs(s - expected) < 1e-12, (m, s)
print("Part 3: roots of unity for n = 3..8 satisfy z^n = 1, are distinct, and sum to 0; B7(d) orthogonality ok for n = 6.")

# %% [markdown]
# ## Part 4: Euler's formula (Theorem 3.13)
#
# Numerically at random angles; via partial sums of the defining series
# (Derivation 1), with the error plotted against the number of terms; and via
# the ODE f' = i f (Derivation 2), checked by a central finite difference.

# %%
def exp_partial_sum(w, N):
    """S_N(w) = sum_{m=0}^{N} w^m / m!, the defining series of e^w truncated after N+1 terms.
    Terms are built by the recurrence term_{m+1} = term_m * w / (m+1)."""
    total, term = 0j, 1 + 0j
    for m in range(N + 1):
        total += term
        term = term * w / (m + 1)
    return total


theta = rng.uniform(-np.pi, np.pi, 1000)
assert np.allclose(np.exp(1j * theta), np.cos(theta) + 1j * np.sin(theta))
assert np.allclose(np.abs(np.exp(1j * theta)), 1)                     # |e^{i theta}| = 1
assert np.allclose(np.conj(np.exp(1j * theta)), np.exp(-1j * theta))  # conjugate flips the sign
assert abs(np.exp(1j * np.pi) + 1) < 1e-15                            # e^{i pi} = -1
# Theorem 3.14: e^{i a} e^{i b} = e^{i(a+b)}, so angles add and moduli multiply.
a, b = rng.uniform(-np.pi, np.pi, (2, 1000))
assert np.allclose(np.exp(1j * a) * np.exp(1j * b), np.exp(1j * (a + b)))
# (3.17): cos and sin as combinations of e^{±i theta}.
assert np.allclose((np.exp(1j * theta) + np.exp(-1j * theta)) / 2, np.cos(theta))
assert np.allclose((np.exp(1j * theta) - np.exp(-1j * theta)) / (2j), np.sin(theta))
# Derivation 2: f' = i f, by central differences.
h = 1e-6
deriv = (np.exp(1j * (theta + h)) - np.exp(1j * (theta - h))) / (2 * h)
assert np.allclose(deriv, 1j * np.exp(1j * theta), atol=1e-8)

Ns = np.arange(0, 31)
fig, ax = plt.subplots(figsize=(6.5, 4.2))
for color, th0, label in [("C0", 1.0, "θ = 1"), ("C1", 2.0, "θ = 2"), ("C2", np.pi, "θ = π")]:
    errs = np.array([abs(exp_partial_sum(1j * th0, N) - np.exp(1j * th0)) for N in Ns])
    assert errs[-1] < 1e-13, "series should have converged by 30 terms"
    assert errs[10] < 1e-2 and errs[20] < 1e-9
    ax.semilogy(Ns, np.maximum(errs, 1e-17), "-o", color=color, ms=3, lw=1.5, label=label)
ax.set_xlabel("N (terms 0..N kept)"); ax.set_ylabel("|S_N(iθ) − e^{iθ}|")
ax.set_title("Euler's formula from the series: partial-sum error")
ax.grid(alpha=.3, which="both"); ax.legend()
fig.tight_layout()
fig.savefig(FIG / "euler_series_error.png", dpi=120)
plt.close(fig)
print("Part 4: Euler's formula ok at 1000 random angles; series error < 1e-13 by N = 30; f' = i f ok.")

# %% [markdown]
# ## Part 5: Complex multiplication as a 2x2 matrix (Theorem 3.15)
#
# M(z) = [[a, -b], [b, a]] is the matrix of w -> z w on R^2. Checks: M(z) w == z w,
# M(z) M(w) == M(z w), det M(z) == |z|^2, M(z) == |z| R_theta, M(i)^2 == -I, and
# Example 4.6 as a matrix power.

# %%
def M(z):
    """Matrix of multiplication by z = a + bi acting on R^2 ≅ C, eq. (3.14)."""
    z = complex(z)
    return np.array([[z.real, -z.imag], [z.imag, z.real]])


def vec(w):
    """The column (Re w, Im w) that represents w in R^2."""
    w = complex(w)
    return np.array([w.real, w.imag])


def rotation(theta):
    return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])


for _ in range(200):
    z, w = (complex(*rng.standard_normal(2)) for _ in range(2))
    assert np.allclose(M(z) @ vec(w), vec(z * w))                 # the matrix does the multiplication
    assert np.allclose(M(z) @ M(w), M(z * w))                     # products compose
    assert np.isclose(np.linalg.det(M(z)), abs(z) ** 2)           # det = |z|^2
    assert np.allclose(M(z), abs(z) * rotation(np.angle(z)))      # rotation by arg z, scaling by |z|
assert np.allclose(M(1j) @ M(1j), -np.eye(2))                     # two quarter turns = -1
# Example 4.6: (1+i)^10 = 32 i, as a matrix power applied to 1 = (1, 0).
v = np.linalg.matrix_power(M(1 + 1j), 10) @ vec(1)
assert np.allclose(v, [0, 32]), v
print("Part 5: M(z) checks ok on 200 random pairs; M(1+i)^10 (1,0) =", np.round(v, 10), "i.e. (1+i)^10 = 32i.")

# %% [markdown]
# ## Part 6: Phasor addition (Theorem 3.20, Example 4.8)
#
# 3 cos(wt) + 4 cos(wt + pi/2) = 5 cos(wt + 0.9273). Add the phasors 3 and 4i as
# vectors, read off modulus and argument, and check against np.cos on a grid.

# %%
omega = 2 * np.pi          # 1 Hz, so one period is one time unit
t = np.linspace(0, 2, 801)
A1, p1, A2, p2 = 3.0, 0.0, 4.0, np.pi / 2
ph1, ph2 = A1 * np.exp(1j * p1), A2 * np.exp(1j * p2)
ph = ph1 + ph2
A, phi = abs(ph), np.angle(ph)
assert np.isclose(A, 5.0) and np.isclose(phi, np.arctan2(4, 3)) and abs(phi - 0.9273) < 5e-5
lhs = A1 * np.cos(omega * t + p1) + A2 * np.cos(omega * t + p2)
rhs = A * np.cos(omega * t + phi)
assert np.allclose(lhs, rhs), "phasor addition disagrees with the trig sum"
assert np.allclose(lhs, np.real(ph * np.exp(1j * omega * t)))   # eq. (3.18) for the sum
assert np.isclose(lhs[0], 3.0) and np.isclose(rhs[0], 3.0)       # the t = 0 check in Example 4.8

fig, (ax_ph, ax_t) = plt.subplots(1, 2, figsize=(11, 4.4))
# Left: the rotating phasors at wt = 0 (solid) and wt = 0.8 rad (faint).
for wt, alpha in [(0.0, 1.0), (0.8, 0.3)]:
    rot = np.exp(1j * wt)
    z1, z2, zs = ph1 * rot, ph2 * rot, ph * rot
    ax_ph.quiver(0, 0, z1.real, z1.imag, angles="xy", scale_units="xy", scale=1, color="C0", alpha=alpha, width=.008)
    ax_ph.quiver(z1.real, z1.imag, z2.real, z2.imag, angles="xy", scale_units="xy", scale=1, color="C1", alpha=alpha, width=.008)
    ax_ph.quiver(0, 0, zs.real, zs.imag, angles="xy", scale_units="xy", scale=1, color="C2", alpha=alpha, width=.008)
ax_ph.plot(5 * np.cos(circle_t), 5 * np.sin(circle_t), color="0.8", lw=1)
ax_ph.plot([], [], color="C0", label="3·e^{i0}"); ax_ph.plot([], [], color="C1", label="4·e^{iπ/2} (head to tail)")
ax_ph.plot([], [], color="C2", label=f"sum: 5·e^{{i{phi:.3f}}}")
ax_ph.set_aspect("equal"); ax_ph.set_xlim(-6, 6); ax_ph.set_ylim(-6, 6); ax_ph.grid(alpha=.3); ax_ph.legend(fontsize=8, loc="lower left")
ax_ph.set_title("Phasors at ωt = 0 (solid) and ωt = 0.8 (faint)"); ax_ph.set_xlabel("Re"); ax_ph.set_ylabel("Im")
# Right: the sinusoids (the real parts of the rotating phasors).
ax_t.plot(t, A1 * np.cos(omega * t + p1), color="C0", lw=1.5, label="3 cos ωt")
ax_t.plot(t, A2 * np.cos(omega * t + p2), color="C1", lw=1.5, label="4 cos(ωt + π/2)")
ax_t.plot(t, lhs, color="C2", lw=2.2, label=f"sum = 5 cos(ωt + {phi:.3f})")
ax_t.axhline(0, color="k", lw=.6)
ax_t.set_xlabel("t (periods)"); ax_t.set_ylabel("displacement"); ax_t.grid(alpha=.3); ax_t.legend(fontsize=8, loc="upper right")
ax_t.set_title("Same frequency in, same frequency out: amplitude 5, not 7")
fig.tight_layout()
fig.savefig(FIG / "phasor_addition.png", dpi=120)
plt.close(fig)
print(f"Part 6: 3cos(wt) + 4cos(wt + pi/2) = {A:.4f} cos(wt + {phi:.4f}); identity holds on {t.size} time samples.")

# %% [markdown]
# ## Part 7: Damped oscillator by complex exponential (Theorem 3.21, Example 4.9)
#
# x'' + x' + 4x = 0, x(0) = 1, x'(0) = 0. Ansatz e^{lambda t} gives
# lambda = -1/2 ± i sqrt(15)/2. The real solution built from Re and Im of e^{lambda t}
# is compared with scipy's solve_ivp.

# %%
gamma, w0 = 0.5, 2.0
lam = np.roots([1, 2 * gamma, w0**2])                   # roots of lambda^2 + 2 gamma lambda + w0^2
assert np.allclose(np.sort(lam.imag), [-np.sqrt(15) / 2, np.sqrt(15) / 2]) and np.allclose(lam.real, -gamma)
wd = np.sqrt(w0**2 - gamma**2)
assert np.isclose(wd, 1.9365, atol=5e-5)

# (3.19): d/dt e^{lambda t} = lambda e^{lambda t} for complex lambda, by central differences.
tt = np.linspace(0, 5, 200)
lam0 = lam[0]
d_num = (np.exp(lam0 * (tt + h)) - np.exp(lam0 * (tt - h))) / (2 * h)
assert np.allclose(d_num, lam0 * np.exp(lam0 * tt), atol=1e-7)


def x_exact(t):
    """Example 4.9: x(t) = e^{-t/2} (cos(w_d t) + (gamma / w_d) sin(w_d t)), the real
    combination of Re e^{lambda t} and Im e^{lambda t} matching x(0) = 1, x'(0) = 0."""
    return np.exp(-gamma * t) * (np.cos(wd * t) + (gamma / wd) * np.sin(wd * t))


# Symbolic check that x_exact satisfies the ODE exactly.
ts = sp.symbols("t", real=True)
xs = sp.exp(-ts / 2) * (sp.cos(sp.sqrt(15) / 2 * ts) + sp.sin(sp.sqrt(15) / 2 * ts) / sp.sqrt(15))
assert sp.simplify(sp.diff(xs, ts, 2) + sp.diff(xs, ts) + 4 * xs) == 0
assert xs.subs(ts, 0) == 1 and sp.simplify(sp.diff(xs, ts).subs(ts, 0)) == 0

sol = solve_ivp(lambda t, y: [y[1], -2 * gamma * y[1] - w0**2 * y[0]], (0, 10), [1.0, 0.0],
                rtol=1e-10, atol=1e-12, dense_output=True)
tg = np.linspace(0, 10, 1001)
x_num = sol.sol(tg)[0]
err = np.max(np.abs(x_num - x_exact(tg)))
assert err < 1e-6, err

fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.plot(tg, x_exact(tg), color="C0", lw=1.8, label="Re/Im of e^{λt} (Example 4.9)")
ax.plot(tg[::40], x_num[::40], "o", color="C1", ms=5, mfc="none", label="scipy solve_ivp")
ax.plot(tg, np.exp(-gamma * tg) * np.sqrt(1 + (gamma / wd) ** 2), "--", color="0.5", lw=1, label="envelope ±e^{−γt}·amplitude")
ax.plot(tg, -np.exp(-gamma * tg) * np.sqrt(1 + (gamma / wd) ** 2), "--", color="0.5", lw=1)
ax.axhline(0, color="k", lw=.6)
ax.set_xlabel("t"); ax.set_ylabel("x(t)"); ax.grid(alpha=.3); ax.legend(fontsize=8)
ax.set_title("x'' + x' + 4x = 0:  λ = −1/2 ± i√15/2  →  decay rate 1/2, frequency √15/2")
fig.tight_layout()
fig.savefig(FIG / "damped_oscillator.png", dpi=120)
plt.close(fig)
print(f"Part 7: lambda = {lam[0]:.4f}, {lam[1]:.4f}; w_d = {wd:.4f}; max |x_exact - solve_ivp| = {err:.1e}.")

# %% [markdown]
# ## Part 8: Injective and surjective tests for finite functions (Definition 3.3)
#
# A finite function is a dict {input: output}. Injective means no output repeats;
# surjective means every element of the codomain appears. Example 4.3 and exercise
# B2 are tested on random functions.
#
# **Build it yourself first:** write `is_injective(f)` and `is_surjective(f, codomain)`.

# %%
def is_injective(f):
    """Definition 3.3: f(a1) = f(a2) implies a1 = a2, i.e. no value of the dict repeats."""
    vals = list(f.values())
    return len(set(vals)) == len(vals)


def is_surjective(f, codomain):
    """Definition 3.3: every b in the codomain equals f(a) for some a."""
    return set(f.values()) == set(codomain)


def is_bijective(f, codomain):
    return is_injective(f) and is_surjective(f, codomain)


def compose(g, f):
    """(g ∘ f)(a) = g(f(a)); needs image(f) ⊆ domain(g)."""
    return {a: g[f[a]] for a in f}


def inverse(f):
    """The reversed dict. It is a function with f ∘ g = id and g ∘ f = id exactly when
    f is bijective (Theorem 3.4)."""
    return {b: a for a, b in f.items()}


def random_function(dom, cod):
    return {a: cod[rng.integers(len(cod))] for a in dom}


def random_injection(dom, cod):
    assert len(dom) <= len(cod)
    return dict(zip(dom, rng.permutation(cod)[: len(dom)].tolist()))


# The three pictures of §3.2.
f_inj = {1: "a", 2: "b"}                   # into {a, b, c}: injective, not surjective
f_sur = {1: "a", 2: "a", 3: "b"}           # onto {a, b}: surjective, not injective
f_bij = {1: "a", 2: "b", 3: "c"}           # onto {a, b, c}: bijective
assert is_injective(f_inj) and not is_surjective(f_inj, "abc")
assert is_surjective(f_sur, "ab") and not is_injective(f_sur)
assert is_bijective(f_bij, "abc")
# Theorem 3.4: the reversed dict is a two-sided inverse exactly for the bijection.
g_inv = inverse(f_bij)
assert compose(g_inv, f_bij) == {a: a for a in f_bij} and compose(f_bij, g_inv) == {b: b for b in "abc"}
assert not is_bijective(f_sur, "ab") and len(inverse(f_sur)) < len(f_sur)   # reversing f_sur loses a key: not a function on {1,2,3}

A_, B_, C_ = list(range(4)), list(range(6)), list(range(5))
for _ in range(200):                                   # Example 4.3: injective ∘ injective is injective
    f, g = random_injection(A_, B_), random_injection(B_, list(range(8)))
    assert is_injective(compose(g, f))
n_sur = 0
for _ in range(2000):                                  # B2(a): surjective ∘ surjective is surjective
    f, g = random_function(A_, A_), random_function(A_, A_)
    if is_surjective(f, A_) and is_surjective(g, A_):
        assert is_surjective(compose(g, f), A_)
        n_sur += 1
n_inj = 0
for _ in range(2000):                                  # B2(b): g∘f injective ⇒ f injective
    f, g = random_function(A_, B_), random_function(B_, C_)
    if is_injective(compose(g, f)):
        assert is_injective(f)
        n_inj += 1
# B2(c): g∘f injective does NOT force g injective.
f_c = {1: "a", 2: "b"}
g_c = {"a": 1, "b": 2, "c": 2}
assert is_injective(compose(g_c, f_c)) and not is_injective(g_c)
print(f"Part 8: finite-function tests ok ({n_sur} random surjective pairs, {n_inj} random injective compositions); B2(c) counterexample confirmed.")

# %% [markdown]
# ## Part 9: Remaining worked-example and exercise numbers
#
# Example 4.7 (cube roots of 8i), Example 4.5 (double sum), exercise A4, A5, B7(b),
# and the forced-oscillator phasor of exercise B8.

# %%
# Example 4.7: cube roots of 8i.
cube_roots = 2 * np.exp(1j * (np.pi / 6 + 2 * np.pi * np.arange(3) / 3))
assert np.allclose(cube_roots**3, 8j)
for expected in [np.sqrt(3) + 1j, -np.sqrt(3) + 1j, -2j]:
    assert np.min(np.abs(cube_roots - expected)) < 1e-12, expected
np_roots = np.roots([1, 0, 0, -8j])                     # numpy agrees
for r in np_roots:
    assert np.min(np.abs(cube_roots - r)) < 1e-9
fig, ax = plt.subplots(figsize=(4.4, 4.4))
ax.plot(2 * np.cos(circle_t), 2 * np.sin(circle_t), color="0.75", lw=1)
ax.plot(np.append(cube_roots, cube_roots[0]).real, np.append(cube_roots, cube_roots[0]).imag, "-", color="C0", lw=1.2)
ax.plot(cube_roots.real, cube_roots.imag, "o", color="C0", ms=7)
ax.plot(0, 8, "s", color="C1", ms=7)
for r, lab in zip(cube_roots, ["√3 + i", "−√3 + i", "−2i"]):
    ax.annotate(lab, (r.real, r.imag), textcoords="offset points", xytext=(8, 6), fontsize=9)
ax.annotate("8i", (0, 8), textcoords="offset points", xytext=(8, 0), fontsize=9)
ax.set_aspect("equal"); ax.set_xlim(-4, 4); ax.set_ylim(-3, 9); ax.grid(alpha=.3)
ax.set_xlabel("Re"); ax.set_ylabel("Im"); ax.set_title("Cube roots of 8i (Example 4.7)")
fig.tight_layout()
fig.savefig(FIG / "cube_roots_8i.png", dpi=120)
plt.close(fig)

# Example 4.5: the double sum both ways.
assert sum(1 for i in range(1, 4) for j in range(1, i + 1)) == 6
assert sum(1 for j in range(1, 4) for i in range(j, 4)) == 6
# Example 4.6 directly in polar form.
assert np.isclose((np.sqrt(2) * np.exp(1j * np.pi / 4)) ** 10, 32j)
# Exercise A4.
assert (2 + 3j) * (1 - 1j) == 5 + 1j and np.isclose((2 + 3j) / (1 - 1j), -0.5 + 2.5j)
assert np.isclose(abs(2 + 3j) ** 2, ((2 + 3j) * np.conj(2 + 3j)).real)
assert np.isclose((1 - np.sqrt(3) * 1j) ** 6, 64)
# Exercise A5.
assert sum(i + j for i in range(1, 5) for j in range(1, 4)) == 54
assert np.isclose(sum(1 / (q * (q + 1)) for q in range(1, 11)), 10 / 11)
# Exercise B7(b): fourth roots of -16.
for r in np.roots([1, 0, 0, 0, 16]):
    assert np.isclose(abs(r), 2) and np.isclose(abs(r.real), np.sqrt(2)) and np.isclose(abs(r.imag), np.sqrt(2))

# Exercise B8: forced oscillator m x'' + c x' + k x = F0 cos(w t), steady state x_p = Re[X e^{i w t}].
m_, c_, k_, F0, w_ = 1.0, 0.5, 4.0, 1.0, 1.5
Xp = F0 / (k_ - m_ * w_**2 + 1j * c_ * w_)
assert abs(abs(Xp) - 0.5252) < 5e-5 and abs(np.angle(Xp) - (-0.4049)) < 5e-5
tf = np.linspace(0, 20, 2001)
E = Xp * np.exp(1j * w_ * tf)
x_p, v_p, a_p = np.real(E), np.real(1j * w_ * E), np.real(-(w_**2) * E)     # derivatives via (3.19)
assert np.allclose(m_ * a_p + c_ * v_p + k_ * x_p, F0 * np.cos(w_ * tf)), "x_p does not satisfy the forced ODE"
print(f"Part 9: cube roots of 8i = {np.round(cube_roots, 6)}; double sum = 6; "
      f"forced oscillator |X| = {abs(Xp):.4f}, arg X = {np.angle(Xp):.4f} rad (lag).")

# %% [markdown]
# ## Your turn
#
# 1. **Binomial sum.** Prove sum_{k=0}^{n} C(n, k) = 2^n two ways: (i) by the binomial
#    theorem with x = y = 1; (ii) by induction using Pascal's rule
#    C(n+1, k) = C(n, k) + C(n, k-1). Then check it here with `sp.binomial` and
#    `sp.summation`, and by brute force for n <= 30 with `math.comb`.
# 2. **Regular polygon.** Show the n-th roots of unity are the vertices of a regular
#    n-gon (exercise C3) and compute its area with the shoelace formula
#    `0.5 * abs(sum(x_k y_{k+1} - x_{k+1} y_k))`. Assert it equals (n/2) sin(2 pi / n)
#    for n = 3..8 and watch it approach pi.
# 3. **A rotation matrix at work.** Write the matrix for multiplication by e^{i pi/4}
#    using `M(...)` from Part 5 and apply it to the four corners of the unit square
#    [0,1]^2 (as columns of a 2x4 array). Plot the square before and after. Apply it
#    eight times: what comes back, and why does M(e^{i pi/4})^8 equal the identity?
# 4. **Exercise B3.** Add `sum k^3 = (n(n+1)/2)^2` to `identities` in Part 1.
# 5. **Exercise C2(c).** Add a brute-force check that every 12 <= n <= 200 is
#    4a + 5b with a, b >= 0, and that 11 is not.

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
