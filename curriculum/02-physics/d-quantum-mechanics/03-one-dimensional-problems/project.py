# %% [markdown]
# # QM-3 Project: one-dimensional problems, solved three ways
# Three tools cover every 1D problem in the notes: (1) a finite-difference
# eigensolver (the Hamiltonian as a symmetric tridiagonal matrix, LA-8's
# spectral theorem, scipy's eigh_tridiagonal); (2) the oscillator's ladder
# operators as matrices in the number basis; (3) a transfer-matrix scattering
# code for any piecewise-constant potential plus delta functions (§3.9). A
# split-operator packet (QM-2 (7.1)) checks stationary scattering against
# time-dependent scattering. Every number quoted in notes §4 and in
# solutions.md is asserted or printed here. Units: hbar = m = 1 unless SI.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
import scipy.constants as C
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import hermite, airy, ai_zeros, factorial
from pathlib import Path
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)

# %% [markdown]
# ## Part 1: the finite-difference eigensolver (§7)
# On a grid x_j = x_0 + j h with psi = 0 at both ends, -psi''/2 + V psi = E psi
# becomes H psi = E psi with H_jj = 1/h^2 + V_j, H_{j,j+1} = -1/(2h^2): a real
# symmetric tridiagonal matrix, so its eigenvalues are real and its
# eigenvectors orthogonal (LA-8). Error O(h^2) for smooth V.

# %%
def fd_eigs(x, V, n_states):
    """Lowest n_states eigenpairs of H = -(1/2) d^2/dx^2 + V on interior grid x (Dirichlet ends)."""
    h = x[1] - x[0]
    d = 1.0 / h**2 + V
    e = -0.5 / h**2 * np.ones(len(x) - 1)
    E, vec = eigh_tridiagonal(d, e, select="i", select_range=(0, n_states - 1))
    return E, vec / np.sqrt(h)


# Infinite well on [0, 1]: E_n = n^2 pi^2 / 2, eq. (3.2)
N = 2000
xw = np.linspace(0, 1, N + 2)[1:-1]
Ew, vw = fd_eigs(xw, np.zeros(N), 6)
n = np.arange(1, 7)
assert np.allclose(Ew, n**2 * np.pi**2 / 2, rtol=1e-5)
# node theorem: psi_n has n-1 interior nodes
for j in range(6):
    s = np.sign(vw[np.abs(vw[:, j]) > 1e-8 * np.abs(vw[:, j]).max(), j])
    assert np.sum(s[1:] != s[:-1]) == j

# Oscillator V = x^2/2: E_n = n + 1/2, eq. (3.12)
xo = np.linspace(-10, 10, 3001)[1:-1]
Eo, vo = fd_eigs(xo, 0.5 * xo**2, 10)
assert np.allclose(Eo, np.arange(10) + 0.5, atol=5e-4)
print("FD: infinite well and oscillator eigenvalues match (3.2), (3.12)")

# %% [markdown]
# ## Part 2: the box with Psi(0) ∝ x(a - x) (Example 4.1)
# c_n = 8 sqrt(15)/(n pi)^3 for odd n, 0 for even n; P(E_1) = 960/pi^6.

# %%
def psi_box(nn, x, a=1.0):
    return np.sqrt(2 / a) * np.sin(nn * np.pi * x / a)


psi0 = lambda x: np.sqrt(30) * x * (1 - x)
cn = np.array([quad(lambda x: psi_box(k, x) * psi0(x), 0, 1)[0] for k in range(1, 40)])
cn_formula = np.array([8 * np.sqrt(15) / (k * np.pi)**3 if k % 2 else 0.0 for k in range(1, 40)])
assert np.allclose(cn, cn_formula, atol=1e-12)
P1 = 960 / np.pi**6
assert np.isclose(cn[0]**2, P1)
print(f"Ex 4.1: c_1 = {cn[0]:.6f}, P(E1) = {P1:.6f}, P(E3) = {cn[2]**2:.6f}, sum|c|^2 (39 terms) = {np.sum(cn**2):.8f}")
kk = np.arange(1, 4000, 2)
H_mean = np.sum((8 * np.sqrt(15) / (kk * np.pi)**3)**2 * kk**2 * np.pi**2 / 2)
assert np.isclose(H_mean, 5.0, rtol=1e-6)
print(f"Ex 4.1: <H> = sum |c_n|^2 E_n = {H_mean:.6f} (exact 5 hbar^2/ma^2); E_1 = {np.pi**2/2:.4f}")
# revival time T_rev = 4 m a^2 / (pi hbar): every phase e^{-i E_n T/hbar} returns to 1
T_rev = 4 / np.pi
ph = np.exp(-1j * kk**2 * np.pi**2 / 2 * T_rev)
assert np.allclose(ph, ph[0]) and np.isclose(abs(ph[0]), 1)
ov = np.sum((8 * np.sqrt(15) / (kk * np.pi)**3)**2 * np.exp(-1j * kk**2 * np.pi**2 / 2 * T_rev / 2))
ov12 = 0.5 * abs(np.exp(-1j * np.pi**2 / 2 * T_rev / 2) + np.exp(-1j * 4 * np.pi**2 / 2 * T_rev / 2))
assert np.isclose(abs(ov), 1) and np.isclose(ov12, 0, atol=1e-12)
print(f"revival: |<Psi(0)|Psi(T_rev)>| = 1; at T_rev/2: odd-n-only x(a-x) state {abs(ov):.4f}, (psi1+psi2)/sqrt2 {ov12:.4f}")

# %% [markdown]
# ## Part 3: the oscillator. Hermite polynomials, ladders, coherent states
# Recurrence H_{n+1} = 2 xi H_n - 2n H_{n-1} (3.17) vs scipy; the n = 2
# termination (Example 4.3) in SymPy; ladder matrices (3.13).

# %%
def hermite_rec(nmax, xi):
    H = [np.ones_like(xi), 2 * xi]
    for k in range(1, nmax):
        H.append(2 * xi * H[k] - 2 * k * H[k - 1])
    return H[: nmax + 1]


xi = np.linspace(-3, 3, 101)
Hr = hermite_rec(10, xi)
for k in range(11):
    assert np.allclose(Hr[k], hermite(k)(xi), rtol=1e-10, atol=1e-8)

xs, K = sp.symbols("xi K")
a0 = sp.Symbol("a0")
coef = {0: a0, 1: 0}
for j in range(0, 6):
    coef[j + 2] = sp.simplify((2 * j + 1 - K) / ((j + 1) * (j + 2)) * coef[j])
h = sum(coef[j] * xs**j for j in range(8))
h2 = sp.expand(h.subs(K, 5))                        # K = 2E/hbar w = 5 kills a_4 onward
assert sp.simplify(h2 - a0 * (1 - 2 * xs**2)) == 0
assert sp.simplify(sp.diff(h2, xs, 2) - 2 * xs * sp.diff(h2, xs) + 4 * h2) == 0
assert sp.expand(-2 * h2.subs(a0, 1)) == sp.expand(sp.hermite(2, xs))
print("Ex 4.3: K = 5 terminates the series at xi^2: h = a0 (1 - 2 xi^2) = -(a0/2) H_2")


def hermite_functions(nmax, x):
    """Normalized psi_n(x) (hbar = m = omega = 1) by the stable recurrence from a_+ (3.13)."""
    psi = [np.pi**-0.25 * np.exp(-x**2 / 2), np.pi**-0.25 * np.sqrt(2) * x * np.exp(-x**2 / 2)]
    for k in range(1, nmax):
        psi.append(np.sqrt(2 / (k + 1)) * x * psi[k] - np.sqrt(k / (k + 1)) * psi[k - 1])
    return np.array(psi[: nmax + 1])


hf = hermite_functions(9, xo)
for k in range(10):
    ref = (2**k * factorial(k) * np.sqrt(np.pi))**-0.5 * hermite(k)(xo) * np.exp(-xo**2 / 2)
    assert np.allclose(hf[k], ref, atol=1e-10)
    assert np.isclose(abs(np.sum(hf[k] * vo[:, k]) * (xo[1] - xo[0])), 1, atol=1e-4)

NB = 60
am = np.diag(np.sqrt(np.arange(1, NB)), 1)          # a_- psi_n = sqrt(n) psi_{n-1}
ap = am.T
comm = am @ ap - ap @ am
assert np.allclose(comm[:-1, :-1], np.eye(NB - 1))   # last entry is a truncation artifact
assert np.isclose(comm[-1, -1], -(NB - 1))
Hm = ap @ am + 0.5 * np.eye(NB)
assert np.allclose(Hm, np.diag(np.arange(NB) + 0.5))
X = (ap + am) / np.sqrt(2)
P = 1j * (ap - am) / np.sqrt(2)
X2, P2 = (X @ X).real, (P @ P).real
nn = np.arange(NB - 2)
assert np.allclose(np.diag(X2)[nn], nn + 0.5) and np.allclose(np.diag(P2)[nn], nn + 0.5)
assert np.allclose(np.diag(X)[nn], 0)
sxsp = np.sqrt(np.diag(X2) * np.diag(P2))[:4]
print(f"Ex 4.2: sigma_x sigma_p / hbar for n = 0..3: {np.round(sxsp, 6)} = n + 1/2; <T> = <V> = (n+1/2)/2")
# B3: (psi0 + psi1)/sqrt2: <x>(t) = sqrt(hbar/2 m w) cos wt
c = np.zeros(NB, complex); c[:2] = 1 / np.sqrt(2)
for t in np.linspace(0, 7, 15):
    ct = c * np.exp(-1j * (np.arange(NB) + 0.5) * t)
    assert np.isclose((ct.conj() @ X @ ct).real, np.cos(t) / np.sqrt(2))

# Coherent state |alpha>, alpha = 2: <x>(t) = sqrt2 Re(alpha e^{-it}), sigma_x = 1/sqrt2 forever
alpha = 2.0
cA = np.exp(-alpha**2 / 2) * alpha**np.arange(NB) / np.sqrt(factorial(np.arange(NB)))
assert np.allclose(am @ cA, alpha * cA, atol=1e-8)
ts = np.linspace(0, 4 * np.pi, 200)
xm, sx = [], []
for t in ts:
    ct = cA * np.exp(-1j * (np.arange(NB) + 0.5) * t)
    m1 = (ct.conj() @ X @ ct).real
    xm.append(m1); sx.append(np.sqrt((ct.conj() @ X @ X @ ct).real - m1**2))
assert np.allclose(xm, np.sqrt(2) * alpha * np.cos(ts), atol=1e-9)
assert np.allclose(sx, 1 / np.sqrt(2), atol=1e-9)
psi_n = hermite_functions(NB - 1, xo)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for k in range(5):
    ax[0].plot(xo, hf[k] + k + 0.5, lw=1)
ax[0].plot(xo, 0.5 * xo**2, "k", lw=0.8); ax[0].set_xlim(-5, 5); ax[0].set_ylim(0, 6)
ax[0].set_title("oscillator eigenfunctions, offset by E_n")
for t in np.linspace(0, np.pi, 5):
    dens = np.abs((cA * np.exp(-1j * (np.arange(NB) + 0.5) * t)) @ psi_n)**2
    ax[1].plot(xo, dens, label=f"t = {t:.2f}")
ax[1].set_xlim(-5, 5); ax[1].legend(); ax[1].set_title("coherent state alpha = 2: rigid oscillation")
fig.savefig(FIG / "oscillator.png", dpi=120); plt.close(fig)
print("ladders: [a-,a+] = 1, H diagonal, coherent state oscillates without spreading")

# %% [markdown]
# ## Part 4: the finite well (§3.6, Example 4.5)
# Width 2a (Griffiths' -a < x < a), z0 = a sqrt(2 m V0)/hbar. Even: tan z = sqrt(z0^2/z^2 - 1),
# odd: -cot z = sqrt(z0^2/z^2 - 1). Written as z sin z - w cos z = 0 and z cos z + w sin z = 0,
# w = sqrt(z0^2 - z^2), so no poles.

# %%
def finite_well_z(z0):
    w = lambda z: np.sqrt(max(z0**2 - z**2, 0.0))
    fe = lambda z: z * np.sin(z) - w(z) * np.cos(z)
    fo = lambda z: z * np.cos(z) + w(z) * np.sin(z)
    zg = np.linspace(1e-9, z0, 20001)
    roots = []
    for f, par in ((fe, "even"), (fo, "odd")):
        vals = np.array([f(z) for z in zg])
        for i in np.nonzero(np.sign(vals[1:]) != np.sign(vals[:-1]))[0]:
            roots.append((brentq(f, zg[i], zg[i + 1], xtol=1e-14), par))
    return sorted(roots)


z0 = 8.0
roots = finite_well_z(z0)
zs = np.array([r[0] for r in roots])
assert len(roots) == int(2 * z0 / np.pi) + 1 == 6
assert [r[1] for r in roots] == ["even", "odd"] * 3
EV = zs**2 / z0**2 - 1                               # E / V0
print("Ex 4.5: z_n =", np.round(zs, 4))
print("        E_n/V0 =", np.round(EV, 4))
print(f"        top state kappa a = {np.sqrt(z0**2 - zs[-1]**2):.4f}; n pi/(2 z_1) = {np.pi/2/zs[0]:.4f}")
print("        infinite-well estimate z = n pi/2:", np.round(np.arange(1, 7) * np.pi / 2, 4))
# FD check with a = 1, V0 = z0^2/2 = 32
V0w = z0**2 / 2
xf = np.linspace(-15, 15, 15001)[1:-1]
Vf = np.where(np.abs(xf) < 1, -V0w, 0.0)
Vf[np.isclose(np.abs(xf), 1)] = -V0w / 2            # wall on a grid point: average, restores O(h^2)
Ef, vf = fd_eigs(xf, Vf, 6)
assert np.allclose(Ef + V0w, zs**2 / 2, rtol=1e-4)
print("        FD energies", np.round(Ef, 3), "vs", np.round(V0w * EV, 3))
for z0s in (0.1, 0.5, 1.5):                          # shallow: exactly one bound state
    assert len(finite_well_z(z0s)) == 1
zz = np.linspace(0.01, z0, 2000)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
tz = np.tan(zz); tz[np.abs(tz) > 12] = np.nan
cz = -1 / np.tan(zz); cz[np.abs(cz) > 12] = np.nan
ax[0].plot(zz, tz, label="tan z"); ax[0].plot(zz, cz, label="-cot z")
ax[0].plot(zz, np.sqrt((z0 / zz)**2 - 1), "k", label="sqrt(z0²/z² - 1)")
ax[0].plot(zs, np.sqrt((z0 / zs)**2 - 1), "ro"); ax[0].set_ylim(0, 12); ax[0].legend()
ax[0].set_title("finite well z0 = 8: six bound states")
ax[1].plot(xf, np.where(np.abs(xf) < 1, -V0w, 0), "k", lw=0.8)
for j in range(6):
    ax[1].plot(xf, Ef[j] + 3 * vf[:, j], lw=1)
ax[1].set_xlim(-3, 3); ax[1].set_title("eigenfunctions offset by E_n")
fig.savefig(FIG / "finite_well.png", dpi=120); plt.close(fig)

# Double well: two wells (width 1, depth 10) separated by barrier b. Splitting ~ e^{-kappa b}
bs = np.linspace(1.0, 3.0, 9)
split, Ebar = [], []
xd = np.linspace(-10, 10, 8001)[1:-1]
for b in bs:
    Vd = np.where((np.abs(xd) > b / 2) & (np.abs(xd) < b / 2 + 1), -10.0, 0.0)
    Ed, vd = fd_eigs(xd, Vd, 2)
    split.append(Ed[1] - Ed[0]); Ebar.append(Ed.mean())
    assert vd[len(xd) // 2, 0] * vd[len(xd) // 2 + 200, 0] > 0      # ground even (no node)
slope = np.polyfit(bs[4:], np.log(split[4:]), 1)[0]
kap = np.sqrt(-2 * np.mean(Ebar[4:]))
assert abs(slope / -kap - 1) < 0.05
print(f"double well: d ln(dE)/db = {slope:.3f}, -kappa = {-kap:.3f}")
print(f"ammonia inversion: nu = 23.87 GHz -> dE = h nu = {C.h*23.87e9/C.e:.3e} eV")

# %% [markdown]
# ## Part 5: the transfer matrix (§3.9)
# Region j: psi = A_j e^{ik_j (x - x_{j-1})} + B_j e^{-ik_j (x - x_{j-1})}, k_j = sqrt(2(E - V_j))
# (imaginary when E < V_j). Interface x_j: (psi, psi') continuous, except psi' jumps by
# 2 alpha_j psi at a delta alpha_j delta(x - x_j). Step matrix S_j = W_{j+1}(0)^{-1} D_j W_j(d_j),
# W(x) = [[e^{ikx}, e^{-ikx}], [ik e^{ikx}, -ik e^{-ikx}]]. M = S_N...S_1 (LA-4's composition).
# With B_N = 0: r = -M21/M22, det M = k_0/k_N, T = (k_0/k_N)/|M22|^2 (3.31).

# %%
def transfer(E, xb, Vs, deltas=None):
    """T(E), R(E) for regions Vs[0..N] separated at xb[0..N-1], optional delta strengths (3.31)."""
    E = np.atleast_1d(np.asarray(E, float))
    deltas = np.zeros(len(xb)) if deltas is None else np.asarray(deltas, float)
    k = [np.sqrt(2 * (E - V) + 0j) for V in Vs]
    M = np.broadcast_to(np.eye(2, dtype=complex), (len(E), 2, 2)).copy()
    for j in range(len(xb)):
        d = 0.0 if j == 0 else xb[j] - xb[j - 1]
        kl, kr = k[j], k[j + 1]
        el, ei = np.exp(1j * kl * d), np.exp(-1j * kl * d)
        Wl = np.stack([np.stack([el, ei], -1), np.stack([1j * kl * el, -1j * kl * ei], -1)], -2)
        Dm = np.array([[1, 0], [2 * deltas[j], 1]], complex)
        Wr_inv = np.stack([np.stack([0.5 * np.ones_like(kr), 0.5 / (1j * kr)], -1),
                           np.stack([0.5 * np.ones_like(kr), -0.5 / (1j * kr)], -1)], -2)
        M = Wr_inv @ (Dm @ (Wl @ M))
    k0, kN = k[0], k[-1]
    T = np.where(np.abs(kN.imag) < 1e-12, (k0.real / np.where(kN.real > 0, kN.real, 1)) / np.abs(M[:, 1, 1])**2, 0.0)
    R = np.abs(M[:, 1, 0] / M[:, 1, 1])**2
    return T, R, M


def T_barrier(E, V0, a):
    """Closed form (3.26): T = [1 + V0^2 sinh^2(kappa a)/(4E(V0 - E))]^{-1}, E < V0."""
    kap = np.sqrt(2 * (V0 - E))
    return 1 / (1 + V0**2 * np.sinh(kap * a)**2 / (4 * E * (V0 - E)))


Es = np.linspace(0.01, 0.99, 300)
Tb, Rb, Mb = transfer(Es, [0, 1.5], [0, 1, 0])
assert np.allclose(Tb, T_barrier(Es, 1, 1.5), rtol=1e-10) and np.allclose(Tb + Rb, 1)
assert np.allclose(np.linalg.det(Mb), 1)                                  # det M = k0/kN = 1
Eh = np.linspace(1.01, 6, 300)                                            # E > V0: sin replaces sinh
Th = transfer(Eh, [0, 1.5], [0, 1, 0])[0]
assert np.allclose(Th, 1 / (1 + np.sin(np.sqrt(2 * (Eh - 1)) * 1.5)**2 / (4 * Eh * (Eh - 1))))

# SymPy: the matching algebra of §3.8 gives F/A = e^{-ika}/[cosh ka + i (kap^2-k^2)/(2 k kap) sinh]
k_, q_, a_ = sp.symbols("k kappa a", positive=True)
B_, C_, D_, F_, x_ = sp.symbols("B C D F x")
L_ = sp.exp(sp.I * k_ * x_) + B_ * sp.exp(-sp.I * k_ * x_)
Mi_ = C_ * sp.exp(q_ * x_) + D_ * sp.exp(-q_ * x_)
R_ = F_ * sp.exp(sp.I * k_ * x_)
sol = sp.solve([(L_ - Mi_).subs(x_, 0), sp.diff(L_ - Mi_, x_).subs(x_, 0),
                (Mi_ - R_).subs(x_, a_), sp.diff(Mi_ - R_, x_).subs(x_, a_)], [B_, C_, D_, F_], dict=True)[0]
target = sp.cosh(q_ * a_) + sp.I * (q_**2 - k_**2) / (2 * k_ * q_) * sp.sinh(q_ * a_)
assert sp.simplify((sp.exp(-sp.I * k_ * a_) / sol[F_] - target).rewrite(sp.exp)) == 0
print("barrier: transfer matrix = closed form (3.26) = SymPy elimination; R + T = 1")

# Step, Example 4.7: E = 2 V0
Ts, Rs, Ms = transfer([2.0], [0], [0, 1])
assert np.isclose(Rs[0], 17 - 12 * np.sqrt(2)) and np.isclose(Ts[0], 12 * np.sqrt(2) - 16)
naive = np.abs(np.linalg.det(Ms[0]) / Ms[0, 1, 1])**2
print(f"Ex 4.7: step E = 2V0: R = {Rs[0]:.4f}, T = {Ts[0]:.4f}; naive |F/A|^2 = {naive:.4f} (> 1!)")
Ts43 = transfer([4 / 3], [0], [0, 1])
assert np.isclose(Ts43[1][0], 1 / 9) and np.isclose(Ts43[0][0], 8 / 9)          # exercise A4
naive43 = np.abs(np.linalg.det(Ts43[2][0]) / Ts43[2][0, 1, 1])**2
assert np.isclose(naive43, 16 / 9)
Tdn = transfer([1.0], [0], [0, -3])                                             # exercise B8: step DOWN, k' = 2k
assert np.isclose(Tdn[1][0], 1 / 9) and np.isclose(Tdn[0][0], 8 / 9)
print(f"A4: E = 4V0/3: R = {Ts43[1][0]:.4f}, T = {Ts43[0][0]:.4f}, naive |F/A|^2 = {naive43:.4f}; B8: step down V = -3E: R = {Tdn[1][0]:.4f}")

# Delta well, Example 4.4: V = -alpha delta, alpha = 1: E_b = -1/2, T(1/2) = 1/2; +alpha same T
Td = transfer([0.5], [0], [0, 0], deltas=[-1.0])[0][0]
Tdp = transfer([0.5], [0], [0, 0], deltas=[1.0])[0][0]
Eg = np.linspace(0.05, 5, 50)
assert np.allclose(transfer(Eg, [0], [0, 0], deltas=[-1.0])[0], 1 / (1 + 1 / (2 * Eg)))
assert np.isclose(Td, 0.5) and np.isclose(Tdp, 0.5)
print(f"Ex 4.4: delta well: T(E = |E_b|) = {Td:.4f}; barrier +alpha: {Tdp:.4f}")

# Finite-well resonances: width w = 2, V0 = 10: T = 1 at E = n^2 pi^2/(2 w^2) - V0 > 0
w_, V0_ = 2.0, 10.0
nres = np.arange(1, 12)
Eres = nres**2 * np.pi**2 / (2 * w_**2) - V0_
Eres = Eres[Eres > 0]
assert np.allclose(transfer(Eres, [0, w_], [0, -V0_, 0])[0], 1, atol=1e-12)
print("well resonances (Ramsauer-Townsend) T = 1 at E =", np.round(Eres, 3))

# ln T vs width: slope -2 kappa (E = 0.5, V0 = 1 -> kappa = 1)
wid = np.linspace(4, 8, 9)
lnT = [np.log(transfer([0.5], [0, w], [0, 1, 0])[0][0]) for w in wid]
sl = np.polyfit(wid, lnT, 1)[0]
assert abs(sl + 2) < 1e-3
print(f"d lnT / d(width) = {sl:.6f} (theory -2 kappa = -2)")

# Double barrier (resonant tunneling): barriers V0 = 5, width 1, well of width 3 between
xb2, V2 = [0, 1, 4, 5], [0, 5, 0, 5, 0]
Eg2 = np.linspace(0.05, 4.95, 20000)
T2 = transfer(Eg2, xb2, V2)[0]
i0 = np.argmax(T2[Eg2 < 1.5])
opt = minimize_scalar(lambda e: -transfer([e], xb2, V2)[0][0], bracket=(Eg2[i0 - 3], Eg2[i0], Eg2[i0 + 3]))
Tsingle = transfer([opt.x], [0, 1], [0, 5, 0])[0][0]
assert -opt.fun > 0.999 and Tsingle < 0.05
print(f"double barrier: resonance at E = {opt.x:.5f}, T = {-opt.fun:.6f}; single barrier there T = {Tsingle:.4f}")

# Dirac comb: V = alpha sum delta(x - j a), a = 1, beta = m alpha a/hbar^2 = 3 (3.35)
beta, Ncell = 3.0, 12
Ec = np.linspace(0.05, 60, 12000)
kc = np.sqrt(2 * Ec)
fB = np.cos(kc) + beta * np.sin(kc) / kc
xcomb = np.arange(Ncell, dtype=float)
Tc = transfer(Ec, xcomb, [0.0] * (Ncell + 1), deltas=[beta] * Ncell)[0]
assert np.all(Tc[np.abs(fB) > 1.1] < 1e-2)
band1 = (np.abs(fB) <= 1) & (Ec < 10)
assert Tc[band1].max() > 0.99
edges = Ec[1:][np.diff((np.abs(fB) <= 1).astype(int)) != 0]
print("Dirac comb beta = 3: band edges (E) =", np.round(edges[:6], 3))
fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
ax[0].plot(Es, Tb, label="barrier V0=1, a=1.5"); ax[0].plot(Eh, Th)
ax[0].plot(Eg2, T2, label="double barrier V0=5 (resonant)", lw=0.8); ax[0].set_xlim(0, 6)
ax[0].legend(); ax[0].set_ylabel("T(E)")
ax[1].plot(Ec, fB); ax[1].axhline(1, c="k", lw=0.5); ax[1].axhline(-1, c="k", lw=0.5)
ax[1].set_ylim(-3, 4); ax[1].set_ylabel("cos ka + beta sin(ka)/ka")
ax[2].plot(Ec, Tc); ax[2].set_ylabel(f"T, {Ncell}-delta comb"); ax[2].set_xlabel("E")
for a_ax in ax[1:]:
    a_ax.set_xlim(0, 60)
fig.savefig(FIG / "transmission.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 6: SI numbers. Example 4.6 (electron barrier), STM, field emission, Gamow (Ex 4.8)

# %%
hb, me, eV = C.hbar, C.m_e, C.e


def T_SI(E_eV, V0_eV, a_m, m=me):
    kap = np.sqrt(2 * m * (V0_eV - E_eV) * eV) / hb
    return 1 / (1 + V0_eV**2 * np.sinh(kap * a_m)**2 / (4 * E_eV * (V0_eV - E_eV))), kap


T05, kap = T_SI(5, 10, 0.5e-9)
T10, _ = T_SI(5, 10, 1.0e-9)
approx = 16 * 0.5 * 0.5 * np.exp(-2 * kap * 0.5e-9)
print(f"Ex 4.6: kappa = {kap:.4e} 1/m, kappa a = {kap*0.5e-9:.3f}; T(0.5 nm) = {T05:.4e} "
      f"(thick approx {approx:.4e}); T(1 nm) = {T10:.4e}; ratio = {T05/T10:.4e}")
assert np.isclose(T05, 4.235e-5, rtol=1e-3) and np.isclose(T10, 4.485e-10, rtol=1e-3)
kST = np.sqrt(2 * me * 4.5 * eV) / hb
print(f"STM: 0.01 A of height changes the current by a factor {np.exp(2*kST*1e-12):.4f}")
print(f"STM: phi = 4.5 eV, kappa = {kST:.4e} 1/m; current drops by e^(2 kappa x 1 A) = {np.exp(2*kST*1e-10):.2f} per angstrom")
for Ffield in (5e9, 2.5e9):
    expo = 4 / 3 * np.sqrt(2 * me) * (4.5 * eV)**1.5 / (hb * eV * Ffield)
    print(f"field emission (Fowler-Nordheim exponent), F = {Ffield:.1e} V/m: T ~ exp(-{expo:.2f}) = {np.exp(-expo):.2e}")

m_al = 6.6447e-27            # alpha mass, kg
MeV, fm = 1e6 * eV, 1e-15
ke = C.e**2 / (4 * np.pi * C.epsilon_0)


def gamow_halflife(E_MeV, Zd, Ad):
    """t_1/2 = ln2 (2 r1/v) e^{2 gamma}, gamma = (1/hbar) int_{r1}^{r2} sqrt(2m(V - E)) dr (3.30)."""
    E = E_MeV * MeV
    r1 = 1.2 * fm * (Ad**(1 / 3) + 4**(1 / 3))      # touching radius: daughter + alpha
    r2 = 2 * Zd * ke / E
    u = r1 / r2
    gamma = np.sqrt(2 * m_al * E) / hb * r2 * (np.arccos(np.sqrt(u)) - np.sqrt(u * (1 - u)))
    v = np.sqrt(2 * E / m_al)
    return np.log(2) * 2 * r1 / v * np.exp(2 * gamma), gamma


# check the closed-form Gamow integral against quadrature
Eu, Zu = 4.27 * MeV, 90
r1u, r2u = 1.2 * fm * (234**(1 / 3) + 4**(1 / 3)), 2 * 90 * ke / Eu
gq = quad(lambda r: np.sqrt(2 * m_al * (2 * Zu * ke / r - Eu)), r1u, r2u)[0] / hb
assert np.isclose(gq, gamow_halflife(4.27, 90, 234)[1], rtol=2e-5)   # quad struggles with the sqrt endpoint
yr = 365.25 * 24 * 3600
tU, gU = gamow_halflife(4.27, 90, 234)
tPo, gPo = gamow_halflife(8.95, 82, 208)
print(f"Ex 4.8: U-238 (E = 4.27 MeV): gamma = {gU:.1f}, t1/2 = {tU/yr:.2e} yr (measured 4.5e9 yr)")
print(f"        Po-212 (E = 8.95 MeV): gamma = {gPo:.1f}, t1/2 = {tPo:.2e} s (measured 3.0e-7 s)")
print(f"        model ratio = 10^{np.log10(tU/tPo):.1f}; measured ratio = 10^{np.log10(4.47e9*yr/2.99e-7):.1f}")
assert 22 < np.log10(tU / tPo) < 26
Egr = np.linspace(4, 9, 100)
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(Egr, np.log10([gamow_halflife(e, 86, 218)[0] for e in Egr]))
ax.set_xlabel("alpha energy (MeV)"); ax.set_ylabel("log10 t1/2 (s), Z_d = 86")
ax.set_title("Gamow: a factor of 2 in E, ~24 decades in t1/2")
fig.savefig(FIG / "gamow.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: a wave packet hits a barrier (QM-2's split-operator solver)
# Prediction: P_transmitted = int |phi(k)|^2 T(k) dk, the stationary T averaged over the packet.

# %%
Ng, Lg = 8192, 384.0                                  # dx = 0.09375 puts both barrier edges on grid points
x = np.linspace(-Lg, Lg, Ng, endpoint=False)
dx = x[1] - x[0]
kg = 2 * np.pi * np.fft.fftfreq(Ng, dx)
Vb = np.where((x > 0) & (x < 1.5), 1.0, 0.0)
Vb[np.isclose(x, 0) | np.isclose(x, 1.5)] = 0.5         # half-weight edges: the grid barrier has width exactly 1.5
k0, sig, x0 = np.sqrt(1.6), 10.0, -100.0
psi = (2 * np.pi * sig**2)**-0.25 * np.exp(-(x - x0)**2 / (4 * sig**2) + 1j * k0 * x)
dt, nst = 0.05, 3200
eV_, eT_ = np.exp(-0.5j * Vb * dt), np.exp(-0.5j * kg**2 * dt)
for _ in range(nst):
    psi = eV_ * np.fft.ifft(eT_ * np.fft.fft(eV_ * psi))
Ptr = np.sum(np.abs(psi[x > 1.5])**2) * dx
Pref = np.sum(np.abs(psi[x < 0])**2) * dx
kq = np.linspace(k0 - 8 / (2 * sig), k0 + 8 / (2 * sig), 4001)
wk = np.exp(-2 * sig**2 * (kq - k0)**2); wk /= np.trapezoid(wk, kq)
Tavg = np.trapezoid(wk * transfer(kq**2 / 2, [0, 1.5], [0, 1, 0])[0], kq)
assert abs(Ptr - Tavg) < 1e-3 and abs(Ptr + Pref - 1) < 1e-6
print(f"packet: P_trans = {Ptr:.4f}, <T> over |phi(k)|^2 = {Tavg:.4f}, T(k0) = {T_barrier(0.8, 1, 1.5):.4f}")
fig, ax = plt.subplots(figsize=(8, 3))
ax.plot(x, np.abs(psi)**2); ax.fill_between(x, 0, Vb * 0.02, color="0.7")
ax.set_xlim(-250, 250); ax.set_title(f"after t = {dt*nst:.0f}: transmitted {Ptr:.3f}, reflected {Pref:.3f}")
fig.savefig(FIG / "packet_barrier.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 8: numbers used in solutions.md

# %%
E1_SI = np.pi**2 * hb**2 / (2 * me * (1e-9)**2) / eV
lam = C.h * C.c / (3 * E1_SI * eV)
print(f"A1: E1 = {E1_SI:.4f} eV, E2 - E1 = {3*E1_SI:.4f} eV, photon lambda = {lam*1e9:.0f} nm")
print(f"A5: P(|x| > 1/kappa) for delta well = e^-2 = {np.exp(-2):.4f}")
# B2: half-box state c_n = 2(1 - cos(n pi/2))/(n pi)
cB = np.array([quad(lambda x: psi_box(k, x) * np.sqrt(2), 0, 0.5)[0] for k in range(1, 9)])
assert np.allclose(cB, 2 * (1 - np.cos(np.arange(1, 9) * np.pi / 2)) / (np.arange(1, 9) * np.pi))
print(f"B2: c_1..c_4 = {np.round(cB[:4], 4)}, P(E1) = 4/pi^2 = {4/np.pi**2:.4f}, P(E2) = {cB[1]**2:.4f}")
for nmax in (10, 100, 1000):
    kk2 = np.arange(1, nmax + 1)
    print(f"    partial <H> up to n = {nmax}: {np.sum((2*(1-np.cos(kk2*np.pi/2))/(kk2*np.pi))**2 * kk2**2*np.pi**2/2):.1f}")
zB5 = 0.5e-9 * np.sqrt(2 * me * 10 * eV) / hb
print(f"B5: z0 = {zB5:.3f}, bound states = {int(2*zB5/np.pi)+1}")
z1B5 = finite_well_z(zB5)[0][0]
print(f"B5: true E1 = {10*(z1B5**2/zB5**2 - 1):.4f} eV (z1 = {z1B5:.4f}); infinite-well estimate E1 = {E1_SI - 10:.4f} eV")
X4 = np.linalg.matrix_power(X, 4).real
assert np.isclose(X4[0, 0], 0.75)
print(f"B1: <0|x^4|0> = {X4[0,0]:.4f} (hbar/m w)^2 = 3 (hbar/2mw)^2")
print(f"C3: 2(gamma_U - gamma_Po) = {2*(gU-gPo):.2f}, e^that = 10^{2*(gU-gPo)/np.log(10):.1f}")
assert int(2 * zB5 / np.pi) + 1 == 6 and len(finite_well_z(zB5)) == 6


def two_delta(alpha, b):
    """kappa for V = -alpha[delta(x-b)+delta(x+b)]: even k(1+tanh kb) = 2 alpha, odd k(1+coth kb) = 2 alpha."""
    ke_ = brentq(lambda k: k * (1 + np.tanh(k * b)) - 2 * alpha, 1e-9, 2 * alpha + 1)
    ko_ = brentq(lambda k: k * (1 + 1 / np.tanh(k * b)) - 2 * alpha, 1e-9, 2 * alpha + 1) if alpha * b > 0.5 else None
    return ke_, ko_


for b in (1.0, 2.0, 3.0):
    ke_, ko_ = two_delta(1.0, b)
    dE = (ke_**2 - ko_**2) / 2
    print(f"C1: alpha = 1, b = {b}: E_even = {-ke_**2/2:.6f}, E_odd = {-ko_**2/2:.6f}, "
          f"dE = {dE:.3e}, 2 alpha^2 e^(-2 alpha b) = {2*np.exp(-2*b):.3e}")
assert two_delta(1.0, 0.4)[1] is None                 # alpha b < 1/2: odd state gone
# C2: triangular well V = F x (x > 0), wall at 0, F = 1: E_n = -a_n (F^2/2)^(1/3)
an = ai_zeros(3)[0]
Etri = -an * 0.5**(1 / 3)
xt = np.linspace(0, 20, 8002)[1:-1]
Et_fd = fd_eigs(xt, xt, 3)[0]
assert np.allclose(Et_fd, Etri, rtol=1e-5) and abs(airy(an[0])[0]) < 1e-10
print(f"C2: Airy zeros {np.round(an, 5)}, E_n = {np.round(Etri, 4)}, FD {np.round(Et_fd, 4)}")

# %% [markdown]
# ## Your turn
# 1. Half-oscillator: V = x^2/2 for x > 0, infinite wall at 0. Predict the spectrum before
#    computing it (which oscillator states vanish at x = 0?), then use fd_eigs on [0, 10].
# 2. Triangular well: Part 8 did F = 1. Plot the FD eigenfunctions against Ai(x - E_n) (scaled)
#    and check the node theorem. Put in SI numbers for an electron in a field of 10^7 V/m
#    (a MOSFET inversion layer) and find E_1 in meV.
# 3. Two-delta double well: extend the C1 cell to plot ln(dE) vs separation 2b and fit the slope.
#    Then build the same potential with transfer(..., deltas=[...]) and find the bound states
#    as the poles of t(E) at negative E (zeros of M22 at k = i kappa).
# 4. Replace the 12-delta comb with 12 rectangular barriers (Kronig-Penney) and compare the gaps.
# 5. Make the Part 7 packet narrower (sig = 2) and watch P_trans drift from T(k0) toward the average.

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
