# %% [markdown]
# # QM-6 Project: the hydrogen atom, numerically and analytically
# A radial Schrodinger solver (QM-3's finite-difference eigensolver on u(r) = r R(r) with the
# effective potential, eq. (3.3)) checked against the analytic hydrogen solution: energies,
# Laguerre radial functions (conventions mapped), node counts, the power-series recursion
# (3.11) terminating or blowing up, expectation values and the virial theorem, Rydberg
# states, potentials that break the l-degeneracy, real orbitals, and the spectrum in SI.
# Every number in notes §4 and solutions.md is asserted or printed here.
# Units: atomic units (hbar = m_e = e^2/4 pi eps0 = 1; lengths in a0, energies in hartree)
# except in Part 8.

# %%
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
import scipy.constants as C
from math import factorial
from pathlib import Path
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import genlaguerre, sph_harm_y, spherical_jn, hyp1f1, gamma
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.time()
HARTREE_EV = 2 * 13.605693122994

# %% [markdown]
# ## Part 1: the radial finite-difference solver (§3.1-3.2, §3.6)
# On r_j = j h, j = 1..N, with u(0) = u(rmax) = 0 (eq. (3.5)), -u''/2 + V_eff u = E u is a
# symmetric tridiagonal matrix problem: QM-3 project Part 1 with V -> V + l(l+1)/2r^2.

# %%
def radial_solve(V, l, h, rmax, k):
    """Lowest k eigenpairs of -u''/2 + [V(r) + l(l+1)/2r^2] u = E u (eq. 3.3, atomic units).
    Returns energies, u normalized so sum u^2 h = 1 (eq. 3.4), and the grid."""
    r = h * np.arange(1, int(round(rmax / h)))
    d = 1.0 / h**2 + V(r) + l * (l + 1) / (2 * r**2)
    e = -0.5 / h**2 * np.ones(len(r) - 1)
    E, u = eigh_tridiagonal(d, e, select="i", select_range=(0, k - 1))
    return E, u / np.sqrt(h), r


coulomb = lambda r: -1.0 / r
H_GRID = dict(h=0.01, rmax=120.0)
levels = {}
for l in range(4):
    E, u, r = radial_solve(coulomb, l, k=4 - l, **H_GRID)
    for i, Ei in enumerate(E):
        n = l + 1 + i
        levels[(n, l)] = (Ei, u[:, i], r)
        assert abs(Ei + 0.5 / n**2) < 1e-3 * 0.5 / n**2, (n, l, Ei)
for n in range(1, 5):  # l-degeneracy: all l < n give the same E_n
    En = [levels[(n, l)][0] for l in range(n)]
    assert np.ptp(En) < 1e-4 * 0.5 / n**2
print("Part 1: E_nl = -1/2n^2 for all n<=4, l<n (rel. err < 1e-3); l-degeneracy holds")
for n in range(1, 5):
    print("   n=%d: " % n + "  ".join(f"l={l}: {levels[(n, l)][0]:.7f}" for l in range(n)))
err = [abs(radial_solve(coulomb, 0, hh, 60.0, 1)[0][0] + 0.5) / 0.5 for hh in (0.01, 0.005)]
print(f"   1s relative error: h=0.01 -> {err[0]:.2e}, h=0.005 -> {err[1]:.2e} (ratio {err[0]/err[1]:.2f})")
assert 3.5 < err[0] / err[1] < 4.5
# Example 4.3: the n = 3 states
n3 = [(l, m) for l in range(3) for m in range(-l, l + 1)]
assert len(n3) == 9 and all(abs(levels[(3, l)][0] + 1 / 18) < 1e-6 for l in range(3))
print(f"Ex 4.3: n=3 has {len(n3)} states; E_3 = {-13.605693/9:.3f} eV; FD: "
      + ", ".join(f"{levels[(3, l)][0]:.7f}" for l in range(3)))
# A1: n = 4
n4 = [(l, m) for l in range(4) for m in range(-l, l + 1)]
print(f"A1: n=4 has {len(n4)} states, {2*len(n4)} with spin")
assert len(n4) == 16

# %% [markdown]
# ## Part 2: analytic radial functions, Laguerre conventions, nodes (§3.7)
# Modern L_k^(alpha) (scipy genlaguerre, sympy assoc_laguerre) vs the old
# L^p_{q-p} = (-1)^p d^p/dx^p [e^x d^q/dx^q (e^-x x^q)] = q! L^(p)_{q-p}.

# %%
def R_nl(n, l, r):
    """Hydrogen radial function, eq. (3.15), modern Laguerre convention, a0 = 1."""
    x = 2 * r / n
    norm = np.sqrt((2 / n)**3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-r / n) * x**l * genlaguerre(n - l - 1, 2 * l + 1)(x)


xs = sp.symbols("x", positive=True)
for q in range(1, 6):
    Lq = sp.exp(xs) * sp.diff(sp.exp(-xs) * xs**q, xs, q)
    for p in range(0, q + 1):
        old = (-1)**p * sp.diff(Lq, xs, p)
        assert sp.simplify(old - factorial(q) * sp.assoc_laguerre(q - p, p, xs)) == 0
print("Part 2: old-convention L^p_{q-p} = q! x modern L^(p)_{q-p} for q <= 5 (SymPy)")
assert np.allclose(genlaguerre(1, 1)(np.array([0.0, 1.0, 3.0])), [2.0, 1.0, -1.0])

rr = np.linspace(0.01, 30, 500)
table = {(1, 0): 2 * np.exp(-rr), (2, 0): (1 - rr / 2) * np.exp(-rr / 2) / np.sqrt(2),
         (2, 1): rr * np.exp(-rr / 2) / np.sqrt(24),
         (3, 0): 2 / np.sqrt(27) * (1 - 2 * rr / 3 + 2 * rr**2 / 27) * np.exp(-rr / 3)}
for key, ref in table.items():  # eq. (3.16)
    assert np.allclose(R_nl(*key, rr), ref, atol=1e-13)
print(f"   (3.16) checked; R_20(0) = {R_nl(2, 0, 0.0):.6f} = 1/sqrt2 = {1/np.sqrt(2):.6f}")

def count_nodes(u):
    big = np.abs(u) > 1e-6 * np.abs(u).max()
    s = np.sign(u[big])
    return int(np.sum(s[1:] != s[:-1]))


fig, axes = plt.subplots(2, 2, figsize=(10, 7))
for (n, l), (E, u, r) in sorted(levels.items()):
    ua = r * R_nl(n, l, r)
    sgn = np.sign(np.dot(u, ua))
    assert np.max(np.abs(sgn * u - ua)) < 2e-3 * np.abs(ua).max(), (n, l)
    assert count_nodes(u) == n - l - 1
    assert abs(quad(lambda x: (x * R_nl(n, l, x))**2, 0, np.inf)[0] - 1) < 1e-10
    ax = axes.flat[n - 1]
    ax.plot(r, ua**2, label=f"{n}{'spdf'[l]} ({n-l-1} nodes)")
    ax.plot(r[::150], u[::150]**2, "k.", ms=3)
for n, ax in enumerate(axes.flat, 1):
    ax.set_xlim(0, 4 * n**2 + 5); ax.set_title(f"n = {n}: r^2|R|^2 (lines analytic, dots FD)")
    ax.set_xlabel("r / a0"); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "radial_densities.png", dpi=110); plt.close(fig)
print("   FD u_nl matches r R_nl; nodes = n - l - 1; int u^2 = 1 for all n <= 4")
# A2: R_21 normalized, peak of r^2 R_21^2 at 4 a0
assert abs(quad(lambda x: x**2 * R_nl(2, 1, x)**2, 0, np.inf)[0] - 1) < 1e-12
rf = np.linspace(0.01, 30, 30000)
rpk = rf[np.argmax(rf**2 * R_nl(2, 1, rf)**2)]
print(f"A2: int r^2 R_21^2 = 1; peak of P_21 at r = {rpk:.2f} a0 (exact 4)")
assert abs(rpk - 4) < 2e-3

# %% [markdown]
# ## Part 3: the series, its termination, and its blow-up (§3.5, Examples 4.1-4.2)

# %%
rho, rho0, l_s = sp.symbols("rho rho0 l", positive=True)
v = sp.Function("v")
u_s = rho**(l_s + 1) * sp.exp(-rho) * v(rho)
lhs = sp.diff(u_s, rho, 2) - (1 - rho0 / rho + l_s * (l_s + 1) / rho**2) * u_s
veq = rho * v(rho).diff(rho, 2) + 2 * (l_s + 1 - rho) * v(rho).diff(rho) + (rho0 - 2 * (l_s + 1)) * v(rho)
assert sp.simplify(sp.expand(lhs / (rho**l_s * sp.exp(-rho))) - veq) == 0
j = sp.symbols("j", integer=True, nonnegative=True)
c = sp.symbols("c0:12")
for lv in range(3):  # recursion (3.11) applied to a truncated series: coefficient of rho^j vanishes
    V = sum(c[k] * rho**k for k in range(12))
    expr = sp.expand(veq.subs(l_s, lv).subs(v(rho), V).doit())
    for jj in range(10):
        coeff = expr.coeff(rho, jj)
        assert sp.simplify(coeff - ((jj + 1) * (jj + 2 * lv + 2) * c[jj + 1]
                                    - (2 * (jj + lv + 1) - rho0) * c[jj])) == 0
print("Part 3: (3.10) and recursion (3.11) verified symbolically")

def series_coeffs(rho_0, l, J):
    """c_0..c_J from (3.11) with c_0 = 1."""
    cs = [1.0]
    for jj in range(J):
        cs.append(cs[-1] * (2 * (jj + l + 1) - rho_0) / ((jj + 1) * (jj + 2 * l + 2)))
    return np.array(cs)


c4 = series_coeffs(4.0, 0, 10)
assert c4[1] == -1.0 and np.all(c4[2:] == 0)
c6 = series_coeffs(6.0, 0, 10)  # B1: n = 3, l = 0
assert np.allclose(c6[:3], [1, -2, 2 / 3]) and np.all(c6[3:] == 0)
print(f"Ex 4.2: rho0=4, l=0: c = {c4[:3]} -> v = 1 - rho;  B1: rho0=6: c = {np.round(c6[:4], 6)}")
# Example 4.1: E = -5 eV
E5 = -5 / HARTREE_EV
kap = np.sqrt(-2 * E5); r0_5 = 2 / kap
c5 = series_coeffs(r0_5, 0, 6)
print(f"Ex 4.1: E=-5 eV = {E5:.5f} Eh, kappa = {kap:.5f}, rho0 = {r0_5:.4f}, c1 = {c5[1]:.4f}, c2 = {c5[2]:.4f}")
assert abs(r0_5 - 3.2992) < 1e-4 and abs(c5[1] + 0.6496) < 1e-4 and abs(c5[2] + 0.0759) < 1e-4
assert np.all(c5[3:] / c5[2:-1] > 0)
# non-terminating rho0 = 4.01: v = 1F1(a; 2; 2 rho), a = 1 - rho0/2; compare series, hyp1f1, asymptotics
cbad = series_coeffs(4.01, 0, 400)
rg = np.linspace(0, 30, 301)
vbad = np.polynomial.polynomial.polyval(rg, cbad)
a_k = 1 - 4.01 / 2
assert np.allclose(vbad, hyp1f1(a_k, 2, 2 * rg), rtol=1e-8)
asym = gamma(2) / gamma(a_k) * np.exp(2 * rg[1:]) * (2 * rg[1:])**(a_k - 2)
ratio_end = vbad[-1] / asym[-1]
print(f"   rho0=4.01: series = 1F1; v(30)/[Gamma(2)/Gamma(a) e^(2rho)(2rho)^(a-2)] = {ratio_end:.3f}")
assert abs(ratio_end - 1) < 0.15
fig, ax = plt.subplots(figsize=(7, 4))
ax.semilogy(rg, np.abs(rg * np.exp(-rg) * np.polynomial.polynomial.polyval(rg, c4)) + 1e-30, label="rho0 = 4 (terminates: R_20)")
ax.semilogy(rg, np.abs(rg * np.exp(-rg) * vbad) + 1e-30, label="rho0 = 4.01 (grows like e^rho)")
ax.set_xlabel("rho"); ax.set_ylabel("|u(rho)|"); ax.legend(); ax.set_title("Theorem 3.1: termination or blow-up")
fig.tight_layout(); fig.savefig(FIG / "series_blowup.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 4: expectation values, the virial theorem, inside the nucleus (§3.8, Ex. 4.4, 4.6)

# %%
def expval(f, n, l):
    return quad(lambda x: f(x) * (x * R_nl(n, l, x))**2, 0, np.inf, limit=200)[0]


r1, r2, rinv = expval(lambda x: x, 1, 0), expval(lambda x: x**2, 1, 0), expval(lambda x: 1 / x, 1, 0)
u10p = lambda x: 2 * (1 - x) * np.exp(-x)
T1 = 0.5 * quad(lambda x: u10p(x)**2, 0, np.inf)[0]
V1 = -rinv
print(f"Ex 4.4: <r> = {r1:.6f}, <r^2> = {r2:.6f}, sigma_r = {np.sqrt(r2 - r1**2):.4f}, <1/r> = {rinv:.6f} "
      f"(a0 = {C.physical_constants['Bohr radius'][0]*1e9:.4f} nm, <r> = {1.5*0.0529177:.4f} nm)")
print(f"   <T> = {T1:.6f} Eh = {T1*HARTREE_EV:.3f} eV, <V> = {V1:.6f} Eh = {V1*HARTREE_EV:.3f} eV")
assert np.isclose(r1, 1.5) and np.isclose(r2, 3.0) and np.isclose(rinv, 1.0)
assert abs(2 * T1 + V1) < 1e-10 and np.isclose(V1, 2 * (-0.5))
E1fd, u1fd, rfd = levels[(1, 0)]
assert abs(np.sum(rfd * u1fd**2) * 0.01 - 1.5) < 1e-3   # same <r> from the FD eigenvector
for n in range(1, 5):            # eq. (3.17) and C2's <1/r^2>
    for l in range(n):
        assert np.isclose(expval(lambda x: x, n, l), 0.5 * (3 * n**2 - l * (l + 1)))
        assert np.isclose(expval(lambda x: x**2, n, l), 0.5 * n**2 * (5 * n**2 + 1 - 3 * l * (l + 1)))
        assert np.isclose(expval(lambda x: 1 / x, n, l), 1 / n**2)
        if l > 0:
            assert np.isclose(expval(lambda x: x**-2, n, l), 1 / (n**3 * (l + 0.5)))
print("   (3.17) <r>, <r^2>, <1/r> and C2's <1/r^2> = 1/(n^3 (l+1/2)) verified for all n <= 4")
r21, r22 = expval(lambda x: x, 2, 1), expval(lambda x: x**2, 2, 1)
print(f"B2: 2p: <r> = {r21:.4f}, <r^2> = {r22:.4f}, sigma_r = {np.sqrt(r22 - r21**2):.4f}")
a0 = C.physical_constants["Bohr radius"][0]
psi0sq = 1 / (np.pi * a0**3)
print(f"A5: |psi_100(0)|^2 = 1/(pi a0^3) = {psi0sq:.3e} m^-3")

def p_inside(x):
    """P(r < x a) for 1s: 1 - e^{-2x}(1 + 2x + 2x^2), by its series 4/3 x^3 - 2 x^4 + 4/3 x^5."""
    return 4 / 3 * x**3 - 2 * x**4 + 4 / 3 * x**5


xb = 1e-15 / a0
P_series = p_inside(xb)
P_quad = quad(lambda s: 4 * s**2 * np.exp(-2 * s), 0, xb, epsabs=0, epsrel=1e-12)[0]
P_naive = 1 - np.exp(-2 * xb) * (1 + 2 * xb + 2 * xb**2)
P_expm1 = -np.expm1(-2 * xb) - np.exp(-2 * xb) * (2 * xb + 2 * xb**2)
print(f"Ex 4.6: x = {xb:.4e}; P series = {P_series:.4e}, quad = {P_quad:.4e}, naive closed form = {P_naive:.3e}")
assert np.isclose(P_series, P_quad, rtol=1e-8) and abs(P_series - 9.00e-15) < 0.01e-15
assert abs(P_naive / P_series - 1) > 0.005  # the cancellation is visible

# %% [markdown]
# ## Part 5: the circular Rydberg state n = 10 (§3.15, Example 4.7)

# %%
E9, u9, r9 = radial_solve(coulomb, 9, 0.05, 400.0, 1)
rpk9 = r9[np.argmax(u9[:, 0]**2)]
rmean9 = np.sum(r9 * u9[:, 0]**2) * 0.05
print(f"Ex 4.7: l=9 lowest E = {E9[0]:.6f} (exact {-1/200}), peak of |u|^2 at {rpk9:.2f} a0 = "
      f"{rpk9*a0*1e9:.2f} nm, <r> = {rmean9:.2f}, sigma/<r> = 1/sqrt(21) = {1/np.sqrt(21):.3f}")
assert abs(E9[0] + 0.005) < 1e-3 * 0.005 and abs(rpk9 - 100) < 0.5 and abs(rmean9 - 105) < 0.1
for n in (1, 10, 100):  # B7: circular-state relative width
    print(f"B7: n = {n}: sigma_r/<r> = 1/sqrt(2n+1) = {1/np.sqrt(2*n+1):.4f}")

# %% [markdown]
# ## Part 6: potentials that break the l-degeneracy (§3.12, §3.3)

# %%
lam = 5.0
yuk = lambda r: -np.exp(-r / lam) / r
Ey0 = radial_solve(yuk, 0, 0.01, 150, 2)[0]
Ey1 = radial_solve(yuk, 1, 0.01, 150, 1)[0]
print(f"Part 6: Yukawa lambda=5: E_1s = {Ey0[0]:.5f}, E_2s = {Ey0[1]:.5f}, E_2p = {Ey1[0]:.5f} (2s below 2p)")
assert Ey0[1] < Ey1[0] < 0
osc = {}
for l in range(4):
    Eo = radial_solve(lambda r: 0.5 * r**2, l, 0.005, 12, 3)[0]
    assert np.allclose(Eo, 2 * np.arange(3) + l + 1.5, atol=1e-3)
    osc[l] = Eo
deg = {N: sum(2 * l + 1 for l in range(N + 1) if (N - l) % 2 == 0) for N in range(5)}
assert [deg[N] for N in range(5)] == [(N + 1) * (N + 2) // 2 for N in range(5)]
print(f"   3D oscillator: E = 2n_r + l + 3/2 (FD), degeneracies by N: {[deg[N] for N in range(5)]}")
def jl_zeros(l, k):
    """First k positive zeros of the spherical Bessel function j_l (eq. 3.6)."""
    xs_, out = np.linspace(0.5, 30, 6000), []
    f = spherical_jn(l, xs_)
    for i in np.nonzero(np.sign(f[:-1]) != np.sign(f[1:]))[0][:k]:
        out.append(brentq(lambda x: spherical_jn(l, x), xs_[i], xs_[i + 1]))
    return np.array(out)


beta = {l: jl_zeros(l, 3) for l in range(4)}
order = sorted(((beta[l][k], f"{k+1}{'spdf'[l]}") for l in range(4) for k in range(3)))[:6]
print("   spherical well: " + ", ".join(f"{name} {b:.4f}" for b, name in order))
assert np.isclose(beta[1][0], 4.4934, atol=1e-4) and np.isclose(beta[0][1], 2 * np.pi)
assert abs(np.tan(beta[1][0]) - beta[1][0]) < 1e-8
E_well = C.hbar**2 * np.pi**2 / (2 * C.m_e * (1e-9)**2) / C.e
print(f"B4: E_11/E_10 = (beta_11/pi)^2 = {(beta[1][0]/np.pi)**2:.4f}; electron, a = 1 nm: E_10 = {E_well:.4f} eV, "
      f"E_11 = {E_well*(beta[1][0]/np.pi)**2:.4f} eV, E_20 = {4*E_well:.4f} eV")
fig, axs = plt.subplots(1, 3, figsize=(11, 4))
for ax, (title, data) in zip(axs, [("Coulomb", {l: [levels[(n, l)][0] for n in range(l + 1, 5)] for l in range(4)}),
                                   ("Yukawa, lambda = 5", {0: Ey0, 1: Ey1}), ("3D oscillator", osc)]):
    for l, Es in data.items():
        for Ei in Es:
            ax.hlines(Ei, l - 0.3, l + 0.3)
    ax.set_xticks(range(4)); ax.set_xticklabels(list("spdf")); ax.set_title(title)
axs[0].set_ylabel("E (hartree)")
fig.tight_layout(); fig.savefig(FIG / "degeneracy_lifted.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 7: real orbitals and cross-sections (§3.9)

# %%
th, ph = rng.uniform(0, np.pi, 50), rng.uniform(0, 2 * np.pi, 50)
Y = lambda l, m: sph_harm_y(l, m, th, ph)
sx, sy, cz = np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)
px = (Y(1, -1) - Y(1, 1)) / np.sqrt(2)
py = 1j * (Y(1, -1) + Y(1, 1)) / np.sqrt(2)
dxy = (Y(2, 2) - Y(2, -2)) / (1j * np.sqrt(2))
dxz = (Y(2, -1) - Y(2, 1)) / np.sqrt(2)  # A6
assert np.allclose(px, np.sqrt(3 / (4 * np.pi)) * sx) and np.allclose(py, np.sqrt(3 / (4 * np.pi)) * sy)
assert np.allclose(dxy, np.sqrt(15 / (4 * np.pi)) * sx * sy)
assert np.allclose(dxz, np.sqrt(15 / (4 * np.pi)) * sx * cz)
assert np.allclose(Y(1, 1), -np.sqrt(3 / (8 * np.pi)) * np.sin(th) * np.exp(1j * ph))
print("Part 7: p_x, p_y (3.19), d_xy, and A6's d_xz = (Y2^-1 - Y2^1)/sqrt2 checked pointwise")

def psi_cart(n, l, angular, x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2) + 1e-12
    t = np.arccos(np.clip(z / r, -1, 1)); p = np.arctan2(y, x)
    return R_nl(n, l, r) * angular(t, p)


g = np.linspace(-1, 1, 241)
cases = [("1s (xz)", 1, 0, lambda t, p: sph_harm_y(0, 0, t, p).real, 5, "xz"),
         ("2p0 (xz)", 2, 1, lambda t, p: sph_harm_y(1, 0, t, p).real, 12, "xz"),
         ("3d0 (xz)", 3, 2, lambda t, p: sph_harm_y(2, 0, t, p).real, 25, "xz"),
         ("2p_x (xz)", 2, 1, lambda t, p: np.sqrt(3 / (4 * np.pi)) * np.sin(t) * np.cos(p), 12, "xz"),
         ("3d_xy (xy)", 3, 2, lambda t, p: np.sqrt(15 / (4 * np.pi)) * np.sin(t)**2 * np.sin(p) * np.cos(p), 25, "xy")]
fig, axs = plt.subplots(1, 5, figsize=(16, 3.4))
for ax, (title, n, l, ang, L, plane) in zip(axs, cases):
    A, B = np.meshgrid(g * L, g * L)
    dens = np.abs(psi_cart(n, l, ang, A, B if plane == "xy" else 0 * A, B if plane == "xz" else 0 * A))**2
    ax.imshow(dens, extent=[-L, L, -L, L], origin="lower", cmap="magma")
    ax.set_title(title); ax.set_xlabel("x / a0")
fig.tight_layout(); fig.savefig(FIG / "orbitals.png", dpi=100); plt.close(fig)

# %% [markdown]
# ## Part 8: SI numbers, the spectrum, hydrogen-like atoms (§3.10-3.14, Examples 4.5, 4.8)

# %%
me, e, hbar, eps0, c, h = C.m_e, C.e, C.hbar, C.epsilon_0, C.c, C.h
mp = C.m_p; md = C.physical_constants["deuteron mass"][0]
ma = C.physical_constants["alpha particle mass"][0]; mmu = C.physical_constants["muon mass"][0]
a0_calc = 4 * np.pi * eps0 * hbar**2 / (me * e**2)
E1 = -me / (2 * hbar**2) * (e**2 / (4 * np.pi * eps0))**2 / e
Rinf = me * e**4 / (8 * eps0**2 * h**3 * c)
RH, RD, RHe = (Rinf / (1 + me / M) for M in (mp, md, ma))
assert np.isclose(a0_calc, a0, rtol=1e-8) and np.isclose(Rinf, C.Rydberg, rtol=1e-8)
print(f"Part 8: a0 = {a0_calc:.5e} m, E1 = {E1:.4f} eV (mu: {E1/(1+me/mp):.4f}), alpha = {C.alpha:.9f} = 1/{1/C.alpha:.3f}")
print(f"   R_inf = {Rinf:.5e} m^-1, R_H = {RH:.5e} m^-1")
lam = lambda R, nf, ni: 1e9 / (R * (1 / nf**2 - 1 / ni**2))
Ha, Ha_inf = lam(RH, 2, 3), lam(Rinf, 2, 3)
print(f"Ex 4.5: H-alpha = {Ha:.2f} nm vac ({Ha/1.000277:.2f} air; R_inf gives {Ha_inf:.2f}), photon {-E1/(1+me/mp)*(1/4-1/9):.3f} eV;"
      f" Lyman limit = {1e9/RH:.2f} nm, Lyman-alpha = {lam(RH,1,2):.2f} nm")
assert abs(Ha - 656.47) < 0.01 and abs(1e9 / RH - 91.18) < 0.01 and abs(lam(RH, 1, 2) - 121.57) < 0.01
print(f"A3: Paschen-alpha = {lam(RH,3,4):.1f} nm, H-beta = {lam(RH,2,4):.2f} nm, H-gamma = {lam(RH,2,5):.2f}, "
      f"H-delta = {lam(RH,2,6):.2f}, Balmer limit = {4e9/RH:.1f} nm")
series = {"Lyman": (1, "tab:purple"), "Balmer": (2, "tab:red"), "Paschen": (3, "tab:brown")}
fig, ax = plt.subplots(figsize=(10, 3))
for name, (nf, col) in series.items():
    wl = [lam(RH, nf, ni) for ni in range(nf + 1, nf + 15)]
    ax.vlines(wl, 0, 1, color=col, lw=1, label=name)
ax.set_xscale("log"); ax.set_xlabel("vacuum wavelength (nm)"); ax.set_yticks([]); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "spectrum.png", dpi=110); plt.close(fig)

# He+ from the solver with Z = 2 (eq. 3.21) and Pickering's line
for l in range(2):
    EZ = radial_solve(lambda r: -2.0 / r, l, 0.005, 60, 3 - l)[0]
    assert np.allclose(EZ, -2.0 / np.arange(l + 1, 4)**2, rtol=1e-3)
He43 = 1e9 / (4 * RHe * (1 / 9 - 1 / 16))
print(f"Ex 4.8 / A4: He+ FD E_n = -2/n^2 Eh; E1(He+) = {4*E1/(1+me/ma):.3f} eV, a = a0/2, "
      f"4->3 line = {He43:.2f} nm, He+ Lyman-alpha = {1e9/(4*RHe*0.75):.2f} nm, R_He = {RHe:.5e}")
assert abs(He43 - 468.72) < 0.01
print(f"   Li2+ E1 = {9*E1:.1f} eV")
# muonic hydrogen (B5) and the proton
mu_mu = mmu * mp / (mmu + mp) / me
amu = a0 / mu_mu
rp = 0.84e-15
PH = p_inside(rp / (a0 * (1 + me / mp)))
Pmu = p_inside(rp / amu)
print(f"B5: m_mu/m_e = {mmu/me:.3f}, mu = {mu_mu:.3f} m_e, E1 = {E1*mu_mu/1e3:.3f} keV, a = {amu*1e15:.1f} fm "
      f"= {amu/rp:.0f} r_p; P(r<r_p): H {PH:.2e}, muH {Pmu:.2e}, ratio {Pmu/PH:.3e}")
# deuterium isotope shift (B6)
dl = lam(RH, 2, 3) - lam(RD, 2, 3)
print(f"B6: R_D = {RD:.5e}, H-alpha(D) = {lam(RD,2,3):.3f} nm, shift H - D = {dl:.4f} nm")
assert abs(dl - 0.1786) < 0.0005
# §3.11 scales
fs = abs(E1) * C.alpha**2 / 16
print(f"   fine structure 2p3/2-2p1/2 = {fs:.3e} eV = {fs*e/h/1e9:.2f} GHz; Lamb 1057.8 MHz = {h*1057.8e6/e:.2e} eV; "
      f"21 cm: {c/1420.405751e6*100:.2f} cm, {h*1420.405751e6/e:.2e} eV; reduced-mass shift of E1 = {abs(E1)*(me/mp)/(1+me/mp)*1e3:.1f} meV")

# C3: potassium, V = -[1 + (Z-1) e^{-r/d}]/r, Z = 19: 4s (4th s state) vs 3d (1st d state)
def k_levels(d, h=0.002):
    V = lambda r: -(1 + 18 * np.exp(-r / d)) / r
    return radial_solve(V, 0, h, 80, 4)[0][3], radial_solve(V, 2, h, 80, 1)[0][0]


E4s, E3d = k_levels(0.5)
print(f"C3: d = 0.5 a0: E_4s = {E4s*HARTREE_EV:.2f} eV, E_3d = {E3d*HARTREE_EV:.2f} eV (measured K: -4.34, -1.67)")
assert E4s < E3d
lo, hi = 0.5, 0.6
assert k_levels(hi)[1] < k_levels(hi)[0]
for _ in range(12):
    mid = 0.5 * (lo + hi)
    a4, a3 = k_levels(mid, 0.004)
    lo, hi = (mid, hi) if a4 < a3 else (lo, mid)
print(f"   ordering flips (3d below 4s) at d = {0.5*(lo+hi):.3f} a0")
E4s_h, E3d_h = -0.5 / 16, -0.5 / 9
print(f"   unscreened hydrogenic: E_4s = {E4s_h*HARTREE_EV:.2f} eV > E_3d = {E3d_h*HARTREE_EV:.2f} eV")

# %% [markdown]
# ## Your turn
# 1. Infinite spherical well: call radial_solve with V = 0 on rmax = 1 (the wall is the grid's
#    end) for l = 0..3, and compare 2E with the squared zeros beta_{nl}^2 from jl_zeros (eq. 3.6).
# 2. Positronium: the reduced mass is m_e/2. Rescale the Coulomb solver (or the formula) and show
#    E_n = -6.80 eV/n^2 and the Lyman-alpha line at 243 nm.
# 3. Z_eff model: with Slater's rules Z_eff(4s, K) = 2.20 and Z_eff(3d, K) = 1.00, and effective
#    principal numbers n* = 3.7 and 3, compare -13.6 Z_eff^2/n*^2 for 4s and 3d. Which is lower?
#    Then redo C3's screened potential for calcium (Z = 20, two valence electrons ignored).

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T0:.1f} s.")
