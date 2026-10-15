# %% [markdown]
# # EM-3 Project: Boundary-Value Problems on a Computer
#
# This script builds the tools of EM-3 and checks the notes against them:
# Legendre polynomials three ways (Bonnet recursion, Rodrigues, Frobenius
# series) with orthogonality and the generating function; the non-terminating
# Frobenius series diverging at x = 1; images for a plane and a sphere (field
# lines, induced charge, force, energy); the slot problem by Fourier series and
# by a red-black SOR relaxation solver; the conducting sphere in a uniform field
# by relaxation on an (ln r, theta) grid; multipole error curves and the
# quadrupole tensor; spherical harmonics; dipoles in external fields. Every
# number quoted in the notes' worked examples and in solutions.md is asserted
# or printed here.
#
# SI units. Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
from math import factorial
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sympy as sp
from scipy import integrate
from scipy.special import roots_legendre, eval_legendre
try:
    from scipy.special import sph_harm_y                        # SciPy >= 1.15: (l, m, polar, azimuth)
except ImportError:                                             # older SciPy: sph_harm(m, l, azimuth, polar)
    from scipy.special import sph_harm
    sph_harm_y = lambda l, m, th, ph: sph_harm(m, l, ph, th)

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.perf_counter()
EPS0 = 8.8541878128e-12             # F/m
K = 1 / (4 * np.pi * EPS0)          # Coulomb's constant
nC = 1e-9

# %% [markdown]
# ## Part 0: Legendre polynomials and the closed-form solutions, symbolically (§3.5–3.7, §4)
# Bonnet (3.18), Rodrigues (3.17), the Frobenius recursion (3.14), Legendre's
# equation (3.13), orthogonality (Theorem 3.5), the generating function (3.19)
# and the series behind the normalization. Then Examples 4.4, 4.5, (3.25), A1, A4.

# %%
x, t = sp.symbols("x t")


def legendre_recurrence(lmax):
    """P_0..P_lmax from Bonnet's recursion (3.18): (l+1) P_{l+1} = (2l+1) x P_l - l P_{l-1}."""
    P = [sp.Integer(1), x]
    for l in range(1, lmax):
        P.append(sp.expand(((2 * l + 1) * x * P[l] - l * P[l - 1]) / (l + 1)))
    return P[: lmax + 1]


def legendre_rodrigues(l):
    """Rodrigues' formula (3.17): P_l = (1 / 2^l l!) d^l/dx^l (x^2 - 1)^l."""
    return sp.expand(sp.diff((x**2 - 1) ** l, x, l) / (2**l * factorial(l)))


def frobenius_series(l, parity, kmax):
    """Series sum_k a_k x^k from the recursion (3.14), a_parity = 1, truncated at x^kmax."""
    a, y = sp.Integer(1), sp.Integer(0)
    for k in range(parity, kmax + 1, 2):
        y += a * x**k
        a = a * sp.Rational((k - l) * (k + l + 1), (k + 1) * (k + 2))
    return y


LMAX = 10
P = legendre_recurrence(LMAX)
for l in range(LMAX + 1):
    fr = frobenius_series(l, l % 2, l + 4)                     # terminates: terms beyond x^l are zero
    assert sp.expand(P[l] - legendre_rodrigues(l)) == 0
    assert sp.expand(P[l] - fr / fr.subs(x, 1)) == 0
    assert sp.expand(P[l] - sp.legendre(l, x)) == 0
    assert sp.expand(sp.diff((1 - x**2) * sp.diff(P[l], x), x) + l * (l + 1) * P[l]) == 0   # (3.13)
    assert P[l].subs(x, 1) == 1 and sp.expand(P[l].subs(x, -x) - (-1) ** l * P[l]) == 0
for l in range(7):
    for m in range(7):
        assert sp.integrate(P[l] * P[m], (x, -1, 1)) == (sp.Rational(2, 2 * l + 1) if l == m else 0)
gen = sp.series(1 / sp.sqrt(1 - 2 * x * t + t**2), t, 0, LMAX + 1).removeO()
assert all(sp.expand(gen.coeff(t, l) - P[l]) == 0 for l in range(LMAX + 1))       # (3.19), (3.28)
G_ = 1 / sp.sqrt(1 - 2 * x * t + t**2)
assert sp.simplify((1 - 2 * x * t + t**2) * sp.diff(G_, t) - (x - t) * G_) == 0                     # C2(a)
assert sp.simplify(sp.diff((1 - x**2) * sp.diff(G_, x), x) + t * sp.diff(t * G_, t, 2)) == 0        # C2(b)
nser = sp.series(sp.log((1 + t) / (1 - t)) / t, t, 0, 12).removeO()               # Theorem 3.5, step 4
assert all(nser.coeff(t, 2 * l) == sp.Rational(2, 2 * l + 1) for l in range(6))
odd1 = frobenius_series(1, 0, 14)                                                    # B4(b): l = 1, even series
assert sp.expand(odd1 - sp.series(1 - x * sp.atanh(x), x, 0, 15).removeO()) == 0
print("P_2 by Rodrigues:", legendre_rodrigues(2), "| P_3 by Frobenius (B4a):", P[3], "| P_4 =", P[4])
print("l = 1 even Frobenius series (B4b):", odd1, "= 1 - x artanh(x) + ...")

r, th, ph, R_, E0_, s0_, V0_, ep, s_, a_ = sp.symbols("r theta phi R E_0 sigma_0 V_0 epsilon_0 s a", positive=True)
lap_sph = lambda f: sp.diff(r**2 * sp.diff(f, r), r) / r**2 + sp.diff(sp.sin(th) * sp.diff(f, th), th) / (r**2 * sp.sin(th))
V44 = -E0_ * (r - R_**3 / r**2) * sp.cos(th)                                        # (4.2)
sig44 = -ep * sp.diff(V44, r).subs(r, R_)
assert sp.simplify(lap_sph(V44)) == 0 and V44.subs(r, R_) == 0 and sp.simplify(sig44 - 3 * ep * E0_ * sp.cos(th)) == 0
assert sp.simplify(sp.integrate(R_ * sp.cos(th) * sig44 * 2 * sp.pi * R_**2 * sp.sin(th), (th, 0, sp.pi)) - 4 * sp.pi * ep * R_**3 * E0_) == 0
Vin, Vout = s0_ / (3 * ep) * r * sp.cos(th), s0_ * R_**3 / (3 * ep) * sp.cos(th) / r**2    # (4.5)
assert sp.simplify(lap_sph(Vin)) == 0 and sp.simplify(lap_sph(Vout)) == 0
assert sp.simplify((Vin - Vout).subs(r, R_)) == 0
assert sp.simplify((sp.diff(Vout, r) - sp.diff(Vin, r)).subs(r, R_) + s0_ * sp.cos(th) / ep) == 0   # EM-2 (3.12)
Vcyl = -E0_ * (s_ - a_**2 / s_) * sp.cos(ph)                                        # (3.25)
assert sp.simplify(sp.diff(s_ * sp.diff(Vcyl, s_), s_) / s_ + sp.diff(Vcyl, ph, 2) / s_**2) == 0
assert sp.simplify(-ep * sp.diff(Vcyl, s_).subs(s_, a_) - 2 * ep * E0_ * sp.cos(ph)) == 0
X0, Y0, Z0, Rs, u_, v_ = sp.symbols("x_0 y_0 z_0 R_s u v", real=True)               # A1: sphere averages
pt = (X0 + Rs * sp.sin(u_) * sp.cos(v_), Y0 + Rs * sp.sin(u_) * sp.sin(v_), Z0 + Rs * sp.cos(u_))
avg = lambda f: sp.simplify(sp.integrate(sp.integrate(f(*pt) * sp.sin(u_), (u_, 0, sp.pi)), (v_, 0, 2 * sp.pi)) / (4 * sp.pi))
assert sp.simplify(avg(lambda X, Y, Z: X**2 + Y**2 - 2 * Z**2) - (X0**2 + Y0**2 - 2 * Z0**2)) == 0
assert sp.simplify(avg(lambda X, Y, Z: X**2 + Y**2 + Z**2) - (X0**2 + Y0**2 + Z0**2 + Rs**2)) == 0
A_l = [sp.Rational(2 * l + 1, 2) * sp.integrate(x**2 * P[l], (x, -1, 1)) for l in range(4)]   # A4: cos^2 = P0/3 + 2P2/3
assert A_l == [sp.Rational(1, 3), 0, sp.Rational(2, 3), 0]
print("Symbolic checks passed: (4.2)-(4.5), (3.25), A1 sphere averages, A4 coefficients", A_l)

# %% [markdown]
# ## Part 1: Why l must be an integer (§3.5, Frobenius derivation steps 3-5)
# For non-integer l neither series of (3.14) terminates; k a_k -> C != 0, so the
# partial sums at x = 1 grow like (C/2) ln N. For integer l one series stops.


# %%
def frob_coeffs(l, parity, kmax):
    """Indices k and coefficients a_k of the series solution (3.14) with a_parity = 1 (floating point)."""
    k = np.arange(parity, kmax + 1, 2, dtype=float)
    ratio = (k[:-1] - l) * (k[:-1] + l + 1) / ((k[:-1] + 1) * (k[:-1] + 2))
    return k, np.concatenate([[1.0], np.cumprod(ratio)])


fig, axs = plt.subplots(1, 2, figsize=(11, 4))
for l in (2.0, 2.5, 3.5):
    k, a = frob_coeffs(l, 0, 200_000)
    S = np.cumsum(a)                                   # partial sums at x = 1
    axs[0].semilogx(k[1:], S[1:], label=f"l = {l}, even series")
    if l == 2.0:
        assert np.all(a[2:] == 0) and S[-1] == -2.0    # 1 - 3x^2: terminates, = -2 P_2
    else:
        C = k[-1] * a[-1]
        sel = k > 1000
        slope = np.polyfit(np.log(k[sel]), S[sel], 1)[0]
        assert abs(slope / (C / 2) - 1) < 2e-3
        print(f"l = {l}: k a_k -> C = {C:.4f}; partial sums at x = 1 grow like {slope:.4f} ln N (C/2 = {C/2:.4f}); S(N=2e5) = {S[-1]:.3f}")
xx = np.linspace(0, 0.999, 400)
k, a = frob_coeffs(2.5, 0, 20_000)
y25 = np.polyval(a[::-1], xx**2)                       # even series in powers of x^2
axs[1].plot(xx, y25 / y25[0], label="l = 2.5 even series (normalized to 1 at x = 0)")
axs[1].plot(xx, eval_legendre(2, xx) / eval_legendre(2, 0), "--", label="l = 2: P_2(x)/P_2(0)")
axs[0].set_xlabel("N (highest power)"); axs[0].set_ylabel("partial sum at x = 1"); axs[0].legend()
axs[1].set_xlabel("x = cos(theta)"); axs[1].legend(); axs[1].set_title("non-integer l blows up at the pole")
fig.tight_layout(); fig.savefig(FIG / "frobenius_divergence.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 2: Images (§3.3, Examples 4.1, 4.2, Exercises A2, B2)


# %%
def V_pts(qs, pos, pts):
    """Potential k sum_i q_i / |r - r_i| (EM-2 (3.8)) at field points pts (M, 3)."""
    return K * (qs[None, :] / np.linalg.norm(pts[:, None, :] - pos[None, :, :], axis=2)).sum(1)


def E_pts(qs, pos, pts):
    """Field k sum_i q_i (r - r_i) / |r - r_i|^3 (EM-2 (3.2)) at field points pts (M, 3)."""
    sep = pts[:, None, :] - pos[None, :, :]
    return K * (qs[None, :, None] * sep / np.linalg.norm(sep, axis=2)[..., None] ** 3).sum(1)


def plane_numbers(q, d):
    """Example 4.1 quantities for q at height d above a grounded plane: (3.3)-(3.5)."""
    sig = lambda s: -q * d / (2 * np.pi * (s**2 + d**2) ** 1.5)
    Qind = integrate.quad(lambda u: sig(d * u) * 2 * np.pi * u * d**2, 0, np.inf)[0]            # s = d u
    Fsurf = integrate.quad(lambda u: sig(d * u) ** 2 / (2 * EPS0) * 2 * np.pi * u * d**2, 0, np.inf)[0]
    Wwork = -integrate.quad(lambda u: K * q**2 / (4 * d * (1 + u) ** 2), 0, np.inf)[0]        # z = d (1 + u)
    return sig, Qind, Fsurf, Wwork


q, d = 1 * nC, 0.01
sig_pl, Qind, Fsurf, Wwork = plane_numbers(q, d)
qs_pl, pos_pl = np.array([q, -q]), np.array([[0, 0, d], [0, 0, -d]])
sv = np.linspace(0, 5 * d, 21); hz = 1e-6
Vz = [V_pts(qs_pl, pos_pl, np.c_[sv, 0 * sv, 0 * sv + j * hz]) for j in range(3)]
sig_fd = -EPS0 * (-3 * Vz[0] + 4 * Vz[1] - Vz[2]) / (2 * hz)           # sigma = -eps0 dV/dz at the plane
assert np.allclose(sig_fd, sig_pl(sv), rtol=1e-6, atol=0) and np.abs(Vz[0]).max() < 1e-9
Eplane = E_pts(qs_pl, pos_pl, np.array([[0, 0, 0.0]]))[0]
assert np.isclose(Qind, -q, rtol=1e-9, atol=0) and np.isclose(Fsurf, K * q**2 / (4 * d**2), rtol=1e-9, atol=0)
assert np.isclose(Wwork, -K * q**2 / (4 * d), rtol=1e-9, atol=0)
print(f"Ex 4.1: E at plane below q = {Eplane[2]:.4e} V/m; sigma(0) = {sig_pl(0)*1e6:.3f}, sigma(d) = {sig_pl(d)*1e6:.3f} uC/m^2; "
      f"Q_ind = {Qind/nC:.6f} nC; Q(<d) = {-q*(1-1/np.sqrt(2))/nC:.3f} nC; half inside s = {np.sqrt(3)*d*100:.2f} cm")
print(f"        F = {K*q**2/(4*d**2):.4e} N (surface-pressure integral {Fsurf:.4e}); W = {Wwork:.4e} J (work integral), "
      f"naive -kq^2/2d = {-K*q**2/(2*d):.4e} J")
q2, d2 = 3 * nC, 0.02
sig2, Q2, F2, W2 = plane_numbers(q2, d2)
print(f"A2: sigma(0) = {sig2(0)*1e6:.4f} uC/m^2, sigma(d) = {sig2(d2)*1e6:.4f}; Q(<2d) = {-q2*(1-1/np.sqrt(5))/nC:.4f} nC; "
      f"F = {F2:.4e} N; work to remove = {-W2:.4e} J")

R, a = 0.1, 0.2
qp, b = -q * R / a, R**2 / a                                            # (3.6)
u = rng.standard_normal((1000, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
qs_sp, pos_sp = np.array([q, qp]), np.array([[0, 0, a], [0, 0, b]])
assert np.abs(V_pts(qs_sp, pos_sp, R * u)).max() < 1e-12 * K * q / R
F_sph = K * q * qp / (a - b) ** 2
assert np.isclose(F_sph, -K * q**2 * R * a / (a**2 - R**2) ** 2, rtol=1e-12, atol=0)        # (3.7)
sig_sph = lambda th_: -q * (a**2 - R**2) / (4 * np.pi * R * (R**2 + a**2 - 2 * a * R * np.cos(th_)) ** 1.5)
thv = np.linspace(0, np.pi, 7); hr = 1e-7
Vr = [V_pts(qs_sp, pos_sp, (R + j * hr) * np.c_[np.sin(thv), 0 * thv, np.cos(thv)]) for j in range(3)]
assert np.allclose(-EPS0 * (-3 * Vr[0] + 4 * Vr[1] - Vr[2]) / (2 * hr), sig_sph(thv), rtol=1e-5, atol=0)
Qsph = integrate.quad(lambda th_: sig_sph(th_) * 2 * np.pi * R**2 * np.sin(th_), 0, np.pi)[0]
assert np.isclose(Qsph, qp, rtol=1e-9, atol=0)
F_neut = K * q * (qp / (a - b) ** 2 + (-qp) / a**2)                    # B2(c): add +qR/a at the centre
print(f"Ex 4.2: q' = {qp/nC:.2f} nC at b = {b*100:.1f} cm; max|V| on sphere = {np.abs(V_pts(qs_sp, pos_sp, R*u)).max():.1e} V; "
      f"F = {F_sph:.4e} N (plane formula, same gap: {-K*q**2/(4*(a-R)**2):.4e} N)")
print(f"        sigma(0) = {sig_sph(0)/nC:.3f} nC/m^2, sigma(pi) = {sig_sph(np.pi)/nC:.4f} nC/m^2, ratio {sig_sph(0)/sig_sph(np.pi):.2f}; "
      f"integral of sigma = {Qsph/nC:.6f} nC")
print(f"B2(c): neutral sphere F = {F_neut:.4e} N; far limit -2kq^2R^3/a^5 = {-2*K*q**2*R**3/a**5:.4e} N")
for aa in (1.001 * R, 1.01 * R, 10 * R):
    Fg = -K * q**2 * R * aa / (aa**2 - R**2) ** 2
    Fn = K * q * (-q * R / aa / (aa - R**2 / aa) ** 2 + q * R / aa / aa**2)
    print(f"B2(d): a = {aa/R:.3f} R: grounded F = {Fg:.4e} N, F/(-kq^2/4(a-R)^2) = {Fg/(-K*q**2/(4*(aa-R)**2)):.4f}; "
          f"neutral F = {Fn:.4e} N, F/(-2kq^2R^3/a^5) = {Fn/(-2*K*q**2*R**3/aa**5):.4f}")
    assert np.isclose(Fn, -K * q**2 * R**3 * (2 * aa**2 - R**2) / (aa**3 * (aa**2 - R**2) ** 2), rtol=1e-12)   # B2(c) closed form

fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
gx, gz = np.meshgrid(np.linspace(-3 * d, 3 * d, 120), np.linspace(1e-5, 4 * d, 90))
Eg = E_pts(qs_pl, pos_pl, np.c_[gx.ravel(), 0 * gx.ravel(), gz.ravel()])
axs[0].streamplot(gx / d, gz / d, Eg[:, 0].reshape(gx.shape), Eg[:, 2].reshape(gx.shape), density=1.4, color="C0")
axs[0].axhspan(-0.3, 0, color="0.6"); axs[0].plot(0, 1, "ro"); axs[0].set_title("q above a grounded plane (z > 0 only)")
gx, gz = np.meshgrid(np.linspace(-3 * R, 3 * R, 150), np.linspace(-2.5 * R, 3 * R, 140))
Eg = E_pts(qs_sp, pos_sp, np.c_[gx.ravel(), 0 * gx.ravel(), gz.ravel()])
inside = (gx**2 + gz**2 < R**2)
Ex_, Ez_ = [np.ma.array(c.reshape(gx.shape), mask=inside) for c in (Eg[:, 0], Eg[:, 2])]
axs[1].streamplot(gx / R, gz / R, Ex_, Ez_, density=1.4, color="C0")
axs[1].add_patch(plt.Circle((0, 0), 1, color="0.6")); axs[1].plot(0, a / R, "ro"); axs[1].set_aspect("equal")
axs[1].set_title("q outside a grounded sphere, a = 2R")
axs[2].plot(sv / d, sig_fd * 1e6, "o", label="-eps0 dV/dz from (3.2)"); axs[2].plot(sv / d, sig_pl(sv) * 1e6, label="(3.3)")
axs[2].set_xlabel("s / d"); axs[2].set_ylabel("sigma (uC/m^2)"); axs[2].legend()
fig.tight_layout(); fig.savefig(FIG / "images.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 3: The slot by Fourier series and by relaxation (§3.4, Example 4.3, B3)
# `sor` is red-black successive over-relaxation (LA-13) for any 5-point stencil
# sum_nb c_nb (V_nb - V) = 0; constant coefficients give the discrete Laplacian.


# %%
def sor(V, fixed, cE, cW, cN, cS, omega, tol=1e-11, maxit=40_000):
    """Red-black SOR for sum_nb c_nb (V_nb - V_ij) = 0 at free nodes; E/W = column j+-1, N/S = row i+-1."""
    ii, jj = np.indices(V.shape)
    csum = cE + cW + cN + cS
    colors = [~fixed & ((ii + jj) % 2 == c) for c in (0, 1)]
    for it in range(1, maxit + 1):
        change = 0.0
        for m in colors:
            Vp = np.pad(V, 1, mode="edge")
            avg = (cE * Vp[1:-1, 2:] + cW * Vp[1:-1, :-2] + cN * Vp[2:, 1:-1] + cS * Vp[:-2, 1:-1]) / csum
            dV = omega * (avg - V)
            V[m] += dV[m]
            change = max(change, np.abs(dV[m]).max())
        if change < tol:
            return V, it
    raise RuntimeError("SOR did not converge")


def slot_series(xv, yv, coef, a=1.0):
    """(3.8): V = sum_n C_n exp(-n pi x / a) sin(n pi y / a), coef[n-1] = C_n."""
    n = np.arange(1, len(coef) + 1)[:, None]
    return (coef[:, None] * np.exp(-n * np.pi * xv / a) * np.sin(n * np.pi * yv / a)).sum(0)


def solve_slot(V0y, h=1 / 40, Lx=4.0):
    """Finite-difference slot 0<y<1, 0<x<Lx: V = V0(y) at x = 0, 0 on the plates and at x = Lx."""
    ny, nx = int(round(1 / h)) + 1, int(round(Lx / h)) + 1
    V = np.zeros((ny, nx)); V[1:-1, 0] = V0y(np.linspace(0, 1, ny)[1:-1])
    fixed = np.zeros_like(V, bool); fixed[[0, -1], :] = True; fixed[:, [0, -1]] = True
    rhoJ = 0.5 * (np.cos(np.pi / (ny - 1)) + np.cos(np.pi / (nx - 1)))
    one = np.ones_like(V)
    return sor(V, fixed, one, one, one, one, 2 / (1 + np.sqrt(1 - rhoJ**2)))


def sine_coeffs(V0y, nmax, a=1.0):
    """(3.9): C_n = (2/a) int_0^a V0(y) sin(n pi y / a) dy, n = 1..nmax (QAWO oscillatory quadrature)."""
    return np.array([2 / a * integrate.quad(V0y, 0, a, weight="sin", wvar=n * np.pi / a)[0] for n in range(1, nmax + 1)])


Cn = sine_coeffs(lambda y: 1.0, 399)                                     # V0 = 1
assert np.allclose(Cn[:3], [4 / np.pi, 0, 4 / (3 * np.pi)], atol=1e-12)
closed = lambda xv, yv: 2 / np.pi * np.arctan(np.sin(np.pi * yv) / np.sinh(np.pi * xv))
partial = [slot_series(0.5, 0.5, Cn[:N])[0] for N in (1, 3, 5)]
assert np.isclose(slot_series(0.5, 0.5, Cn)[0], closed(0.5, 0.5), atol=1e-12)
print(f"Ex 4.3: C1..C3 = {Cn[0]:.4f}, {Cn[1]:.4f}, {Cn[2]:.4f} V0; V(a/2, a/2) with 1, 2, 3 odd terms = "
      + ", ".join(f"{v:.5f}" for v in partial) + f"; closed form {closed(0.5, 0.5):.5f} V0")
for xe in (0.5, 0.1, 0.01):
    full = closed(xe, 0.5)
    Nneed = next(N for N in range(1, 400) if abs(slot_series(xe, 0.5, Cn[:N])[0] - full) < 1e-4)
    print(f"        terms (n <= N) needed for 1e-4 accuracy at x = {xe} a, y = a/2: N = {Nneed}")
Vfd, its = solve_slot(lambda y: np.ones_like(y))
yg, xg = np.meshgrid(np.linspace(0, 1, Vfd.shape[0]), np.linspace(0, 4, Vfd.shape[1]), indexing="ij")
away = xg >= 0.25
err_slot = np.abs(Vfd - closed(np.maximum(xg, 1e-12), yg))[away].max()
assert abs(Vfd[20, 20] - closed(0.5, 0.5)) < 1e-3 and err_slot < 2e-3
print(f"        SOR (h = a/40, {its} sweeps): V(a/2, a/2) = {Vfd[20, 20]:.5f} V0; max |FD - exact| for x >= a/4: {err_slot:.1e} V0")
fig, ax = plt.subplots(figsize=(9, 3.2))
lv = [0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8]
ax.contour(xg, yg, closed(np.maximum(xg, 1e-12), yg), lv, colors="k")
ax.contour(xg, yg, Vfd, lv, colors="r", linestyles="--")
ax.set_xlim(0, 2); ax.set_aspect("equal"); ax.set_title("slot: Fourier series (black) vs SOR (red dashed)")
fig.tight_layout(); fig.savefig(FIG / "slot.png", dpi=120); plt.close(fig)
nb = np.arange(1, 8, 2)                                                  # B3: square box, V0 on x = b = a
terms = 4 / (np.pi * nb) * np.sinh(nb * np.pi / 2) / np.sinh(nb * np.pi) * np.sin(nb * np.pi / 2)
print("B3: box centre partial sums (n = 1; 1,3; 1,3,5; 1..7):", np.round(np.cumsum(terms), 6))
assert np.isclose(np.cumsum(terms)[-1], 0.25, atol=1e-6)

# %% [markdown]
# ## Part 4: Conducting sphere in a uniform field by relaxation (Example 4.4)
# r^2 times (3.11) in rho = ln r is  e^-rho d_rho(e^rho d_rho V) + (1/sin) d_theta(sin d_theta V) = 0.
# Multiplied by e^rho sin(theta) it is a symmetric 5-point stencil; theta sits at
# cell centres, so the pole faces carry sin = 0 and need no special treatment.
# V = 0 on the sphere (r = R = 1), V = -E0 r cos(theta) at r = 50 R; E0 = 1. The
# discretization error is about 1e-3 relative to the local |V|, so it is measured that way.

# %%
Nr, Nt, rmax = 161, 72, 50.0
rho = np.linspace(0, np.log(rmax), Nr); hrho = rho[1] - rho[0]
tht = (np.arange(Nt) + 0.5) * np.pi / Nt; dth = np.pi / Nt
TH, RHO = np.meshgrid(tht, rho, indexing="ij")
cE = np.sin(TH) * np.exp(RHO + hrho / 2) / hrho**2
cW = np.sin(TH) * np.exp(RHO - hrho / 2) / hrho**2
cN = np.exp(RHO) * np.sin(np.minimum(TH + dth / 2, np.pi)) / dth**2
cS = np.exp(RHO) * np.sin(np.maximum(TH - dth / 2, 0)) / dth**2
Vs = np.zeros((Nt, Nr)); Vs[:, -1] = -rmax * np.cos(tht)
fixed = np.zeros_like(Vs, bool); fixed[:, [0, -1]] = True
Vs, its = sor(Vs, fixed, cE, cW, cN, cS, omega=1.95, tol=1e-10)
rr = np.exp(RHO)
Vex = -(rr - 1 / rr**2) * np.cos(TH)                                     # (4.2) with E0 = R = 1
err_V = np.abs(Vs - Vex)[rr <= 3].max()                                # absolute error for R <= r <= 3R
rel_V = (np.abs(Vs - Vex) / np.maximum(np.abs(Vex), 1)).max()           # error relative to max(|V|, E0 R), whole grid
sig_num = -(4 * Vs[:, 1] - Vs[:, 2]) / (2 * hrho)                      # sigma/eps0 = -dV/dr at r = R (V = 0 there)
err_sig = np.abs(sig_num - 3 * np.cos(tht)).max() / 3
p_num = (np.cos(tht) * sig_num * 2 * np.pi * np.sin(tht) * dth).sum() / (4 * np.pi)   # p/(4 pi eps0 R^3 E0)
Q_num = (sig_num * 2 * np.pi * np.sin(tht) * dth).sum()
assert err_V < 5e-3 and rel_V < 2e-3 and err_sig < 2e-3 and abs(p_num - 1) < 2e-3 and abs(Q_num) < 1e-10
print(f"Ex 4.4 relaxation ({its} sweeps, {Nt}x{Nr} grid): max |V - (4.2)| = {err_V:.1e} E0 R for r <= 3R (relative {rel_V:.1e} everywhere); "
      f"max sigma error {err_sig:.1e} of 3 eps0 E0; p = {p_num:.4f} x 4 pi eps0 R^3 E0; Q = {Q_num:.1e}")
Rb, Eb = 0.01, 1e5
print(f"        R = 1 cm, E0 = 1e5 V/m: sigma_max = {3*EPS0*Eb*1e6:.3f} uC/m^2; pole field {3*Eb:.1e} V/m; "
      f"upper-hemisphere charge {3*np.pi*EPS0*Eb*Rb**2/nC:.3f} nC; p = {4*np.pi*EPS0*Rb**3*Eb:.4e} C m")
fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
for sgn in (1, -1):
    Xs, Zs = sgn * rr * np.sin(TH), rr * np.cos(TH)
    axs[0].contour(Xs, Zs, Vs, np.linspace(-3, 3, 13), colors="r", linestyles="--")
    axs[0].contour(Xs, Zs, Vex, np.linspace(-3, 3, 13), colors="k", linewidths=0.6)
axs[0].add_patch(plt.Circle((0, 0), 1, color="0.6")); axs[0].set_xlim(-3, 3); axs[0].set_ylim(-3, 3)
axs[0].set_aspect("equal"); axs[0].set_title("equipotentials: SOR (red) vs (4.2) (black)")
axs[1].plot(np.degrees(tht), sig_num, "o", ms=3, label="relaxation"); axs[1].plot(np.degrees(tht), 3 * np.cos(tht), label="3 cos(theta)")
axs[1].set_xlabel("theta (deg)"); axs[1].set_ylabel("sigma / (eps0 E0)"); axs[1].legend()
fig.tight_layout(); fig.savefig(FIG / "sphere_uniform_field.png", dpi=120); plt.close(fig)
print(f"Ex 4.5: R = 5 cm, sigma0 = 1 uC/m^2: E_in = {1e-6/(3*EPS0):.4e} V/m (along -z); p = {4*np.pi/3*0.05**3*1e-6:.4e} C m")
Qc, Rq, Eq_ = 0.1 * nC, 0.01, 1e5                                         # B5: charged sphere in a field
cth0 = -Qc / (12 * np.pi * EPS0 * Rq**2 * Eq_)
print(f"B5: threshold |Q| = {12*np.pi*EPS0*Rq**2*Eq_/nC:.4f} nC; Q = 0.1 nC: sigma = 0 at theta = {np.degrees(np.arccos(cth0)):.1f} deg; "
      f"field at north pole {3*Eq_ + K*Qc/Rq**2:.4e} V/m, south pole {-3*Eq_ + K*Qc/Rq**2:.4e} V/m (radial)")

# %% [markdown]
# ## Part 5: The multipole expansion (§3.9, Example 4.6, A5, B6)
# Error of the expansion truncated after l = 0, 1, 2 falls like r^-2, r^-3, r^-4.


# %%
def moments(qs, pos):
    """Monopole Q, dipole p (3.29), traceless quadrupole Q_ij (3.30) about the origin."""
    r2 = (pos**2).sum(1)
    return qs.sum(), qs @ pos, 3 * np.einsum("n,ni,nj->ij", qs, pos, pos) - np.eye(3) * (qs @ r2)


def V_multipole(Q, p, Qij, pts, lmax):
    """Truncated expansion (3.28) through order lmax <= 2, in tensor form."""
    r = np.linalg.norm(pts, axis=1); rh = pts / r[:, None]
    V = K * Q / r
    if lmax >= 1:
        V = V + K * (rh @ p) / r**2
    if lmax >= 2:
        V = V + K * np.einsum("mi,ij,mj->m", rh, Qij, rh) / (2 * r**3)
    return V


def fib_sphere(n):
    """n roughly uniform unit vectors (Fibonacci lattice)."""
    i = np.arange(n) + 0.5; thf, phf = np.arccos(1 - 2 * i / n), np.pi * (1 + 5**0.5) * i
    return np.c_[np.sin(thf) * np.cos(phf), np.sin(thf) * np.sin(phf), np.cos(thf)]


qs_c = rng.integers(1, 4, 6) * rng.choice([-1, 1], 6) * nC
pos_c = rng.standard_normal((6, 3)); pos_c *= (rng.uniform(0, 1, 6) ** (1 / 3) / np.linalg.norm(pos_c, axis=1))[:, None]
Qc_, pc_, Qij_c = moments(qs_c, pos_c)
assert abs(np.trace(Qij_c)) < 1e-12 * np.abs(Qij_c).max() and np.allclose(Qij_c, Qij_c.T)
dirs = fib_sphere(200); rs = np.logspace(np.log10(3), 2, 15); errs = np.zeros((3, rs.size))
for j, rv in enumerate(rs):
    Vexact = V_pts(qs_c, pos_c, rv * dirs)
    for L in range(3):
        errs[L, j] = np.sqrt(np.mean((V_multipole(Qc_, pc_, Qij_c, rv * dirs, L) - Vexact) ** 2))
calpha = (dirs @ pos_c.T) / np.linalg.norm(pos_c, axis=1)               # Legendre form of the l = 2 term
V2leg = K * (qs_c * (pos_c**2).sum(1) * eval_legendre(2, calpha)).sum(1) / 10.0**3
assert np.allclose(V2leg, V_multipole(0, 0 * pc_, Qij_c, 10 * dirs, 2), rtol=1e-10, atol=0)
slopes = [np.polyfit(np.log(rs[rs >= 10]), np.log(errs[L, rs >= 10]), 1)[0] for L in range(3)]
assert all(abs(s_ - e_) < 0.05 for s_, e_ in zip(slopes, (-2, -3, -4)))
shift = np.array([0.3, -0.2, 0.5]); assert np.allclose(moments(qs_c, pos_c - shift)[1], pc_ - Qc_ * shift, rtol=1e-12, atol=0)   # (3.31)
assert np.abs(moments(qs_c, pos_c - pc_ / Qc_)[1]).max() < 1e-12 * np.linalg.norm(pc_)      # origin at the centre of charge
print(f"Multipoles: Q = {Qc_/nC:.0f} nC, |p| = {np.linalg.norm(pc_)/nC:.3f} nC m, tr Q_ij = {np.trace(Qij_c):.1e}; "
      f"error slopes after l = 0, 1, 2: " + ", ".join(f"{s_:.3f}" for s_ in slopes))
fig, ax = plt.subplots(figsize=(6, 4.5))
for L in range(3):
    ax.loglog(rs, errs[L], "o-", label=f"truncated after l = {L} (slope {slopes[L]:.2f})")
ax.set_xlabel("r / d"); ax.set_ylabel("RMS error in V (V)"); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "multipole_errors.png", dpi=120); plt.close(fig)


def V_axis_series(qs, zs, r, cth, lmax):
    """(3.28) for charges on the z-axis: k sum_l sum_i q_i z_i^l P_l(cos theta) / r^(l+1)."""
    return sum(K * (qs * zs**l).sum() * eval_legendre(l, cth) / r ** (l + 1) for l in range(lmax + 1))


for name, qs_, zs_ in (("(a) +-q at +-d", np.array([1.0, -1.0]), np.array([1.0, -1.0])),
                       ("(b) linear quadrupole", np.array([1.0, -2.0, 1.0]), np.array([1.0, 0.0, -1.0]))):
    ex = V_pts(qs_, np.c_[0 * zs_, 0 * zs_, zs_], np.array([[0, 0, 5.0]]))[0] / K
    ser = [V_axis_series(qs_, zs_, 5.0, 1.0, L) / K for L in range(7)]
    Qm, pm, Qijm = moments(qs_, np.c_[0 * zs_, 0 * zs_, zs_])
    print(f"Ex 4.6 {name}: Q = {Qm:g}, p_z = {pm[2]:g} qd, Q_zz = {Qijm[2, 2]:g} qd^2, Q_xx = {Qijm[0, 0]:g}; on axis at r = 5d "
          f"(units kq/d): exact {ex:.6f}; through l = 1..6: " + ", ".join(f"{v:.6f}" for v in ser[1:]))
q5 = np.array([3.0, -1.0]); z5 = np.array([1.0, -1.0])                  # A5
Q5, p5, _ = moments(q5, np.c_[0 * z5, 0 * z5, z5])
print(f"A5: Q = {Q5:g} q, p = {p5[2]:g} qa about origin, {moments(q5, np.c_[0*z5, 0*z5, z5 - 1])[1][2]:g} qa about (0,0,a); "
      f"centre of charge z = {p5[2]/Q5:g} a")
qs6 = np.array([1.0, 1.0, -1.0, -1.0]); pos6 = np.array([[1.0, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]])   # B6
Q6, p6, Qij6 = moments(qs6, pos6)
assert np.allclose(Qij6, np.diag([6, -6, 0])) and Q6 == 0 and np.allclose(p6, 0)
pt6 = np.array([[10.0, 0, 0]])
ex6, q6v = V_pts(qs6, pos6, pt6)[0] / K, V_multipole(0, p6, Qij6, pt6, 2)[0] / K
l4_6 = (qs6 * np.linalg.norm(pos6, axis=1) ** 4 * eval_legendre(4, pos6[:, 0] / np.linalg.norm(pos6, axis=1))).sum() / 10.0**5   # l = 4 term of (3.28)
assert np.isclose(l4_6, 1.25e-5) and abs(q6v + l4_6 - ex6) < 4e-7                     # odd l vanish by inversion symmetry
print(f"B6: Q_ij = diag{tuple(int(v) for v in np.diag(Qij6))} qa^2; at r = 10a in-plane, phi = 0: exact {ex6:.7f}, quadrupole {q6v:.7f} kq/a, "
      f"ratio {ex6/q6v:.5f}; l = 4 term {l4_6:.3e}, quadrupole + l = 4 = {q6v + l4_6:.7f}")

# %% [markdown]
# ## Part 6: Spherical harmonics (§3.8): Griffiths' table vs SciPy, orthonormality, plots

# %%
table = {(0, 0): 1 / sp.sqrt(4 * sp.pi) + 0 * th, (1, 0): sp.sqrt(3 / (4 * sp.pi)) * sp.cos(th),
         (1, 1): -sp.sqrt(3 / (8 * sp.pi)) * sp.sin(th) * sp.exp(sp.I * ph),
         (2, 0): sp.sqrt(5 / (16 * sp.pi)) * (3 * sp.cos(th) ** 2 - 1),
         (2, 1): -sp.sqrt(15 / (8 * sp.pi)) * sp.sin(th) * sp.cos(th) * sp.exp(sp.I * ph),
         (2, 2): sp.sqrt(15 / (32 * sp.pi)) * sp.sin(th) ** 2 * sp.exp(2 * sp.I * ph)}
tt, pp = rng.uniform(0, np.pi, 50), rng.uniform(0, 2 * np.pi, 50)
for (l, m), Y in table.items():
    angl = sp.diff(sp.sin(th) * sp.diff(Y, th), th) / sp.sin(th) + sp.diff(Y, ph, 2) / sp.sin(th) ** 2
    assert sp.simplify(angl + l * (l + 1) * Y) == 0                   # eigenvalue -l(l+1)
    f = sp.lambdify((th, ph), Y, "numpy")
    assert np.allclose(f(tt, pp), sph_harm_y(l, m, tt, pp))
    assert np.allclose(np.conj(f(tt, pp)) * (-1) ** m, sph_harm_y(l, -m, tt, pp))   # Y_l^-m = (-1)^m conj(Y_l^m)
xg_, wg_ = roots_legendre(12); phq = np.arange(24) * 2 * np.pi / 24
THq, PHq = np.meshgrid(np.arccos(xg_), phq, indexing="ij"); Wq = np.outer(wg_, np.full(24, 2 * np.pi / 24))
lm = [(l, m) for l in range(3) for m in range(-l, l + 1)]
Ys = np.array([sph_harm_y(l, m, THq, PHq) for l, m in lm])
gram = np.einsum("aij,bij,ij->ab", Ys.conj(), Ys, Wq)
assert np.allclose(gram, np.eye(9), atol=1e-12)
print("Spherical harmonics: Griffiths' table = SciPy (Condon-Shortley), angular-Laplacian eigenvalues -l(l+1), 9x9 Gram matrix = I")
fig = plt.figure(figsize=(12, 7))
TT, PP = np.meshgrid(np.linspace(0, np.pi, 60), np.linspace(0, 2 * np.pi, 80), indexing="ij")
for kk, (l, m) in enumerate([(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)]):
    rad = np.abs(sph_harm_y(l, m, TT, PP)) ** 2
    ax = fig.add_subplot(2, 3, kk + 1, projection="3d")
    ax.plot_surface(rad * np.sin(TT) * np.cos(PP), rad * np.sin(TT) * np.sin(PP), rad * np.cos(TT),
                    facecolors=plt.cm.viridis(rad / rad.max()), linewidth=0)
    lim = rad.max(); ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_title(f"|Y_{l}^{m}|^2"); ax.set_axis_off()
fig.tight_layout(); fig.savefig(FIG / "spherical_harmonics.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 7: Dipoles in external fields (§3.10, A6, B7, C3)
# A physical dipole +-q at r0 +- d/2 in the field of a point charge Q at the origin.
# Energy, torque and force approach (3.32)-(3.34) with errors O(d^2); (3.35) is dipole-dipole.

# %%
Qs, r0 = 1 * nC, np.array([0.03, 0.01, 0.05]); phat = np.array([1.0, 2.0, -0.5]) / np.linalg.norm([1, 2, -0.5])
pmag = 1e-12
rn = np.linalg.norm(r0)
E0v = K * Qs * r0 / rn**3
F_th = K * Qs * (pmag * phat / rn**3 - 3 * (pmag * phat @ r0) * r0 / rn**5)   # (p.grad)E for a point-charge field
prev = None
for dd in (2e-3, 1e-3, 5e-4):
    qd, dv = pmag / dd, dd * phat
    ends = np.array([r0 + dv / 2, r0 - dv / 2]); qq = np.array([qd, -qd])
    Ef = E_pts(np.array([Qs]), np.zeros((1, 3)), ends)
    errs3 = np.array([np.linalg.norm(qq @ Ef - F_th) / np.linalg.norm(F_th),
                      np.linalg.norm(np.cross(dv / 2, qd * Ef[0]) + np.cross(-dv / 2, -qd * Ef[1]) - np.cross(pmag * phat, E0v))
                      / np.linalg.norm(np.cross(pmag * phat, E0v)),
                      abs(qq @ V_pts(np.array([Qs]), np.zeros((1, 3)), ends) + pmag * phat @ E0v) / abs(pmag * phat @ E0v)])
    if prev is not None:
        assert np.all((prev / errs3 > 3.5) & (prev / errs3 < 4.5))   # halving d quarters every error
    prev = errs3
    print(f"Dipole d = {dd*1e3:.1f} mm: relative errors of force, torque, energy = " + ", ".join(f"{e:.2e}" for e in errs3))
p1, p2 = 1e-12 * np.array([0.3, 0.5, 0.8]), 1e-12 * np.array([-0.6, 0.2, 0.4])      # dipole-dipole (3.35)
rv12, dd = np.array([0.02, -0.03, 0.05]), 1e-5
c1 = np.array([-p1, p1]) / np.linalg.norm(p1) * dd / 2; c2 = rv12 + np.array([-p2, p2]) / np.linalg.norm(p2) * dd / 2
q1s, q2s = np.array([-1, 1]) * np.linalg.norm(p1) / dd, np.array([-1, 1]) * np.linalg.norm(p2) / dd
U_pairs = sum(K * q1s[i] * q2s[j] / np.linalg.norm(c2[j] - c1[i]) for i in range(2) for j in range(2))
rh12 = rv12 / np.linalg.norm(rv12)
U_form = K / np.linalg.norm(rv12) ** 3 * (p1 @ p2 - 3 * (p1 @ rh12) * (p2 @ rh12))
assert np.isclose(U_pairs, U_form, rtol=1e-6, atol=0)
kB, eV, pw = 1.380649e-23, 1.602176634e-19, 1.85 * 3.33564e-30
Uw = -2 * K * pw**2 / 0.3e-9**3
print(f"Water p = {pw:.3e} C m; head-to-tail at 0.3 nm: U = {Uw:.3e} J = {Uw/eV:.3f} eV = {Uw/(kB*300):.2f} kT(300 K)")
xA6 = pw * 1e6 / (kB * 300)
print(f"A6: torque max {pw*1e6:.3e} N m; 2pE = {2*pw*1e6:.3e} J; kT = {kB*300:.3e} J; pE/kT = {xA6:.3e}; "
      f"<cos> = coth(x) - 1/x = {1/np.tanh(xA6) - 1/xA6:.3e} (x/3 = {xA6/3:.3e})")
xn, wn = roots_legendre(48); phn = np.arange(64) * 2 * np.pi / 64      # C3: orientational averages
X1, X2, PH = np.meshgrid(xn, xn, phn, indexing="ij"); W3 = np.einsum("i,j,k->ijk", wn, wn, np.full(64, 1 / 64)) / 4
geom = np.sqrt(1 - X1**2) * np.sqrt(1 - X2**2) * np.cos(PH) - 2 * X1 * X2   # (3.35) with r-hat = z
assert np.isclose((W3 * geom**2).sum(), 2 / 3)
for rk in (0.5e-9, 1e-9, 2e-9):
    U0 = K * pw**2 / rk**3; Ug = U0 * geom; bw = np.exp(-Ug / (kB * 300))
    exact = (W3 * Ug * bw).sum() / (W3 * bw).sum(); lead = -2 / 3 * U0**2 / (kB * 300)
    print(f"C3: r = {rk*1e9:.1f} nm: k p^2/r^3 = {U0:.3e} J = {U0/(kB*300):.3f} kT; <U> = {exact:.4e} J, "
          f"leading order -2U0^2/3kT = {lead:.4e} J, ratio {exact/lead:.4f}")
    if rk >= 1e-9:
        assert abs(exact / lead - 1) < 0.1

# %% [markdown]
# ## Your turn
#
# 1. **Charge between two grounded planes** (Exercise C1). q at height z0 between
#    grounded planes z = 0 and z = L. Build the infinite image series (q at
#    z0 + 2nL, -q at -z0 + 2nL), check V = 0 on both planes, and integrate the
#    induced sigma on each plane. Predict the answer first: C1(b) proves
#    Q_bottom = -q (L - z0)/L.
# 2. **A ring of charge** Q, radius a, in the xy-plane. On the axis V = kQ/sqrt(z^2 + a^2).
#    Expand in powers of a/z (or z/a), read off the Legendre coefficients, and
#    compare the off-axis series with a ring of 4000 point charges, outside
#    (r = 1.5a) and inside (r = 0.5a).
# 3. **Water as three point charges.** O at the origin with -2 delta, H's at
#    0.0957 nm, 104.5 degrees apart, with +delta each, delta chosen so p = 1.85 D.
#    How fast does the relative error of the dipole approximation fall with r?
# 4. **Exercise B8.** Repeat Part 3 with V0(y) = 4 V0 (y/a)(1 - y/a). Compare
#    coefficients, convergence, and the SOR solution.
#
# Reference solutions follow. Write your own before reading them.

# %%
Lp, z0p = 1.0, 0.3
nimg = np.arange(-400, 401)
img_q = np.r_[np.full(nimg.size, 1.0), np.full(nimg.size, -1.0)]
img_z = np.r_[z0p + 2 * nimg * Lp, -z0p + 2 * nimg * Lp]
img_pos = np.c_[0 * img_z, 0 * img_z, img_z]
sv1 = np.linspace(0, 3, 7)
for zp in (0.0, Lp):
    assert np.abs(V_pts(img_q, img_pos, np.c_[sv1, 0 * sv1, 0 * sv1 + zp])).max() < 1e-3 * K
Ez_at = lambda s, zp: E_pts(img_q, img_pos, np.array([[s, 0, zp]]))[0, 2]
Qbot = integrate.quad(lambda s: EPS0 * Ez_at(s, 0.0) * 2 * np.pi * s, 0, 15, limit=200)[0]
Qtop = integrate.quad(lambda s: -EPS0 * Ez_at(s, Lp) * 2 * np.pi * s, 0, 15, limit=200)[0]
assert np.isclose(Qbot, -(Lp - z0p) / Lp, atol=1e-4) and np.isclose(Qtop, -z0p / Lp, atol=1e-4)
print(f"Your turn 1 / C1: z0 = 0.3 L: Q_bottom = {Qbot:.5f} q, Q_top = {Qtop:.5f} q (predicted -0.7, -0.3)")
ar, nring = 1.0, 4000
ring = np.c_[np.cos(2 * np.pi * np.arange(nring) / nring), np.sin(2 * np.pi * np.arange(nring) / nring), np.zeros(nring)]
P2n0 = lambda n: np.prod([-(2 * j - 1) / (2 * j) for j in range(1, n + 1)])   # P_2n(0) = (-1)^n (2n)! / (4^n n!^2), in floats
for rr_, tt_ in ((1.5, np.pi / 3), (0.5, np.pi / 6)):
    num = V_pts(np.full(nring, 1.0 / nring), ring, rr_ * np.array([[np.sin(tt_), 0, np.cos(tt_)]]))[0] / K
    ratio_ = min(rr_, 1 / rr_)
    ser = sum(P2n0(n) * ratio_ ** (2 * n) * eval_legendre(2 * n, np.cos(tt_)) for n in range(60)) / max(rr_, 1.0)
    assert np.isclose(num, ser, rtol=1e-10)
    print(f"Your turn 2: ring, r = {rr_} a, theta = {np.degrees(tt_):.0f} deg: point-charge ring {num:.10f}, Legendre series {ser:.10f} (kQ/a)")
print("        on-axis coefficients P_2n(0) = " + ", ".join(f"{P2n0(n):.4f}" for n in range(5)))
bOH, angH = 0.0957e-9, np.radians(104.5)
delta = pw / (2 * bOH * np.cos(angH / 2))
pos_w = np.array([[0, 0, 0], [bOH * np.sin(angH / 2), 0, bOH * np.cos(angH / 2)], [-bOH * np.sin(angH / 2), 0, bOH * np.cos(angH / 2)]])
qs_w = np.array([-2 * delta, delta, delta])
assert np.allclose(moments(qs_w, pos_w)[1] / pw, [0, 0, 1], atol=1e-12)
Qw, pw_vec, Qij_w = moments(qs_w, pos_w)
rw = np.array([0.5, 1, 2, 4, 8, 16, 32]) * 1e-9
ptsw = rw[:, None] * np.array([0, 0, 1.0])                                # on the dipole axis
Vex_w = V_pts(qs_w, pos_w, ptsw)
err_dip = K * (ptsw @ pw_vec) / rw**3 / Vex_w - 1                          # signed relative error of (3.29)
err_quad = V_multipole(0.0, pw_vec, Qij_w, ptsw, 2) / Vex_w - 1            # with the quadrupole term (3.30)
coef_pred = -Qij_w[2, 2] / (2 * pw)                                        # leading error: -Q_zz / (2 p r), from (3.30)
assert abs(err_dip[-1] * rw[-1] / coef_pred - 1) < 0.05                    # r * err -> -Q_zz/(2p) (within 2% at 32 nm)
ratios_q = err_quad[2:-1] / err_quad[3:]
assert np.all((ratios_q > 3.9) & (ratios_q < 4.2))                        # dipole + quadrupole error falls like r^-2
print(f"Your turn 3: delta = {delta/1.602176634e-19:.3f} e; Q_zz = {Qij_w[2, 2]:.3e} C m^2, -Q_zz/(2p) = {coef_pred*1e9:.5f} nm")
print("        r (nm):                 " + "  ".join(f"{v*1e9:9.1f}" for v in rw))
print("        dipole rel. error:      " + "  ".join(f"{v:9.5f}" for v in err_dip))
print("        r x dipole error (nm):  " + "  ".join(f"{v*1e9:9.5f}" for v in err_dip * rw))
print("        dipole+quadrupole err.: " + "  ".join(f"{v:9.2e}" for v in err_quad))
print("        (the octupole cancels the quadrupole near 0.7 nm, so the 1/r law only shows for r >= 2 nm)")
Cp = sine_coeffs(lambda y: 4 * y * (1 - y), 59)
nn = np.arange(1, 60)
assert np.allclose(Cp, np.where(nn % 2, 32 / (nn * np.pi) ** 3, 0), atol=1e-12)
Vfd8, its8 = solve_slot(lambda y: 4 * y * (1 - y))
cen8 = slot_series(0.5, 0.5, Cp)[0]
assert abs(Vfd8[20, 20] - cen8) < 1e-3
print(f"Your turn 4 / B8: C1, C3, C5 = {Cp[0]:.5f}, {Cp[2]:.5f}, {Cp[4]:.6f} V0; V(a/2, a/2) = {cen8:.5f} (one term "
      f"{slot_series(0.5, 0.5, Cp[:1])[0]:.5f}); SOR {Vfd8[20, 20]:.5f} ({its8} sweeps)")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.perf_counter() - T0:.1f} s. Figures in {FIG}")
