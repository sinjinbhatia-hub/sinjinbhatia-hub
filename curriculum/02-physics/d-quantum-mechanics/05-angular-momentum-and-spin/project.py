# %% [markdown]
# # QM-5 Project: Angular Momentum Matrices, Spherical Harmonics, Spin Dynamics, Addition
#
# This script builds the angular-momentum matrices for any j from the ladder
# matrix elements (notes (3.14)), checks the spherical harmonics against the
# differential operators L_z and L^2 (3.17)-(3.19), checks the Pauli algebra,
# runs Larmor precession and Rabi oscillation (Schrodinger and Heisenberg),
# checks the SU(2) -> SO(3) map and the 2pi sign, adds angular momenta with
# np.kron, and generalizes QM-1's Stern-Gerlach simulator to any spin j.
# Every number quoted in the notes' worked examples and in solutions.md is
# asserted or printed here. Units: hbar = 1 except in Part 5 (SI).
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.constants import hbar as HBAR, physical_constants as PC, m_e, c as C_LIGHT
from scipy.integrate import solve_ivp
from scipy.linalg import expm, null_space
from scipy.special import sph_harm_y
from sympy.physics.quantum.cg import CG

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.time()


def comm(A, B):
    return A @ B - B @ A


# %% [markdown]
# ## Part 1: J matrices for any j from the ladder elements (§3.2, eq. (3.14))
# Basis |j, m> ordered m = j, j-1, ..., -j. J_+|j,m> = sqrt(j(j+1) - m(m+1)) |j,m+1>.

# %%
def jmat(j):
    """Return (Jx, Jy, Jz, Jp, Jm) for spin j, hbar = 1, from (3.14)."""
    m = j - np.arange(int(round(2 * j)) + 1)
    Jz = np.diag(m).astype(complex)
    Jp = np.zeros_like(Jz)
    for k in range(1, len(m)):                       # <m+1| J_+ |m>, row k-1, column k
        Jp[k - 1, k] = np.sqrt(j * (j + 1) - m[k] * (m[k] + 1))
    Jm = Jp.conj().T
    return (Jp + Jm) / 2, (Jp - Jm) / 2j, Jz, Jp, Jm


for j in (0.5, 1, 1.5, 2):
    Jx, Jy, Jz, Jp, Jm = jmat(j)
    I = np.eye(len(Jz))
    J2 = Jx @ Jx + Jy @ Jy + Jz @ Jz
    assert np.allclose(comm(Jx, Jy), 1j * Jz) and np.allclose(comm(Jy, Jz), 1j * Jx)
    assert np.allclose(comm(Jz, Jx), 1j * Jy)
    assert np.allclose(J2, j * (j + 1) * I)                          # (3.12)
    for A in (Jx, Jy, Jz):
        assert np.allclose(comm(J2, A), 0)                           # Casimir (3.5)
    assert np.allclose(comm(Jz, Jp), Jp) and np.allclose(comm(Jz, Jm), -Jm)   # (3.7)
    assert np.allclose(Jm @ Jp, J2 - Jz @ Jz - Jz)                   # (3.8)
print("Part 1: commutators, J^2 = j(j+1), Casimir checked for j = 1/2, 1, 3/2, 2")

Sx1, Sy1, Sz1, _, _ = jmat(1)
print("spin-1 Sx*sqrt2 =\n", np.round(np.sqrt(2) * Sx1.real, 12))
assert np.allclose(np.sqrt(2) * Sx1, [[0, 1, 0], [1, 0, 1], [0, 1, 0]])
assert np.allclose(np.sqrt(2) * Sy1, [[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]])

# Example 4.2: the termination inequality <Jx^2 + Jy^2> = lambda - m^2 >= 0 (lambda = 2, l = 1)
lam = 2.0
norms = {mm: lam - mm * (mm + 1) for mm in (1, 2)}
print(f"Ex 4.2: ||L+|1,m>||^2 = lambda - m(m+1): m=1 -> {norms[1]}, m=2 -> {norms[2]}")
assert norms[1] == 0 and norms[2] == -4

# %% [markdown]
# ## Part 2: Spherical harmonics vs L_z, L^2, L_- by finite differences (§3.3)
# SciPy's sph_harm_y(l, m, theta, phi) has theta polar and the Condon-Shortley phase.

# %%
def Y(l, m):
    return lambda th, ph: sph_harm_y(l, m, th, ph)


def Lz_fd(f, th, ph, h=1e-4):
    """L_z = -i d/dphi (3.17)"""
    return -1j * (f(th, ph + h) - f(th, ph - h)) / (2 * h)


def L2_fd(f, th, ph, h=1e-3):
    """L^2 = -[ (1/sin) d_th (sin d_th) + (1/sin^2) d_ph^2 ]  (3.19)"""
    f0 = f(th, ph)
    d1 = (f(th + h, ph) - f(th - h, ph)) / (2 * h)
    d2 = (f(th + h, ph) - 2 * f0 + f(th - h, ph)) / h**2
    dp = (f(th, ph + h) - 2 * f0 + f(th, ph - h)) / h**2
    return -(d2 + d1 / np.tan(th) + dp / np.sin(th)**2)


def Lpm_fd(f, th, ph, s, h=1e-4):
    """L_pm = pm e^{pm i ph} (d_th pm i cot th d_ph)  (3.18)"""
    dth = (f(th + h, ph) - f(th - h, ph)) / (2 * h)
    dph = (f(th, ph + h) - f(th, ph - h)) / (2 * h)
    return s * np.exp(s * 1j * ph) * (dth + s * 1j * dph / np.tan(th))


TH, PH = np.meshgrid(np.linspace(0.3, np.pi - 0.3, 25), np.linspace(0, 2 * np.pi, 31))
for l in range(4):
    for m in range(-l, l + 1):
        f = Y(l, m)
        v = f(TH, PH)
        assert np.allclose(Lz_fd(f, TH, PH), m * v, atol=1e-7)
        assert np.allclose(L2_fd(f, TH, PH), l * (l + 1) * v, atol=2e-5)
        # parity (3.23): Y(pi - th, ph + pi) = (-1)^l Y
        assert np.allclose(f(np.pi - TH, PH + np.pi), (-1)**l * v)
# Condon-Shortley table check
assert np.allclose(Y(1, 1)(TH, PH), -np.sqrt(3 / (8 * np.pi)) * np.sin(TH) * np.exp(1j * PH))
assert np.allclose(Y(2, 1)(TH, PH),
                   -np.sqrt(15 / (8 * np.pi)) * np.sin(TH) * np.cos(TH) * np.exp(1j * PH))
# Example 4.3: L_- Y_1^1 = sqrt(2) Y_1^0, and B3: L_- Y_2^2 = 2 Y_2^1
assert np.allclose(Lpm_fd(Y(1, 1), TH, PH, -1), np.sqrt(2) * Y(1, 0)(TH, PH), atol=1e-7)
assert np.allclose(Lpm_fd(Y(2, 2), TH, PH, -1), 2 * Y(2, 1)(TH, PH), atol=1e-7)
assert np.allclose(Lpm_fd(Y(2, 2), TH, PH, +1), 0, atol=1e-7)          # top rung
print("Part 2: L_z, L^2 eigenvalues, parity, L_- ladder (Ex 4.3) checked for l <= 3")

# orthonormality by Gauss-Legendre (in cos th) x uniform phi quadrature
xg, wg = np.polynomial.legendre.leggauss(12)
phg = np.linspace(0, 2 * np.pi, 16, endpoint=False)
THq, PHq = np.meshgrid(np.arccos(xg), phg, indexing="ij")
W = np.outer(wg, np.full(len(phg), 2 * np.pi / len(phg)))
lm = [(l, m) for l in range(4) for m in range(-l, l + 1)]
Ys = np.array([Y(l, m)(THq, PHq) for l, m in lm])
G = np.einsum("aij,bij,ij->ab", Ys.conj(), Ys, W)
assert np.allclose(G, np.eye(len(lm)), atol=1e-12)
print(f"  Gram matrix of {len(lm)} harmonics = I to {np.abs(G - np.eye(len(lm))).max():.1e}")

# SymPy: angular Laplacian from L_- L_+ + L_z^2 + L_z (3.19), and [Lx, Ly] = i Lz (Ex 4.1)
th, ph = sp.symbols("theta phi", real=True)
F = sp.Function("F")(th, ph)
Lp = lambda g: sp.exp(sp.I * ph) * (sp.diff(g, th) + sp.I * sp.cot(th) * sp.diff(g, ph))
Lm = lambda g: -sp.exp(-sp.I * ph) * (sp.diff(g, th) - sp.I * sp.cot(th) * sp.diff(g, ph))
Lzs = lambda g: -sp.I * sp.diff(g, ph)
L2s = Lm(Lp(F)) + Lzs(Lzs(F)) + Lzs(F)
lap = -(sp.diff(sp.sin(th) * sp.diff(F, th), th) / sp.sin(th) + sp.diff(F, ph, 2) / sp.sin(th)**2)
assert sp.simplify(sp.expand(L2s - lap)) == 0
x, y, z = sp.symbols("x y z", real=True)
g = sp.Function("g")(x, y, z)
Lx = lambda u: -sp.I * (y * sp.diff(u, z) - z * sp.diff(u, y))
Ly = lambda u: -sp.I * (z * sp.diff(u, x) - x * sp.diff(u, z))
Lz = lambda u: -sp.I * (x * sp.diff(u, y) - y * sp.diff(u, x))
assert sp.simplify(Lx(Ly(g)) - Ly(Lx(g)) - sp.I * Lz(g)) == 0
# C1: the would-be l = 1/2 harmonic sqrt(sin th) e^{i ph/2}
fh = sp.sqrt(sp.sin(th)) * sp.exp(sp.I * ph / 2)
assert sp.simplify(Lp(fh)) == 0
gh = sp.simplify(Lm(fh))
assert sp.simplify(gh + sp.cos(th) / sp.sqrt(sp.sin(th)) * sp.exp(-sp.I * ph / 2)) == 0
assert sp.simplify(Lm(gh) + sp.exp(-3 * sp.I * ph / 2) / sp.sin(th)**sp.Rational(3, 2)) == 0
print("  SymPy: (3.19), [Lx,Ly] = i Lz on g(x,y,z), and the l = 1/2 ladder failure (C1)")

fig = plt.figure(figsize=(11, 10))
tg, pg = np.meshgrid(np.linspace(0, np.pi, 60), np.linspace(0, 2 * np.pi, 80))
for l in range(4):
    for m in range(l + 1):
        r = np.abs(Y(l, m)(tg, pg))**2
        ax = fig.add_subplot(4, 4, 4 * l + m + 1, projection="3d")
        ax.plot_surface(r * np.sin(tg) * np.cos(pg), r * np.sin(tg) * np.sin(pg),
                        r * np.cos(tg), cmap="viridis", linewidth=0)
        ax.set_title(rf"$|Y_{l}^{m}|^2$", fontsize=9)
        ax.set_axis_off()
        lim = r.max()
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
fig.savefig(FIG / "spherical_harmonics.png", dpi=110)
plt.close(fig)

# %% [markdown]
# ## Part 3: Pauli algebra and n.sigma eigenvectors (§3.4-3.5, Example 4.4)

# %%
s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)
sig = [sx, sy, sz]
eps = np.zeros((3, 3, 3))
for a, b, cc in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    eps[a, b, cc], eps[b, a, cc] = 1, -1
for a in range(3):
    for b in range(3):
        rhs = (a == b) * s0 + 1j * sum(eps[a, b, k] * sig[k] for k in range(3))
        assert np.allclose(sig[a] @ sig[b], rhs)                        # (3.25)
    assert np.allclose(sig[a].conj().T, sig[a]) and np.isclose(np.trace(sig[a]), 0)
X2, Y2, Z2, _, _ = jmat(0.5)
assert np.allclose(2 * X2, sx) and np.allclose(2 * Y2, sy) and np.allclose(2 * Z2, sz)
print("Part 3: sigma_a sigma_b = delta_ab + i eps_abc sigma_c; matrices from (3.14) = Pauli")


def nvec(theta, phi):
    return np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def ndots(n):
    return n[0] * sx + n[1] * sy + n[2] * sz


def ket_n(theta, phi):
    """|n,+> = cos(th/2)|up> + e^{i ph} sin(th/2)|down>  (3.28)"""
    return np.array([np.cos(theta / 2), np.exp(1j * phi) * np.sin(theta / 2)])


def bloch(psi):
    return np.real([psi.conj() @ s @ psi for s in sig]) / np.vdot(psi, psi).real


th60, ph90 = np.pi / 3, np.pi / 2
kp = ket_n(th60, ph90)
km = np.array([np.sin(th60 / 2), -np.exp(1j * ph90) * np.cos(th60 / 2)])
n = nvec(th60, ph90)
assert np.allclose(ndots(n) @ kp, kp) and np.allclose(ndots(n) @ km, -km)
assert np.allclose(kp, [np.sqrt(3) / 2, 0.5j]) and np.allclose(km, [0.5, -np.sqrt(3) / 2 * 1j])
print(f"Ex 4.4: n = {np.round(n, 4)}, |n+> = {np.round(kp, 4)}, P(Sz=+) = {abs(kp[0])**2:.4f}, "
      f"P(Sz=-) = {abs(kp[1])**2:.4f}; |n-> = {np.round(km, 4)}")
assert np.isclose(abs(kp[0])**2, 0.75) and np.isclose(abs(km[0])**2, 0.25)
assert np.allclose(bloch(kp), n) and np.allclose(bloch(km), -n)        # antipodal (Thm 3.4)
for _ in range(20):                                                    # <S> = n/2 for all n
    t_, p_ = np.arccos(rng.uniform(-1, 1)), rng.uniform(0, 2 * np.pi)
    assert np.allclose(bloch(ket_n(t_, p_)), nvec(t_, p_))
# B4: chi = (3i, 4)/5
chi = np.array([3j, 4]) / 5
print(f"B4: Bloch vector of (3i,4)/5 = {bloch(chi)}, theta = {np.degrees(np.arccos(bloch(chi)[2])):.2f} deg")
assert np.allclose(bloch(chi), [0, -24 / 25, -7 / 25])

# %% [markdown]
# ## Part 4: Rotation operators and SU(2) -> SO(3) (§3.8)

# %%
def U_rot(theta, n):
    """e^{-i theta n.sigma/2} = cos(theta/2) I - i sin(theta/2) n.sigma  (3.36)"""
    return np.cos(theta / 2) * s0 - 1j * np.sin(theta / 2) * ndots(n)


def R_rod(theta, n):
    """Rodrigues: R = cos th I + (1 - cos th) n n^T + sin th [n]_x"""
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.cos(theta) * np.eye(3) + (1 - np.cos(theta)) * np.outer(n, n) + np.sin(theta) * K


for _ in range(10):
    nn = rng.normal(size=3); nn /= np.linalg.norm(nn)
    thr = rng.uniform(0, 4 * np.pi)
    U = U_rot(thr, nn)
    assert np.allclose(U, expm(-1j * thr * ndots(nn) / 2))
    R = R_rod(thr, nn)
    mv = rng.normal(size=3)
    assert np.allclose(U @ ndots(mv) @ U.conj().T, ndots(R @ mv))       # (3.37)
    Rtr = np.real([[np.trace(sig[a] @ U @ sig[b] @ U.conj().T) / 2 for b in range(3)]
                   for a in range(3)])
    assert np.allclose(Rtr, R) and np.isclose(np.linalg.det(Rtr), 1)    # C3
assert np.allclose(U_rot(2 * np.pi, nn), -s0) and np.allclose(U_rot(4 * np.pi, nn), s0)
assert np.allclose(U_rot(np.pi, [1, 0, 0]), -1j * sx)                   # B5
Hd = (sx + sz) / np.sqrt(2)
assert np.allclose(1j * U_rot(np.pi, np.array([1, 0, 1]) / np.sqrt(2)), Hd)
assert np.allclose(U_rot(np.pi / 2, [0, 1, 0]) @ [1, 0], np.array([1, 1]) / np.sqrt(2))
print("Part 4: closed form = expm; U(m.s)U^+ = (Rm).s; R(2pi) = -1, R(4pi) = +1; H = i R_n(pi)")

# %% [markdown]
# ## Part 5: Larmor precession (§3.6, Example 4.5), SI units then hbar = 1
# H = -gamma B S_z. Heisenberg: dS/dt = gamma S x B. Schrodinger: phase e^{-i w t} on |down>.

# %%
gam_p = PC["proton gyromag. ratio"][0]
f1 = gam_p * 1.0 / (2 * np.pi)
print(f"Ex 4.5: proton gamma = {gam_p:.6e} rad/s/T; f(1 T) = {f1/1e6:.3f} MHz, period "
      f"{1e9/f1:.2f} ns; 1.5 T: {1.5*f1/1e6:.2f} MHz; 3 T: {3*f1/1e6:.2f} MHz")
assert abs(f1 / 1e6 - 42.577) < 1e-3
t5 = 5e-9
sx_t5 = 0.5 * np.cos(gam_p * t5)
print(f"  <S_x>(5 ns)/hbar = {sx_t5:.4f}, <S_y>(5 ns)/hbar = {-0.5*np.sin(gam_p*t5):.4f}")
gam_e = PC["electron gyromag. ratio"][0]
print(f"A5: electron |gamma|/2pi at 1 T = {gam_e/2/np.pi/1e9:.2f} GHz; proton at 3 T = "
      f"{3*f1/1e6:.1f} MHz; ratio |gamma_e|/gamma_p = {gam_e/gam_p:.1f}")

# dimensionless: gamma B = w0 = 1, start tilted at theta = 60 deg, phi = 0
w0 = 1.0
psi0 = ket_n(np.pi / 3, 0.0)
H0 = -w0 * sz / 2
ts = np.linspace(0, 2 * np.pi, 201)
sch = np.array([bloch(expm(-1j * H0 * t) @ psi0) for t in ts])
hei = solve_ivp(lambda t, S: w0 * np.cross(S, [0, 0, 1.0]), (0, ts[-1]), nvec(np.pi / 3, 0),
                t_eval=ts, rtol=1e-11, atol=1e-12).y.T
ana = np.c_[np.sin(np.pi / 3) * np.cos(w0 * ts), -np.sin(np.pi / 3) * np.sin(w0 * ts),
            np.full_like(ts, 0.5)]
assert np.allclose(sch, ana, atol=1e-10) and np.allclose(hei, ana, atol=1e-8)
print("  Larmor: Schrodinger = Heisenberg = analytic cone (3.33) to 1e-8")

fig = plt.figure(figsize=(12, 7))
u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 30), np.linspace(0, np.pi, 15))
for k, tk in enumerate(np.linspace(0, 5 * np.pi / 3, 6)):
    ax = fig.add_subplot(2, 3, k + 1, projection="3d")
    ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color="0.85", lw=0.4)
    sel = ts <= tk + 1e-12
    ax.plot(*sch[sel].T, color="C0")
    ax.quiver(0, 0, 0, *sch[sel][-1], color="C3", lw=2)
    ax.quiver(0, 0, 0, 0, 0, 1.3, color="k", lw=1)
    ax.set_title(rf"$\omega_0 t = {tk:.2f}$", fontsize=9)
    ax.set_axis_off()
fig.suptitle("Larmor precession of the Bloch vector about B (black), gamma > 0")
fig.savefig(FIG / "larmor_bloch.png", dpi=110)
plt.close(fig)

# %% [markdown]
# ## Part 6: Rabi oscillations in the lab frame vs the Rabi formula (§3.7, Example 4.6)
# H(t) = -w0 Sz - w1 (cos wt Sx - sin wt Sy), w = w0 - Delta. P_down = (w1/W)^2 sin^2(W t/2).

# %%
def rabi_P(t, w1, D):
    W = np.hypot(w1, D)
    return (w1 / W)**2 * np.sin(W * t / 2)**2


def lab_P(D, w0=20.0, w1=1.0, T=4 * np.pi, n=161):
    w = w0 - D
    def rhs(t, y):
        Hm = -w0 * sz / 2 - w1 * (np.cos(w * t) * sx - np.sin(w * t) * sy) / 2
        return -1j * Hm @ y
    tt = np.linspace(0, T, n)
    sol = solve_ivp(rhs, (0, T), np.array([1, 0], complex), t_eval=tt, method="DOP853",
                    rtol=1e-11, atol=1e-12)
    return tt, np.abs(sol.y[1])**2


fig, ax = plt.subplots(figsize=(7, 4))
for D in (0.0, 0.1, 0.5, 1.0, 2.0):
    tt, P = lab_P(D)
    assert np.allclose(P, rabi_P(tt, 1.0, D), atol=1e-7), D
    ax.plot(tt, P, label=rf"$\Delta/\omega_1 = {D}$")
ax.set_xlabel(r"$\omega_1 t$"); ax.set_ylabel(r"$P_{\uparrow\to\downarrow}$"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "rabi.png", dpi=110); plt.close(fig)
P_res, P_det = rabi_P(np.pi, 1, 0), rabi_P(np.pi, 1, 0.1)
print(f"Ex 4.6: pi-pulse P = {P_res:.6f} at resonance; at Delta = 0.1 w1: P = {P_det:.6f}, "
      f"P_max = {1/1.01:.6f}, best time w1 t = {np.pi/np.sqrt(1.01):.4f}")
assert np.isclose(P_res, 1) and abs(P_det - 0.990038) < 1e-6
print(f"B6: Delta = w1: P_max = {rabi_P(np.pi/np.sqrt(2), 1, 1):.4f} at w1 t = {np.pi/np.sqrt(2):.4f}")
B1 = 10e-6
print(f"  MRI: B1 = 10 uT -> w1/2pi = {gam_p*B1/2/np.pi:.1f} Hz, pi-pulse {np.pi/(gam_p*B1)*1e3:.3f} ms")

# %% [markdown]
# ## Part 7: Adding angular momenta with np.kron (§3.10, Example 4.7)

# %%
def total(j1, j2):
    A, B = jmat(j1), jmat(j2)
    I1, I2 = np.eye(len(A[2])), np.eye(len(B[2]))
    return [np.kron(a, I2) + np.kron(I1, b) for a, b in zip(A, B)]


Jx, Jy, Jz, Jp, Jm = total(0.5, 0.5)
S2 = Jx @ Jx + Jy @ Jy + Jz @ Jz
w, V = np.linalg.eigh(S2)
assert np.allclose(w, [0, 2, 2, 2])
singlet = np.array([0, 1, -1, 0]) / np.sqrt(2)             # order: uu, ud, du, dd
assert np.isclose(abs(V[:, 0] @ singlet), 1)
uu = np.array([1, 0, 0, 0.0])
t0 = Jm @ uu / np.linalg.norm(Jm @ uu)
assert np.allclose(t0, [0, 1, 1, 0] / np.sqrt(2)) and np.allclose(Jm @ Jm @ uu / 2, [0, 0, 0, 1])
for A in (Jx, Jy, Jz):
    assert np.allclose(A @ singlet, 0)                      # B8
s1s2 = sum(np.kron(s, s) for s in sig)
assert np.isclose(singlet @ s1s2 @ singlet, -3) and np.isclose(uu @ s1s2 @ uu, 1)
print("Part 7: two spin-1/2: S^2 eigenvalues {0, 2, 2, 2}; singlet and triplet by lowering")


def multiplets_by_lowering(j1, j2):
    """Coupled states |j, m> as vectors in the product basis (m1 desc, m2 desc), built by
    lowering from the top of each multiplet; new tops orthogonal to previous states with
    Condon-Shortley sign <j1 j1; j2 (j - j1) | j j> > 0."""
    Jx, Jy, Jz, Jp, Jm = total(j1, j2)
    mz = np.real(np.diag(Jz))
    d2 = int(round(2 * j2)) + 1
    states = {}
    for jj in np.arange(j1 + j2, abs(j1 - j2) - 1e-9, -1):
        sub = np.where(np.isclose(mz, jj))[0]
        prev = [states[(jo, jj)] for jo in np.arange(j1 + j2, jj + 1e-9, -1) if (jo, jj) in states]
        P = np.eye(len(sub)) - sum(np.outer(p[sub], p[sub].conj()) for p in prev) if prev \
            else np.eye(len(sub))
        top = np.zeros(len(mz), complex)
        top[sub] = null_space(np.eye(len(sub)) - P)[:, 0] if prev else 1.0
        ref = [i for i in sub if i // d2 == 0]                 # component with m1 = j1
        if ref and abs(top[ref[0]]) > 1e-12:
            top *= np.sign(top[ref[0]].real)
        top /= np.linalg.norm(top)
        vec, mm = top, jj
        while mm >= -jj - 1e-9:
            states[(jj, mm)] = vec
            nv = Jm @ vec
            vec = nv / np.linalg.norm(nv) if np.linalg.norm(nv) > 1e-12 else nv
            mm -= 1
    return states, (Jx, Jy, Jz)


def basis_labels(j1, j2):
    return [(j1 - a, j2 - b) for a in range(int(round(2 * j1)) + 1)
            for b in range(int(round(2 * j2)) + 1)]


def check_cg(j1, j2):
    st, (Jx, Jy, Jz) = multiplets_by_lowering(j1, j2)
    J2 = Jx @ Jx + Jy @ Jy + Jz @ Jz
    labels = basis_labels(j1, j2)
    for (jj, mm), vec in st.items():
        assert np.allclose(J2 @ vec, jj * (jj + 1) * vec) and np.allclose(Jz @ vec, mm * vec)
        for k, (m1, m2) in enumerate(labels):
            ref = float(CG(sp.nsimplify(j1), sp.nsimplify(m1), sp.nsimplify(j2),
                           sp.nsimplify(m2), sp.nsimplify(jj), sp.nsimplify(mm)).doit())
            assert np.isclose(vec[k].real, ref, atol=1e-10), (jj, mm, m1, m2)
    return st


st = check_cg(1, 0.5)
dims = {jj: sum(1 for (a, _) in st if a == jj) for jj in (1.5, 0.5)}
assert dims == {1.5: 4, 0.5: 2} and 3 * 2 == sum(dims.values())
lab = basis_labels(1, 0.5)
v32, v12 = st[(1.5, 0.5)], st[(0.5, 0.5)]
print("Ex 4.7: 1 x 1/2 -> dims", dims, "; |3/2,1/2> =",
      {lab[k]: round(float(v32[k].real), 4) for k in np.nonzero(abs(v32) > 1e-9)[0]},
      "; |1/2,1/2> =", {lab[k]: round(float(v12[k].real), 4) for k in np.nonzero(abs(v12) > 1e-9)[0]})
assert np.isclose(v32[lab.index((0, 0.5))].real, np.sqrt(2 / 3))
assert np.isclose(v12[lab.index((1, -0.5))].real, np.sqrt(2 / 3))
v = st[(0.5, -0.5)]
print("B7/C2: |1/2,-1/2> =", {lab[k]: round(float(v[k].real), 4) for k in np.nonzero(abs(v) > 1e-9)[0]})
check_cg(0.5, 0.5); check_cg(1.5, 1)
for j1, j2 in [(0.5, 0.5), (1, 0.5), (1, 1), (2, 1), (1.5, 1.5)]:
    js = np.arange(abs(j1 - j2), j1 + j2 + 1e-9)
    assert (2 * j1 + 1) * (2 * j2 + 1) == sum(2 * js + 1)   # (3.42)
print("  CG coefficients agree with sympy for 1/2x1/2, 1x1/2, 3/2x1; counting rule checked")

# %% [markdown]
# ## Part 8: Stern-Gerlach cascade for any spin j (QM-1 Part 3 generalized)

# %%
def sg_measure(psi, n, j, rng):
    """Measure n.J on each row of psi; return outcome index (0 = m_n = +j) and collapsed states."""
    Jx, Jy, Jz, _, _ = jmat(j)
    w, U = np.linalg.eigh(n[0] * Jx + n[1] * Jy + n[2] * Jz)
    U = U[:, ::-1]
    P = np.abs(psi @ U.conj())**2
    out = (rng.random(len(psi))[:, None] > np.cumsum(P, axis=1)).sum(axis=1)
    out = np.minimum(out, len(w) - 1)
    return out, U[:, out].T


def binom_check(count, N, p, label):
    sd = np.sqrt(N * p * (1 - p)) if 0 < p < 1 else 0.5
    assert abs(count - N * p) < 4 * sd + 1e-9, label
    print(f"  {label}: {count}/{N} (predicted {N*p:.0f})")


N = 20_000
src = np.tile([1, 0, 0], (N, 1)).astype(complex)               # spin 1, m_z = +1
o, s = sg_measure(src, nvec(np.pi / 2, 0), 1, rng)
for k, p in enumerate([0.25, 0.5, 0.25]):
    binom_check((o == k).sum(), N, p, f"j=1: z+ then x, beam m_x={1-k}")
keep = o == 1
o2, _ = sg_measure(s[keep], nvec(0, 0), 1, rng)
for k, p in enumerate([0.5, 0.0, 0.5]):
    binom_check((o2 == k).sum(), keep.sum(), p, f"j=1: x0 then z, beam m_z={1-k}")
src = np.zeros((N, 4), complex); src[:, 0] = 1                  # spin 3/2 top state
o, _ = sg_measure(src, nvec(np.pi / 3, 0), 1.5, rng)
pb = [0.75**3, 3 * 0.75**2 * 0.25, 3 * 0.75 * 0.25**2, 0.25**3]
for k, p in enumerate(pb):
    binom_check((o == k).sum(), N, p, f"j=3/2: z top then 60deg, beam {k}")
print(f"  j=3/2 binomial probabilities {np.round(pb, 4)} = 27/64, 27/64, 9/64, 1/64")

# %% [markdown]
# ## Part 9: Numbers quoted in the notes and solutions

# %%
r_e = PC["classical electron radius"][0]
v_surf = 1.25 * HBAR / (m_e * r_e)
print(f"§3.4: solid sphere, radius r_e = {r_e:.3e} m, S = hbar/2 -> surface speed "
      f"{v_surf/C_LIGHT:.0f} c")
assert 170 < v_surf / C_LIGHT < 173
# B2: <Lx^2> in |l,m> = (l(l+1) - m^2)/2
Lx1, Ly1, Lz1, _, _ = jmat(2)
for k, mm in enumerate([2, 1, 0, -1, -2]):
    assert np.isclose(Lx1[:, k] @ Lx1[:, k].conj(), (6 - mm**2) / 2)
# B7: spin-1 S_x eigenvectors
w, U = np.linalg.eigh(Sx1)
assert np.allclose(w, [-1, 0, 1])
assert np.isclose(abs(U[:, 2] @ np.array([1, np.sqrt(2), 1]) / 2), 1)
assert np.isclose(abs(U[:, 1] @ np.array([1, 0, -1]) / np.sqrt(2)), 1)
# Example 4.5 extras: rotation angle at 5 ns, MRI polarization, protons in water
print(f"Ex 4.5: w0 t(5 ns) = {gam_p*t5:.4f} rad = {np.degrees(gam_p*t5):.1f} deg; "
      f"hbar w0 = {HBAR*gam_p/PC['elementary charge'][0]:.3e} eV; "
      f"polarization hbar w0 / 2 kT (310 K) = {HBAR*gam_p/(2*PC['Boltzmann constant'][0]*310):.2e}; "
      f"kT = {PC['Boltzmann constant'][0]*310/PC['elementary charge'][0]:.4f} eV; "
      f"water protons/cm^3 = {2*PC['Avogadro constant'][0]/18.015:.2e}")
# Example 4.6 intermediate values
Om = np.sqrt(1.01)
print(f"Ex 4.6: Omega/w1 = {Om:.6f}, Omega t_pi/2 = {Om*np.pi/2:.6f}, sin^2 = {np.sin(Om*np.pi/2)**2:.6f}")
# §3.11: <L.S>/hbar^2 = (j(j+1) - l(l+1) - 3/4)/2 for l = 1
LS = {jj: (jj * (jj + 1) - 2 - 0.75) / 2 for jj in (1.5, 0.5)}
st_ls, (Jx_, Jy_, Jz_) = multiplets_by_lowering(1, 0.5)
Lv, Sv = jmat(1), jmat(0.5)
LdotS = sum(np.kron(Lv[i], Sv[i]) for i in range(3))
for jj in (1.5, 0.5):
    vec = st_ls[(jj, 0.5)]
    assert np.isclose((vec.conj() @ LdotS @ vec).real, LS[jj])
print(f"§3.11: <L.S>/hbar^2 = {LS[1.5]} (j=3/2), {LS[0.5]} (j=1/2)")
# A2: l = 1 Lx eigenvalues; A3: (a.sigma)^2 for a = (1,2,2)
assert np.allclose(np.linalg.eigvalsh(Sx1), [-1, 0, 1])
asig = ndots(np.array([1, 2, 2.0]))
assert np.allclose(asig @ asig, 9 * s0) and np.allclose(np.linalg.eigvalsh(asig), [-3, 3])
assert np.allclose(sx @ sy @ sz, 1j * s0)
# A4: |+y>: P(Sx = +) = 1/2, P(Sz = +) = 1/2, Bloch vector (0,1,0)
assert np.allclose(bloch(np.array([1, 1j]) / np.sqrt(2)), [0, 1, 0])
# B2: minimum cone angle for l = 2, m = 2
print(f"B2: l=2, m=2 cone half-angle = {np.degrees(np.arccos(2/np.sqrt(6))):.2f} deg")
# B4: P(Sx=+) and P(Sy=+) for chi = (3i,4)/5
pxp = abs(np.array([1, 1]) / np.sqrt(2) @ chi)**2
pyp = abs(np.array([1, 1j]).conj() / np.sqrt(2) @ chi)**2
print(f"B4: P(Sx=+) = {pxp:.4f}, P(Sy=+) = {pyp:.4f}, phi = {np.degrees(np.arctan2(-24, 0)):.1f} deg")
assert np.isclose(pxp, 0.5) and np.isclose(pyp, 0.02)
for _ in range(10):                                  # B4(c): |<n'|n>|^2 = (1 + n.n')/2
    a1, b1, a2, b2 = rng.uniform(0, np.pi), rng.uniform(0, 6.3), rng.uniform(0, np.pi), rng.uniform(0, 6.3)
    assert np.isclose(abs(np.vdot(ket_n(a1, b1), ket_n(a2, b2)))**2,
                      (1 + nvec(a1, b1) @ nvec(a2, b2)) / 2)
# B5: H Z H = X
assert np.allclose(Hd @ sz @ Hd, sx)
# B6: line-shape: time average of P is (1/2) w1^2/(w1^2 + D^2)
for D in (0.0, 0.5, 1.0, 2.0):
    tt = np.linspace(0, 400 * np.pi, 200001)
    assert abs(rabi_P(tt, 1, D).mean() - 0.5 / (1 + D**2)) < 2e-3
print(f"B6: detuning that halves P_max at B1 = 10 uT: {gam_p*B1/2/np.pi:.1f} Hz")
# B8: sigma1.sigma2 = 2 S^2 - 3; singlet in the x basis
assert np.allclose(s1s2, 2 * S2 - 3 * np.eye(4))
px_, mx_ = np.array([1, 1]) / np.sqrt(2), np.array([1, -1]) / np.sqrt(2)
sing_x = (np.kron(px_, mx_) - np.kron(mx_, px_)) / np.sqrt(2)
assert np.isclose(abs(sing_x @ singlet), 1)
print(f"B8: singlet in x basis = {np.round(sing_x, 4)} (= -singlet)")
# C3: R for U = R_y(pi/2) is the 90 degree rotation about y
Uy = U_rot(np.pi / 2, [0, 1, 0])
Ry = np.real([[np.trace(sig[a] @ Uy @ sig[b] @ Uy.conj().T) / 2 for b in range(3)] for a in range(3)])
assert np.allclose(Ry, [[0, 0, 1], [0, 1, 0], [-1, 0, 0]])

# %% [markdown]
# ## Your turn
# 1. Spin-1 Larmor precession: evolve |1,1>_x under H = -w0 Sz with jmat(1) and plot <S> and
#    the quadrupole components <Sx^2 - Sy^2>, <{Sx,Sy}>. At what frequency do they rotate?
# 2. Clebsch-Gordan for 1 x 1 by hand (lowering from |2,2>), then call check_cg(1, 1) and read
#    your nine coupled states off multiplets_by_lowering; compare with sympy's CG.
# 3. The line shape: plot the time-averaged flip probability vs Delta/w1 from lab_P and show it
#    is the Lorentzian (1/2) w1^2 / (w1^2 + Delta^2) (F3's Lorentzian, FWHM 2 w1).

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T0:.1f} s.")
