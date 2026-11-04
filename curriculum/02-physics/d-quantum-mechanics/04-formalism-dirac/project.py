# %% [markdown]
# # QM-4 Project: An "Abstract QM" Library and Everything It Can Check
#
# We build two small classes, `State` and `Operator`, over a finite orthonormal
# basis: inner products, adjoints, Hermiticity checks, commutators, Born
# probabilities with collapse (Lüders rule for degenerate eigenvalues), and
# unitary evolution by scipy.linalg.expm. Then we use them on (a) the truncated
# oscillator in the number basis, (b) a position grid and its DFT momentum basis,
# (c) the two-state system and Bloch precession, (d) the uncertainty relation,
# (e) the generalized Ehrenfest theorem, (f) simulated measurements, (g) the
# translation operator, (h) parity and degeneracy. Part 9 checks every number in
# the worked examples and solutions. Units: hbar = m = omega = 1 unless stated.
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.constants import physical_constants as PC, h as H_SI, eV as EV
from scipy.linalg import expm
from scipy.optimize import minimize

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.time()
HBAR = 1.0
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)

# %% [markdown]
# ## Part 0: the library (notes §3.2-3.10)
# A ket is a column of components c_n = <n|psi>; a bra is its conjugate row.
# An operator is its matrix Q_mn = <m|Q|n>. Everything else is LA-8.


# %%
class State:
    """A ket |psi> stored as its components <n|psi> in a fixed orthonormal basis."""

    def __init__(self, vec, normalize=True):
        v = np.asarray(vec, dtype=complex).ravel()
        self.vec = v / np.linalg.norm(v) if normalize else v

    def bra(self):
        """<psi| as a row: the Riesz functional, conjugate transpose (LA-12 Thm 3.9)."""
        return self.vec.conj()

    def inner(self, other):
        """<self|other>, conjugate-linear in the first slot."""
        return self.bra() @ other.vec

    def norm(self):
        return np.sqrt(self.inner(self).real)


class Operator:
    """A linear operator stored as its matrix <m|Q|n>."""

    def __init__(self, mat):
        self.mat = np.asarray(mat, dtype=complex)

    def __matmul__(self, other):
        if isinstance(other, State):
            return State(self.mat @ other.vec, normalize=False)
        return Operator(self.mat @ other.mat)

    def __add__(self, other):
        return Operator(self.mat + other.mat)

    def __sub__(self, other):
        return Operator(self.mat - other.mat)

    def __rmul__(self, c):
        return Operator(c * self.mat)

    def dag(self):
        """Adjoint: <phi|Q psi> = <Q^dag phi|psi>; matrix = conjugate transpose (LA-8 Thm 3.8)."""
        return Operator(self.mat.conj().T)

    def is_hermitian(self, tol=1e-12):
        return np.max(np.abs(self.mat - self.mat.conj().T)) < tol

    def expect(self, psi):
        """<psi|Q|psi> (3.13)."""
        return psi.inner(self @ psi)

    def sigma(self, psi):
        """sigma_Q = ||(Q - <Q>) psi|| (Theorem 3.5)."""
        q = self.expect(psi)
        return np.linalg.norm(self.mat @ psi.vec - q * psi.vec)

    def born(self, psi, tol=1e-9):
        """Distinct eigenvalues q, probabilities <psi|P_q|psi>, projectors P_q (3.11)."""
        assert self.is_hermitian()
        w, V = np.linalg.eigh(self.mat)
        groups = []
        for i, wi in enumerate(w):
            if groups and abs(wi - groups[-1][0]) < tol:
                groups[-1][1].append(i)
            else:
                groups.append([wi, [i]])
        qs, probs, projs = [], [], []
        for q, idx in groups:
            P = V[:, idx] @ V[:, idx].conj().T
            qs.append(q); projs.append(P)
            probs.append(np.real(psi.vec.conj() @ P @ psi.vec))
        return np.array(qs), np.array(probs), projs

    def measure(self, psi, rng):
        """One measurement: outcome q_k with Born probability, then collapse P_k psi/||P_k psi|| (3.12)."""
        qs, probs, projs = self.born(psi)
        k = rng.choice(len(qs), p=probs / probs.sum())
        return qs[k], State(projs[k] @ psi.vec)

    def U(self, t):
        """Propagator exp(-i H t / hbar) (3.22)."""
        return Operator(expm(-1j * self.mat * t / HBAR))


def comm(A, B):
    """[A, B] = AB - BA."""
    return A @ B - B @ A


def close(a, b, tol, label=""):
    err = np.max(np.abs(np.asarray(a) - np.asarray(b)))
    assert err < tol, f"{label}: {a} vs {b} (err {err:.2e})"
    return err


# %% [markdown]
# ## Part 1: Examples 4.1 and 4.3 (box in Dirac notation; a 3x3 observable)

# %%
x_, L_ = sp.symbols("x L", positive=True)
n_ = sp.symbols("n", integer=True, positive=True)
psi_box = sp.sqrt(30 / L_**5) * x_ * (L_ - x_)
cn = sp.simplify(sp.integrate(sp.sqrt(2 / L_) * sp.sin(n_ * sp.pi * x_ / L_) * psi_box, (x_, 0, L_)))
c1 = sp.simplify(cn.subs(n_, 1)); c2 = sp.simplify(cn.subs(n_, 2))
assert sp.simplify(c1 - 8 * sp.sqrt(15) / sp.pi**3) == 0 and c2 == 0
H_direct = sp.integrate(psi_box * (-sp.Rational(1, 2)) * sp.diff(psi_box, x_, 2), (x_, 0, L_))
assert sp.simplify(H_direct - 5 / L_**2) == 0                     # hbar = m = 1: 5 hbar^2/(m L^2)
odd = np.arange(1, 4000, 2)
pn = 960 / (np.pi**6 * odd**6)
assert abs(pn.sum() - 1) < 1e-14
EH_over_E1 = np.sum(pn * odd**2)
close(EH_over_E1, 10 / np.pi**2, 1e-8, "<H>/E1")
print(f"Ex 4.1: |c1|^2 = {pn[0]:.6f}, |c3|^2 = {pn[1]:.3e}, |c5|^2 = {pn[2]:.3e}, "
      f"<H>/E1 = {EH_over_E1:.5f} = 10/pi^2")

Q3 = Operator([[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]])
psi3 = State([1, 1j, 0])
qs, probs, projs = Q3.born(psi3)
close(qs, [-np.sqrt(2), 0, np.sqrt(2)], 1e-12, "eigs")
close(probs, [(3 - 2*np.sqrt(2))/8, 0.25, (3 + 2*np.sqrt(2))/8], 1e-12, "Born")
close(np.sum(qs * probs), 1.0, 1e-12, "<Q> from Born"); close(Q3.expect(psi3), 1.0, 1e-12, "<Q> direct")
close(Q3.sigma(psi3)**2, 0.5, 1e-12, "var Q")
plus = State(projs[2] @ psi3.vec)
target = State([-1j, np.sqrt(2), 1j])
close(abs(target.inner(plus)), 1.0, 1e-12, "collapse")
print(f"Ex 4.3: eigenvalues {np.round(qs, 5)}, probabilities {np.round(probs, 5)}, <Q> = 1, sigma_Q^2 = 0.5")

# %% [markdown]
# ## Part 2 (a): the truncated oscillator, N = 60 (Examples 4.2, 4.6)
# a|n> = sqrt(n)|n-1>; x = (a + a^dag)/sqrt2, p = i(a^dag - a)/sqrt2. In the
# truncated space [a, a^dag] = 1 - N|N-1><N-1|, so [x, p] = i hbar except in the
# last diagonal entry, which is -i hbar (N - 1): the trace of a commutator is 0.

# %%
def ladder(N):
    return Operator(np.diag(np.sqrt(np.arange(1, N)), 1))


for N in (4, 60):
    a = ladder(N); ad = a.dag()
    X = (1 / np.sqrt(2)) * (a + ad); P = (1j / np.sqrt(2)) * (ad - a)
    C = comm(X, P).mat
    expected = 1j * np.eye(N); expected[-1, -1] = -1j * (N - 1)
    close(C, expected, 1e-12, "[x,p] truncated")
    assert X.is_hermitian() and P.is_hermitian() and abs(np.trace(C)) < 1e-10
    if N == 4:
        print("Ex 4.2: x matrix (N=4) * sqrt2 =\n", np.round(X.mat.real * np.sqrt(2), 4))
        print("        diag of [x,p]/i =", np.round(np.diag(C).imag, 4))
print(f"Part a: N=60, [x,p] = i in the bulk, edge entry = {C[-1, -1]:.0f}")

Hosc = Operator(np.diag(np.arange(N) + 0.5))
close(np.linalg.eigvalsh(Hosc.mat), np.arange(N) + 0.5, 1e-12, "spectrum")
Hxp = (0.5) * (P @ P + X @ X)                       # p^2/2 + x^2/2 built from truncated x, p
close(Hxp.mat[:N-1, :N-1], Hosc.mat[:N-1, :N-1], 1e-12, "H bulk")
print(f"Part a: (p^2+x^2)/2 matches (N + 1/2) except the last entry: {Hxp.mat[-1,-1].real:.2f} vs {N-0.5}")
for t in (0.3, 1.0, 2.7):
    U = Hosc.U(t)
    XH = U.dag() @ X @ U
    close(XH.mat, np.cos(t) * X.mat + np.sin(t) * P.mat, 1e-10, "Heisenberg x_H(t)")

def coherent(alpha, N):
    from math import factorial
    c = np.array([alpha**k / np.sqrt(float(factorial(k))) for k in range(N)], complex)
    return State(c * np.exp(-abs(alpha)**2 / 2), normalize=False)

alpha = 1 + 1j
coh = coherent(alpha, N)
close(coh.norm(), 1.0, 1e-12, "coherent norm")
close((a @ coh).vec, alpha * coh.vec, 1e-9, "a|alpha> = alpha|alpha>")
ts = np.linspace(0, 2*np.pi, 9)
xt = [(Hosc.U(t).dag() @ X @ Hosc.U(t)).expect(coh).real for t in ts]
close(xt, 2 * np.cos(ts - np.pi/4), 1e-10, "<x(t)> coherent")
close(X.sigma(coh), np.sqrt(0.5), 1e-10, "sigma_x coherent")
print(f"Ex 4.6: alpha = 1+i, <x>(0) = {xt[0]:.5f}, <x(t)> = 2 cos(t - pi/4), sigma_x = {X.sigma(coh):.5f}")

# %% [markdown]
# ## Part 3 (b): position grid and the DFT as the change of basis to momentum (§3.4)
# On a periodic grid x_j, the unitary DFT matrix F_kj = exp(-i k x_j)/sqrt(Ng) has
# rows <p_k|. Phi = F psi; Parseval is unitarity. Any translation-invariant
# (circulant) operator, such as the central-difference p, is diagonal in this basis.

# %%
Ng, Lbox = 256, 40.0
dx = Lbox / Ng
xg = (np.arange(Ng) - Ng // 2) * dx
kg = 2 * np.pi * np.fft.fftfreq(Ng, d=dx)
F = np.exp(-1j * np.outer(kg, xg)) / np.sqrt(Ng)
close(F.conj().T @ F, np.eye(Ng), 1e-10, "F unitary")
psi_g = State(np.exp(-(xg - 1.0)**2 / 4 + 2j * xg))
Phi = F @ psi_g.vec
close(np.sum(abs(Phi)**2), 1.0, 1e-12, "Parseval")
Phi_cont = (2 / np.pi)**0.25 * np.exp(-(kg - 2)**2 - 1j * (kg - 2) * 1.0) * np.sqrt(2 * np.pi / Lbox)
close(abs(Phi), abs(Phi_cont), 1e-10, "DFT vs continuous FT")
D = (np.roll(np.eye(Ng), 1, axis=1) - np.roll(np.eye(Ng), -1, axis=1)) / (2 * dx)
Pfd = Operator(-1j * HBAR * D)
Pmom = F @ Pfd.mat @ F.conj().T
assert Pfd.is_hermitian()
close(Pmom, np.diag(HBAR * np.sin(kg * dx) / dx), 1e-10, "p diagonal in momentum basis")
Pspec = Operator(F.conj().T @ np.diag(HBAR * kg) @ F)
print(f"Part b: F unitary, Parseval, finite-difference p diagonal in F-basis with eigenvalues sin(k dx)/dx; "
      f"<p> = {Pspec.expect(psi_g).real:.6f}")
close(Pspec.expect(psi_g), 2.0, 1e-8, "<p>")

# %% [markdown]
# ## Part 4 (c): the two-state system H = eps sigma_z + Delta sigma_x (Example 4.5)
# Bloch vector n = <sigma> obeys dn/dt = (2/hbar) h x n (3.29): rotation about h at rate 2|h|/hbar.

# %%
eps, Dl = 0.3, 0.4
H2 = Operator(eps * SZ + Dl * SX)
hvec = np.array([Dl, 0, eps]); hmag = np.linalg.norm(hvec); hhat = hvec / hmag
close(np.linalg.eigvalsh(H2.mat), [-0.5, 0.5], 1e-12, "E+-")
up = State([1, 0])

def bloch(psi):
    return np.array([Operator(S).expect(psi).real for S in (SX, SY, SZ)])

def rotate(v, axis, ang):
    """Rodrigues rotation of v about unit axis by ang."""
    return v*np.cos(ang) + np.cross(axis, v)*np.sin(ang) + axis*np.dot(axis, v)*(1 - np.cos(ang))

tt = np.linspace(0, 2*np.pi, 200)
traj = np.array([bloch(H2.U(t) @ up) for t in tt])
close(traj, [rotate(np.array([0, 0, 1.]), hhat, 2*hmag*t/HBAR) for t in tt], 1e-10, "precession")
Pflip = 1 - np.array([abs(up.inner(H2.U(t) @ up))**2 for t in tt])
rabi = lambda t, e, d: d**2/(e**2+d**2) * np.sin(np.sqrt(e**2+d**2)*t/HBAR)**2
close(Pflip, rabi(tt, eps, Dl), 1e-12, "Rabi formula")
P1 = 1 - abs(up.inner(H2.U(1.0) @ up))**2
close(P1, 0.147103, 1e-6, "P(1)")
print(f"Ex 4.5: E = +-{hmag:.2f}, Pmax = {Dl**2/hmag**2:.2f} at t = pi, P(t=1) = {P1:.6f}, "
      f"Bloch period 2pi/(2|h|) = {np.pi/hmag:.4f}")
eps_scan = np.linspace(-1.5, 1.5, 61)
Pmax_scan = [max(1 - abs(up.inner(Operator(e*SZ + Dl*SX).U(t) @ up))**2
                 for t in np.linspace(0, 20, 801)) for e in eps_scan]
close(Pmax_scan, Dl**2/(eps_scan**2 + Dl**2), 2e-3, "Lorentzian")

fig = plt.figure(figsize=(11, 3.8))
ax = fig.add_subplot(131, projection="3d")
u_, v_ = np.mgrid[0:2*np.pi:30j, 0:np.pi:15j]
ax.plot_wireframe(np.cos(u_)*np.sin(v_), np.sin(u_)*np.sin(v_), np.cos(v_), color="0.85", lw=0.5)
ax.plot(*traj.T, color="C0", lw=2); ax.quiver(0, 0, 0, *hhat, color="C3")
ax.set_title("Bloch vector precesses about h")
ax2 = fig.add_subplot(132); ax2.plot(tt, Pflip, label="numerics"); ax2.plot(tt, rabi(tt, eps, Dl), "--", label="(3.30)")
ax2.set_xlabel("t"); ax2.set_ylabel("P(up -> down)"); ax2.legend()
ax3 = fig.add_subplot(133); ax3.plot(eps_scan, Pmax_scan, "o", ms=3); ax3.plot(eps_scan, Dl**2/(eps_scan**2+Dl**2))
ax3.set_xlabel("detuning eps"); ax3.set_ylabel("max flip probability")
fig.tight_layout(); fig.savefig(FIG / "two_state.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 5 (d): the uncertainty relation (Theorem 3.7) on random Hermitian A, B
# Saturation: eigenvectors of A - i lam B (lam real) satisfy (A - <A>)psi = i lam (B - <B>)psi
# (Exercise B3), so they saturate it exactly. Optimization finds the same floor.

# %%
def rand_herm(n):
    M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return Operator((M + M.conj().T) / 2)

def rand_state(n):
    return State(rng.normal(size=n) + 1j * rng.normal(size=n))

nd = 4
A, B = rand_herm(nd), rand_herm(nd)
Cab = comm(A, B)
ratios = []
for _ in range(2000):
    s = rand_state(nd)
    lhs, rhs = A.sigma(s) * B.sigma(s), abs(Cab.expect(s)) / 2
    assert lhs >= rhs - 1e-12
    ratios.append(lhs / rhs)

def ratio(vr):
    s = State(vr[:nd] + 1j * vr[nd:])
    return A.sigma(s) * B.sigma(s) / (abs(Cab.expect(s)) / 2 + 1e-14)

best = min((minimize(ratio, rng.normal(size=2*nd), method="BFGS") for _ in range(5)), key=lambda r: r.fun)
w, V = np.linalg.eig(A.mat - 1j * 0.7 * B.mat)
sat = State(V[:, 0])
sat_ratio = A.sigma(sat) * B.sigma(sat) / (abs(Cab.expect(sat)) / 2)
close(sat_ratio, 1.0, 1e-9, "eigvec of A - i lam B saturates")
assert best.fun < 1.001
print(f"Part d: 2000 random states obey sigma_A sigma_B >= |<[A,B]>|/2 (min ratio {min(ratios):.3f}); "
      f"BFGS floor {best.fun:.6f}; eigenvector of A - 0.7iB gives {sat_ratio:.12f}")

# Example 4.4: spin-1/2 at theta = pi/3, phi = pi/4 (hbar = 1, S = sigma/2)
th, ph = np.pi/3, np.pi/4
spin = State([np.cos(th/2), np.exp(1j*ph)*np.sin(th/2)])
Sx, Sy, Sz = (Operator(0.5*S) for S in (SX, SY, SZ))
prod = Sx.sigma(spin) * Sy.sigma(spin); bound = 0.5 * abs(Sz.expect(spin))
close(comm(Sx, Sy).mat, 1j * Sz.mat, 1e-14, "[Sx,Sy]")
close([prod, bound], [0.15625, 0.125], 1e-12, "Ex 4.4")
print(f"Ex 4.4: sigma_Sx sigma_Sy = {prod:.5f} hbar^2 >= (hbar/2)|<Sz>| = {bound:.5f} hbar^2")

# %% [markdown]
# ## Part 6 (e): generalized Ehrenfest d<Q>/dt = (i/hbar)<[H, Q]> (3.20)

# %%
def ehrenfest_check(Hop, Q, psi0, t, h=1e-5):
    d_num = (Q.expect(Hop.U(t + h) @ psi0) - Q.expect(Hop.U(t - h) @ psi0)).real / (2*h)
    d_th = ((1j / HBAR) * comm(Hop, Q).expect(Hop.U(t) @ psi0)).real
    return d_num, d_th

for Hop, Q, psi0, lab in [(Hosc, X, coh, "oscillator x"), (Hosc, X @ X, coh, "oscillator x^2"),
                          (H2, Operator(SX), up, "two-state sigma_x"), (H2, Operator(SZ), up, "two-state sigma_z")]:
    dn, dt_ = ehrenfest_check(Hop, Q, psi0, 0.8)
    close(dn, dt_, 1e-7, lab)
print("Part e: Ehrenfest holds for x, x^2 (oscillator) and sigma_x, sigma_z (two-state)")

# %% [markdown]
# ## Part 7 (f): 10^4 measurements of the Example 4.3 observable

# %%
nshots = 10_000
outs = np.array([Q3.measure(psi3, rng)[0] for _ in range(nshots)])
freq = np.array([np.mean(np.isclose(outs, q)) for q in qs])
se = np.sqrt(probs * (1 - probs) / nshots)
assert np.all(np.abs(freq - probs) < 4 * se)
q_, s_ = Q3.measure(psi3, rng)
q2_, s2_ = Q3.measure(s_, rng)
assert q2_ == q_                          # repeat measurement gives the same answer
print(f"Part f: frequencies {np.round(freq, 4)} vs Born {np.round(probs, 4)}; sample mean {outs.mean():.4f} vs <Q> = 1")
fig, ax = plt.subplots(figsize=(5, 3.2))
ax.bar(np.arange(3) - 0.2, freq, 0.4, label="10^4 shots"); ax.bar(np.arange(3) + 0.2, probs, 0.4, label="Born")
ax.set_xticks(range(3), ["-sqrt2", "0", "sqrt2"]); ax.legend(); fig.tight_layout()
fig.savefig(FIG / "measurements.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 8 (g)-(h): translation operator and symmetry
# T(a) = exp(-i p a / hbar) shifts psi(x) -> psi(x - a) (3.31); <psi|T(a)|psi> = exp(-a^2/8s^2).
# Parity commutes with H for even V; a ring has translation AND parity, which do not commute,
# forcing the +-k degeneracy (Theorem 3.15).

# %%
s0, a_sh = 1.0, 2.0
g0 = State((2*np.pi*s0**2)**-0.25 * np.exp(-xg**2 / (4*s0**2)) * np.sqrt(dx), normalize=False)
Ta = Operator(expm(-1j * Pspec.mat * a_sh / HBAR))
close(Ta.mat @ Ta.mat.conj().T, np.eye(Ng), 1e-10, "T unitary")
shifted = Ta @ g0
close(shifted.vec, (2*np.pi*s0**2)**-0.25 * np.exp(-(xg - a_sh)**2/(4*s0**2)) * np.sqrt(dx), 1e-10, "shift")
ov = g0.inner(Ta @ g0).real
close(ov, np.exp(-a_sh**2 / (8*s0**2)), 1e-10, "overlap")
Xg = Operator(np.diag(xg))
close(Xg.expect(shifted).real, a_sh, 1e-10, "<x> shifted")
fd_err = np.max(abs(expm(-1j * Pfd.mat * a_sh) @ g0.vec - shifted.vec))
print(f"Ex 4.7: <psi|T(2)|psi> = {ov:.5f} = e^(-1/2); spectral p shifts exactly, "
      f"finite-difference p gives error {fd_err:.2e} (dispersion)")

# parity on a symmetric double well, V = (x^2 - 1)^2 * 4 on x in [-4, 4]
Nd = 401
xd = np.linspace(-4, 4, Nd); hd = xd[1] - xd[0]
lap = (np.diag(np.ones(Nd-1), 1) + np.diag(np.ones(Nd-1), -1) - 2*np.eye(Nd)) / hd**2
Hdw = Operator(-0.5 * lap + np.diag(4 * (xd**2 - 1)**2))
Par = Operator(np.eye(Nd)[::-1])
close(comm(Par, Hdw).mat, 0, 1e-10, "[P,H]=0")
Ed, Vd = np.linalg.eigh(Hdw.mat)
parities = [Operator(Par.mat).expect(State(Vd[:, k])).real for k in range(6)]
close(parities, [1, -1, 1, -1, 1, -1], 1e-8, "parity labels")
print(f"Part h: double-well levels {np.round(Ed[:4], 4)}, parities {np.round(parities[:4]).astype(int)}, "
      f"tunnel splitting {Ed[1]-Ed[0]:.4e}")

Nr = 64
dr = 2*np.pi / Nr
Tring = Operator(np.roll(np.eye(Nr), 1, axis=0))          # (T psi)_j = psi_{j-1}
lap_r = (Tring.mat + Tring.mat.T - 2*np.eye(Nr)) / dr**2
Hring = Operator(-0.5 * lap_r)
Pring = Operator(np.eye(Nr)[(-np.arange(Nr)) % Nr])      # j -> -j
assert np.allclose(comm(Tring, Hring).mat, 0) and np.allclose(comm(Pring, Hring).mat, 0)
assert not np.allclose(comm(Tring, Pring).mat, 0)
Er = np.linalg.eigvalsh(Hring.mat)
close(Er[1], Er[2], 1e-10, "ring degeneracy"); close(Er[3], Er[4], 1e-10, "ring degeneracy 2")
print(f"Part h: ring levels {np.round(Er[:5], 5)}: pairs degenerate because T and P commute with H but not each other")

# %% [markdown]
# ## Part 9: numbers quoted in solutions.md

# %%
phi_a, psi_a = State([1j, 1]), State([1, 2j])
close(phi_a.inner(psi_a), 1j / np.sqrt(10), 1e-14, "A1")
outer = np.outer(psi_a.vec, phi_a.bra())
close(np.trace(outer), phi_a.inner(psi_a), 1e-14, "A1: Tr |psi><phi| = <phi|psi>")
close(outer.conj().T, np.outer(phi_a.vec, psi_a.bra()), 1e-14, "A1 adjoint of outer product")
Qd = Operator([[2, 1j, 0], [-1j, 2, 0], [0, 0, 3]])
qd, pdg, prj = Qd.born(State([1, 1, 1]))
close(qd, [1, 3], 1e-12, "A5 eigs"); close(pdg, [1/3, 2/3], 1e-12, "A5 probs")
post = State(prj[1] @ State([1, 1, 1]).vec)
close(post.vec, np.array([(1+1j)/2, (1-1j)/2, 1]) / np.sqrt(2), 1e-12, "A5 post state")
close([Qd.expect(State([1, 1, 1])).real, Qd.sigma(State([1, 1, 1]))], [7/3, np.sqrt(8)/3], 1e-12, "A5")
print(f"A5: probabilities {pdg}, post-measurement {np.round(post.vec, 4)}, <Q> = 7/3, sigma = {np.sqrt(8)/3:.4f}")

xs, ps = sp.symbols("x p", real=True); fsym = sp.Function("f")(xs)
pop = lambda g: -sp.I * sp.diff(g, xs)
p3f = pop(pop(pop(fsym)))
assert sp.simplify(xs * p3f - pop(pop(pop(xs * fsym))) - sp.I * 3 * pop(pop(fsym))) == 0
print("A4: [x, p^3] = 3 i hbar p^2 (SymPy)")

# notes numbers: Example 4.4 Bloch vector, Example 4.5 eigenvector, MRI frequency
close(bloch(spin), [np.sqrt(6)/4, np.sqrt(6)/4, 0.5], 1e-12, "Ex 4.4 n")
wE, vE = np.linalg.eigh(H2.mat)
close(abs(vE[:, 1].conj() @ np.array([2, 1]) / np.sqrt(5)), 1.0, 1e-12, "Ex 4.5 |+>")
print(f"Ex 4.4: n = {np.round(bloch(spin), 5)}; Ex 4.5: |+> = (2,1)/sqrt5, theta = {np.degrees(np.arccos(0.6)):.2f} deg")
close(outer @ phi_a.vec, psi_a.vec, 1e-14, "A1: M|phi> = |psi>")

# B1: common eigenbasis of commuting A = diag(1,1,2), B
Ab = Operator(np.diag([1., 1., 2.])); Bb = Operator([[0, 1, 0], [1, 0, 0], [0, 0, 5]])
close(comm(Ab, Bb).mat, 0, 1e-15, "B1 commute")
common = np.array([[1, 1, 0], [1, -1, 0], [0, 0, np.sqrt(2)]]).T / np.sqrt(2)
for v, (aa, bb) in zip(common.T, [(1, 1), (1, -1), (2, 5)]):
    close(Ab.mat @ v, aa * v, 1e-14, "B1 a"); close(Bb.mat @ v, bb * v, 1e-14, "B1 b")
print("B1: common eigenbasis (1,1,0)/sqrt2 (1,1), (1,-1,0)/sqrt2 (1,-1), e3 (2,5)")

# B5: two-flavor oscillation from the mass basis (3.23) vs Rabi (3.30)
th_m, d_m = 0.6, 0.37                               # mixing angle, (E2 - E1)/2
c_, s_m = np.cos(th_m), np.sin(th_m)
Hfl = Operator(d_m * (-np.cos(2*th_m) * SZ + np.sin(2*th_m) * SX))
mu_ = State([0, 1])
for t in (0.5, 2.0, 7.3):
    direct = abs(c_ * s_m * (np.exp(-1j * d_m * t) - np.exp(1j * d_m * t)))**2
    close(abs(State([1, 0]).inner(Hfl.U(t) @ mu_))**2, direct, 1e-12, "B5 direct")
    close(direct, np.sin(2*th_m)**2 * np.sin(d_m * t)**2, 1e-12, "B5 formula")
print("B5: flavor-basis H gives sin^2(2 theta) sin^2(dm^2 c^4 t / 4 hbar E)")
print(f"MRI: proton Larmor frequency at 1.5 T = {1.5 * PC['proton gyromag. ratio in MHz/T'][0]:.2f} MHz")

# B2: chirped Gaussian a = 1, c = -1/4 saturates Schrodinger-Robertson
ac, cc = 1.0, -0.25
chirp = State(np.exp(-xg**2/(4*ac**2) + 1j*cc*xg**2))
sx_, sp_ = Xg.sigma(chirp), Pspec.sigma(chirp)
cov = 0.5 * (Xg @ Pspec + Pspec @ Xg).expect(chirp).real
close([sx_*sp_, cov], [np.sqrt(2)/2, 2*cc*ac**2], 1e-8, "B2")
close((sx_*sp_)**2, cov**2 + 0.25, 1e-8, "B2 saturation")
print(f"B2: sigma_x sigma_p = {sx_*sp_:.5f}, covariance = {cov:.5f}; SR bound saturated")

# B4: free Gaussian a = 1 at t = 2: sigma_x = sqrt2 via Heisenberg picture x + p t / m
Ufree = Operator(F.conj().T @ np.diag(np.exp(-0.5j * kg**2 * 2.0)) @ F)
gfree = State(np.exp(-xg**2 / 4))
close(Xg.sigma(Ufree @ gfree), np.sqrt(2), 1e-8, "B4")
XH = Xg + 2.0 * Pspec
close(XH.sigma(gfree), np.sqrt(2), 1e-8, "B4 Heisenberg")
print(f"B4: sigma_x(2) = {Xg.sigma(Ufree @ gfree):.5f} (Schrodinger) = {XH.sigma(gfree):.5f} (Heisenberg)")

# B5 and notes 3.12: ammonia, neutrinos, Larmor
nu_NH3 = 23.870e9
A_NH3 = H_SI * nu_NH3 / 2 / EV
print(f"NH3: 2A = {2*A_NH3:.4e} eV, transfer time 1/(2 nu) = {1/(2*nu_NH3)*1e12:.2f} ps")
hbarc_eVm = PC["reduced Planck constant times c in MeV fm"][0] * 1e6 * 1e-15
k_osc = 2.5e-3 / (4 * 1e9 * hbarc_eVm) * 1e3          # rad per km at E = 1 GeV
print(f"B5: phase per km = {k_osc:.5e}, first maximum L = {np.pi/2/k_osc:.1f} km; "
      f"coefficient 1/(4 hbar c) = {1e3/(4*1e9*hbarc_eVm):.4f} per (eV^2 km / GeV)")
close(np.pi/2/k_osc, 495.9, 0.1, "B5 L")
f_e = PC["electron gyromag. ratio in MHz/T"][0]; f_p = PC["proton gyromag. ratio in MHz/T"][0]
print(f"Larmor: electron {f_e/1e3:.3f} GHz/T, proton {f_p:.3f} MHz/T")

# B7: Mandelstam-Tamm along the eps=0.3, Delta=0.4 trajectory, and saturation at eps = 0
SZo = Operator(SZ)
close(H2.sigma(up), 0.4, 1e-12, "sigma_H")
for t in np.linspace(0, 10, 101):
    st = H2.U(t) @ up
    dz = ((1j/HBAR) * comm(H2, SZo).expect(st)).real
    assert H2.sigma(st) * SZo.sigma(st) >= 0.5 * HBAR * abs(dz) - 1e-12
H0 = Operator(Dl * SX)
t_orth = np.pi * HBAR / (2 * H0.sigma(up))
close(abs(up.inner(H0.U(t_orth) @ up)), 0, 1e-12, "orthogonal at pi hbar / 2 sigma_H")
print(f"B7: sigma_H = 0.4, Mandelstam-Tamm holds; eps=0 reaches |down> at t = {t_orth:.4f} = pi/(2 sigma_H)")

# C1: driven two-state, RWA prediction sin^2(Omega t / 2)
w0, Om = 1.0, 0.05
from scipy.integrate import solve_ivp
rhs = lambda t, c: -1j * ((w0/2) * SZ + Om*np.cos(w0*t)*SX) @ c
sol = solve_ivp(rhs, [0, np.pi/Om], np.array([0, 1], complex), rtol=1e-9, atol=1e-11, dense_output=True)
Pfull = abs(sol.y[0, -1])**2
assert Pfull > 0.99
print(f"C1: full driven evolution, P(down->up) at t = pi/Omega = {np.pi/Om:.2f}: {Pfull:.5f} (RWA: 1)")

# C3: mixed state purity
rho = 0.5 * np.outer([1, 0], [1, 0]) + 0.5 * np.outer([1, 1], [1, 1]) / 2
r = np.real([np.trace(rho @ S) for S in (SX, SY, SZ)])
close([np.trace(rho @ rho).real, np.linalg.norm(r)], [0.75, np.sqrt(0.5)], 1e-12, "C3")
lam_rho = np.linalg.eigvalsh(rho)
close(lam_rho, [(1 - np.sqrt(0.5))/2, (1 + np.sqrt(0.5))/2], 1e-12, "C3 eigenvalues")
print(f"C3: Tr rho^2 = 0.75, r = {np.round(r, 3)}, |r| = {np.linalg.norm(r):.5f}, eigenvalues {np.round(lam_rho, 5)}")

# %% [markdown]
# ## Your turn
# 1. Interaction picture: for C1's driven spin, compute |c_I(t)> = e^{i H0 t}|psi(t)> and
#    check it obeys i d|c_I>/dt = V_I(t)|c_I>; then drop the counter-rotating terms (RWA) and
#    compare. A reference check follows.
# 2. Density matrix (QM-7 preview): for rho = sum_k w_k |psi_k><psi_k|, verify <Q> = Tr(rho Q)
#    for random states and weights, and that Tr rho^2 < 1 unless one weight is 1.
# 3. Bloch's theorem: build a ring of M = 8 cells with a narrow Gaussian well in each (a
#    smoothed Dirac comb), check [T_cell, H] = 0, diagonalize, and show each eigenvector
#    satisfies psi_{j+n} = e^{i k a} psi_j with k a = 2 pi m / M. A reference check follows.

# %%
tI = 7.0
cS = sol.sol(tI)
cI = expm(1j * (w0/2) * SZ * tI) @ cS
VI = lambda t: expm(1j*(w0/2)*SZ*t) @ (Om*np.cos(w0*t)*SX) @ expm(-1j*(w0/2)*SZ*t)
dcI = (expm(1j*(w0/2)*SZ*(tI+1e-5)) @ sol.sol(tI+1e-5) - expm(1j*(w0/2)*SZ*(tI-1e-5)) @ sol.sol(tI-1e-5)) / 2e-5
close(1j * dcI, VI(tI) @ cI, 1e-5, "interaction picture")

ws = rng.dirichlet(np.ones(3)); kets = [rand_state(4) for _ in range(3)]
rho4 = sum(w_ * np.outer(k_.vec, k_.bra()) for w_, k_ in zip(ws, kets))
close(np.trace(rho4 @ A.mat), sum(w_ * A.expect(k_) for w_, k_ in zip(ws, kets)), 1e-12, "Tr rho Q")
assert np.trace(rho4 @ rho4).real < 1

M, n_cell = 8, 24
Nb = M * n_cell; db = 1.0 / n_cell
xb = np.arange(Nb) * db
dper = lambda y, c: (y - c + M/2) % M - M/2               # periodic distance on the ring
Vb = -40 * sum(np.exp(-(dper(xb, j + 0.5) / 0.05)**2 / 2) for j in range(M))
Tb = np.roll(np.eye(Nb), 1, axis=0)
Hb = -0.5 * (Tb + Tb.T - 2*np.eye(Nb)) / db**2 + np.diag(Vb)
Tcell = np.linalg.matrix_power(Tb, n_cell)
close(Tcell @ Hb - Hb @ Tcell, 0, 1e-8, "Bloch [T,H]")
lam, W = np.linalg.eigh(Hb + 1e-4 * 1j * (Tcell - Tcell.T))  # Hermitian, commutes with H: splits +-k
for k in range(16):
    v = W[:, k]
    j0 = np.argmax(abs(v)); ratio_ = (Tcell @ v)[j0] / v[j0]
    close(Tcell @ v, ratio_ * v, 1e-8, "Bloch eigvec")
    assert abs(abs(ratio_) - 1) < 1e-8 and abs(ratio_**M - 1) < 1e-8
print("Your turn reference: interaction picture, Tr(rho Q), Bloch's theorem on an 8-cell comb all check")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T0:.1f} s.")
