# %% [markdown]
# # QM-2 Project: A Split-Operator Schrödinger Solver and Everything It Can Check
#
# This script builds a 1D time-dependent Schrödinger solver with the
# split-operator (FFT) method and uses it to check the module's results on a
# grid: norm and energy conservation, the analytic spreading of a free
# Gaussian packet, the fixed momentum density, phase vs group velocity, the
# continuity equation, Ehrenfest's theorem (harmonic and quartic), the
# uncertainty product, partial reflection at a step, and the modified
# Hamilton–Jacobi equation with the quantum potential. SymPy checks the
# worked-example integrals and the exact Gaussian packet. A final part checks
# every number quoted in the exercise solutions.
#
# Simulation units: hbar = m = 1 (and omega = 1 for the oscillator).
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy import constants as C
from scipy.integrate import quad, solve_ivp
from scipy.special import erfc

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.time()

# %% [markdown]
# ## Part 1: Scales (§3.12)
# hbar, the Compton wavelength, thermal de Broglie wavelengths vs spacing.

# %%
hbar, h, me, kB, c = C.hbar, C.h, C.m_e, C.k, C.c
lam_C = h / (me * c)
print(f"hbar = {hbar:.6e} J s; electron Compton wavelength h/(m c) = {lam_C:.4e} m, "
      f"reduced hbar/(m c) = {hbar/(me*c):.4e} m")


def lam_th(m, T):
    """Thermal de Broglie wavelength h / sqrt(2 pi m k T), eq. (3.35)."""
    return h / np.sqrt(2 * np.pi * m * kB * T)


lam_e = lam_th(me, 300.0)
m_N2 = 28.0134e-3 / C.N_A
lam_N2 = lam_th(m_N2, 300.0)
n_air = C.atm / (kB * 300.0)
print(f"lambda_th: electron 300 K = {lam_e:.3e} m; N2 300 K = {lam_N2:.3e} m; "
      f"air spacing n^(-1/3) = {n_air**(-1/3):.3e} m; n lambda^3 (N2) = {n_air*lam_N2**3:.2e}")
n_Cu = 8.49e28
print(f"Cu conduction electrons: spacing {n_Cu**(-1/3):.3e} m vs lambda_th {lam_e:.3e} m")
assert abs(lam_C - 2.426e-12) < 1e-15 and abs(lam_e - 4.30e-9) < 1e-11
assert abs(lam_N2 - 1.91e-11) < 1e-13 and n_air * lam_N2**3 < 1e-6

# %% [markdown]
# ## Part 2: SymPy checks of the analytic results
# Example 4.1 (the parabola), the exact Gaussian packet (3.25) solving the
# free equation and having width a sqrt(1 + beta^2) (3.26), and the Madelung
# split (3.32)-(3.33) of §3.11.

# %%
x, t, p = sp.symbols("x t p", real=True)
a, m, hb, k0, L = sp.symbols("a m hbar k_0 L", positive=True)

# Example 4.1: Psi = A (a^2 - x^2) on [-a, a]
A2 = 1 / sp.integrate((a**2 - x**2)**2, (x, -a, a))
psi = sp.sqrt(A2) * (a**2 - x**2)
ex = sp.simplify(sp.integrate(x * psi**2, (x, -a, a)))
ex2 = sp.simplify(sp.integrate(x**2 * psi**2, (x, -a, a)))
ep2 = sp.simplify(sp.integrate(psi * (-hb**2) * sp.diff(psi, x, 2), (x, -a, a)))
ep2_alt = sp.simplify(hb**2 * sp.integrate(sp.diff(psi, x)**2, (x, -a, a)))
assert sp.simplify(A2 - sp.Rational(15, 16) / a**5) == 0
assert ex == 0 and sp.simplify(ex2 - a**2 / 7) == 0
assert sp.simplify(ep2 - sp.Rational(5, 2) * hb**2 / a**2) == 0 and sp.simplify(ep2 - ep2_alt) == 0
prod41 = float(sp.sqrt(sp.Rational(1, 7) * sp.Rational(5, 2)))
print(f"Example 4.1: A^2 = 15/(16a^5), <x^2> = a^2/7, <p^2> = 5hbar^2/(2a^2), "
      f"sigma_x sigma_p = sqrt(5/14) hbar = {prod41:.5f} hbar")
assert abs(prod41 - 0.59761) < 1e-5

# Exact free Gaussian packet, eq. (3.25)
beta = hb * t / (2 * m * a**2)
v = hb * k0 / m
lnPsi = (-sp.log(1 + sp.I * beta) / 2 - (x - v * t)**2 / (4 * a**2 * (1 + sp.I * beta))
         + sp.I * k0 * x - sp.I * hb * k0**2 * t / (2 * m))      # log Psi up to a constant
resid = (sp.I * hb * sp.diff(lnPsi, t)
         + hb**2 / (2 * m) * (sp.diff(lnPsi, x)**2 + sp.diff(lnPsi, x, 2)))   # (i hbar Psi_t - H Psi)/Psi
assert sp.cancel(sp.together(sp.expand(resid))) == 0
expo = -(x - v * t)**2 / (4 * a**2 * (1 + sp.I * beta))
assert sp.cancel(sp.together(expo + expo.subs(sp.I, -sp.I) + (x - v * t)**2 / (2 * a**2 * (1 + beta**2)))) == 0
s2 = sp.symbols("s2", positive=True)
gauss = sp.exp(-x**2 / (2 * s2)) / sp.sqrt(2 * sp.pi * s2)
assert sp.simplify(sp.integrate(gauss, (x, -sp.oo, sp.oo)) - 1) == 0
assert sp.simplify(sp.integrate(x**2 * gauss, (x, -sp.oo, sp.oo)) - s2) == 0   # so sigma^2 = a^2(1+beta^2)
print("Gaussian packet (3.25) solves the free equation; sigma_x(t)^2 = a^2 (1 + (hbar t/2ma^2)^2)")

# Madelung: real part of the equation for Psi = R e^{iS/hbar} is H-J + Q (3.32)
Rf, Sf, Vf = (sp.Function(n)(x, t) for n in ("R", "S", "V"))
Psi_m = Rf * sp.exp(sp.I * Sf / hb)
eq = (sp.I * hb * sp.diff(Psi_m, t) + hb**2 / (2 * m) * sp.diff(Psi_m, x, 2) - Vf * Psi_m)
eq = sp.expand(sp.simplify(eq * sp.exp(-sp.I * Sf / hb)))
parts = sp.collect(eq, sp.I, evaluate=False)
re_part, im_part = parts[sp.S.One], parts[sp.I]
hj = Rf * (sp.diff(Sf, t) + sp.diff(Sf, x)**2 / (2 * m) + Vf - hb**2 / (2 * m) * sp.diff(Rf, x, 2) / Rf)
assert sp.simplify(re_part + hj) == 0
cont = (sp.diff(Rf**2, t) + sp.diff(Rf**2 * sp.diff(Sf, x) / m, x))
assert sp.simplify(im_part * 2 * Rf / hb - cont) == 0
print("Madelung split checked: real part = H-J + quantum potential, imaginary part = continuity")

# %% [markdown]
# ## Part 3: The split-operator solver (§7 and the derivation below)
#
# Over one step dt, exp(-i H dt/hbar) with H = T + V. T = p^2/2m is diagonal
# in momentum space (F3: the derivative becomes multiplication by i k), V is
# diagonal in position space. Strang splitting,
#   exp(-i H dt/hbar) = exp(-i V dt/2hbar) exp(-i T dt/hbar) exp(-i V dt/2hbar) + O(dt^3),
# because the Baker–Campbell–Hausdorff terms of order dt^2 cancel between the
# symmetric halves. Each factor is a unimodular phase multiplying the
# coefficients in its own basis, and the FFT (unitary up to 1/N) moves between
# bases, so the scheme conserves the norm exactly (to rounding). The time step
# must resolve the fastest phase: (hbar k_max^2/2m) dt and V_max dt/hbar should
# be well under pi, or the splitting error dominates.

# %%
class Grid:
    """Periodic grid x_n = x_min + n dx with FFT wavenumbers k_j (hbar = m = 1)."""

    def __init__(self, x_min, x_max, n):
        self.x = np.linspace(x_min, x_max, n, endpoint=False)
        self.dx = self.x[1] - self.x[0]
        self.k = 2 * np.pi * np.fft.fftfreq(n, self.dx)

    def deriv(self, f, order=1):
        """Spectral derivative d^order f/dx^order (F3 rule (3.8) on the DFT)."""
        return np.fft.ifft((1j * self.k)**order * np.fft.fft(f))

    def norm(self, psi):
        return np.sum(np.abs(psi)**2) * self.dx

    def mean_x(self, psi, f=lambda x: x):
        return np.sum(f(self.x) * np.abs(psi)**2) * self.dx

    def p_stats(self, psi):
        """<p>, <p^2> from the momentum-space density |Phi(p)|^2, p = hbar k, eq. (3.18)."""
        w = np.abs(np.fft.fft(psi))**2
        w /= w.sum()
        return np.sum(self.k * w), np.sum(self.k**2 * w)

    def energy(self, psi, V):
        _, p2 = self.p_stats(psi)
        return 0.5 * p2 + self.mean_x(psi, lambda _: V)

    def current(self, psi):
        """J = (hbar/m) Im(psi* dpsi/dx), eq. (3.4)."""
        return np.imag(np.conj(psi) * self.deriv(psi))


def evolve(g, psi, V, dt, n_steps, every=1, callback=None):
    """Strang split-operator steps, eq. (7.1). Returns final psi; calls callback(n, psi)."""
    eV = np.exp(-0.5j * V * dt)
    eT = np.exp(-0.5j * g.k**2 * dt)
    for n in range(1, n_steps + 1):
        psi = eV * np.fft.ifft(eT * np.fft.fft(eV * psi))
        if callback is not None and n % every == 0:
            callback(n, psi)
    return psi


def gaussian(x, x0, k0, a):
    """Normalized Gaussian with sigma_x = a, mean momentum hbar k0, eq. (3.23)."""
    return (2 * np.pi * a**2)**-0.25 * np.exp(-(x - x0)**2 / (4 * a**2) + 1j * k0 * (x - x0))


def gaussian_exact(x, t, x0, k0, a):
    """Exact free evolution of gaussian(), eq. (3.25) (hbar = m = 1, phase origin at x0)."""
    q = 1 + 1j * t / (2 * a**2)
    return ((2 * np.pi * a**2)**-0.25 / np.sqrt(q)
            * np.exp(-(x - x0 - k0 * t)**2 / (4 * a**2 * q) + 1j * k0 * (x - x0) - 0.5j * k0**2 * t))


# %% [markdown]
# ## Part 4: The free Gaussian packet (§3.8, Example 4.4)
# Width vs (3.26); |Phi(p)|^2 fixed; sigma_x sigma_p growing.

# %%
g = Grid(-100, 100, 2048)
a0, kk0, x00 = 1.0, 2.0, -30.0
psi0 = gaussian(g.x, x00, kk0, a0)
dens_p0 = np.abs(np.fft.fft(psi0))**2
dt, nst = 0.01, 2000
rec = {"t": [0.0], "sx": [a0], "prod": [0.5], "norm": [g.norm(psi0)]}
snap = {}


def free_cb(n, psi):
    tt = n * dt
    mx = g.mean_x(psi)
    sx = np.sqrt(g.mean_x(psi, lambda y: (y - mx)**2))
    mp, mp2 = g.p_stats(psi)
    rec["t"].append(tt); rec["sx"].append(sx); rec["norm"].append(g.norm(psi))
    rec["prod"].append(sx * np.sqrt(mp2 - mp**2))
    if n in (500, 1000, 2000):
        snap[n] = psi.copy()


psiT = evolve(g, psi0, np.zeros_like(g.x), dt, nst, callback=free_cb)
tt = np.array(rec["t"]); sx = np.array(rec["sx"])
sx_th = a0 * np.sqrt(1 + (tt / (2 * a0**2))**2)
assert np.max(np.abs(sx - sx_th)) < 1e-8, np.max(np.abs(sx - sx_th))
assert np.max(np.abs(np.array(rec["norm"]) - 1)) < 1e-12
assert np.max(np.abs(np.abs(np.fft.fft(psiT))**2 - dens_p0)) < 1e-9 * dens_p0.max()
assert np.max(np.abs(psiT - gaussian_exact(g.x, nst * dt, x00, kk0, a0))) < 1e-9
prod = np.array(rec["prod"])
assert np.allclose(prod, 0.5 * np.sqrt(1 + (tt / 2)**2), atol=1e-8)
print(f"Free packet: sigma_x(20) = {sx[-1]:.5f} (theory {sx_th[-1]:.5f}); norm drift "
      f"{np.max(np.abs(np.array(rec['norm']) - 1)):.1e}; sigma_x sigma_p(20) = {prod[-1]:.5f}")

# Example 4.4: Phi(p) of the Gaussian, continuous transform from the FFT
pgrid = np.fft.fftshift(g.k)
Phi_num = np.fft.fftshift(g.dx / np.sqrt(2 * np.pi) * np.exp(-1j * g.k * g.x[0]) * np.fft.fft(gaussian(g.x, 0, 0, a0)))
Phi_th = (2 * a0**2 / np.pi)**0.25 * np.exp(-a0**2 * pgrid**2)
assert np.max(np.abs(Phi_num - Phi_th)) < 1e-10
dp = pgrid[1] - pgrid[0]
assert abs(np.sum(np.abs(Phi_num)**2) * dp - 1) < 1e-10        # Plancherel (3.17)
sp_num = np.sqrt(np.sum(pgrid**2 * np.abs(Phi_num)**2) * dp)
print(f"Example 4.4: sigma_p = {sp_num:.6f} = hbar/(2a); product = {a0*sp_num:.6f} hbar")
assert abs(sp_num - 0.5) < 1e-10

fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
for n, col in ((0, "k"), (500, "C0"), (1000, "C1"), (2000, "C2")):
    ps = psi0 if n == 0 else snap[n]
    ax[0].plot(g.x, np.abs(ps)**2, col, label=f"t = {n*dt:.0f}")
ax[0].set(xlim=(-45, 45), xlabel="x", ylabel=r"$|\Psi|^2$", title="Free Gaussian packet")
ax[0].legend()
ax[1].plot(tt, sx, label="numerical")
ax[1].plot(tt, sx_th, "--", label=r"$a\sqrt{1+(\hbar t/2ma^2)^2}$")
ax[1].plot(tt, 2 * prod, ":", label=r"$2\sigma_x\sigma_p/\hbar$")
ax[1].set(xlabel="t", title="Spreading"); ax[1].legend()
fig.tight_layout(); fig.savefig(FIG / "free_packet.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 5: Phase vs group velocity (§3.7, eq. (3.22))
# A wide packet (a = 5): crests move at hbar k0/2m, the envelope at hbar k0/m.

# %%
aw = 5.0
psiw = gaussian(g.x, x00, kk0, aw)
times = [0.0, 4.0, 8.0, 12.0, 16.0, 20.0]
fig, axs = plt.subplots(2, 3, figsize=(12, 5), sharey=True)
vp_meas = []
cur = psiw.copy()
for i, tk in enumerate(times):
    if i > 0:
        cur = evolve(g, cur, np.zeros_like(g.x), 0.01, 400)
    axx = axs.flat[i]
    axx.plot(g.x, cur.real, lw=0.8); axx.plot(g.x, np.abs(cur), "k")
    axx.axvline(x00 + kk0 * tk, color="C3", ls="--", lw=1)          # group
    axx.axvline(x00 + 0.5 * kk0 * tk, color="C2", ls=":", lw=1.5)   # a crest at phase speed
    axx.set(xlim=(x00 - 10 + 0.5 * kk0 * tk, x00 + 10 + kk0 * tk), title=f"t = {tk:.0f}")
axs[0, 0].set_ylabel(r"Re $\Psi$, $|\Psi|$")
fig.suptitle("red dashed: envelope at $v_g=\\hbar k_0/m$; green dotted: crest at $v_p=\\hbar k_0/2m$")
fig.tight_layout(); fig.savefig(FIG / "phase_group.png", dpi=120); plt.close(fig)

# measure v_p = -S_t/S_x at the packet center at t = 10
p_a = evolve(g, psiw, np.zeros_like(g.x), 0.01, 999)
p_b = evolve(g, p_a, np.zeros_like(g.x), 0.01, 1)
p_c = evolve(g, p_b, np.zeros_like(g.x), 0.01, 1)
j = np.argmin(np.abs(g.x - (x00 + kk0 * 10.0)))
S_t = np.angle(p_c[j] * np.conj(p_a[j])) / 0.02
S_x = np.imag(np.conj(p_b) * g.deriv(p_b))[j] / np.abs(p_b[j])**2
v_phase = -S_t / S_x
v_group = (g.mean_x(p_c) - g.mean_x(p_a)) / 0.02
print(f"Phase velocity at packet center = {v_phase:.4f} (k0/2 = 1); group velocity = {v_group:.6f} (k0 = 2)")
assert abs(v_phase - 1.0) < 0.01 and abs(v_group - 2.0) < 1e-6

# %% [markdown]
# ## Part 6: Probability current (3.4) and continuity (3.5) (§3.2, Example 4.2)

# %%
# plane wave and standing wave on the grid (k = 2*pi*5/L exactly periodic)
kp = g.k[5]
assert np.allclose(g.current(np.exp(1j * kp * g.x)), kp)          # J = |A|^2 hbar k / m
assert np.allclose(g.current(np.cos(kp * g.x) + 0j), 0, atol=1e-10)  # standing wave
v_e = hbar * 1e10 / me
print(f"Example 4.2: electron with k = 1e10 /m has J/rho = hbar k/m = {v_e:.4e} m/s; standing wave J = 0")
assert abs(v_e - 1.158e6) < 1e3

# continuity on the grid for a packet in a harmonic well
Vh = 0.5 * g.x**2 * 0.01                     # weak well, omega = 0.1
dtc = 0.005
pm = evolve(g, gaussian(g.x, -20, 1.0, 2.0), Vh, dtc, 400)
pc = evolve(g, pm, Vh, dtc, 1)
pp = evolve(g, pc, Vh, dtc, 1)
drho = (np.abs(pp)**2 - np.abs(pm)**2) / (2 * dtc)
divJ = np.real(g.deriv(g.current(pc)))
res_c = np.max(np.abs(drho + divJ)) / np.max(np.abs(drho))
print(f"Continuity residual max|d rho/dt + dJ/dx| / max|d rho/dt| = {res_c:.2e}")
assert res_c < 1e-4

# %% [markdown]
# ## Part 7: Ehrenfest, energy conservation, minimum uncertainty (§3.4, Example 4.6)

# %%
gh = Grid(-12, 12, 512)


def run_ehrenfest(V, dV, psi, dt, n):
    out = {"t": [0.0], "x": [gh.mean_x(psi)], "p": [gh.p_stats(psi)[0]],
           "F": [-gh.mean_x(psi, dV)], "E": [gh.energy(psi, V(gh.x))], "prod": []}

    def cb(k, ps):
        out["t"].append(k * dt); out["x"].append(gh.mean_x(ps)); out["p"].append(gh.p_stats(ps)[0])
        out["F"].append(-gh.mean_x(ps, dV)); out["E"].append(gh.energy(ps, V(gh.x)))
        mx = out["x"][-1]
        mp, mp2 = gh.p_stats(ps)
        out["prod"].append(np.sqrt(gh.mean_x(ps, lambda y: (y - mx)**2) * (mp2 - mp**2)))

    evolve(gh, psi, V(gh.x), dt, n, callback=cb)
    return {k_: np.array(v_) for k_, v_ in out.items()}


Vho, dVho = (lambda y: 0.5 * y**2), (lambda y: y)
Vq, dVq = (lambda y: 0.25 * y**4), (lambda y: y**3)
coh = np.pi**-0.25 * np.exp(-(gh.x - 2.0)**2 / 2) + 0j       # displaced ground state, x0 = 2
dte = 0.005
ho = run_ehrenfest(Vho, dVho, coh, dte, 2000)
qu = run_ehrenfest(Vq, dVq, coh, dte, 2000)
for name, r in (("harmonic", ho), ("quartic", qu)):
    e1 = np.max(np.abs(np.gradient(r["x"], dte) - r["p"])[1:-1]) / np.max(np.abs(r["p"]))
    e2 = np.max(np.abs(np.gradient(r["p"], dte) - r["F"])[1:-1]) / np.max(np.abs(r["F"]))
    dE = np.max(np.abs(r["E"] - r["E"][0])) / r["E"][0]
    print(f"Ehrenfest ({name}): m d<x>/dt vs <p> rel err {e1:.1e}; d<p>/dt vs <-V'> {e2:.1e}; energy drift {dE:.1e}")
    assert e1 < 1e-4 and e2 < 5e-4 and dE < 1e-4
i90 = np.argmin(np.abs(ho["t"] - np.pi / 2))
assert np.allclose(ho["x"], 2 * np.cos(ho["t"]), atol=1e-4)
assert np.allclose(ho["p"], -2 * np.sin(ho["t"]), atol=1e-4)
assert abs(ho["E"][0] - 2.5) < 1e-10
assert np.allclose(ho["prod"], 0.5, atol=1e-5)      # coherent state stays minimum-uncertainty
print(f"Example 4.6: max|<x> - 2cos t| = {np.max(np.abs(ho['x'] - 2*np.cos(ho['t']))):.1e}, "
      f"<p>(t={ho['t'][i90]:.3f}) = {ho['p'][i90]:.4f}, <H> = {ho['E'][0]:.4f}, "
      f"sigma_x sigma_p in [{ho['prod'].min():.6f}, {ho['prod'].max():.6f}]")

cl = solve_ivp(lambda _, y: [y[1], -y[0]**3], (0, 10), [2.0, 0.0], t_eval=qu["t"], rtol=1e-10, atol=1e-12)
dev = np.max(np.abs(cl.y[0] - qu["x"]))
VpAtMean = -qu["x"]**3
print(f"Quartic: max |<x>_quantum - x_classical| over t<10 = {dev:.3f}; "
      f"max |<-V'> - (-V'(<x>))| = {np.max(np.abs(qu['F'] - VpAtMean)):.3f}")
assert dev > 0.3
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
ax[0].plot(ho["t"], ho["x"], label=r"$\langle x\rangle$"); ax[0].plot(ho["t"], 2 * np.cos(ho["t"]), "--", label="classical")
ax[0].set(title=r"$V=x^2/2$", xlabel="t"); ax[0].legend()
ax[1].plot(qu["t"], qu["x"], label=r"$\langle x\rangle$"); ax[1].plot(qu["t"], cl.y[0], "--", label="classical")
ax[1].set(title=r"$V=x^4/4$", xlabel="t"); ax[1].legend()
fig.tight_layout(); fig.savefig(FIG / "ehrenfest.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 8: A packet hits a step (QM-3 preview)
# E > V0 on average, yet a fraction reflects. Compare with the plane-wave
# R(k) = ((k - q)/(k + q))^2, q = sqrt(k^2 - 2 V0), averaged over |phi(k)|^2.

# %%
gs = Grid(-250, 250, 4096)
V0, ks, aS = 1.5, 2.0, 5.0
Vstep = np.where(gs.x > 0, V0, 0.0)
psiS = evolve(gs, gaussian(gs.x, -100, ks, aS), Vstep, 0.02, 5000)
R_sim = np.sum(np.abs(psiS[gs.x < 0])**2) * gs.dx


def R_plane(k):
    if k**2 <= 2 * V0:
        return 1.0
    q = np.sqrt(k**2 - 2 * V0)
    return ((k - q) / (k + q))**2


wk = lambda k: np.sqrt(2 * aS**2 / np.pi) * np.exp(-2 * aS**2 * (k - ks)**2)   # |phi(k)|^2
R_pred = quad(lambda k: wk(k) * R_plane(k), ks - 1, ks + 1, points=[np.sqrt(2 * V0)])[0]
E_mean = 0.5 * (ks**2 + 1 / (4 * aS**2))
assert abs(E_mean - 2.005) < 1e-12
print(f"Step: <E> = {E_mean:.3f}; reflected probability {R_sim:.5f}; plane-wave average {R_pred:.5f}; R(k0) = {R_plane(ks):.5f}")
assert abs(R_sim - R_pred) < 2e-3 and abs(R_plane(ks) - 1 / 9) < 1e-12
fig, ax = plt.subplots(figsize=(8, 3))
ax.plot(gs.x, np.abs(psiS)**2); ax.axvspan(0, 250, color="0.9"); ax.set(xlabel="x", ylabel=r"$|\Psi|^2$",
                                                                        title=f"after the step: R = {R_sim:.3f}")
fig.tight_layout(); fig.savefig(FIG / "step.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 9: The Hamilton–Jacobi check (§3.11)
# S = hbar arg Psi, rho = |Psi|^2. Check S_t + S_x^2/2m + V + Q = 0 with
# Q = -(hbar^2/2m) (sqrt rho)''/sqrt rho, eq. (3.32), and that dropping Q fails.

# %%
psH = evolve(g, gaussian(g.x, -10, 1.0, 1.0), np.zeros_like(g.x), 0.001, 2999)
psH2 = evolve(g, psH, np.zeros_like(g.x), 0.001, 1)
psH3 = evolve(g, psH2, np.zeros_like(g.x), 0.001, 1)
St = np.angle(psH3 * np.conj(psH)) / 0.002
Sx = np.imag(np.conj(psH2) * g.deriv(psH2)) / np.abs(psH2)**2
sq = np.abs(psH2)
Q = -0.5 * np.real(g.deriv(sq + 0j, 2)) / sq
mx = g.mean_x(psH2); sxx = np.sqrt(g.mean_x(psH2, lambda y: (y - mx)**2))
win = np.abs(g.x - mx) < 2 * sxx
r_q = np.max(np.abs(St + 0.5 * Sx**2 + Q)[win]) / np.max(np.abs(St[win]))
r_c = np.max(np.abs(St + 0.5 * Sx**2)[win]) / np.max(np.abs(St[win]))
print(f"H-J at t = 3: relative residual with Q {r_q:.1e}; without Q {r_c:.3f}")
assert r_q < 1e-4 and r_c > 0.02

# %% [markdown]
# ## Part 10: Worked-example numbers in SI (Examples 4.3, 4.5)

# %%
tau_e = 2 * me * (1e-10)**2 / hbar
tau_g = 2 * 1e-9 * (1e-6)**2 / hbar
tau_g1A = 2 * 1e-9 * (1e-10)**2 / hbar
yr = 365.25 * 86400
print(f"Example 4.3: electron a = 1 A: tau = {tau_e:.4e} s; sigma after 1 s = {hbar*1/(2*me*1e-10)/1e3:.1f} km")
print(f"             1 ug grain, a = 1 um: tau = {tau_g:.4e} s = {tau_g/yr:.3e} yr; a = 1 A: {tau_g1A:.3e} s = {tau_g1A/86400:.2f} d")
assert abs(tau_e - 1.7276e-16) < 1e-19 and abs(tau_g / yr - 6.01e5) < 1e3

L_box = 1e-9
E1 = np.pi**2 * hbar**2 / (2 * me * L_box**2)
eV = C.e
f21 = 3 * E1 / h
print(f"Example 4.5: E1 = {E1/eV:.4f} eV, E2 - E1 = {3*E1/eV:.4f} eV, f = {f21:.4e} Hz, "
      f"lambda = {c/f21*1e6:.3f} um, period = {1/f21:.4e} s, <H> = {2.5*E1/eV:.4f} eV")
assert abs(E1 / eV - 0.3760) < 5e-4 and abs(f21 - 2.727e14) < 1e11


def box_rho(xx, tt_, Lb=1.0):
    """|Psi|^2 for (psi1 + psi2)/sqrt2 in a box, hbar = m = 1, eq. (4.5)."""
    E = lambda n: n**2 * np.pi**2 / (2 * Lb**2)
    s = lambda n: np.sqrt(2 / Lb) * np.sin(n * np.pi * xx / Lb)
    ps = (s(1) * np.exp(-1j * E(1) * tt_) + s(2) * np.exp(-1j * E(2) * tt_)) / np.sqrt(2)
    return np.abs(ps)**2


xb = np.linspace(0, 1, 4001)
Tslosh = 2 * np.pi / (3 * np.pi**2 / 2)
xm = lambda tt_: np.trapezoid(xb * box_rho(xb, tt_), xb)
assert abs(xm(0) - (0.5 - 16 / (9 * np.pi**2))) < 1e-7
assert abs(xm(Tslosh / 2) - (0.5 + 16 / (9 * np.pi**2))) < 1e-7 and abs(xm(Tslosh) - xm(0)) < 1e-9
print(f"Box slosh (hbar=m=L=1): period 4/(3 pi) = {Tslosh:.5f}; <x> amplitude 16/(9 pi^2) = {16/(9*np.pi**2):.5f}")

# %% [markdown]
# ## Part 11: Numbers in the exercise solutions

# %%
lam = 1.0
print(f"A1: P(|x| > sigma) = e^(-sqrt2) = {np.exp(-np.sqrt(2)):.4f}")
assert abs(np.exp(-np.sqrt(2)) - 0.2431) < 1e-4
print(f"A1: Gaussian comparison P(|x - mu| > sigma) = {erfc(1/np.sqrt(2)):.4f}")
k_e = 2 * np.pi / 1e-10
print(f"A2: lambda = 1 A: v_g = {hbar*k_e/me:.4e} m/s, v_p = {hbar*k_e/(2*me):.4e} m/s, E = {(hbar*k_e)**2/(2*me):.4e} J = {(hbar*k_e)**2/(2*me)/eV:.2f} eV")
assert abs(hbar * k_e / me - 7.274e6) < 1e3
m_He = 4.002602e-3 / C.N_A
lHe = lam_th(m_He, 4.2)
nHe = 125.0 / m_He
print(f"A4: He-4 at 4.2 K: lambda_th = {lHe:.3e} m, spacing = {nHe**(-1/3):.3e} m, n lambda^3 = {nHe*lHe**3:.2f}; m = {m_He:.4e} kg, n = {nHe:.3e} /m^3")
assert abs(lHe - 4.26e-10) < 2e-12 and abs(nHe**(-1/3) - 3.76e-10) < 2e-12
# B2: two-sided exponential in momentum space
Phi_exp = lambda pp_: np.sqrt(lam / (2 * np.pi)) * 2 * lam / (lam**2 + pp_**2)
nrm = quad(lambda pp_: Phi_exp(pp_)**2, -np.inf, np.inf)[0]
pp2 = quad(lambda pp_: pp_**2 * Phi_exp(pp_)**2, -np.inf, np.inf)[0]
print(f"B2: int |Phi|^2 = {nrm:.8f}, <p^2> = {pp2:.8f} (hbar^2 lambda^2); product = 1/sqrt2 = {1/np.sqrt(2):.4f}")
assert abs(nrm - 1) < 1e-8 and abs(pp2 - 1) < 1e-7
# B4: box <x>(t), <p>(t)
xs = 16 / (9 * np.pi**2)
pm_amp = 8 / 3          # 2 |x12| * (E2 - E1) = (16/9pi^2) * (3 pi^2/2) in hbar=m=L=1
assert abs(xs * 3 * np.pi**2 / 2 - pm_amp) < 1e-12
print(f"B4: <x> = 1/2 - {xs:.5f} cos(w t), <p> = {pm_amp:.4f} sin(w t) (hbar/L), w = 3pi^2/2 = {3*np.pi**2/2:.4f}")
# B8: Strang splitting is second order: energy error ratio for dt and dt/2
psiq = np.pi**-0.25 * np.exp(-(gh.x - 2.0)**2 / 2) + 0j
ref = evolve(gh, psiq, Vq(gh.x), 0.000625, 3200)
errs = []
for dtt in (0.01, 0.005, 0.0025):
    errs.append(np.sqrt(gh.norm(evolve(gh, psiq, Vq(gh.x), dtt, int(round(2.0 / dtt))) - ref)))
ratios = [errs[0] / errs[1], errs[1] / errs[2]]
print(f"B8: state error at t = 2 for dt = 0.01, 0.005, 0.0025: {[f'{e:.2e}' for e in errs]}; ratios {ratios[0]:.2f}, {ratios[1]:.2f}")
assert all(3.5 < r < 4.5 for r in ratios)
# B6: quartic force deficit for a Gaussian centred at x0 width s: <x^3> = x0^3 + 3 x0 s^2
x0q, sq_ = 2.0, 1 / np.sqrt(2)
assert abs(-qu["F"][0] - (x0q**3 + 3 * x0q * sq_**2)) < 1e-8
print(f"B6: at t = 0, <-V'> = {qu['F'][0]:.4f} vs -V'(<x>) = {-x0q**3:.4f}")
# B5: chirped Gaussian (c = -0.25, sigma_x(0) = 1) focuses to a_min = 1/sqrt2 at t* = 1, product hbar/2 there
chirp = (2 * np.pi)**-0.25 * np.exp(-g.x**2 / 4 - 0.25j * g.x**2)
mp_c, mp2_c = g.p_stats(chirp)
assert abs(np.sqrt(mp2_c - mp_c**2) - np.sqrt(0.25 + 4 * 0.0625)) < 1e-10      # sigma_p^2 = 1/4a^2 + 4c^2a^2
wid = []
evolve(g, chirp, np.zeros_like(g.x), 0.001, 2000, callback=lambda n, ps: wid.append((n * 0.001, np.sqrt(g.mean_x(ps, lambda y: y**2)))))
wid = np.array(wid)
i_min = np.argmin(wid[:, 1])
print(f"B5: chirped packet: sigma_p(0) = {np.sqrt(mp2_c):.5f}, min width {wid[i_min,1]:.5f} at t = {wid[i_min,0]:.3f}")
assert abs(wid[i_min, 0] - 1.0) < 2e-3 and abs(wid[i_min, 1] - 1 / np.sqrt(2)) < 1e-6
# B6: quartic deviation vs packet width: time for |<x> - x_cl| to first exceed 0.1
t_dev = {}
for aw_ in (0.3, 1 / np.sqrt(2), 1.0):
    r_ = run_ehrenfest(Vq, dVq, gaussian(gh.x, 2.0, 0.0, aw_), dte, 2000)
    t_dev[aw_] = r_["t"][np.argmax(np.abs(r_["x"] - cl.y[0]) > 0.1)]
    print(f"B6: a = {aw_:.3f}: |<x> - x_cl| first exceeds 0.1 at t = {t_dev[aw_]:.2f}; "
          f"<-V'>(0) = {r_['F'][0]:.3f} = -(8 + 6a^2)")
    assert abs(r_["F"][0] + 8 + 6 * aw_**2) < 1e-8
# B7: virial for Example 4.6's state at t = 0: d<xp>/dt = 2<T> - <x V'> = 2(1/4) - 9/2 = -4
print(f"B7: d<xp>/dt at t = 0 for Example 4.6 = {2*0.25 - (4 + 0.5):.1f} (classical -2 sin 2t)' = -4)")
# C1: uniform force in momentum space, Phi(k,t) = phi0(k - F t) exp(-i [k^3 - (k - F t)^3] / (6 m hbar F))
F, tC = 0.1, 6.0
phiC = (2 / np.pi)**0.25 * np.exp(-(g.k - F * tC)**2) * np.exp(-1j * (g.k**3 - (g.k - F * tC)**3) / (6 * F))
psiC = np.fft.ifft(phiC * np.exp(-1j * g.k * g.x[0])) * len(g.x) * (g.k[1] - g.k[0]) / np.sqrt(2 * np.pi)
psiN = evolve(g, gaussian(g.x, 0, 0, 1.0), -F * g.x, 0.01, 600)
fid = abs(np.sum(np.conj(psiC) * psiN) * g.dx)
print(f"C1: |<Psi_analytic, Psi_numeric>| at t = 6: {fid:.12f}")
assert abs(fid - 1) < 1e-9
# C2: Bohm trajectories x = v t + y0 sigma(t)/a satisfy m x'' = -dQ/dx (symbolic)
y0 = sp.symbols("y0", real=True)
sig = a * sp.sqrt(1 + beta**2)
traj = v * t + y0 * sig / a
yy = traj - v * t
Qf = hb**2 / (2 * m) * (1 / (2 * sig**2) - x**2 / (4 * sig**4))
force = -sp.diff(Qf, x).subs(x, yy)
assert sp.simplify(m * sp.diff(traj, t, 2) - force) == 0
Sx_field = hb * k0 + hb * (x - v * t) * beta / (2 * a**2 * (1 + beta**2))   # S_x from Im of exponent of (3.25)
assert sp.simplify(sp.diff(traj, t) - Sx_field.subs(x, traj) / m) == 0
print("C2: Bohm trajectories x = v t + y0 sigma(t)/a follow S_x/m and obey m x'' = -dQ/dx")
# C3: box ground state uncertainty product
c3 = np.sqrt(np.pi**2 / 12 - 0.5)
print(f"C3: box ground state sigma_x sigma_p = sqrt(pi^2/12 - 1/2) hbar = {c3:.4f} hbar")
assert abs(c3 - 0.5679) < 1e-4
print(f"B8(c): p_max = pi hbar/dx = {np.pi/gh.dx:.1f} for dx = {gh.dx:.4f}")

# %% [markdown]
# ## Your turn
# 1. **Uniform force.** Set V = -F x (F = 0.1) for the free packet of Part 4
#    (k0 = 0) and check Ehrenfest: <x> = x0 + F t^2/2m, <p> = F t, with
#    sigma_x(t) exactly the free one. (Reference below.)
# 2. **Box slosh with the solver.** Put the n = 1, 2 superposition in a box
#    made of a high wall (V = 1e4 outside [0, 1]) and measure the slosh period;
#    compare with 4/(3 pi) from Part 10. Expect a small shift from wall penetration.
# 3. **Momentum-space Schrödinger equation.** For V = 0 solve
#    i hbar dPhi/dt = (p^2/2m) Phi directly (one line), transform back, and
#    compare with gaussian_exact. (Reference below.)
# 4. Exercise B6: rerun Part 7's quartic case with a = 0.3 and 1.5 and watch how
#    the quantum-classical deviation depends on the packet width.

# %%
F = 0.1
lin = {"t": [], "x": [], "p": [], "sx": []}


def lin_cb(n, ps):
    mx_ = g.mean_x(ps)
    lin["t"].append(n * 0.01); lin["x"].append(mx_); lin["p"].append(g.p_stats(ps)[0])
    lin["sx"].append(np.sqrt(g.mean_x(ps, lambda y: (y - mx_)**2)))


evolve(g, gaussian(g.x, 0, 0, 1.0), -F * g.x, 0.01, 600, every=10, callback=lin_cb)
lin = {k_: np.array(v_) for k_, v_ in lin.items()}
assert np.allclose(lin["x"], 0.5 * F * lin["t"]**2, atol=1e-8)
assert np.allclose(lin["p"], F * lin["t"], atol=1e-8)
assert np.allclose(lin["sx"], np.sqrt(1 + (lin["t"] / 2)**2), atol=1e-8)
print(f"Your turn 1: uniform force, <x>(6) = {lin['x'][-1]:.5f} (F t^2/2 = {0.5*F*36:.5f})")

Phi0 = np.fft.fft(gaussian(g.x, x00, kk0, a0))
psi_from_p = np.fft.ifft(Phi0 * np.exp(-0.5j * g.k**2 * 7.0))
assert np.max(np.abs(psi_from_p - gaussian_exact(g.x, 7.0, x00, kk0, a0))) < 1e-9
print("Your turn 3: momentum-space solution Phi(p,0) e^{-i p^2 t/2m hbar} matches (3.25)")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T0:.1f} s.")
