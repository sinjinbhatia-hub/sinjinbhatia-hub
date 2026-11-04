# %% [markdown]
# # QM-1 Project: The Experiments, a Stern-Gerlach Simulator, and the Postulates in C^2
#
# Parts: (1) blackbody curves (Planck, Rayleigh-Jeans, Wien) and Example 4.1; the photoelectric
# law fitted to synthetic stopping-potential data (Example 4.2); (2) Bohr levels and the spectral
# series on a wavelength axis (Example 4.3), the classical-atom collapse time (EM-8 Example 4.4),
# Compton (Example 4.4), de Broglie and Davisson-Germer (Example 4.5), the Stern-Gerlach deflection
# estimate (3.6), Bohr-Sommerfeld for the oscillator; (3) a Stern-Gerlach simulator: Born
# probabilities and collapse along any axis n, 10^4 atoms through sequential analyzers
# (Example 4.6); (4) the Bloch sphere; (5) the double slit with and without which-path tags;
# (6) commutators, the spin uncertainty product, unitary evolution; (7) Your turn with reference
# solutions; (8) every number in solutions.md.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from pathlib import Path
from scipy.linalg import expm
from scipy.optimize import brentq, minimize_scalar
from scipy.constants import h as H, hbar as HBAR, c as C, e as QE, m_e as ME, k as KB
from scipy.constants import epsilon_0 as EPS0, m_p as MP, m_n as MN, u as AMU
from scipy.constants import physical_constants as PC

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.time()
EV = QE
HC_EVNM = H * C / QE * 1e9          # hc in eV nm, 1239.84
MU_B = PC["Bohr magneton"][0]
M_D = PC["deuteron mass"][0]


def close(a, b, rtol, label=""):
    """assert a ~ b to relative tolerance rtol; return the relative error"""
    err = np.max(np.abs(np.asarray(a) - np.asarray(b))) / np.max(np.abs(b))
    assert err < rtol, f"{label}: {a} vs {b}, rel err {err:.2e}"
    return err

# %% [markdown]
# ## Part 1: Blackbody and photoelectric effect (notes 1.1-1.2, Examples 4.1-4.2)
# Energy density per unit wavelength:
# Planck (3.2) u = 8 pi h c / lam^5 / (e^x - 1), Rayleigh-Jeans (3.1) u = 8 pi k T / lam^4, Wien u = 8 pi h c/lam^5 e^-x,
# with x = h c / (lam k T). Their ratio RJ/Planck = (e^x - 1)/x.

# %%
def u_planck(lam, T):
    """Planck spectral energy density u_lambda (J/m^4), notes (3.2)"""
    x = H * C / (lam * KB * T)
    return 8 * np.pi * H * C / lam**5 / np.expm1(x)


def u_rj(lam, T):
    """Rayleigh-Jeans u_lambda = 8 pi k T / lambda^4, notes (3.1)"""
    return 8 * np.pi * KB * T / lam**4


def u_wien(lam, T):
    """Wien's 1896 law u_lambda = 8 pi h c / lambda^5 exp(-x)"""
    return 8 * np.pi * H * C / lam**5 * np.exp(-H * C / (lam * KB * T))


lam0, T5 = 500e-9, 5000.0
x0 = H * C / (lam0 * KB * T5)
up, ur = u_planck(lam0, T5), u_rj(lam0, T5)
print(f"Ex 4.1: 8 pi h c/lam^5 = {8*np.pi*H*C/lam0**5:.4e}, e^x - 1 = {np.expm1(x0):.2f}, "
      f"x(50 um) = {H*C/(50e-6*KB*T5):.4f}")
print(f"Ex 4.1: x = {x0:.4f}, Planck u = {up:.4e} J/m^4, RJ u = {ur:.4e} J/m^4, RJ/Planck = {ur/up:.2f}")
close(ur / up, np.expm1(x0) / x0, 1e-12, "ratio formula")
assert 54.0 < ur / up < 55.5
# long-wavelength agreement: at 50 um the ratio is within 3%
assert abs(u_rj(50e-6, T5) / u_planck(50e-6, T5) - 1) < 0.03
lam_peak = minimize_scalar(lambda L: -u_planck(L * 1e-9, T5), bounds=(300, 900), method="bounded",
                           options={"xatol": 1e-6}).x * 1e-9
xw = brentq(lambda x: x - 5 * (1 - np.exp(-x)), 1, 10)
b_wien = H * C / (KB * xw)
print(f"Wien peak at 5000 K: {b_wien/T5*1e9:.1f} nm (x = {xw:.6f}, b = {b_wien:.4e} m K)")
close(lam_peak, b_wien / T5, 1e-4, "numerical peak")
# energy in Planck spectrum is finite: integral = (8 pi^5 k^4 / 15 h^3 c^3) T^4
lg = np.linspace(20e-9, 200e-6, 400001)
Uint = np.trapezoid(u_planck(lg, T5), lg) if hasattr(np, "trapezoid") else np.trapz(u_planck(lg, T5), lg)
close(Uint, 8 * np.pi**5 * KB**4 / (15 * H**3 * C**3) * T5**4, 2e-3, "Stefan-Boltzmann")

fig, ax = plt.subplots(figsize=(7, 4.2))
L = np.linspace(100e-9, 3000e-9, 600)
ax.plot(L * 1e9, u_planck(L, T5), "k", label="Planck (1900)")
ax.plot(L * 1e9, u_rj(L, T5), "C3--", label="Rayleigh-Jeans")
ax.plot(L * 1e9, u_wien(L, T5), "C0:", label="Wien (1896)")
ax.set_ylim(0, 1.6 * u_planck(b_wien / T5, T5))
ax.axvline(500, color="0.6", lw=0.8)
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel(r"$u_\lambda$ (J/m$^4$)")
ax.set_title("Blackbody at 5000 K: the ultraviolet catastrophe"); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "blackbody.png", dpi=120); plt.close(fig)

# Photoelectric: sodium, W = 2.3 eV, 400 nm light (Example 4.2)
W_NA = 2.3
E400 = HC_EVNM / 400
V400 = E400 - W_NA
lam_thr = HC_EVNM / W_NA
print(f"Ex 4.2: E(400 nm) = {E400:.4f} eV, V_stop = {V400:.4f} V, threshold {lam_thr:.1f} nm, "
      f"f0 = {W_NA*EV/H:.4e} Hz, v_max = {np.sqrt(2*V400*EV/ME):.4e} m/s")
assert abs(V400 - 0.80) < 0.005
# synthetic Millikan-style data: stopping potential vs frequency with 0.02 V noise
lams = np.array([250, 280, 313, 334, 365, 405, 436, 480, 510]) * 1e-9
f = C / lams
V = (H * f / EV - W_NA) + rng.normal(0, 0.02, f.size)
slope, icpt = np.polyfit(f, V, 1)
h_fit = slope * EV
print(f"fitted h = {h_fit:.4e} J s (true {H:.4e}), fitted W = {-icpt:.3f} eV")
close(h_fit, H, 0.02, "fitted h")
assert abs(-icpt - W_NA) < 0.08
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(f / 1e14, V, "ko", label="synthetic data")
ff = np.linspace(W_NA * EV / H, f.max(), 50)
ax.plot(ff / 1e14, slope * ff + icpt, "C0", label=f"fit: h = {h_fit:.3e} J s")
ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlabel("frequency (10$^{14}$ Hz)"); ax.set_ylabel("stopping potential (V)")
ax.set_title("Photoelectric effect, sodium"); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "photoelectric.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 2: Bohr, Compton, de Broglie, Stern-Gerlach numbers (notes 3.2-3.6, Examples 4.3-4.5)
# Bohr: r_n = n^2 a0, a0 = 4 pi eps0 hbar^2/(m e^2) (3.4); E_n = -E_R/n^2, E_R = m e^4/(2 (4 pi eps0)^2 hbar^2) (3.5).
# Rydberg R = E_R/(h c). Compton (3.7): dlam = (h/mc)(1 - cos theta). de Broglie: lam = h/p.

# %%
K = 1 / (4 * np.pi * EPS0)
a0 = 4 * np.pi * EPS0 * HBAR**2 / (ME * QE**2)
ER = ME * QE**4 / (2 * (4 * np.pi * EPS0)**2 * HBAR**2) / EV
Rinf = ER * EV / (H * C)
close(a0, PC["Bohr radius"][0], 1e-8, "a0")
close(Rinf, PC["Rydberg constant"][0], 1e-8, "Rinf")
v1 = QE**2 * K / HBAR
print(f"Ex 4.3: a0 = {a0*1e9:.6f} nm, E1 = {-ER:.4f} eV, v1/c = {v1/C:.6f} (1/{C/v1:.2f}), R = {Rinf:.6e} 1/m")


def line(ni, nf, R=Rinf, Z=1):
    """wavelength (nm) of the hydrogen-like transition ni -> nf, 1/lam = Z^2 R (1/nf^2 - 1/ni^2)"""
    return 1e9 / (Z**2 * R * (1 / nf**2 - 1 / ni**2))


RH = Rinf / (1 + ME / MP)
RD = Rinf / (1 + ME / M_D)
Ha_inf, Ha_H, Ha_D = line(3, 2), line(3, 2, RH), line(3, 2, RD)
N_AIR = 1.000277
print(f"H-alpha: {Ha_inf:.2f} nm (infinite mass), {Ha_H:.2f} nm (reduced mass, vacuum), "
      f"{Ha_H/N_AIR:.2f} nm (air); D-alpha {Ha_D:.2f} nm, isotope shift {Ha_H-Ha_D:.3f} nm")
assert abs(Ha_inf - 656.1) < 0.05 and abs(Ha_H / N_AIR - 656.28) < 0.03
balmer = [line(n, 2) for n in range(3, 8)]
print("Balmer (nm):", np.round(balmer, 1), " limits: Lyman %.1f, Balmer %.1f, Paschen %.1f" %
      (1e9 / Rinf, 4e9 / Rinf, 9e9 / Rinf))
close(balmer[1], 486.0, 1e-3, "H-beta")

fig, axs = plt.subplots(2, 1, figsize=(8, 5), gridspec_kw={"height_ratios": [2, 1]})
for n in range(1, 8):
    axs[0].hlines(-ER / n**2, 0, 1, color="k", lw=1)
    axs[0].text(1.02, -ER / n**2, f"n={n}", va="center", fontsize=8)
for nf, col, name in [(1, "C4", "Lyman"), (2, "C0", "Balmer"), (3, "C3", "Paschen")]:
    for ni in range(nf + 1, nf + 8):
        lam = line(ni, nf)
        axs[1].vlines(lam, 0, 1, color=col, lw=1)
        x = 0.15 * nf + 0.02 * (ni - nf)
        axs[0].annotate("", xy=(x, -ER / nf**2), xytext=(x, -ER / ni**2),
                        arrowprops=dict(arrowstyle="->", color=col, lw=0.8))
    axs[1].text(line(nf + 1, nf), 1.05, name, color=col, fontsize=8)
axs[0].set_xlim(0, 1.15); axs[0].set_xticks([]); axs[0].set_ylabel("E (eV)")
axs[0].set_title("Bohr levels and the three series")
axs[1].set_xscale("log"); axs[1].set_yticks([]); axs[1].set_ylim(0, 1.3)
axs[1].axvspan(380, 750, color="gold", alpha=0.2)
axs[1].set_xlabel("wavelength (nm), visible band shaded")
fig.tight_layout(); fig.savefig(FIG / "bohr_series.png", dpi=120); plt.close(fig)

# classical atom (EM-8 Example 4.4): t_f = a0^3 / (4 r_e^2 c)
r_e = K * QE**2 / (ME * C**2)
t_f = a0**3 / (4 * r_e**2 * C)
print(f"classical collapse time from a0: {t_f:.3e} s")
assert 1.5e-11 < t_f < 1.6e-11

# correspondence principle (C3): orbital frequency vs transition frequency n -> n-1
for n in (10, 100, 1000):
    f_orb = 2 * ER / n**3
    f_tr = ER * (1 / (n - 1)**2 - 1 / n**2)
    print(f"  n = {n}: f(n->n-1)/f_orbit = {f_tr/f_orb:.5f}")
assert abs(ER * (1 / 81 - 1 / 100) / (2 * ER / 1000) - 1.17284) < 1e-4

# Compton (Example 4.4): Mo K-alpha, 71.07 pm, at 90 degrees
lamC = H / (ME * C)


def compton(lam, theta):
    """scattered wavelength lam' = lam + (h/mc)(1 - cos theta), notes (3.7)"""
    return lam + lamC * (1 - np.cos(theta))


lamX = 71.07e-12
lamXp = compton(lamX, np.pi / 2)
Ein, Eout = HC_EVNM / (lamX * 1e9), HC_EVNM / (lamXp * 1e9)
print(f"Ex 4.4: lambda_C = {lamC*1e12:.4f} pm, lam' = {lamXp*1e12:.3f} pm, E = {Ein/1e3:.3f} -> "
      f"{Eout/1e3:.3f} keV, electron KE = {(Ein-Eout):.0f} eV, fractional shift {lamC/lamX*100:.2f}%")
close(lamC, PC["Compton wavelength"][0], 1e-8, "Compton wavelength")
# sympy check of the Compton derivation: energy + momentum conservation gives (3.7)
lm, lmp, th, m_, c_, h_ = sp.symbols("lambda lambda_p theta m c h", positive=True)
p, pp = h_ / lm, h_ / lmp
pe2 = p**2 + pp**2 - 2 * p * pp * sp.cos(th)                 # electron momentum^2 (law of cosines)
Ee = h_ * c_ / lm - h_ * c_ / lmp + m_ * c_**2              # electron energy
cons = sp.expand(Ee**2 - pe2 * c_**2 - m_**2 * c_**4)       # E^2 = p^2c^2 + m^2c^4 (EM-8 (3.25))
sol = sp.solve(cons, lmp)
assert any(sp.simplify(s - (lm + h_ / (m_ * c_) * (1 - sp.cos(th)))) == 0 for s in sol)

# de Broglie (Example 4.5)
lam_e100 = H / np.sqrt(2 * ME * 100 * EV)
print(f"Ex 4.5: p(100 eV) = {np.sqrt(2*ME*100*EV):.3e} kg m/s")
lam_ball = H / (0.145 * 40.0)
lam_e54 = H / np.sqrt(2 * ME * 54 * EV)
dg = 0.215e-9 * np.sin(np.radians(50))
print(f"Ex 4.5: electron 100 eV lambda = {lam_e100*1e9:.4f} nm, baseball = {lam_ball:.3e} m; "
      f"Davisson-Germer: h/p(54 eV) = {lam_e54*1e9:.4f} nm vs d sin50 = {dg*1e9:.4f} nm")
assert abs(lam_e100 * 1e9 - 0.1226) < 1e-4 and abs(lam_e54 / dg - 1) < 0.02
# relativistic correction at 100 eV is tiny
p_rel = np.sqrt((100 * EV + ME * C**2)**2 - (ME * C**2)**2) / C
assert abs(H / p_rel / lam_e100 - 1) < 1e-4

# Stern-Gerlach deflection estimate (3.6): F = mu_B dB/dz on a silver atom
m_Ag = 107.87 * AMU
F_sg = MU_B * 1e3
t_sg = 0.035 / 500
z_sg = 0.5 * F_sg / m_Ag * t_sg**2
print(f"SG: F = {F_sg:.3e} N, a = {F_sg/m_Ag:.3e} m/s^2, deflection = +-{z_sg*1e3:.3f} mm")
assert 0.12e-3 < z_sg < 0.135e-3

# Bohr-Sommerfeld for the oscillator: J = oint p dq = 2 pi E / omega (CM-8 (3.24))
q_, E_, w_, mm = sp.symbols("q E omega m", positive=True)
A_ = sp.sqrt(2 * E_ / (mm * w_**2))
Jint = 2 * sp.integrate(sp.sqrt(2 * mm * (E_ - mm * w_**2 * q_**2 / 2)), (q_, -A_, A_))
assert sp.simplify(Jint - 2 * sp.pi * E_ / w_) == 0

# %% [markdown]
# ## Part 3: A Stern-Gerlach simulator (notes 3.7-3.8, Example 4.6)
# hbar = 1 from here on. State: a unit 2-vector in the S_z basis. Measuring spin along n:
# outcomes +-1/2, probabilities |<n,+-|psi>|^2 (Born, postulate 3), then the state becomes |n,+-> (postulate 4).

# %%
SX = np.array([[0, 1], [1, 0]], complex) / 2
SY = np.array([[0, -1j], [1j, 0]]) / 2
SZ = np.array([[1, 0], [0, -1]], complex) / 2
UP, DN = np.array([1, 0], complex), np.array([0, 1], complex)


def axis(theta, phi=0.0):
    """unit vector n = (sin th cos ph, sin th sin ph, cos th)"""
    return np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def n_dot_S(n):
    """the observable S_n = n . S = (1/2) n . sigma (LA-8 (3.15))"""
    return n[0] * SX + n[1] * SY + n[2] * SZ


def eigvecs(n):
    """eigenvalues (+1/2, -1/2) and orthonormal eigenvectors of n . S, via eigh (LA-8 Theorem 3.12)"""
    w, U = np.linalg.eigh(n_dot_S(n))
    return w[::-1], U[:, ::-1]


def measure(psi, n, rng):
    """Measure S_n on each row of psi (N x 2). Return outcomes (+1 or -1) and collapsed states.
    P(+) = |<n,+|psi>|^2 (postulate 3); state -> |n,+-> (postulate 4)."""
    _, U = eigvecs(n)
    amp = psi @ U.conj()                  # amp[:, k] = <n,k|psi>
    p_up = np.abs(amp[:, 0])**2
    up = rng.random(len(psi)) < p_up
    out = np.where(up, 1, -1)
    new = np.where(up[:, None], U[:, 0][None, :], U[:, 1][None, :])
    return out, new


# analytic checks: eigenvectors of S_x in the S_z basis, P = 1/2, cos^2(theta/2)
w, U = eigvecs(axis(np.pi / 2))
assert np.allclose(w, [0.5, -0.5])
assert np.isclose(abs(U[:, 0] @ np.array([1, 1]) / np.sqrt(2)), 1)
for th in np.linspace(0, np.pi, 7):
    assert np.isclose(abs(eigvecs(axis(th))[1][:, 0].conj() @ UP)**2, np.cos(th / 2)**2)

N = 10_000
Z, X = axis(0), axis(np.pi / 2)
src = np.tile(UP, (N, 1))                     # atoms already filtered to S_z = +1/2


def check(count, n_tot, p, label):
    """binomial check: count within 4 sigma of n_tot p"""
    sig = np.sqrt(n_tot * p * (1 - p)) if 0 < p < 1 else 0.5
    assert abs(count - n_tot * p) < 4 * sig + 1e-9, f"{label}: {count} vs {n_tot*p}"
    print(f"  {label}: {count}/{n_tot} = {count/n_tot:.4f} (predicted {p:.4f})")


print("Ex 4.6, 10^4 atoms:")
o1, s1 = measure(src, Z, rng)
check((o1 == 1).sum(), N, 1.0, "z+ then z: +")
o2, s2 = measure(src, X, rng)
check((o2 == 1).sum(), N, 0.5, "z+ then x: +")
keep = o2 == 1                                 # block the x- beam
o3, _ = measure(s2[keep], Z, rng)
check((o3 == -1).sum(), keep.sum(), 0.5, "z+ then x+ then z: -")
print(f"  binomial sigmas: {np.sqrt(N*0.25):.1f}, {np.sqrt(keep.sum()*0.25):.1f}")
th60 = axis(np.pi / 3)
o4, s4 = measure(src, th60, rng)
check((o4 == 1).sum(), N, np.cos(np.pi / 6)**2, "z+ then 60deg: +")
o5, _ = measure(s4, Z, rng)
check(((o4 == 1) & (o5 == -1)).sum(), N, np.cos(np.pi / 6)**2 * np.sin(np.pi / 6)**2,
      "z+ -> 60deg+ -> z-  (B5)")
check((o5 == -1).sum(), N, 0.375, "z+ -> 60deg (either) -> z-  (B5)")
# recombined (no measurement at the 60-degree analyzer): both paths kept coherently, so the
# analyzer acts as sum_k |n,k><n,k| = I (completeness, LA-8 3.12) and the state is untouched
_, U60 = eigvecs(th60)
recomb = src @ (U60 @ U60.conj().T).T
assert np.allclose(recomb, src)
o6, _ = measure(recomb, Z, rng)
check((o6 == -1).sum(), N, 0.0, "z+ -> 60deg recombined -> z-")

# %% [markdown]
# ## Part 4: The Bloch sphere (notes 3.9)
# |n,+> = cos(th/2)|up> + e^{i ph} sin(th/2)|down>; its Bloch vector <sigma> = 2<S> equals n.

# %%
def bloch(psi):
    """Bloch vector (<sigma_x>, <sigma_y>, <sigma_z>) = 2 <psi|S|psi>"""
    return 2 * np.real([psi.conj() @ S @ psi for S in (SX, SY, SZ)])


fig = plt.figure(figsize=(5.5, 5.5))
ax = fig.add_subplot(projection="3d")
uu, vv = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
ax.plot_wireframe(np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv), color="0.85", lw=0.5)
for th, ph in [(0, 0), (np.pi / 2, 0), (np.pi / 2, np.pi / 2), (np.pi / 3, np.pi / 4), (2.2, 4.0)]:
    n = axis(th, ph)
    _, U = eigvecs(n)
    explicit = np.array([np.cos(th / 2), np.exp(1j * ph) * np.sin(th / 2)])
    assert np.isclose(abs(explicit.conj() @ U[:, 0]), 1)      # same ray as eigh's vector
    assert np.allclose(bloch(U[:, 0]), n) and np.allclose(bloch(U[:, 1]), -n)
    for v, mk in [(n, "o"), (-n, "x")]:
        ax.plot(*np.c_[np.zeros(3), v], color="C0" if mk == "o" else "C3", lw=1)
        ax.scatter(*v, color="C0" if mk == "o" else "C3", marker=mk)
ax.set_title(r"eigenvectors of $\hat n\cdot\sigma$: $|n,+\rangle$ (o) at $\hat n$, $|n,-\rangle$ (x) at $-\hat n$")
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
fig.tight_layout(); fig.savefig(FIG / "bloch.png", dpi=120); plt.close(fig)
# global phase is invisible; orthogonal kets sit at antipodes
psi = rng.normal(size=2) + 1j * rng.normal(size=2); psi /= np.linalg.norm(psi)
assert np.allclose(bloch(psi), bloch(np.exp(0.7j) * psi))
assert np.isclose(np.linalg.norm(bloch(psi)), 1)

# %% [markdown]
# ## Part 5: The double slit with which-path tags (notes 1.7, 3.10)
# Far-field amplitudes phi_{1,2}(x) = sinc(pi a x/(lam L)) e^{-+ i pi d x/(lam L)}. With tags |t1>, |t2>,
# P(x) = |phi1|^2 + |phi2|^2 + 2 Re(phi1* phi2 <t1|t2>)  (3.16). Visibility = |<t1|t2>|.

# %%
lam_ds, d_ds, a_ds, L_ds = 1.0, 5.0, 1.0, 100.0          # arbitrary units; fringe spacing lam L / d = 20
xs = np.linspace(-90, 90, 3601)          # inside the first envelope zero at |x| = 100


def amps(x):
    """far-field slit amplitudes phi_1, phi_2 (normalized so |phi_j(0)| = 1/sqrt2)"""
    env = np.sinc(a_ds * x / (lam_ds * L_ds)) / np.sqrt(2)
    ph = np.pi * d_ds * x / (lam_ds * L_ds)
    return env * np.exp(-1j * ph), env * np.exp(1j * ph)


def pattern(x, overlap):
    """detection probability density with tag overlap <t1|t2> = overlap, eq. (3.16)"""
    p1, p2 = amps(x)
    return np.abs(p1)**2 + np.abs(p2)**2 + 2 * np.real(np.conj(p1) * p2 * overlap)


def visibility(x, P):
    """(Pmax - Pmin)/(Pmax + Pmin) of the pattern divided by the single-slit envelope |phi1|^2 + |phi2|^2"""
    p1, p2 = amps(x)
    Pn = P / (np.abs(p1)**2 + np.abs(p2)**2)
    return (Pn.max() - Pn.min()) / (Pn.max() + Pn.min())


for ov in (1.0, 0.5, 0.0):
    V_ = visibility(xs, pattern(xs, ov))
    print(f"  tag overlap {ov}: fringe visibility {V_:.3f}")
    assert abs(V_ - ov) < 1e-3
# tagged state built explicitly in C^(screen) x C^2: tracing out the tag gives the same pattern
t1, t2 = np.array([1, 0]), np.array([np.cos(np.pi / 3), np.sin(np.pi / 3)])   # overlap 0.5
p1, p2 = amps(xs)
joint = p1[:, None] * t1[None, :] + p2[:, None] * t2[None, :]
assert np.allclose((np.abs(joint)**2).sum(1), pattern(xs, 0.5))
# electrons one at a time: sample 3000 hits from each pattern
fig, axs = plt.subplots(3, 1, figsize=(7, 6), sharex=True)
for axx, ov, title in zip(axs, (1.0, 0.0, None), ("amplitudes add (no tag)",
                                                  "orthogonal which-path tags", "classical: P1 + P2")):
    P = pattern(xs, ov) if ov is not None else np.abs(p1)**2 + np.abs(p2)**2
    cdf = np.cumsum(P); cdf /= cdf[-1]
    hits = np.interp(rng.random(3000), cdf, xs)
    axx.hist(hits, bins=120, density=True, color="0.75")
    axx.plot(xs, P / np.trapezoid(P, xs) if hasattr(np, "trapezoid") else P / np.trapz(P, xs), "k")
    axx.set_title(title, fontsize=9)
axs[-1].set_xlabel("screen position (units of lambda L / d = 20)")
fig.tight_layout(); fig.savefig(FIG / "double_slit.png", dpi=120); plt.close(fig)
assert np.allclose(pattern(xs, 0.0), np.abs(p1)**2 + np.abs(p2)**2)

# %% [markdown]
# ## Part 6: Commutators, uncertainty, unitary evolution (notes 3.7, 3.9, 3.11)
# [S_x, S_y] = i S_z (hbar = 1). Robertson: Var(A) Var(B) >= (1/4)|<[A,B]>|^2 (3.17).
# Time evolution |psi(t)> = exp(-i H t)|psi(0)> is unitary because H is Hermitian (3.11, Theorem 3.3).

# %%
for A, B, Cc in [(SX, SY, SZ), (SY, SZ, SX), (SZ, SX, SY)]:
    assert np.allclose(A @ B - B @ A, 1j * Cc)
s = sp.Rational(1, 2) * sp.Matrix([[0, 1], [1, 0]]), sp.Rational(1, 2) * sp.Matrix([[0, -sp.I], [sp.I, 0]])
assert s[0] * s[1] - s[1] * s[0] == sp.I * sp.Rational(1, 2) * sp.Matrix([[1, 0], [0, -1]])


def var(Op, psi):
    """Var(Op) = <Op^2> - <Op>^2"""
    m = np.real(psi.conj() @ Op @ psi)
    return np.real(psi.conj() @ Op @ Op @ psi) - m**2


margins = []
for _ in range(2000):
    psi = rng.normal(size=2) + 1j * rng.normal(size=2); psi /= np.linalg.norm(psi)
    lhs = var(SX, psi) * var(SY, psi)
    rhs = 0.25 * abs(psi.conj() @ (SX @ SY - SY @ SX) @ psi)**2
    nx, ny, nz = bloch(psi)
    assert np.isclose(lhs, (1 - nx**2) * (1 - ny**2) / 16) and np.isclose(rhs, nz**2 / 16)
    assert lhs >= rhs - 1e-12
    margins.append(lhs - rhs)
print(f"uncertainty: 2000 random states, min margin {min(margins):.2e} (>= 0), "
      f"|up>: product {var(SX, UP)*var(SY, UP):.4f} = bound {0.25*0.25:.4f}")
# expectation value from Born: sum a_i p_i = <psi|A|psi>
w, U = eigvecs(axis(1.0, 2.0))
psi = (UP + 2j * DN) / np.sqrt(5)
close(np.sum(w * np.abs(U.conj().T @ psi)**2), np.real(psi.conj() @ n_dot_S(axis(1.0, 2.0)) @ psi), 1e-12, "<A>")

# Larmor precession: H = w0 S_z, start in |+x>
w0 = 2.0
plus_x = (UP + DN) / np.sqrt(2)
for t in np.linspace(0, 5, 11):
    Ut = expm(-1j * w0 * SZ * t)
    assert np.allclose(Ut.conj().T @ Ut, np.eye(2))
    pt = Ut @ plus_x
    assert np.allclose(bloch(pt), [np.cos(w0 * t), np.sin(w0 * t), 0])
    assert np.isclose(abs(plus_x.conj() @ pt)**2, np.cos(w0 * t / 2)**2)
g_e = -PC["electron g factor"][0]
w_e = g_e * MU_B * 1.0 / HBAR
print(f"B7: electron in 1 T: omega0 = {w_e:.4e} rad/s, period = {2*np.pi/w_e*1e12:.2f} ps")

# %% [markdown]
# ## Your turn
# 1. **Compton shift vs angle.** Plot lam'(theta) for Mo K-alpha, 0 to 180 degrees, and the electron's
#    kinetic energy. (Reference below.) Where is the shift half its maximum?
# 2. **Hydrogen-like ions.** Redo Part 2 for He+ (Z = 2): E_n = -Z^2 E_R/n^2, r_n = n^2 a0/Z. Show that the
#    He+ n = 4 -> 3 line sits at 468.6 nm and that every *even*-upper-level He+ line n = 2k -> 4 coincides
#    (infinite nuclear mass) with a Balmer line k -> 2. That coincidence is why Pickering's lines were first
#    misassigned to hydrogen; the reduced-mass difference separates them.
# 3. **Polarizers as Stern-Gerlach analyzers (EM-7 3.6).** A Jones vector is a 2-vector; a polarizer at angle
#    a projects onto (cos a, sin a). Rewrite `measure` for photons and run 0, 45, 90 degree polarizers:
#    1/8 of unpolarized light gets through, and 0 without the middle one. Then put the angle-doubling
#    dictionary on paper: polarizer angle a <-> spin axis at theta = 2a.
# 4. **Exercise C2** (hidden-variable toy model) and **B8** (partial tags) extend Parts 3 and 5.

# %%
ths = np.linspace(0, np.pi, 181)
lp = compton(lamX, ths)
KE = HC_EVNM / (lamX * 1e9) - HC_EVNM / (lp * 1e9)
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(np.degrees(ths), (lp - lamX) * 1e12, "k")
ax.set_xlabel("scattering angle (deg)"); ax.set_ylabel(r"$\Delta\lambda$ (pm)")
ax2 = ax.twinx(); ax2.plot(np.degrees(ths), KE / 1e3, "C3--"); ax2.set_ylabel("electron KE (keV)", color="C3")
fig.tight_layout(); fig.savefig(FIG / "compton.png", dpi=120); plt.close(fig)
assert np.isclose((lp[-1] - lamX), 2 * lamC)
print(f"max Compton shift {2*lamC*1e12:.3f} pm; fractional shift at 500 nm {lamC/500e-9:.2e}")
# He+
He43 = line(4, 3, Z=2)
print(f"He+ 4->3 = {He43:.2f} nm; He+ 6->4 = {line(6, 4, Z=2):.3f} nm = H 3->2 {line(3, 2):.3f} nm")
assert abs(He43 - 468.6) < 0.1
for k in range(3, 7):
    assert np.isclose(line(2 * k, 4, Z=2), line(k, 2))
# three polarizers
def polarize(J, a):
    """project Jones vectors (rows of J) onto the axis at angle a; return transmitted intensity, new state"""
    e = np.array([np.cos(a), np.sin(a)])
    amp = J @ e
    return np.abs(amp)**2, np.outer(amp, e)


angles = rng.uniform(0, np.pi, 20000)
J0 = np.c_[np.cos(angles), np.sin(angles)]         # unpolarized: random linear polarizations
I1, J1 = polarize(J0, 0)
I2, J2 = polarize(J1, np.pi / 4)
I3, _ = polarize(J2, np.pi / 2)
I3b, _ = polarize(J1, np.pi / 2)
print(f"polarizers 0,45,90: {I3.mean():.4f} (1/8 = 0.125); 0,90: {I3b.mean():.2e}")
assert abs(I3.mean() - 0.125) < 0.005 and I3b.mean() < 1e-30

# %% [markdown]
# ## Part 8: Numbers quoted in solutions.md

# %%
print("Solutions:")
for lam_nm in (400, 700):
    print(f"  A1: {lam_nm} nm photon = {HC_EVNM/lam_nm:.3f} eV")
E_radio = H * 1e6 / EV
print(f"  A1: 1 MHz photon = {E_radio:.3e} eV; sodium threshold {lam_thr:.1f} nm")
print(f"  A2: E2 = {-ER/4:.3f} eV, E3 = {-ER/9:.3f} eV, Lyman-alpha = {line(2, 1):.2f} nm, "
      f"Hbeta = {line(4, 2):.2f} nm")
lam_n = H / np.sqrt(2 * MN * KB * 300)
print(f"  A5: kT(300 K) = {KB*300/EV:.5f} eV, neutron lambda = {lam_n*1e9:.4f} nm")
E100 = 100e3
lam100 = HC_EVNM / E100 * 1e3                      # pm
lam100p = lam100 + 2 * lamC * 1e12
print(f"  B1: 100 keV: lam = {lam100:.3f} pm, back-scattered {lam100p:.3f} pm, "
      f"E' = {HC_EVNM/(lam100p*1e-3)/1e3:.2f} keV, max KE = {E100/1e3 - HC_EVNM/(lam100p*1e-3)/1e3:.2f} keV")
assert abs(E100 / 1e3 - HC_EVNM / (lam100p * 1e-3) / 1e3 - 28.1) < 0.1
print(f"  B2: R_H = {RH:.6e}, R_D = {RD:.6e} 1/m; H-alpha {Ha_H:.3f}, D-alpha {Ha_D:.3f} nm, "
      f"shift {Ha_H-Ha_D:.3f} nm")
assert abs(Ha_H - Ha_D - 0.179) < 0.002
print(f"  C1: x_peak = {xw:.6f}, b = {b_wien*1e3:.4f} mm K, peak(5000 K) = {b_wien/T5*1e9:.1f} nm, "
      f"Sun 5772 K -> {b_wien/5772*1e9:.1f} nm")
# C2: hidden-variable model, P(+n | z+) = 1 - theta/pi vs quantum cos^2(theta/2)
mvec = rng.normal(size=(200000, 3)); mvec /= np.linalg.norm(mvec, axis=1)[:, None]
mz = mvec[mvec[:, 2] > 0]
for th in (np.pi / 3, np.pi / 2, 2 * np.pi / 3):
    p_hv = (mz @ axis(th) > 0).mean()
    print(f"  C2: theta = {np.degrees(th):.0f} deg: hidden-variable {p_hv:.4f} (1 - th/pi = {1-th/np.pi:.4f}), "
          f"quantum {np.cos(th/2)**2:.4f}")
    assert abs(p_hv - (1 - th / np.pi)) < 0.005
# B4(c), B9: psi = (|up> + 2i|down>)/sqrt5
psi9 = (UP + 2j * DN) / np.sqrt(5)
n9 = bloch(psi9)
th9 = 2 * np.arccos(1 / np.sqrt(5))
assert np.allclose(n9, [0, 0.8, -0.6]) and np.allclose(n9, axis(th9, np.pi / 2))
plus_y, minus_y = (UP + 1j * DN) / np.sqrt(2), (UP - 1j * DN) / np.sqrt(2)
py = abs(np.vdot(plus_y, psi9))**2, abs(np.vdot(minus_y, psi9))**2
assert np.allclose(py, [0.9, 0.1]) and np.isclose(0.5 * (py[0] - py[1]), 0.4)
print(f"  B4/B9: Bloch vector {np.round(n9, 4)}, theta = {np.degrees(th9):.2f} deg, phi = 90 deg; "
      f"P(+y) = {py[0]:.2f}, P(-y) = {py[1]:.2f}, <S_y> = {0.5*(py[0]-py[1]):.2f} hbar")
for beta in (0, np.pi / 6, np.pi / 3, np.pi / 2):          # B8(b) reference
    t1, t2 = np.array([1, 0]), np.array([np.cos(beta), np.exp(0.4j) * np.sin(beta)])
    ov = np.vdot(t1, t2)
    p1, p2 = amps(xs)
    joint = p1[:, None] * t1[None, :] + p2[:, None] * t2[None, :]
    P = (np.abs(joint)**2).sum(1)
    assert np.allclose(P, pattern(xs, ov))
    assert abs(visibility(xs, P) - abs(np.cos(beta))) < 1e-3
print(f"  B8: visibilities |cos beta| = {[float(round(abs(np.cos(b_)), 3)) for b_ in (0, np.pi/6, np.pi/3, np.pi/2)]}")
print(f"  A1: photons needed per electron at 1 MHz: {W_NA/E_radio:.2e}; A5: p_n = {np.sqrt(2*MN*KB*300):.4e} kg m/s")
print(f"  B1: E -> inf back-scatter limit m c^2/2 = {ME*C**2/EV/2e3:.1f} keV")
print(f"  B2: m_p/m = {MP/ME:.2f}, m_d/m = {M_D/ME:.2f}")
print(f"  B5: sigmas {np.sqrt(N*0.1875*0.8125):.1f}, {np.sqrt(N*0.375*0.625):.1f}")
print(f"  B7: precession frequency {w_e/(2*np.pi)/1e9:.2f} GHz")
a_rad = 8 * np.pi**5 * KB**4 / (15 * H**3 * C**3)
print(f"  C1: radiation constant a = {a_rad:.4e} J m^-3 K^-4")
print(f"  C3: n = 10 ratio {ER*(1/81-1/100)/(2*ER/1000):.5f}, n = 100 ratio {(1/99**2-1/100**2)/(2/100**3):.5f}")

# %%
if __name__ == "__main__":
    dt = time.time() - T0
    print(f"All checks passed in {dt:.1f} s.")
    assert dt < 90
