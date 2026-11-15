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


def pin(value, quoted, half_unit, label=""):
    """assert that a computed value rounds to the number quoted in notes.md / solutions.md:
    |value - quoted| <= half a unit in the last quoted digit. If a pin fails, the prose is stale."""
    v, q = np.asarray(value, float), np.asarray(quoted, float)
    assert np.all(np.abs(v - q) <= half_unit * (1 + 1e-9)), f"{label}: computed {value}, prose quotes {quoted}"

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
# the numbers of Example 4.1, as quoted
pin(x0, 5.7551, 5e-5, "x"); pin(8 * np.pi * H * C / lam0**5, 1.5976e8, 5e3, "8 pi h c/lam^5")
pin(np.expm1(x0), 314.8, 0.05, "e^x - 1"); pin(KB * T5, 6.903e-20, 5e-24, "k T")
pin(H * C, 1.9864e-25, 5e-30, "h c in J m"); pin(KB / EV, 8.61733e-5, 5e-11, "k_B in eV/K")
pin(KB * T5 / EV, 0.43087, 5e-6, "k T in eV"); pin(HC_EVNM, 1239.84, 5e-3, "h c in eV nm")
pin(up * 1e-9, 5.07e-4, 5e-7, "Planck u per nm")
pin(up, 5.075e5, 50, "Planck u"); pin(ur, 2.776e7, 5e3, "RJ u"); pin(ur / up, 54.70, 0.005, "RJ/Planck")
# long-wavelength agreement: at 50 um (x = 0.058) the ratio is within 3%
pin(H * C / (50e-6 * KB * T5), 0.058, 5e-4, "x at 50 um")
assert abs(u_rj(50e-6, T5) / u_planck(50e-6, T5) - 1) < 0.03
# Wien is the x >> 1 limit of Planck: Planck/Wien = 1/(1 - e^-x)
close(up / u_wien(lam0, T5), 1 / (1 - np.exp(-x0)), 1e-12, "Wien limit")
lam_peak = minimize_scalar(lambda L: -u_planck(L * 1e-9, T5), bounds=(300, 900), method="bounded",
                           options={"xatol": 1e-6}).x * 1e-9
xw = brentq(lambda x: x - 5 * (1 - np.exp(-x)), 1, 10)
b_wien = H * C / (KB * xw)
print(f"Wien peak at 5000 K: {b_wien/T5*1e9:.1f} nm (x = {xw:.6f}, b = {b_wien:.4e} m K)")
close(lam_peak, b_wien / T5, 1e-4, "numerical peak")
pin(xw, 4.965114, 5e-7, "Wien x"); pin(b_wien, 2.8978e-3, 5e-8, "Wien b"); pin(b_wien / T5 * 1e9, 579.6, 0.05, "peak")
pin(b_wien * 1e3, 2.898, 5e-4, "lambda_max T in mm K (notes 3.1)")
# energy in Planck spectrum is finite: integral = (8 pi^5 k^4 / 15 h^3 c^3) T^4
lg = np.linspace(20e-9, 200e-6, 400001)
Uint = np.trapezoid(u_planck(lg, T5), lg) if hasattr(np, "trapezoid") else np.trapz(u_planck(lg, T5), lg)
close(Uint, 8 * np.pi**5 * KB**4 / (15 * H**3 * C**3) * T5**4, 1e-6, "Stefan-Boltzmann")

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
pin(E400, 3.0996, 5e-5, "E(400 nm)"); pin(V400, 0.7996, 5e-5, "V_stop"); pin(lam_thr, 539.1, 0.05, "threshold")
pin(W_NA * EV / H, 5.56e14, 5e11, "f0"); pin(np.sqrt(2 * V400 * EV / ME), 5.30e5, 5e2, "v_max")
# synthetic Millikan-style data: stopping potential vs frequency with 0.02 V noise
lams = np.array([250, 280, 313, 334, 365, 405, 436, 480, 510]) * 1e-9
assert lams.size == 9 and np.all(HC_EVNM / (lams * 1e9) > W_NA)      # nine wavelengths, all above threshold
f = C / lams
V = (H * f / EV - W_NA) + rng.normal(0, 0.02, f.size)
slope, icpt = np.polyfit(f, V, 1)
h_fit = slope * EV
print(f"fitted h = {h_fit:.4e} J s (true {H:.4e}, off by {abs(h_fit/H-1)*100:.2f}%), fitted W = {-icpt:.3f} eV")
close(h_fit, H, 0.02, "fitted h")
assert abs(-icpt - W_NA) < 0.08
# seeded values quoted in Example 4.2 (seed 0, NumPy's default PCG64 stream): if these fail, the
# data generation changed and the prose must be updated to match
pin(h_fit, 6.611e-34, 5e-38, "fitted h, Example 4.2"); pin(-icpt, 2.287, 5e-4, "fitted W, Example 4.2")
assert abs(h_fit / H - 1) < 0.003
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(f / 1e14, V, "ko", label="synthetic data")
ff = np.linspace(W_NA * EV / H, f.max(), 50)
ax.plot(ff / 1e14, slope * ff + icpt, "C0", label=f"fit: h = {h_fit:.3e} J s")
ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlabel("frequency (10$^{14}$ Hz)"); ax.set_ylabel("stopping potential (V)")
ax.set_title("Photoelectric effect, sodium"); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "photoelectric.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 2: Bohr, Compton, de Broglie, Stern-Gerlach numbers (notes 3.3-3.6, Examples 4.3-4.5)
# Bohr (3.4): r_n = n^2 a0, a0 = 4 pi eps0 hbar^2/(m e^2); E_n = -E_R/n^2, E_R = m e^4/(2 (4 pi eps0)^2 hbar^2).
# Rydberg (3.5): 1/lam = R (1/nf^2 - 1/ni^2), R = E_R/(h c). Stern-Gerlach force (3.6): F_z = mu_z dB_z/dz.
# Compton (3.7): dlam = (h/mc)(1 - cos theta). de Broglie (3.8): lam = h/p.

# %%
K = 1 / (4 * np.pi * EPS0)
a0 = 4 * np.pi * EPS0 * HBAR**2 / (ME * QE**2)
ER = ME * QE**4 / (2 * (4 * np.pi * EPS0)**2 * HBAR**2) / EV
Rinf = ER * EV / (H * C)
close(a0, PC["Bohr radius"][0], 1e-8, "a0")
close(Rinf, PC["Rydberg constant"][0], 1e-8, "Rinf")
v1 = QE**2 * K / HBAR
print(f"Ex 4.3: a0 = {a0*1e9:.6f} nm, E1 = {-ER:.4f} eV, v1/c = {v1/C:.6f} (1/{C/v1:.2f}), R = {Rinf:.6e} 1/m")
# the What-just-happened of notes 3.3: R = m e^4 / (8 eps0^2 h^3 c)
close(Rinf, ME * QE**4 / (8 * EPS0**2 * H**3 * C), 1e-12, "R = m e^4/(8 eps0^2 h^3 c)")
# Example 4.3 as quoted
pin(a0 * 1e9, 0.052918, 5e-7, "a0"); pin(ER, 13.606, 5e-4, "E_R"); pin(C / v1, 137.04, 5e-3, "1/alpha")
pin(Rinf, 1.097373e7, 5, "R_inf"); pin([-ER, -ER / 4, -ER / 9], [-13.606, -3.401, -1.512], 5e-4, "E_1, E_2, E_3")


def line(ni, nf, R=Rinf, Z=1):
    """wavelength (nm) of the hydrogen-like transition ni -> nf, 1/lam = Z^2 R (1/nf^2 - 1/ni^2)"""
    return 1e9 / (Z**2 * R * (1 / nf**2 - 1 / ni**2))


RH = Rinf / (1 + ME / MP)
RD = Rinf / (1 + ME / M_D)
Ha_inf, Ha_H, Ha_D = line(3, 2), line(3, 2, RH), line(3, 2, RD)
N_AIR = 1.000277
print(f"H-alpha: {Ha_inf:.2f} nm (infinite mass), {Ha_H:.2f} nm (reduced mass, vacuum), "
      f"{Ha_H/N_AIR:.2f} nm (air); D-alpha {Ha_D:.2f} nm, isotope shift {Ha_H-Ha_D:.3f} nm")
pin(Ha_inf, 656.11, 5e-3, "H-alpha, infinite mass"); pin(36e9 / (5 * Rinf), 656.11, 5e-3, "36/(5 R)")
pin(RH, 1.096776e7, 5, "R_H"); pin(Ha_H, 656.47, 5e-3, "H-alpha, reduced mass")
pin(Ha_H / N_AIR, 656.29, 5e-3, "H-alpha in air")
assert abs(Ha_H / N_AIR - 656.28) < 0.015          # the measured line in air, to 0.01 nm
pin(ME / MP, 5e-4, 5e-5, "reduced-mass correction to R (notes 3.3)")
balmer = [line(n, 2) for n in range(3, 8)]
print("Balmer (nm):", np.round(balmer, 1), " limits: Lyman %.1f, Balmer %.1f, Paschen %.1f" %
      (1e9 / Rinf, 4e9 / Rinf, 9e9 / Rinf))
pin(balmer, [656.1, 486.0, 433.9, 410.1, 396.9], 0.05, "Balmer series")
pin([1e9 / Rinf, 4e9 / Rinf, 9e9 / Rinf], [91.1, 364.5, 820.1], 0.05, "series limits")
pin(line(2, 1), 121.5, 0.05, "Lyman-alpha")
assert all(line(n, 1) < 380 for n in range(2, 60)) and all(line(n, 3) > 750 for n in range(4, 60))  # UV, IR

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
pin(t_f, 1.556e-11, 5e-15, "collapse time (EM-8 Example 4.4)"); pin(t_f, 1.56e-11, 5e-14, "collapse time (notes 1.3)")
pin(a0 * 1e9, 0.053, 5e-4, "starting radius (notes 1.3)")

# correspondence principle (C3): orbital frequency vs transition frequency n -> n-1
for n in (10, 100, 1000):
    f_orb = 2 * ER / n**3
    f_tr = ER * (1 / (n - 1)**2 - 1 / n**2)
    print(f"  n = {n}: f(n->n-1)/f_orbit = {f_tr/f_orb:.5f}")
    assert np.isclose(f_tr / f_orb, (2 * n - 1) * n / (2 * (n - 1)**2))       # C3(b) closed form
    assert abs(f_tr / f_orb - (1 + 1.5 / n)) < 3 / n**2                       # 1 + 3/(2n) + O(1/n^2)
pin(ER * (1 / 81 - 1 / 100) / (2 * ER / 1000), 1.17284, 5e-6, "C3, n = 10")
pin(ER * (1 / 99**2 - 1 / 100**2) / (2 * ER / 100**3), 1.01520, 5e-6, "C3, n = 100")

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
# Example 4.4 and the What-just-happened of notes 3.5, as quoted
pin(lamC * 1e12, 2.4263, 5e-5, "h/mc"); pin(lamXp * 1e12, 73.496, 5e-4, "lam'")
pin(6.626e-34 / (9.109e-31 * 2.998e8) * 1e12, 2.4263, 5e-5, "h/mc from the constants displayed in Example 4.4")
pin(lamC / lamX * 100, 3.41, 5e-3, "fractional shift at 90 deg, %")
pin([Ein / 1e3, Eout / 1e3], [17.445, 16.869], 5e-4, "photon energies, keV"); pin(Ein - Eout, 576, 0.5, "recoil KE, eV")
pin(2 * lamC * 1e12, 4.85, 5e-3, "maximum shift, pm"); pin(lamC / 500e-9, 5e-6, 5e-7, "fractional shift at 500 nm")
# energy conservation at 90 degrees with the relativistic electron: KE = sqrt(p_e^2 c^2 + m^2 c^4) - m c^2
pe2_num = (H / lamX)**2 + (H / lamXp)**2                       # law of cosines with cos(90 deg) = 0
close((np.sqrt(pe2_num * C**2 + (ME * C**2)**2) - ME * C**2) / EV, Ein - Eout, 1e-4, "Compton recoil energy")
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
# Example 4.5 as quoted
pin(np.sqrt(2 * ME * 100 * EV), 5.403e-24, 5e-28, "p(100 eV)"); pin(lam_e100 * 1e9, 0.1226, 5e-5, "lambda(100 eV)")
pin(lam_e54 * 1e9, 0.1669, 5e-5, "lambda(54 eV)"); pin(np.sin(np.radians(50)), 0.7660, 5e-5, "sin 50")
pin(dg * 1e9, 0.1647, 5e-5, "d sin 50"); pin((lam_e54 / dg - 1) * 100, 1.3, 0.05, "Davisson-Germer mismatch, %")
pin(lam_ball, 1.14e-34, 5e-37, "baseball"); pin(ME * C**2 / EV / 1e3, 511, 0.5, "m c^2, keV")
# relativistic correction at 100 eV is tiny
p_rel = np.sqrt((100 * EV + ME * C**2)**2 - (ME * C**2)**2) / C
assert abs(H / p_rel / lam_e100 - 1) < 1e-4
# notes 3.6: a whole number of de Broglie wavelengths around the orbit, 2 pi r = n h / p, is L = p r = n hbar
r_s, n_s, h_s, p_s = sp.symbols("r n h p", positive=True)
p_orbit = sp.solve(sp.Eq(2 * sp.pi * r_s, n_s * h_s / p_s), p_s)[0]
assert sp.simplify(p_orbit * r_s - n_s * h_s / (2 * sp.pi)) == 0

# Stern-Gerlach deflection estimate (3.6): F = mu_B dB/dz on a silver atom
m_Ag = 107.87 * AMU
F_sg = MU_B * 1e3
t_sg = 0.035 / 500
z_sg = 0.5 * F_sg / m_Ag * t_sg**2
print(f"SG: F = {F_sg:.3e} N, a = {F_sg/m_Ag:.3e} m/s^2, deflection = +-{z_sg*1e3:.3f} mm")
close(MU_B, QE * HBAR / (2 * ME), 1e-8, "mu_B = e hbar / 2m")
pin(MU_B, 9.274e-24, 5e-28, "mu_B"); pin(F_sg, 9.27e-21, 5e-24, "SG force"); pin(F_sg / m_Ag, 5.18e4, 50, "SG acceleration")
pin(t_sg, 7e-5, 5e-7, "time in the magnet"); pin(z_sg * 1e3, 0.127, 5e-4, "SG deflection, mm")
# classical prediction (notes 3.4): moments of random orientation have mu_z = mu cos(alpha) with cos(alpha)
# uniform on [-1, 1], so the deflections fill a flat band between -z_max and +z_max. A separate generator is
# used so that the seeded numbers quoted in the notes (which come from `rng`) do not move.
rng_cl = np.random.default_rng(1)
m_cl = rng_cl.normal(size=(100_000, 3)); m_cl /= np.linalg.norm(m_cl, axis=1)[:, None]
band, _ = np.histogram(m_cl[:, 2], bins=10, range=(-1, 1))
assert np.all(np.abs(band - 10_000) < 4 * np.sqrt(100_000 * 0.1 * 0.9)), band      # flat: no two-spot pattern
# the quantum result puts every atom at one of the two ends, where the classical band has almost nothing:
print(f"SG classical band: fraction of atoms within 5% of either end = {(np.abs(m_cl[:, 2]) > 0.95).mean():.3f}")
assert abs((np.abs(m_cl[:, 2]) > 0.95).mean() - 0.05) < 4 * np.sqrt(0.05 * 0.95 / 100_000)

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


# analytic checks: eigenvectors of S_x in the S_z basis (3.14), P = 1/2, cos^2(theta/2)
w, U = eigvecs(axis(np.pi / 2))
assert np.allclose(w, [0.5, -0.5])
assert np.isclose(abs(np.vdot(np.array([1, 1]) / np.sqrt(2), U[:, 0])), 1)        # |+x> = (1, 1)/sqrt2
assert np.isclose(abs(np.vdot(np.array([1, -1]) / np.sqrt(2), U[:, 1])), 1)       # |-x> = (1, -1)/sqrt2
_, Uy = eigvecs(axis(np.pi / 2, np.pi / 2))
assert np.isclose(abs(np.vdot(np.array([1, 1j]) / np.sqrt(2), Uy[:, 0])), 1)      # |+y> = (1, i)/sqrt2
assert np.isclose(abs(np.vdot(np.array([1, -1j]) / np.sqrt(2), Uy[:, 1])), 1)     # |-y> = (1, -i)/sqrt2
for th in np.linspace(0, np.pi, 7):
    assert np.isclose(abs(eigvecs(axis(th))[1][:, 0].conj() @ UP)**2, np.cos(th / 2)**2)
# the three bases are mutually unbiased (notes 3.8): every overlap between different axes has |.|^2 = 1/2
bases = [np.eye(2, dtype=complex), U, Uy]
for i_b in range(3):
    for j_b in range(3):
        if i_b != j_b:
            assert np.allclose(np.abs(bases[i_b].conj().T @ bases[j_b])**2, 0.5)
# Born probabilities sum to 1 for any state and any axis (Parseval, postulate 3), and <S_n> = sum a_i p_i
for th, ph, a_, b_ in [(0.3, 1.1, 0.6, 0.8j), (1.2, 5.5, 1, 0),
                       (2.0, 4.0, 1 / np.sqrt(3), np.sqrt(2 / 3) * np.exp(0.9j))]:
    psi_t = np.array([a_, b_], complex)
    wn, Un = eigvecs(axis(th, ph))
    probs = np.abs(Un.conj().T @ psi_t)**2
    assert np.isclose(probs.sum(), 1)
    assert np.isclose(wn @ probs, np.real(np.vdot(psi_t, n_dot_S(axis(th, ph)) @ psi_t)))      # (3.12)

N = 10_000
Z, X = axis(0), axis(np.pi / 2)
src = np.tile(UP, (N, 1))                     # atoms already filtered to S_z = +1/2


def check(count, n_tot, p, label, quoted=None):
    """Binomial check: count within 4 sigma of n_tot p. When the theory says p = 0 or p = 1 the count must be
    exact (a repeated measurement, or a recombined beam, has no statistical scatter). `quoted` is the count
    printed in notes.md / solutions.md for seed 0; it must be reproduced exactly, or the prose is stale."""
    if p in (0.0, 1.0):
        assert count == round(n_tot * p), f"{label}: {count} of {n_tot}, but the prediction is exact"
    else:
        sig = np.sqrt(n_tot * p * (1 - p))
        assert abs(count - n_tot * p) < 4 * sig, f"{label}: {count} vs {n_tot*p}"
    if quoted is not None:
        assert count == quoted, f"{label}: seed-0 count is {count}, prose quotes {quoted}"
    print(f"  {label}: {count}/{n_tot} = {count/n_tot:.4f} (predicted {p:.4f})")


print("Ex 4.6, 10^4 atoms:")
o1, s1 = measure(src, Z, rng)
check((o1 == 1).sum(), N, 1.0, "z+ then z: +")
assert np.allclose(s1, UP)                     # collapse (3.10): the state after the outcome is |up>
o2, s2 = measure(src, X, rng)
check((o2 == 1).sum(), N, 0.5, "z+ then x: +", quoted=4946)
keep = o2 == 1                                 # block the x- beam
assert np.allclose(np.abs(s2[keep] @ np.array([1, 1]) / np.sqrt(2)), 1)      # survivors are in |+x>, by (3.10)
o3, _ = measure(s2[keep], Z, rng)
check((o3 == -1).sum(), keep.sum(), 0.5, "z+ then x+ then z: -", quoted=2423)
pin((o3 == -1).mean(), 0.4899, 5e-5, "fraction quoted in Example 4.6")
pin((o3 == -1).sum() / N, 0.25, 0.02, "a quarter of the original beam ends in |down>")
print(f"  binomial sigmas: {np.sqrt(N*0.25):.1f}, {np.sqrt(keep.sum()*0.25):.1f}")
pin([np.sqrt(N * 0.25), np.sqrt(keep.sum() * 0.25)], [50, 35], 0.5, "sigmas quoted in Example 4.6")
th60 = axis(np.pi / 3)
o4, s4 = measure(src, th60, rng)
check((o4 == 1).sum(), N, np.cos(np.pi / 6)**2, "z+ then 60deg: +", quoted=7505)
o5, _ = measure(s4, Z, rng)
check(((o4 == 1) & (o5 == -1)).sum(), N, np.cos(np.pi / 6)**2 * np.sin(np.pi / 6)**2,
      "z+ -> 60deg+ -> z-  (B5)", quoted=1874)
check((o5 == -1).sum(), N, 0.375, "z+ -> 60deg (either) -> z-  (B5)", quoted=3752)


def recombined(psi, n):
    """An SG-n analyzer whose two beams are recombined coherently with no record of the path (notes 3.8 step 4).
    Each path k contributes the amplitude |n,k><n,k|psi>; indistinguishable paths add (Feynman's rule ii).
    Returns the output states and the two path contributions separately."""
    _, Un = eigvecs(n)
    paths = [(psi @ Un[:, k].conj())[:, None] * Un[:, k][None, :] for k in (0, 1)]
    return paths[0] + paths[1], paths


# No measurement happens inside a recombined analyzer: no random number is drawn there. Only the final SGz
# samples. The SGx case uses its own generator so that the draws from `rng`, and with them every seeded
# number quoted in notes.md and solutions.md further down, stay where they were.
rng_rc = np.random.default_rng(2)
for n_rc, name, gen in [(th60, "60deg", rng), (X, "SGx", rng_rc)]:
    out_rc, (path_p, path_m) = recombined(src, n_rc)
    _, Un = eigvecs(n_rc)
    assert np.allclose(Un @ Un.conj().T, np.eye(2))                    # completeness, LA-8 3.12
    assert np.allclose(out_rc, src)                                    # the state is still |up>
    amp_dn = path_p @ DN.conj() + path_m @ DN.conj()                   # <down|n+><n+|up> + <down|n-><n-|up>
    assert np.allclose(amp_dn, 0)                                      # the two paths cancel: P(z-) = 0
    # with the path recorded, the same two paths give probabilities that add (rule iii)
    p_recorded = np.abs(path_p @ DN.conj())**2 + np.abs(path_m @ DN.conj())**2
    assert np.allclose(p_recorded, 0.5 if name == "SGx" else 0.375)
    o6, _ = measure(out_rc, Z, gen)
    check((o6 == -1).sum(), N, 0.0, f"z+ -> {name} recombined -> z-")
# B5(c): the interference term that the which-path record removes is exactly -0.375 at 60 degrees
_, U60 = eigvecs(th60)
a_p, a_m = (DN.conj() @ U60[:, 0]) * (U60[:, 0].conj() @ UP), (DN.conj() @ U60[:, 1]) * (U60[:, 1].conj() @ UP)
assert np.isclose(2 * np.real(np.conj(a_p) * a_m), -0.375) and np.isclose(abs(a_p + a_m), 0)
pin([abs(a_p)**2, abs(a_m)**2], [0.1875, 0.1875], 1e-12, "B5(b) path probabilities")

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


def ket(theta, phi=0.0):
    """the state (3.15): cos(theta/2)|up> + e^{i phi} sin(theta/2)|down>, spin up along axis(theta, phi)"""
    return np.array([np.cos(theta / 2), np.exp(1j * phi) * np.sin(theta / 2)])


# notes 3.9 / Exercise B4: P(+n' | state at n) = cos^2(gamma/2) = (1 + n.n')/2, and the projector identity
# |n',+><n',+| = (I + n'.sigma)/2. A deterministic grid of axis pairs (no random numbers).
grid = [(t_, p_) for t_ in np.linspace(0, np.pi, 5) for p_ in np.linspace(0, 2 * np.pi, 5)[:-1]]
for th1, ph1 in grid:
    assert np.allclose(np.outer(ket(th1, ph1), ket(th1, ph1).conj()), 0.5 * np.eye(2) + n_dot_S(axis(th1, ph1)))
    for th2, ph2 in grid:
        cosg = axis(th1, ph1) @ axis(th2, ph2)
        P_ = abs(np.vdot(ket(th2, ph2), ket(th1, ph1)))**2
        assert np.isclose(P_, 0.5 * (1 + cosg)) and np.isclose(P_, np.cos(np.arccos(np.clip(cosg, -1, 1)) / 2)**2)
assert np.isclose(abs(np.vdot(ket(2.2, 4.0), ket(np.pi - 2.2, 4.0 + np.pi))), 0)      # antipodes are orthogonal
# the conjugate in the first slot (notes section 8): <+y|+y> is 1 with it, 0 without it
plus_y_ = np.array([1, 1j]) / np.sqrt(2)
assert np.isclose(np.vdot(plus_y_, plus_y_), 1) and np.isclose(plus_y_ @ plus_y_, 0)
assert np.isclose(np.vdot(plus_y_, UP), 1 / np.sqrt(2))

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
# A complex overlap <t1|t2> = r e^{i delta} (Exercise B8(a)). This is the only place where the conjugate in
# (3.16) can be tested: with a real overlap, <t1|t2> and <t2|t1> are the same number.
r_ov, delta = 0.6, 0.9
t2c = np.exp(1j * delta) * np.array([r_ov, np.sqrt(1 - r_ov**2)])
ov_c = np.vdot(t1, t2c)                                                       # <t1|t2>, conjugating t1
assert np.isclose(ov_c, r_ov * np.exp(1j * delta))
P_c = (np.abs(p1[:, None] * t1[None, :] + p2[:, None] * t2c[None, :])**2).sum(1)   # explicit joint state
assert np.allclose(P_c, pattern(xs, ov_c))                                    # (3.16)
assert not np.allclose(P_c, pattern(xs, np.conj(ov_c)))                       # the wrong conjugate fails
two_k = 2 * np.pi * d_ds / (lam_ds * L_ds)                                    # phi_{1,2} = A e^{-+ i k x}
A2 = np.abs(p1)**2
assert np.allclose(P_c, 2 * A2 * (1 + r_ov * np.cos(two_k * xs + delta)))     # B8(a)
assert abs(visibility(xs, P_c) - r_ov) < 1e-3                                 # visibility = |<t1|t2>|
central = np.abs(xs) < np.pi / two_k                                          # one fringe spacing around x = 0
x_shift = xs[central][np.argmax((P_c / (2 * A2))[central])]                   # the central fringe moves to -delta/2k
assert abs(x_shift - (-delta / two_k)) < 0.06, x_shift
assert np.isclose(xs[central][np.argmax((pattern(xs, r_ov) / (2 * A2))[central])], 0)   # real overlap: no shift
print(f"  complex overlap {r_ov} e^(i {delta}): visibility {visibility(xs, P_c):.3f}, "
      f"central fringe at x = {x_shift:.2f} (predicted {-delta/two_k:.2f})")
# unequal slits: the visibility is 2|phi1||phi2||<t1|t2>|/(|phi1|^2 + |phi2|^2), below |<t1|t2>|
P_uneq = np.abs(p1)**2 + np.abs(0.5 * p2)**2 + 2 * np.real(np.conj(p1) * 0.5 * p2 * 0.5)
Pn_uneq = P_uneq / (np.abs(p1)**2 + np.abs(0.5 * p2)**2)
assert abs((Pn_uneq.max() - Pn_uneq.min()) / (Pn_uneq.max() + Pn_uneq.min()) - 2 * 0.5 * 0.5 / 1.25) < 1e-3
# electrons one at a time: sample 3000 hits from each pattern
fig, axs = plt.subplots(3, 1, figsize=(7, 6), sharex=True)
for axx, ov, title in zip(axs, (1.0, 0.0, None), ("amplitudes add (no tag)",
                                                  "orthogonal which-path tags", "classical: P1 + P2")):
    P = pattern(xs, ov) if ov is not None else np.abs(p1)**2 + np.abs(p2)**2
    cdf = np.cumsum(P); cdf /= cdf[-1]
    hits = np.interp(rng.random(3000), cdf, xs)
    # Each electron lands at one point; the pattern is in the statistics. Count the hits within an eighth of
    # a fringe spacing of a dark fringe (a quarter of the screen). Amplitudes adding leave it almost empty:
    # for a flat envelope the fraction is 1/4 - sqrt2/(2 pi) = 2.5%. Orthogonal tags, or classical
    # probabilities, put a quarter of the hits there. (The sinc envelope moves both numbers by under 0.5%.)
    spacing = lam_ds * L_ds / d_ds
    near_dark = lambda x: np.abs((x / spacing) % 1 - 0.5) < 0.125
    frac_dark = near_dark(hits).mean()
    expect = np.sum(P[near_dark(xs)]) / np.sum(P)
    assert abs(frac_dark - expect) < 4 * np.sqrt(expect * (1 - expect) / 3000), (title, frac_dark, expect)
    assert abs(expect - (0.25 - np.sqrt(2) / (2 * np.pi) if ov == 1.0 else 0.25)) < 5e-3
    assert (frac_dark < 0.05) if ov == 1.0 else (frac_dark > 0.2)
    print(f"  3000 electrons, {title}: {frac_dark:.3f} of hits near a dark fringe (pattern says {expect:.3f})")
    axx.hist(hits, bins=120, density=True, color="0.75")
    axx.plot(xs, P / np.trapezoid(P, xs) if hasattr(np, "trapezoid") else P / np.trapz(P, xs), "k")
    axx.set_title(title, fontsize=9)
axs[-1].set_xlabel("screen position x (arbitrary units; fringe spacing lambda L / d = 20)")
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
# notes 3.11: in |up>, Delta S_x = Delta S_y = hbar/2 and the bound (hbar/2)|<S_z>| = hbar^2/4 is attained
assert np.isclose(np.sqrt(var(SX, UP)), 0.5) and np.isclose(np.sqrt(var(SY, UP)), 0.5)
assert np.isclose(np.sqrt(var(SX, UP) * var(SY, UP)), 0.5 * abs(UP.conj() @ SZ @ UP))
# B6(c): |+x> also attains the bound, with both sides zero
plus_x = (UP + DN) / np.sqrt(2)
assert np.isclose(var(SX, plus_x), 0) and np.isclose(plus_x.conj() @ SZ @ plus_x, 0)
# expectation value from Born: sum a_i p_i = <psi|A|psi>
w, U = eigvecs(axis(1.0, 2.0))
psi = (UP + 2j * DN) / np.sqrt(5)
close(np.sum(w * np.abs(U.conj().T @ psi)**2), np.real(psi.conj() @ n_dot_S(axis(1.0, 2.0)) @ psi), 1e-12, "<A>")

# Larmor precession: H = w0 S_z, start in |+x>
w0 = 2.0
for t in np.linspace(0, 5, 11):
    Ut = expm(-1j * w0 * SZ * t)
    assert np.allclose(Ut.conj().T @ Ut, np.eye(2))
    pt = Ut @ plus_x
    assert np.allclose(bloch(pt), [np.cos(w0 * t), np.sin(w0 * t), 0])
    assert np.isclose(abs(plus_x.conj() @ pt)**2, np.cos(w0 * t / 2)**2)
# Theorem 3.3 for a Hamiltonian that is not diagonal: U = exp(-i H t) is unitary, it solves the Schrodinger
# equation (3.11), norms and Born probabilities are conserved, and Hermiticity is what does it.
Hn = 0.7 * SX - 0.4 * SY + 0.3 * SZ + 0.2 * np.eye(2)
assert np.allclose(Hn, Hn.conj().T)
psi0 = (UP + 2j * DN) / np.sqrt(5)
for t in (0.3, 1.7, 9.0):
    Ut = expm(-1j * Hn * t)
    assert np.allclose(Ut.conj().T @ Ut, np.eye(2))
    assert np.isclose(np.linalg.norm(Ut @ psi0), 1)
    assert np.isclose((np.abs(eigvecs(axis(1.0, 2.0))[1].conj().T @ (Ut @ psi0))**2).sum(), 1)
    dpsi = (expm(-1j * Hn * (t + 1e-6)) - expm(-1j * Hn * (t - 1e-6))) @ psi0 / 2e-6
    assert np.allclose(1j * dpsi, Hn @ Ut @ psi0, atol=1e-8)                 # i d|psi>/dt = H|psi>
H_bad = Hn + 0.1j * SZ                                                       # not Hermitian
assert not np.allclose(expm(-1j * H_bad).conj().T @ expm(-1j * H_bad), np.eye(2))
# B7(d): electron in B = 1 T with the g = 2.0023 of the exercise
G_B7 = 2.0023
close(G_B7, -PC["electron g factor"][0], 1e-5, "g")
w_e = G_B7 * MU_B * 1.0 / HBAR
print(f"B7: electron in 1 T: omega0 = {w_e:.4e} rad/s, period = {2*np.pi/w_e*1e12:.2f} ps")
pin(w_e, 1.7608e11, 5e6, "B7 omega0"); pin(2 * np.pi / w_e * 1e12, 35.68, 5e-3, "B7 period, ps")
pin(w_e / (2 * np.pi) / 1e9, 28, 0.5, "B7 frequency, GHz")
pin(G_B7 * 9.274e-24 / 1.0546e-34, 1.7608e11, 5e6, "B7 arithmetic as displayed in solutions.md")

# %% [markdown]
# ## Your turn
# 1. **Compton shift vs angle.** Plot lam'(theta) for Mo K-alpha, 0 to 180 degrees, and the electron's
#    kinetic energy. (Reference below.) Where is the shift half its maximum?
# 2. **Hydrogen-like ions.** Redo Part 2 for He+ (Z = 2): E_n = -Z^2 E_R/n^2, r_n = n^2 a0/Z. Show that the
#    He+ n = 4 -> 3 line sits at 468.65 nm for infinite nuclear mass in vacuum (468.6 nm in air once the
#    reduced mass is included) and that every *even*-upper-level He+ line n = 2k -> 4 coincides
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
pin(He43, 468.65, 5e-3, "He+ 4->3, infinite mass, vacuum")
R_He = Rinf / (1 + ME / PC["alpha particle mass"][0])
pin(line(4, 3, R_He, Z=2) / N_AIR, 468.6, 0.05, "He+ 4->3, reduced mass, air")
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
pin([HC_EVNM / 400, HC_EVNM / 700], [3.100, 1.771], 5e-4, "A1 photon energies")
pin(E_radio, 4.136e-9, 5e-13, "A1 radio photon"); pin(lam_thr, 539.1, 0.05, "A1 threshold")
pin(E_radio / W_NA, 2e-9, 5e-10, "A1: fraction of W per radio photon")
pin(W_NA / E_radio, 5.6e8, 5e6, "A1: photons needed")
print(f"  A2: E2 = {-ER/4:.3f} eV, E3 = {-ER/9:.3f} eV, Lyman-alpha = {line(2, 1):.2f} nm, "
      f"Hbeta = {line(4, 2):.2f} nm")
pin([-ER / 4, -ER / 9], [-3.401, -1.512], 5e-4, "A2 energies"); pin(line(2, 1), 121.50, 5e-3, "A2 Lyman-alpha")
pin(line(4, 2), 486.01, 5e-3, "A2 H-beta"); pin(HC_EVNM / (ER / 4), 364.5, 0.05, "A2 ionization photon")
lam_n = H / np.sqrt(2 * MN * KB * 300)
print(f"  A5: kT(300 K) = {KB*300/EV:.5f} eV, neutron lambda = {lam_n*1e9:.4f} nm")
pin(KB * 300, 4.142e-21, 5e-25, "A5 kT in J"); pin(KB * 300 / EV, 0.02585, 5e-6, "A5 kT in eV")
pin(MN, 1.675e-27, 5e-31, "A5 neutron mass"); pin(np.sqrt(2 * MN * KB * 300), 3.725e-24, 5e-28, "A5 p")
pin(lam_n * 1e9, 0.1779, 5e-5, "A5 lambda")
E100 = 100e3
lam100 = HC_EVNM / E100 * 1e3                      # pm
lam100p = lam100 + 2 * lamC * 1e12
print(f"  B1: 100 keV: lam = {lam100:.3f} pm, back-scattered {lam100p:.3f} pm, "
      f"E' = {HC_EVNM/(lam100p*1e-3)/1e3:.2f} keV, max KE = {E100/1e3 - HC_EVNM/(lam100p*1e-3)/1e3:.2f} keV")
E100p = HC_EVNM / (lam100p * 1e-3) / 1e3                # keV
pin(lam100, 12.398, 5e-4, "B1 lambda"); pin(lam100p, 17.251, 5e-4, "B1 lambda'"); pin(E100p, 71.87, 5e-3, "B1 E'")
pin(E100 / 1e3 - E100p, 28.13, 5e-3, "B1 max KE"); pin(ME * C**2 / EV / 2e3, 255.5, 0.05, "B1 m c^2 / 2")
close(1 / E100p, 1 / 100 + 2 / (ME * C**2 / EV / 1e3), 1e-6, "B1(c): 1/E' = 1/E + 2/(m c^2)")
print(f"  B2: R_H = {RH:.6e}, R_D = {RD:.6e} 1/m; H-alpha {Ha_H:.3f}, D-alpha {Ha_D:.3f} nm, "
      f"shift {Ha_H-Ha_D:.3f} nm")
pin([MP / ME, M_D / ME], [1836.15, 3670.48], 5e-3, "B2 mass ratios"); pin(M_D, 3.3436e-27, 5e-32, "deuteron mass")
pin([RH, RD], [1.096776e7, 1.097074e7], 5, "B2 Rydbergs"); pin([Ha_H, Ha_D], [656.470, 656.291], 5e-4, "B2 lines")
pin(Ha_H - Ha_D, 0.179, 5e-4, "B2 isotope shift"); assert Ha_D < Ha_H           # deuterium on the short side
print(f"  C1: x_peak = {xw:.6f}, b = {b_wien*1e3:.4f} mm K, peak(5000 K) = {b_wien/T5*1e9:.1f} nm, "
      f"Sun 5772 K -> {b_wien/5772*1e9:.1f} nm")
# C2: hidden-variable model, P(+n | z+) = 1 - theta/pi vs quantum cos^2(theta/2)
mvec = rng.normal(size=(200000, 3)); mvec /= np.linalg.norm(mvec, axis=1)[:, None]
mz = mvec[mvec[:, 2] > 0]
assert np.all((mz @ axis(0)) > 0)                           # C2(a): repeating SGz gives + again, always
c2_table = []
for th in (np.pi / 3, np.pi / 2, 2 * np.pi / 3):
    p_hv = (mz @ axis(th) > 0).mean()
    print(f"  C2: theta = {np.degrees(th):.0f} deg: hidden-variable {p_hv:.4f} (1 - th/pi = {1-th/np.pi:.4f}), "
          f"quantum {np.cos(th/2)**2:.4f}")
    assert abs(p_hv - (1 - th / np.pi)) < 0.005                # about 3 binomial sigma with 10^5 atoms
    c2_table.append((1 - th / np.pi, p_hv, np.cos(th / 2)**2))
# the table of solutions C2, as quoted (the Monte Carlo column is seed-dependent: seed 0, the draws above)
c2_table = np.array(c2_table)
pin(c2_table[:, 0], [0.6667, 0.5000, 0.3333], 5e-5, "C2 model column")
pin(c2_table[:, 1], [0.6670, 0.4990, 0.3337], 5e-5, "C2 Monte Carlo column (seeded)")
pin(c2_table[:, 2], [0.7500, 0.5000, 0.2500], 5e-5, "C2 quantum column")
pin(c2_table[0, 2] - c2_table[0, 0], 0.08, 5e-3, "C2: 8 percentage points at 60 degrees")
# B4(c), B9: psi = (|up> + 2i|down>)/sqrt5
psi9 = (UP + 2j * DN) / np.sqrt(5)
n9 = bloch(psi9)
th9 = 2 * np.arccos(1 / np.sqrt(5))
assert np.allclose(n9, [0, 0.8, -0.6]) and np.allclose(n9, axis(th9, np.pi / 2))
pin(np.degrees(th9), 126.87, 5e-3, "B4(c) theta")
assert np.allclose(n9 / 2, [0, 0.4, -0.3])                  # B9(a): <S> in units of hbar
plus_y, minus_y = (UP + 1j * DN) / np.sqrt(2), (UP - 1j * DN) / np.sqrt(2)
py = abs(np.vdot(plus_y, psi9))**2, abs(np.vdot(minus_y, psi9))**2
assert np.allclose(py, [0.9, 0.1]) and np.isclose(0.5 * (py[0] - py[1]), 0.4)
print(f"  B4/B9: Bloch vector {np.round(n9, 4)}, theta = {np.degrees(th9):.2f} deg, phi = 90 deg; "
      f"P(+y) = {py[0]:.2f}, P(-y) = {py[1]:.2f}, <S_y> = {0.5*(py[0]-py[1]):.2f} hbar")
gam = 0.4                                                  # B8(b) reference; any phase gamma, it drops out
for beta in (0, np.pi / 6, np.pi / 3, np.pi / 2):
    t1, t2 = np.array([1, 0]), np.array([np.cos(beta), np.exp(1j * gam) * np.sin(beta)])
    ov = np.vdot(t1, t2)                                   # <t1|t2> = cos(beta), conjugating t1
    assert np.isclose(ov, np.cos(beta))
    p1, p2 = amps(xs)
    joint = p1[:, None] * t1[None, :] + p2[:, None] * t2[None, :]
    P = (np.abs(joint)**2).sum(1)
    assert np.allclose(P, pattern(xs, ov))
    assert abs(visibility(xs, P) - abs(np.cos(beta))) < 1e-3
print(f"  B8: visibilities |cos beta| = {[float(round(abs(np.cos(b_)), 3)) for b_ in (0, np.pi/6, np.pi/3, np.pi/2)]}")
pin(np.abs(np.cos([0, np.pi / 6, np.pi / 3, np.pi / 2])), [1, 0.866, 0.5, 0], 5e-4, "B8 visibilities")
print(f"  A1: photons needed per electron at 1 MHz: {W_NA/E_radio:.2e}; A5: p_n = {np.sqrt(2*MN*KB*300):.4e} kg m/s")
print(f"  B1: E -> inf back-scatter limit m c^2/2 = {ME*C**2/EV/2e3:.1f} keV")
print(f"  B2: m_p/m = {MP/ME:.2f}, m_d/m = {M_D/ME:.2f}")
print(f"  B5: sigmas {np.sqrt(N*0.1875*0.8125):.1f}, {np.sqrt(N*0.375*0.625):.1f}")
pin([np.sqrt(N * 0.1875 * 0.8125), np.sqrt(N * 0.375 * 0.625)], [39, 48], 0.5, "B5(d) sigmas")
print(f"  B7: precession frequency {w_e/(2*np.pi)/1e9:.2f} GHz")
a_rad = 8 * np.pi**5 * KB**4 / (15 * H**3 * C**3)
print(f"  C1: radiation constant a = {a_rad:.4e} J m^-3 K^-4")
pin(a_rad, 7.566e-16, 5e-20, "C1 a"); pin(b_wien / 5772 * 1e9, 502.0, 0.05, "C1 solar peak")
x_it, iterates = 5.0, []
for _ in range(3):                                          # C1(a): fixed-point iteration from x = 5
    x_it = 5 * (1 - np.exp(-x_it)); iterates.append(x_it)
pin(iterates, [4.9663, 4.9652, 4.9651], 5e-5, "C1 iterates")
print(f"  C3: n = 10 ratio {ER*(1/81-1/100)/(2*ER/1000):.5f}, n = 100 ratio {(1/99**2-1/100**2)/(2/100**3):.5f}")

# %%
if __name__ == "__main__":
    dt = time.time() - T0
    print(f"All checks passed in {dt:.1f} s.")
    assert dt < 90
