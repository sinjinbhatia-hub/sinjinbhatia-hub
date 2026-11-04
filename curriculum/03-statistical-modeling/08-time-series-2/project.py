# %% [markdown]
# # ST-8 Project: Spectra, the Kalman filter, and GARCH
#
# Builds the three views of a time series from numpy/scipy and checks every
# number in the notes. Parts: (1) ARMA spectra three ways, (2) linear filters,
# (3) the periodogram's inconsistency and the smoothed fix, (4) tapering and a
# hidden periodicity, (5) transfer-function identification of a damped
# oscillator, (6) a generic Kalman filter and RTS smoother, (7) the
# constant-velocity tracker and the local level steady state, (8) ML for the
# local level and local linear trend, (9) ARMA exact likelihood two ways,
# (10) a structural decomposition, (11) an EKF pendulum, (12) GARCH(1,1),
# (13) exercise numbers.  Run: `MPLBACKEND=Agg python3 project.py`.

# %%
import math
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import signal, stats, optimize, linalg

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(8)
T0 = time.time()


def arma_spectrum(om, phi=(), theta=(), s2=1.0):
    """f(w) = s2 |theta(e^{-2 pi i w})|^2 / |phi(e^{-2 pi i w})|^2, eq. (3.9)."""
    z = np.exp(-2j * np.pi * np.asarray(om, float))
    num = 1 + sum(t * z ** (k + 1) for k, t in enumerate(theta))
    den = 1 - sum(p * z ** (k + 1) for k, p in enumerate(phi))
    return s2 * np.abs(num) ** 2 / np.abs(den) ** 2


def periodogram(x):
    """I(w_j) = |sum_t x_t e^{-2 pi i w_j t}|^2 / n at w_j = j/n, j = 0..n//2 (eq. 3.13)."""
    x = np.asarray(x, float)
    n = x.shape[-1]
    X = np.fft.rfft(x - x.mean(axis=-1, keepdims=True), axis=-1)
    return np.arange(n // 2 + 1) / n, np.abs(X) ** 2 / n


def daniell(I, m):
    """Daniell smoother: average of 2m+1 neighbouring ordinates (eq. 3.15), reflected at the ends."""
    Ip = np.concatenate([I[..., 1:m + 1][..., ::-1], I, I[..., -m - 1:-1][..., ::-1]], axis=-1)
    k = np.ones(2 * m + 1) / (2 * m + 1)
    return np.apply_along_axis(lambda r: np.convolve(r, k, mode="valid"), -1, Ip)


def simulate_arma(phi, theta, n, reps=None, burn=500):
    """Simulate ARMA(p,q) with unit-variance Gaussian noise by recursive filtering."""
    shape = (n + burn,) if reps is None else (reps, n + burn)
    w = rng.standard_normal(shape)
    x = signal.lfilter(np.r_[1, theta], np.r_[1, -np.asarray(phi, float)], w, axis=-1)
    return x[..., burn:]


def acf_theory_arma11(phi, theta, hmax, s2=1.0):
    """ARMA(1,1) autocovariance gamma(0..hmax) (ST-7)."""
    g = np.empty(hmax + 1)
    g[0] = s2 * (1 + 2 * phi * theta + theta ** 2) / (1 - phi ** 2)
    g[1] = s2 * (1 + phi * theta) * (phi + theta) / (1 - phi ** 2)
    for h in range(2, hmax + 1):
        g[h] = phi * g[h - 1]
    return g


# %% [markdown]
# ## Part 1: ARMA spectra three ways (notes §3.2-3.3, Examples 4.1-4.2)
# Closed form (3.9) vs the Fourier series of the theoretical ACF (3.4) vs the
# averaged periodogram of 200 long simulations.

# %%
om = np.linspace(0, 0.5, 201)
H = 400
hs = np.arange(-H, H + 1)
gam = {"AR(1) phi=.7": 0.7 ** np.abs(hs) / (1 - 0.49),
       "MA(1) theta=.5": np.where(hs == 0, 1.25, np.where(np.abs(hs) == 1, 0.5, 0.0))}
# AR(2) ACF from the companion matrix: gamma(h) = [T^h Gamma0]_{00}
Tc = np.array([[1.0, -0.5], [1.0, 0.0]])
G0 = linalg.solve_discrete_lyapunov(Tc, np.array([[1.0, 0], [0, 0]]))
g2 = [G0[0, 0]]
M = G0.copy()
for h in range(1, H + 1):
    M = Tc @ M
    g2.append(M[0, 0])
g2 = np.array(g2)
gam["AR(2) 1,-.5"] = g2[np.abs(hs)]
closed = {"AR(1) phi=.7": arma_spectrum(om, [0.7]),
          "MA(1) theta=.5": arma_spectrum(om, (), [0.5]),
          "AR(2) 1,-.5": arma_spectrum(om, [1.0, -0.5])}
for k in closed:
    fs = np.real(np.exp(-2j * np.pi * np.outer(om, hs)) @ gam[k])
    assert np.allclose(fs, closed[k], rtol=1e-8), k
print("Spectral density = Fourier series of the ACF: closed form matches for AR(1), MA(1), AR(2).")

fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
for ax, (k, phi, th) in zip(axes, [("AR(1) phi=.7", [0.7], ()), ("MA(1) theta=.5", (), [0.5]),
                                     ("AR(2) 1,-.5", [1.0, -0.5], ())]):
    w_j, I = periodogram(simulate_arma(phi, th, 1024, reps=200))
    Ibar = I.mean(axis=0)
    sel = (w_j > 0.02) & (w_j < 0.48)
    assert np.allclose(Ibar[sel], arma_spectrum(w_j[sel], phi, th), rtol=0.3), k
    ax.plot(w_j, Ibar, lw=0.6, alpha=0.7, label="mean periodogram (200 runs)")
    ax.plot(om, closed[k], "k", label="closed form")
    ax.set_title(k); ax.set_xlabel("w (cycles/sample)"); ax.set_yscale("log")
axes[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(FIG / "spectra.png", dpi=120); plt.close(fig)

# Example 4.1
for phi in (0.7, -0.7):
    vals = arma_spectrum([0, 0.25, 0.5], [phi])
    print(f"Ex 4.1 AR(1) phi={phi:+.1f}: f(0)={vals[0]:.4f}, f(1/4)={vals[1]:.4f}, f(1/2)={vals[2]:.4f}")
assert np.allclose(arma_spectrum([0, .25, .5], [0.7]), [1 / 0.09, 1 / 1.49, 1 / 2.89])
assert np.allclose(arma_spectrum([0, .25, .5], [-0.7]), [1 / 2.89, 1 / 1.49, 1 / 0.09])
print(f"Ex 4.1 gamma(0)={1/0.51:.3f}, f(0)/f(1/2)={(1/0.09)/(1/2.89):.1f}")

# Example 4.2: AR(2) peak, cos(2 pi w*) = phi1 (phi2 - 1) / (4 phi2), eq. (3.11)
p1, p2 = 1.0, -0.5
c_star = p1 * (p2 - 1) / (4 * p2)
w_star = np.arccos(c_star) / (2 * np.pi)
fine = np.linspace(0, 0.5, 200001)
assert abs(fine[np.argmax(arma_spectrum(fine, [p1, p2]))] - w_star) < 1e-5
roots = np.roots([-p2, -p1, 1])          # roots of 1 - p1 z - p2 z^2
quasi = 2 * np.pi / abs(np.angle(1 / roots[0]))
print(f"Ex 4.2 AR(2): cos(2pi w*)={c_star:.3f}, w*={w_star:.4f}, period={1/w_star:.2f}, "
      f"f(w*)={arma_spectrum(w_star, [p1, p2]):.3f}, f(0)={arma_spectrum(0, [p1, p2]):.3f}; "
      f"root modulus {abs(1/roots[0]):.4f}, quasi-period {quasi:.2f}")
assert abs(c_star - 0.75) < 1e-12 and abs(quasi - 8) < 1e-9

# %% [markdown]
# ## Part 2: linear filters (§3.4, Example 4.3)
# Output spectrum of filtered white noise = |psi(e^{-2 pi i w})|^2 (eq. 3.8).

# %%
filters = {"3-pt MA": np.ones(3) / 3, "difference": np.array([1.0, -1.0]),
           "5-pt MA": np.ones(5) / 5}
wtest = np.array([0, 1 / 6, 0.25, 1 / 3, 0.5])
for name, b in filters.items():
    _, h = signal.freqz(b, 1, worN=2 * np.pi * wtest)
    print(f"Ex 4.3 {name:10s} |psi|^2 at w={np.round(wtest, 3)}: {np.round(np.abs(h) ** 2, 4)}")
_, hma = signal.freqz(np.ones(3) / 3, 1, worN=2 * np.pi * wtest)
assert np.allclose(np.abs(hma) ** 2, ((1 + 2 * np.cos(2 * np.pi * wtest)) / 3) ** 2)
_, hd = signal.freqz([1, -1], 1, worN=2 * np.pi * wtest)
assert np.allclose(np.abs(hd) ** 2, [0, 1, 2, 3, 4])
wn = rng.standard_normal((300, 2048))
fig, ax = plt.subplots(figsize=(6, 3.6))
for name, b in filters.items():
    w_j, I = periodogram(signal.lfilter(b, 1, wn, axis=-1)[:, 10:])
    _, h = signal.freqz(b, 1, worN=2 * np.pi * w_j)
    Ibar = I.mean(axis=0)
    sel = (np.abs(h) ** 2 > 0.05) & (w_j > 0)
    assert np.allclose(Ibar[sel], np.abs(h[sel]) ** 2, rtol=0.3), name
    ax.plot(w_j, Ibar, lw=0.6, label=f"{name} (sim)")
    ax.plot(w_j, np.abs(h) ** 2, "k", lw=0.8)
ax.set_xlabel("w"); ax.set_ylabel("output spectrum"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(FIG / "filters.png", dpi=120); plt.close(fig)
print("Filtered white noise: averaged periodograms match |psi|^2 for all three filters.")
# Slutsky: six 2-point sums then two differences -> peaked spectrum, cos(2 pi w*) = (s-d)/(s+d)
b_sl = np.array([1.0])
for _ in range(6):
    b_sl = np.convolve(b_sl, [1, 1])
for _ in range(2):
    b_sl = np.convolve(b_sl, [1, -1])
_, hs_ = signal.freqz(b_sl, 1, worN=2 * np.pi * fine[::100])
w_sl = fine[::100][np.argmax(np.abs(hs_))]
print(f"Slutsky filter (1+B)^6(1-B)^2: spectral peak at w={w_sl:.4f} (period {1/w_sl:.2f}); "
      f"theory cos(2 pi w*) = (6-2)/(6+2) = 0.5 -> w* = 1/6")
assert abs(w_sl - np.arccos(0.5) / (2 * np.pi)) < 1e-3

# %% [markdown]
# ## Part 3: the periodogram is not consistent; smoothing fixes it (§3.5, Example 4.4)

# %%
for n in (128, 1024, 8192):
    _, I = periodogram(rng.standard_normal((2000, n)))
    j = n // 8                                     # w = 1/8, f = 1
    print(f"white noise n={n:5d}: mean I(1/8)={I[:, j].mean():.3f}, var I(1/8)={I[:, j].var():.3f}")
    assert abs(I[:, j].mean() - 1) < 0.1 and abs(I[:, j].var() - 1) < 0.2
q_lo, q_hi = stats.chi2.ppf([0.025, 0.975], 10)
print(f"Ex 4.4: chi2_10 quantiles {q_lo:.3f}, {q_hi:.3f}; 95% CI = [{10/q_hi:.3f} fhat, {10/q_lo:.3f} fhat]; "
      f"raw periodogram CI from chi2_2: [{2/stats.chi2.ppf(.975,2):.3f} I, {2/stats.chi2.ppf(.025,2):.1f} I]")
assert abs(10 / q_hi - 0.488) < 1e-3 and abs(10 / q_lo - 3.080) < 1e-3
print(f"Ex 4.4: fhat=2.5 -> [{2.5*10/q_hi:.2f}, {2.5*10/q_lo:.2f}]; width ratio {q_hi/q_lo:.2f} vs raw "
      f"{stats.chi2.ppf(.975,2)/stats.chi2.ppf(.025,2):.0f}")

# coverage of the nu = 2(2m+1) interval, AR(1) phi=0.5, w = 1/4
X = simulate_arma([0.5], (), 1024, reps=3000)
w_j, I = periodogram(X)
fhat = daniell(I, 2)
j = 256
f_true = arma_spectrum(0.25, [0.5])
cover = np.mean((10 * fhat[:, j] / q_hi <= f_true) & (f_true <= 10 * fhat[:, j] / q_lo))
print(f"Coverage of the chi2_10 interval at w=1/4 over 3000 AR(1) series: {cover:.3f}")
assert 0.93 < cover < 0.97

# bias-variance in the bandwidth at a sharp AR(2) peak (root modulus 0.95, w = 0.1)
phs = [2 * 0.95 * np.cos(2 * np.pi * 0.1), -0.95 ** 2]
X = simulate_arma(phs, (), 1024, reps=400)
w_j, I = periodogram(X)
ms = [0, 1, 2, 4, 8, 16, 32]
ftrue = arma_spectrum(w_j, phs)
jpk = np.argmax(ftrue)
b2, var = [], []
for m in ms:
    fh = daniell(I, m) if m else I
    b2.append((fh[:, jpk].mean() - ftrue[jpk]) ** 2)
    var.append(fh[:, jpk].var())
b2, var = np.array(b2), np.array(var)
best = ms[np.argmin(b2 + var)]
print(f"Bias-variance at a sharp AR(2) peak: MSE-optimal Daniell m={best}; var(m=0)/var(m=8)={var[0]/var[4]:.1f}; "
      f"bias^2(m=32)/bias^2(m=2)={b2[-1]/b2[2]:.0f}")
assert var[0] > 10 * var[4] and b2[-1] > 10 * b2[2] and 0 < best < 32
fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
ax[0].loglog(np.array(ms[1:]) * 2 + 1, b2[1:], "o-", label="bias$^2$")
ax[0].loglog(np.array(ms[1:]) * 2 + 1, var[1:], "s-", label="variance")
ax[0].loglog(np.array(ms[1:]) * 2 + 1, b2[1:] + var[1:], "k--", label="MSE")
ax[0].set_xlabel("span 2m+1"); ax[0].set_title("at a sharp AR(2) peak"); ax[0].legend()
ax[1].semilogy(w_j, I[0], lw=0.4, color="0.6", label="raw I")
ax[1].semilogy(w_j, daniell(I[:1], 8)[0], label="Daniell m=8")
ax[1].semilogy(w_j, ftrue, "k", label="true f")
ax[1].legend(fontsize=8); ax[1].set_xlabel("w")
fig.tight_layout(); fig.savefig(FIG / "bias_variance.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 4: tapering and a hidden periodicity (§3.6)

# %%
n = 512
t = np.arange(n)
x = 10 * np.cos(2 * np.pi * 0.1037 * t) + 0.01 * np.cos(2 * np.pi * 0.3 * t) + 0.001 * rng.standard_normal(n)
hann = signal.windows.hann(n, sym=False)
I_raw = np.abs(np.fft.rfft(x)) ** 2 / n
I_han = np.abs(np.fft.rfft(x * hann)) ** 2 / np.sum(hann ** 2)
w_j = np.arange(n // 2 + 1) / n
far = np.abs(w_j - 0.3) < 0.004
ratio_raw = I_raw[far].max() / I_raw[np.abs(w_j - 0.25) < 0.01].max()
ratio_han = I_han[far].max() / I_han[np.abs(w_j - 0.25) < 0.01].max()
print(f"Leakage: weak line at 0.3 over background at 0.25: raw {ratio_raw:.1f}x, Hann {ratio_han:.0f}x")
assert ratio_raw < 3 and ratio_han > 50
fig, ax = plt.subplots(figsize=(6, 3.4))
ax.semilogy(w_j, I_raw, lw=0.7, label="no taper"); ax.semilogy(w_j, I_han, lw=0.7, label="Hann")
ax.set_xlabel("w"); ax.legend(); fig.tight_layout(); fig.savefig(FIG / "taper.png", dpi=120); plt.close(fig)


def fisher_g(x):
    """Fisher's g = max_j I(w_j) / sum_j I(w_j), j=1..(n-1)//2, with exact p-value (§3.6)."""
    _, I = periodogram(x)
    m = (len(x) - 1) // 2
    I = I[1:m + 1]
    g = I.max() / I.sum()
    k = np.arange(1, int(1 / g) + 1)
    p = np.sum((-1) ** (k - 1) * np.array([float(math.comb(m, int(kk))) for kk in k])
               * (1 - k * g) ** (m - 1))
    return g, min(max(p, 0.0), 1.0)


x = 0.4 * np.cos(2 * np.pi * 37 / 501 * np.arange(501)) + rng.standard_normal(501)
g, p = fisher_g(x)
gn, pn = fisher_g(rng.standard_normal(501))
print(f"Hidden periodicity (amp 0.4 in unit noise, n=501): Fisher g={g:.4f}, p={p:.2e}; pure noise p={pn:.2f}")
assert p < 1e-3
pvals = np.array([fisher_g(rng.standard_normal(201))[1] for _ in range(1000)])
assert abs(np.mean(pvals < 0.05) - 0.05) < 0.02

# %% [markdown]
# ## Part 5: transfer function of a damped oscillator from the cross-spectrum (§3.7)
# x'' + 2 zeta w0 x' + w0^2 x = u(t), discretized (bilinear), driven by white noise.

# %%
w0, zeta, dt = 2 * np.pi * 0.08, 0.05, 1.0
bd, ad, _ = signal.cont2discrete(([1.0], [1, 2 * zeta * w0, w0 ** 2]), dt, method="bilinear")
bd = bd.ravel()
n = 2 ** 15
u = rng.standard_normal(n)
y = signal.lfilter(bd, ad, u) + 0.5 * rng.standard_normal(n)
U, Y = np.fft.rfft(u), np.fft.rfft(y)
m = 16
fuu = daniell(np.abs(U) ** 2 / n, m)
fyy = daniell(np.abs(Y) ** 2 / n, m)
fyu = daniell((Y * np.conj(U)).real / n, m) + 1j * daniell((Y * np.conj(U)).imag / n, m)
Hhat = fyu / fuu
coh = np.abs(fyu) ** 2 / (fuu * fyy)
w_j = np.arange(n // 2 + 1) / n
_, Htrue = signal.freqz(bd, ad, worN=2 * np.pi * w_j)
band = (w_j > 0.02) & (w_j < 0.2)
relerr = np.abs(np.abs(Hhat[band]) - np.abs(Htrue[band])) / np.abs(Htrue[band])
pk = np.argmax(np.abs(Hhat))
print(f"Transfer function: median rel. gain error {np.median(relerr):.3f} on 0.02<w<0.2; "
      f"estimated resonance w={w_j[pk]:.4f} (true {w_j[np.argmax(np.abs(Htrue))]:.4f}); "
      f"coherence at resonance {coh[pk]:.3f}, at w=0.4 {coh[np.argmin(abs(w_j-0.4))]:.3f}")
assert np.median(relerr) < 0.1 and abs(w_j[pk] - 0.08) < 0.003
fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
ax[0].semilogy(w_j, np.abs(Hhat), lw=0.7, label="|f_yu / f_uu|")
ax[0].semilogy(w_j, np.abs(Htrue), "k", lw=0.8, label="true |H|"); ax[0].legend()
ax[1].plot(w_j, coh); ax[1].set_ylabel("squared coherence")
for a in ax: a.set_xlabel("w")
fig.tight_layout(); fig.savefig(FIG / "transfer.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 6: the Kalman filter and RTS smoother, generic (§3.9-3.11)
# alpha_{t+1} = T alpha_t + R eta_t, y_t = Z alpha_t + eps_t; (a0, P0) are the
# mean and covariance of alpha_1 before y_1 is seen. NaN in y = missing.


# %%
def kalman_filter(y, T, Z, R, Q, H, a0, P0):
    """Kalman filter, eqs. (3.20)-(3.23). Returns dict with filtered a_{t|t}, P_{t|t},
    predicted a_{t|t-1}, P_{t|t-1}, innovations v_t, F_t, gains K_t, and the
    prediction-error-decomposition log-likelihood (3.27)."""
    y = np.atleast_2d(np.asarray(y, float).T).T if np.ndim(y) == 1 else np.asarray(y, float)
    n, p = y.shape
    k = len(a0)
    RQR = R @ Q @ R.T
    out = {key: np.zeros((n, k)) for key in ("af", "ap")}
    out.update(Pf=np.zeros((n, k, k)), Pp=np.zeros((n, k, k)), v=np.full((n, p), np.nan),
               F=np.zeros((n, p, p)), K=np.zeros((n, k, p)))
    a, P, ll = np.asarray(a0, float), np.asarray(P0, float), 0.0
    for t in range(n):
        out["ap"][t], out["Pp"][t] = a, P
        F = Z @ P @ Z.T + H
        out["F"][t] = F
        if not np.any(np.isnan(y[t])):
            v = y[t] - Z @ a
            K = np.linalg.solve(F, Z @ P).T          # P Z' F^{-1}
            a = a + K @ v
            P = P - K @ F @ K.T
            P = 0.5 * (P + P.T)
            out["v"][t], out["K"][t] = v, K
            ll -= 0.5 * (p * np.log(2 * np.pi) + np.linalg.slogdet(F)[1] + v @ np.linalg.solve(F, v))
        out["af"][t], out["Pf"][t] = a, P
        a, P = T @ a, T @ P @ T.T + RQR
    out["ll"] = ll
    return out


def rts_smoother(kf, T):
    """Rauch-Tung-Striebel smoother, eq. (3.26): a_{t|n}, P_{t|n}."""
    af, Pf, ap, Pp = kf["af"], kf["Pf"], kf["ap"], kf["Pp"]
    n = len(af)
    as_, Ps = af.copy(), Pf.copy()
    for t in range(n - 2, -1, -1):
        J = np.linalg.solve(Pp[t + 1].T, (Pf[t] @ T.T).T).T
        as_[t] = af[t] + J @ (as_[t + 1] - ap[t + 1])
        Ps[t] = Pf[t] + J @ (Ps[t + 1] - Pp[t + 1]) @ J.T
    return as_, Ps


def ljung_box(e, m):
    """Ljung-Box Q = n(n+2) sum_{h<=m} rho_h^2/(n-h) ~ chi2_m (ST-7)."""
    e = e - e.mean()
    n = len(e)
    r = np.array([e[h:] @ e[:-h] for h in range(1, m + 1)]) / (e @ e)
    Qs = n * (n + 2) * np.sum(r ** 2 / (n - np.arange(1, m + 1)))
    return Qs, stats.chi2.sf(Qs, m)


# %% [markdown]
# ## Part 7: constant-velocity tracker (Example 4.6) and local-level steady state (Example 4.7)

# %%
Tcv = np.array([[1.0, 1.0], [0.0, 1.0]])
Rcv = np.array([[0.5], [1.0]])
Zcv = np.array([[1.0, 0.0]])
Qcv, Hcv = np.array([[1.0]]), np.array([[1.0]])
a00, P00 = np.array([0.0, 1.0]), np.eye(2)
a10, P10 = Tcv @ a00, Tcv @ P00 @ Tcv.T + Rcv @ Qcv @ Rcv.T
kf1 = kalman_filter(np.array([2.0]), Tcv, Zcv, Rcv, Qcv, Hcv, a10, P10)
print("Ex 4.6 predict: a_{1|0} =", a10, " P_{1|0} =", P10.tolist())
print(f"Ex 4.6 update: F={kf1['F'][0,0,0]:.4f}, K={np.round(kf1['K'][0,:,0],4)}, v={kf1['v'][0,0]:.1f}, "
      f"a_1|1={np.round(kf1['af'][0],4)}, P_1|1={np.round(kf1['Pf'][0],4).tolist()}")
assert np.allclose(P10, [[2.25, 1.5], [1.5, 2.0]]) and abs(kf1["F"][0, 0, 0] - 3.25) < 1e-12
assert np.allclose(kf1["af"][0], [1 + 2.25 / 3.25, 1 + 1.5 / 3.25])

# Example 4.5: AR(2) companion form, eigenvalues and innovation variances with exact observation
T2 = np.array([[1.0, -0.5], [1.0, 0.0]])
ev = np.linalg.eigvals(T2)
P2 = linalg.solve_discrete_lyapunov(T2, np.array([[1.0, 0], [0, 0]]))
kf2s = kalman_filter(np.zeros(6), T2, np.array([[1.0, 0]]), np.array([[1.0], [0]]),
                     np.eye(1), np.zeros((1, 1)), np.zeros(2), P2)
print(f"Ex 4.5: eigenvalues {np.round(ev, 4)}, modulus {abs(ev[0]):.4f}; F_t = {np.round(kf2s['F'][:, 0, 0], 4)}")
assert np.allclose(np.sort_complex(ev), [0.5 - 0.5j, 0.5 + 0.5j]) and np.allclose(kf2s['F'][2:, 0, 0], 1.0)

# simulate a trajectory
n = 300
q_a, sd_y = 0.05, 3.0
Qs_, Hs_ = np.array([[q_a]]), np.array([[sd_y ** 2]])
xs = np.zeros((n, 2)); xs[0] = [0, 1]
for t in range(1, n):
    xs[t] = Tcv @ xs[t - 1] + Rcv[:, 0] * np.sqrt(q_a) * rng.standard_normal()
ys = xs[:, 0] + sd_y * rng.standard_normal(n)
kf = kalman_filter(ys, Tcv, Zcv, Rcv, Qs_, Hs_, np.array([0.0, 0.0]), np.diag([100.0, 10.0]))
asm, Psm = rts_smoother(kf, Tcv)
rmse = lambda e: np.sqrt(np.mean(e ** 2))
r_raw, r_f, r_s = rmse(ys - xs[:, 0]), rmse(kf["af"][20:, 0] - xs[20:, 0]), rmse(asm[20:, 0] - xs[20:, 0])
print(f"Tracker position RMSE: raw {r_raw:.2f}, filtered {r_f:.2f}, smoothed {r_s:.2f}")
assert r_s < r_f < r_raw
zstd = kf["v"][:, 0] / np.sqrt(kf["F"][:, 0, 0])
Qlb, plb = ljung_box(zstd[10:], 10)
print(f"Innovations: mean {zstd.mean():.3f}, var {zstd.var():.3f}, Ljung-Box(10) p={plb:.2f}")
assert plb > 0.01 and abs(zstd.var() - 1) < 0.25
z95 = np.mean(np.abs(xs[20:, 0] - asm[20:, 0]) < 1.96 * np.sqrt(Psm[20:, 0, 0]))
assert 0.88 < z95 < 0.99
# missing data: drop t = 150..199
ym = ys.copy(); ym[150:200] = np.nan
kfm = kalman_filter(ym, Tcv, Zcv, Rcv, Qs_, Hs_, np.array([0.0, 0.0]), np.diag([100.0, 10.0]))
asm_m, Psm_m = rts_smoother(kfm, Tcv)
print(f"Missing 150-199: filtered sd at t=199 {np.sqrt(kfm['Pf'][199,0,0]):.2f} vs {np.sqrt(kf['Pf'][199,0,0]):.2f}"
      f" with data; smoothed sd at t=175 {np.sqrt(Psm_m[175,0,0]):.2f}")
assert kfm["Pf"][199, 0, 0] > 5 * kf["Pf"][199, 0, 0] and Psm_m[175, 0, 0] < kfm["Pf"][175, 0, 0]
fig, ax = plt.subplots(1, 2, figsize=(12, 3.6))
ax[0].plot(ys, ".", ms=2, color="0.6", label="y"); ax[0].plot(xs[:, 0], "k", label="true")
ax[0].plot(kf["af"][:, 0], label="filtered"); ax[0].plot(asm[:, 0], label="smoothed"); ax[0].legend(fontsize=7)
ax[1].plot(xs[:, 0] - asm_m[:, 0], label="smoothed error (gap 150-199)")
ax[1].fill_between(np.arange(n), -1.96 * np.sqrt(Psm_m[:, 0, 0]), 1.96 * np.sqrt(Psm_m[:, 0, 0]), alpha=0.3)
ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(FIG / "tracker.png", dpi=120); plt.close(fig)

# local level, q = 1: Riccati P^2 - qP - q = 0 (eq. 3.24)
one = np.array([[1.0]])
kfl = kalman_filter(rng.standard_normal(60), one, one, one, one, one, np.zeros(1), np.array([[1e6]]))
P_ss = (1 + np.sqrt(5)) / 2
K_ss = P_ss / (P_ss + 1)
print(f"Ex 4.7 local level q=1: P_ss={P_ss:.4f}, K_ss={K_ss:.4f}, P_t|t={P_ss/(P_ss+1):.4f}; "
      f"filter gain at t=1,2,3,4,10: {np.round(kfl['K'][[0,1,2,3,9],0,0],4)}")
assert abs(kfl['K'][3, 0, 0] - 13 / 21) < 1e-5
assert abs(kfl["K"][30, 0, 0] - K_ss) < 1e-10 and abs(kfl["Pp"][30, 0, 0] - P_ss) < 1e-10
# steady-state local level = simple exponential smoothing with alpha = K (Holt-Winters link)
yy = np.cumsum(rng.standard_normal(200)) + rng.standard_normal(200)
kfy = kalman_filter(yy, one, one, one, one, one, np.zeros(1), np.array([[1e6]]))
ses = kfy["af"][29, 0]
for t in range(30, 200):
    ses = ses + K_ss * (yy[t] - ses)
assert abs(ses - kfy["af"][-1, 0]) < 1e-8
print("Steady-state local level filter = simple exponential smoothing with alpha = K_ss.")
fig, ax = plt.subplots(figsize=(5, 3.2))
ax.plot(kfl["K"][:15, 0, 0], "o-"); ax.axhline(K_ss, color="k", ls="--")
ax.set_xlabel("t"); ax.set_ylabel("gain K_t"); fig.tight_layout(); fig.savefig(FIG / "gain.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 8: maximum likelihood for local level and local linear trend (§3.12)


# %%
def fit_ml(negll, x0):
    """Minimize a negative log-likelihood in log-variance coordinates; return estimates and SEs
    (observed information via a finite-difference Hessian, delta method back to variances)."""
    res = optimize.minimize(negll, x0, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-6, maxiter=2000))
    th, d = res.x, 1e-3
    k = len(th)
    Hs = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            ei, ej = np.eye(k)[i] * d, np.eye(k)[j] * d
            Hs[i, j] = (negll(th + ei + ej) - negll(th + ei - ej) - negll(th - ei + ej) + negll(th - ei - ej)) / (4 * d * d)
    se_log = np.sqrt(np.diag(np.linalg.inv(Hs)))
    return np.exp(th), np.exp(th) * se_log


n = 800
s_eta, s_eps = 0.5, 1.0
mu = np.cumsum(np.sqrt(s_eta) * rng.standard_normal(n))
yl = mu + np.sqrt(s_eps) * rng.standard_normal(n)
nll_ll = lambda th: -kalman_filter(yl, one, one, one, np.array([[np.exp(th[0])]]), np.array([[np.exp(th[1])]]),
                                   yl[:1], np.array([[1e7]]))["ll"]
est, se = fit_ml(nll_ll, np.zeros(2))
print(f"Local level ML: sigma_eta^2={est[0]:.3f}+-{se[0]:.3f} (true {s_eta}), "
      f"sigma_eps^2={est[1]:.3f}+-{se[1]:.3f} (true {s_eps})")
assert np.all(np.abs(est - [s_eta, s_eps]) < 3 * se)

n = 500
tv = (0.2, 0.01, 1.0)                                        # level, slope, obs variances
Tlt, Zlt = np.array([[1.0, 1.0], [0.0, 1.0]]), np.array([[1.0, 0.0]])
lv, sl = np.zeros(n), np.zeros(n)
for t in range(1, n):
    sl[t] = sl[t - 1] + np.sqrt(tv[1]) * rng.standard_normal()
    lv[t] = lv[t - 1] + sl[t - 1] + np.sqrt(tv[0]) * rng.standard_normal()
ylt = lv + np.sqrt(tv[2]) * rng.standard_normal(n)
nll_lt = lambda th: -kalman_filter(ylt, Tlt, Zlt, np.eye(2), np.diag(np.exp(th[:2])), np.array([[np.exp(th[2])]]),
                                   np.array([ylt[0], 0.0]), np.diag([1e7, 1e7]))["ll"]
est3, se3 = fit_ml(nll_lt, np.log([0.5, 0.05, 0.5]))
print(f"Local linear trend ML: {np.round(est3,4)} +- {np.round(se3,4)} (true {tv})")
assert np.all(np.abs(est3 - tv) < 3 * se3)

# %% [markdown]
# ## Part 9: ARMA(1,1) exact likelihood: state space vs direct MVN (§3.12)

# %%
phi, th = 0.6, 0.4
n = 200
xa = simulate_arma([phi], [th], n)
Ta, Ra, Za = np.array([[phi, 1.0], [0.0, 0.0]]), np.array([[1.0], [th]]), np.array([[1.0, 0.0]])
P0a = linalg.solve_discrete_lyapunov(Ta, Ra @ Ra.T)
ll_kf = kalman_filter(xa, Ta, Za, Ra, np.eye(1), np.zeros((1, 1)), np.zeros(2), P0a)["ll"]
Sig = linalg.toeplitz(acf_theory_arma11(phi, th, n - 1))
ll_mvn = stats.multivariate_normal(np.zeros(n), Sig).logpdf(xa)
print(f"ARMA(1,1) exact log-lik: Kalman {ll_kf:.6f}, MVN {ll_mvn:.6f}")
assert abs(ll_kf - ll_mvn) < 1e-6
assert np.isclose(P0a[0, 0], acf_theory_arma11(phi, th, 1)[0])

# %% [markdown]
# ## Part 10: structural decomposition: trend + quarterly seasonal + irregular (§3.13)

# %%
s = 4
n = 200
Tst = np.zeros((5, 5)); Tst[:2, :2] = [[1, 1], [0, 1]]; Tst[2, 2:5] = -1; Tst[3, 2] = 1; Tst[4, 3] = 1
Zst = np.array([[1.0, 0, 1.0, 0, 0]])
Rst = np.zeros((5, 3)); Rst[0, 0] = Rst[1, 1] = Rst[2, 2] = 1
Qst = np.diag([0.05, 0.001, 0.01])
st = np.zeros((n, 5)); st[0] = [10, 0.1, 3, -1, -1]
for t in range(1, n):
    st[t] = Tst @ st[t - 1] + Rst @ (np.sqrt(np.diag(Qst)) * rng.standard_normal(3))
yst = Zst[0] @ st.T + 0.7 * rng.standard_normal(n)
kfs = kalman_filter(yst, Tst, Zst, Rst, Qst, np.array([[0.49]]), np.zeros(5), np.eye(5) * 1e6)
ass, _ = rts_smoother(kfs, Tst)
r_tr, r_se = rmse(ass[10:, 0] - st[10:, 0]), rmse(ass[10:, 2] - st[10:, 2])
print(f"Structural model: smoothed trend RMSE {r_tr:.3f}, seasonal RMSE {r_se:.3f} (obs noise sd 0.7)")
assert r_tr < 0.5 and r_se < 0.4
fig, ax = plt.subplots(2, 1, figsize=(8, 4.5), sharex=True)
ax[0].plot(yst, color="0.6", lw=0.7); ax[0].plot(st[:, 0], "k"); ax[0].plot(ass[:, 0], "C1--")
ax[1].plot(st[:, 2], "k"); ax[1].plot(ass[:, 2], "C1--")
fig.tight_layout(); fig.savefig(FIG / "structural.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 11: extended Kalman filter for a large-amplitude pendulum (§3.14)

# %%
kp, dt = 9.81, 0.05
fpend = lambda x: np.array([x[0] + dt * (x[1] - dt * kp * np.sin(x[0])), x[1] - dt * kp * np.sin(x[0])])
Jpend = lambda x: np.array([[1 - dt ** 2 * kp * np.cos(x[0]), dt], [-dt * kp * np.cos(x[0]), 1.0]])
n, sd_m = 400, 0.2
Qp = np.diag([1e-6, 1e-4])
xp = np.zeros((n, 2)); xp[0] = [2.0, 0.0]
for t in range(1, n):
    xp[t] = fpend(xp[t - 1]) + np.sqrt(np.diag(Qp)) * rng.standard_normal(2)
yp = xp[:, 0] + sd_m * rng.standard_normal(n)


def ekf(y, f, J, Q, Hv, a0, P0):
    """Extended Kalman filter, y_t = alpha_t[0] + eps_t, alpha_{t+1} = f(alpha_t) + eta_t (eq. 3.29)."""
    a, P, out = a0.astype(float), P0.astype(float), []
    for yt in y:
        F = P[0, 0] + Hv
        K = P[:, 0] / F
        a = a + K * (yt - a[0]); P = P - np.outer(K, K) * F
        out.append(a.copy())
        A = J(a); a = f(a); P = A @ P @ A.T + Q
    return np.array(out)


ae = ekf(yp, fpend, Jpend, Qp, sd_m ** 2, np.array([1.5, 0.0]), np.eye(2))
Tlin = Jpend(np.zeros(2))
al = ekf(yp, lambda x: Tlin @ x, lambda x: Tlin, Qp, sd_m ** 2, np.array([1.5, 0.0]), np.eye(2))
r_e, r_l = rmse(ae[50:, 0] - xp[50:, 0]), rmse(al[50:, 0] - xp[50:, 0])
print(f"Pendulum (theta0=2 rad): EKF angle RMSE {r_e:.3f}, linearized-at-0 KF {r_l:.3f}, raw {rmse(yp-xp[:,0]):.3f}")
assert r_e < 0.5 * r_l and r_e < sd_m
fig, ax = plt.subplots(figsize=(8, 3.2))
ax.plot(yp, ".", ms=2, color="0.6"); ax.plot(xp[:, 0], "k"); ax.plot(ae[:, 0], "C0", label="EKF")
ax.plot(al[:, 0], "C3", lw=0.7, label="linear KF"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(FIG / "ekf.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 12: GARCH(1,1) (§3.15-3.17, Example 4.8)


# %%
def garch_var(e, om_, a, b, s1):
    """sigma_t^2 = om + a e_{t-1}^2 + b sigma_{t-1}^2 (eq. 3.32), sigma_1^2 = s1, via lfilter."""
    u = om_ + a * e[:-1] ** 2
    rest = signal.lfilter([1.0], [1.0, -b], u, zi=[b * s1])[0]
    return np.r_[s1, rest]


def garch_nll(par, e):
    """Gaussian quasi-negative-log-likelihood (eq. 3.35), up to a constant."""
    om_, a, b = par
    if om_ <= 0 or a < 0 or b < 0 or a + b >= 0.9999:
        return 1e10
    s2 = garch_var(e, om_, a, b, e.var())
    return 0.5 * np.sum(np.log(s2) + e ** 2 / s2)


omg, ag, bg = 0.01, 0.1, 0.85
n = 6000
e = np.zeros(n); s2 = np.zeros(n); s2[0] = omg / (1 - ag - bg)
z = rng.standard_normal(n)
e[0] = np.sqrt(s2[0]) * z[0]
for t in range(1, n):
    s2[t] = omg + ag * e[t - 1] ** 2 + bg * s2[t - 1]
    e[t] = np.sqrt(s2[t]) * z[t]
acf = lambda x, L: np.array([np.corrcoef(x[h:], x[:-h])[0, 1] for h in range(1, L + 1)])
r_e, r_e2 = acf(e, 10), acf(e ** 2, 10)
kurt = stats.kurtosis(e, fisher=False)
kurt_th = 3 * (1 - (ag + bg) ** 2) / (1 - (ag + bg) ** 2 - 2 * ag ** 2)
print(f"Stylized facts: max|ACF(e)| lags 1-10 {np.abs(r_e).max():.3f} (band {2/np.sqrt(n):.3f}); "
      f"ACF(e^2) lag1 {r_e2[0]:.3f}, lag10 {r_e2[9]:.3f}; kurtosis {kurt:.2f} (theory {kurt_th:.3f}); "
      f"var {e.var():.4f} (theory {omg/(1-ag-bg):.2f})")
assert np.abs(r_e).max() < 3 / np.sqrt(n) and r_e2[0] > 0.1 and r_e2[9] > 0.05 and kurt > 3.2
res = optimize.minimize(garch_nll, [0.05, 0.05, 0.5], args=(e,), method="Nelder-Mead",
                        options=dict(xatol=1e-7, fatol=1e-7, maxiter=4000))
gh = res.x
print(f"GARCH QML: omega={gh[0]:.4f}, alpha={gh[1]:.4f}, beta={gh[2]:.4f}, persistence {gh[1]+gh[2]:.4f}")
assert abs(gh[1] - ag) < 0.04 and abs(gh[2] - bg) < 0.06 and abs(gh[1] + gh[2] - 0.95) < 0.03

sbar, pers = omg / (1 - ag - bg), ag + bg
hh = np.arange(1, 61)
fc = sbar + pers ** (hh - 1) * (0.05 - sbar)            # eq. (3.36), sigma^2_{t+1} = 0.05
print(f"Ex 4.8: sigma_bar^2={sbar:.3f}, persistence={pers:.2f}, 10-step forecast={fc[9]:.4f}, "
      f"half-life={np.log(0.5)/np.log(pers):.1f} steps")
assert abs(sbar - 0.2) < 1e-12 and abs(fc[9] - (0.2 - 0.15 * 0.95 ** 9)) < 1e-12
# Monte Carlo check of the forecast: simulate 20000 paths from sigma^2_{t+1}=0.05
S = np.full(20000, 0.05)
mc = []
for h in range(10):
    mc.append(S.mean())
    eps = np.sqrt(S) * rng.standard_normal(S.size)
    S = omg + ag * eps ** 2 + bg * S
assert abs(mc[9] - fc[9]) < 0.005
# ARCH(1) kurtosis 3(1-a^2)/(1-3a^2) for a=0.3
k_arch = 3 * (1 - 0.09) / (1 - 0.27)
print(f"ARCH(1) alpha=0.3 kurtosis {k_arch:.4f}")
fig, ax = plt.subplots(1, 3, figsize=(13, 3.4))
ax[0].plot(e[:1500], lw=0.5); ax[0].set_title("GARCH returns")
ax[1].bar(np.arange(1, 11) - 0.2, r_e, 0.4, label="ACF e"); ax[1].bar(np.arange(1, 11) + 0.2, r_e2, 0.4, label="ACF e^2")
ax[1].axhline(2 / np.sqrt(n), color="k", ls=":"); ax[1].legend()
ax[2].plot(hh, fc, label="forecast"); ax[2].plot(np.arange(1, 11), mc, "o", label="Monte Carlo")
ax[2].axhline(sbar, color="k", ls="--"); ax[2].legend()
fig.tight_layout(); fig.savefig(FIG / "garch.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 13: numbers used in the exercises and solutions

# %%
print("A1 MA(1) theta=-0.5: f(0), f(1/4), f(1/2) =", np.round(arma_spectrum([0, .25, .5], (), [-0.5]), 4))
c2 = 0.5 * (0.3 - 1) / (4 * 0.3)
print(f"A2 AR(2) 0.5,0.3: f(0)={float(arma_spectrum(0, [0.5, 0.3])):.4f}, f(1/2)={float(arma_spectrum(0.5, [0.5, 0.3])):.4f}, "
      f"stationary point cos={c2:.4f} at w={np.arccos(c2)/(2*np.pi):.4f} with f={float(arma_spectrum(np.arccos(c2)/(2*np.pi), [0.5, 0.3])):.4f} (a minimum)")
print("A3 chi2_14 quantiles", np.round(stats.chi2.ppf([.025, .975], 14), 3))
print("A3 Daniell m=3, nu=14: CI factors", round(14 / stats.chi2.ppf(.975, 14), 4), round(14 / stats.chi2.ppf(.025, 14), 4))
P = (0.25 + np.sqrt(0.0625 + 4 * 0.25)) / 2
print(f"A5 local level q=0.25: P_ss={P:.4f}, K_ss={P/(P+1):.4f}")
for q in (0.01, 100):
    P = (q + np.sqrt(q * q + 4 * q)) / 2
    print(f"   q={q}: P_ss={P:.4f}, K_ss={P/(P+1):.4f}")
print("A6 GARCH 0.02,0.08,0.9: var", round(0.02 / 0.02, 3), "half-life", round(np.log(.5) / np.log(.98), 2),
      "5-step from 2.0:", round(1 + 0.98 ** 4 * (2 - 1), 4))
# B5 two-step tracker continuation: y_2 = 2.5
a21, P21 = Tcv @ kf1["af"][0], Tcv @ kf1["Pf"][0] @ Tcv.T + Rcv @ Rcv.T
kf2 = kalman_filter(np.array([2.5]), Tcv, Zcv, Rcv, Qcv, Hcv, a21, P21)
print(f"B5 step 2: a_2|1={np.round(a21,4)}, P_2|1={np.round(P21,4).tolist()}, F={kf2['F'][0,0,0]:.4f}, "
      f"K={np.round(kf2['K'][0,:,0],4)}, a_2|2={np.round(kf2['af'][0],4)}")
kf12 = kalman_filter(np.array([2.0, 2.5]), Tcv, Zcv, Rcv, Qcv, Hcv, a10, P10)
assert np.allclose(kf12["af"][1], kf2["af"][0])
s1, s2_ = 1.0, 4.0
print(f"C1 static fusion: weights {s2_/(s1+s2_):.2f}, {s1/(s1+s2_):.2f}; fused var {1/(1/s1+1/s2_):.2f}")
print(f"B8 ARCH(1) alpha=0.5: kurtosis {3*(1-.25)/(1-.75):.1f}; alpha>=1/sqrt(3)={1/np.sqrt(3):.4f} infinite")

r1 = ag * (1 - ag * bg - bg ** 2) / (1 - 2 * ag * bg - bg ** 2)
print(f"B9 GARCH(0.1,0.85) ACF of e^2: rho(1)={r1:.4f}, rho(2)={r1*pers:.4f}, rho(10)={r1*pers**9:.4f} "
      f"(simulated lag1 {r_e2[0]:.3f}, lag10 {r_e2[9]:.3f})")
r1_s, r10_s = [], []
for _ in range(40):                        # Monte Carlo of the e^2 ACF, shorter paths: heavy tails make it noisy
    zz = rng.standard_normal(4000); ee = np.zeros(4000); ss = sbar
    for t_ in range(4000):
        ee[t_] = np.sqrt(ss) * zz[t_]; ss = omg + ag * ee[t_] ** 2 + bg * ss
    r1_s.append(acf(ee ** 2, 1)[0])
print(f"   mean simulated rho(1) over 40 paths of 4000: {np.mean(r1_s):.3f}; single paths scatter widely because e^4 is heavy-tailed")
# C2: exponential-kernel GP regression = Kalman smoother on an AR(1) state
ell, sg2, nz = 3.0, 1.0, 0.25
yg = rng.standard_normal(50)
Kg = sg2 * np.exp(-np.abs(np.subtract.outer(np.arange(50), np.arange(50))) / ell)
gp_mean = Kg @ np.linalg.solve(Kg + nz * np.eye(50), yg)
phg = np.exp(-1 / ell)
kfg = kalman_filter(yg, np.array([[phg]]), one, one, np.array([[sg2 * (1 - phg ** 2)]]), np.array([[nz]]),
                    np.zeros(1), np.array([[sg2]]))
ks_mean, _ = rts_smoother(kfg, np.array([[phg]]))
print(f"C2 GP (exponential kernel) posterior mean vs RTS smoother: max diff {np.abs(gp_mean - ks_mean[:, 0]).max():.1e}")
assert np.allclose(gp_mean, ks_mean[:, 0], atol=1e-10)

print("A2 roots of 1-0.5z-0.3z^2:", np.round(np.roots([-0.3, -0.5, 1]), 4))
print(f"A3 interval ratio m=3: {stats.chi2.ppf(.975,14)/stats.chi2.ppf(.025,14):.2f} vs m=2: {q_hi/q_lo:.2f}")
k12 = (np.sin(np.pi * 12 / 36) / (12 * np.sin(np.pi / 36))) ** 2
print(f"B2(c) 12-term MA at w=1/36 (3-year cycle in monthly data): |psi|^2={k12:.4f}")
print(f"B5 v_2={kf2['v'][0,0]:.4f}")
pq = (1 + np.sqrt(5)) / 2; th1 = pq / (pq + 1) - 1
print(f"B7 q=1: theta=K-1={th1:.4f}, theta/(1+theta^2)={th1/(1+th1**2):.4f}, -1/(q+2)={-1/3:.4f}")
assert abs(th1 / (1 + th1 ** 2) + 1 / 3) < 1e-12
m8 = 105 * ag ** 4 + 60 * ag ** 3 * bg + 18 * ag ** 2 * bg ** 2 + 4 * ag * bg ** 3 + bg ** 4
print(f"B9(c) fourth-moment condition (a+b)^2+2a^2={pers**2+2*ag**2:.4f}; eighth: E(a z^2+b)^4={m8:.4f}")
Z2 = np.array([[1.0, 0.0], [1.0, 0.0]])
kfc = kalman_filter(np.tile(rng.standard_normal((200, 1)), 2), Tcv, Z2, Rcv, np.eye(1), np.diag([1.0, 16.0]),
                    np.zeros(2), np.eye(2) * 100)
Kc = kfc["K"][-1]
print(f"C1(c) two-sensor steady-state gain: {np.round(Kc, 4).tolist()}; column ratio {Kc[0,0]/Kc[0,1]:.4f}")
assert np.allclose(Kc[:, 0], 16 * Kc[:, 1])
k4 = lambda w: (np.sin(np.pi * 4 * w) / (4 * np.sin(np.pi * w))) ** 2
print(f"B2 4-pt MA |psi|^2 at w=1/8, 1/4, 1/2: {k4(1/8):.4f}, {k4(0.25):.4f}, {k4(0.5):.4f}")


def bootstrap_pf(y, f_prop, loglik_obs, x0, N, rng_):
    """Bootstrap particle filter (§3.14): propagate, weight by p(y_t | x_t), resample.
    Returns the log-likelihood estimate sum_t log mean_i w_t^(i)."""
    x, ll = x0.copy(), 0.0
    for yt in y:
        lw = loglik_obs(yt, x)
        mx = lw.max()
        w = np.exp(lw - mx)
        ll += mx + np.log(w.mean())
        x = x[rng_.choice(N, N, p=w / w.sum())]
        x = f_prop(x, rng_)
    return ll


rpf = np.random.default_rng(3)
ypf = np.cumsum(rpf.standard_normal(100)) + rpf.standard_normal(100)
ll_kf = kalman_filter(ypf, one, one, one, one, one, np.zeros(1), np.array([[4.0]]))["ll"]
lls = [bootstrap_pf(ypf, lambda x, r: x + r.standard_normal(x.size),
                    lambda yt, x: -0.5 * np.log(2 * np.pi) - 0.5 * (yt - x) ** 2,
                    2.0 * rpf.standard_normal(2000), 2000, rpf) for _ in range(5)]
print(f"C3 particle filter on a local level model: log-lik {np.mean(lls):.2f} (sd {np.std(lls):.2f} over 5 runs) "
      f"vs Kalman {ll_kf:.2f}")
assert abs(np.mean(lls) - ll_kf) < 1.0

# %% [markdown]
# ## Your turn
# 1. Coherence: build x_t = c_t + noise, y_t = c_{t-3} filtered + noise with a common AR(1)
#    component c_t. Estimate the squared coherence (Part 5's code) and the phase; the phase
#    slope should be -2 pi * 3 (a pure delay). Where is the coherence high, and why?
# 2. Sensor fusion (exercise C1): observe one constant-velocity object with two sensors,
#    H = diag(1, 16). Stack Z = [[1,0],[1,0]] and run `kalman_filter`. Print the steady-state
#    gain columns: the precise sensor gets exactly 16 times the weight (C1 proves why).
# 3. Stochastic volatility by a particle filter: h_t = mu + phi(h_{t-1}-mu) + eta_t,
#    y_t = exp(h_t/2) z_t. Propagate 2000 particles, weight by N(y_t; 0, e^{h_t}), resample
#    when the effective sample size drops below 1000 (ST-6), and accumulate the log-likelihood
#    log mean(weights), reusing `bootstrap_pf` above. Profile it over phi on a grid.

# %%
if __name__ == "__main__":
    elapsed = time.time() - T0
    print(f"All checks passed in {elapsed:.1f} s. Figures in {FIG}")
