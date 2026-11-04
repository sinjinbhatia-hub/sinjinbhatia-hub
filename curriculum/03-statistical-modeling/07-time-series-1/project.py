# %% [markdown]
# # ST-7 Project: linear time-series models from scratch
#
# Everything here is numpy/scipy: simulation of AR/MA/ARMA/random walks, sample ACF and PACF, the companion
# matrix, Yule-Walker vs conditional least squares vs exact Gaussian ML, Durbin-Levinson, forecasting and
# coverage, ARIMA vs trend regression, Ljung-Box, the Dickey-Fuller distribution and power, spurious regression,
# Engle-Granger cointegration, and simple exponential smoothing as ARIMA(0,1,1). Every number in notes.md
# Section 4 and solutions.md is asserted or printed (Part 0 for the worked examples, Part 11 for the exercises).
# Run: `MPLBACKEND=Agg python3 project.py`.

# %%
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats, optimize, linalg, signal

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(7)
T0 = time.time()


def save(fig, name):
    fig.tight_layout(); fig.savefig(FIG / name, dpi=110); plt.close(fig)


# %% [markdown]
# ## Shared tools (notes 3.3-3.9)
# %%
def psi_weights(phi=(), theta=(), J=400):
    """psi_j of x = sum psi_j w_{t-j} from phi(B) psi(B) = theta(B), eq. (3.20): psi_j = theta_j + sum_k phi_k psi_{j-k}.
    Computed as the impulse response of the filter theta(B)/phi(B) (scipy.signal.lfilter runs that recursion)."""
    imp = np.zeros(J); imp[0] = 1.0
    return signal.lfilter(np.r_[1.0, theta], np.r_[1.0, -np.asarray(phi, float)], imp)


def arma_acov(phi=(), theta=(), sigma2=1.0, H=30):
    """gamma(h) = sigma^2 sum_j psi_j psi_{j+h}, eq. (3.6), for h = 0..H."""
    p = psi_weights(phi, theta, 2000)
    return sigma2 * np.array([p[: len(p) - h] @ p[h:] for h in range(H + 1)])


def simulate(phi=(), theta=(), n=500, burn=500, sigma=1.0, R=1, rng=rng):
    """R independent paths of phi(B) x = theta(B) w, w ~ N(0, sigma^2), started at 0 and burned in."""
    p, q = len(phi), len(theta); N = n + burn; w = sigma * rng.standard_normal((R, N)); x = np.zeros((R, N))
    for t in range(N):
        x[:, t] = w[:, t] + sum(theta[k] * w[:, t - k - 1] for k in range(q) if t - k - 1 >= 0) \
            + sum(phi[k] * x[:, t - k - 1] for k in range(p) if t - k - 1 >= 0)
    return x[:, burn:] if R > 1 else x[0, burn:]


def sacf(x, H):
    """Sample ACF rho_hat(h) = gamma_hat(h)/gamma_hat(0), gamma_hat with divisor n, eq. (3.4)."""
    x = np.asarray(x, float) - np.mean(x); n = len(x)
    g = np.array([x[h:] @ x[: n - h] / n for h in range(H + 1)])
    return g / g[0]


def durbin_levinson(rho, H):
    """PACF phi_hh for h = 1..H from autocorrelations rho(0..H), recursion (3.19)."""
    pacf = np.zeros(H); a = np.zeros(0)
    for h in range(1, H + 1):
        k = (rho[h] - a @ rho[h - 1:0:-1]) / (1 - a @ rho[1:h]) if h > 1 else rho[1]
        a = np.r_[a - k * a[::-1], k]; pacf[h - 1] = k
    return pacf


def pacf_regression(gam, H):
    """PACF by definition: last coefficient of the best linear predictor of x_t from x_{t-1..t-h} (Toeplitz solve)."""
    return np.array([linalg.solve(linalg.toeplitz(gam[:h]), gam[1:h + 1])[-1] for h in range(1, H + 1)])


def companion(phi):
    """Companion matrix of the AR(p) recursion, eq. (3.11)."""
    p = len(phi); A = np.zeros((p, p)); A[0] = phi; A[1:, :-1] += np.eye(p - 1)
    return A


def yule_walker(x, p):
    """phi_hat = Gamma_hat_p^{-1} gamma_hat_p, sigma2_hat = gamma_hat(0) - phi_hat . gamma_hat_p, eq. (3.23)."""
    x = np.asarray(x) - np.mean(x); n = len(x); g = np.array([x[h:] @ x[: n - h] / n for h in range(p + 1)])
    ph = linalg.solve(linalg.toeplitz(g[:p]), g[1:]); return ph, g[0] - ph @ g[1:]


def cls_ar(x, p):
    """Conditional least squares: OLS of x_t on (x_{t-1}, ..., x_{t-p}) (mean removed first)."""
    x = np.asarray(x) - np.mean(x); X = np.column_stack([x[p - k - 1: len(x) - k - 1] for k in range(p)])
    ph, *_ = np.linalg.lstsq(X, x[p:], rcond=None); e = x[p:] - X @ ph
    return ph, e @ e / len(e), e


def ar_loglik_exact(phi, sigma2, x):
    """Exact Gaussian log-likelihood of a zero-mean AR(p): MVN of the first p values (Toeplitz Gamma_p) times the
    conditional densities of the rest, eq. (3.25)."""
    p = len(phi); A = companion(phi); Q = np.zeros((p, p)); Q[0, 0] = sigma2
    Gp = linalg.solve_discrete_lyapunov(A, Q)                     # Gamma_p solves Gamma = A Gamma A^T + sigma^2 e1 e1^T
    c, low = linalg.cho_factor(Gp); x0 = x[:p]
    l0 = -0.5 * p * np.log(2 * np.pi) - np.log(np.diag(c)).sum() - 0.5 * x0 @ linalg.cho_solve((c, low), x0)
    X = np.column_stack([x[p - k - 1: len(x) - k - 1] for k in range(p)]); e = x[p:] - X @ np.asarray(phi)
    return l0 + stats.norm.logpdf(e, scale=np.sqrt(sigma2)).sum()


def ar_loglik_mvn(phi, sigma2, x):
    """The same likelihood as one n-dimensional MVN with Toeplitz covariance, via Cholesky (LA-10 Theorem 3.5)."""
    n = len(x); L = linalg.cholesky(linalg.toeplitz(arma_acov(phi, (), sigma2, n - 1)), lower=True)
    z = linalg.solve_triangular(L, x, lower=True)
    return -0.5 * n * np.log(2 * np.pi) - np.log(np.diag(L)).sum() - 0.5 * z @ z


def ar_mle(x, p):
    """Exact ML for a stationary AR(p): maximize (3.25) over (phi, log sigma2), starting from Yule-Walker."""
    x = np.asarray(x) - np.mean(x); ph0, s0 = yule_walker(x, p)
    def nll(v):
        ph = v[:p]
        if np.max(np.abs(np.linalg.eigvals(companion(ph)))) >= 0.999: return 1e10
        return -ar_loglik_exact(ph, np.exp(v[p]), x)
    r = optimize.minimize(nll, np.r_[ph0, np.log(s0)], method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-8, maxiter=4000))
    return r.x[:p], np.exp(r.x[p])


def ljung_box(e, K, fitted):
    """Q = n(n+2) sum_{h<=K} rho_hat(h)^2/(n-h), compared with chi^2_{K-fitted}, eq. (3.27). Returns (Q, p-value)."""
    n = len(e); r = sacf(e, K)[1:]; Q = n * (n + 2) * np.sum(r**2 / (n - np.arange(1, K + 1)))
    return Q, stats.chi2.sf(Q, K - fitted)


def ols_t(y, X):
    """OLS slope t statistics for many regressions at once: y, X have shape (R, n, ...)."""
    b = np.linalg.solve(np.einsum("rni,rnj->rij", X, X), np.einsum("rni,rn->ri", X, y)[..., None])[..., 0]
    e = y - np.einsum("rni,ri->rn", X, b); s2 = (e**2).sum(1) / (X.shape[1] - X.shape[2])
    V = np.linalg.inv(np.einsum("rni,rnj->rij", X, X)) * s2[:, None, None]
    return b, b / np.sqrt(np.diagonal(V, axis1=1, axis2=2)), e


def df_t(x, const=True):
    """Dickey-Fuller t statistic for rho in dx_t = (alpha) + rho x_{t-1} + e_t, for each row of x, eq. (3.33)."""
    dx, lag = np.diff(x, axis=1), x[:, :-1]
    X = np.stack([np.ones_like(lag), lag], -1) if const else lag[..., None]
    _, t, _ = ols_t(dx, X); return t[:, -1]


# %% [markdown]
# ## Part 0: the worked examples of notes Section 4
# %%
# Ex 4.1 AR(1) phi=0.7 and MA(1) theta=0.5
g = arma_acov((0.7,), (), 1, 5); r1 = g / g[0]
print("Ex 4.1 AR(1) gamma(0) =", round(g[0], 5), " rho(0..5) =", np.round(r1, 5))
assert np.allclose(r1, 0.7 ** np.arange(6)) and np.isclose(g[0], 1 / 0.51)
gm = arma_acov((), (0.5,), 1, 3); print("       MA(1) gamma(0..3) =", gm, " rho(1) =", gm[1] / gm[0])
assert np.allclose(gm, [1.25, 0.5, 0, 0]) and np.isclose(gm[1] / gm[0], 0.4)
assert np.allclose(arma_acov((), (2.0,), 0.25, 3), gm)                         # theta -> 1/theta, sigma2 -> theta^2 sigma2
# Ex 4.2 AR(2) phi=(1,-0.5): roots, companion eigenvalues, quasi-period
z = np.roots([0.5, -1.0, 1.0]); lam = np.linalg.eigvals(companion((1.0, -0.5)))
print("Ex 4.2 roots of phi(z):", z, "|z| =", np.abs(z), " eigenvalues:", lam, "|lam| =", np.abs(lam))
assert np.allclose(np.sort_complex(1 / z), np.sort_complex(lam)) and np.allclose(np.abs(z), np.sqrt(2))
period = 2 * np.pi / np.abs(np.angle(lam[0])); r2 = arma_acov((1.0, -0.5), (), 1, 8); r2n = r2 / r2[0]
print(f"       period = {period:.4f};  gamma(0) = {r2[0]:.4f};  rho(0..8) =", np.round(r2n, 4))
assert np.isclose(period, 8) and np.isclose(r2[0], 2.4) and np.allclose(r2n[1:5], [2 / 3, 1 / 6, -1 / 6, -1 / 4])
hh = np.arange(9); assert np.allclose(r2n, 2 ** (-hh / 2) * (np.cos(np.pi * hh / 4) + np.sin(np.pi * hh / 4) / 3))   # closed form
g0_ar2 = lambda p1, p2: (1 - p2) / ((1 + p2) * ((1 - p2) ** 2 - p1**2))                                   # eq. (3.16)
assert np.isclose(g0_ar2(1.0, -0.5), 2.4) and np.isclose(g0_ar2(1.5, -0.75), arma_acov((1.5, -0.75), (), 1, 0)[0])
# Ex 4.3 Yule-Walker from gamma = (2, 1.2, 0.4)
G = np.array([2.0, 1.2, 0.4]); phYW = linalg.solve(linalg.toeplitz(G[:2]), G[1:]); s2YW = G[0] - phYW @ G[1:]
print("Ex 4.3 YW phi =", phYW, " sigma2 =", s2YW); assert np.allclose(phYW, [0.75, -0.25]) and np.isclose(s2YW, 1.2)
assert np.allclose(arma_acov(phYW, (), s2YW, 2), G)
assert np.allclose(np.abs(np.linalg.eigvals(companion(phYW))), 0.5) and np.isclose(phYW[0]**2 + 4 * phYW[1], -0.4375)   # complex eigenvalues of modulus sqrt(0.25)        # the fitted model reproduces the three autocovariances
# Ex 4.4 PACF of AR(1) phi=0.7 by regression
pa = pacf_regression(arma_acov((0.7,), (), 1, 4), 3); print("Ex 4.4 PACF(1..3) =", np.round(pa, 12))
assert np.allclose(pa, [0.7, 0, 0])
# Ex 4.5 ARMA(1,1) phi=0.5, theta=0.4
ps = psi_weights((0.5,), (0.4,), 5); ga = arma_acov((0.5,), (0.4,), 1, 3)
print("Ex 4.5 psi(0..4) =", ps, " gamma(0..3) =", np.round(ga, 5), " rho(1..3) =", np.round(ga[1:] / ga[0], 5))
print("       AR(1) with the same rho(1) would have rho(2) =", round((1.08 / 1.56) ** 2, 4))
assert np.allclose(ps, [1, 0.9, 0.45, 0.225, 0.1125]) and np.allclose(ga[:2], [2.08, 1.44])
assert np.isclose(ga[1] / ga[0], 1.08 / 1.56) and np.allclose(ga[2:] / ga[1:-1], 0.5)
# Ex 4.6 AR(1) phi=0.8 forecasts from x_n = 2
for h in (1, 5):
    f = 0.8**h * 2; mse = (1 - 0.8 ** (2 * h)) / (1 - 0.64); hw = 1.96 * np.sqrt(mse)
    print(f"Ex 4.6 h={h}: forecast {f:.5f}, MSE {mse:.5f}, 95% PI [{f - hw:.4f}, {f + hw:.4f}]")
assert np.isclose(0.8**5 * 2, 0.65536) and np.isclose((1 - 0.8**10) / 0.36, 2.479516, atol=1e-6)
print(f"       gamma(0) = {1 / 0.36:.4f}; sqrt(P5) = {np.sqrt(2.479516):.5f}; half-width {1.96 * np.sqrt(2.479516):.5f};"
      f" width ratio to the limit {np.sqrt(2.479516 * 0.36):.4f}")
# Ex 4.7 exact AR(1) likelihood for x = (0.5, 1.0, -0.2), phi = 0.6, sigma2 = 1
x3 = np.array([0.5, 1.0, -0.2]); l_fact = ar_loglik_exact((0.6,), 1.0, x3); l_mvn = ar_loglik_mvn((0.6,), 1.0, x3)
terms = [stats.norm.logpdf(0.5, scale=np.sqrt(1 / 0.64)), stats.norm.logpdf(1.0 - 0.3), stats.norm.logpdf(-0.2 - 0.6)]
print("Ex 4.7 log-lik terms", np.round(terms, 5), " total", round(l_fact, 5), " (MVN/Cholesky:", round(l_mvn, 5), ")")
assert np.isclose(l_fact, l_mvn) and np.isclose(l_fact, sum(terms))
print(f"       pieces: {-0.5 * np.log(2 * np.pi * 1.5625):.5f} - 0.08;  {-0.5 * np.log(2 * np.pi):.5f} - 0.245;  ... - 0.32")
# Ex 4.8 SES alpha=0.3 = ARIMA(0,1,1) theta=-0.7
xs = np.array([10.0, 12, 11, 13, 12]); a = 0.3; lev = [xs[0]]
for v in xs[1:]: lev.append(a * v + (1 - a) * lev[-1])
th = a - 1; fc = [xs[0]]                                          # ARIMA(0,1,1): xhat_{t+1} = x_t + theta (x_t - xhat_t)
for t in range(1, len(xs)): fc.append(xs[t] + th * (xs[t] - fc[-1]))
print("Ex 4.8 SES levels", np.round(lev, 5), " ARIMA(0,1,1) forecasts", np.round(fc, 5))
assert np.allclose(lev, fc) and np.isclose(lev[-1], 11.5828)

# %% [markdown]
# ## Part 1: the identification table made visual (notes 3.7)
# Nine processes, each with its sample ACF and PACF and +-2/sqrt(n) bands. Theoretical ACFs of the stationary ones
# are asserted against the sample ACF within 4 Bartlett standard errors, eq. (3.5) generalized (Bartlett's formula).
# %%
n = 2000; H = 20
procs = {"AR(1) 0.7": ((0.7,), ()), "AR(2) real (0.5,0.3)": ((0.5, 0.3), ()), "AR(2) complex (1,-0.5)": ((1.0, -0.5), ()),
         "MA(1) 0.5": ((), (0.5,)), "MA(2) (0.6,0.3)": ((), (0.6, 0.3)), "ARMA(1,1) (0.5,0.4)": ((0.5,), (0.4,))}
series = {k: simulate(*v, n=n) for k, v in procs.items()}
w_ = rng.standard_normal(n); series["random walk"] = np.cumsum(w_); series["RW + drift 0.2"] = np.cumsum(0.2 + w_)
series["white noise"] = rng.standard_normal(n)
fig, ax = plt.subplots(len(series), 3, figsize=(11, 2.0 * len(series)))
for i, (k, x) in enumerate(series.items()):
    r = sacf(x, H); pc = durbin_levinson(r, H); band = 2 / np.sqrt(n)
    ax[i, 0].plot(x[:300], lw=0.7); ax[i, 0].set_ylabel(k, fontsize=8)
    ax[i, 1].bar(range(1, H + 1), r[1:], width=0.4); ax[i, 2].bar(range(1, H + 1), pc, width=0.4, color="C1")
    for j in (1, 2): ax[i, j].axhline(band, ls="--", c="k", lw=0.6); ax[i, j].axhline(-band, ls="--", c="k", lw=0.6)
    if k in procs:
        g = arma_acov(*procs[k], 1, H + 60); rt = g / g[0]; ax[i, 1].plot(range(1, H + 1), rt[1:H + 1], "r.")
        ext = np.r_[rt[::-1], rt[1:]]                                # rho(-H-60..H+60)
        for h in range(1, 11):                                       # Bartlett: var = (1/n) sum_k (rho(k+h)+rho(k-h)-2rho(h)rho(k))^2
            ks = np.arange(1, 50); c = H + 60
            v = np.sum((ext[c + ks + h] + ext[c + ks - h] - 2 * rt[h] * rt[ks]) ** 2) / n
            assert abs(r[h] - rt[h]) < 4 * np.sqrt(v), (k, h)
ax[0, 0].set_title("series (first 300)"); ax[0, 1].set_title("sample ACF (red: theory)"); ax[0, 2].set_title("sample PACF")
save(fig, "identification_table.png")
r_rw = sacf(series["random walk"], 20); print("random walk sample ACF lags 1,5,10,20:", np.round(r_rw[[1, 5, 10, 20]], 3))
assert np.all(r_rw[1:] > 2 / np.sqrt(n))                           # 'significant' at every lag: the slow-decay diagnostic

# %% [markdown]
# ## Part 2: AR(2) stationarity triangle and companion eigenvalues (notes 3.5; LA-7 eq. (3.8))
# %%
for ph in [(0.5, 0.3), (1.0, -0.5), (1.5, -0.75), (1.2, -0.2), (0.3, 0.75)]:
    zz = np.roots(np.r_[-np.array(ph)[::-1], 1.0]); ll = np.linalg.eigvals(companion(ph))
    assert np.allclose(np.sort_complex(1 / zz), np.sort_complex(ll))          # eigenvalues = reciprocal roots
    tri = (ph[0] + ph[1] < 1) and (ph[1] - ph[0] < 1) and abs(ph[1]) < 1
    assert tri == (np.max(np.abs(ll)) < 1 - 1e-12), ph
G1, G2 = np.meshgrid(np.linspace(-2.2, 2.2, 400), np.linspace(-1.2, 1.2, 300))
fig, ax = plt.subplots(figsize=(6, 4))
inside = (G1 + G2 < 1) & (G2 - G1 < 1) & (np.abs(G2) < 1)
ax.contourf(G1, G2, inside & (G1**2 + 4 * G2 < 0), levels=[0.5, 1], colors=["C1"], alpha=0.5)
ax.contourf(G1, G2, inside & (G1**2 + 4 * G2 >= 0), levels=[0.5, 1], colors=["C0"], alpha=0.3)
ax.plot([-2, 0, 2, -2], [-1, 1, -1, -1], "k"); ax.plot(np.linspace(-2, 2), -np.linspace(-2, 2) ** 2 / 4, "k--", lw=0.8)
ax.plot([1.0, 0.5, 1.5], [-0.5, 0.3, -0.75], "ko"); ax.set_xlabel(r"$\phi_1$"); ax.set_ylabel(r"$\phi_2$")
ax.set_title("AR(2) stationary region (orange: complex roots)"); save(fig, "ar2_triangle.png")

# %% [markdown]
# ## Part 3: Yule-Walker vs CLS vs exact ML, 500 replications of AR(2) phi = (1.0, -0.5), n = 200 (notes 3.8)
# The exact likelihood (3.25) is checked against the n-dimensional Toeplitz MVN evaluated by Cholesky.
# %%
phT = np.array([1.0, -0.5]); nR, R = 200, 500; X = simulate(phT, n=nR, R=R)
assert np.isclose(ar_loglik_exact(phT, 1.0, X[0] - X[0].mean()), ar_loglik_mvn(phT, 1.0, X[0] - X[0].mean()))
est = {"YW": [], "CLS": [], "ML": []}
for x in X:
    est["YW"].append(yule_walker(x, 2)[0]); est["CLS"].append(cls_ar(x, 2)[0]); est["ML"].append(ar_mle(x, 2)[0])
est = {k: np.array(v) for k, v in est.items()}
Vasy = 1.0 * linalg.inv(linalg.toeplitz(arma_acov(phT, (), 1, 1)))    # asymptotic cov sigma^2 Gamma_2^{-1} / n, eq. (3.26)
print("asymptotic sd:", np.round(np.sqrt(np.diag(Vasy) / nR), 4))
for k, e in est.items():
    print(f"  {k:3s}: mean {np.round(e.mean(0), 4)}, sd {np.round(e.std(0), 4)}")
    assert np.allclose(e.mean(0), phT, atol=0.03) and np.allclose(e.std(0), np.sqrt(np.diag(Vasy) / nR), rtol=0.15)
    zst = (e[:, 0] - phT[0]) / np.sqrt(Vasy[0, 0] / nR); assert stats.kstest(zst - zst.mean(), "norm").pvalue > 0.01
dML, dYW = np.abs(est["ML"] - est["CLS"]), np.abs(est["YW"] - est["CLS"])
print(f"  |ML - CLS|: mean {dML.mean():.4f}, max {dML.max():.4f};  |YW - CLS|: mean {dYW.mean():.4f}, max {dYW.max():.4f}")
assert dML.mean() < 0.005 and dYW.mean() < 0.03 and dML.max() < 0.05
fig, ax = plt.subplots(1, 2, figsize=(10, 3.5))
for j in range(2):
    for k, e in est.items(): ax[j].hist(e[:, j], bins=30, histtype="step", density=True, label=k)
    gx = np.linspace(*ax[j].get_xlim(), 200); ax[j].plot(gx, stats.norm.pdf(gx, phT[j], np.sqrt(Vasy[j, j] / nR)), "k--")
    ax[j].set_title(rf"$\hat\phi_{j + 1}$ over 500 replications"); ax[j].legend()
save(fig, "estimators.png")

# %% [markdown]
# ## Part 4: PACF by successive regressions = Durbin-Levinson (notes 3.6)
# %%
for ph, th in [((1.0, -0.5), ()), ((0.5, 0.3), ()), ((), (0.5,)), ((0.5,), (0.4,))]:
    g = arma_acov(ph, th, 1, 12); assert np.allclose(pacf_regression(g, 12), durbin_levinson(g / g[0], 12))
g = arma_acov((1.0, -0.5), (), 1, 6); print("AR(2) PACF(1..6):", np.round(durbin_levinson(g / g[0], 6), 6))
assert np.allclose(durbin_levinson(g / g[0], 6)[2:], 0, atol=1e-10) and np.isclose(durbin_levinson(g / g[0], 6)[1], -0.5)
xs2 = X[0]; rs = sacf(xs2, 6)
ols_pacf = [cls_ar(xs2, h)[0][-1] for h in range(1, 7)]
assert np.allclose(durbin_levinson(rs, 6), [yule_walker(xs2, h)[0][-1] for h in range(1, 7)])   # exact on sample too
assert np.max(np.abs(np.array(ols_pacf) - durbin_levinson(rs, 6))) < 0.05                    # OLS version is close

# %% [markdown]
# ## Part 5: forecasting a fitted AR(2) and interval coverage (notes 3.9)
# %%
def ar_forecast(x, ph, s2, Hh):
    """h-step forecasts by iterating the recursion with future shocks set to 0, and MSE sigma^2 sum_{j<h} psi_j^2, eq. (3.30)."""
    hist = list(x[-len(ph):]); f = []
    for _ in range(Hh):
        nxt = sum(ph[k] * hist[-k - 1] for k in range(len(ph))); f.append(nxt); hist.append(nxt)
    ps = psi_weights(ph, (), Hh); return np.array(f), s2 * np.cumsum(ps**2)
Hh = 10; Y = simulate(phT, n=nR + Hh, R=R); cover = np.zeros(Hh)
for y in Y:
    m = y[:nR].mean(); ph, s2, _ = cls_ar(y[:nR], 2); f, mse = ar_forecast(y[:nR] - m, ph, s2, Hh)
    cover += np.abs(y[nR:] - m - f) < 1.96 * np.sqrt(mse)
cover /= R; print("AR(2) forecast coverage h=1..10:", np.round(cover, 3))
assert np.all(np.abs(cover - 0.95) < 0.035)
f, mse = ar_forecast(Y[0][:nR] - Y[0][:nR].mean(), *cls_ar(Y[0][:nR], 2)[:2], 30)
fig, ax = plt.subplots(figsize=(8, 3)); tt = np.arange(nR - 50, nR); tf = np.arange(nR, nR + 30)
ax.plot(tt, Y[0][nR - 50:nR] - Y[0][:nR].mean(), "k"); ax.plot(tf, f, "C0")
ax.fill_between(tf, f - 1.96 * np.sqrt(mse), f + 1.96 * np.sqrt(mse), alpha=0.3)
ax.set_title("AR(2) forecasts decay to the mean; the band widens to the stationary sd"); save(fig, "forecast.png")
assert np.isclose(mse[-1], cls_ar(Y[0][:nR], 2)[1] * np.sum(psi_weights(cls_ar(Y[0][:nR], 2)[0], (), 30) ** 2))

# %% [markdown]
# ## Part 6: ARIMA(0,1,0) with drift vs trend regression on a random walk with drift (notes 3.10)
# %%
nA, HA, RA, dl = 200, 20, 1000, 0.2; Wd = dl + rng.standard_normal((RA, nA + HA)); Xd = np.cumsum(Wd, axis=1)
cov_arima = cov_trend = 0
for x in Xd:
    d = np.diff(x[:nA]); dh, sd = d.mean(), d.std(ddof=1)
    fa = x[nA - 1] + dh * HA; cov_arima += abs(x[-1] - fa) < 1.96 * sd * np.sqrt(HA)
    tt = np.arange(nA); b = np.polyfit(tt, x[:nA], 1); res = x[:nA] - np.polyval(b, tt)
    ft = np.polyval(b, nA - 1 + HA); cov_trend += abs(x[-1] - ft) < 1.96 * res.std(ddof=2)
cov_arima /= RA; cov_trend /= RA
print(f"h={HA} interval coverage: ARIMA(0,1,0)+drift {cov_arima:.3f}, trend regression {cov_trend:.3f}")
assert abs(cov_arima - 0.95) < 0.03 and cov_trend < 0.85

# %% [markdown]
# ## Part 7: Ljung-Box calibration and power (notes 3.8)
# %%
K, RL = 10, 1000; Z = simulate(phT, n=300, R=RL); pv_ok, pv_under = [], []
for z in Z:
    pv_ok.append(ljung_box(cls_ar(z, 2)[2], K, 2)[1]); pv_under.append(ljung_box(cls_ar(z, 1)[2], K, 1)[1])
size, power = np.mean(np.array(pv_ok) < 0.05), np.mean(np.array(pv_under) < 0.05)
print(f"Ljung-Box K={K}: size (correct AR(2)) {size:.3f}, power (AR(1) under-fit) {power:.3f}")
assert 0.03 < size < 0.075 and power > 0.9

# %% [markdown]
# ## Part 8: the Dickey-Fuller distribution and the test's power (notes 3.11)
# %%
RD, nD = 10000, 100; RW = np.cumsum(rng.standard_normal((RD, nD)), axis=1)
t0, tc = df_t(RW, const=False), df_t(RW, const=True)
q0, qc = np.quantile(t0, [0.01, 0.05, 0.10]), np.quantile(tc, [0.01, 0.05, 0.10])
print("DF quantiles (1,5,10%) no constant:", np.round(q0, 3), " with constant:", np.round(qc, 3),
      " t_{98} 5%:", round(stats.t.ppf(0.05, nD - 2), 3))
assert abs(q0[1] + 1.95) < 0.08 and abs(qc[1] + 2.89) < 0.08 and qc[1] < stats.t.ppf(0.05, nD - 2) - 1
print(f"P(chi2_1 < 1) = {stats.chi2.cdf(1, 1):.4f}")
rhohat = (np.sum(RW[:, :-1] * np.diff(RW, axis=1), 1) / np.sum(RW[:, :-1] ** 2, 1))
print("n(phi_hat - 1) for phi = 1: mean", round(np.mean(nD * rhohat), 3), " 5% quantile", round(np.quantile(nD * rhohat, 0.05), 3))
pw = {}
for ph1 in (1.0, 0.99, 0.97, 0.95, 0.9, 0.8):
    xs_ = np.zeros((4000, nD)); e = rng.standard_normal((4000, nD)); xs_[:, 0] = e[:, 0] / np.sqrt(max(1 - ph1**2, 1e-3))
    for t in range(1, nD): xs_[:, t] = ph1 * xs_[:, t - 1] + e[:, t]
    pw[ph1] = np.mean(df_t(xs_, const=True) < qc[1])
print("DF power (5%, constant, n=100):", {k: round(float(v), 3) for k, v in pw.items()})
assert abs(pw[1.0] - 0.05) < 0.015 and pw[0.9] < 0.5 and pw[0.8] > 0.8
fig, ax = plt.subplots(figsize=(7, 3.5)); gx = np.linspace(-5, 3, 300)
ax.hist(t0, bins=80, density=True, histtype="step", label="DF, no constant"); ax.hist(tc, bins=80, density=True, histtype="step", label="DF, constant")
ax.plot(gx, stats.t.pdf(gx, nD - 2), "k--", label="t(98)"); ax.legend(); ax.set_title(r"t statistic for $\phi=1$, 10$^4$ random walks")
save(fig, "dickey_fuller.png")

# %% [markdown]
# ## Part 9: spurious regression and cointegration (notes 3.12)
# %%
RS, nS = 10000, 50
def frac_sig(xa, ya):
    Xm = np.stack([np.ones_like(xa), xa], -1); _, t, _ = ols_t(ya, Xm); return np.mean(np.abs(t[:, 1]) > 2), t[:, 1]
fr_rw, t_rw = frac_sig(np.cumsum(rng.standard_normal((RS, nS)), 1), np.cumsum(rng.standard_normal((RS, nS)), 1))
fr_ar, _ = frac_sig(simulate((0.5,), n=nS, R=RS), simulate((0.5,), n=nS, R=RS))
fr_iid, _ = frac_sig(rng.standard_normal((RS, nS)), rng.standard_normal((RS, nS)))
print(f"fraction |t|>2, n={nS}: random walks {fr_rw:.3f}, AR(1) 0.5 {fr_ar:.3f}, iid {fr_iid:.3f}")
assert fr_rw > 0.5 and 0.04 < fr_iid < 0.07 and fr_iid < fr_ar < fr_rw
# Engle-Granger: DF t on residuals of y on x; critical values simulated under no cointegration (independent RWs)
nE = 100; xa, ya = np.cumsum(rng.standard_normal((RD, nE)), 1), np.cumsum(rng.standard_normal((RD, nE)), 1)
_, _, eres = ols_t(ya, np.stack([np.ones_like(xa), xa], -1)); teg = df_t(eres, const=False); cEG = np.quantile(teg, 0.05)
print(f"Engle-Granger 5% critical value (n={nE}) {cEG:.3f}  vs DF {qc[1]:.3f}")
assert -3.6 < cEG < -3.2
xc = np.cumsum(rng.standard_normal((2000, nE)), 1); u = simulate((0.5,), n=nE, R=2000); yc = 1.0 + 2.0 * xc + u
bc, _, ec = ols_t(yc, np.stack([np.ones_like(xc), xc], -1)); pEG = np.mean(df_t(ec, const=False) < cEG)
print(f"cointegrated pair (y = 1 + 2x + AR(1) 0.5): mean beta_hat {bc[:, 1].mean():.4f}, EG rejection rate {pEG:.3f}")
assert abs(bc[:, 1].mean() - 2) < 0.01 and pEG > 0.9
size_wrong = np.mean(teg < qc[1]); print(f"EG with the plain DF critical value: size {size_wrong:.3f}")
fig, ax = plt.subplots(1, 2, figsize=(10, 3.2))
ax[0].hist(t_rw, bins=80, density=True); ax[0].plot(gx * 3, stats.norm.pdf(gx * 3), "k--")
ax[0].set_title("slope t, independent random walks (n=50)")
ax[1].plot(xc[0] * 2 + 1, label="1 + 2x"); ax[1].plot(yc[0], label="y"); ax[1].plot(ec[0] - 8, label="residual - 8")
ax[1].legend(fontsize=8); ax[1].set_title("cointegrated pair"); save(fig, "spurious_cointegration.png")

# %% [markdown]
# ## Part 10: simple exponential smoothing = ARIMA(0,1,1) on a long series (notes 3.13)
# %%
xi = np.cumsum(rng.standard_normal(300)) + simulate((), (), n=300)
for alpha in (0.1, 0.5, 0.9):
    l = xi[0]; f2 = xi[0]
    for t in range(1, 300):
        l = alpha * xi[t] + (1 - alpha) * l; f2 = xi[t] + (alpha - 1) * (xi[t] - f2)
        assert np.isclose(l, f2)
print("SES == ARIMA(0,1,1) with theta = alpha - 1: checked at alpha = 0.1, 0.5, 0.9")

# %% [markdown]
# ## Part 11: numbers for the exercises
# %%
gA1 = arma_acov((), (0.6, -0.3), 1, 3); print("A1 MA(2) gamma(0..3):", np.round(gA1, 4), " rho(1,2):", np.round(gA1[1:3] / gA1[0], 4))
assert np.allclose(gA1, [1.45, 0.42, -0.3, 0])
zA1 = np.roots([-0.3, 0.6, 1.0]); print("A1 roots of theta(z):", np.round(zA1, 4)); assert np.all(np.abs(zA1) > 1)
for ph in [(1.2, -0.2), (0.4, 0.45), (1.5, -0.75)]:
    print("A2", ph, "|roots| =", np.round(np.abs(np.roots(np.r_[-np.array(ph)[::-1], 1.0])), 4))
lmA2 = np.linalg.eigvals(companion((1.5, -0.75))); assert np.isclose(2 * np.pi / abs(np.angle(lmA2[0])), 12)
print("A3 corr(x_100, x_110) =", round(np.sqrt(100 / 110), 4))
print(f"A5 P(at least one of 20 white-noise lags outside the band) = {1 - 0.95**20:.4f}")
Q = 100 * 102 * (0.15**2 / 99 + 0.1**2 / 98 + 0.05**2 / 97); print(f"A6 Q = {Q:.4f}, chi2_3 95% = {stats.chi2.ppf(0.95, 3):.3f}, p = {stats.chi2.sf(Q, 3):.4f}")
print(f"A6 after an AR(1) fit: chi2_2 95% = {stats.chi2.ppf(0.95, 2):.3f}, p = {stats.chi2.sf(Q, 2):.4f}")
for nn in (4, 5, 6, 10):
    ev = np.linalg.eigvalsh(linalg.toeplitz(np.r_[1, 0.6, np.zeros(nn - 2)])); print(f"B1 n={nn}: min eigenvalue {ev[0]:.4f}")
assert np.linalg.eigvalsh(linalg.toeplitz([1, 0.6, 0, 0]))[0] > 0 and np.linalg.eigvalsh(linalg.toeplitz([1, 0.6, 0, 0, 0]))[0] < 0
gB2 = arma_acov((1.5, -0.75), (), 1, 4); print("B2 gamma(0) =", round(gB2[0], 4), " rho(1..4) =", np.round(gB2[1:] / gB2[0], 4))
rB3 = np.array([1, 0.8, 0.5, 0.3]); print("B3 PACF:", np.round(durbin_levinson(rB3, 3), 4))
assert np.allclose(durbin_levinson(rB3, 3), [0.8, -0.14 / 0.36, (1 / 18) / (11 / 36)])
ph3 = linalg.solve(linalg.toeplitz(rB3[:3]), rB3[1:]); pc3 = durbin_levinson(rB3, 3); v = np.cumprod(1 - pc3**2)
print("B3 AR(3) YW coefficients:", np.round(ph3, 4), " prediction error variances v1..v3 / gamma(0):", np.round(v, 4))
assert np.isclose(ph3[-1], pc3[-1])
fB4, mB4 = ar_forecast(np.array([2.0, 1.0]), (1.0, -0.5), 1.0, 3); print("B4 forecasts", fB4, " MSE", mB4)
assert np.allclose(fB4, [0, -0.5, -0.5]) and np.allclose(mB4, [1, 2, 2.25])
print("B4 95% intervals:", [(round(float(a - 1.96 * np.sqrt(b)), 4), round(float(a + 1.96 * np.sqrt(b)), 4)) for a, b in zip(fB4, mB4)])
print(f"B2 r = {np.sqrt(0.75):.4f}, omega = {np.arccos(1.5 / (2 * np.sqrt(0.75))):.4f} (pi/6 = {np.pi / 6:.4f})")
xb6 = np.array([1.0, 2, 4]); cls6 = (xb6[1:] @ xb6[:-1]) / (xb6[:-1] @ xb6[:-1]); yw6 = (xb6[1:] @ xb6[:-1]) / (xb6 @ xb6)
print(f"B6 x=(1,2,4): CLS {cls6:.4f}, YW {yw6:.4f}"); assert np.isclose(cls6, 2) and np.isclose(yw6, 10 / 21)
print("B5 theta from r=0.4:", (1 - np.sqrt(1 - 4 * 0.16)) / 0.8)
psB8 = np.cumsum(0.5 ** np.arange(3)); print("B8 psi:", psB8, " MSE:", np.cumsum(psB8**2))
assert np.allclose(np.cumsum(psB8**2), [1, 3.25, 6.3125])
yB8 = [9.0, 10.0]
for _ in range(3): yB8.append(yB8[-1] + 0.5 * (yB8[-1] - yB8[-2]))
print("B8 forecasts:", yB8[2:]); assert np.allclose(yB8[2:], [10.5, 10.75, 10.875])
for nn in (50, 100, 200, 400):
    f_, _ = frac_sig(np.cumsum(rng.standard_normal((4000, nn)), 1), np.cumsum(rng.standard_normal((4000, nn)), 1))
    print(f"B7 n={nn}: fraction |t|>2 = {f_:.3f}")
f_d, _ = frac_sig(rng.standard_normal((4000, 400)), rng.standard_normal((4000, 400)))   # differences of RWs are iid
print(f"B7 regress the differences (n=400): fraction |t|>2 = {f_d:.3f}"); assert 0.035 < f_d < 0.07
idw = RW[:, :-1] * np.diff(RW, axis=1); lhs = idw.sum(1); rhs = 0.5 * (RW[:, -1] ** 2 - RW[:, 0] ** 2 - (np.diff(RW, axis=1) ** 2).sum(1))
assert np.allclose(lhs, rhs); print("C1 identity sum x_{t-1} w_t = (x_n^2 - x_1^2 - sum w_t^2)/2 checked")
kap, dt_, sOU = 2.0, 0.1, 1.0; print(f"C3 OU -> AR(1): phi = {np.exp(-kap * dt_):.5f}, sigma2 = {sOU**2 * (1 - np.exp(-2 * kap * dt_)) / (2 * kap):.5f}")

# %% [markdown]
# ## Your turn
# 1. **Seasonal ARIMA.** Simulate monthly data from (1 - B^12) x_t = (1 + Theta B^12) w_t with Theta = -0.6 plus a
#    mean seasonal pattern. Show the sample ACF of x has spikes at 12, 24, 36 that decay slowly; seasonally difference
#    and show the ACF of nabla_12 x cuts off after lag 12 (an MA at lag 12). Fit Theta by conditional least squares.
# 2. **Yule's sunspots.** Simulate Yule's model directly: an AR(2) with complex roots, modulus 0.9 and period 11
#    (phi_1 = 2 (0.9) cos(2 pi/11), phi_2 = -0.81). Fit an AR(2) by Yule-Walker on n = 300 and report the implied
#    period 2 pi / arccos(phi_1 / (2 sqrt(-phi_2))) (notes (3.15)). How variable is it over 200 runs?
# 3. **GARCH(1,1) preview of ST-8.** Simulate r_t = s_t z_t, s_t^2 = 0.05 + 0.1 r_{t-1}^2 + 0.85 s_{t-1}^2. Show the
#    sample ACF of r is inside the bands (white noise!) while the ACF of r^2 is not: the series is uncorrelated but
#    not independent, and no ARMA model of r can see the volatility clustering.
# 4. Change Part 9's n from 50 to 400 (exercise B7) and difference both series before regressing.

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T0:.1f} s.")
