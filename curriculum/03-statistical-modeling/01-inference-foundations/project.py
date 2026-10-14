# %% [markdown]
# # ST-1 Project: Sampling Distributions, Likelihood, and Information
#
# An estimator is a random variable. This script makes that literal: it draws
# 10^4 datasets, applies each estimator to each one, and looks at the resulting
# distribution. Parts: (1) a sampling-distribution simulator for Normal,
# Exponential and Uniform(0, theta), with bias/variance/MSE checked against the
# formulas of notes §3; (2) the bias-variance trade-off for variance estimators;
# (3) the Normal likelihood surface with its Fisher-information ellipse;
# (4) the Gamma MLE from the digamma equation; (5) Cramér-Rao ratios vs n;
# (6) asymptotic normality (QQ plots) and the Bernoulli exact-vs-normal check;
# (7) Wald intervals and the delta method; (8) Rao-Blackwell; (9) exponential
# family identities A' = E[T], A'' = Var[T]; (10) misspecification and the
# KL projection. Every number quoted in the notes' worked examples and in the
# solutions is asserted or printed here.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats, optimize, special, integrate
import sympy as sp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
R = 10_000  # number of simulated datasets


def mc_check(est, truth, bias_th, var_th, z=4.0, label=""):
    """Monte Carlo bias, variance and MSE of the draws `est` (one estimate per dataset),
    asserted equal to theory within z Monte Carlo standard errors. MSE = bias^2 + var, (3.2)."""
    err = est - truth
    bias, var, mse = err.mean(), est.var(), np.mean(err**2)
    se_b = est.std() / np.sqrt(est.size)
    se_v = np.std((est - est.mean())**2) / np.sqrt(est.size)
    se_m = np.std(err**2) / np.sqrt(est.size)
    mse_th = bias_th**2 + var_th
    assert abs(bias - bias_th) < z * se_b, (label, "bias", bias, bias_th)
    assert abs(var - var_th) < z * se_v, (label, "var", var, var_th)
    assert abs(mse - mse_th) < z * se_m, (label, "mse", mse, mse_th)
    return bias, var, mse, mse_th


# %% [markdown]
# ## Part 1: The sampling-distribution simulator
#
# Normal(mu=1, sigma^2=4): estimators of sigma^2 with divisors n-1, n, n+1 (§3.3, Exercise B1).
# Exponential(rate lambda=2): MLE = MoM = n/S and the unbiased (n-1)/S, S = sum X_i (Exercise B3).
# Uniform(0, theta=3): MoM 2*Xbar, MLE max X_i, unbiased (n+1)/n * max (§3.5, Example 4.4).

# %%
MU, SIG2, LAM, THETA = 1.0, 4.0, 2.0, 3.0
TRUTH = {"Normal": SIG2, "Exponential": LAM, "Uniform": THETA}


def draw(family, n, size=R):
    if family == "Normal":
        return rng.normal(MU, np.sqrt(SIG2), (size, n))
    if family == "Exponential":
        return rng.exponential(1 / LAM, (size, n))
    return rng.uniform(0, THETA, (size, n))


def estimators(family, X):
    """Each estimator applied to every row (dataset) of X."""
    n = X.shape[1]
    if family == "Normal":
        ss = ((X - X.mean(1, keepdims=True))**2).sum(1)
        return {"SS/(n-1) = S^2": ss / (n - 1), "SS/n = MLE": ss / n, "SS/(n+1)": ss / (n + 1)}
    if family == "Exponential":
        S = X.sum(1)
        return {"n/S = MLE = MoM": n / S, "(n-1)/S unbiased": (n - 1) / S}
    M = X.max(1)
    return {"2 Xbar = MoM": 2 * X.mean(1), "max = MLE": M, "(n+1)/n max": (n + 1) / n * M}


def theory(family, n):
    """(bias, variance) of each estimator: §3.3 (Normal), Exercise B3 (Exponential), Example 4.4 (Uniform)."""
    if family == "Normal":
        s4 = SIG2**2
        return [(0.0, 2 * s4 / (n - 1)), (-SIG2 / n, 2 * (n - 1) * s4 / n**2),
                (-2 * SIG2 / (n + 1), 2 * (n - 1) * s4 / (n + 1)**2)]
    if family == "Exponential":
        return [(LAM / (n - 1), n**2 * LAM**2 / ((n - 1)**2 * (n - 2))), (0.0, LAM**2 / (n - 2))]
    t2 = THETA**2
    return [(0.0, t2 / (3 * n)), (-THETA / (n + 1), n * t2 / ((n + 1)**2 * (n + 2))), (0.0, t2 / (n * (n + 2)))]


fig, axes = plt.subplots(3, 3, figsize=(13, 10))
print(f"{'family':12s}{'n':>4s}  {'estimator':20s}{'bias MC':>9s}{'theory':>9s}{'var MC':>9s}{'theory':>9s}{'MSE MC':>9s}{'theory':>9s}")
for i, fam in enumerate(TRUTH):
    for j, n in enumerate([5, 20, 100]):
        ests = estimators(fam, draw(fam, n))
        for (name, e), (b_th, v_th) in zip(ests.items(), theory(fam, n)):
            b, v, m, m_th = mc_check(e, TRUTH[fam], b_th, v_th, label=f"{fam} n={n} {name}")
            print(f"{fam:12s}{n:4d}  {name:20s}{b:9.4f}{b_th:9.4f}{v:9.4f}{v_th:9.4f}{m:9.4f}{m_th:9.4f}")
            lo, hi = np.quantile(e, [0.002, 0.995])
            axes[i, j].hist(e, bins=np.linspace(lo, hi, 80), histtype="step", density=True, label=name)
        axes[i, j].axvline(TRUTH[fam], color="k", ls="--", lw=1)
        axes[i, j].set_title(f"{fam}, n = {n}"); axes[i, j].legend(fontsize=7)
fig.suptitle("Sampling distributions of estimators (dashed: true value)")
fig.tight_layout(); fig.savefig(FIG / "sampling_distributions.png", dpi=110); plt.close(fig)
print("Part 1: all bias/variance/MSE values agree with theory within 4 Monte Carlo SEs.")

# %% [markdown]
# ## Part 2: Bias-variance trade-off for the variance estimators (Example 4.1, Exercise B1)
#
# With sigma^2 = 1: MSE(S^2) = 2/(n-1), MSE(MLE) = (2n-1)/n^2, MSE(SS/(n+1)) = 2/(n+1).

# %%
ns = np.arange(2, 51)
mse_S2, mse_mle, mse_np1 = 2 / (ns - 1), (2 * ns - 1) / ns**2, 2 / (ns + 1)
assert np.all(mse_mle < mse_S2) and np.all(mse_np1 <= mse_mle)
n = 10
ex41 = dict(S2_mse=2 / 9, mle_bias=-1 / 10, mle_var=18 / 100, mle_mse=19 / 100, np1_mse=2 / 11)
assert np.isclose(ex41["mle_bias"]**2 + ex41["mle_var"], ex41["mle_mse"])
X = rng.normal(0, 1, (R, n)); ss = ((X - X.mean(1, keepdims=True))**2).sum(1)
mc_check(ss / (n - 1), 1.0, 0.0, 2 / 9, label="Ex4.1 S2")
mc_check(ss / n, 1.0, -0.1, 0.18, label="Ex4.1 MLE")
mc_check(ss / (n + 1), 1.0, -2 / 11, 18 / 121, label="Ex4.1 n+1")
print(f"Example 4.1 (n=10, sigma^2=1): MSE(S^2)={2/9:.4f}, MLE bias={-0.1}, var={0.18}, "
      f"MSE={0.19:.4f}; SS/(n+1) MSE={2/11:.4f} (B1)")
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(ns, mse_S2, label="S^2 (unbiased): 2/(n-1)")
ax.plot(ns, mse_mle, label="MLE SS/n: (2n-1)/n^2")
ax.plot(ns, mse_np1, label="SS/(n+1): 2/(n+1)")
ax.plot(ns, 1 / ns**2, ":", label="MLE bias^2 = 1/n^2")
ax.set_yscale("log"); ax.set_xlabel("n"); ax.set_ylabel("MSE / sigma^4"); ax.legend(); ax.grid(alpha=.3)
ax.set_title("Biased estimators win on MSE"); fig.tight_layout(); fig.savefig(FIG / "variance_estimators_mse.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 3: Normal likelihood surface and the Fisher-information ellipse (Example 4.3)
#
# Data: eight measurements of a nominal 10-ohm resistor. In (mu, v = sigma^2) the per-sample
# Fisher matrix is diag(1/v, 1/(2 v^2)), (3.14). The surface is drawn in (mu, sigma), where the
# observed information at the MLE is diag(n/sigma^2, 2n/sigma^2).

# %%
xr = np.array([9.8, 10.4, 10.1, 9.6, 10.3, 10.0, 9.7, 10.1]); n = xr.size
mu_hat, v_hat = xr.mean(), xr.var()
assert np.isclose(mu_hat, 10.0) and np.isclose(v_hat, 0.07) and np.isclose(xr.var(ddof=1), 0.08)
In = np.diag([n / v_hat, n / (2 * v_hat**2)])
In_inv = np.linalg.inv(In)
assert np.allclose(In, np.diag([114.2857, 816.3265]), atol=1e-3)
assert np.allclose(In_inv, np.diag([0.00875, 0.001225]))
print(f"Example 4.3: mu_hat={mu_hat:.4f}, sigma2_hat={v_hat:.4f}, S^2={xr.var(ddof=1):.4f}; "
      f"n*I = diag({In[0,0]:.2f}, {In[1,1]:.2f}); se(mu)={np.sqrt(In_inv[0,0]):.4f}, se(sigma2)={np.sqrt(In_inv[1,1]):.4f}")


def loglik_mv(mu, v, x=xr):
    """Normal log-likelihood (3.10) in (mu, v = sigma^2)."""
    return -0.5 * x.size * np.log(2 * np.pi * v) - ((x - mu)**2).sum() / (2 * v)


h = 1e-4  # observed information by central differences equals n I_1 at the MLE for this model
H = np.array([[loglik_mv(mu_hat + h, v_hat) - 2 * loglik_mv(mu_hat, v_hat) + loglik_mv(mu_hat - h, v_hat), 0],
              [0, loglik_mv(mu_hat, v_hat + h) - 2 * loglik_mv(mu_hat, v_hat) + loglik_mv(mu_hat, v_hat - h)]]) / h**2
H[0, 1] = H[1, 0] = (loglik_mv(mu_hat + h, v_hat + h) - loglik_mv(mu_hat + h, v_hat - h)
                     - loglik_mv(mu_hat - h, v_hat + h) + loglik_mv(mu_hat - h, v_hat - h)) / (4 * h**2)
assert np.allclose(-H, In, rtol=1e-4, atol=1e-3)

s_hat = np.sqrt(v_hat)
M, S = np.meshgrid(np.linspace(mu_hat - 0.35, mu_hat + 0.35, 300), np.linspace(0.45 * s_hat, 2.4 * s_hat, 300))
L = -n * np.log(S) - (((xr[:, None, None] - M)**2).sum(0)) / (2 * S**2)
L -= L.max()
J = np.diag([n / v_hat, 2 * n / v_hat])  # observed info in (mu, sigma)
Q = -0.5 * (J[0, 0] * (M - mu_hat)**2 + J[1, 1] * (S - s_hat)**2)
c95 = 0.5 * stats.chi2.ppf(0.95, 2)
fig, ax = plt.subplots(figsize=(6.5, 5))
cs = ax.contourf(M, S, L, levels=np.linspace(-12, 0, 25), cmap="viridis")
ax.contour(M, S, L, levels=[-c95, -0.5], colors="w")
ax.contour(M, S, Q, levels=[-c95, -0.5], colors="r", linestyles="--")
ax.plot(mu_hat, s_hat, "r*", ms=12)
ax.set_xlabel("mu"); ax.set_ylabel("sigma"); fig.colorbar(cs, label="log L - max")
ax.set_title("White: exact contours. Red dashed: Fisher ellipse"); fig.tight_layout()
fig.savefig(FIG / "normal_likelihood_surface.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 4: Gamma(alpha, beta) by MoM and MLE (Example 4.5)
#
# Shape alpha, rate beta. MoM: alpha = xbar^2/s2, beta = xbar/s2 (s2 with divisor n). MLE: beta = alpha/xbar
# and log(alpha) - digamma(alpha) = log(xbar) - mean(log x), equation (3.13).

# %%
xg = np.array([0.8, 1.9, 2.4, 0.5, 3.1, 1.2, 2.0, 1.6, 0.9, 2.6])
xbar, s2 = xg.mean(), xg.var()
a_mom, b_mom = xbar**2 / s2, xbar / s2
assert np.isclose(xbar, 1.7) and np.isclose(s2, 0.654)
s_stat = np.log(xbar) - np.log(xg).mean()
a_mle = optimize.brentq(lambda a: np.log(a) - special.digamma(a) - s_stat, 1e-3, 1e3)
b_mle = a_mle / xbar
nll = lambda p: -stats.gamma.logpdf(xg, np.exp(p[0]), scale=np.exp(-p[1])).sum()
p_opt = optimize.minimize(nll, [np.log(a_mom), np.log(b_mom)], method="Nelder-Mead",
                          options=dict(xatol=1e-10, fatol=1e-12, maxiter=5000)).x
a_fit, _, scale_fit = stats.gamma.fit(xg, floc=0)
assert np.allclose([a_mle, b_mle], np.exp(p_opt), rtol=1e-4)
assert np.allclose([a_mle, b_mle], [a_fit, 1 / scale_fit], rtol=1e-3)
print(f"Example 4.5: xbar={xbar:.3f}, s2={s2:.3f}, mean log x={np.log(xg).mean():.4f}, s={s_stat:.4f}")
print(f"  MoM alpha={a_mom:.4f}, beta={b_mom:.4f};  MLE alpha={a_mle:.4f}, beta={b_mle:.4f} "
      f"(minimize: {np.exp(p_opt).round(4)}, scipy fit: {a_fit:.4f}, {1/scale_fit:.4f})")


def gamma_shape_mle(s, iters=30):
    """Vectorized Newton solve of log a - digamma(a) = s, started at Minka's approximation."""
    a = (3 - s + np.sqrt((s - 3)**2 + 24 * s)) / (12 * s)
    for _ in range(iters):
        a = np.maximum(a - (np.log(a) - special.digamma(a) - s) / (1 / a - special.polygamma(1, a)), 1e-8)
    return a


alpha0, n_g = 2.0, 50
Xg = rng.gamma(alpha0, 1.0, (4000, n_g))
a_mom_sim = Xg.mean(1)**2 / Xg.var(1)
a_mle_sim = gamma_shape_mle(np.log(Xg.mean(1)) - np.log(Xg).mean(1))
asym_mle = alpha0 / (alpha0 * special.polygamma(1, alpha0) - 1)   # [I^-1]_{alpha alpha}
asym_mom = 2 * alpha0 * (alpha0 + 1)                                # delta method on (m1, m2)
print(f"  Simulation alpha=2, n=50: n Var(MoM)={n_g*a_mom_sim.var():.2f} (asymptotic {asym_mom:.2f}), "
      f"n Var(MLE)={n_g*a_mle_sim.var():.2f} (asymptotic 1/I = {asym_mle:.2f}); "
      f"asymptotic efficiency of MoM = {asym_mle/asym_mom:.3f}")
assert n_g * a_mle_sim.var() < n_g * a_mom_sim.var()
assert abs(n_g * a_mle_sim.var() / asym_mle - 1) < 0.25 and abs(n_g * a_mom_sim.var() / asym_mom - 1) < 0.25

# %% [markdown]
# ## Part 5: Cramér-Rao ratios Var(MLE) / (1/(n I_1)) vs n (§3.7)
#
# Bernoulli and Poisson: exactly 1. Exponential rate: n^3/((n-1)^2 (n-2)) > 1, tending to 1.
# Normal variance MLE: (n-1)/n < 1 -- allowed, because that MLE is biased (biased bound (3.20)).

# %%
ns5, R5 = np.array([5, 10, 20, 50, 100, 200]), 20_000
ratio_th = {"Bernoulli p=0.3": lambda n: np.ones_like(n, float), "Poisson lam=3": lambda n: np.ones_like(n, float),
            "Exponential lam=2": lambda n: n**3 / ((n - 1.)**2 * (n - 2)), "Normal sigma^2 (MLE)": lambda n: (n - 1) / n}
sims = {"Bernoulli p=0.3": (lambda n: rng.binomial(n, 0.3, R5) / n, lambda n: 0.21 / n),
        "Poisson lam=3": (lambda n: rng.poisson(3 * n, R5) / n, lambda n: 3 / n),
        "Exponential lam=2": (lambda n: n / rng.gamma(n, 0.5, R5), lambda n: 4 / n),
        "Normal sigma^2 (MLE)": (lambda n: rng.chisquare(n - 1, R5) / n, lambda n: 2 / n)}
fig, ax = plt.subplots(figsize=(7, 4))
for k, (name, (sim, crb)) in enumerate(sims.items()):
    rat, se = [], []
    for n in ns5:
        e = sim(n)
        rat.append(e.var() / crb(n)); se.append(np.std((e - e.mean())**2) / np.sqrt(R5) / crb(n))
    rat, se = np.array(rat), np.array(se)
    assert np.all(np.abs(rat - ratio_th[name](ns5)) < 5 * se), (name, rat)
    ax.plot(ns5, rat, "o", color=f"C{k}", label=name)
    ax.plot(ns5, ratio_th[name](ns5), "-", color=f"C{k}", alpha=.5)
    print(f"Part 5 {name:22s} ratio at n=5,10,...,200:", np.round(rat, 3))
ax.axhline(1, color="k", lw=.8); ax.set_xscale("log"); ax.set_xlabel("n"); ax.set_ylabel("Var(MLE) * n I_1")
ax.legend(fontsize=8); ax.set_title("Cramér-Rao ratio (lines: exact theory)"); fig.tight_layout()
fig.savefig(FIG / "cramer_rao_ratios.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 6: Asymptotic normality (§3.8) and the Bernoulli exact check (Example 4.2)
#
# Z = sqrt(n)(lam_hat - lam) sqrt(I_1) for the Exponential rate MLE, against N(0,1) quantiles.

# %%
probs = np.linspace(0.01, 0.99, 99); zq = stats.norm.ppf(probs)
fig, axs = plt.subplots(1, 4, figsize=(14, 3.6), sharey=True)
gaps = []
for ax, n in zip(axs, [5, 20, 100, 500]):
    Z = np.sqrt(n) * (n / rng.gamma(n, 1.0, R) - 1.0)  # lambda = 1, I_1 = 1
    q = np.quantile(Z, probs); gaps.append(np.abs(q - zq)[(probs >= .05) & (probs <= .95)].max())
    ax.plot(zq, q, ".", ms=3); ax.plot(zq, zq, "k-", lw=.8); ax.set_title(f"n = {n}, gap (5-95%) {gaps[-1]:.2f}")
    ax.set_xlabel("N(0,1) quantile")
axs[0].set_ylabel("quantile of sqrt(n)(lam_hat - lam)/lam")
fig.tight_layout(); fig.savefig(FIG / "qq_exponential_mle.png", dpi=110); plt.close(fig)
assert all(np.diff(gaps) < 0) and gaps[-1] < 0.2
print("Part 6: max QQ gap (5%-95% quantiles) for n = 5, 20, 100, 500:", np.round(gaps, 3))

n, p, T_obs = 50, 0.3, 15
p_hat = T_obs / n; I1 = 1 / (p_hat * (1 - p_hat)); se_p = np.sqrt(p_hat * (1 - p_hat) / n)
assert np.isclose(p_hat, 0.3) and np.isclose(I1, 4.7619, atol=1e-4) and np.isclose(se_p**2, 0.0042)
exact = stats.binom.cdf(10, n, p)
approx = stats.norm.cdf((0.2 - p) / np.sqrt(p * (1 - p) / n))
approx_cc = stats.norm.cdf((10.5 - n * p) / np.sqrt(n * p * (1 - p)))
print(f"Example 4.2: p_hat={p_hat}, I1(p_hat)={I1:.4f}, CRB = Var(p_hat) = {se_p**2:.4f}, se={se_p:.4f}, "
      f"95% Wald CI = ({p_hat-1.96*se_p:.3f}, {p_hat+1.96*se_p:.3f})")
print(f"  P(p_hat <= 0.2): exact {exact:.4f}, normal {approx:.4f} (z={(0.2-p)/se_p:.3f}), "
      f"continuity-corrected {approx_cc:.4f}")

# %% [markdown]
# ## Part 7: Poisson Wald interval and the delta method (Example 4.6)
#
# Geiger-counter counts in 20 one-second windows. g(lambda) = exp(-lambda) = P(an empty window).

# %%
counts = np.array([3, 5, 2, 4, 3, 1, 4, 3, 6, 2, 3, 4, 2, 5, 3, 2, 4, 3, 2, 3]); n = counts.size
lam_hat = counts.mean(); se_lam = np.sqrt(lam_hat / n); J_obs = counts.sum() / lam_hat**2
ci = (lam_hat - 1.96 * se_lam, lam_hat + 1.96 * se_lam)
g_hat, se_g = np.exp(-lam_hat), np.exp(-lam_hat) * se_lam
assert counts.sum() == 64 and np.isclose(lam_hat, 3.2) and np.isclose(se_lam, 0.4) and np.isclose(J_obs, 6.25)
assert np.allclose(ci, (2.416, 3.984))
print(f"Example 4.6: lam_hat={lam_hat}, se={se_lam:.3f}, observed info={J_obs:.2f} = n/lam_hat, "
      f"95% CI=({ci[0]:.3f}, {ci[1]:.3f}); exp(-lam_hat)={g_hat:.5f}, delta se={se_g:.5f}")
for n_, tol in [(20, 0.15), (200, 0.03)]:
    Lsim = rng.poisson(3.2 * n_, R) / n_   # sample means at lambda = 3.2
    sd_g = np.exp(-Lsim).std(); delta_sd = np.exp(-3.2) * np.sqrt(3.2 / n_)
    cover = np.mean(np.abs(Lsim - 3.2) <= 1.96 * np.sqrt(Lsim / n_))
    assert abs(sd_g / delta_sd - 1) < tol
    print(f"  simulation n={n_}: sd of exp(-Xbar) = {sd_g:.5f} vs delta {delta_sd:.5f} (ratio {sd_g/delta_sd:.3f}); "
          f"Wald CI coverage = {cover:.3f}")

# %% [markdown]
# ## Part 8: Rao-Blackwell (Example 4.7, Exercise B5)
#
# Poisson(3.2), n = 20, T = sum X_i. E[X_1 | T] = T/n. E[1{X_1 = 0} | T] = (1 - 1/n)^T.

# %%
lam, n = 3.2, 20
Xp = rng.poisson(lam, (R, n)); T = Xp.sum(1)
mc_check(Xp[:, 0].astype(float), lam, 0.0, lam, label="X1")
mc_check(T / n, lam, 0.0, lam / n, label="Xbar")
for t in range(55, 75, 4):  # conditional mean of X_1 given T = t is t/n
    sel = T == t
    assert abs(Xp[sel, 0].mean() - t / n) < 5 * np.sqrt(t / n * (1 - 1 / n) / sel.sum())
g = np.exp(-lam)
v_ind, v_rb, crb_g = g * (1 - g), np.exp(-2 * lam) * (np.exp(lam / n) - 1), lam * np.exp(-2 * lam) / n
mc_check((Xp[:, 0] == 0).astype(float), g, 0.0, v_ind, label="1{X1=0}")
mc_check((1 - 1 / n)**T, g, 0.0, v_rb, label="RB of 1{X1=0}")
e_mle = np.exp(-T / n)
assert v_rb > crb_g
print(f"Example 4.7: Var(X1)={Xp[:,0].var():.3f} (theory {lam}), Var(T/n)={(T/n).var():.4f} (theory {lam/n})")
print(f"B5: Var(1{{X1=0}})={v_ind:.5f}, Var((1-1/n)^T)={v_rb:.3e}, CRB={crb_g:.3e}, ratio={v_rb/crb_g:.4f}; "
      f"MSE of MLE exp(-Xbar) (MC) = {np.mean((e_mle-g)**2):.3e}, bias = {e_mle.mean()-g:.2e}")

# %% [markdown]
# ## Part 9: Exponential-family identities A'(eta) = E[T], A''(eta) = Var[T] (§3.10)

# %%
def A_poisson(eta, K=200):
    """Log-partition log sum_x exp(eta x)/x!, by direct summation (should equal e^eta)."""
    x = np.arange(K)
    return special.logsumexp(eta * x - special.gammaln(x + 1))


def A_bern(eta):
    return np.log1p(np.exp(eta))


for name, A, eta, ET, VT in [("Poisson", A_poisson, np.log(3.2), 3.2, 3.2),
                             ("Bernoulli", A_bern, np.log(0.3 / 0.7), 0.3, 0.21)]:
    d = 1e-4
    A1 = (A(eta + d) - A(eta - d)) / (2 * d)
    A2 = (A(eta + d) - 2 * A(eta) + A(eta - d)) / d**2
    assert np.isclose(A1, ET, rtol=1e-6) and np.isclose(A2, VT, rtol=1e-4)
    print(f"Part 9 {name}: A'(eta)={A1:.6f} (E T={ET}), A''(eta)={A2:.5f} (Var T={VT})")
assert np.isclose(A_poisson(np.log(3.2)), 3.2)
e_ = sp.symbols("eta", real=True)
assert sp.simplify(sp.diff(sp.log(1 + sp.exp(e_)), e_) - 1 / (1 + sp.exp(-e_))) == 0
assert sp.simplify(sp.diff(sp.log(1 + sp.exp(e_)), e_, 2) - sp.exp(e_) / (1 + sp.exp(e_))**2) == 0
# Geometric (Exercise A5): A(eta) = eta - log(1 - e^eta); A' = 1/p, A'' = (1-p)/p^2 with e^eta = 1 - p.
pp = sp.symbols("p", positive=True); Ag = e_ - sp.log(1 - sp.exp(e_))
assert sp.simplify(sp.diff(Ag, e_).subs(e_, sp.log(1 - pp)) - 1 / pp) == 0
assert sp.simplify(sp.diff(Ag, e_, 2).subs(e_, sp.log(1 - pp)) - (1 - pp) / pp**2) == 0

# %% [markdown]
# ## Part 10: Misspecification -- a Normal fitted to Exponential(1) data (§3.12)
#
# The Normal MLE (Xbar, mean squared deviation) converges to the KL projection
# argmin D(q || N(mu, s^2)) = (E X, Var X) = (1, 1), and D at the projection is 0.5 log(2 pi e) - 1.

# %%
def kl_exp_normal(par):
    """D(Exp(1) || N(mu, s^2)) by quadrature, (3.29) of F5."""
    mu, s = par[0], np.exp(par[1])
    f = lambda x: np.exp(-x) * (-x - stats.norm.logpdf(x, mu, s))
    return integrate.quad(f, 0, np.inf)[0]


res = optimize.minimize(kl_exp_normal, [0.0, 0.0], method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-12))
mu_star, s2_star, D_star = res.x[0], np.exp(2 * res.x[1]), res.fun
assert np.allclose([mu_star, s2_star], [1, 1], atol=1e-4) and np.isclose(D_star, 0.5 * np.log(2 * np.pi * np.e) - 1)
nlist = [10, 100, 1000, 10_000, 100_000, 1_000_000]; fits = []
for n in nlist:
    x = rng.exponential(1.0, n); fits.append((x.mean(), x.var()))
fits = np.array(fits)
assert np.allclose(fits[-1], [1, 1], atol=0.01)
print(f"Part 10: KL projection (mu*, sigma*^2) = ({mu_star:.5f}, {s2_star:.5f}), D = {D_star:.4f}; "
      f"MLE at n = 1e6: ({fits[-1,0]:.4f}, {fits[-1,1]:.4f})")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
xx = np.linspace(-2.5, 5, 400)
a1.plot(xx, np.where(xx >= 0, np.exp(-xx), 0), label="truth q = Exp(1)")
a1.plot(xx, stats.norm.pdf(xx, 1, 1), label="KL-closest Normal N(1, 1)"); a1.legend(); a1.set_title("Misspecified fit")
a2.semilogx(nlist, fits[:, 0], "o-", label="mu_hat"); a2.semilogx(nlist, fits[:, 1], "s-", label="sigma2_hat")
a2.axhline(1, color="k", lw=.8); a2.legend(); a2.set_xlabel("n"); a2.set_title("MLE -> KL projection")
fig.tight_layout(); fig.savefig(FIG / "misspecification.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Numbers used in the exercise solutions

# %%
xs = np.array([3, 1, 4, 2, 5])                    # A2: Geometric
pg = xs.size / xs.sum(); se_pg = np.sqrt(pg**2 * (1 - pg) / xs.size)
assert np.isclose(pg, 1 / 3)
xu = np.array([-1.2, 0.4, 2.1, -0.7, 1.5])        # A3: Uniform(-theta, theta)
m2 = np.mean(xu**2); th_mom = np.sqrt(3 * m2)
assert np.isclose(m2, 1.75)
print(f"A2: p_hat={pg:.4f}, mean MLE={xs.mean()}, P(X=1) MLE={pg:.4f}, Wald se={se_pg:.4f}. "
      f"A3: m2={m2}, MoM theta={th_mom:.4f}, MLE max|x|={np.abs(xu).max()}")
n = 10                                            # B2: Uniform estimators, MSE/theta^2
print(f"B2 (n=10): MSE 2Xbar={1/(3*n):.5f}, max={2/((n+1)*(n+2)):.5f}, (n+1)/n max={1/(n*(n+2)):.5f}, "
      f"c* max={1/(n+1)**2:.5f} with c*={(n+2)/(n+1):.4f}")
print(f"B3 (n=10): bias/lam={1/(n-1):.4f}, Var(MLE)/lam^2={n**2/((n-1)**2*(n-2)):.4f}, "
      f"Var(unbiased)/lam^2={1/(n-2):.4f}, CRB/lam^2={1/n:.4f}")
eta_hat = np.log(0.3 / 0.7); se_eta = 1 / np.sqrt(50 * 0.3 * 0.7)
print(f"B4: eta_hat={eta_hat:.4f}, se={se_eta:.4f}")
w = 50 / 52; roots = np.roots([8.08, -8.08, 1.0])  # B6: crossover of MSE(posterior mean) and MSE(p_hat)
for r_ in roots:
    assert np.isclose((1 - w**2) * r_ * (1 - r_) / 50, (1 - w)**2 * (0.5 - r_)**2)
p0 = 0.3; mse_pm = w**2 * p0 * (1 - p0) / 50 + (1 - w)**2 * (0.5 - p0)**2
print(f"B6 / §3.11: posterior mean (T+1)/(n+2)={(15+1)/52:.4f}; at p=0.3 MSE={mse_pm:.6f} vs {0.0042}; "
      f"crossover p = {np.sort(roots).round(4)}")
mu_s, sg = sp.symbols("mu sigma", positive=True)  # A1: MSE of (X1+X2)/3 vs (X1+X2)/2
assert sp.simplify(mu_s**2 / 9 + 2 * sg**2 / 9 - sg**2 / 2 - (mu_s**2 - sp.Rational(5, 2) * sg**2) / 9) == 0
xx_, th_ = sp.symbols("x theta", positive=True)   # A4: Rayleigh, both forms of I equal 1/theta^2
logp = sp.log(xx_ / th_) - xx_**2 / (2 * th_)
dens = sp.exp(logp)
I_a = sp.integrate(sp.diff(logp, th_)**2 * dens, (xx_, 0, sp.oo))
I_b = -sp.integrate(sp.diff(logp, th_, 2) * dens, (xx_, 0, sp.oo))
assert sp.simplify(I_a - 1 / th_**2) == 0 and sp.simplify(I_b - 1 / th_**2) == 0
pairs = rng.normal(rng.normal(0, 5, (20_000, 1)), 1.0, (20_000, 2))  # C1: Neyman-Scott, sigma^2 = 1
ns_mle = ((pairs - pairs.mean(1, keepdims=True))**2).sum() / pairs.size
assert abs(ns_mle - 0.5) < 0.02
print(f"C1: Neyman-Scott MLE of sigma^2 with 20000 pairs = {ns_mle:.4f} (-> sigma^2/2 = 0.5)")

# %% [markdown]
# ## Your turn
#
# 1. (Exercise B8) Mixture weight. For w N(m1, 1) + (1-w) N(m2, 1), show numerically that for fixed
#    (m1, m2) the log-likelihood in w is concave, but with the means free there are two maxima
#    (label switching) and the likelihood at their midpoint is lower. Reference cell below.
# 2. (Exercise C2) Cauchy location. Compute I_1 = 1/2 by quadrature even though Var X does not exist.
#    Simulate n = 100: n Var(MLE) should be near 2, n Var(median) near pi^2/4, and the sample mean's
#    spread should not shrink at all. Reference cell below.
# 3. (ST-2 preview) Bootstrap. For the Gamma data of Part 4, resample xg with replacement 2000 times,
#    recompute alpha_MLE each time with `gamma_shape_mle`, and compare the bootstrap SD with the
#    Wald SE sqrt([I^-1]_{alpha alpha} / n) evaluated at the MLE. With n = 10, expect disagreement. Why?
# 4. In Part 5, add the unbiased Exponential estimator (n-1)/S. Its ratio is n/(n-2): never 1.

# %% [markdown]
# ## Reference solutions for B8 and C2 (write your own first)

# %%
z = rng.random(400) < 0.3
xm = np.where(z, rng.normal(-1, 1, 400), rng.normal(2, 1, 400))
mix_nll = lambda q: -np.log(q[0] * stats.norm.pdf(xm, q[1]) + (1 - q[0]) * stats.norm.pdf(xm, q[2])).sum()
b = [(1e-3, 1 - 1e-3), (-6, 6), (-6, 6)]
fA = optimize.minimize(mix_nll, [0.5, -2, 3], bounds=b); fB = optimize.minimize(mix_nll, [0.5, 3, -2], bounds=b)
assert np.isclose(fA.fun, fB.fun, rtol=1e-6) and np.allclose(fA.x, [1 - fB.x[0], fB.x[2], fB.x[1]], atol=1e-3)
mid = (fA.x + fB.x) / 2
assert mix_nll(mid) > fA.fun + 50                  # midpoint far worse: not concave
ww = np.linspace(0.01, 0.99, 197)
lw = np.array([-mix_nll([w_, fA.x[1], fA.x[2]]) for w_ in ww])
assert np.all(np.diff(lw, 2) < 0)                  # concave in w with the means fixed
print(f"B8: maxima {fA.x.round(3)} and {fB.x.round(3)}, -loglik {fA.fun:.2f}; at midpoint {mix_nll(mid):.2f}")

I_cauchy = integrate.quad(lambda x: (2 * x / (1 + x**2))**2 * stats.cauchy.pdf(x), -np.inf, np.inf)[0]
assert np.isclose(I_cauchy, 0.5)
Rc, n = 2000, 100
Xc = rng.standard_cauchy((Rc, n)); med = np.median(Xc, 1)
grid = med[:, None] + np.linspace(-1, 1, 401)[None, :]
ll = -np.log1p((Xc[:, None, :] - grid[:, :, None])**2).sum(2)
th = grid[np.arange(Rc), ll.argmax(1)]
for _ in range(5):                                  # Newton polish of the grid maximizer
    u = Xc - th[:, None]
    th = th + (2 * u / (1 + u**2)).sum(1) / (2 * ((1 - u**2) / (1 + u**2)**2).sum(1))
iqr = lambda v: np.subtract(*np.quantile(v, [0.75, 0.25]))
print(f"C2: I_1 = {I_cauchy:.4f}; n Var(MLE) = {n*th.var():.3f} (1/I = 2), n Var(median) = {n*med.var():.3f} "
      f"(pi^2/4 = {np.pi**2/4:.3f}); IQR of Xbar = {iqr(Xc.mean(1)):.2f} (one observation: 2)")
assert abs(n * th.var() - 2) < 0.4 and n * med.var() > n * th.var()

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
