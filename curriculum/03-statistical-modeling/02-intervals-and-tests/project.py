# %% [markdown]
# # ST-2 Project: Intervals, Tests, the Bootstrap, and How They Fail
#
# Every procedure in ST-2 makes a promise about long-run frequencies ("covers
# 95% of the time", "rejects a true null 5% of the time"). This script checks
# those promises by brute force and shows where they break: coverage of z, t,
# Wald and Wilson intervals; the chi^2 and t sampling distributions; a small
# bootstrap library and its failure for the maximum; a permutation test; Wilks'
# theorem and the trinity; chi^2 and ANOVA via projections; uniform p-values;
# power and sample size; the winner's curse; p-hacking; Bonferroni, Holm, BH.
# Every number quoted in notes.md, exercises.md or solutions.md is asserted or printed here.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import itertools
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy import integrate, stats

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
Z = stats.norm.ppf(0.975)  # z_{alpha/2} for alpha = 0.05, = 1.95996

# %% [markdown]
# ## Part 1: Sampling distributions via projection (notes §3.1, Example 4.1)
#
# Rotate the data into an orthonormal basis whose first vector is 1/sqrt(n).
# The first coordinate carries the mean; the other n-1 carry (n-1)S^2.

# %%
x3 = np.array([2.0, 5.0, 8.0])
Q = np.column_stack([np.ones(3) / np.sqrt(3), [1 / np.sqrt(2), -1 / np.sqrt(2), 0],
                     np.array([1, 1, -2]) / np.sqrt(6)])
assert np.allclose(Q.T @ Q, np.eye(3))
w = Q.T @ x3
assert np.isclose(w[0], np.sqrt(3) * 5) and np.isclose(w[1] ** 2, 4.5) and np.isclose(w[2] ** 2, 13.5)
assert np.isclose(w[1] ** 2 + w[2] ** 2, ((x3 - x3.mean()) ** 2).sum())  # = 18 = (n-1) S^2
print(f"Ex 4.1: W = {np.round(w, 4)}, W2^2 + W3^2 = {w[1]**2 + w[2]**2:.1f} = (n-1)S^2")

n, R, mu, sig = 5, 10_000, 10.0, 2.0
X = rng.normal(mu, sig, (R, n))
xbar, s2 = X.mean(1), X.var(1, ddof=1)
V = (n - 1) * s2 / sig ** 2                       # Theorem 3.1(b): chi^2_{n-1}
T = (xbar - mu) / np.sqrt(s2 / n)                 # (3.4): t_{n-1}
assert stats.kstest(V, stats.chi2(n - 1).cdf).pvalue > 0.01
assert stats.kstest(T, stats.t(n - 1).cdf).pvalue > 0.01
assert stats.kstest(T, stats.norm.cdf).pvalue < 1e-3   # the t really is not Normal at n = 5
assert abs(np.corrcoef(xbar, s2)[0, 1]) < 0.03          # Theorem 3.1(c): independent
E = rng.exponential(1.0, (R, n))
r_exp = np.corrcoef(E.mean(1), E.var(1, ddof=1))[0, 1]
assert r_exp > 0.3                                      # independence is special to the Normal
print(f"Part 1: KS tests pass; corr(xbar, S^2) Normal = {np.corrcoef(xbar, s2)[0,1]:+.3f}, Exponential = {r_exp:.2f}")

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
g = np.linspace(0.01, 18, 300)
ax[0].hist(V, 80, density=True, alpha=.5); ax[0].plot(g, stats.chi2(n - 1).pdf(g)); ax[0].set_title("(n-1)S²/σ² vs χ²₄")
g = np.linspace(-6, 6, 300)
ax[1].hist(T[np.abs(T) < 6], 120, density=True, alpha=.5)
ax[1].plot(g, stats.t(n - 1).pdf(g), label="t₄"); ax[1].plot(g, stats.norm.pdf(g), "--", label="N(0,1)")
ax[1].legend(); ax[1].set_title("t statistic, n = 5")
fig.tight_layout(); fig.savefig(FIG / "sampling_distributions.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 2: Pivot intervals and their coverage (§3.2, Examples 4.2 and 4.3)


# %%
def t_interval(x, alpha=0.05):
    """(3.6): xbar -/+ t_{n-1, alpha/2} S / sqrt(n)."""
    x = np.asarray(x, float); n = len(x)
    h = stats.t.ppf(1 - alpha / 2, n - 1) * x.std(ddof=1) / np.sqrt(n)
    return x.mean() - h, x.mean() + h


def wald_prop(k, n, z=Z):
    """(3.8): phat -/+ z sqrt(phat(1-phat)/n)."""
    ph = np.asarray(k) / n; h = z * np.sqrt(ph * (1 - ph) / n)
    return ph - h, ph + h


def wilson_prop(k, n, z=Z):
    """(3.9): roots of (phat - p)^2 = z^2 p(1-p)/n."""
    ph = np.asarray(k) / n
    c = (ph + z ** 2 / (2 * n)) / (1 + z ** 2 / n)
    h = z * np.sqrt(ph * (1 - ph) / n + z ** 2 / (4 * n ** 2)) / (1 + z ** 2 / n)
    return c - h, c + h


def exact_coverage(ci, n, p):
    """C(p) = sum_k Bin(k; n, p) 1{p in CI(k)}: coverage by enumeration of all outcomes."""
    k = np.arange(n + 1); lo, hi = ci(k, n)
    return stats.binom.pmf(k, n, p)[(lo <= p) & (p <= hi)].sum()


# SymPy: the Wilson endpoints are the roots of the score-test quadratic.
p_, ph_, z_, n_ = sp.symbols("p phat z n", positive=True)
roots = sp.solve(sp.Eq((ph_ - p_) ** 2, z_ ** 2 * p_ * (1 - p_) / n_), p_)
subs = {ph_: sp.Rational(1, 10), z_: sp.Float(Z), n_: 20}
assert np.allclose(sorted(float(r.subs(subs)) for r in roots), wilson_prop(2, 20))

g_data = np.array([9.79, 9.83, 9.76, 9.85, 9.81, 9.78, 9.88, 9.80, 9.74, 9.82])
lo, hi = t_interval(g_data)
s_g = g_data.std(ddof=1); se_g = s_g / np.sqrt(10)
assert np.isclose(g_data.mean(), 9.806) and np.isclose(((g_data - 9.806) ** 2).sum(), 0.01564)
assert abs(lo - 9.7762) < 1e-4 and abs(hi - 9.8358) < 1e-4, (lo, hi)
c_lo, c_hi = stats.chi2.ppf([0.025, 0.975], 9)
sig_ci = np.sqrt(0.01564 / np.array([c_hi, c_lo]))
print(f"Ex 4.2: S = {s_g:.5f}, se = {se_g:.5f}, t_9 = {stats.t.ppf(.975, 9):.4f}, CI = [{lo:.4f}, {hi:.4f}]; "
      f"chi2_9 quantiles {c_lo:.3f}, {c_hi:.3f}; sigma in [{sig_ci[0]:.4f}, {sig_ci[1]:.4f}]")

wl, wh = wald_prop(2, 20); sl, sh = wilson_prop(2, 20)
cov = {f: [exact_coverage(f, 20, p) for p in (0.10, 0.05)] for f in (wald_prop, wilson_prop)}
print(f"Ex 4.3: Wald [{wl:.4f}, {wh:.4f}], Wilson [{sl:.4f}, {sh:.4f}]; exact coverage at p=0.10, 0.05: "
      f"Wald {cov[wald_prop][0]:.4f}, {cov[wald_prop][1]:.4f}; Wilson {cov[wilson_prop][0]:.4f}, {cov[wilson_prop][1]:.4f}")
assert cov[wald_prop][1] < 0.70 and cov[wilson_prop][1] > 0.90

# Monte Carlo coverage, 10^4 replications per cell.
print("Coverage (10^4 reps)   n:  z(sigma known)  t     z-with-S")
for n in (5, 10, 20, 50):
    X = rng.normal(0, 1, (R, n)); m = X.mean(1); s = X.std(1, ddof=1)
    cz = np.mean(np.abs(m) <= Z / np.sqrt(n))
    ct = np.mean(np.abs(m) <= stats.t.ppf(.975, n - 1) * s / np.sqrt(n))
    cS = np.mean(np.abs(m) <= Z * s / np.sqrt(n))
    print(f"{'':22s}{n:3d}:  {cz:.3f}          {ct:.3f}  {cS:.3f}")
    assert abs(cz - .95) < .01 and abs(ct - .95) < .01
    if n == 5:
        assert cS < 0.90  # plugging S into the z interval undercovers at small n
k_sim = rng.binomial(20, 0.05, R); lo_s, hi_s = wald_prop(k_sim, 20)
mc_wald = np.mean((lo_s <= .05) & (.05 <= hi_s))
assert abs(mc_wald - cov[wald_prop][1]) < 0.015
pgrid = np.linspace(0.01, 0.99, 197)
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
for a, n in zip(ax, (20, 100)):
    for f, lab in ((wald_prop, "Wald"), (wilson_prop, "Wilson")):
        c = np.array([exact_coverage(f, n, p) for p in pgrid]); a.plot(pgrid, c, label=lab, lw=1)
        print(f"n = {n:3d}: mean exact coverage over p in [0.01, 0.99]: {lab} {c.mean():.3f}")
        if lab == "Wilson":
            assert abs(c.mean() - 0.95) < 0.01
    a.axhline(.95, color="k", ls=":"); a.set_title(f"exact coverage, n = {n}"); a.set_xlabel("p"); a.legend()
ax[0].set_ylim(.5, 1); fig.tight_layout(); fig.savefig(FIG / "proportion_coverage.png", dpi=120); plt.close(fig)
print(f"MC check: Wald coverage at n=20, p=0.05 is {mc_wald:.3f} (exact {cov[wald_prop][1]:.3f})")

# %% [markdown]
# ## Part 3: A bootstrap library (§3.4, Example 4.4)
#
# `stat` maps an array of shape (B, n, ...) to shape (B,), reducing axis 1.


# %%
def mean_stat(a): return a.mean(1)
def median_stat(a): return np.median(a, 1)


def corr_stat(a):
    """Pearson r of the pairs a[..., 0], a[..., 1] along axis 1."""
    x = a[..., 0] - a[..., 0].mean(1, keepdims=True); y = a[..., 1] - a[..., 1].mean(1, keepdims=True)
    return (x * y).sum(1) / np.sqrt((x * x).sum(1) * (y * y).sum(1))


def bootstrap_reps(data, stat, B=2000, rng=rng):
    """theta*_b = stat(X*_b), X*_b drawn n times with replacement from the data: draws from F-hat_n."""
    data = np.asarray(data, float)
    return stat(data[rng.integers(0, len(data), (B, len(data)))])


def bootstrap_se(data, stat, B=2000, rng=rng):
    """(3.13): se_boot = sd of the bootstrap replicates."""
    return bootstrap_reps(data, stat, B, rng).std(ddof=1)


def bootstrap_ci(data, stat, B=2000, alpha=0.05, kinds=("normal", "percentile", "basic"), rng=rng):
    """Normal (3.14), percentile (3.15) and basic/pivotal (3.16) intervals from one set of replicates."""
    th = stat(np.asarray(data, float)[None])[0]
    reps = bootstrap_reps(data, stat, B, rng)
    ql, qh = np.quantile(reps, [alpha / 2, 1 - alpha / 2]); h = stats.norm.ppf(1 - alpha / 2) * reps.std(ddof=1)
    out = {"normal": (th - h, th + h), "percentile": (ql, qh), "basic": (2 * th - qh, 2 * th - ql)}
    return {k: out[k] for k in kinds}


med_data = np.array([3.1, 4.7, 2.2, 8.9, 5.0, 3.8, 12.4, 4.1, 6.3, 2.9, 5.6, 4.4, 7.2, 3.5, 19.8])
xs = np.sort(med_data); n = 15
G = stats.binom.sf(7, n, np.arange(n + 1) / n)          # P(median* <= x_(k)) = P(Bin(n, k/n) >= 8)
wk = np.diff(G)                                         # P(median* = x_(k)), k = 1..15
se_exact = np.sqrt(wk @ xs ** 2 - (wk @ xs) ** 2)
reps_med = bootstrap_reps(med_data, median_stat, 10_000)
se_mc = reps_med.std(ddof=1)
assert np.median(med_data) == 4.7 and abs(se_mc - se_exact) / se_exact < 0.03
top = np.argsort(wk)[::-1][:3]
print(f"Ex 4.4: median 4.7; bootstrap SE (B=10^4) = {se_mc:.3f}, ideal (B=inf) = {se_exact:.3f}; "
      f"P(median*=x_(k)) for x = {xs[top]}: {np.round(wk[top], 3)}; s/sqrt(n) = {med_data.std(ddof=1)/np.sqrt(n):.3f}")
fig, ax = plt.subplots(figsize=(5, 3.5))
ax.hist(reps_med, bins=np.arange(2, 9, 0.05)); ax.set_title("bootstrap medians, Example 4.4 (B = 10⁴)")
fig.tight_layout(); fig.savefig(FIG / "bootstrap_median.png", dpi=120); plt.close(fig)

# Coverage of the three intervals: Exponential(1) data for mean/median, bivariate Normal rho=0.5 for r.
Rb, nb, Bb = 400, 30, 1000
truth = {"mean": 1.0, "median": np.log(2), "corr": 0.5}
hits = {(s, k): 0 for s in truth for k in ("normal", "percentile", "basic")}
L = np.linalg.cholesky([[1, .5], [.5, 1]])
for _ in range(Rb):
    e = rng.exponential(1.0, nb); xy = rng.standard_normal((nb, 2)) @ L.T
    for s, d, f in (("mean", e, mean_stat), ("median", e, median_stat), ("corr", xy, corr_stat)):
        for k, (a, b) in bootstrap_ci(d, f, Bb).items():
            hits[s, k] += a <= truth[s] <= b
print("Bootstrap coverage (n=30, 400 reps, B=1000):")
for s in truth:
    c = {k: hits[s, k] / Rb for k in ("normal", "percentile", "basic")}
    print(f"   {s:6s} " + "  ".join(f"{k} {v:.3f}" for k, v in c.items()))
    assert min(c["normal"], c["percentile"]) > 0.87 and c["basic"] > 0.78  # basic is the weakest here

# Failure: the maximum of Uniform(0, theta). P*(max* = max) = 1 - (1 - 1/n)^n -> 1 - 1/e.
nu = 50; u = rng.uniform(0, 1, nu); r_max = bootstrap_reps(u, lambda a: a.max(1), 10_000)
atom = np.mean(r_max == u.max()); atom_th = 1 - (1 - 1 / nu) ** nu
assert abs(atom - atom_th) < 0.015
cov_max = {"percentile": 0, "basic": 0, "exact pivot": 0}
for _ in range(Rb):
    u = rng.uniform(0, 1, nu); ci = bootstrap_ci(u, lambda a: a.max(1), Bb, kinds=("percentile", "basic"))
    for k in ("percentile", "basic"):
        cov_max[k] += ci[k][0] <= 1.0 <= ci[k][1]
    cov_max["exact pivot"] += u.max() <= 1.0 <= u.max() / 0.05 ** (1 / nu)  # X_(n)/theta ~ Beta(n, 1)
cov_max = {k: v / Rb for k, v in cov_max.items()}
print(f"Max of U(0,1), n=50: atom P*(max*=max) = {atom:.3f} (theory {atom_th:.3f}, limit {1-np.exp(-1):.3f}); coverage {cov_max}")
assert cov_max["percentile"] == 0 and cov_max["basic"] < 0.90 and cov_max["exact pivot"] > 0.92

# %% [markdown]
# ## Part 4: A permutation test by full enumeration (§3.4)

# %%
trt = np.array([5.1, 6.3, 7.2, 5.8, 6.9, 7.5]); ctl = np.array([4.2, 5.0, 5.5, 4.8, 6.1, 4.6])
pool = np.r_[trt, ctl]; d_obs = trt.mean() - ctl.mean()
diffs = np.array([pool[list(c)].mean() - np.delete(pool, c).mean() for c in itertools.combinations(range(12), 6)])
assert len(diffs) == 924
p_perm = np.mean(np.abs(diffs) >= abs(d_obs) - 1e-12)
p_t = stats.ttest_ind(trt, ctl).pvalue
print(f"Permutation test: diff {d_obs:.3f}, exact p = {p_perm:.4f} over 924 splits; pooled t p = {p_t:.4f}")
assert abs(p_perm - p_t) < 0.01

# %% [markdown]
# ## Part 5: Neyman–Pearson, power, sample size (§3.5, §3.9, Examples 4.5 and 4.8)

# %%
n, s, a = 16, 2.0, 0.05
c = stats.norm.ppf(1 - a) * s / np.sqrt(n); power_np = stats.norm.sf((c - 1) / (s / np.sqrt(n)))
xb = rng.normal(0, s / np.sqrt(n), 10 ** 5); xb1 = rng.normal(1, s / np.sqrt(n), 10 ** 5)
assert abs(c - 0.8224) < 1e-4 and abs(power_np - 0.6387) < 1e-4
assert abs(np.mean(xb > c) - .05) < .003 and abs(np.mean(xb1 > c) - power_np) < .005
print(f"Ex 4.5: reject if xbar > {c:.4f}; power at mu=1 is {power_np:.4f}")

zb = stats.norm.ppf(0.8); n_z = ((Z + zb) / 0.5) ** 2
pow_z = lambda n, d: stats.norm.sf(Z - d * np.sqrt(n)) + stats.norm.cdf(-Z - d * np.sqrt(n))


def pow_t(n, d, a=0.05):
    """Two-sided one-sample t-test power: noncentral t with delta = d sqrt(n)."""
    tc = stats.t.ppf(1 - a / 2, n - 1); nc = d * np.sqrt(n)
    return stats.nct.sf(tc, n - 1, nc) + stats.nct.cdf(-tc, n - 1, nc)


n_t = next(m for m in range(2, 200) if pow_t(m, 0.5) >= 0.8)
sim = rng.normal(0.5, 1, (20_000, n_t)); p_sim = stats.ttest_1samp(sim, 0, axis=1).pvalue
assert np.isclose(n_z, 31.39, atol=0.01) and pow_z(32, .5) >= .8 and pow_z(31, .5) < .8
assert abs(np.mean(p_sim < .05) - pow_t(n_t, .5)) < .01
print(f"Ex 4.8: n_z = {n_z:.2f} -> 32 (z power {pow_z(32, .5):.3f}); t power at 32 = {pow_t(32, .5):.3f}; "
      f"t needs n = {n_t} (power {pow_t(n_t, .5):.3f}, simulated {np.mean(p_sim < .05):.3f})")
fig, ax = plt.subplots(figsize=(6, 4)); dg = np.linspace(0, 1.2, 100)
for m in (10, 20, 50, 100):
    ax.plot(dg, pow_t(m, dg), label=f"n = {m}")
ax.axhline(.8, ls=":", c="k"); ax.set_xlabel("effect size δ/σ"); ax.set_ylabel("power"); ax.legend()
ax.set_title("two-sided one-sample t test, α = 0.05"); fig.tight_layout(); fig.savefig(FIG / "power_curves.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 6: Likelihood-ratio tests, Wilks, and the trinity (§3.6, Example 4.6)

# %%
n = 10; mu0 = 9.78; xb = g_data.mean()
t_obs = (xb - mu0) / se_g
lam = (np.mean((g_data - xb) ** 2) / np.mean((g_data - mu0) ** 2)) ** (n / 2)   # sup over sigma^2 done analytically
assert np.isclose(lam, (1 + t_obs ** 2 / (n - 1)) ** (-n / 2))                    # (4.2)
m2l = -2 * np.log(lam)
p_exact, p_wilks = 2 * stats.t.sf(abs(t_obs), n - 1), stats.chi2.sf(m2l, 1)
print(f"Ex 4.6: t = {t_obs:.4f}, lambda = {lam:.4f}, -2 log lambda = {m2l:.3f}; p exact t = {p_exact:.4f}, Wilks = {p_wilks:.4f}")

# Wilks: Poisson rate, H0: lambda = 2, n = 30. SymPy: the statistic is n(xbar - l0)^2 / l0 to second order.
xs_, l0_, n_s = sp.symbols("x lambda0 n", positive=True)
lrt_sym = 2 * n_s * (xs_ * sp.log(xs_ / l0_) - (xs_ - l0_))
assert sp.simplify(sp.series(lrt_sym, xs_, l0_, 3).removeO() - n_s * (xs_ - l0_) ** 2 / l0_) == 0


def poisson_trinity(xbar, n, l0):
    """LRT 2n[x log(x/l0) - (x - l0)], Wald n(x - l0)^2/x, score n(x - l0)^2/l0 (§3.6)."""
    xlogx = np.where(xbar > 0, xbar * np.log(np.where(xbar > 0, xbar, 1) / l0), 0.0)
    return 2 * n * (xlogx - (xbar - l0)), n * (xbar - l0) ** 2 / np.maximum(xbar, 1e-300), n * (xbar - l0) ** 2 / l0


n, l0 = 30, 2.0
xbar_null = rng.poisson(l0, (20_000, n)).mean(1)
lr, wa, sc = poisson_trinity(xbar_null, n, l0)
crit = stats.chi2.ppf(.95, 1)
rates = [np.mean(v > crit) for v in (lr, wa, sc)]
qs = np.quantile(lr, [.9, .95, .99]) / stats.chi2.ppf([.9, .95, .99], 1)
print(f"Wilks (Poisson, n=30): null rejection rates LRT/Wald/score = {np.round(rates, 4)}; quantile ratios {np.round(qs, 3)}")
assert all(abs(r - .05) < .01 for r in rates) and np.all(np.abs(qs - 1) < .1)
one = rng.poisson(2.6, n).mean()
print(f"Trinity on one sample (xbar = {one:.3f}, H0 lambda = 2): LRT, Wald, score = {np.round(poisson_trinity(one, n, l0), 3)}")
fig, ax = plt.subplots(figsize=(4.5, 4.5)); pp = (np.arange(1, 2001) - .5) / 2000
ax.plot(stats.chi2.ppf(pp, 1), np.quantile(lr, pp), ".", ms=2); ax.plot([0, 12], [0, 12], "k:")
ax.set_xlabel("χ²₁ quantile"); ax.set_ylabel("simulated −2 log λ"); ax.set_title("Wilks QQ, Poisson n = 30")
fig.tight_layout(); fig.savefig(FIG / "wilks_qq.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: The chi^2 test and ANOVA as projections (§3.7, Example 4.7)

# %%
O = np.array([[30, 20], [18, 32]]); Ex = O.sum(1, keepdims=True) * O.sum(0) / O.sum()
X2 = ((O - Ex) ** 2 / Ex).sum(); Gs = 2 * (O * np.log(O / Ex)).sum()
ref = stats.chi2_contingency(O, correction=False)
assert np.allclose(Ex, [[24, 26], [24, 26]]) and np.isclose(X2, ref[0]) and np.isclose(X2, 5.769, atol=1e-3)
d = .6 - .36; se_d = np.sqrt(.6 * .4 / 50 + .36 * .64 / 50); zp = d / np.sqrt(.48 * .52 * (2 / 50))
assert np.isclose(zp ** 2, X2)
print(f"Ex 4.7: X^2 = {X2:.3f} (p {stats.chi2.sf(X2, 1):.4f}), G = {Gs:.3f} (p {stats.chi2.sf(Gs, 1):.4f}); "
      f"z = {zp:.4f}; diff 0.24, 95% CI [{d - Z*se_d:.3f}, {d + Z*se_d:.3f}]; odds ratio {30*32/(20*18):.3f}")


def anova_projection(groups):
    """F = [||(P1 - P0) y||^2/(k-1)] / [||(I - P1) y||^2/(N-k)]  (3.30), from explicit projection matrices."""
    y = np.concatenate(groups); N, k = len(y), len(groups)
    A = np.zeros((N, k)); A[np.arange(N), np.repeat(np.arange(k), [len(gg) for gg in groups])] = 1
    P1 = A @ np.linalg.solve(A.T @ A, A.T); P0 = np.full((N, N), 1 / N); I = np.eye(N)
    assert np.allclose((P1 - P0) @ (I - P1), 0) and np.allclose(P1 @ P0, P0)        # orthogonal pieces, V0 in V1
    assert np.isclose(np.trace(P1 - P0), k - 1) and np.isclose(np.trace(I - P1), N - k)  # dims = df
    ssb, ssw = y @ (P1 - P0) @ y, y @ (I - P1) @ y
    assert np.isclose(ssb + ssw, ((y - y.mean()) ** 2).sum())
    return (ssb / (k - 1)) / (ssw / (N - k))


grp = [rng.normal(m, 1, s) for m, s in zip((0, .5, 1), (6, 8, 7))]
F_proj = anova_projection(grp)
assert np.isclose(F_proj, stats.f_oneway(*grp).statistic)
two = [rng.normal(0, 1, 9), rng.normal(.7, 1, 11)]
assert np.isclose(anova_projection(two), stats.ttest_ind(*two).statistic ** 2)       # k = 2: F = t^2
Fnull = np.array([anova_projection([rng.normal(0, 1, s) for s in (6, 8, 7)]) for _ in range(1500)])
assert stats.kstest(Fnull, stats.f(2, 18).cdf).pvalue > 0.01
print(f"ANOVA: F by projection = {F_proj:.4f} = scipy f_oneway; two-group F = t^2; null F ~ F(2,18) (KS ok)")

# %% [markdown]
# ## Part 8: p-values, the winner's curse, p-hacking (§3.8, §3.9)

# %%
p_null = stats.ttest_1samp(rng.normal(0, 1, (R, 10)), 0, axis=1).pvalue
p_alt = stats.ttest_1samp(rng.normal(0.5, 1, (R, 10)), 0, axis=1).pvalue
assert stats.kstest(p_null, "uniform").pvalue > 0.01 and abs(np.mean(p_null < .05) - .05) < .006
fig, ax = plt.subplots(figsize=(6, 3.5)); bins = np.linspace(0, 1, 21)
ax.hist(p_null, bins, alpha=.6, label="H0 true"); ax.hist(p_alt, bins, alpha=.6, label="δ = 0.5σ, n = 10")
ax.set_xlabel("p-value"); ax.legend(); fig.tight_layout(); fig.savefig(FIG / "pvalue_histograms.png", dpi=120); plt.close(fig)

n, dlt = 20, 0.2
Xw = rng.normal(dlt, 1, (100_000, n)); tt = stats.ttest_1samp(Xw, 0, axis=1)
pub = (tt.pvalue < .05) & (tt.statistic > 0); est = Xw.mean(1)
se = 1 / np.sqrt(n); cz = Z * se; aa = (cz - dlt) / se
e_pub = dlt + se * stats.norm.pdf(aa) / stats.norm.sf(aa)                 # (3.38), z approximation
print(f"Winner's curse: power {np.mean(tt.pvalue < .05):.3f} (z approx {pow_z(n, dlt):.3f}); mean published estimate "
      f"{est[pub].mean():.3f} = {est[pub].mean()/dlt:.2f} x truth (formula {e_pub:.3f} = {e_pub/dlt:.2f} x); "
      f"wrong-sign share of significant {np.mean(tt.statistic[tt.pvalue < .05] < 0):.3f}")
assert 2.5 < est[pub].mean() / dlt < 3.3 and abs(e_pub / dlt - 2.77) < .01

for m, label in ((20, "20 outcomes"), (5, "5 subgroups")):
    a_, b_ = rng.normal(0, 1, (2, 5000, m, 20))
    fp = np.mean((stats.ttest_ind(a_, b_, axis=2).pvalue < .05).any(1))
    print(f"p-hacking, {label}: false-positive rate {fp:.3f} (1 - 0.95^{m} = {1 - .95 ** m:.3f})")
    assert abs(fp - (1 - .95 ** m)) < .02

# %% [markdown]
# ## Part 9: Multiple testing (§3.10, Example 4.9)


# %%
def bonferroni(p, a=0.05):
    """Reject H_i iff p_i <= a/m (Theorem 3.12)."""
    return np.asarray(p) <= a / len(p)


def holm(p, a=0.05):
    """Step-down: walk up the sorted p's, reject while p_(j) <= a/(m - j + 1), stop at the first failure."""
    p = np.asarray(p); m = len(p); o = np.argsort(p); ok = p[o] <= a / (m - np.arange(m))
    k = m if ok.all() else int(np.argmin(ok)); rej = np.zeros(m, bool); rej[o[:k]] = True
    return rej


def bh(p, q=0.05):
    """Step-up: k = max{j : p_(j) <= j q/m}; reject the k smallest (Theorem 3.13)."""
    p = np.asarray(p); m = len(p); o = np.argsort(p); ok = p[o] <= q * np.arange(1, m + 1) / m
    k = int(np.nonzero(ok)[0].max() + 1) if ok.any() else 0; rej = np.zeros(m, bool); rej[o[:k]] = True
    return rej


p49 = np.array([0.001, 0.004, 0.006, 0.012, 0.019, 0.033, 0.047, 0.21, 0.38, 0.73])
counts = [int(f(p49).sum()) for f in (bonferroni, holm, bh)]
assert counts == [2, 3, 5], counts
print(f"Ex 4.9: rejections Bonferroni/Holm/BH = {counts}")

m, m0, reps = 100, 80, 2000
res = {f.__name__: np.zeros((reps, 3)) for f in (bonferroni, holm, bh)}
for r in range(reps):
    zz = rng.standard_normal(m); zz[m0:] += 3.0; pv = stats.norm.sf(zz)
    for f in (bonferroni, holm, bh):
        rej = f(pv); v = rej[:m0].sum(); res[f.__name__][r] = v, rej.sum(), rej[m0:].sum()
    assert np.all(holm(pv) >= bonferroni(pv))
summ = {k: (np.mean(v[:, 0] > 0), np.mean(v[:, 0] / np.maximum(v[:, 1], 1)), np.mean(v[:, 2]) / (m - m0)) for k, v in res.items()}
for k, (fw, fd, pw) in summ.items():
    print(f"   {k:10s} FWER {fw:.3f}  FDR {fd:.3f}  power {pw:.3f}")
assert summ["bonferroni"][0] <= .065 and summ["holm"][0] <= .065 and summ["bh"][1] <= .05 + .01
assert summ["bh"][0] > .05 and summ["bh"][2] > summ["holm"][2] >= summ["bonferroni"][2]
p3 = stats.norm.sf(3); p5 = stats.norm.sf(5)
print(f"Look-elsewhere: local 3 sigma p = {p3:.2e}; over 100 independent bins global p = {1-(1-p3)**100:.3f}; 5 sigma p = {p5:.2e}")

# %% [markdown]
# ## Part 10: Numbers for the exercises

# %%
h = stats.t.ppf(.975, 15) * 4.8 / 4; cl, ch = stats.chi2.ppf([.025, .975], 15)
print(f"A1: t_15 = {stats.t.ppf(.975, 15):.4f}, t CI [{52.3-h:.3f}, {52.3+h:.3f}]; z CI [{52.3-Z*5/4:.3f}, {52.3+Z*5/4:.3f}]; "
      f"chi2_15 {cl:.3f}, {ch:.3f}; sigma in [{np.sqrt(15*4.8**2/ch):.3f}, {np.sqrt(15*4.8**2/cl):.3f}]")
print(f"A2: Wilson for 0/10 = [{wilson_prop(0, 10)[0]:.4f}, {wilson_prop(0, 10)[1]:.4f}]")
print(f"A3: size P(X>=8|.5) = {stats.binom.sf(7, 10, .5):.4f}, power at .7 = {stats.binom.sf(7, 10, .7):.4f}; "
      f"X>=9: size {stats.binom.sf(8, 10, .5):.4f}, power {stats.binom.sf(8, 10, .7):.4f}")
pa4 = np.array([0.001, 0.012, 0.020, 0.024, 0.040, 0.046, 0.30, 0.70])
assert [int(f(pa4).sum()) for f in (bonferroni, holm, bh)] == [1, 1, 4]
print(f"A6: n = {((Z + stats.norm.ppf(.9)) / .3) ** 2:.2f}")
life = np.array([12.1, 3.4, 27.9, 8.6, 15.2, 5.3, 41.0, 9.5]); S, n = life.sum(), 8
ql, qh = stats.chi2.ppf([.025, .975], 16); thr = stats.chi2.ppf(.95, 16) / (2 * .1)
print(f"B2: sum {S:.1f}, lambda-hat {n/S:.5f}, chi2_16 {ql:.3f}, {qh:.3f}; lambda in [{ql/(2*S):.5f}, {qh/(2*S):.5f}], "
      f"mean life in [{2*S/qh:.2f}, {2*S/ql:.2f}]")
print(f"B3: reject if sum > {thr:.2f}; p = {stats.chi2.sf(2*.1*S, 16):.4f}; power at lambda=0.05: {stats.chi2.sf(2*.05*thr, 16):.4f}")
lh = n / S; LRT = 2 * (n * np.log(lh / .1) - (lh - .1) * S); Wl = n * (1 - .1 / lh) ** 2; Sc = (n / .1 - S) ** 2 / (n / .01)
th = S / n; Wth = n * (th - 10) ** 2 / th ** 2
pex = 2 * min(stats.chi2.sf(2 * .1 * S, 16), stats.chi2.cdf(2 * .1 * S, 16))
print(f"B4: LRT {LRT:.4f} (p {stats.chi2.sf(LRT,1):.3f}), Wald(lambda) {Wl:.4f} (p {stats.chi2.sf(Wl,1):.3f}), "
      f"score {Sc:.4f}, Wald(theta) {Wth:.4f} (p {stats.chi2.sf(Wth,1):.3f}); exact two-sided p {pex:.3f}")
assert np.isclose(Wl, Sc)
xb5 = np.array([np.mean(c) for c in itertools.product([1, 2, 6], repeat=3)])
assert np.isclose(xb5.mean(), 3) and np.isclose(xb5.var(), 14 / 9)
print(f"B5: 27 resamples, var(xbar*) = {xb5.var():.4f} = 14/9, se_boot {np.sqrt(14/9):.4f} vs s/sqrt(n) {np.sqrt(7/3):.4f}")
Dl = 2.5; zu = stats.norm.isf([.05, .04, .001])
print(f"B7: P(p<=.05) = {stats.norm.cdf(Dl - zu[0]):.4f}; likelihood ratio of p at .04: {np.exp(Dl*zu[1]-Dl**2/2):.3f}, at .001: {np.exp(Dl*zu[2]-Dl**2/2):.1f}")
for nn in (10, 100, 1000, 10 ** 4, 10 ** 6):
    print(f"C3: n = {nn:>7d}, z = 2.5: BF01 = {np.sqrt(1 + nn) * np.exp(-2.5 ** 2 / 2 * nn / (1 + nn)):.4f}")
# C1 check: the change-of-variables integral reproduces the t density.
kk = 4
for tv in (0.0, 1.3, 3.0):
    f = lambda v: stats.norm.pdf(tv * np.sqrt(v / kk)) * np.sqrt(v / kk) * stats.chi2.pdf(v, kk)
    assert np.isclose(integrate.quad(f, 0, np.inf)[0], stats.t.pdf(tv, kk))

# %% [markdown]
# ## Your turn
#
# 1. (B8) Bootstrap the slope of a regression y = a + b x + noise with n = 25.
#    Resample (x, y) pairs, refit, and build percentile and basic intervals.
#    Check coverage against the classical t interval for the slope (ST-3 derives it).
#    A reference solution is below; write yours first, then delete it.
# 2. (C2) Optional stopping: under H0, run a z test after every 10 observations
#    up to 100 and stop at the first p < 0.05. Estimate the type I error.
#    Reference below. Then try a Bonferroni-style 0.05/10 per look: too conservative?
# 3. Look-elsewhere: simulate a background-only histogram of 100 Poisson bins
#    with mean 100 each, find the most significant local excess, and estimate the
#    distribution of the minimum local p-value. Compare with 1 - (1 - p)^100.
# 4. Replace Exponential data in Part 3 with Pareto(alpha = 1.5) (infinite variance):
#    what happens to the bootstrap coverage for the mean as n grows?

# %% [markdown]
# ## Reference solutions for B8 and C2 (write your own first)

# %%
def slope(a):
    """Least-squares slope along axis 1 for pairs a[..., 0] = x, a[..., 1] = y."""
    x = a[..., 0] - a[..., 0].mean(1, keepdims=True); y = a[..., 1] - a[..., 1].mean(1, keepdims=True)
    return (x * y).sum(1) / (x * x).sum(1)


hit = {"percentile": 0, "basic": 0, "classical t": 0}
for _ in range(Rb):
    xr = rng.uniform(0, 10, 25); yr = 1 + 0.5 * xr + rng.normal(0, 2, 25); xy = np.c_[xr, yr]
    for k, (a, b) in bootstrap_ci(xy, slope, Bb, kinds=("percentile", "basic")).items():
        hit[k] += a <= 0.5 <= b
    lr_ = stats.linregress(xr, yr); hw = stats.t.ppf(.975, 23) * lr_.stderr
    hit["classical t"] += abs(lr_.slope - 0.5) <= hw
print("B8 slope coverage:", {k: v / Rb for k, v in hit.items()})
assert all(v / Rb > 0.88 for v in hit.values())

Xs = rng.standard_normal((20_000, 100)); looks = np.arange(10, 101, 10)
zs = Xs.cumsum(1)[:, looks - 1] / np.sqrt(looks)
alpha_os = np.mean((np.abs(zs) > Z).any(1)); alpha_bonf = np.mean((np.abs(zs) > stats.norm.isf(.0025)).any(1))
print(f"C2: optional stopping with 10 looks: type I error {alpha_os:.3f}; with 0.05/10 per look: {alpha_bonf:.3f}")
assert 0.17 < alpha_os < 0.22 and alpha_bonf < 0.05

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
