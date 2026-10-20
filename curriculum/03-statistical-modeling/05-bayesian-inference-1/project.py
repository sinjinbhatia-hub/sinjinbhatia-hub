# %% [markdown]
# # ST-5 Project: Bayesian Inference I
#
# Conjugate updating, predictive distributions, priors, credible intervals, Bayesian regression, the Laplace
# approximation and model comparison. Every number in notes.md Section 4 and solutions.md is asserted or printed
# here. Parts: 1 updating; 2 predictive vs plug-in; 3 coverage; 4 priors (SymPy); 5 Normal-Inverse-Gamma;
# 6 exponential-family conjugacy; 7 decisions, MAP = ridge/lasso; 8 Bayesian regression; 9 Laplace; 10 evidence,
# Occam, BIC; 11 Jeffreys-Lindley; 12 posterior predictive checks; 13 stopping rules; 14 numbers for the exercises.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats, special, optimize, integrate
import sympy as sp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
close = lambda a, b, tol=1e-9: bool(np.allclose(a, b, rtol=tol, atol=tol))


def save(fig, name):
    fig.tight_layout(); fig.savefig(FIG / name, dpi=110); plt.close(fig)


# %% [markdown]
# ## Part 1: Conjugate updating as data arrive (notes 3.2, 3.3; Examples 4.1, 4.2)
# %%
def beta_update(a, b, y, n):
    """Beta(a,b) prior, y successes in n trials -> Beta(a+y, b+n-y), eq. (3.4)."""
    return a + y, b + n - y


def normal_update(mu0, tau2, ybar, n, s2):
    """N(mu0,tau2) prior, n observations with mean ybar, known variance s2: precisions add, eq. (3.10)."""
    prec = 1 / tau2 + n / s2
    return (mu0 / tau2 + n * ybar / s2) / prec, 1 / prec


def hpd_beta(a, b, level=0.95):
    """Highest-posterior-density interval of a unimodal Beta: the shortest interval with mass `level` (Theorem 3.9)."""
    d = stats.beta(a, b)
    r = optimize.minimize_scalar(lambda p: d.ppf(p + level) - d.ppf(p), bounds=(0, 1 - level),
                                 method="bounded", options={"xatol": 1e-12})
    return d.ppf(r.x), d.ppf(r.x + level)


a0 = b0 = 2.0                                                    # Example 4.1: Beta(2,2) prior, 7 of 10
y, n = 7, 10
a1, b1 = beta_update(a0, b0, y, n)
post = stats.beta(a1, b1)
w = n / (a0 + b0 + n)
assert (a1, b1) == (9, 5) and close(post.mean(), w * y / n + (1 - w) * a0 / (a0 + b0))        # shrinkage form (3.6)
assert close(post.var(), a1 * b1 / ((a1 + b1) ** 2 * (a1 + b1 + 1)))                          # (3.7)
eq_ci, hpd = post.interval(0.95), hpd_beta(a1, b1)
assert close(post.pdf(hpd[0]), post.pdf(hpd[1]), 1e-5)                                         # HPD: equal density at the ends
log_ev = (special.gammaln(n + 1) - special.gammaln(y + 1) - special.gammaln(n - y + 1)
          + special.betaln(a1, b1) - special.betaln(a0, b0))                                   # evidence (3.5)
assert close(np.exp(log_ev), integrate.quad(lambda p: stats.binom.pmf(y, n, p) * stats.beta.pdf(p, a0, b0), 0, 1)[0], 1e-8)
assert all(close(stats.betabinom.pmf(k, n, 1, 1), 1 / (n + 1)) for k in range(n + 1))          # uniform prior: 1/(n+1)
assert close(integrate.quad(lambda t: t * post.pdf(t), 0, 1)[0], a1 / (a1 + b1))               # P(next success) = mean
assert close(beta_update(1, 1, 7, 10)[0] / sum(beta_update(1, 1, 7, 10)), 8 / 12)             # Laplace's rule (3.9)
print(f"Ex 4.1: Beta({a1:.0f},{b1:.0f}), mean {post.mean():.4f}, sd {post.std():.4f}, w={w:.4f}; 95% equal-tailed "
      f"({eq_ci[0]:.4f}, {eq_ci[1]:.4f}), HPD ({hpd[0]:.4f}, {hpd[1]:.4f}); evidence {np.exp(log_ev):.5f} (uniform {1 / 11:.5f})")

yN = np.array([0.5, 3.1, 1.2, 2.4, -0.3, 2.9, 0.8, 1.4])         # Example 4.2: known sigma^2 = 4, prior N(0,1)
s2, mu0, tau02 = 4.0, 0.0, 1.0
mun, taun2 = normal_update(mu0, tau02, yN.mean(), len(yN), s2)
assert close(yN.mean(), 1.5) and close(1 / taun2, 1 / tau02 + 8 / s2) and close([mun, taun2], [1.0, 1 / 3])
print(f"Ex 4.2: N({mun:.4f}, {taun2:.4f}); predictive variance {s2 + taun2:.4f}; 95% interval "
      f"({mun - 1.96 * taun2 ** .5:.4f}, {mun + 1.96 * taun2 ** .5:.4f})")

p_true, mu_true = 0.7, 1.5                                       # figure grid; sequential updating equals batch
bern = (rng.random(100) < p_true).astype(int)
norm = rng.normal(mu_true, np.sqrt(s2), 100)
marks = [0, 1, 3, 10, 30, 100]
fig, axs = plt.subplots(2, len(marks), figsize=(16, 5.5), sharey="row")
pg, mg = np.linspace(0, 1, 400), np.linspace(-3, 4, 400)
sa, sb, sm, st = 2.0, 2.0, 0.0, 1.0
for i in range(101):
    if i in marks:
        k = marks.index(i)
        assert close((sa, sb), beta_update(2, 2, bern[:i].sum(), i))
        assert close((sm, st), normal_update(0, 1, norm[:i].mean() if i else 0.0, i, s2))
        axs[0, k].plot(pg, stats.beta.pdf(pg, sa, sb)); axs[0, k].axvline(p_true, c="k", ls=":"); axs[0, k].set_title(f"n={i}: Beta({sa:.0f},{sb:.0f})", fontsize=9)
        axs[1, k].plot(mg, stats.norm.pdf(mg, sm, np.sqrt(st))); axs[1, k].axvline(mu_true, c="k", ls=":"); axs[1, k].set_title(f"n={i}: N({sm:.2f},{st:.3f})", fontsize=9)
    if i < 100:
        sa, sb = beta_update(sa, sb, bern[i], 1)
        sm, st = normal_update(sm, st, norm[i], 1, s2)
fig.suptitle("Posteriors sharpen as data arrive (dotted = truth); sequential = batch")
save(fig, "p1_updating.png")

# %% [markdown]
# ## Part 2: Posterior predictive vs plug-in (notes 3.2, 3.3, 3.10)
# %%
m = 10                                                           # ten future trials, Beta-Binomial predictive (3.8)
pred = stats.betabinom(m, a1, b1)
ytil = rng.binomial(m, post.rvs(400_000, random_state=rng))      # sample theta, then y~
assert np.abs(np.bincount(ytil, minlength=m + 1) / ytil.size - pred.pmf(np.arange(m + 1))).max() < 0.004
assert close(pred.var(), m * a1 * b1 * (a1 + b1 + m) / ((a1 + b1) ** 2 * (a1 + b1 + 1)))
plug_mle, plug_mean = stats.binom(m, y / n), stats.binom(m, a1 / (a1 + b1))
assert pred.var() > plug_mean.var() > plug_mle.var()
print(f"Part 2: predictive var of 10 future trials {pred.var():.4f} (MC {ytil.var():.4f}) vs plug-in "
      f"{plug_mle.var():.3f} (MLE), {plug_mean.var():.4f} (posterior mean)")
ynew = rng.normal(rng.normal(mun, np.sqrt(taun2), 400_000), np.sqrt(s2))                       # Normal predictive (3.11)
assert abs(ynew.var() - (s2 + taun2)) < 0.03 and abs(ynew.mean() - mun) < 0.01

fig, ax = plt.subplots(figsize=(6, 3.5)); ks = np.arange(m + 1)
ax.bar(ks - .2, pred.pmf(ks), .4, label="posterior predictive"); ax.bar(ks + .2, plug_mle.pmf(ks), .4, label="plug-in MLE")
ax.legend(); ax.set_title("10 future trials: plug-in is too confident"); save(fig, "p2_predictive.png")

# %% [markdown]
# ## Part 3: Credible intervals and coverage (notes 3.8)
# %%
def eq_tail(a, b, level=0.95):
    lo = (1 - level) / 2
    return stats.beta.ppf(lo, a, b), stats.beta.ppf(1 - lo, a, b)


N, nn = 400_000, 20                                              # Theorem 3.10: theta from the prior -> 95% covers 95%
th = rng.beta(2, 2, N); yy = rng.binomial(nn, th)
lo, hi = eq_tail(2 + yy, 2 + nn - yy)
cov = np.mean((lo <= th) & (th <= hi))
assert abs(cov - 0.95) < 4 * np.sqrt(0.95 * 0.05 / N)
print(f"Part 3a: theta~Beta(2,2), n=20: coverage over prior and data {cov:.4f}")


def freq_coverage(theta, n, a, b):
    """Exact P_theta(theta in equal-tailed credible interval), summing over y = 0..n (no Monte Carlo)."""
    ys = np.arange(n + 1)
    lo, hi = eq_tail(a + ys, b + n - ys)
    return float(np.sum(stats.binom.pmf(ys, n, theta) * ((lo <= theta) & (theta <= hi))))


ns = [5, 10, 20, 50, 100, 200, 400, 1000]
priors = {"Jeffreys Beta(.5,.5)": (.5, .5), "uniform Beta(1,1)": (1, 1), "strong Beta(10,10)": (10, 10)}
cov_tab = {}
fig, axs = plt.subplots(1, 2, figsize=(10, 3.5), sharey=True)
for ax, t0 in zip(axs, [0.5, 0.1]):
    for name, (a, b) in priors.items():
        cov_tab[(name, t0)] = [freq_coverage(t0, k, a, b) for k in ns]
        ax.semilogx(ns, cov_tab[(name, t0)], "o-", label=name)
    ax.axhline(.95, c="k", ls=":"); ax.set_title(f"true theta = {t0}"); ax.set_xlabel("n")
    print(f"Part 3b: theta={t0}, coverage at n=10,100,1000:", {k[0][:7]: [round(v[ns.index(j)], 3) for j in (10, 100, 1000)] for k, v in cov_tab.items() if k[1] == t0})
axs[0].set_ylabel("exact frequentist coverage"); axs[0].legend(fontsize=8)
save(fig, "p3_coverage.png")
sc_ = cov_tab[("strong Beta(10,10)", 0.1)]
assert sc_[1] < 0.6 and cov_tab[("strong Beta(10,10)", 0.5)][1] > 0.97     # strong prior: wrong place under-covers, right place over-covers
assert sc_[1] < sc_[3] < sc_[5] < sc_[-1] < 0.93                           # and its damage fades slowly with n
assert all(abs(cov_tab[(k, 0.1)][-1] - .95) < .03 for k in ("Jeffreys Beta(.5,.5)", "uniform Beta(1,1)"))   # weak priors agree by n=1000
avg_cov = np.mean([freq_coverage(t, 20, 1, 1) for t in np.linspace(0.0005, 0.9995, 4000)])
assert abs(avg_cov - 0.95) < 0.003                                         # coverage averaged over the prior = calibration
print(f"  uniform prior, n=20: coverage averaged over theta~U(0,1) = {avg_cov:.4f}")

# %% [markdown]
# ## Part 4: Priors, Jeffreys, reparametrization (notes 3.7; Example 4.5)
# %%
p_, lam, mu_, sg, rr = sp.symbols("p lambda mu sigma r", positive=True)    # r stands for (y - mu)^2, whose mean is sigma^2
yv = sp.symbols("y", real=True)
fisher = lambda logp, th, repl: sp.simplify(-sp.diff(logp, th, 2).subs(repl))    # -E[d^2 log p] when linear in y or r (ST-1 (3.11))
I_bern = fisher(yv * sp.log(p_) + (1 - yv) * sp.log(1 - p_), p_, {yv: p_})
I_pois = fisher(yv * sp.log(lam) - lam, lam, {yv: lam})
I_mu = fisher(-(yv - mu_) ** 2 / (2 * sg ** 2), mu_, {yv: mu_})
I_sig = fisher(-rr / (2 * sg ** 2) - sp.log(sg), sg, {rr: sg ** 2})
assert sp.simplify(I_bern - 1 / (p_ * (1 - p_))) == 0 and sp.simplify(I_pois - 1 / lam) == 0
assert sp.simplify(I_mu - 1 / sg ** 2) == 0 and sp.simplify(I_sig - 2 / sg ** 2) == 0
print("Part 4: Fisher information  Bernoulli", I_bern, "| Poisson", I_pois, "| Normal mean", I_mu, "| Normal sigma", I_sig)

eta = np.linspace(-8, 8, 801); pp = special.expit(eta)
flat_eta = stats.beta.pdf(pp, 1, 1) * pp * (1 - pp)                        # flat in p, pushed to logit scale (3.21)
jeff_eta = stats.beta.pdf(pp, .5, .5) * pp * (1 - pp)                      # Jeffreys in p, pushed to logit scale
ew = np.linspace(-80, 80, 160001)                                          # tails decay like e^(-|eta|/2): integrate wide
assert close(flat_eta, stats.logistic.pdf(eta)) and close(jeff_eta, 1 / (2 * np.pi * np.cosh(eta / 2)))
assert close(jeff_eta, np.sqrt(pp * (1 - pp)) / np.pi) and close(np.trapezoid(1 / (2 * np.pi * np.cosh(ew / 2)), ew), 1, 1e-6)   # = sqrt(I_eta)/pi: (3.23)
mode_p = (a1 - 1) / (a1 + b1 - 2)                                          # the mode is not invariant
sig_ = special.expit
r = optimize.minimize_scalar(lambda e: -stats.beta.pdf(sig_(e), a1, b1) * sig_(e) * (1 - sig_(e)), bounds=(-5, 5), method="bounded", options={"xatol": 1e-10})
assert close(mode_p, 2 / 3) and close(sig_(r.x), 9 / 14, 1e-6)
print(f"  Beta(9,5): mode in p {mode_p:.4f}; mode on the logit scale maps back to p={sig_(r.x):.4f}")
fl, jf = stats.beta(1, 11), stats.beta(.5, 10.5)                           # Example 4.5: no successes in 10 trials
assert close(stats.gamma(.5, scale=.1).mean(), 0.05)                       # Poisson Jeffreys, zero counts in 10 intervals
print(f"Ex 4.5: y=0, n=10: uniform HPD upper {fl.ppf(.95):.4f}, mean {fl.mean():.4f}; Jeffreys mean {jf.mean():.4f}, "
      f"95% eq-tailed ({jf.ppf(.025):.5f}, {jf.ppf(.975):.4f})")
for sd in (10, 1.5):                                                       # prior predictive check on the log-odds scale
    extreme = np.mean(np.abs(rng.normal(0, sd, 400_000)) > special.logit(.99))
    assert abs(extreme - 2 * stats.norm.cdf(-special.logit(.99) / sd)) < 0.004
    print(f"  prior N(0,{sd}^2) on log-odds: P(p<.01 or p>.99) = {extreme:.4f}")
fig, axs = plt.subplots(1, 3, figsize=(13, 3.5))
pl = np.linspace(0.001, .999, 500)
axs[0].plot(pl, stats.beta.pdf(pl, 1, 1), label="flat"); axs[0].plot(pl, stats.beta.pdf(pl, .5, .5), label="Jeffreys"); axs[0].set(ylim=(0, 4), title="priors on p"); axs[0].legend()
axs[1].plot(eta, flat_eta, label="from flat p (logistic)"); axs[1].plot(eta, jeff_eta, label="from Jeffreys"); axs[1].set_title("same priors on logit(p)"); axs[1].legend(fontsize=8)
for sd in (10, 1.5):
    axs[2].hist(special.expit(rng.normal(0, sd, 20000)), 40, alpha=.5, density=True, label=f"N(0,{sd}^2) on log-odds")
axs[2].set_title("prior predictive: implied p"); axs[2].legend(fontsize=8)
save(fig, "p4_priors.png")

# %% [markdown]
# ## Part 5: Normal-Inverse-Gamma, the t marginal, and ST-2's interval (notes 3.4; Example 4.3)
# %%
g = np.array([9.79, 9.83, 9.76, 9.85, 9.81, 9.78, 9.88, 9.80, 9.74, 9.82])   # ST-2 Example 4.2


def nig_update(mu0, k0, a0, b0, y):
    """NIG hyperparameter update, eq. (3.13): returns (mu_n, kappa_n, alpha_n, beta_n)."""
    y = np.asarray(y, float); n = y.size; yb = y.mean(); S = ((y - yb) ** 2).sum()
    kn = k0 + n
    return (k0 * mu0 + n * yb) / kn, kn, a0 + n / 2, b0 + 0.5 * (S + k0 * n * (yb - mu0) ** 2 / kn)


def nig_t(mun, kn, an, bn, predictive=False):
    """Marginal t for mu (3.14) or for a new observation (3.15)."""
    return stats.t(df=2 * an, loc=mun, scale=np.sqrt(bn / (an * kn) * ((kn + 1) if predictive else 1)))


n_g, S_g = g.size, ((g - g.mean()) ** 2).sum()
ref = nig_update(0.0, 0.0, -0.5, 0.0, g)                                   # reference prior 1/sigma^2: formal limit of the NIG
assert close(ref[:2], (g.mean(), n_g)) and close(ref[2], (n_g - 1) / 2) and close(ref[3], S_g / 2)
ref_ci = nig_t(*ref).interval(0.95)
t_ci = g.mean() + np.array([-1, 1]) * stats.t.ppf(.975, n_g - 1) * g.std(ddof=1) / np.sqrt(n_g)   # ST-2's t interval
assert close(ref_ci, t_ci, 1e-10) and close(np.round(t_ci, 4), [9.7762, 9.8358])
ig = stats.invgamma(ref[2], scale=ref[3]); sig2_ci = (ig.ppf(.025), ig.ppf(.975))
assert close(sig2_ci, (S_g / stats.chi2.ppf(.975, 9), S_g / stats.chi2.ppf(.025, 9)))                # = ST-2 (3.4)
prior = (9.5, 1.0, 2.0, 0.02)                                              # a proper weakly informative NIG
post_n = nig_update(*prior, g)
seq = prior
for gi in g:
    seq = nig_update(*seq, [gi])
assert close(seq, post_n)                                                  # sequential = batch
tm, tp = nig_t(*post_n), nig_t(*post_n, predictive=True)
mus, s2s = np.linspace(9.45, 10.15, 1401), np.geomspace(1e-4, 0.05, 4000)  # brute-force (mu, sigma^2) grid for the t marginal
M, V = np.meshgrid(mus, s2s, indexing="ij")
logp = (-0.5 * n_g * np.log(2 * np.pi * V) - (S_g + n_g * (g.mean() - M) ** 2) / (2 * V)
        + stats.norm.logpdf(M, prior[0], np.sqrt(V / prior[1])) + stats.invgamma.logpdf(V, prior[2], scale=prior[3]))
marg = np.trapezoid(np.exp(logp - logp.max()), s2s, axis=1)
marg /= np.trapezoid(marg, mus)
assert np.abs(marg - tm.pdf(mus)).max() < 0.003 * tm.pdf(mus).max()
jt = stats.t(n_g, g.mean(), np.sqrt(S_g) / n_g)                            # joint Jeffreys prior: t with n (not n-1) df
mu_pts = g.mean() + np.array([0, .03, -.05])
num = [integrate.quad(lambda s: s ** (-n_g - 2) * np.exp(-(S_g + n_g * (m_ - g.mean()) ** 2) / (2 * s * s)), 1e-4, 5)[0] for m_ in mu_pts]
assert close(np.array(num) / num[0], jt.pdf(mu_pts) / jt.pdf(g.mean()), 1e-5)
print(f"Ex 4.3: reference 95% ({ref_ci[0]:.4f}, {ref_ci[1]:.4f}) = ST-2; sigma^2 ({sig2_ci[0]:.5f}, {sig2_ci[1]:.5f}), sigma ({sig2_ci[0] ** .5:.4f}, {sig2_ci[1] ** .5:.4f})")
print(f"  NIG{prior} -> (mu_n,k_n,a_n,b_n)=({post_n[0]:.4f},{post_n[1]:.0f},{post_n[2]:.0f},{post_n[3]:.6f}); mu 95% ({tm.interval(.95)[0]:.4f}, "
      f"{tm.interval(.95)[1]:.4f}); predictive ({tp.interval(.95)[0]:.4f}, {tp.interval(.95)[1]:.4f}); joint-Jeffreys half-width "
      f"{stats.t.ppf(.975, n_g) * np.sqrt(S_g) / n_g:.4f} vs reference {(ref_ci[1] - ref_ci[0]) / 2:.4f}")

# %% [markdown]
# ## Part 6: Exponential-family conjugacy reproduces the three cases (notes 3.5, 3.6; Example 4.4)
# %%
def ef_update(chi, nu, T_sum, n):
    """Conjugate update for the prior exp(eta.chi - nu A(eta)): chi -> chi + sum T, nu -> nu + n, eq. (3.19)."""
    return chi + T_sum, nu + n


yP = np.array([3, 5, 2, 4, 6])                                             # Poisson-Gamma: (chi, nu) = (alpha0, beta0)
chi, nu = ef_update(2.0, 1.0, yP.sum(), yP.size)
assert (chi, nu) == (22.0, 6.0)
lamg = np.linspace(1e-3, 12, 4000)
gp = np.exp(stats.gamma.logpdf(lamg, 2, scale=1) + stats.poisson.logpmf(yP[:, None], lamg).sum(0)); gp /= np.trapezoid(gp, lamg)
assert np.abs(gp - stats.gamma.pdf(lamg, chi, scale=1 / nu)).max() < 1e-4                   # grid posterior = Gamma(22, 6)
Zq = lambda c, v: integrate.quad(lambda e: np.exp(e * c - v * np.exp(e)), -40, 6)[0]        # normalizer of the prior in eta = log(lambda)
assert close(Zq(2, 1), special.gamma(2), 1e-6) and close(integrate.quad(lambda e: np.exp(e) * np.exp(e * chi - nu * np.exp(e)), -40, 6)[0] / Zq(chi, nu), chi / nu, 1e-6)   # (3.20)
nb = stats.nbinom(chi, nu / (nu + 1))                                                       # predictive (3.17)
mc = rng.poisson(rng.gamma(chi, 1 / nu, 400_000))
assert close(nb.mean(), chi / nu) and close(nb.var(), chi * (nu + 1) / nu ** 2)
assert abs(mc.mean() - nb.mean()) < 0.02 and abs(mc.var() - nb.var()) < 0.06
assert np.abs(np.bincount(mc, minlength=16)[:16] / mc.size - nb.pmf(np.arange(16))).max() < 0.002
print(f"Ex 4.4: (chi,nu)=({chi:.0f},{nu:.0f}) -> Gamma(22,6), mean {chi / nu:.4f}; NB predictive mean {nb.mean():.4f}, var {nb.var():.4f}, P(0)={nb.pmf(0):.4f}")
assert close(integrate.quad(lambda e: np.exp(e * 2.0 - 4.0 * np.log1p(np.exp(e))), -40, 40)[0], special.beta(2, 2), 1e-6)   # Bernoulli: chi=2, nu=4 <-> Beta(2,2)
chiN, nuN = ef_update(s2 / tau02 * mu0, s2 / tau02, yN.sum(), yN.size)                       # Normal, sigma known: nu = sigma^2/tau0^2
assert close(chiN / nuN, mun) and close(s2 / nuN, taun2)
c1, c2, nu_ = 3.0, 9.0, 4.0                                                  # N(mu,s2): the exp-family prior is NIG with kappa0 = nu, alpha0 = (nu+3)/2
mg_, vg0 = np.linspace(-1, 3, 7)[:, None], np.linspace(0.3, 3, 9)[None, :]
e1, e2 = mg_ / vg0, -1 / (2 * vg0)
dens = np.exp(e1 * c1 + e2 * c2 - nu_ * (-e1 ** 2 / (4 * e2) - 0.5 * np.log(-2 * e2))) / (2 * vg0 ** 3)   # Jacobian d(eta)/d(mu,s2) = 1/(2 s2^3)
ratio = dens / (stats.norm.pdf(mg_, c1 / nu_, np.sqrt(vg0 / nu_)) * stats.invgamma.pdf(vg0, (nu_ + 3) / 2, scale=(c2 - c1 ** 2 / nu_) / 2))
assert np.ptp(ratio) < 1e-9 * ratio.mean()                                   # constant ratio: same density

# %% [markdown]
# ## Part 7: Decisions, MAP, ridge and lasso (notes 3.9; Example 4.6)
# %%
d = stats.beta(a1, b1)
risk = lambda act, loss: integrate.quad(lambda t: loss(act, t) * d.pdf(t), 0, 1)[0]
best_sq = optimize.minimize_scalar(lambda a: risk(a, lambda a, t: (a - t) ** 2), bounds=(0, 1), method="bounded", options={"xatol": 1e-9}).x
best_ab = optimize.minimize_scalar(lambda a: risk(a, lambda a, t: abs(a - t)), bounds=(0, 1), method="bounded", options={"xatol": 1e-9}).x
assert abs(best_sq - d.mean()) < 1e-5 and abs(best_ab - d.median()) < 1e-4
print(f"Part 7: Beta(9,5): squared-loss action {best_sq:.4f} = mean; absolute-loss action {best_ab:.4f} = median; MAP {mode_p:.4f}")
xc, yc = np.array([-2., -1, 0, 1, 2]), np.array([-2., 0, 1, 0, 1])           # Example 4.6: ST-3's five points, centered
Sxx, Sxy, sig2, tau2 = xc @ xc, xc @ yc, 0.8, 0.5
lam_r = sig2 / tau2                                                          # (3.27)
map_slope = optimize.minimize_scalar(lambda b: ((yc - b * xc) ** 2).sum() / (2 * sig2) + b * b / (2 * tau2)).x
assert close(map_slope, Sxy / (Sxx + lam_r), 1e-6) and close([Sxx, Sxy], [10, 6])
print(f"Ex 4.6: OLS slope {Sxy / Sxx:.4f}, lambda={lam_r:.2f}, MAP = ridge slope {Sxy / (Sxx + lam_r):.4f} (factor {Sxx / (Sxx + lam_r):.4f}); "
      f"posterior var {sig2 / (Sxx + lam_r):.5f} vs OLS var {sig2 / Sxx:.4f}")
soft = lambda z, l: np.sign(z) * max(abs(z) - l, 0)                          # lasso = Laplace-prior MAP (3.28); y|b ~ N(b,1), prior exp(-|b|/sc)
for z, sc in [(1.2, 2.0), (0.3, 2.0)]:
    nlp = lambda b: 0.5 * (z - b) ** 2 + abs(b) / sc
    b_map = optimize.minimize_scalar(nlp, bounds=(-3, 3), method="bounded", options={"xatol": 1e-10}).x
    pm = integrate.quad(lambda b: b * np.exp(-nlp(b)), -12, 12, points=[0])[0] / integrate.quad(lambda b: np.exp(-nlp(b)), -12, 12, points=[0])[0]
    assert abs(b_map - soft(z, 1 / sc)) < 1e-5
    print(f"  Laplace prior scale {sc}, z={z}: MAP {b_map:.4f} = soft-threshold {soft(z, 1 / sc):.4f}; posterior mean {pm:.4f}")

# %% [markdown]
# ## Part 8: Bayesian linear regression and the predictive band (notes 3.11; Example 4.7)
# %%
def blr_posterior(X, y, s2, mu0, S0):
    """beta ~ N(mu0,S0), y|beta ~ N(X beta, s2 I) -> posterior N(mu_n, S_n), eq. (3.30)."""
    P0 = np.linalg.inv(S0)
    Sn = np.linalg.inv(P0 + X.T @ X / s2)
    return Sn @ (P0 @ mu0 + X.T @ y / s2), Sn


def blr_predict(Xs, mun, Sn, s2):
    """Predictive mean x*'mu_n and variance s2 + x*' S_n x*, eq. (3.32)."""
    return Xs @ mun, s2 + np.einsum("ij,jk,ik->i", Xs, Sn, Xs)


def nig_reg(X, y, mu0, V0, a0, b0):
    """beta|s2 ~ N(mu0, s2 V0), s2 ~ IG(a0,b0): posterior (mu_n, V_n, a_n, b_n); predictive is t_{2a_n} (3.33)."""
    Vn = np.linalg.inv(np.linalg.inv(V0) + X.T @ X)
    mn = Vn @ (np.linalg.solve(V0, mu0) + X.T @ y)
    return mn, Vn, a0 + len(y) / 2, b0 + 0.5 * (y @ y + mu0 @ np.linalg.solve(V0, mu0) - mn @ np.linalg.solve(Vn, mn))


def nig_reg_logev(X, y, mu0, V0, a0, b0):
    """log p(y) for the NIG regression model, eq. (3.39): a ratio of normalizers."""
    mn, Vn, an, bn = nig_reg(X, y, mu0, V0, a0, b0)
    return (special.gammaln(an) - special.gammaln(a0) + 0.5 * (np.linalg.slogdet(Vn)[1] - np.linalg.slogdet(V0)[1])
            + a0 * np.log(b0) - an * np.log(bn) - 0.5 * len(y) * np.log(2 * np.pi))


xr, yr = np.arange(1., 6), np.array([2., 4, 5, 4, 5])                        # ST-3 Example 4.1
X = np.column_stack([np.ones(5), xr])
s2r, tau2r = 0.8, 10.0
mun_r, Sn_r = blr_posterior(X, yr, s2r, np.zeros(2), tau2r * np.eye(2))
lam2 = s2r / tau2r
assert close(mun_r, np.linalg.solve(X.T @ X + lam2 * np.eye(2), X.T @ yr))   # ridge, ST-3 (3.39): (3.31)
assert close(Sn_r, s2r * np.linalg.inv(X.T @ X + lam2 * np.eye(2)))
assert close(optimize.minimize(lambda b: ((yr - X @ b) ** 2).sum() / (2 * s2r) + b @ b / (2 * tau2r), np.zeros(2), tol=1e-12).x, mun_r, 1e-5)   # MAP = mean
P0 = tau2r * np.eye(2)                                                       # F5 (3.19) on the joint of (beta, y)
assert close(P0 @ X.T @ np.linalg.solve(X @ P0 @ X.T + s2r * np.eye(5), yr), mun_r)
xs = np.array([[1., 3.], [1., 8.]])
mp, vp = blr_predict(xs, mun_r, Sn_r, s2r)
big_mu, big_S = blr_posterior(X, yr, s2r, np.zeros(2), 1e8 * np.eye(2))      # flat-prior limit
_, vflat = blr_predict(xs, big_mu, big_S, s2r)
h = 1 / 5 + (xs[:, 1] - 3) ** 2 / 10
assert close(vflat, s2r * (1 + h), 1e-6) and np.allclose(big_mu, [2.2, 0.6], atol=1e-5)    # = ST-3 prediction variance
print(f"Ex 4.7: posterior mean {mun_r.round(4)}, sd {np.sqrt(np.diag(Sn_r)).round(4)}; predictive mean/var at x=3: {mp[0]:.4f}/{vp[0]:.4f}, "
      f"at x=8: {mp[1]:.4f}/{vp[1]:.4f}; flat-prior var {vflat.round(4)} = 0.8(1+h), h={h.round(2)}")
V0, a0r, b0r = 10 * np.eye(2), 2.0, 1.0
chain = 0.0                                                                  # chain rule: log p(y) = sum_i log p(y_i | y_<i)
for i in range(5):
    mn_, Vn_, an_, bn_ = (np.zeros(2), V0, a0r, b0r) if i == 0 else nig_reg(X[:i], yr[:i], np.zeros(2), V0, a0r, b0r)
    chain += stats.t(2 * an_, X[i] @ mn_, np.sqrt(bn_ / an_ * (1 + X[i] @ Vn_ @ X[i]))).logpdf(yr[i])
assert close(chain, nig_reg_logev(X, yr, np.zeros(2), V0, a0r, b0r), 1e-9)
fig, axs = plt.subplots(1, 3, figsize=(15, 4))
b0g, b1g = np.meshgrid(np.linspace(-1, 5, 200), np.linspace(-1, 2, 200))
Zc = np.stack([b0g, b1g], -1) - mun_r
axs[0].contour(b0g, b1g, np.einsum("...i,ij,...j", Zc, np.linalg.inv(Sn_r), Zc), [2.30, 6.18, 11.8]); axs[0].plot(*mun_r, "k*")
axs[0].set(xlabel="intercept", ylabel="slope", title="posterior contours (68/95/99.7%)")
xg = np.linspace(-2, 12, 200); Xg = np.column_stack([np.ones_like(xg), xg])
for bs in rng.multivariate_normal(mun_r, Sn_r, 50):
    axs[1].plot(xg, Xg @ bs, c="C0", alpha=.15)
axs[1].plot(xr, yr, "ko"); axs[1].set_title("50 posterior-sampled lines")
mg_, vg_ = blr_predict(Xg, mun_r, Sn_r, s2r)
axs[2].fill_between(xg, mg_ - 2 * np.sqrt(vg_), mg_ + 2 * np.sqrt(vg_), alpha=.3); axs[2].plot(xg, mg_); axs[2].plot(xr, yr, "ko"); axs[2].set_title("predictive band widens away from data")
save(fig, "p8_regression.png")
assert vg_[np.argmin(abs(xg - 3))] < vg_[0] and vg_[-1] > vg_[np.argmin(abs(xg - 8))]         # variance grows outward from the data

# %% [markdown]
# ## Part 9: The Laplace approximation for Bayesian logistic regression (notes 3.11; Example 4.8)
# %%
def laplace_logistic(X, y, s0sq, iters=25):
    """MAP by Newton (ST-4 IRLS with a ridge term) and A = -Hessian at the MAP: posterior ~ N(beta_hat, A^-1), eq. (3.34)."""
    b = np.zeros(X.shape[1])
    for k in range(iters + 1):
        p = special.expit(X @ b)
        H = -(X.T * (p * (1 - p))) @ X - np.eye(len(b)) / s0sq
        if k < iters:
            b = b - np.linalg.solve(H, X.T @ (y - p) - b / s0sq)
    return b, -H


def log_joint(X, y, B, s0sq):
    """log p(y|beta) + log N(beta; 0, s0sq I) for each row of B."""
    e = B @ X.T
    return (y * e - np.logaddexp(0, e)).sum(1) + stats.norm.logpdf(B, 0, np.sqrt(s0sq)).sum(1)


x1, y1 = np.array([-2., -1, .5, 1, 2, 3]), np.array([0., 0, 1, 0, 1, 1])      # Example 4.8: one slope, N(0, 2^2) prior
bh1, A1 = laplace_logistic(x1[:, None], y1, 4.0)
bg = np.linspace(-3, 8, 20001); lj = np.exp(log_joint(x1[:, None], y1, bg[:, None], 4.0))
ev1 = np.trapezoid(lj, bg); pm1 = np.trapezoid(bg * lj, bg) / ev1; sd1 = np.sqrt(np.trapezoid((bg - pm1) ** 2 * lj, bg) / ev1)
sd_l = 1 / np.sqrt(A1[0, 0])
logev_lap = log_joint(x1[:, None], y1, bh1[None], 4.0)[0] + 0.5 * np.log(2 * np.pi) - 0.5 * np.log(A1[0, 0])      # (3.35)
print(f"Ex 4.8: MAP {bh1[0]:.4f}, Laplace sd {sd_l:.4f}; exact mean {pm1:.4f}, sd {sd1:.4f}; log evidence Laplace {logev_lap:.4f} vs exact {np.log(ev1):.4f}")
assert abs(logev_lap - np.log(ev1)) < 0.1 and 0.2 < pm1 - bh1[0] < 0.5 and abs(sd_l / sd1 - 1) < 0.15      # evidence good, skew missed
n2 = 60                                                                        # two parameters: brute-force grid vs Laplace
x2 = rng.normal(size=n2); y2 = (rng.random(n2) < special.expit(-0.5 + 1.5 * x2)).astype(float)
X2 = np.column_stack([np.ones(n2), x2]); s0 = 9.0
bh, A = laplace_logistic(X2, y2, s0)
ga, gb = np.linspace(-2.5, 1.5, 300), np.linspace(-0.5, 4.5, 300)
GA, GB = np.meshgrid(ga, gb, indexing="ij")
LJ = log_joint(X2, y2, np.column_stack([GA.ravel(), GB.ravel()]), s0).reshape(GA.shape)
Wg = np.exp(LJ - LJ.max()); Wn = Wg / Wg.sum()
gm = np.array([(Wn * GA).sum(), (Wn * GB).sum()])
c01 = (Wn * (GA - gm[0]) * (GB - gm[1])).sum()
gcov = np.array([[(Wn * (GA - gm[0]) ** 2).sum(), c01], [c01, (Wn * (GB - gm[1]) ** 2).sum()]])
logev_grid = LJ.max() + np.log(Wg.sum() * (ga[1] - ga[0]) * (gb[1] - gb[0]))
logev_lap2 = log_joint(X2, y2, bh[None], s0)[0] + np.log(2 * np.pi) - 0.5 * np.linalg.slogdet(A)[1]
print(f"Part 9: n=60: MAP {bh.round(3)} vs grid mean {gm.round(3)}; sd Laplace {np.sqrt(np.diag(np.linalg.inv(A))).round(3)} vs grid "
      f"{np.sqrt(np.diag(gcov)).round(3)}; log evidence {logev_lap2:.3f} vs {logev_grid:.3f}")
assert np.all(np.abs(bh - gm) < 0.3 * np.sqrt(np.diag(gcov))) and abs(logev_lap2 - logev_grid) < 0.15
assert np.all(np.abs(np.sqrt(np.diag(np.linalg.inv(A))) / np.sqrt(np.diag(gcov)) - 1) < 0.08)
fig, ax = plt.subplots(figsize=(5, 4.5))
lv = [np.exp(-.5 * c) for c in (6.2, 2.3)]
ax.contour(GA, GB, Wg, lv, colors="C0")
Zl = np.stack([GA - bh[0], GB - bh[1]], -1)
ax.contour(GA, GB, np.exp(-.5 * np.einsum("...i,ij,...j", Zl, A, Zl)), lv, colors="C3", linestyles="--")
ax.set(xlabel="intercept", ylabel="slope", title="grid posterior (blue) vs Laplace (red dashed)")
save(fig, "p9_laplace.png")

# %% [markdown]
# ## Part 10: Evidence, Occam factors, BIC (notes 3.12; Example 4.9)
# %%
def gauss_evidence(Phi, y, s2, alpha):
    """log p(y) for y = Phi beta + eps, beta ~ N(0, I/alpha), noise variance s2: eq. (3.38).
    Returns (log evidence, fit at the posterior mean, log Occam factor)."""
    n_, d_ = Phi.shape
    A_ = alpha * np.eye(d_) + Phi.T @ Phi / s2
    m_ = np.linalg.solve(A_, Phi.T @ y / s2)
    fit = -0.5 * n_ * np.log(2 * np.pi * s2) - ((y - Phi @ m_) ** 2).sum() / (2 * s2)
    occam = -0.5 * alpha * m_ @ m_ + 0.5 * d_ * np.log(alpha) - 0.5 * np.linalg.slogdet(A_)[1]
    return fit + occam, fit, occam


xp = np.linspace(-1, 1, 30); sig_n, alpha_p = 0.3, 0.5
yp = 0.5 + 1.0 * xp - 2.0 * xp ** 2 + rng.normal(0, sig_n, xp.size)
rows = []
for dg in range(8):
    Phi = np.vander(xp, dg + 1, increasing=True)
    le, fit, occam = gauss_evidence(Phi, yp, sig_n ** 2, alpha_p)
    assert close(le, stats.multivariate_normal(np.zeros(xp.size), sig_n ** 2 * np.eye(xp.size) + Phi @ Phi.T / alpha_p).logpdf(yp), 1e-7)   # y ~ N(0, s2 I + Phi Phi'/alpha)
    ll_hat = stats.norm.logpdf(yp, Phi @ np.linalg.lstsq(Phi, yp, rcond=None)[0], sig_n).sum()
    rows.append((dg, le, ll_hat, ll_hat - (dg + 1) * np.log(xp.size) / 2, occam, fit))      # column 3 is -BIC/2
R = np.array(rows)
assert R[:, 1].argmax() == 2 and R[:, 3].argmax() == 2                       # evidence and BIC both peak at the true degree
assert np.all(np.diff(R[:, 2]) >= -1e-9) and np.all(np.diff(R[:, 4]) < 0) and np.allclose(R[:, 1], R[:, 4] + R[:, 5])
print("Ex 4.9: degree | log evidence | max log-lik | -BIC/2 | log Occam | fit at posterior mean")
for r_ in R:
    print(f"  d={int(r_[0])}: " + " ".join(f"{v:8.2f}" for v in r_[1:]))
fig, ax = plt.subplots(figsize=(6.5, 4))
for col, lab, mk in [(2, "max log-likelihood", "o-"), (1, "log evidence", "s-"), (3, "-BIC/2", "^--")]:
    ax.plot(R[:, 0], R[:, col], mk, label=lab)
ax.set(xlabel="polynomial degree", title="evidence peaks at the true degree (2)"); ax.legend()
save(fig, "p10_evidence.png")

# %% [markdown]
# ## Part 11: The Jeffreys-Lindley paradox (notes 3.12; Example 4.10)
# %%
def bf10(xbar, n, sigma, tau):
    """BF10 for H0: mu=0 vs H1: mu~N(0,tau^2), xbar~N(mu, sigma^2/n): ratio of marginal densities of xbar, eq. (3.40)."""
    v0 = sigma ** 2 / n
    return np.exp(stats.norm.logpdf(xbar, 0, np.sqrt(v0 + tau ** 2)) - stats.norm.logpdf(xbar, 0, np.sqrt(v0)))


n_j, xb, sg_j, tau_j = 10 ** 4, 0.02, 1.0, 1.0
z = xb * np.sqrt(n_j) / sg_j; pval = 2 * stats.norm.sf(z); u = n_j * tau_j ** 2 / sg_j ** 2
assert close(bf10(xb, n_j, sg_j, tau_j), (1 + u) ** -0.5 * np.exp(0.5 * z * z * u / (1 + u)), 1e-10)
num_marg = integrate.quad(lambda m_: stats.norm.pdf(xb, m_, sg_j / np.sqrt(n_j)) * stats.norm.pdf(m_, 0, tau_j), -8, 8, points=[0, xb])[0]
assert close(num_marg / stats.norm.pdf(xb, 0, sg_j / np.sqrt(n_j)), bf10(xb, n_j, sg_j, tau_j), 1e-6)         # by quadrature
f_u = lambda uu: (1 + uu) ** -0.5 * np.exp(0.5 * z * z * uu / (1 + uu))
best = optimize.minimize_scalar(lambda lu: -f_u(np.exp(lu)), bounds=(-5, 15), method="bounded")
post_sd = np.sqrt(1 / (n_j / sg_j ** 2 + 1 / tau_j ** 2)); post_m = xb * (n_j / sg_j ** 2) * post_sd ** 2
print(f"Ex 4.10: z={z:.2f}, p={pval:.4f}; BF10={bf10(xb, n_j, sg_j, tau_j):.4f} (BF01={1 / bf10(xb, n_j, sg_j, tau_j):.2f}); best-tau BF10="
      f"{-best.fun:.3f} at tau={np.sqrt(np.exp(best.x) * sg_j ** 2 / n_j):.4f}; e^(z^2/2)={np.exp(z * z / 2):.3f}; "
      f"H1 posterior N({post_m:.4f}, {post_sd:.4f}^2), 95% ({post_m - 1.96 * post_sd:.4f}, {post_m + 1.96 * post_sd:.4f})")
assert abs(pval - 0.0455) < 1e-4 and bf10(xb, n_j, sg_j, tau_j) < 1 / 13 and abs(-best.fun - np.exp((z * z - 1) / 2) / z) < 1e-4
assert -best.fun < np.exp(z * z / 2) and post_m - 1.96 * post_sd > 0         # (3.41); estimation says mu != 0, testing says H0
ns_ = np.logspace(1, 8, 100); xs_ = 2 / np.sqrt(ns_)
fig, ax = plt.subplots(figsize=(6.5, 4))
for tau_ in (1, 0.1):
    ax.loglog(ns_, 1 / bf10(xs_, ns_, 1, tau_), label=f"BF01, tau={tau_}")
ax.axhline(1, c="k", lw=.5); ax.set(xlabel="n (z fixed at 2, p = 0.0455 throughout)", ylabel="BF01 (>1 favors H0)", title="Jeffreys-Lindley"); ax.legend()
save(fig, "p11_jl.png")
assert (1 / bf10(xs_, ns_, 1, 1))[-1] > 1e3 and (1 / bf10(xs_, ns_, 1, 1))[0] < 1

# %% [markdown]
# ## Part 12: Posterior predictive checks (notes 3.10; Example 4.11)
# %%
def ppc_pvalue(data, stat, draws=4000):
    """Fit a Normal (reference prior), simulate y_rep from the posterior predictive, return P(T(y_rep) >= T(y) | y), eq. (3.29)."""
    nn_, s_ = data.size, ((data - data.mean()) ** 2).sum()
    s2d = stats.invgamma((nn_ - 1) / 2, scale=s_ / 2).rvs(draws, random_state=rng)
    mud = data.mean() + np.sqrt(s2d / nn_) * rng.standard_normal(draws)
    rep = mud[:, None] + np.sqrt(s2d)[:, None] * rng.standard_normal((draws, nn_))
    T_rep = np.apply_along_axis(stat, 1, rep)
    return np.mean(T_rep >= stat(data)), T_rep


skewed, normal_data = rng.lognormal(0, 0.8, 60), rng.normal(0, 1, 60)
p_sk, T_sk = ppc_pvalue(skewed, stats.skew)
p_no, _ = ppc_pvalue(normal_data, stats.skew)
print(f"Ex 4.11: lognormal data: observed skewness {stats.skew(skewed):.3f}, p_B={p_sk:.4f} (of 4000); Normal data: p_B={p_no:.3f}")
assert p_sk == 0.0 and 0.05 < p_no < 0.95                                    # none of 4000 replicates is as skewed as the data
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.hist(T_sk, 60, alpha=.6); ax.axvline(stats.skew(skewed), c="r"); ax.set_title("replicated skewness vs observed (red)")
save(fig, "p12_ppc.png")

# %% [markdown]
# ## Part 13: The likelihood principle and stopping rules (notes 3.13; Example 4.12)
# %%
p_bin, p_nb = stats.binom.sf(8, 12, .5), stats.nbinom.sf(8, 3, .5)             # P(Y>=9 | n=12) vs P(>=9 heads before the 3rd tail)
pg_ = np.linspace(1e-4, 1 - 1e-4, 4001)
Lb, Ln = stats.binom.pmf(9, 12, pg_), stats.nbinom.pmf(9, 3, 1 - pg_)
pb_, pn_ = Lb / np.trapezoid(Lb, pg_), Ln / np.trapezoid(Ln, pg_)
assert np.allclose(pb_, pn_) and np.allclose(pb_, stats.beta.pdf(pg_, 10, 4), atol=1e-3)   # same likelihood kernel -> same Beta(10,4)
print(f"Ex 4.12: 9 heads, 3 tails: p-value {p_bin:.4f} (n fixed) vs {p_nb:.4f} (stop at third tail); both posteriors Beta(10,4)")

# %% [markdown]
# ## Part 14: Numbers quoted in solutions.md (exercises A1-C3)
# %%
a_, b_ = beta_update(1, 3, 4, 6)                                              # A1
assert (a_, b_) == (5, 5) and close(a_ * (a_ + 1) / ((a_ + b_) * (a_ + b_ + 1)), 30 / 110) and close(stats.beta(a_, b_).std(), 0.15076, 1e-4)
m2, t2 = normal_update(20, 16, 23, 6, 9)                                      # A2
assert close([m2, t2], [22.742857, 1.371429], 1e-6) and close(stats.norm(m2, t2 ** .5).interval(.95), (20.4475, 25.0382), 1e-5)
assert close(stats.norm(m2, (9 + t2) ** .5).interval(.95), (16.4307, 29.0550), 1e-5)
assert close((5 / 6) ** 9, 0.193807, 1e-5) and close(9 * 6 / 25, 2.16)        # A3: Gamma(9,5) posterior, NB predictive
xe = np.array([2.1, 0.7, 3.4]); ge = stats.gamma(3, scale=1 / xe.sum())      # A4: Jeffreys 1/lambda -> Gamma(n, sum x)
assert close(ge.mean(), 3 / xe.sum()) and close(ge.interval(.95), (0.09979, 1.16527), 1e-4)
assert close([1 - .975 ** (1 / 11), 1 - .025 ** (1 / 11), 1 - .05 ** (1 / 11)], [0.002299, 0.284914, 0.238404], 1e-4)   # A5: Beta(1,11) CDF 1-(1-p)^11
assert close(stats.betabinom(10, 9, 5).pmf(10), 0.038248, 1e-4) and close(np.prod([(9 + j) / (14 + j) for j in range(10)]), 0.038248, 1e-4)   # B1
assert close(stats.binom(10, 9 / 14).pmf(10), 0.012054, 1e-3)
hw = stats.t.ppf(.975, 9) * g.std(ddof=1) * np.sqrt(1 + 1 / 10)               # B2
assert close([g.mean() - hw, g.mean() + hw], [9.70710, 9.90490], 1e-5) and close(g.mean() + np.array([-1, 1]) * 1.96 * g.std(ddof=1), [9.72429, 9.88771], 1e-5)
assert close(stats.gamma(5, scale=1 / 10.2).mean(), 0.490196, 1e-5) and close((10.2 / 11.2) ** 5, 0.626485, 1e-5)   # B3
assert close(10.2 / 4, (1 * 4 + 3 * xe.mean()) / 4) and close(stats.invgamma(5, scale=10.2).mean(), 2.55) and close(np.exp(-3 / 6.2), 0.61649, 1e-4)
iv = stats.invgamma(5, scale=8)                                               # B4
assert close([iv.mean(), 8 / 6], [2.0, 1.33333], 1e-5) and close(iv.interval(.95), (0.78113, 4.92767), 1e-4)
a5, b5 = 2 + 12, 1 + 6.5                                                      # B5
assert close([a5 / b5, a5 * (b5 + 1) / b5 ** 2, (b5 / (b5 + 1)) ** a5, np.exp(-12 / 6.5)], [1.86667, 2.11556, 0.173377, 0.157843], 1e-4)
mm = (np.arange(11) + 1) / 12                                                 # B6: Bayes risk of Laplace's rule vs the MLE, uniform prior, n = 10
assert close(np.mean(mm * (1 - mm) / 13), 1 / 72) and close((1 / 6) / 10, 1 / 60)
nA, nB, nC = rng.beta(121, 881, 400_000), rng.beta(146, 856, 400_000), rng.beta(139, 863, 400_000)   # B7 (third arm 138/1000)
assert close([140 / 1200, 165 / 1200, 121 / 1002, 146 / 1002], [0.116667, 0.1375, 0.120758, 0.145709], 1e-4)
print(f"B7: P(best) A,B,C = {(np.bincount(np.argmax([nA, nB, nC], 0), minlength=3) / 400_000).round(3)}; MC P(B>A) = {np.mean(nB > nA):.4f}")
zz = np.array([1.96, 2.576])                                                  # B8
assert close(np.exp((zz ** 2 - 1) / 2) / zz, [2.11247, 6.49901], 1e-4) and close(np.exp(zz ** 2 / 2), [6.82642, 27.60199], 1e-4)
nn_c = 100; ar = np.sqrt(nn_c) / 2                                            # C1: constant risk of the Beta(sqrt(n)/2, sqrt(n)/2) posterior mean
risk_c = [(nn_c * t * (1 - t) + (ar * (1 - 2 * t)) ** 2) / (nn_c + 2 * ar) ** 2 for t in np.linspace(0, 1, 11)]
r_mm = nn_c / (4 * (nn_c + np.sqrt(nn_c)) ** 2)
assert close(risk_c, r_mm) and close(r_mm, 0.0020661, 1e-4) and close(1 / (4 * nn_c), 0.0025)
assert close((1 + np.array([-1, 1]) * np.sqrt(1 - 4 * nn_c * r_mm)) / 2, [0.29170, 0.70830], 1e-4)   # MLE risk p(1-p)/n < r_mm outside this band
u_star = optimize.brentq(lambda uu: np.sqrt(1 + uu) * np.exp(-0.5 * 4 * uu / (1 + uu)) - 1, 5, 1e6)   # C2: BF01 = 1 at u = n tau^2/sigma^2
assert close(u_star, 49.4353, 1e-4) and close(np.sqrt(u_star) / 100, 0.07031, 1e-3) and abs(np.sqrt(1e4) * np.exp(-2) - 1 / bf10(0.02, 10 ** 4, 1, 1)) < 0.05
bm, _ = laplace_logistic(X2, y2, 1e12)                                        # C3: MLE, -BIC/2, and the terms BIC discards
ll_mle = (y2 * (X2 @ bm) - np.logaddexp(0, X2 @ bm)).sum(); neg_bic2 = ll_mle - np.log(n2)
Jm = (X2.T * (special.expit(X2 @ bm) * (1 - special.expit(X2 @ bm)))) @ X2 / n2
dropped = [stats.norm.logpdf(bm, 0, 3).sum(), np.log(2 * np.pi), -0.5 * np.linalg.slogdet(Jm)[1]]
assert abs(neg_bic2 + sum(dropped) - logev_lap2) < 0.05 and abs(neg_bic2 + sum(dropped) - logev_grid) < 0.02
print(f"C3: log evidence grid {logev_grid:.3f}, Laplace {logev_lap2:.3f}, -BIC/2 {neg_bic2:.3f}, max log-lik {ll_mle:.3f}; dropped terms {np.round(dropped, 3)} sum {sum(dropped):.3f}")
print(f"C2: BF01 = 1 at u = {u_star:.2f} (tau=1: n = {u_star:.1f}; tau=0.1: n = {u_star / 0.01:.0f})")

# %% [markdown]
# ## Your turn
#
# 1. A/B test: with `prob_B_beats_A` below, add a third arm and report the posterior probability that each arm is best.
#    How does it change if the prior is Beta(20, 180) instead of uniform? (Exercise B7.)
# 2. Gamma-Poisson with an exposure offset: y_i ~ Poisson(lambda t_i) gives Gamma(alpha0 + sum y, beta0 + sum t) (Exercise B5).
#    Compare the predictive for the next hour with a Poisson plug-in.
# 3. Sensitivity: plot BF01 for Example 4.10 against tau on a log grid. Where is it smallest? (Exercise C2.)
# 4. Extend Part 3 to the Jeffreys prior at theta = 0.01, n = 20: what is the exact coverage? Why?
# 5. Part 10 with sigma unknown: use `nig_reg_logev` from Part 8 and repeat the degree comparison.
#
# ## Reference solution for B7 (write your own first)

# %%
def prob_B_beats_A(yA, nA, yB, nB, prior=(1, 1)):
    """P(theta_B > theta_A | data) for independent Beta posteriors: the integral of f_A(t) (1 - F_B(t)) dt."""
    A_ = stats.beta(prior[0] + yA, prior[1] + nA - yA); B_ = stats.beta(prior[0] + yB, prior[1] + nB - yB)
    return integrate.quad(lambda t: A_.pdf(t) * B_.sf(t), 0, 1)[0]


pq = prob_B_beats_A(120, 1000, 145, 1000)
assert abs(np.mean(rng.beta(146, 856, 400_000) > rng.beta(121, 881, 400_000)) - pq) < 0.003
print(f"B7 reference: P(B > A) = {pq:.4f} (uniform priors); with Beta(20,180) priors {prob_B_beats_A(120, 1000, 145, 1000, (20, 180)):.4f}")

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
