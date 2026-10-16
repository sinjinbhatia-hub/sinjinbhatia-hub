# %% [markdown]
# # ST-4 Project: GLMs from scratch (IRLS, deviance, logistic, Poisson, overdispersion)
#
# One function, `glm_fit`, fits every model in the notes by iteratively reweighted least squares
# (notes (3.9)-(3.10)). We check it against direct maximum likelihood (scipy.optimize on scipy.stats
# log-densities) and against OLS, print the quadratic convergence of Newton = Fisher scoring, then
# run the logistic, softmax, Poisson-offset and overdispersion experiments of notes sections 3.5-3.9.
# Every number quoted in notes.md, exercises.md and solutions.md is asserted or printed here.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy import integrate, optimize, stats
from scipy.special import digamma, expit, gammaln, xlogy

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
Z = stats.norm.ppf(0.975)

# %% [markdown]
# ## Part 0: The GLM engine
# Links are tuples (g, g^{-1}, g', g''); families carry V(mu), V'(mu) and the unit deviance (3.18).

# %%
EPS = np.finfo(float).eps
_clip = lambda m: np.clip(m, EPS, 1 - EPS)        # as R's binomial(): fitted probabilities never reach exactly 0 or 1
LINKS = {
    "identity": (lambda m: m, lambda e: e, lambda m: np.ones_like(m), lambda m: np.zeros_like(m)),
    "log": (np.log, np.exp, lambda m: 1 / m, lambda m: -1 / m ** 2),
    "logit": (lambda m: np.log(m / (1 - m)), lambda e: _clip(expit(e)), lambda m: 1 / (m * (1 - m)), lambda m: (2 * m - 1) / (m * (1 - m)) ** 2),
    "probit": (stats.norm.ppf, lambda e: _clip(stats.norm.cdf(e)), lambda m: 1 / stats.norm.pdf(stats.norm.ppf(m)),
               lambda m: stats.norm.ppf(m) / stats.norm.pdf(stats.norm.ppf(m)) ** 2),
    "cloglog": (lambda m: np.log(-np.log1p(-m)), lambda e: _clip(-np.expm1(-np.exp(e))),
                lambda m: -1 / ((1 - m) * np.log1p(-m)), lambda m: -(np.log1p(-m) + 1) / ((1 - m) * np.log1p(-m)) ** 2),
    "inverse": (lambda m: 1 / m, lambda e: 1 / e, lambda m: -1 / m ** 2, lambda m: 2 / m ** 3),
}

def _gamma_phi_ml(y, mu, om, dev):
    """ML dispersion for the Gamma: log(nu) - digamma(nu) = D/(2n), phi = 1/nu (Exercise C2)."""
    c = dev / (2 * len(y))
    return 1 / optimize.brentq(lambda v: np.log(v) - digamma(v) - c, 1e-8, 1e12)

FAMILIES = {
    "normal": dict(V=lambda m: np.ones_like(m), dV=lambda m: np.zeros_like(m), dev=lambda y, m: (y - m) ** 2,
                   default="identity", phi_known=False, init=lambda y, om: y.astype(float),
                   ll=lambda y, m, om, ph: stats.norm.logpdf(y, m, np.sqrt(ph / om)).sum(),
                   phi_ml=lambda y, m, om, D: D / len(y)),
    "bernoulli": dict(V=lambda m: m * (1 - m), dV=lambda m: 1 - 2 * m,
                      dev=lambda y, m: 2 * (xlogy(y, y / m) + xlogy(1 - y, (1 - y) / (1 - m))),
                      default="logit", phi_known=True, init=lambda y, om: (om * y + 0.5) / (om + 1),
                      ll=lambda y, m, om, ph: stats.binom.logpmf(np.rint(om * y), om, m).sum()),
    "poisson": dict(V=lambda m: m, dV=lambda m: np.ones_like(m), dev=lambda y, m: 2 * (xlogy(y, y / m) - (y - m)),
                    default="log", phi_known=True, init=lambda y, om: y + 0.5,
                    ll=lambda y, m, om, ph: stats.poisson.logpmf(y, m).sum()),
    "gamma": dict(V=lambda m: m ** 2, dV=lambda m: 2 * m, dev=lambda y, m: 2 * (-np.log(y / m) + (y - m) / m),
                  default="inverse", phi_known=False, init=lambda y, om: y.astype(float),
                  ll=lambda y, m, om, ph: stats.gamma.logpdf(y, a=1 / ph, scale=m * ph).sum(), phi_ml=_gamma_phi_ml),
}

def negbin(r):
    """Negative Binomial with fixed size r: V = mu + mu^2/r (3.31); a GLM for each fixed r."""
    return dict(V=lambda m: m + m ** 2 / r, dV=lambda m: 1 + 2 * m / r,
                dev=lambda y, m: 2 * (xlogy(y, y / m) - (y + r) * np.log((y + r) / (m + r))),
                default="log", phi_known=True, init=lambda y, om: y + 0.5,
                ll=lambda y, m, om, ph: stats.nbinom.logpmf(y, r, r / (r + m)).sum())

def glm_fit(X, y, family="poisson", link=None, offset=None, weights=None, ridge=0.0, observed=False,
            maxit=100, tol=1e-10, beta0=None, est_phi=None, with_null=True):
    """IRLS for g(mu) = X beta + offset. Each pass solves the weighted least squares problem (3.9):
    beta <- (X'WX + ridge*I0)^{-1} X'W z with w_i = om_i/(V g'^2) and z_i = eta_i - o_i + (y_i-mu_i) g'(mu_i) (3.10).
    observed=True takes Newton steps with the observed-information weights (3.11) instead of Fisher scoring.
    The first column is the intercept and is never penalized. Returns a dict with the inference table (3.14)-(3.16)."""
    fam = family if isinstance(family, dict) else FAMILIES[family]
    g, ginv, gp, gpp = LINKS[link or fam["default"]]
    y = np.asarray(y, float); n, p = X.shape
    off = np.zeros(n) if offset is None else offset
    om = np.ones(n) if weights is None else weights
    pen = ridge * np.eye(p); pen[0, 0] = 0.0
    if beta0 is None:
        mu = fam["init"](y, om); eta = g(mu)
    else:
        eta = X @ beta0 + off; mu = ginv(eta)
    obj_old, path, beta = np.inf, [], (None if beta0 is None else np.asarray(beta0, float))
    for it in range(1, maxit + 1):
        gprime, V = gp(mu), fam["V"](mu)
        w = om / (V * gprime ** 2)                                     # expected (Fisher) weights (3.9)
        z = eta - off + (y - mu) * gprime                              # working response (3.10)
        if observed and beta is not None:                              # Newton: beta + (X'W_obs X)^{-1} u, W_obs from (3.11)
            w_obs = w + om * (y - mu) * (fam["dV"](mu) * gprime + V * gpp(mu)) / (V ** 2 * gprime ** 3)
            A = X.T @ (w_obs[:, None] * X) + pen
            beta = beta + np.linalg.solve(A, X.T @ (w * gprime * (y - mu)) - pen @ beta)
        else:                                                          # Fisher scoring = IRLS
            beta = np.linalg.solve(X.T @ (w[:, None] * X) + pen, X.T @ (w * z))
        eta = X @ beta + off; mu = ginv(eta)
        dev = np.sum(om * fam["dev"](y, mu)); obj = dev + ridge * beta[1:] @ beta[1:]
        path.append(beta)
        if abs(obj_old - obj) < tol * (abs(obj) + 0.1):
            break
        obj_old = obj
    gprime, V = gp(mu), fam["V"](mu)
    w = om / (V * gprime ** 2)
    info_inv = np.linalg.inv(X.T @ (w[:, None] * X) + pen)
    pearson = np.sum(om * (y - mu) ** 2 / V)
    estimated = (not fam["phi_known"]) if est_phi is None else est_phi
    phi = pearson / (n - p) if estimated else 1.0
    se = np.sqrt(phi * np.diag(info_inv)); zst = beta / se
    q, pv = (stats.t.ppf(.975, n - p), 2 * stats.t.sf(abs(zst), n - p)) if estimated else (Z, 2 * stats.norm.sf(abs(zst)))
    phi_ml = fam["phi_ml"](y, mu, om, dev) if "phi_ml" in fam else 1.0
    ll = fam["ll"](y, mu, om, phi_ml)
    out = dict(beta=beta, se=se, z=zst, p=pv, ci=np.c_[beta - q * se, beta + q * se], mu=mu, eta=eta, dev=dev, phi=phi,
               w=w, h=w * np.einsum("ij,jk,ik->i", X, info_inv, X), pearson=pearson, cov=phi * info_inv, iters=it,
               path=np.array(path), loglik=ll, aic=-2 * ll + 2 * (p + (1 if "phi_ml" in fam else 0)), df_res=n - p,
               converged=it < maxit)
    out["null_dev"] = glm_fit(X[:, :1], y, family, link, offset, weights, with_null=False)["dev"] if with_null and np.all(X[:, 0] == 1) else np.nan
    return out

def table(fit, names):
    """Print the inference table: coefficient, se, Wald z, p, 95% CI."""
    print("  term        coef      se        z        p      95% CI")
    for k, nm in enumerate(names):
        print(f"  {nm:8s} {fit['beta'][k]:8.4f} {fit['se'][k]:8.4f} {fit['z'][k]:8.3f} {fit['p'][k]:8.4f}  [{fit['ci'][k,0]:.4f}, {fit['ci'][k,1]:.4f}]")
    print(f"  deviance {fit['dev']:.3f} on {fit['df_res']} df (null {fit['null_dev']:.3f}); phi = {fit['phi']:.3f}; AIC = {fit['aic']:.2f}; {fit['iters']} IRLS iterations")

# %% [markdown]
# ## Part 1: Why linear regression fails for binary and count data (notes section 1)

# %%
rng = np.random.default_rng(1)
n = 400
x = rng.uniform(-3, 3, n); Xs = np.column_stack([np.ones(n), x])
yb = rng.binomial(1, expit(1.2 * x))
b_ols = np.linalg.lstsq(Xs, yb, rcond=None)[0]; fit_l = Xs @ b_ols; res_l = yb - fit_l
frac_out = np.mean((fit_l < 0) | (fit_l > 1))
ends, mid = np.abs(x) > 2, np.abs(x) < 0.6
v_ends, v_mid = res_l[ends].var(), res_l[mid].var()
sl_ = x > 1.5; sk, sw = stats.skew(res_l[sl_]), stats.shapiro(res_l[sl_]).pvalue
print(f"Binary: OLS line b = {np.round(b_ols, 3)}; {100*frac_out:.1f}% of fitted values outside [0,1]; residual variance "
      f"{v_ends:.3f} for |x|>2 vs {v_mid:.3f} for |x|<0.6; slice x>1.5 ({sl_.sum()} points): skewness {sk:.2f}, Shapiro p = {sw:.1e}")
assert frac_out > 0.1 and v_mid > 5 * v_ends and sk < -1.5 and sw < 1e-6
xc = rng.uniform(-3, 3, n); Xc = np.column_stack([np.ones(n), xc]); yc = rng.poisson(np.exp(0.3 + 0.6 * xc))
b_c = np.linalg.lstsq(Xc, yc, rcond=None)[0]; fit_c = Xc @ b_c; res_c = yc - fit_c
lo3, hi3 = xc < -1, xc > 1
print(f"Counts: {100*np.mean(fit_c < 0):.1f}% of OLS fitted means are negative; residual variance {res_c[lo3].var():.2f} (x<-1) vs {res_c[hi3].var():.2f} (x>1)")
assert np.mean(fit_c < 0) > 0.05 and res_c[hi3].var() > 10 * res_c[lo3].var()
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
ax[0].scatter(x, yb, s=6, alpha=.4); gx = np.linspace(-3, 3, 100)
ax[0].plot(gx, b_ols[0] + b_ols[1] * gx, "r", label="OLS"); ax[0].plot(gx, expit(1.2 * gx), "k", label="true p(x)")
ax[0].set_ylim(-.4, 1.4); ax[0].legend(); ax[0].set_title("Binary response: line leaves [0, 1]")
ax[1].scatter(xc, yc, s=6, alpha=.4); ax[1].plot(gx, b_c[0] + b_c[1] * gx, "r"); ax[1].plot(gx, np.exp(.3 + .6 * gx), "k")
ax[1].set_title("Counts: negative fitted means, variance grows with the mean")
fig.tight_layout(); fig.savefig(FIG / "ols_fails.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 2: The exponential-dispersion-family table (notes section 3.1, Example 4.1)
# SymPy derives mu = b'(theta), V(mu) = b''(theta(mu)), the canonical link, and the unit deviance 2 int_mu^y (y-t)/V(t) dt.

# %%
th = sp.symbols("theta", real=True); mu_s, y_s, t_s = sp.symbols("mu y t", positive=True)
bfun = {"Normal": (th ** 2 / 2, mu_s, sp.Integer(1)), "Bernoulli": (sp.log(1 + sp.exp(th)), sp.log(mu_s / (1 - mu_s)), mu_s * (1 - mu_s)),
        "Poisson": (sp.exp(th), sp.log(mu_s), mu_s), "Gamma": (-sp.log(-th), -1 / mu_s, mu_s ** 2),
        "InvGauss": (-sp.sqrt(-2 * th), -1 / (2 * mu_s ** 2), mu_s ** 3)}
for name, (b, theta_of_mu, V) in bfun.items():
    assert sp.simplify(sp.diff(b, th).subs(th, theta_of_mu) - mu_s) == 0, name          # b'(theta(mu)) = mu
    assert sp.simplify(sp.diff(b, th, 2).subs(th, theta_of_mu) - V) == 0, name          # b''(theta(mu)) = V(mu)
    d_sym = sp.simplify(2 * sp.integrate((y_s - t_s) / V.subs(mu_s, t_s), (t_s, mu_s, y_s)))   # (3.19)
    print(f"{name:10s} b = {b},  theta(mu) = {theta_of_mu},  V(mu) = {V},  unit deviance = {d_sym}")
# Theorem 3.1 numerically: mean = b'(theta), Var = phi b''(theta), via scipy.stats distributions.
mu0, lam, nu = 3.7, 2.5, 4.0
checks = [(stats.poisson(mu0), mu0, mu0), (stats.bernoulli(0.3), 0.3, 0.3 * 0.7), (stats.gamma(a=nu, scale=mu0 / nu), mu0, mu0 ** 2 / nu),
          (stats.invgauss(mu0 / lam, scale=lam), mu0, mu0 ** 3 / lam), (stats.binom(10, .3), 3.0, 10 * .3 * .7)]
for dist, m_, v_ in checks:
    assert np.isclose(dist.mean(), m_) and np.isclose(dist.var(), v_)
print("Theorem 3.1 checked against scipy.stats: Poisson, Bernoulli, Gamma (phi=1/nu), Inverse Gaussian (phi=1/lambda), Binomial (a=1/m).")
p_ = 0.3; th_b = np.log(p_ / (1 - p_)); b_b = np.log1p(np.exp(th_b))              # Example 4.1
assert np.isclose(th_b, -0.8473, atol=5e-5) and np.isclose(b_b, 0.3567, atol=5e-5) and np.isclose(b_b, -np.log(0.7)) and np.isclose(np.log(10 / 7), b_b)
assert np.isclose(expit(th_b), 0.3) and np.isclose(expit(th_b) * (1 - expit(th_b)), 0.21) and np.isclose(0.21 / 10, 0.021)
assert np.isclose(np.log(3.7), 1.3083, atol=5e-5)
print(f"Ex 4.1: Bernoulli p=0.3: theta = {th_b:.4f}, b = {b_b:.4f}, V = 0.21; Poisson mu=3.7: theta = {np.log(3.7):.4f}")
# Canonical link = (b')^{-1}: logit, log, identity, -1/mu (inverse), mu^{-2}.  Poisson deviance by quadrature equals the closed form:
yv, mv = 7.0, 4.2
quad = 2 * integrate.quad(lambda t: (yv - t) / t, mv, yv)[0]
assert np.isclose(quad, FAMILIES["poisson"]["dev"](np.array(yv), np.array(mv)))
for fn, (yy, mm) in (("bernoulli", (1.0, .3)), ("gamma", (7., 4.2)), ("normal", (7., 4.2))):
    Vf = FAMILIES[fn]["V"]
    assert np.isclose(2 * integrate.quad(lambda t: (yy - t) / Vf(np.array(t)), mm, yy)[0], FAMILIES[fn]["dev"](np.array(yy), np.array(mm)))
print("Unit deviances = 2 int (y-t)/V(t) dt for Normal, Bernoulli, Poisson, Gamma (3.19).")

# %% [markdown]
# ## Part 3: IRLS versus direct maximum likelihood, OLS, and the moment equations (notes 3.3-3.4)

# %%
rng = np.random.default_rng(3)
n = 300
X = np.column_stack([np.ones(n), rng.normal(size=(n, 2))]); bt = np.array([0.4, 0.7, -0.5])
eta_t = X @ bt
data = {("normal", "identity"): eta_t + rng.normal(0, 1.5, n),
        ("bernoulli", "logit"): rng.binomial(1, expit(eta_t)).astype(float),
        ("bernoulli", "probit"): rng.binomial(1, stats.norm.cdf(eta_t)).astype(float),
        ("poisson", "log"): rng.poisson(np.exp(0.3 * eta_t)).astype(float),
        ("gamma", "log"): rng.gamma(4.0, np.exp(0.3 * eta_t) / 4.0)}

def direct_ml(X, y, family, link):
    """Maximize the scipy.stats log-likelihood over beta with BFGS (phi fixed: beta-hat does not depend on it)."""
    ginv = LINKS[link][1]
    def nll(b):
        m = ginv(X @ b)
        return -{"normal": lambda: stats.norm.logpdf(y, m, 1).sum(), "bernoulli": lambda: stats.bernoulli.logpmf(y, np.clip(m, 1e-12, 1 - 1e-12)).sum(),
                 "poisson": lambda: stats.poisson.logpmf(y, m).sum(), "gamma": lambda: stats.gamma.logpdf(y, a=4, scale=m / 4).sum()}[family]()
    r = optimize.minimize(nll, np.zeros(X.shape[1]), method="BFGS", options={"gtol": 1e-9})
    return r.x, -r.fun

for (fam_, lnk_), y_ in data.items():
    f0 = glm_fit(X, y_, fam_, lnk_); f = glm_fit(X, y_, fam_, lnk_, tol=0, maxit=40)   # f: run to machine precision for the score checks
    b_dir, _ = direct_ml(X, y_, fam_, lnk_)
    assert np.max(np.abs(f["beta"] - b_dir)) < 1e-4, (fam_, lnk_, f["beta"], b_dir)
    canon = LINKS[lnk_] is LINKS[FAMILIES[fam_]["default"]]
    score_x = np.abs(X.T @ (y_ - f["mu"])).max()
    print(f"{fam_:9s}/{lnk_:8s}: IRLS {f0['iters']:2d} its, max|beta_IRLS - beta_BFGS| = {np.max(np.abs(f['beta'] - b_dir)):.1e}; max|X'(y-mu)| = {score_x:.1e}")
    if canon:
        assert score_x < 1e-6     # canonical link: X'(y - mu) = 0 (3.7)
    else:
        gpv = LINKS[lnk_][2](f["mu"])
        assert np.abs(X.T @ (f["w"] * gpv * (y_ - f["mu"]))).max() < 1e-6 and score_x > 1e-3   # the real score (3.5) vanishes, X'(y-mu) does not
# Normal/identity versus ST-3's OLS: one IRLS pass, same beta, same se
yn = data[("normal", "identity")]; fn_ = glm_fit(X, yn, "normal", "identity")
b_o = np.linalg.lstsq(X, yn, rcond=None)[0]; s2 = np.sum((yn - X @ b_o) ** 2) / (n - 3)
assert np.allclose(fn_["beta"], b_o) and np.allclose(fn_["se"], np.sqrt(s2 * np.diag(np.linalg.inv(X.T @ X)))) and fn_["iters"] <= 2
assert np.isclose(fn_["dev"], np.sum((yn - X @ b_o) ** 2)) and np.isclose(fn_["phi"], s2)
print(f"Normal/identity: IRLS = OLS in {fn_['iters']} passes; deviance = RSS = {fn_['dev']:.3f}, phi-hat = RSS/(n-p) = {fn_['phi']:.4f}")
# Independent log-likelihood check of the engine's own loglik() for Poisson and Bernoulli
fp = glm_fit(X, data[("poisson", "log")], "poisson")
assert np.isclose(fp["loglik"], stats.poisson.logpmf(data[("poisson", "log")], fp["mu"]).sum())

# Quadratic convergence of Newton = Fisher scoring for the canonical logit link (notes 3.4)
yl = data[("bernoulli", "logit")]
f_ref = glm_fit(X, yl, "bernoulli", tol=1e-15, maxit=30)
trace = glm_fit(X, yl, "bernoulli", beta0=np.zeros(3), tol=0, maxit=8)["path"]
err = np.linalg.norm(trace - f_ref["beta"], axis=1)
dev_tr = [np.sum(FAMILIES["bernoulli"]["dev"](yl, expit(X @ b))) for b in trace]
print("IRLS trace (logit, n=300), start beta=0:")
for t, (e_, d_) in enumerate(zip(err, dev_tr), 1):
    print(f"  iter {t}: ||beta_t - beta_hat|| = {e_:.2e}   deviance = {d_:.6f}")
for e0, e1 in zip(err[:-1], err[1:]):
    if 1e-6 < e0 < 0.3:
        assert e1 < 5 * e0 ** 2         # quadratic: e_{t+1} <= C e_t^2
# non-canonical link (probit): Fisher scoring is not Newton; it still converges, just not quadratically in general
# Gamma/log has constant weights W = I (Exercise B1)
fg = glm_fit(X, data[("gamma", "log")], "gamma", "log")
assert np.allclose(fg["w"], 1.0)
fg1 = glm_fit(np.ones((3, 1)), np.array([1., 2., 6.]), "gamma", "log", beta0=np.zeros(1), maxit=1)
assert np.isclose(fg1["beta"][0], 2.0)

# %% [markdown]
# ## Part 4: Worked example, six data points, one IRLS step by hand (Example 4.2)

# %%
x6 = np.array([-2, -1, 0, 1, 2, 3.]); y6 = np.array([0, 0, 1, 0, 1, 1.]); X6 = np.column_stack([np.ones(6), x6])
mu0 = np.full(6, .5)
u0 = X6.T @ (y6 - mu0); W0 = np.diag(np.full(6, .25)); z0 = (y6 - mu0) / (mu0 * (1 - mu0))
assert np.allclose(u0, [0, 3.5]) and np.allclose(z0, 4 * y6 - 2) and np.allclose(X6.T @ X6, [[6, 3], [3, 19]])
b1 = np.linalg.solve(X6.T @ W0 @ X6, X6.T @ W0 @ z0)
assert np.allclose(b1, [-0.4, 0.8]) and np.isclose(np.linalg.det(X6.T @ X6), 105)
assert np.allclose(b1, 4 * np.linalg.solve(X6.T @ X6, u0))                 # Fisher scoring form beta + I^{-1} u
f6 = glm_fit(X6, y6, "bernoulli", tol=1e-14)
p6_1 = expit(X6 @ b1); w1 = p6_1 * (1 - p6_1)
b2 = np.linalg.solve(X6.T @ (w1[:, None] * X6), X6.T @ (w1 * (X6 @ b1 + (y6 - p6_1) / w1)))
print(f"Ex 4.2: score at 0 = {u0}, W = I/4, z = {z0}, beta1 = {b1}, beta2 = {np.round(b2, 4)}, MLE = {np.round(f6['beta'], 4)} "
      f"(se {np.round(f6['se'], 4)}), {f6['iters']} its; deviance {f6['dev']:.4f}, null deviance {f6['null_dev']:.4f}")
assert np.allclose(glm_fit(X6, y6, "bernoulli", beta0=b1, maxit=1)["beta"], b2)
lrt6 = f6["null_dev"] - f6["dev"]; wz6 = f6["beta"][1] / f6["se"][1]
f6s = glm_fit(X6[:, :1], y6, "bernoulli")       # intercept-only: p-hat = 1/2
assert np.allclose(f6s["mu"], .5) and np.isclose(lrt6, 3.3618, atol=5e-5) and np.isclose(wz6, 1.330, atol=5e-4)
print(f"   LRT for slope {lrt6:.4f} (p = {stats.chi2.sf(lrt6, 1):.4f}); Wald z = {wz6:.3f} (p = {2*stats.norm.sf(wz6):.4f})")
assert np.isclose(f6["null_dev"], 2 * 6 * np.log(2)) and np.allclose(X6.T @ (y6 - f6["mu"]), 0, atol=1e-8)

# %% [markdown]
# ## Part 5: Logistic regression with a known beta: inference, Wald coverage, LRT (notes 3.5, 3.7)
# Story for Example 4.3: y = crack detected, x1 = standardized stress amplitude, x2 = coated surface (synthetic).

# %%
rng = np.random.default_rng(5)
beta_true = np.array([-0.5, 1.0, 0.7])

def sim_logit(n, rng):
    Xl = np.column_stack([np.ones(n), rng.normal(size=n), rng.binomial(1, .5, n)])
    return Xl, rng.binomial(1, expit(Xl @ beta_true)).astype(float)

Xl, yl = sim_logit(200, rng)
fl = glm_fit(Xl, yl, "bernoulli"); table(fl, ["const", "x1", "x2"])
assert np.allclose(fl["cov"], np.linalg.inv(Xl.T @ (fl["w"][:, None] * Xl)))   # se from (X'WX)^{-1}, (3.14)
assert np.isclose(fl["aic"], fl["dev"] + 2 * 3)                                # AIC = D + 2p when ell(saturated) = 0 (Bernoulli)
f_sub = glm_fit(Xl[:, :2], yl, "bernoulli"); lrt = f_sub["dev"] - fl["dev"]
print(f"LRT for x2: D(small) - D(full) = {lrt:.3f}, p = {stats.chi2.sf(lrt, 1):.4f}; Wald z^2 = {fl['z'][2]**2:.3f}")
R, cover = 1000, np.zeros(3)
for _ in range(R):
    Xr, yr = sim_logit(200, rng); fr = glm_fit(Xr, yr, "bernoulli", with_null=False)
    cover += (fr["ci"][:, 0] <= beta_true) & (beta_true <= fr["ci"][:, 1])
cover /= R
print(f"Wald 95% coverage over {R} replications (n=200): {np.round(cover, 3)}")
assert np.all(np.abs(cover - .95) < .025)
print(f"Ex 4.3: x1 OR = {np.exp(fl['beta'][1]):.4f}, 95% CI [{np.exp(fl['ci'][1, 0]):.4f}, {np.exp(fl['ci'][1, 1]):.4f}]")
assert np.allclose(np.round([np.exp(fl["beta"][1]), *np.exp(fl["ci"][1])], 2), [3.34, 2.20, 5.08])
# odds ratio versus relative risk (Example 4.3 and 3.7): RR = OR / (1 - p0 + p0 OR)
b_or = 0.7; OR = np.exp(b_or); print(f"Ex 4.3: OR = e^0.7 = {OR:.4f}")
for p0 in (.01, .1, .5):
    p1 = OR * p0 / (1 - p0) / (1 + OR * p0 / (1 - p0)); rr = p1 / p0
    assert np.isclose(rr, OR / (1 - p0 + p0 * OR))
    print(f"   baseline {p0}: new probability {p1:.4f}, RR = {rr:.4f}")
# latent variable: P(Y=1) = F(x'beta); slope matching gives 4 phi(0) = 1.596, sd matching pi/sqrt(3) = 1.814, LS fit ~ 1.70
c_ls = optimize.minimize_scalar(lambda c: integrate.quad(lambda t: (stats.norm.cdf(t) - expit(c * t)) ** 2, -8, 8)[0], bounds=(1, 2.5), method="bounded").x
print(f"logit/probit scale: slope match 4*phi(0) = {4*stats.norm.pdf(0):.4f}, sd match pi/sqrt(3) = {np.pi/np.sqrt(3):.4f}, least squares c = {c_ls:.4f}")
assert 1.69 < c_ls < 1.71
n_big = 40000; xb = rng.normal(size=n_big); Xb = np.column_stack([np.ones(n_big), xb])
ylat = (0.3 + 1.0 * xb + rng.normal(size=n_big) > 0).astype(float)
fprob = glm_fit(Xb, ylat, "bernoulli", "probit", with_null=False); flogit = glm_fit(Xb, ylat, "bernoulli", "logit", with_null=False)
ratio = flogit["beta"][1] / fprob["beta"][1]
print(f"Latent-normal data: probit slope {fprob['beta'][1]:.3f}, logit slope {flogit['beta'][1]:.3f}, ratio {ratio:.3f}")
assert 1.55 < ratio < 1.85 and abs(fprob["beta"][1] - 1.0) < 0.05

# %% [markdown]
# ## Part 6: ROC, AUC, calibration, Hosmer-Lemeshow, decision boundary (notes 3.7)

# %%
def roc_curve(y, s):
    """Empirical ROC: sweep a threshold down the sorted scores; ties get one step. Returns FPR, TPR, trapezoid AUC."""
    order = np.argsort(-s); s_, y_ = s[order], y[order]
    last = np.r_[np.where(np.diff(s_))[0], len(s_) - 1]
    tpr = np.r_[0, np.cumsum(y_)[last] / y.sum()]; fpr = np.r_[0, np.cumsum(1 - y_)[last] / (1 - y).sum()]
    return fpr, tpr, np.trapezoid(tpr, fpr)

def auc_pairs(y, s):
    """P(S1 > S0) + 0.5 P(S1 = S0) over all positive-negative pairs ((3.27)."""
    d = s[y == 1][:, None] - s[y == 0][None, :]
    return np.mean(d > 0) + 0.5 * np.mean(d == 0)

rng = np.random.default_rng(6)
n2 = 600; X2 = np.column_stack([np.ones(n2), rng.normal(size=(n2, 2))]); b2t = np.array([-0.3, 1.5, -1.0])
y2 = rng.binomial(1, expit(X2 @ b2t)).astype(float); f2 = glm_fit(X2, y2, "bernoulli"); p2 = f2["mu"]
fpr, tpr, auc = roc_curve(y2, p2); auc_p = auc_pairs(y2, p2)
u_stat = stats.mannwhitneyu(p2[y2 == 1], p2[y2 == 0]).statistic / (y2.sum() * (1 - y2).sum())
assert np.isclose(auc, auc_p) and np.isclose(auc, u_stat)
s_tied = np.round(p2, 1); assert np.isclose(roc_curve(y2, s_tied)[2], auc_pairs(y2, s_tied))     # ties count 1/2
assert np.isclose(roc_curve(y2, p2 ** 3)[2], auc)                                                 # monotone transform: AUC unchanged
print(f"AUC = {auc:.4f} (trapezoid) = {auc_p:.4f} (pairs) = {u_stat:.4f} (Mann-Whitney); with scores rounded to 0.1: {auc_pairs(y2, s_tied):.4f}")
ss8 = np.array([.8, .5, .4, .6, .3, .1]); yy8 = np.r_[np.ones(3), np.zeros(3)]       # Example 4.8
assert np.isclose(auc_pairs(yy8, ss8), 7 / 9) and np.isclose(roc_curve(yy8, ss8)[2], 7 / 9)
s1, s0 = np.array([.9, .6, .6, .3]), np.array([.7, .6, .2, .1])
yy = np.r_[np.ones(4), np.zeros(4)]; ss = np.r_[s1, s0]
assert np.isclose(auc_pairs(yy, ss), 11 / 16) and np.isclose(roc_curve(yy, ss)[2], 11 / 16)         # Exercise B2

def hosmer_lemeshow(y, p, g=10):
    """C = sum_g (O_g - E_g)^2 / (n_g pbar_g (1 - pbar_g)), ~ chi^2_{g-2} (heuristic)."""
    idx = np.array_split(np.argsort(p), g); C = 0
    for i in idx:
        E, nn = p[i].sum(), len(i); C += (y[i].sum() - E) ** 2 / (E * (1 - E / nn))
    return C
C_good, C_bad = hosmer_lemeshow(y2, p2), hosmer_lemeshow(y2, p2 ** 3)
print(f"Hosmer-Lemeshow: model {C_good:.2f} (p = {stats.chi2.sf(C_good, 8):.3f}); same ranking after p -> p^3: {C_bad:.1f} (p = {stats.chi2.sf(C_bad, 8):.1e})")
assert stats.chi2.sf(C_good, 8) > 0.05 and stats.chi2.sf(C_bad, 8) < 1e-4
fig, ax = plt.subplots(1, 3, figsize=(14, 4))
ax[0].plot(fpr, tpr); ax[0].plot([0, 1], [0, 1], "k:"); ax[0].set(title=f"ROC, AUC = {auc:.3f}", xlabel="FPR", ylabel="TPR")
for pp, lab in ((p2, "fitted"), (p2 ** 3, "p^3 (same AUC)")):
    ids = np.array_split(np.argsort(pp), 10); ax[1].plot([pp[i].mean() for i in ids], [y2[i].mean() for i in ids], "o-", label=lab)
ax[1].plot([0, 1], [0, 1], "k:"); ax[1].legend(); ax[1].set(title="Calibration (deciles)", xlabel="predicted", ylabel="observed")
gx1, gx2 = np.meshgrid(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100))
ax[2].contourf(gx1, gx2, expit(f2["beta"][0] + f2["beta"][1] * gx1 + f2["beta"][2] * gx2), 20, cmap="RdBu_r", alpha=.6)
ax[2].contour(gx1, gx2, f2["beta"][0] + f2["beta"][1] * gx1 + f2["beta"][2] * gx2, [0], colors="k")
ax[2].scatter(X2[:, 1], X2[:, 2], c=y2, s=6, cmap="coolwarm"); ax[2].set_title("Decision boundary x'beta = 0 (a line)")
fig.tight_layout(); fig.savefig(FIG / "logistic_diagnostics.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 7: Complete separation and the ridge fix (notes 3.4, Example 4.4)

# %%
xs = np.arange(1., 7.); ys = np.array([0, 0, 0, 1, 1, 1.]); Xsep = np.column_stack([np.ones(6), xs])
fsep = glm_fit(Xsep, ys, "bernoulli", beta0=np.zeros(2), tol=0, maxit=25, with_null=False)
dev_sep = lambda b: np.sum(FAMILIES["bernoulli"]["dev"](ys, LINKS["logit"][1](Xsep @ b)))
dv = np.array([dev_sep(b) for b in fsep["path"]]); sl = fsep["path"][:, 1]
print("Separation (x=1..6, y=000111): iteration, intercept, slope, deviance")
for t in (1, 2, 3, 5, 10, 15, 20, 25):
    b_ = fsep["path"][t - 1]; print(f"  {t:2d}  {b_[0]:9.3f} {b_[1]:8.4f}   deviance = {dv[t-1]:.3e}")
print(f"  after 25 passes: slope {sl[-1]:.2f}, se {fsep['se'][1]:.3g}, Wald z = {fsep['z'][1]:.4f}")
assert fsep["se"][1] > 1e4 and abs(fsep["z"][1]) < 0.01                       # the Wald test sees nothing (Hauck-Donner)
# working residual z - eta -> +-1 for the two points next to the boundary, so the slope gains 2/(4-3) = 2 per step
assert np.all(np.diff(sl) > 0) and sl[-1] > 40 and dv[-1] < 1e-6 and np.all(np.abs(np.diff(sl[8:]) - 2) < 0.05)
for c in (1, 10, 100):                                                       # loglik(c beta) increases in c for a separating beta
    bs = np.array([-3.5, 1.]); ll = lambda c: np.sum(ys * (c * (Xsep @ bs)) - np.log1p(np.exp(c * (Xsep @ bs))))
    assert ll(2 * c) > ll(c) and ll(2 * c) < 0
frid = glm_fit(Xsep, ys, "bernoulli", ridge=1.0, with_null=False)
b_step = np.linalg.solve(Xsep.T @ (0.25 * np.eye(6)) @ Xsep + np.diag([0, 1.]), Xsep.T @ (0.25 * (4 * ys - 2)))
assert np.allclose(glm_fit(Xsep, ys, "bernoulli", ridge=1.0, beta0=np.zeros(2), maxit=1)["beta"], b_step)
assert abs(frid["beta"][0] / frid["beta"][1] + 3.5) < 1e-6      # boundary midway between x = 3 and x = 4
print(f"Ridge (lambda = 1 on slope): beta = {np.round(frid['beta'], 4)}, {frid['iters']} iterations; first penalized step from 0: {np.round(b_step, 4)}")
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].plot(range(1, 26), sl, "o-", ms=3); ax[0].axhline(frid["beta"][1], color="r", ls="--", label="ridge, lambda = 1")
ax[0].set(title="Slope per IRLS iteration: no limit", xlabel="iteration"); ax[0].legend()
ax[1].semilogy(range(1, 26), dv); ax[1].set(title="Deviance -> 0 like e^{-t}", xlabel="iteration")
fig.tight_layout(); fig.savefig(FIG / "separation.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 8: Multinomial logistic regression by gradient descent (notes 3.7; LA-14 (3.21))

# %%
def softmax(Z_):
    E = np.exp(Z_ - Z_.max(1, keepdims=True)); return E / E.sum(1, keepdims=True)

def softmax_ce(Wf, Xm, Y):
    """Mean cross-entropy -1/n sum_i sum_k Y_ik log s_ik and its gradient X'(S - Y)/n (LA-14 (3.21))."""
    S = softmax(Xm @ Wf.reshape(Xm.shape[1], -1)); return -np.mean(np.sum(Y * np.log(S), 1)), (Xm.T @ (S - Y) / len(Y)).ravel()

rng = np.random.default_rng(8)
cent = np.array([[-2, 0.], [2, 0.], [0, 2.8]]); lab = rng.integers(0, 3, 600)
P = cent[lab] + 1.3 * rng.normal(size=(600, 2)); Xm = np.column_stack([np.ones(600), P]); Y = np.eye(3)[lab]
W0f = rng.normal(size=9) * .1
g_an = softmax_ce(W0f, Xm, Y)[1]; e9 = np.eye(9) * 1e-6
g_fd = np.array([(softmax_ce(W0f + e, Xm, Y)[0] - softmax_ce(W0f - e, Xm, Y)[0]) / 2e-6 for e in e9])
assert np.allclose(g_an, g_fd, atol=1e-7)                                        # gradient (S - Y) checked by finite differences
assert np.allclose(softmax(np.c_[np.zeros(5), np.linspace(-3, 3, 5)])[:, 1], expit(np.linspace(-3, 3, 5)))   # K=2 softmax is the logistic
assert np.allclose(softmax(Xm[:5] @ W0f.reshape(3, 3)), softmax(Xm[:5] @ W0f.reshape(3, 3) + 7.0))             # shift invariance: K-1 free logits
Wg = np.zeros(9)
for _ in range(6000):
    Wg -= 1.0 * softmax_ce(Wg, Xm, Y)[1]
ref = optimize.minimize(softmax_ce, np.zeros(9), args=(Xm, Y), jac=True, method="L-BFGS-B", options={"maxiter": 2000, "gtol": 1e-10})
acc = np.mean(np.argmax(Xm @ Wg.reshape(3, 3), 1) == lab)
print(f"Softmax: GD loss {softmax_ce(Wg, Xm, Y)[0]:.5f}, L-BFGS loss {ref.fun:.5f}, |grad| {np.linalg.norm(softmax_ce(Wg, Xm, Y)[1]):.1e}, accuracy {acc:.3f}")
assert abs(softmax_ce(Wg, Xm, Y)[0] - ref.fun) < 2e-4 and acc > 0.85
gg1, gg2 = np.meshgrid(np.linspace(-6, 6, 200), np.linspace(-5, 7, 200))
pred = np.argmax(np.c_[np.ones(gg1.size), gg1.ravel(), gg2.ravel()] @ Wg.reshape(3, 3), 1).reshape(gg1.shape)
fig, ax = plt.subplots(figsize=(5, 4.2)); ax.contourf(gg1, gg2, pred, alpha=.25, levels=[-.5, .5, 1.5, 2.5]); ax.scatter(P[:, 0], P[:, 1], c=lab, s=6, cmap="viridis")
ax.set_title("Softmax decision regions: three convex cells"); fig.tight_layout(); fig.savefig(FIG / "softmax_regions.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 9: Poisson regression, offsets, deviance and residuals (notes 3.6, Examples 4.5 and 4.7)

# %%
rng = np.random.default_rng(9)
nn_ = 400; tt = rng.uniform(1, 10, nn_); xx = rng.normal(size=nn_)
Xp = np.column_stack([np.ones(nn_), xx]); yp = rng.poisson(tt * np.exp(-1 + 0.5 * xx)).astype(float)
f_off = glm_fit(Xp, yp, "poisson", offset=np.log(tt))
f_no = glm_fit(Xp, yp, "poisson")
f_free = glm_fit(np.column_stack([Xp, np.log(tt)]), yp, "poisson")
print(f"Offset fit: beta = {np.round(f_off['beta'], 3)} (truth -1, 0.5), deviance {f_off['dev']:.1f}; no offset: deviance {f_no['dev']:.1f}; "
      f"free coefficient on log t: {f_free['beta'][2]:.3f} +/- {f_free['se'][2]:.3f}")
assert np.all(np.abs(f_off["beta"] - [-1, .5]) < 3 * f_off["se"]) and f_free["ci"][2, 0] < 1 < f_free["ci"][2, 1] and f_no["dev"] > f_off["dev"] + 100
# Example 4.5: Geiger-counter counts at dwell times t = 1, 2, 4, 8 s
tg = np.array([1, 2, 4, 8.]); yg = np.array([3, 5, 9, 21.]); Xg = np.ones((4, 1))
fo = glm_fit(Xg, yg, "poisson", offset=np.log(tg)); fn0 = glm_fit(Xg, yg, "poisson")
assert np.isclose(np.exp(fo["beta"][0]), 38 / 15) and np.isclose(fn0["mu"][0], 9.5)
D_off = 2 * np.sum(yg * np.log(yg / fo["mu"])); D_no = 2 * np.sum(yg * np.log(yg / fn0["mu"]))
assert np.isclose(D_off, fo["dev"]) and np.isclose(D_no, fn0["dev"]) and np.isclose(np.sum(yg - fo["mu"]), 0)
dres = np.sign(yg - fn0["mu"]) * np.sqrt(FAMILIES["poisson"]["dev"](yg, fn0["mu"])); pres = (yg - fn0["mu"]) / np.sqrt(fn0["mu"])
terms = yg * np.log(yg / fo["mu"])
print(f"   Ex 4.5 terms y log(y/mu) = {np.round(terms, 4)}, sum {terms.sum():.4f}; p-values: with offset {stats.chi2.sf(D_off, 3):.4f}, "
      f"no offset D {stats.chi2.sf(D_no, 3):.2e}, X2 {stats.chi2.sf(fn0['pearson'], 3):.2e}")
print(f"Ex 4.5: rate {38/15:.4f}; mu with offset {np.round(fo['mu'], 3)}, D = {D_off:.4f}; no offset: mu = 9.5, D = {D_no:.4f}, X^2 = {fn0['pearson']:.3f}")
print(f"Ex 4.7: y = 21, mu = 9.5: unit deviance {FAMILIES['poisson']['dev'](21., 9.5):.4f}, r_D = {dres[3]:.4f}, r_P = {pres[3]:.4f}; all r_D = {np.round(dres, 3)}")
assert np.isclose(dres[3] ** 2, 2 * (21 * np.log(21 / 9.5) - 11.5)) and np.isclose(np.sum(dres ** 2), D_no)

# %% [markdown]
# ## Part 10: Overdispersion: Gamma-mixed Poisson, quasi-Poisson, Negative Binomial (notes 3.6, Example 4.6)

# %%
rng = np.random.default_rng(10)
r0, m0 = 2.0, 4.0                                      # Example 4.6: NB(mean 4, r = 2): Var = mu + mu^2/r = 12
lam_s = rng.gamma(r0, m0 / r0, 1_000_000); ymix = rng.poisson(lam_s)
yv_ = np.arange(0, 6); pmf_mix = np.array([np.mean(ymix == k) for k in yv_])
pmf_nb = stats.nbinom.pmf(yv_, r0, r0 / (r0 + m0))
pmf_int = np.array([integrate.quad(lambda l: stats.poisson.pmf(k, l) * stats.gamma.pdf(l, r0, scale=m0 / r0), 0, np.inf)[0] for k in yv_])
assert np.allclose(pmf_int, pmf_nb) and np.max(np.abs(pmf_mix - pmf_nb)) < 2e-3
print(f"Ex 4.6: mixture mean {ymix.mean():.3f}, variance {ymix.var():.3f} (theory {m0 + m0**2/r0}); pmf(3): quad {pmf_int[3]:.5f} = nbinom {pmf_nb[3]:.5f}; dispersion index {ymix.var()/ymix.mean():.3f}")
assert abs(ymix.var() - 12) < 0.15

def sim_over(n, rng, b=(1.0, 0.5), r=2.0):
    x_ = rng.normal(size=n); Xo = np.column_stack([np.ones(n), x_]); m_ = np.exp(Xo @ b)
    return Xo, rng.poisson(m_ * rng.gamma(r, 1 / r, n)).astype(float)

def fit_nb(Xo, yo):
    """NB GLM with log link; r maximizes the profile log-likelihood over log r (each r is one IRLS fit)."""
    prof = lambda lr: -glm_fit(Xo, yo, negbin(np.exp(lr)), with_null=False, tol=1e-9)["loglik"]
    lr = optimize.minimize_scalar(prof, bounds=(-3, 6), method="bounded", options={"xatol": 1e-2}).x
    return glm_fit(Xo, yo, negbin(np.exp(lr)), with_null=False), np.exp(lr)

b_o = np.array([1.0, 0.5]); R = 300; hit = np.zeros(3); r_hats = []; phis = []; bh = []
for _ in range(R):
    Xo, yo = sim_over(200, rng)
    fp_ = glm_fit(Xo, yo, "poisson", with_null=False); fq = glm_fit(Xo, yo, "poisson", est_phi=True, with_null=False); fnb, rh = fit_nb(Xo, yo)
    bh.append(fp_["beta"][1]); phis.append(fq["phi"]); r_hats.append(rh)
    for k, ff in enumerate((fp_, fq, fnb)):
        hit[k] += ff["ci"][1, 0] <= b_o[1] <= ff["ci"][1, 1]
hit /= R
print(f"Slope 95% coverage ({R} reps, n=200, true r=2): Poisson {hit[0]:.3f}, quasi-Poisson {hit[1]:.3f}, NB {hit[2]:.3f}; "
      f"mean phi-hat {np.mean(phis):.2f}, median r-hat {np.median(r_hats):.2f}; empirical sd of slope {np.std(bh):.4f}")
assert hit[0] < 0.80 and hit[1] > 0.90 and hit[2] > 0.90 and 1.2 < np.median(r_hats) < 3.5
Xo, yo = sim_over(200, np.random.default_rng(11))
fp_ = glm_fit(Xo, yo, "poisson"); fq = glm_fit(Xo, yo, "poisson", est_phi=True); fnb, rh = fit_nb(Xo, yo)
assert np.allclose(fp_["beta"], fq["beta"]) and np.allclose(fq["se"], fp_["se"] * np.sqrt(fq["phi"]))
print(f"One dataset: slope se Poisson {fp_['se'][1]:.4f}, quasi {fq['se'][1]:.4f} (x sqrt(phi) = x{np.sqrt(fq['phi']):.3f}), NB {fnb['se'][1]:.4f} (r-hat {rh:.2f}); "
      f"AIC Poisson {fp_['aic']:.1f}, NB {fnb['aic'] + 2:.1f} (r counted as a parameter)")
assert fnb["aic"] + 2 < fp_["aic"]
fig, ax = plt.subplots(figsize=(5.4, 3.6)); ax.hist(bh, 25, density=True, alpha=.5, label="slope-hat over replications")
gs_ = np.linspace(min(bh), max(bh), 200)
ax.plot(gs_, stats.norm.pdf(gs_, .5, fp_["se"][1]), "r", label="Poisson claim"); ax.plot(gs_, stats.norm.pdf(gs_, .5, fq["se"][1]), "g", label="quasi-Poisson claim")
ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(FIG / "overdispersion.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 11: Analysis of deviance, chi-squared calibration, and where the deviance is not chi-squared (notes 3.8)

# %%
rng = np.random.default_rng(12)
n = 150; Xa = np.column_stack([np.ones(n), rng.normal(size=(n, 3))]); ya = rng.poisson(np.exp(0.8 + 0.4 * Xa[:, 1] + 0.2 * Xa[:, 2])).astype(float)
fits = [glm_fit(Xa[:, :k], ya, "poisson", with_null=False) for k in (1, 2, 3, 4)]
print("Analysis of deviance (Poisson, n=150): model, resid. dev, df, change, p")
for k in range(4):
    dD = fits[k - 1]["dev"] - fits[k]["dev"] if k else np.nan
    print(f"  +x{k}: {fits[k]['dev']:8.3f} {fits[k]['df_res']:4d}   {dD:8.3f}   {stats.chi2.sf(dD, 1) if k else np.nan:.4f}")
print("   p-values:", [f"{stats.chi2.sf(fits[k-1]['dev'] - fits[k]['dev'], 1):.2e}" for k in (1, 2, 3)])
assert stats.chi2.sf(fits[0]["dev"] - fits[1]["dev"], 1) < 1e-3 and stats.chi2.sf(fits[2]["dev"] - fits[3]["dev"], 1) > .001
R = 2000; lr_stat = np.empty(R)
for r_ in range(R):
    Xr = np.column_stack([np.ones(100), rng.normal(size=100)]); yr = rng.poisson(5.0, 100).astype(float)
    lr_stat[r_] = glm_fit(Xr[:, :1], yr, "poisson", with_null=False)["dev"] - glm_fit(Xr, yr, "poisson", with_null=False)["dev"]
rej = np.mean(lr_stat > stats.chi2.ppf(.95, 1)); ks = stats.kstest(lr_stat, stats.chi2(1).cdf).pvalue
print(f"LRT under H0 ({R} reps): rejection rate at 5% = {rej:.4f}, KS p = {ks:.3f}")
assert abs(rej - .05) < .015 and ks > .01
# Bernoulli: D depends on the data only through the fitted values (3.24), so it cannot test fit
Dber, Dgrp, ident = [], [], []
xg = np.linspace(-1.5, 1.5, 8); Xg8 = np.column_stack([np.ones(8), xg])
for _ in range(1000):
    xb_ = rng.normal(size=100); Xb_ = np.column_stack([np.ones(100), xb_]); yb_ = rng.binomial(1, expit(3 * xb_)).astype(float)
    fb = glm_fit(Xb_, yb_, "bernoulli", with_null=False, tol=1e-12); Dber.append(fb["dev"])
    th_ = Xb_ @ fb["beta"]; ident.append(abs(fb["dev"] - (-2) * (fb["beta"] @ Xb_.T @ fb["mu"] - np.sum(np.log1p(np.exp(th_))))))
    mg = 100; ybin = rng.binomial(mg, expit(1.0 * xg)) / mg
    Dgrp.append(glm_fit(Xg8, ybin, "bernoulli", weights=np.full(8, float(mg)), with_null=False)["dev"])
print(f"Bernoulli n=100, beta1=3: mean D = {np.mean(Dber):.1f} vs n-p = 98 (sd {np.std(Dber):.1f}); chi2_98 would have sd 14.0. "
      f"Grouped binomial (8 x 100): mean D = {np.mean(Dgrp):.2f} vs n-p = 6, "
      f"P(D > chi2_6 95% point) = {np.mean(np.array(Dgrp) > stats.chi2.ppf(.95, 6)):.3f}")
assert max(ident) < 1e-6 and np.mean(Dber) < 80 and abs(np.mean(Dgrp) - 6) < 0.3 and abs(np.mean(np.array(Dgrp) > stats.chi2.ppf(.95, 6)) - .05) < .02
# Exercise B6: intercept-only Bernoulli, n = 50, ybar = 0.3: D = 61.09 whatever the order of the data
D6 = -2 * 50 * (.3 * np.log(.3) + .7 * np.log(.7)); assert np.isclose(D6, glm_fit(np.ones((50, 1)), np.r_[np.ones(15), np.zeros(35)], "bernoulli", with_null=False)["dev"])
print(f"Ex B6: D = {D6:.3f}, P(chi2_49 > D) = {stats.chi2.sf(D6, 49):.4f}")

# %% [markdown]
# ## Part 12: Residuals, leverage, Cook's distance (notes 3.8)

# %%
fq = glm_fit(Xa[:, :3], ya, "poisson"); mu_a = fq["mu"]
r_d = np.sign(ya - mu_a) * np.sqrt(FAMILIES["poisson"]["dev"](ya, mu_a)); r_p = (ya - mu_a) / np.sqrt(mu_a)
Wh = np.sqrt(fq["w"])[:, None] * Xa[:, :3]; Hm = Wh @ np.linalg.inv(Wh.T @ Wh) @ Wh.T
assert np.isclose(np.trace(Hm), 3) and np.allclose(Hm, Hm.T) and np.allclose(Hm @ Hm, Hm) and np.allclose(np.diag(Hm), fq["h"])
cook = r_p ** 2 * fq["h"] / (3 * fq["phi"] * (1 - fq["h"]) ** 2); r_std = r_d / np.sqrt(1 - fq["h"])
print(f"Residuals: sum r_D^2 = {np.sum(r_d**2):.3f} = D; sum r_P^2 = {np.sum(r_p**2):.3f} = X^2; sum h = {fq['h'].sum():.3f}; max Cook {cook.max():.3f}; max |std r_D| {np.abs(r_std).max():.2f}")
assert np.isclose(np.sum(r_d ** 2), fq["dev"]) and np.isclose(np.sum(r_p ** 2), fq["pearson"])
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
ax[0].scatter(np.log(mu_a), r_d, s=8); ax[0].scatter(np.log(mu_a), r_p, s=8, alpha=.4); ax[0].axhline(0, color="k", lw=.5); ax[0].set_title("deviance (blue), Pearson (orange)")
stats.probplot(r_d, plot=ax[1]); ax[1].set_title("QQ of deviance residuals"); ax[2].stem(cook, markerfmt=".", basefmt=" "); ax[2].set_title("Cook's distance")
fig.tight_layout(); fig.savefig(FIG / "residuals.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Exercise checks (small numbers quoted in exercises.md and solutions.md)

# %%
assert np.isclose(stats.binom(10, .3).var() / 100, .3 * .7 / 10)                                   # A1: a(phi) = 1/m
p2_, p3_ = expit(-1.4), expit(-0.6)                                                                  # A2
print(f"A2: p(2) = {p2_:.4f}, p(3) = {p3_:.4f}, OR = {np.exp(.8):.4f}, RR = {p3_/p2_:.4f}, ME(3) = {.8*p3_*(1-p3_):.4f}, p=.5 at x = {3/.8}")
assert np.isclose((p3_ / (1 - p3_)) / (p2_ / (1 - p2_)), np.exp(.8))
mu_dry, mu_wet = 3.2 * np.exp(-2.0), 3.2 * np.exp(-1.5)                                              # A3
print(f"A3: dry {mu_dry:.4f}, wet {mu_wet:.4f}, RR = {np.exp(.5):.4f}, P(0 | dry) = {np.exp(-mu_dry):.4f}, P(0 | wet) = {np.exp(-mu_wet):.4f}")
ya4, mu4 = np.array([0, 2, 7.]), np.array([1.2, 3.1, 4.4])                                           # A4
d4 = FAMILIES["poisson"]["dev"](ya4, mu4); rd4 = np.sign(ya4 - mu4) * np.sqrt(d4); rp4 = (ya4 - mu4) / np.sqrt(mu4)
print(f"A4: unit deviances {np.round(d4, 4)}, D = {d4.sum():.4f}; r_D = {np.round(rd4, 4)}; r_P = {np.round(rp4, 4)}; X^2 = {np.sum(rp4**2):.4f}")
yy5, b_ = np.array([2., 4., 6.]), 0.0                                                                # A5
seq = [b_]
for _ in range(7):
    m_ = np.exp(b_); b_ = np.mean(b_ + (yy5 - m_) / m_); seq.append(b_)
print(f"A5: intercept-only Poisson IRLS from 0 on y = (2,4,6): {np.round(seq, 4)}; log 4 = {np.log(4):.4f}")
print(f"   errors |beta_t - log 4| = {np.array2string(np.abs(np.array(seq) - np.log(4)), precision=2)}")
assert np.isclose(seq[1], 3) and abs(seq[-1] - np.log(4)) < 1e-6 and np.isclose(glm_fit(np.ones((3, 1)), yy5, "poisson", beta0=np.array([np.log(4.)]), maxit=1)["beta"][0], np.log(4))
r_ = 10 / 3; assert np.isclose(1 + 5 / r_, 2.5)                                                       # B5: V/mu = 1 + mu/r = 2.5 at mu = 5
print(f"B5: weights at mu = 5, r = {r_:.3f}: quasi {5:.1f}, NB {5/(1+5/r_):.3f}; mu = 50: quasi 50, NB {50/(1+50/r_):.3f}")
print(f"B7: slope match {4*stats.norm.pdf(0):.4f}; Phi(1.5) = {stats.norm.cdf(1.5):.4f}, sigmoid(1.7*1.5) = {expit(1.7*1.5):.4f}")
f_b1 = glm_fit(np.ones((3, 1)), np.array([1., 2., 6.]), "gamma", "log")                          # B1
assert np.isclose(f_b1["beta"][0], np.log(3)) and np.allclose(f_b1["w"], 1)
fpr_b, tpr_b, _ = roc_curve(yy, ss)                                                                  # B2
assert np.allclose(fpr_b, [0, 0, .25, .5, .5, .75, 1]) and np.allclose(tpr_b, [0, .25, .25, .75, 1, 1, 1])
ll1 = lambda b: 3 * b - np.exp(b) - gammaln(4); h_ = 1e-6                                           # B1(c): one observation y=3, mu = e^beta, at beta = 1
true_u = (ll1(1 + h_) - ll1(1 - h_)) / (2 * h_); wrong_u = (3 - np.e) / np.e
assert np.isclose(true_u, 3 - np.e, atol=1e-8) and np.isclose(true_u, 0.2817, atol=5e-5) and np.isclose(wrong_u, 0.1036, atol=5e-5)
assert np.allclose(fsep["path"][0], [-3.6, 108 / 105]) and np.isclose(np.linalg.det(Xsep.T @ Xsep), 105) and np.allclose(Xsep.T @ Xsep, [[6, 21], [21, 91]])   # B3
print(f"B3: pass 1 from 0 = {np.round(fsep['path'][0], 4)}; ridge pass 1 = {np.round(b_step, 4)}")
print(f"B4 (z for offset=1): z = {(0.93-1)/0.06:.3f}, p = {2*stats.norm.sf(abs((0.93-1)/0.06)):.3f}")

# %% [markdown]
# ## Your turn
# (a) Gamma GLM versus OLS on log y; (b) a zero-inflated Poisson; (c) identity-link Poisson, Newton versus Fisher scoring;
# (d) Exercise B8: add the Inverse Gaussian family to `FAMILIES` (reference solution at the bottom). Reference code follows.

# %%
rng = np.random.default_rng(20)                                                      # (a) / Exercise C2
n = 500; xa_ = rng.normal(size=n); Xa_ = np.column_stack([np.ones(n), xa_]); nu_ = 3.0
ya_ = np.exp(1 + 0.5 * xa_) * rng.gamma(nu_, 1 / nu_, n)
fgm = glm_fit(Xa_, ya_, "gamma", "log"); b_lo = np.linalg.lstsq(Xa_, np.log(ya_), rcond=None)[0]
shift = digamma(nu_) - np.log(nu_)
print(f"C2: Gamma GLM beta = {np.round(fgm['beta'], 3)}, phi-hat = {fgm['phi']:.3f} (1/nu = {1/nu_:.3f}); log-OLS beta = {np.round(b_lo, 3)}; "
      f"intercept shift psi(nu) - log(nu) = {shift:.4f}, e^shift = {np.exp(shift):.4f}")
nu_ml = 1 / _gamma_phi_ml(ya_, fgm["mu"], np.ones(n), fgm["dev"])                    # C2(d): log nu - psi(nu) = D / (2n)
assert np.isclose(np.log(nu_ml) - digamma(nu_ml), fgm["dev"] / (2 * n)) and abs(nu_ml - nu_) < 0.5
print(f"   C2(d): Gamma ML shape nu-hat = {nu_ml:.3f} (Pearson 1/phi-hat = {1/fgm['phi']:.3f}), deviance {fgm['dev']:.3f}")
assert abs(fgm["beta"][0] - 1) < 3 * fgm["se"][0] and abs(b_lo[0] - (1 + shift)) < 0.1 and abs(b_lo[1] - .5) < .1
print(f"   C2(c): mean(exp(OLS fit))/mean(y) = {np.exp(Xa_ @ b_lo).mean() / ya_.mean():.4f} vs e^shift = {np.exp(shift):.4f}")
assert np.isclose(fgm["phi"], 1 / nu_, atol=0.05) and abs(np.exp(Xa_ @ b_lo).mean() / ya_.mean() - np.exp(shift)) < 0.08

def zip_nll(par, Xz, yz):
    """Zero-inflated Poisson: P(0) = pi + (1-pi) e^{-mu}, P(k) = (1-pi) Pois(k; mu)  (Exercise C3)."""
    b, pi = par[:-1], expit(par[-1]); m = np.exp(Xz @ b)
    return -np.sum(np.where(yz == 0, np.log(pi + (1 - pi) * np.exp(-m)), np.log(1 - pi) + stats.poisson.logpmf(yz, m)))

k_ = np.arange(0, 400); pk_ = 0.7 * stats.poisson.pmf(k_, 2.0); pk_[0] += 0.3                      # C3(b): ZIP mean and variance
m1_ = np.sum(k_ * pk_); v1_ = np.sum(k_ ** 2 * pk_) - m1_ ** 2
assert np.isclose(m1_, 0.7 * 2) and np.isclose(v1_, 0.7 * 2 * (1 + 0.3 * 2)) and np.isclose(v1_, 2.24)
rng = np.random.default_rng(21)                                                      # (b)
n = 1500; xz = rng.normal(size=n); Xz = np.column_stack([np.ones(n), xz]); mz = np.exp(0.8 + 0.5 * xz)
yz = np.where(rng.random(n) < 0.3, 0, rng.poisson(mz)).astype(float)
rz = optimize.minimize(zip_nll, np.zeros(3), args=(Xz, yz), method="BFGS"); fpz = glm_fit(Xz, yz, "poisson")
print(f"C3: ZIP beta = {np.round(rz.x[:2], 3)}, pi = {expit(rz.x[2]):.3f} (truth 0.8, 0.5, 0.3); zeros observed {np.mean(yz==0):.3f}, "
      f"Poisson-GLM expects {np.mean(np.exp(-fpz['mu'])):.3f}; ZIP logL {-rz.fun:.1f} vs Poisson {fpz['loglik']:.1f}")
assert abs(expit(rz.x[2]) - .3) < .06 and np.all(np.abs(rz.x[:2] - [.8, .5]) < .08) and -rz.fun > fpz["loglik"] + 50

rng = np.random.default_rng(22)                                                      # (c) / Exercise C1
n = 100; xi = rng.uniform(0, 2, n); Xi = np.column_stack([np.ones(n), xi]); yi = rng.poisson(2 + 3 * xi).astype(float)
fF = glm_fit(Xi, yi, "poisson", "identity", tol=0, maxit=40); fN = glm_fit(Xi, yi, "poisson", "identity", observed=True, tol=0, maxit=40)
fF0 = glm_fit(Xi, yi, "poisson", "identity"); fN0 = glm_fit(Xi, yi, "poisson", "identity", observed=True)
tF = glm_fit(Xi, yi, "poisson", "identity", beta0=np.array([1., 1.]), tol=0, maxit=8)["path"]
tN = glm_fit(Xi, yi, "poisson", "identity", beta0=np.array([1., 1.]), observed=True, tol=0, maxit=8)["path"]
eF, eN = np.linalg.norm(tF - fF["beta"], axis=1), np.linalg.norm(tN - fF["beta"], axis=1)
print(f"C1: identity-link Poisson, beta-hat = {np.round(fF['beta'], 4)}; error per iteration from beta=(1,1):\n  Fisher {np.array2string(eF, precision=1)}\n  Newton {np.array2string(eN, precision=1)}")
print(f"   iterations at the default tolerance: Fisher {fF0['iters']}, Newton {fN0['iters']}; se from expected info {np.round(fF['se'], 4)}, from observed info {np.round(np.sqrt(np.diag(np.linalg.inv(Xi.T @ ((yi / fF['mu'] ** 2)[:, None] * Xi)))), 4)}")
gF, vF, dVF, gppF = LINKS["identity"][2](fF["mu"]), FAMILIES["poisson"]["V"](fF["mu"]), FAMILIES["poisson"]["dV"](fF["mu"]), LINKS["identity"][3](fF["mu"])
w_obs_ = 1 / (vF * gF ** 2) + (yi - fF["mu"]) * (dVF * gF + vF * gppF) / (vF ** 2 * gF ** 3)
assert np.allclose(w_obs_, yi / fF["mu"] ** 2)                                      # C1(a): w_obs = y / mu^2 for the identity link
assert np.allclose(fF["beta"], fN["beta"], atol=1e-10) and eN[4] < 10 * eN[3] ** 2 and np.all((eF[1:4] / eF[:3] > .01) & (eF[1:4] / eF[:3] < .1))   # Newton quadratic, Fisher linear (rate ~0.04)

# %% [markdown]
# ## Reference solution to Exercise B8 (write your own first): the Inverse Gaussian family

# %%
FAMILIES["invgauss"] = dict(V=lambda m: m ** 3, dV=lambda m: 3 * m ** 2, dev=lambda y, m: (y - m) ** 2 / (y * m ** 2), default="log",
                            phi_known=False, init=lambda y, om: y.astype(float), ll=lambda y, m, om, ph: stats.invgauss.logpdf(y, m / (1 / ph), scale=1 / ph).sum())
rng = np.random.default_rng(23)
n = 400; xg_ = rng.normal(size=n); Xg_ = np.column_stack([np.ones(n), xg_]); mg_ = np.exp(0.5 + 0.4 * xg_); lam_ = 5.0
yg_ = stats.invgauss.rvs(mg_ / lam_, scale=lam_, random_state=rng)
fig_ = glm_fit(Xg_, yg_, "invgauss", "log")
b_dir = optimize.minimize(lambda b: -stats.invgauss.logpdf(yg_, np.exp(Xg_ @ b) / lam_, scale=lam_).sum(), np.zeros(2), method="BFGS", options={"gtol": 1e-9}).x
print(f"B8: Inverse Gaussian/log IRLS beta = {np.round(fig_['beta'], 4)}, direct ML {np.round(b_dir, 4)}, phi-hat = {fig_['phi']:.3f} (1/lambda = {1/lam_:.3f})")
assert np.allclose(fig_["beta"], b_dir, atol=1e-4) and abs(fig_["phi"] - 1 / lam_) < 0.05

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
