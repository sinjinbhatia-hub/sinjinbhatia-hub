# %% [markdown]
# # ST-3 Project: Least Squares, Inference, Diagnostics, Ridge and Lasso
#
# Everything in the notes is implemented from scratch (no statsmodels, no sklearn) and checked against
# closed forms, brute force or Monte Carlo. Parts: (1) OLS three ways and the normal-equations failure,
# (2) the inference table, Examples 4.1-4.2 (FWL), (3) Gauss-Markov and coverage over 10^4 replications,
# (4) leverage, Cook's distance, the LOO identity, (5) heteroscedasticity: sandwich and WLS,
# (6) multicollinearity: VIF and the SVD, (7) model selection and post-selection, (8) ridge, (9) lasso by
# coordinate descent, (10) numbers used in exercises and solutions, then "Your turn" with reference solutions.
# Every number quoted in notes.md, exercises.md or solutions.md is asserted or printed here.
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg, optimize, stats

T0 = time.time()
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
np.set_printoptions(precision=4, suppress=True)
close = lambda a, b, tol=5e-4: bool(np.all(np.abs(np.asarray(a) - np.asarray(b)) < tol))


def hat_diag(X, lam=0.0):
    """h_ii of H = X (X'X + lam I)^{-1} X' (notes (3.5); lam > 0 is ridge, (3.39))."""
    return np.einsum("ij,ji->i", X, np.linalg.solve(X.T @ X + lam * np.eye(X.shape[1]), X.T))


# %% [markdown]
# ## Part 1: OLS three ways (notes 3.2; LA-5 (3.6), (3.16); LA-9 (3.11); LA-13 Theorem 3.12)

# %%
def ols_normal(X, y):
    """Normal equations X'X b = X'y (3.4). Squares the condition number."""
    return np.linalg.solve(X.T @ X, X.T @ y)


def ols_qr(X, y):
    """X = QR, then R b = Q'y by back substitution (LA-5 (3.16))."""
    Q, R = np.linalg.qr(X)
    return linalg.solve_triangular(R, Q.T @ y)


def ols_svd(X, y):
    """b = X^+ y = V diag(1/sigma) U'y (LA-9 (3.11))."""
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    return Vt.T @ ((U.T @ y) / s)


def make_design(n, s, rng):
    """X = U diag(s) V' with random orthonormal U (n x p), V (p x p): singular values exactly s."""
    U = np.linalg.qr(rng.standard_normal((n, len(s))))[0]
    V = np.linalg.qr(rng.standard_normal((len(s), len(s))))[0]
    return U @ np.diag(s) @ V.T


rng = np.random.default_rng(1)
Xw = rng.standard_normal((60, 5)); yw = Xw @ np.arange(1.0, 6) + 0.1 * rng.standard_normal(60)
for f in (ols_normal, ols_qr, ols_svd):
    assert np.allclose(f(Xw, yw), np.linalg.lstsq(Xw, yw, rcond=None)[0], atol=1e-10)
Xi = make_design(60, np.logspace(0, -8, 5), rng)  # kappa(X) = 1e8, so kappa(X'X) = 1e16
yi = Xi @ np.ones(5)                               # zero residual: the error below is pure conditioning
err = {f.__name__: np.linalg.norm(f(Xi, yi) - np.ones(5)) for f in (ols_normal, ols_qr, ols_svd)}
print("Part 1  kappa(X) = %.1e;  error of beta-hat:" % np.linalg.cond(Xi), {k: "%.1e" % v for k, v in err.items()})
assert err["ols_qr"] < 1e-8 and err["ols_svd"] < 1e-8 and err["ols_normal"] > 0.1  # normal equations: no correct digits

# %% [markdown]
# ## Part 2: the inference table (notes 3.5), Example 4.1, Example 4.2 (FWL)

# %%
def ols_fit(X, y, alpha=0.05):
    """OLS with intercept in column 0: b (3.5), se_j (3.14), t, p, CI (3.15), RSS, TSS, R2, adjR2 (3.8), s^2, overall F (3.18),
    leverages h_ii (3.21) and residuals e."""
    n, p = X.shape
    A = np.linalg.inv(X.T @ X); b = A @ X.T @ y; e = y - X @ b
    rss, tss, df = e @ e, ((y - y.mean()) ** 2).sum(), n - p
    s2 = rss / df; se = np.sqrt(s2 * np.diag(A)); tc = stats.t.ppf(1 - alpha / 2, df); F = ((tss - rss) / (p - 1)) / s2
    return dict(b=b, se=se, t=b / se, p=2 * stats.t.sf(abs(b / se), df), lo=b - tc * se, hi=b + tc * se, rss=rss, tss=tss,
                R2=1 - rss / tss, adjR2=1 - (rss / df) / (tss / (n - 1)), s2=s2, F=F, pF=stats.f.sf(F, p - 1, df),
                h=np.einsum("ij,jk,ik->i", X, A, X), e=e, A=A, df=df, tc=tc)


x41, y41 = np.arange(1.0, 6), np.array([2.0, 4, 5, 4, 5]); X41 = np.c_[np.ones(5), x41]; f41 = ols_fit(X41, y41)
assert (((x41 - 3) ** 2).sum(), ((x41 - 3) * (y41 - 4)).sum()) == (10.0, 6.0) and close(f41["b"], [2.2, 0.6])
assert close(f41["e"], [-0.8, 0.6, 1.0, -0.6, -0.2]) and close([f41["rss"], f41["tss"]], [2.4, 6.0])
assert close([f41["R2"], f41["adjR2"], f41["s2"]], [0.6, 0.46667, 0.8]) and close(f41["se"], [0.93808, 0.28284], 5e-5)
assert close([f41["t"][1], f41["tc"], f41["p"][1]], [2.1213, 3.1824, 0.1240]) and close(f41["F"], f41["t"][1] ** 2, 1e-9) and close(f41["F"], 4.5)
assert close([f41["lo"][1], f41["hi"][1]], [-0.3001, 1.5001]) and close(f41["h"], [0.6, 0.3, 0.2, 0.3, 0.6]) and abs(f41["h"].sum() - 2) < 1e-12
for x0, ch, ph in [(3.0, 1.2730, 3.1181), (5.0, 2.2049, 3.6005)]:  # CI for the mean vs PI for a new y, (3.16)-(3.17)
    v = np.array([1.0, x0]); h0 = v @ f41["A"] @ v; w = f41["tc"] * np.sqrt(f41["s2"])
    assert close([w * np.sqrt(h0), w * np.sqrt(1 + h0)], [ch, ph]) and close([v @ f41["b"], h0], [4.0 if x0 == 3 else 5.2, 0.2 if x0 == 3 else 0.6])
print("Example 4.1: b =", f41["b"], " RSS =", round(f41["rss"], 6), " t = %.4f  p = %.4f  CI = (%.3f, %.3f)" % (f41["t"][1], f41["p"][1], f41["lo"][1], f41["hi"][1]))

# Example 4.2: Frisch-Waugh-Lovell with an intercept and three predictors, both routes.
x1 = np.arange(1.0, 9); x2 = np.array([2, 1, 4, 3, 6, 5, 8, 6.0]); x3 = np.array([1, 0, 2, 1, 3, 1, 2, 4.0])
y42 = np.array([3.1, 2.2, 5.9, 5.0, 8.1, 6.8, 9.9, 9.3]); X42 = np.c_[np.ones(8), x1, x2, x3]
b42 = ols_qr(X42, y42); X_o = X42[:, :3]
M1 = np.eye(8) - X_o @ np.linalg.solve(X_o.T @ X_o, X_o.T)  # annihilator of the other columns (LA-5 Theorem 3.13)
xt, yt = M1 @ x3, M1 @ y42; b_fwl = (xt @ yt) / (xt @ xt); res_full = y42 - X42 @ b42
assert close(b42, [0.9948, 0.2529, 0.7944, 0.3883], 5e-5) and close(b_fwl, b42[3], 1e-12)  # (3.9)
assert close(xt, [0.3727, -0.5206, 0.5674, -0.3258, 0.7622, -1.1311, -1.0431, 1.3184], 5e-5)
assert close([xt @ yt, xt @ xt, res_full @ res_full], [2.1449, 5.5243, 0.12579], 5e-5) and close(yt - b_fwl * xt, res_full, 1e-12)
assert close(np.linalg.inv(X42.T @ X42)[3, 3], 1 / (xt @ xt), 1e-12)  # Var(b_j) = sigma^2/||x~_j||^2 (3.31)
se3 = np.sqrt(res_full @ res_full / 4 / (xt @ xt)); naive = ols_qr(np.c_[np.ones(8), x3], y42)[1]  # df = n - p = 4, not n - 3
assert close(naive, 1.737) and close([se3, b42[3] / se3], [0.0755, 5.146])
print("Example 4.2: full b =", b42, " FWL slope = %.5f; y on x3 alone: %.3f; se = %.4f, t = %.3f" % (b_fwl, naive, se3, b42[3] / se3))

# %% [markdown]
# ## Part 3: Gauss-Markov, E and Var of b, unbiased s^2, coverage over 10^4 replications (notes 3.4-3.5)

# %%
rng = np.random.default_rng(3)
c = np.array([0.4, 0.3, 0.2, 0.1]); d = c - 0.25  # Example 4.3: other linear unbiased weights for the mean, n = 4
assert close([c @ c, 0.25, d @ d, 0.25 + d @ d], [0.30, 0.25, 0.05, 0.30], 1e-12) and abs(d.sum()) < 1e-15
n, p, sig, R = 30, 4, 2.0, 10_000
X = np.c_[np.ones(n), rng.uniform(-1, 1, (n, 3))]; beta = np.array([1.0, 2.0, -1.0, 0.5]); x0 = np.array([1.0, 0.3, -0.2, 0.5])
A = np.linalg.solve(X.T @ X, X.T); H = X @ A; V = np.linalg.inv(X.T @ X)
Y = (X @ beta)[:, None] + sig * rng.standard_normal((n, R)); B = A @ Y
RSS = ((Y - H @ Y) ** 2).sum(0); s2 = RSS / (n - p)
assert np.allclose(np.cov(B), sig**2 * V, rtol=0.1, atol=0.01) and abs(B.mean(1) - beta).max() < 0.05  # (3.11)
assert abs(s2.mean() / sig**2 - 1) < 0.02 and abs((RSS / n).mean() / (sig**2 * (n - p) / n) - 1) < 0.02  # (3.12), (3.19)
tq = stats.t.ppf(0.975, n - p); h0 = x0 @ V @ x0; ynew = x0 @ beta + sig * rng.standard_normal(R)
cov_beta = (np.abs(B - beta[:, None]) <= tq * np.sqrt(s2[None, :] * np.diag(V)[:, None])).mean(1)
cov_pred = (np.abs(ynew - x0 @ B) <= tq * np.sqrt(s2 * (1 + h0))).mean()
cov_mean = (np.abs(x0 @ beta - x0 @ B) <= tq * np.sqrt(s2 * h0)).mean()
print("Part 3  coverage of beta_j intervals:", cov_beta, " mean-response CI: %.4f  prediction interval: %.4f" % (cov_mean, cov_pred))
assert np.all(np.abs(cov_beta - 0.95) < 0.005) and abs(cov_pred - 0.95) < 0.005 and abs(cov_mean - 0.95) < 0.005
assert (np.abs(ynew - x0 @ B) <= tq * np.sqrt(s2 * h0)).mean() < 0.5  # a CI half-width used as a PI covers far too little
Y0 = sig * rng.standard_normal((n, R)); R2 = 1 - ((Y0 - H @ Y0) ** 2).sum(0) / ((Y0 - Y0.mean(0)) ** 2).sum(0)  # all slopes zero
adj = 1 - (1 - R2) * (n - 1) / (n - p)
print("Part 3  null R^2 mean %.4f (theory (p-1)/(n-1) = %.4f), null adjusted R^2 mean %.4f" % (R2.mean(), (p - 1) / (n - 1), adj.mean()))
assert abs(R2.mean() - (p - 1) / (n - 1)) < 0.003 and abs(adj.mean()) < 0.005  # (3.8)
X0 = X[:, :2]; P0 = lambda Z: X0 @ np.linalg.solve(X0.T @ X0, X0.T @ Z)
Fs = (((Y - P0(Y)) ** 2).sum(0) - RSS) / 2 / s2
assert np.allclose(n * np.log(1 + 2 * Fs / (n - p)), n * np.log(((Y - P0(Y)) ** 2).sum(0) / RSS))  # F test = LRT (3.20)
Yn = (X0 @ beta[:2])[:, None] + sig * rng.standard_normal((n, R)); r1 = ((Yn - H @ Yn) ** 2).sum(0); r0 = ((Yn - P0(Yn)) ** 2).sum(0)
rej = (((r0 - r1) / 2) / (r1 / (n - p)) > stats.f.ppf(0.95, 2, n - p)).mean(); assert abs(rej - 0.05) < 0.01
print("Part 3  nested F test (q = 2) size under H0: %.4f" % rej)

# %% [markdown]
# ## Part 4: leverage, studentized residuals, Cook's distance, the LOO identity (notes 3.6, Example 4.4)

# %%
def influence(X, y):
    """h_ii (3.21), r_i = e_i/(s sqrt(1-h_ii)) (3.23), external t_i via RSS_(i) (3.26), Cook D_i = r_i^2 h_ii/(p(1-h_ii)) (3.27)."""
    n, p = X.shape; f = ols_fit(X, y); e, h = f["e"], f["h"]
    r = e / np.sqrt(f["s2"] * (1 - h)); s2_i = (f["rss"] - e**2 / (1 - h)) / (n - p - 1)
    return h, r, e / np.sqrt(s2_i * (1 - h)), r**2 * h / (p * (1 - h)), f


def brute_loo(X, y):
    """Refit n times: Cook's D by its definition, LOO prediction errors (3.25), RSS_(i) (3.26)."""
    n, p = X.shape; b = ols_qr(X, y); s2 = ((y - X @ b) ** 2).sum() / (n - p); D, loo, rss_i = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        m = np.arange(n) != i; bi = ols_qr(X[m], y[m])
        D[i] = np.sum((X @ (b - bi)) ** 2) / (p * s2); loo[i] = y[i] - X[i] @ bi; rss_i[i] = ((y[m] - X[m] @ bi) ** 2).sum()
    return D, loo, rss_i


xl = np.array([1, 2, 3, 4, 10.0])
for y5, D5, r5, hasbig in [(9.9, 0.2163, 0.194, False), (3.0, 16.961, -1.7175, True)]:
    yl = np.array([1.1, 1.9, 3.2, 3.9, y5]); Xl = np.c_[np.ones(5), xl]; h, r, tt, D, fl = influence(Xl, yl); Db, loo, _ = brute_loo(Xl, yl)
    assert close(h, [0.38, 0.28, 0.22, 0.2, 0.92]) and close(D5, D[4], 1e-3) and close(r5, r[4], 1e-3)
    assert close(D, Db, 1e-10) and close(loo, fl["e"] / (1 - h), 1e-10) and (D[4] > 1) == hasbig and (hasbig or D[4] < D[:4].max())
    print("Example 4.4  y5=%.1f  slope=%.3f  e5=%.3f  LOO error=%.2f  r5=%.3f  D5=%.3f  max D(others)=%.3f" % (y5, fl["b"][1], fl["e"][4], loo[4], r[4], D[4], D[:4].max()))

rng = np.random.default_rng(4)
xs = np.r_[rng.uniform(0, 10, 27), 4.0, 5.5, 25.0]; ys = 1 + 0.8 * xs + 0.8 * rng.standard_normal(30)
ys[27] += 6.0; ys[28] -= 5.0                       # two planted outliers in y at ordinary x; index 29: high leverage, on trend
Xs = np.c_[np.ones(30), xs]; h, r, tt, D, fs = influence(Xs, ys); Db, loo, rss_i = brute_loo(Xs, ys)
assert np.isclose(h.sum(), 2) and np.all((h >= 1 / 30 - 1e-12) & (h <= 1)) and np.allclose(h, 1 / 30 + (xs - xs.mean()) ** 2 / ((xs - xs.mean()) ** 2).sum())
assert np.allclose(D, Db) and np.allclose(loo, fs["e"] / (1 - h)) and np.allclose(rss_i, fs["rss"] - fs["e"] ** 2 / (1 - h))
assert set(np.argsort(-np.abs(tt))[:2]) == {27, 28} and h[29] > 2 * 2 / 30 and h[29] == h.max()
print("Part 4  planted outliers |t_ext| =", np.abs(tt[[27, 28]]).round(2), " leverage point h = %.3f (2p/n = %.3f), D = %.3f" % (h[29], 4 / 30, D[29]))
fig, ax = plt.subplots(2, 2, figsize=(10, 7))
ax[0, 0].scatter(Xs @ fs["b"], fs["e"]); ax[0, 0].axhline(0, c="k", lw=.5); ax[0, 0].set(title="residuals vs fitted", xlabel="fitted", ylabel="e")
(q, o), _ = stats.probplot(tt, dist=stats.t, sparams=(27,)); ax[0, 1].plot(q, o, "o"); ax[0, 1].plot(q, q, "k-", lw=.5); ax[0, 1].set(title="QQ plot of studentized residuals", xlabel="t quantile")
ax[1, 0].scatter(h, r**2); ax[1, 0].axvline(4 / 30, ls="--", c="gray"); ax[1, 0].set(title="leverage vs squared std. residual", xlabel="h_ii", ylabel="r_i^2")
for i in (27, 28, 29): ax[1, 0].annotate(str(i), (h[i], r[i] ** 2))
ax[1, 1].stem(D); ax[1, 1].set(title="Cook's distance", xlabel="observation")
fig.tight_layout(); fig.savefig(FIG / "diagnostics.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 5: heteroscedasticity: naive vs sandwich standard errors, and WLS (notes 3.6)

# %%
rng = np.random.default_rng(5)
n, R = 60, 4000
xh = np.sort(rng.uniform(0, 10, n)); Xh = np.c_[np.ones(n), xh]; sd = 0.15 * xh**2 + 0.2  # planted pattern: sd grows like x^2
Ah = np.linalg.solve(Xh.T @ Xh, Xh.T); hh = np.einsum("ij,ji->i", Xh, Ah)
Yh = (Xh @ [1.0, 0.5])[:, None] + sd[:, None] * rng.standard_normal((n, R)); Bh = Ah @ Yh; Eh = Yh - Xh @ Bh
w1 = Ah[1]                                           # slope error = sum_i w1_i eps_i, so Var = sum_i w1_i^2 sigma_i^2 (3.28)
v_naive = (Eh**2).sum(0) / (n - 2) * (w1**2).sum()
v_hc0 = (Eh**2).T @ w1**2                            # sandwich (3.29): replace sigma_i^2 by e_i^2
v_hc3 = ((Eh / (1 - hh)[:, None]) ** 2).T @ w1**2    # HC3: LOO residuals (3.25)
tqh = stats.t.ppf(0.975, n - 2)
cov = {k: float((np.abs(Bh[1] - 0.5) <= tqh * np.sqrt(v)).mean()) for k, v in [("naive", v_naive), ("HC0", v_hc0), ("HC3", v_hc3)]}
true_var = (sd**2 * w1**2).sum(); Wt = 1 / sd**2
Bw = np.linalg.solve(Xh.T @ (Wt[:, None] * Xh), Xh.T * Wt) @ Yh; var_wls = np.linalg.inv(Xh.T @ (Wt[:, None] * Xh))[1, 1]  # WLS (3.30)
print("Part 5  slope coverage:", cov, " mean naive var %.4f vs true %.4f; Var WLS %.4f" % (v_naive.mean(), true_var, var_wls))
assert cov["naive"] < 0.92 and abs(cov["HC3"] - 0.95) < 0.015 and cov["HC0"] > cov["naive"] and var_wls < true_var / 5
assert abs(Bw[1].var() / var_wls - 1) < 0.08 and abs(Bh[1].mean() - 0.5) < 0.01 and abs(Bh[1].var() / true_var - 1) < 0.08  # both unbiased
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6)); ax[0].scatter(Xh @ Bh[:, 0], Eh[:, 0], s=12); ax[0].axhline(0, c="k", lw=.5); ax[0].set(title="residuals vs fitted: the fan", xlabel="fitted", ylabel="e_i")
ax[1].hist(np.sqrt(v_naive), 40, alpha=.6, label="naive se"); ax[1].hist(np.sqrt(v_hc3), 40, alpha=.6, label="HC3 se"); ax[1].axvline(np.sqrt(true_var), c="k", ls="--", label="true sd of slope"); ax[1].legend(fontsize=7); ax[1].set(title="slope standard errors, 4000 datasets")
fig.tight_layout(); fig.savefig(FIG / "heteroscedasticity.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 6: multicollinearity: VIF and the SVD variance decomposition (notes 3.7, Example 4.5)

# %%
rng = np.random.default_rng(6)
rho = 0.99; Xc = np.linalg.qr(rng.standard_normal((100, 2)))[0] @ np.linalg.cholesky(np.array([[1, rho], [rho, 1]])).T  # X'X = [[1, rho], [rho, 1]]
assert np.allclose(Xc.T @ Xc, [[1, rho], [rho, 1]])
U, s, Vt = np.linalg.svd(Xc, full_matrices=False); Vcov = np.linalg.inv(Xc.T @ Xc)  # Var(b)/sigma^2
assert close(s**2, [1.99, 0.01], 1e-12) and close(np.linalg.cond(Xc), np.sqrt(199), 1e-9)
vif = 1 / (1 - rho**2)
assert np.allclose(Vcov, sum(np.outer(Vt[i], Vt[i]) / s[i] ** 2 for i in range(2))) and close(Vcov[0, 0], vif, 1e-9)  # (3.32)
assert close([vif, np.sqrt(vif), Vcov[0, 1] / Vcov[0, 0]], [50.2513, 7.0888, -rho], 5e-5)
vsum, vdiff = np.array([1, 1]) @ Vcov @ [1, 1], np.array([1, -1]) @ Vcov @ [1, -1]
assert close([vsum, vdiff], [2 / (1 + rho), 2 / (1 - rho)], 1e-9) and close([vsum, vdiff, 0.5 / s[0] ** 2, 0.5 / s[1] ** 2], [1.0050, 200.0, 0.2513, 50.0])
Yc = (Xc @ [1.0, 1.0])[:, None] + rng.standard_normal((100, 20_000)); Bc = np.linalg.solve(Xc.T @ Xc, Xc.T @ Yc)
assert abs(Bc[0].var() / vif - 1) < 0.03 and abs(np.corrcoef(Bc)[0, 1] + rho) < 0.005 and abs(Bc.mean(1) - 1).max() < 0.1  # unbiased, huge variance
print("Example 4.5  s =", s.round(4), " VIF = %.4f  se inflation = %.4f  Var(b1+b2)/s2 = %.4f  Var(b1-b2)/s2 = %.1f" % (vif, np.sqrt(vif), vsum, vdiff))
rhos = np.array([0, .5, .9, .99, .999]); kap, var1 = [], []
for r_ in rhos:
    Xr = np.linalg.qr(rng.standard_normal((100, 2)))[0] @ np.linalg.cholesky(np.array([[1, r_], [r_, 1]])).T
    kap.append(np.linalg.cond(Xr)); var1.append(np.linalg.inv(Xr.T @ Xr)[0, 0])
assert np.allclose(var1, 1 / (1 - rhos**2)) and np.allclose(np.array(kap) ** 2, (1 + rhos) / (1 - rhos))  # Var(b_1)/sigma^2 = VIF; kappa^2 = (1+rho)/(1-rho)
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
ax[0].loglog(np.array(kap) ** 2, var1, "o-"); ax[0].set(xlabel="kappa(X)^2", ylabel="Var(b_1)/sigma^2", title="variance vs condition number")
th = np.linspace(0, 2 * np.pi, 200); ell = np.c_[np.cos(th), np.sin(th)] @ (Vt.T / s).T  # {b : b'X'Xb = 1}: semi-axes 1/sigma_i along v_i
ax[1].plot(*ell.T); ax[1].set(aspect="equal", title="confidence ellipse, rho = 0.99", xlabel="b_1", ylabel="b_2")
fig.tight_layout(); fig.savefig(FIG / "collinearity.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 7: model selection: C_p, AIC, BIC, LOO shortcut, post-selection (notes 3.8)

# %%
rng = np.random.default_rng(7)
f_true = lambda x: np.sin(5 * x + 1); dmax, sig7 = 16, 0.5
leg = lambda x, m: np.polynomial.legendre.legvander(x, m - 1)  # columns P_0..P_{m-1}: a well-conditioned polynomial basis


def nested_fits(x, y, dmax):
    """For degrees 0..dmax: RSS_d, leverages h_ii(d), LOO errors (3.25) -- all from ONE QR of the nested Legendre design."""
    Q = np.linalg.qr(leg(x, dmax + 1))[0]; E = y[:, None] - np.cumsum(Q * (Q.T @ y), 1); Hd = np.cumsum(Q**2, 1)
    return (E**2).sum(0), Hd, E / (1 - Hd)


n7 = 30; x7 = np.linspace(-1, 1, n7); mu7 = f_true(x7); Qf = np.linalg.qr(leg(x7, dmax + 1))[0]; pp = np.arange(1, dmax + 2)
bias2 = np.array([((mu7 - (Qf[:, :d + 1] * (Qf[:, :d + 1].T @ mu7)).sum(1)) ** 2).sum() for d in range(dmax + 1)])  # ||(I-H)mu||^2
Rm = 4000; Y7 = mu7[:, None] + sig7 * rng.standard_normal((n7, Rm)); Ynew = mu7[:, None] + sig7 * rng.standard_normal((n7, Rm))
rss_mc, epe_mc = np.zeros(dmax + 1), np.zeros(dmax + 1)
for d in range(dmax + 1):
    Fit = Qf[:, :d + 1] @ (Qf[:, :d + 1].T @ Y7); rss_mc[d] = ((Y7 - Fit) ** 2).sum(0).mean(); epe_mc[d] = ((Ynew - Fit) ** 2).sum(0).mean()
epe = n7 * sig7**2 + bias2 + pp * sig7**2                                           # (3.33)
assert np.allclose(rss_mc, bias2 + (n7 - pp) * sig7**2, rtol=0.03, atol=0.02) and np.allclose(epe_mc, epe, rtol=0.03)
assert np.allclose(rss_mc + 2 * pp * sig7**2, epe, rtol=0.03)                       # optimism 2 p sigma^2 (3.34), C_p unbiased (3.35)
print("Part 7  expected prediction error by degree:", epe.round(2), " best degree", epe.argmin(), "(n = 30)")
yq = Y7[:, 0]; d_ = 4; Xd = leg(x7, d_ + 1)                                          # LOO shortcut = brute force
brute = np.array([yq[i] - leg(x7[i:i + 1], d_ + 1) @ ols_qr(np.delete(Xd, i, 0), np.delete(yq, i)) for i in range(n7)]).ravel()
assert np.allclose(brute, nested_fits(x7, yq, dmax)[2][:, d_])
ns = [20, 40, 80, 160, 320, 640, 1280]; sel = {k: [] for k in ("AIC", "BIC", "LOO")}
for nn in ns:                                                                        # selected degree vs n, 100 datasets each
    xx_ = np.linspace(-1, 1, nn); picks = {k: [] for k in sel}
    for _ in range(100):
        rs, Hd, le = nested_fits(xx_, f_true(xx_) + sig7 * rng.standard_normal(nn), dmax)
        picks["AIC"].append(np.argmin(nn * np.log(rs / nn) + 2 * pp)); picks["BIC"].append(np.argmin(nn * np.log(rs / nn) + pp * np.log(nn)))  # (3.36)
        picks["LOO"].append(np.argmin((le**2).mean(0)))
    for k in sel: sel[k].append(np.mean(picks[k]))
print("Part 7  mean selected degree vs n =", ns, "\n        ", {k: np.round(v, 2) for k, v in sel.items()})
assert sel["LOO"][-1] > sel["LOO"][1] + 0.5 and sel["BIC"][-1] > sel["BIC"][1] + 0.5 and all(b <= a + 1e-9 for a, b in zip(sel["AIC"], sel["BIC"]))
assert sel["AIC"][0] > 12 > sel["LOO"][0]   # p close to n: AIC's penalty 2p is too light and it overfits wildly; LOO does not
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8)); dd = range(dmax + 1)
ax[0].plot(dd, bias2 / n7, label="bias^2 / n"); ax[0].plot(dd, pp * sig7**2 / n7, label="p sigma^2 / n (variance)"); ax[0].plot(dd, epe / n7 - sig7**2, "k", label="reducible error")
ax[0].plot(dd, rss_mc / n7, "--", label="training MSE"); ax[0].set(yscale="log", xlabel="degree", title="bias-variance, n = 30"); ax[0].legend(fontsize=7)
for k, v in sel.items(): ax[1].semilogx(ns, v, "o-", label=k)
ax[1].set(xlabel="n", ylabel="mean selected degree", title="selected degree vs n"); ax[1].legend()
fig.tight_layout(); fig.savefig(FIG / "model_selection.png", dpi=110); plt.close(fig)


def forward_aic(X, y, kmax=15):
    """Greedy forward selection: add the column most correlated with the residual; stop when AIC = n log(RSS/n) + 2k rises."""
    n_ = len(y); chosen, aic_old, e_ = [], n_ * np.log(y @ y / n_), y.copy()
    for k in range(1, kmax + 1):
        j = int(np.argmax([abs(X[:, c] @ e_) / np.linalg.norm(X[:, c]) if c not in chosen else -1 for c in range(X.shape[1])]))
        Xc_ = X[:, chosen + [j]]; e_new = y - Xc_ @ np.linalg.lstsq(Xc_, y, rcond=None)[0]; aic = n_ * np.log(e_new @ e_new / n_) + 2 * k
        if aic >= aic_old: break
        chosen, aic_old, e_ = chosen + [j], aic, e_new
    return chosen


sizes, fr, r2s = [], [], []                          # post-selection: y is pure noise, 50 candidate predictors, n = 100
for _ in range(200):
    Xn = rng.standard_normal((100, 50)); Xn -= Xn.mean(0); yn = rng.standard_normal(100); yn -= yn.mean(); ch = forward_aic(Xn, yn); sizes.append(len(ch))
    if ch: ft = ols_fit(np.c_[np.ones(100), Xn[:, ch]], yn); fr.append((ft["p"][1:] < 0.05).mean()); r2s.append(ft["R2"])
print("Part 7  forward-AIC on pure noise: mean model size %.1f, R^2 %.2f, share of selected coefficients with p < 0.05: %.2f" % (np.mean(sizes), np.mean(r2s), np.mean(fr)))
assert np.mean(sizes) > 5 and np.mean(fr) > 0.4 and np.mean(r2s) > 0.2   # nominal: 0, 0.05, 0

# %% [markdown]
# ## Part 8: ridge (notes 3.9, Example 4.6)

# %%
rng = np.random.default_rng(8)
Xr = make_design(30, np.array([10, 1, 0.1]), rng); U, s, Vt = np.linalg.svd(Xr, full_matrices=False); lam = 1.0
assert close(s, [10, 1, 0.1], 1e-12); fac = s**2 / (s**2 + lam)
assert close(fac, [0.990099, 0.5, 0.009901], 5e-7) and close(fac.sum(), 1.5, 1e-12)                  # Example 4.6: shrinkage factors, (3.41)
assert close(s / (s**2 + lam), [0.0990, 0.5, 0.0990], 5e-5) and close(1 / s, [0.1, 1, 10])           # ridge vs OLS multiplier on u_i'y
yr = Xr @ [1, -1, 2.0] + rng.standard_normal(30)
b_dir = np.linalg.solve(Xr.T @ Xr + lam * np.eye(3), Xr.T @ yr)                                       # (3.39)
b_svd = Vt.T @ (s / (s**2 + lam) * (U.T @ yr))                                                        # (3.40)
b_aug = np.linalg.lstsq(np.r_[Xr, np.sqrt(lam) * np.eye(3)], np.r_[yr, np.zeros(3)], rcond=None)[0]  # ridge = OLS on augmented data
assert np.allclose(b_dir, b_svd) and np.allclose(b_dir, b_aug) and np.isclose(hat_diag(Xr, lam).sum(), 1.5)
alpha, sg = np.array([1.0, 1, 1]), 1.0                                                                # a = V'beta, noise variance 1
mse_ridge = lambda l: np.sum((sg**2 * s**2 + l**2 * alpha**2) / (s**2 + l) ** 2)                      # (3.42)
assert close([mse_ridge(0.0), mse_ridge(1.0)], [101.01, 1.5], 5e-3)
assert mse_ridge(1e-3) < mse_ridge(0.0) and np.isclose((mse_ridge(1e-9) - mse_ridge(0.0)) / 1e-9, -2 * np.sum(1 / s**4), rtol=1e-3)  # MSE'(0) < 0
Bt = np.linalg.solve(Xr.T @ Xr + lam * np.eye(3), Xr.T @ ((Xr @ (Vt.T @ alpha))[:, None] + rng.standard_normal((30, 4000))))
mse_mc = np.mean(np.sum((Bt - (Vt.T @ alpha)[:, None]) ** 2, 0)); assert abs(mse_mc / 1.5 - 1) < 0.05
print("Example 4.6  factors", fac.round(6), " df =", fac.sum(), " MSE(OLS)=%.2f  MSE(ridge, lam=1)=%.4f (MC %.3f)" % (mse_ridge(0.0), mse_ridge(1.0), mse_mc))
pb, bb = 5, np.array([1.0, 0.5, 0.5, 0.2, 0.1]); l_star = pb * 0.25 / (bb @ bb); mse_o = lambda l: (pb * 0.25 + l**2 * (bb @ bb)) / (1 + l) ** 2  # orthonormal design
assert close(mse_o(l_star), pb * 0.25 / (1 + l_star), 1e-12) and mse_o(l_star) < mse_o(0.0) and abs(optimize.minimize_scalar(mse_o, bounds=(0, 10), method="bounded").x - l_star) < 1e-4
n, P = 40, 10                                         # correlated design: coefficient path, LOO-CV via the shortcut, MSE vs lambda
Xp = rng.standard_normal((n, P)) @ np.linalg.cholesky(0.8 ** np.abs(np.subtract.outer(range(P), range(P)))).T; Xp = (Xp - Xp.mean(0)) / Xp.std(0)
bt = np.r_[3, -2, 1.5, np.zeros(7)]; yp = Xp @ bt + 3 * rng.standard_normal(n); yp -= yp.mean()
U, s, Vt = np.linalg.svd(Xp, full_matrices=False); lams = np.logspace(-2, 3, 60)
path = np.array([Vt.T @ (s / (s**2 + l) * (U.T @ yp)) for l in lams])
loo = np.array([np.mean(((yp - Xp @ pth) / (1 - hat_diag(Xp, l))) ** 2) for pth, l in zip(path, lams)])
i0 = 20; Xm = lambda i: np.delete(Xp, i, 0)
brute = np.array([yp[i] - Xp[i] @ np.linalg.solve(Xm(i).T @ Xm(i) + lams[i0] * np.eye(P), Xm(i).T @ np.delete(yp, i)) for i in range(n)])
assert np.allclose(np.mean(brute**2), loo[i0])        # the LOO shortcut is exact for ridge
l_cv = lams[loo.argmin()]; al = Vt @ bt; mse_curve = np.array([np.sum((9 * s**2 + l**2 * al**2) / (s**2 + l) ** 2) for l in lams])  # sigma^2 = 9
assert mse_curve.min() < 0.7 * np.sum(9 / s**2) and 0 < mse_curve.argmin() < len(lams) - 1
mse_cv = np.interp(np.log(l_cv), np.log(lams), mse_curve); assert mse_cv < np.sum(9 / s**2)
print("Part 8  OLS MSE %.2f; ridge min MSE %.2f at lambda %.2f; LOO-CV picks lambda %.2f (MSE there %.2f)" % (np.sum(9 / s**2), mse_curve.min(), lams[mse_curve.argmin()], l_cv, mse_cv))
Xw = rng.standard_normal((10, 25)); yw = rng.standard_normal(10); b_mn = np.linalg.pinv(Xw) @ yw    # p > n: lambda -> 0 is the minimum-norm solution
assert np.allclose(Xw.T @ np.linalg.solve(Xw @ Xw.T + 1e-9 * np.eye(10), yw), b_mn, atol=1e-6) and np.allclose(Xw @ b_mn, yw)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].semilogx(lams, path); ax[0].axvline(l_cv, c="k", ls=":"); ax[0].set(xlabel="lambda", ylabel="coefficient", title="ridge path (dotted: LOO-CV)")
ax[1].loglog(lams, mse_curve); ax[1].axhline(np.sum(9 / s**2), c="r", ls="--", label="OLS"); ax[1].set(xlabel="lambda", ylabel="E||b - beta||^2", title="ridge MSE: minimum below OLS"); ax[1].legend()
fig.tight_layout(); fig.savefig(FIG / "ridge.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 9: the lasso by coordinate descent (notes 3.10, Example 4.7)

# %%
def soft(z, t):
    """S_t(z) = sign(z) (|z| - t)_+, the minimizer of 0.5 (b - z)^2 + t |b| (3.44)."""
    return np.sign(z) * np.maximum(np.abs(z) - t, 0.0)


def lasso_cd(X, y, lam, b=None, tol=1e-9, max_sweeps=5000, l2=0.0):
    """Minimize 0.5||y - Xb||^2 + lam ||b||_1 (+ 0.5 l2 ||b||^2: elastic net) by cycling b_j <- S_lam(x_j'r_j)/(||x_j||^2 + l2) (3.46)."""
    p = X.shape[1]; b = np.zeros(p) if b is None else b.copy(); r = y - X @ b; cn = (X**2).sum(0)
    for _ in range(max_sweeps):
        delta = 0.0
        for j in range(p):
            bj = soft(X[:, j] @ r + cn[j] * b[j], lam) / (cn[j] + l2)
            if bj != b[j]: r += X[:, j] * (b[j] - bj); delta = max(delta, abs(bj - b[j])); b[j] = bj
        if delta < tol: break
    return b


rng = np.random.default_rng(9)
Qo = np.linalg.qr(rng.standard_normal((20, 5)))[0]; z = np.array([3.0, -1.5, 0.4, -0.2, 2.2])  # Example 4.7: orthonormal design
w = rng.standard_normal(20); w -= Qo @ (Qo.T @ w); yo = Qo @ z + w; assert np.allclose(Qo.T @ yo, z)
bl = lasso_cd(Qo, yo, 1.0)
assert close(bl, [2.0, -0.5, 0, 0, 1.2], 1e-7) and np.allclose(bl, soft(z, 1.0)) and close(z / 2, [1.5, -0.75, 0.2, -0.1, 1.1])
print("Example 4.7  soft(z, 1) =", soft(z, 1.0), " ridge(lambda=1) =", z / 2)
n, P = 30, 6                                          # slow reference: split b = b+ - b-, a smooth bound-constrained problem for L-BFGS-B
Xs6 = rng.standard_normal((n, P)); Xs6[:, 1] = Xs6[:, 0] + 0.3 * rng.standard_normal(n); ys6 = Xs6 @ [2, 0, -1, 0, 0, 0.5] + rng.standard_normal(n); lam6 = 8.0
def obj(v): bv = v[:P] - v[P:]; r_ = ys6 - Xs6 @ bv; return 0.5 * r_ @ r_ + lam6 * v.sum(), np.r_[-Xs6.T @ r_ + lam6, Xs6.T @ r_ + lam6]
ref = optimize.minimize(obj, np.zeros(2 * P), jac=True, bounds=[(0, None)] * (2 * P), method="L-BFGS-B", options=dict(ftol=1e-15, gtol=1e-11, maxiter=10_000))
b_ref = ref.x[:P] - ref.x[P:]; b_cd = lasso_cd(Xs6, ys6, lam6); g = Xs6.T @ (ys6 - Xs6 @ b_cd); lmax = np.abs(Xs6.T @ ys6).max()
assert np.allclose(b_cd, b_ref, atol=1e-5) and np.allclose(g[b_cd != 0], lam6 * np.sign(b_cd[b_cd != 0]), atol=1e-6) and np.all(np.abs(g[b_cd == 0]) <= lam6 + 1e-9)  # KKT (3.45)
assert not lasso_cd(Xs6, ys6, lmax * 1.0001).any() and (lasso_cd(Xs6, ys6, lmax * 0.99) != 0).sum() == 1
print("Part 9  CD vs L-BFGS-B reference: max diff %.1e; KKT holds; lambda_max = %.3f; nonzeros at lam=8: %d" % (np.abs(b_cd - b_ref).max(), lmax, (b_cd != 0).sum()))


def lasso_path(X, y, lams, l2=0.0):
    """Warm-started path over a decreasing grid of lambdas (l2 > 0 gives the elastic-net path, exercise C4)."""
    b, out = None, []
    for l in lams: b = lasso_cd(X, y, l, b, tol=1e-7, l2=l2); out.append(b)
    return np.array(out)


def cv_lasso(X, y, lams, K=5):
    idx = np.arange(len(y)) % K; err = np.zeros(len(lams))
    for k in range(K):
        pth = lasso_path(X[idx != k], y[idx != k], lams); err += ((y[idx == k][:, None] - X[idx == k] @ pth.T) ** 2).sum(0)
    return err / len(y)


n, P = 100, 50                                        # sparse recovery: five nonzero coefficients among 50
Xl = rng.standard_normal((n, P)); Xl = (Xl - Xl.mean(0)) / Xl.std(0); bs = np.zeros(P); bs[[3, 11, 20, 31, 44]] = [3, -2, 2.5, -1.5, 2]
yl = Xl @ bs + rng.standard_normal(n); yl -= yl.mean(); lams9 = np.geomspace(np.abs(Xl.T @ yl).max(), 1e-2 * np.abs(Xl.T @ yl).max(), 40)
pth = lasso_path(Xl, yl, lams9); ib = cv_lasso(Xl, yl, lams9).argmin(); b_l = pth[ib]; b_ols = ols_qr(Xl, yl)
supp, true = set(np.flatnonzero(b_l)), set(np.flatnonzero(bs))
print("Part 9  CV-lasso: lambda = %.2f, support size %d, true support recovered: %s, false positives %d; ||b - beta||: lasso %.2f, OLS %.2f" % (
    lams9[ib], len(supp), true <= supp, len(supp - true), np.linalg.norm(b_l - bs), np.linalg.norm(b_ols - bs)))
assert true <= supp and np.linalg.norm(b_l - bs) < 0.5 * np.linalg.norm(b_ols - bs)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8)); zz = np.linspace(-3, 3, 200)
ax[0].semilogx(lams9, pth, lw=1); ax[0].axvline(lams9[ib], c="k", ls=":"); ax[0].set(xlabel="lambda", ylabel="coefficient", title="lasso path, p = 50 (dotted: CV)")
ax[1].plot(zz, zz, "k:", label="OLS"); ax[1].plot(zz, soft(zz, 1), label="lasso (soft)"); ax[1].plot(zz, zz / 2, label="ridge"); ax[1].legend(); ax[1].set(xlabel="z", title="orthonormal design, lambda = 1")
fig.tight_layout(); fig.savefig(FIG / "lasso.png", dpi=110); plt.close(fig)


def compare(kind, reps=12):                           # ridge vs lasso, lambda chosen on a validation set
    out = {"ridge": [], "lasso": [], "missed": [], "false+": []}
    for _ in range(reps):
        Sg = 0.9 ** np.abs(np.subtract.outer(range(P), range(P))) if "correlated" in kind else np.eye(P)
        Xa = rng.standard_normal((300, P)) @ np.linalg.cholesky(Sg).T; Xa = (Xa - Xa.mean(0)) / Xa.std(0)
        bt = 0.45 * rng.standard_normal(P) if kind.startswith("dense") else bs.copy()
        ya = Xa @ bt + 2 * rng.standard_normal(300); ya -= ya.mean(); Xt, yt_, Xv, yv = Xa[:100], ya[:100], Xa[100:200], ya[100:200]
        U_, s_, Vt_ = np.linalg.svd(Xt, full_matrices=False)
        br = [Vt_.T @ (s_ / (s_**2 + l) * (U_.T @ yt_)) for l in np.geomspace(1e-1, 1e4, 40)]
        out["ridge"].append(np.linalg.norm(br[int(np.argmin([np.sum((yv - Xv @ b_) ** 2) for b_ in br]))] - bt))
        bl_ = lasso_path(Xt, yt_, np.geomspace(np.abs(Xt.T @ yt_).max(), 1e-2 * np.abs(Xt.T @ yt_).max(), 20))
        bhat = bl_[int(np.argmin(((yv[:, None] - Xv @ bl_.T) ** 2).sum(0)))]
        out["lasso"].append(np.linalg.norm(bhat - bt)); out["missed"].append(np.sum((bhat == 0) & (bs != 0))); out["false+"].append(np.sum((bhat != 0) & (bs == 0)))
    return {k: float(np.mean(v)) for k, v in out.items()}


res = {k: compare(k) for k in ("sparse", "dense", "sparse-correlated", "dense-correlated")}
print("Part 9  ||b - beta|| (mean of 12 runs) and lasso support errors:")
for k, d in res.items(): print("        %-18s" % k, {m: round(v, 2) for m, v in d.items()})
assert res["sparse"]["lasso"] < 0.7 * res["sparse"]["ridge"] and res["dense"]["ridge"] < res["dense"]["lasso"]

# %% [markdown]
# ## Part 10: numbers used in exercises and solutions (every figure in exercises.md and solutions.md is asserted here)

# %%
xa, ya = np.array([0.0, 1, 2, 3]), np.array([1.0, 2, 4, 5]); Xa = np.c_[np.ones(4), xa]; fa = ols_fit(Xa, ya)     # A1, A2
assert close(fa["b"], [0.9, 1.4]) and close(fa["e"], [0.1, -0.3, 0.3, -0.1]) and close([fa["rss"], fa["tss"], fa["R2"]], [0.2, 10, 0.98])
assert close(fa["h"], [0.7, 0.3, 0.3, 0.7]) and close(fa["e"] / (1 - fa["h"]), [0.3333, -0.4286, 0.4286, -0.3333], 5e-5) and close(Xa.T @ fa["e"], [0, 0], 1e-12)
press = ((fa["e"] / (1 - fa["h"])) ** 2).sum(); assert close([press, press / 4], [0.5896, 0.1474], 5e-5)
assert (((xa - 1.5) ** 2).sum(), ((xa - 1.5) * (ya - 3)).sum()) == (5.0, 7.0) and close(1.4**2 * 5, 9.8) and close(7**2 / (5 * 10), fa["R2"])  # A1
bm1 = ols_qr(Xa[1:], ya[1:]); assert close(bm1, [2 / 3, 1.5], 1e-10) and close(ya[0] - Xa[0] @ bm1, 1 / 3, 1e-10)                                # A2(b): refit without obs 1
assert close(fa["e"] / (1 - fa["h"]) / fa["e"], [10 / 3, 10 / 7, 10 / 7, 10 / 3], 1e-10)                                                          # A2(c): growth factors 1/(1-h)
s2_ = 45 / 17; t17 = stats.t.ppf(0.975, 17)                                                                          # A3
assert close([s2_, np.sqrt(s2_), t17, t17 * np.sqrt(s2_ * 0.1), t17 * np.sqrt(s2_ * 1.1)], [2.6471, 1.6270, 2.1098, 1.0855, 3.6002], 5e-5)
za = np.array([2.5, -0.7, 1.1, -3.2]); assert close(soft(za, 1.0), [1.5, 0, 0.1, -2.2]) and close(za / 2, [1.25, -0.35, 0.55, -1.6])  # A4
sa = np.array([4, 2, 0.5]); assert close(sa**2 / (sa**2 + 1), [0.9412, 0.8, 0.2], 5e-5) and close((sa**2 / (sa**2 + 1)).sum(), 1.9412, 5e-5)  # A5
assert close([1 / 0.5**2, 0.5**2 / (0.5**2 + 1) ** 2], [4.0, 0.16], 1e-12) and close(-1 / (0.5**2 + 1), -0.8, 1e-12) and close(4.0 / 0.16, 25, 1e-12)
bo = np.array([1.0, 0.5, 0.5, 0.2, 0.1]); lo_ = 5 * 0.25 / (bo @ bo); mo = lambda l: (5 * 0.25 + l**2 * (bo @ bo)) / (1 + l) ** 2  # A5(c): orthonormal design, (3.42)
assert close([bo @ bo, lo_, mo(0.0), mo(lo_)], [1.55, 0.8065, 1.25, 0.6920], 5e-5) and np.isclose(mo(lo_), 1.25 / (1 + lo_))
Xd = np.array([[1, 1, 0], [1, 1, 0], [1, 0, 1], [1, 0, 1.0]]); yd = np.array([3.0, 5, 6, 8]); bd = np.linalg.pinv(Xd) @ yd                              # A6
assert np.linalg.matrix_rank(Xd) == 2 and np.allclose(Xd @ [1, -1, -1], 0) and close(bd, [3.6667, 0.3333, 3.3333], 5e-5) and close(Xd @ bd, [4, 4, 7, 7])
assert close(np.linalg.norm(bd), 4.9666, 5e-5) and close(ols_qr(np.c_[np.ones(4), [0, 0, 1, 1]], yd), [4.0, 3.0], 1e-12)
zc = np.array([0.0, 1, -1, 0, 0]); f2 = ols_fit(np.c_[X41, zc], y41)                                                 # B2
assert close([f2["R2"], f2["adjR2"]], [0.6140, 0.2281]) and f41["R2"] < f2["R2"] and f2["adjR2"] < f41["adjR2"]
zt = zc - X41 @ ols_qr(X41, zc); dr = f41["rss"] - f2["rss"]; tz = f2["t"][-1]                                          # B2(c),(d): RSS_old - RSS_new = (z~'e)^2/||z~||^2 = t^2 RSS_new/(n-p-1)
assert close([f2["rss"], f2["s2"], zt @ zt, zt @ f41["e"], dr], [2.3158, 1.1579, 1.9, -0.4, 0.0842], 5e-5) and close(dr, (zt @ f41["e"]) ** 2 / (zt @ zt), 1e-12)
assert close([tz**2, abs(tz)], [0.0727, 0.2697], 5e-5) and close(tz**2, dr / (f2["rss"] / 2), 1e-12) and abs(tz) < 1 and close(zt, [-0.2, 0.9, -1, 0.1, 0.2], 1e-12)                # |t| < 1  <=>  adjusted R^2 falls
assert close(np.corrcoef(y41, X41 @ f41["b"])[0, 1] ** 2, f41["R2"], 1e-12)
bl_, bs_ = ols_qr(X42[:, :3], y42), ols_qr(X42[:, :2], y42); gam = ols_qr(X42[:, :2], x2)[1]                           # B3
assert close(bs_[1], bl_[1] + gam * bl_[2], 1e-12) and close([bs_[1], bl_[1], gam, bl_[2]], [1.0440, 0.3103, 0.8214, 0.8933], 5e-5)
assert close(30 * np.log(1 + 2 * 3.5 / 26), 7.1523) and close(stats.f.sf(3.5, 2, 26), 0.0451) and close(stats.chi2.sf(7.1523, 2), 0.0280)  # B5
x5 = 1 + 2 * 3.5 / 26; assert close([x5 ** -13, x5 ** -15], [stats.f.sf(3.5, 2, 26), stats.chi2.sf(30 * np.log(x5), 2)], 1e-12) and close(x5 ** -2, 0.0280 / 0.0451, 2e-3)  # q = 2: closed forms, exponents (n-p)/2 and n/2
yb = np.array([1.0, 2, 3, 6]); eb = yb - yb.mean()                                                                    # B6
assert close([eb @ eb / 3 / 4, eb @ eb / 16, (eb**2 / 0.5625).sum() / 16], [1.1667, 0.875, 1.5556], 5e-5)
assert eb @ eb == 14 and close(eb @ eb / 16 * 4 / 3, eb @ eb / 3 / 4, 1e-12) and close((eb**2 / 0.5625).sum() / 16, 14 / 9, 1e-12)   # HC1 = naive; HC3 = sum e^2/(n-1)^2
x4 = np.arange(1.0, 6); a4 = np.array([-0.25, 0, 0, 0, 0.25]); X4 = np.c_[np.ones(5), x4]; H4 = X4 @ np.linalg.solve(X4.T @ X4, X4.T)   # B4: (y5 - y1)/4 vs OLS slope
w4 = (x4 - 3) / 10; d4 = a4 - w4
assert close(a4 @ X4, [0, 1], 1e-12) and np.allclose(H4 @ a4, w4) and close(d4, [-0.05, 0.1, 0, -0.1, 0.05], 1e-12) and np.allclose((np.eye(5) - H4) @ a4, d4)
assert close([a4 @ a4, w4 @ w4, d4 @ d4], [0.125, 0.1, 0.025], 1e-12) and close(a4 @ a4, w4 @ w4 + d4 @ d4, 1e-12)
G7 = np.array([[1, 0.9], [0.9, 1]]); ev7 = np.linalg.eigvalsh(G7); Gi7 = np.linalg.inv(G7); lam7 = 0.1                           # B7: rho = 0.9, unit-length columns
assert close(ev7, [0.1, 1.9], 1e-12) and close(np.sqrt(ev7[1] / ev7[0]), 4.3589, 5e-5) and close([Gi7[0, 0], 1 / (1 - 0.81), 0.5 / 1.9 + 0.5 / 0.1], [5.2632] * 3, 5e-5)
assert close([np.ones(2) @ Gi7 @ np.ones(2), np.array([1, -1]) @ Gi7 @ [1, -1], Gi7[0, 1] / Gi7[0, 0]], [1.0526, 20.0, -0.9], 5e-5)
assert close([ev7[0] / (ev7[0] + lam7) ** 2, 1 / ev7[0], ev7[0] / (ev7[0] + lam7), -lam7 / (ev7[0] + lam7), ev7[1] / (ev7[1] + lam7)], [2.5, 10, 0.5, -0.5, 0.95], 1e-12) and close((10 - 2.5) / 0.25, 30)
xc = np.array([1.0, 2, 3, 4]); Xc2 = np.c_[xc, xc]; yc2 = np.array([1.5, 2.5, 3.0, 4.5])                              # C2, C4
b_a, b_b = lasso_cd(Xc2, yc2, 1.0, np.array([5.0, 0.0])), lasso_cd(Xc2, yc2, 1.0, np.array([0.0, 5.0])); lo = lambda b: 0.5 * np.sum((yc2 - Xc2 @ b) ** 2) + np.abs(b).sum()
rb = np.linalg.solve(Xc2.T @ Xc2 + np.eye(2), Xc2.T @ yc2); en = lasso_cd(Xc2, yc2, 1.0, np.array([5.0, 0.0]), l2=1.0); en2 = lasso_cd(Xc2, yc2, 1.0, np.array([0.0, 5.0]), l2=1.0)
assert np.isclose(b_a.sum(), b_b.sum()) and np.isclose(lo(b_a), lo(b_b)) and not np.allclose(b_a, b_b) and np.isclose(rb[0], rb[1]) and np.allclose(en, en2, atol=1e-7) and np.isclose(en[0], en[1])
assert close([b_a.sum(), rb[0], en[0]], [1.0833, 0.5492, 0.5328], 5e-5) and close(b_a, [1.0833, 0], 5e-5)
assert close(Xc2.T @ yc2, [33.5, 33.5]) and close((Xc2**2).sum(0), [30, 30]) and close([32.5 / 30, 33.5 / 61, 32.5 / 61], [b_a.sum(), rb[0], en[0]], 1e-7) and close(b_b, [0, 1.0833], 5e-5)  # C2, C4 by hand
Xm = np.c_[np.c_[np.ones(5), xl], np.eye(5)[:, 4]]; yl4 = np.array([1.1, 1.9, 3.2, 3.9, 3.0]); bm = np.linalg.solve(Xm.T @ Xm, Xm.T @ yl4)   # C3: mean-shift model
em = yl4 - Xm @ bm; tm = bm[2] / np.sqrt(em @ em / 2 * np.linalg.inv(Xm.T @ Xm)[2, 2]); h4, r4, t4, D4, f4 = influence(np.c_[np.ones(5), xl], yl4)
assert close(bm[2], f4["e"][4] / (1 - h4[4]), 1e-10) and close(tm, t4[4], 1e-10) and close([bm[2], tm], [-6.8, -10.8367])
b4 = ols_qr(np.c_[np.ones(4), xl[:4]], yl4[:4]); r4_ = yl4[:4] - np.c_[np.ones(4), xl[:4]] @ b4                                  # C3: the 4-point fit that the dummy leaves behind
assert close(b4, [0.1, 0.97], 1e-10) and close([r4_ @ r4_, em @ em, f4["e"][4], h4[4], b4 @ [1, 10]], [0.063, 0.063, -0.544, 0.92, 9.8], 5e-5)
assert close(tm, f4["e"][4] / np.sqrt(r4_ @ r4_ / 2 * (1 - h4[4])), 1e-10) and close(2 * stats.t.sf(abs(tm), 2), 0.0084, 5e-5) and close(np.sqrt(r4_ @ r4_ / 2), 0.17748, 5e-6)                                                   # t = e_i / (s_(i) sqrt(1 - h_ii))
xsm, ysm, gsm = np.array([1, 2, 3, 6, 7, 8.0]), np.array([8, 9, 10, 1, 2, 3.0]), np.array([0, 0, 0, 1, 1, 1.0])      # Simpson (notes 3.11)
assert close(ols_qr(np.c_[np.ones(6), xsm], ysm), [10.759, -1.1687]) and close(ols_qr(np.c_[np.ones(6), xsm, gsm], ysm), [7, 1, -12], 1e-9)
print("Part 10 exercise numbers verified.")

# %% [markdown]
# ## Your turn
#
# 1. **WLS** (notes 3.6): Part 5 compares OLS and WLS with the true variances. Estimate them instead (fit OLS, regress log e_i^2 on log x_i,
#    use the fitted variances as weights). Does the coverage of the slope interval survive?
# 2. **Interaction and FWL** (exercise B8): reference solution below.
# 3. **Elastic net** (exercise C4): `lasso_cd(..., l2=...)` already has it; reference solution below on a grouped design.
# 4. **Omitted-variable bias** (notes 3.8): reference solution below.
# 5. **Stepwise** (exercise C1): `forward_aic` in Part 7; change the stopping rule to BIC and to a p-value threshold and compare;
#    reference solution below, with a sample-splitting repair.

# %%
rng = np.random.default_rng(10)                       # Reference solution 1 (B8): y ~ 1 + x + g + x*g
n = 80; xg = rng.uniform(0, 5, n); gg = (rng.uniform(size=n) < 0.5).astype(float); yg = 1 + 0.5 * xg + 2 * gg + 0.8 * xg * gg + rng.standard_normal(n)
Xg = np.c_[np.ones(n), xg, gg, xg * gg]; bg = ols_qr(Xg, yg); Mg = np.eye(n) - Xg[:, :3] @ np.linalg.solve(Xg[:, :3].T @ Xg[:, :3], Xg[:, :3].T)
assert np.isclose((Mg @ Xg[:, 3]) @ (Mg @ yg) / ((Mg @ Xg[:, 3]) @ (Mg @ Xg[:, 3])), bg[3])   # FWL slope on the residualized product column
sl0, sl1 = ols_qr(Xg[gg == 0][:, :2], yg[gg == 0])[1], ols_qr(Xg[gg == 1][:, :2], yg[gg == 1])[1]
assert np.isclose(bg[1], sl0) and np.isclose(bg[1] + bg[3], sl1)                               # interaction = difference of within-group slopes
print("B8  interaction coefficient %.4f = slope(g=1) - slope(g=0) = %.4f - %.4f = %.4f (truth 0.8)" % (bg[3], sl1, sl0, sl1 - sl0))
n, Rr = 50, 20_000                                    # Reference solution 2: omitted-variable bias, E[b_short] = b_1 + gamma b_2 (3.37)
xa_ = rng.standard_normal(n); xb_ = 0.6 * xa_ + 0.8 * rng.standard_normal(n); xb_ -= xb_.mean(); xa_ -= xa_.mean()
Xfull = np.c_[np.ones(n), xa_, xb_]; bt = np.array([1.0, 2.0, 1.5]); Yo = (Xfull @ bt)[:, None] + rng.standard_normal((n, Rr))
b_short = np.linalg.solve(Xfull[:, :2].T @ Xfull[:, :2], Xfull[:, :2].T @ Yo)[1]; gamma = (xa_ @ xb_) / (xa_ @ xa_)
print("OVB  mean of short-regression slope %.4f; theory b_1 + gamma b_2 = %.4f (gamma = %.4f)" % (b_short.mean(), bt[1] + gamma * bt[2], gamma))
assert abs(b_short.mean() - (bt[1] + gamma * bt[2])) < 0.01 and abs(gamma * bt[2]) > 0.5

# %%
def forward_select(X, y, rule="aic", kmax=15, alpha=0.05):
    """Reference solution to C1: greedy forward selection (most residual-correlated column first); stop when AIC (penalty 2k),
    BIC (penalty k log n) or the t-test p-value of the newly added coefficient (< alpha) says to."""
    n_ = len(y); chosen, old, e_ = [], n_ * np.log(y @ y / n_), y.copy()
    for k in range(1, kmax + 1):
        j = int(np.argmax([abs(X[:, c] @ e_) / np.linalg.norm(X[:, c]) if c not in chosen else -1 for c in range(X.shape[1])]))
        ft = ols_fit(np.c_[np.ones(n_), X[:, chosen + [j]]], y)
        if rule == "pval":
            if ft["p"][-1] >= alpha: break
        else:
            crit = n_ * np.log(ft["rss"] / n_) + (2 if rule == "aic" else np.log(n_)) * k
            if crit >= old: break
            old = crit
        chosen, e_ = chosen + [j], ft["e"]
    return chosen


rng = np.random.default_rng(11); C1 = {}                 # Reference solution (C1): pure-noise y, n = 100, 50 candidates, 200 datasets
for rule in ("aic", "bic", "pval"):
    sz, sg, r2, sg_split = [], [], [], []
    for _ in range(200):
        Xn = rng.standard_normal((100, 50)); Xn -= Xn.mean(0); yn = rng.standard_normal(100); yn -= yn.mean()
        ch = forward_select(Xn, yn, rule); sz.append(len(ch))
        if ch: ft = ols_fit(np.c_[np.ones(100), Xn[:, ch]], yn); sg.append((ft["p"][1:] < 0.05).mean()); r2.append(ft["R2"])
        ch2 = forward_select(Xn[:50], yn[:50], rule)     # repair: select on the first 50 rows ...
        if ch2: ft2 = ols_fit(np.c_[np.ones(50), Xn[50:, ch2]], yn[50:]); sg_split.append((ft2["p"][1:] < 0.05).mean())  # ... test on the other 50
    C1[rule] = (np.mean(sz), np.mean(r2), np.mean(sg), np.mean(sg_split))
    print("C1  %-4s size %.2f  R^2 %.3f  share p<0.05 %.3f  | select-then-test on separate halves: share p<0.05 %.3f" % ((rule,) + C1[rule]))
assert C1["bic"][0] < 3.5 < C1["aic"][0] and all(v[2] > 0.45 for v in C1.values()) and all(abs(v[3] - 0.05) < 0.03 for v in C1.values())
print("C1  P(at least one of 50 null predictors has p < 0.05) = 1 - 0.95^50 = %.4f" % (1 - 0.95**50)); assert close(1 - 0.95**50, 0.9231, 5e-5)
thr = [np.sqrt(97 * (np.exp(2 / 100) - 1)), np.sqrt(97 * (100 ** (1 / 100) - 1)), stats.t.ppf(0.975, 97)]; print("C1  entry thresholds on |t| (AIC, BIC, p<0.05):", np.round(thr, 3))
assert close(np.square(thr[0]), 1.96, 5e-3) and close(thr, [1.40, 2.14, 1.98], 6e-3)

rng = np.random.default_rng(13)                          # Reference solution (C4): grouped design, 3 groups x 5 correlated columns, beta = 1 on all 15
def grouped_design(n, rng):
    Z = rng.standard_normal((n, 3)); X = rng.standard_normal((n, 50))
    for g in range(3): X[:, 5 * g:5 * g + 5] = Z[:, [g]] + 0.35 * rng.standard_normal((n, 5))
    return (X - X.mean(0)) / X.std(0)
Xcorr = grouped_design(20_000, np.random.default_rng(1)); print("C4  within-group correlation %.3f" % np.corrcoef(Xcorr[:, :5].T)[0, 1]); assert abs(np.corrcoef(Xcorr[:, :5].T)[0, 1] - 0.89) < 0.01
bgrp = np.r_[np.ones(15), np.zeros(35)]; C4 = {m: [] for m in (0, 5, 20, 50)}
for _ in range(12):
    Xa = grouped_design(300, rng); ya = Xa @ bgrp + 2 * rng.standard_normal(300); ya -= ya.mean()
    Xt, yt, Xv, yv = Xa[:100], ya[:100], Xa[100:200], ya[100:200]; lm = np.abs(Xt.T @ yt).max(); lams = np.geomspace(lm, 1e-2 * lm, 20)
    for m in C4:                                         # mu = 0 is the lasso; lambda chosen on the validation half
        pth = lasso_path(Xt, yt, lams, l2=m); C4[m].append(np.linalg.norm(pth[int(np.argmin(((yv[:, None] - Xv @ pth.T) ** 2).sum(0)))] - bgrp))
C4 = {m: float(np.mean(v)) for m, v in C4.items()}
print("C4  mean ||b - beta|| over 12 runs by mu:", {m: round(v, 2) for m, v in C4.items()})
assert C4[20] < 0.6 * C4[0] and C4[5] < C4[0] and close([C4[0], C4[5], C4[20], C4[50]], [2.24, 1.52, 0.92, 0.84], 0.05)

# %%
if __name__ == "__main__":
    print("All checks passed in %.1f s." % (time.time() - T0))
