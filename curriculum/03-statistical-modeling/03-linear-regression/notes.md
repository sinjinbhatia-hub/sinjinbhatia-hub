# ST-3: Linear Regression: Least Squares as Projection, Gauss–Markov, Inference, Diagnostics, Model Selection, Multicollinearity, Ridge and Lasso

**Time:** ~14 hours · **Prerequisites:** ST-2, LA-5, LA-9 · **Unlocks:** ST-4, ST-7, ML-3

> **Conventions in this module.** $n$ observations; $\mathbf y\in\mathbb R^n$ is the response; $X$ is the $n\times p$ design matrix whose first column is $\mathbf 1$ (the intercept), so there are $p-1$ predictors; $\mathbf x_j$ is column $j$, $\mathbf x_i^\top$ is row $i$. Model $\mathbf y = X\boldsymbol\beta+\boldsymbol\varepsilon$ with $\mathbb E\boldsymbol\varepsilon=\mathbf 0$, $\operatorname{Var}\boldsymbol\varepsilon=\sigma^2 I$ ($\sigma$ with no subscript is the noise standard deviation; $\sigma_i$ with a subscript is always a singular value of $X$). $X$ is treated as fixed. Estimates wear hats: $\hat{\boldsymbol\beta}=(X^\top X)^{-1}X^\top\mathbf y$, hat matrix $H=X(X^\top X)^{-1}X^\top$, fitted values $\hat{\mathbf y}=H\mathbf y$, **residuals** $\mathbf e=\mathbf y-\hat{\mathbf y}=(I-H)\mathbf y$ (computed) versus **errors** $\boldsymbol\varepsilon$ (never observed). $\mathrm{RSS}=\|\mathbf e\|^2$, $\mathrm{TSS}=\|\mathbf y-\bar y\mathbf 1\|^2$, $R^2=1-\mathrm{RSS}/\mathrm{TSS}$, $\hat\sigma^2=\mathrm{RSS}/(n-p)$, $h_{ii}=H_{ii}$. Thin SVD $X=U\Sigma V^\top$ with $\sigma_1\ge\dots\ge\sigma_p>0$ (LA-9). Ridge penalty $\lambda\|\mathbf b\|^2$, lasso penalty $\lambda\|\mathbf b\|_1$; both assume $\mathbf y$ centered and columns standardized, intercept unpenalized (§3.9).

## 0. Where this sits

ST-1 gave estimators and the Cramér–Rao bound; ST-2 Theorem 3.1 gave the $t$ and $F$ distributions by rotating Normal data so that one direction carries the mean. This module runs that argument with $p$ directions instead of one. LA-5's projection onto $C(X)$ is the estimator, ST-2's rotation is the distribution theory, and LA-9's SVD explains multicollinearity and ridge. You need LA-5 (normal equations, Pythagoras, QR), LA-9 (SVD, pseudoinverse, condition number), ST-2 ($t$, $F$, Wilks, ANOVA as projections), and F5 ($\operatorname{Cov}(A\mathbf X)=A\Sigma A^\top$, linear maps of Gaussians). It makes possible ST-4 (IRLS is weighted least squares in a loop), ST-5 (ridge is a Gaussian prior), ST-7 (regression with dependent errors), ST-9 (when a coefficient means anything), and ML-3 (regularized linear models).

## 1. The problem that forced this idea

Legendre published least squares in 1805 as a method for fitting comet orbits to more observations than unknowns; Gauss said in 1809 that he had used it since 1795, and justified it by assuming Normal errors (LA-5 tells that story). In 1821 Gauss gave a second justification that needs no Normality: among linear unbiased estimators, least squares has the smallest variance. That is the Gauss–Markov theorem of §3.4.

Regression got its name from data. Galton (1886) measured the heights of parents and adult children and found that children of tall parents were taller than average but *less* tall than their parents, which he called regression towards mediocrity. Pearson and others in the following decades built correlation and multiple regression on top. Fisher (1922–25) worked out the sampling distributions of the coefficients (with Student's $t$) and so made regression a tool for tests and intervals, not only curve fitting.

The forcing problem is the one you meet the first time you have a response and several noisy predictors: how does the response depend on each predictor, how sure are you, and which predictors matter? The geometry of LA-5 answers the first question and the probability of ST-1 and ST-2 answers the second. The third question broke the method. With many or nearly collinear predictors, least squares is unbiased but wildly variable. Hoerl and Kennard (1970) stabilized it by adding $\lambda I$ to $X^\top X$ (ridge), and Tibshirani (1996) changed the penalty from $\|\mathbf b\|^2$ to $\|\mathbf b\|_1$ and found that variable selection came free (the lasso).

## 2. Intuition

**Two pictures of one dataset.** The scatterplot has $n$ points in the plane of $(x,y)$ and the fitted line runs through them. The other picture, the one this module lives in, has *one* point: the vector $\mathbf y\in\mathbb R^n$. The columns $\mathbf 1$ and $\mathbf x$ are two more vectors in $\mathbb R^n$, and their span $C(X)$ is a plane through the origin (for $p=2$).

```
        y  (the data, a point in R^n)
        |\
        | \  e = y - yhat   (perpendicular to the plane)
        |  \
        |   \
        +----yhat------  C(X) = span{1, x}   (the plane of everything the model can say)
```

Least squares drops the perpendicular from $\mathbf y$ onto the plane (LA-5). Everything else is the geometry of that picture. Guess the main facts before they are proved:

1. The residual is perpendicular to every column of $X$, including $\mathbf 1$ (so residuals sum to zero) and every predictor.
2. Pythagoras holds: $\|\mathbf y-\bar y\mathbf 1\|^2=\|\hat{\mathbf y}-\bar y\mathbf 1\|^2+\|\mathbf e\|^2$, so $R^2$ is a squared cosine.
3. The shadow lives in $p$ dimensions and the residual in $n-p$, so noise in the residual should have $n-p$ degrees of freedom (ST-2 Theorem 3.1 with $p$ directions removed instead of one).
4. How well can the data pin down the coefficient on $\mathbf x_j$? Only the part of $\mathbf x_j$ that is *not* already in the span of the other columns carries information about $\beta_j$. If $\mathbf x_j$ lies almost inside that span, a tiny tilt of $\mathbf y$ swings the coefficient: huge variance. This is multicollinearity, and it is LA-9's condition number in disguise.
5. Shrinkage. The contours of $\mathrm{RSS}$ around $\hat{\boldsymbol\beta}$ are ellipses, long in the directions where the data are uninformative (small singular values). Ridge restrains the coefficients to a disc, so it pulls hardest along those long directions. The lasso restrains them to a diamond $\|\mathbf b\|_1\le t$, whose corners sit on the axes; a growing ellipse touches a diamond at a corner for a whole range of centers, and at a corner some coordinates are exactly zero.

## 3. Formal development

### 3.1 Simple regression

Model $y_i=\beta_0+\beta_1x_i+\varepsilon_i$. Write $\bar x,\bar y$ for the means, $S_{xx}=\sum(x_i-\bar x)^2$, $S_{xy}=\sum(x_i-\bar x)(y_i-\bar y)$, $S_{yy}=\sum(y_i-\bar y)^2$. Minimize $Q(b_0,b_1)=\sum_i(y_i-b_0-b_1x_i)^2$.

▶ **Try it first.** Set $\partial Q/\partial b_0=0$ and solve for $b_0$ in terms of $b_1$. Substitute into $\partial Q/\partial b_1=0$.

1. $\partial Q/\partial b_0=-2\sum(y_i-b_0-b_1x_i)=0$ gives $b_0=\bar y-b_1\bar x$. The line passes through $(\bar x,\bar y)$.
2. $\partial Q/\partial b_1=-2\sum x_i(y_i-b_0-b_1x_i)=0$. Substituting step 1: $\sum x_i\big[(y_i-\bar y)-b_1(x_i-\bar x)\big]=0$.
3. Since $\sum(y_i-\bar y)=\sum(x_i-\bar x)=0$, we may replace the leading $x_i$ by $x_i-\bar x$. This gives $S_{xy}=b_1S_{xx}$.

$$\hat\beta_1=\frac{S_{xy}}{S_{xx}},\qquad \hat\beta_0=\bar y-\hat\beta_1\bar x. \tag{3.1}$$

The Hessian of $Q$ is positive definite when $S_{xx}>0$ (not all $x_i$ equal), so this is the unique minimum.

■ **What just happened.** Centering decouples the intercept from the slope. The slope is the covariance of $x$ and $y$ divided by the variance of $x$, and the line is pinned to the center of mass of the data.

**Regression to the mean.** With $r=S_{xy}/\sqrt{S_{xx}S_{yy}}$ and $s_x,s_y$ the sample standard deviations, $\hat\beta_1=r\,s_y/s_x$, so
$$\frac{\hat y-\bar y}{s_y}=r\,\frac{x-\bar x}{s_x}. \tag{3.2}$$
Because $|r|<1$, the prediction is closer to the mean, measured in standard deviations, than the predictor was. This is the sample version of F5 Theorem 3.24: for a bivariate Normal, $\mathbb E[Y\mid X=x]=\mu_Y+\rho(\sigma_Y/\sigma_X)(x-\mu_X)$. Nothing pulls anything. A tall father is partly tall by luck, luck does not inherit, and so the children of the tallest fathers are tall *on average less so*. Regress $x$ on $y$ instead and you get the same effect in reverse: tall children have shorter fathers. Reading (3.2) as a causal force is the regression fallacy.

### 3.2 The model in matrix form: least squares is projection

We assume (A1) $X$ is fixed with $n>p$ and independent columns; (A2) $\mathbb E\boldsymbol\varepsilon=\mathbf 0$; (A3) $\operatorname{Var}\boldsymbol\varepsilon=\sigma^2I$, i.e. errors uncorrelated with equal variance; and later (A4) $\boldsymbol\varepsilon\sim\mathcal N(\mathbf 0,\sigma^2I)$.
$$\mathbf y=X\boldsymbol\beta+\boldsymbol\varepsilon. \tag{3.3}$$
"Linear" means linear in $\boldsymbol\beta$. Columns of $X$ may be $x^2$, $\log x$, dummy variables (§3.6). As $\mathbf b$ ranges over $\mathbb R^p$, $X\mathbf b$ ranges over the column space $C(X)$, a $p$-dimensional subspace of $\mathbb R^n$ that does not contain $\mathbf y$.

**Theorem 3.1 (Normal equations; LA-5 Theorem 3.15).** $\hat{\boldsymbol\beta}$ minimizes $\|\mathbf y-X\mathbf b\|^2$ if and only if
$$X^\top X\hat{\boldsymbol\beta}=X^\top\mathbf y, \tag{3.4}$$
and under (A1) the minimizer is unique:
$$\hat{\boldsymbol\beta}=(X^\top X)^{-1}X^\top\mathbf y. \tag{3.5}$$

*Proof, by geometry.* The point of $C(X)$ closest to $\mathbf y$ is the orthogonal projection (LA-5 Theorem 3.14), so the residual $\mathbf e=\mathbf y-X\hat{\boldsymbol\beta}$ is orthogonal to every column: $X^\top\mathbf e=\mathbf 0$, which is (3.4).
*By calculus.* $f(\mathbf b)=\|\mathbf y-X\mathbf b\|^2$ has gradient $-2X^\top(\mathbf y-X\mathbf b)$ (LA-5 (3.10); LA-14). Setting it to zero is (3.4). The Hessian $2X^\top X$ is positive definite under (A1) (LA-5 Theorem 3.11), so the stationary point is the unique global minimum, and $X^\top X$ is invertible. ∎

■ **What just happened.** "Set the gradient to zero" and "make the error perpendicular" are one sentence: the gradient is $-2X^\top\mathbf e$.

**Theorem 3.2 (The hat matrix).** $H=X(X^\top X)^{-1}X^\top$ satisfies: (a) $H^\top=H$; (b) $H^2=H$; (c) $HX=X$ and $(I-H)X=0$; (d) $\operatorname{tr}H=p$ and $\operatorname{tr}(I-H)=n-p$; (e) $X^\top\mathbf e=\mathbf 0$, so $\mathbf e\perp C(X)$, and $\sum e_i=0$ because $\mathbf 1$ is a column; (f) $\hat{\mathbf y}^\top\mathbf e=0$.

*Proof.* (a) $(X^\top X)^{-1}$ is symmetric. (b) $H^2=X(X^\top X)^{-1}(X^\top X)(X^\top X)^{-1}X^\top=H$. (c) $HX=X(X^\top X)^{-1}X^\top X=X$. (d) By $\operatorname{tr}(AB)=\operatorname{tr}(BA)$, $\operatorname{tr}[X\cdot(X^\top X)^{-1}X^\top]=\operatorname{tr}[(X^\top X)^{-1}X^\top X]=\operatorname{tr}I_p=p$. (e) $X^\top(I-H)\mathbf y=(X^\top-X^\top X(X^\top X)^{-1}X^\top)\mathbf y=\mathbf 0$. (f) $\hat{\mathbf y}^\top\mathbf e=\mathbf y^\top H(I-H)\mathbf y=0$ since $H-H^2=0$. ∎

$H$ is the projection matrix of LA-5 (3.8): a symmetric idempotent matrix has eigenvalues 0 or 1, so $H$ has $p$ ones (the directions of $C(X)$) and $n-p$ zeros.

**Theorem 3.3 (ANOVA decomposition).** If $\mathbf 1$ is a column of $X$, then
$$\|\mathbf y-\bar y\mathbf 1\|^2=\|\hat{\mathbf y}-\bar y\mathbf 1\|^2+\|\mathbf e\|^2,\quad\text{i.e. } \mathrm{TSS}=\mathrm{ESS}+\mathrm{RSS}. \tag{3.6}$$

*Proof.* $\mathbf y-\bar y\mathbf 1=(\hat{\mathbf y}-\bar y\mathbf 1)+\mathbf e$. Since $\mathbf 1\in C(X)$ and $\hat{\mathbf y}\in C(X)$, the first term is in $C(X)$; $\mathbf e\perp C(X)$ by Theorem 3.2(e). Apply Pythagoras (LA-5 (3.1)). ∎

■ **What just happened.** Total variation splits into explained and unexplained because the two pieces are at right angles. Without an intercept $\mathbf 1\notin C(X)$ in general and the split fails; software then silently reports an uncentered $R^2$.

**Geometry of $R^2$.** The vectors $\mathbf y-\bar y\mathbf 1$ (hypotenuse), $\hat{\mathbf y}-\bar y\mathbf 1$ and $\mathbf e$ (legs) form a right triangle. If $\theta$ is the angle between the hypotenuse and the adjacent leg,
$$R^2=\frac{\mathrm{ESS}}{\mathrm{TSS}}=\cos^2\theta=\operatorname{corr}(\mathbf y,\hat{\mathbf y})^2. \tag{3.7}$$
The last equality holds because $(\mathbf y-\bar y\mathbf 1)^\top(\hat{\mathbf y}-\bar y\mathbf 1)=\|\hat{\mathbf y}-\bar y\mathbf 1\|^2$ (the leg $\mathbf e$ is orthogonal to it). In simple regression $R^2=r^2$.

**Adjusted $R^2$.** Adding a column enlarges $C(X)$, and a minimum over a larger set cannot be larger, so $R^2$ never decreases, even if the column is noise. If all slopes are zero and errors are Normal, ESS and RSS are independent multiples of $\chi^2_{p-1}$ and $\chi^2_{n-p}$ (Theorem 3.6 below), so $R^2\sim\mathrm{Beta}(\tfrac{p-1}2,\tfrac{n-p}2)$ with mean $(p-1)/(n-1)$. The fix is to compare unbiased variance estimates (Theorem 3.5(c)) instead of raw sums:
$$R^2_{\rm adj}=1-\frac{\mathrm{RSS}/(n-p)}{\mathrm{TSS}/(n-1)}. \tag{3.8}$$
With $n=30$, $p=4$ and pure noise, the project finds mean $R^2=0.1046$ (theory $3/29=0.1034$) and mean $R^2_{\rm adj}=0.0013$. $R^2_{\rm adj}$ can be negative.

### 3.3 The Frisch–Waugh–Lovell theorem: what "controlling for" means

Partition $X=(X_1\;X_2)$ and let $M_1=I-X_1(X_1^\top X_1)^{-1}X_1^\top$, the projection onto $C(X_1)^\perp$ (LA-5 Theorem 3.13). Write $\hat{\boldsymbol\beta}=(\hat{\boldsymbol\beta}_1,\hat{\boldsymbol\beta}_2)$.

**Theorem 3.4 (Frisch–Waugh–Lovell).** (a) $\hat{\boldsymbol\beta}_2=(X_2^\top M_1X_2)^{-1}X_2^\top M_1\mathbf y$. Equivalently, $\hat{\boldsymbol\beta}_2$ is the coefficient from regressing the residualized $\tilde{\mathbf y}=M_1\mathbf y$ on the residualized $\tilde X_2=M_1X_2$. (b) That regression has the same residual vector $\mathbf e$ as the full one.
$$\hat{\boldsymbol\beta}_2=(\tilde X_2^\top\tilde X_2)^{-1}\tilde X_2^\top\tilde{\mathbf y}=(\tilde X_2^\top\tilde X_2)^{-1}\tilde X_2^\top\mathbf y. \tag{3.9}$$

▶ **Try it first.** Write (3.4) as two block equations, solve the first for $\hat{\boldsymbol\beta}_1$, substitute into the second.

1. Block equations: $X_1^\top X_1\hat{\boldsymbol\beta}_1+X_1^\top X_2\hat{\boldsymbol\beta}_2=X_1^\top\mathbf y$ and $X_2^\top X_1\hat{\boldsymbol\beta}_1+X_2^\top X_2\hat{\boldsymbol\beta}_2=X_2^\top\mathbf y$.
2. From the first: $\hat{\boldsymbol\beta}_1=(X_1^\top X_1)^{-1}X_1^\top(\mathbf y-X_2\hat{\boldsymbol\beta}_2)$.
3. Substitute into the second and collect: $X_2^\top\big[I-X_1(X_1^\top X_1)^{-1}X_1^\top\big](\mathbf y-X_2\hat{\boldsymbol\beta}_2)=\mathbf 0$, i.e. $X_2^\top M_1X_2\hat{\boldsymbol\beta}_2=X_2^\top M_1\mathbf y$. That is (a).
4. $M_1$ is symmetric and idempotent, so $X_2^\top M_1X_2=\tilde X_2^\top\tilde X_2$ and $X_2^\top M_1\mathbf y=\tilde X_2^\top\tilde{\mathbf y}=\tilde X_2^\top\mathbf y$. That is (3.9).
5. For (b), multiply $\mathbf y=X_1\hat{\boldsymbol\beta}_1+X_2\hat{\boldsymbol\beta}_2+\mathbf e$ by $M_1$. Since $M_1X_1=0$ and $M_1\mathbf e=\mathbf e$ ($\mathbf e\perp C(X_1)$), we get $\tilde{\mathbf y}=\tilde X_2\hat{\boldsymbol\beta}_2+\mathbf e$. ∎

■ **What just happened.** The coefficient on a predictor is the simple-regression slope of the *part of $\mathbf y$ the other columns cannot explain* on the *part of the predictor the other columns cannot explain*. That is the whole content of "holding the others fixed": it is a statement about projections in the data you have, not about what would happen if you changed $x_j$. One practical trap: the residualized regression has the right coefficients and residuals but its software output has the wrong degrees of freedom ($n-p_2$ instead of $n-p$, where $p_2$ is the number of columns of $X_2$; in Example 4.2 that is $7$ instead of $4$), so its standard errors are slightly off. Example 4.2 computes both routes.

### 3.4 Gauss–Markov

**Theorem 3.5.** Under (A1)–(A3): (a) $\mathbb E\hat{\boldsymbol\beta}=\boldsymbol\beta$; (b) $\operatorname{Var}(\hat{\boldsymbol\beta})=\sigma^2(X^\top X)^{-1}$; (c) $\mathbb E[\mathrm{RSS}]=(n-p)\sigma^2$, so $\hat\sigma^2$ is unbiased; (d) **(Gauss–Markov)** if $\tilde{\boldsymbol\beta}=C\mathbf y$ is any other linear estimator with $\mathbb E\tilde{\boldsymbol\beta}=\boldsymbol\beta$ for all $\boldsymbol\beta$, then $\operatorname{Var}(\tilde{\boldsymbol\beta})-\operatorname{Var}(\hat{\boldsymbol\beta})$ is positive semidefinite. In particular $\operatorname{Var}(\mathbf c^\top\tilde{\boldsymbol\beta})\ge\operatorname{Var}(\mathbf c^\top\hat{\boldsymbol\beta})$ for every $\mathbf c$: OLS is the **B**est **L**inear **U**nbiased **E**stimator.

▶ **Try it first.** For (d), write $C=(X^\top X)^{-1}X^\top+D$ and ask what unbiasedness forces on $D$.

1. Substitute (3.3): $\hat{\boldsymbol\beta}=(X^\top X)^{-1}X^\top(X\boldsymbol\beta+\boldsymbol\varepsilon)$, so
$$\hat{\boldsymbol\beta}-\boldsymbol\beta=(X^\top X)^{-1}X^\top\boldsymbol\varepsilon. \tag{3.10}$$
2. (a): take expectations in (3.10) and use (A2).
3. (b): with $A=(X^\top X)^{-1}X^\top$, $\operatorname{Cov}(A\boldsymbol\varepsilon)=A(\sigma^2I)A^\top$ (F5 Theorem 3.19) $=\sigma^2(X^\top X)^{-1}X^\top X(X^\top X)^{-1}$, so
$$\operatorname{Var}(\hat{\boldsymbol\beta})=\sigma^2(X^\top X)^{-1}. \tag{3.11}$$
4. (c): $\mathbf e=(I-H)\mathbf y=(I-H)\boldsymbol\varepsilon$ because $(I-H)X=0$. As $I-H$ is symmetric idempotent, $\mathrm{RSS}=\boldsymbol\varepsilon^\top(I-H)\boldsymbol\varepsilon=\operatorname{tr}[(I-H)\boldsymbol\varepsilon\boldsymbol\varepsilon^\top]$, so $\mathbb E[\mathrm{RSS}]=\sigma^2\operatorname{tr}(I-H)=(n-p)\sigma^2$ by Theorem 3.2(d):
$$\mathbb E[\mathrm{RSS}]=(n-p)\sigma^2. \tag{3.12}$$
(For $p=1$, $X=\mathbf 1$, this is ST-1 Theorem 3.6.)
5. (d): $\mathbb E\tilde{\boldsymbol\beta}=CX\boldsymbol\beta=\boldsymbol\beta$ for all $\boldsymbol\beta$ forces $CX=I$. Then $DX=CX-AX=I-I=0$.
6. $\operatorname{Var}(\tilde{\boldsymbol\beta})=\sigma^2CC^\top=\sigma^2(AA^\top+AD^\top+DA^\top+DD^\top)$. The cross term is $AD^\top=(X^\top X)^{-1}(DX)^\top=0$. So $\operatorname{Var}(\tilde{\boldsymbol\beta})=\operatorname{Var}(\hat{\boldsymbol\beta})+\sigma^2DD^\top$, and $DD^\top$ is positive semidefinite. ∎

■ **What just happened.** Any competitor is $\hat{\boldsymbol\beta}$ plus a piece $D\mathbf y=D\boldsymbol\varepsilon$ that is uncorrelated with $\hat{\boldsymbol\beta}$ ($\operatorname{Cov}=\sigma^2AD^\top=0$), and variances of uncorrelated pieces add. Geometrically $D$ has its rows orthogonal to $C(X)$: the competitor adds pure noise from the residual space. Example 4.3 specializes to the sample mean ($X=\mathbf 1$), where the theorem says equal weights minimize the variance of a weighted average.

**What "best" does not mean.** The class is *linear unbiased*. Biased estimators can beat OLS in mean squared error (ridge, §3.9, always does for small $\lambda$). Nonlinear estimators beat it under heavy-tailed errors (the median for a location). Under Normal errors the information matrix for $\boldsymbol\beta$ is $X^\top X/\sigma^2$, so (3.11) attains the Cramér–Rao bound (ST-1 Theorem 3.15, multivariate form) and nothing unbiased beats OLS at all. When (A3) fails the theorem fails: OLS is no longer best (§3.6).

### 3.5 Normal errors: the distribution theory

**Theorem 3.6.** Under (A1)–(A4): (a) $\hat{\boldsymbol\beta}\sim\mathcal N(\boldsymbol\beta,\sigma^2(X^\top X)^{-1})$; (b) $\mathrm{RSS}/\sigma^2\sim\chi^2_{n-p}$; (c) $\hat{\boldsymbol\beta}$ and $\mathrm{RSS}$ are independent.
$$\hat{\boldsymbol\beta}\sim\mathcal N\big(\boldsymbol\beta,\ \sigma^2(X^\top X)^{-1}\big). \tag{3.13}$$

▶ **Try it first.** Redo ST-2 Theorem 3.1 with an orthonormal basis of $\mathbb R^n$ whose first $p$ vectors span $C(X)$. What are $H\boldsymbol\varepsilon$ and $(I-H)\boldsymbol\varepsilon$ in that basis?

1. (a) is F5 Theorem 3.22 applied to (3.10).
2. Take $X=Q_1R$ (QR, LA-5 Theorem 3.20), so $\mathbf q_1,\dots,\mathbf q_p$ span $C(X)$. Extend to an orthonormal basis $\mathbf q_1,\dots,\mathbf q_n$ of $\mathbb R^n$ and let $W_i=\mathbf q_i^\top\boldsymbol\varepsilon/\sigma$. By F5 Theorems 3.22–3.23, $W_1,\dots,W_n$ are iid $\mathcal N(0,1)$.
3. $H\boldsymbol\varepsilon=\sigma\sum_{i\le p}W_i\mathbf q_i$ and $(I-H)\boldsymbol\varepsilon=\sigma\sum_{i>p}W_i\mathbf q_i$.
4. $\mathrm{RSS}=\|(I-H)\boldsymbol\varepsilon\|^2=\sigma^2\sum_{i>p}W_i^2\sim\sigma^2\chi^2_{n-p}$. That is (b).
5. $\hat{\boldsymbol\beta}-\boldsymbol\beta=R^{-1}Q_1^\top\boldsymbol\varepsilon$ depends on $W_1,\dots,W_p$ only; RSS depends on $W_{p+1},\dots,W_n$ only. Independent. ∎

■ **What just happened.** Degrees of freedom are dimensions: the fit uses $p$ directions and the residual keeps $n-p$. Everything below is an application.

**Corollary 3.7 ($t$ statistics, intervals).** With $c_{jj}=[(X^\top X)^{-1}]_{jj}$ and $\mathrm{se}_j=\hat\sigma\sqrt{c_{jj}}$,
$$\frac{\hat\beta_j-\beta_j}{\mathrm{se}_j}\sim t_{n-p},\qquad \hat\beta_j\pm t_{n-p,\,1-\alpha/2}\,\mathrm{se}_j\ \text{ is a }(1-\alpha)\text{ interval}. \tag{3.14, 3.15}$$
*Proof.* The numerator over $\sigma\sqrt{c_{jj}}$ is $\mathcal N(0,1)$; $\hat\sigma^2/\sigma^2=\chi^2_{n-p}/(n-p)$ is independent of it (Theorem 3.6); the ratio is a $t_{n-p}$ by ST-2 Definition 3.2. ∎ The same holds for any $\mathbf c^\top\boldsymbol\beta$ with $c_{jj}$ replaced by $\mathbf c^\top(X^\top X)^{-1}\mathbf c$.

**The mean response and a new observation.** At a point $\mathbf x_0$ (a $p$-vector with a leading 1) let $h_0=\mathbf x_0^\top(X^\top X)^{-1}\mathbf x_0$. The estimate $\mathbf x_0^\top\hat{\boldsymbol\beta}$ has variance $\sigma^2h_0$, so the interval for the *mean response* $\mathbf x_0^\top\boldsymbol\beta$ is
$$\mathbf x_0^\top\hat{\boldsymbol\beta}\pm t_{n-p,\,1-\alpha/2}\,\hat\sigma\sqrt{h_0}. \tag{3.16}$$
A *new observation* is $y_0=\mathbf x_0^\top\boldsymbol\beta+\varepsilon_0$ with $\varepsilon_0$ independent of the training errors. The prediction error $y_0-\mathbf x_0^\top\hat{\boldsymbol\beta}=\varepsilon_0-\mathbf x_0^\top(\hat{\boldsymbol\beta}-\boldsymbol\beta)$ is a sum of two independent Normal terms with mean zero, so its variance is $\sigma^2(1+h_0)$ and the same $t$ argument gives the **prediction interval**
$$\mathbf x_0^\top\hat{\boldsymbol\beta}\pm t_{n-p,\,1-\alpha/2}\,\hat\sigma\sqrt{1+h_0}. \tag{3.17}$$
The extra "1" is the noise of the one new observation, and no amount of data removes it: as $n\to\infty$, $h_0\to0$, the first interval shrinks to a point and the second to $\pm1.96\sigma$. The two are confused constantly (§8). Example 4.1 gives numbers.

**The $F$ test.** Let $X_0$ ($p_0$ columns) span a subspace $C(X_0)\subset C(X)$ and $q=p-p_0$. Write $H_0$ for its projection. Since $\mathbf y-H_0\mathbf y=(\mathbf y-H\mathbf y)+(H-H_0)\mathbf y$ is a sum of orthogonal pieces, $\mathrm{RSS}_0-\mathrm{RSS}=\|(H-H_0)\mathbf y\|^2$, the "extra sum of squares". $H-H_0$ projects onto $C(X)\cap C(X_0)^\perp$, of dimension $q$. If the true mean lies in $C(X_0)$ (the null), then $(H-H_0)\mathbf y=(H-H_0)\boldsymbol\varepsilon$, and the rotation of Theorem 3.6 in a basis adapted to $C(X_0)\subset C(X)\subset\mathbb R^n$ (ST-2 §3.7) gives $(\mathrm{RSS}_0-\mathrm{RSS})/\sigma^2\sim\chi^2_q$, independent of $\mathrm{RSS}/\sigma^2\sim\chi^2_{n-p}$. Hence
$$F=\frac{(\mathrm{RSS}_0-\mathrm{RSS})/q}{\mathrm{RSS}/(n-p)}\sim F_{q,\,n-p}\ \text{ under the null.} \tag{3.18}$$
For $q=1$, $F=t^2$. For the overall test (null: all slopes zero, $q=p-1$), $F=\dfrac{R^2/(p-1)}{(1-R^2)/(n-p)}$.

**Maximum likelihood.** The Normal log-likelihood is $\ell(\boldsymbol\beta,\sigma^2)=-\tfrac n2\log(2\pi\sigma^2)-\|\mathbf y-X\boldsymbol\beta\|^2/(2\sigma^2)$.
1. For any $\sigma^2$, maximizing over $\boldsymbol\beta$ minimizes $\|\mathbf y-X\boldsymbol\beta\|^2$: **OLS is the MLE**.
2. Setting $\partial\ell/\partial\sigma^2=-n/(2\sigma^2)+\mathrm{RSS}/(2\sigma^4)=0$ gives
$$\hat\sigma^2_{\rm MLE}=\mathrm{RSS}/n,\qquad \mathbb E\hat\sigma^2_{\rm MLE}=\tfrac{n-p}{n}\sigma^2. \tag{3.19}$$
The MLE is biased low, by the factor $(n-p)/n$; ST-1 Example 4.1 is the case $p=1$. The unbiased $\hat\sigma^2=\mathrm{RSS}/(n-p)$ is what makes (3.14) exactly $t$.
3. Substituting back, $\max\ell=-\tfrac n2\log(2\pi e\,\mathrm{RSS}/n)$. For nested models the LRT statistic is therefore
$$-2\log\lambda=n\log\frac{\mathrm{RSS}_0}{\mathrm{RSS}}=n\log\Big(1+\frac{qF}{n-p}\Big). \tag{3.20}$$
It is increasing in $F$, so the $F$ test **is** the likelihood-ratio test; for large $n$, $n\log(1+qF/(n-p))\approx qF\to\chi^2_q$, which is Wilks (ST-2 Theorem 3.9). If the errors are not Normal, OLS is a quasi-MLE: by the CLT the $t$ and $F$ results hold approximately for large $n$, but prediction intervals lean on the Normal tails of $\varepsilon_0$ and do not.

### 3.6 Diagnostics and remedies

**Leverage.** Since $\hat y_i=\sum_jH_{ij}y_j$, the diagonal element
$$h_{ii}=\mathbf x_i^\top(X^\top X)^{-1}\mathbf x_i=\frac{\partial\hat y_i}{\partial y_i} \tag{3.21}$$
measures how much $y_i$ pulls its own fit. It depends on $X$ only. By Theorem 3.2, $H_{ii}=\sum_jH_{ij}^2\ge H_{ii}^2$, so $0\le h_{ii}\le1$; $\sum_ih_{ii}=p$; with an intercept $h_{ii}\ge1/n$. A common flag is $h_{ii}>2p/n$. In simple regression $C(X)=\operatorname{span}\{\mathbf 1\}\oplus\operatorname{span}\{\mathbf x-\bar x\mathbf 1\}$ orthogonally, so $H=\mathbf 1\mathbf 1^\top/n+\tilde{\mathbf x}\tilde{\mathbf x}^\top/S_{xx}$ and $h_{ii}=1/n+(x_i-\bar x)^2/S_{xx}$: leverage is distance from the center of the $x$-cloud.

**Residuals are not errors.** $\mathbf e=(I-H)\boldsymbol\varepsilon$ has $\operatorname{Cov}(\mathbf e)=\sigma^2(I-H)(I-H)^\top=\sigma^2(I-H)$, so
$$\operatorname{Var}(e_i)=\sigma^2(1-h_{ii}),\qquad\operatorname{Cov}(e_i,e_j)=-\sigma^2H_{ij}. \tag{3.22}$$
Even under (A3) the raw residuals have unequal variances, and high-leverage points have *small* residuals because they drag the line toward themselves. The **standardized (internally studentized) residual** corrects the scale:
$$r_i=\frac{e_i}{\hat\sigma\sqrt{1-h_{ii}}}. \tag{3.23}$$

**Theorem 3.8 (Leave-one-out identity).** Let $\hat{\boldsymbol\beta}_{(i)}$ be the OLS estimate computed without observation $i$. Then
$$\hat{\boldsymbol\beta}-\hat{\boldsymbol\beta}_{(i)}=(X^\top X)^{-1}\mathbf x_i\frac{e_i}{1-h_{ii}},\tag{3.24}$$
$$y_i-\mathbf x_i^\top\hat{\boldsymbol\beta}_{(i)}=\frac{e_i}{1-h_{ii}},\tag{3.25}$$
$$\mathrm{RSS}_{(i)}=\mathrm{RSS}-\frac{e_i^2}{1-h_{ii}}.\tag{3.26}$$

▶ **Try it first.** $X_{(i)}^\top X_{(i)}=X^\top X-\mathbf x_i\mathbf x_i^\top$ is a rank-one downdate. Recall Sherman–Morrison (LA-2 Exercise C4) and apply it to $A=X^\top X$.

1. Let $A=X^\top X$ and $s=h_{ii}=\mathbf x_i^\top A^{-1}\mathbf x_i$. Sherman–Morrison: $(A-\mathbf x_i\mathbf x_i^\top)^{-1}=A^{-1}+\dfrac{A^{-1}\mathbf x_i\mathbf x_i^\top A^{-1}}{1-s}$. (Check: multiply by $A-\mathbf x_i\mathbf x_i^\top$; the cross terms are $\mathbf x_i\mathbf x_i^\top A^{-1}[-1+\tfrac1{1-s}-\tfrac s{1-s}]=0$.)
2. $X_{(i)}^\top\mathbf y_{(i)}=X^\top\mathbf y-\mathbf x_iy_i$. Multiply: $\hat{\boldsymbol\beta}_{(i)}=\hat{\boldsymbol\beta}-A^{-1}\mathbf x_iy_i+\dfrac{A^{-1}\mathbf x_i\hat y_i}{1-s}-\dfrac{sA^{-1}\mathbf x_iy_i}{1-s}$, using $A^{-1}X^\top\mathbf y=\hat{\boldsymbol\beta}$ and $\mathbf x_i^\top\hat{\boldsymbol\beta}=\hat y_i$.
3. The $y_i$ terms combine: $-y_i(1+\tfrac s{1-s})=-\tfrac{y_i}{1-s}$. So $\hat{\boldsymbol\beta}_{(i)}=\hat{\boldsymbol\beta}-A^{-1}\mathbf x_i\,\dfrac{y_i-\hat y_i}{1-s}$, which is (3.24).
4. Then $y_i-\mathbf x_i^\top\hat{\boldsymbol\beta}_{(i)}=e_i+\mathbf x_i^\top A^{-1}\mathbf x_i\dfrac{e_i}{1-s}=e_i\big(1+\tfrac s{1-s}\big)=\dfrac{e_i}{1-s}$, which is (3.25).
5. For (3.26), by (3.24) the residual of observation $j\ne i$ under $\hat{\boldsymbol\beta}_{(i)}$ is $e_j+H_{ji}e_i/(1-s)$. Square and sum over $j\ne i$, using $\sum_{j\ne i}e_jH_{ji}=-e_is$ (because $H\mathbf e=\mathbf 0$) and $\sum_{j\ne i}H_{ji}^2=s-s^2$:
$\mathrm{RSS}_{(i)}=(\mathrm{RSS}-e_i^2)-\dfrac{2se_i^2}{1-s}+\dfrac{(s-s^2)e_i^2}{(1-s)^2}=\mathrm{RSS}-e_i^2\Big(1+\dfrac s{1-s}\Big)$. ∎

■ **What just happened.** Deleting an observation is a rank-one change, and a rank-one change has a closed-form inverse. One fit gives all $n$ leave-one-out answers, with no refitting. Three uses follow.

*(i) Leave-one-out cross-validation for free:* $\mathrm{LOOCV}=\frac1n\sum_i\big(e_i/(1-h_{ii})\big)^2$ (the PRESS statistic), from (3.25). The factor $1/(1-h_{ii})$ inflates the residuals of high-leverage points, which the fit has bent toward.

*(ii) Externally studentized residuals:* $t_i=e_i/(\hat\sigma_{(i)}\sqrt{1-h_{ii}})$ with $\hat\sigma^2_{(i)}=\mathrm{RSS}_{(i)}/(n-p-1)$ from (3.26). Under (A1)–(A4), $t_i\sim t_{n-p-1}$ (exercise C3). They are the right statistic for "is observation $i$ an outlier?".

*(iii) Cook's distance.* Measure the shift in the fitted values, in units of the estimator's own uncertainty:
$$D_i=\frac{\|X(\hat{\boldsymbol\beta}-\hat{\boldsymbol\beta}_{(i)})\|^2}{p\hat\sigma^2}=\frac{r_i^2}{p}\cdot\frac{h_{ii}}{1-h_{ii}}.\tag{3.27}$$
Derivation: by (3.24), $X(\hat{\boldsymbol\beta}-\hat{\boldsymbol\beta}_{(i)})=XA^{-1}\mathbf x_i\,e_i/(1-h_{ii})$, whose squared norm is $\mathbf x_i^\top A^{-1}X^\top XA^{-1}\mathbf x_i\,e_i^2/(1-h_{ii})^2=h_{ii}e_i^2/(1-h_{ii})^2$. Divide by $p\hat\sigma^2$ and use (3.23). ∎ Influence is the product of *outlyingness in $y$* ($r_i^2$) and *outlyingness in $x$* ($h_{ii}/(1-h_{ii})$). Flags: $D_i>1$, or $>4/n$. **Leverage is about $x$, outlyingness about $y$, influence about both.** Example 4.4 shows the three coming apart.

**Residual plots: what each pattern means.**

| Plot | Pattern | Meaning and response |
|---|---|---|
| $e_i$ or $t_i$ vs $\hat y_i$ and vs each $x_j$ | curve | $E[y\mid x]$ is not linear in that column: add a term or transform |
| same | fan (spread grows with fitted value) | heteroscedasticity: sandwich or WLS (below) |
| $e_i$ vs time or order | runs of the same sign | autocorrelation: Durbin–Watson, ST-7 |
| QQ plot of $t_i$ against $t_{n-p-1}$ quantiles | S-shape, outlying ends | heavy tails: matters for prediction intervals and small-$n$ tests |
| $r_i^2$ vs $h_{ii}$, and $D_i$ | points far up or right | outlier, leverage, influence (above) |

**Heteroscedasticity.** Suppose $\operatorname{Var}(\boldsymbol\varepsilon)=\Sigma=\operatorname{diag}(\sigma_i^2)$. Step 2 of Theorem 3.5 did not use (A3), so $\hat{\boldsymbol\beta}$ is still unbiased. Its variance, however, is
$$\operatorname{Var}(\hat{\boldsymbol\beta})=(X^\top X)^{-1}X^\top\Sigma X(X^\top X)^{-1},\qquad\operatorname{Var}(\hat\beta_j)=\sum_iw_{ji}^2\sigma_i^2, \tag{3.28}$$
where $w_{ji}$ is entry $i$ of row $j$ of $(X^\top X)^{-1}X^\top$, so that (3.10) reads $\hat\beta_j-\beta_j=\sum_iw_{ji}\varepsilon_i$ (the matrix form is F5 Theorem 3.19 applied to (3.10); the scalar form is the variance of a weighted sum of uncorrelated terms). The textbook $\sigma^2(X^\top X)^{-1}$ is the wrong number, so $t$ intervals mis-cover. The **sandwich (Huber–White) estimator** replaces each unknown $\sigma_i^2$ by $e_i^2$:
$$\widehat{\operatorname{Var}}_{\rm HC0}(\hat\beta_j)=\sum_iw_{ji}^2e_i^2.\tag{3.29}$$
It is consistent because $e_i\approx\varepsilon_i$; HC3 uses the leave-one-out residuals $(e_i/(1-h_{ii}))^2$ from (3.25), the finite-sample repair. In the project, the error standard deviation grows like $x^2$ ($n=60$): the naive 95% slope interval covers 0.894 of the time, HC0 covers 0.934, HC3 0.943. The naive formula averages a variance of 0.092 where the truth is 0.134.

**Weighted least squares.** If $\operatorname{Var}(\varepsilon_i)=\sigma^2/w_i$ with known weights, multiply row $i$ by $\sqrt{w_i}$: with $W=\operatorname{diag}(w_i)$, the model $W^{1/2}\mathbf y=W^{1/2}X\boldsymbol\beta+W^{1/2}\boldsymbol\varepsilon$ has $\operatorname{Var}(W^{1/2}\boldsymbol\varepsilon)=\sigma^2I$, so (A3) holds and Theorem 3.5 applies to it:
$$\hat{\boldsymbol\beta}_W=(X^\top WX)^{-1}X^\top W\mathbf y,\qquad\operatorname{Var}(\hat{\boldsymbol\beta}_W)=\sigma^2(X^\top WX)^{-1}.\tag{3.30}$$
It is the BLUE for this $\Sigma$, and geometrically it is orthogonal projection in the weighted inner product $\langle\mathbf u,\mathbf v\rangle_W=\mathbf u^\top W\mathbf v$ (LA-1 Exercise C1). With the true weights in the project, the slope variance falls from 0.134 (OLS) to 0.0117. If the weights must be estimated, you are on your own: the gain is real but so is the risk. ST-4's IRLS is this formula in a loop with weights updated each pass.

**Autocorrelated errors, transformations.** If observations are ordered in time, (A3)'s "uncorrelated" usually fails. OLS stays unbiased but its standard errors are wrong (typically too small). The Durbin–Watson statistic $\sum(e_t-e_{t-1})^2/\sum e_t^2\approx2(1-\hat\rho_1)$ detects it; ST-7 handles it. For multiplicative noise, regress $\log y$ (a slope then reads as a proportional change); the Box–Cox family $(y^\lambda-1)/\lambda$, with $\lambda\to0$ giving $\log y$, chooses a transformation by profile likelihood (named, not developed).

**Categorical predictors, polynomials.** A factor with $k$ levels enters as $k-1$ dummy columns plus the intercept. The coefficient on a dummy is that level's mean minus the *reference* level's mean (other predictors held as in FWL). Using $k$ dummies *and* an intercept makes the columns dependent: the dummies sum to $\mathbf 1$, so $(1,-1,\dots,-1)^\top$ is in the null space of $X$ (LA-3), $X^\top X$ is singular (LA-5 Theorem 3.11), and (3.5) does not exist. Software either drops a column or returns the minimum-norm solution (LA-9 Theorem 3.10); the fitted values are the same, the coefficients are not identified (exercise A6). An interaction column $x\cdot d$ lets the slope differ by group, and its coefficient is the difference between the within-group slopes (exercise B8). Polynomial and spline terms ($x^2$, $(x-k)_+$) are more columns: the design matrix is what matters. Raw powers $1,x,x^2,\dots$ are nearly collinear, which is the next section, so the project uses Legendre polynomials.

### 3.7 Multicollinearity

**Theorem 3.9 (Variance inflation factor).** Let $\tilde{\mathbf x}_j$ be the residual of regressing $\mathbf x_j$ on all other columns of $X$ (including $\mathbf 1$), and $R_j^2$ the $R^2$ of that regression. With $S_{jj}=\sum_i(x_{ij}-\bar x_j)^2$,
$$\operatorname{Var}(\hat\beta_j)=\frac{\sigma^2}{\|\tilde{\mathbf x}_j\|^2}=\frac{\sigma^2}{S_{jj}}\cdot\mathrm{VIF}_j,\qquad\mathrm{VIF}_j=\frac1{1-R_j^2}.\tag{3.31}$$

*Proof.* By (3.9) with $X_2=\mathbf x_j$, $\hat\beta_j=\tilde{\mathbf x}_j^\top\mathbf y/\|\tilde{\mathbf x}_j\|^2$. Its variance is $\tilde{\mathbf x}_j^\top(\sigma^2I)\tilde{\mathbf x}_j/\|\tilde{\mathbf x}_j\|^4=\sigma^2/\|\tilde{\mathbf x}_j\|^2$. By Theorem 3.3 applied to the regression of $\mathbf x_j$ on the others, $\|\tilde{\mathbf x}_j\|^2=\mathrm{RSS}_j=(1-R_j^2)S_{jj}$. ∎

With $\mathbf x_j$ alone the variance would be $\sigma^2/S_{jj}$; $\mathrm{VIF}_j$ is the price of its correlation with the other predictors. This is guess 4 of §2, made exact.

**Theorem 3.10 (Variance in the SVD basis).** With $X=U\Sigma V^\top$ (thin; LA-9), $X^\top X=V\Sigma^2V^\top$ and $(X^\top X)^{-1}=V\Sigma^{-2}V^\top$, so
$$\operatorname{Var}(\hat{\boldsymbol\beta})=\sigma^2\sum_{i=1}^p\frac{\mathbf v_i\mathbf v_i^\top}{\sigma_i^2}.\tag{3.32}$$
Read it off: (i) $\operatorname{Var}(\mathbf v_i^\top\hat{\boldsymbol\beta})=\sigma^2/\sigma_i^2$ and $\operatorname{Cov}(\mathbf v_i^\top\hat{\boldsymbol\beta},\mathbf v_k^\top\hat{\boldsymbol\beta})=0$ for $i\ne k$: in the $V$ basis the estimates are uncorrelated, with variance inversely proportional to the squared singular value. (ii) $\mathbb E\|\hat{\boldsymbol\beta}-\boldsymbol\beta\|^2=\sigma^2\sum_i\sigma_i^{-2}$. (iii) The confidence ellipsoid $\{\mathbf b:\|X(\mathbf b-\hat{\boldsymbol\beta})\|^2\le c\}$ has semi-axes $\sqrt c/\sigma_i$ along $\mathbf v_i$: long and thin, longest along $\mathbf v_p$. (iv) The reason: $X\mathbf v_p=\sigma_p\mathbf u_p$ is tiny, so moving $\boldsymbol\beta$ one unit along $\mathbf v_p$ changes the fit by only $\sigma_p$. The data cannot distinguish $\boldsymbol\beta$ from $\boldsymbol\beta+t\mathbf v_p$ for moderate $t$; noise masquerades as signal in that direction. "Nearly dependent columns" means there is a unit vector $\mathbf v$ with $\|X\mathbf v\|$ small.

**Consequences.** (1) Multicollinearity does *not* bias $\hat{\boldsymbol\beta}$ ((3.10) holds regardless); it inflates variance. (2) Predictions at points resembling the training cloud stay fine ($\operatorname{Var}(\hat y_i)=\sigma^2h_{ii}$, summing to $p\sigma^2$); predictions that break the collinearity pattern are hopeless. (3) The condition number $\kappa(X)=\sigma_1/\sigma_p$ (LA-9 Theorem 3.9) is the diagnostic, and the ratio of the largest to the smallest direction variance in (3.32) is $\kappa^2$, the same $\kappa(X^\top X)=\kappa(X)^2$ that wrecks the normal equations (LA-13 Theorem 3.12). Rules of thumb: $\mathrm{VIF}>10$ ($R_j^2>0.9$), or $\kappa>30$ after scaling the columns to unit length (scale first, or units decide the answer). (4) Remedies: drop or combine predictors, collect data that breaks the pattern, or accept bias for variance (ridge, §3.9). Example 4.5 puts numbers on all of this.

### 3.8 Model selection

**Theorem 3.11 (Prediction error of a fitted subspace).** Let the true mean $\boldsymbol\mu$ be any vector (not necessarily in $C(X)$), $\mathbf y=\boldsymbol\mu+\boldsymbol\varepsilon$ under (A2)–(A3), and let $\mathbf y^*=\boldsymbol\mu+\boldsymbol\varepsilon^*$ be a new response at the same $X$ with $\boldsymbol\varepsilon^*$ independent of $\boldsymbol\varepsilon$. Then
$$\mathbb E\|\mathbf y^*-\hat{\mathbf y}\|^2=n\sigma^2+\|(I-H)\boldsymbol\mu\|^2+p\sigma^2,\qquad\mathbb E[\mathrm{RSS}]=\|(I-H)\boldsymbol\mu\|^2+(n-p)\sigma^2.\tag{3.33}$$

*Proof.* $\hat{\mathbf y}=H\boldsymbol\mu+H\boldsymbol\varepsilon$, so $\mathbf y^*-\hat{\mathbf y}=\boldsymbol\varepsilon^*+(I-H)\boldsymbol\mu-H\boldsymbol\varepsilon$: independent noise, a fixed vector, and a projected noise. Cross terms have zero expectation. $\mathbb E\|\boldsymbol\varepsilon^*\|^2=n\sigma^2$ and $\mathbb E\|H\boldsymbol\varepsilon\|^2=\sigma^2\operatorname{tr}H=p\sigma^2$ (the trace trick of (3.12)). For RSS: $\mathbf e=(I-H)(\boldsymbol\mu+\boldsymbol\varepsilon)$ and $\mathbb E\|(I-H)\boldsymbol\varepsilon\|^2=(n-p)\sigma^2$. ∎

■ **What just happened.** Prediction error = irreducible noise + bias$^2$ (shrinks as $C(X)$ grows) + variance $p\sigma^2$ (grows linearly in $p$). Training error is *optimistic*, by exactly
$$\mathbb E\|\mathbf y^*-\hat{\mathbf y}\|^2-\mathbb E[\mathrm{RSS}]=2p\sigma^2.\tag{3.34}$$
Adding the penalty back gives **Mallows' $C_p$**,
$$C_p=\mathrm{RSS}+2p\hat\sigma^2,\qquad\mathbb E[\mathrm{RSS}+2p\sigma^2]=\mathbb E\|\mathbf y^*-\hat{\mathbf y}\|^2,\tag{3.35}$$
with $\hat\sigma^2$ from the largest model. For a polynomial fit with a nonpolynomial truth ($n=30$, $\sigma=0.5$, $f(x)=\sin(5x+1)$), the project confirms (3.33) and (3.35) by Monte Carlo to within 3% at every degree; the expected prediction error is lowest at degree 6 (9.32); once the bias term has died out (degree 8 on) it rises by $\sigma^2=0.25$ per extra parameter.

**AIC and BIC.** For Gaussian regression, up to constants common to all models,
$$\mathrm{AIC}=n\log(\mathrm{RSS}/n)+2p,\qquad\mathrm{BIC}=n\log(\mathrm{RSS}/n)+p\log n.\tag{3.36}$$
*Sketch of AIC.* The goal is the expected out-of-sample log-likelihood of the fitted model, which (ST-1 Theorem 3.28) is minus the expected KL divergence from the truth up to a constant. The in-sample log-likelihood at the MLE is optimistic, and by the quadratic-approximation argument behind Wilks (ST-2 Theorem 3.9) the optimism is about $k$ for $k$ fitted parameters. So $-2\ell(\hat\theta)+2k$ estimates twice the expected out-of-sample negative log-likelihood. For Gaussian regression it behaves like $C_p/\hat\sigma^2$. BIC comes from a Laplace approximation to the marginal likelihood (ST-5) and penalizes by $\log n$, which exceeds 2 once $n\ge8$. They answer different questions: AIC targets prediction and keeps some extra terms; BIC targets the true finite model, if one exists, and is consistent. When $p/n$ is not small, AIC's penalty is too light: with $n=20$ and 17 candidate columns it picks an average degree of 13.8 out of a maximum of 16 (the corrected AICc, $\mathrm{AIC}+2p(p+1)/(n-p-1)$, is named for this). As $n$ grows to 1280 the mean selected degrees are 7.97 (AIC), 6.65 (BIC) and 7.99 (LOO): more data supports a more complex model.

**Cross-validation.** (3.25) gives leave-one-out CV with one fit, and not only for OLS: the proof of Theorem 3.8 goes through verbatim for any fit of the form $\hat{\boldsymbol\beta}=A^{-1}X^\top\mathbf y$ with $A=X^\top X+\Lambda$, $\Lambda$ a fixed symmetric matrix that makes $A$ positive definite, because deleting row $i$ then still subtracts exactly $\mathbf x_i\mathbf x_i^\top$ from $A$ (the step the proof needs). Ridge (§3.9) is $\Lambda=\lambda I$. $K$-fold needs refits; ML-1 develops it properly, including the optimism of the *winning* model's CV score. The project asserts the shortcut against brute-force refitting.

**Theorem 3.12 (Omitted-variable bias).** Let the truth be $\mathbf y=X_1\boldsymbol\beta_1+X_2\boldsymbol\beta_2+\boldsymbol\varepsilon$ under (A2)–(A3) and regress $\mathbf y$ on $X_1$ alone. Then
$$\mathbb E\hat{\boldsymbol\beta}_1^{\rm short}=\boldsymbol\beta_1+\Gamma\boldsymbol\beta_2,\qquad\Gamma=(X_1^\top X_1)^{-1}X_1^\top X_2.\tag{3.37}$$
*Proof.* $\hat{\boldsymbol\beta}_1^{\rm short}=(X_1^\top X_1)^{-1}X_1^\top(X_1\boldsymbol\beta_1+X_2\boldsymbol\beta_2+\boldsymbol\varepsilon)=\boldsymbol\beta_1+\Gamma\boldsymbol\beta_2+(X_1^\top X_1)^{-1}X_1^\top\boldsymbol\varepsilon$; take expectations. ∎ $\Gamma$ holds the coefficients from regressing each omitted column on the included ones. With one omitted variable the bias is $\gamma\beta_2$, $\gamma=$ slope of $\mathbf x_2$ on $\mathbf x_1$. It vanishes only if $\beta_2=0$ or $X_2\perp X_1$ (after centering), and it does not shrink as $n$ grows. There is an exact in-sample twin, $\hat{\boldsymbol\beta}_1^{\rm short}=\hat{\boldsymbol\beta}_1+\hat\Gamma\hat{\boldsymbol\beta}_2$ (exercise B3), which is FWL read backward.

**Stepwise selection and the post-selection swamp.** Forward, backward and stepwise procedures add or drop columns by $p$-values or AIC. Each step picks the best of many candidates, and the final table then reports $t$ statistics, $p$-values and an $F$ test *as if that model had been chosen in advance*. They are wrong: they condition on a selection event they ignore (the winner's curse of ST-2 §3.9–3.10). The project makes the point with pure noise: $n=100$, 50 candidate predictors, $y$ independent of all of them, forward selection by AIC. It keeps 10.8 predictors on average, reaches $R^2=0.34$, and 52% of the kept coefficients have $p<0.05$ (nominal: 0, 0, 5%). Remedies exist (sample splitting, selective inference) and this module only points at them: if you select, do not read the selected model's $p$-values.

**More predictors than observations ($p>n$).** Then $\operatorname{rank}X\le n<p$, $X^\top X$ is singular, and (3.4) has a $(p-\operatorname{rank}X)$-dimensional affine family of solutions; in the usual case $\operatorname{rank}X=n$ every one of them interpolates the data ($\mathrm{RSS}=0$). The **minimum-norm** solution is $\hat{\boldsymbol\beta}^+=X^+\mathbf y=\sum_{i\le r}(\mathbf u_i^\top\mathbf y/\sigma_i)\mathbf v_i$ (LA-9 Theorem 3.10), equal to $X^\top(XX^\top)^{-1}\mathbf y$ when $\operatorname{rank}X=n$. It is the $\lambda\to0^+$ limit of ridge (§3.9). It interpolates, so $\hat\sigma^2$ and every $t$ and $F$ test are undefined ($n-p<0$): this regime belongs to regularization.

### 3.9 Ridge regression

Penalize the size of the coefficients:
$$\hat{\boldsymbol\beta}_\lambda=\arg\min_{\mathbf b}\ \|\mathbf y-X\mathbf b\|^2+\lambda\|\mathbf b\|^2,\qquad\lambda\ge0.\tag{3.38}$$
The penalty depends on units (multiply $\mathbf x_j$ by 1000 and its coefficient shrinks by 1000, making it nearly free), so first **center $\mathbf y$ and every column and scale columns to unit variance**; the intercept is then $\bar y$, unpenalized.

**Theorem 3.13 (Ridge).** (a) For every $\lambda>0$, even if $X^\top X$ is singular,
$$\hat{\boldsymbol\beta}_\lambda=(X^\top X+\lambda I)^{-1}X^\top\mathbf y.\tag{3.39}$$
(b) In the SVD basis,
$$\hat{\boldsymbol\beta}_\lambda=\sum_{i=1}^p\frac{\sigma_i}{\sigma_i^2+\lambda}\,(\mathbf u_i^\top\mathbf y)\,\mathbf v_i=\sum_i\underbrace{\frac{\sigma_i^2}{\sigma_i^2+\lambda}}_{f_i}\ \frac{\mathbf u_i^\top\mathbf y}{\sigma_i}\,\mathbf v_i.\tag{3.40}$$
(c) $\hat{\mathbf y}_\lambda=\sum_if_i(\mathbf u_i^\top\mathbf y)\mathbf u_i$, and the **effective degrees of freedom** are
$$\mathrm{df}(\lambda)=\operatorname{tr}\big[X(X^\top X+\lambda I)^{-1}X^\top\big]=\sum_{i=1}^p\frac{\sigma_i^2}{\sigma_i^2+\lambda}.\tag{3.41}$$

▶ **Try it first.** Differentiate the objective. Then substitute $X=U\Sigma V^\top$ and see what is diagonal.

1. Gradient: $-2X^\top(\mathbf y-X\mathbf b)+2\lambda\mathbf b=\mathbf 0$, i.e. $(X^\top X+\lambda I)\mathbf b=X^\top\mathbf y$. The Hessian $2(X^\top X+\lambda I)$ has eigenvalues $\ge2\lambda>0$, so the minimizer exists, is unique and is (3.39).
2. $X^\top X+\lambda I=V(\Sigma^2+\lambda I)V^\top$ and $X^\top\mathbf y=V\Sigma U^\top\mathbf y$, so $\hat{\boldsymbol\beta}_\lambda=V(\Sigma^2+\lambda I)^{-1}\Sigma U^\top\mathbf y$. That is (3.40).
3. $X\hat{\boldsymbol\beta}_\lambda=U\Sigma(\Sigma^2+\lambda I)^{-1}\Sigma U^\top\mathbf y=U\operatorname{diag}(f_i)U^\top\mathbf y$, whose trace is $\sum f_i$. ∎

■ **What just happened.** Compare (3.40) with LA-9's pseudoinverse (3.11): OLS multiplies the $i$th coordinate $\mathbf u_i^\top\mathbf y$ by $1/\sigma_i$; ridge multiplies by $\sigma_i/(\sigma_i^2+\lambda)$, which is the same thing times the shrinkage factor $f_i\in(0,1)$. Directions with $\sigma_i^2\gg\lambda$ are untouched; directions with $\sigma_i^2\ll\lambda$ are almost deleted. A small singular value is an *uninformative* direction (§3.7), and ridge declines to chase it. Two more views: ridge is ordinary least squares on the augmented data $\binom{X}{\sqrt\lambda I},\binom{\mathbf y}{\mathbf 0}$ (the squared residual of the augmented problem is exactly (3.38)), which is also how to compute it stably by QR; and the dual form $\hat{\boldsymbol\beta}_\lambda=X^\top(XX^\top+\lambda I)^{-1}\mathbf y$ (the kernel-ridge form ML-4 uses) shows it works with $p>n$.

**Theorem 3.14 (Bias and variance of ridge).** Let $\boldsymbol\alpha=V^\top\boldsymbol\beta$ (the true coefficients in the $V$ basis). Then $\hat\alpha_i=\mathbf v_i^\top\hat{\boldsymbol\beta}_\lambda$ has bias $-\dfrac{\lambda\alpha_i}{\sigma_i^2+\lambda}$ and variance $\dfrac{\sigma^2\sigma_i^2}{(\sigma_i^2+\lambda)^2}$, and
$$\mathrm{MSE}(\lambda)=\mathbb E\|\hat{\boldsymbol\beta}_\lambda-\boldsymbol\beta\|^2=\sum_{i=1}^p\frac{\sigma^2\sigma_i^2+\lambda^2\alpha_i^2}{(\sigma_i^2+\lambda)^2}.\tag{3.42}$$
*Proof.* $\mathbf v_i^\top\hat{\boldsymbol\beta}_\lambda=\frac{\sigma_i}{\sigma_i^2+\lambda}\mathbf u_i^\top\mathbf y$ and $\mathbf u_i^\top\mathbf y=\sigma_i\alpha_i+\mathbf u_i^\top\boldsymbol\varepsilon$ (as $\mathbf u_i^\top X\boldsymbol\beta=\sigma_i\mathbf v_i^\top\boldsymbol\beta$), where $\mathbf u_i^\top\boldsymbol\varepsilon$ has variance $\sigma^2$. The mean is $f_i\alpha_i$, so the bias is $(f_i-1)\alpha_i$; the variance is $\frac{\sigma_i^2}{(\sigma_i^2+\lambda)^2}\sigma^2$. $V$ is orthogonal, so $\|\cdot\|^2$ is the sum over $i$ of bias$^2$ + variance. ∎

**Theorem 3.15 (Hoerl–Kennard, 1970).** For every $\boldsymbol\beta$ and $\sigma^2>0$ there is a $\lambda>0$ with $\mathrm{MSE}(\lambda)<\mathrm{MSE}(0)$.
*Proof.* In (3.42), the bias term $\lambda^2\alpha_i^2/(\sigma_i^2+\lambda)^2$ is $O(\lambda^2)$ and has zero slope at $\lambda=0$. The variance term has slope $-2\sigma^2/\sigma_i^4$ there. So $\mathrm{MSE}'(0)=-2\sigma^2\sum_i\sigma_i^{-4}<0$, and a small enough $\lambda$ lowers the MSE. ∎
For an orthonormal design (all $\sigma_i=1$) (3.42) is $(p\sigma^2+\lambda^2\|\boldsymbol\beta\|^2)/(1+\lambda)^2$, minimized at $\lambda^*=p\sigma^2/\|\boldsymbol\beta\|^2$ with $\mathrm{MSE}(\lambda^*)=p\sigma^2/(1+p\sigma^2/\|\boldsymbol\beta\|^2)<p\sigma^2$. The catch is that $\lambda^*$ depends on the unknown $\boldsymbol\beta$, and the gain $2\sigma^2/\sigma_i^4$ is enormous exactly when some $\sigma_i$ is small.

**Ridge as a Bayesian mode (ST-5 preview).** Let $\boldsymbol\beta\sim\mathcal N(\mathbf 0,\tau^2I)$ and $\mathbf y\mid\boldsymbol\beta\sim\mathcal N(X\boldsymbol\beta,\sigma^2I)$. The log posterior is $-\frac1{2\sigma^2}\big[\|\mathbf y-X\boldsymbol\beta\|^2+(\sigma^2/\tau^2)\|\boldsymbol\beta\|^2\big]+\text{const}$: ridge with $\lambda=\sigma^2/\tau^2$. The posterior is Gaussian with covariance $\sigma^2(X^\top X+\lambda I)^{-1}$, so mean and mode coincide.

**Choosing $\lambda$.** By cross-validation. The derivation of Theorem 3.8 used only that $\hat{\boldsymbol\beta}=A^{-1}X^\top\mathbf y$ and that deleting row $i$ replaces $A$ by $A-\mathbf x_i\mathbf x_i^\top$; both hold for $A=X^\top X+\lambda I$, because the penalty does not depend on which rows are present. With $h_{ii}$ the diagonal of $H_\lambda=X A^{-1}X^\top$, (3.25) holds verbatim. The project verifies it against brute force. Generalized CV replaces each $h_{ii}$ by $\mathrm{df}(\lambda)/n$.

### 3.10 The lasso

$$\hat{\boldsymbol\beta}^{\rm L}_\lambda=\arg\min_{\mathbf b}\ \tfrac12\|\mathbf y-X\mathbf b\|^2+\lambda\|\mathbf b\|_1.\tag{3.43}$$
(The $\tfrac12$ keeps the formulas clean; this $\lambda$ is not the ridge $\lambda$. Same standardization as §3.9.)

**Why exact zeros: geometry.** Equivalently minimize $\|\mathbf y-X\mathbf b\|^2$ subject to $\|\mathbf b\|_1\le t$. The contours are ellipses centered at $\hat{\boldsymbol\beta}$ and the answer is where the first contour touches the constraint set. The $\ell^2$ ball is smooth, so contact is generically at a point with every coordinate nonzero. The $\ell^1$ ball (LA-1 §3.6) is a cross-polytope with corners on the axes, where all but one coordinate vanish; a growing ellipse reaches a corner for a positive-measure set of centers. The analytic version follows.

**Theorem 3.16 (Soft-thresholding).** If $X^\top X=I$ and $\mathbf z=X^\top\mathbf y$ (which is $\hat{\boldsymbol\beta}_{\rm OLS}$), then
$$\hat\beta^{\rm L}_j=S_\lambda(z_j)=\operatorname{sign}(z_j)\,(|z_j|-\lambda)_+.\tag{3.44}$$

▶ **Try it first.** Show $\tfrac12\|\mathbf y-X\mathbf b\|^2=\tfrac12\|\mathbf b-\mathbf z\|^2+\text{const}$. Then minimize $\tfrac12(b-z)^2+\lambda|b|$ for one variable, by cases.

1. With $X^\top X=I$: $\tfrac12\|\mathbf y\|^2-\mathbf b^\top\mathbf z+\tfrac12\|\mathbf b\|^2=\tfrac12\|\mathbf b-\mathbf z\|^2+\big(\tfrac12\|\mathbf y\|^2-\tfrac12\|\mathbf z\|^2\big)$. The problem separates into $p$ scalar problems $\min_b\ g(b)=\tfrac12(b-z)^2+\lambda|b|$, each convex.
2. $|b|$ has a kink at 0. Its *subdifferential* is $\{\operatorname{sign}b\}$ for $b\ne0$ and the whole interval $[-1,1]$ at $b=0$. A convex $g$ is minimized iff $0\in\partial g(b)$, i.e. $0\in b-z+\lambda\,\partial|b|$.
3. Cases. If $b>0$: $b=z-\lambda$, valid when $z>\lambda$. If $b<0$: $b=z+\lambda$, valid when $z<-\lambda$. If $b=0$: need $0\in-z+\lambda[-1,1]$, i.e. $|z|\le\lambda$. These cover every $z$ exactly once. ∎

■ **What just happened.** OLS keeps $z$; ridge returns $z/(1+\lambda)$ (shrinks every coefficient, zeros none); the lasso shifts $z$ toward zero by $\lambda$ and clips at zero. The whole interval $[-\lambda,\lambda]$ maps to 0 because the penalty has slope $\lambda$ on each side of 0 and beats any gradient smaller than that. The exact zeros come from the kink in $|b|$ (the corner of the $\ell^1$ ball), not from magic. Example 4.7.

**Theorem 3.17 (Optimality conditions, general $X$).** $\hat{\mathbf b}$ minimizes (3.43) iff, writing $\mathbf r=\mathbf y-X\hat{\mathbf b}$,
$$\mathbf x_j^\top\mathbf r=\lambda\operatorname{sign}(\hat b_j)\ \text{ if }\hat b_j\ne0,\qquad|\mathbf x_j^\top\mathbf r|\le\lambda\ \text{ if }\hat b_j=0.\tag{3.45}$$
*Proof.* The objective is convex, so $\mathbf 0$ must lie in its subdifferential; coordinatewise that is $-\mathbf x_j^\top\mathbf r+\lambda s_j=0$ with $s_j=\operatorname{sign}\hat b_j$ or $s_j\in[-1,1]$ when $\hat b_j=0$. ∎ Active predictors have residual correlation of magnitude exactly $\lambda$; no inactive predictor exceeds it. Putting $\hat{\mathbf b}=\mathbf 0$ gives $\lambda_{\max}=\max_j|\mathbf x_j^\top\mathbf y|$, the smallest $\lambda$ with every coefficient zero.

**Coordinate descent.** Fix all coordinates but $b_j$ and let $\mathbf r_j=\mathbf y-\sum_{k\ne j}\mathbf x_kb_k$ be the partial residual. As a function of $b_j$ the objective is $\tfrac12\|\mathbf r_j-\mathbf x_jb_j\|^2+\lambda|b_j|+\text{const}=\tfrac{\|\mathbf x_j\|^2}2\big(b_j-\mathbf x_j^\top\mathbf r_j/\|\mathbf x_j\|^2\big)^2+\lambda|b_j|+\text{const}'$, a scalar problem of the form in Theorem 3.16, so
$$b_j\leftarrow\frac{S_\lambda(\mathbf x_j^\top\mathbf r_j)}{\|\mathbf x_j\|^2}.\tag{3.46}$$
Cycle through $j=1,\dots,p$, updating $\mathbf r$ in $O(n)$ per coordinate. Each update is an exact minimization along a coordinate, so the objective never increases; because the nonsmooth part is separable the iterates converge to a global minimum (Tseng, 2001; stated). Run it from $\lambda_{\max}$ downward, warm-starting each $\lambda$ from the last: the **lasso path** is piecewise linear in $\lambda$ (Efron et al., 2004, stated) and the project plots it. ML-2 meets coordinate descent again in general optimization.

**Ridge versus lasso.** On the project's benchmark ($p=50$, $n=100$ training points, error $\|\hat{\mathbf b}-\boldsymbol\beta\|$ averaged over 12 runs, $\lambda$ chosen on a validation set): sparse truth, ridge 1.68 vs lasso 0.81; dense truth (all 50 coefficients small), ridge 1.54 vs lasso 1.63; sparse truth with predictor correlation 0.9, ridge 3.40 vs lasso 1.86. So: sparse truth, lasso; dense truth, ridge, by a smaller margin. What correlation breaks is *selection stability*. From a group of correlated predictors the lasso keeps one, arbitrarily; with two identical columns the solution is not even unique (exercise C2), where ridge splits the weight evenly. The **elastic net** (Zou and Hastie, 2005), penalty $\lambda_1\|\mathbf b\|_1+\lambda_2\|\mathbf b\|^2$, repairs this (exercise C4). Three more facts. The $\lambda$ that minimizes CV prediction error over-selects: on the sparse problem it keeps 27 of 50 predictors, all 5 true ones and 22 false ones, although its error is 0.40 against 1.20 for OLS. The lasso also shrinks large coefficients (it subtracts $\lambda$ from all of them), so people refit OLS on the selected set, at the price of post-selection inference (§3.8). And recovering the exact support needs conditions on the design (the "irrepresentable condition", Zhao and Yu, 2006; named). In Bayesian terms the lasso is the posterior mode under independent Laplace priors $\propto e^{-|\beta_j|/b}$ (ST-5); the posterior mean is not sparse.

### 3.11 Regression is not causation; what the model assumes

A regression coefficient describes how the conditional mean of $y$ moves *across the observations you have*. It is causal only if the omitted-variable bias (3.37) is zero by design (randomization) or by an argument you must defend. If a confounder $Z$ affects both $x$ and $y$ and is omitted, the slope on $x$ estimates $\beta_x+\gamma\beta_z$ with $\gamma$ the slope of $z$ on $x$, which is Theorem 3.12 again.

**Simpson's paradox in regression form.** Take group 0: $x=(1,2,3)$, $y=(8,9,10)$; group 1: $x=(6,7,8)$, $y=(1,2,3)$. The pooled slope is $-1.169$ (intercept 10.76). Add a group dummy $g$: the fit becomes $y=7+x-12g$ exactly, slope $+1$ within each group. By FWL, adding $g$ replaces $x$ and $y$ by their within-group deviations, removing the between-group trend. (This is F4's Simpson trap, now with a coefficient.) Which slope is right depends on the causal story: is $g$ a confounder that should be held fixed, or a mediator of the effect of $x$ that should not? The data cannot say; ST-9's DAGs can.

| Assumption | What it buys | Can you check it? |
|---|---|---|
| $E[y\mid x]$ linear in the columns | unbiased $\hat{\boldsymbol\beta}$ | yes: residual plots vs fitted and each $x_j$; add terms |
| $E[\boldsymbol\varepsilon\mid X]=0$, no omitted confounders | causal reading, unbiasedness under misspecification | **no**, not from the data: design or argument |
| uncorrelated errors | valid standard errors | partly: residual order plot, Durbin–Watson; ST-7 |
| constant variance | valid standard errors and Gauss–Markov | yes: scale–location plot; or use the sandwich |
| Normal errors | exact $t$, $F$, prediction intervals | yes: QQ plot (matters little for CIs at large $n$) |
| $X$ measured without error, full rank | no attenuation, finite variance | rank/VIF yes; measurement error no |

## 4. Worked examples

### Example 4.1: Simple regression by hand, five points

Data $x=(1,2,3,4,5)$, $y=(2,4,5,4,5)$. Then $\bar x=3$, $\bar y=4$, $S_{xx}=4+1+0+1+4=10$, $S_{xy}=(-2)(-2)+(-1)(0)+0+(1)(0)+(2)(1)=6$, $S_{yy}=4+0+1+0+1=6$. By (3.1): $\hat\beta_1=0.6$, $\hat\beta_0=4-0.6\cdot3=2.2$. Fitted values $2.8,3.4,4.0,4.6,5.2$; residuals $\mathbf e=(-0.8,0.6,1.0,-0.6,-0.2)$, which sum to 0 and satisfy $\sum x_ie_i=-0.8+1.2+3.0-2.4-1.0=0$ ✓.

| Source | SS | df | MS | $F$ |
|---|---|---|---|---|
| Regression | $\mathrm{ESS}=3.6$ | 1 | 3.6 | $3.6/0.8=4.5$ |
| Residual | $\mathrm{RSS}=0.64+0.36+1+0.36+0.04=2.4$ | 3 | 0.8 | |
| Total | $\mathrm{TSS}=6$ | 4 | | |

$R^2=3.6/6=0.6$ (=$r^2$), $R^2_{\rm adj}=1-0.8/1.5=0.4667$, $\hat\sigma^2=0.8$. Standard errors: $\mathrm{se}(\hat\beta_1)=\sqrt{0.8/10}=0.2828$ and $\mathrm{se}(\hat\beta_0)=\sqrt{0.8(1/5+9/10)}=0.9381$. Then $t=0.6/0.2828=2.1213$ on 3 df, two-sided $p=0.124$; $F=4.5=t^2$ ✓. With $t_{3,0.975}=3.1824$ the 95% interval for the slope is $0.6\pm0.9001=(-0.300,\,1.500)$, which contains 0, consistent with $p>0.05$. Leverages $h_{ii}=1/5+(x_i-3)^2/10=(0.6,0.3,0.2,0.3,0.6)$, summing to $p=2$ ✓.

Intervals at two new $x_0$, with $w=t_{3,0.975}\hat\sigma=2.8465$:

| $x_0$ | $\hat y_0$ | $h_0$ | CI half-width (3.16) | PI half-width (3.17) |
|---|---|---|---|---|
| 3 | 4.0 | 0.2 | $w\sqrt{0.2}=1.273$ | $w\sqrt{1.2}=3.118$ |
| 5 | 5.2 | 0.6 | $w\sqrt{0.6}=2.205$ | $w\sqrt{1.6}=3.600$ |

The interval for the mean is narrowest at $\bar x$ and widens with leverage; the prediction interval is dominated by the "1" and barely moves.

### Example 4.2: Frisch–Waugh–Lovell, both routes

Eight observations; columns $\mathbf 1$, $\mathbf x_1=(1,\dots,8)$, $\mathbf x_2=(2,1,4,3,6,5,8,6)$, $\mathbf x_3=(1,0,2,1,3,1,2,4)$, $\mathbf y=(3.1,2.2,5.9,5.0,8.1,6.8,9.9,9.3)$.

*Route 1 (full fit, QR):* $\hat{\boldsymbol\beta}=(0.9948,\,0.2529,\,0.7944,\,0.3883)$, $\mathrm{RSS}=0.12579$.

*Route 2 (FWL for $\beta_3$):* regress $\mathbf x_3$ on $(\mathbf 1,\mathbf x_1,\mathbf x_2)$ and keep the residual $\tilde{\mathbf x}_3=(0.3727,-0.5206,0.5674,-0.3258,0.7622,-1.1311,-1.0431,1.3184)$. Then $\tilde{\mathbf x}_3^\top\mathbf y=2.1449$ and $\|\tilde{\mathbf x}_3\|^2=5.5243$, so by (3.9) the slope is $2.1449/5.5243=0.38827$, matching Route 1 to 5 digits. The residuals of $\tilde{\mathbf y}$ on $\tilde{\mathbf x}_3$ equal the full-fit residuals (RSS $0.12579$). Extras: $\operatorname{Var}(\hat\beta_3)=\sigma^2/\|\tilde{\mathbf x}_3\|^2$ (3.31), so $\mathrm{se}=\sqrt{(0.12579/4)/5.5243}=0.0755$ and $t=5.15$ on 4 df. The regression of $\mathbf y$ on $\mathbf x_3$ alone has slope $1.737$: "controlling for" $\mathbf x_1,\mathbf x_2$ moved the coefficient by a factor of 4.5, because $\mathbf x_3$ is correlated with the other predictors.

### Example 4.3: Gauss–Markov for the sample mean

Take $X=\mathbf 1$, $p=1$, $n=4$, so $\hat\beta=\bar y$ with $\operatorname{Var}=\sigma^2/4$. A linear unbiased estimator $\sum c_iy_i$ needs $\sum c_i=1$. Try $\mathbf c=(0.4,0.3,0.2,0.1)$. Write $\mathbf c=\tfrac14\mathbf 1+\mathbf d$ with $\mathbf d=(0.15,0.05,-0.05,-0.15)$, so $\sum d_i=0$ (that is $DX=0$). Then $\operatorname{Var}=\sigma^2\sum c_i^2=\sigma^2(0.25+\|\mathbf d\|^2)=\sigma^2(0.25+0.05)=0.30\sigma^2$, and $\operatorname{Cov}(\bar y,\mathbf d^\top\mathbf y)=\sigma^2\cdot\tfrac14\sum d_i=0$. Gauss–Markov says: among weighted averages, equal weights have the smallest variance, the ordinary result in the language of §3.4 step 6.

### Example 4.4: Leverage, outlyingness and influence

Four points with $x=(1,2,3,4)$, $y=(1.1,1.9,3.2,3.9)$ and a fifth at $x_5=10$, far from the others. Leverages are $(0.38,0.28,0.22,0.20,0.92)$, summing to 2: the fifth point has $h_{55}=0.92\gg2p/n=0.8$ already flagged. Two versions of $y_5$:

| | $y_5=9.9$ (on trend) | $y_5=3.0$ (off trend) |
|---|---|---|
| fitted slope | 0.982 | 0.154 |
| raw residual $e_5$ | 0.008 | $-0.544$ |
| LOO error $e_5/(1-h_{55})$ | 0.10 | $-6.80$ |
| standardized $r_5$ | 0.194 | $-1.717$ |
| Cook's $D_5$ | 0.216 (below its neighbours' maximum, 0.282) | 16.96 (others $\le0.441$) |

In the second case the raw residual $-0.544$ looks harmless, and $|r_5|=1.7$ is unremarkable, yet one point has flipped the slope from 0.98 to 0.15. The point bent the line toward itself, shrinking its own residual; only the leave-one-out error ($-6.80$, the gap between 3.0 and the 9.8 predicted from the other four) and Cook's distance (3.27) expose it. In the first case the same leverage is harmless. Leverage is about $x$, outlyingness about $y$, influence about both.

### Example 4.5: Two nearly collinear predictors, $\rho=0.99$

Take two columns of unit length with $X^\top X=\begin{pmatrix}1&0.99\\0.99&1\end{pmatrix}$. Eigenvalues $1\pm0.99$, so $\sigma_1^2=1.99$, $\sigma_2^2=0.01$, $\sigma_{1,2}=(1.4107,\ 0.1)$, $\kappa(X)=\sqrt{199}=14.1$, with $\mathbf v_1=(1,1)/\sqrt2$ and $\mathbf v_2=(1,-1)/\sqrt2$. By (3.32),
$$\frac{\operatorname{Var}(\hat\beta_1)}{\sigma^2}=\frac{1/2}{1.99}+\frac{1/2}{0.01}=0.2513+50.0=50.2513.$$
By (3.31) the same number is $\mathrm{VIF}=1/(1-0.99^2)=50.2513$ (here $S_{11}=1$): the standard error is inflated by $\sqrt{50.2513}=7.09$. The SVD explains where it comes from: 99.5% of it is the $\mathbf v_2$ term, the direction in which the columns are almost indistinguishable. In that basis $\operatorname{Var}(\hat\beta_1+\hat\beta_2)=2\sigma^2/(1+\rho)=1.005\sigma^2$ (well determined, the direction along which the columns agree) but $\operatorname{Var}(\hat\beta_1-\hat\beta_2)=2\sigma^2/(1-\rho)=200\sigma^2$, and $\operatorname{Corr}(\hat\beta_1,\hat\beta_2)=-0.99$. The data know the sum of the two coefficients and almost nothing about how it splits. Simulation (20,000 replications) gives variance within 3% of 50.25, correlation $-0.99$, and means equal to the truth: unbiased, nearly useless.

### Example 4.6: Ridge shrinkage factors, $\sigma_i=(10,1,0.1)$, $\lambda=1$

$f_i=\sigma_i^2/(\sigma_i^2+1)=100/101,\ 1/2,\ 0.01/1.01=(0.9901,\ 0.5,\ 0.0099)$. Effective df $=\sum f_i=1.5$ (exactly: $100/101+1/101=1$, plus $1/2$). The coefficient multipliers on $\mathbf u_i^\top\mathbf y$ are $\sigma_i/(\sigma_i^2+1)=(0.0990,\ 0.5,\ 0.0990)$ against OLS's $1/\sigma_i=(0.1,\ 1,\ 10)$: the weak direction's multiplier drops from 10 to 0.099. MSE (3.42) with $\sigma^2=1$ and $\boldsymbol\alpha=(1,1,1)$:
$$\mathrm{MSE}(0)=\tfrac1{100}+1+100=101.01,\qquad\mathrm{MSE}(1)=\tfrac{1+100}{101^2}+\tfrac{1+1}{2^2}+\tfrac{0.01+1}{1.01^2}=0.0099+0.5+0.9901=1.5000.$$
Ridge trades a bias of $-0.99\alpha_3$ in the weak direction for removing a variance of 100. Simulation (4000 replications) gives 1.509 for ridge.

### Example 4.7: Lasso on an orthonormal design

Let $X^\top X=I$ with $\mathbf z=X^\top\mathbf y=(3,\,-1.5,\,0.4,\,-0.2,\,2.2)$ and $\lambda=1$. By (3.44): $S_1(\mathbf z)=(2,\,-0.5,\,0,\,0,\,1.2)$ — the two smallest coefficients are zeroed and the others lose exactly 1 in magnitude. Ridge at $\lambda=1$ gives $\mathbf z/2=(1.5,\,-0.75,\,0.2,\,-0.1,\,1.1)$: all five nonzero. Check (3.45): the residual is $\mathbf r=X(\mathbf z-\hat{\mathbf b})+\mathbf r_\perp$, so $\mathbf x_j^\top\mathbf r=z_j-\hat b_j$, which is $\pm1=\lambda\operatorname{sign}\hat b_j$ for the three active predictors, and equals $0.4,\,-0.2$ for the two inactive ones, with $|{\cdot}|\le\lambda$ ✓. The project's coordinate descent returns exactly this vector, and agrees to better than $10^{-7}$ with an independent L-BFGS-B solve of the split problem on a harder $n=30$, $p=6$ design.

## 5. Why it matters and where it connects

- **ST-4 (GLMs):** a GLM is fitted by iteratively reweighted least squares: each step is the WLS (3.30) with weights from the current fit, and Wald $z$ statistics are (3.14) with $(X^\top WX)^{-1}$.
- **ST-5 (Bayesian inference):** ridge is the posterior under a Gaussian prior (§3.9), the lasso under a Laplace prior; Bayesian linear regression reuses (3.39) with $\lambda=\sigma^2/\tau^2$.
- **ST-7 (time series):** the partial autocorrelation at lag $h$ is FWL (3.9) applied to the intermediate lags; Durbin–Watson and autocorrelated-error regression build on §3.6.
- **ST-9 (causal inference):** (3.37) is the confounding formula, "controlling for" is FWL, and the question of *which* variables to include is what DAGs answer.
- **ML-1/ML-2/ML-3:** the bias–variance trade-off (3.33), the LOO shortcut (3.25), regularization paths and coordinate descent, and ridge as weight decay.
- **ML-4:** kernel ridge regression is the dual form $X^\top(XX^\top+\lambda I)^{-1}\mathbf y$; its smoother matrix $K(K+\lambda I)^{-1}$ has effective df of the form (3.41); a Gaussian-process posterior mean is the same object.
- **LA-9 / LA-13 (back-link):** the SVD is not a side remark here: it is what makes (3.32) and (3.40) readable. $\kappa(X)^2$ in the normal equations is a concrete cost, not a footnote (Part 1: at $\kappa=10^8$ the normal equations are wrong by 0.6, QR by $4\times10^{-9}$).

## 6. Exercises

See [exercises.md](exercises.md): 18 problems. A1–A6 are drills: simple regression, leverages and LOO errors, CI vs PI, soft-thresholding, ridge factors, the dummy-variable trap. B1–B8 are the module: hat-matrix facts (B1), $R^2$ and adjusted $R^2$ (B2), the omitted-variable identity (B3), Gauss–Markov by projection (B4), $F$ = LRT (B5), the sandwich (B6), collinearity and the SVD (B7), interaction and FWL in `project.py` (B8). B1, B3, B4, B6 and B8 are essential. C1–C4 go one step further: stepwise post-selection (code), duplicated columns under the lasso, the mean-shift outlier model, and the elastic net (code). Solutions are in [solutions.md](solutions.md); every number in them is asserted in `project.py`.

## 7. Coding project

`project.py` implements everything from scratch (no statsmodels, no sklearn) in ten parts. (1) OLS by normal equations, QR and SVD: agreement on a well-conditioned design and, at $\kappa=10^8$, normal-equations error 0.6 against $4\times10^{-9}$ for QR. (2) The full inference table and Examples 4.1–4.2. (3) $10^4$ replications: coverage of coefficient, mean-response and prediction intervals (all within 0.004 of 0.95), unbiasedness of $\hat\sigma^2$, the null distribution of $R^2$, the $F$ test as an LRT. (4) Leverage, studentized residuals, Cook's distance against brute-force refits, on data with two planted outliers and a high-leverage point. (5) Sandwich versus naive standard errors and WLS. (6) VIF and the SVD decomposition of variance. (7) $C_p$, AIC, BIC and LOO across polynomial degrees, and forward selection on pure noise. (8) Ridge: SVD form against the direct solve, MSE against $\lambda$. (9) Lasso by coordinate descent against an independent solver, the path, sparse recovery, and ridge versus lasso. (10) Every number used in the exercises.

**Build it yourself first:** write `lasso_cd` from (3.46) (about ten lines) and check it on the orthonormal design before reading Part 9; then write the five-line `influence` function from (3.21), (3.23) and (3.27) and compare it with `brute_loo`. "Your turn" asks for WLS with estimated weights, the elastic net, and a stepwise variation; reference solutions for the interaction/FWL (B8), omitted-variable-bias, stepwise (C1) and elastic-net (C4) exercises are at the bottom. You should see: normal-equations error 0.6; naive coverage 0.894 vs HC3 0.943; ridge MSE falling from 12.68 to 7.15 on the correlated design; and a CV-lasso that finds all 5 true predictors among 50.

## 8. Where you'll struggle

- **"Controlling for" is what FWL says and nothing more.** A coefficient is the slope on the residualized predictor. Whether it is a causal effect is a separate claim that needs a causal argument (§3.11, ST-9). Example 4.2: the coefficient moved from 1.74 to 0.39 by adding two columns; neither is "the effect."
- **The normal equations are the wrong way to compute OLS** (LA-5, LA-13), and every student does it anyway because (3.5) is on the page. Use QR or the SVD; `np.linalg.lstsq` does. Run Part 1 and look at the 0.6.
- **$R^2$ is not a measure of model quality.** It rises when you add noise, is high for a wrong model with a trend, and is meaningless without an intercept. It tells you the angle in §3.2 and nothing about whether the model is right.
- **Confidence interval versus prediction interval.** One is for the mean at $\mathbf x_0$, one for a new observation; (3.16) and (3.17) differ by the "1" that never shrinks. If you shrink the first with more data and think you can predict individuals that well, you have made this mistake.
- **Leverage, outliers and influence are three different things.** Leverage is about $x$, outlyingness is about $y$, influence needs both. A high-leverage point can have a tiny residual precisely because it bent the fit (Example 4.4).
- **Multicollinearity does not bias.** It inflates variance, and only the SVD picture (3.32) makes that obvious: a direction in coefficient space the data cannot see. Dropping a collinear predictor to "fix" the standard errors trades variance for omitted-variable bias.
- **Ridge and lasso need standardized predictors,** and everyone forgets once. Without it the penalty is decided by your choice of units. Also center $\mathbf y$ and do not penalize the intercept.
- **Post-selection inference is a swamp this module only points at.** Forward selection on pure noise gave $R^2=0.34$ and half the kept coefficients "significant." CV-chosen lasso keeps many false positives (22 of 45 nulls in the project). Selecting and testing on the same data does not work.
- **The lasso's zeros are the corners of the $\ell^1$ ball** (the kink of $|b|$), not a feature-selection oracle. With correlated or duplicated predictors it picks arbitrarily among them.
- **Harder than it looks:** B8 (interaction and FWL) and C3 (the mean-shift model), because each needs FWL applied to a column that is not a "real" predictor.

## 9. Check yourself

1. State and derive the normal equations in two ways. Why is $\mathbf e\perp C(X)$?
2. Prove $\operatorname{tr}H=p$ and $H^2=H$, and say what each tells you about degrees of freedom.
3. Why is $R^2=\cos^2\theta$? Why can adding a column never lower $R^2$, and why does adjusted $R^2$ not have that flaw?
4. State the Frisch–Waugh–Lovell theorem. What does "holding the others fixed" mean in terms of projections?
5. Prove that $\hat{\boldsymbol\beta}$ is unbiased with variance $\sigma^2(X^\top X)^{-1}$, and that $\hat\sigma^2=\mathrm{RSS}/(n-p)$ is unbiased.
6. State Gauss–Markov. Which assumptions does it use, and name two ways an estimator can beat OLS.
7. Why are $\hat{\boldsymbol\beta}$ and $\mathrm{RSS}$ independent, and why does $\mathrm{RSS}/\sigma^2$ have $n-p$ degrees of freedom? Which two tools from earlier modules give the proof?
8. Write the confidence interval for the mean response and the prediction interval. What does each one's width do as $n\to\infty$?
9. Why is $\operatorname{Var}(e_i)=\sigma^2(1-h_{ii})$? What does that say about residuals at high-leverage points?
10. Give the leave-one-out identity for the prediction error and the formula for Cook's distance. What are its two factors?
11. How does the SVD of $X$ show that nearly collinear predictors have large variance but no bias? What is the variance of $\mathbf v_p^\top\hat{\boldsymbol\beta}$?
12. Write the ridge solution in the SVD basis. By what factor does it shrink the $i$th direction? Prove that some $\lambda>0$ beats OLS in MSE.
13. Derive soft-thresholding. Why does the lasso give exact zeros and ridge does not? Give the optimality conditions (3.45) and the formula for $\lambda_{\max}$.

## 10. Pair with

- **Agresti, *Foundations of Linear and Generalized Linear Models*, ch. 1–3.** The cleanest textbook statement of the projection and distribution theory (Theorems 3.2–3.6) with careful treatment of the $F$ test and of diagnostics; weaker on the SVD and on regularization.
- **Hastie, Tibshirani and Friedman, *The Elements of Statistical Learning* (free), ch. 3.** The chapter on linear methods for regression: best-subset, ridge and lasso with the SVD picture, the lasso path, LARS, and the effective degrees of freedom. Better than these notes on the geometry of the paths and on what happens at scale.
- **Stanford STATS 305A course notes (online).** Linear models from the projection viewpoint, with Gauss–Markov and the distribution theory in the same style as §3.2–3.5.
- **Wasserman, *All of Statistics*, ch. 13.** A brisk, compact treatment of regression, prediction and model selection, with the connection to the bootstrap and to AIC/BIC.
- **Gelman, Hill and Vehtari, *Regression and Other Stories*, ch. 6–12.** The applied judgment this module lacks: interpreting coefficients, residual plots, transformations, interactions, and what "adjusting for" can and cannot do. Read it for §3.6 and §3.11.
- **Strang, MIT 18.06, lectures 15–16.** Projections and least squares in Strang's voice; the best 90 minutes for the picture in §2.
