# ST-3 Solutions

Numbering mirrors `exercises.md`. Every number below is asserted in `project.py` (Part 10, or the reference cells under "Your turn").

## A. Warm-up

**A1.** $\bar x=1.5$, $\bar y=3$. Deviations: $x-\bar x=(-1.5,-0.5,0.5,1.5)$, $y-\bar y=(-2,-1,1,2)$. So $S_{xx}=2.25+0.25+0.25+2.25=5$, $S_{xy}=3+0.5+0.5+3=7$, $S_{yy}=4+1+1+4=10$. By (3.1): $\hat\beta_1=7/5=1.4$, $\hat\beta_0=3-1.4(1.5)=0.9$.

Fitted values $0.9,\,2.3,\,3.7,\,5.1$; residuals $\mathbf e=(0.1,\,-0.3,\,0.3,\,-0.1)$. Check: $\sum e_i=0$; $\sum x_ie_i=0(0.1)+1(-0.3)+2(0.3)+3(-0.1)=0$ ✓.

$\mathrm{TSS}=S_{yy}=10$. $\mathrm{RSS}=0.01+0.09+0.09+0.01=0.2$. $\mathrm{ESS}=\hat\beta_1^2S_{xx}=1.96\cdot5=9.8$, and $9.8+0.2=10$ ✓ (3.6). $R^2=9.8/10=0.98$. Also $r^2=S_{xy}^2/(S_{xx}S_{yy})=49/50=0.98$ ✓.

**A2.** (a) $h_{ii}=\tfrac14+(x_i-1.5)^2/5$: the squared deviations are $(2.25,0.25,0.25,2.25)$, divided by 5 they are $(0.45,0.05,0.05,0.45)$, so $\mathbf h=(0.7,\,0.3,\,0.3,\,0.7)$ and $\sum h_{ii}=2=p$ ✓.

(b) $1-h_{ii}=(0.3,0.7,0.7,0.3)$, so the LOO errors $e_i/(1-h_{ii})$ are $(0.1/0.3,\ -0.3/0.7,\ 0.3/0.7,\ -0.1/0.3)=(0.3333,\,-0.4286,\,0.4286,\,-0.3333)$.

By hand for observation 1: delete $(0,1)$ and fit $(1,2),(2,4),(3,5)$. Then $\bar x=2$, $\bar y=11/3$, $S_{xx}=2$, $S_{xy}=(-1)(2-\tfrac{11}3)+(0)+(1)(5-\tfrac{11}3)=\tfrac53+\tfrac43=3$. Slope $1.5$, intercept $\tfrac{11}3-3=\tfrac23$. The prediction at $x=0$ is $\tfrac23$, and the error is $1-\tfrac23=\tfrac13=0.3333$ ✓.

(c) $\mathrm{PRESS}=2\big(\tfrac19\big)+2\big(\tfrac9{49}\big)=0.2222+0.3673=0.5896$, and PRESS$/n=0.1474$. The in-sample $\mathrm{RSS}/n=0.05$, so the honest error estimate is about three times the training error.

Inflation factors $1/(1-h_{ii})$: $3.33$ for each end point, $1.43$ for each interior point. The end points have leverage 0.7: the line is pulled toward them, so their residuals ($\pm0.1$) say little about how well the *other* data predict them. Dividing by $1-h_{ii}$ undoes the pull. The interior points sit near the pivot of the line and barely move it, so their residuals are already close to their LOO errors.

**A3.** $\hat\sigma^2=45/17=2.6471$, $\hat\sigma=1.6270$, so $t_{17,0.975}\hat\sigma=2.1098\cdot1.6270=3.4326$.

- Confidence interval for the mean (3.16): $3.4326\sqrt{0.1}=3.4326\cdot0.31623=\mathbf{1.0855}$.
- Prediction interval (3.17): $3.4326\sqrt{1.1}=3.4326\cdot1.04881=\mathbf{3.6002}$.

As $n\to\infty$, $h_0=\mathbf x_0^\top(X^\top X)^{-1}\mathbf x_0=O(1/n)\to0$, so the first half-width tends to 0 and the second to $z_{0.975}\sigma=1.96\sigma$ (also $t\to z$ and $\hat\sigma\to\sigma$). The prediction interval keeps the noise of the one new observation, which no amount of data removes.

**A4.** Lasso (3.44), $\lambda=1$: $S_1(2.5)=1.5$; $S_1(-0.7)=0$ because $|-0.7|\le1$; $S_1(1.1)=0.1$; $S_1(-3.2)=-2.2$. So $\hat{\mathbf b}^{\rm L}=(1.5,\,0,\,0.1,\,-2.2)$.

Ridge with objective $\|\mathbf y-X\mathbf b\|^2+\lambda\|\mathbf b\|^2$ on an orthonormal design: $\hat b_j=z_j/(1+\lambda)=z_j/2=(1.25,\,-0.35,\,0.55,\,-1.6)$.

Optimality (3.45): with $X^\top X=I$, $\mathbf x_j^\top\mathbf r=z_j-\hat b_j=(1.0,\,-0.7,\,1.0,\,-1.0)$. For the active coordinates 1, 3, 4 we need $\mathbf x_j^\top\mathbf r=\lambda\operatorname{sign}(\hat b_j)=(+1,+1,-1)$ ✓. For the inactive coordinate 2 we need $|\mathbf x_2^\top\mathbf r|=0.7\le1$ ✓. Only coordinate 2 is zero. Ridge has no zeros.

**A5.** (a) $f_i=\sigma_i^2/(\sigma_i^2+1)=16/17,\ 4/5,\ 0.25/1.25=(0.9412,\ 0.8,\ 0.2)$, and $\mathrm{df}(1)=0.9412+0.8+0.2=1.9412$.

(b) OLS: $\operatorname{Var}(\mathbf v_3^\top\hat{\boldsymbol\beta})=\sigma^2/\sigma_3^2=1/0.25=4$. Ridge (Theorem 3.14): $\sigma^2\sigma_3^2/(\sigma_3^2+\lambda)^2=0.25/1.25^2=0.16$. The variance falls by a factor of 25. The bias is $-\lambda\alpha_3/(\sigma_3^2+\lambda)=-\alpha_3/1.25=-0.8\,\alpha_3$.

(c) $p\sigma^2=1.25$ and $\|\boldsymbol\beta\|^2=1+0.25+0.25+0.04+0.01=1.55$, so $\mathrm{MSE}(\lambda)=(1.25+1.55\lambda^2)/(1+\lambda)^2$. Differentiate:
$$\frac{d\,\mathrm{MSE}}{d\lambda}=\frac{3.1\lambda(1+\lambda)^2-2(1+\lambda)(1.25+1.55\lambda^2)}{(1+\lambda)^4}=\frac{2(1.55\lambda-1.25)}{(1+\lambda)^3}.$$
It is negative for $\lambda<1.25/1.55$ and positive after, so $\lambda^*=1.25/1.55=0.8065$, in general $p\sigma^2/\|\boldsymbol\beta\|^2$. At the minimum, $\lambda^*\|\boldsymbol\beta\|^2=p\sigma^2$, so the numerator is $p\sigma^2(1+\lambda^*)$ and $\mathrm{MSE}(\lambda^*)=p\sigma^2/(1+\lambda^*)=1.25/1.8065=\mathbf{0.6920}$. Against $\mathrm{MSE}(0)=1.25$ that is a 45% reduction. (The catch from §3.9: $\lambda^*$ needs the unknown $\|\boldsymbol\beta\|^2$.)

**A6.** (a) $X(1,-1,-1)^\top=\text{col}_1-\text{col}_2-\text{col}_3=\mathbf 0$ because the two dummy columns sum to the column of ones. So $N(X)\ne\{\mathbf 0\}$, and since $N(X^\top X)=N(X)$ (LA-5 Theorem 3.11), $X^\top X$ is singular and (3.5) does not exist.

(b) $C(X)=\operatorname{span}\{\text{col}_2,\text{col}_3\}$ (the first column is their sum), so $\hat{\mathbf y}$ is the projection onto "one constant per group": the group means $(3+5)/2=4$ and $(6+8)/2=7$, i.e. $\hat{\mathbf y}=(4,4,7,7)$. The coefficient vectors that produce it satisfy $b_0+b_1=4$ and $b_0+b_2=7$: $\mathbf b=(t,\ 4-t,\ 7-t)$, a line parallel to the null vector $(1,-1,-1)$.

(c) Minimize $t^2+(4-t)^2+(7-t)^2$: derivative $6t-22=0$, so $t=11/3$ and $\mathbf b^+=(3.6667,\,0.3333,\,3.3333)$, with $\|\mathbf b^+\|=4.9666$. Check: $\mathbf b^+\cdot(1,-1,-1)=3.6667-0.3333-3.3333=0$, so it is orthogonal to $N(X)$ as LA-9 Theorem 3.10 requires.

(d) Dropping column 2 leaves $(\mathbf 1,\ \text{group-2 indicator})$ and the fit $\hat{\mathbf b}=(4,\ 3)$: the intercept is the reference group's mean (group 1, 4), and the coefficient is the *difference* of means (group 2 minus group 1, $7-4=3$). Same fitted values, but now the coefficients are identified.

## B. Core

**B1.** (a) $H^\top=X[(X^\top X)^{-1}]^\top X^\top=H$ since $(X^\top X)^{-1}$ is symmetric. $H^2=X(X^\top X)^{-1}(X^\top X)(X^\top X)^{-1}X^\top=H$. If $H\mathbf v=\mu\mathbf v$ with $\mathbf v\ne\mathbf 0$ then $\mu\mathbf v=H^2\mathbf v=\mu^2\mathbf v$, so $\mu\in\{0,1\}$. $H$ is symmetric, so it is diagonalizable with real eigenvalues, and $\operatorname{tr}H=\sum\mu=$ the number of unit eigenvalues. By Theorem 3.2(d) $\operatorname{tr}H=p$.

(b) $H_{ii}=(H^2)_{ii}=\sum_jH_{ij}H_{ji}=\sum_jH_{ij}^2=H_{ii}^2+\sum_{j\ne i}H_{ij}^2$. So $H_{ii}\ge0$ (a sum of squares), and $H_{ii}-H_{ii}^2=\sum_{j\ne i}H_{ij}^2\ge0$ gives $H_{ii}\le1$. The argument used only symmetry and idempotence, so it holds for *any* orthogonal projection matrix.

(c) Let $P=\mathbf 1\mathbf 1^\top/n$. Since $\mathbf 1\in C(X)$, $H\mathbf 1=\mathbf 1$, hence $HP=P$ and (transposing) $PH=P$. Then $(H-P)^2=H-HP-PH+P=H-P$, and $H-P$ is symmetric, so by (b) its diagonal is nonnegative: $h_{ii}-1/n\ge0$.

(d) If $h_{ii}=1$ then $\sum_{j\ne i}H_{ij}^2=H_{ii}-H_{ii}^2=0$, so $H_{ij}=0$ for $j\ne i$ and $\hat y_i=\sum_jH_{ij}y_j=y_i$ for every $\mathbf y$. The fit passes exactly through point $i$ whatever its response is, $e_i\equiv0$, and indeed $\operatorname{Var}(e_i)=\sigma^2(1-h_{ii})=0$ (3.22). One degree of freedom of the model is spent on that single observation.

(e) $\sum h_{ii}=\operatorname{tr}H=p$, so the mean leverage is $p/n$; with (b) and (c), $1/n\le h_{ii}\le1$. ∎

**B2.** (a) Let $\tilde{\mathbf z}=(I-H)\mathbf z$ and suppose $\tilde{\mathbf z}\ne\mathbf 0$. Then the enlarged column space is the orthogonal sum $C(X)\oplus\operatorname{span}\{\tilde{\mathbf z}\}$ (the new direction is perpendicular to $C(X)$), so the new hat matrix is $H'=H+\tilde{\mathbf z}\tilde{\mathbf z}^\top/\|\tilde{\mathbf z}\|^2$ and the new residual is $\mathbf e'=\mathbf e-\tilde{\mathbf z}(\tilde{\mathbf z}^\top\mathbf e)/\|\tilde{\mathbf z}\|^2$. Since $\mathbf e'\perp\tilde{\mathbf z}$, Pythagoras gives $\|\mathbf e\|^2=\|\mathbf e'\|^2+(\tilde{\mathbf z}^\top\mathbf e)^2/\|\tilde{\mathbf z}\|^2$, i.e.
$$\mathrm{RSS}_{\rm new}=\mathrm{RSS}-\frac{(\tilde{\mathbf z}^\top\mathbf e)^2}{\|\tilde{\mathbf z}\|^2}\le\mathrm{RSS}.$$
(Also $\tilde{\mathbf z}^\top\mathbf e=\mathbf z^\top(I-H)\mathbf e=\mathbf z^\top\mathbf e$.) TSS does not change, so $R^2$ cannot fall; it stays put only if $\mathbf z^\top\mathbf e=0$.

(b) Because $\mathbf 1^\top\mathbf e=0$, the mean of $\hat{\mathbf y}$ is $\bar y$. So $\operatorname{corr}(\mathbf y,\hat{\mathbf y})=\dfrac{(\mathbf y-\bar y\mathbf 1)^\top(\hat{\mathbf y}-\bar y\mathbf 1)}{\|\mathbf y-\bar y\mathbf 1\|\,\|\hat{\mathbf y}-\bar y\mathbf 1\|}$. Write $\mathbf y-\bar y\mathbf 1=(\hat{\mathbf y}-\bar y\mathbf 1)+\mathbf e$; $\mathbf e$ is orthogonal to $\hat{\mathbf y}-\bar y\mathbf 1\in C(X)$, so the numerator is $\|\hat{\mathbf y}-\bar y\mathbf 1\|^2$ and the correlation is $\|\hat{\mathbf y}-\bar y\mathbf 1\|/\|\mathbf y-\bar y\mathbf 1\|=\sqrt{\mathrm{ESS}/\mathrm{TSS}}$. Squaring gives $R^2$ ∎.

(c) From Example 4.1, $\mathbf e=(-0.8,0.6,1.0,-0.6,-0.2)$, $\mathrm{RSS}=2.4$. $H$ projects onto $\operatorname{span}\{\mathbf 1,\ \mathbf x-3\mathbf 1\}$ with $\mathbf x-3\mathbf 1=(-2,-1,0,1,2)$, $\|\cdot\|^2=10$. For $\mathbf z=(0,1,-1,0,0)$: $\mathbf z^\top\mathbf 1=0$, $\mathbf z^\top(\mathbf x-3\mathbf 1)=-1$, so $H\mathbf z=-\tfrac1{10}(\mathbf x-3\mathbf 1)=(0.2,0.1,0,-0.1,-0.2)$ and $\tilde{\mathbf z}=(-0.2,\,0.9,\,-1,\,0.1,\,0.2)$ with $\|\tilde{\mathbf z}\|^2=0.04+0.81+1+0.01+0.04=1.9$. Also $\mathbf z^\top\mathbf e=0.6-1.0=-0.4$. By (a): $\mathrm{RSS}_{\rm new}=2.4-0.16/1.9=2.4-0.0842=\mathbf{2.3158}$. Now $p=3$, so $\hat\sigma^2=2.3158/2=\mathbf{1.1579}$ (up from $0.8$), $R^2=1-2.3158/6=\mathbf{0.6140}$ (up from $0.6$), and $R^2_{\rm adj}=1-1.1579/1.5=\mathbf{0.2281}$ (down from $0.4667$). The extra column soaked up 0.084 of RSS and cost one degree of freedom.

(d) $R^2_{\rm adj}=1-\dfrac{\mathrm{RSS}/(n-p)}{\mathrm{TSS}/(n-1)}$ rises iff $\mathrm{RSS}_{\rm new}/(n-p-1)<\mathrm{RSS}/(n-p)$. For a single added column $F=t^2$ (after (3.18)), so $t^2=\dfrac{\mathrm{RSS}-\mathrm{RSS}_{\rm new}}{\mathrm{RSS}_{\rm new}/(n-p-1)}$, i.e. $\mathrm{RSS}=\mathrm{RSS}_{\rm new}\big(1+t^2/(n-p-1)\big)$. Substitute:
$$\frac{\mathrm{RSS}_{\rm new}}{n-p-1}<\frac{\mathrm{RSS}_{\rm new}\,(1+t^2/(n-p-1))}{n-p}\iff n-p<(n-p-1)+t^2\iff t^2>1.\ ∎$$
In (c), $t^2=0.0842/(2.3158/2)=0.0727$, $|t|=0.27<1$, so $R^2_{\rm adj}$ fell ✓.

**B3.** (a) The full fit is $\mathbf y=X_1\hat{\boldsymbol\beta}_1+X_2\hat{\boldsymbol\beta}_2+\mathbf e$ with $X_1^\top\mathbf e=\mathbf 0$ (Theorem 3.2(e)). Multiply by $(X_1^\top X_1)^{-1}X_1^\top$:
$$\hat{\boldsymbol\beta}_1^{\rm short}=(X_1^\top X_1)^{-1}X_1^\top\mathbf y=\hat{\boldsymbol\beta}_1+\hat\Gamma\hat{\boldsymbol\beta}_2+(X_1^\top X_1)^{-1}X_1^\top\mathbf e=\hat{\boldsymbol\beta}_1+\hat\Gamma\hat{\boldsymbol\beta}_2.\ ∎$$

(b) Example 4.2 data. Full fit on $(\mathbf 1,\mathbf x_1,\mathbf x_2)$: slope on $\mathbf x_1$ is $0.3103$, on $\mathbf x_2$ is $0.8933$. The slope entry of $\hat\Gamma$ is the slope of $\mathbf x_2$ regressed on $(\mathbf 1,\mathbf x_1)$: $0.8214$. The short fit on $(\mathbf 1,\mathbf x_1)$ gives slope $1.0440$. Check: $0.3103+0.8214\cdot0.8933=0.3103+0.7337=1.0440$ ✓.

(c) With $X$ fixed, $\hat\Gamma=\Gamma$ is not random, and $\mathbb E\hat{\boldsymbol\beta}_1=\boldsymbol\beta_1$, $\mathbb E\hat{\boldsymbol\beta}_2=\boldsymbol\beta_2$ (Theorem 3.5(a)). Taking expectations in (a) gives $\mathbb E\hat{\boldsymbol\beta}_1^{\rm short}=\boldsymbol\beta_1+\Gamma\boldsymbol\beta_2$, which is (3.37). $\Gamma$ depends only on the design. If $\mathbf x_2$ stays correlated with $\mathbf x_1$ as data accumulate, $\Gamma$ tends to a nonzero constant, so the bias $\gamma\beta_2$ does not shrink; only the variance does. The short-model interval therefore concentrates around the wrong number and its coverage tends to 0.

**B4.** (a) $\mathbb E[\mathbf a^\top\mathbf y]=\mathbf a^\top X\boldsymbol\beta=\mathbf c^\top\boldsymbol\beta$ for all $\boldsymbol\beta$ means $(X^\top\mathbf a-\mathbf c)^\top\boldsymbol\beta=0$ for all $\boldsymbol\beta$, hence $X^\top\mathbf a=\mathbf c$.

(b) $H\mathbf a=X(X^\top X)^{-1}X^\top\mathbf a=X(X^\top X)^{-1}\mathbf c$. Then $(H\mathbf a)^\top\mathbf y=\mathbf c^\top(X^\top X)^{-1}X^\top\mathbf y=\mathbf c^\top\hat{\boldsymbol\beta}$. From $\mathbf a=H\mathbf a+(I-H)\mathbf a$: $\mathbf a^\top\mathbf y=\mathbf c^\top\hat{\boldsymbol\beta}+\big((I-H)\mathbf a\big)^\top\mathbf y$.

(c) $\operatorname{Cov}\big(\mathbf c^\top\hat{\boldsymbol\beta},\,((I-H)\mathbf a)^\top\mathbf y\big)=\sigma^2(H\mathbf a)^\top(I-H)\mathbf a=\sigma^2\mathbf a^\top H(I-H)\mathbf a=0$. Variances of uncorrelated terms add, and $\operatorname{Var}(\mathbf u^\top\mathbf y)=\sigma^2\|\mathbf u\|^2$:
$$\operatorname{Var}(\mathbf a^\top\mathbf y)=\operatorname{Var}(\mathbf c^\top\hat{\boldsymbol\beta})+\sigma^2\|(I-H)\mathbf a\|^2\ge\operatorname{Var}(\mathbf c^\top\hat{\boldsymbol\beta}).\ ∎$$

(d) $\mathbf a=(-\tfrac14,0,0,0,\tfrac14)$, $\mathbf c=(0,1)$. $X^\top\mathbf a=(\sum a_i,\ \sum a_ix_i)=(0,\ (5-1)/4)=(0,1)=\mathbf c$ ✓ unbiased for $\beta_1$. Its variance is $\sigma^2\mathbf a^\top\mathbf a=\sigma^2(\tfrac1{16}+\tfrac1{16})=0.125\sigma^2$. OLS: $\sigma^2/S_{xx}=0.1\sigma^2$. Here $H\mathbf a=X(X^\top X)^{-1}\mathbf c=(\mathbf x-3\mathbf 1)/10=(-0.2,-0.1,0,0.1,0.2)$ (the OLS weights), so $(I-H)\mathbf a=\mathbf a-H\mathbf a=(-0.05,\,0.1,\,0,\,-0.1,\,0.05)$ with $\|\cdot\|^2=0.0025+0.01+0.01+0.0025=0.025$. Confirm (c): $0.1+0.025=0.125$ ✓.

(e) Leave the class: (i) biased estimators, such as ridge (§3.9, Theorem 3.15), which trade bias for variance; (ii) nonlinear estimators, such as the median or LAD regression under heavy-tailed errors. Under Normal errors, no *unbiased* estimator beats OLS (Cramér–Rao), but biased ones still can.

**B5.** (a) §3.5 gives $\max\ell=-\tfrac n2\log(2\pi e\,\mathrm{RSS}/n)$ for each model, so $-2\log\lambda=2(\ell_{\rm full}-\ell_0)=n\log(\mathrm{RSS}_0/n)-n\log(\mathrm{RSS}/n)=n\log(\mathrm{RSS}_0/\mathrm{RSS})$. From (3.18), $\mathrm{RSS}_0-\mathrm{RSS}=qF\cdot\mathrm{RSS}/(n-p)$, so $\mathrm{RSS}_0/\mathrm{RSS}=1+qF/(n-p)$ and
$$-2\log\lambda=n\log\Big(1+\frac{qF}{n-p}\Big).\ ∎$$

(b) $qF/(n-p)=7/26=0.26923$, $x:=1.26923$, $\log x=0.23841$. So $-2\log\lambda=30(0.23841)=\mathbf{7.1523}$. For $q=2$ both tails have closed forms: the exact $F_{2,26}$ tail is $x^{-(n-p)/2}=x^{-13}=\mathbf{0.0451}$ and the Wilks tail is $\chi^2_2$: $e^{-7.1523/2}=x^{-n/2}=x^{-15}=\mathbf{0.0280}$.

(c) The $F$ test is exact under (A1)–(A4): it accounts for estimating $\sigma^2$ and uses the heavier-tailed $F_{q,n-p}$. Wilks' $\chi^2_q$ is the $n\to\infty$ limit, so here it is anti-conservative ($0.028<0.045$); the ratio of the two tails is $x^{-p/2}=x^{-2}=0.62$. As $n\to\infty$ with $p,q$ fixed, $n\log(1+qF/(n-p))\to qF$ and $qF_{q,n-p}\to\chi^2_q$, i.e. the exponents $n/2$ and $(n-p)/2$ become indistinguishable relative to their size.

**B6.** (a) For $X=\mathbf 1$: $(X^\top X)^{-1}X^\top=\mathbf 1^\top/n$, so $w_{1i}=1/n$ and (3.29) gives $\sum_i(1/n)^2e_i^2=\sum e_i^2/n^2$. With $e_i=y_i-\bar y$, $\sum e_i^2=(n-1)s^2$, so $\mathrm{HC0}=(n-1)s^2/n^2$ and $\tfrac n{n-1}\mathrm{HC0}=s^2/n$, exactly the naive variance.

(b) $\bar y=3$, $\mathbf e=(-2,-1,0,3)$, $\sum e^2=14$, $s^2=14/3=4.6667$. Naive: $s^2/n=\mathbf{1.1667}$. HC0: $14/16=\mathbf{0.875}$. HC3 with $h_{ii}=1/n=0.25$: each $e_i$ is divided by $0.75$, so $\sum(e_i/0.75)^2/16=(14/0.5625)/16=24.889/16=\mathbf{1.5556}=14/9$. (Check of (a): $0.875\cdot4/3=1.1667$ ✓.)

(c) $\mathbb E\sum e_i^2=(n-1)\sigma^2$ (Theorem 3.5(c), $p=1$). So $\mathbb E[\mathrm{HC0}]=(n-1)\sigma^2/n^2=\tfrac{n-1}n\cdot\tfrac{\sigma^2}n$: biased **low** by the factor $3/4$ at $n=4$. $\mathrm{HC3}=\sum e_i^2/(n-1)^2$ has mean $\sigma^2/(n-1)=\tfrac n{n-1}\cdot\tfrac{\sigma^2}n$: biased **high** by $4/3$. The naive $s^2/n$ is exactly unbiased under homoskedasticity. Both biases vanish as $n\to\infty$; the sandwich accepts a small-sample bias to stay valid when the $\sigma_i^2$ differ.

**B7.** (a) $X^\top X=V\Sigma^2V^\top$ has eigenvalues $1\pm\rho$: $\sigma_1^2=1.9$ with $\mathbf v_1=(1,1)/\sqrt2$ and $\sigma_2^2=0.1$ with $\mathbf v_2=(1,-1)/\sqrt2$. Singular values $(1.3784,\ 0.3162)$, $\kappa(X)=\sqrt{1.9/0.1}=\sqrt{19}=\mathbf{4.3589}$.

(b) (3.32): $(X^\top X)^{-1}=\mathbf v_1\mathbf v_1^\top/1.9+\mathbf v_2\mathbf v_2^\top/0.1$, so the (1,1) entry is $\tfrac12/1.9+\tfrac12/0.1=0.2632+5=\mathbf{5.2632}$. Via the VIF: regressing $\mathbf x_1$ on $\mathbf x_2$ explains $\rho^2=0.81$ of $\|\mathbf x_1\|^2=1$, so $\mathrm{VIF}=1/(1-0.81)=5.2632$ and $\operatorname{Var}(\hat\beta_1)/\sigma^2=\mathrm{VIF}/\|\mathbf x_1\|^2$ ✓. The $\mathbf v_2$ term supplies $5/5.2632=95\%$ of it.

(c) $(1,1)^\top$ is the eigenvector for $1.9$, so $(X^\top X)^{-1}(1,1)^\top=(1,1)^\top/1.9$ and $\operatorname{Var}(\hat\beta_1+\hat\beta_2)/\sigma^2=2/1.9=\mathbf{1.0526}$. $(1,-1)^\top$ has eigenvalue 0.1, so $\operatorname{Var}(\hat\beta_1-\hat\beta_2)/\sigma^2=2/0.1=\mathbf{20}$. $(X^\top X)^{-1}=\frac1{1-\rho^2}\begin{pmatrix}1&-\rho\\-\rho&1\end{pmatrix}$, so $\operatorname{Corr}(\hat\beta_1,\hat\beta_2)=-\rho=\mathbf{-0.9}$. The data determine the **sum** well (the direction in which the columns agree) and the **difference** badly (standard deviation $\sqrt{20}\,\sigma\approx4.5\sigma$): since $\mathbf x_1\approx\mathbf x_2$, only $\beta_1+\beta_2$ moves the fit.

(d) $f_1=1.9/2.0=0.95$, $f_2=0.1/0.2=0.5$. For $\mathbf v_2^\top\hat{\boldsymbol\beta}_\lambda$: variance $\sigma^2\sigma_2^2/(\sigma_2^2+\lambda)^2=0.1\sigma^2/0.04=2.5\sigma^2$ against $\sigma^2/\sigma_2^2=10\sigma^2$ for OLS (four times smaller), and bias $-\lambda\alpha_2/(\sigma_2^2+\lambda)=-0.5\,\alpha_2$ with $\alpha_2=\mathbf v_2^\top\boldsymbol\beta=(\beta_1-\beta_2)/\sqrt2$. The MSE in this direction is $2.5\sigma^2+0.25\alpha_2^2$ against $10\sigma^2$, so ridge wins iff $\alpha_2^2<30\sigma^2$: unless the true difference is large, trading half of $\alpha_2$ for a four-fold variance cut is a good deal. The well-determined direction is nearly untouched ($f_1=0.95$).

**B8.** Reference code (also under "Your turn" in `project.py`):

```python
n = 80; xg = rng.uniform(0, 5, n); gg = (rng.uniform(size=n) < 0.5).astype(float)
yg = 1 + 0.5 * xg + 2 * gg + 0.8 * xg * gg + rng.standard_normal(n)
Xg = np.c_[np.ones(n), xg, gg, xg * gg]; bg = ols_qr(Xg, yg)
M = np.eye(n) - Xg[:, :3] @ np.linalg.solve(Xg[:, :3].T @ Xg[:, :3], Xg[:, :3].T)   # annihilates (1, x, g)
assert np.isclose((M @ Xg[:, 3]) @ (M @ yg) / ((M @ Xg[:, 3]) @ (M @ Xg[:, 3])), bg[3])   # (ii) FWL (3.9)
s0 = ols_qr(Xg[gg == 0][:, :2], yg[gg == 0])[1]; s1 = ols_qr(Xg[gg == 1][:, :2], yg[gg == 1])[1]
assert np.isclose(bg[1], s0) and np.isclose(bg[1] + bg[3], s1)                         # (i) within-group slopes
```
Output: $\hat\beta_3=0.7848=1.2895-0.5047$ (truth $0.8$).

*Proof of (i).* $(\mathbf 1,x,g,xg)$ spans the same space as the four columns $\mathbf 1-g,\ (\mathbf 1-g)\circ x,\ g,\ g\circ x$ (each column of one set is a combination of the other set's columns: $\mathbf 1-g$, $x-xg$, $g$, $xg$). The first pair is supported on the rows with $g=0$ and the second on the rows with $g=1$, so every column of one pair is orthogonal to every column of the other. The projection onto the span is then the sum of the two projections, and $\mathrm{RSS}$ splits into the RSS of a straight-line fit in each group: the fit is two separate regressions. For $g=0$ the fitted line is $\hat\beta_0+\hat\beta_1x$, for $g=1$ it is $(\hat\beta_0+\hat\beta_2)+(\hat\beta_1+\hat\beta_3)x$. Hence $\hat\beta_1$ is the group-0 slope, $\hat\beta_1+\hat\beta_3$ the group-1 slope, and $\hat\beta_3$ their difference ∎.

*What (ii) says.* "Controlling for $x$ and $g$" applied to the product column means: take the part of $x\cdot g$ that is not a linear combination of $\mathbf 1,x,g$, and regress the matching part of $\mathbf y$ on it. That residualized column is exactly "how much the slope differs by group, beyond the main effects", which is why the coefficient is the difference of slopes.

## C. Challenge

**C1.** Code: `forward_select` and the loop under "Your turn" (200 datasets, $n=100$, 50 candidates, pure noise; selection uses the column most correlated with the residual and stops by the chosen rule).

(a)

| Rule | mean size | mean $R^2$ | share of kept coefficients with $p<0.05$ |
|---|---|---|---|
| AIC (penalty $2k$) | 10.64 | 0.327 | 0.500 |
| BIC (penalty $k\log n$) | 2.08 | 0.146 | 1.000 |
| enter while $p<0.05$ | 2.79 | 0.155 | 0.993 |

(b) Selecting on rows 1–50 and refitting on rows 51–100 gives shares $0.050$ (AIC), $0.066$ (BIC), $0.044$ (p-value rule): the nominal 5%.

(c) At every step the procedure picks the best of the remaining candidates, i.e. the largest $|t|$ among up to 50 null statistics, and then reports that winner's $p$-value as if it had been the one predictor chosen in advance. If the 50 null predictors were independent, the chance that at least one has $p<0.05$ is $1-0.95^{50}=0.9231$, so the very first step usually finds a "significant" predictor in pure noise. Adding the $k$th predictor lowers $n\log(\mathrm{RSS}/n)$ by $n\log(1+t^2/(n-k-1))$, where $t$ is the new coefficient's statistic in the enlarged model (intercept plus $k$ predictors, $n-k-1$ residual df; this is (3.18) with $q=1$, $F=t^2$). So with $n=100$ and $k\approx2$, about 97 residual df, the rules are thresholds on $|t|$: BIC needs $|t|>2.14$, the $p$-value rule $|t|>1.98$, AIC only $t^2>1.96$ ($|t|>1.40$). The maximum of 50 null statistics clears the first two routinely, so everything they keep looks significant (1.000, 0.993); AIC's bar is lower, so it keeps many more terms and about half look significant; its $R^2$ of 0.33 is almost entirely fit to noise. In (b) the second half of the data plays no role in the selection event, so the $t$ statistics there have their nominal null distributions. The cost is power and precision: selection sees half the data, estimation sees the other half, the answer depends on the split, and the resulting $p$-values are valid for "the model the first half chose", not for any true model.

**C2.** (a) $X\mathbf b=(b_1+b_2)\mathbf x=s\mathbf x$, so the objective is $\tfrac12\|\mathbf y-s\mathbf x\|^2+\lambda(|b_1|+|b_2|)$. The triangle inequality gives $|b_1|+|b_2|\ge|s|$, with equality iff $b_1b_2\ge0$. Hence for each $s$ the best split has the same signs, and the problem reduces to the one-column lasso $\tfrac12\|\mathbf y-s\mathbf x\|^2+\lambda|s|$ with minimizer $s^*=S_\lambda(\mathbf x^\top\mathbf y)/\|\mathbf x\|^2$. Here $\mathbf x^\top\mathbf y=1.5+5+9+18=33.5$ and $\|\mathbf x\|^2=30$, so $s^*=(33.5-1)/30=\mathbf{1.0833}$. The solution set is the segment $\{b_1,b_2\ge0,\ b_1+b_2=1.0833\}$. Coordinate descent started at $(5,0)$ ends at $(1.0833,0)$ and started at $(0,5)$ ends at $(0,1.0833)$: same objective, different supports.

(b) $X^\top X=30\begin{pmatrix}1&1\\1&1\end{pmatrix}$ and $X^\top\mathbf y=(33.5,33.5)^\top$. $(X^\top X+I)\mathbf b=X^\top\mathbf y$ reads $31b_1+30b_2=33.5$ and $30b_1+31b_2=33.5$. Subtracting, $b_1-b_2=0$; then $61b=33.5$ and $b_1=b_2=\mathbf{0.5492}$ (total $1.0984$). The solution is unique because $X^\top X+I$ is positive definite.

(c) The lasso fixes the *total* weight on the duplicated pair but is indifferent to its split, so which column it "selects" is decided by tie-breaking (here the order in which coordinate descent visits them). With near-duplicates (correlation 0.99 instead of 1) the choice is decided by noise. A statement such as "the lasso chose $x_1$ over $x_2$" carries no information about which one matters; ridge, being strictly convex, shares the credit equally.

**C3.** (a) Apply Theorem 3.4 with $X_1=X$ and $X_2=\mathbf d=\boldsymbol\delta_i$, so $M_1=I-H$. Then $\mathbf d^\top M_1\mathbf d=\boldsymbol\delta_i^\top(I-H)\boldsymbol\delta_i=1-h_{ii}$ and $\mathbf d^\top M_1\mathbf y=\boldsymbol\delta_i^\top(I-H)\mathbf y=e_i$. By (3.9) the coefficient is $\hat c=e_i/(1-h_{ii})$, which is the LOO prediction error (3.25). By the proof of Theorem 3.9 its variance is $\sigma^2/\|M_1\mathbf d\|^2=\sigma^2/(1-h_{ii})$.

(b) For any $(\mathbf b,c)$: $\|\mathbf y-X\mathbf b-c\boldsymbol\delta_i\|^2=\sum_{j\ne i}(y_j-\mathbf x_j^\top\mathbf b)^2+(y_i-\mathbf x_i^\top\mathbf b-c)^2$. Choosing $c$ kills the last term for any $\mathbf b$, so $\mathbf b$ minimizes the first sum, which is the least-squares problem with row $i$ deleted: $\hat{\mathbf b}=\hat{\boldsymbol\beta}_{(i)}$, $\hat c=y_i-\mathbf x_i^\top\hat{\boldsymbol\beta}_{(i)}$ (consistent with (a) and (3.25)), and $\mathrm{RSS}_{\rm aug}=\mathrm{RSS}_{(i)}$ on $n-(p+1)$ degrees of freedom. Therefore
$$t=\frac{\hat c}{\hat\sigma_{(i)}/\sqrt{1-h_{ii}}}=\frac{e_i}{\hat\sigma_{(i)}\sqrt{1-h_{ii}}}=t_i,\qquad\hat\sigma^2_{(i)}=\frac{\mathrm{RSS}_{(i)}}{n-p-1}.$$
This is a single-coefficient $t$ statistic in a regression with $p+1$ columns, so under (A1)–(A4) it is $t_{n-p-1}$ (Corollary 3.7). ∎

(c) Example 4.4 with $y_5=3.0$: $e_5=-0.544$, $h_{55}=0.92$, so $\hat c=-0.544/0.08=-6.8$. Directly: the other four points give $\hat y=0.1+0.97x$ (residuals $0.03,-0.14,0.19,-0.08$, $\mathrm{RSS}_{(5)}=0.063$), which predicts $9.8$ at $x=10$; $3.0-9.8=-6.8$ ✓. Then $\hat\sigma^2_{(5)}=0.063/2=0.0315$, $\hat\sigma_{(5)}=0.17748$, $\mathrm{se}(\hat c)=0.17748/\sqrt{0.08}=0.6275$, $t=-6.8/0.6275=\mathbf{-10.84}$ on 2 df (two-sided $p=1-|t|/\sqrt{t^2+2}=0.0084$, using the closed form for $\nu=2$). Compare the internally standardized $r_5=-1.717$ in Example 4.4: deleting the point before scoring it is what exposes it.

**C4.** (a) Fix all coordinates but $b_j$, with partial residual $\mathbf r_j$. Up to constants the objective in $b_j$ is $\tfrac12\|\mathbf r_j-\mathbf x_jb_j\|^2+\lambda|b_j|+\tfrac12\mu b_j^2=\tfrac12(\|\mathbf x_j\|^2+\mu)b_j^2-(\mathbf x_j^\top\mathbf r_j)b_j+\lambda|b_j|$. Complete the square: $=\tfrac{k}2\big(b_j-\mathbf x_j^\top\mathbf r_j/k\big)^2+\lambda|b_j|$ with $k=\|\mathbf x_j\|^2+\mu$. Dividing by $k$, this is the scalar problem $\tfrac12(b-u)^2+(\lambda/k)|b|$ with $u=\mathbf x_j^\top\mathbf r_j/k$, whose solution (Theorem 3.16) is $S_{\lambda/k}(u)=S_\lambda(\mathbf x_j^\top\mathbf r_j)/k$ (soft-thresholding commutes with positive scaling). That is the update, and `lasso_cd(..., l2=mu)` implements it.

(b) The objective is the sum of a convex function ($\tfrac12\|\mathbf y-X\mathbf b\|^2+\lambda\|\mathbf b\|_1$) and the strongly convex $\tfrac12\mu\|\mathbf b\|^2$, so it is strictly convex and coercive: a minimizer exists and is unique. For $X=(\mathbf x,\mathbf x)$ swapping $b_1\leftrightarrow b_2$ leaves the objective unchanged, so by uniqueness $b_1=b_2=b$. For $b>0$ the objective along that line is $\tfrac12\|\mathbf y-2b\mathbf x\|^2+2\lambda b+\mu b^2$ with derivative $-2\mathbf x^\top\mathbf y+4b\|\mathbf x\|^2+2\lambda+2\mu b=0$, so
$$b=\frac{\mathbf x^\top\mathbf y-\lambda}{2\|\mathbf x\|^2+\mu}=\frac{33.5-1}{60+1}=\frac{32.5}{61}=\mathbf{0.5328}>0.$$
Compare C2: the lasso total $1.0833$ is split arbitrarily; the elastic net splits evenly (total $1.0656$, slightly less because of the $\mu$ shrinkage).

(c) Reference cell under "Your turn" (`grouped_design`, within-group correlation $0.890$, 12 runs): $\|\hat{\mathbf b}-\boldsymbol\beta\|$ is $2.24$ for $\mu=0$ (lasso), $1.52$ for $\mu=5$, $0.92$ for $\mu=20$, $0.84$ for $\mu=50$.

(d) Within a group of five near-duplicates the lasso determines only the *group total* well; how it is divided is poorly determined (the Hessian is nearly flat along the splitting directions), so the validation-chosen fit gives some members large weights and others none, and the error is large even though the totals are right. For a fixed sum, $\sum b_j^2$ is smallest when the weights are equal, so the $\mu$ term makes the splitting directions costly and the group moves together (the *grouping effect*; Zou and Hastie quantify it, and C4(b) is its extreme case). The $\ell^1$ term still shrinks unrelated null coordinates. Returns diminish ($0.92\to0.84$) because a large $\mu$ also shrinks the totals (the "double shrinkage" of the naive elastic net).
