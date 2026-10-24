# Review: ST-3 Linear regression

Reviewer pass started 2026-10-23. Module folder: `curriculum/03-statistical-modeling/03-linear-regression/`.
Baseline: `MPLBACKEND=Agg python3 project.py` exits 0 in 62.7 s; every number printed matches the notes/solutions.
Independent recomputation of all worked-example and solution numbers (scratch: /tmp/claude-0/st3_check.py): see items below.

Format: **ID** | file:location | severity | what is wrong | fix | status

## Findings

Second pass 2026-10-24 (previous pass stopped before recording any items; baseline re-verified: exit 0 in 65.3 s, all printed
numbers match notes/solutions; independent recomputation in /tmp/claude-0/st3_check.py of every worked-example and solution
number, with numpy/scipy/sympy and no module code). Cross-references checked against the cited modules' notes: LA-5 Thm 3.11/3.13/
3.14/3.15/3.20 and (3.1)/(3.6)/(3.8)/(3.10)/(3.16); LA-9 Thm 3.9/3.10 and (3.11); LA-13 Thm 3.12; LA-2 C4; LA-1 C1/§3.6; F5 Thm
3.19/3.22/3.23/3.24; F4 Simpson; ST-1 Thm 3.6/3.15/3.28, Example 4.1; ST-2 Thm 3.1/3.9, Def 3.2, §3.7, §3.9–3.10 — all exist and
say what ST-3 attributes to them. ML-1..4, ST-7, ST-9 are forward links to roadmap entries whose topics match. History (Legendre
1805, Gauss 1809/1795/1821, Galton 1886, Fisher 1922–25, Hoerl–Kennard 1970, Tibshirani 1996, Tseng 2001, Efron et al. 2004,
Zou–Hastie 2005, Zhao–Yu 2006, Huber–White, Durbin–Watson, Box–Cox, Mallows, FWL, Sherman–Morrison, Cook) all confirmed.

**1** | notes.md Example 4.1, intervals table lead-in | low | "$w=t_{3,0.975}\hat\sigma=2.8463$"; actual $3.18245\times\sqrt{0.8}=2.8465$ (the four half-widths derived from it are correct) | change to 2.8465; add assert in project Part 2 | FIXED |
**2** | notes.md §3.3 "What just happened" | low | "wrong degrees of freedom ($n-q$ instead of $n-p$)": $q$ is not defined until §3.5 | define $q$ as the number of columns of $X_2$ at this point | FIXED |
**3** | notes.md §3.8 "More predictors than observations" | low | "(3.4) has a $(p-n)$-dimensional affine family of exact solutions with RSS=0" assumes $\operatorname{rank}X=n$, which is stated only in the next clause | make the rank hypothesis explicit: family is $(p-\operatorname{rank}X)$-dimensional, interpolates when rank is $n$ | FIXED |
**4** | notes.md §3.8 "Cross-validation" and §3.9 "Choosing λ" | medium | claims the LOO identity (3.25) holds "for any linear smoother whose fit is $A^{-1}X^\top y$ with $A$ symmetric positive definite" / "the derivation used only that". The proof of Theorem 3.8 also uses $X_{(i)}^\top X_{(i)}=A-\mathbf x_i\mathbf x_i^\top$, i.e. deleting row $i$ must subtract exactly $\mathbf x_i\mathbf x_i^\top$ from $A$; that needs $A=X^\top X+\Lambda$ with $\Lambda$ fixed (true for ridge, false for a general SPD $A$) | state the hypothesis $A=X^\top X+\Lambda$, $\Lambda$ fixed, in both places | FIXED |
**5** | notes.md §3.8 after (3.35) | low | "lowest at degree 6 (9.32) and rises by $\sigma^2=0.25$ per extra parameter after that": the increments are 0.20, 0.23, 0.25, 0.25… (bias still shrinking at degrees 7–8); exact only from degree 9 on | say "once the bias is negligible (degree 9 on), by exactly σ²" ; assert argmin = 6 and 9.32 in project | FIXED |
**6** | notes.md §3.9 Theorem 3.13 proof step 3 | trivial | "whose trace map is $\sum f_i$" | "whose trace is" | FIXED |
**7** | notes.md Example 4.7 + project.py Part 9 | medium | notes: CD "agrees to $10^{-8}$ with an independent L-BFGS-B solve"; assert uses atol=1e-5 (weaker than the claim); observed 8.1e-9 | assert atol=1e-7 and state "to better than $10^{-7}$" (robust to solver version) | FIXED |
**8** | project.py Part 3 coverage asserts | low | notes §7 claim "all within 0.004 of 0.95"; asserts test 0.005 (observed max deviation 0.0034) | tighten to 0.004 | FIXED |
**9** | project.py Part 5 | low | notes §3.6 quote naive 0.894 / HC0 0.934 / HC3 0.943 coverage, naive var 0.092 vs true 0.134, WLS 0.0117; asserts only test naive<0.92, |HC3−0.95|<0.015, var_wls<true/5 | add asserts on the quoted values (seeded, deterministic) | FIXED |
**10** | project.py Part 7 | low | notes quote best degree 6 (EPE 9.32), AIC mean degree 13.8 at n=20, 7.97/6.65/7.99 at n=1280, forward-AIC 10.8 / R²=0.34 / 52%; asserts only test loose inequalities | add asserts on the quoted values | FIXED |
**11** | project.py Part 8 | low | notes §7 quote "ridge MSE falling from 12.68 to 7.15"; Example 4.6 quotes MC 1.509; asserts test min<0.7·OLS and 5% | add asserts on 12.68, 7.15, 1.509 | FIXED |
**12** | project.py Part 9 | low | notes §3.10/§7/§8 quote ridge/lasso errors 1.68/0.81, 1.54/1.63, 3.40/1.86, and CV-lasso support 27 (22 false positives), error 0.40 vs OLS 1.20; asserted only as inequalities. Also `compare()` reports "missed"/"false+" against the sparse truth `bs` even for the dense truths, so those two printed numbers are meaningless for the dense rows | assert the quoted values; compute support errors against the actual truth `bt` | FIXED |
**13** | project.py "Your turn" (B8, C1 reference cells) | low | solutions quote B8: 0.7848 = 1.2895 − 0.5047; C1 table 10.64/0.327/0.500, 2.08/0.146/1.000, 2.79/0.155/0.993, split shares 0.050/0.066/0.044; not asserted or only loosely | add asserts | FIXED |
**14** | project.py `hat_diag` docstring | trivial | cites "(notes (3.5) …)" for $h_{ii}$; (3.5) is $\hat\beta$ | cite (3.21) and (3.41) | FIXED |
**15** | project.py Example 4.2 comment | trivial | "# df = n - p = 4, not n - 3": the residualized one-column regression would use $n-1=7$ | fix comment | FIXED |
**16** | solutions.md C1(c) | low | "Adding one term to a model with $k\approx2$ others lowers … by $n\log(1+t^2/(n-k-1))$ … about 97 residual df": with $k$ others + the new term + intercept the df is $n-k-2$; 97 corresponds to $k=2$ predictors in the enlarged model | reword so $k$ counts predictors in the enlarged model (thresholds 1.40/2.14/1.98 unchanged) | FIXED |
**17** | project.py runtime | low | 65 s against STYLE's "under about 60 s" (previous pass: 62.7 s); dominated by the pure-Python coordinate-descent loops in Parts 9 and C4 | not changed: reducing the run counts would change every quoted Monte Carlo number; within "about" | NOT FIXED: 61.5 s after this pass (65.3 before); cutting the Monte Carlo run counts would change every quoted number, and "about 60 s" is met in spirit |

## Resolution

Items 1–16 fixed in notes.md (8 edits), solutions.md (1 edit) and project.py (16 edits: 14 new or tightened asserts on numbers quoted in
notes/solutions, 1 docstring, 1 comment). `MPLBACKEND=Agg python3 project.py` re-run after the fixes: exit 0, "All checks passed in 61.5 s";
printed output identical to the pre-fix run except the dense-truth support columns in the ridge-vs-lasso table (now measured against the true
dense β; not quoted anywhere). Item 17 (runtime slightly over 60 s) left as is.

STATUS: review complete — 17 findings (0 high, 2 medium, 15 low/trivial), 16 fixed, 1 not fixed (runtime); project.py exits 0.
