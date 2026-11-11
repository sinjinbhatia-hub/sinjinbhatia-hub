# 00 · Foundations

Six modules, about 62 hours, filling the gaps that the target courses assume. These are not "remedial": F2/F3 (Fourier) and F5 (probability II) are full upper-division material. Take them when the dependency graph says, not all at once.

| # | Module | Hours | Needs | Take it before |
|---|---|---|---|---|
| F1 | [Mathematical language, proof, and complex numbers](01-math-language-proof-complex-numbers/notes.md) | 6 | — | LA-3 (proofs), LA-7 (complex eigenvalues) |
| F2 | [Fourier series](02-fourier-series/notes.md) | 10 | LA-5, LA-8 | EM-3, QM-3, ST-8 |
| F3 | [Fourier transform, delta functions, convolution, DFT and FFT, sampling](03-fourier-transform-dft/notes.md) | 12 | F2 | QM-2, EM-7, ML-8, QC-5 |
| F4 | [Probability I: foundations, random variables, named distributions, expectation](04-probability-1/notes.md) | 14 | F1 | F5, SM-3, ST-1 |
| F5 | [Probability II: joint distributions, multivariate Gaussian, limit theorems, conditional expectation, Markov chains, entropy](05-probability-2/notes.md) | 14 | F4, LA-8 | ST-1, SM-4, ML-1 |
| F6 | [Calculus of variations](06-calculus-of-variations/notes.md) | 6 | multivariable calculus | CM-2 |

## Why these six

- **F1** because the wall in LA-3/LA-4 is a proof-writing wall, not a linear-algebra wall, and because complex numbers appear in LA-7, F2, F3, QM, and QC and must be reflexive.
- **F2/F3** because Fourier analysis is the one topic you flagged as partial, and because it's load-bearing for E&M (boundary-value problems, waves), QM (momentum space, uncertainty), statistics (spectral analysis), ML (convolution), and quantum computing (the QFT is the DFT). Taught *after* inner product spaces so that Fourier coefficients are seen as projections.
- **F4/F5** because stat mech, statistics, and ML each need the same probability and it should be learned once, well. F5 carries the multivariate Gaussian with its conditional-distribution formula, which is linear regression, the Kalman filter, and Gaussian-process regression in one theorem.
- **F6** because Lagrangian mechanics is "minimize a functional" and most engineering curricula never teach how to minimize a functional.

## Order

F1 first (week 1). F4 alongside LA-5 to LA-9. F2 → F3 right after LA-8. F5 after LA-8 and F4. F6 right before CM-2.
