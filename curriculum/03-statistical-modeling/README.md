# 03 · Statistical Modeling

Nine modules, about 108 hours. Inference first (what an estimator is, what a test is), then linear models and their generalization, then the Bayesian route and its computation, then time series, and finally the module that makes the rest honest: causal inference.

| # | Module | Hours | Needs | The one idea |
|---|---|---|---|---|
| 1 | [Inference foundations](01-inference-foundations/notes.md) | 12 | F5, LA-14 | An estimator is a random variable; MLE; Fisher information bounds everything (Cramér–Rao) |
| 2 | [Intervals and tests](02-intervals-and-tests/notes.md) | 12 | ST-1 | A confidence interval is a procedure; a p-value is surprise under the null, not a probability of the null |
| 3 | [Linear regression](03-linear-regression/notes.md) | 14 | ST-2, LA-5, LA-9 | Least squares is projection; Gauss–Markov; ridge shrinks along small singular values; lasso's corners give zeros |
| 4 | [Generalized linear models](04-generalized-linear-models/notes.md) | 12 | ST-3 | One exponential family, one algorithm (IRLS = Newton = weighted least squares) |
| 5 | [Bayesian inference I](05-bayesian-inference-1/notes.md) | 12 | ST-4, F5 | Posterior ∝ likelihood × prior; conjugacy is hyperparameter bookkeeping; evidence penalizes complexity |
| 6 | [Bayesian computation](06-bayesian-computation/notes.md) | 14 | ST-5, F5, CM-8 | Detailed balance designs samplers; HMC is Hamiltonian dynamics on $-\log p$; hierarchical funnels |
| 7 | [Time series I](07-time-series-1/notes.md) | 12 | ST-3, LA-7 | Stationarity; AR stability is an eigenvalue condition; forecasting is conditional expectation; unit roots |
| 8 | [Time series II](08-time-series-2/notes.md) | 12 | ST-7, F3, F5 | Spectrum = Fourier transform of the ACF; the Kalman filter is Gaussian conditioning, recursively |
| 9 | [Causal inference](09-causal-inference/notes.md) | 8 | ST-4, F4 | Association ≠ ATE; randomization buys exchangeability; DAGs tell you what to adjust for and what not to |

## Bridges out of this course

- **ST-1, ST-4 → SM-4, SM-7:** exponential families are maximum-entropy models; $\ln Z$ and its derivatives are Fisher information.
- **ST-3 → ML-3, ML-4:** regression as prediction; ridge is a Gaussian prior; Gaussian processes are Bayesian regression with infinitely many features.
- **ST-6 → SM-6, ML-10:** Metropolis, Gibbs, variational inference, the ELBO.
- **ST-8 → ML-4:** a Gaussian process with a Matérn kernel is a state-space model.
- **ST-9 → everything:** which coefficients mean anything.
