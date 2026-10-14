# 04 · Machine Learning

Ten modules, about 130 hours. The order is the modern consensus (Stanford CS229/CS231n/CS224n, Bishop & Bishop, Murphy): the learning problem, optimization, linear models, kernels, trees, unsupervised methods, then deep learning from scratch through transformers and generative models. Every algorithm is implemented by hand in NumPy on top of LA-14's autodiff engine; libraries appear only as cross-checks.

| # | Module | Hours | Needs | The one idea |
|---|---|---|---|---|
| 1 | [The learning problem](01-the-learning-problem/notes.md) | 10 | ST-1 | Minimize risk you can't see using risk you can; bias–variance; the curse of dimensionality |
| 2 | [Optimization](02-optimization/notes.md) | 14 | LA-14, LA-8, CM-3 | Convexity, $\kappa$, momentum's $\sqrt\kappa$, SGD's noise ball, Adam, KKT and the SVM dual |
| 3 | [Linear models](03-linear-models/notes.md) | 10 | ML-2, ST-4 | The gradient is always "prediction minus target"; softmax is Boltzmann; generative vs discriminative |
| 4 | [Kernels, SVMs, Gaussian processes](04-kernels-svm-gp/notes.md) | 12 | ML-3, LA-8, ST-5 | The representer theorem; support vectors from complementary slackness; a GP is Gaussian conditioning |
| 5 | [Trees and ensembles](05-trees-and-ensembles/notes.md) | 10 | ML-1, F5 | Averaging cuts variance by decorrelation; boosting is gradient descent in function space |
| 6 | [Unsupervised learning](06-unsupervised-learning/notes.md) | 14 | ML-3, LA-9, LA-8, F5 | EM is coordinate ascent on the free energy; PCA is the SVD; spectral clustering is normal-mode analysis of a graph |
| 7 | [Neural networks from scratch](07-neural-networks/notes.md) | 16 | ML-2, ML-3, LA-14 | Backprop is the chain rule backward; initialization is variance propagation; shapes beside every line |
| 8 | [CNNs, RNNs, residual connections](08-architectures-cnn-rnn/notes.md) | 14 | ML-7, F3, LA-7 | Convolution is a circulant matrix the DFT diagonalizes; BPTT is matrix powers; residuals keep the Jacobian near $I$ |
| 9 | [Attention and transformers](09-attention-transformers/notes.md) | 16 | ML-8, LA-12, F3 | Attention is a softmax-weighted lookup (a learned kernel smoother); a GPT built and gradient-checked from scratch |
| 10 | [Generative models and learning theory](10-generative-models-and-theory/notes.md) | 14 | ML-9, SM-7, ST-6, ML-6 | Three ways around $Z$: the ELBO, scores, energies; diffusion is learned Langevin; why deep nets generalize is open |

## Bridges into and out of this course

- **From stat mech (SM-4, SM-7):** softmax, temperature, free energy, Langevin dynamics, Hopfield = attention.
- **From statistics (ST-1 to ST-6):** estimators, exponential families, the ELBO, MCMC, regression as prediction.
- **From linear algebra (LA-9, LA-12, LA-13, LA-14):** SVD/PCA, einsum, conditioning, matrix calculus and backprop.
- **From Fourier (F3):** convolution, circulants, positional encodings, aliasing.
- **From mechanics (CM-3, CM-8):** Lagrange multipliers and KKT; Hamiltonian dynamics in HMC; residual nets as ODE steps.
