# Master Roadmap

**Student:** second-year ME + Applied Math (physics concentration). Strong through multivariable calculus, ODEs, some PDEs. Partial Fourier. Weak linear algebra. Strong Python.

**Goal:** a self-contained, derivation-first curriculum in linear algebra, supporting math, physics (classical mechanics, E&M, thermo/stat mech, quantum mechanics), statistical modeling, machine learning, and quantum computing, sequenced so every prerequisite arrives before it is needed.

This document is the plan. Each course folder holds the actual modules. Read this once, then use [README.md](README.md) as your daily index.

---

## 1. How this roadmap was built

I compared the syllabi and canonical texts that top departments actually use, then resolved the ordering disagreements between them. The sources and what each one settled are in [research/SOURCES.md](research/SOURCES.md). The decisions that shape this curriculum:

| Decision | Source of the idea | Why I adopted it |
|---|---|---|
| Start linear algebra with the **column picture** and column space, not with determinants or abstract axioms | Strang, MIT 18.06 and his 2020 "recommended order" videos | Geometric intuition first is exactly what your teaching standard asks for |
| Delay **determinants** until after subspaces and linear maps; present them as signed volume and as an alternating multilinear form | Axler, *Linear Algebra Done Right* | Determinants taught first produce students who can compute but not think |
| Teach **QR / orthogonality before eigenvalues**, and numerical LA starting from QR, not Gaussian elimination | Trefethen & Bau; Strang | Projections and orthogonality are the single most reused idea in stats, ML, and QM |
| Treat **SVD** as the capstone of core LA and the entry point to data | Strang 18.065; Deisenroth et al. | SVD is what you'll actually use in ML and signal processing |
| Put **Fourier analysis after inner product spaces** | Stein & Shakarchi; Osgood EE261 | Fourier series are orthogonal projection onto a trig basis. Learned that way, it sticks |
| Teach probability from **Blitzstein's story-proof style**, before any statistics | Harvard Stat 110 | You need probability for stat mech, stats, and ML; it should be one shared course, not three partial ones |
| **Calculus of variations** as its own short bridge module before Lagrangian mechanics | Taylor ch. 6, Gelfand & Fomin | It's the prerequisite gap most engineers hit on the way to Lagrange |
| Quantum mechanics with a **strong Dirac-notation / finite-dimensional block** (spin-first material) in addition to wave mechanics | Townsend; Shankar; MIT 8.05 (Zwiebach) | Quantum computing lives entirely in finite-dimensional Hilbert space; wave-mechanics-only QM leaves you unprepared |
| Quantum computing in **Mermin's Qbits-first order**, with Nielsen & Chuang as the reference | Mermin; N&C; MIT 8.370x; Preskill Ph219 | Qubits → gates → protocols → algorithms → noise → error correction is the consensus order |
| Statistics as **inference first, then linear models → GLMs → Bayes → time series** | MIT 18.650; Stanford STATS 305A/B; Agresti; Gelman BDA3 | Exponential families unify GLMs and Bayesian conjugacy; teach them once, use them twice |
| ML as **optimization → linear models → kernels → trees → unsupervised → deep → transformers** | Stanford CS229, CS231n, CS224n; Bishop; Murphy | This is the modern consensus; transformers are the terminal topic of the deep-learning block |
| Explicit **stat-mech ↔ ML bridge** module | Sethna; MacKay; Jaynes | Boltzmann distributions, partition functions, free energy and variational inference are the same mathematics in two costumes |

---

## 2. Gaps I found that you didn't list

You listed linear algebra, Fourier, and probability. Working backward from each target course, these also need to be filled:

| Gap | Where it bites | Where I put it |
|---|---|---|
| **Proof-writing and mathematical language** (quantifiers, induction, "show that") | LA modules 3–4 onward; all of Axler-style LA | Foundations F1 |
| **Complex numbers and complex exponentials** fluency | LA eigenvalues, Fourier, QM, quantum computing | Foundations F1 |
| **Calculus of variations** | Lagrangian mechanics, Fermat's principle, path integrals, variational ML | Foundations F6 |
| **Special functions** (Legendre, Hermite, spherical harmonics, Bessel) | E&M boundary-value problems, hydrogen atom, oscillator | Taught in place inside EM-3 and QM-3/QM-5/QM-6 |
| **Information theory basics** (entropy, KL divergence, mutual information) | Stat mech, Bayesian stats, ML loss functions | Foundations F5 |
| **Group theory, light** (SO(3), SU(2), Pauli group) | Rigid bodies, spin, stabilizer codes | Appendix in LA-12; taught in place in QM-5 and QC-8 |
| **Tensor products and index notation** | Inertia tensor, EM field tensor, multi-qubit systems | LA-12 (full module) |
| **Vector calculus fluency** (div/grad/curl theorems, delta functions, Helmholtz) | All of E&M | EM-1 is a dedicated toolkit module |
| **Numerical methods awareness** (floating point, conditioning, stability) | Every coding project in every course | LA-13 |

Not needed for this path, so deliberately excluded: measure-theoretic probability, real analysis beyond ε–δ awareness, abstract algebra beyond what's listed above, general relativity, QFT. Each gets a one-line pointer where it would naturally come next.

---

## 3. The courses

Folder names are prefixed by the order you'll *start* them. The dependency graph in §5 is what actually governs sequencing.

### 00-python (10 modules, ~87 hours)

Added after the first build round at the student's request ("I want to get really good at Python"). Every module in every other course ends with a Python project, so this course goes first and runs as the applied track beside F1 and early LA. Where a math module derives a theorem, a PY module makes you build the thing yourself before handing you the library version, and every ▶ Try it first asks you to predict output before running.

| # | Module | Hours | Prereqs | Pair with |
|---|---|---|---|---|
| PY-1 | Objects, names, numbers (floating point), control flow, functions, scope | 8 | none | Downey, *Think Python* 3e (free) ch. 1–8; Batchelder, "Facts and myths about Python names and values"; MIT 6.100A lectures 1–4 |
| PY-2 | Collections, iteration, generators, files; the data model | 9 | PY-1 | *Think Python* ch. 9–14; Ramalho, *Fluent Python* 2e ch. 2–3, 17 |
| PY-3 | Functions as values, closures, decorators, recursion, modules, the standard library | 8 | PY-2 | *Fluent Python* ch. 7, 9; Slatkin, *Effective Python* 2e items 19–26 |
| PY-4 | Classes and the object model; dual-number automatic differentiation | 10 | PY-3 | *Fluent Python* ch. 1, 11–13, 16; the Data Model reference |
| PY-5 | Errors, testing (pytest), debugging (pdb), logging, type hints | 8 | PY-4 | Okken, *Python Testing with pytest* 2e ch. 1–5; *Effective Python* ch. 10 |
| PY-6 | NumPy: arrays, memory, broadcasting, vectorization, views vs copies | 10 | PY-5 | NumPy user guide (basics, indexing, broadcasting, copies and views); VanderPlas, *Python Data Science Handbook* 2e ch. 2 |
| PY-7 | Matplotlib and scientific communication | 7 | PY-6 | Matplotlib quick start and Artist tutorial; Rougier, *Scientific Visualization* (free) ch. 1–4 |
| PY-8 | SciPy and SymPy: optimize, integrate, solve_ivp, sparse, special, stats; symbolic → numeric pipelines | 9 | PY-7 | SciPy tutorials; Johansson, *Numerical Python* 3e ch. 3, 5–9; SymPy tutorial |
| PY-9 | pandas: tables, alignment, groupby, merge, time series, tidy data | 8 | PY-6, PY-7 | McKinney, *Python for Data Analysis* 3e (free) ch. 5–8, 10–11; Wickham, "Tidy Data" |
| PY-10 | Performance (profiling, vectorization, multiprocessing, Numba) and real software (packaging, tests, CLI, git) | 10 | PY-8 | Gorelick & Ozsvald, *High Performance Python* 2e ch. 1–2, 6–7, 9; Wilson et al., "Good enough practices in scientific computing" |

### 00-foundations (6 modules, ~62 hours)

| # | Module | Hours | Prereqs | Pair with |
|---|---|---|---|---|
| F1 | Mathematical language, proof, and complex numbers | 6 | none | Hammack, *Book of Proof* (free), ch. 1–10; any complex-numbers refresher |
| F2 | Fourier series | 10 | LA-5, LA-8 | Osgood, *Lectures on the Fourier Transform*, ch. 1; Stein & Shakarchi vol. 1 ch. 1–2 |
| F3 | Fourier transform, delta functions, convolution, DFT and FFT, sampling | 12 | F2 | Osgood ch. 2–6; Stein & Shakarchi ch. 5, 7; Stanford EE261 videos |
| F4 | Probability I: foundations, random variables, named distributions, expectation | 14 | F1 | Blitzstein & Hwang, *Introduction to Probability* (free 2nd ed.) ch. 1–6; Stat 110 videos |
| F5 | Probability II: joint distributions, multivariate Gaussian, limit theorems, conditional expectation, Markov chains, entropy and KL | 14 | F4, LA-8 | Blitzstein & Hwang ch. 7–11; MacKay *Information Theory* ch. 1–2 |
| F6 | Calculus of variations | 6 | multivariable calc | Taylor, *Classical Mechanics* ch. 6; Gelfand & Fomin ch. 1 |

### 01-linear-algebra (14 modules, ~156 hours)

The spine of everything else. Modules 1–9 are the "core" that every other course depends on. Modules 10–14 are the advanced block.

| # | Module | Hours | Prereqs | Pair with |
|---|---|---|---|---|
| LA-1 | Vectors, the two pictures of Ax = b, dot product, norms, span | 8 | none | Strang *Intro to LA* ch. 1; 3Blue1Brown *Essence of LA* 1–4 |
| LA-2 | Solving systems: elimination, echelon forms, LU, inverses, rank | 8 | LA-1 | Strang ch. 2; 18.06 lectures 2–5 |
| LA-3 | Vector spaces, subspaces, independence, basis, dimension, the four fundamental subspaces | 12 | LA-2, F1 | Strang ch. 3; Axler ch. 1–2 |
| LA-4 | Linear maps: matrices as coordinates of a map, change of basis, kernel and image, rank–nullity, isomorphism | 12 | LA-3 | Axler ch. 3; Strang §8.2–8.3; 3B1B 9–13 |
| LA-5 | Orthogonality: inner products, projections, least squares, Gram–Schmidt, QR | 12 | LA-4 | Strang ch. 4; Trefethen & Bau lectures 6–8 |
| LA-6 | Determinants: signed volume, properties, cofactors, Cramer; the alternating-multilinear view | 8 | LA-4 | Strang ch. 5; Axler ch. 9 (4th ed.) |
| LA-7 | Eigenvalues and eigenvectors: diagonalization, powers, matrix exponential and ODE systems, complex eigenvalues, Markov chains | 14 | LA-6, F1 | Strang ch. 6; Axler ch. 5; 18.06 lectures 21–25 |
| LA-8 | Inner product spaces and the spectral theorem: symmetric/Hermitian, orthogonal diagonalization, quadratic forms, positive definiteness, Rayleigh quotient, unitary and normal operators | 14 | LA-5, LA-7 | Axler ch. 6–7; Strang §6.4–6.5 |
| LA-9 | SVD and low-rank approximation: derivation, geometry, Eckart–Young, PCA, pseudoinverse, condition number | 12 | LA-8 | Strang ch. 7; Strang *LA and Learning from Data* I.8–I.9; 18.065 lectures 4–7 |
| LA-10 | Decompositions atlas: LU, Cholesky, QR, Schur, spectral, SVD, polar; which to use when | 6 | LA-9 | Strang *LA & Learning from Data* Part I; Trefethen & Bau lecture 5 |
| LA-11 | Jordan form, generalized eigenvectors, minimal polynomial, Cayley–Hamilton, functions of matrices | 10 | LA-7 | Axler ch. 8; Horn & Johnson ch. 3 |
| LA-12 | Duality, bilinear forms, tensors and tensor products, Kronecker products, index notation; appendix on SO(3), SU(2) | 14 | LA-8 | Axler ch. 9 (4th ed.); Halmos *Finite-Dimensional Vector Spaces* §§13–15, 24–27 |
| LA-13 | Numerical linear algebra: floating point, conditioning, stability, Householder QR, eigenvalue algorithms, Krylov methods, randomized LA | 16 | LA-10 | Trefethen & Bau lectures 12–16, 24–30, 32–38 |
| LA-14 | Matrix calculus: gradients, Jacobians, Hessians, chain rule in matrix form, backpropagation as matrix calculus | 10 | LA-9 | Parr & Howard, *The Matrix Calculus You Need for Deep Learning*; Petersen & Pedersen, *Matrix Cookbook* |

### 02-physics (4 sub-courses, 32 modules, ~360 hours)

#### a-classical-mechanics (8 modules, ~80 hours)

| # | Module | Hours | Prereqs | Pair with |
|---|---|---|---|---|
| CM-1 | Newtonian mechanics revisited: conservation laws, and the problem constraints create | 6 | ODEs | Taylor ch. 1–4 |
| CM-2 | Least action: from Fermat to Hamilton's principle; Euler–Lagrange | 8 | F6 | Taylor ch. 6; Feynman Lectures II-19 |
| CM-3 | Lagrangian mechanics: generalized coordinates, constraints, multipliers, worked systems | 12 | CM-2 | Taylor ch. 7; MIT 8.223 |
| CM-4 | Symmetry and conservation: Noether's theorem, cyclic coordinates, the energy function | 8 | CM-3 | Taylor §7.8; Goldstein §2.6–2.7 |
| CM-5 | Central forces and the two-body problem | 10 | CM-4 | Taylor ch. 8 |
| CM-6 | Small oscillations and normal modes (the generalized eigenproblem) | 12 | CM-3, LA-8 | Taylor ch. 11 |
| CM-7 | Rigid bodies: the inertia tensor, Euler's equations, Euler angles, the spinning top | 12 | CM-3, LA-8, LA-12 | Taylor ch. 10; Goldstein ch. 4–5 |
| CM-8 | Hamiltonian mechanics: Legendre transform, phase space, Poisson brackets, Liouville, canonical transformations, Hamilton–Jacobi | 12 | CM-4 | Taylor ch. 13; Goldstein ch. 8–10; MIT 8.09 |

#### b-electromagnetism (8 modules, ~90 hours)

| # | Module | Hours | Prereqs | Pair with |
|---|---|---|---|---|
| EM-1 | Vector calculus toolkit: grad/div/curl, integral theorems, delta functions, Helmholtz theorem, curvilinear coordinates | 10 | multivariable calc | Griffiths ch. 1 |
| EM-2 | Electrostatics: Coulomb → Gauss → Poisson; potential, energy, conductors | 12 | EM-1 | Griffiths ch. 2; Purcell ch. 1–3 |
| EM-3 | Boundary-value problems: separation of variables, images, multipoles, Legendre polynomials | 14 | EM-2, F2 | Griffiths ch. 3; MIT 8.07 notes |
| EM-4 | Electric fields in matter: polarization, D, dielectrics | 8 | EM-3 | Griffiths ch. 4 |
| EM-5 | Magnetostatics: Biot–Savart, Ampère, the vector potential, magnetic materials | 12 | EM-4 | Griffiths ch. 5–6 |
| EM-6 | Electrodynamics: Faraday, Maxwell's equations, energy and momentum, gauge freedom | 12 | EM-5 | Griffiths ch. 7–8 |
| EM-7 | Electromagnetic waves: wave equation, plane waves, polarization, reflection, waveguides, dispersion | 12 | EM-6, F3 | Griffiths ch. 9 |
| EM-8 | Potentials, radiation, and relativity: retarded potentials, dipole radiation, Lorentz covariance, the field tensor | 10 | EM-7, LA-12 | Griffiths ch. 10–12 |

#### c-thermo-statistical-mechanics (7 modules, ~80 hours)

| # | Module | Hours | Prereqs | Pair with |
|---|---|---|---|---|
| SM-1 | Thermodynamics I: state variables, temperature, first law, ideal gas, heat capacities, processes | 8 | multivariable calc | Schroeder ch. 1 |
| SM-2 | Thermodynamics II: entropy, second law, engines, thermodynamic potentials, Maxwell relations, Legendre transforms | 10 | SM-1 | Schroeder ch. 4–5; Callen ch. 1–7 |
| SM-3 | Statistical foundations: microstates, multiplicity, Boltzmann entropy, Einstein solid, two-state systems, large numbers | 12 | SM-2, F4 | Schroeder ch. 2–3; Kittel & Kroemer ch. 1–2 |
| SM-4 | Canonical ensemble: Boltzmann factor derived, partition function, free energy, equipartition; stat mech as inference (maximum entropy) | 14 | SM-3, F5 | Schroeder ch. 6; Sethna ch. 3, 6; Jaynes (1957) |
| SM-5 | Grand canonical ensemble and quantum statistics: chemical potential, Bose–Einstein, Fermi–Dirac, blackbody, Debye, Fermi gases | 14 | SM-4, QM-3 | Schroeder ch. 7; Kittel & Kroemer ch. 5–7 |
| SM-6 | Interacting systems and phase transitions: Ising model, mean-field theory, Landau theory, Metropolis Monte Carlo | 14 | SM-4 | Schroeder ch. 8; Sethna ch. 8, 12 |
| SM-7 | Bridge: stat mech ↔ probabilistic ML. Softmax = Boltzmann, log-partition = free energy, variational free energy = ELBO, Langevin dynamics = diffusion models, Hopfield and Boltzmann machines | 8 | SM-4, ST-6 | MacKay ch. 31, 41–43; Sethna ch. 8 |

#### d-quantum-mechanics (9 modules, ~110 hours)

| # | Module | Hours | Prereqs | Pair with |
|---|---|---|---|---|
| QM-1 | Why quantum: the experiments that broke classical physics; the postulates previewed as a vector-space theory | 6 | LA-8 | Griffiths & Schroeter preface + ch. 1 intro; Feynman Lectures III-1 |
| QM-2 | The wavefunction and Schrödinger equation: probability, operators, expectation values, uncertainty, momentum space | 12 | QM-1, F3 | Griffiths ch. 1; MIT 8.04 (Zwiebach) |
| QM-3 | One-dimensional problems: infinite well, harmonic oscillator (analytic and ladder operators), free particle, delta and finite wells, tunneling | 16 | QM-2 | Griffiths ch. 2 |
| QM-4 | Formalism: Hilbert space, Dirac notation, Hermitian operators, eigenbases, the generalized statistical interpretation, uncertainty principle derived, time evolution, pictures | 16 | QM-3, LA-8, LA-12 | Griffiths ch. 3; Shankar ch. 1, 4; Townsend ch. 1–4; MIT 8.05 |
| QM-5 | Angular momentum and spin: the commutator algebra, ladder operators, spherical harmonics, spin-1/2, Pauli matrices, Stern–Gerlach, the Bloch sphere | 14 | QM-4 | Griffiths ch. 4.3–4.4; Townsend ch. 1–3; Shankar ch. 12, 14 |
| QM-6 | Three dimensions and hydrogen: separation of variables, the radial equation, spectrum, degeneracy | 12 | QM-5 | Griffiths ch. 4.1–4.2 |
| QM-7 | Composite systems and entanglement: tensor products, identical particles, Bell states, EPR, Bell's inequality, density matrices | 12 | QM-5, LA-12 | Griffiths ch. 5, 12; Townsend ch. 5; Nielsen & Chuang §2.3–2.6 |
| QM-8 | Approximation methods: perturbation theory, variational principle, WKB, time-dependent perturbation, Fermi's golden rule, Rabi oscillations | 14 | QM-6 | Griffiths ch. 7–9, 11 |
| QM-9 | Symmetry, the classical limit, and path integrals: Ehrenfest, Noether in QM, Feynman's sum over paths, decoherence | 8 | QM-8, CM-8 | Shankar ch. 8, 11; Feynman & Hibbs ch. 1–3 |

### 03-statistical-modeling (9 modules, ~108 hours)

| # | Module | Hours | Prereqs | Pair with |
|---|---|---|---|---|
| ST-1 | Inference foundations: models, estimators, bias/variance/MSE, consistency, sufficiency, MLE, method of moments, Fisher information, Cramér–Rao | 12 | F5 | Wasserman *All of Statistics* ch. 6–9; MIT 18.650 ch. 1–3; Casella & Berger ch. 6–7 |
| ST-2 | Intervals and tests: sampling distributions, pivots, the bootstrap, Neyman–Pearson, likelihood-ratio tests, p-values and their abuse, multiple testing | 12 | ST-1 | Wasserman ch. 8, 10; 18.650 ch. 4–5 |
| ST-3 | Linear regression: least squares as projection, Gauss–Markov, inference, diagnostics, model selection, multicollinearity, ridge and lasso | 14 | ST-2, LA-5, LA-9 | Agresti *Foundations of Linear and GLMs* ch. 1–3; ESL ch. 3; Stanford STATS 305A |
| ST-4 | Generalized linear models: exponential families, link functions, logistic and Poisson regression, IRLS, deviance, overdispersion | 12 | ST-3 | Agresti ch. 4–7; McCullagh & Nelder |
| ST-5 | Bayesian inference I: priors, posteriors, conjugacy, credible intervals, predictive distributions, Bayesian regression, model comparison | 12 | ST-4 | Gelman et al. *BDA3* (free) ch. 1–5; McElreath *Statistical Rethinking* ch. 1–4 |
| ST-6 | Bayesian computation: Metropolis–Hastings, Gibbs, Hamiltonian Monte Carlo, diagnostics, hierarchical models, variational inference | 14 | ST-5 | BDA3 ch. 10–13; McElreath ch. 9, 13; Betancourt's HMC intro |
| ST-7 | Time series I: stationarity, autocovariance, AR/MA/ARMA, Yule–Walker, forecasting, ARIMA, unit roots | 12 | ST-3 | Shumway & Stoffer ch. 1–3 |
| ST-8 | Time series II: spectral analysis, state-space models and the Kalman filter, volatility models | 12 | ST-7, F3 | Shumway & Stoffer ch. 4, 6 |
| ST-9 | Causal inference and design: randomization, confounding, DAGs, adjustment, instruments, the limits of regression | 8 | ST-4 | McElreath ch. 5–6; Hernán & Robins *What If* (free) ch. 1–7 |

### 04-machine-learning (10 modules, ~130 hours)

| # | Module | Hours | Prereqs | Pair with |
|---|---|---|---|---|
| ML-1 | The learning problem: ERM, loss functions, generalization, bias–variance, validation, k-NN, the curse of dimensionality | 10 | ST-1 | CS229 notes §1; Bishop & Bishop *Deep Learning* ch. 1; ESL ch. 2 |
| ML-2 | Optimization for ML: convexity, gradient descent with convergence proofs, SGD, momentum, Adam, Newton, Lagrangian duality and KKT | 14 | LA-14, CM-3 | Boyd & Vandenberghe (free) ch. 2–5, 9; Deisenroth *MML* ch. 7 |
| ML-3 | Linear models as ML: logistic/softmax regression, regularization, feature maps, generative vs discriminative, GDA, naive Bayes | 10 | ML-2, ST-4 | CS229 notes §2–4 |
| ML-4 | Kernels, SVMs, and Gaussian processes | 12 | ML-3 | CS229 §6; Bishop PRML ch. 6–7; Rasmussen & Williams ch. 2 |
| ML-5 | Trees and ensembles: decision trees, bagging, random forests, boosting, gradient boosting | 10 | ML-1 | ESL ch. 9–10, 15 |
| ML-6 | Unsupervised learning: k-means, Gaussian mixtures and EM derived, PCA via SVD, spectral clustering, manifold methods | 14 | ML-3, LA-9 | CS229 §10–12; Bishop PRML ch. 9, 12 |
| ML-7 | Neural networks from scratch: MLPs, universal approximation, backprop derived with matrix calculus, autodiff, initialization, normalization, regularization | 16 | ML-2, LA-14 | Bishop & Bishop ch. 6–9; Goodfellow ch. 6–8; CS231n |
| ML-8 | Architectures: convolution (and its link to Fourier), CNNs, recurrent nets, residual connections, training dynamics | 14 | ML-7, F3 | CS231n; Bishop & Bishop ch. 10; Goodfellow ch. 9–10 |
| ML-9 | Attention and transformers: attention as soft lookup, multi-head attention, positional encoding, layer norm, a GPT from scratch, pretraining and finetuning, scaling | 16 | ML-8 | Vaswani et al. 2017; CS224n; Karpathy *nanoGPT*; Bishop & Bishop ch. 12 |
| ML-10 | Generative and probabilistic deep learning: VAEs and the ELBO, diffusion models, energy-based models; learning theory sketch (PAC, VC) | 14 | ML-9, SM-7 | Bishop & Bishop ch. 19–20; Murphy *PML2* ch. 21, 25; Sohl-Dickstein 2015 |

### 05-quantum-computing (9 modules, ~94 hours)

Two entry routes. **Physics route:** after QM-7. **Math route:** after LA-12, using the self-contained postulates in QC-1 (you lose physical intuition but not correctness; most CS quantum courses do this).

| # | Module | Hours | Prereqs | Pair with |
|---|---|---|---|---|
| QC-1 | Qubits: bits vs qubits, state space, measurement, global and relative phase, the Bloch sphere, single-qubit gates as rotations | 8 | LA-8 (QM-5 ideal) | Mermin ch. 1; N&C §1.2–1.3, 4.2 |
| QC-2 | Many qubits: tensor products, entanglement, CNOT, universal gate sets, the circuit model, no-cloning | 10 | QC-1, LA-12 | Mermin ch. 1–2; N&C ch. 4 |
| QC-3 | Protocols: superdense coding, teleportation, BB84, the CHSH game | 8 | QC-2 | Mermin ch. 6; N&C §1.3.7, 12.6 |
| QC-4 | Early algorithms: Deutsch–Jozsa, Bernstein–Vazirani, Simon; query complexity and interference | 10 | QC-2 | Mermin ch. 2; N&C §1.4 |
| QC-5 | The quantum Fourier transform and Shor's algorithm: QFT, phase estimation, order finding, factoring, RSA | 14 | QC-4, F3 | Mermin ch. 3; N&C ch. 5 |
| QC-6 | Grover search and amplitude amplification | 8 | QC-4 | Mermin ch. 4; N&C ch. 6 |
| QC-7 | Noise and open systems: density matrices, quantum channels, Kraus operators, decoherence, fidelity | 10 | QC-2, QM-7 | N&C ch. 8–9; Preskill ch. 3 |
| QC-8 | Quantum error correction: repetition and Shor codes, the stabilizer formalism, Steane code, the surface code, the threshold theorem | 16 | QC-7 | Mermin ch. 5; N&C ch. 10; Gottesman's notes; Fowler et al. 2012 |
| QC-9 | Hamiltonian simulation and variational algorithms: Trotterization, VQE, QAOA, the NISQ landscape | 10 | QC-5, QM-4 | N&C §4.7; Preskill "NISQ" (2018) |

**Total: 80 modules, roughly 910 hours.**

---

## 4. Time, honestly

At 12 hours/week during term and 25 hours/week in summer, this is **about 18 months**. At 6 hours/week it's three years. There is no way to compress 900 hours of real learning; what you *can* do is choose the order so that whatever you finish is useful on its own.

| Milestone | Cumulative hours | What you can do afterward |
|---|---|---|
| LA-1 to LA-9 + F1 + F4 | ~130 | Read any applied math paper's linear algebra; start ML-1, ST-1, CM-1 |
| + F2, F3, F5, F6, LA-10 to LA-14 | ~230 | Full advanced LA; start E&M, QM; ML-2 onward |
| + CM-1 to CM-8, ST-1 to ST-4 | ~360 | Lagrangian/Hamiltonian mechanics; regression and GLMs done |
| + EM-1 to EM-8, QM-1 to QM-7 | ~560 | Undergraduate physics core complete; QC fully unlocked |
| + SM-1 to SM-7, ST-5 to ST-9 | ~690 | Bayesian + time series; stat mech; ML-10 unlocked |
| + ML-1 to ML-10 | ~820 | Build and train a transformer from scratch with understanding |
| + QC-1 to QC-9 | ~910 | Done |

---

## 5. Dependency graph (what must come before what)

The Python course is a chain PY-1 → … → PY-10 that feeds in from the side: PY-5 before LA-2's project in earnest, PY-6/7 before CM-1, PY-8 before CM-3 and ST-1, PY-9 before ST-1, PY-10 before ML-7 and SM-6.

```
F1 ──► LA-3 ──► LA-4 ──► LA-5 ──► LA-8 ──► LA-9 ──► LA-10 ──► LA-13
 │               │         │        │  │      │
 │               └► LA-6 ──┴► LA-7 ─┘  │      └► LA-14 ──► ML-2 ──► ML-7 ──► ML-8 ──► ML-9 ──► ML-10
 │                            │        │                          ▲                         ▲
 │                            └► LA-11 ├──► LA-12 ──► CM-7, EM-8, QM-4, QM-7, QC-2        │
 │                                     │                                                   │
 │                                     ├──► F2 ──► F3 ──► QM-2, EM-7, ST-8, ML-8, QC-5     │
 │                                     │                                                   │
 └──► F4 ──► F5 ──► ST-1 ──► ST-2 ──► ST-3 ──► ST-4 ──► ST-5 ──► ST-6 ──► SM-7 ───────────┘
             │                                   │
             └► SM-3 ──► SM-4 ──► SM-5, SM-6     └► ST-7 ──► ST-8 ; ST-9

F6 ──► CM-2 ──► CM-3 ──► CM-4 ──► CM-5, CM-8 ──► QM-9
                 │
                 └► CM-6, CM-7 (need LA-8, LA-12)

EM-1 ──► EM-2 ──► EM-3 ──► EM-4 ──► EM-5 ──► EM-6 ──► EM-7 ──► EM-8

QM-1 ──► QM-2 ──► QM-3 ──► QM-4 ──► QM-5 ──► QM-6 ──► QM-8 ──► QM-9
                   │                  │
                   └► SM-5            └► QM-7 ──► QC-1 ──► QC-2 ──► QC-3, QC-4 ──► QC-5, QC-6
                                                             │
                                                             └► QC-7 ──► QC-8 ; QC-9
```

---

## 6. Recommended study path

Work **two tracks at once**, never more. One "hard new abstraction" track and one "applied / familiar" track. Switching between them when stuck is more productive than grinding one.

### Phase 0 — Weeks 1–3 (≈70 h)
- PY-1 → PY-5 as the applied track (about 43 h). Light on math, heavy on typing. Do every A exercise.
- F1 (proof, complex numbers). Do every exercise. This is the cheapest insurance in the whole plan.
- LA-1, LA-2. Should feel easy. If it doesn't, slow down here rather than later.
- Then PY-6 and PY-7 during Phase 1 (before any physics project), PY-8 before CM-3, PY-9 before ST-1, PY-10 before ML-7.

### Phase 1 — Weeks 3–14 (≈130 h): Linear algebra core
- **Track A:** LA-3 → LA-4 → LA-5 → LA-6 → LA-7 → LA-8 → LA-9.
- **Track B (from week 6):** F4 Probability I.
- **Checkpoint:** the LA core self-exam in `01-linear-algebra/CHECKPOINT-core.md`. Don't continue until you can do it closed-book.

### Phase 2 — Weeks 15–26 (≈130 h): Advanced LA, Fourier, probability II, mechanics begins
- **Track A:** F2 → F3 (Fourier), then LA-10 → LA-11 → LA-12 → LA-13 → LA-14.
- **Track B:** F5 → F6 → CM-1 → CM-2 → CM-3 → CM-4.

### Phase 3 — Weeks 27–45 (≈230 h): Mechanics, E&M, statistics
- **Track A:** CM-5 → CM-6 → CM-7 → CM-8, then EM-1 → EM-4.
- **Track B:** ST-1 → ST-2 → ST-3 → ST-4, then ML-1.

### Phase 4 — Weeks 46–62 (≈200 h): E&M finish, quantum mechanics, Bayes and time series, stat mech
- **Track A:** EM-5 → EM-8, then QM-1 → QM-5.
- **Track B:** ST-5 → ST-6, SM-1 → SM-4, ST-7 → ST-8.

### Phase 5 — Weeks 63–78 (≈200 h): QM finish, ML, quantum computing, bridges
- **Track A:** QM-6 → QM-9, then QC-1 → QC-9.
- **Track B:** ML-2 → ML-10, SM-5 → SM-7, ST-9.

### If you want ML sooner
After Phase 1 plus F5 and LA-14 (which only needs LA-9), you can run ML-1 → ML-7. Deep learning modules ML-8 to ML-10 wait for F3 and SM-7 respectively. You'll understand the algorithms but will have less of the "why" from physics. That's a legitimate trade.

### If you want quantum computing sooner
After LA-8 and LA-12, take QC-1 → QC-6 via the math route. QC-7 to QC-9 want QM-4 and QM-7 for real understanding of density matrices and Hamiltonians.

---

## 7. Where you will struggle

I'd rather you hear this now than discover it at week 9.

**Linear algebra.**
- **LA-3 and LA-4 are the wall.** Every ME student I've seen hits it. You go from computing to proving in two weeks. The fix is boring: write out a proof every single day for three weeks, even trivial ones. F1 exists to soften this.
- **LA-8:** complex inner products are conjugate-linear in one slot. You *will* get the convention backwards at least once in QM. We adopt the physics convention (conjugate the first slot) throughout, and I flag it every time it appears.
- **LA-11 Jordan form** is fiddly. Most people learn it twice. Learn the *idea* (nilpotent part + diagonal part) and let the algorithm be secondary.
- **LA-12 tensors.** Three definitions circulate (multilinear map, element of a tensor product, "array that transforms like a tensor"). They're equivalent and the module proves it, but the first time through, you'll feel like you're learning three things.
- **LA-13.** Pen-and-paper intuition doesn't transfer to floating point. The distinction between *conditioning* (the problem) and *stability* (the algorithm) takes a few projects to internalize.

**Foundations.**
- **Fourier conventions.** Where the 2π goes differs by field. We fix one convention per module and state it in a box at the top. Expect to translate when reading other books.
- **Probability.** Your intuition for conditional probability is wrong and will stay wrong until you've done fifty problems. Monty Hall, base-rate neglect, and the "two children" problem exist to prove this to you. Also: computing expectations of *functions* of random variables (LOTUS) and Jacobians for transformations are mechanical skills you must drill.

**Physics.**
- **CM-2:** "varying a function" is a new kind of calculus. The δ-versus-d distinction trips everyone.
- **CM-7:** Euler angles are objectively unpleasant. Don't judge yourself by them.
- **CM-8:** canonical transformations and Hamilton–Jacobi theory are abstract and most undergraduates never master them. Aim for *Poisson brackets and Liouville fluent, Hamilton–Jacobi recognizable.*
- **E&M** is a vector-calculus endurance test. EM-3 (boundary-value problems) is pattern recognition plus patience. EM-8 relativity notation is a lot of new symbols at once.
- **Stat mech:** entropy conceptually, the "why the Boltzmann factor" derivation (SM-4 does it two ways), and the chemical potential (notoriously confusing; SM-5 spends real time on it).
- **QM:** the hardest part is emotional, not mathematical: giving up classical pictures. Mathematically, if LA-8 and LA-12 are solid, QM-4 is a *relief* rather than a shock. Angular momentum algebra (QM-5) and perturbation-theory bookkeeping (QM-8) are the technical grinds.

**Statistics.**
- The logic of hypothesis testing is genuinely subtle and most practitioners get it wrong. ST-2 is slow on purpose.
- Frequentist and Bayesian answers to the same question differ in *interpretation*, not just numbers. Holding both in your head is the skill.
- MCMC diagnostics (ST-6) are craft, not theorem. Time-series stationarity and unit roots (ST-7) are where applied people make expensive mistakes.

**Machine learning.**
- The math is the easy part once LA and probability are done. The hard part is engineering judgment: debugging training, knowing when a result is noise. The projects deliberately make you train things that fail.
- Transformers (ML-9) are 90% tensor-shape bookkeeping. Write shapes next to every line.

**Quantum computing.**
- Qubit ordering conventions (little-endian vs big-endian) will burn you. We use the N&C convention: qubit 0 is the leftmost tensor factor. Stated at the top of every module.
- Phase kickback (QC-4) is the single idea that makes all the algorithms click; if it doesn't, stop and redo it.
- The stabilizer formalism (QC-8) is a new algebra, and the surface code is a conceptual leap. Budget the full 16 hours.

---

## 8. Conventions used everywhere

- **Vectors** are column vectors; $\mathbf{x}^\top \mathbf{y}$ is the dot product. Bold lowercase for vectors, uppercase for matrices, $\langle u, v\rangle$ for abstract inner products.
- **Complex inner product** is conjugate-linear in the *first* argument (physics convention), so $\langle u, v\rangle = u^\dagger v$ and bra-ket notation $\langle u | v \rangle$ matches.
- **Fourier transform:** $\hat f(\xi) = \int f(x)\, e^{-2\pi i x \xi}\, dx$ (Osgood/Stein convention) in F3; physics modules use $\int f(x) e^{-ikx}dx/\sqrt{2\pi}$ and say so.
- **Probability:** capital letters for random variables, $\mathbb{E}$, $\mathrm{Var}$, $p(x)$ for densities.
- **Qubits:** $|0\rangle = (1, 0)^\top$, $|1\rangle = (0, 1)^\top$; multi-qubit states ordered with qubit 0 as the leftmost tensor factor.
- **Units:** SI in E&M (Griffiths), $\hbar$ kept explicit in QM, $k_B$ kept explicit in stat mech.
- **Code:** Python 3, NumPy/SciPy/Matplotlib, SymPy for symbolic checks. Every project is a `.py` with `# %%` cell markers; notebooks are generated from them.

---

## 9. What's deliberately *not* here

- Measure-theoretic probability, real/functional analysis, abstract algebra as subjects. Pointers given where they'd come next.
- General relativity, quantum field theory, condensed-matter beyond the Ising model.
- Reinforcement learning (one pointer in ML-10). Worth its own course later.
- Hardware details of quantum computers beyond what QC-9 needs.
