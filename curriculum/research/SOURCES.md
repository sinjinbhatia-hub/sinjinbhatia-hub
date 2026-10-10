# Research Sources and What They Decided

Consulted before writing the roadmap. Each entry says what the source is and what ordering or scope decision it informed.

## Linear algebra

- **MIT 18.06 / 18.06SC (Strang)** — [readings and session list](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/pages/readings), [18.06SC syllabus](https://www.ocw.mit.edu/courses/18-06sc-linear-algebra-fall-2011/pages/syllabus). Three units: Ax = b and the four subspaces; least squares, determinants, eigenvalues; positive-definite matrices and applications. Informed: LA-1 through LA-9 order, orthogonality before eigenvalues, the "four fundamental subspaces" as the organizing picture.
- **Strang's 2020 "recommended order" videos** — [OCW link](https://ocw.mit.edu/courses/18-065-matrix-methods-in-data-analysis-signal-processing-and-machine-learning-spring-2018/external-resources/linear-algebra-teaching-and-mit-opencourseware-youtube_0adc42f4-57dd-4014-8974-4ab0b877d4e9/). Column space first, then the big picture, orthogonality, eigenvalues, singular values. Informed: the column picture opens LA-1.
- **MIT 18.065 (Strang), *Matrix Methods in Data Analysis, Signal Processing, and ML*** — [OCW readings](https://ocw.mit.edu/courses/18-065-matrix-methods-in-data-analysis-signal-processing-and-machine-learning-spring-2018/resources/readings/), [problem sets](https://ocw.mit.edu/courses/mathematics/18-065-matrix-methods-in-data-analysis-signal-processing-and-machine-learning-spring-2018/assignments/MIT18_065S18PSets.pdf). SVD and Eckart–Young early; least squares; eigenvalue computation; randomized methods; gradient descent; backprop; circulants and the Fourier matrix. Informed: LA-9, LA-10, LA-13, LA-14, ML-8's convolution-as-circulant link.
- **Axler, *Linear Algebra Done Right*, 4th ed. (open access)** — [author site](https://linear.axler.net/), [PDF](https://library.oapen.org/bitstream/20.500.12657/85067/1/978-3-031-41026-0.pdf). Vector spaces → finite-dimensional → linear maps → polynomials → eigenvalues → inner product spaces → operators on inner product spaces → operators on complex spaces → multilinear algebra and determinants. Informed: determinants deferred and given the alternating-multilinear treatment; LA-11, LA-12 content; proof style in LA-3/LA-4.
- **Trefethen & Bau, *Numerical Linear Algebra*** — [sample lectures PDF](https://people.maths.ox.ac.uk/trefethen/TB25.pdf), [a syllabus following it](https://staff.washington.edu/rjl/classes/am584a2016/syllabus.html). Forty lectures beginning with QR rather than Gaussian elimination; conditioning and stability; Householder; eigenvalue algorithms; Krylov. Informed: LA-13 entirely; QR early in LA-5.
- **Deisenroth, Faisal, Ong, *Mathematics for Machine Learning*** — [Cambridge listing](https://www.cambridge.org/academic/subjects/computer-science/pattern-recognition-and-machine-learning/mathematics-machine-learning?format=PB). Part I: linear algebra, analytic geometry, matrix decompositions, vector calculus, probability, continuous optimization; Part II: regression, PCA, GMMs, SVMs. Informed: LA-14 scope and the ML-2/ML-6 pairing.

## Fourier and probability

- **Stanford EE261 (Osgood)** — [course page](https://see.stanford.edu/Course/EE261). Fourier series → transform → distributions and delta → convolution, sampling, filters → DFT/FFT → multidimensional. Informed: F2/F3 order and content.
- **Stein & Shakarchi, *Fourier Analysis*, vol. 1** — chapters: genesis; basic properties; convergence; applications; transform on R; on R^d; finite Fourier analysis; Dirichlet's theorem. Informed: Fourier after inner product spaces; finite Fourier analysis as the bridge to the DFT and to QC-5.
- **Harvard Stat 110 (Blitzstein)** — [syllabus](https://stat110.hsites.harvard.edu/resource/syllabus). Sample spaces and conditioning → random variables → named distributions → multivariate → limit laws → Markov chains. Text is Blitzstein & Hwang (free 2nd ed.). Informed: F4/F5 split and content.

## Physics

- **MIT 8.223 Classical Mechanics II** — [syllabus](https://ocw.mit.edu/courses/8-223-classical-mechanics-ii-january-iap-2017/pages/syllabus). Lagrangian mechanics → conserved quantities → orbits and scattering → oscillations → Hamiltonian mechanics. **MIT 8.09 Classical Mechanics III** — [syllabus](https://ocw.mit.edu/courses/8-09-classical-mechanics-iii-fall-2014/pages/syllabus). Constraints, rigid bodies, vibrations, Hamilton–Jacobi, action-angle, perturbation theory, chaos. Informed: CM-2 through CM-8 order; CM-8 scope.
- **Taylor, *Classical Mechanics*** — conservation laws, oscillations, calculus of variations, Lagrange, central forces, non-inertial frames, rigid bodies, normal modes, Hamiltonian mechanics. Informed: primary pairing text.
- **MIT 8.07 Electromagnetism II** — [OCW](https://ocw.mit.edu/courses/8-07-electromagnetism-ii-fall-2012/) (Guth & Chen lecture notes). Electrostatics, magnetostatics, matter, time-dependent fields, Maxwell, waves, radiation, relativistic electrodynamics. Informed: EM-1 through EM-8 order; Griffiths as pairing.
- **MIT 8.044 Statistical Physics I** — confirmed in catalog; content from Schroeder/Kittel pairing and MIT's standard sequence (thermo → probability → microcanonical → canonical → quantum gases). Informed: SM order.
- **MIT 8.04 Quantum Physics I (Zwiebach, 2016)** — [syllabus](https://live.ocw.mit.edu/courses/8-04-quantum-physics-i-spring-2016/pages/syllabus). Part 1 basic concepts (wavefunction, Schrödinger, Hermitian operators, wave packets, uncertainty); Part 2 1D potentials; Part 3 scattering, 3D central potentials, angular momentum, hydrogen. Griffiths required; Shankar reference. Informed: QM-1 through QM-6.
- **Caltech Ph 125** — [catalog](https://catalog.caltech.edu/current/2025-26/ph-125-abc), [a course page](https://sites.astro.caltech.edu/~golwala/ph125ab). Shankar as text; Hilbert spaces, matrix mechanics, angular momentum, symmetries, spin-1/2, approximation methods, identical particles. Informed: QM-4's heavier formalism, QM-7, QM-8.

## Statistics

- **MIT 18.650 Statistics for Applications (2016)** — [syllabus](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/pages/syllabus). Parametric inference, MLE, method of moments, hypothesis testing, goodness of fit, regression, Bayesian statistics, PCA, GLMs. Informed: ST-1 through ST-4 order.
- **Stanford STATS 305A/305B** — [305B site](https://web.stanford.edu/class/stats305b/intro.html), [bulletin](https://bulletin.stanford.edu/courses/1254811). 305A: linear models. 305B: exponential families → GLMs → regularized GLMs → Bayes/empirical Bayes → EM → Gibbs → survival. Informed: exponential families as the unifying device between ST-4 and ST-5.
- **Gelman et al., *Bayesian Data Analysis* 3rd ed.** — Part I fundamentals (single/multiparameter, asymptotics, hierarchical); Part II model checking; Part III computation (MCMC, approximations); Part IV regression. Informed: ST-5/ST-6 split.
- **Shumway & Stoffer, *Time Series Analysis and Its Applications*** — [Springer](https://link.springer.com/book/10.1007/978-3-319-52452-8). Regression and exploratory → ARIMA → spectral → state-space. Informed: ST-7/ST-8.

## Machine learning

- **Stanford CS229** — [course site](https://cs229.stanford.edu/). Linear regression → logistic → GLMs → generative models → kernels/SVMs → neural nets → bias-variance → k-means/EM → PCA → transformers (recent offerings). Informed: ML-1 through ML-6, ML-9 as terminal.
- **Stanford CS224n** — [2026 lecture 1](https://web.stanford.edu/class/cs224n/slides_w26/cs224n-2026-lecture01-intro.pdf). Word vectors → RNNs → attention → transformers → pretraining → post-training. Informed: ML-9 content.
- **Stanford EE364a (Boyd)** — [course page](https://web.stanford.edu/class/ee364a/). Convex sets/functions → problems → duality → unconstrained/equality/interior-point methods. Informed: ML-2.

## Quantum computing

- **MIT 8.370x (edX, three parts)** — [part 1 syllabus](https://openlearninglibrary.mit.edu/courses/course-v1:MITx+8.370.1x+1T2018/syllabus/). Foundations (QM, reversible computation, measurement) → protocols and algorithms (teleportation, superdense coding, Deutsch–Jozsa, Simon, QFT/Shor, Grover, simulation) → noise, channels, QKD. Follow-on 8.371x: codes, fault tolerance, magic states, complexity. Informed: QC order and QC-7/QC-8 placement.
- **Preskill Ph219 notes** — [course page](https://preskill.caltech.edu/ph219/). Ch. 1–6 foundations and information theory; ch. 7 QEC; ch. 8 fault tolerance; ch. 9 topological QC; ch. 10 quantum Shannon theory. Informed: QC-7, QC-8 reference depth.
- **Nielsen & Chuang** — chapters confirmed: 1 intro; 2 QM; 3 CS; 4 circuits; 5 QFT; 6 search; 7 physical realization; 8 noise and operations; 9 distance measures; 10 QEC; 11 entropy; 12 quantum information theory. Primary reference.
- **Mermin, *Quantum Computer Science*** — [Cambridge](https://www.cambridge.org/9780521876582). Cbits and Qbits → general features → RSA/Shor → search → error correction → few-qubit protocols. Informed: Qbits-first pedagogy in QC-1/QC-2.
- **Fowler et al., "Surface codes" (2012)** — [arXiv:1208.0928](https://arxiv.org/pdf/1208.0928). Stabilizers on two qubits → 2D array. Informed: QC-8's surface code treatment.
- **Gottesman, UMD CMSC858G** — [course](https://quics.umd.edu/education/courses/advanced-topics-theory-computing-quantum-error-correction-and-fault-tolerance-cmsc858g-0). Stabilizer codes, QEC conditions, bounds, fault tolerance, threshold, surface codes. Informed: QC-8 scope ceiling.

## Ordering disagreements and how they were resolved

| Disagreement | Resolution |
|---|---|
| Strang (determinants mid-course, computational) vs Axler (determinants last, multilinear) | Determinants in LA-6 after linear maps, taught *both* ways: signed volume and computation first, alternating-form characterization second. |
| Griffiths (wave mechanics first) vs Townsend/Sakurai (spin first) | Wave mechanics first (QM-2/3) because the student has PDE background, then a heavy finite-dimensional formalism block (QM-4/5) before QC. |
| Trefethen & Bau (QR before LU) vs Strang (LU first) | LU in LA-2 because it's how you solve systems by hand; QR in LA-5 as the first "real" decomposition; LA-13 re-derives both numerically in Trefethen's order. |
| Stats-first vs ML-first for regression | Regression derived once in ST-3 (inference view), revisited in ML-3 (prediction view). The duplication is intentional. |
| Probability inside stats vs standalone | Standalone (F4/F5) because stat mech, stats, and ML all need it and it should be learned once. |
