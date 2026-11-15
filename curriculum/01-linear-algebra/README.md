# 01 · Linear Algebra

Fourteen modules, about 156 hours. Modules 1–9 are the **core** that every other course depends on. Modules 10–14 are the **advanced block**. Take the [core checkpoint](CHECKPOINT-core.md) after module 9 before moving on.

Each module folder has `notes.md` (the lesson), `exercises.md`, `solutions.md`, and `project.py` (run with `MPLBACKEND=Agg python3 project.py`; figures land in `figures/`).

## Core

| # | Module | Hours | Needs | The one idea |
|---|---|---|---|---|
| 1 | [Vectors, the two pictures of Ax = b, dot product, norms, span](01-vectors-and-two-pictures/notes.md) | 8 | — | $A\mathbf{x}$ is a combination of the columns of $A$ |
| 2 | [Elimination, echelon forms, LU, inverses, rank](02-elimination-and-lu/notes.md) | 8 | 1 | Row operations are matrices; elimination is a factorization |
| 3 | [Vector spaces, subspaces, basis, dimension, the four fundamental subspaces](03-vector-spaces-and-subspaces/notes.md) | 12 | 2, F1 | Every matrix comes with four subspaces whose dimensions are $r, n-r, r, m-r$ |
| 4 | [Linear maps, change of basis, kernel and image, rank–nullity](04-linear-maps/notes.md) | 12 | 3 | A matrix is what a linear map does to a basis; similar matrices are the same map in different coordinates |
| 5 | [Orthogonality, projections, least squares, Gram–Schmidt, QR](05-orthogonality-and-least-squares/notes.md) | 12 | 4 | The closest point has a perpendicular error |
| 6 | [Determinants](06-determinants/notes.md) | 8 | 4 | The unique signed-volume function |
| 7 | [Eigenvalues and eigenvectors, diagonalization, matrix exponential](07-eigenvalues/notes.md) | 14 | 6, F1 | Directions a matrix only stretches make dynamics one-dimensional |
| 8 | [Inner product spaces and the spectral theorem](08-inner-product-spaces-spectral-theorem/notes.md) | 14 | 5, 7 | Symmetric matrices have orthogonal eigenvectors and real eigenvalues |
| 9 | [SVD and low-rank approximation](09-svd/notes.md) | 12 | 8 | Every matrix is rotate–stretch–rotate; the stretches tell you everything |

**Checkpoint:** [CHECKPOINT-core.md](CHECKPOINT-core.md). Closed book. If you can't do it, go back; the advanced block and every other course assume it.

## Advanced

| # | Module | Hours | Needs | The one idea |
|---|---|---|---|---|
| 10 | [Decompositions atlas](10-decompositions-atlas/notes.md) | 6 | 9 | Factor first, then do the easy thing; Schur is the key you were never given |
| 11 | [Jordan form, generalized eigenvectors, Cayley–Hamilton, matrix functions](11-jordan-form/notes.md) | 10 | 7 | Diagonal plus nilpotent chains; where $t e^{\lambda t}$ comes from |
| 12 | [Duality, bilinear forms, tensors, Kronecker products, index notation](12-duality-and-tensors/notes.md) | 14 | 8 | Three communities, one idea; entanglement is Schmidt rank > 1 |
| 13 | [Numerical linear algebra](13-numerical-linear-algebra/notes.md) | 16 | 10 | Conditioning is the problem's fault; stability is the algorithm's |
| 14 | [Matrix calculus](14-matrix-calculus/notes.md) | 10 | 9 | The derivative is a linear map; backprop is the chain rule in the cheap order |

## Where this course feeds

- **LA-5, LA-8 → F2 (Fourier series):** coefficients are projections onto an orthogonal basis of functions.
- **LA-8 → QM-4, QC-1:** observables are Hermitian operators; the spectral theorem is why measurement outcomes are well defined.
- **LA-9 → ML-6 (PCA), ST-3 (regression), LA-13 (conditioning).**
- **LA-12 → CM-7 (inertia tensor), EM-8 (field tensor), QM-7 / QC-2 (tensor-product state spaces).**
- **LA-14 → ML-2, ML-7 (optimization, backpropagation), ST-1 (Fisher information).**
- **LA-13 → every coding project afterward.**

## Pacing

Core at 12 h/week: about 9 weeks. Add F1 before module 3 and F4 alongside modules 5–9. Advanced block: about 5 weeks. The roadmap's phase plan is in [../ROADMAP.md](../ROADMAP.md#6-recommended-study-path).
