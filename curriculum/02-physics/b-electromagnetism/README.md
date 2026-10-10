# 02b · Electromagnetism

Eight modules, about 90 hours. Griffiths' order, with the vector-calculus toolkit made explicit up front and the tensor formulation at the end building on LA-12.

| # | Module | Hours | Needs | The one idea |
|---|---|---|---|---|
| 1 | [Vector calculus toolkit](01-vector-calculus-toolkit/notes.md) | 10 | multivariable calculus | Div and curl are source and circulation densities; the delta function hides in $\hat{\mathscr r}/\mathscr r^2$ |
| 2 | [Electrostatics](02-electrostatics/notes.md) | 12 | EM-1 | Coulomb → Gauss → Poisson; energy lives in $\tfrac{\varepsilon_0}{2}E^2$ |
| 3 | [Boundary-value problems](03-boundary-value-problems/notes.md) | 14 | EM-2, F2 | Uniqueness justifies images; separation of variables is Fourier series in new clothes; Legendre polynomials |
| 4 | [Electric fields in matter](04-electric-fields-in-matter/notes.md) | 8 | EM-3 | Bound charges $\sigma_b = \mathbf P\cdot\hat n$, $\rho_b = -\nabla\cdot\mathbf P$; $\mathbf D$ is bookkeeping |
| 5 | [Magnetostatics](05-magnetostatics/notes.md) | 12 | EM-4 | $\mathbf B$ curls around currents; the vector potential; bound currents mirror bound charges |
| 6 | [Electrodynamics and Maxwell's equations](06-electrodynamics-maxwell/notes.md) | 12 | EM-5 | Displacement current fixes Ampère and predicts light; Poynting; the stress tensor; gauge freedom |
| 7 | [Electromagnetic waves](07-electromagnetic-waves/notes.md) | 12 | EM-6, F3 | Plane waves, Fresnel, dispersion from driven oscillators, group velocity, waveguides |
| 8 | [Radiation and relativity](08-radiation-and-relativity/notes.md) | 10 | EM-7, LA-12 | Retarded potentials, accelerating charges radiate ($\omega^4$), and $\mathbf E$, $\mathbf B$ are one tensor |

## Bridges out of this course

- **EM-7 → QC-1:** the Jones vector of a polarized wave is a qubit; the Poincaré sphere is the Bloch sphere.
- **EM-8 → QM-1:** Larmor radiation makes the classical atom collapse in $10^{-11}$ s. That is the problem quantum mechanics was forced to solve.
- **EM-6 → QM-4:** gauge invariance becomes the local phase freedom of the wavefunction.
- **EM-3 → QM-5, QM-6:** the Legendre equation and spherical harmonics reappear as angular momentum eigenfunctions and the hydrogen atom.
- **EM-6, EM-8 → LA-12:** the Maxwell stress tensor and the field tensor are where index notation first earns its keep.
