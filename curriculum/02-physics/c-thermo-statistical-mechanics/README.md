# 02c · Thermodynamics and Statistical Mechanics

Seven modules, about 80 hours. Thermodynamics as the phenomenology, then the statistical derivation, then the two places it becomes modern: quantum gases and phase transitions. The last module is the explicit bridge to probabilistic machine learning, so that softmax, free energy, variational inference, and diffusion models arrive as things you already know.

| # | Module | Hours | Needs | The one idea |
|---|---|---|---|---|
| 1 | [Thermodynamics I](01-thermodynamics-1/notes.md) | 8 | multivariable calculus | $U$ is a state function; $Q$ and $W$ are not; equipartition and its first crack |
| 2 | [Thermodynamics II](02-thermodynamics-2/notes.md) | 10 | SM-1 | Entropy exists; potentials are Legendre transforms; Maxwell relations are free |
| 3 | [Statistical foundations](03-statistical-foundations/notes.md) | 12 | SM-2, F4 | $S = k_B\ln\Omega$; the second law is overwhelming probability; Sackur–Tetrode from counting |
| 4 | [Canonical ensemble](04-canonical-ensemble/notes.md) | 14 | SM-3, F5 | $e^{-\beta E}/Z$ derived two ways; $Z$ knows everything; stat mech as inference |
| 5 | [Quantum statistics](05-quantum-statistics/notes.md) | 14 | SM-4, QM-3 | Chemical potential; Bose–Einstein and Fermi–Dirac; Planck's law; Debye |
| 6 | [Phase transitions and the Ising model](06-phase-transitions-ising/notes.md) | 14 | SM-4, F5 | Energy vs entropy tips at $T_c$; mean field, Landau, Metropolis, universality |
| 7 | [Bridge to machine learning](07-stat-mech-ml-bridge/notes.md) | 8 | SM-4, SM-6 | Softmax is Boltzmann; ELBO is minus a free energy; mean field is variational inference |

## Bridges out of this course

- **SM-4, SM-7 → ML-3, ML-10:** softmax, cross-entropy, energy-based and diffusion models.
- **SM-6 → ST-6:** Metropolis is the first MCMC algorithm; detailed balance is the shared foundation.
- **SM-5 → QM-1:** Planck's law is where quantum mechanics started.
- **SM-2 ↔ CM-8:** the Legendre transform appears in both; learn it once.
- **SM-3 ← CM-8:** Liouville's theorem is why the microcanonical ensemble is uniform.
