# 02d · Quantum Mechanics

Nine modules, about 110 hours. Wave mechanics first (because you know PDEs), then a heavy Dirac-notation block (because quantum computing lives there), then spin, hydrogen, entanglement, approximation methods, and the loose ends: symmetry, the classical limit, path integrals, decoherence.

| # | Module | Hours | Needs | The one idea |
|---|---|---|---|---|
| 1 | [Why quantum](01-why-quantum/notes.md) | 6 | LA-8 | Six experiments that forced it; the postulates are the spectral theorem described by nature |
| 2 | [Wavefunction and Schrödinger equation](02-wavefunction-schrodinger/notes.md) | 12 | QM-1, F3 | $|\Psi|^2$ is a probability; $\hat p = -i\hbar\partial_x$ is derived; Fourier uncertainty is Heisenberg uncertainty |
| 3 | [One-dimensional problems](03-one-dimensional-problems/notes.md) | 16 | QM-2 | Confinement quantizes; ladder operators; tunneling is exponential |
| 4 | [Formalism and Dirac notation](04-formalism-dirac/notes.md) | 16 | QM-3, LA-8, LA-12 | A wavefunction is a coefficient list; measurement is projection; $\{A,B\}\to[\hat A,\hat B]/i\hbar$ |
| 5 | [Angular momentum and spin](05-angular-momentum-and-spin/notes.md) | 14 | QM-4 | The commutator algebra alone fixes the spectrum; spin-1/2 is a qubit |
| 6 | [Hydrogen atom](06-hydrogen-atom/notes.md) | 12 | QM-5, CM-5, EM-3 | Kepler's effective potential with $\hbar$; series termination gives Balmer |
| 7 | [Composite systems and entanglement](07-composite-systems-entanglement/notes.md) | 12 | QM-5, LA-12 | Tensor products; Schmidt rank; no signaling; Bell's inequality; density matrices; Pauli exclusion |
| 8 | [Approximation methods](08-approximation-methods/notes.md) | 14 | QM-6 | Perturbation theory, variational principle, WKB, Rabi oscillations, the golden rule |
| 9 | [Symmetry, classical limit, path integrals, decoherence](09-symmetry-classical-limit-path-integrals/notes.md) | 8 | QM-8, CM-8 | Noether with hats; $e^{iS/\hbar}$; imaginary time is statistical mechanics; the environment measures |

## Bridges out of this course

- **QM-5, QM-7 → QC-1, QC-2:** a qubit is a spin-1/2; two qubits live in $\mathbb C^2\otimes\mathbb C^2$; entanglement is Schmidt rank $>1$.
- **QM-8 → QC-1, QC-9:** a single-qubit gate is a Rabi pulse; VQE is the variational principle; adiabatic computing is the adiabatic theorem.
- **QM-7, QM-9 → QC-7, QC-8:** density matrices and decoherence are why error correction exists.
- **QM-3 → SM-5:** box states and oscillators are what Bose–Einstein and Fermi–Dirac statistics count.
- **QM-9 → SM-4:** imaginary-time quantum mechanics is the partition function.
