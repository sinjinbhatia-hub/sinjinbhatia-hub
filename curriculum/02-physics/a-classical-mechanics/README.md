# 02a · Classical Mechanics

Eight modules, about 80 hours. Newton to Hamilton, with the explicit purpose of making quantum mechanics (QM-4), statistical mechanics (SM-3), and the structure of modern physics intelligible. For an ME student, CM-3, CM-6, and CM-7 are also the theory behind dynamics, vibration, and rotating machinery courses.

| # | Module | Hours | Needs | The one idea |
|---|---|---|---|---|
| 1 | [Newtonian mechanics revisited](01-newtonian-mechanics-revisited/notes.md) | 6 | ODEs | Constraints make Newton painful; feel the pain |
| 2 | [Least action](02-least-action/notes.md) | 8 | F6, CM-1 | $\delta\int(T - V)\,dt = 0$ is equivalent to Newton and coordinate-free |
| 3 | [Lagrangian mechanics](03-lagrangian-mechanics/notes.md) | 12 | CM-2 | Choose your coordinates; constraint forces vanish or return as multipliers |
| 4 | [Symmetry and conservation](04-symmetry-and-conservation/notes.md) | 8 | CM-3 | Every continuous symmetry is a conservation law (Noether) |
| 5 | [Central forces](05-central-forces/notes.md) | 10 | CM-4 | Two bodies are one body in an effective potential; Kepler falls out |
| 6 | [Small oscillations and normal modes](06-small-oscillations-normal-modes/notes.md) | 12 | CM-3, LA-8 | $K\mathbf{a} = \omega^2 M\mathbf{a}$: the spectral theorem with a mass-weighted inner product |
| 7 | [Rigid bodies](07-rigid-bodies/notes.md) | 12 | CM-3, LA-8, LA-12 | $\mathbf{L} = I\boldsymbol\omega$ with $I$ a symmetric tensor; principal axes exist because of LA-8 |
| 8 | [Hamiltonian mechanics](08-hamiltonian-mechanics/notes.md) | 12 | CM-4 | Phase space, $\{q, p\} = 1$, incompressible flow; the classical shadow of QM |

## Bridges out of this course

- **CM-8 → QM-4:** $\{A, B\} \to [\hat A, \hat B]/i\hbar$. Heisenberg's equation is $dF/dt = \{F, H\}$ with hats.
- **CM-8 → SM-3:** Liouville's theorem is why the microcanonical ensemble is uniform on the energy surface.
- **CM-6 → QM-3, SM-5:** every normal mode is a harmonic oscillator; quantize each one and you have phonons and the Debye model.
- **CM-5 → QM-6:** the hydrogen atom's radial equation is the Kepler problem's effective potential with $\hbar$.
- **CM-4 → QM-5:** rotational symmetry gives angular momentum conservation classically and the $2\ell + 1$ degeneracy quantum mechanically.
- **CM-3 → ML-2:** Lagrange multipliers for constraints are the same multipliers in constrained optimization and the KKT conditions.
