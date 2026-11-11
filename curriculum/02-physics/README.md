# 02 · Physics

Four sub-courses, 32 modules, about 360 hours. They are listed in the order you'll start them, but they interleave with the other courses; see [../ROADMAP.md §6](../ROADMAP.md#6-recommended-study-path).

| Sub-course | Modules | Hours | Needs first | Feeds |
|---|---|---|---|---|
| [a · Classical mechanics](a-classical-mechanics/) | 8 | ~80 | F6, LA-8, LA-12 (for CM-6/7) | QM-4 (Poisson brackets → commutators), SM-3 (Liouville), QM-6 (central forces) |
| [b · Electromagnetism](b-electromagnetism/) | 8 | ~90 | multivariable calculus, F2, F3, LA-12 (for EM-8) | QM-2 (waves), QC hardware context |
| [c · Thermodynamics and statistical mechanics](c-thermo-statistical-mechanics/) | 7 | ~80 | F4, F5, QM-3 (for SM-5), CM-8 (for SM-3) | ML-10, ST-6 (Metropolis), SM-7 bridge to ML |
| [d · Quantum mechanics](d-quantum-mechanics/) | 9 | ~110 | LA-8, LA-12, F3, CM-8 | QC-1 onward, SM-5 |

## How the four connect to each other

```
CM-8 Hamiltonian ──Poisson brackets──► QM-4 commutators
   │                                        │
   └──Liouville──► SM-3 ensembles           └──► QM-5 spin ──► QC-1 qubits
                     │
EM-6 Maxwell ──► EM-7 waves ──► QM-1 photons, QM-2 wave packets
                     │
                     └──► SM-5 blackbody (the problem that started QM)
```

## Conventions shared across the four

SI units throughout (Griffiths' choice for E&M). $\hbar$ explicit in QM. $k_B$ explicit in stat mech. Lagrangian $L = T - V$. Complex inner products conjugate-linear in the first slot, matching bra-ket. Fourier transforms in physics modules use the symmetric $1/\sqrt{2\pi}$ convention and say so at the top.
