# 05 · Quantum Computing

Nine modules, about 94 hours. Mermin's "qubits first" order with Nielsen & Chuang as the reference, built entirely in NumPy: you write the state-vector simulator, the density-matrix simulator, the stabilizer tableau, the surface-code decoder, and the variational loop yourself.

**Two entry routes.** Physics route: after QM-7. Math route: after LA-8 and LA-12, using QC-1's self-contained postulates (you lose physical intuition, not correctness).

**Conventions used everywhere:** $|0\rangle = (1,0)^\top$; qubit 0 is the **leftmost** tensor factor (big-endian; Qiskit is the opposite); circuits read left to right but matrices multiply right to left; the QFT carries $e^{+2\pi i jk/N}$ (the inverse of NumPy's FFT up to normalization).

| # | Module | Hours | Needs | The one idea |
|---|---|---|---|---|
| 1 | [Qubits and single-qubit gates](01-qubits-and-single-qubit-gates/notes.md) | 8 | LA-8 (QM-5 ideal) | A qubit is a direction on a sphere; gates are rotations; one bit comes out |
| 2 | [Many qubits and universal gates](02-many-qubits-and-universal-gates/notes.md) | 10 | QC-1, LA-12 | Tensor products; CNOT makes entanglement; phase kickback; $\{H, T, \mathrm{CNOT}\}$ is universal; no-cloning |
| 3 | [Protocols](03-protocols/notes.md) | 8 | QC-2, QM-7 | Superdense coding and teleportation trade ebits against bits; BB84; the CHSH game (75% vs 85%) |
| 4 | [Early algorithms](04-early-algorithms/notes.md) | 10 | QC-2 | The Hadamard transform is Fourier analysis on $\mathbb Z_2^n$; Simon's exponential separation |
| 5 | [QFT, phase estimation, Shor](05-qft-phase-estimation-shor/notes.md) | 14 | QC-4, F3 | The QFT is the DFT on amplitudes; phase estimation; order finding factors integers |
| 6 | [Grover and amplitude amplification](06-grover-amplitude-amplification/notes.md) | 8 | QC-4, LA-5 | Two reflections make a rotation; $\sqrt N$ is optimal |
| 7 | [Noise and channels](07-noise-and-channels/notes.md) | 10 | QC-2, QM-7 | CPTP maps, Kraus operators, $T_1$ and $T_2$, fidelity and trace distance |
| 8 | [Quantum error correction](08-quantum-error-correction/notes.md) | 16 | QC-7 | Syndrome measurement discretizes errors; stabilizers are GF(2) linear algebra; the surface code and its threshold |
| 9 | [Simulation, VQE, QAOA](09-simulation-vqe-qaoa/notes.md) | 10 | QC-5, QM-8, CM-8 | Trotter is a symplectic splitting; VQE is the variational principle; the honest NISQ landscape |

## Where the threads meet

- **LA-8 → QM-4 → QC-1:** Hermitian observables, unitary evolution, projective measurement.
- **LA-12 → QM-7 → QC-2, QC-8:** tensor products, Schmidt rank, stabilizer codes.
- **F3 → QC-4, QC-5:** the Hadamard transform and the QFT are Fourier transforms on groups.
- **CM-8 → QC-9:** symplectic integrators and Trotter formulas are one idea.
- **SM-6 → QC-8, QC-9:** the surface-code threshold is an Ising-type phase transition; QAOA cost functions are Ising Hamiltonians.
- **QM-8 → QC-1, QC-9:** gates are Rabi pulses; VQE is the variational principle; adiabatic computing is the adiabatic theorem.
