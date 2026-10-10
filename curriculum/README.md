# Curriculum: Linear Algebra → Physics → Statistics → Machine Learning → Quantum Computing

A self-contained, derivation-first course of study built for one student: a second-year ME + Applied Math undergraduate who needs to rebuild linear algebra from the ground up and take it through to quantum computing and transformers.

**Start here:** [ROADMAP.md](ROADMAP.md) explains every design decision, the dependency graph, time estimates, and where you'll struggle. Then come back to this page and follow the study path.

## Courses

| Folder | Course | Modules | Hours | Status |
|---|---|---|---|---|
| [00-foundations](00-foundations/) | Proof and complex numbers · Fourier series and transforms · Probability I/II · Calculus of variations | 6 | ~62 | building |
| [01-linear-algebra](01-linear-algebra/) | From vectors to SVD, then Jordan form, tensors, numerical LA, matrix calculus | 14 | ~156 | building |
| [02-physics/a-classical-mechanics](02-physics/a-classical-mechanics/) | Newton → least action → Lagrange → Noether → orbits → normal modes → rigid bodies → Hamilton | 8 | ~80 | building |
| [02-physics/b-electromagnetism](02-physics/b-electromagnetism/) | Vector calculus → electrostatics → BVPs → matter → magnetostatics → Maxwell → waves → radiation and relativity | 8 | ~90 | building |
| [02-physics/c-thermo-statistical-mechanics](02-physics/c-thermo-statistical-mechanics/) | Thermodynamics → ensembles → quantum statistics → phase transitions → the ML bridge | 7 | ~80 | building |
| [02-physics/d-quantum-mechanics](02-physics/d-quantum-mechanics/) | Wavefunctions → 1D problems → Dirac formalism → spin → hydrogen → entanglement → approximations | 9 | ~110 | building |
| [03-statistical-modeling](03-statistical-modeling/) | Inference → testing → regression → GLMs → Bayes → MCMC → time series → causal inference | 9 | ~108 | building |
| [04-machine-learning](04-machine-learning/) | Learning problem → optimization → linear → kernels → trees → unsupervised → neural nets → architectures → transformers → generative | 10 | ~130 | building |
| [05-quantum-computing](05-quantum-computing/) | Qubits → gates → protocols → algorithms → QFT/Shor → Grover → noise → error correction → simulation | 9 | ~94 | building |

Every module has the same four files: `notes.md`, `exercises.md`, `solutions.md`, `project.py`. See [TEMPLATE.md](TEMPLATE.md). The structure is explained in [STYLE.md](STYLE.md). The research behind the ordering is in [research/SOURCES.md](research/SOURCES.md).

## The study path, condensed

Work two tracks at once, never more. Full version with week numbers in [ROADMAP.md §6](ROADMAP.md#6-recommended-study-path).

1. **Phase 0 (2 weeks):** F1 → LA-1 → LA-2
2. **Phase 1 (12 weeks):** LA-3 → LA-9 (track A) ‖ F4 (track B). Checkpoint: `01-linear-algebra/CHECKPOINT-core.md`.
3. **Phase 2 (12 weeks):** F2 → F3 → LA-10 → LA-14 ‖ F5 → F6 → CM-1 → CM-4
4. **Phase 3 (19 weeks):** CM-5 → CM-8 → EM-1 → EM-4 ‖ ST-1 → ST-4 → ML-1
5. **Phase 4 (17 weeks):** EM-5 → EM-8 → QM-1 → QM-5 ‖ ST-5 → ST-6 → SM-1 → SM-4 → ST-7 → ST-8
6. **Phase 5 (16 weeks):** QM-6 → QM-9 → QC-1 → QC-9 ‖ ML-2 → ML-10 → SM-5 → SM-7 → ST-9

Shortcuts for ML-first or QC-first are in the roadmap.

## Running the code

```bash
pip install numpy scipy matplotlib sympy
cd curriculum/01-linear-algebra/01-vectors-and-two-pictures
MPLBACKEND=Agg python3 project.py        # figures land in ./figures/
```

Every `project.py` uses `# %%` cell markers, so it opens as a notebook in VS Code or JupyterLab (with Jupytext). Generated `.ipynb` versions sit next to each script once the course is built; `tools/py_to_ipynb.py` regenerates them.

## How to use a module

1. Read §0–2 of `notes.md` (where it sits, the problem, the intuition). Stop. Try to guess the main result.
2. Work §3 with pen and paper. At every **▶ Try it first**, actually try it first.
3. Do the worked examples yourself before reading them.
4. Do all A and B exercises. Check against `solutions.md` only after a real attempt.
5. Run `project.py`, then do the "Your turn" cells.
6. Answer §9 closed-book. If you can't, don't move on.
