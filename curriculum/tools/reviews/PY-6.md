# Review: PY-6 NumPy — arrays, memory, broadcasting, vectorization, views vs copies

Reviewer pass started 2026-10-26. Module folder: `curriculum/00-python/06-numpy/`.
Brief: `curriculum/research/briefs/PY-6.md` (has REVIEW FOCUS). Scratch: /tmp/claude-0/py6/ (harness.py, check.py).

Baseline: `MPLBACKEND=Agg python3 project.py` exits 0 in 8.1 s (Python 3.13.16, numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2).

Method: (1) every fenced `python` block in notes.md (35) and solutions.md (21) extracted and run cumulatively in one
namespace by /tmp/claude-0/py6/harness.py, stdout diffed against the stated `text` block: every non-timing output
reproduces exactly (strides, flags, dtypes, overflow values, broadcast shapes and strides, the axis drill, sort/search
outputs, einsum shapes, residual/error table, float32 summation errors, MSD 992.8 / slope 0.9965, live-pixel counts
213584/137632/67222/51626, convolution shapes, A1-A6, B1-B7, C1-C3 asserts). The only diffs are timings (within 2x).
(2) every hand-derived number recomputed independently in /tmp/claude-0/py6/check.py (sympy/numpy): (3.2) strides on
int64 (2,3), float64 (3,4) and its transpose, float32 (2,3,4), view start bytes 8 and 24, Example 4.1 offset 88; int8/uint8
wraps (-56, 4, [0,100,-56,44,-112], sum -24, 44, 255), 2^24+1 in float32, float32 spacing 0.0625 near 1e6 and 65536 near
1e12, float64 spacing 4 near 3e16, 1/eps = 4.5e15; the arange length rule on all four cases; all 15 broadcast shape pairs
(notes, Example 4.3, §9, A3) and 14 + 300 = 314; (3.7) and Var(x_t^2) = 2t^2 - 2t (sympy), relative SE sqrt(2/1e4) = 1.4%;
B7(b)'s Var(|r_t|^2) = t^2 - t for the 2-D lattice walk (sympy and a 2e5-walker Monte Carlo: 375.8 vs 380 at t = 20) and
5 SE = 0.071; the chain eigenvalues 2 +- sqrt2, 2 and their roots 1.8478/1.4142/0.7654; the inertia tensor diag(8,1,9);
A B = [[2,1],[4,3]], trace 5; the four-point Mandelbrot trace and |1+3i| = sqrt10; softmax of [1000,1001,1002]; isclose
cases; std 1.1180/1.2910; memory sizes 2.2 MB / 2.4 GB / 32 MB / 96 MB / 216 MB / N = 6455 / block 694 = 50 MB; argsort
and inverse permutation; np.partition semantics; the convolve2d signature (no `method`); NEP 50 promotion cases.
(3) cross-references checked against the modules on disk: LA-2/5/6/7/8/9/12/13, CM-1/5/6/7, F3/F5, ST-1, PY-1/2/4/5/7-10,
ML-8/9 briefs, ROADMAP §7 ("90% tensor-shape bookkeeping" exists). (4) REVIEW FOCUS: no `np.random.seed`/`rand`/
`np.matrix` outside the one warning sentence; `inv` only in the why-not demonstration; loops appear only as the measured
baselines, Example 4.6's iteration (justified by §3.10 case i) and solution reference checks; imports numpy/scipy/
matplotlib only (+ stdlib); solutions run; project 8.1 s headless.
Template: §0-§10 with the template headings, conventions box, 8 Build sequences, 6 worked examples (4.4 from CM-6, 4.5
from CM-7), 6 A / 7 B / 3 C, B7 extends project.py, §8 specific and matches the brief's list, §10 matches the brief.
Word count 9,441 (under the PY ceiling of 9,500; above the 6,000-8,000 target, not changed).
History: Higham ch. 4 (summation bounds), Harris et al. Nature 585 (2020), Oliphant *Guide to NumPy*, NEP 50 weak
scalars, Python 3.12 compensated `sum` -- all correct; nothing dubious.

Format: **ID** | file:location | severity | what is wrong | fix | status

## Findings

**1** | solutions.md C1, the timing table after "Cost, for `ij,jk->ik` at n=300" and the conclusion paragraph | medium |
The stated table has `np.einsum, optimize=True` (4.90 ms) beating `P @ R` (7.94 ms) and plain `np.einsum` only 1.5x
slower than `@`, which contradicts the paragraph's conclusion ("`@` wins because it reaches BLAS"; "plain `np.einsum`
loops in C but without BLAS"). A 300x300 `@` at 7.94 ms is itself anomalous (BLAS does it in ~0.4 ms; reproduced 0.38 ms);
the table came from a run with a throttled BLAS. | paste a real run that supports the conclusion (@ 0.38, optimize 0.47,
tensordot 0.63, einsum 6.24, broadcast 94 ms) and quote those numbers in the prose | FIXED

**2** | solutions.md C2, conclusion paragraph after the 3-vector block | medium | The prose quotes numbers that are not in
its own output blocks: "1.6 µs per call", "0.16 µs per element", "crossover is N = 31", "83 ns against 11 ns", "966 ns
against 190 ns", while the stated output says 3.9 µs, 0.15 µs, N = 56, 85/11.0 ns, 892/191 ns. A solution must quote its
own measurements. | paste the real rerun output and write the conclusion from those numbers | FIXED

**3** | solutions.md B5, explanation after the table | medium | "`ma_windows` is off by about 0.1 ... which is as good as
the stored data allows" is wrong: the reference is the float64 mean of the same stored float32 values, and
`sliding_window_view(big, 1000).mean(axis=-1, dtype=np.float64)` reproduces it to 0.0 (verified). The 0.107 comes from
accumulating sums near 1e9 in float32, where the spacing is 64, i.e. 0.064 per rounding after dividing by w = 1000. The
prose also quotes "3.9 ms against 233 ms" while the table says 4.2 and 245.6. | correct the explanation; quote the
table's numbers | FIXED

**4** | notes.md §1 ("the random-walk ensemble from F5"), §3.10 ("The physics check is F5's"), §5 ("F5: Build 7 is the
random-walk MSD") | medium | F5 (00-foundations/05-probability-2) has no random-walk MSD: its random walks are the walk
on a graph (§3.12) and a 2-D return-probability task in "Your turn". The object PY-6 uses is F5 (3.14),
Var(sum a_i X_i) = a^T Sigma a, which for independent +-1 steps gives <x_t^2> = t. | re-target the three links to
F5 (3.14) (and F5's "Your turn" 2-D walk for B7(b)) | FIXED

**5** | notes.md §0 ("LA-13 builds on §3.9"), §3.8 ("LU with partial pivoting is backward stable, LA-13"), §5 ("LA-13:
§3.8's solve versus inv and §3.9's summation orders are the first two experiments of its stability study") | medium |
LA-13 has no solve-versus-inv experiment and no summation-order experiment (notes and project.py grepped). What it has:
§2 "never form A^{-1}", Theorem 3.8 (forward error <~ kappa x backward error), Theorem 3.9 (Wilkinson: pivoted LU is
backward stable *in practice*, growth-factor caveat), §3.1 rounding model. "Backward stable" without the caveat also
overstates LA-13 Theorem 3.9. | name those objects instead | FIXED

**6** | notes.md §5 ("F3: DFT sample points are `linspace(..., endpoint=False)`") | low | F3 builds its sample grid as
`np.arange(N)` times the step (t_n = n/f_s in its notes and project.py); it never uses `linspace(..., endpoint=False)`.
| say what F3 does and that `endpoint=False` is the `linspace` equivalent | FIXED

**7** | notes.md §3.2, paragraph after the promotion block ("A Python scalar never promotes the array: `float32 + 1.0`
stays `float32`") | low | Contradicted by the block's own output: `(np.arange(3) + 0.5).dtype` is `float64`. Under
NEP 50 a Python scalar never widens the array's precision within its kind, but a Python float still turns an integer
array into float64 (verified: int8 + 1.5 -> float64). | state the rule precisely | FIXED

**8** | notes.md §3.7, bullet on `np.partition` ("puts the k smallest values first") | low | `np.partition(v, k)`
places the element belonging at sorted position k there with everything <= before it, so `[:k+1]` are the k+1
smallest; the block demonstrates exactly that with `np.partition(v, 1)[:2]`, so the sentence is off by one against
its own demonstration. | reword | FIXED

**9** | notes.md §3.11, last sentence of Build 8 ("For large kernels SciPy switches to FFT methods, F3") | low |
`scipy.signal.convolve2d` has no `method` argument and is always direct (signature checked); `scipy.signal.convolve`
with `method="auto"` is what switches to FFT. | name the right function | FIXED

**10** | notes.md §7 ("In order it builds and asserts: (1) ... (8)") | low | Not the project's order: Part 1 strides/
dtypes, 2 views, 3 matmul, 4 broadcasting + axis drill + distances, 5 sorting/linalg/einsum/floating point, 6 random
walk, 7 Mandelbrot, 8 windows/convolution, 9 performance/I-O; the notes list distances as item 4 alone and einsum last.
| rewrite the list to match the parts | FIXED

**11** | project.py header ("about 15 s") vs exercises.md intro ("about eight seconds") | low | The script runs in 8.1 s
here; the two statements disagree. | "about 10 s" in both | FIXED

**12** | project.py Part 5, the `make_matrix` loop | medium | The notes' block prints both residuals and *forward
errors*, and the prose claims "its error is only what the conditioning allows, about kappa eps_mach" and that `inv`
loses "more digits"; the project asserts only residuals, so the forward-error claims are unchecked. | compute the
forward errors and assert err_solve < 10 kappa eps for each kappa and err_inv > 10 err_solve at kappa = 1e12 (measured
3e-4 vs 6e-6; seeded, deterministic) | FIXED

**13** | project.py Part 3 (`t_loops[-1] > 30 * t_at[-1]`) and Part 6 (`t_fast < 1.5 * t_big`) | low | Both are far
weaker than the notes' claims (loops/@ ~ 3,900 at n = 200; the narrow-dtype walk is *faster*, 0.064 s vs 0.116 s). A
30x floor passes a `@` 100x slower than it is; "< 1.5x" passes a "cheaper" version that is slower. | 300x and
`t_fast < t_big` (reproduced with >= 2x margin here) | FIXED

**14** | solutions.md B7(b), "For this lattice walk Var(|r_t|^2) = t^2 - t" | low | Stated without derivation and the
assert depends on it. (Verified: |r_t|^2 = t + 2 sum_{s<r} xi_s.xi_r, each product +-1 w.p. 1/4 each and 0 w.p. 1/2,
variance 1/2, the C(t,2) products uncorrelated, so Var = 4 C(t,2) / 2 = t(t-1); Monte Carlo 375.8 vs 380 at t = 20.)
| add the three-line derivation | FIXED

**15** | notes.md Builds 4, 6, 7, 8 | low | Each has a ▶ Try it first but the explanation after its block carries no
■ What just happened marker (9 ▶ against 5 ■); the PY addendum requires one after every ▶. | add the marker to the
existing explanatory paragraph | FIXED

**16** | notes.md cross-links that name a module but not the object: §0 "PY-2 the list and its copying slices", "PY-4
the data model"; Example 4.4 "(CM-6)"; §5 "CM-1, CM-5", "CM-6", "CM-7"; §3.6 "ST-1 derives why" | low | STYLE requires
module code + object. The objects exist (PY-2 §3.3: `a[:]` is a shallow copy; PY-4 §3.4: `__getitem__` receives the
tuple; CM-6 Example 4.5's N-mass chain; CM-1 §3.2's N-particle system and CM-5 §3.13 three bodies; CM-7 §3.3 and LA-12
Example 4.6's einsum; ST-1's unbiasedness of S^2) but are not named. | name them | FIXED

**17** | notes.md §3.2, the stated output of the int8 block ("RuntimeWarning: overflow encountered in scalar add") |
low | The warning goes to stderr with a file:line prefix, so this line is an edited rendering inside a block the
conventions box calls real stdout. | say in the ■ paragraph that the warning line is stderr, shown in place | FIXED

Not changed: timings in every stated block differ from reruns by < 2x (allowed by the conventions box); word count.

## Verification after fixes

- `MPLBACKEND=Agg python3 project.py`: exit 0, 6.2 s and 6.4 s on two runs; the new forward-error asserts pass
  (kappa = 1e12: error solve 6.3e-6, inv 3.2e-4; residual solve 5.0e-16, inv 1.1e-5); the tightened timing floors
  (loops/@ > 300 at n = 200, measured 3,100-3,600; narrow-dtype walk faster than `choice`, 0.04-0.05 s vs 0.10-0.64 s) held.
- harness.py rerun on notes.md (35 blocks) and solutions.md (21 blocks): the same timing-only diffs as before the
  edits; every non-timing output still reproduces exactly, including the three replaced solution tables (B5, C1, C2),
  whose non-timing columns (errors 0.107 / 4.66e-10 / 1.7e4, 216 MB, crossover N = 31) reproduced on a second run.
- No `\"` or other artefacts left by the patch scripts; 9 ▶ and 9 ■ in notes.md.

Status: 17 findings (6 medium, 11 low), 17 fixed, 0 open; project.py exits 0 in ~6.5 s.
