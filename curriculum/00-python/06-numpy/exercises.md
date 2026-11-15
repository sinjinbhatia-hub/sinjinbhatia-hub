# PY-6 Exercises

Conventions as in the notes: `import numpy as np`, `rng = np.random.default_rng(seed)`. Put your answers in a file next to `project.py`; importing `project` re-runs its checks (about ten seconds), so copy the helpers you need instead. **Write your prediction down before you run anything**, and check every shape claim with an `assert` on `.shape`. "Write a function" means: write it, then test it with asserts you chose, including an edge case.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Strides.) Let `a = np.arange(24, dtype=np.float32).reshape(2, 3, 4)`. Predict `shape`, `strides` and whether the result is C-contiguous for `a`, `a.transpose(2, 0, 1)`, `a[:, ::2, 1]` and `a[..., ::-1]`. Then check with the library and with `np.shares_memory(a, ...)`.

**A2.** (Views.) Predict the printed lines, explain the bug, and give two different fixes.

```python
def zero_negatives(x):
    """Return x with negative entries replaced by 0; the caller's array must not change."""
    y = x
    y[y < 0] = 0
    return y


data = np.array([3, -1, 4, -1, 5])
out = zero_negatives(data)
print(data, out is data)
M = np.arange(6).reshape(2, 3)
row = M[1]
row *= 10
print(M)
```

**A3.** (Broadcasting.) For each pair of shapes under `+`, write the result shape or "error". For each error, write the shortest change that makes it work (a `[:, None]`, a `keepdims`, a reshape). Check with `np.broadcast_shapes`.
`(4,3)+(4,)`, `(4,3)+(4,1)`, `(5,1)+(1,6)`, `(2,3,4)+(3,1)`, `(2,3,4)+(3,)`, `(6,)+(6,1)`.

**A4.** (dtypes.) Predict each result, including any warning or exception (NumPy 2), then run:
`np.uint8(250) + np.uint8(10)`; `np.array([250], dtype=np.uint8) + 10`; `np.array([250], dtype=np.uint8) + 300`; `np.arange(5, dtype=np.int8) * 100`; the `.sum()` of that last array; `np.float32(16777216) + np.float32(1)`. For the first and fourth, write code that gives the mathematically correct answer.

**A5.** (Axes.) `X` has shape `(6, 4)`: 6 samples (rows), 4 features (columns). Without a loop, write the expression for each, with its shape: (a) the mean of each feature; (b) the 2-norm of each sample; (c) for each sample, the index of its largest feature; (d) `X` with every sample scaled to unit norm; (e) the feature with the largest variance (`ddof=1`); (f) the running total over samples. Verify against a loop on one case each.

**A6.** (`arange` versus `linspace`.) Predict `len(np.arange(0, 1, 0.1))`, `len(np.arange(0.1, 0.4, 0.1))` and the last element of the second. Explain the difference in one sentence using the length formula. Write `grid(a, b, n)` returning `n` equally spaced points including both ends, and `periodic_grid(a, b, n)` that excludes `b`; assert the spacings.

## B. Core (do all; these are the module)

**B1.** (Strides, Build 1 extended.) Write `strides_for(shape, itemsize, order="C")` for C and Fortran order, and `is_c_contiguous(shape, strides, itemsize)` that ignores axes of length 1 (as NumPy does). Test on at least six arrays: a fresh array in each order, a transpose, `a[::2]`, `a[:, 1]`, and one with a length-1 axis. Compare against `.strides` and `.flags["C_CONTIGUOUS"]`.

**B2.** (Broadcasting.) Write `zscore(X, ddof=1, axis=0)` returning `(Z, mean, std)` where `Z` has zero mean and unit standard deviation along `axis`, without loops. A column with zero variance must come out as zeros, not `nan`. It must work for `axis=0` and `axis=1` (hint: `keepdims`). Test: against a Python loop on a `(50, 4)` array, and with a constant column.

**B3.** (Pairwise distances, the Gram trick.) Squared distances satisfy $D_{ij}^2=\|\mathbf x_i\|^2+\|\mathbf x_j\|^2-2\,\mathbf x_i\cdot\mathbf x_j$. Implement `dist_gram(X)` with one `@`. (a) Check it against `cdist` on random data. (b) Add `1e8` to every coordinate and compare again: report the largest error and whether any squared distance is negative, and explain it in two sentences using cancellation (PY-1, LA-13). (c) Fix it by centering first. (d) For $n=2000$, $d=3$, time `dist_gram`, the broadcast version of Build 4, and `cdist`, and compute the byte size of each version's largest temporary. Write three sentences on which you would use when.

**B4.** (Softmax.) Write `softmax(Z, axis=-1)` that is stable for `Z = [1000, 1001, 1002]`. First show the naive `np.exp(Z) / np.exp(Z).sum()` fails (what does it return, and why?). Test that rows sum to 1, that it is invariant to adding a constant per row, and that it works on a `(B, T, V)` batch with `axis=-1`. Then write `log_softmax` using $\log\sum e^{z}=m+\log\sum e^{z-m}$ and compare with `scipy.special.log_softmax`.

**B5.** (Moving averages.) Implement `moving_average(x, w)` three ways: `sliding_window_view(...).mean(-1)`, `np.convolve`, and the cumulative-sum trick $\bar x_i=(S_{i+w}-S_i)/w$ with $S=\text{cumsum}$ (prepend a 0). Assert they agree on random data. Then repeat with `x = 1e6 + rng.standard_normal(10**6)` in `float32` and compare each with the `float64` answer: which loses digits, and why? Time all three for $N=10^6$, $w=1000$, and say why the windows version does about $w$ times more arithmetic.

**B6.** (Least squares by Vandermonde.) Build the Vandermonde matrix $V_{ij}=x_i^{\,j}$ for `x = np.linspace(0, 1, 40)` by broadcasting (`x[:, None] ** np.arange(deg + 1)`), and fit $y=\sin 2\pi x+0.05\,\eta$ ($\eta$ standard normal) with `np.linalg.lstsq` for degrees 3, 8, 15. Compare the coefficients with `np.polyfit` (note its ordering). For each degree report `np.linalg.cond(V)` and the error of the normal-equations solution `np.linalg.solve(V.T @ V, V.T @ y)` against `lstsq`. Explain what the table says about forming $V^\top V$ (LA-5, LA-13).

**B7.** (Extend `project.py`.) (a) Complete the four "Your turn" items (`assign_vec`, `cov_broadcast`, `life_step`, the fixed `row_normalize`) and make the project's cells accept them. (b) Add a Part 6b: a 2-D random walk with steps chosen uniformly from the four unit directions (draw indices with `rng.integers(0, 4, ...)` and look the steps up in a `(4, 2)` table). Positions have shape `(walkers, time, 2)`. For $5000\times500$, assert the MSD $\langle x^2+y^2\rangle$ has log-log slope $1\pm0.05$ and that it equals $t$ within five standard errors at the last step; save a figure of 10 paths.

## C. Challenge (do at least one)

**C1.** (Write `einsum`.) Write `my_einsum(spec, A, B)` for specs like `"ij,jk->ik"` or `"bij,bjk->bik"` with two operands, no repeated letter inside an operand, and an explicit output. Version 1: give each operand a length-1 axis for every letter it lacks (`None` indexing), multiply by broadcasting, `sum` over the letters missing from the output, `transpose` to the output order. Version 2: the same via `np.tensordot` for the summed letters. Test both on the five cases of §3.8 and on a random spec. For `ij,jk->ik` at $n=300$ report the size of Version 1's intermediate in MB and the times of Version 1, Version 2, `np.einsum` and `@`. Write a short conclusion: why does `@` win, and when is `einsum` the right choice anyway?

**C2.** (Where does vectorizing start to pay?) For $N$ from 1 to $10^6$ (log-spaced), measure the time to compute $\sum_i \sin(x_i)^2$ with (a) a Python loop over a list using `math.sin`, (b) `np.sin(x)**2` then `.sum()` on an array, and (c) `np.dot(np.sin(x), np.sin(x))`. Plot time per element against $N$ on log-log axes, find the crossover $N$ for (a) versus (b), and estimate NumPy's fixed per-call overhead from the small-$N$ plateau. Then answer, with your numbers: in a physics loop that handles one 3-vector per step for $10^6$ steps, would you store the vector as a `tuple`/`math` code or a NumPy array, and where do you vectorize instead (CM-1)?

**C3.** (N-body accelerations.) For $N$ bodies with masses $m_j$ and positions $\mathbf r_j$ (shape `(N, 3)`), the softened gravitational acceleration is $\mathbf a_i=\sum_{j\ne i}G m_j(\mathbf r_j-\mathbf r_i)/(\|\mathbf r_j-\mathbf r_i\|^2+\varepsilon^2)^{3/2}$ (the force law of CM-5, summed over all pairs). Implement it three ways: double loop; full broadcast with an `(N, N, 3)` array; and a chunked broadcast that processes `block` rows of `i` at a time. Assert agreement to $10^{-12}$ relative, and assert Newton's third law as a check: $\sum_i m_i\mathbf a_i\approx0$. Measure time and the byte size of the largest temporary for $N=200, 1000, 3000$ (skip the loop above 1000), and report the $N$ at which the full broadcast would exceed 1 GB. Choose `block` so the temporary stays under 50 MB and show the chunked time is within a small factor of the full broadcast.
