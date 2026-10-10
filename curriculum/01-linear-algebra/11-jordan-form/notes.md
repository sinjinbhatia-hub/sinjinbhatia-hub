# LA-11: Jordan Form, Generalized Eigenvectors, the Minimal Polynomial, Cayley–Hamilton, and Functions of Matrices

**Time:** ~10 hours · **Prerequisites:** LA-7 (LA-10 is used in two proofs; see §0) · **Unlocks:** LA-12 (indirectly), CM-6 (the critically damped oscillator), QM-8 (degenerate perturbation theory intuition)

> **Conventions in this module.** Scalars are complex unless stated otherwise; $A$ is an $n\times n$ matrix, equivalently an operator on $V = \mathbb{C}^n$. $I$ is the identity, $\ker$ the null space, $\operatorname{im}$ the column space (range). The characteristic polynomial is the **monic** one, $p_A(x) = \det(xI - A) = (-1)^n\det(A - xI)$; its roots are the eigenvalues. $N_s$ is the $s\times s$ **shift**: ones on the superdiagonal, zeros elsewhere, so $N_s\mathbf{e}_1 = \mathbf{0}$ and $N_s\mathbf{e}_j = \mathbf{e}_{j-1}$. A **Jordan block** is $J_s(\lambda) = \lambda I_s + N_s$: $\lambda$ on the diagonal, ones on the **superdiagonal** (Axler, Horn & Johnson, Strang and SymPy all do this; some older books put the ones below the diagonal, which is the same thing with the basis reversed). A **chain** $\mathbf{v}_s \to \mathbf{v}_{s-1} \to \cdots \to \mathbf{v}_1 \to \mathbf{0}$ means $N\mathbf{v}_j = \mathbf{v}_{j-1}$ and $N\mathbf{v}_1 = \mathbf{0}$; as a basis it is always written **bottom first**, $(\mathbf{v}_1, \dots, \mathbf{v}_s)$, which is what puts the ones above the diagonal. For an eigenvalue $\lambda$: $d_k(\lambda) = \dim\ker(A - \lambda I)^k$, $a_\lambda$ is the algebraic multiplicity (multiplicity of $\lambda$ as a root of $p_A$), $g_\lambda = d_1(\lambda)$ the geometric multiplicity. Where an inner product appears (Proposition 3.12) it is $\langle\mathbf{u},\mathbf{v}\rangle = \mathbf{u}^*\mathbf{v}$, conjugate-linear in the first slot, the physics convention of LA-8.

## 0. Where this sits

LA-7 ended with a gap. Diagonalization $A = PDP^{-1}$ needs $n$ independent eigenvectors, and you met matrices that don't have them: $\begin{pmatrix}1&1\\0&1\end{pmatrix}$ has eigenvalue $1$ twice and a single eigenvector. LA-7 called such matrices *defective* and moved on. This module closes the gap: every square matrix over $\mathbb{C}$ is similar to a matrix that is diagonal *plus a nilpotent shift*, the Jordan form, and that shift is exactly what produces the $t\,e^{\lambda t}$ solutions you have met in ODEs since your first course.

You need LA-7 cold: eigenvalues as roots of the characteristic polynomial, similarity $P^{-1}AP$, $e^{At}$ for diagonalizable $A$, and the fact that every complex matrix has at least one eigenvalue. You need rank–nullity (LA-4) and direct sums (LA-3). Schur's triangularization (LA-10) is used in two proofs; if you haven't done LA-10, all those proofs need is that every complex matrix is similar to an upper triangular one, and Exercise C1 proves that from scratch in half a page.

What it makes possible: CM-6 (the critically damped oscillator *is* a Jordan block; normal modes with repeated frequencies are *not*, and Proposition 3.12 says why), QM-8 (what a repeated eigenvalue does to eigenvectors under perturbation; the $\varepsilon$-demo of §3.11 is the intuition), ML-2 (heavy-ball momentum at its optimal setting is a Jordan block), ST-7 (companion matrices of AR processes are nonderogatory, Exercise B5), and LA-13 (the "never compute it in floating point" lesson, with the tools to say why).

## 1. The problem that forced this idea

Solve $y'' - 2y' + y = 0$. The characteristic equation $r^2 - 2r + 1 = (r-1)^2 = 0$ has one root, $r = 1$, so $e^t$ is a solution. A second-order equation needs two independent solutions. The other one is $t\,e^t$. You verified this in your ODE course by substitution and were told "for a repeated root, multiply by $t$." Euler had this rule in the 1740s. But *why*? Where does a polynomial factor $t$ come from, in a problem whose only ingredients are exponentials?

Turn it into a system, as LA-7 taught. With $x_1 = y$, $x_2 = y'$:
$$\mathbf{x}' = A\mathbf{x}, \qquad A = \begin{pmatrix} 0 & 1 \\ -1 & 2\end{pmatrix}.$$
LA-7's recipe: find eigenvalues, find eigenvectors, write $\mathbf{x}(t) = \sum_i c_i e^{\lambda_i t}\mathbf{v}_i$. The characteristic polynomial is $p_A(x) = x(x - 2) + 1 = (x-1)^2$, so $\lambda = 1$ twice. Eigenvectors: $A - I = \begin{pmatrix}-1&1\\-1&1\end{pmatrix}$ has rank $1$, so its kernel is a *line*, spanned by $(1,1)^\top$. One eigenvector. The recipe produces $c_1 e^t (1,1)^\top$ and stops. Half the solutions are missing. Diagonalization fails: there is no basis of eigenvectors, so there is no $P$ with $P^{-1}AP$ diagonal.

This is not a rare accident you can ignore. It happens in every critically damped oscillator (CM-6), in companion matrices of ODEs and recurrences with repeated roots (ST-7), in transition matrices with chains of transient states (F5), in momentum methods tuned to their fastest setting (ML-2), and in every matrix close enough to one of those that floating point can't tell the difference (§3.11).

The question "what is the simplest matrix similar to a given $A$, when diagonal is impossible?" was answered by Camille Jordan in his 1870 *Traité des substitutions et des équations algébriques*. Jordan was studying linear substitutions over finite fields for group theory, not ODEs, and the canonical form appears there as a tool. Weierstrass had obtained an equivalent result for complex matrices two years earlier, in his theory of elementary divisors; the name stuck to Jordan. Twelve years before either, Arthur Cayley's 1858 *Memoir on the Theory of Matrices* stated that every matrix satisfies its own characteristic equation, verified it for $2\times 2$ and $3\times 3$, and wrote that he had "not thought it necessary to undertake the labour of a formal proof of the theorem in the general case." Hamilton had the $3\times 3$ case in 1853 from quaternions. Frobenius gave the first general proof in 1878 and in the same paper introduced the minimal polynomial. Sylvester named nilpotent matrices in 1884.

The answer, in one sentence: **after a change of basis, every matrix is diagonal plus a shift, and the shift is what makes $t\,e^{\lambda t}$ appear.** The rest of this module makes each word precise, proves it, shows you how to compute with it, and tells you when you must not.

## 2. Intuition

### Stretch plus shift

Look at the $2\times 2$ Jordan block
$$J = \begin{pmatrix}\lambda & 1\\ 0 & \lambda\end{pmatrix} = \lambda I + N, \qquad N = \begin{pmatrix}0&1\\0&0\end{pmatrix}.$$
It does two things. The $\lambda I$ part stretches every vector by $\lambda$. The $N$ part is a **shift**: $N\mathbf{e}_2 = \mathbf{e}_1$ and $N\mathbf{e}_1 = \mathbf{0}$. It moves $\mathbf{e}_2$ onto $\mathbf{e}_1$ and kills $\mathbf{e}_1$. Apply it twice and everything is dead: $N^2 = 0$. An operator some power of which is zero is called **nilpotent**.

```
   e_2  --N-->  e_1  --N-->  0
```

So $J\mathbf{e}_1 = \lambda\mathbf{e}_1$ ($\mathbf{e}_1$ is an honest eigenvector) but $J\mathbf{e}_2 = \lambda\mathbf{e}_2 + \mathbf{e}_1$: $\mathbf{e}_2$ is stretched *and nudged* toward $\mathbf{e}_1$. It is "almost an eigenvector": $(J - \lambda I)\mathbf{e}_2 \ne \mathbf{0}$ but $(J-\lambda I)^2\mathbf{e}_2 = \mathbf{0}$. That is the whole content of the phrase *generalized eigenvector*.

A $3\times 3$ block is a longer chain,
```
   e_3  --N-->  e_2  --N-->  e_1  --N-->  0        N = J - λI,  N^3 = 0
```
and the general picture of a nilpotent operator is a **collection of chains** of various lengths, each ending in $\mathbf{0}$:
```
   v_3 --> v_2 --> v_1 --> 0          (a chain of length 3)
           w_2 --> w_1 --> 0          (a chain of length 2)
                   u_1 --> 0          (a chain of length 1: an eigenvector)
```
The bottom of each chain is a genuine eigenvector (it is in $\ker N$). The number of chains is $\dim\ker N$. A chain of length $s$ contributes one vector to $\ker N$, two to $\ker N^2$, and so on up to $s$: in the picture, $\ker N^k$ is spanned by the $k$ columns nearest the zero. So if you know $\dim\ker N^k$ for every $k$, you can count the chains of each length. That counting is the easy half of the module (Theorem 3.11). Proving that the chains always exist is the hard half (Theorem 3.6).

### Why a shift makes $t$ appear

Solve $\mathbf{y}' = J\mathbf{y}$ in chain coordinates for the $2\times 2$ block. The equations are
$$y_2' = \lambda y_2, \qquad y_1' = \lambda y_1 + y_2.$$
The top of the chain is on its own: $y_2 = c_2 e^{\lambda t}$. That solution is then *fed into* the equation for $y_1$ as a forcing term, and the forcing is at the natural rate of $y_1$, which is resonance: $y_1' - \lambda y_1 = c_2 e^{\lambda t}$ has the particular solution $c_2\, t\, e^{\lambda t}$. So $y_1 = (c_1 + c_2 t)e^{\lambda t}$. The $t$ is a resonance between a chain vector and the vector below it. A 3-chain feeds twice and gives $t^2 e^{\lambda t}/2$. You can guess the general formula now: $e^{Jt} = e^{\lambda t}(I + tN + \tfrac{t^2}{2}N^2 + \cdots)$, and the series stops because $N$ is nilpotent.

### What the theorem should say

Here is what you should be able to guess before §3. Take any $A$ on $\mathbb{C}^n$.

1. For each eigenvalue $\lambda$, instead of just $\ker(A - \lambda I)$ (the eigenvectors), take the bigger space $G_\lambda$ of vectors killed by *some* power of $A - \lambda I$. These are the generalized eigenvectors.
2. On $G_\lambda$, the operator $A - \lambda I$ is nilpotent, so it is a collection of chains.
3. The $G_\lambda$'s for different $\lambda$ fill the whole space with no overlap: $\mathbb{C}^n = G_{\lambda_1}\oplus\cdots\oplus G_{\lambda_m}$. (For a diagonalizable $A$ this is the eigenvector decomposition of LA-7; now it is true for every $A$.)
4. Put all the chains side by side as a basis. In that basis, $A$ is a block-diagonal matrix of Jordan blocks: **diagonal plus chains.**

Each step is where you'd expect a theorem to be, and each gets proved in §3. The thing to hold onto: an eigenvector is a chain of length 1; a Jordan block is a chain; diagonalizable means all chains have length 1.

### The one thing to be suspicious of

Diagonalizability is a *yes/no* property, and $\begin{pmatrix}\lambda&1\\ \varepsilon&\lambda\end{pmatrix}$ is diagonalizable for every $\varepsilon\ne0$ and not at $\varepsilon = 0$. So the Jordan form jumps discontinuously as the entries move. Floating point cannot see the difference between $\varepsilon = 0$ and $\varepsilon = 10^{-17}$. Keep that in mind throughout: this module is theory that explains ODEs, similarity, and matrix functions. It is not an algorithm you would run on data. §3.11 says what you run instead.

