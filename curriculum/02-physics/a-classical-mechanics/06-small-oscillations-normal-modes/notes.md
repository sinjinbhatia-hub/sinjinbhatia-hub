# CM-6: Small Oscillations and Normal Modes (the Generalized Eigenproblem)

**Time:** ~12 hours · **Prerequisites:** CM-3, LA-8 · **Unlocks:** CM-7, F2 (continuum limit), QM-3 (oscillator), SM-5 (Debye)

> **Conventions in this module.** Generalized coordinates $q = (q_1, \dots, q_n)^\top \in \mathbb{R}^n$ are measured from a stable equilibrium, so the equilibrium is $q = 0$. Lagrangian $L = T - V$ (course convention). The **mass matrix** $M$ is symmetric positive definite; the **stiffness matrix** $K$ is symmetric. The mode ansatz is $q(t) = \mathbf{a}\cos(\omega t - \delta)$ with real mode shape $\mathbf{a}$, angular frequency $\omega \ge 0$ and phase $\delta$. Modes are numbered $\omega_1 \le \omega_2 \le \cdots \le \omega_n$ and normalized by $\mathbf{a}_i^\top M\mathbf{a}_j = \delta_{ij}$, using the **$M$-inner product** $\langle\mathbf{u}, \mathbf{v}\rangle_M = \mathbf{u}^\top M\mathbf{v}$. Normal coordinates are $\xi_j$, with $q = A\xi$, $A = [\mathbf{a}_1 \cdots \mathbf{a}_n]$ and $\Omega^2 = \operatorname{diag}(\omega_1^2, \dots, \omega_n^2)$. Spring constants are $k$; because $k$ is taken, the wavenumber in the string example (4.5) is $\kappa$, and in that example, which has a tension $\tau$ and no springs, $k$ is the bead index. $\dagger$ is conjugate transpose, $\succ 0$ means positive definite. SI units.

## 0. Where this sits

From CM-3 you must be able to write $T$ and $V$ for a constrained system in generalized coordinates and get the equations of motion from Euler–Lagrange; the double pendulum of CM-3 reappears here as Example 4.4. From LA-8 you need the spectral theorem for symmetric matrices, what positive definite means and the second-derivative test, and the three axioms of an inner product. From LA-10 the Cholesky factorization $M = LL^\top$ is used once, at the decisive moment (Theorem 3.4); if you have not done LA-10, the one fact you need is stated there.

What this module makes possible: the inertia tensor of CM-7 is the same kind of object as $M$ here, a symmetric positive definite matrix that defines a kinetic energy, and principal axes are normal modes without the $K$. F2 gets its standing waves and its orthogonal sine basis from the loaded string of Example 4.5. QM-3 quantizes the single oscillator that every normal mode becomes; SM-5 counts those oscillators to get the heat capacity of a solid. For an ME student this module is **modal analysis**: it is what finite-element vibration software computes when you click "solve."

## 1. The problem that forced this idea

Anything with mass and stiffness vibrates: a bridge deck, a turbine blade, an engine on its mounts, the atoms in a molecule or a crystal. A system with $n$ degrees of freedom near equilibrium has $n$ coupled second-order ODEs, and when you write them out they look hopeless. Push on one mass and every other mass responds; the coordinates refuse to separate.

The first people to see through this were working on chains and strings. In 1733 Daniel Bernoulli studied a hanging chain, first as a few heavy beads on a flexible thread and then as a continuous chain, and found special motions in which every bead moves at the same frequency, each at its own fixed amplitude, the whole shape rising and falling in unison; there were as many such motions as beads, each with its own frequency. In the 1750s he claimed something stronger for the vibrating string: that *every* motion is a sum of these special motions. Euler doubted the generality of the claim, d'Alembert had solved the continuous string by a different route (1747), and the three argued for decades. Bernoulli was right, and the argument is the origin of Fourier series: the question "is every vibration a sum of modes?" was on the table half a century before Fourier. In 1759 Lagrange attacked the discrete version, $n$ equal beads on a weightless string (Example 4.5 derives its modes), and recovered the continuous string as a limit (LA-7 §1 tells this part).

Lagrange, in the *Mécanique analytique* (1788), turned the physical insight into algebra for any system. Expand $T$ and $V$ to second order about equilibrium, assume a solution in which every coordinate oscillates at the same frequency, and the differential equations collapse to a polynomial equation of degree $n$ in $\omega^2$ whose roots are the mode frequencies. He argued that for a stable equilibrium the roots are real and positive. He did all this in coordinates, with no matrices, because matrices did not exist yet; Cayley's matrix algebra came in 1858. Lagrange had diagonalized a symmetric matrix seventy years early. His treatment had one gap: he believed that a repeated root produces motions growing like $t$, and only in 1858 did Weierstrass show that with a positive definite kinetic energy it never does. That gap is the completeness-with-degeneracy clause of Theorem 3.4.

Here is why the problem belongs in this course and not only in a dynamics course. Once you see that "find the modes of a vibrating structure" is "find the eigenvectors of a self-adjoint operator in a weighted inner product," three things happen. You stop solving $n$ coupled ODEs and solve one eigenproblem. You get, for free, theorems about reality, stability, orthogonality and completeness that would be painful to prove one system at a time. And you have the template for quantum mechanics, where exactly the same algebra, with $K$ and $M$ replaced by a Hamiltonian, produces energy levels instead of frequencies.

For the engineer: this is modal analysis, and it is not a toy. A finite-element model of a car body has a million degrees of freedom. Nobody integrates those equations. The software assembles $K$ and $M$, computes the lowest few dozen modes, and every subsequent calculation (response to road input, fatigue, noise) is done in the space spanned by those modes. The theorems of §3 are why that works.

## 2. Intuition

### Two pendulums and a weak spring

Hang two identical pendulums side by side, length $\ell$, bob mass $m$, and connect the bobs by a weak spring of stiffness $k$ whose natural length equals the distance between the pivots. Each pendulum alone swings at $\omega_0 = \sqrt{g/\ell}$. Call the angles $\theta_1, \theta_2$.

```
     pivot                 pivot
       |                     |
       | l                   | l
       |                     |
      (m)/\/\/\/\/\/\/\/\/\(m)      spring k, relaxed when theta_1 = theta_2
    theta_1               theta_2
```

Two motions are easy to guess.

**In phase.** Start both bobs at the same angle and release. By symmetry they stay at the same angle forever: $\theta_1(t) = \theta_2(t)$. The spring never changes length, so it exerts no force, so each pendulum swings exactly as if the spring were absent, at $\omega_0$. The spring is invisible to this motion.

**Out of phase.** Start the bobs at opposite angles, $\theta_1 = -\theta_2$, and release. By the left–right symmetry they stay opposite forever. Now the spring is stretched by about $2\ell\theta_1$, so each bob feels the full spring force on top of gravity, and the motion is faster than $\omega_0$. Newton on bob 1 along its arc: $m\ell\ddot\theta_1 = -mg\theta_1 - k(2\ell\theta_1)$, so you would guess $\omega^2 = g/\ell + 2k/m$. (Example 4.1 confirms it.)

Each of these motions is a **normal mode**: every coordinate oscillates at one frequency, and the ratio $\theta_1 : \theta_2$ is fixed in time. The pattern $(1, 1)$ or $(1, -1)$ is the **mode shape**. The in-phase shape is the one the spring cannot see; the out-of-phase shape is the one it works hardest against.

### Every motion is a sum of the two

Now release pendulum 1 from angle $\theta_0$ with pendulum 2 hanging still. This is neither mode. But the initial condition $(\theta_0, 0)$ is half of $(\theta_0, \theta_0)$ plus half of $(\theta_0, -\theta_0)$: half an in-phase mode plus half an out-of-phase mode. Each half evolves at its own frequency, and because the equations are linear the sum of the two evolutions is the motion:

$$
\theta_1(t) = \frac{\theta_0}{2}\left(\cos\omega_1 t + \cos\omega_2 t\right), \qquad
\theta_2(t) = \frac{\theta_0}{2}\left(\cos\omega_1 t - \cos\omega_2 t\right).
$$

Since $\omega_2$ is only slightly bigger than $\omega_1$ for a weak spring, the two cosines drift slowly in and out of step. When they are in step, $\theta_1$ swings at full amplitude and $\theta_2$ is still. Half a drift cycle later they are opposite and all the motion has moved to pendulum 2. Energy sloshes back and forth with a period set by $\omega_2 - \omega_1$. These are **beats**, and you can build this on a desk with two pencils, two washers and a rubber band.

### The guess you should make now

Write the linearized equations of the coupled pendulums as one matrix equation. With $q = (\theta_1, \theta_2)^\top$:

$$
m\ell^2\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}\ddot q = -\begin{pmatrix} mg\ell + k\ell^2 & -k\ell^2\\ -k\ell^2 & mg\ell + k\ell^2\end{pmatrix} q, \qquad\text{i.e.}\qquad M\ddot q = -Kq.
$$

A mode is a motion $q(t) = \mathbf{a}\cos(\omega t - \delta)$ with a fixed shape $\mathbf{a}$. Put it in: $\ddot q = -\omega^2 q$, so $-\omega^2 M\mathbf{a} = -K\mathbf{a}$, so

$$
K\mathbf{a} = \omega^2 M\mathbf{a}.
$$

Here $M$ is a multiple of the identity, so this is an ordinary eigenproblem for $K$: the mode shapes are the eigenvectors of $K$, and the $\omega^2$ are its eigenvalues divided by $m\ell^2$. The eigenvectors of the symmetric matrix $K$ above are $(1, 1)^\top$ and $(1, -1)^\top$ (check: $K(1,1)^\top = mg\ell\,(1,1)^\top$ and $K(1,-1)^\top = (mg\ell + 2k\ell^2)(1,-1)^\top$), exactly the two shapes we guessed, and they are orthogonal, as LA-8 says eigenvectors of a symmetric matrix for distinct eigenvalues must be.

You met this in LA-7: §3.10 turned $M\ddot{\mathbf x} = -K\mathbf x$ into $K\mathbf v = \omega^2M\mathbf v$, and Example 4.3 solved the two-mass, three-spring chain with masses $1$ and $2$ and noticed that its two mode shapes are orthogonal only after weighting by the masses. LA-7 observed that; it did not explain it. The new thing in this module is the explanation, and the cases a $2\times2$ determinant cannot handle: $M$ not even diagonal (the double pendulum), zero modes, repeated frequencies. When $M$ is not a multiple of the identity, $K\mathbf{a} = \omega^2 M\mathbf{a}$ is a **generalized eigenproblem**. The obvious move is to rewrite it as $M^{-1}K\mathbf{a} = \omega^2\mathbf{a}$ and apply the spectral theorem to $M^{-1}K$. Try it on $M = \operatorname{diag}(1, 2)$, $K = \begin{pmatrix}2 & -1\\-1 & 2\end{pmatrix}$: $M^{-1}K = \begin{pmatrix}2 & -1\\-1/2 & 1\end{pmatrix}$, not symmetric. The spectral theorem does not apply to it as a matrix, and its eigenvectors turn out *not* to be orthogonal in the ordinary sense (Example 4.2: their dot product is $1/2$).

So here is the guess to make before §3. The mode shapes are eigenvectors of *something* self-adjoint, and they are orthogonal in *some* sense. What has to change is the ruler. The right inner product on configuration space is not $\mathbf{u}^\top\mathbf{v}$ but $\mathbf{u}^\top M\mathbf{v}$, the one defined by the kinetic energy. In that inner product the mode shapes are orthogonal, $M^{-1}K$ is self-adjoint, and every theorem of LA-8 applies. §3 proves this.

## 3. Formal development

### 3.1 Expanding the Lagrangian about a stable equilibrium

Let the system have $n$ generalized coordinates $q \in \mathbb{R}^n$ and $N$ particles at positions $\mathbf{r}_\alpha(q)$, $\alpha = 1, \dots, N$ (holonomic, time-independent constraints, as in CM-3). The potential is $V(q)$.

**Derivation 1 (the structure of $T$).** ▶ **Try it first.** Differentiate $\mathbf{r}_\alpha(q(t))$ by the chain rule and substitute into $T = \frac12\sum_\alpha m_\alpha|\dot{\mathbf{r}}_\alpha|^2$. What kind of function of $\dot q$ do you get?

1. Chain rule: $\dot{\mathbf{r}}_\alpha = \sum_{i=1}^n \dfrac{\partial\mathbf{r}_\alpha}{\partial q_i}\dot q_i$.
2. Square and sum: $T = \dfrac12\sum_\alpha m_\alpha\,\dot{\mathbf{r}}_\alpha\cdot\dot{\mathbf{r}}_\alpha = \dfrac12\sum_{i,j}\dot q_i\dot q_j\,\underbrace{\sum_\alpha m_\alpha\dfrac{\partial\mathbf{r}_\alpha}{\partial q_i}\cdot\dfrac{\partial\mathbf{r}_\alpha}{\partial q_j}}_{M_{ij}(q)}$.
3. So $T = \frac12\dot q^\top M(q)\dot q$ with $M_{ij}(q)$ as marked. $M$ is symmetric (the dot product is), and $\dot q^\top M\dot q = \sum_\alpha m_\alpha|\dot{\mathbf{r}}_\alpha|^2 \ge 0$, with equality only if every particle is at rest. If the coordinates are independent (the $n$ vectors $\partial\mathbf{r}/\partial q_i \in \mathbb{R}^{3N}$ are linearly independent), then any $\dot q \ne 0$ moves some particle and $T > 0$. So $M(q)$ is **positive definite**. ∎

■ **What just happened.** Kinetic energy is a quadratic form in the velocities, with a configuration-dependent, symmetric, positive definite matrix. Positive definiteness is not an assumption; it is the statement that moving costs energy.

Now suppose $q = 0$ is an **equilibrium**: $\partial V/\partial q_i(0) = 0$ for all $i$. (If the equilibrium is elsewhere, shift coordinates.) Taylor-expand both pieces of $L$ to second order in the small quantities $q$ and $\dot q$:

$$
V(q) = V_0 + \underbrace{\sum_i\frac{\partial V}{\partial q_i}(0)\,q_i}_{=\,0} + \frac12\sum_{i,j}\frac{\partial^2V}{\partial q_i\partial q_j}(0)\,q_iq_j + O(q^3) = V_0 + \frac12 q^\top Kq + O(q^3), \qquad K_{ij} = \frac{\partial^2V}{\partial q_i\partial q_j}(0). \tag{3.1}
$$

$K$ is the **Hessian** of $V$ at the equilibrium (LA-14), symmetric because mixed partials commute. For $T$, write $M(q) = M(0) + O(q)$; then $T = \frac12\dot q^\top M(0)\dot q + O(q\,\dot q^2)$, and the correction is *third* order in small quantities. Dropping everything above second order and the constant $V_0$:

$$
L_2 = \frac12\dot q^\top M\dot q - \frac12 q^\top Kq, \qquad M := M(0). \tag{3.2}
$$

**Derivation 2 (the linearized equations of motion).** ▶ **Try it first.** Apply Euler–Lagrange to (3.2). You will need $\partial(\dot q^\top M\dot q)/\partial\dot q_i = 2\sum_j M_{ij}\dot q_j$; prove that from the symmetry of $M$.

1. $\dfrac{\partial}{\partial\dot q_i}\Big(\dfrac12\sum_{j,l}M_{jl}\dot q_j\dot q_l\Big) = \dfrac12\Big(\sum_l M_{il}\dot q_l + \sum_j M_{ji}\dot q_j\Big) = \sum_j M_{ij}\dot q_j$, using $M_{ji} = M_{ij}$. In vector form, $\partial L_2/\partial\dot q = M\dot q$.
2. The same computation on the potential term gives $\partial L_2/\partial q = -Kq$.
3. Euler–Lagrange, $\dfrac{d}{dt}\dfrac{\partial L_2}{\partial\dot q} = \dfrac{\partial L_2}{\partial q}$, with $M$ constant:
$$M\ddot q = -Kq. \tag{3.3}$$ ∎

■ **What just happened.** Near any equilibrium, any Lagrangian system is a linear system of $n$ coupled oscillators, with the Hessian of $V$ as the stiffness and the kinetic-energy matrix as the mass. Everything that follows uses only (3.3) and the two facts $M = M^\top \succ 0$, $K = K^\top$.

**Where $K$ comes from in a spring network.** When $V$ is a sum of spring energies, $V = \frac12\sum_s c_s e_s^2$ with $e_s$ the elongation of spring $s$, and the elongations are linear in $q$ to first order, $e = Bq$ for a matrix $B$ with one row per spring, then
$$V = \tfrac12 q^\top B^\top CBq, \qquad K = B^\top CB, \quad C = \operatorname{diag}(c_1, c_2, \dots).$$
This is the $A^\top CA$ structure of Strang's applied-mathematics framework (§10): $K$ is automatically symmetric and positive semidefinite ($q^\top Kq = \sum_s c_s e_s^2 \ge 0$), and it is positive definite exactly when $Bq = 0$ forces $q = 0$, that is, when no nonzero displacement leaves every spring unstretched. A displacement in the null space $N(B)$ (LA-3) is a **rigid motion** or a **mechanism**, and §3.5 shows it is a zero-frequency mode. Examples 4.3 and 4.6 and the project build $K$ this way.

**Stability.** If $K$ is positive definite, $V$ has a strict local minimum at $q = 0$ (LA-8's second-derivative test), and the equilibrium is stable: the exact energy $E = T + V$ is conserved and $\ge V_0$, and the region $\{E \le V_0 + \epsilon\}$ of phase space shrinks onto the point $(q, \dot q) = (0, 0)$ as $\epsilon \to 0$, so a motion that starts near equilibrium with little energy stays near it. (This is the Lagrange–Dirichlet theorem; the argument is exact, not linearized.) If $K$ has a negative eigenvalue, Theorem 3.5 produces a solution of (3.3) that grows exponentially, and the equilibrium is unstable. If $K$ is only positive semidefinite, the linear theory cannot decide; the higher-order terms of $V$ do (cubic, quartic, or none at all when the flat direction is a symmetry, §3.5).

### 3.2 The ansatz and the generalized eigenproblem

Look for solutions of (3.3) in which every coordinate oscillates at the same frequency with a fixed shape:

$$
q(t) = \mathbf{a}\cos(\omega t - \delta), \qquad \mathbf{a} \in \mathbb{R}^n \text{ constant}. \tag{3.4}
$$

Then $\ddot q = -\omega^2\mathbf{a}\cos(\omega t - \delta)$, and (3.3) becomes $-\omega^2M\mathbf{a}\cos(\omega t - \delta) = -K\mathbf{a}\cos(\omega t - \delta)$. The cosine is not identically zero, so

$$
K\mathbf{a} = \omega^2 M\mathbf{a}. \tag{3.5}
$$

This is the **generalized eigenvalue problem** for the pair $(K, M)$: $\lambda = \omega^2$ is a **generalized eigenvalue** and $\mathbf{a}$ a **generalized eigenvector**. It has a nonzero solution $\mathbf{a}$ exactly when $K - \omega^2M$ is singular:

$$
\det(K - \omega^2 M) = 0. \tag{3.6}
$$

The left side is a polynomial of degree $n$ in $\omega^2$ (the coefficient of $(\omega^2)^n$ is $(-1)^n\det M \ne 0$), so there are $n$ roots counted with multiplicity. A solution (3.4) with $\omega^2$ a root and $\mathbf{a}$ a corresponding null vector is a **normal mode**; $\omega$ is its **natural frequency** and $\mathbf{a}$ its **mode shape**.

Using the complex ansatz $q = \mathbf{a}e^{i\omega t}$ gives the same (3.5), and the real part of $\mathbf{a}e^{i(\omega t - \delta)}$ is (3.4). Nothing in (3.5) yet says that $\omega^2$ is real or that $\mathbf{a}$ can be chosen real. That is the first theorem.

### 3.3 Five theorems

**Theorem 3.1 (Real frequencies).** Let $M = M^\top \succ 0$ and $K = K^\top$. Every generalized eigenvalue $\omega^2$ of (3.5) is real, and every eigenvalue has a real eigenvector.

▶ **Try it first.** Allow $\mathbf{a} \in \mathbb{C}^n$ and $\omega^2 \in \mathbb{C}$. Multiply (3.5) on the left by $\mathbf{a}^\dagger$. Show that both $\mathbf{a}^\dagger K\mathbf{a}$ and $\mathbf{a}^\dagger M\mathbf{a}$ are real numbers.

*Proof.*
1. From (3.5), $\mathbf{a}^\dagger K\mathbf{a} = \omega^2\,\mathbf{a}^\dagger M\mathbf{a}$.
2. For any real symmetric $S$ and complex $\mathbf{a}$, the number $z = \mathbf{a}^\dagger S\mathbf{a}$ is real: $\bar z = \mathbf{a}^\top S\bar{\mathbf{a}}$, and a $1\times 1$ matrix equals its transpose, so $\bar z = (\mathbf{a}^\top S\bar{\mathbf{a}})^\top = \bar{\mathbf{a}}^\top S^\top\mathbf{a} = \mathbf{a}^\dagger S\mathbf{a} = z$.
3. Writing $\mathbf{a} = \mathbf{u} + i\mathbf{v}$ with $\mathbf{u}, \mathbf{v}$ real, $\mathbf{a}^\dagger M\mathbf{a} = \mathbf{u}^\top M\mathbf{u} + \mathbf{v}^\top M\mathbf{v} + i(\mathbf{u}^\top M\mathbf{v} - \mathbf{v}^\top M\mathbf{u})$. The imaginary part vanishes by symmetry of $M$, and the real part is $> 0$ for $\mathbf{a} \ne 0$ because $M$ is positive definite and at least one of $\mathbf{u}, \mathbf{v}$ is nonzero.
4. So $\omega^2 = \mathbf{a}^\dagger K\mathbf{a}\,/\,\mathbf{a}^\dagger M\mathbf{a}$ is a real number divided by a positive one: real.
5. With $\omega^2$ real, $K - \omega^2M$ is a real singular matrix, so it has a nonzero real null vector. ∎

■ **What just happened.** The ratio in step 4 is the **generalized Rayleigh quotient** $R(\mathbf{a}) = \mathbf{a}^\top K\mathbf{a}/\mathbf{a}^\top M\mathbf{a}$: twice the potential energy stored by the shape, divided by twice the kinetic energy the shape would have at unit rate. Every mode frequency is this ratio at its own shape, which is why a stiff shape (large $\mathbf{a}^\top K\mathbf{a}$) that moves little mass (small $\mathbf{a}^\top M\mathbf{a}$) is fast. It is LA-8's Rayleigh quotient with the ruler changed.

**Theorem 3.2 ($M$-orthogonality).** If $\mathbf{a}_i$ and $\mathbf{a}_j$ are eigenvectors of (3.5) with $\omega_i^2 \ne \omega_j^2$, then
$$\mathbf{a}_i^\top M\mathbf{a}_j = 0 \quad\text{and}\quad \mathbf{a}_i^\top K\mathbf{a}_j = 0. \tag{3.7}$$

▶ **Try it first.** This is LA-8's proof that eigenvectors of a symmetric matrix for distinct eigenvalues are orthogonal. Write that proof out, then insert $M$ where it belongs.

*Proof.*
1. $K\mathbf{a}_i = \omega_i^2 M\mathbf{a}_i$. Multiply on the left by $\mathbf{a}_j^\top$: $\mathbf{a}_j^\top K\mathbf{a}_i = \omega_i^2\,\mathbf{a}_j^\top M\mathbf{a}_i$.
2. $K\mathbf{a}_j = \omega_j^2M\mathbf{a}_j$. Multiply by $\mathbf{a}_i^\top$: $\mathbf{a}_i^\top K\mathbf{a}_j = \omega_j^2\,\mathbf{a}_i^\top M\mathbf{a}_j$.
3. Both $K$ and $M$ are symmetric, so $\mathbf{a}_j^\top K\mathbf{a}_i = \mathbf{a}_i^\top K\mathbf{a}_j$ and $\mathbf{a}_j^\top M\mathbf{a}_i = \mathbf{a}_i^\top M\mathbf{a}_j$ (a $1 \times 1$ matrix equals its transpose). Subtract step 2 from step 1: $0 = (\omega_i^2 - \omega_j^2)\,\mathbf{a}_i^\top M\mathbf{a}_j$.
4. Since $\omega_i^2 \ne \omega_j^2$, $\mathbf{a}_i^\top M\mathbf{a}_j = 0$. Then step 2 gives $\mathbf{a}_i^\top K\mathbf{a}_j = 0$ too. ∎

■ **What just happened.** The symmetric-matrix orthogonality proof goes through word for word, except that the thing that vanishes is $\mathbf{a}_i^\top M\mathbf{a}_j$, not $\mathbf{a}_i^\top\mathbf{a}_j$. Mode shapes for different frequencies are orthogonal *when you measure with $M$*. Physically: the kinetic energy of a superposition of two modes has no cross term. For a free molecule (Example 4.3) with the translational zero mode $\mathbf{a}_0 = (1, \dots, 1)^\top$ and diagonal $M$, (3.7) reads $\sum_i m_i a_i = 0$ for every vibrational mode: a vibration never moves the center of mass.

**The inner-product view.** Define, for $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$,
$$\langle\mathbf{u}, \mathbf{v}\rangle_M = \mathbf{u}^\top M\mathbf{v}. \tag{3.8}$$

**Lemma 3.3.** $\langle\cdot, \cdot\rangle_M$ is an inner product on $\mathbb{R}^n$, in LA-8's sense, if and only if $M$ is symmetric positive definite.

*Proof.* Check the three axioms. *Bilinear*: $\mathbf{u}^\top M\mathbf{v}$ is linear in $\mathbf{u}$ and in $\mathbf{v}$ by the rules of matrix multiplication, for any $M$. *Symmetric*: $\langle\mathbf{v}, \mathbf{u}\rangle_M = \mathbf{v}^\top M\mathbf{u} = (\mathbf{v}^\top M\mathbf{u})^\top = \mathbf{u}^\top M^\top\mathbf{v}$, which equals $\langle\mathbf{u}, \mathbf{v}\rangle_M$ for all $\mathbf{u}, \mathbf{v}$ iff $M^\top = M$ (take $\mathbf{u} = \mathbf{e}_i$, $\mathbf{v} = \mathbf{e}_j$ to see that the entries must agree). *Positive definite*: $\langle\mathbf{v}, \mathbf{v}\rangle_M = \mathbf{v}^\top M\mathbf{v} > 0$ for all $\mathbf{v} \ne 0$ is exactly the definition of $M \succ 0$. ∎

So by Derivation 1 the kinetic energy hands us a legitimate inner product: $T = \frac12\langle\dot q, \dot q\rangle_M$, and $\|\mathbf{v}\|_M = \sqrt{\langle\mathbf{v},\mathbf{v}\rangle_M}$ is a norm (LA-1's C1 with a full matrix of weights). Now rewrite (3.5) as $M^{-1}K\mathbf{a} = \omega^2\mathbf{a}$ and compute, for any $\mathbf{u}, \mathbf{v}$:
$$\langle\mathbf{u}, M^{-1}K\mathbf{v}\rangle_M = \mathbf{u}^\top M M^{-1}K\mathbf{v} = \mathbf{u}^\top K\mathbf{v} = (K\mathbf{u})^\top\mathbf{v} = (M M^{-1}K\mathbf{u})^\top\mathbf{v} = (M^{-1}K\mathbf{u})^\top M\mathbf{v} = \langle M^{-1}K\mathbf{u}, \mathbf{v}\rangle_M.$$
The operator $M^{-1}K$ is **self-adjoint with respect to $\langle\cdot,\cdot\rangle_M$**, even though its matrix in the standard basis is not symmetric. The spectral theorem is at heart a statement about self-adjoint operators on an inner product space, not about symmetric arrays of numbers. LA-8's Theorems 3.9 and 3.10 are already stated that way, and the proof of its Theorem 3.11 uses symmetry only to let the operator hop across the inner product, which $M^{-1}K$ does in $\langle\cdot,\cdot\rangle_M$ (LA-8 Exercise B1(d)). Applied here, it says that $M^{-1}K$ has real eigenvalues and an eigenbasis orthonormal in $\langle\cdot,\cdot\rangle_M$. Theorems 3.1 and 3.2 are what that abstract statement says in coordinates. The next theorem makes the connection to a plain symmetric matrix concrete, which is also how the computer does it.

**Theorem 3.4 (Reduction to a standard symmetric eigenproblem; completeness).** Let $M = LL^\top$ be the Cholesky factorization of $M$ ($L$ lower triangular with positive diagonal; it exists and is unique because $M \succ 0$, LA-10). Define
$$\tilde K = L^{-1}KL^{-\top}, \qquad L^{-\top} := (L^{-1})^\top = (L^\top)^{-1}. \tag{3.9}$$
Then (i) $\tilde K$ is symmetric; (ii) $\omega^2$ is a generalized eigenvalue of $(K, M)$ with eigenvector $\mathbf{a}$ if and only if $\omega^2$ is an ordinary eigenvalue of $\tilde K$ with eigenvector $\mathbf{y} = L^\top\mathbf{a}$; (iii) consequently there exist $n$ real $\omega_1^2 \le \dots \le \omega_n^2$ (with multiplicity) and eigenvectors $\mathbf{a}_1, \dots, \mathbf{a}_n$ with
$$\mathbf{a}_i^\top M\mathbf{a}_j = \delta_{ij}, \tag{3.10}$$
and they form a basis of $\mathbb{R}^n$. This holds whether or not the $\omega_j^2$ are distinct.

▶ **Try it first.** Substitute $M = LL^\top$ into (3.5), multiply by $L^{-1}$, and insert $I = L^{-\top}L^\top$ between $K$ and $\mathbf{a}$.

*Proof.*
1. (i) $\tilde K^\top = (L^{-1}KL^{-\top})^\top = (L^{-\top})^\top K^\top (L^{-1})^\top = L^{-1}KL^{-\top} = \tilde K$, using $K^\top = K$.
2. (ii) $K\mathbf{a} = \omega^2LL^\top\mathbf{a}$. Multiply on the left by $L^{-1}$ and write $\mathbf{a} = L^{-\top}(L^\top\mathbf{a}) = L^{-\top}\mathbf{y}$ on the left side: $L^{-1}KL^{-\top}\mathbf{y} = \omega^2L^\top\mathbf{a} = \omega^2\mathbf{y}$, i.e. $\tilde K\mathbf{y} = \omega^2\mathbf{y}$. Each step is reversible because $L$ is invertible, and $\mathbf{y} \ne 0 \iff \mathbf{a} \ne 0$.
3. (iii) $\tilde K$ is real symmetric, so by the spectral theorem (LA-8) it has real eigenvalues $\omega_1^2 \le \dots \le \omega_n^2$ and an orthonormal eigenbasis $\mathbf{y}_1, \dots, \mathbf{y}_n$, $\mathbf{y}_i^\top\mathbf{y}_j = \delta_{ij}$, regardless of repeated eigenvalues. Set $\mathbf{a}_j = L^{-\top}\mathbf{y}_j$. By (ii) these solve (3.5), and
$$\mathbf{a}_i^\top M\mathbf{a}_j = \mathbf{y}_i^\top L^{-1}(LL^\top)L^{-\top}\mathbf{y}_j = \mathbf{y}_i^\top\mathbf{y}_j = \delta_{ij}.$$
4. The $\mathbf{a}_j$ are linearly independent: if $\sum_j c_j\mathbf{a}_j = 0$, apply $\langle\mathbf{a}_i, \cdot\rangle_M$ to both sides and use (3.10) to get $c_i = 0$. So they are $n$ independent vectors in $\mathbb{R}^n$: a basis. ∎

■ **What just happened.** The change of variable $\mathbf{y} = L^\top\mathbf{a}$ is a change of coordinates in which the kinetic energy becomes the ordinary Euclidean one: $\dot q^\top M\dot q = \dot q^\top LL^\top\dot q = |L^\top\dot q|^2$. In those coordinates the mass matrix is the identity, the problem is a standard symmetric eigenproblem, and orthogonality is ordinary orthogonality. Translating back, ordinary orthogonality of the $\mathbf{y}$'s is $M$-orthogonality of the $\mathbf{a}$'s. This is exactly what `scipy.linalg.eigh(K, M)` does (LAPACK's `sygv` family of drivers: Cholesky, transform, symmetric eigensolver, transform back), and it is why LA-10's Cholesky factorization matters for anything that vibrates. Any factorization $M = RR^\top$ works, including the symmetric square root $M^{1/2}$ from LA-8; Cholesky is the cheapest.

Completeness with degeneracies is the part Theorem 3.2 could not give you. Theorem 3.2 only speaks about distinct frequencies; Theorem 3.4 says that even when $\omega_i^2 = \omega_j^2$ an $M$-orthonormal pair of shapes exists in that eigenspace (Gram–Schmidt in the $M$-inner product produces one). Every motion is a sum of modes because the modes span $\mathbb{R}^n$. Bernoulli's claim is proved, for finitely many degrees of freedom.

**Theorem 3.5 (Stability).** With $M$, $K$ as above and $\omega_j^2$ the generalized eigenvalues: all $\omega_j^2 \ge 0$ if and only if $K$ is positive semidefinite, and all $\omega_j^2 > 0$ if and only if $K$ is positive definite.

▶ **Try it first.** One direction is a line, using a formula from the proof of Theorem 3.1. For the other, expand an arbitrary $\mathbf{v} \ne 0$ in the modes of Theorem 3.4 and compute $\mathbf{v}^\top K\mathbf{v}$.

*Proof.* ($\Leftarrow$) Each $\omega_j^2 = \mathbf{a}_j^\top K\mathbf{a}_j/\mathbf{a}_j^\top M\mathbf{a}_j$ by Theorem 3.1 step 4; the denominator is positive, and the numerator is $\ge 0$ ($> 0$) when $K$ is PSD (PD). ($\Rightarrow$) Take any $\mathbf{v} \ne 0$ and expand it in the basis of Theorem 3.4: $\mathbf{v} = \sum_j c_j\mathbf{a}_j$ with $c_j = \langle\mathbf{a}_j, \mathbf{v}\rangle_M$, not all zero. Using $K\mathbf{a}_j = \omega_j^2M\mathbf{a}_j$ and (3.10),
$$\mathbf{v}^\top K\mathbf{v} = \sum_{i,j}c_ic_j\,\mathbf{a}_i^\top K\mathbf{a}_j = \sum_{i,j}c_ic_j\,\omega_j^2\,\mathbf{a}_i^\top M\mathbf{a}_j = \sum_j\omega_j^2c_j^2,$$
which is $\ge 0$ if all $\omega_j^2 \ge 0$, and $> 0$ if all $\omega_j^2 > 0$. ∎

When some $\omega_j^2 < 0$, write $\omega_j^2 = -\gamma_j^2$ with $\gamma_j > 0$. The ansatz $q = \mathbf{a}_je^{st}$ in (3.3) gives $s^2M\mathbf{a}_j = -K\mathbf{a}_j = \gamma_j^2M\mathbf{a}_j$, so $s = \pm\gamma_j$ and $q = \mathbf{a}_j(c_1e^{\gamma_jt} + c_2e^{-\gamma_jt})$: almost every initial condition grows exponentially. That is the linear-theory statement that a saddle of $V$ is unstable.

■ **What just happened.** In the mode basis, $\mathbf{v}^\top K\mathbf{v} = \sum_j\omega_j^2c_j^2$: the potential energy of any shape is a weighted sum of squares with the $\omega_j^2$ as weights. So the sign pattern of the $\omega_j^2$ is the sign pattern of $K$, and $M$, being positive definite, never changes a sign. A negative $\omega_j^2$ is a mode shape along which $V$ curves down, and the motion runs away along it.

### 3.4 Normal coordinates and the general solution

Collect the $M$-orthonormal modes of Theorem 3.4 into the **modal matrix** $A = [\mathbf{a}_1 \cdots \mathbf{a}_n]$ and let $\Omega^2 = \operatorname{diag}(\omega_1^2, \dots, \omega_n^2)$. The $n$ equations (3.5) stack into $KA = MA\Omega^2$, and (3.10) is
$$A^\top MA = I, \qquad A^\top KA = A^\top MA\Omega^2 = \Omega^2. \tag{3.11}$$
From $A^\top MA = I$ we read off the inverse: $A^{-1} = A^\top M$. (Not $A^\top$. $A$ is not an orthogonal matrix unless $M = I$.)

**Definition.** The **normal coordinates** $\xi \in \mathbb{R}^n$ are defined by $q = A\xi$, i.e. $\xi = A^\top Mq$, $\xi_j = \langle\mathbf{a}_j, q\rangle_M$.

**Derivation 3 (decoupling).** ▶ **Try it first.** Substitute $q = A\xi$ into (3.2) and use (3.11).

1. $T = \frac12\dot q^\top M\dot q = \frac12\dot\xi^\top A^\top MA\dot\xi = \frac12\dot\xi^\top\dot\xi = \frac12\sum_j\dot\xi_j^2$.
2. $V_2 = \frac12 q^\top Kq = \frac12\xi^\top A^\top KA\xi = \frac12\xi^\top\Omega^2\xi = \frac12\sum_j\omega_j^2\xi_j^2$.
3. So
$$L_2 = \frac12\sum_{j=1}^n\left(\dot\xi_j^2 - \omega_j^2\xi_j^2\right), \tag{3.12}$$
and Euler–Lagrange gives $n$ independent equations $\ddot\xi_j = -\omega_j^2\xi_j$. ∎

■ **What just happened.** In normal coordinates the system is $n$ uncoupled unit-mass oscillators. The coupling was never physical; it was an artifact of describing the motion in coordinates that were not the eigenvectors. The energy $E = \frac12\sum_j(\dot\xi_j^2 + \omega_j^2\xi_j^2)$ is a sum of $n$ separately conserved modal energies (Exercise B4).

**Derivation 4 (general solution and initial conditions).** ▶ **Try it first.** Solve $\ddot\xi_j = -\omega_j^2\xi_j$ with $\xi_j(0), \dot\xi_j(0)$ given; then express $\xi_j(0)$ and $\dot\xi_j(0)$ in terms of $q(0) = q_0$ and $\dot q(0) = v_0$.

1. For $\omega_j > 0$: $\xi_j(t) = \xi_j(0)\cos\omega_jt + \dfrac{\dot\xi_j(0)}{\omega_j}\sin\omega_jt$. For $\omega_j = 0$: $\ddot\xi_j = 0$, so $\xi_j(t) = \xi_j(0) + \dot\xi_j(0)\,t$.
2. From $\xi = A^\top Mq$ at $t = 0$: $\xi_j(0) = \mathbf{a}_j^\top Mq_0 = \langle\mathbf{a}_j, q_0\rangle_M$ and $\dot\xi_j(0) = \mathbf{a}_j^\top Mv_0 = \langle\mathbf{a}_j, v_0\rangle_M$. These are the components of $q_0$ and $v_0$ along an orthonormal basis: projections, exactly as in LA-5, with the $M$-inner product.
3. Assemble $q = A\xi$:
$$q(t) = \sum_{j=1}^n\mathbf{a}_j\Big[\langle\mathbf{a}_j, q_0\rangle_M\cos\omega_jt + \langle\mathbf{a}_j, v_0\rangle_M\,\frac{\sin\omega_jt}{\omega_j}\Big], \tag{3.13}$$
with $\frac{\sin\omega_jt}{\omega_j}$ read as $t$ when $\omega_j = 0$ (its limit). Equivalently each mode contributes $\mathbf{a}_jC_j\cos(\omega_jt - \delta_j)$ with $C_j = \sqrt{\xi_j(0)^2 + \dot\xi_j(0)^2/\omega_j^2}$ and $\tan\delta_j = \dot\xi_j(0)/(\omega_j\xi_j(0))$: the form (3.4), one amplitude and one phase per mode, $2n$ constants for $2n$ initial conditions. ∎

■ **What just happened.** To predict the motion from any start: project the initial position and velocity onto each mode *with the $M$-inner product*, let each projection oscillate at its own frequency, add. If you project with the ordinary dot product you get the wrong answer whenever $M \ne cI$. (3.13) is the modal superposition formula that vibration software evaluates.

### 3.5 Zero modes

If $K$ is only positive semidefinite, its null space $N(K)$ is nontrivial, and any $\mathbf{a}_0 \in N(K)$ satisfies $K\mathbf{a}_0 = 0 = 0\cdot M\mathbf{a}_0$: it is a mode with $\omega = 0$, a **zero mode**. Its "oscillation" is $\xi_0(t) = \xi_0(0) + \dot\xi_0(0)\,t$: uniform drift, no restoring force.

Where they come from. A zero mode is a direction in configuration space along which $V$ is flat to second order. For a free system (a molecule, a satellite, the triangle of Example 4.6) the obvious ones are rigid translations and rotations: $V$ depends only on relative positions, so moving everything together costs nothing. In three dimensions a free body of $N$ atoms has 3 translational and 3 rotational zero modes (2 rotational if the atoms are collinear, since rotation about the axis moves nothing), leaving $3N - 6$ (or $3N - 5$) vibrational modes. CO₂ has $N = 3$, linear: $3\cdot3 - 5 = 4$ vibrations, two stretches and a doubly degenerate bend. In spring-network language zero modes are $N(B)$: displacements that stretch no spring. A structure with a zero mode that is not a rigid motion is a mechanism, and it folds.

Noether's theorem (CM-4) is the other face of the same fact. Translation invariance of $V$ gives a zero mode $\mathbf{a}_0$ and also a conserved momentum; the conserved quantity is $\dot\xi_0 = \langle\mathbf{a}_0, \dot q\rangle_M$, the velocity of the zero mode. For $\mathbf{a}_0 \propto (1, \dots, 1)^\top$ and diagonal $M$ this is $\sum_im_i\dot q_i$, the total momentum.

Two practical consequences. The general solution (3.13) contains terms linear in $t$: a free molecule given a push drifts away while vibrating. And every formula that divides by $\omega_j$ breaks on a zero mode; the project's `fit_initial_conditions` handles the case explicitly, and §8 says more.

### 3.6 Degenerate modes and symmetry

If $\omega_i^2 = \omega_j^2$, then any combination $c_i\mathbf{a}_i + c_j\mathbf{a}_j$ is also an eigenvector with that eigenvalue (linearity of (3.5)), so "the mode shape" is not defined: the whole two-dimensional eigenspace is. Theorem 3.4 guarantees an $M$-orthonormal basis of it; the basis is not unique. Pick the one that makes physical sense, and expect the computer to pick a different one.

Degeneracy is almost always caused by symmetry. Suppose an invertible matrix $S$ is a **symmetry** of the system in the sense that $S^\top MS = M$ and $S^\top KS = K$ (relabeling the atoms of a symmetric molecule, rotating the triangle of Example 4.6 by $120°$, reflecting the coupled pendulums left to right). Then $S$ preserves $\langle\cdot,\cdot\rangle_M$ and $L_2$, and if $\mathbf{a}$ is a mode with frequency $\omega$ then so is $S\mathbf{a}$: from $S^\top KS = K$ we get $KS = S^{-\top}K$, and likewise $MS = S^{-\top}M$, so
$$K(S\mathbf{a}) = S^{-\top}K\mathbf{a} = \omega^2S^{-\top}M\mathbf{a} = \omega^2M(S\mathbf{a}).$$
The symmetry maps each eigenspace into itself. If the eigenspace is one-dimensional, $S\mathbf{a} = \pm\mathbf{a}$ (for a reflection, $S^2 = I$): the mode is itself symmetric or antisymmetric, like the in-phase and out-of-phase pendulum modes under left–right reflection, and the stretches of CO₂. If $S\mathbf{a}$ is not a multiple of $\mathbf{a}$, the eigenspace has dimension at least 2 and the frequency is degenerate. In Example 4.6, rotating a shear mode of the triangle by $120°$ gives a different shear mode at the same frequency, and that is the whole reason the pair is degenerate. Choosing a basis of the eigenspace adapted to the symmetry (eigenvectors of $S$ within it) is what chemists do when they label molecular vibrations by symmetry species, and it is the same move as CM-4's cyclic coordinates: find the operator that commutes with the dynamics and diagonalize it first. In QM-5 the identical argument, with $[H, S] = 0$, explains why the energy levels of a symmetric system come in degenerate multiplets.

### 3.7 Damping and driving

Real structures dissipate. The standard linear model adds a damping matrix $C$ and an external force $F(t)$:
$$M\ddot q + C\dot q + Kq = F(t). \tag{3.14}$$
In general $C$ wrecks the decoupling: $A^\top CA$ need not be diagonal, and the modes of the damped system are complex and must be found from a $2n$-dimensional first-order problem (set $\mathbf{u} = (q, \dot q)$ and solve $\dot{\mathbf{u}} = A\mathbf{u}$ as in LA-7 §3.8 and Example 4.5). Engineers almost always avoid this with an assumption.

**Theorem 3.6 (Rayleigh damping decouples).** If $C = \alpha M + \beta K$ for scalars $\alpha, \beta \ge 0$, then in normal coordinates (3.14) becomes $n$ independent damped driven oscillators
$$\ddot\xi_j + 2\zeta_j\omega_j\dot\xi_j + \omega_j^2\xi_j = f_j(t), \qquad 2\zeta_j\omega_j = \alpha + \beta\omega_j^2, \qquad f = A^\top F. \tag{3.15}$$

▶ **Try it first.** Substitute $q = A\xi$ into (3.14), multiply on the left by $A^\top$, and ask which of the three resulting matrices (3.11) already makes diagonal. What must $C$ look like for the third one to be diagonal too?

*Proof.* Substitute $q = A\xi$ into (3.14) and multiply on the left by $A^\top$: $A^\top MA\ddot\xi + A^\top CA\dot\xi + A^\top KA\xi = A^\top F$. By (3.11) the first and third matrices are $I$ and $\Omega^2$. The middle one is $A^\top(\alpha M + \beta K)A = \alpha I + \beta\Omega^2$, diagonal with entries $\alpha + \beta\omega_j^2$. ∎

■ **What just happened.** The modes diagonalize $M$ and $K$ and nothing else, so damping decouples exactly when $A^\top CA$ is diagonal, and a combination of $M$ and $K$ is the simplest $C$ for which that is automatic. Real damping (joints, friction, material loss) is not of this form; engineers assume it anyway because it keeps $n$ scalar equations instead of $n$ coupled ones.

The modal damping ratio is $\zeta_j = \alpha/(2\omega_j) + \beta\omega_j/2$: mass-proportional damping ($\alpha$) damps low modes most, stiffness-proportional damping ($\beta$) damps high modes most. Two parameters cannot match measured damping in more than two modes, but the decoupling is worth more than the accuracy, and in practice one fits $\alpha, \beta$ to two measured modes. (Exercise C3 generalizes to the Caughey series $C = M\sum_lc_l(M^{-1}K)^l$, which also decouples and has as many parameters as you like.)

**Frequency response.** Drive with $F(t) = \operatorname{Re}[\mathbf{F}e^{i\omega t}]$ and look for the steady state $\xi_j = \operatorname{Re}[\hat\xi_je^{i\omega t}]$. From (3.15), $(\omega_j^2 - \omega^2 + 2i\zeta_j\omega_j\omega)\hat\xi_j = \mathbf{a}_j^\top\mathbf{F}$, so
$$\hat q(\omega) = \sum_{j=1}^n\frac{\mathbf{a}_j\,(\mathbf{a}_j^\top\mathbf{F})}{\omega_j^2 - \omega^2 + 2i\zeta_j\omega_j\omega}. \tag{3.16}$$
Each term is a single-oscillator resonance: the modal amplitude $|\hat\xi_j|$ peaks near $\omega = \omega_j$ at height $\approx |\mathbf{a}_j^\top\mathbf{F}|/(2\zeta_j\omega_j^2)$, and coordinate $i$ sees that peak multiplied by $|a_{ij}|$, the $i$-th entry of $\mathbf{a}_j$. A sweep of the drive frequency shows one peak per mode (the project does this for a 3-DOF chain); this is how modal frequencies and damping ratios are measured on a real structure with a shaker and accelerometers. Two readable features of (3.16): a mode is not excited by a force applied where its shape has a node ($\mathbf{a}_j^\top\mathbf{F} = 0$), and it is not seen by a sensor placed at a node of its shape (the $j$-th term has no component there).

**Modal superposition as the engineering method.** A finite-element model replaces a continuous structure by $n \sim 10^5$–$10^7$ nodal displacements, assembles $K = B^\top CB$ and $M$ from element matrices, and solves (3.5) for the lowest $r \ll n$ modes with a Lanczos or subspace iteration (LA-13; the Cholesky reduction of Theorem 3.4 is inside). Every later analysis is done with $q \approx \sum_{j\le r}\mathbf{a}_j\xi_j$: $r$ scalar equations (3.15) instead of $n$ coupled ones, with Rayleigh or modal damping assumed so that they decouple. Truncation to $r$ modes is accurate when the loads have no frequency content near the discarded modes, which is why the first question in a vibration analysis is always "what is the frequency range of the excitation?" The exaggerated shapes you see in a modal-analysis animation are the $\mathbf{a}_j$.

### 3.8 The quantum preview

Each normal coordinate obeys $\ddot\xi_j = -\omega_j^2\xi_j$, with Hamiltonian $\frac12(p_j^2 + \omega_j^2\xi_j^2)$: the harmonic oscillator of QM-3. Quantizing one mode gives energy levels $\hbar\omega_j(n_j + \frac12)$; quantizing all of them gives $E = \sum_j\hbar\omega_j(n_j + \frac12)$, and a quantum of mode $j$ in a crystal is a **phonon**. The Debye model (SM-5) takes the dispersion relation $\omega(\kappa)$ of the loaded string in its long-wavelength form $\omega = c\kappa$, counts the modes of a block of solid, and gets the $T^3$ heat capacity; the finite number of atoms sets the maximum frequency, exactly as $N$ masses give $N$ modes in Example 4.5. The continuum limit of Example 4.5 is a field, a displacement $y(x, t)$ at every point, and it is a system of infinitely many oscillators, one per wavenumber. Quantizing those is quantum field theory, and the photon is the phonon of the electromagnetic field.

## 4. Worked examples

Systems (1)–(6) of the module are the examples. Every number below is reproduced by `project.py`.

### Example 4.1: Two equal pendulums coupled by a spring, with beats

$m = 1$ kg, $\ell = 1$ m, $g = 9.81$ m/s², $k = 1$ N/m; the spring joins the bobs and its natural length equals the pivot separation $D$.

*Lagrangian.* Bob $i$ is at $(\pm D/2 + \ell\sin\theta_i, -\ell\cos\theta_i)$, so $T = \frac12m\ell^2(\dot\theta_1^2 + \dot\theta_2^2)$ exactly, and $M = m\ell^2I$. Gravity: $V_g = mg\ell(1 - \cos\theta_1) + mg\ell(1 - \cos\theta_2) = \frac12mg\ell(\theta_1^2 + \theta_2^2) + O(\theta^4)$. Spring: the bobs' horizontal separation is $D + \ell(\sin\theta_2 - \sin\theta_1)$ and their vertical offset $\ell(\cos\theta_1 - \cos\theta_2)$ is $O(\theta^2)$, so the extension is $\ell(\theta_2 - \theta_1) + O(\theta^2)$ and $V_s = \frac12k\ell^2(\theta_2 - \theta_1)^2 + O(\theta^3)$. The Hessian (3.1):
$$K = \begin{pmatrix}mg\ell + k\ell^2 & -k\ell^2\\-k\ell^2 & mg\ell + k\ell^2\end{pmatrix} = \begin{pmatrix}10.81 & -1\\-1 & 10.81\end{pmatrix}\ \mathrm{J}, \qquad M = \begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}\ \mathrm{kg\,m^2}.$$

*Modes.* (3.6): $(mg\ell + k\ell^2 - m\ell^2\omega^2)^2 - k^2\ell^4 = 0$, so $m\ell^2\omega^2 = mg\ell + k\ell^2 \mp k\ell^2$:
$$\omega_1^2 = \frac{g}{\ell} = 9.81\ \mathrm{s^{-2}},\ \ \omega_1 = 3.1321\ \mathrm{s^{-1}}; \qquad \omega_2^2 = \frac{g}{\ell} + \frac{2k}{m} = 11.81\ \mathrm{s^{-2}},\ \ \omega_2 = 3.4366\ \mathrm{s^{-1}}.$$
Shapes from $(K - \omega^2M)\mathbf{a} = 0$: for $\omega_1^2$ the first row reads $k\ell^2a_1 - k\ell^2a_2 = 0$, so $\mathbf{a}_1 \propto (1, 1)^\top$; for $\omega_2^2$ it reads $-k\ell^2a_1 - k\ell^2a_2 = 0$, so $\mathbf{a}_2 \propto (1, -1)^\top$. The §2 guesses, confirmed. $M$-normalize: $(1,1)M(1,1)^\top = 2m\ell^2 = 2$, so $\mathbf{a}_1 = (1,1)^\top/\sqrt2$, $\mathbf{a}_2 = (1,-1)^\top/\sqrt2$ (in units $\mathrm{kg^{-1/2}m^{-1}}$). Because $M$ is numerically the identity here, $M$-orthonormal and orthonormal coincide. In Examples 4.2–4.4 they do not; in 4.5 and 4.6, where $M = mI$ again, they differ only by the factor $\sqrt m$.

*Fitting initial conditions, and beats.* Release pendulum 1 from $\theta_0$ with pendulum 2 at rest at $0$: $q_0 = (\theta_0, 0)^\top$, $v_0 = 0$. The projections in (3.13) are $\langle\mathbf{a}_1, q_0\rangle_M = \mathbf{a}_1^\top Mq_0 = m\ell^2\theta_0/\sqrt2 = \theta_0/\sqrt2$ and $\langle\mathbf{a}_2, q_0\rangle_M = \theta_0/\sqrt2$, so
$$q(t) = \frac{\theta_0}{\sqrt2}\frac{(1,1)^\top}{\sqrt2}\cos\omega_1t + \frac{\theta_0}{\sqrt2}\frac{(1,-1)^\top}{\sqrt2}\cos\omega_2t,$$
$$\theta_1(t) = \frac{\theta_0}{2}(\cos\omega_1t + \cos\omega_2t) = \theta_0\cos(\varepsilon t)\cos(\bar\omega t), \qquad \theta_2(t) = \frac{\theta_0}{2}(\cos\omega_1t - \cos\omega_2t) = \theta_0\sin(\varepsilon t)\sin(\bar\omega t), \tag{4.1}$$
with $\bar\omega = \frac{\omega_1 + \omega_2}{2} = 3.2843\ \mathrm{s^{-1}}$ and $\varepsilon = \frac{\omega_2 - \omega_1}{2} = 0.15224\ \mathrm{s^{-1}}$. The identities used: $\cos a + \cos b = 2\cos\frac{a+b}{2}\cos\frac{a-b}{2}$ and $\cos a - \cos b = -2\sin\frac{a+b}{2}\sin\frac{a-b}{2}$, with $a = \omega_1t$, $b = \omega_2t$, $\frac{a-b}{2} = -\varepsilon t$ (the minus sign and the odd sine cancel). Pendulum 1 swings at the carrier frequency $\bar\omega$ inside the slow envelope $\theta_0|\cos\varepsilon t|$; pendulum 2 inside $\theta_0|\sin\varepsilon t|$. All the energy has moved to pendulum 2 when $\varepsilon t = \pi/2$, and the envelope repeats when $\varepsilon t = \pi$:
$$T_{\text{transfer}} = \frac{\pi}{\omega_2 - \omega_1} = \frac{\pi}{0.30448} = 10.318\ \mathrm{s}, \qquad T_{\text{beat}} = \frac{2\pi}{\omega_2 - \omega_1} = 20.636\ \mathrm{s}. \tag{4.2}$$
For weak coupling, $\omega_2 = \omega_1\sqrt{1 + 2k/(m\omega_1^2)} \approx \omega_1 + k/(m\omega_1)$, so $T_{\text{beat}} \approx 2\pi m\omega_1/k = 19.68$ s; the 4.6% gap is the next term of the square root. Two pencils and a rubber band give a beat period you can time with a watch.

The project integrates the *exact* nonlinear equations of this system (full $\sin\theta$ gravity, exact spring geometry) and overlays (4.1). At $\theta_0 = 0.05$ rad they agree to 2% of $\theta_0$ over two beat periods. At $\theta_0 = 1$ rad the overlay is meaningless. Section 8 explains why the failure comes earlier than "$\theta$ small" suggests.

### Example 4.2: Two masses, three springs, unequal masses — orthogonality done wrong and right

Wall–$k$–$m_1$–$k$–$m_2$–$k$–wall with $m_1 = 1$, $m_2 = 2$, $k = 1$; displacements $x_1, x_2$. $V = \frac12x_1^2 + \frac12(x_2 - x_1)^2 + \frac12x_2^2$, so
$$M = \begin{pmatrix}1 & 0\\0 & 2\end{pmatrix}, \qquad K = \begin{pmatrix}2 & -1\\-1 & 2\end{pmatrix}.$$
This is LA-7's Example 4.3. LA-7 found the numbers below and noticed the $M$-orthogonality; here it is an instance of Theorem 3.2. (With $m_1 = m_2 = m$ the frequencies would be $\omega^2 = k/m$ and $3k/m$ and $M$-orthogonality would reduce to ordinary orthogonality; the unequal masses are the point.)

(3.6) with $\lambda = \omega^2$: $\det\begin{pmatrix}2 - \lambda & -1\\-1 & 2 - 2\lambda\end{pmatrix} = (2 - \lambda)(2 - 2\lambda) - 1 = 2\lambda^2 - 6\lambda + 3 = 0$, so
$$\lambda_{1,2} = \frac{3 \mp\sqrt3}{2}: \qquad \omega_1^2 = 0.63397,\ \omega_1 = 0.79623; \qquad \omega_2^2 = 2.36603,\ \omega_2 = 1.53819.$$
Shapes from the first row, $(2 - \lambda)a_1 - a_2 = 0$, so $a_2 = (2 - \lambda)a_1$:
$$\mathbf{a}_1 \propto \begin{pmatrix}1\\ (1+\sqrt3)/2\end{pmatrix} = \begin{pmatrix}1\\1.36603\end{pmatrix}, \qquad \mathbf{a}_2 \propto \begin{pmatrix}1\\(1-\sqrt3)/2\end{pmatrix} = \begin{pmatrix}1\\-0.36603\end{pmatrix}.$$
Slow mode: both masses move together, the heavy one further. Fast mode: opposed, the heavy one barely moving.

*Ordinary dot product:* $1\cdot1 + 1.36603\cdot(-0.36603) = 1 - 0.5 = 0.5 \ne 0$. **Not orthogonal.**
*$M$-inner product (3.8):* $1\cdot1\cdot1 + 2\cdot1.36603\cdot(-0.36603) = 1 - 1 = 0$. **Orthogonal.** (Exactly: $(1 + \sqrt3)(1 - \sqrt3)/4 = -1/2$, times $m_2 = 2$ gives $-1$.)

Theorem 3.2, in two lines of arithmetic. Do this one by hand; it is the computation that makes the module stick. For the record, $M^{-1}K = \begin{pmatrix}2 & -1\\-1/2 & 1\end{pmatrix}$ is not symmetric and has these same non-orthogonal eigenvectors with the same real eigenvalues; Exercise B3 shows why its eigenvalues are real anyway (it is similar to the symmetric $\tilde K$ of (3.9)). Normalizing, $\mathbf{a}_1^\top M\mathbf{a}_1 = 1 + 2(1.36603)^2 = 4.73205$ and $\mathbf{a}_2^\top M\mathbf{a}_2 = 1 + 2(0.36603)^2 = 1.26795$, so
$$\mathbf{a}_1 = \begin{pmatrix}0.45970\\0.62796\end{pmatrix}, \qquad \mathbf{a}_2 = \begin{pmatrix}0.88807\\-0.32506\end{pmatrix},$$
and the project checks $A^\top MA = I$, $A^\top KA = \operatorname{diag}(0.63397, 2.36603)$.

### Example 4.3: The linear triatomic molecule (CO₂ model)

Three masses on a line, $m$–$M_c$–$m$ (O–C–O), joined by two identical springs $k$; displacements $x_1, x_2, x_3$ along the line. $V = \frac12k(x_2 - x_1)^2 + \frac12k(x_3 - x_2)^2$, so
$$K = k\begin{pmatrix}1 & -1 & 0\\-1 & 2 & -1\\0 & -1 & 1\end{pmatrix}, \qquad M = \begin{pmatrix}m & 0 & 0\\0 & M_c & 0\\0 & 0 & m\end{pmatrix}.$$
$K = B^\top CB$ with $B = \begin{pmatrix}-1 & 1 & 0\\0 & -1 & 1\end{pmatrix}$ (two springs, two rows), and $N(B) = \operatorname{span}\{(1,1,1)^\top\}$: the translation. So $K$ is only PSD and there is exactly one zero mode.

Use the reflection through the carbon (§3.6): the two oxygens trade places and every displacement reverses direction, $S(x_1, x_2, x_3)^\top = (-x_3, -x_2, -x_1)^\top$. It satisfies $S^\top MS = M$, $S^\top KS = K$ and $S^2 = I$. Since the three frequencies turn out distinct, each mode is symmetric ($S\mathbf{a} = \mathbf{a}$) or antisymmetric ($S\mathbf{a} = -\mathbf{a}$) under it, and that is where the names below come from. (The bare relabeling $x_1 \leftrightarrow x_3$ is also a symmetry, but it is $-S$, so it swaps the two labels.) Guess, then verify in (3.5):

- **Zero mode** $\mathbf{a}_0 = (1, 1, 1)^\top$: $K\mathbf{a}_0 = 0$, $\omega_0 = 0$. Rigid translation.
- **Symmetric stretch** $\mathbf{a}_s = (1, 0, -1)^\top$ (oxygens out together, carbon still): $K\mathbf{a}_s = k(1, 0, -1)^\top = k\,\mathbf{a}_s$ and $M\mathbf{a}_s = m\,\mathbf{a}_s$, so $\omega_s^2 = k/m$.
- **Antisymmetric stretch** $\mathbf{a}_a = (1, -2m/M_c, 1)^\top$ (oxygens one way, carbon the other, chosen so that $\sum_im_ia_i = 0$): $K\mathbf{a}_a = k\big(1 + \tfrac{2m}{M_c},\ -2 - \tfrac{4m}{M_c},\ 1 + \tfrac{2m}{M_c}\big)^\top = k\big(1 + \tfrac{2m}{M_c}\big)(1, -2, 1)^\top$, and $M\mathbf{a}_a = (m, -2m, m)^\top = m(1, -2, 1)^\top$. So (3.5) holds with
$$\omega_a^2 = \frac{k}{m}\Big(1 + \frac{2m}{M_c}\Big) = \frac{k}{m} + \frac{2k}{M_c}, \qquad \frac{\omega_a}{\omega_s} = \sqrt{1 + \frac{2m}{M_c}}. \tag{4.3}$$

Three mutually $M$-orthogonal vectors (check $\mathbf{a}_0^\top M\mathbf{a}_a = m - 2m + m = 0$: the center of mass does not move) in a 3-dimensional space: complete by Theorem 3.4, so these are all the modes. The antisymmetric stretch is faster because the light carbon does most of the moving: its Rayleigh quotient has a large numerator and a small denominator.

*Numbers* with atomic masses $m = 16$, $M_c = 12$ and $k = 1$: $\omega_s^2 = 1/16 = 0.0625$, $\omega_s = 0.25$; $\omega_a^2 = 1/16 + 2/12 = 0.229167$, $\omega_a = 0.47871$; $\omega_a/\omega_s = \sqrt{11/3} = 1.9149$. $\mathbf{a}_a = (1, -8/3, 1)^\top$, $\mathbf{a}_a^\top M\mathbf{a}_a = 16 + 12\cdot\frac{64}{9} + 16 = 117.33$, so the normalized shape is $(1, -2.6667, 1)^\top/10.832$. The measured CO₂ stretches are about $1388\ \mathrm{cm^{-1}}$ (symmetric; Raman-active, infrared-dark because the dipole does not change; a Fermi resonance with the bend overtone splits this line into $1285$ and $1388\ \mathrm{cm^{-1}}$, and $1388$ is the value usually quoted) and $2349\ \mathrm{cm^{-1}}$ (antisymmetric; the strong infrared band at $4.3\ \mu\mathrm{m}$), ratio $1.69$. Two independent springs get the ordering and the symmetry classification right and the ratio 13% high; the real potential also couples the two bonds directly. The bending modes are the degenerate pair of transverse motions that this one-dimensional model cannot see. At $667\ \mathrm{cm^{-1}}$ ($15\ \mu\mathrm{m}$, near the peak of Earth's thermal emission) the bend, not the antisymmetric stretch, is the band that makes CO₂ a greenhouse gas.

### Example 4.4: The double pendulum, linearized

From CM-3, with equal masses $m$ and lengths $\ell$ and angles $\theta_1$ (upper), $\theta_2$ (lower) from the vertical: $\mathbf{r}_1 = \ell(\sin\theta_1, -\cos\theta_1)$, $\mathbf{r}_2 = \mathbf{r}_1 + \ell(\sin\theta_2, -\cos\theta_2)$. Then $|\dot{\mathbf{r}}_1|^2 = \ell^2\dot\theta_1^2$ and $|\dot{\mathbf{r}}_2|^2 = \ell^2\dot\theta_1^2 + \ell^2\dot\theta_2^2 + 2\ell^2\dot\theta_1\dot\theta_2(\cos\theta_1\cos\theta_2 + \sin\theta_1\sin\theta_2)$, so
$$T = \frac12m\ell^2\big[2\dot\theta_1^2 + \dot\theta_2^2 + 2\dot\theta_1\dot\theta_2\cos(\theta_1 - \theta_2)\big], \qquad V = -mg\ell(2\cos\theta_1 + \cos\theta_2).$$
Second order: $\cos(\theta_1 - \theta_2) \to 1$ in $T$ (the correction is $O(\theta^2\dot\theta^2)$) and $-\cos\theta \to -1 + \theta^2/2$ in $V$. Reading off $M_{ij} = \partial^2T/\partial\dot\theta_i\partial\dot\theta_j$ and $K_{ij} = \partial^2V/\partial\theta_i\partial\theta_j$ at $0$:
$$M = m\ell^2\begin{pmatrix}2 & 1\\1 & 1\end{pmatrix}, \qquad K = mg\ell\begin{pmatrix}2 & 0\\0 & 1\end{pmatrix}.$$
$M$ is not diagonal: the coordinates are coupled through the kinetic energy, not the potential. This is the case the ordinary spectral theorem cannot touch.

Set $\lambda = \omega^2\ell/g$. (3.6): $\det\begin{pmatrix}2 - 2\lambda & -\lambda\\-\lambda & 1 - \lambda\end{pmatrix} = (2 - 2\lambda)(1 - \lambda) - \lambda^2 = \lambda^2 - 4\lambda + 2 = 0$,
$$\lambda = 2 \mp\sqrt2: \qquad \omega_1^2 = (2 - \sqrt2)\frac{g}{\ell} = 0.58579\,\frac{g}{\ell}, \qquad \omega_2^2 = (2 + \sqrt2)\frac{g}{\ell} = 3.41421\,\frac{g}{\ell}. \tag{4.4}$$
With $g/\ell = 9.81\ \mathrm{s^{-2}}$: $\omega_1 = 2.3972\ \mathrm{s^{-1}}$, $\omega_2 = 5.7874\ \mathrm{s^{-1}}$, and $\omega_2/\omega_1 = \sqrt{(2 + \sqrt2)/(2 - \sqrt2)} = 1 + \sqrt2 = 2.4142$.

Shapes from the first row, $(2 - 2\lambda)a_1 = \lambda a_2$, i.e. $a_2/a_1 = 2(1 - \lambda)/\lambda$. For $\lambda_1 = 2 - \sqrt2$: $a_2/a_1 = \dfrac{2(\sqrt2 - 1)}{2 - \sqrt2} = \dfrac{2(\sqrt2 - 1)}{\sqrt2(\sqrt2 - 1)} = \sqrt2$. For $\lambda_2 = 2 + \sqrt2$: $a_2/a_1 = \dfrac{-2(1 + \sqrt2)}{2 + \sqrt2} = \dfrac{-2(1+\sqrt2)}{\sqrt2(\sqrt2 + 1)} = -\sqrt2$.
$$\mathbf{a}_1 \propto \begin{pmatrix}1\\\sqrt2\end{pmatrix}\ \text{(in phase; the lower bob swings }\sqrt2\text{ times further)}, \qquad \mathbf{a}_2 \propto \begin{pmatrix}1\\-\sqrt2\end{pmatrix}\ \text{(opposed)}.$$
Check (3.7): $\mathbf{a}_1^\top M\mathbf{a}_2 = m\ell^2\big[1\cdot(2\cdot1 + 1\cdot(-\sqrt2)) + \sqrt2\,(1\cdot1 + 1\cdot(-\sqrt2))\big] = m\ell^2[2 - \sqrt2 + \sqrt2 - 2] = 0$. ✓ Ordinary dot product: $1 - 2 = -1 \ne 0$. Normalization: $\mathbf{a}_1^\top M\mathbf{a}_1 = m\ell^2(2 + 2\sqrt2 + 2) = (4 + 2\sqrt2)m\ell^2 = 6.8284\,m\ell^2$ and $\mathbf{a}_2^\top M\mathbf{a}_2 = (4 - 2\sqrt2)m\ell^2 = 1.1716\,m\ell^2$. (A coincidence worth noticing, not a connection: $\lambda^2 - 4\lambda + 2$ is also satisfied by $\omega^2/\omega_0^2$ for the $j = 1$ and $j = 3$ modes of the $N = 3$ loaded string below.)

### Example 4.5: $N$ masses on a string, and the continuum limit

$N$ equal masses $m$ at spacing $d$ on a massless string under tension $\tau$, ends fixed at $x_0 = 0$ and $x_{N+1} = L = (N+1)d$. Transverse displacements $y_1, \dots, y_N$, with $y_0 = y_{N+1} = 0$.

**Derivation 5 (the equations).** ▶ **Try it first.** Newton's law on mass $k$, keeping only first order in the slopes.

1. The segment from mass $k$ to mass $k+1$ has slope $(y_{k+1} - y_k)/d$, and the transverse component of its tension on mass $k$ is $\tau\sin\phi \approx \tau(y_{k+1} - y_k)/d$ (the tension's magnitude changes only at second order). The segment on the left contributes $\tau(y_{k-1} - y_k)/d$.
2. $m\ddot y_k = \dfrac{\tau}{d}\,(y_{k+1} - 2y_k + y_{k-1})$, $k = 1, \dots, N$, with $y_0 = y_{N+1} = 0$.
3. In the form (3.3): $M = mI$ and $K = \dfrac{\tau}{d}\begin{pmatrix}2 & -1 & & \\ -1 & 2 & -1 & \\ & \ddots & \ddots & \ddots\\ & & -1 & 2\end{pmatrix}$, tridiagonal. ($V = \frac{\tau}{2d}\sum_{k=0}^N(y_{k+1} - y_k)^2$, the tension times the extra length of string to second order; this is $B^\top CB$ with $B$ the $(N+1)\times N$ difference matrix.) Write $\omega_0^2 = \tau/(md)$. The same matrices describe $N$ masses joined by $N+1$ springs of stiffness $\tau/d$ moving longitudinally, so everything below applies to a one-dimensional crystal. ∎

**Derivation 6 (the modes).** ▶ **Try it first.** Try the shape $a_k = \sin(k\phi)$ for some angle $\phi$. Use $\sin((k\pm1)\phi) = \sin k\phi\cos\phi \pm \cos k\phi\sin\phi$.

1. $(K\mathbf{a})_k = \dfrac{\tau}{d}(2a_k - a_{k+1} - a_{k-1})$. With $a_k = \sin k\phi$: $a_{k+1} + a_{k-1} = 2\sin k\phi\cos\phi$ (the cosine terms cancel), so $(K\mathbf{a})_k = \dfrac{\tau}{d}(2 - 2\cos\phi)\sin k\phi = \dfrac{\tau}{d}\,4\sin^2\!\frac{\phi}{2}\,\sin k\phi$, using $1 - \cos\phi = 2\sin^2(\phi/2)$.
2. So $K\mathbf{a} = m\cdot4\omega_0^2\sin^2\!\frac{\phi}{2}\,\mathbf{a} = \omega^2M\mathbf{a}$ with $\omega^2 = 4\omega_0^2\sin^2(\phi/2)$: the recurrence is satisfied in the interior for *every* $\phi$. The boundary conditions select $\phi$.
3. $a_0 = \sin 0 = 0$ holds automatically. $a_{N+1} = \sin((N+1)\phi) = 0$ requires $(N+1)\phi = j\pi$ for an integer $j$: $\phi_j = j\pi/(N+1)$.
4. $j = 1, \dots, N$ give $N$ values $\phi_j \in (0, \pi)$, hence $N$ distinct frequencies ($\sin(\phi/2)$ is increasing on $(0,\pi)$). $j = 0$ and $j = N+1$ give $\mathbf{a} = 0$; $j \to 2(N+1) - j$ flips the sign of the shape, so other $j$ repeat these. By Theorem 3.4 there are exactly $N$ modes, so this is all of them:
$$a^{(j)}_k = \sin\frac{jk\pi}{N+1}, \qquad \omega_j = 2\omega_0\sin\frac{j\pi}{2(N+1)}, \qquad j = 1, \dots, N. \tag{4.5}$$ ∎

■ **What just happened.** Mode $j$ has $j - 1$ interior nodes: it is the discrete version of the $j$-th harmonic of a string. The frequencies are *not* in the ratio $1 : 2 : 3$; they bend over because of the sine. For $N = 5$: $\omega_j/\omega_0 = 2\sin(j\pi/12) = 0.5176,\ 1,\ 1.4142,\ 1.7321,\ 1.9319$. For $N = 3$: $\omega_j^2/\omega_0^2 = 2 - \sqrt2,\ 2,\ 2 + \sqrt2$. And Theorem 3.2 with $M = mI$ says $\sum_{k=1}^N\sin\frac{jk\pi}{N+1}\sin\frac{j'k\pi}{N+1} = 0$ for $j \ne j'$: a discrete orthogonality of sines, obtained with no trigonometry at all. (Exercise B5 computes the normalization, $(N+1)/2$.)

**Derivation 7 (the continuum limit).** ▶ **Try it first.** Fix the length $L$ and the mass per unit length $\mu = m/d$, and let $N \to \infty$. Write $y_k(t) = y(x_k, t)$ with $x_k = kd$ and Taylor-expand $y(x \pm d)$ in step 2 of Derivation 5.

1. As $N \to \infty$ with $L$ and $\mu$ fixed, $d = L/(N+1) \to 0$ and $m = \mu d \to 0$.
2. Suppose the displacements are samples of a smooth function, $y_k(t) = y(x_k, t)$. Taylor: $y(x \pm d) = y \pm d\,y_x + \frac12d^2y_{xx} \pm \frac16d^3y_{xxx} + O(d^4)$, so $y_{k+1} - 2y_k + y_{k-1} = d^2\,y_{xx}(x_k) + O(d^4)$. The second difference is $d^2$ times the second derivative.
3. Substitute into $m\ddot y_k = \frac{\tau}{d}(y_{k+1} - 2y_k + y_{k-1})$ with $m = \mu d$: $\mu d\,y_{tt} = \frac{\tau}{d}\big(d^2y_{xx} + O(d^4)\big)$, so $y_{tt} = \frac{\tau}{\mu}y_{xx} + O(d^2)$. Letting $d \to 0$:
$$\frac{\partial^2y}{\partial t^2} = c^2\frac{\partial^2y}{\partial x^2}, \qquad c = \sqrt{\tau/\mu}, \qquad y(0, t) = y(L, t) = 0. \tag{4.6}$$
The wave equation. A system of ODEs has become a PDE: one unknown $y(x, t)$ at each of a continuum of points.
4. The mode shapes (4.5) are $a^{(j)}_k = \sin\frac{jk\pi}{N+1} = \sin\frac{j\pi x_k}{L}$: the function $\sin(j\pi x/L)$ sampled at the masses. In the limit the modes are the standing waves $\sin(j\pi x/L)$, the eigenfunctions of $d^2/dx^2$ on $[0, L]$ with zero boundary values, and the discrete orthogonality noted after Derivation 6 becomes $\int_0^L\sin\frac{j\pi x}{L}\sin\frac{j'\pi x}{L}\,dx = \frac L2\delta_{jj'}$. This is the sine basis of F2, and Bernoulli's "every motion is a sum of modes" becomes "every initial shape of the string is a Fourier sine series."
5. Frequencies. Write $\kappa_j = j\pi/L$ for the wavenumber of $\sin\kappa_jx$. Then $\phi_j = j\pi/(N+1) = \kappa_jd$, and $\omega_0 = \sqrt{\tau/(md)} = \sqrt{\tau/(\mu d^2)} = c/d$, so (4.5) reads
$$\omega(\kappa) = \frac{2c}{d}\sin\frac{\kappa d}{2}. \tag{4.7}$$
This is the **dispersion relation** of the chain: the frequency as a function of the wavenumber.
6. Long wavelengths, $\kappa d \ll 1$ (wavelength $2\pi/\kappa$ much larger than the spacing): $\sin(\kappa d/2) \approx \kappa d/2$ and
$$\omega \approx c\kappa.$$
Linear: every wavelength travels at the same speed $c$, which is what the wave equation (4.6) says ($y = f(x - ct)$ solves it for any $f$). The $j$-th mode has $\omega_j \approx c\,j\pi/L$, frequency $f_j = jc/(2L)$: the harmonic series of a string.
7. Short wavelengths: $\omega(\kappa)$ bends over and reaches its maximum $\omega_{\max} = 2c/d = 2\omega_0$ at $\kappa = \pi/d$, where the wavelength is $2d$ and neighboring masses move oppositely. That is the fastest pattern a lattice can hold; there is no mode above it, and there are only $N$ modes in total. The group velocity $d\omega/d\kappa = c\cos(\kappa d/2)$ falls to zero there: a disturbance at that wavelength does not propagate. The chain is *dispersive* (different wavelengths, different speeds); the continuum it approximates is not. ∎

■ **What just happened.** The discrete chain and the continuous string agree for wavelengths long compared with the spacing, and disagree at wavelengths comparable to it, where the lattice runs out of modes. This is exactly the structure of the Debye model (SM-5): treat the crystal as a continuum with $\omega = c\kappa$, but stop counting modes at a cutoff chosen so that there are $3N$ of them. It is also your first look at how a PDE emerges from a large system of ODEs, and it is worth reading twice.

### Example 4.6: Three masses on an equilateral triangle (degeneracy)

Three equal masses $m$ at the vertices of an equilateral triangle of circumradius $R = 1$, positions $\mathbf{r}_i = (\cos\theta_i, \sin\theta_i)$ with $\theta_i = 90°, 210°, 330°$, joined by three identical springs $k$ along the sides (natural length $\sqrt3$). Each mass moves in the plane: $q = (\mathbf{u}_1, \mathbf{u}_2, \mathbf{u}_3) \in \mathbb{R}^6$ and $M = mI_6$.

*Stiffness.* To first order the extension of the spring from $i$ to $j$ is the projection of the relative displacement onto the spring's direction: with $\hat{\mathbf{n}}_{ij} = (\mathbf{r}_j - \mathbf{r}_i)/|\mathbf{r}_j - \mathbf{r}_i|$,
$$e_{ij} = |\mathbf{r}_j + \mathbf{u}_j - \mathbf{r}_i - \mathbf{u}_i| - |\mathbf{r}_j - \mathbf{r}_i| = \hat{\mathbf{n}}_{ij}\cdot(\mathbf{u}_j - \mathbf{u}_i) + O(u^2)$$
(the first-order change in a vector's length is the change projected onto its direction). So $e = Bq$ with $B$ the $3 \times 6$ matrix whose row for spring $ij$ holds $-\hat{\mathbf{n}}_{ij}$ in slot $i$ and $+\hat{\mathbf{n}}_{ij}$ in slot $j$, and $K = kB^\top B$. The project assembles this numerically.

*Zero modes.* $B$ has 3 rows, so $\dim N(B) \ge 3$: at least three zero modes. Two translations, $(1,0,1,0,1,0)^\top$ and $(0,1,0,1,0,1)^\top$ (every $\mathbf{u}_j - \mathbf{u}_i = 0$), and the rigid rotation $\mathbf{u}_i = \hat{\mathbf{z}}\times\mathbf{r}_i = (-\sin\theta_i, \cos\theta_i)$, for which $(\mathbf{u}_j - \mathbf{u}_i)\cdot(\mathbf{r}_j - \mathbf{r}_i) = [\hat{\mathbf{z}}\times(\mathbf{r}_j - \mathbf{r}_i)]\cdot(\mathbf{r}_j - \mathbf{r}_i) = 0$. A triangle of bars is a rigid structure, so these are the only zero modes ($\operatorname{rank}B = 3$; the project checks it). The eigenvalue $0$ has a 3-dimensional eigenspace; the computer returns some arbitrary $M$-orthonormal basis of it, mixing translation and rotation.

*Breathing mode.* $\mathbf{u}_i = \epsilon\,\mathbf{r}_i$ (every mass straight out). Each spring extends by $e_{ij} = \epsilon\,\hat{\mathbf{n}}_{ij}\cdot(\mathbf{r}_j - \mathbf{r}_i) = \epsilon|\mathbf{r}_j - \mathbf{r}_i| = \epsilon\sqrt3$. The two springs at vertex $i$ pull it toward the other two vertices with force $k\epsilon\sqrt3$ each, along unit vectors at $\pm30°$ from the inward direction $-\mathbf{r}_i$ (the interior angle is $60°$), so the net force is $k\epsilon\sqrt3\cdot2\cos30°\cdot(-\mathbf{r}_i) = -3k\epsilon\,\mathbf{r}_i = -3k\,\mathbf{u}_i$. Thus $m\ddot{\mathbf{u}}_i = -3k\mathbf{u}_i$ and
$$\omega_b^2 = \frac{3k}{m}.$$

*The degenerate pair.* Four modes found; two remain, $M$-orthogonal to the four. The trace of $K$ is the sum of its eigenvalues (LA-7), and $\operatorname{tr}K = k\sum_{\text{springs}}|\mathbf{b}_{ij}|^2 = k\cdot3\cdot2 = 6k$ (each row of $B$ holds two unit vectors). So $0 + 0 + 0 + 3k + \lambda_5 + \lambda_6 = 6k$, i.e. $\lambda_5 + \lambda_6 = 3k$. Now use §3.6 with $S$ the $120°$ rotation of the triangle (which permutes the masses and rotates each displacement; $S^\top S = I$, $S^\top KS = K$). $S$ maps the remaining 2-dimensional subspace to itself. If $\lambda_5 \ne \lambda_6$, each of the two eigenvectors would be mapped to a real multiple of itself; with $S^3 = I$ that multiple is $1$, so the mode would be rotation-invariant. But a rotation-invariant displacement gives every vertex the same radial and the same tangential component, and those two patterns are the breathing mode and the rigid rotation, already used. So no vector in the remaining subspace is $S$-invariant, hence $\lambda_5 = \lambda_6$:
$$\omega_5^2 = \omega_6^2 = \frac{3k}{2m}.$$
These are the two **shear** modes, which deform the triangle into an isosceles shape (one vertex out, the other two in and sideways); the whole 2-dimensional space of them is the eigenspace, and the rotated copy of any shear mode is another. The project computes the spectrum $(0, 0, 0, 1.5, 1.5, 3)\,k/m$ numerically, draws the shapes, and shows that the pair it returns is some rotation of the "textbook" pair. Nothing is wrong; §8.

## 5. Why it matters and where it connects

- **CM-7.** The kinetic energy of a rigid body is $\frac12\boldsymbol\omega^\top I\boldsymbol\omega$ with $I$ the inertia tensor: a symmetric positive definite matrix playing the role of $M$. Principal axes are its eigenvectors; the spectral theorem of LA-8 is why they exist and are orthogonal, exactly as here. Small oscillations of a rigid body about a stable orientation (the spinning top of CM-7, the car body of C2) are this module with that $M$.
- **F2 (Fourier series).** Derivation 7 is where the sine basis comes from, and Theorem 3.2 is where its orthogonality comes from. The Fourier coefficients of an initial shape are the projections $\langle\mathbf{a}_j, q_0\rangle_M$ of (3.13), with $M$ the identity (a sum over masses) becoming the integral $\int_0^L$.
- **QM-3, QM-4.** Each normal coordinate is a harmonic oscillator; quantize it and you get $\hbar\omega_j(n_j + \frac12)$. The algebra of this module, with $(K, M)$ replaced by $(H, \text{identity})$ and real symmetric replaced by Hermitian, is QM-4's: eigenvalues are what you measure, eigenvectors are orthogonal, completeness is why you can expand any state. The generalized problem returns in QM-4 as the Schrödinger equation in a non-orthogonal basis, $H\mathbf{c} = ES\mathbf{c}$ with $S$ the overlap matrix (the Roothaan equations of quantum chemistry are exactly $K\mathbf{a} = \omega^2M\mathbf{a}$).
- **SM-5.** The Debye model is Derivation 7's dispersion relation, counted: $\omega = c\kappa$ up to a cutoff fixed by the number of atoms. Phonons are quantized normal modes; the $T^3$ law comes from counting modes with $\omega < k_BT/\hbar$.
- **CM-4.** Zero modes are Noether's theorem in linearized form (§3.5); degenerate modes are symmetry in linearized form (§3.6). QM-5 uses the same argument for degenerate multiplets.
- **LA-10, LA-13.** Theorem 3.4 is why Cholesky is in the decompositions atlas. Industrial modal analysis computes the lowest $r$ eigenpairs of a sparse $10^6 \times 10^6$ pair $(K, M)$ with Lanczos on the Cholesky-reduced problem; LA-13's Krylov methods are the engine.
- **ML-6.** Spectral clustering solves the generalized eigenproblem $L\mathbf{v} = \lambda D\mathbf{v}$ with $L$ the graph Laplacian (a $B^\top B$ matrix, exactly the stiffness of a spring network on the graph) and $D$ the diagonal degree matrix in the role of $M$. The zero mode is the constant vector; the first nonzero mode (the Fiedler vector) cuts the graph. Same theorems, same code.
- **Engineering.** Modal analysis, modal testing (§3.7), tuned mass dampers (B2(d)), why a helicopter rotor's blade frequencies must avoid multiples of the rotation rate: all are this module. The Tacoma Narrows collapse (1940) is often told as resonance; it was flutter, a torsional mode of the deck made unstable by the wind (an effective negative damping), which needs this module's modes plus an aerodynamic force that depends on the motion.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills: the generalized eigenproblem by hand, $M$-normalization, when $\langle\cdot,\cdot\rangle_M$ fails to be an inner product, a Hessian with a stable and an unstable equilibrium, the $N = 3$ string, and a beat period. B1–B8 are the module: B1 (Rayleigh damping decouples, general damping does not), B2 (double pendulum with unequal masses), B3 (why $M^{-1}K$ has real eigenvalues although it is not symmetric), B4 (modal energies), B5 (discrete sine orthogonality), B6 (zero modes and momentum), B7 (extend `project.py`), B8 (the receptance formula and nodes). B3, B6 and B7 are essential. C1–C3: the 2D lattice and its separable modes, a four-degree-of-freedom car model, and the Caughey series.

## 7. Coding project

`project.py` builds, in order:

1. **`normal_modes(K, M)`** by the Cholesky reduction of Theorem 3.4 and `np.linalg.eigh`, asserted against `scipy.linalg.eigh(K, M)`; $A^\top MA = I$ and $A^\top KA = \Omega^2$ asserted (3.11), on random $(K, M)$ pairs including indefinite $K$.
2. **Examples 4.1, 4.2, 4.4** with every printed number from §4 asserted, the Hessians computed symbolically from the exact $T$ and $V$, and the non-orthogonality of the eigenvectors of $M^{-1}K$ shown explicitly on Example 4.2.
3. **Beats.** The exact nonlinear coupled pendulums integrated with `solve_ivp` and overlaid with (4.1), at $\theta_0 = 0.05$ and $1.0$ rad; agreement asserted small for the first and large for the second.
4. **CO₂** (Example 4.3): frequencies, the zero mode, mode shapes drawn as arrows on the three atoms.
5. **The loaded string** for $N = 5, 20, 100$: numerical frequencies against (4.5) (asserted), mode shapes against sines, the dispersion relation (4.7) against $\omega = c\kappa$, and a grid of snapshots of the first three modes.
6. **The triangle** (Example 4.6): $K = kB^\top B$ assembled from geometry, $\operatorname{rank}B = 3$, spectrum $(0,0,0,1.5,1.5,3)k/m$ asserted, the degenerate pair identified by near-equal eigenvalues, mode shapes plotted.
7. **Frequency response** of a damped, driven 3-DOF chain: (3.16) against a direct solve of $(K - \omega^2M + i\omega C)\hat q = \mathbf{F}$, three resonance peaks located at the $\omega_j$.
8. **`fit_initial_conditions(q0, v0, modes)`** by $M$-projection (3.13), zero modes included, asserted against direct integration of (3.3) for CO₂ with a drifting initial condition.
9. **Exercise checks.** Every number quoted in `solutions.md`, including the reference solution to B7 (the exact double pendulum against its modal solution) and the car model of C2.

**Build it yourself first:** before reading Part 1, write `normal_modes(K, M)` (Cholesky, form $\tilde K$, `eigh`, map back with `solve_triangular`, sort) and `fit_initial_conditions`. The two things you will get wrong the first time are forgetting to $M$-normalize and dividing by $\omega = 0$.

**Your turn** (at the end of the script): a $5\times5$ lattice of masses and its lowest modes, with a degenerate pair split by a 1% heavier bead; Rayleigh damping checked for decoupling on the triangle, and the damped free response by modal superposition; the bounce–pitch–axle car of Exercise C2 (whose ride and wheel-hop modes Part 9 already computes) driven by a road profile, with body response plotted against speed.

## 8. Where you'll struggle

- **$M^{-1}K$ is not the matrix you think it is.** Every student rewrites $K\mathbf{a} = \omega^2M\mathbf{a}$ as $M^{-1}K\mathbf{a} = \omega^2\mathbf{a}$, calls `eig`, and then wonders why the eigenvectors are not orthogonal and why the normalization "doesn't work." $M^{-1}K$ is not symmetric (Example 4.2: $\begin{pmatrix}2&-1\\-1/2&1\end{pmatrix}$). Its eigenvalues are real and its eigenvectors are complete, but they are orthogonal only in the $M$-inner product, and that fact takes about a day to accept. The fix is to stop thinking of $M^{-1}K$ as a matrix and think of it as an operator on the inner product space $(\mathbb{R}^n, \langle\cdot,\cdot\rangle_M)$, where it is self-adjoint. Operationally: never form $M^{-1}K$; use the Cholesky reduction or `eigh(K, M)`, and check $A^\top MA = I$, not $A^\top A = I$.
- **Linearization fails earlier than you expect, and it fails in phase.** "Small angle" is not enough. For beats the relevant comparison is between the anharmonic frequency shift of a pendulum, about $\omega_0\theta_0^2/16$, and the mode splitting $\omega_2 - \omega_1$ that sets the beat. At $\theta_0 = 0.3$ rad the shift is 6% of the splitting; after two beat periods the carrier has drifted by most of a radian and the linear and nonlinear curves are out of phase, though the envelopes still roughly match. Students routinely run the linear model far into this regime and read the output as physics. Always integrate the nonlinear system once and look.
- **Degenerate modes: the computer's eigenvectors are not the textbook's.** For the triangle, `eigh` returns some $M$-orthonormal basis of the 2-dimensional shear eigenspace and some basis of the 3-dimensional zero-mode space, and neither looks like "translation in $x$" or "one vertex out." Nothing is wrong. Any basis of the eigenspace is correct, and floating point decides which one you get. If you want the symmetric ones, project onto the eigenspace and diagonalize the symmetry operator there, or add a tiny symmetry-breaking perturbation and watch the degenerate pair split.
- **The continuum limit deserves slow reading.** Derivation 7 is the first time in this curriculum that a PDE comes out of a system of ODEs. Three separate things are happening: the second difference becomes a second derivative (Taylor), the sampled sine becomes a function (the mode shapes), and the dispersion relation straightens out at long wavelength (the frequencies). Students who read it fast come away believing that the continuum is "more exact," when it is the discrete chain that is the physical model of a crystal and the continuum that is the approximation.
- **Zero modes crash naive code.** Any formula with $\omega_j$ in a denominator ((3.13)'s $\sin\omega_jt/\omega_j$, the damping ratio $\zeta_j = \alpha/(2\omega_j) + \dots$, a frequency response at $\omega = 0$) fails on a zero mode, and `eigh` returns the zero eigenvalue as $\pm10^{-16}$, so `np.sqrt` produces `nan` from the negative one. Threshold explicitly: treat $|\omega^2| < \epsilon\cdot\max\omega^2$ as zero and use $\xi_j(t) = \xi_j(0) + \dot\xi_j(0)t$. The project does this; check that yours does.
- **Conventions.** $\lambda$ is $\omega^2$, not $\omega$; `eigh` returns $\omega^2$ ascending; the mass-normalized $\mathbf{a}_j$ have units $\mathrm{kg^{-1/2}}$ (or $\mathrm{kg^{-1/2}m^{-1}}$ for angles), so normal coordinates have units $\mathrm{kg^{1/2}m}$. Some books normalize $\mathbf{a}_j^\top M\mathbf{a}_j = 1$, some $\max_k|a_k| = 1$, some $\mathbf{a}_j^\top\mathbf{a}_j = 1$; only the first makes (3.11) hold without extra factors.

## 9. Check yourself

1. Starting from $T = \frac12\sum_\alpha m_\alpha|\dot{\mathbf{r}}_\alpha|^2$, what is $M_{ij}$, and why is $M$ positive definite?
2. What is $K$ in terms of $V$? What does $K \succ 0$ tell you about the equilibrium?
3. Substitute $q = \mathbf{a}\cos(\omega t - \delta)$ into $M\ddot q = -Kq$. What equation results? Why is it called *generalized*?
4. Prove in four lines that $\omega^2$ is real.
5. Prove in four lines that $\mathbf{a}_i^\top M\mathbf{a}_j = 0$ for $\omega_i \ne \omega_j$. Where exactly does the proof differ from LA-8's?
6. What are the three axioms of an inner product, and which property of $M$ makes $\langle\mathbf{u},\mathbf{v}\rangle_M$ satisfy each?
7. Write $\tilde K$ in terms of $K$ and the Cholesky factor of $M$. Why is it symmetric? How do its eigenvectors map to mode shapes?
8. What is $A^{-1}$ when $A^\top MA = I$? Give the formula for the normal coordinates in terms of $q$.
9. Write down the general solution of $M\ddot q = -Kq$ fitted to $q(0) = q_0$, $\dot q(0) = v_0$, including a zero mode.
10. For CO₂ with $m = 16$, $M_c = 12$: the three mode shapes and the ratio of the two stretch frequencies.
11. For $N$ beads on a string, what is $\omega_j$? What does it become as $N \to \infty$ at fixed $L$, $\mu$, and for which $j$ is that a good approximation?
12. Why do the triangle's two shear modes have the same frequency? Why does $C = \alpha M + \beta K$ decouple?

## 10. Pair with

- **Taylor, *Classical Mechanics*, ch. 11.** The standard undergraduate treatment: two masses and three springs, weak coupling and beats, the double pendulum, several coupled pendulums, and normal coordinates, with more pictures than these notes. Taylor writes the eigenproblem as $(K - \omega^2M)\mathbf{a} = 0$ and solves by determinants; read him for the physics of each example, these notes for why the theorems hold.
- **Goldstein, Poole & Safko, *Classical Mechanics*, ch. 6.** The congruence transformation $A^\top MA = I$, $A^\top KA = \Omega^2$ done properly, and the forced and damped cases. Terser and more general than Taylor; read after §3 here.
- **French, *Vibrations and Waves*, ch. 5.** The best intuition for coupled oscillators and beats, with experiments you can do. Read it alongside §2 here.
- **Strang, MIT 18.06, lecture 27 "Positive definite matrices and minima"**, and for the $A^\top CA$ framework his *Introduction to Applied Mathematics*, ch. 1–2 (also the core of his MIT 18.085 lectures). Where $K = B^\top CB$ comes from, why it is positive definite, and the mechanical meaning of every factor. Strang's springs-and-masses example is Example 4.2 with equal masses.
- **MIT 8.03 (Lee)** on OCW: the lectures on coupled oscillators, normal modes and the beaded string, with live demonstrations of beats and the dispersion relation.
