# CM-6: Small Oscillations and Normal Modes (the Generalized Eigenproblem)

**Time:** ~12 hours · **Prerequisites:** CM-3, LA-8 · **Unlocks:** CM-7, F2 (continuum limit), QM-3 (oscillator), SM-5 (Debye)

> **Conventions in this module.** Generalized coordinates $q = (q_1, \dots, q_n)^\top \in \mathbb{R}^n$ are measured from a stable equilibrium, so the equilibrium is $q = 0$. Lagrangian $L = T - V$ (course convention). The **mass matrix** $M$ is symmetric positive definite; the **stiffness matrix** $K$ is symmetric. The mode ansatz is $q(t) = \mathbf{a}\cos(\omega t - \delta)$ with real mode shape $\mathbf{a}$, angular frequency $\omega \ge 0$ and phase $\delta$. Modes are numbered $\omega_1 \le \omega_2 \le \cdots \le \omega_n$ and normalized by $\mathbf{a}_i^\top M\mathbf{a}_j = \delta_{ij}$, using the **$M$-inner product** $\langle\mathbf{u}, \mathbf{v}\rangle_M = \mathbf{u}^\top M\mathbf{v}$. Normal coordinates are $\xi_j$, with $q = A\xi$, $A = [\mathbf{a}_1 \cdots \mathbf{a}_n]$ and $\Omega^2 = \operatorname{diag}(\omega_1^2, \dots, \omega_n^2)$. Spring constants are $k$; because $k$ is taken, the wavenumber in the string example is $\kappa$. $\dagger$ is conjugate transpose, $\succ 0$ means positive definite. SI units.

## 0. Where this sits

From CM-3 you must be able to write $T$ and $V$ for a constrained system in generalized coordinates and get the equations of motion from Euler–Lagrange; the double pendulum of CM-3 reappears here as Example 4.4. From LA-8 you need the spectral theorem for symmetric matrices, what positive definite means and the second-derivative test, and the three axioms of an inner product. From LA-10 the Cholesky factorization $M = LL^\top$ is used once, at the decisive moment (Theorem 3.4); if you have not done LA-10, the one fact you need is stated there.

What this module makes possible: the inertia tensor of CM-7 is the same kind of object as $M$ here, a symmetric positive definite matrix that defines a kinetic energy, and principal axes are normal modes without the $K$. F2 gets its standing waves and its orthogonal sine basis from the loaded string of Example 4.5. QM-3 quantizes the single oscillator that every normal mode becomes; SM-5 counts those oscillators to get the heat capacity of a solid. For an ME student this module is **modal analysis**: it is what finite-element vibration software computes when you click "solve."

## 1. The problem that forced this idea

Anything with mass and stiffness vibrates: a bridge deck, a turbine blade, an engine on its mounts, the atoms in a molecule or a crystal. A system with $n$ degrees of freedom near equilibrium has $n$ coupled second-order ODEs, and when you write them out they look hopeless. Push on one mass and every other mass responds; the coordinates refuse to separate.

The first people to see through this were working on strings and chains. In the 1730s Daniel Bernoulli studied the oscillations of a hanging chain and of a string loaded with beads, and found that these systems have special motions in which every point moves at the same frequency, each at its own fixed amplitude, the whole shape rising and falling in unison. For a string carrying $n$ beads he found $n$ such motions with $n$ distinct frequencies (Example 4.5 derives them). He then claimed something stronger: that *every* motion of the system is a sum of these special motions. Euler doubted the generality of the claim, d'Alembert had solved the continuous string by a different route (1747), and the three argued for decades. Bernoulli was right, and the argument is the origin of Fourier series: the question "is every vibration a sum of modes?" was on the table a century before Fourier.

Lagrange, in the *Mécanique analytique* (1788), turned the physical insight into algebra. Expand $T$ and $V$ to second order about equilibrium, assume a solution in which every coordinate oscillates at the same frequency, and the differential equations collapse to a polynomial equation of degree $n$ in $\omega^2$ whose roots are the mode frequencies. He showed the roots are real and, for a stable equilibrium, positive, and that the amplitude patterns belonging to different roots satisfy a relation we would now call orthogonality. He did all this in coordinates, with no matrices, because matrices did not exist yet; Cayley's matrix algebra came in 1858. Lagrange had diagonalized a symmetric matrix seventy years early, while working out how a string sounds.

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

You met this in LA-7, where the two-mass, three-spring chain gave $\ddot{\mathbf x} = -(K/m)\mathbf x$ and you solved it by diagonalizing $K$. The new thing in this module is what happens when $M$ is *not* a multiple of the identity: unequal masses, or the double pendulum, where $M$ is not even diagonal. Then $K\mathbf{a} = \omega^2 M\mathbf{a}$ is a **generalized eigenproblem**. The obvious move is to rewrite it as $M^{-1}K\mathbf{a} = \omega^2\mathbf{a}$ and apply the spectral theorem to $M^{-1}K$. Try it on $M = \operatorname{diag}(1, 2)$, $K = \begin{pmatrix}2 & -1\\-1 & 2\end{pmatrix}$: $M^{-1}K = \begin{pmatrix}2 & -1\\-1/2 & 1\end{pmatrix}$, not symmetric. The spectral theorem does not apply to it as a matrix, and its eigenvectors turn out *not* to be orthogonal in the ordinary sense (Example 4.2: their dot product is $1/2$).

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
This is the $A^\top CA$ structure of Strang's lecture on positive definite matrices: $K$ is automatically symmetric and positive semidefinite ($q^\top Kq = \sum_s c_s e_s^2 \ge 0$), and it is positive definite exactly when $Bq = 0$ forces $q = 0$, that is, when no nonzero displacement leaves every spring unstretched. A displacement in the null space $N(B)$ (LA-3) is a **rigid motion** or a **mechanism**, and §3.5 shows it is a zero-frequency mode. Examples 4.3 and 4.6 and the project build $K$ this way.

**Stability.** If $K$ is positive definite, $V$ has a strict local minimum at $q = 0$ (LA-8's second-derivative test), and the equilibrium is stable: the exact energy $E = T + V$ is conserved and $\ge V_0$, and the region $\{E \le V_0 + \epsilon\}$ of phase space shrinks onto the point $(q, \dot q) = (0, 0)$ as $\epsilon \to 0$, so a motion that starts near equilibrium with little energy stays near it. (This is the Lagrange–Dirichlet theorem; the argument is exact, not linearized.) If $K$ has a negative eigenvalue, Theorem 3.5 produces a solution of (3.3) that grows exponentially, and the equilibrium is unstable. If $K$ is only positive semidefinite, the linear theory cannot decide; the cubic terms do.

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
The operator $M^{-1}K$ is **self-adjoint with respect to $\langle\cdot,\cdot\rangle_M$**, even though its matrix in the standard basis is not symmetric. LA-8's spectral theorem is a theorem about self-adjoint operators on inner product spaces, not about symmetric arrays of numbers; applied here it says that $M^{-1}K$ has real eigenvalues and an eigenbasis orthonormal in $\langle\cdot,\cdot\rangle_M$. Theorems 3.1 and 3.2 are what that abstract statement says in coordinates. The next theorem makes the connection to a plain symmetric matrix concrete, which is also how the computer does it.

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

■ **What just happened.** The change of variable $\mathbf{y} = L^\top\mathbf{a}$ is a change of coordinates in which the kinetic energy becomes the ordinary Euclidean one: $\dot q^\top M\dot q = \dot q^\top LL^\top\dot q = |L^\top\dot q|^2$. In those coordinates the mass matrix is the identity, the problem is a standard symmetric eigenproblem, and orthogonality is ordinary orthogonality. Translating back, ordinary orthogonality of the $\mathbf{y}$'s is $M$-orthogonality of the $\mathbf{a}$'s. This is exactly what `scipy.linalg.eigh(K, M)` does (LAPACK `sygv`: Cholesky, transform, symmetric eigensolver, transform back), and it is why LA-10's Cholesky factorization matters for anything that vibrates. Any factorization $M = RR^\top$ works, including the symmetric square root $M^{1/2}$ from LA-8; Cholesky is the cheapest.

Completeness with degeneracies is the part Theorem 3.2 could not give you. Theorem 3.2 only speaks about distinct frequencies; Theorem 3.4 says that even when $\omega_i^2 = \omega_j^2$ an $M$-orthonormal pair of shapes exists in that eigenspace (Gram–Schmidt in the $M$-inner product produces one). Every motion is a sum of modes because the modes span $\mathbb{R}^n$. Bernoulli's claim is proved, for finitely many degrees of freedom.

**Theorem 3.5 (Stability).** With $M$, $K$ as above and $\omega_j^2$ the generalized eigenvalues: all $\omega_j^2 \ge 0$ if and only if $K$ is positive semidefinite, and all $\omega_j^2 > 0$ if and only if $K$ is positive definite.

*Proof.* ($\Leftarrow$) Each $\omega_j^2 = \mathbf{a}_j^\top K\mathbf{a}_j/\mathbf{a}_j^\top M\mathbf{a}_j$ by Theorem 3.1 step 4; the denominator is positive, and the numerator is $\ge 0$ ($> 0$) when $K$ is PSD (PD). ($\Rightarrow$) Take any $\mathbf{v} \ne 0$ and expand it in the basis of Theorem 3.4: $\mathbf{v} = \sum_j c_j\mathbf{a}_j$ with $c_j = \langle\mathbf{a}_j, \mathbf{v}\rangle_M$, not all zero. Using $K\mathbf{a}_j = \omega_j^2M\mathbf{a}_j$ and (3.10),
$$\mathbf{v}^\top K\mathbf{v} = \sum_{i,j}c_ic_j\,\mathbf{a}_i^\top K\mathbf{a}_j = \sum_{i,j}c_ic_j\,\omega_j^2\,\mathbf{a}_i^\top M\mathbf{a}_j = \sum_j\omega_j^2c_j^2,$$
which is $\ge 0$ if all $\omega_j^2 \ge 0$, and $> 0$ if all $\omega_j^2 > 0$. ∎

When some $\omega_j^2 < 0$, write $\omega_j^2 = -\gamma_j^2$; the ansatz $q = \mathbf{a}_je^{st}$ gives $s^2 = -\gamma_j^2\cdot(-1)$, i.e. $q = \mathbf{a}_j(c_1e^{\gamma_jt} + c_2e^{-\gamma_jt})$, and almost every initial condition grows. That is the linear-theory statement that a saddle of $V$ is unstable.

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
In general $C$ wrecks the decoupling: $A^\top CA$ need not be diagonal, and the modes of the damped system are complex and must be found from a $2n$-dimensional first-order problem (LA-7's companion form). Engineers almost always avoid this with an assumption.

**Theorem 3.6 (Rayleigh damping decouples).** If $C = \alpha M + \beta K$ for scalars $\alpha, \beta \ge 0$, then in normal coordinates (3.14) becomes $n$ independent damped driven oscillators
$$\ddot\xi_j + 2\zeta_j\omega_j\dot\xi_j + \omega_j^2\xi_j = f_j(t), \qquad 2\zeta_j\omega_j = \alpha + \beta\omega_j^2, \qquad f = A^\top F. \tag{3.15}$$

*Proof.* Substitute $q = A\xi$ into (3.14) and multiply on the left by $A^\top$: $A^\top MA\ddot\xi + A^\top CA\dot\xi + A^\top KA\xi = A^\top F$. By (3.11) the first and third matrices are $I$ and $\Omega^2$. The middle one is $A^\top(\alpha M + \beta K)A = \alpha I + \beta\Omega^2$, diagonal with entries $\alpha + \beta\omega_j^2$. ∎

The modal damping ratio is $\zeta_j = \alpha/(2\omega_j) + \beta\omega_j/2$: mass-proportional damping ($\alpha$) damps low modes most, stiffness-proportional damping ($\beta$) damps high modes most. Two parameters cannot match measured damping in more than two modes, but the decoupling is worth more than the accuracy, and in practice one fits $\alpha, \beta$ to two measured modes. (Exercise C3 generalizes to the Caughey series $C = M\sum_la_l(M^{-1}K)^l$, which also decouples and has as many parameters as you like.)

**Frequency response.** Drive with $F(t) = \operatorname{Re}[\mathbf{F}e^{i\omega t}]$ and look for the steady state $\xi_j = \operatorname{Re}[\hat\xi_je^{i\omega t}]$. From (3.15), $(\omega_j^2 - \omega^2 + 2i\zeta_j\omega_j\omega)\hat\xi_j = \mathbf{a}_j^\top\mathbf{F}$, so
$$\hat q(\omega) = \sum_{j=1}^n\frac{\mathbf{a}_j\,(\mathbf{a}_j^\top\mathbf{F})}{\omega_j^2 - \omega^2 + 2i\zeta_j\omega_j\omega}. \tag{3.16}$$
Each term is a single-oscillator resonance with its peak near $\omega = \omega_j$ and height $\approx |\mathbf{a}_j^\top\mathbf{F}|/(2\zeta_j\omega_j^2)$. A sweep of the drive frequency shows one peak per mode (the project does this for a 3-DOF chain); this is how modal frequencies and damping ratios are measured on a real structure with a shaker and accelerometers. Two readable features of (3.16): a mode is not excited by a force applied where its shape has a node ($\mathbf{a}_j^\top\mathbf{F} = 0$), and it is not seen by a sensor placed at a node of its shape (the $j$-th term has no component there).

**Modal superposition as the engineering method.** A finite-element model replaces a continuous structure by $n \sim 10^5$–$10^7$ nodal displacements, assembles $K = B^\top CB$ and $M$ from element matrices, and solves (3.5) for the lowest $r \ll n$ modes with a Lanczos or subspace iteration (LA-13; the Cholesky reduction of Theorem 3.4 is inside). Every later analysis is done with $q \approx \sum_{j\le r}\mathbf{a}_j\xi_j$: $r$ scalar equations (3.15) instead of $n$ coupled ones, with Rayleigh or modal damping assumed so that they decouple. Truncation to $r$ modes is accurate when the loads have no frequency content near the discarded modes, which is why the first question in a vibration analysis is always "what is the frequency range of the excitation?" The exaggerated shapes you see in a modal-analysis animation are the $\mathbf{a}_j$.

### 3.8 The quantum preview

Each normal coordinate obeys $\ddot\xi_j = -\omega_j^2\xi_j$, with Hamiltonian $\frac12(p_j^2 + \omega_j^2\xi_j^2)$: the harmonic oscillator of QM-3. Quantizing one mode gives energy levels $\hbar\omega_j(n_j + \frac12)$; quantizing all of them gives $E = \sum_j\hbar\omega_j(n_j + \frac12)$, and a quantum of mode $j$ in a crystal is a **phonon**. The Debye model (SM-5) takes the dispersion relation $\omega(\kappa)$ of the loaded string in its long-wavelength form $\omega = c\kappa$, counts the modes of a block of solid, and gets the $T^3$ heat capacity; the finite number of atoms sets the maximum frequency, exactly as $N$ masses give $N$ modes in Example 4.5. The continuum limit of Example 4.5 is a field, a displacement $y(x, t)$ at every point, and it is a system of infinitely many oscillators, one per wavenumber. Quantizing those is quantum field theory, and the photon is the phonon of the electromagnetic field.
