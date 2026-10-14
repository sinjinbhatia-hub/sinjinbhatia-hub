# CM-2: Least Action: From Fermat to Hamilton's Principle; Euler–Lagrange

**Time:** ~8 hours · **Prerequisites:** F6 (Euler–Lagrange equation, fundamental lemma, second variation), CM-1 (Newton's laws; the constraint problems of §3.10) · **Unlocks:** CM-3, and through it CM-4 to CM-8; the same ideas return in EM-6 and QM-9

> **Conventions in this module.** SI units; $g = 9.81\ \mathrm{m/s^2}$; the vertical coordinate ($y$ or $z$) points **up**, so $V = mgy$. Generalized coordinates are $q = (q_1, \dots, q_n)$; a dot is $d/dt$. The **Lagrangian** $L(q, \dot q, t) = T - V$ is a function of $2n + 1$ independent slots: $\partial L/\partial\dot q_i$ treats $q$, $\dot q$, $t$ as independent, while $d/dt$ is the total derivative along a path (F6's conventions with $x \to t$, $y \to q$). The **action** is $S[q] = \int_{t_1}^{t_2}L\,dt$. A **variation** is $q + \varepsilon\eta$ with $\eta(t_1) = \eta(t_2) = 0$; **stationary** means $\frac{d}{d\varepsilon}S[q + \varepsilon\eta]\big|_{\varepsilon = 0} = 0$ for every such $\eta$, which is weaker than "minimum". The **Euler–Lagrange expression** is $\mathcal E_i := \frac{d}{dt}\frac{\partial L}{\partial\dot q_i} - \frac{\partial L}{\partial q_i}$, signed so that in Cartesian coordinates $\mathcal E = m\ddot x - F$; the equations of motion are $\mathcal E_i = 0$. **Generalized momentum** $p_i := \partial L/\partial\dot q_i$; **generalized force** $\mathcal F_i := \partial L/\partial q_i$ (not CM-3's $Q_i$). New coordinates are $Q_j$, with Jacobian $J_{ij} = \partial q_i/\partial Q_j$. Charge is $q_e$, so it is never confused with a coordinate; $\Phi$ and $\mathbf A$ are the scalar and vector potentials; bold $\mathbf E$ is the electric field, plain $E$ is energy. No summation convention.

## 0. Where this sits

F6 built the machine: a functional $\int L(x, y, y')\,dx$ is stationary exactly when the Euler–Lagrange equation holds. CM-1 made you solve three constrained problems in Newtonian form and watch unknown forces eat the page. This module joins the two. It proves that Newton's second law, for any system whose forces come from a potential, is *equivalent* to one statement about whole paths: the true motion makes $S = \int(T - V)\,dt$ stationary. Then it proves the property that makes this more than a curiosity: the Euler–Lagrange equations take the same form in every coordinate system.

It unlocks CM-3 (the coordinate freedom becomes a method for constrained systems), CM-4 (conservation laws from $p_i$ and the total-derivative theorem), CM-8 (Hamilton–Jacobi from §3.8), EM-6 (§3.9's charged-particle Lagrangian), and QM-9 ($e^{iS/\hbar}$).

## 1. The problem that forced this idea

In the first century, Hero of Alexandria showed in his *Catoptrica* that a light ray reflected by a flat mirror from one point to another takes the shortest available path, and that the equal-angle law follows. A law of motion phrased as an optimization over paths: for sixteen centuries it stayed a curiosity.

Fermat made it a principle. In 1662 he argued that light crossing into a new medium takes the path of least *time*, and derived the law of refraction, provided light is slower in the denser medium. Descartes had derived the same sine law in 1637 from a mechanical model that needed light to be *faster* there. The Cartesians objected to more than the speed. Claude Clerselier, defending Descartes, wrote to Fermat that least time was a moral principle, not a physical one: how could light know in advance which path would be quickest? The speed question was settled in 1850, when Foucault measured light travelling more slowly in water than in air. The objection about foresight never quite went away; §2 and §8 answer it.

Maupertuis carried the idea to matter. In a 1744 memoir to the Paris Academy he proposed that nature minimizes a quantity he called **action** (for light, speed times distance along the path), and in 1746 applied it to collisions as evidence of a Creator who economizes. Also in 1744, Euler published a precise version in an appendix to his *Methodus inveniendi*: for a body under central forces, $\int v\,ds$ along the orbit is stationary. Euler's version was the sound one, and he credited Maupertuis. When Samuel König claimed in 1751 that Leibniz had stated the principle first, citing a copy of a letter whose original was never produced, Euler led the Berlin Academy's verdict of forgery (1752), and Voltaire answered with a satire that Frederick the Great had burned. The quarrel left the principle a reputation for mysticism it did not deserve.

Lagrange stripped the metaphysics out. His papers of 1760–61 applied his $\delta$-method (F6 §1) to dynamics, and his *Mécanique analytique* (1788) reduced mechanics to analysis; its preface warns that the book contains no figures, only algebraic operations. Hamilton, in *On a General Method in Dynamics* (1834) and its *Second Essay* (1835), gave the principle its modern form: with the start and end **times** fixed, the motion makes $\int(T - V)\,dt$ stationary. He came from optics, where he had packed ray optics into one function of a ray's endpoints, and he did the same for mechanics (§3.8).

The twentieth century is why you learn this. When Newton's laws were replaced, the action principle survived: a relativistic free particle's action is $-mc^2$ times its elapsed proper time; Hilbert obtained Einstein's field equations from an action in 1915; Noether proved in 1918 that symmetries of an action give conservation laws (CM-4); and Dirac (1933) and Feynman (1940s) found that quantum mechanics weights each path by $e^{iS/\hbar}$, with classical mechanics surviving the interference (§3.10).

The working problem is narrower, and you met it in CM-1 §3.10. $\mathbf F = m\ddot{\mathbf r}$ is a vector equation tied to Cartesian components: in polar coordinates the acceleration is not $(\ddot r, \ddot\theta)$, and getting it right means differentiating rotating unit vectors. For a pendulum you want $\theta$; Newton hands you $x$, $y$, an unknown tension, and seven steps of algebra to make the tension cancel. What was needed is a formulation that is **the same in every coordinate system**, so you can choose coordinates that fit the problem and never see the forces enforcing the constraints. A number attached to a whole path does not care which coordinates describe the path. That is the idea.

## 2. Intuition

### Competing paths

Fix two events: the system is at configuration $A$ at time $t_1$ and at $B$ at time $t_2$. Every curve in the $(t, q)$ plane from $(t_1, A)$ to $(t_2, B)$ is a candidate history. All of them leave at the same time and arrive at the same time.

```
   q
   ^                          B
   |        .-------.        *      each curve is a whole history q(t);
   |      /    ___    \     /|      S assigns each one a number
   |     /  .-'   '-.  \__ / |
   |    / /    true    \  /  |      the true one: S does not change
   |   //      path     \/   |      to first order when you wiggle it
   |  *--------------------- |
   |  A                      |
   +--+----------------------+----> t
      t1                     t2
```

Hamilton's principle says: compute $S = \int(T - V)\,dt$ for each, and the true history is one where $S$ is stationary. Note what is fixed. The question is never "how long does it take to get from $A$ to $B$?" It is "given that the system left $A$ at $t_1$ and arrived at $B$ at $t_2$, what did it do in between?"

### Chop time into pieces (F6's picture, with $t$)

Take one coordinate $q$ with $T = \frac12 m\dot q^2$ and potential $V(q)$. Cut $[t_1, t_2]$ into $N + 1$ steps of length $h$, hold $q_0 = A$ and $q_{N+1} = B$ fixed, and let $q_1, \dots, q_N$ be free. Approximate the velocity on each step by $(q_{k+1} - q_k)/h$ and the potential by its average at the two ends of the step:
$$
S_h(q_1, \dots, q_N) = \sum_{k=0}^{N}h\left[\frac m2\left(\frac{q_{k+1} - q_k}{h}\right)^2 - \frac{V(q_k) + V(q_{k+1})}{2}\right]. \tag{2.1}
$$
This is an ordinary function of $N$ variables, and its stationary points are where its gradient vanishes.

▶ **Try it first.** $q_k$ appears in two kinetic terms (steps $k - 1$ and $k$) and in two potential terms. Compute $\partial S_h/\partial q_k$ and set it to zero. Which equation do you recognize?

The kinetic terms give $m(q_k - q_{k-1})/h - m(q_{k+1} - q_k)/h$; each potential term gives $-\frac h2V'(q_k)$. So
$$
\frac{\partial S_h}{\partial q_k} = -h\left[m\,\frac{q_{k+1} - 2q_k + q_{k-1}}{h^2} + V'(q_k)\right]. \tag{2.2}
$$
The fraction is the central-difference approximation to $\ddot q(t_k)$, and $-V'(q_k)$ is the force. The gradient of the discretized action is **$-h$ times "mass times acceleration minus force"** at each node; setting it to zero gives $m\ddot q = -V'(q)$ as $h \to 0$. You have guessed Theorem 3.2. (The discrete equation is also the Störmer–Verlet integrator of CM-8 §3.11.) The project minimizes (2.1) with L-BFGS-B and lands on the ODE solution.

Component $k$ of (2.2) involves only $q_{k-1}, q_k, q_{k+1}$: a condition on the whole path breaks into one condition per instant. That is why "the path makes $S$ stationary" is not teleology. It is a list of local equations written as one sentence.

### Why $T - V$: fast versus high

Throw a ball straight up so it leaves the ground ($y = 0$) at $t = 0$ and lands at $t = \tau$. Compare the candidate paths $y_h(t) = 4h\,t(\tau - t)/\tau^2$, parabolas of peak height $h$ ($h = 0$ is the ball sitting on the ground). With $T = \frac12 m\dot y^2$ and $V = mgy$, direct integration gives
$$
S(h) = \int_0^\tau\left(\tfrac12 m\dot y_h^2 - mgy_h\right)dt = \frac{8mh^2}{3\tau} - \frac23 mgh\tau. \tag{2.3}
$$
Read the two terms. The kinetic term grows like $h^2$: to go higher and still land on time you must move faster, and speed is penalized quadratically. The potential term *rewards* height linearly, because $T - V$ is smaller where $V$ is larger. The action balances "spend time up high" against "don't rush". Setting $dS/dh = 0$ gives $h^* = g\tau^2/8$. Check against the real throw: landing at $\tau$ needs launch speed $g\tau/2$, so the real peak is $(g\tau/2)^2/(2g) = g\tau^2/8$. The family contains the true path, and stationary action picked it out.

Flip the sign to $T + V$ and the potential term also penalizes height; the "best" path dips to $h = -g\tau^2/8$, the flight of a ball whose gravity points up. Section 3.3 turns this into a theorem.

### What to expect

Two guesses to carry into §3. First, stationary action is *equivalent* to Newton's second law, not a deeper law that implies it. Second, $S$ is a number attached to a path, and a path is the same path whatever coordinates describe it. So the stationarity condition, the Euler–Lagrange equation, should look the same in every coordinate system. The second guess is the important one, and its proof (Theorem 3.5) is the center of the module.

## 3. Formal development

### 3.1 The action and its first variation

A **path** is a $C^2$ function $q: [t_1, t_2] \to \mathbb R^n$, and the Lagrangian $L(q, \dot q, t)$ is a $C^2$ function of its $2n + 1$ slots. The **action** is
$$
S[q] = \int_{t_1}^{t_2}L\big(q(t), \dot q(t), t\big)\,dt. \tag{3.1}
$$
F6 (3.3)–(3.5) computed the derivative of such a functional in a direction $\eta$: differentiate under the integral with the chain rule, then integrate the $\dot\eta$ term by parts. With several coordinates (F6 §3.7), and keeping the boundary term instead of discarding it, the result for any $C^1$ functions $\eta_i$ is
$$
\frac{d}{d\varepsilon}S[q + \varepsilon\eta]\Big|_{\varepsilon = 0} = -\int_{t_1}^{t_2}\sum_{i=1}^n\mathcal E_i[q](t)\,\eta_i(t)\,dt + \Big[\sum_{i=1}^n\frac{\partial L}{\partial\dot q_i}\,\eta_i\Big]_{t_1}^{t_2}, \tag{3.2}
$$
with the Euler–Lagrange expressions
$$
\mathcal E_i[q] := \frac{d}{dt}\frac{\partial L}{\partial\dot q_i} - \frac{\partial L}{\partial q_i}. \tag{3.3}
$$
We keep the boundary term because §3.8 needs it.

**Lemma 3.1 (Stationary ⟺ Euler–Lagrange).** A $C^2$ path $q$ makes $S$ stationary if and only if $\mathcal E_i[q](t) = 0$ for every $i$ and every $t \in [t_1, t_2]$.

*Proof.* (⇐) If every $\mathcal E_i \equiv 0$, the integral in (3.2) vanishes, and the boundary term vanishes because $\eta$ does at both ends. (⇒) Vary one coordinate at a time: $\eta_j = 0$ for $j \ne i$. Then (3.2) says $\int\mathcal E_i\,\eta_i\,dt = 0$ for every admissible $\eta_i$. $\mathcal E_i$ is continuous, so by the fundamental lemma (F6 Theorem 3.4) $\mathcal E_i \equiv 0$. ∎

F6 Theorem 3.5 proved (⇒), which is all you need to *find* candidate paths. The equivalence with Newton below needs (⇐) as well.

### 3.2 Hamilton's principle, and Newton from it

**Hamilton's principle.** For a system whose forces derive from a potential $V(q, t)$, the motion from configuration $A$ at $t_1$ to $B$ at $t_2$ is a path that makes $S = \int_{t_1}^{t_2}(T - V)\,dt$ stationary among all $C^2$ paths with those endpoints.

In this module the principle is not an axiom. It is a theorem about Newtonian mechanics.

**Theorem 3.2 (Hamilton ⟺ Newton).** Let $N$ particles with constant masses $m_a$ and positions $\mathbf r_a = (x_a, y_a, z_a)$ move in a potential $V(\mathbf r_1, \dots, \mathbf r_N, t)$, and let
$$
L = \sum_{a=1}^N\tfrac12 m_a\left(\dot x_a^2 + \dot y_a^2 + \dot z_a^2\right) - V(\mathbf r_1, \dots, \mathbf r_N, t).
$$
A $C^2$ path makes $S$ stationary if and only if, for every $a$,
$$
m_a\ddot{\mathbf r}_a = -\nabla_aV, \tag{3.4}
$$
where $\nabla_a = (\partial/\partial x_a, \partial/\partial y_a, \partial/\partial z_a)$: Newton's second law with force $\mathbf F_a = -\nabla_aV$.

▶ **Try it first.** There are $3N$ coordinates. Write the E–L equation for $x_1$ alone. Which terms of $L$ contain $\dot x_1$? Which contain $x_1$?

**Derivation.**

1. The coordinates are the $3N$ numbers $x_a, y_a, z_a$. By Lemma 3.1, stationarity is equivalent to $\mathcal E = 0$ for each of them.
2. **$x$-components.** Only $\frac12 m_a\dot x_a^2$ contains $\dot x_a$, and $V$ contains no velocities, so $\partial L/\partial\dot x_a = m_a\dot x_a$. $T$ contains no positions, so $\partial L/\partial x_a = -\partial V/\partial x_a$. Hence $\mathcal E_{x_a} = \frac{d}{dt}(m_a\dot x_a) + \partial V/\partial x_a = m_a\ddot x_a + \partial V/\partial x_a$ (mass constant), and $\mathcal E_{x_a} = 0 \iff m_a\ddot x_a = -\partial V/\partial x_a$.
3. **$y$-components.** $\partial L/\partial\dot y_a = m_a\dot y_a$ and $\partial L/\partial y_a = -\partial V/\partial y_a$, so $\mathcal E_{y_a} = m_a\ddot y_a + \partial V/\partial y_a$, and $\mathcal E_{y_a} = 0 \iff m_a\ddot y_a = -\partial V/\partial y_a$.
4. **$z$-components.** $\partial L/\partial\dot z_a = m_a\dot z_a$ and $\partial L/\partial z_a = -\partial V/\partial z_a$, so $\mathcal E_{z_a} = m_a\ddot z_a + \partial V/\partial z_a$, and $\mathcal E_{z_a} = 0 \iff m_a\ddot z_a = -\partial V/\partial z_a$.
5. Steps 2–4 for every $a$ are the components of (3.4). ∎

■ **What just happened.** In Cartesian coordinates the E–L expression is literally $m_a\ddot{\mathbf r}_a - \mathbf F_a$, component by component. "Stationary action" and "$\mathbf F = m\mathbf a$" are the same equations for conservative forces, time-dependent ones included. Nothing is gained yet; the gain is Theorem 3.5. (If the particles interact through $V = \sum_{a<b}U(|\mathbf r_a - \mathbf r_b|)$, then $\nabla_aU = -\nabla_bU$: Newton's third law, CM-1 §3.1, comes free.)

### 3.3 Why $T - V$ and not something else

**Theorem 3.3 (The relative sign is forced).** For one particle on a line with $T = \frac12 m\dot x^2$, let $L = \alpha T + \beta V(x)$ with constants $\alpha \ne 0$ and $\beta$. The Euler–Lagrange equation is
$$
m\ddot x = \frac\beta\alpha\,V'(x), \tag{3.5}
$$
so it reproduces Newton's $m\ddot x = -V'(x)$ for every potential $V$ if and only if $\beta/\alpha = -1$.

▶ **Try it first.** Compute $\partial L/\partial\dot x$ and $\partial L/\partial x$, form $\mathcal E$, and ask what $\beta/\alpha$ must be for the potential $V(x) = x$.

**Derivation.**

1. $\partial L/\partial\dot x = \alpha m\dot x$ and $\partial L/\partial x = \beta V'(x)$.
2. $\mathcal E = \alpha m\ddot x - \beta V'(x)$. Setting $\mathcal E = 0$ and dividing by $\alpha$ gives (3.5).
3. If $\beta/\alpha = -1$, (3.5) is Newton's law.
4. Conversely, suppose (3.5) agrees with Newton for every $V$. For $V(x) = x$ (a constant force $-1$), (3.5) says $m\ddot x = \beta/\alpha$ and Newton says $m\ddot x = -1$. So $\beta/\alpha = -1$.
5. If $\alpha = 0$, then $\mathcal E = -\beta V'(x)$ contains no acceleration: it is an algebraic condition on $x$, not an equation of motion. ∎

With $\beta = -\alpha$, $L = \alpha(T - V)$, and by Theorem 3.4(b) the factor $\alpha$ changes nothing, so take $\alpha = 1$. With many particles in three dimensions the same computation runs component by component (steps 2–4 of Theorem 3.2 with the constants carried along).

■ **What just happened.** The minus sign is what Newton's law requires, nothing more mysterious. Why the kinetic term is $\frac12 m\dot x^2$ and not another function of $\dot x$ is a separate question: Landau & Lifshitz §3–4 answer it with Galilean relativity, and Exercise B1 shows why something like that is needed. Note also that $T - V$ is not conserved; $T + V$ is, and it is the Beltrami quantity (F6 (3.11)) of $L = T - V$ up to sign.

### 3.4 The Lagrangian is not unique

**Theorem 3.4 (Total derivatives and rescaling).** Let $L(q, \dot q, t)$ be $C^2$.

(a) Let $F(q, t)$ be $C^3$ and
$$
L' = L + \frac{dF}{dt}, \qquad \frac{dF}{dt} = \sum_{j=1}^n\frac{\partial F}{\partial q_j}\dot q_j + \frac{\partial F}{\partial t}. \tag{3.6}
$$
Then $\mathcal E_i[L'] = \mathcal E_i[L]$ for every path, so $L$ and $L'$ have the same equations of motion.

(b) For a constant $c \ne 0$, $\mathcal E_i[cL] = c\,\mathcal E_i[L]$, so $cL$ and $L$ have the same equations of motion.

▶ **Try it first.** For (a), compute $S'[q] - S[q]$ for a path with fixed endpoints. What does the difference depend on?

**Proof via the action.** $S'[q] = S[q] + \int_{t_1}^{t_2}\frac{dF}{dt}\,dt = S[q] + F(q(t_2), t_2) - F(q(t_1), t_1)$. All competing paths share their endpoints, so the added term is the same number for every one of them. Adding a constant to a function does not move its stationary points, so $S$ and $S'$ are stationary on the same paths, and Lemma 3.1 turns this into equal sets of solutions. ∎

That shows the two Lagrangians have the same *solutions*. A direct computation shows more: the extra term contributes nothing to $\mathcal E_i$ on *any* path.

**Proof by computation.** Let $G = dF/dt$ as in (3.6) and $\partial_j = \partial/\partial q_j$.

1. $G$ contains $\dot q_i$ only in the term $\partial_iF\,\dot q_i$, and $\partial_iF$ is a function of $(q, t)$, so $\partial G/\partial\dot q_i = \partial_iF$.
2. Along a path, $\frac{d}{dt}\partial_iF\big(q(t), t\big) = \sum_j\partial_j\partial_iF\,\dot q_j + \partial_t\partial_iF$.
3. Differentiating (3.6) with respect to the slot $q_i$: $\partial G/\partial q_i = \sum_j\partial_i\partial_jF\,\dot q_j + \partial_i\partial_tF$.
4. Mixed partials of the $C^2$ function $F$ commute, so steps 2 and 3 agree and $\mathcal E_i[G] = 0$ identically. $\mathcal E_i$ is linear in $L$, so $\mathcal E_i[L'] = \mathcal E_i[L] + \mathcal E_i[G] = \mathcal E_i[L]$.
5. (b) $\mathcal E_i$ is built from derivatives of $L$, so $\mathcal E_i[cL] = c\,\mathcal E_i[L]$, which vanishes exactly when $\mathcal E_i[L]$ does, since $c \ne 0$. ∎

■ **What just happened.** Different-looking Lagrangians can describe the same motion, so compare Lagrangians by their E–L expressions, not their formulas. You will meet $dF/dt$ again in a Galilean boost (CM-4 §3.3), in a gauge transformation (§3.9), and in CM-8's generating functions. $F$ may depend on $q$ and $t$ but not on $\dot q$; otherwise $dF/dt$ contains $\ddot q$. Exercise A4 tests the boundary.

### 3.5 Covariance: the Euler–Lagrange equations in any coordinates

A **point transformation** replaces coordinates $q$ by new ones $Q$ via
$$
q_i = q_i(Q_1, \dots, Q_n, t), \qquad i = 1, \dots, n,
$$
a smooth map ($C^3$ is enough), invertible for each $t$, with nonsingular Jacobian $J_{ij} = \partial q_i/\partial Q_j$. Time dependence is allowed, which covers moving frames such as a rotating one. A path $Q(t)$ gives a path $q(t) = q(Q(t), t)$ with velocity
$$
\dot q_i = \sum_{k=1}^n\frac{\partial q_i}{\partial Q_k}\dot Q_k + \frac{\partial q_i}{\partial t}. \tag{3.7}
$$
The **transformed Lagrangian** is the same number written in the new variables, $\tilde L(Q, \dot Q, t) := L\big(q(Q, t), \dot q(Q, \dot Q, t), t\big)$ with $\dot q$ from (3.7). Write $\tilde{\mathcal E}_j$ for its E–L expressions.

**Theorem 3.5 (Covariance of the Euler–Lagrange equations).** For every $C^2$ path $Q(t)$, with $q(t) = q(Q(t), t)$,
$$
\tilde{\mathcal E}_j[Q](t) = \sum_{i=1}^nJ_{ij}\,\mathcal E_i[q](t), \qquad j = 1, \dots, n. \tag{3.8}
$$
Since $J$ is invertible, $\tilde{\mathcal E}_j \equiv 0$ for all $j$ if and only if $\mathcal E_i \equiv 0$ for all $i$: a path solves the E–L equations of $L$ in the old coordinates exactly when it solves those of $\tilde L$ in the new ones.

Two facts about (3.7) do all the work. Both concern the coordinate change only, not $L$.

**Lemma A (dots cancel).** $\dfrac{\partial\dot q_i}{\partial\dot Q_j} = \dfrac{\partial q_i}{\partial Q_j}$.

*Proof.* In (3.7), $\dot q_i$ is linear in the slots $\dot Q_k$, with coefficients $\partial q_i/\partial Q_k$ that depend only on $(Q, t)$. Differentiating with respect to $\dot Q_j$ picks out the coefficient of $\dot Q_j$. ∎

**Lemma B ($d/dt$ commutes with $\partial/\partial Q$).** Along any path, $\dfrac{d}{dt}\dfrac{\partial q_i}{\partial Q_j} = \dfrac{\partial\dot q_i}{\partial Q_j}$.

*Proof.* The left side is the total time derivative of a function of $(Q, t)$: $\sum_k\frac{\partial^2q_i}{\partial Q_k\partial Q_j}\dot Q_k + \frac{\partial^2q_i}{\partial t\,\partial Q_j}$. The right side differentiates (3.7) with respect to the slot $Q_j$ with $\dot Q$ held fixed: $\sum_k\frac{\partial^2q_i}{\partial Q_j\partial Q_k}\dot Q_k + \frac{\partial^2q_i}{\partial Q_j\,\partial t}$. They agree term by term because mixed partials commute. ∎

▶ **Try it first.** Compute $\partial\tilde L/\partial\dot Q_j$ by the chain rule (through $\dot q$ only) and simplify with Lemma A. Take $d/dt$. Separately compute $\partial\tilde L/\partial Q_j$ (through both $q$ and $\dot q$). Subtract, and watch which terms cancel.

**Proof of Theorem 3.5.** Write $p_i = \partial L/\partial\dot q_i$, evaluated at the point $(q, \dot q, t)$ that corresponds to $(Q, \dot Q, t)$.

1. $\tilde L$ depends on $\dot Q_j$ only through the $\dot q_i$, because the $q_i$ depend on $(Q, t)$ alone. Chain rule, then Lemma A:
$$\frac{\partial\tilde L}{\partial\dot Q_j} = \sum_ip_i\,\frac{\partial\dot q_i}{\partial\dot Q_j} = \sum_ip_i\,J_{ij}.$$
2. Differentiate along the path with the product rule, and use Lemma B on the second term:
$$\frac{d}{dt}\frac{\partial\tilde L}{\partial\dot Q_j} = \sum_i\dot p_i\,J_{ij} + \sum_ip_i\,\frac{\partial\dot q_i}{\partial Q_j}.$$
3. $\tilde L$ depends on $Q_j$ through *both* $q_i$ and $\dot q_i$. Chain rule:
$$\frac{\partial\tilde L}{\partial Q_j} = \sum_i\frac{\partial L}{\partial q_i}\,J_{ij} + \sum_ip_i\,\frac{\partial\dot q_i}{\partial Q_j}.$$
4. Subtract step 3 from step 2. The terms $\sum_ip_i\,\partial\dot q_i/\partial Q_j$ cancel exactly, leaving
$$\tilde{\mathcal E}_j = \sum_i\Big(\dot p_i - \frac{\partial L}{\partial q_i}\Big)J_{ij} = \sum_iJ_{ij}\,\mathcal E_i,$$
which is (3.8).
5. In matrix form (3.8) reads $\tilde{\boldsymbol{\mathcal E}} = J^\top\boldsymbol{\mathcal E}$. $J^\top$ is invertible because $J$ is, so $\tilde{\boldsymbol{\mathcal E}} = \mathbf 0 \iff \boldsymbol{\mathcal E} = \mathbf 0$ at each instant. ∎

**Why it had to work.** By construction $\tilde S[Q] = S[q]$ for every path, and fixed endpoints in one set of coordinates are fixed endpoints in the other. So the two actions are stationary on the same paths, which is the guess of §2. Making this a proof needs (3.2) for the nonlinear family $q(Q + \varepsilon\eta, t)$ plus the fundamental lemma; the chain-rule proof above avoids that and gives the stronger identity (3.8), valid off solutions too.

■ **What just happened.** The E–L expressions transform with the transposed Jacobian, the rule obeyed by the components of a gradient. In Cartesian coordinates $\boldsymbol{\mathcal E}$ is $m\mathbf a - \mathbf F$; in other coordinates $\tilde{\mathcal E}_j$ is that same vector dotted into the coordinate direction $\partial\mathbf r/\partial Q_j$. So you may write $L = T - V$ in *whatever* coordinates fit the problem, apply (3.3), and the equations come out right, with no unit-vector derivatives and no fictitious-force bookkeeping. The proof looks like index shuffling; it is the license for all of Lagrangian mechanics. CM-3 §3.4 extends it to coordinates on a constraint surface, where the map is no longer invertible on the whole space.

### 3.6 Generalized momentum and generalized force

**Definition.** For a Lagrangian $L(q, \dot q, t)$, the **generalized momentum** conjugate to $q_i$ and the **generalized force** are
$$
p_i := \frac{\partial L}{\partial\dot q_i}, \qquad \mathcal F_i := \frac{\partial L}{\partial q_i}, \qquad\text{so that the E–L equation reads}\quad \dot p_i = \mathcal F_i. \tag{3.9}
$$
In Cartesian coordinates with velocity-independent $V$, $p_i = m\dot x_i$ and $\mathcal F_i = -\partial V/\partial x_i$: momentum and force. In other coordinates neither is what its name suggests. In polar coordinates (Example 4.2), $p_\theta = mr^2\dot\theta$ is an angular momentum, with units kg·m²/s, and $\mathcal F_r = mr\dot\theta^2 - V'(r)$ contains a term from the $r$-dependence of $T$ that no agent applies. For a charged particle (§3.9), $\mathbf p = m\mathbf v + q_e\mathbf A$. Step 1 of Theorem 3.5's proof showed how momenta transform: $\tilde p_j = \sum_iJ_{ij}p_i$, the same rule as (3.8).

**Cyclic coordinates.** If $L$ does not contain $q_k$, then $\mathcal F_k = 0$ and (3.9) says $p_k$ is constant. This is the first conservation law the Lagrangian formalism gives for free; CM-4 builds Noether's theorem from it.

### 3.7 Is it a minimum?

The name says "least action". Here is what is true, for $L = \frac12 m\dot q^2 - V(q, t)$ in one coordinate (several coordinates change nothing essential).

▶ **Try it first.** Expand $S[q + \varepsilon\eta]$ to second order in $\varepsilon$, using a Taylor expansion of $V$ about $q(t)$. What multiplies $\varepsilon^2$? When is the expansion exact?

**Derivation.**

1. Kinetic part, exactly: $\frac12 m(\dot q + \varepsilon\dot\eta)^2 = \frac12 m\dot q^2 + \varepsilon m\dot q\dot\eta + \frac12\varepsilon^2m\dot\eta^2$.
2. Potential part: $V(q + \varepsilon\eta) = V(q) + \varepsilon V'(q)\eta + \frac12\varepsilon^2V''(q)\eta^2 + O(\varepsilon^3)$, primes meaning $\partial/\partial q$.
3. Integrate. The $\varepsilon^1$ term is the first variation, which vanishes on a solution. So on a solution $q$,
$$
S[q + \varepsilon\eta] = S[q] + \frac{\varepsilon^2}{2}\,\delta^2S[\eta] + O(\varepsilon^3), \qquad \delta^2S[\eta] = \int_{t_1}^{t_2}\Big(m\dot\eta^2 - V''\big(q(t), t\big)\eta^2\Big)\,dt. \tag{3.10}
$$
4. If $V$ is at most quadratic in $q$, the series in step 2 stops at $\varepsilon^2$, and (3.10) holds exactly, with no $O(\varepsilon^3)$, for perturbations of any size. ∎

This is F6 (3.17) with $L_{qq} = -V''$, $L_{q\dot q} = 0$, $L_{\dot q\dot q} = m$.

**No force, or a constant force ($V'' = 0$).** Then $S[q + \eta] = S[q] + \frac12 m\int\dot\eta^2\,dt$ exactly. The integral is positive unless $\dot\eta \equiv 0$, which with $\eta(t_1) = 0$ means $\eta \equiv 0$. The true path is the **unique global minimum**. Examples 4.1 and 4.4 put numbers on this.

**The harmonic oscillator.** $V = \frac12 kq^2$, $\omega = \sqrt{k/m}$, period $T_{\text{osc}} = 2\pi/\omega$, interval length $\tau = t_2 - t_1$. Now $V'' = k > 0$, and the two terms of (3.10) compete.

**Theorem 3.6 (A minimum only for short times).** For the oscillator, the true path between fixed endpoints is a strict minimum of $S$ if $\tau < T_{\text{osc}}/2 = \pi/\omega$. If $\tau > \pi/\omega$, then $\eta(t) = \sin\big(\pi(t - t_1)/\tau\big)$ gives $S[q + \varepsilon\eta] < S[q]$ for every $\varepsilon \ne 0$, while $\eta = \sin\big(n\pi(t - t_1)/\tau\big)$ with $n\pi/\tau > \omega$ gives $S[q + \varepsilon\eta] > S[q]$: the path is a saddle, not a maximum.

▶ **Try it first.** Put $\eta = \sin(\pi s/\tau)$, $s = t - t_1$, into (3.10) and compute both integrals. When is the result negative?

**Derivation.**

1. By step 4 above, $S[q + \eta] - S[q] = \frac12\mathcal Q[\eta]$ exactly, with $\mathcal Q[\eta] = \int_0^\tau(m\dot\eta^2 - k\eta^2)\,ds$. (Any admissible perturbation is $\varepsilon\eta$ with $\varepsilon = 1$.)
2. Expand $\eta$ in its Fourier sine series on $[0, \tau]$: $\eta(s) = \sum_{n\ge1}b_n\sin(n\pi s/\tau)$. Every term vanishes at both ends, as $\eta$ does.
3. The cosine coefficients of $\dot\eta$ come from integration by parts, with boundary terms killed by $\eta(0) = \eta(\tau) = 0$: $\frac2\tau\int_0^\tau\dot\eta\cos(n\pi s/\tau)\,ds = \frac{n\pi}\tau b_n$. The mean of $\dot\eta$ is $(\eta(\tau) - \eta(0))/\tau = 0$.
4. Parseval's identity (F2 Theorem 3.8, applied to the odd extension of $\eta$ and the even extension of $\dot\eta$) gives $\int_0^\tau\eta^2\,ds = \frac\tau2\sum_nb_n^2$ and $\int_0^\tau\dot\eta^2\,ds = \frac\tau2\sum_n(n\pi/\tau)^2b_n^2$. Therefore
$$
\mathcal Q[\eta] = \frac{m\tau}2\sum_{n=1}^\infty b_n^2\left[\left(\frac{n\pi}\tau\right)^2 - \omega^2\right]. \tag{3.11}
$$
5. If $\tau < \pi/\omega$, then $(n\pi/\tau)^2 \ge (\pi/\tau)^2 > \omega^2$ for every $n$: every bracket is positive, and $\mathcal Q[\eta] > 0$ unless all $b_n = 0$, i.e. $\eta \equiv 0$. Strict minimum.
6. If $\tau > \pi/\omega$, the $n = 1$ bracket is negative. For $\eta = \sin(\pi s/\tau)$ (only $b_1 = 1$), $\mathcal Q = \frac{m\tau}2\big((\pi/\tau)^2 - \omega^2\big) < 0$, so $S$ decreases along $\varepsilon\eta$ for every $\varepsilon \ne 0$. For $\eta = \sin(n\pi s/\tau)$ with $n\pi/\tau > \omega$ the bracket is positive and $S$ increases. ∎

■ **What just happened.** Kinetic energy penalizes a wiggle in proportion to its frequency squared; the potential rewards displacement at the fixed rate $k$. On a short interval even the slowest admissible wiggle, half a sine wave, is too fast to win; past half a period it wins. At exactly $\tau = \pi/\omega$ the $n = 1$ bracket is zero: every solution starting at $q_1$ reaches $-q_1$ half a period later whatever its amplitude, so a whole family of true paths shares those endpoints. That endpoint is a **conjugate point** (F6 §3.11). Each further multiple of $\pi/\omega$ adds one more descending direction (Exercise B4).

**In general.** If $\partial^2L/\partial\dot q^2 > 0$, as it is for any genuine kinetic energy, F6's Legendre argument (F6 Theorem 3.8) shows the true path is **never a maximum**: fast wiggles always raise $S$. Before the first conjugate point it is a local minimum; past it, a saddle (Jacobi's condition: F6 §3.11, proved in Gelfand & Fomin ch. 5). So the accurate name is the principle of **stationary** action, and the equations of motion come from the first variation alone.

### 3.8 The action as a function of its endpoints

Suppose that for endpoints near $(q_1, t_1)$ and $(q_2, t_2)$ there is exactly one true path $q_{\text{cl}}$, depending smoothly on the endpoints (true away from conjugate points). Its action is then a function of the endpoints, $S_{\text{cl}}(q_2, t_2; q_1, t_1) := S[q_{\text{cl}}]$. Hamilton called it the **principal function**.

▶ **Try it first.** Move $q_2$ slightly. The new true path differs from the old by some $\xi(t)$ with $\xi(t_1) = 0$ but $\xi(t_2) \ne 0$. Use (3.2) with its boundary term.

**Derivation** (one coordinate; with several, repeat for each).

1. Replace $q_2$ by $q_2 + \varepsilon$. The true path becomes $q_{\text{cl}} + \varepsilon\xi + O(\varepsilon^2)$ with $\xi(t_1) = 0$, $\xi(t_2) = 1$.
2. By (3.2), $\frac{d}{d\varepsilon}S_{\text{cl}}\big|_{\varepsilon = 0} = -\int\mathcal E[q_{\text{cl}}]\,\xi\,dt + \big[p\,\xi\big]_{t_1}^{t_2}$. The integral vanishes because $q_{\text{cl}}$ solves E–L; the boundary term is $p(t_2)\cdot1 - p(t_1)\cdot0$.
3. So, with the mirror-image argument at the other end ($\xi(t_1) = 1$, $\xi(t_2) = 0$, which flips the sign),
$$
\frac{\partial S_{\text{cl}}}{\partial q_2} = p_2, \qquad \frac{\partial S_{\text{cl}}}{\partial q_1} = -p_1. \tag{3.12}
$$
4. For the time dependence, let the endpoint slide along the true path: $S_{\text{cl}}\big(q_{\text{cl}}(t), t; q_1, t_1\big) = \int_{t_1}^tL\,dt'$, since a piece of a true path is the true path between its own ends. Differentiate in $t$. The right side gives $L$; the chain rule on the left gives $\partial S_{\text{cl}}/\partial t_2 + (\partial S_{\text{cl}}/\partial q_2)\,\dot q$. With (3.12),
$$
\frac{\partial S_{\text{cl}}}{\partial t_2} = L - p_2\dot q_2 = -E, \tag{3.13}
$$
where $E = p\dot q - L$ is the energy function; for $L = \frac12 m\dot q^2 - V$ it is $T + V$. ∎

**Free particle.** On the straight line, $S_{\text{cl}} = \frac{m(q_2 - q_1)^2}{2(t_2 - t_1)}$, so $\partial S_{\text{cl}}/\partial q_2 = m(q_2 - q_1)/(t_2 - t_1) = mv = p_2$ and $\partial S_{\text{cl}}/\partial t_2 = -\frac12 mv^2 = -E$. With $m = 1$, $q_1 = t_1 = 0$, $q_2 = t_2 = 1$: $S_{\text{cl}} = 0.5$, $\partial S/\partial q_2 = 1$, $\partial S/\partial t_2 = -0.5$ (project, Part 3).

■ **What just happened.** The slope of the action with respect to where you end is the momentum you arrive with; the slope with respect to when you end is minus the energy. Substituting (3.12) into (3.13) gives a PDE for $S$ alone, $\partial S/\partial t + E(q, \partial S/\partial q, t) = 0$: the Hamilton–Jacobi equation of CM-8 §3.10. In QM-8 it is the leading order of the Schrödinger equation for $\psi \approx e^{iS/\hbar}$.

### 3.9 Velocity-dependent potentials and the Lorentz force

Some velocity-dependent forces still fit the scheme. Suppose $U(q, \dot q, t)$ is a function with $L = T - U$. In Cartesian coordinates, steps 2–4 of Theorem 3.2 with $V$ replaced by $U$ give $m\ddot x_i = F_i$ with
$$
F_i = -\frac{\partial U}{\partial x_i} + \frac{d}{dt}\frac{\partial U}{\partial\dot x_i}. \tag{3.14}
$$
The case that matters is a charge $q_e$ in electric and magnetic fields described by potentials $\Phi(\mathbf r, t)$ and $\mathbf A(\mathbf r, t)$, with $\mathbf E = -\nabla\Phi - \partial_t\mathbf A$ and $\mathbf B = \nabla\times\mathbf A$ (EM-6 derives these; here they are definitions).

**Theorem 3.7 (Lorentz force from a Lagrangian).** With
$$
L = \tfrac12 m|\dot{\mathbf r}|^2 - q_e\Phi(\mathbf r, t) + q_e\,\dot{\mathbf r}\cdot\mathbf A(\mathbf r, t), \tag{3.15}
$$
the Euler–Lagrange equations are equivalent to
$$
m\ddot{\mathbf r} = q_e\left(\mathbf E + \dot{\mathbf r}\times\mathbf B\right). \tag{3.16}
$$

▶ **Try it first.** Compute $p_x = \partial L/\partial\dot x$. Its time derivative contains $dA_x/dt$ along the path; expand that with the chain rule. Collect the velocity terms and compare with $(\dot{\mathbf r}\times\mathbf B)_x$ written out using $\mathbf B = \nabla\times\mathbf A$.

**Derivation.** Write $(x_1, x_2, x_3) = (x, y, z)$, $v_j = \dot x_j$, $\partial_j = \partial/\partial x_j$.

1. In components, $L = \frac12 m\sum_jv_j^2 - q_e\Phi + q_e\sum_jv_jA_j$.
2. $p_i = \partial L/\partial v_i = mv_i + q_eA_i(\mathbf r, t)$.
3. Along the path $A_i$ changes because time passes and because the particle moves: $\frac{d}{dt}A_i = \partial_tA_i + \sum_jv_j\partial_jA_i$. So $\dot p_i = m\ddot x_i + q_e\partial_tA_i + q_e\sum_jv_j\partial_jA_i$.
4. $\partial L/\partial x_i = -q_e\partial_i\Phi + q_e\sum_jv_j\,\partial_iA_j$.
5. $\mathcal E_i = \dot p_i - \partial L/\partial x_i = 0$ rearranges to
$$m\ddot x_i = q_e\left(-\partial_i\Phi - \partial_tA_i\right) + q_e\sum_jv_j\left(\partial_iA_j - \partial_jA_i\right).$$
6. The first bracket is the field component $(\mathbf E)_i$.
7. The second term is $q_e(\mathbf v\times\mathbf B)_i$. For $i = x$: $(\mathbf v\times\mathbf B)_x = v_yB_z - v_zB_y$ with $B_z = \partial_xA_y - \partial_yA_x$ and $B_y = \partial_zA_x - \partial_xA_z$, so $(\mathbf v\times\mathbf B)_x = v_y\partial_xA_y + v_z\partial_xA_z - v_y\partial_yA_x - v_z\partial_zA_x$. Add and subtract $v_x\partial_xA_x$: this is $\sum_jv_j\partial_xA_j - \sum_jv_j\partial_jA_x$, the $i = x$ case of step 5. The relabeling $x \to y \to z \to x$ maps the formulas for $\mathbf v\times\mathbf B$ and $\nabla\times\mathbf A$ to themselves, so it carries this computation to the $y$ and $z$ components.
8. Steps 5–7 give (3.16), and every step is reversible. ∎

■ **What just happened.** The magnetic force depends on velocity and does no work, so no ordinary $V(\mathbf r)$ can produce it; the term $q_e\dot{\mathbf r}\cdot\mathbf A$ can. In the language of (3.14), $U = q_e\Phi - q_e\dot{\mathbf r}\cdot\mathbf A$. The price is that the canonical momentum is $\mathbf p = m\mathbf v + q_e\mathbf A$, not $m\mathbf v$. QM-2 quantizes $\mathbf p$, which is why $\mathbf A$ appears in the Schrödinger equation.

**Gauge freedom is Theorem 3.4.** The potentials are not unique: for any $\chi(\mathbf r, t)$, $\mathbf A' = \mathbf A + \nabla\chi$ and $\Phi' = \Phi - \partial_t\chi$ give the same $\mathbf E$ and $\mathbf B$ (Exercise B6). The Lagrangian changes by
$$
L' - L = q_e\left(\partial_t\chi + \dot{\mathbf r}\cdot\nabla\chi\right) = q_e\frac{d\chi}{dt}, \tag{3.17}
$$
a total derivative, so by Theorem 3.4 the equations of motion cannot change. Example 4.6 uses this.

### 3.10 Feynman's picture (a preview, stated not derived)

In Feynman's formulation of quantum mechanics (1948, building on Dirac 1933), the amplitude to go from $A$ at $t_1$ to $B$ at $t_2$ is a sum over *all* paths of $e^{iS[\text{path}]/\hbar}$, each with the same magnitude; $\hbar = 1.055\times10^{-34}$ J·s. Near a stationary path $S$ changes only at second order (3.10), so neighboring paths have nearly equal phases and add up. Elsewhere $S$ changes at first order, neighboring phases rotate rapidly, and the contributions cancel. For the ball of Example 4.1, $|S_{\text{cl}}|/\hbar \approx 3.8\times10^{34}$ radians: only the stationary path survives. That is the classical limit, and it is why the principle is about stationarity, not minimality; saddles survive the interference too, with the conjugate points of §3.7 contributing phase factors. QM-9 does this properly.

### 3.11 What the principle does not cover

Hamilton's principle as stated needs forces from a potential $V(q, t)$ or a generalized potential $U(q, \dot q, t)$. **Friction and drag have no such potential** in general: they remove energy, while the E–L expression of $T - V$ contains no damping term and the velocity-dependent forces produced by a generalized potential linear in the velocities, like the magnetic one, do no work: for $U = U_0(\mathbf r, t) - \sum_jA_j(\mathbf r, t)v_j$, (3.14) gives them as $\sum_jv_j(\partial_iA_j - \partial_jA_i)$, as in step 5 of Theorem 3.7, and their power $\sum_{i,j}v_iv_j(\partial_iA_j - \partial_jA_i)$ vanishes because the bracket is antisymmetric in $i, j$. CM-3 §3.9 handles dissipation by adding Rayleigh's dissipation function to the E–L equations by hand. Exercise C3 shows a trick that works for one linearly damped oscillator and why it does not generalize. **Rolling without slipping** and other non-holonomic constraints also fall outside this form; CM-3 §3.2 says what to do.

## 4. Worked examples

### Example 4.1: Newton from the principle, for a thrown ball

A ball of mass $m$ moves in three dimensions with $L = \frac12 m(\dot x^2 + \dot y^2 + \dot z^2) - mgz$. Theorem 3.2's components: $\partial L/\partial\dot x = m\dot x$ and $\partial L/\partial x = 0$ give $m\ddot x = 0$; the same for $y$; $\partial L/\partial\dot z = m\dot z$ and $\partial L/\partial z = -mg$ give $m\ddot z = -mg$. A parabola.

Now the vertical throw of §2 with $m = 1$ kg, $\tau = 1$ s, $z(0) = z(\tau) = 0$. The true path is $z_{\text{cl}} = \frac12 gt(\tau - t)$, with peak $g\tau^2/8 = 1.22625$ m, the $h^*$ of (2.3). Its action:
$$
\int_0^\tau\tfrac12 m\dot z_{\text{cl}}^2\,dt = \frac{mg^2}2\int_0^\tau\Big(\frac\tau2 - t\Big)^2dt = \frac{mg^2\tau^3}{24}, \qquad \int_0^\tau mgz_{\text{cl}}\,dt = \frac{mg^2}2\int_0^\tau t(\tau - t)\,dt = \frac{mg^2\tau^3}{12},
$$
so $S_{\text{cl}} = mg^2\tau^3\big(\tfrac1{24} - \tfrac1{12}\big) = -mg^2\tau^3/24 = -4.00984$ J·s. The "sit on the ground" path has $S = 0$, which is larger. Because $V$ is linear, (3.10) is exact with $V'' = 0$: for $\eta = \sin(\pi t)$,
$$S[z_{\text{cl}} + \eta] = S_{\text{cl}} + \tfrac12 m\int_0^1\pi^2\cos^2(\pi t)\,dt = -4.00984 + \frac{\pi^2}4 = -4.00984 + 2.46740.$$
The true path is the global minimum. The project's L-BFGS-B minimizer of (2.1) reproduces $z_{\text{cl}}$ at every node to better than $10^{-6}$: for a constant force the discrete equation from (2.2) is exact on parabolas, so only the optimizer's tolerance remains.

### Example 4.2: Polar coordinates, and the centripetal term for free

A particle in a plane with a central potential $V(r)$. With $x = r\cos\theta$, $y = r\sin\theta$, (3.7) gives $\dot x = \dot r\cos\theta - r\dot\theta\sin\theta$ and $\dot y = \dot r\sin\theta + r\dot\theta\cos\theta$, so $\dot x^2 + \dot y^2 = \dot r^2 + r^2\dot\theta^2$ and
$$
L = \tfrac12 m\left(\dot r^2 + r^2\dot\theta^2\right) - V(r). \tag{4.1}
$$
By Theorem 3.5 we may apply (3.3) to (4.1) directly.

- $r$: $p_r = m\dot r$, $\mathcal F_r = mr\dot\theta^2 - V'(r)$. E–L: $m(\ddot r - r\dot\theta^2) = -V'(r)$.
- $\theta$: $p_\theta = mr^2\dot\theta$, $\mathcal F_\theta = 0$. E–L: $\frac{d}{dt}(mr^2\dot\theta) = 0$.

Newton in polar form needs the rotating unit vectors differentiated twice to get $\mathbf a = (\ddot r - r\dot\theta^2)\hat{\mathbf r} + (r\ddot\theta + 2\dot r\dot\theta)\hat{\boldsymbol\theta}$. The $r$ equation is the $\hat{\mathbf r}$ component, its centripetal term $-r\dot\theta^2$ arriving from $\partial T/\partial r$ with no unit vectors in sight. The $\theta$ equation, $mr(r\ddot\theta + 2\dot r\dot\theta) = 0$, is $r$ times the $\hat{\boldsymbol\theta}$ component. Theorem 3.5 predicts these factors: with $J = \begin{pmatrix}\cos\theta & -r\sin\theta\\ \sin\theta & r\cos\theta\end{pmatrix}$, (3.8) gives $\tilde{\mathcal E}_r = \cos\theta\,\mathcal E_x + \sin\theta\,\mathcal E_y$ and $\tilde{\mathcal E}_\theta = -r\sin\theta\,\mathcal E_x + r\cos\theta\,\mathcal E_y$. The project checks this in SymPy for an arbitrary $V(x, y)$; Exercise B2 does it by hand.

$\theta$ is cyclic, so $p_\theta = mr^2\dot\theta$ (angular momentum, not "mass times $\dot\theta$") is conserved: Kepler's second law in one line. For a circular orbit in $V = -k/r$, $\ddot r = 0$ gives $mr\dot\theta^2 = k/r^2$, so $\dot\theta^2 = k/(mr^3)$, Kepler's third law; with $m = k = r = 1$, $\dot\theta = 1$ and the period is $2\pi$ (Exercise C2).

### Example 4.3: The pendulum in one line, the inclined plane in another

CM-1 §3.10(a) took seven steps, an unknown tension, and a mid-course switch to the angle to reach $\ddot\theta = -(g/\ell)\sin\theta$ (CM-1 (3.13)). Here: the bob sits at $(\ell\sin\theta, -\ell\cos\theta)$, its speed is $\ell\dot\theta$, and its height is $-\ell\cos\theta$, so
$$
L = \tfrac12 m\ell^2\dot\theta^2 + mg\ell\cos\theta \quad\Longrightarrow\quad \frac{d}{dt}\left(m\ell^2\dot\theta\right) = -mg\ell\sin\theta.
$$
That is the entire derivation. The tension never appears because $\theta$ describes only configurations the rod allows (CM-3 proves that dropping the constraint force is legitimate, and §3.8 there recovers the tension as a multiplier). The project integrates CM-1's Cartesian equations, with tension $mv^2/\ell - mgy/\ell$ from CM-1 (3.12), against this $\theta$ equation from $\theta_0 = 1$ rad, $\ell = 1$ m: the angle histories agree to about $10^{-11}$ rad over 5 s.

**Inclined plane.** A block slides on a frictionless incline at angle $\alpha$; let $s$ be the distance travelled down the slope. $T = \frac12 m\dot s^2$; the height drops by $s\sin\alpha$, so $V = -mgs\sin\alpha$. $L = \frac12 m\dot s^2 + mgs\sin\alpha$ gives $m\ddot s = mg\sin\alpha$, with no normal force anywhere. At $\alpha = 30°$, $\ddot s = 4.905$ m/s².

### Example 4.4: Free-particle action along two paths

Take $m = 1$ with $x(0) = 0$ and $x(1) = 1$. Compare the straight line $x = t$ with a bent path that goes out to $x = 2$ at $t = \frac12$ and back to $x = 1$ at $t = 1$, uniformly on each half.

*Straight:* $\dot x = 1$, so $S = \frac12\cdot1\cdot1^2\cdot1 = 0.5$.

*Bent:* first half $\dot x = 2/0.5 = 4$, second half $\dot x = (1 - 2)/0.5 = -2$. $S = \frac12(4^2)(0.5) + \frac12(2^2)(0.5) = 4 + 1 = 5$.

*Check with §3.7.* The perturbation $\eta = x_{\text{bent}} - x_{\text{straight}}$ rises from 0 to 1.5 at $t = \frac12$ and returns to 0, so $\dot\eta = \pm3$ and $\frac12 m\int_0^1\dot\eta^2\,dt = \frac12\cdot9 = 4.5$. So $S_{\text{bent}} = 0.5 + 4.5 = 5$ exactly. (The bent path has a kink, so it is not $C^2$, but the free-particle identity needs only a piecewise $C^1$ $\eta$: the cross term $\int m\dot x_{\text{cl}}\dot\eta\,dt = m\cdot1\cdot[\eta]_0^1 = 0$ regardless.) The detour costs ten times the action, because it forces speed and kinetic energy grows like speed squared.

### Example 4.5: The oscillator's second variation and the $T/2$ threshold

Take $m = k = 1$: $\omega = 1$, $T_{\text{osc}} = 2\pi$, threshold $\tau = \pi \approx 3.1416$. Endpoints $q(0) = 0$, $q(\tau) = 1$. The true path is $q_{\text{cl}} = \sin t/\sin\tau$ (it solves $\ddot q = -q$ and fits both ends), with
$$
S_{\text{cl}} = \int_0^\tau\frac{\cos^2t - \sin^2t}{2\sin^2\tau}\,dt = \frac{\sin2\tau}{4\sin^2\tau} = \tfrac12\cot\tau.
$$
Perturb by $\varepsilon\sin(\pi t/\tau)$. By (3.11) with $b_1 = \varepsilon$, and the exactness in §3.7 step 4,
$$
S[q_{\text{cl}} + \varepsilon\eta] = S_{\text{cl}} + c(\tau)\,\varepsilon^2, \qquad c(\tau) = \frac\tau4\left[\Big(\frac\pi\tau\Big)^2 - 1\right]. \tag{4.2}
$$

- $\tau = 2 < \pi$: $S_{\text{cl}} = \frac12\cot2 = -0.22883$ and $c = \frac12(\pi^2/4 - 1) = \frac12(2.46740 - 1) = 0.73370$. A minimum.
- $\tau = 4 > \pi$: $S_{\text{cl}} = \frac12\cot4 = 0.43185$ and $c = \pi^2/16 - 1 = 0.61685 - 1 = -0.38315$. Along $\sin(\pi t/4)$ the action goes *down*, quadratically and without bound, since the expansion is exact.

The project confirms both coefficients by quadrature for $\varepsilon$ up to 2. The Hessian of the discrete action (2.1), here $(m/h)\,\mathrm{tridiag}(-1, 2, -1) - hk\,I$, has no negative eigenvalue at $\tau = 2$ and exactly one at $\tau = 4$, with eigenvector matching $\sin(\pi t/4)$ to overlap above $0.9999$. On the $\tau = 4$ problem L-BFGS-B runs off toward $S_h \to -\infty$: there is no minimum to find. The discretized action is a convex quadratic below a threshold and an unbounded nonconvex one above it. The threshold is where the lowest Hessian eigenvalue, $(4m/h)\sin^2\!\big(\pi/(2(N+1))\big) - hk$, crosses zero: $\tau_h = \frac{2(N+1)}{\omega}\sin\frac{\pi}{2(N+1)} = \frac\pi\omega\big(1 - \frac{\pi^2}{24(N+1)^2} + \dots\big)$. That is $\pi$ up to $O(h^2)$, and $\pi - 8.0\times10^{-6}$ for the project's $N = 400$.

### Example 4.6: Charged particle in a uniform magnetic field

Uniform $\mathbf B = B\hat{\mathbf z}$, no electric field. Take the **symmetric gauge** $\mathbf A = \frac B2(-y, x, 0)$: $(\nabla\times\mathbf A)_z = \partial_x(Bx/2) - \partial_y(-By/2) = B$, and the other components vanish since $A_z = 0$ and nothing depends on $z$. With $\Phi = 0$, (3.15) is
$$
L = \tfrac12 m(\dot x^2 + \dot y^2 + \dot z^2) + \tfrac12 q_eB\,(x\dot y - y\dot x).
$$

- $x$: $p_x = m\dot x - \frac12 q_eBy$, so $\dot p_x = m\ddot x - \frac12 q_eB\dot y$; $\partial L/\partial x = \frac12 q_eB\dot y$. E–L: $m\ddot x = q_eB\dot y$.
- $y$: $p_y = m\dot y + \frac12 q_eBx$, $\dot p_y = m\ddot y + \frac12 q_eB\dot x$; $\partial L/\partial y = -\frac12 q_eB\dot x$. E–L: $m\ddot y = -q_eB\dot x$.
- $z$: $m\ddot z = 0$.

Check with (3.16): $\dot{\mathbf r}\times B\hat{\mathbf z} = B(\dot y, -\dot x, 0)$. ✓ Note that $p_x \ne m\dot x$.

Solve with $u = \dot x + i\dot y$: $\dot u = \frac{q_eB}m(\dot y - i\dot x) = -i\omega_cu$, where $\omega_c = q_eB/m$ is the **cyclotron frequency**. So $u = u_0e^{-i\omega_ct}$: the velocity turns clockwise (seen from $+z$, for $q_eB > 0$) at constant speed, and the orbit is a circle of radius $R = |u_0|/\omega_c = mv/(q_eB)$.

Numbers: $m = 2$, $q_e = 1$, $B = 1$, start at the origin with $\mathbf v = (1, 0, 0)$. Then $\omega_c = 0.5$, $T_c = 2\pi/\omega_c = 4\pi \approx 12.566$, $R = 2$, and integrating $u$ gives $x + iy = \frac{iu_0}{\omega_c}(e^{-i\omega_ct} - 1) = 2i(e^{-it/2} - 1)$: a circle centered at $(0, -2)$, through $(2, -2)$ after a quarter period. The project asserts each number.

**Gauge check.** The Landau gauge $\mathbf A_L = (0, Bx, 0)$ of CM-3 Example 4.8 gives the same $\mathbf B$, and $\mathbf A_L - \mathbf A = \frac B2(y, x, 0) = \nabla\chi$ with $\chi = \frac12 Bxy$. By (3.17) and Theorem 3.4 the E–L equations agree; the project confirms it term for term, although $p_x$ and $p_y$ differ between the gauges. Canonical momentum depends on the gauge; the motion does not.

## 5. Why it matters and where it connects

- **CM-3** is Theorem 3.5 put to work: coordinates on the constraint surface, $T - V$ in them, then (3.3). CM-3 §3.4 gives the d'Alembert derivation; §3.8 recovers constraint forces as multipliers.
- **CM-4** uses §3.6 (a cyclic coordinate conserves its $p_k$) and Theorem 3.4 (a symmetry of $L$ up to $dF/dt$ is a symmetry of the motion; the Galilean boost needs it).
- **CM-5** starts from Example 4.2's $L$ and conserved $p_\theta$.
- **CM-8** builds $H = \sum_ip_i\dot q_i - L$ from this module's $p_i$, turns (3.12)–(3.13) into the Hamilton–Jacobi equation (§3.10 there), and its symplectic integrators include the discrete E–L equation of §2.
- **EM-6**: the Lagrangian (3.15) and gauge freedom as a total derivative (3.17).
- **QM-2, QM-9**: the momentum that becomes $-i\hbar\nabla$ is $p_i = \partial L/\partial\dot q_i$ (so $m\mathbf v + q_e\mathbf A$ for a charge), and the path integral weights paths by $e^{iS/\hbar}$.
- **ML-2**: minimizing (2.1) is smooth unconstrained optimization, and for the oscillator Theorem 3.6 says when it is convex (up to Example 4.5's $O(h^2)$ shift of the threshold); Example 4.5's negative eigenvalue is the saddle that convexity assumptions exclude.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are drills (actions along explicit paths, CM-1's Atwood machine in one line, the wrong sign, harmless added terms, spherical momenta). B1–B7 are the module; B2 (covariance by hand), B4 (counting descending directions), B5 (the oscillator's principal function), and B7 (extend `project.py`) are essential. C1–C3: the brachistochrone as mechanics, the Kepler circle against an ellipse, a Lagrangian for damping.

## 7. Coding project

`project.py`, in seven parts: (1) the discrete action (2.1) with gradient (2.2), minimized by L-BFGS-B for a projectile, an oscillator ($\tau = 2$) and a pendulum ($\tau = 0.8$ s), asserted against shooting solutions of the ODE, plus (2.3) and Example 4.1; (2) $S[q_{\text{cl}} + \varepsilon\eta]$ plotted and its slope at $0$ asserted below $10^{-6}$ on the true path and above $10^{-2}$ on a straight line, plus Example 4.4; (3) Example 4.5, the discrete Hessian's negative eigenvalues for $\tau = 2, 4, 5, 7$ and its convexity threshold $\tau_h$, L-BFGS-B running away past $T/2$, and the endpoint derivatives (3.12)–(3.13); (4) a SymPy `euler_lagrange` verifying Theorems 3.2–3.5 and 3.7 with arbitrary $V$, $\Phi$, $\mathbf A$, including Theorem 3.5 for the time-dependent rotating frame of Exercise B3; (5) cyclotron motion, integrating the SymPy E–L equations of Example 4.6; (6) CM-1's Cartesian pendulum against the one-line $\theta$ equation; (7) the numbers quoted in the solutions.

**Build it yourself:** write `discrete_action` (value and gradient) and `euler_lagrange` before reading them. The gradient is the test of whether you understood (2.2): it is three lines.

## 8. Where you'll struggle

- **"Nature minimizes action" sounds as if the particle looks ahead.** It does not. Theorem 3.2 makes the principle *equivalent* to $\mathbf F = m\mathbf a$ for conservative forces, and (2.2) shows why: the whole-path condition is one local condition per instant. The global phrasing is valuable because it is coordinate-free (Theorem 3.5), not because it adds physics; only in quantum mechanics (§3.10) does the whole-path view carry physics of its own. If you catch yourself asking how the particle knows where it will end up, reread (2.2).
- **Stationary, not least.** The oscillator past half a period is a saddle (Theorem 3.6); the Kepler circle over a full period has a flat direction (Exercise C2). Never argue "this path has lower action, so it is the true one." Say "stationary", and when minimality matters, compute (3.10).
- **$\delta$ versus $d$.** $\varepsilon\eta(t)$ compares two paths at the same instant; $d/dt$ moves along one path. Take $\partial/\partial\dot q$ with $q$, $\dot q$ as independent slots first, then substitute the path and take $d/dt$. Most wrong E–L equations apply $d/dt$ too early. The project's `total_d` makes the order mechanical.
- **The covariance proof looks like bookkeeping.** It is the most important proof here, because it licenses everything from CM-3 on. Remember "dots cancel", "$d/dt$ commutes with $\partial/\partial Q$", and the cancellation in step 4; do Exercise B2 by hand once and watch $J^\top$ appear.
- **$p_i \ne m\dot q_i$.** $p_\theta = mr^2\dot\theta$ has units of angular momentum; $m\mathbf v + q_e\mathbf A$ is not even gauge-invariant. Compute $p_i = \partial L/\partial\dot q_i$ every time and check units.
- **Endpoint derivatives (§3.8).** Moving $q_2$ changes the whole true path; only the boundary term of (3.2) survives because the interior change multiplies $\mathcal E = 0$. Getting $\partial S/\partial q_2 = 0$ means you dropped that term.
- **$e^{iS/\hbar}$ is a preview.** §3.10 derives nothing; do not build on it before QM-9.

## 9. Check yourself

1. State Hamilton's principle precisely: what is held fixed, and what is varied?
2. Derive $m\ddot z = -\partial V/\partial z$ from $L = \frac12 m|\dot{\mathbf r}|^2 - V(\mathbf r)$ in two lines.
3. Why is the Lagrangian $T - V$ and not $T + V$? Give the one-line computation.
4. Show in one line, using the action, that $L$ and $L + dF(q, t)/dt$ give the same motion.
5. State Theorem 3.5 in words. Which matrix relates $\tilde{\boldsymbol{\mathcal E}}$ and $\boldsymbol{\mathcal E}$?
6. Write $L$ for a particle in polar coordinates with $V(r)$. Which momentum is conserved, and what are its units?
7. Why is the free particle's true path an exact global minimum of $S$?
8. For an oscillator of period $T$, for which interval lengths is the true path a minimum? Name a descending direction past the threshold.
9. What are $\partial S_{\text{cl}}/\partial q_2$ and $\partial S_{\text{cl}}/\partial t_2$?
10. Write the Lagrangian of a charged particle and its canonical momentum. Why is that momentum not $m\mathbf v$?
11. What does a gauge transformation do to $L$, and why does the motion not change?
12. Why does the discrete gradient (2.2) show that stationary action is not teleological?

## 10. Pair with

- **Taylor, *Classical Mechanics*, ch. 6 and §7.1.** Ch. 6 is F6 at a gentler pace; §7.1 does Hamilton's principle for unconstrained motion with many worked problems. Use it for practice volume.
- **Feynman, *Lectures on Physics* Vol. II, ch. 19 ("The Principle of Least Action").** A transcribed lecture with the best intuition anywhere for why $T - V$, how the classical path emerges from quantum interference, and what "least" means. Read it after §2.
- **Landau & Lifshitz, *Mechanics*, §1–2 (and §3–4).** Terse and complete. §2 states the principle; §3–4 derive the free-particle Lagrangian from Galilean relativity, which this module only cites (Exercise B1).
- **Susskind & Hrabovsky, *The Theoretical Minimum*, lectures 6–7.** A conversational pass through action, Lagrangians, and the first conservation laws; good for a second reading once the formalism is familiar.
