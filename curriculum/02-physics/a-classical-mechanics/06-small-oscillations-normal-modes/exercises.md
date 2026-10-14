# CM-6 Exercises

Conventions as in the notes: $q$ is measured from equilibrium, $M\ddot q = -Kq$ (3.3), modes solve $K\mathbf{a} = \omega^2M\mathbf{a}$ (3.5), $\langle\mathbf{u}, \mathbf{v}\rangle_M = \mathbf{u}^\top M\mathbf{v}$, modes are $M$-normalized ($\mathbf{a}_i^\top M\mathbf{a}_j = \delta_{ij}$) and ordered $\omega_1 \le \omega_2 \le \cdots$. Every number the solutions quote is reproduced in Part 9 of `project.py`.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Wall–$k_1$–$m_1$–$k_2$–$m_2$–$k_3$–wall on a line, with $m_1 = 1$, $m_2 = 2$, $k_1 = 1$, $k_2 = k_3 = 2$ (consistent units). Write $M$ and $K$. Solve $\det(K - \omega^2M) = 0$ by hand, find both mode shapes, and $M$-normalize them. Check $\mathbf{a}_1^\top M\mathbf{a}_2 = 0$, then compute the ordinary dot product $\mathbf{a}_1^\top\mathbf{a}_2$ of the normalized shapes.

**A2.** For the double pendulum of Example 4.4 with $m = \ell = 1$, $g = 9.81$: $M$-normalize the shapes $(1, \sqrt2)^\top$ and $(1, -\sqrt2)^\top$ (five significant figures), assemble $A$, and check $A^\top MA = I$ and $A^\top KA = \operatorname{diag}(\omega_1^2, \omega_2^2)$. Then show that $AA^\top = M^{-1}$ follows from $A^\top MA = I$, and check it numerically.

**A3.** (When $\langle\cdot,\cdot\rangle_M$ is not an inner product.)
(a) $M = \begin{pmatrix}1 & 2\\2 & 1\end{pmatrix}$ is symmetric. Find $\mathbf{v} \ne 0$ with $\mathbf{v}^\top M\mathbf{v} < 0$. Which axiom fails?
(b) $M = \begin{pmatrix}1 & 1\\0 & 1\end{pmatrix}$ has $\mathbf{v}^\top M\mathbf{v} > 0$ for every $\mathbf{v} \ne 0$. Prove that, then find the axiom it fails.
(c) Could either matrix be the mass matrix of a mechanical system? Answer using Derivation 1.

**A4.** A particle of unit mass moves in the plane in the potential $V(x, y) = \frac14(x^2 - 1)^2 + (y - x)^2$. Find every equilibrium. At each one compute the Hessian $K$ and decide stability with Theorem 3.5. Give the two normal-mode frequencies at a stable equilibrium, and the exponential growth rate at an unstable one.

**A5.** For $N = 3$ beads on a string (Example 4.5), check by direct multiplication that $(1/\sqrt2, 1, 1/\sqrt2)^\top$, $(1, 0, -1)^\top$ and $(1/\sqrt2, -1, 1/\sqrt2)^\top$ are mode shapes, and find each $\omega^2/\omega_0^2$. Count the interior nodes of each shape and compare with (4.5).

**A6.** You want to build the coupled pendulums of Example 4.1 with $\ell = 0.6$ m and bobs of mass $m = 0.2$ kg, so that the energy moves from one pendulum to the other in 6 s ($T_{\text{beat}} = 12$ s). What spring constant do you need? Compute it exactly from $\omega_2^2 = g/\ell + 2k/m$, and again from the weak-coupling approximation $\omega_2 - \omega_1 \approx k/(m\omega_1)$. How big is the discrepancy?

## B. Core (do all; these are the module)

**B1.** (Damping.)
(a) Prove Theorem 3.6. Write every step, including where (3.11) is used.
(b) Add one dashpot $c = 1$ between mass 1 of Example 4.2 and the wall, so $C = \operatorname{diag}(1, 0)$. Compute $A^\top CA$ with the $M$-normalized modes of Example 4.2 and show it is not diagonal. Show that no $\alpha, \beta$ give $C = \alpha M + \beta K$.
(c) Use Rayleigh damping $\alpha = 0.1$, $\beta = 0.05$ instead. Compute $\zeta_1$ and $\zeta_2$.
(d) Find $\alpha, \beta$ that give both modes $\zeta = 0.02$. Derive the general formulas for equal damping $\zeta$ in two modes $\omega_1, \omega_2$.
(e) In two sentences, say what the off-diagonal entry in (b) does physically.

**B2.** (Double pendulum with unequal masses.) The upper bob has mass $m_1$, the lower $m_2$, and both strings have length $\ell$.
(a) Linearize CM-3's Lagrangian to get $M = \ell^2\begin{pmatrix}m_1 + m_2 & m_2\\m_2 & m_2\end{pmatrix}$ and $K = g\ell\operatorname{diag}(m_1 + m_2,\ m_2)$.
(b) With $\mu = m_2/(m_1 + m_2)$, show that $\omega^2 = \dfrac{g}{\ell}\,\dfrac{1}{1 \pm \sqrt\mu}$ and that the shapes have $a_2/a_1 = \pm1/\sqrt\mu$.
(c) Check the formulas against Example 4.4, then evaluate them for $m_1 = 3$, $m_2 = 1$.
(d) Interpret the limits $m_2 \to 0$ and $m_1 \to 0$ physically.

**B3.** (Why $M^{-1}K$ has real eigenvalues.)
(a) With $M = LL^\top$ and $\tilde K = L^{-1}KL^{-\top}$, show $M^{-1}K = L^{-\top}\tilde KL^\top$. Conclude that $M^{-1}K$ is diagonalizable with real eigenvalues even though it is not symmetric.
(b) For Example 4.2, compute $L$ and $\tilde K$. Find the eigenvalues of $\tilde K$ and its orthonormal eigenvectors $\mathbf{y}_j$, and check that $L^{-\top}\mathbf{y}_j$ reproduces the normalized modes of the notes.
(c) Show that the matrix of the operator $M^{-1}K$ in the basis $\{\mathbf{a}_j\}$ is $\Omega^2$. Explain why `np.linalg.eig(np.linalg.solve(M, K))` returns the right frequencies but eigenvectors you should not use directly in (3.13).

**B4.** (Modal energies.)
(a) Show that $E = \frac12\dot q^\top M\dot q + \frac12q^\top Kq$ equals $\sum_jE_j$ with $E_j = \frac12(\dot\xi_j^2 + \omega_j^2\xi_j^2)$, and that each $E_j$ is separately constant.
(b) For Example 4.1 released from $(\theta_0, 0)$ at rest, compute $E_1$ and $E_2$ and check that their sum is $\frac12q_0^\top Kq_0$.
(c) At $t = T_{\text{transfer}}$ show $\theta_1 = 0$ and $|\dot\theta_1| \le \varepsilon\theta_0$. Pendulum 1 has almost no energy then, yet by (a) nothing has changed. Reconcile the two statements.

**B5.** (Discrete sine orthogonality.)
(a) Prove $\displaystyle\sum_{k=1}^N\sin\frac{jk\pi}{N+1}\sin\frac{j'k\pi}{N+1} = \frac{N+1}{2}\,\delta_{jj'}$ for $1 \le j, j' \le N$. Do the case $j \ne j'$ two ways: from Theorem 3.2 with no trigonometry, and by direct summation.
(b) Write the $M$-normalized modes of the loaded string ($M = mI$).
(c) For $N = 5$, expand $q_0 = (1, 2, 3, 2, 1)^\top$ in the sine modes. Explain without computing why the even-$j$ coefficients vanish.
(d) Write $y_k(t)$ for the string released from rest in the shape $q_0$.

**B6.** (Zero modes and momentum.) In the CO₂ model of Example 4.3 ($m = 16$, $M_c = 12$, $k = 1$), all atoms start at equilibrium, and oxygen 1 is struck so that it starts with velocity 1 while the others are at rest.
(a) Compute $\xi_j(0)$ and $\dot\xi_j(0)$ for the three normalized modes.
(b) Write $x_1(t)$, $x_2(t)$ and $x_3(t)$ from (3.13).
(c) Show that the center of mass moves at speed $4/11$, and that this motion is exactly the zero-mode term.
(d) Prove for any free chain with diagonal $M$ and zero mode $\propto(1, \dots, 1)^\top$: the zero-mode velocity $\dot\xi_0$ is proportional to the total momentum and is conserved, and every vibrational mode carries zero total momentum.
(e) What does code that evaluates $\sin(\omega_jt)/\omega_j$ naively do on this problem? How does `fit_initial_conditions` avoid it?

**B7.** (Extend `project.py`: when does the linear picture fail?)
(a) From the Lagrangian of Example 4.4 ($m = \ell = 1$), derive the exact equations of motion of the double pendulum in the form $\begin{pmatrix}2 & \cos\Delta\\\cos\Delta & 1\end{pmatrix}\ddot\theta = \mathbf{f}(\theta, \dot\theta)$ with $\Delta = \theta_1 - \theta_2$, and code them as `double_pendulum_rhs`.
(b) Release from $(\theta_0, 0)$ at rest for $\theta_0 = 0.1$ and $1.0$ rad. Integrate for 10 s with `solve_ivp` and compare with `fit_initial_conditions`. Report $\max_t|q_{\text{exact}} - q_{\text{modal}}|/\theta_0$.
(c) Project the exact motion onto the modes, $\xi = A^\top Mq$, $\dot\xi = A^\top M\dot q$, and track $E_1(t)$ and $E_2(t)$ from B4. Report the relative swing $(\max E_1 - \min E_1)/E_1(0)$.
(d) By hand, show that this initial condition gives $E_1 = E_2 = \frac12g\theta_0^2$ exactly.
(e) At 0.1 rad the position error is large while $E_1$ barely moves. Explain why the two diagnostics disagree, and which one tells you whether the mode picture is valid.

**B8.** (Receptance and nodes.) Take the chain of project Part 7: $M = I$, $K = \operatorname{tridiag}(-1, 2, -1)$ ($3\times3$), $C = 0.02M + 0.01K$.
(a) Show $K^{-1} = \sum_j\mathbf{a}_j\mathbf{a}_j^\top/\omega_j^2 = A\Omega^{-2}A^\top$ in two ways: from (3.16) at $\omega = 0$, and directly from (3.11).
(b) Compute the static displacement under a unit force on mass 1, both directly and as a modal sum, and give each mode's contribution. Which mode dominates, and why?
(c) Show that mode 2 has a node at mass 2. Use (3.16) to explain why a force on mass 2 gives two resonance peaks in every response, and why a force on mass 1 gives two peaks in the response of mass 2.
(d) For a force on mass 1, estimate $|\hat q_1(\omega_j)|$ at each resonance from the $j$-th term of (3.16) alone, and compare with the full sum computed in the script.

## C. Challenge (do at least one)

**C1.** (The discrete drum.) Put 25 beads on a $5\times5$ grid. Each moves transversely and is joined to its four neighbors (or to the fixed frame) by string under tension. In units $\tau/d = m = 1$, $K = K_1\otimes I_5 + I_5\otimes K_1$, where $K_1 = \operatorname{tridiag}(-1, 2, -1)$ is the $N = 5$ string matrix.
(a) Derive this $K$ from $V = \frac12\sum(\text{segment height differences})^2$ by ordering the beads row by row.
(b) Prove that if $K_1\mathbf{u} = s\mathbf{u}$ and $K_1\mathbf{v} = t\mathbf{v}$, then $K(\mathbf{u}\otimes\mathbf{v}) = (s + t)(\mathbf{u}\otimes\mathbf{v})$. Deduce all 25 frequencies, $\omega_{jl}^2 = s_j + s_l$ with $s_j = 4\sin^2(j\pi/12)$.
(c) Show that the 25 product vectors are orthonormal, so these are all the modes.
(d) Count the distinct frequencies and account for every degeneracy. Which ones come from the symmetry of the square, and which do not?
(e) Do "Your turn" 1 in `project.py`.

**C2.** (A four-degree-of-freedom car.) The body has mass $m_b = 1200$ kg and pitch moment of inertia $I = 1800\ \mathrm{kg\,m^2}$. Its center of mass is $a = 1.1$ m behind the front axle and $b = 1.5$ m ahead of the rear axle. Use body heave $z$ (up), pitch $\theta$ (nose up) and axle heights $z_f$, $z_r$ as coordinates. Suspension springs: $k_{sf} = 45$ kN/m and $k_{sr} = 35$ kN/m (each lumps both sides of its axle). Unsprung masses: 90 kg front, 70 kg rear. Tyre stiffness: 200 kN/m per axle. The road is flat at zero.
(a) Write the four spring elongations as $e = Bq$, then form $K = B^\top CB$ and $M$. Which entry couples bounce to pitch, and what condition makes it vanish?
(b) Compute the four modes with `normal_modes` and give the frequencies in Hz.
(c) Compare with one-degree-of-freedom estimates: bounce and pitch with the suspension and tyre springs in series, and each axle hopping between its suspension spring and tyre with the body held fixed.
(d) For the two body modes, find the point along the body that does not move, and name each mode. The coupling entry is small, yet the shapes mix bounce and pitch heavily. Why?

**C3.** (The Caughey series.)
(a) For any polynomial $p$, show that $C = M\,p(M^{-1}K)$ is symmetric and that $A^\top CA = p(\Omega^2)$.
(b) Show that a polynomial of degree $n - 1$ can give each of $n$ modes any damping ratio you choose, and that Rayleigh damping is the degree-1 case.
(c) For the chain of B8, find $c_0, c_1, c_2$ with $p(x) = c_0 + c_1x + c_2x^2$ so that $\zeta = (0.02, 0.03, 0.05)$, and verify with the script.
(d) Give two reasons engineers rarely go beyond degree 1.
