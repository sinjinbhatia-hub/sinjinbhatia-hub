# CM-2 Exercises

Conventions as in the notes: $L = T - V$; $\mathcal E_i = \frac{d}{dt}\frac{\partial L}{\partial\dot q_i} - \frac{\partial L}{\partial q_i}$; $p_i = \partial L/\partial\dot q_i$; $\mathcal F_i = \partial L/\partial q_i$; vertical coordinate up, $V = mgy$; charge $q_e$.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** A free particle of mass $m = 2$ must be at $x = 0$ at $t = 0$ and at $x = 3$ at $t = 2$.
(a) Compute $S$ along the true path.
(b) Compute $S$ along $x(t) = 3t^2/4$, which has the same endpoints.
(c) Write the second path as $x_{\text{cl}} + \eta$ and verify $S = S_{\text{cl}} + \frac12 m\int_0^2\dot\eta^2\,dt$. Why is there no cross term?

**A2.** (CM-1's Atwood machine in one line.) Masses $m_1$ and $m_2$ hang from a massless string over a massless, frictionless pulley; heights $y_1$, $y_2$ are measured upward and $y_1 + y_2$ is constant. Using the single coordinate $y = y_1$, write $T$, $V$ and $L$, and derive $\ddot y$. Compare with CM-1 (3.16). Where did the tension go?

**A3.** (The wrong sign.) For a mass on a spring, $V = \frac12 kx^2$. Derive the E–L equation from $L = T + V$, solve it with $x(0) = x_0$, $\dot x(0) = 0$, and describe the motion. Then show that $L = 3(T - V)$ gives the correct motion.

**A4.** To $L = \frac12 m\dot q^2$ add each term below. Which ones change the equation of motion, and to what? For each that does not, find $F(q, t)$ with the added term equal to $dF/dt$.
(a) $3q\dot q$  (b) $t\dot q$  (c) $\dot q\cos q$  (d) $q + t\dot q$

**A5.** A particle in spherical coordinates has $L = \frac12 m(\dot r^2 + r^2\dot\theta^2 + r^2\sin^2\theta\,\dot\phi^2) - V(r)$. Compute $p_r, p_\theta, p_\phi$ and $\mathcal F_r, \mathcal F_\theta, \mathcal F_\phi$. Which momentum is conserved? Give the SI units of each $p$.

## B. Core (do all)

**B1.** (Why not some other kinetic term?)
(a) Show that every $L = f(\dot q)$ with $f'' \ne 0$ gives $\ddot q = 0$.
(b) Now take $L = f(\dot q) - V(q)$. Show the E–L equation is $f''(\dot q)\,\ddot q = -V'(q)$, and conclude that it agrees with Newton's $m\ddot q = -V'(q)$ for every $V$ and every initial condition only if $f(v) = \frac12 mv^2 + av + b$.
(c) What do the terms $av$ and $b$ do, and why?
(d) Part (a) says the free particle alone cannot pin down $f$. In two or three sentences, what extra physical requirement do Landau & Lifshitz (§4) use, and how does $\frac12 mv^2$ satisfy it while $v^4$ does not?

**B2.** (Covariance by hand.) For a particle in the plane with potential $V(x, y)$, let $\mathcal E_x, \mathcal E_y$ be the Cartesian E–L expressions and $\tilde{\mathcal E}_r, \tilde{\mathcal E}_\theta$ those of $\tilde L(r, \theta, \dot r, \dot\theta)$ obtained by substituting $x = r\cos\theta$, $y = r\sin\theta$.
(a) Compute $\tilde{\mathcal E}_r$ and $\tilde{\mathcal E}_\theta$ directly, using $\tilde V(r, \theta) := V(r\cos\theta, r\sin\theta)$, $\partial\tilde V/\partial r = \cos\theta\,V_x + \sin\theta\,V_y$ and $\partial\tilde V/\partial\theta = -r\sin\theta\,V_x + r\cos\theta\,V_y$.
(b) Express $\ddot x$ and $\ddot y$ in terms of $r, \theta$ and their derivatives.
(c) Verify (3.8): $\tilde{\mathcal E}_r = \cos\theta\,\mathcal E_x + \sin\theta\,\mathcal E_y$ and $\tilde{\mathcal E}_\theta = -r\sin\theta\,\mathcal E_x + r\cos\theta\,\mathcal E_y$.

**B3.** (A time-dependent point transformation.) A free particle moves in a plane, $L = \frac12 m(\dot x^2 + \dot y^2)$. Describe it in axes rotating at constant rate $\Omega$: $x = X\cos\Omega t - Y\sin\Omega t$, $y = X\sin\Omega t + Y\cos\Omega t$.
(a) Show $\tilde L = \frac12 m\big[(\dot X - \Omega Y)^2 + (\dot Y + \Omega X)^2\big]$.
(b) Derive the E–L equations and identify the centrifugal and Coriolis terms.
(c) Compute $p_X$. Why is it not $m\dot X$? Relate it to the inertial momentum.
(d) Which theorem guarantees, without solving anything, that the solutions of (b) are exactly the inertial motions (straight lines at constant speed in the $(x, y)$ plane) written in rotating coordinates? (CM-3 Example 4.9 does the general 3D case.)

**B4.** (Counting descending directions.) For the oscillator with $m = k = 1$ and endpoints fixed at $t = 0$ and $t = \tau$ (not a multiple of $\pi$):
(a) Using (3.11), show that the number of independent descending directions equals the number of points $\pi, 2\pi, 3\pi, \dots$ strictly inside $(0, \tau)$.
(b) Give the count for $\tau = 2, 4, 5, 7$ and compare with the project's discrete Hessian.
(c) At $\tau = \pi$ exactly, find all solutions of $\ddot q = -q$ with $q(0) = 0$, $q(\pi) = 0$, and compute their actions. Relate this family to the vanishing $n = 1$ bracket in (3.11).

**B5.** (The oscillator's principal function.) For $L = \frac12 m\dot q^2 - \frac12 m\omega^2q^2$ with $\sin\omega\tau \ne 0$, the true path from $q_1$ at $t = 0$ to $q_2$ at $t = \tau$ is $q(t) = [q_1\sin\omega(\tau - t) + q_2\sin\omega t]/\sin\omega\tau$.
(a) Show that $S_{\text{cl}} = \dfrac{m\omega\big[(q_1^2 + q_2^2)\cos\omega\tau - 2q_1q_2\big]}{2\sin\omega\tau}$. (Hint: integrate $\dot q^2$ by parts.)
(b) Verify $\partial S_{\text{cl}}/\partial q_2 = m\dot q(\tau)$ and $\partial S_{\text{cl}}/\partial q_1 = -m\dot q(0)$, as (3.12) predicts.
(c) Show $\partial S_{\text{cl}}/\partial\tau = -E$. For $m = \omega = 1$, $q_1 = 0$, $q_2 = 1$, $\tau = 2$, give $S_{\text{cl}}$, $p(0)$, $p(\tau)$ and $\partial S_{\text{cl}}/\partial\tau$.
(d) What happens to $S_{\text{cl}}$ as $\tau \to \pi/\omega$? Connect this with Theorem 3.6.

**B6.** (Gauge.)
(a) Show that $\mathbf A \to \mathbf A + \nabla\chi$, $\Phi \to \Phi - \partial_t\chi$ leaves $\mathbf E$ and $\mathbf B$ unchanged and changes (3.15) by $q_e\,d\chi/dt$.
(b) For uniform $\mathbf B = B\hat{\mathbf z}$, show that the symmetric gauge $\mathbf A_s = \frac B2(-y, x, 0)$ and the Landau gauge $\mathbf A_L = (0, Bx, 0)$ differ by a gradient, and find $\chi$.
(c) Compute $p_x$ and $p_y$ in both gauges. Which are conserved in which gauge, and why? Does the motion depend on the gauge?

**B7.** (Code.) In `project.py`, take the oscillator with $m = k = 1$, $q(0) = 0$, $q(2) = 1$, and add $F(q, t) = q^2\sin t$: $L' = L + dF/dt$.
(a) Write out $dF/dt$.
(b) Predict, then compute, $S' - S$ along three paths with these endpoints: the true path $\sin t/\sin 2$, the line $t/2$, and $t^2/4$. Assert that all three differences agree.
(c) In two or three sentences: why must the minimizer of a discretized $S'$ equal that of $S$ if you discretize the added term as $\sum_k\big[F(q_{k+1}, t_{k+1}) - F(q_k, t_k)\big]$? What changes if you instead discretize $dF/dt = 2q\dot q\sin t + q^2\cos t$ with a left-point rule?

## C. Challenge (do at least one)

**C1.** (The brachistochrone as mechanics.) The descent time along a wire $y(x)$, with $y$ measured *down* from the start (F6 Example 4.2), is $\int L\,dx$ with $L(y, y') = \sqrt{(1 + y'^2)/(2gy)}$. Read $x$ as "time" and $y$ as "coordinate".
(a) Compute the "momentum" $p = \partial L/\partial y'$ and the "energy" $h = py' - L$.
(b) Which of the two is conserved along the optimal curve, and why?
(c) With $y' = \tan\alpha$ ($\alpha$ the wire's angle below horizontal) and $v = \sqrt{2gy}$, show the conserved quantity says $\sin\beta/v = \text{const}$, where $\beta$ is the angle from the vertical. Which law of optics is this?
(d) Derive $y(1 + y'^2) = \text{const}$ and check that the cycloid $x = a(\theta - \sin\theta)$, $y = a(1 - \cos\theta)$ satisfies it.

**C2.** (Kepler: circle against ellipse.) Take $m = k = 1$, $V = -1/r$.
(a) The circular orbit $r = 1$ has $\dot\theta = 1$ (Example 4.2). Compute its action over $t \in [0, 2\pi]$.
(b) Launch from $(1, 0)$ with speed 1 in a direction tilted by $\gamma$ from the circle's velocity. Using $E = -k/(2a)$ and $P = 2\pi a^{3/2}$ (CM-5), show the orbit has $a = 1$, returns to $(1, 0)$ at $t = 2\pi$, and has eccentricity $e = \sin\gamma$ (use $\ell^2 = a(1 - e^2)$ for the angular momentum $\ell$).
(c) Prove the virial relation $2\langle T\rangle = -\langle V\rangle$ over one period by averaging $\frac{d}{dt}(\mathbf r\cdot\mathbf p)$. Deduce that every orbit in (b) has the same action as the circle.
(d) For the non-Keplerian path $r = 1 + \varepsilon\sin(t/2)$, $\theta = t$ (same endpoints, same duration), show $S = 3\pi + \frac{13\pi}8\varepsilon^2 + O(\varepsilon^3)$.
(e) Use (c) and (d) with §3.7: is the circle over one full period a strict minimum of $S$? What is $t = 2\pi$ relative to $t = 0$?

**C3.** (A Lagrangian for a damped oscillator.)
(a) Show that $L = e^{\gamma t}\big(\frac12 m\dot q^2 - \frac12 kq^2\big)$ gives $m\ddot q + m\gamma\dot q + kq = 0$.
(b) Compute $p = \partial L/\partial\dot q$. Is it the physical momentum?
(c) Compute $h = p\dot q - L$. Is it conserved? Is it the mechanical energy? Find $dh/dt$.
(d) Two coordinates are coupled through $V(q_1, q_2)$ with $\partial V/\partial q_1$ and $\partial V/\partial q_2$ both not identically zero, and you want damping rates $\gamma_1 \ne \gamma_2$. Show that no Lagrangian of the form $\frac12 m\big(a_1(t)\dot q_1^2 + a_2(t)\dot q_2^2\big) - b(t)V(q_1, q_2)$ produces $m\ddot q_i + m\gamma_i\dot q_i + \partial V/\partial q_i = 0$. What does this say about §3.11's claim?
