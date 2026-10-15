# F6: Calculus of Variations

**Time:** ~6 hours · **Prerequisites:** multivariable calculus (partial derivatives, chain rule, differentiating under an integral sign, Lagrange multipliers in $\mathbb{R}^n$) and ODEs (separable equations, first integrals); no linear algebra needed · **Unlocks:** CM-2 (and through it CM-3, CM-4); indirectly QM-9 (path integrals) and ML-10 (variational inference)

> **Conventions in this module.** The independent variable is $x$ on a fixed interval $[a, b]$; the unknown is a function $y(x)$, and $y' = dy/dx$. A **Lagrangian** $L(x, y, y')$ is an ordinary function of *three independent slots*; $\partial L/\partial y$ and $\partial L/\partial y'$ are partial derivatives that treat the three slots as independent variables, exactly as $\partial f/\partial v$ treats $u, v, w$ as independent in $f(u, v, w)$. The symbol $\frac{d}{dx}$ is the *total* derivative along a curve: $\frac{d}{dx}F(x, y(x), y'(x)) = F_x + F_y y' + F_{y'} y''$. A **functional** is $J[y] = \int_a^b L(x, y, y')\,dx$; square brackets remind you that the input is a function. A **variation** of $y$ is $y_\varepsilon = y + \varepsilon\eta$ with $\eta(a) = \eta(b) = 0$; $\varepsilon$ is a real number, $\eta$ is a fixed function. The **first variation** is $\delta J := \frac{d}{d\varepsilon} J[y + \varepsilon\eta]\big|_{\varepsilon = 0}$, and physicists write $\delta y := \varepsilon\eta$. "**Stationary**" means $\delta J = 0$ for every admissible $\eta$; it does *not* mean minimum. Subscripts on $L$ are partial derivatives: $L_y = \partial L/\partial y$, $L_{y'} = \partial L/\partial y'$, $L_{y'y'} = \partial^2 L/\partial y'^2$. In mechanics the independent variable is time $t$, the function is $q(t)$, and $L = T - V$. **Sign convention for the brachistochrone only:** $y$ is measured *downward* from the starting point, so $y > 0$ below the start; everywhere else $y$ points up. $g = 9.81\ \mathrm{m/s^2}$.

## 0. Where this sits

You know how to minimize a function of several variables: set the gradient to zero. This module does the same thing when the unknown is not a point in $\mathbb{R}^n$ but a *function*, and the thing being minimized is an integral that depends on that function. The only prerequisites are multivariable calculus and ODEs; the whole module is one idea (the Euler–Lagrange equation) plus the craft of solving the ODEs it produces.

Everything in CM-2 onward is "minimize a functional": Hamilton's principle (CM-2), constraints by multipliers (CM-3), conservation laws from symmetries (CM-4). The same move reappears as Fermat's principle in optics, as the stationary-phase path that dominates Feynman's path integral in QM-9, as the Dirichlet principle behind Laplace's equation in EM-2, and as variational inference in ML-10, where the "function" being optimized is a probability density. Most engineering curricula never teach this module and then assume it. This one spends six hours on it.

## 1. The problem that forced this idea

In June 1696 Johann Bernoulli published a challenge in the *Acta Eruditorum*; when he reissued it in January 1697 he addressed it "to the most acute mathematicians of the whole world." A bead slides without friction along a wire from a point $A$ to a lower point $B$, not directly below $A$. What shape of wire gets the bead from $A$ to $B$ in the least time? He called the curve the **brachistochrone** (Greek: shortest time). Galileo had asked the question in 1638 and guessed a circular arc; he was wrong. Bernoulli knew the answer (a cycloid) and gave the world six months. Leibniz solved it quickly and asked for the deadline to be extended so that others could try. Solutions came from Jakob Bernoulli (Johann's brother), Leibniz, Tschirnhaus, l'Hôpital, and an anonymous one from England that Johann recognized at once, "as the lion by his claw." It was Newton's. The May 1697 issue of the *Acta* printed several of them alongside Johann's own. According to an account from Newton's niece, he came home from the Mint at four in the afternoon, very tired, found the problem waiting, and did not go to bed until he had solved it, at four the next morning.

Why was this hard for the people who had just invented calculus? Because calculus finds the best **number**. You write $f(x)$, set $f'(x) = 0$, and solve. The brachistochrone asks for the best **function**: among all curves $y(x)$ joining $A$ to $B$, which one makes the descent time $T[y]$ smallest? There are infinitely many directions in which to perturb a curve, and no obvious meaning for "$dT/dy$."

The problem had precedents. Hero of Alexandria (first century) had shown that reflected light takes the shortest path. Fermat, in 1662, argued that refracted light takes the path of *least time*, and derived Snell's law from it; nobody could say what physical mechanism made light "choose" the fastest route, and Descartes's followers attacked the idea as teleological. Jakob Bernoulli, replying to his brother in 1697, posed the harder **isoperimetric** problems (fixed perimeter, maximal area; the version in Virgil's *Aeneid*, where Dido buys as much land as an ox-hide can enclose, is where the name "Dido's problem" comes from), which need a constraint handled inside the minimization. Each of these was solved by a trick specific to the problem.

Euler turned the tricks into a method. In his 1744 book *Methodus inveniendi lineas curvas maximi minimive proprietate gaudentes* ("a method for finding curves enjoying a maximum or minimum property") he replaced the curve by a polygon with many vertices, moved one vertex at a time, and took the limit as the number of vertices went to infinity. Out came a differential equation that every minimizing curve must satisfy. That is exactly the route §2 and §3.5 take, and it is the reason the project's numerical minimizer works.

In August 1755 a 19-year-old from Turin, Lagrange, wrote to Euler with a cleaner derivation. Instead of polygons, he perturbed the whole curve at once by an arbitrary small function and introduced the symbol $\delta$ for this new kind of differential. The derivation became a few lines of algebra. Euler adopted the method, named the subject the *calculus of variations*, and, by the traditional account, held back his own write-up so that Lagrange could publish first. In 1834–35 Hamilton showed that the whole of Newtonian mechanics is one such problem: the trajectory a system actually follows is the one that makes a certain integral (the action) stationary. Fermat's "least time" for light and Newton's $F = ma$ for matter turned out to be the same kind of statement, and that is what CM-2 is about.

## 2. Intuition

### A function is a point; a functional is a landscape

Take the set of all smooth curves $y(x)$ on $[a, b]$ with $y(a) = A$ and $y(b) = B$. Think of each curve as a single **point** in an enormous space. A functional $J$ assigns a number to each point. For instance, arc length
$$
J[y] = \int_a^b \sqrt{1 + y'(x)^2}\,dx
$$
assigns to each curve its length. So $J$ is a **landscape** over function space: hills and valleys, with each location a whole curve. Minimizing $J$ means finding the lowest point of the landscape.

In $\mathbb{R}^n$ you find the bottom of a valley by checking that every directional derivative is zero: for every direction $\mathbf{h}$, $\frac{d}{d\varepsilon} f(\mathbf{x} + \varepsilon\mathbf{h})\big|_{\varepsilon=0} = \nabla f \cdot \mathbf{h} = 0$. In function space a **direction** is a function $\eta(x)$, and moving in that direction means replacing $y$ by $y + \varepsilon\eta$:

```
   y + εη  ___
          /   \__           original y
    ____/        \___   ___/
   A                 \_/        B
   |  <-- η is a bump that is zero at both ends -->  |
```

Because the endpoints are pinned, $\eta(a) = \eta(b) = 0$. The directional derivative of $J$ in direction $\eta$ is $\frac{d}{d\varepsilon} J[y + \varepsilon\eta]$ at $\varepsilon = 0$. At the bottom of the valley this is zero *for every* $\eta$. The Euler–Lagrange equation of §3 is "gradient $= 0$" written out for this landscape. The project plots $J[y + \varepsilon\eta]$ against $\varepsilon$: for the true minimizer the curve is flat at $\varepsilon = 0$; for any other $y$ it is sloped.

### Euler's picture: discretize, then let $N \to \infty$

Here is a way to make "gradient $= 0$" literal. Replace the unknown curve by its values $y_1, \dots, y_N$ at equally spaced interior points $x_i = a + ih$, $h = (b - a)/(N + 1)$, with $y_0 = A$ and $y_{N+1} = B$ fixed. Connect the dots with straight segments. The functional becomes an ordinary function of $N$ variables:
$$
J_h(y_1, \dots, y_N) = \sum_{i=0}^{N} L\!\left(x_i,\ y_i,\ \frac{y_{i+1} - y_i}{h}\right) h .
$$
Now minimize it the ordinary way: set all $N$ partial derivatives to zero.

**One interior point.** Arc length from $(0, 0)$ to $(2, 2)$ through a single movable point $(1, y_1)$:
$$
J_h(y_1) = \sqrt{1 + y_1^2} + \sqrt{1 + (2 - y_1)^2}, \qquad
\frac{dJ_h}{dy_1} = \frac{y_1}{\sqrt{1 + y_1^2}} - \frac{2 - y_1}{\sqrt{1 + (2 - y_1)^2}} .
$$
Setting this to zero says $\frac{p_0}{\sqrt{1 + p_0^2}} = \frac{p_1}{\sqrt{1 + p_1^2}}$ where $p_0 = y_1$ and $p_1 = 2 - y_1$ are the slopes of the two segments. The function $p \mapsto p/\sqrt{1 + p^2}$ is strictly increasing, so $p_0 = p_1$, so $y_1 = 1$. The slope does not change at the vertex.

**Many interior points.** The same computation at vertex $i$ gives $\frac{p_{i-1}}{\sqrt{1 + p_{i-1}^2}} = \frac{p_i}{\sqrt{1 + p_i^2}}$: the slope does not change at *any* vertex. The polygon is a straight line. Nothing surprising, but look at the *structure* of the condition: $\partial J_h/\partial y_i$ involves only the two segments touching vertex $i$, and it says that a certain quantity ($\partial L/\partial y'$, evaluated on each segment) is the same on both sides. The condition is **local**, and locality is what differential equations express.

**General $L$.** Now take any Lagrangian. The variable $y_i$ appears in two terms of the sum: term $i$, through the slot $y$ and through the slope $(y_{i+1} - y_i)/h$; and term $i - 1$, through the slope $(y_i - y_{i-1})/h$. Differentiating (the chain rule, done in full in §3.5) gives
$$
\frac{\partial J_h}{\partial y_i} = h\left[\ \frac{\partial L}{\partial y}\Big|_{\text{node } i} \;-\; \frac{\left.\frac{\partial L}{\partial y'}\right|_{\text{segment } i} - \left.\frac{\partial L}{\partial y'}\right|_{\text{segment } i - 1}}{h}\ \right].
$$
The second term inside the bracket is a difference quotient: the change in $\partial L/\partial y'$ from one segment to the next, divided by $h$. As $N \to \infty$ it becomes a derivative. So the guess for the continuum condition is
$$
\frac{\partial L}{\partial y} - \frac{d}{dx}\left(\frac{\partial L}{\partial y'}\right) = 0 .
$$
That is the Euler–Lagrange equation. You have now guessed the main theorem of the module. Two things are worth noticing before we prove it. First, "gradient $= 0$" in $N$ variables became a single *differential equation* in the limit, because each component of the gradient only talks to its neighbors. Second, each gradient component carries a factor $h$: the gradient of the discretized functional is $h$ times the Euler–Lagrange expression (§3.5 makes this exact, with the error term). That factor is why the project's optimizer and the ODE solver agree.

### What to expect from the answer

For arc length, $\partial L/\partial y = 0$ and $\partial L/\partial y' = y'/\sqrt{1 + y'^2}$, so the equation says $y'/\sqrt{1 + y'^2}$ is constant, so $y'$ is constant: a straight line. Good. For the brachistochrone the equation is a nonlinear second-order ODE whose solution is a cycloid. For the hanging chain it is $\cosh$. In each case the calculus of variations does not hand you the curve; it hands you an ODE, and you still have to solve it. Hence the ODE prerequisite.

## 3. Formal development

### 3.1 Functionals

**Definition 3.1.** Let $L(x, y, p)$ be a function of three real variables with continuous second partial derivatives. The **functional** with Lagrangian $L$ is
$$
J[y] = \int_a^b L\big(x, y(x), y'(x)\big)\,dx, \tag{3.1}
$$
defined for every continuously differentiable $y$ on $[a, b]$. (We write the third slot as $y'$ rather than $p$ from now on, but remember it is an independent slot.)

Four functionals you will meet:

| Problem | Lagrangian $L(x, y, y')$ | Where |
|---|---|---|
| Arc length of the graph of $y$ | $\sqrt{1 + y'^2}$ | Example 4.1 |
| Descent time of a bead ($y$ down, start at rest) | $\sqrt{1 + y'^2}\big/\sqrt{2gy}$ | Example 4.2 |
| Area of the surface of revolution of $y$ about the $x$-axis | $2\pi\,y\sqrt{1 + y'^2}$ | Example 4.6 |
| Action of a particle of mass $m$ in potential $V$ (independent variable $t$) | $\tfrac12 m\dot q^2 - V(q)$ | §3.13, CM-2 |

Where these come from: arc length is $\int ds$ with $ds = \sqrt{dx^2 + dy^2} = \sqrt{1 + y'^2}\,dx$. Descent time is $\int ds/v$ with $v = \sqrt{2gy}$ from energy conservation. Surface area is $\int 2\pi y\, ds$ (a ring of radius $y$ and width $ds$). The action is new; §3.13 shows that making it stationary gives Newton's second law.

### 3.2 Admissible functions, variations, and the first variation

**Definition 3.2.** The **admissible functions** for the fixed-endpoint problem are the $y \in C^2[a, b]$ with $y(a) = A$ and $y(b) = B$. An **admissible variation** (or test function) is any $\eta \in C^1[a, b]$ with $\eta(a) = \eta(b) = 0$. For admissible $y$ and $\eta$ and real $\varepsilon$, the function $y + \varepsilon\eta$ is a $C^1$ curve with the right endpoints, because $\eta$ vanishes there. (It need not be $C^2$, since $\eta$ is only $C^1$; that is why the competitors in the next definition are taken in $C^1$, where $J$ is defined. The $C^2$ assumption on $y$ itself is used only in Theorem 3.5, and Exercise B3 removes it.)

**Definition 3.3.** $J$ has a **(weak) local minimum** at $y$ if $J[\tilde y] \ge J[y]$ for every $\tilde y \in C^1[a, b]$ with $\tilde y(a) = A$, $\tilde y(b) = B$ and $\max|\tilde y - y| + \max|\tilde y' - y'|$ small enough. $J$ is **stationary** at $y$ if
$$
\delta J := \frac{d}{d\varepsilon} J[y + \varepsilon\eta]\Big|_{\varepsilon = 0} = 0 \quad\text{for every admissible } \eta. \tag{3.2}
$$

If $y$ is a local minimum then for each fixed $\eta$ the ordinary function $\Phi(\varepsilon) = J[y + \varepsilon\eta]$ has a local minimum at $\varepsilon = 0$, so $\Phi'(0) = 0$. So **a minimum is stationary**. The converse is false, and §3.11 says what more is needed.

▶ **Try it first.** Write $\Phi(\varepsilon) = \int_a^b L(x, y + \varepsilon\eta, y' + \varepsilon\eta')\,dx$ and differentiate with respect to $\varepsilon$ under the integral sign using the chain rule. What two terms do you get?

**Derivation (the first variation).**

1. $\Phi(\varepsilon) = \int_a^b L\big(x,\ y(x) + \varepsilon\eta(x),\ y'(x) + \varepsilon\eta'(x)\big)\,dx$. The integrand and its partial derivative with respect to $\varepsilon$ (computed in step 2) are continuous functions of $(x, \varepsilon)$ on the compact set $[a, b] \times [-1, 1]$, because $L$ is $C^2$ and $y, y', \eta, \eta'$ are continuous; so we may differentiate under the integral sign.
2. By the chain rule applied to the second and third slots of $L$ (the first slot, $x$, does not depend on $\varepsilon$):
$$
\Phi'(\varepsilon) = \int_a^b \Big[ L_y\big(x, y + \varepsilon\eta, y' + \varepsilon\eta'\big)\,\eta(x) + L_{y'}\big(x, y + \varepsilon\eta, y' + \varepsilon\eta'\big)\,\eta'(x) \Big]\,dx .
$$
3. Set $\varepsilon = 0$:
$$
\delta J = \int_a^b \Big[ L_y\big(x, y, y'\big)\,\eta + L_{y'}\big(x, y, y'\big)\,\eta' \Big]\,dx . \tag{3.3}
$$
∎

■ **What just happened.** The first variation is the directional derivative of $J$ in direction $\eta$. It is linear in $\eta$, as a directional derivative should be. It contains both $\eta$ and $\eta'$, and the whole art of the next step is to get rid of the $\eta'$.

### 3.3 The fundamental lemma

Before we can conclude anything from "$\int(\ldots)\eta\,dx = 0$ for all $\eta$," we need the one genuinely rigorous step of the module.

**Theorem 3.4 (Fundamental lemma of the calculus of variations).** Let $f$ be continuous on $[a, b]$. If
$$
\int_a^b f(x)\,\eta(x)\,dx = 0
$$
for every $\eta \in C^1[a, b]$ with $\eta(a) = \eta(b) = 0$, then $f(x) = 0$ for every $x \in [a, b]$.

▶ **Try it first.** Suppose $f(x_0) > 0$ at some interior point. Continuity gives you a small interval around $x_0$ on which $f$ stays positive. Build an $\eta$ that is positive on that interval and zero everywhere else, and compute the sign of $\int f\eta$.

**Proof.**

1. Suppose the conclusion fails: $f(x_0) \ne 0$ for some $x_0$. If $x_0$ is an endpoint, continuity gives a nearby interior point where $f \ne 0$, so we may take $x_0 \in (a, b)$. Replacing $f$ by $-f$ if needed, assume $f(x_0) > 0$.
2. By continuity there is a closed interval $[c, d] \subset (a, b)$ containing $x_0$ on which $f(x) > \tfrac12 f(x_0)$.
3. Define the **bump function**
$$
\eta(x) = \begin{cases} (x - c)^2 (d - x)^2, & c \le x \le d, \\ 0, & \text{otherwise.} \end{cases}
$$
On $[c, d]$ it is a polynomial; at $x = c$ and $x = d$ both $\eta$ and $\eta'$ are zero (each factor is squared), so $\eta$ is $C^1$ on all of $[a, b]$. It vanishes at $a$ and $b$ because $[c, d]$ is interior. And $\eta > 0$ on $(c, d)$. So $\eta$ is admissible.
4. Then
$$
\int_a^b f\eta\,dx = \int_c^d f\eta\,dx \ \ge\ \tfrac12 f(x_0) \int_c^d \eta\,dx \ >\ 0,
$$
because the integrand is at least $\tfrac12 f(x_0)\eta \ge 0$ on $[c, d]$ and $\eta$ is strictly positive on the open interval. This contradicts the hypothesis. ∎

■ **What just happened.** "Orthogonal to every test function" forces "zero." The bump is the whole trick: it tests $f$ on an interval so small that $f$ cannot change sign there. If you need test functions that are infinitely differentiable (field theory does), use $\exp\!\big(-1/((x - c)(d - x))\big)$ on $(c, d)$ instead of the polynomial; the proof does not change. This lemma is what most physics courses skip, and it is where the rigor of the Euler–Lagrange equation lives.

### 3.4 The Euler–Lagrange equation

**Theorem 3.5 (Euler–Lagrange).** Let $L$ be $C^2$ and let $y \in C^2[a, b]$ with $y(a) = A$, $y(b) = B$. If $J$ is stationary at $y$, then for every $x \in [a, b]$,
$$
\frac{\partial L}{\partial y}\big(x, y(x), y'(x)\big) - \frac{d}{dx}\left[\frac{\partial L}{\partial y'}\big(x, y(x), y'(x)\big)\right] = 0 . \tag{3.4}
$$

▶ **Try it first.** In (3.3), integrate the $L_{y'}\eta'$ term by parts. What happens to the boundary term, and why?

**Derivation.**

1. Start from (3.3): $\delta J = \int_a^b L_y\,\eta\,dx + \int_a^b L_{y'}\,\eta'\,dx = 0$ for every admissible $\eta$.
2. Integrate the second term by parts with $u = L_{y'}(x, y(x), y'(x))$ and $dv = \eta'\,dx$. Since $L$ is $C^2$ and $y$ is $C^2$, $u$ is $C^1$ in $x$, so this is legitimate:
$$
\int_a^b L_{y'}\,\eta'\,dx = \Big[L_{y'}\,\eta\Big]_a^b - \int_a^b \frac{d}{dx}\big(L_{y'}\big)\,\eta\,dx .
$$
3. The boundary term is $L_{y'}(b)\,\eta(b) - L_{y'}(a)\,\eta(a) = 0$ because $\eta(a) = \eta(b) = 0$. (This is the only place the endpoint conditions on $\eta$ are used, and §3.9 is about what happens when they are not available.)
4. Therefore
$$
\delta J = \int_a^b \left[ L_y - \frac{d}{dx} L_{y'} \right] \eta\,dx = 0 \quad\text{for every admissible } \eta. \tag{3.5}
$$
5. The bracket is a continuous function of $x$ (it involves $L$'s second derivatives and $y''$, all continuous). By the fundamental lemma it is identically zero. That is (3.4). ∎

■ **What just happened.** Integration by parts moved the derivative off $\eta$ and onto $L_{y'}$, at the price of a boundary term that the pinned endpoints killed. Then the lemma converted "zero for every $\eta$" into "zero pointwise." Compare the $\mathbb{R}^n$ statement: $\nabla f \cdot \mathbf{h} = 0$ for all $\mathbf{h}$ forces $\nabla f = \mathbf{0}$; the lemma is the function-space version of "pick $\mathbf{h} = \mathbf{e}_i$."

Written out with the chain rule, (3.4) is a second-order ODE for $y$:
$$
L_{y'y'}\,y'' + L_{y'y}\,y' + L_{y'x} - L_y = 0 . \tag{3.6}
$$
Two boundary conditions, $y(a) = A$ and $y(b) = B$, go with a second-order equation. Solutions of (3.4) are called **extremals** (not "minimizers"), for the reason that follows.

**Stationary is not minimum.** Theorem 3.5 is a *necessary* condition for a minimum, and it is the same necessary condition for a maximum or for a saddle. It is $f'(x) = 0$, not $f''(x) > 0$. A great circle on a sphere is an extremal of arc length between any two of its points, including along the long way round, which is certainly not the shortest route. When this module says "the minimizer satisfies E–L," read it as "if a minimizer exists and is $C^2$, it satisfies E–L"; existence is a separate question (Example 4.6 shows a case where no extremal exists), and sufficiency is §3.11.

### 3.5 The discrete derivation, done properly

This section makes §2's picture exact. It explains why the project's numerical minimizer works, and it is the bridge to the gradient in LA-14 and to the optimization in ML-2.

Grid: $x_i = a + ih$ for $i = 0, \dots, N + 1$, with $h = (b - a)/(N + 1)$; unknowns $y_1, \dots, y_N$; fixed $y_0 = A$, $y_{N+1} = B$. Slopes $p_i = (y_{i+1} - y_i)/h$ for $i = 0, \dots, N$. Discretized functional (left Riemann sum with forward-difference slopes):
$$
J_h(y_1, \dots, y_N) = h \sum_{i=0}^{N} L(x_i, y_i, p_i) . \tag{3.7}
$$

▶ **Try it first.** Which terms of the sum contain $y_i$? Differentiate each with respect to $y_i$, remembering that $p_i$ contains $-y_i/h$ and $p_{i-1}$ contains $+y_i/h$.

**Derivation.**

1. $y_i$ appears in exactly two terms: term $i$, where it is the second slot and also enters $p_i = (y_{i+1} - y_i)/h$ with $\partial p_i/\partial y_i = -1/h$; and term $i - 1$, where it enters $p_{i-1} = (y_i - y_{i-1})/h$ with $\partial p_{i-1}/\partial y_i = +1/h$.
2. Chain rule on each:
$$
\frac{\partial J_h}{\partial y_i} = h\left[ L_y(x_i, y_i, p_i) - \frac{1}{h} L_{y'}(x_i, y_i, p_i) \right] + h\left[ \frac{1}{h} L_{y'}(x_{i-1}, y_{i-1}, p_{i-1}) \right].
$$
3. Regroup:
$$
\frac{\partial J_h}{\partial y_i} = h\left[ L_y(x_i, y_i, p_i) - \frac{L_{y'}(x_i, y_i, p_i) - L_{y'}(x_{i-1}, y_{i-1}, p_{i-1})}{h} \right]. \tag{3.8}
$$
4. Now suppose the $y_i$ are samples $y(x_i)$ of a smooth function. Then $p_i = y'(x_i) + \tfrac12 h\,y''(x_i) + O(h^2)$ (forward difference), so $L_y(x_i, y_i, p_i) = L_y(x_i, y, y') + O(h)$ and likewise $L_{y'}(x_i, y_i, p_i) = Q(x_i) + O(h)$ where $Q(x) := L_{y'}(x, y(x), y'(x))$ is smooth. The difference quotient of a smooth function with $O(h)$ perturbations that vary smoothly from node to node is $Q'(x_i) + O(h)$. (The $O(h)$ perturbation is itself a smooth function of $x_i$ times $h$, so its difference quotient is $O(h)$, not $O(1)$.)
5. Therefore the bracket in (3.8) equals $L_y - \frac{d}{dx}L_{y'}$ at $x_i$, up to $O(h)$, and
$$
\frac{\partial J_h}{\partial y_i} = h\left[\frac{\partial L}{\partial y} - \frac{d}{dx}\frac{\partial L}{\partial y'}\right]_{x = x_i} + O(h^2) . \tag{3.9}
$$
∎

■ **What just happened.** The gradient of the discretized functional is $h$ times the Euler–Lagrange expression sampled at the nodes, with an error one order smaller. So (i) the stationary points of $J_h$ converge to solutions of E–L; (ii) a gradient-based minimizer like L-BFGS (ML-2) applied to $J_h$ is, in effect, solving the E–L equation; (iii) the Euler–Lagrange expression *is* the gradient of $J$ with respect to $y$, up to the volume element $h$ that converts a sum into an integral. Section 3.12 gives that last statement a name. The project checks (3.9) numerically and watches the error drop by a factor of four when $h$ is halved. (For Lagrangians of the special form $L = c\,y'^2 + f(x, y)$ with $c$ constant, $L_y = f_y(x_i, y_i)$ does not see the slope and $(L_{y'}|_i - L_{y'}|_{i-1})/h = 2c\,(y_{i+1} - 2y_i + y_{i-1})/h^2$ is a centered second difference, so the $O(h)$ terms in step 4 vanish and the agreement is even better. A cross term such as $y\,y'$ already breaks this. The project uses a generic $L$ so you see the honest rate.)

### 3.6 Three special cases and their first integrals

Each of these turns the second-order ODE (3.6) into something easier. Check them by *looking at the formula for $L$*, not by thinking about whether $y$ depends on $x$ (of course it does; that is not the question).

**Case A: $L$ does not contain $y$.** Then $L_y = 0$ and (3.4) says $\frac{d}{dx} L_{y'} = 0$, so
$$
\frac{\partial L}{\partial y'} = \text{const}. \tag{3.10}
$$
A first-order ODE. Arc length is this case. In mechanics this is conservation of momentum when the Lagrangian does not contain the coordinate (CM-4 calls such a coordinate *cyclic*).

**Case B: $L$ does not contain $y'$.** Then $L_{y'} = 0$ and (3.4) says $L_y(x, y) = 0$: an *algebraic* equation for $y$ at each $x$, no ODE at all, and the boundary conditions generally cannot be imposed. The variational problems in ML-10 and SM-4 are of this type: the unknown is a probability density $q$ and $L$ depends on $q$ but not $q'$, so the stationarity condition is algebraic and solves to $q \propto \exp(\ldots)$ (Exercise C4).

**Case C: $L$ does not contain $x$ explicitly.** This is the useful one.

**Theorem 3.6 (Beltrami identity).** If $L = L(y, y')$ has no explicit $x$-dependence, then along every solution of the Euler–Lagrange equation
$$
L - y'\,\frac{\partial L}{\partial y'} = \text{const}. \tag{3.11}
$$

▶ **Try it first.** Compute $\frac{d}{dx}\big(L - y' L_{y'}\big)$ along a curve with the chain rule, keeping every term. Then substitute the E–L equation.

**Derivation.**

1. Along a curve $y(x)$, with $L = L(x, y, y')$ for the moment general, the total derivative is
$$
\frac{d}{dx}\big(L - y' L_{y'}\big) = \underbrace{L_x + L_y y' + L_{y'} y''}_{dL/dx} - \underbrace{y'' L_{y'}}_{\frac{d}{dx}(y')\cdot L_{y'}} - \underbrace{y' \frac{d}{dx}L_{y'}}_{y' \cdot \frac{d}{dx}(L_{y'})} .
$$
2. The two $y'' L_{y'}$ terms cancel:
$$
\frac{d}{dx}\big(L - y' L_{y'}\big) = L_x + y'\left(L_y - \frac{d}{dx}L_{y'}\right). \tag{3.12}
$$
3. On a solution of E–L the parenthesis is zero, so $\frac{d}{dx}(L - y' L_{y'}) = L_x$. If $L$ has no explicit $x$, $L_x = 0$ and the left side is constant. ∎

■ **What just happened.** Identity (3.12) holds for *every* curve; it is pure chain rule. Only in step 3 do we use that $y$ solves E–L. Two consequences: (a) on a solution, the Beltrami quantity is constant whenever $L_x = 0$, which drops the order of the ODE by one; (b) on a curve that is *not* a solution, the Beltrami quantity is *not* constant (unless $y'$ happens to vanish wherever the E–L residual does not, since by (3.12) its derivative is $y'$ times that residual), so checking its constancy is a test of whether you have the right curve. The project does exactly that test on the brachistochrone. In mechanics, with $x \to t$ and $L = T - V$, the quantity $\dot q\,\partial L/\partial\dot q - L$ is the energy, and (3.11) is conservation of energy when $L$ does not depend on $t$; CM-4 derives this as the time-translation case of Noether's theorem, and CM-8 renames the same quantity the Hamiltonian.

### 3.7 Several dependent variables

Suppose the unknown is a curve in the plane, $(y_1(x), y_2(x))$, or $n$ functions $y_1, \dots, y_n$ of one variable, with
$$
J[y_1, \dots, y_n] = \int_a^b L(x, y_1, \dots, y_n, y_1', \dots, y_n')\,dx .
$$

▶ **Try it first.** Vary only $y_k$: replace $y_k$ by $y_k + \varepsilon\eta_k$ and leave the other $y_j$ alone. What does the argument of §3.4 give?

**Derivation.** With only $y_k$ varied, $\Phi(\varepsilon) = J[y_1, \dots, y_k + \varepsilon\eta_k, \dots, y_n]$ and the chain rule gives $\Phi'(0) = \int_a^b (L_{y_k}\eta_k + L_{y_k'}\eta_k')\,dx$, exactly (3.3) with $y$ replaced by $y_k$. Integrate by parts, kill the boundary term with $\eta_k(a) = \eta_k(b) = 0$, apply the fundamental lemma:
$$
\frac{\partial L}{\partial y_k} - \frac{d}{dx}\frac{\partial L}{\partial y_k'} = 0, \qquad k = 1, \dots, n . \tag{3.13}
$$
One Euler–Lagrange equation per unknown function, a system of $n$ second-order ODEs. ∎

■ **What just happened.** Varying one function at a time turned the $n$-function problem into $n$ copies of §3.4. The equations are coupled only because each $L_{y_k}$ and $L_{y_k'}$ may contain all the $y_j$; the derivation itself never mixes them.

The Beltrami identity generalizes by the same computation as (3.12): if $L$ has no explicit $x$, then $L - \sum_k y_k' L_{y_k'}$ is constant along solutions. (Compute $\frac{d}{dx}$ of it; the $y_k'' L_{y_k'}$ terms cancel pairwise and what remains is $L_x + \sum_k y_k'(\text{E–L}_k)$.) In CM-3 this is a particle with $n$ generalized coordinates; the project's geodesics on a sphere use $n = 2$.

### 3.8 Higher derivatives

If $L = L(x, y, y', y'')$, integrate by parts twice. The admissible $\eta$ must now have $\eta = \eta' = 0$ at both ends (both $y$ and $y'$ are prescribed there), both boundary terms vanish, and the result is the **Euler–Poisson equation**
$$
\frac{\partial L}{\partial y} - \frac{d}{dx}\frac{\partial L}{\partial y'} + \frac{d^2}{dx^2}\frac{\partial L}{\partial y''} = 0 . \tag{3.14}
$$
The sign alternates because each integration by parts contributes a minus. The engineering example: a beam's elastic energy is $\int \tfrac12 EI\,(y'')^2\,dx$ and the work of a distributed load $q(x)$ is $\int q\,y\,dx$; the Euler–Poisson equation of $\tfrac12 EI y''^2 - qy$ is $EI\,y'''' = q$, the Euler–Bernoulli beam equation you know from solid mechanics.

### 3.9 Natural boundary conditions

What if one endpoint is free, say $y(a) = A$ is fixed but $y(b)$ may be anything? Then an admissible $\eta$ must still satisfy $\eta(a) = 0$, but $\eta(b)$ is arbitrary.

▶ **Try it first.** Redo step 3 of §3.4 without $\eta(b) = 0$. Which boundary term survives, and what must be true of it?

**Derivation.**

1. Integration by parts as before gives
$$
\delta J = L_{y'}(b)\,\eta(b) - \underbrace{L_{y'}(a)\,\eta(a)}_{= 0} + \int_a^b \Big[L_y - \tfrac{d}{dx}L_{y'}\Big]\eta\,dx .
$$
2. The variations with $\eta(b) = 0$ are among the admissible ones, and for them $\delta J$ reduces to the integral. So by §3.4 the Euler–Lagrange equation holds on $[a, b]$.
3. With the integral now zero for *every* admissible $\eta$, stationarity requires $L_{y'}(b)\,\eta(b) = 0$ for every value of $\eta(b)$, hence
$$
\frac{\partial L}{\partial y'}\Big|_{x = b} = 0 . \tag{3.15}
$$
∎

■ **What just happened.** When you do not impose a boundary condition, the variational principle imposes one for you; (3.15) is called the **natural boundary condition**. For arc length, $L_{y'} = y'/\sqrt{1 + y'^2} = 0$ means $y'(b) = 0$: the shortest path from a point to a vertical line meets the line at a right angle. For a string, $L_{y'} \propto$ the transverse force, and (3.15) says a free end carries no transverse force. Exercise B7 uses (3.15) to find the fastest wire from a point to a vertical wall.

### 3.10 Constraints

**Isoperimetric constraints.** Minimize $J[y] = \int_a^b L\,dx$ subject to $K[y] = \int_a^b G(x, y, y')\,dx = \ell$ (a fixed-length chain, a fixed-perimeter fence, a normalized probability density). Write $E_J := L_y - \frac{d}{dx}L_{y'}$ and $E_K := G_y - \frac{d}{dx}G_{y'}$ for the two Euler–Lagrange expressions.

**Theorem 3.7 (Lagrange multiplier rule).** Let $y$ be admissible and $C^2$, with $K[y] = \ell$, and suppose $y$ is not itself an extremal of $K$ (that is, $E_K \not\equiv 0$ along $y$). If $J$ is stationary at $y$ among admissible functions satisfying $K = \ell$, then there is a constant $\lambda$ such that $y$ satisfies the Euler–Lagrange equation of $L + \lambda G$:
$$
E_J + \lambda E_K = 0 . \tag{3.16}
$$

▶ **Try it first.** You cannot use a one-parameter variation $y + \varepsilon\eta$: it generally breaks the constraint. Use two directions, $y + \varepsilon_1\eta_1 + \varepsilon_2\eta_2$, so the constraint becomes a single equation in the two numbers $(\varepsilon_1, \varepsilon_2)$, and apply the ordinary Lagrange multiplier rule in $\mathbb{R}^2$.

**Derivation.**

1. Since $E_K \not\equiv 0$, the fundamental lemma (in contrapositive form) gives an admissible $\eta_2$ with $\int_a^b E_K\,\eta_2\,dx \ne 0$. Fix it. Let $\eta_1$ be any admissible variation.
2. Define two ordinary functions of two real variables: $\Phi(\varepsilon_1, \varepsilon_2) = J[y + \varepsilon_1\eta_1 + \varepsilon_2\eta_2]$ and $\Psi(\varepsilon_1, \varepsilon_2) = K[y + \varepsilon_1\eta_1 + \varepsilon_2\eta_2] - \ell$. Then $\Psi(0, 0) = 0$. By the computation of §3.2 and §3.4 (first variation, then integration by parts with vanishing boundary terms), the partial derivatives at the origin are
$$
\Phi_i := \frac{\partial\Phi}{\partial\varepsilon_i}(0,0) = \int_a^b E_J\,\eta_i\,dx, \qquad \Psi_i := \frac{\partial\Psi}{\partial\varepsilon_i}(0,0) = \int_a^b E_K\,\eta_i\,dx, \qquad i = 1, 2 .
$$
3. $\Psi_2 \ne 0$ by step 1. By the implicit function theorem, near the origin the constraint $\Psi(\varepsilon_1, \varepsilon_2) = 0$ is solved by a $C^1$ function $\varepsilon_2 = \varphi(\varepsilon_1)$ with $\varphi(0) = 0$ and $\varphi'(0) = -\Psi_1/\Psi_2$.
4. The functions $y + \varepsilon_1\eta_1 + \varphi(\varepsilon_1)\eta_2$ are admissible and satisfy the constraint for all small $\varepsilon_1$. So $f(\varepsilon_1) := \Phi(\varepsilon_1, \varphi(\varepsilon_1))$ is stationary at $0$:
$$
0 = f'(0) = \Phi_1 + \Phi_2\,\varphi'(0) = \Phi_1 - \frac{\Phi_2}{\Psi_2}\Psi_1 .
$$
5. Define $\lambda := -\Phi_2/\Psi_2$. It depends only on $\eta_2$, which was fixed in step 1, and not on $\eta_1$. Step 4 reads $\Phi_1 + \lambda\Psi_1 = 0$, i.e.
$$
\int_a^b \big(E_J + \lambda E_K\big)\,\eta_1\,dx = 0 \quad\text{for every admissible } \eta_1 .
$$
6. By the fundamental lemma, $E_J + \lambda E_K = 0$ on $[a, b]$. ∎

■ **What just happened.** Exactly the finite-dimensional multiplier rule, $\nabla f = -\lambda\nabla g$, run in the two-dimensional slice of function space spanned by $\eta_1$ and $\eta_2$. The hypothesis "$y$ is not an extremal of $K$" is the function-space version of "$\nabla g \ne 0$" (the constraint qualification in ML-2's KKT conditions). In practice: form $L^* = L + \lambda G$, write its E–L equation, solve with $\lambda$ as an unknown constant, and fix $\lambda$ from the constraint at the end. Example 4.4 does this for the hanging chain.

**Dido's problem.** Among curves $y(x) \ge 0$ on $[-w, w]$ with $y(\pm w) = 0$ and fixed length $\ell$, $2w < \ell \le \pi w$ (beyond $\pi w$ the answer is no longer a graph over $[-w, w]$), which encloses the most area with the $x$-axis? Maximize $\int y\,dx$ subject to $\int\sqrt{1 + y'^2}\,dx = \ell$. With $L^* = y + \lambda\sqrt{1 + y'^2}$: $L^*_y = 1$ and $L^*_{y'} = \lambda y'/\sqrt{1 + y'^2}$, so E–L says $\frac{d}{dx}\big[\lambda y'/\sqrt{1 + y'^2}\big] = 1$, hence $\lambda y'/\sqrt{1 + y'^2} = x - x_0$. Solve for $y'$: $\lambda^2 y'^2 = (x - x_0)^2(1 + y'^2)$, so $|y'| = |x - x_0|/\sqrt{\lambda^2 - (x - x_0)^2}$. The sign: $y'/\sqrt{1 + y'^2} = (x - x_0)/\lambda$, and for an arc that bulges *upward* (rising to the left of $x_0$, falling to the right) we need $\lambda < 0$, which gives $y' = -(x - x_0)/\sqrt{\lambda^2 - (x - x_0)^2}$. Integrating, $y - y_0 = +\sqrt{\lambda^2 - (x - x_0)^2}$, i.e. $(x - x_0)^2 + (y - y_0)^2 = \lambda^2$. **A circular arc of radius $|\lambda|$.** Symmetry gives $x_0 = 0$; the two remaining constants $y_0, \lambda$ come from $y(w) = 0$ and the length. When the "shoreline" is straight and the perimeter is free to be placed, the arc is a semicircle, the answer Dido is said to have found.

**Holonomic constraints.** A different kind of constraint is pointwise: $g(x, y_1, y_2) = 0$ for *every* $x$ (a bead on a wire, a pendulum bob on a rod). The rule is: form $L^* = L + \lambda(x)\,g$ with a multiplier that is a **function** of $x$, write the E–L equations (3.13) for $y_1$ and $y_2$ treating $\lambda(x)$ as given, and solve together with $g = 0$. Why a function? Because a pointwise constraint is infinitely many constraints, one per $x$, and each gets its own multiplier. We state this without derivation; CM-3 derives it and shows that $\lambda(x)\,\partial g/\partial y_k$ is the constraint force, which is why engineers care about the multiplier and not only the motion.

### 3.11 The second variation and the Legendre condition

Return to $\Phi(\varepsilon) = J[y + \varepsilon\eta]$. The first variation is $\Phi'(0)$; the **second variation** is
$$
\delta^2 J := \Phi''(0) = \int_a^b \Big[ L_{yy}\,\eta^2 + 2L_{yy'}\,\eta\eta' + L_{y'y'}\,\eta'^2 \Big]\,dx, \tag{3.17}
$$
obtained by differentiating under the integral once more (chain rule on each of $L_y\eta$ and $L_{y'}\eta'$). If $y$ is a local minimum, $\Phi$ has a minimum at $0$ for every $\eta$, so $\delta^2 J \ge 0$ for every admissible $\eta$: the analogue of "Hessian positive semidefinite." Which term dominates?

▶ **Try it first.** Take $\eta$ small in size but wildly oscillating, so that $\eta'$ is large. Which of the three terms in (3.17) wins?

**Derivation (Legendre's necessary condition).**

1. Fix $x_0 \in (a, b)$ and suppose $L_{y'y'}(x_0, y(x_0), y'(x_0)) < 0$. By continuity there is an interval $[c, d] \ni x_0$ inside $(a, b)$ on which $L_{y'y'} < -m$ for some $m > 0$, while $|L_{yy}| \le M$ and $|L_{yy'}| \le M$ there.
2. Choose $\eta(x) = \sin^2\!\big(k(x - c)\big)$ on $[c, d]$ and $0$ elsewhere, with $k = n\pi/(d - c)$ for an integer $n \ge 1$ so that $\eta$ and $\eta' = k\sin(2k(x - c))$ both vanish at $c$ and $d$. Then $\eta$ is admissible, $0 \le \eta \le 1$, and $|\eta'| \le k$.
3. Over a whole number of periods, $\sin^4$ averages to $3/8$ and $\sin^2$ averages to $1/2$, so $\int_c^d \eta^2\,dx = \tfrac38(d - c)$ and $\int_c^d \eta'^2\,dx = \tfrac12 k^2 (d - c)$, while $\big|\int_c^d 2L_{yy'}\eta\eta'\,dx\big| \le 2Mk(d - c)$.
4. Therefore
$$
\delta^2 J \le M\cdot\tfrac38(d - c) + 2Mk(d - c) - m\cdot\tfrac12 k^2(d - c) .
$$
The right side is a quadratic in $k$ with negative leading coefficient, so for $n$ (hence $k$) large enough, $\delta^2 J < 0$, and $y$ is not a minimum. ∎

**Theorem 3.8 (Legendre condition).** If $y$ is a weak local minimum of $J$, then $L_{y'y'}(x, y(x), y'(x)) \ge 0$ for every $x \in [a, b]$. (For a maximum, $\le 0$.)

■ **What just happened.** The $\eta'^2$ term wins because differentiating a rapid oscillation multiplies it by $k$. So the sign of $L_{y'y'}$ along the extremal is the first thing to check. For arc length, $L_{y'y'} = (1 + y'^2)^{-3/2} > 0$: consistent with a minimum. For the brachistochrone, $L_{y'y'} = (1 + y'^2)^{-3/2}/\sqrt{2gy} > 0$. For the action $\tfrac12 m\dot q^2 - V$, $L_{\dot q\dot q} = m > 0$, which is why Hamilton's principle is usually (not always) "least" action.

The Legendre condition is necessary, not sufficient. The sufficient condition for a weak local minimum is strict Legendre ($L_{y'y'} > 0$) plus **Jacobi's condition**: the interval $(a, b]$ contains no point conjugate to $a$, where a conjugate point is where a neighboring extremal through $(a, A)$ crosses $y$ again. Geodesics on a sphere show what goes wrong: all great circles through a point $P$ meet again at the antipode $P^*$, so an arc from $P$ that goes past $P^*$ contains a conjugate point and is not a minimum, even though it satisfies E–L and Legendre. We state Jacobi's condition without proof; Gelfand & Fomin ch. 5 prove it, and you do not need it for CM-2.

### 3.12 The functional derivative

Equation (3.5) says $\delta J = \int_a^b \big[L_y - \frac{d}{dx}L_{y'}\big]\eta\,dx$. Compare with $\mathbb{R}^n$, where $\frac{d}{d\varepsilon}f(\mathbf{x} + \varepsilon\mathbf{h}) = \sum_i \frac{\partial f}{\partial x_i}h_i$. The sum over $i$ has become an integral over $x$, the components $h_i$ have become the values $\eta(x)$, and the role of $\partial f/\partial x_i$ is played by the bracket.

**Definition 3.9.** The **functional derivative** of $J$ at $y$ is the function $\dfrac{\delta J}{\delta y}(x)$ such that
$$
\delta J = \int_a^b \frac{\delta J}{\delta y}(x)\,\eta(x)\,dx \quad\text{for every admissible } \eta. \tag{3.18}
$$
For $J = \int L\,dx$ with fixed endpoints,
$$
\frac{\delta J}{\delta y}(x) = \frac{\partial L}{\partial y} - \frac{d}{dx}\frac{\partial L}{\partial y'}, \tag{3.19}
$$
and the Euler–Lagrange equation is $\delta J/\delta y = 0$: **the gradient of the landscape is zero.** The discrete statement (3.9) is $\partial J_h/\partial y_i \approx h\,\frac{\delta J}{\delta y}(x_i)$; the functional derivative is the gradient per unit length of $x$. The same symbol appears in field theory as $\delta S/\delta\phi(x)$, and in ML-10, where the functional is the ELBO or a KL divergence and the function is a density $q(z)$: setting $\delta F/\delta q = 0$ with a multiplier for $\int q = 1$ is how the optimal variational distribution is derived. There $L$ contains $q$ but not $q'$, so by Case B of §3.6 the equation is algebraic and gives $q \propto \exp(\ldots)$. SM-4's maximum-entropy derivation of the Boltzmann distribution is the same computation with different constraints.

### 3.13 Hamilton's principle, one paragraph

Let a particle of mass $m$ move on a line with position $q(t)$ in a potential $V(q)$. Define the **action** $S[q] = \int_{t_1}^{t_2}\big(\tfrac12 m\dot q^2 - V(q)\big)\,dt$. Its Lagrangian $L = T - V$ has $L_q = -V'(q)$ and $L_{\dot q} = m\dot q$, so the Euler–Lagrange equation (3.4) with $x \to t$ is $-V'(q) - \frac{d}{dt}(m\dot q) = 0$, i.e.
$$
m\ddot q = -\frac{dV}{dq} = F(q) .
$$
Newton's second law is the statement that the actual motion makes the action stationary. Since $L$ has no explicit $t$, the Beltrami quantity $L - \dot q L_{\dot q} = \tfrac12 m\dot q^2 - V - m\dot q^2 = -(T + V)$ is constant: energy conservation falls out of §3.6 (Example 4.5 checks this by hand). That is all we say here. CM-2 asks the questions this paragraph skips: why $T - V$ and not $T + V$, what happens with several particles and constraints, and why "stationary" rather than "least."

### 3.14 How to solve a variational problem

1. Write the quantity to be extremized as $\int_a^b L(x, y, y')\,dx$ and identify $L$. Decide what is fixed at the endpoints.
2. If there is a constraint $\int G\,dx = \ell$, replace $L$ by $L + \lambda G$ with $\lambda$ an unknown constant.
3. Look at the formula for $L$. No $y$? Use (3.10). No $y'$? The condition is algebraic. No explicit $x$? Use Beltrami (3.11). Otherwise write (3.4) in full.
4. Solve the resulting first- or second-order ODE. Expect a substitution.
5. Impose the boundary conditions; if an endpoint is free, impose the natural condition (3.15) there. Impose the constraint to find $\lambda$.
6. Check the Legendre sign $L_{y'y'}$ to know whether you could be looking at a minimum, a maximum, or neither. If the physics says there must be a minimum and your extremal is unique, you have it.

## 4. Worked examples

### Example 4.1: The shortest path in the plane is a straight line

Minimize $J[y] = \int_a^b\sqrt{1 + y'^2}\,dx$ with $y(a) = A$, $y(b) = B$.

$L = \sqrt{1 + y'^2}$ contains neither $y$ nor $x$. Use Case A: $L_{y'} = y'/\sqrt{1 + y'^2} = C$. Square and solve: $y'^2 = C^2(1 + y'^2)$, so $y'^2 = C^2/(1 - C^2)$, so $y' = m$ for a constant $m$ (we need $|C| < 1$, which is automatic since $|y'| < \sqrt{1 + y'^2}$). Integrate: $y = mx + c$, and the two endpoint conditions fix $m = (B - A)/(b - a)$ and $c$. The full E–L equation (3.6) says the same thing with more work: $L_{y'y'} = (1 + y'^2)^{-3/2}$, $L_{y'y} = L_{y'x} = L_y = 0$, so (3.6) is $y''/(1 + y'^2)^{3/2} = 0$, i.e. $y'' = 0$. Legendre: $L_{y'y'} > 0$. Exercise C1 proves this extremal is the *global* minimum. The project's minimizer, started from a wiggly curve, returns to the line to within $10^{-6}$.

### Example 4.2: The brachistochrone

**Setup.** Measure $y$ *downward* from the start, so the bead begins at $(0, 0)$ and ends at $(x_1, y_1)$ with $x_1, y_1 > 0$. It starts from rest; friction is absent; so energy conservation gives $\tfrac12 v^2 = gy$, i.e. $v = \sqrt{2gy}$. The time to traverse an arc element $ds = \sqrt{1 + y'^2}\,dx$ is $ds/v$, so
$$
T[y] = \int_0^{x_1}\frac{\sqrt{1 + y'^2}}{\sqrt{2gy}}\,dx, \qquad L = \frac{\sqrt{1 + y'^2}}{\sqrt{2gy}} .
$$
**Beltrami.** $L$ has no explicit $x$. Compute $L_{y'} = \dfrac{y'}{\sqrt{1 + y'^2}\sqrt{2gy}}$, so
$$
L - y'L_{y'} = \frac{1 + y'^2 - y'^2}{\sqrt{1 + y'^2}\sqrt{2gy}} = \frac{1}{\sqrt{2gy\,(1 + y'^2)}} = \text{const}.
$$
Hence $y(1 + y'^2) = 2a$ for some constant $a > 0$ (the factor 2 is hindsight). This is the first integral; the project checks it is constant along the solution and not along the straight line.

**Solve the ODE.** $y'^2 = (2a - y)/y$, so on the part of the curve where $y$ increases with $x$,
$$
\frac{dy}{dx} = \sqrt{\frac{2a - y}{y}} \quad\Longrightarrow\quad dx = \sqrt{\frac{y}{2a - y}}\,dy .
$$
The integral $\int\sqrt{y/(2a - y)}\,dy$ needs a substitution, and the one that works is
$$
y = a(1 - \cos\theta) = 2a\sin^2(\theta/2), \qquad dy = a\sin\theta\,d\theta .
$$
Then $2a - y = a(1 + \cos\theta) = 2a\cos^2(\theta/2)$, so $\sqrt{y/(2a - y)} = \tan(\theta/2)$, and
$$
dx = \tan(\theta/2)\cdot a\sin\theta\,d\theta = \tan(\theta/2)\cdot 2a\sin(\theta/2)\cos(\theta/2)\,d\theta = 2a\sin^2(\theta/2)\,d\theta = a(1 - \cos\theta)\,d\theta .
$$
Integrate: $x = a(\theta - \sin\theta) + C$. At the start $y = 0$, so $\theta = 0$, and $x = 0$ forces $C = 0$. The solution is the **cycloid**
$$
x = a(\theta - \sin\theta), \qquad y = a(1 - \cos\theta), \tag{4.1}
$$
the path traced by a point on the rim of a wheel of radius $a$ rolling along the underside of the $x$-axis. (Where $y$ decreases with $x$, past the bottom of the arch, the other sign of the square root applies; the same parametrization covers it with $\theta > \pi$.)

**The start is vertical.** $\dfrac{dy}{dx} = \dfrac{\sin\theta}{1 - \cos\theta} = \cot(\theta/2) \to \infty$ as $\theta \to 0$. The fastest wire leaves the start pointing straight down, however far to the side the endpoint is. This is physically sensible (get speed first, then spend it on horizontal distance) and it is what makes the brachistochrone awkward numerically (§7).

**Fit the endpoint.** Two unknowns, $a$ and the final angle $\theta_1$, and two equations: $a(\theta_1 - \sin\theta_1) = x_1$, $a(1 - \cos\theta_1) = y_1$. Divide to eliminate $a$:
$$
\frac{\theta_1 - \sin\theta_1}{1 - \cos\theta_1} = \frac{x_1}{y_1}, \tag{4.2}
$$
a transcendental equation solved numerically (the left side increases from $0$ at $\theta_1 = 0$ to $\infty$ at $\theta_1 = 2\pi$, so there is exactly one root). Then $a = y_1/(1 - \cos\theta_1)$.

**The descent time.** Along the cycloid, $ds = \sqrt{(dx/d\theta)^2 + (dy/d\theta)^2}\,d\theta = a\sqrt{(1 - \cos\theta)^2 + \sin^2\theta}\,d\theta = a\sqrt{2 - 2\cos\theta}\,d\theta = 2a\sin(\theta/2)\,d\theta$, and $v = \sqrt{2gy} = \sqrt{2ga(1 - \cos\theta)} = 2\sqrt{ga}\,\sin(\theta/2)$. So $dt = ds/v = \sqrt{a/g}\,d\theta$: the bead advances at a *constant rate in $\theta$*, and
$$
T_{\text{cyc}} = \theta_1\sqrt{a/g} . \tag{4.3}
$$

**Numbers.** Descend from $(0, 0)$ to a point $1\ \mathrm{m}$ to the right and $1\ \mathrm{m}$ lower. (With $y$ measured upward this endpoint is $(1, -1)$; in our downward convention it is $(x_1, y_1) = (1, 1)$.) Equation (4.2) with $x_1/y_1 = 1$ gives $\theta_1 = 2.4120\ \mathrm{rad}$ ($138.2°$), then $a = 1/(1 - \cos 2.4120) = 0.57292\ \mathrm{m}$, and
$$
T_{\text{cyc}} = 2.4120\sqrt{0.57292/9.81} = 0.5829\ \mathrm{s}.
$$
Compare the straight wire $y = x$: it is an incline of length $\sqrt2$ at $45°$, acceleration $g\sin 45° = g/\sqrt2$, so $\sqrt2 = \tfrac12(g/\sqrt2)T^2$ gives $T_{\text{line}} = 2/\sqrt{g} = 0.6386\ \mathrm{s}$. (Or integrate $T[y]$ with $y' = 1$: $\int_0^1\sqrt2/\sqrt{2gx}\,dx = 2/\sqrt g$, same thing.) The cycloid wins by $8.7\%$. Not dramatic for this geometry, but it is the best possible, and it beats the circular arc Galileo proposed too: the quarter circle through the same two points that starts vertically, $(x, y) = (1 - \cos\varphi, \sin\varphi)$, takes $\int_0^{\pi/2} d\varphi/\sqrt{2g\sin\varphi} = 0.5920\ \mathrm{s}$, between the cycloid and the line. The cycloid arrives at the endpoint with slope $\cot(\theta_1/2) = 0.382$, i.e. $20.9°$ below horizontal. (All of these numbers are asserted in Part 2 of the project.)

**Dipping below the endpoint.** The cycloid's lowest point is at $\theta = \pi$, depth $2a$. If $\theta_1 > \pi$ the wire goes *below* the endpoint and climbs back up to it. From (4.2), $\theta_1 = \pi$ corresponds to $x_1/y_1 = \pi/2$, so the fastest path dips below its destination whenever the endpoint is more than $\pi/2 \approx 1.57$ times as far sideways as it is down. For $(x_1, y_1) = (3, 1)$: $\theta_1 = 4.052 > \pi$, $a = 0.6197$, and the wire descends to depth $2a = 1.239$ before rising to $1$. Students find this hard to believe. The project draws it.

### Example 4.3: Fermat's principle and Snell's law

Light travels at speed $c/n$ in a medium of refractive index $n$. Fermat's principle: between two points, light follows the path of least time. If $n = n(y)$ depends only on height, and the path is a graph $y(x)$, the travel time is
$$
T[y] = \frac1c\int_a^b n(y)\sqrt{1 + y'^2}\,dx, \qquad L = n(y)\sqrt{1 + y'^2}
$$
(the factor $1/c$ is irrelevant to the minimization). $L$ has no explicit $x$, so Beltrami:
$$
L - y'L_{y'} = n(y)\sqrt{1 + y'^2} - \frac{n(y)\,y'^2}{\sqrt{1 + y'^2}} = \frac{n(y)}{\sqrt{1 + y'^2}} = \text{const}. \tag{4.4}
$$
Geometry: let $\theta$ be the angle between the ray and the **vertical** (the normal to the layers of constant $n$). The ray direction is $(1, y')/\sqrt{1 + y'^2}$, and its horizontal component is $\sin\theta = 1/\sqrt{1 + y'^2}$. So (4.4) reads
$$
n(y)\sin\theta(y) = \text{const} \tag{4.5}
$$
along the whole ray. For a continuously varying $n$ this is the law of the mirage and of gradient-index fibers (Exercise B1). For **piecewise-constant** $n$, say $n_1$ above the line $y = 0$ and $n_2$ below, the ray is straight in each half (there $n$ is constant, so Example 4.1 applies), and (4.5) evaluated just above and just below the interface gives
$$
n_1\sin\theta_1 = n_2\sin\theta_2,
$$
**Snell's law**, with $\theta_1, \theta_2$ the angles of incidence and refraction measured from the normal. (At the interface $y'$ jumps, so strictly our $C^2$ hypothesis fails there; the honest version minimizes over the crossing point $x_c$ directly, as the project does. With the start at $(0, z_1)$, a height $z_1$ above the interface, and the end at $(D, -z_2)$, a depth $z_2$ below it, $c\,T(x_c) = n_1\sqrt{x_c^2 + z_1^2} + n_2\sqrt{(D - x_c)^2 + z_2^2}$, and $T'(x_c) = 0$ reads $n_1 x_c/\sqrt{x_c^2 + z_1^2} = n_2(D - x_c)/\sqrt{(D - x_c)^2 + z_2^2}$, which is Snell's law in one line.) **Numbers:** from $A = (0, 1)$ in air ($n_1 = 1$) to $B = (2, -1)$ in glass ($n_2 = 1.5$), the time-minimizing crossing point is $x_c = 1.3627$, giving $\theta_1 = 53.73°$, $\theta_2 = 32.51°$, and $n_1\sin\theta_1 = n_2\sin\theta_2 = 0.8062$.

### Example 4.4: The catenary, with the multiplier

A uniform chain of linear density $\rho$ and length $\ell$ hangs between $(x_0', y_0')$ and $(x_1', y_1')$ with $\ell$ greater than the distance between them. Its potential energy is $\int\rho g\,y\,ds$, so it minimizes
$$
V[y] = \rho g\int y\sqrt{1 + y'^2}\,dx \quad\text{subject to}\quad K[y] = \int\sqrt{1 + y'^2}\,dx = \ell .
$$
**Multiplier.** By Theorem 3.7, drop the constant $\rho g$ and use $L^* = (y + \lambda)\sqrt{1 + y'^2}$. It has no explicit $x$, so Beltrami:
$$
L^* - y'L^*_{y'} = (y + \lambda)\sqrt{1 + y'^2} - \frac{(y + \lambda)y'^2}{\sqrt{1 + y'^2}} = \frac{y + \lambda}{\sqrt{1 + y'^2}} = c . \tag{4.6}
$$
**Solve.** $(y + \lambda)^2 = c^2(1 + y'^2)$, so $y' = \sqrt{\big((y + \lambda)/c\big)^2 - 1}$. Substitute $y + \lambda = c\cosh u$: then $y' = c\sinh u\cdot u'$ and $\sqrt{\cosh^2 u - 1} = \sinh u$, so $c\sinh u\,u' = \sinh u$, i.e. $u' = 1/c$, $u = (x - x_0)/c$. Hence
$$
y(x) = -\lambda + c\cosh\!\left(\frac{x - x_0}{c}\right). \tag{4.7}
$$
Three constants: $c$, $x_0$, $\lambda$. Three conditions: the two endpoints and the length. If you measure $y$ from the level $y = -\lambda$ (the catenary's *directrix*), the shape is $y = c\cosh((x - x_0)/c)$, which is what tables quote. The multiplier has a meaning: the tension in the chain is $\rho g(y + \lambda)$, and its horizontal component is the constant $\rho g c$, so $\lambda$ is the depth of the directrix below the level $y = 0$ and $c$ is the horizontal tension in units of weight per length.

**Numbers.** Chain of length $\ell = 3\ \mathrm{m}$ between $(-1, 0)$ and $(1, 0)$. Symmetry gives $x_0 = 0$. The endpoint condition $y(\pm1) = 0$ gives $\lambda = c\cosh(1/c)$. The length of (4.7) is $\int_{-1}^{1}\sqrt{1 + \sinh^2(x/c)}\,dx = \int_{-1}^1\cosh(x/c)\,dx = 2c\sinh(1/c)$, so $c$ solves
$$
2c\sinh(1/c) = 3 \quad\Longrightarrow\quad c = 0.61647\ \mathrm{m},
$$
by bisection. Then $\lambda = 0.61647\cosh(1.62213) = 1.6217\ \mathrm{m}$ and the lowest point is $y(0) = c - \lambda = -1.0053\ \mathrm{m}$: a three-metre chain between points two metres apart sags one metre. The project minimizes the discretized energy with a length constraint and recovers this curve to $2.3\times10^{-4}$ (40 interior nodes); a second constrained optimizer also reports the multiplier it used, $1.6216$, matching $\lambda$ to $1.5\times10^{-4}$.

### Example 4.5: Beltrami by hand for a falling particle

Take $L = \tfrac12 y'^2 - y$ (unit mass, $g = 1$, $y$ up, independent variable $x$ playing the role of time). E–L: $L_y = -1$, $L_{y'} = y'$, so $-1 - y'' = 0$, i.e. $y'' = -1$: free fall. Beltrami quantity:
$$
H := L - y'L_{y'} = \tfrac12 y'^2 - y - y'^2 = -\big(\tfrac12 y'^2 + y\big).
$$
Check (3.12) directly: $\frac{dH}{dx} = -(y'y'' + y') = y'\,(-1 - y'')$, which is $y'$ times the E–L expression $L_y - \frac{d}{dx}L_{y'} = -1 - y''$, exactly as (3.12) says, and it vanishes exactly when $y'' = -1$. So $\tfrac12 y'^2 + y$ is constant along solutions: kinetic plus potential energy. On the solution $y = -\tfrac12 x^2 + v_0 x + y_0$: $\tfrac12 y'^2 + y = \tfrac12(v_0 - x)^2 - \tfrac12 x^2 + v_0 x + y_0 = \tfrac12 v_0^2 + y_0$. Constant, as promised. On a non-solution such as $y = x^3$, $H = -(\tfrac92 x^4 + x^3)$ is not constant. The project verifies this symbolically.

### Example 4.6: The catenoid, briefly

Rotate the graph of $y > 0$ about the $x$-axis; the area is $2\pi\int y\sqrt{1 + y'^2}\,dx$. This is the catenary Lagrangian with $\lambda = 0$, so the minimal surface of revolution is $y = c\cosh((x - x_0)/c)$ rotated: the **catenoid**, the shape of a soap film between two coaxial rings. Two rings of radius $R$ at $x = \pm d$ require $c\cosh(d/c) = R$, which has a solution only when $d/R \le 0.6627$; pull the rings further apart and no catenoid exists, the film breaks, and the area-minimizing configuration is two flat disks (the Goldschmidt solution). So E–L can have *no* admissible solution, and the true minimum then lives outside the class of smooth graphs. The project computes the critical ratio.

## 5. Why it matters and where it connects

- **CM-2** takes §3.13 and asks why $L = T - V$. Hamilton's principle, Fermat's principle (Example 4.3) and the brachistochrone are the same mathematics with different $L$; CM-2 is pairing Taylor ch. 6 and Feynman II-19 for exactly this reason.
- **CM-3** uses §3.7 (one E–L equation per generalized coordinate) and the holonomic-constraint rule of §3.10 with $\lambda(t)$; the multiplier there is the constraint force.
- **CM-4** proves Noether's theorem; the Beltrami identity (3.11) is its time-translation case (energy), and Case A of §3.6 is its space-translation case (momentum).
- **CM-8** calls $y'L_{y'} - L$ the Hamiltonian and $L_{y'}$ the canonical momentum; the Legendre transform that takes you from $L$ to $H$ is the same one that takes you from internal energy to free energy in SM-2.
- **EM-2, EM-3**: the Dirichlet principle (Exercise C3: Laplace's equation is the E–L equation of $\int|\nabla u|^2$; EM-2 §3.12 points back to it) turns the electrostatic boundary-value problems that EM-3 solves into minimization problems. The beam equation of §3.8 does the same for structures in ME, and minimizing a discretized energy is what finite-element software actually does.
- **CM-6**: the loaded string's potential energy $V = \frac{\tau}{2d}\sum_k (y_{k+1} - y_k)^2$ is a discretized $\int\tfrac12\tau\,y'^2\,dx$ of exactly the form (3.7), and its stiffness matrix is the Hessian of that $J_h$.
- **QM-8**: the variational principle for the ground-state energy, $E_0 \le \langle\psi|H|\psi\rangle/\langle\psi|\psi\rangle$, is a functional of $\psi$ whose stationary points are the eigenstates; its E–L equation is the Schrödinger equation.
- **QM-9**: Feynman's path integral weights every path by $e^{iS/\hbar}$; as $\hbar \to 0$ the paths near the stationary action reinforce and all others cancel, which is why the classical path is the one where $\delta S = 0$. This module's $\delta$ is that $\delta$.
- **LA-14** defines the gradient as the vector $\nabla f$ with $df = \nabla f\cdot d\mathbf{x}$; §3.12's functional derivative is defined by the same sentence with an integral in place of the dot product, and (3.9) is the exact relation between them.
- **ML-2**: the discretized problem (3.7) is a smooth unconstrained optimization, which is why L-BFGS solves it in the project; Theorem 3.7 is the Lagrange multiplier rule whose finite-dimensional version becomes KKT there.
- **ML-10, SM-4, F5**: variational inference maximizes the ELBO over densities $q$; maximum entropy maximizes $-\int q\ln q$ subject to moment constraints. Both are §3.6 Case B with multipliers, and Exercise C4 does the Gaussian case. The word "variational" in "variational autoencoder" is this module's word.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are drills: computing E–L equations, a first variation by hand, recognizing which first integral applies, and a natural boundary condition. B1–B8 are the module: B1 (light in a graded medium), B3 (the du Bois-Reymond lemma), B4 (the discrete derivation with a different quadrature), B6 (an isoperimetric problem with a computable multiplier), B7 (the brachistochrone to a vertical wall, using the natural boundary condition), and B8 (extending `project.py` to the catenoid) are essential. C1–C4 go beyond: a global-minimality proof, the brachistochrone that returns to its starting height, the Euler–Lagrange equation for a field, and the maximum-entropy derivation of the Gaussian.

## 7. Coding project

`project.py` is a discretized functional minimizer plus the checks that connect it to the theory. Its parts are numbered as below; Exercises B4 and B8 refer to these part numbers.

1. **Part 1, arc length → straight line.** A function $y$ is stored as its values at $N$ interior nodes and the functional is evaluated on the polygon through them (for arc length the trapezoid rule is exact: $J_h = \sum_i\sqrt{h^2 + \Delta y_i^2}$). L-BFGS-B, started from a wiggly curve, returns the straight line to a maximum deviation asserted $< 10^{-3}$ (it comes out near $2\times10^{-7}$). Along the iterations the Euler–Lagrange expression $-y''/(1 + y'^2)^{3/2}$, evaluated by central differences at the nodes, is tracked and plotted collapsing from about $800$ to $4\times10^{-6}$.
2. **Part 2, the brachistochrone.** The polygon functional is the *exact* sliding time along each straight segment, $2\Delta s/(v_i + v_{i+1})$ with $v = \sqrt{2gy}$. The minimizer recovers the cycloid with a descent time within $0.5\%$ of (4.3) (it comes out $0.06\%$ high, asserted positive and below $0.1\%$, as a faithful discretization must: it cannot beat the true minimum). Then **a discretization that lies**: the same problem with the midpoint rule for $\int ds/\sqrt{2gy}$ finds a polygon with a *smaller* discrete time than the true minimum, impossible for a faithful discretization, by diving almost vertically at the first node, where the integrand is singular. Watch for this whenever an integrand blows up. The part ends by drawing the $(3, 1)$ case, where the wire dips below its destination.
3. **Part 3, shooting and the Beltrami test.** The brachistochrone's E–L equation, reduced by Beltrami to $dx/dy = \sqrt{y/(2a - y)}$, is integrated with `solve_ivp` from the start, and `brentq` adjusts the constant $2a$ until the path hits the endpoint; the constant matches the cycloid's $2a$ to $10^{-9}$. (The usual shooting parameter, the initial slope, is infinite here, so the first-integral constant takes its place.) Then $H = L - y'L_{y'}$ is evaluated along the cycloid (constant to $10^{-15}$), along the numerical polygon from Part 2 (relative spread $5\times10^{-4}$ away from the vertical start, asserted below $6\times10^{-4}$), and along the straight line, where it varies by a factor of $9$ (asserted above $8.5$).
4. **Part 4, the catenary.** The chain's energy $\sum_i \Delta s_i\,(y_i + y_{i+1})/2$ is minimized with the length constraint $\sum_i\Delta s_i = \ell$ through SLSQP's equality-constraint interface and recovers (4.7) to $2.3\times10^{-4}$ (asserted below $3\times10^{-4}$); a second solver (`trust-constr`) reports the multiplier it used, $1.6216$, which matches $\lambda$ from Example 4.4 to $1.5\times10^{-4}$ (asserted below $2\times10^{-4}$).
5. **Part 5, Snell's law** from Fermat's principle with two media: minimize the optical path over the crossing point, compare $n_1\sin\theta_1$ with $n_2\sin\theta_2$.
6. **Part 6, geodesics on a sphere** with two dependent variables $(\theta(t), \varphi(t))$, by minimizing the energy $\int(\dot\theta^2 + \sin^2\theta\,\dot\varphi^2)\,dt$ (whose minimizer is the great circle at constant speed; minimizing length directly is degenerate because any reparametrization has the same length). Compared pointwise with the great circle.
7. **Part 7, $J[y + \varepsilon\eta]$ versus $\varepsilon$** for the straight line (flat at $\varepsilon = 0$) and for $y = x^2$ (sloped), the picture of §2; the finite-difference slopes at $\varepsilon = 0$ are checked against (3.3).
8. **Part 8, gradient versus E–L.** The gradient (3.8) of the discretized functional (3.7) is compared with $h$ times the continuum E–L expression for a generic Lagrangian on a smooth non-solution curve; the error drops by a factor of $4$ each time $h$ is halved, which is (3.9). **Part 8b** repeats this for the trapezoid-rule discretization of Exercise B4, where the factor is $8$.
9. **Part 9, SymPy.** Identity (3.12) on four Lagrangians, Example 4.5's energy, and the fact that the line, the cycloid and $\cosh$ satisfy their E–L equations are checked symbolically, as is the sign of $\lambda$ for Dido's upward arc in §3.10; the catenoid's critical ratio $0.6627$ of Example 4.6 is computed, and the headline numbers of Exercises B6, B7 and C2 are asserted.
10. **Part 10, exercise references.** Every other number quoted in `solutions.md` is computed or asserted here, and the reference implementation for Exercise B8 (the catenoid between two rings, including the case where no catenoid exists) is at its end. Do the exercises before reading it.

**Build it yourself** before reading the code: `polygon_length(y_int, x, ya, yb)` (the trapezoid-rule arc length of the polygon through the nodes), `polygon_descent_time(y_int, x, ya, yb)` (sum of $2\Delta s/(v_i + v_{i+1})$ with $v = \sqrt{2gy}$), and `el_residual_arclength(y, h)` (the E–L expression $-y''/(1 + y'^2)^{3/2}$ by central differences). Then minimize the first two with L-BFGS-B and compare.

**Your turn** (at the end of the script): the catenoid between two rings, including the case where no catenoid exists (this is Exercise B8; a reference implementation sits at the end of Part 10); the brachistochrone with Coulomb friction $\mu$ (replace the speed law by $v^2 = 2g(y - \mu x)$, which follows from taking the normal force to be the weight's component $mg\cos\alpha$ and neglecting its centripetal part $mv^2\kappa$, the usual simplification: the friction work over $ds$ is then $\mu mg\cos\alpha\,ds = \mu mg\,dx$); and a one-dimensional particle with $q(0) = q(T) = 0$ whose discretized action $\sum_i h\big[\tfrac12(\Delta q_i/h)^2 - g\,(q_i + q_{i+1})/2\big]$ is minimized numerically, with the discrete trajectory compared against the parabola $q(t) = \tfrac g2\,t(T - t)$ and its second difference against $-g$, to recover Newton's law from the optimizer.

## 8. Where you'll struggle

- **$\delta$ versus $d$.** $\delta y$ changes the *function* at fixed $x$ (move the curve up or down); $dy$ moves *along* the curve to a neighboring $x$. Physics books write $\delta L = L_y\,\delta y + L_{y'}\,\delta y'$ and then "integrate by parts," which hides that $\delta y'$ means $(\delta y)'$ and that $\delta$ commutes with $d/dx$ for a reason (it is $\varepsilon\eta'$). Every time you are confused, go back to the definition: $\delta J$ is $\frac{d}{d\varepsilon}J[y + \varepsilon\eta]$ at $\varepsilon = 0$, an ordinary derivative of an ordinary function of one real variable. That definition is rigorous, it is what we used in every derivation, and nothing in the $\delta$ notation adds to it except brevity.
- **Partial versus total derivatives of $L$.** In $\partial L/\partial y'$, treat $y'$ as a letter, not as the derivative of $y$. Students compute $\partial L/\partial y$ for $L = \sqrt{1 + y'^2}$ and write something nonzero because "$y'$ depends on $y$." It does not, in the slot sense. Only the $\frac{d}{dx}$ in (3.4) sees the curve, and that is where the chain rule (3.6) lives.
- **Stationary points are not minima**, and you will call them minima anyway. E–L is $f' = 0$. The Legendre sign is the cheapest sanity check; the sphere's long great-circle arc is the counterexample to keep in mind. In mechanics the "least" in "least action" is a historical misnomer that CM-2 discusses.
- **The brachistochrone substitution.** Nobody finds $y = a(1 - \cos\theta)$ unaided; it is found by knowing the answer. What you *can* do alone is recognize $\int\sqrt{y/(2a - y)}\,dy$ as calling for a substitution that rationalizes the square root of a ratio, and $y = 2a\sin^2\phi$ is the standard one for $\sqrt{y/(2a - y)}$. Learn it as a move, not a miracle.
- **"$L$ has no explicit $x$" is read off the formula.** For $L = \sqrt{1 + y'^2}/\sqrt{2gy}$ there is no $x$ in the formula, so Beltrami applies; the fact that $y$ is a function of $x$ is irrelevant to this test. Conversely $L = x\,y'^2$ *does* contain $x$ and Beltrami is false for it, even though $L_y = 0$ (use Case A instead).
- **The fundamental lemma gets skipped**, and then "$\int(\ldots)\eta = 0$ for all $\eta$, so $(\ldots) = 0$" feels like a leap. It is not a leap; it is Theorem 3.4, and the bump function is the whole proof. Reproduce it once closed-book.
- **Three constants in the catenary** (and in most isoperimetric problems): $c$, $x_0$, and the multiplier $\lambda$. Students fit two and wonder why the length is wrong. Count your constants against your conditions before solving anything.
- **The cycloid dips below the endpoint** when $x_1/y_1 > \pi/2$, and it leaves vertically no matter what. Both offend intuition. Both are correct. The project draws both.
- **Numerics can cheat.** The midpoint-rule brachistochrone in the project finds a "better than optimal" time by exploiting the quadrature near the singular start. The discretized functional must be a faithful evaluation of the real functional on the discrete curve; when the integrand is singular, use the exact expression for each segment.

## 9. Check yourself

1. Define $\delta J$ in one line using $\varepsilon$ and $\eta$. Why must $\eta$ vanish at the endpoints?
2. State the fundamental lemma and sketch the bump-function proof in three sentences.
3. Derive the Euler–Lagrange equation from (3.3) in four lines. Where exactly is $\eta(a) = \eta(b) = 0$ used?
4. Write E–L for $L = \sqrt{1 + y'^2}$ and solve it.
5. Which first integral applies to $L = y^2\sqrt{1 + y'^2}$? To $L = x^2 y'^2$? To $L = y'^2 + x\,y$? Write it for each.
6. Derive the Beltrami identity. What is the Beltrami quantity for $L = \tfrac12 m\dot q^2 - V(q)$, and what is it called?
7. What is the natural boundary condition at a free endpoint, and what does it say for arc length?
8. State the Lagrange multiplier rule for an isoperimetric constraint. Is $\lambda$ a number or a function? And for a holonomic constraint?
9. What does $\partial J_h/\partial y_i$ equal, to leading order in $h$, when $J_h$ is the discretized functional? Why does this explain the numerical minimizer?
10. Why is the brachistochrone's start vertical? When does the solution dip below the endpoint?
11. State Legendre's condition. What is $L_{y'y'}$ for the brachistochrone, and what does its sign tell you?
12. What is the functional derivative $\delta J/\delta y$ for $J = \int L(x, y, y')\,dx$, and what is the E–L equation in that notation?

## 10. Pair with

- **Taylor, *Classical Mechanics*, ch. 6.** The cleanest short textbook treatment; §6.1–6.3 do the derivation and the brachistochrone at the level of this module, and §6.4 is the start of CM-2. Better than these notes on the mechanics motivation; read it before CM-2.
- **Gelfand & Fomin, *Calculus of Variations*, ch. 1.** The rigorous version: the fundamental lemma (and the du Bois-Reymond lemma of Exercise B3), the Euler equation, natural boundary conditions, multiple integrals (Exercise C3) and the variational derivative, all proved. Dense but short. Ch. 2 continues with several unknown functions, higher derivatives and the isoperimetric rule of §3.10. Better than these notes on precision; ch. 5 has the Jacobi condition we stated without proof.
- **Feynman Lectures vol. II, ch. 19, "The Principle of Least Action."** A transcribed special lecture, outside the course's regular sequence. No rigor, maximum physical intuition about why "stationary" and what the action is doing. Better than anything else on the *feel* of the subject; read it in one sitting after §3.13.
