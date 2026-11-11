# QM-6: Three Dimensions and the Hydrogen Atom: Separation of Variables, the Radial Equation, the Spectrum, Degeneracy

**Time:** ~12 hours · **Prerequisites:** QM-5, CM-5 (the effective potential), EM-3 (Legendre, separation of variables) · **Unlocks:** QM-8, SM-5 (context), atomic-physics pointers

> **Conventions in this module.** $\hbar$ explicit; `project.py` uses **atomic units** $\hbar = m_e = e = 4\pi\varepsilon_0 = 1$ (lengths in $a_0$, energies in hartree $E_h = 27.211$ eV $= 2\times13.606$ eV) outside its SI part. $m$ is the mass in the Schrödinger equation: strictly the reduced mass $\mu = m_em_p/(m_e+m_p)$ of CM-5 (3.2)–(3.3), replaced by $m_e$ where the $0.05\%$ difference does not matter; where it does (spectroscopy, muonic atoms) we write $\mu$. The magnetic quantum number is $m_\ell$ when it could be confused with a mass, $m$ otherwise. Bohr radius $a_0 = 4\pi\varepsilon_0\hbar^2/m_ee^2 = 0.052918$ nm; Bohr energies $E_n = E_1/n^2$, $E_1 = -13.606$ eV (infinite nuclear mass) or $-13.598$ eV (with $\mu$). Wavefunctions $\psi_{n\ell m}(r,\theta,\phi) = R_{n\ell}(r)Y_\ell^m(\theta,\phi)$ with QM-5's Condon–Shortley $Y_\ell^m$ (SciPy `sph_harm_y`). $u(r) = rR(r)$, normalized by $\int_0^\infty|u|^2dr = 1$. Radial quantum number $n_r$ = number of radial nodes; principal quantum number $n = n_r + \ell + 1$. **Associated Laguerre polynomials** $L_k^{(\alpha)}(x)$ in the modern convention used by Griffiths & Schroeter (3rd ed.), Zwiebach, and `scipy.special.genlaguerre`: $L_k^{(\alpha)}(0) = \binom{k+\alpha}{k}$. Older texts (Griffiths 2nd ed., many chemistry books) use a polynomial that is $(-1)^\alpha(k+\alpha)!$ times this one; §3.7 and Project Part 2 map them. CM-5 wrote the classical angular momentum as $\ell$; here $\ell$ is the quantum number and the classical magnitude is $L$.

## 0. Where this sits

You need QM-5's spherical harmonics ($\hat L^2Y_\ell^m = \hbar^2\ell(\ell+1)Y_\ell^m$ and the Laplacian split QM-5 (3.20)), QM-3's two techniques for 1D bound states (the series-termination method of QM-3 §3.3 and the finite-difference solver of QM-3 project Part 1), CM-5's reduction of a two-body central problem to one radial equation with an effective potential (CM-5 (3.8)), and EM-3's Frobenius argument for Legendre's equation (EM-3 §3.5, steps 3–5).

This module solves the only realistic atom with an exact solution. QM-8 perturbs it (fine structure, Zeeman and Stark effects, and the selection rules quoted in §3.11); QM-7 adds the exclusion principle that turns hydrogen's orbitals into §3.14's periodic table; SM-5 uses the hydrogenic level structure as context for ionization equilibrium.

## 1. The problem that forced this idea

By 1885 hydrogen's visible spectrum was four sharp lines at 656, 486, 434 and 410 nm. Johann Balmer, a Basel schoolteacher, found that all four fit $\lambda = B\,n^2/(n^2-4)$ for $n = 3, 4, 5, 6$, and Rydberg rewrote it as $1/\lambda = R(1/n_f^2 - 1/n_i^2)$. Nothing in nineteenth-century physics produced integers from an atom. Worse, after Rutherford (1911) the atom was a nucleus with an electron around it, and EM-8 Example 4.4 showed that such an electron radiates and spirals in within about $1.6\times10^{-11}$ s. A theory of the atom had to explain two facts at once: atoms are stable, and they emit only certain frequencies.

Bohr's 1913 model (QM-1 Theorem 3.1) got the numbers by decree: allowed circular orbits with $L = n\hbar$ that do not radiate. It produced $R$ from $m_e$, $e$, $h$ to within experimental error, and it was wrong about almost everything else. The ground state has $L = 0$ in reality, which no Bohr orbit has; the model could not say why some lines are bright and others absent; it failed outright for helium.

In early 1926 Pauli derived the hydrogen spectrum from Heisenberg's matrix mechanics using the quantum version of the Laplace–Runge–Lenz vector (CM-4 §3.8 gives the history). Weeks later Schrödinger's first paper on wave mechanics solved his new equation for the Coulomb potential by the method of this module and recovered Balmer's formula. His route made the integers come from something physicists already understood: a wave equation with boundary conditions, quantized the way a drumhead is quantized. That is what got wave mechanics taken seriously within months.

The solution gives more than the energies. Each level has a set of states labeled $(n, \ell, m)$, and the count of them, $n^2$ per level, $2n^2$ with spin, is the length of the rows of the periodic table (§3.14). The chemistry of carbon and the color of a neon sign sit downstream of the series termination in §3.5.

## 2. Intuition

**It is CM-5 with $\hbar$.** Two bodies attract by a central force. In CM-5 you separated the center of mass, used conservation of $\mathbf L$ to put the motion in a plane, and were left with one radial equation in an effective potential $V(r) + L^2/2\mu r^2$. Quantum mechanics does the same three moves. The center of mass separates exactly. Angular momentum is conserved ($[\hat{\mathbf L}, \hat H] = 0$, QM-5 B1), so states can carry definite $\hat L^2$ and $\hat L_z$; their angular dependence is a spherical harmonic $Y_\ell^m$, and $L^2$ becomes $\hbar^2\ell(\ell+1)$. What is left is a one-dimensional Schrödinger equation on the half-line $r > 0$, for the function $u = rR$, in the potential
$$V_{\rm eff}(r) = -\frac{e^2}{4\pi\varepsilon_0r} + \frac{\hbar^2\ell(\ell+1)}{2mr^2}.$$
This is QM-3 again: a 1D well, bound states where standing waves fit.

**Sketch the well.** For $\ell = 0$ the well is the bare Coulomb $-1/r$, bottomless at the origin. For $\ell \ge 1$ the centrifugal $+1/r^2$ wins at small $r$ and builds a wall, so the well has a minimum at some radius that grows like $\ell(\ell+1)$:
```
 V_eff      l=2 \
   |          \     l=1
   |   l=0     \____...........  0
   |    |   ___/       ---
   |    |  /   ---
   |    | /
   |    |/   (l=0 has no wall; -1/r all the way down)
   +--------------------------> r
```
States with larger $\ell$ live farther out.

**Count nodes.** In any 1D well the $k$-th state has $k-1$ nodes, and more nodes mean more kinetic energy (QM-3's node theorem). In each $\ell$-well the radial states are labeled by their node count $n_r = 0, 1, 2, \ldots$, and the energy rises with $n_r$. It also rises with $\ell$, because the centrifugal wall squeezes the state. So guess: the energy depends on a combination of $n_r$ and $\ell$.

**Guess the formula.** Bohr's levels $E \propto -1/n^2$ are right (Balmer says so). The ground state has no nodes and no angular momentum: $n_r = 0$, $\ell = 0$, $n = 1$. The simplest guess consistent with "each node costs about as much as each unit of $\ell$" is $n = n_r + \ell + 1$, which makes $E$ depend on $n_r + \ell$ only. Then a level $n$ contains $\ell = 0, 1, \ldots, n-1$, each with $2\ell+1$ values of $m$: $\sum_{\ell=0}^{n-1}(2\ell+1) = n^2$ states. That is the theorem. The surprise is that it is exact: a state with one radial node and no angular momentum has *exactly* the energy of a state with no node and $\ell = 1$.

**Why that is a surprise.** In CM-5 (3.15) the energy of a Kepler orbit depends only on the semi-major axis $a$, not on the eccentricity: a circle and a needle-thin ellipse with the same $a$ have the same energy. Eccentricity is the classical shadow of $\ell$ at fixed $n$. The Laplace–Runge–Lenz vector (CM-4 §3.8) is the extra conserved quantity that makes this happen, and only the $1/r$ potential has it. Any other central potential has a conserved $\mathbf L$ but no LRL vector, and its levels depend on $n_r$ and $\ell$ separately (§3.12).

**Where the integers come from.** At large $r$ a bound state decays like $e^{-\kappa r}$ with $\hbar^2\kappa^2/2m = -E$. Near $r = 0$ it behaves like $r^{\ell+1}$. In between it is a polynomial. A polynomial times $e^{-\kappa r}$ satisfies the equation only for special $\kappa$; for any other $\kappa$ the series that replaces the polynomial secretly grows like $e^{+2\kappa r}$. This is QM-3's oscillator argument and EM-3's Legendre argument, run a third time.

## 3. Formal development

### 3.1 Separation of variables for any central potential

The time-independent Schrödinger equation for one particle of mass $m$ in a potential $V(r)$ depending only on $r = |\mathbf r|$ is $-\frac{\hbar^2}{2m}\nabla^2\psi + V(r)\psi = E\psi$. With QM-5 (3.20), $\nabla^2 = \frac{1}{r^2}\partial_r(r^2\partial_r) - \hat L^2/\hbar^2r^2$:
$$
-\frac{\hbar^2}{2m}\left[\frac{1}{r^2}\frac{\partial}{\partial r}\Big(r^2\frac{\partial\psi}{\partial r}\Big)\right] + \frac{\hat L^2}{2mr^2}\psi + V(r)\psi = E\psi. \tag{3.1}
$$

▶ **Try it first.** Put $\psi = R(r)Y(\theta,\phi)$ into (3.1), multiply by $2mr^2/\hbar^2RY$, and collect everything that depends on $r$ on one side.

**Derivation.**
1. Substituting, $\hat L^2$ acts only on $Y$ and $\partial_r$ only on $R$. Multiplying by $2mr^2/(\hbar^2RY)$:
$$\underbrace{\frac{1}{R}\frac{d}{dr}\Big(r^2\frac{dR}{dr}\Big) - \frac{2mr^2}{\hbar^2}\big(V(r) - E\big)}_{\text{function of } r} = \underbrace{\frac{1}{Y}\frac{\hat L^2Y}{\hbar^2}}_{\text{function of } \theta, \phi}.$$
2. A function of $r$ alone equal to a function of angles alone is a constant. Call it $\lambda$. The angular equation is $\hat L^2Y = \hbar^2\lambda Y$.
3. QM-5 §3.3 solved this: the single-valued solutions regular at the poles are $Y = Y_\ell^m$ with $\lambda = \ell(\ell+1)$, $\ell = 0, 1, 2, \ldots$, $m = -\ell, \ldots, \ell$. (EM-3 §3.5 did the $m = 0$ case by Frobenius.)
4. The radial equation is
$$
\frac{d}{dr}\Big(r^2\frac{dR}{dr}\Big) - \frac{2mr^2}{\hbar^2}\big(V(r) - E\big)R = \ell(\ell+1)R. \tag{3.2}
$$
5. Substitute $R = u/r$. Then $R' = u'/r - u/r^2$, $r^2R' = ru' - u$, and $(r^2R')' = u' + ru'' - u' = ru''$. Dividing (3.2) by $r$ and multiplying by $-\hbar^2/2m$:
$$
\boxed{-\frac{\hbar^2}{2m}\frac{d^2u}{dr^2} + \left[V(r) + \frac{\hbar^2}{2m}\frac{\ell(\ell+1)}{r^2}\right]u = Eu.} \tag{3.3}
$$
6. Normalization: $\int|\psi|^2d^3r = \int_0^\infty|R|^2r^2dr\int|Y_\ell^m|^2d\Omega$. The $Y_\ell^m$ are normalized on the sphere, and $|R|^2r^2 = |u|^2$, so
$$
\int_0^\infty|u(r)|^2\,dr = 1. \tag{3.4}
$$
∎

■ **What just happened.** For *every* central potential, the 3D problem is a 1D problem (3.3) on $r > 0$ with the effective potential $V_{\rm eff} = V + \hbar^2\ell(\ell+1)/2mr^2$, which is CM-5 (3.8) with $L^2\to\hbar^2\ell(\ell+1)$. The quantity $|u|^2dr$ is the probability of finding the particle in the shell between $r$ and $r + dr$, whatever the direction. The energy cannot depend on $m$, since $m$ does not appear in (3.3): that is the $(2\ell+1)$-fold degeneracy forced by rotational symmetry (QM-5 §3.2). Only the potential can make $E$ depend on $\ell$.

**$\ell(\ell+1)$, not $\ell^2$.** Classically the centrifugal term is $L^2/2mr^2$ and a state with $L_z = L$ is a planar orbit. Quantum mechanically the maximum $L_z^2$ is $\hbar^2\ell^2$ but $\langle L^2\rangle = \hbar^2\ell(\ell+1)$: even in the state $m = \ell$, $\langle L_x^2 + L_y^2\rangle = \hbar^2\ell$ (QM-5 B2), because the components you did not fix still fluctuate. The extra $\hbar^2\ell$ is real rotational kinetic energy, it is in the Hamiltonian, and it changes the levels. For $\ell = 0$ all three components are sharp at zero and there is no centrifugal term at all; this is the state Bohr's model had no room for.

### 3.2 The boundary condition at the origin

▶ **Try it first.** Near $r = 0$, assume $V$ is less singular than $1/r^2$ (true for Coulomb). Try $u \approx r^s$ in (3.3) and find the allowed $s$.

**Derivation.**
1. For small $r$ the centrifugal term dominates $V$ and $E$, so $u'' \approx \ell(\ell+1)u/r^2$. With $u = r^s$: $s(s-1) = \ell(\ell+1)$, so $s = \ell+1$ or $s = -\ell$.
2. $s = -\ell$ with $\ell \ge 1$: $\int_0|u|^2dr = \int_0r^{-2\ell}dr$ diverges. Not normalizable; rejected.
3. $s = -\ell = 0$, i.e. $\ell = 0$ and $u(0) \neq 0$: then $R \approx u(0)/r$, which is normalizable. But $\nabla^2(1/r) = -4\pi\delta^3(\mathbf r)$ (EM-1/EM-2's point charge), so $-\frac{\hbar^2}{2m}\nabla^2\psi$ contains a delta function at the origin that nothing on the right of the Schrödinger equation can match. Rejected.
4. Therefore
$$
u(0) = 0, \qquad u(r) \propto r^{\ell+1}\ \text{ as } r\to0. \tag{3.5}
$$
∎

■ **What just happened.** The half-line problem (3.3) behaves as if there were an infinite wall at $r = 0$. That is why every $\ell = 0$ radial function $u$ looks like an *odd* 1D wavefunction, and why the finite-difference solver of Project Part 1 is QM-3's solver with Dirichlet ends at $r = 0$ and $r = r_{\max}$. Note that $R = u/r \to$ a nonzero constant for $\ell = 0$: s-electrons are found at the nucleus (§4.6), and for $\ell \ge 1$, $R\propto r^\ell$ vanishes there.

### 3.3 Warm-up: the infinite spherical well

$V = 0$ for $r < a$, $\infty$ for $r > a$. Set $k = \sqrt{2mE}/\hbar$.

**$\ell = 0$.** (3.3) is $u'' = -k^2u$ with $u(0) = u(a) = 0$. This is QM-3 (3.2) on $[0, a]$: $u = \sqrt{2/a}\,\sin(n\pi r/a)$, $E_{n0} = n^2\pi^2\hbar^2/2ma^2$.

**$\ell \ge 1$.** The regular solution of (3.2) with $V = 0$ is the **spherical Bessel function** $R = j_\ell(kr)$ (stated, not derived; $j_0(x) = \sin x/x$, $j_1(x) = \sin x/x^2 - \cos x/x$, and the irregular $n_\ell(kr)\sim r^{-\ell-1}$ is excluded by §3.2). The wall requires $j_\ell(ka) = 0$. If $\beta_{n\ell}$ is the $n$-th positive zero of $j_\ell$,
$$
E_{n\ell} = \frac{\hbar^2\beta_{n\ell}^2}{2ma^2}. \tag{3.6}
$$
For $\ell = 1$, $j_1(x) = 0$ means $\tan x = x$, first root $\beta_{11} = 4.4934$. The ordering of $\beta$'s is $1s\,(3.1416) < 1p\,(4.4934) < 1d\,(5.7635) < 2s\,(6.2832) < 1f\,(6.9879) < 2p\,(7.7253)$ (Project Part 6 computes them with `scipy.special.spherical_jn`). Here the label $n$ counts zeros within each $\ell$, as nuclear physicists do. The levels depend on $n$ and $\ell$ separately: no accidental degeneracy. This ordering, plus spin–orbit coupling, is the start of the nuclear shell model.

### 3.4 Hydrogen: the two-body reduction and the dimensionless equation

The electron (mass $m_e$, charge $-e$) and proton ($m_p$, $+e$) interact through $V = -e^2/(4\pi\varepsilon_0|\mathbf r_e - \mathbf r_p|)$. With CM-5's coordinates (3.2), the kinetic energy $\frac{p_e^2}{2m_e} + \frac{p_p^2}{2m_p} = \frac{P^2}{2M} + \frac{p^2}{2\mu}$, where $\mathbf P$ is conjugate to the center of mass $\mathbf R$ and $\mathbf p$ to $\mathbf r = \mathbf r_e - \mathbf r_p$. The Hamiltonian is a sum of a free particle in $\mathbf R$ and a one-body problem in $\mathbf r$, so $\Psi = e^{i\mathbf K\cdot\mathbf R}\psi(\mathbf r)$ and $\psi$ obeys (3.1) with $m\to\mu$ and
$$
V(r) = -\frac{e^2}{4\pi\varepsilon_0}\frac{1}{r}. \tag{3.7}
$$
Since $\mu = m_e/(1 + m_e/m_p) = 0.99946\,m_e$, we write $m$ and remember the correction when comparing with spectra (§4.5, Exercise B6).

▶ **Try it first.** For a bound state ($E < 0$) define $\kappa = \sqrt{-2mE}/\hbar$ and $\rho = \kappa r$. Divide (3.3) by $E$ and find what dimensionless constant multiplies $1/\rho$.

**Derivation.**
1. Divide (3.3) by $-E = \hbar^2\kappa^2/2m$:
$$\frac{1}{\kappa^2}\frac{d^2u}{dr^2} = \left[1 - \frac{me^2}{2\pi\varepsilon_0\hbar^2\kappa}\frac{1}{\kappa r} + \frac{\ell(\ell+1)}{(\kappa r)^2}\right]u.$$
2. With $\rho = \kappa r$ and $\rho_0 \equiv me^2/(2\pi\varepsilon_0\hbar^2\kappa)$:
$$
\frac{d^2u}{d\rho^2} = \left[1 - \frac{\rho_0}{\rho} + \frac{\ell(\ell+1)}{\rho^2}\right]u, \qquad \rho_0 = \frac{me^2}{2\pi\varepsilon_0\hbar^2\kappa} = \frac{2}{\kappa a_0}, \qquad a_0 = \frac{4\pi\varepsilon_0\hbar^2}{me^2}. \tag{3.8}
$$
∎

■ **What just happened.** All the physics constants collapsed into one number, $\rho_0$, which is a disguised energy: $\rho_0 = 2/(\kappa a_0)$ and $E = -\hbar^2\kappa^2/2m$, so $E = -\frac{\hbar^2}{2ma_0^2}\cdot\frac{4}{\rho_0^2}$. Finding the allowed $\rho_0$ is finding the spectrum. The length $a_0$ appeared without being put in: it is the only length you can build from $\hbar$, $m$, and $e^2/4\pi\varepsilon_0$.

### 3.5 Asymptotics, the ansatz, the series, and termination

▶ **Try it first.** Find the behavior of $u$ as $\rho\to\infty$ and as $\rho\to0$ from (3.8). Then write $u = \rho^{\ell+1}e^{-\rho}v(\rho)$ and try to guess what equation $v$ obeys.

**Derivation (asymptotics).**
1. $\rho\to\infty$: the bracket in (3.8) tends to 1, so $u'' \approx u$, $u \approx Ae^{-\rho} + Be^{\rho}$. Normalizability kills $B$: $u \sim e^{-\rho}$.
2. $\rho\to0$: the $\ell(\ell+1)/\rho^2$ term dominates and §3.2 gives $u\sim\rho^{\ell+1}$.
3. Factor both out:
$$
u(\rho) = \rho^{\ell+1}e^{-\rho}v(\rho). \tag{3.9}
$$
No information has been lost; $v$ is just $u$ divided by known functions. The hope is that $v$ is simple.

**Derivation (the equation for $v$).**
4. $\dfrac{du}{d\rho} = \rho^\ell e^{-\rho}\big[(\ell + 1 - \rho)v + \rho v'\big]$.
5. Differentiate again and collect:
$$\frac{d^2u}{d\rho^2} = \rho^\ell e^{-\rho}\Big[\rho v'' + 2(\ell + 1 - \rho)v' + \Big(\rho - 2(\ell+1) + \frac{\ell(\ell+1)}{\rho}\Big)v\Big].$$
6. The right side of (3.8) is $\rho^\ell e^{-\rho}\big[\rho - \rho_0 + \ell(\ell+1)/\rho\big]v$. Set them equal; the $\rho v$ and $\ell(\ell+1)v/\rho$ terms cancel:
$$
\rho\frac{d^2v}{d\rho^2} + 2(\ell + 1 - \rho)\frac{dv}{d\rho} + \big[\rho_0 - 2(\ell+1)\big]v = 0. \tag{3.10}
$$
(Project Part 3 checks this with SymPy.)

**Derivation (the recursion).**
7. Try $v = \sum_{j\ge0}c_j\rho^j$. Term by term, re-indexing so every sum runs over $\rho^j$:
   - $\rho v'' = \sum_j j(j+1)c_{j+1}\rho^j$,
   - $2(\ell+1)v' = \sum_j 2(\ell+1)(j+1)c_{j+1}\rho^j$,
   - $-2\rho v' = \sum_j(-2j)c_j\rho^j$,
   - $[\rho_0 - 2(\ell+1)]v = \sum_j[\rho_0 - 2(\ell+1)]c_j\rho^j$.
8. The coefficient of $\rho^j$ must vanish: $(j+1)(j + 2\ell + 2)c_{j+1} = [2(j + \ell + 1) - \rho_0]c_j$, so
$$
\boxed{c_{j+1} = \frac{2(j + \ell + 1) - \rho_0}{(j+1)(j + 2\ell + 2)}\,c_j.} \tag{3.11}
$$
∎

**Theorem 3.1 (Termination).** The solution (3.9) is normalizable if and only if the series for $v$ terminates, which happens if and only if $\rho_0 = 2n$ with $n = j_{\max} + \ell + 1$ a positive integer and $j_{\max}\ge0$ the degree of $v$.

*Proof.*
1. Suppose the series does not terminate. For large $j$, (3.11) gives $c_{j+1}/c_j \to 2/(j+1)$, the ratio of coefficients of $e^{2\rho} = \sum(2\rho)^j/j!$.
2. More carefully: for $j$ beyond some $J$, $c_{j+1}/c_j > (2-\epsilon)/(j+1)$ for any fixed $\epsilon > 0$, so the tail of $v$ dominates a constant times $e^{(2-\epsilon)\rho}$. Then $u = \rho^{\ell+1}e^{-\rho}v$ grows like $e^{(1-\epsilon)\rho}$. That is the rejected $Be^{\rho}$ of step 1, sneaking back in. Not normalizable.
3. So the series must stop: some $c_{j_{\max}+1} = 0$ with $c_{j_{\max}}\ne0$. By (3.11) that requires $2(j_{\max} + \ell + 1) = \rho_0$. Conversely, if $\rho_0 = 2(j_{\max}+\ell+1)$ every later coefficient vanishes, $v$ is a polynomial, and $u = \rho^{\ell+1}e^{-\rho}\times$polynomial is normalizable. ∎

■ **What just happened.** Three moves, the same three as QM-3 §3.3 (oscillator: Gaussian asymptotics, Hermite recursion (3.16), termination gives $E = (n+\frac12)\hbar\omega$) and EM-3 §3.5 (Legendre: regularity at the poles, recursion (3.14), termination gives integer $\ell$). Asymptotics peel off the known behavior; the recursion is mechanical; the non-terminating series is a disguised version of the solution you threw away; termination is quantization. Recognize this pattern and the rest of the hydrogen solution is bookkeeping. Project Part 3 runs (3.11) at $\rho_0 = 4$ and at $\rho_0 = 4.01$ and watches the second one track $e^{2\rho}$.

### 3.6 The Bohr formula, the quantum numbers, and the degeneracy

**Derivation.**
1. Theorem 3.1: $\rho_0 = 2n$, and by (3.8) $\rho_0 = 2/(\kappa a_0)$, so
$$
\kappa_n = \frac{1}{na_0}. \tag{3.12}
$$
2. $E = -\hbar^2\kappa^2/2m$:
$$
\boxed{E_n = -\frac{\hbar^2}{2ma_0^2}\frac{1}{n^2} = -\frac{m}{2\hbar^2}\Big(\frac{e^2}{4\pi\varepsilon_0}\Big)^2\frac{1}{n^2} = \frac{E_1}{n^2}, \qquad E_1 = -13.606\ \text{eV}.} \tag{3.13}
$$
3. Numbers: $a_0 = 4\pi\varepsilon_0\hbar^2/m_ee^2 = 5.2918\times10^{-11}$ m, and $\hbar^2/(2m_ea_0^2) = 13.606$ eV. With $\mu$ instead of $m_e$: $E_1 = -13.598$ eV (Project Part 8). A useful form: $E_1 = -\frac12\alpha^2m_ec^2$ with the fine-structure constant $\alpha = e^2/(4\pi\varepsilon_0\hbar c) = 1/137.036$, and $a_0 = \hbar/(\alpha m_ec)$.
4. **Allowed quantum numbers.** $n = j_{\max} + \ell + 1$ with $j_{\max}\ge0$, so for given $n$: $\ell = 0, 1, \ldots, n - 1$, and for each $\ell$, $m = -\ell, \ldots, \ell$. The radial node count is $n_r = j_{\max} = n - \ell - 1$ (Project Part 2 counts them: a degree-$j_{\max}$ polynomial with all roots real and positive, a Laguerre property).
5. **Degeneracy.**
$$
d_n = \sum_{\ell=0}^{n-1}(2\ell + 1) = 2\cdot\frac{(n-1)n}{2} + n = n^2, \qquad 2n^2 \text{ with spin.} \tag{3.14}
$$
∎

■ **What just happened.** $E$ depends on $n = n_r + \ell + 1$ only: §2's guess, now exact. Of the $n^2$-fold degeneracy, the factor $2\ell + 1$ is rotational symmetry (any $V(r)$), and the sum over $\ell$ is special to $1/r$ (§3.12). Spectroscopic notation: $\ell = 0, 1, 2, 3$ is s, p, d, f, so the states are $1s$; $2s, 2p$; $3s, 3p, 3d$; and so on.

**The quantum numbers as a complete set of commuting observables.** QM-4 §3.7 said a set of commuting observables is complete when their joint eigenvalues label a basis without ambiguity. Here $\hat H$, $\hat L^2$, $\hat L_z$ commute ($[\hat{\mathbf L}, \hat H] = 0$ and $[\hat L^2, \hat L_z] = 0$, QM-5). $\hat H$ alone is not complete (its eigenvalue $E_n$ has $n^2$ states); adding $\hat L^2$ and $\hat L_z$ makes it so for spinless hydrogen: $|n\ell m\rangle$ is the unique simultaneous eigenvector. With spin, add $\hat S_z$ and label $|n\ell m_\ell m_s\rangle$; QM-8's fine structure prefers the coupled set $\hat H, \hat L^2, \hat J^2, \hat J_z$ (QM-5 §3.10).

### 3.7 The radial functions and Laguerre polynomials

▶ **Try it first.** Run (3.11) for $n = 1$, $\ell = 0$ ($\rho_0 = 2$) and for $n = 2$, $\ell = 1$ ($\rho_0 = 4$). What is $v$ in each case?

In both cases $j_{\max} = 0$, so $v = c_0$: the states with $\ell = n - 1$ (no radial nodes) have $u \propto \rho^ne^{-\rho}$, $\rho = r/na_0$. For general $(n, \ell)$, Eq. (3.10) with $\rho_0 = 2n$ and $x = 2\rho$ becomes $xv'' + (2\ell + 2 - x)v' + (n - \ell - 1)v = 0$, which is the associated Laguerre equation $xy'' + (\alpha + 1 - x)y' + ky = 0$ with $\alpha = 2\ell+1$ and $k = n - \ell - 1$. So $v\propto L^{(2\ell+1)}_{n-\ell-1}(2r/na_0)$ and
$$
R_{n\ell}(r) = \sqrt{\Big(\frac{2}{na_0}\Big)^3\frac{(n-\ell-1)!}{2n\,(n+\ell)!}}\;e^{-r/na_0}\Big(\frac{2r}{na_0}\Big)^\ell L^{(2\ell+1)}_{n-\ell-1}\Big(\frac{2r}{na_0}\Big). \tag{3.15}
$$
The normalization constant is quoted (it follows from a standard Laguerre integral, done in Griffiths & Schroeter §4.2); Project Part 2 checks (3.15) against the finite-difference eigenvectors and checks $\int u^2 = 1$ by quadrature for every $n\le4$.

**The convention trap.** In (3.15) $L^{(\alpha)}_k$ is the modern polynomial: `scipy.special.genlaguerre(k, alpha)`, with $L^{(1)}_1(x) = 2 - x$. Older books define $L_q(x) = e^x\frac{d^q}{dx^q}(e^{-x}x^q)$ and $L^p_{q-p}(x) = (-1)^p\frac{d^p}{dx^p}L_q(x)$, which equals $q!\,L^{(p)}_{q-p}(x)$ in the modern notation. With that polynomial the prefactor carries $[(n+\ell)!]^3$ instead of $(n+\ell)!$; the product is the same function. Mix a modern polynomial with an old prefactor and every normalization check fails by factorials. Project Part 2 builds the old polynomial in SymPy and asserts the factor $q!$.

The first few, with $a = a_0$:
$$
\begin{aligned}
R_{10} &= 2a^{-3/2}e^{-r/a}, \\
R_{20} &= \tfrac{1}{\sqrt2}a^{-3/2}\Big(1 - \frac{r}{2a}\Big)e^{-r/2a}, \qquad
R_{21} = \tfrac{1}{\sqrt{24}}a^{-3/2}\frac{r}{a}e^{-r/2a}, \\
R_{30} &= \tfrac{2}{\sqrt{27}}a^{-3/2}\Big(1 - \frac{2r}{3a} + \frac{2r^2}{27a^2}\Big)e^{-r/3a}.
\end{aligned} \tag{3.16}
$$

**Radial probability density.** $P(r) = r^2|R_{n\ell}|^2 = |u_{n\ell}|^2$ has $n - \ell - 1$ zeros at $r > 0$ and $n - \ell$ humps. Plotted (Project figure `radial_densities.png`): $1s$ is one hump peaked at $a_0$; $2s$ has a small inner hump and a big outer one with a node at $2a_0$; $2p$ is one hump peaked at $4a_0$; $3s$ has three humps, $3p$ two, $3d$ one, peaked at $9a_0$. The outermost hump carries most of the probability and sits near $n^2a_0$, but the s-states keep a small amount of probability close to the nucleus. That inner hump is what §3.14's screening argument is about.

**Expectation values (quoted general formulas).** For the state $n\ell$,
$$
\langle r\rangle = \frac{a_0}{2}\big[3n^2 - \ell(\ell+1)\big], \qquad \langle r^2\rangle = \frac{a_0^2n^2}{2}\big[5n^2 + 1 - 3\ell(\ell+1)\big], \qquad \Big\langle\frac1r\Big\rangle = \frac{1}{n^2a_0}. \tag{3.17}
$$
§3.8 derives the $n = 1$ values; Exercise C2 derives $\langle1/r\rangle$ for every state by the Feynman–Hellmann theorem; the other two follow from Kramers' relation (in the problems of Griffiths & Schroeter §4.2). Project Part 4 checks all three by quadrature for every $n\le4$. Note $\langle1/r\rangle \neq 1/\langle r\rangle$.

### 3.8 The ground state

$$
\psi_{100}(r) = \frac{1}{\sqrt{\pi a_0^3}}\,e^{-r/a_0}, \qquad E_1 = -13.6\ \text{eV}. \tag{3.18}
$$
(From $R_{10}Y_0^0 = 2a^{-3/2}e^{-r/a}\cdot\frac{1}{\sqrt{4\pi}}$.) Spherically symmetric, no nodes, largest at $r = 0$.

▶ **Try it first.** With $P(r) = 4r^2e^{-2r/a}/a^3$ and $\int_0^\infty r^ke^{-2r/a}dr = k!(a/2)^{k+1}$, compute $\langle r\rangle$, $\langle1/r\rangle$, and the most probable $r$.

**Derivation.**
1. $\langle r^k\rangle = \frac{4}{a^3}\int_0^\infty r^{k+2}e^{-2r/a}dr = \frac{4}{a^3}(k+2)!\Big(\frac a2\Big)^{k+3} = \frac{(k+2)!}{2^{k+1}}a^k$.
2. $k = 1$: $\langle r\rangle = \frac{3!}{4}a = \frac32a_0$. $k = 2$: $\langle r^2\rangle = \frac{4!}{8}a^2 = 3a_0^2$, so $\sigma_r = \sqrt{3 - 9/4}\,a_0 = 0.866\,a_0$. $k = -1$: $\langle1/r\rangle = \frac{1!}{1}a^{-1} = 1/a_0$.
3. Most probable: $\frac{d}{dr}(r^2e^{-2r/a}) = (2r - 2r^2/a)e^{-2r/a} = 0$ at $r = a_0$. Bohr's radius is the *peak* of $P(r)$; the mean is 50% farther out because the tail is long.
4. Potential energy: $\langle V\rangle = -\frac{e^2}{4\pi\varepsilon_0}\langle\frac1r\rangle = -\frac{e^2}{4\pi\varepsilon_0a_0} = 2E_1 = -27.2$ eV, using $E_1 = -e^2/(8\pi\varepsilon_0a_0)$, which is (3.13) with one factor of $a_0$ expanded.
5. Kinetic energy, directly: $u = 2a^{-3/2}re^{-r/a}$, $u' = 2a^{-3/2}(1 - r/a)e^{-r/a}$, and with no centrifugal term $\langle T\rangle = \frac{\hbar^2}{2m}\int_0^\infty|u'|^2dr$ (integrate $-u\,u''$ by parts; the boundary term vanishes). $\int_0^\infty(1-x)^2e^{-2x}dx = \frac12 - \frac12 + \frac14 = \frac14$, so $\int|u'|^2dr = 4a^{-3}\cdot\frac a4 = a^{-2}$ and $\langle T\rangle = \hbar^2/2ma_0^2 = -E_1 = +13.6$ eV.
6. Check: $\langle T\rangle + \langle V\rangle = 13.6 - 27.2 = -13.6$ eV $= E_1$. ∎

**Theorem 3.2 (Virial theorem).** In any stationary state of $\hat H = \hat p^2/2m + V(\mathbf r)$, $2\langle T\rangle = \langle\mathbf r\cdot\nabla V\rangle$. For $V\propto r^k$, $2\langle T\rangle = k\langle V\rangle$; for Coulomb ($k = -1$), $\langle T\rangle = -E$ and $\langle V\rangle = 2E$.

*Proof.* QM-4 (3.20) with $\hat Q = \hat{\mathbf r}\cdot\hat{\mathbf p}$: in a stationary state $\frac{d}{dt}\langle Q\rangle = 0$, so $\langle[\hat H, \hat{\mathbf r}\cdot\hat{\mathbf p}]\rangle = 0$. Using $[\hat p^2, \hat x_i] = -2i\hbar\hat p_i$ and $[V, \hat p_i] = i\hbar\,\partial_iV$, $[\hat H, \hat x_i\hat p_i] = [\hat H, \hat x_i]\hat p_i + \hat x_i[\hat H, \hat p_i] = -\frac{i\hbar}{m}\hat p_i^2 + i\hbar\,\hat x_i\partial_iV$. Summing over $i$, $0 = -2i\hbar\langle T\rangle + i\hbar\langle\mathbf r\cdot\nabla V\rangle$. For $V\propto r^k$, $\mathbf r\cdot\nabla V = kV$ (Euler's theorem for homogeneous functions). With $E = \langle T\rangle + \langle V\rangle$ and $k = -1$: $\langle T\rangle = -\frac12\langle V\rangle$, so $E = \frac12\langle V\rangle$. ∎

This is CM-5's (and CM-4's) virial theorem for time-averages over a Kepler orbit, with the time average replaced by the stationary-state average. It holds for every $n\ell m$, so $\langle V\rangle_n = 2E_n$, which is (3.17)'s $\langle1/r\rangle = 1/n^2a_0$ in another guise.

### 3.9 Orbital shapes, and real orbitals

$|\psi_{n\ell m}|^2 = |R_{n\ell}(r)|^2|Y_\ell^m(\theta,\phi)|^2$, and $|Y_\ell^m|^2 \propto |P_\ell^m(\cos\theta)|^2$ does not depend on $\phi$: every complex orbital is a figure of revolution about $z$. $2p_0$ ($m = 0$) is two lobes along $\pm z$ ($|Y_1^0|^2\propto\cos^2\theta$); $2p_{\pm1}$ is a doughnut in the $xy$-plane ($\propto\sin^2\theta$); $3d_0$ is two lobes along $z$ plus a ring ($\propto(3\cos^2\theta - 1)^2$). Project figure `orbitals.png` shows cross-sections.

Chemists want orbitals pointing along $x$ and $y$ because bonds point in definite directions. Since $E$ does not depend on $m$, any combination of $m = \pm1$ states with the same $n, \ell$ is also an energy eigenstate.

▶ **Try it first.** With $Y_1^{\pm1} = \mp\sqrt{3/8\pi}\,\sin\theta\,e^{\pm i\phi}$ (Condon–Shortley, QM-5 (3.22)), find combinations proportional to $x/r = \sin\theta\cos\phi$ and $y/r = \sin\theta\sin\phi$.

**Derivation.**
1. $Y_1^{-1} - Y_1^1 = \sqrt{3/8\pi}\,\sin\theta\,(e^{-i\phi} + e^{i\phi}) = 2\sqrt{3/8\pi}\,\sin\theta\cos\phi$.
2. $Y_1^{-1} + Y_1^1 = \sqrt{3/8\pi}\,\sin\theta\,(e^{-i\phi} - e^{i\phi}) = -2i\sqrt{3/8\pi}\,\sin\theta\sin\phi$.
3. Divide by $\sqrt2$ to normalize (the $Y$'s are orthonormal, so $\|Y_1^{-1}\mp Y_1^1\|^2 = 2$):
$$
p_x \equiv \frac{Y_1^{-1} - Y_1^1}{\sqrt2} = \sqrt{\frac{3}{4\pi}}\frac{x}{r}, \qquad p_y \equiv \frac{i(Y_1^{-1} + Y_1^1)}{\sqrt2} = \sqrt{\frac{3}{4\pi}}\frac{y}{r}, \qquad p_z \equiv Y_1^0 = \sqrt{\frac{3}{4\pi}}\frac zr. \tag{3.19}
$$
4. The same with $m = \pm2$: $d_{xy} \equiv \frac{Y_2^2 - Y_2^{-2}}{i\sqrt2} = \sqrt{\frac{15}{4\pi}}\frac{xy}{r^2}$, using $Y_2^{\pm2} = \sqrt{15/32\pi}\sin^2\theta\,e^{\pm2i\phi}$ and $\sin^2\theta\sin2\phi = 2xy/r^2$. ∎

■ **What just happened.** Real orbitals are a change of basis inside a degenerate eigenspace (LA-8: any orthonormal basis of an eigenspace is an eigenbasis). They are eigenstates of $\hat H$ and $\hat L^2$ but **not** of $\hat L_z$: $p_x$ is an equal superposition of $m = +1$ and $m = -1$, so measuring $L_z$ in $p_x$ gives $\pm\hbar$ with probability $\frac12$ each. The complete set of commuting observables for the real basis is not $\hat H, \hat L^2, \hat L_z$; it includes parity-like operators instead. Project Part 7 asserts (3.19) pointwise against `sph_harm_y`.

### 3.10 The spectrum: the Rydberg formula

A transition $n_i\to n_f$ ($n_i > n_f$) emits a photon of energy $E_{n_i} - E_{n_f}$. (Why the atom radiates at all, given that $|n\ell m\rangle$ is stationary, is the coupling to the quantized field; QM-8 computes the rates.) With $E = hc/\lambda$:

**Derivation.**
1. $\dfrac{hc}{\lambda} = E_{n_i} - E_{n_f} = |E_1|\Big(\dfrac{1}{n_f^2} - \dfrac{1}{n_i^2}\Big)$.
2. Divide by $hc$:
$$
\frac1\lambda = R\Big(\frac{1}{n_f^2} - \frac{1}{n_i^2}\Big), \qquad R_\infty = \frac{|E_1|}{hc} = \frac{m_ee^4}{8\varepsilon_0^2h^3c} = 1.09737\times10^7\ \text{m}^{-1}. \tag{3.20}
$$
3. With the reduced mass, $R_{\rm H} = R_\infty/(1 + m_e/m_p) = 1.09678\times10^7$ m$^{-1}$. ∎

Rydberg had fit this constant from wavelengths. Here it is computed from $m_e$, $e$, $h$, $c$, $\varepsilon_0$. **Series:** $n_f = 1$ is Lyman (ultraviolet, 121.6 nm down to the series limit 91.18 nm), $n_f = 2$ is Balmer (visible, 656.5 nm down to 364.7 nm), $n_f = 3$ is Paschen (infrared, from 1875.6 nm). The wavelengths quoted are vacuum values; tables of visible lines usually give wavelengths in air, shorter by the refractive index $1.00028$ (656.3 nm for H$\alpha$).

**Selection rules (stated; QM-8 derives them).** Electric-dipole transitions, which carry nearly all the light, connect $|n\ell m\rangle$ to $|n'\ell'm'\rangle$ only if $\Delta\ell = \pm1$ and $\Delta m = 0, \pm1$; $\Delta n$ is unrestricted. The reason: the photon carries one unit of angular momentum (QM-5 §3.9), and the dipole operator $\mathbf r$ is odd under parity while $Y_\ell^m$ has parity $(-1)^\ell$ (QM-5 (3.23)). So $2s\to1s$ is forbidden as a one-photon dipole transition: the $2s$ state lives about $0.12$ s (it decays by emitting two photons) instead of the nanoseconds of $2p$.

### 3.11 The real hydrogen atom: the hierarchy of corrections

The Hamiltonian (3.7) leaves things out. Their sizes, in units of the Bohr energy scale $\alpha^2m_ec^2 = 27.2$ eV:

| Effect | Physics | Size | Number |
|---|---|---|---|
| Bohr levels | (3.13) | $\alpha^2m_ec^2$ | $E_1 = -13.6$ eV |
| Reduced mass | CM-5 two-body | $\times m_e/m_p$ | $7.4$ meV shift of $E_1$ |
| Fine structure | relativistic $p^4$ term, spin–orbit, Darwin term | $\alpha^4m_ec^2$ | $2p_{3/2}$–$2p_{1/2}$: $4.53\times10^{-5}$ eV, $10.9$ GHz |
| Lamb shift | QED: vacuum fluctuations of the field | $\sim\alpha^5m_ec^2\ln(1/\alpha)$ | $2s_{1/2}$–$2p_{1/2}$: $1058$ MHz, $4.4\times10^{-6}$ eV |
| Hyperfine | electron spin–proton spin magnetic coupling | $\alpha^4m_ec^2\cdot m_e/m_p$ | $1s$: $1420.4$ MHz, $5.9\times10^{-6}$ eV, $\lambda = 21.1$ cm |

Fine structure (QM-8 derives it) gives $E_{nj} = E_n\big[1 + \frac{\alpha^2}{n^2}\big(\frac{n}{j + 1/2} - \frac34\big)\big]$, which depends on $n$ and the total angular momentum $j = \ell\pm\frac12$ but still not on $\ell$ separately: $2s_{1/2}$ and $2p_{1/2}$ stay degenerate in Dirac's theory. Willis Lamb and Robert Retherford measured them to differ by about 1 GHz in 1947, and explaining it was a founding problem of quantum electrodynamics. The hyperfine line was predicted for interstellar hydrogen by Hendrik van de Hulst in 1944 and detected by Harold Ewen and Edward Purcell in 1951; the 21 cm line maps the hydrogen of the Milky Way. Each correction is smaller than the last by roughly $\alpha\approx1/137$ or $m_e/m_p\approx1/1836$, which is why perturbation theory (QM-8) works so well here.

### 3.12 Why the $\ell$-degeneracy is "accidental": $SO(4)$

Rotational symmetry explains $2\ell+1$. Nothing in rotations explains why $2s$ and $2p$ have the same energy. For the Coulomb potential there is an extra conserved vector, CM-4 §3.8's Laplace–Runge–Lenz vector $\mathbf A = \mathbf p\times\mathbf L - mk\hat{\mathbf r}$ (with $k = e^2/4\pi\varepsilon_0$), which points at perihelion and has length $mk\,e$ (CM-5 §3.12). Its existence is why classical Kepler orbits close and why their energy depends only on $a$.

Pauli's 1926 calculation, in outline (a pointer, not a derivation; Zwiebach ch. 12 gives it in full):
1. The quantum operator $\hat{\mathbf A} = \frac12(\hat{\mathbf p}\times\hat{\mathbf L} - \hat{\mathbf L}\times\hat{\mathbf p}) - mk\hat{\mathbf r}/r$ (symmetrized to be Hermitian) commutes with $\hat H$.
2. On a bound eigenspace of energy $E$, the rescaled $\hat{\mathbf M} = \hat{\mathbf A}/\sqrt{-2mE}$ and $\hat{\mathbf L}$ close under commutators into the Lie algebra of $SO(4)$, rotations in four dimensions. The combinations $\hat{\mathbf J}_{1,2} = \frac12(\hat{\mathbf L}\pm\hat{\mathbf M})$ are two commuting angular momenta.
3. $\hat{\mathbf L}\cdot\hat{\mathbf A} = 0$ forces $j_1 = j_2 \equiv j$, and a Casimir identity gives $E = -mk^2/[2\hbar^2(2j+1)^2]$. With $n = 2j + 1$, that is (3.13), and the degeneracy is $(2j_1+1)(2j_2+1) = n^2$.

So the degeneracy is not accidental; it is a symmetry larger than the obvious one. "Accidental" means only "not explained by rotations." **Any perturbation of $1/r$ breaks it.** Project Part 6 shows this for a screened (Yukawa) potential $-e^{-r/\lambda}/r$, where $2s$ falls below $2p$, and for the 3D oscillator, which has its own different accidental degeneracy (its hidden symmetry is $SU(3)$): $E = (2n_r + \ell + \frac32)\hbar\omega$, with $(N+1)(N+2)/2$ states at $N = 2n_r + \ell$, so $s$ and $d$ are degenerate but $s$ and $p$ are not.

### 3.13 Hydrogen-like atoms

One electron around a nucleus of charge $Ze$ and mass $M$: replace $e^2\to Ze^2$ and $m\to\mu$ in (3.7)–(3.13).
$$
E_n^{(Z)} = Z^2\frac{\mu}{m_e}\frac{E_1}{n^2}, \qquad a^{(Z)} = \frac{m_e}{\mu}\frac{a_0}{Z}, \qquad R_M = \frac{R_\infty}{1 + m_e/M}. \tag{3.21}
$$
$\text{He}^+$: $E_1 = -54.42$ eV, size $a_0/2$, every wavelength divided by $4$ (Example 4.8). $\text{Li}^{2+}$: $E_1 = -9\times13.6 = -122.4$ eV. **Muonic hydrogen** (a muon, $m_\mu = 206.77m_e$, bound to a proton): $\mu = 185.84\,m_e$, $E_1 = -2.53$ keV (X-rays), size $a_0/185.84 = 285$ fm, only about 340 proton radii. The muon spends enough time inside the proton that its levels measure the proton's charge radius; the 2010 measurement at the Paul Scherrer Institute gave $0.84$ fm, smaller than the then-accepted value from ordinary hydrogen, and started the "proton radius puzzle" (Exercise B5 shows why the muon is $\sim\!6\times10^6$ times more sensitive). **Positronium** ($e^+e^-$): $\mu = m_e/2$, so every level is halved (project "Your turn").

### 3.14 Many-electron atoms, qualitatively

This page is qualitative and says so; QM-7 supplies the exclusion principle and QM-8 the variational method that make it quantitative.

**The filling picture.** Treat each electron as moving in an average central potential produced by the nucleus and the other electrons (the "central-field approximation"). Each one-electron state is still labeled $n\ell m_\ell m_s$, because rotational symmetry survives the averaging. The Pauli principle (QM-7) allows one electron per state, so each $(n, \ell)$ subshell holds $2(2\ell+1)$ electrons: s holds 2, p holds 6, d holds 10. Fill from the lowest energy up.

**Screening breaks the $\ell$-degeneracy.** The averaged potential is not $1/r$. Far from the nucleus an outer electron sees charge $+e$ (nucleus $Z$ screened by $Z - 1$ others); close in it sees nearly $Ze$. So $V(r)$ interpolates from $-Ze^2/4\pi\varepsilon_0r$ to $-e^2/4\pi\varepsilon_0r$, and §3.12 says the $\ell$-degeneracy goes. Which way? An s-electron's $R$ is nonzero at the origin and its inner humps (§3.7) penetrate the core, where it sees the large charge; a d-electron is held out by its centrifugal wall. So at fixed $n$, $E_{ns} < E_{np} < E_{nd}$. The effect can beat the $1/n^2$ spacing: for potassium ($Z = 19$) the 19th electron goes into $4s$, not $3d$. The empirical **Madelung rule**, fill in order of increasing $n + \ell$ and then increasing $n$, gives $1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, \ldots$ and the shape of the periodic table: rows of length 2, 8, 8, 18, 18, 32, the $2n^2$ of (3.14) rearranged. Chemistry is the outer subshell; the transition metals are the filling of $3d$ after $4s$. Exercise C3 reproduces the $4s$/$3d$ ordering for potassium with a one-parameter screened potential and the project's solver.

**Hund's rules** (named, not derived; Griffiths & Schroeter §5.2) decide how electrons fill a partly full subshell: maximize total spin $S$, then total orbital $L$, then $J = |L - S|$ for less than half full and $L + S$ for more. They come from the Coulomb repulsion between electrons plus exchange symmetry (QM-7).

### 3.15 Correspondence: Rydberg states and circular orbits

For large $n$ the quantum atom should look like CM-5's Kepler problem. The circular orbit corresponds to $\ell = n - 1$ (maximum angular momentum for the energy, as a circle has maximum $L$ for given $a$), $n_r = 0$.

▶ **Try it first.** For $\ell = n - 1$, the recursion stops at $j_{\max} = 0$. Write $P(r) = |u|^2$ and find its peak.

**Derivation.**
1. $u = c_0\rho^ne^{-\rho}$ with $\rho = r/na_0$, so $P(r) \propto r^{2n}e^{-2r/na_0}$.
2. $\frac{d}{dr}\ln P = \frac{2n}{r} - \frac{2}{na_0} = 0$ at $r_{\rm peak} = n^2a_0$: exactly Bohr's radius $r_n = n^2a_0$ (QM-1 Theorem 3.1).
3. $P(r)$ is a Gamma density in $r$ with shape $2n + 1$ and scale $na_0/2$, so $\langle r\rangle = (2n+1)\frac{na_0}{2} = n^2a_0 + \frac{na_0}{2}$ (matching (3.17) with $\ell = n-1$) and $\sigma_r/\langle r\rangle = 1/\sqrt{2n+1}$. ∎

■ **What just happened.** Bohr's orbit radius is the peak of the circular state's radial density, and the shell sharpens like $1/\sqrt{2n+1}$ as $n$ grows: at $n = 100$ the relative thickness is 7%. Rydberg atoms with $n\sim100$ are micrometers across and are real laboratory objects. Even then, the state $|n, n-1, n-1\rangle$ is a doughnut, not a localized electron going around; a localized orbiting packet is a superposition of many $n$ (a "Rydberg wave packet"), which QM-8's time-dependent tools can build. Also $\langle r\rangle\propto n^2$ for every $\ell$ by (3.17), and $L = \hbar\sqrt{\ell(\ell+1)}\to\hbar n$ as $n\to\infty$, which is Bohr's $L = n\hbar$ recovered only in the limit.

## 4. Worked examples

### Example 4.1: the radial equation in atomic units, and an energy that is not allowed

Atomic units set $\hbar = m_e = e^2/4\pi\varepsilon_0 = 1$, so $a_0 = 1$ and $\hbar^2/m_ea_0^2 = 1$ hartree $= 27.211$ eV. Then (3.3) is
$$-\tfrac12u'' + \Big[-\frac1r + \frac{\ell(\ell+1)}{2r^2}\Big]u = Eu,$$
$\kappa = \sqrt{-2E}$, and $\rho_0 = 2/\kappa$. For $E_1 = -\frac12$ hartree ($-13.606$ eV): $\kappa = 1$, $\rho_0 = 2$, $n = 1$. For $E = -5$ eV $= -0.18375$ hartree: $\kappa = \sqrt{0.36750} = 0.60622$, $\rho_0 = 3.2991$. Theorem 3.1 needs $\rho_0$ to be an even integer, and $3.2991$ is not: no bound state at $-5$ eV. The nearest allowed values are $\rho_0 = 2$ ($-13.6$ eV) and $\rho_0 = 4$ ($-3.40$ eV). With $\rho_0 = 3.2991$ and $\ell = 0$, (3.11) gives $c_1 = (2 - 3.2991)/2\,c_0 = -0.6496\,c_0$, $c_2 = (4 - 3.2991)/6\,c_1 = -0.0759\,c_0$, and every later ratio is positive and tends to $2/(j+1)$: the series never stops and $u$ grows like $e^{\rho}$ (Project Part 3 prints these).

### Example 4.2: $R_{20}$ from the recursion

$n = 2$, $\ell = 0$: $\rho_0 = 4$, $\kappa = 1/2a_0$, $\rho = r/2a_0$.
1. (3.11) with $j = 0$: $c_1 = \frac{2(0 + 0 + 1) - 4}{(1)(2)}c_0 = -c_0$.
2. $j = 1$: $c_2 = \frac{2(1 + 1) - 4}{(2)(3)}c_1 = 0$. The series stops: $v = c_0(1 - \rho)$, $j_{\max} = 1$, $n = 1 + 0 + 1 = 2$ ✓.
3. $u = c_0\rho e^{-\rho}(1 - \rho)$, so $R = u/r = \frac{c_0}{2a_0}\Big(1 - \frac{r}{2a_0}\Big)e^{-r/2a_0}$.
4. Normalize: with $x = r/a_0$, $\int_0^\infty r^2R^2dr = A^2a_0^3\int_0^\infty x^2\big(1 - x + \frac{x^2}{4}\big)e^{-x}dx = A^2a_0^3\big(2 - 6 + \frac{24}{4}\big) = 2A^2a_0^3$, where $A = c_0/2a_0$. So $A = 1/\sqrt{2a_0^3}$:
$$R_{20} = \frac{1}{\sqrt2}a_0^{-3/2}\Big(1 - \frac{r}{2a_0}\Big)e^{-r/2a_0},$$
matching (3.16). One radial node, at $r = 2a_0$ ($n - \ell - 1 = 1$ ✓). Check against (3.15): $L^{(1)}_1(x) = 2 - x$ at $x = r/a_0$ gives $2(1 - r/2a_0)$, and the prefactor $\sqrt{1\cdot1!/(4\cdot2!)}\,a_0^{-3/2} = a_0^{-3/2}/\sqrt8$, so $R_{20} = \frac{2}{\sqrt8}a_0^{-3/2}(1 - r/2a_0)e^{-r/2a_0}$ ✓. The project asserts this and Exercise B1 repeats the procedure for $R_{30}$.

### Example 4.3: the nine states of $n = 3$

$\ell\in\{0, 1, 2\}$ and $n_r = n - \ell - 1$:

| $\ell$ | name | $n_r$ (radial nodes) | $m$ values | count |
|---|---|---|---|---|
| 0 | $3s$ | 2 | $0$ | 1 |
| 1 | $3p$ | 1 | $-1, 0, 1$ | 3 |
| 2 | $3d$ | 0 | $-2, -1, 0, 1, 2$ | 5 |

Total $1 + 3 + 5 = 9 = 3^2$; with spin, 18, the length of the fourth and fifth rows of the periodic table. All nine have $E_3 = -13.606/9 = -1.512$ eV. In the finite-difference solver (Project Part 1) the third $\ell = 0$ eigenvalue, the second $\ell = 1$, and the first $\ell = 2$ all agree with $-1/18 = -0.0555556$ hartree to better than $10^{-6}$ hartree.

### Example 4.4: $\langle r\rangle$, $\langle r^2\rangle$, $\langle1/r\rangle$ in the ground state

From §3.8: $\langle r\rangle = 1.5\,a_0 = 0.0794$ nm, $\langle r^2\rangle = 3a_0^2$, $\sigma_r = 0.866\,a_0$, $\langle1/r\rangle = 1/a_0$, most probable $r = a_0 = 0.0529$ nm. The spread $\sigma_r$ is more than half the mean: there is no "orbit radius" in any sharp sense. Energies: $\langle V\rangle = -27.211$ eV, $\langle T\rangle = +13.606$ eV. Project Part 4 computes all of these by `scipy.integrate.quad` on $u_{10}$ and also on the finite-difference eigenvector, and asserts the virial relation $2\langle T\rangle = -\langle V\rangle$ to $10^{-10}$ on the analytic function.

### Example 4.5: Balmer-$\alpha$ and the Lyman limit

Balmer-$\alpha$ is $3\to2$: $\frac1\lambda = R_{\rm H}\big(\frac14 - \frac19\big) = \frac{5R_{\rm H}}{36}$, so $\lambda = \frac{36}{5\times1.096776\times10^7}$ m $= 656.47$ nm (vacuum; $656.29$ nm in air). Red: the color of an emission nebula. With $R_\infty$ instead you get $656.11$ nm, off by $0.36$ nm, much more than spectroscopic precision, which is why the reduced mass matters here. Photon energy $hc/\lambda = 13.598\,(1/4 - 1/9) = 1.889$ eV.

The Lyman limit is $\infty\to1$: $\lambda = 1/R_{\rm H} = 91.18$ nm, photon energy $13.598$ eV. Light at shorter wavelength ionizes ground-state hydrogen; that is why interstellar hydrogen makes the universe nearly opaque below 91 nm, and why the Lyman limit is a sharp edge in the spectra of distant galaxies. Lyman-$\alpha$ ($2\to1$): $\lambda = 4/(3R_{\rm H}) = 121.57$ nm.

### Example 4.6: the probability that a 1s electron is inside the nucleus

Take the nucleus as a sphere of radius $b = 10^{-15}$ m. $P(r < b) = \int_0^b\frac{4}{a^3}r^2e^{-2r/a}dr$. With $s = 2r/a$:
$$P(r<b) = \frac12\int_0^{2x}s^2e^{-s}ds = 1 - e^{-2x}(1 + 2x + 2x^2), \qquad x = \frac ba.$$
Here $x = 10^{-15}/5.2918\times10^{-11} = 1.8897\times10^{-5}$. Do **not** evaluate the closed form in floating point: it is $1$ minus a number within $10^{-14}$ of $1$, and double precision keeps only about two digits of the difference. Expand instead: $e^{-2x}(1 + 2x + 2x^2) = 1 - \frac43x^3 + 2x^4 - \ldots$, so
$$P \approx \tfrac43x^3\big(1 - \tfrac32x\big) = \tfrac43(1.8897\times10^{-5})^3 = 9.00\times10^{-15}.$$
(Small $x$ means $e^{-2r/a}\approx1$ inside, so $P\approx|\psi(0)|^2\cdot\frac43\pi b^3$, the same thing.) The project computes it three ways: the series, `quad`, and the naive closed form, which comes out as $9.10\times10^{-15}$, wrong in the second digit.

A probability of $10^{-14}$ sounds like nothing. It is not zero, and that is the point: a classical orbit at $a_0$ never comes near the nucleus, while the real $1s$ electron has its *largest* density there. This overlap is what makes electron capture possible (a proton-rich nucleus absorbs an inner s-electron) and what produces the hyperfine 21 cm line, whose size is set by $|\psi(0)|^2$.

### Example 4.7: the circular state $n = 10$ against Bohr's orbit

$n = 10$, $\ell = 9$, $n_r = 0$: $P(r)\propto r^{20}e^{-r/5a_0}$. By §3.15 the peak is at $n^2a_0 = 100a_0 = 5.29$ nm, which is Bohr's $r_{10} = 100a_0$ exactly (CM-5's circular orbit (3.9) with $L = 10\hbar$ gives the same $r_0 = L^2/\mu k = 100a_0$). Mean $\langle r\rangle = \frac{a_0}{2}(300 - 90) = 105a_0$, $\sigma_r/\langle r\rangle = 1/\sqrt{21} = 0.218$. Compare $1s$: $\sigma_r/\langle r\rangle = 0.866/1.5 = 0.577$. Project Part 5 solves the $\ell = 9$ radial equation numerically on a grid to $r = 400a_0$ and finds the peak of $|u|^2$ at $100.0a_0$ and $E = -0.005000$ hartree $= -1/(2\cdot10^2)$.

### Example 4.8: Pickering's lines are ionized helium

In 1896 Edward Pickering found a series of lines in the star $\zeta$ Puppis that interleaved with the Balmer lines. Bohr pointed out in 1913 that they come from $\text{He}^+$. By (3.21) with $Z = 2$, $R_{\rm He} = R_\infty/(1 + m_e/m_\alpha) = 1.09722\times10^7$ m$^{-1}$ and $\frac1\lambda = 4R_{\rm He}\big(\frac{1}{n_f^2} - \frac{1}{n_i^2}\big)$. For $n_f = 4$, $4R(1/16 - 1/n_i^2) = R(1/4 - 4/n_i^2) = R\big(\frac{1}{2^2} - \frac{1}{(n_i/2)^2}\big)$: even $n_i$ reproduce (nearly) the Balmer lines, odd $n_i$ fall between them. The $\text{He}^+$ $4\to3$ line: $\frac1\lambda = 4R_{\rm He}\cdot\frac{7}{144}$, $\lambda = 468.72$ nm (vacuum), the strong He II line used in stellar classification. Ground state: $E_1 = 4\times(-13.606)/(1 + m_e/m_\alpha) = -54.415$ eV.

## 5. Why it matters and where it connects

- **QM-8.** Every perturbation in the table of §3.11 is a matrix element between the $|n\ell m\rangle$ of this module: fine structure needs $\langle1/r^3\rangle$ and $\langle1/r^2\rangle$, the Stark effect mixes the degenerate $2s$ and $2p_0$ (possible only because of §3.12's degeneracy), and the selection rules of §3.10 come from $\langle n'\ell'm'|\mathbf r|n\ell m\rangle$. The variational method's first test case is helium with hydrogenic trial functions.
- **QM-7.** Pauli plus §3.14's filling gives the periodic table; the helium ground state is two $1s$ electrons in a spin singlet.
- **CM-4, CM-5.** The two-body reduction, the effective potential, Bohr's circular orbit as the peak of $|u_{n,n-1}|^2$ (§3.15), the virial theorem, and the LRL vector as the source of the degeneracy (§3.12).
- **EM-3, QM-3.** The Frobenius-termination argument a third time; the finite-difference solver applied to the half-line.
- **SM-5.** Hydrogen's levels and their $2n^2$ degeneracies enter the Saha equation for the ionization of stellar atmospheres; the $n^2$ growth makes the partition function formally divergent, a famous subtlety resolved by the finite size of real atoms in a plasma.
- **Astronomy and chemistry.** The Lyman limit, the 21 cm line, He II 468.6 nm; the real orbitals of §3.9 are the starting point for molecular orbitals and hybridization ($sp^3$ is a change of basis among $2s$ and the three $2p$'s).

## 6. Exercises

See `exercises.md`. A (6): degeneracy counting, normalizing $R_{21}$, series wavelengths, $\text{He}^+$, $|\psi(0)|^2$, real d-orbitals. B (7): the $v$-equation and $R_{30}$ by recursion, $\langle r\rangle$ and $\langle r^2\rangle$ for $2p$, the virial theorem for power laws, the spherical well, muonic hydrogen inside the proton, the deuterium isotope shift, Rydberg-state widths. C (3): the LRL degeneracy count, $\langle1/r\rangle$ and $\langle1/r^2\rangle$ by Feynman–Hellmann, the $4s$/$3d$ ordering in potassium by code. B1, B2, B4, and B6 are essential; C3 extends `project.py`.

## 7. Coding project

`project.py` (about 10 s). Part 1: `radial_solve(V, l, h, rmax, k)`, QM-3's tridiagonal finite-difference solver applied to (3.3) with $u(0) = u(r_{\max}) = 0$; hydrogen energies for every $n\le4$ and $\ell < n$, asserted to $-1/2n^2$ within $10^{-3}$ relative (actually about $3\times10^{-5}$ at $h = 0.01$), and the $\ell$-degeneracy asserted. Part 2: the analytic $R_{n\ell}$ (3.15) from `genlaguerre`, the old-convention polynomial built in SymPy and the $q!$ factor asserted, (3.16) checked, numerical vs analytic radial functions, node counts $n - \ell - 1$, and `radial_densities.png`. Part 3: SymPy derivation of (3.10) and (3.11); the series at $\rho_0 = 4$ (terminates, gives $R_{20}$) and $\rho_0 = 4.01$ (tracks $e^{2\rho}$), Example 4.1's coefficients, `series_blowup.png`. Part 4: ground-state expectation values and virial (Example 4.4), (3.17) for all $n\le4$, the inside-the-nucleus probability three ways (Example 4.6). Part 5: the $n = 10$ circular state (Example 4.7). Part 6: non-Coulomb potentials (Yukawa splitting, 3D oscillator degeneracies, spherical well zeros) and `degeneracy_lifted.png`. Part 7: real orbitals (3.19) checked pointwise, cross-sections of $1s$, $2p_0$, $3d_0$, $2p_x$, $3d_{xy}$ in `orbitals.png`. Part 8: SI constants, $R_\infty$, the Lyman/Balmer/Paschen wavelengths (Example 4.5), `spectrum.png`, $\text{He}^+$ from the solver with $Z = 2$ (Example 4.8), muonic hydrogen, deuterium, the §3.11 scales, and the potassium screening model of Exercise C3.

**Build it yourself first:** write `radial_solve` from QM-3's `fd_eigs` by adding the centrifugal term and starting the grid at $r = h$, and confirm $-1/2n^2$ before reading Part 1; write the recursion (3.11) as a loop and find by experiment what happens when $\rho_0$ is not an even integer. Then the "Your turn" items: the infinite spherical well from the solver against `spherical_jn` zeros, positronium, and a $Z_{\rm eff}$ model of the $4s$/$3d$ ordering.

## 8. Where you'll struggle

- **The solution is long, and most of it is three moves you have already made.** People drown in the algebra of §3.5 because they read it as new. It is not: asymptotics at both ends, factor them out, power series, termination. Write the four headings on a card, and when you redo the derivation (you should, closed-book, twice), check off each move. The only hydrogen-specific step is that $\rho_0$, a disguised energy, sits in the recursion's numerator.
- **$\ell(\ell+1)$ vs $\ell^2$.** You will write $\hbar^2\ell^2/2mr^2$ for the centrifugal term at least once, from CM-5 habit. The fix is to remember where the term comes from: $\hat L^2/2mr^2$, and $\hat L^2$'s eigenvalue is $\hbar^2\ell(\ell+1)$. The difference is physical, not cosmetic: $\langle L_x^2 + L_y^2\rangle = \hbar^2\ell$ even when $m = \ell$. The $\ell = 0$ case is the opposite surprise: no centrifugal barrier and no angular kinetic energy, which a classical orbit cannot have.
- **Laguerre conventions differ by factorials.** Griffiths' 2nd edition, Griffiths & Schroeter's 3rd, Arfken, and SciPy do not all mean the same polynomial by "associated Laguerre." The prefactor in (3.15) is only right with the modern polynomial. The fix is mechanical: before using anyone's $R_{n\ell}$ formula, evaluate it for $n = 2$, $\ell = 0$ at $r = 0$ and compare with $R_{20}(0) = a_0^{-3/2}/\sqrt2$. Project Part 2 maps the two conventions with an assert.
- **The orbit picture.** Bohr's model is in your head and it will keep producing wrong intuitions: that the ground state has angular momentum, that the electron has a radius, that it avoids the nucleus. The $1s$ density is largest *at* the nucleus, $\sigma_r$ is 58% of $\langle r\rangle$, and Example 4.6's $9\times10^{-15}$ is a nonzero probability to be inside the proton. Only the circular Rydberg states (§3.15) resemble orbits, and only in $r$, not in angle.
- **Units.** SI hydrogen formulas carry $4\pi\varepsilon_0$ in places that are easy to lose. Do every derivation in atomic units, then restore SI at the end by $E\to E\times27.211$ eV and $r\to r\times0.0529$ nm. Check every SI result against $13.6$ eV and $0.0529$ nm.
- **Catastrophic cancellation.** Example 4.6's closed form gives the wrong answer in floating point. Whenever you compute "1 minus something close to 1," expand analytically first, or use `np.expm1`.
- **The SO(4) explanation is a pointer, not a derivation.** §3.12 tells you what is true and where the proof lives. If you want it, Zwiebach ch. 12 takes about six pages; do it after QM-8, not now.
- **The periodic-table page is qualitative.** §3.14 is an honest sketch; the Madelung rule has exceptions (chromium and copper fill $3d$ early), and the real theory (Hartree–Fock) is beyond this curriculum. Do not try to derive the rule; know why screening favors s over d.
- **The finite-difference solver lies near $r = 0$ for s-states if the grid is coarse.** The Coulomb singularity makes the $1s$ error the largest. Halve $h$ and watch the error drop by 4 before trusting a digit (Project Part 1 prints the $1s$ error at two step sizes).

## 9. Check yourself

1. Derive (3.3) from (3.2) by substituting $R = u/r$.
2. Why must $u(0) = 0$ for $\ell = 0$, given that $R\sim1/r$ would be normalizable?
3. What are the asymptotic forms of $u$ at $\rho\to0$ and $\rho\to\infty$, and why is $e^{+\rho}$ rejected?
4. State the recursion (3.11) and explain in two sentences why a non-terminating series is unacceptable.
5. Derive $E_n$ from $\rho_0 = 2n$ and $\rho_0 = 2/\kappa a_0$.
6. List the states with $n = 2$ and count them, with and without spin.
7. Compute $\langle r\rangle$ for the ground state and the most probable $r$. Why are they different?
8. State the virial theorem for Coulomb and use it to get $\langle T\rangle$ and $\langle V\rangle$ in any state from $E_n$.
9. Write $p_x$ in terms of $Y_1^{\pm1}$. Is it an eigenstate of $\hat L_z$?
10. Compute the Balmer-$\alpha$ wavelength from $R$. Why does the reduced mass matter here?
11. Which central potentials have energies independent of $\ell$ at fixed $n$, and what conserved quantity is responsible?
12. Why does the 19th electron of potassium go into $4s$ rather than $3d$?

## 10. Pair with

- **Griffiths & Schroeter (3rd ed.), §4.1–4.2.** The standard derivation in the same order as §3.1–3.8, with the modern Laguerre convention; its end-of-chapter problems do the normalization and Kramers' relation that §3.7 quotes.
- **MIT 8.04, lectures 20–24.** Central potentials, the effective potential, and hydrogen on video; good for seeing the radial wells drawn.
- **Zwiebach, *Mastering Quantum Mechanics*, ch. 11–12.** The cleanest account of the hydrogen-atom algebra and the full $SO(4)$/LRL derivation that §3.12 only points to.
- **Shankar, *Principles of Quantum Mechanics*, ch. 13.** The hydrogen atom with the most careful discussion of numerical estimates and of the real atom's corrections.
- **Feynman Lectures on Physics, vol. III, ch. 19.** The periodic table from hydrogen orbitals, qualitatively, in the spirit of §3.14; read it for the physical picture.
