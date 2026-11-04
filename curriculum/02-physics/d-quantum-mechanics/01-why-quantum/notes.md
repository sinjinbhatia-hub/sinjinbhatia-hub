# QM-1: Why Quantum: The Experiments That Broke Classical Physics; the Postulates Previewed as a Vector-Space Theory

**Time:** ~6 hours · **Prerequisites:** LA-8 (spectral theorem, Hermitian operators, Dirac-notation preview); EM-8 and SM-5 helpful · **Unlocks:** QM-2

> **Conventions in this module.** SI units, with the electron-volt ($1\ \text{eV} = 1.602\times10^{-19}$ J) as the practical energy unit; $hc = 1239.84$ eV nm. $h$ is Planck's constant and $\hbar = h/2\pi$ is written explicitly everywhere (except in `project.py` Parts 3–6, which set $\hbar = 1$). $f$ is frequency in Hz, $\omega = 2\pi f$, $\lambda$ wavelength, $k = 2\pi/\lambda$. $e > 0$ is the elementary charge, $m$ the electron mass, $k_B$ Boltzmann's constant. **Dirac notation (LA-8 §3.12):** a ket $|\psi\rangle$ is a column vector, the bra $\langle\phi|$ is its conjugate transpose, and $\langle\phi|\psi\rangle = \boldsymbol\phi^\dagger\boldsymbol\psi$ is **conjugate-linear in the first slot**, linear in the second (LA-8's physics convention). Observables are Hermitian operators, written with a hat when the distinction matters ($\hat A$). Spin basis: $|{\uparrow}\rangle = |{+z}\rangle = (1, 0)^\top$, $|{\downarrow}\rangle = |{-z}\rangle = (0, 1)^\top$, which are LA-8's $\mathbf e_1, \mathbf e_2$ and QC-1's $|0\rangle, |1\rangle$. Pauli matrices as in LA-8. A unit vector $\hat{\mathbf n} = (\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$. $[A, B] = AB - BA$ is the commutator. $\Delta A = \sqrt{\langle A^2\rangle - \langle A\rangle^2}$ is the standard deviation.

## 0. Where this sits

You need LA-8: complex inner products, Hermitian matrices have real eigenvalues and orthonormal eigenbases (Theorems 3.9, 3.10, 3.12), Parseval (Theorem 3.6), projectors $|q\rangle\langle q|$ (Theorem 3.13), unitary matrices (Theorem 3.20), the Pauli matrices and the eigenvectors of $\hat{\mathbf n}\cdot\boldsymbol\sigma$ (Proposition 3.22), and the Born-rule preview (3.16). LA-7's matrix exponential (§3.9) gives time evolution. EM-8's Example 4.4 (the classical atom collapses in $10^{-11}$ s) and the energy–momentum relation EM-8 (3.25) are used once each. SM-5 (later in the roadmap) derives Planck's law from Bose–Einstein statistics; here we derive it the way Planck could. This module makes QM-2 possible: QM-2 runs the same postulates on wavefunctions.

## 1. The problem that forced this idea

Between 1900 and 1927 physics was pushed, one experiment at a time, into a theory nobody wanted. Each step below was forced: a measurement that classical physics predicted confidently and got wrong. Read this section as a sequence of moves in a game where the board kept changing. The formulas are derived in §3; here are the facts.

**1.1 Blackbody radiation (1900).** A hot cavity with a small hole emits a spectrum that depends only on temperature. Classical physics treats each standing-wave mode of the cavity as an oscillator with average energy $k_BT$ (equipartition). The number of modes per unit wavelength grows as $\lambda^{-4}$, so the predicted energy density, the Rayleigh–Jeans law, diverges at short wavelength and the total energy is infinite. Paul Ehrenfest named this the "ultraviolet catastrophe" in 1911. Measurements by Lummer, Pringsheim, Rubens and Kurlbaum in Berlin showed a peaked curve. Wilhelm Wien had a formula (1896) that fit short wavelengths but failed at long ones. In October 1900 Max Planck found an interpolating formula that fit everything; in December he derived it by assuming that an oscillator of frequency $f$ can only hold energy in multiples of $hf$. He regarded this as a calculational device. At 500 nm and 5000 K the classical prediction is 55 times too large (Example 4.1).

**1.2 The photoelectric effect (1905).** Light falling on a metal ejects electrons. Classically, a light wave's energy is spread across its wavefront and grows with intensity, so brighter light should give faster electrons, and very dim light should need time to accumulate enough energy. Philipp Lenard (1902) found the opposite: the electrons' maximum energy depends on the light's frequency, not its intensity, and there is a threshold frequency below which nothing comes out however bright the light. Albert Einstein (1905) proposed that light itself comes in quanta of energy $hf$, each absorbed by one electron. That predicts a straight line, stopping potential versus frequency, with slope $h/e$ for every metal. Robert Millikan, who expected to refute it, measured the line in 1914–1916 and confirmed it, with $h$ to within about half a percent. Example 4.2 does sodium.

**1.3 The stability of atoms and the spectral lines (1885–1913).** Rutherford's 1911 scattering experiments showed that an atom is a tiny positive nucleus with electrons around it. An orbiting electron accelerates, and by Larmor's formula it radiates. EM-8 Example 4.4 integrated the energy loss: a classical hydrogen atom spirals in from $0.053$ nm in $1.56\times10^{-11}$ s, emitting a continuous chirp (`project.py` recomputes it). Real atoms last forever and emit sharp lines at fixed wavelengths. Johann Balmer (1885) found that hydrogen's visible lines fit $\lambda = B\,n^2/(n^2 - 4)$; Johannes Rydberg (1888–1890) wrote the general form $1/\lambda = R(1/n_f^2 - 1/n_i^2)$. Nobody could say why. Niels Bohr (1913) postulated that only orbits with angular momentum $L = n\hbar$ are allowed, that electrons in them do not radiate, and that a jump between orbits emits one photon carrying the energy difference. That gives Rydberg's formula with $R$ computed from $m$, $e$, $h$ (§3.3, Example 4.3). The postulates were contradictory with classical electrodynamics and Bohr knew it. The numbers were right.

**1.4 The Stern–Gerlach experiment (1922).** Otto Stern and Walther Gerlach, in Frankfurt, sent a beam of silver atoms through a magnet whose field was strongly non-uniform. A magnetic moment in a field gradient feels a force proportional to its component along the field (§3.4). Classically the moments point every which way, so the beam should smear into a continuous band. It split into two distinct spots. They interpreted this, wrongly, as Bohr–Sommerfeld "space quantization" of orbital angular momentum; the right explanation, an intrinsic angular momentum of the electron with only two possible components, came from George Uhlenbeck and Samuel Goudsmit in 1925 and was named spin. Silver's single outer electron has zero orbital angular momentum, so the experiment was measuring the spin of one electron all along. The decisive pattern, what happens when you put analyzers in sequence along different axes, is the cleanest entry to the whole theory (Townsend's book starts there, and so does §2).

**1.5 The Compton effect (1923).** Arthur Compton scattered molybdenum X-rays (wavelength 71 pm) off graphite and found that the scattered rays had a longer wavelength, by an amount that depended only on the scattering angle. Classically, a wave scattering off a free electron makes the electron oscillate at the incident frequency and re-radiate at that same frequency: no shift. Compton treated the X-ray as a particle with energy $hf$ and momentum $h/\lambda$ colliding with an electron, and the shift came out exactly (§3.5, Example 4.4). After Compton, the photon's momentum was not optional.

**1.6 Matter waves (1924–1927).** Louis de Broglie, in his 1924 thesis, turned Einstein's relation around: if light waves carry momentum $p = h/\lambda$, then particles with momentum $p$ should have a wavelength $\lambda = h/p$. Bohr's rule then becomes a standing-wave condition: a whole number of wavelengths fits around the orbit (§3.6). In 1927 Clinton Davisson and Lester Germer at Bell Labs, scattering 54 eV electrons off a nickel crystal, saw a diffraction peak at $50°$, where the crystal's atomic spacing and de Broglie's wavelength agree to about one percent (Example 4.5). George Paget Thomson saw electron diffraction rings through thin films in the same year.

**1.7 The double slit, one electron at a time.** Feynman (*Lectures* III, ch. 1) calls the two-slit experiment the phenomenon that "has in it the heart of quantum mechanics" and its only mystery. Fire electrons at a wall with two slits and record where each lands on a screen. Each electron arrives at one point, like a particle. After many arrivals, the points build up an interference pattern, like a wave: there are places on the screen where electrons arrive with only slit 1 open, and with only slit 2 open, but almost never with both open. The experiment was done with genuine single electrons by Merli, Missiroli and Pozzi (Bologna, 1974) and by Tonomura and co-workers (Hitachi, 1989), whose film of the pattern building up dot by dot is worth watching. And if you add any apparatus that records which slit each electron went through, however gently, the interference disappears. Feynman's three rules summarize the data: (i) an event's probability is the squared magnitude of a complex amplitude; (ii) when an event can happen in several indistinguishable ways, the amplitudes add; (iii) when the ways can be distinguished, even in principle, the probabilities add.

**1.8 Two theories, then one (1925–1927).** In 1925 Werner Heisenberg wrote down rules for arrays of transition amplitudes; Max Born and Pascual Jordan recognized the arrays as matrices and, with Heisenberg, built **matrix mechanics**, in which position and momentum are matrices with $xp - px = i\hbar$. Wolfgang Pauli derived the hydrogen spectrum from it in early 1926. Also in 1926 Erwin Schrödinger, starting from de Broglie, wrote a wave equation, solved it for hydrogen and got Bohr's levels: **wave mechanics**. The two looked nothing alike. Schrödinger showed later in 1926 that they give the same predictions (Pauli and Carl Eckart independently showed the same); Dirac's and Jordan's transformation theory (1926–27) and von Neumann's Hilbert-space book (1932) explained why: both are representations of one structure, vectors in a complex inner product space and Hermitian operators on it. That structure is what LA-8 taught you, and it is the subject of the rest of this module.

## 2. Intuition

Take the Stern–Gerlach data first, because they are cleaner than anything else in §1. Call an analyzer that splits the beam along axis $\hat{\mathbf n}$ "SG$\hat{\mathbf n}$." Each atom leaves through one of two ports, $+$ or $-$. Block one port and you have a filter.

```
 oven ─► SGz ─(+)─► SGz ─► 100% (+), 0% (−)          experiment 1
 oven ─► SGz ─(+)─► SGx ─►  50% (+), 50% (−)          experiment 2
 oven ─► SGz ─(+)─► SGx ─(+)─► SGz ─► 50% (+), 50% (−)    experiment 3
```

Experiment 1: a measurement repeated gives the same answer. Experiment 2: an atom with a definite $z$-answer has a perfectly random $x$-answer. Experiment 3 is the one that breaks the classical picture: after the $x$-measurement, the $z$-answer is random again. The $x$-filter did not merely reveal something; it erased the $z$-information. And if you build an "SGx" that splits the beam and then recombines both halves without recording which way each atom went, the atoms come out still all $+z$ (Feynman III ch. 5 builds this apparatus; `project.py` Part 3 simulates it).

Now look for a mathematical object that does this. A classical spinning top has a definite vector $\mathbf S$, and measuring $S_x$ then $S_z$ would return the same $S_z$ as before. That fails experiment 3. Try instead: **the state is a unit vector in a two-dimensional space**. The two $z$-outcomes are two orthonormal vectors $|{\uparrow}\rangle$, $|{\downarrow}\rangle$. The two $x$-outcomes are another orthonormal basis $|{+x}\rangle$, $|{-x}\rangle$ of the same space, at "45 degrees" to the first, so that each $z$-vector has equal-sized components along both $x$-vectors. Rule: the probability of an outcome is the squared length of the state's component along that outcome's vector, and after the measurement the state *is* that vector. Run the three experiments in your head with this rule. Experiment 1: $|{\uparrow}\rangle$ has component 1 along itself. Experiment 2: $|{\uparrow}\rangle$ has components of size $1/\sqrt2$ along $|{\pm x}\rangle$, so probabilities $\frac12$. Experiment 3: the state is now $|{+x}\rangle$, whose components along $|{\uparrow}\rangle$ and $|{\downarrow}\rangle$ are each $1/\sqrt2$. The rule reproduces everything.

You have seen this picture before. A polarizer is a measurement with two outcomes (EM-7 §3.6): unpolarized light through a vertical polarizer then a horizontal one gives nothing, but insert a $45°$ polarizer between them and an eighth of the light gets through. The middle polarizer erased the "vertical" information exactly as SGx erased the $z$-information. The Jones vector is a 2-vector; the Stern–Gerlach state is a 2-vector; the mathematics is the same, with one twist (§3.9): spin "orthogonal" axes are $180°$ apart in real space, polarizer ones are $90°$ apart.

Why must the space be *complex*? Because there is a third axis. The $y$-outcomes must also be "45 degrees" from both the $z$-basis and the $x$-basis. In $\mathbb R^2$, a basis at $45°$ to both $\{\mathbf e_1, \mathbf e_2\}$ and $\{(1, 1)/\sqrt2, (1, -1)/\sqrt2\}$ does not exist. In $\mathbb C^2$ it does: $(1, \pm i)/\sqrt2$, the eigenvectors of $\sigma_y$ in LA-8 Proposition 3.22.

The double slit says the same thing in another way. Each route to a point on the screen contributes a complex number; the numbers add before squaring, so they can cancel. Vectors add; probabilities are squared lengths; squared lengths of sums are not sums of squared lengths. That last fact *is* interference, and it is Pythagoras failing for non-orthogonal vectors. If the electron leaves a record of its path, the two routes end in orthogonal states of the electron-plus-record, Pythagoras holds, and the interference term is gone.

Collect the pattern from §1. In every experiment:

1. **States add like vectors** (amplitudes add: double slit, Stern–Gerlach recombination).
2. **Measurement outcomes are discrete** (two spots; spectral lines; quanta $hf$).
3. **Probabilities are squared lengths of components** (Born).
4. **Measuring changes the state** to the vector of the outcome (experiment 3).

Guess the theorem. If outcomes are labeled by an orthonormal basis and probabilities are squared components, the probabilities sum to 1 by Parseval. If the outcomes are numbers attached to basis vectors, then "the quantity being measured" is a matrix with those numbers as eigenvalues and those vectors as eigenvectors. Real numbers and an orthonormal eigenbasis: a Hermitian matrix, by the spectral theorem. That is postulates 1–4 of §3.7. You could have written them down from LA-8 and these four facts.

## 3. Formal development

Sections 3.1–3.6 derive the formula each experiment of §1 hinges on. Sections 3.7–3.12 state the postulates and work them out completely for spin-1/2.

### 3.1 Blackbody radiation: Rayleigh–Jeans and Planck

The energy density per unit wavelength in a cavity at temperature $T$ is $u_\lambda = (\text{modes per unit volume per unit wavelength})\times(\text{average energy per mode})$.

▶ **Try it first.** Count the standing electromagnetic waves in a cube of side $\ell$ with wavelength between $\lambda$ and $\lambda + d\lambda$. The allowed wavevectors are $\mathbf k = (\pi/\ell)(n_x, n_y, n_z)$ with positive integers $n_i$, and each $\mathbf k$ has two polarizations.

**Mode count.**
1. The number of allowed $\mathbf k$ with $|\mathbf k| < k$ is the volume of one octant of a ball of radius $k$, divided by the volume $(\pi/\ell)^3$ per lattice point: $\frac18\cdot\frac43\pi k^3\cdot\ell^3/\pi^3$. With two polarizations, $N(k) = k^3V/(3\pi^2)$, $V = \ell^3$.
2. $dN/dk = k^2V/\pi^2$. With $k = 2\pi/\lambda$, $|dk| = (2\pi/\lambda^2)\,d\lambda$, so $dN = \frac{4\pi^2}{\lambda^2}\frac{V}{\pi^2}\frac{2\pi}{\lambda^2}d\lambda = \frac{8\pi V}{\lambda^4}d\lambda$.

**Classical average energy.** Each mode is a harmonic oscillator; equipartition (SM-4) gives $\langle E\rangle = k_BT$. Hence the **Rayleigh–Jeans law**
$$u_\lambda^{\rm RJ} = \frac{8\pi k_BT}{\lambda^4}. \tag{3.1}$$
$\int_0^\infty u^{\rm RJ}_\lambda d\lambda$ diverges at $\lambda\to0$. Nothing in the derivation is shaky; the classical theory really predicts this.

▶ **Try it first.** Now let a mode of frequency $f$ hold only energies $E_n = nhf$, $n = 0, 1, 2, \dots$, with Boltzmann probabilities $\propto e^{-E_n/k_BT}$. Compute the average energy.

**Planck's average energy.** Write $\beta = 1/k_BT$ and $\epsilon = hf$.
1. The partition function is a geometric series: $Z = \sum_{n\ge0}e^{-n\beta\epsilon} = 1/(1 - e^{-\beta\epsilon})$.
2. $\langle E\rangle = \sum_n E_ne^{-\beta E_n}/Z = -\partial_\beta\ln Z = \frac{\epsilon e^{-\beta\epsilon}}{1 - e^{-\beta\epsilon}} = \frac{hf}{e^{hf/k_BT} - 1}$.
3. With $f = c/\lambda$ and the mode count:
$$u_\lambda = \frac{8\pi hc}{\lambda^5}\frac{1}{e^{x} - 1},\qquad x = \frac{hc}{\lambda k_BT}. \tag{3.2}$$
4. Limits: for $x\ll1$, $e^x - 1\approx x$ and (3.2) becomes (3.1). For $x\gg1$, (3.2) becomes Wien's $8\pi hc\lambda^{-5}e^{-x}$. The ratio of classical to quantum is exactly $(e^x - 1)/x$.

■ **What just happened.** Discreteness suppresses high-frequency modes: when $hf\gg k_BT$, the smallest allowed packet of energy is so expensive that the mode is almost never excited. Equipartition is right for modes with $hf\ll k_BT$ and wrong for the rest. The total energy is now finite, $\propto T^4$ (Stefan–Boltzmann; `project.py` checks the integral), and the peak sits at $\lambda_{\max}T = 2.898$ mm K (Exercise C1). SM-5 rederives (3.2) by counting photons as bosons.

### 3.2 The photoelectric law

A photon of energy $hf$ is absorbed by one electron. Removing the electron from the metal costs at least the **work function** $W$, so the fastest electrons leave with kinetic energy $hf - W$. A reverse voltage $V_s$ stops them when $eV_s$ equals that kinetic energy:
$$eV_s = hf - W. \tag{3.3}$$
Three predictions, all confirmed by Millikan: $V_s$ is linear in $f$; the slope is $h/e$, the same for every metal; there is a threshold $f_0 = W/h$ below which no electrons emerge, whatever the intensity. Intensity changes the *number* of photons and so the current, not $V_s$.

### 3.3 The Bohr model

▶ **Try it first.** An electron of mass $m$ circles a proton at radius $r$ with speed $v$. Write Newton's second law with the Coulomb force, then impose $L = mvr = n\hbar$. Solve for $r$ and the total energy.

**Theorem 3.1 (Bohr levels).** With $L = n\hbar$, $n = 1, 2, \dots$ and an infinitely heavy nucleus,
$$r_n = n^2a_0,\quad a_0 = \frac{4\pi\varepsilon_0\hbar^2}{me^2};\qquad E_n = -\frac{E_R}{n^2},\quad E_R = \frac{me^4}{2(4\pi\varepsilon_0)^2\hbar^2}. \tag{3.4}$$
A transition $n_i\to n_f$ emits a photon with
$$\frac1\lambda = R\Big(\frac1{n_f^2} - \frac1{n_i^2}\Big),\qquad R = \frac{E_R}{hc}. \tag{3.5}$$

*Derivation.* Let $\kappa = e^2/4\pi\varepsilon_0$.
1. Circular orbit: $mv^2/r = \kappa/r^2$, so $mv^2 = \kappa/r$.
2. Quantization: $v = n\hbar/(mr)$. Substitute: $m\,n^2\hbar^2/(m^2r^2) = \kappa/r$, so $r = n^2\hbar^2/(m\kappa) = n^2a_0$.
3. Energy: $E = \frac12mv^2 - \kappa/r = \frac12\kappa/r - \kappa/r = -\kappa/(2r)$ (the virial theorem for $1/r$, CM-4). With step 2, $E_n = -m\kappa^2/(2\hbar^2n^2)$, which is (3.4).
4. Photon energy $hc/\lambda = E_{n_i} - E_{n_f} = E_R(1/n_f^2 - 1/n_i^2)$, which is (3.5). ∎

The speed in the ground state is $v_1 = \kappa/\hbar = \alpha c$ with $\alpha = e^2/(4\pi\varepsilon_0\hbar c)\approx1/137$, the fine-structure constant: nonrelativistic, barely. A finite nuclear mass $M$ replaces $m$ by the reduced mass $\mu = mM/(m + M)$ (CM-5's two-body reduction), so $R_M = R_\infty/(1 + m/M)$ (Exercise B2).

■ **What just happened.** One rule, $L = n\hbar$, turned Rydberg's empirical constant into $me^4/(8\varepsilon_0^2h^3c)$, correct to the precision of 1913 spectroscopy. The model asserts, without reason, that allowed orbits do not radiate.

**The Bohr–Sommerfeld rule and why it is not the whole story.** Arnold Sommerfeld (1916) and William Wilson generalized Bohr's rule to $\oint p\,dq = nh$ for each periodic coordinate, where $\oint p\,dq$ is CM-8's action variable $J$ (CM-8 (3.24)). For the circular orbit, $q = \varphi$, $p = L$ constant, $\oint L\,d\varphi = 2\pi L = nh$, which is $L = n\hbar$. For a harmonic oscillator, $J = 2\pi E/\omega$ (`project.py` checks the integral), so $E_n = n\hbar\omega$. CM-8 §3.9 explains why action variables were the natural things to quantize (they are adiabatic invariants). And the rule fails, in ways that told people a new mechanics was needed: the oscillator's true levels are $(n + \frac12)\hbar\omega$ (QM-3); hydrogen's true ground state has zero angular momentum, which no Bohr orbit has (QM-6); helium, with two electrons, could not be made to work at all; the rule says nothing about how bright each line is; and it needs periodic orbits, which chaotic systems do not have (Einstein pointed this out in 1917). QM-8's WKB method recovers $\oint p\,dq = (n + \frac12)h$ as an approximation to the real theory.

### 3.4 The Stern–Gerlach force

A magnetic dipole $\boldsymbol\mu$ in a field $\mathbf B$ has energy $U = -\boldsymbol\mu\cdot\mathbf B$ (EM-5's $\mathbf m$; we write $\boldsymbol\mu$ as atomic physics does). The force is $-\nabla U$:
$$\mathbf F = \nabla(\boldsymbol\mu\cdot\mathbf B),\qquad F_z\approx\mu_z\frac{\partial B_z}{\partial z}. \tag{3.6}$$
The approximation: the magnet is built so that $\mathbf B\approx B_z(z)\hat{\mathbf z}$ along the beam, with a large $\partial B_z/\partial z$. ($\nabla\cdot\mathbf B = 0$ forces some $\partial B_x/\partial x$ as well, but the moment precesses rapidly about the strong $z$ field, so $\mu_x$ and $\mu_y$ average to zero and only $\mu_z\,\partial_zB_z$ survives.) A uniform field exerts a torque but no force; the gradient is the point.

**Classical prediction.** If the atoms' moments have magnitude $\mu$ and random orientation, then $\mu_z = \mu\cos\alpha$ with $\cos\alpha$ uniformly distributed on $[-1, 1]$ (the area element on a sphere is $d(\cos\alpha)\,d\varphi$). The deflection is proportional to $\mu_z$, so the beam spreads into a uniform band between $\pm z_{\max}$.

**Observation.** Two spots, at $\pm z_{\max}$. Each atom has $\mu_z = +\mu_B$ or $-\mu_B$, where $\mu_B = e\hbar/2m = 9.274\times10^{-24}$ J/T is the Bohr magneton. For an electron's spin $\boldsymbol\mu\approx-(e/m)\mathbf S$ (gyromagnetic factor $g\approx2$, QM-5), so $\mu_z = \mp\mu_B$ means $S_z = \pm\hbar/2$.

**Size of the effect.** With illustrative values of the right order (gradient $10^3$ T/m, a 3.5 cm magnet, atoms at 500 m/s): $F = \mu_B\times10^3 = 9.27\times10^{-21}$ N; a silver atom ($m_{\rm Ag} = 107.87$ u) accelerates at $5.18\times10^4$ m/s$^2$ for $7\times10^{-5}$ s and is deflected by $\pm0.127$ mm. A fraction of a millimetre, which is why the experiment was hard.

### 3.5 The Compton shift

▶ **Try it first.** A photon of wavelength $\lambda$ (energy $hc/\lambda$, momentum $h/\lambda$) hits an electron at rest and leaves at angle $\theta$ with wavelength $\lambda'$. Write energy conservation and momentum conservation. Eliminate the electron's final momentum using $E^2 = p^2c^2 + m^2c^4$ (EM-8 (3.25)).

**Theorem 3.2 (Compton).**
$$\lambda' - \lambda = \frac{h}{mc}(1 - \cos\theta). \tag{3.7}$$

*Derivation.* Let $p = h/\lambda$, $p' = h/\lambda'$, and let the electron leave with momentum $\mathbf p_e$ and energy $E_e$.
1. Momentum: $\mathbf p_e = \mathbf p - \mathbf p'$, so by the law of cosines $p_e^2 = p^2 + p'^2 - 2pp'\cos\theta$.
2. Energy: $pc + mc^2 = p'c + E_e$, so $E_e = (p - p')c + mc^2$.
3. Square step 2 and use $E_e^2 = p_e^2c^2 + m^2c^4$: $(p - p')^2c^2 + 2(p - p')mc^3 + m^2c^4 = (p^2 + p'^2 - 2pp'\cos\theta)c^2 + m^2c^4$.
4. Cancel $m^2c^4$, expand $(p - p')^2 = p^2 - 2pp' + p'^2$, cancel $p^2 + p'^2$, divide by $2c^2$: $-pp' + (p - p')mc = -pp'\cos\theta$, i.e. $(p - p')mc = pp'(1 - \cos\theta)$.
5. Divide by $pp'mc$: $\frac1{p'} - \frac1p = \frac{1 - \cos\theta}{mc}$. Multiply by $h$: (3.7). ∎ (`project.py` checks this with SymPy.)

■ **What just happened.** Treating light as a particle with momentum $h/\lambda$ and the electron relativistically gives a shift set by one length, the **Compton wavelength** $h/mc = 2.426$ pm, independent of the incident wavelength. That is why Compton needed X-rays: for visible light ($\sim500$ nm) the fractional shift is $5\times10^{-6}$; for 71 pm X-rays it is 3.4%.

### 3.6 de Broglie waves and Bohr's rule

de Broglie's hypothesis: a particle of momentum $p$ has wavelength
$$\lambda = \frac{h}{p}. \tag{3.8}$$
For a circular orbit, require the wave to fit around the circumference without mismatch: $2\pi r = n\lambda = nh/p$. Then $L = pr = nh/2\pi = n\hbar$, Bohr's condition. Quantization becomes the same thing as the discrete overtones of a guitar string: a boundary condition on a wave. That is exactly how Schrödinger's equation produces discrete levels (QM-3), and why wave mechanics came first. For a nonrelativistic particle of kinetic energy $E$, $p = \sqrt{2mE}$; Example 4.5 puts in numbers.

### 3.7 The postulates in vector-space language

Everything in §§3.1–3.6 is a patch. Here is the theory that replaced the patches, stated for a finite-dimensional state space $\mathbb C^N$. QM-2 through QM-4 extend it to infinite dimensions; nothing below changes except that sums become integrals and some care is needed about which vectors are allowed.

**Postulate 1 (states).** The state of an isolated system is a unit vector $|\psi\rangle$ in a complex inner product space $\mathcal H$ (LA-8 §3.1), with $\langle\psi|\psi\rangle = 1$. The vectors $|\psi\rangle$ and $e^{i\alpha}|\psi\rangle$ represent the same physical state (a **global phase** is unobservable; Exercise A4 proves it from Postulate 3). *Relative* phases between components are observable.

**Postulate 2 (observables).** Each measurable quantity is a Hermitian operator $\hat A = \hat A^\dagger$ on $\mathcal H$. By the spectral theorem (LA-8 Theorem 3.12) it has real eigenvalues $a_i$ and an orthonormal eigenbasis $|a_i\rangle$, and $\hat A = \sum_ia_i|a_i\rangle\langle a_i|$ (LA-8 (3.10) in Dirac form).

**Postulate 3 (Born rule).** A measurement of $\hat A$ gives one of the eigenvalues $a_i$, with probability
$$p_i = |\langle a_i|\psi\rangle|^2. \tag{3.9}$$
These are a probability distribution: $p_i\ge0$, and $\sum_ip_i = \sum_i|\langle a_i|\psi\rangle|^2 = \|\psi\|^2 = 1$ by Parseval (LA-8 Theorem 3.6). This is why the eigenbasis must be orthonormal and complete, and why observables must be Hermitian: the spectral theorem is the consistency condition for the Born rule. (If an eigenvalue is repeated, $p = \langle\psi|P|\psi\rangle$ with $P$ the projector onto its eigenspace.)

**Postulate 4 (collapse).** Immediately after the outcome $a_i$, the state is
$$|\psi\rangle\ \longmapsto\ \frac{|a_i\rangle\langle a_i|\psi\rangle}{|\langle a_i|\psi\rangle|} \tag{3.10}$$
which is $|a_i\rangle$ up to a phase: the projection (LA-8 Theorem 3.13) of $|\psi\rangle$ onto the eigenvector, renormalized. A repeated measurement then gives $a_i$ with probability $|\langle a_i|a_i\rangle|^2 = 1$ (experiment 1 of §2).

**Postulate 5 (time evolution).** Between measurements, the state evolves by the **Schrödinger equation**, generated by the Hermitian energy operator $\hat H$ (the Hamiltonian):
$$i\hbar\frac{d}{dt}|\psi(t)\rangle = \hat H|\psi(t)\rangle,\qquad |\psi(t)\rangle = e^{-i\hat Ht/\hbar}|\psi(0)\rangle\quad(\hat H\text{ time-independent}). \tag{3.11}$$
The solution is LA-7's matrix exponential (LA-7 §3.9) for the linear system $\dot{\mathbf x} = A\mathbf x$ with $A = -i\hat H/\hbar$.

**Postulate 6 (composite systems).** The state space of a system made of parts with spaces $\mathcal H_A$ and $\mathcal H_B$ is the tensor product $\mathcal H_A\otimes\mathcal H_B$ (LA-12), of dimension $\dim\mathcal H_A\cdot\dim\mathcal H_B$. QM-7 develops this; §3.10 uses it once, for which-path tags.

▶ **Try it first.** Show that $U(t) = e^{-i\hat Ht/\hbar}$ is unitary. You need $(e^{M})^\dagger = e^{M^\dagger}$ (take the dagger of the power series term by term) and that $e^{M}e^{-M} = I$.

**Theorem 3.3 (evolution is unitary).** If $\hat H^\dagger = \hat H$, then $U(t)^\dagger U(t) = I$; hence $\|\psi(t)\| = \|\psi(0)\|$ and total probability is conserved.

*Proof.* $U^\dagger = (e^{-i\hat Ht/\hbar})^\dagger = e^{+i\hat H^\dagger t/\hbar} = e^{i\hat Ht/\hbar}$, and since $i\hat Ht/\hbar$ commutes with $-i\hat Ht/\hbar$, $U^\dagger U = e^{0} = I$. ∎

That one line is why Hermiticity of $\hat H$ matters twice: real energies (Postulate 2) and probability conservation (Theorem 3.3). By LA-8 Theorem 3.20, $U$ preserves every inner product, so probabilities computed from $\langle a_i|\psi(t)\rangle$ always sum to 1.

▶ **Try it first.** The average outcome of many measurements of $\hat A$ is $\sum_ia_ip_i$. Write $p_i$ by (3.9) as $\langle\psi|a_i\rangle\langle a_i|\psi\rangle$ and pull the sum inside.

**Theorem 3.4 (expectation value).** In the state $|\psi\rangle$, the mean outcome of $\hat A$ is
$$\langle\hat A\rangle = \sum_ia_i|\langle a_i|\psi\rangle|^2 = \langle\psi|\hat A|\psi\rangle. \tag{3.12}$$

*Proof.* $\sum_ia_i|\langle a_i|\psi\rangle|^2 = \sum_ia_i\langle\psi|a_i\rangle\langle a_i|\psi\rangle = \langle\psi|\big(\sum_ia_i|a_i\rangle\langle a_i|\big)|\psi\rangle = \langle\psi|\hat A|\psi\rangle$, using $\overline{\langle a_i|\psi\rangle} = \langle\psi|a_i\rangle$ and the spectral decomposition of Postulate 2. ∎

■ **What just happened.** The six postulates are LA-8 with physical labels. Hermitian ↔ real outcomes and orthonormal eigenbasis; Parseval ↔ probabilities sum to 1; projector ↔ collapse; Hermitian generator ↔ unitary evolution. The physics is in two places: which operator corresponds to which quantity, and which $\hat H$ governs the system. For spin, experiment fixes both (next section). For a particle on a line, QM-2 has to work them out.

### 3.8 Spin-1/2: the complete worked instance

The state space is $\mathbb C^2$. The spin observables are
$$\hat S_x = \frac\hbar2\sigma_x = \frac\hbar2\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad\hat S_y = \frac\hbar2\sigma_y = \frac\hbar2\begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad\hat S_z = \frac\hbar2\sigma_z = \frac\hbar2\begin{pmatrix}1&0\\0&-1\end{pmatrix}. \tag{3.13}$$
Why these: each must be Hermitian with eigenvalues $\pm\hbar/2$ (the two spots), and the three eigenbases must be mutually "unbiased" (every eigenvector of one has probability $\frac12$ along every eigenvector of another, experiment 2 for every pair of perpendicular axes). LA-8 Proposition 3.22 shows the Pauli matrices have these properties. QM-5 derives (3.13) from the commutation relations alone; here, the matrices are justified by reproducing the data.

**Eigenvectors.** By Proposition 3.22, $\hat S_z$ has $|{\uparrow}\rangle$ ($+\hbar/2$) and $|{\downarrow}\rangle$ ($-\hbar/2$), and
$$|{\pm x}\rangle = \frac{1}{\sqrt2}\big(|{\uparrow}\rangle\pm|{\downarrow}\rangle\big),\qquad|{\pm y}\rangle = \frac1{\sqrt2}\big(|{\uparrow}\rangle\pm i|{\downarrow}\rangle\big). \tag{3.14}$$
Check one: $\hat S_x|{+x}\rangle = \frac\hbar2\frac1{\sqrt2}\sigma_x(1, 1)^\top = \frac\hbar2\frac1{\sqrt2}(1, 1)^\top$. ✓

▶ **Try it first.** Predict experiments 1–3 of §2 from (3.9), (3.10), (3.14). Each probability is one inner product.

**The data, derived.**
1. *Experiment 1.* The SGz$+$ filter leaves $|{\uparrow}\rangle$ (Postulate 4). $P(S_z = +\hbar/2) = |\langle{\uparrow}|{\uparrow}\rangle|^2 = 1$.
2. *Experiment 2.* $\langle{+x}|{\uparrow}\rangle = \frac1{\sqrt2}(1\cdot1 + 1\cdot0) = \frac1{\sqrt2}$ (conjugating the first slot changes nothing here, the entries are real). $P(S_x = +\hbar/2) = \frac12$, and $P(-) = |\langle{-x}|{\uparrow}\rangle|^2 = \frac12$. The sum is 1, by Parseval.
3. *Experiment 3.* After the SGx$+$ filter the state is $|{+x}\rangle$ by (3.10). $P(S_z = -\hbar/2) = |\langle{\downarrow}|{+x}\rangle|^2 = |1/\sqrt2|^2 = \frac12$. The $z$-information is gone because $|{+x}\rangle$ does not remember which $z$-state it was projected from: $|{+x}\rangle$ is the same vector whether it came from $|{\uparrow}\rangle$ or from $|{\downarrow}\rangle$.
4. *Recombined analyzer.* If both SGx beams are kept and recombined coherently with nothing recording the path, no measurement occurred. The apparatus acts as $|{+x}\rangle\langle{+x}| + |{-x}\rangle\langle{-x}| = I$ (completeness, LA-8 §3.12), the state stays $|{\uparrow}\rangle$, and $P(S_z = -\hbar/2) = 0$.
5. *Any angle.* An SG$\hat{\mathbf n}$ with $\hat{\mathbf n}$ at angle $\theta$ from $z$ in the $xz$-plane measures $\hat S_n = \hat{\mathbf n}\cdot\hat{\mathbf S} = \frac\hbar2\hat{\mathbf n}\cdot\boldsymbol\sigma$, with $|{+n}\rangle = \cos\frac\theta2|{\uparrow}\rangle + \sin\frac\theta2|{\downarrow}\rangle$ (Proposition 3.22 with $\phi = 0$). So $P(+n\mid{\uparrow}) = \cos^2(\theta/2)$: 1 at $\theta = 0$, $\frac12$ at $90°$, 0 at $180°$.

■ **What just happened.** Two-by-two Hermitian matrices, one inner product per probability, reproduce every Stern–Gerlach result, including the one with no classical explanation. "A qubit is a spin-1/2" is the literal statement that QC-1's state space and gates are this $\mathbb C^2$ with these operators; QC-1 starts from here. Example 4.6 runs the numbers and `project.py` Part 3 runs $10^4$ simulated atoms.

### 3.9 The general spin state and the Bloch sphere

▶ **Try it first.** A general vector in $\mathbb C^2$ has four real parameters. Normalization removes one and the global phase another. Show that every state can be written with two angles.

**Theorem 3.5 (Bloch parametrization).** Every spin-1/2 state is, up to global phase,
$$|\psi\rangle = \cos\frac\theta2|{\uparrow}\rangle + e^{i\phi}\sin\frac\theta2|{\downarrow}\rangle,\qquad0\le\theta\le\pi,\ 0\le\phi<2\pi, \tag{3.15}$$
and this is the $+\hbar/2$ eigenvector of $\hat{\mathbf n}\cdot\hat{\mathbf S}$ with $\hat{\mathbf n} = (\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$. Moreover $\langle\psi|\hat{\mathbf S}|\psi\rangle = \frac\hbar2\hat{\mathbf n}$.

*Derivation.*
1. Write $|\psi\rangle = a|{\uparrow}\rangle + b|{\downarrow}\rangle$ with $|a|^2 + |b|^2 = 1$. Multiply by the global phase $e^{-i\arg a}$ so that $a\ge0$ real.
2. Since $a\in[0, 1]$ and $|b| = \sqrt{1 - a^2}$, set $a = \cos(\theta/2)$, $|b| = \sin(\theta/2)$ with $\theta\in[0, \pi]$, and $b = e^{i\phi}\sin(\theta/2)$.
3. That (3.15) is the $+$ eigenvector of $\hat{\mathbf n}\cdot\boldsymbol\sigma$ is LA-8 Proposition 3.22, step 5.
4. Expectations by (3.12): $\langle\sigma_z\rangle = \cos^2\frac\theta2 - \sin^2\frac\theta2 = \cos\theta$. $\langle\sigma_x\rangle = \bar ab + \bar ba = 2\operatorname{Re}(\bar ab) = 2\cos\frac\theta2\sin\frac\theta2\cos\phi = \sin\theta\cos\phi$. $\langle\sigma_y\rangle = -i\bar ab + i\bar ba = 2\operatorname{Im}(\bar ab) = \sin\theta\sin\phi$. So $\langle\boldsymbol\sigma\rangle = \hat{\mathbf n}$. ∎

■ **What just happened.** Every spin-1/2 state is "spin up along some axis." The state space, rays in $\mathbb C^2$, is a sphere, the **Bloch sphere**: the state (3.15) sits at the point $\hat{\mathbf n}$. Orthogonal states sit at antipodes ($|{\uparrow}\rangle$ at the north pole, $|{\downarrow}\rangle$ at the south; $|{\pm x}\rangle$ at $\pm\hat{\mathbf x}$), which is why the half-angle appears. The probability of $+n'$ for a state at $\hat{\mathbf n}$ is $\cos^2(\gamma/2) = \frac12(1 + \hat{\mathbf n}\cdot\hat{\mathbf n}')$, with $\gamma$ the angle between the axes (Exercise B4). `project.py` Part 4 draws the eigenvectors of several $\hat{\mathbf n}\cdot\boldsymbol\sigma$ on the sphere and checks $\langle\boldsymbol\sigma\rangle = \hat{\mathbf n}$, completing LA-8's picture. For polarization the same sphere is EM-7's Poincaré sphere, with polarizer angle $\alpha$ at sphere angle $2\alpha$.

**Precession, as an example of Postulate 5.** An electron spin in a field $B\hat{\mathbf z}$ has $\hat H = \omega_0\hat S_z$ with $\omega_0 = g\mu_BB/\hbar$. Then $e^{-i\hat Ht/\hbar} = \operatorname{diag}(e^{-i\omega_0t/2}, e^{i\omega_0t/2})$ (diagonal, so LA-7's exponential is entrywise), and (3.15) becomes $\cos\frac\theta2|{\uparrow}\rangle + e^{i(\phi + \omega_0t)}\sin\frac\theta2|{\downarrow}\rangle$ after removing a global phase: $\phi\to\phi + \omega_0t$. The Bloch vector precesses about $\hat{\mathbf z}$ at $\omega_0$, which is the classical Larmor precession of a magnetic moment. Exercise B7 puts in numbers; QM-8's Rabi oscillations and QC-1's gates are this computation with other $\hat H$.

### 3.10 The double slit in vector language

Let $\phi_1(x)$ and $\phi_2(x)$ be the amplitudes to arrive at screen position $x$ via slit 1 and slit 2 (complex numbers, from the wave optics of EM-7). Rule (ii) of §1.7: with no record of the path, the amplitude is $\phi_1 + \phi_2$ and the detection probability density is $|\phi_1 + \phi_2|^2$.

▶ **Try it first.** Now suppose each electron also flips a "tag" system (a photon scattered near the slit, a spin in a detector) into state $|t_1\rangle$ or $|t_2\rangle$ depending on its path. By Postulate 6 the joint state at position $x$ is $\phi_1(x)|t_1\rangle + \phi_2(x)|t_2\rangle$. The screen does not look at the tag, so add the probabilities of both tag outcomes. Compute.

**Theorem 3.6 (which-path information).** With tags $|t_1\rangle$, $|t_2\rangle$ (unit vectors),
$$P(x) = |\phi_1(x)|^2 + |\phi_2(x)|^2 + 2\operatorname{Re}\big[\overline{\phi_1(x)}\phi_2(x)\langle t_1|t_2\rangle\big]. \tag{3.16}$$

*Derivation.*
1. Take any orthonormal basis $|k\rangle$ of the tag space. The probability of position $x$ and tag outcome $k$ is $|\langle k|\Phi_x\rangle|^2$ with $|\Phi_x\rangle = \phi_1(x)|t_1\rangle + \phi_2(x)|t_2\rangle$.
2. Summing over $k$: $P(x) = \sum_k|\langle k|\Phi_x\rangle|^2 = \|\Phi_x\|^2$ (Parseval again).
3. Expand: $\|\Phi_x\|^2 = |\phi_1|^2\langle t_1|t_1\rangle + |\phi_2|^2\langle t_2|t_2\rangle + \bar\phi_1\phi_2\langle t_1|t_2\rangle + \bar\phi_2\phi_1\langle t_2|t_1\rangle$, and the last two terms are complex conjugates. ∎

■ **What just happened.** Identical tags ($\langle t_1|t_2\rangle = 1$, no path information) give $|\phi_1 + \phi_2|^2$: full interference. Orthogonal tags (perfect path information) give $|\phi_1|^2 + |\phi_2|^2$: the classical sum, rule (iii). In between, the fringe visibility is exactly $|\langle t_1|t_2\rangle|$ (Exercise B8). Nothing "disturbed" the electron's motion; the interference term was multiplied by the overlap of two states of *another* system. This is the simplest instance of decoherence, which QM-7 and QM-9 develop: the environment acquiring path information is what makes quantum superpositions look classical. `project.py` Part 5 builds the joint state explicitly and samples electrons one at a time.

### 3.11 Commutators and incompatibility

▶ **Try it first.** Multiply $\sigma_x\sigma_y$ and $\sigma_y\sigma_x$. What is the difference?

**Proposition 3.7.** $[\hat S_x, \hat S_y] = i\hbar\hat S_z$, and cyclically $[\hat S_y, \hat S_z] = i\hbar\hat S_x$, $[\hat S_z, \hat S_x] = i\hbar\hat S_y$.

*Derivation.* $\sigma_x\sigma_y = \begin{pmatrix}0&1\\1&0\end{pmatrix}\begin{pmatrix}0&-i\\i&0\end{pmatrix} = \begin{pmatrix}i&0\\0&-i\end{pmatrix} = i\sigma_z$, and $\sigma_y\sigma_x = \begin{pmatrix}0&-i\\i&0\end{pmatrix}\begin{pmatrix}0&1\\1&0\end{pmatrix} = \begin{pmatrix}-i&0\\0&i\end{pmatrix} = -i\sigma_z$. So $[\sigma_x, \sigma_y] = 2i\sigma_z$ and $[\hat S_x, \hat S_y] = \frac{\hbar^2}{4}2i\sigma_z = i\hbar\hat S_z$. The other two follow by the same two multiplications with the indices relabeled (each is written out in Exercise B6). ∎

**Why the commutator measures incompatibility.** If $[\hat A, \hat B] = 0$ for Hermitian $\hat A$, $\hat B$, they have a common orthonormal eigenbasis (each preserves the other's eigenspaces; QM-4 proves it), so a state can have definite values of both and measuring one does not disturb the other. Experiment 3 says $\hat S_x$ and $\hat S_z$ have no common eigenvector, and Proposition 3.7 says they do not commute. The quantitative statement is the uncertainty principle.

**Theorem 3.8 (Robertson uncertainty relation; derived carefully in QM-4).** For Hermitian $\hat A$, $\hat B$ and any state,
$$\Delta A\,\Delta B\ \ge\ \tfrac12\big|\langle[\hat A, \hat B]\rangle\big|. \tag{3.17}$$

*Sketch.* Let $|f\rangle = (\hat A - \langle\hat A\rangle)|\psi\rangle$ and $|g\rangle = (\hat B - \langle\hat B\rangle)|\psi\rangle$, so $\|f\|^2 = (\Delta A)^2$ and $\|g\|^2 = (\Delta B)^2$. Cauchy–Schwarz (LA-8 Theorem 3.4) gives $\|f\|^2\|g\|^2\ge|\langle f|g\rangle|^2\ge(\operatorname{Im}\langle f|g\rangle)^2$, and $\langle f|g\rangle - \langle g|f\rangle = \langle\psi|[\hat A, \hat B]|\psi\rangle$, so $\operatorname{Im}\langle f|g\rangle = \frac1{2i}\langle[\hat A, \hat B]\rangle$. QM-4 fills in the steps and adds the Schrödinger term from the real part.

For spin, (3.17) with Proposition 3.7 reads $\Delta S_x\Delta S_y\ge\frac\hbar2|\langle\hat S_z\rangle|$. In $|{\uparrow}\rangle$: $\langle\hat S_x^2\rangle = \hbar^2/4$ (since $\sigma_x^2 = I$) and $\langle\hat S_x\rangle = 0$, so $\Delta S_x = \hbar/2$; likewise $\Delta S_y = \hbar/2$. The product is $\hbar^2/4$ and the bound is $\frac\hbar2\cdot\frac\hbar2 = \hbar^2/4$: equality. Exercise B6 shows the inequality for every state reduces to $n_x^2n_y^2\ge0$; `project.py` Part 6 checks it on 2000 random states.

**Correspondence with CM-8.** The classical angular momentum components satisfy the Poisson bracket relations $\{L_x, L_y\} = L_z$ (CM-8 §3.6). Dirac's 1925 rule is
$$\{A, B\}\ \longrightarrow\ \frac{1}{i\hbar}[\hat A, \hat B], \tag{3.18}$$
and it turns $\{L_x, L_y\} = L_z$ into Proposition 3.7 and $\{q, p\} = 1$ into $[\hat x, \hat p] = i\hbar$ (QM-2, QM-4). Spin has no classical counterpart, but it obeys the same algebra as orbital angular momentum, which is why QM-5 treats both at once.

### 3.12 What comes next: the same postulates on functions

QM-2 and QM-3 apply §3.7 to a particle on a line. The state space becomes $L^2(\mathbb R)$, square-integrable complex functions with LA-8's inner product $\langle\phi|\psi\rangle = \int\bar\phi\psi\,dx$ (LA-8 §3.1's function-space example). The "components" of $|\psi\rangle$ in the position "basis" are the values $\psi(x)$, and $|\psi(x)|^2dx$ is the Born probability of finding the particle in $[x, x + dx]$. The Hamiltonian is a differential operator, $\hat H = -\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V(x)$, and (3.11) becomes a PDE, the Schrödinger equation you will solve with separation of variables and boundary conditions.

Why is this order historical? Schrödinger's wave equation was something physicists in 1926 knew how to solve, it gave hydrogen directly, and spin was not even recognized as a separate two-dimensional degree of freedom until Pauli introduced his matrices in 1927. Why do we do wave mechanics next, when spin is simpler? Because positions, momenta, atoms, tunnelling and bound states, the things QM is used for in mechanical and materials physics, live there; because you are strong at PDEs; and because the infinite-dimensional case has genuine new features ($\hat x$ and $\hat p$ have no normalizable eigenvectors; $[\hat x, \hat p] = i\hbar$ cannot hold for finite matrices, since the trace of a commutator is 0 and the trace of $i\hbar I$ is not). This module's payoff comes in QM-4: when the formalism is stated in general, you will already know it.

### 3.13 What the postulates say and do not say

The postulates are a complete recipe for predictions: given $|\psi\rangle$, $\hat H$ and the observable, they give every probability, and those probabilities have never failed a test. What they do not say is what a measurement is. Postulate 5 says isolated systems evolve unitarily and continuously. Postulate 4 says measurement causes a discontinuous, non-unitary jump. A measuring apparatus is made of atoms, so it should obey Postulate 5 too; then where does the jump happen, and why? That is the **measurement problem**, and it is unresolved in the sense that physicists who agree on every prediction disagree about it. The modern partial answer is **decoherence** (QM-9): when a system interacts with an environment, the environment records which-path information, and by Theorem 3.6 the interference terms between macroscopically different outcomes get multiplied by overlaps that are astronomically small. Decoherence explains why you never see a superposition of a dead and live cat, and why the Born-rule probabilities behave classically. It does not, on its own, explain why one outcome occurs. Interpretations (Copenhagen, many-worlds, Bohmian mechanics, objective-collapse models) differ on that point and agree on the numbers in this module.

The advice: compute first, philosophize later. Every question in QM-1 to QM-8 can be answered without taking a position. QM-9 returns to it with the tools to say something precise.

## 4. Worked examples

### Example 4.1: Rayleigh–Jeans versus Planck at 500 nm and 5000 K (project Part 1)

$x = hc/(\lambda k_BT) = (6.626\times10^{-34})(2.998\times10^8)/\big[(500\times10^{-9})(1.381\times10^{-23})(5000)\big] = 5.7551$.

Planck (3.2): $8\pi hc/\lambda^5 = 8\pi(1.986\times10^{-25})/(3.125\times10^{-32}) = 1.598\times10^{8}$ J/m$^4$; $e^{5.7551} - 1 = 314.8$; $u_\lambda = 5.075\times10^5$ J/m$^4$ (that is, $5.08\times10^{-4}$ J/m$^3$ per nm of bandwidth).

Rayleigh–Jeans (3.1): $8\pi k_BT/\lambda^4 = 8\pi(6.904\times10^{-20})/(6.25\times10^{-26}) = 2.776\times10^7$ J/m$^4$.

Ratio: $2.776\times10^7/5.075\times10^5 = 54.70 = (e^x - 1)/x$. At green light and the temperature of the Sun's surface, classical physics is off by a factor of 55, and the error grows exponentially toward the ultraviolet. At 50 μm ($x = 0.058$) the two agree within 3%. The Planck curve at 5000 K peaks at 579.6 nm. Figure `blackbody.png` shows all three laws.

### Example 4.2: The photoelectric stopping potential for sodium under 400 nm light (project Part 1)

Photon energy: $hc/\lambda = 1239.84/400 = 3.0996$ eV. Work function $W = 2.3$ eV. By (3.3), $eV_s = 3.0996 - 2.3 = 0.7996$ eV, so $V_s = 0.80$ V. The fastest electrons have speed $\sqrt{2(0.7996)(1.602\times10^{-19})/(9.109\times10^{-31})} = 5.30\times10^5$ m/s. Threshold: $f_0 = W/h = 5.56\times10^{14}$ Hz, $\lambda_0 = 1239.84/2.3 = 539.1$ nm, so green light at 540 nm, however intense, ejects nothing.

`project.py` generates stopping potentials at nine wavelengths from 250 to 510 nm with 0.02 V noise and fits a line: slope$\times e$ gives $h = 6.61\times10^{-34}$ J s (true $6.626\times10^{-34}$, within 0.3%) and intercept $W = 2.29$ eV. That is Millikan's measurement in miniature (`photoelectric.png`).

### Example 4.3: Bohr levels and the Balmer series (project Part 2)

From (3.4) with CODATA constants: $a_0 = 0.052918$ nm, $E_R = 13.606$ eV, $v_1/c = \alpha = 1/137.04$, $R_\infty = E_R/hc = 1.097373\times10^7$ m$^{-1}$.

Levels: $E_1 = -13.606$, $E_2 = -3.401$, $E_3 = -1.512$ eV. $H\alpha$ is $3\to2$: $1/\lambda = R_\infty(\frac14 - \frac19) = \frac{5}{36}R_\infty$, $\lambda = 36/(5R_\infty) = 656.11$ nm.

Refinements, to see how good the model is: the reduced-mass Rydberg $R_H = R_\infty/(1 + m/m_p) = 1.096776\times10^7$ m$^{-1}$ gives $656.47$ nm in vacuum; dividing by the refractive index of air, 1.000277, gives 656.29 nm. The measured $H\alpha$ line in air is at 656.28 nm. The Balmer series $n\to2$ for $n = 3, \dots, 7$: 656.1, 486.0, 433.9, 410.1, 396.9 nm, converging on the series limit $4/R_\infty = 364.5$ nm. Lyman ($\to1$) lies entirely in the ultraviolet (121.5 nm down to 91.1 nm); Paschen ($\to3$) in the infrared (limit 820.1 nm). `bohr_series.png` plots the levels and the three series on a wavelength axis.

### Example 4.4: The Compton shift at 90° for molybdenum X-rays (project Part 2)

$h/mc = 6.626\times10^{-34}/\big[(9.109\times10^{-31})(2.998\times10^8)\big] = 2.4263$ pm. At $\theta = 90°$, $1 - \cos\theta = 1$, so by (3.7) $\lambda' = 71.07 + 2.4263 = 73.496$ pm, a 3.41% shift. Photon energies: $hc/\lambda = 1239.84/0.07107 = 17.445$ keV before, $16.869$ keV after. The electron recoils with kinetic energy $17.445 - 16.869 = 0.576$ keV $= 576$ eV, far above the few-eV binding of graphite's outer electrons, which is why treating them as free works. The maximum shift, at $180°$, is $2h/mc = 4.85$ pm (Your turn 1 plots $\Delta\lambda(\theta)$).

### Example 4.5: de Broglie wavelengths: an electron at 100 eV and a baseball (project Part 2)

Electron: $p = \sqrt{2mE} = \sqrt{2(9.109\times10^{-31})(100)(1.602\times10^{-19})} = 5.40\times10^{-24}$ kg m/s; $\lambda = h/p = 0.1226$ nm. That is the spacing of atoms in a crystal, so crystals diffract 100 eV electrons. (The relativistic correction at 100 eV changes $\lambda$ by less than $10^{-4}$, since $E\ll mc^2 = 511$ keV.)

Davisson–Germer check: at 54 eV, $\lambda = h/\sqrt{2mE} = 0.1669$ nm. Nickel's surface rows are spaced $d = 0.215$ nm, and the grating condition $d\sin\theta = \lambda$ at the observed $\theta = 50°$ gives $0.215\times0.7660 = 0.1647$ nm. Agreement to 1.3%.

Baseball, 0.145 kg at 40 m/s: $\lambda = 6.626\times10^{-34}/5.8 = 1.14\times10^{-34}$ m, about $10^{-19}$ times a proton's diameter. No slit could ever show it. This is why everyday objects look classical: not because quantum mechanics stops applying, but because $\lambda$ is absurdly small compared to every length in the experiment.

### Example 4.6: The Stern–Gerlach sequence $S_z\to S_x\to S_z$, all probabilities (project Part 3)

Start: $10^4$ atoms leave an SGz with the $-$ port blocked, so each is in $|{\uparrow}\rangle = (1, 0)^\top$.

*Stage 1, SGx.* Eigenvectors from (3.14): $|{+x}\rangle = (1, 1)^\top/\sqrt2$, $|{-x}\rangle = (1, -1)^\top/\sqrt2$, eigenvalues $\pm\hbar/2$. Amplitudes: $\langle{+x}|{\uparrow}\rangle = \frac1{\sqrt2}$, $\langle{-x}|{\uparrow}\rangle = \frac1{\sqrt2}$. Probabilities $\frac12$, $\frac12$; sum 1. Expected: 5000 atoms through $+x$. Block $-x$. The survivors are in $|{+x}\rangle$ by (3.10).

*Stage 2, SGz.* $\langle{\uparrow}|{+x}\rangle = \frac1{\sqrt2}$, $\langle{\downarrow}|{+x}\rangle = \frac1{\sqrt2}$: probabilities $\frac12$, $\frac12$. Expected: 2500 atoms with $S_z = +\hbar/2$ and 2500 with $-\hbar/2$.

*Overall:* $P(\text{end in }{\downarrow}\mid\text{start }{\uparrow}) = \frac12\cdot\frac12 = \frac14$ of the original beam, from atoms that were filtered to be "up" at the start. Without the SGx stage, the same probability is $|\langle{\downarrow}|{\uparrow}\rangle|^2 = 0$. With the SGx stage present but both of its beams recombined and nothing recorded, it is also 0 (§3.8 step 4).

*Expectation values* by (3.12): in $|{\uparrow}\rangle$, $\langle\hat S_x\rangle = \frac\hbar2(1, 0)\sigma_x(1, 0)^\top = 0$; in $|{+x}\rangle$, $\langle\hat S_z\rangle = \frac\hbar2\cdot\frac12(1\cdot1 + 1\cdot(-1)) = 0$. Both agree with $\sum a_ip_i = \frac\hbar2(\frac12) - \frac\hbar2(\frac12) = 0$.

*Simulation* (one run, seed 0): 4946 of 10,000 pass SGx$+$; of those, 2423 (0.4899) exit SGz$-$. Both within 4 binomial standard deviations ($\sigma = 50$ for the first, 35 for the second). At $60°$ instead of $90°$: $P(+n\mid{\uparrow}) = \cos^2 30° = 0.75$, simulated 0.7505.

## 5. Why it matters and where it connects

- **QM-2, QM-3.** The postulates on $L^2(\mathbb R)$: Born's rule becomes $|\psi(x)|^2$; de Broglie's $\lambda = h/p$ becomes $\hat p = -i\hbar\,d/dx$ acting on $e^{ipx/\hbar}$; §3.6's standing waves become the infinite well's $\sin(n\pi x/L)$; §3.3's failed oscillator gets its $\frac12\hbar\omega$.
- **QM-4.** States the postulates in general, derives (3.17), adds the Heisenberg picture. Theorem 3.4 and §3.11 are its first pages.
- **QM-5.** Derives (3.13) from $[\hat S_x, \hat S_y] = i\hbar\hat S_z$ alone, and spin-1, spin-3/2, ... .
- **QM-6.** Solves hydrogen properly and recovers (3.4) exactly, with the right degeneracies and $\ell = 0$ ground state that Bohr lacked.
- **QM-7, QM-9.** Postulate 6, entanglement and Bell's inequality (Exercise C2's hidden-variable model is a one-particle warm-up); decoherence as the generalization of Theorem 3.6.
- **SM-5.** Rederives (3.2) by Bose–Einstein counting; the photon gas.
- **QC-1.** A qubit is this module's $\mathbb C^2$: $|0\rangle = |{\uparrow}\rangle$, measurement in the computational basis is SGz, the Bloch sphere is §3.9, a single-qubit gate is $e^{-i\hat Ht/\hbar}$ for some $2\times2$ $\hat H$ (§3.9's precession).
- **EM-7.** Jones vectors and polarizers are the same mathematics as Stern–Gerlach (Your turn 3); EM-8's radiating atom is the contradiction of §1.3.
- **CM-8.** Action variables are what Bohr–Sommerfeld quantized; Poisson brackets become commutators by (3.18).

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are drills with the experimental formulas and two-component spinors. B1–B9 are the module: B3 (Bohr–Sommerfeld's successes and failures), B5 (sequential analyzers at $60°$), B6 (the spin uncertainty product) and B8 (partial which-path information, extends `project.py`) are essential. C1 derives Wien's displacement law, C2 builds a classical hidden-variable model and watches it fail, C3 is Bohr's correspondence principle.

## 7. Coding project

`project.py` runs in about 2 seconds; figures go to `figures/`.

1. **Part 1:** Planck, Rayleigh–Jeans and Wien curves at 5000 K with Example 4.1's numbers and the Stefan–Boltzmann integral; a fitted $h$ from synthetic photoelectric data (Example 4.2).
2. **Part 2:** Bohr radius, energies, Rydberg constant from constants; the Lyman, Balmer, Paschen series on a log-wavelength axis; the classical collapse time; Compton (with a SymPy check of Theorem 3.2); de Broglie and Davisson–Germer; the Stern–Gerlach deflection; Bohr–Sommerfeld for the oscillator.
3. **Part 3, the Stern–Gerlach simulator:** `measure(psi, n, rng)` takes an array of states, computes Born probabilities onto the eigenvectors of $\hat{\mathbf n}\cdot\hat{\mathbf S}$, samples outcomes and collapses. $10^4$ atoms through every sequence of §3.8, with binomial asserts.
4. **Part 4:** the Bloch sphere with $|{\pm n}\rangle$ for five axes; $\langle\boldsymbol\sigma\rangle = \pm\hat{\mathbf n}$; global phase invisible.
5. **Part 5:** the double slit by (3.16) with tag overlaps 1, 0.5, 0 (visibilities 1.000, 0.500, 0.000), the joint state built explicitly, and 3000 electrons sampled one at a time for each case (`double_slit.png`).
6. **Part 6:** commutators numerically and in SymPy; the uncertainty product on 2000 random states; Larmor precession by `expm`, unitarity asserted.

**Your turn** (reference solutions included): the Compton shift versus angle; He$^+$ and the hydrogen-like $Z^2$ scaling, including the Pickering-line coincidence; three polarizers as a Stern–Gerlach sequence. Part 8 prints every number in `solutions.md`.

**Build it yourself:** before reading Part 3, write `measure` yourself: given a state and an axis, build $\hat{\mathbf n}\cdot\hat{\mathbf S}$, diagonalize it with `eigh`, compute $|\langle n_\pm|\psi\rangle|^2$ (watch the conjugate), draw the outcome, return the eigenvector. Then reproduce experiment 3 of §2 with $10^4$ atoms.

## 8. Where you'll struggle

- **The hardest part is not mathematical.** You will keep wanting the electron to *have* a position, or a spin direction, before it is measured, with the measurement merely revealing it. Experiment 3 of §2 and Exercise C2 are the evidence against the simplest version of that picture, and Bell's theorem (QM-7) against every local version. The fix is operational: say "the state is $|{+x}\rangle$" and "the probability of finding $S_z = +\hbar/2$ is $\frac12$," and stop saying "the spin is pointing along $x$, so its $z$-component is..." The second sentence has no meaning the theory supports.
- **"The electron is a wave" and "the electron is a particle" are both wrong.** The correct statement: the electron's state is a vector; the things you observe (a dot on a screen, a spot from a magnet) are measurement outcomes, eigenvalues, each occurring with a Born probability. The wave behaviour is the vector's components adding with phases; the particle behaviour is that each measurement gives one definite outcome. Neither picture is the object.
- **The Bohr model is wrong, and we taught it anyway.** It gets hydrogen's energies exactly and nearly everything else wrong (angular momentum of the ground state, helium, line intensities, the oscillator's zero-point energy, chaotic systems). It is in this module because it got the numbers right in 1913 and because you should see exactly how a theory that fits the data can still be wrong. Do not use its orbit picture for intuition about real atoms; QM-6's orbitals replace it.
- **Postulates before wavefunctions is unusual.** Griffiths and most courses start with $\Psi(x, t)$ and state the vector-space version in chapter 3. We do it the other way because you have LA-8, and the payoff is that QM-4 will be a relief instead of a shock. The cost is that the postulates here may feel abstract until QM-2 gives them a second example. Accept the abstraction for now and run `project.py` Part 3 until the rule "probability = squared inner product, then collapse" is automatic.
- **The conjugate in the first slot.** $\langle{+y}|{\uparrow}\rangle$: the bra conjugates, so $\langle{+y}| = \frac1{\sqrt2}(1, -i)$, and the inner product with $(1, 0)^\top$ is $\frac1{\sqrt2}$. Forget the conjugate and $\langle{+y}|{+y}\rangle$ comes out 0. In NumPy, `np.vdot(a, b)` conjugates `a`; `a @ b` does not.
- **Half-angles.** Orthogonal spin states are $180°$ apart on the Bloch sphere, not $90°$. $|{+x}\rangle$ and $|{\uparrow}\rangle$ are $90°$ apart on the sphere and have overlap $\cos 45°$, not 0. The formula $P = \cos^2(\theta/2)$ has the half-angle for this reason; polarizers have $\cos^2\alpha$ because polarization angle $\alpha$ sits at sphere angle $2\alpha$.
- **Global versus relative phase.** $-|{\uparrow}\rangle$ is the same state as $|{\uparrow}\rangle$. $|{\uparrow}\rangle + |{\downarrow}\rangle$ and $|{\uparrow}\rangle - |{\downarrow}\rangle$ are different states (they are $|{\pm x}\rangle$, orthogonal). Exercise A4 makes you prove both.
- **Interpretation will nag.** Section 3.13 is all you need until QM-9. If you find yourself spending an evening on many-worlds, go back to computing.

## 9. Check yourself

1. State the Rayleigh–Jeans and Planck laws and say which assumption in the derivation differs.
2. Why does the photoelectric stopping potential not depend on light intensity?
3. Derive $r_n = n^2a_0$ and $E_n = -13.6\text{ eV}/n^2$ from $L = n\hbar$.
4. Name three things the Bohr–Sommerfeld rule gets wrong.
5. Derive the Compton shift from energy and momentum conservation.
6. Show that de Broglie's standing-wave condition on a circle is $L = n\hbar$.
7. What does classical physics predict for the Stern–Gerlach pattern, and why?
8. State the six postulates, and for each of 2, 3, 4, 5 name the LA-8 theorem it relies on.
9. Prove that $e^{-i\hat Ht/\hbar}$ is unitary, and derive $\langle\hat A\rangle = \langle\psi|\hat A|\psi\rangle$ from the Born rule.
10. Compute $P(S_z = -\hbar/2)$ after SGz$+$, SGx$+$; and after SGz$+$ followed by an SGx whose beams are recombined.
11. Where on the Bloch sphere is $\frac1{\sqrt2}(|{\uparrow}\rangle + i|{\downarrow}\rangle)$?
12. Compute $[\hat S_x, \hat S_y]$ and state what it implies about $\Delta S_x\Delta S_y$ in the state $|{\uparrow}\rangle$.

## 10. Pair with

- **Feynman, *Lectures on Physics* III, ch. 1–3 (required).** The two-slit experiment with bullets, waves and electrons, the three rules, and probability amplitudes built from scratch; nobody explains interference and which-path information better. Ch. 5–6 do Stern–Gerlach with Feynman's recombining analyzers.
- **Townsend, *A Modern Approach to Quantum Mechanics*, ch. 1.** Stern–Gerlach first, exactly as in §2 and §3.8, with more sequential-analyzer experiments and the step from data to kets done slowly.
- **Griffiths & Schroeter, *Introduction to Quantum Mechanics* (3rd ed.), preface and ch. 1 opening.** The standard wavefunction-first route you will take in QM-2, and a short, sensible discussion of what the wavefunction means.
- **Eisberg & Resnick, *Quantum Physics of Atoms, Molecules, Solids, Nuclei, and Particles*, ch. 1–4.** The experiments of §§1.1–1.6 in historical and experimental detail: real data, apparatus, the Franck–Hertz experiment, more Bohr model than you need.
- **Zwiebach, *Mastering Quantum Mechanics*, ch. 1–2.** A modern account of the experiments and the linearity of quantum mechanics, with Mach–Zehnder interferometers and the Elitzur–Vaidman bomb tester; matches MIT 8.04's opening lectures.
