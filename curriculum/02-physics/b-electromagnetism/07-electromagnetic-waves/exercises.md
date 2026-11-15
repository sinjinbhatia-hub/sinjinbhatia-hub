# EM-7 Exercises

Conventions as in the notes: physics sign $e^{i(\mathbf k\cdot\mathbf r - \omega t)}$, Griffiths' sign for $r_p$, $\mu = \mu_0$ everywhere unless stated. Numbers you need: $c = 2.998\times10^8$ m/s, $\varepsilon_0 = 8.854\times10^{-12}$ F/m, $\mu_0 = 4\pi\times10^{-7}$ H/m, $e = 1.602\times10^{-19}$ C, $m_e = 9.109\times10^{-31}$ kg, $\hbar = 1.055\times10^{-34}$ J s.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** In vacuum, $\mathbf E(x, t) = (300\ \text{V/m})\,\hat{\mathbf y}\cos(kx + \omega t)$ with $\lambda = 500$ nm.
(a) Which way does the wave travel? Find $f$, $\omega$, $k$.
(b) Write $\mathbf B(x, t)$ as a vector, using (3.6).
(c) Find the intensity.

**A2.** A 5.0 mW laser pointer makes a beam of 1.0 mm diameter, uniform across the spot. Find the intensity, the amplitude $E_0$ and the rms field. What force does the beam exert on a mirror at normal incidence?

**A3.** (a) Unpolarized light passes through two polarizers with axes at $0°$ and $90°$. What fraction gets through? Insert a third polarizer at $45°$ between them. What fraction gets through now? Use the Jones matrices (3.10).
(b) Show that a half-wave plate with fast axis at $\theta$ has Jones matrix $\begin{pmatrix}\cos2\theta & \sin2\theta\\\sin2\theta & -\cos2\theta\end{pmatrix}$ (up to an overall phase) and that it turns linear polarization at angle $\varphi$ into linear polarization at $2\theta - \varphi$. What does $\theta = 22.5°$ do to horizontal light?

**A4.** For water, $n = 1.33$: find (a) the reflectance at normal incidence from air, (b) Brewster's angle from air, (c) the critical angle from inside the water. (d) A diver looking up sees the whole sky compressed into a bright circle overhead ("Snell's window"). Which of (a)–(c) sets the circle's angular radius, and why?

**A5.** Seawater has $\sigma\approx4$ S/m and $\varepsilon_r\approx81$. Check that it is a good conductor at 10 kHz, and find the skin depth at 10 kHz and at 1 MHz. Why do navies talk to submerged submarines at very low frequencies?

## B. Core (do all)

**B1. Energy at a string junction.** A wave $A\cos(kz - \omega t)$ on a string of tension $T$ and mass per length $\mu$ carries average power $P = \tfrac12\sqrt{T\mu}\,\omega^2A^2$ (derive it from the rate of work $-T\,\partial_zf\,\partial_tf$ done by the left side on the right).
(a) Using (3.4), show that the reflected plus transmitted power equals the incident power.
(b) Define the string's impedance $Z = \sqrt{T\mu}$. Rewrite (3.4) in terms of $Z_1$, $Z_2$, and compare with (3.13) written in terms of $Z = Z_0/n$. What plays the role of impedance for light?

**B2. The s-polarization Fresnel equations.** Derive (3.17) in full. Take $\mathbf E = \tilde E_0\hat{\mathbf y}$ for all three waves, compute each $\mathbf B$ from (3.6), write all four boundary conditions, show which are trivial and which repeat others, and solve. Then show $R_s + T_s = 1$ and that, for $\mu_1 = \mu_2$ and $n_1\ne n_2$, $r_s$ never vanishes.

**B3. Polarization by reflection.** Unpolarized light falls on glass ($n = 1.5$) at $45°$.
(a) Using Example 4.4, what fraction of the incident power is reflected?
(b) The reflected beam is partially polarized. Its degree of polarization is $(I_{\max} - I_{\min})/(I_{\max} + I_{\min})$ as an analyser is rotated. Compute it.
(c) Explain, with a sketch of the plane of incidence, why polarized sunglasses have vertical transmission axes, and why they are most effective for glare arriving at about $53°$ from the vertical off water.

**B4. The evanescent wave carries no energy across.** s-polarized light is totally internally reflected at an interface from $n_1$ to $n_2$.
(a) Using (3.6) with the complex transmitted wave vector $(k_{Tx}, 0, i\kappa)$, write the real fields $\mathbf E_T$ and $\mathbf B_T$ in medium 2.
(b) Show that the time-averaged Poynting vector in medium 2 has $\langle S_z\rangle = 0$ but $\langle S_x\rangle\ne0$.
(c) Glass ($n_1 = 1.5$) to air at $60°$, $\lambda_0 = 633$ nm: find the decay length, $r_s$ as a complex number, $|r_s|$, and its phase.

**B5. Plasmas and the ionosphere.**
(a) Derive $\varepsilon_r = 1 - \omega_p^2/\omega^2$ from (3.23) with $\omega_0 = 0$, $\gamma = 0$, and hence $\omega^2 = \omega_p^2 + c^2k^2$.
(b) Show $v_pv_g = c^2$, with $v_p > c > v_g$ above $\omega_p$.
(c) For the ionosphere with $N = 10^{12}$ m⁻³, find $f_p$. For a 20 MHz signal, find $n$, $v_g/c$, and the extra delay (relative to vacuum) of a pulse crossing 100 km of such plasma.
(d) The extra delay scales as $1/f^2$ for $f\gg f_p$. Show this, and explain how a receiver using two frequencies can measure and remove it (GPS does).

**B6. Rectangular waveguide modes.** For the guide of Example 4.6 ($a = 2.3$ cm, $b = 1.0$ cm):
(a) List the five lowest cutoff frequencies with their mode labels. Over what band is the guide single-mode?
(b) Show that TM modes need $m, n\ge1$, so TM$_{10}$ and TM$_{01}$ do not exist.
(c) Write the TE$_{10}$ fields explicitly from (3.31)–(3.32). Show that they are a superposition of two plane waves bouncing between the walls $x = 0$ and $x = a$ at angle $\theta$ to the axis with $\cos\theta = ck/\omega$, and show that $v_g = c\cos\theta$.
(d) How long does a 10 GHz pulse take to travel 1 m of this guide?

**B7. Coaxial cable impedance.** A coax has inner radius $a$, outer radius $b$, filled with a dielectric of index $n$ ($\mu = \mu_0$).
(a) Generalize (3.34) to the filled cable. Compute the voltage $V_0$ between the conductors and the current $I_0$ on the inner conductor (use Ampère's law around a circle between them).
(b) Show $Z\equiv V_0/I_0 = \frac{Z_0}{2\pi n}\ln\frac ba$ and that the time-averaged power $\int\langle\mathbf S\rangle\cdot d\mathbf a$ equals $\tfrac12V_0I_0$.
(c) Evaluate $Z$ for $b/a = 3.5$ in vacuum and with polyethylene ($\varepsilon_r = 2.25$).

**B8. Anti-reflection coating (requires `project.py`).** A layer of index $n_c$ and thickness $d$ sits on a substrate of index $n_s$, in air.
(a) Let $r_1$, $r_2$ be the normal-incidence amplitude reflection coefficients of the air–layer and layer–substrate interfaces, and $\delta = n_c\omega d/c$. Sum the multiply reflected waves and use the Stokes relations $t_1t_1' = 1 - r_1^2$, $r_1' = -r_1$ (check them from (3.13)) to show the total reflection coefficient is $r = (r_1 + r_2e^{2i\delta})/(1 + r_1r_2e^{2i\delta})$.
(b) Show that $r = 0$ at the vacuum wavelength $\lambda_0$ if $d = \lambda_0/(4n_c)$ and $n_c = \sqrt{n_s}$.
(c) MgF₂ ($n_c = 1.38$) is the usual coating for glass ($n_s = 1.5$). Find $R$ at $\lambda_0$ for a quarter-wave layer and its thickness for $\lambda_0 = 550$ nm.
(d) In `project.py` (Your turn 1), use `fdtd` to send a long pulse with carrier wavelength 120 cells onto a substrate of $n_s = 2.25$, first bare, then with a 20-cell layer of $n_c = 1.5$. Measure the reflected energy fraction in both cases and compare with (b).

## C. Challenge (do at least one)

**C1. Group velocity faster than light.** In a dilute gas with one resonance, write $F = Nq^2f/(m\varepsilon_0)$.
(a) From (3.25), show that at $\omega = \omega_0$, $n = 1$, $dn/d\omega = -F/(\gamma^2\omega_0)$, so $n_g = 1 - F/\gamma^2$; and that $\alpha(\omega_0) = F/(c\gamma)$.
(b) For $\omega_0 = 1$, $\gamma = 0.01$ (units with $c = 1$), find $v_g(\omega_0)$ for $F = 5\times10^{-5}$ and $F = 2\times10^{-4}$, and the intensity absorption length $1/\alpha$ in each case. Confirm both with a few lines of code using the exact $\tilde\varepsilon$ of (3.24).
(c) What does a negative group velocity mean for a pulse entering and leaving a slab? Explain, using (b), why no signal can be sent faster than $c$ this way.

**C2. Fresnel's rhomb.** Beyond the critical angle both $r_s$ and $r_p$ have modulus 1.
(a) Using the right-handed convention ($-r_p$ in our notation), show that the phase difference $\Delta$ between the p and s reflected components satisfies $\tan(\Delta/2) = \cos\theta\sqrt{\sin^2\theta - n^2}/\sin^2\theta$, with $n = n_2/n_1$.
(b) For glass with $n_1 = 1.51$ in air, find the angles at which $\Delta = 45°$ and the maximum $\Delta$.
(c) Explain how two total internal reflections at such an angle turn $45°$ linear light into circular light, and why this "wave plate" works over a much wider range of wavelengths than a birefringent quarter-wave plate.

**C3. The QM-2 wave packet.** A free particle starts as $\psi(x, 0) = (2\pi\sigma_0^2)^{-1/4}e^{-x^2/(4\sigma_0^2)}e^{ik_0x}$, so $|\psi|^2$ has standard deviation $\sigma_0$, and evolves as $\psi(x, t) = \int\phi(k)e^{i(kx - \omega(k)t)}dk$ with $\omega = \hbar k^2/2m$.
(a) Repeat the four steps of (3.30) with $x$ and $t$ in the roles of $t$ and $z$ to show $\sigma(t) = \sigma_0\sqrt{1 + (\hbar t/2m\sigma_0^2)^2}$, and that the centre moves at $\hbar k_0/m$.
(b) How long does it take the width to double for an electron with $\sigma_0 = 1$ Å, and for a 1 g marble with $\sigma_0 = 1\ \mu$m?
(c) Which quantity plays the role of $k''$, and why is the result exact here but only approximate for light in glass?
