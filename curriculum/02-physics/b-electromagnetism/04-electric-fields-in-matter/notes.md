# EM-4: Electric Fields in Matter: Polarization, $\mathbf D$, Dielectrics

**Time:** ~8 hours · **Prerequisites:** EM-3 · **Unlocks:** EM-5 (bound currents and $\mathbf H$ copy bound charges and $\mathbf D$ line by line), EM-7 (the refractive index is $\sqrt{\varepsilon_r(\omega)}$)

> **Conventions in this module.** SI units; $k \equiv 1/4\pi\varepsilon_0$. Field point $\mathbf r$, source point $\mathbf r'$, separation $\boldsymbol{\mathscr r} = \mathbf r - \mathbf r'$ (EM-2). **Polarization** $\mathbf P$ is dipole moment per unit volume, in C/m². **Bound charge**: $\sigma_b = \mathbf P\cdot\hat{\mathbf n}$ with $\hat{\mathbf n}$ the outward normal of the dielectric, $\rho_b = -\nabla\cdot\mathbf P$. **Free charge** $\rho_f$, $\sigma_f$ is every charge that is not bound: charge on conductors, embedded ions, anything you put there. $\mathbf D = \varepsilon_0\mathbf E + \mathbf P$, also in C/m². **Linear media**: $\mathbf P = \varepsilon_0\chi_e\mathbf E$, $\varepsilon = \varepsilon_0(1 + \chi_e)$, $\varepsilon_r = \varepsilon/\varepsilon_0 = 1 + \chi_e$ (the dielectric constant, dimensionless). Atomic polarizability $\alpha$ is defined by $\mathbf p = \alpha\mathbf E$; tables quote the "polarizability volume" $\alpha/4\pi\varepsilon_0$ in Å³. $\mathbf E$ means the **macroscopic** field unless "microscopic" is said. Boundary conditions are "above minus below" with $\hat{\mathbf n}$ pointing from below to above (EM-2 §3.7). Spherical $\theta$ is measured from the polarization or applied-field axis, called $z$.

## 0. Where this sits

You need from EM-2: Gauss's law in both forms (3.4)–(3.5), the potential of a charge distribution (3.8), the boundary conditions at a surface charge (3.10)–(3.12), the field energy (3.16), and the capacitor formulas (3.18)–(3.19). From EM-3: the uniqueness theorems (§3.2), the $\cos\theta$ shell of Example 4.5 (its field is the field of a uniformly polarized sphere), the separation-of-variables solution (3.15), and the dipole in an external field (3.32)–(3.34). From EM-1: the product rules (§3.5) and the divergence theorem (3.18).

This module puts matter between the plates. The payoff is a pair of ideas, bound charge and the auxiliary field $\mathbf D$, that EM-5 copies for magnetism (bound currents, $\mathbf H$) and that EM-6 needs to write Maxwell's equations in matter. The frequency-dependent $\varepsilon_r(\omega)$ previewed in §3.15 is the whole of EM-7's dispersion chapter.

## 1. The problem that forced this idea

In 1837 Faraday built two identical spherical capacitors, charged one, and shared its charge with the other. When the second was empty the charge divided equally. When the second was filled with shellac, sulphur, or glass, it took more than half: the insulator had increased the capacitor's capacity to hold charge at a given potential. Faraday called the ratio the *specific inductive capacity* (our $\varepsilon_r$) and, with Whewell's help, coined the word *dielectric* for an insulator that transmits electric influence through itself. He had no model of why.

The model took forty years. Mossotti (1850) and, independently and more completely, Clausius (1879) supposed that a dielectric is made of small conducting or polarizable spheres in an insulating matrix, each of which acquires a dipole moment in a field, and asked how the measurable constant $\varepsilon_r$ depends on how polarizable each sphere is and how many there are per unit volume. Their answer, with Lorentz's 1880 correction for the field that one molecule actually feels, is the Clausius–Mossotti relation of §3.14. Debye (1912) added molecules with permanent dipole moments and explained why water's $\varepsilon_r$ is 80 while glass's is 5.

Here is the problem that forces the whole apparatus. A cubic metre of glass contains about $10^{29}$ electrons and nuclei. Coulomb's law applied literally, (EM-2 (3.2)) with $N = 10^{29}$ terms, is not merely impractical; it answers the wrong question. The microscopic field inside glass swings between $10^{11}$ V/m next to a nucleus and nearly zero between atoms, over distances of $10^{-10}$ m, and no capacitor measurement sees any of that. What is measured is an average over regions containing millions of atoms. So two things are needed: a definition of the **macroscopic field** as a spatial average, with a proof that the average obeys the same electrostatics as the microscopic field; and a way to separate the charge **you control** (free charge: what you put on the plates) from the charge **the material produces in response** (bound charge: the displaced electron clouds). The free charge is the input. The bound charge is part of the answer, like the induced charge on a conductor in EM-3, and must be solved for.

## 2. Intuition

### A stretched atom is a dipole

Put a neutral atom in a field $\mathbf E$. The nucleus is pulled along $\mathbf E$, the electron cloud against it, until the attraction between them balances the pull. The atom is now a tiny dipole $\mathbf p$ pointing along $\mathbf E$, with $\mathbf p\propto\mathbf E$ for any field you can make in a laboratory (the internal field holding the atom together is $10^{11}$ V/m; you are perturbing it by one part in $10^5$). A polar molecule like water already has a $\mathbf p$; the field exerts a torque that turns it toward $\mathbf E$, against thermal jostling that tries to randomize it. Either way, a block of material in a field becomes a block of aligned dipoles, and $\mathbf P$ is their density.

### Where the charge shows up

Line up dipoles head to tail:

```
   (-)(+) (-)(+) (-)(+) (-)(+) (-)(+)
    ^                              ^
  uncancelled -                uncancelled +
```

Inside the chain every $+$ sits next to a $-$ and they cancel on average. Only the ends are left. A block of uniformly polarized material therefore carries charge only on the two faces that $\mathbf P$ pokes through: negative where $\mathbf P$ enters, positive where it leaves. On a face tilted at angle $\theta$ to $\mathbf P$, fewer chain ends per unit area poke through, by a factor $\cos\theta$. So you should guess $\sigma_b = P\cos\theta = \mathbf P\cdot\hat{\mathbf n}$.

If $\mathbf P$ is *not* uniform the cancellation inside fails: where the chains get stronger from left to right, each $+$ is a little bigger than the $-$ to its right, and net charge is left behind in the volume. More dipole moment leaving a region than entering it means the region has been left negative. Guess: $\rho_b = -\nabla\cdot\mathbf P$. And the total bound charge of any isolated object must be zero, because you only moved charge around; you did not create any.

### What the bound charge does to the field

The bound charges on a uniformly polarized slab are $+\sigma_b$ on one face and $-\sigma_b$ on the other: a parallel-plate capacitor's charges, producing a field $\sigma_b/\varepsilon_0 = P/\varepsilon_0$ inside, pointing from $+$ to $-$, which is *against* $\mathbf P$. A polarized object fights its own polarization. For a sphere the end charges are spread over the curved surface and produce less: a third of the slab's value, $-\mathbf P/3\varepsilon_0$, which is the one number from this module you must carry around. For a long needle polarized along its axis, the end charges are far away and the interior field is nearly zero.

Now the capacitor. Plates with fixed free charge $\pm\sigma_f$ make a field $\sigma_f/\varepsilon_0$. Fill the gap with a dielectric: it polarizes along the field, its faces acquire $\mp\sigma_b$ against the plates, and the net charge the field sees is $\sigma_f - \sigma_b$. The field drops, $V = Ed$ drops, $C = Q/V$ rises. That is Faraday's observation, and $\varepsilon_r = \sigma_f/(\sigma_f - \sigma_b)$.

### Why anyone would invent $\mathbf D$

Gauss's law counts all charge: $\varepsilon_0\nabla\cdot\mathbf E = \rho_f + \rho_b = \rho_f - \nabla\cdot\mathbf P$. Move the $\mathbf P$ to the left: $\nabla\cdot(\varepsilon_0\mathbf E + \mathbf P) = \rho_f$. The combination $\mathbf D = \varepsilon_0\mathbf E + \mathbf P$ has a divergence that sees only the charge you put there. It is a **bookkeeping field**: it lets you use Gauss's law with free charge alone, without first solving for the bound charge you do not yet know. In a linear material $\mathbf P\propto\mathbf E$, so $\mathbf D\propto\mathbf E$ and knowing $\mathbf D$ gives $\mathbf E$ at once. The guess for the module's main theorems: $\nabla\cdot\mathbf D = \rho_f$; at a surface the normal component of $\mathbf D$ jumps by $\sigma_f$ only; a dielectric with no free charge inside it has no bound volume charge either.

## 3. Formal development

### 3.1 Induced dipoles and the atomic polarizability

**Definition 3.1.** For an atom or nonpolar molecule in a field $\mathbf E$ weak compared with the internal atomic fields, the induced dipole moment is
$$\mathbf p = \alpha\mathbf E, \tag{3.1}$$
and $\alpha$ is the **atomic polarizability** (units C m²/V $=$ F m²). Tables quote $\alpha/4\pi\varepsilon_0$, which has units of volume: hydrogen $0.667$ Å³, helium $0.205$, argon $1.64$, N$_2$ $1.74$, water $1.45$ (electronic part only). For molecules of low symmetry $\alpha$ is a symmetric tensor, $p_i = \alpha_{ij}E_j$ (LA-12); CO$_2$ is nearly twice as polarizable along its axis as across it.

**A crude model that gets the size right.** Model the atom as a point nucleus $+q$ at the centre of a uniformly charged ball of radius $a$ with total charge $-q$ (the electron cloud).

▶ **Try it first.** Displace the nucleus a distance $d$ from the cloud's centre. What field does the cloud exert on it (EM-2 Example 4.1)? Set that equal to the external field.

1. The field of the cloud at distance $d\le a$ from its centre is, by EM-2 Example 4.1 with $Q = -q$, $R = a$, $E_{\text{cloud}} = \dfrac{qd}{4\pi\varepsilon_0a^3}$ directed toward the centre: a linear restoring field.
2. Equilibrium: the external field pulls the nucleus with $qE$, the cloud pulls it back with $q\cdot qd/(4\pi\varepsilon_0a^3)$. So $d = 4\pi\varepsilon_0a^3E/q$.
3. The dipole moment is $p = qd$ (nucleus displaced by $d$ relative to the cloud's centre of charge):
   $$\alpha = 4\pi\varepsilon_0a^3,\qquad\frac{\alpha}{4\pi\varepsilon_0} = a^3. \tag{3.2}$$
∎

■ **What just happened.** The polarizability volume is the atom's volume, roughly. For hydrogen the measured $0.667$ Å³ corresponds to $a = 0.874$ Å (Exercise A1), not far from the Bohr radius $0.529$ Å; the exact quantum answer is $\tfrac92a_0^3 = 0.667$ Å³ (QM-8). The restoring field in step 1 is linear in $d$, which is *why* $\mathbf p\propto\mathbf E$: a harmonic spring. That spring reappears as the Lorentz oscillator in §3.15.

### 3.2 Polar molecules: torque, alignment, force

A molecule with a permanent moment $\mathbf p$ in a field $\mathbf E$ has, from EM-3 (3.32)–(3.34),
$$U = -\mathbf p\cdot\mathbf E,\qquad\mathbf N = \mathbf p\times\mathbf E,\qquad\mathbf F = (\mathbf p\cdot\nabla)\mathbf E. \tag{3.3}$$
The torque turns $\mathbf p$ toward $\mathbf E$; the force is nonzero only in a *non-uniform* field. In a gas or liquid the torque competes with thermal agitation, and the energy difference between aligned and anti-aligned, $2pE$, is tiny compared with $k_BT$ for any laboratory field: for water in $10^5$ V/m, $2pE/k_BT = 3\times10^{-4}$ (Exercise A3). The thermal average of the component along $\mathbf E$ is $\langle p_z\rangle \approx p^2E/3k_BT$ (SM-3 derives this Langevin average; EM-3 Exercise A6 did the integral), so a polar gas behaves like a nonpolar one with an *orientational* polarizability $p^2/3k_BT$ added to its electronic $\alpha$. This is Debye's term, and it is why $\varepsilon_r$ of polar substances falls with temperature.

For later use, the field of a point dipole $\mathbf p$ at the origin (EM-2 (4.3)):
$$\mathbf E_{\text{dip}}(\mathbf r) = \frac{1}{4\pi\varepsilon_0r^3}\left[3(\mathbf p\cdot\hat{\mathbf r})\hat{\mathbf r} - \mathbf p\right],\qquad V_{\text{dip}} = \frac{\mathbf p\cdot\hat{\mathbf r}}{4\pi\varepsilon_0r^2}. \tag{3.4}$$

### 3.3 Polarization

**Definition 3.2.** The **polarization** $\mathbf P(\mathbf r)$ is the dipole moment per unit volume: a volume $d\tau'$ at $\mathbf r'$ carries dipole moment $\mathbf P(\mathbf r')\,d\tau'$. If there are $N$ molecules per unit volume with average moment $\langle\mathbf p\rangle$,
$$\mathbf P = N\langle\mathbf p\rangle. \tag{3.5}$$

$\mathbf P$ is a macroscopic field: it is defined by averaging over a region containing many molecules, exactly as $\rho$ in EM-2 was. How it depends on $\mathbf E$ is a property of the material (§3.10). For now $\mathbf P$ is given, and the question is what field it produces.

### 3.4 The field of a polarized object: bound charges

This is the module's central result.

▶ **Try it first.** Write the potential at $\mathbf r$ of the dipole $\mathbf P\,d\tau'$ at $\mathbf r'$ using (3.4), integrate over the object, then try to turn the integrand into something that looks like $\rho/\mathscr r$ or $\sigma/\mathscr r$ by integrating by parts.

**Theorem 3.3 (Bound charges).** The potential of a polarized object is
$$V(\mathbf r) = \frac{1}{4\pi\varepsilon_0}\int_{\mathcal V}\frac{\mathbf P(\mathbf r')\cdot\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\,d\tau' \tag{3.6}$$
and equals the potential of a surface charge plus a volume charge,
$$V(\mathbf r) = \frac{1}{4\pi\varepsilon_0}\oint_{\mathcal S}\frac{\sigma_b}{\mathscr r}\,da' + \frac{1}{4\pi\varepsilon_0}\int_{\mathcal V}\frac{\rho_b}{\mathscr r}\,d\tau', \tag{3.7}$$
with
$$\sigma_b = \mathbf P\cdot\hat{\mathbf n},\qquad\rho_b = -\nabla\cdot\mathbf P. \tag{3.8}$$

**Derivation.**

1. The potential of the single dipole $\mathbf P(\mathbf r')\,d\tau'$ is, by (3.4) with $\mathbf r$ replaced by the separation $\boldsymbol{\mathscr r} = \mathbf r - \mathbf r'$, $dV = \dfrac{1}{4\pi\varepsilon_0}\dfrac{\mathbf P(\mathbf r')\cdot\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\,d\tau'$. Superpose: (3.6).
2. **The key identity.** Differentiating $1/\mathscr r$ with respect to the *source* coordinates,
   $$\nabla'\left(\frac1{\mathscr r}\right) = \frac{\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}. \tag{3.9}$$
   Check one component: $\partial_{x'}\left[(x-x')^2 + (y-y')^2 + (z-z')^2\right]^{-1/2} = (x - x')/\mathscr r^3$, which is the $x$ component of $\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$. The sign is opposite to the field-point gradient $\nabla(1/\mathscr r) = -\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$ of EM-2 §3.5, because $\mathscr r$ grows when $\mathbf r$ moves away and shrinks when $\mathbf r'$ moves toward $\mathbf r$. (Project Part 1 checks (3.9) with SymPy.)
3. Substitute: $V = \dfrac{1}{4\pi\varepsilon_0}\displaystyle\int\mathbf P\cdot\nabla'\left(\frac1{\mathscr r}\right)d\tau'$.
4. Product rule (EM-1 §3.5) for $\nabla'\cdot(f\mathbf A)$ with $f = 1/\mathscr r$, $\mathbf A = \mathbf P$:
   $$\nabla'\cdot\left(\frac{\mathbf P}{\mathscr r}\right) = \frac1{\mathscr r}\nabla'\cdot\mathbf P + \mathbf P\cdot\nabla'\left(\frac1{\mathscr r}\right),\quad\text{so}\quad\mathbf P\cdot\nabla'\left(\frac1{\mathscr r}\right) = \nabla'\cdot\left(\frac{\mathbf P}{\mathscr r}\right) - \frac{\nabla'\cdot\mathbf P}{\mathscr r}.$$
5. Integrate the first term with the divergence theorem (EM-1 (3.18)) over the object's volume $\mathcal V$ with boundary $\mathcal S$ and outward normal $\hat{\mathbf n}$:
   $$V = \frac{1}{4\pi\varepsilon_0}\oint_{\mathcal S}\frac{\mathbf P\cdot\hat{\mathbf n}}{\mathscr r}\,da' + \frac{1}{4\pi\varepsilon_0}\int_{\mathcal V}\frac{-\nabla'\cdot\mathbf P}{\mathscr r}\,d\tau'.$$
6. Compare with EM-2 (3.8), the potential of a surface charge $\sigma$ and a volume charge $\rho$: the first term is the potential of $\sigma_b = \mathbf P\cdot\hat{\mathbf n}$ on the surface, the second of $\rho_b = -\nabla\cdot\mathbf P$ in the volume. ∎

■ **What just happened.** One integration by parts turned "a continuum of dipoles" into "a charge distribution." Nothing approximate happened, except in step 1, where we used the point-dipole potential (3.4) for a volume element; that is exact for field points far from $d\tau'$ compared with molecular size and is the subject of §3.7 for field points inside the material. From here on, a polarized object *is* the charge distribution (3.8), and every tool of EM-2 and EM-3 applies to it.

**Corollary 3.4 (Total bound charge is zero).** For any bounded object,
$$Q_b = \int_{\mathcal V}\rho_b\,d\tau + \oint_{\mathcal S}\sigma_b\,da = -\int_{\mathcal V}\nabla\cdot\mathbf P\,d\tau + \oint_{\mathcal S}\mathbf P\cdot d\mathbf a = 0 \tag{3.10}$$
by the divergence theorem. ∎ Polarizing an object rearranges charge; it creates none. Use (3.10) as a check on every bound-charge calculation (Example 4.2, Exercises B1, B3, B5; project Part 1 checks it on a grid for an irregular $\mathbf P$).

If $\mathbf P$ jumps across an *internal* surface, (3.8) applies on each side with its own outward normal: $\sigma_b = \mathbf P_1\cdot\hat{\mathbf n}_1 + \mathbf P_2\cdot\hat{\mathbf n}_2 = (\mathbf P_1 - \mathbf P_2)\cdot\hat{\mathbf n}_{12}$, the discontinuity of the normal component (Exercise B5).

### 3.5 Bound charge is real charge: counting at a tilted surface

The derivation above is algebra. Here is the same result by counting electrons, so that you believe $\sigma_b$ is charge you could measure.

▶ **Try it first.** A material has $N$ molecules per unit volume, each a pair $\pm q$ separated by $\mathbf d$, so $\mathbf P = Nq\mathbf d$. Cut a surface through it with normal $\hat{\mathbf n}$. Which molecules have their $+q$ outside the surface and their $-q$ inside?

1. Think of polarizing as taking a neutral medium and displacing all positive charge by $\mathbf d$ relative to the negative charge. Consider a patch $da$ of the surface. The positive charges that end up *outside* the patch, having started inside, are those that began within the slanted prism of base $da$ and height $\mathbf d\cdot\hat{\mathbf n} = d\cos\theta$ behind the patch, volume $d\cos\theta\,da$.
2. That prism contains $N\,d\cos\theta\,da$ molecules, so the charge pushed out through $da$ is $Nqd\cos\theta\,da = P\cos\theta\,da = \mathbf P\cdot\hat{\mathbf n}\,da$. If the surface is the object's boundary, this charge is stuck there: $\sigma_b = \mathbf P\cdot\hat{\mathbf n}$. If $\theta > 90°$, $\cos\theta<0$ and negative charge is left exposed.
3. For a closed surface inside the material, the net charge pushed out of the enclosed volume is $\oint\mathbf P\cdot d\mathbf a$. The volume was neutral before, so it is now left with $Q_{\text{inside}} = -\oint\mathbf P\cdot d\mathbf a = -\int\nabla\cdot\mathbf P\,d\tau$. Since this holds for every volume, $\rho_b = -\nabla\cdot\mathbf P$. ∎

■ **What just happened.** $\sigma_b$ and $\rho_b$ are not mathematical fictions produced by integration by parts; they are the electrons that were pushed across the surface, and (3.5)'s density times (3.1)'s displacement gives the same formula as the vector calculus did.

### 3.6 The uniformly polarized sphere

A sphere of radius $R$ with uniform $\mathbf P = P\hat{\mathbf z}$. By (3.8), $\rho_b = 0$ and $\sigma_b = \mathbf P\cdot\hat{\mathbf r} = P\cos\theta$.

▶ **Try it first.** This surface charge is EM-3 Example 4.5 with $\sigma_0 = P$. Write down the fields.

**Route 1 (EM-3's $\cos\theta$ shell).** EM-3 (4.5) with $\sigma_0 = P$:
$$V_{\text{in}} = \frac{P}{3\varepsilon_0}\,r\cos\theta = \frac{Pz}{3\varepsilon_0},\qquad\mathbf E_{\text{in}} = -\frac{\mathbf P}{3\varepsilon_0};\qquad V_{\text{out}} = \frac{PR^3}{3\varepsilon_0}\frac{\cos\theta}{r^2} = \frac{\mathbf p\cdot\hat{\mathbf r}}{4\pi\varepsilon_0r^2},\quad\mathbf p = \tfrac43\pi R^3\mathbf P. \tag{3.11}$$
Inside, a **uniform** field, antiparallel to $\mathbf P$, one third of the slab's $P/\varepsilon_0$. Outside, a **pure dipole** with the sphere's total moment, exactly, not just to leading order.

**Route 2 (two offset balls), which also explains the $\tfrac13$.** Superpose a ball of uniform $+\rho$ centred at $+\mathbf d/2$ and a ball of $-\rho$ centred at $-\mathbf d/2$, with $d\ll R$. In the overlap the charges cancel; a thin crescent of thickness $d\cos\theta$ at the surface is left with one sign, which is the $\sigma_b = \rho d\cos\theta$ of §3.5 with $P = \rho d$.

1. Inside a uniform ball the field is $\rho(\mathbf r - \mathbf r_c)/3\varepsilon_0$ (EM-2 Example 4.1). At a point inside both balls,
   $$\mathbf E = \frac{\rho}{3\varepsilon_0}\left(\mathbf r - \frac{\mathbf d}2\right) - \frac{\rho}{3\varepsilon_0}\left(\mathbf r + \frac{\mathbf d}2\right) = -\frac{\rho\mathbf d}{3\varepsilon_0} = -\frac{\mathbf P}{3\varepsilon_0}.$$
   Uniform, with no approximation at all.
2. Outside, each ball is a point charge $\pm Q = \pm\tfrac43\pi R^3\rho$ at $\pm\mathbf d/2$: a physical dipole with $\mathbf p = Q\mathbf d = \tfrac43\pi R^3\mathbf P$, which becomes the pure dipole field as $d\to0$ at fixed $\rho d$. ∎

So
$$\boxed{\ \mathbf E_{\text{in}} = -\frac{\mathbf P}{3\varepsilon_0}\ }\qquad\text{(uniformly polarized sphere)}. \tag{3.12}$$

■ **What just happened.** The $\tfrac13$ is the $\tfrac13$ in the interior field $\rho\mathbf r/3\varepsilon_0$ of a uniform ball, which is the $3$ in "$\nabla\cdot\mathbf r = 3$": the field of a ball spreads its flux over three dimensions. A slab spreads over one and gives $P/\varepsilon_0$; a long cylinder polarized across its axis spreads over two and gives $P/2\varepsilon_0$ (Exercise B7). The project computes (3.12) three ways: Coulomb sum over surface patches ($10^{-10}$), two offset balls ($10^{-9}$), and a lattice of $2109$ point dipoles, where the result is different and instructive (§3.7).

### 3.7 What "the field inside a dielectric" means

Theorem 3.3 used the point-dipole potential for each $d\tau'$. For a field point inside the material that is wrong for the molecules nearby: a few ångströms from a molecule its field is not a dipole field, and the microscopic field varies by orders of magnitude over atomic distances. Yet (3.7) is used for interior points without apology. Here is why it is right.

**Definition 3.5.** The **macroscopic field** at $\mathbf r$ is the average of the microscopic field over a sphere centred at $\mathbf r$ of radius $R_s$ large compared with molecular spacing (so it contains thousands of molecules) and small compared with any distance over which $\mathbf P$ varies.

**Lemma 3.6 (Averaging lemma).** Let $\mathbf E_{\text{avg}}$ be the average over a ball of radius $R_s$ of the field produced by a set of charges.
(a) For charges **outside** the ball, $\mathbf E_{\text{avg}}$ equals the field they produce at the ball's centre.
(b) For charges **inside** the ball, with total dipole moment $\mathbf p$ about the centre,
$$\mathbf E_{\text{avg}} = -\frac{\mathbf p}{4\pi\varepsilon_0R_s^3}. \tag{3.13}$$

▶ **Try it first.** For (b): the average over the ball of the field of a point charge $q$ at $\mathbf r'$ is $\frac1{\mathcal V}\int_{\text{ball}}\frac{q}{4\pi\varepsilon_0}\frac{\mathbf r - \mathbf r'}{|\mathbf r - \mathbf r'|^3}d\tau$. Read that integral as the field *at* $\mathbf r'$ produced by a uniformly charged ball.

**Derivation.**

1. (a) Each component of the field of charges outside the ball is harmonic inside it (EM-2 (3.9) with $\rho = 0$, applied to $\partial V/\partial x_i$). The mean-value property EM-3 (3.1) gives the centre value as the average over every spherical *surface* inside the ball; averaging those over radius gives the centre value as the *volume* average. ∎(a)
2. (b) For one charge $q$ at $\mathbf r'$ inside, with $\mathcal V = \tfrac43\pi R_s^3$:
   $$\mathbf E_{\text{avg}} = \frac1{\mathcal V}\int_{\text{ball}}\frac{q}{4\pi\varepsilon_0}\frac{\mathbf r - \mathbf r'}{|\mathbf r - \mathbf r'|^3}\,d\tau = -\int_{\text{ball}}\frac{(q/\mathcal V)}{4\pi\varepsilon_0}\frac{\mathbf r' - \mathbf r}{|\mathbf r' - \mathbf r|^3}\,d\tau.$$
   The right side is minus the field at $\mathbf r'$ of a ball of uniform density $\rho_u = q/\mathcal V$ (the integration variable $\mathbf r$ is now the source point).
3. By EM-2 Example 4.1, that field at an interior point is $\rho_u\mathbf r'/3\varepsilon_0 = q\mathbf r'/(4\pi\varepsilon_0R_s^3)$. So $\mathbf E_{\text{avg}} = -q\mathbf r'/(4\pi\varepsilon_0R_s^3)$, and $q\mathbf r'$ is the dipole moment of $q$ about the centre. Sum over all charges inside: (3.13). ∎
4. (Consistency.) If $\mathbf r'$ is *outside*, step 3 gives instead $\rho_u\mathcal V\hat{\mathbf r}'/(4\pi\varepsilon_0r'^2)$, so $\mathbf E_{\text{avg}} = -q\hat{\mathbf r}'/(4\pi\varepsilon_0r'^2)$, which is the field at the centre due to $q$: (a) again. The project checks both cases by quadrature.

■ **What just happened.** "Average of the field over the ball due to $q$" and "field at $q$ due to a uniformly charged ball" are the same integral read two ways, and the second is a problem you solved in EM-2.

**Theorem 3.7 (The macroscopic field is what (3.7) computes).** Now compute the macroscopic field at a point inside a dielectric. Split the material into the molecules outside the averaging sphere and those inside.

1. *Outside.* These molecules are at least $R_s$ away, many molecular diameters, so the dipole approximation (3.4) is accurate for each, and by Lemma 3.6(a) their average field equals the field they produce at the centre: exactly (3.6) integrated over the material minus the sphere.
2. *Inside.* The sphere contains dipole moment $\mathbf p = \tfrac43\pi R_s^3\mathbf P$ (since $\mathbf P$ is uniform over so small a region). By (3.13) their average field is $-\mathbf P/3\varepsilon_0$, *whatever the detailed microscopic arrangement*.
3. But $-\mathbf P/3\varepsilon_0$ is exactly the field that (3.6), applied to a uniformly polarized sphere, produces at interior points: that is (3.12). So the inside contribution to the macroscopic field is also what (3.6) gives.
4. Adding: the macroscopic field is (3.6) integrated over the whole material, which by Theorem 3.3 is the field of $\sigma_b$ and $\rho_b$. ∎

■ **What just happened.** The point-dipole approximation is wrong for the nearby molecules, but the *average* of their true field equals the *average* of their dipole-approximated field, both being $-\mathbf P/3\varepsilon_0$. Hence (3.7) and everything built on it is correct for the macroscopic field. The project's lattice makes the distinction vivid: $2109$ point dipoles on a cubic lattice filling a sphere give $|\mathbf E| \le 0.006\,P/\varepsilon_0$ *at the lattice sites*, not $0.333\,P/\varepsilon_0$; the macroscopic value is an average over the cell, dominated by the huge fields near each dipole, not the value at any particular point. The value at a site is what a molecule actually feels, and §3.14 needs it.

### 3.8 Gauss's law in dielectrics

▶ **Try it first.** Put $\rho = \rho_f + \rho_b$ into EM-2 (3.4) and move the $\mathbf P$ term across.

1. The total charge density is $\rho = \rho_f + \rho_b = \rho_f - \nabla\cdot\mathbf P$, so $\varepsilon_0\nabla\cdot\mathbf E = \rho_f - \nabla\cdot\mathbf P$.
2. Collect divergences: $\nabla\cdot(\varepsilon_0\mathbf E + \mathbf P) = \rho_f$. Define the **electric displacement** $\mathbf D \equiv \varepsilon_0\mathbf E + \mathbf P$:
   $$\nabla\cdot\mathbf D = \rho_f. \tag{3.14}$$
3. Integrate over a volume and apply the divergence theorem:
   $$\oint_{\mathcal S}\mathbf D\cdot d\mathbf a = Q_{f,\text{enc}}. \tag{3.15}$$
∎

■ **What just happened.** Gauss's law for $\mathbf D$ has only the free charge on the right. In any geometry where symmetry fixes the direction of $\mathbf D$ (spherical, cylindrical, planar), (3.15) gives $\mathbf D$ from the free charge alone, before you know anything about the material's response. Then, for a linear material, $\mathbf E = \mathbf D/\varepsilon$ (§3.10).

**What $\mathbf D$ is not.** $\mathbf D$ is not "the electric field in the material." $\mathbf E$ is. A test charge in a cavity, a force on a bound charge, the energy density: all are set by $\mathbf E$. $\mathbf D$ is a bookkeeping combination that happens to have a convenient divergence. Specifically:

- **Its curl need not vanish.** $\nabla\times\mathbf D = \varepsilon_0\nabla\times\mathbf E + \nabla\times\mathbf P = \nabla\times\mathbf P$, and $\nabla\times\mathbf P\ne0$ at the edges of a polarized object in general. So $\mathbf D$ has no scalar potential, and
- **there is no Coulomb's law for $\mathbf D$.** You cannot compute $\mathbf D$ from $\rho_f$ by $\int\rho_f\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$. A bar electret (a permanently polarized cylinder) has $\rho_f = 0$ everywhere, yet $\mathbf D\ne0$: $\mathbf D = \varepsilon_0\mathbf E + \mathbf P$ with $\mathbf E$ the field of the bound end charges. Only $\nabla\cdot\mathbf D = 0$ is true, and a divergence-free field is far from unique (EM-1 Helmholtz, §3.12). $\mathbf D$ is determined by $\rho_f$ alone only when symmetry kills the curl part.
- **The parallel with $\mathbf H$.** In EM-5, magnetized matter has bound currents $\mathbf J_b = \nabla\times\mathbf M$, $\mathbf K_b = \mathbf M\times\hat{\mathbf n}$, and $\mathbf H = \mathbf B/\mu_0 - \mathbf M$ obeys $\nabla\times\mathbf H = \mathbf J_f$ with free current only, while $\nabla\cdot\mathbf H = -\nabla\cdot\mathbf M$ need not vanish. Every sentence of this subsection has a twin there.

### 3.9 Boundary conditions

▶ **Try it first.** Repeat EM-2 §3.7 with (3.15) for the pillbox. Then ask what, if anything, changed for the loop.

**Theorem 3.8.** At a surface carrying free charge $\sigma_f$, with $\hat{\mathbf n}$ from "below" to "above",
$$D^\perp_{\text{above}} - D^\perp_{\text{below}} = \sigma_f, \tag{3.16}$$
$$\mathbf E^\parallel_{\text{above}} = \mathbf E^\parallel_{\text{below}},\qquad\text{equivalently}\qquad\mathbf D^\parallel_{\text{above}} - \mathbf D^\parallel_{\text{below}} = \mathbf P^\parallel_{\text{above}} - \mathbf P^\parallel_{\text{below}}. \tag{3.17}$$

*Proof.* A pillbox straddling the surface, faces of area $A$, thickness $\to0$: by (3.15) the flux is $(D^\perp_{\text{above}} - D^\perp_{\text{below}})A$ and the enclosed free charge is $\sigma_fA$, giving (3.16). The loop argument of EM-2 §3.7 used only $\nabla\times\mathbf E = 0$, which is unchanged by matter (the bound charges are static charges like any other), so $\mathbf E^\parallel$ is continuous; the $\mathbf D^\parallel$ statement follows from $\mathbf D = \varepsilon_0\mathbf E + \mathbf P$. ∎

■ **What just happened.** The normal component of $\mathbf D$ is continuous across a dielectric interface with no free charge on it, while $E^\perp$ jumps by $\sigma_b/\varepsilon_0$; the tangential component of $\mathbf E$ is continuous, while $\mathbf D^\parallel$ jumps. Field lines refract at an interface: Exercise A5 derives $\tan\theta_2/\tan\theta_1 = \varepsilon_2/\varepsilon_1$. The project's half-filled capacitor plots $D_y$ (flat) against $E_y$ (a step) through the interface.

### 3.10 Linear dielectrics

**Definition 3.9.** A **linear dielectric** is one in which, for the fields of interest,
$$\mathbf P = \varepsilon_0\chi_e\mathbf E,\qquad\text{so}\qquad\mathbf D = \varepsilon_0(1 + \chi_e)\mathbf E = \varepsilon\mathbf E,\qquad\varepsilon_r \equiv \frac{\varepsilon}{\varepsilon_0} = 1 + \chi_e. \tag{3.18}$$
$\chi_e$ is the **electric susceptibility** (dimensionless), $\varepsilon$ the **permittivity**, $\varepsilon_r$ the **dielectric constant**. $\mathbf E$ here is the *total* macroscopic field, including the part produced by the material's own bound charge; (3.18) is therefore an implicit equation, which is why "dielectric sphere in a field" (§3.11) needs a self-consistency argument. For an anisotropic crystal $\chi_e$ is a symmetric tensor and $\mathbf P\nparallel\mathbf E$ (LA-12; EM-7 birefringence).

| Material (20 °C, 1 atm unless noted) | $\varepsilon_r$ | Material | $\varepsilon_r$ |
|---|---|---|---|
| vacuum | 1 (exactly) | polyethylene, Teflon | 2.3, 2.1 |
| helium | 1.000065 | benzene | 2.28 |
| hydrogen | 1.00025 | fused quartz, glass | 3.8, 5–10 |
| dry air | 1.00054 | mica | 5.4 |
| nitrogen | 1.00055 | diamond, NaCl | 5.7, 5.9 |
| water vapour (100 °C) | 1.00587 | silicon | 11.8 |
| liquid argon (87 K) | 1.5 | methanol | 33 |
| | | water | 80.1 |
| | | ice (−30 °C) | 99 |
| | | BaTiO$_3$ (ceramic, ferroelectric) | $10^3$–$10^4$ |

Gases: $\chi_e\sim10^{-4}$–$10^{-3}$, proportional to density. Nonpolar liquids and solids: $2$–$12$. Polar liquids: tens, falling with temperature. Ferroelectrics: thousands, nonlinear, hysteretic (§3.15). The first three columns are from Griffiths Table 4.2 (CRC Handbook).

**Theorem 3.10 (Bound charge in a homogeneous linear dielectric).** Where $\chi_e$ is uniform,
$$\rho_b = -\frac{\chi_e}{1 + \chi_e}\,\rho_f. \tag{3.19}$$
In particular, a homogeneous linear dielectric containing no free charge has no bound volume charge: all its bound charge sits on its surfaces.

*Proof.* From (3.18), $\mathbf P = \varepsilon_0\chi_e\mathbf E = \varepsilon_0\chi_e\mathbf D/\varepsilon = \dfrac{\chi_e}{1 + \chi_e}\mathbf D$. Take the divergence, with $\chi_e$ constant, and use (3.14): $\rho_b = -\nabla\cdot\mathbf P = -\dfrac{\chi_e}{1+\chi_e}\nabla\cdot\mathbf D = -\dfrac{\chi_e}{1 + \chi_e}\rho_f$. ∎

■ **What just happened.** Bound charge in the bulk appears only where there is free charge in the bulk, and it partly cancels it: a free point charge $q$ embedded in a dielectric is surrounded by bound charge $-\dfrac{\chi_e}{\varepsilon_r}q$ (apply (3.19) to a delta function), leaving net charge $q/\varepsilon_r$ at the centre. That is why the field of a point charge inside a dielectric is
$$\mathbf E = \frac{1}{4\pi\varepsilon}\frac{q}{r^2}\hat{\mathbf r} = \frac{1}{\varepsilon_r}\,\frac{q}{4\pi\varepsilon_0r^2}\hat{\mathbf r}: \tag{3.20}$$
Gauss for $\mathbf D$ gives $D = q/4\pi r^2$, and $\mathbf E = \mathbf D/\varepsilon$. The material screens the charge by the factor $\varepsilon_r$ (Exercise B4 does this inside a finite sphere, where the missing bound charge reappears on the surface).

**Capacitors.** For a capacitor completely filled with a linear dielectric, (3.15) gives the same $\mathbf D$ as the vacuum $\mathbf E$ times $\varepsilon_0$ for the same free charge, so $\mathbf E$ is reduced by $\varepsilon_r$ everywhere, $V$ is reduced by $\varepsilon_r$, and
$$C = \varepsilon_rC_{\text{vac}}. \tag{3.21}$$
This is Faraday's measurement. Partial filling is Example 4.3.

### 3.11 A dielectric sphere in a uniform field

A sphere of linear dielectric $\varepsilon_r$, radius $R$, is placed in a uniform field $\mathbf E_0 = E_0\hat{\mathbf z}$. Find the field inside.

▶ **Try it first.** Guess that $\mathbf P$ is uniform. Then the sphere's own field inside is (3.12). Write the total interior field and demand that it produce exactly that $\mathbf P$ through (3.18).

**Route 1: self-consistency.**

1. *Guess* a uniform $\mathbf P$. Its field inside the sphere is $-\mathbf P/3\varepsilon_0$ by (3.12), so the total interior field is $\mathbf E_{\text{in}} = \mathbf E_0 - \mathbf P/3\varepsilon_0$, which is uniform. A uniform $\mathbf E_{\text{in}}$ produces a uniform $\mathbf P = \varepsilon_0\chi_e\mathbf E_{\text{in}}$: the guess is consistent with itself.
2. Solve: $\mathbf P = \varepsilon_0\chi_e(\mathbf E_0 - \mathbf P/3\varepsilon_0)\Rightarrow\mathbf P\left(1 + \dfrac{\chi_e}3\right) = \varepsilon_0\chi_e\mathbf E_0\Rightarrow\mathbf P = \dfrac{3\varepsilon_0\chi_e}{3 + \chi_e}\mathbf E_0 = 3\varepsilon_0\dfrac{\varepsilon_r - 1}{\varepsilon_r + 2}\mathbf E_0$.
3. Then $\mathbf E_{\text{in}} = \mathbf E_0 - \dfrac{\mathbf P}{3\varepsilon_0} = \left(1 - \dfrac{\varepsilon_r - 1}{\varepsilon_r + 2}\right)\mathbf E_0$:
   $$\mathbf E_{\text{in}} = \frac{3}{\varepsilon_r + 2}\,\mathbf E_0,\qquad\mathbf P = 3\varepsilon_0\frac{\varepsilon_r - 1}{\varepsilon_r + 2}\mathbf E_0,\qquad\mathbf p = 4\pi\varepsilon_0R^3\frac{\varepsilon_r - 1}{\varepsilon_r + 2}\mathbf E_0. \tag{3.22}$$
4. *Why the guess is the answer.* Outside, the field is $\mathbf E_0$ plus the dipole field of $\mathbf p$ (3.11); inside it is (3.22). This pair satisfies Laplace's equation on both sides, tends to $\mathbf E_0$ at infinity, and (check) has continuous $V$ and continuous $D^\perp$ at $r = R$. By the uniqueness theorem (EM-3 §3.2, extended to dielectrics in Exercise C3 of EM-3's style; the proof is the same Green's-identity argument with $\varepsilon|\nabla V|^2$), it is *the* solution. ∎

**Route 2: separation of variables (the check).** By EM-3 (3.15), with only the $\ell = 1$ terms surviving (the applied field is $\ell = 1$, and for every other $\ell$ the two matching conditions below are homogeneous with nonzero determinant, forcing zero coefficients):
$$V_{\text{in}} = -A\,r\cos\theta,\qquad V_{\text{out}} = -E_0\,r\cos\theta + \frac{B\cos\theta}{r^2}.$$
Continuity of $V$ at $R$: $-AR = -E_0R + B/R^2$. Continuity of $D^\perp$ (3.16) with no free charge: $\varepsilon_r\,\partial V_{\text{in}}/\partial r = \partial V_{\text{out}}/\partial r$, i.e. $\varepsilon_rA = E_0 + 2B/R^3$. From the first, $B = R^3(E_0 - A)$; substitute: $\varepsilon_rA = E_0 + 2E_0 - 2A$, so
$$A = \frac{3E_0}{\varepsilon_r + 2},\qquad B = R^3E_0\frac{\varepsilon_r - 1}{\varepsilon_r + 2},$$
and $\mathbf E_{\text{in}} = A\hat{\mathbf z}$ is (3.22). The project solves the pair with SymPy. ∎

■ **What just happened.** Route 1 looks circular: it assumes the field is uniform to prove the field is uniform. It is not circular; it is a fixed-point argument whose legitimacy rests on uniqueness. Route 2 does not assume anything and gets the same answer; that agreement is what should convince you. Limits: $\varepsilon_r\to1$ gives $\mathbf E_{\text{in}}\to\mathbf E_0$; $\varepsilon_r\to\infty$ gives $\mathbf E_{\text{in}}\to0$ and $\mathbf p\to4\pi\varepsilon_0R^3\mathbf E_0$, the conducting sphere of EM-3 (4.3). The field is *reduced* inside (the bound charges on the surface oppose $\mathbf E_0$) but by less than $\varepsilon_r$: for $\varepsilon_r = 4$, $E_{\text{in}} = E_0/2$, not $E_0/4$. The same argument for a sphere of $\varepsilon_{\text{in}}$ in a medium $\varepsilon_{\text{out}}$ gives $\mathbf E_{\text{in}} = 3\varepsilon_{\text{out}}\mathbf E_0/(\varepsilon_{\text{in}} + 2\varepsilon_{\text{out}})$, and for a cylinder across its axis $2\varepsilon_{\text{out}}\mathbf E_0/(\varepsilon_{\text{in}} + \varepsilon_{\text{out}})$ (Exercise B7; project Part 4 solves it on a grid).

### 3.12 Energy in dielectric systems

▶ **Try it first.** Build the free charge up slowly from zero, $\delta\rho_f$ at a time. The work for each increment is $\int V\delta\rho_f\,d\tau$ (EM-2 (3.13)). Replace $\delta\rho_f$ by $\nabla\cdot\delta\mathbf D$ and integrate by parts.

**Theorem 3.11.** For a system of linear dielectrics and free charges,
$$W = \frac12\int\mathbf D\cdot\mathbf E\,d\tau. \tag{3.23}$$

**Derivation.**

1. Increase the free charge by $\delta\rho_f(\mathbf r)$ everywhere while the dielectric responds. The work done on the free charge is $\delta W = \int V\,\delta\rho_f\,d\tau$, with $V$ the current potential.
2. By (3.14), $\delta\rho_f = \nabla\cdot(\delta\mathbf D)$. Product rule: $V\nabla\cdot(\delta\mathbf D) = \nabla\cdot(V\delta\mathbf D) - \delta\mathbf D\cdot\nabla V = \nabla\cdot(V\delta\mathbf D) + \delta\mathbf D\cdot\mathbf E$.
3. The divergence term becomes a surface integral at infinity, where $V\sim1/r$ and $\delta D\sim1/r^2$, so it vanishes. Hence $\delta W = \int\mathbf E\cdot\delta\mathbf D\,d\tau$. This step is general: no linearity used.
4. For a linear medium $\mathbf D = \varepsilon\mathbf E$, so $\mathbf E\cdot\delta\mathbf D = \varepsilon\mathbf E\cdot\delta\mathbf E = \delta\left(\tfrac12\varepsilon E^2\right) = \delta\left(\tfrac12\mathbf D\cdot\mathbf E\right)$. Integrate from zero field to the final state: (3.23). ∎

■ **What just happened.** The $\tfrac12$ comes from integrating a linear response ($\int_0^E\varepsilon E'\,dE' = \tfrac12\varepsilon E^2$), not from the double-counting argument of EM-2 (3.14); for a nonlinear material you must keep $\int\mathbf E\cdot\delta\mathbf D$. For a capacitor, (3.23) reduces to $\tfrac12QV = \tfrac12CV^2$ (Example 4.2 and project Part 3 check this).

**What (3.23) includes.** EM-2 (3.16), $W_0 = \frac{\varepsilon_0}2\int E^2d\tau$, is the work to assemble *all* the charge, free and bound, treating the bound charges as if they were free and could be placed by hand. (3.23) is the work *you* do on the free charge while the material polarizes by itself. The difference,
$$W - \frac{\varepsilon_0}{2}\int E^2\,d\tau = \frac12\int\mathbf P\cdot\mathbf E\,d\tau, \tag{3.24}$$
is the energy stored in the stretched molecular springs: for $\mathbf p = \alpha\mathbf E$, each molecule stores $\tfrac12\alpha E^2 = \tfrac12\mathbf p\cdot\mathbf E$, and $N$ per volume gives $\tfrac12\mathbf P\cdot\mathbf E$ (Exercise B6 computes the split for a half-filled capacitor: $15\%$ of the total for $\varepsilon_r = 4$). Both $W_0$ and $W$ are "the energy"; they answer different questions. $W$ is the one that gives forces on free charges and dielectrics by differentiation, which is what you want.

### 3.13 Forces on dielectrics: the slab and the paradox

A dielectric slab is pulled *into* a capacitor. The force comes from the fringing field at the plate edges acting on the slab's bound charges, which is hard to compute directly. The energy method sidesteps it.

**Setup.** Plates of length $\ell$ (along $x$) and width $w$, gap $d$, slab of thickness $d$ and $\varepsilon_r$ inserted a distance $x$ from one end. The filled part and the empty part are capacitors in parallel:
$$C(x) = \frac{\varepsilon_0w}{d}\left[\varepsilon_rx + (\ell - x)\right] = \frac{\varepsilon_0w}{d}\left(\ell + \chi_ex\right),\qquad\frac{dC}{dx} = \frac{\varepsilon_0\chi_ew}{d}.$$

▶ **Try it first.** Compute $F = -dW/dx$ holding $Q$ fixed. Then hold $V$ fixed instead and see what goes wrong.

1. **Fixed $Q$** (isolated capacitor). $W = Q^2/2C$, and the only energy in the problem is the field's, so
   $$F = -\frac{dW}{dx}\bigg|_Q = \frac{Q^2}{2C^2}\frac{dC}{dx} = \frac12V^2\frac{dC}{dx} = \frac{\varepsilon_0\chi_ewV^2}{2d}, \tag{3.25}$$
   positive: toward larger $x$, into the capacitor, where $V = Q/C$ is the instantaneous voltage.
2. **Fixed $V$** (battery attached). Now $W = \tfrac12CV^2$ and $dW/dx = \tfrac12V^2\,dC/dx > 0$: the stored energy *increases* as the slab goes in. Naively $F = -dW/dx$ would push the slab *out*. The paradox.
3. **Resolution.** With $V$ fixed, the battery is part of the system. As $C$ increases by $dC$, charge $dQ = V\,dC$ flows onto the plates, and the battery does work $V\,dQ = V^2\,dC$. Energy balance for a displacement $dx$: battery's work $=$ increase of field energy $+$ work done by the force on the slab, $V^2dC = \tfrac12V^2dC + F\,dx$. So
   $$F = \frac12V^2\frac{dC}{dx},$$
   the same as (3.25). The battery supplies twice the energy increase; half goes into the field, half into pulling the slab. ∎

■ **What just happened.** The force is a property of the instantaneous state (positions, charges, fields), so it cannot depend on what you *later* hold fixed. $F = -\partial W/\partial x$ is correct only when $W$ includes every energy reservoir that exchanges energy during the displacement; at fixed $Q$ there is only the field, at fixed $V$ there is also the battery. (The same trap, with the sign reversed, appears for magnetic forces at fixed current in EM-6, and in thermodynamics as the difference between $\partial U/\partial V$ at fixed $S$ and $\partial F/\partial V$ at fixed $T$, SM-2.) Project Part 6 solves the finite-plate problem with fringing fields on a grid: $C(x)$ is linear with slope $3.000\,\varepsilon_0w/d$ for $\chi_e = 3$, and the forces computed at fixed $Q$ and at fixed $V$ agree.

### 3.14 From $\alpha$ to $\chi_e$: the Clausius–Mossotti relation

If $N$ molecules per unit volume each have polarizability $\alpha$, is $\chi_e = N\alpha/\varepsilon_0$? Only for a dilute gas. The field that polarizes a molecule is not the macroscopic $\mathbf E$ (an average that includes the molecule's own neighbours *and itself*) but the **local field** $\mathbf E_{\text{loc}}$: the field at the molecule's site due to everything except the molecule.

▶ **Try it first.** Carve out a sphere around the molecule, large enough to contain many neighbours. Use Theorem 3.7's split: the macroscopic field is (field at the centre from everything outside the sphere) $+$ (average field from the dipoles inside, $-\mathbf P/3\varepsilon_0$). What does the molecule at the centre actually feel from the dipoles inside the sphere?

1. By §3.7, $\mathbf E = \mathbf E_{\text{outside}} + \left(-\mathbf P/3\varepsilon_0\right)$, where $\mathbf E_{\text{outside}}$ is the field at the centre from all matter outside the sphere. So $\mathbf E_{\text{outside}} = \mathbf E + \mathbf P/3\varepsilon_0$.
2. The molecule feels $\mathbf E_{\text{outside}}$ plus the field at its site from the *other* dipoles inside the sphere, $\mathbf E_{\text{near}}$. For a cubic lattice, or a random isotropic arrangement, $\mathbf E_{\text{near}} = 0$ by symmetry: the dipole field (3.4) summed over a cubic shell of neighbours cancels. (Project Part 2's lattice: the field at an interior site is $\le0.006\,P/\varepsilon_0$, against the $-0.333\,P/\varepsilon_0$ of the macroscopic field. The two differ by exactly the $+P/3\varepsilon_0$ of this step.)
3. So the **Lorentz local field** is
   $$\mathbf E_{\text{loc}} = \mathbf E + \frac{\mathbf P}{3\varepsilon_0}. \tag{3.26}$$
4. Each molecule acquires $\mathbf p = \alpha\mathbf E_{\text{loc}}$, so $\mathbf P = N\alpha\left(\mathbf E + \mathbf P/3\varepsilon_0\right)$. Solve for $\mathbf P$: $\mathbf P\left(1 - \dfrac{N\alpha}{3\varepsilon_0}\right) = N\alpha\mathbf E$, so with $x\equiv N\alpha/3\varepsilon_0$,
   $$\chi_e = \frac{N\alpha/\varepsilon_0}{1 - N\alpha/3\varepsilon_0} = \frac{3x}{1 - x},\qquad\text{equivalently}\qquad\frac{\varepsilon_r - 1}{\varepsilon_r + 2} = \frac{N\alpha}{3\varepsilon_0}. \tag{3.27}$$
   (From $\varepsilon_r - 1 = 3x/(1-x)$: $(\varepsilon_r - 1)(1 - x) = 3x$, so $\varepsilon_r - 1 = x(\varepsilon_r + 2)$.) ∎

■ **What just happened.** The local field exceeds the macroscopic field because the molecule does not feel its own depolarizing contribution. The result is that $\chi_e$ grows *faster* than linearly with density, and (3.27) predicts a divergence at $N\alpha = 3\varepsilon_0$: the polarization catastrophe, in which the dipoles sustain each other with no applied field. Real ferroelectrics (§3.15) are roughly this. For gases $x\sim10^{-4}$ and the correction is negligible (Example 4.5); for liquid argon it is a $20\%$ effect and (3.27) works; for water it fails for the static $\varepsilon_r$ (Exercise C1), because a polar liquid's orientational response is not a fixed $\alpha$, the neighbours are hydrogen-bonded rather than isotropic, and $\mathbf E_{\text{near}}\ne0$. Debye's extension adds $p^2/3k_BT$ to $\alpha$; Onsager's (1936) fixes the local field. With $\alpha$ the optical (electronic) polarizability and $\varepsilon_r = n^2$, (3.27) is the **Lorentz–Lorenz** relation for the refractive index.

### 3.15 Three previews in one subsection

**Why $\varepsilon_r$ depends on frequency.** Make the spring of §3.1 explicit: an electron of mass $m$, charge $-e$, bound by a restoring force $-m\omega_0^2x$ and damped by $-m\gamma\dot x$, is driven by $E_0e^{-i\omega t}$. The steady state of $m\ddot x = -m\omega_0^2x - m\gamma\dot x - eE$ is $x = -\dfrac{e/m}{\omega_0^2 - \omega^2 - i\gamma\omega}E$, so $p = -ex$ gives the **Lorentz oscillator** polarizability
$$\alpha(\omega) = \frac{e^2/m}{\omega_0^2 - \omega^2 - i\gamma\omega}, \tag{3.28}$$
complex and frequency-dependent. At $\omega\ll\omega_0$ it is the static $e^2/m\omega_0^2$; near $\omega_0$ it is large and imaginary (absorption); far above it is negative and small (free-electron-like). Then $\varepsilon_r(\omega) = 1 + N\alpha(\omega)/\varepsilon_0$ (dilute) or via (3.27), and the refractive index $n = \sqrt{\varepsilon_r}$ varies with colour: dispersion. EM-7 derives everything that follows from (3.28), including why glass is transparent and why water's $\varepsilon_r$ drops from $80$ to $n^2 = 1.78$ at optical frequencies (the molecules cannot rotate at $10^{15}$ Hz; only the electrons respond, with $\alpha/4\pi\varepsilon_0 = 1.45$ Å³, and (3.27) then gives $1.76$).

**Ferroelectrics.** In BaTiO$_3$ below $120$ °C the Ti$^{4+}$ ion sits off-centre in its oxygen cage, every unit cell carries a permanent dipole, and the crystal has a spontaneous $\mathbf P$ with no applied field: the polarization catastrophe realized. $\mathbf P$ versus $\mathbf E$ traces a hysteresis loop (remanent polarization, coercive field), the crystal breaks into domains of different $\mathbf P$ direction, and $\varepsilon_r$ is huge, nonlinear, and temperature-dependent. Everything here has a magnetic twin in EM-5's ferromagnetism, with the same loop.

**The engineering angle.** (i) Every dielectric breaks down above a **dielectric strength** $E_{\max}$: about $3$ MV/m for air, tens of MV/m for polymer films, over $100$ MV/m for mica. The maximum energy density a capacitor can store is $\tfrac12\varepsilon E_{\max}^2$, which is why high-energy-density capacitors use thin polymer films (high $E_{\max}$) and ceramics (high $\varepsilon_r$) and why air gaps spark at $3$ kV/mm. (ii) In power electronics, film capacitors ($\varepsilon_r\approx2$–$3$, linear, low loss) do the filtering and snubbing, while BaTiO$_3$-based ceramic capacitors pack $\mu$F into millimetres at the price of capacitance that changes with voltage and temperature: (3.18) with a field-dependent $\varepsilon_r$. (iii) **Piezoelectricity**: in crystals without a centre of symmetry (quartz, PZT), stress produces polarization, $P_i = d_{ijk}\sigma_{jk}$, and conversely a field produces strain, $\epsilon_{jk} = d_{ijk}E_i$. The piezoelectric tensor $d_{ijk}$ is a third-rank tensor (LA-12) coupling the stress tensor of CM-7's continuum mechanics to $\mathbf P$; for quartz $d_{11} = 2.3$ pC/N, for PZT several hundred. A stressed plate acquires $\sigma_b = \mathbf P\cdot\hat{\mathbf n}$ on its faces and an open-circuit voltage (Exercise C2); this is the pressure sensor, the ultrasound transducer, and the quartz clock.

## 4. Worked examples

### Example 4.1: The uniformly polarized sphere, with numbers (Part 2)

$R = 1$ cm, $P = 1\ \mu$C/m² along $\hat{\mathbf z}$. By §3.6: $\rho_b = 0$, $\sigma_b = P\cos\theta$, with maximum $10^{-6}$ C/m² at the poles and zero at the equator. Total bound charge: $\int P\cos\theta\,R^2\sin\theta\,d\theta\,d\phi = 0$ ✓ (3.10).

*Inside:* $E_{\text{in}} = -P/3\varepsilon_0 = -10^{-6}/(3\times8.854\times10^{-12}) = -3.7647\times10^4$ V/m, uniform, along $-\hat{\mathbf z}$.

*Outside:* pure dipole with $p = \tfrac43\pi R^3P = \tfrac43\pi(10^{-6})(10^{-6}) = 4.1888\times10^{-12}$ C m. Just outside the north pole ($r = R$, $\theta = 0$), (3.4) gives $E_r = 2p/(4\pi\varepsilon_0R^3) = 2P/3\varepsilon_0 = +7.5294\times10^4$ V/m. The jump across the pole is $7.5294\times10^4 - (-3.7647\times10^4) = 1.1294\times10^5 = P/\varepsilon_0 = \sigma_b(0)/\varepsilon_0$ ✓ (EM-2 (3.10)).

*Potential at the pole* (continuous): $V = PR/3\varepsilon_0 = 376.47$ V. Project Part 2 sums the Coulomb fields of the $\cos\theta$ surface charge on $4608$ patches: interior field uniform to $2\times10^{-10}$, exterior dipole to $2\times10^{-14}$.

### Example 4.2: A metal sphere in a dielectric shell (Part 5)

A metal sphere of radius $a = 1$ cm carries $Q = 1$ nC and is surrounded by a concentric shell of linear dielectric, $\varepsilon_r = 4$, out to $b = 2$ cm. Find $\mathbf D$, $\mathbf E$, $\mathbf P$, all bound charges, the potential of the sphere, and the energy.

1. *$\mathbf D$ by (3.15).* Spherical symmetry puts $\mathbf D$ radial; a Gaussian sphere of radius $r > a$ encloses free charge $Q$: $\mathbf D = \dfrac{Q}{4\pi r^2}\hat{\mathbf r}$ for all $r > a$ (and $0$ inside the metal). The dielectric does not appear.
2. *$\mathbf E = \mathbf D/\varepsilon$.* In the shell, $\mathbf E = \dfrac{Q}{4\pi\varepsilon_0\varepsilon_rr^2}\hat{\mathbf r}$; outside, $\dfrac{Q}{4\pi\varepsilon_0r^2}\hat{\mathbf r}$. Numbers: $E(a^+) = 8.988\times10^9\times10^{-9}/(4\times10^{-4}) = 22469$ V/m; $E(b^-) = 5617$ V/m; $E(b^+) = 22469$ V/m. $E$ jumps by a factor $\varepsilon_r$ at $b$; $D$ does not (3.16).
3. *$\mathbf P = \mathbf D - \varepsilon_0\mathbf E = \left(1 - \dfrac1{\varepsilon_r}\right)\mathbf D = \dfrac{\chi_e}{\varepsilon_r}\dfrac{Q}{4\pi r^2}\hat{\mathbf r}$* in the shell. $\rho_b = -\nabla\cdot\mathbf P = 0$ (the $1/r^2$ radial field is divergence-free, as (3.19) promised with $\rho_f = 0$ in the shell).
4. *Bound surface charges (3.8).* Inner surface $r = a$, outward normal of the dielectric $-\hat{\mathbf r}$: $\sigma_b(a) = -\dfrac{\chi_e}{\varepsilon_r}\dfrac{Q}{4\pi a^2} = -\dfrac34\times\dfrac{10^{-9}}{4\pi\times10^{-4}} = -5.968\times10^{-7}$ C/m². Outer surface, normal $+\hat{\mathbf r}$: $\sigma_b(b) = +\dfrac{\chi_e}{\varepsilon_r}\dfrac{Q}{4\pi b^2} = +1.492\times10^{-7}$ C/m². Totals: $Q_b(a) = -\dfrac{\chi_e}{\varepsilon_r}Q = -0.75$ nC, $Q_b(b) = +0.75$ nC, sum zero ✓ (3.10).
5. *Check by Gauss for $\mathbf E$ with total charge.* Just outside the metal the enclosed charge is $Q + Q_b(a) = Q/\varepsilon_r = 0.25$ nC, giving $E = kQ/(\varepsilon_ra^2)$ ✓. Outside $b$, $Q + Q_b(a) + Q_b(b) = Q$ ✓.
6. *Potential of the sphere.* $V(a) = \displaystyle\int_a^\infty E\,dr = \frac{Q}{4\pi\varepsilon_0}\left[\frac1{\varepsilon_r}\left(\frac1a - \frac1b\right) + \frac1b\right] = 8.988\left[\frac14(100 - 50) + 50\right] = 561.72$ V, against $898.76$ V with no dielectric. Capacitance $C = Q/V(a) = 1.7802$ pF.
7. *Energy (3.23).* $W = \tfrac12\int\mathbf D\cdot\mathbf E\,d\tau = \tfrac12QV(a) = 2.8086\times10^{-7}$ J. The project does the radial integral $\int\tfrac12\varepsilon E^2\,4\pi r^2dr$ by quadrature and matches to $10^{-6}$.

The pattern "$\mathbf D$ from free charge and symmetry, then $\mathbf E$, then $\mathbf P$, then bound charges, then check they sum to zero" is the module's standard procedure.

### Example 4.3: A capacitor half-filled, two ways (Part 3)

Plates of area $A = 100$ cm², gap $d = 2$ mm, $C_0 = \varepsilon_0A/d = 44.27$ pF. A slab of $\varepsilon_r = 4$ fills half the volume. Two geometries.

**(a) Series: a slab of thickness $d/2$ against one plate.** Free charge $\pm\sigma_f$ on the plates. By (3.16) at each plate, $D = \sigma_f$ in both the slab and the gap, uniform (no free charge inside, planar symmetry). So $E_{\text{slab}} = \sigma_f/\varepsilon_r\varepsilon_0$, $E_{\text{gap}} = \sigma_f/\varepsilon_0$, and
$$V = E_{\text{slab}}\frac d2 + E_{\text{gap}}\frac d2 = \frac{\sigma_fd}{2\varepsilon_0}\left(\frac1{\varepsilon_r} + 1\right),\qquad C = \frac{\sigma_fA}{V} = \frac{2\varepsilon_r}{\varepsilon_r + 1}C_0 = 1.6\,C_0 = 70.83\text{ pF}.$$
(Two capacitors in series, EM-2 §3.10: $1/C = 1/(2\varepsilon_rC_0) + 1/(2C_0)$.) At $V = 100$ V: $\sigma_f = D = CV/A = 7.0834\times10^{-7}$ C/m²; $E_{\text{gap}} = D/\varepsilon_0 = 8.0\times10^4$ V/m; $E_{\text{slab}} = 2.0\times10^4$ V/m; check $2\times10^4\times10^{-3} + 8\times10^4\times10^{-3} = 100$ V ✓. Bound charge on the slab's faces: $P = D(1 - 1/\varepsilon_r) = 5.3125\times10^{-7}$ C/m², $\sigma_b = -P$ on the face touching the positive plate and $+P$ on the face toward the gap. The net charge the gap field sees at the slab's far face is $\sigma_b = +P$, and at the positive plate $\sigma_f - P = \sigma_f/\varepsilon_r$, which is why $E_{\text{slab}}$ is a quarter of $E_{\text{gap}}$.

**(b) Parallel: a slab of full thickness $d$ covering half the area.** Both halves see the same $V$, so $E = V/d$ everywhere (tangential continuity (3.17) at the slab's vertical edge is automatic), and $D$ is $\varepsilon_rε_0V/d$ in the slab and $\varepsilon_0V/d$ in the gap. The free charge is $Q = \varepsilon_r\varepsilon_0\frac Vd\frac A2 + \varepsilon_0\frac Vd\frac A2$, so
$$C = \frac{\varepsilon_r + 1}{2}C_0 = 2.5\,C_0 = 110.68\text{ pF}.$$
Half the plate carries four times the charge density of the other half.

The project solves both on a $40\times40$ grid with harmonic-mean face permittivities: $C/C_0 = 1.6000$ and $2.5000$, $D_y$ constant across the series interface to machine precision, $E_y$ jumping by exactly $4$, and $V$ exactly linear in (b). `slab_capacitor_DE.png` is the $D$-versus-$E$ plot through the interface.

### Example 4.4: A dielectric sphere in a uniform field, $\varepsilon_r = 4$ (Part 4)

$E_0 = 10^5$ V/m, $R = 2$ cm, $\varepsilon_r = 4$. By (3.22): $E_{\text{in}} = 3E_0/(\varepsilon_r + 2) = \tfrac36\times10^5 = 5.000\times10^4$ V/m. $P = \varepsilon_0\chi_eE_{\text{in}} = 3\varepsilon_0\times5\times10^4 = 1.3281\times10^{-6}$ C/m² (equivalently $3\varepsilon_0\tfrac{3}{6}E_0 = 1.5\varepsilon_0E_0$). $p = \tfrac43\pi R^3P = 4.4506\times10^{-11}$ C m (equivalently $4\pi\varepsilon_0R^3\tfrac{\varepsilon_r-1}{\varepsilon_r+2}E_0 = 2\pi\varepsilon_0R^3E_0$).

*Field just outside the pole:* $E_0 + 2p/(4\pi\varepsilon_0R^3) = E_0 + 2\times\tfrac12E_0 = 2E_0 = 2.000\times10^5$ V/m. *Jump at the pole:* $2\times10^5 - 5\times10^4 = 1.5\times10^5 = P/\varepsilon_0 = \sigma_b(0)/\varepsilon_0$ ✓. *$D^\perp$ continuity (3.16):* inside, $\varepsilon_r\varepsilon_0E_{\text{in}} = 4\varepsilon_0\times5\times10^4 = 2\varepsilon_0\times10^5$; outside, $\varepsilon_0\times2\times10^5$ ✓. *Self-consistency check:* $E_0 - P/3\varepsilon_0 = 10^5 - 5\times10^4 = 5\times10^4 = E_{\text{in}}$ ✓.

The project confirms (3.22) symbolically by Route 2 and solves the 2D analogue (a cylinder) on a $241^2$ grid: $E_{\text{in}}/E_0 = 0.407$ against $2/(\varepsilon_r + 1) = 0.400$, with the $1.7\%$ error from staircasing the circle and shrinking as the grid is refined; for $\varepsilon_r = 80$ the cylinder is nearly a conductor, $E_{\text{in}}/E_0 = 0.025$.

### Example 4.5: Clausius–Mossotti for nitrogen at STP (Part 7; a TSM example)

Estimate $\chi_e$ of N$_2$ at $0$ °C and $1$ atm from its molecular polarizability and compare with the tabulated value.

1. *Number density* from the ideal gas law (SM-1): $N = p/k_BT = 101325/(1.3807\times10^{-23}\times273.15) = 2.6868\times10^{25}$ m$^{-3}$.
2. *Polarizability.* Tabulated $\alpha/4\pi\varepsilon_0 = 1.74$ Å³ for N$_2$. The crude model (3.2) with $a = 1.2$ Å (about the N–N bond length plus the electron cloud) gives $a^3 = 1.73$ Å³: the right size.
3. *(3.27):* $x = N\alpha/3\varepsilon_0 = \tfrac{4\pi}{3}N(\alpha/4\pi\varepsilon_0) = \tfrac{4\pi}3\times2.6868\times10^{25}\times1.74\times10^{-30} = 1.958\times10^{-4}$. Then $\chi_e = 3x/(1 - x) = 5.876\times10^{-4}$; the dilute formula $3x$ gives $5.875\times10^{-4}$, so the local-field correction is one part in $5000$.
4. *Compare.* Griffiths' table gives $\varepsilon_r = 1.00055$ at $20$ °C, i.e. $\chi_e = 5.5\times10^{-4}$; since $\chi_e\propto N\propto1/T$ at fixed pressure, at $0$ °C that is $5.5\times10^{-4}\times293.15/273.15 = 5.90\times10^{-4}$. Our estimate is within $0.5\%$ of it, and within $7\%$ of the raw $20$ °C figure. ∎

The lesson: for a gas, $\chi_e = N\alpha/\varepsilon_0$ to excellent accuracy, and a one-line atom model gets $\alpha$ to $1\%$ once $a$ is sensible. The project also runs (3.27) for liquid argon ($\varepsilon_r = 1.509$ against measured $\approx1.5$: works) and water (Exercise C1: fails), and plots $\chi_e$ against $N\alpha/\varepsilon_0$ to show the superlinear rise.

## 5. Why it matters and where it connects

- **EM-5.** Replace $\mathbf P\to\mathbf M$, $\sigma_b = \mathbf P\cdot\hat{\mathbf n}\to\mathbf K_b = \mathbf M\times\hat{\mathbf n}$, $\rho_b = -\nabla\cdot\mathbf P\to\mathbf J_b = \nabla\times\mathbf M$, $\mathbf D\to\mathbf H$, and Theorem 3.3's derivation (one integration by parts of the dipole potential) is repeated for the vector potential. The uniformly magnetized sphere has $\mathbf B_{\text{in}} = \tfrac23\mu_0\mathbf M$, the twin of (3.12). The "no Coulomb's law for $\mathbf D$" warning becomes "no Biot–Savart law for $\mathbf H$."
- **EM-6.** Maxwell's equations in matter are (3.14), $\nabla\times\mathbf H = \mathbf J_f + \partial\mathbf D/\partial t$, and the vacuum pair. The $\partial\mathbf D/\partial t$ is Maxwell's displacement current, and the energy density $\tfrac12(\mathbf D\cdot\mathbf E + \mathbf B\cdot\mathbf H)$ generalizes (3.23).
- **EM-7.** Plane waves in a linear medium travel at $c/\sqrt{\varepsilon_r\mu_r}$; the refractive index is $n = \sqrt{\varepsilon_r}$; reflection at an interface is the boundary conditions (3.16)–(3.17) applied to oscillating fields (Fresnel); dispersion is (3.28) put into (3.27).
- **EM-3.** The $\cos\theta$ shell (Example 4.5 there) and the uniform-field sphere (Example 4.4 there) are the two EM-3 problems this module reuses; §3.11 is the dielectric version of the conducting sphere, and uniqueness is what makes Route 1 legitimate.
- **QM-8.** The static polarizability of hydrogen, $\alpha/4\pi\varepsilon_0 = \tfrac92a_0^3$, is a second-order perturbation-theory result; (3.1) with $U = -\tfrac12\alpha E^2$ is the quadratic Stark effect.
- **SM-3, SM-6.** The Langevin average $\langle p_z\rangle = p^2E/3k_BT$ of §3.2 is the first Boltzmann-weighted orientation average you meet; Debye's $\alpha + p^2/3k_BT$ in (3.27) is why $\varepsilon_r$ of polar liquids falls with $T$. The induced-dipole interaction $\mathbf p = \alpha\mathbf E$ with a neighbour's dipole field gives the Debye and London $1/r^6$ van der Waals terms used in SM-6.
- **LA-12, CM-7.** The polarizability and susceptibility of a crystal are symmetric second-rank tensors; the piezoelectric coefficient is a third-rank tensor coupling the stress tensor to $\mathbf P$. The stress tensor itself is CM-7's continuum object.
- **LA-13.** Project Part 3's solver is the five-point finite-difference Laplacian with variable coefficients, assembled sparse and solved directly; the harmonic-mean face permittivity is the standard treatment of discontinuous coefficients in finite-volume methods.
- **Engineering.** Capacitor design is (3.21) plus the breakdown field; the force on a dielectric (3.25) is the principle of electrostatic actuators and of the dielectrophoretic manipulation of cells; piezoelectric transducers are §3.15(iii).

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 (warm-up) drill the polarizability model, bound charges of a simple $\mathbf P$, torque against $k_BT$, a filled capacitor's charges, and refraction of field lines. B1–B8 are the module: B1 and B3 (radial $\mathbf P$ in a sphere, fields two ways), B2 (total bound charge and total dipole moment, proofs), B4 (point charge inside a dielectric sphere), B5 (two-layer capacitor, interface bound charge), B6 (the energy split (3.24)), B7 (cylinder in a uniform field by separation of variables, then `project.py`), B8 (force on a slab, both ways, with numbers). B1, B4, B5, B7, B8 are essential. C1 (Clausius–Mossotti where it works and where it fails), C2 (a piezoelectric plate), C3 (prove the averaging lemma).

## 7. Coding project

`project.py` (about 9 s) checks every number above and builds the module's pictures:

1. **Bound charges, symbolically and on a grid.** SymPy verifies (3.2), (3.9), and $\rho_b$, $\sigma_b$, and (3.10) for $\mathbf P = kr\hat{\mathbf r}$ and $kr^2\hat{\mathbf r}$; `bound_charge_on_grid` computes $-\nabla\cdot\mathbf P$ by central differences for an irregular bounded $\mathbf P$ and finds $Q_b = 10^{-17}$ against $\int|\rho_b| = 1.8$.
2. **The polarized sphere three ways** (`coulomb_field`, `dipole_field`, `ball_field`): surface patches of $\sigma_b = P\cos\theta$ (Gauss–Legendre in $\cos\theta$), two offset balls, and a cubic lattice of $2109$ dipoles. The first two give (3.12) and the exterior dipole; the lattice gives $\approx0$ at the sites and the dipole field outside. Then Lemma 3.6 by quadrature. `polarized_sphere.png`.
3. **A finite-difference solver** for $\nabla\cdot(\varepsilon_r\nabla V) = -\rho_f/\varepsilon_0$: `solve_poisson2d` assembles the flux-balance equations $\sum_f\varepsilon_f(V_{\text{nb}} - V_i) = -\rho_ih^2/\varepsilon_0$ on each node's dual cell as a sparse matrix and solves directly; `node_charges` returns $\oint\mathbf D\cdot d\mathbf a$ around each node (so plate charge, hence $C$); `field_energy` returns $\tfrac12\sum_f\varepsilon_f(\Delta V)^2$, the discrete (3.23); `harmonic_faces` gives the face permittivity $2\varepsilon_1\varepsilon_2/(\varepsilon_1 + \varepsilon_2)$, exact for an interface midway between nodes because the two half-cells are in series for the normal flux. Example 4.3 both ways; $D_y$ continuous, $E_y$ jumping; $W = \tfrac12CV^2$ to $10^{-10}$. `slab_capacitor_DE.png`.
4. **Cylinder in a uniform field** (`cylinder_in_field`) for five $(\varepsilon_{\text{in}}, \varepsilon_{\text{out}})$ pairs against $2\varepsilon_{\text{out}}/(\varepsilon_{\text{in}} + \varepsilon_{\text{out}})$, a convergence check, and the sphere (3.22) by SymPy. `cylinder_uniform_field.png`.
5. **Energy.** Example 4.2 by quadrature of $\tfrac12\varepsilon E^2$; Examples 4.1, 4.3, 4.4 numbers.
6. **Force on a slab.** `slab_capacitor` solves finite plates of length $8d$ inside a grounded box with the slab inserted to $x$; $C(x)$ from the plate charge is linear with slope $\chi_e\varepsilon_0w/d$ (fit $3.000$), and $\tfrac12V^2\,dC/dx$ at fixed $V$ agrees with $-d(Q^2/2C)/dx$ at fixed $Q$. `slab_force.png`.
7. **Clausius–Mossotti** (`clausius_mossotti`): Example 4.5, argon, water, and every exercise number. `clausius_mossotti.png`.

**Build it yourself** before reading Part 3: write `solve_poisson2d(eps_xf, eps_yf, rho, h, fixed, fixed_val)` and `node_charges`. The one idea is that each face between two nodes carries one permittivity and one flux $\varepsilon_f(V_{\text{nb}} - V_i)$; summing over a node's four faces gives $\oint\mathbf D\cdot d\mathbf a$ over its dual cell, which must equal the free charge there. Dirichlet nodes get the trivial row $V_i = V_{\text{fixed}}$. Then compare to the reference, which is about $30$ lines. **Your turn** (end of the file): a coaxial cable with two concentric dielectrics ($\mathbf D$ from (3.15), then $C'$, then $\sigma_b$ at the interface and the zero total), Clausius–Mossotti for water (where it fails and why), and a piezoelectric quartz plate (induced $\sigma_b$, charge, open-circuit voltage). Reference solutions follow the prompts; write yours first.

## 8. Where you'll struggle

- **You will want $\mathbf D$ to be "the field in the material." It is not.** $\mathbf E$ is the field: it is what a test charge feels, what appears in the energy density (3.23) together with $\mathbf D$, what has a potential. $\mathbf D$ is $\varepsilon_0\mathbf E$ plus the material's dipole density, defined because its divergence is convenient. Symptoms of the confusion: writing "$\mathbf D = \varepsilon_0\mathbf E$ outside the dielectric, so $\mathbf D$ is the vacuum field" (true, and irrelevant inside); trying to compute $\mathbf D$ from free charge in a geometry without symmetry (§3.8: no Coulomb's law for $\mathbf D$); forgetting that $\nabla\times\mathbf D\ne0$ at the edge of a slab. The fix: in every problem, say which field each equation is about, and get $\mathbf E$ before you say anything physical.
- **The averaging argument of §3.7 is genuinely subtle and most courses hand-wave it.** Griffiths puts it in §4.2.3 and most readers skip it. The point is not that the macroscopic field is "approximately" the dipole-sum field; it is that the *average* over a sphere of the true microscopic field of the nearby molecules equals $-\mathbf P/3\varepsilon_0$ exactly, by Lemma 3.6, which is the same number the dipole approximation gives. If you cannot reconstruct Lemma 3.6(b) (field averaged over a ball equals minus the field of a uniform ball at the charge's position), reread §3.7; the project's lattice, where the field at a site is $0$ while the macroscopic field is $-P/3\varepsilon_0$, is the picture to hold.
- **Which energy, and why the $\tfrac12$.** $\tfrac12\int\mathbf D\cdot\mathbf E$ is the work done on free charges with the material responding (it includes the spring energy of the dipoles); $\tfrac{\varepsilon_0}2\int E^2$ is the work to assemble all charges as if free (it does not). They differ by $\tfrac12\int\mathbf P\cdot\mathbf E$ (3.24). The $\tfrac12$ in (3.23) is from linear response, $\int\varepsilon E\,dE$, not double counting, and for a nonlinear dielectric you must integrate $\mathbf E\cdot\delta\mathbf D$ along the actual path. Use (3.23) for forces on free charges and dielectrics; it is the one whose derivative is the force.
- **The force-on-a-slab paradox catches everyone.** At fixed $V$ the field energy *increases* as the slab goes in, and $-dW/dx$ gives the wrong sign. The force is a property of the state and cannot depend on what you hold fixed; the battery's work $V^2dC$ is twice the field-energy increase. Always ask "what reservoirs exchange energy during this virtual displacement?" before writing $F = -dW/dx$. (Expect the same trap, with the opposite sign of the battery's contribution, for forces at fixed current in EM-6.)
- **Signs of bound charges.** $\sigma_b = \mathbf P\cdot\hat{\mathbf n}$ with $\hat{\mathbf n}$ pointing *out of the dielectric*. At the inner surface of a shell that is $-\hat{\mathbf r}$ (Example 4.2 step 4). At an interface between two dielectrics use both outward normals: $\sigma_b = (\mathbf P_1 - \mathbf P_2)\cdot\hat{\mathbf n}_{12}$ (Exercise B5). The check that catches every sign error: the total bound charge of an isolated object is zero (3.10), and the bound charge adjacent to a free charge has the opposite sign.
- **The dielectric-sphere argument looks circular.** "Assume uniform $\mathbf P$, find a uniform field, conclude uniform $\mathbf P$" proves only that uniform $\mathbf P$ is *consistent*. Uniqueness (EM-3 §3.2) is what upgrades consistent to correct, and Route 2 of §3.11 is the independent check. Do both once; afterwards use Route 1, which is three lines.
- **$\mathbf E$ in (3.18) is the total field.** Students write $\mathbf P = \varepsilon_0\chi_e\mathbf E_0$ with the *applied* field and get the sphere wrong by the factor $3/(\varepsilon_r + 2)$. The material responds to the field that is actually there, including its own.
- **Bound volume charge in a linear dielectric.** By (3.19) there is none unless free charge is embedded, *and the medium is homogeneous*. A dielectric whose $\varepsilon_r$ varies in space does have $\rho_b = -\nabla\cdot\mathbf P\ne0$ even with no free charge (Exercise B3 asks why the radial-$\mathbf P$ sphere, which is not a linear dielectric, has $\rho_b\ne0$).
- **Exercise B7 is harder than it looks**: the cylindrical separation of variables is EM-3 (3.24) with only the $n = 1$ term, but the matching condition is $\varepsilon_r\partial V_{\text{in}}/\partial s = \partial V_{\text{out}}/\partial s$, not continuity of $\partial V/\partial s$. Write the two matching conditions before touching the algebra.

## 9. Check yourself

1. Define $\mathbf P$, $\sigma_b$, $\rho_b$, and state the sign convention for $\hat{\mathbf n}$.
2. Derive $\sigma_b = \mathbf P\cdot\hat{\mathbf n}$ and $\rho_b = -\nabla\cdot\mathbf P$ from the dipole potential in five lines. Which identity and which theorem are used?
3. Why is the total bound charge of an isolated polarized object zero?
4. A uniformly polarized sphere: what are $\sigma_b$, the interior field, and the exterior field? Why is the interior field a third of the slab's $P/\varepsilon_0$?
5. State the averaging lemma for charges inside a sphere, and explain in two sentences why it makes the macroscopic field computable from $\sigma_b$ and $\rho_b$.
6. Define $\mathbf D$ and derive $\nabla\cdot\mathbf D = \rho_f$. Is $\nabla\times\mathbf D = 0$? Can $\mathbf D$ be computed from $\rho_f$ by a Coulomb-type integral?
7. Write the two boundary conditions at a dielectric interface with no free charge. Which component of which field jumps, and by how much?
8. For a homogeneous linear dielectric, what is $\rho_b$ in terms of $\rho_f$? What follows for a dielectric with no free charge inside?
9. A dielectric sphere in a uniform field: write $\mathbf E_{\text{in}}$ and give the self-consistency argument in three lines. What theorem makes it legitimate?
10. Write the energy of a linear dielectric system. What does it include that $\tfrac{\varepsilon_0}2\int E^2$ does not?
11. A slab is pulled into a capacitor at fixed $V$. The field energy increases. Where does the energy come from, and what is the force?
12. State the Clausius–Mossotti relation and the local-field assumption behind it. For which substances does it fail, and why?
13. Why does $\varepsilon_r$ of water drop from $80$ to about $1.8$ at optical frequencies?

## 10. Pair with

- **Griffiths, *Introduction to Electrodynamics*, ch. 4.** Same notation, same order, same conventions as these notes. §4.2.3 is the averaging argument of §3.7 in his words; §4.4.1 has two more Gauss-for-$\mathbf D$ examples; §4.4.3–4.4.4 has the energy and force discussions with the paradox stated plainly. Problems 4.10–4.20 (bound charges and $\mathbf D$), 4.23–4.26 (linear dielectrics, the sphere), 4.28 (force), 4.37–4.40 (Clausius–Mossotti, Debye) are the set to do.
- **Purcell & Morin, *Electricity and Magnetism*, ch. 10.** Starts from the atoms and builds up, with more physics about *why* molecules polarize and excellent numerical estimates (atomic polarizabilities, the water molecule, the Lorentz field). Its §10.11 on the field inside a dielectric is a good second reading of §3.7; its treatment of the local field is more careful than Griffiths'. Gaussian units: multiply by $4\pi\varepsilon_0$ as needed.
- **Feynman, *Lectures on Physics* II, ch. 10–11.** Ch. 10 is the clearest short account of bound charge and of why $\mathbf D$ is bookkeeping; ch. 11 ("Inside dielectrics") does the Clausius–Mossotti and Debye arguments, the polarization catastrophe, ferroelectricity, and piezoelectricity in Feynman's style, which is the best place to read §3.14–3.15 slowly.
