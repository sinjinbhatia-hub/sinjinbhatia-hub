# SM-2: Thermodynamics II: Entropy, the Second Law, Engines, Thermodynamic Potentials, Maxwell Relations, Legendre Transforms

**Time:** ~10 hours · **Prerequisites:** SM-1 · **Unlocks:** SM-3, CM-8 (the Legendre transform parallel)

> **Conventions in this module.** As in SM-1: $N$ particles or $n$ moles, $k_B$ Boltzmann's constant, $R = N_Ak_B = 8.314\ \mathrm{J\,mol^{-1}K^{-1}}$, ideal gas $PV = Nk_BT = nRT$, $U = \tfrac f2 Nk_BT$ for $f$ quadratic degrees of freedom, $\gamma = C_P/C_V$, adiabats $TV^{\gamma-1} = $ const. **Sign convention for the first law:** $dU = đQ - đW$ with $đQ$ the heat *into* the system and $đW$ the work done *by* it; $đ$ marks an inexact differential (SM-1). In cycle bookkeeping $Q_h > 0$ is the heat taken from the hot reservoir and $Q_c > 0$ the heat dumped into the cold one, except in (3.4) and §3.4, where every heat is signed as heat *into* the system. Entropy $S$ is defined by $dS = đQ_{\text{rev}}/T$. Potentials and their natural variables: $U(S,V,N)$, $H(S,P,N) = U + PV$, $F(T,V,N) = U - TS$, $G(T,P,N) = U - TS + PV$. Chemical potential $\mu = (\partial U/\partial N)_{S,V}$; thermodynamic identity $dU = T\,dS - P\,dV + \mu\,dN$. $(\partial x/\partial y)_z$ means the derivative with $z$ held fixed, and the subscript is never optional. Thermal expansion $\alpha = \frac1V(\partial V/\partial T)_P$, isothermal compressibility $\kappa_T = -\frac1V(\partial V/\partial P)_T$. Temperatures are in kelvin unless marked °C.

## 0. Where this sits

SM-1 gave you the first law, $U$ as a state function, and the ideal-gas processes (isothermal, adiabatic, isochoric, isobaric). It could tell you how much heat a process takes but not whether the process can happen at all. This module adds the second law and the state function it implies, entropy, then builds the machinery that makes thermodynamics usable: the potentials $H, F, G$ as Legendre transforms of $U$, the Maxwell relations, the chemical potential, and phase coexistence.

SM-3 explains *why* entropy exists ($S = k_B\ln\Omega$) and fixes the constant this module leaves undetermined. CM-8's Legendre transform (Definition 3.1 there) is the same mathematics as §3.8 here; whichever you do second, you should find it familiar. SM-6 returns to the van der Waals loop of §3.11 as a mean-field theory.

## 1. The problem that forced this idea

Steam engines came first and theory second. In 1824 Sadi Carnot, a French engineer, published *Réflexions sur la puissance motrice du feu*. He asked a practical question: is there a limit to how much work you can get out of a given amount of fuel, and does it depend on the working substance (steam, air, alcohol)? His answer was that the limit depends only on the temperatures between which the engine runs, and that a reversible engine attains it. He argued this from the impossibility of perpetual motion. He did it while believing heat was a conserved fluid, "caloric", that falls from hot to cold the way water falls through a mill wheel. The argument was right; the picture underneath it was wrong. Émile Clapeyron put Carnot's argument into equations and drew it on a $PV$ diagram in 1834, which is how it reached later physicists.

By the late 1840s Joule's experiments had shown heat is not conserved: work turns into heat at a fixed exchange rate, which is SM-1's first law. That seemed to destroy Carnot's argument, which needed something to be conserved through the engine. William Thomson (later Lord Kelvin) saw in 1848 that Carnot's result gave a temperature scale independent of any substance, and in 1851 stated the second law in the form now named after him. Rudolf Clausius (1850) reconciled Carnot with Joule: energy is conserved, *and* there is a second, independent principle, that heat does not flow by itself from cold to hot. In 1865 he identified the conserved-in-reversible-processes quantity hidden in Carnot's argument, $\int đQ/T$, showed it defines a state function, and named it *entropy* from the Greek for "transformation". He summarized: the energy of the world is constant; the entropy of the world tends to a maximum.

Josiah Willard Gibbs, at Yale, then turned the subject into geometry. In papers from 1873 to 1878 (the long one is *On the Equilibrium of Heterogeneous Substances*) he treated $U(S, V, N)$ as a surface, introduced the potentials and the chemical potential, and derived the conditions for phase equilibrium, the phase rule, and the foundations of chemical thermodynamics. Maxwell was impressed enough to sculpt a plaster model of Gibbs's thermodynamic surface for water. Van der Waals's 1873 thesis produced the first equation of state with a liquid–gas transition, and Maxwell (1875) supplied the equal-area rule that fixes where the transition happens.

The forcing problem in one line: **energy conservation cannot explain why heat flows hot to cold, why engines have a maximum efficiency, or why anything is irreversible.** A cup of coffee cooling and a cup spontaneously heating from the room both conserve energy. Something else decides. That something is a second law, and its content is a new state function.

## 2. Intuition

**Entropy as a count.** SM-3 will show that $S = k_B\ln\Omega$, where $\Omega$ is the number of microscopic arrangements compatible with the macroscopic state. Take that on credit for now. Ten moles of gas have about $10^{25}$ particles, so $\Omega$ is something like $10^{10^{25}}$, and a change in macroscopic state that multiplies $\Omega$ by even $e^{100}$ makes the new state overwhelmingly more likely. The second law, "entropy of an isolated system does not decrease," is then the statement that systems wander into the macrostates with the most microstates and stay there.

**Why heat flows hot to cold.** Put two bodies in contact, sharing a fixed total energy. Moving energy $dU$ from body 1 to body 2 changes the total count by
$$d\ln\Omega_{\text{tot}} = \left(-\frac{\partial \ln\Omega_1}{\partial U_1} + \frac{\partial \ln\Omega_2}{\partial U_2}\right)dU.$$
Energy flows the way that increases the count, and stops when the two slopes are equal. The quantity that equalizes in thermal contact is what we call temperature, so the slope must be a function of temperature: $\partial S/\partial U = 1/T$. Why $1/T$ and not $T$? Because the slope must be larger for the colder body (energy flows *into* the colder one), and adding energy to a cold system multiplies its count of arrangements more than adding it to a hot one. Guess now: when a small heat $dQ$ leaves a body at $T_h$ and enters one at $T_c$, the total entropy changes by $dQ(1/T_c - 1/T_h)$, positive exactly when $T_h > T_c$.

**Why engines have a maximum efficiency.** An engine runs in a cycle, so its own entropy returns to where it started. It takes heat $Q_h$ from a hot reservoir, which loses entropy $Q_h/T_h$. If the total entropy is not to decrease, the engine must dump at least that much entropy somewhere, and the only place is the cold reservoir: $Q_c/T_c \ge Q_h/T_h$. Work carries no entropy. So the work $W = Q_h - Q_c$ is at most $Q_h(1 - T_c/T_h)$. You can guess Carnot's efficiency from bookkeeping alone. The real derivation (§3.1 to §3.4) runs the other way, from engines to entropy, because historically entropy had to be *found*.

**Free energy.** A system sitting in a room at temperature $T$ can trade heat with the room. If you want work out of it, some of its energy is locked up: to lower the system's entropy by $\Delta S$ you must hand the room at least $T\Delta S$ of heat. What is left to extract is $U - TS$, the Helmholtz free energy $F$. "Free" means available. A chemist at constant temperature *and* atmospheric pressure also has to push back the atmosphere, which adds $PV$: the Gibbs free energy $G = U - TS + PV$. Systems at fixed $T$ and $P$ roll downhill in $G$ the way a ball rolls downhill in potential energy.

**Legendre transforms.** $U(S, V)$ describes a system by its entropy and volume. But no lab instrument controls entropy; you control temperature by putting the system in a bath. The Legendre transform is the move "stop describing the system by its entropy and start describing it by the temperature of the bath it is in." The curve $U(S)$ can be described by its points or by its tangent lines; a tangent line has slope $T = \partial U/\partial S$ and intercept $U - TS = F$. Describing the curve by "intercept as a function of slope" loses nothing if the curve is convex. CM-8 does the same thing to the Lagrangian with velocity traded for momentum.

**Guess the Maxwell relations.** If $F(T, V)$ is a function with $\partial F/\partial T = -S$ and $\partial F/\partial V = -P$, what does equality of mixed partials $\partial^2F/\partial V\partial T = \partial^2 F/\partial T\partial V$ give you? Try it before §3.9.

## 3. Formal development

### 3.1 The Carnot cycle

A **heat engine** is a device that runs in a cycle, takes heat $Q_h$ from a hot reservoir at $T_h$, rejects $Q_c$ to a cold reservoir at $T_c$, and delivers work $W$. Over a cycle $\Delta U = 0$, so $W = Q_h - Q_c$ and the **efficiency** is
$$\eta = \frac{W}{Q_h} = 1 - \frac{Q_c}{Q_h}. \tag{3.1}$$
A **reservoir** is a body so large that heat flowing in or out does not change its temperature.

The **Carnot cycle** uses only the two reservoirs and no other heat exchange, so heat can only enter or leave isothermally at $T_h$ or $T_c$, and every other leg must be adiabatic. For an ideal gas the four legs are: (1→2) isothermal expansion at $T_h$ from $V_1$ to $V_2$; (2→3) adiabatic expansion to $V_3$, cooling to $T_c$; (3→4) isothermal compression at $T_c$ to $V_4$; (4→1) adiabatic compression back to $V_1$, warming to $T_h$.

▶ **Try it first.** Use SM-1's two results: on an isotherm $\Delta U = 0$ so $Q = W = Nk_BT\ln(V_f/V_i)$; on an adiabat $TV^{\gamma-1}$ is constant. Compute $Q_c/Q_h$.

**Derivation.**

1. Leg 1→2: $\Delta U = 0$, so $Q_h = W_{12} = Nk_BT_h\ln(V_2/V_1)$.
2. Leg 3→4: heat leaves the gas, $Q_c = Nk_BT_c\ln(V_3/V_4)$ (positive, since $V_3 > V_4$).
3. Adiabats: $T_hV_2^{\gamma-1} = T_cV_3^{\gamma-1}$ and $T_hV_1^{\gamma-1} = T_cV_4^{\gamma-1}$.
4. Divide the two adiabat equations: $(V_2/V_1)^{\gamma-1} = (V_3/V_4)^{\gamma-1}$, so $V_3/V_4 = V_2/V_1$.
5. Therefore the logarithms in steps 1 and 2 are equal, and
$$\frac{Q_c}{Q_h} = \frac{T_c}{T_h}, \tag{3.2}$$
$$\eta_{\text{Carnot}} = 1 - \frac{T_c}{T_h}. \tag{3.3}$$
∎

■ **What just happened.** The ratio of heats equals the ratio of temperatures; the volumes and $\gamma$ dropped out. Every step was reversible: run the cycle backwards and it becomes a refrigerator with the same $Q$s and $W$ reversed. For now (3.3) is a fact about one ideal-gas engine. §3.3 shows it is a fact about every reversible engine.

### 3.2 The second law: Kelvin and Clausius

**Kelvin statement.** No process is possible whose *sole* result is the extraction of heat from a single reservoir and its complete conversion into work.

**Clausius statement.** No process is possible whose *sole* result is the transfer of heat from a colder body to a hotter one.

"Sole" is the load-bearing word. An isothermal ideal-gas expansion converts heat from one reservoir entirely into work, but the gas ends up in a different state (larger volume), so that is not the sole result. A refrigerator moves heat from cold to hot, but uses work. The laws forbid doing these things *with nothing else changing*.

**Theorem 3.1 (Equivalence).** The Kelvin and Clausius statements are equivalent.

▶ **Try it first.** Prove the contrapositives: a device that violates one statement, coupled to an ordinary engine or refrigerator, gives a device that violates the other.

**Proof.**

1. *Not Kelvin ⇒ not Clausius.* Suppose device $K$ takes heat $Q$ from a hot reservoir and converts it all to work $W = Q$. Use $W$ to drive an ordinary refrigerator between the cold and hot reservoirs: it extracts $Q_c$ from the cold one and delivers $Q_c + W$ to the hot one. Net effect of $K$ plus refrigerator: the hot reservoir receives $Q_c + W - Q = Q_c$, the cold one loses $Q_c$, no net work. Heat $Q_c$ has moved from cold to hot as the sole result. Clausius is violated.
2. *Not Clausius ⇒ not Kelvin.* Suppose device $C$ moves heat $Q_c$ from cold to hot with no other effect. Run an ordinary engine between the same reservoirs, adjusted so it takes $Q_h$ from the hot one and rejects exactly $Q_c$ to the cold one, producing $W = Q_h - Q_c > 0$. Net effect: the cold reservoir is unchanged, the hot reservoir has lost $Q_h - Q_c$, and that amount is now work. Heat from a single reservoir has been converted entirely into work. Kelvin is violated. ∎

■ **What just happened.** Two statements that sound different (one about engines, one about refrigerators) forbid the same thing. From here on "the second law" means either. Textbooks list a dozen more forms (Carathéodory's, "entropy of an isolated system never decreases," "no engine beats Carnot"). They are all one law; §8 has advice on not treating them as twelve.

### 3.3 Carnot's theorem

**Theorem 3.2 (Carnot).** (a) No engine operating between two reservoirs is more efficient than a reversible engine operating between the same two reservoirs. (b) All reversible engines operating between the same two reservoirs have the same efficiency.

▶ **Try it first.** Assume an engine $X$ beats a reversible engine $R$. Run $R$ backwards as a refrigerator, sized so that it returns to the hot reservoir exactly the heat $X$ takes from it. What is the net effect?

**Proof.**

1. Let $X$ have efficiency $\eta_X$ and $R$ be reversible with efficiency $\eta_R$, and suppose $\eta_X > \eta_R$.
2. Run $X$ forward: it takes $Q_h$ from the hot reservoir and produces $W_X = \eta_XQ_h$, rejecting $Q_h - W_X$ to the cold one.
3. Run $R$ backward (allowed because it is reversible). Size it so it delivers exactly $Q_h$ to the hot reservoir. Reversed, every heat and work flips sign, so it needs work input $W_R = \eta_RQ_h$ and takes $Q_h - W_R$ from the cold reservoir.
4. Net effect per cycle of the pair: the hot reservoir is unchanged. Net work produced: $W_X - W_R = (\eta_X - \eta_R)Q_h > 0$. Net heat taken from the cold reservoir: $(Q_h - W_R) - (Q_h - W_X) = W_X - W_R$.
5. So heat from a single reservoir (the cold one) has been converted entirely into work, with nothing else changing. That violates Kelvin. Hence $\eta_X \le \eta_R$, which is (a).
6. If $X$ is also reversible, swap the roles of $X$ and $R$: $\eta_R \le \eta_X$. Hence $\eta_X = \eta_R$, which is (b). ∎

■ **What just happened.** The efficiency of a reversible engine cannot depend on the working substance, the volumes, or the engine design: only on the two reservoirs, i.e. on their temperatures. We computed it once, for the ideal gas, in (3.3), so $1 - T_c/T_h$ is the efficiency of *every* reversible engine. Kelvin (1848) turned this around: since $Q_c/Q_h$ for a reversible engine is a universal function of the two reservoirs, it can *define* a temperature ratio without reference to any thermometer substance. That is the **thermodynamic temperature**, and (3.2) says the ideal-gas temperature agrees with it.

For the next step write the heats as signed heats *into the engine*, $Q_1 = Q_h$ at $T_1 = T_h$ and $Q_2 = -Q_c$ at $T_2 = T_c$. Then (3.2) reads, for every reversible two-reservoir cycle,
$$\frac{Q_1}{T_1} + \frac{Q_2}{T_2} = 0. \tag{3.4}$$
For an irreversible engine, part (a) gives $1 + Q_2/Q_1 \le 1 - T_2/T_1$, so $Q_1/T_1 + Q_2/T_2 \le 0$.

### 3.4 The Clausius inequality and the existence of entropy

**Theorem 3.3 (Clausius inequality).** For any cyclic process in which a system receives heats $đQ$ from reservoirs at temperatures $T$,
$$\oint \frac{đQ}{T} \le 0, \tag{3.5}$$
with equality if the cycle is reversible. Here $T$ is the temperature of the reservoir supplying $đQ$. During an irreversible process the system itself may not even have a temperature.

▶ **Try it first.** Replace each reservoir the system touches with a small reversible Carnot engine fed from one master reservoir at $T_0$. Use (3.4) for each small engine, then apply Kelvin to everything taken together.

**Derivation.**

1. Split the cycle into steps $i = 1, \dots, n$; in step $i$ the system receives heat $đQ_i$ (signed) from a reservoir at $T_i$.
2. Introduce a master reservoir at $T_0$. For each $i$, place a reversible Carnot engine $C_i$ between $T_0$ and $T_i$, run so that it delivers exactly $đQ_i$ to reservoir $i$ during step $i$. Reservoir $i$ then ends with no net change.
3. By (3.4) applied to $C_i$ (heat $đQ_i^{(0)}$ in from $T_0$, heat $đQ_i$ out to $T_i$): $đQ_i^{(0)}/T_0 - đQ_i/T_i = 0$, so $đQ_i^{(0)} = T_0\,đQ_i/T_i$.
4. After one full cycle, the system is back in its initial state, each $C_i$ has completed whole cycles, and each reservoir $i$ is unchanged. The only net effect: heat $Q_0 = \sum_i đQ_i^{(0)} = T_0\sum_i đQ_i/T_i$ has left the master reservoir, and by the first law (everything else is cyclic) it has all become work.
5. Kelvin forbids $Q_0 > 0$. So $T_0\sum_i đQ_i/T_i \le 0$, and since $T_0 > 0$, $\sum_i đQ_i/T_i \le 0$. In the limit of infinitely many steps this is (3.5).
6. If the cycle is reversible, run it backwards: every $đQ_i$ flips sign, and step 5 applied to the reversed cycle gives $-\oint đQ/T \le 0$. Both inequalities hold, so $\oint đQ_{\text{rev}}/T = 0$. ∎

**Corollary 3.4 (Entropy exists).** There is a state function $S$ with
$$dS = \frac{đQ_{\text{rev}}}{T}, \qquad S(B) - S(A) = \int_A^B \frac{đQ_{\text{rev}}}{T}\ \text{along any reversible path}. \tag{3.6}$$

*Proof.* Take two reversible paths $\mathcal{P}_1, \mathcal{P}_2$ from $A$ to $B$. Going out along $\mathcal{P}_1$ and back along $\mathcal{P}_2$ reversed is a reversible cycle, so by step 6 the integral of $đQ_{\text{rev}}/T$ around it vanishes, i.e. $\int_{\mathcal{P}_1} = \int_{\mathcal{P}_2}$. An integral that depends only on its endpoints defines a function of state, exactly as path-independence of $\int\mathbf F\cdot d\mathbf r$ defines a potential. Fix $S(A)$ arbitrarily; this leaves an additive constant, which SM-3 will fix. ∎

**Corollary 3.5 (Irreversible processes).** For any process from $A$ to $B$,
$$S(B) - S(A) \ge \int_A^B \frac{đQ}{T}, \tag{3.7}$$
with $T$ the reservoir temperature. In particular, for a thermally isolated system $\Delta S \ge 0$; and since the system plus every reservoir it touches form an isolated whole,
$$\Delta S_{\text{universe}} \ge 0. \tag{3.8}$$

*Proof.* Go from $A$ to $B$ along the actual (possibly irreversible) path and return along a reversible one. By (3.5), $\int_{A\to B}đQ/T + \int_{B\to A}đQ_{\text{rev}}/T \le 0$. The second integral is $S(A) - S(B)$ by (3.6). Rearrange. For an isolated system $đQ = 0$ throughout. ∎

■ **What just happened.** Carnot's engine theorem, applied to a fleet of tiny engines, produced a new state function. $1/T$ is an integrating factor: $đQ_{\text{rev}}$ is inexact (SM-1), but $đQ_{\text{rev}}/T$ is exact. The second law in its most useful form is (3.8), and it makes no reference to engines at all.

### 3.5 Entropy changes you should be able to compute

The recipe is always the same: $S$ is a state function, so to find $\Delta S$ for *any* process, real or idealized, invent a *reversible* path between the same end states and integrate $đQ_{\text{rev}}/T$ along it.

**Heating at constant volume.** $đQ_{\text{rev}} = C_V\,dT$, so
$$\Delta S = \int_{T_1}^{T_2}\frac{C_V\,dT}{T} = C_V\ln\frac{T_2}{T_1}\quad(C_V\text{ constant}). \tag{3.9}$$
Same with $C_P$ at constant pressure.

**Isothermal expansion of an ideal gas.** $\Delta U = 0$, so $đQ_{\text{rev}} = P\,dV = Nk_BT\,dV/V$, and
$$\Delta S = Nk_B\ln\frac{V_f}{V_i}. \tag{3.10}$$
The reservoir supplies $Q = T\Delta S$ and loses entropy $Q/T$; for the reversible process $\Delta S_{\text{universe}} = 0$.

**Free expansion.** A gas in an insulated box expands into vacuum when a partition is removed. $Q = 0$ and $W = 0$, so $\Delta U = 0$, and for an ideal gas $T$ is unchanged. The end states are the same as for the isothermal expansion, so $\Delta S$ is the same, (3.10). But now there is no reservoir and no heat. **Entropy increased with no heat flow at all.** That is not a contradiction of $dS = đQ_{\text{rev}}/T$: the actual process is irreversible, so its $đQ = 0$ is not $đQ_{\text{rev}}$, and (3.7) says $\Delta S \ge 0$, which is what we got. This is the cleanest test of whether you understand that $S$ belongs to the state, not to the process.

**Heat flow between two bodies.** Bodies with constant heat capacities $C_1$, $C_2$ at $T_1 > T_2$ are put in contact and isolated from everything else.

▶ **Try it first.** Show $\Delta S_{\text{total}} > 0$ without computing the final temperature, by looking at each small heat transfer.

**Derivation.**

1. At some instant the bodies are at $T_1' > T_2'$ (heat flows from 1 to 2 until they are equal, so this ordering persists). A small heat $dQ > 0$ leaves body 1 and enters body 2.
2. Each body separately undergoes a change it could undergo reversibly (exchange heat at its own temperature), so $dS_1 = -dQ/T_1'$ and $dS_2 = +dQ/T_2'$.
3. Total:
$$dS_{\text{tot}} = dQ\left(\frac{1}{T_2'} - \frac{1}{T_1'}\right) > 0. \tag{3.11}$$
4. Every increment is positive, so the integrated $\Delta S_{\text{tot}} > 0$. Explicitly, energy conservation gives $T_f = (C_1T_1 + C_2T_2)/(C_1+C_2)$ and $\Delta S_{\text{tot}} = C_1\ln(T_f/T_1) + C_2\ln(T_f/T_2)$. For $C_1 = C_2 = C$ this is $C\ln\big(T_f^2/T_1T_2\big) > 0$ because the arithmetic mean $T_f$ exceeds the geometric mean $\sqrt{T_1T_2}$. ∎

■ **What just happened.** Clausius's statement and "entropy increases" are one statement: heat flowing hot to cold increases total entropy, cold to hot would decrease it. The project's Part 3 plots $S_1 + S_2$ against the energy split and shows the maximum sits at equal temperatures.

**Mixing.** $N_A$ molecules of gas A in volume $V_A$ and $N_B$ of gas B in $V_B$, same $T$ and $P$, partition removed. Ideal gases do not interact, so each gas undergoes its own free expansion into $V = V_A + V_B$:
$$\Delta S_{\text{mix}} = N_Ak_B\ln\frac{V}{V_A} + N_Bk_B\ln\frac{V}{V_B} = -Nk_B\sum_i x_i\ln x_i, \tag{3.12}$$
with $x_i = N_i/N = V_i/V$ (equal $T$, $P$). For one mole each, $\Delta S = 2R\ln 2 = 11.53$ J/K.

*The Gibbs paradox.* If A and B are the *same* gas, removing the partition changes nothing you can measure, and sliding it back restores the initial state. So $\Delta S$ must be zero. But (3.12) does not know how different A and B are; it gives $2R\ln2$ for nitrogen and oxygen and seems to give it for nitrogen and nitrogen. The resolution needs entropy to scale as $\ln(V/N)$, not $\ln V$, for identical particles; the classical argument here supplies it only by demanding extensivity (see (3.16) below), and SM-3 explains it by counting indistinguishable particles.

**Melting.** At the melting point, at constant $T_m$ and $P$, the latent heat $L$ is absorbed reversibly, so
$$\Delta S = \frac{L}{T_m}. \tag{3.13}$$
For 1 kg of ice ($L = 334$ kJ/kg) at 273.15 K: $\Delta S = 1222.8$ J/K.

### 3.6 The thermodynamic identity, and the entropy of an ideal gas

For a reversible change of a simple system, $đQ = T\,dS$ and $đW = P\,dV$, so the first law becomes
$$dU = T\,dS - P\,dV + \mu\,dN, \tag{3.14}$$
where the last term (allowing particle number to change) defines the chemical potential $\mu$ (§3.10). Every term in (3.14) is a differential of state functions, so it holds between any two nearby equilibrium states however the system got there. What does *not* hold for an irreversible change is the identification $T\,dS = đQ$. Reading off the partial derivatives of $U(S, V, N)$:
$$T = \left(\frac{\partial U}{\partial S}\right)_{V,N}, \qquad P = -\left(\frac{\partial U}{\partial V}\right)_{S,N}, \qquad \mu = \left(\frac{\partial U}{\partial N}\right)_{S,V}. \tag{3.15}$$
Solved for $dS$: $dS = dU/T + (P/T)\,dV - (\mu/T)\,dN$, so $1/T = (\partial S/\partial U)_{V,N}$, the relation guessed in §2.

**Entropy of an ideal gas.**

▶ **Try it first.** Write $1/T$ and $P/T$ as functions of $U$ and $V$ for the ideal gas, integrate $dS$, and then ask what extensivity demands of the $N$-dependent "constant".

**Derivation.**

1. From $U = \tfrac f2Nk_BT$: $1/T = \tfrac f2 Nk_B/U$. From $PV = Nk_BT$: $P/T = Nk_B/V$.
2. At fixed $N$: $dS = \tfrac f2 Nk_B\,dU/U + Nk_B\,dV/V$, which integrates to $S = Nk_B\big[\ln V + \tfrac f2\ln U\big] + \phi(N)$.
3. Entropy is extensive: doubling the system ($U, V, N \to \lambda U, \lambda V, \lambda N$) multiplies $S$ by $\lambda$. Substituting, this requires $\phi(\lambda N) = \lambda\phi(N) - \lambda Nk_B(1 + \tfrac f2)\ln\lambda$, which is satisfied by $\phi(N) = Nk_B\big[\sigma - (1 + \tfrac f2)\ln N\big]$ with $\sigma$ a constant. (Check: $\phi(\lambda N) = \lambda Nk_B[\sigma - (1+\tfrac f2)(\ln\lambda + \ln N)]$.)
4. Hence
$$S(U, V, N) = Nk_B\left[\ln\frac{V}{N} + \frac f2\ln\frac{U}{N} + \sigma\right]. \tag{3.16}$$
∎

■ **What just happened.** Thermodynamics determines $S$ up to the constant $\sigma$. SM-3 computes $\sigma$ by counting quantum states (the Sackur–Tetrode equation); it involves Planck's constant, which no amount of classical thermodynamics could produce. The $V/N$ in (3.16) is exactly what makes the Gibbs paradox go away: mixing a gas with itself leaves $V/N$ unchanged.

**Third law.** Nernst (1906) found experimentally that entropy changes in reactions vanish as $T \to 0$; Planck (1911) sharpened this to: the entropy of a perfect crystal tends to zero as $T \to 0$. A consequence you can check: since $S(T) - S(0) = \int_0^T C(T')\,dT'/T'$ must be finite, heat capacities must vanish as $T \to 0$. Equation (3.16) has $S \to -\infty$ as $U \to 0$, which says the classical ideal gas fails at low temperature. SM-5 (quantum gases) is where it is replaced.

### 3.7 Engines, refrigerators, heat pumps

The entropy argument of §2 now runs in full rigor, and gives all three bounds at once. A cyclic device exchanges $Q_h$ with $T_h$ and $Q_c$ with $T_c$, and its own $\Delta S$ per cycle is zero, so by (3.8) the reservoirs' entropy changes must sum to $\ge 0$.

- **Engine** ($Q_h$ out of hot, $Q_c$ into cold): $-Q_h/T_h + Q_c/T_c \ge 0$, so $\eta = 1 - Q_c/Q_h \le 1 - T_c/T_h$, which is (3.3) again.
- **Refrigerator** (extract $Q_c$ from the cold interior using work $W$, dump $Q_h = Q_c + W$ to the room): $Q_h/T_h - Q_c/T_c \ge 0$, so the coefficient of performance is
$$\mathrm{COP}_{\text{fridge}} = \frac{Q_c}{W} \le \frac{T_c}{T_h - T_c}. \tag{3.17}$$
- **Heat pump** (same machine; what you want is the $Q_h$ delivered indoors):
$$\mathrm{COP}_{\text{hp}} = \frac{Q_h}{W} \le \frac{T_h}{T_h - T_c}. \tag{3.18}$$

COPs exceed 1 and that is not a violation of anything: the work is not turned into the heat, it *moves* heat that was already there. A heat pump keeping a house at 20 °C with the outdoor air at 0 °C has a Carnot COP of $293.15/20 = 14.66$. Real units achieve about 3 to 4 in those conditions. A resistive heater has COP exactly 1. So a real heat pump delivers three or four times the heat per joule of electricity, which is why they beat resistive heating even after you account for how the electricity was made (Exercise C3).

**The Otto cycle** (the idealized gasoline engine; SM-1 set it up and stopped at the heats). Air-standard: (1→2) adiabatic compression from $V_1$ to $V_1/r$, where $r$ is the compression ratio; (2→3) heat input at constant volume (combustion, fast); (3→4) adiabatic expansion back to $V_1$; (4→1) heat rejection at constant volume (exhaust).

**Derivation.**

1. Heat in: $Q_h = C_V(T_3 - T_2)$. Heat out: $Q_c = C_V(T_4 - T_1)$. So $\eta = 1 - (T_4 - T_1)/(T_3 - T_2)$.
2. Adiabats with volume ratio $r$: $T_2 = T_1r^{\gamma-1}$ and $T_3 = T_4r^{\gamma-1}$.
3. Subtract: $T_3 - T_2 = r^{\gamma-1}(T_4 - T_1)$. Therefore
$$\eta_{\text{Otto}} = 1 - r^{1-\gamma}. \tag{3.19}$$
∎

For $r = 10$, $\gamma = 1.4$: $\eta = 0.602$. That is below Carnot between the hottest and coldest temperatures of the cycle, because heat enters over a range of temperatures, not all at $T_3$.

**The Diesel cycle** replaces the constant-volume heat input with a constant-pressure one (fuel injected and burned as the piston moves out) from $V_2$ to $V_3 = r_cV_2$; $r_c$ is the cutoff ratio, $r = V_1/V_2$ the compression ratio.

**Derivation.**

1. Heat in at constant pressure: $Q_h = C_P(T_3 - T_2)$. Heat out at constant volume: $Q_c = C_V(T_4 - T_1)$. So $\eta = 1 - (T_4 - T_1)/[\gamma(T_3 - T_2)]$.
2. Adiabatic compression: $T_2 = T_1r^{\gamma-1}$. Isobaric: $V \propto T$, so $T_3 = r_cT_2 = r_cr^{\gamma-1}T_1$.
3. Adiabatic expansion from $V_3$ to $V_1$, volume ratio $V_1/V_3 = r/r_c$: $T_4 = T_3(r_c/r)^{\gamma-1} = r_cr^{\gamma-1}T_1\,r_c^{\gamma-1}r^{1-\gamma} = r_c^{\gamma}T_1$.
4. Substitute:
$$\eta_{\text{Diesel}} = 1 - \frac{1}{r^{\gamma-1}}\,\frac{r_c^{\gamma} - 1}{\gamma(r_c - 1)}. \tag{3.20}$$
∎

For $r = 18$, $r_c = 2$: $\eta = 0.632$. At the same compression ratio Otto would give $0.685$; the factor $(r_c^\gamma - 1)/[\gamma(r_c - 1)] > 1$ is the Diesel penalty. Diesels win in practice because they can run at $r \approx 15$ to $20$ without knock (the fuel is injected after compression), while gasoline engines are limited to about $r \approx 10$.

**Real-engine losses.** A real car engine converts roughly 25–35% of the fuel's energy to shaft work at a good operating point, well below (3.19). The losses: heat conducted through cylinder walls during combustion and expansion; combustion that takes finite time rather than happening at constant volume; friction; pumping losses moving air in and exhaust out; exhaust leaving hotter than ambient; and the fact that the working gas is not air with constant $\gamma$. Every one of these is an irreversibility, and every one generates entropy. The project's particle engine shows the simplest version: moving the piston ten times faster costs efficiency.

### 3.8 Thermodynamic potentials as Legendre transforms

**Convexity first.** CM-8 Definition 3.1 defines the Legendre transform of a convex $f$ as $g(p) = \max_v[pv - f(v)]$, with $p = f'(v)$ and $v = g'(p)$, and shows nothing is lost passing from $f$ to $g$ (Theorem 3.1 there). It requires convexity. So we need: is $U(S)$ convex at fixed $V$, $N$?

**Theorem 3.6.** At fixed $V, N$: $S(U)$ is concave $\iff$ $U(S)$ is convex $\iff$ $C_V > 0$. And equilibrium stability requires it.

▶ **Try it first.** Compute $\partial^2 S/\partial U^2$ using $\partial S/\partial U = 1/T$ and $C_V = (\partial U/\partial T)_V$.

**Derivation.**

1. $\dfrac{\partial^2 S}{\partial U^2} = \dfrac{\partial}{\partial U}\dfrac1T = -\dfrac{1}{T^2}\dfrac{\partial T}{\partial U} = -\dfrac{1}{T^2C_V}$. $\tag{3.21}$
   So $S(U)$ is concave exactly when $C_V > 0$.
2. Likewise $\partial^2U/\partial S^2 = (\partial T/\partial S)_V = T/C_V$ (using $C_V = T(\partial S/\partial T)_V$), positive exactly when $C_V > 0$.
3. *Why stability demands it.* Split a system into two equal halves, each with energy $U/2$. A fluctuation moves $\delta$ from one half to the other. The total entropy becomes $S_{1/2}(U/2 + \delta) + S_{1/2}(U/2 - \delta) = 2S_{1/2}(U/2) + S_{1/2}''\,\delta^2 + O(\delta^4)$. If $S'' > 0$, the fluctuation *raises* the entropy and grows: energy would pile up on one side spontaneously, and the uniform state would not be an equilibrium. So a stable uniform equilibrium requires $S'' < 0$, i.e. $C_V > 0$. ∎

■ **What just happened.** Positive heat capacity is not an empirical accident; it is the condition for the entropy maximum to be a maximum. (Self-gravitating systems such as stars have negative heat capacity, and they are not in thermodynamic equilibrium in this sense.) Convexity is also exactly what licenses the Legendre transform.

**The construction.** The sign convention in thermodynamics is reversed from CM-8: define
$$F(T) = \min_S\big[U(S) - TS\big] = -\max_S\big[TS - U(S)\big],$$
so $F = -g$ with $g$ the CM-8 transform of $U$. The minimizer satisfies $T = \partial U/\partial S$ (3.15). The tangent line to $U(S)$ of slope $T$ meets the $U$-axis at $U - TS = F$: the free energy is the intercept. Project Part 5 does the minimization numerically on a table of $U(S)$ and recovers the closed form.

▶ **Try it first.** Take the differentials of $U + PV$, $U - TS$, and $U - TS + PV$ and substitute (3.14). Which differentials survive?

**Derivation of the three potentials.**

1. **Enthalpy**, transform of $U$ in $V$: $H = U + PV$ (the minimizer of $U + PV$ over $V$ satisfies $\partial U/\partial V = -P$). Then $dH = dU + P\,dV + V\,dP = T\,dS + V\,dP + \mu\,dN$:
$$H(S, P, N),\qquad dH = T\,dS + V\,dP + \mu\,dN. \tag{3.22}$$
2. **Helmholtz free energy**, transform in $S$: $F = U - TS$. $dF = dU - T\,dS - S\,dT$:
$$F(T, V, N),\qquad dF = -S\,dT - P\,dV + \mu\,dN. \tag{3.23}$$
3. **Gibbs free energy**, transform in both: $G = U - TS + PV$. $dG = dU - T\,dS - S\,dT + P\,dV + V\,dP$:
$$G(T, P, N),\qquad dG = -S\,dT + V\,dP + \mu\,dN. \tag{3.24}$$
∎

■ **What just happened.** Each transform cancels one term of (3.14) and replaces it with its conjugate: $T\,dS \to -S\,dT$, $-P\,dV \to V\,dP$. The variable whose differential appears is the *natural variable*; in it, the potential contains everything about the system, because its partial derivatives return all the others ($S = -(\partial F/\partial T)_V$, $P = -(\partial F/\partial V)_T$). Written in the wrong variables, say $F(T, S)$, information is lost. This is CM-8's (3.3), $g'(p) = v$, four times over.

**Why each one is useful.**

**$F$ is the maximum work at constant temperature.** A system in contact with a reservoir at $T$ goes from equilibrium state $A$ to equilibrium state $B$, both at temperature $T$, doing work $W$ of any kind.

1. First law: $\Delta U = Q - W$, where $Q$ is the heat from the reservoir.
2. Second law (3.8): $\Delta S + \Delta S_{\text{res}} \ge 0$ with $\Delta S_{\text{res}} = -Q/T$. So $Q \le T\Delta S$.
3. Therefore $W = Q - \Delta U \le T\Delta S - \Delta U = -\Delta(U - TS)$, using that $T$ is the same at $A$ and $B$:
$$W \le -\Delta F, \tag{3.25}$$
with equality for a reversible process. ∎

The maximum work you can extract is the decrease of $F$, not of $U$; the difference $T\Delta S$ must go to the reservoir as heat. If the volume is fixed and no work is done, (3.25) says $\Delta F \le 0$: at fixed $T$ and $V$, a system evolves to minimize $F$.

**$G$ is minimized at constant temperature and pressure.** This is the argument that makes chemistry and phase transitions go. A system exchanges heat with, and pushes against, a reservoir (the atmosphere) at fixed $T$ and $P$.

1. The reservoir loses heat $Q$ to the system and has volume $\Delta V$ taken from it. Its energy changes by $\Delta U_{\text{res}} = -Q + P\Delta V$ (the system pushes back the atmosphere, doing work on it). If the system does no other work, $\Delta U = Q - P\Delta V$, so $\Delta U_{\text{res}} = -\Delta U$.
2. The reservoir is internally in equilibrium at fixed $T$, $P$, so (3.14) applies to it: $\Delta S_{\text{res}} = (\Delta U_{\text{res}} + P\Delta V_{\text{res}})/T = (-\Delta U - P\Delta V)/T$, since $\Delta V_{\text{res}} = -\Delta V$.
3. Total entropy change:
$$\Delta S_{\text{tot}} = \Delta S - \frac{\Delta U + P\Delta V}{T} = -\frac{\Delta(U - TS + PV)}{T} = -\frac{\Delta G}{T} \ge 0. \tag{3.26}$$
∎

■ **What just happened.** "Total entropy of system plus surroundings increases" has been rewritten entirely in terms of the system: $\Delta G \le 0$. You never need to track the atmosphere's entropy again. Every chemical equilibrium constant and every phase boundary in §3.11 comes from minimizing $G$. With other work allowed, the same steps give $W_{\text{other}} \le -\Delta G$: the maximum electrical work from a battery is $-\Delta G$ of its reaction.

**$H$ is the heat at constant pressure.** At constant $P$ with only $P\,dV$ work, $Q = \Delta U + P\Delta V = \Delta H$. Calorimetry at atmospheric pressure measures enthalpy changes; latent heats and heats of reaction are tabulated as $\Delta H$; and $C_P = (\partial H/\partial T)_P$. Throttling (below) conserves $H$.

### 3.9 Maxwell relations and what they buy

**Theorem 3.7 (Maxwell relations).**
$$\left(\frac{\partial T}{\partial V}\right)_S = -\left(\frac{\partial P}{\partial S}\right)_V,\quad \left(\frac{\partial T}{\partial P}\right)_S = \left(\frac{\partial V}{\partial S}\right)_P,\quad \left(\frac{\partial S}{\partial V}\right)_T = \left(\frac{\partial P}{\partial T}\right)_V,\quad \left(\frac{\partial S}{\partial P}\right)_T = -\left(\frac{\partial V}{\partial T}\right)_P. \tag{3.27}$$

▶ **Try it first.** You did this at the end of §2 for $F$. Do it for $U$, $H$, $G$.

**Derivation.** For a smooth function $\Phi(x, y)$, $\partial^2\Phi/\partial y\partial x = \partial^2\Phi/\partial x\partial y$ (Clairaut).

1. $U(S, V)$: $\partial U/\partial S = T$, $\partial U/\partial V = -P$. Mixed partials: $(\partial T/\partial V)_S = -(\partial P/\partial S)_V$.
2. $H(S, P)$: $\partial H/\partial S = T$, $\partial H/\partial P = V$. So $(\partial T/\partial P)_S = (\partial V/\partial S)_P$.
3. $F(T, V)$: $\partial F/\partial T = -S$, $\partial F/\partial V = -P$. So $-(\partial S/\partial V)_T = -(\partial P/\partial T)_V$.
4. $G(T, P)$: $\partial G/\partial T = -S$, $\partial G/\partial P = V$. So $-(\partial S/\partial P)_T = (\partial V/\partial T)_P$. ∎

■ **What just happened.** Four relations, each one line, each from one potential. The useful ones are the last two: they convert derivatives of entropy, which you cannot measure directly, into derivatives of $P$, $V$, $T$, which you can (the equation of state). There is a mnemonic (the Born square: write V F T across the top, U G down the middle, S H P across the bottom, potentials flanked by their natural variables), but deriving a relation from the right potential takes ten seconds and never gives you a wrong sign. Derive; do not memorize.

Two calculus tools we need. With $\alpha$ and $\kappa_T$ as in the conventions box:

- **Cyclic rule.** For three variables related by one equation, $(\partial x/\partial y)_z(\partial y/\partial z)_x(\partial z/\partial x)_y = -1$. Proof: write $dx = (\partial x/\partial y)_z\,dy + (\partial x/\partial z)_y\,dz$ and set $dx = 0$; this gives $(\partial y/\partial z)_x = -(\partial x/\partial z)_y/(\partial x/\partial y)_z$. Applied to $P, V, T$: $(\partial P/\partial T)_V = -(\partial V/\partial T)_P/(\partial V/\partial P)_T = \alpha/\kappa_T$.
- **Changing what is held fixed.** For $S(T, V)$: $(\partial S/\partial T)_P = (\partial S/\partial T)_V + (\partial S/\partial V)_T(\partial V/\partial T)_P$ (chain rule along a path at constant $P$).

**Application 1: entropy and energy of a van der Waals gas.** Per mole, $P = RT/(v - b) - a/v^2$ (SM-1), where $a$ measures attraction and $b$ the excluded volume.

1. Maxwell (third of (3.27)): $(\partial s/\partial v)_T = (\partial P/\partial T)_v = R/(v - b)$.
2. Is $c_v$ a function of $T$ only? $(\partial c_v/\partial v)_T = T\,\partial^2 s/\partial v\,\partial T = T(\partial^2P/\partial T^2)_v = 0$, since $P$ is linear in $T$.
3. Integrating, with $c_v$ taken constant:
$$s(T, v) = c_v\ln T + R\ln(v - b) + \text{const}. \tag{3.28}$$
4. Energy: from $du = T\,ds - P\,dv$, $(\partial u/\partial v)_T = T(\partial s/\partial v)_T - P = RT/(v - b) - P = a/v^2$. So $u = c_vT - a/v$: attractions lower the energy, and a vdW gas cools on free expansion (Exercise B4). ∎

We got the entropy of a non-ideal gas from its equation of state and one heat capacity, without any microscopic model. That is the Maxwell relations' main job.

**Application 2: $C_P - C_V$ in general.** SM-1 showed $C_P - C_V = Nk_B$ for the ideal gas. Here is the general result, and it is the showpiece calculation of the module.

▶ **Try it first.** Multiply the "changing what is held fixed" identity by $T$, then use a Maxwell relation and the cyclic rule.

**Derivation.**

1. $C_P = T(\partial S/\partial T)_P$ and $C_V = T(\partial S/\partial T)_V$. Multiply the identity for $(\partial S/\partial T)_P$ by $T$:
   $C_P = C_V + T(\partial S/\partial V)_T(\partial V/\partial T)_P$.
2. Maxwell: $(\partial S/\partial V)_T = (\partial P/\partial T)_V$.
3. Cyclic rule: $(\partial P/\partial T)_V = \alpha/\kappa_T$. And $(\partial V/\partial T)_P = V\alpha$.
4. Therefore
$$C_P - C_V = \frac{TV\alpha^2}{\kappa_T} = -T\frac{(\partial P/\partial T)_V^2}{(\partial P/\partial V)_T}. \tag{3.29}$$
∎

■ **What just happened.** Every quantity on the right is measured from the equation of state. Checks: for the ideal gas $\alpha = 1/T$, $\kappa_T = 1/P$, giving $TV/(T^2)\cdot P = PV/T = Nk_B$. Consequences: since $\kappa_T > 0$ for any stable substance (the volume analogue of Theorem 3.6), $C_P \ge C_V$ always; and $C_P = C_V$ exactly where $\alpha = 0$, e.g. water at 4 °C. For solids and liquids $C_V$ is hard to measure (you would have to stop the sample from expanding); (3.29) gets it from $C_P$. Example 4.5 does it for water; project Part 4 checks (3.29) symbolically for the vdW gas.

**Application 3: the Joule–Thomson effect.** Gas is pushed steadily through a porous plug or valve from pressure $P_1$ to $P_2 < P_1$, in an insulated pipe.

1. Follow a parcel of gas of volume $V_1$ upstream, $V_2$ downstream. The upstream gas pushes it in, doing work $P_1V_1$ on it; it pushes the downstream gas, doing work $P_2V_2$. No heat flows. So $U_2 - U_1 = P_1V_1 - P_2V_2$, i.e. $H_1 = H_2$. **Throttling conserves enthalpy.**
2. The Joule–Thomson coefficient is $\mu_{JT} = (\partial T/\partial P)_H$. By the cyclic rule, $\mu_{JT} = -(\partial H/\partial P)_T/(\partial H/\partial T)_P = -(\partial H/\partial P)_T/C_P$.
3. From (3.22) at fixed $T$, $N$: $(\partial H/\partial P)_T = T(\partial S/\partial P)_T + V = -T(\partial V/\partial T)_P + V$ (fourth Maxwell relation).
4. Therefore
$$\mu_{JT} = \frac{V}{C_P}\,(T\alpha - 1). \tag{3.30}$$
∎

For an ideal gas $T\alpha = 1$ and $\mu_{JT} = 0$: no temperature change. For a real gas at low density write the vdW equation as $Pv \approx RT[1 + B(T)/v]$ with $B = b - a/RT$ (expand $1/(v - b) \approx (1 + b/v)/v$), so $v \approx RT/P + B(T)$. Then $T(\partial v/\partial T)_P - v = TB' - B = 2a/RT - b$, and
$$\mu_{JT} \approx \frac{1}{c_P}\left(\frac{2a}{RT} - b\right).$$
Below the **inversion temperature** $T_{\text{inv}} = 2a/(Rb) = \tfrac{27}{4}T_c$, attractions win: the expanding gas does work against its own intermolecular attraction, its kinetic energy drops, and it cools ($\mu_{JT} > 0$ means $T$ falls when $P$ falls). Above it, the repulsive core wins and the gas warms. Linde's liquefier (1895) throttles compressed gas and uses the cooled gas to precool the incoming stream through a counter-current heat exchanger, so the temperature ratchets down until liquid drips out. Hydrogen and helium have inversion temperatures (about 200 K and 40 K) below room temperature and must be precooled before throttling cools them at all.

At finite density the full vdW inversion curve, where $T(\partial v/\partial T)_P = v$, follows from $(\partial v/\partial T)_P = -(\partial P/\partial T)_v/(\partial P/\partial v)_T$. The condition becomes $RT\,v/(v-b)^2 - RT/(v - b) = 2a/v^2$, i.e. $RT = (2a/b)(1 - b/v)^2$, and then $P = RT/(v - b) - a/v^2$. In reduced variables this is $P_r = 24\sqrt{3T_r} - 12T_r - 27$: cooling happens inside the curve, which spans $\tfrac34 < T_r < \tfrac{27}{4}$ and peaks at $P_r = 9$, $T_r = 3$. Project Part 9 plots it and checks the sign of $\mu_{JT}$ on both sides.

### 3.10 The chemical potential

$\mu = (\partial U/\partial N)_{S,V}$ is the energy cost of adding one particle while holding entropy and volume fixed. That definition is correct and nearly useless for intuition, because nobody adds a particle at fixed entropy. Equivalent forms from (3.23), (3.24): $\mu = (\partial F/\partial N)_{T,V} = (\partial G/\partial N)_{T,P}$.

**Particles flow toward lower $\mu$.** Two systems at equal $T$ and fixed volumes exchange particles; $dN_A = -dN_B$. From $dS = -(\mu/T)\,dN$ at fixed $U$, $V$ (and with no net energy flow), $dS_{\text{tot}} = (\mu_B - \mu_A)\,dN_A/T \ge 0$, so $dN_A > 0$ exactly when $\mu_B > \mu_A$. Particles move from high $\mu$ to low $\mu$ until the two are equal, just as energy moves from high $T$ to low $T$.

**Theorem 3.8 (Euler relation, $G = \mu N$, Gibbs–Duhem).** For a simple single-component system,
$$U = TS - PV + \mu N, \tag{3.31}$$
$$G = \mu N, \tag{3.32}$$
$$S\,dT - V\,dP + N\,d\mu = 0. \tag{3.33}$$

▶ **Try it first.** $U(S, V, N)$ is extensive: $U(\lambda S, \lambda V, \lambda N) = \lambda U(S, V, N)$. Apply CM-4's Euler theorem on homogeneous functions, eq. (3.10) there.

**Derivation.**

1. $U$ is homogeneous of degree 1 in $(S, V, N)$. CM-4 (3.10) with $k = 1$: $S\,\partial U/\partial S + V\,\partial U/\partial V + N\,\partial U/\partial N = U$. Substitute (3.15): $TS - PV + \mu N = U$. That is (3.31).
2. Then $G = U - TS + PV = \mu N$. That is (3.32): **for a pure substance, $\mu$ is the Gibbs free energy per particle.**
3. Differentiate (3.31): $dU = T\,dS + S\,dT - P\,dV - V\,dP + \mu\,dN + N\,d\mu$. Subtract (3.14): $0 = S\,dT - V\,dP + N\,d\mu$. That is (3.33). ∎

■ **What just happened.** $\mu$ is not a new independent variable: by Gibbs–Duhem, $T$, $P$, $\mu$ cannot all be varied independently. Dividing (3.33) by $N$: $d\mu = -s\,dT + v\,dP$ with $s, v$ per particle. The Clausius–Clapeyron equation is this line applied twice.

**Chemical potential of an ideal gas.** At fixed $T$, (3.33) gives $d\mu = (V/N)\,dP = k_BT\,dP/P$. Integrating from a reference pressure $P_0$:
$$\mu(T, P) = \mu_0(T) + k_BT\ln\frac{P}{P_0}. \tag{3.34}$$
Thermodynamics leaves $\mu_0(T)$ undetermined; SM-4 and SM-5 compute it ($\mu = k_BT\ln(n\lambda^3)$ with $\lambda$ the thermal wavelength). Squeeze a gas and its $\mu$ rises logarithmically: more "escaping tendency". A warning you should take seriously: $\mu$ will not feel intuitive until SM-5, where it controls Bose–Einstein and Fermi–Dirac occupation. That is normal.

### 3.11 Phase equilibrium

**Coexistence requires equal $\mu$.** A pure substance at fixed $T$, $P$ is split between phases 1 and 2 (say liquid and gas), $N_1 + N_2 = N$ fixed. $G = N_1g_1(T, P) + N_2g_2(T, P)$ with $g_i = \mu_i$ the Gibbs energy per particle of phase $i$ (3.32). Moving $dN$ from phase 2 to phase 1: $dG = (\mu_1 - \mu_2)\,dN$. By (3.26) $G$ decreases to a minimum, so if $\mu_1 < \mu_2$ everything converts to phase 1. Both phases can be present in equilibrium only if
$$\mu_1(T, P) = \mu_2(T, P). \tag{3.35}$$
One equation in two unknowns: coexistence happens on a **curve** $P(T)$ in the $(T, P)$ plane. Three phases coexist only where two such curves cross: a **triple point**.

**Theorem 3.9 (Clausius–Clapeyron).** Along a coexistence curve,
$$\frac{dP}{dT} = \frac{L}{T\,\Delta v}, \tag{3.36}$$
where $L$ is the latent heat per particle (or per mole, with $\Delta v$ per mole) absorbed going from phase 1 to phase 2, and $\Delta v = v_2 - v_1$ the volume change in the same direction.

▶ **Try it first.** Move a small step along the coexistence curve. Both $\mu$s change, but they must stay equal. Use $d\mu = -s\,dT + v\,dP$ for each phase.

**Derivation.**

1. On the curve, $\mu_1 = \mu_2$ at both ends of a small step, so $d\mu_1 = d\mu_2$.
2. Gibbs–Duhem per particle for each phase: $-s_1\,dT + v_1\,dP = -s_2\,dT + v_2\,dP$.
3. Rearrange: $dP/dT = (s_2 - s_1)/(v_2 - v_1)$.
4. The transition happens at constant $T$ reversibly, so $s_2 - s_1 = L/T$ by (3.13). ∎

**Sign conventions.** $L$ and $\Delta v$ must refer to the *same direction* of transition. Liquid → gas: $L > 0$, $\Delta v > 0$, slope positive. Ice → water: $L > 0$ but water is denser, so $\Delta v < 0$ and the melting curve slopes *backward*: pressure lowers the melting point (Exercise B7 gives about −135 bar per kelvin).

**Integrated form for vaporization.** Neglect $v_{\text{liquid}}$ next to $v_{\text{gas}}$ and treat the vapour as ideal, $v_g = RT/P$ per mole. Then $dP/dT = LP/(RT^2)$, i.e. $d\ln P/dT = L/(RT^2)$, and with $L$ constant
$$\ln\frac{P_2}{P_1} = -\frac{L}{R}\left(\frac1{T_2} - \frac1{T_1}\right). \tag{3.37}$$
Example 4.6 uses it for boiling at altitude. Project Part 8 shows that constant $L$ overestimates water's vapour pressure by 18% at 25 °C, and that letting $L$ vary as $L_0 + \Delta C_P(T - T_0)$ (Kirchhoff's relation, from $dL/dT \approx C_{P,\text{gas}} - C_{P,\text{liq}}$) brings the error down to 3%.

**The phase diagram of water, in words.** Three regions (ice, liquid, vapour) meeting at the triple point, 273.16 K and 611.657 Pa. The vaporization curve rises from there and ends at the critical point (647.1 K, 22.06 MPa), beyond which liquid and gas are no longer distinct: you can go from one to the other around the critical point without crossing any boundary. The sublimation curve runs below the triple point. The melting curve leaves the triple point tilted slightly to the left (negative slope, from (3.36)) and runs up to about 200 MPa, where other crystal structures of ice (ice III, V, VI, VII, and more) take over; water has more than a dozen solid phases.

**The van der Waals liquid–gas transition and Maxwell's construction.** Below $T_c$ a vdW isotherm $P(v)$ has a loop: a stretch where $(\partial P/\partial v)_T > 0$, i.e. $\kappa_T < 0$, which is mechanically unstable (squeeze it and the pressure drops, so it collapses further). The critical point, where the loop first appears, is where $\partial P/\partial v$ and $\partial^2P/\partial v^2$ both vanish: $RT/(v - b)^2 = 2a/v^3$ and $RT/(v - b)^3 = 3a/v^4$. Dividing gives $v - b = 2v/3$, so $v_c = 3b$, then $T_c = 8a/(27Rb)$ and $P_c = a/(27b^2)$. For CO₂ ($a = 0.364$ Pa m⁶ mol⁻², $b = 4.267\times10^{-5}$ m³/mol): $T_c = 304.0$ K, $P_c = 7.40$ MPa, close to the measured 304.1 K and 7.38 MPa because $a$ and $b$ are fitted to the critical point.

The real system does not follow the loop. At some pressure $P^*$ it splits into liquid at $v_l$ and gas at $v_g$, both on the isotherm at $P^*$, and the isotherm is flat between them. Which $P^*$?

▶ **Try it first.** Coexistence needs $\mu_l = \mu_g$ at the same $T$, $P^*$. At fixed $T$, $d\mu = v\,dP$ along the isotherm. Integrate from the liquid end to the gas end along the vdW curve.

**Derivation.**

1. The vdW free energy $f(T, v) = -RT\ln(v - b) - a/v + f_0(T)$ (from integrating $P = -(\partial f/\partial v)_T$) defines a smooth $\mu = f + Pv$ at every point of the isotherm, including the unstable part. At fixed $T$, $d\mu = v\,dP$ (Gibbs–Duhem with $dT = 0$).
2. Equal chemical potentials at the two ends: $0 = \mu_g - \mu_l = \int_{l}^{g} v\,dP$, integrated along the isotherm.
3. Integrate by parts: $\int_l^g v\,dP = [Pv]_l^g - \int_{v_l}^{v_g}P\,dv = P^*(v_g - v_l) - \int_{v_l}^{v_g}P\,dv$.
4. So
$$\int_{v_l}^{v_g}\big[P(v) - P^*\big]\,dv = 0: \tag{3.38}$$
the horizontal line at $P^*$ cuts the loop into two lobes of equal area. ∎

■ **What just happened.** The equal-area rule is equal $\mu$ in geometric form. Project Part 7 finds $P^*$ by root-finding on (3.38) for CO₂ at temperatures from 230 K to 303 K, checks $\mu_l = \mu_g$ directly, and plots the coexistence dome. At 280 K it gives $P^* = 5.28$ MPa against a measured vapour pressure of about 4.2 MPa: the vdW equation is qualitatively right and quantitatively mediocre away from the critical point. SM-6 revisits this as a mean-field theory and explains why it gets the behaviour near $T_c$ wrong too.

### 3.12 What thermodynamics cannot tell you

Everything above follows from two laws plus an equation of state and a heat capacity. Thermodynamics cannot supply those inputs, and it cannot justify itself:

- **The equation of state** ($PV = Nk_BT$, van der Waals, water's steam tables). Measured, or derived from a microscopic model in SM-4.
- **Heat capacities.** Why $C_V = \tfrac32Nk_B$ for helium, why it freezes out for molecular vibrations (SM-1's "first crack" in equipartition), why it vanishes as $T \to 0$. SM-4 and SM-5.
- **The additive constant in $S$**, $\sigma$ in (3.16), and $\mu_0(T)$ in (3.34). Classical thermodynamics is blind to it; it contains Planck's constant. SM-3.
- **Why the second law holds at all.** Mechanics is time-reversible; nothing in Newton's laws prefers heat flowing hot to cold. The answer is statistical: the second law is overwhelmingly probable rather than certain. SM-3.
- **Fluctuations.** Thermodynamics speaks of exact values of $U$, $P$; real systems fluctuate around them by $\sim 1/\sqrt N$. The project's particle engine fluctuates visibly at $N = 40\,000$. SM-4.

That list is the motivation for statistical mechanics.

## 4. Worked examples

### Example 4.1: Carnot engine between 500 K and 300 K, with entropy bookkeeping

One mole of ideal gas runs a Carnot cycle between $T_h = 500$ K and $T_c = 300$ K; on the hot isotherm the volume doubles.

- Efficiency (3.3): $\eta = 1 - 300/500 = 0.400$.
- Hot isotherm: $Q_h = RT_h\ln 2 = 8.314 \times 500 \times 0.6931 = 2881.6$ J.
- Step 4 of §3.1: the cold isotherm compresses by the same factor 2, so $Q_c = RT_c\ln2 = 1728.9$ J.
- Work: $W = Q_h - Q_c = 1152.6$ J, and $W/Q_h = 0.400$. ✓
- Entropy: on the hot leg the gas gains $Q_h/T_h = R\ln2 = 5.763$ J/K and the hot reservoir loses the same. On the adiabats nothing changes. On the cold leg the gas loses $Q_c/T_c = 5.763$ J/K and the cold reservoir gains it. Per cycle: $\Delta S_{\text{gas}} = 0$ (state function, closed cycle), $\Delta S_{\text{hot}} + \Delta S_{\text{cold}} = -5.763 + 5.763 = 0$. The universe's entropy is unchanged: that is what reversible means.

Entropy flows *through* the engine from the hot reservoir to the cold one, unchanged in amount; energy flows through too, but only $Q_c$ of it has to go along with the entropy. The rest, $W$, leaves as work, which carries no entropy. That is the picture to keep.

### Example 4.2: Free expansion vs reversible isothermal expansion, 1 mol, factor 2

One mole of ideal gas at 300 K doubles its volume.

*Reversibly, isothermally* (in contact with a 300 K reservoir, piston withdrawn slowly): $\Delta S_{\text{gas}} = R\ln2 = 5.763$ J/K by (3.10). The gas does work $W = RT\ln 2 = 1728.9$ J and absorbs $Q = 1728.9$ J from the reservoir. $\Delta S_{\text{res}} = -Q/T = -5.763$ J/K. $\Delta S_{\text{universe}} = 0$.

*Freely* (insulated box, partition removed): $Q = 0$, $W = 0$, $\Delta U = 0$, final $T = 300$ K, final volume doubled. Same end state, so $\Delta S_{\text{gas}} = 5.763$ J/K. No reservoir is involved, so $\Delta S_{\text{universe}} = +5.763$ J/K. The entropy increase is the same; what differs is that in the reversible case it was paid for by the reservoir, and in the free expansion it was created from nothing. The 1728.9 J of work you could have extracted is gone; $T\,\Delta S_{\text{universe}}$ is exactly the lost work.

### Example 4.3: Two bricks at 400 K and 200 K

Two identical 2 kg bricks ($c = 840$ J kg⁻¹ K⁻¹, so $C = 1680$ J/K each) at 400 K and 200 K are pressed together and insulated.

- Energy: $C(T_f - 400) + C(T_f - 200) = 0$, so $T_f = 300$ K.
- Hot brick: $\Delta S = C\ln(300/400) = 1680\times(-0.2877) = -483.3$ J/K.
- Cold brick: $\Delta S = C\ln(300/200) = 1680\times0.4055 = +681.2$ J/K.
- Total: $\Delta S_{\text{tot}} = C\ln(9/8) = 1680 \times 0.1178 = +197.9$ J/K.

The cold brick gains more entropy than the hot one loses, because the same heat (168 kJ) enters at lower temperature. If instead you ran a reversible engine between the bricks until they reached a common temperature, $\Delta S_{\text{tot}} = 0$ forces $T_f = \sqrt{400\times200} = 282.8$ K, and the energy difference becomes work: Exercise C1 finds $57.6$ kJ.

### Example 4.4: The Legendre transform from $U(S, V)$ to $F(T, V)$ for an ideal gas

Take one mole of monatomic ideal gas ($f = 3$). Fix a reference state $(T_0, V_0)$ with $U_0 = \tfrac32RT_0$, and write (3.16) relative to it:
$$S(U, V) = S_0 + \tfrac32R\ln\frac{U}{U_0} + R\ln\frac{V}{V_0}.$$

**Invert to get $U(S, V)$:**
$$U(S, V) = U_0\left(\frac{V_0}{V}\right)^{2/3}\exp\!\left[\frac{2(S - S_0)}{3R}\right].$$
Check (3.15): $T = \partial U/\partial S = 2U/(3R)$, so $U = \tfrac32RT$ ✓; $P = -\partial U/\partial V = \tfrac{2}{3}U/V = RT/V$ ✓. And $\partial^2U/\partial S^2 = 4U/(9R^2) > 0$: convex, so the transform is legitimate.

**Transform.** Minimize $U(S) - TS$ over $S$: the minimizer has $\partial U/\partial S = T$, i.e. $U = \tfrac32RT$, and at that point $S = S_0 + \tfrac32R\ln(T/T_0) + R\ln(V/V_0)$. So
$$F(T, V) = \tfrac32RT - T\left[S_0 + \tfrac32R\ln\frac{T}{T_0} + R\ln\frac{V}{V_0}\right].$$

**Check the derivatives (3.23):** $-\partial F/\partial V = RT/V = P$ ✓. $-\partial F/\partial T = -\tfrac32R + S_0 + \tfrac32R\ln(T/T_0) + R\ln(V/V_0) + \tfrac32R = S$ ✓.

**Numbers.** With $T_0 = 300$ K, $V = V_0$, and $S_0 = 0$ (the constant thermodynamics cannot fix; it shifts $F$ by $-TS_0$): $F(300\text{ K}) = \tfrac32R(300) = 3741.5$ J; $F(600\text{ K}) = \tfrac32R(600) - 600\cdot\tfrac32R\ln2 = 7483.0 - 5186.8 = 2296.2$ J; $F(150\text{ K}) = 3167.5$ J. Project Part 5 tabulates $U(S)$ on a grid, computes $\min_S[U(S) - TS]$ numerically, and matches all three to better than $10^{-3}$, then recovers $S(600\text{ K}) = -\partial F/\partial T = \tfrac32R\ln2 = 8.645$ J/K by finite differences. This is CM-8's $\tfrac12mv^2 \to p^2/2m$ example with a different convex function.

### Example 4.5: $C_P - C_V$ for liquid water at 25 °C

Data at 25 °C, 1 atm: density 997.05 kg/m³, so molar volume $v = 0.018015/997.05 = 1.8068\times10^{-5}$ m³/mol; $\alpha = 2.572\times10^{-4}$ K⁻¹; $\kappa_T = 4.524\times10^{-10}$ Pa⁻¹.

By (3.29):
$$c_P - c_v = \frac{Tv\alpha^2}{\kappa_T} = \frac{298.15 \times 1.8068\times10^{-5}\times(2.572\times10^{-4})^2}{4.524\times10^{-10}} = 0.788\ \mathrm{J\,mol^{-1}K^{-1}},$$
or $43.7$ J kg⁻¹ K⁻¹. With $c_P = 4181.3$ J kg⁻¹ K⁻¹ that gives $c_v \approx 4137.6$ J kg⁻¹ K⁻¹ and $\gamma = c_P/c_v = 1.0106$. Tabulated values (IAPWS-95 formulation) at 25 °C give $c_P = 4.1813$ and $c_v \approx 4.138$ kJ kg⁻¹ K⁻¹, a difference of about 44 J kg⁻¹ K⁻¹, matching. For a liquid $C_P$ and $C_V$ differ by about 1%; for a diatomic gas by 40%. That is why "the heat capacity of water" is one number for most purposes, and why it is not quite.

### Example 4.6: Boiling point of water at 3000 m

The standard atmosphere gives $P = 101325\,(1 - 2.25577\times10^{-5}\,h)^{5.25588}$ Pa; at $h = 3000$ m, $P = 70.11$ kPa. Take $L = 40.66$ kJ/mol at the normal boiling point $T_1 = 373.15$ K, $P_1 = 101.325$ kPa. Solve (3.37) for $T_2$:
$$\frac{1}{T_2} = \frac{1}{T_1} - \frac{R}{L}\ln\frac{P_2}{P_1} = \frac{1}{373.15} - \frac{8.314}{40660}\ln\frac{70.11}{101.325} = 2.6799\times10^{-3} + 7.53\times10^{-5},$$
$T_2 = 362.95$ K $= 89.8$ °C. Letting $L$ vary with $T$ (project Part 8) gives 89.9 °C; the correction is tiny over a 10 K range. Your pasta takes longer at altitude because the water boils about 10 °C cooler, not because anything is wrong with the stove. A pressure cooker works the same way in reverse (Exercise B7).

## 5. Why it matters and where it connects

- **SM-3** derives $S = k_B\ln\Omega$, recovers (3.16) with the constant $\sigma$ (Sackur–Tetrode), resolves the Gibbs paradox of §3.5, and explains why (3.8) holds.
- **SM-4** computes $F = -k_BT\ln Z$ directly from a partition function; every relation in §3.8–3.9 then applies to it. The canonical ensemble is the microscopic version of "system in a bath at $T$," which is exactly what the Legendre transform to $F$ describes.
- **SM-5** needs $\mu$ from §3.10 to state Bose–Einstein and Fermi–Dirac statistics; $G = \mu N$ and (3.34) are where you first meet it.
- **SM-6** revisits the van der Waals loop and the equal-area rule of §3.11 as mean-field theory, and the convexity of Theorem 3.6 becomes the reason free energies have flat pieces at first-order transitions.
- **CM-8.** Definition 3.1 and Theorem 3.1 there are the same Legendre transform as §3.8 here, with $L(\dot q) \to H(p)$ in place of $U(S) \to F(T)$. The "derivative of the new function is the old variable" identity, CM-8 (3.3), is $(\partial F/\partial T)_V = -S$. The Maxwell relations are the thermodynamic analogue of the symmetry of mixed partials that makes Hamilton's equations consistent.
- **CM-4.** Euler's theorem on homogeneous functions, (3.10) there, is what gives $G = \mu N$ here.
- **SM-7 and ML.** In SM-7, the variational free energy of machine learning (the negative ELBO) is $\langle E\rangle - TS$ with $T = 1$: the same $F = U - TS$, minimized for the same reason.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are mechanical: COPs, entropy changes for heating, gas processes and melting, the differentials of the potentials, an Otto efficiency. B1–B9 are the module: B1 (a second proof that reversible engines tie), B2 (entropy generated by an irreversible compression), B3 (stability and concavity), B4 (vdW free expansion), B5 ($C_P - C_V$ for a vdW gas), B6 (Joule–Thomson inversion), B7 (Clausius–Clapeyron: pressure cooker and ice), B8 (Gibbs–Duhem and $\mu$ of the ideal gas), and B9 (`project.py` extension: latent heat of vdW CO₂ two ways). Essential: B2, B3, B5, B7, B9. C1–C3: maximum work from two finite bodies, thermodynamic stability conditions, and a heat pump fed by a power plant.

## 7. Coding project

`project.py` builds and checks:

1. **A Carnot engine made of particles.** 40,000 non-interacting particles on a line between a wall and a piston. During isothermal legs the wall is a thermal wall that re-emits particles with speeds drawn from the flux distribution at the bath temperature; during adiabatic legs it reflects. The piston moves slowly. Heat and work are counted particle by particle, the first law is asserted to machine precision, and the measured efficiency is compared with $1 - T_c/T_h$: 0.485 against 0.500 and 0.284 against 0.300, both below Carnot as they must be. The temperature after an adiabat is checked against $TL^2 = $ const ($\gamma = 3$ in one dimension). A run with the piston ten times faster gives 0.405: finite speed is irreversibility, and irreversibility costs work.
2. **Worked-example numbers** for Examples 4.1–4.3, melting, mixing, Otto, Diesel, heat pump.
3. **Entropy production.** Two bodies with $S_i = C_i\ln U_i$; $S_1 + S_2$ is maximal at equal temperatures, and 1000 random initial splits all gain entropy when heat flows hot to cold.
4. **SymPy.** The Maxwell relation $(\partial S/\partial V)_T = (\partial P/\partial T)_V$ for ideal and vdW gases from $F$; $C_V$ independent of $V$; $C_P - C_V$ from (3.29) for both (ideal gives $nR$); a Maxwell relation from $U(S,V)$; $\mu_{JT} = 0$ for the ideal gas.
5. **Numerical Legendre transform** of a tabulated $U(S)$ against the closed-form $F(T)$ (Example 4.4).
6. **Example 4.5.**
7. **Maxwell's equal-area construction** for vdW CO₂, with coexistence pressures by root-finding, a check that $\mu_l = \mu_g$, and the coexistence dome.
8. **Clausius–Clapeyron for water**, constant $L$ vs $L(T)$ integrated with `solve_ivp`, against IAPWS saturation pressures, and Example 4.6.
9. **The vdW Joule–Thomson inversion curve**, checked against $P_r = 24\sqrt{3T_r} - 12T_r - 27$.
10. **Numbers for the exercise solutions**, then "Your turn" and the reference solution to B9.

**Build it yourself:** before reading Part 1, write `run_leg` yourself: positions and velocities as arrays, a time step, a reflecting or thermal wall at 0, a moving piston at $L(t)$, and two accumulators (heat at the wall, work at the piston). The thermal wall is the subtle part: a wall at temperature $T$ does not emit Maxwellian velocities but flux-weighted ones, $p(v) \propto v\,e^{-v^2/2T}$, because fast particles hit it more often. Sample with $v = \sqrt{-2T\ln u}$. Then write `coexistence(T)` for Part 7 before reading it.

**Your turn** (in the script): a Stirling cycle with and without regeneration; the Gibbs paradox numerically from (3.16); heat-pump COP against outdoor temperature and the break-even point with a gas furnace.

## 8. Where you'll struggle

- **Treating the dozen statements of the second law as different laws.** Kelvin, Clausius, Carnot's theorem, the Clausius inequality, $\Delta S_{\text{universe}} \ge 0$, "$G$ decreases at fixed $T, P$": one law. The fix is to be able to derive each from one other in a few lines, which §3.2–3.4 and (3.26) do. If you can go Kelvin → Carnot → Clausius inequality → $S$ exists → $\Delta S \ge 0$ → $\Delta G \le 0$ closed-book, you understand the second law. If you are memorizing them separately, you don't yet.
- **"Entropy is disorder."** It is a slogan that will mislead you. A supercooled liquid freezing spontaneously looks more ordered and the total entropy still rises (the latent heat raises the surroundings' entropy). Entropy is a count of microstates (SM-3), or here, $\int đQ_{\text{rev}}/T$. Use those.
- **Free expansion.** The test of whether you understand that $S$ is a state function. If you find yourself saying "$Q = 0$ so $\Delta S = 0$," stop: $dS = đQ/T$ holds only for reversible heat. Find a reversible path between the same end states (§3.5).
- **Which reservoir's $T$ goes in the Clausius inequality.** The reservoir's. During an irreversible step the system may not have a temperature.
- **Legendre direction.** Which variable is traded for which, and the sign. In $F = U - TS$ the new variable is $T = \partial U/\partial S$; the old one comes back as $S = -\partial F/\partial T$. It is the same confusion as $L \to H$ in CM-8. The fix: always write the differential ($dF = -S\,dT - P\,dV$) before using a potential; it tells you the natural variables and the signs at once.
- **Maxwell relations.** They are trivial to derive (one line each from the right potential) and students memorize them instead, then get signs wrong under exam pressure. Never write one without the potential it came from. And never drop the subscripts: $(\partial S/\partial V)_T$ and $(\partial S/\partial V)_P$ are different quantities.
- **The chemical potential** will confuse you until SM-5. That's normal. For now hold on to three facts: particles flow from high $\mu$ to low $\mu$; $\mu = G/N$ for a pure substance; coexisting phases have equal $\mu$.
- **Clausius–Clapeyron signs.** $L$ and $\Delta v$ must describe the same direction of the transition. Ice → water has $L > 0$ and $\Delta v < 0$. Mixing directions gives you a melting curve with the wrong slope.
- **The particle engine's efficiency is noisy.** With $N = 40\,000$, energy fluctuates by about 1%, so the measured efficiency wanders by about 0.01 run to run. That is SM-4's $1/\sqrt N$ showing up, not a bug.

## 9. Check yourself

1. State Kelvin's and Clausius's forms of the second law. Prove that violating Kelvin lets you violate Clausius.
2. Why can the efficiency of a reversible engine depend only on the two reservoir temperatures?
3. In $\oint đQ/T \le 0$, whose temperature is $T$?
4. A gas undergoes free expansion. $Q = 0$. Why is $\Delta S \ne 0$, and how do you compute it?
5. Show in one line that heat flowing from hot to cold increases total entropy.
6. Write the COP bound for a heat pump. Why can a COP exceed 1?
7. Derive $\eta_{\text{Otto}} = 1 - r^{1-\gamma}$ from the two adiabats.
8. Show $\partial^2S/\partial U^2 = -1/(T^2C_V)$. Why does stability require $C_V > 0$?
9. Write $dF$ and $dG$. From $dG$, derive a Maxwell relation.
10. Why is $G$ minimized at fixed $T$ and $P$? (Three lines, from system plus reservoir.)
11. Derive $C_P - C_V = TV\alpha^2/\kappa_T$. When is it zero?
12. Why does throttling conserve $H$? Why does an ideal gas not cool on throttling?
13. Derive $G = \mu N$ from extensivity. Then derive Clausius–Clapeyron from Gibbs–Duhem.

## 10. Pair with

- **Schroeder, *An Introduction to Thermal Physics*, ch. 4–5.** Engines, refrigerators, free energy, and phase transformations with excellent physical examples and numbers. Schroeder introduces entropy statistically first (his ch. 2–3, our SM-3), so read his ch. 4–5 for the applications rather than the logic.
- **Callen, *Thermodynamics and an Introduction to Thermostatistics*, ch. 1–7.** The axiomatic version: entropy postulated as a concave, extensive function, everything else derived as geometry. Read it for the structure of §3.8–3.10; nothing is cleaner on Legendre transforms, Euler and Gibbs–Duhem, and stability.
- **Fermi, *Thermodynamics*, ch. 3–5.** The engine-first route of §3.2–3.4 done briefly and rigorously by a master; his proof of the Clausius inequality is the one used here.
- **Feynman Lectures on Physics, Vol. I, ch. 44–46.** Ch. 44 runs Carnot's argument with Feynman's weight-lifting machines; ch. 45 derives the Clausius–Clapeyron equation and other illustrations; ch. 46 (the ratchet and pawl) is the best short explanation of why you cannot beat the second law with a clever gadget.
