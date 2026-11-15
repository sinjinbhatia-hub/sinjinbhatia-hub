# SM-2 Solutions

Every number below is printed by Part 10 of `project.py` (lines labelled with the exercise number) or by the B9 reference solution at the end of the script. Equation numbers refer to `notes.md`.

## A. Warm-up

**A1.** (a) By (3.17), $\mathrm{COP}_{\text{fridge}} \le T_c/(T_h - T_c) = 275/20 = 13.75$.

(b) $\mathrm{COP} = Q_c/W$, so $W = Q_c/\mathrm{COP} = 300/3 = 100$ W of electrical power. The fridge is cyclic, so by the first law the heat into the kitchen is $Q_h = Q_c + W = 300 + 100 = 400$ W. A running fridge with its door open heats the kitchen: it moves 300 W from inside to outside the box and adds its own 100 W.

**A2.** The end states are equilibrium states at 1 atm, so use a reversible path at constant pressure: $đQ_{\text{rev}} = mc_P\,dT$ and, by (3.9) with $C_P$ in place of $C_V$,
$$\Delta S = mc_P\ln\frac{T_2}{T_1} = 4186\ln\frac{353.15}{293.15} = 779.5\ \mathrm{J/K}.$$
It does not depend on how you heat it: $S$ is a state function and the water's end states are the same. What does depend on the method is the entropy change of the *heat source*. A hotplate barely warmer than the water at each moment loses almost exactly 779.5 J/K; a red-hot rock at ~1000 K loses much less (it gives up the same heat at higher temperature), so $\Delta S_{\text{universe}}$ is larger.

**A3.** Choose a reversible path: heat at constant volume from $T_1$ to $2T_1$, then expand isothermally from $V_1$ to $3V_1$. Monatomic means $c_v = \tfrac32R$ per mole. By (3.9) and (3.10),
$$\Delta S = n c_v\ln 2 + nR\ln 3 = 2R\left(\tfrac32\ln2 + \ln3\right) = 2(8.314)(1.0397 + 1.0986) = 35.56\ \mathrm{J/K}.$$
The actual process does not matter because $S$ is a function of state: (3.16) gives $S$ as a function of $(T, V)$ at fixed $N$, and any path between the same end states gives the same $\Delta S$ of the gas. (The universe's $\Delta S$ does depend on the process.)

**A4.**
1. $H = U + PV$: $dH = dU + P\,dV + V\,dP = (T\,dS - P\,dV + \mu\,dN) + P\,dV + V\,dP = T\,dS + V\,dP + \mu\,dN$.
2. $F = U - TS$: $dF = dU - T\,dS - S\,dT = -S\,dT - P\,dV + \mu\,dN$.
3. $G = U - TS + PV$: $dG = dU - T\,dS - S\,dT + P\,dV + V\,dP = -S\,dT + V\,dP + \mu\,dN$.

From step 1, $H$ is naturally a function of $(S, P, N)$ with $(\partial H/\partial S)_{P,N} = T$ and $(\partial H/\partial P)_{S,N} = V$. Equality of mixed partials, $\partial^2H/\partial P\,\partial S = \partial^2 H/\partial S\,\partial P$, gives
$$\left(\frac{\partial T}{\partial P}\right)_{S} = \left(\frac{\partial V}{\partial S}\right)_{P},$$
the second relation in (3.27).

**A5.** By (3.19), $\eta = 1 - r^{1-\gamma} = 1 - 8^{-0.4} = 0.565$. On the adiabatic compression $TV^{\gamma-1}$ is constant, so $T_2 = T_1r^{\gamma-1} = 300\times8^{0.4} = 689.2$ K. (Real engines are lower in efficiency for the reasons in §3.7.)

**A6.** The ice melts at constant $T_m = 273.15$ K, absorbing $Q = mL = 0.05\times334\times10^3 = 16\,700$ J reversibly from its own point of view, so by (3.13) $\Delta S_{\text{ice}} = 16\,700/273.15 = +61.14$ J/K. The room is a reservoir at 293.15 K that gives up the same heat: $\Delta S_{\text{room}} = -16\,700/293.15 = -56.97$ J/K. Total: $\Delta S_{\text{universe}} = 16\,700\,(1/273.15 - 1/293.15) = +4.17$ J/K. The positive total is what says the process runs spontaneously, by (3.8); the reverse (water at 0 °C freezing in a 20 °C room) would have $-4.17$ J/K and does not happen. Note the ice became *less* ordered and the room's entropy fell; only the sum matters.

## B. Core

**B1.** Let $X$ have efficiency $\eta_X$ and $R$ be reversible with $\eta_R$, and suppose $\eta_X > \eta_R$.

1. Run $X$ forward: it takes $Q_h$ from the hot reservoir, delivers work $W = \eta_XQ_h$, and rejects $Q_h - W$ to the cold reservoir.
2. Run $R$ backward as a refrigerator driven by exactly that work $W$. Forward, $R$ would take $Q_h'$ from the hot reservoir and produce $W = \eta_RQ_h'$; reversed, it consumes $W$, takes $Q_h' - W$ from the cold reservoir, and delivers $Q_h' = W/\eta_R$ to the hot one.
3. Because $\eta_R < \eta_X$, $Q_h' = W/\eta_R > W/\eta_X = Q_h$.
4. Net effect of the pair per cycle: work produced minus work consumed is $W - W = 0$. The hot reservoir gains $Q_h' - Q_h > 0$. The cold reservoir loses $(Q_h' - W) - (Q_h - W) = Q_h' - Q_h$.
5. So heat $Q_h' - Q_h > 0$ has moved from the cold reservoir to the hot one with no other change anywhere (both engines are back in their initial states, no net work). That violates Clausius. Hence $\eta_X \le \eta_R$. ∎

This is part (a) of Theorem 3.2; part (b) follows exactly as in the notes, by exchanging the roles of two reversible engines.

**B2.** (a) The gas ends at 300 K, so $\Delta U = 0$ (ideal gas). The external pressure is constant at $2P_1$ and the volume falls by $V_1/2$, so the work done *on* the gas is
$$W_{\text{on}} = P_{\text{ext}}\,|\Delta V| = 2P_1\cdot\tfrac12V_1 = P_1V_1 = RT = 8.314\times300 = 2494.3\ \mathrm{J}.$$
By the first law, $Q = \Delta U - W_{\text{on}} = -2494.3$ J: the gas gives 2494.3 J to the reservoir.

Gas: same end states as a reversible isothermal compression to half the volume, so by (3.10) $\Delta S_{\text{gas}} = R\ln\tfrac12 = -5.763$ J/K.

Reservoir: it receives 2494.3 J at 300 K, $\Delta S_{\text{res}} = 2494.3/300 = +8.314$ J/K ($= R$ exactly, since $Q = RT$).

Universe: $\Delta S = R(1 - \ln 2) = 8.314 - 5.763 = +2.551$ J/K $> 0$, as (3.8) requires.

(b) A reversible isothermal compression needs $W_{\text{rev}} = RT\ln2 = 1728.9$ J. The excess is $2494.3 - 1728.9 = 765.4$ J, and $T\,\Delta S_{\text{universe}} = 300\times2.551 = 765.4$ J. They are equal: $W_{\text{on}} - W_{\text{rev}} = RT - RT\ln 2 = T\cdot R(1 - \ln 2)$. The extra work you spent was degraded into extra heat dumped into the reservoir, and that extra heat divided by $T$ is precisely the entropy generated.

**B3.** (a) $1/T = dS/dU = 2k_BU/\varepsilon^2$, so $T = \varepsilon^2/(2k_BU)$. Inverting, $U = \varepsilon^2/(2k_BT)$ and
$$C_V = \frac{dU}{dT} = -\frac{\varepsilon^2}{2k_BT^2} < 0.$$
Adding energy *lowers* the temperature. This matches (3.21): $d^2S/dU^2 = 2k_B/\varepsilon^2 > 0$, so $S$ is convex, and $-1/(T^2C_V) > 0$ forces $C_V < 0$.

(b) Let each half have entropy $S_h(x)$ with the same convex form, $S_h'' > 0$. Move energy $\delta$ from one half to the other:
$$S_{\text{tot}}(\delta) = S_h(U/2 + \delta) + S_h(U/2 - \delta) = 2S_h(U/2) + S_h''(U/2)\,\delta^2 + O(\delta^4).$$
The $\delta^2$ coefficient is positive, so $\delta = 0$ is a minimum of $S_{\text{tot}}$: any fluctuation raises the entropy and is amplified. Physically, the half that gains a little energy gets *colder* (because $C_V < 0$), so heat flows *into* it from the now-hotter other half, and the imbalance runs away until one half holds essentially all the energy. The uniform state is not an equilibrium.

(c) $U(S) = \varepsilon\sqrt{S/k_B}$, and $U''(S) = -\tfrac14\varepsilon k_B^{-1/2}S^{-3/2} < 0$: concave. Then $U(S) - TS$ has a stationary point where $U' = T$, but since $U - TS$ is concave in $S$ that point is a *maximum*. The minimum over $S$ runs off to $S \to \infty$, where $U - TS \to -\infty$. So $\min_S[U - TS] = -\infty$ for every $T$: the transform carries no information, and the map $S \leftrightarrow T$ is not one-to-one in the useful direction. The Legendre transform of §3.8 needs convexity, and convexity is $C_V > 0$ (Theorem 3.6). A system in a heat bath with $C_V < 0$ has no stable equilibrium at the bath temperature.

**B4.** Free expansion: $Q = 0$, $W = 0$, so $u$ is conserved.

1. $c_vT_i - a/v_1 = c_vT_f - a/v_2$, so $T_f = T_i - \dfrac{a}{c_v}\left(\dfrac1{v_1} - \dfrac1{v_2}\right)$.
2. $c_v = \tfrac52R = 20.786$ J mol⁻¹ K⁻¹; $1/v_1 - 1/v_2 = 1000 - 500 = 500$ m⁻³ mol. So $T_f = 300 - (0.364/20.786)(500) = 300 - 8.756 = 291.24$ K.
3. Entropy from (3.28), a state function, between $(300\ \text{K}, 1\ \text{L})$ and $(291.24\ \text{K}, 2\ \text{L})$:
$$\Delta s = c_v\ln\frac{291.24}{300} + R\ln\frac{2\times10^{-3} - b}{1\times10^{-3} - b} = -0.616 + 8.314\ln\frac{1.95733}{0.95733} = -0.616 + 5.946 = 5.331\ \mathrm{J/K}.$$
4. Ideal gas: $R\ln 2 = 5.763$ J/K, with no temperature change.

The vdW gas cools because $u = c_vT - a/v$: as $v$ grows the (negative) attraction energy $-a/v$ rises toward zero, and with $u$ fixed that energy has to come from the kinetic part $c_vT$. Molecules climbing out of each other's attractive wells slow down. $\Delta s$ is still positive (the process is irreversible and isolated, (3.7)), but smaller than ideal: the cooling term is negative, and excluded volume makes the effective volume ratio $(v_2 - b)/(v_1 - b) = 2.045$ slightly more than 2, which partly offsets it.

**B5.** (a) The second form of (3.29), per mole, is $c_P - c_v = -T(\partial P/\partial T)_v^2/(\partial P/\partial v)_T$. For $P = RT/(v - b) - a/v^2$:
$$\left(\frac{\partial P}{\partial T}\right)_v = \frac{R}{v - b},\qquad \left(\frac{\partial P}{\partial v}\right)_T = -\frac{RT}{(v - b)^2} + \frac{2a}{v^3}.$$
Then
$$c_P - c_v = \frac{TR^2/(v - b)^2}{RT/(v - b)^2 - 2a/v^3}.$$
Divide numerator and denominator by $RT/(v - b)^2$:
$$c_P - c_v = \frac{R}{1 - \dfrac{2a(v - b)^2}{RTv^3}}.$$
Project Part 4 derives the same expression with SymPy and Part 10 checks the numerical value against it.

(b) $a \to 0$: the denominator is 1, so $c_P - c_v = R$, SM-1's ideal-gas result. (With $a = 0$ but $b \ne 0$ it is also exactly $R$.)

(c) $v = 10^{-3}$ m³/mol, $T = 300$ K: $2a(v - b)^2/(RTv^3) = 2(0.364)(9.5733\times10^{-4})^2/(8.314\times300\times10^{-9}) = 0.2675$, so $c_P - c_v = 8.314/0.7325 = 11.35$ J mol⁻¹ K⁻¹, 37% above $R$.

(d) At the critical point $(\partial P/\partial v)_T = 0$, so the denominator in (a) goes to zero and $c_P - c_v \to \infty$. In terms of (3.29): $\kappa_T \to \infty$ (the fluid is infinitely compressible there), so at constant pressure a tiny heat input produces a huge expansion, and almost all the heat goes into it. $c_P$ diverges at the critical point; that is a real effect (it is measured), though the exponent vdW predicts is wrong (SM-6).

**B6.** (a) Per mole, (3.30) reads $\mu_{JT} = [T(\partial v/\partial T)_P - v]/c_P$. At low density expand the vdW equation in $1/v$:
1. $\dfrac{RT}{v - b} = \dfrac{RT}{v}\left(1 - \dfrac bv\right)^{-1} \approx \dfrac{RT}{v}\left(1 + \dfrac bv\right)$, so $P \approx \dfrac{RT}{v}\left[1 + \dfrac{B(T)}{v}\right]$ with $B(T) = b - \dfrac{a}{RT}$.
2. Solve for $v$ to first order in $B$: $v \approx \dfrac{RT}{P} + B(T)$ (replace $B/v$ by its leading value $BP/RT$).
3. $T(\partial v/\partial T)_P = \dfrac{RT}{P} + TB'(T)$, with $B' = a/(RT^2)$. So $T(\partial v/\partial T)_P - v = TB' - B = \dfrac{a}{RT} - b + \dfrac{a}{RT} = \dfrac{2a}{RT} - b$.
4. Therefore $\mu_{JT} \approx \dfrac{1}{c_P}\left(\dfrac{2a}{RT} - b\right)$, which changes sign at $T_{\text{inv}} = 2a/(Rb)$.

(b) N₂: $T_{\text{inv}} = 2(0.1370)/(8.314\times3.87\times10^{-5}) = 852$ K. H₂: $2(0.02476)/(8.314\times2.661\times10^{-5}) = 224$ K.

(c) With $c_P = \tfrac72R = 29.10$ J mol⁻¹ K⁻¹ and $1\ \text{bar} = 10^5$ Pa:
N₂ at 300 K: $2a/RT = 1.098\times10^{-4}$, minus $b$ gives $7.12\times10^{-5}$ m³/mol; divided by $c_P$: $2.44\times10^{-6}$ K/Pa $= 0.244$ K/bar.
H₂ at 300 K: $2a/RT = 1.985\times10^{-5}$, minus $b$ gives $-6.76\times10^{-6}$; $\mu_{JT} = -0.023$ K/bar.

(d) Nitrogen at room temperature is far below its inversion temperature, so every throttling step cools it ($\mu_{JT} > 0$ and $\Delta P < 0$ give $\Delta T < 0$), and Linde's counter-current exchanger compounds the cooling. Hydrogen at 300 K is above its inversion temperature: throttling *warms* it, and the regenerative loop would run the wrong way. It must first be precooled below its inversion temperature (in practice with liquid nitrogen) before throttling cools it. The vdW number 224 K is rough; the criterion, not the precise value, is what this calculation buys.

**B7.** (a) From (3.37) with $P_2/P_1 = 2$:
$$\frac{1}{T_2} = \frac{1}{373.15} - \frac{8.314}{40\,660}\ln 2 = 2.6799\times10^{-3} - 1.417\times10^{-4} = 2.5382\times10^{-3}\ \mathrm{K^{-1}},$$
so $T_2 = 394.0$ K $= 120.8$ °C. Food cooks faster because the water is 20 °C hotter.

(b) Take the transition ice → water. Then $L = +333.6$ kJ/kg (heat absorbed melting) and $\Delta v = v_{\text{water}} - v_{\text{ice}} = 1.0002\times10^{-3} - 1.0908\times10^{-3} = -9.06\times10^{-5}$ m³/kg, both in the same direction. By (3.36),
$$\frac{dP}{dT} = \frac{333.6\times10^3}{273.15\times(-9.06\times10^{-5})} = -1.348\times10^{7}\ \mathrm{Pa/K} = -134.8\ \mathrm{bar/K}.$$
Negative: raising the pressure by 135 bar lowers the melting point by 1 K. Had you taken $L$ for melting but $\Delta v$ for freezing, you would get a positive slope, which is wrong for water.

(c) A 70 kg skater on a blade contact of 3 mm × 10 cm exerts $P = 70\times9.81/(3\times10^{-4}) = 22.9$ bar, which lowers the melting point by $22.9/134.8 = 0.170$ K. Skating works on ice at −5 °C or −10 °C, so pressure melting cannot be the main mechanism. The usual explanations are frictional heating of the contact and a thin, liquid-like layer on the ice surface; Clausius–Clapeyron rules out the textbook story by two orders of magnitude.

**B8.** (a) From (3.16) with $U = \tfrac f2Nk_BT$:
$$F = U - TS = \tfrac f2Nk_BT - Nk_BT\left[\ln\frac VN + \frac f2\ln\!\left(\tfrac f2k_BT\right) + \sigma\right].$$
At fixed $T$, $V$, only the explicit $N$s vary. $\partial(N\ln(V/N))/\partial N = \ln(V/N) - 1$, so
$$\mu = \left(\frac{\partial F}{\partial N}\right)_{T,V} = \tfrac f2k_BT - k_BT\left[\ln\frac VN + \frac f2\ln\!\left(\tfrac f2k_BT\right) + \sigma\right] + k_BT.$$
Substitute $V/N = k_BT/P$ and split off a reference pressure $P_0$:
$$\mu = k_BT\ln\frac{P}{P_0} + \mu_0(T),\qquad \mu_0(T) = k_BT\left[1 + \tfrac f2 - \sigma - \ln\frac{k_BT}{P_0} - \tfrac f2\ln\!\left(\tfrac f2k_BT\right)\right],$$
the form (3.34). (The logarithms of dimensional quantities are made dimensionless by the units hidden in $\sigma$; SM-3 makes that explicit.)

(b) $G = F + PV = F + Nk_BT$. Multiply the $\mu$ of (a) by $N$: $N\mu = \tfrac f2Nk_BT - Nk_BT[\cdots] + Nk_BT = F + Nk_BT = G$. ✓

(c) At fixed $T$, from (a), $d\mu = k_BT\,dP/P$. Gibbs–Duhem (3.33) at $dT = 0$ says $N\,d\mu = V\,dP$, i.e. $d\mu = (V/N)\,dP = (k_BT/P)\,dP$. ✓

(d) Per mole, $\Delta\mu = RT\ln(P_2/P_1) = 8.314\times300\times\ln 10 = 5743$ J/mol. If the 10 atm and 1 atm samples are connected, particles flow from high $\mu$ to low $\mu$ (§3.10): from the compressed sample to the dilute one, which is also what pressure equalization says. For an ideal gas at one temperature "flows to lower $\mu$" and "flows to lower $P$" agree; with different temperatures or different substances they do not, and $\mu$ is the one that is right.

**B9.** The reference solution is at the end of `project.py`.

1. Call `coexistence(T)` on a grid $T = 240, 245, \dots, 300$ K; each call returns $(P^*, v_l, v_g)$.
2. Clausius–Clapeyron: $dP^*/dT \approx [P^*(T + 5) - P^*(T - 5)]/10$ at interior points, and $L_{cc} = T(v_g - v_l)\,dP^*/dT$.
3. Entropy route: from (3.28) at fixed $T$, $s_g - s_l = R\ln[(v_g - b)/(v_l - b)]$ (the $c_v\ln T$ terms cancel), so $L_s = RT\ln[(v_g - b)/(v_l - b)]$.
4. The script asserts $|L_{cc}/L_s - 1| < 2\times10^{-3}$ at every interior point. At 270 K, $P^* = 4.526$ MPa and $L = 4.81$ kJ/mol. At 240 K, $L = 6.27$ kJ/mol; at 300 K, $L = 1.73$ kJ/mol.

Why they must agree exactly: the equal-area rule is exact $\mu_l = \mu_g$ within the vdW model (§3.11). The derivation of (3.36) from $d\mu_l = d\mu_g$ and Gibbs–Duhem uses nothing but that equality, with the exact $v_l$, $v_g$. No ideal-vapour or negligible-liquid-volume approximation enters (those are only in (3.37)). And (3.28) is the exact entropy of the model. So the two routes are the same identity; the residual disagreement is the $O(h^2)$ error of the central difference.

As $T \to T_c$, $v_l$ and $v_g$ merge at $v_c = 3b$, so $\ln[(v_g - b)/(v_l - b)] \to 0$ and $L \to 0$: the script asserts $L$ decreases monotonically over the grid. The latent heat vanishes at the critical point because the two phases become the same phase. Treat the absolute values as qualitative: the vdW constants are fitted to the critical point, not to the liquid.

## C. Challenge

**C1.** (a) The engine and the bricks together form an isolated system, and the engine is reversible and cyclic, so the total entropy is unchanged:
$$C\ln\frac{T_f}{T_1} + C\ln\frac{T_f}{T_2} = 0 \implies T_f = \sqrt{T_1T_2} = \sqrt{400\times200} = 282.8\ \mathrm{K}.$$
Energy: the hot brick gives up $C(T_1 - T_f)$, the cold one receives $C(T_f - T_2)$, and the difference is the work:
$$W = C(T_1 + T_2 - 2T_f) = 1680\,(600 - 565.69) = 57\,648\ \mathrm{J}.$$

(b) $T_1 + T_2 - 2\sqrt{T_1T_2} = (\sqrt{T_1} - \sqrt{T_2})^2$, so $W_{\max} = C(\sqrt{T_1} - \sqrt{T_2})^2$, which is positive for $T_1 \ne T_2$ and zero only when there is no temperature difference to exploit. The script checks the two forms agree (57 648 J).

(c) Warming both bricks from $\sqrt{T_1T_2}$ to $(T_1 + T_2)/2 = 300$ K takes $2C[(T_1 + T_2)/2 - \sqrt{T_1T_2}] = C(T_1 + T_2 - 2\sqrt{T_1T_2}) = W_{\max} = 57\,648$ J. If that heat is supplied by dissipating the work (stirring, a resistor), each brick's entropy rises by $C\ln[300/282.8]$, total
$$2C\ln\frac{(T_1 + T_2)/2}{\sqrt{T_1T_2}} = C\ln\frac{(T_1 + T_2)^2}{4T_1T_2} = C\ln\frac{360\,000}{320\,000} = C\ln\frac98 = 197.9\ \mathrm{J/K},$$
exactly the entropy generated in Example 4.3. So direct contact is equivalent to "run the reversible engine, then throw its work back into the bricks as heat": the work that could have been extracted is what was lost, and its dissipation is where the entropy came from.

**C2.** (a) Use the energy form. At fixed total $S$, $V$, $N$, a stable equilibrium has minimum $U$. (If some rearrangement at fixed $S$ lowered $U$, you could extract the difference as work and return it as heat, ending at the original $U$ with higher $S$, contradicting the entropy maximum.) Split the system into two equal halves and let a fluctuation move entropy $\sigma$ and volume $\nu$ from one half to the other. With $U_h$ the energy function of a half,
$$U_h\!\left(\tfrac S2 + \sigma, \tfrac V2 + \nu\right) + U_h\!\left(\tfrac S2 - \sigma, \tfrac V2 - \nu\right) = 2U_h + \begin{pmatrix}\sigma & \nu\end{pmatrix}\mathbf H\begin{pmatrix}\sigma\\ \nu\end{pmatrix} + \cdots,$$
where $\mathbf H$ is the Hessian of $U_h$ (the linear terms cancel). The uniform state is a minimum for every direction $(\sigma, \nu)$ only if $\mathbf H$ is positive semidefinite: $U(S, V)$ is convex, equivalently $S(U, V)$ concave (the same argument in entropy form gives a negative semidefinite Hessian of $S$). A $2\times2$ symmetric matrix is positive semidefinite iff its diagonal entries and its determinant are $\ge 0$ (LA-8, positive semidefiniteness of a quadratic form: both eigenvalues nonnegative; their product is the determinant and their sum the trace).

(b) Diagonal entry: $U_{SS} = (\partial T/\partial S)_V = T/C_V \ge 0$, so $C_V > 0$ (as in Theorem 3.6).

Determinant: $\mathbf H$ is the derivative matrix of the map $(S, V) \mapsto (U_S, U_V) = (T, -P)$, so $\det\mathbf H = \partial(T, -P)/\partial(S, V) = -\partial(T, P)/\partial(S, V)$, where $\partial(f, g)/\partial(x, y)$ is the Jacobian determinant. Factor the map through $(T, V)$: $(S, V) \mapsto (T, V) \mapsto (T, P)$. Determinants of composed maps multiply (chain rule: the derivative matrix of a composition is the product of derivative matrices). The two factors are triangular:
$$\frac{\partial(T, V)}{\partial(S, V)} = \det\begin{pmatrix}(\partial T/\partial S)_V & (\partial T/\partial V)_S\\ 0 & 1\end{pmatrix} = \frac{T}{C_V},\qquad \frac{\partial(T, P)}{\partial(T, V)} = \det\begin{pmatrix}1 & 0\\ (\partial P/\partial T)_V & (\partial P/\partial V)_T\end{pmatrix} = \left(\frac{\partial P}{\partial V}\right)_T.$$
So $\det\mathbf H = -\dfrac{T}{C_V}\left(\dfrac{\partial P}{\partial V}\right)_T = \dfrac{T}{C_V\,V\kappa_T}$. With $C_V > 0$, $\det\mathbf H \ge 0$ requires $\kappa_T > 0$.

(c) With $\kappa_T > 0$, (3.29) gives $C_P - C_V = TV\alpha^2/\kappa_T \ge 0$, with equality only where $\alpha = 0$. On the part of a vdW isotherm where $(\partial P/\partial v)_T > 0$, $\kappa_T < 0$, so $\det\mathbf H < 0$: there is a direction of fluctuation that lowers $U$ at fixed total $S$, $V$, and a uniform fluid there splits into two phases spontaneously. That is the mechanical instability the Maxwell construction of §3.11 replaces with the flat coexistence line.

**C3.** (a) Furnace: 0.95 J of heat per joule of fuel. Heat pump: the plant turns 1 J of fuel into 0.40 J of electricity, and the heat pump turns that into $3\times0.40 = 1.20$ J of heat. The heat pump wins by about 26% even after paying for a fossil-fuelled plant, because it moves outdoor heat rather than making heat.

(b) Against an ideal furnace (1 J per J), break-even is $0.40\times\mathrm{COP} = 1$, so $\mathrm{COP} = 1/0.40 = 2.50$.

(c) By (3.18), $\mathrm{COP}_{\text{hp}} \le 293.15/(293.15 - 263.15) = 293.15/30 = 9.77$. Real units deliver far less because: the refrigerant must be colder than the outdoor air to absorb heat and hotter than the indoor air to release it, so the cycle runs across a temperature span perhaps twice the 30 K you computed (and COP falls roughly as the inverse of that span); the compressor is irreversible; the expansion valve is a Joule–Thomson throttle, irreversible by design; and below freezing the outdoor coil frosts and must be periodically defrosted with heat from indoors.

(d) Set $0.45\,T_{\text{in}}/(T_{\text{in}} - T_{\text{out}}) = 2.5$ with $T_{\text{in}} = 293.15$ K: $T_{\text{in}} - T_{\text{out}} = 0.45\times293.15/2.5 = 52.77$ K, so $T_{\text{out}} = 240.4$ K $= -32.8$ °C. In the script, extend Your-turn item 3 with `brentq(lambda To: 0.45*293.15/(293.15 - To) - 2.5, 200, 290)` and mark the root on the plot. On this crude model, a heat pump on fossil electricity beats even a perfect furnace down to about −33 °C; with low-carbon electricity the comparison is not close.
