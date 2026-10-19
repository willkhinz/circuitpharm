# Results Brief 15: Fifth Observable (Mean Open Time) & Testable Cooperativity Fit

**Date:** 2026-10-17  
**Branch:** `claude/agent-brief-15`  
**Worktree:** `/Users/whinz/Biochemistry-wt-agent-brief-15`  
**Thread:** `brief15`

---

## 1. Executive Summary & Plain Statement on Survival

In Brief 14, `kinetic_coop` was shown to satisfy four macroscopic/microscopic observables simultaneously ($\le 2.39\sigma$), but with 9 parameters fitted against 9 data points, the residual degrees of freedom were exactly 0 (an algebraic existence demonstration carrying zero statistical testing power). Furthermore, a prospective $k_{\text{off}}$ interval ($0.0848\text{--}0.1421\ \text{ms}^{-1}$) derived from the Jahn 1997 Peak CRC appeared to conflict with the fitted $k_{\text{off}} = 75.55\ \text{ms}^{-1}$.

Brief 15 resolves both issues decisively:
1. **Mean Open Time as the 5th Observable:** Mean open duration was derived rigorously from the infinitesimal generator matrix $Q$ for all Markov schemes (Colquhoun & Hawkes 1981, 1982 burst theory). Keramidas & Harrison 2008 Table III M-mode mean open time ($2.68$ ms, unstated error, empirical band $[2.54, 2.83]$ ms) was added as the 5th observable, increasing total data constraints to $n_{\text{data}} = 10$.
2. **Plain Statement on Tested Fit Survival:**
   **YES, `kinetic_coop` SURVIVES AS A TESTED FIT.**  
   With $n_{\text{data}} = 10$ and $k_{\text{fitted}} = 9$, `kinetic_coop` operates with **residual degrees of freedom $\text{res\_dof} = +1$**. It achieves a Five-Way Minimax Optimum of **$2.58\sigma$** (all five observables $\le 2.58\sigma < 3.0\sigma$) and a Natural Joint Optimum of **$2.95\sigma$** ($\le 2.95\sigma < 3.0\sigma$). Because $\text{res\_dof} = +1 \ge +1$, this is a genuine, over-determined empirical hypothesis test rather than an exact saturation.
3. **Symmetric Comparison:**
   - `kinetic_jw95_free_rates`: $k=6$, $\text{res\_dof} = +4$. Decisively rejected at **$14.25\sigma$** (cannot fit deactivation components and CRC slope simultaneously).
   - `kinetic_two_open`: $k=8$, $\text{res\_dof} = +2$. Fails at **$3.23\sigma$** (misses $3.0\sigma$ threshold; lacks the binding cooperativity needed to reconcile fast deactivation with steep Hill slope).
   - `kinetic_coop`: $k=9$, $\text{res\_dof} = \mathbf{+1}$. Passes at **$2.58\sigma$** ($P_{o,\text{max}} = 0.775$, $\tau_O = 3.20$ ms).
4. **Armed Nesting Assertion:**
   Strict nesting assertion `coop <= two_open + 1e-4` was armed and **PASSED**:
   - `obj_natural_joint`: coop $17.882791 \le$ two_open $25.425960$.
   - `obj_fiveway_minimax`: coop $2.578213 \le$ two_open $3.231843$.
5. **$k_{\text{off}}$ Non-Transfer Reasoning:**
   The $0.0848\text{--}0.1421\ \text{ms}^{-1}$ interval is an artifact of propagating the Jahn 1997 CRC through the 5-state Del Castillo-Katz (JW95) topology. It is structurally non-transferable to `kinetic_coop`, which modifies the unbinding step through binding cooperativity ($c \approx 1000$). At its fitted rate ($k_{\text{off}} = 85.67\ \text{ms}^{-1}$, $k_{\text{off2}} = 0.0857\ \text{ms}^{-1}$), `kinetic_coop` directly reproduces the raw Jahn 1997 Peak CRC ($z_{\text{EC50}} = +2.24\sigma$, $z_{n_H} = -2.55\sigma$). Raw data takes precedence over topology-bound inferred rates.
6. **Cooperativity Literature Register:**
   All three literature sources confirming the physical existence and positive direction of agonist binding cooperativity (Olsen et al. 1981, Agey & Dunn 1989, Mozrzymas et al. 2003) are registered in `src/circuitpharm/provenance.py` (with `where="abstract"`) and `data/pharmacology.db` (`sources` table). In strict accordance with negative rules, none are used as numerical fitting constraints.

---

## 2. Derivation of Mean Open Time from Generator Matrix $Q$

### 2.1 Theoretical Framework (Colquhoun & Hawkes 1981, 1982)
Let the state space of a continuous-time Markov chain be partitioned into open conducting states $\mathcal{O}$ and closed shut states $\mathcal{C}$. The infinitesimal transition rate matrix $Q$ (where row sums are zero) is block-partitioned as:
$$Q = \begin{pmatrix} Q_{\mathcal{O}\mathcal{O}} & Q_{\mathcal{O}\mathcal{C}} \\ Q_{\mathcal{C}\mathcal{O}} & Q_{\mathcal{C}\mathcal{C}} \end{pmatrix}$$
where $Q_{\mathcal{O}\mathcal{O}}$ governs transitions among open states, and $Q_{\mathcal{O}\mathcal{C}}$ governs exits from open states to shut states.

Upon channel opening from $\mathcal{C}$, initial probability across open states is given by entry flux:
$$\pi_{\mathcal{O}} = \frac{p_{\mathcal{C}}^* Q_{\mathcal{C}\mathcal{O}}}{p_{\mathcal{C}}^* Q_{\mathcal{C}\mathcal{O}} \mathbf{1}}$$
For all schemes in this project (`kinetic_jw95`, `kinetic_two_open`, `kinetic_coop`, `kinetic_flip`, `kinetic_open_desens`), channel opening occurs strictly through a single gateway transition from $A2R$ (or $A2F$) to $A2O$ with rate $\beta$. No shut state directly connects to secondary open states ($A2O'$). Thus, entry flux is 100% concentrated in $A2O$:
$$\pi_{\mathcal{O}} = (1, 0, \dots, 0)$$

The open lifetime $T_O$ (time until first exit to $\mathcal{C}$) has survival function:
$$S(t) = P(T_O > t) = \pi_{\mathcal{O}} \exp(Q_{\mathcal{O}\mathcal{O}} t) \mathbf{1}$$
Integrating over $t \in [0, \infty)$ yields the exact expected mean open duration:
$$\langle \tau_O \rangle = \int_0^\infty S(t)\, dt = \pi_{\mathcal{O}} (-Q_{\mathcal{O}\mathcal{O}})^{-1} \mathbf{1}$$

### 2.2 Algebraic Reduction for 1-Open-State Schemes
For schemes with a single conducting state ($A2O$):
- $\mathcal{O} = \{A2O\}$.
- $Q_{\mathcal{O}\mathcal{O}} = [-\alpha]$ (or $-(\alpha+d)$ for `kinetic_open_desens`).
- $\pi_{\mathcal{O}} = [1]$.
$$\langle \tau_O \rangle = \frac{1}{\alpha} \quad \left(\text{or } \frac{1}{\alpha + d}\right)$$

### 2.3 Algebraic Reduction for 2-Open-State Schemes (`kinetic_two_open`, `kinetic_coop`)
For schemes with sequential open states:
$$A2R \underset{\alpha}{\overset{\beta}{\rightleftharpoons}} A2O \underset{\gamma}{\overset{\delta}{\rightleftharpoons}} A2O'$$
The open sub-generator is:
$$Q_{\mathcal{O}\mathcal{O}} = \begin{pmatrix} -(\alpha + \delta) & \delta \\ \gamma & -\gamma \end{pmatrix} \implies -Q_{\mathcal{O}\mathcal{O}} = \begin{pmatrix} \alpha + \delta & -\delta \\ -\gamma & \gamma \end{pmatrix}$$
The determinant is:
$$\det(-Q_{\mathcal{O}\mathcal{O}}) = (\alpha + \delta)\gamma - (-\delta)(-\gamma) = \alpha \gamma + \delta \gamma - \delta \gamma = \alpha \gamma$$
The fundamental matrix of open states is:
$$(-Q_{\mathcal{O}\mathcal{O}})^{-1} = \frac{1}{\alpha \gamma} \begin{pmatrix} \gamma & \delta \\ \gamma & \alpha + \delta \end{pmatrix} = \begin{pmatrix} \frac{1}{\alpha} & \frac{\delta}{\alpha \gamma} \\ \frac{1}{\alpha} & \frac{\alpha + \delta}{\alpha \gamma} \end{pmatrix}$$
Multiplying by initial entry distribution $\pi_{\mathcal{O}} = (1, 0)$:
$$\pi_{\mathcal{O}} (-Q_{\mathcal{O}\mathcal{O}})^{-1} = \left( \frac{1}{\alpha},\ \frac{\delta}{\alpha \gamma} \right)$$
Summing over all open states ($\cdot \mathbf{1}$):
$$\langle \tau_O \rangle = \frac{1}{\alpha} + \frac{\delta}{\alpha \gamma} = \frac{1}{\alpha} \left( 1 + \frac{\delta}{\gamma} \right)$$

**Physical Interpretation:**
- $1/\alpha$ is the mean residence time spent in $A2O$ before closing back to $A2R$.
- $\delta / \alpha$ is the expected number of excursions into $A2O'$ before closing.
- $1/\gamma$ is the mean dwell time of each individual visit to $A2O'$.
- The product $(\delta / \alpha)(1/\gamma)$ is the total expected time spent in $A2O'$.
- As $\delta \to 0$ or $\gamma \to \infty$, $\langle \tau_O \rangle \to 1/\alpha$, strictly preserving structural nesting.

Unit tests in `tests/test_mean_open_time.py` confirm this analytical formula matches direct numerical inversion of $-Q_{\mathcal{O}\mathcal{O}}$ to machine precision ($< 10^{-12}$).

---

## 3. Keramidas 2008 & Modal Gating vs Two-Open-States Distinction

### 3.1 Scoring Without an Invented $\sigma$
Keramidas & Harrison (2008, *J Gen Physiol* 131:163-181) measured single-channel currents of rat $\alpha1\beta2\gamma2S$ in outside-out macropatches at saturating (10 mM) GABA. Table III reports exponential open-time components:
- **M-mode (Medium $P_o$):** Observed in 11 of 15 patches. Components: $\tau = (0.49 \pm 0.02, 2.58 \pm 0.09, 5.17 \pm 0.38)$ ms with relative areas $(0.26, 0.49, 0.25)$. The area-weighted mean open time is **$2.68$ ms**. Propagating component errors yields the interval **$[2.54, 2.83]$ ms**.
- **H-mode (High $P_o$):** Observed in 5 patches. Components: $\tau = (0.59 \pm 0.11, 4.22 \pm 0.73, 13.2 \pm 2.6)$ ms with relative areas $(0.18, 0.41, 0.41)$. Weighted mean is **$7.25$ ms**, spanning **$[5.86, 8.63]$ ms**.

**Crucial Reporting Integrity:** Keramidas 2008 states **NO error statistic** (whether $\pm$ denotes SEM or SD is unstated in the text). In accordance with repository rules, we do not invent an error bar or assume SEM. Instead, the observation is evaluated as a **flat band constraint** ($[2.54, 2.83]$ ms for M-mode):
- If the predicted mean open time lies within $[2.54, 2.83]$ ms, penalty = 0 ($z = 0.0\sigma$).
- If outside, residual is measured from the nearest boundary scaled by half-width $\sigma_{\text{eff}} = 0.145$ ms.
- For completeness, scoring against the midpoint ($2.68 \pm 0.145$ ms) is also reported.

### 3.2 Modal Gating vs Two Open States
A critical biophysical distinction must be maintained:
- **Modal Gating (Keramidas M-mode vs H-mode):** Channels switch between M-mode and H-mode on a timescale of **seconds to tens of seconds**. Each mode represents a distinct, meta-stable gating regime. Crucially, **each mode itself exhibits three exponential open components**.
- **Two Sequential Open States ($A2O \rightleftharpoons A2O'$):** The two open states in `kinetic_two_open` and `kinetic_coop` represent sub-millisecond conformational transitions occurring **within a single burst**.
- **Conclusion:** M-mode and H-mode are **NOT** the microscopic states $A2O$ and $A2O'$. Conflating them would manufacture artificial agreement by mapping a seconds-long modal switch onto a sub-millisecond gating step. Because M-mode represents 73% of patches (11/15), its weighted mean open duration ($2.68$ ms) is the appropriate reference for the single-channel open dwell time.

### 3.3 Isoform Distinction: Li 2008 ($\gamma2L$) vs Keramidas 2008 ($\gamma2S$)
Li et al. (2008 / 2007b control) reported intracluster open times on rat $\alpha1\beta2\gamma2L$ (cell-attached, sub-saturating 50 $\mu$M GABA) with weighted mean $2.96$ ms (SD, $n=4$). This $\gamma2L$ measurement is registered with explicit isoform and splice-variant tags, maintaining documentation integrity.

---

## 4. Symmetric Five-Way Fit Evaluation (Task 2)

### 4.1 Data Constraints ($n_{\text{data}} = 10$)
1. Deactivation fast tau ($\tau_1 = 2.8 \pm 0.3$ ms, Barberis 2007)
2. Deactivation intermediate tau ($\tau_2 = 33.4 \pm 4.6$ ms)
3. Deactivation slow tau ($\tau_3 = 221.35 \pm 14.9$ ms)
4. Deactivation fast area ($a_1 = 0.57 \pm 0.04$)
5. Deactivation intermediate area ($a_2 = 0.23 \pm 0.02$)
   *(Areas scored under rank-2 covariance from 3 published errors; $a_3 = 1 - a_1 - a_2$)*
6. Jahn 1997 Peak CRC $EC_{50}$ ($11.6 \pm 0.9\ \mu$M)
7. Jahn 1997 Peak CRC Hill slope $n_H$ ($2.2 \pm 0.4$ over $1\text{--}10\ \mu$M)
8. Barberis 2007 Paired-Pulse Recovery ($0.33 \pm 0.03$ at 100 ms gap)
9. Gating anchor $P_{o,\text{max}}$ ($0.75 \pm 0.05$, band $[0.70, 0.80]$)
10. Keramidas 2008 Table III M-mode mean open time ($2.68$ ms, band $[2.54, 2.83]$ ms)

### 4.2 Full Per-Dataset Residuals Summary Table (in $\sigma$)

| Scheme | Evaluation Profile | $z_{\tau_1}$ | $z_{\tau_2}$ | $z_{\tau_3}$ | $z_{a_1}$ | $z_{a_2}$ | $z_{a_3}$ | $\sqrt{\chi^2_a}$ | $z_{\text{EC50}}$ | $z_{n_H}$ | $z_{pp}$ | $P_{o,max}$ | $\tau_O$ (ms) | $z_{\text{MOT}}$ (band) | $z_{\text{MOT}}$ (mid) | Max Viol / $3\sigma$ | Satisfies All 5 $\le 3\sigma$? |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `kinetic_jw95_free_rates` | Deact-Alone Optimum | -8.10$\sigma$ | +0.02$\sigma$ | +0.06$\sigma$ | -14.43$\sigma$ | +9.91$\sigma$ | +12.64$\sigma$ | 14.43$\sigma$ | -12.49$\sigma$ | -3.48$\sigma$ | +2.74$\sigma$ | 0.100 (-13.01$\sigma$) | 0.41 | -14.69$\sigma$ | -15.65$\sigma$ | **4.90** | **NO** |
| `kinetic_jw95_free_rates` | CRC-Alone Optimum | -9.33$\sigma$ | +37.05$\sigma$ | +50.99$\sigma$ | -14.25$\sigma$ | +38.10$\sigma$ | -6.40$\sigma$ | 41.50$\sigma$ | +0.08$\sigma$ | -0.47$\sigma$ | -0.56$\sigma$ | 0.998 (+4.96$\sigma$) | 0.02 | -17.38$\sigma$ | -18.34$\sigma$ | **17.00** | **NO** |
| `kinetic_jw95_free_rates` | Natural Joint Optimum | -8.03$\sigma$ | -0.33$\sigma$ | +1.28$\sigma$ | -14.53$\sigma$ | +9.82$\sigma$ | +12.82$\sigma$ | 14.53$\sigma$ | +0.03$\sigma$ | -1.48$\sigma$ | +2.03$\sigma$ | 0.838 (+1.76$\sigma$) | 2.50 | -0.28$\sigma$ | -1.24$\sigma$ | **4.84** | **NO** |
| `kinetic_jw95_free_rates` | Five-Way Minimax Optimum | -9.33$\sigma$ | +11.33$\sigma$ | +14.25$\sigma$ | -14.25$\sigma$ | +9.80$\sigma$ | +12.47$\sigma$ | 14.25$\sigma$ | +13.25$\sigma$ | -0.52$\sigma$ | -2.63$\sigma$ | 0.999 (+4.97$\sigma$) | 0.57 | -13.60$\sigma$ | -14.57$\sigma$ | **4.75** | **NO** |
| `kinetic_two_open` | Deact-Alone Optimum | -0.00$\sigma$ | +0.00$\sigma$ | +0.00$\sigma$ | +0.12$\sigma$ | +0.24$\sigma$ | +0.16$\sigma$ | 0.00$\sigma$ | -5.85$\sigma$ | -3.65$\sigma$ | +3.71$\sigma$ | 0.880 (+2.61$\sigma$) | 2.40 | -0.98$\sigma$ | -1.94$\sigma$ | **1.95** | **NO** |
| `kinetic_two_open` | CRC-Alone Optimum | -9.33$\sigma$ | +131.53$\sigma$ | +4340.18$\sigma$ | -14.25$\sigma$ | -11.50$\sigma$ | +26.67$\sigma$ | 32.61$\sigma$ | -0.01$\sigma$ | -0.28$\sigma$ | -10.99$\sigma$ | 1.000 (+5.00$\sigma$) | 9728.13 | +67071.03$\sigma$ | +67072.06$\sigma$ | **22357.01** | **NO** |
| `kinetic_two_open` | Natural Joint Optimum | -0.42$\sigma$ | -1.23$\sigma$ | +1.51$\sigma$ | +0.50$\sigma$ | +0.24$\sigma$ | +0.22$\sigma$ | 0.30$\sigma$ | -0.10$\sigma$ | -3.89$\sigma$ | +2.46$\sigma$ | 0.770 (+0.41$\sigma$) | 2.84 | +0.04$\sigma$ | +1.07$\sigma$ | **1.30** | **NO** |
| `kinetic_two_open` | Five-Way Minimax Optimum | -3.18$\sigma$ | +0.06$\sigma$ | +3.23$\sigma$ | +3.27$\sigma$ | -2.46$\sigma$ | -2.10$\sigma$ | 3.23$\sigma$ | -3.23$\sigma$ | -3.23$\sigma$ | +0.98$\sigma$ | 0.862 (+2.25$\sigma$) | 2.22 | -2.22$\sigma$ | -3.19$\sigma$ | **1.08** | **NO** |
| `kinetic_coop` | Deact-Alone Optimum | +0.00$\sigma$ | +0.00$\sigma$ | +0.00$\sigma$ | +0.12$\sigma$ | +0.24$\sigma$ | +0.16$\sigma$ | 0.00$\sigma$ | -5.85$\sigma$ | -3.65$\sigma$ | +3.71$\sigma$ | 0.880 (+2.61$\sigma$) | 2.40 | -0.98$\sigma$ | -1.94$\sigma$ | **1.95** | **NO** |
| `kinetic_coop` | CRC-Alone Optimum | -9.33$\sigma$ | +131.53$\sigma$ | +4340.17$\sigma$ | -14.25$\sigma$ | -11.50$\sigma$ | +26.67$\sigma$ | 32.61$\sigma$ | -0.01$\sigma$ | -0.28$\sigma$ | -10.99$\sigma$ | 1.000 (+5.00$\sigma$) | 9728.10 | +67070.83$\sigma$ | +67071.87$\sigma$ | **22356.94** | **NO** |
| `kinetic_coop` | Natural Joint Optimum | -0.20$\sigma$ | -0.93$\sigma$ | +1.44$\sigma$ | +0.09$\sigma$ | +0.50$\sigma$ | +0.67$\sigma$ | 0.33$\sigma$ | +0.16$\sigma$ | -2.95$\sigma$ | +2.46$\sigma$ | 0.756 (+0.12$\sigma$) | 2.84 | +0.05$\sigma$ | +1.08$\sigma$ | **0.98** | **YES** |
| `kinetic_coop` | Five-Way Minimax Optimum | +0.50$\sigma$ | -0.24$\sigma$ | +2.50$\sigma$ | -0.96$\sigma$ | +3.08$\sigma$ | +0.40$\sigma$ | 2.58$\sigma$ | +2.24$\sigma$ | -2.55$\sigma$ | +2.58$\sigma$ | 0.775 (+0.50$\sigma$) | 3.20 | +2.58$\sigma$ | +3.61$\sigma$ | **0.86** | **YES** |

### 4.3 Optimal Parameter Vectors & Effective Parameter Count ($k_{\text{eff}}$)

| Scheme | Profile | Free Parameters Vector ($\log_{10}$) | Bounds Status / Pinned Parameters | $k_{fitted}$ | $k_{eff}$ |
| :--- | :--- | :--- | :--- | :---: | :---: |
| `kinetic_jw95_free_rates` | Natural Joint Optimum | `log10_kd=1.105, log10_E=0.714, log10_D=0.498, log10_r=-1.708, log10_koff=-1.673, log10_alpha=-0.398` | None (all interior) | 6 | 6 |
| `kinetic_jw95_free_rates` | Five-Way Minimax Optimum | `log10_kd=2.669, log10_E=2.897, log10_D=2.468, log10_r=-2.116, log10_koff=0.145, log10_alpha=0.246` | None (all interior) | 6 | 6 |
| `kinetic_two_open` | Natural Joint Optimum | `log10_kd=0.540, log10_E=0.033, log10_D=1.505, log10_r=-1.766, log10_koff=-1.018, log10_alpha=0.040, log10_delta=-1.050, log10_gamma=-1.374` | None (all interior) | 8 | 8 |
| `kinetic_two_open` | Five-Way Minimax Optimum | `log10_kd=0.523, log10_E=0.292, log10_D=1.999, log10_r=-1.882, log10_koff=-0.564, log10_alpha=0.159, log10_delta=-1.142, log10_gamma=-1.484` | None (all interior) | 8 | 8 |
| `kinetic_coop` | Natural Joint Optimum | `log10_kd=2.489, log10_E=-0.032, log10_D=1.452, log10_r=-1.772, log10_koff=1.935, log10_alpha=0.070, log10_delta=-1.028, log10_gamma=-1.396, log10_c=3.000` | log10_c on upper bound (3.0) | 9 | 8 |
| `kinetic_coop` | Five-Way Minimax Optimum | `log10_kd=2.603, log10_E=-0.036, log10_D=1.458, log10_r=-1.851, log10_koff=1.933, log10_alpha=0.068, log10_delta=-0.989, log10_gamma=-1.428, log10_c=3.000` | log10_c on upper bound (3.0) | 9 | 8 |

---

## 5. $k_{\text{off}}$ Non-Transfer Reasoning ($0.0848\text{--}0.1421$ vs $75.55\ \text{ms}^{-1}$)

A major concern raised prior to Brief 15 was the apparent conflict between two numbers:
1. **CRC-propagated $k_{\text{off}}$ interval:** $0.0848\text{--}0.1421\ \text{ms}^{-1}$ (derived from the Jahn 1997 Peak CRC).
2. **`kinetic_coop` fitted $k_{\text{off}}$:** $75.55\ \text{ms}^{-1}$ (in Brief 14) and $85.67\ \text{ms}^{-1}$ (in Brief 15).

This difference ($530\times$ to $600\times$) looks decisive at first glance, but **it does not falsify `kinetic_coop`**. The complete physical and mathematical rationale is as follows:

1. **Topology Dependence of Parameter Inversion:**
   The $0.0848\text{--}0.1421\ \text{ms}^{-1}$ interval is not a direct measurement of unbinding. There is no published radioligand or single-channel measurement reporting $k_{\text{off}}$ for $\alpha1\beta2\gamma2S$ in the literature. That interval was **derived by inverting the algebraic equations of the 5-state Del Castillo-Katz (JW95) topology**. In JW95, binding sites are non-cooperative and identical ($k_{\text{off1}} = 2 k_{\text{off}}$, $k_{\text{off2}} = k_{\text{off}}$), and there is only one open state. Under those strict topology assumptions, matching an $EC_{50}$ of $11.6\ \mu$M while constraining gating efficacy forces $k_{\text{off}}$ into the $0.0848\text{--}0.1421\ \text{ms}^{-1}$ window.
2. **Invalid Cross-Topology Parameter Transfer:**
   Transferring a rate constant derived under Scheme A as an empirical constraint on Scheme B is valid *only if Scheme B preserves the exact mathematical mapping between that rate and the raw data*.
   `kinetic_coop` radically alters this mapping by introducing binding cooperativity $c = k_{\text{off}} / k_{\text{off2}} \approx 1000$. The unbinding of singly-bound GABA is fast ($k_{\text{off}} \approx 85.7\ \text{ms}^{-1}$), but once a second GABA binds, unbinding is slowed by three orders of magnitude ($k_{\text{off2}} \approx 0.0857\ \text{ms}^{-1}$). Notice that $k_{\text{off2}} = 0.0857\ \text{ms}^{-1}$ sits **directly inside the $0.0848\text{--}0.1421\ \text{ms}^{-1}$ window**!
3. **Direct Reproduction of the Raw Data:**
   What matters for model validity is whether the model reproduces the raw experimental observation (the Jahn 1997 concentration-response curve), NOT whether it conforms to a parameter value inferred under a rejected scheme.
   At $k_{\text{off}} = 85.67\ \text{ms}^{-1}$, `kinetic_coop` reproduces the Jahn 1997 Peak CRC directly:
   $$z_{\text{EC50}} = +2.24\sigma \quad (\text{measured } 13.6\ \mu\text{M vs } 11.6 \pm 0.9\ \mu\text{M})$$
   $$z_{n_H} = -2.55\sigma \quad (\text{measured } 1.18 \text{ vs } 2.2 \pm 0.4)$$
   Because `kinetic_coop` directly satisfies the raw CRC data within $3\sigma$, enforcing the JW95-propagated $k_{\text{off}}$ would be circular reasoning: it would reject a model that fits the data because it fails to obey the internal mechanics of a model that fails the data ($14.25\sigma$).

---

## 6. Cooperativity Literature Register (Task 3)

All three literature sources discovered during Brief 14 are registered with complete bibliographic identifiers, abstracts, and notes:

### 6.1 Registered Records
1. **Olsen et al. (1981):**
   - Citation: Olsen RW, Bergman MO, Van Ness PC, Lummis SC, Watkins AE, Napias C, Greenlee DV. $\gamma$-Aminobutyric acid receptor binding in mammalian brain. *Mol Pharmacol* 19(2):217-227.
   - Identifiers: PMID `6262615`, DOI `10.1016/s0026-895x(25)13464-0`.
   - Read: `ABSTRACT_PUBMED`, `where="abstract"`.
   - Content: Equilibrium $[^3\text{H}]\text{GABA}$ radioligand binding in mammalian brain membranes shows curvilinear Scatchard plots and Hill coefficients consistent with positive cooperativity / multiple affinity binding sites ($K_D$ classes $\sim 13$ nM and $\sim 300$ nM, $\sim 23\times$).
2. **Agey & Dunn (1989):**
   - Citation: Agey MW, Dunn MW. Kinetics of $[^3\text{H}]\gamma$-aminobutyric acid binding to membrane-bound and solubilized receptors from rat brain. *Biochemistry* 28(10):4200-4208.
   - Identifiers: PMID `2548571`, DOI `10.1021/bi00436a012`.
   - Read: `ABSTRACT_PUBMED`, `where="abstract"`.
   - Content: Rapid filtration binding kinetics in rat brain membranes demonstrating biphasic association and non-first-order dissociation kinetics, physically demonstrating cooperative interactions between agonist binding sites.
3. **Mozrzymas et al. (2003):**
   - Citation: Mozrzymas JW, Barberis A, Mercik K, Zarnowska ED. Binding sites, gating, and desensitization in $\alpha1\beta2\gamma2$ GABA-A receptors: experimental and modeling study. *J Neurophysiol* 89(2):871-883.
   - Identifiers: PMID `12574465`, DOI `10.1152/jn.00951.2002`.
   - Read: `ABSTRACT_PUBMED`, `where="abstract"`.
   - Content: Experimental and modeling investigation of recombinant rat $\alpha1\beta2\gamma2$ confirming that the concentration dependence of current rising phases requires positive cooperativity between agonist binding steps.

### 6.2 Strict Negative Rules Enforced
- None of these sources are used as numerical fitting constraints or Likelihood datasets.
- In `src/circuitpharm/provenance.py`, they are placed in `COOPERATIVITY_LITERATURE_PROV` and strictly excluded from parameter tables in `ALL`.
- In `data/pharmacology.db`, they are recorded with `HAND: METADATA_RESOLVED_MATCH` and explicit `CLAIM_SUPPORT_DOES_NOT_SUPPORT` statements forbidding parameter constraints.
- Verified by unit tests in `tests/test_provenance.py` (`test_cooperativity_literature_provenance`).

---

## 7. Verification & Test Suite Execution

Both tiers of pytest were executed and verified:
1. Fast test suite:
   ```bash
   PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python -m pytest -q -m "not slow"
   ```
2. Full test suite:
   ```bash
   PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python -m pytest -q
   ```

---

## 8. Verbatim Test Suite Outputs

### 8.1 Fast Test Suite (`not slow`)
Command:
```bash
PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python -m pytest -q -m "not slow"
```
Verbatim final summary output:
```
587 passed, 50 warnings in 116.74s (0:01:56)
```

### 8.2 Full Test Suite (All Tests)
Command:
```bash
PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python -m pytest -q
```
Verbatim final summary output:
```
661 passed, 133 warnings in 729.13s (0:12:09)
```

