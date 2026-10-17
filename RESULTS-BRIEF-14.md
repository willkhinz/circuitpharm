# Results Brief 14 — Putting Recovery in the Objective, Sourced Anchor Resolution, and Nesting Restoration

**Worktree:** `/Users/whinz/Biochemistry-wt-agent-brief-14`  
**Branch:** `claude/agent-brief-14`  
**Base Commit on `main`:** `918cea9`  
**Task 3 Commit:** `e67c1f1` (`Task 3: Update manuscript and tests to sourced 11.6 uM anchor with dated note in S3 and derived falsification threshold`)  
**Task 1 Commit:** `be8bba2` (`Task 1: Add recovery to joint objective, true four-way minimax, and enforce nesting warm-start and assertion`)  
**Date:** 2026-10-17  
**Author:** Antigravity Worker (Brief 14)  
**Reporting to:** Claude Opus 5 (Planner) via Agent Bridge (`thread: brief14`)  
**Reference Brief:** `AGENT-BRIEF-14.md`  

---

## 1. Executive Summary & Answers to the Four Core Questions

In Briefs 12 and 13, the joint optimization procedure fitted only three observables (deactivation, CRC, and peak open probability) and evaluated paired-pulse recovery only *downstream* after optimization. Consequently, the reported 4.34σ for `kinetic_two_open` was an upper bound, and the reported 17.69σ for `kinetic_coop` was an artifact of an unconstrained recovery explosion.

In Brief 14, paired-pulse recovery was integrated directly into `obj_natural_joint` using `BARBERIS2007_PAIRED_PULSE` (0.33 ± 0.03, no literals), a true four-way minimax objective `obj_fourway_minimax` was implemented and optimized directly across all four observable families, and nesting warm-starts and assertions were strictly enforced.

Here are the answers to the four questions posed in `AGENT-BRIEF-14.md`:

### 1.1 The True Four-Way Minimax per Scheme vs Suspended Values
| Scheme | Free Params ($k$) | Obs ($n_{\text{data}}$) | Res. DOF | Suspended Value (Brief 12/13) | **True Four-Way Minimax (Brief 14)** | **Omission Effect (Omission Contribution)** | Worst Observable Family |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `kinetic_jw95_free_rates` | 6 | 9 | +3 | 14.25σ | **14.25σ** | **0.00σ (0.0%)** | Deactivation ($z_{a_1} = -14.25\sigma, \sqrt{\chi^2_a} = 14.25\sigma$) |
| `kinetic_two_open` | 8 | 9 | +1 | 4.34σ | **3.10σ** | **1.24σ (28.6% of apparent gap)** | Deact vs CRC trade-off ($\sqrt{\chi^2_a} = 3.10\sigma, z_{n_H} = -3.10\sigma$) |
| `kinetic_coop` | 9 | 9 | **0** | 17.69σ | **2.39σ** | **15.30σ (86.5% of apparent gap)** | Balanced minimax across all 4 ($\le 2.39\sigma$) |

* **`kinetic_jw95_free_rates`:** 0.00σ of the 14.25σ was omission. The failure is completely dominated by structural mode cancellation on deactivation ($a_1 \le 0$). Recovery was not active at the bottleneck.
* **`kinetic_two_open`:** **1.24σ** of the 4.34σ (28.6%) was artificial inflation from omitting recovery during optimization. When the optimizer is explicitly asked to satisfy all four, the worst-case residual drops from 4.34σ to **3.10σ** ($1.03\times$ beyond $3\sigma$).
* **`kinetic_coop`:** **15.30σ** of the 17.69σ (86.5%) was pure optimization failure. Omitting recovery allowed the unconstrained optimizer to choose $c = 184$ with $k_{\text{off}} = 79.3\ \text{ms}^{-1}$, which blew out recovery to 17.69σ. When warm-started from `kinetic_two_open` and optimized under the true 4-way objective, the minimax residual drops to **2.39σ**.

### 1.2 Is there a vector satisfying all four within 3σ?
* **`kinetic_jw95_free_rates`:** **NO.** (Minimax 14.25σ, 4.75× beyond 3σ).
* **`kinetic_two_open`:** **NO.** (Minimax 3.10σ, 1.03× beyond 3σ; misses the 3σ feasibility threshold by only 0.10σ).
* **`kinetic_coop`:** **YES.** Under an objective that contains every observable, `kinetic_coop` achieves:
  * **Four-Way Minimax Optimum:** **2.39σ** (all 4 observable families $\le 2.39\sigma$, factor of **0.80×** of 3σ).
  * **Natural Joint Optimum:** **2.93σ** (all 4 observable families $\le 2.93\sigma$, factor of **0.98×** of 3σ).

### 1.3 If no: which pair is irreconcilable under an objective that includes both?
* **`kinetic_jw95_free_rates`:** Deactivation fast amplitude $a_1$ ($z_{a_1} = -14.25\sigma$) vs EC50 ($z_{\text{EC50}} = +14.23\sigma$). Single-open topology cannot produce positive $a_1$ while fitting the concentration-response curve.
* **`kinetic_two_open`:** **Deactivation vs CRC Slope ($n_H$).**
  Under the four-way objective, recovery is successfully brought inside 3σ ($z_{pp} = +2.81\sigma$) and $P_{o,max} = 0.790$ ($z_{po} = +0.81\sigma$). But the model is tightly pinched at **3.10σ** between deactivation ($\sqrt{\chi^2_a} = 3.10\sigma, z_{\tau_2} = -3.10\sigma$) and CRC ($z_{\text{EC50}} = -3.10\sigma, z_{n_H} = -3.10\sigma$).

### 1.4 Whether the nesting assertion held, and for which schemes
* **YES, strictly held and verified for `kinetic_coop` nested from `kinetic_two_open`.**
  * `obj_natural_joint`: coop $8.951122 \le$ two_open $15.006455$ (passed).
  * `obj_fourway_minimax`: coop $2.386095 \le$ two_open $3.099302$ (passed).
  Both objectives strictly obeyed monotonicity ($\text{coop} \le \text{two\_open} + 10^{-4}$), verified by automated assertion.

---

## 2. Status of Previous Findings & Brief Verdicts

1. **Brief 12 Verdict:**
   * **SUPERSEDED.** The reported 4.34σ was an upper bound resulting from omitting recovery from the objective. The true four-way minimax is **3.10σ**. Superseded header added to `RESULTS-BRIEF-12.md`.
2. **Brief 13 Verdict:**
   * **SUPERSEDED.** The reported 17.69σ was an optimizer failure caused by omitting recovery and omitting warm-start from `kinetic_two_open`. The true four-way minimax is **2.39σ**. Superseded header added to `RESULTS-BRIEF-13.md`.
3. **`FINDING-minimax-conflict.md`:**
   * **WITHDRAWN.** The headline claim that *"no scheme satisfies all four observables"* is withdrawn because `kinetic_coop` achieves 2.39σ (all $\le 3\sigma$). Furthermore, for `kinetic_two_open`, the gap shrank from 4.34σ to 3.10σ ($1.03\times$), demonstrating that most of the apparent conflict was scoring and optimization omission rather than an insurmountable physical barrier.

---

## 3. Parametric Accounting & Physical Interpretation

### 3.1 Degrees of Freedom & The Exact Saturation Caveat
The observable count across all four targets is:
* Deactivation: 3 taus + 2 independent areas (rank-2 simplex) = 5 terms
* Peak CRC: EC50 + Hill slope $n_H$ = 2 terms
* Paired-pulse recovery: 1 term ($0.33 \pm 0.03$)
* Peak open probability: 1 term ($0.75 \pm 0.05$)
* **Total $n_{\text{data}} = 9$**

Degrees of freedom per scheme:
* `kinetic_jw95_free_rates`: $k=6 \implies \text{res\_dof} = 9 - 6 = +3$. Over-determined, decisively rejected (14.25σ).
* `kinetic_two_open`: $k=8 \implies \text{res\_dof} = 9 - 8 = +1$. Operates with 1 residual degree of freedom; misses 3σ by 0.10σ (3.10σ).
* `kinetic_coop`: $k=9 \implies \text{res\_dof} = 9 - 9 = \mathbf{0}$.

> [!WARNING]
> **THE EXACT SATURATION CAVEAT ($n_{\text{data}} = 9, k_{\text{fitted}} = 9$, Residual DOF = 0):**  
> `kinetic_coop` fits 9 free parameters against exactly 9 data points. Consequently, its residual degrees of freedom is **0**.  
> The 2.39σ result is an **existence demonstration of algebraic flexibility**, NOT a statistical hypothesis test.

### 3.2 Physical Plausibility of Cooperativity Parameter $c$
At the Four-Way Minimax Optimum for `kinetic_coop`:
* Cooperativity factor: $c = 815.94$ ($\log_{10} c = 2.9117$).
* First unbinding: $k_{\text{off}} = 75.55\ \text{ms}^{-1}$ (lifetime $\tau_1 = 13.2\ \mu\text{s}$)
* Second unbinding: $k_{\text{off2}} = 0.0926\ \text{ms}^{-1}$ (lifetime $\tau_2 = 10.8\ \text{ms}$)
* Microscopic ratio: $c = k_{\text{off}} / k_{\text{off2}} = 815.94$.

While mathematically admissible and detailed-balance compliant (the state graph is an acyclic tree), an 800-fold difference in unbinding rates between the first and second agonist binding sites represents extreme positive cooperativity. Singly bound GABA dissociates in 13 microseconds, whereas doubly bound GABA is stabilized by three orders of magnitude.

---

## 4. Full Per-Dataset Residuals Summary Table (in $\sigma$)

| Scheme | Evaluation Profile | $z_{\tau_1}$ | $z_{\tau_2}$ | $z_{\tau_3}$ | $z_{a_1}$ | $z_{a_2}$ | $z_{a_3}$ | $\sqrt{\chi^2_a}$ | $z_{\text{EC50}}$ | $z_{n_H}$ | $z_{pp}$ (Recovery) | $P_{o,max}$ ($z_{po}$) | CRC RMSE Ratio (Rep.) | Max Viol / $3\sigma$ | Satisfies All 4 $\le 3\sigma$? |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `kinetic_jw95_free_rates` | Deact-Alone Optimum | -8.10$\sigma$ | +0.02$\sigma$ | +0.06$\sigma$ | -14.43$\sigma$ | +9.91$\sigma$ | +12.64$\sigma$ | 14.43$\sigma$ | -12.49$\sigma$ | -3.48$\sigma$ | +2.74$\sigma$ | 0.100 (-13.01$\sigma$) | 81.99$\times$ | **4.81** | **NO** |
| `kinetic_jw95_free_rates` | CRC-Alone Optimum | -9.33$\sigma$ | +37.05$\sigma$ | +50.99$\sigma$ | -14.25$\sigma$ | +38.10$\sigma$ | -6.40$\sigma$ | 41.50$\sigma$ | +0.08$\sigma$ | -0.47$\sigma$ | -0.56$\sigma$ | 0.998 (+4.96$\sigma$) | 2.58$\times$ | **17.00** | **NO** |
| `kinetic_jw95_free_rates` | Two-Observable Pareto Trade-off | -9.33$\sigma$ | -5.85$\sigma$ | +5.65$\sigma$ | -14.25$\sigma$ | +9.80$\sigma$ | +12.47$\sigma$ | 14.25$\sigma$ | -12.81$\sigma$ | +0.10$\sigma$ | -11.00$\sigma$ | 0.999 (+4.98$\sigma$) | 111.83$\times$ | **4.75** | **NO** |
| `kinetic_jw95_free_rates` | Natural Joint Optimum | -8.39$\sigma$ | -0.14$\sigma$ | +0.60$\sigma$ | -14.38$\sigma$ | +9.82$\sigma$ | +12.62$\sigma$ | 14.38$\sigma$ | +0.06$\sigma$ | -1.08$\sigma$ | +0.78$\sigma$ | 0.750 (-0.00$\sigma$) | 14.55$\times$ | **4.79** | **NO** |
| `kinetic_jw95_free_rates` | Four-Way Minimax Optimum | -9.33$\sigma$ | -5.73$\sigma$ | +14.25$\sigma$ | -14.25$\sigma$ | +9.80$\sigma$ | +12.47$\sigma$ | 14.25$\sigma$ | +14.23$\sigma$ | -0.75$\sigma$ | -6.17$\sigma$ | 0.999 (+4.98$\sigma$) | 33.11$\times$ | **4.75** | **NO** |
| `kinetic_two_open` | Deact-Alone Optimum | +0.00$\sigma$ | +0.00$\sigma$ | -0.00$\sigma$ | +0.12$\sigma$ | +0.25$\sigma$ | +0.17$\sigma$ | 0.00$\sigma$ | -2.66$\sigma$ | -3.80$\sigma$ | +3.91$\sigma$ | 0.858 (+2.16$\sigma$) | 24.23$\times$ | **1.30** | **NO** |
| `kinetic_two_open` | CRC-Alone Optimum | -9.33$\sigma$ | +131.53$\sigma$ | +4340.18$\sigma$ | -14.25$\sigma$ | -11.50$\sigma$ | +26.67$\sigma$ | 32.61$\sigma$ | -0.01$\sigma$ | -0.28$\sigma$ | -10.99$\sigma$ | 1.000 (+5.00$\sigma$) | 6.81$\times$ | **1446.73** | **NO** |
| `kinetic_two_open` | Two-Observable Pareto Trade-off | +0.70$\sigma$ | -3.11$\sigma$ | -3.09$\sigma$ | +2.86$\sigma$ | -0.48$\sigma$ | -2.91$\sigma$ | 3.11$\sigma$ | -3.11$\sigma$ | -3.11$\sigma$ | +6.58$\sigma$ | 0.756 (+0.12$\sigma$) | 19.49$\times$ | **2.19** | **NO** |
| `kinetic_two_open` | Natural Joint Optimum | +0.05$\sigma$ | -0.20$\sigma$ | -0.02$\sigma$ | +0.48$\sigma$ | -0.13$\sigma$ | -0.52$\sigma$ | 0.53$\sigma$ | -0.07$\sigma$ | -3.83$\sigma$ | -0.00$\sigma$ | 0.754 (+0.08$\sigma$) | 23.76$\times$ | **1.28** | **NO** |
| `kinetic_two_open` | Four-Way Minimax Optimum | -2.08$\sigma$ | -3.10$\sigma$ | +0.35$\sigma$ | +2.98$\sigma$ | -0.78$\sigma$ | -2.90$\sigma$ | 3.10$\sigma$ | -3.10$\sigma$ | -3.10$\sigma$ | +2.81$\sigma$ | 0.790 (+0.81$\sigma$) | 19.44$\times$ | **1.03** | **NO** |
| `kinetic_coop` | Deact-Alone Optimum | +0.00$\sigma$ | +0.00$\sigma$ | -0.00$\sigma$ | +0.12$\sigma$ | +0.25$\sigma$ | +0.17$\sigma$ | 0.00$\sigma$ | -2.66$\sigma$ | -3.80$\sigma$ | +3.91$\sigma$ | 0.858 (+2.16$\sigma$) | 24.23$\times$ | **1.30** | **NO** |
| `kinetic_coop` | CRC-Alone Optimum | -9.33$\sigma$ | +131.53$\sigma$ | +4340.17$\sigma$ | -14.25$\sigma$ | -11.50$\sigma$ | +26.67$\sigma$ | 32.61$\sigma$ | -0.01$\sigma$ | -0.28$\sigma$ | -10.99$\sigma$ | 1.000 (+5.00$\sigma$) | 6.81$\times$ | **1446.72** | **NO** |
| `kinetic_coop` | Two-Observable Pareto Trade-off | -0.57$\sigma$ | -3.10$\sigma$ | -2.43$\sigma$ | +3.09$\sigma$ | -1.01$\sigma$ | -2.86$\sigma$ | 3.10$\sigma$ | -3.10$\sigma$ | -3.10$\sigma$ | +5.48$\sigma$ | 0.770 (+0.40$\sigma$) | 19.33$\times$ | **1.83** | **NO** |
| `kinetic_coop` | Natural Joint Optimum | +0.43$\sigma$ | -0.07$\sigma$ | -0.31$\sigma$ | -0.02$\sigma$ | +0.20$\sigma$ | -0.03$\sigma$ | 0.19$\sigma$ | +0.18$\sigma$ | -2.93$\sigma$ | -0.07$\sigma$ | 0.749 (-0.01$\sigma$) | 12.01$\times$ | **0.98** | **YES** |
| `kinetic_coop` | Four-Way Minimax Optimum | +2.39$\sigma$ | -0.34$\sigma$ | -2.39$\sigma$ | +0.09$\sigma$ | +1.71$\sigma$ | -1.20$\sigma$ | 2.39$\sigma$ | +2.39$\sigma$ | -2.39$\sigma$ | +1.02$\sigma$ | 0.721 (-0.58$\sigma$) | 11.59$\times$ | **0.80** | **YES** |

---

## 5. Optimal Parameter Vectors & Pinned Bounds

| Scheme | Profile | Free Parameters Vector ($\log_{10}$) | Bounds Status / Pinned Parameters | $k_{fitted}$ | $k_{eff}$ |
| :--- | :--- | :--- | :--- | :---: | :---: |
| `kinetic_jw95_free_rates` | Two-Observable Pareto Trade-off | `log10_kd=-0.126, log10_E=3.000, log10_D=2.965, log10_r=-1.104, log10_koff=0.507, log10_alpha=2.000` | log10_E on upper bound (3.0); log10_alpha on upper bound (2.0) | 6 | 4 |
| `kinetic_jw95_free_rates` | Natural Joint Optimum | `log10_kd=1.799, log10_E=0.477, log10_D=2.137, log10_r=-1.718, log10_koff=-0.404, log10_alpha=-1.510` | None (all interior) | 6 | 6 |
| `kinetic_jw95_free_rates` | Four-Way Minimax Optimum | `log10_kd=2.425, log10_E=3.000, log10_D=2.963, log10_r=-1.138, log10_koff=0.352, log10_alpha=2.000` | log10_E on upper bound (3.0); log10_alpha on upper bound (2.0) | 6 | 4 |
| `kinetic_two_open` | Two-Observable Pareto Trade-off | `log10_kd=0.519, log10_E=0.221, log10_D=1.507, log10_r=-1.781, log10_koff=-0.819, log10_alpha=0.251, log10_delta=-1.295, log10_gamma=-1.230` | None (all interior) | 8 | 8 |
| `kinetic_two_open` | Natural Joint Optimum | `log10_kd=0.454, log10_E=-0.042, log10_D=1.423, log10_r=-1.832, log10_koff=-1.015, log10_alpha=1.531, log10_delta=-1.100, log10_gamma=-1.476` | None (all interior) | 8 | 8 |
| `kinetic_two_open` | Four-Way Minimax Optimum | `log10_kd=0.459, log10_E=0.257, log10_D=1.773, log10_r=-1.795, log10_koff=-0.718, log10_alpha=0.228, log10_delta=-1.200, log10_gamma=-1.235` | None (all interior) | 8 | 8 |
| `kinetic_coop` | Two-Observable Pareto Trade-off | `log10_kd=0.502, log10_E=0.238, log10_D=1.602, log10_r=-1.772, log10_koff=-0.767, log10_alpha=0.238, log10_delta=-1.262, log10_gamma=-1.234, log10_c=0.000` | None (all interior) | 9 | 9 |
| `kinetic_coop` | Natural Joint Optimum | `log10_kd=2.438, log10_E=-0.071, log10_D=1.381, log10_r=-1.825, log10_koff=1.937, log10_alpha=1.333, log10_delta=-1.082, log10_gamma=-1.483, log10_c=2.981` | None (all interior) | 9 | 9 |
| `kinetic_coop` | Four-Way Minimax Optimum | `log10_kd=2.530, log10_E=-0.073, log10_D=1.271, log10_r=-1.858, log10_koff=1.878, log10_alpha=1.413, log10_delta=-1.146, log10_gamma=-1.460, log10_c=2.912` | None (all interior) | 9 | 9 |

---

## 6. Summary of Task 3 Resolution (Anchor & Falsification Derivation)

Per Task 3 of `AGENT-BRIEF-14.md`, the manuscript and tests were updated to the sourced $11.6 \pm 0.9\ \mu\text{M}$ anchor:
1. **Manuscript (`knowledge/08-manuscript.md`):**
   - All tables, rates, and derived quantities updated from $20.0\ \mu\text{M} \to 11.6\ \mu\text{M}$ ($EC_{50}$ headroom $123.3\times \to 43.3\times$; contrast ratio $120\times \to 42\times$).
   - A prominent, dated note added at the top of Section §S3 preserving all previous legacy values and detailing the exact movement of every rate constant and asymptotic headroom factor.
2. **Falsification Threshold Derivation (`tests/test_review_regressions.py`):**
   - The threshold in `test_falsification_intervals_are_set_from_reachable_gain_not_the_asymptote` was derived directly from the physical divergence bounds between the exact reachable gain $G_{\text{reach}}(c) = (1 + \Delta_{\text{act}} P_o) / (1 + \Delta_{\text{aff}} \frac{c}{c + K_d})$ and asymptotic headroom $H_{\text{head}}(c)$, rather than widening a numerical tolerance.
   - All 16 manuscript consistency tests and 57 review regression tests pass.
