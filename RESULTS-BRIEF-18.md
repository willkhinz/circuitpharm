# Results Brief 18: Three Open States, Hill Slope Invariance, and Parameter Saturation

**Date:** 2026-10-18  
**Branch:** `claude/agent-brief-18`  
**Worktree:** `/Users/whinz/Biochemistry-wt-agent-brief-18`  
**Thread:** `brief18`  

---

## 1. Executive Summary & Plain Answers to Brief 18 Questions

Following the full-text reading of Keramidas & Harrison 2008 (PMC2213567), which reported that GABA-gated bursts of $\alpha1\beta2\gamma2S$ receptors require **three open components** to describe single-channel activation in both M-mode and H-mode, we developed, implemented, and evaluated 7-state cooperative models incorporating a third conducting open state (`KineticThreeOpenModel` in series and `KineticThreeOpenParallelModel` in parallel).

Here are the direct, plain answers to the four questions posed in Brief 18 §3:

### 1.1 What is the maximum achievable Hill slope for the three-open scheme?
* **At the Five-Way Minimax Optimum:** $n_H = 1.184$ ($z_{n_H} = -2.54\sigma$ against the published target $2.2 \pm 0.4$, Jahn et al. 1997). This is virtually identical to `kinetic_coop`'s $n_H = 1.178$ ($z_{n_H} = -2.55\sigma$).
* **Under Peak CRC Optimization Alone:** $n_H = 1.068$ ($z_{n_H} = -2.83\sigma$) for sequential and $1.064$ ($z_{n_H} = -2.84\sigma$) for parallel.
* **Comparison to Target & Ceiling:** The third open state **DOES NOT** relieve the $1.25$ ceiling. Both topologies remain strictly locked beneath $n_H \le 1.25$ and more than $2.5\sigma$ below the published $2.2 \pm 0.4$ target.
* **Mathematical Mechanism:** All open states ($A_2O, A_2O', A_2O''$) branch from the same doubly-bound shut state ($A_2R$) without binding or unbinding agonist. Consequently, all open states collapse into a single effective gating efficacy:
  $$E_{\text{eff}} = \frac{[O_{\text{total}}]}{[A_2R]} = E\left(1 + \frac{\delta}{\gamma} + \frac{\delta \delta_2}{\gamma \gamma_2}\right)$$
  The concentration-dependence of activation is governed **entirely by the binding steps** ($R \rightleftharpoons AR \rightleftharpoons A_2R$). Gating transitions among open states alter dwell times, but leave the shape and steepness of the concentration-response curve mathematically invariant.

### 1.2 Does the five-way minimax improve on `kinetic_coop`'s $2.58\sigma$?
* **Numerically:** The worst residual shifts marginally from $2.58\sigma$ to **$2.54\sigma$** (a negligible $0.04\sigma$ difference).
* **Parameter Behavior:** At the minimax optimum, $\log_{10} \delta_2$ pins hard to its lower bound ($-5.0$), where $\delta_2 \to 0$ collapses the third open state back into `kinetic_coop`. The optimizer explicitly shuts off the third open state.

### 1.3 Watch the degrees of freedom: Is this a tested fit?
* **NO. `kinetic_three_open` CEASES TO BE A TESTED HYPOTHESIS FIT.**
* **Degrees of Freedom Accounting (Reported Before Minimax):**
  - Total independent data constraints: $n_{\text{data}} = 10$.
  - Free fitted parameters: $k_{\text{fitted}} = 11$.
  - **Residual Degrees of Freedom:** $\mathbf{\text{res\_dof} = 10 - 11 = -1 \le 0}$.
* **Verdict:** Because $\text{res\_dof} < 0$, the model is overparameterized. Any residual improvement is an algebraic property of parameter saturation, not empirical evidence in favor of the scheme. `kinetic_coop` ($k=9, \text{res\_dof} = +1$) remains the sole surviving tested model.

### 1.4 Does the cooperativity requirement change from $[35, 2500]$ toward 5–23$\times$?
* **NO.** At the five-way minimax optimum, $c = 999.7$ ($\log_{10} c = 3.000$, pinned to the search box upper bound), with $k_{\text{off}} = 85.53\ \text{ms}^{-1}$ and $k_{\text{off2}} = 0.0856\ \text{ms}^{-1}$.
* Adding a third open state provides zero relief to the cooperativity requirement, because the demand for $c \gg 1$ arises from the tension between fast macroscopic deactivation ($\tau_1 = 2.8$ ms) and the slow dissociation of dually-bound receptor required to maintain nanomolar/micromolar sensitivity without premature deactivation.

---

## 2. Model Architecture, Connectivity Justification & Verification (Task 1)

### 2.1 Sequential vs Parallel Topology: Decision & Justification
Keramidas & Harrison (2008) established the existence of three open dwell-time components in single-channel burst records, but did not specify the microscopic transition topology. We considered two primary configurations:

1. **Sequential / Series Topology (`KineticThreeOpenModel`, Default):**
   $$R \underset{k_{\text{off}}}{\overset{2 k_{\text{on}} x}{\rightleftharpoons}} AR \underset{2 k_{\text{off2}}}{\overset{k_{\text{on}} x}{\rightleftharpoons}} A_2R \underset{\alpha}{\overset{\beta}{\rightleftharpoons}} A_2O \underset{\gamma}{\overset{\delta}{\rightleftharpoons}} A_2O' \underset{\gamma_2}{\overset{\delta_2}{\rightleftharpoons}} A_2O''$$
   with desensitization branching from $A_2R$: $A_2R \underset{r}{\overset{d}{\rightleftharpoons}} A_2D$.
   - **Justification:** Sequential opening represents progressive conformational relaxation or dilation of the pore domain within the open manifold, extending the natural two-open cascade ($A_2R \rightleftharpoons A_2O \rightleftharpoons A_2O'$).
   - **Analytical Properties:**
     - Maximum open probability:
       $$P_{o,\text{max}} = \frac{E(1 + K_2 + K_2 K_3)}{1 + E(1 + K_2 + K_2 K_3)}$$
       where $E = \beta/\alpha$, $K_2 = \delta/\gamma$, and $K_3 = \delta_2/\gamma_2$.
     - Mean open duration:
       $$\langle \tau_O \rangle = \frac{1}{\alpha} (1 + K_2 + K_2 K_3)$$
     - Direct nested reduction to `kinetic_coop` as $\delta_2 \to 0$.

2. **Parallel Topology (`KineticThreeOpenParallelModel`, Alternative):**
   $$A_2O' \underset{\delta}{\overset{\gamma}{\rightleftharpoons}} A_2O \underset{\alpha}{\overset{\beta}{\rightleftharpoons}} A_2R \underset{\alpha_p}{\overset{\beta_p}{\rightleftharpoons}} A_2O''$$
   - **Alternative Prediction:** In parallel topology, the receptor can transition directly between the resting shut state $A_2R$ and $A_2O''$ with independent efficacy $E_p = \beta_p / \alpha_p$. Exits from $A_2O''$ close directly to $A_2R$ without passing through $A_2O$, producing a split entry flux $\pi_{\mathcal{O}} = \left(\frac{\beta}{\beta + \beta_p}, 0, \frac{\beta_p}{\beta + \beta_p}\right)$.
   - Both topologies were implemented, unit-tested, and evaluated. As demonstrated in Section 4, both topologies exhibit identical Hill slope ceilings.

### 2.2 Microscopic Reversibility & Explicit Zero Assertions
* **Acyclic Tree Structure:** Both schemes contain $V = 7$ states and $E = 6$ bidirectional transitions. The fundamental cycle rank is:
  $$\text{Cycle Rank} = E - V + 1 = 6 - 7 + 1 = 0$$
  Because there are **no closed loops or cycles**, microscopic reversibility (detailed balance) is identically satisfied for all parameter sets without requiring cycle constraints ($P_{\text{cycle}} = 1.0$).
* **Explicit Zero Connectivity:** Transitions between non-adjacent states are strictly zero:
  - Series: $Q(A_2R, A_2O') = 0$, $Q(A_2R, A_2O'') = 0$, $Q(A_2O, A_2O'') = 0$, $Q(A_2D, \mathcal{O}) = 0$.
  - Parallel: $Q(A_2O, A_2O'') = 0$, $Q(A_2O', A_2O'') = 0$, $Q(A_2D, \mathcal{O}) = 0$.
  - These conditions are strictly verified by unit tests in `tests/test_three_open.py`.

### 2.3 Nesting & Likelihood Invariance
* In `comparison.py`, `kinetic_three_open` declares `nests=("kinetic_coop",)`.
* In `_model_three_open`, when $\delta_2 \le 10^{-5}\ \text{ms}^{-1}$, $\delta_2$ is clamped to $0.0$, collapsing $A_2O''$.
* Verified in `tests/test_three_open.py`: Embedding a `kinetic_coop` parameter vector via `nesting_embedding` reproduces the exact log-likelihood of `kinetic_coop` to numerical precision (absolute difference $< 10^{-12}$).
* Armed with an optimizer refusal assertion: `assert three_open <= coop + 1e-4`, which passed unconditionally.

---

## 3. Full-Text Literature Update & Critical Caveats (Task 3)

### 3.1 Keramidas & Harrison (2008) Full-Text Read
In `circuitpharm.literature`, `SOURCES["keramidas2008"]` was updated following full-text examination (PMC2213567):
* `where="full_text_pmc"`
* Recorded verbatim findings:
  - M-Mode bursts require 3 shut and 3 open components (10 patches).
  - H-Mode bursts require 2 shut and 3 open states (5 patches).
  - Component dwell times and fractional areas documented in `OpenTimeFit`.
  - Verified absence of stated error statistics: the $\pm$ values have no stated error type in the paper. The flat band constraint $[2.54, 2.83]$ ms remains the sole source-mandated treatment.

### 3.2 Burst Dwell-Time Components vs Macroscopic Scheme States (MANDATORY CAVEAT)
As required by Brief 18 §4, we explicitly state:
> **Assumption Warning:** The three open components reported by Keramidas & Harrison (2008) are exponential time constants extracted from dwell-time histograms of single-channel openings *within bursts* under continuous, saturating GABA (10 mM). Our kinetic scheme states ($A_2O, A_2O', A_2O''$) are states within a whole-receptor Markov model intended to reproduce non-equilibrium macroscopic synaptic currents (deactivation, recovery).
>
> **These are not automatically identical objects.** While a Markov scheme with $N$ open states yields an open dwell-time distribution with at most $N$ exponential components, the mapping assumes that intra-burst gating is governed by the same stationary Markov generator that drives fast macroscopic activation. We treat this mapping as an explicit modeling assumption, not an established biological identity.

### 3.3 Modal Gating is Completely Unrepresented
Keramidas & Harrison observed two distinct gating regimes:
* **M-mode (Medium $P_o \approx 0.7$):** 3 shut components, 3 open components.
* **H-mode (High $P_o \approx 0.9$):** 2 shut components, 3 open components.
Channels switch between M-mode and H-mode on a timescale of **seconds to tens of seconds**.
**Nothing in this project represents modal gating.** Neither `kinetic_two_open`, `kinetic_coop`, nor `kinetic_three_open` accounts for modal switching. All models assume a single stationary gating regime corresponding to M-mode (the predominant mode in 11/15 patches).

---

## 4. Symmetric Five-Way Fit Evaluation (Task 2)

### 4.1 Degrees of Freedom & Fit Status
| Scheme | Model Parameters ($k_{\text{fitted}}$) | Data Points ($n_{\text{data}}$) | Residual D.O.F. | Statistical Fit Status |
| :--- | :---: | :---: | :---: | :--- |
| `kinetic_two_open` | 8 | 10 | **+2** | Genuine Over-Determined Test (FAILS at $3.23\sigma$) |
| `kinetic_coop` | 9 | 10 | **+1** | Genuine Over-Determined Test (PASSES at $2.58\sigma$) |
| `kinetic_three_open` | 11 | 10 | **-1** | **OVERPARAMETERIZED / SATURATION (res_dof $\le 0$)** |

> [!IMPORTANT]
> Because `kinetic_three_open` operates with $\text{res\_dof} = -1 < 0$, it has more tunable parameters than independent data constraints. Its fit is not an empirical test of a hypothesis; it is an overparameterized system.

### 4.2 Per-Dataset Residuals Summary Table (in $\sigma$)

| Scheme | Evaluation Profile | $z_{\tau_1}$ | $z_{\tau_2}$ | $z_{\tau_3}$ | $z_{a_1}$ | $z_{a_2}$ | $z_{a_3}$ | $\sqrt{\chi^2_a}$ | $z_{\text{EC50}}$ | $z_{n_H}$ | $z_{pp}$ | $P_{o,\text{max}}$ | $\tau_O$ (ms) | $z_{\text{MOT}}$ (band) | Max Viol / $3\sigma$ | Satisfies All 5 $\le 3\sigma$? |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `kinetic_two_open` | Deact-Alone | -0.00$\sigma$ | +0.00$\sigma$ | +0.00$\sigma$ | +0.12$\sigma$ | +0.24$\sigma$ | +0.16$\sigma$ | 0.00$\sigma$ | -5.85$\sigma$ | -3.65$\sigma$ | +3.71$\sigma$ | 0.880 (+2.61$\sigma$) | 2.40 | -0.98$\sigma$ | **1.95** | **NO** |
| `kinetic_two_open` | Natural Joint | -0.42$\sigma$ | -1.23$\sigma$ | +1.51$\sigma$ | +0.50$\sigma$ | +0.24$\sigma$ | +0.22$\sigma$ | 0.30$\sigma$ | -0.10$\sigma$ | -3.89$\sigma$ | +2.46$\sigma$ | 0.770 (+0.41$\sigma$) | 2.84 | +0.04$\sigma$ | **1.30** | **NO** |
| `kinetic_two_open` | Five-Way Minimax | -3.18$\sigma$ | +0.06$\sigma$ | +3.23$\sigma$ | +3.27$\sigma$ | -2.46$\sigma$ | -2.10$\sigma$ | 3.23$\sigma$ | -3.23$\sigma$ | -3.23$\sigma$ | +0.98$\sigma$ | 0.862 (+2.25$\sigma$) | 2.22 | -2.22$\sigma$ | **1.08** | **NO** |
| `kinetic_coop` | Deact-Alone | +0.00$\sigma$ | +0.00$\sigma$ | +0.00$\sigma$ | +0.12$\sigma$ | +0.24$\sigma$ | +0.16$\sigma$ | 0.00$\sigma$ | -5.85$\sigma$ | -3.65$\sigma$ | +3.71$\sigma$ | 0.880 (+2.61$\sigma$) | 2.40 | -0.98$\sigma$ | **1.95** | **NO** |
| `kinetic_coop` | Natural Joint | -0.20$\sigma$ | -0.93$\sigma$ | +1.44$\sigma$ | +0.09$\sigma$ | +0.50$\sigma$ | +0.67$\sigma$ | 0.33$\sigma$ | +0.16$\sigma$ | -2.95$\sigma$ | +2.46$\sigma$ | 0.756 (+0.12$\sigma$) | 2.84 | +0.05$\sigma$ | **0.98** | **YES** |
| `kinetic_coop` | Five-Way Minimax | +0.50$\sigma$ | -0.24$\sigma$ | +2.50$\sigma$ | -0.96$\sigma$ | +3.08$\sigma$ | +0.40$\sigma$ | 2.58$\sigma$ | +2.24$\sigma$ | -2.55$\sigma$ | +2.58$\sigma$ | 0.775 (+0.50$\sigma$) | 3.20 | +2.58$\sigma$ | **0.86** | **YES** |
| `kinetic_three_open` | Deact-Alone | +0.00$\sigma$ | +0.00$\sigma$ | -0.00$\sigma$ | +0.12$\sigma$ | +0.24$\sigma$ | +0.16$\sigma$ | 0.00$\sigma$ | -5.85$\sigma$ | -3.65$\sigma$ | +3.71$\sigma$ | 0.880 (+2.61$\sigma$) | 2.40 | -0.98$\sigma$ | **1.95** | **NO** |
| `kinetic_three_open` | Natural Joint | -0.20$\sigma$ | -0.93$\sigma$ | +1.44$\sigma$ | +0.08$\sigma$ | +0.51$\sigma$ | +0.67$\sigma$ | 0.33$\sigma$ | +0.15$\sigma$ | -2.95$\sigma$ | +2.46$\sigma$ | 0.756 (+0.12$\sigma$) | 2.84 | +0.06$\sigma$ | **0.98** | **YES** |
| `kinetic_three_open` | Five-Way Minimax | +0.47$\sigma$ | -0.26$\sigma$ | +2.54$\sigma$ | -0.95$\sigma$ | +3.04$\sigma$ | +0.40$\sigma$ | 2.54$\sigma$ | +2.54$\sigma$ | -2.54$\sigma$ | +2.54$\sigma$ | 0.775 (+0.50$\sigma$) | 3.18 | +2.39$\sigma$ | **0.85** | **YES** |

### 4.3 Optimal Parameters & Bound Status
* **`kinetic_two_open` Minimax Optimum ($k_{\text{eff}} = 8$):**
  `log10_kd=0.523, log10_E=0.292, log10_D=1.999, log10_r=-1.882, log10_koff=-0.564, log10_alpha=0.159, log10_delta=-1.142, log10_gamma=-1.484` (All interior).
* **`kinetic_coop` Minimax Optimum ($k_{\text{eff}} = 8$):**
  `log10_kd=2.603, log10_E=-0.036, log10_D=1.458, log10_r=-1.851, log10_koff=1.933, log10_alpha=0.068, log10_delta=-0.989, log10_gamma=-1.428, log10_c=3.000` (`log10_c` pinned on upper bound 3.0).
* **`kinetic_three_open` Minimax Optimum ($k_{\text{eff}} = 9$):**
  `log10_kd=2.610, log10_E=-0.036, log10_D=1.459, log10_r=-1.849, log10_koff=1.932, log10_alpha=0.071, log10_delta=-0.989, log10_gamma=-1.427, log10_c=3.000, log10_delta2=-5.000, log10_gamma2=0.000` (`log10_c` on upper bound 3.0; `log10_delta2` pinned to lower bound -5.0).

---

## 5. Maximum Achievable Hill Slope Analysis (Task 5)

| Scheme | Architecture | Max $n_H$ ($EC_{50} \in [5, 25]\ \mu$M) | Sourced Peak CRC Optimum ($n_H$) | Five-Way Minimax Optimum ($n_H$) | Residual vs $2.2 \pm 0.4$ Target |
| :--- | :--- | :---: | :---: | :---: | :---: |
| `kinetic_two_open` | 2 Open (Series) | 0.976 | 0.976 | 0.907 | $-3.23\sigma$ |
| `kinetic_coop` | 2 Open (Series) + Coop | 1.074 | 1.074 | 1.178 | $-2.55\sigma$ |
| `kinetic_three_open` | 3 Open (Series) + Coop | 1.068 | 1.068 | 1.184 | $-2.54\sigma$ |
| `kinetic_three_open_parallel` | 3 Open (Parallel) + Coop | 1.064 | 1.064 | 1.175 | $-2.56\sigma$ |

### Mathematical Proof of the Hill-Slope Invariance
Why does adding open states fail to raise the Hill slope?
1. **Binding Polynomial Factorization:**
   In any scheme where all open states branch exclusively from the doubly-liganded state $A_2R$:
   $$[O_{\text{total}}] = \sum_{i} [O_i] = E_{\text{eff}} [A_2R]$$
   For an agonist concentration $x$, the fractional open response normalized to peak saturation is:
   $$Y(x) = \frac{[O_{\text{total}}](x)}{[O_{\text{total}}]_{\text{max}}} = \frac{[A_2R](x)}{[A_2R]_{\text{max}}}$$
   Notice that $E_{\text{eff}}$ cancels out identically in the normalized concentration-response relation!
2. **Independence from Open State Count:**
   The normalized response $Y(x)$ is governed **exclusively by the binding transitions** that generate $A_2R$:
   $$Y(x) = \frac{\frac{x^2}{K_1 K_2}}{1 + \frac{x}{K_1} + \frac{x^2}{K_1 K_2} (1 + D + E_{\text{eff}})}$$
3. **The Independent-Binding Ceiling ($4 - 2\sqrt{2} \approx 1.1716$):**
   Without cooperativity ($c=1$), an independent two-site binding receptor has an analytical maximum logit slope at $EC_{50}$ of:
   $$n_H(EC_{50}) = 4 - 2\sqrt{2} \approx 1.1716$$
   While positive cooperativity ($c \to \infty$) can increase the theoretical asymptote to $2.0$, when the scheme is jointly constrained by macroscopic deactivation ($\tau_1 = 2.8$ ms), $c \approx 1000$ yields an operational Hill slope over the measured $1\text{--}10\ \mu$M window of **$1.18$**.
4. **Conclusion:** Adding a third open state (series or parallel) redistributes lifetimes within the open manifold, but has **zero mathematical effect on the agonist-binding cooperativity**. It cannot raise the Hill slope toward 2.2.

---

## 6. Verification & Test Suite Execution (Task 5)

Both tiers of pytest were executed in the branch worktree:

### 6.1 Fast Test Suite (`not slow`)
Command:
```bash
PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python -m pytest -q -m "not slow"
```
Verbatim final summary line:
```
594 passed, 41 warnings in 154.03s (0:02:34)
```
*(Baseline on main: 588 passed; net +6 tests added and passing).*

### 6.2 Full Test Suite (All Tests)
Command:
```bash
PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python -m pytest -q
```
Verbatim final summary line:
```
668 passed, 116 warnings in 938.30s (0:15:38)
```
*(Baseline on main: 662 passed; net +6 tests added and passing).*

---

## 7. Conclusions & Recommendations for Manuscript

1. **Structural Corroboration Confirmed:** Keramidas & Harrison 2008 provides strong external biological corroboration that $\alpha1\beta2\gamma2S$ receptors possess multiple open states.
2. **Model Selection Remains `kinetic_coop`:** Adding a third open state pushes residual degrees of freedom to $-1 \le 0$ (saturation) without improving the fit ($2.54\sigma$ vs $2.58\sigma$, with $\delta_2 \to 0$). `kinetic_coop` ($k=9, \text{res\_dof} = +1$) remains the sole parsimonious, tested model.
3. **Hill Slope Ceiling is Fundamental to Gating Architecture:** The gap between the model's $n_H \approx 1.18$ and the published $2.2 \pm 0.4$ cannot be resolved by adding open states. If the published Hill slope is biological truth, it requires either higher-order binding cooperativity ($>2$ binding steps), subunit-level concerted allostery (MWC model), or state-dependent agonist trapping, rather than additional post-binding open states.
