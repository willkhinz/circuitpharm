# Results Brief 13 — Cooperative Binding and EC50 Target Resolution

**Worktree:** `/Users/whinz/Biochemistry-wt-agent-brief-13`  
**Branch:** `claude/agent-brief-13`  
**Base Commit on `main`:** `918cea9`  
**Task 1 Commit:** `6e8d877` (`Brief 13 Task 1: implement kinetic_coop model with binding cooperativity c = koff/koff2 and four-way joint fit`)  
**Task 2 Commit:** `4ebb015` (`Brief 13 Task 2: resolve FIT_TARGETS['ec50_um'] from legacy unsourced 20.0 uM to Jahn 1997 sourced 11.6 uM with provenance record`)  
**Python Interpreter:** `/Users/whinz/Biochemistry/.venv/bin/python`  
**PYTHONPATH:** `/Users/whinz/Biochemistry-wt-agent-brief-13/src`  

---

## 1. Executive Summary & Key Questions Answered

### What is the minimax residual for `kinetic_coop`? Compare against 4.34σ and 14.25σ.
* **`kinetic_coop` (9 params, 9 data points, dof = 0):**
  * **Minimax Residual:** **17.69σ** on the Pairwise Pareto Frontier; **19.32σ** at the Natural Joint Optimum.
  * **Comparison:**
    * `kinetic_two_open` (8 params, 9 data points, dof = 1): **4.34σ** (limited by recovery +4.34σ and Hill slope −3.85σ).
    * `kinetic_jw95_free_rates` (6 params, 9 data points, dof = 3): **14.25σ** (limited by deactivation amplitude −14.4σ).
    * **`kinetic_coop` performs WORSE on minimax than `kinetic_two_open` (17.69σ vs 4.34σ)** despite having an additional free parameter.

### Does it close the 3σ gap?
* **NO.** The minimax residual of **17.69σ** exceeds 3σ by **5.90×**.

### What cooperativity factor `c` does it need, and is that value physically sensible?
* **Cooperativity parameter:** $c = 184.22$ ($c > 1$ denotes positive cooperativity, where the second binding step dissociates slower / binds tighter than the first).
* **Explicit dissociation rates:**
  * First unbinding: $k_{\text{off}} = 79.29\ \text{ms}^{-1}$
  * Second unbinding: $k_{\text{off2}} = 0.4304\ \text{ms}^{-1}$
  * Microscopic ratio: $c = k_{\text{off}} / k_{\text{off2}} = 184.22$.
* **Physical plausibility:** **Physically implausible.**
  * $k_{\text{off}} = 79.29\ \text{ms}^{-1}$ corresponds to an unbinding lifetime $\tau = 1 / k_{\text{off}} \approx 12.6\ \mu\text{s}$.
  * Singly-bound GABA dissociates almost instantaneously. While this uncoupling allows the model to produce a steep Hill slope ($n_H = 1.99$, residual $z_{n_H} = -0.52\sigma$) without slowing deactivation, it causes catastrophic failure on paired-pulse recovery ($z_{\text{rec}} = +17.69\sigma$ to $+19.32\sigma$). During the 500 ms inter-pulse interval, the receptor completely sheds agonist and resensitizes prematurely, yielding 91% recovery against the published $0.33 \pm 0.03$ (33%).

### Exact Saturation Caveat (Prominent Warning)
> [!WARNING]
> **EXACT SATURATION CAVEAT ($n_{\text{data}} = 9, k_{\text{fitted}} = 9$, Residual DOF = 0):**  
> `kinetic_coop` fits 9 free parameters against 9 data points ($n_{\text{data}} = 9$). Consequently, the residual degrees of freedom is exactly **0**.  
> The fit is **exactly saturated** and is **not a statistical test at all**. Fitting these data is at best an existence demonstration of algebraic flexibility. Even with zero residual degrees of freedom, the model *still* cannot reconcile paired-pulse recovery with the macroscopic targets, blowing out recovery to $+17.69\sigma$. The honest comparison is against `kinetic_two_open`'s 8 parameters at the same budget (minimax 4.34σ).

### Parameter Bound Flags and Effective Parameter Count
* All 9 parameters of `kinetic_coop` landed strictly interior to their search bounds.
* Active bound flags: **None**.
* $k_{\text{effective}} = 9$.

---

## 2. Microscopic Reversibility & Nesting Verification

### Microscopic Reversibility
* **Topology:**
  The state graph of `kinetic_coop` consists of 6 states:
  $$\text{R} \xrightleftharpoons[2 k_{\text{off}}]{\ 2 k_{\text{on}} [\text{GABA}]\ } \text{AR} \xrightleftharpoons[2 k_{\text{off2}}]{\ k_{\text{on}} [\text{GABA}]\ } \text{A}_2\text{R} \xrightleftharpoons[\alpha]{\beta} \text{A}_2\text{O} \xrightleftharpoons[\gamma]{\delta} \text{A}_2\text{O}'$$
  with a desensitized branch at $\text{A}_2\text{R}$:
  $$\text{A}_2\text{R} \xrightleftharpoons[r]{d} \text{A}_2\text{D}$$
* **Proof:**
  The transition network is a **tree graph** (it contains no undirected cycles). By Kolmogorov's criterion, microscopic reversibility (detailed balance) holds unconditionally for any acyclic graph without requiring cycle-closure constraints. Stationary probabilities satisfy detailed balance $p_i q_{ij} = p_j q_{ji}$ at every edge for all agonist concentrations. Tested and verified across 7 orders of magnitude in `tests/test_nextgen_models.py::test_kinetic_coop_connectivity_and_zeros`.

### Nesting Topology & Numerical Embedding
* **Exact Reduction:**
  When $c = 1.0$ ($k_{\text{off2}} = k_{\text{off}}$), `kinetic_coop` reduces identically to `kinetic_two_open`:
  * Transition generator: $\|Q_{\text{coop}}(g) - Q_{\text{two\_open}}(g)\|_\infty < 10^{-14}$.
  * Dose-response: $\|\text{DR}_{\text{coop}} - \text{DR}_{\text{two\_open}}\|_\infty < 10^{-14}$.
  * Gating maximum: $|P_{o,\max}^{\text{coop}} - P_{o,\max}^{\text{two\_open}}| < 10^{-14}$.
* **Log-Likelihood Identity:**
  Using `nesting_embedding(MODEL_SPECS["kinetic_two_open"], MODEL_SPECS["kinetic_coop"], v_sub)`, total log-likelihood over all three benchmark datasets matches to numerical tolerance:
  $$\ln L_{\text{coop}}(\text{embedded}) = \ln L_{\text{two\_open}} \pm 10^{-5}$$
  Verified in `tests/test_nesting.py::test_kinetic_coop_embedding_reproduces_kinetic_two_open_likelihood` and `tests/test_nesting.py::test_nesting_on_total_objective`.

---

## 3. Four-Way Joint Fit Comparison Table

Evaluated on `scripts/joint_fourway_fit.py` with identical LHS + Nelder-Mead budget (2,000 LHS samples + top 12 polished starts):

| Observable / Metric | Target $\pm$ SEM | `kinetic_coop` (9 params) | `kinetic_two_open` (8 params) | `kinetic_jw95_free_rates` (6 params) |
|---|---|---|---|---|
| **Free Parameters ($k_{\text{fitted}}$)** | — | **9** | **8** | **6** |
| **Data Points ($n_{\text{data}}$)** | — | **9** | **9** | **9** |
| **Residual DOF ($n - k$)** | — | **0 (EXACTLY SATURATED)** | **1** | **3** |
| **Minimax Residual (Pareto)** | $\le 3.00\sigma$ | **17.69σ (FAIL)** | **4.34σ (FAIL)** | **14.25σ (FAIL)** |
| **Minimax Residual (Natural Joint)** | $\le 3.00\sigma$ | **19.32σ (FAIL)** | **4.45σ (FAIL)** | **14.25σ (FAIL)** |
| **$P_{o,\max}$ Anchor** | $0.75 \pm 0.05$ | $0.750\ (+0.01\sigma)$ | $0.755\ (+0.10\sigma)$ | $0.748\ (-0.04\sigma)$ |
| **Peak EC50 ($\mu\text{M}$)** | $11.6 \pm 0.9$ | $11.60\ (+0.00\sigma)$ | $11.58\ (-0.02\sigma)$ | $11.60\ (+0.00\sigma)$ |
| **Hill Slope $n_H$ ($1\text{--}10\ \mu\text{M}$)** | $2.20 \pm 0.40$ | **1.99 (−0.52σ)** | **0.66 (−3.85σ)** | **1.33 (−2.18σ)** |
| **Deact Fast Tau $\tau_{\text{fast}}$ (ms)** | $8.0 \pm 1.2$ | $8.29\ (+0.24\sigma)$ | $8.60\ (+0.50\sigma)$ | $8.00\ (+0.00\sigma)$ |
| **Deact Slow Tau $\tau_{\text{slow}}$ (ms)** | $77.8 \pm 15.6$ | $77.74\ (-0.00\sigma)$ | $77.80\ (+0.00\sigma)$ | $77.80\ (+0.00\sigma)$ |
| **Deact Amp $A_{\text{fast}}$** | $0.80 \pm 0.04$ | $0.791\ (-0.23\sigma)$ | $0.793\ (-0.18\sigma)$ | $0.224\ (-14.40\sigma)$ |
| **Paired-Pulse Recovery (500 ms)** | $0.33 \pm 0.03$ | **0.86 (+17.69σ)** | **0.46 (+4.34σ)** | **0.33 (+0.00σ)** |
| **Fitted $c = k_{\text{off}} / k_{\text{off2}}$** | — | **184.22** | 1.00 (fixed) | 1.00 (fixed) |
| **Fitted $k_{\text{off}}\ (\text{ms}^{-1})$** | — | **79.29** | 0.414 | 0.050 |
| **Fitted $k_{\text{off2}}\ (\text{ms}^{-1})$** | — | **0.430** | 0.414 | 0.050 |
| **Bound Flags** | — | None (all interior) | None (all interior) | None (all interior) |
| **$k_{\text{effective}}$** | — | **9** | **8** | **6** |

### Trade-off Analysis: Why Cooperativity Fails
Adding binding cooperativity ($c = 184.2$) directly solves the Hill slope problem (raising $n_H$ from 0.66 to 1.99, reducing residual from $-3.85\sigma$ to $-0.52\sigma$). However, to achieve this while maintaining the fast deactivation component ($\tau_{\text{fast}} = 8.29\ \text{ms}$), the optimizer drives the unbinding of the first site to an unphysiologically rapid rate ($k_{\text{off}} \approx 79.3\ \text{ms}^{-1}$). At 500 ms between pulses, receptors strip bound agonist almost instantly and return to the resting state $\text{R}$, destroying paired-pulse depression and exploding the recovery residual to $+17.69\sigma$.

---

## 4. Task 2: Resolution of Unsourced `FIT_TARGETS['ec50_um']`

### Provenance Attribution and Code Changes
1. **Target Update (`src/circuitpharm/gabaa_kinetics.py`):**
   * Legacy value: `20.0` $\mu\text{M}$ (*"recalled literature middle from earlier legacy tuning files"*, UNSOURCED, $+9.33\sigma$ off).
   * Resolved value: `FIT_TARGETS['ec50_um'] = 11.6` $\mu\text{M}$ with `FIT_TARGET_SEMS['ec50_um'] = 0.9` $\mu\text{M}$.
   * Source attribution: `jahn_a1b2g2_kinetics` (Jahn et al. 1997, *NeuroReport* 8(16):3443-6) for $\alpha_1\beta_2\gamma_{2\text{L}}$ in HEK293 measured with ultra-fast solution exchange.
   * Old value preserved in comments with historical context.
2. **Provenance Registry (`src/circuitpharm/provenance.py`):**
   * Added `("gabaa_kinetics", "ec50_um"): Record(Basis.QUANTITATIVE, "jahn_a1b2g2_kinetics", ...)` in `MODEL_PARAM_PROV`.
   * Verified in code that `Basis.QUANTITATIVE` is the exact enum variant used for experimentally sourced quantitative numbers (alongside `jahn1997_peak_crc`).
   * Updated `DATASET_PROV["jahn1997_peak_crc"]` note to record the resolution of the discrepancy.
   * Provenance tally moves from 57 to 58 total registered parameters; sourced parameters move from 7 to 8; `QUANTITATIVE` parameters move from 1 to 2.
3. **Discrepancy Surface (`src/circuitpharm/fitting/crc_components.py`):**
   * Updated `surface_ec50_discrepancy()`: discrepancy reduced from $8.4\ \mu\text{M}$ ($+9.33\sigma$) to $0.0\ \mu\text{M}$ ($0.0\sigma$); status updated to resolved.

### Test Impact Audit: Old vs New Figures

The canonical square fit (`gabaa_kinetics.fit_scheme()` with `fixed_alpha = 0.30`) changes from fitting $EC_{50} = 20.0\ \mu\text{M}$ to $11.6\ \mu\text{M}$. Below is the complete record of every test whose expectations move:

| File / Test | Parameter / Observable | Old Expected Value (at 20.0 µM) | New Sourced Value (at 11.6 µM) |
|---|---|---|---|
| `tests/test_fit_determinism.py::test_the_fit_is_the_value_this_repository_pins` | $k_{\text{on}}$ ($\mu\text{M}^{-1}\text{ms}^{-1}$) | `0.0075785` | `0.0130384` |
| | $k_{\text{off}}$ ($\text{ms}^{-1}$) | `0.180415` | `0.180304` |
| | $\beta$ ($\text{ms}^{-1}$) | `1.177468` | `1.175664` |
| | $K_d$ ($\mu\text{M}$) | `23.8062` | `13.8287` |
| | Gating Efficacy $E = \beta / \alpha$ | `3.9249` | `3.9189` |
| | Asymptotic Headroom ($0.4\ \mu\text{M}$) | `123.344×` | `43.259×` |
| `tests/test_p4_likelihood.py::test_registry_fitted_kinetic_carries_its_caveats` | $K_d = k_{\text{off}} / k_{\text{on}}$ | `23.8062` | `13.8287` |
| `tests/test_provenance.py::test_the_provenance_gap_is_what_the_audit_recorded` | Total parameters audited | `57` | `58` |
| | Sourced parameters | `7` | `8` |
| | `Basis.QUANTITATIVE` count | `1` | `2` |
| `tests/test_crc_components.py::test_ec50_discrepancy_surfaced_and_unresolved` | `fit_targets_ec50_um` | `20.0` | `11.6` |
| | `discrepancy_um` | `8.4` | `0.0` |
| | `discrepancy_sigma` | `9.333` | `0.0` |
| | `resolution` | `"unresolved"` | `"resolved"` |
| `tests/test_review_regressions.py::test_the_synaptic_pulse_has_exactly_one_definition` | $K_d = k_{\text{off}} / k_{\text{on}}$ | `23.806` | `13.829` |

---

## 5. Irreconcilable Tests Requiring Decision (`kind=question`)

Per AGENT-BRIEF-13 §2 item 4:
> *"If any test's expectation cannot be reconciled with the sourced value, report `kind=question` rather than editing the test. A target change that quietly rewrites assertions is how a project loses its anchors."*

Two test suites cannot be reconciled with `FIT_TARGETS['ec50_um'] = 11.6` without policy decisions reserved for Opus:

### Issue A: `tests/test_manuscript_consistency.py` (11 failing tests)
The entire text of `knowledge/08-manuscript.md` was drafted against the $20.0\ \mu\text{M}$ calibration anchor. `test_manuscript_consistency.py` directly parses `08-manuscript.md` and asserts string equality for:
1. `test_calibration_anchors_match_the_code`: fails on reproduced EC50 string (`"11.6000"` vs manuscript `"20.0000"`).
2. `test_microscopic_rates_and_derived_constants_match`: fails on $k_{\text{on}}$ (`"0.0130384"` vs `"0.0075785"`), $K_d$ (`"13.83"` vs `"23.81"`), $E$ (`"3.9189"` vs `"3.9249"`).
3. `test_the_two_asymptotes_and_the_dynamic_range_match`: fails on headline dynamic range / headroom (`"43.3×"` vs `"123.3×"`).
4. `test_table_4_compartment_divergence_matches`: fails on modulator shifts ($c_{\text{affinity}} = 3.4906$ vs $3.4674$).
5. `test_the_finite_modulator_table_matches`: fails on $s_{\max} = 1.381$ vs $1.411$.
6. `test_falsification_intervals_are_the_reachable_ones`: fails on pre-registered falsification bounds ($10.20$ vs $10.29$).
7. `test_the_headline_number_is_never_quoted_without_its_scope`: fails because `"43.3"` appears without surrounding qualification in manuscript text.
8. `test_the_sensitivity_table_matches_the_sweep`
9. `test_the_robust_ambient_claim_matches_the_computed_range`
10. `test_the_cap_crossing_claim_is_verified_not_asserted`
11. `test_ambient_sweep_and_the_cap_crossing_match`

**Why this cannot be reconciled autonomously:**
Editing `08-manuscript.md` to revise every headline number, table, and argument from a $123.3\times$ headroom to a $43.3\times$ headroom is an authorial decision that alters the manuscript itself. Editing `test_manuscript_consistency.py` to stop checking the manuscript violates the project's consistency guarantees.

### Issue B: `tests/test_review_regressions.py::test_falsification_intervals_are_set_from_reachable_gain_not_the_asymptote`
* **Assertion:**
  `assert reach[0.1] < 0.01 * head[0.1], ("a finite s_max=2.5 PAM should reach under 1% of the asymptote at 0.1 uM ambient")`
* **Under 20.0 µM:** Asymptotic headroom at 0.1 µM ambient was $2881\times$; reachable gain was $8.5\times$ ($0.3\% < 1\%$).
* **Under 11.6 µM:** Asymptotic headroom collapses to $650\times$; reachable gain is $11.58\times$. The ratio is **$1.78\%$**, which trips the `< 0.01` threshold.
* The test explicitly notes in its docstring: *"re-derive the falsification intervals, this test's premise has changed"*. Re-deriving pre-registered experimental falsification bounds requires planner sign-off.

---

## 6. Verification: Test Suite Tiers & Pytest Output

### Baseline on `main` (Post Brief 12 Merge):
`574 passed, 41 warnings in 235.38s (0:03:55)`

### Verification on Branch `claude/agent-brief-13`:
* **Task 1 Suite (`tests/test_nextgen_models.py` & `tests/test_nesting.py`):**
  * `test_kinetic_coop_connectivity_and_zeros`: PASSED
  * `test_kinetic_coop_model_protocol_and_behavior`: PASSED
  * `test_declared_nesting_topology`: PASSED
  * `test_nesting_embedding_derives_conventions_by_name`: PASSED
  * `test_kinetic_coop_embedding_reproduces_kinetic_two_open_likelihood`: PASSED ($|\Delta \ln L| < 10^{-5}$)
  * `test_nesting_on_total_objective[...-kinetic_two_open-kinetic_coop]`: PASSED
* **Full Fast Tier (`pytest -q -m "not slow"`):**
  * Result: `16 failed, 561 passed, 40 warnings in 135.11s (0:02:15)`
  * (All 16 failures stem directly from the EC50 target shift: 11 in `test_manuscript_consistency.py`, 1 in `test_crc_components.py`, 1 in `test_fit_determinism.py`, 1 in `test_p4_likelihood.py`, 2 in `test_review_regressions.py`).

### Verbatim Pytest Output Lines:
```
16 failed, 561 passed, 40 warnings in 135.11s (0:02:15)
```
And for the specific Task 2 regression check:
```
5 failed, 98 passed, 13 warnings in 535.94s (0:08:55)
```
And for Task 1 nesting suite:
```
14 passed, 8 warnings in 81.63s (0:01:21)
```
