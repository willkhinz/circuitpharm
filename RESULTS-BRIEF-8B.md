# RESULTS-BRIEF-8B: Does a Pre-Open Conformational Step Reproduce Measured Deactivation?

**Agent:** 8B (Pre-Open / Flip Scheme Investigation)  
**Worktree:** `/Users/whinz/Biochemistry-wt-agent-brief-8b`  
**Branch:** `claude/agent-brief-8b`  
**Committed Search Script:** `scripts/flip_structural_search.py`  
**Execution Command:**  
```bash
PYTHONPATH=/Users/whinz/Biochemistry-wt-agent-brief-8b/src /Users/whinz/Biochemistry/.venv/bin/python scripts/flip_structural_search.py
```
**Runtime:** 1.97 s total (JW95: 0.78 s; Flip: 1.13 s)  
**Search Budget:** 1500 Latin Hypercube samples, top 12 polished via bounded Nelder–Mead, random seed 42.

---

## 1. Executive Summary and Direct Answer to the Brief's Question

> **Question:** Does adding a pre-open conformational step — a shut, fully-liganded intermediate `A2F` between `A2R` and `A2O` (the "flip" mechanism of Lape et al. 2008) — produce a fast, large-amplitude component and so reproduce the published three-component deactivation of $\alpha1\beta2\gamma2$ (Barberis et al. 2007)?

### Plain Answer:
**NO. The pre-open flip mechanism FAILS to reproduce the published deactivation, displaying the exact same structural signature of failure observed in the existing schemes.**

Specifically:
1. **No meaningful improvement in $\chi^2$:** Over all six published deactivation numbers (read directly from `circuitpharm.literature.DEACTIVATIONS[1]` with no literals), the best $\chi^2$ for `kinetic_flip` is **527.4** (against **527.9** for `kinetic_jw95_free_rates`, against $\sim 6$ for a genuine fit).
2. **Persistent amplitude collapse on the fastest mode:** In `kinetic_flip`, the predicted fastest mode amplitude collapses to $\text{area}_1 = -0.002906$, missing the published $0.57 \pm 0.04$ by **$-14.32\sigma$**. In both schemes, the published data require the fastest component ($\tau = 2.8\text{ ms}$) to carry the dominant amplitude ($0.57$), but the Markov relaxation topology forces this mode's amplitude to virtually zero (or slightly negative).
3. **Only two appreciable modes:** In `kinetic_flip`, only **2 of 5** relaxation modes carry $|\text{amplitude}| \ge 0.10$ (mode 4: $40.6\%$ at $\tau = 33.26\text{ ms}$; mode 5: $59.7\%$ at $\tau = 222.49\text{ ms}$). The fastest mode carries only $-0.29\%$.
4. **Three positive amplitudes are not obtained:** The fastest amplitude is negative ($-0.0029$), failing the physical requirement of a pure multi-exponential decay.
5. **Cross-checks fail decisively:** Even at this deactivation optimum, `kinetic_flip` misses the Jahn 1997 peak CRC baseline RMSE by **$115.0\times$** ($0.5655$ vs $0.004917$), and its maximum open probability is severely depressed ($P_{o,\max} = 0.2366$, far below the physiological band $[0.70, 0.80]$).

This establishes a sharp negative result: inserting a pre-open conformational equilibrium between agonist binding and channel opening **does not resolve the deactivation shape failure**.

---

## 2. Cross-Scheme Search Comparison (Same Code Path & Budget)

Evaluated under the exact same optimization protocol: 1500 Latin Hypercube starts, best 12 polished via bounded Nelder–Mead, minimizing $\chi^2$ over the six published deactivation parameters from Barberis 2007 Table 1 (`DEACTIVATIONS[1]`):

$$\chi^2 = \sum_{i=1}^3 \left(\frac{\tau_i - \tau_i^{\text{pub}}}{\sigma_{\tau, i}}\right)^2 + \sum_{i=1}^3 \left(\frac{\text{area}_i - \text{area}_i^{\text{pub}}}{\sigma_{\text{area}, i}}\right)^2$$

| Scheme | Free Params | Best $\chi^2$ (6 items) | Modes $\ge 0.10$ | 3 Positive Amplitudes? | Pinned at Bound (`at_bound`) |
|---|---|---|---|---|---|
| `kinetic_jw95_free_rates` | 6 | **527.9** | 2 of 4 | no | `log10_kd (L)` |
| `kinetic_flip` | 8 | **527.4** | 2 of 5 | no | `log10_D (U)` |

---

## 3. Detailed Component Breakdown and Residuals

Published reference targets: `DEACTIVATIONS[1]` (Barberis et al. 2007, 2 ms pulse of 10 mM GABA in $\alpha1\beta2\gamma2$):
* $\tau = (2.80 \pm 0.30, 33.40 \pm 4.60, 221.35 \pm 14.90)\text{ ms}$
* $\text{areas} = (0.57 \pm 0.04, 0.23 \pm 0.02, 0.20 \pm 0.03)$

### Model 1: `kinetic_jw95_free_rates` (Best $\chi^2 = 527.86$)

| Quantity | Published Target | Published Error | Predicted Value | Residual ($\sigma$) |
|---|---|---|---|---|
| $\tau_1$ (ms) | 2.80 | $\pm 0.30$ | 0.2970 | -8.34$\sigma$ |
| $\tau_2$ (ms) | 33.40 | $\pm 4.60$ | 33.5386 | +0.03$\sigma$ |
| $\tau_3$ (ms) | 221.35 | $\pm 14.90$ | 222.8265 | +0.10$\sigma$ |
| $\text{area}_1$ | 0.57 | $\pm 0.04$ | -0.003296 | **-14.33$\sigma$** |
| $\text{area}_2$ | 0.23 | $\pm 0.02$ | 0.405744 | +8.79$\sigma$ |
| $\text{area}_3$ | 0.20 | $\pm 0.03$ | 0.597551 | +13.25$\sigma$ |

* **All 4 relaxation modes:**
  * $\tau$: `[0.2970, 33.54, 222.2, 222.8]` ms
  * $\text{areas}$: `[-0.00330, 0.40574, -0.00000, 0.59755]`
* **Fitted parameters:**
  * `log10_kd` = -1.0000 (value: 0.1000 $\mu$M) $\leftarrow$ **PINNED AT LOWER BOUND**
  * `log10_E` = -0.5270 (value: 0.2971)
  * `log10_D` = -0.2650 (value: 0.5433)
  * `log10_r` = -1.7142 (value: 0.01931 $\text{ms}^{-1}$)
  * `log10_koff` = -2.3468 (value: 0.00450 $\text{ms}^{-1}$)
  * `log10_alpha` = 0.4137 (value: 2.592 $\text{ms}^{-1}$)

---

### Model 2: `kinetic_flip` (Best $\chi^2 = 527.41$)

| Quantity | Published Target | Published Error | Predicted Value | Residual ($\sigma$) |
|---|---|---|---|---|
| $\tau_1$ (ms) | 2.80 | $\pm 0.30$ | 0.2934 | -8.36$\sigma$ |
| $\tau_2$ (ms) | 33.40 | $\pm 4.60$ | 33.2583 | -0.03$\sigma$ |
| $\tau_3$ (ms) | 221.35 | $\pm 14.90$ | 222.4881 | +0.08$\sigma$ |
| $\text{area}_1$ | 0.57 | $\pm 0.04$ | -0.002906 | **-14.32$\sigma$** |
| $\text{area}_2$ | 0.23 | $\pm 0.02$ | 0.406070 | +8.80$\sigma$ |
| $\text{area}_3$ | 0.20 | $\pm 0.03$ | 0.596811 | +13.23$\sigma$ |

* **All 5 relaxation modes:**
  * $\tau$: `[0.03716, 0.2934, 0.3288, 33.26, 222.5]` ms
  * $\text{areas}$: `[+0.00003, -0.00291, -0.00000, +0.40607, +0.59681]`
* **Fitted parameters:**
  * `log10_kd` = 0.1094 (value: 1.2866 $\mu$M)
  * `log10_E` = -0.4959 (value: 0.3193)
  * `log10_D` = 3.0000 (value: 1000.0) $\leftarrow$ **PINNED AT UPPER BOUND**
  * `log10_r` = -1.7110 (value: 0.01945 $\text{ms}^{-1}$)
  * `log10_koff` = 0.4830 (value: 3.0410 $\text{ms}^{-1}$)
  * `log10_alpha` = 0.4110 (value: 2.5765 $\text{ms}^{-1}$)
  * `log10_delta` = 0.1336 (value: 1.3602 $\text{ms}^{-1}$)
  * `log10_gamma` = -1.3908 (value: 0.04067 $\text{ms}^{-1}$)

---

## 4. Task 3: Independent Evidence Cross-Checks

Even though deactivation is NOT fitted ($\chi^2 = 527.4$, far from $\sim 6$), Brief 8B requires evaluating the optimum against the other project anchors without cherry-picking:

| Scheme | Paired-Pulse (100 ms) | PP $\sigma$ Miss | Jahn Peak CRC RMSE | Jahn RMSE Ratio | $P_{o,\max}$ | $P_{o,\max} \in [0.70, 0.80]$? |
|---|---|---|---|---|---|---|
| `kinetic_jw95_free_rates` | 0.4050 (pub 0.33 $\pm$ 0.03) | +2.50$\sigma$ | 0.416746 | 84.8$\times$ | 0.2291 | **NO** |
| `kinetic_flip` | 0.3625 (pub 0.33 $\pm$ 0.03) | +1.08$\sigma$ | 0.565454 | **115.0$\times$** | 0.2366 | **NO** |

### Observations:
1. **Paired-pulse recovery:** `kinetic_flip` reaches $R = 0.3625$ ($+1.08\sigma$), close to the published $0.33 \pm 0.03$. However, as established in Brief 6, matching recovery at a deactivation optimum is an artifact of trading away trace shape.
2. **Jahn 1997 peak CRC:** `kinetic_flip` has an RMSE of **$0.565454$**, which is **$115.0\times$** higher than the baseline single-dataset fit RMSE of $0.004917$. This matches the failure mode in Brief 6 where reconciling vectors missed Jahn by $97.6\times$ and $124.4\times$.
3. **Open probability ceiling:** $P_{o,\max} = 0.2366$, which fails the physiological calibration band $[0.70, 0.80]$ by nearly a factor of 3.

---

## 5. Mathematical Mechanics and Nesting Analysis

### Why the Pre-Open State Fails to Supply a Fast Large-Amplitude Component:
In the sequential scheme:
$$R \underset{k_{\text{off}}}{\overset{2 k_{\text{on}} g}{\rightleftharpoons}} AR \underset{2 k_{\text{off}}}{\overset{k_{\text{on}} g}{\rightleftharpoons}} A2R \underset{\gamma}{\overset{\delta}{\rightleftharpoons}} A2F \underset{\alpha}{\overset{\beta}{\rightleftharpoons}} A2O$$
with desensitisation $A2R \underset{r}{\overset{d}{\rightleftharpoons}} A2D$.

During deactivation at zero agonist ($g = 0$), receptors in $A2O$ must relax through $A2F \to A2R$ before unbinding ($A2R \to AR \to R$).
* The transition $A2F \leftrightarrow A2O$ is governed by $\beta$ and $\alpha$. Because channel closing is fast ($\alpha \approx 2.58\text{ ms}^{-1}$, $\tau \sim 0.3\text{ ms}$), the $A2O \to A2F$ transition quickly reaches quasi-equilibrium.
* However, the projection of this rapid relaxation onto the observable $A2O$ current depends on the initial occupancy of $A2F$ relative to $A2O$ at the end of the 2 ms pulse. Because the pulse is brief (2 ms) and $P_{o,\max}$ is constrained, $A2F$ cannot store a massive reservoir of shut receptors that re-feed $A2O$ with large positive net amplitude without simultaneously corrupting the intermediate and slow phases.
* Consequently, the optimizer drives the fast mode's amplitude to $-0.0029$ and allocates all observable deactivation amplitude ($40.6\%$ and $59.7\%$) to the two slower modes ($\tau \approx 33\text{ ms}$ and $\tau \approx 222\text{ ms}$). The fast $2.8\text{ ms}$ mode carrying $57\%$ of the area cannot be produced by this topology.

### Nesting Against `kinetic_jw95`:
In `src/circuitpharm/fitting/comparison.py`, `kinetic_flip` is declared with `nests = ()` (excluded from strictly nesting `kinetic_jw95`):
* **Finite search box:** Collapsing $A2F$ onto $A2R$ requires $\delta \to \infty$ with $\gamma$ fixed. In any bounded search box ($\log_{10}\delta \le 3.0$), $\delta$ sits on the boundary ($1000\text{ ms}^{-1}$) rather than infinity.
* **Desensitisation bypass:** Because desensitisation enters strictly from $A2R$ ($A2R \leftrightarrow A2D$), rapid flipping ($\delta \gg \gamma$) shifts occupancy away from $A2R$ into $A2F$, artificially suppressing desensitisation entry.
* **Measured discrepancy:** At $\delta = 1000$, $\gamma = 1$ with JW95 conventions ($r = 0.002, k_{\text{off}} = 0.180415, \alpha = 0.30$), the measured discrepancy on Deact+Recovery is **$161.60\text{ nats}$** (and **$3.37\text{ nats}$** on Jahn peak). This is verified by `test_kinetic_flip_exclusion_with_measured_discrepancy` in `tests/test_nesting.py`.

---

## 6. Acceptance Criteria Verification

- [x] `src/circuitpharm/models/kinetic_flip.py` implemented with the stated 6-state topology and documented state order (`open_index = 4`).
- [x] Connectivity test `test_kinetic_flip_connectivity_and_zeros` in `tests/test_nextgen_models.py` with explicit negative assertions (all zero transitions asserted).
- [x] `kinetic_flip` entry in `MODEL_SPECS` in `src/circuitpharm/fitting/comparison.py` (additive edit); nesting against `kinetic_jw95` excluded honestly with measured discrepancy ($161.60\text{ nats}$) recorded and tested in `tests/test_nesting.py`.
- [x] `scripts/flip_structural_search.py` committed; no published value or error hardcoded (all read from `DEACTIVATIONS[1]`).
- [x] `RESULTS-BRIEF-8B.md` reports $\chi^2$ for `kinetic_flip` ($527.4$) and `kinetic_jw95_free_rates` ($527.9$) from the same code path and budget, with per-quantity $\sigma$, mode counts, and bound flags.
- [x] Task 3 cross-checks fully evaluated and reported: paired-pulse ($+1.08\sigma$), Jahn CRC ($115.0\times$ RMSE), and $P_{o,\max}$ ($0.2366$, failing $[0.70, 0.80]$).
- [x] Both test tiers verified with `PYTHONPATH` set, counts labelled.

---

## 7. Test Suite Verification

Both test tiers were executed with `PYTHONPATH=/Users/whinz/Biochemistry-wt-agent-brief-8b/src` and interpreter `/Users/whinz/Biochemistry/.venv/bin/python`:

1. **Fast Test Tier (`pytest -m "not slow" -q`):**
   * **Result:** `544 passed, 41 warnings in 244.71s (0:04:04)`
   * **Comparison:** Baseline was 524 on `main` and ~535 on the `7A` branch. All 544 tests passed with zero failures.

2. **Full Test Tier (`pytest -q`):**
   * **Result:** `616 passed, 124 warnings in 1064.29s (0:17:44)`
   * **Comparison:** Baseline was 596 on `main`. All 616 tests passed with zero failures.
   * **Discipline:** The known Apple Silicon split-R̂ artefact in `test_posterior_recovers_known_parameters` was not touched or gated.
