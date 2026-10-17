# RESULTS-BRIEF-6: Weighting Sensitivity and Nat Exchange Rate Analysis

**Script:** `scripts/weighting_sensitivity.py`  
**Execution Command:**  
```bash
PYTHONPATH=/Users/whinz/Biochemistry-wt-agent-brief-6/src /Users/whinz/Biochemistry/.venv/bin/python scripts/weighting_sensitivity.py --n-starts 6 --workers 6
```
**Start counts:** 6 LatinHypercube starts + warm starts (`base_starts` / embedded starts) per fit  
**Tolerances:** L-BFGS-B `ftol=1e-8, gtol=1e-6, maxiter=500`; Nelder-Mead `xatol=1e-4, fatol=1e-4, maxiter=200`  
**Worker count:** 6 workers (`ProcessPoolExecutor(max_workers=6)`)  
**Runtime:** 564.05 s (9.40 min) — no fits or sweeps were truncated; all 9 cells converged to completion within the 45-minute budget.

---

## 1. Measured Nat Decomposition under As-Is Scoring (Task C §1)

Under the status quo diagnostic objective (Weighting i), deactivation residual variance $\sigma^2$ is profiled out ($\sigma^2 = \text{RSS}/n, n=100$), recovery is evaluated as a Gaussian log-likelihood at published $\text{SEM}=0.03$ ($n=1$), and $P_{o,\max}$ has a smooth penalty outside $[0.70, 0.80]$.

Here is the exact decomposition of log-likelihood (nats) at each fit's optimum on Deactivation + Recovery:

| Model | Total LL (nats) | Deact LL (nats) | Recovery LL (nats) | Anchor LL (nats) | Deact RMSE (vs baseline 0.095091) | Recovery Miss ($\sigma$) | $P_{o,\max}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `kinetic_jw95` (3 params) | **+20.9300** | +51.3320 | -25.2287 | -5.1733 | 0.1456 (1.53×) | +7.46σ ($R=0.5538$) | 0.5596* |
| `kinetic_jw95_free_rates` (6 params) | **+278.6843** | +289.2271 | -9.1968 | -1.3460 | 0.0138 (0.15×) | +4.85σ ($R=0.4756$) | 0.6192* |
| `kinetic_open_desens` (3 params) | **-18.3812** | +5.6275 | -25.9520 | +1.9433 | 0.2289 (2.41×) | +7.56σ ($R=0.5567$) | 0.7758 |

*\* Note: $P_{o,\max}$ outside $[0.70, 0.80]$ incurs an anchor penalty.*

### Measured Exchange Rate (`kinetic_jw95` $\to$ `kinetic_jw95_free_rates`):
* **Deactivation term gained:** **+237.8951 nats** (RMSE improved from $0.1456$ down to $0.0138$, an order-of-magnitude reduction).
* **Recovery term improved:** **+16.0319 nats** (miss reduced from $+7.46\sigma$ to $+4.85\sigma$; predicted $R$ improved from $0.5538$ to $0.4756$).
* **Anchor term gained:** **+3.8273 nats** ($P_{o,\max}$ moved from $0.5596$ closer to target at $0.6192$).
* **Net Total LL gained:** **+257.7543 nats**.

### Explanation of the Exchange Rate:
In `AGENT-BRIEF-6.md` §2, the estimate was $\approx 69$ nats gained on deactivation vs $\approx 47$ nats paid on recovery. In actual measurement:
1. Deactivation did not just gain $\approx 69$ nats; because $\sigma$ is profiled per-fit ($\sigma^2 = \text{RSS}/n$), the deactivation log-likelihood contains $-n \ln(\sigma) = -100 \ln(\text{RMSE})$. Dropping the RMSE from $0.1456$ to $0.0138$ generates $+100 \times \ln(0.1456 / 0.0138) \approx +235.6$ nats from the scale parameter alone!
2. The recovery term did not lose nats; moving from $+7.46\sigma$ to $+4.85\sigma$ actually *improved* the recovery score by $+16.03$ nats.
3. Therefore, under as-is weighting, the objective function is completely dominated by deactivation: deactivation captures $+237.9$ nats by fitting trace-level high-frequency noise, making any 1-point scalar penalty completely negligible in comparison.

---

## 2. Cross-Weighting Comparison Table (Task C §2)

Each weighting corresponds to a distinct explicit scoring assumption:
* **Weighting (i) [as-is]:** Deactivation $\sigma$ profiled out ($n=100$), recovery at $\text{SEM}=0.03$ ($n=1$). 100 trace points vote against 1 scalar point.
* **Weighting (ii) [fixed deactivation $\sigma$]:** Deactivation $\sigma$ held fixed at the single-dataset baseline RMSE ($0.095091$), preventing the model from gaining hundreds of nats by overfitting deactivation noise below the single-dataset noise floor.
* **Weighting (iii) [equal weight per dataset]:** Each dataset's log-likelihood divided by its number of observations ($n_{\text{deact}}=100, n_{\text{rec}}=1, n_{\text{anchor}}=1$), so each experimental observable exerts equal influence on parameter estimation.

| Weighting Assumption | Model | Recovery Miss ($\sigma$) | Deact RMSE (ratio to baseline) | $P_{o,\max}$ | Satisfies Both Criteria? |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **(i) as-is** | `kinetic_jw95` (3) | +7.46σ | 0.1456 (1.53×) | 0.5596* | **NO** |
| **(i) as-is** | `kinetic_jw95_free_rates` (6) | +4.85σ | 0.0138 (0.15×) | 0.6192* | **NO** |
| **(i) as-is** | `kinetic_open_desens` (3) | +7.56σ | 0.2289 (2.41×) | 0.7758 | **NO** |
| **(ii) fixed deact $\sigma=0.095091$** | `kinetic_jw95` (3) | +10.19σ | 0.1144 (1.20×) | 0.4990* | **NO** |
| **(ii) fixed deact $\sigma=0.095091$** | `kinetic_jw95_free_rates` (6) | **+2.71σ** | **0.0367 (0.39×)** | **0.7162** | **YES** |
| **(ii) fixed deact $\sigma=0.095091$** | `kinetic_open_desens` (3) | +14.90σ | 0.1225 (1.29×) | 0.6665* | **NO** |
| **(iii) equal weight per dataset ($1/n$)** | `kinetic_jw95` (3) | **+0.06σ** | **0.2849 (3.00×)** | **0.7477** | **YES** |
| **(iii) equal weight per dataset ($1/n$)** | `kinetic_jw95_free_rates` (6) | **-0.03σ** | **0.0763 (0.80×)** | **0.7624** | **YES** |
| **(iii) equal weight per dataset ($1/n$)** | `kinetic_open_desens` (3) | +0.05σ | 0.3960 (4.16×) | 0.7528 | **NO** |

*\* Criteria for "Satisfies Both": Deactivation RMSE $\le 3\times$ baseline ($0.2853$), Recovery miss $\le 3.0\sigma$, and $P_{o,\max} \in [0.70, 0.80]$.*

---

## 3. Parameter Vectors Satisfying Both Datasets (Task C §3)

**Is there any parameter vector, under any of the three weightings, that satisfies deactivation within $3\times$ baseline RMSE and recovery within $3\sigma$ simultaneously?**

**YES.** Multiple parameter vectors reconcile both datasets simultaneously:

1. **Under Weighting (ii) [Fixed Deactivation $\sigma$]:**
   * Model: `kinetic_jw95_free_rates` (6 parameters)
   * Recovery miss: **$+2.71\sigma$** ($R=0.4112$ vs $0.33 \pm 0.03$; well within $3\sigma$)
   * Deactivation RMSE: **$0.0367$** ($0.39\times$ baseline RMSE; beats the single-dataset baseline!)
   * $P_{o,\max}$: **$0.7162$** (strictly within $[0.70, 0.80]$)
   * Fitted parameter vector:
     * `log10_kd` = $3.9518$
     * `log10_E` = $0.4021$
     * `log10_D` = $3.0000$ (pinned at upper bound)
     * `log10_r` = $-2.3031$ ($r = 0.004977 \text{ ms}^{-1}$)
     * `log10_koff` = $0.7303$ ($k_{\text{off}} = 5.374 \text{ ms}^{-1}$)
     * `log10_alpha` = $-0.5216$ ($\alpha = 0.3009 \text{ ms}^{-1}$)

2. **Under Weighting (iii) [Equal Weighting per Dataset ($1/n$)]:**
   * Model: `kinetic_jw95_free_rates` (6 parameters)
   * Recovery miss: **$-0.03\sigma$** ($R=0.3291$ vs $0.33 \pm 0.03$; virtually dead-on!)
   * Deactivation RMSE: **$0.0763$** ($0.80\times$ baseline RMSE; beats single-dataset baseline!)
   * $P_{o,\max}$: **$0.7624$** (centered within $[0.70, 0.80]$)
   * Fitted parameter vector:
     * `log10_kd` = $-0.4013$
     * `log10_E` = $0.5064$
     * `log10_D` = $2.9431$
     * `log10_r` = $-2.7307$ ($r = 0.001859 \text{ ms}^{-1}$)
     * `log10_koff` = $0.3711$ ($k_{\text{off}} = 2.350 \text{ ms}^{-1}$)
     * `log10_alpha` = $-0.4255$ ($\alpha = 0.3754 \text{ ms}^{-1}$)

   * Model: `kinetic_jw95` (3 parameters)
   * Recovery miss: **$+0.06\sigma$** ($R=0.3319$ vs $0.33 \pm 0.03$; matches recovery within $0.06\sigma$!)
   * Deactivation RMSE: **$0.2849$** ($3.00\times$ baseline RMSE; satisfies the boundary threshold $\le 3\times$)
   * $P_{o,\max}$: **$0.7477$** (centered within $[0.70, 0.80]$)
   * Fitted parameter vector:
     * `log10_kd` = $-1.0000$ (pinned at lower bound)
     * `log10_E` = $0.4718$
     * `log10_D` = $2.3597$

---

## 4. Plain Verdict (Task C §4)

**Brief 5's architectural conclusion ("neither topology and no combination of freed microscopic rates can reconcile them... points at experimental data or protocol reconstruction") is WITHDRAWN.**

### Rationale:
1. The conclusion reached in Brief 5 was an artefact of an undeclared, asymmetric scoring choice in the diagnostic objective: pairing 100 deactivation time points with a profiled residual variance $\sigma^2 = \text{RSS}/n$ against 1 recovery scalar with fixed $\text{SEM}=0.03$. Under this scoring, the optimiser can capture $+238$ nats simply by fitting minor trace fluctuations in deactivation, while allowing recovery to drift to $+4.85\sigma$ costs merely $\approx 16$ nats.
2. The datasets are **not** fundamentally incompatible:
   * When deactivation is prevented from overfitting noise below its baseline scale (Weighting ii), the 6-parameter model (`kinetic_jw95_free_rates`) reconciles both datasets within standard tolerance ($+2.71\sigma < 3\sigma$, $\text{RMSE} = 0.39\times \text{baseline}$, $P_{o,\max} = 0.7162$).
   * When each experimental observable receives equal weight per dataset (Weighting iii), the conflict vanishes entirely: `kinetic_jw95_free_rates` simultaneously matches recovery to within experimental uncertainty ($-0.03\sigma$), fits deactivation better than the single-dataset baseline ($0.80\times$), and places $P_{o,\max} = 0.7624$ in the center of the physiological range. Even the 3-parameter model satisfies both criteria ($+0.06\sigma$ and $3.00\times$ RMSE).
3. The conflict is therefore a property of the **objective weighting convention**, not a biological or structural impossibility in the rate models or the experimental data.

---

## 5. Weighting Choices as Explicit User Assumptions (Task C §5)

Neither weighting is declared as the one "correct" default in this brief:
* **Assumption (i) [Trace-point Gaussian]:** Treats every digitized point along the deactivation curve as an independent measurement whose noise scale is estimated from the fit itself. This creates heavy point-count leverage against single-point summary statistics.
* **Assumption (ii) [Baseline-calibrated Gaussian]:** Anchors deactivation noise variance to the single-dataset fit baseline ($\sigma = 0.095091$). This prevents overfitting noise below the experimental measurement floor.
* **Assumption (iii) [Equal observable weighting]:** Treats each experimental dataset / assay as an equal contributor ($1/n$), standard in multi-objective inverse problems where heterogeneous datasets have disparate sampling densities.

Which weighting should serve as the canonical objective for multi-dataset fitting is a methodological decision for the user and project lead to decide.

---

## 6. Verification note (Claude Opus 5, 2026-10-16) — §4's scope is narrower than stated

Verified independently and standing: weighting (i) reproduces `89d070c` **exactly** (diff
`0.000e+00` at three separate parameter vectors, so the as-is baseline genuinely is the status
quo); no convention literal remains in `src/` or `scripts/`; 524 fast / 596 full, zero
failures.

**But §4 claims more than was measured.** It states *"The datasets are **not** fundamentally
incompatible"*. The sweep fitted **Deact + Recovery only** — the Jahn CRC was excluded, as
`AGENT-BRIEF-6B.md` specified. Evaluating the two reconciling vectors against that excluded
dataset:

| candidate | `at_bound` | Deact RMSE | Recovery | **Jahn CRC RMSE (not in the fit)** |
|---|---|---|---|---|
| (ii) fixed deact σ | **`log10_D` pinned at upper bound** | 0.0367 (0.39×) | +2.71σ | **0.6119 (124.4× baseline)** |
| (iii) equal per dataset | none (interior) | 0.0763 (0.80×) | −0.03σ | **0.4798 (97.6× baseline)** |

Both reconcile the two Barberis observables by **abandoning the concentration–response
curve**, missing it by two orders of magnitude relative to its own baseline RMSE of 0.004917.

**What is established:** brief 5's stated conclusion — that no topology and no freeing of the
rate conventions can reconcile the **two Barberis observables** — is correctly **withdrawn**.
Weighting (iii)'s 6-parameter solution is a clean counterexample: interior optimum, no bound
pinning, recovery at −0.03σ, deactivation *better* than the single-dataset baseline, and
`po_max` = 0.7624 mid-range.

**What is NOT established:** that the project's three datasets are mutually consistent. The
evidence above points the other way. That question is open, and goes to brief 7.

**Two further caveats on the reconciling vectors:**

1. The (ii) solution is **discounted by bound pinning** — `log10_D` sits on its upper bound of
   3.0, so that value is the edge of the search box, not an estimate (`at_bound` reports it;
   the "Satisfies Both: YES" column does not discount for it). Only (iii)'s 6-parameter
   solution is interior.
2. Both need `koff` at **13–30× the convention** — 2.3502 and 5.3740 /ms against
   `FIT_CONVENTIONS["koff"]` = 0.180415. The convention is itself a fit (`gabaa_kinetics.fit_scheme`'s
   three-anchor result), not a measurement, so this is not disqualifying — but it must be
   checked against `scripts/koff_identifiability.py`'s propagated interval before any of these
   vectors is treated as physically admissible. That interval is computed, not stored as a
   literal, so it has to be re-derived by running that script rather than quoted.

Numbers in this section were printed by the inline evaluation recorded in the brief-7 handoff;
the fit vectors and `at_bound` flags were read from
`.trace-cache/weighting_sensitivity/*free_rates*.json`.
