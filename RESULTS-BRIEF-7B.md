# RESULTS-BRIEF-7B: Deactivation Components by Waveform Fitting

**Agent:** Gemini 3.8 Flash (High) (Worker 7B)  
**Assigned Method:** Empirical waveform simulation and triexponential fitting (reproducing Barberis et al. 2007 experimental procedure)  
**Execution Command:**  
```bash
PYTHONPATH=/Users/whinz/Biochemistry-wt-agent-brief-7b/src /Users/whinz/Biochemistry/.venv/bin/python scripts/deact_components_brief7b.py
```
**Runtime:** 119.96 s (2.00 min)  
**Interpreter:** `/Users/whinz/Biochemistry/.venv/bin/python`  
**Git Branch:** `claude/agent-brief-7b` in `/Users/whinz/Biochemistry-wt-agent-brief-7b`  

---

## 1. Executive Summary & The Conflict Verdict

### The Question Brief 7 Was Built to Answer:
> *Scored against the published deactivation parameters at their published errors, rather than against a reconstructed curve, does the Barberis deactivation / paired-pulse recovery conflict exist at all?*

### Plain Answer:
**YES, the conflict persists and is dramatically magnified.**

In Brief 6, Weighting (iii) (`equal_per_dataset`, $1/n$ weighting per dataset) appeared to reconcile the Barberis deactivation trace and paired-pulse recovery under the 6-parameter model (`kinetic_jw95_free_rates`), matching recovery at $-0.03\sigma$ and yielding a deactivation curve RMSE of $0.0763$ ($0.80\times$ baseline).

However, when scored under **Weighting (iv)** directly against the attributable measurements and their published error bars (Barberis 2007 Table 1):
1. **The slow deactivation component is severely violated:** At the Weighting (iii) optimum, the model produces $\tau_3 = 743.89\text{ ms}$ against the published $221.35 \pm 14.9\text{ ms}$ — an enormous **$+35.07\sigma$ miss**.
2. **The intermediate component is absent:** The model has no intermediate relaxation mode near $33.4\text{ ms}$; the triexponential fit forces $\tau_1$ and $\tau_2$ to collapse together at $3.48\text{ ms}$, missing published $\tau_2 = 33.4 \pm 4.6\text{ ms}$ by **$-6.50\sigma$**.
3. **The likelihood penalty is immense:** The predicted deactivation components incur a log-likelihood penalty of **$-684.64\text{ nats}$**. The $+2.59\text{ nats}$ gained by matching recovery ($R = 0.3291$) is utterly dwarfed by the hundreds of nats lost on deactivation kinetics.
4. **Conclusion:** The apparent reconciliation under curve-level RMSE in Brief 6 was an artifact of averaging squared residuals over a $1000\text{ ms}$ window where the slow tail has tiny absolute amplitude ($< 0.05$). When scored honestly against the published parameters with their stated errors ($n=6$, no pseudo-data expansion), the 5-state Markov models cannot simultaneously reconcile deactivation and paired-pulse recovery.

---

## 2. Deliverable 3: Ground-Truth Recovery Benchmark

To verify that our predictor is exact, `fit_triexponential_decay` was evaluated directly on `BARBERIS2007_DEACT.normalized_current` (generated from the six published numbers):

| Quantity | Published Target | Recovered Value | Error / Discrepancy |
| :--- | :--- | :--- | :--- |
| $\tau_1$ | $2.80\text{ ms}$ | $2.8000\text{ ms}$ | $< 10^{-12}\text{ ms}$ |
| $\tau_2$ | $33.40\text{ ms}$ | $33.4000\text{ ms}$ | $< 10^{-12}\text{ ms}$ |
| $\tau_3$ | $221.35\text{ ms}$ | $221.3500\text{ ms}$ | $< 10^{-12}\text{ ms}$ |
| $a_1$ | $0.570$ | $0.5700$ | $< 10^{-12}$ |
| $a_2$ | $0.230$ | $0.2300$ | $< 10^{-12}$ |
| $a_3$ | $0.200$ | $0.2000$ | $< 10^{-12}$ |
| $\tau_w$ (derived) | $52.5 \pm 2.9\text{ ms}$ | $53.5480\text{ ms}$ | $+0.36\sigma$ from published |
| Residual Cost | — | $7.78 \times 10^{-32}$ | Machine precision zero ($R^2 = 1.00000000$) |
| Condition Number | — | $39.01$ | Healthy, well-conditioned Jacobian |
| Singular Values | — | $(5.4678, 1.9620, 0.9188, 0.3719, 0.1402)$ | All 5 modes active |
| Component LL | — | $+2.10\text{ nats}$ | Maximal ($z_{\tau} = 0, z_{a} = 0$) |

**Verdict:** The triexponential fitting engine recovers the true parameters to machine precision with a well-conditioned Jacobian ($s_{\min} = 0.1402, \kappa = 39.01$). Any discrepancy observed on model waveforms is a property of the model's kinetics, not optimization failure.

---

## 3. Comprehensive Evaluation Across All 9 Brief-6 Optima

### Summary Comparison Table:

| Model | Weighting | $\tau_1$ (ms) | $\tau_2$ (ms) | $\tau_3$ (ms) | $a_1$ | $a_2$ | $a_3$ | Rec Miss ($\sigma$) | LL Comp (nats) | LL Rec (nats) | Total (iv) LL (nats) | Well-Cond? |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `kinetic_jw95` (3) | (i) as-is | $8.57$ | $8.57$ | $821.22$ | $0.378$ | $0.292$ | $0.330$ | $+7.46\sigma$ | $-1103.3$ | $-25.2$ | $-1128.6$ | **False** |
| `kinetic_jw95` (3) | (iii) equal $1/n$ | $8.83$ | $374.72$ | $1165.89$ | $0.432$ | $0.009$ | $0.559$ | $+0.06\sigma$ | $-5105.1$ | $+2.6$ | $-5102.5$ | **False** |
| `kinetic_jw95` (3) | (ii) fixed $\sigma$ | $8.58$ | $8.58$ | $771.93$ | $0.662$ | $0.081$ | $0.256$ | $+10.19\sigma$ | $-913.0$ | $-49.3$ | $-962.3$ | **False** |
| `kinetic_jw95_free_rates` (6) | (i) as-is | $3.33$ | $169.21$ | $179.04$ | $0.651$ | $1.000$ | $-0.651$ | $+4.85\sigma$ | $-1585.2$ | $-9.2$ | $-1594.4$ | **False** |
| `kinetic_jw95_free_rates` (6) | (iii) equal $1/n$ | $3.48$ | $3.48$ | $743.89$ | $0.355$ | $0.390$ | $0.254$ | $-0.03\sigma$ | $-684.6$ | $+2.6$ | $-682.0$ | **False** |
| `kinetic_jw95_free_rates` (6) | (ii) fixed $\sigma$ | $3.62$ | $18.60$ | $296.68$ | $0.682$ | $0.000$ | $0.318$ | $+2.71\sigma$ | $-97.4$ | $-1.1$ | $-98.5$ | **False** |
| `kinetic_open_desens` (3) | (i) as-is | $9.34$ | $630.96$ | $986.32$ | $0.527$ | $0.000$ | $0.473$ | $+7.56\sigma$ | $-10099.4$ | $-26.0$ | $-10125.4$ | **False** |
| `kinetic_open_desens` (3) | (iii) equal $1/n$ | $5.45$ | $5.45$ | $1764.47$ | $0.174$ | $0.119$ | $0.707$ | $+0.05\sigma$ | $-5625.5$ | $+2.6$ | $-5622.9$ | **False** |
| `kinetic_open_desens` (3) | (ii) fixed $\sigma$ | $9.81$ | $9.88$ | $748.96$ | $0.000$ | $0.744$ | $0.256$ | $+14.90\sigma$ | $-1343.5$ | $-108.4$ | $-1451.9$ | **False** |

*Published target reference: $\tau = (2.8 \pm 0.3, 33.4 \pm 4.6, 221.35 \pm 14.9)\text{ ms}$, areas $= (0.57 \pm 0.04, 0.23 \pm 0.02, 0.20 \pm 0.03)$, recovery $= 0.33 \pm 0.03$.*

---

## 4. Why the Triexponential Fit is Ill-Conditioned Across All Models

A central finding of Brief 7B is that **every single fit across all 9 cells is ill-conditioned (`Well-Cond: False`)**:

### The Mathematical Mechanism:
The JW95 scheme has 5 states, which under zero GABA (washout) yields **4 non-zero relaxation eigenvalues**. However:
1. **Eigenvalue clustering:** In `kinetic_jw95`, the analytical washout eigenvalues yield relaxation time constants of approximately $[0.93, 5.54, 6.53, 762.1]\text{ ms}$.
   - The fastest mode ($0.93\text{ ms}$) decays almost instantly within $< 2\text{ ms}$.
   - The two middle modes ($5.54\text{ ms}$ and $6.53\text{ ms}$) are separated by only $15\%$, merging macroscopically into a single apparent decay near $8.5\text{ ms}$.
   - The slow mode ($762\text{ ms}$) represents resensitisation from $A_2D$.
2. **Missing Intermediate Timescale:** There is **no eigenvalue anywhere near $33.4\text{ ms}$** in any of the three models. The model has only two macroscopic timescales: an unbinding/deactivation phase ($\sim 3.5 - 8.5\text{ ms}$) and a desensitisation recovery phase ($\sim 150 - 1000\text{ ms}$).
3. **Rank Deficiency:** When fitting 3 exponentials to a waveform that only contains 2 distinct macroscopic decay rates:
   - Either $\tau_1$ and $\tau_2$ merge to the exact same value ($\tau_1 \approx \tau_2 \approx 8.57\text{ ms}$ in `jw95`, $3.48\text{ ms}$ in `free_rates`), creating an infinite 1D continuum of ties where $a_1 + a_2 = \text{const}$ with zero change in residual.
   - Or an amplitude collapses to zero ($a_2 \to 0$), leaving the corresponding $\tau_2$ completely unidentifiable.
4. **Jacobian Spectrum:** While the true Barberis curve has a condition number of $39.0$ with five non-zero singular values, the model waveforms produce condition numbers ranging from $2.5 \times 10^5$ to $4.9 \times 10^8$ and $\infty$, with 1 or 2 singular values dropping to zero.

The optimiser did not fail; it accurately discovered that a 3-exponential model is structurally over-parameterised on a 5-state Markov chain.

---

## 5. Grid Sensitivity Check: Hybrid vs. Uniform Grid

As required by the brief, every model fit was evaluated on both the dataset's hybrid grid (`_T_BARBERIS`, 101 points) and an identical-count uniform grid ($10.1\text{ ms}$ spacing over $0 - 1000\text{ ms}$):

| Model & Weighting | Hybrid Grid $\tau_1$ (ms) | Uniform Grid $\tau_1$ (ms) | Ratio / Shift | Sensitivity Verdict |
| :--- | :---: | :---: | :---: | :--- |
| **Published Truth + 0.005 noise** | **$2.78\text{ ms}$** | **$1.60\text{ ms}$** | **$0.58\times$** | Uniform grid scatters $43\%$ due to single early sample |
| `kinetic_jw95` (as-is) | $8.57\text{ ms}$ | $7.03\text{ ms}$ | $0.82\times$ | Unstable; collapses with $a_1 \to 0$ |
| `kinetic_jw95` (equal $1/n$) | $8.83\text{ ms}$ | $7.99\text{ ms}$ | $0.90\times$ | Unstable; shifts slow modes |
| `kinetic_jw95` (fixed $\sigma$) | $8.58\text{ ms}$ | $7.03\text{ ms}$ | $0.82\times$ | Unstable; early tail unconstrained |
| `kinetic_jw95_free_rates` (as-is) | **$3.33\text{ ms}$** | **$3.19\text{ ms}$** | **$0.96\times$** | **Stable** ($4\%$ shift; fast mode has large footprint) |
| `kinetic_jw95_free_rates` (equal $1/n$) | **$3.48\text{ ms}$** | **$3.24\text{ ms}$** | **$0.93\times$** | **Stable** ($7\%$ shift; fast mode has large footprint) |
| `kinetic_jw95_free_rates` (fixed $\sigma$) | **$3.62\text{ ms}$** | **$3.51\text{ ms}$** | **$0.97\times$** | **Stable** ($3\%$ shift) |
| `kinetic_open_desens` (as-is) | $9.34\text{ ms}$ | $8.40\text{ ms}$ | $0.90\times$ | Unstable |
| `kinetic_open_desens` (equal $1/n$) | $5.45\text{ ms}$ | $4.43\text{ ms}$ | $0.81\times$ | Unstable |

### Key Takeaway:
For models where $k_{\text{off}}$ is freed (`kinetic_jw95_free_rates`), $\tau_{\text{fast}}$ is around $3.3 - 3.6\text{ ms}$ with large amplitude ($a \approx 0.65 - 0.70$), giving it sufficient presence across multiple $10\text{ ms}$ bins to remain relatively stable ($< 7\%$ shift). For fixed-rate models where the decay is fast and steep, the uniform grid severely under-samples the initial decay and drifts significantly. The hybrid grid remains strictly necessary for unbiased identification.

---

## 6. Shared Acceptance Criteria Checklist

- [x] **A function that, given a ReceptorModel, returns predicted ($\tau$, area) components for the Barberis protocol, by your assigned method:**  
  Implemented as pure function `predict_components_fit(model, *, pulse_ms=2.0, conc_mm=10.0)` in `src/circuitpharm/fitting/deact_components.py`.
- [x] **A likelihood term scoring those predictions against the table above at the published errors — not at a chosen sigma, and not at the reconstructed curve:**  
  Implemented as `deact_components_log_likelihood(taus, areas, *, deact_row=None)` in `src/circuitpharm/fitting/deact_components.py`, reading published values and errors from `circuitpharm.literature.DEACTIVATIONS[1]` and its note.
- [x] **A test that your predictor reproduces the published components when given a model built to have them, or, if that is not constructible, a stated reason and a weaker check that is honest about what it proves:**  
  Implemented in `tests/test_deact_components.py` (`test_recovers_published_components_from_published_curve` and `test_predict_components_fit_on_dataset_and_surrogate_model`).
- [x] **`RESULTS-BRIEF-7B.md` with the command, runtime, your predicted components at each of the three models' brief-6 optima, and a plain answer to the question above:**  
  Documented in this file.
- [x] **Both test tiers pass with PYTHONPATH set; counts stated separately and labelled:**  
  - Fast test suite (`pytest -q -m "not slow"`): **532 passed**, 40 warnings in 297.30 s (baseline on main was 524; +8 new tests).  
  - Full test suite (`pytest -q`): **604 passed**, 131 warnings in 991.57 s (16 min 31 s; baseline on main was 596; +8 new tests). Verified 100% green.
- [x] **Do not touch TRAIN, HOLDOUT, holdout_guard, the datasets, the sampling grids, provenance.py, or the sources table. Weighting (iv) is an additional objective; leave (i)–(iii) in place:**  
  None of those structures were touched or modified. Weighting (iv) is implemented purely as an additional objective.
- [x] **Do not delete or modify scripts/weighting_sensitivity.py's existing behaviour:**  
  `scripts/weighting_sensitivity.py` remains completely untouched.

---

## 7. Recommendation for the Planning Agent (Opus)

Per the discipline rules (§0 and line 68 of AGENT-BRIEF-7B):
- We do **not** unilaterally adopt Weighting (iv) as the project default.
- We report our findings to Opus and the user to decide how Weighting (iv) should be integrated into future model ranking.
