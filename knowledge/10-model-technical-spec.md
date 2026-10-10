# Technical Specification: Next-Generation GABA-A Kinetic & Neuropharmacology Model

**Scope:** Architectural and algorithmic design for the computational model in `src/circuitpharm/`.  
**Consensus Date:** 2026-10-08 (via structured design tree alignment).  
**Revised:** 2026-10-09, after roadmap phases P0–P7 were built.

> **This document used to describe things that do not exist, and that is a defect of the
> same kind as a mis-cited number.** A specification is read as a description of the
> system; one that overstates it misleads exactly as a wrong figure does, and is harder to
> catch because nothing recomputes it. Every section below has been checked against the
> implementation. Where the design was not built, it now says so, in a **NOT BUILT** note
> that names what exists instead — rather than being quietly deleted, because the gap is
> the information.
>
> The four that mattered most, all now corrected in place: Model C's modulator does **not**
> act through a thermodynamic cycle; **no dataset here is digitised**; profile likelihood
> over six microscopic rates is **VOID**, not a result; and the attribution in
> `decomposition.py` is **not** a Shapley value.

---

## 1. Architectural Overview & Modular Layout

To move from an isolated 5-state Markov calculation to a rigorously tested, data-constrained, multi-model neuropharmacology platform without breaking existing tests (277 passing at the time of writing; 464 passing and 19 skipped after P0-P7), the codebase is organized into clean modular subpackages under `src/circuitpharm/`:

```
src/circuitpharm/
├── models/
│   ├── base.py                 # Abstract ReceptorModel protocol & base interface
│   ├── operational.py          # Model A: Phenomenological Operational Scalar Ceiling
│   ├── kinetic_jw95.py         # Model B: 5-state Jones & Westbrook Kinetic Allostery
│   └── extended_desens.py      # Model C: Dual-pathway desensitization & state-dependent PAM
├── fitting/
│   ├── data.py                 # Benchmark datasets; NONE digitised -- see §4
│   ├── reparam.py              # The 3 combinations an equilibrium CRC determines
│   ├── likelihood.py           # Per-observable likelihood, sigma profiled out
│   ├── posterior.py            # Posterior over the identifiable combinations
│   ├── comparison.py           # Blocked-CV model comparison (P5)
│   ├── identifiability.py      # Profile likelihood & Fisher Information Matrix
│   └── mcmc.py                 # Pure NumPy/SciPy affine-invariant ensemble MCMC
├── protocols/
│   ├── waveforms.py            # Transient, pulse-train & ambient spillover generators
│   ├── decomposition.py        # Log-additive & Shapley mechanistic attribution
│   └── oed.py                  # Optimal experimental design & falsification protocol
└── gabaa_kinetics.py           # Preserved backwards-compatible façade wrapping Model B
```

---

## 2. Model Hierarchy (Models A, B, C)

All models implement the polymorphic `ReceptorModel` protocol (`models/base.py`):
```python
class ReceptorModel(Protocol):
    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float: ...
    def simulate_waveform(self, t: np.ndarray, gaba_t: np.ndarray, pam_factor: float = 1.0) -> WaveformResult: ...
    def dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0) -> np.ndarray: ...
```

### Model A: Operational Scalar Ceiling (`models/operational.py`)
* Phenomenological operational model of allosterism (Black & Leff / modern operational Hill).
* Agonist response: $E([G]) = E_{\max} \frac{[G]^n}{[G]^n + EC_{50}^n}$.
* Modulator acts via empirical potency-shift factor $s \le s_{\max}$ and scaling factor $\alpha_{\max}$, imposing an extrinsic efficacy ceiling blind to microscopic state flux.

### Model B: Kinetic Allostery Scheme (`models/kinetic_jw95.py`)
* 5-state Del Castillo–Katz / Jones–Westbrook (1995) Markov scheme:
  $$R \underset{k_{\text{off}}}{\overset{2 k_{\text{on}} [G]}{\rightleftharpoons}} AR \underset{2 k_{\text{off}}}{\overset{k_{\text{on}} [G]}{\rightleftharpoons}} A_2R \underset{\alpha}{\overset{\beta}{\rightleftharpoons}} A_2O$$
  with desensitization from the shut doubly-bound state: $A_2R \underset{r}{\overset{d}{\rightleftharpoons}} A_2D$.
* PAM acts via accelerated binding / slowed unbinding ($k_{\text{off}} \leftarrow k_{\text{off}} / s$) and/or gating acceleration ($\beta \leftarrow \beta \cdot g$), with zero ad-hoc gain caps.

### Model C: Extended Desensitization & State-Dependent Affinity (`models/extended_desens.py`)
* Motivated by the *Nature Communications* GABAA desensitization pathway.
* Dual desensitization branches: fast ($D_{\text{fast}}$) and slow ($D_{\text{slow}}$) conformational sinks.
* **NOT BUILT: "state-dependent affinity via a thermodynamic cycle".** The implementation
  scales `d_fast` by a single linear `pam_desens_factor`. That is a one-parameter
  perturbation of one rate, not a cycle: a thermodynamic cycle would require the modulator's
  binding constants to each of $R$, $O$ and $D$ and would have to satisfy detailed balance
  around the $R \to O \to D \to R$ loop, which scaling one rate does not. The design is
  still worth building — it is the structural difference that would make Model C testable
  against Model B — but it is not what the code does, and the code's `pam_desens_factor` is
  additionally unconstrained by any dataset here (P5: its posterior fails every
  convergence diagnostic, `knowledge/13-preregistration.md` §2).

---

## 3. Mechanistic Decomposition Engine (`protocols/decomposition.py`)

Isolates the three physical contributors to the phasic vs. tonic modulation difference:
1. **Agonist Occupancy ($\Delta \ln \text{Occ}$):** Fraction of receptors bound by GABA.
2. **Kinetic Redistribution ($\Delta \ln \text{Kinet}$):** Conformational equilibrium between bound-shut, open, and desensitized states ($\beta/\alpha$, $d/r$).
3. **Operating Point Headroom ($\Delta \ln \text{Headroom}$):** Proximity to the maximal open probability ceiling $P_{\text{o,max}} = \beta/(\alpha+\beta)$.

**Algorithm:**
* Exact factorial grid decomposing total logarithmic gain:
  $$\ln G = \Delta \ln(\text{Occupancy}) + \Delta \ln(\text{Kinetic Redistribution}) + \Delta \ln(\text{Operating Point Headroom})$$
* **NOT A SHAPLEY VALUE, and the fields are no longer named as one.** The computation is
  $|\Delta_i| / \sum_j |\Delta_j|$: a normalised absolute-log share. It is not
  permutation-averaged, it discards sign, and it does not sum to the total log gain. The
  fields are `log_share_pct_*` (roadmap P4-3). What *is* exact is the identity above —
  `FactorialDecompositionResult.check()` verifies it to `rtol=1e-10` and raises rather than
  returning a result that does not satisfy it. The earlier definition of
  `headroom_factor = (P_o,max − P_o,base)/(P_o,max − P_o,pam)` did **not** satisfy it, and
  the docstring asserting the identity was wrong for as long as it stood.

---

## 4. Parameter Estimation, Identifiability & Uncertainty (`fitting/`)

* **Datasets (`fitting/data.py`): NONE OF THEM ARE DIGITISED.** Three categories are
  distinguished in code, because collapsing the middle one either way misleads:
  * `digitised` — individual points read off a published figure. **Count: zero.**
  * `parametric` — generated from a paper's own *published fitted parameters*. `JAHN1997_PEAK_CRC`
    is this: a Hill curve from Jahn et al. 1997's EC₅₀ = 11.6 ± 0.9 µM and n_H = 2.2 ± 0.4.
    The parameters are real measurements and attributable; the nine points are the paper's
    model of its data. It carries `synthetic=True`, so it taints any fit it touches.
  * `synthetic` — invented, for exercising machinery. `DOSE_RESPONSE_BENCHMARK` and
    `DEACTIVATION_BENCHMARK` are this. The latter's dominant 15 ms component was for three
    commits attributed to a paper that measured 76.1 ms.
  * **NOT BUILT: single-channel open probability and burst durations, and a real rapid-perfusion
    deactivation trace.** `fitting.data.MISSING_DATASETS` names all three gaps and why they are
    blocked (full texts unreachable from this environment). The mean open time is the
    highest-value acquisition in the project: it is what would turn
    `gabaa_kinetics.FIT_FIXED_ALPHA` from a convention into a measurement.
* **Identifiability (`fitting/identifiability.py`, `fitting/reparam.py`):**
  * **An equilibrium concentration-response determines THREE combinations, not six rates.**
    `K_d = k_off/k_on`, `E = β/α` and `D = d/r`; the six rates enter only through them, so
    the likelihood is invariant under a three-parameter scaling group. Proved, and verified
    numerically to 2 × 10⁻¹² along each generator.
  * **Profile likelihood over the six microscopic rates returns `Tier.VOID`** — it profiles
    directions the objective cannot see, and its intervals are prior widths. The machinery
    stays callable behind `ProvisionalResultWarning` rather than being deleted.
  * Profiles over the three identifiable combinations are real and are classified
    IDENTIFIABLE / PRACTICALLY NON-IDENTIFIABLE / STRUCTURALLY NON-IDENTIFIABLE. At 0.2%
    noise `E` and `D` each span more than two decades.
  * **NOT BUILT, and provably not buildable from this data: a Fisher information matrix.**
    What `cost_hessian_and_spectrum` returns is the Hessian of a chi-squared at whatever
    parameters it was handed, with its FULL spectrum including negatives. Filtering to
    positive eigenvalues is how the previous version turned "not at a minimum, and
    rank-deficient" into a finite, plausible condition number.
* **Affine-Invariant MCMC Ensemble (`fitting/mcmc.py`, `fitting/posterior.py`):**
  * Lightweight pure NumPy/SciPy Goodman & Weare ensemble sampler (no heavy PPL dependencies).
    Sampling is in log₁₀, where a flat prior IS log-uniform — the previous version returned
    0.0 inside a box on the rates and called itself log-uniform, which for parameters
    spanning decades is a materially different prior.
  * **The stretch scale is a MIXTURE (2, 20, 200), not the canonical 2.** A single `a = 2`
    cannot propose a contraction below one half, which froze 2 of 24 walkers 2.5 decades
    from the mode for an entire run while split-R̂ read 62. Per-walker acceptance is now a
    mandatory diagnostic.
  * **Convergence diagnostics gate the result**: integrated autocorrelation time, chain
    length in multiples of it, split-R̂, acceptance fraction and frozen walkers. A run that
    fails any of them returns VOID quantities naming the failure — reading one raises.
  * **NOT BUILT: validation that the headroom divergence holds across the 95% credible
    interval.** `dynamic_range_posterior` propagates draws through all three tiers and
    reports median + [2.5, 97.5] per tier, which is the input to such a claim; the claim
    itself has not been made, because no posterior here rests on a digitised dataset.

---

## 5. Realistic Waveform & Electrophysiology Engine (`protocols/waveforms.py`)

Provides standardized physiological stimulation protocols:
1. **Biexponential Synaptic Pulse:** Variable rise ($\tau_{\text{rise}} = 0.1\text{–}0.5\text{ ms}$) and biexponential clearance ($\tau_{\text{fast}} = 1\text{ ms}, \tau_{\text{slow}} = 10\text{–}30\text{ ms}$).
2. **High-Frequency Pulse Trains:** 10 Hz, 50 Hz, 100 Hz bursts to quantify frequency-dependent desensitization accumulation ($P_D$) and paired-pulse depression.
3. **Activity-Dependent Spillover:** Baseline ambient GABA with transient spillover elevations.

**Standard Output Metrics:**
* Peak open probability and peak current: $P_{\text{peak}}$, $I_{\text{peak}} = g_{\max} P_{\text{peak}} (V - E_{\text{Cl}})$.
* Integrated inhibitory charge: $Q = \int I(t)\,dt$.
* Biexponential decay time constants: $\tau_{\text{fast}}, \tau_{\text{slow}}, \tau_{\text{weighted}}$.
* Paired-pulse ratio $I_2/I_1$, measured **from the open probability at pulse onset, not
  from zero**. At 100 Hz the channel has not closed between pulses — deactivation τ is
  3.3 ms against a 10 ms period — so a peak measured from zero is mostly residual from the
  preceding pulse. Measured from zero, this scheme reports paired-pulse *facilitation*
  deepening with frequency (0.863 → 0.900 → 0.958) while in fact depressing threefold
  (0.857 → 0.552 → 0.282). Both numbers are returned; `ppr_2_over_1_from_zero` is the
  artefact, kept so it stays demonstrable.

---

## 6. Three-Tiered Dynamic Range Metric (`DynamicRangeEvaluation`)

Replaces single headline point estimates with a structured metric dataclass:
1. **Theoretical Asymptotic Headroom:** $P_{\text{open},\infty} / P_0$ as $k_{\text{off}} \to 0^+$ at stationary ambient $[G]$.
2. **Pharmacologically Reachable Gain:** Gain evaluated at ligand-specific operational efficacy ceilings ($s_{\max} = 1.25, 1.50, 2.50$).
3. **Physiological Charge Transfer Ratio:** Realized ratio $\Delta Q_{\text{tonic}} / \Delta Q_{\text{phasic}}$ under realistic driving forces and time-varying waveforms.
* `dynamic_range_posterior` reports median and [2.5%, 97.5%] intervals per tier; the
  point-estimate entry point `evaluate_dynamic_range` is kept unchanged so existing results
  stay reproducible, and there is deliberately no way to call the interval version and get a
  number without an interval. The collapse boundary is returned as a **curve**
  (`collapse_scan`), not the boolean `headroom_collapses` alone, which threw away everything
  the scan computed. Note which tier is worst determined: the asymptotic headroom carries
  the uncertainty in both ratios including the `D` no anchor constrains, and spans 2.99-fold
  where the reachable gains span 1.05-fold.

---

## 7. Model Discrimination & Optimal Experimental Design (`protocols/oed.py`)

* **Model Selection (`fitting/comparison.py`).** Each model is fitted independently to its
  own MLE first — an information criterion evaluated anywhere else ranks starting guesses.
  `k` counts **estimated identifiable parameters plus estimated noise scales**, not
  `len(param_names)`. AIC, BIC and **AICc** (n is 9–30 here, so the second-order correction
  is not a refinement). Cross-validation is by **contiguous concentration block**: random
  point-wise folds leak across a smooth curve and flatter the flexible model. A CV winner is
  named only if it clears both the paired per-fold standard error and a 2-nat floor;
  otherwise the tie goes to the model with fewest estimated parameters. Akaike weights are
  reported **only if the CV ranking agrees**, and are VOID otherwise.
  `oed.evaluate_model_fit` stays callable and VOID.
* **Optimal Experimental Design (`protocols/design.py`).**
  * **NOT BUILT: Jensen-Shannon divergence.** The score is a posterior-predictive
    separation, `E|y_A − y_B| / sqrt(var(y_A) + var(y_B) + σ_meas²)`, which falls when a
    posterior widens — the property the previous `|Δy|/σ_noise` lacked by construction. The
    search ranges over agonist, modulator **and the observable itself**, because what to
    measure is part of the design.
  * The multiple-comparison cost of having searched the protocol space is **carried, not
    mentioned**: the reported score is corrected for the number of protocols tried.
  * **NOT BUILT, and the finding is that it cannot be: "non-overlapping falsification
    boundaries".** At conventional measurement noise the best of 54 protocols scores 1.71
    raw and **0.00** corrected, and the two models' predictive intervals overlap over 60% of
    their width. The actionable requirement is a ~4× reduction in σ_meas, not a different
    concentration. See `knowledge/13-preregistration.md`.
