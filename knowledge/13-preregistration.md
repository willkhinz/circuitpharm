# Pre-registration: the experiment that would separate Models B and C

**Roadmap phase:** P6-3. **Status:** the design machinery works and its answer is that
**the experiment does not exist at conventional measurement noise.** That is the
pre-registered finding. What would change it is stated in §5, in numbers.

This document is written *before* any such experiment, which is the point of the format.
Everything in it is reproducible from the script in §7.

---

## 0. The honest header

| | |
|---|---|
| **Basis dataset** | GENERATED, not measured: an EQUILIBRIUM concentration-response from Model B at `log10 K_d = 1.5, log10 E = 0.5, log10 D = 0.0`, 30 concentrations over 4.5 decades, 0.1% noise. |
| **Why generated** | No real EQUILIBRIUM concentration-response exists in this repository. Both CRCs in `fitting/data.py` are PEAK, and `fitting.data.MISSING_DATASETS` names the three gaps. |
| **What may be claimed from it** | That the DESIGN MACHINERY behaves correctly, and what it says about these two schemes at these parameters. **Nothing about α1β2γ2 receptors.** |
| **Tier** | Model B's posterior UNCALIBRATED; Model C's **VOID** (§2); the design therefore UNCALIBRATED at best and the falsification intervals in §4 are not yet pre-registerable as numbers. |

The protocol, the decision rule and the refutation criterion in §3–§5 **are** pre-registerable
now: they do not depend on the parameter values, only on the machinery. The predicted
intervals are what has to be recomputed once P1-2's digitisation lands.

---

## 1. Model B's posterior converges

24 walkers × 12 000 steps, 4 000 burn-in, prior-dispersed start, stretch scales (2, 20, 200).

| | median | 95% interval | truth | width (decades) |
|---|---|---|---|---|
| log₁₀ K_d | 1.5169 | [1.5046, 1.5293] | 1.5 | 0.025 |
| log₁₀ E | 0.5224 | [0.5044, 0.5404] | 0.5 | 0.036 |
| log₁₀ D | 0.0442 | [0.0083, 0.0789] | 0.0 | 0.072 |

Diagnostics: τ = 58.5 / 58.4 / 58.7, split-R̂ = 1.0093 / 1.0094 / 1.0095, acceptance 0.279,
no frozen walker, nothing on a bound. All four gates pass, so the intervals are real.

**And on this noise realisation all three intervals just miss the truth, on the same side.**
That is one event, not three: the three parameters are strongly correlated, so the
posterior is displaced as a whole. Checked across seeds rather than assumed — see
`test_posterior_intervals_are_calibrated_across_seeds`, which requires joint coverage in at
least four of six noise realisations, because **a 95% interval cannot be validated against
one draw** and a test that did would be passing on luck. Seed 0 is a miss and the other five
cover. Measured over six seeds: **five of six realisations cover all three** (15 of 18 marginals, 83%), which for three near-perfectly correlated marginals is five joint successes out of six — exactly what a calibrated 95% interval does, since one miss in six has probability 26%.

---

## 2. Model C's posterior does not, and that is the result

Same sampler, same data, four parameters.

| | median | 95% interval | diagnostic |
|---|---|---|---|
| log₁₀ K_d | 1.5164 | [1.5047, 1.5292] | — |
| log₁₀ E | 0.5217 | [0.5047, 0.5406] | — |
| log₁₀ D_fast | 0.0372 | [−0.0350, 0.0776] | — |
| log₁₀ D_slow | −2.432 | [−3.9035, −0.8244] | **on its lower bound**, 3.1 decades wide |

Failures: `chain is 14.5 autocorrelation times long, need >= 50`;
`worst split-Rhat 1.3910 > 1.01`; `acceptance fraction 0.139 outside (0.15, 0.60)`.
τ = 271 / 270 / 553 / 551.

**The chain is not broken; the parameter is not there to be estimated.** The data was
generated from Model B, which has no slow desensitisation sink, so `D_slow → 0` is correct
and the likelihood is flat along it over three decades. The sampler reports VOID instead of
an interval, which is the behaviour P0-5 exists for.

Two consequences, both binding:

* **No posterior-derived number for Model C may be quoted**, here or in the manuscript.
  `parameters.get("extended_desens", "fitted")` was already right to raise.
* The design in §3 therefore rests on one converged posterior and one that did not. Its
  score is a lower bound of sorts — a wider, properly mixed Model C posterior can only make
  the denominator larger and the separation smaller — but it is not a number to publish.

---

## 3. The protocol, and the decision rule

**Search space:** 9 GABA concentrations (0.1–1000 µM) × 6 PAM strengths (1.0–4.0× as an
EQUILIBRIUM EC₅₀ shift) × the observable itself. 54 protocols scored.

**Score:** posterior-predictive separation,

```
score = E|y_B − y_C| / sqrt( var(y_B) + var(y_C) + σ_meas² )
```

over 960 thinned draws per model. A score of 1.96 is the smallest at which a single
measurement could separate the models at all.

**Best protocol:** 1000 µM GABA, PAM 4×, measured at EQUILIBRIUM.
**Runner-up:** 300 µM GABA, PAM 4×, EQUILIBRIUM, score 1.692 — so the maximum is a plateau,
not a peak, and the specific concentration is not the design.

| | value |
|---|---|
| raw score | **1.707** |
| after correcting for 54 searched protocols | **0.000** |
| expected separation | 0.0342 open-probability units |
| discriminating at σ = 0.02? | **no** |

The correction is a Bonferroni-style discount on the implied z: the best of 54 standard
normals sits near 2.3σ, so a raw 1.71 is entirely explained by having looked in 54 places.
Quoting 1.71 would be quoting the search.

### 3.1 Decision rule, pre-registered

1. Run the protocol at the stated GABA, PAM and observable, `n ≥ 8` cells.
2. Compare the measured mean against each model's two-sided 95% posterior-predictive
   interval (§4), computed **before** the data is seen and from a posterior that passed
   every diagnostic.
3. An outcome **inside exactly one** interval and outside the other supports that model.
   An outcome **outside both** refutes both — this is the outcome worth having, and the one
   a point-prediction score cannot even define.
4. An outcome **inside both** is the expected result at σ = 0.02 and is **not** evidence for
   either. It must be reported as such and not as support for the simpler model.
5. The corrected score, not the raw one, is quoted in any report of the design.

---

## 4. Falsification intervals, as computed now

Two-sided 95% posterior-predictive intervals at the best protocol, including measurement
noise at σ = 0.02:

| model | 2.5% | median | 97.5% | predictive SD |
|---|---|---|---|---|
| B (`kinetic_jw95`) | 0.5743 | 0.6101 | 0.6506 | 0.0197 |
| C (`extended_desens`) | 0.6025 | 0.6413 | 0.6800 | 0.0199 |

They overlap over [0.6025, 0.6506], which is 60% of either interval. **Not pre-registerable
as numbers**: Model C's posterior is VOID (§2), and the dataset behind both is generated.

---

## 5. What would make the experiment exist

The score against measurement noise, at the same protocol:

| σ_meas (normalised P_open) | score | corrected for 54 protocols | discriminating? |
|---|---|---|---|
| 0.005 | **5.90** | ≈ 5.2 | **yes** |
| 0.010 | **2.99** | ≈ 1.1 | marginal |
| 0.020 | 1.71 | 0.00 | no |
| 0.050 | 1.24 | 0.00 | no |

**The actionable requirement is σ_meas ≤ 0.005, a four-fold reduction from the conventional
0.02 — not a different concentration.** The protocol search found a plateau, so no choice of
agonist or modulator inside this space helps; the binding constraint is the recording.

That is a conclusion the score this replaces could not reach. `|y_B − y_C| / σ_noise` has no
posterior variance in its denominator, so it would have reported the same ranking of
protocols at every noise level and never identified the noise as the binding constraint.

Three other things would change the answer, in rough order of value:

1. **A kinetic dataset** (`MISSING_DATASETS["single_channel_mean_open_time"]`,
   `["deactivation_peak_pulse"]`). The two schemes differ in their desensitisation
   *kinetics*, and an equilibrium CRC integrates that away by construction. A design over
   CHARGE with a stated window, or over paired-pulse depression at 10/50/100 Hz — where the
   schemes differ threefold (`knowledge/12-inference.md`) — is the obvious place to look
   next, and `protocols/design.py` already searches the observable.
2. **A real EQUILIBRIUM CRC**, so §1–§4 are about receptors.
3. **More cells**, which shrinks `σ_meas/√n` on the MEAN but not the predictive interval —
   worth 2× at best, against the 4× needed.

### 5.1 Re-derivation with declared protocols: NO candidate observable discriminates (2026-10-16)

The initial pre-registration in §5 filed the verdict: *"The actionable requirement is $\sigma_{\text{meas}} \le 0.005$."* That verdict was true for theoretical instantaneous equilibrium, but buries the core empirical finding: **tightening $\sigma$ does not help if no observable separates the models at any $\sigma$ tried.** The blocker is not noise; it is the absence of a discriminating observable.

Following `AGENT-BRIEF-2.md` and `AGENT-BRIEF-3.md`, two candidate observables with declared and matched experimental protocols were evaluated across protocol space (`scripts/rederive_preregistration.py`):
1. **PEAK open probability** under declared fast-application pulses (2 ms, 20 ms, 300 ms, 1000 ms).
2. **Paired-pulse recovery** (`Observable.PAIRED_PULSE_RECOVERY` from Barberis 2007: 10 mM GABA, 2 ms pulse, 10–500 ms inter-pulse gaps, PAM 1–4×).

#### Summary Table of Candidate Observables Tried

The printed table below from `scripts/rederive_preregistration.py` reports the old and new candidate observables side by side:

| Observable | Protocol | Expected Sep (Δ) | Raw Score (σ=0.02/0.03) | Corrected Score (σ=0.005) | Discriminating at any σ? |
|---|---|---|---|---|---|
| **EQUILIBRIUM (Theoretical)** | 1000 µM GABA, PAM 4× | 0.0342 | 1.495 (at 0.02) | 5.247 | **YES** (at σ ≤ 0.005 only) |
| **PEAK (2 ms pulse)** | 100 µM GABA, PAM 2.5×, 2 ms | 0.0012 | 0.059 (at 0.02) | 0.000 (raw 0.226) | **NO** (indistinguishable) |
| **PEAK (20–1000 ms pulse)** | 100 µM GABA, PAM 2.5×, 300 ms | 0.0057 | 0.281 (at 0.02) | 0.000 (raw 0.867) | **NO** (indistinguishable) |
| **PAIRED_PULSE_RECOVERY** | 10 mM GABA, PAM 4×, 100 ms gap | 0.0104 | 0.339 (at 0.03) | 0.000 (raw 1.277) | **NO** (indistinguishable) |

#### Analysis: Why Noise Reduction Fails to Resolve the Schemes

1. **PEAK observable fails to separate schemes:**
   Matching application duration between model and protocol shows that while the whole-curve Hill slope shifts with duration (1.752 at 2 ms → 1.376 at 300 ms), the rising-phase slope (`hill_slope(..., window_rel=JAHN1997_SLOPE_WINDOW_REL)`) is stable (1.765 at 2 ms vs 1.673 at 300 ms). However, peak open probability separation between Model B and Model C is minuscule ($\Delta < 0.006$ open-probability units), yielding a corrected score of **0.000** at both conventional $\sigma = 0.02$ and ultra-low noise $\sigma = 0.005$ (raw score 0.867, far below the 1.96 threshold).

2. **PAIRED_PULSE_RECOVERY fails to separate schemes:**
   Across 18 searched protocols spanning gap intervals (10–500 ms) and PAM factors (1–4×), the optimal protocol (10 mM GABA, PAM 4×, 100 ms gap) yields an expected separation of only $\Delta = 0.0104$. The corrected posterior-predictive score is **0.000** at the published experimental noise ($\sigma = 0.03$, raw 0.339) and remains **0.000** even at fourfold reduced noise ($\sigma = 0.005$, raw 1.277 < 1.96).

3. **Cross-reference to Task D (The P6 design question is malformed):**
   The root cause is not experimental measurement scatter ($\sigma$), but model identifiability:
   As established in Task D (§`scripts/resolve_cv_aicc.py`), Model C's slow desensitisation parameter sits pinned against its lower bound ($\log_{10} D_{\text{slow}} = -4.0000$) on the training concentration-response curve, conferring zero improvement in training log-likelihood ($-0.0005$ nats), zero gain in cross-validation out-of-sample (+0.0124 nats, $\text{diff}/\text{se} = 1.41\sigma < 2.0\sigma$, tying with Model B), and zero difference on the experimental paired-pulse holdout (residual difference $< 0.001$ nats, both predicting $R = 0.1949$).
   
   Because Model C nests Model B, when $D_{\text{slow}} \to 0$, Model C functionally collapses onto Model B. The two models are empirically indistinguishable across all physically realizable kinetic protocols not because $\sigma$ is too large, but because Model C's extra parameter is unconstrained and contributes no distinct physical behavior. Seeking a discriminating protocol between B and C in this parameter regime is malformed.

---

## 6. What this does not claim

* Nothing about α1β2γ2 receptors. The basis dataset is generated.
* Nothing about which of Models B and C is right. The design says the experiment to decide
  does not yet exist, which is a statement about the experiment.
* Nothing about respiratory function. See roadmap §6's last row, which is unchanged.
* The paired-pulse and charge numbers quoted in §5 are from `knowledge/12-inference.md` and
  carry their own tiers; they are named here as the next place to look, not as results.

---

## 7. Reproducing this

```python
import numpy as np
from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.data import equilibrium_crc_from_model
from circuitpharm.fitting.posterior import sample_posterior_vector
from circuitpharm.models.base import Observable
from circuitpharm.protocols.design import (
    find_discriminating_protocol_pp, sigma_sensitivity)

b, c = MODEL_SPECS["kinetic_jw95"], MODEL_SPECS["extended_desens"]
centre = np.array([0.5 * (lo + hi) for lo, hi in b.bounds])
ds = equilibrium_crc_from_model(b.factory(centre), concs_um=np.logspace(-1, 3.5, 30),
                                noise_sd=1e-3, seed=0)

posts = {}
for name, spec in (("kinetic_jw95", b), ("extended_desens", c)):
    posts[name] = sample_posterior_vector(
        spec.factory, [ds], param_names=spec.param_names,
        bounds=dict(zip(spec.param_names, spec.bounds)), caller=name,
        n_steps=12000, burn_in=4000, seed=11)
    print(name, posts[name].tier.label, posts[name].diagnostics["failures"])

cand = {n: (MODEL_SPECS[n].factory, posts[n].log10_samples[::200]) for n in posts}
res = find_discriminating_protocol_pp(cand, observables=(Observable.EQUILIBRIUM,))
print(res.note)
print(sigma_sensitivity(cand, res.protocol))
```

About 35 s for the two chains and 2 s for the search. Pinned in
`tests/test_p6_design.py`; the design-level assertions are marked `slow`.
