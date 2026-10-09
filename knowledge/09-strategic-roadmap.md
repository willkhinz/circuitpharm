# Strategic Research Manifesto: From Numerical Phenomenon to Rigorous Computational Neuropharmacology

**Target Directive for All Agents:**
> *Move from "a Markov model produces a surprising difference between phasic and tonic modulation" to "a rigorously tested model explains when that difference should occur, identifies the kinetic mechanisms responsible, and makes predictions that competing models cannot."*

---

## 1. Executive Summary & Vision

The core intuition developed in `knowledge/08-manuscript.md` (and underlying `src/circuitpharm/gabaa_kinetics.py`) is compelling:
**Phasic and tonic inhibition operate in fundamentally distinct agonist concentration regimes, and therefore possess radically different headroom for positive allosteric modulation.**

However, a model that simply displays a striking numerical result (e.g., the 184.6-fold asymptotic dynamic range under Jones & Westbrook 1995 parameters) remains vulnerable if:
1. The parameters are under-identified (3 calibration anchors cannot uniquely constrain 6 microscopic rates).
2. The dynamic range calculation relies on an unachievable infinite-affinity limit ($k_{\text{off}} \to 0^+$) or a hand-selected resting baseline.
3. The model does not separate the physical drivers of the effect (occupancy vs. kinetic redistribution vs. operating point).
4. Competing models are not formulated, benchmarked, or pitted against prospective experiments.

To transform this manuscript from an exploratory computational paper into a foundational contribution in computational neuropharmacology, every future agent working in this repository must align work with the **Six Strategic Advancements** detailed below.

---

## 2. The Six Prioritized Advancements

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                               THE STRATEGIC PYRAMID                                    │
│                                                                                        │
│   [4] Discriminating Experiment  ───> Falsifiable, prospective, hypothesis-testing     │
│   [3] Competing Models (A, B, C) ───> Scalar ceiling vs Kinetic vs Extended mechanism │
│   [2] Mechanistic Decomposition  ───> Occupancy vs Kinetic redistribution vs Headroom  │
│   [1] Empirical Data Fitting     ───> Likelihood, Bayesian MCMC, identifiability        │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### Advancement 1: Make the Mechanism Identifiable, Not Just the Output
**Priority:** Essential / Highest Priority
* **The Problem:** The current model shows *that* tonic modulation exhibits greater fold-change than phasic modulation, but does not isolate *why*.
* **Actionable Requirement:** Decompose the overall modulation response into three mathematically distinct, separable factors:
  1. **Agonist Occupancy:** The fraction of receptors occupying ligand-bound states ($R \leftrightarrow RG \leftrightarrow RG_2$) across steady-state vs. transient pulses.
  2. **Kinetic Redistribution:** How allosteric modulation alters flux between closed, open, and desensitized conformations ($\beta/\alpha$, $d/r$) holding occupancy fixed.
  3. **Operating Point / Saturation Proximity:** How close baseline open probability $P_{\text{open},0}$ already is to the gating ceiling $P_{\text{o,max}} = \beta/(\alpha+\beta)$, dictating theoretical headroom independently of channel kinetics.
* **Analysis:** Systematically vary each factor while holding the other two constant. Determine which factor is necessary and which is sufficient to generate the compartment divergence.

---

### Advancement 2: Fit to Empirical Data & Quantify Parameter Uncertainty
**Priority:** Essential / Highest Priority
* **The Problem:** The nominal Jones & Westbrook rate set is one point in an under-constrained parameter manifold. Reviewers will rightly object that alternative parameterizations could erase the phasic–tonic divergence.
* **Actionable Requirement:**
  1. Fit kinetic parameters to published patch-clamp data:
     * Concentration–response curves (EC₅₀, Hill slope).
     * Single-channel kinetics (mean open time, burst duration, open probability).
     * Deactivation and desensitization time courses (rapid agonist application/removal).
  2. Perform formal identifiability analysis (Profile Likelihood or Bayesian MCMC).
  3. Calculate posterior predictive distributions for all quantities of interest.
  4. Establish that the phasic–tonic headroom divergence holds across the entire data-consistent parameter ensemble, not merely at the nominal point.

---

### Advancement 3: Simulate Realistic, Time-Varying GABA Waveforms
**Priority:** High
* **The Problem:** Real synapses do not see instantaneous step pulses or permanent square waves; tonic extrasynaptic space experiences fluctuating spillover. Furthermore, steady-state $P_{\text{open}}$ ratios do not map 1:1 to integrated synaptic charge transfer $\int I(t)\,dt$.
* **Actionable Requirement:** Drive the kinetic scheme with realistic temporal waveforms:
  1. Synaptic transients with physiological rise ($\tau_{\text{rise}} \approx 0.1\text{–}0.5\text{ ms}$) and biexponential clearance ($\tau_{\text{fast}} \approx 1\text{ ms}, \tau_{\text{slow}} \approx 10\text{–}30\text{ ms}$).
  2. Repetitive pulse trains (e.g., 10 Hz, 50 Hz, 100 Hz) to assess accumulation of desensitization ($P_D$) and short-term synaptic plasticity under PAM exposure.
  3. Sustained extrasynaptic baselines ($[G] \approx 0.2\text{–}0.8\ \mu\text{M}$) modulated by activity-dependent spillover.
* **Metrics:** Report peak current ($I_{\text{peak}}$), total charge transfer ($Q = \int I\,dt$), deactivation kinetics ($\tau_{\text{decay}}$), and steady-state tonic current ($I_{\text{tonic}}$).

---

### Advancement 4: Formulate Competing Models & Predict Divergent Outcomes
**Priority:** Essential
* **The Problem:** Without explicit competitors, the paper cannot establish what makes its kinetic formulation uniquely necessary.
* **Actionable Requirement:** Implement and contrast three distinct model classes:
  * **Model A (Scalar Ceiling):** Phenomenological operational model where PAM efficacy is governed by an empirical Michaelis-Menten / Hill maximal fold-increase, blind to microscopic kinetics.
  * **Model B (Kinetic Allostery):** Microscopic Monod-Wyman-Changeux or Del Castillo-Katz scheme where modulation accelerates channel opening ($\beta$) or slows unbinding ($k_{\text{off}}$) with zero ad-hoc gain limits.
  * **Model C (Extended Kinetic Scheme):** Dual-pathway scheme incorporating subunit-dependent desensitization transitions or state-dependent modulator affinity (e.g., preferential PAM binding to closed vs. desensitized states; ref. *Nature Communications* desensitization pathway).
* **Evaluation:** Benchmark using BIC/AIC and out-of-sample cross-validation. Identify specific agonist/modulator concentration protocols where Models A, B, and C make diametrically opposed predictions.

---

### Advancement 5: Replace Single Headline Numbers with Robust Predictive Ensembles
**Priority:** High
* **The Problem:** The headline "$184.6\times$ dynamic range" is an asymptotic limit as $k_{\text{off}} \to 0^+$ at a specific baseline ($[G] = 0.4\ \mu\text{M}$), which can distract from the biological claim.
* **Actionable Requirement:** Clearly separate three distinct claims:
  1. **Theoretical Maximum Headroom:** The asymptotic limit of the Markov model ($P_{\text{open},\infty} / P_{\text{open},0}$).
  2. **Pharmacologically Reachable Gain:** The gain attained under clinically/experimentally achievable concentrations of specific PAMs (e.g., $s_{\text{max}} \le 2.5\times$ for diazepam/alogabat).
  3. **Physiological Current & Charge Impact:** The resulting change in inhibitory charge transfer under realistic neuronal driving forces.
* Propagate full parameter uncertainty to show credible intervals (95% CI) rather than point estimates.

---

### Advancement 6: Design the Optimal Discriminating Experiment
**Priority:** Potentially Transformative
* **The Problem:** Descriptive models explain the past; transformative models prescribe future experiments.
* **Actionable Requirement:** Use Bayesian Optimal Experimental Design (OED) to calculate the precise experimental protocol that maximizes information gain and discriminates between Models A, B, and C:
  * Concentration points ($[G]$ concentrations and PAM titration steps).
  * Fast-application pulse lengths (transient 1 ms vs. stepped 100 ms vs. sustained).
  * Falsification criteria: specific, quantitatively bounded outcomes that would definitively refute our kinetic hypothesis.

---

## 3. The Four Core Results of the Ambitious Manuscript

1. **Result 1 (Data-Constrained Kinetic Family):**
   The kinetic scheme reproduces multiple independent single-channel and macroscopic datasets; parameter uncertainty is bounded and quantified via MCMC.
2. **Result 2 (Mechanistic Decomposition):**
   A rigorous breakdown demonstrating that operating point saturation ($P_{\text{o,max}} - P_0$) and desensitization trapping govern the phasic ceiling, while low baseline occupancy dictates tonic expansiveness.
3. **Result 3 (Generalizable & Bounded Predictions):**
   A complete phase map across $[G]$, receptor subtype kinetics, and modulator affinity/efficacy showing exactly where the compartment divergence thrives and where it collapses.
4. **Result 4 (Prospective Experimental Test):**
   A falsifiable protocol predicting out-of-sample current waveforms that distinguish kinetic allostery from phenomenological models.

---

## 4. Architectural Guardrail: Anti-Complexity Rule

> **Never add biological complexity solely for sophistication.**
> Do not introduce more receptor states, unanchored conductance parameters, or larger circuit models unless:
> 1. It explains empirical data that the simpler model structurally fails to explain, OR
> 2. It generates a distinct, testable, and falsifiable prediction.

All agents must preserve this roadmap across all future implementation sprints.
