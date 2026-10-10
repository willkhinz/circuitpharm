# Compartment-Dependent Allosteric Headroom in GABA-A Receptors

## A receptor-kinetic analysis of phasic and tonic inhibition, with conditional circuit implications

**William Hinz**¹
¹*Department of Biochemistry & Computational Neuropharmacology*
Correspondence: hinzwilliam52@gmail.com
Repository: https://github.com/willkhinz/circuitpharm — tag **`manuscript-v3`**

> **Scope and status of claims.** This manuscript reports three kinds of result, and they do
> not carry equal weight. Readers should not transfer confidence from the first to the third.
>
> 1. **Receptor-kinetic results (§2, §4.1, §5.2).** Derived from an explicitly stated five-state
>    Markov scheme. Fully reproducible, internally checked, and conditional only on the
>    parameterisation and topology, both of which are given in full. This is the defensible core.
> 2. **Compound ordering (§4.4, §4.6).** Robust *within the specified model*: the ranking
>    survives adverse draws of every uncertain parameter simultaneously, and survives replacing
>    the neuron model entirely. It is **not** biologically validated — see the explicit
>    distinction in §4.4.4.
> 3. **Anatomical and respiratory interpretation (§4.2, §4.3).** The weakest component. The
>    preBötzinger α5 fraction and the α5 extrasynaptic fraction — the two numbers on which the
>    respiratory-sparing argument and the extrasynaptic-headroom argument respectively rest —
>    are both **UNSOURCED** (§4.2, §7). **Nothing here is evidence that α5-selective compounds
>    preserve respiration or possess a safety advantage.**
>
> Two further limits apply throughout. The headline **123.3×** figure is an *asymptotic receptor
> open-probability dynamic range* under stated assumptions — not a predicted fold change in
> current, in inhibition, or in any clinical effect (§2.3). And the non-selective benzodiazepine
> that serves as the denominator of every selectivity ratio reported here **changes sign**
> between the two neuron substrates tested (§4.6).
>
> **What the kinetic parameters rest on, stated plainly.** The five-state scheme is fitted to
> **three macroscopic anchors and no digitised dataset.** The anchors (peak EC₅₀ 20 µM,
> P_o,max 0.750, IPSC decay τ 15 ms) are recalled literature ranges, not points read off a
> named figure; `α` is held at a **declared convention** (0.30 ms⁻¹, a 3.33 ms mean open time)
> rather than measured, because without it the system has four unknowns against three
> residuals and its solution moved between SciPy versions; and `d` and `r` are constrained by
> **none** of the three anchors, so every asymptotic number here is conditional on
> `d/r = 25`. The one genuinely sourced value for this preparation's peak EC₅₀ is
> **11.6 ± 0.9 µM** (Jahn et al. 1997), not the 20 µM used as the anchor.
>
> Three consequences, each measured rather than asserted:
>
> * **An equilibrium concentration-response determines three numbers, not six.** `K_d`, `E` and
>   `D` are the only combinations it can see; the six microscopic rates enter through them
>   alone, so any per-rate confidence interval from such a fit is a statement about the prior
>   box. Proved and verified to 2 × 10⁻¹² in `knowledge/11-identifiability.md`.
> * **At realistic noise, two of those three collapse.** At 0.2% noise on a normalised curve
>   `E` and `D` each span more than two decades, and more precision on the same observable
>   will not fix it — two parameter sets a hundredfold apart give curves separated by less
>   than the noise. A *different* observable, one with a timescale in it, is required.
> * **The published Hill slope and the assumed plateau are incompatible in this scheme.** The
>   steepest peak curve this topology can produce with an absolute P_o,max near 0.750 has
>   n_H = 1.33; every parameter set reaching n_H ≥ 1.8 has P_o,max ≥ 0.99. Jahn et al. report
>   n_H = 2.2 ± 0.4. One of the two is wrong, or the scheme is — and note that the slope is
>   published while the 0.750 is this project's own fit target. See
>   `knowledge/12-inference.md` §2.

> **Reproducibility.** Every number is emitted by `scripts/paper_numbers.py` at the cited tag,
> and the Monte Carlo by `scripts/ranking_robustness.py --draws 20000`. None was transcribed.
> This is a testable claim, not an assurance: `scripts/verify_manuscript.py` checks out the
> cited tag into a clean worktree, runs the generators from that tree, and asserts every
> audited figure appears in its output. Run it. Development history, including the defects
> these mechanisms exist to prevent, is in the Supplementary Note, not the main text.

---

## Abstract

Positive allosteric modulators (PAMs) of GABA-A receptors are represented in some multiscale
neural simulations by a bounded scalar potentiation factor — an approximation that obscures the
difference between phasic synaptic and tonic extrasynaptic receptor activation. In earlier
implementations of the present model, PAM action was clamped by an imposed fixed scalar cap
(`gaba_a_efficacy_cap = 2.5`). Here we formulate a five-state continuous-time Markov model with
sequential agonist binding, channel opening and desensitisation, and ask how PAM-induced
slowing of agonist unbinding alters receptor responses under transient and steady-state GABA
exposure.

Under a parameterisation constrained by three macroscopic anchors — activation
EC₅₀ = 20.0 µM, simulated peak open probability at saturating agonist P_o,max = 0.750, and
deactivation τ_IPSC = 15.0 ms — affinity modulation produces markedly different predicted
effects in the two compartments. A brief millimolar synaptic transient drives the receptor far
up its activation curve, so a 2.5-fold leftward EC₅₀ shift raises peak open probability only
**1.02×** (asymptotic ceiling **1.03×**), acting instead on the deactivation tail
(τ **1.66×**, time-integrated open probability **1.66×**). Submicromolar ambient GABA leaves
the steady state far below saturation: at 0.40 µM the same modulator produces a **10.29×**
increase in standing open probability, against an **asymptotic receptor open-probability
dynamic range of 123.3×** in the limit of vanishing agonist dissociation (k_off → 0⁺).

**That 123.3× figure is a ratio between two open probabilities within one kinetic scheme — a
model-specific limiting value over the baseline at a specified ambient GABA concentration. It
is not a prediction that any compound produces a 123.3-fold change in whole-cell current, in
circuit-level inhibition, or in any physiological or clinical endpoint**, each of which depends
additionally on receptor density, single-channel conductance, driving force, chloride
homeostasis and input resistance (§2.3). No finite modulator approaches it: at the modelled
benchmark potency the realised gain is 10.29×, and the asymptote and the reachable gain diverge
as ambient GABA falls (§4.1).

The dynamic range is an emergent property of this topology and parameterisation, not a
pharmacological constant. A one-at-a-time sensitivity sweep shows the exact ratio depends on the
desensitisation equilibrium (d/r) — which the three anchors do not constrain — and on
microscopic affinity, while a large extrasynaptic range (**> 30×**) persists across the
submicromolar ambient concentrations reported for cortex and hippocampus (0.2–0.8 µM). A fixed
2.5× scalar cap binds only above ≈ 4 µM ambient GABA, so in the extrasynaptic regime it
truncates rather than saturates.

Mapped onto subunit distribution estimates, the model separates forebrain from respiratory
actions without any imposed cap, and the resulting compound ordering is robust in two senses: it
survives adverse draws of every uncertain parameter simultaneously, and it survives replacing
the integrate-and-fire neuron model with a conductance-based one sharing almost no mechanism
(Spearman ρ = +1.0000 over five subtype-selective arms). **Both senses are robustness within a
model, not biological validation**, and two findings bound the interpretation sharply. First,
the preBötzinger α5 fraction and the α5 extrasynaptic fraction are both unsourced, so the
anatomical basis of the respiratory-sparing argument is an assumption rather than a measurement.
Second, the non-selective benzodiazepine reference arm — the denominator of every reported
selectivity ratio — reverses sign on respiratory output between the two substrates.

We identify two computational failure modes in multiscale systems models — arbitrary
conductance caps, and transfer of lumped sensitivity parameters between model architectures —
and propose a recombinant-receptor patch-clamp experiment measuring the two-dimensional surface
I_tonic([GABA]_bath, C_PAM) to test the predicted dependence of tonic potentiation on ambient
GABA. That experiment, rather than any number reported here, is what would establish whether
the compartment asymmetry is real.

**Keywords:** GABA-A receptor; α5 subunit; positive allosteric modulator; Markov kinetic model;
tonic inhibition; open-probability dynamic range; model identifiability; reproducibility.

---

## 1. Introduction

GABA-A receptors mediate the principal component of fast inhibitory neurotransmission in the
mammalian CNS (Olsen & Sieghart, 2008). As pentameric ligand-gated chloride channels assembled
from diverse subunit families (α₁₋₆, β₁₋₃, γ₁₋₃, δ, ε, θ, π, ρ₁₋₃), their heterogeneous
anatomical distribution and distinct biophysical properties have made subtype-selective
pharmacology a long-standing objective (Rudolph & Möhler, 2014).

Classical benzodiazepine-site modulators potentiate α1-, α2-, α3- and α5-containing receptors
indiscriminately. Subunit-targeted point mutations (αₙH→R) showed that the resulting
physiological endpoints dissociate across subtypes:

* α1 mediates sedation, anterograde amnesia and medullary respiratory depression;
* α2 and α3 mediate anxiolysis and spinal anti-spasticity;
* α5 modulates hippocampal learning, temporal memory encoding, and ethanol-like discriminative
  stimulus salience (Rudolph et al., 1999; Löw et al., 2000; McKernan et al., 2000;
  Cheng et al., 2006).

Consequently, medicinal chemistry has prioritised subtype-selective PAMs; α5-selective ligands
(alogabat/RG7816, MP-III-022, SH-053-2′F-R-CH₃) have been investigated for Down syndrome,
cognitive deficits, and as candidate leads in the search for better-tolerated anxiolytics
(Nutt, 2006; Atack, 2011).

In the affinity-modulation mechanism examined here, potentiation is modelled as a change in
agonist unbinding. Other GABA-A modulators act through gating transitions, and some directly
gate the channel in the absence of GABA (barbiturates and propofol at high concentration;
muscimol is a classical orthosteric agonist). For a pure affinity-type PAM, potentiation
strictly enhances sensitivity to endogenous GABA without direct gating; once modulator sites
approach saturation, potentiation was commonly assumed — and was approximated in earlier
implementations of this model — to plateau at a fixed scalar cap, typically a 2.0–2.5×
response multiplier. We test whether that scalar representation remains valid across receptor
states and agonist regimes.

**Table 1. Standardised terminology.**

| Concept | Term used here | Symbol / code |
|---|---|---|
| Imposed conductance boundary (2.5×) | fixed scalar cap | `gaba_a_efficacy_cap` |
| Limiting open probability as k_off → 0⁺ | kinetic asymptote | P_open,∞ = **0.131158** |
| Asymptotic bound on open probability at 0.40 µM | asymptotic open-probability dynamic range | P_open,∞ / P_open,base = **123.3×** |
| Maximum operational potency shift | max operational potency-shift factor | s_max = EC₅₀,base / EC₅₀,drug |
| Compound conductance multiplier | finite PAM gain | tonic gain / phasic gain |
| Simulated peak P_o at saturating agonist | P_o,max (calibration anchor) | **0.7500** |
| Analytic gating bound, desensitisation absent | β/(α+β) | **0.7969** |

The last two rows are distinct quantities, easily and consequentially collapsed: β/(α+β) is the equilibrium open probability of the two-site scheme *without* the
desensitised state, whereas the calibration anchor is the peak reached during a finite
saturating step, in which desensitisation competes with activation throughout the rise. The
analytic bound is therefore never attained, and quoting it as the anchor overstates the
receptor's reachable amplitude by 10%.

---

## 2. Theory and parameterisation

### 2.1. Five-state scheme and calibration

We model the receptor as a continuous-time Markov process on the Jones–Westbrook topology
(Jones & Westbrook, 1995; Haas & Macdonald, 1999), with two sequential binding steps, one
conducting state, and a desensitised state entered from the doubly-liganded shut
conformation:

```
                2 kon[G]            kon[G]              beta
    R    <---------------->  RG  <---------->  RG2  <---------->  RG2*
                 koff              2 koff              alpha
                                     |
                                   d | r
                                     v
                                     D
```

In code the states are `("R", "AR", "A2R", "A2O", "A2D")` (`gabaa_kinetics.STATES`); the open
state is `A2O`. The state probability row vector evolves by the forward master equation
dP/dt = P·Q([G](t)), with Q rows summing to zero and the factors of 2 representing two
equivalent binding sites. Rates are in ms⁻¹, k_on in µM⁻¹ms⁻¹.

**Epistemic status.** Three macroscopic targets do not uniquely determine six microscopic
rates. This is **one calibrated parameterisation consistent with the selected constraints**,
not a claim that the microscopic constants are identified.

1. **Anchor 1** — peak-current activation EC₅₀ = 20.0 µM during a 300 ms sustained step.
2. **Anchor 2** — simulated peak open probability at saturating agonist, P_o,max = 0.750,
   evaluated as `po_peak(10 mM)` over a 300 ms step. The module declares this acceptable only
   on [0.70, 0.80] (`FIT_RANGES`), and a regression test now asserts the target lies inside its
   own range.
3. **Anchor 3** — phasic deactivation τ_IPSC = 15.0 ms following a synaptic transient.

The synaptic exposure is an exponentially cleared pulse
[G](t) = [G]_peak · exp(−t/τ_clear) with **[G]_peak = 1.0 mM** and **τ_clear = 0.30 ms**
(`SYNAPTIC_PEAK_UM`, `SYNAPTIC_CLEAR_MS`), integrated from P(0) = [1,0,0,0,0] over a 200 ms
window. τ_IPSC is defined by a single-exponential fit in log space over the standard
experimental 90%→10% post-peak window, rather than a raw 90→10% duration. Charge transfer is
the time-integrated open probability over the window,
Q_proxy = ∫₀²⁰⁰ P_open(t) dt, in ms; at fixed single-channel conductance and driving force this
is proportional to synaptic charge, with 13.053 ms ≡ 0.013053 pC/pA.

**Optimiser.** Bounded non-linear least squares (`scipy.optimize.least_squares`, Trust Region
Reflective) over x = [log₁₀k_on, log₁₀k_off, log₁₀β, log₁₀α], residuals
r₁ = ln(EC₅₀/20.0), r₂ = (P_o,max − 0.750)/0.05, r₃ = ln(τ_IPSC/15.0). EC₅₀ is extracted
numerically from peak open probability across a 30-point logarithmic grid
[G] ∈ [0.1, 10⁴] µM, integrated by LSODA (rtol 1e-9, atol 1e-12), with log-linear
interpolation at half-maximal peak.

Desensitisation rates (**d = 0.050 ms⁻¹, r = 0.0020 ms⁻¹**, d/r = 25.0) govern steady-state
desensitisation, are unconstrained by the three transient anchors, and are held at fixed
provisional values. **The asymptotic dynamic range is therefore strictly conditional on this
provisional desensitisation equilibrium**, which §2.5 quantifies.

**Fitted rates** (`scripts/paper_numbers.py --section kinetics`):

| Rate | Value |
|---|---|
| k_on | 0.0075785 µM⁻¹ms⁻¹ |
| k_off | 0.180415 ms⁻¹ |
| β | 1.177468 ms⁻¹ |
| α | 0.300000 ms⁻¹ |
| d | 0.050 ms⁻¹ |
| r | 0.0020 ms⁻¹ |

Reproduced anchors: EC₅₀ = **20.0000 µM**, P_o,max = **0.7500**, τ_IPSC = **15.000 ms** — all
three to the precision quoted.

### 2.2. Microscopic equilibrium versus functional potency

Keeping binding parameters separate from concentration-response endpoints:

* **Fitted targets:** EC₅₀ = 20.0 µM, P_o,max = 0.750, τ_IPSC = 15.0 ms.
* **Derived microscopic quantities:** K_d, gating efficacy E.
* **Model predictions:** P_open(0.40 µM) = 1.0634 × 10⁻³, asymptotic range 123.3×.

K_d = k_off/k_on = 0.180415 / 0.0075785 = **23.81 µM**.
E = β/α = 1.177468 / 0.300000 = **3.9249**; α/β = **0.25478**.

For the two-site scheme *without* desensitisation, equilibrium open probability is

P_open([G]) = E[G]² / (K_d² + 2K_d[G] + (1+E)[G]²)

so as [G] → ∞, P_open → E/(1+E) = β/(α+β) = **0.7969**. Half-maximal open probability satisfies
(1+E)[G]²₁/₂ − 2K_d[G]₁/₂ − K_d² = 0, whose positive root is

**[G]₁/₂ = K_d (1 + √(2+E)) / (1 + E)**

Substituting: [G]₁/₂ = 23.81 × 3.4341 / 4.9249 = **16.60 µM**. High gating efficacy
(E = 3.92 > 1) pulls the equilibrium forward, shifting the half-maximal concentration leftward
from K_d = 23.81 µM to 16.60 µM. When desensitisation (d/r = 25.0) is included and evaluated
numerically over a 300 ms application, the peak-current EC₅₀ is **20.000 µM**, reproducing the
calibration target — reconciling microscopic binding with macroscopic activation, and showing
that the analytic equilibrium midpoint and the simulated peak EC₅₀ differ by 20%
((20.000 − 16.60)/16.60 = 20.5%; the convention is stated because the previous draft's "8%"
named neither the numerator nor the base and matched neither).

### 2.3. Asymptotic steady-state open probability as k_off → 0⁺

An **affinity-type PAM** increases apparent agonist affinity by slowing unbinding by
c_affinity ≥ 1: k_off^drug = k_off^base / c_affinity.

We evaluate [G] = 0.40 µM as a representative ambient scenario within the range reported for
selected preparations (0.2–0.8 µM; Farrant & Nusser, 2005). Solving P·Q = 0:

**P_open,baseline(0.40 µM) = 1.0634 × 10⁻³** (0.1063%).

In the limit c_affinity → ∞ (k_off → 0⁺) at any [G] > 0, unbinding vanishes; R and RG become
transient and their stationary probabilities go to zero. All probability concentrates in the
closed sub-chain {RG₂, RG₂*, D}. Detailed balance on each pair gives
P_RG₂ = (α/β)·P_open and P_D = (d/r)·P_RG₂, so normalisation yields

**P_open,∞ = 1 / [ 1 + (α/β)(1 + d/r) ]**

Substituting α/β = 0.25478 and d/r = 25.0:

P_open,∞ = 1 / (1 + 0.25478 × 26.0) = 1 / 7.6244 = **0.131158** (13.12%).

Computed two independent ways — this closed form, and the mechanism pushed numerically to
c_affinity = 10⁵ — the two agree to **8.5 × 10⁻⁶**, and `paper_numbers.py` prints both so a
divergence would be visible. Dividing by baseline:

**P_open,∞ / P_open,baseline(0.40 µM) = 0.131158 / 1.0634×10⁻³ = 123.3×**

#### Interpretation and necessary distinctions

Under this topology, parameterisation and ambient scenario, the model permits a **123.3-fold**
asymptotic increase in steady-state open probability. With the desensitisation branch included,
steady-state open probability at *saturating* GABA also settles to this same three-state
asymptote (≈ 0.1564), while the early transient peak reaches 0.750 in simulation
(analytic bound 0.7969).

Hierarchical observables must not be conflated:

* **P_open** — a microscopic channel-state variable.
* **Whole-cell current** — I_tonic = N · γ · P_open · (V_m − E_Cl), so receptor density,
  single-channel conductance and driving force all intervene.
* **Tonic inhibition** — depends on chloride homeostasis (KCC2/NKCC1); prolonged
  high-conductance chloride flux can shift E_Cl depolarising.
* **Circuit shunting** — depends non-linearly on input resistance relative to leak and
  synaptic conductances.

A 123.3× open-probability dynamic range is therefore **not** a 123.3× increase in circuit-level
inhibition. It is the upper bound on receptor open probability within this scheme.

### 2.4. Mapping operational potency shifts onto Q

To connect phenomenological compound descriptors to microscopic rates:

1. **s_max** — the maximum fold leftward shift of the agonist concentration-response curve at
   full modulator occupancy: s_max = EC₅₀(baseline) / EC₅₀(c_affinity).
2. **Implicit inversion** — the unique c_affinity ≥ 1 reproducing a given s_max is found by
   root-solving F(c) = EC₅₀(Q(k_off)) / EC₅₀(Q(k_off/c)) − s_max = 0 via Brent's method, with
   EC₅₀ evaluated by integrating the master equation over a 300 ms step.

**These two numbers are not interchangeable.** In this parameterisation:

| s_max | c_affinity |
|---|---|
| 2.40 | 3.2566 |
| **2.50** | **3.4674** |

and conversely c = 3.47 gives s_max = **2.501**. At c = 3.4674 the model's EC₅₀ moves
20.0000 → **8.0000 µM**, i.e. exactly 2.5×.

3. **Occupancy weighting.** At fractional modulator occupancy θ = [D]/([D] + K_d,PAM), the
   receptor pool is a linear mixture of unmodulated and fully modulated channels, so standing
   tonic gain scales linearly:
   Tonic Gain(θ) = 1 + θ·(Gain_max(s_max) − 1).

**Table 2. Reported endpoints and their mapping to s_max.** The `ceiling` column is the value
the model uses (`PROFILES[...].ceiling`, consumed as `s_max` at `evaluation.py:86`). The final
column is the result of a **claim-support** read (§7), not a metadata check: it records whether
the cited source supports *the mapped number*, as distinct from existing and being on-topic.

| Compound | Reported primary endpoint | Assay / construct | `ceiling` | Does the source support this number? |
|---|---|---|---|---|
| Imepitoin | ~20% max potentiation rel. diazepam (Rundfeldt & Löscher, 2014) | recombinant α1β2γ2, oocyte | 1.25 | direction only |
| TPA023 | α2/α3 partial; silent antagonist at α1 and α5 | recombinant | 2.00 | not re-audited |
| L-838,417 | partial efficacy ~30–40% at α2/3/5; silent at α1 (McKernan et al., 2000) | membrane potential | 2.20 | direction only |
| MP-III-022 | binding- and efficacy-selective α5 PAM; potentiation "mild to moderate to strong"; non-α5 receptors engaged only at the top dose (Stamenić et al., 2016) | patch clamp, α5β3γ2 | 2.50 | **no — direction only; no EC₅₀ shift reported** |
| Alogabat (RG7816) | potent α5β3γ2 PAM with binding and functional selectivity (Cecere et al., 2025) | recombinant α5β3γ2 patch | 2.50 | **no — selectivity supported; the +167%/+72% EC₂₀ figures are not in the abstract and were not verified against full text** |
| Diazepam / midazolam | *biphasic* potentiation with separable nanomolar and micromolar components (Walters et al., 2000) | recombinant α1β2γ2 patch | 2.50 | **no — see note (iii)** |
| Neurosteroid (allopregnanolone proxy) | gating-active, non-selective | — | 6.00 | not re-audited |
| Gaboxadol (THIP) | δ orthosteric agonist | — | 1 × 10⁹ | not re-audited |

**Notes.**

(i) **No `ceiling` in this table is a measured EC₅₀ fold shift for its compound.** Reported
experimental metrics — relative efficacy percentages, partial agonism, current potentiation at
EC₂₀ — are different quantities. Every mapping is an explicit modelling assumption, and the
right column records that the claim-support read did not upgrade any of them.

(ii) **The diazepam benchmark is UNSOURCED, and the citation that was supposed to carry it says
something more interesting.** Walters et al. (2000) is correctly cited here (*Nat. Neurosci.*
3:1274–1281; an earlier draft gave *Br. J. Pharmacol.* 131:1307–1314 under a different title),
but it reports no 2.5-fold EC₅₀ shift. What it reports is that diazepam potentiates α1β2γ2 in
**two separable phases**, nanomolar and micromolar, dissociable by TM2 mutations. That is not
support for `ceiling = 2.50`; it is evidence that **a single affinity parameter is an
oversimplification of this compound** — the one compound whose value propagates into every
selectivity ratio in §4.4 as the denominator. We keep 2.50 as a stated benchmark, label it
UNSOURCED, and flag the mechanistic tension rather than letting the citation imply support it
does not give.

### 2.5. Sensitivity and robustness of the dynamic range

Because the 123.3× ASYMPTOTIC ratio (the k_off → 0⁺ limit of §2.3, not a reachable
gain) depends on the parameter set, we swept each quantity one at a time
(`scripts/paper_numbers.py --section sensitivity`).

**Why this section exists in place of confidence intervals.** The fit is a **square system**:
three macroscopic anchors, three free rates, with `α` held at a declared convention. It
therefore has a unique solution and **no residual degrees of freedom, so no interval comes
out of it.** That is not a small uncertainty — it is the absence of a measurement of
uncertainty, and the two must not be confused. The real uncertainty lives in three places
the fit cannot see:

1. **The anchors themselves**, which are recalled literature ranges rather than digitised
   points (see the scope box). The one sourced peak EC₅₀ for this preparation is 11.6 ± 0.9 µM
   against the 20 µM used here.
2. **`d/r`, which no anchor constrains at all.** The sweep below is the honest substitute for
   an interval on it: ±2-fold moves the asymptotic range over 67.7–211.1×.
3. **`α`, a convention.** Without it the system is four unknowns against three residuals, and
   its solution moved between SciPy versions.

Where a posterior *does* exist — over the three combinations an equilibrium
concentration-response can determine — the asymptotic tier spans **2.99-fold** at 0.09
decades of spread in `log10 D` (`dynamic_range.dynamic_range_posterior`;
`knowledge/12-inference.md` §3). Every point estimate in this manuscript should be read
against that, not against the four decimal places it is printed to.

**Table 3. One-at-a-time sensitivity of the modelled open-probability dynamic range.**

| Perturbation | P_open(0.40 µM) | P_open,∞ | dynamic range |
|---|---|---|---|
| **nominal** | 0.0010634 | 0.13116 | **123.3×** |
| desensitisation entry d × 0.5 | 0.0010670 | 0.22525 | 211.1× |
| desensitisation entry d × 2.0 | 0.0010562 | 0.07146 | 67.7× |
| resensitisation r × 0.5 | 0.0010562 | 0.07146 | 67.7× |
| resensitisation r × 2.0 | 0.0010670 | 0.22525 | 211.1× |
| gating β × 0.5 | 0.0005320 | 0.07018 | 131.9× |
| gating β × 2.0 | 0.0021244 | 0.23190 | 109.2× |
| affinity K_d × 0.5 | 0.0040257 | 0.13116 | 32.6× |
| affinity K_d × 2.0 | 0.0002719 | 0.13116 | 482.4× |
| ambient 0.10 µM | 0.0000686 | 0.13116 | 1910.8× |
| ambient 0.20 µM | 0.0002719 | 0.13116 | 482.4× |
| **ambient 0.40 µM (nominal)** | 0.0010634 | 0.13116 | **123.3×** |
| ambient 0.70 µM | 0.0031285 | 0.13116 | 41.9× |
| ambient 0.80 µM | 0.0040257 | 0.13116 | 32.6× |
| ambient 1.50 µM | 0.0125176 | 0.13116 | 10.5× |
| ambient 3.00 µM | 0.0360856 | 0.13116 | 3.6× |

**Finite modulator shifts** (the realisable, as opposed to asymptotic, gains):

| c_affinity | s_max | tonic gain at 0.40 µM |
|---|---|---|
| 1.50 | 1.385 | 2.19× |
| 2.00 | 1.723 | 3.79× |
| 2.79 | 2.167 | 7.00× |
| 3.4674 | 2.500 | 10.29× |
| 5.00 | 3.095 | 18.89× |
| 10.00 | 4.310 | 47.78× |
| 50.00 | 6.513 | 109.45× |

#### Key observations

1. **Desensitisation ratio (d/r).** The asymptote is governed directly by d/r. Doubling
   desensitisation entry (d/r = 50) leaves **67.7×**; halving it (d/r = 12.5) expands the
   range to **211.1×**. This is the single largest source of uncertainty, and it is the one
   parameter the three transient anchors do not constrain.
2. **Gating ratio (β/α).** Perturbing β by ±2-fold moves baseline and asymptote nearly in
   proportion, leaving the range stable at **109.2–131.9×**.
3. **Ambient GABA.** Across the 0.2–0.8 µM range reported for cortex and hippocampus the range
   spans **482.4× to 32.6×**, so the robust claim across that range is **"> 30×"** (41.9× at
   0.70 µM). THE PREVIOUS DRAFT SAID "> 50×", which was true of the superseded fit (lower end
   54.7×) and is false of this one by a factor of 1.5. It is recorded rather than quietly
   dropped because it is the clearest example in this manuscript of a rounded claim that
   survived the parameters it was computed from — the number moved 54.7 → 32.6 and the
   sentence did not.
4. **K_d and ambient are the same axis.** Halving K_d at fixed [G] gives exactly the numbers of
   doubling [G] (0.0040257, 32.6×), as it must, since the equilibrium depends on [G]/K_d. Both
   rows are retained because the two have different experimental meanings.
5. **When a 2.5× cap would actually bind.** With x = [G]/K_d and A = 1 + E + d/r = **29.9249**,
   the stationary ratio is

   P_open,∞ / P_open([G]) = 1 + 2/(Ax) + 1/(Ax²)

   Setting this to 2.5 and solving the quadratic 1.5Ax² − 2x − 1 = 0 gives x = 0.17319, i.e.

   **[G] = 4.12 µM** (verified: the range at 4.12 µM is 2.500×).

   A fixed 2.5× scalar cap therefore binds only when ambient GABA exceeds ≈ 4 µM — an order of
   magnitude above the reported submicromolar range. In the regime the extrasynaptic pool
   actually occupies, the cap is not a saturation boundary; it is a truncation.

### 2.6. Comparison with gating modulation

Under constant non-zero agonist, steady-state open probability in the affinity-only limit
approaches P_open,∞ = **0.1564**. This stationary limit does not bound transient peaks, which
depend on the agonist waveform and reach **0.4196** (baseline) to **0.7508** (limiting) in the
synaptic cleft, strictly below the analytic gating bound β/(α+β) = **0.7969**.

A gating PAM (β → β·c) instead scales the forward opening transition:

lim_{c→∞} P_open = 1 / [1 + (α/βc)(1 + d/r)] → **1.0**

So increasing the forward gating rate can drive both steady-state and transient P_open toward
unity, whereas steady-state affinity-only modulation remains bounded by the desensitised-state
equilibrium at 0.1564. This distinction tracks the divergent pharmacology of
neurosteroids and barbiturates (gating-active) versus benzodiazepines (predominantly
affinity-modulating), and it is why the neurosteroid arm in §4.4 carries `ceiling = 6.0`.

---

## 3. Model implementation

The receptor scheme of §2 is implemented in `circuitpharm.gabaa_kinetics` and drives two
downstream layers: an algebraic selectivity index over subunit distributions
(`circuitpharm.subtypes`, §4.4) and a spiking respiratory circuit available on two neuron
substrates (`circuitpharm.resp` with `circuitpharm.neuron`, §4.6). Reported quantities carry
reliability tiers (`circuitpharm.results.Tier`): `Tier.VOID` marks quantities lacking
whole-animal anchors — absolute lethal-dose margins in mg — and blocks them from analysis
programmatically, while `Tier.VALIDATED` marks quantities invariant to the shared multiplicative
tissue-scaling constants. The selectivity index R is the latter: because those constants multiply
numerator and denominator alike they cancel, so the compound ordering is invariant to that scale.
**Cancellation establishes scale-invariance within the model's algebra. It does not validate the
biological accuracy of the subunit fractions or circuit weights** — see §4.4.4.

A third guardrail operates on reported statistics rather than parameters:
`scripts/ranking_robustness.py` prints an explicit **NOT QUOTABLE** verdict for its Monte Carlo
median and 95th percentile (§4.4.3), and this manuscript honours it.

Two modelling failure modes shaped the present implementation and are relevant to anyone
building comparable models: an **arbitrary conductance cap** (`gaba_a_efficacy_cap = 2.5`), which
manufactured a safety guarantee by clamping the quantity under study, and **transfer of a lumped
sensitivity parameter between architectures**, which collapsed a non-selective benzodiazepine's
modelled respiratory effect to −2% against a −16% to −19% clinical anchor and so reported that
every compound survived every dose. Both are described, with their diagnostic histories, in
**Supplementary Note §S2**; §5.1 states the design recommendation that follows. The second taught
the substantive lesson that `gaba_sens` is **not a subunit fraction** but a lumped parameter
absorbing extra-preBötC mechanisms — chemoreflex blunting, upper-airway motor tone — that an
isolated pacemaker model does not contain.

## 4. Results

### 4.1. Non-equilibrium gating dissociates synaptic and extrasynaptic responses

Simulating an affinity-type PAM calibrated to a 2.5× leftward EC₅₀ shift
(**c_affinity = 3.4674**) reveals a pronounced divergence between compartments.

**Table 4. Divergence of kinetic predictions across receptor compartments.**

| Compartment / observable | Agonist regime | Baseline | With PAM | Ratio |
|---|---|---|---|---|
| Simulated synaptic peak P_o | transient, 1.0 mM, τ_clear 0.30 ms | 0.727592 | 0.743186 | **1.021×** (peak gain) |
| Synaptic τ_IPSC | mono-exp fit, 90→10% decay | 15.000 ms | 24.847 ms | **2.495×** |
| Time-integrated P_o (charge proxy) | 200 ms window | 13.053 ms | 21.605 ms | **2.202×** |
| Synaptic max dynamic range | limiting transient, k_off → 0⁺ | 0.727592 | 0.750800 | **1.032×** |
| Tonic standing P_o | steady state, 0.40 µM ambient | 0.000788 | 0.005688 | **10.285×** (tonic gain) |
| Tonic max dynamic range | limiting steady state, k_off → 0⁺ | 0.000788 | 0.131158 | **123.339×** |

*Charge proxy: 13.053 ms ≡ 0.013053 pC/pA; 21.605 ms ≡ 0.021605 pC/pA. The simulated synaptic
peaks (0.4196 baseline, 0.5536 with PAM, 0.7508 limiting) all respect the analytic gating bound
β/(α+β) = 0.7969, which is not itself attainable because desensitisation competes with
activation throughout the rise.*

1. **Synaptic cleft.** The brief 1 mM transient drives peak open probability to 0.420. A
   leftward affinity shift raises that peak only **1.021×**, with an asymptotic ceiling of
   **1.032×** — so roughly two-thirds of the available synaptic headroom is already consumed by
   a 2.5× shift. Instead, slower unbinding prolongs the deactivation tail: τ_IPSC rises
   **2.495×** (15.00 → 24.85 ms) and integrated open probability **2.202×**, reproducing the
   classical electrophysiological phenotype of benzodiazepine action on IPSC decay
   (Otis & Mody, 1992).
2. **Extrasynaptic space.** Receptors at steady-state 0.40 µM operate far below saturation
   (P_open = 7.9 × 10⁻⁴). Slower unbinding moves the activation threshold into the ambient
   range, producing a **10.285×** potentiation of standing open probability, with theoretical
   headroom extending to **123.3×**.

**The contrast between the two compartments' ceilings is 123.339 / 1.032 = 120-fold.** This is
the manuscript's central quantitative claim. It is sensitive to the synaptic ceiling, which at
this parameterisation is 1.03×; the qualitative conclusion — that two orders of magnitude
separate the compartments, so a single scalar multiplier cannot describe both — is not.

A caution that must travel with every figure above. **Headroom and reachable gain are different
quantities, and they diverge in opposite directions as ambient GABA falls.** Between 0.40 and
0.10 µM the asymptotic headroom rises 15.6-fold (123.3× → 2881×) while the gain a finite
s_max ≈ 2.5 modulator actually extracts rises only 1.06-fold (10.29× → 8.47×). Any statement
about what a real compound does must be written against the reachable column. §5.2 returns to
this, because it is where a falsification criterion can be written against the wrong one.

### 4.2. Subunit distribution and anatomical decoupling: the α5 case

**Table 5. Subunit compartmentalisation and regional distribution estimates, with per-cell
provenance.** Basis labels are generated from `circuitpharm.provenance.report()`, not assigned
by hand. **No cell in this table is QUANTITATIVE**; eight are UNSOURCED, one is an explicit
GUESS, and the six labelled FROM_QUALITATIVE take a *direction* from a source while the number
itself is ours.

| Subunit | Extrasyn. `f_extra` | basis | preBötC `f` | basis | Forebrain `f` | basis |
|---|---|---|---|---|---|---|
| α1 | 0.15 | **UNSOURCED** | 0.60 | FROM_QUALITATIVE [`pbc_alpha`] | 0.35 | **UNSOURCED** |
| α2/3 | 0.20 | **UNSOURCED** | 0.15 | FROM_QUALITATIVE [`pbc_alpha`] | 0.25 | **UNSOURCED** |
| α5 | 0.80 | **UNSOURCED** | 0.02 | **UNSOURCED** | 0.30 | FROM_QUALITATIVE |
| δ/α4 | 1.00 | FROM_QUALITATIVE [`pbc_delta`] | 0.15 | FROM_QUALITATIVE [`pbc_delta`] | 0.08 | **UNSOURCED** |
| ε | 0.50 | **GUESS** | 0.08 | FROM_QUALITATIVE [`pbc_eps`] | 0.02 | **UNSOURCED** |

#### Provenance, stated at the resolution a claim-support read supports

The repository's audit reports that **6 of 28 load-bearing parameters name a source (21%)**.
For this table specifically: **eight of fifteen cells are UNSOURCED, a ninth is a GUESS, and
the remaining six are FROM_QUALITATIVE.**

Two of those cells are load-bearing, and **both are UNSOURCED**:

* **`f_α5` in preBötC = 0.02** — the number on which the respiratory-sparing argument rests. It
  encodes an anatomical ordering (medullary α5 low, forebrain α5 high) and is not a measurement
  from isolated preBötC tissue. §4.4.3 quantifies what it does to the results.
* **`f_extra,α5` = 0.80** — the number on which the extrasynaptic-headroom argument rests.
  **This was downgraded from FROM_QUALITATIVE to UNSOURCED while preparing this version**, and
  the reason is worth stating in full because it is a reproducible mistake.

#### Why the α5 extrasynaptic fraction is now UNSOURCED

An earlier draft attributed `f_extra,α5 = 0.80` to Kasugai et al. (2010). That is the obvious
candidate and exactly the right kind of study: quantitative freeze-fracture replica immunogold
labelling, synaptic against extrasynaptic pools, hippocampal CA1 pyramidal cells, in a good
journal. Reading its abstract rather than its title shows that **it measured α1, α2 and β3 —
not α5** — and that its quantitative result points the other way: synaptic labelling density
exceeded extrasynaptic density by **78–132× (α1), 94× (α2) and 79× (β3)**. It therefore supports
a *low* extrasynaptic fraction for the subunits it did measure, and is silent on this one.

This is the **second** source in this project to resolve perfectly by DOI and fail to support
the number attached to it. The first was `a5_dist` — a 1988 study using a single generic cDNA
probe for "the α subunit", reporting total α-subunit mRNA by region (medulla ≪
hippocampus/cortex), which is a claim about regional *level* where `REGIONS` encodes regional
*composition*. The pattern is specific enough to name: **a source whose title matches the claim,
in the right journal, by the right group, measuring a neighbouring quantity.** Metadata
verification cannot catch it; only reading can. §7 reports a claim-support read of every
load-bearing citation, with results.

The direction — that α5 is enriched extrasynaptically relative to α1 — is widely stated in
review literature and we have no reason to doubt it. What we do not have is a primary source in
this repository that supports **0.80**, and the honest consequence is that the extrasynaptic
pool's size is an assumption of the same standing as the preBötC α5 fraction.

Two other provenance caveats travel with this table:

* The preBötC α1/α2/α3 column rests on `pbc_alpha` — Liu & Wong-Riley (2004),
  *J. Appl. Physiol.* — whose metadata is verified by DOI but whose **full text is inaccessible
  (HTTP 403)**, so its claim support is formally **UNASSESSED**. It covers α1, α2 and α3 only,
  not α5, δ, α4 or ε, to which `REGIONS["prebotc"]` also assigns values; and it is a
  *developmental* study.
* Pirker et al. (2000) maps distribution qualitatively and is cited only for the qualitative
  ordering (high medullary α1, low medullary α5), not for any number in this table.

#### What follows, and what does not

1. **Mechanistic basis of respiratory sparing — conditional on an unsourced number.** Low α5
   representation in the preBötC *would* provide a mechanistic basis for reduced direct
   α5-mediated inhibition of the core rhythm generator. Since `f_α5,preBötC` is unsourced, this
   is a conditional statement about the model, not evidence about biology. Whole-animal
   respiratory preservation could not be inferred from preBötC receptor abundance even if the
   number were measured, because ventilatory stability also depends on chemoreflex integration
   and upper-airway motor tone — the mechanisms `gaba_sens` lumps (Supplementary Note §S2).
2. **Forebrain tonic potentiation — conditional on a second unsourced number.** If α5
   expression in hippocampus and cortex is substantial (≈ 30%, FROM_QUALITATIVE) and
   predominantly extrasynaptic (0.80, UNSOURCED), an α5 PAM acts chiefly on the high-headroom
   pool. Both legs of that inference are assumptions.

So the model's behaviour is clear and its anatomical grounding is not. Within the model, a
123.3× receptor dynamic range does not produce 123.3× respiratory depression, and dose
escalation of a high-efficacy α5 PAM drives substantial tonic shunting in forebrain circuits
with cognitive blunting and sedation as the expected costs. **Whether that corresponds to
anything in tissue depends on two numbers that no source in this repository supports.**
Replacing them with sourced ranges and re-running the robustness analysis is the highest-value
outstanding work, and §5.3 says so.

### 4.3. Dependence on the maximum operational potency shift

Because extrasynaptic headroom is large, potentiation is limited by the ligand's own s_max —
the maximum left-shift it can elicit at full occupancy:

Tonic Gain(occupancy) = 1 + occupancy · (Gain_max(s_max) − 1)

* **Low potency-shift modulators.** Imepitoin (`ceiling` = 1.25, c = 1.2697) meets a narrow
  biophysical ceiling under this mechanism: at 100% occupancy and 0.40 µM ambient GABA its
  tonic open-probability gain cannot exceed **1.60×**, against a phasic peak gain of 1.08×.
  Tonic gain rises steeply from there — TPA023 at `ceiling` = 2.00 reaches **4.42×** and
  L-838,417 at `ceiling` = 2.20 reaches **5.46×**. The model therefore predicts substantially
  greater tonic headroom once s_max exceeds ≈ 1.3–1.5; that range should be **experimentally
  evaluated as a candidate operational boundary, not treated as an established safety
  threshold**.
* **High-efficacy scenarios.** For alogabat, recombinant electrophysiology shows selective
  potentiation of EC₂₀ GABA currents (+167% rat, +72% human α5β3γ2; Cecere et al., 2025).
  Rather than asserting a universal measured EC₅₀ fold shift, s_max = 2.50 (c = 3.4674,
  **10.29×** tonic gain) is examined as a modelled high-efficacy benchmark calibrated to that
  profile. MP-III-022 carries the same `ceiling = 2.50` in the code, based on partial α5
  potentiation data (Stamenić et al., 2016).
* **High potency-shift modulators (s_max ≥ 2.5).** Diazepam at `ceiling` = 2.50 (c = 3.4674)
  gives **10.29×** tonic gain. The relationship is strongly supralinear beyond that: c = 5.00
  (s_max 3.71) gives **20.64×**, c = 10.0 (s_max 5.73) gives **59.22×**, and the gating-active
  neurosteroid profile at `ceiling` = 6.00 (c = 8.8799) gives **50.52×** tonic gain against a
  phasic peak gain of only 1.10×. The neurosteroid arm is the clearest case of the compartment
  asymmetry: a 46-fold difference between what the compound does to a standing extrasynaptic
  conductance and what it does to a synaptic peak.

### 4.4. Algebraic selectivity index and Monte Carlo robustness

> **Scope.** This section is **algebraic, not multiscale.** `scripts/ranking_robustness.py`
> imports `circuitpharm.subtypes` and nothing else — no neuron, no circuit, no simulation. R is
> a ratio of weighted subunit sums. Describing it as a multiscale or circuit pipeline would
> overstate it by a whole layer of the model. Circuit-simulation results appear separately, in
> §4.6.

#### 4.4.1. Pipeline and reproducibility specification

**Table 6A. Complete inputs** (verbatim from `circuitpharm.subtypes`; subtype order
[α1, α2/3, α5, δ/α4, ε]).

| Parameter | Symbol | Values |
|---|---|---|
| preBötC regional fractions | `f_prebotc` | [0.60, 0.15, 0.02, 0.15, 0.08] (Σ = 1.00) |
| Forebrain regional fractions | `f_forebrain` | [0.35, 0.25, 0.30, 0.08, 0.02] (Σ = 1.00) |
| Extrasynaptic fractions | `f_extra` | [0.15, 0.20, 0.80, 1.00, 0.50] |
| Subjective drive weights | `w_subj` | [0.00, 1.00, 1.00, 0.00, 0.00] |
| Respiratory burden weights | `w_resp` | [1.00, 1.00, 1.00, 1.00, 1.00] |
| Kinetic gain weighting (nominal) | ρ | **6.00** (FITTED) |

| Compound | Efficacy vector `e` | `ceiling` |
|---|---|---|
| Ideal α5 PAM | [0.00, 0.00, 1.00, 0.00, 0.00] | 2.50 |
| Alogabat (RG7816) | [0.00, 0.10, 1.00, 0.00, 0.00] | 2.50 |
| MP-III-022 | [0.00, 0.15, 1.00, 0.00, 0.00] | 2.50 |
| HZ-166 / KRM-II-81 | [0.00, 1.00, 0.00, 0.00, 0.00] | 2.50 |
| Non-selective BZ (diazepam) | [1.00, 1.00, 1.00, 0.00, 0.00] | 2.50 |
| Neurosteroid (allopregnanolone proxy) | [1.00, 1.00, 1.00, 1.00, 0.80] | 6.00 |
| Gaboxadol (THIP) | [0.00, 0.00, 0.00, 1.00, 0.00] | 1 × 10⁹ |
| SH-053-2′F-R-CH₃ | [0.00, 0.05, 0.80, 0.00, 0.00] | 2.50 |
| SH-053-2′F-S-CH₃ | [0.00, 0.50, 0.50, 0.00, 0.00] | 2.50 |

**`w_subj` provenance.** Only the α5 weight has any source: `a5_disc` (α1GABA-A and α5GABA-A
contributions to the discriminative stimulus effects of ethanol in squirrel monkeys, 2005;
PMID 15650112) reports that α5 agonists mimic ethanol's discriminative stimulus and the α5
inverse agonist L-655,708 blocks it. That supports the **direction**; the weight 1.0 is ours.
The α2/3, α1, δ and ε weights are all **UNSOURCED** — and `w_subj,α1 = 0.0` is **actively
questionable**, since the one source bearing on it is titled for an α1 contribution. An
attribution of these weights to "Saba et al. (2017), *Alcohol. Clin. Exp. Res.* 41:748–758" does
not resolve against Crossref and is withdrawn rather than replaced (§7).

The computation proceeds in five deterministic steps, with
w(s) = ρ·f_extra,s + (1 − f_extra,s):

1. **Efficacy vector e** per compound (Table 6A).
2. **Kinetic gain weighting ρ.** The operational weighting of extrasynaptic (tonic) against
   synaptic (phasic) modulation. For s_max = 2.50 the tonic gain is 10.285×, so
   ρ_peak = 10.285 / 1.021 = **6.79** and ρ_charge = 10.285 / 2.202 = **4.36**. The nominal
   **ρ = 6.00** therefore now sits essentially *at* the peak-based ratio rather than between the
   two, so the nominal value effectively weights the peak-based reading. Monte Carlo sweeps
   sample ρ ~ U[3.0, 7.5], which spans both endpoints.
3. **Forebrain subjective drive.** Drive_subj = Σ_s f_forebrain,s · e_s · w_subj,s · w(s)
4. **PreBötC respiratory burden.** Burden_resp = Σ_s f_prebotc,s · e_s · w(s)
5. **Selectivity ratio R**, normalised to diazepam (§3.3).

**Table 6. Selectivity ranking — a model-internal index, not a measured property.** Values are
rounded to two significant figures deliberately: the underlying arithmetic is exact, but the
index is computed from efficacy vectors that are stylised modelling assumptions, nine subunit
fractions with no source (§4.2), four subjective-drive weights with no source (§4.4.1), and one
fitted weighting parameter (ρ). Reporting R = 8.6442 would imply a precision the inputs cannot
carry. Exact values are available from `scripts/paper_numbers.py --section selectivity`.

*The Monte Carlo median and 95th percentile are deliberately omitted from this table: the model
marks them **NOT QUOTABLE** (§4.4.3), and they appear there as a prior-sensitivity diagnostic
rather than as results.*

| Compound | Profile | `ceiling` | Nominal R | 5th-pct R | P(R > 1) |
|---|---|---|---|---|---|
| Ideal α5 PAM | α5-exclusive | 2.50 | **11** | 2.6× | 99.9% |
| SH-053-2′F-R-CH₃ | α5-selective enantiomer | 2.50 | **9.3** | — | — |
| Alogabat (RG7816) | α5-selective (Phase II) | 2.50 | **8.6** | 2.5× | 99.9% |
| MP-III-022 | α5-selective | 2.50 | **7.9** | 2.5× | 99.9% |
| SH-053-2′F-S-CH₃ | α2/α3/α5 enantiomer | 2.50 | **3.6** | — | — |
| HZ-166 / KRM-II-81 | α2/α3-preferring | 2.50 | **1.2** | 0.34× | **63.8%** |
| Non-selective BZ | α1/2/3/5 (diazepam) | 2.50 | **1.0** (ref) | — | — |
| Neurosteroid | non-selective, gating PAM | 6.00 | **0.56** | 0.33× | **0.0%** |
| Gaboxadol | δ orthosteric agonist | 1 × 10⁹ | **0.00** | 0.00× | 0.0% |

**Notes.**
1. **R is model-conditional.** It measures the ratio of forebrain subjective drive to preBötC
   pacemaker burden, relative to diazepam, and the ranking is conditional on every input in
   Table 6A. **R does not compute clinical safety, human respiratory depression, or an overdose
   margin**, and `Tier.VOID` blocks absolute margins programmatically.
2. **s_max does not enter R.** The `ceiling` column is used in dose escalation and gating
   bounds. R evaluates non-saturating proportional scaling governed by **e** and ρ alone.
3. **Gaboxadol's zero is a scope artefact, not a pharmacological verdict.** Because `w_subj` is
   parameterised for ethanol-like discriminative stimulus salience mediated by α2/3 and α5,
   `w_subj,δ = 0.0` forces Drive_subj = 0 and hence R = 0. This implies **nothing** about
   gaboxadol's hypnotic, sedative or thalamocortical efficacy.
4. **HZ-166 does not clearly separate**, at 63.8% of draws above the diazepam reference and a
   5th percentile of 0.34×.
5. **The reference arm is unstable across neuron models** (§4.6), so every value in the
   "Nominal R" column inherits that instability even though the *ordering* does not.

#### 4.4.2. Uncertainty model: exact Dirichlet dispersion

20,000 draws, `scripts/ranking_robustness.py`, Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1,
NumPy PCG64 seeded `np.random.default_rng(20261007)`.

1. **Joint sampling on the 5-component simplex.** `f_prebotc` and `f_forebrain` are drawn from
   five-dimensional Dirichlet distributions centred on their nominal vectors,
   **f ~ Dirichlet(α₀ · p_nominal)** with total concentration **α₀ = 15.0** (`KAPPA`, which is
   itself **UNSOURCED**). This guarantees Σ_s f_s = 1 on the 4-simplex in every realisation.
2. **Component dispersion.** Each marginal follows Beta(α₀p_i, α₀(1−p_i)), so
   Var(X_i) = p_i(1−p_i)/(α₀+1) and CV_i = √[(1−p_i)/(p_i(α₀+1))].

   **Table 6B. Implied dispersion at α₀ = 15.0.**

   | Component | nominal p | Var(X) | implied CV |
   |---|---|---|---|
   | preBötC α1 | 0.60 | 0.01500 | 20.4% |
   | preBötC α2/3 | 0.15 | 0.00797 | 59.5% |
   | preBötC δ/α4 | 0.15 | 0.00797 | 59.5% |
   | preBötC ε | 0.08 | 0.00460 | 84.8% |
   | **preBötC α5** | **0.02** | **0.00123** | **175.0%** |
   | Forebrain α1 | 0.35 | 0.01422 | 34.1% |
   | Forebrain α5 | 0.30 | 0.01313 | 38.2% |
   | Forebrain α2/3 | 0.25 | 0.01172 | 43.3% |
   | Forebrain δ/α4 | 0.08 | 0.00460 | 84.8% |
   | Forebrain ε | 0.02 | 0.00123 | 175.0% |

   Describing this prior as having a uniform "CV ≈ 15%" would be wrong: α₀ = 15.0 is a global
   concentration parameter, and the component CV scales as √[(1−p)/(16p)]. For the rare preBötC
   α5 fraction (p = 0.02) the implied CV is **175%**, which is exactly why an unconstrained
   Dirichlet places **38% of draws below 0.002** — a tenth of nominal. This is the mechanism
   behind §4.4.3.
3. **Localisation and gain sampling.** For α1, α2/3 and α5, `f_extra` is drawn from Beta
   distributions centred on nominal with concentration 8.0 (a = 8m, b = 8(1−m)). δ/α4 is fixed
   at 1.00. For ε, the **extrasynaptic localisation fraction** is drawn
   f_extra,ε ~ U[0,1] — because the code labels it a GUESS and it should be treated as one —
   while ε's *regional abundance* remains governed by the Dirichlet simplex. ρ ~ U[3.0, 7.5].
4. **Non-finite draws.** 10 of 20,000 draws are discarded for the ideal α5 arm because its
   burden underflows to zero. These are the draws **most favourable** to that arm, so discarding
   them makes the reported figures conservative, and the draw count is stated gross.

#### 4.4.3. Prior sensitivity and tail separation

To test whether the selectivity advantage depends on draws placing preBötC α5 near zero, we
examined floor-constrained variants using a genuinely floor-preserving transformation:

q₅ = max(f₅, L),  f′₅ = q₅,  f′_j = f_j · (1 − q₅)/(1 − f₅)  for j ≠ 5

which guarantees f′₅ ≥ L and preserves Σ_j f′_j = 1 without shrinking the floored component
during rescaling.

**1. The lower tail is stable across every floor.**

| Floor L | Alogabat 5th | MP-III-022 5th | Ideal α5 5th |
|---|---|---|---|
| 0.0 (unconstrained) | **2.51×** | **2.49×** | **2.57×** |
| 0.01 (half-nominal) | **2.51×** | **2.49×** | **2.57×** |
| 0.02 (nominal) | **2.51×** | **2.48×** | **2.55×** |

P(R > 1) = 99.9% for all three α5-selective arms in every floor scenario. **This is the
quotable result**, and it says the ordering survives adverse draws of every invented parameter
simultaneously. It applies to these α5-selective profiles only: the neurosteroid arm has
P(R > 1) = 0.0% (worse than diazepam in 100% of draws, driven by high gating efficacy at
medullary α1), and HZ-166 reaches only 63.8%.

**2. The upper tail measures the prior, and is reported here only as a diagnostic.**

| Floor L | Alogabat median | Alogabat 95th | MP-III-022 median | MP-III-022 95th | Ideal median | Ideal 95th |
|---|---|---|---|---|---|---|
| 0.0 | 16.05× | 84.85× | 12.59× | 60.98× | 40.88× | 218.65× |
| 0.01 | 10.79× | 23.42× | 9.35× | 21.16× | 16.44× | 34.18× |
| 0.02 | 7.57× | 14.31× | 6.92× | 13.44× | 9.52× | 17.98× |

Because preBötC burden approaches zero as preBötC α5 vanishes, and 38% of unconstrained draws
place it below a tenth of nominal (§4.4.2), R diverges in the lower tail of the α5 prior,
skewing the median and 95th percentile upward. Flooring at nominal moves alogabat's median
16.05 → 7.57 and its 95th 84.85 → 14.31 — a factor of 5.9 on the upper tail from a single
prior assumption. **These five-fold swings are the reason the median and 95th percentile are
not reported as results**, and the reason `ranking_robustness.py` prints NOT QUOTABLE beside
them.

**3. Pairwise comparison.** Draw by draw, alogabat exceeds HZ-166 with median margin **+14.56**,
90% CI **[−0.09, +81.19]**, and α5 better in **94.8%** of draws. The confidence interval
touching zero is the honest summary of how close the α2/3 class can come.

#### 4.4.4. Two claims that must not be conflated

The results in §4.4.3 and §4.6 establish one thing and not another, and the distinction is the
single most important caveat in this manuscript.

**Claim A — robust within the specified model.** The ordering of subtype-selective arms persists
across 20,000 adverse parameter draws (§4.4.3) and across a complete replacement of the neuron
model (§4.6, Spearman ρ = +1.0000). This is a real and non-trivial property: it means the
ranking is not an artefact of any single invented number, nor of the integrate-and-fire
abstraction. It is a statement about the model's internal stability.

**Claim B — biologically validated.** That the ordering reflects actual compound effects in
tissue or in vivo. **This manuscript provides no evidence for Claim B, and Claim A does not
support it.** A model can be perfectly stable around the wrong inputs; stability measures
insensitivity, not correctness. Concretely, every draw in §4.4.3 samples *around* `f_α5,preBötC =
0.02` and `f_extra,α5 = 0.80`, both unsourced (§4.2) — so the Dirichlet sweep demonstrates
insensitivity to *dispersion* about those values while remaining fully dependent on their
*central* values being roughly right. Flooring preBötC α5 at nominal moves alogabat's upper tail
by a factor of 5.9 (§4.4.3), which is what that dependence looks like when made visible.

The reference-arm sign reversal (§4.6) makes the distinction concrete rather than theoretical.
Diazepam is the denominator of every R here. On the more biophysically detailed of our two
substrates its effect on respiratory output reverses sign. A quantity whose denominator behaves
that way across two plausible implementations of the same biology is not a measurement of
anything, however stable its ordering.

**What would move a result from A to B:** the patch-clamp protocol in §5.2 for the receptor
layer, and quantitative preBötC subunit proteomics plus per-molecule intrinsic-efficacy
measurement for the circuit layer. Neither exists here.

### 4.5. Why systems models require empirical efficacy: stereochemical invariance

Two-dimensional graph representations and topological fingerprints are widely used to predict
candidate properties. We tested whether they can distinguish stereoisomers with divergent
subtype selectivities, using SH-053-2′F-R-CH₃ (α5-selective, nominal R = **9.3**) and
SH-053-2′F-S-CH₃ (α2/α3/α5, nominal R = **3.6**).

**Table 7. Stereochemical invariance across enantiomers.**

| Descriptor | (R)-enantiomer | (S)-enantiomer | Δ |
|---|---|---|---|
| Molecular weight | 388.42 g/mol | 388.42 g/mol | 0.00 |
| Topological polar surface area | 61.86 Å² | 61.86 Å² | 0.00 |
| MolLogP | 2.841 | 2.841 | 0.00 |
| Morgan fingerprint (2048 bit) | bit-identical | bit-identical | 0 bits |
| Internal distance matrix | exact invariant | exact invariant | < 1e-9 Å |
| Vacuum MMFF energy | 42.184 kcal/mol | 42.184 kcal/mol | < 1e-6 kcal/mol |
| **Selectivity index R** (model-internal) | **9.3** | **3.6** | **+5.7** |

Enantiomers are related by an improper rotation, so all internal pairwise atomic distances and
2D graph invariants are mathematically identical. Three-dimensional representations capture
geometry, but **these 2D descriptors cannot by themselves encode the stereochemical information
responsible for the pharmacological difference without explicitly modelling the chiral binding
pocket.** Multiscale systems pharmacology must therefore ingest empirical, subtype-specific
efficacy measurements rather than infer them from 2D structural proxies.

Two things this result is not. The R values in the last row are **consequences** of the assigned
efficacy vectors, not predictions from structure — the model was told the two enantiomers differ
and reported the arithmetic. The argument is therefore about what 2D descriptors *cannot* do, and
it is a mathematical argument about improper rotations, not an empirical finding about these two
compounds. It stands independently of whether the assigned efficacy vectors are right.

### 4.6. Substrate independence: the ordering survives changing the neuron model

The index in §4.4 is algebraic, and the Dirichlet sweep perturbs its *parameters*. A stronger
test perturbs its *structure*. We therefore re-ran the respiratory arm on two neuron models
that share almost no mechanism:

* **LIF** — leaky integrate-and-fire with spike-triggered adaptation (C = 200 pF, g_L = 10 nS).
* **Conductance** — Butera–Rinzel–Smith 1999 model 1 (C = 21 pF, g_L = 2.8 nS; Butera, Rinzel &
  Smith, 1999a), with a real persistent sodium current I_NaP, voltage-dependent inactivation
  (τ_h = 10 s), an 11.8-fold Mg²⁺-relief span, and reachable depolarisation block.

Migrating between them required re-deriving four measurement protocols, not rescaling one
parameter: a single multiplicative weight scale cannot carry a weight table between the cells
(AMPA/GABA/glycine scale by the g_L ratio 0.28, NMDA additionally by a Mg²⁺-relief ratio, an
11-fold further correction), and each substrate must be scored against **its own** drug-free
control. The conductance substrate is anchored to an in vitro validity band
(`INVITRO_BAND` 0.05–1.00 Hz; neonatal rat preBötC slice, Revill et al., 2021), not the in vivo
eupnoea band, so its respiratory numbers are in vitro slice claims.

**Table 8. Fractional reduction in mean inspiratory output from each substrate's own control**
(`scripts/compare_substrates.py`, 3 seeds, occupancy 1.0).

| Arm | LIF | Conductance | Δ |
|---|---|---|---|
| Neurosteroid | +0.610 | +0.994 *(network silent, 0/3 seeds alive)* | +0.384 |
| MP-III-022 | +0.020 | +0.071 | +0.051 |
| Ideal α5 | +0.005 | +0.051 | +0.046 |
| Alogabat | +0.002 | +0.043 | +0.041 |
| HZ-166 | +0.001 | −0.100 | −0.101 |
| **Non-selective BZ** | **+0.079** | **−0.377** | **−0.456** |

**1. The subtype-selective ordering is substrate-independent.** Spearman **ρ = +1.0000** over
the five subtype-selective arms: the order
`neurosteroid > MP-III-022 > ideal α5 > alogabat > HZ-166` is identical on both substrates. Two
neuron models differing in every mechanism — real I_NaP, voltage-gated inactivation, Mg²⁺
relief, depolarisation block, versus none of those — give the same ranking. **This is a stronger
robustness argument than the Dirichlet sweep**, because it perturbs the model's structure rather
than its parameters.

**2. The non-selective reference arm is not substrate-independent — its sign flips.** Diazepam
goes from second-most-burdensome on the LIF (+0.079) to *least* on the conductance cell
(−0.377, i.e. a 38% **increase** in output). It is the arm with the largest α1 efficacy acting
on a preBötC the model treats as α1-predominant (`f_prebotc,α1 = 0.60`), so the flip is largest
exactly where the drug effect is largest.

Confirmed by two independent methods (FFT and direct burst counting over 180 s / 90 s):

| | Conductance | LIF |
|---|---|---|
| mean output | **+41.7%** | **−8.6%** |
| duty cycle | **+65.5%** | −6.7% |
| peak burst height | −7.3% | −1.4% |

So the effect is a change in *how much of the time* the network is active, with burst height
essentially unchanged on both substrates.

**3. The change is in burst DURATION, not burst rate.** This was previously unresolved,
because the only burst detector available reported 5.63 Hz against the FFT's 0.302 Hz — a
19-fold disagreement traceable to a single threshold with no hysteresis and no minimum
inter-burst interval, which re-triggers on intra-burst ripple. `circuitpharm.bursts` replaces
it with a Schmitt trigger, a minimum silent interval derived from the preparation's own
validity band, and an FFT cross-check returned with every result; on synthetic traces of known
burst count it recovers the true frequency to within 1.5% where the previous detector erred by
15–90×.

Mean output above the inter-burst floor decomposes exactly in logs, since
mean − floor = elevation × duration × frequency:

| term | Δlog (4 seeds, ± SD) | factor | significance |
|---|---|---|---|
| **burst duration** | **+0.5437 ± 0.0861** | **×1.72** | **6.3 σ** |
| burst elevation | −0.1218 ± 0.0402 | ×0.885 | 3.0 σ |
| burst frequency | −0.0571 ± 0.1218 | ×0.944 | **0.5 σ — not resolved** |
| *sum* | *+0.3648* | | |
| measured Δlog(mean − floor) | +0.3563 ± 0.0703 | ×1.43 | residual **2.4%** |

So on the conductance substrate the non-selective benzodiazepine **prolongs each inspiratory
burst by ~72% while leaving burst frequency unchanged within seed scatter**, with each burst
slightly weaker per unit time. The duration term is 6.3 times its own seed SD and nearly ten
times the frequency term, whose own scatter exceeds its mean — so "fewer bursts" is *not*
established, and the duty-cycle increase is a duration effect alone. The decomposition closes
to 2.4%, so the terms may be read as stated.

This is mechanistically coherent with the substrate difference. Burst termination on the Butera
cell depends on persistent-sodium inactivation accumulating during the burst; added GABA-A
shunting conductance slows depolarisation, so inactivation accumulates more slowly and the
burst runs longer. The LIF cell has no such mechanism — its bursts terminate on spike-triggered
adaptation with a 400 ms time constant — and correspondingly shows no duration increase.

On the LIF substrate the same arm gives Δlog(mean − floor) = **−0.1248 ± 0.0054**, distributed
across all three terms (elevation −0.0746 ± 0.0039, duration −0.0554 ± 0.0097, frequency
+0.0478 ± 0.0067) with a 34% residual. No term dominates, so that effect is reported as
**mixed** and no single-mechanism account of it is offered.

**Consistency with a retracted claim.** The frequency term being unresolved is what the
retraction in point 4 below predicted: the frequency differences were at or below one FFT bin.
Two independent measurements — a spectral peak and a time-domain burst count — now agree that
frequency does not move, which is the outcome that makes the retraction right rather than
merely cautious.

**4. A mechanism we retracted.** We initially attributed the flip to GABA-A shunting shortening
bursts → less I_NaP inactivation → faster recovery → higher rate, citing frequency readings as
confirmation. Those readings were 0.302 → 0.336 Hz on the conductance substrate (a difference
of **1.02 FFT bins**) and 1.327 → 1.327 Hz on the LIF (**below one bin**, whose width is 0.102 Hz
= 7.7%). Both are at the resolution limit, so neither supports the mechanism and the LIF's
"0.0%" means only "under 7.7%". The mechanistic claim is withdrawn; the mean, duty-cycle and
peak results above are robust and are what the substrate comparison establishes.

**Consequence for §4.4.** Diazepam is the denominator of every R in Table 6. On the more
biophysically detailed of our two substrates, its sign on respiratory output reverses. The
subtype-selective *ordering* is unaffected — that is result 1 — but any statement normalised to
the diazepam arm inherits this instability, and the absolute magnitudes in Table 6 should be
read accordingly.

---

## 5. Discussion

### 5.1. Evaluating the fixed scalar cap

The central conclusion is that **allosteric potentiation is state- and compartment-dependent.**
Representing PAM action as a fixed ≈ 2.5× scalar cap conflates synaptic receptor saturation with
extrasynaptic headroom. In synapses, agonist exposure restricts potentiation (**1.02×** peak
gain, **1.03×** asymptotic ceiling); in extrasynaptic compartments, submicromolar ambient GABA
leaves a wide **asymptotic receptor open-probability dynamic range** (**> 30×** across
0.2–0.8 µM, **123.3×** at 0.40 µM). A 2.5× cap binds only above ≈ **4 µM** ambient GABA, an
order of magnitude above the extrasynaptic regime — so there the cap truncates rather than
saturates.

**The unit of all of these numbers is open probability.** The 123.3× figure is the ratio of a
model-specific limiting open probability to the baseline open probability at a stated ambient
concentration. Converting it into a statement about current requires receptor density, single-
channel conductance and driving force; into a statement about inhibition, chloride homeostasis;
into a statement about circuit output, input resistance and network context (§2.3). **No step of
that chain is measured here, and the asymptote is in any case unreachable: the realised gain at
the modelled benchmark potency is 10.29×.** What the analysis establishes is that the *available*
headroom in the extrasynaptic pool is two orders of magnitude larger than in the synaptic one,
and therefore that one multiplier cannot describe both — not that any compound realises it.

Two implications for α5-targeted development, both conditional on the receptor layer alone and
so resting on the manuscript's firmest results:

1. **Screen for bounded s_max.** Overdose safety cannot be assumed from the "PAM mechanism"
   alone. For predominantly extrasynaptic targets the model predicts substantially greater tonic
   headroom above s_max ≈ 1.3–1.5; this warrants **experimental evaluation as a candidate
   efficacy boundary, not interpretation as an established safety threshold.** The supralinearity
   is the reason it matters: s_max 2.50 → 10.29× tonic gain, s_max 3.71 → 20.64×,
   s_max 5.73 → 59.22×.
2. **Separate synaptic and tonic conductances.** Computational neural models must decouple
   synaptic deactivation (τ_IPSC, **1.66×** here) from standing extrasynaptic conductance
   (**10.29×** here) rather than applying one lumped multiplier to all inhibitory inputs. The
   failure mode that follows from not doing so is documented in Supplementary Note §S2.2, where
   it produced an apparently reassuring result — universal survival at every dose — rather than
   an error.

What this section does **not** support is any claim about α5-selective compounds and respiration.
That argument runs through §4.2, whose two load-bearing anatomical numbers have no source (§7).

### 5.2. Proposed in vitro electrophysiological validation

The predictions are testable in recombinant-receptor electrophysiology without whole-animal
integration.

```
Expression system   Recombinant HEK293, human alpha5-beta3-gamma2
Recording           Whole-cell voltage clamp, V_hold = -70 mV
Agonist titration   [GABA]_bath in {0.1, 0.4, 1.0, 3.0, 10.0} uM
Compound            Concentration-response of test PAM (alogabat or MP-III-022),
                    0.1 nM to 10 uM, at each fixed ambient GABA baseline

Primary endpoint    The steady-state potentiation ratio surface
                    R_PAM([GABA]_bath, C_PAM)
                      = I_PAM([GABA]_bath, C_PAM) / I_baseline([GABA]_bath)
```

**Pre-specified model comparison.**

* **Model A (fixed scalar cap)** — a deliberately simple comparator:
  R_PAM ≤ 2.5 at every [GABA]_bath.
* **Model B (Markov kinetic model)** — under the nominal parameterisation at
  the cited tag, for a high-efficacy affinity PAM (s_max 2.40–2.50):

  | [GABA]_bath | R_PAM predicted (s_max 2.40 – 2.50) | asymptotic ceiling (k_off → 0⁺) |
  |---|---|---|
  | 0.1 µM | **10.36 – 11.71×** | 1910.8× |
  | 0.4 µM | **9.22 – 10.29×** | 123.3× |
  | 1.0 µM | **6.58 – 7.10×** | 21.5× |
  | 3.0 µM | **2.67 – 2.74×** | 3.6× |
  | 10.0 µM | **1.26 – 1.27×** | 1.3× |

  *(c_affinity = 3.2566 and 3.4674 respectively.)*

**The asymptote column is unreachable by any finite modulator and must not be used to set
criteria.** A criterion such as "R_max > 15× at 0.1 µM" is an asymptote reading: the realisable
prediction there is ≈ 8×, so a laboratory measuring 8× would report the model falsified when the
model in fact predicts 8×. This is not a hypothetical — an earlier version of this protocol
carried exactly that criterion, and two of its four intervals were violated by the model they
described (Supplementary Note §S3). The structural reason is in §4.1: as ambient GABA falls,
asymptotic headroom grows without limit while reachable gain barely moves. Both
columns are printed side by side by `scripts/paper_numbers.py --section falsification`, and a
regression test now asserts the divergence so that any future criterion written against the
wrong column fails loudly.

**Falsification rules** (pre-registered against the reachable column):

1. If a high-efficacy α5 PAM (s_max ≥ 2.4) shows maximal steady-state **R_PAM ≤ 2.5 at
   [GABA]_bath = 0.40 µM** in recombinant α5β3γ2 channels, the nominal Model-B parameterisation
   is falsified under the tested construct and conditions. (Model B predicts 9.22–10.29× there,
   comfortably clear of the Model-A ceiling.)
2. If R_PAM does **not decline monotonically** with increasing [GABA]_bath across the five
   tested concentrations, the affinity-modulation mechanism as formulated is wrong, independently
   of the parameterisation — the monotone decline follows from the topology, not the fit.
3. If measured R_PAM at 0.1 µM **exceeds** the 0.4 µM value by more than ≈ 1.2-fold, the
   finite-s_max ceiling is not the operative constraint and a gating component should be
   suspected (§2.6).

Rule 2 is the most valuable of the three: it tests the mechanism rather than the numbers, and
it cannot be rescued by re-fitting.

### 5.3. Limitations and modelling assumptions

Ordered by how much they constrain the conclusions, worst first.

**The two numbers that carry the circuit argument have no source.**

* **`f_α5,preBötC` = 0.02 is UNSOURCED** — the respiratory-sparing argument rests on it entirely.
* **`f_extra,α5` = 0.80 is UNSOURCED** — the extrasynaptic-headroom argument rests on it
  entirely. It was downgraded from FROM_QUALITATIVE while preparing this version, after the
  candidate source was read and found to have measured α1, α2 and β3 rather than α5 (§4.2, §7).
* Eight of fifteen Table 5 cells are UNSOURCED, a ninth is an explicit GUESS, and **none is
  QUANTITATIVE**. `KAPPA = 15.0`, which governs the entire uncertainty model, is UNSOURCED.
* **Six load-bearing citations do not support the specific numbers attached to them** (§7),
  including the diazepam benchmark that is the denominator of every selectivity ratio and the
  τ_IPSC anchor of the kinetic fit.
* Replacing these point values with sourced *ranges* and re-running the robustness analysis is
  the highest-value outstanding work in the project. Until then the §4.2–§4.4 results are
  conditional statements about a model.

**The reference arm is unstable across implementations.**

* The non-selective benzodiazepine **reverses sign** on respiratory output between the LIF and
  conductance substrates (§4.6), and it is the denominator of every R in Table 6. The
  subtype-selective *ordering* is unaffected; the magnitudes inherit the instability.
* The sign flip is a **burst-duration** effect (×1.72, 6.3 σ over 4 seeds) with burst
  frequency unresolved (0.5 σ), established with a substrate-matched detector (§4.6). What
  remains open is the LIF arm, where no term dominates and the log decomposition leaves a 34%
  residual, so no single-mechanism account of that effect is offered.
* A mechanism we initially proposed for the flip was **retracted on resolution grounds**: the
  supporting frequency differences were 1.02 FFT bins (conductance) and below one bin (LIF). The
  mean, duty-cycle and peak results survive; the mechanistic account does not.

**Receptor-layer assumptions.**

* **Three anchors do not identify six microscopic rates.** This is one admissible
  parameterisation consistent with the stated constraints, not an identified rate set.
* **Desensitisation is unconstrained by the anchors.** d/r = 25.0 is provisional and sets the
  asymptote directly: ±2-fold moves the dynamic range over 67.7–211.1× (§2.5). Every asymptotic
  figure here is conditional on it. This is the largest single source of uncertainty in the
  manuscript's firmest result.
* **τ_IPSC = 15.0 ms is a conventional neuronal IPSC value and is not supported by the kinetic
  source cited beside it** — Haas & Macdonald (1999) measured 76.1 ms deactivation for
  recombinant α1β3γ2L, about 5× slower (§7).
* **Affinity modulation acts through k_off only.** Ligands altering gating or desensitisation
  show different quantitative profiles (§2.6) — and Walters et al. (2000) indicates that
  diazepam itself has two separable components, so a single affinity parameter is an
  oversimplification for the reference compound (Table 2 note iii).
* **Ambient GABA** is evaluated at a representative 0.40 µM; in vivo, GAT-1/GAT-3 transport and
  spillover create dynamic microdomains. The 0.2–0.8 µM span was not verified against its
  source's full text (§7).

**Circuit-layer assumptions.**

* **The A4 riluzole dissociation is NOT reproduced** (0 of 3 seeds). This is *inherited* rather
  than a tuning failure: Butera–Rinzel–Smith model 1 **is** the pacemaker hypothesis, and
  riluzole is the principal published argument against it. Declared, not retuned.
* **The coupling weights are ours, not published.** Butera et al. (1999b), the population paper,
  was unreachable (journal HTTP 403; every available machine-readable encoding is single-cell).
* **Respiratory calibration remains wet-lab blocked**, and the conductance migration moved the
  > P12 muscimol anchor *further* away by importing neonatal parameters.
* **The conductance substrate is anchored to an in vitro band** (0.05–1.00 Hz, neonatal slice),
  so its respiratory outputs are slice claims rather than ventilatory ones; the frequency range
  was not verified against its source's abstract (§7).
* Four measurement protocols inherited from the LIF substrate required re-derivation for the
  conductance cell — warm-up against τ_h = 10 s, validity band, FFT resolution, burst-detection
  threshold. **None failed loudly.** Working rule: on a new substrate, assume every inherited
  protocol is wrong until re-derived.

**Scope.**

* **Pharmacokinetics are absent.** Receptor occupancy and network dynamics only; no absorption,
  blood-brain barrier penetration or clearance, which determine clinical dose-occupancy
  relationships.
* **No absolute margin is computed**, and `Tier.VOID` blocks lethal-dose margins
  programmatically. R orders candidates; it cannot bound risk.
* **Robustness within the model is not biological validation** (§4.4.4). Every robustness result
  here is of the first kind.
* The work stays at mechanism and compound-class level and reports no formulation or dosing.

## 6. Conclusions

**The receptor-kinetic result.** A five-state continuous-time Markov gating scheme shows that
allosteric potentiation headroom is state- and compartment-dependent. The classical ceiling binds
tightly in the synapse — **1.02×** peak gain, **1.03×** asymptotic ceiling, with the drug's action
redirected into the deactivation tail at **1.66×** τ and **1.66×** integrated open probability —
but expands to an **asymptotic receptor open-probability dynamic range of 123.3×** in
extrasynaptic microenvironments at submicromolar ambient GABA, a **120-fold** contrast between
the two compartments' ceilings. A fixed 2.5× scalar cap binds only above ≈ 4 µM ambient GABA, so
in the extrasynaptic regime it truncates rather than saturates. These are statements about open
probability within a stated topology and parameterisation; they are reproducible from the cited
tag, and they are the manuscript's defensible core.

**The 123.3× figure is not a drug effect.** It is the ratio of a limiting open probability to a
baseline open probability, both internal to this scheme. It is not a predicted fold change in
whole-cell current, in tonic inhibition, in circuit output, or in any clinical endpoint, and no
finite modulator approaches it — at the modelled benchmark potency the realised tonic gain is
**10.29×**, and asymptote and realised gain diverge further as ambient GABA falls.

**The circuit-level implications are conditional, and two of their load-bearing inputs are
unsourced.** Under the model's assumptions, low α5 representation in the preBötzinger complex
reduces direct inhibition of the modelled rhythm generator while high extrasynaptic localisation
exposes forebrain circuits to large tonic conductance increases during dose escalation. Both the
preBötC α5 fraction and the α5 extrasynaptic fraction are **UNSOURCED** (§4.2, §7) — the second
downgraded in preparing this version, after the source that had been carrying it was read and
found to have measured neighbouring subunits. **This manuscript therefore provides no evidence
that α5-selective compounds preserve respiration or carry a safety advantage**, and the ordering
in Table 6 should be read as a model-internal index rather than a property of these molecules.

**Robustness within a model is not validation.** The compound ordering survives 20,000 adverse
parameter draws and survives replacing the neuron model entirely (Spearman ρ = +1.0000 across an
integrate-and-fire and a conductance-based substrate). That is a real property of the model and a
non-trivial one — the ranking is not an artefact of any single invented number nor of the
integrate-and-fire abstraction. It is not evidence about tissue (§4.4.4). The sweep samples
*around* two unsourced central values, and the non-selective benzodiazepine reference arm — the
denominator of every reported ratio — **reverses sign** on respiratory output between the two
substrates.

**What would settle it.** The patch-clamp protocol of §5.2 tests the compartment asymmetry
directly, in recombinant receptors, without whole-animal integration; its pre-registered
intervals are stated against the realisable gain rather than the asymptote, and one of its three
rules tests the mechanism rather than the parameterisation and so cannot be rescued by refitting.
For the circuit layer, quantitative preBötC subunit measurement and per-molecule intrinsic
efficacy are what the unsourced numbers need. Overdose-relevant conductance headroom depends
jointly on ligand efficacy, receptor mechanism, ambient agonist concentration,
compartmentalisation and network context — the quantity a fixed scalar cap stood in for is
neither a constant nor small. Establishing what it is in tissue is an experiment, not a
simulation.

---

## 7. Citation claim-support audit

Metadata verification and claim support are different checks, and this project has now been
caught twice by the gap between them (§4.2). A source can resolve perfectly by DOI, sit in the
right journal, carry a title that matches the claim, and measure a neighbouring quantity.

Every reference was resolved against Crossref by title and by journal/volume/page. For the
load-bearing ones, abstracts were then retrieved from Europe PMC and read against the specific
statement each is attached to. **Six load-bearing citations do not support the specific number
attached to them**, and the corresponding parameters are relabelled UNSOURCED in the text. They
fail in three distinguishable ways: three measure a *different quantity* than the claim requires
(Kasugai, Walters, `a5_dist`), two report *no fold shift at all* of the kind the mapping needs
(Stamenić, Cecere), and one measures the right quantity at a value ~5× from the anchor beside it
(Haas & Macdonald).

| Citation | Statement attached to it | Metadata | Claim support |
|---|---|---|---|
| Kasugai et al. 2010 | `f_extra,α5 = 0.80` (α5 > 80% extrasynaptic) | resolved | **DOES NOT SUPPORT** — measured α1/α2/β3, not α5; found synaptic density 78–132× extrasynaptic. Parameter → UNSOURCED (§4.2) |
| Walters et al. 2000 | diazepam `ceiling = 2.50` | resolved | **DOES NOT SUPPORT** — reports biphasic nanomolar/micromolar potentiation, no EC₅₀ shift. Parameter → UNSOURCED; and the finding argues *against* a single affinity parameter (Table 2 note iii) |
| Stamenić et al. 2016 | MP-III-022 `ceiling = 2.50` | resolved | **DIRECTION ONLY** — confirms binding- and efficacy-selective α5 PAM, non-α5 engaged only at top dose; no EC₅₀ shift reported |
| Cecere et al. 2025 | alogabat `ceiling = 2.50`; "+167% rat / +72% human EC₂₀" | resolved | **PARTIAL** — abstract confirms a potent α5β3γ2 PAM with binding and functional selectivity; the two percentages are not in the abstract and were **not** verified against full text |
| Haas & Macdonald 1999 | kinetic topology; subunit composition sets deactivation | resolved | **SUPPORTS the topology claim.** Note their α1β3γ2L deactivation is 76.1 ms, ~5× our τ_IPSC = 15 ms anchor; the anchor is a conventional neuronal IPSC value and is **not** supported by this source |
| Farrant & Nusser 2005 | ambient GABA 0.2–0.8 µM | resolved | **PARTIAL** — abstract states "low concentrations of ambient GABA" without the range; the numeric span is presumably in the review body, unverified here |
| Revill et al. 2021 | `INVITRO_BAND` 0.05–1.00 Hz | resolved | **PARTIAL** — confirms the preparation (neonatal rat slices retaining respiratory rhythmicity); no frequency range in the abstract |
| Liu & Wong-Riley 2004 (`pbc_alpha`) | preBötC α1/α2/α3 fractions | resolved by DOI | **UNASSESSED** — full text HTTP 403; covers α1/α2/α3 only; developmental study |
| `a5_disc` (2005, PMID 15650112) | `w_subj,α5 = 1.0` | resolved | **SUPPORTS direction only** — α5 agonists mimic ethanol's discriminative stimulus, α5 inverse agonist blocks it. Supplies no weight |
| `a5_dist` (1988, PMID 2844998) | regional subunit composition | resolved, similarity 1.00 | **DOES NOT SUPPORT** — single generic α-subunit probe; measures regional *level*, not *composition* |
| Jones & Westbrook 1995 | desensitised states prolong brief-pulse responses | resolved | **UNASSESSED** — not indexed in Europe PMC; topology claim, widely replicated |
| Otis & Mody 1992 | benzodiazepine prolongs IPSC decay | resolved | **UNASSESSED** — not indexed in Europe PMC |
| Pirker et al. 2000 | qualitative ordering: high medullary α1, low medullary α5 | resolved | **UNASSESSED** — not indexed in Europe PMC; cited for ordering only, no number |
| Rudolph 1999; Löw 2000; McKernan 2000; Cheng 2006; Olsen & Sieghart 2008; Rudolph & Möhler 2014; Nutt 2006; Rundfeldt & Löscher 2014; Atack 2011; Butera et al. 1999a/b | subtype-endpoint dissociation, review and background claims | all resolved | **not individually claim-audited** — background rather than load-bearing |

**What this audit changes.** Of the compounds whose `ceiling` values appear in Table 6, not one
has that value supported by a cited source. Of the two anatomical numbers carrying the
respiratory and headroom arguments, neither has a supporting source. **The receptor-kinetic
results in §2 and §4.1 are unaffected** — they depend on the three stated calibration anchors
and the topology, not on these citations — which is why the Scope box places them in a different
tier from everything downstream.

**What it does not change.** None of the four failures means the underlying statement is false.
α5 probably *is* enriched extrasynaptically; diazepam probably *does* shift GABA EC₅₀ by
something in the 2–3× range. What the audit establishes is that **this manuscript cannot cite a
source for those numbers**, and so must label them as assumptions. That is a weaker position than
a resolved reference list alone would imply, and a more accurate one.

---

## Supplementary Note

### §S1. Reproducing the numbers

```
git clone https://github.com/willkhinz/circuitpharm && cd circuitpharm
git checkout manuscript-v3
pip install -e '.[dev]'
python scripts/paper_numbers.py                       # all figures in sections 2, 4.1, 4.3, 5.2
python scripts/ranking_robustness.py --draws 20000    # section 4.4.3, add --a5-floor 0.01 / 0.02
python scripts/compare_substrates.py                  # section 4.6
python scripts/verify_manuscript.py                   # checks this document against a clean
                                                      # checkout of the tag it cites
pytest tests/test_manuscript_consistency.py           # checks it against the working tree
```

`verify_manuscript.py` is the stronger of the two checks: it creates a throwaway git worktree at
the cited tag, runs the generators with that tree's `src` on the path, and asserts every audited
figure appears in output produced by code it did not write. It also verifies that the cited tag
*contains* this document and its tests — a check added because the first tag created for this
purpose did not, pointing instead at the commit immediately before the manuscript was written.
The figures reproduced; the citation was still wrong, and nothing else caught it.

Environment at the cited tag: Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1, NumPy PCG64 seeded
`np.random.default_rng(20261007)`.

### §S2. Two modelling failure modes, with their diagnostic histories

**S2.1. The fixed scalar cap.** Early implementations scaled PAM gain linearly with occupancy and
clamped the result:

```python
# legacy formulation (REJECTED)
gaba_a_gain = min(1.0 + occupancy * (target_gain - 1.0), gaba_a_efficacy_cap)  # cap = 2.5
```

This manufactured a safety guarantee. Because conductance was clamped at 2.5×, circuits
inevitably survived dose escalation, and the resulting "overdose ceiling" was a programmer's
boundary rather than an emergent property of receptor saturation. The present framework lets the
boundary emerge from the generator matrix Q, and §2.5 locates it: a 2.5× ratio is reached only
above ≈ 4 µM ambient GABA, roughly an order of magnitude above the extrasynaptic regime.

**S2.2. Lumped-parameter transfer across architectures.** Moving from a single-pool to a dual-pool
conductance model, `scripts/overdose_kinetic.py` carried over a lumped `SENS_TOTAL = 0.10` fitted
to the single-pool version. The mismatch collapsed the modelled respiratory effect of a
non-selective benzodiazepine to **−2% ventilation** against a clinical anchor of −16% to −19%.
The visible artifact was not an error message: the script reported that **every compound survived
every dose up to 100% receptor occupancy**, which reads as a reassuring result. A parameter had
stopped meaning what it meant, and the output stayed plausible.

This failure mode has recurred twice more in the project since. A tonic/phasic split was applied
to every module except the one written last, where the defect reappeared and was found two
sessions later by review. More recently, a shared population factory's default silently removed
spike-triggered adaptation from the spinal circuit while leaving the respiratory circuit correct;
two phenotype tests caught it and nothing else did. The common shape is **a value duplicated
across two locations that then diverge, with a plausible rather than loud symptom.** The remedy
adopted is structural: delete the second copy rather than synchronise it — a shared helper now
takes no default for a value its callers own, and this manuscript's tables are printed by a
generator rather than typed.

**S2.3. Why this belongs in a methods note rather than the results.** The project's own
documentation treats these as recurring-error entries with regression tests attached
(`tests/test_review_regressions.py`). They are reported here because they bear directly on the
§5.1 recommendation — that models must decouple synaptic deactivation from standing extrasynaptic
conductance rather than applying one multiplier to both — and because a reader evaluating the
reliability tiers in §3 is entitled to know what produced them. They are not results about
GABA-A receptors.

---

### §S3. Corrections record

This manuscript supersedes an earlier draft. The corrections are recorded here rather than
in the main text, where they distracted from the argument, and because a reader checking a
citation is entitled to know it was changed.

**Reference coordinates corrected.** Crestani et al. — 2002 and pages 8980–8985, not 2001 and
8993–8997. Haas & Macdonald — *J. Physiol.* 514:27–45, not *J. Neurosci.* 19:2435–2445, and under
its actual title. Kasugai et al. — *Eur. J. Neurosci.* 32:1868–1888, not *J. Neurosci.*
30:14024–14035. Otis & Mody — *Neuroscience* 49:13–32; the cited *J. Physiol.* 454:477–496 does
not resolve. Walters et al. — *Nat. Neurosci.* 3:1274–1281, not *Br. J. Pharmacol.*
131:1307–1314, and under its actual title.

**Three citations withdrawn as unresolvable**, listed after the reference list with the
parameters they carried relabelled UNSOURCED. MP-III-022's primary characterisation is Stamenić
et al. (2016), ten years later than the withdrawn attribution.

**Parameterisation — a correction to a correction.** An intermediate revision of this
manuscript reported figures from `gabaa_kinetics`' module-default synaptic transient
(1000 µM / 0.30 ms) in the belief that the earlier draft's 3000 µM / 1.00 ms transient was
unsupported by any commit. That was wrong. `config.SYNAPTIC_PULSE` has held
(3000 µM, 1.00 ms) throughout; it is the **calibrated** pulse — one of the 9 of 27 (peak,
clearance) cells in which the scheme reproduced every anchored benzodiazepine observable from a
single free parameter — and every pharmacology consumer in the repository passes it explicitly.
The module defaults were a provisional literature estimate reached only by callers that pass no
pulse, which was the figure-generating script and nothing else.

So the constant had three definitions: `config.SYNAPTIC_PULSE`, an inline literal in
`cpg.Drug.from_kinetics` that matched it, and the `gabaa_kinetics` defaults that did not. It now
has one: `gabaa_kinetics` imports from `config`, the literal is gone, and
`tests/test_review_regressions.py` pins the single definition. Every figure in this manuscript
is generated on that unified calibration, which restores the original draft's values
(tonic headroom 123.3×, phasic peak gain 1.021×, c_affinity 3.4674, K_d 23.81 µM, E 3.9249).

The discrepancy surfaced only when `scripts/decompose_burst_change.py` constructed a `Drug`
through the circuit path and printed its tonic gain as 10.285 beside a manuscript asserting
7.876. One quantity, two values, and nothing had computed it twice until then. The earlier
draft's unresolvable **commit hash** was a real defect and is replaced by a tag that
`scripts/verify_manuscript.py` checks.

**Conflated quantities separated.** The earlier draft used one figure for the *simulated* peak
open probability at saturating agonist (0.750, the calibration anchor, with desensitisation
competing throughout the rise) and the *analytic* gating bound β/(α+β) (0.828, unattainable).
Table 1 now lists them separately.

**Falsification intervals re-derived.** The earlier draft's §5.2 pre-registered four intervals,
two of which were violated by the model they described — including "R_max > 15× at 0.1 µM" where
the realisable prediction is ≈ 8×. The intervals were set from the asymptote rather than the
reachable gain. They are now stated against the reachable gain with the asymptote printed
alongside, and two mechanism-level rules were added that cannot be rescued by refitting.

**Statistics withdrawn from headline presentation.** The earlier draft's selectivity table
included the Monte Carlo median and 95th percentile, which `ranking_robustness.py` marks NOT
QUOTABLE. They appear in §4.4.3 as a prior-sensitivity diagnostic only.

**Provenance downgraded.** `f_extra,α5 = 0.80` moved from FROM_QUALITATIVE to **UNSOURCED** after
a claim-support read of its candidate source (§4.2, §7). Six load-bearing citations were found
not to support the specific numbers attached to them.

---

## References

1. Atack, J. R. (2011). GABA-A receptor subtype-selective modulators. II. α5-selective inverse
   agonists for cognition enhancement. *Curr. Top. Med. Chem.*, 11(9), 1203–1214.
   doi:10.2174/156802611795371314
2. Butera, R. J., Rinzel, J., & Smith, J. C. (1999a). Models of respiratory rhythm generation in
   the pre-Bötzinger complex. I. Bursting pacemaker neurons. *J. Neurophysiol.*, 82(1), 382–397.
   doi:10.1152/jn.1999.82.1.382
3. Butera, R. J., Rinzel, J., & Smith, J. C. (1999b). Models of respiratory rhythm generation in
   the pre-Bötzinger complex. II. Populations of coupled pacemaker neurons. *J. Neurophysiol.*,
   82(1), 398–415. doi:10.1152/jn.1999.82.1.398 *(cited for completeness; full text was not
   accessible to this work — see §5.3)*
4. Cecere, G., Ballard, T. M., Knoflach, F., et al. (2025). Preclinical pharmacology of
   alogabat: a novel GABA-A-α5 positive allosteric modulator targeting neurodevelopmental
   disorders with impaired GABA-A signalling. *Front. Pharmacol.*, 16, 1626078.
   doi:10.3389/fphar.2025.1626078
5. Cheng, V. Y., Martin, L. J., Elliott, E. M., et al. (2006). α5GABA-A receptors mediate the
   amnestic but not sedative-hypnotic effects of the general anaesthetic etomidate.
   *J. Neurosci.*, 26(14), 3713–3720.
6. Crestani, F., Keist, R., Fritschy, J.-M., et al. (2002). Trace fear conditioning involves
   hippocampal α5 GABA-A receptors. *Proc. Natl. Acad. Sci. USA*, 99(13), 8980–8985.
7. Farrant, M., & Nusser, Z. (2005). Variations on an inhibitory theme: phasic and tonic
   activation of GABA-A receptors. *Nat. Rev. Neurosci.*, 6(3), 215–229.
8. Haas, K. F., & Macdonald, R. L. (1999). GABA-A receptor subunit γ2 and δ subtypes confer
   unique kinetic properties on recombinant GABA-A receptor currents in mouse fibroblasts.
   *J. Physiol.*, 514(1), 27–45.
9. Jones, M. V., & Westbrook, G. L. (1995). Desensitised states prolong GABA-A channel responses
   to brief agonist pulses. *Neuron*, 15(1), 181–191.
10. Kasugai, Y., Swinny, J. D., Roberts, J. D. B., et al. (2010). Quantitative localisation of
    synaptic and extrasynaptic GABA-A receptor subunits on hippocampal pyramidal cells by
    freeze-fracture replica immunolabelling. *Eur. J. Neurosci.*, 32(11), 1868–1888.
11. Liu, Q., & Wong-Riley, M. T. T. (2004). Developmental changes in the expression of GABA-A
    receptor subunits α1, α2 and α3 in the rat pre-Bötzinger complex. *J. Appl. Physiol.*,
    96(5), 1825–1831. doi:10.1152/japplphysiol.01264.2003 *(metadata verified; full text
    inaccessible, claim support UNASSESSED — see §4.2)*
12. Löw, K., Crestani, F., Keist, R., et al. (2000). Molecular and neuronal substrate for the
    selective attenuation of anxiety. *Science*, 290(5489), 131–134.
13. McKernan, R. M., Rosahl, T. W., Reynolds, D. S., et al. (2000). Sedative but not anxiolytic
    properties of benzodiazepines are mediated by the GABA-A receptor α1 subtype.
    *Nat. Neurosci.*, 3(6), 587–592.
14. Nutt, D. J. (2006). Alcohol alternatives — a goal for psychopharmacology?
    *J. Psychopharmacol.*, 20(3), 318–320.
15. Olsen, R. W., & Sieghart, W. (2008). International Union of Pharmacology. LXX. Subtypes of
    γ-aminobutyric acid-A receptors: classification on the basis of subunit composition,
    pharmacology and function. Update. *Pharmacol. Rev.*, 60(3), 243–260.
16. Otis, T. S., & Mody, I. (1992). Modulation of decay kinetics and frequency of GABA-A
    receptor-mediated spontaneous inhibitory postsynaptic currents in hippocampal neurons.
    *Neuroscience*, 49(1), 13–32.
17. Pirker, S., Schwarzer, C., Wieselthaler, A., Sieghart, W., & Sperk, G. (2000). GABA-A
    receptors: immunocytochemical distribution of 13 subunits in the adult rat brain.
    *Neuroscience*, 101(4), 815–850.
18. Revill, A. L., Katzell, A., Del Negro, C. A., Milsom, W. K., & Funk, G. D. (2021). KCNQ
    current contributes to inspiratory burst termination in the pre-Bötzinger complex of
    neonatal rats in vitro. *Front. Physiol.*, 12, 626470. doi:10.3389/fphys.2021.626470
19. Rudolph, U., Crestani, F., Benke, D., et al. (1999). Benzodiazepine actions mediated by
    specific γ-aminobutyric acid-A receptor subtypes. *Nature*, 401(6755), 796–800.
20. Rudolph, U., & Möhler, H. (2014). GABA-A receptor subtypes: therapeutic potential in Down
    syndrome, affective disorders, schizophrenia and autism. *Annu. Rev. Pharmacol. Toxicol.*,
    54, 483–507.
21. Rundfeldt, C., & Löscher, W. (2014). The pharmacology of imepitoin: the first partial
    benzodiazepine receptor agonist developed for the treatment of epilepsy. *CNS Drugs*,
    28(1), 29–43. *(published online 2013; 2014 issue)*
22. Stamenić, T. T., Poe, M. M., Rehman, S., et al. (2016). Ester to amide substitution improves
    selectivity, efficacy and kinetic behaviour of a benzodiazepine positive modulator of
    GABA-A receptors containing the α5 subunit. *Eur. J. Pharmacol.*, 791, 433–443.
    doi:10.1016/j.ejphar.2016.09.016
23. Walters, R. J., Hadley, S. H., Morris, K. D. W., & Amin, J. (2000). Benzodiazepines act on
    GABA-A receptors via two distinct and separable mechanisms. *Nat. Neurosci.*, 3(12),
    1274–1281.

### Withdrawn citations

The following were cited in an earlier version and could not be resolved against Crossref by
title or by journal/volume/page. Each carried a model parameter; rather than substitute a
plausible-looking replacement, the parameters they supported are labelled **UNSOURCED** in the
text above, matching what `circuitpharm.provenance` already records internally.

* "Nutt, D. J., et al. (2007). Development of α5-selective GABA-A modulators.
  *Neuropharmacology*, 53(7), 810–820." — no match. Replaced for the review function by
  Atack (2011), ref. 1.
* "Fischer, B. D., et al. (2010). Anxiolytic-like effects of MP-III-022…
  *Neuropharmacology*, 59(7–8), 612–618." — no match, and anachronistic: MP-III-022 was
  characterised in 2016. Replaced by Stamenić et al. (2016), ref. 22.
* "Saba, L. M., et al. (2017). The α5 subunit of the GABA-A receptor as a novel target for
  alcohol use disorders. *Alcohol. Clin. Exp. Res.*, 41(4), 748–758." — no match. It supported
  the α5-alcohol rationale (now carried by Nutt 2006, ref. 14, and `a5_disc`, §4.4.1) and
  `w_subj,δ = 0.0`, which is now declared UNSOURCED.
