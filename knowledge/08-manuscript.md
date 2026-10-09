# Compartment-Dependent Allosteric Headroom in GABA-A Receptors Revealed by Markov-Kinetic Modelling of Phasic and Tonic Inhibition

**William Hinz**¹
¹*Department of Biochemistry & Computational Neuropharmacology*
Correspondence: hinzwilliam52@gmail.com
Repository: https://github.com/willkhinz/circuitpharm — tag **`manuscript-v2`** (commit `1c0faa5`)

> **Reproducibility.** Every number in this manuscript is emitted by
> `scripts/paper_numbers.py` at the tagged commit, and the Monte Carlo by
> `scripts/ranking_robustness.py --draws 20000`. None was transcribed. The previous draft
> quoted a parameterisation no commit produces and cited a commit hash that exists in no ref;
> see `knowledge/07-paper-review.md` for the audit and `tests/test_review_regressions.py`
> (E21) for the regression tests that now prevent it.

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
**1.32×** (asymptotic ceiling **1.67×**), acting instead on the deactivation tail
(τ **1.90×**, time-integrated open probability **2.42×**). Submicromolar ambient GABA leaves
the steady state far below saturation: at 0.40 µM the same modulator produces a **7.88×**
increase in standing open probability, against an asymptotic bound of **184.6×** in the limit
of vanishing agonist dissociation (k_off → 0⁺). That bound is a model-specific
open-probability dynamic range, not a claim of equivalent whole-cell or circuit-level
potentiation.

The dynamic range is an emergent property of this topology and parameterisation, not a
pharmacological constant. A one-at-a-time sensitivity sweep shows the exact ratio depends on
the desensitisation equilibrium (d/r) and microscopic affinity (K_d), while a large
extrasynaptic range (**> 48×**) persists across the submicromolar ambient concentrations
reported for cortex and hippocampus (0.2–0.8 µM). Incorporating receptor compartmentalisation
and α5-subunit distribution, subtype-selective PAMs separate forebrain from respiratory
actions without any imposed cap; under the model's assumed low α5 representation in the
preBötzinger complex, α5-selective modulation produces less direct inhibition of the modelled
core respiratory rhythm generator — though whole-animal ventilatory preservation cannot be
inferred from receptor abundance alone.

Two results constrain how far this should be taken. First, the compound ordering is
**substrate-independent**: migrating the respiratory circuit from an integrate-and-fire cell to
a Butera–Rinzel–Smith conductance cell with a real persistent sodium current leaves the
ordering of five subtype-selective arms unchanged (Spearman ρ = +1.0000). Second, the
**non-selective benzodiazepine reference arm is not** substrate-independent: its sign on
respiratory output reverses between the two cell models, and that arm is the denominator of
every selectivity ratio reported here.

We identify two computational failure modes in multiscale systems models — arbitrary
conductance caps, and transfer of lumped sensitivity parameters between model architectures
(recurring error E12) — and propose a recombinant-receptor patch-clamp experiment measuring the
two-dimensional surface I_tonic([GABA]_bath, C_PAM) to test the predicted dependence of tonic
potentiation on ambient GABA.

**Keywords:** GABA-A receptor; α5 subunit; positive allosteric modulator; Markov kinetic model;
tonic inhibition; preBötzinger complex; open-probability dynamic range; model identifiability.

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
| Limiting open probability as k_off → 0⁺ | kinetic asymptote | P_open,∞ = **0.156377** |
| Asymptotic bound on open probability at 0.40 µM | asymptotic open-probability dynamic range | P_open,∞ / P_open,base = **184.6×** |
| Maximum operational potency shift | max operational potency-shift factor | s_max = EC₅₀,base / EC₅₀,drug |
| Compound conductance multiplier | finite PAM gain | tonic gain / phasic gain |
| Simulated peak P_o at saturating agonist | P_o,max (calibration anchor) | **0.7500** |
| Analytic gating bound, desensitisation absent | β/(α+β) | **0.8282** |

The last two rows are a distinction the previous draft collapsed. They are not the same
quantity: β/(α+β) is the equilibrium open probability of the two-site scheme *without* the
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
is proportional to synaptic charge, with 7.454 ms ≡ 0.007454 pC/pA.

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
| k_on | 0.0112168 µM⁻¹ms⁻¹ |
| k_off | 0.333062 ms⁻¹ |
| β | 0.646493 ms⁻¹ |
| α | 0.134142 ms⁻¹ |
| d | 0.050 ms⁻¹ |
| r | 0.0020 ms⁻¹ |

Reproduced anchors: EC₅₀ = **20.0000 µM**, P_o,max = **0.7500**, τ_IPSC = **15.000 ms** — all
three to the precision quoted.

### 2.2. Microscopic equilibrium versus functional potency

Keeping binding parameters separate from concentration-response endpoints:

* **Fitted targets:** EC₅₀ = 20.0 µM, P_o,max = 0.750, τ_IPSC = 15.0 ms.
* **Derived microscopic quantities:** K_d, gating efficacy E.
* **Model predictions:** P_open(0.40 µM) = 8.4704 × 10⁻⁴, asymptotic range 184.6×.

K_d = k_off/k_on = 0.333062 / 0.0112168 = **29.69 µM**.
E = β/α = 0.646493 / 0.134142 = **4.8195**; α/β = **0.20749**.

For the two-site scheme *without* desensitisation, equilibrium open probability is

P_open([G]) = E[G]² / (K_d² + 2K_d[G] + (1+E)[G]²)

so as [G] → ∞, P_open → E/(1+E) = β/(α+β) = **0.8282**. Half-maximal open probability satisfies
(1+E)[G]²₁/₂ − 2K_d[G]₁/₂ − K_d² = 0, whose positive root is

**[G]₁/₂ = K_d (1 + √(2+E)) / (1 + E)**

Substituting: [G]₁/₂ = 29.69 × 3.6114 / 5.8195 = **18.43 µM**. High gating efficacy
(E = 4.82 > 1) pulls the equilibrium forward, shifting the half-maximal concentration leftward
from K_d = 29.69 µM to 18.43 µM. When desensitisation (d/r = 25.0) is included and evaluated
numerically over a 300 ms application, the peak-current EC₅₀ is **20.000 µM**, reproducing the
calibration target — reconciling microscopic binding with macroscopic activation, and showing
that the analytic equilibrium midpoint and the simulated peak EC₅₀ differ by 8%.

### 2.3. Asymptotic steady-state open probability as k_off → 0⁺

An **affinity-type PAM** increases apparent agonist affinity by slowing unbinding by
c_affinity ≥ 1: k_off^drug = k_off^base / c_affinity.

We evaluate [G] = 0.40 µM as a representative ambient scenario within the range reported for
selected preparations (0.2–0.8 µM; Farrant & Nusser, 2005). Solving P·Q = 0:

**P_open,baseline(0.40 µM) = 8.4704 × 10⁻⁴** (0.0847%).

In the limit c_affinity → ∞ (k_off → 0⁺) at any [G] > 0, unbinding vanishes; R and RG become
transient and their stationary probabilities go to zero. All probability concentrates in the
closed sub-chain {RG₂, RG₂*, D}. Detailed balance on each pair gives
P_RG₂ = (α/β)·P_open and P_D = (d/r)·P_RG₂, so normalisation yields

**P_open,∞ = 1 / [ 1 + (α/β)(1 + d/r) ]**

Substituting α/β = 0.20749 and d/r = 25.0:

P_open,∞ = 1 / (1 + 0.20749 × 26.0) = 1 / 6.3947 = **0.156377** (15.64%).

Computed two independent ways — this closed form, and the mechanism pushed numerically to
c_affinity = 10⁵ — the two agree to **7.5 × 10⁻⁶**, and `paper_numbers.py` prints both so a
divergence would be visible. Dividing by baseline:

**P_open,∞ / P_open,baseline(0.40 µM) = 0.156377 / 8.4704×10⁻⁴ = 184.6×**

#### Interpretation and necessary distinctions

Under this topology, parameterisation and ambient scenario, the model permits a **184.6-fold**
asymptotic increase in steady-state open probability. With the desensitisation branch included,
steady-state open probability at *saturating* GABA also settles to this same three-state
asymptote (≈ 0.1564), while the early transient peak reaches 0.750 in simulation
(analytic bound 0.8282).

Hierarchical observables must not be conflated:

* **P_open** — a microscopic channel-state variable.
* **Whole-cell current** — I_tonic = N · γ · P_open · (V_m − E_Cl), so receptor density,
  single-channel conductance and driving force all intervene.
* **Tonic inhibition** — depends on chloride homeostasis (KCC2/NKCC1); prolonged
  high-conductance chloride flux can shift E_Cl depolarising.
* **Circuit shunting** — depends non-linearly on input resistance relative to leak and
  synaptic conductances.

A 184.6× open-probability dynamic range is therefore **not** a 184.6× increase in circuit-level
inhibition. It is the upper bound on receptor open probability within this scheme.

### 2.4. Mapping operational potency shifts onto Q

To connect phenomenological compound descriptors to microscopic rates:

1. **s_max** — the maximum fold leftward shift of the agonist concentration-response curve at
   full modulator occupancy: s_max = EC₅₀(baseline) / EC₅₀(c_affinity).
2. **Implicit inversion** — the unique c_affinity ≥ 1 reproducing a given s_max is found by
   root-solving F(c) = EC₅₀(Q(k_off)) / EC₅₀(Q(k_off/c)) − s_max = 0 via Brent's method, with
   EC₅₀ evaluated by integrating the master equation over a 300 ms step.

**These two numbers are not interchangeable**, and the previous draft's Table 3 listed
c = 2.79 as s_max = 2.50. In this parameterisation:

| s_max | c_affinity |
|---|---|
| 2.40 | 2.7882 |
| **2.50** | **2.9321** |

and conversely c = 2.79 gives s_max = **2.401**. At c = 2.9321 the model's EC₅₀ moves
20.0000 → **8.0000 µM**, i.e. exactly 2.5×.

3. **Occupancy weighting.** At fractional modulator occupancy θ = [D]/([D] + K_d,PAM), the
   receptor pool is a linear mixture of unmodulated and fully modulated channels, so standing
   tonic gain scales linearly:
   Tonic Gain(θ) = 1 + θ·(Gain_max(s_max) − 1).

**Table 2. Reported endpoints and their mapping to s_max.**

| Compound | Reported primary endpoint | Assay / construct | Model `ceiling` | Citation status |
|---|---|---|---|---|
| Imepitoin | ~20% max potentiation rel. diazepam (Rundfeldt & Löscher, 2014) | recombinant α1β2γ2, oocyte | **1.25** | resolved |
| TPA023 | α2/α3 partial; silent antagonist at α1 and α5 | recombinant | **2.00** | not re-audited |
| L-838,417 | partial efficacy ~30–40% at α2/3/5; silent at α1 (McKernan et al., 2000) | membrane potential | **2.20** | resolved |
| MP-III-022 | α5-selective potentiation, weaker partial modulation at α2/α3 (Stamenić et al., 2016) | patch clamp, α5β3γ2 | **2.50** | resolved |
| Alogabat (RG7816) | +167% (rat), +72% (human) potentiation at EC₂₀ GABA; no potentiation at α1/α2/α3 (Cecere et al., 2025) | recombinant α5β3γ2 patch | **2.50** | resolved |
| Diazepam / midazolam | two distinct and separable mechanisms of potentiation (Walters et al., 2000) | recombinant α1β2γ2 patch | **2.50** | **corrected** — see note |
| Neurosteroid (allopregnanolone proxy) | gating-active, non-selective | — | **6.00** | not re-audited |
| Gaboxadol (THIP) | δ orthosteric agonist | — | **1 × 10⁹** | not re-audited |

**Notes.** (i) The `ceiling` column is the repository's own `PROFILES[...].ceiling`, i.e. the
value the model actually uses. The previous draft's Table 6A listed alogabat at 2.40 and
MP-III-022 and HZ-166 at 2.20; all three are **2.50** in the code. The selectivity index R is
unaffected (s_max does not enter it, §4.4 note 2) but the table misreported its own inputs.
(ii) Reported experimental metrics (relative efficacy %, partial agonism, current potentiation
at EC₂₀) are **not** direct EC₅₀ fold shifts. Each mapping to a `ceiling` is an explicit
modelling assumption, and none of these values is a measured EC₅₀ shift for that compound.
(iii) The diazepam citation is corrected from the previous draft, which attributed it to
*Br. J. Pharmacol.* 131:1307–1314 under a different title; the paper is Walters et al. (2000),
*Nat. Neurosci.* 3:1274–1281. This matters because diazepam is the reference arm of every
selectivity ratio reported in §4.4. (iv) The MP-III-022 citation is likewise corrected: the
previous draft cited "Fischer et al. (2010), *Neuropharmacology* 59:612–618", which does not
resolve and is anachronistic — MP-III-022 was first characterised in Stamenić et al. (2016).

### 2.5. Sensitivity and robustness of the dynamic range

Because the 184.6× ratio depends on the parameter set, we swept each quantity one at a time
(`scripts/paper_numbers.py --section sensitivity`).

**Table 3. One-at-a-time sensitivity of the modelled open-probability dynamic range.**

| Perturbation | P_open(0.40 µM) | P_open,∞ | dynamic range |
|---|---|---|---|
| **nominal** | 0.0008470 | 0.15638 | **184.6×** |
| desensitisation entry d × 0.5 | 0.0008489 | 0.26308 | 309.9× |
| desensitisation entry d × 2.0 | 0.0008433 | 0.08634 | 102.4× |
| resensitisation r × 0.5 | 0.0008433 | 0.08634 | 102.4× |
| resensitisation r × 2.0 | 0.0008489 | 0.26308 | 309.9× |
| gating β × 0.5 | 0.0004237 | 0.08482 | 200.2× |
| gating β × 2.0 | 0.0016926 | 0.27046 | 159.8× |
| affinity K_d × 0.5 | 0.0032505 | 0.15638 | 48.1× |
| affinity K_d × 2.0 | 0.0002154 | 0.15638 | 725.8× |
| ambient 0.10 µM | 0.0000543 | 0.15638 | 2881.1× |
| ambient 0.20 µM | 0.0002154 | 0.15638 | 725.8× |
| **ambient 0.40 µM (nominal)** | 0.0008470 | 0.15638 | **184.6×** |
| ambient 0.70 µM | — | 0.15638 | 62.1× |
| ambient 0.80 µM | 0.0032505 | 0.15638 | 48.1× |
| ambient 1.50 µM | 0.0104257 | 0.15638 | 15.0× |
| ambient 3.00 µM | 0.0324370 | 0.15638 | 4.8× |

**Finite modulator shifts** (the realisable, as opposed to asymptotic, gains):

| c_affinity | s_max | tonic gain at 0.40 µM |
|---|---|---|
| 1.50 | 1.433 | 2.21× |
| 2.00 | 1.834 | 3.84× |
| 2.79 | 2.401 | 7.18× |
| 2.9321 | 2.500 | 7.88× |
| 5.00 | 3.708 | 20.25× |
| 10.00 | 5.725 | 56.46× |
| 50.00 | 10.748 | 158.08× |

#### Key observations

1. **Desensitisation ratio (d/r).** The asymptote is governed directly by d/r. Doubling
   desensitisation entry (d/r = 50) still leaves **102.4×**; halving it (d/r = 12.5) expands the
   range to **309.9×**. This is the single largest source of uncertainty, and it is the one
   parameter the three transient anchors do not constrain.
2. **Gating ratio (β/α).** Perturbing β by ±2-fold moves baseline and asymptote nearly in
   proportion, leaving the range stable at **159.8–200.2×**.
3. **Ambient GABA.** Across the 0.2–0.8 µM range reported for cortex and hippocampus the range
   spans **725.8× to 48.1×**. Note the lower end falls marginally *below* 50, so the robust
   claim is **"> 48×"**, or "> 50× across 0.2–0.7 µM" (62.1× at 0.70 µM). The previous draft's
   "> 50×" sat one unit on the wrong side of the computed value.
4. **K_d and ambient are the same axis.** Halving K_d at fixed [G] gives exactly the numbers of
   doubling [G] (0.0032505, 48.1×), as it must, since the equilibrium depends on [G]/K_d. Both
   rows are retained because the two have different experimental meanings.
5. **When a 2.5× cap would actually bind.** With x = [G]/K_d and A = 1 + E + d/r = **30.8195**,
   the stationary ratio is

   P_open,∞ / P_open([G]) = 1 + 2/(Ax) + 1/(Ax²)

   Setting this to 2.5 and solving the quadratic 1.5Ax² − 2x − 1 = 0 gives x = 0.17045, i.e.

   **[G] = 5.06 µM** (verified: the range at 5.06 µM is 2.500×).

   A fixed 2.5× scalar cap therefore binds only when ambient GABA exceeds ≈ 5 µM — an order of
   magnitude above the reported submicromolar range. In the regime the extrasynaptic pool
   actually occupies, the cap is not a saturation boundary; it is a truncation.

### 2.6. Comparison with gating modulation

Under constant non-zero agonist, steady-state open probability in the affinity-only limit
approaches P_open,∞ = **0.1564**. This stationary limit does not bound transient peaks, which
depend on the agonist waveform and reach **0.4196** (baseline) to **0.6998** (limiting) in the
synaptic cleft, strictly below the analytic gating bound β/(α+β) = **0.8282**.

A gating PAM (β → β·c) instead scales the forward opening transition:

lim_{c→∞} P_open = 1 / [1 + (α/βc)(1 + d/r)] → **1.0**

So increasing the forward gating rate can drive both steady-state and transient P_open toward
unity, whereas steady-state affinity-only modulation remains bounded by the desensitised-state
equilibrium at 0.1564. This distinction tracks the divergent pharmacology of
neurosteroids and barbiturates (gating-active) versus benzodiazepines (predominantly
affinity-modulating), and it is why the neurosteroid arm in §4.4 carries `ceiling = 6.0`.

---

## 3. Computational architecture and self-correction

Two modelling traps distorted earlier PAM safety analyses in this project.

```
[Legacy platform (simulator.py)]
  |-- conductance scaled linearly: g_eff = g_base * (1 + occ * (gain - 1))
  `-- clamped at a hard-coded boundary: min(g_eff, 2.5)   <-- TRAP 1: fixed scalar cap
         |
[Intermediate recalibration (overdose_kinetic.py, early)]
  |-- split into synaptic vs extrasynaptic pools
  `-- SENS_TOTAL = 0.10 carried over                     <-- TRAP 2: calibration collapse (E12)
      `-- drug effect collapsed to -2% ventilation; artifact: 100% survival at every dose
         |
[Audited engine (circuitpharm, tag manuscript-v2)]
  |-- 5-state Markov generator (Q matrix)
  |-- dual-pool compartmentalisation: synaptic 1.32x vs extrasynaptic 184.6x range
  |-- subunit distribution estimates (preBotC vs forebrain), with per-cell provenance
  |-- two neuron substrates (LIF and Butera-Rinzel-Smith conductance)
  `-- epistemic guardrails: Tier.VOID for unanchored margins, Tier.VALIDATED for R
```

### 3.1. Trap 1 — the fixed scalar cap

Early implementations scaled PAM gain linearly with occupancy and clamped conductance:

```python
# legacy formulation (REJECTED)
gaba_a_gain = min(1.0 + occupancy * (target_gain - 1.0), gaba_a_efficacy_cap)  # cap = 2.5
```

This manufactured a safety guarantee: because conductance was clamped at 2.5×, circuits
inevitably survived dose escalation. The "overdose ceiling" was a programmer's boundary, not an
emergent property of receptor saturation. In the revised framework the boundary emerges from
Q — and §2.5 shows where it actually lies (≈ 5 µM ambient GABA, not at any fixed multiplier).

### 3.2. Trap 2 — the calibration collapse (E12)

Moving from a single-pool to a dual-pool conductance model, `scripts/overdose_kinetic.py`
carried over the lumped `SENS_TOTAL = 0.10` from the single-pool fit. The mismatch collapsed
the drug's respiratory effect to **−2% ventilation** for a non-selective benzodiazepine, against
a clinical anchor of −16% to −19%. The resulting artifact: the script reported that *every
compound survived every dose up to 100% occupancy* — a false impression of universal overdose
safety produced by a parameter that no longer meant what it had meant.

The insight that followed: **tissue-level sensitivity (`gaba_sens`) is not a subunit fraction.**
It is a lumped parameter absorbing extra-preBötC mechanisms — chemoreflex blunting, upper-airway
motor tone — that an isolated pacemaker model does not contain.

E12 has since recurred in two further forms, both documented: a tonic/phasic split applied
everywhere except the module written last, and (most recently) a shared population factory
whose default silently removed spike-triggered adaptation from the spinal circuit while leaving
the respiratory one correct. The general shape is a value duplicated across two locations that
then diverge, with a *plausible* rather than loud symptom. The structural remedy adopted here is
to delete the second copy rather than synchronise it.

### 3.3. Epistemic reliability tiers

`circuitpharm.results.Tier` marks every reported quantity:

* **`Tier.VOID`** — quantities lacking whole-animal anchors (absolute lethal-dose margins in
  mg) are blocked from analysis programmatically.
* **`Tier.VALIDATED` (scale-invariant)** — quantities independent of the shared multiplicative
  tissue-scaling constants, such as the selectivity index

  **R = [Drive_subj(compound) / Burden_resp(compound)] / [Drive_subj(diazepam) / Burden_resp(diazepam)]**

Because the tissue-scaling constants multiply numerator and denominator alike, they cancel, so
the compound ordering is invariant to that scale. **Cancellation establishes mathematical
scale-invariance within the model's algebra; it does not validate the biological accuracy of
the underlying subunit fractions or circuit weights, and it does not compute clinical safety,
human ventilatory depression, or an overdose margin.**

A third guardrail operates at the level of reported statistics: `ranking_robustness.py` prints
an explicit **NOT QUOTABLE** verdict for the Monte Carlo median and 95th percentile (§4.4.3).
The previous draft placed both in a headline table. This manuscript does not.

---

## 4. Results

### 4.1. Non-equilibrium gating dissociates synaptic and extrasynaptic responses

Simulating an affinity-type PAM calibrated to a 2.5× leftward EC₅₀ shift
(**c_affinity = 2.9321**) reveals a pronounced divergence between compartments.

**Table 4. Divergence of kinetic predictions across receptor compartments.**

| Compartment / observable | Agonist regime | Baseline | With PAM | Ratio |
|---|---|---|---|---|
| Simulated synaptic peak P_o | transient, 1.0 mM, τ_clear 0.30 ms | 0.419609 | 0.553593 | **1.319×** (peak gain) |
| Synaptic τ_IPSC | mono-exp fit, 90→10% decay | 15.000 ms | 28.426 ms | **1.895×** |
| Time-integrated P_o (charge proxy) | 200 ms window | 7.454 ms | 18.030 ms | **2.419×** |
| Synaptic max dynamic range | limiting transient, k_off → 0⁺ | 0.419609 | 0.699800 | **1.668×** |
| Tonic standing P_o | steady state, 0.40 µM ambient | 0.000847 | 0.006671 | **7.876×** (tonic gain) |
| Tonic max dynamic range | limiting steady state, k_off → 0⁺ | 0.000847 | 0.156377 | **184.608×** |

*Charge proxy: 7.454 ms ≡ 0.007454 pC/pA; 18.030 ms ≡ 0.018030 pC/pA. The simulated synaptic
peaks (0.4196 baseline, 0.5536 with PAM, 0.6998 limiting) all respect the analytic gating bound
β/(α+β) = 0.8282, which is not itself attainable because desensitisation competes with
activation throughout the rise.*

1. **Synaptic cleft.** The brief 1 mM transient drives peak open probability to 0.420. A
   leftward affinity shift raises that peak only **1.319×**, with an asymptotic ceiling of
   **1.668×** — so roughly two-thirds of the available synaptic headroom is already consumed by
   a 2.5× shift. Instead, slower unbinding prolongs the deactivation tail: τ_IPSC rises
   **1.895×** (15.00 → 28.43 ms) and integrated open probability **2.419×**, reproducing the
   classical electrophysiological phenotype of benzodiazepine action on IPSC decay
   (Otis & Mody, 1992).
2. **Extrasynaptic space.** Receptors at steady-state 0.40 µM operate far below saturation
   (P_open = 8.5 × 10⁻⁴). Slower unbinding moves the activation threshold into the ambient
   range, producing a **7.876×** potentiation of standing open probability, with theoretical
   headroom extending to **184.6×**.

**The contrast between the two compartments' ceilings is 184.608 / 1.668 = 111-fold.** This is
the manuscript's central quantitative claim, and it is weaker than the previous draft's
asserted 187-fold — because the synaptic ceiling in this parameterisation is 1.67×, not 1.12×.
The qualitative conclusion is unchanged: two orders of magnitude separate the compartments, and
a single scalar multiplier cannot describe both.

A caution the previous draft did not state. **Headroom and reachable gain are different
quantities, and they diverge in opposite directions as ambient GABA falls.** Between 0.40 and
0.10 µM the asymptotic headroom rises 15.6-fold (184.6× → 2881×) while the gain a finite
s_max ≈ 2.5 modulator actually extracts rises only 1.08-fold (7.88× → 8.47×). Any statement
about what a real compound does must be written against the reachable column. §5.2 returns to
this, because the previous draft's falsification criteria were set from the wrong one.

### 4.2. Subunit distribution and anatomical decoupling: the α5 case

**Table 5. Subunit compartmentalisation and regional distribution estimates, with per-cell
provenance** (basis from `circuitpharm.provenance.report()`; `python -c "from
circuitpharm.provenance import report; print(report())"`).

| Subunit | Extrasyn. share `f_extra` | basis | preBötC `f` | basis | Forebrain `f` | basis |
|---|---|---|---|---|---|---|
| α1 | 0.15 | **UNSOURCED** | 0.60 | FROM_QUALITATIVE [`pbc_alpha`] | 0.35 | **UNSOURCED** |
| α2/3 | 0.20 | **UNSOURCED** | 0.15 | FROM_QUALITATIVE [`pbc_alpha`] | 0.25 | **UNSOURCED** |
| α5 | 0.80 | FROM_QUALITATIVE | 0.02 | **UNSOURCED** | 0.30 | FROM_QUALITATIVE |
| δ/α4 | 1.00 | FROM_QUALITATIVE [`pbc_delta`] | 0.15 | FROM_QUALITATIVE [`pbc_delta`] | 0.08 | **UNSOURCED** |
| ε | 0.50 | **GUESS** | 0.08 | FROM_QUALITATIVE [`pbc_eps`] | 0.02 | **UNSOURCED** |

#### Provenance, stated at the resolution the audit supports

The repository's audit reports that **6 of 28 load-bearing parameters name a source (21%)**.
That is the honest figure for this table, and it is substantially weaker than the previous
draft's attribution of Table 5 to "regional expression literature (Pirker et al., 2000;
Kasugai et al., 2010)". Specifically:

* **`f_α5` in preBötC = 0.02 is UNSOURCED.** It is a modelling assumption encoding an
  anatomical ordering (medullary α5 low, forebrain α5 high), not a measurement from isolated
  preBötC tissue. §4.4.3 shows the consequences of that assumption in full.
* **Seven of the fifteen cells above are UNSOURCED, an eighth (`f_extra,ε`) is labelled an
  outright GUESS in the source code, and the remaining seven are FROM_QUALITATIVE** — i.e.
  a direction taken from a source, with the number itself ours. No cell in this table is
  QUANTITATIVE.
* The preBötC α1/α2/α3 column rests on `pbc_alpha` — Liu & Wong-Riley (2004),
  *J. Appl. Physiol.*, metadata verified by DOI but **full text inaccessible (HTTP 403)**, so
  claim support is formally **UNASSESSED**. Two cautions travel with it: it covers α1, α2 and α3
  only (not α5, δ, α4 or ε, to which `REGIONS["prebotc"]` also assigns values), and it is a
  *developmental* study — an axis this project has been bitten by twice.
* The repository's own internal source for regional distribution, `a5_dist`, **resolves
  perfectly by DOI and does not support the numbers attributed to it**: it is a 1988 study using
  a single generic cDNA probe for "the α subunit", reporting total α-subunit mRNA by region
  (medulla ≪ hippocampus/cortex). That is a claim about regional *level*; `REGIONS` encodes
  regional *composition* (rows sum to 1 by construction). The source is not merely coarse — it
  concerns a different quantity.
* Kasugai et al. (2010) **is** the better source for the α5 extrasynaptic fraction, and is cited
  here at its correct coordinates (*Eur. J. Neurosci.* 32:1868–1888; the previous draft gave
  *J. Neurosci.* 30:14024–14035). It quantitatively establishes that α5 is predominantly
  extrasynaptic in hippocampus, supporting `f_extra,α5 = 0.80` as a direction. Pirker et al.
  (2000) maps distribution qualitatively, showing high medullary α1 and low medullary α5.

Two consequences follow, both weaker than the previous draft's:

1. **Mechanistic basis of respiratory sparing.** Low α5 representation in the preBötC provides
   a mechanistic basis for reduced *direct* α5-mediated inhibition of the core rhythm generator.
   Whole-animal respiratory preservation cannot be inferred from preBötC receptor abundance
   alone, since ventilatory stability also depends on chemoreflex integration and upper-airway
   motor tone — the very mechanisms `gaba_sens` lumps (§3.2).
2. **Forebrain tonic potentiation.** In hippocampus and cortex, α5 expression is substantial
   (≈ 30%) and predominantly extrasynaptic (> 80%), so an α5 PAM acts directly on the
   high-headroom pool.

Consequently, while a 184.6× receptor dynamic range does not produce 184.6× respiratory
depression, it does alter expected forebrain behaviour: dose escalation of a high-efficacy α5
PAM can drive substantial tonic shunting conductances, with cognitive blunting, memory
disruption and sedation as the expected costs.

### 4.3. Dependence on the maximum operational potency shift

Because extrasynaptic headroom is large, potentiation is limited by the ligand's own s_max —
the maximum left-shift it can elicit at full occupancy:

Tonic Gain(occupancy) = 1 + occupancy · (Gain_max(s_max) − 1)

* **Low potency-shift modulators.** Imepitoin (`ceiling` = 1.25, c = 1.2804) meets a narrow
  biophysical ceiling under this mechanism: at 100% occupancy and 0.40 µM ambient GABA its
  tonic open-probability gain cannot exceed **1.62×**, against a phasic peak gain of 1.08×.
  Tonic gain rises steeply from there — TPA023 at `ceiling` = 2.00 reaches **4.70×** and
  L-838,417 at `ceiling` = 2.20 reaches **5.87×**. The model therefore predicts substantially
  greater tonic headroom once s_max exceeds ≈ 1.3–1.5; that range should be **experimentally
  evaluated as a candidate operational boundary, not treated as an established safety
  threshold**. (The previous draft assigned L-838,417 s_max = 1.50 in its Table 2 where the
  code says 2.20, moving this compound out of the low-shift group entirely.)
* **High-efficacy scenarios.** For alogabat, recombinant electrophysiology shows selective
  potentiation of EC₂₀ GABA currents (+167% rat, +72% human α5β3γ2; Cecere et al., 2025).
  Rather than asserting a universal measured EC₅₀ fold shift, s_max = 2.50 (c = 2.9321,
  **7.88×** tonic gain) is examined as a modelled high-efficacy benchmark calibrated to that
  profile. MP-III-022 carries the same `ceiling = 2.50` in the code, based on partial α5
  potentiation data (Stamenić et al., 2016).
* **High potency-shift modulators (s_max ≥ 2.5).** Diazepam at `ceiling` = 2.50 (c = 2.9321)
  gives **7.88×** tonic gain. The relationship is strongly supralinear beyond that: c = 5.00
  (s_max 3.71) gives **20.25×**, c = 10.0 (s_max 5.73) gives **56.46×**, and the gating-active
  neurosteroid profile at `ceiling` = 6.00 (c = 10.9059) gives **62.69×** tonic gain against a
  phasic peak gain of only 1.55×. The neurosteroid arm is the clearest case of the compartment
  asymmetry: a 40-fold difference between what the compound does to a standing extrasynaptic
  conductance and what it does to a synaptic peak.

### 4.4. Algebraic selectivity index and Monte Carlo robustness

> **Scope.** This section is **algebraic, not multiscale.** `scripts/ranking_robustness.py`
> imports `circuitpharm.subtypes` and nothing else — no neuron, no circuit, no simulation. R is
> a ratio of weighted subunit sums. The previous draft titled this section "Multiscale Circuit
> Selectivity Pipeline", which overstates it by a whole layer of the model. Circuit-simulation
> results appear separately, in §4.6.

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
questionable**, since the one source bearing on it is titled for an α1 contribution. The
previous draft attributed these weights to "Saba et al. (2017), *Alcohol. Clin. Exp. Res.*
41:748–758", which does not resolve; that attribution is withdrawn rather than replaced.

The computation proceeds in five deterministic steps, with
w(s) = ρ·f_extra,s + (1 − f_extra,s):

1. **Efficacy vector e** per compound (Table 6A).
2. **Kinetic gain weighting ρ.** The operational weighting of extrasynaptic (tonic) against
   synaptic (phasic) modulation. For s_max = 2.50 the tonic gain is 7.876×, so
   ρ_peak = 7.876 / 1.319 = **5.97** and ρ_charge = 7.876 / 2.419 = **3.26**. The nominal
   **ρ = 6.00** therefore now sits essentially *at* the peak-based ratio rather than between the
   two; under the previous draft's parameterisation it lay between 6.79 and 4.36. Monte Carlo
   sweeps sample ρ ~ U[3.0, 7.5], which still spans both endpoints.
3. **Forebrain subjective drive.** Drive_subj = Σ_s f_forebrain,s · e_s · w_subj,s · w(s)
4. **PreBötC respiratory burden.** Burden_resp = Σ_s f_prebotc,s · e_s · w(s)
5. **Selectivity ratio R**, normalised to diazepam (§3.3).

**Table 6. Selectivity ranking.** Nominal R, and the 5th-percentile R over 20,000 Monte Carlo
draws. *The Monte Carlo median and 95th percentile are deliberately omitted: the model marks
them NOT QUOTABLE (§4.4.3), and they are reported there as a prior-sensitivity diagnostic
rather than as results.*

| Compound | Profile | `ceiling` | Nominal R | 5th-pct R | P(R > 1) |
|---|---|---|---|---|---|
| Ideal α5 PAM | α5-exclusive | 2.50 | **10.8750×** | 2.57× | 99.9% |
| Alogabat (RG7816) | α5-selective (Phase II) | 2.50 | **8.6442×** | 2.51× | 99.9% |
| SH-053-2′F-R-CH₃ | α5-selective enantiomer | 2.50 | **9.3487×** | — | — |
| MP-III-022 | α5-selective | 2.50 | **7.8750×** | 2.49× | 99.9% |
| SH-053-2′F-S-CH₃ | α2/α3/α5 enantiomer | 2.50 | **3.6250×** | — | — |
| HZ-166 / KRM-II-81 | α2/α3-preferring | 2.50 | **1.2083×** | 0.34× | **63.8%** |
| Non-selective BZ | α1/2/3/5 (diazepam) | 2.50 | **1.0000×** (ref) | — | — |
| Neurosteroid | non-selective, gating PAM | 6.00 | **0.5633×** | 0.33× | **0.0%** |
| Gaboxadol | δ orthosteric agonist | 1 × 10⁹ | **0.0000×** | 0.00× | 0.0% |

**Notes.**
1. **R is model-conditional.** It measures the ratio of forebrain subjective drive to preBötC
   pacemaker burden, relative to diazepam. The ranking is conditional on the stated parameters:
   ρ = 6.00 is a fitted weighting; `w_subj` assigns zero contribution to α1, δ and ε *by
   design*, and eight of the fifteen Table 5 cells carry no source at all. **R does not compute clinical
   safety, human respiratory depression, or an overdose margin.**
2. **s_max does not enter R.** The `ceiling` column is each compound's maximum operational
   potency shift, used in dose escalation and gating bounds. R evaluates non-saturating
   proportional scaling governed by **e** and ρ alone. The previous draft's note said the same
   while its Table 6A listed s_max values inconsistent with the code.
3. **Gaboxadol's zero is a scope artefact, not a pharmacological verdict.** Gaboxadol is a δ
   orthosteric agonist. Because `w_subj` is parameterised for ethanol-like discriminative
   stimulus salience mediated by α2/3 and α5, `w_subj,δ = 0.0` forces Drive_subj = 0 and hence
   R = 0. This reflects the domain-specific design of the drive metric and implies **nothing**
   about gaboxadol's hypnotic, sedative or thalamocortical efficacy.
4. **HZ-166 is the arm that does not clearly separate**, at 63.8% of draws above the diazepam
   reference and a 5th percentile of 0.34×.

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
   them makes the reported figures conservative. The previous draft reported "20,000 draws"
   without this disclosure.

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

### 4.5. Why systems models require empirical efficacy: stereochemical invariance

Two-dimensional graph representations and topological fingerprints are widely used to predict
candidate properties. We tested whether they can distinguish stereoisomers with divergent
subtype selectivities, using SH-053-2′F-R-CH₃ (α5-selective, nominal R = **9.3487×**) and
SH-053-2′F-S-CH₃ (α2/α3/α5, nominal R = **3.6250×**).

**Table 7. Stereochemical invariance across enantiomers.**

| Descriptor | (R)-enantiomer | (S)-enantiomer | Δ |
|---|---|---|---|
| Molecular weight | 388.42 g/mol | 388.42 g/mol | 0.00 |
| Topological polar surface area | 61.86 Å² | 61.86 Å² | 0.00 |
| MolLogP | 2.841 | 2.841 | 0.00 |
| Morgan fingerprint (2048 bit) | bit-identical | bit-identical | 0 bits |
| Internal distance matrix | exact invariant | exact invariant | < 1e-9 Å |
| Vacuum MMFF energy | 42.184 kcal/mol | 42.184 kcal/mol | < 1e-6 kcal/mol |
| **Selectivity index R** | **9.3487×** | **3.6250×** | **+5.7237×** |

Enantiomers are related by an improper rotation, so all internal pairwise atomic distances and
2D graph invariants are mathematically identical. Three-dimensional representations capture
geometry, but **these 2D descriptors cannot by themselves encode the stereochemical information
responsible for the pharmacological difference without explicitly modelling the chiral binding
pocket.** Multiscale systems pharmacology must therefore ingest empirical, subtype-specific
efficacy measurements rather than infer them from 2D structural proxies. Note that the R values
on the right-hand column are *consequences* of the measured efficacy vectors, not predictions
from structure — which is precisely the point.

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

**3. What is not established.** Whether that is *more* bursts or *longer* bursts. The burst
detector reported 5.63 Hz on the conductance substrate against the FFT's 0.302 Hz — a 19-fold
discrepancy — because a 0.35 × max threshold applied after a 20 ms rate low-pass fires
repeatedly *within* a single burst. A burst-level detector matched to this substrate is
outstanding work.

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
Representing PAM action as a fixed ≈ 2.5× scalar cap conflates synaptic receptor saturation
with extrasynaptic headroom. In synapses, agonist exposure restricts potentiation (**1.32×**
peak gain, **1.67×** asymptotic ceiling); in extrasynaptic compartments, submicromolar ambient
GABA leaves a broad range (**> 48×** across 0.2–0.8 µM, **184.6×** at 0.40 µM). A 2.5× cap
binds only above ≈ **5 µM** ambient GABA (§2.5), an order of magnitude above the extrasynaptic
regime — so in that regime the cap is a truncation, not a saturation boundary.

Two implications for α5-targeted development:

1. **Screen for bounded s_max.** Overdose safety cannot be assumed from the "PAM mechanism"
   alone. For predominantly extrasynaptic targets the model predicts substantially greater
   tonic headroom above s_max ≈ 1.3–1.5; this warrants **experimental evaluation as a candidate
   efficacy boundary, not interpretation as an established safety threshold.** The supralinearity
   matters here: s_max 2.50 → 7.88× tonic gain, but s_max 3.71 → 20.25× and s_max 5.73 → 56.46×.
2. **Separate synaptic and tonic conductances.** Computational neural models must decouple
   synaptic deactivation (τ_IPSC, **1.90×** here) from standing extrasynaptic conductance
   (**7.88×** here) rather than applying one lumped multiplier to all inhibitory inputs. Doing
   otherwise is error E12 (§3.2), which has now recurred three times in this project.

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
  tag `manuscript-v2`, for a high-efficacy affinity PAM (s_max 2.40–2.50):

  | [GABA]_bath | R_PAM predicted (s_max 2.40 – 2.50) | asymptotic ceiling (k_off → 0⁺) |
  |---|---|---|
  | 0.1 µM | **7.66 – 8.47×** | 2881.1× |
  | 0.4 µM | **7.18 – 7.88×** | 184.6× |
  | 1.0 µM | **5.87 – 6.33×** | 31.5× |
  | 3.0 µM | **2.94 – 3.03×** | 4.8× |
  | 10.0 µM | **1.34 – 1.35×** | 1.5× |

  *(c_affinity = 2.7882 and 2.9321 respectively.)*

**The asymptote column is unreachable by any finite modulator and must not be used to set
criteria.** The previous draft pre-registered "R_max > 15× at 0.1 µM", which is an asymptote
reading; the realisable prediction there is ≈ 8×. A laboratory measuring 8× would have reported
the model falsified when the model in fact predicts 8×. Two of that draft's four intervals were
violated by the model they claimed to describe. The structural reason is in §4.1: as ambient
GABA falls, asymptotic headroom grows without limit while reachable gain barely moves. Both
columns are printed side by side by `scripts/paper_numbers.py --section falsification`, and a
regression test now asserts the divergence so that any future criterion written against the
wrong column fails loudly.

**Falsification rules** (pre-registered against the reachable column):

1. If a high-efficacy α5 PAM (s_max ≥ 2.4) shows maximal steady-state **R_PAM ≤ 2.5 at
   [GABA]_bath = 0.40 µM** in recombinant α5β3γ2 channels, the nominal Model-B parameterisation
   is falsified under the tested construct and conditions. (Model B predicts 7.18–7.88× there,
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

**Receptor layer**

* **Kinetic topology.** Five states, with affinity modulation acting through k_off. Ligands
  altering gating or desensitisation kinetics will show different quantitative profiles (§2.6).
* **Desensitisation is unconstrained by the anchors.** d/r = 25.0 is provisional, and it sets
  the asymptote directly: ±2-fold moves the dynamic range over 102.4–309.9× (§2.5). Every
  asymptotic figure in this manuscript is conditional on it.
* **Three anchors do not identify six rates.** This is one admissible parameterisation.
* **Ambient GABA.** Evaluated at a representative 0.40 µM. In vivo, GAT-1/GAT-3 transport and
  synaptic spillover create dynamic microdomains of varying concentration.

**Distribution layer**

* **Seven of fifteen Table 5 cells are UNSOURCED, an eighth is a GUESS, the other seven are
  FROM_QUALITATIVE, and preBötC α5 — the single most consequential value — has no source**
  (§4.2). None is QUANTITATIVE. Replacing these point values with sourced
  *ranges* and re-running the robustness analysis is the highest-value outstanding work.
* `KAPPA = 15.0`, the Dirichlet concentration governing the entire uncertainty model, is
  **UNSOURCED**.
* The repository's internal regional-distribution source does not support the quantity it is
  attributed to (§4.2), and `pbc_alpha`'s claim support is formally UNASSESSED (full text 403).

**Circuit layer**

* **The A4 riluzole dissociation is NOT reproduced** (0 of 3 seeds). This is an *inherited*
  limitation rather than a tuning failure: Butera–Rinzel–Smith model 1 **is** the pacemaker
  hypothesis, and riluzole is the principal published argument against it. It is declared, not
  retuned.
* **The coupling weights are ours, not published.** Butera, Rinzel & Smith (1999b) — the
  population paper — was unreachable (journal HTTP 403; every available machine-readable
  encoding is single-cell), so the network weights were anchored by us.
* **Respiratory calibration remains wet-lab blocked**, and the conductance migration moved the
  > P12 muscimol anchor *further* away by importing neonatal parameters.
* **The conductance substrate is anchored to an in vitro band** (0.05–1.00 Hz neonatal slice),
  so its respiratory outputs are slice claims, not ventilatory ones.
* **The non-selective benzodiazepine reference arm changes sign between substrates** (§4.6), and
  it is the denominator of every R in Table 6.
* Four measurement protocols inherited from the LIF substrate had to be re-derived for the
  conductance cell (warm-up against τ_h = 10 s, validity band, FFT resolution, burst detection
  threshold). **None failed loudly.** The working rule adopted: on a new substrate, assume every
  inherited protocol is wrong until re-derived.

**Scope**

* **Pharmacokinetics are absent.** The model evaluates receptor occupancy and network dynamics,
  not absorption, blood-brain barrier penetration or clearance, which determine clinical
  dose-occupancy relationships.
* **No absolute margin is computed.** `Tier.VOID` blocks lethal-dose margins programmatically.
  R orders candidates; it cannot bound risk. Bounding risk needs the in vitro anchor and
  per-molecule measurement of intrinsic allosteric efficacy.
* **The work stays at mechanism and compound-class level**, and reports no formulation or dosing.

---

## 6. Conclusions

A five-state continuous-time Markov gating scheme shows that allosteric potentiation headroom is
state- and compartment-dependent. The classical ceiling binds tightly in the synapse
(**1.32×** peak gain, **1.67×** asymptotic ceiling, with the drug's action redirected into the
deactivation tail at **1.90×** τ and **2.42×** charge), but expands to a **184.6-fold**
asymptotic open-probability dynamic range in extrasynaptic microenvironments at submicromolar
ambient GABA — an **111-fold** contrast between the two compartments' ceilings. A fixed 2.5×
scalar cap binds only above ≈ 5 µM ambient GABA, so in the extrasynaptic regime it truncates
rather than saturates.

Under the model's assumption of low α5 representation in the preBötzinger complex — a value with
no source — α5-selective PAMs reduce direct inhibition of the modelled core rhythm generator,
while their high extrasynaptic localisation exposes forebrain circuits to substantial tonic
conductance increases during dose escalation. The resulting compound ordering is robust in two
independent senses: it survives adverse draws of every invented parameter simultaneously
(5th-percentile R ≥ 2.48× for all three α5-selective arms, in every prior-floor scenario), and
it survives replacing the neuron model entirely (Spearman ρ = +1.0000 across an
integrate-and-fire and a conductance-based substrate). The non-selective benzodiazepine
reference arm does **not** survive the latter test — its sign on respiratory output reverses —
and that arm is the denominator of every ratio reported here.

Overdose-relevant conductance headroom therefore depends jointly on ligand efficacy, receptor
mechanism, ambient agonist concentration, compartmentalisation and network context. The
quantity that a fixed scalar cap was standing in for is not a constant, and it is not small.

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
   *(Corrected from the previous draft, which gave 2001 and pages 8993–8997.)*
7. Farrant, M., & Nusser, Z. (2005). Variations on an inhibitory theme: phasic and tonic
   activation of GABA-A receptors. *Nat. Rev. Neurosci.*, 6(3), 215–229.
8. Haas, K. F., & Macdonald, R. L. (1999). GABA-A receptor subunit γ2 and δ subtypes confer
   unique kinetic properties on recombinant GABA-A receptor currents in mouse fibroblasts.
   *J. Physiol.*, 514(1), 27–45. *(Corrected from the previous draft, which gave*
   J. Neurosci. *19(7):2435–2445 under a different title.)*
9. Jones, M. V., & Westbrook, G. L. (1995). Desensitised states prolong GABA-A channel responses
   to brief agonist pulses. *Neuron*, 15(1), 181–191.
10. Kasugai, Y., Swinny, J. D., Roberts, J. D. B., et al. (2010). Quantitative localisation of
    synaptic and extrasynaptic GABA-A receptor subunits on hippocampal pyramidal cells by
    freeze-fracture replica immunolabelling. *Eur. J. Neurosci.*, 32(11), 1868–1888.
    *(Corrected from the previous draft, which gave* J. Neurosci. *30(42):14024–14035.)*
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
    *Neuroscience*, 49(1), 13–32. *(Corrected from the previous draft, which gave*
    J. Physiol. *454(1):477–496, a coordinate that does not resolve.)*
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
    doi:10.1016/j.ejphar.2016.09.016 *(replaces the previous draft's unresolvable "Fischer et
    al. 2010" attribution for MP-III-022)*
23. Walters, R. J., Hadley, S. H., Morris, K. D. W., & Amin, J. (2000). Benzodiazepines act on
    GABA-A receptors via two distinct and separable mechanisms. *Nat. Neurosci.*, 3(12),
    1274–1281. *(Corrected from the previous draft, which gave* Br. J. Pharmacol.
    *131(7):1307–1314 under a different title. This is the reference arm of every selectivity
    ratio in §4.4.)*

### Withdrawn citations

The following appeared in the previous draft and could not be resolved against Crossref by
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
