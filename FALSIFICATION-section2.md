# Falsification — §2's structural claim as stated is BROKEN

**Found by:** Claude Opus 5, 2026-10-19, under the instruction that nothing publishable may be
disprovable. Recorded immediately, before the follow-up test completes, because the break is
already established and the manuscript must not sit in its current state.

## What breaks

`MANUSCRIPT-DRAFT.md` §2 claims that five-state single-open schemes cannot supply the measured
fast deactivation amplitude of **+0.57**, and §2.3 claims that observing three positive
amplitudes **necessitates** two conducting open states.

**Both claims compare a proven statement about eigenmodes with a measured quantity that is a
fit.** The published amplitudes (0.57 ± 0.04, 0.23 ± 0.02, 0.20 ± 0.03) are the output of a
**triexponential fit** to a recorded current. The theorems constrain the **eigenmode**
amplitudes of the generator. These are not the same quantity, and the difference is not small.

### Evidence 1 — a fit can return three positive amplitudes where the eigendecomposition has one

Searching 600 single-open parameter vectors, fitting a three-exponential to each model's own
waveform (as the experimenters did):

| fitted fractional amplitudes | fit RMS | positive **eigenmode** amplitudes |
|---|---|---|
| 0.4089, 0.2683, 0.3228 | **4.6 × 10⁻⁵** | **1** |
| 0.6591, 0.2917, 0.0492 | 7.6 × 10⁻⁴ | 2 |
| 0.0194, 0.6269, 0.3537 | 2.1 × 10⁻⁶ | 2 |

19 of 600 vectors gave three positive fitted amplitudes at RMS < 10⁻³. **So three positive
fitted components do not require two open states**, and §2.3's necessity argument does not
transfer from eigenmodes to the published data.

### Evidence 2 — an *admissible* single-open fit reaches a fast amplitude of 0.92

Using the project's **own canonical route** (`predict_components_fit`) and its **own
admissibility criterion** (well-conditioned: Jacobian condition ≤ 10⁴, no collapsed τ, no
vanishing amplitude), three single-open vectors pass admissibility with large positive fitted
fast amplitudes:

| fitted amplitudes | condition number | R² |
|---|---|---|
| **0.7150**, −0.0287, 0.3137 | 149 | 0.99818 |
| **0.9224**, −0.0155, 0.0931 | 5.19 × 10³ | 0.9999997 |
| **0.4576**, −0.0113, 0.5537 | 2.55 × 10³ | 0.9999990 |

The published value is 0.57. **It is bracketed.** The claim that single-open schemes cannot put
amplitude on the fast deactivation component is false as stated for the fitted quantity the
paper actually compares against.

## What may survive, and it is a different claim

In all three admissible counterexamples the **middle** amplitude is *negative* (−0.029, −0.016,
−0.011). The published data are **all three positive**. So the surviving candidate claim is not
about the magnitude of a₁ at all:

> An admissible single-open triexponential fit cannot return **three strictly positive**
> amplitudes.

That is the necessity theorem's content, restricted to admissible fits — and it is **not** what
§2 argues. A direct search for an admissible all-positive single-open fit is running; Evidence 1
suggests such fits exist but did not test them for admissibility, so the question is open.

## What the theorems still prove

Unaffected, because they are statements about the generator and were verified as such:

* `w_A2O/|w_A2R| = α/(|λ₁| − α) < 1/E`, unconditionally.
* The extremal **eigenmode** amplitude is non-positive on the reachable cone, with the Duhamel
  and secular identities verified numerically (0 violations in 2000 + 400 checks).

These remain correct. **What was wrong was the bridge from them to the measurement**, and that
bridge is mine: `DECISION-deactivation-objective.md` ruled the fit route canonical *precisely
because* the published numbers come from a fit, and I then allowed both theorems to be stated
about eigenmodes and compared against fitted values anyway. The decision document and the
manuscript contradicted each other for five rounds and I did not notice.

## Required manuscript action

1. **§2's headline claim must be withdrawn** in its current form. The eigenmode theorem cannot
   be used to exclude a fitted amplitude of 0.57.
2. **§2.3's necessity claim must be withdrawn** in its current form.
3. Whether a restricted claim survives — admissible single-open fits cannot be all-positive —
   depends on the running test. If it does, it must be stated as a claim about **admissible
   fitted decompositions**, proven or searched *in that space*, not inherited from the
   eigenmode theorems.
4. Every σ figure in §3 and §5 that was computed against **eigen-route** components needs
   recomputing on the canonical fit route, or relabelling.

Until 1–3 are resolved the paper has **no established structural claim**, and the two results
that remain intact are the methodological leverage finding and the Hill-slope bracket.

---

## Second falsification, same error class — the Hill-slope boundary (brief 22)

**The n_H > 2.000 falsification boundary is broken**, and for exactly the reason that broke §2:
a quantity derived in one definition, compared against a measurement made in another.

The bracket (4 − 2√2, 2.000) is the **local logarithmic derivative at EC50, at equilibrium**.
The paper reports — and Jahn 1997 measured — a slope fitted over a **finite 1–10 µM window on
peak currents**. Those are different quantities, and the difference is not marginal. Under the
windowed definition, unconstrained two-site models reach

* `kinetic_jw95`: n_H = **2.9325**
* `kinetic_coop`: n_H = **11.7284**

both far above the supposed 2.000 ceiling. Jahn's measurement also violates **both** hypotheses
of the bracket: peak currents are not equilibrium occupancies, and a finite-window regression is
not a local derivative.

**Consequence:** §5's claim that "a Hill slope reliably above 2.000 would exclude the entire
two-site class" is withdrawn. So is the associated precision requirement (σ_nH ≤ 0.067). The
bracket remains true of the quantity it describes — the equilibrium slope at EC50 — and that is
all it may be used for.

I reported this boundary as a headline falsification criterion. It was wrong.

## Third defect — the abstract contradicts the body

The abstract still asserts that the 10–90% risetime discriminates across the cooperativity
plateau. §5.1 and §8 withdrew that claim after the intervals were found to overlap completely
(minimum resolvable difference 0.00 ms). **I withdrew it in the body and left it in the
abstract.** Corrected.

## The pattern, stated plainly

Three claims have now broken, and all three broke the same way:

| claim | proven/derived about | compared against |
|---|---|---|
| §2 amplitude constraint | eigenmode amplitudes | a **fitted** amplitude |
| §2.3 necessity | eigenmode count | **fitted** component count |
| §5 Hill boundary | local derivative at EC50, at equilibrium | **windowed regression** on peak currents |

This is one error made three times: **establishing a result about one mathematical quantity and
testing it against a different measured quantity.** It is the same failure the paper's own
methodological section documents — scoring a reconstruction as data — in a more abstract
costume. The project's own `DECISION-deactivation-objective.md` names the rule that would have
caught all three: *the model must pass through the same pipeline as the measurement.* I wrote
that rule and then broke it three times in the same manuscript.

## What survives the audit

* **Theorem 2.2 as a statement about eigenmodes** — survives corner attacks: degenerate
  eigenvalues, vanishing off-diagonals, extreme ratios (10^±6), leaf and non-leaf topologies.
* **Theorem 2.3's necessity, restricted to admissible fits** — brief 22 reports 0 of 96 single-open
  configurations yielding three positive well-conditioned amplitudes; every attempt failed by
  amplitude collapse or τ coalescence. My own 3000-vector-per-scheme sweep is still running and
  is the larger test.
* **The fifth observable is robust to band choice** — `kinetic_coop` passes under the M-mode
  band (2.54σ), the H-mode band (2.39σ) **and** the full modal span (2.25σ). The modal span is
  the physically defensible choice for an ensemble measurement and removes the boundary pinning,
  so this is now *stronger* than when it rested on the narrowest band.
* **The leverage figures** are exact objective-function precision weights (n/σ²), not artefacts.
* **The 203.0625 floor** is the exact profiled minimum Mahalanobis distance on the rank-2
  2-simplex covariance — a live result, not a residue of the withdrawn sum-to-one objective.
* A §3 discrepancy (3.28σ in text vs 3.23σ in table) is resolved: different priors, both correct.
