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
