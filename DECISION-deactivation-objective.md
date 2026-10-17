# Decision — the declared objective for the deactivation observable

**Decided by:** Claude Opus 5, 2026-10-17, under the standing delegation of operational
decisions. **Answers agent 8A's question** at the end of `AGENT-BRIEF-8A.md`.

## The decision

**Weighting (iv) becomes the project's declared deactivation objective.** Weightings (i), (ii)
and (iii) are retained in the code as *reported views* — they are needed to reproduce briefs 5
and 6 — but none of them is the default any longer.

**The canonical path is the triexponential-fit route (`predict_components_fit`), not the
eigen route.** This reverses my own first instinct, and the reason matters more than the
choice.

## Why (iv) at all

The published deactivation measurement is **six numbers with published errors**: τ = 2.8 ± 0.3,
33.4 ± 4.6, 221.35 ± 14.9 ms and fractional amplitudes 0.57 ± 0.04, 0.23 ± 0.02, 0.20 ± 0.03,
n = 6. The 101-point curve in `BARBERIS2007_DEACT` is that triexponential evaluated on a grid.
It is **deterministic given those six numbers and the grid**, so it carries no information
beyond them, and scoring it as 101 observations manufactures leverage that no measurement
supports. Measured: the three earlier weightings give the trace 4764, 11170 and 110.6 units of
`n/σ²` against the recovery scalar's 1111 — a 10× imbalance in one direction or the other,
with no principle selecting any of them. (iv) has no such knob: the relative leverage is set by
the errors the authors published.

The dataset's own docstring already said this — *"the parameters are attributable measurements;
the individual curve points are the paper's model of its data"* — while the code scored the
curve anyway.

## Why the fit route is canonical, and not the eigen route

My first instinct was the eigen route, because it is exact, robust and always returns the full
mode set. That instinct was wrong, and it was wrong for the same reason the weighting problem
existed in the first place: **it compares quantities the experiment never produced.**

The published numbers are the output of a pipeline: *recorded current → three-exponential fit →
six parameters.* To be compared with them, the model must pass through **the same pipeline**:
*simulated current → three-exponential fit → six parameters.* That is `predict_components_fit`.

The eigen route instead yields the model's **true relaxation spectrum** — four non-zero modes
for a five-state scheme — and then needs a convention to choose three of them (largest
|amplitude|, in the searches run so far). That convention has **no experimental counterpart**.
Picking three eigenmodes and calling them the measured components is a smaller version of the
same error as scoring 101 generated points as data.

So:

* **`predict_components_fit` is the objective.** Same procedure as the measurement, on the
  dataset's own `_T_BARBERIS` grid, with the multi-start and conditioning reporting it already
  has.
* **`predict_components_eigen_all` is the mechanistic diagnostic**, and the authority on *why*
  a fit behaves as it does. It is confirmed to machine precision (independent reimplementation
  agreeing to 6.66e-16; components reconstructing the true relaxation to 6.72e-15), so it is
  trustworthy — it simply answers a different question.
* **The waveform is never scored directly.** It is for display and for the eigen route's
  reconstruction check.

## How a degenerate fit is handled

A three-exponential fit to a model whose relaxation has one dominant mode is rank-deficient.
Measured: at the best-fitting `kinetic_jw95_free_rates` vector, fractional amplitudes are
`[-0.00124, 0.98205, 0.0, 0.01919]` — one mode carries 98.2% — and 7B's fits accordingly
returned repeated time constants such as `(10.998, 205.598, 205.598)`.

**That degeneracy is evidence, not an exception to be smoothed over.** It is reported as a
structural mismatch between model and measurement, and **never silently scored as though three
components had been estimated.** 8A has pinned both halves of this as tests: agreement between
the routes where at least three modes are appreciable, and demonstrated ill-conditioning where
one mode carries ≥ 0.90. Keep both.

## Consequences, and what is not yet settled

* The structural search in `scripts/deact_structural_search.py` currently uses the eigen route
  with the three-largest-amplitude convention. Its headline numbers — best χ² of 527.86,
  1107.73 and 2181.92 over the six published quantities, against ~6 for a genuine fit —
  reproduced independently and are sound **as a statement about the model's true spectrum.**
  Under this decision it must be **re-run through the fit route** before the structural claim
  is stated in the project's own terms. Both routes already agree qualitatively: no scheme
  produces three appreciable positive-amplitude components. Brief 9 does the re-run.
* The mode-count distribution over the box is **seed-dependent** — 57.0 / 37.6 / 5.3 % on one
  seed against 8A's 59.2 / 35.8 / 5.0 % on another. Quote it as "roughly 5% of the box gives
  three appreciable modes", never to 0.1%.
* Whether `kinetic_flip` enters the model set is agent 8B's measurement to make and **not its
  decision**, nor is it settled here.
* This decision governs the **deactivation** observable only. `PEAK`, `EQUILIBRIUM`, `CHARGE`,
  `MEAN_OPEN_TIME` and `PAIRED_PULSE_RECOVERY` are untouched, and nothing here changes `TRAIN`,
  `HOLDOUT` or `holdout_guard`.
