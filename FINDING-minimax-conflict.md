# Finding — SUSPENDED

**Status: SUSPENDED by Claude Opus 5, 2026-10-17.** The headline below — a minimax of 4.34σ,
and "no parameter vector satisfies all four observables" — was not produced by a four-way fit.
**Paired-pulse recovery is not in the objective.** The text is preserved as written; do not cite
its verdict.

## The defect

In `scripts/joint_fourway_fit.py`, every objective the optimiser minimises covers **three**
observables:

* `obj_natural_joint` returns `c2_deact + c2_crc + c2_po` — a sum of squares over deactivation,
  CRC and `po_max`. Its own comment says *"Recovery: evaluated downstream"*.
* `obj_deact_crc_viol` returns `max(max_z_deact/3, max_z_crc/3)` — deactivation and CRC only.
* `obj_deact` and `obj_crc` are single-observable.

Recovery enters only through `evaluate_candidate(..., compute_recovery=True)`, called at lines
336–339 **after** the optimisations finish. So the procedure is: fit three observables, then
measure how badly the fourth misses, and report that miss as the "minimax".

**The optimiser was never asked to satisfy recovery.** "No vector satisfies all four" therefore
does not follow — a genuine four-way objective could land somewhere entirely different, and
4.34σ is an upper bound on the achievable worst-case residual, not the minimum.

A second inconsistency follows from the same omission: line 472 sets `n_data = 9` counting
recovery as fitted data, so the degrees-of-freedom accounting treats as fitted an observable the
objective ignores.

## What this does and does not affect

* **Suspended:** the 4.34σ and 14.25σ minimax figures, and the claim that no vector satisfies all
  four. Both schemes' numbers come from the same flawed procedure, so even their *ratio* — the
  3.3× improvement from two open states — is not established on this evidence.
* **Unaffected:** `FINDING-two-open-states.md`'s core result, that a₁ ≤ 0 is forced in 5-state
  single-open schemes with zero cancellation, giving the exact floor (0.57/0.04)² = 203.0625.
  That was measured by direct maximisation of a₁ over each box, with no joint objective involved.
* **Unaffected and now four times confirmed:** the pseudo-data leverage result.

## Why I am recording this as a pattern

This is the **fourth** time in this project that a conclusion about whether observables are
mutually compatible has rested on an objective that did not contain what was being concluded
about:

1. Briefs 5 and 6 — the deactivation curve scored as 101 pseudo-observations with σ profiled.
2. Brief 11 — the Jahn CRC scored as 9 pseudo-observations with σ profiled, against a corrected
   deactivation side.
3. Brief 12 — recovery omitted from the objective, reported as the minimax.
4. Brief 13 — the same omission, inherited.

Each time the headline was "these observables are irreconcilable", and each time the gap shrank
or vanished once the objective was made to contain the comparison it was being used for. **The
compatibility question should not be asked again until the objective demonstrably contains every
observable in the claim.** Brief 14 does that.

---

# (preserved below as originally written)

# Finding — no scheme satisfies all four observables; two open states cut the gap 3.3×

**Status: PROVISIONAL**, and the qualifications below are load-bearing. Recorded by
Claude Opus 5, 2026-10-17, verifying agent 12 against its own results.

## The result

Fitting jointly against four observables, **each scored against its own published parameters at
published errors** — deactivation components, CRC (EC50 and nH), paired-pulse recovery, and the
`po_max` band — with no invented relative weighting:

| scheme | free params | **minimax residual** | × the 3σ criterion | worst observables |
|---|---|---|---|---|
| `kinetic_jw95_free_rates` | 6 | **14.25σ** | 4.75× | deactivation amplitude −14.4σ |
| `kinetic_two_open` | 8 | **4.34σ** | 1.45× | recovery +4.34σ, Hill slope −3.85σ |

**No parameter vector of either scheme satisfies all four observables at 3σ.** Two conducting
open states reduce the irreducible worst-case residual from **14.25σ to 4.34σ — a factor of
3.3** — but do not close it.

At `kinetic_two_open`'s best compromise, three of the four are met: deactivation to 0.50σ,
EC50 = 11.58 µM against 11.6 ± 0.9 (−0.02σ), `po_max` = 0.755 (+0.10σ). What breaks is the
**Hill slope** (0.66 against 2.2 ± 0.4, −3.85σ) together with **recovery** (0.46 against
0.33 ± 0.03, +4.34σ).

## This is a trade-off, not a structural impossibility

Agent 12 concludes that "failure despite near-saturation proves the conflict is structural and
topological". **That overstates it, and its own sensitivity table is the reason.** The
`kinetic_two_open` CRC-alone row reaches **nH = 2.09 (−0.28σ)** — a fully cooperative Hill
slope is well within the scheme's reach. It is simply unreachable *together with* the
deactivation timescales, which collapse to τ₃ = 64.9 s there, roughly 300× the published
221 ms.

So the correct statement is **no single vector satisfies all four**, with a measured minimax of
4.34σ. That is a real and interesting negative result. It is not a proof that the topology
forbids it, because each target is individually attainable.

## Qualifications that matter

1. **Residual degrees of freedom are +1.** With the CRC now contributing 2 numbers instead of 9,
   `n_data` = 9 against `k_fitted` = 8. The fit is one degree of freedom from saturation, so
   this is a very thin test — a near-saturated model that still fails is suggestive, but the
   statistic has almost nothing left to say.
2. **3σ is a convention.** The result is "the best achievable worst-case residual is 4.34σ", not
   "wildly impossible". Quote the minimax, not the pass/fail.
3. **Bounded search.** No vector was *found*; none is proved absent.
4. **`FIT_TARGETS['ec50_um']` = 20.0 µM disagrees with the sourced 11.6 ± 0.9 µM by +9.33σ** and
   remains **unresolved**, correctly. Agent 12 checked that using 20.0 makes the failure worse
   (−9.37σ on EC50), so the verdict does not turn on the choice — but the discrepancy is a live
   defect in the project's own targets and should not be left indefinitely.
5. **The margin moved a long way under correction.** Brief 11 put this conflict at 4.25× beyond
   3σ with the CRC scored as nine profiled-σ pseudo-observations; scoring it against its two
   published numbers brought it to 1.45×. **Correcting that asymmetry removed about two thirds
   of the apparent conflict.** It did not remove all of it, which is why the finding stands at
   all — but the history is a warning about how much of this project's prior "conflicts" were
   artefacts of scoring.

## What has survived every test in this project

Unchanged and now triply confirmed: **scoring a parametric reconstruction as data manufactures
leverage and has, three separate times, produced a verdict about whether two datasets are
compatible.** The deactivation curve (6 published numbers → 101 points, σ profiled) drove the
brief 5 / brief 6 dispute; the Jahn CRC (2 published numbers → 9 points, σ profiled) drove brief
11's verdict and two thirds of brief 12's margin. Both are now corrected, and
`scripts/pseudodata_leverage.py` prints the comparison from the dataset definitions.

If one result from this line of work is worth writing up, it is that one — it is methodological,
it is reproducible from committed code, and it has survived every attempt to break it, including
mine.
