# Finding — the fast deactivation amplitude requires two conducting open states

**Status: PROVISIONAL.** Measured by Claude Opus 5, 2026-10-17, verifying agent 10 against its
own code. This is the **third** version of an amplitude-obstruction claim in this project; the
first two were artefacts and are withdrawn in `FINDING-amplitude-obstruction.md`. What is
different this time is stated under "Why this one is not the previous two".

## The result

Published target: the fastest deactivation component carries fractional amplitude
**a₁ = 0.57 ± 0.04** at τ₁ = 2.8 ± 0.3 ms (Barberis 2007, α1β2γ2S).

Maximising a₁ over each scheme's full parameter box, requiring the decomposition to be
**admissible** — `Σ|aᵢ| / |Σaᵢ| ≤ 1.5`, i.e. at most 50% cancellation between modes, which
excludes the pathologies that produced spurious values of 194, 9203 and 10.2 in earlier
attempts:

| scheme | states | open states | free params | **max admissible a₁** | cancellation | `po_max` there |
|---|---|---|---|---|---|---|
| `kinetic_jw95_free_rates` | 5 | 1 | 6 | **−0.00000** | **1.000** (none) | 0.9990 |
| `kinetic_open_desens_free_rates` | 5 | 1 | 6 | **−0.00000** | **1.000** (none) | 0.9990 |
| `kinetic_flip` | 6 | 1 | 8 | **+0.24515** | 1.500 (at the limit) | 0.8976 |
| `kinetic_two_open` | 6 | **2** | 8 | **+1.11574** | 1.500 (at the limit) | 0.0118 |

Three things follow, and they are different claims with different strengths.

1. **In both 5-state single-open schemes, a₁ ≤ 0 is forced**, with no cancellation at all
   (`Σ|aᵢ|/|Σaᵢ|` = 1.000 exactly). The published value is +0.57. This is a hard structural
   obstruction, and it produces an exact arithmetic floor on the amplitude objective of
   **(0.57/0.04)² = 203.0625** — measured as 203.06 in every single-open scheme.
2. **Adding a shut pre-open state relieves it only partially.** `kinetic_flip` is *not* forced
   to a₁ ≤ 0 — it reaches +0.245 — but that is still well short of +0.57 and requires 50%
   amplitude cancellation to achieve. Its best all-six χ² (273.86) is therefore
   indistinguishable from the 5-state schemes' (273.92, 273.95) **for a different reason than
   theirs**, which agent 10's report gets wrong — see "Corrections" below.
3. **Two conducting open states break the obstruction.** `kinetic_two_open` reaches a₁ = +1.12
   and fits all six published components **exactly**.

## The exact fit is not viable, and that is the headline caveat

At the vector where `kinetic_two_open` reproduces the six components exactly, agent 10's own
Task-4 cross-checks give:

| check | value | verdict |
|---|---|---|
| paired-pulse recovery | 0.3350 vs 0.33 ± 0.03 | **+0.17σ — matches** |
| Jahn 1997 peak CRC RMSE | 0.724604 vs baseline 0.004917 | **147.37× — fails** |
| `po_max` | 0.9355 | outside [0.70, 0.80] — **fails** |
| parameters on bounds | none | interior |

So two open states buy the deactivation **shape** and paired-pulse recovery simultaneously, and
then lose the concentration–response curve by two orders of magnitude. **This is the same wall
every previous reconciliation in this project has hit**, and it is now the central open
question rather than a footnote.

## Corrections to agent 10's report

* **The degrees-of-freedom accounting is invalid.** It reports χ² = 0.00 "against 5 degrees of
  freedom (p = 1.000)". `kinetic_two_open` fits **8 free parameters** to **5 independent data**
  (3 time constants + 2 independent amplitudes under the sum constraint), so the residual
  degrees of freedom are **5 − 8 = −3**. An exact fit is guaranteed by parameter count and
  **cannot** be a test of the model; p = 1.000 is meaningless here. What the result does
  establish is an **existence** claim — that such a vector exists in this scheme and does not
  in the single-open ones — and the contrast with `kinetic_flip`, which has the **same 8
  parameters** and still fails, is what makes that non-vacuous. That contrast is the real
  evidence, not the zero.
* **The stated mechanism is wrong for `kinetic_flip`.** The report says single-open schemes
  share "the exact same 203.06 arithmetic amplitude floor" because "the fastest observable
  relaxation mode remains forced to zero or negative amplitude". Measured above, flip reaches
  a₁ = +0.245, so it is **not** forced. It fails quantitatively, not structurally, and the two
  cases should not be collapsed.
* **Prose error:** §3 describes the two open states as "AR\* and A2R\*", implying a
  mono-liganded conducting state. The implementation is correct — both are doubly-liganded
  (`A2O`, `A2O'`), confirmed in `kinetic_two_open.py`. Only the prose is wrong.
* Reported test counts were "565 fast / 637 full". Independent measurement of the previous
  branch gave 561 fast where the agent had claimed 538 passed / 7 skipped; counts here are
  pending my own re-run. Three agents in a row have misreported counts while the code was
  green.

## Why this one is not the previous two

The first version rested on χ² = 452.99, which turned out to be the arithmetic minimum of a
sum-to-one normalisation and was identical across schemes for that reason. The second rested on
"the time constants fit exactly", which held only on the non-canonical route. This one has:

* **a conditioning criterion** that excludes the cancelling-amplitude regimes which generated
  every spurious number so far, applied uniformly;
* **a positive control** — two open states break the obstruction;
* **a negative control with matched parameter count** — `kinetic_flip`, 8 parameters, does not;
* **an exact arithmetic value** for the floor, (0.57/0.04)², derived rather than fitted;
* **zero cancellation** at the 5-state optima, so the forced a₁ ≤ 0 is not a numerical artefact.

It is still provisional, for these reasons:

1. **Bounded search.** "Forced" means no admissible vector was found in the declared box at
   this budget. A proof from the generator's structure would be worth far more than another
   search, and none has been attempted.
2. **The 1.5 cancellation bound is a convention I chose**, not a measured quantity. The 5-state
   results are insensitive to it (they sit at 1.000, nowhere near the limit) but `flip` and
   `two_open` both press exactly against it, so **their numbers depend on it** and would change
   if it were set elsewhere.
3. **The eigen route's three-of-four selection rule still applies** to the a₁ figures above.
4. **No scheme reconciles deactivation with the CRC**, so nothing here is yet a receptor model.

Numbers produced by `scripts/deact_obstruction_decomposition.py` (agent 10) and by the
admissibility scan that brief 11 is tasked with committing as
`scripts/amplitude_admissibility.py`.
