# Agent brief 16 — is c = 1000 required, or just where the optimiser landed?

**Read first:** `FINDING-coop-five-way.md`, `FINDING-two-open-states.md` with its addendum,
`RESULTS-BRIEF-15.md`, `DECISION-deactivation-objective.md` with its amendment,
`AGENT-BRIEF-6.md` §0. `FINDING-amplitude-obstruction.md` and `FINDING-minimax-conflict.md` are
withdrawn. `PYTHONPATH=<your worktree>/src`.

**Your worktree:** branch `claude/agent-brief-16`, based on `claude/agent-brief-15`. One agent.

Brief 15 is accepted in full. The mean-open-time derivation from the generator, the band
treatment of an unstated error without inventing a σ, the explicit modal-gating caveat, the
registration of three sources as abstracts, and the observation that **`koff2` = 0.0857 ms⁻¹
falls inside the JW95-propagated window** — that last one I had missed, and it is the best
corroboration this line of work has produced. Good.

## 1. The question

`kinetic_coop` satisfies all five observables at 2.58σ. Its single biggest liability is that it
needs **c = koff/koff2 = 1000**, against published site-to-site affinity ratios of ~23×
(Olsen 1981) and ~5× (Agey & Dunn 1989) — 43–200× beyond anything published.

**Nobody has asked whether c = 1000 is required.** It is where a minimax optimiser landed. The
question that decides whether this liability is real:

> How does the five-way minimax vary with `c` held fixed? Is there a value of `c` inside the
> published range that still satisfies all five observables?

If the fit is nearly as good at `c` = 20–50, the extreme value is an artefact of an
under-constrained optimum and **the liability largely evaporates**. If the minimax rises steeply
away from ~1000, then the model makes a **strong falsifiable prediction** — that α1β2γ2S has
site-to-site cooperativity two orders of magnitude beyond published estimates — and that is a
far more interesting result than a fit.

**Both outcomes are good. Neither is expected.**

## 2. Task 1 — profile `c`

Use the project's existing profile-likelihood machinery (see `scripts/koff_identifiability.py`
and the profile tooling it relies on; `fit_mle_vector` reports top-optima spread). For `c` on a
grid spanning at least **1 to 10,000** — log-spaced, at least 20 points, including the published
anchors **5 and 23** and the fitted **1000** — hold `log10(c)` fixed, re-optimise the remaining
8 parameters under the **five-way minimax** objective from brief 15, and report:

* the five-way minimax at each `c`;
* **which observable binds** at each `c`, since that is how you see what the cooperativity is
  buying;
* the full per-observable residual row at `c` = 5, 23, 100, 1000;
* whether any `c` ≤ 50 satisfies all five at ≤3σ;
* the profile's shape: is there a plateau, a sharp minimum, or a monotone ramp?

**Warm-start each grid point from the neighbouring solution**, not from scratch — a profile with
independent restarts at each point produces a ragged curve that is optimiser noise, and this
project has mistaken that for structure before. Report any non-monotonicity that survives, and
re-run those points from both directions before believing it.

**Report the runtime and budget per point.** The five-way objective is expensive; if the grid
would exceed ~40 minutes, reduce the number of points rather than the restarts, and say what you
used.

## 3. Task 2 — what does c = 1000 actually predict?

A cooperativity of 1000 is a physical claim with consequences beyond the five fitted
observables. At `koff` = 85.67 ms⁻¹ the singly-bound receptor dissociates in **~12 µs**, so it is
essentially never occupied at equilibrium. Derive and report, from the fitted model:

1. The **fractional occupancy of the singly-bound state** at its peak across the concentration
   range. If it is negligible everywhere, say so — that is the physical content of c = 1000.
2. The predicted **Hill slope as a function of c** along the profile. The CRC is where
   cooperativity shows up most directly, so this says whether c is doing the work the brief 13
   reasoning predicted.
3. Whether the model predicts anything **qualitatively falsifiable** — a measurement that would
   distinguish c = 1000 from c = 23. Name it concretely, with the observable and the expected
   difference. Do **not** search for or invent literature values for it; just state the
   prediction.

## 4. Task 3 — consolidate the record

This project now has findings in four states and a reader coming to it cold cannot tell what is
live. Write **`STATE-OF-RESULTS.md`** at the repo root, one page, stating for each claim: what it
is, its status, the number that supports it, the script that prints that number, and its main
caveat. Cover at minimum:

* the pseudo-data leverage result (**live**, `scripts/pseudodata_leverage.py`);
* the forced a₁ ≤ 0 in 5-state single-open schemes and the (0.57/0.04)² floor (**live**);
* `kinetic_coop`'s five-way fit (**provisional**, convention-dependent — state both the 2.58σ
  band figure and the 3.56σ point figure);
* `FINDING-amplitude-obstruction.md` (**withdrawn**, arithmetic identity);
* `FINDING-minimax-conflict.md` (**withdrawn**, recovery omitted from the objective);
* the `FIT_TARGETS['ec50_um']` correction (**resolved**, 20.0 → 11.6 ± 0.9 µM);
* what this brief's profile concludes about `c`.

**Do not restate claims more confidently than their source documents do.** Where a finding was
withdrawn, say what was wrong in one line. This file is the index, not a new argument.

## 5. Verify

Both tiers, **final pytest line pasted verbatim**. Brief 15 measured 587 fast / 661 full. Do not
touch or gate the Apple Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] `c` profiled over ≥20 log-spaced points spanning 1–10,000, including 5, 23 and 1000, warm-started along the grid.
- [ ] Five-way minimax and the binding observable reported at every point; full residual rows at c = 5, 23, 100, 1000.
- [ ] A plain answer: does any `c` ≤ 50 satisfy all five at ≤3σ?
- [ ] Profile shape characterised; any surviving non-monotonicity re-run from both directions.
- [ ] Singly-bound occupancy reported; Hill slope vs `c` reported; one concrete falsifiable prediction distinguishing c = 1000 from c = 23.
- [ ] `STATE-OF-RESULTS.md` covering all seven items, with status, number, script and caveat each.
- [ ] Runtime and budget per profile point stated.
- [ ] Both tiers pass; pytest lines verbatim.
- [ ] Commit per task; `kind=blocked` with what landed rather than running to the timeout.

## Not yours to decide

`kind=question`: adoption of any scheme into the canonical set; promotion of any finding to
settled; whether the project adopts a point or band treatment for unstated-error sources as its
default — report the consequence either way and leave the choice.
