# Agent brief 11 — can two open states satisfy deactivation, recovery AND the CRC together?

**Read first:** `FINDING-two-open-states.md` (the current state, provisional),
`FINDING-amplitude-obstruction.md` (withdrawn — do not cite its claims),
`DECISION-deactivation-objective.md`, `BRIEF-7-SHARED.md`, `AGENT-BRIEF-6.md` §0. Discipline
rules unchanged. `PYTHONPATH=<your worktree>/src`. Interpreter
`/Users/whinz/Biochemistry/.venv/bin/python`.

**Your worktree:** branch `claude/agent-brief-11`, based on `claude/agent-brief-10`. One agent.

Your brief-10 work was good and most of it stands. Four corrections are in
`FINDING-two-open-states.md` under "Corrections to agent 10" — read them; tasks 1 and 2 below
are those corrections, and they are small. **Task 3 is the real work and the question this
whole project has been circling.**

## Task 1 — fix the degrees-of-freedom accounting everywhere

`RESULTS-BRIEF-10.md` reports `kinetic_two_open` at χ² = 0.00 "against 5 degrees of freedom
(p = 1.000)". It fits **8 free parameters** to **5 independent data** (3 time constants + 2
independent amplitudes under the sum constraint), so residual dof = **5 − 8 = −3**. An exact
fit is guaranteed by parameter count and cannot be a test; the p-value is meaningless.

Report, for every scheme and objective: **n_data, k_fitted, residual dof**, and a p-value
**only where residual dof > 0**. Where it is ≤ 0, say "not testable — k ≥ n" and leave it. The
project already has `k_effective` and `at_bound` for discounting bound-pinned parameters; use
them.

**The real evidence is the contrast, not the zero:** `kinetic_flip` has the same 8 parameters
and still fails. State it that way.

## Task 2 — two corrections to the mechanism, and make admissibility part of the objective

**(a)** Your §3 says single-open schemes share the 203.06 floor because "the fastest observable
relaxation mode remains forced to zero or negative amplitude". That is true for the two 5-state
schemes and **false for `kinetic_flip`**, which reaches a₁ = **+0.24515**. Flip fails
quantitatively (0.245 against 0.57 required, and only with 50% cancellation), not structurally.
Correct it; do not collapse the two cases.

**(b)** Your §3 describes the two open states as "AR\* and A2R\*". The implementation is right —
both are doubly-liganded `A2O`, `A2O'`. Prose only.

**(c)** Promote the admissibility criterion into the objective, as a committed script
`scripts/amplitude_admissibility.py` and a gate in `deact_components`: a decomposition is
admissible only if `Σ|aᵢ| / |Σaᵢ| ≤ L1_MAX`, default 1.5. Every earlier spurious number in this
project — 194, 9203, 10.2, and my own 1.24 — came from cancelling-amplitude regimes that passed
the existing `Σa ≈ P_o(0)` check, because the sum can be 1.0 while individual amplitudes are
±10⁴.

**Report the sensitivity**, because this bound is a convention I chose and not a measurement.
Sweep `L1_MAX` over at least {1.0, 1.25, 1.5, 2.0, 5.0} and tabulate max admissible a₁ per
scheme. The 5-state results sit at 1.000 and should be insensitive; `kinetic_flip` and
`kinetic_two_open` press exactly against 1.5 and their numbers **will** move. Say by how much.
If the qualitative conclusion — 5-state forced ≤ 0, two-open exceeds 0.57 — does not hold
across that whole sweep, **say so plainly**; it would mean the finding depends on my
convention.

## Task 3 — the question that matters

`kinetic_two_open` fits the six deactivation components exactly and matches paired-pulse
recovery at +0.17σ, and then misses the Jahn CRC by **147.37×** its baseline with
`po_max` = 0.9355. Every reconciliation attempted in this project has died exactly there.

**So: is there any parameter vector of `kinetic_two_open` that satisfies all three
simultaneously?**

Fit jointly against:

1. the six published deactivation components, at their published errors, under the corrected
   2-simplex objective and the admissibility gate;
2. paired-pulse recovery, 0.33 ± 0.03 (`BARBERIS2007_PAIRED_PULSE`);
3. the Jahn 1997 peak CRC (`JAHN1997_PEAK_CRC`), with its baseline RMSE of 0.004917;
4. the `po_max` gating anchor, band [0.70, 0.80].

**Weighting is the hazard here and this project has already been burned by it twice.** The CRC
is a 9-point parametric reconstruction; the deactivation components are 6 published numbers
with published errors; recovery is 1 scalar with a published SEM. **Do not invent a weighting.**
Score each dataset by its own published errors where they exist, exactly as
`DECISION-deactivation-objective.md` requires, and **report the per-dataset residuals
separately in σ** rather than only a total. If you find yourself choosing a relative weight,
stop and report `kind=question`.

Then answer in `RESULTS-BRIEF-11.md`, per scheme and with every residual shown:

* Is there a vector satisfying all four within 3σ (CRC within 3× baseline RMSE)? **Yes or no.**
* If yes: the full vector, every bound flag, and `k_effective`. Then **hold** — adoption is not
  yours.
* If no: **which pair is irreconcilable, and by how much.** That is the project's result, and it
  is a perfectly good one.

Also run the same joint fit for `kinetic_jw95_free_rates` as the reference, so the comparison
comes from one code path and one budget.

## Task 4 — verify

Both tiers, `PYTHONPATH` set, **paste the final pytest line verbatim, do not transcribe it**.
Measured baseline on `claude/agent-brief-9`: **561 fast, 0 skipped**. You reported 565 fast /
637 full on brief 10; I have not re-measured those. Do not touch or gate the Apple Silicon
split-R̂ artefact.

## Acceptance criteria

- [ ] dof reported as n_data, k_fitted, residual dof everywhere; p-values only where residual dof > 0.
- [ ] The flip mechanism corrected; the AR\*/A2R\* prose corrected.
- [ ] `scripts/amplitude_admissibility.py` committed; L1 gate in `deact_components`; sensitivity sweep over at least 5 values of `L1_MAX` tabulated per scheme.
- [ ] Explicit statement of whether the qualitative conclusion survives the whole sweep.
- [ ] Joint four-way fit for `kinetic_two_open` and `kinetic_jw95_free_rates`, same code path and budget, **per-dataset residuals in σ**, no invented weighting.
- [ ] `RESULTS-BRIEF-11.md` answers Task 3's question yes or no, and if no, names the irreconcilable pair and the margin.
- [ ] Both tiers pass; final pytest lines pasted verbatim.
- [ ] Commit per task; `kind=blocked` with what landed rather than running to the timeout.

## Not yours to decide

`kind=question`: adoption of `kinetic_two_open` into the canonical model set; promotion of
`FINDING-two-open-states.md` from provisional to settled; any choice of relative weighting
between datasets.
