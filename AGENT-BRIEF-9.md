# Agent brief 9 — re-derive the obstruction under the canonical objective, then test the two untested mechanisms

**Read first, in full:** `BRIEF-7-SHARED.md`, `AGENT-BRIEF-6.md` §0,
`DECISION-deactivation-objective.md`, and `FINDING-amplitude-obstruction.md`. The discipline
rules are unchanged. Run tests with `PYTHONPATH=<your worktree>/src`. Interpreter
`/Users/whinz/Biochemistry/.venv/bin/python`.

**Your worktree:** already created, branch `claude/agent-brief-9`, based on
`claude/agent-brief-8a`. **One agent this time** — the tasks are strictly ordered and parallel
work would collide on `comparison.py` and `literature.py`.

---

## Where things stand

`FINDING-amplitude-obstruction.md` records the current result and is marked **PROVISIONAL**.
The short version: fitting the published deactivation **time constants alone** gives χ² = 0.00
in both `kinetic_jw95_free_rates` and `kinetic_flip`, while fitting the **amplitudes alone**
gives χ² = 452.99 in both, and all six together gives 527.86 in both. The obstruction is in the
amplitudes, and adding a pre-open state buys 0.45 in 527.86 for two extra parameters. At 6×
search budget and a different seed the numbers did not move.

**Why it is provisional, and it is my reservation about my own result:** those numbers come
from the eigen route, which yields four non-zero modes for a five-state scheme and then selects
three by largest |amplitude|. That selection rule has no experimental counterpart, and
`DECISION-deactivation-objective.md` has already ruled the **triexponential-fit route**
canonical, because the published numbers are the output of *current → three-exponential fit →
six parameters* and the model must pass through the same pipeline. **The identical χ² across
schemes is suggestive of a real constraint, but it is also exactly what a binding selection
rule would produce.** Task 2 settles which.

## Tasks, in this order

### Task 1 — merge 8B, and reconcile the literature field name

`git merge claude/agent-brief-8b`. It will conflict: 8A renamed the field to
`component_errors` with a `component_errs` alias in `__post_init__`, while 8B is based on 7A
and uses `component_errs`. Keep **8A's** choice (`component_errors` canonical, `component_errs`
aliased) and migrate 8B's call sites. Both will also touch `comparison.py` — 8B added the
`kinetic_flip` spec; keep it.

After merging, confirm both field names resolve and that the `note` prose is intact, and that
`scripts/flip_structural_search.py` and `scripts/deact_structural_search.py` both still run.

### Task 2 — re-derive the obstruction through the fit route

Write `scripts/deact_obstruction_decomposition.py`, committed, reporting for **each** of
`kinetic_jw95_free_rates`, `kinetic_flip`, `extended_desens` and `kinetic_open_desens`:

* best χ² fitting **all six** published quantities,
* best χ² fitting the **three time constants only**,
* best χ² fitting the **three amplitudes only**,

each computed **twice** — once via `predict_components_eigen_all` with the three-largest-
|amplitude| rule, and once via `predict_components_fit`, the canonical route. Read every
published value and error from `circuitpharm.literature.DEACTIVATIONS[1]`; **no literals**.

**Where the fit route degenerates, report it as degenerate** — `is_well_conditioned` false, or
repeated τ — and do **not** substitute an eigen-route number for it. A scheme whose waveform
cannot support a three-exponential fit has failed to reproduce the measurement, and that is the
finding for that cell, not a missing value to be filled.

Then answer plainly, in `RESULTS-BRIEF-9.md`:

> Under the canonical fit route, is the obstruction still in the amplitudes rather than the
> time constants, and is it still unrelieved by adding states?

**If the fit route disagrees with the eigen route, say so prominently and do not reconcile them
by preferring the one that matches the existing finding.** The disagreement would mean the
provisional finding was an artefact of the selection rule, and
`FINDING-amplitude-obstruction.md` must then be withdrawn — update it, in place, saying so.

### Task 3 — test the two mechanisms that have NOT been tested

The amplitude reading in `FINDING-amplitude-obstruction.md` predicts that what is needed is a
mechanism that **removes open channels quickly**, not another shut state upstream. Two
candidates, neither yet tested:

1. **Open-state desensitisation with its rates actually freed.** `kinetic_open_desens` exists
   (`A2R ↔ A2O ↔ A2D`) but its spec has only **3** free parameters — `koff`, `alpha` and `r`
   are still pinned at `FIT_CONVENTIONS`, so the topology has never had a fair test. Add
   `kinetic_open_desens_free_rates` with the same 6 free parameters as
   `kinetic_jw95_free_rates`, additive, and run it through Task 2's decomposition.
2. **A second open state.** Add `src/circuitpharm/models/kinetic_two_open.py`:
   `R ↔ AR ↔ A2R ↔ A2O ↔ A2O'`, with `A2R ↔ A2D` desensitisation as in JW95, where both
   `A2O` and `A2O'` are conducting. The observable is the **sum** of the two open-state
   occupancies — be explicit about this in the docstring and make sure the component predictors
   receive it correctly, since `open_index` is a single integer and this scheme needs two.
   **If that requires extending the predictor's interface, say so and do it cleanly** rather
   than hard-coding an index. Connectivity test with explicit zero assertions, as in
   `tests/test_nextgen_models.py`.

Declare nesting honestly for both, using `nests` / `nesting_embedding`. If a limit sits on a
bound or does not hold exactly, record the measured discrepancy and exclude it — as was done
for `extended_desens` (16.66 nats) and `kinetic_flip` (161.60 nats).

### Task 4 — if anything fits, check it against the rest of the evidence

Any scheme reaching χ² of order tens rather than hundreds on the six published numbers must
then report, at that same optimum:

1. Paired-pulse recovery against 0.33 ± 0.03 (`BARBERIS2007_PAIRED_PULSE`), in σ.
2. The Jahn 1997 peak CRC RMSE against its baseline of 0.004917. **This is where every previous
   reconciliation died** — brief 6's missed by 97.6× and 124.4×, and `kinetic_flip`'s
   deactivation optimum by 115×.
3. `po_max` against `FIT_RANGES` [0.70, 0.80]. **0.8004 is outside the band**; an earlier report
   called it "aligned with" it, and `kinetic_flip` gave 0.2366.
4. Every parameter sitting on a bound, flagged. A pinned parameter is the edge of the search
   box, not an estimate.

Report all four whichever way they fall. **Do not report a favourable subset** — that error has
now been made twice in this project, in both directions.

### Task 5 — verify

Both tiers, `PYTHONPATH` set, counts stated separately and labelled. Baselines: 546 fast / 618
full on `claude/agent-brief-8a`; 8B measured 544 fast / 616 full on its own branch, so expect
roughly 546–560 fast after the merge plus your additions, and **account for the number** rather
than just reporting it. Do not touch or gate the known Apple Silicon split-R̂ artefact in
`test_posterior_recovers_known_parameters`.

## Acceptance criteria

- [ ] 8B merged; `component_errors` canonical with `component_errs` aliased; `note` prose intact; both existing search scripts still run.
- [ ] `scripts/deact_obstruction_decomposition.py` committed, covering 4 schemes × 3 objectives × 2 routes, with no published value hardcoded.
- [ ] Degenerate fit-route cells reported as degenerate, never backfilled from the eigen route.
- [ ] `RESULTS-BRIEF-9.md` answers the Task 2 question plainly, and states explicitly whether `FINDING-amplitude-obstruction.md` stands or is withdrawn.
- [ ] If withdrawn, `FINDING-amplitude-obstruction.md` is updated in place to say so, preserving what it originally claimed.
- [ ] `kinetic_open_desens_free_rates` spec added; `kinetic_two_open.py` added with connectivity test and summed-open observable handled without a hard-coded index.
- [ ] Nesting declared or excluded with measured discrepancies for both new entries.
- [ ] Task 4's four checks reported for any scheme that fits, all of them.
- [ ] Both tiers pass, counts labelled and accounted for.

## Not yours to decide

Report `kind=question` rather than deciding: whether any new scheme is **adopted** into the
project's model set, and whether `FINDING-amplitude-obstruction.md` is promoted from provisional
to settled. Measuring is yours; both of those are the planner's and the user's.
