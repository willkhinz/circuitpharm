# Agent brief 10 — the amplitude objective is misspecified; rebuild it, then re-ask everything

**Read first:** `BRIEF-7-SHARED.md`, `AGENT-BRIEF-6.md` §0, `DECISION-deactivation-objective.md`,
and `FINDING-amplitude-obstruction.md` — **which is now WITHDRAWN**, including its own
"what survives" section. Discipline rules unchanged. `PYTHONPATH=<your worktree>/src` for all
test runs. Interpreter `/Users/whinz/Biochemistry/.venv/bin/python`.

**Your worktree:** branch `claude/agent-brief-10`, based on `claude/agent-brief-9`, so agent 9's
warm-start work, conditioning checks and `kinetic_two_open` are all present. **One agent.**

## 1. What went wrong, so you do not rebuild it

Agent 9's work was sound; the **objective it was measuring with is not**. Two defects, both
mine, found after its report:

**(a) The amplitudes-only χ² has an arithmetic floor.** Fractional amplitudes sum to 1 by
construction. Any model that zeroes its fastest selected amplitude is left with `(0, x, 1-x)`,
and minimising the published χ² over that single number gives **452.9856 at x = 0.4054**, of
which `(0.57/0.04)² = 203.06` is the zeroed fast term alone. Measured: **452.99 in every
scheme**, matching to 0.0044. The invariance across topologies that was read as evidence of a
structural constraint is just that floor.

**(b) The amplitude term is misspecified.** The published amplitudes 0.57, 0.23, 0.20 **sum to
exactly 1.0**. They carry **two** independent degrees of freedom, not three. The objective
scores all three as independent measurements with independent errors, and the published errors
(0.04, 0.02, 0.03) cannot all be independent either. Scoring a constrained composition as
unconstrained double-counts and inflates χ².

Consequence: of the headline all-six χ² of 527.86, **452.99 is floor and only 74.87 is
anything else**. Every χ² in `RESULTS-BRIEF-9.md` and `RESULTS-BRIEF-8A.md` inherits this.

**Agent 9's verdict that the finding stands was rejected**, for the same reason 7A's was: it
summarised per-scheme numbers as ranges ("amplitudes-only 71.06 to 1031.15; taus-only 106.77 to
123.19"), which hides that for the reference scheme the amplitude number is the *lower* one.
**State every decomposition result per scheme. Never as a range across schemes.**

Also, for the record and not a criticism of the code: your reported counts were **538 passed,
7 skipped**. Measured here on your own branch with `PYTHONPATH` set: **561 passed, 0 skipped,
0 failed**. The code was green and better than claimed. Three agents in a row have now
misreported test counts while the code was fine — run the tier, paste the final line verbatim,
do not transcribe.

## 2. Task 1 — specify the amplitude comparison correctly

This is the whole brief. Decide and implement a defensible treatment of a **constrained
composition** measured with stated per-component errors, and justify it in the code and in your
report. Options, and you may propose better:

* **Drop one component.** Score `(a₁, a₂)` only, with `a₃ = 1 - a₁ - a₂` implied. Simple,
  honest, but the choice of which to drop is arbitrary and the result must be shown to be
  invariant to it — **check all three droppings and report whether χ² depends on the choice. If
  it does, this option is wrong.**
* **Compositional transform.** Score additive log-ratios, e.g. `log(a₁/a₃)` and `log(a₂/a₃)`,
  propagating the published errors through the transform. Handles the constraint properly but
  is undefined for non-positive amplitudes, which these models produce — **say what you do
  then, and do not silently drop those vectors.**
* **Full covariance.** Treat the three as a degenerate Gaussian with a rank-2 covariance
  implied by the sum constraint, and use a pseudo-inverse. Most principled; state the
  covariance you assume, because the paper does not publish one.

Whatever you choose: **it must have no arithmetic floor.** Demonstrate that explicitly —
construct the amplitude triple the objective most prefers and show its χ² is of order the
degrees of freedom, not 450. Add a test pinning that.

## 3. Task 2 — re-derive everything under the corrected objective

Per scheme, both routes, all three objectives (all six / τ only / amplitudes only), at a budget
you state and defend. The warm start you built gives ~0.9 s per fit evaluation, so this is now
affordable — **use a real budget this time**, not 50 LHS with 2 polished in a six-dimensional
box.

Schemes, in priority order: `kinetic_jw95_free_rates`, `kinetic_open_desens_free_rates`,
`kinetic_flip`, `kinetic_two_open`, `extended_desens`.

Then answer, in `RESULTS-BRIEF-10.md`, plainly and per scheme:

1. Under the corrected amplitude objective, **does any scheme reproduce the published
   three-component deactivation?** Give χ² against its degrees of freedom.
2. Is the failure, where it occurs, in the time constants or the amplitudes? **Per scheme.**
3. Does `kinetic_flip`'s non-improvement over `kinetic_jw95_free_rates` survive? It was 527.41
   against 527.86 — both dominated by the floor, so that comparison has to be redone.
4. Does `kinetic_two_open` — two conducting open states — differ from the others? It has never
   been scored on a sound objective.

**Supersede, do not silently replace.** `RESULTS-BRIEF-8A.md` and `RESULTS-BRIEF-9.md` keep
their numbers; add a header to each saying the objective they used was later found to have an
arithmetic floor, and pointing at this brief's results. Preserving what was claimed is the
house rule and I have just applied it to my own finding.

## 4. Task 3 — the one claim that has survived everything

Untouched by any of this: scoring a 101-point reconstruction as data, when it is deterministic
given six published numbers, manufactured a 10× leverage imbalance (measured: 4764, 11170 and
110.6 units of `n/σ²` against the recovery scalar's 1111) and drove the entire brief 5 / brief 6
dispute about whether two datasets contradict each other.

Write `scripts/pseudodata_leverage.py`, committed, that prints that comparison from the dataset
definitions themselves — the point count, the σ under each weighting, the resulting `n/σ²`, and
the ratio — reading everything from `circuitpharm.fitting.data` and
`circuitpharm.literature`, **no literals**. This is the project's most durable result and it
currently exists only in prose and in my commit messages.

## 5. Verify

Both tiers, `PYTHONPATH` set, **final pytest line pasted verbatim**, counts labelled and
accounted for against **561 fast** (measured by me on `claude/agent-brief-9`) and the full tier
you measure. Do not touch or gate the Apple Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] A corrected amplitude objective, with the constraint handled, justified in code and report.
- [ ] Explicit demonstration that it has **no arithmetic floor**, plus a test pinning that.
- [ ] If the drop-one option is used, invariance to which component is dropped is checked and reported.
- [ ] Decomposition re-derived per scheme, both routes, three objectives, at a stated real budget.
- [ ] `RESULTS-BRIEF-10.md` answers §3's four questions plainly and **per scheme, never as ranges**.
- [ ] `RESULTS-BRIEF-8A.md` and `RESULTS-BRIEF-9.md` given superseded headers, numbers preserved.
- [ ] `scripts/pseudodata_leverage.py` committed, no literals.
- [ ] Both tiers pass; final pytest lines pasted verbatim.
- [ ] Commit per task. Report `kind=blocked` with what landed rather than running to the timeout.

## Not yours to decide

`kind=question`: adoption of `kinetic_two_open` or `kinetic_flip` into the canonical model set,
and any promotion of a finding from provisional to settled. Measuring is yours.
