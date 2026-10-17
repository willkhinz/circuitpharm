# Agent brief 14 — put recovery in the objective; and the two EC50 decisions, answered

**Read first:** `FINDING-minimax-conflict.md` — **now SUSPENDED**, and the reason is Task 1.
Also `FINDING-two-open-states.md`, `DECISION-deactivation-objective.md` with its amendment,
`AGENT-BRIEF-6.md` §0. `FINDING-amplitude-obstruction.md` is withdrawn.
`PYTHONPATH=<your worktree>/src`.

**Your worktree:** branch `claude/agent-brief-14`, based on `claude/agent-brief-13`. One agent.

Your brief-13 work is accepted and both your questions are answered in Task 3. Your `kinetic_coop`
implementation, the reversibility argument (acyclic tree, Kolmogorov) and the nesting check at
c = 1 to 1e-5 are all good. **But its headline comparison, and brief 12's, rest on a defect I
missed when I wrote both briefs.**

## Task 1 — the four-way fit has never included the fourth observable

In `scripts/joint_fourway_fit.py`, every objective the optimiser minimises covers **three**
observables:

* `obj_natural_joint` → `c2_deact + c2_crc + c2_po`, with the comment *"Recovery: evaluated
  downstream"*.
* `obj_deact_crc_viol` → `max(max_z_deact/3, max_z_crc/3)`, deactivation and CRC only.

Recovery enters only via `evaluate_candidate(..., compute_recovery=True)` at lines 336–339,
**after** optimisation. So the procedure fits three observables and then reports how badly the
fourth misses — and that miss is what both briefs called the "minimax". The optimiser was never
asked to satisfy recovery, so **4.34σ is an upper bound, not a minimum**, and "no vector
satisfies all four" does not follow. `kinetic_coop`'s 17.69σ has the same status.

Note also line 472: `n_data = 9` counts recovery as fitted data while the objective ignores it.

**Fix it:**

1. Add recovery to the joint objective — `z_pp² = ((R − 0.33)/0.03)²`, read from
   `BARBERIS2007_PAIRED_PULSE`, no literals — so `obj_natural_joint` covers all four.
2. Add a **true four-way minimax** objective: `max` over the worst standardised residual of
   *every* observable, including recovery, and optimise **that** directly. The existing Pareto
   objective is a two-observable minimax and should be kept and relabelled as such, not deleted.
3. **Assert the nesting property**, since it now applies and this project has twice been misled
   by its absence: `kinetic_coop` reduces to `kinetic_two_open` at `c = 1`, which you verified
   to 1e-5, so on **any** objective optimised over both boxes, `kinetic_coop`'s optimum must be
   **≤** `kinetic_two_open`'s. Warm-start `kinetic_coop` from `kinetic_two_open`'s optimum
   embedded at `c = 1` via `nesting_embedding`, and **refuse to report a worse value** — raise,
   as `dataset_consistency.py` already does. A nested model scoring worse is a failed
   optimisation, not a finding.
4. Fix `n_data` to count exactly what the objective contains, and report `n_data`, `k_fitted`,
   residual dof per scheme as you now do.

Then re-run for `kinetic_jw95_free_rates`, `kinetic_two_open` and `kinetic_coop`, same code path
and budget, and answer in `RESULTS-BRIEF-14.md`, per scheme, every residual in σ:

* **The true four-way minimax** for each scheme. Compare to the suspended 14.25σ / 4.34σ /
  17.69σ and **say how much of each was the omission.**
* Is there a vector satisfying all four within 3σ? **Yes or no**, and this time the objective
  contains the question.
* If no: which pair is irreconcilable **under an objective that includes both of them**.
* Whether the nesting assertion held, and for which schemes.

**Do not invent a weighting.** Every quantity has a published error.

## Task 2 — expect the answer to change, and say so plainly

If including recovery substantially lowers the minimax, then `FINDING-minimax-conflict.md` is
**withdrawn rather than merely suspended** — update it in place, preserving what it claimed. If
the gap survives, it becomes a much stronger result, because it will then be the first version
of this claim where the objective contained every observable in it.

Either way: state in `RESULTS-BRIEF-14.md` whether brief 12's and brief 13's verdicts stand, and
add superseded headers to `RESULTS-BRIEF-12.md` and `RESULTS-BRIEF-13.md` pointing here. Numbers
preserved; this project supersedes rather than deletes, and I have applied that to two findings
of my own.

## Task 3 — your two EC50 questions, answered

**Question 1 — `tests/test_manuscript_consistency.py`, 11 tests pinned to the 20 µM anchor.**

**Update the manuscript. Do not decouple it from `FIT_TARGETS`.**

The manuscript's numbers were derived from an **unsourced** value — your own register records
20.0 µM as *"not stated / recalled literature middle from earlier legacy tuning files"*.
Decoupling manuscript reproduction from `FIT_TARGETS` would let the manuscript keep quoting
figures from an anchor the project has rejected, and hide the dependency that makes those
figures wrong. The coupling is doing its job: it caught this.

So: recompute the manuscript's affected numbers from the sourced 11.6 ± 0.9 µM, update
`knowledge/08-manuscript.md`, and **add a prominent, dated note in the manuscript recording that
the previous version was anchored on an unsourced 20 µM, what changed, and by how much** —
including the headroom moving 123.3× → 43.3× and the rate constants you listed. Keep the old
figures visible in that note. A manuscript that silently changes its numbers is worse than one
that shows its corrections.

**Question 2 — `test_falsification_intervals_are_set_from_reachable_gain_not_the_asymptote`,
asserting `reach[0.1] < 0.01 * head[0.1]`, now 1.78%.**

**Re-derive the interval; do not widen the threshold.** But read this carefully, because the two
can look identical in a diff and only one is legitimate.

The 1% constant was calibrated when headroom was 2881×. The sourced anchor collapses it to 650×,
so the constant was fitted to a quantity derived from the wrong EC50. Re-deriving the
falsification interval from the sourced anchor is correcting a derived value after fixing its
input. **Changing `0.01` to `0.02` to make the assertion pass is widening a gate, and the
standing rule forbids it.**

Concretely: work out what the test is actually asserting — that falsification intervals come from
*reachable* gain rather than the asymptote — and express the threshold in terms of the quantity
that claim depends on, so that it is derived from the anchor rather than hard-coded beside it. If
you cannot do that, **report `kind=question` with the derivation you attempted.** Do not ship a
changed constant with a comment saying the anchor moved.

## Task 4 — verify

Both tiers, **final pytest line pasted verbatim**. Baselines measured by me on `main`:
**574 fast / 646 full**, zero failures. Do not touch or gate the Apple Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] Recovery in `obj_natural_joint`; a true four-way minimax objective added and optimised; the two-observable Pareto objective kept and relabelled.
- [ ] Nesting warm-start from `kinetic_two_open` at `c = 1`, with a refusal to report a worse nested optimum.
- [ ] `n_data` counts exactly what the objective contains; dof reported per scheme.
- [ ] `RESULTS-BRIEF-14.md` gives the true four-way minimax per scheme, says how much of the old figures was the omission, and answers the all-four question with an objective that contains it.
- [ ] `FINDING-minimax-conflict.md` withdrawn or strengthened, updated in place either way.
- [ ] Superseded headers on `RESULTS-BRIEF-12.md` and `RESULTS-BRIEF-13.md`.
- [ ] Manuscript updated from the sourced anchor, with a dated note preserving the old figures.
- [ ] Falsification threshold **derived**, not widened — or `kind=question` with the attempted derivation.
- [ ] Both tiers pass; pytest lines verbatim.
- [ ] Commit per task; `kind=blocked` with what landed rather than running to the timeout.
