# Agent brief 20 — prove the positive direction, and retire the last two conventions

**Read first:** `knowledge/15-amplitude-constraint-theorem.md` (now a closed theorem),
`VERIFICATION-theorem-17.md`, `MANUSCRIPT-DRAFT.md` §2.2, §3.1, §5 and §8,
`RESULTS-BRIEF-18.md` §5. `PYTHONPATH=<your worktree>/src`.

**Your worktree:** branch `claude/agent-brief-20`, based on `main` (briefs 17, 18, 19 all
merged). One agent — no parallel work this round, the tasks touch overlapping files.

Briefs 17, 18 and 19 are all accepted. Between them they closed four of the manuscript's soft
spots. Verified independently before merging: the Duhamel cone identity (0 violations in 2000
checks), the secular identity (0 in 400), the analytic Hill-slope bracket (4 − 2√2 = 1.17157 and
the cooperative limit 2.000, both reproduced), and the 608-passed count that reconciled brief
19's apparent regression as xdist flakiness.

**Three gaps remain. This brief closes them or classifies them honestly.**

## 1. Task 1 — the positive direction is still a conjecture

The theorem proves single-open five-state schemes **cannot** produce a positive fast amplitude.
The complementary claim — that **two conducting open states can** — rests on a measurement
(a₁ = +1.116) and a narrative about "constructive interference" that §4 of the theorem document
now correctly labels a conjecture.

**Prove it, or find its boundary.** The machinery is already in hand and the asymmetry is
suggestive: in the single-open case the open state is a **leaf** of the tree, and the proof turns
on the gateway weight w_A2R being strictly negative. With two conducting states the observable is
`w^(O1) + w^(O2)`, a sum over two components of the same eigenvector, and the leaf argument no
longer applies to the summed projection.

Deliverables, ranked, and **say which you delivered**:

1. A theorem giving a **sufficient condition** on the rates under which a two-open scheme admits
   a₁ > 0 — ideally with the published 0.57 reachable inside it.
2. A **necessary** condition, i.e. a proof that a₁ > 0 requires at least two conducting states in
   the open manifold. This is arguably the more valuable direction: it would make "two open
   states are necessary" a theorem rather than an exhaustive-search result.
3. The precise inequality that remains open, written in the rate constants.

Do **not** re-label a numerical sweep as a proof. Brief 17 did that and it cost a round.

## 2. Task 2 — the band-to-σ convention, and what it does to the degrees of freedom

Two of the ten data are **bands** (`po_max` and mean open time), scored by treating the
half-width as a σ. The manuscript concedes this inflates `res_dof`, which is the difference
between `kinetic_coop` being a test (+1) and not being one.

* State what the half-width-as-σ convention assumes, and what a defensible alternative would be
  (a uniform prior over the band and its implied effective σ = half-width/√3 is the obvious
  candidate — **check whether that is right rather than taking my word**).
* Recompute the five-way minimax and `res_dof` under the alternative. **Does `kinetic_coop`
  still pass, and is it still a test?** If the answer changes, that is a material result and the
  manuscript needs it.
* Report both, and say which the project should adopt — then **stop**, because that choice is
  not yours.

## 3. Task 3 — the risetime prediction needs an error bar

`MANUSCRIPT-DRAFT.md` §5.1 offers the discriminating experiment as two bare numbers: 10–90%
risetime of 3.54 ms at c = 23 against 9.19 ms at c = 1000. **A prediction offered to
experimentalists without an uncertainty is not usable.**

Propagate the uncertainty properly: across the admissible parameter region at each `c` — not a
single optimum — report the predicted risetime as an interval. Then state the **minimum
resolvable difference**: given the spread, how precisely would an experimenter need to measure
risetime to distinguish c = 23 from c = 1000, and is that within reach of the ultrafast
application methods the source papers already use (10–90% exchange in 60–100 µs)?

If the intervals overlap, **say so** — the prediction would then not discriminate, and the
manuscript's §5.1 claim would have to be withdrawn.

## 4. Task 4 — a falsification statement the manuscript should carry

From brief 18's analysis, any two-site scheme in this family is confined to
n_H ∈ (4 − 2√2, 2.000) ≈ (1.172, 2.000), with open states unable to move it and cooperativity
the only parameter that does.

That is a **falsifiable boundary**: a Hill slope reliably measured above 2.0 for α1β2γ2S would
exclude the entire model class, not merely a parameter choice. The published 2.2 ± 0.4 straddles
it.

State this explicitly in `RESULTS-BRIEF-20.md`, with the bracket derived rather than asserted,
and note what measurement precision on n_H would be needed to make the test decisive.

## 5. Verify

Both tiers, **final pytest line pasted verbatim**. Note: `pytest-xdist` is active and has
produced inconsistent totals across runs — **run the tier twice and report both** if they
disagree, rather than reporting one. Baseline on `main` after the three merges: I measure 608
fast; confirm it. Do not touch or gate the Apple Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] Positive-direction result delivered at one of the three ranked levels, with the level stated explicitly.
- [ ] Band-to-σ alternative analysed; five-way minimax and `res_dof` recomputed under it; both reported.
- [ ] Risetime prediction given as intervals over the admissible region, with the minimum resolvable difference stated and compared to achievable exchange times.
- [ ] The n_H ∈ (1.172, 2.000) falsification boundary derived and stated.
- [ ] Both tiers pass; pytest lines verbatim; run twice if xdist disagrees.
- [ ] Commit per task; `kind=blocked` with what landed rather than running to the timeout.

## Not yours to decide

`kind=question`: which band convention the project adopts; any rewrite of the manuscript's
claims; adoption of any scheme; promotion of any finding to settled.
