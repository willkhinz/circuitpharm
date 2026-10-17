# Agent brief 12 — score the CRC against its two published numbers, then re-ask the four-way question

**Read first:** `DECISION-deactivation-objective.md` **including the 2026-10-17 amendment at
the end**, `FINDING-two-open-states.md`, `RESULTS-BRIEF-11.md`. `FINDING-amplitude-obstruction.md`
is withdrawn — do not cite it. Discipline rules per `AGENT-BRIEF-6.md` §0.
`PYTHONPATH=<your worktree>/src`.

**Your worktree:** branch `claude/agent-brief-12`, based on `claude/agent-brief-11`. One agent.

Your brief-11 work stands and the corrections you made are accepted. **But your "no" answer is
suspended**, for a reason that is my fault, not yours.

## 1. The asymmetry

You answered Task 3 "no": deactivation and the Jahn CRC are irreconcilable by at least 4.25×
beyond 3σ, with compromises blowing recovery to +14.4σ–+19.5σ. You reached that with the
**deactivation** side scored against its six published parameters at published errors — the
corrected objective — and the **CRC** side scored as **nine parametric pseudo-observations with
σ profiled out**:

```
JAHN1997_PEAK_CRC:  kind=parametric, synthetic=True, n=9, noise_model='estimated'
```

generated in `fitting/data.py` from **two** published numbers:

```
_JAHN_EC50_UM = 11.6   _JAHN_EC50_SEM = 0.9
_JAHN_HILL    = 2.2    _JAHN_HILL_SEM = 0.4
```

So the verdict compares a corrected dataset against an uncorrected one, and the uncorrected one
is the one that wins the trade. `DECISION-deactivation-objective.md` originally scoped weighting
(iv) to deactivation alone; **that scoping was my error** and the amendment extends it to every
parametric dataset.

## 2. Task 1 — weighting (iv) for the CRC

Score the model against **EC50 = 11.6 ± 0.9 µM** and **nH = 2.2 ± 0.4**, not against the nine
generated points. That means extracting an EC50 and a Hill slope **from the model**, by the same
procedure the source used, and comparing at the published errors. Read both values and both
errors from the project's own records — **no literals** — and say where you read them from.

Three things to handle explicitly rather than quietly:

1. **`FIT_TARGETS['ec50_um']` is not the sourced value.** `data.py` line ~515 says so. Surface
   the discrepancy, report both numbers, and **do not resolve it in either direction** — report
   `kind=question` if the choice would change your answer.
2. **How you extract nH from the model matters.** The published statement covers "between 0.001
   and 0.01 mM GABA", i.e. a **1–10 µM** window, while the reconstruction spans 0.3–3000 µM with
   four decades extrapolated. `models.base.hill_slope` exists and
   `JAHN1997_PEAK_CRC_RESTRICTED` covers the measured window. **Fit the slope over the measured
   window**, state that you did, and report what you would get over the full range as a
   sensitivity check.
3. **Normalisation.** The reconstruction is `X/(1+X)`, which saturates at 1. State how you
   normalise the model's peak response before extracting EC50, since `po_max` is itself an
   observable with a declared band.

Keep the nine-point objective available as a reported view, exactly as weightings (i)–(iii)
were kept. **Do not delete it.**

## 3. Task 2 — re-run the four-way fit, symmetrically

Same four observables, same code path and budget as `scripts/joint_fourway_fit.py`, now with
**every parametric dataset scored against its published parameters**:

| observable | scored against |
|---|---|
| deactivation | 6 published components at published errors (already done) |
| **peak CRC** | **EC50 = 11.6 ± 0.9 µM, nH = 2.2 ± 0.4** (new) |
| paired-pulse recovery | 0.33 ± 0.03, digitised scalar, stated SEM (unchanged) |
| `po_max` | band [0.70, 0.80] (unchanged) |

Schemes: `kinetic_two_open` and `kinetic_jw95_free_rates`, as before.

**Do not invent a relative weighting.** Every quantity now has a published error; use it. If you
find yourself needing a weight, stop and report `kind=question` — that failure produced briefs 5
and 6 and then produced this one.

Answer in `RESULTS-BRIEF-12.md`, with every residual in σ, per scheme:

* Is there a vector satisfying all four within 3σ? **Yes or no.**
* **Does your brief-11 "no" survive the symmetric treatment?** If it does, say so — it is then a
  much stronger result, because it no longer depends on which dataset was corrected. If it does
  not, say that plainly and state that `RESULTS-BRIEF-11.md`'s verdict is withdrawn, updating
  that file in place while preserving what it claimed.
* If irreconcilable, **name the pair and the margin**, in the corrected units.

Report dof honestly as you now do: n_data, k_fitted, residual dof, p-values only where residual
dof > 0. Note that n_data **drops** under this correction — the CRC contributes 2 numbers, not 9
— so residual dof gets worse, and `kinetic_two_open`'s 8 parameters against roughly 9 data is
close to saturation. **Say what that does to the testability of any fit you find.**

## 4. Verify

Both tiers, **final pytest line pasted verbatim**. Baseline you reported: 568 fast / 640 full.
Do not touch or gate the Apple Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] CRC scored against EC50 and nH at published errors; values and errors read from project records, source named, no literals.
- [ ] The `FIT_TARGETS['ec50_um']` discrepancy surfaced with both numbers, unresolved.
- [ ] Hill slope extracted over the **measured** 1–10 µM window, with a full-range sensitivity check reported.
- [ ] Normalisation convention stated.
- [ ] Nine-point objective retained as a reported view.
- [ ] Four-way fit re-run symmetrically, same code path and budget, per-dataset residuals in σ, no invented weighting.
- [ ] `RESULTS-BRIEF-12.md` states whether brief 11's "no" survives; if not, `RESULTS-BRIEF-11.md` updated in place.
- [ ] dof reported honestly, with explicit comment on near-saturation.
- [ ] Both tiers pass; pytest lines verbatim.
- [ ] Commit per task; `kind=blocked` with what landed rather than running to the timeout.

## Not yours to decide

`kind=question`: the `FIT_TARGETS` vs sourced EC50 choice if it changes your answer; adoption of
any scheme; promotion of any finding to settled.
