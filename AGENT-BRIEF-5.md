# Agent brief 5 — the conflict was measured with half the model pinned

**Written by:** Claude Opus 5, 2026-10-16, after verifying brief 4's execution.
**Self-contained: assume you know nothing.** One agent, one commit per task.
**Push at the end of EVERY task.**

```bash
git worktree add ../cp-brief5 -b claude/agent-brief-5 claude/agent-brief-4
cd ../cp-brief5
python3.11 -m venv .venv && source .venv/bin/activate   # 3.13 does NOT work
pip install -e ".[dev]"
pytest -q -m "not slow" -p no:randomly
```

**Branch from `claude/agent-brief-4`, not from the trunk** — this brief builds directly on it.
The rules in `AGENT-BRIEF.md` §Rules apply in full: never skip a test; verify with **bare
`pytest`** too; nothing in `src/` imports from `scripts/`; numbers printed by code, never
typed into prose; **record your own errors in the commit message**; cost-model before you
profile; do not loosen a bound to stop a fit pinning against it.

---

## 0. Brief 4 was good work, and I verified it

580 tests pass. `holdout_guard` is intact (`data.py` and `likelihood.py` untouched, guard
still fires at both fitting entry points), `TRAIN`/`HOLDOUT` are unchanged, and the diagnostic
is a genuinely separate path. The square-table methodology is sound — ratios against each
dataset's *best achievable* RMSE, published SEM for the holdout.

**The central finding reproduces.** I checked it independently rather than accepting the
table: at the parameters that satisfy recovery (R = 0.3695, **+1.32σ**), deactivation RMSE is
**0.3423** against a best achievable **0.0951**. The reverse holds too — fitting deactivation
alone misses recovery by **+16.79σ**. The two Barberis observables, from the same paper,
preparation and cells, are in real tension under the model as fitted. **And the Jahn CRC is
exonerated**, which was not the expected answer.

Two things are nevertheless wrong, and they are the subject of this brief.

---

## 1. Task A — the stated mechanism is wrong. Correct it.

**Size:** S. **Do this first; it is currently printed as a finding.**

`scripts/isolate_dataset_pull.py` prints, and the conclusion rests on:

> *"In the Jones-Westbrook scheme (Model B), desensitization states are accessed from the open
> channel. Slow recovery from desensitization (slow r) necessarily feeds back into the open
> state, producing a long resensitization tail that slows deactivation."*

**That is not this scheme's wiring.** Read `KineticAllosteryModel.q_matrix`. States are
`[R, AR, A2R, A2O, A2D]`, and the two rows that matter are:

| row | transitions |
|---|---|
| `A2R` | `→ AR` at `2·koff`, `→ A2O` at `beta`, **`→ A2D` at `d`** |
| `A2D` | **`→ A2R` at `r`**, and nothing else |

So desensitisation **enters from the doubly-bound CLOSED state and returns to it**. `A2D`
connects only to `A2R`. It is not accessed from the open channel and it does not feed back
through the open state. The classic JW95 side-branch.

### Steps

1. Delete the incorrect mechanism from `scripts/isolate_dataset_pull.py` and from any
   knowledge document that repeats it. **Keep the finding** — the tension is real and
   reproduces — but state it as *observed* rather than explained, until Task B says whether
   there is anything left to explain.
2. Add a test asserting the actual connectivity of `A2D` (that its only non-diagonal entry is
   to `A2R`), so a future claim about this topology has something to fail against. A
   mechanism asserted in prose and contradicted by the Q matrix is exactly what this project
   builds guards for.
3. In the commit, say plainly that the mechanism was wrong and the finding stands.

---

## 2. Task B — redo the consistency verdict with the conventions freed

**Size:** M. **The substance of this brief.**

### Why the verdict is premature

`MODEL_SPECS["kinetic_jw95"]` fits **three** parameters:

```
param_names = ("log10_kd", "log10_E", "log10_D")
conventions = ("koff", "alpha", "r")        # held FIXED, not fitted
```

Those three are the **equilibrium-identifiable ratios**. `koff`, `alpha` and **`r`** are
conventions. And `r` is the resensitisation rate — **the rate that governs recovery from
desensitisation**, which is precisely what the paired-pulse holdout measures. It is pinned at
`0.002` and was never fitted.

So "the three real datasets are NOT mutually consistent" was established while **half the
scheme's rates were held at conventions**, including the one that controls the observable
being missed. That is a statement about a three-parameter sub-model, not about the scheme and
not about the data.

This is also exactly what the project's own ladder predicts.
`knowledge/11-identifiability.md` §1.2 says an equilibrium observable fixes three
combinations and that **absolute rates require kinetic observables**. There are now two
kinetic observables in the repo — the deactivation and the paired-pulse recovery — and
neither is being used to determine an absolute rate. `comparison.py:569` already says so in
as many words: *"the absolute rates are NOT fitted: koff, alpha and the resensitisation rates
are …"*.

### A negative I already checked, so you do not repeat it

Freeing `r` **alone**, holding `D = d/r` at the fitted 25, does **not** resolve it. Predicted
recovery stays at 0.748–0.982 across `r ∈ [2e-4, 5e-2]` — a miss of +13.9σ to +21.7σ
throughout — while deactivation RMSE drifts from 0.174 to 0.222. One direction is not the
test; a joint fit over the enlarged space is.

### Steps

1. Add a model spec that frees the conventions — at minimum `r`, and `koff` and `alpha` too
   if the cost allows. Bounds in log10, as the existing specs do. **Do not delete or change
   the existing three-parameter spec**: it is what every earlier result was computed with, and
   both must coexist so the comparison is like-for-like.
2. Re-run brief 4's square consistency table under the enlarged spec. Same script, same
   denominators, one extra row-set.
3. **Report whether the conflict survives.** Three outcomes, all publishable:
   * **It dissolves** — then the "datasets are inconsistent" verdict was an artifact of the
     conventions, and the real finding is that the paired-pulse recovery *determines `r`*,
     which is the first absolute rate the project has ever had two independent kinetic
     observables for. Say so, and say what `r` comes out as.
   * **It survives** — then the inconsistency is a property of the topology rather than of
     the parameterisation, and Task C becomes worth doing.
   * **It partly dissolves** — report the residuals, do not round the verdict either way.
4. Watch for the failure modes this branch has already hit: parameters pinning at bounds
   (`ModelFit.at_bound` exists — use it), `po_max` leaving `FIT_RANGES` (0.70, 0.80), and
   `k_effective` discounting bound-pinned parameters in any information criterion.
5. **Cost-model first.** A PEAK prediction integrates a waveform per concentration; going
   from 3 to 6 free parameters multiplies the fit cost. Time one evaluation, multiply by
   grid × restarts, then choose.

**Done when** the square table exists for both the three-parameter and the enlarged spec, and
the commit states whether freeing the conventions resolves the conflict.

---

## 3. Task C — only if the conflict survives Task B

**Size:** M. **Conditional. Do not start it if Task B dissolves the conflict.**

If the tension is still there with the rates free, then it is topological, and the question
becomes which topology. The honest next test is an **alternative desensitisation wiring** —
for instance a desensitised state reachable from the **open** state rather than from `A2R`
(which is what brief 4 wrongly believed Model B already was), or a second desensitised state
that does not return through `A2R`.

1. Implement **one** alternative wiring as a new model, alongside the existing ones. Do not
   modify Model B.
2. Fit it to the same datasets with the same diagnostic and print the same square table.
3. Report whether the alternative reconciles the two Barberis observables. **If it does not
   either, that is the stronger result**: two topologies and a free rate set cannot fit two
   observables from the same cells, which would point at the data or at the protocol
   reconstruction rather than at the model.
4. **Do not adopt a new model as the project's model.** This is a test of whether the
   conflict is topological, not a model-selection exercise. `TRAIN`, `HOLDOUT` and the
   existing specs stay as they are.

---

## 4. Do not do these

* **Do not fit against the holdout** or weaken `holdout_guard`. Brief 4's diagnostic path is
  the pattern: a separate, clearly-named entry point whose output is never adopted.
* **Do not change `TRAIN` or `HOLDOUT`.**
* **Do not delete the three-parameter spec** or change what earlier results were computed
  with. Add alongside.
* **Do not loosen a bound** to stop a fit pinning against it, and do not widen
  `FIT_RANGES["po_max"]`.
* **Do not revisit** settled questions: `po_max`'s basis (`FITTED`), `DataKind` subdivision,
  Model C's support (unsupported), the Jahn CRC option (Option 3, and brief 4 exonerated it).
* **Do not widen the split-R̂ gate** in `test_posterior_recovers_known_parameters`.
* **Do not chase the blocked literature** or invent a mean open time for α1β2γ2L.

---

## 5. Summary

| # | task | depends on | size |
|---|---|---|---|
| A | correct the wrong mechanism; keep the finding; guard the topology with a test | — | S |
| B | **redo the consistency verdict with `r` (and ideally `koff`, `alpha`) freed** | A | M |
| C | alternative desensitisation topology — **only if B's conflict survives** | B | M |

Task B is the brief. The conflict may be real, but it has not yet been measured with the
model's own rates free — and the rate that governs the missed observable is currently a
convention.
