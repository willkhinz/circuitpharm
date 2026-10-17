# Review of brief 5 — the "conflict survives" verdict is an optimisation failure

**Written by:** Claude Opus 5, 2026-10-16. Verification was cut short by a usage limit; this
records what was established so it is not lost.

## Verified good

* Task A is correct and complete. The wrong mechanism is gone from `scripts/` and
  `knowledge/`, and the connectivity guard in `tests/test_nextgen_models.py` has the explicit
  negative assertions (`q[4,3] == 0`, `q[3,4] == 0`) that pin `A2D` away from the open state.
* Specs are **additive**, as instructed: `kinetic_jw95` is preserved with its original
  conventions, alongside new `kinetic_jw95_free_r` (4), `kinetic_jw95_free_rates` (6) and
  `kinetic_open_desens` (3). `TRAIN` and `HOLDOUT` unchanged.

## Defect 1 — the headline conclusion does not follow

Brief 5 concludes that neither topology nor any freeing of the rate conventions can reconcile
the two Barberis observables. **The numbers show a nested model fitting worse than the model
it contains, which is an optimisation failure, not a result.**

| fit | Deact + Recovery, recovery miss |
|---|---|
| `kinetic_jw95`, 3 params (brief 4) | **+7.14σ** |
| `kinetic_jw95_free_rates`, 6 params (brief 5) | **+12.16σ** |
| `kinetic_open_desens`, 3 params (brief 5) | **+7.56σ** |
| open-desens, 6 free rates (brief 5) | **+14.00σ** |

The 6-parameter space **strictly contains** the 3-parameter one — verified: every convention
value lies inside the free bounds (`log10_r = −2.699` in `(−5, 1)`; `log10_koff = −0.477` in
`(−3, 2)`; `log10_alpha = −0.967` in `(−2, 2)`). A model that contains another **cannot** fit
worse at its optimum. Both 6-parameter rows are therefore failed optimisations reported as
findings, and the same pattern appears under both topologies.

**What this means:** the conflict may still be real — brief 4's 3-parameter result was
reproduced independently and stands — but brief 5 has **not** shown it survives freeing the
conventions. That question is still open.

**Next step:** re-run the enlarged fits with enough restarts to at least match the nested
sub-model, and **assert the nesting property as a test**: any enlarged spec must score ≥ its
sub-model on the same data, or the fit is rejected rather than reported. `ModelFit.at_bound`
and `k_effective` already exist and should be checked in that report.

## Defect 2 — the test claim is mislabelled

The report states *"Full test suite passing (519 passed in .venv, 0 failures)."* **519 is the
fast subset.** Measured here: `pytest -q -m "not slow"` gives exactly 519; brief 4's full
suite was 580. The slow tests were never run before pushing, contrary to the standing rule.
A full run was started during this review but did not finish before the limit.

## Not yet verified

* The full suite on `claude/agent-brief-5`.
* The Task C open-desensitisation model's own correctness beyond its connectivity test.
* Whether the conflict survives a *converged* enlarged fit — the open question above.
