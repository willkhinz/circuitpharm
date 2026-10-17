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

## Defect 1 — CORRECTION: the nesting-violation diagnosis was wrong

The original review argued below that the 6-parameter enlarged fits were failed optimisations because their recovery miss (+12.16σ) was worse than the 3-parameter model (+7.14σ). **That inference was invalid.**

### What was corrected (see `scripts/nesting_check.py` and `AGENT-BRIEF-6.md` §1–§2)

1. **Nesting constrains the total objective, not single-dataset residuals.** An enlarged model is free to trade off error between datasets. Under the diagnostic objective:
   * `BARBERIS2007_DEACT` has ~100 points with σ profiled out (`sigma2 = rss/n`). Halving RMSE gains $n \ln 2 \approx 69$ nats.
   * `BARBERIS2007_PAIRED_PULSE` is 1 scalar with fixed SEM = 0.03. Moving from 7σ to 12σ costs only $(144 - 49)/2 \approx 47$ nats.
   * Because 69 > 47, the enlarged model correctly buys trace fit by sacrificing the recovery point. The larger recovery miss is an exchange rate effect under an undeclared weighting, not an optimisation failure.

2. **Nesting holds exactly.** Measured in a clean process via `scripts/nesting_check.py` at `(log10_kd, log10_E, log10_D) = (-0.428, 0.104, 1.9721)` on Deact + Recovery:
   * `kinetic_jw95` (3 params): total log-likelihood = **-2.435202**
   * embedded in `kinetic_jw95_free_rates` (6 params): total log-likelihood = **-2.435225**
   * Agreement: $2.29 \times 10^{-5}$ nats; rate constants agree to $10^{-16}$; Q matrix max abs difference is $2.84 \times 10^{-14}$.

3. **Correction to convention values.** The original review erroneously quoted `log10_koff = -0.477` and `log10_alpha = -0.967`. The actual convention values in `FIT_CONVENTIONS` are `log10_koff = -0.743727` (`koff = 0.180415`) and `log10_alpha = -0.522879` (`alpha = 0.30`).

### Original text (preserved for history)

> Brief 5 concludes that neither topology nor any freeing of the rate conventions can reconcile
> the two Barberis observables. **The numbers show a nested model fitting worse than the model
> it contains, which is an optimisation failure, not a result.**
>
> | fit | Deact + Recovery, recovery miss |
> |---|---|
> | `kinetic_jw95`, 3 params (brief 4) | **+7.14σ** |
> | `kinetic_jw95_free_rates`, 6 params (brief 5) | **+12.16σ** |
> | `kinetic_open_desens`, 3 params (brief 5) | **+7.56σ** |
> | open-desens, 6 free rates (brief 5) | **+14.00σ** |
>
> The 6-parameter space **strictly contains** the 3-parameter one — verified: every convention
> value lies inside the free bounds (`log10_r = −2.699` in `(−5, 1)`; `log10_koff = −0.477` in
> `(−3, 2)`; `log10_alpha = −0.967` in `(−2, 2)`). A model that contains another **cannot** fit
> worse at its optimum. Both 6-parameter rows are therefore failed optimisations reported as
> findings, and the same pattern appears under both topologies.

## Defect 2 — the test claim is mislabelled

The report states *"Full test suite passing (519 passed in .venv, 0 failures)."* **519 is the
fast subset.** Measured here: `pytest -q -m "not slow"` gives exactly 519; brief 4's full
suite was 580. The slow tests were never run before pushing, contrary to the standing rule.
A full run was started during this review but did not finish before the limit.

## Status of verification

* **The full suite on `claude/agent-brief-5` is now verified: 585 passed in 689.94 s.** (Closes the gap left open in this review.)
* The Task C open-desensitisation model's own correctness beyond its connectivity test.
* Whether the conflict survives under declared alternative weightings — addressed in Brief 6.
