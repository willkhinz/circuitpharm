# Agent brief 7B — predicted deactivation components by fitting the simulated waveform

**Read `/Users/whinz/Biochemistry/BRIEF-7-SHARED.md` first, in full.** It holds the discipline
rules, the published measurement table, why weighting (iv) exists, and the shared acceptance
criteria. This file only adds what is specific to your method.

**Your worktree:** you are already in it, on branch `claude/agent-brief-7b`.
**Another agent (7A) is solving the same problem by a different method, in a different
worktree.** Do not look for its work, do not coordinate, do not try to match it.

---

## Your method: do what the experimenters did

Barberis 2007 did not report eigenvalues. They applied a 2 ms pulse of 10 mM GABA, recorded
the current, and **fitted a three-exponential to the decay**. Your job is to reproduce that
procedure against the model's own simulated current, so that what is compared to the published
numbers is the same kind of quantity the authors produced.

1. **Simulate the protocol.** 2 ms pulse of 10 mM GABA from rest, then washout, over the
   dataset's own fit window (`_BARBERIS_FIT_WINDOW_MS = 1000.0` in
   `src/circuitpharm/fitting/data.py`). Use the existing waveform machinery rather than
   writing a new integrator; say which function you used.
2. **Fit a triexponential** `Σ aᵢ exp(−t/τᵢ)` to the decay, normalised the way
   `BARBERIS2007_DEACT.normalized_current` is normalised — confirm that before choosing, do
   not guess.
3. **Report (τ₁, τ₂, τ₃) and (a₁, a₂, a₃)** ordered fastest-first, matching the published
   table's convention.

**The sampling grid is part of the measurement, and this has already bitten this project
once.** `data.py` records that a uniform 10.1 ms grid left τ_fast unidentified — 2.60 ± 2.43
ms, 98% scatter — against 2.81 ± 0.18 ms on the hybrid log/linear grid, measured over 400
replicates by `scripts/deact_sampling.py`. **Fit on the dataset's own `_T_BARBERIS` grid**,
and state what you would get on a uniform grid as a sensitivity check. If your fitted τ_fast
is not stable across both, your predictor is reporting the grid and not the model.

**Your fit is itself an optimisation and must be shown to have converged.** Multi-start it,
report the spread of the top optima the way `fit_mle_vector` does, and refuse to report a
single winner from a tie. A three-exponential fit to a four-exponential truth is
ill-conditioned by construction — expect trouble and report it rather than smoothing it over.

## Specific deliverables

1. `src/circuitpharm/fitting/deact_components.py` (new) with
   `predict_components_fit(model, *, pulse_ms=2.0, conc_mm=10.0) -> tuple[taus, areas]`.
   Pure function, no global state. Document the grid, the normalisation and the start strategy.
2. A likelihood term — your call where, and say where — scoring predicted components against
   `circuitpharm.literature.DEACTIVATIONS[1]`'s published values and errors. Read the errors
   from that record and from its `note` field, which carries the component and area errors;
   **do not retype them as literals.** If they exist only as prose in `note`, parse them or
   promote them to structured fields — and if you promote them, that is a change to a
   literature record, so state it prominently and keep the prose intact.
3. A test that your predictor **recovers the published components from the published curve**:
   run it against `BARBERIS2007_DEACT.normalized_current` itself, which was generated from
   those six numbers, and assert it returns them to a stated tolerance. This is the strongest
   check available to your method and it is cheap — if it fails, your predictor is wrong and
   nothing downstream means anything.
4. `RESULTS-BRIEF-7B.md`, per the shared criteria, additionally reporting: your components at
   each of the three models' brief-6 optima, the top-optima spread for each fit, the
   uniform-grid sensitivity check, and whether the triexponential fit was well-conditioned at
   each optimum.

## What is yours to decide and what is not

Yours: the fit's start strategy and tolerances, where the likelihood term lives, how you
report conditioning.

**Not yours:** whether weighting (iv) becomes the project default. Report `kind=question`
rather than adopting it.
