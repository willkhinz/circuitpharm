# Brief 7 — shared context for agents 7A and 7B

**Read `AGENT-BRIEF-6.md` §0 in full first.** Its discipline rules apply unchanged: a failed
optimisation is not a result; never fabricate a citation or a number; numbers are generated,
not transcribed; never skip, delete or `xfail` a test to get green; record your own errors in
the commit message. Baselines on `main`: **524** fast (`pytest -q -m "not slow"`), **596**
full. State which tier you ran.

**You are already in the git worktree you must work in.** It was created for you. Do not
create another, do not switch branch, do not touch `/Users/whinz/Biochemistry` — the planner
reads from there.

**Run tests with `PYTHONPATH=<your worktree>/src`.** The shared `.venv` has one editable
install pinned to the main checkout, which sits on `main`. Without the env var, pytest
collects your tests but imports `circuitpharm` from `main`, and you will see ~15 failures and
~15 collection errors that are not yours. Measured; do not spend time on it.

Interpreter: `/Users/whinz/Biochemistry/.venv/bin/python`.

---

## Why this brief exists

Brief 6 asked whether the Barberis deactivation / paired-pulse conflict was real or an
artefact of how the objective weighs a 100-point trace against a 1-point scalar. It measured
three weightings. Under the as-is weighting the datasets look irreconcilable; under two
alternatives they reconcile. The verdict therefore depended entirely on a weighting
convention nobody had chosen — recorded in `RESULTS-BRIEF-6.md`.

**The decision taken (Claude Opus 5, 2026-10-16): none of the three is adopted.** All three
are conventions for weighting *pseudo-observations*, and the measured leverage each dataset
carries, as `n/sigma^2`, shows why:

| weighting | deactivation | recovery | ratio |
|---|---|---|---|
| (i) as-is, sigma profiled | 4764.3 | 1111.1 | **10.05× too much trace** (at the 3-param optimum) |
| (ii) sigma fixed at 0.095091 | 11169.7 | 1111.1 | **10.05× too much trace** |
| (iii) equal weight per dataset, 1/n | 110.6 | 1111.1 | **10.0× too much scalar** |

(ii) and (iii) are symmetric overcorrections, an order of magnitude off in opposite
directions. Choosing between them is choosing which answer you get.

**The reason none of them is right is already written in the repo.** The docstring of
`BARBERIS2007_DEACT` in `src/circuitpharm/fitting/data.py` says:

> Declared `kind="parametric"`, `synthetic=True`: **the parameters are attributable
> measurements; the individual curve points are the paper's model of its data.**

The 101 curve points are not 101 measurements. They are a triexponential evaluated on a grid,
generated from six published numbers. Scoring them as independent observations invents ~100
pseudo-data from 6 real ones, and then the weighting argues about how loudly the invention
should speak.

**The actual measurements, with their published errors** (Barberis 2007, Eur J Neurosci
25(9):2726-2740, PMC1950087, Table 1; rat α1β2γ2S, 2 ms pulse of 10 mM GABA, 22–24 °C,
−70 mV, ultrafast exchange; γ2S substituted for γ2L, recorded in
`circuitpharm.literature.DEACTIVATIONS[1]`):

| quantity | value | published error | n |
|---|---|---|---|
| τ₁ | 2.8 ms | ± 0.3 | 6 |
| τ₂ | 33.4 ms | ± 4.6 | 6 |
| τ₃ | 221.35 ms | ± 14.9 | 6 |
| area₁ | 0.57 | ± 0.04 | 6 |
| area₂ | 0.23 | ± 0.02 | 6 |
| area₃ | 0.20 | ± 0.03 | 6 |
| τ_w (derived) | 52.5 ms | ± 2.9 (SEM) | 6 |
| paired-pulse recovery at 100 ms | 0.33 | ± 0.03 | 8 |

**Weighting (iv), which this brief builds:** score the model against those published
parameters at their published errors, instead of against the reconstructed curve. Then every
datum is a real measurement with a stated error bar, `n` is honest, and **there is no
weighting convention left to choose** — the relative leverage is set by the errors the authors
published. That is the objective this project's discipline actually implies.

## The question (iv) is built to answer

> Scored against the published deactivation parameters at their published errors, rather than
> against a reconstructed curve, does the Barberis deactivation / paired-pulse recovery
> conflict exist at all?

Both answers are acceptable. Neither is expected.

## Why there are two of you

Predicting published time constants from a kinetic model can be done two ways, and **they are
not obviously equivalent**. Agent **7A** does it analytically from the Q matrix eigenvalues;
agent **7B** does it by fitting a triexponential to the simulated waveform, as the
experimenters did. You are each doing it one way, in your own worktree, without coordinating.

Afterwards your predictions are compared at the same parameter vector. **If they disagree,
that disagreement is the finding** — it would mean the paper's own three-exponential fit is
lossy with respect to the model's true relaxation spectrum, which bears directly on whether
the published parameters are the right target. Do not try to make your answer match the other
agent's. Report what your method gives.

A known subtlety you will meet: the JW95 scheme has **5 states, so 4 non-zero relaxation
eigenvalues**, while the paper reports **3** exponential components. Do not paper over the
mismatch. State explicitly how you handle it and what is lost.

## Shared acceptance criteria (both agents)

- [ ] A function that, given a `ReceptorModel`, returns predicted (τ, area) components for the Barberis protocol, by your assigned method.
- [ ] A likelihood term scoring those predictions against the table above at the **published** errors — not at a chosen sigma, and not at the reconstructed curve.
- [ ] A test that your predictor reproduces the published components **when given a model built to have them**, or, if that is not constructible, a stated reason and a weaker check that is honest about what it proves.
- [ ] `RESULTS-BRIEF-7<your letter>.md` with the command, runtime, your predicted components at each of the three models' brief-6 optima, and a plain answer to the question above.
- [ ] Both test tiers pass with `PYTHONPATH` set; counts stated separately and labelled.
- [ ] Do not touch `TRAIN`, `HOLDOUT`, `holdout_guard`, the datasets, the sampling grids, `provenance.py`, or the `sources` table. Weighting (iv) is an **additional** objective; leave (i)–(iii) in place.
- [ ] Do not delete or modify `scripts/weighting_sensitivity.py`'s existing behaviour.

## Report

`kind=done` only when every criterion is met. `kind=question` if you need a decision — in
particular, **do not** decide by yourself that weighting (iv) replaces the project default;
report your result and let the planner and user decide.
