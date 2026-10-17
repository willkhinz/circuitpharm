# Agent brief 9C — make the canonical route affordable, then answer the question

**Continuing `AGENT-BRIEF-9.md` and `AGENT-BRIEF-9B.md`, same worktree, same branch
(`claude/agent-brief-9`).** Both still apply except where this file overrides them. **I stopped
your last run after 21 minutes** — not because of anything you did wrong. You were committing
per task as instructed (`93c990c`, `ad2ef0b`, both kept), and the anti-polling change worked.
I stopped it because the task as written was infeasible and I have a fix.

## 1. The measurement that changes the plan

I ran your `scripts/deact_obstruction_decomposition.py` at a deliberately tiny budget
(`--n-samples-fit 12 --n-polish-fit 2`):

| route | all six | τ only | amplitudes only | wall time |
|---|---|---|---|---|
| eigen, 150 samples / 3 polish | 528.07 | **0.00** | 453.01 | 0.46 s total |
| **fit (canonical), 12 samples / 2 polish** | 1679.46 | **205.68** | 1473.79 | **1540 s total** |

Two things follow.

**(a) The cost is prohibitive as written.** About **37 s per fit-route objective evaluation**
against roughly **0.001 s** for the eigen route — four orders of magnitude. At the brief's
default of 100 samples that is ~67 minutes per cell, so two schemes × three objectives is near
seven hours. You would have hit the timeout again, and that would have been my fault, not
yours.

**(b) The headline result is in question, and the budget is why we cannot yet tell.** Under the
canonical route **τ-only is 205.68, not 0.00**. If that survives a proper search, then "the
time constants fit exactly and the obstruction is purely in the amplitudes" is **false under
the canonical route**, and `FINDING-amplitude-obstruction.md` falls. But 12 Latin-hypercube
samples in a six-dimensional box is far too few to distinguish a real result from a failed
optimisation. **Do not report either way from a budget like that.**

## 2. The fix: warm-start the triexponential fit from the eigen components

The fit route is slow because every objective evaluation multi-starts a triexponential fit from
scratch. But we already have excellent starting values: the eigen decomposition gives the
model's true dominant components for that very waveform, at ~0.001 s.

So in `fit_triexponential_decay`, accept explicit start values, and in `predict_components_fit`
default them to the three largest-|amplitude| eigen components from
`predict_components_eigen_all`. A fit started essentially at its answer should converge in a
few iterations.

**This must not become a way of smuggling the eigen answer into the canonical route**, so guard
it:

* Keep **random/multi-start fits as well** as the warm start — fewer of them is fine.
* **Assert the warm-started optimum is no worse than the best random-start optimum** on the
  fit's own cost. If a random start ever beats the warm start, the warm start has biased the
  answer: report that, with both numbers, and treat the better one as the fit.
* Report, per cell, how often the warm start won and the spread of the top optima, the way
  `fit_mle_vector` does.
* Where the fit is genuinely degenerate — one eigen mode carrying ≥0.90 — it will stay
  degenerate with a warm start too. **That is still reported as degenerate, never scored.**

Measure and report the speed-up. If it is not at least ~20×, say so, and fall back to §4.

## 3. Then answer the question properly

With the warm start in place, run the decomposition at a budget you can defend — state it —
for `kinetic_jw95_free_rates` and `kinetic_open_desens_free_rates` first, both routes, all
three objectives. Then answer in `RESULTS-BRIEF-9.md`:

> Under the canonical fit route, at an adequate search budget, is the obstruction still in the
> amplitudes rather than the time constants?

and state explicitly whether `FINDING-amplitude-obstruction.md` **stands or is withdrawn**. If
τ-only remains of order 200 rather than 0, it is withdrawn — **say so plainly and update that
file in place, preserving what it originally claimed.** A withdrawal is a result, and this
project has already recorded two of mine.

Keep 9B's conditioning requirement in force: raw `Σaⱼ` and `P_o(0)` at every optimum,
ill-conditioned cells rejected and labelled with the stated threshold, and an explicit statement
of whether the amplitudes-only figure survives it.

## 4. If the warm start does not rescue it

Fall back to a declared **two-stage** procedure and label it as such: **search** with the cheap
eigen objective, then **score** the resulting optimum canonically with the fit route. That is
legitimate — searching cheaply and scoring canonically is normal practice — but it is a weaker
claim, because the search may never visit the region the canonical objective prefers. **If you
use it, say so in the results file in those terms.**

## 5. Discipline reminders that killed previous runs

* **Foreground, and commit per task.** Both your last-run commits survived because of this.
* **Report `kind=blocked` with what landed** rather than running to the timeout. That is a
  correct outcome.
* Do not resurrect either of my withdrawn mechanisms for the obstruction — see
  `FINDING-amplitude-obstruction.md`. The measurements stand; the "why" is open and is not your
  task.

## Acceptance criteria

- [ ] `fit_triexponential_decay` accepts start values; `predict_components_fit` warm-starts from the eigen components by default.
- [ ] Random starts retained; warm-start-no-worse assertion in place and its outcome reported per cell.
- [ ] Measured speed-up reported.
- [ ] Decomposition for `kinetic_jw95_free_rates` and `kinetic_open_desens_free_rates`, both routes, three objectives, at a stated and defensible budget.
- [ ] Raw `Σaⱼ` and `P_o(0)` at every optimum; ill-conditioned cells rejected with the threshold stated.
- [ ] Degenerate cells labelled, never scored or backfilled from the eigen route.
- [ ] `RESULTS-BRIEF-9.md` answers §3's question plainly and states whether the finding stands or is withdrawn; if withdrawn, that file is updated in place.
- [ ] Both tiers run with `PYTHONPATH` set; counts stated, labelled, accounted for against 546 fast / 618 full.

## Not yours to decide

`kind=question`: adoption of any new scheme; promotion of the finding from provisional to
settled.
