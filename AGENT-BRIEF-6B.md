# Agent brief 6B — finish brief 6 inside a compute budget

**For:** the same Gemini worker, same thread, continuing brief 6.
**Worktree:** `/Users/whinz/Biochemistry-wt-agent-brief-6`, branch `claude/agent-brief-6`.
**You are already in it.** It was created for you. Do not create a worktree, do not switch
branch, do not touch `/Users/whinz/Biochemistry` — the planner reads from there.
**Written by:** Claude Opus 5, 2026-10-16 21:00.

---

## 0. What happened, and why it was my fault not yours

I cancelled your brief-6 run at 20:56. **Your work is safe** — I committed it as `dadcda2`
(WIP, clearly labelled as yours and as unverified).

`AGENT-BRIEF-6.md` Task C asked for three weightings × three models, each with **≥24
Latin-hypercube starts** *and* tolerances tightened to `ftol=1e-14` / `xatol=1e-9` over stiff
ODE solves. That is roughly **nine times** brief 5's diagnostic cost. Measured from your run:
the cheapest of the nine fits (3-parameter, as-is) took **~25 minutes**, and the 6-parameter
fit had run **27.5 minutes** without finishing. One of nine fits was cached after 55 minutes.
The run would have hit its 4-hour timeout with **no comparison to show** and tasks D and E
unstarted.

**That was my sizing error.** You implemented what I specified, including the per-fit caching
I asked for, and the caching is why nothing was lost. This brief replaces the compute plan
with one that fits in a budget. Everything else about brief 6 still applies — **read
`AGENT-BRIEF-6.md` §0 (discipline), §1 and §2 in full**; they are still the context, and §2 is
still the question.

### Already done — do not redo

* **Task A** — `20fcd78`, `REVIEW-BRIEF-5.md` corrected. Complete.
* **Task B** — `89d070c`, `ModelSpec.nests`, `nesting_embedding`, `extra_starts`,
  `tests/test_nesting.py`. Complete. **You will now use it** — see §1.
* Your in-flight Task C is `dadcda2`: `scripts/weighting_sensitivity.py` plus edits to
  `diagnostic.py` and `dataset_consistency.py`. **I have not verified any of it.**

---

## 1. The re-scoped compute plan

**Stop brute-forcing the enlarged fits. Use the nesting machinery you just built.**

24 blind Latin-hypercube starts in a 6-dimensional box is the expensive way to find a basin
you already know the location of. `nesting_embedding` maps the 3-parameter optimum into
6-parameter coordinates, and from that point the enlarged fit **cannot go backwards** — that
is exactly what `tests/test_nesting.py` now asserts. So:

1. **Fit the 3-parameter spec first**, per weighting. Modest starts are fine here: the box is
   3-dimensional.
2. **Warm-start every enlarged fit** from the embedded 3-parameter optimum, plus a handful of
   random starts for basin coverage — **6**, not 24.
3. **Relax the tolerances back** to the diagnostic's own values: `ftol=1e-8, gtol=1e-6,
   maxiter=500`, Nelder–Mead `xatol=1e-4, fatol=1e-4, maxiter=200`. The tightened values were
   my instruction and they are not what is protecting these fits — **the nesting assertion
   is.** Say so in a comment where you set them.

**Hard budget, and it is a stop condition, not a target:**

* **≤ 5 minutes per fit.** If a single fit exceeds that, stop it, record the start count and
  elapsed time, and move on with fewer starts for that spec. A fit you had to truncate is
  reported as truncated.
* **≤ 45 minutes for the whole sweep.** If you reach 45 minutes, stop, write up what you
  have, and report what is missing. **An incomplete comparison reported honestly is worth
  more than a complete one nobody is awake to read.**
* Keep the per-fit cache. **Make every source edit before you launch the sweep** — your cache
  key hashes `diagnostic.py`, `comparison.py` and `likelihood.py`, so editing any of them
  after caching silently invalidates every entry. Your one cached fit will be invalidated by
  the start-count change anyway; that is expected, not a bug.
* `--workers 8` on a 10-core machine drove the load average to 28. **Use 6.** Say what you
  used.

**Resolution is being traded for coverage, deliberately.** The deliverable is the
*comparison across three weightings*, not a high-precision optimum in any one of them. If a
number is too coarse to support a conclusion, **say it is too coarse** rather than quoting it
with false confidence.

---

## 2. Tasks

### Task C1 — verify your own WIP before building on it

`dadcda2` is unverified. Before running anything, check and state:

1. That each of the three weightings in `diagnostic_log_likelihood` does what it claims:
   * **(i) as-is** — deactivation σ profiled out (`sigma2 = rss/n`), recovery at stated SEM.
   * **(ii) deactivation σ fixed** at the single-dataset baseline RMSE `0.095091`.
   * **(iii) equal weight per dataset** — each dataset's contribution divided by its `n`.
2. That **(i) reproduces the existing behaviour exactly.** Evaluate one model under the
   as-is weighting and against `diagnostic_log_likelihood` as it was before `dadcda2`
   (`git show 89d070c:src/circuitpharm/fitting/diagnostic.py`). They must agree to numerical
   tolerance. **If they do not, weighting (i) is not the status quo and every comparison
   against it is meaningless.** Fix it before going further.
3. That the cache key covers everything that changes a fit: spec, weighting, fixed σ, start
   count, seed, and the source hash. Name anything it misses.

### Task C2 — run the sweep and report

Three weightings × three models (`kinetic_jw95`, `kinetic_jw95_free_rates`,
`kinetic_open_desens`) on Deact + Recovery, under §1's budget.

Write `RESULTS-BRIEF-6.md` at the worktree root, with the exact command, start counts,
tolerances, worker count and runtime at the top, and then:

1. **The measured nat decomposition** at each fit's optimum under the as-is weighting:
   deactivation term, recovery term and `po_max` anchor term separately, and how many nats
   the enlarged model gains on deactivation versus pays on recovery. My ≈69-vs-≈47 estimate
   in `AGENT-BRIEF-6.md` §2 was arithmetic from `n` and the SEM — **measure it, do not quote
   mine.**
2. **The recovery σ miss under each of the three weightings**, for each model, in one table.
3. **Is there any parameter vector, under any of the three, that satisfies deactivation
   within 3× its baseline RMSE *and* recovery within 3σ simultaneously?**
4. **The plain verdict.** If the conflict survives all three, brief 5's architectural
   conclusion is strengthened and should be restated **with the weighting declared**. If it
   dissolves under (ii) or (iii), brief 5's conclusion is **withdrawn** and replaced by a
   statement about the weighting. Both outcomes are acceptable. Neither is expected.
5. Label each weighting as **an assumption**. Do not pick a favourite and present its number
   as the answer — that choice is the user's (§4).

### Task C3 — the E12 fix and the guard

Carry these over from `AGENT-BRIEF-6.md` Task C steps 4–6; they are cheap and independent of
the sweep:

* The warm starts hardcode **`alpha` = 0.36** (`log10` = −0.443697) against
  `FIT_CONVENTIONS["alpha"][0]` = **0.30** (`log10` = −0.522879), and `koff` as `log10(0.18)`
  = −0.744727 against −0.743727. Sites: `scripts/dataset_consistency.py` `starts_6` (both
  entries) and `src/circuitpharm/fitting/diagnostic.py` `fit_diagnostic_consistency`, the
  `elif len(names) == 6:` branch. Replace all of them with `nesting_embedding`.
* **No convention literal may survive** in `src/` or `scripts/`. Prove it: paste the output of
  `grep -rn "0\.36\|log10(0\.18)\|-0\.44\|-0\.74" src/ scripts/`, empty or every hit
  explained.
* `scripts/dataset_consistency.py` must **raise**, not print, if a nested fit scores worse in
  **total** log-likelihood than its embedded sub-model.

### Task D — rewrite the three tests that hardcode fitted vectors

Unchanged from `AGENT-BRIEF-6.md` Task D. In `tests/test_dataset_consistency.py`:
`test_conflict_survives_freeing_rate_conventions`,
`test_conflict_survives_under_alternative_open_desensitisation_topology`, and
`test_datasets_are_not_mutually_consistent`. Each must fit, then assert a property of the fit
it just computed. No hardcoded fitted vectors. **Do not delete, skip or `xfail` any of them.**
If the sweep shows the conflict does not survive some weighting, rename each test to what the
evidence supports and assert that, and say in the commit message that the previous assertion
was reached under an undeclared weighting.

### Task E — verify

```bash
cd /Users/whinz/Biochemistry-wt-agent-brief-6
/Users/whinz/Biochemistry/.venv/bin/python -m pytest -q -m "not slow" -p no:randomly
/Users/whinz/Biochemistry/.venv/bin/python -m pytest -q -p no:randomly          # full, ~11.5 min
```

Baselines: **519** fast, **585** full, plus what Task B and this brief add. **State both
counts separately and label which is which.**

**Known platform difference — do not touch:** a marginal split-R̂ failure in
`test_posterior_recovers_known_parameters` on Apple Silicon is a known platform artefact. **Do
not widen the gate**, do not mark it. Report it if you see it.

Then prove the import path still works from outside the repo:

```bash
cd /tmp && /Users/whinz/Biochemistry/.venv/bin/python -c "import circuitpharm; print(circuitpharm.__file__)"
```

---

## 3. Acceptance criteria — `kind=done` only when every one is met

- [ ] Task C1's three checks stated, with weighting (i) shown to reproduce pre-`dadcda2` behaviour to numerical tolerance.
- [ ] Enlarged fits warm-started via `nesting_embedding`; ≤6 random starts; tolerances back to the diagnostic's own values with a comment saying why.
- [ ] Every source edit made **before** the sweep launched.
- [ ] `RESULTS-BRIEF-6.md` contains all five items in Task C2, the command, start counts, tolerances, worker count and runtime, and a plain verdict.
- [ ] Any fit or sweep truncated by the budget is **reported as truncated**, with its elapsed time and start count.
- [ ] No convention literal in `src/` or `scripts/`; `grep` output pasted.
- [ ] `scripts/dataset_consistency.py` raises rather than printing a table whose nested fit scores worse on total log-likelihood.
- [ ] All three tests in Task D rewritten with no hardcoded fitted vectors; none skipped, deleted or `xfail`ed.
- [ ] `tests/test_nesting.py` and `tests/test_fit_determinism.py` both still pass.
- [ ] Both pytest tiers pass, **counts stated separately and labelled**.
- [ ] Every number was printed by a committed script, and you name the script.
- [ ] Commit messages record any error you made.

## 4. Out of scope, and what is the user's call

Out of scope: `TRAIN`/`HOLDOUT`, `holdout_guard`, the datasets, the sampling grids,
`BARBERIS2007_PAIRED_PULSE.sem`, `provenance.py`, the `sources` table, the P6 σ requirement,
Model C's adoption. The alternative weightings are **additional reported views, not edits to
the data.**

Report `kind=question` rather than deciding, if you reach either:

1. **Which weighting becomes the project's declared default.** Report all three; the user
   chooses. Silently adopting one is the exact failure this brief exists to fix.
2. Any need for more than the 45-minute sweep budget. Report what another hour would buy.

**Separately noted, not your task:** your Task B measurement that `extended_desens` does not
nest `kinetic_jw95` — 16.66 nats on Deact+Recovery, because `r_fast` = 0.005 against
`r` = 0.002 — contradicts that spec's own note claiming "NESTS MODEL B … so C cannot fit worse
than B at its own optimum". That is a good catch and it reaches past brief 6. It goes in
brief 7. Do not act on it here.

## 5. If you get blocked

`kind=blocked` with the exact command, the exact output, and what you have ruled out.
