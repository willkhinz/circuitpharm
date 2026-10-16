# Agent brief 3 — five tasks, in order

**Written by:** Claude Opus 5, 2026-10-16, after verifying `AGENT-BRIEF-2.md`'s execution.
**Self-contained: assume you know nothing.** One agent, one worktree, one commit per task.

```bash
git worktree add ../cp-brief3 -b claude/agent-brief-3 claude/quirky-sagan-edtdet
cd ../cp-brief3
python3.11 -m venv .venv && source .venv/bin/activate   # 3.13 does NOT work
pip install -e ".[dev]"
pytest -q -m "not slow" -p no:randomly      # should be green before you start
```

The rules in `AGENT-BRIEF.md` §Rules apply in full — **re-read them.** Especially: never skip
a test; verify with **bare `pytest`** too; nothing in `src/` imports from `scripts/`; numbers
are printed by code, never typed into prose; absence of a declaration is an error, not a
default; **record your own errors in the commit message**; cost-model before you profile.

---

## 0. What brief 2 achieved, verified

All five tasks landed (`8d6c15b`, `08f3799`, `fed8f03`, `9f5b160`, `1641147`). Verified
independently rather than taken from the commit messages:

**Task A — the SEM fix works, exactly.** `noise_model` is now a declared per-dataset field
(`"stated_sem"` for the holdout, `"estimated"` for the four training sets), and `DatasetFit`
reports `noise_model`, `stated_sem` and `sigma_hat` side by side. An 18.6σ miss now scores
**−171.202**, matching the Gaussian with σ = 0.03 to three decimals. The degenerate **+0.577**
is gone. The fitted-model holdout score is **−7.663**, matching the −7.66 predicted in brief 2.

**Task B** — E27 is in `WORKLOG.md`. **Task E** — σ ≤ 0.005 survives, recorded side by side
with the old numbers.

**Three findings worth more than the tasks that produced them**, all in §1–§3 below.

---

## 1. Task A — ~~push~~ DONE during review; check CI instead

**Size:** XS.

Brief 2's five commits were **local only** when this brief was written
(`rev-list` returned `0 5`) — the SEM fix that every quantitative claim now rests on existed
in one working copy, unbacked-up and unseen by CI. **Pushing this brief carried them**, so
they are now on `origin/claude/quirky-sagan-edtdet` (`b81b0e9` and its six parents).

What remains of this task:

1. **Check CI on that push.** The SEM change alters scoring, so expect movement in any test
   asserting a log-likelihood. 504 fast tests pass locally; CI also runs a lean 3.12 matrix
   this session has not exercised.
2. **Push at the end of every task from now on**, not at the end of the brief. Five commits
   of unbacked-up work is the failure this records.

---

## 2. Task B — the `po_max` anchor is missing from the comparison path

**Size:** M. **This is an E12-class defect and it CONFOUNDS Task D of brief 2.**

### The defect

`identifiability.py:157` anchors absolute open probability:

```python
chi2_pomax = float(((model.po_max - 0.75) / 0.05) ** 2)
```

**`comparison.py` has no such term.** And it cannot get the information from the data,
because **both fitted datasets are `normalisation="fraction_of_max"`** — `jahn1997_peak_crc`
and `barberis2007_deact`. The codebase already states the consequence, in `_predict`:

> *"a normalised dataset then constrains only the SHAPE of the curve … Absolute open
> probability has to come from somewhere else."*

In `compare_models` / `fit_one` it comes from nowhere. So `log10_E`, which sets `po_max`,
wanders. Brief 2's Task D fitted `MODEL_SPECS["kinetic_jw95"]`, whose free parameters are
`(log10_kd, log10_E, log10_D)` with `po_max` **derived**, and produced:

| Task D option | `log10_D` | derived `po_max` | PPR prediction |
|---|---|---|---|
| 1, keep 9 points, + deactivation | **−3.0000 = its lower bound** | **0.0100** | 0.9967 (+22.2σ) |
| 2, restricted, + deactivation | −2.9730 (interior by 0.027) | 0.5960 | 1.0000 (+22.3σ) |
| 1, CRC alone, no deactivation | **−3.0000 = its lower bound** | 0.9939 | 0.1942 (−4.5σ) |

`FIT_RANGES["po_max"]` is **(0.70, 0.80)**. Two of three fits are outside it, one by 70×.
Rows 1 and 3 are the **same dataset** hitting the **same bound** with `po_max` differing
**100-fold** — which is only possible if nothing constrains `E`.

### Why this matters more than it looks

**Brief 2's Task D cannot support a recommendation**, and its headline — that both Jahn
options miss the holdout by ~+22σ — is a consequence of `po_max` collapsing, not of the
option choice. Note the tell: the one row with a plausible `po_max` (0.9939) is also the one
with the **best** holdout score (−4.5σ rather than +22σ).

### Steps

1. Add an absolute `po_max` anchor to the comparison/fitting path, **or** state in code why
   the comparison path should not have one. If you add it, use the existing term rather than
   writing a second copy — the anchor living in one of two places is exactly error **E12**,
   and duplicating it would be **E22**.
2. **Make it visible in the result.** Any fit whose derived `po_max` falls outside
   `FIT_RANGES` should say so loudly; a fit at a parameter bound is a failed optimisation
   reported as a number, and nothing currently flags it.
3. Add a guard: **no fit may be reported as a result with a free parameter sitting at its
   bound** without that being stated. `log10_D = −3.0000` appeared three times in brief 2 and
   was reported as a finding each time.
4. Re-run the three Task D fits and print the before/after `po_max`. **If anchoring changes
   the option ranking, say so plainly** — that is the point of the task.

**Done when** `po_max` cannot silently leave `FIT_RANGES`, bound-pinned fits are flagged, and
the three fits are re-reported with the anchor in place.

---

## 3. Task C — redo the `JAHN1997_PEAK_CRC` comparison, and recommend

**Size:** M. **Depends on Task B.** **ENDS IN A RECOMMENDATION TO THE USER, NOT A CHANGE TO
`TRAIN`.**

Brief 2's Task D built the machinery correctly — `get_jahn_crc()`,
`JAHN1997_PEAK_CRC_RESTRICTED` (4 points: 1, 3, 10 µM in the measured window plus the 3000 µM
saturation anchor), `scripts/jahn_options.py`, `TRAIN` left untouched, tests. **Reuse all of
it.** Only the numbers are confounded.

1. Re-run `scripts/jahn_options.py` with Task B's anchor in place.
2. Print per option: does `log10_D` still hit its bound? Is derived `po_max` inside
   `FIT_RANGES`? **What is the paired-pulse holdout residual in units of the published SEM?**
3. Recommend, with the printed table, and **state the cost.** If restricting the CRC to 4
   points makes most downstream fits uninformative, say that — `HANDOFF.md` explicitly allows
   that this may be the honest state.
4. Keep the Mortensen 2012 EC₅₀ cross-check as a **consistency check only**: α1β2γ2**S**
   6.6 µM (pEC₅₀ 5.180 ± 0.0593, n = 34) against Jahn's 11.6 ± 0.9 µM, **not the same
   observable** (20–30 ms delivery latency, so not a fast-application peak) and **no Hill
   coefficient reported**, so it cannot constrain shape. Never pool it.
5. **Stop.** Do not touch `TRAIN`, the P5 comparison or the P6 design.

---

## 4. Task D — Model C's extra parameter buys nothing. Settle it.

**Size:** M. **Depends on Task B.**

Brief 2's Task C reported "the three criteria do not agree" and left it. The numbers say
something sharper than that, and it is worth stating properly:

* **AICc prefers Model B by 12.00 points**, with a train log-likelihood gain for C of
  **0.0005 nats** (35.0668 vs 35.0663) — because C's slow desensitisation parameter sits at
  its lower bound (−4.0).
* **The holdout cannot tell them apart**: −7.6631 vs −7.6638, a difference of **0.0008 nats**,
  both predicting R = 0.1942.
* **CV prefers Model C by 5.6948 nats against a paired SE of 5.6905** — and
  `resolve_cv`'s test is `threshold = max(min_nats, se)` with `if diff < threshold`, i.e.
  **`diff > se`, a 1σ bar.** The margin clears it by **0.0043 nats = 0.08% of its own stated
  uncertainty**, which brief 2's agent traced to a single fold's `sigma_train` numerical floor.

So this is not three criteria disagreeing on equal terms. **Model C's extra parameter sits at
a bound and buys 0.0005 nats in training and 0.0008 nats out of sample; the sole criterion
preferring it does so by 0.08% of a 1σ threshold.**

### Steps

1. **Strengthen `resolve_cv`'s paired-SE test.** A 1σ bar admits roughly a one-in-three chance
   of noise, which is far too permissive for a tie-break the project calls load-bearing. Its
   own docstring already explains why totals-without-SE is wrong; the same reasoning applies
   to SE-at-1σ. Report `diff / se` explicitly in the verdict so a razor-thin margin is
   visible instead of inferred.
2. **Fix or floor the `sigma_train` numerical floor** that drove the fold. A selection margin
   that turns on a numerical floor is not a measurement.
3. Re-run with Task B's anchor. Report whether CV still prefers C.
4. **State the conclusion about Model C plainly.** If its extra parameter remains at a bound
   with no measurable gain anywhere, say that the data does not support it — and note that
   this is the scenario `resolve_cv`'s docstring was written for, where C *nests* B and "best
   total alone always drifts upward in complexity".

**Done when** the verdict prints `diff / se`, a sub-2σ margin cannot name a winner, and the
commit states whether Model C is supported.

---

## 5. Task E — record what P6 actually found: no discriminating observable

**Size:** S. **Mostly writing, and the writing is the finding.**

Brief 2's Task E is honest and its conclusion is stronger than how it was filed. It reports:

* **PEAK does not discriminate** with the application duration declared — peak separation
  under **0.006** in open probability, score 0.688 at σ = 0.005.
* **Paired-pulse recovery does not discriminate either** — corrected score **0.000** at the
  published σ = 0.03 *and* at σ = 0.005.
* Verdict filed: "σ ≤ 0.005 requirement SURVIVES".

That verdict is true but buries the result. **Tightening σ does not help if no observable
separates the models at any σ tried.** The honest statement is: *two* candidate observables
have now been tested with their protocols declared, and **neither discriminates** — so the
pre-registration's blocker is no longer noise, it is the absence of a discriminating
observable.

### Steps

1. Rewrite `knowledge/13-preregistration.md` §5.1 to lead with that. Keep the old and new
   numbers side by side — a pre-registered conclusion that changes *is* the result.
2. Make the list of observables tried, and their outcomes, a **printed table** from a script,
   not prose.
3. Cross-reference Task D: if Model C's extra parameter sits at a bound with no measurable
   gain, then **B and C may be empirically indistinguishable for reasons that have nothing to
   do with σ**, and the P6 design question may be malformed rather than merely hard. Say so if
   the numbers support it. Do not overclaim if they do not.

---

## 6. Do not do these

* **Do not fit against the holdout.** `holdout_guard` raises on the fitting path and that is
  correct — I confirmed it fires. Score with `score_holdout`.
* **Do not change `TRAIN`, the P5 comparison or the P6 design** in Task C. Recommend.
* **Do not loosen a bound to stop a fit pinning against it.** The bound is not the problem;
  an unconstrained direction is. Task B is the fix.
* **Do not regrade `po_max`'s basis** (closed `FITTED`, P0-13) or **subdivide `DataKind`**
  (closed; `source_measured_domain` carries coverage).
* **Do not widen the split-R̂ gate** in `test_posterior_recovers_known_parameters` — a
  marginal failure on Apple Silicon is a known arm64-vs-x86 difference.
* **Do not chase the blocked literature** (Jahn full text, the one-page Hamill comment,
  Barberis Figures 4C/6C) or **invent a mean open time for α1β2γ2L** — no source states one.

---

## 7. Summary

| # | task | depends on | size |
|---|---|---|---|
| A | ~~push~~ done during review — check CI on `b81b0e9` | — | XS |
| B | **`po_max` anchor missing from the comparison path** (confounds brief 2's Task D) | A | M |
| C | redo the Jahn option comparison, recommend to the user | B | M |
| D | strengthen `resolve_cv`'s 1σ bar; settle Model C | B | M |
| E | record that no observable discriminates, not just that σ survives | — | S |

Task B first now, because it is the measurement instrument
for C and D — as the SEM fix was for brief 2.
