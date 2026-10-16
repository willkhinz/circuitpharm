# Agent brief 2 — five tasks, in order

**Written by:** Claude Opus 5, 2026-10-15, after verifying `AGENT-BRIEF.md`'s execution.
**Self-contained: assume you know nothing.** One agent, one worktree, one commit per task.

```bash
git worktree add ../cp-brief2 -b claude/agent-brief-2 claude/quirky-sagan-edtdet
cd ../cp-brief2
python3.11 -m venv .venv && source .venv/bin/activate   # 3.13 does NOT work
pip install -e ".[dev]"
pytest -q -m "not slow" -p no:randomly     # ~75 s, 498 tests — should be green before you start
```

Read first: `knowledge/14-literature.md`, then `knowledge/11-identifiability.md` §1.2.
`PLAN-phase3.md` is the wider plan. The rules in `AGENT-BRIEF.md` §Rules still apply in full
— **re-read them**, especially: never skip a test; verify with **bare `pytest`** as well as
`python -m pytest`; nothing in `src/` imports from `scripts/`; numbers are printed by code,
never typed into prose; one canonical name per quantity; absence of a declaration is an error,
not a default; a test must check the invariant its name claims; **record your own errors in
the commit message**; cost-model before you profile.

---

## 0. What the last brief achieved, verified

All four tasks landed (`538c879`, `1648538`, `f54716a`, `64ee237`, `640bf65`). Tree clean,
498 fast tests green. Verified independently, not taken on trust:

| | |
|---|---|
| Task A | Two commits — E27 substantive, lint separate. Tree clean. ✓ |
| Task B | `source_measured_domain` added **and enforced**. Jahn `(1.0, 10.0)` → **6 of 9** points outside (0.667); Barberis `(0.0, 1000.0)` → **0 of 101**. ✓ |
| Task C | P0-13 closed as DECIDED with its reasoning. ✓ |
| Task D | `PAIRED_PULSE_RECOVERY` + `PairedPulseDataset`; `BARBERIS2007_PAIRED_PULSE` registered `role="holdout"`, `kind="digitised"`, **`synthetic=False`** — the project's **first non-synthetic dataset**. `HOLDOUT` non-empty. ✓ |

**The scientific finding is real and it is the most important result on this branch.** Models
fitted only to the peak CRC **do not generalise out of sample**, and the mechanism is
identified: Jahn's steep nH = 2.2 drives desensitisation `D` to its lower bound, and
paired-pulse recovery is *governed* by desensitisation. Reproduced here independently:

| model | predicted R | observed | miss |
|---|---|---|---|
| Model B, unfitted defaults | 0.8893 | 0.33 ± 0.03 | **18.6× SEM** |
| Model C, unfitted defaults | 0.8529 | 0.33 ± 0.03 | **17.4× SEM** |
| B and C fitted to Jahn CRC | 0.1942 | 0.33 ± 0.03 | **4.5× SEM** |

So the Jahn reconstruction artefact now has a **measured out-of-sample consequence**, which
is exactly what the holdout was for.

---

## 1. Task A — the holdout score ignores the published SEM. Fix it.

**Size:** M. **Do first: every other quantitative claim about the holdout depends on it.**

### The defect

`likelihood._score` estimates the noise scale from the residuals, for every dataset,
unconditionally:

```python
sigma2 = max(rss / n, 1e-12)
ll = -0.5 * n * (np.log(2.0 * np.pi) + 1.0 + np.log(sigma2))
```

For a **single-point** holdout that means estimating a variance from **one** residual, and
the consequence is not subtle. With R = 0.1942 against 0.33 ± 0.03:

| σ used | holdout log-likelihood |
|---|---|
| σ estimated from the residual (**what the code does**) | **+0.577** |
| the published SEM, 0.03 (**what the data says**) | **−7.66** |

Verified to four decimals against the committed value. **A 4.5σ miss is being scored as
`+0.577`, which reads as a good fit.** And it is why Models B and C score 0.5774 and 0.5773 —
indistinguishable, because `sigma_hat` absorbs the residual rather than measuring it.

`BARBERIS2007_PAIRED_PULSE` **already records `sem=0.03` and `n_cells=8`**. The published
uncertainty is present and unused.

### Why this is not just a holdout problem

Estimating σ from residuals is defensible for a multi-point reconstructed CRC, where the
project may not trust the error bars. It is **not** defensible for a stated scalar with a
published SEM from n = 8 patches, and it silently makes holdout log-likelihoods
**incomparable across datasets of different precision** — which is the whole point of a
holdout score.

### Steps

1. Give `_score` a per-dataset choice of noise model, explicit rather than implicit:
   * **use the dataset's stated `sem`** when it has one and it is trustworthy (the new
     holdout);
   * **estimate from residuals** otherwise (the existing CRCs).
   Make the choice a **declared field on the dataset**, not an inference from which
   attributes happen to exist — rule 6. Suggested: `noise_model: Literal["stated_sem",
   "estimated"]`, required where a `sem` exists so no dataset gets one by accident.
2. **Do not silently change the existing CRC scores.** Those feed `compare_models`, the CV
   margin and the AICc. Keep them on `"estimated"` and say so in the commit; if any number
   moves, print the before/after.
3. Report `sigma_hat` *and* the stated `sem` on `DatasetFit` so a reader can see which was
   used and how far apart they are.
4. Re-run the holdout evaluation. **The 4.5× SEM miss should now show up as a strongly
   negative log-likelihood.** Print it; do not retype it.
5. Tests: a stated-SEM dataset scores `−7.66` (within tolerance) rather than `+0.577`; a
   dataset declaring `"stated_sem"` with no `sem` raises; the existing CRC path is unchanged.

**Done when** the holdout log-likelihood reflects the published measurement error, the choice
is declared per dataset, and no existing score moved without being printed.

---

## 2. Task B — put E27 in the register

**Size:** XS.

`tests/test_review_regressions.py:827` registers **E27** (deactivation protocol-dependence)
and its test is green. But the register itself is **`WORKLOG.md`**, which holds E1–E26, and
**E27 is absent from it.** The registration is half-done: the guard exists, the record does
not.

### Steps

1. Add E27 to `WORKLOG.md` in the format E26 uses (read it first — it is a narrative entry
   with a verification line, not a table row).
2. Content: deactivation is protocol-dependent, and this is a **published measurement**, not
   a model artefact — Barberis 2007 measured the same receptor in the same study at
   **52.5 ms (2 ms pulse)** and **364 ms (3 s pulse)**, 6.9× apart. Dixon 2014's 5.9 ms
   differs from Barberis's 2 ms-pulse value by 8.9× at **nearly the same pulse duration**,
   because Dixon fitted two exponentials where Barberis fitted three and Barberis's 221 ms
   component carries 83% of its weighted τ. So a weighted deactivation τ is undefined until
   **pulse duration, component count and fit window** are all declared.
3. Cross-reference `knowledge/14-literature.md` §3 and the guard test.

**Done when** `grep -c E27 WORKLOG.md` is non-zero and the entry says what the error was.

---

## 3. Task C — resolve CV versus AICc, or prove it cannot be resolved

**Size:** M. **This is now load-bearing, because a holdout exists.**

The Task D agent reported that after fitting to the Jahn CRC, **cross-validation prefers
Model C while AICc prefers Model B** — unresolved. **That was reported, not independently
reproduced; verify it first.**

This matters more than it looks. `HANDOFF.md` §6 records that the CV tie-break is
load-bearing: Model C *nests* Model B and once scored 0.84 nats better out of sample **on
noise**, and `resolve_cv` requires the margin to clear both a paired per-fold SE and a 2-nat
floor before breaking ties on fewest parameters. Remove that and the protocol selects the
larger model on noise.

### Steps

1. Reproduce the disagreement and print both criteria with their margins. **Cost-model
   first** — a PEAK prediction integrates a waveform per concentration (~40× a deactivation
   score), so a naive sweep is hours.
2. Establish whether the disagreement survives `resolve_cv`'s own guards (paired SE, 2-nat
   floor). If it does not, there is no disagreement to resolve and that is the answer.
3. Now that a **real holdout** exists, score both models on it and report whether the holdout
   agrees with CV, with AICc, or with neither. Use Task A's corrected noise model — this is
   why Task A goes first.
4. If all three disagree, **say so and stop.** Do not invent a tie-break. Three criteria
   disagreeing on nested models fitted to one parametric CRC is a finding about the data, not
   a problem to be optimised away.

**Done when** a script prints CV margin, AICc margin and holdout score for B and C side by
side, and the commit states plainly whether they agree.

---

## 4. Task D — the `JAHN1997_PEAK_CRC` decision, now with evidence

**Size:** M. **Depends on Tasks A and C.** **ENDS WITH A RECOMMENDATION TO THE USER, NOT A
CHANGE TO `TRAIN`.**

This is Phase 3 of `PLAN-phase3.md` and it is **much more decidable than when it was
written**, because two things now exist that did not:

* **Task B of the last brief made the artefact measurable:** 6 of Jahn's 9 points lie outside
  the 1–10 µM window its source measured.
* **The holdout made its consequence measurable:** fitting to that CRC drives `D` to its
  lower bound and produces a **4.5σ** out-of-sample failure on paired-pulse recovery.

That is a far stronger argument for restricting the dataset than anything available before.
**It is still not your decision to take.**

### The three options (`HANDOFF.md` §5.1)

1. **Keep as is** — only sourced CRC in the repo; artefact documented and now measured.
2. **Restrict to the measured window plus a saturation anchor** — asserts only what the
   source supports; leaves ~3 points, and most downstream fits may become uninformative,
   *which may be the honest state.*
3. **Split in two** — a defensible restricted dataset plus the full reconstruction retained
   to *demonstrate* the artefact.

### Steps

1. Implement all three behind a flag.
2. Fit under each, **with `barberis2007_deact` included** (an absolute rate is available now —
   that was the whole reason this phase was moved after the deactivation work).
3. For each option print: does `D` still hit its bound? Is the plateau still 0.99? **What is
   the paired-pulse holdout residual, in units of the published SEM?** That last column is
   the one that did not exist before and it is the most informative.
4. Note the independent EC₅₀ check: Mortensen 2012 gives α1β2γ2**S** EC₅₀ = **6.6 µM**
   (pEC₅₀ 5.180 ± 0.0593, n = 34) against Jahn's 11.6 ± 0.9 µM. **Not the same observable** —
   Mortensen delivers GABA with a 20–30 ms latency, so it is not a fast-application peak, and
   6.6 < 11.6 is the direction this project's own PEAK protocol-dependence predicts. A
   consistency check on the EC₅₀ anchor, **never a second point to pool.** Mortensen reports
   **no Hill coefficient**, so it cannot constrain shape.
5. Recommend, with the printed table, and **state the cost plainly.** If the honest answer is
   "restricting it makes most fits uninformative", say that — `HANDOFF.md` explicitly allows
   that this may be correct.
6. **Stop.** Do not touch `TRAIN`, the P5 comparison or the P6 design.

---

## 5. Task E — re-derive the P6 σ requirement

**Size:** M–L (MCMC). **Depends on Tasks A, C, D.** Attempt only if the earlier tasks are
landed and green.

`knowledge/13-preregistration.md` concludes the discriminating experiment does not exist at
σ = 0.02 and needs **σ ≤ 0.005**. That rests on posteriors fitted partly through `PEAK`,
which was later found protocol-dependent: application length alone moves peak EC₅₀
**sevenfold** (299 µM at 2 ms → 39 µM at 1000 ms), and `PEAK_APPLICATION_MS` is fixed at
300 ms while no dataset declares its own.

### Steps

1. Re-run with the application duration **declared and matched** between model and dataset.
   Use `models.base.hill_slope(..., window_rel=...)` for any slope comparison.
2. Report whether σ ≤ 0.005 **survives.** It may soften — some apparent indistinguishability
   may be protocol smearing rather than a property of the schemes.
3. If it softens, **do not quietly relax the pre-registration.** A pre-registered conclusion
   that changes *is* the result: record old and new side by side with the protocol stated for
   each.
4. **Re-target the discriminating observable onto the paired-pulse recovery.** P6 originally
   chose Barberis's desensitisation onset, but that **cannot be fully reconstructed** (τ₂/τ₃
   appear only in Figure 4C) and the source is **internally inconsistent about its n** (text
   says 7, the Fig. 4 caption says 3). The paired-pulse recovery is fully stated, is already
   registered, and — given it separates fitted from unfitted models by 4.5σ versus 18.6σ —
   may discriminate better anyway.
5. Apple Silicon: a marginal split-R̂ failure in `test_posterior_recovers_known_parameters` is
   a known arm64-vs-x86 difference through the matrix exponentials. **Do not widen that gate.**

---

## 6. Do not do these

* **Do not change `TRAIN`, the P5 comparison or the P6 design** as part of Task D. The
  recommendation goes to the user.
* **Do not "fix" the out-of-sample failure by refitting against the holdout.** It is a
  holdout. Fitting to it destroys the only out-of-sample evidence the project has, and
  `fitting/data.py` already guards this path — read the message at its `"Pass
  fitting.data.TRAIN"` check before going near it.
* **Do not regrade `po_max`'s basis.** Closed as `FITTED` in P0-13 (`64ee237`), with reasons.
* **Do not subdivide `DataKind`.** Closed; `source_measured_domain` carries coverage instead.
* **Do not chase the blocked literature.** Jahn 1997's full text and the one-page Hamill
  comment (NeuroReport 8(16):iv) need institutional Ovid access, not more searching. Barberis
  Figures 4C/6C would need figure digitisation — a *lower* evidence grade than reading a
  table, and must be graded as such if ever done.
* **Do not invent a mean open time for α1β2γ2L.** No source states one; four readings span
  5.1× and the spread is *definition*, not disagreement. `FIT_FIXED_ALPHA` stays a
  `CONVENTION`.

---

## 7. Summary

| # | task | depends on | size |
|---|---|---|---|
| A | **holdout score must use the published SEM** (+0.577 should be −7.66) | — | M |
| B | E27 into `WORKLOG.md` | — | XS |
| C | resolve CV vs AICc, or prove it cannot be | A | M |
| D | `JAHN1997_PEAK_CRC` recommendation, with holdout evidence | A, C | M |
| E | re-derive P6 σ, re-targeted onto paired-pulse | A, C, D | M/L |

Task A first: it is the measurement instrument for everything after it.
