# Agent brief 4 — are the datasets mutually consistent?

**Written by:** Claude Opus 5, 2026-10-16, after verifying brief 3's execution.
**Self-contained: assume you know nothing.** One agent, one worktree, one commit per task.
**Push at the end of EVERY task, not at the end of the brief.**

```bash
git worktree add ../cp-brief4 -b claude/agent-brief-4 claude/quirky-sagan-edtdet
cd ../cp-brief4
python3.11 -m venv .venv && source .venv/bin/activate   # 3.13 does NOT work
pip install -e ".[dev]"
pytest -q -m "not slow" -p no:randomly      # green before you start
```

The rules in `AGENT-BRIEF.md` §Rules apply in full — **re-read them.** Especially: never skip
a test; verify with **bare `pytest`** too; nothing in `src/` imports from `scripts/`; numbers
printed by code, never typed into prose; **record your own errors in the commit message**;
cost-model before you profile; **do not loosen a bound to stop a fit pinning against it.**

---

## 0. Brief 3 was done well. Verified.

All four tasks landed and pushed. Checked independently, not read off the commit messages:

* **Task B's anchor works, and its own self-caught error is the kind worth copying** — the
  fixed stated-SEM `po_max` anchor was being charged as an *estimated* noise parameter in `k`,
  penalising anchored models; found and fixed by the agent, recorded in the message.
  `ModelFit` now carries `at_bound`, `po_max`, `po_max_in_range` and a `k_effective` that
  discounts bound-pinned parameters. **That last idea is better than what I asked for**: a
  parameter stopped by its box was never estimated and should not be charged to `k`.
* **Task D vindicated the diagnosis.** With the anchor in and `sigma_train` floored, CV's
  margin for Model C collapsed from **+5.6948 nats to +0.0124** (diff/se = 1.41σ against the
  new 2σ bar). Model C is settled: **not supported** — its nested parameter is pinned at
  `-4.0000` with −0.0005 nats of train gain and no out-of-sample gain.
* **Task E's conclusion is the right one** and states it plainly: no candidate observable
  discriminates B from C at any σ tried, because **when `D_slow → 0` Model C collapses onto
  Model B.** The blocker is identifiability, not measurement scatter, "rendering the P6
  design question malformed in this regime."
* **Task C recommended Option 3 and left `TRAIN` untouched**, as instructed. Verified: `TRAIN`
  still has 4 datasets with the 9-point Jahn CRC; `HOLDOUT` still has 1.

**Two of my own suspicions were wrong and I checked them rather than filing them.** I thought
predicted recovery might be inverted in desensitisation — it is correctly monotone (R falls
from 0.995 to 0.033 as `d` goes 1e-4 → 5). I then thought the CRC-alone fit's R = 0.1942 at
`log10_D = −3` was inconsistent — it reproduces to three decimals at
`(log10_kd, log10_E, log10_D) = (2.0, 2.21, −3.0)`, where low recovery comes from slow
*reset*, not from desensitisation. No bug either time.

---

## 1. The question this brief exists to answer

The holdout misses by **~17σ** when the model is fitted to `TRAIN`, and by **−4.5σ** when
fitted to the CRC alone. Both are bad. The obvious hypothesis is that the scheme simply cannot
produce the observed recovery — **and that is false.** I checked it on a coarse
11 × 11 × 13 grid, restricted to `po_max` inside `FIT_RANGES` (0.70, 0.80):

| | |
|---|---|
| grid points with `po_max` in range | 143 |
| achievable R over those | **0.0643 to 0.9998** |
| observed R | 0.33 ± 0.03 |
| closest grid point | **1.32σ**, R = 0.3695, `po_max` = 0.7597, at `(log10_kd, log10_E, log10_D) = (3.00, 0.50, 2.50)` |

**So a region of parameter space satisfies the holdout to ~1.3σ with `po_max` in range, and
the fit to `TRAIN` does not find it.** That is not a falsification and not an optimiser
failure to assume — it is evidence that **the training data and the holdout may be pulling in
different directions**, i.e. that the datasets are mutually inconsistent.

That matters more than any remaining fit, because if `TRAIN` and the holdout cannot be
satisfied together, then **no posterior fitted on `TRAIN` means anything about this receptor**,
and the project needs to know which dataset is the problem. There are strong candidates: the
Jahn CRC is a **parametric reconstruction with 6 of its 9 points outside the 1–10 µM window
its source measured**, and every real dataset in the repo is **γ2S** while the project's
receptor is **γ2L**.

---

## 2. Task A — the dataset-consistency diagnostic

**Size:** M. **The core of this brief.**

### What to build

Find the parameter vector that best satisfies **all three** real datasets *jointly* — the Jahn
CRC, the Barberis deactivation and the Barberis paired-pulse recovery — and report **how much
each one has to be compromised** to accommodate the others.

### THIS IS A DIAGNOSTIC, NOT A FIT. Read this twice.

* **`TRAIN` does not change. `HOLDOUT` does not change. No result from this task is adopted
  as a fitted model.** The output is a table of residuals and the words "these datasets are /
  are not mutually consistent".
* It is **not** training on the holdout. The distinction is that nothing here is kept: you are
  asking *"does any parameter set satisfy everything?"*, which is a property of the data, not
  a model you then use. State that in the script's docstring and in the commit.
* `holdout_guard` will raise if you route the holdout through the fitting path — **correctly**,
  I confirmed it fires. Build the diagnostic objective explicitly rather than defeating the
  guard, and **do not weaken or bypass `holdout_guard`.** If you need a separate entry point,
  add one that is clearly named as a diagnostic and that `TRAIN`-path callers cannot reach.

### Steps

1. Build a joint objective over the three real datasets, each scored with **its own declared
   `noise_model`** (the holdout has `stated_sem = 0.03`; the other two are `"estimated"`), plus
   the `po_max` anchor from brief 3 Task B.
2. Minimise it. **Cost-model first** — a PEAK prediction integrates a waveform per
   concentration, so time one evaluation and multiply before choosing a grid or restart count.
   Start from the grid point above, which is known to satisfy the holdout.
3. Print, for the joint optimum and for each single-dataset optimum, a square table: **rows =
   which datasets were fitted, columns = residual on each dataset in units of its own σ.** The
   diagonal is the fit quality; the off-diagonal is the conflict.
4. State the verdict: is there a parameter set that satisfies all three within, say, 3σ each?
   **If not, which pair is irreconcilable, and by how much?**
5. Report `at_bound` and `po_max_in_range` for every fit in the table. A bound-pinned joint
   optimum means the conflict is being absorbed by the box, not resolved.

**Done when** a script prints that square table and the commit states whether the three real
datasets can be satisfied together.

---

## 3. Task B — isolate which dataset pulls

**Size:** S–M. **Depends on Task A.** Mostly reading Task A's table, then one more run.

Two numbers already point somewhere: fitted to the **CRC alone** the holdout misses by
**−4.5σ** (under-recovers, R = 0.1942); fitted to **CRC + deactivation** it misses by **+17σ**
(over-recovers, R ≈ 0.84). Adding the deactivation **flips the sign** of the holdout error and
makes it four times worse.

1. Score the holdout after fitting to: CRC alone; deactivation alone; both. Three rows.
2. Report whether the deactivation and the holdout are individually compatible — they come
   from **the same paper, the same preparation and the same cells** (Barberis 2007, rat
   α1β2γ2S, 22–24 °C), so if *those two* cannot be satisfied together that is a much stronger
   statement than any disagreement involving the Jahn reconstruction.
3. If instead the Jahn CRC is the outlier, say so — and note that brief 3 Task C recommended
   keeping it (Option 3) because restricting it made `log10_D` indeterminate over 4.4 decades.
   **A dataset that cannot be reconciled with two others from a real preparation is a
   different problem from an indeterminate one**, and would change that recommendation.

**Done when** the three-row table is printed and the commit names which dataset is in tension
with which.

---

## 4. Task C — write the γ2S substitution into the conclusions

**Size:** S. **Independent; can be done any time.**

Every real dataset in the repo is **α1β2γ2S**. The project's receptor is **α1β2γ2L**. The
substitution is recorded per dataset, but it is **not** carried into the conclusions, and after
Task A it may be the explanation for whatever is found.

1. From `knowledge/14-literature.md`: only **two** sources in the corpus of seven are γ2L (Li
   2008, rat; Dixon 2014, human), and **neither** yields a usable time course — Dixon publishes
   only a weighted scalar with no components, so it cannot be reconstructed at all.
2. Add a short, explicit limitation section to `knowledge/08-manuscript.md` and the roadmap:
   every quantitative result on this branch rests on **γ2S data applied to a γ2L model**, and
   the two splice variants are not assumed equivalent anywhere — it is a recorded substitution,
   not an identity.
3. **Do not** try to correct for it. There is no basis to.

---

## 5. Task D — land the settled conclusions

**Size:** S. **Do last.**

Three questions are now closed and should read as closed rather than open:

1. **Model C is not supported.** Its nested parameter is pinned at its bound with no gain in
   training, CV (now 1.41σ against a 2σ bar) or out of sample, and it collapses onto Model B as
   `D_slow → 0`.
2. **The P6 design question is malformed in this regime** — no observable can discriminate two
   models that are not distinguishable. Record it as malformed, not as "needs lower σ".
3. **`JAHN1997_PEAK_CRC` keeps its 9 points** (Option 3), with the restricted variant retained
   to demonstrate the artefact. **Subject to Task B**: if the CRC turns out to be irreconcilable
   with the two Barberis datasets, reopen this and say so.

Update `knowledge/09-strategic-roadmap.md` and `knowledge/08-manuscript.md`, and regenerate the
manuscript's numbers with `scripts/paper_numbers.py` rather than editing them by hand.

---

## 6. Do not do these

* **Do not fit against the holdout, and do not weaken `holdout_guard`.** Task A is a
  diagnostic whose output is discarded; if that distinction ever becomes unclear in the code,
  you have built the wrong thing.
* **Do not change `TRAIN` or `HOLDOUT`.**
* **Do not loosen a bound** to stop a fit pinning against it, and do not widen
  `FIT_RANGES["po_max"]` to make `po_max = 0.59` acceptable. Brief 3's anchored fits land at
  0.59, outside (0.70, 0.80); that is a finding to report, not a range to adjust.
* **Do not widen the split-R̂ gate** in `test_posterior_recovers_known_parameters` — a marginal
  failure on Apple Silicon is a known arm64-vs-x86 difference.
* **Do not revisit** `po_max`'s basis (closed `FITTED`), `DataKind` subdivision (closed), or
  Model C's support (closed by brief 3 Task D).
* **Do not chase the blocked literature** (Jahn full text, the one-page Hamill comment,
  Barberis Figures 4C/6C) or **invent a mean open time for α1β2γ2L** — no source states one.

---

## 7. Summary

| # | task | depends on | size |
|---|---|---|---|
| A | **can the three real datasets be satisfied together?** (diagnostic, nothing adopted) | — | M |
| B | isolate which dataset pulls — deactivation vs holdout is the sharp pair | A | S–M |
| C | write the γ2S→γ2L substitution into the conclusions | — | S |
| D | land the settled conclusions (Model C, P6 malformed, Jahn Option 3) | B | S |

Task A first. If the three real datasets cannot be satisfied together, that outranks every
remaining fit on this branch — and if they can, the ~17σ holdout miss is a fitting problem with
a known target to aim at.
