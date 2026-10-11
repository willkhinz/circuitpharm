# Plan for the next upgrade — written for delegated agents

**Author:** Claude Opus 5, 2026-10-10, after reading the literature (commit `48da27c`).
**Status: NOT STARTED. Awaiting the user's go-ahead.**

Read `HANDOFF.md`, then `knowledge/14-literature.md`, then this file. The literature pass
**changed the right order of work**, so this plan supersedes `HANDOFF.md` §5 where they
disagree. The reason is in §0.

---

## 0. Why the order changed

`HANDOFF.md` §5 said: decide the `JAHN1997_PEAK_CRC` question first. That was correct when
the only dataset in the repo was a parametric reconstruction of a four-decade curve from a
one-decade slope — the question was how much to *subtract*.

Reading the literature did **not** produce a digitisable concentration–response curve (no
source gives one; Mortensen 2012 gives an EC₅₀ with no Hill coefficient and no curve
points), so **the Jahn question does not dissolve**. But it produced something better:

> **Barberis 2007's deactivation after a 2 ms pulse is fully specified by the source** —
> three time constants, three areas, all with errors, in the functional form the authors
> themselves fitted. Recomputing Σ A_i τ_i gives 53.55 ms against a published 52.5 ± 2.9 ms.
> It can be reconstructed with **no extrapolation beyond the measured window**, which is
> exactly what the Jahn CRC cannot claim.

Deactivation is the observable that makes an **absolute rate (`k_off`) identifiable**;
without it the fit determines only ratios (`knowledge/11-identifiability.md`). The
`D`-on-its-bound, 0.99-plateau artefact that motivated the Jahn question is plausibly a
*symptom of having no absolute-rate information at all*.

**So: add the information before deciding what to subtract.** Phase 3 (the Jahn decision)
must be run *after* Phase 2, and may well get a different answer than it would today. If
Phase 2 is skipped, Phase 3 is decided on strictly less evidence for no reason.

---

## 1. Rules every agent must follow

These are not style preferences; each one was learned from a failure recorded in this repo.

1. **Never skip, disable, `xfail` or quarantine a test to get green.** If a test fails,
   either the code is wrong or the test's invariant has changed — and if the latter, rewrite
   the invariant *and say so in the commit*. See `tests/test_data_safety.py`'s own docstring
   for a worked example of doing this correctly.
2. **Numbers are printed by code, never typed into prose.** Any figure quoted in a markdown
   document must be emitted by a script in `scripts/`. `scripts/paper_numbers.py` exists
   because a draft once quoted a self-consistent set of numbers that no commit produced.
3. **`[LEAD]` ≠ `[VERIFIED]`.** A lead is a place to look, never a number to use. Promote
   only after reading the full text yourself. **Never fabricate a citation or a number.** An
   honest "the full text does not state this" is worth more here than a plausible value.
4. **A resolving DOI proves existence, not support.** Any source row asserting
   `CLAIM_SUPPORT_SUPPORTS` needs a `HAND:` verdict saying who assessed it and which
   specific claims it supports. `tests/test_data_safety.py` enforces this.
5. **State the protocol, always.** `PEAK` is not one observable; neither, as of this pass,
   are `DEACTIVATION` or `MEAN_OPEN_TIME`. Use
   `models.base.hill_slope(..., window_rel=...)` for any slope comparison — the window is an
   argument precisely because a note in prose failed to constrain the computation once.
6. **Record your own errors in the commit message** rather than quietly fixing them. This is
   load-bearing convention in this repo.
7. **Absence of a string is not absence of a value, and presence is not presence.** A
   case-insensitive grep for `nH` in Mortensen 2012 returns 34 hits, all author initials.
   Compute a quantity two ways before concluding it is missing or present.
8. **Error statistics are not assumed.** Keramidas 2008 and Dixon 2014 declare neither SEM
   nor SD. They are recorded `UNSTATED`. Do not "tidy" these into SEM.

### Environment

```bash
python3.11 -m venv .venv && source .venv/bin/activate   # 3.13 does NOT work
pip install -e ".[dev]"
pytest -q -m "not slow" -p no:randomly    # ~2 min, 473 tests
pytest -q -p no:randomly                  # full, ~20 min; run before any push
```

* **`-p no:randomly` whenever comparing runs.** Test order is randomised by default.
* **Apple Silicon:** a marginal split-R̂ failure in
  `test_posterior_recovers_known_parameters` is a known arm64/Accelerate vs x86/OpenBLAS
  difference in the third decimal through the matrix exponentials. The chains were lengthened
  (12000/4000 → 24000/8000) rather than the gate loosened. **Do not widen the gate.**
* `ProvisionalResultWarning` when fitting a `synthetic=True` dataset is the machinery
  working, not a bug.
* Cache simulation traces on the model source hash — analysis iteration then costs seconds.

### ⚠️ Concurrency — read this before starting

On 2026-10-10 debugging agents modified 37 files in this working tree *while another agent
was mid-task*, including an unused-import sweep and a `pythonpath` addition to
`pyproject.toml`. Nothing broke, but the commit had to be assembled hunk-by-hunk to avoid
committing another agent's in-flight work.

**Therefore: one agent per phase, each in its own git worktree.**

```bash
git worktree add ../cp-phase2 -b claude/phase2-deactivation claude/quirky-sagan-edtdet
```

Phases 1–4 are a dependency chain and must be **sequential**. Phase 5 is independent and may
run in parallel in its own worktree. Do not run two agents in one checkout.

---

## 2. Phase 1 — Declare the protocol on the two new protocol-dependent observables

**Depends on:** nothing. **Blocks:** Phases 2, 3, 4. **Size:** small. **Risk:** low.

The project already learned this lesson once for `PEAK` and fixed it with
`PEAK_APPLICATION_MS` and `hill_slope(window_rel=)`. Do the same for deactivation and mean
open time, where the literature has now shown the same failure.

### Steps

1. Read `knowledge/14-literature.md` §3 and §2, and `scripts/literature_gamma2.py`.
2. Find how `PEAK_APPLICATION_MS` is declared and threaded (`src/circuitpharm/config.py`,
   `protocols/waveforms.py`, `fitting/data.py`). **Follow that existing pattern** — do not
   invent a second mechanism. Note commit `6e15266` ("Extend the E22 guard to subpackages;
   it found a fourth copy of the pulse"): check for duplicated constants before adding one.
3. Add to the dataset/observable description, as **arguments or dataclass fields, not
   comments**:
   * deactivation: `pulse_ms`, `conc_mm`, `n_components`, `fit_window_ms`
   * mean open time: `conditioning` — one of a small closed set, e.g.
     `INTRABURST`, `INTRACLUSTER`, `NONSTATIONARY_WINDOW`, `POPULATION`
4. Make the absence of a declaration an **error, not a default**. The `PEAK` lesson was that
   a silent default is what allowed the mismatch. A dataset tagged `DEACTIVATION` with no
   `pulse_ms` must fail loudly.
5. Add a guard test in the style of the existing E22 guard: no module may compare two
   deactivation values whose declarations differ.

### Acceptance

* `pytest -q -m "not slow" -p no:randomly` green.
* A new test demonstrates that comparing Dixon's 5.9 ms with Barberis's 52.5 ms **raises or
  is flagged**, citing their differing `pulse_ms`/`n_components`.
* No new copy of an existing constant (grep for duplicates first).

---

## 3. Phase 2 — Enter the Barberis deactivation as the first non-synthetic dataset

**Depends on:** Phase 1. **Blocks:** Phases 3, 4. **Size:** medium. **Risk:** medium —
**this is the highest-value step in the plan.**

This is what makes `k_off` identifiable and would let the first
`ProvisionalResultWarning` stop on its own merits.

### The data (all from `scripts/literature_gamma2.py`; do not retype — import it)

Barberis et al. 2007, Eur J Neurosci 25(9):2726–2740, PMC1950087, rat **α1β2γ2S**, HEK293
outside-out pooled with small lifted whole cells, **22–24 °C**, −70 mV, ultrafast exchange
60–100 µs, **2 ms pulse of 10 mM GABA**:

| component | τ (ms) | area |
|---|---|---|
| fast | 2.8 ± 0.3 | 0.57 ± 0.04 |
| middle | 33.4 ± 4.6 | 0.23 ± 0.02 |
| slow | 221.35 ± 14.9 | 0.20 ± 0.03 |

τ_w = 52.5 ± 2.9 ms (SEM, n = 6). Σ A_i = 1.000. Σ A_i τ_i = 53.55 ms.

### Steps

1. Build the time course as a sum of three exponentials from **those** parameters, at a
   declared `fit_window_ms`. Generate it in code from the imported constants.
2. **Decide and document the `synthetic` flag.** This is the judgement call of the phase and
   it must be made explicitly, not by default:
   * **For `synthetic=False`:** every parameter of the functional form is stated by the
     source, over the window the source measured, with no extrapolation. The read-back check
     (53.55 vs 52.5 ± 2.9, 0.36 SEM) confirms the set is complete and self-consistent.
   * **Against:** it is a *curve reconstructed from fitted parameters*, not digitised points
     — structurally the same operation as `JAHN1997_PEAK_CRC`, which this project calls a
     reconstruction artefact. And it is **γ2S**, not γ2L.
   * **Recommended:** a **third grade** between the two, e.g.
     `provenance.Basis.SOURCE_PARAMETRIC` or a `reconstruction=` field — "every parameter
     stated by the source, within the measured window" is genuinely different from both
     "digitised from a figure" and "generated from our formula". If a third grade is added,
     `ProvisionalResultWarning` must state which grade it is complaining about, and the Jahn
     CRC must be regraded at the same time, since the distinction is precisely what separates
     them: Jahn extrapolates four decades from one, Barberis extrapolates nothing.
     **Ask the user before inventing a new grade** — it changes what every warning means.
3. Register in `provenance.py` and point `source_key` at `barberis2007` (already in the
   `sources` table with a `HAND:` verdict).
4. Record the **γ2S substitution** explicitly, in the dataset, not only in prose.
5. Re-run the identifiability analysis: does adding this dataset make `k_off` identifiable,
   as `knowledge/11-identifiability.md` predicts? **Report the answer either way.** A
   negative result here is a real finding and must not be buried.

### Acceptance

* Full suite green.
* A script prints the before/after identifiability result for `k_off`, and
  `knowledge/11-identifiability.md` is updated with whatever it actually printed.
* The `synthetic`/grade decision is written down with its reasoning, and the Jahn CRC is
  regraded consistently if a new grade was introduced.

---

## 4. Phase 3 — Decide what `JAHN1997_PEAK_CRC` should be

**Depends on:** Phase 2 (**do not start before it — see §0**). **Size:** medium.

`HANDOFF.md` §5.1 and `knowledge/12-inference.md` §2.4 state the problem: nine points of a
global Hill curve with nH = 2.2 spanning 0.3–3000 µM, four decades of shape extrapolated
from a slope measured over one (1–10 µM). The source constrains three things only: the slope
over 1–10 µM, the EC₅₀, and saturation by 3 mM.

### New evidence from the literature pass

* **No digitisable CRC exists in the corpus.** The Jahn CRC remains the only one. The
  question does **not** dissolve.
* **An independent EC₅₀ now exists:** Mortensen 2012 gives α1β2γ2**S** EC₅₀ = **6.6 µM**
  (pEC₅₀ 5.180 ± 0.0593, n = 34), HEK293 whole-cell — against Jahn's 11.6 ± 0.9 µM. **These
  are not the same observable:** Mortensen delivers GABA with a 20–30 ms latency, so it is
  not a fast-application peak. The direction (6.6 < 11.6) is what the project's own PEAK
  protocol-dependence table predicts for a longer effective application. Use it as a
  *consistency check on the EC₅₀ anchor*, **not** as a second point to pool.
* Mortensen reports **no Hill coefficient**, so it cannot constrain the curve's shape.

### Steps

1. Implement all three options from `HANDOFF.md` §5.1 behind a flag: keep-as-is; restrict to
   the measured window plus a saturation anchor; split into a defensible restricted dataset
   plus the full reconstruction retained to *demonstrate* the artefact.
2. **Fit under each, with the Phase 2 deactivation dataset included**, and report what
   changes: does `D` still land on its bound? Is the plateau still 0.99? Does the restricted
   3-point version remain uninformative once an absolute rate is available?
3. Make a recommendation **with the printed evidence**, and **state plainly what it costs.**
   If the honest answer is "most downstream fits become uninformative", say that — the
   handoff explicitly allows that this may be the correct state.
4. **Do not pick silently.** Bring the comparison back to the user before changing `TRAIN`,
   the P5 comparison or the P6 design.

### Acceptance

* A script prints the fit outcome under all three options, side by side.
* A written recommendation with its cost stated, and the user consulted before `TRAIN` moves.

---

## 5. Phase 4 — Re-derive the P6 σ requirement with the protocol declared

**Depends on:** Phases 1–3. **Size:** medium-large (MCMC). **Risk:** low.

`knowledge/13-preregistration.md` concludes the discriminating experiment does not exist at
σ = 0.02 and needs **σ ≤ 0.005** — a 4× noise reduction. That rests on posteriors fitted
partly through `PEAK`, and `PEAK` was later found protocol-dependent: application length
alone moves peak EC₅₀ **sevenfold** (299 µM at 2 ms → 39 µM at 1000 ms).
`PEAK_APPLICATION_MS` is fixed at 300 ms and no dataset declares its own.

### Steps

1. Re-run the P6 derivation with the application duration **declared and matched** between
   model and dataset.
2. Report whether the pre-registered σ ≤ 0.005 conclusion **survives**. It may soften: some
   of the apparent indistinguishability may be protocol smearing rather than a real property
   of the two schemes.
3. If it softens, do **not** quietly relax the pre-registration. A pre-registered conclusion
   that changes is itself the result; record the old and new values side by side with the
   reason.
4. Barberis's desensitisation onset was P6's candidate discriminating observable. Note from
   §7 of `knowledge/14-literature.md` that it **cannot be fully reconstructed** (τ₂/τ₃ are
   only in Figure 4C) and that its n is ambiguous in the source (text 7, caption 3). The
   **paired-pulse recovery at a 100 ms gap, 0.33 ± 0.03 (n = 8)**, is fully stated and is the
   better holdout. Consider re-targeting the design onto it.

### Acceptance

* Full suite green; `-p no:randomly` used for all comparisons.
* `knowledge/13-preregistration.md` states both the original and re-derived σ, with the
  protocol declared in each.

---

## 6. Phase 5 — Land the `P_o,max = 0.750` verdict

**Depends on:** nothing. **May run in parallel, in its own worktree.** **Size:** small.

Task 4 is **already answered** in `knowledge/14-literature.md` §5 and
`scripts/literature_gamma2.py --section po`. This phase only records it in the code.

### The finding

Five of six candidate values are **intraburst** or **intracluster**. The only genuinely
macroscopic reading is **0.56** (nonstationary fluctuation analysis, Keramidas 2008) — which
is **below** the 0.750 convention, not above it. The handoff's "0.69 by nonstationary
variance analysis" conflates two numbers in the same paper.

### Steps

1. Change the `parameters._NOMINAL_DEFECT` provenance note from "a convention pending
   measurement" to: **a convention that the available measurements cannot replace**, because
   the quantity nearest to it is conditional and the one unconditional reading is 0.56 at an
   incomparable protocol (1–2 ms application; γ2S; errorless; n-less; marked *unpublished
   data* by its own authors; and outside `FIT_RANGES`' [0.70, 0.80] in any case).
2. Keep `Basis.CONVENTION`. **Do not change 0.750 to 0.56.**
3. Update roadmap **P0-13** to say the question is *resolved as unanswerable from the current
   literature*, rather than open pending a search.
4. Point both at `knowledge/14-literature.md` §5 and the `po` section of the script.

### Acceptance

* `pytest -q tests/test_provenance.py tests/test_literature.py -p no:randomly` green.
* `_NOMINAL_DEFECT` is still 0.750 and still `CONVENTION`.

---

## 7. What is still genuinely blocked, and needs a human

No agent can close these; they need institutional access or a judgement call.

| item | what it needs |
|---|---|
| **Jahn 1997 full text** — does it state a mean open time at all? | NeuroReport via Ovid. No PMC copy. **The last open question about the project's own primary source.** |
| **Hamill 1997, NeuroReport 8(16):iv** — one page | Ovid. Highest interest-per-page in the project; bears directly on PR #2. |
| **Barberis Figures 4C and 6C** | Figure digitisation — a *lower* evidence grade than reading a table, and must be graded as such if ever done. |
| **An error statistic for Keramidas 2008 / Dixon 2014** | Neither paper declares one. Probably unobtainable; leave `UNSTATED`. |
| **The new `synthetic`/grade decision in Phase 2** | **Ask the user.** It changes what every `ProvisionalResultWarning` in the project means. |

---

## 8. Summary for the user

| phase | what | depends on | size | parallel? |
|---|---|---|---|---|
| 1 | declare the protocol on deactivation + mean open time | — | S | no |
| 2 | **enter the Barberis deactivation; make `k_off` identifiable** | 1 | M | no |
| 3 | decide `JAHN1997_PEAK_CRC`, now with an absolute rate in hand | 2 | M | no |
| 4 | re-derive the P6 σ requirement with the protocol declared | 1–3 | M/L | no |
| 5 | record the `P_o,max` verdict | — | S | **yes** |

**Two decisions are yours, not an agent's:** the evidence grade in Phase 2, and the
`JAHN1997_PEAK_CRC` option in Phase 3.
