# Plan — Phase 1b through 4, written for delegated agents

**Author:** Claude Opus 5, 2026-10-10, after reviewing Phase 5 (`3a8f962`).
**Status: NOT STARTED. Awaiting go-ahead.**
**Supersedes** `PLAN-next-upgrade.md` §2–§5. Read `knowledge/14-literature.md` first.

---

## 0. Where things actually stand

| phase | state |
|---|---|
| **1** — protocol declarations | **implemented, UNCOMMITTED, and partly unlanded.** `DeactivationDataset`, `Observable`, `ObservableMismatch` are in HEAD; `MeanOpenTimeDataset` and `check_deactivation_protocols` are **not**. Two real defects (§2). |
| **5** — `P_o,max` verdict | **done, reviewed, fixed, pushed** (`3a8f962`). Four defects found and corrected; see that commit. |
| **2, 3, 4** | not started. |

**One decision I flagged for the user has dissolved.** `PLAN-next-upgrade.md` §3 step 2 said a
third evidence grade might need inventing for "a curve reconstructed from source-stated
parameters", and told the agent to ask first. **It already exists.** `fitting.data.DataKind`
is `Literal["digitised", "parametric", "synthetic"]`, and `parametric` is documented as:

> generated from the FITTED PARAMETERS a paper published (EC50, Hill slope, Emax). The
> parameters are real measurements and are attributable; the individual points are the
> paper's own model of its data, not its data. Weaker than a digitisation and far stronger
> than a guess.

That is exactly the category. **Use `kind="parametric"`. Do not invent a grade, and do not
ask about it.** Note what follows from it: `_validate_dataset` requires
`parametric ⇒ synthetic=True`, so **entering Barberis does NOT stop
`ProvisionalResultWarning`** and does not make any result VALIDATED. It makes `k_off`
identifiable, which is a different and still worthwhile thing. Say so plainly rather than
implying the warning should go away.

---

## 1. Rules

Everything in `PLAN-next-upgrade.md` §1 still applies. **Six additions, every one of them
from a defect found in the last two days — four of them in work done by agents following
the previous version of this plan.**

1. **Verify with bare `pytest`, not `python -m pytest`.** `python -m pytest` prepends the
   current directory to `sys.path`; bare `pytest` does not, and **CI runs bare `pytest`.**
   `tests/test_literature.py` imports from `scripts/` and so passed in two local sessions and
   died in CI with `ModuleNotFoundError: No module named 'scripts'`. Fixed in `5643f74` by
   `pythonpath = ["."]`. My own fault: I verified exclusively with `python -m pytest`. **Run
   both, and treat bare `pytest` as the authority.**
2. **Never add a constant that duplicates a value which already has a home — derive it.**
   Phase 5 added `parameters.NOMINAL_DEFECT = 0.750` and `PO_MAX_CONVENTION = 0.750` as
   literals, when `gabaa_kinetics.FIT_TARGETS["po_max"]` already existed *and was already
   imported in the same file*. Because they were literals, the canonical value could have
   moved and the test would still have passed. `dynamic_range.py` warns about this in prose;
   it is recurring error E12/E22.
3. **Never give one quantity two field names.** Phase 1 added `pulse_ms` *and*
   `pulse_duration_ms`, `conc_mm` *and* `gaba_conc_um`, cross-filled only when one is
   `None` — so contradictory pairs are **silently accepted** (§2.1). E22 reproduced inside
   the module built to prevent E22.
4. **A test must check the invariant its name claims.** Phase 5's
   `test_nominal_defect_stays_convention_and_unreplaced` never looked at a `Basis` — and had
   it done so faithfully it would have failed, because the code records `FITTED` while the
   roadmap claimed `CONVENTION`.
5. **Absence of a declaration is an error, not a default.** This was step 4 of Phase 1 and
   Phase 1 violated it: `MeanOpenTimeDataset.conc_mm` defaults to `10.0` (§2.2). The default
   happens to be *correct* for Barberis, which is worse than wrong — a wrong one would have
   been noticed.
6. **Do not assert in prose what the code does not say.** If a markdown file makes a claim
   about a code value, a test must pin the two together. Phase 5's roadmap paragraph asserted
   a `Basis` the code did not have, and nothing compared them.

**Worktrees.** One agent per phase: `git worktree add ../cp-phaseN -b claude/phaseN <base>`.
Phases are sequential. Before staging, diff against the session-start state and commit **only
your own hunks** — the tree routinely carries 50+ files of other agents' in-flight work, and
two commits in the last two days had to be assembled hunk-by-hunk with
`git apply --cached` on a filtered patch because of it.

---

## 2. Phase 1b — fix Phase 1's two defects, then land it

**Size:** small. **Blocks:** Phase 2. **Do this first; it is the prerequisite.**

### 2.1 Contradictory aliases are accepted

Reproduce before fixing:

```python
d = DeactivationDataset(
    citation_label="contradiction", time_ms=t, normalized_current=y,
    observable=Observable.DEACTIVATION,
    pulse_ms=1.0, pulse_duration_ms=999.0,    # one quantity, two values
    conc_mm=3.0,  gaba_conc_um=1.0,           # one quantity, two values (3000x apart)
    n_components=2, fit_window_ms=50.0,
    synthetic=True, kind="synthetic", origin="probe")
# ACCEPTED. d.pulse_ms == 1.0 and d.pulse_duration_ms == 999.0
```

**Preferred fix:** one canonical field per quantity — keep `pulse_ms` and `conc_mm`, delete
`pulse_duration_ms` and `gaba_conc_um`, and update the callers (`DEACTIVATION_BENCHMARK` sets
all four). **If the aliases must stay** for compatibility, `__post_init__` must *raise* on a
contradiction rather than cross-fill only when one side is `None`. Add a test asserting the
contradiction raises — the code above is the test case.

### 2.2 `MeanOpenTimeDataset.conc_mm` has a silent default

`conc_mm: float = 10.0`. Make it required, or default to `None` and raise in
`_validate_dataset` when the observable is `MEAN_OPEN_TIME`. Follow how Phase 1 already made
`pulse_ms` required for `Observable.DEACTIVATION` — that part was done right.

### 2.3 Not a defect, do not "fix" it

`kind="digitised"`/`synthetic=False` are the **dataclass defaults**, which looks backwards for
this project. But the pre-existing `DoseResponseDataset` has exactly the same defaults, so it
is codebase-wide convention and **out of scope here.** Changing it touches every dataset.
Raise it as a separate question if you think it matters; do not change it inside this phase.

### 2.4 Then land it

Phase 1's implementation and its test in `tests/test_literature.py` are uncommitted. Commit
them **together** — the test references `MeanOpenTimeDataset` and
`check_deactivation_protocols`, so a commit with one and not the other produces a tree whose
tests reference code it does not contain. Verify by building the staged tree
(`git archive $(git write-tree) | tar -x -C <tmp>`) and running it there.

### Acceptance

* Bare `pytest -q -m "not slow" -p no:randomly` green, **and** `python -m pytest` green.
* A test proves the contradictory-alias case raises.
* A test proves a `MEAN_OPEN_TIME` dataset with no stated concentration raises.
* `git grep pulse_duration_ms` and `git grep gaba_conc_um` return either nothing or only
  code that raises on contradiction.

---

## 3. Phase 2 — enter the Barberis deactivation; make `k_off` identifiable

**Depends on:** 1b. **Size:** medium. **This is the highest-value step in the plan.**

### The data — import it, do not retype it

Everything is in `scripts/literature_gamma2.py`. `DEACTIVATIONS[1]` is the row
(`label="Barberis 2007, brief pulse"`), with `.components == ((2.8, 33.4, 221.35),
(0.57, 0.23, 0.20))`. **Import those; a literal copy here would be the fourth copy of a
number this project already owns.**

Barberis et al. 2007, Eur J Neurosci 25(9):2726–2740, PMC1950087, source_key
**`barberis2007`** (already in the `sources` table with a `HAND:` verdict):

* rat **α1β2γ2S** — the Methods say γ2S though the title says α1β2γ2
* HEK293 outside-out pooled with small lifted whole cells, **22–24 °C**, −70 mV
* ultrafast exchange, 10–90% in 60–100 µs
* **2 ms pulse of 10 mM GABA**
* τ = 2.8 ± 0.3 / 33.4 ± 4.6 / 221.35 ± 14.9 ms; areas 0.57 ± 0.04 / 0.23 ± 0.02 / 0.20 ± 0.03
* SEM, **n = 6**; Σ A_i = 1.000; τ_w = 52.5 ± 2.9 ms

### Steps

1. Build the trace as Σ A_i·exp(−t/τ_i) from the imported parameters, over a declared
   `fit_window_ms`. **Choose the window deliberately and say why**: the slowest component is
   221 ms, so a window shorter than ~1 s truncates the component that carries 44.3 of the
   52.5 ms weighted τ — which is precisely the mechanism behind the 8.9× Dixon/Barberis gap
   (`knowledge/14-literature.md` §3). Getting this wrong reproduces the error the dataset
   exists to document.
2. Declare it:
   ```
   kind="parametric", synthetic=True, source_key="barberis2007",
   observable=Observable.DEACTIVATION,
   pulse_ms=2.0, conc_mm=10.0, n_components=3, fit_window_ms=<chosen>,
   normalisation=<see below>, n_cells=6, origin="<which published parameters, and how>"
   ```
   `origin` must name the three τ, the three areas, and the pulse — `_validate_dataset`
   requires a parametric dataset to state "which published parameters were used and how the
   curve was generated", and that check will reject a vague string.
3. **Decide `normalisation` explicitly.** Σ A_i = 1 means the trace starts at 1.0, so this is
   `"fraction_of_max"`, not `"absolute"`. Getting it wrong "does not bias a fit, it breaks
   it" — the module says so, and `Normalisation` exists because an optimiser once closed a
   0.25 plateau gap by deleting desensitisation.
4. **Record the γ2S substitution in the dataset itself**, not only in prose. The project's own
   receptor is γ2L and only two sources in the corpus are γ2L (Li 2008, Dixon 2014) — neither
   of which gives a usable deactivation time course.
5. Register in `provenance.py`. Add to `ALL`; decide `role="train"` vs `"holdout"` and say
   why. Note `HOLDOUT` is currently **empty**, so if this becomes the holdout, P5's
   cross-validation and P6's falsification bound get their first out-of-sample score — but
   then it cannot also make `k_off` identifiable in the fit. **That trade is the real
   decision in this phase.** Report it rather than resolving it silently; the honest default
   is `train`, with the Barberis paired-pulse recovery (0.33 ± 0.03 at a 100 ms gap) as the
   holdout instead.
6. **Re-run the identifiability analysis** (`knowledge/11-identifiability.md`,
   `fitting/identifiability.py`): does adding this make `k_off` identifiable? **Report the
   answer either way.** A negative result is a real finding — it would mean the limit is the
   parametric reconstruction rather than the absence of a deactivation observable — and must
   not be buried.
7. Print the before/after from a script. Update `knowledge/11-identifiability.md` with what it
   actually printed, not with what was expected.

### Acceptance

* Full suite green under **bare** `pytest`.
* A script prints the `k_off` identifiability result before and after.
* `ProvisionalResultWarning` still fires (parametric ⇒ synthetic), and the commit says so
  rather than implying progress toward VALIDATED that did not happen.
* `check_deactivation_protocols(BARBERIS, DEACTIVATION_BENCHMARK)` raises — different pulse,
  concentration, component count and window.

---

## 4. Phase 3 — decide `JAHN1997_PEAK_CRC`

**Depends on:** Phase 2. **Do not start before it** — the whole reason for the reordering is
that this decision should be taken *with* an absolute rate available. **Size:** medium.

Problem statement: `HANDOFF.md` §5.1 and `knowledge/12-inference.md` §2.4. Nine points of a
global Hill curve, nH = 2.2, spanning 0.3–3000 µM — four decades of shape extrapolated from a
slope measured over one (1–10 µM). The source constrains three things: the slope over
1–10 µM, the EC₅₀, and saturation by 3 mM.

### What the literature pass changed

* **No digitisable CRC exists in the corpus.** Jahn's remains the only one. The question does
  **not** dissolve — that was the condition `HANDOFF.md` set, and it was not met.
* **An independent EC₅₀ exists:** Mortensen 2012, α1β2γ2**S**, **6.6 µM** (pEC₅₀
  5.180 ± 0.0593, n = 34), HEK293 whole-cell — against Jahn's 11.6 ± 0.9 µM. **Not the same
  observable:** Mortensen delivers GABA with a 20–30 ms latency, so it is not a
  fast-application peak, and the direction (6.6 < 11.6) is what this project's own PEAK
  protocol-dependence table predicts for a longer effective application. Use it as a
  **consistency check on the EC₅₀ anchor, never as a second point to pool.**
* Mortensen reports **no Hill coefficient**, so it cannot constrain curve shape.
* **Both Jahn and Barberis are `kind="parametric"`** — but they are not equally so, and the
  taxonomy does not currently distinguish them: Jahn extrapolates four decades from one,
  Barberis extrapolates nothing beyond its measured window. Whether `parametric` should be
  subdivided, or a `within_measured_window` flag added, is a **question for the user** — it is
  the sharpened remnant of the grade question that otherwise dissolved. Raise it; do not
  implement it unprompted.

### Steps

1. Implement the three `HANDOFF.md` §5.1 options behind a flag: keep as is; restrict to the
   measured window plus a saturation anchor; split in two (defensible restricted dataset, plus
   the full reconstruction retained to *demonstrate* the artefact).
2. Fit under each **with the Phase 2 deactivation included**, and print: does `D` still land
   on its bound? Is the plateau still 0.99? Does the restricted ~3-point version stay
   uninformative once an absolute rate is available?
3. Recommend, **with the printed evidence**, and state what it costs. If the honest answer is
   "most downstream fits become uninformative", say it — `HANDOFF.md` explicitly allows that
   this may be the correct state.
4. **Do not pick silently.** Bring the comparison to the user before changing `TRAIN`, the P5
   comparison or the P6 design.

---

## 5. Phase 4 — re-derive the P6 σ requirement with the protocol declared

**Depends on:** 1b–3. **Size:** medium–large (MCMC). **Risk:** low.

`knowledge/13-preregistration.md` concludes the discriminating experiment does not exist at
σ = 0.02 and needs **σ ≤ 0.005**, a 4× noise reduction. That rests on posteriors fitted partly
through `PEAK`, which was later found protocol-dependent: application length alone moves peak
EC₅₀ **sevenfold** (299 µM at 2 ms → 39 µM at 1000 ms). `PEAK_APPLICATION_MS` is fixed at
300 ms and no dataset declares its own.

### Steps

1. Re-run the derivation with the application duration **declared and matched** between model
   and dataset. Use `models.base.hill_slope(..., window_rel=...)` for any slope comparison.
2. Report whether σ ≤ 0.005 **survives**. It may soften — some of the apparent
   indistinguishability may be protocol smearing rather than a property of the two schemes.
3. If it softens, **do not quietly relax the pre-registration.** A pre-registered conclusion
   that changes *is* the result: record old and new side by side with the protocol stated for
   each.
4. **Re-target the discriminating observable.** P6 chose Barberis's desensitisation onset, but
   `knowledge/14-literature.md` §4 shows it **cannot be fully reconstructed** (τ₂/τ₃ only in
   Figure 4C) and its n is ambiguous in the source (text 7, Fig. 4 caption 3). The
   **paired-pulse recovery at a 100 ms gap, 0.33 ± 0.03 (n = 8)** is fully stated. Prefer it,
   and declare its pulse duration (2 ms) and gap.
5. Apple Silicon: a marginal split-R̂ failure in `test_posterior_recovers_known_parameters` is
   a known arm64-vs-x86 difference through the matrix exponentials. **Do not widen the gate.**

---

## 6. Still blocked on a human

| item | needs |
|---|---|
| **Jahn 1997 full text** — does it state a mean open time at all? | NeuroReport via Ovid. No PMC copy. The last open question about the project's own primary source. |
| **Hamill 1997, NeuroReport 8(16):iv** — one page | Ovid. Highest interest-per-page in the project; bears on PR #2. |
| **Barberis Figures 4C / 6C** | Figure digitisation — a *lower* grade than reading a table, and must be graded as such. |
| **Should `parametric` be subdivided?** (§4) | User. Jahn and Barberis are both `parametric` and should not be equally trusted. |
| **Should `po_max`'s basis be `CONVENTION` rather than `FITTED`?** | User. Argued in roadmap P0-13's OPEN paragraph; regrading moves provenance counts that `tests/test_provenance.py` defends. |
| **Train or holdout for Barberis?** (§3 step 5) | User, after seeing step 6's identifiability result. |

---

## 7. Summary

| phase | what | depends on | size |
|---|---|---|---|
| **1b** | fix Phase 1's aliases + silent default, then land Phase 1 | — | S |
| **2** | **Barberis deactivation; `k_off` identifiable** | 1b | M |
| **3** | decide `JAHN1997_PEAK_CRC`, with an absolute rate in hand | 2 | M |
| **4** | re-derive P6 σ with the protocol declared | 1b–3 | M/L |

Phase 5 is done. The grade question that previously blocked Phase 2 is answered by
`kind="parametric"`, which already exists.
