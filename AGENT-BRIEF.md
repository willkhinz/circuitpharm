# Agent brief — four tasks, in order

**Written by:** Claude Opus 5, 2026-10-14. **Self-contained: assume you know nothing.**
One agent, one worktree, four tasks, A → D. Each has its own commit.

```bash
git worktree add ../cp-brief -b claude/agent-brief claude/quirky-sagan-edtdet
cd ../cp-brief
python3.11 -m venv .venv && source .venv/bin/activate   # 3.13 does NOT work
pip install -e ".[dev]"
```

Read before you start: `knowledge/14-literature.md`, then `knowledge/11-identifiability.md`
§1.2.1. `PLAN-phase3.md` is the wider plan; this brief is tasks 0, (a), (b) and 2b of it,
with the two decisions already made for you.

---

## Rules — every one of these is from a real defect in this branch

1. **Never skip, disable, `xfail` or quarantine a test.** If one fails, either the code is
   wrong or the invariant changed — and if the latter, rewrite the invariant *and say so in
   the commit*.
2. **Verify with bare `pytest`, not just `python -m pytest`.** The latter prepends cwd to
   `sys.path`; **CI runs bare `pytest`.** A test here once passed in two local sessions and
   died in CI for exactly this reason.
3. **Nothing in `src/` may import from `scripts/`.** It broke `import circuitpharm` from
   every directory except the repo root, invisibly to CI. Shared data lives in
   `circuitpharm.literature`; `scripts/` imports *from* the package.
4. **Numbers are printed by code, never typed into prose.** Any figure in a markdown file
   must be emitted by a script in `scripts/`.
5. **One canonical name per quantity.** Do not add aliased fields. `pulse_ms` *and*
   `pulse_duration_ms` shipped here once and silently accepted contradictory pairs.
6. **Absence of a declaration is an error, not a default.** `MeanOpenTimeDataset.conc_mm`
   shipped with a silent `10.0` that happened to be *correct*, which is worse than wrong.
7. **A test must check the invariant its name claims.**
8. **Record your own errors in the commit message** rather than quietly fixing them. This is
   load-bearing convention here.
9. **Commit only your own hunks.** The shared tree carries ~50 unrelated dirty files; three
   commits on this branch had to be assembled with `git apply --cached` on a filtered patch
   because of it. Task A clears this, so do it first.
10. **Cost-model before you profile.** Time one objective evaluation, multiply by
    grid × restarts × `maxiter`, then choose. A profile here was once sized at 2.5 hours
    because tight Nelder-Mead tolerances made it exhaust `maxiter` instead of converging.

```bash
pytest -q -m "not slow" -p no:randomly   # ~75 s, 496 tests
pytest -q -p no:randomly                 # full, ~11 min, 561 tests — before every push
```

`-p no:randomly` whenever comparing runs. On Apple Silicon a marginal split-R̂ failure in
`test_posterior_recovers_known_parameters` is a known arm64-vs-x86 difference through the
matrix exponentials — **do not widen that gate.** `ProvisionalResultWarning` on a
`synthetic=True` dataset is the machinery working.

---

## Task A — land the uncommitted tree

**Size:** S. **Do first; it unblocks the rest in practice.**

~50 files have been dirty since 2026-10-10 19:22. **It is not a pure lint sweep**, so read it
rather than discarding or bulk-adding it.

* **Substantive:** `tests/test_review_regressions.py` gains
  `test_deactivation_comparisons_require_identical_protocols`, registered as **E27** — the
  protocol-dependence of deactivation entered in the project's recurring-error register. That
  is real work and has never been committed.
* **Lint:** unused imports, `f`-strings with no placeholders, and a dead
  `floors = [r[5] for r in rows]` (verified assigned and never read).
* A small `tests/test_provisional.py` change.

### Steps

1. Read the whole diff. **Do not `git add -A` on faith.**
2. Two commits: **(1)** E27 and any other substantive change, **(2)** the lint. One half is a
   finding, the other is housekeeping, and they should not share a message.
3. `git grep` every removed name to confirm the dead-code removals are dead.
4. Full suite green under bare `pytest`, then push.

**Done when** `git status --short` is empty or near it, and E27 is in the register.

---

## Task B — make the Jahn extrapolation machine-visible

**Size:** S–M. **This is a decision already taken; implement it, do not re-litigate it.**

### The decision, and why

`jahn1997_peak_crc` and `barberis2007_deact` are both `kind="parametric"` and both in `TRAIN`,
so the taxonomy cannot express a difference it should:

* **Jahn** spans **0.3 → 3000 µM** (9 points) while its source measured the Hill slope over
  **1–10 µM only**. Six of its nine points lie outside the window the source measured.
* **Barberis** spans **0 → 1000 ms** and its source fitted the decay over that window.
  **Nothing is outside.**

**DECIDED: do NOT subdivide `DataKind`.** `DataKind` records *how the points came to be*
(digitised / parametric / synthetic), and both of these came to be the same way — generated
from published fitted parameters. Splitting the enum would conflate **provenance** with
**coverage**, and would force a regrade of existing datasets, change what
`ProvisionalResultWarning` means, and touch `_validate_dataset`'s `kind`↔`synthetic` coupling
for no gain.

**DECIDED: add a coverage field instead.** Coverage is a per-dataset fact with a crisp
definition, and making it a field rather than a comment is this project's standing rule — *a
note in prose does not constrain a computation; only a parameter does.*

### Steps

1. Add to both dataset classes:

   ```python
   #: For kind="parametric" ONLY: the domain the SOURCE actually measured, in this dataset's
   #: own x units and as (lo, hi). REQUIRED for parametric data, because the entire
   #: difference between a faithful reconstruction and an extrapolation is whether the
   #: dataset's own points lie inside it.
   source_measured_domain: tuple[float, float] | None = None
   ```

2. Add a derived, axis-agnostic property. **Do not use decades** — one axis is concentration
   (log-natural) and the other is time starting at 0, so a decade measure breaks on the
   second. Count points instead:

   ```python
   @property
   def n_points_outside_source_domain(self) -> int: ...
   @property
   def fraction_outside_source_domain(self) -> float: ...
   ```

   Use `concs_um` for a `DoseResponseDataset` and `time_ms` for a `DeactivationDataset`.
3. In `_validate_dataset`, require `source_measured_domain` when `kind == "parametric"`, with
   a message explaining why (rule 6 — it must be an error, not a default). Leave
   `digitised` and `synthetic` unaffected.
4. Declare the two existing values:
   * `jahn1997_peak_crc`: `source_measured_domain=(1.0, 10.0)` — the 1–10 µM slope window.
     In `origin`, note that the EC₅₀ (11.6 µM) and saturation (3 mM) are *separately*
     constrained by the source even though they sit outside that window, so the field records
     the **slope** window and is not the whole of what the source says.
   * `barberis2007_deact`: `source_measured_domain=(0.0, 1000.0)`.
5. Tests: Jahn's `n_points_outside_source_domain` is **6 of 9**; Barberis's is **0**; a
   parametric dataset built without the field **raises**. Assert the asymmetry directly —
   that is the finding.
6. Print it. Add a section to `scripts/literature_gamma2.py` (or a short new script) that
   tabulates every registered dataset with its kind, its source domain and its fraction
   outside, so the Jahn artefact is visible from a command rather than from prose. Update
   `knowledge/12-inference.md` §2.4 and `MISSING_DATASETS` to point at it.

**Done when** `fraction_outside_source_domain` separates the two datasets, a parametric
dataset cannot be constructed without declaring its source domain, and a script prints the
table.

---

## Task C — close the `po_max` basis question

**Size:** XS. **Also already decided.**

Roadmap **P0-13** carries an OPEN paragraph asking whether
`provenance.MODEL_PARAM_PROV[("operational", "po_max")]` should be `Basis.CONVENTION` rather
than `Basis.FITTED`.

**DECIDED: keep `Basis.FITTED`. Close the question.** The reasoning, which belongs in the
document:

* `Basis` records **how a number was derived**, not whether it *could be replaced*. The
  2026-10-10 literature audit was a finding about replaceability, and it already lives in
  `_NOMINAL_DEFECT["operational"]` where a reader of the defaults meets it.
* `Basis.FITTED`'s own definition carries the warning this number most needs: *"it means the
  number carries that anchor's assumptions and **cannot be quoted independently**."*
  `po_max = 0.75` **is** the fit target, and `provenance.py` already calls it "circular if
  used as an anchor". `Basis.CONVENTION` — "a structural choice, not an empirical quantity" —
  would **lose** that circularity warning, which is the opposite of the intent.
* Regrading moves counts that `tests/test_provenance.py` defends, for no informational gain.

### Steps

1. Rewrite P0-13's **OPEN** paragraph as **DECIDED**, with the three reasons above.
2. Do **not** change any code. `tests/test_literature.py::`
   `test_po_max_convention_matches_what_the_code_records` already pins the basis to `FITTED`
   and tells a future regrader to update the roadmap in the same commit — check that its
   message still reads correctly against the new wording.

**Done when** P0-13 records a decision with its reasoning and no code changed.

---

## Task D — fill `HOLDOUT`

**Size:** M. **The payoff. Highest value per unit of work in the plan.**

`HOLDOUT` is **empty** and every dataset is `train`. P5's cross-validation and P6's
falsification bound both consume a holdout, so two machines that exist cannot be run.

### The data — already read and recorded

Barberis et al. 2007, *Eur J Neurosci* 25(9):2726–2740, PMC1950087, `source_key`
**`barberis2007`** (already in the `sources` table with a `HAND:` verdict):

> **Paired-pulse recovery at a 100 ms gap = 0.33 ± 0.03** (SEM; Fig. 5 caption gives n = 8),
> **2 ms pulses** of saturating (**10 mM**) GABA, rat **α1β2γ2S**, HEK293 outside-out pooled
> with small lifted whole cells, **22–24 °C**, −70 mV, ultrafast exchange 60–100 µs.

It is a **fully stated scalar** — no reconstruction, no extrapolation, no figure
digitisation. That makes it stronger evidence than either `parametric` dataset in `TRAIN`, and
it will be the project's **first non-synthetic dataset**.

### The trap, checked for you

`protocols/waveforms.paired_pulse_ratio(res, gaba, *, n_expected)` **already exists — reuse
it, do not write a second one** (recurring error E12/E22 is exactly that). The conventions
line up, and this is the non-obvious part:

* The existing helper measures amplitude **from the value at pulse onset, not from zero**,
  because at short gaps the channel has not closed and the raw second peak is mostly residual
  from the first.
* **Barberis defines recovery the same way**: its `R` is built from `I1` (first peak), `I2`
  (second peak) and `I3`, *"the current value immediately before the application"* —
  onset-referenced, not zero-referenced.

So `recovery = peak2 / peak1` measured from zero is the **wrong quantity** and would disagree
with the source for reasons that look like a model failure. The helper is written for a pulse
*train*, so confirm it behaves at `n_expected=2` before relying on it.

### Steps

1. Add `Observable.PAIRED_PULSE_RECOVERY` and a scalar dataset type with **required**
   `pulse_ms`, `conc_mm` and `gap_ms`. Follow `MeanOpenTimeDataset` (the existing
   scalar-observation precedent) and rules 5 and 6: one canonical name per quantity, and no
   silent defaults.
2. Prediction path in `likelihood._predict`, guarded on **`q_matrix`, not
   `simulate_waveform`** — every model here has a `simulate_waveform`, and
   `OperationalScalarModel`'s relaxes toward a Hill target at an assumed time constant, so a
   `hasattr(model, "simulate_waveform")` guard would score that assumption as if it were a
   rate. That mistake cost an hour in Phase 2; the message should say the model has no states
   to relax between.
3. Register it with `role="holdout"`. **Justify the grade explicitly in `origin`.** This is
   the first dataset that is not `synthetic=True`, so `_validate_dataset` will require a
   `source_key`, and `ProvisionalResultWarning` behaviour will change — **find out what that
   unblocks or breaks and report it**, do not just make it green.
4. Record the **γ2S substitution** in the dataset, not only in prose: the project's receptor
   is γ2L.
5. Re-run P5's cross-validation and P6's falsification bound now that a holdout exists.
   **Report what the out-of-sample score actually is, including if it is bad.** A negative
   result here is a real finding — it would say the models do not generalise — and must not be
   buried. Update `knowledge/11-identifiability.md` and the roadmap with what was printed.

**Done when** `HOLDOUT` is non-empty, `compare_models` runs against a real holdout, the full
suite is green under bare `pytest`, and the commit states plainly whether the models survive
out-of-sample scoring.

---

## Do not do these

* **Do not decide `JAHN1997_PEAK_CRC`** (Phase 3 of `PLAN-phase3.md`). Task B makes its
  artefact measurable, which is its prerequisite; the choice between keeping, restricting and
  splitting the dataset is the user's and needs the three-way fit comparison first.
* **Do not re-derive the P6 σ requirement** (Phase 4). It depends on Task D.
* **Do not chase the blocked literature.** The Jahn 1997 full text and the one-page Hamill
  comment (NeuroReport 8(16):iv) need institutional Ovid access, not more searching. Barberis
  Figures 4C/6C would need figure digitisation, which is a *lower* evidence grade than
  reading a table and must be graded as such if ever done.
* **Do not invent a mean open time for α1β2γ2L.** No source states one. Four readings exist,
  they span 5.1×, and the spread is *definition* (intraburst vs intracluster vs a
  non-stationary window) rather than disagreement. `FIT_FIXED_ALPHA` stays a `CONVENTION`.
