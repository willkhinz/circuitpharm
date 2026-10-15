# Plan — Phase 0 through 4, written for delegated agents

**Author:** Claude Opus 5, 2026-10-14. **Status: NOT STARTED.**
**Supersedes** `PLAN-phase2.md` for everything still open. Read in this order:
`knowledge/14-literature.md` → `knowledge/11-identifiability.md` §1.2.1 → this file.

---

## 0. Where things stand, 2026-10-14

| phase | state |
|---|---|
| **1** protocol declarations | DONE (`b0c547d`) |
| **2** Barberis deactivation, `k_off` | DONE (`b0c547d`, `a18805a`) |
| **5** `P_o,max` verdict | DONE, reviewed and fixed (`3a8f962`) |
| **0** land the uncommitted tree | **NOT STARTED — do this first, see §2** |
| **2b** fill `HOLDOUT` | **NOT STARTED — new, see §4** |
| **3** `JAHN1997_PEAK_CRC` decision | **NOT STARTED — unblocked, see §5** |
| **4** P6 σ re-derivation | **NOT STARTED — see §6** |

Current dataset registry:

| key | observable | kind | role |
|---|---|---|---|
| `dose_response_peak` | PEAK | synthetic | train |
| `deactivation_charge` | CHARGE | synthetic | train |
| `jahn1997_peak_crc` | PEAK | **parametric** | train |
| `barberis2007_deact` | DEACTIVATION | **parametric** | train |

**`HOLDOUT` is empty.** Every dataset is `train`. P5's cross-validation and P6's falsification
bound both consume a holdout, so neither can be scored out-of-sample today.

---

## 1. Rules

Everything in `PLAN-next-upgrade.md` §1 and `PLAN-phase2.md` §1 still applies. **Five
additions, all from defects found during phases 1–2 and 5.**

1. **Nothing in `src/` may import from `scripts/`.** `fitting/data.py` did it to avoid
   retyping the Barberis components — right instinct, wrong direction — and it broke
   `import circuitpharm` from every directory except the repo root. CI never caught it
   because CI runs pytest from the root. Shared data belongs in the package
   (`circuitpharm.literature`); `scripts/` imports *from* it.
2. **Verify with bare `pytest` as well as `python -m pytest`.** The latter prepends cwd to
   `sys.path`; CI runs the former. Treat bare `pytest` as the authority.
3. **Cost-model an objective before profiling it.** A PEAK prediction integrates a waveform
   per concentration and is ~40× dearer than a deactivation score. One profile was sized at
   13 points × 10 restarts × `maxiter=4000` × 17 ms = **2.5 hours**, because
   `xatol=1e-6/fatol=1e-9` made Nelder-Mead exhaust its iteration cap instead of converging.
   Time one evaluation, multiply, then choose.
4. **A profile width over a noise-free parametric dataset is NOT an uncertainty.** Residuals
   are the model's systematic misfit, `sigma_hat` is tiny, curvature is enormous. The `k_off`
   profile's 95% width came out *narrower than its own grid step* and would have been quoted
   as precision. Propagate the **source's stated errors** instead — for `k_off` that was 1.6×
   wide, dominating the profile by more than an order of magnitude. See
   `propagate_published_errors` in `scripts/koff_identifiability.py` for the pattern.
5. **A dataset's SAMPLING GRID is part of its information content.** The Barberis trace
   reproduced its source's parameters to the digit on a uniform 10.1 ms grid and still could
   not constrain the component carrying 57% of its amplitude. **A noiseless refit proves
   nothing** — check recovery at the project's own σ. `scripts/deact_sampling.py` is the
   pattern.

**Worktrees, and why this matters more than it sounds.** One agent per phase:
`git worktree add ../cp-phaseN -b claude/phaseN <base>`. The shared tree has carried ~50
dirty files for four days, which forced **three** separate commits to be assembled
hunk-by-hunk with `git apply --cached` on a filtered patch. Before staging, diff against the
session-start state and commit only your own hunks.

---

## 2. Phase 0 — land the uncommitted tree

**Depends on:** nothing. **Blocks:** everything, in practice. **Size:** S. **Do it first.**

50 files have been dirty since 2026-10-10 19:22. This is **not** a pure lint sweep, which is
why it must be read rather than discarded wholesale:

* **Real work:** `tests/test_review_regressions.py` gains
  `test_deactivation_comparisons_require_identical_protocols`, registered as **E27** — the
  protocol-dependence of deactivation added to the project's recurring-error register. That is
  Phase 1's work and it has never been committed.
* **Genuine lint:** unused imports, `f`-strings with no placeholders, and at least one dead
  assignment (`floors = [r[5] for r in rows]`, assigned and never read — verified).
* A small `tests/test_provisional.py` change.

### Steps

1. Read the whole diff. Do not `git add -A` on faith.
2. Split it: **(a)** the E27 test and any other substantive change, **(b)** the lint. Two
   commits, because one of them is a finding and the other is housekeeping.
3. Confirm the dead-code removals really are dead (`git grep` each removed name).
4. Full suite green under bare `pytest`, then push.

### Acceptance

`git status --short` is empty or near it, and E27 is in the register with its test.

---

## 3. Two decisions that are the user's, not an agent's

Do not implement either unprompted. Both are now *more* pressing than when first raised.

**(a) Should `parametric` be subdivided?** `jahn1997_peak_crc` and `barberis2007_deact` are
both `kind="parametric"` and sit in `TRAIN` together, so the taxonomy cannot express a trust
difference it clearly should: **Jahn extrapolates four decades of curve shape from a slope
measured over one; Barberis extrapolates nothing beyond its measured window** (recomputing
Σ A_i τ_i gives 53.55 ms against a published 52.5 ± 2.9). A `within_measured_window` field
would capture it without multiplying grades. This bears directly on Phase 3.

**(b) Should `po_max`'s basis be `CONVENTION` rather than `FITTED`?** Argued in roadmap
P0-13's OPEN paragraph. `FITTED` means "set to make a *different* anchored quantity come out
right", but `po_max` **is** the anchor. Against the change: `FITTED`'s docstring carries the
"cannot be quoted independently" warning, which is exactly the circularity caveat `po_max`
needs and `CONVENTION` would lose. Regrading moves counts that `tests/test_provenance.py`
defends. **Left as `FITTED` deliberately**, with the argument recorded.

---

## 4. Phase 2b — fill `HOLDOUT`

**Depends on:** Phase 0. **Size:** M. **Highest value per unit of work in this plan.**

`HOLDOUT` is empty, so P5's CV and P6's falsification bound are both unscoreable
out-of-sample. The cheapest honest fill is **already read and recorded**:

> Barberis et al. 2007, Table 1 (source_key `barberis2007`): **paired-pulse recovery at a
> 100 ms gap = 0.33 ± 0.03** (SEM; Fig. 5 caption gives n = 8), with **2 ms pulses of
> saturating (10 mM) GABA**, rat α1β2γ2S, HEK293, 22–24 °C.

It is a **fully stated scalar** — no reconstruction, no extrapolation, no figure
digitisation. That makes it stronger evidence than either `parametric` dataset in `TRAIN`.

### Steps

1. Add `Observable.PAIRED_PULSE_RECOVERY` and a dataset type for a scalar observation with
   **required** `pulse_ms` and `gap_ms`. Follow `MeanOpenTimeDataset`, which is the existing
   scalar-observation precedent, and follow Phase 1's rule: **absence of a declaration is an
   error, not a default.** `MeanOpenTimeDataset.conc_mm` shipped with a silent `10.0` default
   that happened to be right, which is worse than wrong.
2. **Do not add aliased field names.** Phase 1 shipped `pulse_ms`/`pulse_duration_ms` and
   `conc_mm`/`gaba_conc_um`, and contradictory pairs were silently accepted until another
   agent made them raise. One canonical name per quantity.
3. Give the model a prediction path. **`protocols/waveforms.paired_pulse_ratio` already
   exists — reuse it, do not write a second one** (recurring error E12/E22 is exactly that).
   Checked, and the conventions line up, which is the non-obvious part:

   * The existing helper measures amplitude **from the value at pulse onset, not from zero**,
     because at short gaps the channel has not closed and the raw second peak is mostly
     residual from the first.
   * **Barberis defines recovery the same way**: its `R` is built from `I1` (first peak),
     `I2` (second peak) and `I3`, "the current value immediately before the application" —
     i.e. onset-referenced, not zero-referenced.

   So `recovery = peak2 / peak1` measured from zero would be the **wrong** quantity and
   would disagree with the source for reasons that look like a model failure. Its signature
   is `paired_pulse_ratio(res, gaba, *, n_expected)` and it is written for a pulse *train*,
   so confirm it behaves for `n_expected=2` before relying on it.
4. Score it in `likelihood._predict`, with the `q_matrix` guard: a model with no states has
   no paired-pulse recovery. `OperationalScalarModel` has a `simulate_waveform` that relaxes
   toward a Hill target, so **`hasattr(model, "simulate_waveform")` is not a sufficient
   guard** — that mistake cost an hour in Phase 2.
5. Register with `role="holdout"`, `kind="digitised"` (it is a stated value, not a
   reconstruction — but **justify the grade explicitly in `origin`**, since `_validate_dataset`
   requires `synthetic=False` to name a `source_key` and this will be the project's **first
   non-synthetic dataset**; expect `ProvisionalResultWarning` behaviour to change and check
   what that unblocks or breaks).
6. Re-run P5's CV and P6's falisification bound now that a holdout exists. **Report what the
   out-of-sample score actually is**, including if it is bad.

### Acceptance

`HOLDOUT` is non-empty; `compare_models` runs with a real holdout; the full suite is green;
and the commit says plainly whether the models survive out-of-sample scoring.

---

## 5. Phase 3 — decide `JAHN1997_PEAK_CRC`

**Depends on:** Phase 0; **prerequisite already met** (Phase 2 landed an absolute rate).
**Size:** M.

Problem: `HANDOFF.md` §5.1, `knowledge/12-inference.md` §2.4. Nine points of a global Hill
curve, nH = 2.2, spanning 0.3–3000 µM — **four decades extrapolated from a slope measured
over one** (1–10 µM). The source constrains three things: the slope over 1–10 µM, the EC₅₀,
and saturation by 3 mM.

### What has changed since the question was posed

* **It does not dissolve.** No digitisable CRC exists anywhere in the corpus of seven papers
  read; Jahn's remains the only one. That was `HANDOFF.md`'s stated condition for dissolution
  and it was not met.
* **An absolute rate is now available.** `k_off` is constrained by the Barberis trace, so the
  three options can be fitted *with* absolute-rate information — which is the whole reason
  this phase was moved after Phase 2. The `D`-on-its-bound / 0.99-plateau artefact may be a
  symptom of having had none.
* **An independent EC₅₀ exists:** Mortensen 2012 gives α1β2γ2**S** EC₅₀ = **6.6 µM**
  (pEC₅₀ 5.180 ± 0.0593, n = 34), HEK293 whole-cell, against Jahn's 11.6 ± 0.9 µM. **These
  are not the same observable** — Mortensen delivers GABA with a 20–30 ms latency, so it is
  not a fast-application peak, and 6.6 < 11.6 is the direction this project's own PEAK
  protocol-dependence table predicts. Use it as a **consistency check on the EC₅₀ anchor,
  never as a second point to pool.** Mortensen reports **no Hill coefficient**, so it cannot
  constrain shape.

### Steps

1. Implement the three `HANDOFF.md` §5.1 options behind a flag: keep as is; restrict to the
   measured window plus a saturation anchor; split in two (a defensible restricted dataset
   plus the full reconstruction kept to *demonstrate* the artefact).
2. Fit under each **with `barberis2007_deact` included**. Print: does `D` still land on its
   bound? Is the plateau still 0.99? Does the restricted ~3-point version stay uninformative
   now that an absolute rate is available?
3. **Cost-model first** (rule 3). Three fits against a PEAK dataset is affordable; a profile
   per option is not. If you need intervals, propagate source errors (rule 4) rather than
   widening the profile budget.
4. Recommend **with printed evidence**, and state the cost plainly. If the honest answer is
   "most downstream fits become uninformative", say so — `HANDOFF.md` explicitly allows that
   this may be the correct state.
5. **Do not pick.** Bring the comparison to the user before touching `TRAIN`, the P5
   comparison or the P6 design. Flag §3(a) at the same time: the right option may depend on
   whether `parametric` gets subdivided.

---

## 6. Phase 4 — re-derive the P6 σ requirement with the protocol declared

**Depends on:** Phases 0, 2b, 3. **Size:** M–L (MCMC). **Risk:** low.

`knowledge/13-preregistration.md` concludes the discriminating experiment does not exist at
σ = 0.02 and needs **σ ≤ 0.005** — a 4× noise reduction. That rests on posteriors fitted
partly through `PEAK`, which was later found protocol-dependent: application length alone
moves peak EC₅₀ **sevenfold** (299 µM at 2 ms → 39 µM at 1000 ms). `PEAK_APPLICATION_MS` is
fixed at 300 ms and no dataset declares its own.

### Steps

1. Re-run the derivation with the application duration **declared and matched** between model
   and dataset. Use `models.base.hill_slope(..., window_rel=...)` for any slope comparison.
2. Report whether σ ≤ 0.005 **survives**. It may soften: some apparent indistinguishability
   may be protocol smearing rather than a property of the two schemes.
3. If it softens, **do not quietly relax the pre-registration.** A pre-registered conclusion
   that changes *is* the result — record old and new side by side with the protocol stated
   for each.
4. **Re-target the discriminating observable.** P6 chose Barberis's desensitisation onset, but
   `knowledge/14-literature.md` §4 shows it **cannot be fully reconstructed** (τ₂/τ₃ appear
   only in Figure 4C) and the source is **internally inconsistent about its n** (text says 7,
   the Fig. 4 caption says 3). The paired-pulse recovery from Phase 2b is fully stated — prefer
   it, and declare its pulse duration and gap.
5. Apple Silicon: a marginal split-R̂ failure in `test_posterior_recovers_known_parameters` is
   a known arm64-vs-x86 difference through the matrix exponentials. **Do not widen the gate.**

---

## 7. Blocked on a human, not an agent

| item | needs |
|---|---|
| **Jahn 1997 full text** — does it state a mean open time at all? | NeuroReport via Ovid. No PMC copy. The last open question about the project's own primary source. |
| **Hamill 1997, NeuroReport 8(16):iv** — one page | Ovid. Highest interest-per-page in the project; bears on PR #2's withdrawal. |
| **Barberis Figures 4C / 6C** | Figure digitisation — a *lower* grade than reading a table, and must be graded as such. |
| **A mean open time for α1β2γ2L** | Does not exist in the corpus. Four readings span 5.1× and the spread is *definition*, not disagreement. `FIT_FIXED_ALPHA` stays a `CONVENTION` until one is measured at a declared conditioning. **Still the highest-value acquisition in the project.** |
| §3(a) subdivide `parametric`? · §3(b) regrade `po_max`? | User. |

---

## 8. Summary

| # | phase | depends on | size |
|---|---|---|---|
| 0 | **land the uncommitted tree** (contains E27, not just lint) | — | S |
| 2b | **fill `HOLDOUT`** with Barberis paired-pulse recovery | 0 | M |
| 3 | decide `JAHN1997_PEAK_CRC`, with an absolute rate in hand | 0 | M |
| 4 | re-derive P6 σ with the protocol declared | 0, 2b, 3 | M/L |

Phase 0 is small and unblocks the rest in practice. **Phase 2b is the highest value per unit
of work**: it is the only step that adds a dataset needing no reconstruction, and it turns two
existing but unscoreable machines (P5's CV, P6's bound) on.
