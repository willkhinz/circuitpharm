# Handoff — state of play, 2026-10-11

> **Amended 2026-10-11 after a local session read the literature.** The egress block was
> environment-specific — PMC is reachable from outside the cloud container, and all seven
> open-access targets were read in full. §4 and §5.3 are corrected below, including one
> substantive error of mine (the "0.69 macroscopic P_o" was a conflation) and four wrong
> citations. **`knowledge/14-literature.md` and `scripts/literature_gamma2.py` now
> supersede §4 of this document.**

Written for someone picking this up in a fresh local session. It says what was just done,
what is true now, what is blocked, and what I would do next and in what order.

---

## 0. Thirty-second version

* **PR #1 is MERGED.** The whole P0–P7 roadmap is on `main` (`6f9fcfa`). CI green: 499
  passed, 0 skipped, on lean 3.11, lean 3.12 and full.
* **PR #2 is OPEN as a draft** on branch `claude/quirky-sagan-edtdet`, head `a8d9019`. It
  **withdraws** what P4 had recorded as its main result. CI was still running when this was
  written — check it before trusting anything below about test counts.
* **The repo is now PUBLIC.** That was the fix for CI dying in 3 seconds since 2026-10-09:
  it was a private-repo Actions minutes/billing refusal, not a code problem.
* **The project's blocker is no longer scientific. It is an egress allowlist.** See §4.

---

## 1. What PR #2 changes, and why it matters more than it looks

P4 recorded, as its headline finding, that the published Hill slope and this project's
assumed plateau **"cannot both hold"**. That is **withdrawn**. It was the wrong comparison.

The Jahn et al. 1997 abstract says, verbatim (confirmed by independent search, not taken on
trust from an agent):

> "The slope between 0.001 and 0.01 mM GABA was 2.2 ± 0.4, indicating at least three
> binding sites for GABA."

0.001–0.01 mM is **1–10 µM**, against their EC₅₀ of 11.6 µM. So 2.2 is a **local slope on
the rising phase** (0.086–0.862 × EC₅₀), not a whole-curve Hill coefficient. The project
compared it against a regression over the **entire** curve.

For a Hill curve the two are identical — a Hill curve is a straight line in logit space —
which is exactly why this survived review for so long. For a receptor scheme they differ by
~0.4, because its logit curve **bends** between a low-agonist limiting slope equal to the
binding-site count (exactly 2 here) and a flatter slope through EC₅₀.

Redone with the measurement matched (33 × 31 grid):

| | whole-curve (as compared before) | rising phase (as Jahn measured) |
|---|---|---|
| steepest slope at `P_o,max ∈ [0.73, 0.77]` | 1.328 | **1.689** |
| lowest `P_o,max` reaching slope ≥ 1.8 | 0.9898 | **0.9084** |
| gap to the published 2.2 ± 0.4 | 2.18 σ | **1.28 σ** |

**1.28 σ is ordinary agreement.** There is nothing to explain away and no case, on this
evidence, for a different state diagram.

### Corroborations

1. Published **whole-curve** Hill fits for α1β2γ2 peak currents sit at **1.3–1.6**
   (nH = 1.5 ± 0.09, EC₅₀ 36 ± 6 µM, HEK293). That is where this scheme's whole-curve slope
   (1.294) already sits — the scheme agrees with the literature on the measurement the
   literature actually reports.
2. A **two-site** scheme reaches a rising-phase slope of **~2.0**, within 0.5 σ of 2.2. So
   the slope does *not* establish "at least three binding sites", and cryo-EM of the
   synaptic α1β2γ2 receptor shows **two**. The source's own inference does not follow.
3. **O. P. Hamill published a comment on exactly that inference** in the same NeuroReport
   issue — PMID 9480006. **Not read** (egress-blocked). This is the single most interesting
   unread document in the project and the first thing to fetch once §4 is unblocked.

### What replaced it as a genuine finding

**The PEAK observable is protocol-dependent.** Varying only the application length:

| application | 2 ms | 5 ms | 20 ms | 100 ms | 1000 ms |
|---|---|---|---|---|---|
| peak EC₅₀ | 299 µM | 113 µM | 45.6 µM | 39.2 µM | 39.2 µM |
| rising-phase slope | 1.77 | 1.79 | 1.76 | 1.69 | 1.69 |

A **sevenfold** EC₅₀ shift from the protocol alone. `PEAK_APPLICATION_MS` is fixed at 300 ms
and no dataset declares its own. This is the **third** observable-mismatch class, after
PEAK-vs-EQUILIBRIUM (P1) and absolute-vs-normalised (P4): same tag, same units, different
protocol. Jahn's 11.6 µM and this scheme's peak EC₅₀ are not strictly comparable either.

### The consequence for the dataset — read this before fitting anything

`JAHN1997_PEAK_CRC`'s nine points are a **global** Hill curve with nH = 2.2 spanning
0.3–3000 µM: **four decades of curve shape extrapolated from a slope measured over one.**
The source constrains the curve in three places only — the slope over 1–10 µM, the EC₅₀, and
saturation by 3 mM.

`sem` now widens outside the measured window. That does **not** move the fit, because the
*central values* carry the extrapolation, not their error bars — so the `D`-on-its-bound,
0.99-plateau fit is now **recorded as the reconstruction's artefact** rather than patched
over or quietly deleted.

**This is an open decision, not a solved problem.** See §5.1.

### The process lesson, which is the sharpest part

`fitting/data.py` recorded `"(slope over 1-10 uM)"` **in the same comment block as the
comparison**, and the comparison was made anyway. The qualifier was written down and not
acted on.

So the window is now an argument to `models.base.hill_slope` rather than a remark beside it.
**A note in prose does not constrain a computation; only a parameter does.** `hill_slope`
recovers 2.2 under both definitions on a true Hill curve — that is the control proving the
difference is a property of the scheme's curve, not of the measurement code.

---

## 2. Where the code is

New in PR #2:

| thing | where |
|---|---|
| `hill_slope(concs, response, *, window_rel=, band=)` | `src/circuitpharm/models/base.py` |
| `JAHN1997_SLOPE_WINDOW_REL` = `(1/11.6, 10/11.6)` | same, exported from `models/` |
| corrected dataset uncertainty | `src/circuitpharm/fitting/data.py` (`_jahn_spread`) |
| the whole corrected argument | `knowledge/12-inference.md` §2.0–2.6 |

Tests that pin the correction (all currently passing locally):

* `tests/test_p4_likelihood.py::test_the_slope_and_the_plateau_are_compatible_once_measured_the_same_way`
* `tests/test_p4_likelihood.py::test_two_binding_sites_can_produce_the_published_slope`
* `tests/test_p4_likelihood.py::test_the_peak_observable_is_protocol_dependent`
* `tests/test_p4_likelihood.py::test_mle_reproduces_what_the_source_actually_measured`
* `tests/test_provisional.py::test_the_published_slope_is_a_local_one_and_the_scheme_matches_it`

---

## 3. How to run it

```bash
python3.11 -m venv .venv && source .venv/bin/activate   # 3.13 does NOT work
pip install -e ".[dev]"

pytest -q -m "not slow"        # fast subset
pytest -q                      # everything; ~20 min
pytest -q tests/test_p4_likelihood.py -p no:randomly
```

Timing and gotchas worth knowing before you lose an hour to them:

* **The full suite takes ~20 minutes**, and the `full` CI job ~57 min against a 90 min
  timeout. That timeout was raised from 60 after a run was killed at exactly 60:00.
* **`-p no:randomly`** when comparing runs; test order is randomised by default.
* **R̂ is platform-dependent.** `test_posterior_recovers_known_parameters` was tuned because
  arm64/Accelerate and x86/OpenBLAS differ in the third decimal of split-R̂ through the
  matrix exponentials. Chains were lengthened (12000/4000 → 24000/8000) rather than the gate
  loosened. **If you are on an Apple Silicon Mac and see a marginal R̂ failure, that is this
  — do not "fix" it by widening the gate.**
* Fitting against any `synthetic=True` dataset raises a loud `ProvisionalResultWarning` on
  purpose. That is the machinery working, not a bug.

---

## 4. The data problem — THIS IS THE REAL BLOCKER

**No dataset in the repo is digitised, so no posterior is VALIDATED.** Everything built is
machinery demonstrated on generated or parametrically-reconstructed data.

**UPDATE 2026-10-11 — RESOLVED. A local session read all seven open-access targets.**
`pmc.ncbi.nlm.nih.gov` is reachable from outside the cloud container and serves full text,
so the block below was environment-specific, not a property of the literature. The findings
are in `knowledge/14-literature.md` and `scripts/literature_gamma2.py`. **Read those rather
than the table below**, which is kept because it records what was believed before, and
because several of its citations were WRONG — see the corrections after it.

The original (cloud-session) state, for the record:

```
pmc.ncbi.nlm.nih.gov   CONNECT tunnel failed, response 403
europepmc.org          CONNECT tunnel failed, response 403
rupress.org            CONNECT tunnel failed, response 403
discovery.ucl.ac.uk    CONNECT tunnel failed, response 403
WebFetch               getaddrinfo ENOTFOUND  (every host, incl. example.com)
```

Web **search** works, which is how the abstracts were read.

### The exact targets, in value order

| gap | target | what it contains | status |
|---|---|---|---|
| **mean open time** — would turn `FIT_FIXED_ALPHA` from a CONVENTION into a measurement and make the whole kinetic fit data-determined | Keramidas & Harrison 2008, *J Gen Physiol* 131(2):163–181, **PMC2213567**; and 2010, 135(1):59–75, **PMC2806416** | explicit open-time decompositions (τ + areas), intraburst P_o ~0.7 / ~0.9 | **α1β2γ2S — the SHORT splice variant.** Not this project's receptor. Usable only as a recorded substitution. |
| same | Li et al. 2008, *Br J Pharmacol*, **PMC2241790** | per-patch open-time components for α1β2γ2**L**, Table 1 has the averages | 50 µM GABA = sub-saturating (~EC₄₀) |
| same | Jahn et al. 1997, NeuroReport 8(16):3443–6 | **may not contain a mean open time at all** — the abstract gives BURST duration (10.3 ± 3.0 ms) | no PMC copy; paywalled via Ovid. Establish whether the number exists before chasing it. |
| **deactivation** — makes an absolute rate (`k_off`) identifiable; without it the fit determines only ratios | Dixon et al. 2014, *J Biol Chem* 289(8):5399–5411, **PMC3937617**; and 2015, *Br J Pharmacol*, **PMC4507157** | weighted deactivation τ for α1β2γ2L, plus pulse duration / concentration / temperature / components | a τ of 5.9 ± 0.5 ms (n = 10) is recorded **as a LEAD only** — a targeted search returned the citation **without** the number |
| **holdout** — P5's CV and P6's falsification bound both consume it | Barberis et al. 2007, ***Eur J Neurosci* 25(9):2726–2740**, **PMC1950087** | desensitisation onset | ⚠️ I cited 26(7) — wrong. ⚠️ Title says α1β2γ2 but the **Methods say γ2S**. The onset **cannot be reconstructed** (τ₂/τ₃ live only in Fig. 4C) and the paper is internally inconsistent on n (text 7, caption 3). **A better holdout is in the same paper**: paired-pulse recovery at 100 ms, 0.33 ± 0.03 (n = 8). |
| EC₅₀ reference table | Mortensen, Patel & Smart 2012, *Front Cell Neurosci* 6:1 | ⚠️ **I called this the authoritative "EC₅₀/nH table". It contains NO Hill coefficients.** The nH half does not exist. | open access |
| **the Hamill comment** | PMID 9480006, NeuroReport 8(16):**iv** | ⚠️ **One page, no abstract, no DOI, no PMC copy — its TITLE is all that exists.** That the title asks how many binding steps are involved is still independent evidence the three-site inference was contested in the same issue, but there is no argument to read. |

### Two ways to unblock

1. **Allowlist those hosts** in the cloud environment's network policy. Then one agent pass
   closes most of the table.
2. **Drop the PDFs in** from institutional access. Faster and more reliable.

**Grading convention, which is load-bearing here.** `MISSING_DATASETS` marks every claim
`[VERIFIED]` (read verbatim here) or `[LEAD]` (reported by a literature-search agent and
*not* independently confirmed). **No unconfirmed number was entered as data.** Keep that
discipline — a mis-citation in this project is worse than a gap, because a gap is visible.

### Mis-citation risks already spotted

* **Mortensen et al. 2010** (*J Physiol* 588:1251) is widely cited for α1β2γ2 but used
  **α1β3γ2** — β3, not β2. Check before citing.
* An Ovid page renders the Jahn title with **garbled subunit names** ("α1δ2γ2L"). PubMed has
  it right. Do not let a tool ingest the Ovid rendering.
* Oocyte EC₅₀s (41–51 µM) run ~4× higher than HEK293 (~5–16 µM). **Do not pool across
  expression systems.**

---

## 5. What I would do next, in order

### 5.1 Decide what `JAHN1997_PEAK_CRC` should be — an open design question

It currently keeps all nine points, with widened `sem` outside the measured window and the
artefact documented. That is honest but unsatisfying: a fit still draws shape information
from four decades the source never measured.

Three options, none obviously right:

| option | for | against |
|---|---|---|
| **keep as is** (current) | only sourced CRC in the repo; artefact documented and tested | the fit still lands on `D` at its bound and a 0.99 plateau |
| **restrict to the measured window** + a saturation anchor | asserts only what the source supports | leaves ~3 points; most downstream fits become uninformative — which may be the *honest* state |
| **split in two** — a defensible restricted dataset plus the full reconstruction kept to *demonstrate* the artefact | both uses served explicitly | more machinery; touches `TRAIN`, P5 comparison, P6 design |

I did not pick one unilaterally because it changes what every downstream fit means. **This
is the first thing to decide.**

### 5.2 Re-examine the P6 pre-registration in light of §1

`knowledge/13-preregistration.md` concludes the discriminating experiment does not exist at
σ = 0.02 and needs σ ≤ 0.005. That conclusion rests on posteriors fitted partly through the
PEAK observable — whose protocol dependence (§1) was not known when it was written. **The
σ ≤ 0.005 requirement should be re-derived with the application duration declared.** It may
soften, because some of the apparent indistinguishability may be protocol smearing.

### 5.3 ~~Decide whether `P_o,max = 0.750` means what the project thinks~~ — ANSWERED, and I had it wrong

**My original §5.3 was built on a conflation, and the recommendation in it was wrong.**

I wrote that "the one genuinely macroscopic estimate is the **0.69** by nonstationary
variance analysis". **That number does not exist as described.** 0.69 ± 0.02 (n = 11) is
Keramidas 2008's Table II **intraburst P_o for the M-mode**. The paper's actual
nonstationary-fluctuation-analysis value, quoted verbatim in
`scripts/literature_gamma2.py`, is **0.56** — and the authors mark it *unpublished data*,
with no error bar and no n.

So the finding is the opposite of what I implied: the only genuinely macroscopic reading is
**below** the 0.750 convention, not above it.

**Conclusion (local session, and I agree): `P_o,max = 0.750` stays a `CONVENTION`.** The 0.56
is γ2S rather than this project's γ2L, errorless, n-less, self-declared unpublished,
measured at a 1–2 ms application — the regime where this project's own PEAK EC₅₀ is ~299 µM
rather than 39 µM — and outside `FIT_RANGES`' [0.70, 0.80] in any case. Five of the six
candidate values are **intraburst** or **intracluster**, which are conditional on being
inside a burst and are not the population quantity a macroscopic scaling convention needs.

The `0.8 intraburst P_o` I also listed was **not found as a stated value in any of the seven
texts**; the nearest real numbers are Keramidas 2010's 0.81 ± 0.01 at 5 mM and Keramidas
2008's pair.

### 5.4 Then, and only then, re-fit

With a real mean open time, `FIT_FIXED_ALPHA` stops being a convention and the three-anchor
square system (3 anchors, 3 free rates, no residual DOF, **no confidence interval**) can
become an actual estimation problem.

---

## 6. Things that will bite you

* **`PEAK` is not one observable.** Always state the application duration. §1.
* **Normalised vs absolute** is a separate trap from PEAK-vs-EQUILIBRIUM; `Normalisation`
  exists because an optimiser closed a 0.25 plateau gap by deleting desensitisation.
* **The CV tie-break is load-bearing.** Model C *nests* Model B and scored 0.84 nats better
  out of sample on noise. `resolve_cv` needs the margin to clear both a paired per-fold SE
  and a 2-nat floor, then breaks ties on fewest parameters. Remove it and the protocol
  selects the larger model on noise.
* **The three-anchor fit has no confidence interval** — a square system, not a precise one.
* **Nothing here bears on the respiratory question.** Roadmap §6's last row is unchanged.

---

## 7. My own errors in this round, recorded so they are not repeated

* Quoted the rising-phase slope as **1.71** (Jahn's absolute 1–10 µM) where the helper
  computes **1.59** (same *relative* part of a curve whose EC₅₀ is not his). Both readings
  are now given rather than the flattering one.
* Pinned **1.69** from a 33 × 31 grid into a test that runs an **11 × 9** one, which does not
  contain that grid's best point. The test now asserts what the coarse grid produces (1.58).
* Set a `> 0.3` threshold that the measured **0.300** missed by nothing at all.
* Delegated the literature hunt and got back a report with real findings **and** numbers I
  could not confirm. The Jahn quote verified; the Dixon 5.9 ms did not. Only the verified
  one was recorded as data. **Treat agent output as a lead list, never as a source.**
* **And then broke that rule in the same document.** I labelled the P_o values LEAD-grade
  and *then built §5.3's recommendation on one of them* — "the one genuinely macroscopic
  estimate is the 0.69". It is not a macroscopic estimate at all; it is an intraburst
  M-mode P_o, and the real fluctuation-analysis number (0.56) points the other way. A
  grading discipline only works if it also stops you reasoning from the ungraded side.
  Corrected in §5.3.
* Four citations entered wrong from search metadata rather than the articles: Dixon 2014's
  issue (289(8) → **289(9)**), Barberis 2007's volume (26(7) → **25(9):2726–2740**), Dixon
  2015's splice variant (γ2L → **γ2S**), and calling Mortensen 2012 an "EC₅₀/nH table" when
  **it contains no Hill coefficients**. All four are in §4, flagged in place.

---

## 8. Pointers

| document | what it is |
|---|---|
| `knowledge/09-strategic-roadmap.md` | the plan; §8 has every phase CLOSED with a "what it actually found" block |
| `knowledge/12-inference.md` §2 | the full corrected slope argument — **start here** |
| `knowledge/11-identifiability.md` | what the data can determine, asked before anything was fitted |
| `knowledge/13-preregistration.md` | the P6 design, subject to §5.2 |
| `knowledge/08-manuscript.md` | the manuscript; numbers regenerated by `scripts/paper_numbers.py` |
| `fitting.data.MISSING_DATASETS` | the three gaps, with targets and `[VERIFIED]`/`[LEAD]` grading |

PRs: [#1 merged](https://github.com/willkhinz/circuitpharm/pull/1) ·
[#2 open](https://github.com/willkhinz/circuitpharm/pull/2)
