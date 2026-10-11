# 14 — The literature, read

Status: **the egress block is gone.** On 2026-10-10 `pmc.ncbi.nlm.nih.gov` was reachable
from a local session and **all seven open-access targets in `fitting.data.MISSING_DATASETS`
were read in full text.** `europepmc.org`, `rupress.org` and `discovery.ucl.ac.uk` still
return 403 — confirmed with a browser user agent, so a real block rather than UA sniffing —
but every target had a PMC copy, so it no longer matters.

Every number below is in `scripts/literature_gamma2.py`, which holds it with its conditions
and **prints** the derived quantities. Run it rather than trusting this prose:

```bash
python scripts/literature_gamma2.py
```

`tests/test_literature.py` pins the readings that other modules will consume.

---

## 0. The thirty-second version

1. **The Dixon 5.9 ms deactivation [LEAD] is CONFIRMED verbatim.** Promoted to [VERIFIED].
2. **No source states a mean open time for α1β2γ2L.** Four readings exist, spanning 5.1×,
   and the spread is **definition**, not disagreement.
3. **Deactivation is protocol-dependent, and this is now a published measurement rather
   than a model prediction.** Barberis 2007 measured the same receptor in the same study at
   **52.5 ms (2 ms pulse)** and **364 ms (3 s pulse)** — 6.9× from pulse duration alone.
4. **`P_o,max = 0.750` is being compared against the wrong quantity.** Five of the six
   candidate values are *intraburst* or *intracluster*. The only genuinely macroscopic
   reading is **0.56**, which is **below** 0.750, not above it.
5. **The handoff's "0.69 by nonstationary variance analysis" is a conflation of two
   different numbers in the same paper.** The 0.69 is an intraburst M-mode P_o; the actual
   nonstationary-fluctuation-analysis value is 0.56, and its own authors mark it
   *unpublished data*.
6. **The Hamill comment is one page long** (NeuroReport 8(16):**iv**), has no abstract, no
   DOI and no PMC copy. Its **title** is all that is available — and the title alone is
   informative.
7. **Three citation errors** in `MISSING_DATASETS`/`HANDOFF.md` are corrected, and
   **Mortensen 2012 contains no Hill coefficients at all.**

---

## 1. What was read

| source | read | receptor *as stated* | prep | temp |
|---|---|---|---|---|
| Keramidas & Harrison 2008, J Gen Physiol 131(2):163–181, PMC2213567 | full text | α1β2γ2**S** | HEK293 excised outside-out | 21 ± 1 °C |
| Keramidas & Harrison 2010, J Gen Physiol 135(1):59–75, PMC2806416 | full text | α1β2γ2**S** | HEK293 excised outside-out | 21 ± 1 °C |
| Li et al. 2008, Br J Pharmacol 153(3):598–608, PMC2241790 | full text | rat α1β2γ2**L** | HEK293 **cell-attached** | room temp (no value) |
| Dixon et al. 2014, J Biol Chem 289(**9**):5399–5411, PMC3937617 | full text | human α1β2γ2**L** | HEK293 outside-out / whole cell | room temp (no value) |
| Dixon et al. 2015, Br J Pharmacol 172(14):3522–3536, PMC4507157 | full text | human α1β2γ2**S** | HEK293 outside-out | 22 ± 1 °C |
| Barberis et al. 2007, Eur J Neurosci **25(9)**:2726–2740, PMC1950087 | full text | rat α1β2γ2**S** | HEK293 outside-out + lifted cells | 22–24 °C |
| Mortensen, Patel & Smart 2012, Front Cell Neurosci 6:1, PMC3262152 | full text | α1β2γ2**S** | HEK293 whole-cell | not stated |
| Jahn et al. 1997, NeuroReport 8(16):3443–3446, PMID 9427304 | **abstract only** | α1β2γ2**L** | HEK293, ultrafast exchange | not stated |
| Hamill 1997, NeuroReport 8(16):**iv**, PMID 9480006 | **unread** | — | — | — |

**Only two sources are γ2L — the project's receptor: Li 2008 (rat) and Dixon 2014 (human).**
Everything else is a splice-variant substitution and must be recorded as one.

### Citation errors corrected

* **Dixon 2014 is 289(9)**, not 289(8).
* **Barberis 2007 is Eur J Neurosci 25(9):2726–2740**, not 26(7).
* **Dixon 2015 is γ2S**, not a second γ2L deactivation source as `MISSING_DATASETS` claimed.
* **Mortensen 2012 reports no Hill coefficients.** It is described in `MISSING_DATASETS` and
  `HANDOFF.md` as "the most authoritative same-conditions EC₅₀/nH table". The nH half does
  not exist: Table 1 carries pEC₅₀ ± SEM (n), EC₅₀ and I_max only, and "Hill coefficient"
  appears exactly once, in the Methods, defining *n* in the fitted equation. A
  case-insensitive grep for `nH` returns 34 hits, **all of them author initials**
  ("Schindelin H.", "Shin H. S."). Checked because a string count is not a value.
* Also: Barberis's **title and abstract** say "α1β2γ2" while its **Methods** say γ2S. The
  same trap as Mortensen et al. 2010 (β3 behind an "α1β2γ2" citation), one level down.

---

## 2. Mean open time — the quantity does not exist as a single number

| reading | receptor | [GABA] | mean open time | conditional on |
|---|---|---|---|---|
| Barberis 2007, **stated** | α1β2γ2S | 10 mM | **1.42 ± 0.05 ms** (SEM, n=6) | a 4 s window **after** a 2 ms pulse — non-stationary |
| Keramidas 2008 M-mode, *derived* | α1β2γ2S | 10 mM | **2.68 ms** (n=10) | being inside a burst, medium gating mode |
| Li 2008 / Li 2007b, *derived* | α1β2γ2**L** | 50 µM | **2.96 ms** (n=4) | being inside a cluster; 50 µM is sub-saturating |
| Keramidas 2008 H-mode, *derived* | α1β2γ2S | 10 mM | **7.25 ms** (n=5) | being inside a burst, high gating mode |

**Spread 5.1×.** Barberis is the only source that *states* a mean open time; the other three
are **area-weighted means computed by this project** from the components the sources report:

* Keramidas 2008 Table III, M-mode: τ_O = 0.49 ± 0.02 / 2.58 ± 0.09 / 5.17 ± 0.38 ms,
  areas 0.26 / 0.49 / 0.25. **Error statistic not declared by the paper.**
* Keramidas 2008 Table III, H-mode: τ_O = 0.59 ± 0.11 / 4.22 ± 0.73 / 13.2 ± 2.6 ms,
  areas 0.18 / 0.41 / 0.41.
* Li 2008 Table 1: OT = 0.28 ± 0.05 / 3.0 ± 0.7 / 7.3 ± 3.2 ms (**s.d.**, n = 4),
  fractions 0.22 / 0.65 / 0.13.

No source publishes a covariance between τ_i and a_i, so **no error bar on the derived means
is computable** and none is quoted. The script prints a ±1-error sensitivity range instead,
labelled as not an error bar.

**The useful fact:** the two intraburst/intracluster readings nearest the project's regime
agree at **2.68 and 2.96 ms across both splice variants** — a γ2S/γ2L cross-check that
happens to be reassuring. Barberis's 1.42 ms is lower because it is non-stationary and
includes late brief singly-bound openings; H-mode's 7.25 ms is higher because it conditions
on the high mode. Neither is in conflict with the other two.

### Jahn 1997, settled as far as it can be

The abstract (re-read verbatim from PubMed) contains **no mean open time**. It reports:

> EC₅₀ 11.6 ± 0.9 µM · saturates with 3 mM GABA · slope over 0.001–0.01 mM = 2.2 ± 0.4
> · rise time ~120 ms at 0.001 mM → ~0.8 ms at 10 mM · **burst** duration 10.3 ± 3.0 ms
> · single channel slope conductance ~29 pS for >95% of the current

Burst duration is not mean open time. **Whether the full text states one remains unknown**
and cannot be settled from the abstract: NeuroReport 1997 has no PMC copy and is paywalled.
That is the honest end of this thread. (Jahn's ~29 pS agrees well with Barberis's
27.7 ± 0.8 pS, n = 5 — a small independent corroboration of the preparation.)

---

## 3. Deactivation — protocol-dependent, and now measured

| study | receptor | pulse | [GABA] | τ_w | n |
|---|---|---|---|---|---|
| Dixon 2014 | α1β2γ2**L** | ≤1 ms | 3 mM | **5.9 ± 0.5 ms** | 10 |
| Barberis 2007 | α1β2γ2S | 2 ms | 10 mM | **52.5 ± 2.9 ms** | 6 |
| Barberis 2007 | α1β2γ2S | **3 s** | 10 mM | **364 ± 40 ms** | 7 |
| Dixon 2015 | α1β2γ2S | ≤1 ms | 5 mM | *9.0 ms — **simulated**, not measured* | — |

Two things follow, and they are different.

**(a) Pulse duration alone moves τ_w 6.9×.** The two Barberis rows are the same receptor,
same patches, same study, same analysis — only the pulse differs. This is the deactivation
analogue of the PEAK application-length result in `12-inference.md` §2.3, except that it is
a *published measurement* rather than a property of this project's scheme.

**(b) The Dixon/Barberis 8.9× gap is NOT pulse duration** — ≤1 ms versus 2 ms is nearly the
same protocol. It is the **fit**:

* Dixon fitted **two** exponentials; Barberis fitted **three**.
* Barberis's slow component (221.35 ms, area 0.20) contributes **44.3 of its 52.5 ms**.
* Dixon's weighted 5.9 ms sits beside Barberis's **fast component, τ₁ = 2.8 ± 0.3 ms**.

So the two are plausibly compatible *on the fast component* and differ on whether the slow
tail was captured and fitted at all. A weighted deactivation τ is therefore only defined
once **pulse duration, component count and fit window** are all declared. The project's
dataset tags currently carry none of the three.

### The one reconstructable observable in the whole corpus

Barberis's 2 ms-pulse deactivation is **fully specified by the source**: three time
constants, three areas, all with errors, and the functional form the authors themselves
fitted. Recomputing Σ A_i τ_i gives **53.55 ms** against a published τ_w of 52.5 ± 2.9 —
**0.36 SEM**, which both validates this project's transcription and confirms the parameter
set is complete and self-consistent.

**This matters for the `JAHN1997_PEAK_CRC` question (§5.1 of the handoff, Task 2).** Jahn's
CRC is a reconstruction that extrapolates four decades of curve shape from a slope measured
over one. Barberis's deactivation is the opposite case: a curve reconstructed from
parameters the source states, over the window the source measured, with **no extrapolation**.
It is the first observable in this project that could carry `synthetic=False` honestly —
subject to the splice-variant substitution (γ2S) and the protocol declaration.

---

## 4. The holdout — partly unusable, with a better alternative in the same paper

Barberis 2007, desensitisation onset, 3 s pulse of 10 mM GABA, rat α1β2γ2S:

* τ₁ = **2.9 ± 0.1 ms**, A₁ = **0.56 ± 0.025**
* steady-state component weight **0.076 ± 0.013**
* steady-state:peak **measured at 200 ms** = **0.21 ± 0.02** (Table 1)
* "in ~3 ms the current is reduced by more than one half"

**τ₂ and τ₃ exist but are reported only in Figure 4C**, so the onset time course **cannot be
reconstructed** — the opposite of the deactivation case above. And the source is
**internally inconsistent about n**: the text says **n = 7**, the Figure 4 caption says
**three patches** for α1β2γ2. Recorded as an ambiguity, not resolved. (The Fig. 3 captions
*do* agree with the text at n = 6/9, so this is specific to Fig. 4.)

**A better holdout sits in Table 1 of the same paper:** paired-pulse recovery at a **100 ms
gap = 0.33 ± 0.03** (Fig. 5 caption: n = 8), with 2 ms pulses of saturating GABA. A stated
interval, a stated value, no reconstruction needed. This is the cheapest honest holdout
available and it should be preferred over the desensitisation onset.

Other fully-stated Barberis scalars, usable the same way: 10–90% rise time at 10 mM GABA
= **0.29 ± 0.02 ms**; single-channel conductance **27.7 ± 0.8 pS** (n = 5).

---

## 5. What `P_o,max = 0.750` is — Task 4, answered

| value | receptor | n | quantity | conditional on |
|---|---|---|---|---|
| 0.87 ± 0.02 | α1β2γ2S | 5 | intraburst P_o, **H-mode** | inside a burst **and** the high mode |
| 0.81 ± 0.01 | α1β2γ2S | 7 | intraburst P_o, modes pooled | inside a burst |
| 0.69 ± 0.02 | α1β2γ2S | 11 | intraburst P_o, **M-mode** | inside a burst **and** the medium mode |
| **0.56** | α1β2γ2S | not stated | **P_o macropatch, nonstationary fluctuation analysis** | **nothing — this is the population quantity** |
| 0.56 ± 0.04 | α1β2γ2**L** | 3–8 | intraburst P_o | inside a burst |
| 0.42 ± 0.04 | α1β2γ2**L** | 4 | cluster open probability | inside a cluster; sub-saturating |

**Five of the six are conditional.** A macroscopic model's peak-scaling convention needs the
**population peak**, and the only reading of that quantity is **0.56** — which is **0.19
below** the 0.750 convention, not above it. Every value at or above 0.750 is conditioned on
being inside a burst or cluster.

### Two corrections to the handoff's §5.3

1. **"0.69 by nonstationary variance analysis" does not exist.** It conflates two different
   numbers in the *same paper*. The 0.69 ± 0.02 (n = 11) is Keramidas 2008 **Table II's
   intraburst M-mode P_o**. The actual nonstationary-fluctuation-analysis value in that
   paper is **0.56**, stated verbatim as:

   > "we estimated the channel open probability (P_O macropatch) in response to 1–2-ms
   > exposure to agonist using nonstationary fluctuation analysis. The P_O macropatch
   > values we obtained were … for α1β2γ2S GABA_A Rs, GABA–0.56 and THIP–0.42
   > **(unpublished data)**."

   The handoff concluded "the one genuinely macroscopic estimate is the 0.69" and built its
   Task-4 recommendation on that. The macroscopic estimate is 0.56, and **its own authors
   mark it unpublished** — no error bar, no n.

2. **Beware a second 0.69.** Keramidas 2010 also reports 0.69 ± 0.02 (n = 17) — for
   **α3β3γ2S at 200 µM GABA**, a different receptor entirely. Two unrelated 0.69s across the
   two papers.

3. The **"0.8 as an intraburst P_o"** lead was **not located** in any of the seven texts. The
   nearest readings are Keramidas 2010's 0.81 ± 0.01 at 5 mM and Keramidas 2008's H-mode
   0.87. Recorded as UNLOCATED rather than silently matched to either.

### What this does *not* settle

It settles **which quantity** 0.750 should be (the population peak) and **which readings are
not that quantity** (five of six). It does **not** establish that 0.750 should become 0.56:

* 0.56 is γ2**S**, unpublished-by-its-authors, errorless and n-less.
* It was measured at a **1–2 ms** application — by this project's own protocol-dependence
  table, the regime where peak EC₅₀ is ~299 µM, not the 39 µM of the 300 ms
  `PEAK_APPLICATION_MS`. So it is not even comparable to the project's PEAK observable.
* `FIT_RANGES` admits `po_max` only on [0.70, 0.80], so 0.56 would be **rejected** by the
  module's own validation — the same wall the `initial_paper.md` 0.84 hit from the other
  side.

So the finding is **a mismatch to record, not a value to swap in.** The `_NOMINAL_DEFECT`
basis should stay `CONVENTION`, with the note changed from "a convention pending
measurement" to "a convention that the available measurements **cannot** replace, because
the measured quantity nearest to it is conditional and the one unconditional reading is 0.56
at an incomparable protocol." That is a stronger and more honest statement than the roadmap's
P0-13 currently makes.

---

## 6. The Hamill comment — what is actually available

PubMed (via NCBI eutils, PMID 9480006):

> Neuroreport. 1997 Nov 10;8(16):**iv**.
> **How many transmitter binding steps are involved in opening fast receptor-gated
> channels?** Hamill OP.
> Comment on Neuroreport. 1997 Nov 10;8(16):3443-6.

Page **iv** — a single page in the issue's front matter. **No abstract, no DOI, no PMC
copy.** The comment relationship to Jahn 1997 is confirmed; the **content is unread** and
cannot be obtained from any reachable host.

The **title is itself the finding.** It asks how many binding steps are involved — which is
precisely the inference Jahn drew from the Hill slope ("the slope … was 2.2 ± 0.4, indicating
at least three binding sites for GABA"). A contemporaneous one-page commentary in the same
issue, titled as a question about exactly that inference, is **independent corroboration
that the three-binding-site claim was contested at publication**. PR #2 withdrew this
project's "structural conflict" on the grounds that a two-site scheme reaches a rising-phase
slope of ~2.0; Hamill evidently questioned the inference in 1997 on grounds unknown.

**This remains the highest interest-per-page unread document in the project.** It needs
institutional access to NeuroReport (Ovid) — one page. Nothing else blocks it.

---

## 7. EC₅₀, and a new observable mismatch

| source | receptor | EC₅₀ | prep | application |
|---|---|---|---|---|
| Jahn 1997 | α1β2γ2L | 11.6 ± 0.9 µM | HEK293 | **ultrafast exchange** (peak) |
| Mortensen 2012 | α1β2γ2S | **6.6 µM** (pEC₅₀ 5.180 ± 0.0593, n = 34) | HEK293 whole-cell | **20–30 ms delivery latency** |

Both HEK293, and they differ 1.76×. These are **not the same observable**: Mortensen states
GABA is "delivered … with a latency of 20–30 ms", which is slow enough that the measured
current is not a fast-application peak. The project's own PEAK table makes the direction
predictable — longer effective application lowers peak EC₅₀ — and 6.6 < 11.6 is consistent
with that. **Do not pool them,** and note this is a *third* axis beyond the oocyte/HEK293
split already recorded (41–51 µM vs ~5–16 µM).

Mortensen's I_max for α1β2γ2S is 2230 ± 193 pA (n = 18).

---

## 8. What is still missing after this pass

1. **A mean open time for α1β2γ2L stated by a source.** Does not exist in the corpus read.
   Li 2008's components are the closest, are derived, and are at a sub-saturating 50 µM in
   the **cell-attached** configuration.
2. **The Jahn full text** — to settle whether a mean open time is in it. Paywalled, no PMC.
3. **The Hamill comment's content** — one page, Ovid.
4. **Barberis Figures 4C and 6C** — τ₂/τ₃ of the desensitisation onset, and the steady-state
   intraburst open-time distribution. Both are plots; the numbers are not in the text. These
   would need figure digitisation, which is a different and lower-grade activity than
   reading a table, and must be graded accordingly if ever done.
5. **An error statistic for Keramidas 2008 and Dixon 2014.** Neither paper declares SEM or
   SD. Their error bars are recorded as UNSTATED rather than assumed to be SEM.
6. **A resolution of Barberis's n = 7 / n = 3 inconsistency** for the desensitisation onset.

None of these blocks the next step. What blocks the next step is a **decision**, not a
document: which protocol each observable declares. See `MISSING_DATASETS` and the roadmap.
