# Review of `initial_paper.md` against the model it claims to describe

Reviewed draft: `~/.gemini/antigravity/brain/e36f005a-bd4f-4eac-b808-352a4cf19933/initial_paper.md`
(699 lines, 80,686 bytes, modified 2026-10-08 19:00).
Reviewed against: working tree at `2ff9c76` plus the uncommitted `substrate.py` migration.
Every number below was re-derived by running the repository, not read off the draft.

Reproduce the whole verification with:

    python scripts/paper_numbers.py
    python scripts/ranking_robustness.py --draws 20000 [--a5-floor 0.01|0.02]

## RETRACTION, 2026-10-08 — §2 and §1b of this review are withdrawn

**The central finding below is wrong, and the draft was right.** This notice is first because
anyone reading the review needs it before anything else in it.

§2 claimed the draft's kinetic numbers were "not reproducible from this repository" and that
its parameterisation existed in no commit. In fact the draft was quoting the repository's **own
pharmacology calibration**, and the audit that "refuted" it made three errors:

1. **It grepped for the wrong symbol.** `git log -S"SYNAPTIC_PEAK_UM = 3000"` finds nothing
   because the value lives in `config.SYNAPTIC_PULSE`, which has held
   `(peak_um=3000.0, clear_ms=1.00)` all along — exactly the draft's transient. A third copy
   sat as an inline literal in `cpg.Drug.from_kinetics`. The absence of a string was read as
   the absence of a value.
2. **It compared against the wrong one of two calibrations.** `gabaa_kinetics` carried module
   defaults of 1000 µM / 0.30 ms, reached only by callers that pass no pulse. Every
   pharmacology consumer — `cpg.py`, `evaluation.py`, and four calibration scripts — passes
   `config.SYNAPTIC_PULSE` explicitly. So `scripts/paper_numbers.py`, which I wrote, was the
   only thing in the project using the module defaults, and the manuscript I "corrected" was
   the only document describing a calibration nothing else used.
3. **The `FIT_RANGES` argument was wrong too.** I claimed `Po_max = 0.84` "would have been
   rejected" by the declared range [0.70, 0.80]. Both calibrations hit the fit target
   `po_max = 0.7500` exactly. The draft's 0.8382 is **β/(α+β)**, correctly computed for the
   config-pulse parameters — a different quantity, not an out-of-range anchor.

On the config pulse the draft reproduces to every digit: `kon` 0.0146842, `koff` 0.469177,
`beta` 0.559023, `alpha` 0.107891, `Kd` 31.95, `E` 5.1814, `c_affinity` 2.7854, tonic gain
7.2142, phasic gain 1.0623, tonic headroom 210.74, τ ratio 1.6564, charge ratio 1.6552, and
every cell of its sensitivity table (351.1, 117.4, 229.8, 180.9, 54.7, 829.7, 3295.5, cap
crossing 5.40 µM). Its `c = 2.79 → s_max = 2.50` was right, and so was its `> 50×` across
0.2–0.8 µM, which I had "corrected" to `> 48×`.

**How this was found.** Not by re-reading the review. By writing
`scripts/decompose_burst_change.py`, which built a `Drug` through the circuit path and printed
its tonic gain as **7.214** — the draft's number — beside a manuscript claiming 7.876. One
quantity, two values, visible only because something finally computed it twice.

**The fix is structural.** `gabaa_kinetics.SYNAPTIC_PEAK_UM`/`SYNAPTIC_CLEAR_MS` now import
from `config.SYNAPTIC_PULSE`, and `cpg.py`'s literal is gone, so the constant has **one**
definition. `tests/test_review_regressions.py` pins that. The manuscript is regenerated on the
unified calibration, which restores the draft's figures.

**What this review got right, and which still stands:**

* The non-existent commit hash `c8f17a9` (§2).
* The **`P_o,max` conflation** (§Table 1) — the draft labelled β/(α+β) as the calibration
  anchor. The anchor is 0.750; β/(α+β) = 0.8382 is unattainable because desensitisation
  competes during the rise. Real finding, wrong supporting argument.
* **The falsification intervals** (§7b) — two of four are violated by the model *on the draft's
  own calibration*: 0.1 µM predicts 6.98–7.66× against a pre-registered > 15×, and 3.0 µM gives
  3.00–3.10× against < 3.0×. This was the review's most consequential finding and it is
  independent of the pulse confusion.
* **The citation errors** (§3): five wrong coordinates, three unresolvable.
* **The claim-support failures** (§7 of the manuscript): six load-bearing citations that do not
  support their number, including Kasugai on `f_extra,α5` and Walters on the diazepam benchmark.
* **The disowned statistics** (§5): Table 6 carried the median and 95th percentile that
  `ranking_robustness.py` prints NOT QUOTABLE.
* **The "multiscale" mislabel** (§6).
* **The missing substrate-independence result** (§7).

**The lesson, which is the project's own and which I committed anyway.** A value duplicated
across locations diverges, and the symptom is plausible rather than loud. I had written that
sentence into `substrate.py`, into a memory file, and into this repository's worklog, and then
diagnosed a three-way duplication as fabrication. The check that would have caught it is the
one I eventually ran by accident: compute the quantity twice, by different paths, and compare.

---

## Verdict in one line

The **algebraic half of the paper is sound and reproduces exactly**; the **kinetic half quotes
numbers that no commit of this repository produces**, and **8 of its 20 references are wrong or
unresolvable** — including four that carry load-bearing parameters.

---

## 1. What verified exactly

These were recomputed from the repo and matched the draft to the digit:

| Draft claim | Source | Status |
|---|---|---|
| Nominal R: ideal 10.88, alogabat 8.64, MP-III-022 7.88, HZ-166 1.21, BZ 1.00, neurosteroid 0.56, gaboxadol 0.00 | `subtypes.py` at ρ=6.0 | **exact** (10.8750, 8.6442, 7.8750, 1.2083, 1.0000, 0.5633, 0.0000) |
| Stereoisomer pair R = 9.35 / 3.63 | `PROFILES['sh053_R'/'sh053_S']` | **exact** (9.3487, 3.6250) |
| Table 6A input vectors (`f_prebotc`, `f_forebrain`, `f_extra`, `w_subj`, all seven efficacy vectors) | `subtypes.py` | **exact** |
| MC 5th percentiles 2.51 / 2.49 / 2.57; medians 16.05 / 12.59 / 40.88; 95ths 84.85 / 60.98 / 218.65 | `ranking_robustness.py`, seed 20261007 | **exact** |
| Floored variants: alogabat median 16.05→10.79→7.57, 95th 84.85→23.42→14.31 | `--a5-floor 0.01 / 0.02` | **exact** |
| "38% of draws put preBötC α5 below a tenth of nominal" | same | **exact** |
| Neurosteroid P(R>1) = 0.0% | same | **exact** |
| Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1, seed `20261007` | the venv itself | **exact** |
| Dirichlet CV formula √((1−p)/(p(α₀+1))) and Table 6B | derived | **correct** |
| Analytic cap-crossing identity, ratio = 1 + 2/(Ax) + 1/(Ax²) | derived | **correct form** (see §2 for the value) |
| Stereochemical invariance argument (2D descriptors cannot separate enantiomers) | `chirality.py` | **correct, and the right conclusion** |

§4.2–4.5 and §3 are, with the exceptions noted below, publishable as they stand. The Dirichlet
dispersion treatment in §4.4.2 is better than the repo's own documentation of it: the draft
derived the per-component CV and found the 175% CV on preBötC α5, which explains the 38% tail
the script reports but never explained.

---

## 1b. The two halves were produced differently, and that is diagnostic

§4.4.3 describes the Monte Carlo floor as "a genuinely floor-preserving transformation" and
writes out

    q5 = max(f5, L),  f5' = q5,  f_j' = f_j (1 - q5)/(1 - f5)

which is, line for line, the **uncommitted** fix in the working tree's
`scripts/ranking_robustness.py` — and it quotes 7.57 / 14.31, the numbers that fix produces,
not the 7.63 / 14.5 that `HEAD` produces. The draft therefore had this working tree and
**genuinely ran** the Monte Carlo against it.

Which means the kinetic numbers in §2 and Table 4 were **not** run against anything. They were
not produced by an older commit — `git log -S` rules that out — and they were not produced by
the tree the author demonstrably had. An author with the code in hand, who ran one section and
reported figures for the other, produced a section that is internally consistent and
externally unmoored. Internal consistency is why reading it does not catch the problem: Kd,
E, P_o,max, the asymptote and the ratio are all mutually coherent at `Po_max = 0.84`.

**Practical consequence for the revision:** trust §3 and §4.2–4.5 as run. Re-derive all of §2,
Table 4 and §5.2 from `scripts/paper_numbers.py`. Do not spot-check §2 — replace it.

---

## 2. The kinetic numbers are not reproducible from this repository

The draft states its numbers were "directly generated from the verified rate implementation."
They were not generated from *this* implementation. Two constants differ at the root:

| Quantity | Draft | `gabaa_kinetics.py` |
|---|---|---|
| `P_o,max` calibration anchor | 0.84 | **0.75** (`FIT_TARGETS`) |
| Synaptic peak `[G]` | 3.0 mM | **1.0 mM** (`SYNAPTIC_PEAK_UM`) |
| Synaptic `τ_clear` | 1.0 ms | **0.30 ms** (`SYNAPTIC_CLEAR_MS`) |

Three facts make this a defect in the draft rather than a stale repo:

1. `git log -S"po_max=0.84"` and `git log -S"SYNAPTIC_PEAK_UM = 3000"` return **nothing**.
   Those values have never existed in this repository's history.
2. `FIT_RANGES` declares `po_max` acceptable only on **[0.70, 0.80]**. The draft's 0.84 would
   have been **rejected by the module's own validation**.
3. The draft cites commit **`c8f17a9`** twice as its reproducibility anchor. `git cat-file -t
   c8f17a9` → *"Not a valid object name"*. It exists in no ref.

The downstream numbers all move, and the ones that move most are the ones the argument rests on:

| Draft | Repo | |
|---|---|---|
| Kd = 31.95 µM | **29.69 µM** | |
| E = β/α = 5.181 | **4.819** | |
| P_o,max = 0.8382 | **0.7500** | |
| P_open,∞ = 0.166169 | **0.156377** | closed form and mechanism-pushed agree to 7.5×10⁻⁶ |
| P_open,base(0.40 µM) = 7.8847×10⁻⁴ | **8.4704×10⁻⁴** | |
| **tonic headroom 210.7×** | **184.6×** | |
| c_affinity for s_max = 2.5: 2.785 | **2.9321** | |
| **phasic peak gain 1.062×** | **1.319×** | |
| **phasic max headroom 1.124×** | **1.668×** | |
| charge ratio 1.655× | **2.419×** | |
| τ ratio 1.656× | **1.895×** | |
| tonic gain 7.214× | **7.876×** | |
| synaptic peak P_o 0.6683 | **0.4196** | |
| synaptic charge 13.05 ms | **7.454 ms** | |
| 2.5× cap binds above 5.4 µM | **5.06 µM** | |

**Why this is not a rounding complaint.** The paper's thesis is the contrast between a
near-saturated synapse and a wide-open extrasynaptic space. The abstract, §4.1, §5.1 and the
Conclusions all lead with "1.062× / 1.124×" as the synaptic figure. At the repo's
parameterisation the synaptic peak gain is **1.319×** and its ceiling **1.668×** — in
deviation-from-unity terms, **5× and 5.4× larger**. The headline ratio of the two compartments
falls from the draft's 187× to **111×**.

The qualitative conclusion survives: two orders of magnitude still separate the compartments,
and the fixed 2.5× scalar cap is still wrong for the extrasynaptic pool. But a reader
reproducing the draft from this repository finds nothing matching, and the specific claim that
affinity modulation is "nearly inert" at the synapse is substantially weaker at 1.32× than at
1.06×. **Every figure in §2 and Table 4 must be regenerated.** `scripts/paper_numbers.py` now
emits them so they cannot drift again.

Also: the draft's Table 3 lists `c = 2.79` as `s_max = 2.50`. In this parameterisation
c = 2.79 gives **s_max = 2.401**; s_max = 2.50 needs c = 2.9321. The draft conflated the
microscopic k_off multiplier with the operational EC50 shift in the one table whose purpose is
to keep them apart.

### The sensitivity sweep's conclusion holds, with different numbers

| Sweep | Draft range | Repo range |
|---|---|---|
| d/r halved / doubled | 351.1× / 117.4× | **309.9× / 102.4×** |
| β ±2-fold | 180.9–229.8× | **159.8–200.2×** |
| ambient 0.2–0.8 µM | 829.7× – 54.7× | **725.8× – 48.1×** |

The claim "a large extrasynaptic dynamic range (>50×) robustly persists across 0.2–0.8 µM"
**survives** — the repo's worst case in that span is 48.1×, which is marginally *below* 50.
State it as **>48×**, or as ">50× across 0.2–0.7 µM", rather than letting a round number sit
one unit on the wrong side of the computed value.

---

## 3. Citations: 8 of 20 defective

Resolved against Crossref by title, then by journal/volume/page. Four of the eight carry
load-bearing parameters.

### Wrong coordinates

| Draft | Actual | Carries |
|---|---|---|
| Crestani et al. (**2001**) PNAS 99(13) **8993–8997** | **2002**, PNAS 99 **8980–8985** | α5 trace-fear rationale |
| Walters et al. (2000) ***Br. J. Pharmacol.* 131(7) 1307–1314**, "Evaluation of the allosteric interaction…" | **Nat. Neurosci. 3:1274–1281**, "Benzodiazepines act on GABA-A receptors via two distinct and separable mechanisms" | **diazepam s_max = 2.50 — the reference arm of the entire index** |
| Haas & Macdonald (1999) ***J. Neurosci.* 19(7) 2435–2445**, "…composition determines deactivation kinetics" | ***J. Physiol.* 514:27–45**, "…γ2 and δ subtypes confer unique kinetic properties…" | **the kinetic topology itself** |
| Kasugai et al. (2010) ***J. Neurosci.* 30(42) 14024–14035** | ***Eur. J. Neurosci.* 32:1868–1888** | **`EXTRASYN['a5'] = 0.80` — the whole headroom argument** |
| Otis & Mody (1992) ***J. Physiol.* 454(1) 477–496** | no such paper; the decay-kinetics result is **Otis & Mody 1992, *Neuroscience* 49:13–32** | the τ_IPSC-prolongation phenotype |

The Kasugai error is the one already on record: `provenance.report()` says *"a5_dist resolves
perfectly by DOI and does not support the numbers it is the obvious candidate for."* The draft
reached the same source independently and cited it to the wrong journal.

### Unresolvable

| Draft | Finding |
|---|---|
| Nutt et al. (2007) "Development of α5-selective GABA-A modulators", *Neuropharmacology* 53(7) 810–820 | **no match** by title or by journal/volume/page |
| Fischer et al. (2010) "Anxiolytic-like effects of MP-III-022", *Neuropharmacology* 59(7–8) 612–618 | **no match**, and **anachronistic** — the earliest primary MP-III-022 literature Crossref indexes is 2024–2026 (Lyu et al.; Mirković et al. on the related GL-II-73). Carries **MP-III-022's s_max = 2.20** in Table 2 and §4.3. |
| Saba et al. (2017) *Alcohol. Clin. Exp. Res.* 41(4) 748–758 | **no match**. Carries the α5-alcohol rationale **and `w_subj,δ = 0.0`**, the assignment that produces gaboxadol's R = 0.00 |

Verified correct: Cecere 2025, Cheng 2006, Farrant & Nusser 2005, Jones & Westbrook 1995,
Löw 2000, McKernan 2000 (Nat. Neurosci. 3:587–592 ✓), Nutt 2006 (J. Psychopharmacol. 20:318–320 ✓),
Olsen & Sieghart 2008, Pirker 2000, Rudolph 1999, Rudolph & Möhler 2014, Rundfeldt & Löscher
(online 2013, 2014 issue — acceptable as cited).

**Three unresolvable citations supporting three model parameters is the same failure the project
already documented for its own source registry.** The fix is the same: resolve or delete. A
parameter whose only support is an unresolvable citation should be labelled UNSOURCED, not
footnoted.

---

## 4. The draft overstates provenance that the repo itself rates lower

`provenance.report()` says **6 of 28 load-bearing parameters name a source (21%)**. The draft
attributes Table 5 to "regional expression literature (Pirker et al., 2000; Kasugai et al.,
2010)". The repo's audit disagrees about which cells those sources reach:

| Table 5 / 6A cell | Draft attribution | `provenance.py` |
|---|---|---|
| `prebotc.a5 = 0.02` | "model assumption reflecting this anatomical ordering" | **UNSOURCED** |
| `forebrain.a1, a23, d_a4, eps` | Pirker et al. 2000 | **all four UNSOURCED** |
| `EXTRASYN.a1 = 0.15`, `a23 = 0.20` | implied Kasugai | **both UNSOURCED** |
| `EXTRASYN.eps = 0.50` | "50% (Assumed)" ✓ | **GUESS** ✓ agrees |
| `w_subj` for a1, a23, d_a4, eps | Saba et al. 2017 | **all four UNSOURCED** |
| `KAPPA = 15.0` | presented with full derivation, no justification of the value | **UNSOURCED** |

The draft is honest about `prebotc.a5` and `eps` in prose. The problem is the table headers and
the §4.2 provenance note, which name two real papers over a column that is mostly unsourced.
**Table 5 should carry a per-cell provenance column generated from `provenance.py`**, so the
manuscript cannot claim more support than the code records.

---

## 5. The draft quotes the exact statistics the model marks NOT QUOTABLE

`ranking_robustness.py` prints, verbatim:

> NOT QUOTABLE: the median and 95th percentile. 38% of draws put preBötC a5 below a tenth of
> nominal … so the upper tail measures the prior.

The draft's **Table 6 makes "MC Median R" and "95th Percentile R" two of its six columns**, and
§4.4.3 leads with the medians. §4.4.3 item 2 does then explain the prior-sensitivity correctly
and quantitatively — so the draft knows. But a reader who reads Table 6 and stops has taken two
numbers the model explicitly refuses to stand behind.

**Fix:** drop the median and 95th columns from Table 6, keep the 5th percentile and the
ordering, and move the medians into §4.4.3 as a stated prior-sensitivity diagnostic.

Also: the script discards 10 of 20,000 draws for the ideal arm as non-finite (burden → 0, the
draws most favourable to it). The draft says "20,000-draw" without the discard. Minor, but it is
the kind of omission that makes a reproduction look like a mismatch.

Smaller: HZ-166's P(R>1) is **63.8%** at seed 20261007; the draft says 64.8%.

---

## 6. §4.4 is called multiscale and is not

The draft titles §4.4 "**Multiscale Circuit** Selectivity Pipeline" and attributes ρ = 6.00 to
`circuitpharm.evaluation`. `ranking_robustness.py` imports:

    from circuitpharm.subtypes import (REGIONS, SUBTYPES, SUBJECTIVE_WEIGHT, EXTRASYN, PROFILES)

and nothing else. No neuron, no circuit, no simulation — the script's own docstring says *"This
is analytic -- no circuit simulation."* R is a ratio of weighted subunit sums.

This matters because the project *does* have a circuit layer, and the paper's credibility
depends on not blurring which claims came from it. **Retitle to "Algebraic Selectivity Index"**
(the draft's own §4.4.1 heading already says "Algebraic Pipeline" — the section heading simply
contradicts its first subheading).

---

## 7. What the draft is missing: the project's newest and strongest result

The draft contains **no circuit-simulation result at all**, and the project's best evidence that
the ranking is not an artefact of one neuron model is absent.

**The link-5 result** (`scripts/compare_substrates.py`, 3 seeds, occupancy 1.0, each arm as a
fractional reduction in mean inspiratory output from *its own* substrate's control):

| arm | LIF | conductance | Δ |
|---|---|---|---|
| neurosteroid | +0.610 | +0.994 *(alive 0/3)* | +0.384 |
| mp_iii_022 | +0.020 | +0.071 | +0.051 |
| ideal_a5 | +0.005 | +0.051 | +0.046 |
| alogabat | +0.002 | +0.043 | +0.041 |
| hz_166 | +0.001 | −0.100 | −0.101 |
| **nonselective_bz** | **+0.079** | **−0.377** | **−0.456** |

* **Spearman ρ = +1.0000 over the five subtype-selective arms.** Identical ordering —
  `neurosteroid > mp_iii_022 > ideal_a5 > alogabat > hz_166` — on an integrate-and-fire cell and
  on a Butera–Rinzel–Smith conductance cell with a real I_NaP, voltage-gated inactivation, an
  11.8× Mg²⁺-relief span and reachable depolarisation block. The two substrates share almost no
  mechanism. **This is a far better robustness argument than the Dirichlet sweep**, because it
  perturbs the model's *structure* rather than its parameters, and it belongs in §4.
* **Exactly one arm moves, and only by changing sign.** The non-selective BZ goes from
  2nd-most-burdensome (+0.079) to least (−0.377 — a 38% *increase* in output). It is the arm
  with the largest α1 efficacy acting on a preBötC the model treats as α1-predominant
  (`REGIONS` α1 = 0.60), so the flip is largest exactly where the drug effect is largest.
  Confirmed by two independent methods (FFT and direct burst counting): mean output +41.7%
  cond vs −8.6% LIF, duty cycle **+65.5%** cond vs −6.7% LIF, peak essentially unchanged on both.
* **Not established:** whether that is more bursts or longer bursts. `burst_metrics` reported
  5.63 Hz against the FFT's 0.302 Hz (19×) because a 0.35×max threshold fires repeatedly within
  one burst. A burst-level detector for the conductance substrate is outstanding.

**This cuts directly against the draft's own reference arm.** Diazepam is the denominator of
every R in Table 6, and on the more biophysically detailed substrate its sign on respiratory
output reverses. The draft should say so.

### Other declared limitations the draft omits

* **A4 (riluzole dissociation) is NOT reproduced** — 0/3 seeds. This is *inherited*: BRS model 1
  **is** the pacemaker hypothesis, and riluzole is the principal published argument against it.
  Declared, not retuned.
* **The coupling weights are ours, not published.** BRS part II was unreachable (journal 403;
  every available encoding is single-cell).
* **Respiratory calibration remains wet-lab blocked**, and link 4 moved the >P12 muscimol anchor
  *further* away by importing neonatal parameters.
* **The conductance substrate is anchored to an in vitro band** (`INVITRO_BAND` 0.05–1.00 Hz,
  neonatal rat, Revill et al. 2021), not the in vivo `EUPNOEA_BAND`. Any respiratory number from
  it is an in vitro slice claim.

---

## 7b. The pre-registered falsification intervals do not match the model

This is the most consequential single finding, because §5.2 is the manuscript's falsification
test — the part a wet lab would actually run. The draft pre-registers four intervals for a
high-efficacy affinity PAM (s_max ≈ 2.4–2.5). Regenerated from the repo:

| `[GABA]_bath` | draft pre-registers | repo at s_max 2.40–2.50 | asymptote | verdict |
|---|---|---|---|---|
| 0.1 µM | **R_max > 15×** | **7.66 – 8.47×** | 2881× | **VIOLATED, ~2× low** |
| 0.4 µM | R_max ∈ [7, 15] | 7.18 – 7.88× | 184.6× | holds, at the bottom edge |
| 3.0 µM | **R_max < 3.0×** | **2.94 – 3.03×** | 4.8× | **VIOLATED at s_max 2.50 (3.03)** |
| 10.0 µM | R_max < 1.5× | 1.34 – 1.35× | 1.5× | holds |

(`c_affinity` = 2.7882 for s_max 2.40, 2.9321 for s_max 2.50.)

**Two of four intervals are violated by the very model they claim to predict.** A lab running
this protocol and measuring ~8× at 0.1 µM would report the model falsified, when the model
predicts ~8×. The draft's intervals look as though they were set from the asymptote column
rather than from a finite-s_max calculation — the asymptote at 0.1 µM really is >15× (2881×),
but no affinity PAM with s_max ≈ 2.5 reaches it.

The underlying reason is structural and worth stating in the paper: at low ambient GABA the
*asymptotic* headroom grows without limit, but the gain a **finite** s_max can extract does
not. Between 0.4 and 0.1 µM the asymptote rises 15.6× while the realised gain rises only 1.08×.
The headroom and the reachable gain diverge, and the falsification test must be written against
the reachable gain.

**Fix:** replace the four intervals with the computed ones, widen them to a defensible
tolerance, and print the asymptote alongside so the distinction is explicit. The falsification
rule as written ("R_PAM ≤ 2.5 at 0.4 µM falsifies Model B") **is** sound and worth keeping — the
repo predicts 7.18–7.88× there, comfortably clear of 2.5.

---

## 8. Required corrections, in priority order

1. **Delete the commit hash or replace it with a real one.** `c8f17a9` does not exist. Tag the
   commit the numbers come from and cite the tag.
2. **Regenerate every number in §2, §2.5 and Table 4** from `scripts/paper_numbers.py`. Paste,
   never retype. Soften "1.062× peak gain" to the computed **1.319×** throughout the abstract,
   §4.1, §5.1 and the Conclusions, and change the headline to **184.6×** / **111× contrast**.
3. **Fix 5 miscited references; resolve or delete 3 unresolvable ones.** The three unresolvable
   ones each carry a parameter; if they cannot be resolved, mark those parameters UNSOURCED in
   the manuscript as the repo already does internally.
4. **Drop the MC median and 95th-percentile columns from Table 6.** The model refuses to stand
   behind them in its own output.
5. **Retitle §4.4** to drop "Multiscale Circuit".
6. **Add Table 5 provenance per cell** from `provenance.py`, and stop attributing unsourced cells
   to Pirker and Kasugai.
7. **Regenerate the §5.2 falsification intervals** — two of four are violated by the model (§7b). This is the paper's falsification test; it must not be wrong.
8. **Add §4.6: substrate independence** (the link-5 table above), and the BZ sign flip as a
   limitation on the reference arm.
9. **Add the four omitted limitations** to §5.3.
10. Fix `po_max` either way *deliberately*: if 0.84 is the better anchor, change `FIT_TARGETS` and
   widen `FIT_RANGES` **in the code**, re-run, and let the paper quote the result. Do not leave
   the paper and the code disagreeing in either direction.
11. Minor: HZ-166 P(R>1) 63.8% not 64.8%; disclose the 10 discarded draws; ">48×" not ">50×" for
    the 0.2–0.8 µM span.

## What should not change

The thesis. Compartment-dependent allosteric headroom is real in this model, the fixed 2.5×
scalar cap is wrong for the extrasynaptic pool, the two traps in §3 are correctly diagnosed, the
scale-invariance argument for R is valid, the stereochemical-invariance result is correct and
well argued, and the proposed patch-clamp protocol in §5.2 is a genuine falsification test with
pre-registered intervals. Those intervals need regenerating from the repo's parameterisation,
but the design is sound and is the most valuable part of the manuscript.
