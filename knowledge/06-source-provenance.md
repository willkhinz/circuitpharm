# Source provenance audit

**Status: in progress, 2026-10-07.** Mechanical metadata resolution is automated
(`scripts/verify_sources.py`); claim-support assessment is a reading task and is done by hand
below, source by source. Nothing here is marked verified because a DOI resolved.

---

## 1. The hole, stated precisely

The recorded gap was "all 90 sources are marked UNVERIFIED." That understates it. The gap is
one layer lower:

**The numbers that drive the project's one VALIDATED result are not unverified — they are
unsourced.**

`scripts/ranking_robustness.py` computes its score from exactly four inputs:
`PROFILES[key].eff()`, `REGIONS` (per-region subunit fractions), `SUBJECTIVE_WEIGHT`, and
`EXTRASYN`. Of these:

| input | where it lives | source link |
|---|---|---|
| `REGIONS` — per-region α1/α23/α5/δα4/ε fractions | `subtypes.py:34` | **none** |
| `SUBJECTIVE_WEIGHT` — which subtypes carry ethanol's stimulus | `subtypes.py:50` | **none in code** (but see `a5_disc`, §3) |
| `EXTRASYN` — synaptic vs extrasynaptic split | `subtypes.py:73` | **none**; `eps` is labelled `# eps = GUESS` in the source |
| `PROFILES[...].eff()` — per-compound subtype efficacies | `subtypes.py:140+` | per-compound source keys exist |

`subtypes.py` *appears* to reference five source keys. All five —
`alogabat`, `gl_ii_73`, `imepitoin`, `sh053`, `tpa023` — are **compound names that collide
with source keys**. They appear as `PROFILES` dictionary keys, not as citations. The true
count of citations in the module that holds the project's central numbers is **zero**.

The database cannot supply them either. `subunit_expression` has 8 rows and every one is
**qualitative**: "predominant", "present", "high", "enriched", "restricted". There are no
numbers in the knowledge base from which `a1=0.60, a23=0.15, a5=0.02` could have been read.
They were written separately from the sources.

This is the same failure already in the error catalogue — *"scored an invented range as if it
were data, from the phrase 'strongly potentiated'"* — sitting under the one result the project
calls VALIDATED, and never caught by eight review passes because the reviews read code, not
provenance.

---

## 2. What the robustness analysis does and does not rescue

Measured 2026-10-07, before writing anything else here, because if the headline were an
artifact of the prior then nothing else in this document would matter.

The Dirichlet prior is centred on the (unsourced) nominal fractions with `KAPPA = 15`. For
preBötC α5 that is `alpha = 0.3`, which puts **38% of draws below 0.002** — a tenth of
nominal — and gives the forebrain:preBötC α5 ratio a p99 of ~3.5×10⁷. Since an α5-selective
compound's respiratory burden is roughly proportional to that fraction, the score should
diverge, and the whole result could have been prior shape rather than pharmacology.

**It is not, for the statistic the claim rests on.** Flooring the drawn preBötC α5 fraction:

| preBötC α5 floor | alogabat p5 | median | p95 | P(>1) |
|---|---|---|---|---|
| none (as shipped) | **2.51** | 16.05 | 84.9 | 99.9% |
| 0.005 | **2.51** | 13.56 | 36.0 | 99.9% |
| 0.01 | **2.51** | 10.83 | 23.5 | 99.9% |
| 0.02 (= nominal) | **2.51** | 7.63 | 14.5 | 99.9% |

`corr(score, 1/α5) = 0.02`. The 5th-percentile floor and the ordering are robust; my
suspicion that they were artifacts was wrong.

**What IS an artifact is the optimistic end.** The median halves and p95 falls 6× as the floor
rises to nominal. So:

> **QUOTABLE:** the 5th-percentile floor (≥2.49×) and the ordering.
> **NOT QUOTABLE:** the median and 95th percentile. The previously-reported *"median 33.76×,
> p95 9264.80×"* for the ideal α5 arm are properties of the prior's tail. A p95 of 9,264× was
> never a pharmacological claim.

Reproduce with `scripts/ranking_robustness.py --a5-floor 0.02`.

**A defect found in passing.** For a perfectly α5-selective compound the burden can reach
zero, giving an undefined score. The script's `isfinite` filter discards those draws — which
are **exactly the draws most favourable to that arm**. That makes the reported numbers
conservative, the honest direction, but the discard was never counted. It is now reported.

**What remains unrescued:** the prior's *centre* and its *width* (`KAPPA = 15`,
`EXTRASYN_CONC = 8`) are both unsourced choices of mine. "Robust across 20,000 draws" means
robust across a prior I invented. That is a weaker claim than it reads as, and the floor's
insensitivity to the α5 floor does not speak to it.

---

## 2b. Mechanical metadata resolution: all 90 sources

`scripts/verify_sources.py`, run 2026-10-07, written to the `sources.verification` column.
Each row now reads `METADATA_<status> (date; how); CLAIM_SUPPORT_UNASSESSED` — the two are
separate fields and are never collapsed.

| status | n | meaning |
|---|---|---|
| `RESOLVED_MATCH` | **62** | work resolves by DOI/PMID/title and the recorded citation matches |
| `RESOLVED_MISMATCH` | 5 | a work resolves but the recorded title differs — see below |
| `UNRESOLVED` | 12 | nothing resolved; may still be a correct citation |
| `NON_LITERATURE` | 6 | trial registry, encyclopedia, software resource — no DOI by nature |
| `INTERNAL` | 5 | this project's own output |

48 abstracts retrieved and available for the reading half.

**The 5 mismatches are citation-hygiene problems, not fabrications.** In each case the
`citation` field holds a truncation or a paraphrase rather than the verbatim title:
`bk` (prefix of the real title), `girk` (abbreviated), `ganaxolone` (different subtitle —
worth confirming it is the same work), `gl_ii_73` and `alogabat` (a summary of the finding
used in place of a title, similarity 0.35). Real, because a citation that is not the title
cannot be checked mechanically by anyone else — but not invented sources.

**The 12 unresolved** are mostly publisher pages without an embedded identifier
(ScienceDirect), plus one method reference with no URL (`matsuoka`). They need resolving by
hand or by adding DOIs.

### Two bugs in the verification tool itself, both found the same way

Both were caught because sources this audit had *just cited as provenance* came back
UNRESOLVED:

1. **URL-encoded DOIs were invisible.** PLOS writes `?id=10.1371%2Fjournal.pone.0030608`;
   the DOI regex needs a literal slash. `pbc_eps` went UNRESOLVED → RESOLVED_MATCH.
2. **Nature URLs carry the DOI suffix without the prefix.** `/articles/s41598-017-17379-x`
   is `10.1038/s41598-017-17379-x`. `pbc_delta` and `vrodent` recovered.

A verification tool that *under*-reports resolution is not a safe failure. It makes real
citations look missing and so overstates the provenance gap — misleading in the opposite
direction to the one everyone guards against.

---

## 2c. Per-parameter provenance, enforced in code

`src/circuitpharm/provenance.py`, with `tests/test_provenance.py` (12 tests) failing if a
parameter is added without a record, if a record names a source absent from the database, or
if the headline count moves without being updated.

**6 of 28 load-bearing parameters name a source (21%).** 18 UNSOURCED, 8 FROM_QUALITATIVE,
1 FITTED, 1 GUESS.

The worst single entry is `REGIONS["prebotc"]["a5"] = 0.02`. An α5-selective compound's
modelled respiratory burden is roughly proportional to it, so it sets the entire safety
margin — and no source in the knowledge base states it.

---

## 3. Sources assessed by hand

Metadata status comes from `scripts/verify_sources.py`; the verdict column is a reading
judgement and is mine.

### `a5_dist` — **RESOLVES PERFECTLY, DOES NOT SUPPORT THE NUMBERS**

* Recorded: *"Regional distribution of the GABA-A/benzodiazepine receptor (alpha subunit) mRNA
  in rat brain"*, 1988, PMID 2844998.
* Metadata: `RESOLVED_MATCH`, title similarity 1.00. The citation is exactly right.
* Abstract retrieved (1,194 chars) via Europe PMC.

What it actually reports: a **single cDNA probe for "the alpha subunit"** — generic, 1988,
before the α subtypes were separable — used to quantify *total* α-subunit mRNA by region.
Result: *"highest in the cerebellum followed by the thalamus = frontal cortex = hippocampus =
parietal cortex = hypothalamus much greater than pons = striatum = medulla."*

**Verdict: INSUFFICIENT for the attributed use.** It does not resolve α1 vs α2 vs α3 vs α5 at
all, so it cannot support any per-subtype regional fraction. It supports only "total GABA-A α
expression is much lower in medulla than in hippocampus/cortex" — which is a claim about
regional *level*, whereas `REGIONS` encodes regional *composition* (its rows sum to 1 by
construction). The source is therefore not merely too coarse; it is about a different
quantity.

This is the first source audited in the project and it is the reason this document separates
metadata from claim support. A pipeline that wrote "VERIFIED" on a resolved DOI would have
marked this one verified.

### `pbc_alpha` — metadata verified, content inaccessible

* Recorded: *"Developmental changes in expression of GABA-A receptor subunits alpha1, alpha2,
  alpha3 in the rat pre-Botzinger complex"*, 2004.
* Crossref on `10.1152/japplphysiol.01264.2003` returns *"Developmental changes in the
  expression of GABA<sub>A</sub> receptor subunits α1, α2, and α3 in the rat pre-Bötzinger
  complex"*, Liu & Wong-Riley, J Appl Physiol 2004. **Exact match.**
* Full text: HTTP 403. Abstract not in Europe PMC under this identifier.

**Verdict: CLAIM SUPPORT UNASSESSED.** Relevant by title to the preBötC α1/α2/α3 fractions,
and the correct paper to use. Two cautions that must travel with it: it covers α1, α2 and α3
only (not α5, δ, α4 or ε, which `REGIONS["prebotc"]` also assigns numbers to), and it is a
**developmental** study — an axis this project has already been bitten by twice (the
NKCC1/KCC2 chloride switch, and Butera's neonatal parameters).

### `a5_disc` — supports its claim, qualitatively

* *"Contribution of α1GABA-A and α5GABA-A receptor subtypes to the discriminative stimulus
  effects of ethanol in squirrel monkeys"*, 2005, PMID 15650112. `RESOLVED_MATCH`, sim 0.96,
  abstract retrieved.
* Recorded finding: α5 agonists (QH-ii-066, panadiplon) mimic ethanol's discriminative
  stimulus; the α5 inverse agonist L-655,708 blocks it.

**Verdict: SUPPORTS, qualitatively.** This backs `SUBJECTIVE_WEIGHT[a5] = 1.0` as a
*direction* — α5 carries the stimulus. It does not supply a weight, and it says nothing about
`SUBJECTIVE_WEIGHT[a23] = 1.0`, which remains unsourced. Nor does it license zeroing α1,
given the paper's title names an α1 contribution.

**Open question this raises:** `SUBJECTIVE_WEIGHT` sets `a1 = 0.0`, and the one source that
bears on it is a paper about the contribution *of α1 and α5*. That needs the paper read
before the zero stands.

---

## 4. What filling the hole properly requires

1. **Per-parameter provenance, enforced.** Every number in `REGIONS`, `EXTRASYN` and
   `SUBJECTIVE_WEIGHT` carries a provenance record naming its source and the derivation, or is
   explicitly marked unsourced. A test asserts completeness so a new number cannot be added
   without one.
2. **Point values become ranges where only a qualitative statement exists.** "α1-predominant"
   does not mean 0.60; it means something like 0.4–0.8. The ranking already propagates
   uncertainty, so the honest move is to widen the prior to what the literature actually
   licenses rather than to keep an invented centre with an invented width.
3. **Read `pbc_alpha` and `a5_disc` properly** — the two sources that genuinely bear on the
   load-bearing numbers. Both are behind paywalls for full text; abstracts are available for
   one.
4. **Re-run the robustness analysis** under the widened prior and report what survives. The
   floor's insensitivity to the α5 floor is encouraging but does not predict this.

Until (1) and (2) are done, the selectivity ranking's tier should read: *ordering robust to
its own prior's tail; prior itself unsourced.* That is not VALIDATED in the sense
`results.Tier` defines it.

---

## Claim-support pass, 2026-10-08

Prompted by a manuscript review asking that each citation support *the exact statement attached
to it*, not merely a related topic. Method: resolve against Crossref by title and by
journal/volume/page (as before), then retrieve the abstract from Europe PMC and read it against
the specific claim. Reproduce with the queries in `scripts/verify_sources.py`; the readings
below are mine.

**Six load-bearing citations do not support the specific number attached to them.** They fail
in three distinguishable ways.

### Measuring a different quantity than the claim requires

**`EXTRASYN[a5] = 0.80` ← Kasugai et al. 2010.** *Eur. J. Neurosci.* 32:1868–1888. The obvious
citation and the right kind of study: quantitative freeze-fracture replica immunogold labelling,
synaptic against extrasynaptic pools, hippocampal CA1 pyramidal cells. **Its abstract states it
measured α1, α2 and β3.** It does not measure α5, and its quantitative result points the other
way — synaptic labelling density exceeded extrasynaptic by **78–132× (α1), 94× (α2), 79× (β3)**.

→ `EXTRASYN_PROV["a5"]` downgraded **FROM_QUALITATIVE → UNSOURCED**. This widens the audit gap
to 19 UNSOURCED / 7 FROM_QUALITATIVE, pinned in `tests/test_provenance.py`.

**This is the second source to resolve perfectly and fail to support its number** (`a5_dist`
was the first, §3). The pattern is specific enough to name and watch for:

> **A source whose title matches the claim, in the right journal, by the right group, measuring
> a neighbouring quantity.**

Metadata verification is structurally incapable of catching it — `a5_dist` resolved at
similarity 1.00 — and so is a plausibility check, because the citation *is* plausible. Only
reading catches it. Worth noting too that `EXTRASYN[a5]` never carried a source key: it was
FROM_QUALITATIVE with an empty key, so **a basis label was asserting more confidence than any
recorded citation supported.** The audit should probably not permit that state.

**diazepam `ceiling = 2.50` ← Walters et al. 2000.** *Nat. Neurosci.* 3:1274–1281 (itself a
correction; an earlier draft cited *Br. J. Pharmacol.* 131:1307–1314 under a different title).
Reports that diazepam potentiates α1β2γ2 in **two separable phases**, nanomolar and micromolar,
dissociable by TM2 mutations. No EC₅₀ fold shift. This is not weak support for 2.50 — it is
evidence that **a single affinity parameter is an oversimplification for this compound**, which
is the reference arm and therefore the denominator of every selectivity ratio the project
reports.

**regional composition ← `a5_dist` (1988, PMID 2844998).** Previously documented in §3; included
here for completeness.

### Reporting no fold shift of the kind the mapping needs

**MP-III-022 `ceiling` ← Stamenić et al. 2016.** *Eur. J. Pharmacol.* 791:433–443. Strongly
supports the *direction*: a binding- and efficacy-selective α5 PAM, potentiation "from mild to
moderate to strong" over 1–10 mg/kg, non-α5 receptors engaged only at the top dose, no ataxia or
sedation at 1–10 mg/kg. Supplies no EC₅₀ shift. **SUPPORTS DIRECTION ONLY.**

**alogabat `ceiling` ← Cecere et al. 2025.** *Front. Pharmacol.* 16:1626078. Abstract confirms a
potent α5β3γ2 PAM with binding and functional selectivity, plus receptor-occupancy, pMRI and EEG
target engagement. The **+167% rat / +72% human EC₂₀** figures are **not in the abstract** and
were not verified against full text. **PARTIAL.**

### Right quantity, different value

**τ_IPSC = 15.0 ms ← Haas & Macdonald 1999.** *J. Physiol.* 514:27–45 (also a corrected
citation). Supports the topology claim well — deactivation rate is set by subunit composition —
but the measured value for α1β3γ2L is **76.1 ms**, about 5× slower than our anchor. The 15 ms
anchor is a conventional neuronal IPSC decay figure and **is not supported by the source cited
beside it.**

### Partial, not load-bearing on a number

* **ambient GABA 0.2–0.8 µM ← Farrant & Nusser 2005.** Abstract says "low concentrations of
  ambient GABA" without the range; presumably in the review body, unverified here.
* **`INVITRO_BAND` 0.05–1.00 Hz ← Revill et al. 2021.** Confirms the preparation (neonatal rat
  slices retaining respiratory rhythmicity); no frequency range in the abstract.

### Not assessable

Jones & Westbrook 1995, Otis & Mody 1992 and Pirker et al. 2000 are not indexed in Europe PMC
(all pre-date reliable abstract coverage). Metadata verified by Crossref; claim support
**UNASSESSED**. Pirker is cited only for a qualitative ordering and carries no number.
`pbc_alpha` remains UNASSESSED on full-text 403.

### What this does and does not mean

None of these findings makes the underlying statement false. α5 probably *is* enriched
extrasynaptically; diazepam probably *does* shift GABA EC₅₀ by something in the 2–3× range. What
the pass establishes is that **this project cannot cite a source for those numbers** and must
carry them as assumptions. That is a weaker position than a resolved reference list implies, and
an accurate one.

The receptor-kinetic results are unaffected: they rest on the three stated calibration anchors
and the topology, not on these citations. The τ_IPSC finding is the one exception, and it bears
on an anchor rather than a conclusion.
