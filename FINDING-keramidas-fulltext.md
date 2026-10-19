# Finding — Keramidas & Harrison 2008 read in full text: two checks, both consequential

**Read by:** Claude Opus 5, 2026-10-18, **full text** via PubMed Central (PMC2213567).
Previously this project held only the extracted parameter values with
`err_kind='UNSTATED'`. Source: Keramidas A, Harrison NL, *J Gen Physiol* 131(2):163–81, 2008.
PMID 18227274, DOI [10.1085/jgp.200709871](https://doi.org/10.1085/jgp.200709871). Retrieved
from PubMed.

## 1. The unstated error kind is VERIFIED, not assumed

`MANUSCRIPT-DRAFT.md` §3.1 rests the mean-open-time band treatment on the source not stating
whether its errors are SD or SEM. **That is now checked against the full text and holds.**

Searched the complete article for every conventional form: `s.e.m` / `S.E.M.`, "standard error",
"standard deviation", `SD`, `± SE`, `mean ±`, `+/-`. Result: **zero genuine hits.** The four
matches on `s.?e.?m` are the words "assemble" and "homologous"; "standard error", "standard
deviation", `SD`, `mean ±` return **nothing at all**.

**Consequence for the manuscript.** The band treatment is no longer a conservative choice in the
face of ambiguity — it is the only treatment the source supports. Assigning a Gaussian σ would
require inventing an error convention the authors never published. §3.1 should be strengthened
accordingly: the point-treatment alternative (which fails at 3.56σ) is not a legitimate rival
reading, because it requires a σ that does not exist in the source.

This removes one of the manuscript's three stated soft spots.

## 2. The same paper reports THREE open states for α1β2γ2S — our best model has two

Verbatim from the full text:

> "GABA-gated bursts required three shut and three open components to adequately describe M-Mode
> activation (… = 10 patches). This is the same number of components reported for αβγGABARs."

and for the high-P_o mode, from the abstract:

> "the 'high' P_o mode (∼0.9) was characterized by two shut and three open states"

So an independent single-channel study of **the same subunit combination** (α1β2γ2S, outside-out
patches) concludes **three open states**, in both gating modes, and notes that other studies of
αβγ receptors report the same number.

**This cuts two ways and both matter.**

* **Corroboration.** The manuscript's central structural claim is that one conducting open state
  is insufficient and a second is required. That claim was derived here purely from macroscopic
  deactivation amplitudes. An independent line of evidence — single-channel burst analysis, a
  different experiment on the same receptor — independently requires more than one open state.
  **This is the first external corroboration the structural result has had**, and it did not come
  from the data we fitted.
* **Our model is probably still too simple.** We have two conducting open states; this source
  requires three. That is a live candidate explanation for the one tension the manuscript cannot
  resolve — the Hill slope ceiling of 1.25 against a published 2.2 ± 0.4 (§5). Adding open states
  monotonically increased the achievable Hill slope in every test so far.

**What this does not establish.** Keramidas's three open components are *within bursts* at
saturating GABA, obtained by fitting dwell-time histograms; our open states are states in a
macroscopic gating scheme. They are not automatically the same objects, and the mapping between
a burst-level dwell-time decomposition and a whole-receptor kinetic scheme needs stating rather
than assuming. Note also that the two modes differ in their **shut** state count (three for
M-mode, two for H-mode), which is modal gating and is not something any scheme in this project
represents at all.

## 3. Consequences

1. `MANUSCRIPT-DRAFT.md` §3.1 and §8 to be revised: the band treatment is source-mandated, not
   a convention. One soft spot removed.
2. A three-open-state scheme should be built and tested. If the Hill slope ceiling rises toward
   2.2, that resolves the manuscript's residual tension and converts an unexplained deficit into
   a prediction that was met. Brief 18.
3. `circuitpharm.literature` should record that the full text has now been read, and that the
   source reports three open components — the record currently carries only the parameter values.
