"""Benchmark datasets. Two are currently SYNTHETIC -- generated from a formula for
mechanics testing, not measured. See `synthetic` and `origin` on each, and roadmap P1-2
for the digitisation that replaces them.

WHAT WENT WRONG HERE, recorded because the defect is subtle and was shipped three commits
after the audit that should have prevented it (roadmap P0-6):

`DEACTIVATION_BENCHMARK` is `0.70*exp(-t/15) + 0.30*exp(-t/70)`, a formula, and it carried
`citation="Haas & Macdonald 1999 ..."`. `knowledge/07-paper-review.md` records that paper
as measuring **76.1 ms** for this quantity -- it is one of the six load-bearing citations
that claim-support audit found did not support the number attached to them. So the
generated trace's dominant 15 ms component was attributed to a paper that measured
something five times slower.

`DOSE_RESPONSE_BENCHMARK` has the same problem in a subtler form: its plateau is exactly
0.750, which is `gabaa_kinetics.FIT_TARGETS["po_max"]` -- the model's own fit target. So
"fitting the model to the data" was partly circular, and its `sem` values have no stated
origin at all.

THE FIX IS NOT A RENAME. Both datasets keep their names and stay usable, so no caller
breaks and the mechanics tests still have something to run against. What changed is that

  * `synthetic=True` and `origin` say what they are, in the object itself;
  * the citation text survives under `motivated_by`, a field whose name cannot be mistaken
    for provenance, while `source_key` stays empty;
  * `assert_real_data()` makes any fit that consumes them warn and taint its conclusions
    automatically -- so when P1-2 lands real digitisations the flag flips to False, the
    warning stops on its own, and no call site changes.

A synthetic trace may exist and may be used. It may never carry provenance it does not
have, and it may never be the basis of a claim about receptors (roadmap §5.4).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Sequence

import numpy as np

from ..models.base import Observable


#: HOW A DATASET'S NUMBERS CAME TO EXIST. Three categories, not two, because the honest
#: middle case is the common one and collapsing it either way misleads:
#:
#:   "digitised"  individual data points read off a published figure. The gold standard,
#:                and what roadmap P1-2 asks for.
#:   "parametric" generated from the FITTED PARAMETERS a paper published (EC50, Hill slope,
#:                Emax). The parameters are real measurements and are attributable; the
#:                individual points are the paper's own model of its data, not its data.
#:                Weaker than a digitisation and far stronger than a guess.
#:   "synthetic"  invented. For exercising machinery only; no claim about receptors may
#:                rest on it (roadmap §5.4).
DataKind = Literal["digitised", "parametric", "synthetic"]

#: WHAT THE RESPONSE COLUMN IS, which is part of the observable and was missing from it.
#:
#:   "absolute"        an open probability, or a current already on an absolute scale. A
#:                     model's prediction is comparable to it directly.
#:   "fraction_of_max" I/I_max -- the near-universal convention for a published
#:                     concentration-response, and NOT the same quantity. Its asymptote is
#:                     1 by construction, while this scheme's absolute peak open
#:                     probability saturates near 0.75, so comparing the two directly asks
#:                     the model to reach a maximum it cannot have. Measured: fitting the
#:                     Jahn curve that way sent D to 1e-11 -- the optimiser deleted
#:                     desensitisation, because raising the absolute plateau to 1 is the
#:                     only way to meet a normalised one. That is the P1 observable
#:                     mismatch (roadmap §2.3) in a second guise: PEAK against EQUILIBRIUM
#:                     was caught, PEAK against NORMALISED PEAK was not.
#:
#: A normalised dataset carries NO information about absolute open probability, so a fit
#: against one must not be allowed to claim any. `likelihood._predict` normalises the
#: prediction the same way the data was normalised, which confines the dataset's evidence
#: to the curve's shape -- its EC50 and slope -- where it belongs.
Normalisation = Literal["absolute", "fraction_of_max"]


@dataclass(frozen=True)
class DoseResponseDataset:
    """Concentration-response measurements, or a stand-in for them. See `kind`."""

    citation_label: str          # display text only; see `source_key` for provenance
    preparation: str
    concs_um: np.ndarray
    mean_response: np.ndarray
    sem: np.ndarray
    observable: Observable = Observable.PEAK
    #: What the response column IS: an absolute open probability, or I/I_max. See
    #: `Normalisation` -- getting this wrong does not bias a fit, it breaks it.
    normalisation: Normalisation = "absolute"
    #: How these numbers came to exist. See `DataKind`.
    kind: DataKind = "digitised"
    #: True = NOT measured point-by-point, so it taints any fit that uses it. Kept as a
    #: field (rather than derived from `kind`) because it is what the taint check reads and
    #: what every existing test and caller asks for; `__post_init__` keeps the two
    #: consistent so they cannot disagree.
    synthetic: bool = False
    #: For a synthetic dataset: the exact expression and why it exists.
    origin: str = ""
    #: A paper that INSPIRED the shape. NOT provenance. Never set this for real data --
    #: real data gets `source_key`.
    motivated_by: str = ""
    #: A key into `provenance`. Set ONLY for genuinely digitised or published data.
    source_key: str = ""
    #: Figure/table the numbers came from, e.g. "Fig 2B, filled circles".
    figure: str = ""
    #: Tool, method, date and who, for a real digitisation.
    digitisation: str = ""
    n_cells: int | None = None
    role: Literal["train", "holdout"] = "train"

    def __post_init__(self):
        _validate_dataset(self)

    @property
    def key(self) -> str:
        return self.citation_label


@dataclass(frozen=True)
class DeactivationDataset:
    """A macroscopic deactivation current trace, or a synthetic stand-in for one."""

    citation_label: str
    pulse_duration_ms: float
    gaba_conc_um: float
    time_ms: np.ndarray
    normalized_current: np.ndarray
    observable: Observable = Observable.CHARGE
    #: What the response column IS: an absolute open probability, or I/I_max. See
    #: `Normalisation` -- getting this wrong does not bias a fit, it breaks it.
    normalisation: Normalisation = "absolute"
    kind: DataKind = "digitised"
    synthetic: bool = False
    origin: str = ""
    motivated_by: str = ""
    source_key: str = ""
    figure: str = ""
    digitisation: str = ""
    n_cells: int | None = None
    role: Literal["train", "holdout"] = "train"

    def __post_init__(self):
        _validate_dataset(self)

    @property
    def key(self) -> str:
        return self.citation_label


def _validate_dataset(ds) -> None:
    """A dataset cannot be both generated and sourced, and must say which it is."""
    if ds.kind not in ("digitised", "parametric", "synthetic"):
        raise ValueError(f"{ds.citation_label!r}: unknown kind {ds.kind!r}")

    if ds.normalisation not in ("absolute", "fraction_of_max"):
        raise ValueError(
            f"{ds.citation_label!r}: unknown normalisation {ds.normalisation!r}; see "
            f"`Normalisation`. It must be stated, because a model prediction is comparable "
            f"to one of the two and not the other.")
    peak_response = float(np.max(np.asarray(
        getattr(ds, "mean_response", getattr(ds, "normalized_current", [0.0])),
        dtype=float)))
    if ds.normalisation == "fraction_of_max" and not (0.9 <= peak_response <= 1.0 + 1e-9):
        raise ValueError(
            f"{ds.citation_label!r} declares normalisation='fraction_of_max' but its "
            f"largest response is {peak_response:.4g}. An I/I_max column reaches 1 at the "
            f"concentration it was normalised by; if this one does not, either it is not "
            f"normalised or the normalising point is missing from the dataset, and a fit "
            f"would divide the prediction by the wrong concentration.")

    # `kind` and `synthetic` must agree. "parametric" counts as synthetic for the purpose
    # of tainting a fit -- the individual points are a model, not measurements -- while
    # still carrying a real source_key, which is why the two fields both exist.
    expect_synthetic = ds.kind in ("synthetic", "parametric")
    if bool(ds.synthetic) != expect_synthetic:
        raise ValueError(
            f"{ds.citation_label!r}: kind={ds.kind!r} implies synthetic="
            f"{expect_synthetic} but synthetic={ds.synthetic!r}. A parametric "
            f"reconstruction taints a fit exactly as an invented trace does -- its points "
            f"are the paper's model of its data -- while still naming a real source.")

    if ds.kind == "parametric":
        if not ds.source_key:
            raise ValueError(
                f"{ds.citation_label!r} is parametric but names no `source_key`. The whole "
                f"point of the category is that the PARAMETERS are attributable; without a "
                f"source it is just a synthetic trace.")
        if not ds.origin.strip():
            raise ValueError(
                f"{ds.citation_label!r} is parametric but has no `origin` stating which "
                f"published parameters were used and how the curve was generated.")
        return

    if ds.synthetic:
        if not ds.origin.strip():
            raise ValueError(
                f"{ds.citation_label!r} is marked synthetic but has no `origin`. State the "
                f"expression that generated it and why it exists, or it is indistinguishable "
                f"from data.")
        if ds.source_key:
            raise ValueError(
                f"{ds.citation_label!r} is marked synthetic AND carries source_key="
                f"{ds.source_key!r}. A wholly generated trace has no provenance; put the "
                f"paper that inspired it in `motivated_by`, or -- if its PARAMETERS are "
                f"genuinely published -- declare kind='parametric' instead.")
    else:
        if not ds.source_key:
            raise ValueError(
                f"{ds.citation_label!r} claims to be real data (synthetic=False) but names "
                f"no `source_key`. Register it in `provenance` first; roadmap §2.2.")
        if not ds.figure.strip() or not ds.digitisation.strip():
            raise ValueError(
                f"{ds.citation_label!r} is real data but does not say which `figure` it came "
                f"from or how it was obtained (`digitisation`). Roadmap P1-2 requires both.")


# ---------------------------------------------------------------------------------------
# 1. Concentration-response, PEAK observable.
#
# SYNTHETIC. Shaped to resemble a published alpha1beta2gamma2 macroscopic peak-current
# curve (EC50 ~20-30 uM, Hill ~1.3-1.5) but generated, not digitised. Its plateau is
# exactly the project's own po_max fit target, which is why a fit against it is partly
# circular -- see the module docstring.
# ---------------------------------------------------------------------------------------
_CONCS = np.array([0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0], dtype=float)
_RESP = np.array([0.005, 0.021, 0.075, 0.245, 0.540, 0.710, 0.745, 0.750, 0.750],
                 dtype=float)
_SEM = np.array([0.002, 0.005, 0.012, 0.025, 0.030, 0.020, 0.015, 0.010, 0.010],
                dtype=float)

DOSE_RESPONSE_BENCHMARK = DoseResponseDataset(
    citation_label="synthetic alpha1beta2gamma2 peak concentration-response",
    preparation="nominally recombinant alpha1beta2gamma2, whole-cell; NOT a real recording",
    concs_um=_CONCS,
    mean_response=_RESP,
    sem=_SEM,
    observable=Observable.PEAK,
    kind="synthetic",
    synthetic=True,
    origin=("hand-written response and SEM arrays shaped to resemble a sigmoid with "
            "EC50 ~20-30 uM and a 0.750 plateau. The plateau is exactly "
            "gabaa_kinetics.FIT_TARGETS['po_max'], so fitting the model to this curve is "
            "partly circular; the SEM values have no origin at all. Exists so the "
            "likelihood and identifiability machinery can be exercised before P1-2 "
            "supplies a digitisation."),
    motivated_by=("Mortensen et al. 2012 (J Physiol 590:69-78) and Sigel & Steinmann 2012 "
                  "report curves of this general shape. NEITHER was digitised and neither "
                  "supports these specific numbers."),
)

# ---------------------------------------------------------------------------------------
# 2. Macroscopic deactivation trace, CHARGE observable.
#
# SYNTHETIC, and the one that carried a mis-citation. See the module docstring.
# ---------------------------------------------------------------------------------------
_T_DEACT = np.linspace(0.0, 100.0, 50, dtype=float)
_I_DEACT = 0.70 * np.exp(-_T_DEACT / 15.0) + 0.30 * np.exp(-_T_DEACT / 70.0)

DEACTIVATION_BENCHMARK = DeactivationDataset(
    citation_label="synthetic biexponential deactivation trace",
    pulse_duration_ms=1.0,
    gaba_conc_um=1000.0,
    time_ms=_T_DEACT,
    normalized_current=_I_DEACT,
    observable=Observable.CHARGE,
    kind="synthetic",
    synthetic=True,
    origin=("0.70*exp(-t/15) + 0.30*exp(-t/70), evaluated on 50 points over 100 ms. A "
            "formula, not a digitisation. Exists so the kinetic-timescale half of the "
            "likelihood can be exercised before P1-2; it cannot make any rate "
            "identifiable in a way that means anything, because the timescales in it were "
            "chosen rather than measured."),
    motivated_by=("Haas & Macdonald 1999 (J Physiol 514:27-45) and Jones & Westbrook 1995 "
                  "measure macroscopic deactivation. IMPORTANT: "
                  "knowledge/07-paper-review.md records Haas & Macdonald 1999 as measuring "
                  "76.1 ms for this quantity, so the 15 ms component above is NOT that "
                  "paper's number and must not be attributed to it."),
)


# ---------------------------------------------------------------------------------------
# 3. PARAMETRIC concentration-response from a named source. PEAK observable.
#
# THE FIRST SOURCED KINETIC DATASET IN THIS PROJECT. Until now every target in
# `gabaa_kinetics.FIT_TARGETS` came from what `config.CALIBRATIONS` calls "recalled
# literature RANGES, not digitised from specific published figures", which is why
# `gabaa_kinetics_baseline` is rated UNCALIBRATED and why "no quantitative claim here is
# citable".
#
# Jahn et al. 1997 measured exactly this preparation -- alpha1beta2gamma2L transiently
# expressed in HEK293, patch clamp with ultra-fast solution exchange -- and published the
# fitted Hill parameters. The curve below is generated FROM THOSE PARAMETERS, so it is
# `kind="parametric"`: the EC50 and the slope are real attributable measurements, the
# individual points are the paper's own model of its data. That is weaker than a
# digitisation and far stronger than the synthetic curve above, and it taints a fit exactly
# as a synthetic one does.
#
# TWO THINGS THIS SURFACES, both recorded rather than smoothed over:
#
#   1. THE PROJECT'S EC50 ANCHOR IS NOT THE SOURCED VALUE. FIT_TARGETS['ec50_um'] is
#      20.0 uM, the middle of a recalled 10-30 uM range. The one sourced value for this
#      preparation is 11.6 +/- 0.9 uM. 20.0 is inside FIT_RANGES' accepted [10, 30] band,
#      so nothing was violated -- but refitting to 11.6 would move every derived kinetic
#      number again, and that is P4's job, not a quiet substitution here.
#
#   2. THE ABSTRACT DOES NOT GIVE A MEAN OPEN TIME. It reports a BURST duration of
#      10.3 +/- 3.0 ms, and a burst contains several openings separated by brief closures.
#      So this source does NOT determine `gabaa_kinetics.FIT_FIXED_ALPHA`, which stays a
#      CONVENTION. Getting the mean open time needs the full text, which this environment
#      cannot reach (pubmed, PMC and EuropePMC are all egress-blocked).
#
# THE HILL SLOPE IS NOT WHAT THIS PROJECT FIRST TOOK IT FOR, and getting that wrong
# produced a "structural conflict" that has now been withdrawn. The abstract says, verbatim:
#
#     "The slope between 0.001 and 0.01 mM GABA was 2.2 +/- 0.4, indicating at least
#      three binding sites for GABA."
#
# So 2.2 is a LOCAL slope over 1-10 uM -- the rising phase, 0.086-0.862 x their EC50 of
# 11.6 uM -- and NOT the Hill coefficient of a fit to the whole curve. The distinction is
# invisible for a Hill curve, which is a straight line in logit space and so has the same
# slope everywhere, and it is large for a receptor scheme, whose logit curve BENDS between
# a low-agonist limiting slope equal to the number of binding sites and a flatter slope
# through EC50. Measured on the current fit:
#
#     measurement                                      EC50        slope    gap to 2.2
#     scheme, PEAK, whole-curve regression            19.946 uM     1.294     2.3 sigma
#     scheme, PEAK, rising phase, same RELATIVE part  19.946 uM     1.594     1.5 sigma
#     scheme, PEAK, rising phase, same ABSOLUTE 1-10  19.946 uM     1.711     1.2 sigma
#     scheme, EQUILIBRIUM, whole-curve                 5.222 uM     1.654
#     Jahn 1997, peak current, over 1-10 uM           11.6   uM     2.2 +/- 0.4
#
# Two readings of "the same way" are possible because the scheme's peak EC50 is not Jahn's
# -- the same relative part of the curve, or the same absolute concentrations -- and both
# are given rather than the flattering one. Either way the gap is INSIDE the published
# error bar where the whole-curve comparison put it outside. Three independent checks say
# the same thing:
#
#   * published WHOLE-CURVE Hill fits for a1b2g2 peak currents sit at 1.3-1.6 (e.g.
#     1.5 +/- 0.09, EC50 36 +/- 6 uM, HEK293) -- which is where this scheme sits, not
#     where 2.2 is;
#   * a two-site scheme's limiting log-log slope at low agonist is exactly 2, and its
#     steepest rising-phase fit is ~2.0, within 0.5 sigma of the published 2.2 -- so the
#     slope does NOT establish "at least three binding sites", and cryo-EM of the synaptic
#     a1b2g2 receptor shows TWO GABA sites (Nature 2018);
#   * O. P. Hamill published a comment on exactly this inference in the same NeuroReport
#     issue ("How many transmitter binding steps are involved in opening fast
#     receptor-gated channels?", PMID 9480006), which is not read here but is the right
#     place to look next.
#
# See `knowledge/12-inference.md` Section 2 and `models.base.hill_slope`, which takes the
# measurement window as an argument so this comparison cannot be made loosely again.
#
# WHAT THAT COSTS THIS DATASET, and it is not small. The nine points below are a GLOBAL
# Hill curve with nH = 2.2 spanning 0.3-3000 uM: four decades of curve shape extrapolated
# from a slope the paper measured over one decade. The source constrains the curve in three
# places -- the slope over 1-10 uM, the EC50, and saturation by 3 mM -- and says nothing
# about its shape between 10 uM and 3 mM. Fitting a model to the extrapolated region is
# fitting to an assumption, and that is what drove `log10_D` onto its bound and the plateau
# to 0.99 in `knowledge/12-inference.md` Section 2.1. The `sem` below now says so: inside
# the measured window it propagates the published parameter errors, and OUTSIDE it widens
# to the spread over every slope consistent with the published whole-curve fits for this
# receptor. The dataset stays usable -- it is still the only sourced concentration-response
# here -- and it no longer asserts precision the abstract does not carry.
#
# The EC50 must still not be quietly substituted for the 20 uM target (see point 1 above).
# ---------------------------------------------------------------------------------------
_JAHN_EC50_UM = 11.6
_JAHN_EC50_SEM = 0.9
_JAHN_HILL = 2.2
_JAHN_HILL_SEM = 0.4
_JAHN_CONCS = np.array([0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0],
                       dtype=float)
_JAHN_X = (_JAHN_CONCS / _JAHN_EC50_UM) ** _JAHN_HILL
# A pure Hill curve, asymptote 1. At the paper's 3 mM saturating concentration it
# reaches 0.999995 -- approached, not attained -- which is what 'saturates' means and
# which keeps both published parameters exactly recoverable from the points by a
# logit-space regression (asserted in tests/test_provisional.py).
_JAHN_RESP = _JAHN_X / (1.0 + _JAHN_X)

# Propagated 1-sigma spread from the published parameter errors INSIDE the window the
# slope was measured over, and the spread over admissible curve shapes OUTSIDE it. NOT a
# measured SEM per point -- the paper's per-point errors are in the figure, which was not
# digitised -- so it is an uncertainty on the CURVE, floored at 0.01 so no point claims
# more precision than a normalised macroscopic current can carry.
#
# The two regimes exist because the source constrains the two regions differently. Over
# 1-10 uM it reports a slope and an error, so the curve there is known to +/- that error.
# Everywhere else the abstract gives only the EC50 and saturation by 3 mM, so the honest
# spread is over every slope consistent with published whole-curve Hill fits for this
# receptor (1.3-1.6) through the steepest the local measurement allows (2.6). A fit then
# draws its shape information from the decade that was actually measured.
_JAHN_SLOPE_WINDOW_UM = (1.0, 10.0)      # "between 0.001 and 0.01 mM GABA", verbatim
_JAHN_OUTSIDE_HILL_RANGE = (1.3, 2.6)    # published whole-curve fits .. local upper bound


def _jahn_spread() -> np.ndarray:
    def half_range(hills) -> np.ndarray:
        lo_hi = []
        for ec50 in (_JAHN_EC50_UM - _JAHN_EC50_SEM, _JAHN_EC50_UM + _JAHN_EC50_SEM):
            for nh in hills:
                x = (_JAHN_CONCS / ec50) ** nh
                lo_hi.append(x / (1.0 + x))
        stack = np.vstack(lo_hi)
        return (stack.max(axis=0) - stack.min(axis=0)) / 2.0

    inside = half_range((_JAHN_HILL - _JAHN_HILL_SEM, _JAHN_HILL + _JAHN_HILL_SEM))
    outside = half_range(_JAHN_OUTSIDE_HILL_RANGE)
    lo, hi = _JAHN_SLOPE_WINDOW_UM
    measured = (_JAHN_CONCS >= lo) & (_JAHN_CONCS <= hi)
    return np.maximum(np.where(measured, inside, outside), 0.01)


JAHN1997_PEAK_CRC = DoseResponseDataset(
    citation_label="Jahn 1997 alpha1beta2gamma2L peak GABA concentration-response",
    preparation=("alpha1beta2gamma2L transiently expressed in HEK293; patch clamp with "
                 "ultra-fast solution exchange; peak current"),
    concs_um=_JAHN_CONCS,
    mean_response=_JAHN_RESP,
    sem=_jahn_spread(),
    observable=Observable.PEAK,
    # I/I_max. The paper's own Hill fit has a unit asymptote, which is what a published
    # concentration-response reports; this scheme's ABSOLUTE peak open probability
    # saturates near 0.75. See `Normalisation` for what comparing them directly did.
    normalisation="fraction_of_max",
    kind="parametric",
    synthetic=True,          # see DataKind: parametric taints a fit as synthetic does
    source_key="jahn_a1b2g2_kinetics",
    figure="dose-response relation reported in the abstract (EC50 and Hill slope)",
    digitisation=("NOT digitised. Generated from the paper's own fitted Hill parameters "
                  "EC50 = 11.6 +/- 0.9 uM and nH = 2.2 +/- 0.4 (slope measured over "
                  "1-10 uM). A pure Hill curve, asymptote 1, reaching 0.999995 at the "
                  "3 mM where the paper reports saturation. `sem` is the half-range of the curve over those parameter "
                  "errors, floored at 0.01 -- an uncertainty on the CURVE, not a measured "
                  "per-point SEM. Values taken from the abstract; the full text is "
                  "unreachable from this environment (pubmed/PMC/EuropePMC egress-blocked) "
                  "so the figure was never seen. 2026-10-09."),
    origin=("Hill curve from the published EC50 and slope; see `digitisation`. The paper's "
            "measured quantities are the two parameters, not these nine points."),
    motivated_by="",
    n_cells=None,
    role="train",
)

#: What is still missing, named so it cannot be forgotten (roadmap P1-2, P3).
#:
#: P3's second anchor test -- "k_off becomes identifiable once a kinetic dataset is in the
#: likelihood" -- needs a real macroscopic DEACTIVATION trace, and there is none here. A
#: search surfaced candidate values for alpha1beta2gamma2 (a weighted deactivation tau of
#: ~5.9 ms in one rapid-application macropatch study; ~9 ms fast / ~148 ms slow in
#: another) but none could be verified: every full text is behind this environment's
#: egress proxy, and this project has twice been bitten by a source that resolves
#: perfectly and measures a neighbouring quantity (knowledge/06-source-provenance.md).
#: Attaching one of those numbers to a citation I have not read against the claim is
#: exactly forbidden pattern 5.5, so the gap is recorded instead of filled.
# WHY THESE ARE STILL MISSING, re-checked 2026-10-10 and the answer has CHANGED AGAIN.
#
# 2026-10-11 (cloud): every document below was located and is OPEN ACCESS; what blocked
# them was an egress policy. pmc.ncbi.nlm.nih.gov, europepmc.org, rupress.org and
# discovery.ucl.ac.uk all returned "CONNECT tunnel failed, response 403", and WebFetch
# failed DNS for every host. Only abstracts could be read, via web search.
#
# 2026-10-10 (local): pmc.ncbi.nlm.nih.gov IS reachable and serves full text. All SEVEN
# open-access targets were read in full. europepmc.org, rupress.org and discovery.ucl.ac.uk
# still return 403 (confirmed with a browser user agent, so it is a real block, not UA
# sniffing) -- but every target had a PMC copy, so that no longer matters.
#
# THE GAPS ARE THEREFORE NO LONGER "UNREAD". They are now characterised, and two of the
# three turn out to be harder than "find the number", because the number is not unique:
#
#   * DEACTIVATION is protocol-dependent, measured. Barberis 2007 reports the SAME receptor
#     in the SAME study at 52.5 ms (2 ms pulse) and 364 ms (3 s pulse) -- 6.9x from pulse
#     duration alone. Dixon 2014's 5.9 ms differs from Barberis's 2 ms-pulse 52.5 ms by
#     8.9x at nearly the same pulse duration, because Dixon fitted TWO exponentials where
#     Barberis fitted THREE and Barberis's 221 ms component carries 44.3 of its 52.5 ms.
#   * MEAN OPEN TIME is definition-dependent, spanning 1.42-7.25 ms across sources, and no
#     source states one for alpha1beta2gamma2L at all.
#
# This is the same failure class as PEAK-vs-EQUILIBRIUM (P1), absolute-vs-normalised (P4)
# and PEAK's application dependence (P6): same tag, same units, different protocol. It is
# recorded rather than averaged away. scripts/literature_gamma2.py holds every reading with
# its conditions and PRINTS the derived quantities; knowledge/14-literature.md is the
# reading. Nothing below has been entered as a dataset yet -- see each entry for why.
#
# GRADING, because the distinction is the whole point of this project:
#   [VERIFIED]   read here in the FULL TEXT at the stated URL, with conditions attached.
#   [ABSTRACT]   read here, but only the PubMed abstract exists / was reachable.
#   [LEAD]       reported by a literature-search agent and NOT independently confirmed.
#   [UNREAD]     located, content not obtained.
#
# TWO CITATION ERRORS IN THE PREVIOUS VERSION OF THIS BLOCK, corrected below: Dixon 2014 is
# J Biol Chem 289(9), not 289(8); Barberis 2007 is Eur J Neurosci 25(9):2726-2740, not
# 26(7). Both were carried from agent-supplied metadata.
MISSING_DATASETS = {
    "deactivation_peak_pulse": (
        "A macroscopic deactivation time course for alpha1beta2gamma2 after a brief "
        "saturating GABA pulse. This is the dataset that makes an absolute rate (k_off) "
        "identifiable at all; without it the fit determines only ratios. "
        "STATUS: the numbers now exist and the gap has changed character -- what is "
        "missing is a TIME COURSE, and what was found is weighted scalars plus one full "
        "component set. "
        "[VERIFIED] Dixon, Sah, Lynch & Keramidas 2014, J Biol Chem 289(9):5399-5411 "
        "(PMID 24425869, PMC3937617): weighted deactivation tau = 5.9 +/- 0.5 ms (n = 10) "
        "for human alpha1beta2gamma2L, HEK293 outside-out, room temperature (no value "
        "given), <=1 ms application of 3 mM GABA, -70 mV. The previously [LEAD] value is "
        "CONFIRMED VERBATIM. Error statistic is NOT declared by the paper. The individual "
        "components behind that weighted tau are NOT in the full text -- the only 'two "
        "exponential' fit mentioned is of the SIMULATED current. "
        "[VERIFIED] Barberis et al. 2007 (PMC1950087) gives a FULLY SPECIFIED triple "
        "exponential for rat alpha1beta2gamma2S: tau = 2.8 +/- 0.3, 33.4 +/- 4.6, 221.35 "
        "+/- 14.9 ms with areas 0.57 +/- 0.04, 0.23 +/- 0.02, 0.20 +/- 0.03 (SEM, n = 6), "
        "2 ms pulse of 10 mM GABA, 22-24 C, tau_w = 52.5 +/- 2.9 ms. Recomputing "
        "sum(A_i tau_i) gives 53.55 ms, 0.36 SEM from the published tau_w -- so this is "
        "the ONE observable found that can be reconstructed from parameters the source "
        "itself fitted, with no extrapolation. It is gamma2S, and it disagrees with Dixon "
        "by 8.9x for the reasons in the header. "
        "NOT ENTERED because using either requires deciding WHICH quantity the project's "
        "deactivation observable is: pulse duration, component count and fit window must "
        "all be declared first, exactly as PEAK_APPLICATION_MS had to be."),
    "single_channel_mean_open_time": (
        "Mean open time for alpha1beta2gamma2L. Would replace "
        "gabaa_kinetics.FIT_FIXED_ALPHA -- currently a CONVENTION -- with a measurement. "
        "STATUS: NO SOURCE STATES ONE FOR gamma2L. Four readings exist, spanning 5.1x, and "
        "the spread is DEFINITION rather than disagreement. "
        "[ABSTRACT] Jahn et al. 1997, NeuroReport 8(16):3443-3446 (PMID 9427304, no PMC "
        "copy, paywalled): the abstract contains NO mean open time, only BURST duration "
        "10.3 +/- 3.0 ms. Whether the full text states one REMAINS UNKNOWN -- this is as "
        "far as the abstract can settle it. "
        "[VERIFIED, DERIVED] Keramidas & Harrison 2008 (PMC2213567), rat "
        "alpha1beta2gamma2S, 10 mM GABA, HEK293 excised outside-out at 21 +/- 1 C: Table "
        "III gives tau_O = 0.49/2.58/5.17 ms with areas 0.26/0.49/0.25 (M-mode, n = 10) "
        "and 0.59/4.22/13.2 ms with areas 0.18/0.41/0.41 (H-mode, n = 5). Area-weighted "
        "means 2.68 and 7.25 ms are COMPUTED, not stated. Error statistic NOT declared. "
        "[VERIFIED, DERIVED] Li et al. 2008 (PMC2241790), rat alpha1beta2gamma2L -- the "
        "right splice variant -- 50 uM GABA (SUB-SATURATING), HEK293 CELL-ATTACHED, room "
        "temperature: Table 1 gives OT = 0.28 +/- 0.05 / 3.0 +/- 0.7 / 7.3 +/- 3.2 ms "
        "(s.d., n = 4) with fractions 0.22/0.65/0.13; area-weighted mean 2.96 ms, COMPUTED. "
        "NOTE the control rows are reproduced from Li et al. 2007b (J Physiol 584:789-800), "
        "so cite THAT for these open times. "
        "[VERIFIED] Barberis et al. 2007 (PMC1950087) is the only source to STATE a mean "
        "open time for alpha1beta2gamma2: 1.42 +/- 0.05 ms (SEM, n = 6) -- but over a 4 s "
        "window AFTER a 2 ms pulse, so NON-STATIONARY, mixing early high-occupancy "
        "openings with late brief singly-bound ones. Its steady-state intraburst "
        "distribution is in Figure 6C and could not be read. "
        "NOT ENTERED because 'mean open time' is not one quantity: 1.42 "
        "(non-stationary), 2.68 (intraburst M-mode), 2.96 (intracluster, sub-saturating), "
        "7.25 (intraburst H-mode). FIT_FIXED_ALPHA needs a stated choice among these, and "
        "the two intraburst/intracluster readings nearest the project's regime agree at "
        "2.68 and 2.96 ms across BOTH splice variants, which is the useful fact."),
    "holdout": (
        "Any observable not used in fitting. P5's cross-validation and P6's falsification "
        "bound both consume it, so until one exists neither can be scored out-of-sample. "
        "STATUS: the target was read and is PARTLY unusable, with one usable alternative "
        "inside the same paper. "
        "[VERIFIED] Barberis et al. 2007, Eur J Neurosci 25(9):2726-2740 (PMID 17561840, "
        "PMC1950087), rat alpha1beta2gamma2S, HEK293 outside-out pooled with small lifted "
        "whole cells, 22-24 C, ultrafast exchange (60-100 us), 3 s pulses of 10 mM GABA. "
        "Desensitisation onset: tau_1 = 2.9 +/- 0.1 ms, A_1 = 0.56 +/- 0.025, steady-state "
        "weight 0.076 +/- 0.013, and steady-state:peak measured at 200 ms = 0.21 +/- 0.02. "
        "BUT tau_2 and tau_3 are reported ONLY in Figure 4C, so the onset time course "
        "CANNOT be reconstructed -- and the source is INTERNALLY INCONSISTENT about n for "
        "this measurement: the text says n = 7, the Figure 4 caption says three patches "
        "for alpha1beta2gamma2. That ambiguity is recorded, not resolved. "
        "[VERIFIED] BETTER HOLDOUT, same paper, fully stated as a scalar: paired-pulse "
        "recovery at a 100 ms gap = 0.33 +/- 0.03 (Table 1; Fig. 5 caption gives n = 8). "
        "This is a stated interval with a stated value and needs no reconstruction, which "
        "makes it the cheapest honest holdout available. "
        "NOT ENTERED pending the protocol-declaration decision above, since a paired-pulse "
        "observable must declare its pulse duration (2 ms here) and gap."),
}


#: Every registered dataset, for the provenance test to iterate.
ALL: dict[str, object] = {
    "dose_response_peak": DOSE_RESPONSE_BENCHMARK,
    "deactivation_charge": DEACTIVATION_BENCHMARK,
    "jahn1997_peak_crc": JAHN1997_PEAK_CRC,
}

TRAIN = tuple(d for d in ALL.values() if getattr(d, "role", "train") == "train")
HOLDOUT = tuple(d for d in ALL.values() if getattr(d, "role", "train") == "holdout")


def assert_real_data(datasets: Sequence[object], *, caller: str) -> list[str]:
    """Warn, loudly and by name, for every synthetic dataset a fit is about to consume.

    Returns the list of synthetic dataset labels, so the caller can decide what to taint.
    A helper rather than a copied five-line block, because "forgot it in one of five
    places" is this project's recurring error E12 (roadmap §5.3).
    """
    from ..provisional import warn_provisional

    fake = [getattr(d, "citation_label", repr(d)) for d in datasets
            if getattr(d, "synthetic", False)]
    if fake:
        warn_provisional(
            what=f"{caller} (fitting against {len(fake)} synthetic dataset(s))",
            defect=("the following datasets are generated from formulae, not measured: "
                    + "; ".join(fake)
                    + ". A fit against them measures the formula, not a receptor."),
            register_item="P0-6 / P1-2",
            promote_by=("digitise the concentration-response and deactivation curves from "
                        "named published figures, register them in provenance.py, and set "
                        "synthetic=False -- this warning then stops on its own"),
        )
    return fake


def holdout_guard(datasets: Sequence[object], *, caller: str) -> None:
    """Refuse to fit against a held-out dataset.

    P5's cross-validation and P6's falsification both consume the holdout, and a fit that
    has seen it makes both meaningless. Enforced here rather than remembered.
    """
    leaked = [getattr(d, "citation_label", repr(d)) for d in datasets
              if getattr(d, "role", "train") == "holdout"]
    if leaked:
        raise ValueError(
            f"{caller} was passed held-out dataset(s) {leaked}. The holdout exists so P5 "
            f"can score out-of-sample and P6 can state a falsification bound; fitting "
            f"against it destroys both. Pass `fitting.data.TRAIN`.")


def peak_crc_from_model(model, *, concs_um=None, noise_sd: float = 0.0, seed: int = 0,
                        label: str = "method check",
                        normalisation: str = "absolute") -> DoseResponseDataset:
    """A PEAK concentration-response generated from a model at known parameters.

    The PEAK twin of `equilibrium_crc_from_model`, and the same single justification: a
    method check, never evidence about receptors. PEAK rather than EQUILIBRIUM because it is
    the one axis all three models can produce -- Model A has no desensitisation and so no
    separate equilibrium, and `OperationalScalarModel.dose_response` raises rather than
    pretending otherwise. A model-comparison method check therefore has to be on this
    observable or it cannot include Model A at all.
    """
    concs = (np.asarray(concs_um, dtype=float) if concs_um is not None
             else np.logspace(-1.0, 3.5, 14))
    y = np.asarray(model.peak_dose_response(concs, pam_factor=1.0), dtype=float)
    if normalisation == "fraction_of_max":
        y = y / float(y[int(np.argmax(concs))])
    if noise_sd > 0:
        y = y + np.random.default_rng(seed).normal(0.0, noise_sd, y.shape)
        if normalisation == "fraction_of_max":
            # the validator requires the top point to be 1: it IS the normalising point,
            # and adding noise to it would mean the column was normalised by something the
            # dataset does not contain.
            y[int(np.argmax(concs))] = 1.0
    sem = np.full(y.shape, noise_sd if noise_sd > 0 else 1e-3)
    return DoseResponseDataset(
        citation_label=f"generated PEAK CRC ({label})",
        preparation="none -- generated from a model at known parameters",
        concs_um=concs, mean_response=y, sem=sem,
        observable=Observable.PEAK, normalisation=normalisation,
        kind="synthetic", synthetic=True,
        origin=(f"model.peak_dose_response() at {len(concs)} concentrations"
                + (f" plus N(0, {noise_sd}) noise, seed {seed}" if noise_sd else "")
                + (", normalised to its highest concentration"
                   if normalisation == "fraction_of_max" else "")
                + ". A METHOD CHECK for roadmap P5: it exists to verify that the model "
                  "comparison recovers a model it generated the data from. No claim about "
                  "receptors may rest on it."),
        motivated_by="", role="train")


def equilibrium_crc_from_model(model, *, concs_um=None, noise_sd: float = 0.0,
                               seed: int = 0, label: str = "method check") -> DoseResponseDataset:
    """An EQUILIBRIUM concentration-response generated from a model at known parameters.

    THE ONE LEGITIMATE USE OF GENERATED DATA (roadmap P5's precedent): checking that an
    estimation method recovers an answer it is known to have. Nothing about receptors may
    rest on it, and it is labelled `synthetic` so it taints anything it touches.

    It exists because P3's equilibrium identifiability analysis has NO dataset to run on:
    both concentration-response curves in this module are PEAK, and feeding one to an
    equilibrium objective drives D to its bound by deleting desensitisation
    (`identifiability.equilibrium_chi2_identifiable` now refuses). Until a real
    EQUILIBRIUM curve is digitised -- see MISSING_DATASETS -- the method can be VALIDATED
    here and applied nowhere.
    """
    concs = (np.asarray(concs_um, dtype=float) if concs_um is not None
             else np.logspace(-1.0, 3.5, 14))
    y = np.asarray(model.dose_response(concs, pam_factor=1.0), dtype=float)
    if noise_sd > 0:
        y = y + np.random.default_rng(seed).normal(0.0, noise_sd, y.shape)
    # The declared sem must BE the noise that was added, or the two uncertainty machineries
    # disagree for a reason that has nothing to do with either of them: a chi-squared
    # profile divides by this sem while the posterior estimates the scale from the
    # residuals, so a sem floored five times too high widened the profile intervals by five
    # times (measured at noise_sd = 2e-4 against a 1e-3 floor). A NOISELESS curve still
    # needs a nonzero scale for chi-squared to be finite, so that one case keeps a declared
    # placeholder -- and it is a placeholder, not a measurement.
    sem = np.full(y.shape, noise_sd if noise_sd > 0 else 1e-3)
    return DoseResponseDataset(
        citation_label=f"generated EQUILIBRIUM CRC ({label})",
        preparation="none -- generated from a model at known parameters",
        concs_um=concs, mean_response=y, sem=sem,
        observable=Observable.EQUILIBRIUM, kind="synthetic", synthetic=True,
        origin=(f"model.dose_response() evaluated at {len(concs)} concentrations"
                + (f" plus N(0, {noise_sd}) noise, seed {seed}" if noise_sd else "")
                + ". A METHOD CHECK: it exists to verify that estimation recovers known "
                  "parameters. No claim about receptors may rest on it."),
        motivated_by="", role="train")
