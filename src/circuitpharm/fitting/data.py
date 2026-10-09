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


@dataclass(frozen=True)
class DoseResponseDataset:
    """Concentration-response measurements, or a stand-in for them. See `kind`."""

    citation_label: str          # display text only; see `source_key` for provenance
    preparation: str
    concs_um: np.ndarray
    mean_response: np.ndarray
    sem: np.ndarray
    observable: Observable = Observable.PEAK
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
# THE HILL SLOPE IS A STRUCTURAL PROBLEM FOR THE 5-STATE SCHEME, and it is the most
# interesting thing this source brings. Measured on the current fit, over the 10-90% band:
#
#     observable                 EC50       Hill slope
#     scheme, PEAK              19.946 uM     1.294
#     scheme, EQUILIBRIUM        5.222 uM     1.654
#     Jahn 1997, peak current   11.6   uM     2.2 +/- 0.4   (slope over 1-10 uM)
#
# The scheme's peak curve is markedly SHALLOWER than the measurement (1.29 against 2.2) as
# well as sitting at twice the concentration -- and the slope gap is structural, not a
# fitting failure: the scheme has two equivalent binding sites, and Jahn et al. read their
# slope of 2.2 as evidence for at least three. Fitting to this dataset therefore cannot
# succeed by moving rates; it needs a different state diagram.
#
# That is exactly what P5's model comparison exists to adjudicate, and it is the first
# falsifiable mismatch in this project between the kinetic scheme and a sourced
# measurement. It must NOT be fitted away, and the EC50 must not be quietly substituted
# for the 20 uM target either (see point 1 above).
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

# Propagated 1-sigma spread from the published parameter errors, by evaluating the Hill
# curve at EC50 +/- SEM and nH +/- SEM and taking the half-range. NOT a measured SEM per
# point -- the paper's per-point errors are in the figure, which was not digitised -- so it
# is an uncertainty on the CURVE, and it is floored at 0.01 so no point claims more
# precision than a normalised macroscopic current can carry.
def _jahn_spread() -> np.ndarray:
    lo_hi = []
    for ec50 in (_JAHN_EC50_UM - _JAHN_EC50_SEM, _JAHN_EC50_UM + _JAHN_EC50_SEM):
        for nh in (_JAHN_HILL - _JAHN_HILL_SEM, _JAHN_HILL + _JAHN_HILL_SEM):
            x = (_JAHN_CONCS / ec50) ** nh
            lo_hi.append(x / (1.0 + x))
    stack = np.vstack(lo_hi)
    return np.maximum((stack.max(axis=0) - stack.min(axis=0)) / 2.0, 0.01)


JAHN1997_PEAK_CRC = DoseResponseDataset(
    citation_label="Jahn 1997 alpha1beta2gamma2L peak GABA concentration-response",
    preparation=("alpha1beta2gamma2L transiently expressed in HEK293; patch clamp with "
                 "ultra-fast solution exchange; peak current"),
    concs_um=_JAHN_CONCS,
    mean_response=_JAHN_RESP,
    sem=_jahn_spread(),
    observable=Observable.PEAK,
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
MISSING_DATASETS = {
    "deactivation_peak_pulse": (
        "A macroscopic deactivation time course for alpha1beta2gamma2 after a brief "
        "saturating GABA pulse, with the pulse duration, concentration, temperature and n "
        "stated. This is the dataset that makes an absolute rate (k_off) identifiable at "
        "all; without it the fit determines only ratios. BLOCKED: full texts unreachable "
        "from this environment."),
    "single_channel_mean_open_time": (
        "Mean open time for alpha1beta2gamma2L. Would replace "
        "gabaa_kinetics.FIT_FIXED_ALPHA -- currently a CONVENTION -- with a measurement "
        "and make the whole kinetic fit data-determined. Jahn et al. 1997 is the right "
        "paper and the right preparation, but its abstract reports BURST duration "
        "(10.3 +/- 3.0 ms), which is a different quantity. BLOCKED: needs the full text."),
    "holdout": (
        "Any observable not used in fitting -- paired-pulse recovery at a stated interval, "
        "a desensitisation onset time course, or a PAM concentration-response at a "
        "different ambient GABA. P5's cross-validation and P6's falsification bound both "
        "consume it, so until one exists neither can be scored out-of-sample."),
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
