"""The single source of model parameters (roadmap §2.5, P4-3).

WHY THIS EXISTS. The three receptor models in `models/` carry dataclass defaults that are
not calibrations. `KineticAllosteryModel`'s are a chimera of two different fits --
beta/alpha from the manuscript's config-pulse fit, k_on/k_off from the superseded
1000 uM / 0.30 ms fit -- giving K_d = 29.62 uM, a value in no fit, no commit and no
document, and reproducing neither anchor (EC50 6.34 uM against 20, P_o,max 0.8382 against
0.750). `ExtendedDesensitizationModel`'s are round numbers with no recorded origin. Both
entry points to the three-tiered headroom hierarchy used to default to the first set, so
the project's headline numbers were being produced at parameters nothing anchored.

So production code asks here, and `get()` returns the parameters together with a
`Quantity` that says what they rest on. The quantity is the point: there is no way to take
the numbers without also being handed their provenance.

WHAT IS AVAILABLE, honestly:

  kinetic_jw95 / "fitted"   The `gabaa_kinetics.fit_scheme` result. Three macroscopic
                            anchors (PEAK EC50, P_o,max, IPSC decay tau) with alpha held at
                            a declared convention so the system is square. UNCALIBRATED:
                            the anchors come from recalled literature ranges, not
                            digitised figures, and d/r is unconstrained by any of them.
                            Reproducible across SciPy versions to 6 significant figures.

  operational / "declared"  Model A's own stated Hill parameters. UNCALIBRATED, and note
                            that two of them (po_max 0.75, tau_deact 15 ms) ARE the
                            project's fit targets rather than anything fitted, so using
                            them as anchors is circular. The literature cannot replace
                            po_max either -- see _NOMINAL_DEFECT["operational"].

  extended_desens           NOTHING. There is no fit for Model C, and its defaults have no
                            origin. `get()` raises rather than handing back round numbers
                            with a tier on them, because P5's comparison would then be
                            ranking a fitted model against an unfitted one and calling the
                            difference evidence about mechanism.

  any model / "nominal"     The dataclass defaults, returned with a VOID quantity. For
                            reproducing a historical number or exercising algebra only.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Literal

from .provisional import void
from .results import Quantity, Tier

Which = Literal["fitted", "declared", "nominal"]

#: Models this registry knows about, and what it can supply for each.
AVAILABLE = {
    "kinetic_jw95": ("fitted", "nominal"),
    "operational": ("declared", "nominal"),
    "extended_desens": ("nominal",),
}

_NOMINAL_DEFECT = {
    "kinetic_jw95": (
        "a chimera of two fits: beta/alpha from the manuscript's config-pulse fit, "
        "kon/koff from the superseded 1000 uM / 0.30 ms fit, giving K_d = 29.62 uM, a "
        "value in no fit, no commit and no document. Reproduces neither anchor: EC50 "
        "6.34 uM against 20, P_o,max 0.8382 against 0.750."),
    # THE SINGLE HOME of the P_o,max verdict. It was briefly duplicated into a
    # PO_MAX_PROVENANCE module constant beside a NOMINAL_DEFECT = 0.750 literal; both were
    # removed. 0.750 has exactly one definition, gabaa_kinetics.FIT_TARGETS["po_max"], and
    # a second literal copy would be E12/E22 re-forming (see dynamic_range.py's note).
    "operational": (
        "declared round numbers; tau_deact 15 ms is the project's own fit TARGET rather "
        "than a fitted value; po_max is FIT_TARGETS['po_max'], a convention that the "
        "available measurements CANNOT REPLACE -- the quantity nearest to it is "
        "conditional (5 of 6 literature readings are intraburst or intracluster) and the "
        "one unconditional reading is 0.56 at an incomparable protocol: 1-2 ms "
        "application, gamma2S not gamma2L, no error bar, no n, marked UNPUBLISHED DATA by "
        "its own authors, and outside FIT_RANGES' [0.70, 0.80] in any case. So the "
        "literature audit of 2026-10-10 did NOT license changing it, and recorded why. "
        "See knowledge/14-literature.md section 5 and "
        "`scripts/literature_gamma2.py --section po`."),
    "extended_desens": (
        "round numbers with no recorded origin (kon 0.012, koff 0.35, beta 0.60, "
        "alpha 0.10). Reproduces neither anchor: EC50 5.18 uM, P_o,max 0.8571."),
}


@lru_cache(maxsize=8)
def _fitted_kinetic() -> tuple:
    """The fit_scheme result, cached: it costs a few seconds of ODE solves."""
    from .config import SYNAPTIC_PULSE
    from .gabaa_kinetics import FIT_FIXED_ALPHA, FIT_TARGETS, fit_scheme

    s = fit_scheme(verbose=False, pulse=dict(SYNAPTIC_PULSE))
    rates = dict(kon=s.kon, koff=s.koff, beta=s.beta, alpha=s.alpha, d=s.d, r=s.r)
    prov = (
        f"gabaa_kinetics.fit_scheme on config.SYNAPTIC_PULSE: three macroscopic anchors "
        f"(PEAK EC50 {FIT_TARGETS['ec50_um']} uM, P_o,max {FIT_TARGETS['po_max']}, IPSC "
        f"decay tau {FIT_TARGETS['tau_ms']} ms) with alpha held at FIT_FIXED_ALPHA = "
        f"{FIT_FIXED_ALPHA} ms^-1 (mean open time {1 / FIT_FIXED_ALPHA:.2f} ms) so the "
        f"system is square and the solution unique. K_d = {s.koff / s.kon:.4f} uM, "
        f"E = {s.beta / s.alpha:.4f}.")
    return rates, prov


def get(model: str, which: Which = "fitted") -> tuple[dict, Quantity]:
    """Parameters for `model`, with a `Quantity` describing what they rest on.

    Returns `(rates, provenance_quantity)`. The quantity's value is the parameter dict
    itself, so a caller that wants the numbers without the provenance has to go out of its
    way -- and for `which="nominal"` it is VOID, so reading it raises.
    """
    if model not in AVAILABLE:
        raise KeyError(
            f"unknown model {model!r}; this registry knows {sorted(AVAILABLE)}")
    if which not in ("fitted", "declared", "nominal"):
        raise ValueError(f"`which` must be fitted, declared or nominal; got {which!r}")

    if which == "nominal":
        rates = _nominal(model)
        return rates, void(
            f"{model}_nominal_parameters", rates,
            defect=_NOMINAL_DEFECT[model], register_item="P0-13",
            promote_by=("fit the model to sourced data on each dataset's own observable "
                        "(P4), or for Model A digitise the Hill parameters"),
            caveats=("for reproducing a historical number or exercising algebra only; not "
                     "for any quantity a reader might quote",))

    if which not in AVAILABLE[model]:
        raise NotImplementedError(
            f"there is no {which!r} parameter set for {model!r}; this registry can supply "
            f"{AVAILABLE[model]}. "
            + ("Model C has never been fitted and its defaults have no recorded origin, so "
               "returning them with a tier on them would let P5 rank a fitted model "
               "against an unfitted one and call the difference evidence about mechanism. "
               "Fit it first (P4) or compare on `nominal` for both and say so."
               if model == "extended_desens" else
               f"Use one of {AVAILABLE[model]}."))

    if model == "kinetic_jw95":
        rates, prov = _fitted_kinetic()
        return dict(rates), Quantity(
            name="kinetic_jw95_fitted_parameters", _value=dict(rates),
            tier=Tier.UNCALIBRATED, provenance=prov,
            promote_by=("digitise the three anchor observables from named figures and "
                        "replace FIT_FIXED_ALPHA with a measured mean open time -- see "
                        "fitting.data.MISSING_DATASETS -- which makes the fit "
                        "data-determined rather than convention-determined"),
            caveats=("the three anchors come from recalled literature ranges, not "
                     "digitised figures (config.CALIBRATIONS rates this UNCALIBRATED)",
                     "d and r are NOT fitted by any anchor, and every asymptotic headroom "
                     "number is conditional on d/r = 25",
                     "alpha is a CONVENTION, not a measurement; without it the fit is not "
                     "reproducible across SciPy versions (roadmap P0-7)",
                     "the one sourced value for this preparation's EC50 is 11.6 uM "
                     "(Jahn 1997), not the 20 uM anchor used here"))

    # operational / declared
    rates = _nominal("operational")
    return rates, Quantity(
        name="operational_declared_parameters", _value=rates, tier=Tier.UNCALIBRATED,
        provenance=("Model A's own declared Hill parameters. EC50 25 uM and nH 1.4 are "
                    "stated round numbers; po_max 0.75 and tau_deact 15 ms are the "
                    "project's fit TARGETS. The 2026-10-10 literature audit found po_max "
                    "unreplaceable from the literature rather than merely unmeasured; the "
                    "reasoning is in _NOMINAL_DEFECT[\"operational\"], not repeated here."),
        promote_by="digitise a concentration-response and fit the Hill parameters to it",
        caveats=("two of these are the project's own fit targets, so using them as "
                 "independent anchors is circular",
                 "the one sourced slope for this preparation, 2.2 +/- 0.4 (Jahn 1997), "
                 "is a LOCAL slope over 1-10 uM and is NOT comparable to this whole-curve "
                 "nH of 1.4; published whole-curve fits for a1b2g2 are 1.3-1.6, so 1.4 is "
                 "in range. See knowledge/12-inference.md section 2"))


def _nominal(model: str) -> dict:
    from .models.extended_desens import ExtendedDesensitizationModel
    from .models.kinetic_jw95 import KineticAllosteryModel
    from .models.operational import OperationalScalarModel

    cls = {"kinetic_jw95": KineticAllosteryModel,
           "operational": OperationalScalarModel,
           "extended_desens": ExtendedDesensitizationModel}[model]
    return cls().get_params()


def build(model: str, which: Which = "fitted"):
    """`get()` plus construction: returns `(model_instance, provenance_quantity)`."""
    from .models.extended_desens import ExtendedDesensitizationModel
    from .models.kinetic_jw95 import KineticAllosteryModel
    from .models.operational import OperationalScalarModel

    rates, prov = get(model, which)
    cls = {"kinetic_jw95": KineticAllosteryModel,
           "operational": OperationalScalarModel,
           "extended_desens": ExtendedDesensitizationModel}[model]
    return cls(**rates), prov


def report() -> str:
    """What the registry can supply, and what each set rests on."""
    out = ["=" * 78, "PARAMETER REGISTRY", "=" * 78]
    for model, kinds in AVAILABLE.items():
        out.append(f"\n{model}: {', '.join(kinds)}")
        for kind in kinds:
            if kind == "nominal":
                out.append(f"  {kind:<9} VOID -- {_NOMINAL_DEFECT[model][:70]}...")
                continue
            _, q = get(model, kind)
            out.append(f"  {kind:<9} {q.tier.label}")
            out.append(f"            {q.provenance[:100]}...")
    missing = [m for m, k in AVAILABLE.items() if "fitted" not in k]
    out.append(f"\nNO FIT EXISTS FOR: {', '.join(missing)}. P5 cannot rank these against a "
               f"fitted model.")
    return "\n".join(out)
