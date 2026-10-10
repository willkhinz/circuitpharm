"""Three-tiered Dynamic Range Evaluation and Headroom Hierarchy.

Distinguishes:
1. Theoretical Asymptotic Headroom (infinite affinity limit koff -> 0+).
2. Pharmacologically Reachable Gain (bounded by finite ligand operational efficacy s_max).
3. Physiological Charge Transfer Ratio (realised integral I dt under time-varying waveforms).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .models.base import Observable, ObservedQuantity
from .provisional import void
from .models.kinetic_jw95 import KineticAllosteryModel
from .config import SYNAPTIC_PULSE
from .protocols.waveforms import synaptic_transient
from .results import Quantity, ResultSet, Tier


@dataclass(frozen=True)
class DynamicRangeEvaluation:
    """Structured three-tiered headroom evaluation."""
    ambient_gaba_um: float
    theoretical_asymptotic_headroom: float
    reachable_gains_by_smax: dict[float, float]
    physiological_charge_ratio: float
    phasic_charge_gain: float
    tonic_charge_gain: float
    headroom_collapses: bool  # True if tonic advantage drops <= 1.0
    #: The window the CHARGE tier was integrated over. A charge without its window is not
    #: a quantity (models.base.ObservedQuantity), and the ratio depends on it.
    charge_window_ms: float = float("nan")
    #: Where the tonic advantage reaches 1.0, as a curve rather than the boolean above.
    #: Maps ambient GABA (uM) -> the charge ratio there. A boolean throws away everything
    #: the scan computed (roadmap P6-2).
    collapse_scan: tuple = ()


# THE CLEFT PEAK DEFAULT DISAGREES WITH config.SYNAPTIC_PULSE, DELIBERATELY LEFT FOR A
# DECISION. `cleft_peak_um` defaults to 1000.0 in both public entry points below, while
# config.SYNAPTIC_PULSE["peak_um"] is 3000.0. So this module evaluates the dynamic range on a
# transient that matches neither of the project's two historical calibrations:
#
#     config.SYNAPTIC_PULSE          3000 uM / 1.00 ms   <- what the manuscript quotes
#     superseded gabaa_kinetics      1000 uM / 0.30 ms   <- removed by E22
#     this module's defaults         1000 uM / 1.00 ms   <- neither
#
# Nothing published is affected: the manuscript's figures come from gabaa_kinetics via
# scripts/paper_numbers.py, and no script calls this module -- its only callers are tests,
# one of which (test_nextgen_models.py) pins cleft_peak_um=1000.0 explicitly. But this is the
# newer tiered/posterior implementation of the SAME quantity, so the moment the manuscript
# sources its dynamic range from here the two will disagree on the input.
#
# That is recurring error E12/E22 re-forming in the replacement layer. The clearance literal
# has been pointed at config (behaviour-identical, 1.0 == 1.0, which is why it was invisible);
# the PEAK is a real numerical choice and changing it would move every dynamic-range figure
# this module produces, so it is flagged rather than silently switched.
def evaluate_dynamic_range(
    model: KineticAllosteryModel,
    ambient_gaba_um: float = 0.40,
    cleft_peak_um: float = 1000.0,
    smax_values: tuple[float, ...] = (1.25, 1.50, 2.40, 2.50),
    transient_duration_ms: float = 50.0,
) -> DynamicRangeEvaluation:
    """Evaluate the three distinct tiers of dynamic range and headroom.

    Parameters
    ----------
    model : KineticAllosteryModel
        REQUIRED. There is deliberately no default: this function used to fall back to
        `KineticAllosteryModel()`, whose parameters are a chimera of two fits reproducing
        neither anchor (EC50 6.34 uM against 20, P_o,max 0.8382 against 0.750), so the
        three-tiered headroom hierarchy -- the whole point of Advancement 5 -- was being
        computed at parameters nothing anchored. Obtain one from
        `circuitpharm.parameters.get(...)` (roadmap P0-13).
    ambient_gaba_um : float
        Basal ambient extrasynaptic GABA (uM).
    cleft_peak_um : float
        Peak synaptic cleft GABA (uM).
    smax_values : tuple[float, ...]
        Ligand operational potency shift factors to test.
    transient_duration_ms : float
        Window length for dynamic charge integration (ms).
    """
    if model is None:
        raise TypeError(
            "evaluate_dynamic_range() requires an explicit `model`. The former fallback to "
            "KineticAllosteryModel() produced headline headroom numbers at uncalibrated "
            "chimeric defaults; see roadmap P0-13.")
    m = model
    amb = max(float(ambient_gaba_um), 1e-4)

    # 1. Tier 1: Theoretical asymptotic headroom
    po_inf = m.po_inf
    po_base = m.steady_state(amb, pam_factor=1.0)
    asymptotic_headroom = float(po_inf / max(po_base, 1e-12))

    # 2. Tier 2: Reachable gains at specific s_max limits
    reachable: dict[float, float] = {}
    for sm in smax_values:
        po_pam = m.steady_state(amb, pam_factor=sm)
        reachable[sm] = float(po_pam / max(po_base, 1e-12))

    # 3. Tier 3: Physiological charge transfer ratio (using default s_max = 2.50)
    t = np.linspace(0.0, transient_duration_ms, 200)
    # Phasic waveform
    gaba_phasic = synaptic_transient(t, peak_um=cleft_peak_um, rise_ms=0.1,
                                     clear_ms=SYNAPTIC_PULSE["clear_ms"])
    res_phasic_ctrl = m.simulate_waveform(t, gaba_phasic, pam_factor=1.0)
    res_phasic_pam = m.simulate_waveform(t, gaba_phasic, pam_factor=2.50)
    
    phasic_gain = float(res_phasic_pam.charge_integral / max(res_phasic_ctrl.charge_integral, 1e-12))

    # Tonic waveform (sustained ambient)
    gaba_tonic = np.full_like(t, amb)
    res_tonic_ctrl = m.simulate_waveform(t, gaba_tonic, pam_factor=1.0)
    res_tonic_pam = m.simulate_waveform(t, gaba_tonic, pam_factor=2.50)
    
    tonic_gain = float(res_tonic_pam.charge_integral / max(res_tonic_ctrl.charge_integral, 1e-12))

    charge_ratio = float(tonic_gain / max(phasic_gain, 1e-12))
    collapses = charge_ratio <= 1.0

    # THE COLLAPSE BOUNDARY AS A CURVE. `headroom_collapses` answers "does the tonic
    # advantage survive at this one ambient concentration?"; the scan answers "where does
    # it stop surviving?", which is the question a reader actually has and which the
    # boolean discarded (roadmap P6-2). Equilibrium gains are used here, not the
    # charge integrals: this is a boundary location over a sweep, and a 200-point ODE
    # solve per grid point buys no accuracy in the answer.
    # `scan_amb`, NOT `amb`: the loop variable must not shadow the ambient concentration
    # this evaluation is ABOUT. It did, and the returned dataclass reported
    # ambient_gaba_um = 30.0 (the last grid point) for a run at 0.4 -- the same class of
    # defect as resp_metrics rebinding its `band` parameter to the FFT mask (P0-12), found
    # the same way, by reading the output rather than the code.
    scan = []
    for scan_amb in (0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 3.0, 10.0, 30.0):
        base_t = m.steady_state(scan_amb, pam_factor=1.0)
        pam_t = m.steady_state(scan_amb, pam_factor=2.50)
        if base_t <= 0.0:
            continue
        g_t = pam_t / base_t
        scan.append((float(scan_amb), float(g_t / max(phasic_gain, 1e-12))))

    return DynamicRangeEvaluation(
        ambient_gaba_um=amb,
        theoretical_asymptotic_headroom=asymptotic_headroom,
        reachable_gains_by_smax=reachable,
        physiological_charge_ratio=charge_ratio,
        phasic_charge_gain=phasic_gain,
        tonic_charge_gain=tonic_gain,
        headroom_collapses=collapses,
        charge_window_ms=float(transient_duration_ms),
        collapse_scan=tuple(scan),
    )


def dynamic_range_report(model: KineticAllosteryModel, **kw) -> ResultSet:
    """The three tiers as tiered quantities (roadmap P2).

    THE TIERS OF THE METRIC AND THE TIERS OF THE EVIDENCE ARE DIFFERENT THINGS, and
    conflating them is how "184.6x" ended up quoted as a result. `DynamicRangeEvaluation`'s
    three tiers say WHAT is being measured -- asymptotic, reachable, realised. The
    `results.Tier` on each quantity says how much it can be trusted. Every one here is
    UNCALIBRATED: the mechanism is sound, the magnitudes rest on three anchors from recalled
    literature ranges, one declared convention (`gabaa_kinetics.FIT_FIXED_ALPHA`), and d/r,
    which no anchor constrains.
    """
    ev = evaluate_dynamic_range(model, **kw)
    rs = ResultSet(f"dynamic range at ambient {ev.ambient_gaba_um} uM")
    promote = ("fit the rates to sourced data on their own observables and propagate the "
               "posterior (P4); quote a median with a credible interval, never a point")

    rs.add(Quantity(
        "theoretical_asymptotic_headroom", ev.theoretical_asymptotic_headroom,
        Tier.UNCALIBRATED, "x",
        provenance=("P_open,inf / P_open(ambient), where P_open,inf = 1/(1 + (alpha/beta)"
                    "(1 + d/r)) is the k_off -> 0+ limit. A formula, not a measurement: "
                    "its inputs are the fitted E and the UNFITTED d/r."),
        promote_by=promote,
        caveats=("asymptotic: k_off -> 0+ is reached by no ligand, so this is a ceiling "
                 "the mechanism cannot deliver, not a predicted effect",
                 "strictly conditional on d/r = 25, which no anchor constrains")))

    for sm, gain in sorted(ev.reachable_gains_by_smax.items()):
        rs.add(Quantity(
            f"reachable_gain_at_smax_{sm:g}", gain, Tier.UNCALIBRATED, "x",
            provenance=(f"equilibrium open-probability gain at an EQUILIBRIUM EC50 shift "
                        f"of {sm:g}x, which for this scheme is exactly the factor dividing "
                        f"K_d"),
            promote_by=promote,
            caveats=("s_max here is an EQUILIBRIUM EC50 shift; gabaa_kinetics.calibrate_pam "
                     "calibrates against the PEAK shift and returns 2.945 for a nominal "
                     "2.5, so the two must not be interchanged",)))

    rs.add(Quantity(
        "physiological_charge_ratio", ev.physiological_charge_ratio, Tier.UNCALIBRATED, "x",
        provenance=(f"(tonic charge gain) / (phasic charge gain), both integrated over "
                    f"{ev.charge_window_ms:g} ms at PAM 2.50x. Tonic gain "
                    f"{ev.tonic_charge_gain:.4f}x, phasic {ev.phasic_charge_gain:.4f}x."),
        promote_by=promote,
        caveats=(f"CHARGE depends on its window; this is {ev.charge_window_ms:g} ms and a "
                 f"different window gives a different number",
                 "the phasic arm is a synaptic transient and the tonic arm a sustained "
                 "step, so the two charges are not the same kind of integral -- only their "
                 "RATIO OF GAINS is dimensionless and comparable")))

    if ev.collapse_scan:
        below = [a for a, x in ev.collapse_scan if x <= 1.0]
        boundary = min(below) if below else float("nan")
        rs.add(Quantity(
            "collapse_ambient_um", boundary, Tier.UNCALIBRATED, "uM",
            provenance=("lowest ambient GABA in the scan at which the tonic advantage "
                        "reaches 1.0, i.e. where the compartment divergence stops. Scan: "
                        + ", ".join(f"{a:g}->{x:.3f}x" for a, x in ev.collapse_scan)),
            promote_by=promote,
            caveats=("a grid bound, not an interpolated crossing",
                     "AMBIENT_GABA_UM itself is UNSOURCED over a 0.1-1 uM literature "
                     "range, and this boundary is what makes that range matter")))
    return rs


def charge_quantity(ev: DynamicRangeEvaluation) -> ObservedQuantity:
    """The realised charge ratio as an `ObservedQuantity`, window attached."""
    return ObservedQuantity(ev.physiological_charge_ratio, Observable.CHARGE,
                            window_ms=ev.charge_window_ms)


# ===================================================================== P6-2: intervals
#
# WHY A SECOND ENTRY POINT RATHER THAN A FLAG ON THE FIRST. `evaluate_dynamic_range` and
# `dynamic_range_report` evaluate the three tiers AT ONE PARAMETER SET. That is the right
# thing when the question is "what does this model say", and it is what the manuscript's
# existing numbers came from, so it stays as it is and keeps reproducing them. It is the
# wrong thing when the question is "what does the DATA say", because a single parameter set
# has no uncertainty and a fold-change quoted without one is the defect the whole roadmap
# is about. So the interval version takes posterior DRAWS and has no point-estimate mode:
# you cannot call it and accidentally get a number without an interval.
@dataclass(frozen=True)
class TierPosterior:
    """One tier evaluated over posterior draws."""

    name: str
    median: float
    ci_95: tuple[float, float]
    n_draws: int
    #: Draws that could not be evaluated (a non-finite tier, a model that refused to
    #: build). Reported, because a median over the survivors of a 40%-failure rate is a
    #: median over a subset nobody chose.
    n_failed: int

    @property
    def spread_factor(self) -> float:
        """How many fold wide the interval is. 1.0 would be a point estimate."""
        lo, hi = self.ci_95
        return float(hi / lo) if lo > 0 else float("inf")


def _percentiles(values: np.ndarray) -> tuple[float, tuple[float, float]]:
    return (float(np.median(values)),
            (float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))))


def dynamic_range_posterior(
    draws,
    *,
    conventions: dict[str, float] | None = None,
    ambient_gaba_um: float = 0.40,
    smax_values: tuple[float, ...] = (1.25, 1.50, 2.50),
    cleft_peak_um: float = 1000.0,
    transient_duration_ms: float = 50.0,
    charge_smax: float = 2.50,
    n_draws: int | None = 200,
    seed: int = 0,
    include_charge: bool = True,
) -> dict[str, TierPosterior]:
    """The three tiers over posterior draws of `(log10 K_d, log10 E, log10 D)`.

    `draws` is `(n, 3)` in that order -- `PosteriorResult.log10_samples`, or any array of
    them. `conventions` supplies `koff`, `alpha` and `r`, which no dataset in this
    repository constrains (`fitting.data.MISSING_DATASETS`); it is required in the sense
    that the defaults are `comparison.FIT_CONVENTIONS` and every absolute rate is
    conditional on them.

    `n_draws` THINS THE CHAIN. Posterior draws are autocorrelated, so 200 taken evenly
    across the chain carry nearly as much information as 20,000 taken consecutively and
    cost a hundredth as much -- which matters because tier 3 is four ODE integrations per
    draw. Thinning is by stride, not by taking a prefix: a prefix is the start of the chain
    and may be burn-in that was not discarded.

    `include_charge=False` skips tier 3, for the two analytic tiers at a thousand draws.
    """
    from .fitting.comparison import FIT_CONVENTIONS

    conv = dict(conventions or {k: v for k, (v, _) in FIT_CONVENTIONS.items()
                                if k in ("koff", "alpha", "r")})
    for key in ("koff", "alpha", "r"):
        if key not in conv:
            raise ValueError(
                f"`conventions` must supply {key!r}: no dataset here carries a timescale, "
                f"so the absolute rates are not identifiable and this function will not "
                f"invent one. See fitting.comparison.FIT_CONVENTIONS.")

    d = np.atleast_2d(np.asarray(draws, dtype=float))
    if d.shape[1] != 3:
        raise ValueError(
            f"draws must be (n, 3) for (log10_kd, log10_E, log10_D); got {d.shape}")
    if n_draws is not None and d.shape[0] > n_draws:
        d = d[:: max(1, d.shape[0] // n_draws)][:n_draws]
    if d.shape[0] < 2:
        raise ValueError(
            f"{d.shape[0]} draw(s) is not a posterior. An interval needs a distribution; "
            f"for a point evaluation call `evaluate_dynamic_range`, which says it is one.")

    amb = max(float(ambient_gaba_um), 1e-4)
    collected: dict[str, list[float]] = {"theoretical_asymptotic_headroom": []}
    for sm in smax_values:
        collected[f"reachable_gain_at_smax_{sm:g}"] = []
    if include_charge:
        collected["physiological_charge_ratio"] = []
        collected["phasic_charge_gain"] = []
        collected["tonic_charge_gain"] = []
    failed = 0

    t = np.linspace(0.0, transient_duration_ms, 200)
    gaba_phasic = synaptic_transient(t, peak_um=cleft_peak_um, rise_ms=0.1,
                                     clear_ms=SYNAPTIC_PULSE["clear_ms"])
    gaba_tonic = np.full_like(t, amb)

    for row in d:
        kd, e_ratio, d_ratio = (10.0 ** float(x) for x in row)
        try:
            m = KineticAllosteryModel(
                kon=conv["koff"] / kd, koff=conv["koff"],
                beta=conv["alpha"] * e_ratio, alpha=conv["alpha"],
                d=conv["r"] * d_ratio, r=conv["r"])
            po_base = m.steady_state(amb, pam_factor=1.0)
            if not (np.isfinite(po_base) and po_base > 0.0):
                raise FloatingPointError("zero baseline open probability")
            row_vals = {"theoretical_asymptotic_headroom": float(m.po_inf / po_base)}
            for sm in smax_values:
                row_vals[f"reachable_gain_at_smax_{sm:g}"] = float(
                    m.steady_state(amb, pam_factor=sm) / po_base)
            if include_charge:
                q_pc = m.simulate_waveform(t, gaba_phasic, pam_factor=1.0).charge_integral
                q_pp = m.simulate_waveform(t, gaba_phasic,
                                           pam_factor=charge_smax).charge_integral
                q_tc = m.simulate_waveform(t, gaba_tonic, pam_factor=1.0).charge_integral
                q_tp = m.simulate_waveform(t, gaba_tonic,
                                           pam_factor=charge_smax).charge_integral
                if min(q_pc, q_tc) <= 0:
                    raise FloatingPointError("zero control charge")
                phasic = float(q_pp / q_pc)
                tonic = float(q_tp / q_tc)
                row_vals["phasic_charge_gain"] = phasic
                row_vals["tonic_charge_gain"] = tonic
                row_vals["physiological_charge_ratio"] = float(tonic / phasic)
            if not all(np.isfinite(v) for v in row_vals.values()):
                raise FloatingPointError("non-finite tier")
        except Exception:
            failed += 1
            continue
        for k, v in row_vals.items():
            collected[k].append(v)

    kept = len(collected["theoretical_asymptotic_headroom"])
    if kept < 2:
        raise RuntimeError(
            f"only {kept} of {d.shape[0]} draws evaluated. A median over that is not a "
            f"posterior summary; check the conventions and the ambient concentration.")

    out = {}
    for name, vals in collected.items():
        med, ci = _percentiles(np.asarray(vals, dtype=float))
        out[name] = TierPosterior(name=name, median=med, ci_95=ci, n_draws=kept,
                                  n_failed=failed)
    return out


def dynamic_range_posterior_report(
    draws, *, posterior_tier: Tier = Tier.UNCALIBRATED, posterior_note: str = "",
    **kw) -> ResultSet:
    """The interval version as tiered quantities: median, 95% interval, window, observable.

    `posterior_tier` IS THE CALLER'S TO SUPPLY AND CANNOT BE BETTER THAN THE POSTERIOR'S.
    An interval from a chain that failed its diagnostics is VOID however carefully the
    tiers are then computed, and this function has no way to check the chain it was handed
    an array from -- so the caller passes the tier its `PosteriorResult` carried, and
    passing a better one is the kind of thing a reviewer can see.
    """
    tiers = dynamic_range_posterior(draws, **kw)
    amb = kw.get("ambient_gaba_um", 0.40)
    window = kw.get("transient_duration_ms", 50.0)
    rs = ResultSet(f"dynamic range at ambient {amb} uM, over "
                   f"{tiers['theoretical_asymptotic_headroom'].n_draws} posterior draws")

    def quantity(name: str, units: str, provenance: str, caveats: tuple) -> Quantity:
        tp = tiers[name]
        full_prov = (f"{provenance} Posterior median over {tp.n_draws} draws with a "
                     f"[2.5, 97.5] percentile interval of "
                     f"[{tp.ci_95[0]:.4g}, {tp.ci_95[1]:.4g}] {units} -- a "
                     f"{tp.spread_factor:.3g}-fold spread."
                     + (f" {tp.n_failed} draw(s) could not be evaluated and are excluded."
                        if tp.n_failed else "")
                     + (f" {posterior_note}" if posterior_note else ""))
        if posterior_tier is Tier.VOID:
            return void(name, (tp.median, tp.ci_95),
                        defect=("summarised from a posterior the caller declared VOID: "
                                + (posterior_note or "its chain failed a diagnostic")),
                        register_item="P6-2",
                        promote_by="run the chain to convergence (see fitting/posterior.py)",
                        caveats=caveats)
        return Quantity(name=name, _value=tp.median, tier=posterior_tier, units=units,
                        provenance=full_prov,
                        promote_by=("digitise the anchor observables and acquire a kinetic "
                                    "dataset so d/r and the absolute rates are fitted "
                                    "rather than conventional"),
                        caveats=caveats + (
                            f"the 95% interval is [{tp.ci_95[0]:.4g}, {tp.ci_95[1]:.4g}], "
                            f"a {tp.spread_factor:.3g}-fold spread; the median alone is not "
                            f"the result",))

    rs.add(quantity(
        "theoretical_asymptotic_headroom", "x",
        ("P_open,inf / P_open(ambient), with P_open,inf = E/(1 + E + D) -- the k_off -> 0+ "
         "limit, which is a FORMULA in the two fitted ratios and not a measurement."),
        ("asymptotic: k_off -> 0+ is reached by no ligand, so this is a ceiling the "
         "mechanism cannot deliver, not a predicted effect",
         "D is sampled here rather than fixed at 25, which is why this interval is wide")))

    for sm in kw.get("smax_values", (1.25, 1.50, 2.50)):
        rs.add(quantity(
            f"reachable_gain_at_smax_{sm:g}", "x",
            (f"equilibrium open-probability gain at an EQUILIBRIUM EC50 shift of {sm:g}x, "
             f"which for this scheme is exactly the factor dividing K_d."),
            ("s_max here is an EQUILIBRIUM EC50 shift; gabaa_kinetics.calibrate_pam "
             "calibrates against the PEAK shift and returns 2.945 for a nominal 2.5, so "
             "the two must not be interchanged",)))

    if "physiological_charge_ratio" in tiers:
        rs.add(quantity(
            "physiological_charge_ratio", "x",
            (f"(tonic charge gain)/(phasic charge gain), both integrated over {window:g} ms "
             f"at PAM {kw.get('charge_smax', 2.50):g}x. Observable: CHARGE."),
            (f"CHARGE depends on its window; this is {window:g} ms and a different window "
             f"gives a different number",
             "the phasic arm is a synaptic transient and the tonic arm a sustained step, so "
             "only their RATIO OF GAINS is dimensionless and comparable")))
    return rs
