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
from .models.kinetic_jw95 import KineticAllosteryModel
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
    gaba_phasic = synaptic_transient(t, peak_um=cleft_peak_um, rise_ms=0.1, clear_ms=1.0)
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
