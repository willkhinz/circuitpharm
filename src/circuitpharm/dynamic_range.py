"""Three-tiered Dynamic Range Evaluation and Headroom Hierarchy.

Distinguishes:
1. Theoretical Asymptotic Headroom (infinite affinity limit koff -> 0+).
2. Pharmacologically Reachable Gain (bounded by finite ligand operational efficacy s_max).
3. Physiological Charge Transfer Ratio (realised integral I dt under time-varying waveforms).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .models.kinetic_jw95 import KineticAllosteryModel
from .protocols.waveforms import synaptic_transient


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

    return DynamicRangeEvaluation(
        ambient_gaba_um=amb,
        theoretical_asymptotic_headroom=asymptotic_headroom,
        reachable_gains_by_smax=reachable,
        physiological_charge_ratio=charge_ratio,
        phasic_charge_gain=phasic_gain,
        tonic_charge_gain=tonic_gain,
        headroom_collapses=collapses,
    )
