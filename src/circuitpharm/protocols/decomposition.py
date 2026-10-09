"""Mechanistic decomposition of allosteric modulation into physical drivers.

Disentangles:
1. Occupancy: changes in ligand-bound fraction.
2. Kinetic Redistribution: conformational shifts between bound-shut, open, and desensitised states.
3. Operating Point Headroom: proximity to the maximal gating ceiling Po_max.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from ..models.kinetic_jw95 import KineticAllosteryModel


@dataclass(frozen=True)
class FactorialDecompositionResult:
    """Quantitative attribution of positive allosteric modulation to its 3 mechanistic drivers."""
    condition: str                     # "phasic" or "tonic"
    gaba_conc_um: float                # exposure concentration or peak (uM)
    pam_factor: float                  # modulation multiplier applied
    total_fold_gain: float             # G = P_open,pam / P_open,base
    log_gain: float                    # ln(G)
    delta_log_occupancy: float         # ln(Occ_pam / Occ_base)
    delta_log_gating: float            # ln(P(O|bound)_pam / P(O|bound)_base)
    delta_log_headroom: float          # ln(Headroom ratio)
    shapley_pct_occupancy: float       # % of log-gain attributed to occupancy
    shapley_pct_gating: float          # % of log-gain attributed to kinetic redistribution
    shapley_pct_headroom: float        # % of log-gain attributed to headroom saturation proximity


def decompose_modulation_gain(
    model: KineticAllosteryModel,
    gaba_um: float,
    pam_factor: float = 2.50,
    condition_label: str = "custom",
) -> FactorialDecompositionResult:
    """Decompose modulation gain into Occupancy, Gating, and Headroom factors.

    Parameters
    ----------
    model : KineticAllosteryModel
        Baseline 5-state kinetic model.
    gaba_um : float
        GABA concentration in uM (e.g. 0.40 for tonic, 1000.0 for synaptic peak).
    pam_factor : float
        Allosteric affinity shift factor.
    condition_label : str
        Label for the exposure regime ("tonic", "phasic", etc.).
    """
    g = max(float(gaba_um), 1e-6)
    pf = max(float(pam_factor), 1.0)

    # 1. Baseline states
    dist_base = model.state_distribution(g, pam_factor=1.0)
    po_base = float(dist_base[3])  # P(A2O)
    occ_base = float(dist_base[1] + dist_base[2] + dist_base[3] + dist_base[4])  # AR + A2R + A2O + A2D
    # Conditional gating efficiency given doubly-bound: P(A2O) / (P(A2R) + P(A2O) + P(A2D))
    a2_total_base = float(dist_base[2] + dist_base[3] + dist_base[4])
    cond_gate_base = (po_base / max(a2_total_base, 1e-12)) if a2_total_base > 0 else 0.0

    # 2. Modulated states
    dist_pam = model.state_distribution(g, pam_factor=pf)
    po_pam = float(dist_pam[3])
    occ_pam = float(dist_pam[1] + dist_pam[2] + dist_pam[3] + dist_pam[4])
    a2_total_pam = float(dist_pam[2] + dist_pam[3] + dist_pam[4])
    cond_gate_pam = (po_pam / max(a2_total_pam, 1e-12)) if a2_total_pam > 0 else 0.0

    # 3. Headroom proximity: 1 - P / Po_max
    po_max = model.po_max
    headroom_base = max(po_max - po_base, 1e-6)
    headroom_pam = max(po_max - po_pam, 1e-6)
    headroom_factor = headroom_base / headroom_pam

    # Total fold gain
    gain = max(po_pam / max(po_base, 1e-12), 1e-12)
    log_gain = float(np.log(gain))

    # Factorial log differences
    delta_occ = float(np.log(max(occ_pam / max(occ_base, 1e-12), 1e-12)))
    delta_gate = float(np.log(max(cond_gate_pam / max(cond_gate_base, 1e-12), 1e-12)))
    delta_head = float(np.log(headroom_factor))

    # Shapley attribution across the 3 terms
    total_abs = abs(delta_occ) + abs(delta_gate) + abs(delta_head)
    if total_abs > 1e-9:
        pct_occ = 100.0 * abs(delta_occ) / total_abs
        pct_gate = 100.0 * abs(delta_gate) / total_abs
        pct_head = 100.0 * abs(delta_head) / total_abs
    else:
        pct_occ = 33.33
        pct_gate = 33.33
        pct_head = 33.34

    return FactorialDecompositionResult(
        condition=condition_label,
        gaba_conc_um=g,
        pam_factor=pf,
        total_fold_gain=gain,
        log_gain=log_gain,
        delta_log_occupancy=delta_occ,
        delta_log_gating=delta_gate,
        delta_log_headroom=delta_head,
        shapley_pct_occupancy=pct_occ,
        shapley_pct_gating=pct_gate,
        shapley_pct_headroom=pct_head,
    )


def compare_phasic_tonic_mechanisms(
    model: KineticAllosteryModel | None = None,
    tonic_gaba_um: float = 0.40,
    phasic_gaba_um: float = 1000.0,
    pam_factor: float = 2.50,
) -> dict[str, FactorialDecompositionResult]:
    """Run comparative mechanistic decomposition between Phasic and Tonic conditions."""
    m = model or KineticAllosteryModel()
    tonic_res = decompose_modulation_gain(m, tonic_gaba_um, pam_factor=pam_factor, condition_label="tonic")
    phasic_res = decompose_modulation_gain(m, phasic_gaba_um, pam_factor=pam_factor, condition_label="phasic")
    return {"tonic": tonic_res, "phasic": phasic_res}
