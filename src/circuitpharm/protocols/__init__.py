"""Electrophysiological protocols, decomposition, and experimental design."""
from .waveforms import (
    synaptic_transient,
    synaptic_transient_biexp_clearance,
    pulse_train,
    paired_pulse_ratio,
    ambient_with_spillover,
    extract_electrophys_metrics,
)
from .decomposition import (
    FactorialDecompositionResult,
    decompose_modulation_gain,
    compare_phasic_tonic_mechanisms,
)
from .oed import (
    ModelComparisonResult,
    DiscriminatingProtocol,
    evaluate_model_fit,
    find_discriminating_protocol,
)
from .design import (
    DEFAULT_SIGMA_MEAS,
    DesignResult,
    ProtocolPoint,
    design_report,
    find_discriminating_protocol_pp,
    posterior_predictive_score,
    sigma_sensitivity,
)

__all__ = [
    "synaptic_transient",
    "synaptic_transient_biexp_clearance",
    "pulse_train",
    "paired_pulse_ratio",
    "ambient_with_spillover",
    "extract_electrophys_metrics",
    "FactorialDecompositionResult",
    "decompose_modulation_gain",
    "compare_phasic_tonic_mechanisms",
    "ModelComparisonResult",
    "DiscriminatingProtocol",
    "evaluate_model_fit",
    "find_discriminating_protocol",
    # P6: the posterior-predictive design
    "ProtocolPoint",
    "DesignResult",
    "DEFAULT_SIGMA_MEAS",
    "posterior_predictive_score",
    "find_discriminating_protocol_pp",
    "sigma_sensitivity",
    "design_report",
]
