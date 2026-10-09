"""Electrophysiological protocols, decomposition, and experimental design."""
from .waveforms import (
    synaptic_transient,
    pulse_train,
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

__all__ = [
    "synaptic_transient",
    "pulse_train",
    "ambient_with_spillover",
    "extract_electrophys_metrics",
    "FactorialDecompositionResult",
    "decompose_modulation_gain",
    "compare_phasic_tonic_mechanisms",
    "ModelComparisonResult",
    "DiscriminatingProtocol",
    "evaluate_model_fit",
    "find_discriminating_protocol",
]
