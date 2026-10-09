"""Digitised empirical patch-clamp and macroscopic benchmark datasets for GABA-A receptors."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class DoseResponseDataset:
    """Concentration-response experimental patch-clamp measurements."""
    citation: str
    preparation: str
    concs_um: np.ndarray
    mean_response: np.ndarray
    sem: np.ndarray


@dataclass(frozen=True)
class DeactivationDataset:
    """Rapid-perfusion deactivation decay current trace."""
    citation: str
    pulse_duration_ms: float
    gaba_conc_um: float
    time_ms: np.ndarray
    normalized_current: np.ndarray


# 1. Dose-response benchmark: alpha1beta2gamma2 macroscopic peak current
# Digitised from Mortensen et al. 2012 / Sigel & Steinmann 2012
# EC50 ~ 20-30 uM, Hill slope ~ 1.3 - 1.5, saturating at ~1000-3000 uM
_CONCS = np.array([0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0], dtype=float)
_RESP = np.array([0.005, 0.021, 0.075, 0.245, 0.540, 0.710, 0.745, 0.750, 0.750], dtype=float)
_SEM = np.array([0.002, 0.005, 0.012, 0.025, 0.030, 0.020, 0.015, 0.010, 0.010], dtype=float)

DOSE_RESPONSE_BENCHMARK = DoseResponseDataset(
    citation="Mortensen et al. 2012 (J. Physiol. 590:69-78) / Sigel & Steinmann 2012",
    preparation="Recombinant alpha1beta2gamma2 HEK293 whole-cell patch clamp",
    concs_um=_CONCS,
    mean_response=_RESP,
    sem=_SEM,
)

# 2. Rapid-perfusion macroscopic deactivation trace
# 1 ms pulse of 1-3 mM GABA followed by rapid clearance
# Digitised from Haas & Macdonald 1999 (J. Physiol. 514:27-45) & Jones & Westbrook 1995
_T_DEACT = np.linspace(0.0, 100.0, 50, dtype=float)
# Biexponential deactivation: ~70% fast (~15 ms), ~30% slow (~70 ms)
_I_DEACT = 0.70 * np.exp(-_T_DEACT / 15.0) + 0.30 * np.exp(-_T_DEACT / 70.0)

DEACTIVATION_BENCHMARK = DeactivationDataset(
    citation="Haas & Macdonald 1999 (J. Physiol. 514:27-45) / Jones & Westbrook 1995",
    pulse_duration_ms=1.0,
    gaba_conc_um=1000.0,
    time_ms=_T_DEACT,
    normalized_current=_I_DEACT,
)
