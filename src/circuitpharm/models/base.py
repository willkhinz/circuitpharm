"""Abstract protocol and base result structures for GABA-A receptor models."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable, Any
import numpy as np


@dataclass(frozen=True)
class WaveformResult:
    """Standardized electrophysiological output from a receptor model simulation."""
    t: np.ndarray             # time array in ms
    gaba: np.ndarray          # agonist concentration array in uM
    p_open: np.ndarray        # open probability over time
    peak_p_open: float        # maximum open probability attained
    charge_integral: float    # integral of open probability over time (ms)
    decay_tau_ms: float       # estimated weighted decay constant after agonist peak (ms)
    steady_state_tail: float  # final open probability in the trace
    states: np.ndarray | None = None  # optional full state trajectory (shape: n_times x n_states)

    def current_pA(self, g_max_ns: float = 1.0, v_hold_mv: float = -70.0, e_rev_mv: float = -70.0) -> np.ndarray:
        """Calculate macro current trace in pA given driving force."""
        # I = g * (V - E_rev); g in nS, V in mV -> I in pA
        return g_max_ns * self.p_open * (v_hold_mv - e_rev_mv)


@runtime_checkable
class ReceptorModel(Protocol):
    """Protocol satisfied by all competing GABA-A kinetic and operational models."""

    @property
    def name(self) -> str:
        """Human-readable identifier of the model."""
        ...

    @property
    def param_names(self) -> tuple[str, ...]:
        """Names of the model's tunable parameters."""
        ...

    def get_params(self) -> dict[str, float]:
        """Return a dictionary of current parameters."""
        ...

    def with_params(self, **kwargs: float) -> "ReceptorModel":
        """Return a copy of the model with updated parameter values."""
        ...

    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float:
        """Calculate equilibrium open probability at constant GABA concentration.

        Parameters
        ----------
        gaba_um : float
            Agonist concentration in micromolar.
        pam_factor : float
            Allosteric modulation multiplier (>= 1.0 for positive modulation).
        """
        ...

    def dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0) -> np.ndarray:
        """Evaluate steady-state open probability across an array of GABA concentrations."""
        ...

    def simulate_waveform(
        self,
        t: np.ndarray,
        gaba_t: np.ndarray,
        pam_factor: float = 1.0,
        initial_state: np.ndarray | None = None,
    ) -> WaveformResult:
        """Simulate time-domain open probability in response to a time-varying GABA transient.

        Parameters
        ----------
        t : np.ndarray
            1D array of timepoints in ms.
        gaba_t : np.ndarray
            1D array of GABA concentrations in uM matching `t`.
        pam_factor : float
            Allosteric modulation factor.
        initial_state : np.ndarray, optional
            Initial state vector or initial open probability.
        """
        ...
