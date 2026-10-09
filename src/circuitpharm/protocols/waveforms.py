"""Electrophysiological GABA exposure waveforms and metric extraction."""
from __future__ import annotations

from typing import Sequence
import numpy as np

from ..models.base import WaveformResult


def synaptic_transient(
    t: np.ndarray,
    peak_um: float = 1000.0,
    rise_ms: float = 0.10,
    clear_ms: float = 1.00,
    t_start: float = 0.0,
) -> np.ndarray:
    """Generate a biexponential synaptic GABA transient in the cleft.

    Parameters
    ----------
    t : np.ndarray
        Timepoints in ms.
    peak_um : float
        Peak cleft GABA concentration in uM (typically 1000 - 3000 uM).
    rise_ms : float
        Activation/rise time constant in ms (typically 0.05 - 0.3 ms).
    clear_ms : float
        Clearance/unbinding time constant in ms (typically 0.3 - 2.0 ms).
    t_start : float
        Onset time of the release event in ms.
    """
    t_rel = np.maximum(t - t_start, 0.0)
    mask = t >= t_start

    if rise_ms <= 1e-4:
        # Monoexponential instantaneous jump and decay
        gaba = np.where(mask, peak_um * np.exp(-t_rel / max(clear_ms, 1e-4)), 0.0)
    else:
        # Biexponential profile: exp(-t/clear) - exp(-t/rise)
        r = rise_ms
        c = clear_ms
        t_peak = (r * c / (c - r)) * np.log(c / r) if c != r else r
        norm_factor = np.exp(-t_peak / c) - np.exp(-t_peak / r)
        norm_factor = max(norm_factor, 1e-6)
        
        raw = np.exp(-t_rel / c) - np.exp(-t_rel / r)
        gaba = np.where(mask, peak_um * (raw / norm_factor), 0.0)

    return np.asarray(gaba, dtype=float)


def pulse_train(
    t: np.ndarray,
    freq_hz: float = 20.0,
    n_pulses: int = 5,
    peak_um: float = 1000.0,
    rise_ms: float = 0.10,
    clear_ms: float = 1.00,
    start_ms: float = 10.0,
) -> np.ndarray:
    """Generate a train of synaptic GABA pulses to probe desensitisation accumulation.

    Parameters
    ----------
    t : np.ndarray
        Time vector in ms.
    freq_hz : float
        Repetition frequency in Hz (e.g. 10 Hz, 50 Hz, 100 Hz).
    n_pulses : int
        Number of consecutive synaptic events.
    peak_um : float
        Peak GABA concentration per pulse in uM.
    rise_ms : float
        Rise time constant (ms).
    clear_ms : float
        Clearance time constant (ms).
    start_ms : float
        Onset of first pulse in ms.
    """
    period_ms = 1000.0 / max(freq_hz, 0.1)
    gaba = np.zeros_like(t, dtype=float)
    for i in range(n_pulses):
        t_pulse = start_ms + i * period_ms
        gaba += synaptic_transient(t, peak_um=peak_um, rise_ms=rise_ms, clear_ms=clear_ms, t_start=t_pulse)
    return gaba


def ambient_with_spillover(
    t: np.ndarray,
    ambient_um: float = 0.40,
    spillover_peak_um: float = 2.0,
    event_times_ms: Sequence[float] | None = None,
    tau_spillover_ms: float = 30.0,
) -> np.ndarray:
    """Generate sustained ambient GABA with activity-dependent spillover transients.

    Parameters
    ----------
    t : np.ndarray
        Timepoints in ms.
    ambient_um : float
        Basal tonic extrasynaptic GABA in uM.
    spillover_peak_um : float
        Peak magnitude of spillover transient above baseline.
    event_times_ms : Sequence[float], optional
        Onset times of synaptic spillover events.
    tau_spillover_ms : float
        Decay constant of extrasynaptic clearance.
    """
    gaba = np.full_like(t, float(ambient_um), dtype=float)
    if event_times_ms is not None:
        for t_ev in event_times_ms:
            t_rel = np.maximum(t - t_ev, 0.0)
            mask = t >= t_ev
            gaba += np.where(mask, spillover_peak_um * np.exp(-t_rel / max(tau_spillover_ms, 1e-3)), 0.0)
    return gaba


def extract_electrophys_metrics(
    res: WaveformResult,
    g_max_ns: float = 1.0,
    driving_force_mv: float = 30.0,
) -> dict[str, float]:
    """Compute standard electrophysiological summary metrics from a WaveformResult.

    driving_force_mv: |V_hold - E_Cl| in mV.
    """
    current_pA = g_max_ns * res.p_open * driving_force_mv
    peak_current = float(np.max(current_pA))

    _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
    # Charge Q = integral I dt (pA * ms = fC)
    charge_fC = float(_trapz(current_pA, res.t))

    return {
        "peak_open_prob": res.peak_p_open,
        "peak_current_pA": peak_current,
        "charge_integral_fC": charge_fC,
        "decay_tau_ms": res.decay_tau_ms,
        "steady_state_tail": res.steady_state_tail,
    }
