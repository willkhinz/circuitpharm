"""Electrophysiological GABA exposure waveforms and metric extraction."""
from __future__ import annotations

from typing import Sequence
import numpy as np

from ..models.base import (DEFAULT_E_CL_MV, DEFAULT_V_HOLD_MV, WaveformResult)


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
    r, c = float(rise_ms), float(clear_ms)

    if r <= 1e-4:
        # Monoexponential instantaneous jump and decay
        gaba = np.where(mask, peak_um * np.exp(-t_rel / max(c, 1e-4)), 0.0)
    elif abs(c - r) < 1e-9 * max(c, r):
        # THE tau_rise -> tau_clear LIMIT, which previously returned an all-zero trace.
        #
        # exp(-t/c) - exp(-t/r) is identically zero when c == r, and the old code then
        # clamped its normaliser to 1e-6 and emitted `peak_um * 0/1e-6` = 0 at every
        # sample -- measured: 0 uM for a requested 1000 uM. A receptor model handed an
        # empty agonist trace returns a flat resting trace, a near-zero charge and (with
        # the old estimator) a decay tau of exactly the fit target: three plausible
        # numbers from an experiment that never happened.
        #
        # The analytic limit is the alpha function, (t/r)*exp(-t/r), whose peak is at
        # t = r with value exp(-1). Note the NEAR-equal case was always fine (c = 0.999,
        # r = 1.0 gives the full 1000 uM), so the guard is a relative tolerance rather
        # than an equality test -- do not "simplify" it back to `c == r` (roadmap P0-2).
        raw = (t_rel / r) * np.exp(-t_rel / r)
        gaba = np.where(mask, peak_um * raw / np.exp(-1.0), 0.0)
    else:
        # Biexponential profile: exp(-t/clear) - exp(-t/rise)
        t_peak = (r * c / (c - r)) * np.log(c / r)
        norm_factor = np.exp(-t_peak / c) - np.exp(-t_peak / r)
        if norm_factor <= 1e-12:
            raise ValueError(
                f"the biexponential normaliser came out as {norm_factor!r} for "
                f"rise_ms={r!r}, clear_ms={c!r}. Clamping it (the previous behaviour) "
                f"scales the whole waveform by an arbitrary factor; these time constants "
                f"are degenerate and the alpha-function limit above should have caught "
                f"them.")
        raw = np.exp(-t_rel / c) - np.exp(-t_rel / r)
        gaba = np.where(mask, peak_um * (raw / norm_factor), 0.0)

    return np.asarray(gaba, dtype=float)


def synaptic_transient_biexp_clearance(
    t: np.ndarray,
    peak_um: float = 1000.0,
    rise_ms: float = 0.10,
    clear_fast_ms: float = 1.00,
    clear_slow_ms: float = 20.0,
    weight_fast: float = 0.80,
    t_start: float = 0.0,
) -> np.ndarray:
    """Synaptic transient with TWO clearance components, normalised to `peak_um`.

    The technical spec (§5) and roadmap P6-1 both call for biexponential clearance
    (tau_fast ~1 ms, tau_slow ~10-30 ms); `synaptic_transient` has a single `clear_ms` and
    is kept unchanged as the documented special case, so every existing result stays
    reproducible (roadmap §2.6).

    Implemented as a weighted sum of two single-clearance transients sharing one rise, then
    rescaled so the realised peak is exactly `peak_um` -- the weights shape the decay, they
    do not silently change the amplitude.

    NORMALISATION IS TO THE SAMPLED PEAK, not an analytic one: the mixture's peak has no
    closed form. So the amplitude is exact to the resolution of `t`, and this function and
    `synaptic_transient` agree in SHAPE rather than bit-for-bit when the mixture collapses
    to one component (they differ by the ratio of sampled to analytic peak, ~0.02% on a
    0.025 ms grid). Compare shapes, not absolute samples, when checking the reduction.
    """
    if not (0.0 <= weight_fast <= 1.0):
        raise ValueError(f"weight_fast must lie in [0, 1], got {weight_fast!r}")
    if clear_slow_ms < clear_fast_ms:
        raise ValueError(
            f"clear_slow_ms ({clear_slow_ms}) is faster than clear_fast_ms "
            f"({clear_fast_ms}); swap them rather than relying on the labels being "
            f"ignored downstream.")
    fast = synaptic_transient(t, peak_um=1.0, rise_ms=rise_ms,
                              clear_ms=clear_fast_ms, t_start=t_start)
    slow = synaptic_transient(t, peak_um=1.0, rise_ms=rise_ms,
                              clear_ms=clear_slow_ms, t_start=t_start)
    mixed = weight_fast * fast + (1.0 - weight_fast) * slow
    pk = float(np.max(mixed)) if mixed.size else 0.0
    if pk <= 1e-12:
        raise ValueError(
            f"the mixed transient has no positive peak for rise_ms={rise_ms!r}, "
            f"clear_fast_ms={clear_fast_ms!r}, clear_slow_ms={clear_slow_ms!r}")
    return np.asarray(peak_um * mixed / pk, dtype=float)


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
    *,
    v_hold_mv: float = DEFAULT_V_HOLD_MV,
    e_cl_mv: float = DEFAULT_E_CL_MV,
) -> dict[str, float]:
    """Standard electrophysiological summary metrics from a `WaveformResult`.

    ONE DRIVING-FORCE CONVENTION, shared with `WaveformResult.current_pA`. This function
    used to take its own `driving_force_mv=30.0` because `current_pA`'s defaults made it
    return zeros, leaving two conventions for one quantity; both now come from
    `models.base` (roadmap P0-11).

    The charge is tagged `CHARGE` with the integration window it was taken over, because a
    charge without its window is not a quantity (`ObservedQuantity`).
    """
    current = res.current_pA(g_max_ns, v_hold_mv=v_hold_mv, e_cl_mv=e_cl_mv)
    peak_current = float(np.max(np.abs(current)))

    _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
    # Charge Q = integral I dt (pA * ms = fC)
    charge_fC = float(_trapz(current, res.t))
    window = float(res.t[-1] - res.t[0]) if res.t.size > 1 else 0.0

    decay = res.decay
    return {
        "peak_open_prob": res.peak_p_open,
        "peak_current_pA": peak_current,
        "charge_integral_fC": charge_fC,
        "charge_window_ms": window,
        "decay_tau_ms": res.decay_tau_ms,
        "decay_tau_fast_ms": decay.tau_fast_ms if decay else float("nan"),
        "decay_tau_slow_ms": decay.tau_slow_ms if decay else float("nan"),
        "decay_weight_fast": decay.weight_fast if decay else float("nan"),
        "decay_r_squared": decay.r_squared if decay else 0.0,
        "steady_state_tail": res.steady_state_tail,
        "driving_force_mv": v_hold_mv - e_cl_mv,
    }


def paired_pulse_ratio(res: WaveformResult, gaba: np.ndarray, *,
                       n_expected: int) -> dict[str, float]:
    """Paired-pulse and steady-state ratios from a pulse-train response.

    THE AMPLITUDE IS MEASURED FROM THE VALUE AT PULSE ONSET, NOT FROM ZERO, and the
    difference is not cosmetic. At 100 Hz with these rates the channel has not closed
    between pulses (deactivation tau is 3.3 ms against a 10 ms period), so the open
    probability at the second onset is 0.477 and the raw peak is mostly residual from the
    first pulse. Measured on the fitted model:

        frequency   peak/peak (from zero)   amplitude/amplitude (from onset)
         10 Hz              0.863                      0.857
         50 Hz              0.900                      0.552
        100 Hz              0.958                      0.282

    The peak-to-zero ratio therefore reports paired-pulse FACILITATION deepening with
    frequency -- 0.86 to 0.96 -- for a scheme that is in fact depressing threefold over the
    same range. It is measuring temporal summation and calling it a response. Experimental
    practice is to subtract the residual; this does the same thing, by onset baseline.

    Both are returned. `ppr_2_over_1` is the amplitude ratio and is the quantity to use;
    `ppr_2_over_1_from_zero` is the old peak-to-zero number, kept so the artefact can be
    demonstrated rather than merely described, and named so it cannot be mistaken for the
    response.

    Pulse onsets are located on the AGONIST trace rather than on the response, because
    under heavy desensitisation a later response peak can vanish entirely and peak-finding
    on `p_open` would then silently return the wrong number of pulses.
    """
    g = np.asarray(gaba, dtype=float)
    p = np.asarray(res.p_open, dtype=float)
    if g.shape != p.shape:
        raise ValueError(f"gaba shape {g.shape} does not match p_open {p.shape}")
    thr = 0.5 * float(np.max(g)) if g.size else 0.0
    if thr <= 0:
        raise ValueError("the agonist trace has no positive amplitude")
    above = g > thr
    # A RISING-EDGE SEARCH CANNOT SEE A PULSE THAT IS ALREADY HIGH AT t=0, so t=0 is
    # treated as an onset when the trace starts above threshold.
    #
    # HARDENING, NOT A LIVE BUG FIX. This cannot be triggered through `pulse_train`, which
    # was the suspected route: even at `start_ms=0.0` the 0.1 ms rise leaves the first
    # sample below the 5% threshold, so `above[0]` is False and both edges are found --
    # verified. It bites only a hand-built waveform, or a train whose rise is shortened to
    # zero. Left unguarded, the symptom is not an exception but a SILENT off-by-one: with
    # `n_expected` matching it raises, but a longer train quietly compares pulse 2 against
    # pulse 3 and reports the wrong ratio.
    onsets = np.flatnonzero(np.diff(above.astype(int)) == 1) + 1
    if above.size and above[0]:
        onsets = np.insert(onsets, 0, 0)
    if onsets.size < 2:
        raise ValueError(
            f"found {onsets.size} pulse onset(s) in the agonist trace, expected "
            f"{n_expected}; the train and the time vector disagree.")
    bounds = list(onsets) + [p.size]
    peaks, amps, baselines = [], [], []
    for i in range(len(bounds) - 1):
        seg = p[bounds[i]:bounds[i + 1]]
        base = float(p[bounds[i]])
        peaks.append(float(np.max(seg)))
        baselines.append(base)
        amps.append(float(np.max(seg)) - base)

    first_peak = peaks[0]
    first_amp = amps[0]
    if first_peak <= 1e-15 or first_amp <= 1e-15:
        raise ValueError(
            f"the first pulse produced no measurable response (peak {first_peak:.3g}, "
            f"amplitude above onset {first_amp:.3g})")
    return {
        "n_pulses_found": float(len(peaks)),
        # the response, measured as an experimenter would
        "ppr_2_over_1": amps[1] / first_amp,
        "steady_over_first": amps[-1] / first_amp,
        "amp_1": first_amp,
        "amp_last": amps[-1],
        "baseline_at_pulse_2": baselines[1],
        "baseline_at_last_pulse": baselines[-1],
        # the summation artefact, kept and labelled
        "ppr_2_over_1_from_zero": peaks[1] / first_peak,
        "steady_over_first_from_zero": peaks[-1] / first_peak,
        "peak_1": first_peak,
        "peak_last": peaks[-1],
    }
