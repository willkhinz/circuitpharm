"""Burst detection that survives a change of neuron substrate.

WHY THIS MODULE EXISTS. `cpg.burst_metrics` reported **5.63 Hz** on a conductance-substrate
preBotC trace whose FFT peak was **0.302 Hz** -- a 19x disagreement -- and it did so silently,
returning a plausible period and duty cycle. Three properties of that detector caused it, and
all three are fine on the LIF cell and wrong on a conductance cell:

  1. ONE THRESHOLD, NO HYSTERESIS. `rate > 0.35 * rate.max()` re-triggers every time the
     intra-burst rate dips. A 21 pF cell with real Na dynamics has far more intra-burst
     structure than a LIF population, and a 20 ms rate low-pass leaves it in.
  2. THRESHOLD REFERENCED TO `max()`. One outlier sample sets the scale for the whole trace.
  3. NO MINIMUM INTER-BURST INTERVAL AND NO MINIMUM DURATION. Nothing stops one burst being
     counted as a dozen.

So the detector measured ripple and called it rhythm. It is the fourth measurement protocol
inherited from the LIF substrate to fail on the conductance cell, after the 4 s warm-up against
tau_h = 10 s, the in vivo validity band, and the FFT search floor. None of the four failed
loudly, which is why the rule in knowledge/05-design-conductance-substrate.md is to assume every
inherited protocol is wrong until re-derived.

WHAT IS DERIVED AND WHAT IS CONVENTIONAL, kept separate on purpose.

DERIVED from the preparation's own validity band, so the detector cannot be tuned per-result:

    T_min      = 1000 / band_hi            shortest period the preparation can plausibly show
    min_gap_ms = (1 - DUTY_MAX) * T_min    shortest plausible SILENT interval; an apparent gap
                                           shorter than this is intra-burst ripple, so merge
    min_dur_ms = DUTY_MIN * T_min          shortest plausible burst; shorter events are noise
    smooth_ms  = SMOOTH_FRAC * T_min       low-pass well inside one burst

CONVENTIONAL, and declared as such because nothing in this project derives them:

    HI_FRAC / LO_FRAC = 0.50 / 0.20 of the p5->p95 range (Schmitt-trigger thresholds)
    DUTY_MIN / DUTY_MAX = 0.05 / 0.70

The conventional four are the ones a result must not depend on, so
`tests/test_bursts.py::test_the_detector_is_insensitive_to_its_conventional_parameters` sweeps
them and requires the burst count to hold.

THE CROSS-CHECK IS PART OF THE RESULT, not an optional extra. Every call returns the FFT peak
beside the detector's own frequency and a `consistent` flag comparing them against the FFT's
resolution. The 19x failure was invisible because nothing compared the two numbers; here they
arrive together and a caller that ignores `consistent` has to ignore it on purpose.
"""
from __future__ import annotations

import numpy as np

#: Schmitt-trigger thresholds, as fractions of the p5->p95 amplitude range. CONVENTIONAL.
HI_FRAC = 0.50
LO_FRAC = 0.20

#: Plausible duty-cycle bounds, used to turn a frequency band into time constants.
#: CONVENTIONAL. In vitro preBotC inspiratory duty cycle is well inside (0.05, 0.70).
DUTY_MIN = 0.05
DUTY_MAX = 0.70

#: Low-pass width as a fraction of the shortest plausible period. CONVENTIONAL.
SMOOTH_FRAC = 0.05

#: Agreement tolerance between the detector and the FFT: the larger of this many FFT bins
#: and this relative difference. Two bins because a peak can sit either side of a bin edge.
FFT_TOL_BINS = 2.0
FFT_TOL_REL = 0.20


def _percentile_range(x):
    """Amplitude reference from percentiles, never from max(). See failure 2 above."""
    p5, p95 = float(np.percentile(x, 5)), float(np.percentile(x, 95))
    return p5, p95


def _fft_peak(t_ms, x, f_lo, f_hi_search):
    """Dominant frequency in Hz, amplitude-invariant, with the bin width that found it."""
    if len(x) < 8:
        return np.nan, np.nan
    dt_s = float(t_ms[1] - t_ms[0]) / 1000.0
    y = (x - x.mean()) * np.hanning(len(x))
    F = np.abs(np.fft.rfft(y))
    fr = np.fft.rfftfreq(len(x), dt_s)
    m = (fr > f_lo) & (fr < f_hi_search)
    if not m.any():
        return np.nan, float(fr[1] - fr[0]) if len(fr) > 1 else np.nan
    return float(fr[m][np.argmax(F[m])]), float(fr[1] - fr[0])


def detect_bursts(t, rate, band, *, hi_frac=HI_FRAC, lo_frac=LO_FRAC,
                  duty_min=DUTY_MIN, duty_max=DUTY_MAX, smooth_frac=SMOOTH_FRAC):
    """Burst times, period, duty cycle and an FFT cross-check.

    `t` in ms, uniformly sampled. `band` is the preparation's (lo, hi) validity band in Hz --
    `config.EUPNOEA_BAND` in vivo, `config.INVITRO_BAND` for the neonatal slice. The band is
    not a filter on the answer; it sets the detector's time constants, which is why a
    substrate change requires passing its own band rather than reusing the default.

    Returns a dict with:
        n_bursts, period_ms, period_cv, duty, onsets, offsets, durations
        freq_hz        detector frequency, 1000 / period_ms
        fft_hz         FFT peak over the same window
        fft_bin_hz     FFT resolution, so the caller can judge a disagreement
        consistent     do the two agree within tolerance?
        reason         why not, when something is wrong
    """
    t = np.asarray(t, float)
    r = np.asarray(rate, float)
    nan = dict(n_bursts=0, period_ms=np.nan, period_cv=np.nan, duty=np.nan,
               onsets=[], offsets=[], durations=[], freq_hz=np.nan, fft_hz=np.nan,
               fft_bin_hz=np.nan, consistent=False)
    if t.size != r.size:
        return dict(nan, reason=f"t and rate differ in length: {t.size} vs {r.size}")
    if t.size < 16:
        return dict(nan, reason=f"trace too short to analyse: {t.size} samples")
    if not np.isfinite(r).all() or not np.isfinite(t).all():
        # Declared, not worked around. A NaN that reaches here means an upstream integration
        # diverged, and reporting a period for it is how a divergence becomes a result.
        return dict(nan, reason="trace contains non-finite samples")

    band_lo, band_hi = float(band[0]), float(band[1])
    if not (0.0 < band_lo < band_hi):
        raise ValueError(f"band must be 0 < lo < hi, got {band!r}")

    dt_ms = float(t[1] - t[0])
    t_min_ms = 1000.0 / band_hi                 # shortest plausible period
    min_gap_ms = (1.0 - duty_max) * t_min_ms    # shortest plausible silence
    min_dur_ms = duty_min * t_min_ms            # shortest plausible burst
    smooth_ms = max(dt_ms, smooth_frac * t_min_ms)

    w = max(1, int(round(smooth_ms / dt_ms)))
    sm = np.convolve(r, np.ones(w) / w, mode="same")
    if sm.size > 3 * w:                          # drop convolution edges, as resp_metrics does
        sm, t = sm[w:-w], t[w:-w]
    if sm.size < 16:
        return dict(nan, reason="trace too short after edge trimming")

    p5, p95 = _percentile_range(sm)
    span = p95 - p5
    mean = float(sm.mean())
    fft_hz, fft_bin = _fft_peak(t, sm, 0.5 * band_lo, max(5.0, 4.0 * band_hi))

    # A flat trace has no bursts, and a relative threshold on one finds spurious ones. The
    # comparison is to the MEAN, so this is scale-free: a quiescent network and a tonically
    # firing one both land here, correctly, with n_bursts = 0.
    if span <= 1e-9 or span < 0.05 * max(1e-9, abs(mean)):
        return dict(nan, fft_hz=fft_hz, fft_bin_hz=fft_bin, consistent=False,
                    reason=f"no amplitude modulation (p95-p5 = {span:.3g}, mean = {mean:.3g})")

    hi_thr, lo_thr = p5 + hi_frac * span, p5 + lo_frac * span

    # ---- Schmitt trigger: rise through hi_thr opens a burst, fall through lo_thr closes it
    above_hi, below_lo = sm > hi_thr, sm < lo_thr
    onsets, offsets, inside = [], [], False
    for i in range(1, len(sm)):
        if not inside and above_hi[i] and not above_hi[i - 1]:
            onsets.append(t[i]); inside = True
        elif inside and below_lo[i] and not below_lo[i - 1]:
            offsets.append(t[i]); inside = False

    if not onsets:
        return dict(nan, fft_hz=fft_hz, fft_bin_hz=fft_bin, consistent=False,
                    reason="no burst onset crossed the upper threshold")

    # ---- merge events separated by less than one plausible silent interval
    # THIS IS THE STEP THAT FIXES THE 19x ERROR. Without it, intra-burst ripple that dips
    # through lo_thr and back ends and restarts a burst.
    merged_on, merged_off = [onsets[0]], []
    for k in range(len(onsets)):
        end = offsets[k] if k < len(offsets) else t[-1]
        nxt = onsets[k + 1] if k + 1 < len(onsets) else None
        if nxt is not None and (nxt - end) < min_gap_ms:
            continue                            # ripple, not a real gap: keep the burst open
        merged_off.append(end)
        if nxt is not None:
            merged_on.append(nxt)
    merged_on = merged_on[:len(merged_off)]

    # ---- discard events shorter than one plausible burst
    keep = [(a, b) for a, b in zip(merged_on, merged_off) if (b - a) >= min_dur_ms]
    if not keep:
        return dict(nan, fft_hz=fft_hz, fft_bin_hz=fft_bin, consistent=False,
                    reason=(f"every candidate burst was shorter than {min_dur_ms:.1f} ms, "
                            f"the shortest plausible burst in a {band_lo}-{band_hi} Hz band"))
    on = np.array([a for a, _ in keep])
    off = np.array([b for _, b in keep])
    durations = off - on

    # The LAST burst may be truncated by the window, which shortens the mean duration and so
    # biases duty DOWN. Dropped from the duration statistics, kept for the period, because an
    # onset is a valid period marker whether or not its offset fits in the window.
    complete = durations[off < t[-1] - 0.5 * dt_ms] if len(durations) else durations
    dur_for_duty = complete if complete.size else durations

    if len(on) > 2:
        iv = np.diff(on)
        period = float(iv.mean())
        cv = float(iv.std() / period) if period > 0 else np.nan
    else:
        period, cv = np.nan, np.nan

    duty = (float(dur_for_duty.mean() / period)
            if dur_for_duty.size and np.isfinite(period) and period > 0 else np.nan)
    freq = 1000.0 / period if np.isfinite(period) and period > 0 else np.nan

    # ---- the cross-check, returned whether or not anyone looks at it
    tol = max(FFT_TOL_BINS * fft_bin, FFT_TOL_REL * fft_hz) if np.isfinite(fft_hz) else np.nan
    consistent = bool(np.isfinite(freq) and np.isfinite(fft_hz) and abs(freq - fft_hz) <= tol)
    reason = "ok" if consistent else (
        "fewer than 3 bursts, so the detector has no period to compare"
        if not np.isfinite(freq) else
        f"detector {freq:.3f} Hz vs FFT {fft_hz:.3f} Hz differs by "
        f"{abs(freq - fft_hz):.3f} Hz, tolerance {tol:.3f} Hz "
        f"(bin {fft_bin:.4f} Hz) -- one of them is measuring something else")

    return dict(n_bursts=int(len(on)), period_ms=period, period_cv=cv, duty=duty,
                onsets=on.tolist(), offsets=off.tolist(), durations=durations.tolist(),
                freq_hz=freq, fft_hz=fft_hz, fft_bin_hz=fft_bin,
                consistent=consistent, reason=reason,
                thresholds=dict(p5=p5, p95=p95, hi=hi_thr, lo=lo_thr),
                derived=dict(t_min_ms=t_min_ms, min_gap_ms=min_gap_ms,
                             min_dur_ms=min_dur_ms, smooth_ms=smooth_ms, window=w))
