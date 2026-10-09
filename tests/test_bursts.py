"""Acceptance tests for `circuitpharm.bursts`, against synthetic ground truth.

WHY SYNTHETIC. The detector exists because the previous one reported 5.63 Hz on a trace whose
FFT peak was 0.302 Hz. Validating a replacement on simulation output would mean comparing two
estimates with no ground truth between them; these traces have a known burst count, frequency
and duty cycle by construction, so "correct" is defined before the detector runs.
"""
import numpy as np
import pytest

from circuitpharm.bursts import (DUTY_MAX, DUTY_MIN, HI_FRAC, LO_FRAC, SMOOTH_FRAC,
                                 detect_bursts)

INVITRO = (0.05, 1.00)
EUPNOEA = (0.30, 2.50)


def burst_train(f=0.30, duty=0.35, ripple=0.45, noise=1.2, T=120000.0, dt=1.0,
                seed=0, hi=30.0, lo=2.0, phase0=0.0, ripple_ms=80.0):
    """A rectangular burst train with intra-burst ripple and additive noise.

    `ripple` is the failure mode: a within-burst oscillation that dips the rate back through
    a single threshold and so re-triggers a naive detector.
    """
    t = np.arange(0.0, T, dt)
    ph = (t / 1000.0 * f + phase0) % 1.0
    base = np.where(ph < duty, hi, lo)
    rng = np.random.default_rng(seed)
    x = base * (1.0 + ripple * np.sin(2 * np.pi * t / ripple_ms))
    x = x + rng.normal(0.0, noise, t.size)
    return t, np.clip(x, 0.0, None)


@pytest.mark.parametrize("f,duty", [(0.08, 0.35), (0.30, 0.15), (0.30, 0.35),
                                    (0.30, 0.60), (0.90, 0.35)])
def test_recovers_known_frequency_and_duty(f, duty):
    t, x = burst_train(f=f, duty=duty)
    d = detect_bursts(t, x, INVITRO)
    assert d["consistent"], d["reason"]
    assert abs(d["freq_hz"] - f) < 0.05 * f, f"frequency {d['freq_hz']:.4f} vs truth {f}"
    assert abs(d["duty"] - duty) < 0.05, f"duty {d['duty']:.3f} vs truth {duty}"


@pytest.mark.parametrize("ripple", [0.0, 0.3, 0.45, 0.6, 0.8])
def test_intra_burst_ripple_does_not_multiply_the_burst_count(ripple):
    """THE REGRESSION. Ripple is what produced the 19x error, so the count must not move
    with it at all."""
    t, x = burst_train(f=0.30, duty=0.35, ripple=ripple)
    d = detect_bursts(t, x, INVITRO)
    assert d["consistent"], f"ripple {ripple}: {d['reason']}"
    assert abs(d["freq_hz"] - 0.30) < 0.015, (
        f"ripple {ripple} moved the detected frequency to {d['freq_hz']:.4f} Hz")


def test_the_old_detector_fails_where_this_one_does_not():
    """Documents why this module exists, and will fail if `burst_metrics` is ever fixed --
    at which point this test should be deleted, not relaxed."""
    from circuitpharm.cpg import burst_metrics
    t, x = burst_train(f=0.30, duty=0.35, ripple=0.45)
    old = burst_metrics(t, x)
    old_hz = 1000.0 / old["period"]
    new = detect_bursts(t, x, INVITRO)
    assert abs(new["freq_hz"] - 0.30) < 0.015, "the new detector should be right here"
    assert old_hz > 3.0, (
        f"burst_metrics now reports {old_hz:.3f} Hz on the ripple case, not a spurious "
        "high frequency -- if it has been fixed, delete this test rather than loosening it")
    # And the specific reason the duty-cycle result from the old detector was nonetheless
    # usable: period and duration inflate together, so their RATIO survives.
    assert abs(old["duty"] - 0.35) < 0.10, (
        "the old detector's duty cycle was approximately right despite a ~40x period error, "
        "because duration and period are both inflated; that is why the duty-cycle finding "
        "was kept when the frequency finding was retracted")


def test_the_detector_is_insensitive_to_its_conventional_parameters():
    """HI_FRAC, LO_FRAC, DUTY_MIN, DUTY_MAX and SMOOTH_FRAC are conventional, not derived.
    A result that depends on them is a result about them."""
    t, x = burst_train(f=0.30, duty=0.35, ripple=0.45)
    base = detect_bursts(t, x, INVITRO)
    assert base["consistent"]
    for kw in (dict(hi_frac=0.40), dict(hi_frac=0.65), dict(lo_frac=0.10),
               dict(lo_frac=0.30), dict(duty_max=0.55), dict(duty_max=0.80),
               dict(duty_min=0.02), dict(duty_min=0.10), dict(smooth_frac=0.02),
               dict(smooth_frac=0.10)):
        d = detect_bursts(t, x, INVITRO, **kw)
        assert d["consistent"], f"{kw}: {d['reason']}"
        assert abs(d["freq_hz"] - base["freq_hz"]) < 0.02 * base["freq_hz"], (
            f"{kw} moved the frequency from {base['freq_hz']:.4f} to {d['freq_hz']:.4f}; "
            "the result depends on a conventional parameter")
        assert abs(d["duty"] - base["duty"]) < 0.06, (
            f"{kw} moved duty from {base['duty']:.3f} to {d['duty']:.3f}")


def test_the_time_constants_are_derived_from_the_band_not_fixed():
    """The band is a property of the preparation. If it stops driving the constants, the
    detector has been decoupled from the substrate and E17 is reachable again."""
    t, x = burst_train(f=0.30)
    slow = detect_bursts(t, x, (0.05, 1.00))["derived"]
    fast = detect_bursts(t, x, (0.30, 2.50))["derived"]
    assert fast["t_min_ms"] < slow["t_min_ms"]
    for k in ("min_gap_ms", "min_dur_ms", "smooth_ms"):
        assert fast[k] < slow[k], f"{k} did not scale with the band"
    # and the documented relationships hold exactly
    assert slow["t_min_ms"] == pytest.approx(1000.0 / 1.00)
    assert slow["min_gap_ms"] == pytest.approx((1.0 - DUTY_MAX) * slow["t_min_ms"])
    assert slow["min_dur_ms"] == pytest.approx(DUTY_MIN * slow["t_min_ms"])
    assert slow["smooth_ms"] == pytest.approx(SMOOTH_FRAC * slow["t_min_ms"])


def test_a_trace_starting_mid_burst_is_not_phase_biased():
    """cpg.burst_metrics returned NaN duty for exactly this case before a prior fix; an
    analysis window opens at an arbitrary phase after settling, so this is the common case,
    not an edge case."""
    out = [detect_bursts(*burst_train(f=0.30, duty=0.35, phase0=p), band=INVITRO)
           for p in (0.0, 0.2, 0.5, 0.8)]
    for d in out:
        assert d["consistent"], d["reason"]
        assert abs(d["duty"] - 0.35) < 0.05, f"duty {d['duty']:.3f} at a shifted phase"
    duties = [d["duty"] for d in out]
    assert max(duties) - min(duties) < 0.04, f"duty varies with window phase: {duties}"


def test_a_truncated_final_burst_does_not_drag_duty_down():
    """A window that closes mid-burst gives one short duration. Including it biases duty
    DOWN, which on a drug arm reads as a shortened burst -- a plausible wrong answer."""
    # choose T so the window closes part-way through a burst
    t, x = burst_train(f=0.30, duty=0.35, T=121200.0)
    d = detect_bursts(t, x, INVITRO)
    assert d["consistent"], d["reason"]
    assert abs(d["duty"] - 0.35) < 0.05, f"duty {d['duty']:.3f} with a truncated last burst"


@pytest.mark.parametrize("name,maker", [
    ("silent", lambda t: np.zeros_like(t)),
    ("tonic constant", lambda t: np.full_like(t, 20.0)),
    ("tonic with small noise",
     lambda t: 20.0 + np.random.default_rng(1).normal(0, 0.3, t.size)),
])
def test_non_bursting_traces_report_no_bursts_with_a_reason(name, maker):
    t = np.arange(0.0, 60000.0, 1.0)
    d = detect_bursts(t, maker(t), INVITRO)
    assert d["n_bursts"] == 0, f"{name}: invented {d['n_bursts']} bursts"
    assert not d["consistent"]
    assert "modulation" in d["reason"], f"{name}: unhelpful reason {d['reason']!r}"


def test_non_finite_input_is_declared_not_worked_around():
    """A NaN here means an upstream integration diverged. Reporting a period for it is how a
    divergence becomes a result -- the project's signature failure shape."""
    t = np.arange(0.0, 60000.0, 1.0)
    x = np.where(t > 30000.0, np.nan, 5.0)
    d = detect_bursts(t, x, INVITRO)
    assert d["n_bursts"] == 0 and not d["consistent"]
    assert "non-finite" in d["reason"]
    assert not np.isfinite(d["period_ms"])


def test_the_fft_cross_check_catches_a_planted_disagreement():
    """The check must be able to fail, or it is decoration. A band whose upper bound is far
    below the true rhythm makes min_gap_ms longer than the real period, so the detector
    merges genuine bursts and must report the disagreement rather than the merge."""
    t, x = burst_train(f=0.90, duty=0.35, ripple=0.0, noise=0.2)
    d = detect_bursts(t, x, (0.02, 0.10))      # wrong band on purpose
    assert not d["consistent"], (
        f"a 0.90 Hz rhythm analysed in a 0.02-0.10 Hz band reported consistent: {d}")
    assert "differs by" in d["reason"] or "fewer than 3" in d["reason"], d["reason"]


def test_period_cv_reports_irregularity():
    """Regularity is reported so a 'rhythm' that is really scattered events cannot pass as
    one on frequency alone."""
    t, x = burst_train(f=0.30, duty=0.35, ripple=0.0, noise=0.2)
    regular = detect_bursts(t, x, INVITRO)
    assert regular["period_cv"] < 0.05, f"a clean train gave cv {regular['period_cv']:.3f}"

    rng = np.random.default_rng(3)
    tt = np.arange(0.0, 240000.0, 1.0)
    y = np.full_like(tt, 2.0)
    centre = 0.0
    while centre < tt[-1] - 4000.0:                 # jittered inter-burst intervals
        centre += rng.uniform(1500.0, 6000.0)
        y[(tt > centre) & (tt < centre + 900.0)] = 30.0
    jittered = detect_bursts(tt, y, INVITRO)
    assert jittered["period_cv"] > regular["period_cv"] * 4, (
        f"jittered train cv {jittered['period_cv']:.3f} vs regular "
        f"{regular['period_cv']:.3f}")
