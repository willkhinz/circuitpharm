"""P6: waveforms, charge with its window, and the discriminating experiment.

Three ANCHORs here, each tying new machinery to something that is true independently of it:

  * `test_charge_ratio_recovers_steady_state_limit` -- as the integration window grows past
    the slowest relaxation, the CHARGE gain must converge to the EQUILIBRIUM gain. A
    time-domain quantity meeting an algebraic one.
  * `test_paired_pulse_depression_sign` -- a desensitisation scheme must depress, and
    deeper at higher frequency. A scheme that does not is mis-wired, and no amount of
    fitting would reveal it.
  * `test_discrimination_score_falls_when_posterior_widens` -- the regression guard for C6.
    The score this replaces was insensitive to parameter uncertainty BY CONSTRUCTION, so a
    test that it is now sensitive is the only thing that stops the old behaviour returning.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from circuitpharm.dynamic_range import (
    dynamic_range_posterior,
    dynamic_range_posterior_report,
    evaluate_dynamic_range,
)
from circuitpharm.fitting.comparison import FIT_CONVENTIONS, MODEL_SPECS
from circuitpharm.models import KineticAllosteryModel
from circuitpharm.models.base import Observable, ObservableMismatch
from circuitpharm.protocols.design import (
    DEFAULT_SIGMA_MEAS,
    ProtocolPoint,
    _penalise,
    design_report,
    find_discriminating_protocol_pp,
    posterior_predictive_score,
    sigma_sensitivity,
)
from circuitpharm.protocols.waveforms import (
    paired_pulse_ratio,
    pulse_train,
    synaptic_transient,
    synaptic_transient_biexp_clearance,
)
from circuitpharm.results import Tier, VoidQuantityError

CONV = {k: v for k, (v, _) in FIT_CONVENTIONS.items() if k in ("koff", "alpha", "r")}


def _model(log10_kd=1.4, log10_e=0.6, log10_d=1.4) -> KineticAllosteryModel:
    return MODEL_SPECS["kinetic_jw95"].factory(
        np.array([log10_kd, log10_e, log10_d], dtype=float))


@pytest.fixture(scope="module")
def draws():
    """A stand-in posterior: tight enough to be informative, wide enough to have a spread.

    NOT a real posterior. These tests are about the machinery that consumes draws, and
    running a 20 s chain in each of them would buy nothing -- `tests/test_p4_likelihood.py`
    is where the sampler itself is checked.
    """
    rng = np.random.default_rng(0)
    return np.column_stack([rng.normal(1.4, 0.06, 150),
                            rng.normal(0.6, 0.08, 150),
                            rng.normal(1.4, 0.09, 150)])


# ================================================================= P6-1 waveforms
def test_biexponential_clearance_reduces_to_single_tau():
    """With both taus equal the mixture must be the single-exponential transient.

    SHAPE, not absolute samples: `synaptic_transient_biexp_clearance` normalises to the
    SAMPLED peak because the mixture's peak has no closed form, so the two differ by the
    ratio of sampled to analytic peak (~0.02% on this grid). That is documented in the
    function and is what this test asserts, rather than a bit-for-bit equality that would
    be asserting the normalisation is absent.
    """
    t = np.arange(0.0, 60.0, 0.025)
    single = synaptic_transient(t, peak_um=1000.0, rise_ms=0.1, clear_ms=5.0)
    mixed = synaptic_transient_biexp_clearance(
        t, peak_um=1000.0, rise_ms=0.1, clear_fast_ms=5.0, clear_slow_ms=5.0,
        weight_fast=0.8)
    assert mixed / mixed.max() == pytest.approx(single / single.max(), abs=1e-9)
    assert mixed.max() == pytest.approx(1000.0, rel=1e-9)
    # and the weight must not change the amplitude, only the decay
    for w in (0.0, 0.3, 1.0):
        m = synaptic_transient_biexp_clearance(t, peak_um=1000.0, weight_fast=w)
        assert m.max() == pytest.approx(1000.0, rel=1e-9), w


def test_biexponential_clearance_refuses_inverted_taus():
    t = np.arange(0.0, 20.0, 0.05)
    with pytest.raises(ValueError, match="faster than"):
        synaptic_transient_biexp_clearance(t, clear_fast_ms=20.0, clear_slow_ms=1.0)
    with pytest.raises(ValueError, match="weight_fast"):
        synaptic_transient_biexp_clearance(t, weight_fast=1.5)


def test_a_slower_clearance_component_carries_more_charge():
    """The point of the second component: it is what a single tau cannot represent."""
    t = np.arange(0.0, 200.0, 0.05)
    trapz = getattr(np, "trapezoid", None) or np.trapz
    fast_only = synaptic_transient_biexp_clearance(t, clear_slow_ms=1.0, clear_fast_ms=1.0)
    with_slow = synaptic_transient_biexp_clearance(t, clear_fast_ms=1.0,
                                                   clear_slow_ms=25.0, weight_fast=0.8)
    assert trapz(with_slow, t) > 2.0 * trapz(fast_only, t)


@pytest.mark.slow
def test_paired_pulse_depression_sign():
    """ANCHOR. A desensitisation scheme depresses, and deeper at higher frequency.

    I2/I1 < 1 at 100 Hz, and the depression must deepen MONOTONICALLY across 10, 50 and
    100 Hz: at higher frequency the receptor has less time to resensitise between pulses,
    so more of the population is in A2D when the next pulse arrives. If this ever comes out
    non-monotone, the desensitisation arm is wired backwards, and that is a defect no
    fitting procedure would surface -- it would simply fit `d` and `r` to compensate.

    WRITING THIS TEST FOUND A MEASUREMENT DEFECT, and the second half of it is the
    regression guard for that. `paired_pulse_ratio` used to measure each peak FROM ZERO. At
    100 Hz the channel has not closed between pulses -- deactivation tau is 3.3 ms against
    a 10 ms period -- so the open probability at the second onset is 0.477 and the raw peak
    is mostly residual. Measured from zero, this scheme reports paired-pulse FACILITATION
    deepening with frequency (0.863 -> 0.900 -> 0.958) while actually depressing threefold
    (0.857 -> 0.552 -> 0.282). The amplitude is now taken from the onset baseline, as an
    experimenter would, and both numbers are returned so the artefact stays demonstrable.
    """
    m = _model()
    ratios = {}
    for freq in (10.0, 50.0, 100.0):
        period = 1000.0 / freq
        t = np.arange(0.0, 10.0 + period * 21.0, 0.05)
        g = pulse_train(t, freq_hz=freq, n_pulses=20, peak_um=1000.0, rise_ms=0.1,
                        clear_ms=1.0, start_ms=10.0)
        res = m.simulate_waveform(t, g, initial_state=m.state_distribution(0.0))
        pp = paired_pulse_ratio(res, g, n_expected=20)
        ratios[freq] = pp
        assert pp["ppr_2_over_1"] < 1.0, (
            f"at {freq:g} Hz the second pulse is not depressed (I2/I1 = "
            f"{pp['ppr_2_over_1']:.4f})")
        assert pp["steady_over_first"] <= pp["ppr_2_over_1"] + 1e-9, (
            f"at {freq:g} Hz depression does not accumulate across the train")

    ppr = [ratios[f]["ppr_2_over_1"] for f in (10.0, 50.0, 100.0)]
    assert ppr[0] > ppr[1] > ppr[2], (
        f"paired-pulse depression is not monotone in frequency: {dict(zip((10, 50, 100), ppr))}")
    assert ppr[2] < 0.5, (
        f"at 100 Hz the depression is only to {ppr[2]:.3f}; it was 0.282 when measured")

    # the artefact, pinned: measured from zero the SAME traces report facilitation
    from_zero = [ratios[f]["ppr_2_over_1_from_zero"] for f in (10.0, 50.0, 100.0)]
    assert from_zero[0] < from_zero[1] < from_zero[2], (
        f"the peak-to-zero ratio no longer rises with frequency ({from_zero}), so the "
        f"summation artefact this measurement exists to avoid has gone -- check whether "
        f"the model's deactivation got faster before relaxing the test")
    assert ratios[100.0]["baseline_at_pulse_2"] > 0.3, (
        "the residual at the second onset is small, so there is nothing for the "
        "peak-to-zero measure to inflate and the comparison above is vacuous")


@pytest.mark.slow
def test_charge_ratio_recovers_steady_state_limit():
    """ANCHOR. As the window grows past the slowest relaxation, CHARGE gain -> EQUILIBRIUM gain.

    For a SUSTAINED application the time-averaged open probability converges to the
    stationary one, so the ratio of charges with and without a PAM must converge to the
    ratio of stationary open probabilities -- a number computed from algebra, not from an
    integration. This is the test that ties the whole CHARGE tier to a result that does not
    depend on the integrator.

    The slowest relaxation here is 1/r = 500 ms, so the window is taken to 5/r = 2500 ms
    and the tolerance is 1%, as the roadmap specifies.
    """
    m = _model()
    pam = 2.5
    expected = (m.steady_state(50.0, pam_factor=pam) / m.steady_state(50.0, pam_factor=1.0))
    window = 5.0 / m.r
    t = np.linspace(0.0, window, 3000)
    g = np.full_like(t, 50.0)
    resting = m.state_distribution(0.0)
    q_ctrl = m.simulate_waveform(t, g, pam_factor=1.0,
                                 initial_state=resting).charge_integral
    q_pam = m.simulate_waveform(t, g, pam_factor=pam,
                                initial_state=resting).charge_integral
    assert q_ctrl > 0
    assert q_pam / q_ctrl == pytest.approx(expected, rel=0.01), (
        f"charge gain {q_pam / q_ctrl:.6f} against equilibrium gain {expected:.6f} over a "
        f"{window:.0f} ms window (5/r)")


# ========================================================= P6-2 tiers with intervals
def test_dynamic_range_posterior_reports_an_interval_per_tier(draws):
    tiers = dynamic_range_posterior(draws, n_draws=30, include_charge=False)
    assert "theoretical_asymptotic_headroom" in tiers
    for name, tp in tiers.items():
        lo, hi = tp.ci_95
        assert lo <= tp.median <= hi, name
        assert tp.n_draws >= 25, (name, tp.n_draws)
        assert tp.spread_factor > 1.0, (
            f"{name} has a {tp.spread_factor} spread, i.e. no interval at all")


def test_the_asymptotic_tier_is_the_one_with_a_wide_interval(draws):
    """The headline number is the least determined one, which is the point of P6-2.

    `theoretical_asymptotic_headroom` is `E/(1+E+D)` over the ambient open probability, so
    it carries the uncertainty in BOTH ratios and `D` is the one no anchor constrains. The
    reachable gains at fixed `s_max` are far better determined because they are ratios of
    the same curve at two K_d values and the E and D dependence largely cancels. A single
    "123x" headline hides exactly this.
    """
    tiers = dynamic_range_posterior(draws, n_draws=30, include_charge=False)
    asymptotic = tiers["theoretical_asymptotic_headroom"].spread_factor
    reachable = tiers["reachable_gain_at_smax_2.5"].spread_factor
    assert asymptotic > 2.0, f"asymptotic spread is only {asymptotic:.3f}-fold"
    assert asymptotic > 10.0 * (reachable - 1.0) + 1.0, (
        f"asymptotic {asymptotic:.3f}-fold against reachable {reachable:.3f}-fold: the "
        f"headline tier is no longer the badly determined one and the P6-2 argument needs "
        f"rewriting")


def test_a_point_estimate_cannot_be_obtained_from_the_interval_function(draws):
    """There is deliberately no way to call this and get a number without an interval."""
    with pytest.raises(ValueError, match="not a posterior"):
        dynamic_range_posterior(draws[:1], n_draws=None, include_charge=False)
    with pytest.raises(ValueError, match=r"\(n, 3\)"):
        dynamic_range_posterior(np.zeros((10, 4)), include_charge=False)


def test_the_conventions_are_required_and_named(draws):
    with pytest.raises(ValueError, match="koff"):
        dynamic_range_posterior(draws, conventions={"alpha": 0.3, "r": 0.002},
                                include_charge=False)
    # and the defaults are the declared fit conventions, not invented here
    assert set(CONV) == {"koff", "alpha", "r"}
    assert CONV["alpha"] == pytest.approx(0.30)


def test_a_void_posterior_makes_every_tier_void(draws):
    rs = dynamic_range_posterior_report(
        draws, posterior_tier=Tier.VOID, posterior_note="chain failed split-Rhat",
        n_draws=20, include_charge=False)
    assert rs.items
    for q in rs.items:
        assert q.tier is Tier.VOID, q.name
        with pytest.raises(VoidQuantityError):
            q.value
    text = str(rs)
    assert "VOID" in text


def test_thinning_takes_a_stride_not_a_prefix():
    """A prefix of a chain may be burn-in that was not discarded."""
    # a chain that drifts: the first 50 draws are nothing like the last 50
    drifting = np.column_stack([np.linspace(0.5, 2.0, 100),
                                np.full(100, 0.6), np.full(100, 1.4)])
    thinned = dynamic_range_posterior(drifting, n_draws=10, include_charge=False)
    prefix = dynamic_range_posterior(drifting[:10], n_draws=None, include_charge=False)
    a = thinned["theoretical_asymptotic_headroom"]
    b = prefix["theoretical_asymptotic_headroom"]
    assert a.spread_factor > 2.0 * b.spread_factor, (
        f"thinned spread {a.spread_factor:.3f} against prefix {b.spread_factor:.3f}: the "
        f"thinning is taking a prefix, which would hide an undiscarded burn-in")


def test_the_point_estimate_function_still_reproduces_its_own_number(draws):
    """`evaluate_dynamic_range` stays as it was, so existing results stay reproducible."""
    m = _model()
    ev = evaluate_dynamic_range(m, ambient_gaba_um=0.40)
    assert ev.ambient_gaba_um == pytest.approx(0.40)
    assert ev.theoretical_asymptotic_headroom > 1.0
    assert ev.collapse_scan, "the collapse boundary must still be a curve"
    # the point value must sit inside the posterior interval when the draws are centred
    # on the same parameters -- if it does not, one of the two is computing a different
    # quantity and the P6-2 rewrite silently changed the metric.
    tiers = dynamic_range_posterior(draws, n_draws=40, include_charge=False)
    lo, hi = tiers["theoretical_asymptotic_headroom"].ci_95
    assert lo <= ev.theoretical_asymptotic_headroom <= hi, (
        f"point estimate {ev.theoretical_asymptotic_headroom:.4f} is outside the "
        f"posterior interval [{lo:.4f}, {hi:.4f}] computed at the same centre")


# ======================================================= P6-3 the discriminating protocol
@pytest.fixture(scope="module")
def candidates(draws):
    rng = np.random.default_rng(7)
    c_draws = np.column_stack([rng.normal(1.4, 0.06, 150), rng.normal(0.6, 0.08, 150),
                               rng.normal(1.0, 0.12, 150), rng.normal(0.5, 0.15, 150)])
    return {"kinetic_jw95": (MODEL_SPECS["kinetic_jw95"].factory, draws),
            "extended_desens": (MODEL_SPECS["extended_desens"].factory, c_draws)}


def test_discrimination_score_falls_when_posterior_widens(candidates, draws):
    """ANCHOR, and the regression guard for C6.

    The score this replaces was `|y_A - y_B| / sigma_noise`: a difference of POINT
    predictions over a fixed measurement noise, which cannot see parameter uncertainty at
    all and would return exactly the same number for both arms of this test. Inflating both
    posteriors 10x about their means must drop the score, because a protocol only
    discriminates if the models disagree by more than either of them is unsure.
    """
    point = ProtocolPoint(gaba_um=100.0, pam_factor=2.5,
                          observable=Observable.EQUILIBRIUM)
    narrow, _, sd_narrow, _ = posterior_predictive_score(candidates, point)

    wide = {}
    for name, (factory, d) in candidates.items():
        d = np.atleast_2d(d)
        wide[name] = (factory, d.mean(axis=0) + (d - d.mean(axis=0)) * 10.0)
    broad, _, sd_broad, _ = posterior_predictive_score(wide, point)

    assert broad < narrow, f"score rose from {narrow:.4f} to {broad:.4f} on a wider posterior"
    assert broad < 0.6 * narrow, (
        f"score only fell from {narrow:.4f} to {broad:.4f} on a 10x wider posterior, which "
        f"is too little to be seeing parameter uncertainty")
    for name in sd_narrow:
        assert sd_broad[name] > sd_narrow[name], name


def test_the_score_needs_exactly_two_models(candidates):
    with pytest.raises(ValueError, match="pairwise"):
        posterior_predictive_score(
            {k: v for i, (k, v) in enumerate(candidates.items()) if i == 0},
            ProtocolPoint(1.0, 1.0, Observable.EQUILIBRIUM))
    three = dict(candidates)
    three["another"] = next(iter(candidates.values()))
    with pytest.raises(ValueError, match="pairwise"):
        posterior_predictive_score(three, ProtocolPoint(1.0, 1.0, Observable.EQUILIBRIUM))


def test_zero_measurement_noise_is_refused(candidates):
    """It would make every protocol infinitely discriminating."""
    with pytest.raises(ValueError, match="sigma_meas"):
        posterior_predictive_score(candidates,
                                   ProtocolPoint(1.0, 1.0, Observable.EQUILIBRIUM),
                                   sigma_meas=0.0)


def test_the_search_carries_its_multiple_comparison_cost():
    """The best of n protocols is biased upward, and the correction must be applied."""
    assert _penalise(3.0, 1) == pytest.approx(3.0)
    assert _penalise(3.0, 100) < 3.0
    assert _penalise(3.0, 10000) < _penalise(3.0, 100)
    # a score the search can fully explain corrects to zero, which is the finding
    assert _penalise(1.5, 500) == 0.0
    # and a genuinely large score survives a large search
    assert _penalise(8.0, 1000) > 5.0


def test_the_design_reports_an_undiscriminating_best_as_such(candidates):
    """If the best protocol in the space cannot separate the models, that is the result."""
    res = find_discriminating_protocol_pp(candidates,
                                          observables=(Observable.EQUILIBRIUM,))
    assert res.n_searched > 10
    assert res.score_penalised <= res.score
    if not res.discriminating:
        assert "DOES NOT CLEAR" in res.note
    assert res.runner_up is not None
    # the falsification intervals are two-sided and include the measurement noise
    for name, (lo, med, hi) in res.predictive_intervals.items():
        assert lo < med < hi, name
        assert hi - lo > 2.0 * DEFAULT_SIGMA_MEAS, (
            f"{name}'s interval is narrower than the measurement noise alone, so it cannot "
            f"be including parameter uncertainty")


def test_an_equilibrium_design_refuses_a_peak_native_model(draws):
    """Model A has no equilibrium response; a design must say so rather than use its peak."""
    rng = np.random.default_rng(3)
    a_draws = np.column_stack([rng.normal(1.1, 0.05, 40), rng.normal(0.15, 0.05, 40),
                               rng.normal(-0.12, 0.03, 40)])
    cand = {"operational": (MODEL_SPECS["operational"].factory, a_draws),
            "kinetic_jw95": (MODEL_SPECS["kinetic_jw95"].factory, draws)}
    with pytest.raises(ObservableMismatch, match="no equilibrium response"):
        posterior_predictive_score(cand, ProtocolPoint(10.0, 2.5,
                                                       Observable.EQUILIBRIUM))
    # and the search skips that observable rather than scoring it as zero separation
    res = find_discriminating_protocol_pp(
        cand, observables=(Observable.EQUILIBRIUM, Observable.PEAK),
        gaba_grid=(1.0, 10.0, 100.0), pam_grid=(1.0, 2.5))
    assert res.protocol.observable is Observable.PEAK
    assert "skipped" in res.note


def test_a_charge_design_is_refused_without_a_window(candidates):
    with pytest.raises(ObservableMismatch, match="window"):
        posterior_predictive_score(candidates,
                                   ProtocolPoint(10.0, 2.5, Observable.CHARGE))


def test_sigma_sensitivity_is_reported_not_assumed(candidates):
    """`sigma_meas` is the only assumed input, so the score is reported against a sweep."""
    point = ProtocolPoint(gaba_um=100.0, pam_factor=2.5,
                          observable=Observable.EQUILIBRIUM)
    sweep = sigma_sensitivity(candidates, point)
    assert len(sweep) >= 3
    ordered = [sweep[s] for s in sorted(sweep)]
    assert ordered == sorted(ordered, reverse=True), (
        f"the score must fall as measurement noise rises: {sweep}")


def test_a_void_posterior_makes_the_design_void(candidates):
    res = find_discriminating_protocol_pp(candidates,
                                          observables=(Observable.EQUILIBRIUM,),
                                          gaba_grid=(1.0, 10.0), pam_grid=(1.0, 2.5))
    rs = design_report(res, posterior_tier=Tier.VOID,
                       posterior_note="chain failed split-Rhat")
    for q in rs.items:
        assert q.tier is Tier.VOID, q.name
        with pytest.raises(VoidQuantityError):
            q.value

    ok = design_report(res, posterior_tier=Tier.UNCALIBRATED)
    for q in ok.items:
        assert q.tier is Tier.UNCALIBRATED, q.name
    joined = " ".join(c for q in ok.items for c in q.caveats)
    assert "biased upward" in joined
    assert "sigma_meas" in joined
    assert "not interchangeable" in joined


# ============================================ the pre-registration's own numbers (P6-4)
@pytest.fixture(scope="module")
def real_posteriors():
    """The two posteriors `knowledge/13-preregistration.md` is computed from.

    ~35 s for both chains. Module-scoped, because three tests below read them and each one
    sampling its own would triple that for no extra assurance.
    """
    from circuitpharm.fitting.data import equilibrium_crc_from_model
    from circuitpharm.fitting.posterior import sample_posterior_vector

    spec_b = MODEL_SPECS["kinetic_jw95"]
    centre = np.array([0.5 * (lo + hi) for lo, hi in spec_b.bounds])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = equilibrium_crc_from_model(
            spec_b.factory(centre), concs_um=np.logspace(-1.0, 3.5, 30), noise_sd=1e-3,
            seed=0, label="P6 design basis, truth = Model B")
        out = {}
        for name in ("kinetic_jw95", "extended_desens"):
            spec = MODEL_SPECS[name]
            out[name] = sample_posterior_vector(
                spec.factory, [ds], param_names=spec.param_names,
                bounds=dict(zip(spec.param_names, spec.bounds)), caller=f"prereg/{name}",
                n_steps=12000, burn_in=4000, seed=11)
    return ds, out


@pytest.mark.slow
def test_model_b_posterior_converges_and_model_c_does_not(real_posteriors):
    """ANCHOR for §1 and §2 of the pre-registration, and the point of the whole document.

    Model C's fourth parameter is not identifiable on an equilibrium CRC generated from
    Model B -- `D_slow` is genuinely zero there -- so its chain cannot mix along that
    direction and the diagnostics must refuse to summarise it. The chain is not broken; the
    parameter is not there. If this ever passes its diagnostics, either the data changed or
    the gate stopped working, and §2's conclusion has to be rewritten either way.
    """
    _, posts = real_posteriors
    b, c = posts["kinetic_jw95"], posts["extended_desens"]

    assert b.converged, b.diagnostics["failures"]
    assert b.tier is Tier.UNCALIBRATED
    assert not any(any(p) for p in b.pinned_at_bound.values())

    assert not c.converged
    assert c.tier is Tier.VOID
    for name in c.param_names:
        with pytest.raises(VoidQuantityError):
            c.intervals_95[name].value
    assert c.pinned_at_bound["log10_D_slow"][0], (
        "log10_D_slow is no longer on its lower bound; the data was generated without a "
        "slow sink, so if the fit now finds one, check the generator before the sampler")
    lo, hi = c.intervals_95["log10_D_slow"].get(acknowledge_void=True)
    assert hi - lo > 2.0, (
        f"log10_D_slow's interval is {hi - lo:.3f} decades; it was 3.1, and a narrow one "
        f"would mean the unidentifiable direction had become identifiable")
    # the shared parameters must nevertheless agree between the two models: they are the
    # same ratios and the data is the same, so a disagreement would be a wiring error
    for shared in ("log10_kd", "log10_E"):
        assert (b.medians[shared].value
                == pytest.approx(c.medians[shared].get(acknowledge_void=True), abs=0.01)), shared


@pytest.mark.slow
def test_the_designed_experiment_does_not_exist_at_conventional_noise(real_posteriors):
    """ANCHOR for §3. The corrected score must be below 1.96, and the document says so.

    This is the pre-registered FINDING, not a limitation of the code: 54 protocols were
    searched, the best raw score is 1.71, and the best of 54 standard normals sits near
    2.3 sigma, so the whole of it is explained by the search. A version of this test that
    asserted the opposite would be asserting that the project has a discriminating
    experiment, which it does not.
    """
    _, posts = real_posteriors
    cand = {n: (MODEL_SPECS[n].factory, posts[n].log10_samples[::200]) for n in posts}
    res = find_discriminating_protocol_pp(cand, observables=(Observable.EQUILIBRIUM,))

    assert res.n_searched == 54, f"the search space changed ({res.n_searched} protocols)"
    assert 1.0 < res.score < 2.5, res.score
    assert res.score_penalised < 1.96
    assert not res.discriminating
    assert "DOES NOT CLEAR" in res.note
    # a plateau, not a peak: the specific concentration is not the design
    assert res.runner_up is not None
    assert res.runner_up[1] > 0.9 * res.score, (
        f"the runner-up scores {res.runner_up[1]:.3f} against {res.score:.3f}; the maximum "
        f"is now a peak rather than a plateau and §3's caveat needs rewriting")
    # and the two predictive intervals must overlap substantially -- that IS the finding
    (lo_b, _, hi_b) = res.predictive_intervals["kinetic_jw95"]
    (lo_c, _, hi_c) = res.predictive_intervals["extended_desens"]
    overlap = max(0.0, min(hi_b, hi_c) - max(lo_b, lo_c))
    assert overlap > 0.4 * min(hi_b - lo_b, hi_c - lo_c), (
        f"the intervals overlap by only {overlap:.4f}; if they have separated, the "
        f"experiment now exists and §3 is wrong")


@pytest.mark.slow
def test_lower_measurement_noise_is_what_would_make_it_exist(real_posteriors):
    """ANCHOR for §5. The binding constraint is sigma_meas, not the concentration.

    The score runs 5.90 / 2.99 / 1.71 / 1.24 at sigma = 0.005 / 0.01 / 0.02 / 0.05, so a
    four-fold noise reduction crosses the threshold while no choice of agonist or modulator
    in the searched space does. THIS IS THE CONCLUSION THE OLD SCORE COULD NOT REACH:
    `|y_B - y_C| / sigma_noise` has no posterior variance in its denominator, so it ranks
    protocols identically at every noise level and can never identify the noise as the
    constraint.
    """
    from circuitpharm.protocols.design import _penalise

    _, posts = real_posteriors
    cand = {n: (MODEL_SPECS[n].factory, posts[n].log10_samples[::200]) for n in posts}
    res = find_discriminating_protocol_pp(cand, observables=(Observable.EQUILIBRIUM,))
    sweep = sigma_sensitivity(cand, res.protocol, sigmas=(0.005, 0.01, 0.02, 0.05))

    assert sweep[0.005] > 4.0, sweep
    assert sweep[0.05] < 2.0, sweep
    assert _penalise(sweep[0.005], 54) >= 1.96, (
        f"a score of {sweep[0.005]:.3f} at sigma 0.005 no longer survives the "
        f"multiple-comparison correction, so §5's actionable requirement is wrong")
    assert _penalise(sweep[0.02], 54) < 1.96


@pytest.mark.slow
def test_posterior_intervals_are_calibrated_across_seeds():
    """A 95% interval cannot be validated against one draw, so this uses six.

    THE SINGLE-SEED VERSION OF THIS TEST WAS PASSING ON LUCK AND IS RECORDED SO IT IS NOT
    RESTORED. On the pre-registration's own noise realisation (seed 0) all three intervals
    miss the truth, on the same side -- which is ONE event, not three, because the three
    ratios are strongly correlated and the posterior is displaced as a whole. Measured
    across six seeds: five cover all three, 15 of 18 marginals, 83%. One miss in six has
    probability 26% at a true 95%, so that is a calibrated interval, and a test asserting
    six-for-six would fail a quarter of the time for no reason.

    The threshold is four of six: loose enough not to be flaky, tight enough that a
    sampler whose intervals were half the right width (and would cover ~2 of 6) fails.
    """
    from circuitpharm.fitting.data import equilibrium_crc_from_model
    from circuitpharm.fitting.posterior import sample_posterior_vector

    spec = MODEL_SPECS["kinetic_jw95"]
    truth = np.array([0.5 * (lo + hi) for lo, hi in spec.bounds])
    joint_hits, converged = 0, 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for seed in range(6):
            ds = equilibrium_crc_from_model(
                spec.factory(truth), concs_um=np.logspace(-1.0, 3.5, 30), noise_sd=1e-3,
                seed=seed)
            p = sample_posterior_vector(
                spec.factory, [ds], param_names=spec.param_names,
                bounds=dict(zip(spec.param_names, spec.bounds)), caller="calibration",
                n_steps=12000, burn_in=4000, seed=11)
            converged += bool(p.converged)
            hits = []
            for i, name in enumerate(spec.param_names):
                lo, hi = p.intervals_95[name].get(acknowledge_void=True)
                hits.append(lo <= truth[i] <= hi)
            joint_hits += all(hits)

    assert converged >= 5, (
        f"only {converged} of 6 chains passed their diagnostics at 12000 steps; the "
        f"calibration result below is then about chains that did not converge")
    assert joint_hits >= 4, (
        f"only {joint_hits} of 6 noise realisations had all three 95% intervals cover the "
        f"truth; 5 of 6 was measured, and 4 is the floor below which the intervals are too "
        f"narrow rather than unlucky")
