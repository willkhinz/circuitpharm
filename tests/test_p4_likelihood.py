"""P4: a likelihood each dataset's own observable can be scored against, and two
independent uncertainty machineries that have to agree.

The three ANCHOR tests here are the ones that would catch a silent regression in the
inference layer:

  * `test_posterior_recovers_known_parameters` -- the sampler, run on data generated at
    known parameters, must cover them. A posterior that does not recover a truth it was
    given cannot be trusted on a truth it was not.
  * `test_posterior_median_in_profile_ci` -- profile likelihood and MCMC share no
    optimiser, no grid and no convergence criterion. If they disagree about where a
    parameter is, one of them is wrong, and neither passing its own internal checks would
    reveal it.
  * `test_mcmc_diagnostics_gate_the_result` -- a deliberately under-run chain must return
    VOID, not a narrow interval. This is the behaviour P0-5 exists for.

The chains take ~20 s each, so the sampling tests are marked `slow`. They are not
optional: `-m "not slow"` runs the cheap half, and CI runs everything.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from circuitpharm.fitting.data import (
    JAHN1997_PEAK_CRC,
    equilibrium_crc_from_model,
)
from circuitpharm.fitting.identifiability import IDENTIFIABLE_BOUNDS
from circuitpharm.fitting.likelihood import (
    concentrated_log_likelihood,
    fit_mle,
    make_log_likelihood,
)
from circuitpharm.fitting.mcmc import sample_ensemble
from circuitpharm.fitting.posterior import (
    agreement_with_profiles,
    posterior_report,
    sample_identifiable_posterior,
)
from circuitpharm.fitting.reparam import EQUILIBRIUM_IDENTIFIABLE, IdentifiableParams
from circuitpharm.models import (
    ExtendedDesensitizationModel,
    KineticAllosteryModel,
    OperationalScalarModel,
)
from circuitpharm.models.base import (
    PEAK_APPLICATION_MS,
    Observable,
    ObservableMismatch,
)
from circuitpharm.results import Tier, VoidQuantityError

#: Known-truth model for the method checks. K_d = 25 uM, E = 4, D = 25.
TRUTH = KineticAllosteryModel(kon=0.01, koff=0.25, beta=0.8, alpha=0.2, d=0.05, r=0.002)
TRUTH_VEC = IdentifiableParams.from_microscopic(**TRUTH.get_params()).as_vector(
    with_timescale=False)

#: The conventions the three identifiable ratios have to be combined with to get rates.
#: They are the TRUTH's own values here, which is legitimate for a method check and
#: nothing else: `to_microscopic` refuses to invent them precisely so that a call site has
#: to say where they came from.
CONV = dict(koff=TRUTH.koff, alpha=TRUTH.alpha, r=TRUTH.r)


def factory(p: IdentifiableParams) -> KineticAllosteryModel:
    return KineticAllosteryModel(**p.to_microscopic(**CONV))


@pytest.fixture(scope="module")
def quiet():
    """Silence the synthetic-data warning. Asserted separately in `test_taint_warns`."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


@pytest.fixture(scope="module")
def dataset():
    """0.1% noise on 30 concentrations, with `sem` EQUAL to the noise that was added.

    The sem matters: a chi-squared profile divides by it while the posterior estimates the
    scale from the residuals, so a sem that misstates the noise makes the two machineries
    disagree for a reason that is about neither of them. At a 1e-3 sem against 2e-4 of
    actual noise the profile intervals came out five times too wide.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return equilibrium_crc_from_model(
            TRUTH, concs_um=np.logspace(-1.0, 3.5, 30), noise_sd=0.001, seed=0)


@pytest.fixture(scope="module")
def posterior(dataset):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sample_identifiable_posterior(factory, [dataset], seed=11)


# ===================================================== the exact peak propagator (P4-0)
def test_scaled_expm_matches_scipy():
    """The hand-rolled scaling-and-squaring must match the library it replaced.

    `scipy.linalg.expm` is correct and costs 8 ms on a 5x5, essentially all overhead, which
    is why the inner loop does not call it. That trade is only acceptable with this test:
    checked on the generators this project builds, at agonist from 0 to 100 mM and over the
    three step sizes `peak_dose_response` produces. 1e-9 absolute against a measured
    worst case of 2.6e-10.
    """
    from scipy.linalg import expm

    from circuitpharm.models.base import _expm_scaled

    for model in (KineticAllosteryModel(), ExtendedDesensitizationModel(), TRUTH):
        for g in (0.0, 0.3, 10.0, 3000.0, 1e5):
            for dt in (0.5, 5.0, 300.0):
                a = np.asarray(model.q_matrix(g), dtype=float).T * dt
                assert _expm_scaled(a) == pytest.approx(expm(a), abs=1e-9), (
                    f"{type(model).__name__} at {g} uM, dt {dt} ms")


def test_exact_peak_agrees_with_the_ode_it_replaced():
    """ANCHOR. The matrix-exponential peak must match `solve_ivp` to the solver's own error.

    This is the one test that licenses a 300x speedup: if the fast path and the slow path
    ever disagree by more than the integrator's tolerance, every PEAK number in the project
    moved for a reason that is not physics. 1e-6 absolute is chosen against the Radau
    settings (rtol 1e-5, atol 1e-7), and the measured difference is 2.2e-7.
    """
    concs = np.array([0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0])
    for model in (KineticAllosteryModel(), ExtendedDesensitizationModel(), TRUTH):
        fast = np.asarray(model.peak_dose_response(concs), dtype=float)
        slow = []
        for c in concs:
            t = np.linspace(0.0, PEAK_APPLICATION_MS, 600)
            g = np.full_like(t, float(c))
            slow.append(model.simulate_waveform(
                t, g, initial_state=model.state_distribution(0.0)).peak_p_open)
        assert fast == pytest.approx(np.array(slow), abs=1e-6), type(model).__name__


def test_the_six_state_generator_matches_its_own_rhs():
    """`ExtendedDesensitizationModel.q_matrix` duplicates `simulate_waveform`'s `rhs`.

    The duplication is deliberate (building a 6x6 per ODE step costs more than it saves),
    so it needs a test that it cannot drift. Checked as `dP/dt = P Q` against a one-step
    finite difference of the integrator at several concentrations.
    """
    m = ExtendedDesensitizationModel()
    for g in (0.0, 1.0, 50.0, 3000.0):
        q = m.q_matrix(g)
        assert q.sum(axis=1) == pytest.approx(np.zeros(6), abs=1e-12), (
            f"rows of a generator must sum to zero; at {g} uM they do not")
        p = m.state_distribution(5.0)
        dt = 1e-6
        t = np.array([0.0, dt])
        res = m.simulate_waveform(t, np.full(2, float(g)), initial_state=p)
        numeric = (res.states[1] - res.states[0]) / dt
        assert numeric == pytest.approx(p @ q, abs=2e-3), f"at {g} uM"


# ========================================================== the likelihood itself (P4-1)
def test_each_dataset_is_scored_on_its_own_observable(quiet):
    """A PEAK dataset must be compared against the PEAK curve, not the equilibrium one.

    The two differ by a factor of ~3 in EC50 at identical parameters, so getting this wrong
    is not a rounding error -- it drives D to its bound by deleting desensitisation (see
    knowledge/11-identifiability.md §2.2). Verified by scoring the SAME model against the
    same concentrations tagged two different ways and requiring the residuals to differ.
    """
    concs = JAHN1997_PEAK_CRC.concs_um
    peak_ll = concentrated_log_likelihood(TRUTH, [JAHN1997_PEAK_CRC], n_params=3)
    eq_version = equilibrium_crc_from_model(TRUTH, concs_um=concs, noise_sd=0.0)
    eq_ll = concentrated_log_likelihood(TRUTH, [eq_version], n_params=3)

    assert peak_ll.per_dataset[0].observable is Observable.PEAK
    assert eq_ll.per_dataset[0].observable is Observable.EQUILIBRIUM
    assert peak_ll.per_dataset[0].rss != pytest.approx(eq_ll.per_dataset[0].rss, rel=1e-6)


def test_a_charge_dataset_is_refused_not_approximated(quiet):
    """A time-course observable needs a waveform and a window; it is not a CRC."""

    class FakeCharge:
        citation_label = "a charge dataset"
        observable = Observable.CHARGE
        concs_um = np.array([1.0, 10.0])
        mean_response = np.array([0.1, 0.5])
        synthetic = True
        role = "train"

    with pytest.raises(ObservableMismatch, match="CHARGE"):
        concentrated_log_likelihood(TRUTH, [FakeCharge()], n_params=3)


def test_sigma_is_estimated_so_no_result_depends_on_a_hand_set_noise_level(quiet):
    """The profiled scale must equal sqrt(RSS/n) exactly, and the constants must be there.

    If `sigma_hat` were a fixed input, every AIC difference and Akaike weight in P5 would
    inherit it (roadmap C5). Checked against the closed form rather than a regression value.
    """
    ds = equilibrium_crc_from_model(TRUTH, noise_sd=0.01, seed=3)
    off = KineticAllosteryModel(kon=0.02, koff=0.25, beta=0.8, alpha=0.2, d=0.05, r=0.002)
    r = concentrated_log_likelihood(off, [ds], n_params=3)
    fit = r.per_dataset[0]
    assert fit.sigma_hat == pytest.approx(np.sqrt(fit.rss / fit.n), rel=1e-12)
    expected = -0.5 * fit.n * (np.log(2.0 * np.pi) + 1.0 + np.log(fit.rss / fit.n))
    assert fit.log_likelihood == pytest.approx(expected, rel=1e-12)
    # n_estimated counts the three parameters AND the scale that was estimated for the
    # dataset. P5's k is this number, not the length of a parameter dataclass.
    assert r.n_estimated == 4


def test_the_holdout_is_unreachable_from_a_fit():
    """P5's out-of-sample score and P6's falsification bound both die to a leak."""

    class Held:
        citation_label = "held out"
        observable = Observable.EQUILIBRIUM
        concs_um = np.array([1.0, 10.0])
        mean_response = np.array([0.1, 0.5])
        synthetic = False
        role = "holdout"

    with pytest.raises(ValueError, match="holdout"):
        concentrated_log_likelihood(TRUTH, [Held()], n_params=3)
    with pytest.raises(ValueError, match="holdout"):
        make_log_likelihood([Held()], caller="t", n_params=3)


def test_taint_warns_once_per_fit_not_once_per_iteration():
    """The guards run when the closure is built, so an optimiser does not bury the warning.

    `provisional.py` sets `simplefilter("always")` for its category precisely so the
    warning cannot be shown once and suppressed. That makes it the inner loop's problem:
    at one warning per likelihood evaluation a 10,000-evaluation fit emits 10,000 of them
    and the one that matters is unreadable. One per fit is the fix, and it is a behaviour
    worth pinning.
    """
    from circuitpharm.provisional import ProvisionalResultWarning

    ds = equilibrium_crc_from_model(TRUTH, noise_sd=0.0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        score = make_log_likelihood([ds], caller="test", n_params=3)
        for _ in range(25):
            score(TRUTH)
    provisional = [w for w in caught if issubclass(w.category, ProvisionalResultWarning)]
    assert len(provisional) == 1, (
        f"{len(provisional)} warnings for one fit of 25 evaluations")
    assert "generated EQUILIBRIUM CRC" in str(provisional[0].message)


# =================================================================== the MLE (P4-2)
def test_mle_recovers_known_parameters_without_noise(quiet):
    """ANCHOR, and the estimator's method check. Noiseless data pins the truth exactly."""
    ds = equilibrium_crc_from_model(TRUTH, noise_sd=0.0)
    r = fit_mle(factory, [ds], n_starts=16, seed=0)
    assert r.params.as_vector(with_timescale=False) == pytest.approx(TRUTH_VEC, abs=1e-3)
    assert r.converged_uniquely, r.note


def test_mle_reports_a_spread_rather_than_picking_a_winner(quiet):
    """The top-5 spread must be reported, and `converged_uniquely` must reflect it."""
    ds = equilibrium_crc_from_model(TRUTH, noise_sd=0.0)
    r = fit_mle(factory, [ds], n_starts=12, seed=1)
    assert set(r.top_spread) == set(EQUILIBRIUM_IDENTIFIABLE)
    assert r.n_converged >= 1
    assert r.converged_uniquely == all(v <= 0.05 for v in r.top_spread.values())
    assert ("agree to within" in r.note) == r.converged_uniquely


@pytest.mark.slow
def test_mle_reproduces_what_the_source_actually_measured(quiet):
    """ANCHOR. Fitted to the Jahn PEAK curve, Model B reproduces the two MEASURED
    quantities -- and the things it does to achieve the third are the reconstruction's
    artefact, not a result.

    WHAT THE SOURCE MEASURES is an EC50 of 11.6 +/- 0.9 uM and a slope of 2.2 +/- 0.4
    "between 0.001 and 0.01 mM GABA" -- a LOCAL slope over 1-10 uM. It does not measure a
    whole-curve Hill coefficient, so this test no longer asks for one. The earlier version
    did, and in passing it was certifying an artefact: see below.

    Measured here: EC50 11.629 uM, rising-phase slope 2.007 -- both inside the published
    errors, and the rising-phase slope is the like-for-like comparison.

    THE ARTEFACT, asserted so it stays visible rather than being read as a finding. The
    nine points of this dataset are a GLOBAL Hill curve with nH = 2.2 spanning four
    decades, reconstructed from a slope measured over one (see `fitting/data.py`). A curve
    that steep that far out can only be matched by a scheme that has deleted
    desensitisation, so the fit puts `log10_D` on its lower bound and the absolute plateau
    at 0.995 against this project's 0.750 anchor. That is a property of the
    RECONSTRUCTION, not of the receptor and not of the scheme -- which is why
    `knowledge/12-inference.md` Section 2's "structural conflict" was withdrawn.

    THE CONVENTIONS ARE THE TRUTH'S OWN HERE. `koff`, `alpha` and `r` are held at `CONV`,
    so this fit determines three ratios and not six rates.
    """
    from circuitpharm.models.base import JAHN1997_SLOPE_WINDOW_REL, hill_slope

    r = fit_mle(factory, [JAHN1997_PEAK_CRC], n_starts=8, seed=0,
                bounds=[IDENTIFIABLE_BOUNDS[k] for k in EQUILIBRIUM_IDENTIFIABLE])
    fitted = factory(r.params)
    concs = np.logspace(-2.0, 4.5, 300)
    y = np.asarray(fitted.peak_dose_response(concs), dtype=float)
    ec50 = float(np.interp(0.5, y / y.max(), concs))

    # the two MEASURED quantities, each compared the way it was measured
    assert ec50 == pytest.approx(11.6, abs=1.8), (
        f"fitted PEAK EC50 {ec50:.3f} uM against Jahn's 11.6 +/- 0.9")
    rising = hill_slope(concs, y, window_rel=JAHN1997_SLOPE_WINDOW_REL)
    assert rising == pytest.approx(2.2, abs=0.4), (
        f"fitted rising-phase slope {rising:.3f} against Jahn's 2.2 +/- 0.4 over 1-10 uM")

    # and the price of matching a four-decade reconstruction of a one-decade measurement
    assert r.at_bound["log10_D"][0], (
        "log10_D used to sit on its lower bound here. If it no longer does, the "
        "reconstruction in fitting/data.py has stopped forcing desensitisation out and "
        "the artefact described in this docstring needs re-examining")
    assert float(y.max()) > 0.95, (
        f"the fit's absolute plateau is {y.max():.4f}; it used to be driven to ~0.995 by "
        f"the global-Hill reconstruction, which is the artefact this test documents")
    assert "ON a bound" in r.note


def test_a_normalised_dataset_is_compared_against_a_normalised_prediction(quiet):
    """The fix for the second observable mismatch (roadmap §2.3, data.Normalisation).

    The Jahn curve is I/I_max with a unit asymptote; this scheme's absolute peak open
    probability saturates near 0.75 at the three-anchor parameters. Compared directly, the
    residual at saturation is ~0.25 and the only way to close it is to delete
    desensitisation. Checked by scoring the same model against the same numbers tagged both
    ways and requiring the normalised one to fit far better.
    """
    import dataclasses

    three_anchor = KineticAllosteryModel(
        kon=0.0075785, koff=0.180415, beta=1.177468, alpha=0.30, d=0.05, r=0.002)
    as_absolute = dataclasses.replace(JAHN1997_PEAK_CRC, normalisation="absolute")
    norm = concentrated_log_likelihood(three_anchor, [JAHN1997_PEAK_CRC], n_params=3)
    absol = concentrated_log_likelihood(three_anchor, [as_absolute], n_params=3)
    assert norm.per_dataset[0].rss < 0.25 * absol.per_dataset[0].rss, (
        f"normalised RSS {norm.per_dataset[0].rss:.4f} against absolute "
        f"{absol.per_dataset[0].rss:.4f} -- the normalisation is not being applied")


def test_a_dataset_that_claims_normalisation_must_actually_be_normalised():
    """An I/I_max column reaches 1; one that does not would be divided at the wrong point."""
    import dataclasses

    from circuitpharm.fitting.data import DOSE_RESPONSE_BENCHMARK

    with pytest.raises(ValueError, match="largest response"):
        dataclasses.replace(DOSE_RESPONSE_BENCHMARK, normalisation="fraction_of_max")
    with pytest.raises(ValueError, match="unknown normalisation"):
        dataclasses.replace(DOSE_RESPONSE_BENCHMARK, normalisation="percent")


def test_the_slope_and_the_plateau_are_compatible_once_measured_the_same_way():
    """ANCHOR, and the CORRECTION of a result this project previously recorded.

    THE OLD CLAIM, now withdrawn: "the measured slope and the assumed plateau cannot both
    hold -- they are more than two sigma apart at every parameter set satisfying the
    other". That was an artefact of the measurement, not a property of the scheme.

    Jahn et al. report "a Hill-type slope of 2.2 +/- 0.4 OVER 1-10 uM" against an EC50 of
    11.6 uM -- a fit to the RISING PHASE, 0.086-0.862 x EC50, entirely below EC50. The old
    sweep instead regressed over the 2-98% band of the WHOLE curve. For a Hill curve those
    two agree exactly, which is why the distinction went unnoticed; for THIS scheme they do
    not, because its logit curve bends between a low-agonist limiting slope of 2 (the
    number of binding sites) and a flatter slope through EC50.

    The window was recorded in `fitting/data.py` beside the comparison all along. It was
    written down and not acted on, which is the more useful half of this lesson.

    Measured the way the source measured it, at the assumed plateau:

      * steepest rising-phase slope with P_o,max in [0.73, 0.77] is 1.69 on a 33 x 31
        grid, against 1.33 for the whole-curve regression on the same grid;
      * the gap to the published 2.2 +/- 0.4 falls from 2.2 sigma to 1.3 sigma, which is
        ordinary agreement, not a contradiction;
      * slope >= 1.8 at the plateau needs P_o,max >= ~0.91, not >= 0.99.

    So the scheme and the measurement are compatible, and nothing needs to be explained
    away. What survives is softer and still worth having: this scheme is a little flatter
    on the rising phase than the central published value, so a genuinely steeper curve
    would still be evidence against it -- see the companion assertion on the two-site
    maximum, which is what the slope CAN rule on.

    Grid coarsened to 11 x 9 as before, which costs some of the maximum: the fine grid's
    best point (E = 3.16) is not on the coarse one, so the figures below are 1.58 rather
    than 1.69 and 1.55 sigma rather than 1.28. The conclusion is unchanged, and the
    assertions are written against what the coarse grid actually produces rather than
    against the headline numbers, so this test cannot pass by quoting the docstring.
    """
    from circuitpharm.models.base import JAHN1997_SLOPE_WINDOW_REL, hill_slope

    koff, alpha, r_rate = 0.180415, 0.30, 0.002
    concs = np.logspace(-2.0, 4.5, 300)

    def measured(e_val, d_val):
        m = KineticAllosteryModel(kon=koff / 50.0, koff=koff, beta=alpha * e_val,
                                  alpha=alpha, d=r_rate * d_val, r=r_rate)
        y = np.asarray(m.peak_dose_response(concs), dtype=float)
        return (hill_slope(concs, y),                                   # whole curve
                hill_slope(concs, y, window_rel=JAHN1997_SLOPE_WINDOW_REL),  # as Jahn did
                float(y.max()))

    steepest_whole = steepest_rising = 0.0
    lowest_plateau_when_steep = 1.0
    for log_e in np.linspace(-1.0, 3.0, 11):
        for log_d in np.linspace(-3.0, 3.0, 9):
            whole, rising, p_max = measured(10.0 ** log_e, 10.0 ** log_d)
            if 0.73 <= p_max <= 0.77:
                steepest_whole = max(steepest_whole, whole)
                steepest_rising = max(steepest_rising, rising)
            if rising >= 1.8:
                lowest_plateau_when_steep = min(lowest_plateau_when_steep, p_max)

    # the two measurements of the SAME curves differ by about 0.35 -- this is the finding
    assert steepest_rising > steepest_whole + 0.25, (
        f"rising-phase {steepest_rising:.3f} vs whole-curve {steepest_whole:.3f}: the two "
        f"measurements have converged, so the correction recorded in "
        f"knowledge/12-inference.md Section 2 no longer has a cause and needs re-examining")

    # and measured consistently, the published value is NOT out of reach
    gap_sigma = (2.2 - steepest_rising) / 0.4
    assert gap_sigma < 2.0, (
        f"the published slope is {gap_sigma:.2f} sigma above the steepest rising-phase "
        f"slope at the plateau; the withdrawn claim was that this exceeds 2 sigma")
    assert steepest_rising == pytest.approx(1.58, abs=0.08)   # 1.69 on the fine grid

    # slope >= 1.8 no longer demands a plateau of 0.99
    assert lowest_plateau_when_steep < 0.95, (
        f"slope >= 1.8 still needs P_o,max {lowest_plateau_when_steep:.4f}")


def test_two_binding_sites_can_produce_the_published_slope():
    """What the slope CAN rule on, which is not what Jahn et al. used it for.

    Jahn et al. read 2.2 +/- 0.4 as evidence for at least THREE binding sites. It is not:
    a two-site scheme has a limiting log-log slope of exactly 2 at low agonist, and a Hill
    fit over the rising phase samples that limit. Scanned over E and D, the steepest
    rising-phase slope this TWO-site scheme can produce is ~2.0 -- within half a sigma of
    the published value.

    So the published slope does not discriminate two sites from three, and the inference
    drawn from it in the source does not follow. This is a statement about the measurement,
    not a correction to the receptor: three sites may well be right for other reasons.
    """
    from circuitpharm.models.base import JAHN1997_SLOPE_WINDOW_REL, hill_slope

    koff, alpha, r_rate = 0.180415, 0.30, 0.002
    concs = np.logspace(-2.0, 4.5, 300)
    best = 0.0
    for log_e in np.linspace(-1.0, 3.0, 11):
        for log_d in np.linspace(-3.0, 3.0, 9):
            m = KineticAllosteryModel(kon=koff / 50.0, koff=koff, beta=alpha * 10.0 ** log_e,
                                      alpha=alpha, d=r_rate * 10.0 ** log_d, r=r_rate)
            y = np.asarray(m.peak_dose_response(concs), dtype=float)
            best = max(best, hill_slope(concs, y, window_rel=JAHN1997_SLOPE_WINDOW_REL))

    assert best == pytest.approx(2.0, abs=0.15)
    assert (2.2 - best) / 0.4 < 1.0, (
        f"a two-site scheme reaches a rising-phase slope of {best:.3f}; if this ever falls "
        f"far below the published 2.2 - 1 sigma, the slope WOULD discriminate site count "
        f"and knowledge/12-inference.md Section 2 needs rewriting the other way")


def test_the_peak_observable_is_protocol_dependent():
    """PEAK is not one observable unless the application duration is stated.

    The peak of a DESENSITISING response depends on how long the agonist is applied: at low
    agonist the response is still rising when a short application ends, so the measured
    peak and the concentration at which it half-saturates both move. Measured on this
    scheme at the plateau-respecting parameters, varying only the application length:

        2 ms   -> EC50 ~300 uM      1000 ms -> EC50 ~39 uM

    a SEVENFOLD shift from the protocol alone, with the rising-phase slope moving 1.69-1.79
    over the same range. `PEAK_APPLICATION_MS` is 300 ms here; Jahn et al. used ultra-fast
    exchange and the abstract does not state the application length, so their EC50 of
    11.6 uM and this scheme's are not strictly comparable either.

    This is the third observable-mismatch class in the project, after PEAK-vs-EQUILIBRIUM
    (P1) and absolute-vs-normalised (P4): same tag, same units, different protocol. It is
    pinned because a future change to PEAK_APPLICATION_MS would silently move every peak
    EC50 the project reports.
    """
    from circuitpharm.models.base import JAHN1997_SLOPE_WINDOW_REL, hill_slope

    koff, alpha, r_rate = 0.180415, 0.30, 0.002
    m = KineticAllosteryModel(kon=koff / 50.0, koff=koff, beta=alpha * 3.16, alpha=alpha,
                              d=r_rate * 1e-3, r=r_rate)
    concs = np.logspace(-2.0, 4.0, 140)

    def ec50_for(ms):
        y = np.asarray(m.peak_dose_response(concs, application_ms=ms), dtype=float)
        return float(np.interp(0.5, y / y.max(), concs)), hill_slope(
            concs, y, window_rel=JAHN1997_SLOPE_WINDOW_REL)

    ec50_short, slope_short = ec50_for(2.0)
    ec50_long, slope_long = ec50_for(1000.0)

    assert ec50_short / ec50_long > 5.0, (
        f"peak EC50 moved only {ec50_short / ec50_long:.1f}x between a 2 ms and a 1000 ms "
        f"application ({ec50_short:.1f} -> {ec50_long:.1f} uM); the protocol dependence "
        f"this test exists to pin has gone")
    assert slope_short == pytest.approx(1.77, abs=0.10)
    assert slope_long == pytest.approx(1.69, abs=0.10)


def _peak_ec50_and_slope(model) -> tuple[float, float]:
    """EC50 and Hill slope of a PEAK curve, by regression in logit space.

    NOT by interpolating a decade-spaced grid, which is how a previous version of this
    measurement reported 12.08 uM for a curve whose EC50 was 11.6 (roadmap §5.3 E9).
    """
    c = np.logspace(-1.0, 4.0, 400)
    y = np.asarray(model.peak_dose_response(c), dtype=float)
    y = y / y.max()
    sel = (y > 0.02) & (y < 0.98)
    a = np.vstack([np.log10(c[sel]), np.ones(int(sel.sum()))]).T
    b = np.log(y[sel] / (1.0 - y[sel]))
    coef, *_ = np.linalg.lstsq(a, b, rcond=None)
    return float(10.0 ** (-coef[1] / coef[0])), float(coef[0] / np.log(10.0))


# ============================================================= the sampler (P4-4)
@pytest.mark.slow
def test_sampler_recovers_a_known_gaussian():
    """ANCHOR, and the sampler's own method check, independent of any receptor.

    A 3-D correlated Gaussian has a closed-form answer, so this separates "the sampler
    works" from "the posterior is well determined" -- which the receptor problem below
    cannot do, since a failure there could be either. Tolerances are set from the Monte
    Carlo error at this chain length (~2500 effective samples), not from what passed.
    """
    mean = np.array([1.0, -0.5, 2.0])
    cov = np.array([[1.0, 0.6, 0.2], [0.6, 1.0, 0.4], [0.2, 0.4, 1.0]])
    prec = np.linalg.inv(cov)
    lo, hi = mean - 12.0, mean + 12.0

    def log_prob(v):
        if np.any(v <= lo) or np.any(v >= hi):
            return -np.inf
        d = v - mean
        return float(-0.5 * d @ prec @ d)

    run = sample_ensemble(log_prob, lo=lo, hi=hi, names=("a", "b", "c"),
                          n_walkers=30, n_steps=12000, burn_in=3000, seed=5,
                          a=(2.0, 20.0))
    assert run.converged, run.diagnostics["failures"]
    # 9000 kept steps at tau ~ 44 is ~200 effective samples per walker. 4000 steps was
    # NOT enough -- split-Rhat came out 1.0167 against a 1.01 gate -- and the gate is
    # right: lengthening the chain is the fix, loosening it would be the bug.
    assert run.samples.mean(axis=0) == pytest.approx(mean, abs=0.1)
    assert np.cov(run.samples.T) == pytest.approx(cov, abs=0.1)
    assert not run.diagnostics["frozen_walkers"]


def test_a_frozen_walker_is_a_diagnostic_failure():
    """A bimodal target with a deep valley freezes walkers at a=2; the run must say so.

    This is the trap that cost a day of P4: 2 of 24 walkers sat 2.5 decades from the mode
    for an entire 6000-step run while the other 22 reported the right answer, and the only
    sign was a split-Rhat of 62. The stretch move cannot propose a contraction below 1/a,
    so a walker across a valley cannot reach the ensemble. Detected now rather than noticed.
    """
    lo, hi = np.array([-12.0, -12.0]), np.array([12.0, 12.0])

    def log_prob(v):
        if np.any(v <= lo) or np.any(v >= hi):
            return -np.inf
        near = -0.5 * float((v - 0.0) @ (v - 0.0)) / 0.04
        far = -0.5 * float((v - 10.0) @ (v - 10.0)) / 0.04 - 60.0
        return float(max(near, far))

    run = sample_ensemble(log_prob, lo=lo, hi=hi, names=("x", "y"), n_walkers=20,
                          n_steps=600, burn_in=100, seed=2, a=2.0, init="prior")
    assert not run.converged
    assert any("frozen" in f for f in run.diagnostics["failures"]), (
        run.diagnostics["failures"])
    assert run.diagnostics["frozen_walkers"]


def test_mixing_stretch_scales_is_what_fixes_it():
    """The same target, the same seed, with a large scale added: no frozen walker.

    Pins the claim in `sample_ensemble`'s docstring. If this ever fails while the test
    above passes, the mixture stopped being applied.
    """
    lo, hi = np.array([-12.0, -12.0]), np.array([12.0, 12.0])

    def log_prob(v):
        if np.any(v <= lo) or np.any(v >= hi):
            return -np.inf
        near = -0.5 * float((v - 0.0) @ (v - 0.0)) / 0.04
        far = -0.5 * float((v - 10.0) @ (v - 10.0)) / 0.04 - 60.0
        return float(max(near, far))

    run = sample_ensemble(log_prob, lo=lo, hi=hi, names=("x", "y"), n_walkers=20,
                          n_steps=600, burn_in=100, seed=2, a=(2.0, 20.0, 200.0),
                          init="prior")
    assert not run.diagnostics["frozen_walkers"], run.diagnostics["accepts_per_walker"]


def test_a_single_stretch_scale_below_one_is_refused():
    lo, hi = np.zeros(2), np.ones(2) * 3.0
    with pytest.raises(ValueError, match="stretch scale"):
        sample_ensemble(lambda v: 0.0, lo=lo, hi=hi, names=("x", "y"), a=1.0,
                        n_walkers=8, n_steps=20, burn_in=5)
    with pytest.raises(ValueError, match="stretch scale"):
        sample_ensemble(lambda v: 0.0, lo=lo, hi=hi, names=("x", "y"), a=(2.0, 0.5),
                        n_walkers=8, n_steps=20, burn_in=5)


# ======================================================= the posterior (P4-4), slow
@pytest.mark.slow
def test_the_convergence_gate_is_met_with_margin_not_on_the_nose(posterior):
    """The gate must be met by more than the diagnostic varies BETWEEN PLATFORMS.

    E23. At the former default of 12000 steps / 4000 burn-in this posterior reached split-Rhat
    = 1.01089 against `max_rhat = 1.01`, i.e. it failed by 0.0009. The chain is seeded, so it
    is deterministic on any one machine -- and that is what made it nasty. x86/OpenBLAS and
    arm64/Accelerate sum in different orders, that propagates through the matrix exponentials
    and 24 walkers into the third decimal of Rhat, and the result was CI permanently GREEN
    while a macOS checkout was permanently RED on identical code. A test that passes on the
    build machine and fails on every developer's is worse than one that fails everywhere.

    `converged` decides whether percentiles are reported as credible intervals or marked VOID,
    so the honest fix was a longer chain (24000/8000, measured Rhat 1.00377), not a looser
    threshold -- 1.01 is the Vehtari et al. recommendation and worth keeping.

    This test pins the MARGIN rather than the pass, because a pass alone cannot distinguish
    "converged" from "converged by 0.0001 on this machine only". If a future scheme or sampler
    change erodes the margin, this fails while `converged` is still True, which is the early
    warning the previous arrangement did not have.
    """
    rhat = (posterior.diagnostics.get("split_rhat")
            or posterior.diagnostics.get("rhat") or {})
    assert rhat, f"no split-Rhat in diagnostics: {sorted(posterior.diagnostics)}"
    worst = max(rhat.values())

    #: Half the distance between the achieved 1.00377 and the 1.01 gate. Any platform
    #: difference large enough to cross this is large enough to want investigating.
    MARGIN_FLOOR = 1.006
    assert worst <= MARGIN_FLOOR, (
        f"worst split-Rhat is {worst:.5f}, inside {MARGIN_FLOOR} of the 1.01 gate. The chain "
        f"converges on this machine but the margin is now thinner than the arithmetic "
        f"difference between platforms, so CI and a local checkout can disagree. Lengthen the "
        f"chain rather than moving the gate.")


def test_posterior_recovers_known_parameters(posterior):
    """ANCHOR. Every 95% interval must cover the truth it was generated from.

    Three intervals at 95% is a weak coverage test on its own -- that is why the Gaussian
    check above exists -- but a posterior that misses a truth it was handed is
    disqualifying, and this is the cheapest way to notice.
    """
    assert posterior.converged, posterior.diagnostics["failures"]
    assert posterior.tier is Tier.UNCALIBRATED, (
        "the data is generated, so the interval cannot be VALIDATED however well the "
        "chain converged")
    for i, name in enumerate(EQUILIBRIUM_IDENTIFIABLE):
        lo, hi = posterior.intervals_95[name].value
        assert lo <= TRUTH_VEC[i] <= hi, (
            f"{name}: 95% interval [{lo:.4f}, {hi:.4f}] excludes the truth "
            f"{TRUTH_VEC[i]:.4f}")
        assert not any(posterior.pinned_at_bound[name]), (
            f"{name}'s interval reaches its prior bound, so it is reporting the box")


@pytest.mark.slow
def test_posterior_median_in_profile_ci(posterior, dataset, quiet):
    """ANCHOR. Two independent uncertainty machineries must agree.

    Profile likelihood re-optimises the other parameters at each node of a grid; MCMC
    integrates over them. They share no optimiser, no grid and no convergence criterion, so
    agreement is the strongest check available here and disagreement means one is wrong.
    Both the medians and the intervals are checked: a median inside a very wide profile
    interval would pass while the two still disagreed about the uncertainty.
    """
    report = agreement_with_profiles(posterior, dataset)
    for name, r in report.items():
        assert r.median_inside_profile_ci, (
            f"{name}: posterior median {r.posterior_median:.4f} is outside the profile "
            f"95% CI {r.profile_ci_95}")
        assert r.interval_overlap > 0.8, (
            f"{name}: the two 95% intervals overlap by only "
            f"{r.interval_overlap:.2f} of the narrower one -- posterior "
            f"{r.posterior_ci_95} against profile {r.profile_ci_95}")


@pytest.mark.slow
def test_mcmc_diagnostics_gate_the_result(dataset, quiet):
    """ANCHOR. An under-run chain returns VOID, not a narrow interval.

    10 steps cannot estimate an autocorrelation time, let alone a credible interval. The
    percentiles still exist -- they are where the walkers happened to be -- and the point
    is that reading one RAISES rather than returning a plausible number.
    """
    post = sample_identifiable_posterior(factory, [dataset], n_steps=10, burn_in=2,
                                         seed=4)
    assert not post.converged
    assert post.tier is Tier.VOID
    for name in EQUILIBRIUM_IDENTIFIABLE:
        with pytest.raises(VoidQuantityError):
            post.intervals_95[name].value
        with pytest.raises(VoidQuantityError):
            post.medians[name].value
        # the number is still reachable for debugging, at a call site a reviewer can see
        assert len(post.intervals_95[name].get(acknowledge_void=True)) == 2


@pytest.mark.slow
def test_posterior_report_never_prints_a_void_value(dataset, quiet):
    post = sample_identifiable_posterior(factory, [dataset], n_steps=10, burn_in=2,
                                         seed=4)
    text = str(posterior_report(post))
    assert "VOID" in text
    for name in EQUILIBRIUM_IDENTIFIABLE:
        raw = post.medians[name].get(acknowledge_void=True)
        assert f"{raw:.4f}" not in text, (
            f"the report printed {name}'s VOID value")


# =================================================== the parameter registry (P4-3)
def test_registry_refuses_to_invent_a_fit_for_model_c():
    """P5 cannot rank a fitted model against an unfitted one and call it evidence."""
    from circuitpharm import parameters

    with pytest.raises(NotImplementedError, match="never been fitted"):
        parameters.get("extended_desens", "fitted")
    with pytest.raises(KeyError, match="unknown model"):
        parameters.get("not_a_model")
    with pytest.raises(ValueError, match="fitted, declared or nominal"):
        parameters.get("kinetic_jw95", "sensible")


def test_registry_nominal_parameters_are_void():
    """The dataclass defaults are reachable and unquotable, which is the P0-13 fix."""
    from circuitpharm import parameters

    for model in ("kinetic_jw95", "operational", "extended_desens"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rates, q = parameters.get(model, "nominal")
        assert q.tier is Tier.VOID
        with pytest.raises(VoidQuantityError):
            q.value
        assert rates, model
        assert q.get(acknowledge_void=True) == rates


def test_registry_fitted_kinetic_carries_its_caveats():
    """The fitted set must arrive with d/r, FIT_FIXED_ALPHA and the Jahn EC50 named."""
    from circuitpharm import parameters

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rates, q = parameters.get("kinetic_jw95", "fitted")
    assert q.tier is Tier.UNCALIBRATED
    joined = " ".join(q.caveats).lower()
    for needed in ("d and r", "convention", "11.6"):
        assert needed.lower() in joined, f"the caveat about {needed!r} is missing"
    assert rates["koff"] / rates["kon"] == pytest.approx(23.8062, rel=1e-3)
    assert rates["alpha"] == pytest.approx(0.30, rel=1e-12)


def test_registry_build_returns_a_usable_model():
    from circuitpharm import parameters

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model, q = parameters.build("kinetic_jw95", "fitted")
    assert isinstance(model, KineticAllosteryModel)
    assert q.tier is Tier.UNCALIBRATED
    assert model.peak_dose_response(np.array([1000.0]))[0] > 0.5

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model_a, _ = parameters.build("operational", "declared")
    assert isinstance(model_a, OperationalScalarModel)


def test_no_model_is_instantiated_from_its_defaults_in_production_code():
    """ANCHOR. Every model in `src/` must get its parameters from the registry.

    AST-walked rather than grepped, so a call split across lines or wrapped in another
    expression cannot hide. `parameters.py` is exempt because it IS the registry, and
    `models/` is exempt because a class may construct its own kind (`apply_pam` returns a
    modified copy). Everything else that writes `KineticAllosteryModel()` is reintroducing
    the defect P0-13 records: the project's headline numbers were produced at a chimeric
    parameter set that reproduced no anchor.
    """
    import ast
    import pathlib

    targets = {"KineticAllosteryModel", "ExtendedDesensitizationModel",
               "OperationalScalarModel"}
    root = pathlib.Path(__file__).resolve().parent.parent / "src" / "circuitpharm"
    exempt = {"parameters.py"}
    offenders = []
    for path in sorted(root.rglob("*.py")):
        if path.name in exempt or path.parent.name == "models":
            continue
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            fn = node.func
            name = (fn.id if isinstance(fn, ast.Name)
                    else fn.attr if isinstance(fn, ast.Attribute) else None)
            if name in targets and not node.args and not node.keywords:
                offenders.append(f"{path.relative_to(root)}:{node.lineno} {name}()")
    assert not offenders, (
        "these call sites build a model from its dataclass defaults instead of asking "
        "`parameters.get`: " + "; ".join(offenders))


def test_registry_report_names_what_has_no_fit():
    from circuitpharm import parameters

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        text = parameters.report()
    assert "extended_desens" in text
    assert "NO FIT EXISTS FOR" in text
    assert "extended_desens" in text.split("NO FIT EXISTS FOR")[1]


@pytest.mark.slow
def test_a_saved_chain_names_its_inputs(dataset, quiet, tmp_path):
    """A chain archive must carry the sha, the datasets and the dependency versions.

    Roadmap P4-2. Six months on, "the posterior" in a figure caption has to resolve to a
    file that says what produced it; a bare array of numbers cannot.
    """
    import json

    from circuitpharm.fitting.posterior import save_chain

    post = sample_identifiable_posterior(factory, [dataset], n_steps=400, burn_in=100,
                                         seed=7)
    path = save_chain(post, tmp_path / "chain.npz", [dataset], label="a method check")
    with np.load(path, allow_pickle=False) as z:
        assert z["log10_samples"].shape == (post.log10_samples.shape[0], 3)
        meta = json.loads(str(z["metadata"]))
    assert meta["git_sha"]
    assert meta["label"] == "a method check"
    assert meta["param_names"] == list(EQUILIBRIUM_IDENTIFIABLE)
    assert meta["datasets"][0]["label"] == dataset.citation_label
    assert meta["datasets"][0]["synthetic"] is True
    assert meta["datasets"][0]["normalisation"] == "absolute"
    assert meta["versions"]["scipy"]
    # an under-run chain is VOID and the archive must say so rather than looking clean
    assert meta["tier"] == "VOID"
    assert meta["diagnostics"]["failures"]
