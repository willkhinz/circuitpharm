"""P5: model comparison that means something.

The comparison this replaces (`protocols/oed.py::evaluate_model_fit`) ranked models at
unfitted parameters, with `k` counting declared parameters and `sigma` hand-set. Those are
three defects, not three approximations, and each has a test here that it has been fixed.

The ANCHOR is `test_nested_model_recovers_the_truth`: generate from Model B, and the
protocol must select Model B over Model A and must NOT select Model C, which nests B. A
protocol that cannot recover a truth it was handed cannot adjudicate one it was not.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from circuitpharm.fitting.comparison import (
    FIT_CONVENTIONS,
    MODEL_SPECS,
    blocked_cv,
    compare_models,
    comparison_report,
    fit_one,
    recovery_check,
    resolve_cv,
)
from circuitpharm.fitting.data import equilibrium_crc_from_model, peak_crc_from_model
from circuitpharm.models.base import ObservableMismatch
from circuitpharm.results import Tier, VoidQuantityError

ALL_THREE = ("operational", "kinetic_jw95", "extended_desens")


@pytest.fixture(scope="module")
def quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


@pytest.fixture(scope="module")
def peak_from_b():
    """PEAK data generated from Model B at the midpoint of its own fitting box.

    The midpoint, so the truth is inside every model's search space and a failure to
    recover it cannot be blamed on the bounds. PEAK, because Model A has no equilibrium and
    refuses to produce one, so it is the only axis all three can be compared on.
    """
    spec = MODEL_SPECS["kinetic_jw95"]
    centre = np.array([0.5 * (lo + hi) for lo, hi in spec.bounds])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return peak_crc_from_model(
            spec.factory(centre), concs_um=np.logspace(-1.0, 3.5, 24), noise_sd=0.005,
            seed=0, label="P5 test fixture, truth = Model B")


# ======================================================== the three fixed defects
def test_k_counts_identifiable_params_only(quiet, peak_from_b):
    """ANCHOR. `k` is estimated identifiable parameters plus estimated noise scales.

    The old count was `len(m.param_names)`: 6 for Model A, 6 for Model B, 9 for Model C.
    That charges Model C for a `pam_desens_factor` the pam_factor = 1 curve never touches
    and charges Model B for `d` and `r` separately when only `d/r` enters the likelihood --
    a penalty wrong in the direction that systematically disfavours the detailed scheme.
    """
    from circuitpharm.models import (
        ExtendedDesensitizationModel,
        KineticAllosteryModel,
        OperationalScalarModel,
    )

    declared = {"operational": len(OperationalScalarModel().param_names),
                "kinetic_jw95": len(KineticAllosteryModel().param_names),
                "extended_desens": len(ExtendedDesensitizationModel().param_names)}
    for name in ALL_THREE:
        fit = fit_one(MODEL_SPECS[name], [peak_from_b], n_starts=3)
        # one scale per dataset, profiled out of the likelihood but still estimated
        assert fit.n_sigma == 1
        assert fit.k == fit.n_free + 1
        assert fit.k < declared[name], (
            f"{name}: k = {fit.k} is not below the declared count {declared[name]}, so "
            f"the fix has been undone")
    # and Model C is charged exactly one more than Model B, which it nests
    b = fit_one(MODEL_SPECS["kinetic_jw95"], [peak_from_b], n_starts=3)
    c = fit_one(MODEL_SPECS["extended_desens"], [peak_from_b], n_starts=3)
    assert c.k == b.k + 1


def test_the_ranking_does_not_depend_on_an_assumed_sigma(quiet, peak_from_b):
    """ANCHOR. Multiplying the declared `sem` by 100 must change nothing at all.

    The old likelihood took `measurement_noise_std=0.03` and every AIC difference scaled
    with it: at 0.03 the Akaike weights were near-uniform, at 0.003 a delta. A ranking that
    moves that far on an unjustified knob is not a ranking. Here the scale is PROFILED OUT
    -- estimated from the residuals -- so the declared column is not an input at all, and
    the fits, the criteria and the verdict must be bit-for-bit identical.

    A WEAKER VERSION OF THIS TEST WAS WRONG AND IS RECORDED HERE SO IT IS NOT REWRITTEN.
    It rescaled the DATA by 0.5 and asserted the fitted `log10_kd` barely moved. It moved
    from 1.512 to 1.236, and correctly: halving the responses asks for a curve of half the
    amplitude, which this scheme reaches by changing `E`, which moves the EC50 with it.
    That is the model working, not the likelihood failing. The claim being tested is
    invariance to the DECLARED uncertainty, not to the data.
    """
    import dataclasses

    loud = dataclasses.replace(peak_from_b, sem=peak_from_b.sem * 100.0)
    for name in ("kinetic_jw95", "extended_desens"):
        a = fit_one(MODEL_SPECS[name], [peak_from_b], n_starts=3)
        b = fit_one(MODEL_SPECS[name], [loud], n_starts=3)
        assert a.log_likelihood == pytest.approx(b.log_likelihood, rel=1e-12), name
        assert a.aicc == pytest.approx(b.aicc, rel=1e-12), name
        for k in a.values:
            assert a.values[k] == pytest.approx(b.values[k], rel=1e-12), f"{name}.{k}"


def test_every_model_is_fitted_before_it_is_ranked(quiet, peak_from_b):
    """An information criterion is a statement about a model AT its maximum.

    Checked by requiring each model's fitted log-likelihood to beat its own value at the
    centre of its box. If a model were being scored at handed-in parameters this would
    fail, which is precisely what the old `evaluate_model_fit` did.
    """
    from circuitpharm.fitting.likelihood import make_log_likelihood

    for name in ALL_THREE:
        spec = MODEL_SPECS[name]
        fit = fit_one(spec, [peak_from_b], n_starts=3)
        centre = np.array([0.5 * (lo + hi) for lo, hi in spec.bounds])
        at_centre = make_log_likelihood(
            [peak_from_b], caller="test", n_params=spec.n_free
        )(spec.factory(centre)).total_log_likelihood
        assert fit.log_likelihood > at_centre, name


# ================================================================ blocked CV
def test_cv_folds_are_blocked_not_pointwise(quiet, peak_from_b):
    """Folds must be contiguous in concentration, and must partition the foldable points.

    Random point-wise folds leak across a smooth curve -- a point's neighbours are almost
    its answer -- so they score interpolation and flatter the most flexible model.
    """
    from circuitpharm.fitting.comparison import _blocks

    concs = np.asarray(peak_from_b.concs_um, dtype=float)
    folds = _blocks(peak_from_b, 4)
    assert len(folds) == 4
    # a partition: every foldable index exactly once
    flat = np.concatenate(folds)
    assert flat.size == len(set(flat.tolist()))
    # contiguous in sorted concentration: each block's range must not straddle another's
    spans = [(concs[b].min(), concs[b].max()) for b in folds]
    for (lo1, hi1), (lo2, hi2) in zip(spans, spans[1:]):
        assert hi1 <= lo2, f"blocks {spans} overlap in concentration, so they are not blocks"


def test_the_normalising_concentration_is_never_held_out():
    """For an I/I_max dataset the top point DEFINES the scale, so it always trains.

    A fold that held it out would score the model against a normalisation nobody gave it,
    and the prediction would be normalised at a different concentration than the data.
    """
    from circuitpharm.fitting.comparison import _blocks

    spec = MODEL_SPECS["kinetic_jw95"]
    centre = np.array([0.5 * (lo + hi) for lo, hi in spec.bounds])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        norm = peak_crc_from_model(spec.factory(centre),
                                   concs_um=np.logspace(-1.0, 3.5, 12),
                                   normalisation="fraction_of_max")
        absolute = peak_crc_from_model(spec.factory(centre),
                                       concs_um=np.logspace(-1.0, 3.5, 12))
    top = int(np.argmax(np.asarray(norm.concs_um)))
    assert top not in np.concatenate(_blocks(norm, 3)).tolist()
    assert top in np.concatenate(_blocks(absolute, 3)).tolist()


def test_cv_scores_the_heldout_block_at_the_training_sigma(quiet, peak_from_b):
    """Re-estimating sigma on the held-out block is grading your own exam.

    A model that predicts the block badly could then declare the noise large and score as
    well as one that predicts it well. Checked by requiring Model A -- which cannot
    reproduce the 5-state peak shape and therefore extrapolates badly -- to score far
    worse out of sample than in sample.
    """
    a_fit = fit_one(MODEL_SPECS["operational"], [peak_from_b], n_starts=3)
    a_cv = blocked_cv(MODEL_SPECS["operational"], [peak_from_b], n_folds=4, n_starts=2)
    assert a_cv.total_lpd < a_fit.log_likelihood, (
        "Model A's held-out score is not worse than its in-sample fit, which means the "
        "held-out points are being scored at a scale fitted to themselves")
    assert a_cv.n_failed_folds == 0


# ============================================================= the tie-break
def test_a_sub_two_nat_cv_margin_is_a_tie_and_parsimony_wins():
    """The measured case: Model C beat Model B out of sample by 0.84 nats, on nothing.

    Per-fold differences were [0.0, 0.7, 0.0, 0.1]. Calling that a win selects the larger
    model over the one the data came from. `resolve_cv` requires the margin to clear both
    the paired standard error and a 2-nat floor, and breaks a tie on fewest estimated
    parameters -- the one-standard-error rule.
    """
    from circuitpharm.fitting.comparison import CVScore, ModelFit

    def fit(name, k):
        return ModelFit(name=name, values={}, log_likelihood=0.0, n_free=k - 1, n_sigma=1,
                        n_data=24, aic=0.0, aicc=0.0, bic=0.0, at_bound=(),
                        converged_uniquely=True, note="")

    b_folds = (27.2, 18.8, 17.7, 23.1)
    c_folds = (27.2, 19.5, 17.7, 23.2)
    cv = {"kinetic_jw95": CVScore("kinetic_jw95", float(sum(b_folds)), b_folds, 4, 0),
          "extended_desens": CVScore("extended_desens", float(sum(c_folds)), c_folds,
                                     4, 0)}
    fits = {"kinetic_jw95": fit("kinetic_jw95", 4),
            "extended_desens": fit("extended_desens", 5)}
    v = resolve_cv(cv, fits)
    assert set(v.tied) == {"kinetic_jw95", "extended_desens"}
    assert v.winner == "kinetic_jw95", "parsimony must break the tie"
    assert v.margin_nats == pytest.approx(sum(c_folds) - sum(b_folds), abs=1e-9)
    assert 0.5 < v.margin_nats < 1.0
    assert "fewest estimated parameters" in v.note


def test_a_large_cv_margin_is_not_a_tie():
    """The rule must still be able to name a winner, or it is not a rule."""
    from circuitpharm.fitting.comparison import CVScore, ModelFit

    def fit(name, k):
        return ModelFit(name=name, values={}, log_likelihood=0.0, n_free=k - 1, n_sigma=1,
                        n_data=24, aic=0.0, aicc=0.0, bic=0.0, at_bound=(),
                        converged_uniquely=True, note="")

    cv = {"big": CVScore("big", 80.0, (20.0, 20.0, 20.0, 20.0), 4, 0),
          "small": CVScore("small", 40.0, (10.0, 10.0, 10.0, 10.0), 4, 0)}
    fits = {"big": fit("big", 6), "small": fit("small", 3)}
    v = resolve_cv(cv, fits)
    assert v.winner == "big"
    assert v.tied == ("big",)
    assert v.margin_nats == pytest.approx(40.0)


# ============================================================= the observable contract
def test_comparison_requires_same_observable(quiet):
    """Model A has no equilibrium, so it cannot be ranked on one (P1 contract)."""
    spec = MODEL_SPECS["kinetic_jw95"]
    centre = np.array([0.5 * (lo + hi) for lo, hi in spec.bounds])
    eq = equilibrium_crc_from_model(spec.factory(centre),
                                    concs_um=np.logspace(-1.0, 3.0, 12), noise_sd=0.005)
    with pytest.raises(ObservableMismatch, match="cannot produce it|no equilibrium curve"):
        compare_models(["operational", "kinetic_jw95"], [eq], n_folds=3, n_starts=2)
    # and B against C on equilibrium is fine, since both have one
    res = compare_models(["kinetic_jw95", "extended_desens"], [eq], n_folds=3,
                         n_starts=2, cv_starts=2)
    assert res.cv_winner in ("kinetic_jw95", "extended_desens")


def test_datasets_spanning_two_observables_are_refused(quiet, peak_from_b):
    spec = MODEL_SPECS["kinetic_jw95"]
    centre = np.array([0.5 * (lo + hi) for lo, hi in spec.bounds])
    eq = equilibrium_crc_from_model(spec.factory(centre), noise_sd=0.005)
    with pytest.raises(ObservableMismatch, match="span observables"):
        compare_models(["kinetic_jw95", "extended_desens"], [peak_from_b, eq])


def test_a_cross_native_comparison_says_so(quiet, peak_from_b):
    """B and C are EQUILIBRIUM-native; a PEAK comparison must carry that caveat."""
    res = compare_models(ALL_THREE, [peak_from_b], n_folds=4, n_starts=3, cv_starts=2)
    joined = " ".join(res.verdict.caveats)
    assert "native observable" in joined
    assert "EQUILIBRIUM-native" in joined


# ============================================================= the holdout
def test_the_holdout_cannot_be_trained_on(quiet, peak_from_b):
    import dataclasses

    held = dataclasses.replace(peak_from_b, role="holdout")
    with pytest.raises(ValueError, match="holdout"):
        compare_models(["kinetic_jw95", "extended_desens"], [held])


def test_the_holdout_is_scored_once_and_reported_separately(quiet, peak_from_b):
    import dataclasses

    spec = MODEL_SPECS["kinetic_jw95"]
    centre = np.array([0.5 * (lo + hi) for lo, hi in spec.bounds])
    held = dataclasses.replace(
        peak_crc_from_model(spec.factory(centre),
                            concs_um=np.logspace(-0.5, 3.0, 10), noise_sd=0.005, seed=9),
        role="holdout")
    res = compare_models(["kinetic_jw95", "extended_desens"], [peak_from_b], [held],
                         n_folds=4, n_starts=3, cv_starts=2)
    assert set(res.holdout_log_likelihood) == {"kinetic_jw95", "extended_desens"}
    # the holdout must not be what chose the winner: the verdict's provenance cites CV
    assert "holdout" not in res.verdict.provenance.lower()


# ============================================================= the verdict
def test_disagreeing_rankings_produce_a_void_verdict():
    """When CV and AICc disagree the comparison is unresolved, and saying so is the result.

    Constructed directly rather than hunting for a dataset that produces disagreement: the
    behaviour under disagreement is what matters and it must not depend on finding an
    example.
    """
    from circuitpharm.fitting.comparison import CVScore, ModelFit, _caveats

    # exercised through the same path compare_models uses, via resolve_cv's output
    cv = {"x": CVScore("x", 80.0, (20.0,) * 4, 4, 0),
          "y": CVScore("y", 40.0, (10.0,) * 4, 4, 0)}
    fits = {"x": ModelFit("x", {}, 0.0, 3, 1, 24, 0.0, 10.0, 0.0, (), True, ""),
            "y": ModelFit("y", {}, 0.0, 3, 1, 24, 0.0, 1.0, 0.0, (), True, "")}
    v = resolve_cv(cv, fits)
    assert v.winner == "x"
    assert min(fits, key=lambda n: fits[n].aicc) == "y"
    assert _caveats(fits, [], cv) != ()


@pytest.mark.slow
def test_the_verdict_is_void_when_the_two_rankings_disagree(quiet, peak_from_b):
    """End to end: a disagreement must make the Akaike weights unreadable, not just noted."""
    # min_nats=0 disables the parsimony tie-break, which is what produced the measured
    # disagreement in the first place: CV then prefers Model C by 0.84 nats while AICc
    # prefers Model B. Forced here on purpose -- the point is the VOIDing, not the gap.
    res = compare_models(["kinetic_jw95", "extended_desens"], [peak_from_b], n_folds=4,
                         n_starts=4, cv_starts=2, min_nats=0.0)
    # NOT A CONDITIONAL SKIP. The disagreement is deterministic at this seed and this is
    # the measured case -- CV prefers Model C by 0.84 nats while AICc prefers Model B. If it
    # ever stops happening that is a finding about the comparison, not a reason to pass
    # quietly; and the FULL-install CI job fails on any skip, so a skip here would read as a
    # missing optional dependency.
    assert not res.rankings_agree, (
        f"with the parsimony floor disabled the two rankings now agree on "
        f"{res.cv_winner}. That is the measured disagreement gone, and the VOID path below "
        f"is then untested -- construct a disagreement explicitly rather than skipping")
    assert res.verdict.tier is Tier.VOID
    with pytest.raises(VoidQuantityError):
        res.verdict.value
    for name in res.akaike_weights:
        with pytest.raises(VoidQuantityError):
            res.akaike_weights[name].value
    assert "UNRESOLVED" in res.note


@pytest.mark.slow
def test_nested_model_recovers_the_truth(quiet):
    """ANCHOR. Generate from Model B; select Model B.

    Over Model A, which has the wrong shape, and NOT Model C, which nests B and therefore
    can never fit worse -- so the only thing that stops it winning is the complexity
    penalty and the tie-break actually working. Measured: AICc -194.18 for B against
    -190.98 for C, and a CV difference of 0.84 nats that `resolve_cv` correctly calls a
    tie.

    If this fails, nothing the comparison says about Models A, B and C on real data means
    anything, and that is the whole reason the test exists.
    """
    r = recovery_check("kinetic_jw95", noise_sd=0.005, n_points=24, n_starts=4,
                       n_folds=4, seed=0)
    assert r.recovered, r.note
    assert r.selected == "kinetic_jw95"
    assert r.rankings_agree
    # and Model A must be decisively worse, not merely second
    assert r.cv["kinetic_jw95"] - r.cv["operational"] > 20.0, r.cv


@pytest.mark.slow
def test_the_report_names_the_conventions_that_were_not_fitted(quiet, peak_from_b):
    """Every absolute rate here is conditional on `FIT_CONVENTIONS`, and the result says so."""
    res = compare_models(ALL_THREE, [peak_from_b], n_folds=4, n_starts=3, cv_starts=2)
    joined = " ".join(res.verdict.caveats)
    assert "FIT_CONVENTIONS" in joined
    assert "only the ratios" in joined
    assert "synthetic" in joined or "formula" in joined
    text = str(comparison_report(res))
    assert res.cv_winner in text
    # two of the five conventions have no recorded origin and must say so
    unsourced = [k for k, (_, why) in FIT_CONVENTIONS.items() if "no recorded origin" in why]
    assert set(unsourced) == {"r", "r_fast", "r_slow"}
