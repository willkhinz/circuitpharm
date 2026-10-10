"""P3: identifiability as a question, answered before anything is fitted.

The full derivation and the measured numbers are in knowledge/11-identifiability.md. These
tests pin both halves: the STRUCTURAL result (what an equilibrium observable can determine
at all) and the PRACTICAL one (what it determines at realistic noise, which is much less).
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from circuitpharm.fitting.data import equilibrium_crc_from_model
from circuitpharm.fitting.identifiability import (
    EQUILIBRIUM_IDENTIFIABLE,
    IDENTIFIABLE_BOUNDS,
    IdentifiabilityClass,
    equilibrium_chi2_identifiable,
    fit_identifiable,
    identifiability_report,
    invariance_report,
    profile_likelihood,
)
from circuitpharm.fitting.reparam import UNLOCKED_BY, IdentifiableParams
from circuitpharm.models import KineticAllosteryModel
from circuitpharm.models.base import ObservableMismatch
from circuitpharm.results import Tier

#: Known-truth model for the method check. K_d = 25 uM, E = 4, D = 25.
TRUTH = KineticAllosteryModel(kon=0.01, koff=0.25, beta=0.8, alpha=0.2, d=0.05, r=0.002)
TRUTH_VEC = IdentifiableParams.from_microscopic(**TRUTH.get_params()).as_vector(
    with_timescale=False)


@pytest.fixture(scope="module")
def clean():
    return equilibrium_crc_from_model(TRUTH, noise_sd=0.0)


@pytest.fixture(scope="module")
def noisy():
    return equilibrium_crc_from_model(TRUTH, noise_sd=0.002)


# ============================================================ §1 the structural result
def test_the_reparameterisation_round_trips():
    """ANCHOR. from_microscopic -> to_microscopic must be exact in the ratios.

    rtol 1e-12: this is a pair of logarithms and their inverse, so only rounding is
    admissible.
    """
    p = IdentifiableParams.from_microscopic(**TRUTH.get_params())
    back = p.to_microscopic(alpha=TRUTH.alpha, r=TRUTH.r)
    for k, v in TRUTH.get_params().items():
        assert back[k] == pytest.approx(v, rel=1e-12), k


def test_the_map_is_many_to_one_and_says_so():
    """Scaling a pair must land on the same point: that IS the structural result."""
    a = IdentifiableParams.from_microscopic(**TRUTH.get_params(), with_timescale=False)
    scaled = dict(TRUTH.get_params())
    scaled["kon"] *= 7.0
    scaled["koff"] *= 7.0
    b = IdentifiableParams.from_microscopic(**scaled, with_timescale=False)
    assert b.as_vector(with_timescale=False) == pytest.approx(
        a.as_vector(with_timescale=False), rel=1e-12)

    # and recovering rates without the missing conventions must refuse
    with pytest.raises(ValueError, match="no absolute timescale"):
        a.to_microscopic(alpha=1.0, r=1.0)
    with pytest.raises(ValueError, match="alpha and r must be supplied"):
        IdentifiableParams.from_microscopic(**TRUTH.get_params()).to_microscopic()


def test_the_three_scaling_directions_are_flat():
    """ANCHOR. The invariance, measured. Spread < 1e-9 on an objective of order 1e1-1e4."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rep = invariance_report(TRUTH)
    assert set(rep["generators"]) == {"kon/koff", "beta/alpha", "d/r"}
    for name, info in rep["generators"].items():
        assert info["objective_spread"] < 1e-9, f"{name} is not flat"


def test_the_unlock_table_names_a_measurement_for_each_absolute_rate():
    """The ladder in §1.2 must be in code, not only in prose."""
    assert set(UNLOCKED_BY) == {"log10_koff", "log10_alpha", "log10_d"}
    assert "TIME COURSE" in UNLOCKED_BY["log10_koff"]
    assert "mean open time" in UNLOCKED_BY["log10_alpha"]
    assert "desensitisation" in UNLOCKED_BY["log10_d"]


def test_profiling_an_unidentifiable_rate_is_refused():
    with pytest.raises(ValueError, match="not an equilibrium-identifiable"):
        profile_likelihood("log10_koff")


# ============================================================ §2.2 the observable trap
def test_an_equilibrium_objective_refuses_a_peak_dataset():
    """ANCHOR on the refusal. Fitting equilibrium to peak data deletes desensitisation.

    The equilibrium plateau is E/(1+E+D) and a peak curve plateaus near 0.75, so the only
    way to reconcile them is D -> 0. Measured before the guard: log10_D -> -15.
    """
    from circuitpharm.fitting import DOSE_RESPONSE_BENCHMARK, JAHN1997_PEAK_CRC

    for ds in (DOSE_RESPONSE_BENCHMARK, JAHN1997_PEAK_CRC):
        with pytest.raises(ObservableMismatch, match="EQUILIBRIUM"):
            equilibrium_chi2_identifiable(TRUTH_VEC, ds)
        with pytest.raises(ObservableMismatch):
            fit_identifiable(ds)


# ============================================================ §2 the practical result
def test_noiseless_recovery_is_exact(clean):
    """ANCHOR, and the method validation. Without noise the MLE IS the truth.

    rtol 1e-4 on each log10 coordinate and a cost below 1e-12: if this fails the estimator
    is broken, and no conclusion downstream of it means anything. This is the test that
    earns the right to interpret the noisy case below.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p, cost = fit_identifiable(clean)
    assert cost < 1e-12
    assert p.as_vector(with_timescale=False) == pytest.approx(TRUTH_VEC, rel=1e-4)


def test_noiseless_profiles_are_identifiable_and_narrow(clean):
    """All three combinations bounded on both sides, each inside half a decade."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i, name in enumerate(EQUILIBRIUM_IDENTIFIABLE):
            r = profile_likelihood(name, clean, n_grid=15)
            assert r.classification == IdentifiabilityClass.IDENTIFIABLE, r.reason
            width = r.ci_95[1] - r.ci_95[0]
            assert width < 0.6, f"{name} CI is {width:.3f} decades"
            assert r.ci_95[0] <= TRUTH_VEC[i] <= r.ci_95[1], (
                f"{name}'s interval excludes the truth")


def test_two_tenths_of_a_percent_of_noise_destroys_E_and_D(noisy):
    """ANCHOR, and the result P3 exists to produce.

    At 0.2% noise on a NORMALISED curve, E and D span more than two decades -- a factor of
    ~130 -- and are classified practically non-identifiable. K_d survives, more loosely.

    The thresholds are deliberately coarse (>1.5 decades for E and D, <1.5 for K_d): the
    point is the order of magnitude of the collapse, not its third digit, and the exact
    widths depend on the noise realisation. The measured values are 2.115, 2.137 and 1.335.

    THE CONSEQUENCE, which is the actionable part: more precision on an equilibrium
    concentration-response cannot fix this, because two parameter sets a hundredfold apart
    produce curves separated by less than the noise. A different OBSERVABLE is needed --
    see reparam.UNLOCKED_BY and fitting.data.MISSING_DATASETS.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        widths, classes = {}, {}
        for name in EQUILIBRIUM_IDENTIFIABLE:
            r = profile_likelihood(name, noisy, n_grid=15)
            widths[name] = r.ci_95[1] - r.ci_95[0]
            classes[name] = r.classification

    assert widths["log10_E"] > 1.5, widths
    assert widths["log10_D"] > 1.5, widths
    assert classes["log10_E"] == IdentifiabilityClass.PRACTICALLY_NON_IDENTIFIABLE
    assert classes["log10_D"] == IdentifiabilityClass.PRACTICALLY_NON_IDENTIFIABLE
    assert widths["log10_kd"] < 1.5, widths


def test_a_wrong_parameter_set_fits_the_noise_better_than_the_truth(noisy):
    """ANCHOR. The sharpest statement of practical non-identifiability.

    The MLE is 2-4x away from the truth in every coordinate and fits the noise realisation
    BETTER than the truth does (measured: 4.96 against 11.41), while the two predicted
    curves differ by at most 2.4e-3 -- the noise. So the data do not merely fail to pin the
    parameters down; they actively prefer the wrong ones.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p, cost = fit_identifiable(noisy)
    cost_truth = equilibrium_chi2_identifiable(TRUTH_VEC, noisy)
    assert cost < cost_truth

    fitted = KineticAllosteryModel(**p.to_microscopic(koff=1.0, alpha=1.0, r=1e-3))
    y_t = TRUTH.dose_response(noisy.concs_um)
    y_f = fitted.dose_response(noisy.concs_um)
    assert float(np.abs(y_t - y_f).max()) < 5e-3, (
        "the two curves are distinguishable, so this is not a flat-valley result")
    assert float(np.abs(p.as_vector(with_timescale=False) - TRUTH_VEC).max()) > 0.3


# ============================================================ reporting
def test_the_report_is_tiered_and_the_structural_count_is_validated():
    """The one VALIDATED claim in this layer: three of six, proved not fitted."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rs = identifiability_report(equilibrium_crc_from_model(TRUTH, noise_sd=0.0))
    q = rs.quantity("n_identifiable_of_six")
    assert q.tier is Tier.VALIDATED
    assert q.value == 3
    assert "provable" in q.provenance
    for name in EQUILIBRIUM_IDENTIFIABLE:
        qq = rs.quantity(name)
        assert qq.tier is Tier.UNCALIBRATED, (
            "the method is verified but the dataset is generated, so nothing here is a "
            "measurement")
        assert any("not a digitisation" in c for c in qq.caveats)


def test_bounds_are_physical_and_stated():
    for name in EQUILIBRIUM_IDENTIFIABLE:
        lo, hi = IDENTIFIABLE_BOUNDS[name]
        assert lo < hi and hi - lo >= 4.0, f"{name}'s bounds are too narrow to be honest"
