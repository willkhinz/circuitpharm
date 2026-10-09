"""P0-5 and P0-6: defective machinery stays callable, loudly, with VOID conclusions.

The contract these tests pin (roadmap §2.4's last paragraph):

  * every call emits `ProvisionalResultWarning`, on EVERY call and not just the first;
  * raw computed arrays stay plain and usable;
  * inferential summaries are `Tier.VOID`, so `.value` raises and
    `.get(acknowledge_void=True)` works at a call site a reviewer can see;
  * a synthetic dataset taints whatever consumes it, automatically.

`test_warning_fires_on_every_call_not_just_the_first` is the load-bearing one: without the
`simplefilter("always")` in `provisional.py`, Python shows a UserWarning once per code
location per process, so the banner reaches whoever makes the first call and nobody after.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from circuitpharm.fitting import (
    DEACTIVATION_BENCHMARK,
    DOSE_RESPONSE_BENCHMARK,
    HOLDOUT,
    TRAIN,
    assert_real_data,
    compute_profile_likelihood,
    cost_hessian_and_spectrum,
    equilibrium_dose_response_chi2,
    holdout_guard,
    run_ensemble_mcmc,
)
from circuitpharm.fitting.data import DoseResponseDataset, _validate_dataset
from circuitpharm.fitting.identifiability import invariance_report
from circuitpharm.models import KineticAllosteryModel
from circuitpharm.models.base import Observable
from circuitpharm.provisional import ProvisionalResultWarning, warn_provisional
from circuitpharm.results import Quantity, Tier, VoidQuantityError

MODEL = KineticAllosteryModel()
X0 = np.array([MODEL.kon, MODEL.koff, MODEL.beta, MODEL.alpha, MODEL.d, MODEL.r])


# ======================================================================== P0-5 layer 1
def test_warning_fires_on_every_call_not_just_the_first():
    """PERMANENT RECORD. The `simplefilter("always")` in provisional.py.

    Python's default filter shows a UserWarning once per unique location per process, so
    the banner would reach the first caller and be invisible to the loop that generates the
    table that goes in the write-up. Two calls must produce two records.
    """
    with warnings.catch_warnings(record=True) as rec:
        warnings.resetwarnings()                     # drop any pre-set filters
        import circuitpharm.provisional as prov
        warnings.simplefilter("always", prov.ProvisionalResultWarning)
        for _ in range(3):
            warn_provisional(what="X", defect="d", register_item="P0-0", promote_by="p")
    fired = [r for r in rec if issubclass(r.category, ProvisionalResultWarning)]
    assert len(fired) == 3, f"expected 3 warnings, got {len(fired)}"


def test_warn_provisional_refuses_a_vague_banner():
    """A banner saying "results may be unreliable" trains people to ignore the category."""
    for kwargs in (dict(what="", defect="d", register_item="P0", promote_by="p"),
                   dict(what="w", defect="   ", register_item="P0", promote_by="p"),
                   dict(what="w", defect="d", register_item="", promote_by="p"),
                   dict(what="w", defect="d", register_item="P0", promote_by="")):
        with pytest.raises(ValueError, match="specific"):
            warn_provisional(**kwargs)


def test_banner_names_the_defect_and_the_promotion_path():
    with pytest.warns(ProvisionalResultWarning) as rec:
        compute_profile_likelihood("kon", MODEL, scan_factors=(0.8, 1.0, 1.25))
    text = "\n".join(str(r.message) for r in rec)
    assert "K_d = k_off/k_on" in text
    assert "P0-5" in text
    assert "reparameterise" in text


# ======================================================================== P0-5 layer 2
def test_provisional_conclusions_are_void_and_raw_arrays_are_not():
    """ANCHOR for the VOID/plain split in the register's table.

    Both halves must hold. Marking everything VOID passes a one-sided test trivially, and
    so does marking nothing; the split is the content.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        prof = compute_profile_likelihood("koff", MODEL, scan_factors=(0.5, 1.0, 2.0))
        spec = cost_hessian_and_spectrum(MODEL)
        chain = run_ensemble_mcmc(X0, n_walkers=12, n_steps=60, burn_in=20, seed=7)

    # --- plain: real computations, readable without acknowledgement
    assert isinstance(prof.grid_values, np.ndarray) and np.all(np.isfinite(prof.grid_values))
    assert isinstance(prof.profile_costs, np.ndarray)
    assert np.all(np.isfinite(prof.profile_costs))
    assert spec.hessian.shape == (6, 6) and np.all(np.isfinite(spec.hessian))
    assert spec.eigenvalues.shape == (6,)
    assert chain.flat_samples.shape[1] == 6
    assert 0.0 <= chain.acceptance_fraction <= 1.0
    assert isinstance(chain.diagnostics, dict)

    # --- VOID: inferences the objective cannot support
    voids = [prof.ci_95_bounds, prof.is_identifiable, prof.diagnostic_reason,
             spec.condition_number]
    voids += list(chain.posterior_means.values())
    voids += list(chain.posterior_std.values())
    voids += list(chain.credible_intervals_95.values())
    for q in voids:
        assert isinstance(q, Quantity)
        assert q.tier is Tier.VOID
        with pytest.raises(VoidQuantityError):
            q.value
        q.get(acknowledge_void=True)        # must not raise


def test_void_messages_name_the_reason():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = cost_hessian_and_spectrum(MODEL)
    with pytest.raises(VoidQuantityError) as exc:
        spec.condition_number.value
    msg = str(exc.value)
    assert "K_d = k_off/k_on" in msg
    assert "P0-5" in msg


def test_full_eigenvalue_spectrum_is_returned():
    """Negatives are the finding, not noise to filter.

    The previous implementation took the condition number over positive eigenvalues only,
    which hid both the flat subspace and the fact that the parameters are not at a minimum.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = cost_hessian_and_spectrum(MODEL)
    assert spec.eigenvalues.size == 6
    assert np.any(spec.eigenvalues < 0), "a negative eigenvalue was expected here"
    # the caveat must say so, since the number itself looks like a normal condition number
    caveats = " ".join(spec.condition_number.caveats)
    assert "NEGATIVE" in caveats and "not a Fisher information matrix" in caveats
    assert len(spec.near_null_directions) == 3


# ======================================================================== P0-5 the reason
def test_objective_is_invariant_to_uniform_rate_scaling():
    """ANCHOR, PERMANENT RECORD. The structural reason the conclusions are VOID.

    rtol=1e-12: the invariance is exact algebra, not an approximation, so the only
    admissible difference is float rounding. Amend this test when P3 lands a likelihood
    with a kinetic observable -- do NOT delete it; it is the only durable statement of why
    P3 exists.
    """
    base = equilibrium_dose_response_chi2(X0)
    for factor in (2.0, 10.0, 100.0, 1000.0):
        assert equilibrium_dose_response_chi2(X0 * factor) == pytest.approx(base, rel=1e-12)


def test_the_three_invariance_generators_are_flat():
    """ANCHOR. Each PAIR scaling independently leaves the objective unchanged.

    This is stronger and more specific than the uniform-scaling test: it identifies WHICH
    three directions carry no information, which is what P3's reparameterisation is built
    from. Spread tolerance 1e-9 absolute on an objective of order 1e4.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rep = invariance_report(MODEL)
    assert set(rep["generators"]) == {"kon/koff", "beta/alpha", "d/r"}
    for name, info in rep["generators"].items():
        assert info["objective_spread"] < 1e-9, f"{name} is not flat"


def test_zero_curvature_subspace_is_three_dimensional():
    """ANCHOR. v^T H v / max|H| ~ 0 for all three generators.

    1e-8 relative: the Hessian is finite-difference, so an exact zero is not available;
    the measured values are ~1e-12 to 1e-11, four orders inside this bound. Note this is a
    CURVATURE statement, not a rank statement -- see `invariance_report`'s docstring for
    why the roadmap's earlier "rank <= 3" claim was withdrawn.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = cost_hessian_and_spectrum(MODEL)
    keys = list(spec.param_names)
    theta = np.array([getattr(MODEL, k) for k in keys])
    scale = float(np.abs(spec.hessian).max())
    for a, b in (("kon", "koff"), ("beta", "alpha"), ("d", "r")):
        v = np.zeros(6)
        v[keys.index(a)] = theta[keys.index(a)]
        v[keys.index(b)] = theta[keys.index(b)]
        u = v / np.linalg.norm(v)
        assert abs(float(u @ spec.hessian @ u)) / scale < 1e-8


# ======================================================================== P0-12 fitting
def test_optimiser_failure_raises_rather_than_scoring_as_identifiable():
    """Sign-inverted defect: `base_cost + 10.0` is ABOVE the 3.84 threshold."""
    import circuitpharm.fitting.identifiability as idn

    class _Failed:
        success = False
        message = "SENTINEL_MINIMIZE_FAILURE"
        fun = 1.0
        x = np.ones(5)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        orig = idn.minimize
        try:
            idn.minimize = lambda *a, **k: _Failed()
            with pytest.raises(RuntimeError, match="SENTINEL_MINIMIZE_FAILURE"):
                idn.compute_profile_likelihood("kon", MODEL, scan_factors=(1.0,))
        finally:
            idn.minimize = orig


def test_profile_baseline_is_a_real_optimum():
    """`base_cost` must come from re-optimising, or delta can be negative and the CI logic
    is not a profile likelihood at all."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        prof = compute_profile_likelihood("kon", MODEL, scan_factors=(0.5, 1.0, 2.0))
    nominal = equilibrium_dose_response_chi2(X0)
    assert float(prof.profile_costs.min()) < nominal, (
        "the profile never beat the cost at the nominal parameters, so the baseline was "
        "probably not re-optimised")


def test_log_prior_is_log_uniform_in_the_sampling_space():
    """The old `log_prior` claimed log-uniform and was uniform on the rates."""
    from circuitpharm.fitting import LOG10_BOUNDS, log_prior_log10
    keys = ("kon", "koff", "beta", "alpha", "d", "r")
    inside = np.array([(lo + hi) / 2 for lo, hi in (LOG10_BOUNDS[k] for k in keys)])
    assert log_prior_log10(inside, keys) == 0.0
    outside = inside.copy()
    outside[0] = LOG10_BOUNDS["kon"][1] + 1.0
    assert log_prior_log10(outside, keys) == -np.inf


def test_mcmc_diagnostics_gate_the_result():
    """ANCHOR. A deliberately under-run chain must report itself unconverged.

    60 steps over 12 walkers in 6 dimensions cannot support an autocorrelation estimate,
    let alone an interval. The old defaults (16 x 150) were in the same regime and nothing
    said so.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        chain = run_ensemble_mcmc(X0, n_walkers=12, n_steps=60, burn_in=20, seed=3)
    assert chain.converged is False
    assert chain.diagnostics["failures"], "an under-run chain reported no failures"
    joined = " ".join(" ".join(q.caveats) for q in chain.credible_intervals_95.values())
    assert "diagnostics failed" in joined


def test_mcmc_refuses_a_degenerate_start():
    bad = X0.copy()
    bad[0] = 1e3            # far outside LOG10_BOUNDS["kon"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match="outside LOG10_BOUNDS"):
            run_ensemble_mcmc(bad, n_walkers=12, n_steps=10, burn_in=2)


# ======================================================================== P0-6
def test_synthetic_datasets_declare_themselves():
    for ds in (DOSE_RESPONSE_BENCHMARK, DEACTIVATION_BENCHMARK):
        assert ds.synthetic is True
        assert ds.origin.strip()
        assert ds.source_key == "", "a generated trace must not carry a provenance key"


def test_no_citation_field_on_synthetic_data():
    """PERMANENT RECORD. `citation=` on a formula asserts provenance that does not exist.

    The deactivation trace is `0.70*exp(-t/15) + 0.30*exp(-t/70)` and carried
    `citation="Haas & Macdonald 1999 ..."` -- the paper knowledge/07-paper-review.md
    records as measuring 76.1 ms for that quantity. The text survives under
    `motivated_by`; the field name must not come back by habit.
    """
    for ds in (DOSE_RESPONSE_BENCHMARK, DEACTIVATION_BENCHMARK):
        assert not hasattr(ds, "citation")
        assert ds.motivated_by.strip()
    assert "76.1 ms" in DEACTIVATION_BENCHMARK.motivated_by


def test_module_docstring_declares_synthetic():
    """Blunt on purpose -- the same technique test_manuscript_consistency.py uses."""
    import circuitpharm.fitting.data as data
    first = (data.__doc__ or "").strip().splitlines()[0]
    assert "SYNTHETIC" in first


def test_a_wholly_invented_dataset_cannot_carry_a_source():
    """`kind="synthetic"` means invented, so a source_key on it asserts provenance that
    does not exist. A PARAMETRIC dataset may carry one -- its parameters are real -- which
    is why the three-way `DataKind` exists rather than a synthetic/not-synthetic flag."""
    with pytest.raises(ValueError, match="wholly generated trace has no provenance"):
        DoseResponseDataset(
            citation_label="x", preparation="p",
            concs_um=np.array([1.0]), mean_response=np.array([0.1]),
            sem=np.array([0.01]), kind="synthetic", synthetic=True, origin="formula",
            source_key="some_key")


def test_kind_and_synthetic_cannot_disagree():
    """Two fields describing one property must not be settable into contradiction."""
    with pytest.raises(ValueError, match="implies synthetic"):
        DoseResponseDataset(
            citation_label="x", preparation="p", concs_um=np.array([1.0]),
            mean_response=np.array([0.1]), sem=np.array([0.01]),
            kind="digitised", synthetic=True, origin="o", source_key="k",
            figure="f", digitisation="d")
    with pytest.raises(ValueError, match="implies synthetic"):
        DoseResponseDataset(
            citation_label="x", preparation="p", concs_um=np.array([1.0]),
            mean_response=np.array([0.1]), sem=np.array([0.01]),
            kind="parametric", synthetic=False, origin="o", source_key="k")


def test_parametric_data_must_name_its_source_and_its_construction():
    for kw, match in (
        (dict(kind="parametric", synthetic=True, origin="o"), "names no `source_key`"),
        (dict(kind="parametric", synthetic=True, source_key="k"), "no `origin` stating"),
    ):
        with pytest.raises(ValueError, match=match):
            DoseResponseDataset(
                citation_label="x", preparation="p", concs_um=np.array([1.0]),
                mean_response=np.array([0.1]), sem=np.array([0.01]), **kw)


def test_real_data_must_name_its_source_and_figure():
    with pytest.raises(ValueError, match="names no `source_key`"):
        DoseResponseDataset(
            citation_label="x", preparation="p", concs_um=np.array([1.0]),
            mean_response=np.array([0.1]), sem=np.array([0.01]), synthetic=False)
    with pytest.raises(ValueError, match="which `figure`"):
        DoseResponseDataset(
            citation_label="x", preparation="p", concs_um=np.array([1.0]),
            mean_response=np.array([0.1]), sem=np.array([0.01]), synthetic=False,
            source_key="k")


def test_synthetic_datasets_taint_their_results():
    """ANCHOR. Consuming a synthetic dataset warns and names it."""
    with pytest.warns(ProvisionalResultWarning) as rec:
        fake = assert_real_data([DOSE_RESPONSE_BENCHMARK, DEACTIVATION_BENCHMARK],
                                caller="test")
    assert len(fake) == 2
    text = "\n".join(str(r.message) for r in rec)
    assert "synthetic alpha1beta2gamma2 peak concentration-response" in text
    assert "synthetic biexponential deactivation trace" in text
    assert "P0-6" in text


def test_real_data_does_not_warn():
    """The taint must switch itself off when P1-2 lands, with no call site changed."""
    real = DoseResponseDataset(
        citation_label="real", preparation="p", concs_um=np.array([1.0]),
        mean_response=np.array([0.1]), sem=np.array([0.01]),
        synthetic=False, source_key="k", figure="Fig 1A",
        digitisation="by hand, 2026-10-09")
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        assert assert_real_data([real], caller="test") == []
    assert not [r for r in rec if issubclass(r.category, ProvisionalResultWarning)]


def test_holdout_is_unreachable_from_fitting():
    holdout = DoseResponseDataset(
        citation_label="held out", preparation="p", concs_um=np.array([1.0]),
        mean_response=np.array([0.1]), sem=np.array([0.01]),
        synthetic=False, source_key="k", figure="Fig 2", digitisation="x", role="holdout")
    with pytest.raises(ValueError, match="destroys both"):
        holdout_guard([holdout], caller="test")
    holdout_guard(list(TRAIN), caller="test")        # train data passes


def test_datasets_declare_an_observable():
    assert DOSE_RESPONSE_BENCHMARK.observable is Observable.PEAK
    assert DEACTIVATION_BENCHMARK.observable is Observable.CHARGE


# ======================================================================== P1-2 real data
def test_the_jahn_dataset_is_parametric_and_sourced():
    """ANCHOR. The project's first quantitatively sourced kinetic dataset.

    The published parameters are EC50 = 11.6 +/- 0.9 uM and nH = 2.2 +/- 0.4 (Jahn et al.
    1997, NeuroReport 8(16):3443-6). The generated curve must reproduce that EC50 exactly
    at its half-maximum -- it is constructed from it, so anything but equality means the
    construction is wrong.
    """
    from circuitpharm.fitting import JAHN1997_PEAK_CRC as ds

    assert ds.kind == "parametric"
    assert ds.synthetic is True, "a parametric curve must still taint a fit"
    assert ds.source_key == "jahn_a1b2g2_kinetics"
    assert ds.observable is Observable.PEAK
    assert ds.role == "train"
    assert "NOT digitised" in ds.digitisation
    assert ds.motivated_by == "", "a sourced dataset uses source_key, not motivated_by"

    # The curve is a Hill function with the published parameters, so recover BOTH of them
    # by regression in logit space, where log(y/(1-y)) = nH * (log c - log EC50) is exactly
    # linear. Doing it this way rather than by interpolating the half-maximum on the
    # concentration grid matters: the grid is decade-spaced, and a log-linear interpolation
    # between 3 and 30 uM returns 12.08 uM for a curve whose EC50 is exactly 11.6 -- a 4%
    # artefact of the interpolation that would have looked like a 4% construction error.
    c, y = ds.concs_um, ds.mean_response
    m = (y > 1e-6) & (y < 1.0 - 1e-6)
    slope, intercept = np.polyfit(np.log10(c[m]), np.log10(y[m] / (1.0 - y[m])), 1)
    ec50 = float(10.0 ** (-intercept / slope))
    assert slope == pytest.approx(2.2, rel=1e-9), "the published Hill slope is not recovered"
    assert ec50 == pytest.approx(11.6, rel=1e-9), "the published EC50 is not recovered"
    # The Hill asymptote is 1 and is APPROACHED, not reached: at the paper's 3 mM
    # saturating concentration the curve is 0.999995. "Saturates" means the last decade
    # buys nothing measurable, not that the formula equals its limit.
    assert y.max() == pytest.approx(1.0, abs=1e-4)
    assert y.max() < 1.0
    assert np.all(np.diff(y) > 0), "a concentration-response curve must be monotone"
    assert np.all(ds.sem >= 0.01), "the SEM floor is not being applied"


def test_the_jahn_source_resolves_in_the_knowledge_base():
    """A source_key that names nothing is worse than no key at all."""
    import sqlite3

    con = sqlite3.connect("data/pharmacology.db")
    rows = list(con.execute("select citation, verification from sources where key = ?",
                            ("jahn_a1b2g2_kinetics",)))
    assert rows, "jahn_a1b2g2_kinetics is not in the sources table"
    citation, verification = rows[0]
    assert "NeuroReport" in citation and "3443" in citation
    # the record must say plainly what was and was not read
    assert "FULL TEXT NOT READ" in verification
    assert "CLAIM_SUPPORT_DOES_NOT_SUPPORT for MEAN OPEN TIME" in verification


def test_the_scheme_cannot_reproduce_the_measured_hill_slope():
    """ANCHOR, and the first falsifiable mismatch against a sourced measurement.

    Jahn et al. report a Hill-type slope of 2.2 +/- 0.4 and read it as evidence for at
    least three binding sites. The 5-state scheme has TWO, and its peak curve comes out at
    1.294 -- below even the lower error bound of 1.8. This is structural: no choice of
    rates makes a two-site scheme that steep, so it cannot be fitted away, and P5's model
    comparison is where it gets adjudicated.

    Tolerance: asserted as an inequality against the measured lower bound rather than as a
    value, because the point is the direction and the size of the gap, not the third digit.
    """
    from circuitpharm.config import SYNAPTIC_PULSE
    from circuitpharm.gabaa_kinetics import fit_scheme

    s = fit_scheme(verbose=False, pulse=dict(SYNAPTIC_PULSE))
    c = np.logspace(-1, 4, 60)
    y = np.array([s.po_peak(x) for x in c])
    f = y / y.max()
    m = (f > 0.1) & (f < 0.9)
    slope = float(np.polyfit(np.log10(c[m]), np.log10(f[m] / (1 - f[m])), 1)[0])

    assert slope == pytest.approx(1.294, rel=0.05)
    assert slope < 2.2 - 0.4, (
        f"the scheme's peak Hill slope is {slope:.3f}; Jahn et al. measured 2.2 +/- 0.4. "
        f"If this ever rises above 1.8 the structural argument in fitting/data.py needs "
        f"re-examining.")


def test_the_missing_datasets_are_named_not_forgotten():
    """The gap is the finding, so it is in the code rather than only in a commit message."""
    from circuitpharm.fitting import MISSING_DATASETS

    assert set(MISSING_DATASETS) == {"deactivation_peak_pulse",
                                     "single_channel_mean_open_time", "holdout"}
    for key, note in MISSING_DATASETS.items():
        assert len(note) > 100, f"{key}'s note does not say what is needed or why"
    assert "k_off" in MISSING_DATASETS["deactivation_peak_pulse"]
    assert "FIT_FIXED_ALPHA" in MISSING_DATASETS["single_channel_mean_open_time"]
