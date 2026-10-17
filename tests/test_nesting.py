"""Machine-checked nesting tests (Brief 6 Task B).

Verifies that for every declared (sub, big) pair in MODEL_SPECS:
  1. The embedding preserves total log-likelihood: ll_embedded == approx(ll_sub).
  2. The enlarged model scores at least as well on total log-likelihood: ll_big >= ll_sub - tol.
  3. No assertions are made on per-dataset residuals (AGENT-BRIEF-6 §2).
"""
import numpy as np
import pytest

from circuitpharm.fitting.comparison import FIT_CONVENTIONS, MODEL_SPECS, nesting_embedding
from circuitpharm.fitting.data import (
    BARBERIS2007_DEACT,
    BARBERIS2007_PAIRED_PULSE,
    JAHN1997_PEAK_CRC,
)
from circuitpharm.fitting.diagnostic import (
    diagnostic_log_likelihood,
    fit_diagnostic_consistency,
)


def test_declared_nesting_topology():
    """Verify declared nesting relations in MODEL_SPECS."""
    assert MODEL_SPECS["kinetic_jw95_free_r"].nests == ("kinetic_jw95",)
    assert MODEL_SPECS["kinetic_jw95_free_rates"].nests == ("kinetic_jw95", "kinetic_jw95_free_r")
    assert MODEL_SPECS["kinetic_jw95"].nests == ()
    assert MODEL_SPECS["operational"].nests == ()
    assert MODEL_SPECS["kinetic_open_desens"].nests == ()
    assert MODEL_SPECS["kinetic_open_desens_free_rates"].nests == ("kinetic_open_desens",)
    assert MODEL_SPECS["extended_desens"].nests == ()
    assert MODEL_SPECS["kinetic_flip"].nests == ()
    assert MODEL_SPECS["kinetic_two_open"].nests == ()


def test_extended_desens_exclusion_with_measured_discrepancy():
    """Verify extended_desens is excluded from nesting kinetic_jw95 with measured discrepancy.

    Per AGENT-BRIEF-6 §3 Task B.2: extended_desens does NOT strictly nest kinetic_jw95 under
    FIT_CONVENTIONS because:
      1. r_fast convention is 0.0050 vs r=0.0020 (differs by 2.5x).
      2. D_slow -> 0 sits on the search bound (log10_D_slow = -4.0).
    The measured discrepancy is recorded here and verified to be non-zero.
    """
    v3 = np.array([-0.428, 0.104, 1.9721])
    m3 = MODEL_SPECS["kinetic_jw95"].factory(v3)
    v_c = np.array([-0.428, 0.104, 1.9721, -4.0])
    mc = MODEL_SPECS["extended_desens"].factory(v_c)

    dsets_dr = [BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]
    ll3_dr = diagnostic_log_likelihood(m3, dsets_dr, n_params=3, anchor_po_max=True).total_log_likelihood
    llc_dr = diagnostic_log_likelihood(mc, dsets_dr, n_params=4, anchor_po_max=True).total_log_likelihood
    discrepancy_dr = abs(ll3_dr - llc_dr)

    dsets_j = [JAHN1997_PEAK_CRC]
    ll3_j = diagnostic_log_likelihood(m3, dsets_j, n_params=3).total_log_likelihood
    llc_j = diagnostic_log_likelihood(mc, dsets_j, n_params=4).total_log_likelihood
    discrepancy_j = abs(ll3_j - llc_j)

    # Record measured discrepancies
    assert discrepancy_dr > 10.0, f"Discrepancy on Deact+Recovery should be substantial (~16.66 nats), got {discrepancy_dr:.2f}"
    assert discrepancy_j > 0.1, f"Discrepancy on Jahn peak should be non-zero (~0.28 nats), got {discrepancy_j:.2f}"
    assert MODEL_SPECS["extended_desens"].nests == ()


def test_kinetic_flip_exclusion_with_measured_discrepancy():
    """Verify kinetic_flip is excluded from nesting kinetic_jw95 with measured discrepancy.

    Per AGENT-BRIEF-8B Task 2: kinetic_flip does NOT strictly nest kinetic_jw95 under
    finite parameter bounds because:
      1. The theoretical limit collapsing A2F onto A2R requires delta -> inf with gamma fixed.
         In any bounded parameter box (e.g. log10_delta <= 3.0), delta sits on the upper bound.
      2. Rapid flipping (delta >> gamma) shifts occupancy to A2F, bypassing desensitisation
         from A2R. At delta=1000, gamma=1 (log10_delta=3.0, log10_gamma=0.0) with conventions
         r=0.002, koff=0.180415, alpha=0.30, the measured discrepancy is ~161.60 nats on
         Deact+Recovery and ~3.37 nats on Jahn peak.
    The measured discrepancy is recorded here and verified to be non-zero.
    """
    v3 = np.array([-0.428, 0.104, 1.9721])
    m3 = MODEL_SPECS["kinetic_jw95"].factory(v3)

    v_flip = np.array([
        -0.428, 0.104, 1.9721,
        float(np.log10(FIT_CONVENTIONS["r"][0])),
        float(np.log10(FIT_CONVENTIONS["koff"][0])),
        float(np.log10(FIT_CONVENTIONS["alpha"][0])),
        3.0, 0.0,
    ])
    mf = MODEL_SPECS["kinetic_flip"].factory(v_flip)

    dsets_dr = [BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]
    ll3_dr = diagnostic_log_likelihood(m3, dsets_dr, n_params=3, anchor_po_max=True).total_log_likelihood
    llf_dr = diagnostic_log_likelihood(mf, dsets_dr, n_params=8, anchor_po_max=True).total_log_likelihood
    discrepancy_dr = abs(ll3_dr - llf_dr)

    dsets_j = [JAHN1997_PEAK_CRC]
    ll3_j = diagnostic_log_likelihood(m3, dsets_j, n_params=3).total_log_likelihood
    llf_j = diagnostic_log_likelihood(mf, dsets_j, n_params=8).total_log_likelihood
    discrepancy_j = abs(ll3_j - llf_j)

    # Record measured discrepancies
    assert discrepancy_dr > 10.0, f"Discrepancy on Deact+Recovery should be substantial (~161.60 nats), got {discrepancy_dr:.2f}"
    assert discrepancy_j > 0.1, f"Discrepancy on Jahn peak should be non-zero (~3.37 nats), got {discrepancy_j:.2f}"
    assert MODEL_SPECS["kinetic_flip"].nests == ()


def test_kinetic_two_open_exclusion_with_measured_discrepancy():
    """Verify kinetic_two_open is excluded from nesting kinetic_jw95 with measured discrepancy.

    Per AGENT-BRIEF-9 Task 3: kinetic_two_open does NOT strictly nest kinetic_jw95 under
    finite parameter bounds because:
      1. The theoretical limit collapsing A2O' requires delta -> 0 (log10_delta -> -inf).
         In any bounded parameter box (e.g. log10_delta >= -3.0), delta sits on the lower bound.
      2. At delta=1e-3, gamma=1.0 (log10_delta=-3.0, log10_gamma=0.0) with conventions
         r=0.002, koff=0.180415, alpha=0.30, the measured discrepancy is ~0.005 nats on
         Deact+Recovery and ~0.018 nats on Jahn peak.
    The measured discrepancy is recorded here and verified to be non-zero, and nests is empty.
    """
    v3 = np.array([-0.428, 0.104, 1.9721])
    m3 = MODEL_SPECS["kinetic_jw95"].factory(v3)

    v_2o = np.array([
        -0.428, 0.104, 1.9721,
        float(np.log10(FIT_CONVENTIONS["r"][0])),
        float(np.log10(FIT_CONVENTIONS["koff"][0])),
        float(np.log10(FIT_CONVENTIONS["alpha"][0])),
        -3.0, 0.0,
    ])
    m_2o = MODEL_SPECS["kinetic_two_open"].factory(v_2o)

    dsets_dr = [BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]
    ll3_dr = diagnostic_log_likelihood(m3, dsets_dr, n_params=3, anchor_po_max=True).total_log_likelihood
    ll2o_dr = diagnostic_log_likelihood(m_2o, dsets_dr, n_params=8, anchor_po_max=True).total_log_likelihood
    discrepancy_dr = abs(ll3_dr - ll2o_dr)

    dsets_j = [JAHN1997_PEAK_CRC]
    ll3_j = diagnostic_log_likelihood(m3, dsets_j, n_params=3).total_log_likelihood
    ll2o_j = diagnostic_log_likelihood(m_2o, dsets_j, n_params=8).total_log_likelihood
    discrepancy_j = abs(ll3_j - ll2o_j)

    # Record measured discrepancies
    assert discrepancy_dr > 0.001, f"Discrepancy on Deact+Recovery should be non-zero (~0.005 nats), got {discrepancy_dr:.4f}"
    assert discrepancy_j > 0.001, f"Discrepancy on Jahn peak should be non-zero (~0.018 nats), got {discrepancy_j:.4f}"
    assert MODEL_SPECS["kinetic_two_open"].nests == ()


def test_nesting_embedding_derives_conventions_by_name():
    """Verify nesting_embedding correctly maps parameters from FIT_CONVENTIONS by name."""
    sub = MODEL_SPECS["kinetic_jw95"]
    big = MODEL_SPECS["kinetic_jw95_free_rates"]
    v_sub = np.array([-0.5, 0.2, 1.5])
    v_big = nesting_embedding(sub, big, v_sub)

    assert v_big.shape == (6,)
    assert v_big[0] == pytest.approx(-0.5)
    assert v_big[1] == pytest.approx(0.2)
    assert v_big[2] == pytest.approx(1.5)
    # Rate slots derived from FIT_CONVENTIONS by name, never literals
    assert v_big[3] == pytest.approx(float(np.log10(FIT_CONVENTIONS["r"][0])))
    assert v_big[4] == pytest.approx(float(np.log10(FIT_CONVENTIONS["koff"][0])))
    assert v_big[5] == pytest.approx(float(np.log10(FIT_CONVENTIONS["alpha"][0])))

    # Also check kinetic_open_desens into kinetic_open_desens_free_rates
    sub_od = MODEL_SPECS["kinetic_open_desens"]
    big_od = MODEL_SPECS["kinetic_open_desens_free_rates"]
    v_big_od = nesting_embedding(sub_od, big_od, v_sub)
    assert v_big_od.shape == (6,)
    assert v_big_od[0] == pytest.approx(-0.5)
    assert v_big_od[1] == pytest.approx(0.2)
    assert v_big_od[2] == pytest.approx(1.5)
    assert v_big_od[3] == pytest.approx(float(np.log10(FIT_CONVENTIONS["r"][0])))
    assert v_big_od[4] == pytest.approx(float(np.log10(FIT_CONVENTIONS["koff"][0])))
    assert v_big_od[5] == pytest.approx(float(np.log10(FIT_CONVENTIONS["alpha"][0])))


@pytest.mark.slow
@pytest.mark.parametrize("sub_name,big_name", [
    ("kinetic_jw95", "kinetic_jw95_free_r"),
    ("kinetic_jw95", "kinetic_jw95_free_rates"),
    ("kinetic_jw95_free_r", "kinetic_jw95_free_rates"),
])
@pytest.mark.parametrize("dset_name,dsets", [
    ("deact_recovery", [BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]),
    ("joint_all_three", [JAHN1997_PEAK_CRC, BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]),
])
def test_nesting_on_total_objective(sub_name, big_name, dset_name, dsets):
    """Machine-check that enlarged models fit >= sub-models on total log-likelihood.

    IMPORTANT: Per AGENT-BRIEF-6 §2, nesting constrains ONLY total log-likelihood.
    Nesting does NOT constrain per-dataset residuals. An enlarged model is free
    to reallocate error between datasets. Under this diagnostic's scoring (profiled
    sigma over 100 points vs fixed SEM on 1 point), buying trace fit by sacrificing
    the recovery point is mathematically favoured (69 nats gained vs 47 nats paid).
    Therefore, this test deliberately asserts NOTHING about per-dataset residuals.
    """
    spec_jw95 = MODEL_SPECS["kinetic_jw95"]
    spec_free_r = MODEL_SPECS["kinetic_jw95_free_r"]
    sub = MODEL_SPECS[sub_name]
    big = MODEL_SPECS[big_name]

    # Reference point from AGENT-BRIEF-6 §1 (clean benchmark point)
    v3 = np.array([-0.428, 0.104, 1.9721])
    if sub_name == "kinetic_jw95":
        v_sub = v3
    elif sub_name == "kinetic_jw95_free_r":
        v_sub = nesting_embedding(spec_jw95, spec_free_r, v3)
    else:
        raise ValueError(f"Unexpected sub_name: {sub_name}")

    m_sub = sub.factory(v_sub)
    ll_sub = float(diagnostic_log_likelihood(m_sub, dsets, n_params=sub.n_free, anchor_po_max=True).total_log_likelihood)

    v_emb = nesting_embedding(sub, big, v_sub)
    m_emb = big.factory(v_emb)
    ll_emb = float(diagnostic_log_likelihood(m_emb, dsets, n_params=big.n_free, anchor_po_max=True).total_log_likelihood)

    # 1. Embedding identity: same physical model evaluates to the same log-likelihood
    assert ll_emb == pytest.approx(ll_sub, abs=1e-3, rel=1e-4)

    # 2. Optimization starting from embedded sub-model: enlarged model cannot fit worse
    fit_big = fit_diagnostic_consistency(
        big.factory, dsets, param_names=big.param_names, bounds=big.bounds,
        starts=[v_emb], anchor_po_max=True, n_random_starts=0,
    )
    # Optimizer tolerance tol=1e-4: the start vector IS v_emb, so Nelder-Mead / L-BFGS-B
    # can only improve upon or maintain ll_emb up to numerical precision.
    assert fit_big.total_log_likelihood >= ll_sub - 1e-4
