"""P2/P4-3: the mechanistic decomposition, and the identity it must satisfy.

The previous version asserted an identity its code did not satisfy -- the third factor was
a proximity-to-ceiling ratio, which is not a factor of the gain -- and formed its
percentages from absolute values, discarding sign. See the module docstring of
`protocols/decomposition.py` for the replacement and roadmap §5.9.
"""
from __future__ import annotations

import numpy as np
import pytest

from circuitpharm.models import KineticAllosteryModel
from circuitpharm.protocols.decomposition import (
    FactorialDecompositionResult,
    _factors,
    compare_phasic_tonic_mechanisms,
    decompose_modulation_gain,
    decomposition_report,
)
from circuitpharm.results import Quantity, Tier

MODEL = KineticAllosteryModel()


def test_decomposition_identity_is_exact():
    """ANCHOR. ln G equals the sum of the three terms, to rtol 1e-10.

    Over 20 random parameter draws and both regimes. The tolerance is this tight because
    the decomposition is an IDENTITY -- Po = Occ x (P_A2/Occ) x (Po/P_A2) by construction --
    so the only admissible error is float rounding. If this ever fails, the state partition
    is wrong; do not widen it.
    """
    rng = np.random.default_rng(20261009)
    for _ in range(20):
        m = KineticAllosteryModel(
            kon=float(10 ** rng.uniform(-3, -1.5)),
            koff=float(10 ** rng.uniform(-1.5, 0.5)),
            beta=float(10 ** rng.uniform(-1, 0.5)),
            alpha=float(10 ** rng.uniform(-1.5, 0.0)),
            d=float(10 ** rng.uniform(-2, -0.5)),
            r=float(10 ** rng.uniform(-4, -2)),
        )
        for gaba in (0.1, 0.4, 10.0, 1000.0):
            for pf in (1.0, 1.5, 2.5, 8.0):
                r = decompose_modulation_gain(m, gaba, pam_factor=pf)
                total = (r.delta_log_occupancy + r.delta_log_double_occupancy
                         + r.delta_log_gating)
                assert total == pytest.approx(r.log_gain, rel=1e-10, abs=1e-12)
                assert abs(r.residual) < 1e-10 * max(abs(r.log_gain), 1.0)


def test_the_shares_sum_to_one_hundred_percent():
    """Signed shares of an exactly additive decomposition must sum to 100%."""
    for gaba in (0.4, 1000.0):
        r = decompose_modulation_gain(MODEL, gaba, pam_factor=2.5)
        total = (r.log_share_pct_occupancy + r.log_share_pct_double_occupancy
                 + r.log_share_pct_gating)
        assert total == pytest.approx(100.0, abs=1e-9)


def test_no_modulation_gives_unit_gain_and_zero_shares():
    """The old code returned 33.33/33.33/33.34 here, inventing structure from nothing."""
    r = decompose_modulation_gain(MODEL, 0.4, pam_factor=1.0)
    assert r.total_fold_gain == pytest.approx(1.0, rel=1e-12)
    assert r.log_share_pct_occupancy == 0.0
    assert r.log_share_pct_double_occupancy == 0.0
    assert r.log_share_pct_gating == 0.0


def test_the_factors_reconstruct_the_open_probability():
    """ANCHOR. The product of the three factors IS Po, exactly."""
    for gaba in (0.05, 0.4, 20.0, 3000.0):
        dist = MODEL.state_distribution(gaba)
        occ, dbl, gate = _factors(dist)
        assert occ * dbl * gate == pytest.approx(float(dist[3]), rel=1e-12)


def test_an_affinity_modulator_leaves_the_gating_share_at_its_ceiling():
    """ANCHOR, and the finding that overturns the roadmap's expectation.

    An affinity-type PAM divides k_off and touches neither beta/alpha nor d/r, so the
    equilibrium distribution WITHIN the doubly-liganded states cannot move:
    P(open | doubly bound) is pinned at E/(1+E+D) regardless of agonist concentration or
    modulator strength. So there is no "operating point headroom" to spend at equilibrium,
    and the whole phasic/tonic divergence is binding-side.

    Exact to 1e-12: this is algebra, not an approximation.
    """
    ceiling = MODEL.gating_efficacy / (1.0 + MODEL.gating_efficacy + MODEL.desens_ratio)
    for gaba in (0.1, 0.4, 10.0, 1000.0):
        for pf in (1.0, 2.5, 100.0):
            r = decompose_modulation_gain(MODEL, gaba, pam_factor=pf)
            assert r.delta_log_gating == pytest.approx(0.0, abs=1e-12)
            assert r.gating_share_baseline == pytest.approx(ceiling, rel=1e-12)
            assert r.gating_share_ceiling == pytest.approx(ceiling, rel=1e-12)


def test_a_gating_modulator_moves_the_gating_share_and_an_affinity_one_does_not():
    """ANCHOR. The decomposition is a mechanism assay: the two PAM types are separable.

    Measured at 0.4 uM: a 2.5x gating modulator roughly doubles the gating share
    (0.166169 -> 0.332536, tracking its own raised ceiling) while barely moving occupancy
    (3.166e-2 -> 3.298e-2). A 2.5x affinity modulator does the opposite -- occupancy
    triples and the gating share does not move at all.
    """
    base = _factors(MODEL.state_distribution(0.4))
    aff = _factors(MODEL.apply_pam(affinity_factor=2.5).state_distribution(0.4))
    gat = _factors(MODEL.apply_pam(gating_factor=2.5).state_distribution(0.4))

    # affinity: occupancy moves a lot, gating not at all
    assert aff[0] / base[0] > 2.5
    assert aff[2] == pytest.approx(base[2], rel=1e-12)
    # gating: gating share moves a lot, occupancy barely
    assert gat[2] / base[2] == pytest.approx(2.0, rel=0.05)
    assert gat[0] / base[0] < 1.1


def test_the_tonic_gain_exceeds_the_phasic_gain():
    """The project's central qualitative claim, on this instrument.

    Asserted as an inequality with a wide margin rather than as a value: the magnitudes
    depend on the fitted rates (UNCALIBRATED), the ordering does not.
    """
    both = compare_phasic_tonic_mechanisms(MODEL, pam_factor=2.5)
    assert both["tonic"].total_fold_gain > 3.0
    assert both["phasic"].total_fold_gain < 1.1
    assert both["tonic"].total_fold_gain > 4.0 * both["phasic"].total_fold_gain


def test_phasic_gain_is_dominated_by_the_double_occupancy_term():
    """At a saturating transient, occupancy has nowhere to go, so term 2 carries the gain."""
    r = decompose_modulation_gain(MODEL, 1000.0, pam_factor=2.5)
    assert r.log_share_pct_double_occupancy > 90.0
    assert r.log_share_pct_occupancy < 10.0


def test_an_explicit_model_is_required():
    with pytest.raises(TypeError, match="requires an explicit"):
        decompose_modulation_gain(None, 0.4)
    with pytest.raises(TypeError, match="requires an explicit"):
        compare_phasic_tonic_mechanisms(None)


def test_negative_modulation_is_out_of_domain():
    with pytest.raises(ValueError, match="below 1.0"):
        decompose_modulation_gain(MODEL, 0.4, pam_factor=0.5)


def test_the_identity_check_raises_rather_than_returning_a_bad_result():
    """`check()` is called inside the constructor path, so a broken partition cannot ship."""
    bad = FactorialDecompositionResult(
        condition="x", gaba_conc_um=1.0, pam_factor=2.0, total_fold_gain=2.0,
        log_gain=np.log(2.0), delta_log_occupancy=0.1, delta_log_double_occupancy=0.1,
        delta_log_gating=0.1, log_share_pct_occupancy=33.3,
        log_share_pct_double_occupancy=33.3, log_share_pct_gating=33.3,
        gating_share_baseline=0.1, gating_share_ceiling=0.2,
        residual=np.log(2.0) - 0.3)
    with pytest.raises(AssertionError, match="identity failed"):
        bad.check()


def test_the_report_is_tiered():
    """P2: no bare quoteable floats out of this module."""
    rs = decomposition_report(MODEL)
    assert rs.items
    for q in rs.items:
        assert isinstance(q, Quantity)
        assert q.tier is Tier.UNCALIBRATED, (
            f"{q.name} is {q.tier.label}; the identity is exact but the magnitudes rest on "
            f"recalled anchors plus a convention, so nothing here is VALIDATED")
        assert q.provenance and q.promote_by
    names = {q.name for q in rs.items}
    assert "tonic_fold_gain" in names and "phasic_share_gating" in names
