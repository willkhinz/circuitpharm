"""P2: every quoteable number out of the new layer carries a tier.

The package's whole scientific identity is `results.Tier` and `provenance.Basis`
(`results.py`'s docstring records the reporting failure they exist to prevent), and the
3,177 lines added on 2026-10-08 bypassed both. This file is the standing check that they
do not drift back out.
"""
from __future__ import annotations

import inspect
import warnings

import numpy as np
import pytest

from circuitpharm import dynamic_range as dr
from circuitpharm.models import (
    ExtendedDesensitizationModel,
    KineticAllosteryModel,
    OperationalScalarModel,
)
from circuitpharm.models.base import Observable
from circuitpharm.protocols import decomposition, oed
from circuitpharm.results import Quantity, ResultSet, Tier, VoidQuantityError

MODEL = KineticAllosteryModel()
MODELS = [OperationalScalarModel(), KineticAllosteryModel(), ExtendedDesensitizationModel()]


def test_dynamic_range_report_is_tiered():
    rs = dr.dynamic_range_report(MODEL)
    assert isinstance(rs, ResultSet) and rs.items
    for q in rs.items:
        assert q.tier is Tier.UNCALIBRATED, f"{q.name} is {q.tier.label}"
        assert q.provenance and q.promote_by
        assert q.caveats, f"{q.name} has no caveats; every tier here is conditional"


def test_the_asymptotic_tier_says_it_is_unreachable():
    """ANCHOR on the wording, because this is the number that got quoted as a result.

    The asymptotic headroom is a k_off -> 0+ limit. Its caveats must say so, and must say
    it is conditional on d/r, which no anchor constrains.
    """
    rs = dr.dynamic_range_report(MODEL)
    q = rs.quantity("theoretical_asymptotic_headroom")
    joined = " ".join(q.caveats).lower()
    assert "no ligand" in joined or "unreachab" in joined
    assert "d/r" in joined


def test_the_charge_tier_carries_its_window():
    rs = dr.dynamic_range_report(MODEL)
    q = rs.quantity("physiological_charge_ratio")
    assert "ms" in q.provenance
    assert any("window" in c for c in q.caveats)
    # and the ObservedQuantity form refuses to be compared across windows
    ev = dr.evaluate_dynamic_range(MODEL)
    oq = dr.charge_quantity(ev)
    assert oq.observable is Observable.CHARGE
    assert oq.window_ms == ev.charge_window_ms


def test_the_collapse_boundary_is_a_curve_not_a_boolean():
    """ANCHOR. The scan must locate where the tonic advantage ends.

    Measured: the ratio falls monotonically with ambient GABA and crosses 1.0 between 3 and
    10 uM. A boolean `headroom_collapses` discarded all of that. Monotonicity is asserted
    rather than the crossing value, because the crossing moves with the fitted rates while
    the direction is structural -- more ambient agonist means less room.
    """
    ev = dr.evaluate_dynamic_range(MODEL)
    assert len(ev.collapse_scan) >= 6
    ratios = [x for _, x in ev.collapse_scan]
    assert all(b <= a + 1e-9 for a, b in zip(ratios, ratios[1:])), (
        f"the collapse scan is not monotone in ambient GABA: {ratios}")
    assert ratios[0] > 1.0 and ratios[-1] < 1.0, "the scan does not bracket the crossing"


def test_oed_metrics_are_void_until_fitted():
    """ANCHOR. Reading an AIC raises, and the message names the unfitted-parameter reason."""
    from circuitpharm.fitting import JAHN1997_PEAK_CRC as ds

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = oed.evaluate_model_fit(MODELS, ds.concs_um, ds.mean_response,
                                         observable=Observable.PEAK)
    assert results
    for r in results:
        for q in (r.aic, r.bic, r.akaike_weight):
            assert isinstance(q, Quantity) and q.tier is Tier.VOID
            with pytest.raises(VoidQuantityError) as exc:
                q.value
            assert "UNFITTED" in str(exc.value)
            assert np.isfinite(q.get(acknowledge_void=True))
        # the raw fit statistics stay plain and usable
        assert np.isfinite(r.rss) and np.isfinite(r.log_likelihood)
        assert r.observable is Observable.PEAK


def test_oed_refuses_to_compare_across_observables():
    """Model A is PEAK-native, Models B and C EQUILIBRIUM-native; one axis is required."""
    from circuitpharm.models.base import ObservableMismatch
    from circuitpharm.fitting import JAHN1997_PEAK_CRC as ds

    with pytest.raises(ObservableMismatch, match="different observables"):
        oed.evaluate_model_fit(MODELS, ds.concs_um, ds.mean_response)
    with pytest.raises(ObservableMismatch, match="different observables"):
        oed.find_discriminating_protocol(MODELS)


def test_the_discrimination_score_is_void():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p = oed.find_discriminating_protocol(MODELS, observable=Observable.PEAK)
    for q in (p.discrimination_score, p.falsification_boundaries):
        assert q.tier is Tier.VOID
        with pytest.raises(VoidQuantityError):
            q.value
    # the protocol itself is plain -- it is a choice of conditions, not an inference
    assert p.gaba_um > 0 and p.pam_factor >= 1.0
    assert set(p.predicted_responses) == {m.name for m in MODELS}
    # the defect text (shared by both quantities) is where the parameter-uncertainty
    # reason lives; the per-quantity caveat adds the multiple-comparison cost
    assert "parameter uncertainty" in p.discrimination_score.provenance
    assert any("multiple-comparison" in c for c in p.discrimination_score.caveats)


def test_void_quantities_are_withheld_from_printed_output():
    """`results.ResultSet.__str__` must not print a VOID value -- the original defect."""
    from circuitpharm.fitting import JAHN1997_PEAK_CRC as ds

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = oed.evaluate_model_fit(MODELS, ds.concs_um, ds.mean_response,
                                         observable=Observable.PEAK)
    rs = ResultSet("oed")
    for r in results:
        rs.add(r.aic)
    text = str(rs)
    assert "VOID" in text and "not quotable" in text
    for r in results:
        assert f"{r.aic.get(acknowledge_void=True):.4g}" not in text


# ------------------------------------------------------------------ the standing sweep
#: Public callables that legitimately return a bare number: pure algebra on a model's own
#: parameters, where there is nothing to attach evidence to. Everything else must return a
#: Quantity/ResultSet or a container whose inferential fields are Quantities.
_ALGEBRA_ALLOWLIST = {
    "kd_um", "gating_efficacy", "desens_ratio", "desens_fast_ratio",
    "desens_slow_ratio", "total_desens_ratio", "po_max", "po_inf", "effective_ec50",
    "steady_state", "name", "param_names", "native_observable", "get_params",
}


@pytest.mark.parametrize("model", MODELS, ids=lambda m: m.name)
def test_model_properties_that_return_bare_floats_are_algebra_only(model):
    """A derived model property returning a float must be algebra on its own parameters.

    Not a tier check -- a model's K_d is not a model OUTPUT, it is a restatement of its
    inputs. DATACLASS FIELDS are skipped for the same reason, only more so: kon and alpha
    are the inputs themselves. What is swept is the DERIVED properties, and the allowlist
    is deliberately short so a new float-returning property has to be justified by adding
    its name here.
    """
    fields = set(getattr(type(model), "__dataclass_fields__", {}))
    for name, attr in inspect.getmembers(type(model)):
        if name.startswith("_") or name in fields:
            continue
        if not isinstance(attr, property):
            continue
        try:
            val = getattr(model, name)
        except Exception:
            continue
        if isinstance(val, float):
            assert name in _ALGEBRA_ALLOWLIST, (
                f"{type(model).__name__}.{name} returns a bare float. If it is a model "
                f"OUTPUT it must carry a Tier (roadmap §2.1); if it is algebra on the "
                f"model's own parameters, add it to _ALGEBRA_ALLOWLIST with a reason.")


def test_the_report_entry_points_all_return_resultsets():
    """The three places a reader gets numbers from must all be tiered."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert isinstance(dr.dynamic_range_report(MODEL), ResultSet)
        assert isinstance(decomposition.decomposition_report(MODEL), ResultSet)
