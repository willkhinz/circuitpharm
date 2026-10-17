"""Unit tests for KineticFlipModel, MODEL_SPECS integration, and flip structural search (Brief 8B)."""
from __future__ import annotations

import numpy as np
import pytest

from circuitpharm.fitting.comparison import FIT_CONVENTIONS, MODEL_SPECS
from circuitpharm.fitting.data import BARBERIS2007_PAIRED_PULSE, JAHN1997_PEAK_CRC
from circuitpharm.fitting.deact_components import (
    parse_deactivation_errors,
    predict_components_eigen_all,
    verify_eigen_reconstruction,
)
from circuitpharm.literature import DEACTIVATIONS
from circuitpharm.models import KineticFlipModel, ReceptorModel
from circuitpharm.models.base import Observable
from circuitpharm.models.kinetic_flip import OPEN_IDX_6, STATES_6
from scripts.flip_structural_search import (
    check_at_bound,
    evaluate_spec_vector,
    run_structural_search,
)


def test_kinetic_flip_model_instantiation_and_properties():
    """Verify state order, open index, rate accessors, and protocol properties."""
    model = KineticFlipModel()
    assert isinstance(model, ReceptorModel)
    assert model.name == "KineticFlip"
    assert model.native_observable == Observable.EQUILIBRIUM
    assert model.open_index == OPEN_IDX_6 == 4
    assert STATES_6[model.open_index] == "A2O"
    assert STATES_6 == ("R", "AR", "A2R", "A2F", "A2O", "A2D")

    params = model.get_params()
    assert len(params) == 8
    for p in ("kon", "koff", "delta", "gamma", "beta", "alpha", "d", "r"):
        assert p in params
        assert params[p] > 0.0

    # Dimensionless ratios
    assert np.isclose(model.kd_um, model.koff / model.kon)
    assert np.isclose(model.flipping_ratio, model.delta / model.gamma)
    assert np.isclose(model.gating_efficacy, model.beta / model.alpha)
    assert np.isclose(model.desens_ratio, model.d / model.r)

    # po_max formula: (F * E) / (1 + F * (1 + E))
    f = model.flipping_ratio
    e = model.gating_efficacy
    expected_po_max = (f * e) / (1.0 + f * (1.0 + e))
    assert np.isclose(model.po_max, expected_po_max)

    # In delta -> inf limit, po_max -> E / (1 + E) = beta / (beta + alpha)
    large_flip = model.with_params(delta=1e8, gamma=1.0)
    expected_jw95_po_max = large_flip.beta / (large_flip.beta + large_flip.alpha)
    assert np.isclose(large_flip.po_max, expected_jw95_po_max, atol=1e-5)


def test_kinetic_flip_eigen_reconstruction_agrees_with_simulation():
    """Verify that 5-component eigen-decomposition reproduces simulate_waveform to < 1e-6."""
    model = KineticFlipModel()
    dev = verify_eigen_reconstruction(model, open_index=model.open_index, atol=1e-6)
    assert dev < 1e-6


def test_kinetic_flip_spec_in_model_specs():
    """Verify MODEL_SPECS entry for kinetic_flip and its nesting exclusion."""
    assert "kinetic_flip" in MODEL_SPECS
    spec = MODEL_SPECS["kinetic_flip"]
    assert spec.name == "kinetic_flip"
    assert len(spec.param_names) == 8
    assert len(spec.bounds) == 8
    assert spec.conventions == ()
    assert spec.nests == ()  # excluded honestly per Brief 8B Task 2

    # Factory round-trip
    vec = np.array([0.0, 0.0, 1.0, -2.0, 0.0, 0.0, 1.0, 0.0])
    m = spec.factory(vec)
    assert isinstance(m, KineticFlipModel)
    assert np.isclose(m.kd_um, 1.0)
    assert np.isclose(m.flipping_ratio, 10.0)


def test_check_at_bound_detects_pinned_parameters():
    """Verify check_at_bound correctly identifies bound-pinned parameters."""
    bounds = ((-1.0, 4.0), (-2.0, 3.0))
    names = ("log10_kd", "log10_E")

    # Interior
    assert check_at_bound(np.array([1.0, 0.0]), bounds, names) == []

    # Lower bound on first
    pinned_l = check_at_bound(np.array([-0.9999, 0.0]), bounds, names)
    assert "log10_kd (L)" in pinned_l

    # Upper bound on second
    pinned_u = check_at_bound(np.array([1.0, 2.9999]), bounds, names)
    assert "log10_E (U)" in pinned_u


def test_structural_search_no_literals_and_runs_cleanly():
    """Verify that structural search runs cleanly with no literals and returns valid structure."""
    b = DEACTIVATIONS[1]
    target_taus = np.asarray(b.components[0], dtype=float)
    target_areas = np.asarray(b.components[1], dtype=float)
    tau_errs, area_errs = (np.array(x, dtype=float) for x in parse_deactivation_errors(b))

    spec = MODEL_SPECS["kinetic_flip"]
    res = evaluate_spec_vector(spec, np.zeros(8), target_taus, target_areas, tau_errs, area_errs)
    assert res["valid"] is True
    assert res["chi2"] > 0.0
    assert len(res["pred_taus"]) == 3
    assert len(res["pred_areas"]) == 3

    # Fast mini-search
    mini = run_structural_search("kinetic_flip", n_samples=20, n_polish=2, seed=1)
    assert mini["best_chi2"] > 0.0
    assert len(mini["best_vec"]) == 8
    assert mini["elapsed_s"] > 0.0
    assert np.isfinite(mini["pp_pred"])
    assert np.isfinite(mini["crc_rmse"])
