"""Comprehensive bug checks and unit tests for next-generation GABA-A models and protocols."""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from circuitpharm.models import (
    ReceptorModel,
    WaveformResult,
    OperationalScalarModel,
    KineticAllosteryModel,
    ExtendedDesensitizationModel,
)
from circuitpharm.protocols import (
    synaptic_transient,
    pulse_train,
    ambient_with_spillover,
    extract_electrophys_metrics,
    compare_phasic_tonic_mechanisms,
    evaluate_model_fit,
    find_discriminating_protocol,
)
from circuitpharm.fitting import (
    compute_kinetic_objective,
    compute_profile_likelihood,
    compute_fisher_information_matrix,
    equilibrium_dose_response_chi2,
    run_ensemble_mcmc,
)
from circuitpharm.results import Tier, VoidQuantityError
from circuitpharm.dynamic_range import evaluate_dynamic_range
from circuitpharm.gabaa_kinetics import Scheme


# ==============================================================================
# 1. Model Protocol & Individual Implementations
# ==============================================================================

def test_operational_scalar_model_protocol_and_behavior():
    model = OperationalScalarModel(ec50_um=20.0, hill_n=1.5, po_max=0.75, s_max=2.5)
    assert isinstance(model, ReceptorModel)
    assert model.name == "OperationalScalar"
    assert "ec50_um" in model.param_names

    # Zero agonist -> zero open prob
    assert model.steady_state(0.0) == 0.0

    # At EC50 -> exactly half of po_max
    p_half = model.steady_state(20.0, pam_factor=1.0)
    assert abs(p_half - 0.375) < 1e-4

    # PAM leftward shift capped at s_max (2.5x)
    p_pam_2x = model.steady_state(10.0, pam_factor=2.0)
    p_pam_5x = model.steady_state(10.0, pam_factor=5.0)  # should cap at 2.5x
    p_pam_2_5x = model.steady_state(10.0, pam_factor=2.5)
    assert abs(p_pam_5x - p_pam_2_5x) < 1e-9
    assert p_pam_2x < p_pam_2_5x

    # Dose response monotonic
    concs = np.logspace(-1, 3, 20)
    dr = model.dose_response(concs)
    assert np.all(np.diff(dr) >= 0.0)
    assert dr[-1] <= 0.75

    # Waveform simulation
    t = np.linspace(0.0, 50.0, 100)
    gaba_t = np.where((t >= 5.0) & (t <= 6.0), 1000.0, 0.0)
    res = model.simulate_waveform(t, gaba_t)
    assert isinstance(res, WaveformResult)
    assert res.peak_p_open > 0.5
    assert res.charge_integral > 0.0
    assert np.all(res.p_open >= 0.0) and np.all(res.p_open <= 1.0)


def test_kinetic_jw95_model_protocol_and_behavior():
    model = KineticAllosteryModel()
    assert isinstance(model, ReceptorModel)
    assert model.name == "KineticAllostery_JW95"

    # Microscopic constants
    assert abs(model.kd_um - (model.koff / model.kon)) < 1e-9
    assert abs(model.gating_efficacy - (model.beta / model.alpha)) < 1e-9
    assert abs(model.desens_ratio - (model.d / model.r)) < 1e-9

    # State distribution simplex constraint
    dist = model.state_distribution(0.40)
    assert len(dist) == 5
    assert abs(np.sum(dist) - 1.0) < 1e-6
    assert np.all(dist >= 0.0)

    # Steady state matches open probability in distribution
    po = model.steady_state(0.40)
    assert abs(po - dist[3]) < 1e-9

    # Infinite agonist limit matches po_inf
    po_high = model.steady_state(1e7)
    assert abs(po_high - model.po_inf) < 1e-4

    # Scheme compatibility wrapper
    s = Scheme(kon=model.kon, koff=model.koff, beta=model.beta, alpha=model.alpha, d=model.d, r=model.r)
    model_from_s = s.to_model()
    assert abs(model_from_s.steady_state(0.40) - po) < 1e-9

    # Dynamic waveform integration
    t = np.linspace(0.0, 30.0, 100)
    gaba = synaptic_transient(t, peak_um=1000.0, rise_ms=0.1, clear_ms=1.0)
    res = model.simulate_waveform(t, gaba)
    assert res.peak_p_open > 0.4
    assert res.states is not None
    assert res.states.shape == (100, 5)
    # Check probability conservation at all time points
    sums = np.sum(res.states, axis=1)
    assert np.all(np.abs(sums - 1.0) < 1e-4)


def test_kinetic_jw95_a2d_connectivity_is_strictly_to_a2r():
    """Assert that in Model B (JW95), A2D connects strictly to A2R, not to A2O.

    States: 0: R, 1: AR, 2: A2R, 3: A2O, 4: A2D.
    Desensitisation enters from doubly-bound closed A2R (rate d) and resensitises
    exclusively to A2R (rate r). It does NOT connect to or from open state A2O.
    """
    model = KineticAllosteryModel()
    q = model.q_matrix(gaba_um=100.0)

    # Row 4 (A2D outgoing transitions): only non-diagonal entry is to column 2 (A2R)
    for col in range(5):
        if col == 4:
            assert q[4, col] == -model.r, f"A2D diagonal should be -r, got {q[4, col]}"
        elif col == 2:
            assert q[4, col] == model.r, f"A2D -> A2R rate should be r, got {q[4, col]}"
        else:
            assert q[4, col] == 0.0, f"A2D transition to state {col} must be 0, got {q[4, col]}"

    # Column 4 (incoming transitions into A2D): only non-diagonal entry is from row 2 (A2R)
    for row in range(5):
        if row == 4:
            continue
        elif row == 2:
            assert q[row, 4] == model.d, f"A2R -> A2D rate should be d, got {q[row, 4]}"
        else:
            assert q[row, 4] == 0.0, f"Transition from state {row} to A2D must be 0, got {q[row, 4]}"

    # Explicit negative assertions: no direct connection with open state A2O (index 3)
    assert q[4, 3] == 0.0, "A2D must not transition to A2O (open state)"
    assert q[3, 4] == 0.0, "A2O (open state) must not transition to A2D"


def test_kinetic_open_desens_connectivity_is_strictly_to_a2o():
    """Assert that in KineticOpenDesensModel, A2D connects strictly to A2O, not to A2R.

    States: 0: R, 1: AR, 2: A2R, 3: A2O, 4: A2D.
    Desensitisation enters from open state A2O (rate d) and resensitises
    exclusively to A2O (rate r). It does NOT connect to or from closed state A2R.
    """
    from circuitpharm.models.kinetic_open_desens import KineticOpenDesensModel

    model = KineticOpenDesensModel()
    q = model.q_matrix(gaba_um=100.0)

    # All rows sum to zero
    assert np.allclose(q.sum(axis=1), 0.0, atol=1e-12)

    # Row 4 (A2D outgoing transitions): only non-diagonal entry is to column 3 (A2O)
    for col in range(5):
        if col == 4:
            assert q[4, col] == -model.r, f"A2D diagonal should be -r, got {q[4, col]}"
        elif col == 3:
            assert q[4, col] == model.r, f"A2D -> A2O rate should be r, got {q[4, col]}"
        else:
            assert q[4, col] == 0.0, f"A2D transition to state {col} must be 0, got {q[4, col]}"

    # Column 4 (incoming transitions into A2D): only non-diagonal entry is from row 3 (A2O)
    for row in range(5):
        if row == 4:
            continue
        elif row == 3:
            assert q[row, 4] == model.d, f"A2O -> A2D rate should be d, got {q[row, 4]}"
        else:
            assert q[row, 4] == 0.0, f"Transition from state {row} to A2D must be 0, got {q[row, 4]}"

    # Explicit negative assertions: no direct connection with closed state A2R (index 2)
    assert q[4, 2] == 0.0, "A2D must not transition to A2R (closed state)"
    assert q[2, 4] == 0.0, "A2R (closed state) must not transition to A2D"

    # Equilibrium state distribution sums to 1.0
    for gaba in [0.0, 1.0, 10.0, 100.0, 10000.0]:
        dist = model.state_distribution(gaba)
        assert np.isclose(np.sum(dist), 1.0, atol=1e-12)
        assert np.all(dist >= 0.0)


def test_kinetic_flip_connectivity_and_zeros():
    """Assert strictly verified connectivity for KineticFlipModel with explicit negative assertions.

    Topology:
        R (0) <-> AR (1) <-> A2R (2) <-> A2F (3) <-> A2O (4)
                               |
                             A2D (5)

    Explicit negative assertions required by Agent Brief 8B Task 1:
      1. A2F (state 3) connects strictly to A2R (2) and A2O (4), and nowhere else.
      2. A2D (state 5) connects strictly to A2R (2), and nowhere else.
      3. A2O (state 4) connects strictly to A2F (3), and nowhere else.
    Assert all zeros, not just the non-zeros.
    """
    from circuitpharm.models.kinetic_flip import KineticFlipModel, OPEN_IDX_6, STATES_6

    assert OPEN_IDX_6 == 4
    assert STATES_6 == ("R", "AR", "A2R", "A2F", "A2O", "A2D")

    model = KineticFlipModel()
    q = model.q_matrix(gaba_um=100.0)

    # 1. Generator rows sum to zero
    assert np.allclose(q.sum(axis=1), 0.0, atol=1e-12)

    # 2. A2D (state 5): connects only to A2R (state 2)
    # Row 5 (outgoing transitions from A2D):
    for col in range(6):
        if col == 5:
            assert q[5, col] == -model.r, f"A2D diagonal should be -r, got {q[5, col]}"
        elif col == 2:
            assert q[5, col] == model.r, f"A2D -> A2R rate should be r, got {q[5, col]}"
        else:
            assert q[5, col] == 0.0, f"A2D transition to state {col} ({STATES_6[col]}) must be 0, got {q[5, col]}"

    # Column 5 (incoming transitions into A2D):
    for row in range(6):
        if row == 5:
            continue
        elif row == 2:
            assert q[row, 5] == model.d, f"A2R -> A2D rate should be d, got {q[row, 5]}"
        else:
            assert q[row, 5] == 0.0, f"Transition from state {row} ({STATES_6[row]}) to A2D must be 0, got {q[row, 5]}"

    # Explicit negative assertions for A2D:
    assert q[5, 4] == 0.0, "A2D must not transition to A2O (open state)"
    assert q[4, 5] == 0.0, "A2O must not transition to A2D"
    assert q[5, 3] == 0.0, "A2D must not transition to A2F (flip state)"
    assert q[3, 5] == 0.0, "A2F must not transition to A2D"
    assert q[5, 1] == 0.0, "A2D must not transition to AR"
    assert q[1, 5] == 0.0, "AR must not transition to A2D"
    assert q[5, 0] == 0.0, "A2D must not transition to R"
    assert q[0, 5] == 0.0, "R must not transition to A2D"

    # 3. A2F (state 3): connects only to A2R (state 2) and A2O (state 4)
    # Row 3 (outgoing transitions from A2F):
    for col in range(6):
        if col == 3:
            assert q[3, col] == -(model.gamma + model.beta), f"A2F diagonal should be -(gamma+beta), got {q[3, col]}"
        elif col == 2:
            assert q[3, col] == model.gamma, f"A2F -> A2R rate should be gamma, got {q[3, col]}"
        elif col == 4:
            assert q[3, col] == model.beta, f"A2F -> A2O rate should be beta, got {q[3, col]}"
        else:
            assert q[3, col] == 0.0, f"A2F transition to state {col} ({STATES_6[col]}) must be 0, got {q[3, col]}"

    # Column 3 (incoming transitions into A2F):
    for row in range(6):
        if row == 3:
            continue
        elif row == 2:
            assert q[row, 3] == model.delta, f"A2R -> A2F rate should be delta, got {q[row, 3]}"
        elif row == 4:
            assert q[row, 3] == model.alpha, f"A2O -> A2F rate should be alpha, got {q[row, 3]}"
        else:
            assert q[row, 3] == 0.0, f"Transition from state {row} ({STATES_6[row]}) to A2F must be 0, got {q[row, 3]}"

    # Explicit negative assertions for A2F:
    assert q[3, 0] == 0.0, "A2F must not transition to R"
    assert q[0, 3] == 0.0, "R must not transition to A2F"
    assert q[3, 1] == 0.0, "A2F must not transition to AR"
    assert q[1, 3] == 0.0, "AR must not transition to A2F"
    assert q[3, 5] == 0.0, "A2F must not transition to A2D"
    assert q[5, 3] == 0.0, "A2D must not transition to A2F"

    # 4. A2O (state 4): connects strictly to A2F (state 3)
    # Row 4 (outgoing transitions from A2O):
    for col in range(6):
        if col == 4:
            assert q[4, col] == -model.alpha, f"A2O diagonal should be -alpha, got {q[4, col]}"
        elif col == 3:
            assert q[4, col] == model.alpha, f"A2O -> A2F rate should be alpha, got {q[4, col]}"
        else:
            assert q[4, col] == 0.0, f"A2O transition to state {col} ({STATES_6[col]}) must be 0, got {q[4, col]}"

    # Column 4 (incoming transitions into A2O):
    for row in range(6):
        if row == 4:
            continue
        elif row == 3:
            assert q[row, 4] == model.beta, f"A2F -> A2O rate should be beta, got {q[row, 4]}"
        else:
            assert q[row, 4] == 0.0, f"Transition from state {row} ({STATES_6[row]}) to A2O must be 0, got {q[row, 4]}"

    # Explicit negative assertions for A2O:
    assert q[4, 0] == 0.0, "A2O must not transition to R"
    assert q[0, 4] == 0.0, "R must not transition to A2O"
    assert q[4, 1] == 0.0, "A2O must not transition to AR"
    assert q[1, 4] == 0.0, "AR must not transition to A2O"
    assert q[4, 2] == 0.0, "A2O must not transition to A2R"
    assert q[2, 4] == 0.0, "A2R must not transition to A2O directly (must pass through A2F)"
    assert q[4, 5] == 0.0, "A2O must not transition to A2D"
    assert q[5, 4] == 0.0, "A2D must not transition to A2O"

    # 5. Detailed balance check: stationary distribution P @ Q == 0
    for gaba in [0.0, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0]:
        dist = model.state_distribution(gaba)
        assert np.isclose(np.sum(dist), 1.0, atol=1e-12)
        assert np.all(dist >= 0.0)
        q_g = model.q_matrix(gaba)
        assert np.allclose(dist @ q_g, 0.0, atol=1e-12)


def test_kinetic_flip_model_protocol_and_behavior():
    """Verify protocol conformance, po_max, steady state, and waveform simulation for KineticFlipModel."""
    from circuitpharm.models.kinetic_flip import KineticFlipModel, OPEN_IDX_6

    model = KineticFlipModel()
    assert isinstance(model, ReceptorModel)
    assert model.name == "KineticFlip"
    assert model.open_index == OPEN_IDX_6 == 4
    assert len(model.param_names) == 8

    # Zero agonist -> open prob 0 at resting
    p0 = model.steady_state(0.0)
    assert np.isclose(p0, 0.0, atol=1e-12)

    # State distribution simplex constraint (6 states)
    dist = model.state_distribution(10.0)
    assert len(dist) == 6
    assert np.isclose(np.sum(dist), 1.0, atol=1e-12)
    assert np.all(dist >= 0.0)

    # High agonist steady state converges to po_inf
    po_high = model.steady_state(1e8)
    assert np.isclose(po_high, model.po_inf, atol=1e-4)

    # po_max formula check
    f = model.flipping_ratio
    e = model.gating_efficacy
    expected_po_max = (f * e) / (1.0 + f * (1.0 + e))
    assert np.isclose(model.po_max, expected_po_max)

    # Transient waveform simulation
    t = np.linspace(0.0, 50.0, 100)
    gaba = synaptic_transient(t, peak_um=2000.0, rise_ms=0.1, clear_ms=1.0)
    res = model.simulate_waveform(t, gaba)
    assert res.peak_p_open > 0.0
    assert res.states is not None
    assert res.states.shape == (100, 6)
    assert np.all(np.abs(np.sum(res.states, axis=1) - 1.0) < 1e-4)


def test_kinetic_two_open_connectivity_and_zeros():
    """Assert strictly verified connectivity for KineticTwoOpenModel with explicit negative assertions.

    Topology:
        R (0) <-> AR (1) <-> A2R (2) <-> A2O (3) <-> A2O' (4)
                               |
                             A2D (5)

    Explicit negative assertions required by Agent Brief 9 Task 3:
      1. A2O' (state 4) connects strictly to A2O (3), and nowhere else.
      2. A2D (state 5) connects strictly to A2R (2), and nowhere else.
      3. A2O (state 3) connects strictly to A2R (2) and A2O' (4), and nowhere else.
    Assert all zeros, not just the non-zeros.
    """
    from circuitpharm.models.kinetic_two_open import KineticTwoOpenModel, OPEN_INDICES_6, STATES_6

    assert OPEN_INDICES_6 == (3, 4)
    assert STATES_6 == ("R", "AR", "A2R", "A2O", "A2O'", "A2D")

    model = KineticTwoOpenModel()
    q = model.q_matrix(gaba_um=100.0)

    # 1. Generator rows sum to zero
    assert np.allclose(q.sum(axis=1), 0.0, atol=1e-12)

    # 2. A2D (state 5): connects only to A2R (state 2)
    # Row 5 (outgoing transitions from A2D):
    for col in range(6):
        if col == 5:
            assert q[5, col] == -model.r, f"A2D diagonal should be -r, got {q[5, col]}"
        elif col == 2:
            assert q[5, col] == model.r, f"A2D -> A2R rate should be r, got {q[5, col]}"
        else:
            assert q[5, col] == 0.0, f"A2D transition to state {col} ({STATES_6[col]}) must be 0, got {q[5, col]}"

    # Column 5 (incoming transitions into A2D):
    for row in range(6):
        if row == 5:
            continue
        elif row == 2:
            assert q[row, 5] == model.d, f"A2R -> A2D rate should be d, got {q[row, 5]}"
        else:
            assert q[row, 5] == 0.0, f"Transition from state {row} ({STATES_6[row]}) to A2D must be 0, got {q[row, 5]}"

    # Explicit negative assertions for A2D:
    assert q[5, 4] == 0.0, "A2D must not transition to A2O' (second open state)"
    assert q[4, 5] == 0.0, "A2O' must not transition to A2D"
    assert q[5, 3] == 0.0, "A2D must not transition to A2O (first open state)"
    assert q[3, 5] == 0.0, "A2O must not transition to A2D"
    assert q[5, 1] == 0.0, "A2D must not transition to AR"
    assert q[1, 5] == 0.0, "AR must not transition to A2D"
    assert q[5, 0] == 0.0, "A2D must not transition to R"
    assert q[0, 5] == 0.0, "R must not transition to A2D"

    # 3. A2O' (state 4): connects strictly to A2O (state 3)
    # Row 4 (outgoing transitions from A2O'):
    for col in range(6):
        if col == 4:
            assert q[4, col] == -model.gamma, f"A2O' diagonal should be -gamma, got {q[4, col]}"
        elif col == 3:
            assert q[4, col] == model.gamma, f"A2O' -> A2O rate should be gamma, got {q[4, col]}"
        else:
            assert q[4, col] == 0.0, f"A2O' transition to state {col} ({STATES_6[col]}) must be 0, got {q[4, col]}"

    # Column 4 (incoming transitions into A2O'):
    for row in range(6):
        if row == 4:
            continue
        elif row == 3:
            assert q[row, 4] == model.delta, f"A2O -> A2O' rate should be delta, got {q[row, 4]}"
        else:
            assert q[row, 4] == 0.0, f"Transition from state {row} ({STATES_6[row]}) to A2O' must be 0, got {q[row, 4]}"

    # Explicit negative assertions for A2O':
    assert q[4, 0] == 0.0, "A2O' must not transition to R"
    assert q[0, 4] == 0.0, "R must not transition to A2O'"
    assert q[4, 1] == 0.0, "A2O' must not transition to AR"
    assert q[1, 4] == 0.0, "AR must not transition to A2O'"
    assert q[4, 2] == 0.0, "A2O' must not transition to A2R directly"
    assert q[2, 4] == 0.0, "A2R must not transition to A2O' directly"
    assert q[4, 5] == 0.0, "A2O' must not transition to A2D"
    assert q[5, 4] == 0.0, "A2D must not transition to A2O'"

    # 4. A2O (state 3): connects strictly to A2R (state 2) and A2O' (state 4)
    # Row 3 (outgoing transitions from A2O):
    for col in range(6):
        if col == 3:
            assert q[3, col] == -(model.alpha + model.delta), f"A2O diagonal should be -(alpha+delta), got {q[3, col]}"
        elif col == 2:
            assert q[3, col] == model.alpha, f"A2O -> A2R rate should be alpha, got {q[3, col]}"
        elif col == 4:
            assert q[3, col] == model.delta, f"A2O -> A2O' rate should be delta, got {q[3, col]}"
        else:
            assert q[3, col] == 0.0, f"A2O transition to state {col} ({STATES_6[col]}) must be 0, got {q[3, col]}"

    # Column 3 (incoming transitions into A2O):
    for row in range(6):
        if row == 3:
            continue
        elif row == 2:
            assert q[row, 3] == model.beta, f"A2R -> A2O rate should be beta, got {q[row, 3]}"
        elif row == 4:
            assert q[row, 3] == model.gamma, f"A2O' -> A2O rate should be gamma, got {q[row, 3]}"
        else:
            assert q[row, 3] == 0.0, f"Transition from state {row} ({STATES_6[row]}) to A2O must be 0, got {q[row, 3]}"

    # Explicit negative assertions for A2O:
    assert q[3, 0] == 0.0, "A2O must not transition to R"
    assert q[0, 3] == 0.0, "R must not transition to A2O"
    assert q[3, 1] == 0.0, "A2O must not transition to AR"
    assert q[1, 3] == 0.0, "AR must not transition to A2O"
    assert q[3, 5] == 0.0, "A2O must not transition to A2D"
    assert q[5, 3] == 0.0, "A2D must not transition to A2O"

    # 5. Detailed balance check: stationary distribution P @ Q == 0
    for gaba in [0.0, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0]:
        dist = model.state_distribution(gaba)
        assert np.isclose(np.sum(dist), 1.0, atol=1e-12)
        assert np.all(dist >= 0.0)
        q_g = model.q_matrix(gaba)
        assert np.allclose(dist @ q_g, 0.0, atol=1e-12)


def test_kinetic_two_open_model_protocol_and_behavior():
    """Verify protocol conformance, po_max, steady state, and waveform simulation for KineticTwoOpenModel."""
    from circuitpharm.models.kinetic_two_open import KineticTwoOpenModel, OPEN_INDICES_6

    model = KineticTwoOpenModel()
    assert isinstance(model, ReceptorModel)
    assert model.name == "KineticTwoOpen"
    assert model.open_indices == OPEN_INDICES_6 == (3, 4)
    assert len(model.param_names) == 8

    # Zero agonist -> open prob 0 at resting
    p0 = model.steady_state(0.0)
    assert np.isclose(p0, 0.0, atol=1e-12)

    # State distribution simplex constraint (6 states)
    dist = model.state_distribution(10.0)
    assert len(dist) == 6
    assert np.isclose(np.sum(dist), 1.0, atol=1e-12)
    assert np.all(dist >= 0.0)

    # Steady state is the sum of both open states
    ss = model.steady_state(10.0)
    assert np.isclose(ss, dist[3] + dist[4])

    # High agonist steady state converges to po_inf
    po_high = model.steady_state(1e8)
    assert np.isclose(po_high, model.po_inf, atol=1e-4)

    # po_max formula check
    e = model.gating_efficacy
    k2 = model.two_open_ratio
    eff_e = e * (1.0 + k2)
    expected_po_max = eff_e / (1.0 + eff_e)
    assert np.isclose(model.po_max, expected_po_max)

    # Transient waveform simulation: p_open must be sum of states 3 and 4
    t = np.linspace(0.0, 50.0, 100)
    gaba = synaptic_transient(t, peak_um=2000.0, rise_ms=0.1, clear_ms=1.0)
    res = model.simulate_waveform(t, gaba)
    assert res.peak_p_open > 0.0
    assert res.states is not None
    assert res.states.shape == (100, 6)
    assert np.all(np.abs(np.sum(res.states, axis=1) - 1.0) < 1e-4)
    assert np.allclose(res.p_open, res.states[:, 3] + res.states[:, 4], atol=1e-12)


def test_extended_desensitization_model_protocol_and_behavior():
    model = ExtendedDesensitizationModel()
    assert isinstance(model, ReceptorModel)
    assert model.name == "ExtendedDesensitization_DualD"

    # Dual-pathway ratios
    assert abs(model.total_desens_ratio - (model.desens_fast_ratio + model.desens_slow_ratio)) < 1e-9

    # State distribution simplex constraint (6 states)
    dist = model.state_distribution(1.0)
    assert len(dist) == 6
    assert abs(np.sum(dist) - 1.0) < 1e-6
    assert np.all(dist >= 0.0)

    # Infinite agonist limit
    po_high = model.steady_state(1e7)
    assert abs(po_high - model.po_inf) < 1e-4

    # Transient waveform simulation
    t = np.linspace(0.0, 50.0, 100)
    gaba = synaptic_transient(t, peak_um=2000.0, rise_ms=0.1, clear_ms=1.0)
    res = model.simulate_waveform(t, gaba)
    assert res.peak_p_open > 0.3
    assert res.states is not None
    assert res.states.shape == (100, 6)
    assert np.all(np.abs(np.sum(res.states, axis=1) - 1.0) < 1e-4)


# ==============================================================================
# 2. Waveforms & Protocols
# ==============================================================================

def test_waveform_generators():
    t = np.linspace(0.0, 100.0, 1000)

    # Synaptic transient
    tr = synaptic_transient(t, peak_um=1000.0, rise_ms=0.2, clear_ms=1.5, t_start=5.0)
    assert tr[0] == 0.0
    assert abs(np.max(tr) - 1000.0) < 5.0
    assert tr[-1] < 1.0

    # Pulse train
    train = pulse_train(t, freq_hz=50.0, n_pulses=3, peak_um=800.0, clear_ms=1.0, start_ms=10.0)
    # Peak should approach 800 uM
    assert np.max(train) >= 700.0
    # Between pulses it should drop
    assert np.min(train[t > 15.0]) < 100.0

    # Ambient spillover
    amb = ambient_with_spillover(t, ambient_um=0.40, spillover_peak_um=3.0, event_times_ms=[20.0, 60.0])
    assert amb[0] == 0.40
    assert np.max(amb) > 3.0


def test_electrophys_metrics_extraction():
    # 400 samples, not 100: the biexponential decay fit (P0-4) needs at least 6 points
    # inside the 90%->10% window, and a 1 ms clearance on a 50 ms span at 100 samples
    # does not supply them. A coarse grid now yields an honest NaN rather than 15.0.
    t = np.linspace(0.0, 50.0, 400)
    gaba = synaptic_transient(t, peak_um=1000.0, rise_ms=0.1, clear_ms=1.0)
    model = KineticAllosteryModel()
    res = model.simulate_waveform(t, gaba)
    # ONE driving-force convention, from models.base (P0-11); there is no
    # `driving_force_mv` argument any more.
    metrics = extract_electrophys_metrics(res, g_max_ns=2.0, v_hold_mv=-60.0,
                                          e_cl_mv=-75.0)

    assert metrics["peak_current_pA"] > 0.0
    assert metrics["charge_integral_fC"] > 0.0
    assert metrics["driving_force_mv"] == pytest.approx(15.0)
    assert metrics["charge_window_ms"] == pytest.approx(50.0)
    # decay may legitimately not fit on a given grid; if it does, it must be positive
    assert np.isnan(metrics["decay_tau_ms"]) or metrics["decay_tau_ms"] > 0.0


# ==============================================================================
# 3. Mechanistic Decomposition
# ==============================================================================

def test_mechanistic_decomposition_attribution():
    """Mechanics. The identity and the mechanism findings live in test_decomposition.py.

    The three factors are now occupancy / double-occupancy / gating and they satisfy
    ln G = sum exactly; the old occupancy / gating / HEADROOM triple did not, because the
    proximity-to-ceiling term is not a factor of the gain (roadmap §5.9).
    """
    model = KineticAllosteryModel()
    decomp = compare_phasic_tonic_mechanisms(model, tonic_gaba_um=0.40,
                                             phasic_gaba_um=1000.0, pam_factor=2.50)

    tonic = decomp["tonic"]
    phasic = decomp["phasic"]

    for r in (tonic, phasic):
        total = (r.log_share_pct_occupancy + r.log_share_pct_double_occupancy
                 + r.log_share_pct_gating)
        assert abs(total - 100.0) < 1e-3
        assert abs(r.residual) < 1e-10 * max(abs(r.log_gain), 1.0)

    assert tonic.total_fold_gain > phasic.total_fold_gain
    assert tonic.total_fold_gain > 2.0
    assert phasic.total_fold_gain < 1.5


# ==============================================================================
# 4. Model Comparison & Optimal Experimental Design
# ==============================================================================

def test_model_comparison_and_oed():
    """Mechanics. Every CONCLUSION from this module is VOID until P5/P6.

    Also note the explicit `observable=`: Model A is PEAK-native and Models B and C are
    EQUILIBRIUM-native, so there is no default axis and omitting it now raises.
    """
    from circuitpharm.models.base import Observable, ObservableMismatch

    m_a = OperationalScalarModel()
    m_b = KineticAllosteryModel()
    m_c = ExtendedDesensitizationModel()

    concs = np.array([0.1, 1.0, 10.0, 30.0, 100.0, 1000.0])
    y_true = m_b.dose_response(concs)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        comp = evaluate_model_fit([m_a, m_b, m_c], concs, y_true,
                                  measurement_noise_std=0.02,
                                  observable=Observable.PEAK)
        assert len(comp) == 3
        total_w = sum(c.akaike_weight.get(acknowledge_void=True) for c in comp)
        assert abs(total_w - 1.0) < 1e-4
        for c in comp:
            assert c.aic.tier is Tier.VOID

        protocol = find_discriminating_protocol([m_a, m_b, m_c], noise_sigma=0.02,
                                                observable=Observable.PEAK)
    assert protocol.discrimination_score.get(acknowledge_void=True) > 0.0
    assert len(protocol.predicted_responses) == 3
    assert len(protocol.falsification_boundaries.get(acknowledge_void=True)) == 3

    # OperationalScalarModel is PEAK-native; passing EQUILIBRIUM must raise ObservableMismatch
    with pytest.raises(ObservableMismatch, match="equilibrium"):
        evaluate_model_fit([m_a, m_b, m_c], concs, y_true,
                           measurement_noise_std=0.02,
                           observable=Observable.EQUILIBRIUM)
    with pytest.raises(ObservableMismatch, match="equilibrium"):
        find_discriminating_protocol([m_a, m_b, m_c], noise_sigma=0.02,
                                    observable=Observable.EQUILIBRIUM)


# ==============================================================================
# 5. Parameter Estimation & Identifiability
# ==============================================================================

def test_fitting_objective_and_identifiability():
    """Mechanics only. The SCIENCE of this layer is pinned in tests/test_provisional.py.

    This test used to assert `cost >= 0.0` and `np.isfinite(cost)` on a chi-squared of
    9290, which is why nothing noticed (roadmap C9). The value assertions now live in
    test_provisional.py; what remains here is that the shapes and types are right.
    """
    model = KineticAllosteryModel()
    x0 = np.array([model.kon, model.koff, model.beta, model.alpha, model.d, model.r])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # deprecated alias still works; the real name says what it constrains
        cost = compute_kinetic_objective(x0)
        assert cost == pytest.approx(equilibrium_dose_response_chi2(x0))

        prof = compute_profile_likelihood("kon", model, scan_factors=(0.8, 1.0, 1.2))
        assert prof.param_name == "kon"
        assert len(prof.profile_costs) == 3
        # the CI is a VOID Quantity now, not a tuple -- reading it must raise
        with pytest.raises(VoidQuantityError):
            prof.ci_95_bounds.value
        lo, hi = prof.ci_95_bounds.get(acknowledge_void=True)
        assert lo <= hi

        # the old FIM name is deprecated; it returns the Hessian and a VOID condition number
        fim, cond = compute_fisher_information_matrix(model)
        assert fim.shape == (6, 6)
        assert cond.tier is Tier.VOID
        assert np.isfinite(cond.get(acknowledge_void=True))


def test_lightweight_ensemble_mcmc():
    model = KineticAllosteryModel()
    x0 = np.array([model.kon, model.koff, model.beta, model.alpha, model.d, model.r])

    # Short MCMC run to verify sampling mechanics
    res = run_ensemble_mcmc(x0, n_walkers=14, n_steps=15, burn_in=5, seed=123)
    assert isinstance(res.flat_samples, np.ndarray)
    assert res.flat_samples.shape[1] == 6
    assert len(res.posterior_means) == 6
    assert 0.0 <= res.acceptance_fraction <= 1.0


# ==============================================================================
# 6. Three-Tiered Dynamic Range Hierarchy
# ==============================================================================

def test_dynamic_range_hierarchy():
    model = KineticAllosteryModel()
    eval_dr = evaluate_dynamic_range(model, ambient_gaba_um=0.40, cleft_peak_um=1000.0)

    # Tier 1: Asymptotic headroom must be large (> 100x)
    assert eval_dr.theoretical_asymptotic_headroom > 100.0

    # Tier 2: Reachable gains at finite s_max
    assert 2.50 in eval_dr.reachable_gains_by_smax
    assert eval_dr.reachable_gains_by_smax[2.50] > 1.5
    # Order monotonic with s_max
    assert eval_dr.reachable_gains_by_smax[1.25] < eval_dr.reachable_gains_by_smax[2.50]

    # Tier 3: Physiological charge ratio under waveforms
    assert eval_dr.physiological_charge_ratio > 1.0  # tonic advantage preserved
    assert not eval_dr.headroom_collapses
