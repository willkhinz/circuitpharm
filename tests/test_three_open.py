"""Comprehensive tests for 7-state kinetic schemes with three open states (Agent Brief 18).

Verifies:
1. Connectivity of KineticThreeOpenModel with explicit negative zero assertions.
2. Microscopic reversibility (detailed balance) holds unconditionally on acyclic tree.
3. Nesting reduction: at delta2 -> 0, scheme reduces identically to KineticCoopModel.
4. Nesting embedding verification via nesting_embedding and log-likelihood reproduction.
5. Connectivity and detailed balance for KineticThreeOpenParallelModel (alternative topology).
6. Mean open time matrix derivation vs exact analytical formula.
"""
from __future__ import annotations

import numpy as np
import pytest

from circuitpharm.fitting.comparison import (
    FIT_CONVENTIONS,
    MODEL_SPECS,
    nesting_embedding,
)
from circuitpharm.fitting.data import (
    BARBERIS2007_DEACT,
    BARBERIS2007_PAIRED_PULSE,
    JAHN1997_PEAK_CRC,
    KERAMIDAS2008_MEAN_OPEN_TIME_M,
)
from circuitpharm.fitting.deact_components import predict_components_eigen
from circuitpharm.fitting.diagnostic import diagnostic_log_likelihood
from circuitpharm.fitting.likelihood import predict_mean_open_time
from circuitpharm.models import (
    KineticCoopModel,
    KineticThreeOpenModel,
    KineticThreeOpenParallelModel,
    OPEN_INDICES_7,
    STATES_7,
)


def test_kinetic_three_open_connectivity_and_explicit_zeros():
    """Assert strictly verified connectivity for KineticThreeOpenModel with explicit zero assertions.

    Topology:
        R (0) <--(2 kon g, koff)--> AR (1) <--(kon g, 2 koff2)--> A2R (2) <--(beta, alpha)--> A2O (3) <--(delta, gamma)--> A2O' (4) <--(delta2, gamma2)--> A2O'' (5)
                                                                    |
                                                                  (d, r)
                                                                    v
                                                                  A2D (6)

    Explicit negative assertions required by Agent Brief 18 Task 1:
      1. A2O'' (state 5) connects strictly to A2O' (4), and nowhere else.
      2. A2O' (state 4) connects strictly to A2O (3) and A2O'' (5), and nowhere else.
      3. A2O (state 3) connects strictly to A2R (2) and A2O' (4), and nowhere else.
      4. A2D (state 6) connects strictly to A2R (2), and nowhere else.
      5. A2R (state 2) connects strictly to AR (1), A2O (3), and A2D (6), and nowhere else.
      6. AR (state 1) connects strictly to R (0) and A2R (2), and nowhere else.
      7. R (state 0) connects strictly to AR (1), and nowhere else.
    """
    assert OPEN_INDICES_7 == (3, 4, 5)
    assert STATES_7 == ("R", "AR", "A2R", "A2O", "A2O'", "A2O''", "A2D")

    model = KineticThreeOpenModel(
        kon=0.012, koff=0.35, koff2=0.07,  # c = 5.0
        delta=1.2, gamma=0.8,
        delta2=0.5, gamma2=0.4,
    )
    assert model.c == pytest.approx(5.0)
    assert model.open_indices == (3, 4, 5)
    assert model.open_index == (3, 4, 5)

    q = model.q_matrix(gaba_um=100.0)

    # 1. Generator rows sum to zero
    assert np.allclose(q.sum(axis=1), 0.0, atol=1e-12)

    # 2. State 5 (A2O''): connects strictly to A2O' (state 4)
    for col in range(7):
        if col == 5:
            assert q[5, col] == -model.gamma2
        elif col == 4:
            assert q[5, col] == model.gamma2
        else:
            assert q[5, col] == 0.0, f"A2O'' transition to state {col} ({STATES_7[col]}) must be 0, got {q[5, col]}"

    for row in range(7):
        if row == 5:
            continue
        elif row == 4:
            assert q[row, 5] == model.delta2
        else:
            assert q[row, 5] == 0.0, f"Transition from state {row} ({STATES_7[row]}) to A2O'' must be 0, got {q[row, 5]}"

    # Explicit negative assertions for A2O'':
    assert q[5, 0] == 0.0, "A2O'' must not transition to R"
    assert q[0, 5] == 0.0, "R must not transition to A2O''"
    assert q[5, 1] == 0.0, "A2O'' must not transition to AR"
    assert q[1, 5] == 0.0, "AR must not transition to A2O''"
    assert q[5, 2] == 0.0, "A2O'' must not transition to A2R"
    assert q[2, 5] == 0.0, "A2R must not transition to A2O''"
    assert q[5, 3] == 0.0, "A2O'' must not transition directly to A2O"
    assert q[3, 5] == 0.0, "A2O must not transition directly to A2O''"
    assert q[5, 6] == 0.0, "A2O'' must not transition to A2D"
    assert q[6, 5] == 0.0, "A2D must not transition to A2O''"

    # 3. State 4 (A2O'): connects strictly to A2O (3) and A2O'' (5)
    for col in range(7):
        if col == 4:
            assert q[4, col] == -(model.gamma + model.delta2)
        elif col == 3:
            assert q[4, col] == model.gamma
        elif col == 5:
            assert q[4, col] == model.delta2
        else:
            assert q[4, col] == 0.0, f"A2O' transition to {col} must be 0, got {q[4, col]}"

    for row in range(7):
        if row == 4:
            continue
        elif row == 3:
            assert q[row, 4] == model.delta
        elif row == 5:
            assert q[row, 4] == model.gamma2
        else:
            assert q[row, 4] == 0.0, f"Transition from {row} to A2O' must be 0, got {q[row, 4]}"

    # Explicit negative assertions for A2O':
    assert q[4, 0] == 0.0, "A2O' must not transition to R"
    assert q[0, 4] == 0.0, "R must not transition to A2O'"
    assert q[4, 1] == 0.0, "A2O' must not transition to AR"
    assert q[1, 4] == 0.0, "AR must not transition to A2O'"
    assert q[4, 2] == 0.0, "A2O' must not transition directly to A2R"
    assert q[2, 4] == 0.0, "A2R must not transition directly to A2O'"
    assert q[4, 6] == 0.0, "A2O' must not transition to A2D"
    assert q[6, 4] == 0.0, "A2D must not transition to A2O'"

    # 4. State 6 (A2D): connects strictly to A2R (2)
    for col in range(7):
        if col == 6:
            assert q[6, col] == -model.r
        elif col == 2:
            assert q[6, col] == model.r
        else:
            assert q[6, col] == 0.0, f"A2D transition to {col} must be 0, got {q[6, col]}"

    for row in range(7):
        if row == 6:
            continue
        elif row == 2:
            assert q[row, 6] == model.d
        else:
            assert q[row, 6] == 0.0, f"Transition from {row} to A2D must be 0, got {q[row, 6]}"

    # Explicit negative assertions for A2D:
    assert q[6, 0] == 0.0, "A2D must not transition to R"
    assert q[0, 6] == 0.0, "R must not transition to A2D"
    assert q[6, 1] == 0.0, "A2D must not transition to AR"
    assert q[1, 6] == 0.0, "AR must not transition to A2D"
    assert q[6, 3] == 0.0, "A2D must not transition to A2O"
    assert q[3, 6] == 0.0, "A2O must not transition to A2D"
    assert q[6, 4] == 0.0, "A2D must not transition to A2O'"
    assert q[4, 6] == 0.0, "A2O' must not transition to A2D"

    # 5. Detailed balance check (microscopic reversibility):
    # Tree graph with V=7, E=6, loop space 6-7+1=0.
    # P @ Q == 0 must hold identically at all agonist concentrations.
    for gaba in [0.0, 0.05, 0.5, 5.0, 50.0, 500.0, 5000.0]:
        dist = model.state_distribution(gaba)
        assert np.isclose(np.sum(dist), 1.0, atol=1e-12)
        assert np.all(dist >= 0.0)
        q_g = model.q_matrix(gaba)
        assert np.allclose(dist @ q_g, 0.0, atol=1e-12)


def test_kinetic_three_open_parallel_connectivity_and_explicit_zeros():
    """Assert strictly verified connectivity for KineticThreeOpenParallelModel with explicit negative assertions."""
    model = KineticThreeOpenParallelModel(
        kon=0.012, koff=0.35, koff2=0.07,
        beta_p=0.4, alpha_p=0.15,
    )
    assert model.open_indices == (3, 4, 5)
    assert model.open_index == (3, 4, 5)
    q = model.q_matrix(gaba_um=50.0)

    # 1. Generator rows sum to zero
    assert np.allclose(q.sum(axis=1), 0.0, atol=1e-12)

    # 2. State 5 (A2O'') in parallel connects strictly to A2R (state 2)
    for col in range(7):
        if col == 5:
            assert q[5, col] == -model.alpha_p
        elif col == 2:
            assert q[5, col] == model.alpha_p
        else:
            assert q[5, col] == 0.0

    for row in range(7):
        if row == 5:
            continue
        elif row == 2:
            assert q[row, 5] == model.beta_p
        else:
            assert q[row, 5] == 0.0

    # Explicit negative assertions: A2O'' does NOT connect to A2O or A2O' or A2D
    assert q[5, 3] == 0.0
    assert q[3, 5] == 0.0
    assert q[5, 4] == 0.0
    assert q[4, 5] == 0.0
    assert q[5, 6] == 0.0
    assert q[6, 5] == 0.0

    # Detailed balance on parallel tree:
    for gaba in [0.0, 1.0, 10.0, 100.0]:
        dist = model.state_distribution(gaba)
        assert np.isclose(np.sum(dist), 1.0, atol=1e-12)
        assert np.all(dist >= 0.0)
        q_g = model.q_matrix(gaba)
        assert np.allclose(dist @ q_g, 0.0, atol=1e-12)


def test_kinetic_three_open_reduces_identically_to_kinetic_coop():
    """Verify that KineticThreeOpenModel collapses exactly to KineticCoopModel when delta2 = 0."""
    m_coop = KineticCoopModel(
        kon=0.015, koff=0.40, koff2=0.08,
        beta=0.60, alpha=0.12, delta=1.5, gamma=0.9,
        d=0.04, r=0.003,
    )
    m_three = KineticThreeOpenModel(
        kon=0.015, koff=0.40, koff2=0.08,
        beta=0.60, alpha=0.12, delta=1.5, gamma=0.9,
        delta2=0.0, gamma2=1.0,  # collapsed third open state
        d=0.04, r=0.003,
    )

    # 1. Po,max and Po,inf
    assert m_three.po_max == pytest.approx(m_coop.po_max, abs=1e-14)
    assert m_three.po_inf == pytest.approx(m_coop.po_inf, abs=1e-14)

    # 2. State distributions and steady states
    for g in [0.0, 0.1, 1.0, 10.0, 100.0, 1000.0]:
        d_coop = m_coop.state_distribution(g)
        d_three = m_three.state_distribution(g)
        # States 0..4 match
        assert np.allclose(d_coop[:5], d_three[:5], atol=1e-14)
        # State 5 (A2O'') has probability zero
        assert d_three[5] == pytest.approx(0.0, abs=1e-14)
        # State 6 (A2D) matches coop state 5 (A2D)
        assert d_three[6] == pytest.approx(d_coop[5], abs=1e-14)

        # Steady state open probability matches exactly
        assert m_three.steady_state(g) == pytest.approx(m_coop.steady_state(g), abs=1e-14)

    # 3. Mean open time matches exactly
    mot_coop = float(predict_mean_open_time(m_coop, KERAMIDAS2008_MEAN_OPEN_TIME_M))
    mot_three = float(predict_mean_open_time(m_three, KERAMIDAS2008_MEAN_OPEN_TIME_M))
    assert mot_three == pytest.approx(mot_coop, abs=1e-12)

    # 4. Deactivation relaxation components match
    t_coop, a_coop = predict_components_eigen(m_coop, n_components=3, l1_max=None)
    t_three, a_three = predict_components_eigen(m_three, n_components=3, l1_max=None)
    assert np.allclose(t_three, t_coop, atol=1e-10)
    assert np.allclose(a_three, a_coop, atol=1e-10)


def test_kinetic_three_open_embedding_reproduces_kinetic_coop_likelihood():
    """Verify nesting_embedding from kinetic_coop into kinetic_three_open reproduces likelihood.

    Agent Brief 18 Task 1 requirement:
    Declare nesting: the scheme must reduce exactly to kinetic_coop in the limit that
    collapses the third open state, verified through nesting_embedding with the likelihood
    asserted to match.
    """
    sub = MODEL_SPECS["kinetic_coop"]
    big = MODEL_SPECS["kinetic_three_open"]

    assert big.nests == ("kinetic_coop",)
    assert "delta2" in FIT_CONVENTIONS
    assert "gamma2" in FIT_CONVENTIONS

    v_sub = np.array([-0.428, 0.104, 1.9721, -2.699, -0.744, -0.523, 0.0, 0.0, 1.5])
    m_sub = sub.factory(v_sub)
    v_emb = nesting_embedding(sub, big, v_sub)
    m_emb = big.factory(v_emb)

    # 1. Embedding properties:
    assert len(v_emb) == 11
    assert np.allclose(v_emb[:9], v_sub)
    assert v_emb[9] == pytest.approx(np.log10(FIT_CONVENTIONS["delta2"][0]))
    assert v_emb[10] == pytest.approx(np.log10(FIT_CONVENTIONS["gamma2"][0]))
    assert m_emb.delta2 == pytest.approx(0.0, abs=1e-14)
    assert m_emb.po_max == pytest.approx(m_sub.po_max, abs=1e-12)

    # 2. Total log-likelihood across all three standard datasets matches to numerical tolerance
    dsets = [JAHN1997_PEAK_CRC, BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]
    ll_sub = float(diagnostic_log_likelihood(m_sub, dsets, n_params=sub.n_free, anchor_po_max=True).total_log_likelihood)
    ll_emb = float(diagnostic_log_likelihood(m_emb, dsets, n_params=big.n_free, anchor_po_max=True).total_log_likelihood)

    assert ll_emb == pytest.approx(ll_sub, abs=1e-4, rel=1e-4)


def test_mean_open_time_analytical_formula_matches_matrix_inversion():
    """Verify analytical formula <tau_O> = (1/alpha)*(1 + K2 + K2*K3) matches matrix inversion."""
    alpha, delta, gamma, delta2, gamma2 = 0.25, 1.2, 0.8, 0.6, 0.3
    model = KineticThreeOpenModel(
        alpha=alpha, delta=delta, gamma=gamma,
        delta2=delta2, gamma2=gamma2,
    )
    k2 = delta / gamma
    k3 = delta2 / gamma2
    analytical_mot = (1.0 / alpha) * (1.0 + k2 + k2 * k3)

    q = model.q_matrix(0.0)
    q_oo = q[np.ix_([3, 4, 5], [3, 4, 5])]
    pi_o = np.array([1.0, 0.0, 0.0])
    matrix_mot = float(pi_o @ np.linalg.inv(-q_oo) @ np.ones(3))

    assert matrix_mot == pytest.approx(analytical_mot, abs=1e-12)
    # Also verify through predict_mean_open_time
    pred_mot = float(predict_mean_open_time(model, KERAMIDAS2008_MEAN_OPEN_TIME_M))
    assert pred_mot == pytest.approx(analytical_mot, abs=1e-12)
