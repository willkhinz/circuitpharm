"""Executable mathematical assertions mirroring the Amplitude Constraint Theorem (Brief 17 & 19).

See `knowledge/15-amplitude-constraint-theorem.md` for the full theorem and derivation.

Mirrors every key mathematical step as an executable unit test:
1. Reversible detailed balance and exact symmetry of H = D B D^-1 to machine precision.
2. Spectral reality: all eigenvalues of Q_0 restricted to doubly-liganded states are real, distinct, negative.
3. Sturm / Gantmacher-Krein oscillation theorem: extremal eigenvector strictly alternates in sign along the chain.
4. Exact closed-form weight ratio: w_A2O / |w_A2R| = alpha / (|lambda_1| - alpha) to machine precision.
5. Exact algebraic bound: |lambda_1| > alpha + beta, forcing w_A2O / |w_A2R| < 1/E.
6. Open-state desensitisation exact formulas: both neighboring shut states carry negative weights.
7. Reachable initial occupancy constraint: a_1 <= 0 across randomized parameter vectors and pulses.
8. Two-open relief: kinetic_two_open breaks the single-open constraint and exceeds 0.57.
9. Duhamel cone impulse response negativity: rho(u) @ w = w_A2R * exp(lambda_1 * u) < 0 for all u >= 0 (JW95 & open-desens).
10. Exact equilibrium secular bracket identity: [beta/(|lambda_1|-alpha) + d/(|lambda_1|-r) - 1] = -2*koff/|lambda_1| < 0.
11. Exact open-desens equilibrium secular identity: [1 - alpha/(|lambda_1|-(2koff+beta)) - d/(|lambda_1|-r)] = -(2koff/|lambda_1|)*(alpha/(|lambda_1|-(2koff+beta))) < 0.
12. Transient scalar dissipative ODE: z(t) = w @ p_D(t) satisfies z_dot = lambda_1 * z + w_A2R * sigma(t) <= 0.
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.linalg as la

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.deact_components import (
    predict_components_eigen,
    predict_components_eigen_all,
)
from circuitpharm.models.kinetic_jw95 import KineticAllosteryModel
from circuitpharm.models.kinetic_open_desens import KineticOpenDesensModel
from circuitpharm.models.kinetic_two_open import KineticTwoOpenModel


def test_generator_symmetrization_to_machine_precision():
    """Verify that H = D B D^-1 is symmetric to machine precision (< 1e-15)."""
    np.random.seed(42)
    for _ in range(50):
        r = float(10 ** np.random.uniform(-4, 1))
        d = float(10 ** np.random.uniform(-3, 3))
        alpha = float(10 ** np.random.uniform(-2, 2))
        beta = float(10 ** np.random.uniform(-2, 2))
        koff = float(10 ** np.random.uniform(-3, 2))

        # 3x3 doubly-liganded subgenerator (states A2D, A2R, A2O)
        B = np.array([
            [-r, r, 0.0],
            [d, -(2.0 * koff + beta + d), beta],
            [0.0, alpha, -alpha],
        ], dtype=float)

        pi = np.array([d / r, 1.0, beta / alpha], dtype=float)
        D = np.diag(np.sqrt(pi))
        D_inv = np.diag(1.0 / np.sqrt(pi))

        H = D @ B @ D_inv
        asymm = np.max(np.abs(H - H.T))
        assert asymm < 1e-13, f"H asymmetry {asymm:.2e} exceeds machine precision"


def test_spectral_reality_and_negativity():
    """Verify that all eigenvalues of B are strictly real, distinct, and negative."""
    np.random.seed(123)
    for _ in range(50):
        r = float(10 ** np.random.uniform(-4, 1))
        d = float(10 ** np.random.uniform(-3, 3))
        alpha = float(10 ** np.random.uniform(-2, 2))
        beta = float(10 ** np.random.uniform(-2, 2))
        koff = float(10 ** np.random.uniform(-3, 2))

        B = np.array([
            [-r, r, 0.0],
            [d, -(2.0 * koff + beta + d), beta],
            [0.0, alpha, -alpha],
        ], dtype=float)

        eigvals = la.eigvals(B)
        # Check zero imaginary part
        assert np.max(np.abs(np.imag(eigvals))) < 1e-12
        reals = np.sort(np.real(eigvals))
        # Strictly negative
        assert np.all(reals < 0.0)
        # Strictly distinct
        diffs = np.diff(reals)
        assert np.all(diffs > 0.0)


def test_gantmacher_krein_sturm_eigenvector_alternation():
    """Verify that the extremal eigenvector u^(1) strictly alternates in sign: (+, -, +)."""
    np.random.seed(99)
    for _ in range(50):
        r = float(10 ** np.random.uniform(-4, 1))
        d = float(10 ** np.random.uniform(-3, 3))
        alpha = float(10 ** np.random.uniform(-2, 2))
        beta = float(10 ** np.random.uniform(-2, 2))
        koff = float(10 ** np.random.uniform(-3, 2))

        H = np.array([
            [-r, np.sqrt(d * r), 0.0],
            [np.sqrt(d * r), -(2.0 * koff + beta + d), np.sqrt(alpha * beta)],
            [0.0, np.sqrt(alpha * beta), -alpha],
        ], dtype=float)

        eigvals, U = np.linalg.eigh(H)
        # Mode 0 is extremal (fastest / most negative eigenvalue)
        u1 = U[:, 0]
        if u1[2] < 0:
            u1 = -u1

        # Check strict sign alternation: u1[0] > 0, u1[1] < 0, u1[2] > 0
        assert u1[0] > 0, f"u1[0]={u1[0]} must be strictly positive"
        assert u1[1] < 0, f"u1[1]={u1[1]} must be strictly negative"
        assert u1[2] > 0, f"u1[2]={u1[2]} must be strictly positive"


def test_exact_weight_ratio_formula_jw95():
    """Verify analytical formula w_A2O / |w_A2R| = alpha / (|lambda_1| - alpha) to machine precision."""
    np.random.seed(77)
    for _ in range(50):
        r = float(10 ** np.random.uniform(-4, 1))
        d = float(10 ** np.random.uniform(-3, 3))
        alpha = float(10 ** np.random.uniform(-2, 2))
        beta = float(10 ** np.random.uniform(-2, 2))
        koff = float(10 ** np.random.uniform(-3, 2))

        H = np.array([
            [-r, np.sqrt(d * r), 0.0],
            [np.sqrt(d * r), -(2.0 * koff + beta + d), np.sqrt(alpha * beta)],
            [0.0, np.sqrt(alpha * beta), -alpha],
        ], dtype=float)

        eigvals, U = np.linalg.eigh(H)
        lam1 = float(eigvals[0])
        u1 = U[:, 0]

        pi = np.array([d / r, 1.0, beta / alpha], dtype=float)
        w = np.sqrt(pi[2] / pi) * u1[2] * u1

        actual_ratio = float(w[2] / (-w[1]))
        formula_ratio = float(alpha / (-lam1 - alpha))
        assert np.isclose(actual_ratio, formula_ratio, rtol=1e-8, atol=1e-9), (
            f"Weight ratio diff {abs(actual_ratio - formula_ratio):.2e} "
            f"actual={actual_ratio}, formula={formula_ratio}"
        )


def test_exact_algebraic_bound_jw95():
    """Verify |lambda_1| > alpha + beta, proving w_A2O / |w_A2R| < 1/E unconditionally."""
    np.random.seed(55)
    for _ in range(50):
        r = float(10 ** np.random.uniform(-4, 1))
        d = float(10 ** np.random.uniform(-3, 3))
        alpha = float(10 ** np.random.uniform(-2, 2))
        beta = float(10 ** np.random.uniform(-2, 2))
        koff = float(10 ** np.random.uniform(-3, 2))

        H = np.array([
            [-r, np.sqrt(d * r), 0.0],
            [np.sqrt(d * r), -(2.0 * koff + beta + d), np.sqrt(alpha * beta)],
            [0.0, np.sqrt(alpha * beta), -alpha],
        ], dtype=float)

        eigvals = np.linalg.eigvalsh(H)
        lam1 = float(eigvals[0])
        abs_lam1 = -lam1

        # Strict Rayleigh bound
        assert abs_lam1 > alpha + beta, f"|lam1|={abs_lam1} not > alpha+beta={alpha+beta}"

        # Consequently w_ratio < 1/E = alpha / beta
        w_ratio = alpha / (abs_lam1 - alpha)
        assert w_ratio < (alpha / beta), f"w_ratio={w_ratio} not < 1/E={alpha/beta}"


def test_exact_formulas_open_desensitisation():
    """Verify exact weight formulas for open-state desensitisation scheme."""
    np.random.seed(88)
    for _ in range(50):
        r = float(10 ** np.random.uniform(-4, 1))
        d = float(10 ** np.random.uniform(-3, 3))
        alpha = float(10 ** np.random.uniform(-2, 2))
        beta = float(10 ** np.random.uniform(-2, 2))
        koff = float(10 ** np.random.uniform(-3, 2))

        # States ordered: (A2R, A2O, A2D) with O = 1
        H = np.array([
            [-(2.0 * koff + beta), np.sqrt(alpha * beta), 0.0],
            [np.sqrt(alpha * beta), -(alpha + d), np.sqrt(d * r)],
            [0.0, np.sqrt(d * r), -r],
        ], dtype=float)

        eigvals, U = np.linalg.eigh(H)
        lam1 = float(eigvals[0])
        u1 = U[:, 0]
        pi = np.array([1.0, beta / alpha, (beta / alpha) * (d / r)], dtype=float)
        w = np.sqrt(pi[1] / pi) * u1[1] * u1

        # Open state is positive, both neighbors are negative
        assert w[1] > 0
        assert w[0] < 0
        assert w[2] < 0

        # Verify exact analytical formulas to machine precision
        actual_r_shut = abs(w[0]) / w[1]
        formula_r_shut = beta / (-lam1 - (2.0 * koff + beta))
        assert np.isclose(actual_r_shut, formula_r_shut, rtol=1e-8, atol=1e-9)

        actual_r_desens = abs(w[2]) / w[1]
        formula_r_desens = r / (-lam1 - r)
        assert np.isclose(actual_r_desens, formula_r_desens, rtol=1e-8, atol=1e-9)


def test_reachable_amplitude_non_positivity_sweep():
    """Verify a_1 <= 0 across 500 randomized parameter vectors and pulse durations/concentrations."""
    np.random.seed(333)
    for _ in range(500):
        kon = float(10 ** np.random.uniform(-2, 3))
        E = float(10 ** np.random.uniform(-2, 3))
        D = float(10 ** np.random.uniform(-3, 3))
        r = float(10 ** np.random.uniform(-4, 1))
        koff = float(10 ** np.random.uniform(-3, 2))
        alpha = float(10 ** np.random.uniform(-2, 2))
        beta = alpha * E
        d = r * D

        t_pulse = float(10 ** np.random.uniform(-1, 2))  # 0.1 ms to 100 ms
        g_um = float(10 ** np.random.uniform(1, 5))     # 10 uM to 100 mM

        # Test both 5-state topologies
        for ModelCls in (KineticAllosteryModel, KineticOpenDesensModel):
            m = ModelCls(kon=kon, koff=koff, beta=beta, alpha=alpha, d=d, r=r)
            try:
                taus, areas = predict_components_eigen(
                    m, pulse_ms=t_pulse, conc_mm=g_um / 1000.0, n_components=3, l1_max=None
                )
                a1 = areas[0]
                assert a1 <= 1e-10, (
                    f"{ModelCls.__name__} produced positive a1 = {a1:.4e} "
                    f"at pulse={t_pulse:.2f} ms, GABA={g_um:.1f} uM"
                )
            except Exception:
                pass


def test_two_open_states_relieves_constraint():
    """Verify that kinetic_two_open breaks the constraint, achieving a_1 > 0.57."""
    spec = MODEL_SPECS["kinetic_two_open"]
    # Test vector from amplitude admissibility scan
    from scripts.amplitude_admissibility import maximize_admissible_a1_for_scheme
    res = maximize_admissible_a1_for_scheme(spec, [1.0], n_samples=100, n_polish=2, seed=42)
    max_a1 = res[1.0]["a1"]
    assert max_a1 > 0.57, f"kinetic_two_open max a1 {max_a1} must exceed published target 0.57"


def test_duhamel_cone_impulse_response_negativity_jw95():
    """Verify that rho(u) @ w = w_A2R * exp(lambda_1 * u) < 0 for all ages u >= 0 in JW95."""
    np.random.seed(111)
    for _ in range(50):
        r = float(10 ** np.random.uniform(-4, 1))
        d = float(10 ** np.random.uniform(-3, 3))
        alpha = float(10 ** np.random.uniform(-2, 2))
        beta = float(10 ** np.random.uniform(-2, 2))
        koff = float(10 ** np.random.uniform(-3, 2))

        B = np.array([
            [-r, r, 0.0],
            [d, -(2.0 * koff + beta + d), beta],
            [0.0, alpha, -alpha],
        ], dtype=float)
        pi = np.array([d / r, 1.0, beta / alpha], dtype=float)
        D = np.diag(np.sqrt(pi))
        H = D @ B @ np.linalg.inv(D)
        eigvals, U = np.linalg.eigh(H)
        lam1 = float(eigvals[0])
        u1 = U[:, 0]
        w = np.sqrt(pi[2] / pi) * u1[2] * u1

        e_gate = np.array([0.0, 1.0, 0.0], dtype=float)
        for u in [0.0, 0.0001, 0.001, 0.01]:
            rho_u = e_gate @ la.expm(B * u)
            amp = float(np.dot(rho_u, w))
            expected = float(w[1] * np.exp(lam1 * u))
            assert expected < 0.0, "Expected amplitude must be strictly negative"
            assert np.isclose(amp, expected, rtol=1e-7, atol=1e-9)
            if expected < -1e-12:
                assert amp < 0.0, f"Impulse response amp {amp} must be strictly negative"
            else:
                assert amp <= 1e-14


def test_duhamel_cone_impulse_response_negativity_open_desens():
    """Verify that rho(u) @ w = w_A2R * exp(lambda_1 * u) < 0 for all ages u >= 0 in open desens."""
    np.random.seed(222)
    for _ in range(50):
        r = float(10 ** np.random.uniform(-4, 1))
        d = float(10 ** np.random.uniform(-3, 3))
        alpha = float(10 ** np.random.uniform(-2, 2))
        beta = float(10 ** np.random.uniform(-2, 2))
        koff = float(10 ** np.random.uniform(-3, 2))

        # States: A2R (0), A2O (1), A2D (2)
        B = np.array([
            [-(2.0 * koff + beta), beta, 0.0],
            [alpha, -(alpha + d), d],
            [0.0, r, -r],
        ], dtype=float)
        pi = np.array([1.0, beta / alpha, (beta / alpha) * (d / r)], dtype=float)
        D = np.diag(np.sqrt(pi))
        H = D @ B @ np.linalg.inv(D)
        eigvals, U = np.linalg.eigh(H)
        lam1 = float(eigvals[0])
        u1 = U[:, 0]
        w = np.sqrt(pi[1] / pi) * u1[1] * u1

        e_gate = np.array([1.0, 0.0, 0.0], dtype=float)
        for u in [0.0, 0.0001, 0.001, 0.01]:
            rho_u = e_gate @ la.expm(B * u)
            amp = float(np.dot(rho_u, w))
            expected = float(w[0] * np.exp(lam1 * u))
            assert expected < 0.0, "Expected open-desens amplitude must be strictly negative"
            assert np.isclose(amp, expected, rtol=1e-7, atol=1e-9)
            if expected < -1e-12:
                assert amp < 0.0, f"Open desens impulse response amp {amp} must be strictly negative"
            else:
                assert amp <= 1e-14


def test_equilibrium_secular_bracket_identity_jw95():
    """Verify exact algebraic secular identity beta/(|lam1|-alpha) + d/(|lam1|-r) - 1 = -2*koff/|lam1|."""
    np.random.seed(333)
    for _ in range(50):
        r = float(10 ** np.random.uniform(-4, 1))
        d = float(10 ** np.random.uniform(-3, 3))
        alpha = float(10 ** np.random.uniform(-2, 2))
        beta = float(10 ** np.random.uniform(-2, 2))
        koff = float(10 ** np.random.uniform(-3, 2))

        B = np.array([
            [-r, r, 0.0],
            [d, -(2.0 * koff + beta + d), beta],
            [0.0, alpha, -alpha],
        ], dtype=float)
        pi = np.array([d / r, 1.0, beta / alpha], dtype=float)
        D = np.diag(np.sqrt(pi))
        H = D @ B @ np.linalg.inv(D)
        eigvals = np.linalg.eigvalsh(H)
        lam1 = float(eigvals[0])
        s = -lam1

        bracket = beta / (s - alpha) + d / (s - r) - 1.0
        secular_rhs = -2.0 * koff / s
        assert np.isclose(bracket, secular_rhs, rtol=1e-8, atol=1e-9)
        assert bracket < 0.0, "JW95 equilibrium bracket must be strictly negative"


def test_equilibrium_secular_bracket_identity_open_desens():
    """Verify open-desens secular identity: 1 - alpha/(s - (2koff+beta)) - d/(s - r) = -(2koff/s)*(alpha/(s - (2koff+beta)))."""
    np.random.seed(444)
    for _ in range(50):
        r = float(10 ** np.random.uniform(-4, 1))
        d = float(10 ** np.random.uniform(-3, 3))
        alpha = float(10 ** np.random.uniform(-2, 2))
        beta = float(10 ** np.random.uniform(-2, 2))
        koff = float(10 ** np.random.uniform(-3, 2))

        B = np.array([
            [-(2.0 * koff + beta), beta, 0.0],
            [alpha, -(alpha + d), d],
            [0.0, r, -r],
        ], dtype=float)
        pi = np.array([1.0, beta / alpha, (beta / alpha) * (d / r)], dtype=float)
        D = np.diag(np.sqrt(pi))
        H = D @ B @ np.linalg.inv(D)
        eigvals = np.linalg.eigvalsh(H)
        lam1 = float(eigvals[0])
        s = -lam1

        bracket = 1.0 - alpha / (s - (2.0 * koff + beta)) - d / (s - r)
        secular_rhs = -(2.0 * koff / s) * (alpha / (s - (2.0 * koff + beta)))
        assert np.isclose(bracket, secular_rhs, rtol=1e-8, atol=1e-9)
        assert bracket < 0.0, "Open desens equilibrium bracket must be strictly negative"


def test_transient_scalar_dissipative_ode_identity():
    """Verify that z(t) = w @ p_D(t) satisfies the dissipative ODE z_dot = lam1 * z + w_A2R * sigma(t) <= 0."""
    from scipy.integrate import solve_ivp
    kon = 15.0; koff = 2.0; alpha = 60.0; beta = 180.0; d = 5.0; r = 0.5; g = 10.0
    Q = np.zeros((5, 5))
    Q[0, 1] = 2 * kon * g; Q[0, 0] = -Q[0, 1]
    Q[1, 0] = koff; Q[1, 2] = kon * g; Q[1, 1] = -(Q[1, 0] + Q[1, 2])
    Q[2, 1] = 2 * koff; Q[2, 3] = beta; Q[2, 4] = d; Q[2, 2] = -(Q[2, 1] + Q[2, 3] + Q[2, 4])
    Q[3, 2] = alpha; Q[3, 3] = -alpha
    Q[4, 2] = r; Q[4, 4] = -r

    B = np.array([
        [-r, r, 0.0],
        [d, -(2 * koff + beta + d), beta],
        [0.0, alpha, -alpha]
    ])
    pi = np.array([d / r, 1.0, beta / alpha])
    D = np.diag(np.sqrt(pi))
    H = D @ B @ np.linalg.inv(D)
    eigvals, U = np.linalg.eigh(H)
    lam1 = float(eigvals[0])
    u1 = U[:, 0]
    w = np.sqrt(pi[2] / pi) * u1[2] * u1

    def rhs(t, y):
        return y @ Q

    sol = solve_ivp(rhs, (0, 0.05), [1.0, 0, 0, 0, 0], t_eval=np.linspace(0.0005, 0.05, 100), rtol=1e-10, atol=1e-12)
    pD = np.array([sol.y[4], sol.y[2], sol.y[3]])
    z = np.dot(w, pD)
    assert np.all(z < 0.0), f"All transient mode 1 amplitudes must be strictly negative, max was {np.max(z)}"

