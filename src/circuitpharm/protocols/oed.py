"""Statistical Model Comparison and Optimal Experimental Design (OED).

Compares Models A (Scalar), B (5-State Kinetic), and C (Extended Dual-Desensitisation),
evaluates AIC/BIC out-of-sample goodness-of-fit, and discovers experimental stimulation
protocols that maximally discriminate between candidate mechanisms.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence
import numpy as np

from ..models.base import ReceptorModel


@dataclass(frozen=True)
class ModelComparisonResult:
    """Statistical model selection metrics across candidate models."""
    model_name: str
    n_params: int
    n_data: int
    rss: float
    log_likelihood: float
    aic: float
    bic: float
    akaike_weight: float


@dataclass(frozen=True)
class DiscriminatingProtocol:
    """An experimentally realizable protocol that maximally separates candidate models."""
    gaba_um: float
    pam_factor: float
    pulse_duration_ms: float
    predicted_responses: dict[str, float]
    falsification_boundaries: dict[str, tuple[float, float]]
    discrimination_score: float  # Normalized separation between predictions


def evaluate_model_fit(
    models: Sequence[ReceptorModel],
    concs_um: np.ndarray,
    observed_responses: np.ndarray,
    measurement_noise_std: float = 0.03,
) -> list[ModelComparisonResult]:
    """Calculate log-likelihood, AIC, BIC, and Akaike weights across candidate models.

    Parameters
    ----------
    models : Sequence[ReceptorModel]
        List of competing models (e.g. Model A, Model B, Model C).
    concs_um : np.ndarray
        Agonist concentrations of the test data.
    observed_responses : np.ndarray
        Measured open probabilities or normalized current amplitudes.
    measurement_noise_std : float
        Estimated standard deviation of measurement noise.
    """
    y_obs = np.asarray(observed_responses, dtype=float)
    concs = np.asarray(concs_um, dtype=float)
    n_data = len(y_obs)
    sigma = max(measurement_noise_std, 1e-6)

    raw_results = []
    aics = []

    for m in models:
        y_pred = m.dose_response(concs, pam_factor=1.0)
        residuals = y_obs - y_pred
        rss = float(np.sum(residuals ** 2))
        
        # Gaussian log-likelihood
        log_lik = float(-0.5 * np.sum((residuals / sigma) ** 2) - n_data * np.log(sigma * np.sqrt(2.0 * np.pi)))
        k = len(m.param_names)
        aic = 2.0 * k - 2.0 * log_lik
        bic = k * np.log(max(n_data, 1)) - 2.0 * log_lik
        
        aics.append(aic)
        raw_results.append((m.name, k, rss, log_lik, aic, bic))

    # Akaike weights: exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
    min_aic = min(aics)
    delta_aics = np.array([a - min_aic for a in aics])
    w = np.exp(-0.5 * delta_aics)
    weights = w / np.sum(w)

    final_results = []
    for (name, k, rss, log_lik, aic, bic), weight in zip(raw_results, weights):
        final_results.append(ModelComparisonResult(
            model_name=name,
            n_params=k,
            n_data=n_data,
            rss=rss,
            log_likelihood=log_lik,
            aic=aic,
            bic=bic,
            akaike_weight=float(weight),
        ))

    return final_results


def find_discriminating_protocol(
    models: Sequence[ReceptorModel],
    gaba_candidates_um: Sequence[float] | None = None,
    pam_candidates: Sequence[float] | None = None,
    noise_sigma: float = 0.02,
) -> DiscriminatingProtocol:
    """Find the GABA and PAM exposure combination that maximizes separation between models.

    Parameters
    ----------
    models : Sequence[ReceptorModel]
        The candidate models to discriminate between.
    gaba_candidates_um : Sequence[float], optional
        Candidate GABA concentrations to test (defaults to 0.1 to 1000 uM sweep).
    pam_candidates : Sequence[float], optional
        Candidate PAM multipliers (defaults to 1.5, 2.5, 4.0).
    noise_sigma : float
        Typical experimental patch-clamp noise level.
    """
    gaba_grid = gaba_candidates_um or [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 1000.0]
    pam_grid = pam_candidates or [1.5, 2.0, 2.5, 3.0, 4.0]

    best_score = -1.0
    best_gaba = 1.0
    best_pam = 2.5
    best_preds: dict[str, float] = {}

    for g in gaba_grid:
        for p in pam_grid:
            preds = {m.name: m.steady_state(g, pam_factor=p) for m in models}
            vals = list(preds.values())
            # Spread = max pairwise distance normalized by noise
            max_dist = max(vals) - min(vals)
            score = max_dist / max(noise_sigma, 1e-4)

            if score > best_score:
                best_score = score
                best_gaba = g
                best_pam = p
                best_preds = preds

    # Generate 95% confidence tolerance intervals around each prediction (2 * sigma)
    boundaries = {
        name: (max(0.0, val - 2.0 * noise_sigma), min(1.0, val + 2.0 * noise_sigma))
        for name, val in best_preds.items()
    }

    return DiscriminatingProtocol(
        gaba_um=best_gaba,
        pam_factor=best_pam,
        pulse_duration_ms=100.0,
        predicted_responses=best_preds,
        falsification_boundaries=boundaries,
        discrimination_score=float(best_score),
    )
