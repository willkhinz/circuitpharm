"""Lightweight pure NumPy/SciPy Affine-Invariant Ensemble MCMC Sampler.

Implements the Goodman & Weare (2010) stretch move to estimate joint posterior parameter
distributions and credible intervals without external PPL dependencies.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence
import numpy as np

from .identifiability import compute_kinetic_objective


@dataclass(frozen=True)
class MCMCChainResult:
    """Posterior sampling results from affine-invariant MCMC."""
    param_names: tuple[str, ...]
    flat_samples: np.ndarray          # shape: (n_samples, n_params)
    posterior_means: dict[str, float]
    posterior_std: dict[str, float]
    credible_intervals_95: dict[str, tuple[float, float]]
    acceptance_fraction: float


def log_prior(theta: np.ndarray) -> float:
    """Log-prior on kinetic parameters: log-uniform across physiological rates."""
    # Bounds: kon in [1e-4, 1.0], koff in [1e-3, 10.0], beta in [0.01, 10.0],
    # alpha in [0.001, 5.0], d in [0.001, 2.0], r in [1e-5, 0.5]
    bounds = [
        (1e-4, 1.0),    # kon
        (1e-3, 10.0),   # koff
        (0.01, 10.0),   # beta
        (0.001, 5.0),   # alpha
        (0.001, 2.0),   # d
        (1e-5, 0.5),    # r
    ]
    for val, (lo, hi) in zip(theta, bounds):
        if val <= lo or val >= hi:
            return -np.inf
    return 0.0


def log_posterior(theta: np.ndarray, param_keys: Sequence[str]) -> float:
    """Unnormalized log-posterior = log-prior + log-likelihood."""
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    # log_lik = -0.5 * chi2
    chi2 = compute_kinetic_objective(theta, param_keys)
    return float(lp - 0.5 * chi2)


def run_ensemble_mcmc(
    x0: np.ndarray,
    param_keys: Sequence[str] = ("kon", "koff", "beta", "alpha", "d", "r"),
    n_walkers: int = 16,
    n_steps: int = 150,
    burn_in: int = 50,
    a: float = 2.0,
    seed: int = 42,
) -> MCMCChainResult:
    """Run affine-invariant ensemble MCMC over kinetic parameters.

    Parameters
    ----------
    x0 : np.ndarray
        Central starting parameter vector (length D).
    param_keys : Sequence[str]
        Parameter names.
    n_walkers : int
        Number of ensemble walkers (must be >= 2 * D).
    n_steps : int
        Number of MCMC steps per walker.
    burn_in : int
        Number of initial steps to discard.
    a : float
        Scale parameter for the stretch move (typically 2.0).
    seed : int
        Random seed for reproducibility.
    """
    rng = np.random.default_rng(seed)
    dim = len(x0)
    assert n_walkers >= 2 * dim, f"n_walkers must be >= 2*dim ({2*dim})"

    # Initialize walker swarm in a small ball around x0
    walkers = np.zeros((n_walkers, dim), dtype=float)
    for i in range(n_walkers):
        perturb = rng.normal(0.0, 0.05, dim)
        walkers[i] = np.maximum(x0 * (1.0 + perturb), 1e-4)

    # Evaluate initial log-posteriors
    log_probs = np.array([log_posterior(w, param_keys) for w in walkers])
    # Fallback if any started in invalid region
    for i in range(n_walkers):
        if not np.isfinite(log_probs[i]):
            walkers[i] = x0.copy()
            log_probs[i] = log_posterior(walkers[i], param_keys)

    chain = np.zeros((n_steps, n_walkers, dim), dtype=float)
    n_accepted = 0
    total_proposals = 0

    for step in range(n_steps):
        # Stretch move: update each walker against a random complementary walker
        for k in range(n_walkers):
            # Pick a complementary walker j != k
            j = int(rng.choice([idx for idx in range(n_walkers) if idx != k]))
            # Sample z ~ g(z) = 1 / (2*(sqrt(a) - 1/sqrt(a)) * sqrt(z))
            u = rng.uniform(0.0, 1.0)
            z = ((u * (np.sqrt(a) - 1.0 / np.sqrt(a)) + 1.0 / np.sqrt(a)) ** 2)
            
            y = walkers[j] + z * (walkers[k] - walkers[j])
            new_lp = log_posterior(y, param_keys)
            
            # Acceptance prob: min(1, z^(D-1) * P(Y) / P(X))
            log_accept = (dim - 1) * np.log(z) + (new_lp - log_probs[k])
            if np.log(rng.uniform(0.0, 1.0)) < log_accept:
                walkers[k] = y
                log_probs[k] = new_lp
                n_accepted += 1
            total_proposals += 1

        chain[step] = walkers.copy()

    # Discard burn-in
    samples = chain[burn_in:].reshape(-1, dim)
    acc_rate = float(n_accepted / max(total_proposals, 1))

    # Compute posterior stats
    means = {k: float(np.mean(samples[:, i])) for i, k in enumerate(param_keys)}
    stds = {k: float(np.std(samples[:, i])) for i, k in enumerate(param_keys)}
    cis = {
        k: (float(np.percentile(samples[:, i], 2.5)), float(np.percentile(samples[:, i], 97.5)))
        for i, k in enumerate(param_keys)
    }

    return MCMCChainResult(
        param_names=tuple(param_keys),
        flat_samples=samples,
        posterior_means=means,
        posterior_std=stds,
        credible_intervals_95=cis,
        acceptance_fraction=acc_rate,
    )
