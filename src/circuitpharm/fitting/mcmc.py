"""Affine-invariant ensemble MCMC (Goodman & Weare 2010 stretch move), pure NumPy.

READ `fitting/identifiability.py`'s module docstring FIRST. This samples the same
objective, which is a function of (K_d, E, D) only, so the posterior is flat along three
of six directions and bounded there only by the prior box. Its credible intervals are
therefore PRIOR WIDTHS, not data-driven intervals, and they come back `Tier.VOID`
(roadmap P0-5). The chain itself is a real sample from a real distribution and stays plain.

Three things fixed here beyond the VOID marking:

  * THE PRIOR MATCHED ITS NAME. `log_prior` was documented "log-uniform" and returned 0.0
    inside the box, which is uniform ON THE RATES. For rate constants spanning decades
    that is a materially different prior -- it puts most of its mass at the top of each
    range. Sampling now happens in log10 space, where a flat prior IS log-uniform, which
    also removes the positivity constraint for free.
  * CONVERGENCE DIAGNOSTICS GATE THE RESULT. Integrated autocorrelation time, chain length
    in multiples of it, split-Rhat and acceptance fraction are computed and a run that
    fails any of them returns a result whose summaries name the failing diagnostic. The
    previous defaults -- 16 walkers x 150 steps for 6 dimensions -- are far short of what
    is needed and nothing said so.
  * THE STRETCH MOVE'S JACOBIAN IS APPLIED IN THE SAMPLING SPACE, which is now log10.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from ..provisional import void, warn_provisional
from ..results import Quantity
from .data import DOSE_RESPONSE_BENCHMARK, assert_real_data
from .identifiability import equilibrium_dose_response_chi2

_DEFECT = ("sampled from an objective that is a function of (K_d, E, D) only, so the "
           "posterior is flat along three of six directions and bounded there by the "
           "prior box -- the reported intervals are prior widths, not data")
_PROMOTE = ("P3/P4: sample the identifiable reparameterisation against a likelihood that "
            "includes a kinetic observable, with each dataset on its own observable")

#: log10 bounds, one per rate. A FLAT prior in this space is log-uniform in the rate,
#: which is what the docstring always claimed. Each is a stated range, not a measurement;
#: `provenance.PRIOR_PROV` carries the bases.
LOG10_BOUNDS = {
    "kon": (-4.0, 0.0),      # uM^-1 ms^-1
    "koff": (-3.0, 1.0),     # ms^-1
    "beta": (-2.0, 1.0),     # ms^-1
    "alpha": (-3.0, 0.7),    # ms^-1
    "d": (-3.0, 0.3),        # ms^-1
    "r": (-5.0, -0.3),       # ms^-1
}


@dataclass(frozen=True)
class MCMCChainResult:
    """Posterior sample, with its summaries marked VOID and its diagnostics plain.

    `flat_samples` (in the NATIVE rate space), `log10_samples`, `acceptance_fraction` and
    `diagnostics` are real computations. `posterior_means`, `posterior_std` and
    `credible_intervals_95` are inferences from a structurally flat posterior.
    """

    param_names: tuple[str, ...]
    flat_samples: np.ndarray
    log10_samples: np.ndarray
    acceptance_fraction: float
    diagnostics: dict
    converged: bool
    posterior_means: dict[str, Quantity]
    posterior_std: dict[str, Quantity]
    credible_intervals_95: dict[str, Quantity]
    defect: str = _DEFECT


def log_prior_log10(theta_log10: np.ndarray, param_keys: Sequence[str]) -> float:
    """Flat prior inside `LOG10_BOUNDS`, which is LOG-UNIFORM in the rate itself.

    The old `log_prior` returned 0.0 inside a box on the RATES and called itself
    log-uniform. Those are different priors: uniform on [1e-3, 10] puts ~99.9% of its mass
    above 0.01, while log-uniform spreads it evenly over the decades. For parameters whose
    plausible values span four orders of magnitude that changes the posterior, so the
    sampler now works in the space where the claim is true.
    """
    for val, key in zip(theta_log10, param_keys):
        lo, hi = LOG10_BOUNDS[key]
        if not (lo < float(val) < hi):
            return -np.inf
    return 0.0


def log_posterior_log10(theta_log10: np.ndarray, param_keys: Sequence[str]) -> float:
    lp = log_prior_log10(theta_log10, param_keys)
    if not np.isfinite(lp):
        return -np.inf
    chi2 = equilibrium_dose_response_chi2(10.0 ** np.asarray(theta_log10, dtype=float),
                                          param_keys)
    if not np.isfinite(chi2):
        return -np.inf
    return float(lp - 0.5 * chi2)


# -------------------------------------------------------------------- diagnostics
def integrated_autocorr_time(chain: np.ndarray, c: float = 5.0) -> float:
    """Integrated autocorrelation time of one parameter, with Sokal's adaptive window.

    `chain` is (n_steps, n_walkers). Returns NaN when the chain is too short for the
    window to close, which is itself the diagnostic -- it must not be reported as a small
    number.
    """
    x = np.asarray(chain, dtype=float)
    if x.ndim != 2 or x.shape[0] < 10:
        return float("nan")
    n = x.shape[0]
    y = x - x.mean(axis=0, keepdims=True)
    # mean autocovariance over walkers, via FFT
    f = np.fft.rfft(y, n=2 * n, axis=0)
    acf = np.fft.irfft(f * np.conjugate(f), n=2 * n, axis=0)[:n].real.mean(axis=1)
    if acf[0] <= 0:
        return float("nan")
    acf /= acf[0]
    taus = 2.0 * np.cumsum(acf) - 1.0
    window = np.flatnonzero(np.arange(len(taus)) < c * taus)
    if window.size == 0:
        return float("nan")
    return float(taus[window[-1]])


def split_rhat(chain: np.ndarray) -> float:
    """Split-Rhat for one parameter. `chain` is (n_steps, n_walkers)."""
    x = np.asarray(chain, dtype=float)
    n_steps, n_walkers = x.shape
    half = n_steps // 2
    if half < 4:
        return float("nan")
    segs = np.concatenate([x[:half], x[half:2 * half]], axis=1)   # 2 * n_walkers segments
    m = segs.shape[1]
    means = segs.mean(axis=0)
    variances = segs.var(axis=0, ddof=1)
    w = float(variances.mean())
    b = float(half * means.var(ddof=1))
    if w <= 0:
        return float("nan")
    var_hat = ((half - 1) / half) * w + b / half
    return float(np.sqrt(var_hat / w))


def run_ensemble_mcmc(
    x0: np.ndarray,
    param_keys: Sequence[str] = ("kon", "koff", "beta", "alpha", "d", "r"),
    n_walkers: int = 32,
    n_steps: int = 2000,
    burn_in: int = 500,
    a: float = 2.0,
    seed: int = 42,
    min_tau_multiples: float = 50.0,
    max_rhat: float = 1.01,
    accept_range: tuple[float, float] = (0.15, 0.60),
) -> MCMCChainResult:
    """Sample the posterior in log10 space. Summaries are VOID; see the module docstring.

    Defaults are 32 walkers x 2000 steps rather than the previous 16 x 150, which for six
    dimensions is nowhere near enough to estimate an autocorrelation time, let alone a
    credible interval. A run that fails a diagnostic still returns -- the chain is real --
    but `converged` is False and the summaries say which diagnostic failed.
    """
    param_keys = tuple(param_keys)
    assert_real_data([DOSE_RESPONSE_BENCHMARK], caller="run_ensemble_mcmc")
    warn_provisional(what="the MCMC posterior", defect=_DEFECT,
                     register_item="P0-5", promote_by=_PROMOTE)

    dim = len(param_keys)
    if n_walkers < 2 * dim:
        raise ValueError(f"n_walkers must be >= 2*dim ({2 * dim}), got {n_walkers}")
    if burn_in >= n_steps:
        raise ValueError(f"burn_in ({burn_in}) must be < n_steps ({n_steps})")

    rng = np.random.default_rng(seed)
    centre = np.log10(np.maximum(np.asarray(x0, dtype=float), 1e-12))
    lo = np.array([LOG10_BOUNDS[k][0] for k in param_keys])
    hi = np.array([LOG10_BOUNDS[k][1] for k in param_keys])
    if np.any(centre <= lo) or np.any(centre >= hi):
        raise ValueError(
            f"x0 maps to log10 {centre} which is outside LOG10_BOUNDS "
            f"[{lo}, {hi}] for {param_keys}. The walkers would all start at -inf "
            f"log-posterior and the old code silently reset them to x0, which collapses "
            f"the ensemble to a point.")

    walkers = np.clip(centre + rng.normal(0.0, 0.05, (n_walkers, dim)),
                      lo + 1e-9, hi - 1e-9)
    log_probs = np.array([log_posterior_log10(w, param_keys) for w in walkers])
    if not np.all(np.isfinite(log_probs)):
        raise RuntimeError(
            "some walkers started with a non-finite log-posterior. Resetting them to the "
            "centre (the previous behaviour) collapses the ensemble and makes the stretch "
            "move degenerate; widen LOG10_BOUNDS or move x0 instead.")

    chain = np.zeros((n_steps, n_walkers, dim), dtype=float)
    n_accepted = 0
    total = 0
    sqrt_a = np.sqrt(a)

    for step in range(n_steps):
        for k in range(n_walkers):
            j = int(rng.integers(0, n_walkers - 1))
            if j >= k:
                j += 1                       # uniform over the complementary walkers
            z = (rng.uniform(0.0, 1.0) * (sqrt_a - 1.0 / sqrt_a) + 1.0 / sqrt_a) ** 2
            y = walkers[j] + z * (walkers[k] - walkers[j])
            new_lp = log_posterior_log10(y, param_keys)
            log_accept = (dim - 1) * np.log(z) + (new_lp - log_probs[k])
            if np.log(rng.uniform(0.0, 1.0)) < log_accept:
                walkers[k] = y
                log_probs[k] = new_lp
                n_accepted += 1
            total += 1
        chain[step] = walkers

    post = chain[burn_in:]
    log10_samples = post.reshape(-1, dim)
    samples = 10.0 ** log10_samples
    acc = float(n_accepted / max(total, 1))

    taus = {k: integrated_autocorr_time(post[:, :, i]) for i, k in enumerate(param_keys)}
    rhats = {k: split_rhat(post[:, :, i]) for i, k in enumerate(param_keys)}
    n_kept = post.shape[0]
    worst_tau = max((t for t in taus.values() if np.isfinite(t)), default=float("nan"))
    worst_rhat = max((r for r in rhats.values() if np.isfinite(r)), default=float("nan"))

    failures = []
    if not np.isfinite(worst_tau):
        failures.append("integrated autocorrelation time could not be estimated "
                        "(chain too short for the window to close)")
    elif n_kept < min_tau_multiples * worst_tau:
        failures.append(f"chain is {n_kept / worst_tau:.1f} autocorrelation times long, "
                        f"need >= {min_tau_multiples:g}")
    if not np.isfinite(worst_rhat):
        failures.append("split-Rhat could not be estimated")
    elif worst_rhat > max_rhat:
        failures.append(f"worst split-Rhat {worst_rhat:.4f} > {max_rhat}")
    if not (accept_range[0] <= acc <= accept_range[1]):
        failures.append(f"acceptance fraction {acc:.3f} outside {accept_range}")

    converged = not failures
    diagnostics = dict(tau=taus, split_rhat=rhats, n_kept=n_kept,
                       acceptance_fraction=acc, failures=tuple(failures),
                       converged=converged)

    # A failed diagnostic is an ADDITIONAL reason the summaries are unusable; they are VOID
    # either way because of the flat objective.
    extra = (("diagnostics failed: " + "; ".join(failures)),) if failures else ()
    def _v(name, value, note=()):
        return void(name, value, defect=_DEFECT, register_item="P0-5",
                    promote_by=_PROMOTE, caveats=extra + tuple(note))

    means = {k: _v(f"posterior_mean_{k}", float(np.mean(samples[:, i])))
             for i, k in enumerate(param_keys)}
    stds = {k: _v(f"posterior_std_{k}", float(np.std(samples[:, i])))
            for i, k in enumerate(param_keys)}
    cis = {k: _v(f"posterior_ci95_{k}",
                 (float(np.percentile(samples[:, i], 2.5)),
                  float(np.percentile(samples[:, i], 97.5))),
                 note=("for a flat direction this interval is the PRIOR range, not an "
                       "inference",))
           for i, k in enumerate(param_keys)}

    return MCMCChainResult(
        param_names=param_keys,
        flat_samples=samples,
        log10_samples=log10_samples,
        acceptance_fraction=acc,
        diagnostics=diagnostics,
        converged=converged,
        posterior_means=means,
        posterior_std=stds,
        credible_intervals_95=cis,
    )


def log_prior(theta: np.ndarray) -> float:
    """Deprecated. The old rate-space box prior, kept so nothing breaks.

    It is UNIFORM ON THE RATES despite its former docstring saying log-uniform. Use
    `log_prior_log10`, where a flat prior genuinely is log-uniform.
    """
    bounds = [LOG10_BOUNDS[k] for k in ("kon", "koff", "beta", "alpha", "d", "r")]
    for val, (lo, hi) in zip(theta, bounds):
        if not (10.0 ** lo < float(val) < 10.0 ** hi):
            return -np.inf
    return 0.0
