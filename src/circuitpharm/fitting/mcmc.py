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
    means = segs.mean(axis=0)
    variances = segs.var(axis=0, ddof=1)
    w = float(variances.mean())
    b = float(half * means.var(ddof=1))
    if w <= 0:
        return float("nan")
    var_hat = ((half - 1) / half) * w + b / half
    return float(np.sqrt(var_hat / w))


@dataclass(frozen=True)
class EnsembleRun:
    """A finished ensemble run in whatever space it was asked to sample.

    Nothing here is VOID: a chain is a chain, and its diagnostics are measurements of that
    chain. Whether the SUMMARIES mean anything depends on the objective, which is the
    caller's business -- `run_ensemble_mcmc` below wraps this one in VOID quantities
    because its objective is flat in three directions, while `posterior.py` wraps the same
    machinery in real intervals because its objective is not.
    """

    param_names: tuple[str, ...]
    #: (n_kept, n_walkers, dim), post burn-in, in the sampling space.
    chain: np.ndarray
    #: (n_kept * n_walkers, dim), the same thing flattened.
    samples: np.ndarray
    acceptance_fraction: float
    diagnostics: dict
    converged: bool


def sample_ensemble(
    log_prob,
    *,
    lo: np.ndarray,
    hi: np.ndarray,
    names: Sequence[str],
    x0: np.ndarray | None = None,
    init: str = "prior",
    n_walkers: int = 32,
    n_steps: int = 2000,
    burn_in: int = 500,
    a: float | Sequence[float] = 2.0,
    seed: int = 42,
    jitter: float = 0.05,
    min_tau_multiples: float = 50.0,
    max_rhat: float = 1.01,
    accept_range: tuple[float, float] = (0.15, 0.60),
    min_walker_share: float = 0.1,
) -> EnsembleRun:
    """Goodman & Weare stretch-move ensemble sampler over an arbitrary log-density.

    `x0`, `lo` and `hi` are in the SAMPLING space already -- this function does no
    transformation, so the caller decides whether it is sampling rates or their logarithms
    and the Jacobian is the caller's to get right. The stretch move's own `(dim-1)log z`
    Jacobian is applied here, in the sampling space, which is the only place it is correct.

    Extracted from `run_ensemble_mcmc` so P4's posterior over the identifiable combinations
    runs the same sampler and the same diagnostics rather than a second copy of them: "the
    fix was applied in one of several places" is this project's recurring error E12.

    `init="prior"` spreads the walkers uniformly over the box, which is what makes
    split-Rhat mean anything: the statistic compares between-chain to within-chain variance
    and is near 1 by construction if every walker starts in the same small ball.
    `init="ball"` is the alternative (Gaussian of width `jitter` about `x0`), which is
    standard practice around an MLE but only diagnoses mixing within the mode it starts in.

    ON `a`, AND A FROZEN-WALKER TRAP FOUND HERE. The stretch move proposes
    `y = x_j + z(x_k - x_j)` with `z` in `[1/a, a]`, so with the canonical `a = 2` it CANNOT
    propose a contraction tighter than one half. A walker sitting in a low-probability
    region separated from the ensemble by a valley is then stuck: the nearest point it can
    reach is the midpoint, and if that midpoint is worse than where it is, every proposal is
    rejected for the whole run. Measured on P4's posterior (`fitting/posterior.py`), from a
    prior-dispersed start: 2 of 24 walkers accepted 94 and 405 moves against a median of
    3350, froze 2.5 decades from the mode, and inflated split-Rhat to 62 while the other 22
    walkers sat on the right answer. Along the line from one of them to the mode the density
    fell by 16 log units before rising, so the move it needed was `z = 0.14` and the move it
    could make was `z = 0.5`.

    `a` MAY THEREFORE BE SEVERAL SCALES, sampled uniformly per proposal. A single large `a`
    crosses the valley but wastes most proposals: on that posterior `a = 200` gave
    acceptance 0.045 and tau 175, against 0.49 and tau 57 for `a = 2`, which cannot cross
    it at all. Mixing `(2, 20, 200)` gives acceptance 0.28, tau 59, split-Rhat 1.006 and no
    frozen walker -- better than either alone on every diagnostic. Each fixed-`a` stretch
    move satisfies detailed balance with respect to the target, and choosing between them
    independently of the current state leaves that intact, so the mixture is a valid kernel
    rather than a heuristic.

    The freezing is now DETECTED rather than relied on being noticed: per-walker acceptance
    is reported, and a walker that moved less than `min_walker_share` of the median is a
    diagnostic failure, because an ensemble with a frozen member has not sampled anything.
    """
    names = tuple(names)
    dim = len(names)
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    if lo.shape != (dim,) or hi.shape != (dim,):
        raise ValueError(
            f"lo and hi must each have one entry per name ({dim}); got "
            f"{lo.shape}, {hi.shape}")
    if np.any(hi <= lo):
        raise ValueError(f"every bound must have hi > lo; got lo={lo}, hi={hi}")
    if n_walkers < 2 * dim:
        raise ValueError(f"n_walkers must be >= 2*dim ({2 * dim}), got {n_walkers}")
    if burn_in >= n_steps:
        raise ValueError(f"burn_in ({burn_in}) must be < n_steps ({n_steps})")
    if init not in ("prior", "ball"):
        raise ValueError(f"init must be 'prior' or 'ball', got {init!r}")
    scales = (float(a),) if np.isscalar(a) else tuple(float(x) for x in a)
    if not scales or any(x <= 1.0 for x in scales):
        raise ValueError(
            f"every stretch scale must be > 1 (a scale bounds the stretch factor to "
            f"z in [1/a, a]); got {a!r}. a = 1 proposes y = x_k and never moves.")

    rng = np.random.default_rng(seed)
    if init == "ball":
        if x0 is None:
            raise ValueError("init='ball' needs x0")
        centre = np.asarray(x0, dtype=float)
        if centre.shape != (dim,):
            raise ValueError(f"x0 must have one entry per name ({dim}); got {centre.shape}")
        if np.any(centre <= lo) or np.any(centre >= hi):
            raise ValueError(
                f"x0 {centre} is outside the box [{lo}, {hi}] for {names}. Every walker "
                f"would start at -inf log-density, and resetting them to x0 -- the previous "
                f"behaviour -- collapses the ensemble to a point, which makes the stretch "
                f"move degenerate rather than merely inefficient.")
        walkers = np.clip(centre + rng.normal(0.0, jitter, (n_walkers, dim)),
                          lo + 1e-9, hi - 1e-9)
    else:
        walkers = lo + rng.uniform(0.0, 1.0, (n_walkers, dim)) * (hi - lo)

    log_probs = np.array([log_prob(w) for w in walkers], dtype=float)
    if not np.all(np.isfinite(log_probs)):
        bad = int(np.sum(~np.isfinite(log_probs)))
        raise RuntimeError(
            f"{bad} of {n_walkers} walkers started with a non-finite log-density. Widen "
            f"the box, move x0, or -- for init='prior' -- give the density a finite value "
            f"everywhere inside the box; resetting them collapses the ensemble.")

    chain = np.zeros((n_steps, n_walkers, dim), dtype=float)
    accepts = np.zeros(n_walkers, dtype=np.int64)
    sqrt_scales = np.sqrt(np.array(scales, dtype=float))
    n_scales = len(scales)

    for step in range(n_steps):
        for k in range(n_walkers):
            j = int(rng.integers(0, n_walkers - 1))
            if j >= k:
                j += 1                       # uniform over the complementary walkers
            # one scale per proposal, chosen independently of the state: a mixture of
            # valid kernels is a valid kernel, and no single scale does both jobs.
            sqrt_a = (sqrt_scales[0] if n_scales == 1
                      else sqrt_scales[int(rng.integers(0, n_scales))])
            z = (rng.uniform(0.0, 1.0) * (sqrt_a - 1.0 / sqrt_a) + 1.0 / sqrt_a) ** 2
            y = walkers[j] + z * (walkers[k] - walkers[j])
            new_lp = log_prob(y)
            log_accept = (dim - 1) * np.log(z) + (new_lp - log_probs[k])
            if np.log(rng.uniform(0.0, 1.0)) < log_accept:
                walkers[k] = y
                log_probs[k] = new_lp
                accepts[k] += 1
        chain[step] = walkers

    post = chain[burn_in:]
    total = n_steps * n_walkers
    acc = float(accepts.sum() / max(total, 1))
    taus = {k: integrated_autocorr_time(post[:, :, i]) for i, k in enumerate(names)}
    rhats = {k: split_rhat(post[:, :, i]) for i, k in enumerate(names)}
    n_kept = post.shape[0]
    worst_tau = max((t for t in taus.values() if np.isfinite(t)), default=float("nan"))
    worst_rhat = max((r for r in rhats.values() if np.isfinite(r)), default=float("nan"))
    median_accepts = float(np.median(accepts))
    frozen = [int(i) for i in np.flatnonzero(accepts < min_walker_share * median_accepts)]

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
    if frozen:
        failures.append(
            f"walker(s) {frozen} accepted {[int(accepts[i]) for i in frozen]} moves "
            f"against a median of {median_accepts:.0f} -- they are frozen, not sampling. "
            f"The stretch move cannot contract below 1/max(a) = "
            f"{1 / max(scales):.3g}, so a walker across a valley from the ensemble can be "
            f"unable to reach it; add a larger scale to `a`")

    diagnostics = dict(tau=taus, split_rhat=rhats, n_kept=n_kept,
                       acceptance_fraction=acc, accepts_per_walker=accepts.copy(),
                       frozen_walkers=tuple(frozen), stretch_a=scales, init=init,
                       failures=tuple(failures), converged=not failures)
    return EnsembleRun(param_names=names, chain=post,
                       samples=post.reshape(-1, dim), acceptance_fraction=acc,
                       diagnostics=diagnostics, converged=not failures)


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

    centre = np.log10(np.maximum(np.asarray(x0, dtype=float), 1e-12))
    try:
        run = sample_ensemble(
            lambda v: log_posterior_log10(v, param_keys),
            x0=centre, init="ball",
            lo=np.array([LOG10_BOUNDS[k][0] for k in param_keys]),
            hi=np.array([LOG10_BOUNDS[k][1] for k in param_keys]),
            names=param_keys, n_walkers=n_walkers, n_steps=n_steps, burn_in=burn_in,
            a=a, seed=seed, min_tau_multiples=min_tau_multiples, max_rhat=max_rhat,
            accept_range=accept_range)
    except ValueError as exc:
        # The box check speaks in the sampling space; say which space that is, since the
        # caller passed rates.
        if "outside the box" in str(exc):
            raise ValueError(
                f"x0 maps to log10 {centre}, which is outside LOG10_BOUNDS for "
                f"{param_keys}. {exc}") from None
        raise

    log10_samples = run.samples
    samples = 10.0 ** log10_samples
    acc = run.acceptance_fraction
    diagnostics = run.diagnostics
    converged = run.converged
    failures = list(diagnostics["failures"])

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
