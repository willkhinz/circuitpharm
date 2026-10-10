"""A posterior over the combinations the data can determine (roadmap P4-4).

WHY THIS EXISTS ALONGSIDE `mcmc.py`. That module samples six microscopic rates against an
equilibrium objective that is a function of three ratios, so three of its six directions
are flat and bounded only by the prior box. Its credible intervals are prior widths, and
they come back `Tier.VOID` saying so. The machinery is not wrong -- the *question* was.
This module asks the answerable version: sample `(log10 K_d, log10 E, log10 D)`, which are
exactly the combinations `reparam.EQUILIBRIUM_IDENTIFIABLE` proves an equilibrium
observable determines, against the per-observable likelihood in `likelihood.py`.

FOUR THINGS THAT MAKE THE INTERVALS MEAN SOMETHING:

  * THE SAMPLED SPACE IS THE IDENTIFIABLE ONE. No direction of this posterior is
    structurally flat, so no interval here is a restatement of its bounds. The bounds are
    still reported, and `sampler_pinned_at_bound` says if the chain spent its time against
    one -- which is PRACTICAL non-identifiability and a result, not a failure.
  * A FLAT PRIOR IN log10 IS LOG-UNIFORM IN THE RATIO, which is the right uninformative
    choice for a quantity spanning decades and is what sampling in this space buys.
  * THE DIAGNOSTICS GATE THE RESULT. A chain that fails tau, split-Rhat or acceptance
    returns intervals that are VOID -- reading one raises -- rather than a number with a
    footnote. An unconverged chain's percentiles are not an interval; they are where the
    walkers happened to be.
  * TWO INDEPENDENT MACHINERIES MUST AGREE. `agreement_with_profiles` checks each
    posterior median against the P3 profile-likelihood 95% CI for the same parameter.
    Profile likelihood and MCMC share no code and almost no assumptions, so agreement is
    the strongest check available here and disagreement means one of them is wrong. The
    check is a test, not decoration.

ON PROFILING SIGMA OUT AND THEN SAMPLING. `likelihood.concentrated_log_likelihood`
substitutes the MLE `sigma^2 = RSS/n`, giving `const - (n/2) log RSS`. Profiling a nuisance
parameter is not in general the same as marginalising it -- but here it happens to be, up
to a constant: marginalising `sigma` under the scale-invariant (Jeffreys) prior `1/sigma`
gives `Gamma(n/2) 2^(n/2-1) RSS^(-n/2)`, whose logarithm is `const - (n/2) log RSS` with
the same parameter dependence. So this posterior IS the marginal posterior over the three
ratios with `sigma` integrated out under a Jeffreys prior, and the constant it drops
depends only on `n`, which is fixed. Said explicitly because "profile == marginal" is
false in general and would be a real defect if assumed rather than checked.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..provisional import void
from ..results import Quantity, ResultSet, Tier
from .identifiability import IDENTIFIABLE_BOUNDS
from .likelihood import make_log_likelihood
from .mcmc import sample_ensemble
from .reparam import EQUILIBRIUM_IDENTIFIABLE, IdentifiableParams

_UNCONVERGED_PROMOTE = (
    "run the chain longer (n_steps >= 50 integrated autocorrelation times after burn-in) "
    "and re-check split-Rhat; the sampler reports both, so the target is measurable rather "
    "than guessed")


@dataclass(frozen=True)
class PosteriorResult:
    """A posterior over the identifiable combinations, with its summaries gated.

    `log10_samples`, `acceptance_fraction` and `diagnostics` are measurements of the chain
    and are always plain. `medians` and `intervals_95` are inferences and are VOID when a
    diagnostic failed, UNCALIBRATED when the data behind them is parametric or synthetic,
    and VALIDATED only when the chain converged against measured data.
    """

    param_names: tuple[str, ...]
    log10_samples: np.ndarray
    acceptance_fraction: float
    diagnostics: dict
    converged: bool
    medians: dict[str, Quantity]
    intervals_95: dict[str, Quantity]
    #: Parameters at the posterior median, as an `IdentifiableParams`, or `None` for a
    #: parameterisation that does not fit in that dataclass (Model C has four). Plain: it
    #: is a coordinate of the sample, and the tier lives on `medians`.
    median_params: IdentifiableParams | None
    #: Per parameter: did the 95% interval run into its prior bound? True means practical
    #: non-identifiability on that side -- the data did not close the interval, the box did.
    pinned_at_bound: dict[str, tuple[bool, bool]]
    bounds: dict[str, tuple[float, float]]
    tier: Tier
    note: str


def sample_identifiable_posterior(model_factory, datasets, **kw) -> PosteriorResult:
    """Sample `(log10 K_d, log10 E, log10 D)`; see `sample_posterior_vector`.

    `model_factory` takes an `IdentifiableParams`. A thin wrapper, because P5's Model C has
    four identifiable combinations and does not fit in that dataclass -- the general form
    takes a plain vector.
    """
    names = tuple(EQUILIBRIUM_IDENTIFIABLE)
    return sample_posterior_vector(
        lambda v: model_factory(IdentifiableParams.from_vector(v)), datasets,
        param_names=names,
        bounds={k: IDENTIFIABLE_BOUNDS[k] for k in names} | dict(kw.pop("bounds", {})
                                                                 or {}),
        caller="sample_identifiable_posterior", **kw)


def sample_posterior_vector(
    model_factory,
    datasets,
    *,
    param_names,
    bounds,
    caller: str = "sample_posterior_vector",
    x0=None,
    n_walkers: int = 24,
    # 24000/8000, RAISED FROM 12000/4000 BECAUSE THE OLD DEFAULT SAT ON THE GATE.
    # At 12000/4000 the P4 anchor posterior reaches split-Rhat = 1.01089 against
    # max_rhat = 1.01 -- it FAILS by 0.0009. And it fails only on some machines: the chain is
    # seeded, so it is deterministic per platform, but arm64/Accelerate and x86/OpenBLAS sum
    # in different orders and that propagates through the matrix exponentials and 24 walkers
    # into the third decimal of Rhat. CI (x86) was permanently green and a macOS checkout
    # permanently red, on identical code.
    #
    # A gate that discriminates finer than the diagnostic's own platform reproducibility is
    # not measuring convergence, it is measuring the BLAS. The gate is right -- 1.01 is the
    # Vehtari et al. recommendation, and `converged` decides whether percentiles are reported
    # as credible intervals or marked VOID -- so the chain was lengthened to meet it with
    # margin rather than the threshold loosened to admit it.
    #
    # Measured on this scheme (worst split-Rhat, 24 walkers, seed 11):
    #     12000/4000  1.01089  FAIL      3.4 s
    #     18000/6000  1.00521            5.2 s
    #     24000/8000  1.00377            7.1 s   <- chosen: 0.0062 of margin, ~3x the
    #     36000/12000 1.00240           10.7 s      shortfall it had to cover
    #     48000/16000 1.00257           14.2 s   <- plateau; more steps buy nothing
    n_steps: int = 24000,
    burn_in: int = 8000,
    seed: int = 11,
    stretch_a=(2.0, 20.0, 200.0),
    min_tau_multiples: float = 50.0,
    max_rhat: float = 1.01,
    accept_range: tuple[float, float] = (0.15, 0.60),
) -> PosteriorResult:
    """Sample a log10 parameter vector against the per-observable likelihood.

    The caller supplies `model_factory`, for the same reason `fit_mle` requires one:
    turning three ratios into six rates needs conventions for `alpha` and `r`, and this
    module will not invent them (see `reparam.IdentifiableParams.to_microscopic`).

    THE WALKERS START FROM THE PRIOR, not from an MLE. A chain started in a small ball at
    the optimum reports a split-Rhat near 1 whether or not it mixed, because the statistic
    compares between-walker to within-walker variance and the walkers began identical. Over-
    dispersed starts are what make the diagnostic a test. Pass `x0` only to reproduce a
    particular run.

    `stretch_a` IS A MIXTURE, not the canonical 2. With `a = 2` alone the stretch move
    cannot propose a contraction below one half, and on this posterior that froze 2 of 24
    walkers 2.5 decades from the mode for an entire 6000-step run; with `a = 200` alone
    acceptance falls to 0.045 and tau rises to 175. The mixture gets acceptance 0.28, tau 59
    and split-Rhat 1.006 -- see `mcmc.sample_ensemble`'s docstring for the measurements and
    for why mixing scales is exact rather than a tuning trick.

    `model_factory` takes a plain `np.ndarray` in `param_names` order, every entry in
    log10. `bounds` maps each name to `(lo, hi)` and is required: a flat prior needs a box
    to be a prior at all, and leaving it implicit is how an interval becomes a restatement
    of a default nobody chose.
    """
    names = tuple(param_names)
    box = dict(bounds)
    missing = [n for n in names if n not in box]
    if missing:
        raise ValueError(f"no bounds given for {missing}; a flat prior needs a box")
    lo = np.array([box[k][0] for k in names], dtype=float)
    hi = np.array([box[k][1] for k in names], dtype=float)

    datasets = list(datasets)
    # ONCE, before the chain: the taint and holdout checks, and the closure the chain uses.
    score = make_log_likelihood(datasets, caller=caller, n_params=len(names))
    synthetic = [getattr(d, "synthetic", False) for d in datasets]

    def log_prob(v: np.ndarray) -> float:
        if np.any(v <= lo) or np.any(v >= hi):
            return -np.inf                      # flat prior == log-uniform in the ratio
        try:
            ll = score(model_factory(np.asarray(v, dtype=float))).total_log_likelihood
        except Exception:
            return -np.inf
        return ll if np.isfinite(ll) else -np.inf

    run = sample_ensemble(log_prob, lo=lo, hi=hi, names=names,
                          x0=None if x0 is None else np.asarray(x0, dtype=float),
                          init="prior" if x0 is None else "ball",
                          n_walkers=n_walkers, n_steps=n_steps, burn_in=burn_in,
                          seed=seed, a=stretch_a,
                          min_tau_multiples=min_tau_multiples, max_rhat=max_rhat,
                          accept_range=accept_range)

    s = run.samples
    med = {k: float(np.median(s[:, i])) for i, k in enumerate(names)}
    ci = {k: (float(np.percentile(s[:, i], 2.5)), float(np.percentile(s[:, i], 97.5)))
          for i, k in enumerate(names)}
    # "within 2% of the box width of a bound" -- a sampled percentile never lands exactly
    # on a bound, so an exact comparison would never fire and the check would be decorative.
    width = hi - lo
    pinned = {k: (bool(ci[k][0] - lo[i] < 0.02 * width[i]),
                  bool(hi[i] - ci[k][1] < 0.02 * width[i]))
              for i, k in enumerate(names)}

    failures = tuple(run.diagnostics["failures"])
    if failures:
        tier = Tier.VOID
        note = ("the chain failed " + "; ".join(failures)
                + ". Its percentiles are not credible intervals -- they are where the "
                  "walkers happened to be -- so they are VOID.")
    elif any(synthetic):
        tier = Tier.UNCALIBRATED
        note = ("the chain converged on every diagnostic, so the intervals are real "
                "statements about the likelihood; the DATA behind them is parametric or "
                "synthetic, so they are not statements about a receptor.")
    else:
        tier = Tier.VALIDATED
        note = ("the chain converged on every diagnostic and every dataset is measured.")

    pin_notes = tuple(
        f"{k}: the 95% interval reaches its prior bound "
        f"({'lower' if p[0] else ''}{' and ' if all(p) else ''}{'upper' if p[1] else ''}"
        f"), so the DATA did not close it on that side -- practically non-identifiable"
        for k, p in pinned.items() if any(p))

    def _wrap(name: str, value, extra: tuple[str, ...] = ()) -> Quantity:
        if tier is Tier.VOID:
            return void(name, value,
                        defect=("sampled from a chain that failed its convergence "
                                "diagnostics: " + "; ".join(failures)),
                        register_item="P4-4", promote_by=_UNCONVERGED_PROMOTE,
                        caveats=pin_notes + extra)
        return Quantity(
            name=name, _value=value, tier=tier,
            provenance=(f"posterior over {names} from the per-observable concentrated "
                        f"likelihood (sigma marginalised under a Jeffreys prior, see the "
                        f"module docstring), {run.diagnostics['n_kept']} kept steps x "
                        f"{n_walkers} walkers, seed {seed}; worst split-Rhat "
                        f"{max(v for v in run.diagnostics['split_rhat'].values()):.4f}, "
                        f"acceptance {run.acceptance_fraction:.3f}"),
            promote_by=("digitise a real EQUILIBRIUM concentration-response -- see "
                        "fitting.data.MISSING_DATASETS -- so the interval is about "
                        "receptors and not about a generated curve"
                        if tier is Tier.UNCALIBRATED else ""),
            caveats=pin_notes + extra
            + (("the prior is flat over " + repr(box) + " in log10, i.e. log-uniform in "
                "the ratio; an interval close to a bound is reporting the bound",)
               if any(any(p) for p in pinned.values()) else ()))

    medians = {k: _wrap(f"posterior_median_{k}", v) for k, v in med.items()}
    intervals = {
        k: _wrap(f"posterior_ci95_{k}", v,
                 extra=(f"width {v[1] - v[0]:.3f} decades, a factor of "
                        f"{10.0 ** (v[1] - v[0]):.3g}",))
        for k, v in ci.items()}

    median_params = (
        IdentifiableParams(log10_kd=med["log10_kd"], log10_E=med["log10_E"],
                           log10_D=med["log10_D"])
        if names == tuple(EQUILIBRIUM_IDENTIFIABLE) else None)
    return PosteriorResult(
        param_names=names, log10_samples=s,
        acceptance_fraction=run.acceptance_fraction, diagnostics=run.diagnostics,
        converged=run.converged, medians=medians, intervals_95=intervals,
        median_params=median_params,
        pinned_at_bound=pinned, bounds={k: (float(lo[i]), float(hi[i]))
                                        for i, k in enumerate(names)},
        tier=tier, note=note)


@dataclass(frozen=True)
class AgreementResult:
    """Posterior median against profile-likelihood CI, per parameter."""

    param_name: str
    posterior_median: float
    profile_ci_95: tuple[float, float]
    posterior_ci_95: tuple[float, float]
    median_inside_profile_ci: bool
    #: Overlap of the two 95% intervals as a fraction of the narrower one. 1.0 means one
    #: contains the other; 0.0 means they are disjoint, which means one machinery is wrong.
    interval_overlap: float


def agreement_with_profiles(post: PosteriorResult, dataset) -> dict[str, AgreementResult]:
    """Check the posterior against P3's profile likelihoods on the same dataset.

    THE POINT OF A SECOND MACHINERY. Profile likelihood re-optimises the other parameters
    at each node of a grid; MCMC integrates over them. They share no optimiser, no grid and
    no convergence criterion, so if they disagree about where the parameter is, at least one
    is broken -- and that is a far stronger statement than either passing its own internal
    checks. `tests/test_p4_likelihood.py` asserts on this.
    """
    from .identifiability import profile_likelihood

    out = {}
    for name in post.param_names:
        prof = profile_likelihood(name, dataset)
        pci = tuple(float(x) for x in prof.ci_95)
        median = post.medians[name].get(acknowledge_void=True)
        qci = post.intervals_95[name].get(acknowledge_void=True)
        inner_lo, inner_hi = max(pci[0], qci[0]), min(pci[1], qci[1])
        overlap = max(0.0, inner_hi - inner_lo)
        narrower = min(pci[1] - pci[0], qci[1] - qci[0])
        out[name] = AgreementResult(
            param_name=name, posterior_median=median, profile_ci_95=pci,
            posterior_ci_95=qci,
            median_inside_profile_ci=bool(pci[0] <= median <= pci[1]),
            interval_overlap=float(overlap / narrower) if narrower > 0 else 0.0)
    return out


def save_chain(post: PosteriorResult, path, datasets, *, label: str = "") -> str:
    """Archive a chain with everything needed to say what produced it (roadmap P4-2).

    A `.npz` holding the samples plus the git sha, the dataset keys, the dependency
    versions and the diagnostics. A chain without those is not reproducible and not worth
    keeping: six months on, "the posterior" in a figure caption has to resolve to a file
    that names its inputs, and a bare array of numbers cannot.

    The git sha is recorded as "unknown" rather than omitted if `git` is unavailable, so
    the field is always present and its absence is never mistaken for a clean tree.
    """
    import json
    import pathlib
    import platform
    import subprocess

    import scipy

    out = pathlib.Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)

    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                             timeout=10, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], capture_output=True,
                               text=True, timeout=10, check=True).stdout.strip()
        if dirty:
            sha += "-dirty"
    except Exception:
        sha = "unknown"

    meta = dict(
        git_sha=sha,
        label=label,
        param_names=list(post.param_names),
        tier=post.tier.name,
        note=post.note,
        converged=bool(post.converged),
        acceptance_fraction=float(post.acceptance_fraction),
        bounds={k: list(v) for k, v in post.bounds.items()},
        pinned_at_bound={k: list(v) for k, v in post.pinned_at_bound.items()},
        datasets=[dict(label=getattr(d, "citation_label", repr(d)),
                       observable=getattr(getattr(d, "observable", None), "value", None),
                       normalisation=getattr(d, "normalisation", None),
                       kind=getattr(d, "kind", None),
                       synthetic=bool(getattr(d, "synthetic", False)),
                       source_key=getattr(d, "source_key", "")) for d in datasets],
        diagnostics=dict(
            tau={k: float(v) for k, v in post.diagnostics["tau"].items()},
            split_rhat={k: float(v) for k, v in post.diagnostics["split_rhat"].items()},
            n_kept=int(post.diagnostics["n_kept"]),
            frozen_walkers=list(post.diagnostics["frozen_walkers"]),
            stretch_a=list(post.diagnostics["stretch_a"]),
            init=post.diagnostics["init"],
            failures=list(post.diagnostics["failures"]),
            accepts_per_walker=[int(x)
                                for x in post.diagnostics["accepts_per_walker"]],
        ),
        versions=dict(python=platform.python_version(), numpy=np.__version__,
                      scipy=scipy.__version__),
    )
    np.savez_compressed(out, log10_samples=post.log10_samples,
                        metadata=np.array(json.dumps(meta, indent=2)))
    return str(out)


def posterior_report(post: PosteriorResult) -> ResultSet:
    """The posterior as a tiered result set, for `audit` and the manuscript tables."""
    rs = ResultSet(f"posterior over {', '.join(post.param_names)} "
                   f"({post.tier.label}: {post.note})")
    for name in post.param_names:
        rs.add(post.medians[name])
        rs.add(post.intervals_95[name])
    return rs
