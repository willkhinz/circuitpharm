"""Dataset consistency diagnostic: asks whether datasets are mutually reconcilable.

THIS IS A DIAGNOSTIC, NOT A FIT. No parameter set found here is adopted as a fitted model,
and TRAIN / HOLDOUT roles are unchanged.

The question answered here is:
    "Does any parameter set in the model's feasible space satisfy all three real
     datasets simultaneously (Jahn CRC, Barberis deactivation, Barberis paired-pulse)?"

This is a property of the data and model structure, not a fitted receptor model.
Standard fitting paths (e.g. `make_log_likelihood`, `fit_one`, `compare_models`) enforce
`holdout_guard` to prevent data leakage into estimation. This diagnostic entry point is
explicitly separate and unreachable from TRAIN-path estimators.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence
import numpy as np
from scipy.optimize import minimize

from .data import (
    NoiseModel,
    assert_real_data,
)
from .likelihood import (
    DatasetFit,
    LikelihoodResult,
    _observed,
    _predict,
)
from ..models.base import Observable
from ..gabaa_kinetics import (
    FIT_RANGES,
    FIT_TARGETS,
    po_max_anchor_log_likelihood,
)
from scipy.stats import qmc


@dataclass(frozen=True)
class DiagnosticFitResult:
    """Result of a consistency diagnostic optimization.

    Explicitly labeled as diagnostic to prevent confusion with `ModelFit` or `MLEResult`.
    """
    values: dict[str, float]
    vector: np.ndarray
    total_log_likelihood: float
    per_dataset: tuple[DatasetFit, ...]
    po_max: float | None
    po_max_in_range: bool
    at_bound: tuple[str, ...]
    converged: bool
    note: str = ""


def diagnostic_log_likelihood(
    model,
    datasets: Sequence[object],
    *,
    n_params: int = 0,
    anchor_po_max: bool = True,
    weighting: str = "as_is",
    fixed_deact_sigma: float = 0.095091,
) -> LikelihoodResult:
    """Evaluate joint log-likelihood over arbitrary datasets for consistency diagnostics.

    Declared alternative weightings (Task C):
      - "as_is": Concentrated Gaussian for estimated noise (sigma profiled out), stated SEM
        for stated_sem datasets. Total log-likelihood is unweighted sum.
      - "fixed_deact_sigma": Deactivation sigma fixed at baseline RMSE (0.095091) so trace
        cannot buy likelihood by shrinking its own sigma.
      - "equal_per_dataset": Each dataset's log-likelihood divided by its point count n (1/n),
        so a 100-point trace does not outvote a 1-point scalar 100 to 1. Anchor has n=1.

    If `anchor_po_max` is True, adds the canonical P_o,max gating anchor log-likelihood.
    Does NOT invoke `holdout_guard` because this function is exclusively for consistency
    diagnostics across data sources, not for training models.
    """
    valid_weightings = ("as_is", "fixed_deact_sigma", "equal_per_dataset")
    if weighting not in valid_weightings:
        raise ValueError(f"Unknown weighting {weighting!r}. Must be one of {valid_weightings}")

    datasets = list(datasets)
    fits: list[DatasetFit] = []
    total_ll = 0.0

    for ds in datasets:
        y = _observed(ds)
        pred = _predict(model, ds)
        if pred.shape != y.shape:
            raise ValueError(
                f"{ds.citation_label!r}: prediction shape {pred.shape} does not match "
                f"observation shape {y.shape}"
            )
        n = int(y.size)
        rss = float(np.sum((pred - y) ** 2))
        if not np.isfinite(rss):
            raise ValueError(
                f"{ds.citation_label!r}: residual sum of squares is {rss!r}."
            )
        sigma2_hat = max(rss / n, 1e-12)
        sigma_hat = float(np.sqrt(sigma2_hat))

        noise_model = getattr(ds, "noise_model", "estimated")
        raw_sem = getattr(ds, "sem", getattr(ds, "sem_ms", None))
        if raw_sem is None:
            stated_sem = None
        elif isinstance(raw_sem, (int, float)):
            stated_sem = float(raw_sem)
        else:
            arr = np.asarray(raw_sem, dtype=float)
            stated_sem = float(arr.item()) if arr.size == 1 else tuple(float(x) for x in arr)

        if weighting == "fixed_deact_sigma" and getattr(ds, "observable", None) is Observable.DEACTIVATION:
            # (ii) Deactivation sigma fixed at baseline RMSE (0.095091)
            sigma2 = float(fixed_deact_sigma) ** 2
            ll = -0.5 * (rss / sigma2 + n * np.log(2.0 * np.pi * sigma2))
            noise_model = "fixed_sigma"
            stated_sem = float(fixed_deact_sigma)
        elif noise_model == "stated_sem":
            if raw_sem is None:
                raise ValueError(f"{ds.citation_label!r}: declares noise_model='stated_sem' but has no sem")
            if isinstance(raw_sem, (int, float)):
                sem_val = float(raw_sem)
                if sem_val <= 0.0 or not np.isfinite(sem_val):
                    raise ValueError(f"{ds.citation_label!r}: stated sem must be positive and finite")
                sigma2 = sem_val ** 2
                ll = -0.5 * (rss / sigma2 + n * np.log(2.0 * np.pi * sigma2))
            else:
                sem_arr = np.asarray(raw_sem, dtype=float)
                if sem_arr.size == 1:
                    sem_val = float(sem_arr.item())
                    sigma2 = sem_val ** 2
                    ll = -0.5 * (rss / sigma2 + n * np.log(2.0 * np.pi * sigma2))
                else:
                    sigma2 = sem_arr ** 2
                    ll = -0.5 * float(np.sum((pred - y) ** 2 / sigma2 + np.log(2.0 * np.pi * sigma2)))
        elif noise_model == "estimated":
            sigma2 = sigma2_hat
            ll = -0.5 * n * (np.log(2.0 * np.pi) + 1.0 + np.log(sigma2))
        else:
            raise ValueError(f"{ds.citation_label!r}: unknown noise_model {noise_model!r}")

        if weighting == "equal_per_dataset":
            # (iii) Each dataset contribution divided by point count n
            contrib_ll = ll / n
        else:
            contrib_ll = ll

        fits.append(DatasetFit(
            label=ds.citation_label,
            observable=ds.observable,
            n=n,
            rss=rss,
            sigma_hat=sigma_hat,
            log_likelihood=float(contrib_ll),
            noise_model=noise_model,
            stated_sem=stated_sem,
        ))
        total_ll += contrib_ll

    if anchor_po_max and hasattr(model, "po_max"):
        po_max_val = float(model.po_max)
        ll_anc = po_max_anchor_log_likelihood(po_max_val)
        # po_max gating anchor has n=1, so under equal_per_dataset ll_anc / 1 = ll_anc
        total_ll += float(ll_anc)
        target = FIT_TARGETS["po_max"]
        sigma = (FIT_RANGES["po_max"][1] - FIT_RANGES["po_max"][0]) / 2.0
        rss_anc = float((po_max_val - target) ** 2)
        fits.append(DatasetFit(
            label="po_max gating anchor",
            observable=model.native_observable,
            n=1,
            rss=rss_anc,
            sigma_hat=sigma,
            log_likelihood=float(ll_anc),
            noise_model="stated_sem",
            stated_sem=sigma,
        ))

    return LikelihoodResult(
        total_log_likelihood=total_ll,
        per_dataset=tuple(fits),
        n_estimated=n_params + sum(1 for f in fits if f.noise_model == "estimated"),
    )


def _optimize_diagnostic_start(args: tuple) -> object:
    model_factory, datasets, param_names, bounds, anchor_po_max, weighting, fixed_deact_sigma, s = args
    bnds = list(bounds)
    lo = np.array([b[0] for b in bnds])
    hi = np.array([b[1] for b in bnds])

    def objective(v: np.ndarray) -> float:
        try:
            m = model_factory(np.asarray(v, dtype=float))
            sc = diagnostic_log_likelihood(
                m, datasets, n_params=len(param_names), anchor_po_max=anchor_po_max,
                weighting=weighting, fixed_deact_sigma=fixed_deact_sigma,
            )
            return -float(sc.total_log_likelihood)
        except Exception:
            return 1e12

    # Standard diagnostic tolerances: nesting assertion (rather than tight optimizer tolerances)
    # protects against regression (AGENT-BRIEF-6B §1).
    s_clamped = np.clip(s, lo, hi)
    r = minimize(
        objective, s_clamped, method="L-BFGS-B", bounds=bnds,
        options=dict(ftol=1e-8, gtol=1e-6, maxiter=500),
    )
    start2 = r.x if np.all(np.isfinite(r.x)) else s_clamped
    r2 = minimize(
        objective, start2, method="Nelder-Mead", bounds=bnds,
        options=dict(xatol=1e-4, fatol=1e-4, maxiter=200),
    )
    cand = r2 if (np.isfinite(r2.fun) and r2.fun <= r.fun) else r
    return cand


def fit_diagnostic_consistency(
    model_factory: Callable[[np.ndarray], object],
    datasets: Sequence[object],
    *,
    param_names: Sequence[str],
    bounds: Sequence[tuple[float, float]],
    starts: Sequence[np.ndarray] | None = None,
    anchor_po_max: bool = True,
    n_random_starts: int = 6,
    seed: int = 42,
    weighting: str = "as_is",
    fixed_deact_sigma: float = 0.095091,
    n_workers: int = 1,
) -> DiagnosticFitResult:
    """Find the parameter vector that best satisfies `datasets` jointly.

    Multi-start optimization: L-BFGS-B followed by Nelder-Mead polish.
    Standard diagnostic tolerances: nesting assertion (rather than tight optimizer tolerances)
    protects against regression (AGENT-BRIEF-6B §1).
    """
    datasets = list(datasets)
    names = list(param_names)
    bnds = list(bounds)
    lo = np.array([b[0] for b in bnds])
    hi = np.array([b[1] for b in bnds])
    width = hi - lo

    all_starts = []
    if starts is not None:
        all_starts.extend([np.asarray(s, dtype=float) for s in starts])
    else:
        # Default starts informed by prior knowledge:
        base_3_starts = [
            np.array([3.00, 0.50, 2.50]),
            np.array([2.02, 2.21, -3.00]),
            np.array([3.06, 0.16, 1.41]),
            (lo[:3] + hi[:3]) / 2.0,
        ]
        if len(names) == 3:
            all_starts.extend(base_3_starts)
        elif len(names) in (4, 6):
            from .comparison import MODEL_SPECS, nesting_embedding
            sub_spec = MODEL_SPECS["kinetic_jw95"]
            spec_key = "kinetic_jw95_free_r" if len(names) == 4 else "kinetic_jw95_free_rates"
            big_spec = MODEL_SPECS[spec_key]
            for s3 in base_3_starts:
                all_starts.append(nesting_embedding(sub_spec, big_spec, s3))
        else:
            all_starts.append((lo + hi) / 2.0)

    if n_random_starts > 0:
        sampler = qmc.LatinHypercube(d=len(bnds), seed=seed)
        random_starts = lo + sampler.random(n_random_starts) * width
        all_starts.extend(random_starts)

    task_args = [
        (model_factory, datasets, names, bnds, anchor_po_max, weighting, fixed_deact_sigma, s)
        for s in all_starts
    ]

    best_cand = None
    best_cost = 1e12

    if n_workers > 1 and len(all_starts) > 1:
        import multiprocessing as mp
        ctx = mp.get_context("fork")
        with ctx.Pool(min(n_workers, len(all_starts))) as pool:
            cands = pool.map(_optimize_diagnostic_start, task_args)
    else:
        cands = [_optimize_diagnostic_start(arg) for arg in task_args]

    for cand in cands:
        if np.isfinite(cand.fun) and cand.fun < best_cost:
            best_cost = float(cand.fun)
            best_cand = cand

    if best_cand is None:
        raise RuntimeError("No optimization start converged to a finite loss.")

    best_v = np.asarray(best_cand.x, dtype=float)
    best_m = model_factory(best_v)
    sc = diagnostic_log_likelihood(
        best_m, datasets, n_params=len(names), anchor_po_max=anchor_po_max,
        weighting=weighting, fixed_deact_sigma=fixed_deact_sigma,
    )

    at_bound = tuple(
        names[i]
        for i in range(len(names))
        if (best_v[i] - lo[i] <= 1e-3 * width[i]) or (hi[i] - best_v[i] <= 1e-3 * width[i])
    )

    po_max = float(best_m.po_max) if hasattr(best_m, "po_max") else None
    po_max_in_range = (
        bool(FIT_RANGES["po_max"][0] <= po_max <= FIT_RANGES["po_max"][1])
        if po_max is not None
        else True
    )

    return DiagnosticFitResult(
        values={names[i]: float(best_v[i]) for i in range(len(names))},
        vector=best_v,
        total_log_likelihood=float(sc.total_log_likelihood),
        per_dataset=sc.per_dataset,
        po_max=po_max,
        po_max_in_range=po_max_in_range,
        at_bound=at_bound,
        converged=bool(best_cand.success),
        note=(
            f"{'WARNING: pinned at bound(s) ' + str(at_bound) if at_bound else ''}"
            f"{'; po_max out of range' if not po_max_in_range else ''}"
        ).strip("; "),
    )
