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
from ..gabaa_kinetics import (
    FIT_RANGES,
    FIT_TARGETS,
    po_max_anchor_log_likelihood,
)


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
) -> LikelihoodResult:
    """Evaluate joint log-likelihood over arbitrary datasets for consistency diagnostics.

    Each dataset is scored according to its own declared `noise_model`:
      - "stated_sem": Gaussian log-likelihood using the dataset's declared SEM.
      - "estimated": concentrated Gaussian log-likelihood with sigma profiled out.

    If `anchor_po_max` is True, adds the canonical P_o,max gating anchor log-likelihood.
    Does NOT invoke `holdout_guard` because this function is exclusively for consistency
    diagnostics across data sources, not for training models.
    """
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

        if noise_model == "stated_sem":
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

        fits.append(DatasetFit(
            label=ds.citation_label,
            observable=ds.observable,
            n=n,
            rss=rss,
            sigma_hat=sigma_hat,
            log_likelihood=float(ll),
            noise_model=noise_model,
            stated_sem=stated_sem,
        ))
        total_ll += ll

    if anchor_po_max and hasattr(model, "po_max"):
        po_max_val = float(model.po_max)
        ll_anc = po_max_anchor_log_likelihood(po_max_val)
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


def fit_diagnostic_consistency(
    model_factory: Callable[[np.ndarray], object],
    datasets: Sequence[object],
    *,
    param_names: Sequence[str],
    bounds: Sequence[tuple[float, float]],
    starts: Sequence[np.ndarray] | None = None,
    anchor_po_max: bool = True,
    n_random_starts: int = 0,
    seed: int = 42,
) -> DiagnosticFitResult:
    """Find the parameter vector that best satisfies `datasets` jointly.

    Cost-model-conscious multi-start optimization: L-BFGS-B followed by short Nelder-Mead polish.
    """
    datasets = list(datasets)
    names = list(param_names)
    bnds = list(bounds)
    lo = np.array([b[0] for b in bnds])
    hi = np.array([b[1] for b in bnds])
    width = hi - lo

    def objective(v: np.ndarray) -> float:
        try:
            m = model_factory(np.asarray(v, dtype=float))
            sc = diagnostic_log_likelihood(m, datasets, n_params=len(names), anchor_po_max=anchor_po_max)
            return -float(sc.total_log_likelihood)
        except Exception:
            return 1e12

    all_starts = []
    if starts is not None:
        all_starts.extend([np.asarray(s, dtype=float) for s in starts])
    else:
        # Default starts informed by prior knowledge:
        if len(names) == 3:
            # 1. Known grid point satisfying holdout
            all_starts.append(np.array([3.00, 0.50, 2.50]))
            # 2. Known CRC alone optimum
            all_starts.append(np.array([2.02, 2.21, -3.00]))
            # 3. Known Deact alone optimum
            all_starts.append(np.array([3.06, 0.16, 1.41]))
            # 4. Box center
            all_starts.append((lo + hi) / 2.0)
        elif len(names) == 4:
            # log10_kd, log10_E, log10_D, log10_r
            r_conv = float(np.log10(0.002))
            all_starts.append(np.array([3.00, 0.50, 2.50, r_conv]))
            all_starts.append(np.array([2.02, 2.21, -3.00, r_conv]))
            all_starts.append(np.array([3.06, 0.16, 1.41, r_conv]))
            all_starts.append((lo + hi) / 2.0)
        elif len(names) == 6:
            # log10_kd, log10_E, log10_D, log10_r, log10_koff, log10_alpha
            r_conv = float(np.log10(0.002))
            koff_conv = float(np.log10(0.18))
            alpha_conv = float(np.log10(0.36))
            all_starts.append(np.array([3.00, 0.50, 2.50, r_conv, koff_conv, alpha_conv]))
            all_starts.append(np.array([2.02, 2.21, -3.00, r_conv, koff_conv, alpha_conv]))
            all_starts.append(np.array([3.06, 0.16, 1.41, r_conv, koff_conv, alpha_conv]))
            all_starts.append((lo + hi) / 2.0)
        else:
            all_starts.append((lo + hi) / 2.0)

    if n_random_starts > 0:
        rng = np.random.default_rng(seed)
        for _ in range(n_random_starts):
            all_starts.append(lo + rng.uniform(0.0, 1.0, len(bnds)) * width)

    best_cand = None
    best_cost = 1e12

    for s in all_starts:
        s_clamped = np.clip(s, lo, hi)
        r = minimize(
            objective, s_clamped, method="L-BFGS-B", bounds=bnds,
            options=dict(ftol=1e-8, gtol=1e-6, maxiter=500)
        )
        start2 = r.x if np.all(np.isfinite(r.x)) else s_clamped
        r2 = minimize(
            objective, start2, method="Nelder-Mead", bounds=bnds,
            options=dict(xatol=1e-4, fatol=1e-4, maxiter=200)
        )
        cand = r2 if (np.isfinite(r2.fun) and r2.fun <= r.fun) else r
        if np.isfinite(cand.fun) and cand.fun < best_cost:
            best_cost = float(cand.fun)
            best_cand = cand

    if best_cand is None:
        raise RuntimeError("No optimization start converged to a finite loss.")

    best_v = np.asarray(best_cand.x, dtype=float)
    best_m = model_factory(best_v)
    sc = diagnostic_log_likelihood(best_m, datasets, n_params=len(names), anchor_po_max=anchor_po_max)

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
