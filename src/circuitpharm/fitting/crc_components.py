"""CRC component prediction, extraction, and likelihood scoring (Weighting iv for CRC).

AGENT BRIEF 12 (Task 1):
Weighting (iv) declared objective extended from deactivation to every parametric dataset
per DECISION-deactivation-objective.md amendment (2026-10-17).

The Jahn 1997 peak concentration-response curve (JAHN1997_PEAK_CRC) is a parametric
reconstruction generated from two published numbers:
    EC50 = 11.6 +/- 0.9 uM
    nH   = 2.2 +/- 0.4 (measured over 1-10 uM)

Scoring the model against nine generated pseudo-observations with profiled sigma
manufactured artificial leverage and biased joint fits against deactivation (which was
scored against its six published parameters). Under Weighting (iv), the model is scored
directly against the two published numbers at their published errors:
    z_ec50 = (ec50_model - 11.6) / 0.9
    z_hill = (nH_model - 2.2) / 0.4
    chi2_crc = z_ec50^2 + z_hill^2

DISCREPANCIES AND CONVENTIONS HANDLED EXPLICITLY:
1. Discrepancy between FIT_TARGETS['ec50_um'] (20.0 uM, recalled literature middle) and
   JAHN_EC50_UM (11.6 uM, sourced measurement). Both numbers are reported and surfaced,
   unresolved.
2. Hill slope extraction window: The source reports "The slope between 0.001 and 0.01 mM
   GABA was 2.2 +/- 0.4", i.e. an absolute window of 1-10 uM. The slope is extracted over
   this measured window by default. Sensitivity checks over the full curve (2-98% band)
   and relative window (JAHN1997_SLOPE_WINDOW_REL) are also computed and reported.
3. Normalisation convention: The model's peak response is normalised to its own maximum
   (fraction_of_max, y / y_max), matching the experimental I / I_max protocol. Absolute
   open probability is separately constrained by the P_o,max anchor [0.70, 0.80].
4. The 9-point objective (RMSE and ratio to baseline 0.004917) is retained as a reported
   view and never deleted.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

from ..gabaa_kinetics import FIT_TARGETS
from ..models.base import (
    JAHN1997_SLOPE_WINDOW_REL,
    PEAK_APPLICATION_MS,
    ReceptorModel,
    hill_slope,
)
from .data import (
    JAHN1997_PEAK_CRC,
    JAHN1997_PEAK_CRC_RESTRICTED,
    JAHN_EC50_SEM,
    JAHN_EC50_UM,
    JAHN_HILL,
    JAHN_HILL_SEM,
    JAHN_SLOPE_WINDOW_UM,
)

#: Baseline RMSE of the 9-point parametric reconstruction (reported view)
CRC_BASELINE_RMSE = 0.004917


@dataclass(frozen=True)
class CrcParameters:
    """Extracted concentration-response curve parameters from a model prediction.

    Attributes:
        ec50_um: Concentration at half-maximal normalized peak response (uM).
        hill_measured: Hill coefficient fitted over the measured 1-10 uM window.
        hill_rel: Hill coefficient fitted over the relative rising-phase window.
        hill_full: Hill coefficient fitted over the full 2-98% response band (sensitivity check).
        y_max: Maximum peak open probability attained at saturation.
        normalisation: Normalisation convention applied ('fraction_of_max').
    """

    ec50_um: float
    hill_measured: float
    hill_rel: float
    hill_full: float
    y_max: float
    normalisation: str = "fraction_of_max"


@dataclass(frozen=True)
class CrcComponentsLogLikelihood:
    """Scoring result of model CRC parameters against published targets at published errors.

    Attributes:
        ec50_um: Model predicted EC50 in uM.
        hill_slope: Model predicted Hill slope over measured window (1-10 uM).
        z_ec50: Residual in sigma for EC50 against published 11.6 +/- 0.9 uM.
        z_hill: Residual in sigma for Hill slope against published 2.2 +/- 0.4.
        chi2_ec50: z_ec50^2.
        chi2_hill: z_hill^2.
        chi2_total: chi2_ec50 + chi2_hill.
        log_likelihood: Gaussian log-likelihood term.
        max_z: Maximum absolute residual in sigma: max(|z_ec50|, |z_hill|).
        hill_slope_full: Full-range Hill slope sensitivity check.
        z_hill_full: Residual in sigma for full-range Hill slope.
        hill_slope_rel: Relative-window Hill slope sensitivity check.
        z_hill_rel: Residual in sigma for relative-window Hill slope.
        crc_rmse: 9-point reconstruction RMSE (reported view).
        crc_ratio: 9-point reconstruction RMSE ratio to baseline 0.004917 (reported view).
    """

    ec50_um: float
    hill_slope: float
    z_ec50: float
    z_hill: float
    chi2_ec50: float
    chi2_hill: float
    chi2_total: float
    log_likelihood: float
    max_z: float
    hill_slope_full: float
    z_hill_full: float
    hill_slope_rel: float
    z_hill_rel: float
    crc_rmse: float
    crc_ratio: float


def get_jahn_published_records() -> dict[str, Any]:
    """Return published numbers and errors read directly from project records (no literals)."""
    return {
        "ec50_um": float(JAHN_EC50_UM),
        "ec50_sem": float(JAHN_EC50_SEM),
        "hill": float(JAHN_HILL),
        "hill_sem": float(JAHN_HILL_SEM),
        "measured_window_um": tuple(float(x) for x in JAHN_SLOPE_WINDOW_UM),
        "source": "Jahn et al. 1997 NeuroReport 8(16):3443-6 via circuitpharm.fitting.data",
        "dataset_key": JAHN1997_PEAK_CRC.citation_label,
        "restricted_dataset_key": JAHN1997_PEAK_CRC_RESTRICTED.citation_label,
    }


def surface_ec50_discrepancy() -> dict[str, Any]:
    """Surface the discrepancy between FIT_TARGETS['ec50_um'] and the sourced value.

    Per AGENT-BRIEF-12 Task 1 item 1:
    FIT_TARGETS['ec50_um'] is 20.0 uM (the middle of a recalled 10-30 uM literature range).
    The one sourced value for alpha1beta2gamma2L in HEK293 is 11.6 +/- 0.9 uM (Jahn 1997).
    Both numbers are reported, unresolved.
    """
    recalled_target = float(FIT_TARGETS["ec50_um"])
    sourced_val = float(JAHN_EC50_UM)
    sourced_sem = float(JAHN_EC50_SEM)
    diff_um = recalled_target - sourced_val
    diff_sigma = diff_um / sourced_sem
    return {
        "fit_targets_ec50_um": recalled_target,
        "fit_targets_source": "circuitpharm.gabaa_kinetics.FIT_TARGETS['ec50_um'] (recalled literature middle)",
        "jahn_sourced_ec50_um": sourced_val,
        "jahn_sourced_ec50_sem": sourced_sem,
        "jahn_source": "Jahn et al. 1997 (NeuroReport 8(16):3443-6) via circuitpharm.fitting.data",
        "discrepancy_um": diff_um,
        "discrepancy_sigma": diff_sigma,
        "resolution": "unresolved (reported explicitly, neither substituted quietly)",
    }


def extract_crc_parameters(
    model: ReceptorModel,
    *,
    concs_um: np.ndarray | None = None,
    slope_window_um: tuple[float, float] = JAHN_SLOPE_WINDOW_UM,
    full_range_band: tuple[float, float] = (0.02, 0.98),
    n_points_curve: int = 400,
    application_ms: float = PEAK_APPLICATION_MS,
) -> CrcParameters:
    """Extract EC50 and Hill slopes from model peak concentration-response curve.

    Normalisation:
        The peak open probability curve across concentrations is divided by its own
        maximum response (y_norm = y / y_max), matching the experimental I / I_max
        protocol. EC50 is extracted by linear interpolation in log-concentration space
        at y_norm = 0.5.

    Hill Slope:
        Fitted over the measured window (slope_window_um, default 1-10 uM) by logit-space
        linear regression. Sensitivity checks over the full range (full_range_band) and
        the relative window (JAHN1997_SLOPE_WINDOW_REL) are also extracted.
    """
    if concs_um is None:
        c = np.logspace(-2.0, 4.5, n_points_curve)
    else:
        c = np.asarray(concs_um, dtype=float)

    # Simulate peak dose response
    y = np.asarray(model.peak_dose_response(c, application_ms=application_ms), dtype=float)
    y_max = float(np.max(y))
    if not (np.isfinite(y_max) and y_max > 0.0):
        raise ValueError(f"Model {model.name} returned non-positive peak open probability: {y_max!r}")

    y_norm = y / y_max

    # Extract EC50 via interpolation on log-concentration grid
    log_c = np.log10(c)
    if y_norm.min() >= 0.5 or y_norm.max() <= 0.5:
        # Curve does not span the half-maximum
        log_ec50 = float("nan")
        ec50_um = float("nan")
    else:
        log_ec50 = float(np.interp(0.5, y_norm, log_c))
        ec50_um = float(10.0 ** log_ec50)

    # 1. Slope over measured window (1-10 uM)
    try:
        hill_meas = hill_slope(c, y, window_um=slope_window_um)
    except Exception:
        hill_meas = float("nan")

    # 2. Slope over relative window
    try:
        hill_rel = hill_slope(c, y, window_rel=JAHN1997_SLOPE_WINDOW_REL)
    except Exception:
        hill_rel = float("nan")

    # 3. Slope over full range (sensitivity check)
    try:
        hill_full = hill_slope(c, y, band=full_range_band)
    except Exception:
        hill_full = float("nan")

    return CrcParameters(
        ec50_um=ec50_um,
        hill_measured=hill_meas,
        hill_rel=hill_rel,
        hill_full=hill_full,
        y_max=y_max,
        normalisation="fraction_of_max",
    )


def crc_components_log_likelihood(
    params_or_model: CrcParameters | ReceptorModel,
    *,
    target_ec50: float = JAHN_EC50_UM,
    sem_ec50: float = JAHN_EC50_SEM,
    target_hill: float = JAHN_HILL,
    sem_hill: float = JAHN_HILL_SEM,
    compute_9point_view: bool = True,
) -> CrcComponentsLogLikelihood:
    """Compute Gaussian log-likelihood and residuals in sigma for CRC parameters.

    Scores model EC50 against target_ec50 +/- sem_ec50 and Hill slope against
    target_hill +/- sem_hill (Weighting iv).
    """
    if isinstance(params_or_model, CrcParameters):
        p = params_or_model
        model = None
    else:
        model = params_or_model
        p = extract_crc_parameters(model)

    ec50 = p.ec50_um
    hill_val = p.hill_measured

    z_ec50 = float((ec50 - target_ec50) / sem_ec50) if np.isfinite(ec50) else 999.0
    z_hill = float((hill_val - target_hill) / sem_hill) if np.isfinite(hill_val) else 999.0

    chi2_ec50 = float(z_ec50 ** 2)
    chi2_hill = float(z_hill ** 2)
    chi2_total = chi2_ec50 + chi2_hill

    # Uncorrelated bivariate Gaussian log-likelihood
    # log L = -0.5 * chi2 - log(2 * pi * sigma_1 * sigma_2)
    ll = float(-0.5 * chi2_total - np.log(2.0 * np.pi * sem_ec50 * sem_hill))

    max_z = float(max(abs(z_ec50), abs(z_hill)))

    z_hill_full = float((p.hill_full - target_hill) / sem_hill) if np.isfinite(p.hill_full) else 999.0
    z_hill_rel = float((p.hill_rel - target_hill) / sem_hill) if np.isfinite(p.hill_rel) else 999.0

    # 9-point reported view
    if compute_9point_view and model is not None:
        try:
            from .likelihood import _predict
            pred_9 = _predict(model, JAHN1997_PEAK_CRC)
            obs_9 = np.asarray(JAHN1997_PEAK_CRC.mean_response, dtype=float)
            crc_rmse = float(np.sqrt(np.mean((pred_9 - obs_9) ** 2)))
            crc_ratio = float(crc_rmse / CRC_BASELINE_RMSE)
        except Exception:
            crc_rmse = float("nan")
            crc_ratio = float("nan")
    else:
        crc_rmse = float("nan")
        crc_ratio = float("nan")

    return CrcComponentsLogLikelihood(
        ec50_um=ec50,
        hill_slope=hill_val,
        z_ec50=z_ec50,
        z_hill=z_hill,
        chi2_ec50=chi2_ec50,
        chi2_hill=chi2_hill,
        chi2_total=chi2_total,
        log_likelihood=ll,
        max_z=max_z,
        hill_slope_full=p.hill_full,
        z_hill_full=z_hill_full,
        hill_slope_rel=p.hill_rel,
        z_hill_rel=z_hill_rel,
        crc_rmse=crc_rmse,
        crc_ratio=crc_ratio,
    )
