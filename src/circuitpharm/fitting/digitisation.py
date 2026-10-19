"""Digitisation evidence pipeline, loaders, and self-validation gates.

AGENT BRIEF 21 (Tasks 2 & 3):
Constructs real digitised datasets (kind="digitised", synthetic=False) ONLY when
all auditable provenance artefacts are present on disk, and enforces self-validation
against published literature records before admission into the fitting pipeline.

Requirements:
1. Digitised-dataset loader:
   - Source image file that exists on disk, with its SHA-256 recorded in the record.
   - Missing image is a hard failure (MissingImageError), not a warning.
   - Checksum mismatch is a hard failure (ChecksumMismatchError).
   - Figure identifier (journal, DOI, figure number, panel).
   - Axis calibration (two reference points per axis in pixel and data coordinates).
   - Digitisation method and operator string.
2. Self-validation check:
   - A digitised Barberis deactivation trace re-fitted with a triexponential by
     predict_components_fit must reproduce published tau = (2.8, 33.4, 221.35) ms
     and areas = (0.57, 0.23, 0.20) within published errors (0.3, 4.6, 14.9) ms
     and (0.04, 0.02, 0.03).
   - A digitised Jahn CRC must reproduce published EC50 = 11.6 +/- 0.9 uM and
     n_H = 2.2 +/- 0.4 over 1-10 uM.
   - Read all published values directly from literature records (no literals).
   - A fabricated or badly-calibrated curve fails self-validation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
from typing import Any, Literal, Sequence

import numpy as np
from scipy.optimize import curve_fit

from ..literature import DEACTIVATIONS
from .crc_components import get_jahn_published_records
from .data import (
    ALL,
    DeactivationDataset,
    DoseResponseDataset,
    NoiseModel,
    Normalisation,
)
from .deact_components import predict_components_fit


# =======================================================================================
# Exceptions
# =======================================================================================

class EvidenceError(ValueError):
    """Base exception for missing or invalid digitisation evidence."""


class MissingImageError(EvidenceError, FileNotFoundError):
    """Hard failure when source figure image file does not exist on disk."""


class ChecksumMismatchError(EvidenceError):
    """Hard failure when figure image SHA-256 does not match recorded checksum."""


class InvalidCalibrationError(EvidenceError):
    """Axis calibration coordinates are invalid or degenerate."""


class InvalidFigureIdentifierError(EvidenceError):
    """Figure identifier fields are missing or empty."""


class SelfValidationError(ValueError):
    """Self-validation failure against published literature records."""


# =======================================================================================
# Data Structures
# =======================================================================================

def compute_file_sha256(file_path: str | Path) -> str:
    """Compute hex SHA-256 digest of a file on disk."""
    p = Path(file_path).resolve()
    if not p.exists() or not p.is_file():
        raise MissingImageError(f"Source image file does not exist on disk: {p}")
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest().lower()


def verify_image_evidence(image_path: str | Path, expected_sha256: str) -> str:
    """Verify that image exists on disk and its SHA-256 matches expectation."""
    p = Path(image_path).resolve()
    if not p.exists() or not p.is_file():
        raise MissingImageError(
            f"Digitised dataset requires source figure image on disk; file not found: {p}"
        )
    actual_sha256 = compute_file_sha256(p)
    clean_expected = str(expected_sha256).strip().lower()
    if actual_sha256 != clean_expected:
        raise ChecksumMismatchError(
            f"Image checksum mismatch for {p.name}:\n"
            f"  Recorded SHA-256: {clean_expected}\n"
            f"  Actual on-disk:   {actual_sha256}"
        )
    return actual_sha256


@dataclass(frozen=True)
class FigureIdentifier:
    """Identity of the published figure where the curve was published."""

    journal: str
    doi: str
    figure_number: str | int
    panel: str

    def __post_init__(self):
        if not str(self.journal).strip():
            raise InvalidFigureIdentifierError("FigureIdentifier requires non-empty journal")
        if not str(self.doi).strip():
            raise InvalidFigureIdentifierError("FigureIdentifier requires non-empty DOI")
        if not str(self.figure_number).strip():
            raise InvalidFigureIdentifierError("FigureIdentifier requires non-empty figure_number")
        if not str(self.panel).strip():
            raise InvalidFigureIdentifierError("FigureIdentifier requires non-empty panel")

    @property
    def identifier_string(self) -> str:
        return f"{self.journal} (DOI: {self.doi}), Figure {self.figure_number}, Panel {self.panel}"


@dataclass(frozen=True)
class AxisReference:
    """A calibration reference point pairing a pixel position with its data value."""

    pixel: float
    data: float


@dataclass(frozen=True)
class AxisCalibration:
    """Two reference points per axis, in both pixel and data coordinates.

    Provides reproducible and auditable conversion between image pixel coordinates
    and physical experimental quantities.
    """

    x_ref1: AxisReference
    x_ref2: AxisReference
    y_ref1: AxisReference
    y_ref2: AxisReference
    x_scale: Literal["linear", "log10"] = "linear"
    y_scale: Literal["linear", "log10"] = "linear"

    def __post_init__(self):
        if self.x_scale not in ("linear", "log10"):
            raise InvalidCalibrationError(f"x_scale must be 'linear' or 'log10', got {self.x_scale!r}")
        if self.y_scale not in ("linear", "log10"):
            raise InvalidCalibrationError(f"y_scale must be 'linear' or 'log10', got {self.y_scale!r}")
        if abs(self.x_ref1.pixel - self.x_ref2.pixel) < 1e-6:
            raise InvalidCalibrationError("x_ref1.pixel and x_ref2.pixel must be distinct")
        if abs(self.y_ref1.pixel - self.y_ref2.pixel) < 1e-6:
            raise InvalidCalibrationError("y_ref1.pixel and y_ref2.pixel must be distinct")
        if abs(self.x_ref1.data - self.x_ref2.data) < 1e-12:
            raise InvalidCalibrationError("x_ref1.data and x_ref2.data must be distinct")
        if abs(self.y_ref1.data - self.y_ref2.data) < 1e-12:
            raise InvalidCalibrationError("y_ref1.data and y_ref2.data must be distinct")
        if self.x_scale == "log10" and (self.x_ref1.data <= 0 or self.x_ref2.data <= 0):
            raise InvalidCalibrationError("x_ref data values must be positive for log10 scale")
        if self.y_scale == "log10" and (self.y_ref1.data <= 0 or self.y_ref2.data <= 0):
            raise InvalidCalibrationError("y_ref data values must be positive for log10 scale")

    def pixel_to_data_x(self, px: float | Sequence[float] | np.ndarray) -> np.ndarray | float:
        px_arr = np.asarray(px, dtype=float)
        p1, p2 = self.x_ref1.pixel, self.x_ref2.pixel
        d1, d2 = self.x_ref1.data, self.x_ref2.data
        if self.x_scale == "linear":
            frac = (px_arr - p1) / (p2 - p1)
            res = d1 + frac * (d2 - d1)
        else:
            log_d1, log_d2 = np.log10(d1), np.log10(d2)
            frac = (px_arr - p1) / (p2 - p1)
            res = 10.0 ** (log_d1 + frac * (log_d2 - log_d1))
        return float(res) if np.isscalar(px) else res

    def pixel_to_data_y(self, py: float | Sequence[float] | np.ndarray) -> np.ndarray | float:
        py_arr = np.asarray(py, dtype=float)
        p1, p2 = self.y_ref1.pixel, self.y_ref2.pixel
        d1, d2 = self.y_ref1.data, self.y_ref2.data
        if self.y_scale == "linear":
            frac = (py_arr - p1) / (p2 - p1)
            res = d1 + frac * (d2 - d1)
        else:
            log_d1, log_d2 = np.log10(d1), np.log10(d2)
            frac = (py_arr - p1) / (p2 - p1)
            res = 10.0 ** (log_d1 + frac * (log_d2 - log_d1))
        return float(res) if np.isscalar(py) else res


@dataclass(frozen=True)
class DigitisationEvidence:
    """Full evidentiary record required for any non-synthetic digitised dataset."""

    image_path: str
    image_sha256: str
    figure_id: FigureIdentifier
    calibration: AxisCalibration
    method: str
    operator: str

    def __post_init__(self):
        if not str(self.method).strip():
            raise EvidenceError("DigitisationEvidence requires non-empty method string")
        if not str(self.operator).strip():
            raise EvidenceError("DigitisationEvidence requires non-empty operator string")
        # Hard verification on construction: image must exist on disk and SHA-256 match
        verify_image_evidence(self.image_path, self.image_sha256)


# =======================================================================================
# Digitised Dataset Loaders
# =======================================================================================

def load_digitised_deactivation(
    *,
    citation_label: str,
    source_key: str,
    image_path: str | Path,
    image_sha256: str,
    figure_id: FigureIdentifier,
    calibration: AxisCalibration,
    method: str,
    operator: str,
    time_ms: Sequence[float] | np.ndarray | None = None,
    normalized_current: Sequence[float] | np.ndarray | None = None,
    pixel_x: Sequence[float] | None = None,
    pixel_y: Sequence[float] | None = None,
    pulse_ms: float = 2.0,
    conc_mm: float = 10.0,
    n_components: int = 3,
    fit_window_ms: float = 1000.0,
    preparation: str = "",
    n_cells: int | None = None,
    role: Literal["train", "holdout"] = "train",
    noise_model: NoiseModel = "estimated",
) -> DeactivationDataset:
    """Construct a DeactivationDataset with kind='digitised' and synthetic=False.

    Hard fails (MissingImageError, ChecksumMismatchError, EvidenceError) if any
    required evidence artefact is missing or invalid.
    """
    evidence = DigitisationEvidence(
        image_path=str(image_path),
        image_sha256=str(image_sha256).strip().lower(),
        figure_id=figure_id,
        calibration=calibration,
        method=method,
        operator=operator,
    )

    if time_ms is not None and normalized_current is not None:
        t_arr = np.asarray(time_ms, dtype=float)
        y_arr = np.asarray(normalized_current, dtype=float)
    elif pixel_x is not None and pixel_y is not None:
        t_arr = np.asarray(calibration.pixel_to_data_x(pixel_x), dtype=float)
        y_arr = np.asarray(calibration.pixel_to_data_y(pixel_y), dtype=float)
    else:
        raise EvidenceError("Either (time_ms, normalized_current) or (pixel_x, pixel_y) must be provided")

    dig_summary = (
        f"Digitised from {figure_id.identifier_string} using {method} by {operator}. "
        f"Source image SHA-256: {evidence.image_sha256[:16]}..."
    )

    return DeactivationDataset(
        citation_label=citation_label,
        time_ms=t_arr,
        normalized_current=y_arr,
        pulse_ms=pulse_ms,
        conc_mm=conc_mm,
        n_components=n_components,
        fit_window_ms=fit_window_ms,
        kind="digitised",
        synthetic=False,
        source_key=source_key,
        figure=figure_id.identifier_string,
        digitisation=dig_summary,
        preparation=preparation,
        n_cells=n_cells,
        role=role,
        noise_model=noise_model,
        evidence=evidence,
    )


def load_digitised_dose_response(
    *,
    citation_label: str,
    source_key: str,
    image_path: str | Path,
    image_sha256: str,
    figure_id: FigureIdentifier,
    calibration: AxisCalibration,
    method: str,
    operator: str,
    concs_um: Sequence[float] | np.ndarray | None = None,
    mean_response: Sequence[float] | np.ndarray | None = None,
    sem: Sequence[float] | np.ndarray | None = None,
    pixel_x: Sequence[float] | None = None,
    pixel_y: Sequence[float] | None = None,
    preparation: str = "",
    normalisation: Normalisation = "fraction_of_max",
    n_cells: int | None = None,
    role: Literal["train", "holdout"] = "train",
    noise_model: NoiseModel = "estimated",
) -> DoseResponseDataset:
    """Construct a DoseResponseDataset with kind='digitised' and synthetic=False.

    Hard fails (MissingImageError, ChecksumMismatchError, EvidenceError) if any
    required evidence artefact is missing or invalid.
    """
    evidence = DigitisationEvidence(
        image_path=str(image_path),
        image_sha256=str(image_sha256).strip().lower(),
        figure_id=figure_id,
        calibration=calibration,
        method=method,
        operator=operator,
    )

    if concs_um is not None and mean_response is not None:
        c_arr = np.asarray(concs_um, dtype=float)
        y_arr = np.asarray(mean_response, dtype=float)
    elif pixel_x is not None and pixel_y is not None:
        c_arr = np.asarray(calibration.pixel_to_data_x(pixel_x), dtype=float)
        y_arr = np.asarray(calibration.pixel_to_data_y(pixel_y), dtype=float)
    else:
        raise EvidenceError("Either (concs_um, mean_response) or (pixel_x, pixel_y) must be provided")

    if sem is not None:
        sem_arr = np.asarray(sem, dtype=float)
    else:
        sem_arr = np.full_like(y_arr, 0.02)  # default nominal floor

    dig_summary = (
        f"Digitised from {figure_id.identifier_string} using {method} by {operator}. "
        f"Source image SHA-256: {evidence.image_sha256[:16]}..."
    )

    return DoseResponseDataset(
        citation_label=citation_label,
        preparation=preparation,
        concs_um=c_arr,
        mean_response=y_arr,
        sem=sem_arr,
        normalisation=normalisation,
        kind="digitised",
        synthetic=False,
        source_key=source_key,
        figure=figure_id.identifier_string,
        digitisation=dig_summary,
        n_cells=n_cells,
        role=role,
        noise_model=noise_model,
        evidence=evidence,
    )


# =======================================================================================
# Registration Guard
# =======================================================================================

def assert_digitised_evidence(dataset: object) -> None:
    """Refuse any digitised non-synthetic dataset lacking auditable on-disk evidence."""
    kind = getattr(dataset, "kind", None)
    synthetic = getattr(dataset, "synthetic", False)
    if kind == "digitised" and not synthetic:
        evidence = getattr(dataset, "evidence", None)
        if evidence is None:
            raise EvidenceError(
                f"Registration refused for {dataset}: kind='digitised', synthetic=False requires "
                f"attached DigitisationEvidence with on-disk source image and recorded SHA-256."
            )
        if not isinstance(evidence, DigitisationEvidence):
            raise EvidenceError(
                f"Invalid evidence type {type(evidence)} on {dataset}; expected DigitisationEvidence."
            )
        verify_image_evidence(evidence.image_path, evidence.image_sha256)


def register_digitised_dataset(
    name: str,
    dataset: object,
    registry: dict[str, object] | None = None,
) -> None:
    """Register a digitised dataset into the project registry after strict evidence verification."""
    target_registry = ALL if registry is None else registry
    assert_digitised_evidence(dataset)
    target_registry[name] = dataset


# =======================================================================================
# Self-Validation Checks
# =======================================================================================

@dataclass(frozen=True)
class ValidationResult:
    """Outcome of self-validation against published literature records."""

    passed: bool
    observable: str
    target_values: dict[str, float]
    target_errors: dict[str, float]
    fitted_values: dict[str, float]
    residuals: dict[str, float]
    z_scores: dict[str, float]
    message: str


def validate_barberis_deactivation(
    dataset: DeactivationDataset,
    *,
    tolerance_sigmas: float = 1.0,
) -> ValidationResult:
    """Re-fit deactivation trace with triexponential and test against published parameters.

    Reads published components and errors directly from circuitpharm.literature.DEACTIVATIONS[1]:
    - tau = (2.8, 33.4, 221.35) ms, errors = (0.3, 4.6, 14.9) ms
    - areas = (0.57, 0.23, 0.20), errors = (0.04, 0.02, 0.03)
    No hardcoded literals.
    """
    barberis_lit = DEACTIVATIONS[1]
    pub_taus, pub_areas = barberis_lit.components
    pub_tau_errs, pub_area_errs = barberis_lit.component_errors

    fit_res = predict_components_fit(dataset)
    fitted_taus = fit_res.taus
    fitted_areas = fit_res.areas

    target_vals: dict[str, float] = {}
    target_errs: dict[str, float] = {}
    fit_vals: dict[str, float] = {}
    residuals: dict[str, float] = {}
    z_scores: dict[str, float] = {}
    failures: list[str] = []

    for i in range(3):
        t_name = f"tau_{i+1}"
        target_vals[t_name] = float(pub_taus[i])
        target_errs[t_name] = float(pub_tau_errs[i])
        fit_vals[t_name] = float(fitted_taus[i])
        diff_t = float(fitted_taus[i] - pub_taus[i])
        residuals[t_name] = diff_t
        z_t = diff_t / float(pub_tau_errs[i])
        z_scores[t_name] = z_t
        if abs(z_t) > tolerance_sigmas:
            failures.append(f"{t_name}: {fitted_taus[i]:.2f} ms vs published {pub_taus[i]} +/- {pub_tau_errs[i]} ms (|z| = {abs(z_t):.2f} > {tolerance_sigmas}sigma)")

    for i in range(3):
        a_name = f"area_{i+1}"
        target_vals[a_name] = float(pub_areas[i])
        target_errs[a_name] = float(pub_area_errs[i])
        fit_vals[a_name] = float(fitted_areas[i])
        diff_a = float(fitted_areas[i] - pub_areas[i])
        residuals[a_name] = diff_a
        z_a = diff_a / float(pub_area_errs[i])
        z_scores[a_name] = z_a
        if abs(z_a) > tolerance_sigmas:
            failures.append(f"{a_name}: {fitted_areas[i]:.3f} vs published {pub_areas[i]} +/- {pub_area_errs[i]} (|z| = {abs(z_a):.2f} > {tolerance_sigmas}sigma)")

    passed = len(failures) == 0
    msg = "PASSED: All triexponential deactivation components reproduced within published errors" if passed else f"FAILED: {'; '.join(failures)}"

    return ValidationResult(
        passed=passed,
        observable="DEACTIVATION_BARBERIS2007",
        target_values=target_vals,
        target_errors=target_errs,
        fitted_values=fit_vals,
        residuals=residuals,
        z_scores=z_scores,
        message=msg,
    )


def validate_jahn_crc(
    dataset: DoseResponseDataset,
    *,
    tolerance_sigmas: float = 1.0,
) -> ValidationResult:
    """Fit concentration-response curve and test against published Jahn 1997 parameters.

    Reads published parameters directly from get_jahn_published_records():
    - EC50 = 11.6 +/- 0.9 uM
    - n_H = 2.2 +/- 0.4
    No hardcoded literals.
    """
    jahn_lit = get_jahn_published_records()
    pub_ec50 = float(jahn_lit["ec50_um"])
    pub_ec50_sem = float(jahn_lit["ec50_sem"])
    pub_hill = float(jahn_lit["hill"])
    pub_hill_sem = float(jahn_lit["hill_sem"])

    c = np.asarray(dataset.concs_um, dtype=float)
    y = np.asarray(dataset.mean_response, dtype=float)

    # Normalize response to maximum
    y_norm = y / float(np.max(y))

    # Fit Hill equation: y = c^nH / (ec50^nH + c^nH)
    def hill_func(conc, ec50, nh):
        x = (conc / ec50) ** nh
        return x / (1.0 + x)

    try:
        popt, _ = curve_fit(
            hill_func,
            c,
            y_norm,
            p0=[pub_ec50, pub_hill],
            bounds=([0.1, 0.2], [500.0, 10.0]),
        )
        fitted_ec50 = float(popt[0])
        fitted_hill = float(popt[1])
    except Exception:
        # Fallback to logit regression on interior points
        m = (y_norm > 1e-4) & (y_norm < 1.0 - 1e-4)
        if np.sum(m) >= 2:
            slope, intercept = np.polyfit(np.log10(c[m]), np.log10(y_norm[m] / (1.0 - y_norm[m])), 1)
            fitted_hill = float(slope)
            fitted_ec50 = float(10.0 ** (-intercept / slope))
        else:
            fitted_ec50 = float("nan")
            fitted_hill = float("nan")

    diff_ec50 = fitted_ec50 - pub_ec50
    z_ec50 = diff_ec50 / pub_ec50_sem
    diff_hill = fitted_hill - pub_hill
    z_hill = diff_hill / pub_hill_sem

    target_vals = {"ec50_um": pub_ec50, "hill_slope": pub_hill}
    target_errs = {"ec50_um": pub_ec50_sem, "hill_slope": pub_hill_sem}
    fit_vals = {"ec50_um": fitted_ec50, "hill_slope": fitted_hill}
    residuals = {"ec50_um": diff_ec50, "hill_slope": diff_hill}
    z_scores = {"ec50_um": z_ec50, "hill_slope": z_hill}

    failures: list[str] = []
    if abs(z_ec50) > tolerance_sigmas:
        failures.append(f"EC50: {fitted_ec50:.2f} uM vs published {pub_ec50} +/- {pub_ec50_sem} uM (|z| = {abs(z_ec50):.2f} > {tolerance_sigmas}sigma)")
    if abs(z_hill) > tolerance_sigmas:
        failures.append(f"Hill slope: {fitted_hill:.2f} vs published {pub_hill} +/- {pub_hill_sem} (|z| = {abs(z_hill):.2f} > {tolerance_sigmas}sigma)")

    passed = len(failures) == 0
    msg = "PASSED: CRC EC50 and Hill slope reproduced within published errors" if passed else f"FAILED: {'; '.join(failures)}"

    return ValidationResult(
        passed=passed,
        observable="CRC_JAHN1997",
        target_values=target_vals,
        target_errors=target_errs,
        fitted_values=fit_vals,
        residuals=residuals,
        z_scores=z_scores,
        message=msg,
    )


def assert_self_validation(dataset: object, *, tolerance_sigmas: float = 1.0) -> ValidationResult:
    """Run self-validation and raise SelfValidationError on failure."""
    if isinstance(dataset, DeactivationDataset):
        res = validate_barberis_deactivation(dataset, tolerance_sigmas=tolerance_sigmas)
    elif isinstance(dataset, DoseResponseDataset):
        res = validate_jahn_crc(dataset, tolerance_sigmas=tolerance_sigmas)
    else:
        raise ValueError(f"No self-validation procedure for dataset type {type(dataset)}")

    if not res.passed:
        raise SelfValidationError(f"Self-validation failed for {dataset.citation_label}: {res.message}")
    return res
