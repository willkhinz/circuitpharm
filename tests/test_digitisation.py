"""Tests for digitisation evidence verification, loaders, and self-validation.

AGENT BRIEF 21 (Task 2 & 3):
Verifies that:
1. The loader refuses without an image (hard failure MissingImageError).
2. The loader refuses on checksum mismatch (hard failure ChecksumMismatchError).
3. The loader accepts a correctly-formed record with on-disk image, matching SHA-256,
   figure identifier, axis calibration, and method/operator.
4. Registration is impossible without complete evidence.
5. Axis calibration accurately transforms pixel coordinates to data units.
6. Self-validation verifies that unperturbed curves reproduce published parameters
   within their published errors (read from literature records, no literals).
7. Self-validation catches deliberately perturbed curves (test fixtures, not data).
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import numpy as np
import pytest

from circuitpharm.fitting import (
    ALL,
    BARBERIS2007_DEACT,
    JAHN1997_PEAK_CRC,
    AxisCalibration,
    AxisReference,
    ChecksumMismatchError,
    DigitisationEvidence,
    EvidenceError,
    FigureIdentifier,
    InvalidCalibrationError,
    InvalidFigureIdentifierError,
    MissingImageError,
    SelfValidationError,
    ValidationResult,
    assert_digitised_evidence,
    assert_self_validation,
    compute_file_sha256,
    load_digitised_deactivation,
    load_digitised_dose_response,
    register_digitised_dataset,
    validate_barberis_deactivation,
    validate_jahn_crc,
    verify_image_evidence,
)
from circuitpharm.fitting.crc_components import get_jahn_published_records
from circuitpharm.literature import DEACTIVATIONS


# =======================================================================================
# Test Fixtures (explicitly labelled as test fixtures, NOT real data)
# =======================================================================================

@pytest.fixture
def dummy_image_file(tmp_path: Path) -> tuple[Path, str]:
    """Create a temporary image file on disk and return its path and SHA-256."""
    img_path = tmp_path / "barberis_fig1.png"
    # Write a deterministic sequence of bytes simulating an image file
    content = b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x10\x00\x00\x00\x10\x08\x06\x00\x00\x00\x1f\xf3\xffa"
    img_path.write_bytes(content)
    sha256 = hashlib.sha256(content).hexdigest()
    return img_path, sha256


@pytest.fixture
def valid_figure_id() -> FigureIdentifier:
    return FigureIdentifier(
        journal="Eur J Neurosci",
        doi="10.1111/j.1460-9568.2007.05530.x",
        figure_number=1,
        panel="C",
    )


@pytest.fixture
def valid_linear_calibration() -> AxisCalibration:
    return AxisCalibration(
        x_ref1=AxisReference(pixel=100.0, data=0.0),
        x_ref2=AxisReference(pixel=600.0, data=1000.0),
        y_ref1=AxisReference(pixel=500.0, data=0.0),
        y_ref2=AxisReference(pixel=100.0, data=1.0),
        x_scale="linear",
        y_scale="linear",
    )


@pytest.fixture
def valid_log_calibration() -> AxisCalibration:
    return AxisCalibration(
        x_ref1=AxisReference(pixel=50.0, data=0.1),
        x_ref2=AxisReference(pixel=650.0, data=1000.0),
        y_ref1=AxisReference(pixel=500.0, data=0.0),
        y_ref2=AxisReference(pixel=100.0, data=1.0),
        x_scale="log10",
        y_scale="linear",
    )


# =======================================================================================
# 1. Loader Evidence & Checksum Verification Tests
# =======================================================================================

def test_loader_refuses_without_image_file(tmp_path: Path, valid_figure_id, valid_linear_calibration):
    """Hard failure (MissingImageError) if image file does not exist on disk."""
    non_existent = tmp_path / "ghost_figure.png"
    assert not non_existent.exists()

    with pytest.raises(MissingImageError, match="source figure image on disk"):
        load_digitised_deactivation(
            citation_label="test deact",
            source_key="barberis2007",
            image_path=non_existent,
            image_sha256="0" * 64,
            figure_id=valid_figure_id,
            calibration=valid_linear_calibration,
            method="WebPlotDigitizer",
            operator="agent-brief-21",
            time_ms=np.linspace(0, 1000, 50),
            normalized_current=np.exp(-np.linspace(0, 1000, 50) / 50.0),
        )


def test_loader_refuses_on_checksum_mismatch(dummy_image_file, valid_figure_id, valid_linear_calibration):
    """Hard failure (ChecksumMismatchError) if on-disk SHA-256 does not match recorded checksum."""
    img_path, true_sha256 = dummy_image_file
    tampered_sha256 = "f" * 64  # deliberate mismatch

    with pytest.raises(ChecksumMismatchError, match="Image checksum mismatch"):
        load_digitised_deactivation(
            citation_label="test deact",
            source_key="barberis2007",
            image_path=img_path,
            image_sha256=tampered_sha256,
            figure_id=valid_figure_id,
            calibration=valid_linear_calibration,
            method="WebPlotDigitizer",
            operator="agent-brief-21",
            time_ms=np.linspace(0, 1000, 50),
            normalized_current=np.exp(-np.linspace(0, 1000, 50) / 50.0),
        )


def test_loader_refuses_empty_metadata(dummy_image_file, valid_figure_id, valid_linear_calibration):
    """Hard failure if method or operator is empty."""
    img_path, sha256 = dummy_image_file
    with pytest.raises(EvidenceError, match="non-empty method"):
        load_digitised_deactivation(
            citation_label="test deact",
            source_key="barberis2007",
            image_path=img_path,
            image_sha256=sha256,
            figure_id=valid_figure_id,
            calibration=valid_linear_calibration,
            method="   ",
            operator="operator",
            time_ms=[0.0, 10.0],
            normalized_current=[1.0, 0.5],
        )

    with pytest.raises(EvidenceError, match="non-empty operator"):
        load_digitised_deactivation(
            citation_label="test deact",
            source_key="barberis2007",
            image_path=img_path,
            image_sha256=sha256,
            figure_id=valid_figure_id,
            calibration=valid_linear_calibration,
            method="method",
            operator="",
            time_ms=[0.0, 10.0],
            normalized_current=[1.0, 0.5],
        )


def test_loader_refuses_invalid_figure_identifier():
    """FigureIdentifier requires non-empty journal, DOI, figure_number, and panel."""
    with pytest.raises(InvalidFigureIdentifierError, match="journal"):
        FigureIdentifier(journal="", doi="10.1234/test", figure_number=1, panel="A")
    with pytest.raises(InvalidFigureIdentifierError, match="DOI"):
        FigureIdentifier(journal="J Physiol", doi="", figure_number=1, panel="A")
    with pytest.raises(InvalidFigureIdentifierError, match="figure_number"):
        FigureIdentifier(journal="J Physiol", doi="10.1234/test", figure_number="", panel="A")
    with pytest.raises(InvalidFigureIdentifierError, match="panel"):
        FigureIdentifier(journal="J Physiol", doi="10.1234/test", figure_number=1, panel="")


def test_loader_refuses_degenerate_calibration():
    """AxisCalibration requires non-coincident reference points."""
    with pytest.raises(InvalidCalibrationError, match="x_ref1.pixel and x_ref2.pixel must be distinct"):
        AxisCalibration(
            x_ref1=AxisReference(pixel=100.0, data=0.0),
            x_ref2=AxisReference(pixel=100.0, data=10.0),
            y_ref1=AxisReference(pixel=0.0, data=0.0),
            y_ref2=AxisReference(pixel=100.0, data=1.0),
        )

    with pytest.raises(InvalidCalibrationError, match="x_ref data values must be positive for log10"):
        AxisCalibration(
            x_ref1=AxisReference(pixel=10.0, data=0.0),  # non-positive in log scale
            x_ref2=AxisReference(pixel=100.0, data=10.0),
            y_ref1=AxisReference(pixel=0.0, data=0.0),
            y_ref2=AxisReference(pixel=100.0, data=1.0),
            x_scale="log10",
        )


def test_axis_calibration_transformations(valid_linear_calibration, valid_log_calibration):
    """Test linear and log10 pixel-to-data conversions."""
    # Linear: x_ref1=(100, 0), x_ref2=(600, 1000) => 2 data units per pixel
    assert valid_linear_calibration.pixel_to_data_x(100.0) == pytest.approx(0.0)
    assert valid_linear_calibration.pixel_to_data_x(350.0) == pytest.approx(500.0)
    assert valid_linear_calibration.pixel_to_data_x(600.0) == pytest.approx(1000.0)

    # Log10: x_ref1=(50, 0.1), x_ref2=(650, 1000) => 4 decades across 600 pixels (150 px/decade)
    assert valid_log_calibration.pixel_to_data_x(50.0) == pytest.approx(0.1)
    assert valid_log_calibration.pixel_to_data_x(200.0) == pytest.approx(1.0)
    assert valid_log_calibration.pixel_to_data_x(350.0) == pytest.approx(10.0)
    assert valid_log_calibration.pixel_to_data_x(650.0) == pytest.approx(1000.0)


# =======================================================================================
# 2. Loader Acceptance of Correct Records
# =======================================================================================

def test_loader_accepts_valid_deactivation_record(dummy_image_file, valid_figure_id, valid_linear_calibration):
    """Loader successfully constructs a DeactivationDataset when all evidence is present."""
    img_path, sha256 = dummy_image_file
    t = np.array([0.0, 2.0, 10.0, 50.0, 200.0, 1000.0])
    y = np.array([1.0, 0.8, 0.5, 0.25, 0.1, 0.01])

    ds = load_digitised_deactivation(
        citation_label="Barberis 2007 test digitised trace",
        source_key="barberis2007",
        image_path=img_path,
        image_sha256=sha256,
        figure_id=valid_figure_id,
        calibration=valid_linear_calibration,
        method="WebPlotDigitizer manual extraction",
        operator="agent-brief-21",
        time_ms=t,
        normalized_current=y,
        pulse_ms=2.0,
        conc_mm=10.0,
    )

    assert ds.kind == "digitised"
    assert ds.synthetic is False
    assert ds.evidence is not None
    assert ds.evidence.image_sha256 == sha256
    assert ds.evidence.figure_id == valid_figure_id
    assert ds.evidence.calibration == valid_linear_calibration
    assert np.allclose(ds.time_ms, t)
    assert np.allclose(ds.normalized_current, y)

    # Test registration succeeds with valid evidence
    mock_registry: dict[str, object] = {}
    register_digitised_dataset("test_deact", ds, registry=mock_registry)
    assert "test_deact" in mock_registry
    assert mock_registry["test_deact"] is ds


def test_loader_accepts_pixel_coordinates(dummy_image_file, valid_figure_id, valid_linear_calibration):
    """Loader accurately maps pixel coordinates through calibration to data units."""
    img_path, sha256 = dummy_image_file
    px = [100.0, 350.0, 600.0]  # maps to [0.0, 500.0, 1000.0] ms
    py = [100.0, 300.0, 500.0]  # maps to [1.0, 0.5, 0.0] normalized current

    ds = load_digitised_deactivation(
        citation_label="Barberis pixel test",
        source_key="barberis2007",
        image_path=img_path,
        image_sha256=sha256,
        figure_id=valid_figure_id,
        calibration=valid_linear_calibration,
        method="automated pixel edge detection",
        operator="agent-brief-21",
        pixel_x=px,
        pixel_y=py,
    )

    assert np.allclose(ds.time_ms, [0.0, 500.0, 1000.0])
    assert np.allclose(ds.normalized_current, [1.0, 0.5, 0.0])


def test_loader_accepts_valid_dose_response_record(dummy_image_file, valid_figure_id, valid_log_calibration):
    """Loader successfully constructs a DoseResponseDataset when all evidence is present."""
    img_path, sha256 = dummy_image_file
    concs = np.array([0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 1000.0])
    resp = np.array([0.01, 0.05, 0.15, 0.45, 0.75, 0.92, 1.0])

    ds = load_digitised_dose_response(
        citation_label="Jahn 1997 test digitised CRC",
        source_key="jahn_a1b2g2_kinetics",
        image_path=img_path,
        image_sha256=sha256,
        figure_id=valid_figure_id,
        calibration=valid_log_calibration,
        method="WebPlotDigitizer",
        operator="agent-brief-21",
        concs_um=concs,
        mean_response=resp,
    )

    assert ds.kind == "digitised"
    assert ds.synthetic is False
    assert ds.evidence is not None
    assert np.allclose(ds.concs_um, concs)
    assert np.allclose(ds.mean_response, resp)


def test_registration_refuses_without_evidence():
    """Registration refused if kind='digitised' and synthetic=False but evidence is missing."""
    from circuitpharm.fitting.data import DoseResponseDataset

    fake_unverified = DoseResponseDataset(
        citation_label="unverified real claim",
        preparation="prep",
        concs_um=np.array([1.0, 10.0]),
        mean_response=np.array([0.1, 0.9]),
        sem=np.array([0.01, 0.05]),
        kind="digitised",
        synthetic=False,
        source_key="jahn_a1b2g2_kinetics",
        figure="Fig 1",
        digitisation="claimed by hand",
    )

    with pytest.raises(EvidenceError, match="requires attached DigitisationEvidence"):
        register_digitised_dataset("fake", fake_unverified)


# =======================================================================================
# 3. Self-Validation Check Tests
# =======================================================================================

def test_barberis_self_validation_passes_on_published_components(dummy_image_file, valid_figure_id, valid_linear_calibration):
    """Unperturbed Barberis deactivation curve reproduces published components within published errors.

    Published components read directly from circuitpharm.literature.DEACTIVATIONS[1]:
      tau = (2.8, 33.4, 221.35) ms, errors = (0.3, 4.6, 14.9) ms
      areas = (0.57, 0.23, 0.20), errors = (0.04, 0.02, 0.03)
    """
    img_path, sha256 = dummy_image_file
    barberis_row = DEACTIVATIONS[1]
    pub_taus, pub_areas = barberis_row.components

    # Construct trace from published components
    t = BARBERIS2007_DEACT.time_ms
    y = np.zeros_like(t)
    for tau, area in zip(pub_taus, pub_areas):
        y += area * np.exp(-t / tau)

    ds = load_digitised_deactivation(
        citation_label="Barberis exact unperturbed reconstruction (test fixture)",
        source_key="barberis2007",
        image_path=img_path,
        image_sha256=sha256,
        figure_id=valid_figure_id,
        calibration=valid_linear_calibration,
        method="test fixture",
        operator="test",
        time_ms=t,
        normalized_current=y,
    )

    res = validate_barberis_deactivation(ds, tolerance_sigmas=1.0)
    assert res.passed is True, f"Unperturbed Barberis trace failed self-validation: {res.message}"
    assert all(abs(z) <= 1.0 for z in res.z_scores.values())


def test_barberis_self_validation_catches_deliberately_perturbed_tau(dummy_image_file, valid_figure_id, valid_linear_calibration):
    """Self-validation detects a curve whose fast tau was deliberately perturbed beyond published error.

    Published tau_1 is 2.8 +/- 0.3 ms. Perturbing tau_1 to 4.5 ms (+5.67 sigma) must fail.
    Explicitly labeled as a TEST FIXTURE, NOT REAL DATA.
    """
    img_path, sha256 = dummy_image_file
    barberis_row = DEACTIVATIONS[1]
    pub_taus, pub_areas = barberis_row.components

    # Perturb tau_1 from 2.8 ms to 4.5 ms (well outside 2.8 +/- 0.3 ms)
    perturbed_taus = (4.5, pub_taus[1], pub_taus[2])
    t = BARBERIS2007_DEACT.time_ms
    y_perturbed = np.zeros_like(t)
    for tau, area in zip(perturbed_taus, pub_areas):
        y_perturbed += area * np.exp(-t / tau)

    perturbed_fixture = load_digitised_deactivation(
        citation_label="Deliberately perturbed tau_1 test fixture (NOT DATA)",
        source_key="barberis2007",
        image_path=img_path,
        image_sha256=sha256,
        figure_id=valid_figure_id,
        calibration=valid_linear_calibration,
        method="deliberate perturbation test fixture",
        operator="test",
        time_ms=t,
        normalized_current=y_perturbed,
    )

    res = validate_barberis_deactivation(perturbed_fixture, tolerance_sigmas=1.0)
    assert res.passed is False, "Deliberately perturbed tau_1 was NOT caught by self-validation!"
    assert "tau_1" in res.message
    assert abs(res.z_scores["tau_1"]) > 1.0

    with pytest.raises(SelfValidationError, match="Self-validation failed"):
        assert_self_validation(perturbed_fixture)


def test_barberis_self_validation_catches_deliberately_perturbed_area(dummy_image_file, valid_figure_id, valid_linear_calibration):
    """Self-validation detects a curve whose fast area was deliberately perturbed beyond published error.

    Published area_1 is 0.57 +/- 0.04. Perturbing area_1 to 0.35 (-5.5 sigma) must fail.
    Explicitly labeled as a TEST FIXTURE, NOT REAL DATA.
    """
    img_path, sha256 = dummy_image_file
    barberis_row = DEACTIVATIONS[1]
    pub_taus, pub_areas = barberis_row.components

    # Perturb area_1 from 0.57 to 0.35, compensating area_2 to maintain sum=1
    perturbed_areas = (0.35, 0.45, 0.20)
    t = BARBERIS2007_DEACT.time_ms
    y_perturbed = np.zeros_like(t)
    for tau, area in zip(pub_taus, perturbed_areas):
        y_perturbed += area * np.exp(-t / tau)

    perturbed_fixture = load_digitised_deactivation(
        citation_label="Deliberately perturbed area_1 test fixture (NOT DATA)",
        source_key="barberis2007",
        image_path=img_path,
        image_sha256=sha256,
        figure_id=valid_figure_id,
        calibration=valid_linear_calibration,
        method="deliberate perturbation test fixture",
        operator="test",
        time_ms=t,
        normalized_current=y_perturbed,
    )

    res = validate_barberis_deactivation(perturbed_fixture, tolerance_sigmas=1.0)
    assert res.passed is False, "Deliberately perturbed area_1 was NOT caught by self-validation!"
    assert "area_1" in res.message
    assert abs(res.z_scores["area_1"]) > 1.0


def test_jahn_self_validation_passes_on_published_parameters(dummy_image_file, valid_figure_id, valid_log_calibration):
    """Unperturbed Jahn CRC reproduces EC50 = 11.6 +/- 0.9 uM and nH = 2.2 +/- 0.4.

    Published parameters read directly from get_jahn_published_records().
    """
    img_path, sha256 = dummy_image_file
    jahn_rec = get_jahn_published_records()
    pub_ec50 = jahn_rec["ec50_um"]
    pub_hill = jahn_rec["hill"]

    c = JAHN1997_PEAK_CRC.concs_um
    x = (c / pub_ec50) ** pub_hill
    y = x / (1.0 + x)

    ds = load_digitised_dose_response(
        citation_label="Jahn unperturbed test fixture",
        source_key="jahn_a1b2g2_kinetics",
        image_path=img_path,
        image_sha256=sha256,
        figure_id=valid_figure_id,
        calibration=valid_log_calibration,
        method="test fixture",
        operator="test",
        concs_um=c,
        mean_response=y,
    )

    res = validate_jahn_crc(ds, tolerance_sigmas=1.0)
    assert res.passed is True, f"Unperturbed Jahn CRC failed self-validation: {res.message}"
    assert abs(res.z_scores["ec50_um"]) <= 1.0
    assert abs(res.z_scores["hill_slope"]) <= 1.0


def test_jahn_self_validation_catches_deliberately_perturbed_ec50(dummy_image_file, valid_figure_id, valid_log_calibration):
    """Self-validation detects a CRC whose EC50 was deliberately perturbed beyond published error.

    Published EC50 is 11.6 +/- 0.9 uM. Perturbing EC50 to 20.0 uM (+9.33 sigma) must fail.
    Explicitly labeled as a TEST FIXTURE, NOT REAL DATA.
    """
    img_path, sha256 = dummy_image_file
    jahn_rec = get_jahn_published_records()
    pub_hill = jahn_rec["hill"]

    # Perturb EC50 from 11.6 uM to 20.0 uM
    c = JAHN1997_PEAK_CRC.concs_um
    x_perturbed = (c / 20.0) ** pub_hill
    y_perturbed = x_perturbed / (1.0 + x_perturbed)

    perturbed_fixture = load_digitised_dose_response(
        citation_label="Deliberately perturbed EC50 test fixture (NOT DATA)",
        source_key="jahn_a1b2g2_kinetics",
        image_path=img_path,
        image_sha256=sha256,
        figure_id=valid_figure_id,
        calibration=valid_log_calibration,
        method="deliberate perturbation test fixture",
        operator="test",
        concs_um=c,
        mean_response=y_perturbed,
    )

    res = validate_jahn_crc(perturbed_fixture, tolerance_sigmas=1.0)
    assert res.passed is False, "Deliberately perturbed EC50 was NOT caught by self-validation!"
    assert "EC50" in res.message
    assert abs(res.z_scores["ec50_um"]) > 1.0

    with pytest.raises(SelfValidationError, match="Self-validation failed"):
        assert_self_validation(perturbed_fixture)


def test_jahn_self_validation_catches_deliberately_perturbed_hill_slope(dummy_image_file, valid_figure_id, valid_log_calibration):
    """Self-validation detects a CRC whose Hill slope was deliberately perturbed beyond published error.

    Published Hill slope is 2.2 +/- 0.4. Perturbing Hill slope to 1.2 (-2.5 sigma) must fail.
    Explicitly labeled as a TEST FIXTURE, NOT REAL DATA.
    """
    img_path, sha256 = dummy_image_file
    jahn_rec = get_jahn_published_records()
    pub_ec50 = jahn_rec["ec50_um"]

    # Perturb Hill slope from 2.2 to 1.2
    c = JAHN1997_PEAK_CRC.concs_um
    x_perturbed = (c / pub_ec50) ** 1.2
    y_perturbed = x_perturbed / (1.0 + x_perturbed)

    perturbed_fixture = load_digitised_dose_response(
        citation_label="Deliberately perturbed Hill slope test fixture (NOT DATA)",
        source_key="jahn_a1b2g2_kinetics",
        image_path=img_path,
        image_sha256=sha256,
        figure_id=valid_figure_id,
        calibration=valid_log_calibration,
        method="deliberate perturbation test fixture",
        operator="test",
        concs_um=c,
        mean_response=y_perturbed,
    )

    res = validate_jahn_crc(perturbed_fixture, tolerance_sigmas=1.0)
    assert res.passed is False, "Deliberately perturbed Hill slope was NOT caught by self-validation!"
    assert "Hill slope" in res.message
    assert abs(res.z_scores["hill_slope"]) > 1.0
