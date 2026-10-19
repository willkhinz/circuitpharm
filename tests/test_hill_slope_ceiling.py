"""Executable tests for Unconstrained Hill Slope Ceilings (Brief 17 Task 3).

Distinguishes:
- Statement A: 'The scheme cannot produce n_H = 2.2 at all' (FALSE)
- Statement B: 'The scheme cannot produce n_H = 2.2 while satisfying the other observables' (TRUE)
"""
from __future__ import annotations

import pathlib
import sys

_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

import numpy as np
import pytest

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.crc_components import extract_crc_parameters


def test_jw95_unconstrained_hill_slope_exceeds_target():
    """Verify that even 5-state JW95 with 3 free parameters reaches n_H > 2.2 unconstrained."""
    spec = MODEL_SPECS["kinetic_jw95"]
    # Optimum vector found by unconstrained differential evolution
    v = np.array([1.2881204286837873, 3.0, -2.6805855906866864], dtype=float)

    m = spec.factory(v)
    crc = extract_crc_parameters(m)

    # Measured Hill slope over 1-10 uM window (matching Jahn 1997)
    assert crc.hill_measured > 2.2, f"Expected n_H > 2.2, got {crc.hill_measured:.3f}"
    assert crc.hill_measured > 2.8, f"Expected n_H > 2.8, got {crc.hill_measured:.3f}"


def test_statement_a_is_false_statement_b_is_true():
    """Verify distinction between Statement A (false) and Statement B (true)."""
    spec = MODEL_SPECS["kinetic_jw95"]
    v = np.array([1.2881204286837873, 3.0, -2.6805855906866864], dtype=float)

    m = spec.factory(v)
    crc = extract_crc_parameters(m)

    # Statement A is FALSE: The model can easily produce n_H = 2.2 (in fact reaching ~2.93)
    statement_a_cannot_produce_22 = bool(crc.hill_measured < 2.2)
    assert not statement_a_cannot_produce_22, "Statement A should be FALSE"

    # Statement B is TRUE: At this high-Hill regime, other observables are severely violated:
    # 1. EC50 is 1.82 uM vs published 11.6 +/- 0.9 uM (> 10 sigma deviation)
    # 2. Po_max is ~0.999 vs accepted range [0.70, 0.80]
    ec50_distorted = abs(crc.ec50_um - 11.6) > 5.0
    po_max_distorted = float(m.po_max) > 0.90
    statement_b_joint_conflict = ec50_distorted and po_max_distorted
    assert statement_b_joint_conflict, "Statement B should be TRUE (other observables violated)"
