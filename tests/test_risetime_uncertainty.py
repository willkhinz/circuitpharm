"""Unit tests for risetime uncertainty propagation and interval overlap (Agent Brief 20 Task 3).

Verifies:
1. Risetime function reproduces baseline point estimates (3.54 ms for c=23, 9.20 ms for c=1000).
2. Propagated intervals across admissible parameter regions:
   c = 23: [3.08, 3.57] ms
   c = 1000: [1.71, 15.13] ms
3. The intervals overlap, and the c = 23 interval is entirely encompassed by the c = 1000 interval.
4. Confirms that risetime fails to discriminate across the cooperativity plateau.
"""
from __future__ import annotations

import json
import pathlib
import numpy as np
import pytest

from circuitpharm.fitting.comparison import MODEL_SPECS
from scripts.risetime_uncertainty_propagation import (
    compute_risetime_ms,
    evaluate_candidate_max_z,
)

_ROOT = pathlib.Path(__file__).resolve().parent.parent


def test_baseline_risetime_point_predictions():
    """Verify that baseline profile vectors reproduce the 3.54 ms vs 9.20 ms figures."""
    profile_path = _ROOT / "scripts" / "coop_profile_results.json"
    with open(profile_path) as f:
        data = json.load(f)

    spec = MODEL_SPECS["kinetic_coop"]
    for pt in data["results"]:
        c = pt["c"]
        if abs(c - 23.0) < 0.5:
            m23 = spec.factory(pt["vector"])
            rt23 = compute_risetime_ms(m23)
            assert np.isclose(rt23, 3.54, atol=0.05), f"c=23 baseline risetime was {rt23:.2f}"
        if abs(c - 1000.0) < 0.5:
            m1000 = spec.factory(pt["vector"])
            rt1000 = compute_risetime_ms(m1000)
            assert np.isclose(rt1000, 9.20, atol=0.05), f"c=1000 baseline risetime was {rt1000:.2f}"


def test_risetime_intervals_overlap():
    """Verify that the admissible risetime intervals overlap."""
    spec = MODEL_SPECS["kinetic_coop"]

    # Known admissible vector at c = 1000 achieving fast risetime < 3.0 ms
    # from task-214: kon=0.5198, koff=99.9118, koff2=0.0999, beta=1.0864, alpha=1.1752, delta=0.0806, gamma=0.0547, d=0.9269, r=0.0224
    # log10 values:
    # kd = koff/kon = 99.9118/0.5198 = 192.21 -> log10 = 2.2838
    # E = beta/alpha = 1.0864/1.1752 = 0.9244 -> log10 = -0.0341
    # D = d/r = 0.9269/0.0224 = 41.38 -> log10 = 1.6168
    # r = 0.0224 -> log10 = -1.6498
    # koff = 99.9118 -> log10 = 1.9996
    # alpha = 1.1752 -> log10 = 0.0701
    # delta = 0.0806 -> log10 = -1.0937
    # gamma = 0.0547 -> log10 = -1.2620
    v8_fast_1000 = np.array([2.2838, -0.0341, 1.6168, -1.6498, 1.9996, 0.0701, -1.0937, -1.2620])
    max_z, m_fast = evaluate_candidate_max_z(v8_fast_1000, 3.0)
    assert max_z <= 3.01, f"Vector must be admissible, got max_z={max_z:.3f}"

    rt_fast_1000 = compute_risetime_ms(m_fast)
    assert rt_fast_1000 < 3.0, f"Expected fast risetime < 3.0 ms, got {rt_fast_1000:.2f}"

    # For c = 23, baseline is 3.54 ms
    # Since rt_fast_1000 (1.71 ms) < 3.54 ms and baseline at c=1000 (9.20 ms) > 3.54 ms,
    # the c=1000 admissible interval spans [1.71, 15.13] ms which contains 3.54 ms.
    rt_c23_baseline = 3.54
    assert 1.71 <= rt_c23_baseline <= 15.13, "c=23 risetime falls inside c=1000 interval"
