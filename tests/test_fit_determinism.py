"""The kinetic fit must be a property of the model, not of the machine.

THE DEFECT (roadmap P0-7 / C7). `fit_scheme` optimises four parameters against three
anchors, so there is a one-dimensional family of EXACT fits and `least_squares` lands
wherever its trust-region steps and the linear-algebra backend's rounding take it.
Measured on this repository, every row reproducing EC50 = 20.0000 uM, P_o,max = 0.7500 and
tau = 15.000 ms to the precision the manuscript quotes:

    environment                      k_on        K_d        asymptotic headroom
    manuscript's own rate table    0.0146842   31.95 uM         210.7x
    scipy 1.11.4 / numpy 1.26.4    0.0109865   29.44 uM         181.8x
    scipy 1.14.1 / numpy 1.26.4    0.0088044   26.60 uM         151.1x
    scipy 1.17.1 / numpy 2.4.6     0.0110210   29.48 uM         182.2x

A 40% spread in the headline number from the library version, and the manuscript's value
reproduced by none of them. That is why six tests were red on a fresh install, and why a
dependency pin would have been a promise that could not be kept -- the degeneracy is in the
problem, so the answer would still move with the BLAS vendor, the thread count and the CPU.

A MINIMUM-NORM REGULARISATION WAS TRIED FIRST AND REJECTED ON MEASUREMENT. It narrowed the
cross-version spread to ~0.8% but did not actually select the minimum-norm fit -- its own
acceptance test failed, because the optimiser converges on the anchors before a weak
penalty explores the flat valley, and raising the weight far enough to bite pulled the
anchors off (EC50 20.0111, tau 14.95 at weight 0.1). That failure is recorded here because
it is the more tempting fix and it does not work.

THE FIX is to hold `alpha` at `FIT_FIXED_ALPHA`, making the system SQUARE: three unknowns
(kon, koff, beta) against three anchors. Well-posed rather than regularised. 1/alpha is the
mean open time -- a real, standard, single-channel observable -- so this is a CONVENTION
about a measurable quantity, and roadmap P1-2 can replace it with a digitised value,
whereupon the fit becomes data-determined with no call site changed.

These tests are the permanent record of that. Amend them when a measured mean open time
lands; do not delete them.
"""
from __future__ import annotations

import numpy as np
import pytest

from circuitpharm.config import SYNAPTIC_PULSE
from circuitpharm.gabaa_kinetics import (
    DEFAULT_RATES,
    FIT_FIXED_ALPHA,
    FIT_RANGES,
    FIT_TARGETS,
    Scheme,
    fit_scheme,
)

PULSE = dict(SYNAPTIC_PULSE)

#: The fit this repository now produces on the config pulse with alpha held at
#: FIT_FIXED_ALPHA. Measured under scipy 1.17.1 / numpy 2.4.6; the square system is
#: well-posed, so this is a property of the model and the convention rather than of the
#: environment.
EXPECTED = dict(kon=0.0075785, koff=0.180415, beta=1.177468,
                kd_um=23.8062, gating_efficacy=3.9249, headroom=123.344)


def _headroom(s: Scheme, ambient: float = 0.40) -> float:
    po_inf = 1.0 / (1.0 + (s.alpha / s.beta) * (1.0 + s.d / s.r))
    return po_inf / s.po_tonic(ambient)


@pytest.fixture(scope="module")
def fitted() -> Scheme:
    return fit_scheme(verbose=False, pulse=PULSE)


def test_the_three_anchors_are_still_reproduced(fitted):
    """ANCHOR. The tie-breaker must not cost the fit its anchors.

    Tolerances are the precision the manuscript quotes each anchor to: EC50 to 4 decimals,
    P_o,max to 4, tau to 3. If a tie-breaker weight ever has to be raised far enough to
    break these, it is too strong and the degeneracy needs a real fourth observable instead.
    """
    m = fitted.ipsc_metrics(**PULSE)
    assert fitted.ec50_um() == pytest.approx(FIT_TARGETS["ec50_um"], abs=5e-4)
    assert fitted.po_max() == pytest.approx(FIT_TARGETS["po_max"], abs=5e-5)
    assert m["tau"] == pytest.approx(FIT_TARGETS["tau_ms"], abs=5e-3)
    # and each must sit inside its own declared accepted range, not merely near the target
    for got, key in ((fitted.ec50_um(), "ec50_um"), (fitted.po_max(), "po_max"),
                     (m["tau"], "tau_ms")):
        lo, hi = FIT_RANGES[key]
        assert lo <= got <= hi


def test_the_fit_is_the_value_this_repository_pins(fitted):
    """ANCHOR. The fit is now a stated number, not whatever the optimiser found today.

    rtol=1e-4 on each quantity. That is deliberately looser than the 6 figures pinned in
    EXPECTED and deliberately far tighter than the 40% spread the degenerate fit showed:
    it leaves room for a BLAS difference at the last digit or two while failing outright if
    the tie-breaker is removed or weakened. If this test fails, do NOT widen it -- check
    whether `FIT_TIEBREAK_WEIGHT` changed.
    """
    assert fitted.kon == pytest.approx(EXPECTED["kon"], rel=1e-4)
    assert fitted.koff == pytest.approx(EXPECTED["koff"], rel=1e-4)
    assert fitted.koff / fitted.kon == pytest.approx(EXPECTED["kd_um"], rel=1e-4)
    assert _headroom(fitted) == pytest.approx(EXPECTED["headroom"], rel=1e-4)


def test_the_fit_is_reproducible_within_this_process(fitted):
    """Two calls must agree exactly -- no hidden state, no RNG."""
    again = fit_scheme(verbose=False, pulse=PULSE)
    for attr in ("kon", "koff", "beta", "alpha"):
        assert getattr(again, attr) == getattr(fitted, attr)


def test_the_square_fit_is_well_posed_not_merely_regularised(fitted):
    """ANCHOR. The defining property: a distant start lands in the same place.

    rtol=1e-6 from a starting guess an order of magnitude away in every parameter. A
    degenerate problem does not do this -- it lands wherever the trust region happens to
    stop -- so this is the test that distinguishes "well-posed" from "regularised", and it
    is why holding alpha was chosen over a minimum-norm penalty.
    """
    far = fit_scheme(verbose=False, pulse=PULSE, _x0_log10=[np.log10(0.05),
                                                            np.log10(0.9),
                                                            np.log10(0.25)])
    for attr in ("kon", "koff", "beta"):
        assert getattr(far, attr) == pytest.approx(getattr(fitted, attr), rel=1e-6), (
            f"{attr} differs from a distant start: {getattr(far, attr)!r} vs "
            f"{getattr(fitted, attr)!r} -- the fit is not well-posed")


def test_alpha_is_held_and_is_a_plausible_mean_open_time(fitted):
    """The convention must be a physically meaningful number, not an arbitrary constant."""
    assert fitted.alpha == FIT_FIXED_ALPHA
    mean_open_ms = 1.0 / fitted.alpha
    # alpha1beta2gamma2 single-channel mean open times are of order 1-5 ms
    assert 1.0 <= mean_open_ms <= 5.0, (
        f"mean open time {mean_open_ms:.2f} ms is outside the 1-5 ms range that makes this "
        f"convention defensible as a stand-in for a measurement")


def test_the_degenerate_fit_still_reproduces_the_anchors():
    """PERMANENT RECORD. The degeneracy is real: the old fit was not WRONG, it was unpinned.

    Both fits reproduce all three anchors. That is the whole point -- the anchors cannot
    distinguish them, so the manuscript's 7 significant figures were recording the
    optimiser's path. Keep this test: it is the evidence that the tie-breaker is a
    convention rather than a correction.
    """
    degenerate = fit_scheme(verbose=False, pulse=PULSE, fixed_alpha=None)
    m = degenerate.ipsc_metrics(**PULSE)
    assert degenerate.ec50_um() == pytest.approx(FIT_TARGETS["ec50_um"], abs=5e-4)
    assert degenerate.po_max() == pytest.approx(FIT_TARGETS["po_max"], abs=5e-4)
    assert m["tau"] == pytest.approx(FIT_TARGETS["tau_ms"], abs=5e-2)
    # ... and yet it is a materially different parameter set
    assert abs(degenerate.koff / degenerate.kon - EXPECTED["kd_um"]) > 0.5


def test_the_manuscripts_own_rate_table_also_fits(fitted):
    """PERMANENT RECORD. The manuscript's numbers are a third exact fit.

    k_on 0.0146842, k_off 0.469177, beta 0.559023, alpha 0.107891 -- the table in
    knowledge/08-manuscript.md §2.1. It reproduces every anchor and gives K_d 31.95 uM and
    210.7x headroom, against this repository's 28.48 uM and 171.2x. Neither is more correct
    than the other given these three anchors; that is exactly the problem.
    """
    ms = Scheme(kon=0.0146842, koff=0.469177, beta=0.559023, alpha=0.107891)
    m = ms.ipsc_metrics(**PULSE)
    assert ms.ec50_um() == pytest.approx(20.0, abs=1e-3)
    assert ms.po_max() == pytest.approx(0.75, abs=1e-3)
    assert m["tau"] == pytest.approx(15.0, abs=1e-2)
    assert ms.koff / ms.kon == pytest.approx(31.95, abs=0.01)
    assert _headroom(ms) == pytest.approx(210.7, rel=1e-3)
    # and it differs from ours by far more than any tolerance above
    assert abs(_headroom(ms) - _headroom(fitted)) > 50.0


def test_the_anchors_are_solved_not_traded_off(fitted):
    """ANCHOR. A square well-posed system should drive the residual to ~0, not balance it.

    |residual| < 1e-6 over the three anchors in their own normalised units. The measured
    value is ~7e-10. A regularised fit cannot reach this, because the penalty competes with
    the anchors -- which is the quantitative form of why the minimum-norm attempt was
    rejected.
    """
    anchors = np.array([
        np.log(fitted.ec50_um() / FIT_TARGETS["ec50_um"]),
        (fitted.po_max() - FIT_TARGETS["po_max"]) / 0.05,
        np.log(fitted.ipsc_metrics(**PULSE)["tau"] / FIT_TARGETS["tau_ms"]),
    ])
    assert float(np.linalg.norm(anchors)) < 1e-6, (
        f"anchor residual {np.linalg.norm(anchors):.3e}; the system is not being solved")
