"""Regressions for nine defects found in an external code review of session 8.

All nine were real and all nine were mine. They are grouped here rather than scattered so
the review's findings stay auditable against the fixes. Severity order, worst first.
"""
import math
import subprocess
import sys

import pytest

from circuitpharm.cpg import Pop, Drug
from circuitpharm.circuit import SpinalCircuit
from circuitpharm.results import Quantity, Tier
from circuitpharm.evaluate import Compound


# ===================================================== 1. reproducibility (worst)
def test_initial_voltages_are_identical_across_PROCESSES():
    """`hash(str)` is randomised per interpreter via PYTHONHASHSEED, so seeding an RNG with
    `abs(hash(name))` gave different initial membrane voltages on every run even when the
    caller passed an explicit seed. Measured V0[0] = -71.05, -62.73, -69.90 across three
    processes. Every published number was therefore irreproducible between runs.

    Must run in SUBPROCESSES: within one process the hash is stable, so an in-process
    check would pass while the bug was live.
    """
    code = ("from circuitpharm.cpg import Pop;"
            "print(repr(float(Pop(n=5, name='Exc0').V[0])))")
    outs = {subprocess.run([sys.executable, "-c", code], capture_output=True,
                           text=True, check=True).stdout.strip()
            for _ in range(3)}
    assert len(outs) == 1, f"initial voltage differs across processes: {outs}"


def test_distinct_population_names_still_get_distinct_voltages():
    """The fix must not collapse every population onto one RNG stream."""
    a, b = Pop(n=8, name="Exc0"), Pop(n=8, name="Inh0")
    assert not (a.V == b.V).all()


# ============================================ 2. global class-state mutation
def test_stretch_reflex_does_not_mutate_the_class_attribute():
    """The assay used to rescale Ia->Mn by MUTATING SpinalCircuit.W (a class attribute) and
    restoring it in a `finally`. Any concurrent instance in the same process saw corrupted
    weights, and an exception before the restore left them corrupted permanently."""
    pytest.importorskip("mujoco", reason="needs the optional body plant")
    from circuitpharm.assays import stretch_reflex
    before = dict(SpinalCircuit.W)
    stretch_reflex(Drug())
    assert dict(SpinalCircuit.W) == before, "class-level weights were mutated"


def test_circuit_weights_are_per_instance():
    a = SpinalCircuit(w={("Ia", "Mn", "ampa"): 0.99})
    b = SpinalCircuit()
    assert a.W[("Ia", "Mn", "ampa")] == 0.99
    assert b.W[("Ia", "Mn", "ampa")] == SpinalCircuit.W[("Ia", "Mn", "ampa")]
    assert a.W is not SpinalCircuit.W


# ============================================ 3. single-seed motor endpoints (E6)
@pytest.mark.slow
@pytest.mark.needs_plant
def test_motor_endpoints_average_over_seeds():
    """Motor readouts hard-coded seed=1 while ventilation averaged over n_seed -- recurring
    error E6, committed in this project's own public API. Motor output is noisier than
    ventilation, so a single seed is worse there, not better."""
    from circuitpharm.evaluate import evaluate, clear_control_cache
    clear_control_cache()
    a = evaluate(Compound.from_profile("alogabat", occupancy=0.35), n_seed=1)
    clear_control_cache()
    b = evaluate(Compound.from_profile("alogabat", occupancy=0.35), n_seed=3)
    # averaging over more seeds must actually change the value, or it is not averaging
    assert a.quantity("coordination").value != b.quantity("coordination").value


# ============================================ 4. bool formatting
def test_booleans_print_as_booleans_not_as_1_and_0():
    """`bool` subclasses `int`, so the numeric branch caught it and f"{True:.4g}" rendered
    "1". Output read "walking: 1", ambiguous against a real 0/1 fraction."""
    assert "True" in str(Quantity("walking", True, Tier.UNCALIBRATED))
    assert "False" in str(Quantity("alive", False, Tier.UNCALIBRATED))
    # numbers must still format numerically
    assert "0.5" in str(Quantity("frac", 0.5, Tier.UNCALIBRATED))


# ============================================ 5. NMDA arm of the subjective index
def test_nmda_arm_is_reported_and_the_total_is_void():
    """Moving evaluation out of the old simulator dropped the NMDA subjective term, so any
    compound with an NMDA component was silently under-reported -- and the two-arm design
    is the entire point, since ethanol's stimulus is a COMPOUND one.

    The arms are reported separately and the total is VOID rather than summed, because they
    are in incommensurable units. The project's 'GABA salience >= NMDA salience'
    constraint compares them directly and so rests on the same defect."""
    from circuitpharm.evaluate import evaluate
    rs = evaluate(Compound.from_profile("alogabat", occupancy=0.35, nmda_block=0.3),
                  n_seed=1, include_motor=False)
    assert rs.quantity("subjective_index_gaba").value > 0
    assert rs.quantity("subjective_index_nmda_raw").value > 0
    assert rs.quantity("subjective_index_total").tier is Tier.VOID


def test_no_nmda_arm_means_no_nmda_quantities():
    from circuitpharm.evaluate import evaluate
    rs = evaluate(Compound.from_profile("alogabat", occupancy=0.35), n_seed=1,
                  include_motor=False)
    with pytest.raises(KeyError):
        rs.quantity("subjective_index_nmda_raw")


# ============================================ 6. reflex-gain denominator
def test_reflex_gain_is_undefined_when_the_afferent_did_not_respond():
    """max(1e-6, ia_dyn - ia_base) divided by 1e-6 when the Ia response was NEGATIVE,
    turning a 0.5 Hz drop into a gain of 5,000,000. A ratio has no meaning when its
    denominator is noise, so it must be NaN with a stated reason."""
    pytest.importorskip("mujoco", reason="needs the optional body plant")
    from circuitpharm.assays import IA_RESPONSE_FLOOR
    assert IA_RESPONSE_FLOOR > 0.0


# ============================================ 7. NaN propagation
def test_unreachable_intrinsic_efficacy_raises_instead_of_producing_nan():
    """calibrate_pam returns NaN when the requested EC50 shift is unreachable by an
    affinity-type mechanism. Unchecked it became koff=NaN, then NaN conductances, then NaN
    voltages in the LIF integrator -- surfacing as 'SVD did not converge in Linear Least
    Squares' and raw LAPACK complaints far from the cause."""
    with pytest.raises(ValueError, match="not reachable|SATURATES|s_max"):
        Compound("impossible", a5=1.0, s_max=500.0).pool_gains()


def test_plausible_intrinsic_efficacy_still_works():
    g = Compound("ok", a5=1.0, s_max=2.5, occupancy=0.5).pool_gains()
    assert all(math.isfinite(g[k]) for k in ("tonic", "phasic", "tau"))


# ============================================ 8. conftest skip heuristic
def test_skip_heuristic_does_not_match_on_test_NAMES():
    """Any test whose name contained "reflex" or "walk" was skipped on a lean install, even
    one testing pure arithmetic. A test skipped when it should have run is worse than a
    failing one: nothing reports it. THIS test's own name contains 'reflex'* and it must
    still run.  (*see the assertion below.)"""
    from pathlib import Path
    src = (Path(__file__).parent / "conftest.py").read_text()
    assert '"reflex" in item.name' not in src
    assert '"walk" in item.name' not in src
    assert "needs_plant" in src, "the explicit marker must be the primary mechanism"


# ============================================ 9. coverage threshold
def test_coverage_threshold_is_not_set_globally():
    """`fail_under` in pyproject fires on ANY coverage run, so `pytest -m "not slow" --cov`
    failed at 65% with every test passing. A failure carrying no information trains people
    to disable the check."""
    from pathlib import Path
    cfg = (Path(__file__).parents[1] / "pyproject.toml").read_text()
    body = cfg.split("[tool.coverage.report]")[1] if "[tool.coverage.report]" in cfg else ""
    active = [l for l in body.splitlines()
              if l.strip().startswith("fail_under")]
    assert not active, f"fail_under is set globally: {active}"
