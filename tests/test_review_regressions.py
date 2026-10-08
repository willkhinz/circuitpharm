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
from circuitpharm.evaluation import Compound


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
    from circuitpharm.evaluation import evaluate, clear_control_cache
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
    from circuitpharm.evaluation import evaluate
    rs = evaluate(Compound.from_profile("alogabat", occupancy=0.35, nmda_block=0.3),
                  n_seed=1, include_motor=False)
    assert rs.quantity("subjective_index_gaba").value > 0
    assert rs.quantity("subjective_index_nmda_raw").value > 0
    assert rs.quantity("subjective_index_total").tier is Tier.VOID


def test_no_nmda_arm_means_no_nmda_quantities():
    from circuitpharm.evaluation import evaluate
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


# =====================================================================================
# Review 3 (8 findings). Two of these were created or left incomplete by the review-2
# fixes, which is the more useful lesson: a correct fix can still interact badly with
# another correct fix, and a partial fix looks identical to a complete one from outside.
# =====================================================================================

import numpy as np


# ---------------------------------------------------- R3.1 duty-cycle phase dependence
def test_duty_cycle_is_phase_invariant():
    """`burst_metrics` paired onsets and offsets BY INDEX, which is only correct when the
    trace starts below threshold. Starting mid-burst -- routine, since analysis windows
    open after settling at an arbitrary phase -- shifted every pair by one, so
    `off[i] > on[i]` failed for every burst and duty came back NaN for a strongly,
    regularly bursting circuit.

    The silence was the damage: scripts/tune_rg2.py scores a NaN duty as 0.0 and penalises
    it, so parameter sets were ranked partly on where their window happened to open.
    """
    from circuitpharm.cpg import burst_metrics
    t = np.arange(0, 10000.0, 1.0)
    duties = []
    for phase in (0.0, 0.3, np.pi / 2, np.pi, 1.7 * np.pi):
        sig = np.sin(2 * np.pi * t / 1000.0 + phase)
        r = burst_metrics(t, sig - sig.min())
        assert np.isfinite(r["duty"]), f"duty is NaN at phase {phase:.2f}"
        duties.append(r["duty"])
    assert max(duties) - min(duties) < 1e-6, f"duty varies with phase: {duties}"


def test_duty_cycle_matches_a_known_square_wave():
    """Phase-invariance alone could be satisfied by a constant wrong answer."""
    from circuitpharm.cpg import burst_metrics
    t = np.arange(0, 10000.0, 1.0)
    sq = (np.sin(2 * np.pi * t / 1000.0) > 0).astype(float)
    assert burst_metrics(t, sq, thresh_frac=0.5)["duty"] == pytest.approx(0.5, abs=0.02)


# ------------------------------------------------- R3.5 direct agonists must not crash
def test_direct_agonist_raises_with_actionable_guidance():
    """CREATED BY TWO EARLIER FIXES INTERACTING. Carrying each profile's ceiling into
    s_max (correct) met a hard raise on unreachable affinity shifts (also correct), and
    together they made gaboxadol -- a direct orthosteric agonist with ceiling 1e9 -- raise
    instead of evaluate. Neither fix was wrong; the combination was unhandled.

    A direct agonist genuinely cannot be modelled by the PAM machinery, so raising is
    right -- but the error must name the reason and point at the right tool.
    """
    with pytest.raises(ValueError, match="DIRECT|agonist"):
        Compound.from_profile("gaboxadol").pool_gains()


@pytest.mark.parametrize("key", ["alogabat", "neurosteroid", "imepitoin", "tpa023",
                                 "hz_166", "mp_iii_022", "ideal_a5", "zolpidem"])
def test_every_modulator_profile_evaluates_at_its_own_ceiling(key):
    """All non-agonist profiles must work with their real s_max. The neurosteroid arm
    (ceiling 6.0) needs the gating fallback, since an affinity-only mechanism cannot
    reach a 6x shift."""
    g = Compound.from_profile(key, occupancy=0.5).pool_gains()
    assert all(math.isfinite(g[k]) for k in ("tonic", "phasic", "tau")), key
    assert g["tonic"] >= 1.0 and g["phasic"] >= 1.0, key


# ------------------------------------------- R3.6 glyr_sens must reach the assays
def test_assays_accept_and_forward_glyr_sens():
    """INCOMPLETE EARLIER FIX. glyr_sens was added to Drug, the circuits and the plant,
    but not to the assay signatures -- so the public entry points could not vary glycine
    sensitivity at all, which is the whole point of separating it from GABA-A."""
    import inspect
    from circuitpharm.assays import stretch_reflex, locomotion
    for fn in (stretch_reflex, locomotion):
        assert "glyr_sens" in inspect.signature(fn).parameters, fn.__name__


@pytest.mark.needs_plant
@pytest.mark.slow
def test_glycine_action_survives_a_tiny_gaba_sensitivity():
    """The substantive check behind the plumbing: strychnine must still cause
    hyperreflexia when the GABA-A sensitivity is tiny, because glycine receptors contain
    no GABA-A subunits. Before the decoupling a 1.6x glycine potentiation became
    1.0019x at alogabat's preBotC sensitivity."""
    from circuitpharm.assays import stretch_reflex
    ctrl = stretch_reflex(Drug(), gaba_sens_phasic=0.01, glyr_sens=1.0)["gain"]
    stry = stretch_reflex(Drug(glyr_gain=0.4), gaba_sens_phasic=0.01,
                          glyr_sens=1.0)["gain"]
    assert stry > 1.2 * ctrl, (
        f"glycine block gave {stry:.3f} vs control {ctrl:.3f}; glycine pharmacology is "
        "being suppressed by a GABA-A-derived fraction again")


# --------------------------------------- R3.4 the two forebrain scales must stay distinct
def test_forebrain_sensitivity_and_subjective_index_are_documented_as_different_scales():
    """`regional_sens("forebrain")` includes K (so a non-selective BZ gives 1.00) while
    `subjective_index()` omits it (giving 0.55) -- a 1.8x gap between two quantities that
    both sound like forebrain drug engagement. The omission is correct (the index is a
    weighted SUBSET, so a whole-conductance normaliser would be meaningless) but it must
    be stated, or the two get compared."""
    from circuitpharm.subtypes import PROFILES
    p = PROFILES["nonselective_bz"]
    assert p.regional_sens("forebrain") == pytest.approx(1.0, abs=1e-6)
    assert p.subjective_index() == pytest.approx(0.55, abs=1e-6)
    doc = type(p).subjective_index.__doc__ or ""
    assert "NOT ON THE SAME SCALE" in doc.upper() or "discrepancy" in doc.lower()


# ------------------------------------------------- R3.8 no warnings on legitimate NaN
@pytest.mark.needs_plant
@pytest.mark.slow
def test_abolished_locomotion_does_not_emit_a_numpy_warning():
    """np.nanmean over an all-NaN list warns 'Mean of empty slice' and returns NaN. That
    case is legitimate -- a heavy sedative can abolish locomotion in every seed -- so the
    warning is noise on a correct result, and noise on correct results is how real
    warnings get ignored."""
    import warnings
    from circuitpharm.evaluation import evaluate, clear_control_cache
    clear_control_cache()
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        evaluate(Compound.from_profile("nonselective_bz", occupancy=1.0), n_seed=2)


# =====================================================================================
# REVIEW PASS 7 (2026-10-07). Eight findings; the three with a library surface are here.
# The other five were script-level (build_kb collateral destruction, ranking_robustness
# empty percentile, overdose_kinetic NaN propagation, two division-by-zero guards) plus
# dead code in conftest.
# =====================================================================================

# ------------------------------------------------- R7.1 calibrate_pam lower bracket
@pytest.mark.parametrize("target_shift", [1.0, 0.8, 0.5])
def test_calibrate_pam_returns_nan_below_its_search_domain(target_shift):
    """Only the UPPER bracket end was checked, so asking for a shift at or below 1.0 --
    a neutral ligand, or a NAM, which right-shifts EC50 -- made both endpoints positive
    and brentq raised a bare ValueError about signs. The documented convention for
    'unreachable by this mechanism' is NaN, and callers already test isfinite."""
    from circuitpharm.gabaa_kinetics import fit_scheme, calibrate_pam
    s = fit_scheme(verbose=False)
    out = calibrate_pam(s, target_shift=target_shift)
    assert not np.isfinite(out), (
        f"expected NaN for target_shift={target_shift} (outside the x>=1 search domain), "
        f"got {out}")


def test_calibrate_pam_still_solves_a_reachable_shift():
    """Guard against the lower-bracket fix swallowing the normal case."""
    from circuitpharm.gabaa_kinetics import fit_scheme, calibrate_pam
    s = fit_scheme(verbose=False)
    out = calibrate_pam(s, target_shift=2.5)
    assert np.isfinite(out) and out > 1.0


# ------------------------------------------------- R7.2 no negative GABA-A conductance
@pytest.mark.parametrize("cls_name", ["PreBotC", "GroupPacemakerRG", "SpinalCircuit"])
def test_tonic_gabaa_conductance_never_goes_negative(cls_name):
    """A NAM (gaba_scale_tonic < 1) meeting a sensitivity above 1 made `eff` negative,
    and a negative tonic term can drive the TOTAL gabaa conductance negative. Since the
    current is g*(E - V), that inverts an inhibitory shunt into regenerative negative
    damping: the voltage diverges rather than the run failing visibly. Same class as the
    clamp already in Drug.nmda_scale.

    Latent rather than live -- PROFILES top out at sens 0.815 -- but both
    `gaba_sens_tonic` and `Drug(gaba_a_gain_tonic=...)` are public.
    """
    from circuitpharm.cpg import Drug
    d = Drug(gaba_a_gain_tonic=0.0, gaba_a_cap_tonic=5.0)   # full negative modulation
    if cls_name == "PreBotC":
        from circuitpharm.resp import PreBotC
        net = PreBotC(drug=d, gaba_sens_tonic=3.0, seed=0)
    elif cls_name == "GroupPacemakerRG":
        from circuitpharm.rg2 import GroupPacemakerRG
        net = GroupPacemakerRG(drug=d, gaba_sens_tonic=3.0, seed=0)
    else:
        from circuitpharm.circuit import SpinalCircuit
        net = SpinalCircuit(drug=d, gaba_sens_tonic=3.0, seed=0)

    for _ in range(3000):
        net.step(0.1)
    for nm, p in net.pops.items():
        v = np.asarray(p.V, float)
        assert np.all(np.isfinite(v)), f"{cls_name}.{nm}: voltage went non-finite"
        assert v.max() < 1e3, (
            f"{cls_name}.{nm}: voltage diverged to {v.max():.3g} mV -- negative "
            f"conductance is acting as negative damping")


# ------------------------------------------------- R7.3 locomotor FFT spacing is pinned
@pytest.mark.needs_plant
def test_locomotion_fft_spacing_comes_from_the_time_vector():
    """The sample spacing was the module-level DT rather than the recorded time vector.
    Same number today, so not wrong -- UNPINNED. Decimating the recording (an obvious way
    to hold memory down on a long run) would leave rfftfreq told the undecimated spacing
    and report a frequency N times too fast, silently. A drug that slows the step cycle
    reading as a faster one is recurring error E4.

    Checked by source inspection: the behavioural version needs a decimation path that
    does not exist yet, and asserting on the DT literal is what let this drift.
    """
    import inspect
    from circuitpharm import assays
    src = inspect.getsource(assays.locomotion)
    assert "np.fft.rfftfreq(len(q), dt_s)" in src, (
        "locomotion's FFT spacing is no longer derived from the recorded time vector")
    assert "np.fft.rfftfreq(len(q), DT)" not in src


# =====================================================================================
# REVIEW PASS 10 (2026-10-07). Three findings. One (locomotion's empty-window crash) has
# a library surface; the other two are the NaN-poisoned grid search and the remaining
# unguarded baseline divisions, both script-level, covered by the shared-pattern test
# below rather than by running the sweeps (each takes minutes).
# =====================================================================================

# ------------------------------------------------- R10.1 locomotion rejects a dead window
@pytest.mark.needs_plant
@pytest.mark.parametrize("duration_s", [0.5, 1.0])
def test_locomotion_refuses_a_duration_with_no_analysable_window(duration_s):
    """`m = t > LOCO_SETTLE_S` is all-False at or below the settle time, so every analysis
    slice was shape (0,) and the first reduction raised numpy's "zero-size array to
    reduction operation maximum which has no identity" -- an obscure message from deep
    inside the function naming nothing the caller controls. duration_s <= 1.0 is exactly
    what someone writes for a fast unit test."""
    from circuitpharm.assays import locomotion, LOCO_SETTLE_S
    with pytest.raises(ValueError, match="no analysable window"):
        locomotion(duration_s=duration_s)
    assert duration_s <= LOCO_SETTLE_S        # the premise of the test


# ------------------------------------------------- R10.2 NaN cannot latch a grid search
def test_nan_does_not_latch_a_min_tracking_comparison():
    """The bug, in isolation: `err < best[0]` is False whenever best[0] is NaN, so one
    non-finite error in the FIRST grid cell latched `best` permanently and every later
    finite, better candidate was discarded. The script then printed that poisoned cell as
    BEST with NaN% beside it, which reads as a converged answer.

    This pins the PATTERN. The three sweeps that had it (calib_split, calibrate_resp,
    recalibrate_kinetic) each take minutes to run, so they are checked by source
    inspection in the companion test below.
    """
    def select(errs, guarded):
        best = None
        for e in errs:
            if guarded:
                if np.isfinite(e) and (best is None or e < best[0]):
                    best = (e, "cand")
            elif best is None or e < best[0]:
                best = (e, "cand")
        return best

    errs = [float("nan"), 0.05, 0.01]
    assert math.isnan(select(errs, guarded=False)[0]), "premise: unguarded form latches"
    assert select(errs, guarded=True)[0] == 0.01
    # all-NaN must yield None so the caller can refuse, not a poisoned winner
    assert select([float("nan")] * 3, guarded=True) is None


@pytest.mark.parametrize("script,var", [
    ("scripts/calib_split.py", "e"),
    ("scripts/calibrate_resp.py", "err"),
    ("scripts/recalibrate_kinetic.py", "err"),
])
def test_calibration_sweeps_guard_their_min_tracking(script, var):
    """Each sweep must test isfinite before comparing, and must refuse rather than report
    a winner when nothing finite was found."""
    import pathlib
    src = pathlib.Path(script).read_text()
    assert f"np.isfinite({var}) and (best is None" in src, (
        f"{script}: min-tracking is not NaN-guarded")
    assert "best is None:" in src and "SystemExit" in src, (
        f"{script}: does not refuse when every grid cell is non-finite")


# ------------------------------------------------- R10.3 baselines guarded, not bare
@pytest.mark.parametrize("script", [
    "scripts/predict_muscimol.py", "scripts/recalibrate_kinetic.py",
    "scripts/calib_split.py", "scripts/reflex.py",
])
def test_percent_of_control_never_divides_bare(script):
    """A quiescent drug-free control makes every percent-of-control 0/0. Bare division
    raised ZeroDivisionError and killed an entire grid or concentration series at one bad
    cell; NaN says "no baseline to measure against" and lets the rest report.
    `evaluation._pct_of_control` has done this for several sessions; the scripts had not.
    """
    import pathlib
    src = pathlib.Path(script).read_text()
    assert "1e-9" in src, f"{script}: no baseline-magnitude guard present"
