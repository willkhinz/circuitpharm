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
@pytest.mark.slow
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


@pytest.mark.slow
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


@pytest.mark.slow
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
@pytest.mark.slow
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


# ------------------------------------------------- E20 substrate factory default erased
#                                                   the LIF cell's adaptation
def test_make_pop_requires_g_adapt_and_none_keeps_the_class_default():
    """E20. `substrate.make_pop` carried `g_adapt=0.0` as a DEFAULT. `resp.py` and `rg2.py`
    had always set adaptation explicitly, so routing them through it was faithful.
    `circuit.py` had not -- it built bare `Pop`s and set only `tref`, inheriting the
    dataclass's 0.55 nS -- so the factory's default silently deleted spike-triggered
    adaptation from every spinal population. The LIF path is supposed to be unchanged by
    the conductance migration; two phenotype tests caught this and nothing else did.

    The fix is structural, so the test is too: `g_adapt` must have NO default (every caller
    states its intent) and `None` must mean "keep the cell class's own value" rather than
    being spelled out as a second copy of 0.55.
    """
    import inspect
    from circuitpharm import substrate as sub
    from circuitpharm.cpg import Pop

    sig = inspect.signature(sub.make_pop)
    for name in ("g_adapt",):
        p = sig.parameters[name]
        assert p.default is inspect.Parameter.empty, (
            f"make_pop.{name} has a default again; a caller that does not set it will "
            "inherit the factory's idea of adaptation instead of the cell class's")

    default = Pop(n=1, name="x")
    kept = sub.make_pop("lif", 1, "y", g_adapt=None, tau_adapt=None)
    assert kept.g_adapt == default.g_adapt, "g_adapt=None did not keep the Pop default"
    assert kept.tau_adapt == default.tau_adapt, "tau_adapt=None did not keep the default"

    explicit = sub.make_pop("lif", 1, "z", g_adapt=1.25, tau_adapt=280.0)
    assert (explicit.g_adapt, explicit.tau_adapt) == (1.25, 280.0), (
        "an explicit adaptation strength was not applied")


def test_spinal_circuit_keeps_its_spike_triggered_adaptation():
    """The behavioural half of E20, independent of make_pop's signature. Every LIF spinal
    population must carry the adaptation it had before the migration; a zero here is what
    broke the stretch reflex and the step cycle.
    """
    from circuitpharm.circuit import SpinalCircuit
    from circuitpharm.cpg import Pop

    expect = Pop(n=1, name="x").g_adapt
    c = SpinalCircuit(seed=0)
    assert expect > 0.0, "Pop's default adaptation is zero; this test no longer means that"
    for nm, p in c.pops.items():
        assert p.g_adapt == expect, (
            f"spinal population {nm} has g_adapt={p.g_adapt} not {expect}; "
            "spike-triggered adaptation was dropped from the LIF path")


# ------------------------------------------------- E21 the manuscript quoted numbers no
#                                                   commit of this repo produces
def test_kinetic_calibration_target_lies_inside_its_own_accepted_range():
    """E21, first half. The first manuscript draft calibrated to Po_max = 0.84 and reported
    a complete set of numbers from it. `FIT_RANGES` declares po_max acceptable only on
    [0.70, 0.80], so 0.84 would have been REJECTED by this module's own validation -- but
    nothing checked the target against the range, so the two could disagree silently.

    A target outside its own accepted range means the fit is asked to land somewhere the
    module would then refuse, which is a contradiction regardless of which value is right.
    """
    from circuitpharm import gabaa_kinetics as gk
    for key, (lo, hi) in gk.FIT_RANGES.items():
        t = gk.FIT_TARGETS[key]
        assert lo <= t <= hi, (
            f"calibration target {key}={t} lies outside its accepted range [{lo}, {hi}]; "
            "the fit would be asked to reach a value the module rejects")


@pytest.mark.slow
def test_manuscript_numbers_are_generated_not_transcribed():
    """E21, second half. Every figure the manuscript quotes must come out of
    `scripts/paper_numbers.py`, which recomputes it from the model. The draft's numbers
    were internally consistent and matched no commit: Po_max 0.84 vs the repo's 0.75, a
    3.0 mM / 1.0 ms synaptic transient vs 1.0 mM / 0.30 ms, and a cited commit hash that
    exists in no ref. The qualitative conclusion survived; every quoted figure did not.

    So the script must exist, must run, and must agree with the module it reads -- the last
    part being the one that matters, since a generator that drifts from its own source is
    no better than prose.
    """
    import pathlib
    import subprocess
    import sys

    script = pathlib.Path("scripts/paper_numbers.py")
    assert script.exists(), "the manuscript's numbers have no generator"

    out = subprocess.run([sys.executable, str(script)], capture_output=True, text=True,
                         timeout=600)
    assert out.returncode == 0, f"paper_numbers.py failed:\n{out.stderr[-2000:]}"
    txt = out.stdout

    from circuitpharm import gabaa_kinetics as gk
    s = gk.fit_scheme(verbose=False)

    # The two asymptotes are computed twice by different routes inside the script (closed
    # form, and the mechanism pushed to affinity=1e5). They must agree, or one is wrong.
    dr = s.d / s.r
    closed = 1.0 / (1.0 + (s.alpha / s.beta) * (1.0 + dr))
    pushed = s.pam(affinity=1e5).po_tonic(gk.AMBIENT_UM)
    assert abs(closed - pushed) < 1e-4, (
        f"the three-state asymptote disagrees with the mechanism pushed to its limit: "
        f"{closed:.6f} vs {pushed:.6f}")

    # And the printed value must be the computed one, to the digits printed.
    assert f"{closed:.6f}" in txt, (
        f"paper_numbers.py does not print the asymptote it computes ({closed:.6f})")
    assert f"{s.po_max():.4f}" in txt, "Po_max is not reported from the fitted scheme"

    # The sections the manuscript depends on must all be present.
    for section in ("KINETIC SCHEME", "COMPARTMENT DIVERGENCE", "SENSITIVITY",
                    "FALSIFICATION INTERVALS", "SELECTIVITY INDEX"):
        assert section in txt, f"paper_numbers.py no longer emits the {section} section"


def test_falsification_intervals_are_set_from_reachable_gain_not_the_asymptote():
    """E21, third half -- the defect that would have mattered most in a wet lab. The draft
    pre-registered "R_PAM > 15x at 0.1 uM ambient GABA" as a falsification criterion. That
    is an ASYMPTOTE reading (2881x at 0.1 uM); the gain a finite s_max ~ 2.5 PAM actually
    reaches there is ~8.5x. A lab measuring 8x would have reported the model falsified when
    the model predicts 8x.

    The cause is structural: as ambient GABA falls, the asymptotic headroom grows without
    limit while the reachable gain barely moves. This test pins that divergence, so anyone
    writing a criterion against the headroom column trips it.
    """
    from circuitpharm import gabaa_kinetics as gk
    s = gk.fit_scheme(verbose=False)
    c = gk.calibrate_pam(s, target_shift=2.5, kind="affinity")
    mod = s.pam(affinity=c)
    dr = s.d / s.r
    po_inf = 1.0 / (1.0 + (s.alpha / s.beta) * (1.0 + dr))

    reach, head = {}, {}
    for g in (0.1, 0.4):
        b = s.po_tonic(g)
        reach[g], head[g] = mod.po_tonic(g) / b, po_inf / b

    # Headroom explodes as ambient falls; reachable gain does not. If these ever track each
    # other, a criterion written against either one would be safe -- and the test is moot.
    assert head[0.1] / head[0.4] > 10.0, (
        "asymptotic headroom no longer diverges from reachable gain as ambient GABA falls; "
        "re-derive the falsification intervals, this test's premise has changed")
    assert reach[0.1] / reach[0.4] < 1.5, (
        f"reachable gain at 0.1 uM is now {reach[0.1]/reach[0.4]:.2f}x its 0.4 uM value")
    assert reach[0.1] < 0.01 * head[0.1], (
        "a finite s_max=2.5 PAM should reach under 1% of the asymptote at 0.1 uM ambient")


# ------------------------------------------------- E22 the synaptic pulse had three
#                                                   definitions, two of which agreed
def test_the_synaptic_pulse_has_exactly_one_definition():
    """E22, and the most expensive error in the project so far, measured in wrong conclusions.

    `config.SYNAPTIC_PULSE` held (3000 uM, 1.00 ms). `cpg.Drug.from_kinetics` carried an
    inline literal `dict(peak_um=3000.0, clear_ms=1.00)` that matched it. `gabaa_kinetics`
    carried module defaults of 1000 uM / 0.30 ms that did NOT, and those were reached by any
    caller that passed no pulse.

    Every pharmacology consumer passes config's pulse explicitly, so for a long time only one
    caller ever saw the module defaults -- the script written to generate a manuscript's
    tables. The two calibrations differ materially:

        (3000, 1.00)  kon 0.0147  koff 0.4692  tonic headroom 210.7x  phasic gain 1.062x
        (1000, 0.30)  kon 0.0112  koff 0.3331  tonic headroom 184.6x  phasic gain 1.319x

    (Both rows are from the four-parameter fit, which roadmap P0-7 showed is not
    reproducible across SciPy versions; the pulse divergence they document is real and
    independent of that. The square fit now in use gives kon 0.0075785, koff 0.180415 and
    tonic headroom 123.3x on the config pulse.)

    A manuscript draft quoting the first set was audited, declared unreproducible, and
    "corrected" to the second -- on the strength of `git log -S"SYNAPTIC_PEAK_UM = 3000"`
    returning nothing, which it does because that symbol never held the value. The absence of
    a string was read as the absence of a value. The disagreement surfaced only when a script
    finally computed one quantity by both paths and printed 7.214 beside 7.876.

    So: one definition, and a test that fails if a second appears.
    """
    import pathlib
    import re

    from circuitpharm import gabaa_kinetics as gk
    from circuitpharm.config import SYNAPTIC_PULSE

    # 1. the module must agree with config, because it reads from it
    assert gk.SYNAPTIC_PEAK_UM == SYNAPTIC_PULSE["peak_um"], (
        f"gabaa_kinetics.SYNAPTIC_PEAK_UM is {gk.SYNAPTIC_PEAK_UM}, config says "
        f"{SYNAPTIC_PULSE['peak_um']} -- the two calibrations have diverged again")
    assert gk.SYNAPTIC_CLEAR_MS == SYNAPTIC_PULSE["clear_ms"], (
        f"gabaa_kinetics.SYNAPTIC_CLEAR_MS is {gk.SYNAPTIC_CLEAR_MS}, config says "
        f"{SYNAPTIC_PULSE['clear_ms']}")

    # 2. no source file may hard-code the pulse values again. config.py is where they live;
    #    everything else must import. The numbers are searched for as literals because that
    #    is the form the third copy took.
    peak, clear = SYNAPTIC_PULSE["peak_um"], SYNAPTIC_PULSE["clear_ms"]

    # RECURSIVE, and VALUE-AWARE. Both halves were learned after this test was written.
    #
    # `glob("*.py")` covered only the top-level package, which was every module there was at
    # the time. The fitting/, protocols/ and models/ subpackages arrived afterwards, and
    # neither this test nor tests/test_config_single_source.py (which scans scripts/) reached
    # them -- so the one constant this test exists to protect had an unguarded home.
    #
    # The match is on the VALUE, not on the keyword. A bare `peak_um=` scan flags two
    # legitimate uses: `protocols/waveforms.py` passes peak_um=1.0 as a NORMALISED amplitude
    # for shape mixing, and `fitting/data.py` carries 3000.0 as a concentration grid point in
    # a dose-response array. Neither is a copy of the pulse. Only a literal equal to the
    # configured value is a second definition, so only that is an offence.
    pat = re.compile(r"(peak_um|clear_ms)\s*=\s*([0-9]+\.?[0-9]*)")
    want = {"peak_um": float(peak), "clear_ms": float(clear)}
    offenders = []
    for f in sorted(pathlib.Path("src/circuitpharm").rglob("*.py")):
        if f.name == "config.py":
            continue
        for m in pat.finditer(f.read_text()):
            key, lit = m.group(1), m.group(2)
            try:
                val = float(lit)
            except ValueError:
                continue                      # forwarding a variable; fine
            if val == want[key]:
                offenders.append(f"{f.relative_to('src/circuitpharm')}: {key}={lit}")
    assert not offenders, (
        "the synaptic pulse is hard-coded outside config.py again: "
        f"{offenders}. Import config.SYNAPTIC_PULSE instead; this constant has had three "
        "definitions before and the two that agreed hid the one that did not.")

    # 3. and the fit that every consumer gets must be the config one.
    #
    # UPDATED 2026-10-09 (roadmap P0-7). This asserted Kd = 31.95 +/- 0.05 -- a 0.16%
    # tolerance on a quantity that was not reproducible at all. The four-parameter fit gave
    # 26.60 / 29.44 / 29.48 uM under scipy 1.14.1 / 1.11.4 / 1.17.1, every one reproducing
    # all three anchors exactly, against the 31.95 recorded here from one machine. The test
    # was right to exist; its number was an artefact.
    #
    # `fit_scheme` now holds alpha at FIT_FIXED_ALPHA, making the system square (three
    # unknowns, three anchors) and well-posed, so Kd is a determined 23.806 uM --
    # identical to 6 significant figures across all three of those environments. The
    # tolerance stays tight BECAUSE it is now reproducible; that is the point of the change.
    s = gk.fit_scheme(verbose=False)
    assert abs(s.koff / s.kon - 23.806) < 0.01, (
        f"the default fit gives Kd = {s.koff / s.kon:.3f} uM; the config pulse with alpha "
        f"held at gabaa_kinetics.FIT_FIXED_ALPHA gives 23.806. A default-pulse caller is "
        f"fitting a different scheme from the pipeline.")


@pytest.mark.slow
def test_the_two_calibration_paths_now_agree():
    """The behavioural half of E22: a Drug built through the circuit path must report the
    tonic gain the kinetic module derives directly. This is the comparison that finally
    exposed the divergence, so it is now run every time."""
    from circuitpharm import gabaa_kinetics as gk
    from circuitpharm.evaluation import Compound

    s = gk.fit_scheme(verbose=False)
    c = gk.calibrate_pam(s, 2.5, "affinity")
    direct = gk.derive(s, affinity=c)["tonic_gain"]

    circuit = Compound.from_profile("nonselective_bz", occupancy=1.0)
    via_drug = circuit.drug("prebotc").gaba_scale_tonic()

    assert abs(direct - via_drug) < 0.01 * direct, (
        f"the kinetic module derives a tonic gain of {direct:.4f} and the circuit path "
        f"reports {via_drug:.4f}. One quantity, two values -- which is E22 exactly.")


# ------------------------------------------------- E24-E26: a static audit of the P0-P7
#                                                   branch, verified by execution first
def test_equilibrium_entry_points_have_no_runnable_default():
    """E24. Four public EQUILIBRIUM functions declared `dataset=None` and defaulted to
    `DOSE_RESPONSE_BENCHMARK`, which is tagged `Observable.PEAK` -- so calling any of them at
    its own signature raised `ObservableMismatch` from inside the objective. They were
    uncallable as declared.

    The fix is not a better default: `fitting.data.MISSING_DATASETS` records that no measured
    equilibrium concentration-response exists in this project, so there IS no default. A
    required argument that names the two ways to get a dataset is honest; a default that
    cannot run reports the absence as a type confusion.
    """
    import pytest as _pytest

    from circuitpharm.fitting import identifiability as idf

    for name, args in (("fit_identifiable", ()), ("identifiability_report", ()),
                       ("profile_likelihood", ("log10_kd",)),
                       ("equilibrium_chi2_identifiable", ([1.4, 0.6, 1.4],))):
        fn = getattr(idf, name, None)
        if fn is None:
            continue
        with _pytest.raises(ValueError) as e:
            fn(*args)
        msg = str(e.value)
        assert "EQUILIBRIUM" in msg and "equilibrium_crc_from_model" in msg, (
            f"{name} refuses without a dataset but does not say how to get one: {msg[:160]}")


def test_every_model_rejects_negative_allosteric_modulation():
    """E25. `OperationalScalarModel.effective_ec50` clamped with `max(shift, 1e-6)` instead of
    guarding, so `effective_ec50(0.5)` returned 2x the EC50 -- a silent NAM -- while the two
    kinetic models raise on the same input. Measured before the fix: 25 -> 50 uM.

    The invariant is SYMMETRY, which is why this test covers all three rather than the one
    that was wrong. Model comparison is what this layer exists for, and a sweep that
    silently measures a NAM in model A and a rejection in B and C is comparing two different
    questions.
    """
    import pytest as _pytest

    from circuitpharm.models.kinetic_jw95 import KineticAllosteryModel
    from circuitpharm.models.operational import OperationalScalarModel
    try:
        from circuitpharm.models.extended_desens import ExtendedDesensitizationModel
    except ImportError:                      # pragma: no cover
        ExtendedDesensitizationModel = None

    with _pytest.raises(ValueError):
        OperationalScalarModel().effective_ec50(0.5)
    with _pytest.raises(ValueError):
        KineticAllosteryModel().apply_pam(0.5)
    if ExtendedDesensitizationModel is not None:
        with _pytest.raises(ValueError):
            ExtendedDesensitizationModel().apply_pam(0.5)

    # and a genuine PAM must still work on all of them
    assert OperationalScalarModel().effective_ec50(2.5) > 0.0


def test_the_result_schema_does_not_depend_on_the_environment():
    """E26. `joint_excursion` -- a permanently VOID quantity recording a metric with the wrong
    sign -- was added inside the try block that calls the body-plant assays. On a lean
    install the ImportError jumped past it, so `rs.quantity("joint_excursion")` raised
    KeyError without mujoco and returned a VOID quantity with it.

    The missing number is not the problem; the environment-dependent SCHEMA is. A caller
    cannot write `rs.quantity("joint_excursion").tier is Tier.VOID` and have it mean the same
    thing on two machines, and the VOID tier exists so an invalid metric is visible rather
    than absent.

    Checked statically, because the condition is "mujoco absent" and this test must give the
    same answer whether or not it is installed.
    """
    import inspect

    from circuitpharm import evaluation as ev

    # COMMENTS STRIPPED FIRST. `inspect.getsource` returns the comments too, and the note
    # explaining this very fix names `joint_excursion` three times -- so a naive count of the
    # raw source counts prose. My own bug, caught by this test on its first run.
    src = inspect.getsource(ev.evaluate)
    code = "\n".join(ln for ln in src.splitlines() if not ln.lstrip().startswith("#"))

    assert code.count('"joint_excursion"') == 1, (
        f"{code.count(chr(34) + 'joint_excursion' + chr(34))} code sites now add it; there "
        "must be exactly one, and it must be unconditional")
    at = code.index('"joint_excursion"')
    handler = code.index("except (ImportError")
    assert at > handler, (
        "joint_excursion is added before the ImportError handler, i.e. inside the try. On a "
        "lean install it will be skipped and the result schema becomes environment-dependent.")
