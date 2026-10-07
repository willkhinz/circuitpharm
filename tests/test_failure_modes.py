"""Regression tests for the project's catalogued failure modes (E1-E12).

WHY THIS FILE IS THE HEART OF THE TEST SUITE. Each of these errors was made at least once
during development, several were made twice, and every one of them produced a RESULT THAT
LOOKED FINE -- a flat parameter sweep, a plausible percentage, a reassuring safety margin.
None announced themselves as bugs. Documenting them in a worklog was not enough, because a
documented failure mode can still silently return.

So each error is an executable test here. Where possible the test is a fast pure-unit check
rather than a full circuit integration, and where a two-sided test is available (the
mechanism is engaged AND would disengage if the error returned) it is preferred, because a
one-sided assertion can pass for the wrong reason.

The full catalogue with symptoms and fixes is in WORKLOG.md.
"""
import numpy as np
import pytest

from circuitpharm.cpg import Drug, Pop, Syn, TAU
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.config import RESP_OP, EUPNOEA_BAND, CALIBRATIONS
from circuitpharm.results import Tier


def _run(duration_ms=8000.0, warm_ms=3000.0, **kw):
    b = PreBotC(seed=0, **{**RESP_OP, **kw})
    for i in range(int(duration_ms / 0.1)):
        b.step(0.1)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > warm_ms
    return b, resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])


# ===================================================================== E1
# Suprathreshold tonic drive in a network that must burst.
# SYMPTOM: population fires tonically; adaptation and mutual inhibition have no effect;
# every weight sweep comes out flat. Made TWICE (locomotor RG, then preBotC).
@pytest.mark.slow
def test_E1_operating_point_produces_a_real_rhythm():
    """Two-sided: the configured drive must BURST, and a supra-rheobase drive must not.

    A one-sided check that the rhythm exists would still pass if the drive crept upward,
    because a tonically firing population has a nonzero mean output.
    """
    _, r = _run()
    assert r["alive"], f"rhythm absent at the configured operating point: {r['reason']}"
    assert r["mod"] > 1.0, f"modulation {r['mod']:.2f} too low to be bursting"


@pytest.mark.slow
def test_E1_suprarheobase_drive_destroys_bursting():
    """The other half of E1: if drive is pushed well above rheobase the rhythm must die.

    This is what makes the test above meaningful -- it proves the rhythm depends on the
    drive being sub-rheobase rather than being an artefact of the metric.
    """
    _, r = _run(drive=600.0)
    assert not r["alive"] or r["mod"] < 1.0, (
        "a strongly supra-rheobase drive still reported a healthy rhythm; the burst "
        f"mechanism is not load-bearing (mod={r['mod']:.2f})")


# ===================================================================== E3
# Interneurons left without tonic drive sit just below threshold and never fire, so the
# pathway they mediate silently does nothing and its weight sweep looks flat.
@pytest.mark.slow
def test_E3_inhibitory_population_actually_fires():
    b, _ = _run(duration_ms=4000.0, warm_ms=1000.0)
    assert b.pops["Inh"].rate >= 0.0
    # the diagnostic that matters: Inh must have spiked at all during the run
    assert b.trace["Inh"] and max(b.trace["Inh"]) > 0.0, (
        "the inhibitory population never fired, so inhibition is not engaged and any "
        "sweep over inhibitory weights would be meaningless")


# ===================================================================== E4
# A relative threshold on a degrading signal reports respiratory DEPRESSION as
# tachypnoea. First form: frac*max() crossings gave 13 Hz (780 breaths/min). Second form
# (found later): a genuinely fragmented burst train has real spectral power at ~4 Hz, which
# the FFT reports correctly and which `mod > 0.8 and n >= 3` wrongly passed as ALIVE.
# This mattered because `alive` is what the overdose scan keys on.
def test_E4_fragmented_fast_rhythm_is_not_alive():
    """Synthetic 4 Hz train with good modulation must be rejected as non-physiological."""
    t = np.arange(0, 10000.0, 1.0)
    out = 20.0 * (1.0 + 0.9 * np.sin(2 * np.pi * 4.0 * t / 1000.0))
    r = resp_metrics(t, out, out)
    assert r["freq"] > EUPNOEA_BAND[1]
    assert r["frag"] is True
    assert not r["alive"], (
        "a 4 Hz (240 breaths/min) train passed as a live respiratory rhythm; this is the "
        "bug that made overdose indices optimistic")
    assert "eupnoea" in r["reason"] or "frequency" in r["reason"]


def test_E4_physiological_rhythm_is_alive():
    """The complement: a clean 1.3 Hz rhythm must pass. Guards against a gate so strict
    it rejects everything, which would hide drug effects instead of revealing them."""
    t = np.arange(0, 10000.0, 1.0)
    out = 20.0 * (1.0 + 0.9 * np.sin(2 * np.pi * 1.3 * t / 1000.0))
    r = resp_metrics(t, out, out)
    assert r["alive"], r["reason"]
    assert r["freq"] == pytest.approx(1.3, abs=0.15)
    assert r["frag"] is False


def test_E4_apnoea_is_detected_as_mean_collapse():
    t = np.arange(0, 10000.0, 1.0)
    r = resp_metrics(t, np.zeros_like(t), np.zeros_like(t))
    assert not r["alive"]
    assert "collapsed" in r["reason"]


# ===================================================================== E7
# Voltage-gated gates do not work on LIF neurons: V is reset at threshold and never
# reaches spike voltages, so a slow voltage-gated gate equilibrates at the subthreshold
# mean and never moves (observed h-gate span 0.01 against a needed 0.7).
def test_E7_lif_voltage_never_reaches_spike_range():
    """Documents the structural reason burst termination uses spike-triggered adaptation
    rather than voltage-gated inactivation. If this ever fails, voltage-gated mechanisms
    become available and the modelling choice should be revisited."""
    p = Pop(n=20, name="E7probe", nap=True)
    rng = np.random.default_rng(0)
    vmax = -1e9
    for _ in range(20000):
        p.step(0.1, {}, {}, 180.0, rng)
        vmax = max(vmax, float(p.V.max()))
    assert vmax <= p.Vth + 1e-9, (
        f"LIF voltage reached {vmax:.2f} mV, above threshold {p.Vth}; the reset is not "
        "being applied as assumed")


def test_E7_spike_triggered_adaptation_accumulates_and_decays():
    p = Pop(n=10, name="E7adapt")
    p.g_adapt, p.tau_adapt = 1.0, 100.0
    rng = np.random.default_rng(0)
    for _ in range(2000):
        p.step(0.1, {}, {}, 400.0, rng)
    assert p.a.max() > 0.0, "adaptation never accumulated despite sustained firing"
    peak = float(p.a.max())
    for _ in range(5000):          # no drive: adaptation must decay away
        p.step(0.1, {}, {}, 0.0, rng)
    assert float(p.a.max()) < 0.5 * peak


# ===================================================================== E11
# Scaling a drug effect on conductance but NOT on kinetics floors the achievable
# sensitivity: ventilation stuck at -25% at any gaba_sens, because tau stayed fully
# prolonged however small the conductance share.
def test_E11_sens_scales_both_conductance_and_tau():
    """A PAM changes peak conductance AND decay tau; `sens` must scale both identically."""
    # cap raised out of the way: this test is about `sens` scaling, and the efficacy cap
    # (default 2.5) would otherwise clip the gain and mask it. Capping is tested separately
    # in test_identities::test_efficacy_caps_bind_per_pool.
    drug = Drug(gaba_a_gain=3.0, gaba_a_tau=2.0, gaba_a_efficacy_cap=1e9)
    full = Syn(5, "gabaa", drug, sens=1.0)
    half = Syn(5, "gabaa", drug, sens=0.5)
    none = Syn(5, "gabaa", drug, sens=0.0)

    # conductance scaling
    assert full.w_scale == pytest.approx(3.0)
    assert half.w_scale == pytest.approx(1.0 + 0.5 * (3.0 - 1.0))
    assert none.w_scale == pytest.approx(1.0)
    # tau scaling must follow the SAME rule, not stay at the full drug value
    assert full.tau == pytest.approx(TAU["gabaa"] * 2.0)
    assert half.tau == pytest.approx(TAU["gabaa"] * (1.0 + 0.5 * (2.0 - 1.0)))
    assert none.tau == pytest.approx(TAU["gabaa"]), (
        "with sens=0 the decay tau was still prolonged; this is E11, which floored the "
        "achievable sensitivity at any gaba_sens")


def test_E11_zero_sens_is_a_complete_no_op():
    drug = Drug(gaba_a_gain=5.0, gaba_a_tau=3.0, glyr_gain=4.0, nmda_block=0.8,
                gaba_a_efficacy_cap=1e9)
    for kind in ("gabaa", "gly", "nmda"):
        s = Syn(5, kind, drug, sens=0.0)
        assert s.w_scale == pytest.approx(1.0), kind
        assert s.tau == pytest.approx(TAU[kind]), kind


# ===================================================================== E12
# A calibration constant carried across a change of parameterisation. Made twice: across
# regions, then across the tonic/phasic pool split -- the second time one worklog entry
# after writing down that the migration would be needed. Symptom both times: a
# suspiciously REASSURING result (every reflex effect abolished; every overdose survived).
def test_E12_invalidated_anchor_is_marked_void():
    """The knowledge that the preBotC anchor is invalid must live in code, not prose.

    Human whole-body ventilation is the wrong observable for an isolated preBotC, so this
    parameter was absorbing chemoreflex, airway and cortical mechanisms the circuit does
    not contain. Marking it VOID is what stops it being quoted or silently reused.
    """
    c = CALIBRATIONS["prebotc_gaba_sens"]
    assert c.tier is Tier.VOID
    assert c.promote_by, "a VOID calibration must say what would promote it"
    assert "P12" in c.promote_by, (
        "the promotion path must retain the >P12 age requirement; the preBotC NKCC1/KCC2 "
        "crossover is why a newborn preparation's GABA-A sign is uninformative")


def test_E12_uncalibrated_axes_are_not_claimed_as_validated():
    for key in ("spinal_gaba_sens", "forebrain_gaba_sens"):
        assert CALIBRATIONS[key].tier is not Tier.VALIDATED, key


def test_E12_only_calibration_independent_result_is_validated():
    """Exactly one entry should be VALIDATED: the ranking, which is a ratio and therefore
    independent of the lumped scale. If something else becomes VALIDATED, it must have
    earned it with an anchor."""
    validated = [k for k, c in CALIBRATIONS.items() if c.tier is Tier.VALIDATED]
    assert validated == ["ranking"], (
        f"unexpected VALIDATED calibrations: {validated}. Promotion requires a "
        "structure-matched anchor, not a refit.")


# ============================================================ general diagnostic
# "A parameter sweep that produces identical results across a wide range means the
# mechanism is not engaged at all." Encoded as a test: the drug must actually do something.
@pytest.mark.slow
def test_gaba_pam_has_a_monotone_effect_on_ventilation():
    """Guards the general heuristic. If a PAM sweep came out flat, a population upstream
    would be silent and every pharmacology result would be vacuous."""
    means = []
    for gain in (1.0, 2.0, 4.0):
        _, r = _run(duration_ms=8000.0,
                    drug=Drug(gaba_a_gain=gain, gaba_a_tau=1.0 + 0.6 * (gain - 1.0)),
                    gaba_sens=0.5)
        means.append(r["mean"])
    assert means[0] > means[-1], (
        f"increasing GABA-A potentiation did not reduce ventilation ({means}); the drug "
        "node is not engaged")
    assert abs(means[0] - means[-1]) / means[0] > 0.05, (
        f"PAM sweep was effectively flat ({means}) -- suspect a silent population")
