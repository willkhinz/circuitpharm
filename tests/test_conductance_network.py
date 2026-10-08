"""Link 5: the conductance cell wired into the preBotC network.

These tests cover the PLUMBING, which is where this migration's predicted defects live --
every one of them silent. The scientific questions (does the rhythm survive I_NaP removal,
does the compound ordering survive the substrate change) are answered by scripts, because
they are measurements whose outcome is not known in advance and a test that asserts a
hypothesis is not a test.
"""
import hashlib

import numpy as np
import pytest

from circuitpharm.config import RESP_OP
from circuitpharm.resp import PreBotC, REQUIRED_W, resp_metrics

# A known-alive conductance operating point (stage-2 output of scripts/anchor_cond_resp.py:
# weight-block scale 0.9 x the LIF table). Used so these tests work whether or not
# config.COND_RESP_OP has been registered yet, and so they do not silently start testing a
# different operating point when it is.
OP = dict(
    drive=20.0, drive_other=30.0, gaba_tonic=1.2,
    w={"ee_ampa": 0.405, "ee_nmda": 0.22275, "ei_ampa": 0.495, "ei_nmda": 0.0,
       "ie_gaba": 0.405, "ie_gly": 0.315, "eo_ampa": 0.63, "eo_nmda": 0.27},
)


def _run(net, seconds=6.0, warm_ms=2000.0):
    for i in range(int(round(seconds * 1000.0 / 0.1))):
        net.step(0.1)
        if i % 10 == 0:
            net.record()
    A = net.arrays()
    m = A["t"] > warm_ms
    return resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])


# ===================================================== the LIF path must not have moved
def test_lif_substrate_is_byte_identical_after_the_refactor():
    """Adding the substrate switch restructured population construction and the step loop.
    The LIF path must be untouched, because every published result in this project came from
    it AND because the substrate comparison is meaningless against a drifting baseline.

    Hash recorded against the committed pre-refactor module, which produced the same digest.
    """
    net = PreBotC(**RESP_OP, seed=0)
    tr = []
    for i in range(40000):
        net.step(0.1)
        if i % 10 == 0:
            tr.append([net.pops[k].rate for k in ("Exc", "Inh", "Out")])
    a = np.round(np.asarray(tr), 9)
    assert hashlib.sha256(a.tobytes()).hexdigest()[:20] == "cb808accc11db8462b7f"


# ===================================================== E14: refuse LIF-scaled quantities
def test_cond_substrate_refuses_to_run_without_its_own_operating_point():
    """RESP_OP's `drive` is in pA and its weights in nS relative to the LIF cell
    (C=200 pF, g_L=10 nS). The Butera cell is C=21 pF, g_L=2.8 nS. Reusing those numbers
    does not raise on its own -- it produces a network that does not oscillate, or
    oscillates for the wrong reason. So the refusal has to be explicit."""
    with pytest.raises(ValueError, match="needs its own operating point"):
        PreBotC(substrate="cond", seed=0)


@pytest.mark.parametrize("drop", ["drive", "drive_other", "gaba_tonic", "ee_ampa",
                                  "ie_gaba", "eo_nmda"])
def test_a_partial_operating_point_is_refused(drop):
    """THE DANGEROUS CASE is a partial operating point, not a missing one. Supply `drive`
    and the excitatory weights while `ie_gaba`, `ie_gly` and `gaba_tonic` quietly keep their
    LIF values, and the network is internally inconsistent by ~an order of magnitude on the
    INHIBITORY arm -- which is the arm the drug acts through. It would run, produce a
    rhythm, and be wrong about pharmacology."""
    op = dict(OP, w=dict(OP["w"]))
    if drop in op:
        del op[drop]
    else:
        del op["w"][drop]
    with pytest.raises(ValueError, match="incomplete"):
        PreBotC(substrate="cond", op=op, seed=0)


def test_an_unknown_substrate_is_refused():
    with pytest.raises(ValueError, match="must be 'lif' or 'cond'"):
        PreBotC(substrate="definitely-not-a-substrate", seed=0)


# ========================================== E18: no double-counted burst termination
def test_cond_populations_have_no_spike_triggered_adaptation():
    """The LIF assigned `g_adapt` AFTER construction, so pointing it at a conductance cell
    would have switched adaptation back on at its LIF-tuned 2.5 nS/spike on top of a now-real
    I_NaP inactivation. Burst termination would be double-counted and the symptom would be
    bursts ending slightly early -- a plausible duty cycle, which is this project's signature
    failure shape."""
    net = PreBotC(substrate="cond", op=OP, seed=0)
    for nm, p in net.pops.items():
        assert p.g_adapt == 0.0, f"{nm} carries g_adapt={p.g_adapt} on the cond substrate"
    lif = PreBotC(**RESP_OP, seed=0)
    assert lif.pops["Exc"].g_adapt == RESP_OP["g_adapt"], (
        "the LIF path must keep its adaptation; it is that substrate's burst terminator")


# ================================= E15: NMDA handed over unevaluated on the cond substrate
def test_cond_substrate_passes_nmda_raw_so_the_block_tracks_voltage():
    """`Syn.conductance(V)` applies the Mg2+ block at the PRE-STEP voltage, which discards
    the entire voltage-dependent relief on a cell that sub-steps through a 70 mV spike. The
    symptom would be "the substrate change did not move the NMDA result" -- a reassuring
    robustness check that is in fact the bug.

    Checked at the seam: the cell must RECEIVE the unattenuated conductance, so what it is
    handed is strictly larger than what the LIF path would hand it at the same voltage.
    """
    net = PreBotC(substrate="cond", op=OP, seed=0)
    syn = net.syn[("Exc", "nmda")]
    syn.g = np.full(net.pops["Exc"].n, 3.0)
    p = net.pops["Exc"]
    raw = syn.g
    evaluated = syn.conductance(p.V)
    assert np.all(raw > evaluated), (
        "the raw conductance is not larger than the pre-evaluated one, so Syn.conductance "
        "is not attenuating and this test cannot detect the defect it is for")

    import inspect
    src = inspect.getsource(PreBotC.step)
    assert 'raw = ("nmda",) if self.substrate == "cond" else ()' in src
    assert "s.g if rec in raw else s.conductance(p.V)" in src


# ============================= the control cache must not mix substrates
def test_control_cache_is_keyed_on_the_substrate():
    """`_CTRL_CACHE`'s safety argument covers the DRUG parameters only: a drug-free run gives
    eff = 1 for any sensitivity. The neuron model is a different axis, and it changes the
    drug-free control completely -- the conductance network is anchored separately and runs
    at a different frequency and mean output by construction.

    Keyed on seed alone, a conductance arm would have been normalised against a LIF control,
    making every reported percent-of-control silently nonsense. The substrate comparison is
    precisely a comparison of those percentages, so the one number this upgrade exists to
    produce would have been the one corrupted.
    """
    import inspect
    from circuitpharm import evaluation
    src = inspect.getsource(evaluation.evaluate)
    assert '("resp", s, substrate)' in src, (
        "the control-cache key does not include the substrate")
    assert "substrate=substrate" in src


def test_evaluate_and_simulate_resp_accept_a_substrate():
    import inspect
    from circuitpharm.evaluation import _simulate_resp, evaluate
    assert inspect.signature(evaluate).parameters["substrate"].default == "lif"
    assert inspect.signature(_simulate_resp).parameters["substrate"].default == "lif"


# ============================================================ the network actually runs
@pytest.mark.slow
def test_cond_network_produces_a_living_rhythm_in_band():
    """The operating point above must give a control rhythm inside EUPNOEA_BAND. This is the
    precondition for everything else: a drug arm compared against a dead control is not a
    measurement."""
    m = _run(PreBotC(substrate="cond", op=OP, seed=0), seconds=12.0)
    assert m["alive"], f"cond control is not alive: {m['reason']}"
    assert 0.30 <= m["freq"] <= 2.50


@pytest.mark.slow
def test_cond_network_voltage_actually_spikes():
    """The whole point of the substrate: V must leave the LIF's [-71, -44] mV window inside
    the NETWORK, not just in an isolated cell. If synaptic load keeps the population
    subthreshold, the upgrade has bought nothing where it matters."""
    net = PreBotC(substrate="cond", op=OP, seed=0)
    vmax = -1e9
    for _ in range(60000):
        net.step(0.1)
        vmax = max(vmax, float(net.pops["Exc"].V.max()))
    assert vmax > 0.0, (
        f"peak V in the coupled network only reached {vmax:+.2f} mV; the conductance cell "
        f"is not spiking under synaptic load")
