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
def test_cond_substrate_refuses_to_run_without_its_own_operating_point(monkeypatch):
    """RESP_OP's `drive` is in pA and its weights in nS relative to the LIF cell
    (C=200 pF, g_L=10 nS). The Butera cell is C=21 pF, g_L=2.8 nS. Reusing those numbers
    does not raise on its own -- it produces a network that does not oscillate, or
    oscillates for the wrong reason. So the refusal has to be explicit.

    THE GUARD IS EXERCISED, not inferred from an empty registry. This test used to assert
    that `PreBotC(substrate="cond")` raises, which held only while COND_RESP_OP was None.
    Registering an anchored operating point made it fail -- correctly, since the guard's
    precondition no longer existed, but it meant the test had been checking the registry's
    state rather than the guard. COND_RESP_OP is now blanked for the duration so the refusal
    path itself is tested, whether or not an operating point happens to be registered.
    """
    import circuitpharm.config as cfg
    monkeypatch.setattr(cfg, "COND_RESP_OP", None, raising=False)
    with pytest.raises(ValueError, match="needs its own operating point"):
        PreBotC(substrate="cond", seed=0)


def test_cond_substrate_runs_from_the_registered_operating_point():
    """The other side of the same guard: with one registered, no explicit `op` is needed."""
    from circuitpharm.config import COND_RESP_OP
    assert COND_RESP_OP is not None, "no conductance operating point is registered"
    net = PreBotC(substrate="cond", seed=0)
    assert net.drive == COND_RESP_OP["drive"]
    assert net.drive_other == COND_RESP_OP["drive_other"]
    assert net.gaba_tonic == COND_RESP_OP["gaba_tonic"]
    # the NMDA weights must stay ~11x smaller than AMPA relative to the LIF table: that
    # ratio is the Mg-relief correction, and losing it reintroduces the runaway into
    # depolarisation block that made the network unanchorable.
    w = COND_RESP_OP["w"]
    assert w["ee_nmda"] / w["ee_ampa"] < 0.10, (
        f"ee_nmda/ee_ampa is {w['ee_nmda']/w['ee_ampa']:.3f}; the NMDA weights are no "
        f"longer scaled down by the Mg2+ relief ratio and the network will run away into "
        f"depolarisation block")


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
def test_cond_substrate_passes_nmda_raw_so_the_block_tracks_voltage(monkeypatch):
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

    # CHECKED AT THE SEAM, BEHAVIOURALLY. This used to assert on the source text of
    # PreBotC.step and broke the moment that logic moved into the shared substrate helper --
    # a legitimate refactor failing a test that was pinning an implementation rather than a
    # behaviour. What matters is what the CELL receives, so capture it.
    captured = {}
    orig = type(p).step

    def spy(self, dt, g, E, Idrive, rng, raw=()):
        # Keyed by population name. Guarding on "nmda" not in `captured` while storing
        # under "g_nmda" let every population overwrite the last, so this captured Out
        # (g=0) instead of Exc (g=3) -- my own bug, and the test failed for it rather than
        # for the defect it is about.
        if "nmda" in g:
            captured[self.name] = (float(np.asarray(g["nmda"]).ravel()[0]), tuple(raw))
        return orig(self, dt, g, E, Idrive, rng, raw=raw)

    monkeypatch.setattr(type(p), "step", spy)
    net.step(0.1)
    g_exc, raw_exc = captured["Exc0"]
    assert raw_exc == ("nmda",), (
        f"the cell was handed raw={raw_exc}; NMDA is not being passed unevaluated, so the "
        f"Mg block stays frozen at the pre-step voltage")
    assert g_exc == pytest.approx(3.0, rel=1e-9), (
        f"the Exc cell received g_nmda={g_exc:.4f} rather than the raw 3.0, so "
        f"Syn.conductance attenuated it before handover")


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


# ================= the three substrate-coupled settings must not drift apart
def test_substrate_table_couples_operating_point_band_and_settling():
    """Three things change with the substrate and each is silent if wrong.

    * operating point -- pA/nS are relative to the cell, so RESP_OP must not be reused
    * validity band   -- EUPNOEA_BAND is IN VIVO rat (1-2 Hz); the Butera cell is neonatal
                         rodent IN VITRO (~0.11-0.24 Hz control burst frequency), so gating
                         a cond run against the in vivo band reports every healthy rhythm
                         as dead
    * settling time   -- tau_h is 10 s, 25x the LIF's tau_adapt. A 4 s warm-up is 0.4 time
                         constants and measures a decaying transient; that produced a
                         confident, seed-verified, WRONG operating point once already.

    Kept in one table so they cannot drift apart, and pinned here because each was wrong
    independently at some point during the migration.
    """
    from circuitpharm.config import INVITRO_BAND
    from circuitpharm.evaluation import _SUBSTRATE
    lif, cond = _SUBSTRATE["lif"], _SUBSTRATE["cond"]

    assert lif["band"] is None, "the LIF must keep the in vivo default band"
    assert lif["warm_ms"] == 4000.0 and lif["duration_ms"] == 14000.0, (
        "the LIF timings must not change; every published result used them")

    assert cond["band"] == "invitro"
    # at least 3 tau_h of settling
    assert cond["warm_ms"] >= 3 * 10000.0, (
        f"cond warm-up {cond['warm_ms']} ms is under 3 tau_h (30 s); the slow I_NaP "
        f"inactivation gate has not equilibrated and the measurement is of a transient")
    assert cond["duration_ms"] - cond["warm_ms"] >= 20000.0, (
        "less than 20 s of analysable window after settling")
    assert INVITRO_BAND[0] < 0.30 <= INVITRO_BAND[1], (
        "the in vitro band must reach below the in vivo floor but still bound the top")


def test_the_in_vitro_band_still_excludes_fragmented_rhythms():
    """The band exists to catch disintegrated rhythms (recurring error E4, second form:
    under heavy block the burst train fragments and the FFT correctly reports ~4 Hz, which
    is not tachypnoea). Lowering the floor for the in vitro preparation must not raise the
    ceiling past that."""
    from circuitpharm.config import EUPNOEA_BAND, INVITRO_BAND
    assert INVITRO_BAND[1] <= EUPNOEA_BAND[1], (
        "the in vitro band's ceiling exceeds the in vivo one, so it would admit the "
        "fragmented high-frequency rhythms the gate exists to reject")
    for frag_hz in (3.8, 4.9):
        assert not (INVITRO_BAND[0] <= frag_hz <= INVITRO_BAND[1])


def test_the_in_vitro_band_is_sourced_not_convenient():
    """A band introduced to make a failing point qualify would be indistinguishable in code
    from one introduced because the preparation differs. The difference is the record: the
    source must be in the knowledge base with a hand-read verdict."""
    import sqlite3
    import pathlib
    db = pathlib.Path(__file__).resolve().parent.parent / "data" / "pharmacology.db"
    if not db.exists():
        pytest.skip("no database in this checkout")
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    row = con.execute("select citation, verification from sources "
                      "where key='pbc_invitro_freq'").fetchone()
    con.close()
    assert row, "the in vitro frequency source is not in the knowledge base"
    assert "HAND:" in row[1] and "SUPPORTS" in row[1], (
        f"the in vitro band's source carries no hand-read claim-support verdict: {row[1][:90]}")


# =============== the shared substrate helper, and the two newly-wired circuits
def test_derive_weights_separates_nmda_from_everything_else():
    """The correction that took four failed sweeps to find. A single scale cannot carry a
    weight table between these cells: NMDA needs the Mg-relief ratio on top of the g_L
    ratio, making it ~11x smaller than AMPA.

    tau_nmda is 100 ms, 20x the AMPA tau, so lumped population NMDA accumulates 20x more
    conductance per unit firing rate. The LIF's weights were tuned where Mg block held NMDA
    at a measured mean relief of 0.063; on a cell that relieves to ~0.70 the same weight is
    ~11x too strong -- and since relief RISES with depolarisation it is positive feedback
    into depolarisation block, not a scale error.
    """
    from circuitpharm.substrate import GL_RATIO, RELIEF_RATIO, derive_weights
    lif = dict(ee_ampa=1.0, ee_nmda=1.0, ie_gaba=1.0, ie_gly=1.0)
    w = derive_weights(lif)
    assert w["ee_ampa"] == pytest.approx(GL_RATIO)
    assert w["ie_gaba"] == pytest.approx(GL_RATIO)
    assert w["ie_gly"] == pytest.approx(GL_RATIO)
    assert w["ee_nmda"] == pytest.approx(GL_RATIO * RELIEF_RATIO)
    ratio = w["ee_ampa"] / w["ee_nmda"]
    assert 9.0 < ratio < 13.0, (
        f"AMPA/NMDA scale ratio is {ratio:.1f}, not ~11. Losing this separation is what "
        f"made the network unanchorable through three consecutive sweeps.")


def test_minimum_settle_is_three_tau_h():
    from circuitpharm.substrate import MIN_SETTLE_MS, TAU_H_MS
    assert TAU_H_MS == 10000.0
    assert MIN_SETTLE_MS >= 3 * TAU_H_MS


@pytest.mark.parametrize("ctor,kw", [
    ("GroupPacemakerRG", {}),
    ("SpinalCircuit", {}),
])
def test_newly_wired_circuits_refuse_the_cond_substrate_unanchored(ctor, kw):
    """`rg2` and `circuit` now take `substrate=`, and both must refuse 'cond' while their
    operating points are unanchored rather than silently inheriting LIF-scaled pA and nS."""
    import circuitpharm.circuit as circuit_mod
    import circuitpharm.rg2 as rg2_mod
    cls = {"GroupPacemakerRG": rg2_mod.GroupPacemakerRG,
           "SpinalCircuit": circuit_mod.SpinalCircuit}[ctor]
    with pytest.raises(ValueError, match="needs its own operating point"):
        cls(substrate="cond", seed=0, **kw)


@pytest.mark.parametrize("ctor", ["GroupPacemakerRG", "SpinalCircuit"])
def test_newly_wired_circuits_reject_an_unknown_substrate(ctor):
    import circuitpharm.circuit as circuit_mod
    import circuitpharm.rg2 as rg2_mod
    cls = {"GroupPacemakerRG": rg2_mod.GroupPacemakerRG,
           "SpinalCircuit": circuit_mod.SpinalCircuit}[ctor]
    with pytest.raises(ValueError, match="must be 'lif' or 'cond'"):
        cls(substrate="definitely-not-a-substrate", seed=0)


def test_locomotor_interneuron_bias_is_no_longer_a_literal():
    """`drive_in` was the literal 70.0 in rg2's step loop -- pA against the LIF cell, and
    the third LIF-scaled quantity found hiding in a step loop after the drive and the weight
    table had both been moved out. resp.py had the identical trap at 60.0."""
    import inspect
    from circuitpharm.rg2 import GroupPacemakerRG
    src = inspect.getsource(GroupPacemakerRG.step)
    assert "else 70.0" not in src
    assert "self.drive_in" in src
    assert GroupPacemakerRG(seed=0).drive_in == 70.0, "the LIF default must not change"


def test_spinal_circuit_forwards_the_substrate_to_its_rhythm_generator():
    """Review 4 found SpinalCircuit failing to forward sensitivity into its RG. The same
    omission for `substrate` would leave a conductance spinal circuit driven by a LIF
    rhythm generator -- a mixed-substrate network, which is worse than either."""
    import inspect
    from circuitpharm.circuit import SpinalCircuit
    src = inspect.getsource(SpinalCircuit.__init__)
    assert "substrate=substrate" in src, "the substrate is not forwarded into the RG"
    assert SpinalCircuit(seed=0).rg.substrate == "lif"


def test_all_three_circuits_share_one_substrate_implementation():
    """Recurring error E12 was a fix applied everywhere except the module written later,
    where the same defect reappeared. Three circuits now route through `substrate.py` rather
    than carrying three copies of the logic."""
    import inspect
    from circuitpharm import circuit, resp, rg2
    for mod in (resp, rg2, circuit):
        src = inspect.getsource(mod)
        assert "substrate as sub" in src, f"{mod.__name__} does not use the shared helper"
        assert "sub.make_pop(" in src, f"{mod.__name__} builds populations directly"
        assert "sub.raw_receptors(" in src, f"{mod.__name__} does not use shared raw-NMDA"
        assert "sub.tonic_gaba(" in src, f"{mod.__name__} inlines the tonic GABA clamp"
