"""Acceptance tests for the conductance-based cell (roadmap link 4).

These were specified in knowledge/05-design-conductance-substrate.md BEFORE the module was
written, as sections 4 (A1-A7) and 8 (E13-E18), and the names below carry those labels. The
reason for writing them first is that this project's defects are overwhelmingly "plausible
and wrong" -- 49 of them across eight review passes, none found by my own scrutiny -- and the
only defence that has worked is a test that fails loudly.

Each A-test is an observable the cell was NOT fitted to. The parameters came from the
published model; nothing here was tuned to make these pass.
"""
import hashlib

import numpy as np
import pytest
from dataclasses import replace

from circuitpharm.cpg import Pop
from circuitpharm.neuron import (
    BRS1999_MODEL1, CELLS, CondPop, CellParams, ParamSet, ProvenanceMismatch,
    isolated, run_isolated, spike_regime,
)

P = BRS1999_MODEL1


def _mg_relief(V):
    """The Mg2+ block expression shipped in cpg.Syn.conductance, so these tests measure the
    real one rather than a re-derivation."""
    return 1.0 / (1.0 + 0.28 * np.exp(-0.062 * np.asarray(V, float)))


# ===================================================================== A1: E7, inverted
def test_A1_voltage_reaches_spike_range_unlike_the_LIF():
    """E7 asserts the LIF's V never exceeds its threshold, and its docstring says: "If this
    ever fails, voltage-gated mechanisms become available and the modelling choice should be
    revisited." This is that revisit, from the other side -- the conductance cell must do
    what E7 proves the LIF cannot."""
    r = run_isolated(seconds=10.0)
    assert r["v_max"] > 0.0, (
        f"peak V {r['v_max']:+.2f} mV never crossed 0; a conductance cell that does not "
        f"produce a real spike has no advantage over the LIF it replaces")


def test_A1_slow_voltage_gated_gate_moves():
    """The h gate (I_NaP inactivation) must actually traverse its range.

    ON THE THRESHOLD, HONESTLY. The E7 note records "observed h-gate span 0.01 against a
    needed 0.7", and the design document carried that 0.7 forward as >0.5. The 0.7 was never
    sourced -- it was my own estimate of what inactivation would need, not a measurement or a
    published figure. The cell's actual span at the published operating point is ~0.45, so
    the pre-committed number would have failed a cell that is behaving correctly.

    Rather than quietly relax the bar, the bar is replaced by a FUNCTIONAL one: the next test
    shows that freezing h abolishes bursting, which is the claim the 0.7 was a proxy for and a
    far stronger check. The span assertion below is kept only as a coarse floor at 30x the
    LIF's measured 0.01, with the measured value recorded so drift is visible.
    """
    r = run_isolated(seconds=20.0)
    assert r["h_span"] >= 0.30, (
        f"h span {r['h_span']:.4f} is below 30x the LIF's measured 0.01; the slow gate is "
        f"not traversing its range")
    assert r["h_span"] == pytest.approx(0.4544, abs=0.05), (
        f"h span has drifted from the value measured when this cell was adopted "
        f"(0.4544); got {r['h_span']:.4f}")


def test_A1_h_is_what_terminates_the_burst():
    """THE STRONG FORM, and the test that supersedes the h-span threshold above.

    Clamping h at a constant abolishes bursting in every direction: held at its minimum the
    cell is quiescent (I_NaP too inactivated to fire at all), held higher it spikes TONICALLY
    at 69-119 Hz with no burst structure. So the rhythm comes from the mechanism claimed --
    slow voltage-dependent inactivation of a persistent sodium current -- and not from
    something incidental to the implementation.

    The LIF cannot pass this test, or fail it: it has no h gate to clamp. That is what makes
    this the acceptance criterion for link 4.
    """
    free = run_isolated(seconds=20.0)
    assert free["regime"] == "bursting", f"baseline cell is {free['regime']}, not bursting"

    quiet = run_isolated(seconds=20.0, freeze_h=0.46)
    assert quiet["regime"] == "quiescent", (
        f"h frozen at its minimum gave {quiet['regime']}; I_NaP should be too inactivated "
        f"to support firing")

    for hv in (0.60, 0.90):
        tonic = run_isolated(seconds=20.0, freeze_h=hv)
        assert tonic["regime"] == "tonic", (
            f"h frozen at {hv} gave {tonic['regime']}, expected tonic spiking: with no slow "
            f"inactivation there is nothing left to terminate a burst")
        assert tonic["rate_hz"] > 30.0


# ============================================== A2: bursting without lumped adaptation
def test_A2_bursts_with_no_spike_triggered_adaptation():
    """The LIF needs spike-triggered adaptation to burst at all; it was a lumped stand-in for
    I_NaP inactivation plus Ca-dependent K current plus synaptic depression. This cell must
    burst with adaptation entirely absent, which is what separates "rhythm from I_NaP" from
    "rhythm from a stand-in"."""
    c = isolated()
    assert c.g_adapt == 0.0, "adaptation must be OFF by default on this substrate (E18)"
    r = run_isolated(seconds=20.0)
    assert r["regime"] == "bursting"
    assert r["n_bursts"] >= 3 and r["spikes_per_burst"] > 2.0


def test_A2_bursting_requires_the_persistent_sodium_current():
    """Removing I_NaP must abolish ISOLATED-cell bursting -- the mechanism is load-bearing,
    not decorative.

    The other half of this prediction (A4: the NETWORK rhythm should PERSIST without I_NaP,
    the pacemaker-versus-network-rhythm controversy) needs the coupled population and is
    deferred to link 5. It is deliberately not asserted here, because asserting half a
    prediction and calling it validated is how this project produced its retracted Matsuoka
    result.
    """
    r = run_isolated(seconds=20.0, params=replace(P, g_nap=0.0))
    assert r["regime"] == "quiescent", (
        f"g_NaP = 0 gave {r['regime']} at {r['rate_hz']:.1f} Hz; bursting that survives "
        f"removal of the current it is supposed to depend on is coming from somewhere else")


# ============================================ A3 / A5: the published excitability sequence
def test_A3_excitability_sequence_quiescent_bursting_tonic():
    """Raising tonic drive must move the cell quiescent -> bursting -> tonic spiking. This is
    published behaviour of the model and was not fitted by us, so reproducing it is
    tier-promoting evidence rather than a tautology."""
    seq = {ia: run_isolated(seconds=20.0, i_app=ia)["regime"]
           for ia in (-30.0, -15.0, -5.0, 0.0, 10.0, 25.0, 50.0)}
    assert seq[-30.0] == "quiescent" and seq[-15.0] == "quiescent", seq
    assert seq[-5.0] == "bursting" and seq[0.0] == "bursting", seq
    assert seq[10.0] == "tonic" and seq[25.0] == "tonic" and seq[50.0] == "tonic", seq


def test_A5_firing_rate_increases_monotonically_with_drive():
    """f-I monotonicity across the tonic regime."""
    rates = [run_isolated(seconds=10.0, i_app=ia)["rate_hz"]
             for ia in (10.0, 25.0, 50.0, 100.0)]
    assert all(b > a for a, b in zip(rates, rates[1:])), f"f-I not monotonic: {rates}"


# =========================================================== E13: integration convergence
def test_E13_burst_period_is_converged_at_the_shipped_step():
    """An HH cell integrated with too large a step does not blow up -- it produces a
    plausible but WRONG firing rate. That is this project's signature failure mode, so the
    shipped `dt_max` must be shown converged rather than assumed.

    Measured: 0.02% period error against a step four times finer. The design document asked
    for under 2%.
    """
    fine = run_isolated(seconds=20.0, params=replace(P, dt_max=P.dt_max / 4.0))
    ship = run_isolated(seconds=20.0, params=replace(P, dt_max=P.dt_max))
    err = abs(ship["burst_period_ms"] - fine["burst_period_ms"]) / fine["burst_period_ms"]
    assert err < 0.02, (
        f"burst period moved {100*err:.2f}% between dt_max={P.dt_max} and "
        f"{P.dt_max/4}; the shipped step is not converged")
    assert ship["regime"] == fine["regime"] == "bursting"


def test_E13_substep_count_follows_dt_max():
    """The cell subdivides the caller's dt internally, so a circuit can swap substrates
    without touching its stepping loop. Two dt_max values that both yield one substep must
    therefore give IDENTICAL results -- which is why dt_max=0.2 and dt_max=0.1 agree exactly
    at an external dt of 0.1, rather than that agreement indicating dt-independence."""
    a = run_isolated(seconds=5.0, dt=0.1, params=replace(P, dt_max=0.2))
    b = run_isolated(seconds=5.0, dt=0.1, params=replace(P, dt_max=0.1))
    assert np.array_equal(a["spikes_ms"], b["spikes_ms"])
    c = run_isolated(seconds=5.0, dt=0.1, params=replace(P, dt_max=0.05))
    assert not np.array_equal(a["spikes_ms"], c["spikes_ms"]), (
        "halving the substep changed nothing at all, which means dt_max is not reaching "
        "the integrator")


# ================================================= E14: provenance locking on parameters
def test_E14_weights_from_a_different_cell_are_refused():
    """The predicted most-likely defect of this migration. The LIF cell is C=200 pF,
    g_L=10 nS; this one is C=21 pF, g_L=2.8 nS. Every synaptic weight in the package is in
    nS, hand-tuned against the LIF. Pairing them silently is an order-of-magnitude error that
    makes the network fail to oscillate, or oscillate for the wrong reason, with nothing
    raising."""
    good = ParamSet(cell=P, weights={"ee_ampa": 0.5},
                    weights_provenance=f"tuned against {P.name}")
    assert good.weights["ee_ampa"] == 0.5

    with pytest.raises(ProvenanceMismatch, match="does not name cell"):
        ParamSet(cell=P, weights={"ee_ampa": 0.55},
                 weights_provenance="hand-tuned against the LIF Pop, session 4")

    with pytest.raises(ProvenanceMismatch, match="no provenance"):
        ParamSet(cell=P, weights={"ee_ampa": 0.55})


def test_E14_a_cell_with_no_weights_needs_no_provenance():
    """Link 4 ships the cell alone; the weight table arrives with link 5. An empty table must
    not be blocked by the guard that protects a populated one."""
    assert ParamSet(cell=P).weights == {}


def test_E14_registry_entry_carries_its_source_and_preparation():
    """A parameter set without its preparation is not citable, and the preparation is
    load-bearing here: these are NEONATAL (~P0-P4) values, which is why adopting them moves
    the >P12 muscimol anchor further away rather than closer."""
    p = CELLS["butera_rinzel_smith_1999_model1"]
    for token in ("Butera", "1999", "J Neurophysiol", "P0-P4"):
        assert token in p.provenance, f"provenance does not record {token!r}"
    assert "not from recall" in p.provenance


def test_E14_the_one_parameter_that_was_misremembered_is_right():
    """E_L is the bifurcation parameter of this model -- the knob that moves the cell
    quiescent -> bursting -> tonic. The recalled value was -65 mV; the published value read
    off a machine-readable encoding is -57.5 mV. Starting from -65 would have put the cell in
    the wrong regime, and the rhythm would then have been retuned around a wrong resting
    drive. Pinned because it is the single value most expensive to get wrong."""
    assert P.e_l == -57.5
    assert P.C == 21.0 and P.g_na == 28.0 and P.g_k == 11.2
    assert P.g_nap == 2.8 and P.g_l == 2.8
    assert P.e_na == 50.0 and P.e_k == -85.0
    assert P.m == (-34.0, -5.0) and P.n == (-29.0, -4.0)
    assert P.mp == (-40.0, -6.0) and P.h == (-48.0, 6.0)
    assert P.tau_n == 10.0 and P.tau_h == 10000.0

    # The recalled value does not merely shift the regime: it gives a SILENT cell. And the
    # bursting band is only a few mV wide, with tonic spiking on both sides -- which is what
    # makes E_L the single parameter here most expensive to get wrong, and why reading it off
    # a machine-readable source rather than reciting it was worth the round trip.
    #
    #   E_L   -65.0  -60.0  -57.5  -55.0  -52.5
    #         silent tonic  BURST  tonic  tonic
    wrong = run_isolated(seconds=20.0, params=replace(P, e_l=-65.0))
    assert wrong["regime"] == "quiescent" and wrong["n_spikes"] == 0, (
        f"the misremembered E_L gives {wrong['regime']} at {wrong['rate_hz']:.1f} Hz rather "
        f"than silence; recheck that this pin is testing what it claims")
    for el in (-60.0, -55.0):
        near = run_isolated(seconds=20.0, params=replace(P, e_l=el))
        assert near["regime"] != "bursting", (
            f"E_L={el} also bursts, so the bursting window is wider than measured when this "
            f"cell was adopted and this pin is weaker than it reads")


# ====================================== E15: Mg2+ relief is evaluated at the cell's own V
def test_E15_nmda_block_relieves_during_a_spike():
    """The main NMDA benefit of this upgrade, and the way it would most plausibly be thrown
    away silently.

    The LIF interface hands over a conductance already evaluated at the PRE-STEP voltage.
    That is harmless when V moves a few mV per step, and it discards the entire
    voltage-dependent relief when V sweeps ~70 mV through a spike. The symptom would be "the
    upgrade did not change the NMDA result", which reads as a reassuring robustness check and
    is in fact the bug.

    Measured on this cell: relief ranges 0.071 .. 0.842 over the trace, an 11.8x span against
    the 4.50x measured on the LIF preBotC (scripts/diag_substrate_limits.py).
    """
    r = run_isolated(seconds=20.0)
    rel = _mg_relief(r["V"])
    span = float(rel.max() / rel.min())
    assert span > 8.0, (
        f"Mg2+ relief spans only {span:.2f}x; the LIF substrate already reached 4.50x, so "
        f"this would not be an improvement worth the cost")
    peak_over_rest = float(rel.max() / _mg_relief(np.percentile(r["V"], 5)))
    assert peak_over_rest > 5.0, f"peak/rest relief ratio only {peak_over_rest:.2f}x"


def test_E15_raw_nmda_conductance_is_voltage_scaled_inside_the_cell():
    """`raw=("nmda",)` must make the cell apply the block itself. Checked by the effect on
    the membrane: an identical raw conductance delivers MORE current when the cell is
    depolarised, which is the state-dependence the LIF cannot express."""
    g, E = {"nmda": np.array([5.0])}, {"nmda": 0.0}
    frac = {}
    for V0 in (-65.0, -20.0):
        c = isolated()
        c.V[:] = V0
        intrinsic = c._currents(c.V, {}, {}, ())
        i_without = c._currents(c.V, g, E, raw=()) - intrinsic
        i_with = c._currents(c.V, g, E, raw=("nmda",)) - intrinsic
        # the raw path must be a FRACTION of the unscaled one (the block is attenuation)
        frac[V0] = float((i_with / i_without)[0])
    assert frac[-65.0] < frac[-20.0], (
        f"NMDA attenuation did not weaken with depolarisation ({frac[-65.0]:.4f} at -65 mV "
        f"vs {frac[-20.0]:.4f} at -20 mV); the block is not evaluated at the cell's own V")
    assert frac[-65.0] == pytest.approx(float(_mg_relief(-65.0)), rel=1e-9)
    assert frac[-20.0] == pytest.approx(float(_mg_relief(-20.0)), rel=1e-9)


# ================================================================= E16: spike detection
def test_E16_spike_count_is_invariant_to_the_detection_threshold():
    """A bare upward crossing counts one broad spike several times as V jitters across the
    line, inflating every firing rate -- and an inflated rate looks like a plausible number.
    With a re-arm lockout, the count must not depend on where in the spike's rising phase the
    detector sits."""
    counts = {}
    for vd in (-30.0, -20.0, -10.0):
        r = run_isolated(seconds=10.0,
                         params=replace(P, v_detect=vd, v_reset_lock=vd - 10.0))
        counts[vd] = int(r["spikes_ms"].size)
    lo, hi = min(counts.values()), max(counts.values())
    assert hi - lo <= 2, f"spike count varies with the detection threshold: {counts}"


def test_E16_a_cell_held_depolarised_does_not_spike_repeatedly():
    """The lockout's job, in isolation: with V clamped above the detection voltage and never
    returning below the re-arm level, exactly one detection may occur."""
    c = isolated()
    rng = np.random.default_rng(0)
    n = 0
    for _ in range(500):
        c.V[:] = 0.0                       # above v_detect, never below v_reset_lock
        n += int(c.step(0.1, {}, {}, 0.0, rng)[0])
    assert n == 1, f"detector fired {n} times while V was held depolarised; expected 1"


# ============================================================ A7: the LIF path is intact
def test_A7_lif_population_trace_is_unchanged():
    """Everything above is additive: `cpg.Pop` was not touched. This pins its trace so that a
    later refactor cannot alter the LIF substrate while the comparison between substrates is
    the deliverable -- a silently drifting baseline would invalidate it.

    The hash is a going-forward regression guard, recorded when the conductance backend was
    added; it is not evidence about the past.
    """
    p = Pop(n=8, name="A7probe", nap=True)
    rng = np.random.default_rng(12345)
    V = []
    for _ in range(2000):
        p.step(0.1, {}, {}, 180.0, rng)
        V.append(p.V.copy())
    digest = hashlib.sha256(np.round(np.asarray(V), 6).tobytes()).hexdigest()[:16]
    assert digest == "66cc15b79cc4c4e7", (
        f"LIF trace hash changed to {digest}; cpg.Pop must stay byte-identical while the "
        f"substrate comparison depends on it")


# ================================================================= E18 / interface shims
def test_E18_adaptation_is_off_by_default_but_still_available():
    """Spike-triggered adaptation was a stand-in for I_NaP inactivation. Now that h is real,
    leaving it at its LIF-tuned strength double-counts burst termination; the symptom would be
    bursts terminating slightly early, i.e. a plausible duty cycle. Off by default, and the
    machinery is kept so the double-counting can be measured rather than argued about."""
    assert CondPop(n=1, name="d").g_adapt == 0.0
    c = CondPop(n=1, name="e", g_adapt=1.0, sigma_i=0.0)
    rng = np.random.default_rng(0)
    for _ in range(3000):
        c.step(0.1, {}, {}, 0.0, rng)
    assert c.a.max() > 0.0, "adaptation did not accumulate when explicitly enabled"


def test_interface_matches_the_LIF_population():
    """A circuit must be able to swap substrates without touching its stepping loop; that is
    what makes the LIF-vs-conductance comparison practical."""
    c, p = CondPop(n=4, name="iface"), Pop(n=4, name="iface")
    for attr in ("V", "spk", "rate", "a", "g_adapt", "tau_adapt", "EK", "EL", "C", "Vth"):
        assert hasattr(c, attr), f"CondPop is missing {attr!r}, which cpg.Pop exposes"
    rng = np.random.default_rng(0)
    out = c.step(0.1, {"ampa": np.zeros(4)}, {"ampa": 0.0}, 0.0, rng)
    assert out.shape == p.spk.shape and out.dtype == bool


def test_seeding_is_stable_across_processes():
    """Same lesson as cpg.Pop: `hash(str)` is randomised per process via PYTHONHASHSEED, so
    initial voltages differed between runs even under an explicit seed, making every result
    irreproducible. crc32 is stable."""
    a, b = CondPop(n=6, name="Exc0"), CondPop(n=6, name="Exc0")
    assert np.array_equal(a.V, b.V)
    assert not np.array_equal(a.V, CondPop(n=6, name="Inh0").V)


# ========================= the substrate has an UPPER bound on drive; the LIF has none
def test_conductance_cell_shows_depolarisation_block_and_the_LIF_does_not():
    """A substrate difference with no LIF analogue, and the one that broke the network.

    Measured on an isolated cell with g_NaP = 0 (as the Inh and Out populations are built):
    39 Hz at 40 pA rising to 136 Hz at 150 pA, then SILENT at 300 pA. At that potential the
    Na inactivation term (1-n) is nearly zero, so no spike can be produced however much
    current is injected.

    The LIF cannot do this. More drive means a higher rate, capped only by tref. So any
    reasoning of the form "output collapsed, therefore increase the drive" is valid on the
    LIF substrate and WRONG on this one -- which is exactly the mistake made when first
    deriving a conductance operating point, where `drive_other` was raised to 40-70 pA on a
    threshold calculation while the Out population was already in block at 25 pA plus
    synaptic input.
    """
    no_nap = replace(P, g_nap=0.0)
    firing = run_isolated(seconds=6.0, params=no_nap, i_app=150.0)
    blocked = run_isolated(seconds=6.0, params=no_nap, i_app=300.0)
    assert firing["rate_hz"] > 50.0, (
        f"the no-I_NaP cell should fire briskly at 150 pA, got {firing['rate_hz']:.1f} Hz")
    assert blocked["rate_hz"] == 0.0, (
        f"expected depolarisation block at 300 pA, got {blocked['rate_hz']:.1f} Hz. If this "
        f"cell no longer blocks, the Out-population diagnosis needs revisiting.")

    # The LIF, driven far harder, keeps firing: it has no inactivation to lose.
    lif = Pop(n=4, name="blockprobe")
    rng = np.random.default_rng(0)
    n_spk = 0
    for _ in range(20000):
        n_spk += int(Pop.step(lif, 0.1, {}, {}, 2000.0, rng).sum())
    assert n_spk > 0, (
        "the LIF stopped firing at very high drive, so it has acquired a block mechanism "
        "and no longer provides the contrast this test documents")
