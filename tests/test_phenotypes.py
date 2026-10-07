"""Published phenotypes the model reproduces WITHOUT having been fitted to them.

This file is the only place in the suite that tests biology rather than code. Everything
else asserts identities or guards failure modes; these assert that the circuits reproduce
experimental results they were not tuned against, which is the one thing that can promote a
quantity to VALIDATED.

There are few of them, and that is the honest state of the model.
"""
import numpy as np
import pytest

from circuitpharm.rg2 import GroupPacemakerRG
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.cpg import Drug, burst_metrics
from circuitpharm.config import RESP_OP

T, DT, WARM = 20000.0, 0.1, 6000.0


def _rg(ie_gly=None, drug=None, seed=1):
    w = None if ie_gly is None else dict(ie_gly=ie_gly)
    r = GroupPacemakerRG(drug=drug, w=w, seed=seed)
    for i in range(int(T / DT)):
        r.step(DT)
        if i % 10 == 0:
            r.record()
    A = r.arrays(); m = A["t"] > WARM
    k = np.ones(50) / 50.0
    sm = lambda v: np.convolve(v, k, mode="same")
    f, e = sm(A["RG_F"][m]), sm(A["RG_E"][m])
    return dict(f=f, e=e, t=A["t"][m],
                corr=float(np.corrcoef(f, e)[0, 1]) if f.std() > 1e-9 and e.std() > 1e-9
                else float("nan"),
                nb_F=burst_metrics(A["t"][m], f)["n_bursts"],
                nb_E=burst_metrics(A["t"][m], e)["n_bursts"],
                per_F=burst_metrics(A["t"][m], f)["period"],
                per_E=burst_metrics(A["t"][m], e)["period"])


# ============================================================ locomotor RG: strychnine
# PUBLISHED TARGET (lamprey fictive locomotion): strychnine, i.e. glycine receptor block,
# ELIMINATES left-right alternation while robust rhythmic activity PERSISTS. Phase and
# rhythm are therefore separable, and glycinergic inhibition sets phase only.
#
# This is the entire reason the group-pacemaker architecture replaced a Matsuoka
# oscillator: in Matsuoka, mutual inhibition IS the oscillator, so removing it stops the
# rhythm, contradicting the data. An earlier prediction of this project rested on the
# Matsuoka pathway and was retracted.
@pytest.mark.slow
def test_module_defaults_alternate_at_all():
    """Guards a real shipped bug: the module's defaults were once the untuned initial
    guesses while the validated set lived only in a past run's stdout. One half-centre was
    silent entirely, so anything importing the module inherited a broken rhythm."""
    r = _rg()
    assert r["f"].max() > 1.0, "flexor half-centre silent"
    assert r["e"].max() > 1.0, "extensor half-centre silent at module defaults"
    assert r["corr"] < -0.3, f"not alternating at defaults (corr {r['corr']:+.2f})"


@pytest.mark.slow
def test_strychnine_phenotype_glycine_removal_abolishes_alternation():
    """Glycine removed COMPLETELY -> alternation lost, rhythm intact."""
    coupled = _rg()
    removed = _rg(ie_gly=0.0)

    # alternation must weaken substantially
    assert removed["corr"] - coupled["corr"] > 0.25, (
        f"removing glycine barely changed alternation "
        f"({coupled['corr']:+.2f} -> {removed['corr']:+.2f}); the coupling is not "
        "load-bearing for phase")
    # ...and the rhythm must SURVIVE. This is the half that distinguishes a group
    # pacemaker from a Matsuoka oscillator.
    assert removed["nb_F"] >= 8 and removed["nb_E"] >= 8, (
        "removing glycine abolished the rhythm, which is the Matsuoka failure mode this "
        "architecture exists to avoid")


@pytest.mark.slow
def test_E10_partial_block_cannot_detect_the_coupling():
    """Recurring error E10, as an executable demonstration.

    Entrainment is cheap: a 95% conductance reduction leaves alternation fully intact, so a
    partial block is blind to whether a coupling matters. This test asserts the error is
    real, which is what justifies testing couplings only by complete removal -- including
    in scripts/tune_rg2.py, which still uses the 95% form and is therefore uninformative
    about this phenotype.
    """
    coupled = _rg()
    blocked_95 = _rg(drug=Drug(glyr_gain=0.05))
    removed = _rg(ie_gly=0.0)

    assert abs(blocked_95["corr"] - coupled["corr"]) < 0.15, (
        "95% glycine block changed alternation more than expected; if partial block has "
        "become informative, E10 no longer applies here and the guidance should change")
    assert abs(removed["corr"] - coupled["corr"]) > 2 * abs(
        blocked_95["corr"] - coupled["corr"]), (
        "complete removal was not clearly more informative than 95% block")


# ================================================ preBotC: inhibition vs adaptation roles
@pytest.mark.slow
def test_adaptation_not_inhibition_terminates_prebotc_bursts():
    """Removing adaptation must abolish the rhythm (tonic firing); removing synaptic
    inhibition must NOT. Establishes which mechanism is rhythmogenic, and it is the result
    that made the NMDA-disinhibition negative finding interpretable."""
    def run(**kw):
        b = PreBotC(seed=0, **{**RESP_OP, **kw})
        for i in range(int(14000 / DT)):
            b.step(DT)
            if i % 10 == 0:
                b.record()
        A = b.arrays(); m = A["t"] > 4000
        return resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])

    intact = run()
    no_adapt = run(g_adapt=0.0)
    no_inhib = run(w=dict(ie_gaba=0.0, ie_gly=0.0, **dict(RESP_OP["w"])))

    assert intact["alive"], intact["reason"]
    assert not no_adapt["alive"], (
        "removing spike-triggered adaptation left a rhythm; adaptation is supposed to be "
        "what terminates the burst")
    assert no_inhib["alive"], (
        "removing synaptic inhibition abolished the rhythm; in a group pacemaker "
        "inhibition should shape output, not generate the rhythm")
    # inhibition is load-bearing for AMPLITUDE even though it is not rhythmogenic --
    # this is what made the 'disinhibition cannot rescue NMDA block' result non-trivial
    assert no_inhib["mean"] > 1.2 * intact["mean"], (
        "removing inhibition did not raise output; then there would be no inhibition for "
        "an NMDA antagonist to release, and that negative result would be vacuous")


# ============================================================ spinal reflex phenotypes
# These are reproduced by a circuit that was NOT fitted to them: the reflex arc's
# connectivity was hand-built from anatomy, and the drug effects follow from receptor
# placement. Strychnine hyperreflexia in particular is a textbook phenotype.
@pytest.mark.slow
def test_strychnine_causes_hyperreflexia():
    """Glycine block must INCREASE reflex gain. Glycinergic Renshaw and reciprocal Ia
    inhibition restrain the motoneuron pool, so removing it releases the reflex -- the
    classic strychnine phenotype."""
    from circuitpharm.assays import stretch_reflex
    from circuitpharm.cpg import Drug
    ctrl = stretch_reflex(Drug())["gain"]
    stry = stretch_reflex(Drug(glyr_gain=0.4))["gain"]
    assert stry > 1.2 * ctrl, (
        f"glycine block gave gain {stry:.3f} vs control {ctrl:.3f}; strychnine must "
        "produce hyperreflexia")


@pytest.mark.slow
def test_benzodiazepine_depresses_the_reflex():
    """Diazepam depresses the stretch reflex; the model must too."""
    from circuitpharm.assays import stretch_reflex
    from circuitpharm.cpg import Drug
    ctrl = stretch_reflex(Drug())["gain"]
    bz = stretch_reflex(Drug(gaba_a_gain=2.0, gaba_a_tau=1.6))["gain"]
    assert bz < 0.9 * ctrl, f"benzodiazepine gave gain {bz:.3f} vs control {ctrl:.3f}"


@pytest.mark.slow
def test_glun2b_selectivity_window_spares_the_spinal_reflex():
    """THE MECHANISM THE WHOLE NMDA ARM RESTS ON. At the same 60% block, a GluN2B-selective
    antagonist must spare the spinal reflex (low local GluN2B fraction) where a
    non-selective one depresses it. If this collapses, there is no selectivity window."""
    from circuitpharm.assays import stretch_reflex
    from circuitpharm.cpg import Drug
    ctrl = stretch_reflex(Drug())["gain"]
    nonsel = stretch_reflex(Drug(nmda_block=0.6, glun2b_selectivity=0.0))["gain"]
    sel = stretch_reflex(Drug(nmda_block=0.6, glun2b_selectivity=1.0,
                              glun2b_fraction=0.15))["gain"]
    assert nonsel < 0.92 * ctrl, "non-selective NMDA block should depress the reflex"
    assert sel > 0.95 * ctrl, (
        f"GluN2B-selective block depressed the reflex to {100*sel/ctrl:.0f}% of control; "
        "the selectivity window has collapsed")
    assert sel > nonsel


@pytest.mark.slow
def test_E5_motoneuron_peak_is_not_pinned_at_the_firing_ceiling():
    """Recurring error E5: with the full Ia->Mn weight the motoneuron dynamic peak sits at
    the tref=8 ms cap (125 Hz) and EVERY drug reads ~100% of control. The assay rescales
    the weight so controls sit mid-range; this asserts it worked."""
    from circuitpharm.assays import stretch_reflex
    from circuitpharm.cpg import Drug
    r = stretch_reflex(Drug())
    cap_hz = 1000.0 / 8.0        # tref = 8 ms
    assert r["mn_dyn"] < 0.8 * cap_hz, (
        f"motoneuron dynamic peak {r['mn_dyn']:.1f} Hz is at the {cap_hz:.0f} Hz "
        "structural cap; drug effects will be masked")
