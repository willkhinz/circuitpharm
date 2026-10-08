"""The public evaluation API, including that it refuses to hand out VOID numbers."""
import numpy as np
import pytest

from circuitpharm.evaluate import Compound, evaluate
from circuitpharm.results import Tier, VoidQuantityError
from circuitpharm.subtypes import PROFILES


# ------------------------------------------------------------------ construction
def test_from_profile_matches_the_named_profile():
    c = Compound.from_profile("alogabat")
    p = PROFILES["alogabat"]
    for s in ("a1", "a23", "a5", "d_a4", "eps"):
        assert getattr(c, s) == getattr(p, s)


def test_zero_occupancy_is_drug_free():
    """Occupancy 0 must give identity gains, or every control is contaminated."""
    g = Compound("x", a5=1.0, occupancy=0.0).pool_gains()
    assert g["tonic"] == pytest.approx(1.0)
    assert g["phasic"] == pytest.approx(1.0)
    assert g["tau"] == pytest.approx(1.0)


def test_pool_gains_are_linear_in_occupancy():
    """A partially occupied receptor population is a MIXTURE of modulated and unmodulated
    receptors, so each pool's conductance gain is exactly linear in occupancy. That
    linearity is why the gains saturate at full occupancy, and that saturation is the only
    honest ceiling."""
    c = Compound("x", a5=1.0, s_max=2.5)
    c.occupancy = 1.0
    full = c.pool_gains()
    for occ in (0.25, 0.5, 0.75):
        c.occupancy = occ
        g = c.pool_gains()
        assert g["tonic"] == pytest.approx(1.0 + occ * (full["tonic"] - 1.0), rel=1e-9)
        assert g["phasic"] == pytest.approx(1.0 + occ * (full["phasic"] - 1.0), rel=1e-9)


def test_tonic_gain_greatly_exceeds_phasic():
    """The central kinetic result: an affinity-type PAM barely moves a near-saturated
    synaptic peak while strongly potentiating sub-saturating tonic current."""
    g = Compound("x", a5=1.0, occupancy=1.0, s_max=2.5).pool_gains()
    assert g["tonic"] > 3.0 * g["phasic"]


def test_gains_saturate_at_full_occupancy():
    """Beyond full occupancy there is no mechanism left -- this is the real ceiling."""
    c = Compound("x", a5=1.0, s_max=2.5, occupancy=1.0)
    g = c.pool_gains()
    assert g["tonic"] == pytest.approx(g["tonic_at_full"], rel=1e-9)


def test_higher_intrinsic_efficacy_gives_a_higher_ceiling():
    """`s_max` is a LIGAND property and is what actually sets overdose protection."""
    lo = Compound("lo", a5=1.0, s_max=2.0, occupancy=1.0).pool_gains()["tonic"]
    hi = Compound("hi", a5=1.0, s_max=5.0, occupancy=1.0).pool_gains()["tonic"]
    assert hi > lo


# ------------------------------------------------- the calibration-independent ranking
def test_selectivity_ratio_of_reference_is_one():
    assert Compound.from_profile("nonselective_bz").selectivity_ratio() == \
        pytest.approx(1.0)


def test_a5_arms_beat_nonselective_and_neurosteroid_loses():
    """Reproduces the 20,000-draw result's central tendency. The neurosteroid arm coming
    out WORSE than a plain benzodiazepine independently reproduces an earlier retraction."""
    a5 = Compound.from_profile("alogabat").selectivity_ratio()
    ns = Compound.from_profile("neurosteroid").selectivity_ratio()
    a23 = Compound.from_profile("hz_166").selectivity_ratio()
    assert a5 > 3.0, f"a5 arm ratio {a5:.2f} unexpectedly low"
    assert a5 > a23 > ns
    assert ns < 1.0, "neurosteroid should be worse than a non-selective benzodiazepine"


def test_zolpidem_is_bad_because_a1_carries_no_subjective_weight():
    """a1-preferring: all respiratory burden, no subjective drive. A sanity check that the
    subtype weighting is wired the right way round."""
    assert Compound.from_profile("zolpidem").selectivity_ratio() < 0.5


def test_selectivity_ratio_is_scale_invariant():
    """The whole point of the ratio: it must not change when the unknown lumped
    sensitivity factor changes. Scaling every efficacy scales numerator and denominator
    together, so the ratio is fixed."""
    a = Compound.from_profile("alogabat")
    b = Compound.from_profile("alogabat")
    b.a1, b.a23, b.a5 = 2 * a.a1, 2 * a.a23, 2 * a.a5
    assert b.selectivity_ratio() == pytest.approx(a.selectivity_ratio(), rel=1e-9)


def test_enantiomer_pair_is_separated():
    """The project's central chiral pair: one stereocentre switches which GABA arm the
    molecule serves, and the ranking must see it."""
    r = Compound.from_profile("sh053_R").selectivity_ratio()
    s = Compound.from_profile("sh053_S").selectivity_ratio()
    assert r > 2.0 * s


# ------------------------------------------------------------------- tiered results
@pytest.mark.slow
def test_evaluate_returns_tiered_results_and_withholds_void():
    rs = evaluate(Compound.from_profile("alogabat", occupancy=0.35), n_seed=1)

    assert rs.quantity("selectivity_ratio").tier is Tier.VALIDATED
    assert rs.quantity("ventilation").tier is Tier.UNCALIBRATED
    assert rs.quantity("overdose_index").tier is Tier.VOID

    # the refusal must actually fire
    with pytest.raises(VoidQuantityError):
        rs.quantity("overdose_index").value
    # ...and the printed form must not leak it
    assert "overdose_index" in str(rs)
    assert "withheld" in str(rs)


@pytest.mark.slow
def test_nmda_arm_adds_its_own_void_quantity():
    """With an NMDA component the respiratory contribution is a known negative result and
    must be withheld rather than reported."""
    rs = evaluate(Compound.from_profile("alogabat", occupancy=0.35, nmda_block=0.3),
                  n_seed=1)
    q = rs.quantity("nmda_respiratory_contribution")
    assert q.tier is Tier.VOID
    assert "preserve" in q.provenance.lower() or "stimulate" in q.provenance.lower()


@pytest.mark.slow
def test_every_uncalibrated_quantity_states_how_to_promote_it_or_why_not():
    rs = evaluate(Compound.from_profile("alogabat", occupancy=0.35), n_seed=1)
    for q in rs.by_tier(Tier.UNCALIBRATED):
        assert q.provenance, f"{q.name} has no provenance"
    # the two headline ones must name their promotion path explicitly
    assert "P12" in rs.quantity("ventilation").promote_by
    assert "learns" in rs.quantity("subjective_index_gaba").promote_by


# ------------------------------------------------------------------ motor endpoints
@pytest.mark.slow
@pytest.mark.needs_plant
def test_evaluate_includes_all_three_endpoints():
    """The package has respiratory, reflex and locomotor endpoints; the public API must
    expose all three, or users will reach past it into the modules."""
    rs = evaluate(Compound.from_profile("alogabat", occupancy=0.35), n_seed=1)
    for name in ("ventilation", "reflex_gain", "step_period", "coordination"):
        assert rs.quantity(name) is not None, name


@pytest.mark.slow
@pytest.mark.needs_plant
def test_joint_excursion_is_void_in_the_api_too():
    """The invalid metric must be VOID at the API boundary, not just renamed in the assay.
    Renaming protects a careful reader; the tier protects everyone else."""
    rs = evaluate(Compound.from_profile("alogabat", occupancy=0.35), n_seed=1)
    q = rs.quantity("joint_excursion")
    assert q.tier is Tier.VOID
    with pytest.raises(VoidQuantityError):
        q.value
    assert "WRONG SIGN" in q.provenance


@pytest.mark.slow
def test_motor_endpoints_can_be_skipped():
    """They need the optional body plant, so evaluation must work without them."""
    rs = evaluate(Compound.from_profile("alogabat", occupancy=0.35), n_seed=1,
                  include_motor=False)
    assert rs.quantity("ventilation") is not None
    with pytest.raises(KeyError):
        rs.quantity("reflex_gain")


def test_drugfree_control_is_independent_of_sensitivity():
    """PINS THE ASSUMPTION THE CONTROL CACHE RESTS ON.

    Controls are cached on seed alone, which is only valid because a drug-free Drug()
    gives eff = 1 + sens*(x-1) with x == 1, i.e. 1 for ANY sensitivity. If a future change
    makes a drug-free control depend on sensitivity, the cache silently returns the wrong
    normaliser and every percentage shifts. This test fails first.
    """
    from circuitpharm.cpg import Drug, Syn, TAU
    d = Drug()
    for sens in (0.0, 0.05, 0.5, 1.0):
        for kind in ("gabaa", "gly", "nmda", "ampa"):
            syn = Syn(4, kind, d, sens=sens)
            assert syn.w_scale == pytest.approx(1.0), (kind, sens)
            assert syn.tau == pytest.approx(TAU[kind]), (kind, sens)


def test_control_cache_can_be_cleared():
    from circuitpharm.evaluate import _CTRL_CACHE, clear_control_cache
    _CTRL_CACHE[("probe", 0)] = "x"
    clear_control_cache()
    assert not _CTRL_CACHE


def test_evaluate_degrades_gracefully_without_the_body_plant():
    """A lean install must still get the respiratory and receptor-level results.

    REGRESSION TEST for a bug that passed locally and would have failed CI. The guard in
    evaluate() originally wrapped only the IMPORT of the assay functions -- but `assays`
    imports mujoco lazily INSIDE each function, so on a lean install that import succeeds
    and the ImportError surfaces on the first CALL, outside the guard. evaluate() therefore
    crashed on exactly the install CI uses (`.[dev]`, no plant extra) while passing on this
    machine, where mujoco happens to be present.

    Simulated here by hiding the module, so the test is meaningful even WITH mujoco
    installed -- otherwise it would silently pass for the wrong reason on a full install.
    """
    import builtins
    import sys

    real_import = builtins.__import__

    def fake(name, *a, **k):
        if name.split(".")[0] in ("mujoco", "dm_control"):
            raise ImportError(f"No module named {name!r} (simulated lean install)")
        return real_import(name, *a, **k)

    saved = {m: sys.modules[m] for m in list(sys.modules)
             if m.split(".")[0] in ("mujoco", "dm_control", "circuitpharm")}
    try:
        builtins.__import__ = fake
        for m in saved:
            del sys.modules[m]
        from circuitpharm.evaluate import Compound as C, evaluate as ev
        rs = ev(C.from_profile("alogabat", occupancy=0.35), n_seed=1)
        names = [q.name for q in rs.items]
        assert "ventilation" in names, "respiratory endpoint lost on a lean install"
        assert "selectivity_ratio" in names
        assert "motor_endpoints" in names, "no notice that motor endpoints are unavailable"
        assert not {"reflex_gain", "step_period"} & set(names)
        assert "plant" in rs.quantity("motor_endpoints").promote_by
    finally:
        builtins.__import__ = real_import
        for m in list(sys.modules):
            if m.split(".")[0] in ("mujoco", "dm_control", "circuitpharm"):
                del sys.modules[m]
        sys.modules.update(saved)
