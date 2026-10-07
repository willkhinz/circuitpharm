"""Exact identities the model must satisfy. These are not approximate comparisons.

An identity test is the strongest kind available to a model with no validated calibration:
it does not ask whether the model is right about biology, it asks whether the code computes
what it claims to. Every assertion here should hold to numerical precision, so any
tolerance is either zero or set by float round-off -- never by "close enough".
"""
import numpy as np
import pytest

from circuitpharm.cpg import Drug
from circuitpharm.subtypes import PROFILES, REGIONS, SUBTYPES, EXTRASYN


# --------------------------------------------------------------- drug-free identity
def test_drugfree_drug_is_identity():
    """A default Drug must perturb nothing. If this breaks, every control is wrong."""
    d = Drug()
    assert d.gaba_scale() == 1.0
    assert d.gaba_scale_tonic() == 1.0
    assert d.nmda_scale() == 1.0
    assert d.gaba_a_tau == 1.0
    assert d.glyr_gain == 1.0


def test_tonic_gain_defaults_to_phasic():
    """`gaba_a_gain_tonic=None` must reproduce the legacy single-pool behaviour EXACTLY.

    This backward-compatibility guarantee is what allows pre-split results to remain
    reproducible after the tonic/phasic split was introduced.
    """
    d = Drug(gaba_a_gain=2.0, gaba_a_tau=1.6)
    assert d.gaba_scale() == 2.0
    assert d.gaba_scale_tonic() == 2.0


def test_tonic_and_phasic_are_independently_settable():
    d = Drug(gaba_a_gain=1.06, gaba_a_gain_tonic=7.21,
             gaba_a_efficacy_cap=1e9, gaba_a_cap_tonic=1e9)
    assert d.gaba_scale() == pytest.approx(1.06)
    assert d.gaba_scale_tonic() == pytest.approx(7.21)


def test_efficacy_caps_bind_per_pool():
    d = Drug(gaba_a_gain=10.0, gaba_a_gain_tonic=100.0,
             gaba_a_efficacy_cap=2.5, gaba_a_cap_tonic=50.0)
    assert d.gaba_scale() == 2.5
    assert d.gaba_scale_tonic() == 50.0


# ------------------------------------------------------- NMDA selectivity window
def test_glun2b_selectivity_bounds_the_block():
    """A fully GluN2B-selective antagonist cannot block more than the GluN2B pool.

    This bound IS the selectivity window the whole NMDA arm rests on: brainstem is
    GluN2D-dominant (low GluN2B fraction) so a selective antagonist should barely touch it.
    """
    full = Drug(nmda_block=1.0, glun2b_selectivity=1.0, glun2b_fraction=0.15)
    assert full.nmda_scale() == pytest.approx(1.0 - 0.15)

    nonsel = Drug(nmda_block=1.0, glun2b_selectivity=0.0, glun2b_fraction=0.15)
    assert nonsel.nmda_scale() == pytest.approx(0.0)


def test_nmda_block_zero_is_identity():
    for sel in (0.0, 0.5, 1.0):
        assert Drug(nmda_block=0.0, glun2b_selectivity=sel).nmda_scale() == 1.0


# ------------------------------------------------- tonic/phasic split is a partition
@pytest.mark.parametrize("key", sorted(PROFILES))
def test_regional_split_sums_to_total(key):
    """tonic + phasic must equal the unsplit total, for every profile and region.

    The split is a PARTITION of the same conductance by localisation. If it stopped summing
    the per-region calibration would silently change meaning.
    """
    p = PROFILES[key]
    for region in REGIONS:
        t, ph = p.regional_sens_split(region)
        assert t + ph == pytest.approx(p.regional_sens(region), abs=1e-12)


@pytest.mark.parametrize("key", sorted(PROFILES))
def test_subjective_split_sums_to_total(key):
    p = PROFILES[key]
    t, ph = p.subjective_index_split()
    assert t + ph == pytest.approx(p.subjective_index(), abs=1e-12)


def test_split_components_are_non_negative():
    for key, p in PROFILES.items():
        for region in REGIONS:
            t, ph = p.regional_sens_split(region)
            assert t >= -1e-15 and ph >= -1e-15, key


def test_extrasyn_fractions_are_fractions():
    for s in SUBTYPES:
        assert 0.0 <= EXTRASYN[s] <= 1.0, s
    # delta/alpha4 receptors are exclusively extrasynaptic -- this one is solid biology,
    # not an estimate, so it is pinned.
    assert EXTRASYN["d_a4"] == 1.0


def test_regional_fractions_sum_to_one():
    for region, f in REGIONS.items():
        assert sum(f[s] for s in SUBTYPES) == pytest.approx(1.0, abs=1e-9), region


# ------------------------------------------------------- per-region calibration anchors
def test_nonselective_bz_reproduces_each_regional_anchor():
    """The per-region constants are DEFINED so a non-selective BZ reproduces each anchor.

    Error E12 was committed by applying one region's constant to another. This asserts the
    construction still holds, which is what makes selectivity move sensitivity relative to
    an anchor rather than by refitting.
    """
    from circuitpharm.subtypes import ANCHOR
    ns = PROFILES["nonselective_bz"]
    for region, anchor in ANCHOR.items():
        assert ns.regional_sens(region) == pytest.approx(anchor, rel=1e-9), region


# ----------------------------------------------------------------- alpha5 localisation
def test_a5_is_predominantly_extrasynaptic():
    """Load-bearing for the central result: a5 carries ethanol's discriminative stimulus
    AND sits mostly extrasynaptically, where a PAM has ~200x headroom rather than ~1.1x.
    That combination is why the efficacy-ceiling safety argument fails exactly where the
    a5 strategy places the drug."""
    assert EXTRASYN["a5"] > 0.5
    assert EXTRASYN["a5"] > EXTRASYN["a1"]
