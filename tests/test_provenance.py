"""Every load-bearing number must declare where it came from.

`results.Tier` makes a model OUTPUT carry its reliability. These tests do the same for model
INPUTS, and they exist because an audit found the harder problem: the numbers driving the
project's one VALIDATED result were not unverified, they were UNSOURCED, and eight review
passes had not caught it because reviews read code rather than provenance.

The point is not to make unsourced numbers look sourced. It is to make the count visible and
to stop a new number being added without a record.
"""
import pathlib
import sqlite3

import pytest

from circuitpharm.provenance import (
    ALL, Basis, EXTRASYN_PROV, PRIOR_PROV, REGIONS_PROV, SUBJECTIVE_PROV, audit, report,
)
from circuitpharm.subtypes import EXTRASYN, REGIONS, SUBJECTIVE_WEIGHT, SUBTYPES

DB = pathlib.Path(__file__).resolve().parent.parent / "data" / "pharmacology.db"


# ======================================================= completeness: no silent additions
def test_every_regional_subunit_fraction_has_provenance():
    missing = [(r, s) for r in REGIONS for s in SUBTYPES if (r, s) not in REGIONS_PROV]
    assert not missing, (
        f"REGIONS entries with no provenance record: {missing}. These numbers are the "
        f"direct input to the selectivity ranking; a new one must not be addable without "
        f"saying where it came from.")


def test_every_extrasynaptic_fraction_has_provenance():
    assert not [s for s in EXTRASYN if s not in EXTRASYN_PROV]


def test_every_subjective_weight_has_provenance():
    assert not [s for s in SUBJECTIVE_WEIGHT if s not in SUBJECTIVE_PROV]


def test_provenance_does_not_describe_parameters_that_no_longer_exist():
    """The reverse direction. A record for a deleted parameter is stale documentation, which
    is how this project's worklog drifted from its code before."""
    stale = [k for k in REGIONS_PROV if k[0] not in REGIONS or k[1] not in SUBTYPES]
    assert not stale, f"provenance records for non-existent REGIONS entries: {stale}"
    assert not [s for s in EXTRASYN_PROV if s not in EXTRASYN]
    assert not [s for s in SUBJECTIVE_PROV if s not in SUBJECTIVE_WEIGHT]


# ============================================================== the records must be honest
def test_sourced_records_actually_name_a_source():
    """`is_sourced` must not be satisfiable by a basis alone. A record claiming a traceable
    derivation with no source key is exactly the failure this module exists to prevent."""
    for name, table in ALL.items():
        for key, rec in table.items():
            if rec.is_sourced:
                assert rec.source_key, f"{name}.{key} claims to be sourced with no key"


def test_every_named_source_exists_in_the_knowledge_base():
    """A provenance record pointing at a source key that is not in the database is a
    dangling citation -- worse than no citation, because it reads as traceable."""
    if not DB.exists():
        pytest.skip("pharmacology.db not present")
    con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    keys = {k for (k,) in con.execute("select key from sources")}
    con.close()
    dangling = sorted({rec.source_key for t in ALL.values() for rec in t.values()
                       if rec.source_key and rec.source_key not in keys})
    assert not dangling, f"provenance names sources absent from the database: {dangling}"


def test_every_record_explains_itself():
    """An UNSOURCED record with no note is indistinguishable from an oversight."""
    for name, table in ALL.items():
        for key, rec in table.items():
            assert rec.note.strip(), f"{name}.{key} has no note explaining its basis"


# ================================================= the headline count, pinned so it moves visibly
def test_the_provenance_gap_is_what_the_audit_recorded():
    """Pinned deliberately. The gap is the finding, so it should be impossible to close it
    quietly OR to widen it without the number changing in a diff.

    As audited 2026-10-07: 6 of 28 named a source (21%), 18 UNSOURCED, 8 FROM_QUALITATIVE,
    1 FITTED, 1 GUESS.

    UPDATED 2026-10-08, and the gap WIDENED: 19 UNSOURCED, 7 FROM_QUALITATIVE. `EXTRASYN[a5]`
    = 0.80 was downgraded FROM_QUALITATIVE -> UNSOURCED after a claim-support read of its
    candidate source. Kasugai et al. 2010 is the obvious citation for "a5 is predominantly
    extrasynaptic" -- right journal, right method (quantitative freeze-fracture replica
    immunogold), right preparation -- and its abstract shows it measured a1, a2 and b3, with
    synaptic labelling density 78-132x extrasynaptic. It is silent on a5.

    That is the SECOND source in this project to resolve perfectly by DOI and fail to support
    the number attached to it (`a5_dist` was the first). The pattern: a source whose title
    matches the claim, in the right journal, by the right group, measuring a neighbouring
    quantity. Metadata verification cannot catch it.

    The `sourced` count is unchanged at 6 because `EXTRASYN[a5]` never carried a source KEY --
    it was FROM_QUALITATIVE with an empty key, which is itself worth noting: a basis label can
    assert more confidence than any recorded citation supports.

    UPDATED 2026-10-09 (roadmap P0-5/P0-6/P0-13), and the gap widened again, from 28
    parameters to 43. Nothing was demoted this time; 15 numbers that were ALREADY
    load-bearing were brought under the audit for the first time:

      * MODEL_PARAMS (14). The three competing receptor models in `models/` carry
        dataclass defaults that both entry points to the three-tiered headroom hierarchy
        used to default to. `KineticAllosteryModel`'s are a chimera of two different fits
        -- beta/alpha from the manuscript's config-pulse fit, kon/koff from the superseded
        1000 uM / 0.30 ms fit -- giving K_d = 29.62 uM, a value in no fit, no commit and no
        document. Two of Model A's defaults are the project's own FIT TARGETS (po_max 0.75,
        tau_deact 15.0 ms) rather than anything fitted, so they are recorded FITTED and are
        circular if used as anchors.
      * DATASETS (2). Both `fitting/data.py` benchmarks are synthetic. The deactivation
        trace carried a citation to Haas & Macdonald 1999, which this project's own
        claim-support audit records as measuring 76.1 ms for that quantity against the
        trace's dominant 15 ms component.

    One of the 14 is `gabaa_kinetics.FIT_FIXED_ALPHA`, the first CONVENTION-rated entry in
    the audit: the channel closing rate held fixed so the kinetic fit is well-posed. It is
    deliberately NOT in RANKING_INPUT_TABLES -- it is a kinetic fit convention, not an
    input to the selectivity ranking, and putting it there would have moved the figure the
    manuscript quotes.

    So the sourced FRACTION fell from 21% to 13% without a single citation being lost: the
    denominator grew because the audit's coverage grew. That is the honest direction, and
    it is the reason this pin names the composition rather than only the total.
    """
    from circuitpharm.provenance import RANKING_INPUT_TABLES

    # TWO SCOPES, both pinned. The ranking's inputs are what the manuscript quotes and
    # what the project's one VALIDATED result rests on; the full audit additionally covers
    # the receptor-model defaults and the benchmark datasets, which feed the
    # receptor-kinetic results instead. Conflating them would let either drift behind the
    # other.
    rank = audit(RANKING_INPUT_TABLES)
    assert (rank["sourced"], rank["total"]) == (6, 28), (
        f"the ranking-input scope is now {rank['sourced']} of {rank['total']}, not 6 of "
        f"28. knowledge/08-manuscript.md quotes this figure, so a change here needs the "
        f"manuscript updated too.")

    a = audit()
    assert a["total"] == 44, f"parameter count changed to {a['total']}; re-run the audit"
    assert a["sourced"] == 6, (
        f"{a['sourced']} parameters now name a source, not 6. If this went UP, update this "
        f"pin and knowledge/06-source-provenance.md. If it went DOWN, something lost its "
        f"citation.")
    assert a["by_basis"][Basis.UNSOURCED.value] == 31, (
        f"{a['by_basis'][Basis.UNSOURCED.value]} UNSOURCED, not 31. Widening this gap is a "
        f"legitimate finding -- it means a claim-support read demoted something, or the "
        f"audit's coverage grew -- but it must be recorded in "
        f"knowledge/06-source-provenance.md and in this docstring, not just absorbed by "
        f"the pin.")
    assert a["by_basis"][Basis.FROM_QUALITATIVE.value] == 7
    assert a["by_basis"][Basis.GUESS.value] == 2, (
        "2 GUESS entries: EXTRASYN['eps'] (labelled a guess in subtypes.py itself) and "
        "ExtendedDesensitizationModel.pam_desens_factor (neither the 0.1 coefficient nor "
        "the linear form has a source).")
    assert a["by_basis"][Basis.FITTED.value] == 3, (
        "3 FITTED: the GAIN_RATIO prior, plus Model A's po_max and tau_deact_ms, which are "
        "the project's own fit targets rather than independent measurements.")
    assert a["by_basis"][Basis.CONVENTION.value] == 1, (
        "1 CONVENTION: gabaa_kinetics.FIT_FIXED_ALPHA. A convention is a legitimate basis "
        "-- it is a declared choice rather than an unsourced guess -- but it must stay "
        "countable, because this one determines every derived kinetic number.")


def test_the_most_consequential_number_is_flagged_unsourced():
    """preBotC a5. An a5-selective compound's modelled respiratory burden is roughly
    proportional to it, so it sets the entire safety margin -- and no source in the
    knowledge base states it. If this ever becomes sourced, that is a real advance and this
    test should be updated to say so."""
    rec = REGIONS_PROV[("prebotc", "a5")]
    assert rec.basis is Basis.UNSOURCED
    assert not rec.source_key


def test_a5_dist_is_not_cited_for_the_regional_fractions():
    """a5_dist is the obvious candidate source for REGIONS and it does NOT support them: a
    1988 study with a single probe for 'the alpha subunit', reporting TOTAL alpha mRNA by
    region, which neither resolves subtypes nor measures composition. It resolves perfectly
    by DOI, so a metadata-only check would have marked it verified.

    Citing it would be worse than citing nothing, so this pins that it is not cited.
    """
    cited = {rec.source_key for rec in REGIONS_PROV.values() if rec.source_key}
    assert "a5_dist" not in cited


def test_the_robustness_priors_carry_provenance_too():
    """KAPPA and EXTRASYN_CONC decide what 'robust across 20,000 draws' means. They are
    arguably more load-bearing than the point values, since the draws are what the headline
    claim rests on -- and both are unsourced choices."""
    assert PRIOR_PROV["KAPPA"].basis is Basis.UNSOURCED
    assert PRIOR_PROV["EXTRASYN_CONC"].basis is Basis.UNSOURCED
    assert PRIOR_PROV["GAIN_RATIO"].basis is Basis.FITTED


def test_report_renders():
    txt = report()
    assert "PARAMETER PROVENANCE" in txt and "UNSOURCED" in txt
