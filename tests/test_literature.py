"""The literature readings, pinned.

WHY THESE TESTS EXIST. `scripts/literature_gamma2.py` holds numbers transcribed by hand
from seven full texts. Hand transcription is exactly what this project distrusts elsewhere,
so the readings are pinned here and the DERIVED quantities are recomputed rather than
restated.

SCOPE, deliberately narrow. These tests assert:
  * that the derived means are what the published components actually give,
  * that the structural facts other modules will rely on hold (which quantity is
    macroscopic, which sources are gamma2L, that no mean open time is stated for gamma2L),
  * that every source carries the provenance fields the project's rules require.

They do NOT assert the prose of knowledge/14-literature.md, and they do NOT forbid a reading
from being corrected -- if a full text is re-read and a number moves, the fix belongs in the
script and these values move with it. What they forbid is a number moving SILENTLY, and a
new source being added without its conditions.
"""
import math
import sqlite3

import pytest

from scripts.literature_gamma2 import (
    DEACTIVATIONS, NOMINAL_DEFECT, OPEN_TIME_FITS, OTHER, PO_CANDIDATES, SOURCES,
    STATED_MEAN_OPEN_TIME, _weighted,
)

GAMMA2L_SOURCES = {"li2008", "dixon2014", "jahn_a1b2g2_kinetics"}


# ---------------------------------------------------------------------------------------
# Provenance completeness. The rule this project runs on: a number without its conditions
# is not a reading.
# ---------------------------------------------------------------------------------------

def test_every_source_identifies_itself():
    for key, s in SOURCES.items():
        assert s.key == key
        assert s.url.startswith("https://"), key
        assert s.doi or s.pmid, f"{key} has neither DOI nor PMID"
        assert s.citation.strip(), key
        assert s.read in {"FULL_TEXT_PMC", "ABSTRACT_PUBMED", "NOT_READ"}, key


def test_every_source_states_its_conditions_and_error_convention():
    """A source read in full text must say what receptor, prep, temperature and error
    statistic it reports -- including saying 'UNSTATED' where the paper does not."""
    for key, s in SOURCES.items():
        if s.read != "FULL_TEXT_PMC":
            continue
        assert s.receptor_verbatim, key
        assert s.prep, key
        assert s.temp_c, key
        assert s.err_kind in {"SEM", "SD", "UNSTATED"}, f"{key}: {s.err_kind!r}"


def test_the_unstated_error_conventions_are_recorded_not_assumed():
    """Keramidas 2008 and Dixon 2014 declare no error statistic. Assuming SEM would be the
    easy move and would be unsupported, so it is not made."""
    assert SOURCES["keramidas2008"].err_kind == "UNSTATED"
    assert SOURCES["dixon2014"].err_kind == "UNSTATED"
    # ...while the ones that DO declare are recorded as declared.
    assert SOURCES["keramidas2010"].err_kind == "SEM"
    assert SOURCES["barberis2007"].err_kind == "SEM"
    assert SOURCES["li2008"].err_kind == "SD"


def test_sources_are_registered_in_the_database():
    c = sqlite3.connect("data/pharmacology.db")
    present = {k for (k,) in c.execute("select key from sources")}
    for key in SOURCES:
        assert key in present, f"{key} read but not registered in the sources table"


def test_claim_support_is_asserted_only_where_something_was_actually_read():
    """A source that was READ -- in full text, or as an abstract -- may assert
    CLAIM_SUPPORT_SUPPORTS for the claims that text states. An UNREAD source may not.

    This assertion was initially written as "only full text may assert claim support",
    which was wrong: an abstract genuinely does support the claims it itself states, and
    Jahn 1997's EC50 and Hill slope are read verbatim from its abstract. The real
    distinction is read versus unread, not full text versus abstract. What an abstract
    cannot support is a claim it does not make -- which is why the Jahn row carries an
    explicit CLAIM_SUPPORT_DOES_NOT_SUPPORT for mean open time.
    """
    c = sqlite3.connect("data/pharmacology.db")
    for key, src in SOURCES.items():
        row = c.execute("select verification from sources where key=?", (key,)).fetchone()
        assert row, f"{key} is not in the sources table"
        (verif,) = row
        assert "HAND:" in verif, f"{key} has no hand verdict"
        if src.read == "NOT_READ":
            assert "CLAIM_SUPPORT_SUPPORTS" not in verif, (
                f"{key} was never read but claims claim-support")
            assert "CLAIM_SUPPORT_UNASSESSED" in verif, key
        else:
            assert "CLAIM_SUPPORT_SUPPORTS" in verif, key


def test_sources_whose_limits_were_found_record_what_they_do_not_support():
    """The negatives are findings. Each of these sources was read and turned out NOT to
    support something the project expected of it, and the row has to say so."""
    c = sqlite3.connect("data/pharmacology.db")
    for key in ("jahn_a1b2g2_kinetics", "keramidas2008", "dixon2014", "dixon2015",
                "li2008", "mortensen2012", "barberis2007"):
        (verif,) = c.execute("select verification from sources where key=?", (key,)).fetchone()
        assert "CLAIM_SUPPORT_DOES_NOT_SUPPORT" in verif, (
            f"{key} was read and its limits were found, but the row does not record them")


# ---------------------------------------------------------------------------------------
# Open time. The derived means are recomputed, never restated.
# ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("fit", OPEN_TIME_FITS, ids=lambda f: f.label)
def test_open_time_areas_are_a_partition(fit):
    assert fit.areas == pytest.approx([a for a in fit.areas])
    assert sum(fit.areas) == pytest.approx(1.0, abs=5e-3), (
        f"{fit.label}: published areas sum to {sum(fit.areas)}, which should be 1")
    assert len(fit.taus_ms) == len(fit.areas) == len(fit.tau_errs)
    assert all(t > 0 for t in fit.taus_ms)
    assert list(fit.taus_ms) == sorted(fit.taus_ms), "components should be ordered fast->slow"


def test_derived_mean_open_times_are_what_the_components_give():
    got = {f.label: _weighted(f.taus_ms, f.areas) for f in OPEN_TIME_FITS}
    assert got["Keramidas 2008 M-Mode"] == pytest.approx(2.684, abs=5e-3)
    assert got["Keramidas 2008 H-Mode"] == pytest.approx(7.249, abs=5e-3)
    assert got["Li 2008 / Li 2007b control"] == pytest.approx(2.961, abs=5e-3)


def test_no_source_states_a_mean_open_time_for_gamma2L():
    """The headline negative result of the literature pass. Barberis is the only source that
    STATES a mean open time, and it is gamma2S."""
    for key, d in STATED_MEAN_OPEN_TIME.items():
        assert not d["receptor"].endswith("gamma2L"), (
            f"{key} now states a gamma2L mean open time -- if a new source was read, "
            "MISSING_DATASETS['single_channel_mean_open_time'] must be updated too")
    assert set(STATED_MEAN_OPEN_TIME) == {"barberis2007"}


def test_the_two_intraburst_readings_agree_across_splice_variants():
    """The one reassuring cross-check in the corpus: gamma2S at 10 mM (M-mode, intraburst)
    and gamma2L at 50 uM (intracluster) land within 0.3 ms of each other."""
    by_label = {f.label: f for f in OPEN_TIME_FITS}
    m = _weighted(*[getattr(by_label["Keramidas 2008 M-Mode"], a) for a in ("taus_ms", "areas")])
    li = _weighted(*[getattr(by_label["Li 2008 / Li 2007b control"], a)
                     for a in ("taus_ms", "areas")])
    assert abs(m - li) < 0.35, (m, li)
    # and they are NOT the same receptor, which is the point
    assert by_label["Keramidas 2008 M-Mode"].receptor.endswith("gamma2S")
    assert by_label["Li 2008 / Li 2007b control"].receptor.endswith("gamma2L")


def test_mean_open_time_spread_is_large_enough_to_matter():
    """If this ever collapses, the 'mean open time is not one quantity' claim in
    knowledge/14-literature.md and MISSING_DATASETS has stopped being true."""
    vals = [_weighted(f.taus_ms, f.areas) for f in OPEN_TIME_FITS]
    vals += [d["value_ms"] for d in STATED_MEAN_OPEN_TIME.values()]
    assert max(vals) / min(vals) > 4.0, vals


# ---------------------------------------------------------------------------------------
# Deactivation. The protocol dependence, and the read-back check on the transcription.
# ---------------------------------------------------------------------------------------

def test_the_dixon_lead_is_confirmed_verbatim():
    """This value was [LEAD]-graded and unusable until the full text was read."""
    d = next(d for d in DEACTIVATIONS if d.label == "Dixon 2014")
    assert d.tau_w_ms == 5.9 and d.err_ms == 0.5 and d.n == "10"
    assert d.receptor.endswith("gamma2L")
    assert d.pulse == "<=1 ms" and d.conc_mm == 3.0
    assert not d.components, (
        "Dixon 2014 does not publish the components behind its weighted tau; if that "
        "changed, the full text was re-read and MISSING_DATASETS should say so")


def test_barberis_components_reproduce_its_published_weighted_tau():
    """The read-back check. If the transcription of either the components or the weighted
    value were wrong, these would not agree."""
    d = next(d for d in DEACTIVATIONS if d.label == "Barberis 2007, brief pulse")
    taus, areas = d.components
    assert sum(areas) == pytest.approx(1.0, abs=5e-3)
    recomputed = _weighted(taus, areas)
    assert recomputed == pytest.approx(53.55, abs=0.05)
    assert abs(recomputed - d.tau_w_ms) < d.err_ms, (
        f"recomputed {recomputed} vs published {d.tau_w_ms} +/- {d.err_ms}")


def test_deactivation_is_protocol_dependent_within_a_single_study():
    """Same receptor, same patches, same study -- only the pulse duration differs. This is a
    published measurement, which is why it is stronger evidence than the PEAK result."""
    brief = next(d for d in DEACTIVATIONS if d.label == "Barberis 2007, brief pulse")
    long = next(d for d in DEACTIVATIONS if d.label == "Barberis 2007, long pulse")
    assert brief.source_key == long.source_key == "barberis2007"
    assert brief.receptor == long.receptor
    assert brief.conc_mm == long.conc_mm
    assert long.tau_w_ms / brief.tau_w_ms > 5.0, (brief.tau_w_ms, long.tau_w_ms)


def test_the_dixon_barberis_gap_is_not_explained_by_pulse_duration():
    """<=1 ms vs 2 ms is nearly the same protocol, yet the weighted taus differ ~9x. The
    explanation is the component count, and Dixon's weighted value sits near Barberis's
    FAST component rather than its weighted one."""
    dixon = next(d for d in DEACTIVATIONS if d.label == "Dixon 2014")
    barb = next(d for d in DEACTIVATIONS if d.label == "Barberis 2007, brief pulse")
    assert barb.tau_w_ms / dixon.tau_w_ms > 8.0
    tau_fast = barb.components[0][0]
    assert abs(dixon.tau_w_ms - tau_fast) < abs(dixon.tau_w_ms - barb.tau_w_ms)
    # the slow component carries most of Barberis's weighted value
    slow_contrib = barb.components[0][2] * barb.components[1][2]
    assert slow_contrib / barb.tau_w_ms > 0.75, slow_contrib


def test_the_simulated_row_is_marked_as_not_a_measurement():
    d = next(d for d in DEACTIVATIONS if "SIMULATED" in d.label)
    assert "SIMULATION" in d.where.upper() or "SIMULAT" in d.where.upper()
    assert math.isnan(d.err_ms), "a simulated point value must not carry an error bar"
    assert "Do not use as data" in d.note


# ---------------------------------------------------------------------------------------
# P_o. Task 4: which quantity is 0.750 meant to be.
# ---------------------------------------------------------------------------------------

def test_exactly_one_po_candidate_is_macroscopic():
    macro = [c for c in PO_CANDIDATES if "MACROPATCH" in c["quantity"]]
    assert len(macro) == 1, [c["quantity"] for c in PO_CANDIDATES]
    assert macro[0]["value"] == 0.56
    assert "nothing" in macro[0]["conditions_on"]
    assert "UNPUBLISHED" in macro[0]["note"].upper()


def test_every_other_po_candidate_declares_what_it_conditions_on():
    for c in PO_CANDIDATES:
        assert c["conditions_on"], c
        assert c["quantity"], c
        assert c["source_key"] in SOURCES, c
        if "MACROPATCH" not in c["quantity"]:
            assert ("INTRABURST" in c["quantity"] or "CLUSTER" in c["quantity"]), c


def test_the_only_macroscopic_reading_is_below_the_convention():
    """The substantive Task-4 result. Every reading at or above 0.750 is conditional; the
    unconditional one is well below it. So the measurements do not support RAISING the
    convention, and the one that could lower it is at an incomparable protocol."""
    macro = next(c for c in PO_CANDIDATES if "MACROPATCH" in c["quantity"])
    assert macro["value"] < NOMINAL_DEFECT
    above = [c for c in PO_CANDIDATES if c["value"] >= NOMINAL_DEFECT]
    assert above, "expected at least one candidate at or above the convention"
    assert all("INTRABURST" in c["quantity"] for c in above), (
        "a non-conditional reading now sits at or above the convention -- Task 4's "
        "conclusion in knowledge/14-literature.md would need revisiting")


def test_the_0_69_is_attributed_to_the_right_measurement():
    """The handoff recorded 0.69 as 'nonstationary variance analysis'. It is not: it is an
    intraburst M-mode P_o, and the actual fluctuation-analysis value is 0.56. This test
    exists so the conflation cannot come back."""
    c = next(c for c in PO_CANDIDATES
             if c["value"] == 0.69 and c["source_key"] == "keramidas2008")
    assert "INTRABURST" in c["quantity"]
    assert "M-mode" in c["quantity"]
    assert "NOT" in c["note"]


def test_gamma2L_po_readings_are_both_conditional():
    g = [c for c in PO_CANDIDATES if c["receptor"].endswith("gamma2L")]
    assert len(g) == 2
    assert all("INTRABURST" in c["quantity"] or "CLUSTER" in c["quantity"] for c in g)


# ---------------------------------------------------------------------------------------
# The negatives. These are findings, and they are the ones most likely to be quietly lost.
# ---------------------------------------------------------------------------------------

def test_jahn_is_recorded_as_abstract_only_with_no_mean_open_time():
    s = SOURCES["jahn_a1b2g2_kinetics"]
    assert s.read == "ABSTRACT_PUBMED"
    assert "NO mean open time" in s.note
    assert "STILL UNKNOWN" in s.note
    burst = next(o for o in OTHER
                 if o["source_key"] == "jahn_a1b2g2_kinetics" and "burst" in o["quantity"])
    assert burst["where"] == "ABSTRACT"
    assert "NOT a mean open time" in burst["note"]


def test_hamill_is_recorded_as_unread():
    s = SOURCES["hamill1997"]
    assert s.read == "NOT_READ"
    assert s.pmid == "9480006"
    assert not s.doi, "the Hamill comment has no DOI; inventing one would be a fabrication"
    assert "UNREAD" in s.note.upper()


def test_mortensen_is_recorded_as_having_no_hill_coefficients():
    """MISSING_DATASETS and HANDOFF.md called it an 'EC50/nH table'. The nH half does not
    exist, and a case-insensitive grep for 'nH' hits only author initials."""
    assert "NO Hill coefficients" in SOURCES["mortensen2012"].note
    ec50 = next(o for o in OTHER if o["source_key"] == "mortensen2012")
    assert "No Hill coefficient is reported" in ec50["note"]


def test_the_holdout_records_the_sources_internal_inconsistency():
    """Barberis's text says n = 7 for the desensitisation onset; its Figure 4 caption says
    three patches. Recorded as an ambiguity rather than resolved by picking one."""
    o = next(o for o in OTHER if "desensitisation onset" in o["quantity"])
    assert "7" in o["n"] and "3" in o["n"]
    assert "INCONSISTENT" in o["n"].upper()


def test_splice_variant_substitutions_are_visible():
    """Only two of the full texts are the project's receptor. If that is ever miscounted the
    substitutions stop being recorded, which is the specific failure MISSING_DATASETS was
    written to prevent."""
    gamma2l = {k for k, s in SOURCES.items() if s.receptor_verbatim.endswith("gamma2L")}
    assert gamma2l == GAMMA2L_SOURCES, gamma2l
    full_text_gamma2l = {k for k in gamma2l if SOURCES[k].read == "FULL_TEXT_PMC"}
    assert full_text_gamma2l == {"li2008", "dixon2014"}, full_text_gamma2l


def test_the_corrected_citations_stay_corrected():
    """Three citation errors came from agent-supplied metadata. Pinned so a later edit
    cannot reintroduce them."""
    assert "289(9)" in SOURCES["dixon2014"].journal_ref
    assert "289(8)" not in SOURCES["dixon2014"].journal_ref
    assert "25(9):2726-2740" in SOURCES["barberis2007"].journal_ref
    assert "26(7)" not in SOURCES["barberis2007"].journal_ref
    assert SOURCES["dixon2015"].receptor_verbatim.endswith("gamma2S")
    assert "iv" in SOURCES["hamill1997"].journal_ref


def test_li_control_data_is_attributed_to_its_actual_origin():
    """Li 2008's GABA-alone open times are reproduced from Li et al. 2007b. Citing 2008 for
    them would be a mis-citation of exactly the kind this project tracks."""
    assert "Li et al. (2007b)" in SOURCES["li2008"].note
    assert "jphysiol.2007.142794" in SOURCES["li2008"].note
