"""The literature readings, pinned.

WHY THESE TESTS EXIST. `circuitpharm.literature` holds numbers transcribed by hand
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

# From the PACKAGE, not from scripts/. These constants used to live in
# scripts/literature_gamma2.py, which made this module's import depend on the repo root being
# on sys.path -- it passed under `python -m pytest` and died in CI under bare `pytest`
# (fixed once in 5643f74 by pyproject's pythonpath, and fixed at the root here by moving the
# data into the package). `fitting.data` imported them from scripts/ too, which broke
# `import circuitpharm` from any directory other than the repo root.
from circuitpharm.literature import (
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


def test_deactivation_and_open_time_datasets_enforce_protocol_declarations():
    """Phase 1: datasets enforce positive protocol declarations and valid conditioning."""
    import numpy as np
    from circuitpharm.models.base import Observable, ObservableMismatch
    from circuitpharm.fitting.data import DeactivationDataset, MeanOpenTimeDataset
    from circuitpharm.fitting.comparison import check_deactivation_protocols

    # DeactivationDataset missing pulse_ms raises ValueError
    with pytest.raises(ValueError, match="deactivation requires positive protocol declarations"):
        DeactivationDataset(
            citation_label="invalid deact",
            time_ms=np.linspace(0, 100, 10),
            normalized_current=np.linspace(1, 0, 10),
            observable=Observable.DEACTIVATION,
            conc_mm=1.0,
            n_components=2,
            fit_window_ms=100.0,
            # missing pulse_ms
        )

    # MeanOpenTimeDataset with missing or invalid conditioning raises ValueError
    with pytest.raises(ValueError, match="mean open time requires `conditioning`"):
        MeanOpenTimeDataset(
            citation_label="invalid open time",
            mean_open_time_ms=2.5,
            conditioning="INVALID",  # type: ignore
        )

    # Comparing Dixon (5.9 ms, 1 ms pulse, 2 components) with Barberis (52.5 ms, 2 ms pulse, 3 components)
    # raises ObservableMismatch citing differing pulse_ms and n_components
    dixon_ds = DeactivationDataset(
        citation_label="Dixon 2014",
        pulse_ms=1.0,
        conc_mm=3.0,
        n_components=2,
        fit_window_ms=50.0,
        time_ms=np.linspace(0, 50, 10),
        normalized_current=np.linspace(1, 0, 10),
        observable=Observable.DEACTIVATION,
        synthetic=True,
        kind="synthetic",
        origin="test fixture",
    )
    barberis_ds = DeactivationDataset(
        citation_label="Barberis 2007",
        pulse_ms=2.0,
        conc_mm=10.0,
        n_components=3,
        fit_window_ms=1000.0,
        time_ms=np.linspace(0, 1000, 10),
        normalized_current=np.linspace(1, 0, 10),
        observable=Observable.DEACTIVATION,
        synthetic=True,
        kind="synthetic",
        origin="test fixture",
    )
    with pytest.raises(ObservableMismatch) as exc_info:
        check_deactivation_protocols(dixon_ds, barberis_ds)
    err = str(exc_info.value)
    assert "pulse duration" in err
    assert "n_components" in err


def test_parametric_datasets_declare_source_measured_domain_and_quantify_extrapolation():
    """Task B: parametric datasets must declare source_measured_domain.

    Jahn spans 0.3 -> 3000 uM while measured only over 1-10 uM slope window:
    n_points_outside_source_domain is 6 of 9 (fraction 6/9).
    Barberis spans 0 -> 1000 ms and fits within that window:
    n_points_outside_source_domain is 0 of 101 (fraction 0.0).

    The finding is the asymmetry: fraction_outside_source_domain separates them.
    A parametric dataset constructed without source_measured_domain raises ValueError.
    """
    import numpy as np
    from circuitpharm.fitting.data import (
        BARBERIS2007_DEACT,
        JAHN1997_PEAK_CRC,
        DeactivationDataset,
        DoseResponseDataset,
    )
    from circuitpharm.models.base import Observable

    # 1. Jahn: 6 of 9 points outside (1.0, 10.0)
    assert JAHN1997_PEAK_CRC.source_measured_domain == (1.0, 10.0)
    assert len(JAHN1997_PEAK_CRC.concs_um) == 9
    assert JAHN1997_PEAK_CRC.n_points_outside_source_domain == 6
    assert JAHN1997_PEAK_CRC.fraction_outside_source_domain == pytest.approx(6 / 9)

    # 2. Barberis: 0 of 101 points outside (0.0, 1000.0)
    assert BARBERIS2007_DEACT.source_measured_domain == (0.0, 1000.0)
    assert len(BARBERIS2007_DEACT.time_ms) == 101
    assert BARBERIS2007_DEACT.n_points_outside_source_domain == 0
    assert BARBERIS2007_DEACT.fraction_outside_source_domain == 0.0

    # 3. Assert asymmetry directly
    assert JAHN1997_PEAK_CRC.fraction_outside_source_domain > 0.6
    assert BARBERIS2007_DEACT.fraction_outside_source_domain == 0.0

    # 4. Parametric datasets built without source_measured_domain raise ValueError
    with pytest.raises(ValueError, match="source_measured_domain"):
        DoseResponseDataset(
            citation_label="unconstrained parametric dose-response",
            preparation="HEK293",
            concs_um=np.array([1.0, 10.0]),
            mean_response=np.array([0.1, 0.9]),
            sem=np.array([0.01, 0.01]),
            observable=Observable.PEAK,
            kind="parametric",
            synthetic=True,
            source_key="jahn_a1b2g2_kinetics",
            origin="test fixture",
            # source_measured_domain omitted
        )

    with pytest.raises(ValueError, match="source_measured_domain"):
        DeactivationDataset(
            citation_label="unconstrained parametric deactivation",
            time_ms=np.array([0.0, 10.0]),
            normalized_current=np.array([1.0, 0.5]),
            pulse_ms=2.0,
            conc_mm=10.0,
            n_components=1,
            fit_window_ms=10.0,
            observable=Observable.DEACTIVATION,
            kind="parametric",
            synthetic=True,
            source_key="barberis2007",
            origin="test fixture",
            # source_measured_domain omitted
        )


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


def test_po_max_convention_matches_what_the_code_records():
    """Phase 5: P_o,max stays 0.750 and the project records WHY the literature cannot
    replace it.

    THIS TEST REPLACES ONE THAT WAS NAMED FOR AN INVARIANT IT DID NOT CHECK, which is worth
    recording. The first version was called
    `test_nominal_defect_stays_convention_and_unreplaced` and asserted a float and a
    substring against two module constants -- `parameters.NOMINAL_DEFECT` and
    `parameters.PO_MAX_CONVENTION` -- that the same change had just added, that nothing in
    the model consumed, and that were hard-coded 0.750 literals rather than derived from
    `FIT_TARGETS`. So it would have passed unchanged if the real value had moved, and it
    never checked the `Basis` its own name claimed. The roadmap meanwhile asserted
    `Basis.CONVENTION` while the code recorded `Basis.FITTED`.

    What is asserted instead: the value has ONE definition, the recorded basis is whatever
    the code actually says, and the verdict's reasoning is present where a reader will find
    it. The basis is pinned rather than demanded -- if it is deliberately regraded to
    CONVENTION (see roadmap P0-13), this test and the roadmap must change together, which is
    the point.
    """
    from circuitpharm import parameters
    from circuitpharm.gabaa_kinetics import FIT_RANGES, FIT_TARGETS
    from circuitpharm.provenance import Basis, MODEL_PARAM_PROV

    # one definition, and it is the canonical one
    assert FIT_TARGETS["po_max"] == 0.75
    for dead in ("NOMINAL_DEFECT", "PO_MAX_CONVENTION", "PO_MAX_PROVENANCE"):
        assert not hasattr(parameters, dead), (
            f"parameters.{dead} is a second home for a value that already has one at "
            "gabaa_kinetics.FIT_TARGETS['po_max'] -- that is E12/E22 re-forming")

    # the convention was NOT replaced by the 0.56 reading
    assert FIT_TARGETS["po_max"] != 0.56
    lo, hi = FIT_RANGES["po_max"]
    assert not (lo <= 0.56 <= hi), (
        "0.56 is now inside FIT_RANGES, which removes one of the four stated reasons the "
        "literature cannot replace the convention -- roadmap P0-13 needs revisiting")

    # the recorded basis, pinned to what the code says rather than to what prose claims
    rec = MODEL_PARAM_PROV[("operational", "po_max")]
    assert rec.basis is Basis.FITTED, (
        f"po_max's basis is now {rec.basis.value}. If that was deliberate, update roadmap "
        "P0-13's 'OPEN' paragraph in the same commit -- it currently records FITTED and "
        "argues CONVENTION may be more honest.")

    # and the verdict's reasoning lives where a reader of the defaults will meet it
    note = parameters._NOMINAL_DEFECT["operational"]
    assert "CANNOT REPLACE" in note
    assert "0.56" in note and "UNPUBLISHED DATA" in note
    assert "14-literature" in note


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


def test_paired_pulse_holdout_registration_and_validation():
    """Task D: Barberis 2007 paired-pulse recovery is registered as the first real holdout.

    Observable.PAIRED_PULSE_RECOVERY requires pulse_ms, conc_mm and gap_ms.
    The dataset is scalar (0.33 +/- 0.03 at 100 ms gap, 2 ms pulses, 10 mM GABA, n=8).
    It is digitised, synthetic=False, role='holdout', and records gamma2S substitution.
    OperationalScalarModel raises ObservableMismatch (no q_matrix).
    Markov schemes predict finite recoveries and score against holdout without taint.
    """
    import warnings
    import numpy as np
    from circuitpharm.fitting.data import (
        ALL,
        HOLDOUT,
        BARBERIS2007_PAIRED_PULSE,
        PairedPulseDataset,
        assert_real_data,
    )
    from circuitpharm.fitting.likelihood import score_holdout, _predict
    from circuitpharm.models import (
        ExtendedDesensitizationModel,
        KineticAllosteryModel,
        OperationalScalarModel,
    )
    from circuitpharm.models.base import Observable, ObservableMismatch
    from circuitpharm.provisional import ProvisionalResultWarning

    # 1. Registered in ALL and HOLDOUT
    assert "barberis2007_paired_pulse" in ALL
    assert ALL["barberis2007_paired_pulse"] is BARBERIS2007_PAIRED_PULSE
    assert len(HOLDOUT) == 1
    assert HOLDOUT[0] is BARBERIS2007_PAIRED_PULSE

    # 2. Dataset properties and gamma2S recording
    ds = BARBERIS2007_PAIRED_PULSE
    assert ds.observable is Observable.PAIRED_PULSE_RECOVERY
    assert ds.kind == "digitised"
    assert ds.synthetic is False
    assert ds.role == "holdout"
    assert ds.source_key == "barberis2007"
    assert ds.mean_response == pytest.approx(0.33)
    assert ds.sem == pytest.approx(0.03)
    assert ds.noise_model == "stated_sem"
    assert ds.pulse_ms == 2.0
    assert ds.conc_mm == 10.0
    assert ds.gap_ms == 100.0
    assert ds.n_cells == 8
    assert "gamma2S" in ds.receptor
    assert "gamma2L" in ds.receptor
    assert "gamma2S" in ds.preparation

    # 3. Protocol declaration requirements: missing or non-positive fields raise ValueError
    for kw, match in (
        (dict(mean_response=0.33, conc_mm=10.0, gap_ms=100.0), "pulse_ms"),
        (dict(mean_response=0.33, pulse_ms=2.0, gap_ms=100.0), "conc_mm"),
        (dict(mean_response=0.33, pulse_ms=2.0, conc_mm=10.0), "gap_ms"),
        (dict(mean_response=0.33, pulse_ms=0.0, conc_mm=10.0, gap_ms=100.0), "pulse_ms"),
        (dict(mean_response=0.33, pulse_ms=2.0, conc_mm=10.0, gap_ms=100.0, noise_model="invalid"), "noise_model"),
        (dict(mean_response=0.33, pulse_ms=2.0, conc_mm=10.0, gap_ms=100.0, sem=None, noise_model="stated_sem"), "no sem"),
    ):
        with pytest.raises(ValueError, match=match):
            PairedPulseDataset(
                citation_label="invalid protocol",
                kind="synthetic",
                synthetic=True,
                origin="test fixture",
                **kw,
            )

    # 4. Guard on q_matrix: OperationalScalarModel raises ObservableMismatch
    op = OperationalScalarModel()
    with pytest.raises(ObservableMismatch, match="no `q_matrix`|no states to relax between"):
        _predict(op, ds)

    # 5. First non-synthetic dataset does not taint with ProvisionalResultWarning
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fake = assert_real_data([ds], caller="test")
        assert fake == []
    prov = [w for w in caught if issubclass(w.category, ProvisionalResultWarning)]
    assert len(prov) == 0

    # 6. Scoring against holdout: Markov models predict finite recoveries and score cleanly
    b = KineticAllosteryModel()
    pred_b = float(_predict(b, ds))
    assert 0.0 < pred_b < 1.0

    c = ExtendedDesensitizationModel()
    pred_c = float(_predict(c, ds))
    assert 0.0 < pred_c < 1.0

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res_b = score_holdout(b, HOLDOUT, caller="test")
        res_c = score_holdout(c, HOLDOUT, caller="test")
    prov = [w for w in caught if issubclass(w.category, ProvisionalResultWarning)]
    assert len(prov) == 0
    assert np.isfinite(res_b.total_log_likelihood)
    assert np.isfinite(res_c.total_log_likelihood)
    # Stated-SEM scoring: nominal defaults miss by >17x SEM, scoring strongly negative
    assert res_b.total_log_likelihood < -100.0
    assert res_c.total_log_likelihood < -100.0
    fit_b = res_b.per_dataset[0]
    assert fit_b.noise_model == "stated_sem"
    assert fit_b.stated_sem == 0.03
    assert fit_b.sigma_hat == pytest.approx(0.5593, abs=1e-3)


def test_holdout_noise_model_uses_published_sem_scoring():
    """Task A: the holdout score must use the published SEM, scoring -7.66 rather than +0.577.

    For a single-point holdout, estimating sigma from 1 residual gave +0.577 (absorbing the miss).
    Using the published SEM (0.03) gives -7.66 (scoring the 4.5x SEM miss as negative log-likelihood).
    Both sigma_hat and stated_sem are reported on DatasetFit.
    """
    import numpy as np
    from circuitpharm.fitting.data import (
        BARBERIS2007_PAIRED_PULSE,
        DOSE_RESPONSE_BENCHMARK,
        HOLDOUT,
        JAHN1997_PEAK_CRC,
        PairedPulseDataset,
        DoseResponseDataset,
    )
    from circuitpharm.fitting.comparison import MODEL_SPECS
    from circuitpharm.fitting.likelihood import score_holdout, concentrated_log_likelihood

    # 1. Fitted Model B parameters against Jahn CRC
    spec = MODEL_SPECS["kinetic_jw95"]
    v = np.array([2.0187309337948527, 2.214634354318136, -3.0])
    m_fit = spec.factory(v)

    res = score_holdout(m_fit, HOLDOUT, caller="test")
    fit = res.per_dataset[0]
    assert fit.noise_model == "stated_sem"
    assert fit.stated_sem == 0.03
    assert fit.sigma_hat == pytest.approx(0.1358, abs=1e-3)
    # Crucial assertion: -7.66 rather than +0.577
    assert res.total_log_likelihood == pytest.approx(-7.66, abs=0.02)
    assert fit.log_likelihood == pytest.approx(-7.66, abs=0.02)
    assert res.n_estimated == 0  # holdout has n_params=0 and stated_sem has no estimated noise scale

    # 2. A dataset declaring "stated_sem" with no sem raises
    with pytest.raises(ValueError, match="no sem"):
        DoseResponseDataset(
            citation_label="no sem stated",
            preparation="p",
            concs_um=np.array([1.0]),
            mean_response=np.array([0.5]),
            sem=None,
            kind="synthetic",
            synthetic=True,
            origin="o",
            noise_model="stated_sem",
        )

    # 3. Existing CRC path remains on "estimated" with unchanged behavior
    assert JAHN1997_PEAK_CRC.noise_model == "estimated"
    assert DOSE_RESPONSE_BENCHMARK.noise_model == "estimated"
    res_crc = concentrated_log_likelihood(m_fit, [JAHN1997_PEAK_CRC], n_params=3)
    fit_crc = res_crc.per_dataset[0]
    assert fit_crc.noise_model == "estimated"
    assert fit_crc.sigma_hat == pytest.approx(np.sqrt(fit_crc.rss / fit_crc.n), rel=1e-12)
    assert res_crc.n_estimated == 4  # 3 params + 1 estimated noise scale

