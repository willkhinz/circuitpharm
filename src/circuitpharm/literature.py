"""Every number read out of the alpha1beta2gamma2 full texts, with its conditions attached.

WHY THIS MODULE IS IN THE PACKAGE AND NOT IN scripts/. It began life as
`scripts/literature_gamma2.py`, and `fitting.data` imported the Barberis deactivation
components from there so they would not be retyped. That made the LIBRARY depend on the
scripts directory, which is not part of the installed distribution -- so
`import circuitpharm` succeeded from the repository root (where pyproject's
`pythonpath = ["."]` happens to put scripts/ in reach) and failed from anywhere else with
`ModuleNotFoundError: No module named 'scripts'`. CI never caught it because CI runs pytest
from the root.

The dependency now runs the correct way: the data lives here, inside the package, and
`scripts/literature_gamma2.py` is a presentation layer that imports it and prints the derived
quantities. Nothing in `src/` imports from `scripts/`.

WHAT IS IN HERE. The egress block that stopped the previous session was specific to the cloud
container; on 2026-10-10 pmc.ncbi.nlm.nih.gov was reachable and all seven open-access targets
in `fitting.data.MISSING_DATASETS` were read in full text. This module holds what those texts
say -- value, error AND ITS DECLARED STATISTIC, n, units, temperature, preparation, the
subunit combination exactly as the paper states it, whether it came from text, table, figure
or abstract, and the DOI/PMID/URL actually read.

THREE OF THE QUANTITIES THE PROJECT USES ARE DERIVED HERE, NOT READ:

  * NO source gives a single "mean open time" for alpha1beta2gamma2L. Keramidas & Harrison
    2008 and Li et al. 2008 give exponential DECOMPOSITIONS (tau_i with areas a_i). The
    area-weighted mean sum(a_i * tau_i) is computed, and it is ours, not theirs.
  * Barberis et al. 2007 report both the three components of their deactivation fit AND the
    weighted tau, so recomputing the latter from the former is a read-back check on this
    module's transcription.
  * `_NOMINAL_DEFECT` = 0.750 was compared against four P_o values that turned out not to be
    the same quantity; `PO_CANDIDATES` records what each one conditions on.

THE FINDING THAT MATTERS MOST. Every observable this project wanted to anchor on is
protocol-dependent, and the sources do not declare the protocol in comparable terms. PEAK was
already known (knowledge/12-inference.md section 2.3). DEACTIVATION: Barberis measured the
SAME receptor in the SAME study at 52.5 ms (2 ms pulse) and 364 ms (3 s pulse) -- 6.9x from
the protocol alone, a published measurement rather than a model prediction. MEAN OPEN TIME:
1.42 to 7.25 ms across sources, and the spread is mostly DEFINITION, not disagreement.

See knowledge/14-literature.md. `scripts/literature_gamma2.py` prints all of it.
"""
from dataclasses import dataclass, field  # noqa: F401

from dataclasses import dataclass

# =======================================================================================
# SOURCES. `where` records what was actually read -- full text at PMC, or only the PubMed
# abstract. `verification` uses the vocabulary already in data/pharmacology.db.
# =======================================================================================


@dataclass(frozen=True)
class Source:
    key: str
    citation: str
    journal_ref: str
    doi: str
    pmid: str
    url: str
    year: int
    read: str                    # FULL_TEXT_PMC | ABSTRACT_PUBMED | NOT_READ
    receptor_verbatim: str       # the subunit combination EXACTLY as the paper states it
    species: str
    prep: str
    temp_c: str
    err_kind: str                # SEM | SD | UNSTATED -- as the paper declares it
    note: str = ""


SOURCES = {
    "keramidas2008": Source(
        key="keramidas2008",
        citation=("Keramidas A, Harrison NL. Agonist-dependent single channel current and "
                  "gating in alpha4beta2delta and alpha1beta2gamma2S GABA-A receptors."),
        journal_ref="J Gen Physiol 131(2):163-181",
        doi="10.1085/jgp.200709871", pmid="18227271",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC2213567/",
        year=2008, read="FULL_TEXT_PMC",
        receptor_verbatim="alpha1beta2gamma2S",
        species="rat (subunit source not restated in Methods)",
        prep="HEK293, excised outside-out patch, -70 and +70 mV",
        temp_c="21 +/- 1",
        err_kind="UNSTATED",
        note=("The paper states NO error convention. SEM is the journal's usual and the "
              "same group's 2010 paper declares SEM, but this one does not, so every "
              "error bar below is recorded UNSTATED rather than assumed."),
    ),
    "keramidas2010": Source(
        key="keramidas2010",
        citation=("Keramidas A, Harrison NL. The activation mechanism of alpha1beta2gamma2S "
                  "and alpha3beta3gamma2S GABA-A receptors."),
        journal_ref="J Gen Physiol 135(1):59-75",
        doi="10.1085/jgp.200910317", pmid="20038526",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC2806416/",
        year=2010, read="FULL_TEXT_PMC",
        receptor_verbatim="alpha1beta2gamma2S",
        species="rat (not restated)",
        prep="HEK293, excised outside-out patch, -70 mV",
        temp_c="21 +/- 1",
        err_kind="SEM",
        note="Declares: 'Data are stated as mean +/- SEM for (n) separate experiments'.",
    ),
    "li2008": Source(
        key="li2008",
        citation=("Li P, Reichert DE, Rodriguez AD, Manion BD, Evers AS, Eterovic VA, "
                  "Steinbach JH, Akk G. Mechanisms of potentiation of the mammalian GABA-A "
                  "receptor by the marine cembranoid eupalmerin acetate."),
        journal_ref="Br J Pharmacol 153(3):598-608",
        doi="10.1038/sj.bjp.0707597", pmid="18037908",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC2241790/",
        year=2008, read="FULL_TEXT_PMC",
        receptor_verbatim="alpha1beta2gamma2L",
        species="rat",
        prep="HEK293 (FLAG-tagged), CELL-ATTACHED single-channel patch",
        temp_c="room temperature (no value given)",
        err_kind="SD",
        note=("Table 1 declares mean +/- s.d. The GABA-alone control rows are NOT this "
              "paper's own data: it states they are 'from Li et al. (2007b)' = Li P et al., "
              "J Physiol 584:789-800, doi 10.1113/jphysiol.2007.142794. Cite that for the "
              "control open times, not this paper. Cell-attached, not outside-out."),
    ),
    "dixon2014": Source(
        key="dixon2014",
        citation=("Dixon C, Sah P, Lynch JW, Keramidas A. GABA-A receptor alpha and gamma "
                  "subunits shape synaptic currents via different mechanisms."),
        journal_ref="J Biol Chem 289(9):5399-5411",
        doi="10.1074/jbc.M113.514695", pmid="24425869",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC3937617/",
        year=2014, read="FULL_TEXT_PMC",
        receptor_verbatim="alpha1beta2gamma2L",
        species="human (alpha1 pCIS2, beta2 pcDNA3.1+, gamma2L pcDNA3.1+; 1:1:3)",
        prep="HEK293, outside-out patch / whole cell, -70 mV",
        temp_c="room temperature (no value given)",
        err_kind="UNSTATED",
        note=("MISCITED in MISSING_DATASETS as 289(8); the article is 289(9):5399-5411. "
              "One of only two gamma2L sources found."),
    ),
    "dixon2015": Source(
        key="dixon2015",
        citation=("Dixon CL, Harrison NL, Lynch JW, Keramidas A. Zolpidem and eszopiclone "
                  "prime alpha1beta2gamma2 GABA-A receptors for longer duration of activity."),
        journal_ref="Br J Pharmacol 172(14):3522-3536",
        doi="10.1111/bph.13142", pmid="25817537",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC4507157/",
        year=2015, read="FULL_TEXT_PMC",
        receptor_verbatim="alpha1beta2gamma2S",
        species="human (alpha1 pCIS2, beta2 pcDNA3.1+, gamma2S pcDNA3.1+; 1:1:3)",
        prep="HEK293, outside-out patch (macropatch/single channel), -70 mV",
        temp_c="22 +/- 1",
        err_kind="UNSTATED",
        note=("gamma2S, NOT gamma2L -- MISSING_DATASETS lists it as a second gamma2L "
              "deactivation source, which it is not. Its MEASURED control deactivation tau "
              "appears only in Figure 8A; the text states only the SIMULATED value."),
    ),
    "barberis2007": Source(
        key="barberis2007",
        citation=("Barberis A, Mozrzymas JW, Ortinski PI, Vicini S. Desensitization and "
                  "binding properties determine distinct alpha1beta2gamma2 and "
                  "alpha3beta2gamma2 GABA-A receptor-channel kinetic behavior."),
        journal_ref="Eur J Neurosci 25(9):2726-2740",
        doi="10.1111/j.1460-9568.2007.05530.x", pmid="17561840",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC1950087/",
        year=2007, read="FULL_TEXT_PMC",
        receptor_verbatim="alpha1beta2gamma2 (Methods: 'Rat alpha1, beta2 and gamma2S')",
        species="rat",
        prep=("HEK293, outside-out patch (>=100 uM GABA) pooled with small lifted whole "
              "cells (<=30 uM); ultrafast perfusion, 10-90% exchange 60-100 us; -70 mV "
              "macroscopic, -100 mV single channel"),
        temp_c="22-24",
        err_kind="SEM",
        note=("Title and abstract say 'alpha1beta2gamma2'; the Methods say gamma2S. Treat "
              "as gamma2S. MISCITED in MISSING_DATASETS as Eur J Neurosci 26(7); the "
              "article is 25(9):2726-2740."),
    ),
    "mortensen2012": Source(
        key="mortensen2012",
        citation=("Mortensen M, Patel B, Smart TG. GABA potency at GABA-A receptors found "
                  "in synaptic and extrasynaptic zones."),
        journal_ref="Front Cell Neurosci 6:1",
        doi="10.3389/fncel.2012.00001", pmid="22319471",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC3262152/",
        year=2012, read="FULL_TEXT_PMC",
        receptor_verbatim="alpha1beta2gamma2S",
        species="not stated for the authors' own HEK293 data",
        prep="HEK293, whole-cell voltage clamp; GABA delivered with 20-30 ms latency",
        temp_c="not stated",
        err_kind="SEM",
        note=("A REVIEW presenting the authors' own dataset. Reports pEC50, EC50 and I_max "
              "only -- it contains NO Hill coefficients. MISSING_DATASETS and HANDOFF.md "
              "describe it as an 'EC50/nH table'; the nH half does not exist. The slow "
              "20-30 ms application means this is not a fast-application peak measurement."),
    ),
    "jahn_a1b2g2_kinetics": Source(
        key="jahn_a1b2g2_kinetics",
        citation=("Jahn K, Hertle I, Bufler J, Adelsberger H, Pestel E, Zieglgaensberger W, "
                  "Dudel J, Franke C. Activation kinetics and single channel properties of "
                  "recombinant alpha1beta2gamma2L GABA-A receptor channels."),
        journal_ref="NeuroReport 8(16):3443-3446",
        doi="10.1097/00001756-199711100-00006", pmid="9427304",
        url="https://pubmed.ncbi.nlm.nih.gov/9427304/",
        year=1997, read="ABSTRACT_PUBMED",
        receptor_verbatim="alpha1beta2gamma2L",
        species="not stated in abstract",
        prep="HEK293, patch clamp with ultra-fast solution exchange",
        temp_c="not stated in abstract",
        err_kind="UNSTATED (abstract uses +/- without naming the statistic)",
        note=("No PMC copy; paywalled. The ABSTRACT contains NO mean open time -- it "
              "reports BURST duration 10.3 +/- 3.0 ms. Whether the full text states a mean "
              "open time is STILL UNKNOWN and cannot be settled from the abstract. PubMed "
              "renders the subunits correctly as alpha1beta2gamma2L; an Ovid page garbles "
              "them."),
    ),
    "hamill1997": Source(
        key="hamill1997",
        citation=("Hamill OP. How many transmitter binding steps are involved in opening "
                  "fast receptor-gated channels?"),
        journal_ref="NeuroReport 8(16):iv",
        doi="", pmid="9480006",
        url="https://pubmed.ncbi.nlm.nih.gov/9480006/",
        year=1997, read="NOT_READ",
        receptor_verbatim="n/a (commentary)",
        species="n/a", prep="n/a", temp_c="n/a", err_kind="n/a",
        note=("PubMed confirms it is a Comment on Jahn 1997 (NeuroReport 8(16):3443-6) and "
              "that it occupies page 'iv' -- a SINGLE page in the issue's front matter. No "
              "abstract, no DOI, no PMC copy. The TITLE is the only content available, and "
              "it is itself informative: it questions the number of binding steps, which "
              "is the inference Jahn drew from the Hill slope. CONTENT STILL UNREAD."),
    ),
}

# =======================================================================================
# OPEN TIME. Exponential decompositions, verbatim. `areas` are the fractions the source
# reports; the area-weighted mean is DERIVED below and is not a number any source states.
# =======================================================================================


@dataclass(frozen=True)
class OpenTimeFit:
    label: str
    source_key: str
    receptor: str
    agonist_um: float
    conditioning: str        # what the measurement is conditional on -- the crux
    taus_ms: tuple
    tau_errs: tuple
    areas: tuple
    n_patches: str
    err_kind: str


OPEN_TIME_FITS = (
    OpenTimeFit(
        label="Keramidas 2008 M-Mode", source_key="keramidas2008",
        receptor="alpha1beta2gamma2S", agonist_um=10_000.0,
        conditioning="WITHIN bursts, medium-P_o gating mode (11 of ~15 patches)",
        taus_ms=(0.49, 2.58, 5.17), tau_errs=(0.02, 0.09, 0.38),
        areas=(0.26, 0.49, 0.25), n_patches="10", err_kind="UNSTATED",
    ),
    OpenTimeFit(
        label="Keramidas 2008 H-Mode", source_key="keramidas2008",
        receptor="alpha1beta2gamma2S", agonist_um=10_000.0,
        conditioning="WITHIN bursts, high-P_o gating mode (5 patches)",
        taus_ms=(0.59, 4.22, 13.2), tau_errs=(0.11, 0.73, 2.6),
        areas=(0.18, 0.41, 0.41), n_patches="5", err_kind="UNSTATED",
    ),
    OpenTimeFit(
        label="Li 2008 / Li 2007b control", source_key="li2008",
        receptor="alpha1beta2gamma2L", agonist_um=50.0,
        conditioning="WITHIN clusters (intracluster), 50 uM GABA is SUB-SATURATING",
        taus_ms=(0.28, 3.0, 7.3), tau_errs=(0.05, 0.7, 3.2),
        areas=(0.22, 0.65, 0.13), n_patches="4", err_kind="SD",
    ),
)

#: A mean open time STATED as such by a source, rather than derived. Barberis is the only
#: one, and it is a non-stationary measurement over a 4 s window after a brief pulse, so it
#: mixes the early high-occupancy openings with the late brief singly-bound ones.
STATED_MEAN_OPEN_TIME = {
    "barberis2007": dict(
        receptor="alpha1beta2gamma2 (gamma2S per Methods)",
        value_ms=1.42, err_ms=0.05, err_kind="SEM", n="6",
        conditioning=("mean open time over a 4 s window AFTER a 2 ms pulse of saturating "
                      "(10 mM) GABA, in patches with few channels -- NON-STATIONARY"),
        note=("The same paper's STEADY-STATE intraburst open-time distribution (two "
              "exponentials, with its own mean open time) is reported only in Figure 6C "
              "and so could not be read. That figure, not this number, is what would be "
              "comparable to the Keramidas and Li decompositions."),
    ),
}

# =======================================================================================
# DEACTIVATION. The protocol column is the finding: two of these rows are the same
# receptor in the same study, and they differ 6.9-fold.
# =======================================================================================


@dataclass(frozen=True)
class Deactivation:
    label: str
    source_key: str
    receptor: str
    pulse: str
    conc_mm: float
    tau_w_ms: float
    err_ms: float
    err_kind: str
    n: str
    components: tuple = ()      # (taus, areas) where the source states them
    component_errors: tuple = ()  # (tau_errs, area_errs) where the source states them (promoted for weighting iv)
    where: str = "full text"
    note: str = ""


DEACTIVATIONS = (
    Deactivation(
        label="Dixon 2014", source_key="dixon2014", receptor="alpha1beta2gamma2L",
        pulse="<=1 ms", conc_mm=3.0, tau_w_ms=5.9, err_ms=0.5, err_kind="UNSTATED", n="10",
        note=("VERBATIM: 'The weighted deactivation time constants for alpha1beta2gamma2L "
              "and alpha1beta2gamma1 GABA-A Rs were 5.9 +/- 0.5 (n = 10) and 9.1 +/- 0.9 "
              "ms (n = 6), respectively.' The [LEAD] is CONFIRMED. The individual "
              "components behind this weighted tau are NOT stated anywhere in the full "
              "text -- the only 'two exponential' fit mentioned is of the SIMULATED "
              "current, not the measured one."),
    ),
    Deactivation(
        label="Barberis 2007, brief pulse", source_key="barberis2007",
        receptor="alpha1beta2gamma2 (gamma2S)",
        pulse="2 ms", conc_mm=10.0, tau_w_ms=52.5, err_ms=2.9, err_kind="SEM", n="6",
        components=((2.8, 33.4, 221.35), (0.57, 0.23, 0.20)),
        component_errors=((0.3, 4.6, 14.9), (0.04, 0.02, 0.03)),
        note=("Triple exponential; component errors are 0.3, 4.6, 14.9 ms and area errors "
              "0.04, 0.02, 0.03. Fully specified -- the only observable found that can be "
              "reconstructed from parameters the source itself fitted, with no "
              "extrapolation beyond the measured window."),
    ),
    Deactivation(
        label="Barberis 2007, long pulse", source_key="barberis2007",
        receptor="alpha1beta2gamma2 (gamma2S)",
        pulse="3 s", conc_mm=10.0, tau_w_ms=364.0, err_ms=40.0, err_kind="SEM", n="7",
        note=("Same receptor, same patches, same study as the row above. VERBATIM: 'we "
              "studied the kinetics of the current relaxation after a long (3 s) pulse. "
              "The deactivation process was markedly slower ... with weighted time "
              "constants of 743.3 +/- 73.7 and 364 +/- 40 ms, respectively (n = 7)' "
              "-- 743.3 is alpha3beta2gamma2, 364 is alpha1beta2gamma2."),
    ),
    Deactivation(
        label="Dixon 2015 (SIMULATED, not measured)", source_key="dixon2015",
        receptor="alpha1beta2gamma2S",
        pulse="<=1 ms", conc_mm=5.0, tau_w_ms=9.0, err_ms=float("nan"),
        err_kind="n/a (single simulated value)", n="n/a",
        where="full text, but of a SIMULATION",
        note=("This is the ensemble current simulated from their fitted mechanism, with "
              "two exponential components, NOT a measurement. The measured control value "
              "is in Figure 8A and was not readable. Do not use as data."),
    ),
)

# =======================================================================================
# P_o CANDIDATES. parameters._NOMINAL_DEFECT = 0.750. The four values that surfaced as
# [LEAD]s are different quantities; `conditions_on` is why they are not interchangeable.
# =======================================================================================

NOMINAL_DEFECT = 0.750

PO_CANDIDATES = (
    dict(value=0.69, err=0.02, n="11", source_key="keramidas2008",
         receptor="alpha1beta2gamma2S", agonist="10 mM GABA",
         quantity="INTRABURST P_o, M-mode",
         conditions_on="being inside a burst AND in the medium gating mode",
         note=("This is the 0.69 that MISSING_DATASETS/HANDOFF.md attribute to "
               "'nonstationary variance analysis'. IT IS NOT. It is Table II's "
               "intraburst P_o for the M-mode. The conflation matters because the "
               "handoff concluded '[t]he one genuinely macroscopic estimate is the "
               "0.69' and built the Task-4 recommendation on that.")),
    dict(value=0.87, err=0.02, n="5", source_key="keramidas2008",
         receptor="alpha1beta2gamma2S", agonist="10 mM GABA",
         quantity="INTRABURST P_o, H-mode",
         conditions_on="being inside a burst AND in the high gating mode",
         note="Source of the '~0.9' lead. Keramidas 2010 restates the pair as '~0.7 and ~0.9'."),
    dict(value=0.56, err=float("nan"), n="not stated", source_key="keramidas2008",
         receptor="alpha1beta2gamma2S", agonist="1-2 ms exposure",
         quantity="P_o MACROPATCH by nonstationary fluctuation analysis",
         conditions_on="nothing -- this IS the population/macroscopic quantity",
         note=("THE actual nonstationary-fluctuation-analysis number, and the only "
               "macroscopic P_o found. VERBATIM: 'we estimated the channel open "
               "probability (P_O macropatch) in response to 1-2-ms exposure to agonist "
               "using nonstationary fluctuation analysis. The P_O macropatch values we "
               "obtained were ... for alpha1beta2gamma2S GABA-A Rs, GABA-0.56 and "
               "THIP-0.42 (UNPUBLISHED DATA).' Marked unpublished by its own authors, "
               "with no error bar and no n, and measured at a 1-2 ms application -- the "
               "regime where this project's own PEAK EC50 is ~299 uM, not 39 uM.")),
    dict(value=0.81, err=0.01, n="7", source_key="keramidas2010",
         receptor="alpha1beta2gamma2S", agonist="5 mM GABA",
         quantity="INTRABURST P_o, modes pooled",
         conditions_on="being inside a burst",
         note=("Concentration series: 0.81 (5 mM), 0.61 +/- 0.04 (n=11, 200 uM), 0.33 "
               "+/- 0.06 (n=4, 20 uM), 0.15 +/- 0.01 (n=5, 2 uM). Methods are explicit "
               "that P_o here means intraburst.")),
    dict(value=0.56, err=0.04, n="3 to 8", source_key="dixon2014",
         receptor="alpha1beta2gamma2L", agonist="3 mM GABA",
         quantity="INTRABURST P_o",
         conditions_on="being inside a burst",
         note=("The ONLY P_o found for gamma2L, the project's actual receptor. Table 1; "
               "0.37 +/- 0.01 at 2 uM GABA. n given only as 'averages from 3 to 8 "
               "patches' for the whole table.")),
    dict(value=0.42, err=0.04, n="4", source_key="li2008",
         receptor="alpha1beta2gamma2L", agonist="50 uM GABA",
         quantity="CLUSTER open probability",
         conditions_on="being inside a cluster; 50 uM is sub-saturating",
         note=("Defined in the paper as mean open time / (mean open time + mean closed "
               "time). Control data from Li et al. 2007b.")),
)

#: The 0.8 'intraburst P_o' lead in HANDOFF.md section 5.3 was NOT found as a stated value
#: for alpha1beta2gamma2 in any of the seven full texts. The nearest readings are
#: Keramidas 2010's 0.81 +/- 0.01 at 5 mM (intraburst, modes pooled) and Keramidas 2008's
#: H-mode 0.87. Recorded as UNLOCATED rather than silently matched to either.
PO_LEAD_UNLOCATED = "0.8 as an intraburst P_o -- no stated source found in the seven texts"

# =======================================================================================
# OTHER VERIFIED READINGS, kept because they bear on datasets the project already holds.
# =======================================================================================

OTHER = (
    dict(quantity="EC50, peak current", value="11.6 +/- 0.9 uM", n="not stated",
         source_key="jahn_a1b2g2_kinetics", receptor="alpha1beta2gamma2L",
         prep="HEK293, ultra-fast exchange", where="ABSTRACT",
         note="Saturates with 3 mM GABA. The anchor JAHN1997_PEAK_CRC is built on."),
    dict(quantity="Hill slope over 0.001-0.01 mM", value="2.2 +/- 0.4", n="not stated",
         source_key="jahn_a1b2g2_kinetics", receptor="alpha1beta2gamma2L",
         prep="HEK293, ultra-fast exchange", where="ABSTRACT",
         note=("Re-confirmed verbatim from PubMed. 0.001-0.01 mM = 1-10 uM = "
               "0.086-0.862 x EC50: a LOCAL rising-phase slope, as PR #2 established.")),
    dict(quantity="burst duration", value="10.3 +/- 3.0 ms", n="not stated",
         source_key="jahn_a1b2g2_kinetics", receptor="alpha1beta2gamma2L",
         prep="HEK293", where="ABSTRACT",
         note="NOT a mean open time. The abstract gives no mean open time."),
    dict(quantity="single channel slope conductance", value="~29 pS", n="not stated",
         source_key="jahn_a1b2g2_kinetics", receptor="alpha1beta2gamma2L", prep="HEK293",
         where="ABSTRACT", note=">95% of current. Compare Barberis 27.7 +/- 0.8 pS."),
    dict(quantity="EC50", value="6.6 uM (pEC50 5.180 +/- 0.0593)", n="34",
         source_key="mortensen2012", receptor="alpha1beta2gamma2S",
         prep="HEK293 whole-cell, 20-30 ms application latency", where="Table 1",
         note=("I_max 2230 +/- 193 pA (n=18). NOT comparable to Jahn's 11.6 uM: slow "
               "application, whole-cell, gamma2S. No Hill coefficient is reported.")),
    dict(quantity="desensitisation onset, fast component",
         value="tau_1 = 2.9 +/- 0.1 ms, A_1 = 0.56 +/- 0.025",
         n="7 in text, 3 in the Fig. 4 caption -- THE SOURCE IS INCONSISTENT",
         source_key="barberis2007", receptor="alpha1beta2gamma2 (gamma2S)",
         prep="HEK293, 3 s pulse of 10 mM GABA", where="full text + Fig. 4 caption",
         note=("THE HOLDOUT OBSERVABLE. tau_2 and tau_3 exist but are reported only in "
               "Figure 4C, so the onset time course CANNOT be fully reconstructed. "
               "Steady-state weight 0.076 +/- 0.013 (n=7); steady-state:peak measured at "
               "200 ms is 0.21 +/- 0.02. 'in ~3 ms the current is reduced by more than "
               "one half'.")),
    dict(quantity="paired-pulse recovery, 100 ms gap", value="0.33 +/- 0.03",
         n="8 (Fig. 5 caption)", source_key="barberis2007",
         receptor="alpha1beta2gamma2 (gamma2S)",
         prep="HEK293, 2 ms pulses of saturating GABA", where="Table 1",
         note="A second candidate holdout, fully stated as a scalar."),
    dict(quantity="10-90% rise time at 10 mM GABA", value="0.29 +/- 0.02 ms", n="not stated",
         source_key="barberis2007", receptor="alpha1beta2gamma2 (gamma2S)",
         prep="HEK293, ultrafast exchange", where="Table 1", note=""),
    dict(quantity="single channel conductance", value="27.7 +/- 0.8 pS", n="5",
         source_key="barberis2007", receptor="alpha1beta2gamma2 (gamma2S)",
         prep="HEK293, outside-out", where="full text", note=""),
    dict(quantity="mean burst length", value="148 +/- 16 ms (3 mM), 23 +/- 2 ms (2 uM)",
         n="3 to 8", source_key="dixon2014", receptor="alpha1beta2gamma2L",
         prep="HEK293, outside-out", where="Table 1",
         note="The 2 uM figure counts only bursts with >=2 events."),
)



def _weighted(taus, areas):
    """Area-weighted mean of an exponential decomposition. DERIVED, never published."""
    return sum(a * t for a, t in zip(areas, taus))
