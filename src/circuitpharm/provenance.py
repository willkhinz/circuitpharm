"""Where every load-bearing number came from, enforced in code.

WHY THIS MODULE EXISTS. `results.Tier` makes a model OUTPUT carry its reliability. This does
the same for model INPUTS, and it was added because an audit found the harder problem:
the numbers driving the project's one VALIDATED result were not unverified, they were
UNSOURCED. `subtypes.REGIONS`, `subtypes.EXTRASYN` and `subtypes.SUBJECTIVE_WEIGHT` between
them hold 20 numbers that `scripts/ranking_robustness.py` scores on, and the module that
holds them contains zero citations. (It appears to contain five; all five are compound names
that collide with source keys.)

The knowledge base could not have supplied them either: its `subunit_expression` table holds
eight rows and every one is qualitative -- "predominant", "present", "high", "enriched". There
were no numbers to read.

So each number now carries a `Basis` saying what kind of thing it is. The point is not to
make the unsourced ones look sourced; it is to make the count visible and to stop a new
number being added without one. `tests/test_provenance.py` fails if any parameter is missing
a record.

A BASIS IS A CLAIM ABOUT DERIVATION, NEVER ABOUT CORRECTNESS. An UNSOURCED number may well
be right; a SOURCED one may be misread. See knowledge/06-source-provenance.md for the
per-source reading, including `a5_dist`, which resolves perfectly by DOI and does NOT support
the numbers attributed to it.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Basis(str, Enum):
    """How a number came to have the value it has."""

    #: read as a number from a named source
    QUANTITATIVE = "QUANTITATIVE"
    #: a named source makes a QUALITATIVE statement ("predominant", "minimal") that was
    #: converted to a number here. The conversion rule is the provenance, and it is ours.
    FROM_QUALITATIVE = "FROM_QUALITATIVE"
    #: set to make a different, anchored quantity come out right. Legitimate, but it means
    #: the number carries that anchor's assumptions and cannot be quoted independently.
    FITTED = "FITTED"
    #: a structural choice, not an empirical quantity (e.g. a reference value fixed at 1.0)
    CONVENTION = "CONVENTION"
    #: no source identified. Stated so, rather than dressed up.
    UNSOURCED = "UNSOURCED"
    #: explicitly a guess in the original code
    GUESS = "GUESS"


@dataclass(frozen=True)
class Record:
    basis: Basis
    source_key: str = ""      # key into the `sources` table, where one exists
    note: str = ""            # the derivation, or why none exists

    @property
    def is_sourced(self) -> bool:
        return self.basis in (Basis.QUANTITATIVE, Basis.FROM_QUALITATIVE) and bool(
            self.source_key)


# ---------------------------------------------------------------------------------------
# REGIONS -- per-region subunit composition. These 15 numbers, with EXTRASYN and
# SUBJECTIVE_WEIGHT, are the entire input to the selectivity ranking.
#
# `a5_dist` is NOT cited for any of them, deliberately. It is the obvious candidate and it
# was assessed by hand: a 1988 study using a single probe for "the alpha subunit", reporting
# TOTAL alpha mRNA by region (cerebellum > thalamus = cortex = hippocampus >> pons =
# striatum = medulla). It does not resolve alpha subtypes, so it cannot support a
# per-subtype fraction -- and it is about regional LEVEL, where REGIONS encodes regional
# COMPOSITION (its rows sum to 1 by construction). Citing it would be worse than citing
# nothing.
# ---------------------------------------------------------------------------------------
_PBC_QUAL = ("Liu & Wong-Riley 2004 report developmental expression of alpha1/2/3 in rat "
             "preBotC. The paper is paywalled (HTTP 403) and has not been read; its title "
             "and DOI are confirmed. It covers alpha1/2/3 only -- not a5, delta/a4 or "
             "epsilon, which this row also assigns. The specific fraction is ours.")

REGIONS_PROV = {
    ("prebotc", "a1"):   Record(Basis.FROM_QUALITATIVE, "pbc_alpha", _PBC_QUAL),
    ("prebotc", "a23"):  Record(Basis.FROM_QUALITATIVE, "pbc_alpha", _PBC_QUAL),
    ("prebotc", "a5"):   Record(Basis.UNSOURCED, "", (
        "0.02 encodes 'a5 minimal in medulla'. No source in the knowledge base states a "
        "preBotC a5 fraction. This is the single most consequential number in the package: "
        "an a5-selective compound's modelled respiratory burden is roughly proportional to "
        "it, so it sets the whole safety margin.")),
    ("prebotc", "d_a4"): Record(Basis.FROM_QUALITATIVE, "pbc_delta", (
        "pbc_delta and pbc_a4 establish that delta- and alpha4-containing receptors are "
        "PRESENT in respiratory networks and that alpha4 loss causes respiratory "
        "dysfunction. Neither gives a fraction; 0.15 is ours.")),
    ("prebotc", "eps"):  Record(Basis.FROM_QUALITATIVE, "pbc_eps", (
        "pbc_eps reports epsilon ENRICHED on ventral respiratory column neurons. "
        "'Enriched' is not 0.08; the number is ours.")),

    ("spinal", "a1"):    Record(Basis.UNSOURCED, "", "no spinal subunit source recorded"),
    ("spinal", "a23"):   Record(Basis.UNSOURCED, "", "no spinal subunit source recorded"),
    ("spinal", "a5"):    Record(Basis.UNSOURCED, "", "no spinal subunit source recorded"),
    ("spinal", "d_a4"):  Record(Basis.UNSOURCED, "", "no spinal subunit source recorded"),
    ("spinal", "eps"):   Record(Basis.UNSOURCED, "", "no spinal subunit source recorded"),

    ("forebrain", "a1"):   Record(Basis.UNSOURCED, "", "no forebrain subunit source recorded"),
    ("forebrain", "a23"):  Record(Basis.UNSOURCED, "", "no forebrain subunit source recorded"),
    ("forebrain", "a5"):   Record(Basis.FROM_QUALITATIVE, "", (
        "0.30 encodes the in-code comment '>25% of CA1/CA3 neurons express a5, primarily "
        "extrasynaptic'. That statement has NO source key attached anywhere in the "
        "repository, so the figure it rests on is untraceable. It is also a statement about "
        "the fraction of NEURONS expressing a5, which is not the fraction of RECEPTORS that "
        "are a5 -- the quantity REGIONS actually holds.")),
    ("forebrain", "d_a4"): Record(Basis.UNSOURCED, "", "no forebrain subunit source recorded"),
    ("forebrain", "eps"):  Record(Basis.UNSOURCED, "", "no forebrain subunit source recorded"),
}

# ---------------------------------------------------------------------------------------
# EXTRASYN -- synaptic vs extrasynaptic split. Decides how much headroom a PAM has in each
# pool (gabaa_kinetics: ~1.1x synaptic against ~211x extrasynaptic), so these numbers carry
# the project's central safety argument.
# ---------------------------------------------------------------------------------------
EXTRASYN_PROV = {
    "a1":   Record(Basis.UNSOURCED, "", "0.15; no source recorded"),
    "a23":  Record(Basis.UNSOURCED, "", "0.20; no source recorded"),
    # DOWNGRADED from FROM_QUALITATIVE after a claim-support read. A manuscript draft
    # attributed this to Kasugai et al. 2010 (Eur J Neurosci 32:1868-1888), which is the
    # obvious candidate and the right kind of study -- quantitative freeze-fracture replica
    # immunogold, synaptic vs extrasynaptic pools, hippocampal CA1 pyramidal cells. Its
    # abstract says it measured **a1, a2 and b3**. It does not measure a5 at all, and its
    # quantitative result runs the other way: synaptic labelling density exceeded
    # extrasynaptic by 78-132x (a1), 94x (a2) and 79x (b3). So it supports a LOW
    # extrasynaptic fraction for the subunits it did measure, and is silent on this one.
    #
    # This is the second source to resolve perfectly and fail to support the number it was
    # attached to; `a5_dist` was the first. The pattern is specific enough to name: a source
    # whose title matches the claim, in the right journal, by the right group, measuring a
    # neighbouring quantity.
    "a5":   Record(Basis.UNSOURCED, "", (
        "0.80 encodes the widely repeated statement that hippocampal a5 is predominantly "
        "extrasynaptic. The direction is well supported in review literature, but no primary "
        "source in this repository supports it, and the best candidate (Kasugai et al. 2010) "
        "measures a1/a2/b3 and not a5. Load-bearing: this is the whole extrasynaptic-headroom "
        "argument, so it is recorded as UNSOURCED rather than FROM_QUALITATIVE.")),
    "d_a4": Record(Basis.FROM_QUALITATIVE, "pbc_delta", (
        "1.00: delta-containing receptors are exclusively extrasynaptic. This one is a "
        "structural fact rather than a measured fraction, and is the most defensible entry "
        "in the table.")),
    "eps":  Record(Basis.GUESS, "", (
        "0.50, labelled '# eps = GUESS' in subtypes.py itself. The robustness analysis "
        "treats it as one, drawing uniform(0,1) rather than a concentrated Beta.")),
}

# ---------------------------------------------------------------------------------------
# SUBJECTIVE_WEIGHT -- which subtypes carry ethanol's discriminative stimulus.
# ---------------------------------------------------------------------------------------
SUBJECTIVE_PROV = {
    "a5":   Record(Basis.FROM_QUALITATIVE, "a5_disc", (
        "a5_disc (Contribution of a1GABA-A and a5GABA-A receptor subtypes to the "
        "discriminative stimulus effects of ethanol in squirrel monkeys, 2005) reports that "
        "a5 agonists QH-ii-066 and panadiplon mimic ethanol's stimulus and the a5 inverse "
        "agonist L-655,708 blocks it. Supports the DIRECTION. The weight 1.0 is ours.")),
    "a23":  Record(Basis.UNSOURCED, "", (
        "1.0, equal to a5, with no source. gabox_disc bears on a2/a3 discriminative effects "
        "but has not been read.")),
    "a1":   Record(Basis.UNSOURCED, "", (
        "0.0 -- and this one is actively questionable. The single source bearing on it, "
        "a5_disc, is titled 'Contribution of a1GABA-A AND a5GABA-A receptor subtypes to the "
        "discriminative stimulus effects of ethanol'. Zeroing a1 on the strength of a paper "
        "whose title names an a1 contribution needs that paper read before it stands.")),
    "d_a4": Record(Basis.UNSOURCED, "", "0.0; no source recorded"),
    "eps":  Record(Basis.UNSOURCED, "", "0.0; no source recorded"),
}

# ---------------------------------------------------------------------------------------
# Priors used by the robustness analysis. These are not pharmacological quantities but they
# decide what "robust across 20,000 draws" means, so they need provenance too -- arguably
# more than the point values, since the draws are what the headline claim rests on.
# ---------------------------------------------------------------------------------------
PRIOR_PROV = {
    "KAPPA": Record(Basis.UNSOURCED, "", (
        "15.0, the Dirichlet concentration on REGIONS. Sets how far the draws explore from "
        "an unsourced centre. At the nominal preBotC a5 of 0.02 it gives alpha=0.3, which "
        "puts 38% of draws below a tenth of nominal.")),
    "EXTRASYN_CONC": Record(Basis.UNSOURCED, "", (
        "8.0, the Beta concentration on EXTRASYN. Unsourced width around unsourced "
        "centres.")),
    "GAIN_RATIO": Record(Basis.FITTED, "", (
        "(3.0, 7.5) tonic:phasic, taken from this project's own mechanism-mix sweep "
        "(session 7d), not from literature. Internally derived, so it cannot corroborate "
        "a result the same model produces.")),
}

# ---------------------------------------------------------------------------------------
# MODEL PARAMETERS for the three competing receptor models (models/). Added with P0-13.
#
# All of these are PROVISIONAL dataclass defaults, not calibrations, and the register entry
# records why that matters: both entry points to the three-tiered headroom hierarchy used
# to default to them.
# ---------------------------------------------------------------------------------------
_CHIMERA = (
    "KineticAllosteryModel's defaults are a CHIMERA of two fits: beta/alpha are the "
    "manuscript's current config-pulse fit (0.559023 / 0.107891) while kon/koff are "
    "approximately the superseded 1000 uM / 0.30 ms fit. The resulting K_d = 29.62 uM "
    "appears in no fit, no commit and no document, and the set reproduces neither anchor "
    "(measured EC50 6.34 uM against 20, P_o,max 0.8382 against 0.750).")

MODEL_PARAM_PROV = {
    ("kinetic_jw95", "kon"):   Record(Basis.UNSOURCED, "", _CHIMERA),
    ("kinetic_jw95", "koff"):  Record(Basis.UNSOURCED, "", _CHIMERA),
    ("kinetic_jw95", "beta"):  Record(Basis.UNSOURCED, "", _CHIMERA),
    ("kinetic_jw95", "alpha"): Record(Basis.UNSOURCED, "", _CHIMERA),
    ("kinetic_jw95", "d"):     Record(Basis.UNSOURCED, "", (
        "0.050 ms^-1, carried over from gabaa_kinetics.DEFAULT_RATES, where it is "
        "explicitly NOT fitted -- the three macroscopic anchors do not constrain "
        "desensitisation. Every asymptotic headroom number is conditional on it.")),
    ("kinetic_jw95", "r"):     Record(Basis.UNSOURCED, "", (
        "0.0020 ms^-1, same status as d. d/r = 25.0 sets the asymptote directly.")),
    ("extended_desens", "all"): Record(Basis.UNSOURCED, "", (
        "kon 0.012, koff 0.35, beta 0.60, alpha 0.10 and four desensitisation rates, all "
        "round numbers with no recorded origin. Measured EC50 5.18 uM against the 20 uM "
        "anchor and P_o,max 0.8571 against 0.750.")),
    ("extended_desens", "pam_desens_factor"): Record(Basis.GUESS, "", (
        "1.0, with the modulation entering as d_fast/(1 + 0.1*(pf-1)*f). Neither the 0.1 "
        "coefficient nor the linear form has a source; a thermodynamic cycle over the "
        "R/O/D loop (which the technical spec claims) would constrain it by detailed "
        "balance instead.")),
    ("operational", "ec50_um"): Record(Basis.UNSOURCED, "",
                                       "25.0 uM, a declared round number."),
    ("operational", "hill_n"):  Record(Basis.UNSOURCED, "", "1.4, no source."),
    ("operational", "po_max"):  Record(Basis.FITTED, "", (
        "0.75 -- this is gabaa_kinetics.FIT_TARGETS['po_max'], i.e. the project's own fit "
        "TARGET rather than something fitted to data. Circular if used as an anchor.")),
    ("operational", "s_max"):   Record(Basis.UNSOURCED, "", (
        "2.50, the classical BZ-site intrinsic efficacy range. Carried as a ceiling, not "
        "measured here.")),
    ("gabaa_kinetics", "FIT_FIXED_ALPHA"): Record(Basis.CONVENTION, "", (
        "0.30 ms^-1, the channel closing rate HELD FIXED in gabaa_kinetics.fit_scheme so "
        "the fit is well-posed. 1/alpha is the mean open time, 3.33 ms, which is inside "
        "the 1-5 ms range reported for alpha1beta2gamma2 single channels. "
        "WHY: four parameters against three anchors left the optimiser's landing point "
        "dependent on the environment -- measured k_on 0.0088044 / 0.0109865 / 0.0110210 "
        "and asymptotic headroom 151.1x / 181.8x / 182.2x across scipy 1.14.1 / 1.11.4 / "
        "1.17.1, every one reproducing all three anchors, against the manuscript's own "
        "0.0146842 and 210.7x which none of them reproduces. Holding alpha makes the "
        "system square (3 unknowns, 3 anchors); verified identical to 6 significant "
        "figures across all three environments and from a distant starting guess. "
        "A minimum-norm regularisation was tried first and rejected on measurement: it "
        "narrowed the spread to ~0.8% without actually selecting the minimum-norm fit. "
        "This is a CONVENTION about a REAL single-channel observable, so roadmap P1-2 can "
        "replace it with a digitised mean open time and the fit becomes data-determined; "
        "it makes the rates REPRODUCIBLE, not identifiable.")),
    ("operational", "tau_deact_ms"): Record(Basis.FITTED, "", (
        "15.0 ms, again the project's own tau fit target -- and the value the old decay "
        "estimator returned on failure, which is what made a failed fit look perfect "
        "(P0-4).")),
}

# ---------------------------------------------------------------------------------------
# DATASETS in fitting/data.py. Both currently synthetic (P0-6).
# ---------------------------------------------------------------------------------------
DATASET_PROV = {
    "dose_response_peak": Record(Basis.UNSOURCED, "", (
        "SYNTHETIC. Hand-written response and SEM arrays shaped to resemble an "
        "alpha1beta2gamma2 peak concentration-response. Its plateau is exactly 0.750, "
        "which is the model's own po_max fit target, so a fit against it is partly "
        "circular; the SEM values have no origin. Mortensen et al. 2012 and Sigel & "
        "Steinmann 2012 are recorded as MOTIVATION only -- neither was digitised.")),
    "jahn1997_peak_crc": Record(Basis.QUANTITATIVE, "jahn_a1b2g2_kinetics", (
        "PARAMETRIC, and the project's FIRST sourced kinetic dataset. Jahn et al. 1997 "
        "(NeuroReport 8(16):3443-6, PMID 9427304) measured alpha1beta2gamma2L in HEK293 "
        "with ultra-fast solution exchange and published EC50 = 11.6 +/- 0.9 uM and a "
        "Hill-type slope of 2.2 +/- 0.4 over 1-10 uM, saturating at 3 mM. The nine points "
        "in the dataset are generated from those two parameters, so the PARAMETERS are the "
        "measurement and the points are the paper's own model of its data -- recorded "
        "QUANTITATIVE on that basis, with kind='parametric' and the fit tainted as for a "
        "synthetic trace. The full text could not be read (pubmed/PMC/EuropePMC are "
        "egress-blocked here), so the figure was never digitised. "
        "ONE DISCREPANCY IT SURFACES, unresolved: FIT_TARGETS['ec50_um'] is 20.0 uM "
        "against this 11.6. A SECOND, the scheme's 1.294 peak Hill slope against this 2.2, "
        "was recorded here as structural and is WITHDRAWN: the 2.2 is the slope 'between "
        "0.001 and 0.01 mM GABA', a LOCAL measurement on the rising phase, and 1.294 was a "
        "whole-curve regression. Measured alike the scheme gives 1.59-1.71, inside the "
        "published error bar; published whole-curve fits for this receptor are 1.3-1.6. "
        "The authors' inference of at-least-three binding sites does not follow either -- "
        "a two-site scheme reaches ~2.0 on that window and cryo-EM shows two sites. "
        "See knowledge/12-inference.md section 2.")),
    "deactivation_charge": Record(Basis.UNSOURCED, "", (
        "SYNTHETIC. 0.70*exp(-t/15) + 0.30*exp(-t/70). It previously carried "
        "citation='Haas & Macdonald 1999 / Jones & Westbrook 1995'; "
        "knowledge/06-source-provenance.md and 07-paper-review.md record Haas & Macdonald "
        "1999 as measuring 76.1 ms for this quantity, so the dominant 15 ms component is "
        "NOT that paper's number. The reference survives under `motivated_by`.")),
}

#: WAVEFORM SHAPE PARAMETERS (roadmap P6-1). Every number that shapes an agonist
#: time course, with what it rests on -- because a CHARGE quantity is a property of the
#: waveform as much as of the receptor, and a spillover amplitude nobody sourced propagates
#: into every charge ratio downstream of it without appearing in any parameter table.
#:
#: NOT ONE OF THESE IS QUANTITATIVE. The two that matter most are the spillover amplitude
#: and its tau: they set how much agonist an extrasynaptic receptor actually sees, which is
#: the whole mechanism of the tonic arm.
WAVEFORM_PROV = {
    ("synaptic_transient", "peak_um"): Record(Basis.FROM_QUALITATIVE, "", (
        "1000 uM cleft peak. The 0.3-1 mM range is standard for a central synapse and is "
        "qualitatively well established, but no single measurement is cited for 1000 "
        "specifically. The fit's anchors are macroscopic and this value enters through "
        "them, so changing it moves the fitted rates.")),
    ("synaptic_transient", "rise_ms"): Record(Basis.CONVENTION, "", (
        "0.10 ms. A declared convention standing in for the solution-exchange time of a "
        "fast-perfusion experiment; the tau_rise -> tau_clear limit is handled explicitly "
        "(P0-2) because the closed form is singular there.")),
    ("synaptic_transient", "clear_ms"): Record(Basis.FROM_QUALITATIVE, "", (
        "1.00 ms single clearance tau. Kept as the documented special case of the "
        "biexponential form so existing results reproduce.")),
    ("synaptic_transient_biexp_clearance", "clear_fast_ms"): Record(
        Basis.FROM_QUALITATIVE, "", (
        "1.0 ms, the fast component the technical spec calls for. Diffusion out of the "
        "cleft is fast and this is the right order; the digit is not measured.")),
    ("synaptic_transient_biexp_clearance", "clear_slow_ms"): Record(
        Basis.FROM_QUALITATIVE, "", (
        "20.0 ms, mid-range of the 10-30 ms the spec gives for the slow component. The "
        "slow tail is what carries most of the charge -- a transient with it carries more "
        "than twice the integral of one without -- so this is load-bearing for every "
        "CHARGE number and is UNSOURCED in its digit.")),
    ("synaptic_transient_biexp_clearance", "weight_fast"): Record(Basis.GUESS, "", (
        "0.80. No source at all. It trades the fast component against the slow one and "
        "therefore sets the charge directly; 0.8 was chosen because it makes the fast "
        "component dominant, which is qualitatively right and quantitatively arbitrary.")),
    ("pulse_train", "freq_hz"): Record(Basis.CONVENTION, "", (
        "10 / 50 / 100 Hz, the three frequencies the design asks for. A protocol choice, "
        "not a measurement -- and the three are chosen to span the range over which this "
        "scheme's paired-pulse depression changes threefold.")),
    ("ambient_with_spillover", "ambient_um"): Record(Basis.FROM_QUALITATIVE, "", (
        "0.40 uM, the midpoint of the 0.2-0.8 uM range reported for cortex and "
        "hippocampus. See config.AMBIENT_GABA_UM: the range is attributable, the midpoint "
        "is a choice, and the asymptotic headroom moves from 482x to 33x across it -- so "
        "this one number spans more than a decade of the headline result.")),
    ("ambient_with_spillover", "spillover_peak_um"): Record(Basis.GUESS, "", (
        "2.0 uM above baseline. NO SOURCE. Spillover amplitude at an extrasynaptic "
        "receptor depends on release probability, uptake, geometry and distance, none of "
        "which this model represents; 2.0 uM is five times ambient and was chosen to be "
        "visible rather than measured. Any charge ratio computed with spillover on "
        "inherits it.")),
    ("ambient_with_spillover", "tau_spillover_ms"): Record(Basis.GUESS, "", (
        "30.0 ms. NO SOURCE. Longer than the synaptic clearance because extrasynaptic "
        "clearance is diffusion- and uptake-limited rather than cleft-limited, which is "
        "the correct direction and no more than that.")),
    ("peak_dose_response", "application_ms"): Record(Basis.CONVENTION, "", (
        "300 ms, shared with gabaa_kinetics.Scheme.po_peak. NOT arbitrary and NOT to be "
        "changed: a 1 ms step is binding-rate-limited at low agonist, which inflates the "
        "fitted EC50, forces an implausibly high microscopic affinity and compresses every "
        "derived PAM gain. Published concentration-response curves come from applications "
        "of hundreds of ms.")),
    ("peak_dose_response", "n_grid"): Record(Basis.CONVENTION, "", (
        "600 time samples, held at the value the solve_ivp implementation used so that "
        "replacing it with an exact propagator did not move any PEAK number for a reason "
        "unrelated to the model. The peak is a maximum over these samples.")),
}


ALL = {"REGIONS": REGIONS_PROV, "EXTRASYN": EXTRASYN_PROV,
       "SUBJECTIVE_WEIGHT": SUBJECTIVE_PROV, "PRIORS": PRIOR_PROV,
       "MODEL_PARAMS": MODEL_PARAM_PROV, "DATASETS": DATASET_PROV,
       "WAVEFORMS": WAVEFORM_PROV}


#: The tables that feed the SELECTIVITY RANKING -- the project's one VALIDATED result.
#:
#: This subset exists because "6 of 28 parameters name a source" is a claim about the
#: ranking's inputs, and the audit later grew to cover the receptor-model defaults and the
#: benchmark datasets as well (43 records). Those are a different layer: they feed the
#: receptor-kinetic results, not the ranking. Folding them into one number would silently
#: change what the manuscript's sentence means, which is the same scoping error as
#: comparing `subjective_index` with `regional_sens` -- see `subtypes.subjective_index`.
RANKING_INPUT_TABLES = ("REGIONS", "EXTRASYN", "SUBJECTIVE_WEIGHT", "PRIORS")


def audit(tables: tuple[str, ...] | None = None) -> dict:
    """Count the bases across registered parameters.

    `tables` selects a scope by name; `None` audits everything. Pass
    `RANKING_INPUT_TABLES` for the figure the manuscript quotes.
    """
    names = tuple(ALL) if tables is None else tuple(tables)
    unknown = [n for n in names if n not in ALL]
    if unknown:
        raise KeyError(f"unknown provenance table(s) {unknown}; have {list(ALL)}")

    tally: dict = {}
    for name in names:
        for rec in ALL[name].values():
            tally[rec.basis.value] = tally.get(rec.basis.value, 0) + 1
    total = sum(tally.values())
    sourced = sum(1 for n in names for r in ALL[n].values() if r.is_sourced)
    return dict(total=total, sourced=sourced, by_basis=tally, tables=names,
                sourced_fraction=sourced / total if total else 0.0)


def report() -> str:
    a = audit()
    rank = audit(RANKING_INPUT_TABLES)
    out = ["=" * 78,
           "PARAMETER PROVENANCE",
           "=" * 78,
           f"SELECTIVITY-RANKING INPUTS: {rank['sourced']} of {rank['total']} name a "
           f"source ({100*rank['sourced_fraction']:.0f}%). This is the figure the "
           f"manuscript quotes;",
           "it is the scope of the project's one VALIDATED result.",
           f"ALL REGISTERED PARAMETERS: {a['sourced']} of {a['total']} "
           f"({100*a['sourced_fraction']:.0f}%), adding the receptor-model defaults and "
           f"the benchmark datasets.",
           "A named source means the derivation is traceable, NOT that the number is right:",
           "a5_dist resolves perfectly by DOI and does not support the numbers it is the",
           "obvious candidate for. See knowledge/06-source-provenance.md.",
           ""]
    for basis in Basis:
        n = a["by_basis"].get(basis.value, 0)
        if n:
            out.append(f"  {basis.value:<18} {n:3d}")
    for name, table in ALL.items():
        out.append(f"\n-- {name} " + "-" * max(0, 72 - len(name)))
        for k, rec in table.items():
            key = ".".join(k) if isinstance(k, tuple) else k
            src = f" [{rec.source_key}]" if rec.source_key else ""
            out.append(f"  {key:<22} {rec.basis.value:<18}{src}")
    return "\n".join(out)


if __name__ == "__main__":
    print(report())
