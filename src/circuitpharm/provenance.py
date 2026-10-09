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

ALL = {"REGIONS": REGIONS_PROV, "EXTRASYN": EXTRASYN_PROV,
       "SUBJECTIVE_WEIGHT": SUBJECTIVE_PROV, "PRIORS": PRIOR_PROV}


def audit() -> dict:
    """Count the bases across every registered parameter."""
    tally: dict = {}
    for table in ALL.values():
        for rec in table.values():
            tally[rec.basis.value] = tally.get(rec.basis.value, 0) + 1
    total = sum(tally.values())
    sourced = sum(1 for t in ALL.values() for r in t.values() if r.is_sourced)
    return dict(total=total, sourced=sourced, by_basis=tally,
                sourced_fraction=sourced / total if total else 0.0)


def report() -> str:
    a = audit()
    out = ["=" * 78,
           "PARAMETER PROVENANCE",
           "=" * 78,
           f"{a['sourced']} of {a['total']} load-bearing parameters name a source "
           f"({100*a['sourced_fraction']:.0f}%).",
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
