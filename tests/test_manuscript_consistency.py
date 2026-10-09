"""Every quantitative claim in knowledge/08-manuscript.md must match the code.

WHY THIS FILE EXISTS. The previous manuscript draft quoted a complete, internally consistent
kinetic parameterisation that no commit of this repository produces -- Po_max 0.84 against the
code's 0.75, a 3.0 mM / 1.0 ms synaptic transient against 1.0 mM / 0.30 ms -- and cited a commit
hash present in no ref. Internal consistency is exactly why reading it did not catch the
problem: Kd, E, Po_max, the asymptote and the ratio were all mutually coherent at the wrong
anchor. The only thing that catches it is recomputing each figure and looking for the string.

So this is a string test on purpose. It greps the manuscript for the decimal representation of
each number it should contain, computed fresh. If someone retypes a table, moves an anchor, or
edits a figure by hand, the number stops appearing and this fails. It is deliberately blunt:
a test that parsed the tables could be fooled by the same drift it is meant to catch.

LIMITATION, stated so nobody over-trusts it. This proves the manuscript's numbers are the
model's numbers. It says nothing about whether the model is right, whether a citation supports
the claim attached to it, or whether a sentence of prose describes the number beside it.
Those are read, not tested; `knowledge/07-paper-review.md` is the record of that reading.
"""
import pathlib

import pytest

MS = pathlib.Path("knowledge/08-manuscript.md")


@pytest.fixture(scope="module")
def text():
    assert MS.exists(), f"{MS} is missing"
    return MS.read_text()


@pytest.fixture(scope="module")
def main_text(text):
    """The manuscript WITHOUT its Supplementary Note.

    Superseded values -- the rejected 0.84 anchor, the asymptote-derived falsification
    criterion, the four-decimal selectivity indices -- appear legitimately in the
    Supplementary Note, which records what changed and why. They must not appear in the main
    text, where they would read as live claims. Several assertions below were written against
    the whole document and failed on exactly that distinction, so the split is explicit.
    """
    i = text.find("## Supplementary Note")
    return text if i == -1 else text[:i]


def _table_rows(text, caption, end):
    """Just the pipe-delimited rows of one table, excluding its caption prose."""
    block = text[text.index(caption):text.index(end, text.index(caption))]
    return [ln for ln in block.split("\n") if ln.strip().startswith("|")]


@pytest.fixture(scope="module")
def scheme():
    from circuitpharm import gabaa_kinetics as gk
    return gk.fit_scheme(verbose=False)


def _require(text, s, what):
    assert s in text, f"{what}: the manuscript does not contain {s!r}"


def test_the_cited_commit_and_tag_exist(text):
    """The previous draft's reproducibility anchor was a hash present in no ref."""
    import re
    import subprocess

    refs = re.findall(r"`(manuscript-v\d+)`", text) + re.findall(r"commit `([0-9a-f]{7,40})`", text)
    assert refs, "the manuscript cites no commit or tag to reproduce it from"
    for ref in set(refs):
        r = subprocess.run(["git", "cat-file", "-t", ref], capture_output=True, text=True)
        assert r.returncode == 0, (
            f"the manuscript cites {ref!r}, which git cannot resolve -- the same defect as the "
            f"previous draft's c8f17a9")


def test_calibration_anchors_match_the_code(text, main_text, scheme):
    from circuitpharm import gabaa_kinetics as gk
    _require(text, f"{gk.FIT_TARGETS['ec50_um']:.1f}", "EC50 anchor")
    _require(text, f"{gk.FIT_TARGETS['po_max']:.3f}", "Po_max anchor")
    _require(text, f"{scheme.ec50_um():.4f}", "reproduced EC50")
    _require(text, f"{scheme.po_max():.4f}", "reproduced Po_max")
    # the pulse constants, which is where the previous draft diverged
    _require(text, f"{gk.SYNAPTIC_CLEAR_MS}", "synaptic clearance tau")
    assert "1.0 mM" in text, "the manuscript does not state the 1.0 mM synaptic peak"
    # 0.84 was a rejected calibration anchor. It may be discussed in the Supplementary Note;
    # in the main text it would read as a live value. 0.8470/0.8489 are P_open figures and
    # 0.8282 is the analytic gating bound, all legitimate.
    stripped = main_text
    for ok in ("0.8470", "0.8489", "0.8282"):
        stripped = stripped.replace(ok, "")
    assert "0.84" not in stripped, (
        "0.84 appears in the MAIN TEXT -- that is the rejected Po_max anchor, which belongs "
        "only in the corrections record")


def test_microscopic_rates_and_derived_constants_match(text, scheme):
    s = scheme
    for val, what in ((f"{s.kon:.7f}", "k_on"), (f"{s.koff:.6f}", "k_off"),
                      (f"{s.beta:.6f}", "beta"), (f"{s.alpha:.6f}", "alpha")):
        _require(text, val, what)
    _require(text, f"{s.koff / s.kon:.2f}", "K_d")
    _require(text, f"{s.beta / s.alpha:.4f}", "E = beta/alpha")
    _require(text, f"{s.alpha / s.beta:.5f}", "alpha/beta")
    # the analytic gating bound, which must appear SEPARATELY from the 0.75 anchor
    _require(text, f"{(s.beta / s.alpha) / (1 + s.beta / s.alpha):.4f}", "beta/(alpha+beta)")


def test_the_two_asymptotes_and_the_dynamic_range_match(text, scheme):
    import numpy as np
    from circuitpharm import gabaa_kinetics as gk
    s = scheme
    po_inf = 1.0 / (1.0 + (s.alpha / s.beta) * (1.0 + s.d / s.r))
    po_base = s.po_tonic(gk.AMBIENT_UM)
    _require(text, f"{po_inf:.6f}", "P_open,inf")
    _require(text, f"{po_inf / po_base:.1f}", "asymptotic dynamic range")
    # Po_base appears in scientific notation in the text
    mant = f"{po_base:.4e}".split("e")[0]
    _require(text, mant, "P_open,base mantissa")
    # the closed-form [G]_1/2
    E, Kd = s.beta / s.alpha, s.koff / s.kon
    _require(text, f"{Kd * (1 + np.sqrt(2 + E)) / (1 + E):.2f}", "[G]_1/2")


def test_table_4_compartment_divergence_matches(text, scheme):
    from circuitpharm import gabaa_kinetics as gk
    s = scheme
    c = gk.calibrate_pam(s, target_shift=2.5, kind="affinity")
    _require(text, f"{c:.4f}", "c_affinity for s_max 2.5")
    d = gk.derive(s, affinity=c)
    b, m = s.ipsc_metrics(), s.pam(affinity=c).ipsc_metrics()
    for val, what in ((f"{b['peak']:.6f}", "baseline synaptic peak"),
                      (f"{m['peak']:.6f}", "PAM synaptic peak"),
                      (f"{d['phasic_gain']:.3f}", "phasic gain"),
                      (f"{d['tau_ratio']:.3f}", "tau ratio"),
                      (f"{d['charge_ratio']:.3f}", "charge ratio"),
                      (f"{d['phasic_headroom']:.3f}", "phasic headroom"),
                      (f"{d['tonic_gain']:.3f}", "tonic gain"),
                      (f"{d['tonic_headroom']:.3f}", "tonic headroom")):
        _require(text, val, what)
    # the headline contrast between the two ceilings
    _require(text, f"{d['tonic_headroom'] / d['phasic_headroom']:.0f}", "ceiling contrast")


def test_falsification_intervals_are_the_reachable_ones(text, scheme):
    """The defect that would have mattered in a wet lab: the previous draft pre-registered
    asymptote values. Each reachable bound must appear, and the discredited '> 15x at 0.1 uM'
    must not."""
    from circuitpharm import gabaa_kinetics as gk
    s = scheme
    cs = [gk.calibrate_pam(s, target_shift=sm, kind="affinity") for sm in (2.40, 2.50)]
    for g in (0.1, 0.4, 1.0, 3.0, 10.0):
        b = s.po_tonic(g)
        for c in cs:
            _require(text, f"{s.pam(affinity=c).po_tonic(g) / b:.2f}",
                     f"falsification bound at {g} uM")
    # The discredited criterion may be QUOTED while being explained -- that is §5.2's job, and
    # it is the clearest way to make the asymptote/reachable distinction concrete. But every
    # occurrence must sit beside that explanation, and the pre-registered rules must be free
    # of it.
    i = 0
    while (i := text.find("R_max > 15", i)) != -1:
        near = text[max(0, i - 500):i + 500]
        assert "asymptote" in near, (
            "'R_max > 15' appears without the asymptote explanation nearby -- it reads as a "
            "live criterion")
        i += 1
    rules = text[text.index("**Falsification rules**"):text.index("### 5.3.")]
    assert "R_max > 15" not in rules, (
        "the pre-registered falsification rules still contain the asymptote-derived bound")


def test_ambient_sweep_and_the_cap_crossing_match(text, scheme):
    import numpy as np
    from circuitpharm import gabaa_kinetics as gk
    s = scheme
    po_inf = 1.0 / (1.0 + (s.alpha / s.beta) * (1.0 + s.d / s.r))
    for g in (0.10, 0.20, 0.40, 0.70, 0.80, 1.50, 3.00):
        _require(text, f"{po_inf / s.po_tonic(g):.1f}", f"dynamic range at {g} uM")
    A = 1.0 + s.beta / s.alpha + s.d / s.r
    _require(text, f"{A:.4f}", "A = 1 + E + d/r")
    x = (2.0 + np.sqrt(4.0 + 6.0 * A)) / (3.0 * A)
    _require(text, f"{x * (s.koff / s.kon):.2f}", "ambient GABA at which a 2.5x cap binds")


def test_selectivity_table_matches_subtypes_module(text):
    from circuitpharm.subtypes import (EXTRASYN, PROFILES, REGIONS, SUBJECTIVE_WEIGHT,
                                       SUBTYPES)
    RHO = 6.0
    w = lambda x: RHO * EXTRASYN[x] + (1.0 - EXTRASYN[x])
    drive = lambda p: sum(REGIONS["forebrain"][x] * getattr(p, x) * SUBJECTIVE_WEIGHT[x] * w(x)
                          for x in SUBTYPES)
    burden = lambda p: sum(REGIONS["prebotc"][x] * getattr(p, x) * w(x) for x in SUBTYPES)
    ref = PROFILES["nonselective_bz"]
    r0 = drive(ref) / burden(ref)
    # TWO SIGNIFICANT FIGURES, deliberately. The arithmetic is exact, but R is computed from
    # stylised efficacy vectors, nine unsourced subunit fractions, four unsourced subjective
    # weights and one fitted parameter. Printing 8.6442 implies a precision the inputs cannot
    # carry, so the manuscript rounds and this test enforces BOTH halves of that: the rounded
    # value must appear, and the 4-decimal value must NOT, so false precision cannot creep
    # back in via a later edit. Exact values stay available from paper_numbers.py.
    def two_sf(x):
        from decimal import Decimal
        if x == 0:
            return "0.00"
        import math
        exp = math.floor(math.log10(abs(x)))
        q = round(x, 1 - exp)
        return f"{q:.{max(0, 1 - exp)}f}"

    for k in ("ideal_a5", "alogabat", "mp_iii_022", "hz_166", "neurosteroid",
              "sh053_R", "sh053_S"):
        p = PROFILES[k]
        bd = burden(p)
        assert bd > 1e-15, f"{k} has zero burden; the nominal R is undefined"
        r = (drive(p) / bd) / r0
        _require(text, two_sf(r), f"nominal R for {k} at two significant figures")
        # Checked against the TABLE ROWS only. The captions and §4.4.4 cite a four-decimal
        # value precisely to explain why the table does not carry one.
        rows = "\n".join(_table_rows(text, "**Table 6. Selectivity ranking", "#### 4.4.2."))
        assert f"{r:.4f}" not in rows, (
            f"Table 6 reports R = {r:.4f} for {k} to four decimals; the inputs are stylised "
            "assumptions and unsourced fractions, so that precision is false")
    # and the input vectors the table reprints
    for region, key in (("prebotc", "preBötC"), ("forebrain", "Forebrain")):
        vals = ", ".join(f"{REGIONS[region][x]:.2f}" for x in SUBTYPES)
        _require(text, f"[{vals}]", f"{key} regional fractions")
    _require(text, ", ".join(f"{EXTRASYN[x]:.2f}" for x in SUBTYPES), "extrasynaptic fractions")


def test_the_disowned_statistics_are_not_presented_as_results(text):
    """`ranking_robustness.py` prints NOT QUOTABLE for the MC median and 95th percentile. The
    previous draft made both columns of its headline table. They may appear as a stated
    prior-sensitivity diagnostic -- which is where 4.4.3 puts them -- but not in Table 6."""
    start = text.index("**Table 6. Selectivity ranking")
    end = text.index("#### 4.4.2")
    table = text[start:end]
    assert "NOT QUOTABLE" in table, (
        "Table 6 does not carry the model's own quotability caveat")
    table = "\n".join([table] + _table_rows(text, "**Table 6. Selectivity ranking",
                                             "#### 4.4.2."))
    for banned in ("16.05", "84.85", "40.88", "218.65", "12.59", "60.98"):
        assert banned not in table, (
            f"Table 6 presents {banned}, a median or 95th percentile the model marks NOT "
            "QUOTABLE, as a result")


def test_provenance_claims_match_the_audit(text):
    """Table 5 must not attribute an UNSOURCED cell to a citation. The previous draft
    attributed the whole table to Pirker and Kasugai; the audit rates 9 of 15 cells
    UNSOURCED and one an outright GUESS."""
    from circuitpharm.provenance import ALL, RANKING_INPUT_TABLES, Basis, audit

    # SCOPED to the ranking's inputs, which is what the manuscript's sentence is about.
    # The audit later grew to cover the receptor-model defaults and the benchmark datasets
    # (43 records); those feed the receptor-kinetic results, not the ranking, and folding
    # them in would change what this claim means without the sentence changing. See
    # provenance.RANKING_INPUT_TABLES.
    a = audit(RANKING_INPUT_TABLES)
    # the headline provenance fraction the manuscript quotes, and its percentage
    _require(text, f"{a['sourced']} of {a['total']}", "the sourced-parameter fraction")
    _require(text, f"({100 * a['sourced_fraction']:.0f}%)", "the sourced-parameter percentage")
    assert a["by_basis"].get(Basis.GUESS.value, 0) >= 1, (
        "no GUESS-rated parameter remains; §4.2's wording assumes one")
    assert "UNSOURCED" in text and "GUESS" in text, (
        "the manuscript does not carry the audit's own labels")

    # Table 5 reprints a basis per cell. Every cell the audit rates UNSOURCED must be
    # labelled as such in the manuscript rather than attributed to a citation -- which is
    # exactly what the previous draft did with Pirker and Kasugai.
    table = "\n".join(_table_rows(text, "**Table 5.", "#### Provenance, stated"))
    # Derive the expected counts from the audit rather than asserting a remembered number.
    # The first version of this test asserted "at least 9 UNSOURCED" from the manuscript's own
    # prose, and the prose was wrong: the audit gives 7 UNSOURCED + 1 GUESS + 7
    # FROM_QUALITATIVE over these 15 cells. A test that trusts the document it is checking
    # checks nothing.
    from circuitpharm.provenance import EXTRASYN_PROV, REGIONS_PROV
    from circuitpharm.subtypes import SUBTYPES

    want: dict = {}
    for sub in SUBTYPES:
        want[EXTRASYN_PROV[sub].basis.value] = want.get(EXTRASYN_PROV[sub].basis.value, 0) + 1
    for region in ("prebotc", "forebrain"):
        for sub in SUBTYPES:
            b = REGIONS_PROV[(region, sub)].basis.value
            want[b] = want.get(b, 0) + 1
    assert sum(want.values()) == 15, f"Table 5 should describe 15 cells, audit has {want}"

    for basis, n in want.items():
        got = sum(1 for cell in table.split("|") if basis in cell)
        assert got == n, (
            f"Table 5 labels {got} cells {basis} where the audit rates {n}. The audit is "
            f"authoritative; the table is reprinting it.")
    assert Basis.QUANTITATIVE.value not in table, (
        "a Table 5 cell claims QUANTITATIVE provenance, which the audit does not grant any "
        "of them")


# ===================================================================== P7 additions
#
# The quantities P1-P6 added. Same blunt technique, same reason: these are the numbers a
# reader would quote, and the manuscript has already once carried an internally consistent
# set that no commit produced.
def test_the_sensitivity_table_matches_the_sweep(text, scheme):
    """Every row of Table 3, recomputed. The whole table was stale after P0-7's refit.

    It was stale in a way reading could not catch: the RATIO column had been updated and the
    two probability columns had not, so each row was internally inconsistent while the table
    as a whole still looked like a sweep. Recomputing every cell is the only check that sees
    it.
    """
    from dataclasses import replace

    s = scheme

    def pinf(sc):
        return 1.0 / (1.0 + (sc.alpha / sc.beta) * (1.0 + sc.d / sc.r))

    rows = {
        "nominal": (s, 0.40),
        "desensitisation entry d × 0.5": (replace(s, d=s.d * 0.5), 0.40),
        "desensitisation entry d × 2.0": (replace(s, d=s.d * 2.0), 0.40),
        "resensitisation r × 0.5": (replace(s, r=s.r * 0.5), 0.40),
        "resensitisation r × 2.0": (replace(s, r=s.r * 2.0), 0.40),
        "gating β × 0.5": (replace(s, beta=s.beta * 0.5), 0.40),
        "gating β × 2.0": (replace(s, beta=s.beta * 2.0), 0.40),
        "affinity K_d × 0.5": (replace(s, koff=s.koff * 0.5), 0.40),
        "affinity K_d × 2.0": (replace(s, koff=s.koff * 2.0), 0.40),
    }
    for g in (0.10, 0.20, 0.40, 0.70, 0.80, 1.50, 3.00):
        rows[f"ambient {g:.2f} µM"] = (s, g)

    table = text[text.index("**Table 3."):text.index("**Finite modulator shifts**")]
    for label, (sc, g) in rows.items():
        base = sc.po_tonic(g)
        for value, what in ((f"{base:.7f}", "P_open(ambient)"),
                            (f"{pinf(sc):.5f}", "P_open,inf"),
                            (f"{pinf(sc) / base:.1f}×", "dynamic range")):
            assert value in table, (
                f"Table 3, row {label!r}: {what} should be {value!r} and is not in the table")


def test_the_finite_modulator_table_matches(text, scheme):
    """Table 4's s_max and tonic gain per c_affinity. These moved MORE than the headline.

    At c = 50 the realised gain fell from 177.95x to 109.45x on the refit -- a 1.6-fold
    change in a row a reader might quote as "the ceiling" -- while the headline asymptote
    moved 1.7-fold. Neither is more important than the other and both are now recomputed.
    """
    s = scheme
    ec0 = s.ec50_um()
    base = s.po_tonic(0.40)
    table = text[text.index("**Finite modulator shifts**"):text.index("#### Key observations")]
    for c in (1.50, 2.00, 2.79, 5.00, 10.00, 50.00):
        mod = s.pam(affinity=c)
        for value, what in ((f"{ec0 / mod.ec50_um():.3f}", "s_max"),
                            (f"{mod.po_tonic(0.40) / base:.2f}×", "tonic gain")):
            assert value in table, (
                f"Table 4, c_affinity {c}: {what} should be {value!r}")


def test_the_robust_ambient_claim_matches_the_computed_range(main_text, scheme):
    """The rounded claim over 0.2-0.8 uM must not exceed what the sweep supports.

    THIS IS THE ONE A STRING TEST ALMOST MISSED. The previous draft said "> 50x" across
    0.2-0.8 uM, which was true of the superseded fit (lower end 54.7x) and is false of this
    one (32.6x) by a factor of 1.5. Every individual number in the sweep was updated and the
    SENTENCE was not, because no test asserted the relationship between them. This one
    does: whatever threshold the manuscript claims, it must be at or below the minimum of
    the computed range.
    """
    import re

    s = scheme
    po_inf = 1.0 / (1.0 + (s.alpha / s.beta) * (1.0 + s.d / s.r))
    worst = min(po_inf / s.po_tonic(g) for g in (0.2, 0.4, 0.7, 0.8))

    claims = [float(m) for m in re.findall(r'\*\*> (\d+)×\*\*', main_text)]
    assert claims, "the main text makes no '> Nx' extrasynaptic range claim at all"
    for claimed in claims:
        assert claimed <= worst, (
            f"the main text claims '> {claimed:g}x' but the range falls to {worst:.1f}x "
            f"within 0.2-0.8 uM, so the claim is false at the top of the ambient range")
    # and it should not be needlessly weak either -- within a factor of 2 of the truth
    assert max(claims) >= worst / 2.0, (
        f"the strongest claim is '> {max(claims):g}x' against a computed floor of "
        f"{worst:.1f}x, which understates the result by more than twofold")


def test_the_cap_crossing_claim_is_verified_not_asserted(main_text, scheme):
    """The ambient concentration where a 2.5x cap binds, and the rounded prose beside it."""
    import re

    import numpy as np

    s = scheme
    E, Kd = s.beta / s.alpha, s.koff / s.kon
    A = 1.0 + E + s.d / s.r
    x = (2.0 + np.sqrt(4.0 + 6.0 * A)) / (3.0 * A)
    exact = x * Kd
    _require(main_text, f"{exact:.2f}", "the cap-crossing concentration")

    # the rounded version in prose must round to the same place
    rounded = sorted({int(m) for m in re.findall(r'≈ \*?\*?(\d+) µM\*?\*? ambient', main_text)})
    assert rounded, "the main text does not state a rounded cap-crossing concentration"
    for r in rounded:
        assert abs(r - exact) < 1.0, (
            f"prose says the cap binds above ~{r} uM; the computed crossing is "
            f"{exact:.2f} uM")


def test_no_superseded_fit_value_survives_in_the_main_text(main_text):
    """The pre-P0-7 fit's numbers must not appear as live claims.

    The old fit was `k_on = 0.0146842, k_off = 0.469177, beta = 0.559023, alpha = 0.107891`,
    giving `K_d = 31.95`, `P_open,inf = 0.166169` and the famous `210.7x`. It was not wrong
    arithmetic -- it was the solution of an under-determined system that moved between SciPy
    versions (roadmap P0-7), which is why it was replaced rather than corrected. Any of its
    numbers still in the main text is a sentence that was not reread.

    The Supplementary Note may discuss them, which is why this runs on `main_text`.
    """
    superseded = {
        "0.0146842": "old k_on", "0.469177": "old k_off", "0.559023": "old beta",
        "0.107891": "old alpha", "31.95": "old K_d", "0.166169": "old P_open,inf",
        "210.7": "the old headline dynamic range", "5.1814": "old E",
        "0.8382": "old gating bound", "2.7854": "old c_affinity at s_max 2.5",
        "7.8847": "old P_open,base mantissa", "19.02": "old [G]_1/2",
        "3295.5": "old range at 0.1 uM", "351.1": "old d x 0.5 range",
    }
    found = {v: k for k, v in superseded.items() if k in main_text}
    assert not found, f"superseded fit values still in the main text: {found}"


def test_the_headline_number_is_never_quoted_without_its_scope(main_text, scheme):
    """Every occurrence of the headline ratio must sit near the word that bounds it.

    `123.3x` is an ASYMPTOTIC ratio of open probabilities at k_off -> 0+, which no ligand
    reaches, and the reachable gain at benchmark potency is 10.29x. Quoting the headline
    without that distinction nearby is the single most consequential way this manuscript
    could mislead, and it is a thing a string test can actually check.
    """
    s = scheme
    po_inf = 1.0 / (1.0 + (s.alpha / s.beta) * (1.0 + s.d / s.r))
    headline = f"{po_inf / s.po_tonic(0.40):.1f}"

    i, checked = 0, 0
    while (i := main_text.find(headline, i)) != -1:
        near = main_text[max(0, i - 700):i + 700].lower()
        assert any(w in near for w in ("asymptot", "k_off → 0", "ceiling", "not a prediction",
                                       "reachable", "limit")), (
            f"'{headline}' at offset {i} appears with no asymptote/reachable qualifier within "
            f"700 characters; it reads as a predicted effect size")
        checked += 1
        i += 1
    assert checked >= 5, (
        f"only {checked} occurrences of the headline checked; the manuscript used to quote "
        f"it far more often and a drop this large means the search string is wrong")
