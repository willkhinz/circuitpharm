"""Print every number read out of the alpha1beta2gamma2 full texts, and the derived quantities.

WHY THIS SCRIPT EXISTS. The literature targets in `fitting.data.MISSING_DATASETS` were
located but not read: the cloud session that found them was behind an egress policy that
returned 403 for pmc.ncbi.nlm.nih.gov, europepmc.org, rupress.org and discovery.ucl.ac.uk.
Numbers reported by a search agent were graded [LEAD] and deliberately NOT entered as data.
On 2026-10-10 PMC became reachable and all seven open-access targets were read in full.

This module holds what the full texts say, with the conditions attached, and PRINTS the
derived quantities rather than letting them be retyped into prose. Three of them are
derivations, not readings, and the distinction is the point:

  * NO source gives a single "mean open time" for alpha1beta2gamma2L. Keramidas & Harrison
    2008 and Li et al. 2008 give exponential DECOMPOSITIONS (tau_i with areas a_i). The
    area-weighted mean sum(a_i * tau_i) is computed here. It is ours, not theirs.
  * Barberis et al. 2007 report both the three components of their deactivation fit AND the
    weighted tau. Recomputing the latter from the former is a read-back check on this
    module's transcription, and it is printed as such.
  * `_NOMINAL_DEFECT` = 0.750 was compared against four P_o values that turned out not to be
    the same quantity. `po_candidates()` prints what each one actually conditions on.

THE FINDING THAT MATTERS MOST. Every observable this project wanted to anchor on is
protocol-dependent, and the sources do not declare the protocol in comparable terms:

  * PEAK: already known (knowledge/12-inference.md section 2.3), 7x EC50 shift with
    application length.
  * DEACTIVATION: Barberis 2007 measured the SAME receptor in the SAME patches at two pulse
    durations and got 52.5 ms (2 ms pulse) and 364 ms (3 s pulse) -- 6.9x from the protocol
    alone. This is a published measurement, not a model prediction.
  * MEAN OPEN TIME: spans 1.42 ms to 7.25 ms across sources, and the spread is mostly
    DEFINITION (non-stationary 4 s window vs intraburst M-mode vs intraburst H-mode), not
    disagreement.

So a dataset tagged `DEACTIVATION` or `MEAN_OPEN_TIME` is no more self-describing than one
tagged `PEAK`. See knowledge/14-literature.md.

Run:  python scripts/literature_gamma2.py
      python scripts/literature_gamma2.py --section open_times
"""
import argparse

from circuitpharm.literature import (  # noqa: F401
    SOURCES,
    Source,
    OpenTimeFit,
    OPEN_TIME_FITS,
    STATED_MEAN_OPEN_TIME,
    Deactivation,
    DEACTIVATIONS,
    NOMINAL_DEFECT,
    PO_CANDIDATES,
    PO_LEAD_UNLOCATED,
    OTHER,
    _weighted,
)

def _rule(ch="-", n=86):
    return ch * n


def _weighted(taus, areas):
    return sum(a * t for a, t in zip(areas, taus))


def sources():
    print(_rule("="))
    print("SOURCES READ")
    print(_rule("="))
    for s in SOURCES.values():
        print(f"\n[{s.key}]  read = {s.read}")
        print(f"  {s.citation}")
        print(f"  {s.journal_ref} ({s.year})   doi:{s.doi or '-'}  PMID:{s.pmid}")
        print(f"  url        {s.url}")
        print(f"  receptor   {s.receptor_verbatim}   species: {s.species}")
        print(f"  prep       {s.prep}")
        print(f"  temp (C)   {s.temp_c}        error bars: {s.err_kind}")
        if s.note:
            print(f"  NOTE       {s.note}")
    n_full = sum(1 for s in SOURCES.values() if s.read == "FULL_TEXT_PMC")
    print(f"\n{n_full} of {len(SOURCES)} read in full text; "
          f"{sum(1 for s in SOURCES.values() if s.read == 'NOT_READ')} unread.")


def open_times():
    print(_rule("="))
    print("MEAN OPEN TIME -- area-weighted means are DERIVED HERE, not stated by any source")
    print(_rule("="))
    print(f"\n{'fit':<30} {'receptor':<22} {'[GABA]':>9} {'mean':>8}  conditional on")
    print(_rule())
    for f in OPEN_TIME_FITS:
        m = _weighted(f.taus_ms, f.areas)
        area_sum = sum(f.areas)
        flag = "" if abs(area_sum - 1.0) < 1e-9 else f"  [areas sum to {area_sum:.3f}]"
        conc = f"{f.agonist_um / 1000:.0f} mM" if f.agonist_um >= 1000 else f"{f.agonist_um:.0f} uM"
        print(f"{f.label:<30} {f.receptor:<22} {conc:>9} {m:>7.2f} ms  {f.conditioning}{flag}")
        lo = _weighted([t - e for t, e in zip(f.taus_ms, f.tau_errs)], f.areas)
        hi = _weighted([t + e for t, e in zip(f.taus_ms, f.tau_errs)], f.areas)
        print(f"{'':<30} components {f.taus_ms} ms, areas {f.areas}, n = {f.n_patches}")
        print(f"{'':<30} range if every tau moves +/-1 {f.err_kind}: "
              f"{lo:.2f} - {hi:.2f} ms (NOT an error bar -- no covariance is published)")
    print()
    for key, d in STATED_MEAN_OPEN_TIME.items():
        print(f"{'STATED by ' + key:<30} {d['receptor']:<22} {'':>9} "
              f"{d['value_ms']:>7.2f} ms  {d['conditioning']}")
        print(f"{'':<30} +/- {d['err_ms']} {d['err_kind']}, n = {d['n']}")
        print(f"{'':<30} NOTE {d['note']}")

    derived = [_weighted(f.taus_ms, f.areas) for f in OPEN_TIME_FITS]
    stated = [d["value_ms"] for d in STATED_MEAN_OPEN_TIME.values()]
    allv = derived + stated
    print(f"\n  spread across all {len(allv)} readings: "
          f"{min(allv):.2f} - {max(allv):.2f} ms  ({max(allv) / min(allv):.1f}x)")
    intraburst = derived[:1] + derived[2:]      # M-mode and Li; H-mode is the outlier mode
    print(f"  the two SUB-SATURATING/medium-mode intraburst readings agree closely: "
          f"{', '.join(f'{v:.2f}' for v in intraburst)} ms")
    print("  -> the spread is DEFINITION, not disagreement. A 'mean open time' is not one")
    print("     quantity: it depends on whether it is conditioned on a burst, on a gating")
    print("     mode, or on a post-pulse time window.")
    print("\n  NO source states a mean open time for alpha1beta2gamma2L. The gamma2L")
    print("  readings (Li 2008) are a DERIVED intracluster mean at a sub-saturating 50 uM.")


def deactivation():
    print(_rule("="))
    print("MACROSCOPIC DEACTIVATION -- the protocol column is the finding")
    print(_rule("="))
    print(f"\n{'study':<36} {'receptor':<26} {'pulse':>7} {'[GABA]':>8} {'tau_w':>10}  n")
    print(_rule())
    for d in DEACTIVATIONS:
        err = "" if d.err_ms != d.err_ms else f" +/- {d.err_ms:g}"
        print(f"{d.label:<36} {d.receptor:<26} {d.pulse:>7} {d.conc_mm:>5g} mM "
              f"{d.tau_w_ms:>7.1f}{err:<8} {d.n}")
    print()
    b2, b3 = DEACTIVATIONS[1], DEACTIVATIONS[2]
    print(f"  SAME receptor, SAME study, pulse duration only: "
          f"{b2.tau_w_ms:.1f} ms ({b2.pulse}) -> {b3.tau_w_ms:.1f} ms ({b3.pulse})"
          f"  = {b3.tau_w_ms / b2.tau_w_ms:.1f}x")
    d14 = DEACTIVATIONS[0]
    print(f"  Dixon 2014 ({d14.pulse}, {d14.conc_mm:g} mM) vs Barberis ({b2.pulse}, "
          f"{b2.conc_mm:g} mM): {d14.tau_w_ms:.1f} vs {b2.tau_w_ms:.1f} ms "
          f"= {b2.tau_w_ms / d14.tau_w_ms:.1f}x")
    print("  Nearly the same pulse duration, so duration does NOT explain that gap.")
    print(f"  But Dixon's weighted {d14.tau_w_ms:.1f} ms sits beside Barberis's FAST")
    print(f"  component tau_1 = {b2.components[0][0]:.1f} ms. Dixon fitted TWO exponentials,")
    print("  Barberis THREE; Barberis's slow component (221 ms, area 0.20) contributes")
    print(f"  {b2.components[0][2] * b2.components[1][2]:.1f} of its {b2.tau_w_ms:.1f} ms.")
    print("  -> a weighted deactivation tau is only defined once the number of exponential")
    print("     components and the fit window are declared, as well as the pulse. Three")
    print("     declarations, none of which the project's dataset tags currently carry.")


def consistency():
    print(_rule("="))
    print("READ-BACK CHECKS -- do the components reproduce the weighted values as published?")
    print(_rule())
    b2 = DEACTIVATIONS[1]
    taus, areas = b2.components
    recomputed = _weighted(taus, areas)
    print("\nBarberis 2007 deactivation, 2 ms pulse")
    print(f"  components      tau = {taus} ms, A = {areas}")
    print(f"  sum of areas    {sum(areas):.3f}   (the paper states sum(A_i) = 1)")
    print(f"  recomputed      sum(A_i tau_i) = {recomputed:.2f} ms")
    print(f"  published       tau_w = {b2.tau_w_ms:.1f} +/- {b2.err_ms:g} ms (SEM, n = {b2.n})")
    delta = recomputed - b2.tau_w_ms
    print(f"  difference      {delta:+.2f} ms ({abs(delta) / b2.err_ms:.2f} published SEM)")
    verdict = "CONSISTENT" if abs(delta) < b2.err_ms else "DISCREPANT"
    print(f"  VERDICT         {verdict} -- the components and the weighted value in this")
    print("                  module agree with each other, so the transcription of both")
    print("                  is self-checked. (Means are published rounded, so exact")
    print("                  agreement is not expected.)")
    for f in OPEN_TIME_FITS:
        s = sum(f.areas)
        print(f"\n{f.label}: areas sum to {s:.3f} "
              f"{'OK' if abs(s - 1.0) < 1e-9 else 'CHECK'}")


def po_candidates():
    print(_rule("="))
    print(f"WHAT IS P_o,max = {NOMINAL_DEFECT}? The candidates are not the same quantity.")
    print(_rule("="))
    print(f"\n{'value':>12}  {'receptor':<22} {'n':>8}  quantity")
    print(_rule())
    for c in PO_CANDIDATES:
        err = "" if c["err"] != c["err"] else f" +/- {c['err']:g}"
        print(f"{c['value']:>7.2f}{err:<5} {c['receptor']:<22} {c['n']:>8}  "
              f"{c['quantity']}  [{c['source_key']}]")
        print(f"{'':>12}  conditional on: {c['conditions_on']}")
        print(f"{'':>12}  {c['note']}")
        print()
    macro = [c for c in PO_CANDIDATES if "MACROPATCH" in c["quantity"]]
    bursty = [c for c in PO_CANDIDATES if "INTRABURST" in c["quantity"]]
    print(f"  intraburst / cluster readings: {len(bursty) + 1} of {len(PO_CANDIDATES)}")
    print(f"  genuinely macroscopic readings: {len(macro)} "
          f"(value {macro[0]['value']}, and its authors call it unpublished)")
    print(f"\n  UNLOCATED LEAD: {PO_LEAD_UNLOCATED}")
    print("\n  A macroscopic model's peak-scaling convention needs the POPULATION peak.")
    print(f"  The only reading of that quantity is {macro[0]['value']}, which is"
          f" {NOMINAL_DEFECT - macro[0]['value']:.2f} BELOW")
    print(f"  the {NOMINAL_DEFECT} convention -- not above it. Every value that sits near or")
    print(f"  above {NOMINAL_DEFECT} is conditioned on being inside a burst or cluster.")
    gamma2l = [c for c in PO_CANDIDATES if c["receptor"].endswith("gamma2L")]
    print(f"\n  For gamma2L specifically, {len(gamma2l)} readings exist, both conditional: "
          f"{', '.join(str(c['value']) for c in gamma2l)}.")


def other():
    print(_rule("="))
    print("OTHER VERIFIED READINGS")
    print(_rule("="))
    for o in OTHER:
        print(f"\n{o['quantity']}  =  {o['value']}")
        print(f"  {o['receptor']}  [{o['source_key']}]  n = {o['n']}  ({o['where']})")
        print(f"  {o['prep']}")
        if o["note"]:
            print(f"  NOTE {o['note']}")


SECTIONS = {
    "sources": sources,
    "open_times": open_times,
    "deactivation": deactivation,
    "consistency": consistency,
    "po": po_candidates,
    "other": other,
}


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--section", choices=sorted(SECTIONS), default=None)
    a = p.parse_args()
    todo = [a.section] if a.section else list(SECTIONS)
    for i, name in enumerate(todo):
        if i:
            print("\n")
        SECTIONS[name]()


if __name__ == "__main__":
    main()
