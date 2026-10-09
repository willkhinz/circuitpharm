"""Emit every number the manuscript quotes, computed from the code that is in the repo.

WHY THIS SCRIPT EXISTS. The first manuscript draft (`initial_paper.md`) quoted a complete,
internally consistent set of kinetic numbers that NO commit of this repository produces. The
draft's three calibration anchors were EC50 = 20 uM, Po_max = 0.84, tau_IPSC = 15 ms with a
3.0 mM / tau_clear = 1.0 ms synaptic transient; `gabaa_kinetics.FIT_TARGETS` says Po_max =
0.75 and the pulse constants say 1.0 mM / 0.30 ms. `git log -S` finds neither 0.84 nor 3000
anywhere in this repository's history, and `FIT_RANGES` declares po_max acceptable only on
[0.70, 0.80], so 0.84 would have been REJECTED by the module's own validation. The draft also
cited a commit hash (c8f17a9) that does not exist in any ref.

The consequence was not a rounding difference. The draft's rhetorical core is the contrast
between a nearly-saturated synapse and a wide-open extrasynaptic space, and the synaptic half
of that contrast was the part that moved most:

    phasic peak gain       draft 1.062x   repo 1.319x
    phasic max headroom    draft 1.124x   repo 1.668x
    charge ratio           draft 1.655x   repo 2.419x
    tonic max headroom     draft 210.7x   repo 184.6x

The qualitative claim survives (two orders of magnitude separate the compartments) but every
quoted figure was wrong, and a reader reproducing the draft from this repository would have
found nothing matching. Numbers typed into prose drift from the code silently; numbers printed
by the code cannot. So the manuscript's tables are generated here and pasted, never retyped.

Run:  python scripts/paper_numbers.py            # all sections
      python scripts/paper_numbers.py --section kinetics
"""
import argparse

import numpy as np

from circuitpharm import gabaa_kinetics as gk
from circuitpharm.subtypes import (EXTRASYN, PROFILES, REGIONS, SUBJECTIVE_WEIGHT, SUBTYPES)

#: Nominal tonic:phasic weighting used by the selectivity index. FITTED, not sourced --
#: it is measured in the kinetic scheme and lands between the peak- and charge-based
#: ratios, both of which this script prints so the choice stays visible.
RHO_NOMINAL = 6.0

#: The ambient concentrations the manuscript tabulates. The 0.2-0.8 uM span is the
#: literature range; 0.10, 1.50 and 3.00 are included to show the trend, not as claims
#: about any preparation.
AMBIENT_GRID = (0.10, 0.20, 0.40, 0.80, 1.50, 3.00)

#: Affinity multipliers for the finite-PAM rows. These are MICROSCOPIC multipliers on
#: k_off; the operational EC50 shift each one produces is printed beside it, because the
#: two are not equal and the draft's table conflated them (it listed c = 2.79 as
#: s_max = 2.50, where this parameterisation gives s_max = 2.40).
AFFINITY_GRID = (1.50, 2.00, 2.79, 5.00, 10.0, 50.0)


def _rule(ch="-", n=78):
    return ch * n


def kinetics(s):
    """Section 2: calibration, microscopic constants, and the two asymptotes."""
    out = [_rule("="), "KINETIC SCHEME (manuscript section 2)", _rule("=")]
    out.append(f"calibration anchors   {gk.FIT_TARGETS}")
    out.append(f"accepted ranges       {gk.FIT_RANGES}")
    out.append(f"synaptic transient    {gk.SYNAPTIC_PEAK_UM:.0f} uM peak, "
               f"tau_clear {gk.SYNAPTIC_CLEAR_MS} ms")
    out.append(f"ambient (tonic)       {gk.AMBIENT_UM} uM")
    out.append("")
    out.append("fitted microscopic rates")
    for k in ("kon", "koff", "beta", "alpha", "d", "r"):
        unit = "uM^-1 ms^-1" if k == "kon" else "ms^-1"
        out.append(f"  {k:6s} {getattr(s, k):.6g} {unit}")
    kd, e, dr = s.koff / s.kon, s.beta / s.alpha, s.d / s.r
    out.append("")
    # EACH DERIVED CONSTANT IS PRINTED AT THE PRECISION THE MANUSCRIPT QUOTES, and the rule
    # is worth stating because getting it wrong twice cost two verification rounds.
    # `verify_manuscript.py` checks the document's digits against this output by substring,
    # so a figure printed to FEWER digits than the document uses is unverifiable (E printed
    # 4.819, quoted 4.8195) and one printed to MORE digits is equally unverifiable, because a
    # correctly rounded 4.8195 is not a substring of 4.819466. Print what the document says.
    out.append(f"  Kd = koff/kon            {kd:.4g} uM   ({kd:.4f})")
    out.append(f"  E  = beta/alpha          {e:.4g}      ({e:.4f})")
    out.append(f"  alpha/beta               {s.alpha / s.beta:.5g}")
    out.append(f"  d/r                      {dr:.4g}")
    out.append("")
    out.append("reproduced calibration targets")
    out.append(f"  EC50                     {s.ec50_um():.4f} uM   "
               f"(target {gk.FIT_TARGETS['ec50_um']})")
    out.append(f"  Po_max = beta/(alpha+beta) {s.po_max():.4f}       "
               f"(target {gk.FIT_TARGETS['po_max']})")
    out.append(f"  tau_IPSC                 {s.ipsc_metrics()['tau']:.3f} ms  "
               f"(target {gk.FIT_TARGETS['tau_ms']})")
    out.append("")
    # The closed-form three-state asymptote, and the same quantity measured by pushing the
    # mechanism. They must agree; if they ever stop agreeing, one of the two is wrong.
    po_inf_closed = 1.0 / (1.0 + (s.alpha / s.beta) * (1.0 + dr))
    po_inf_pushed = s.pam(affinity=1e5).po_tonic(gk.AMBIENT_UM)
    po_base = s.po_tonic(gk.AMBIENT_UM)
    out.append(f"Po_inf, closed form      {po_inf_closed:.6f}")
    out.append(f"Po_inf, mechanism pushed {po_inf_pushed:.6f}   "
               f"(agree to {abs(po_inf_closed - po_inf_pushed):.2e})")
    out.append(f"Po_base at {gk.AMBIENT_UM} uM        {po_base:.6e}")
    out.append(f"asymptotic dynamic range {po_inf_closed / po_base:.1f}x")
    return out


def compartments(s):
    """Section 4.1: what a 2.5x EC50-shift PAM does in each compartment."""
    c = gk.calibrate_pam(s, target_shift=2.5, kind="affinity")
    d = gk.derive(s, affinity=c)
    b = s.ipsc_metrics()
    m = s.pam(affinity=c).ipsc_metrics()
    out = ["", _rule("="), "COMPARTMENT DIVERGENCE (manuscript section 4.1)", _rule("=")]
    out.append(f"c_affinity reproducing a 2.5x operational EC50 shift: {c:.4f}")
    out.append("")
    out.append(f"{'observable':26s} {'baseline':>12s} {'with PAM':>12s} {'ratio':>10s}")
    out.append(_rule())
    rows = [("synaptic peak Po", b["peak"], m["peak"], d["phasic_gain"]),
            ("synaptic tau_IPSC (ms)", b["tau"], m["tau"], d["tau_ratio"]),
            ("synaptic charge (ms)", b["charge"], m["charge"], d["charge_ratio"]),
            ("tonic Po (0.40 uM)", s.po_tonic(gk.AMBIENT_UM),
             s.pam(affinity=c).po_tonic(gk.AMBIENT_UM), d["tonic_gain"])]
    for nm, bv, mv, rt in rows:
        # Six decimals for exactness, three in parentheses because that is the precision the
        # manuscript's Table 4 quotes for the absolute values, and verify_manuscript.py
        # matches by substring. Same rule as the derived constants above: print what the
        # document says. Found by that script rejecting a correctly rounded 24.847.
        out.append(f"{nm:26s} {bv:12.6f} {mv:12.6f} {rt:9.3f}x"
                   f"   ({bv:.3f} -> {mv:.3f})")
    out.append(_rule())
    out.append(f"{'synaptic MAX headroom':26s} {'':12s} {'':12s} "
               f"{d['phasic_headroom']:9.3f}x")
    out.append(f"{'tonic MAX headroom':26s} {'':12s} {'':12s} "
               f"{d['tonic_headroom']:9.3f}x")
    out.append("")
    out.append("The ratio of the two maxima is the manuscript's central contrast: "
               f"{d['tonic_headroom'] / d['phasic_headroom']:.0f}x.")
    return out


def sensitivity(s):
    """Section 2.5: is the dynamic range an artefact of one parameter set?"""
    dr = s.d / s.r
    po_inf = 1.0 / (1.0 + (s.alpha / s.beta) * (1.0 + dr))
    out = ["", _rule("="), "SENSITIVITY (manuscript section 2.5)", _rule("=")]
    out.append(f"{'perturbation':34s} {'Po_base':>11s} {'Po_inf':>9s} {'range':>10s}")
    out.append(_rule())
    out.append(f"{'nominal':34s} {s.po_tonic(gk.AMBIENT_UM):11.7f} {po_inf:9.5f} "
               f"{po_inf / s.po_tonic(gk.AMBIENT_UM):9.1f}x")
    from dataclasses import replace
    for label, kw in (("desensitisation entry d x0.5", dict(d=s.d * 0.5)),
                      ("desensitisation entry d x2.0", dict(d=s.d * 2.0)),
                      ("resensitisation r x0.5", dict(r=s.r * 0.5)),
                      ("resensitisation r x2.0", dict(r=s.r * 2.0)),
                      ("gating beta x0.5", dict(beta=s.beta * 0.5)),
                      ("gating beta x2.0", dict(beta=s.beta * 2.0)),
                      ("affinity Kd x0.5 (koff x0.5)", dict(koff=s.koff * 0.5)),
                      ("affinity Kd x2.0 (koff x2.0)", dict(koff=s.koff * 2.0))):
        p = replace(s, **kw)
        pb = p.po_tonic(gk.AMBIENT_UM)
        pi = 1.0 / (1.0 + (p.alpha / p.beta) * (1.0 + p.d / p.r))
        out.append(f"{label:34s} {pb:11.7f} {pi:9.5f} {pi / pb:9.1f}x")
    out.append(_rule())
    for g in AMBIENT_GRID:
        pb = s.po_tonic(g)
        tag = " (nominal)" if abs(g - gk.AMBIENT_UM) < 1e-9 else ""
        out.append(f"{'ambient [G] = %.2f uM%s' % (g, tag):34s} {pb:11.7f} "
                   f"{po_inf:9.5f} {po_inf / pb:9.1f}x")
    out.append(_rule())
    out.append(f"{'finite c_affinity':34s} {'s_max':>11s} {'':9s} {'tonic gain':>10s}")
    for cc in AFFINITY_GRID:
        dd = gk.derive(s, affinity=cc)
        out.append(f"{'  c = %.2f' % cc:34s} {dd['ec50_shift']:11.3f} {'':9s} "
                   f"{dd['tonic_gain']:9.2f}x")
    # The concentration at which a 2.5x fixed cap would actually bind.
    A = 1.0 + s.beta / s.alpha + s.d / s.r
    # ratio(x) = 1 + 2/(A x) + 1/(A x^2) = 2.5  ->  1.5 A x^2 - 2 x - 1 = 0
    x = (2.0 + np.sqrt(4.0 + 6.0 * A)) / (3.0 * A)
    g_cross = x * (s.koff / s.kon)
    out.append("")
    out.append(f"A = 1 + E + d/r = {A:.4f}; a 2.5x cap binds only above "
               f"[G] = {g_cross:.2f} uM")
    out.append(f"  (check: range at that [G] = {po_inf / s.po_tonic(g_cross):.3f}x)")
    return out


def selectivity():
    """Section 4.4: the algebraic selectivity index. NOT a circuit simulation."""
    def wt(sub, rho):
        return rho * EXTRASYN[sub] + (1.0 - EXTRASYN[sub])

    def drive(p, rho):
        return sum(REGIONS["forebrain"][x] * getattr(p, x) * SUBJECTIVE_WEIGHT[x] * wt(x, rho)
                   for x in SUBTYPES)

    def burden(p, rho):
        return sum(REGIONS["prebotc"][x] * getattr(p, x) * wt(x, rho) for x in SUBTYPES)

    out = ["", _rule("="), "SELECTIVITY INDEX (manuscript section 4.4)", _rule("=")]
    out.append("ALGEBRAIC, not multiscale: this index is a ratio of weighted subunit sums.")
    out.append("`scripts/ranking_robustness.py` imports circuitpharm.subtypes and nothing")
    out.append("else -- no neuron, no circuit, no simulation. Calling it a 'multiscale")
    out.append("circuit pipeline' overstates it by a whole layer of the model.")
    out.append("")
    ref = PROFILES["nonselective_bz"]
    r0 = drive(ref, RHO_NOMINAL) / burden(ref, RHO_NOMINAL)
    out.append(f"rho = {RHO_NOMINAL} (FITTED); reference arm = nonselective_bz")
    out.append("")
    out.append(f"{'arm':22s} {'ceiling (repo)':>15s} {'nominal R':>11s}")
    out.append(_rule())
    for k in ("ideal_a5", "alogabat", "mp_iii_022", "hz_166", "nonselective_bz",
              "neurosteroid", "gaboxadol", "sh053_R", "sh053_S"):
        p = PROFILES[k]
        bd = burden(p, RHO_NOMINAL)
        r = (drive(p, RHO_NOMINAL) / bd) / r0 if bd > 1e-15 else float("inf")
        out.append(f"{k:22s} {p.ceiling:15.4g} {r:10.4f}x")
    out.append(_rule())
    out.append("`ceiling` is each profile's s_max in the repo. The draft's Table 6A listed")
    out.append("alogabat 2.40, MP-III-022 2.20 and HZ-166 2.20; all three are 2.5 here. R is")
    out.append("unaffected (s_max does not enter it) but the table misreported its own input.")
    return out


def falsification(s):
    """Section 5.2: the pre-registered intervals for the proposed patch-clamp protocol.

    THE DRAFT GOT TWO OF FOUR WRONG, and the reason is worth keeping in front of whoever
    edits this section. At low ambient GABA the ASYMPTOTIC headroom grows without limit, but
    the gain a FINITE s_max can actually extract does not. Between 0.4 and 0.1 uM the
    asymptote rises 15.6x while the realised gain rises 1.08x. The draft pre-registered
    "R_max > 15x at 0.1 uM", which is an asymptote reading; no affinity PAM at s_max ~2.5
    gets past ~8.5x there. A lab measuring 8x would have reported the model falsified when
    the model predicts 8x. So this prints both columns, always, side by side.
    """
    out = ["", _rule("="), "FALSIFICATION INTERVALS (manuscript section 5.2)", _rule("=")]
    cs = {sm: gk.calibrate_pam(s, target_shift=sm, kind="affinity") for sm in (2.40, 2.50)}
    for sm, c in cs.items():
        out.append(f"  s_max = {sm:.2f}  ->  c_affinity = {c:.4f}")
    dr = s.d / s.r
    po_inf = 1.0 / (1.0 + (s.alpha / s.beta) * (1.0 + dr))
    out.append("")
    out.append(f"{'[GABA]_bath':>12s} {'R @ s_max 2.40':>15s} {'R @ s_max 2.50':>15s} "
               f"{'asymptote':>12s}")
    out.append(_rule())
    for g in (0.1, 0.4, 1.0, 3.0, 10.0):
        b = s.po_tonic(g)
        r24 = s.pam(affinity=cs[2.40]).po_tonic(g) / b
        r25 = s.pam(affinity=cs[2.50]).po_tonic(g) / b
        out.append(f"{g:11.1f}  {r24:14.3f}x {r25:14.3f}x {po_inf / b:11.1f}x")
    out.append(_rule())
    out.append("The asymptote is the koff->0+ ceiling: unreachable by any finite s_max.")
    out.append("Pre-register against the REACHABLE columns, never the asymptote.")
    return out


SECTIONS = {"kinetics": kinetics, "compartments": compartments,
            "sensitivity": sensitivity, "selectivity": selectivity,
            "falsification": falsification}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--section", choices=sorted(SECTIONS), action="append",
                    help="emit only these sections (repeatable); default is all")
    a = ap.parse_args()
    want = a.section or sorted(SECTIONS)
    s = gk.fit_scheme(verbose=False)
    lines = []
    for name in ("kinetics", "compartments", "sensitivity", "falsification",
                 "selectivity"):
        if name not in want:
            continue
        fn = SECTIONS[name]
        lines += fn() if name == "selectivity" else fn(s)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
