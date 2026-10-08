#!/usr/bin/env python3
"""CLI for circuitpharm — evaluate a compound's receptor profile with reliability tiers.

    python simulator.py calibration            # what is anchored and what is not
    python simulator.py subtypes               # regional subunit composition
    python simulator.py compound alogabat      # evaluate a named profile
    python simulator.py profile --a5 1.0 --a23 0.1 --occupancy 0.35
    python simulator.py rank                   # the calibration-independent ranking

This is deliberately a THIN wrapper: all evaluation logic lives in
`circuitpharm.evaluate`, so the library is usable without the CLI and the CLI cannot drift
away from the library. The previous version carried its own evaluation code and had drifted
in three ways that all flattered the model -- a single sensitivity applied to both receptor
pools (wrong by ~7x on one of them), dose escalation clipped at a hand-set ceiling, and an
overdose index printed as a plain number with a warning underneath that was duly ignored.

WHAT THE NUMBERS MEAN. Every quantity carries a tier:
    VALIDATED     reproduces something it was not fitted to, or is an identity, or
                  survives propagation of every parameter it depends on
    UNCALIBRATED  mechanism sound, no valid quantitative anchor -- ordering yes, scale no
    VOID          rests on something known to be invalid; its value is not printed

"Dose" here is modulator site OCCUPANCY. There is no pharmacokinetics in this package:
no mg/kg, no brain:plasma ratio, no time course.
"""
import argparse

from circuitpharm.evaluate import Compound, evaluate
from circuitpharm.config import calibration_report
from circuitpharm.subtypes import (PROFILES, REGIONS, SUBTYPES, EXTRASYN,
                                   SUBJECTIVE_WEIGHT)


def cmd_calibration(a):
    print(calibration_report())


def cmd_subtypes(a):
    print(f"{'region':<12}" + "".join(f"{s:>9}" for s in SUBTYPES))
    for r, f in REGIONS.items():
        print(f"{r:<12}" + "".join(f"{f[s]:9.2f}" for s in SUBTYPES))
    print(f"\n{'extrasyn':<12}" + "".join(f"{EXTRASYN[s]:9.2f}" for s in SUBTYPES))
    print(f"{'subj weight':<12}"
          + "".join(f"{SUBJECTIVE_WEIGHT[s]:9.1f}" for s in SUBTYPES))
    print("""
a5 and a2/a3 carry ethanol's discriminative stimulus; a1 and delta do not.
a1 dominates the medulla (respiratory); a5 is forebrain-restricted AND predominantly
extrasynaptic, which is why an affinity-type PAM has ~200x headroom on its pool and only
~1.1x on the synaptic one -- and therefore why the efficacy-ceiling safety argument is
weakest exactly where an a5-selective strategy places the drug.""")


def cmd_rank(a):
    rows = []
    for key in sorted(PROFILES):
        c = Compound.from_profile(key)
        r = c.selectivity_ratio()
        if r == r:                      # drop NaN (no respiratory burden at all)
            rows.append((r, key, PROFILES[key].name))
    rows.sort(reverse=True)
    print("SELECTIVITY RANKING — subjective drive per unit preBotC respiratory burden,")
    print("relative to a non-selective benzodiazepine (= 1.00).\n")
    print(f"{'ratio':>8}  {'key':<18} compound")
    print("-" * 74)
    for ratio, key, name in rows:
        print(f"{ratio:8.2f}  {key:<18} {name}")
    print("""
This is the one safety-relevant quantity here that is calibration-INDEPENDENT: it is a
ratio, so the unknown lumped sensitivity factor cancels. It survives 20,000-draw
propagation over every estimated parameter at the 5th percentile.

It RANKS compounds. It bounds nothing -- no safety margin, no overdose multiple, no dose.""")


def cmd_compound(a):
    # s_max omitted unless the user gave one, so the profile's ceiling is used
    kw = dict(occupancy=a.occupancy, nmda_block=a.nmda_block,
              glun2b_sel=a.glun2b_sel, glyr=a.glyr)
    if a.s_max is not None:
        kw["s_max"] = a.s_max
    print(evaluate(Compound.from_profile(a.key, **kw), n_seed=a.seeds))


def cmd_profile(a):
    # an arbitrary profile has no ceiling of its own, so fall back to the BZ-site range
    print(evaluate(Compound(
        name=a.name, a1=a.a1, a23=a.a23, a5=a.a5, d_a4=a.d_a4, eps=a.eps,
        occupancy=a.occupancy, s_max=(2.5 if a.s_max is None else a.s_max),
        nmda_block=a.nmda_block, glun2b_sel=a.glun2b_sel, glyr=a.glyr), n_seed=a.seeds))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add_dose_args(p):
        p.add_argument("--occupancy", type=float, default=0.35,
                       help="fraction of modulator sites bound (NOT a mass dose)")
        # DEFAULT None, not 2.5. A concrete default here silently overwrote each
        # profile's own measured ceiling: `simulator.py compound neurosteroid` reported a
        # 6.0-ceiling compound as 2.5.
        p.add_argument("--s-max", dest="s_max", type=float, default=None,
                       help="ligand intrinsic allosteric efficacy = max GABA EC50 shift "
                            "at full occupancy (~2-3x for classical BZ-site ligands). "
                            "THIS is what determines overdose protection; it is a property "
                            "of the individual molecule and must be measured")
        p.add_argument("--nmda-block", type=float, default=0.0)
        p.add_argument("--glun2b-sel", dest="glun2b_sel", type=float, default=1.0,
                       help="1 = GluN2B-selective, 0 = non-selective")
        p.add_argument("--glyr", type=float, default=1.0)
        p.add_argument("--seeds", type=int, default=4)

    sub.add_parser("calibration", help="what is anchored and what is not"
                   ).set_defaults(fn=cmd_calibration)
    sub.add_parser("subtypes", help="regional subunit composition and localisation"
                   ).set_defaults(fn=cmd_subtypes)
    sub.add_parser("rank", help="the calibration-independent selectivity ranking"
                   ).set_defaults(fn=cmd_rank)

    p1 = sub.add_parser("compound", help="evaluate a named profile")
    p1.add_argument("key", choices=sorted(PROFILES))
    add_dose_args(p1)
    p1.set_defaults(fn=cmd_compound)

    p2 = sub.add_parser("profile", help="evaluate an arbitrary receptor profile")
    for s in ("a1", "a23", "a5", "d_a4", "eps"):
        p2.add_argument(f"--{s.replace('_', '-')}", dest=s, type=float, default=0.0)
    p2.add_argument("--name", default="novel compound")
    add_dose_args(p2)
    p2.set_defaults(fn=cmd_profile)

    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
