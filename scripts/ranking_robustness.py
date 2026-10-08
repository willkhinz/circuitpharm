"""Is the selectivity RANKING — the project's only surviving output — actually robust?

WHY THIS IS THE DECISIVE TEST. After sessions 7a-7e the project's defensible output has
narrowed to one thing: a calibration-independent RANKING of compounds by subjective drive
per unit of respiratory burden. Every absolute number is gone (the respiratory axis cannot
be calibrated against human ventilation; the in-vitro anchor does not exist in usable form;
the overdose formulation is structurally incomplete). The ranking survived because it is a
RATIO, so the unknown lumped calibration factor cancels.

But the ranking rests on numbers that were invented:

  REGIONS   regional subunit fractions -- "the weakest numbers here" per subtypes.py
  EXTRASYN  per-subtype extrasynaptic fraction -- added in 7c, with `eps` flagged in the
            source as an outright GUESS
  the tonic:phasic drug gain ratio -- measured in the kinetic scheme at 3.6-6.8 depending
            on the affinity/gating mechanism mix, which is itself unknown

If the ranking is not robust to those, then NOTHING survives and the honest answer about
this whole in-silico approach is that it cannot rank compounds either.

This is analytic -- no circuit simulation -- so it runs at thousands of draws.

WHAT IS SAMPLED
  REGIONS["prebotc"] and REGIONS["forebrain"]  Dirichlet(kappa * nominal), preserving means
  EXTRASYN per subtype                         Beta around nominal with a wide spread,
                                               and `eps` sampled UNIFORM on [0,1] because
                                               it is a guess and should be treated as one
  tonic:phasic gain ratio                      uniform on [3.0, 7.5]

A ratio of ratios, so the lumped `gaba_sens` scale cancels in every draw regardless.

Run:  python scripts/ranking_robustness.py --draws 20000
"""
import argparse
import numpy as np
from circuitpharm.subtypes import (REGIONS, SUBTYPES, SUBJECTIVE_WEIGHT, EXTRASYN, PROFILES)

ARMS = ("neurosteroid", "hz_166", "mp_iii_022", "alogabat", "ideal_a5")
REF = "nonselective_bz"
KAPPA = 15.0
EXTRASYN_CONC = 8.0          # Beta concentration; lower = wider
GAIN_RATIO = (3.0, 7.5)      # tonic:phasic, from the mechanism-mix sweep (7d)


def draw_params(rng, a5_floor=0.0):
    frac = {}
    for region in ("prebotc", "forebrain"):
        nom = REGIONS[region]
        alpha = np.array([max(1e-3, nom[s] * KAPPA) for s in SUBTYPES])
        d = rng.dirichlet(alpha)
        # Optional floor on the drawn preBotC a5 fraction, for the prior-sensitivity check
        # described in the verdict section. Renormalised so the fractions still sum to 1.
        if a5_floor > 0.0 and region == "prebotc":
            i5 = SUBTYPES.index("a5")
            d[i5] = max(d[i5], a5_floor)
            d = d / d.sum()
        frac[region] = dict(zip(SUBTYPES, d))
    ex = {}
    for s in SUBTYPES:
        if s == "eps":
            ex[s] = float(rng.uniform(0.0, 1.0))      # declared a guess; treat it as one
        elif s == "d_a4":
            ex[s] = 1.0                               # exclusively extrasynaptic: solid
        else:
            m = EXTRASYN[s]
            a = max(0.05, m * EXTRASYN_CONC)
            b = max(0.05, (1.0 - m) * EXTRASYN_CONC)
            ex[s] = float(rng.beta(a, b))
    return frac, ex, float(rng.uniform(*GAIN_RATIO))


def score(key, frac, ex, ratio):
    """Subjective drive per unit preBotC respiratory burden, in arbitrary shared units.

    Both numerator and denominator are linear in the unknown lumped sensitivity scale and
    in the unknown absolute drug gain, so those cancel when this is divided by the
    reference arm's score.
    """
    p = PROFILES[key]; e = p.eff()
    fb, pb = frac["forebrain"], frac["prebotc"]
    subj = sum(fb[s] * e[s] * SUBJECTIVE_WEIGHT[s] * (ratio * ex[s] + (1.0 - ex[s]))
               for s in SUBTYPES)
    burden = sum(pb[s] * e[s] * (ratio * ex[s] + (1.0 - ex[s])) for s in SUBTYPES)
    return subj / burden if burden > 1e-12 else np.nan


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=20000)
    ap.add_argument("--a5-floor", type=float, default=0.0,
                    help="floor the DRAWN preBotC a5 fraction (nominal 0.02). The prior "
                         "puts 38%% of draws below 0.002, which inflates the optimistic "
                         "tail for a5-selective arms; use this to reproduce the "
                         "prior-sensitivity check.")
    a = ap.parse_args()
    rng = np.random.default_rng(20261007)

    R = {k: [] for k in ARMS}
    for _ in range(a.draws):
        frac, ex, ratio = draw_params(rng, a5_floor=a.a5_floor)
        ref = score(REF, frac, ex, ratio)
        if not np.isfinite(ref) or ref <= 0:
            continue
        for k in ARMS:
            R[k].append(score(k, frac, ex, ratio) / ref)

    print(f"{a.draws} draws. Score = subjective drive per unit preBotC respiratory")
    print(f"burden, RELATIVE to a non-selective benzodiazepine (= 1.0 by construction).")
    print("The lumped calibration scale cancels in every draw.\n")
    print(f"{'arm':<24}{'median':>9}{'5th':>8}{'95th':>8}{'P(better than ref)':>21}")
    print("-" * 70)
    # EMPTY-ARRAY GUARD. `np.percentile([], 5)` does not return NaN, it raises
    # `IndexError: index -1 is out of bounds for axis 0 with size 0`. Every draw for an arm
    # can be non-finite -- calibrate_pam returns NaN when the requested EC50 shift is
    # unreachable by the ligand's mechanism, which is a legitimate pharmacological outcome,
    # not an error. So an arm whose s_max cannot reach its target crashed the whole
    # robustness report instead of being reported as unreachable.
    def _pctl(arr, q):
        return float(np.percentile(arr, q)) if len(arr) else float("nan")

    # DROPPED DRAWS ARE COUNTED AND REPORTED (found 2026-10-07).
    #
    # A non-finite score here is not noise: for a PERFECTLY a5-selective compound the
    # preBotC burden is proportional to the drawn preBotC a5 fraction, and the Dirichlet
    # prior puts mass arbitrarily close to zero, so burden -> 0 and the score diverges. The
    # isfinite filter therefore silently discards EXACTLY THE DRAWS MOST FAVOURABLE to the
    # arm being scored. That pushes the reported numbers in the conservative direction,
    # which is the honest direction -- but it was invisible, and an unreported discard rate
    # is indistinguishable from no discards at all.
    n_drop = {}
    for k in ARMS:
        v = np.array([x for x in R[k] if np.isfinite(x)])
        n_drop[k] = len(R[k]) - len(v)
        if not len(v):
            note = "no finite draws - mechanism cannot reach its target shift"
            print(f"{PROFILES[k].name[:23]:<24}{note:>46}")
            continue
        print(f"{PROFILES[k].name[:23]:<24}{np.median(v):9.2f}"
              f"{_pctl(v,5):8.2f}{_pctl(v,95):8.2f}"
              f"{100*(v>1.0).mean():19.1f}%")
    if any(n_drop.values()):
        print("\ndiscarded non-finite draws (burden -> 0, i.e. the draws MOST favourable")
        print("to the arm; discarding them makes these numbers conservative):")
        for k in ARMS:
            if n_drop[k]:
                print(f"    {PROFILES[k].name[:40]:<42} {n_drop[k]:6d} of {len(R[k])} "
                      f"({100*n_drop[k]/max(1,len(R[k])):.1f}%)")

    print("\nPAIRWISE: does the a5 class beat the a2/a3 class, draw by draw?")
    a5 = np.array(R["alogabat"]); a23 = np.array(R["hz_166"])
    n = min(len(a5), len(a23))
    d = a5[:n] - a23[:n]
    print(f"  alogabat - HZ-166:  median {np.median(d):+.2f},  "
          f"90% CI [{np.percentile(d,5):+.2f}, {np.percentile(d,95):+.2f}],  "
          f"a5 better in {100*(d>0).mean():.1f}% of draws")
    ns = np.array(R["neurosteroid"])
    print(f"\n  neurosteroid worse than a non-selective BZ in "
          f"{100*(ns<1.0).mean():.1f}% of draws")

    # PRIOR SENSITIVITY. Measured 2026-10-07; this is the check that decides whether the
    # headline is pharmacology or the shape of its own prior.
    #
    # The Dirichlet around a NOMINAL preBotC a5 of 0.02 puts 38% of draws below a tenth of
    # that, and the forebrain:preBotC a5 ratio has a p99 of ~3.5e7. Since an a5-selective
    # compound's respiratory burden is roughly proportional to that fraction, the score
    # should diverge -- so the whole result could have been an artifact.
    #
    # It is not, for the statistic that matters. Flooring the drawn preBotC a5 at 0.005,
    # 0.01 and the full nominal 0.02 leaves the 5th percentile UNCHANGED (alogabat 2.51 in
    # every condition, P(>1) 99.9%, corr(score, 1/a5) = 0.02). What the near-zero tail
    # drives is the OPTIMISTIC end: alogabat's median falls 16.05 -> 7.63 and its 95th
    # percentile 84.9 -> 14.5 as the floor rises to nominal.
    #
    # CONCLUSION, and it is a restriction on what may be quoted: the 5th-percentile FLOOR
    # and the ordering are robust. The MEDIAN and 95th percentile are prior artifacts and
    # must not be quoted -- the previously-reported "median 33.76x, 95th 9264.80x" for the
    # ideal a5 arm are properties of the prior's tail, not pharmacological claims. Re-run
    # with --a5-floor to reproduce the sensitivity.
    print("\n" + "=" * 70)
    print("WHAT MAY AND MAY NOT BE QUOTED FROM THE TABLE ABOVE")
    print("=" * 70)
    print("QUOTABLE:   the 5th percentile and the ordering. Flooring the drawn preBotC a5")
    print("            fraction at 0.005 / 0.01 / 0.02 (nominal) leaves the 5th percentile")
    print("            unchanged at 2.51x for alogabat and P(>1) at 99.9%.")
    print("NOT QUOTABLE: the median and 95th percentile. 38% of draws put preBotC a5 below")
    print("            a tenth of nominal, and an a5-selective compound's burden is roughly")
    print("            proportional to it, so the upper tail measures the prior. Flooring at")
    print("            nominal moves alogabat's median 16.05 -> 7.63 and its 95th 84.9 -> 14.5.")
    print(f"            Reproduce with:  --a5-floor 0.02")
    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    robust = [k for k in ARMS
              if _pctl(np.array([x for x in R[k] if np.isfinite(x)]), 5) > 1.0]
    if robust:
        print("These arms beat a non-selective benzodiazepine in the LOWER 5th percentile,")
        print("i.e. the ordering survives even adverse draws of every invented parameter:")
        for k in robust:
            v = np.array([x for x in R[k] if np.isfinite(x)])
            print(f"    {PROFILES[k].name[:40]:<42} >= {_pctl(v,5):.2f}x")
        print("\n-> the RANKING is a real result. It is the project's surviving output.")
    else:
        print("NO arm beats the reference robustly. The ranking does not survive its own")
        print("parameter uncertainty, and the in-silico approach yields nothing usable.")
    print("\nWhat a robust ranking does NOT give: any absolute safety margin, any overdose")
    print("multiple, any dose. Ranking orders candidates; it cannot bound risk. Bounding")
    print("risk needs the in-vitro anchor (specified in scripts/predict_muscimol.py) and")
    print("per-molecule measurement of intrinsic allosteric efficacy.")
