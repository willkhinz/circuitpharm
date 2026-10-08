"""Subunit-resolved GABA-A pharmacology.

WHY THIS MODULE EXISTS. The model previously treated GABA-A as ONE lumped conductance with
a single regional sensitivity (`gaba_sens=0.15`, fitted to midazolam). That abstraction
hid the most important fact in the project:

    alpha5 mediates ethanol's DISCRIMINATIVE STIMULUS (a5 agonists mimic it, a5 inverse
    agonists block it, a1 antagonists do not attenuate it, delta is not required)
    alpha1 mediates SEDATION, RESPIRATORY DEPRESSION, ataxia and motor impairment

Those are different subtypes, and alpha5 is forebrain-restricted (hippocampus and
olfactory bulb highest; pons and medulla much lower) while alpha1 dominates the medulla
where the respiratory network lives. So subtype selectivity is not a refinement -- it is
the mechanism by which the subjective effect can be separated from the harm.

CALIBRATION IS PRESERVED. A region's effective drug sensitivity is

    sens(region, drug) = K * sum_over_subtypes( fraction[region][s] * efficacy[drug][s] )

with K fixed so that a NON-SELECTIVE benzodiazepine in the preBotC reproduces the fitted
value of 0.15 (which was anchored on human midazolam ventilation data). Selectivity then
moves sensitivity away from that anchor in a principled way rather than by refitting.

Subunit fractions are ESTIMATES assembled from regional expression literature, not
measured proteomics. They are the weakest numbers here and are the first thing to replace
with real data. Their ORDERING (a1 dominant in medulla, a5 forebrain-restricted) is solid;
their precise values are not.
"""
from dataclasses import dataclass, field

SUBTYPES = ("a1", "a23", "a5", "d_a4", "eps")

# fraction of regional GABA-A conductance carried by each subtype (estimates)
REGIONS = {
    # medulla / preBotC: a1-predominant in adult, plus a2/a3 on specific neurons,
    # extrasynaptic delta, a4, and epsilon enriched on NK1R+ rhythm neurons. a5 minimal.
    "prebotc":  dict(a1=0.60, a23=0.15, a5=0.02, d_a4=0.15, eps=0.08),
    # spinal cord: a1 and a2/a3 both substantial, a5 low
    "spinal":   dict(a1=0.45, a23=0.35, a5=0.05, d_a4=0.10, eps=0.05),
    # forebrain circuits carrying the subjective effect (hippocampus-weighted: >25% of
    # CA1/CA3 neurons express a5, primarily extrasynaptic)
    "forebrain": dict(a1=0.35, a23=0.25, a5=0.30, d_a4=0.08, eps=0.02),
}

# which subtypes carry ethanol's discriminative stimulus, and with what weight.
# a5: a5 agonists (QH-ii-066, panadiplon) mimic ethanol; L-655,708 blocks it.
# a2/a3: HZ-166 and YT-III-31 substituted FULLY for ethanol in rhesus without rate change.
# a1: a1 antagonists do NOT attenuate ethanol's stimulus -> weight 0.
# delta: not necessary for ethanol's stimulus -> weight 0.
SUBJECTIVE_WEIGHT = dict(a1=0.0, a23=1.0, a5=1.0, d_a4=0.0, eps=0.0)

# ---------------------------------------------------------------------------------------
# SYNAPTIC vs EXTRASYNAPTIC LOCALISATION. Added 2026-10-07 (session 7c) and it is not a
# refinement -- it decides the project's central safety argument.
#
# circuitpharm/gabaa_kinetics.py shows that an affinity-type PAM has utterly different headroom
# in the two pools, because they see different agonist concentrations:
#     SYNAPTIC       near-saturating cleft transient (~mM)  -> headroom ~1.1x
#     EXTRASYNAPTIC  ambient ~0.4 uM vs EC50 ~20 uM         -> headroom ~211x
# So WHERE a subtype sits determines how much a PAM can do to it, and a single
# `gaba_a_efficacy_cap` cannot express that.
#
# Consequence, verified numerically: the a5 arm needs gain 2.54 for a 0.50 subjective
# target. Against the hand-set ceiling of 2.5 (a synaptic-flavoured number) it is
# UNREACHABLE; against extrasynaptic headroom it is trivially reachable, and a single
# clinical-potency ligand already delivers 7.2x. The project assumed BOTH a tight ceiling
# (its only source of GABA-arm overdose protection) AND full subjective delivery by an
# a5-selective PAM. Those are the same parameter and cannot both hold.
#
# Values: a5 is predominantly extrasynaptic in hippocampus; delta/a4 are exclusively
# extrasynaptic; a1/a2/a3 are predominantly synaptic. ORDERING is solid, precise values are
# not -- and `eps` is a guess with no good basis, flagged accordingly.
EXTRASYN = dict(a1=0.15, a23=0.20, a5=0.80, d_a4=1.00, eps=0.50)   # eps = GUESS

# PER-REGION calibration constants. One global constant is wrong: each region has its own
# empirical anchor, and they differ by an order of magnitude.
#   prebotc -> 0.15  : fitted to human midazolam ventilation (-14 to -19% at a sedative dose)
#   spinal  -> 1.00  : the value at which the stretch-reflex validations PASSED
#                      (benzodiazepine 2.0x gave reflex 67-79%, i.e. a 21-33% reduction,
#                       matching diazepam's known reflex depression)
#   forebrain -> 1.00: reference region for the subjective index
# Applying the preBotC constant to the spinal cord silently abolished every reflex drug
# effect (all arms came out at ~100% of control), which is how this was caught.
K_REGION = {"prebotc": None, "spinal": None, "forebrain": 1.0}
ANCHOR = {"prebotc": 0.15, "spinal": 1.00, "forebrain": 1.00}


@dataclass
class GabaProfile:
    """A compound's positive-modulatory efficacy at each GABA-A subtype (0 = none)."""
    name: str
    a1: float = 0.0
    a23: float = 0.0
    a5: float = 0.0
    d_a4: float = 0.0
    eps: float = 0.0
    ceiling: float = 2.5          # efficacy cap (overdose protection); high = no ceiling
    hepatotox_flag: str = "unknown"

    def eff(self) -> dict:
        return {s: getattr(self, s) for s in SUBTYPES}

    def regional_sens(self, region: str) -> float:
        """Total drug-modulatable fraction. Unchanged — the sum of the two pools below,
        kept so existing callers and the per-region calibration are untouched."""
        f = REGIONS[region]; e = self.eff()
        return K_REGION[region] * sum(f[s] * e[s] for s in SUBTYPES)

    def regional_sens_split(self, region: str) -> tuple:
        """(extrasynaptic/tonic, synaptic/phasic) sensitivity in this region.

        Partitions regional_sens by subtype localisation (EXTRASYN). The two must be driven
        by DIFFERENT drug gains, because a PAM's headroom differs ~200-fold between them.
        Sums to regional_sens() by construction.
        """
        f = REGIONS[region]; e = self.eff(); K = K_REGION[region]
        tonic = K * sum(f[s] * e[s] * EXTRASYN[s] for s in SUBTYPES)
        phasic = K * sum(f[s] * e[s] * (1.0 - EXTRASYN[s]) for s in SUBTYPES)
        return tonic, phasic

    def subjective_index(self) -> float:
        """Forebrain subjective drive, weighted by which subtypes carry ethanol's stimulus.

        DELIBERATELY NOT ON THE SAME SCALE AS `regional_sens("forebrain")`, and the
        difference is a real trap worth stating. `regional_sens` multiplies by
        `K_REGION["forebrain"]` (1.111) so a non-selective benzodiazepine reproduces the
        regional anchor of 1.00. This function omits K, so the same compound yields 0.55 --
        a 1.8x discrepancy between two quantities that both sound like "forebrain drug
        engagement".

        The omission is intentional: this is NOT a sensitivity. It is a weighted sum over
        only those subtypes that carry ethanol's discriminative stimulus
        (`SUBJECTIVE_WEIGHT` zeroes alpha1 and delta), so applying a constant fitted to
        TOTAL regional conductance would be meaningless -- it would scale a subset by a
        normaliser derived from the whole.

        What follows is that the two must never be compared or combined, and that the
        subjective index has NO absolute scale at all: the forebrain anchor is a unit
        convention (see `circuitpharm.config` CALIBRATIONS, where it is recorded as
        UNCALIBRATED for exactly this reason), so only RATIOS between compounds mean
        anything. A "subjective target of 0.50" is a number in invented units.
        """
        f = REGIONS["forebrain"]; e = self.eff()
        return sum(f[s] * e[s] * SUBJECTIVE_WEIGHT[s] for s in SUBTYPES)

    def subjective_index_split(self) -> tuple:
        """(tonic, phasic) components of the forebrain subjective drive.

        Matters because a5 — the subtype carrying ethanol's discriminative stimulus — is
        predominantly EXTRASYNAPTIC. So the subjective effect is largely a TONIC-current
        phenomenon, and must be scaled by the tonic gain (large headroom), not the phasic
        one. Using a single gain for both is what produced the spurious reachability wall.
        """
        f = REGIONS["forebrain"]; e = self.eff()
        tonic = sum(f[s] * e[s] * SUBJECTIVE_WEIGHT[s] * EXTRASYN[s] for s in SUBTYPES)
        phasic = sum(f[s] * e[s] * SUBJECTIVE_WEIGHT[s] * (1.0 - EXTRASYN[s])
                     for s in SUBTYPES)
        return tonic, phasic


# ---- compound profiles -------------------------------------------------------------
# BZ-site ligands need gamma2, so they are inactive at delta- and epsilon-containing
# receptors. Neurosteroids act at a distinct transmembrane site and are subtype
# non-selective, so they engage delta and epsilon too.
PROFILES = {
    "nonselective_bz":  GabaProfile("non-selective BZ (diazepam/midazolam)",
                                    a1=1.0, a23=1.0, a5=1.0, ceiling=2.5),
    "zolpidem":         GabaProfile("zolpidem (a1-preferring)",
                                    a1=1.0, a23=0.3, a5=0.0, ceiling=2.5),
    "mp_iii_022":       GabaProfile("MP-III-022 (a5-selective, nonmodulatory at a1)",
                                    a1=0.0, a23=0.15, a5=1.0, ceiling=2.5),
    "hz_166":           GabaProfile("HZ-166 / KRM-II-81 (a2/a3)",
                                    a1=0.0, a23=1.0, a5=0.0, ceiling=2.5),
    "yt_iii_31":        GabaProfile("YT-III-31 (a3)",
                                    a1=0.0, a23=0.7, a5=0.0, ceiling=2.5),
    "tpa023":           GabaProfile("TPA023 (a2/a3 partial; silent antagonist a1 AND a5)",
                                    a1=0.0, a23=0.6, a5=0.0, ceiling=2.0),
    "l_838417":         GabaProfile("L-838417 (a2/a3/a5 partial; antagonist a1)",
                                    a1=0.0, a23=0.6, a5=0.6, ceiling=2.2),
    "neurosteroid":     GabaProfile("neurosteroid (non-selective, distinct site)",
                                    a1=1.0, a23=1.0, a5=1.0, d_a4=1.0, eps=0.8,
                                    ceiling=6.0),
    "gaboxadol":        GabaProfile("gaboxadol (delta direct agonist)",
                                    a1=0.0, a23=0.0, a5=0.0, d_a4=1.0, eps=0.0,
                                    ceiling=1e9),
    # THE CENTRAL CHIRAL PAIR OF THIS PROJECT. One stereocentre switches which GABA arm
    # the molecule serves. MP-III-022, the lead a5 candidate, derives from the R enantiomer.
    "sh053_R":          GabaProfile("SH-053-2'F-R-CH3 (a5-selective enantiomer)",
                                    a1=0.0, a23=0.05, a5=0.80, ceiling=2.5),
    "sh053_S":          GabaProfile("SH-053-2'F-S-CH3 (a2/a3/a5 enantiomer)",
                                    a1=0.0, a23=0.50, a5=0.50, ceiling=2.5),
    "alogabat":         GabaProfile("alogabat RG7816 (a5-selective, Phase 2)",
                                    a1=0.0, a23=0.10, a5=1.00, ceiling=2.5),
    "gl_ii_73":         GabaProfile("GL-II-73 (a5-preferring, no a1 affinity)",
                                    a1=0.0, a23=0.30, a5=1.00, ceiling=2.5),
    "imepitoin":        GabaProfile("imepitoin (non-selective, 20% of diazepam efficacy)",
                                    a1=0.20, a23=0.20, a5=0.20, ceiling=1.25),
    # hypothetical best-case: a5 only, zero elsewhere
    "ideal_a5":         GabaProfile("IDEAL a5-only PAM",
                                    a1=0.0, a23=0.0, a5=1.0, ceiling=2.5),
}

# fix each region's constant so a NON-SELECTIVE benzodiazepine reproduces that region's
# empirical anchor. Selectivity then moves sensitivity relative to that anchor.
_ns = PROFILES["nonselective_bz"]
for _r in ("prebotc", "spinal", "forebrain"):
    _raw = sum(REGIONS[_r][s] * _ns.eff()[s] for s in SUBTYPES)
    K_REGION[_r] = ANCHOR[_r] / _raw


def table():
    rows = []
    for key, p in PROFILES.items():
        subj = p.subjective_index()
        resp = p.regional_sens("prebotc")
        mot = p.regional_sens("spinal")
        # RATIOS ARE NORMALISED to a non-selective benzodiazepine. The raw quotient
        # subj/resp mixes two different scales: `subjective_index` omits
        # K_REGION["forebrain"] (it is a weighted SUBSET of subtypes) while `regional_sens`
        # includes it, so the bare number is in invented units and was printed as though
        # it meant something -- directly contradicting the warning in
        # `subjective_index`'s own docstring. Dividing by the reference compound's quotient
        # cancels both the K factor and the lumped sensitivity scale, which is the same
        # construction `Compound.selectivity_ratio` uses and the only form that is
        # interpretable.
        rows.append(dict(key=key, name=p.name, subj=subj, resp=resp, motor=mot,
                         raw_resp=(subj / resp if resp > 1e-9 else float("inf")),
                         raw_motor=(subj / mot if mot > 1e-9 else float("inf")),
                         ceiling=p.ceiling))
    ref = next(r for r in rows if r["key"] == "nonselective_bz")
    for r in rows:
        for k, rk in (("ratio_resp", "raw_resp"), ("ratio_motor", "raw_motor")):
            r[k] = (r[rk] / ref[rk] if ref[rk] not in (0.0, float("inf"))
                    and r[rk] != float("inf") else float("inf"))
    return rows


if __name__ == "__main__":
    print("per-region calibration (non-selective BZ reproduces each anchor):")
    for r in ("prebotc", "spinal", "forebrain"):
        print(f"  {r:<10} K={K_REGION[r]:7.4f}  anchor={ANCHOR[r]}")
    print()
    print("ratios are NORMALISED to a non-selective benzodiazepine (= 1.0), because the")
    print("raw quotient mixes two scales and is in invented units. See subjective_index().")
    print()
    print(f"{'compound':<46} {'subj':>6} {'resp':>6} {'motor':>6} "
          f"{'vs BZ resp':>11} {'vs BZ motor':>12}")
    print("-" * 92)
    for r in sorted(table(), key=lambda r: -r["ratio_resp"]):
        rr = f"{r['ratio_resp']:10.1f}" if r['ratio_resp'] != float("inf") else f"{'inf':>10}"
        rm = f"{r['ratio_motor']:11.1f}" if r['ratio_motor'] != float("inf") else f"{'inf':>11}"
        print(f"{r['name']:<46} {r['subj']:6.3f} {r['resp']:6.4f} {r['motor']:6.4f} {rr} {rm}")
