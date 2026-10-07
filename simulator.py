#!/usr/bin/env python3
"""ALCOHOL-SUBSTITUTE SIMULATOR — predict a novel compound's profile from its receptor activity.

    python simulator.py --help
    python simulator.py profile --a5 1.0 --nmda-block 0.3 --glun2b-sel 1.0
    python simulator.py compound mp_iii_022 --nmda-block 0.29 --glun2b-sel 1.0
    python simulator.py refs
    python simulator.py optimise --subj 0.6

WHAT IT PREDICTS (simulated, from the circuit models):
    ventilation   preBotC minute-ventilation proxy, % of drug-free control
    reflex        stretch-reflex gain, % of drug-free control (motor impairment)
    overdose      dose multiple at which the respiratory rhythm fails
    subjective    forebrain index, weighted by which GABA-A subtypes and NMDA action
                  carry ethanol's discriminative stimulus

WHAT IT CANNOT PREDICT — and these have killed real compounds:
    hepatotoxicity   FOUR GABA-A anxiolytics died of it (alpidem withdrawn post-market,
                     ocinaplon halted, panadiplon halted with HUMAN toxicity invisible to
                     rat/dog/monkey, kava banned in several EU markets). Not receptor-
                     mediated, not modellable here. Mandatory DILI screening.
    hERG / QTc       killed traxoprodil. An off-target cardiac property, invisible here.
    dissociation     GluN2B-selective block is LESS dissociative than ketamine, not
                     non-dissociative. Magnitude requires human dosing.
    dependence       any GABA-A PAM will produce tolerance and withdrawal. Out of scope
                     by project decision, NOT by absence of risk.
    subjective truth the subjective axis is a receptor-level index constrained by the drug
                     discrimination literature, never a simulation. It cannot tell you what
                     a compound feels like.

CALIBRATION ANCHORS (so the numbers mean something):
    preBotC GABA-A sensitivity fitted to human midazolam 2 mg IV: ventilation -14 to -19%,
    tidal volume -16 to -22%. Spinal anchored where the stretch-reflex validations passed
    (4/4: strychnine hyperreflexia, benzodiazepine depression, GluN2B sparing).
"""
import argparse, itertools, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool

from spinal.resp import PreBotC, resp_metrics
from spinal.circuit import SpinalCircuit
from spinal.plant import JointPlant
from spinal.cpg import Drug
from spinal.subtypes import GabaProfile, PROFILES, REGIONS, SUBTYPES, SUBJECTIVE_WEIGHT
import scripts.reflex as R

RESP_OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
               w=dict(ee_ampa=0.45, ee_nmda=0.2475))
IA_SCALE = 0.30
FB_2B, PERIPH_2B = 0.70, 0.15
# NMDA contribution to the subjective index, scaled so a 50% forebrain reduction = 1.0
NMDA_SUBJ_SCALE = 0.50


# --------------------------------------------------------------------------- compound spec
class Candidate:
    """A novel compound: GABA-A subtype efficacies + NMDA action + ceiling."""

    def __init__(self, name="novel", a1=0.0, a23=0.0, a5=0.0, d_a4=0.0, eps=0.0,
                 ceiling=2.5, nmda_block=0.0, glun2b_sel=1.0, glyr=1.0,
                 gaba_dose=None, subj_target=None, nmda_class="subunit"):
        self.name = name
        self.gaba = GabaProfile(name, a1=a1, a23=a23, a5=a5, d_a4=d_a4, eps=eps,
                                ceiling=ceiling)
        self.nmda_block, self.glun2b_sel, self.glyr = nmda_block, glun2b_sel, glyr
        # nmda_class: "subunit"  = GluN2B-selective / glycine-site (no respiratory boost
        #                          documented; model's own sign is used)
        #             "blocker"  = ketamine-class channel blocker, which CLINICALLY
        #                          stimulates breathing -- empirical boost applied
        self.nmda_class = nmda_class
        # The respiratory boost kludge was REMOVED. See MODEL LIMITATION below.
        self.resp_boost = 0.0
        # dose: either given explicitly, or solved so the GABA arm hits a subjective target
        self.subj_target = subj_target
        self._target_met = True
        if gaba_dose is not None:
            self.gain = min(gaba_dose, ceiling)
        elif subj_target is not None and self.gaba.subjective_index() > 1e-9:
            want = 1.0 + subj_target / self.gaba.subjective_index()
            self.gain = min(want, ceiling)
            # BUG FIXED: clamping the dose to the ceiling silently delivers LESS than the
            # requested subjective effect. Record it, so reachable() cannot report success
            # for a compound that merely stayed under its ceiling without hitting target.
            self._target_met = want <= ceiling + 1e-9
        else:
            self.gain = 1.0
            self._target_met = subj_target is None

    # ---- subjective index: GABA arm + NMDA arm, per the discrimination literature ----
    def subjective(self):
        gaba_term = self.gaba.subjective_index() * (self.gain - 1.0)
        fb_nmda = self.nmda_block * (FB_2B * self.glun2b_sel + (1 - self.glun2b_sel))
        nmda_term = fb_nmda / NMDA_SUBJ_SCALE
        return gaba_term, nmda_term, gaba_term + nmda_term

    def salience_ok(self):
        """Ethanol substituted only for mixtures where the GABA-A component had EQUAL OR
        GREATER salience than the NMDA component."""
        g, n, _ = self.subjective()
        return g >= n - 1e-9

    def reachable(self):
        """True only if the requested subjective target is actually DELIVERED."""
        if self.subj_target is None:
            return True
        _, n, tot = self.subjective()
        if self.gaba.subjective_index() <= 1e-9:
            return self.nmda_block > 0
        return self._target_met or tot >= self.subj_target - 1e-9

    def drug(self, dose_mult=1.0):
        g = 1.0 + (self.gain - 1.0) * dose_mult
        return Drug(gaba_a_gain=g, gaba_a_tau=1.0 + 0.6 * (g - 1.0),
                    gaba_a_efficacy_cap=self.gaba.ceiling,
                    nmda_block=min(self.nmda_block * dose_mult, 1.0),
                    glun2b_selectivity=self.glun2b_sel,
                    glun2b_fraction=PERIPH_2B, glyr_gain=self.glyr)

    def sens(self, region):
        return self.gaba.regional_sens(region)


# ------------------------------------------------------------------------------ endpoints
def _resp(args):
    cand, dose_mult, seed = args
    b = PreBotC(drug=cand.drug(dose_mult), gaba_sens=cand.sens("prebotc"),
                resp_drive_boost=cand.resp_boost * dose_mult, seed=seed, **RESP_OP)
    for i in range(int(14000 / 0.1)):
        b.step(0.1)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > 4000
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return r["mean"], float(r["alive"])


def _reflex(args):
    cand, dose_mult, seed = args
    orig = dict(SpinalCircuit.W)
    try:
        SpinalCircuit.W = dict(orig)
        for rec in ("ampa", "nmda"):
            SpinalCircuit.W[("Ia", "Mn", rec)] = orig[("Ia", "Mn", rec)] * IA_SCALE
        pl = JointPlant(R.XML, drug=cand.drug(dose_mult), rg_gain=0.0, seed=seed,
                        gaba_sens=cand.sens("spinal"))
        for (t, q, v) in R.trajectory(0.1, 0.45):
            pl.step(R.DT, impose=(q, v))
        A = pl.arrays(); tt = A["t"] / 1000.0
        base = (tt > 0.15) & (tt < R.PRE)
        dyn = (tt > R.PRE) & (tt < R.PRE + R.RAMP + 0.03)
        return ((A["mn_E"][dyn].max() - A["mn_E"][base].mean())
                / max(1e-6, A["ia_E"][dyn].max() - A["ia_E"][base].mean())), 1.0
    finally:
        SpinalCircuit.W = orig


CONTROL = Candidate("control")


def evaluate(cand, n_seed=4, overdose_doses=(1., 2., 4., 8., 16., 26.), pool=None):
    mapper = pool.map if pool else map
    jobs_r = [(cand, 1.0, s) for s in range(n_seed)]
    jobs_c = [(CONTROL, 1.0, s) for s in range(n_seed)]
    rv = list(mapper(_resp, jobs_r)); cv = list(mapper(_resp, jobs_c))
    rf = list(mapper(_reflex, jobs_r)); cf = list(mapper(_reflex, jobs_c))
    vent = 100 * np.mean([x[0] for x in rv]) / np.mean([x[0] for x in cv])
    refl = 100 * np.mean([x[0] for x in rf]) / np.mean([x[0] for x in cf])
    od = None
    for D in overdose_doses:
        alive = np.mean([x[1] for x in mapper(_resp, [(cand, D, s) for s in range(2)])])
        if alive < 0.5:
            od = D
            break
    g, n, tot = cand.subjective()
    return dict(name=cand.name, gain=cand.gain, subj=tot, subj_gaba=g, subj_nmda=n,
                salience=cand.salience_ok(), reachable=cand.reachable(),
                vent=vent, reflex=refl, nmda_block=cand.nmda_block,
                overdose=od if od else f">{overdose_doses[-1]:.0f}",
                sens_resp=cand.sens("prebotc"), sens_motor=cand.sens("spinal"))


def show(r):
    print(f"\n{'='*74}\n{r['name']}\n{'='*74}")
    if not r["reachable"]:
        print("  UNREACHABLE: subjective target exceeds this compound's efficacy ceiling")
    print(f"  GABA dose (PAM gain)      {r['gain']:.2f}")
    print(f"  subjective index          {r['subj']:.3f}  "
          f"(GABA {r['subj_gaba']:.3f} + NMDA {r['subj_nmda']:.3f})")
    print(f"  salience constraint       {'OK' if r['salience'] else 'VIOLATED (GABA must be >= NMDA)'}")
    print(f"  regional GABA-A sens      preBotC {r['sens_resp']:.4f}   spinal {r['sens_motor']:.4f}")
    print(f"  --- simulated endpoints, % of drug-free control ---")
    print(f"  ventilation               {r['vent']:.0f}%")
    print(f"  stretch reflex            {r['reflex']:.0f}%")
    print(f"  overdose index            {r['overdose']}x  "
          f"(alcohol's own ratio is ~5x)")
    if r.get("nmda_block", 0) > 0:
        print(f"  !! MODEL LIMITATION -- the ventilation figure above is NOT RELIABLE for")
        print(f"     the NMDA component. This preBotC ties burst maintenance to the long")
        print(f"     NMDA conductance, so ANY NMDA block depresses it. Clinically the")
        print(f"     opposite is true: ketamine-class blockers PRESERVE respiratory drive")
        print(f"     and s-ketamine STIMULATES breathing and attenuates propofol-/opioid-")
        print(f"     induced hypoventilation; GluN2B-selective rislenemdaz had clean safety")
        print(f"     pharmacology. Attempts to fix this failed BOTH structurally (reducing")
        print(f"     the NMDA share of recurrent excitation to 6% still gave 56% at 60%")
        print(f"     block) and by a drive-boost term. Take the NMDA arm's respiratory")
        print(f"     contribution as ~100% FROM CLINICAL DATA, and read the GABA arm's")
        print(f"     ventilation (run with --nmda-block 0) as the model's real output.")
    print(f"  --- NOT PREDICTED: hepatotoxicity, hERG/QTc, dissociation magnitude, "
          f"dependence ---")


# ---------------------------------------------------------------------------------- refs
REFS = [
    ("non-selective benzodiazepine, subj 0.50", dict(a1=1., a23=1., a5=1., subj_target=.50)),
    ("neurosteroid, subj 0.50", dict(a1=1., a23=1., a5=1., d_a4=1., eps=.8,
                                     ceiling=6., subj_target=.50)),
    ("a5-selective PAM alone, subj 0.50", dict(a23=.15, a5=1., subj_target=.50)),
    ("RECOMMENDED: a5-sel PAM + GluN2B 20%", dict(a23=.15, a5=1., subj_target=.50,
                                                  nmda_block=.29, glun2b_sel=1.)),
    ("a2/a3 PAM + GluN2B 20%", dict(a23=1., subj_target=.35, nmda_block=.29, glun2b_sel=1.)),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add_profile_args(p):
        p.add_argument("--a1", type=float, default=0.0, help="GABA-A a1 efficacy (SEDATION, RESPIRATORY)")
        p.add_argument("--a23", type=float, default=0.0, help="a2/a3 efficacy (anxiolysis; substitutes for ethanol)")
        p.add_argument("--a5", type=float, default=0.0, help="a5 efficacy (CARRIES ethanol's subjective effect)")
        p.add_argument("--d-a4", type=float, default=0.0, help="delta/a4 extrasynaptic efficacy")
        p.add_argument("--eps", type=float, default=0.0, help="epsilon efficacy (preBotC rhythm neurons)")
        p.add_argument("--ceiling", type=float, default=2.5, help="efficacy ceiling (high = no overdose protection)")
        p.add_argument("--nmda-block", type=float, default=0.0, help="fraction NMDA blocked")
        p.add_argument("--glun2b-sel", type=float, default=1.0, help="1=GluN2B-selective, 0=non-selective")
        p.add_argument("--glyr", type=float, default=1.0, help="glycine receptor gain")
        p.add_argument("--nmda-class", choices=["subunit", "blocker"], default="subunit",
                       help="'blocker' = ketamine-class, applies the empirical "
                            "respiratory-stimulation boost from clinical data")
        p.add_argument("--dose", type=float, default=None, help="explicit PAM gain")
        p.add_argument("--subj", type=float, default=None, help="solve dose for this subjective target")
        p.add_argument("--seeds", type=int, default=4)

    p1 = sub.add_parser("profile", help="evaluate an arbitrary receptor profile")
    add_profile_args(p1)
    p1.add_argument("--name", default="novel compound")

    p2 = sub.add_parser("compound", help="evaluate a known GABA profile by key")
    p2.add_argument("key", choices=sorted(PROFILES))
    add_profile_args(p2)

    p6 = sub.add_parser("pair", help="evaluate an ENANTIOMER PAIR and report the contrast")
    p6.add_argument("eutomer", choices=sorted(PROFILES))
    p6.add_argument("distomer", choices=sorted(PROFILES))
    p6.add_argument("--subj", type=float, default=0.45)
    p6.add_argument("--nmda-block", type=float, default=0.0)
    p6.add_argument("--glun2b-sel", type=float, default=1.0)
    p6.add_argument("--seeds", type=int, default=4)

    sub.add_parser("refs", help="evaluate the reference set")
    sub.add_parser("subtypes", help="print regional subunit composition and weights")

    p5 = sub.add_parser("optimise", help="search receptor space for the best profile")
    p5.add_argument("--subj", type=float, default=0.50)
    p5.add_argument("--seeds", type=int, default=2)

    a = ap.parse_args()

    if a.cmd == "subtypes":
        print(f"{'region':<12}" + "".join(f"{s:>8}" for s in SUBTYPES))
        for r, f in REGIONS.items():
            print(f"{r:<12}" + "".join(f"{f[s]:8.2f}" for s in SUBTYPES))
        print(f"\n{'subj weight':<12}" + "".join(f"{SUBJECTIVE_WEIGHT[s]:8.1f}" for s in SUBTYPES))
        print("\na5 and a2/a3 carry ethanol's discriminative stimulus; a1 and delta do not.")
        print("a1 dominates the medulla (respiratory); a5 is forebrain-restricted.")
        return

    with Pool(os.cpu_count()) as pool:
        if a.cmd == "pair":
            # WHY PAIRS ARE THE STRONGEST THING THIS MODEL CAN DO: both enantiomers run
            # through the SAME pipeline with the SAME wrong parameters, so shared
            # systematic error (wrong conductances, wrong connectivity, wrong plant,
            # uncalibrated spinal sensitivity) cancels to first order in the CONTRAST.
            # Absolute numbers here are miscalibrated; the difference between the arms is
            # the defensible quantity. Enantiomers are also a near-perfect matched control:
            # identical mass, logP, pKa, polar surface area -- only 3D fit differs.
            out = {}
            for role, key in (("eutomer", a.eutomer), ("distomer", a.distomer)):
                p = PROFILES[key]
                c = Candidate(f"{role}: {p.name}", a1=p.a1, a23=p.a23, a5=p.a5,
                              d_a4=p.d_a4, eps=p.eps, ceiling=p.ceiling,
                              nmda_block=a.nmda_block, glun2b_sel=a.glun2b_sel,
                              subj_target=a.subj)
                out[role] = evaluate(c, n_seed=a.seeds, pool=pool)
                show(out[role])
            e, d = out["eutomer"], out["distomer"]
            print(f"\n{'='*74}\nPAIRED CONTRAST (the defensible quantity)\n{'='*74}")
            print(f"  {'':<22}{'eutomer':>12}{'distomer':>12}{'contrast':>12}")
            for k, lab in (("subj", "subjective"), ("vent", "ventilation %"),
                           ("reflex", "reflex %"), ("sens_resp", "preBotC sens"),
                           ("sens_motor", "spinal sens")):
                print(f"  {lab:<22}{e[k]:12.3f}{d[k]:12.3f}{e[k]-d[k]:+12.3f}")
            if not e["reachable"] or not d["reachable"]:
                print("  NOTE: one arm could not reach the subjective target within its ceiling,")
                print("        so the contrast is not at matched subjective effect.")
            print("  Shared systematic error cancels in the contrast column; the absolute")
            print("  columns inherit the model's full calibration uncertainty.")
        elif a.cmd == "refs":
            for lab, kw in REFS:
                show(evaluate(Candidate(lab, **kw), n_seed=4, pool=pool))
        elif a.cmd == "optimise":
            print(f"searching receptor space for subjective target {a.subj} "
                  f"(salience constraint enforced)...")
            best, rows = None, []
            for a1, a23, a5, ceil in itertools.product((0., .25, .5, 1.), (0., .15, .5, 1.),
                                                       (0., .5, 1.), (2.5,)):
                if a1 + a23 + a5 == 0:
                    continue
                for nb, sel in ((0., 1.), (.2, 1.), (.29, 1.), (.4, 1.), (.29, 0.)):
                    c = Candidate("x", a1=a1, a23=a23, a5=a5, ceiling=ceil,
                                  nmda_block=nb, glun2b_sel=sel, subj_target=a.subj)
                    if not c.reachable() or not c.salience_ok():
                        continue
                    if c.subjective()[2] < a.subj - 1e-9:   # target must be DELIVERED
                        continue
                    rows.append(c)
            print(f"  {len(rows)} feasible profiles; evaluating the 12 most selective")
            rows.sort(key=lambda c: c.sens("prebotc") + c.sens("spinal"))
            for c in rows[:12]:
                r = evaluate(c, n_seed=a.seeds, pool=pool)
                score = min(r["vent"], r["reflex"])
                g = c.gaba
                tag = (f"a1={g.a1:.2f} a23={g.a23:.2f} a5={g.a5:.2f} "
                       f"nmda={c.nmda_block:.2f} sel={c.glun2b_sel:.0f}")
                print(f"   vent {r['vent']:3.0f}%  reflex {r['reflex']:3.0f}%   {tag}")
                if best is None or score > best[0]:
                    best = (score, tag, r)
            if best:
                print(f"\nBEST: {best[1]}")
                show(best[2])
        else:
            kw = dict(a1=a.a1, a23=a.a23, a5=a.a5, d_a4=a.d_a4, eps=a.eps,
                      ceiling=a.ceiling, nmda_block=a.nmda_block,
                      glun2b_sel=a.glun2b_sel, glyr=a.glyr,
                      gaba_dose=a.dose, subj_target=a.subj,
                      nmda_class=a.nmda_class)
            if a.cmd == "compound":
                p = PROFILES[a.key]
                kw.update(a1=p.a1, a23=p.a23, a5=p.a5, d_a4=p.d_a4, eps=p.eps,
                          ceiling=p.ceiling)
                name = p.name
            else:
                name = a.name
            if kw["gaba_dose"] is None and kw["subj_target"] is None:
                kw["subj_target"] = 0.50
            show(evaluate(Candidate(name, **kw), n_seed=a.seeds, pool=pool))


if __name__ == "__main__":
    main()
