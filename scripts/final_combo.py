"""FINAL: choose the GABA-A component and report the recommended combination.

Two options for the GABA arm, with a real trade-off:

  BZ-SITE PAM (e.g. a benzodiazepine-site partial agonist)
    + ceiling-limited: requires endogenous GABA, saturates -> cannot cause apnoea
      (validated: plateaus at ~-44% ventilation even at 24x; an uncapped modulator
       reaches -85%)
    - lower subjective efficiency: ethanol generalises to benzodiazepines, but
      benzodiazepines do NOT generalise back to ethanol (one-way asymmetry), i.e. the
      BZ stimulus is only PART of ethanol's

  NEUROSTEROID (allopregnanolone-like, distinct transmembrane site)
    + higher subjective efficiency: ethanol is MORE likely to substitute for neuroactive
      steroids than for benzodiazepines
    + no extra respiratory penalty per unit effect: ganaxolone's preclinical sedation is
      "comparable to the benzodiazepine midazolam"
    - HIGH-EFFICACY neurosteroids have direct agonist activity at higher concentrations,
      so they may LOSE the overdose ceiling. Brexanolone carries a boxed warning for
      excessive sedation and SUDDEN LOSS OF CONSCIOUSNESS (4% of patients) requiring
      continuous pulse oximetry -- consistent with ceiling loss at high exposure.

`subj_eff` (subjective efficiency of the GABA arm) is an ASSUMPTION with only qualitative
literature support, so the recommendation is reported across a range of it rather than at
one value.
"""
import sys, os, itertools; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.circuit import SpinalCircuit
from circuitpharm.cpg import Drug
import scripts.reflex as R

FOREBRAIN_2B, PERIPHERAL_2B, GABA_SENS, IA_SCALE = 0.70, 0.15, 0.15, 0.30
N_SEED = 3
RESP_OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
               w=dict(ee_ampa=0.45, ee_nmda=0.2475))
# (label, efficacy ceiling, subjective efficiency)
GABA_ARMS = [("BZ-PAM", 2.5, 1.0), ("neurosteroid", 6.0, 1.5)]
PAMS = (1.0, 1.4, 1.8, 2.2, 2.6, 3.2)
FBS = (0.0, 0.10, 0.20, 0.30, 0.40)
SUBJ_TARGET = 1.20


def make(pam, fb, sel, cap):
    kw = dict(gaba_a_gain=pam, gaba_a_tau=1.0 + 0.6 * (pam - 1.0),
              gaba_a_efficacy_cap=cap)
    if fb > 0:
        b = fb / (FOREBRAIN_2B if sel else 1.0)
        kw.update(nmda_block=min(b, 1.0),
                  glun2b_selectivity=1.0 if sel else 0.0,
                  glun2b_fraction=PERIPHERAL_2B)
    return Drug(**kw)


def feasible(fb, sel):
    return fb / (FOREBRAIN_2B if sel else 1.0) <= 1.0


def dispatch(job):
    kind, pam, fb, sel, cap, seed = job
    d = make(pam, fb, sel, cap)
    if kind == "resp":
        b = PreBotC(drug=d, gaba_sens=GABA_SENS, seed=seed, **RESP_OP)
        for i in range(int(14000 / 0.1)):
            b.step(0.1)
            if i % 10 == 0:
                b.record()
        A = b.arrays(); m = A["t"] > 4000
        r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
        return job, r["mean"], r["amp"], float(r["alive"])
    orig = dict(SpinalCircuit.W)
    try:
        SpinalCircuit.W = dict(orig)
        for rec in ("ampa", "nmda"):
            SpinalCircuit.W[("Ia", "Mn", rec)] = orig[("Ia", "Mn", rec)] * IA_SCALE
        r = R.run(d, seed=seed)
        g = (r["mn_dyn"] - r["mn_base"]) / max(1e-6, r["ia_dyn"] - r["ia_base"])
        return job, g, 0.0, 1.0
    finally:
        SpinalCircuit.W = orig


def build():
    jobs, seen = [], set()
    for (lab, cap, eff) in GABA_ARMS:
        for pam, fb, sel in itertools.product(PAMS, FBS, (True, False)):
            if fb == 0.0 and not sel:
                continue
            if not feasible(fb, sel):
                continue
            key = (pam, fb, sel, cap)
            if key in seen:
                continue
            seen.add(key)
            for k in ("resp", "reflex"):
                jobs += [(k, pam, fb, sel, cap, s) for s in range(N_SEED)]
    return jobs


if __name__ == "__main__":
    jobs = build()
    print(f"{len(jobs)} sims on {os.cpu_count()} cores\n")
    with Pool(os.cpu_count()) as p:
        out = p.map(dispatch, jobs)
    A = {}
    for job, v, amp, alive in out:
        kind, pam, fb, sel, cap, seed = job
        A.setdefault((kind, pam, fb, sel, cap), []).append((v, amp, alive))
    cv = np.mean([x[0] for x in A[("resp", 1.0, 0.0, True, 2.5)]])
    cr = np.mean([x[0] for x in A[("reflex", 1.0, 0.0, True, 2.5)]])
    print(f"controls: ventilation {cv:.1f}, reflex gain {cr:.3f}")

    for eff_assume in (1.0, 1.3, 1.5, 1.8):
        print(f"\n{'='*78}\nsubjective efficiency of the neurosteroid arm = {eff_assume}x "
              f"(BZ-PAM fixed at 1.0x); subjective target {SUBJ_TARGET}")
        print(f"{'GABA arm':<13} {'PAM':>5} {'fb':>5} {'NMDA arm':>14} {'subj':>5} "
              f"{'vent%':>6} {'refl%':>6} {'apnoea':>7}")
        best = {}
        for (lab, cap, _) in GABA_ARMS:
            eff = eff_assume if lab == "neurosteroid" else 1.0
            cand = []
            for (kind, pam, fb, sel, c), v in A.items():
                if kind != "resp" or c != cap:
                    continue
                if ("reflex", pam, fb, sel, c) not in A:
                    continue
                gaba = eff * (min(pam, cap) - 1.0) / 1.5
                nmda = fb / 0.50
                if gaba < nmda - 1e-9:          # salience constraint
                    continue
                if gaba + nmda < SUBJ_TARGET - 1e-9:
                    continue
                vent = 100 * np.mean([x[0] for x in v]) / cv
                refl = 100 * np.mean([x[0] for x in A[("reflex", pam, fb, sel, c)]]) / cr
                apn = np.mean([x[2] for x in v]) < 0.5
                if apn:
                    continue
                cand.append((min(vent, refl), lab, pam, fb, sel, gaba + nmda,
                             vent, refl, apn))
            if cand:
                best[lab] = max(cand)
        for lab in ("BZ-PAM", "neurosteroid"):
            if lab not in best:
                print(f"{lab:<13} -- no feasible combination at this target")
                continue
            _, l, pam, fb, sel, subj, vent, refl, apn = best[lab]
            print(f"{l:<13} {pam:5.1f} {fb:5.0%} "
                  f"{'GluN2B-sel' if sel else 'non-selective':>14} {subj:5.2f} "
                  f"{vent:6.0f} {refl:6.0f} {'YES' if apn else 'no':>7}")
