"""THE deliverable: search the (GABA-A PAM, NMDA block, selectivity) space for a
combination that reaches an alcohol-like subjective intensity at acceptable motor and
respiratory cost.

CALIBRATED: preBotC gaba_sens = 0.15, fitted against human midazolam 2 mg IV
(ventilation -14 to -19%, tidal volume -16 to -22%, compensatory rate rise) and
cross-checked against propofol modelled as a non-ceiling-limited modulator (~-60% TV).

SUBJECTIVE AXIS -- an explicit, transparent assumption, not a simulation.
Drug discrimination says ethanol's stimulus is GABA-A positive modulation + NMDA
antagonism, and ethanol substituted only for mixtures where the GABA-A component had
EQUAL OR GREATER salience than the NMDA component. Operationalised as:

    gaba_term = (gaba_scale - 1) / 1.5          # PAM 2.5x (the ceiling) -> 1.0
    nmda_term = forebrain_NMDA_reduction / 0.50 # 50% reduction          -> 1.0
    subjective = gaba_term + nmda_term
    CONSTRAINT: gaba_term >= nmda_term          # the salience requirement

The normalisation constants are arbitrary in scale; what is NOT arbitrary is the
constraint's direction and the fact that the GABA arm must carry at least half the effect.
Conclusions are reported as a function of the subjective threshold so the arbitrary part
is visible rather than hidden.

HARM AXES (simulated): preBotC ventilation (minute-ventilation proxy) and stretch-reflex
gain, both with peripheral GluN2B fraction 0.15.
"""
import sys, os, itertools; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.circuit import SpinalCircuit
from circuitpharm.cpg import Drug
import scripts.reflex as R

FOREBRAIN_2B, PERIPHERAL_2B = 0.70, 0.15
GABA_SENS = 0.15
IA_SCALE = 0.30
N_SEED = 3
RESP_OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
               w=dict(ee_ampa=0.45, ee_nmda=0.2475))

PAM_LEVELS = (1.0, 1.4, 1.8, 2.2, 2.5)
FB_TARGETS = (0.0, 0.10, 0.20, 0.30, 0.40)
ARMS = (True, False)                      # GluN2B-selective, non-selective


def terms(pam, fb):
    gaba_scale = min(pam, 2.5)
    return (gaba_scale - 1.0) / 1.5, fb / 0.50


def make_drug(pam, fb, selective):
    if fb <= 0:
        return Drug(gaba_a_gain=pam, gaba_a_tau=1.0 + 0.6 * (pam - 1.0))
    b = fb / (FOREBRAIN_2B if selective else 1.0)
    return Drug(gaba_a_gain=pam, gaba_a_tau=1.0 + 0.6 * (pam - 1.0),
                nmda_block=min(b, 1.0),
                glun2b_selectivity=1.0 if selective else 0.0,
                glun2b_fraction=PERIPHERAL_2B)


def feasible(pam, fb, selective):
    b = fb / (FOREBRAIN_2B if selective else 1.0)
    return b <= 1.0


def _resp(drug, seed):
    b = PreBotC(drug=drug, gaba_sens=GABA_SENS, seed=seed, **RESP_OP)
    for i in range(int(14000 / 0.1)):
        b.step(0.1)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > 4000
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return r["mean"], r["amp"], float(r["alive"])


def _reflex(drug, seed):
    orig = dict(SpinalCircuit.W)
    try:
        SpinalCircuit.W = dict(orig)
        SpinalCircuit.W[("Ia", "Mn", "ampa")] = orig[("Ia", "Mn", "ampa")] * IA_SCALE
        SpinalCircuit.W[("Ia", "Mn", "nmda")] = orig[("Ia", "Mn", "nmda")] * IA_SCALE
        r = R.run(drug, seed=seed)
        return (r["mn_dyn"] - r["mn_base"]) / max(1e-6, r["ia_dyn"] - r["ia_base"])
    finally:
        SpinalCircuit.W = orig


def dispatch(job):
    kind, pam, fb, sel, seed = job
    d = make_drug(pam, fb, sel)
    if kind == "resp":
        mn, amp, alive = _resp(d, seed)
        return kind, pam, fb, sel, seed, mn, amp, alive
    return kind, pam, fb, sel, seed, _reflex(d, seed), 0.0, 1.0


def build():
    jobs = []
    for pam, fb, sel in itertools.product(PAM_LEVELS, FB_TARGETS, ARMS):
        if fb == 0.0 and not sel:
            continue                      # fb=0 is arm-independent; keep one copy
        if not feasible(pam, fb, sel):
            continue
        for k in ("resp", "reflex"):
            jobs += [(k, pam, fb, sel, s) for s in range(N_SEED)]
    return jobs


if __name__ == "__main__":
    jobs = build()
    print(f"{len(jobs)} sims on {os.cpu_count()} cores")
    with Pool(os.cpu_count()) as p:
        out = p.map(dispatch, jobs)
    A = {}
    for kind, pam, fb, sel, seed, v, amp, alive in out:
        A.setdefault((kind, pam, fb, sel), []).append((v, amp, alive))
    c_v = np.mean([x[0] for x in A[("resp", 1.0, 0.0, True)]])
    c_r = np.mean([x[0] for x in A[("reflex", 1.0, 0.0, True)]])
    print(f"controls: ventilation {c_v:.1f}, reflex gain {c_r:.3f}\n")

    rows = []
    for (kind, pam, fb, sel), v in A.items():
        if kind != "resp":
            continue
        key = (pam, fb, sel)
        if ("reflex", pam, fb, sel) not in A:
            continue
        g, n = terms(pam, fb)
        vent = 100 * np.mean([x[0] for x in v]) / c_v
        refl = 100 * np.mean([x[0] for x in A[("reflex", pam, fb, sel)]]) / c_r
        apn = np.mean([x[2] for x in v]) < 0.5
        rows.append(dict(pam=pam, fb=fb, sel=sel, gaba=g, nmda=n, subj=g + n,
                         vent=vent, refl=refl, apn=apn,
                         ok=(g >= n - 1e-9)))

    print(f"{'PAM':>5} {'fb':>5} {'arm':>5} {'gaba':>5} {'nmda':>5} {'subj':>5} "
          f"{'vent%':>6} {'refl%':>6} {'salience':>9}")
    print("-" * 64)
    for r in sorted(rows, key=lambda r: (-r["subj"], r["pam"])):
        print(f"{r['pam']:5.1f} {r['fb']:5.0%} {'2B' if r['sel'] else 'ns':>5} "
              f"{r['gaba']:5.2f} {r['nmda']:5.2f} {r['subj']:5.2f} "
              f"{r['vent']:6.0f} {r['refl']:6.0f} "
              f"{'OK' if r['ok'] else 'violates':>9}")

    print("\n=== best combination at each subjective threshold ===")
    print("(salience constraint satisfied; maximise the worse of the two harm axes)")
    for thr in (0.6, 0.8, 1.0, 1.2, 1.4):
        cand = [r for r in rows if r["ok"] and r["subj"] >= thr and not r["apn"]]
        if not cand:
            print(f"  subj>={thr:.1f}: NO FEASIBLE COMBINATION")
            continue
        best = max(cand, key=lambda r: min(r["vent"], r["refl"]))
        print(f"  subj>={thr:.1f}: PAM {best['pam']:.1f}x + {best['fb']:.0%} forebrain "
              f"NMDA ({'GluN2B-sel' if best['sel'] else 'non-selective'}) -> "
              f"ventilation {best['vent']:.0f}%, reflex {best['refl']:.0f}%")
