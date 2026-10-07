"""THE project deliverable: motor and respiratory cost at MATCHED subjective effect.

Comparing at matched dose is meaningless -- a selective and a non-selective antagonist at
the same occupancy produce different forebrain effects. The right comparison is: fix the
forebrain (subjective) effect, then ask what each option costs on the two harm axes.

Subjective effect is NOT simulated. Per the project architecture it is an EMPIRICAL
CONSTRAINT from the drug-discrimination literature: ethanol's discriminative stimulus is a
compound stimulus = GABA-A positive modulation + NMDA antagonism, and ethanol substituted
only for mixtures where the GABA-A component had EQUAL OR GREATER salience than the NMDA
component.

Axes:
  subjective  (constraint) forebrain NMDA conductance reduction, GluN2B fraction 0.70
  motor       (simulated)  stretch-reflex gain, peripheral GluN2B fraction 0.15
  respiratory (simulated)  preBotC mean output = minute-ventilation proxy, 2B 0.15

WEAKEST LINK, stated plainly: the subjective axis is a linear receptor-level index, not a
circuit model. Drug discrimination requires learning and cannot be simulated with a frozen
policy. All numbers below are therefore "cost at a *nominal* matched subjective effect".

NOTE (bug fixed): every worker function must be at MODULE level. macOS multiprocessing
uses spawn, so a function defined inside `if __name__ == "__main__"` is not picklable.
"""
import sys, os; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.circuit import SpinalCircuit
from circuitpharm.cpg import Drug
import scripts.reflex as R

FOREBRAIN_2B, PERIPHERAL_2B = 0.70, 0.15
PAM_GAIN, PAM_TAU = 2.0, 1.6           # fixed GABA-A component in every arm
IA_SCALE = 0.30                         # reflex operating point (off the ceiling)
TARGETS = (0.15, 0.25, 0.35, 0.45)      # forebrain NMDA conductance reduction
N_SEED = 4
RESP_OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
               w=dict(ee_ampa=0.45, ee_nmda=0.2475))


def block_for(target, selective):
    """nmda_block needed to achieve `target` forebrain NMDA reduction."""
    accessible = FOREBRAIN_2B if selective else 1.0
    b = target / accessible
    return b, b <= 1.0


def drug_for(target, selective):
    if target <= 0.0:
        return Drug()
    b, _ = block_for(target, selective)
    return Drug(gaba_a_gain=PAM_GAIN, gaba_a_tau=PAM_TAU,
                nmda_block=min(b, 1.0),
                glun2b_selectivity=1.0 if selective else 0.0,
                glun2b_fraction=PERIPHERAL_2B)


def _resp(drug, seed):
    b = PreBotC(drug=drug, seed=seed, **RESP_OP)
    for i in range(int(14000 / 0.1)):
        b.step(0.1)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > 4000
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return r["mean"], float(r["alive"])


def _reflex(drug, seed):
    orig = dict(SpinalCircuit.W)
    try:
        SpinalCircuit.W = dict(orig)
        SpinalCircuit.W[("Ia", "Mn", "ampa")] = orig[("Ia", "Mn", "ampa")] * IA_SCALE
        SpinalCircuit.W[("Ia", "Mn", "nmda")] = orig[("Ia", "Mn", "nmda")] * IA_SCALE
        r = R.run(drug, seed=seed)
        g = (r["mn_dyn"] - r["mn_base"]) / max(1e-6, r["ia_dyn"] - r["ia_base"])
        return g, 1.0
    finally:
        SpinalCircuit.W = orig


def dispatch(job):
    """MODULE-LEVEL so it is picklable under spawn."""
    kind, target, selective, seed = job
    drug = drug_for(target, selective)
    v, alive = _resp(drug, seed) if kind == "resp" else _reflex(drug, seed)
    return kind, target, selective, seed, v, alive


def build_jobs():
    jobs = [(k, 0.0, False, s) for k in ("resp", "reflex") for s in range(N_SEED)]
    for t in TARGETS:
        for sel in (True, False):
            for k in ("resp", "reflex"):
                jobs += [(k, t, sel, s) for s in range(N_SEED)]
    return jobs


if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        out = p.map(dispatch, build_jobs())
    agg = {}
    for kind, t, sel, s, v, alive in out:
        agg.setdefault((kind, t, sel), []).append((v, alive))
    c_resp = np.mean([a[0] for a in agg[("resp", 0.0, False)]])
    c_refl = np.mean([a[0] for a in agg[("reflex", 0.0, False)]])
    print(f"controls: ventilation {c_resp:.1f} Hz, reflex gain {c_refl:.3f}")
    print(f"GABA-A component fixed at PAM {PAM_GAIN}x (ceiling 2.5) in every arm")
    print(f"forebrain GluN2B fraction {FOREBRAIN_2B}, peripheral {PERIPHERAL_2B}\n")
    print(f"{'forebrain':>9} {'arm':>14} {'nmda':>6} {'periph':>7} | "
          f"{'vent %':>8} {'reflex %':>9} {'apnoea':>7}")
    print(f"{'target':>9} {'':>14} {'block':>6} {'cut':>7} |")
    print("-" * 72)
    rows = []
    for t in TARGETS:
        for sel in (True, False):
            b, ok = block_for(t, sel)
            periph = min(b, 1.0) * (PERIPHERAL_2B if sel else 1.0)
            rv = np.array([a[0] for a in agg[("resp", t, sel)]])
            al = np.array([a[1] for a in agg[("resp", t, sel)]])
            rf = np.array([a[0] for a in agg[("reflex", t, sel)]])
            vent = 100 * rv.mean() / c_resp
            refl = 100 * rf.mean() / c_refl
            arm = "GluN2B-sel" if sel else "non-selective"
            print(f"{t:9.0%} {arm:>14} {min(b,1.0):6.2f} {periph:7.1%} | "
                  f"{vent:8.0f} {refl:9.0f} {'YES' if al.mean()<0.5 else 'no':>7}"
                  f"{'' if ok else '  (block>1, infeasible)'}")
            rows.append((t, sel, vent, refl, ok))
    print("\nselective advantage at matched subjective effect:")
    for t in TARGETS:
        s = next(r for r in rows if r[0] == t and r[1])
        n = next(r for r in rows if r[0] == t and not r[1])
        print(f"  {t:.0%} forebrain: ventilation {s[2]:.0f}% vs {n[2]:.0f}% "
              f"(+{s[2]-n[2]:.0f} pts), reflex {s[3]:.0f}% vs {n[3]:.0f}% "
              f"(+{s[3]-n[3]:.0f} pts)" + ("" if s[4] else "   [selective infeasible]"))
