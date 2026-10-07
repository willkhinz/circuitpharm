"""Is synaptic inhibition load-bearing for preBotC burst amplitude at all?

WHY. Adding an NMDA component to Exc->Inh (scripts/calib_ei_nmda.py) did NOT change the
sign of the ventilation response to NMDA block -- ventilation still fell to ~40% at 90%
block regardless of the share. The hypothesis for why: in this circuit bursts are
TERMINATED BY SPIKE-TRIGGERED ADAPTATION, not by synaptic inhibition, so disinhibition has
nothing to release. If that is right, removing inhibition entirely should barely move
ventilation.

Per recurring error E10, a coupling is tested by REMOVING IT COMPLETELY, never by partial
block -- a 95% reduction previously left entrainment fully intact elsewhere in this project.

Four conditions, all drug-free unless stated:
  intact            reference
  no inhibition     ie_gaba = ie_gly = 0        (is inhibition load-bearing?)
  no adaptation     g_adapt = 0                 (is adaptation what terminates bursts?)
  no inh, NMDA 90%  does removing inhibition rescue NMDA block?

Run:  python scripts/diag_inhib_load.py
"""
import sys, os; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from spinal.resp import PreBotC, resp_metrics
from spinal.cpg import Drug

EE_AMPA, EE_NMDA = 0.45, 0.2475
N_SEED = 4
DUR_MS, DT = 14000.0, 0.1

CONDS = {
    "intact (reference)":        dict(),
    "inhibition REMOVED":        dict(w=dict(ie_gaba=0.0, ie_gly=0.0)),
    "adaptation REMOVED":        dict(g_adapt=0.0),
    "NMDA 90% block":            dict(blk=0.9),
    "NMDA 90% + inhib REMOVED":  dict(blk=0.9, w=dict(ie_gaba=0.0, ie_gly=0.0)),
}


def job(a):
    label, seed = a
    cfg = dict(CONDS[label])
    blk = cfg.pop("blk", 0.0)
    wextra = cfg.pop("w", {})
    g_adapt = cfg.pop("g_adapt", 2.5)
    w = dict(ee_ampa=EE_AMPA, ee_nmda=EE_NMDA); w.update(wextra)
    b = PreBotC(drug=Drug(nmda_block=blk, glun2b_selectivity=0.0),
                gaba_sens=0.15, seed=seed,
                drive=170.0, g_adapt=g_adapt, tau_adapt=400.0, w=w)
    for i in range(int(DUR_MS / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > 4000
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return label, r["mean"], r["mod"], r["freq"], float(r["alive"])


if __name__ == "__main__":
    jobs = [(k, s) for k in CONDS for s in range(N_SEED)]
    print(f"{len(jobs)} sims on {os.cpu_count()} cores")
    with Pool(os.cpu_count()) as p:
        out = p.map(job, jobs)
    A = {}
    for lab, mn, md, fq, al in out:
        A.setdefault(lab, []).append((mn, md, fq, al))
    ref = np.mean([x[0] for x in A["intact (reference)"]])
    print(f"\n{'condition':<28}{'vent %':>9}{'mod':>8}{'freq Hz':>9}{'alive':>7}")
    print("-" * 61)
    for k in CONDS:
        mn = np.mean([x[0] for x in A[k]])
        md = np.mean([x[1] for x in A[k]])
        fq = np.mean([x[2] for x in A[k]])
        al = np.mean([x[3] for x in A[k]])
        print(f"{k:<28}{100*mn/ref:8.0f}%{md:8.2f}{fq:9.2f}{al:7.2f}")
    print("\nRESULT (2026-10-07), and it REFUTED the hypothesis this script was written")
    print("to test. The hypothesis was that inhibition is not load-bearing, leaving")
    print("nothing for disinhibition to release. The data says otherwise:")
    print("  inhibition REMOVED        -> 167%. Inhibition IS load-bearing; it holds")
    print("                               ventilation to ~60% of its uninhibited value,")
    print("                               so there IS inhibition available to release.")
    print("  adaptation REMOVED        -> 689%, mod 0.00, not alive. Tonic firing, no")
    print("                               rhythm: adaptation is what terminates bursts.")
    print("  NMDA 90%                  -> 45%")
    print("  NMDA 90% + inhib REMOVED  -> 55%")
    print("\nSo the conclusion holds for a STRONGER reason than hypothesised. Removing")
    print("100% of synaptic inhibition -- far more disinhibition than any drug could")
    print("produce -- recovers only ~10 points of a ~55 point deficit, whereas the same")
    print("removal is worth 67 points drug-free. The loss of the long-tau NMDA")
    print("conductance degrades the CHARACTER of the rhythm (note freq 4.06 Hz, i.e.")
    print("fragmented fast bursting, not physiological eupnoea), and disinhibition")
    print("cannot restore synchronous bursting it can only amplify what remains.")
    print("\nCAVEAT ON THE FREQ COLUMN: 4.06 Hz at 90% block is ~240 breaths/min, which")
    print("is not physiological for a depressed rat. The ventilation conclusion uses")
    print("'mean', not 'freq', so it stands -- but the frequency readout under severe")
    print("rhythm degradation is not trustworthy and should not be quoted. Related to")
    print("recurring error E4, in a form the FFT fix does not fully solve.")
