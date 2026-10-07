"""Find an Ia->Mn gain that puts the dynamic reflex OFF the refractory ceiling,
so drug modulation is visible. Target dynamic peak 40-70 Hz (ceiling is 125 Hz)."""
import sys, os; sys.path.insert(0,".")
import numpy as np
from multiprocessing import Pool
from circuitpharm.circuit import SpinalCircuit
from circuitpharm.cpg import Drug
import scripts.reflex as R

def probe(scale):
    orig = dict(SpinalCircuit.W)
    try:
        SpinalCircuit.W = dict(orig)
        SpinalCircuit.W[("Ia","Mn","ampa")] = orig[("Ia","Mn","ampa")]*scale
        SpinalCircuit.W[("Ia","Mn","nmda")] = orig[("Ia","Mn","nmda")]*scale
        c = R.run(Drug())
        g = R.run(Drug(gaba_a_gain=2.0, gaba_a_tau=1.6))
        n = R.run(Drug(nmda_block=0.6, glun2b_selectivity=0.0))
        s = R.run(Drug(glyr_gain=0.4))
        return dict(scale=scale, dyn=c['mn_dyn'], sta=c['mn_sta'], base=c['mn_base'],
                    pam=100*g['mn_dyn']/max(1e-9,c['mn_dyn']),
                    nmda=100*n['mn_dyn']/max(1e-9,c['mn_dyn']),
                    stry=100*s['mn_dyn']/max(1e-9,c['mn_dyn']))
    finally:
        SpinalCircuit.W = orig

if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        res = p.map(probe, [0.06,0.10,0.15,0.22,0.32,0.50,1.0])
    print(f"{'scale':>6} {'Mn dyn':>7} {'Mn sta':>7} {'Mn base':>8} | "
          f"{'PAM %':>7} {'NMDA %':>7} {'stry %':>7}   (ceiling = 125 Hz)")
    for r in res:
        print(f"{r['scale']:6.2f} {r['dyn']:7.1f} {r['sta']:7.1f} {r['base']:8.1f} | "
              f"{r['pam']:7.0f} {r['nmda']:7.0f} {r['stry']:7.0f}")
