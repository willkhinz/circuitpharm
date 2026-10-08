"""Verify the spinal CPG: does it produce alternating flexor/extensor rhythm?"""
import time
import numpy as np
from half_centre_cpg import HalfCentreCPG
from circuitpharm.cpg import Drug, burst_metrics
# SIDE-EFFECT GUARD added 2026-10-07. This script previously did its work at MODULE level,
# so merely importing it ran it. For `port_muscle.py` that regenerated the MuJoCo body
# model, and for `build_kb.py`/`build_compounds.py` it rewrote the pharmacology database --
# both destructive, both triggered by any tool that imports or scans the package. The body
# is unchanged; it is now reached only when the file is run directly.


def main():

    def run(drug=None, T=4000.0, dt=0.1, rec_every=10, **kw):
        c = HalfCentreCPG(drug=drug, seed=1, **kw)
        n = int(T / dt)
        for i in range(n):
            c.step(dt)
            if i % rec_every == 0: c.record()
        return c.arrays()

    t0 = time.time(); A = run(); el = time.time() - t0
    t = A["t"]; mask = t > 500           # discard transient
    bf = burst_metrics(t[mask], A["Mn_F"][mask]); be = burst_metrics(t[mask], A["Mn_E"][mask])
    print(f"drug-free: {el:.1f}s wall for 4.0 s sim ({4.0/el:.2f}x realtime)")
    print(f"  Mn_F mean={A['Mn_F'][mask].mean():6.2f} Hz  bursts={bf['n_bursts']:2d} "
          f"period={bf['period']:.0f} ms duty={bf['duty']:.2f}")
    print(f"  Mn_E mean={A['Mn_E'][mask].mean():6.2f} Hz  bursts={be['n_bursts']:2d} "
          f"period={be['period']:.0f} ms duty={be['duty']:.2f}")
    print(f"  RG_F mean={A['RG_F'][mask].mean():6.2f} Hz   RG_E mean={A['RG_E'][mask].mean():6.2f} Hz")
    f, e = A["Mn_F"][mask], A["Mn_E"][mask]
    if f.std() > 1e-6 and e.std() > 1e-6:
        print(f"  flexor/extensor correlation = {np.corrcoef(f,e)[0,1]:+.3f}  "
              f"(negative => alternation)")
    else:
        print("  NO RHYTHM (one or both pools silent/tonic)")


if __name__ == "__main__":
    main()
