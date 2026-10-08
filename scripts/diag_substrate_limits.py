#!/usr/bin/env python3
"""What the LIF substrate does to the voltage-dependent mechanisms the pharmacology needs.

DIAGNOSTIC ONLY -- changes no model state and fits nothing. It exists because the roadmap's
stated reason for moving to conductance-based neurons was "LIF is not the biophysics," which
argues from principle, and this project has been wrong from principle before (the Matsuoka
locomotor oscillator, retracted). So: measure it.

Produces the table in knowledge/05-design-conductance-substrate.md section 1. Three
quantities, all of which a LIF truncates because V is reset at threshold and never reaches
spike voltages (recurring error E7):

  1. NMDA Mg2+ relief. cpg.py says "voltage dependence is why NMDA block is
     state-dependent." Measure whether this substrate can express that.
  2. GABA-A driving force (E_GABA - V). A PAM raises conductance; the CURRENT it buys is
     g*(E - V), so the driving force is what converts the kinetic finding into behaviour.
  3. Glutamate driving force, which in a real cell collapses at the burst peak and so
     self-limits excitation.

WHY THE RESULT IS NOT JUST "THE WEIGHTS ABSORBED IT". The weights were hand-tuned on this
substrate, so absolute magnitudes were compensated -- a scalar weight rescales a mean. What a
scalar weight cannot restore is a voltage-DEPENDENCE: it cannot widen a relief range, and it
cannot make one receptor pool's driving force swing while another's stays flat, because both
pools share a single weight-independent E_rev and a single V. So read the RANGES and RATIOS
below, not the means.

Run:  python scripts/diag_substrate_limits.py [--seconds 15]
"""
import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from circuitpharm.cpg import E_REV                      # noqa: E402
from circuitpharm.config import RESP_OP                 # noqa: E402
from circuitpharm.resp import PreBotC                   # noqa: E402

DT = 0.1            # ms, the integration step every circuit in this package uses
SETTLE_MS = 5000.0  # discard; the network starts at an arbitrary phase
SAMPLE_EVERY = 20   # thin the trace; 250k samples is already far past convergence

# Mg2+ block of the NMDA conductance, copied from cpg.Syn.conductance so this script
# measures the shipped expression rather than a re-derivation of it.
def mg_relief(V):
    return 1.0 / (1.0 + 0.28 * np.exp(-0.062 * V))


# The comparison column: what a spiking cell's voltage range makes available. Not a
# simulation -- the same algebra evaluated over a real cell's excursion.
SPIKING_RANGE = (-65.0, 20.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seconds", type=float, default=15.0)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    net = PreBotC(**RESP_OP, seed=a.seed)
    n_steps = int(a.seconds * 1000.0 / DT)
    V = []
    for i in range(n_steps):
        net.step(DT)
        if i * DT > SETTLE_MS and i % SAMPLE_EVERY == 0:
            V.append(net.pops["Exc"].V.copy())

    if not V:
        raise SystemExit(f"no samples: --seconds {a.seconds} is shorter than the "
                         f"{SETTLE_MS/1000.0:.0f} s settling period")
    V = np.concatenate(V)
    relief = mg_relief(V)
    df_gaba = np.abs(E_REV["gabaa"] - V)
    df_glut = np.abs(E_REV["ampa"] - V)

    lo, hi = SPIKING_RANGE
    r_lo, r_hi = mg_relief(lo), mg_relief(hi)

    print(f"\npreBotC Exc population, {a.seconds:.0f} s at dt={DT} ms, "
          f"{V.size} samples after settling\n")
    print(f"  V range               {V.min():+.2f} .. {V.max():+.2f} mV")
    print(f"  V p1/p50/p99          {np.percentile(V,1):+.2f} / "
          f"{np.percentile(V,50):+.2f} / {np.percentile(V,99):+.2f} mV")
    print(f"  (Vth = {net.pops['Exc'].Vth:+.0f} mV, so nothing above it but noise)\n")

    print("  NMDA Mg2+ relief      mean {:.4f}   range {:.4f} .. {:.4f}".format(
        relief.mean(), relief.min(), relief.max()))
    print(f"    dynamic range       {relief.max()/relief.min():.2f}x"
          f"   (a cell spanning {lo:+.0f}..{hi:+.0f} mV: {r_hi/r_lo:.1f}x)")
    print(f"    relief at {hi:+.0f} mV      {r_hi:.4f}   -> unreachable on this substrate\n")

    print(f"  GABA-A driving force  mean {df_gaba.mean():.1f} mV   max {df_gaba.max():.1f} mV")
    print(f"    at V = {hi:+.0f} mV        {abs(E_REV['gabaa']-hi):.1f} mV"
          f"   = {abs(E_REV['gabaa']-hi)/df_gaba.mean():.2f}x the measured mean")
    print(f"  glutamate driving f.  mean {df_glut.mean():.1f} mV   "
          f"min {df_glut.min():.1f} mV   (never collapses toward 0)\n")

    print("  Reading: the tonic GABA-A pool is always on at the RESTING driving force,")
    print("  which this substrate represents correctly. The phasic pool arrives correlated")
    print("  with the burst -- when a real cell would be depolarised and the driving force")
    print("  up to ~9x larger. That asymmetry is in the TIMING, so no scalar weight fixes")
    print("  it. See knowledge/05-design-conductance-substrate.md section 1.\n")


# __main__ GUARD IS LOAD-BEARING (recurring error E8, hit twice). Without it, any tool that
# imports or scans scripts/ runs this; and a heredoc-fed version of this file made __main__
# resolve to "<string>", so macOS spawn workers could not re-import it -- 13,334 tracebacks
# and 144 MB of output. Keep this a real file with this guard.
if __name__ == "__main__":
    main()
