#!/usr/bin/env python3
"""A4: does the NETWORK rhythm survive removal of I_NaP, when the isolated cell's does not?

ROADMAP LINK 5, and the one genuinely falsifiable prediction this upgrade makes that needs no
new parameters from anybody's paper.

THE PREDICTION. In vitro, blocking the persistent sodium current (riluzole) abolishes
pacemaking in isolated preBotC neurons while the network rhythm persists -- the
pacemaker-versus-network-rhythm controversy. So a model whose cells are real pacemakers should
show a DISSOCIATION:

    isolated cell, g_NaP -> 0   :  bursting must STOP     (already shown, link 4)
    network,       g_NaP -> 0   :  rhythm should PERSIST

The LIF model cannot pose this test at all. It has no I_NaP to block: its "I_NaP" is an
instantaneous m-gate with no inactivation, and burst termination is a spike-triggered
adaptation conductance standing in for it. Removing that stand-in removes the rhythm by
construction, which would tell us nothing about the biology.

WHAT COUNTS AS A RESULT, PRE-COMMITTED. This script MEASURES; it does not assert. Either
outcome is informative and both get reported:

  * rhythm persists  -> the network is a network oscillator, matching the preparation, and
    the conductance substrate has bought a published, independently-testable behaviour the
    LIF substrate could not express.
  * rhythm dies      -> the model is a pacemaker-driven rhythm, which does NOT match the
    preparation. That is a real limitation of the network and must be recorded as one. It
    would NOT be a reason to retune until the rhythm survives: that would be fitting the
    model to a prediction and then claiming the prediction as validation, which is exactly
    how this project's retracted Matsuoka result happened.

Run:  python scripts/probe_inap_dissociation.py [--seconds 20] [--seeds 3]
"""
import argparse
import pathlib
import sys
from dataclasses import replace

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from circuitpharm.neuron import BRS1999_MODEL1, run_isolated   # noqa: E402
from circuitpharm.resp import PreBotC, resp_metrics            # noqa: E402

NO_NAP = replace(BRS1999_MODEL1, g_nap=0.0)


def network(cell, seconds, seed):
    net = PreBotC(substrate="cond", cell=cell, seed=seed)
    for i in range(int(round(seconds * 1000.0 / 0.1))):
        net.step(0.1)
        if i % 10 == 0:
            net.record()
    A = net.arrays()
    m = A["t"] > 4000.0                      # discard the settling transient
    return resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seconds", type=float, default=20.0)
    ap.add_argument("--seeds", type=int, default=3)
    a = ap.parse_args()

    print("=" * 78)
    print("ISOLATED CELL  (expected: bursting with I_NaP, silent without)")
    print("=" * 78)
    for label, cell in (("g_NaP = 2.8 nS", BRS1999_MODEL1), ("g_NaP = 0", NO_NAP)):
        r = run_isolated(seconds=a.seconds, params=cell)
        print(f"  {label:16s} {r['regime']:10s} bursts={r['n_bursts']:3d} "
              f"rate={r['rate_hz']:6.2f} Hz  h-span={r['h_span']:.4f}")

    print()
    print("=" * 78)
    print("NETWORK  (the prediction: rhythm PERSISTS without I_NaP)")
    print("=" * 78)
    out = {}
    for label, cell in (("g_NaP = 2.8 nS", BRS1999_MODEL1), ("g_NaP = 0", NO_NAP)):
        rows = [network(cell, a.seconds, s) for s in range(a.seeds)]
        alive = sum(r["alive"] for r in rows)
        fq = float(np.mean([r["freq"] for r in rows]))
        mn = float(np.mean([r["mean"] for r in rows]))
        md = float(np.mean([r["mod"] for r in rows]))
        out[label] = dict(alive=alive, freq=fq, mean=mn, mod=md,
                          reasons=[r["reason"] for r in rows])
        print(f"  {label:16s} alive {alive}/{a.seeds}   freq {fq:6.3f} Hz   "
              f"mean {mn:6.2f}   mod {md:5.2f}")
        for i, r in enumerate(rows):
            print(f"       seed {i}: {r['reason']}")

    print()
    print("=" * 78)
    intact, blocked = out["g_NaP = 2.8 nS"], out["g_NaP = 0"]
    if not intact["alive"]:
        print("INCONCLUSIVE: the control network is not alive, so there is no rhythm whose")
        print("survival could be tested. Fix the operating point first.")
    elif blocked["alive"]:
        print("DISSOCIATION REPRODUCED. The isolated cell stops bursting without I_NaP while")
        print("the network rhythm persists -- the published in vitro behaviour, on a test the")
        print("LIF substrate cannot pose. Frequency moved "
              f"{intact['freq']:.3f} -> {blocked['freq']:.3f} Hz "
              f"({100*(blocked['freq']-intact['freq'])/intact['freq']:+.1f}%).")
    else:
        print("DISSOCIATION NOT REPRODUCED. Removing I_NaP abolishes the NETWORK rhythm too,")
        print("so this network is pacemaker-driven and does not match the preparation, where")
        print("riluzole leaves the rhythm intact. Record as a limitation of the network.")
        print("Do NOT retune the coupling until the rhythm survives: fitting the model to a")
        print("prediction and then claiming the prediction as validation is how the Matsuoka")
        print("locomotor result had to be retracted.")
    print("=" * 78)


# __main__ guard: E8 bit this project twice, and the second form was a file whose __main__
# resolved to "<string>" so spawn workers could not re-import it.
if __name__ == "__main__":
    main()
