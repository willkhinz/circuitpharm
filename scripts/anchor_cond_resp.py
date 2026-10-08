#!/usr/bin/env python3
"""Anchor the conductance-substrate preBotC to the SAME OBSERVABLE as the LIF one.

ROADMAP LINK 5. `PreBotC(substrate="cond")` refuses to run without an operating point,
because `drive` is in pA and every weight in nS, both relative to the LIF cell (C=200 pF,
g_L=10 nS) rather than the Butera cell (C=21 pF, g_L=2.8 nS). This script finds that
operating point; the result is written into `config.COND_RESP_OP`.

WHY ANCHORED AND NOT INHERITED. Butera-Rinzel-Smith part II supplies network coupling
conductances for a population of exactly these cells, which would have made these weights
published rather than ours -- the actual point of link 5. They are not retrievable: the
journal full text returns HTTP 403 and every accessible encoding (the curated CellML, ModelDB
247647) is single-cell only. Inventing them was not an option. So the CELL is published and
the COUPLING is ours, and link 5 is partial for that reason. Said plainly here because a
reader could otherwise reasonably assume the whole network came from the paper.

MATCHED OPERATING POINT, NOT MATCHED PARAMETERS. This is the part that makes the
substrate comparison mean anything. The same nS weight does different things on the two
cells, so comparing substrates at equal weights compares two differently-broken networks.
Each substrate is anchored independently to the same observable -- control burst frequency
inside EUPNOEA_BAND, with a comparable modulation depth -- and only then is the drug applied
and the fractional change from each substrate's own control compared.

The target is the LIF's measured control, not a literature value: freq 1.271 Hz, mod 4.74,
mean 29.3 (RESP_OP, seed 0). That is deliberate. The question this enables is "does the same
drug produce the same FRACTIONAL change on a different substrate", which needs the two
controls to agree on the observable, and says nothing about either being the right absolute
frequency. The respiratory axis has no valid absolute anchor -- that is unchanged by this
script and is a problem about the observable, not the neurons.

Run:  python scripts/anchor_cond_resp.py [--seconds 8] [--procs 8] [--stage 1|2]
"""
import argparse
import itertools
import multiprocessing as mp
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from circuitpharm.config import RESP_OP                      # noqa: E402
from circuitpharm.resp import PreBotC, REQUIRED_W, resp_metrics  # noqa: E402

# The LIF control, measured. Re-measured by --verify rather than trusted.
TARGET = dict(freq=1.271, mod=4.736, mean=29.303)

# The LIF weight table, as PreBotC assembles it. The search scales this as a block, so the
# RATIOS between pathways are held at the LIF's and only the overall scale moves. Holding the
# ratios is a choice, not a necessity: it keeps the search small and means the conductance
# network is the same circuit at a different scale rather than a differently-wired one, so a
# substrate difference cannot be a wiring difference in disguise.
LIF_W = dict(RESP_OP["w"])
LIF_W.update(ei_ampa=0.55, ei_nmda=0.0, ie_gaba=0.45, ie_gly=0.35,
             eo_ampa=0.70, eo_nmda=0.30)


def make_op(scale, drive, drive_other, gaba_tonic):
    return dict(drive=float(drive), drive_other=float(drive_other),
                gaba_tonic=float(gaba_tonic),
                w={k: LIF_W[k] * scale for k in REQUIRED_W})


# tau_h is 10 s, the slowest timescale in the Butera cell and 25x the LIF's tau_adapt.
# A warm-up shorter than a few multiples of it measures a transient. The first version of
# this script used 4 s and produced a confident, verified, wrong operating point.
SETTLE_MS = 30000.0          # 3 x tau_h
WINDOW_MS = 15000.0          # each of two consecutive late windows


def evaluate(args):
    """One operating point, measured OUTSIDE the transient and checked for drift.

    Two consecutive late windows are measured separately. Both must be alive and their
    frequencies must AGREE, because a network whose rhythm is still decaying can read
    perfectly well over any single window -- which is exactly how the first anchoring
    attempt passed a three-stage search with a seed-robustness check.
    """
    scale, drive, drive_other, gaba_tonic, _seconds, seed = args
    op = make_op(scale, drive, drive_other, gaba_tonic)
    total_ms = SETTLE_MS + 2.0 * WINDOW_MS
    base = dict(scale=scale, drive=drive, drive_other=drive_other, gaba_tonic=gaba_tonic)
    try:
        net = PreBotC(substrate="cond", op=op, seed=seed)
        for i in range(int(round(total_ms / 0.1))):
            net.step(0.1)
            if i % 10 == 0:
                net.record()
        A = net.arrays()
        halves = []
        for k in (0, 1):
            lo = SETTLE_MS + k * WINDOW_MS
            m = (A["t"] >= lo) & (A["t"] < lo + WINDOW_MS)
            halves.append(resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m]))
    except Exception as exc:                      # a diverged cell must not kill the sweep
        return dict(**base, score=float("inf"), freq=float("nan"), mod=float("nan"),
                    mean=float("nan"), alive=False, drift=float("nan"),
                    reason=f"{type(exc).__name__}: {exc}"[:50])

    f1, f2 = halves[0]["freq"], halves[1]["freq"]
    alive = bool(halves[0]["alive"] and halves[1]["alive"])
    freq = float(np.mean([f1, f2]))
    mod = float(np.mean([h["mod"] for h in halves]))
    mean = float(np.mean([h["mean"] for h in halves]))
    drift = abs(f1 - f2) / max(freq, 1e-9)

    # DRIFT IS A HARD GATE, not a penalty term. A penalty lets a drifting point win a
    # sparse grid on the strength of its frequency being close on average.
    if not alive:
        score, reason = float("inf"), (halves[0]["reason"] if not halves[0]["alive"]
                                       else halves[1]["reason"])[:50]
    elif drift > 0.20:
        score, reason = float("inf"), f"drifting: {f1:.2f} then {f2:.2f} Hz"
    else:
        score = (abs(freq - TARGET["freq"]) / TARGET["freq"]
                 + 0.20 * abs(np.log(max(mod, 1e-6) / TARGET["mod"])))
        reason = "alive, stable"
    return dict(**base, score=float(score), freq=freq, mod=mod, mean=mean,
                alive=alive, drift=float(drift), reason=reason)


def sweep(grid, seconds, procs, seed=0):
    jobs = [(s, d, do, gt, seconds, seed) for s, d, do, gt in grid]
    if procs > 1:
        with mp.Pool(procs) as pool:
            return pool.map(evaluate, jobs)
    return [evaluate(j) for j in jobs]


def report(rows, n=12):
    rows = sorted(rows, key=lambda r: r["score"])
    print(f"\n{'scale':>6}{'drive':>7}{'d_oth':>7}{'gtonic':>7} | "
          f"{'freq':>7}{'mod':>7}{'mean':>8}{'drift':>7} {'alive':>6}{'score':>8}  reason")
    print("-" * 104)
    for r in rows[:n]:
        sc = "inf" if not np.isfinite(r["score"]) else f"{r['score']:.4f}"
        print(f"{r['scale']:6.2f}{r['drive']:7.1f}{r['drive_other']:7.1f}"
              f"{r['gaba_tonic']:7.2f} | {r['freq']:7.3f}{r['mod']:7.2f}{r['mean']:8.1f}"
              f"{100*r['drift']:6.1f}% {str(r['alive']):>6}{sc:>8}  {r['reason']}")
    alive = [r for r in rows if r["alive"]]
    print(f"\n{len(alive)}/{len(rows)} points produced a living rhythm")
    return rows[0] if rows and np.isfinite(rows[0]["score"]) else None


# Stage-2 finalists, carried forward verbatim from its output. Short sweeps are not
# trustworthy for frequency: the stage-1 winner read 1.026 Hz at 8 s and 1.102 Hz at 12 s,
# a 7% move from duration alone. So the winner of a sweep is a CANDIDATE, and stage 3 is
# what decides -- longer runs, several seeds, and a requirement that every seed be alive.
# Picking on a single short run would be choosing an operating point partly on seed noise.
FINALISTS = [
    #  scale  drive  d_other  gaba_tonic      (stage-2 12 s reading)
    (0.90,   20.0,   30.0,    1.20),          # 1.356 Hz, mod 6.29
    (1.00,   20.0,   20.0,    1.50),          # 1.102 Hz, mod 5.87
    (1.00,   15.0,   15.0,    1.20),          # 1.017 Hz, mod 4.39
    (1.10,   20.0,   30.0,    1.20),          # 1.441 Hz, mod 3.50
]


def verify(a):
    """Re-measure the finalists at length, across seeds. Every seed must be alive."""
    print(f"stage 3: {len(FINALISTS)} finalists x {a.seeds} seeds, {a.seconds:.0f} s each")
    print(f"target (the LIF's measured control): freq {TARGET['freq']} Hz, "
          f"mod {TARGET['mod']}\n")
    print(f"{'scale':>6}{'drive':>7}{'d_oth':>7}{'gtonic':>7} | {'alive':>7}"
          f"{'freq':>9}{'sd':>7}{'mod':>7}{'mean':>8}{'|dfreq|':>9}")
    print("-" * 82)
    best = None
    for scale, drive, do, gt in FINALISTS:
        jobs = [(scale, drive, do, gt, a.seconds, s) for s in range(a.seeds)]
        if a.procs > 1:
            with mp.Pool(min(a.procs, len(jobs))) as pool:
                rows = pool.map(evaluate, jobs)
        else:
            rows = [evaluate(j) for j in jobs]
        n_alive = sum(r["alive"] for r in rows)
        fq = np.array([r["freq"] for r in rows], float)
        md = float(np.mean([r["mod"] for r in rows]))
        mn = float(np.mean([r["mean"] for r in rows]))
        dfreq = abs(float(fq.mean()) - TARGET["freq"]) / TARGET["freq"]
        print(f"{scale:6.2f}{drive:7.1f}{do:7.1f}{gt:7.2f} | {n_alive:4d}/{a.seeds}"
              f"{fq.mean():9.3f}{fq.std():7.3f}{md:7.2f}{mn:8.1f}{100*dfreq:8.1f}%")
        # EVERY SEED MUST BE ALIVE. An operating point that is alive in 3 of 4 seeds is
        # not a usable control: the drug arms would then be compared against a baseline
        # that sometimes does not exist, and a drug would get credit for a seed that was
        # already dead.
        if n_alive == a.seeds and (best is None or dfreq < best[0]):
            best = (dfreq, scale, drive, do, gt, float(fq.mean()), md, mn)

    if best is None:
        print("\nNO finalist was alive in every seed. Do not register any of them; widen "
              "the stage-2 grid instead. Do NOT lower the seed requirement.")
        return
    dfreq, scale, drive, do, gt, fq, md, mn = best
    print("\nREGISTER THIS:\n")
    print("COND_RESP_OP = MappingProxyType(dict(")
    print(f"    drive={drive}, drive_other={do}, gaba_tonic={gt},")
    print("    w=MappingProxyType(dict("
          + ", ".join(f"{k}={LIF_W[k]*scale:.6g}" for k in REQUIRED_W) + ")),")
    print("))")
    print(f"\n# weight-block scale {scale} x the LIF table; control {fq:.3f} Hz vs the "
          f"LIF's {TARGET['freq']} ({100*dfreq:.1f}% off), mod {md:.2f}, mean {mn:.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seconds", type=float, default=8.0)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--stage", type=int, default=1, choices=(1, 2, 3))
    ap.add_argument("--seeds", type=int, default=4,
                    help="stage 3 only: how many seeds to verify each candidate across")
    a = ap.parse_args()

    if a.stage == 3:
        return verify(a)
    if a.stage == 1:
        # COARSE. Scale the whole weight block and move the two drives. gaba_tonic tracks
        # the scale, since it is an inhibitory conductance in the same units.
        # EXC DRIVE IS SEARCHED WHERE THE ISOLATED CELL BURSTS: -5..0 pA, verified with
        # neuron.run_isolated (quiescent at -10, bursting at -5 and 0, TONIC from +5 up).
        # The first version swept {5..30} pA -- entirely inside the tonic regime -- which is
        # why it found almost nothing alive and what it did find was a transient.
        #
        # The range extends below -5 because network cells also receive recurrent AMPA and
        # NMDA excitation on top of `drive`, so the external drive that puts a COUPLED cell
        # in the bursting regime is lower than for an isolated one.
        grid = [(s, d, do, 1.5 * s)
                for s, d, do in itertools.product(
                    (0.5, 0.75, 1.0, 1.5),            # weight-block scale
                    (-12.0, -8.0, -5.0, -2.0, 0.0, 3.0),   # Exc drive, pA
                    (10.0, 25.0))]                    # Inh/Out drive, pA (no I_NaP there)
    else:
        # REFINE around the stage-1 winner, which was scale 1.00 / drive 20 / drive_other 20
        # / gaba_tonic 1.50, giving 1.026 Hz against the 1.271 Hz target.
        #
        # Stage 1 found only 3 of 50 points alive, all at scale 1.00. That narrowness is a
        # result, not a search artifact: a population of intrinsically bursting Butera cells
        # has a much smaller region of synchronised, in-band behaviour than the LIF network
        # did, because each cell is already an oscillator and the coupling has to entrain
        # rather than create the rhythm. Reported rather than widened away.
        grid = [(s, d, do, gt)
                for s, d, do, gt in itertools.product(
                    (0.9, 1.0, 1.1),
                    (15.0, 20.0, 25.0),
                    (15.0, 20.0, 30.0),
                    (1.2, 1.5, 1.8))]

    print(f"stage {a.stage}: {len(grid)} operating points, {a.seconds:.0f} s each, "
          f"{a.procs} procs")
    print(f"target (the LIF's measured control): freq {TARGET['freq']} Hz, "
          f"mod {TARGET['mod']}, mean {TARGET['mean']}")
    best = report(sweep(grid, a.seconds, a.procs))
    if best is None:
        print("\nNO operating point produced a living rhythm in band. Widen the grid; do "
              "NOT relax EUPNOEA_BAND to make a point qualify.")
        return
    print("\nBEST:")
    print(f"  COND_RESP_OP = dict(")
    print(f"      drive={best['drive']}, drive_other={best['drive_other']},")
    print(f"      gaba_tonic={best['gaba_tonic']},")
    print(f"      w={{" + ", ".join(f"{k!r}: {LIF_W[k]*best['scale']:.6g}"
                                    for k in REQUIRED_W) + "},")
    print(f"  )   # freq {best['freq']:.3f} Hz vs LIF {TARGET['freq']}, "
          f"mod {best['mod']:.2f} vs {TARGET['mod']}, scale {best['scale']}")
    print("\nVerify at a longer duration and several seeds before registering it.")


# __main__ GUARD IS LOAD-BEARING. E8 bit this project twice: once because scripts did
# destructive work at import, and once because a heredoc-fed file made __main__ resolve to
# "<string>", so macOS spawn workers could not re-import it -- 13,334 tracebacks and 144 MB
# of output. multiprocessing here makes that failure mode live again.
if __name__ == "__main__":
    main()
