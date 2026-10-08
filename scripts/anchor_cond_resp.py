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

from circuitpharm.config import INVITRO_BAND, RESP_OP        # noqa: E402
from circuitpharm.resp import PreBotC, REQUIRED_W, resp_metrics  # noqa: E402

# The LIF control, measured. Kept for reporting the contrast, NOT as a target any more.
TARGET = dict(freq=1.271, mod=4.736, mean=29.303)

# THE CONDUCTANCE ARM IS ANCHORED TO ITS OWN PREPARATION, not to the LIF's frequency.
#
# Three attempts tried to force it to 1.271 Hz and all three failed the same way: the only
# stable, deeply modulated points (drift 0%, modulation 2.5-9.9) sat at 0.27-0.34 Hz, while
# every point reaching the in vivo band had modulation collapsed to 0.2-0.9 and drift of
# 40-164%. The target was wrong, not the network.
#
# EUPNOEA_BAND is an IN VIVO rat band (1-2 Hz). The Butera-Rinzel-Smith cell is neonatal
# rodent IN VITRO, where control inspiratory burst frequency is 6.6 +/- 3.1 to 14.6 +/- 2.0
# bursts/min = ~0.11-0.24 Hz (source pbc_invitro_freq: Revill et al. 2021,
# Front Physiol 12:626470). Forcing it 5-12x above that is what destroyed the rhythm.
#
# So frequency is now a RANGE the rhythm must fall in, taken from the source, and the thing
# being optimised is rhythm QUALITY -- modulation depth and stationarity. The two substrates
# therefore do not share a frequency, which is correct: only fractional change from each
# substrate's own control is comparable between them (design doc section 6).
INVITRO_TARGET = (0.11, 0.24)       # sourced control baseline range, Hz

# MINIMUM CONTROL OUTPUT, derived from resp_metrics' own gates rather than chosen.
#
# resp_metrics declares an arm dead on EITHER an absolute floor (mean < 1.0, "mean output
# collapsed") or a relative one (mean < 0.2 * control). For the relative gate to be the
# binding one -- i.e. for there to be a GRADED range between full control and dead -- we
# need 0.2 * ctrl > 1.0, so ctrl mean > 5.0.
#
# Why it matters here: the best-scoring point of the first successful sweep had control mean
# 4.3, which puts "dead" at the absolute floor of 1.0 and leaves just 3.3 of range for a
# drug effect to live in. The LIF's 29.3 leaves 23.4. A deep, stationary, correctly-banded
# rhythm that no drug effect can be measured against is not a usable operating point, and
# optimising modulation alone selected exactly that.
#
# THE GATE IS THE DERIVED 5.0. The margin is a PREFERENCE, not a veto, and that
# distinction was corrected after it did damage.
#
# This was first set to 10.0 "for margin". At verification the best rhythm in the finalist
# set -- modulation 4.44, closest to the LIF's 4.736, and the tightest frequency SD of the
# four (0.015 Hz) -- came in at mean 9.5 and was DISQUALIFIED for being 0.5 below a number
# I had chosen, while satisfying the derived criterion comfortably (7.6 of graded range
# against the 1.0 absolute floor). The margin was overriding the thing it existed to
# protect, and the selection then fell to a point with modulation 2.43.
#
# So: hard gate at the derived 5.0; shortfall against 10 costs a little score instead of
# vetoing. Recorded in full because adjusting a threshold after seeing which candidate it
# excludes is exactly the move that needs a reason on the record -- and because the reason
# has to be that 5.0 is derived and 10.0 was not, rather than that I preferred the outcome.
MIN_CTRL_MEAN = 5.0
PREFERRED_CTRL_MEAN = 10.0

# The LIF weight table, as PreBotC assembles it. The search scales this as a block, so the
# RATIOS between pathways are held at the LIF's and only the overall scale moves. Holding the
# ratios is a choice, not a necessity: it keeps the search small and means the conductance
# network is the same circuit at a different scale rather than a differently-wired one, so a
# substrate difference cannot be a wiring difference in disguise.
LIF_W = dict(RESP_OP["w"])
LIF_W.update(ei_ampa=0.55, ei_nmda=0.0, ie_gaba=0.45, ie_gly=0.35,
             eo_ampa=0.70, eo_nmda=0.30)


# RECEPTOR TYPES NEED DIFFERENT SCALES, and this is the finding that made the network work.
#
# tau_nmda is 100 ms, 20x the AMPA tau, so with lumped population weights the NMDA
# conductance accumulates 20x more per unit firing rate. The LIF's weights were tuned on a
# substrate where the Mg2+ block held NMDA at a MEASURED mean relief of 0.063
# (scripts/diag_substrate_limits.py). A depolarised conductance cell reaches ~0.70, so the
# same weight delivers ~11x more NMDA -- and relief RISES with depolarisation, making it
# positive feedback into depolarisation block rather than a simple scale error.
#
# So:   AMPA, GABA, glycine  ->  x (g_L ratio) = 0.28
#       NMDA                 ->  x (g_L ratio) x (0.063/0.70) = 0.025
#
# Derived, then confirmed by prediction: of twelve points probed across both scales, the
# single one that produced a living rhythm was exactly (0.280, 0.025).
GL_RATIO = 2.8 / 10.0
RELIEF_RATIO = 0.063 / 0.70


def make_op(s_ampa, s_nmda, s_inh, drive, drive_other, gaba_tonic):
    w = {}
    for k in REQUIRED_W:
        if k.endswith("_nmda"):
            w[k] = LIF_W[k] * s_nmda
        elif k.endswith("_ampa"):
            w[k] = LIF_W[k] * s_ampa
        else:                                   # ie_gaba, ie_gly
            w[k] = LIF_W[k] * s_inh
    return dict(drive=float(drive), drive_other=float(drive_other),
                gaba_tonic=float(gaba_tonic), w=w)


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
    s_ampa, s_nmda, s_inh, drive, drive_other, gaba_tonic, _seconds, seed = args
    op = make_op(s_ampa, s_nmda, s_inh, drive, drive_other, gaba_tonic)
    total_ms = SETTLE_MS + 2.0 * WINDOW_MS
    base = dict(s_ampa=s_ampa, s_nmda=s_nmda, s_inh=s_inh, drive=drive,
                drive_other=drive_other, gaba_tonic=gaba_tonic)
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
            halves.append(resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m],
                                       band=INVITRO_BAND))
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
    elif mean < MIN_CTRL_MEAN:
        # Hard gate, not a penalty: a rhythm with too little output cannot host a graded
        # drug measurement however good it looks otherwise, and a penalty term would let a
        # beautifully modulated but unmeasurable point win.
        score, reason = float("inf"), f"output too low to measure against: {mean:.1f} < {MIN_CTRL_MEAN}"
    else:
        # Frequency must land in the sourced in vitro baseline range; outside it (but still
        # inside INVITRO_BAND) costs proportionally. Within it, what is optimised is
        # MODULATION DEPTH -- a rhythm that only just clears the 0.8 modulation gate is one
        # noise realisation from being reported dead, and leaves no range to measure a drug
        # against.
        lo, hi = INVITRO_TARGET
        off = 0.0 if lo <= freq <= hi else (lo - freq) / lo if freq < lo else (freq - hi) / hi
        short = max(0.0, (PREFERRED_CTRL_MEAN - mean) / PREFERRED_CTRL_MEAN)
        score = off + 2.0 / max(mod, 1e-6) + 0.5 * short
        reason = "alive, stable"
    return dict(**base, score=float(score), freq=freq, mod=mod, mean=mean,
                alive=alive, drift=float(drift), reason=reason)


def sweep(grid, seconds, procs, seed=0):
    jobs = [(sa, sn, si, d, do, gt, seconds, seed) for sa, sn, si, d, do, gt in grid]
    if procs > 1:
        with mp.Pool(procs) as pool:
            return pool.map(evaluate, jobs)
    return [evaluate(j) for j in jobs]


def report(rows, n=12):
    rows = sorted(rows, key=lambda r: r["score"])
    print(f"\n{'ampa':>6}{'nmda':>7}{'inh':>6}{'drive':>7}{'d_oth':>7}{'gtonic':>7} | "
          f"{'freq':>7}{'mod':>7}{'mean':>8}{'drift':>7} {'alive':>6}{'score':>8}  reason")
    print("-" * 112)
    for r in rows[:n]:
        sc = "inf" if not np.isfinite(r["score"]) else f"{r['score']:.4f}"
        print(f"{r['s_ampa']:6.3f}{r['s_nmda']:7.4f}{r['s_inh']:6.3f}{r['drive']:7.1f}"
              f"{r['drive_other']:7.1f}{r['gaba_tonic']:7.2f} | {r['freq']:7.3f}"
              f"{r['mod']:7.2f}{r['mean']:8.1f}{100*r['drift']:6.1f}% "
              f"{str(r['alive']):>6}{sc:>8}  {r['reason']}")
    alive = [r for r in rows if r["alive"]]
    print(f"\n{len(alive)}/{len(rows)} points produced a living rhythm")
    return rows[0] if rows and np.isfinite(rows[0]["score"]) else None


# Stage-2 finalists, carried forward verbatim from its output. Short sweeps are not
# trustworthy for frequency: the stage-1 winner read 1.026 Hz at 8 s and 1.102 Hz at 12 s,
# a 7% move from duration alone. So the winner of a sweep is a CANDIDATE, and stage 3 is
# what decides -- longer runs, several seeds, and a requirement that every seed be alive.
# Picking on a single short run would be choosing an operating point partly on seed noise.
FINALISTS = [
    # (s_ampa, s_nmda, s_inh, drive, drive_other, gaba_tonic)
    # From the stage-2 output (41/48 alive). Only points that passed every hard gate --
    # alive in the in vitro band, drift <= 20%, control mean >= 10 -- are carried forward.
    (0.40, 0.025, 0.21, 2.0,  5.0, 0.21),   # 0.270 Hz, mod 4.11, mean 10.3, drift 0.0%
    (0.28, 0.025, 0.21, 0.0, 12.0, 0.21),   # 0.270 Hz, mod 2.42, mean 15.5, drift 0.0%
    (0.28, 0.025, 0.21, 0.0,  5.0, 0.21),   # 0.270 Hz, mod 2.39, mean 13.3, drift 0.0%
    (0.28, 0.025, 0.21, 2.0,  5.0, 0.21),   # 0.270 Hz, mod 2.03, mean 16.2, drift 0.0%
]


def verify(a):
    """Re-measure the finalists at length, across seeds. Every seed must be alive."""
    print(f"stage 3: {len(FINALISTS)} finalists x {a.seeds} seeds, {a.seconds:.0f} s each")
    print(f"target: frequency in {INVITRO_TARGET[0]}-{INVITRO_TARGET[1]} Hz (band "
          f"{INVITRO_BAND[0]}-{INVITRO_BAND[1]}), control mean >= {MIN_CTRL_MEAN}, "
          f"maximising modulation.\n")
    print(f"{'ampa':>6}{'nmda':>7}{'inh':>6}{'drive':>7}{'d_oth':>7}{'gtonic':>7} | {'alive':>7}"
          f"{'freq':>9}{'sd':>7}{'mod':>7}{'mean':>8}{'|dfreq|':>9}")
    print("-" * 82)
    best = None
    for sa, sn, si, drive, do, gt in FINALISTS:
        jobs = [(sa, sn, si, drive, do, gt, a.seconds, s) for s in range(a.seeds)]
        if a.procs > 1:
            with mp.Pool(min(a.procs, len(jobs))) as pool:
                rows = pool.map(evaluate, jobs)
        else:
            rows = [evaluate(j) for j in jobs]
        n_alive = sum(r["alive"] for r in rows)
        fq = np.array([r["freq"] for r in rows], float)
        md = float(np.mean([r["mod"] for r in rows]))
        mn = float(np.mean([r["mean"] for r in rows]))
        # Scored against the SOURCED in vitro range, not the LIF's frequency: the two
        # substrates represent different preparations and must not share a frequency target.
        lo, hi = INVITRO_TARGET
        fmean = float(fq.mean())
        dfreq = (0.0 if lo <= fmean <= hi
                 else (lo - fmean) / lo if fmean < lo else (fmean - hi) / hi)
        print(f"{sa:6.3f}{sn:7.4f}{si:6.3f}{drive:7.1f}{do:7.1f}{gt:7.2f} | "
              f"{n_alive:4d}/{a.seeds}"
              f"{fq.mean():9.3f}{fq.std():7.3f}{md:7.2f}{mn:8.1f}{100*dfreq:8.1f}%")
        # EVERY SEED MUST BE ALIVE. An operating point that is alive in 3 of 4 seeds is
        # not a usable control: the drug arms would then be compared against a baseline
        # that sometimes does not exist, and a drug would get credit for a seed that was
        # already dead.
        # EVERY hard gate must hold in EVERY seed, then modulation decides. Ranking on
        # frequency proximity alone would reselect the deep-but-unmeasurable family.
        ok = (n_alive == a.seeds and mn >= MIN_CTRL_MEAN
              and all(abs(x - fmean) / max(fmean, 1e-9) <= 0.20 for x in fq))
        short = max(0.0, (PREFERRED_CTRL_MEAN - mn) / PREFERRED_CTRL_MEAN)
        rank = dfreq + 2.0 / max(md, 1e-6) + 0.5 * short
        if ok and (best is None or rank < best[0]):
            best = (rank, sa, sn, si, drive, do, gt, fmean, md, mn)

    if best is None:
        print("\nNO finalist passed every gate in every seed. Do not register any of them; "
              "widen the grid instead. Do NOT lower the seed requirement or MIN_CTRL_MEAN "
              "-- that value is derived from resp_metrics' own dead-arm criteria, not "
              "chosen, and is the one threshold here that must not move.")
        return
    rank, sa, sn, si, drive, do, gt, fq, md, mn = best
    print("\nREGISTER THIS:\n")
    print("COND_RESP_OP = MappingProxyType(dict(")
    print(f"    drive={drive}, drive_other={do}, gaba_tonic={gt},")
    wb = make_op(sa, sn, si, drive, do, gt)["w"]
    print("    w=MappingProxyType(dict("
          + ", ".join(f"{k}={v:.6g}" for k, v in wb.items()) + ")),")
    print("))")
    print(f"\n# s_ampa {sa}, s_nmda {sn}, s_inh {si}; control {fq:.3f} Hz "
          f"({60*fq:.1f} bursts/min, sourced in vitro range {60*INVITRO_TARGET[0]:.1f}-"
          f"{60*INVITRO_TARGET[1]:.1f}/min), mod {md:.2f} (LIF {TARGET['mod']}), "
          f"mean {mn:.1f}. The LIF's {TARGET['freq']} Hz is a DIFFERENT preparation and is "
          f"not a target.")


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
        # THE DEEP-MODULATION REGION. At s_inh 0.14 the network gave modulation 7.2-9.9
        # (better than the LIF's 4.74) with zero drift, at 0.30-0.34 Hz -- which read as
        # dead only because it was being gated against an in vivo band. Searching around
        # there now, with the correct band and quality-based scoring.
        grid = [(0.28, 0.025, si, d, 5.0, gt)
                for si, d, gt in itertools.product(
                    (0.10, 0.14, 0.17),             # inhibitory scale
                    (-2.0, 0.0, 2.0, 4.0),          # Exc drive, pA
                    (0.14, 0.21, 0.30))]            # gaba_tonic, nS
    else:
        # STAGE 2: find a point that is deep AND measurable. The first successful sweep gave
        # two families, and neither passed both requirements:
        #
        #   mod 7.0-9.7, mean 1.8-6.5  at 0.135 Hz  -> deep but unmeasurable
        #   mod 2.1-2.5, mean 18.8-22.5 at 0.270 Hz -> measurable but shallow
        #
        # `eo_ampa`/`eo_nmda` set the Out population's output and are scaled by s_ampa, so
        # s_ampa is the lever on mean. It is bounded above by depolarisation block (Out sat
        # at -19.9 mV firing 0 Hz at s_ampa 0.28 before the NMDA correction), so this
        # searches upward carefully and the drift and block gates do the rejecting.
        grid = [(sa, 0.025, si, d, do, gt)
                for sa, si, d, do, gt in itertools.product(
                    (0.28, 0.40, 0.55),             # AMPA scale -- the lever on mean output
                    (0.17, 0.21),                   # inhibitory scale
                    (0.0, 2.0),                     # Exc drive, pA
                    (5.0, 12.0),                    # Inh/Out bias, pA
                    (0.21, 0.30))]                  # gaba_tonic, nS

    print(f"stage {a.stage}: {len(grid)} operating points, {a.seconds:.0f} s each, "
          f"{a.procs} procs")
    print(f"target: frequency in the SOURCED in vitro baseline range "
          f"{INVITRO_TARGET[0]}-{INVITRO_TARGET[1]} Hz (validity band "
          f"{INVITRO_BAND[0]}-{INVITRO_BAND[1]}), maximising modulation depth.")
    print(f"the LIF control, for contrast only: {TARGET['freq']} Hz, mod {TARGET['mod']}, "
          f"mean {TARGET['mean']}  -- NOT a target; different preparation.")
    best = report(sweep(grid, a.seconds, a.procs))
    if best is None:
        print("\nNO operating point produced a living rhythm in band. Widen the grid; do "
              "NOT relax EUPNOEA_BAND to make a point qualify.")
        return
    print("\nBEST:")
    print("  COND_RESP_OP = dict(")
    print(f"      drive={best['drive']}, drive_other={best['drive_other']},")
    print(f"      gaba_tonic={best['gaba_tonic']},")
    wb = make_op(best["s_ampa"], best["s_nmda"], best["s_inh"],
                 best["drive"], best["drive_other"], best["gaba_tonic"])["w"]
    print("      w={" + ", ".join(f"{k!r}: {v:.6g}" for k, v in wb.items()) + "},")
    print("  )")
    print(f"  # s_ampa {best['s_ampa']}, s_nmda {best['s_nmda']}, s_inh {best['s_inh']}; "
          f"freq {best['freq']:.3f} Hz vs LIF {TARGET['freq']}, mod {best['mod']:.2f} "
          f"vs {TARGET['mod']}, drift {100*best['drift']:.1f}%")
    print("\nVerify at a longer duration and several seeds before registering it.")


# __main__ GUARD IS LOAD-BEARING. E8 bit this project twice: once because scripts did
# destructive work at import, and once because a heredoc-fed file made __main__ resolve to
# "<string>", so macOS spawn workers could not re-import it -- 13,334 tracebacks and 144 MB
# of output. multiprocessing here makes that failure mode live again.
if __name__ == "__main__":
    main()
