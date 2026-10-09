#!/usr/bin/env python3
"""Is a drug's change in respiratory output MORE bursts or LONGER bursts?

THE OPEN QUESTION THIS CLOSES. `scripts/compare_substrates.py` found that the non-selective
benzodiazepine's effect on modelled respiratory output REVERSES SIGN between the two neuron
substrates: a 7.9% reduction on the LIF cell, a 37.7% increase on the Butera conductance cell.
Two independent measures agreed that the change is in how much of the TIME the network is
active -- duty cycle +65.5% (cond) against -6.7% (LIF), with peak burst height essentially
unchanged on both. What could not be said was whether a higher duty cycle meant more bursts or
longer ones, because the only burst detector available reported 5.63 Hz against the FFT's
0.302 Hz, a 19-fold disagreement (see circuitpharm.bursts for why).

`circuitpharm.bursts.detect_bursts` resolves that, and the decomposition is exact in logs.
Mean output over a window is, to the accuracy of a rectangular burst approximation,

    mean  ~  amplitude  x  duty  =  amplitude  x  duration  x  frequency

so, taking logs, the fractional change in mean output SEPARATES additively:

    dlog(mean)  =  dlog(amplitude)  +  dlog(duration)  +  dlog(frequency)
                   -------------      --------------     ---------------
                   stronger bursts    longer bursts      more bursts

Each term is measured independently, and the residual between dlog(mean) and the sum of the
three is reported. A large residual means the rectangular approximation has failed -- bursts
with non-trivial shape -- and the decomposition should not be quoted. It is printed rather
than hidden because that is the check that tells you whether to believe the rest.

EVERY RUN'S FFT CROSS-CHECK IS REPORTED. `detect_bursts` returns `consistent`, comparing its
own frequency against the FFT peak within the FFT's own resolution. Any run where that is
False is excluded from the aggregate and listed, because it means the two measures of the same
quantity disagree and at least one is wrong. The previous detector had no such check, which is
how a 19x error survived into a worklog.

Run:  python scripts/decompose_burst_change.py [--seeds 3] [--arms nonselective_bz,alogabat]
"""
import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from circuitpharm.bursts import detect_bursts                      # noqa: E402
from circuitpharm.config import EUPNOEA_BAND, INVITRO_BAND         # noqa: E402
from circuitpharm.evaluation import Compound, _SUBSTRATE, RESP_OP  # noqa: E402

DEFAULT_ARMS = ("nonselective_bz", "alogabat", "hz_166", "neurosteroid")

#: The validity band belongs to the PREPARATION, and it sets the detector's time constants.
#: Reusing the in vivo band on the in vitro substrate is recurring error E17 and would make
#: min_gap_ms 2.5x too small -- i.e. it would reintroduce exactly the ripple-counting failure
#: this script exists to avoid.
BANDS = {"lif": EUPNOEA_BAND, "cond": INVITRO_BAND}


def run(arm, substrate, seed):
    """One simulation, returning the burst decomposition plus the raw window mean."""
    from circuitpharm.resp import PreBotC
    # `Compound(arm, ...)` sets NAME ONLY and leaves every subtype efficacy at 0.0, i.e. it
    # builds a drug-free compound wearing a drug's name. `from_profile` is the constructor
    # that carries the efficacy vector AND the profile's own s_max ceiling. Writing the
    # former produced a decomposition in which every term was exactly 0.0000 -- caught only
    # because identical-to-control is implausible for a modulator at full occupancy.
    cand = (Compound("control", occupancy=0.0) if arm == "control"
            else Compound.from_profile(arm, occupancy=1.0))
    cfg = _SUBSTRATE[substrate]
    st, sp = cand.sens("prebotc")
    kw = dict(**RESP_OP) if substrate == "lif" else dict(substrate="cond")
    b = PreBotC(drug=cand.drug("prebotc"), gaba_sens_tonic=st, gaba_sens_phasic=sp,
                seed=seed, **kw)
    for i in range(int(cfg["duration_ms"] / 0.1)):
        b.step(0.1)
        if i % 10 == 0:
            b.record()
    A = b.arrays()
    m = A["t"] > cfg["warm_ms"]
    t, out = A["t"][m], A["Out"][m]
    d = detect_bursts(t, out, BANDS[substrate])
    d["mean"] = float(out.mean())

    # THE DECOMPOSITION MUST BE EXACT, OR ITS RESIDUAL EATS THE RESULT. A first version used
    # (p95 - p5) as "amplitude" and decomposed log(mean). Those do not compose: `mean`
    # includes the inter-burst floor, which no burst term can account for, so the residual
    # ran to 55% of the effect and the script's own guard refused to quote a decomposition
    # whose direction was in fact unambiguous.
    #
    # Over a window containing whole cycles, with `base` the inter-burst floor,
    #
    #     mean - base = (total area above base) / T = elevation x duration / period
    #                 = elevation x duration x frequency
    #
    # which is exact up to the floor being steady. So the quantity decomposed is
    # log(mean - base), and the floor is reported as its own term, because a modulator that
    # raises or lowers tonic drive changes mean output WITHOUT touching the bursts -- the
    # likeliest single confound in this measurement.
    th = d.get("thresholds")
    if th is None or not d["onsets"]:
        d.update(base=np.nan, elev=np.nan, mean_above=np.nan)
        return d
    base = float(th["p5"])
    inside = np.zeros(t.size, bool)
    for a_, b_ in zip(d["onsets"], d["offsets"]):
        inside |= (t >= a_) & (t <= b_)
    d["base"] = base
    d["elev"] = float(np.mean(out[inside] - base)) if inside.any() else np.nan
    d["mean_above"] = float(out.mean() - base)
    return d


def dlog(new, old):
    return (float(np.log(new / old)) if np.isfinite(new) and np.isfinite(old)
            and new > 0 and old > 0 else np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--arms", default=",".join(DEFAULT_ARMS))
    a = ap.parse_args()
    seeds = list(range(a.seeds))
    arms = [x.strip() for x in a.arms.split(",") if x.strip()]

    for substrate in ("lif", "cond"):
        band = BANDS[substrate]
        print("=" * 100)
        print(f"SUBSTRATE {substrate.upper()}   band {band[0]}-{band[1]} Hz   "
              f"settle {_SUBSTRATE[substrate]['warm_ms']/1000:.0f} s   "
              f"window {(_SUBSTRATE[substrate]['duration_ms']-_SUBSTRATE[substrate]['warm_ms'])/1000:.0f} s")
        print("=" * 100)

        ctrl = [run("control", substrate, s) for s in seeds]
        bad = [(s, c["reason"]) for s, c in zip(seeds, ctrl) if not c["consistent"]]
        for s, why in bad:
            print(f"  CONTROL seed {s} EXCLUDED: {why}")
        ok_seeds = [s for s, c in zip(seeds, ctrl) if c["consistent"]]
        if not ok_seeds:
            print("  no usable control on this substrate; nothing to decompose\n")
            continue
        cmap = {s: c for s, c in zip(seeds, ctrl) if c["consistent"]}
        c0 = cmap[ok_seeds[0]]
        print(f"  control: {np.mean([cmap[s]['freq_hz'] for s in ok_seeds]):.4f} Hz "
              f"(FFT {np.mean([cmap[s]['fft_hz'] for s in ok_seeds]):.4f}, "
              f"bin {c0['fft_bin_hz']:.4f}), "
              f"duty {np.mean([cmap[s]['duty'] for s in ok_seeds]):.3f}, "
              f"{np.mean([cmap[s]['n_bursts'] for s in ok_seeds]):.1f} bursts, "
              f"mean {np.mean([cmap[s]['mean'] for s in ok_seeds]):.2f}")
        print(f"  derived: min_gap {c0['derived']['min_gap_ms']:.0f} ms, "
              f"min_dur {c0['derived']['min_dur_ms']:.0f} ms, "
              f"smooth {c0['derived']['smooth_ms']:.0f} ms")
        print()
        print(f"  {'arm':17s} {'dlog(m-b)':>10s} = {'elev':>8s} + {'duration':>9s} + "
              f"{'freq':>8s} | {'resid':>7s}  verdict")
        print("  " + "-" * 94)
        for arm in arms:
            rows = []
            for s in ok_seeds:
                d = run(arm, substrate, s)
                if not d["consistent"]:
                    print(f"  {arm:17s} seed {s} EXCLUDED: {d['reason'][:66]}")
                    continue
                c = cmap[s]
                rows.append((dlog(d["mean_above"], c["mean_above"]),
                             dlog(d["elev"], c["elev"]),
                             dlog(np.mean(d["durations"]), np.mean(c["durations"])),
                             dlog(d["freq_hz"], c["freq_hz"]),
                             dlog(d["mean"], c["mean"]),
                             dlog(d["base"], c["base"])))
            if not rows:
                print(f"  {arm:17s} no usable seed")
                continue
            M, A_, D, F, MT, B = (float(np.nanmean([r[i] for r in rows]))
                                  for i in range(6))
            resid = M - (A_ + D + F)
            # A modulator at full occupancy producing EXACTLY no change is a construction
            # error, not a result -- it is what a drug-free `Compound` looks like. Say so,
            # rather than letting the residual guard swallow it as "0 is not < 0.35*0".
            if abs(M) < 1e-12 and abs(A_) < 1e-12 and abs(D) < 1e-12 and abs(F) < 1e-12:
                print(f"  {arm:17s} IDENTICAL TO CONTROL in every term -- the compound was "
                      f"probably built without its efficacy vector")
                continue
            # The verdict names the dominant term WITH ITS SIGN. A first version took the
            # largest term by absolute value and labelled it "longer bursts" even when the
            # term was negative, so a shortening read as a lengthening.
            terms = {"bursts": (D, "longer", "shorter"),
                     "rate": (F, "more frequent", "less frequent"),
                     "strength": (A_, "stronger", "weaker")}
            lead = max(terms, key=lambda k: abs(terms[k][0]))
            val, up, down = terms[lead]
            second = sorted((abs(v[0]) for v in terms.values()), reverse=True)[1]
            dominates = abs(val) > 2.0 * second
            closes = abs(resid) < 0.35 * abs(M) if abs(M) > 1e-9 else False
            verdict = ("DO NOT QUOTE: residual too large" if not closes else
                       f"{up if val > 0 else down} bursts" if dominates else "mixed")
            print(f"  {arm:17s} {M:10.4f} = {A_:8.4f} + {D:9.4f} + {F:8.4f} | "
                  f"{resid:7.4f}  {verdict}")
            print(f"  {'':17s} {'(raw mean':>10s} {MT:8.4f}{', floor':>9s} {B:8.4f})")
        print()

    print("=" * 100)
    print("Reading the table: dlog(mean) is the natural log of the fractional change in mean")
    print("inspiratory output. Positive = more output. The three terms are independent")
    print("measurements and must sum to it; `resid` is how badly the rectangular-burst")
    print("approximation fails, and a verdict is withheld when it exceeds 35% of the total.")


if __name__ == "__main__":
    main()
