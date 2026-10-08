#!/usr/bin/env python3
"""Does the ORDERING of compounds by simulated respiratory burden survive the substrate?

ROADMAP LINKS 4-5, THE DELIVERABLE -- and a corrected one. The design document originally
said this step would test whether the SELECTIVITY RANKING is substrate-independent. That was
wrong: `scripts/ranking_robustness.py` never imports a circuit module (its own docstring says
"This is analytic -- no circuit simulation"), so the ranking cannot depend on the neuron model
and the VALIDATED headline was never at risk. See section 6 of
knowledge/05-design-conductance-substrate.md.

What DOES run through the circuit is everything `evaluation.evaluate()` simulates: the
UNCALIBRATED respiratory and motor endpoints, where shape and ORDERING between compounds are
usable and absolute scale is not. Those go through current = g*(E - V), which is exactly what
`scripts/diag_substrate_limits.py` measured the LIF truncating -- mean GABA-A driving force
10.4 mV against up to 95 mV in a real burst, and NMDA held in near-permanent Mg2+ block. So
the ordering is genuinely at risk, and this script measures whether it moves.

MATCHED OPERATING POINTS, NOT MATCHED PARAMETERS. The same nS weight does different things on
the two cells, so each substrate is anchored independently to the same observable (control
burst frequency; scripts/anchor_cond_resp.py) and every compound is then reported as a
FRACTIONAL CHANGE FROM ITS OWN SUBSTRATE'S CONTROL. Comparing absolute outputs across
substrates would be meaningless -- the two controls differ several-fold in mean output by
construction.

SAME SEEDS ON BOTH ARMS, or the substrate gets credit for seed noise.

WHAT THIS CANNOT SHOW. The absolute respiratory numbers stay VOID on both substrates: the
respiratory axis has no valid quantitative anchor, which is a problem about the OBSERVABLE
(whole-body ventilation is wrong for an isolated preBotC) and is untouched by what the neurons
are made of. And the VALIDATED selectivity ranking is unaffected either way. This script's
reach is the UNCALIBRATED tier only.

Run:  python scripts/compare_substrates.py [--seeds 3] [--arms alogabat,hz_166,...]
"""
import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from circuitpharm.evaluation import Compound, _simulate_resp, clear_control_cache  # noqa: E402

DEFAULT_ARMS = ("nonselective_bz", "alogabat", "mp_iii_022", "hz_166",
                "ideal_a5", "neurosteroid")


def controls(substrate, seeds):
    """Drug-free controls for one substrate, computed ONCE.

    These were originally recomputed inside `burden()`, i.e. once per arm -- six times the
    necessary work on the arm list's default, and the conductance substrate is ~4x the LIF's
    cost per simulation. Controls do not depend on the arm, so this was pure waste, and it is
    what made the first runs look hung rather than slow.
    """
    return [_simulate_resp(Compound("control", occupancy=0.0), seed=s, substrate=substrate)
            for s in seeds]


def burden(arm, substrate, seeds, occupancy, ctrl):
    """Fractional reduction in mean inspiratory output vs this substrate's own control.

    Mean output rather than frequency: a GABA-A PAM in this circuit reduces output amplitude
    well before it moves the cycle period, so frequency is the less sensitive readout and
    would compress the ordering into ties.
    """
    drug = [_simulate_resp(Compound.from_profile(arm, occupancy=occupancy),
                           seed=s, substrate=substrate) for s in seeds]
    # Paired by seed, then averaged -- not a ratio of averages, which would let one noisy
    # control seed dominate.
    frac = [1.0 - d["mean"] / c["mean"] if c["mean"] > 1e-9 else np.nan
            for c, d in zip(ctrl, drug)]
    alive = sum(d["alive"] for d in drug)
    return dict(frac=float(np.nanmean(frac)), sd=float(np.nanstd(frac)),
                alive=alive, n=len(seeds),
                ctrl_mean=float(np.mean([c["mean"] for c in ctrl])),
                drug_mean=float(np.mean([d["mean"] for d in drug])),
                ctrl_freq=float(np.mean([c["freq"] for c in ctrl])),
                drug_freq=float(np.mean([d["freq"] for d in drug])))


def spearman(a, b):
    """Rank correlation without scipy (it is a dependency here, but this keeps the
    comparison's arithmetic visible rather than delegated)."""
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--occupancy", type=float, default=1.0)
    ap.add_argument("--arms", type=str, default=",".join(DEFAULT_ARMS))
    a = ap.parse_args()
    arms = tuple(x.strip() for x in a.arms.split(",") if x.strip())
    seeds = list(range(a.seeds))

    rows = {}
    for sub in ("lif", "cond"):
        clear_control_cache()          # the key includes the substrate, but be explicit
        print(f"[{sub}] controls...", flush=True)
        ctrl = controls(sub, seeds)
        rows[sub] = {}
        for arm in arms:
            rows[sub][arm] = burden(arm, sub, seeds, a.occupancy, ctrl)
            r = rows[sub][arm]
            # Printed as it goes, flushed. A twelve-minute run that prints only at the end
            # is indistinguishable from a hung one, and partial results are worth keeping.
            print(f"[{sub}] {arm:<18} frac={r['frac']:+.4f} alive={r['alive']}/{r['n']}",
                  flush=True)

    print(f"\nRespiratory burden at occupancy {a.occupancy}, {a.seeds} seeds, paired by seed.")
    print("Each arm is reported as a FRACTIONAL reduction in mean inspiratory output from")
    print("ITS OWN substrate's drug-free control. Absolute values are not comparable across")
    print("substrates and are shown only to make the two controls' difference visible.\n")
    print(f"{'arm':<20}{'LIF frac':>11}{'sd':>7}{'COND frac':>11}{'sd':>7}"
          f"{'d(frac)':>9}{'LIF alive':>11}{'COND alive':>12}")
    print("-" * 90)
    for arm in arms:
        L, C = rows["lif"][arm], rows["cond"][arm]
        print(f"{arm:<20}{L['frac']:11.4f}{L['sd']:7.3f}{C['frac']:11.4f}{C['sd']:7.3f}"
              f"{C['frac']-L['frac']:+9.4f}{L['alive']:8d}/{L['n']}{C['alive']:9d}/{C['n']}")

    cl = rows["lif"][arms[0]]["ctrl_mean"]; cc = rows["cond"][arms[0]]["ctrl_mean"]
    fl = rows["lif"][arms[0]]["ctrl_freq"]; fc = rows["cond"][arms[0]]["ctrl_freq"]
    print(f"\ncontrols: LIF mean {cl:.2f} @ {fl:.3f} Hz   |   COND mean {cc:.2f} @ {fc:.3f} Hz")
    print(f"          mean output differs {cl/max(cc,1e-9):.1f}x -- which is why only the")
    print(f"          FRACTIONAL change is compared.")

    vl = np.array([rows["lif"][x]["frac"] for x in arms], float)
    vc = np.array([rows["cond"][x]["frac"] for x in arms], float)
    ok = np.isfinite(vl) & np.isfinite(vc)
    print("\n" + "=" * 90)
    if ok.sum() < 3:
        print("INCONCLUSIVE: fewer than three arms produced a finite burden on both")
        print("substrates, so there is no ordering to compare.")
        return
    rho = spearman(vl[ok], vc[ok])
    order_l = [arms[i] for i in np.argsort(-vl[ok])]
    order_c = [arms[i] for i in np.argsort(-vc[ok])]
    print(f"ORDERING (most to least respiratory burden)")
    print(f"  LIF : {' > '.join(order_l)}")
    print(f"  COND: {' > '.join(order_c)}")
    print(f"\nSpearman rho = {rho:+.4f} over {int(ok.sum())} arms")
    if order_l == order_c:
        print("\nORDERING PRESERVED. The simulated respiratory ordering is robust to the")
        print("neuron model, which strengthens the UNCALIBRATED tier -- it does not touch")
        print("the VALIDATED selectivity ranking (substrate-independent by construction) or")
        print("the VOID absolute numbers (blocked on an observable, not on the neurons).")
    else:
        print("\nORDERING MOVED. The LIF's truncated phasic driving force was shaping the")
        print("simulated pharmacology, so every simulated endpoint in the package inherits")
        print("that. This does NOT overturn the selectivity ranking, which never runs a")
        print("neuron -- but it does mean the simulated respiratory ordering cannot be")
        print("quoted without naming its substrate.")
    print("=" * 90)


# __main__ guard: E8 bit this project twice.
if __name__ == "__main__":
    main()
