"""UNCERTAINTY PROPAGATION — can this model actually tell two compounds apart?

THE QUESTION THIS ANSWERS, which nothing else in the project does. Every number the
simulator reports is a POINT ESTIMATE computed from parameters that are estimates. The
module docstring of spinal/subtypes.py says it outright: the subunit fractions "are the
weakest numbers here... their ORDERING is solid; their precise values are not." If the
spread those estimates induce is wider than the gap between two compounds, then the model
CANNOT RANK THEM, and every comparison in the project so far is inside its own noise.

That is a question about honesty, not resolution, and no amount of added physical detail
answers it. It is also cheap to answer, so it should have been answered first.

WHAT IS SAMPLED (the parameters the project itself flags as weak):
  1. regional subunit fractions -- Dirichlet around the nominal values, preserving their
     mean and ordering but admitting real spread. These are literature estimates, not
     proteomics.
  2. the preBotC calibration anchor -- lognormal around 0.15. The underlying clinical data
     is midazolam 2 mg IV at -14.3 +/- 5.9% and -19 +/- 7% minute ventilation, and the
     local calibration table maps that spread to roughly sens in [0.09, 0.25].

WHAT IS HELD FIXED, and why: the efficacy ceiling (2.5) and the circuit parameters. Those
carry their own uncertainty, so every interval below is a LOWER BOUND on the true spread.

PAIRING. Both arms are evaluated under the SAME parameter draw and the SAME seeds, so the
per-draw DIFFERENCE cancels shared systematic error -- the same logic as the `pair` command.
Two intervals are therefore reported and they answer different questions:
  marginal  : how well do we know this compound's ventilation?        (wide)
  paired    : how well do we know which compound is better?           (narrower)
A model can be badly miscalibrated in absolute terms and still rank reliably. If even the
PAIRED interval spans zero, it cannot rank.

Run:  python scripts/uncertainty.py --draws 150
"""
import sys, os, argparse; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from spinal.resp import PreBotC, resp_metrics
from spinal.cpg import Drug
from spinal.subtypes import REGIONS, SUBTYPES, SUBJECTIVE_WEIGHT

RESP_OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
               w=dict(ee_ampa=0.45, ee_nmda=0.2475))
DUR_MS, DT, WARM = 14000.0, 0.1, 4000.0
CEILING = 2.5
SUBJ_TARGET = 0.50
KAPPA = 15.0          # Dirichlet concentration; lower = wider spread on the fractions
ANCHOR_MED, ANCHOR_SIG = 0.15, 0.31      # lognormal: 90% CI ~ [0.09, 0.25]

# the two arms being compared, as GABA-A subtype efficacies
ARMS = {
    "a5-selective (alogabat-like)": dict(a1=0.0, a23=0.10, a5=1.00, d_a4=0.0, eps=0.0),
    "non-selective BZ":             dict(a1=1.0, a23=1.00, a5=1.00, d_a4=0.0, eps=0.0),
}


def sample_params(rng):
    """One draw of the weak parameters. Returns (fractions_by_region, prebotc_anchor)."""
    frac = {}
    for region, nominal in REGIONS.items():
        alpha = np.array([max(1e-3, nominal[s] * KAPPA) for s in SUBTYPES])
        frac[region] = dict(zip(SUBTYPES, rng.dirichlet(alpha)))
    anchor = float(np.exp(rng.normal(np.log(ANCHOR_MED), ANCHOR_SIG)))
    return frac, anchor


def regional_sens(eff, frac, region, anchor):
    """sens = anchor * sum(f*e) / sum(f over BZ-sensitive subtypes).

    The denominator is what a NON-SELECTIVE BZ would produce in this region, so the anchor
    is reproduced by construction for that reference compound -- the same per-region
    calibration used in spinal/subtypes.py, recomputed here from SAMPLED fractions instead
    of mutating module globals (multiprocessing uses spawn; mutating globals in a worker is
    error E8 territory).
    """
    f = frac[region]
    raw = sum(f[s] * eff.get(s, 0.0) for s in SUBTYPES)
    ref = f["a1"] + f["a23"] + f["a5"]          # non-selective BZ reaches gamma2 subtypes
    return anchor * raw / max(1e-9, ref)


def subjective_index(eff, frac):
    f = frac["forebrain"]
    return sum(f[s] * eff.get(s, 0.0) * SUBJECTIVE_WEIGHT[s] for s in SUBTYPES)


def job(a):
    draw, arm, seed = a
    rng = np.random.default_rng(1000 + draw)        # same draw -> same params for all arms
    frac, anchor = sample_params(rng)
    eff = ARMS[arm]
    si = subjective_index(eff, frac)
    if si <= 1e-9:
        return draw, arm, np.nan, np.nan, False
    gain = min(1.0 + SUBJ_TARGET / si, CEILING)
    reached = (1.0 + SUBJ_TARGET / si) <= CEILING + 1e-9
    sens = regional_sens(eff, frac, "prebotc", anchor)
    d = Drug(gaba_a_gain=gain, gaba_a_tau=1.0 + 0.6 * (gain - 1.0),
             gaba_a_efficacy_cap=CEILING)
    b = PreBotC(drug=d, gaba_sens=sens, seed=seed, **RESP_OP)
    for i in range(int(DUR_MS / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > WARM
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return draw, arm, r["mean"], sens, reached


def ctrl_job(seed):
    b = PreBotC(drug=Drug(), gaba_sens=ANCHOR_MED, seed=seed, **RESP_OP)
    for i in range(int(DUR_MS / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > WARM
    return resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])["mean"]


def pct(x, q):
    return float(np.percentile(x, q))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=150)
    ap.add_argument("--seeds", type=int, default=2)
    a = ap.parse_args()

    jobs = [(d, arm, s) for d in range(a.draws) for arm in ARMS
            for s in range(a.seeds)]
    print(f"{len(jobs)} sims + {a.seeds} controls on {os.cpu_count()} cores "
          f"({a.draws} parameter draws x {len(ARMS)} arms x {a.seeds} seeds)")
    with Pool(os.cpu_count()) as p:
        ctrl = np.mean(p.map(ctrl_job, range(a.seeds)))
        out = p.map(job, jobs)

    V = {}      # (draw, arm) -> ventilation % of control
    S = {arm: [] for arm in ARMS}
    unreached = {arm: 0 for arm in ARMS}
    bad_draw = set()        # (draw) where ANY arm failed to reach target
    for draw, arm, mn, sens, reached in out:
        if np.isnan(mn):
            continue
        V.setdefault((draw, arm), []).append(100.0 * mn / ctrl)
        S[arm].append(sens)
        if not reached:
            unreached[arm] += 1
            bad_draw.add(draw)
    reached_all = set(range(a.draws)) - bad_draw
    per = {arm: np.array([np.mean(V[(d, arm)]) for d in range(a.draws)
                          if (d, arm) in V]) for arm in ARMS}

    print(f"\ndrug-free control ventilation: {ctrl:.2f} (model units)")
    print(f"subjective target {SUBJ_TARGET}, ceiling {CEILING}, "
          f"Dirichlet kappa {KAPPA}, anchor lognormal({ANCHOR_MED}, {ANCHOR_SIG})\n")

    # analytic knife-edge check -- no simulation needed, and it explains the reachability
    print("CEILING MARGIN at the NOMINAL parameters (analytic, no simulation):")
    for arm, eff in ARMS.items():
        si = subjective_index(eff, REGIONS)
        want = 1.0 + SUBJ_TARGET / si if si > 1e-9 else float("inf")
        marg = 100.0 * (CEILING - want) / CEILING
        flag = "REACHABLE" if want <= CEILING else "UNREACHABLE"
        print(f"  {arm:<32}needs gain {want:5.2f} vs ceiling {CEILING}  "
              f"margin {marg:+6.1f}%  {flag}")
    print("  A margin of a few percent means reachability is decided by the third")
    print("  significant figure of a literature-estimated subunit fraction.\n")

    print("MARGINAL — how well do we know each compound's ventilation?")
    print("  (includes under-dosed draws where the target was unreachable, so a good")
    print("   number here can mean 'safe' OR merely 'not actually delivering the effect')")
    print(f"  {'arm':<32}{'median':>9}{'5th':>8}{'95th':>8}{'width':>8}")
    print("  " + "-" * 63)
    for arm in ARMS:
        v = per[arm]
        print(f"  {arm:<32}{pct(v,50):8.0f}%{pct(v,5):7.0f}%{pct(v,95):7.0f}%"
              f"{pct(v,95)-pct(v,5):7.0f}pp")
        print(f"  {'  preBotC sens':<32}{np.median(S[arm]):8.3f}"
              f"{pct(S[arm],5):7.3f} {pct(S[arm],95):6.3f}")
        if unreached[arm]:
            print(f"  {'  NOTE':<32}subjective target unreachable in "
                  f"{100*unreached[arm]/(a.draws*a.seeds):.0f}% of draws")

    # REACHABILITY IS A FIRST-CLASS RESULT, not a footnote. A compound clamped to its
    # ceiling delivers LESS than the requested subjective effect, so its ventilation looks
    # good for the wrong reason -- it is being under-dosed. Comparing a full dose of one arm
    # against a subtherapeutic dose of another is not a comparison. The simulator already
    # guards this via Candidate._target_met; the first version of THIS script did not, and
    # the smoke test duly reported "better in 100% of draws" partly on that artefact.
    print("\nREACHABILITY — can the arm deliver the target subjective effect at all?")
    for arm in ARMS:
        n_un = unreached[arm] / max(1, a.seeds)
        print(f"  {arm:<32}target unreachable in {100*n_un/a.draws:5.0f}% of draws")

    arms = list(ARMS)
    common = [d for d in range(a.draws)
              if all((d, arm) in V for arm in arms) and d in reached_all]
    if not common:
        print("\n  NO DRAW lets both arms reach the target. The paired comparison is")
        print("  undefined: there is no matched-subjective-effect condition to compare in.")
        sys.exit(0)
    print(f"\n  matched-subjective draws usable for the paired contrast: "
          f"{len(common)}/{a.draws}")
    diff = np.array([np.mean(V[(d, arms[0])]) - np.mean(V[(d, arms[1])])
                     for d in common])
    print(f"\nPAIRED — how well do we know WHICH IS BETTER? ({len(common)} common draws)")
    print(f"  {arms[0]} minus {arms[1]}, per draw:")
    print(f"    median  {np.median(diff):+.1f} pp")
    print(f"    90% CI  [{pct(diff,5):+.1f}, {pct(diff,95):+.1f}] pp")
    frac_pos = float((diff > 0).mean())
    print(f"    the a5 arm is better in {100*frac_pos:.0f}% of draws")
    spans_zero = pct(diff, 5) < 0 < pct(diff, 95)
    print(f"\n  VERDICT: paired interval {'SPANS ZERO' if spans_zero else 'excludes zero'}"
          f" -> the model {'CANNOT' if spans_zero else 'CAN'} rank these two arms")
    print("           under its own parameter uncertainty.")
    print("\n  Every interval here is a LOWER BOUND: the ceiling and all circuit")
    print("  parameters were held fixed, and they carry uncertainty too.")
