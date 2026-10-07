"""Calibrate preBotC GABA-A sensitivity against measured benzodiazepine respiratory data.

ANCHOR 1 (primary). Human midazolam 2 mg IV -- a clearly sedative/anxiolytic dose, in the
range that substitutes for ethanol in drug discrimination -- produces:
    minute ventilation  -14.3 +/- 5.9%  (second cohort -19 +/- 7%; younger -6 to -9%)
    tidal volume        -16 to -22%
    respiratory rate    +10% (COMPENSATORY INCREASE)
So PAM 2.0x in model units is defined as "clinically sedative benzodiazepine dose".

ANCHOR 2 (slope check). Propofol at a sedative dose reduces tidal volume ~60%. Propofol is
far more efficacious than a BZ-site PAM, so it sits past the BZ ceiling (modelled 5x).

Readout is decomposed into RATE and AMPLITUDE because the clinical signature is
amplitude-dominant with rate preserved or slightly raised -- a pattern the model must
reproduce, not just the scalar ventilation number.

gaba_sens = fraction of preBotC GABA-A conductance that is drug-modulatable. Biologically
motivated: BZ-site PAMs need gamma2, and the respiratory network also expresses delta,
alpha4 and epsilon subunits, with epsilon (enriched on NK1R+ rhythm neurons) BZ-insensitive.
"""
import sys, os; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.cpg import Drug

OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
          w=dict(ee_ampa=0.45, ee_nmda=0.2475))
N_SEED, T, DT, WARM = 4, 14000.0, 0.1, 4000.0
SENS = (0.40, 0.25, 0.15, 0.10, 0.06, 0.03)
PAMS = (1.0, 2.0, 5.0)

TARGET_VENT_PAM2 = -16.5      # %
TARGET_AMP_PAM5 = -60.0       # %
TARGET_RATE_PAM2 = +10.0      # %


def job(a):
    sens, pamx, seed = a
    d = Drug(gaba_a_gain=pamx, gaba_a_tau=1.0 + 0.6 * (pamx - 1.0))
    b = PreBotC(drug=d, gaba_sens=sens, seed=seed, **OP)
    for i in range(int(T / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > WARM
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    span_s = (A["t"][m][-1] - A["t"][m][0]) / 1000.0
    rate = r["n"] / span_s * 60.0                     # breaths per minute
    return sens, pamx, seed, r["mean"], r["amp"], rate, float(r["alive"])


if __name__ == "__main__":
    jobs = [(s, p, k) for s in SENS for p in PAMS for k in range(N_SEED)]
    with Pool(os.cpu_count()) as pool:
        out = pool.map(job, jobs)
    agg = {}
    for s, p, k, mn, amp, rate, al in out:
        agg.setdefault((s, p), []).append((mn, amp, rate, al))
    print(f"targets: PAM2 vent {TARGET_VENT_PAM2:+.0f}%, PAM2 rate "
          f"{TARGET_RATE_PAM2:+.0f}%, PAM5 amplitude {TARGET_AMP_PAM5:+.0f}%\n")
    print(f"{'sens':>5} | {'ctrl bpm':>8} | {'PAM2 vent':>9} {'PAM2 amp':>8} "
          f"{'PAM2 rate':>9} | {'PAM5 vent':>9} {'PAM5 amp':>8} | {'err':>5}")
    best = None
    for s in SENS:
        g = lambda p, i: np.mean([x[i] for x in agg[(s, p)]])
        c_v, c_a, c_r = g(1.0, 0), g(1.0, 1), g(1.0, 2)
        pct = lambda p, i, c: 100 * g(p, i) / max(1e-9, c) - 100
        v2, a2, r2 = pct(2.0, 0, c_v), pct(2.0, 1, c_a), pct(2.0, 2, c_r)
        v5, a5 = pct(5.0, 0, c_v), pct(5.0, 1, c_a)
        err = (abs(v2 - TARGET_VENT_PAM2) / 16.5
               + abs(a5 - TARGET_AMP_PAM5) / 60.0
               + 0.5 * abs(r2 - TARGET_RATE_PAM2) / 25.0)
        print(f"{s:5.2f} | {c_r:8.0f} | {v2:8.1f}% {a2:7.1f}% {r2:8.1f}% | "
              f"{v5:8.1f}% {a5:7.1f}% | {err:5.3f}")
        if best is None or err < best[0]:
            best = (err, s, v2, a2, r2, a5)
    print(f"\nBEST: gaba_sens={best[1]:.2f} -> PAM2 vent {best[2]:+.1f}% "
          f"(target {TARGET_VENT_PAM2:+.0f}), PAM2 amp {best[3]:+.1f}%, "
          f"PAM2 rate {best[4]:+.1f}%, PAM5 amp {best[5]:+.1f}% "
          f"(target {TARGET_AMP_PAM5:+.0f})")
