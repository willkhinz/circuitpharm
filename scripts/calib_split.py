"""Calibrate the SPLIT (tonic/phasic) preBotC sensitivity against the midazolam anchor.

WHY THIS IS NECESSARY, AND IT IS ERROR E12 AGAIN. `scripts/overdose_kinetic.py` first ran
with the split sensitivities scaled so their TOTAL equalled 0.10 -- the value fitted for
the ONE-POOL drug. That is a calibration constant carried across a change of
parameterisation, which is exactly recurring error E12, and the worklog warned about this
specific migration one entry before I made it.

The damage: the split assigns the non-selective arm tonic 0.0177 / phasic 0.0823, so tonic
sensitivity fell 5.6x versus the 0.10 it was fitted at. Because the kinetic drug acts
almost ENTIRELY on the tonic pool (tonic gain ~7-11x vs phasic ~1.06x), the drug effect
collapsed: the non-selective arm came out at -2% ventilation where the anchor demands
-16 to -19%, and every arm then "survived" every dose. A reassuring result produced by an
uncalibrated parameter.

WHAT IS FITTED HERE. One scalar `k` mapping `GabaProfile.regional_sens_split()` output to
actual sensitivities, applied identically to every arm so each compound's SELECTIVITY
RATIO is preserved and only the overall scale is set by data.

ANCHOR. Human midazolam 2 mg IV, a clearly sedative dose and in the range that substitutes
for ethanol in drug discrimination: minute ventilation -14.3 +/- 5.9% (second cohort
-19 +/- 7%). The dose is defined as the one delivering the matched subjective target
(0.50) for the non-selective arm, which replaces the old "PAM 2.0x = sedative" fiat with
something tied to the discrimination literature.

Rate is NOT a target: compensatory tachypnoea is chemoreflex-mediated and this model has no
CO2 loop, so it is unreachable in principle (session 7c).

Run:  python scripts/calib_split.py
"""
import os
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.cpg import Drug
from circuitpharm.subtypes import PROFILES
from circuitpharm.gabaa_kinetics import fit_scheme, derive, calibrate_pam
from circuitpharm.config import RESP_OP, SYNAPTIC_PULSE, AMBIENT_GABA_UM

# Operating point, cleft pulse and ambient GABA come from circuitpharm.config --
# the single source. These were previously duplicated as literals in each script,
# which is HOW recurring error E12 (a calibration carried across a change of
# parameterisation) happened twice. dict() copies keep local mutation safe.
OP = dict(RESP_OP); OP["w"] = dict(RESP_OP["w"])
PULSE = dict(SYNAPTIC_PULSE)
AMBIENT = AMBIENT_GABA_UM
T, DT, WARM = 14000.0, 0.1, 4000.0
N_SEED = 3
SUBJ_TARGET = 0.50
S_MAX = 2.5                      # classical BZ intrinsic allosteric efficacy
TARGET_VENT = -16.5
VENT_RANGE = (-26.0, -8.0)       # clinical spread across both cohorts

KS = (0.3, 1.0, 3.0, 6.0, 10.0, 20.0, 40.0, 80.0)

SCHEME = fit_scheme(verbose=False, pulse=PULSE)
AFF = calibrate_pam(SCHEME, S_MAX, "affinity", pulse=PULSE, ambient_um=AMBIENT)
G = derive(SCHEME, affinity=AFF, ambient_um=AMBIENT, pulse=PULSE)


def ref_occupancy(key):
    """Occupancy delivering the matched subjective target. Linear in occ (mixture model)."""
    st, sp = PROFILES[key].subjective_index_split()
    per_occ = st * (G["tonic_gain"] - 1.0) + sp * (G["phasic_gain"] - 1.0)
    o = SUBJ_TARGET / per_occ if per_occ > 1e-12 else None
    return o if (o is not None and o <= 1.0) else None


def job(a):
    k, seed = a
    if k is None:
        b = PreBotC(drug=Drug(), gaba_sens=0.0, seed=seed, **OP)
    else:
        occ = ref_occupancy("nonselective_bz")
        lin = lambda g: 1.0 + occ * (g - 1.0)
        t, ph = PROFILES["nonselective_bz"].regional_sens_split("prebotc")
        b = PreBotC(drug=Drug(gaba_a_gain=lin(G["phasic_gain"]),
                              gaba_a_gain_tonic=lin(G["tonic_gain"]),
                              gaba_a_tau=lin(G["tau_ratio"]),
                              gaba_a_efficacy_cap=1e9, gaba_a_cap_tonic=1e9),
                    gaba_sens_tonic=min(1.0, t * k),
                    gaba_sens_phasic=min(1.0, ph * k), seed=seed, **OP)
    for i in range(int(T / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > WARM
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return k, r["mean"], r["amp"], float(r["alive"])


if __name__ == "__main__":
    occ = ref_occupancy("nonselective_bz")
    t, ph = PROFILES["nonselective_bz"].regional_sens_split("prebotc")
    print(f"kinetic drug at S_max {S_MAX}: tonic {G['tonic_gain']:.2f}x  "
          f"phasic {G['phasic_gain']:.3f}x  tau {G['tau_ratio']:.2f}x")
    print(f"reference occupancy for subj {SUBJ_TARGET}: {occ:.3f}")
    print(f"  -> delivered gains at that dose: tonic "
          f"{1+occ*(G['tonic_gain']-1):.2f}x  phasic {1+occ*(G['phasic_gain']-1):.3f}x")
    print(f"non-selective BZ split sens (unscaled): tonic {t:.4f}  phasic {ph:.4f}")
    print(f"\nfitting ONE scalar k against ventilation {TARGET_VENT}% "
          f"(clinical spread {VENT_RANGE[0]} to {VENT_RANGE[1]}%)\n")

    jobs = [(None, s) for s in range(N_SEED)] + [(k, s) for k in KS
                                                 for s in range(N_SEED)]
    with Pool(os.cpu_count()) as p:
        out = p.map(job, jobs)
    A = {}
    for k, mn, amp, al in out:
        A.setdefault(k, []).append((mn, amp, al))
    cv = np.mean([x[0] for x in A[None]])
    print(f"drug-free control ventilation {cv:.2f}")
    print(f"\n{'k':>7}{'tonic sens':>12}{'phasic sens':>13}{'ventilation':>13}"
          f"{'alive':>7}{'':>4}")
    print("-" * 60)
    best = None
    for k in KS:
        mn = np.mean([x[0] for x in A[k]])
        al = np.mean([x[2] for x in A[k]])
        # A quiescent control makes this 0/0. Reporting NaN says "no baseline to
        # measure against"; ZeroDivisionError kills the whole sweep at one bad cell.
        dv = 100 * (mn - cv) / cv if abs(cv) > 1e-9 else float("nan")
        ok = VENT_RANGE[0] <= dv <= VENT_RANGE[1]
        print(f"{k:7.1f}{min(1.0,t*k):12.4f}{min(1.0,ph*k):13.4f}"
              f"{dv:12.1f}%{al:7.2f}{'  OK' if ok else '':>4}")
        e = abs(dv - TARGET_VENT)
        # NaN-SAFE SELECTION. `err < best[0]` is False whenever best[0] is NaN, so a
        # single non-finite error in the FIRST grid cell latched `best` permanently and
        # every later finite, better candidate was silently discarded -- the script then
        # printed that poisoned cell as BEST. A grid search that reports the first point
        # it tried, with NaN% beside it, looks like a converged answer.
        #
        # Note this script is where I MADE the path reachable: the review-7 fix above
        # turned a ZeroDivisionError on a quiescent control into a NaN `dv`, which lands
        # straight here. A fix that converts a loud failure into a quiet wrong answer is
        # worse than no fix, and that is the second time in this project.
        if np.isfinite(e) and (best is None or e < best[0]):
            best = (e, k, dv)
    if best is None:
        raise SystemExit("no grid point produced a finite error -- every control was "
                         "quiescent or every run failed. Nothing to calibrate against.")
    print(f"\nBEST k = {best[1]:.1f}  -> ventilation {best[2]:+.1f}% "
          f"(target {TARGET_VENT}%)")
    print(f"  non-selective BZ sensitivities at that k: "
          f"tonic {min(1.0,t*best[1]):.4f}  phasic {min(1.0,ph*best[1]):.4f}")
    for key in ("alogabat", "mp_iii_022", "ideal_a5", "hz_166"):
        tt, pp = PROFILES[key].regional_sens_split("prebotc")
        print(f"  {key:<16} tonic {min(1.0,tt*best[1]):.4f}  "
              f"phasic {min(1.0,pp*best[1]):.4f}")
    print("\nUse this k in scripts/overdose_kinetic.py. The earlier run of that script")
    print("used a total-matched 0.10 instead and was therefore UNCALIBRATED -- it")
    print("reported -2% for the non-selective arm and 'everything survives'.")
