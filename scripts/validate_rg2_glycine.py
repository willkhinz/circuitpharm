"""Does the locomotor RG reproduce the strychnine phenotype? Tested by COMPLETE removal.

THE PUBLISHED VALIDATION TARGET. In lamprey fictive locomotion, strychnine (glycine
receptor block) ELIMINATES left-right alternation while robust rhythmic activity PERSISTS
-- same burst proportion, same rostro-caudal coordination, just co-activation instead of
alternation. That phenotype is the whole reason circuitpharm/rg2.py replaced the Matsuoka
oscillator: in Matsuoka, mutual inhibition IS the oscillator, so removing it stops the
rhythm, which contradicts the data.

So rg2 must show: glycine removed -> alternation lost, rhythm intact.

WHY THIS SCRIPT EXISTS SEPARATELY FROM scripts/tune_rg2.py. That tuning script tests the
phenotype with `Drug(glyr_gain=0.05)` -- a 95% block. That is recurring error E10: a 95%
conductance reduction previously left entrainment fully intact elsewhere in this project
(correlation -0.49), because entrainment is cheap. Partial block CANNOT detect whether a
coupling is load-bearing. The coupling must be removed COMPLETELY (ie_gly = 0).

Running the tuning sweep showed exactly the signature of that error: control and
"blocked" conditions came out identical in all 12 top rows (corr -0.29 vs -0.28; period
1242 vs 1242). That is uninformative, not negative.

ALSO FIXED HERE: error E8. This must be a FILE with a __main__ guard, not a heredoc piped
to stdin. macOS multiprocessing uses spawn, so workers re-import __main__; when __main__ is
"<string>" they cannot, each worker dies, the Pool respawns it forever, and the run emits
tracebacks until the disk fills (observed: 13,334 identical tracebacks, 144 MB).

Run:  python scripts/validate_rg2_glycine.py
"""
import os
import numpy as np
from multiprocessing import Pool
from circuitpharm.rg2 import GroupPacemakerRG
from circuitpharm.cpg import Drug, burst_metrics

# best-scoring cell from scripts/tune_rg2.py (drive 260, g_adapt 1.2, tau_adapt 280,
# ie_gly 3.0, ee_ampa 0.55) -- NOTE these are NOT the module defaults, which is itself a
# defect: the tuned values were only ever printed to stdout in a past session and never
# written back into rg2.py, so the module as shipped does not alternate at all.
TUNED = dict(drive=260.0, g_adapt=1.2, tau_adapt=280.0)
WEE = 0.55
T, DT, WARM = 30000.0, 0.1, 6000.0


def run(a):
    label, wgly, drug = a
    r = GroupPacemakerRG(drug=drug,
                         w=dict(ie_gly=wgly, ee_ampa=WEE, ee_nmda=WEE * 0.55),
                         seed=1, **TUNED)
    for i in range(int(T / DT)):
        r.step(DT)
        if i % 10 == 0:
            r.record()
    A = r.arrays(); m = A["t"] > WARM
    k = np.ones(50) / 50.0
    sm = lambda v: np.convolve(v, k, mode="same")
    f, e = sm(A["RG_F"][m]), sm(A["RG_E"][m])
    if f.std() < 1e-6 or e.std() < 1e-6:
        return label, None
    bf, be = burst_metrics(A["t"][m], f), burst_metrics(A["t"][m], e)
    return label, dict(corr=float(np.corrcoef(f, e)[0, 1]),
                       per_F=bf["period"], per_E=be["period"],
                       nb_F=bf["n_bursts"], nb_E=be["n_bursts"],
                       inrg=float(sm(A["InRG_F"][m]).mean()))


JOBS = [("coupled (ie_gly=3.0)",        3.0, Drug()),
        ("95% block (E10, wrong test)", 3.0, Drug(glyr_gain=0.05)),
        ("COMPLETE removal (ie_gly=0)", 0.0, Drug()),
        ("coupling x3 (ie_gly=9.0)",    9.0, Drug())]

if __name__ == "__main__":
    print(f"{len(JOBS)} sims on {min(os.cpu_count(), len(JOBS))} cores\n")
    with Pool(min(os.cpu_count(), len(JOBS))) as p:
        out = p.map(run, JOBS)

    print(f"{'condition':<30}{'corr':>8}{'per_F':>8}{'per_E':>8}"
          f"{'nb_F':>6}{'nb_E':>6}{'InRG Hz':>9}")
    print("-" * 75)
    res = {}
    for label, r in out:
        if r is None:
            print(f"{label:<30}{'ARRHYTHMIC':>8}")
            continue
        res[label] = r
        pf = r["per_F"] if r["per_F"] == r["per_F"] else float("nan")
        pe = r["per_E"] if r["per_E"] == r["per_E"] else float("nan")
        print(f"{label:<30}{r['corr']:+8.2f}{pf:8.0f}{pe:8.0f}"
              f"{r['nb_F']:6d}{r['nb_E']:6d}{r['inrg']:9.1f}")

    co = res.get("coupled (ie_gly=3.0)")
    rm = res.get("COMPLETE removal (ie_gly=0)")
    print("\n" + "=" * 75)
    print("VERDICT on the strychnine phenotype")
    print("=" * 75)
    if not co or not rm:
        print("  cannot judge: a required condition was arrhythmic")
    else:
        lost = rm["corr"] - co["corr"]        # positive => alternation weakened
        alive = rm["nb_F"] >= 8 and rm["nb_E"] >= 8
        print(f"  coupled correlation          {co['corr']:+.2f}")
        print(f"  correlation without glycine  {rm['corr']:+.2f}   "
              f"(change {lost:+.2f})")
        print(f"  rhythm survives removal      {'YES' if alive else 'NO'} "
              f"({rm['nb_F']} / {rm['nb_E']} bursts)")
        if lost > 0.25 and alive:
            print("\n  -> PHENOTYPE REPRODUCED: removing glycine weakened alternation while")
            print("     the rhythm persisted. The coupling is load-bearing for PHASE only,")
            print("     which is the published result and the reason for this architecture.")
        elif alive:
            print("\n  -> NOT REPRODUCED. The rhythm persists (correct) but alternation is")
            print("     UNCHANGED by removing glycine completely, so the glycinergic")
            print("     pathway is NOT what sets the phase in this model. Any apparent")
            print("     anti-phase correlation comes from the intrinsic frequency mismatch")
            print("     (`asym`), not from the coupling. The module must NOT be described")
            print("     as reproducing the strychnine phenotype until this is fixed.")
        else:
            print("\n  -> NOT REPRODUCED, and worse: removing glycine ABOLISHED the rhythm,")
            print("     which is the Matsuoka failure mode this architecture was built to")
            print("     avoid.")
