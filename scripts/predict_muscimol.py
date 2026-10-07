"""The model's muscimol prediction — a ZERO-FREE-PARAMETER test of the respiratory axis.

WHY THIS IS THE RIGHT EXPERIMENT TO PREDICT. Sessions 7c/7d established that the
respiratory axis cannot be calibrated against human whole-body ventilation: two unknowns
(the ligand's effective potency and `gaba_sens`) against effectively one constraint, with
the amplitude/rate split set in vivo by a chemoreflex the model lacks. Fitting a
benzodiazepine requires knowing simultaneously what the PAM does to the conductance
(tonic/phasic split, intrinsic allosteric efficacy, and the bell-shaped concentration-
response found in 7e) AND how sensitive the network is to a conductance change. Those are
entangled, which is why the fit was unidentifiable.

A DIRECT ORTHOSTERIC AGONIST breaks the entanglement completely, and better than first
realised:

  * it bypasses the PAM mechanism entirely -- no affinity/gating split, no bell-shaped
    curve, no intrinsic-efficacy unknown;
  * and `gaba_sens` is not merely better constrained, it is FIXED AT 1.0. The only reason
    `gaba_sens` is below 1 is that benzodiazepine-site ligands require a gamma2 subunit and
    are inactive at delta-, alpha4- and epsilon-containing receptors. GABA and muscimol
    bind the ORTHOSTERIC site, which every GABA-A receptor has. So a bath agonist reaches
    the whole population.

That leaves NO free parameter. The curve below is therefore a genuine prediction of the
circuit as built, not a fit waiting for data, and it is falsifiable by one experiment.

MODELLING CHOICE. Bath-applied muscimol produces a persistent, non-desensitising standing
conductance, so it is modelled as a multiplier on the TONIC extrasynaptic GABA-A
conductance (`gaba_tonic`, 1.5 nS baseline) with sens = 1.0. It is NOT applied as a phasic
PAM: a bath agonist does not prolong synaptic decay the way a benzodiazepine does.

THE EXPERIMENT THIS PREDICTS, stated so it can be run and the model falsified:
    muscimol (or GABA) concentration-response on inspiratory burst FREQUENCY and AMPLITUDE
    in a rhythmic preBotC slice or arterially perfused preparation, spanning partial to
    complete suppression (roughly 0.1-10 uM, given that 5 uM gives full-or-partial
    suppression), in animals OLDER THAN P12.

The P12 requirement is not pedantry. In the preBotC, NKCC1 falls precipitously at P12 and
KCC2 rises to dominance, the curves crossing near P11, so in younger animals GABA-A is
weakly inhibitory or frankly depolarising -- which is why diazepam INCREASES burst
frequency in newborn en-bloc preparations and why that observation says nothing about
mature pharmacology. Twice in this project a sign conclusion was drawn from a preparation
without first checking its age. Don't make it three times.

Run:  python scripts/predict_muscimol.py
"""
import os
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.cpg import Drug

OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
          w=dict(ee_ampa=0.45, ee_nmda=0.2475))
BASE_TONIC = 1.5          # nS, the model's drug-free standing GABA-A conductance
MULTS = (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0)
N_SEED, T, DT, WARM = 4, 14000.0, 0.1, 4000.0


def job(a):
    mult, seed = a
    # direct orthosteric agonist: a standing conductance on the WHOLE receptor population.
    # sens = 1.0 because GABA/muscimol bind the orthosteric site present on every GABA-A
    # receptor, including the BZ-insensitive delta/alpha4/epsilon populations.
    b = PreBotC(drug=Drug(), gaba_tonic=BASE_TONIC * mult,
                gaba_sens_tonic=1.0, gaba_sens_phasic=1.0,
                seed=seed, **OP)
    for i in range(int(T / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > WARM
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return mult, r["mean"], r["amp"], r["freq"], r["mod"], float(r["alive"])


if __name__ == "__main__":
    jobs = [(m, s) for m in MULTS for s in range(N_SEED)]
    print(f"{len(jobs)} sims on {os.cpu_count()} cores")
    print(f"direct GABA-A agonist modelled as a standing conductance, "
          f"baseline {BASE_TONIC} nS, sens fixed at 1.0 (NO free parameters)\n")
    with Pool(os.cpu_count()) as p:
        out = p.map(job, jobs)
    A = {}
    for mult, mn, amp, fq, mod, al in out:
        A.setdefault(mult, []).append((mn, amp, fq, mod, al))

    c = np.mean([x[0] for x in A[1.0]])
    ca = np.mean([x[1] for x in A[1.0]])
    cf = np.mean([x[2] for x in A[1.0]])
    print(f"drug-free: ventilation {c:.2f}, amplitude {ca:.1f}, frequency {cf:.2f} Hz\n")
    print("PREDICTED muscimol-equivalent concentration-response")
    print(f"{'g_tonic x':>10}{'g_tonic nS':>12}{'ventilation':>13}{'amplitude':>11}"
          f"{'frequency':>11}{'mod':>7}{'rhythm':>9}")
    print("-" * 74)
    partial, full = None, None
    for m in MULTS:
        mn = np.mean([x[0] for x in A[m]])
        amp = np.mean([x[1] for x in A[m]])
        fq = np.mean([x[2] for x in A[m]])
        mod = np.mean([x[3] for x in A[m]])
        al = np.mean([x[4] for x in A[m]])
        state = "intact" if al >= 0.99 else ("FAILED" if al < 0.5 else "partial")
        print(f"{m:10.1f}{BASE_TONIC*m:12.2f}{100*mn/c:12.0f}%{100*amp/ca:10.0f}%"
              f"{100*fq/cf:10.0f}%{mod:7.2f}{state:>9}")
        if partial is None and 100 * mn / c < 70:
            partial = m
        if full is None and al < 0.5:
            full = m

    print("\n" + "=" * 74)
    print("FALSIFIABLE PREDICTIONS")
    print("=" * 74)
    print(f"  partial suppression (ventilation <70% of control) at g_tonic "
          f"x{partial if partial else '>'}{'' if partial else MULTS[-1]}")
    print(f"  rhythm failure at g_tonic "
          f"{'x%g' % full if full else '> x%g (not reached)' % MULTS[-1]}")
    if full and partial:
        print(f"  -> the window between partial suppression and failure spans only "
              f"{full/partial:.1f}x in conductance")
        print("     A measured muscimol curve should show the same NARROWNESS. If the real")
        print("     curve is much shallower, this circuit is over-sensitive to tonic")
        print("     inhibition and `gaba_sens` has been absorbing that error.")
    print("\n  The SHAPE is the test, not the absolute conductance: the experiment gives")
    print("  micromolar muscimol, the model gives nS, and the two are linked by a receptor")
    print("  density this model does not contain. So compare NORMALISED curves -- the ratio")
    print("  of the concentration causing failure to the concentration causing 50%")
    print("  suppression. That ratio is parameter-free in both.")
    print("\n  If the model's curve is too steep, the likely culprit is the single lumped")
    print("  tonic conductance: real tonic inhibition is distributed over subtypes with")
    print("  different affinities, which would broaden the response.")
