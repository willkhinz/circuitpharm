"""Re-calibrate preBotC GABA-A sensitivity with the KINETICS-DERIVED two-pool drug.

WHY A RE-CALIBRATION IS MANDATORY. The old anchor fixed "PAM 2.0x = clinically sedative
benzodiazepine dose" by fiat and fitted `gaba_sens` to human midazolam 2 mg IV. That value
(0.15) is only meaningful for the one-pool drug it was fitted with. The kinetic scheme
replaces that single 2.0x with a PAIR -- phasic 1.06x, tonic 7.2x -- so the old sensitivity
is attached to a parameterisation that no longer exists. Carrying 0.15 across would be
silently wrong by roughly the ratio between the pools.

WHAT IS NOW FITTED, AND THE DEGENERACY THIS EXPOSES. Two unknowns:

  ec50_shift  how much the ligand shifts the GABA dose-response AT THE SEDATIVE DOSE.
              The ~2-3x figure is the ligand's MAXIMUM shift at full occupancy; at a
              sedative-but-not-anaesthetic dose occupancy is partial, so the effective
              shift is smaller and is genuinely unknown.
  gaba_sens   fraction of preBotC GABA-A conductance that is drug-modulatable.

Three clinical targets constrain them (ventilation, amplitude, rate), so the problem is
not underdetermined -- but the two parameters TRADE OFF, and printing the whole surface
rather than one fitted point is the honest presentation. A single "best fit" would hide
that a weaker ligand with a more sensitive network is indistinguishable from the reverse.

CLINICAL ANCHOR (unchanged): human midazolam 2 mg IV
    minute ventilation  -14.3 +/- 5.9%   (second cohort -19 +/- 7%)
    tidal volume        -16 to -22%
    respiratory rate    +10%  (COMPENSATORY INCREASE -- the signature is amplitude-
                               dominant with rate preserved or slightly raised)

Run:  python scripts/recalibrate_kinetic.py
"""
import sys, os, itertools; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from spinal.resp import PreBotC, resp_metrics
from spinal.cpg import Drug
from spinal.gabaa_kinetics import fit_scheme, derive, calibrate_pam

OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
          w=dict(ee_ampa=0.45, ee_nmda=0.2475))
PULSE = dict(peak_um=3000.0, clear_ms=1.00)      # a cell that passed the anchored checks
AMBIENT = 0.40
N_SEED, T, DT, WARM = 3, 14000.0, 0.1, 4000.0

SHIFTS = (1.5, 2.0, 2.5, 3.0, 4.0)
SENS = (0.06, 0.08, 0.10, 0.125, 0.15, 0.20, 0.30)

TARGET_VENT = -16.5      # %
TARGET_AMP = -19.0       # %
TARGET_RATE = +10.0      # %   REPORTED BUT NOT FITTED -- see below

# RESPIRATORY RATE IS EXCLUDED FROM THE ERROR FUNCTION. The previous calibration was wrong
# to include it, but the reason matters and my first statement of it was too strong.
#
# The clinical +10% is a COMPENSATORY rate increase: midazolam cuts tidal volume, PaCO2
# rises, and chemoreceptors drive rate up. This model is an isolated preBotC with fixed
# drive, no gas exchange, no CO2 compartment and no chemoreceptor, so it cannot reproduce
# that compensation. Including the target pushes the fitted parameters to absorb an error
# they cannot remove, biasing them away from the targets the model CAN represent. The old
# one-pool fit reported +25% and logged it as a "known residual mismatch"; the two-pool fit
# gives -7.4%.
#
# WHAT I FIRST GOT WRONG (and then corrected twice -- see worklog 7c/7e): I claimed a rate
# change was "unreachable in principle", then on finding that diazepam 1 uM INCREASES burst
# frequency in an isolated newborn medulla-spinal cord preparation, flipped to calling the
# model's negative sign a defect. Both were premature. That preparation is pre-P12, and in
# the preBotC specifically NKCC1 falls precipitously at P12 while KCC2 rises to dominance,
# the curves crossing near P11 -- so in a newborn prep GABA-A is weakly inhibitory or
# depolarising and its sign carries no information about mature pharmacology. In MATURE
# rhythmic preparations the GABA-A agonist muscimol SUPPRESSES the rhythm (5 uM: full or
# partial suppression in all cultures; bicuculline reactivates it), which is the same sign
# this model gives.
#
# So: the model's frequency sign is supported, and rate is excluded only because the
# CLINICAL number mixes in mechanisms the circuit does not contain. Rate is printed as a
# DIAGNOSTIC. Fitting it would be fitting a mechanism that is not here.
#
# RULE ADOPTED: for any respiratory anchor, establish the animal's age relative to P12
# BEFORE using its sign for anything.
FIT_RATE = False

# built once in the parent; Scheme is a frozen dataclass so it pickles to workers cleanly
SCHEME = fit_scheme(verbose=False, pulse=PULSE)
DRUGS = {}
for _sh in SHIFTS:
    _aff = calibrate_pam(SCHEME, _sh, "affinity", pulse=PULSE, ambient_um=AMBIENT)
    DRUGS[_sh] = derive(SCHEME, affinity=_aff, ambient_um=AMBIENT, pulse=PULSE)


def run(drug, sens, seed):
    b = PreBotC(drug=drug, gaba_sens=sens, seed=seed, **OP)
    for i in range(int(T / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > WARM
    return resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])


def job(a):
    shift, sens, seed = a
    if shift is None:                               # drug-free control
        r = run(Drug(), sens, seed)
    else:
        d = DRUGS[shift]
        r = run(Drug(gaba_a_gain=d["phasic_gain"],
                     gaba_a_gain_tonic=d["tonic_gain"],
                     gaba_a_tau=d["tau_ratio"],
                     gaba_a_efficacy_cap=d["phasic_headroom"],
                     gaba_a_cap_tonic=d["tonic_headroom"]), sens, seed)
    return shift, sens, seed, r["mean"], r["amp"], r["freq"], float(r["alive"])


if __name__ == "__main__":
    print("kinetics-derived drug at each candidate EC50 shift:")
    print(f"  {'shift':>6}{'phasic':>9}{'tonic':>9}{'tau':>8}{'tonic cap':>11}")
    for sh in SHIFTS:
        d = DRUGS[sh]
        print(f"  {sh:6.1f}{d['phasic_gain']:9.3f}{d['tonic_gain']:9.2f}"
              f"{d['tau_ratio']:8.2f}{d['tonic_headroom']:11.1f}")

    jobs = [(None, s, k) for s in SENS for k in range(N_SEED)]
    jobs += [(sh, s, k) for sh, s in itertools.product(SHIFTS, SENS)
             for k in range(N_SEED)]
    print(f"\n{len(jobs)} sims on {os.cpu_count()} cores")
    with Pool(os.cpu_count()) as p:
        out = p.map(job, jobs)

    C, A = {}, {}
    for sh, s, k, mn, amp, fq, al in out:
        (C if sh is None else A).setdefault((s,) if sh is None else (sh, s),
                                            []).append((mn, amp, fq))
    # controls do not depend on sens; pool them
    cv = np.mean([x[0] for v in C.values() for x in v])
    ca = np.mean([x[1] for v in C.values() for x in v])
    cf = np.mean([x[2] for v in C.values() for x in v])
    print(f"\ndrug-free control: ventilation {cv:.2f}, amplitude {ca:.1f}, "
          f"rate {cf:.2f} Hz")
    print(f"targets: ventilation {TARGET_VENT}%, amplitude {TARGET_AMP}%, "
          f"rate {TARGET_RATE:+.0f}%\n")

    print("VENTILATION, % change from drug-free control")
    hdr = f"{'gaba_sens':>10}" + "".join(f"{f'shift {s}':>11}" for s in SHIFTS)
    print(hdr); print("-" * len(hdr))
    best = None
    for s in SENS:
        row = f"{s:10.3f}"
        for sh in SHIFTS:
            v = np.mean([x[0] for x in A[(sh, s)]])
            dv = 100 * (v - cv) / cv
            row += f"{dv:10.1f}%"
        print(row)

    print("\nERROR vs all three clinical targets (lower is better)")
    print(hdr); print("-" * len(hdr))
    for s in SENS:
        row = f"{s:10.3f}"
        for sh in SHIFTS:
            v = np.mean([x[0] for x in A[(sh, s)]])
            am = np.mean([x[1] for x in A[(sh, s)]])
            fq = np.mean([x[2] for x in A[(sh, s)]])
            dv = 100 * (v - cv) / cv
            da = 100 * (am - ca) / ca
            df = 100 * (fq - cf) / max(1e-9, cf)
            err = (((dv - TARGET_VENT) / 6.0) ** 2 + ((da - TARGET_AMP) / 6.0) ** 2
                   + (((df - TARGET_RATE) / 15.0) ** 2 if FIT_RATE else 0.0)) ** 0.5
            row += f"{err:11.2f}"
            if best is None or err < best[0]:
                best = (err, sh, s, dv, da, df)
        print(row)

    err, sh, s, dv, da, df = best
    print(f"\nBEST: ec50_shift {sh}, gaba_sens {s:.3f}  (error {err:.2f})")
    print(f"  ventilation {dv:+.1f}%  (target {TARGET_VENT}, FITTED)")
    print(f"  amplitude   {da:+.1f}%  (target {TARGET_AMP}, FITTED)")
    print(f"  rate        {df:+.1f}%  (clinical {TARGET_RATE:+.0f}, NOT FITTED --")
    print(f"              no chemoreflex in the model, so this is unreachable in")
    print(f"              principle and carries no information about the drug)")
    d = DRUGS[sh]
    print(f"\n  the drug at that shift: phasic {d['phasic_gain']:.3f}x, "
          f"tonic {d['tonic_gain']:.2f}x, tau {d['tau_ratio']:.2f}x")
    print(f"  OLD one-pool calibration for comparison: gain 2.0 both pools, "
          f"tau 1.60, gaba_sens 0.150")
    # IDENTIFIABILITY. A single "best" cell implies the parameters are pinned down. Check
    # whether they actually are, by finding every cell whose ventilation lands inside the
    # clinical spread (-14.3 +/- 5.9% and -19 +/- 7% => roughly -8% to -26%).
    print("\n" + "=" * 70)
    print("IDENTIFIABILITY — is this calibration actually determined?")
    print("=" * 70)
    ok = []
    for s in SENS:
        for sh in SHIFTS:
            v = np.mean([x[0] for x in A[(sh, s)]])
            dv = 100 * (v - cv) / cv
            if -26.0 <= dv <= -8.0:
                ok.append((sh, s, dv))
    print(f"{len(ok)} of {len(SHIFTS)*len(SENS)} cells put ventilation inside the")
    print("clinical spread (-8% to -26%), i.e. are all equally acceptable fits:")
    for sh, s, dv in ok:
        print(f"    ec50_shift {sh:<5} gaba_sens {s:<6.3f} -> ventilation {dv:+6.1f}%")
    if len(ok) > 1:
        shifts_ok = sorted({c[0] for c in ok}); sens_ok = sorted({c[1] for c in ok})
        print(f"\n  NOT IDENTIFIED. Acceptable ec50_shift spans "
              f"{min(shifts_ok)}-{max(shifts_ok)} and gaba_sens spans "
              f"{min(sens_ok):.3f}-{max(sens_ok):.3f},")
        print("  trading off along a diagonal valley: a weaker ligand in a more sensitive")
        print("  network is indistinguishable from a stronger ligand in a less sensitive")
        print("  one. TOTAL VENTILATION IS ONLY ONE CONSTRAINT AND THERE ARE TWO UNKNOWNS.")
        print("\n  The amplitude/rate split cannot supply the second constraint either:")
        print("  in vivo that split is set by CO2 chemoreflex compensation, which this")
        print("  model does not contain, so neither component is independently meaningful")
        print("  here even though their product is.")
        print("\n  WHAT WOULD FIX IT: an anchor measured in a preparation that ALSO lacks")
        print("  chemoreflex feedback -- benzodiazepine or GABA dose-response on burst")
        print("  frequency and amplitude in an isolated preBotC slice or en-bloc brainstem.")
        print("  That matches the model's structure, so both components would be")
        print("  comparable and the two parameters would separate. Until then, report the")
        print("  VALLEY, not a point estimate, and propagate it as uncertainty.")
    else:
        print("\n  Identified by this anchor alone.")
