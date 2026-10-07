"""Overdose protection as a LIGAND property: does the rhythm survive full occupancy?

!! DO NOT QUOTE THIS SCRIPT'S NUMBERS (as of session 7d). The FORM is correct -- the
!! occupancy-mixture formulation below is right and worth keeping -- but the output is VOID
!! because the split tonic/phasic sensitivity it needs has no defensible calibration:
!!
!!   * the first run used SENS_TOTAL = 0.10, carried over from the ONE-POOL drug. That is
!!     recurring error E12. It collapsed the drug effect to -2% ventilation for a
!!     non-selective benzodiazepine (anchor: -16 to -19%) and consequently reported that
!!     every arm survives every dose up to full occupancy. That result was an artefact.
!!   * re-calibrating (scripts/calib_split.py) requires k=10, which CLIPS phasic
!!     sensitivity at 1.0 -- i.e. 100% of preBotC phasic GABA-A modulatable, contradicting
!!     the subunit argument that justifies gaba_sens being below 1 in the first place.
!!   * no affinity/gating mechanism mix reproduces the one-pool 2.0x on both pools.
!!
!! Conclusion recorded in the worklog: `gaba_sens` is not the subunit fraction its
!! docstring claims, but a lumped factor absorbing extra-preBotC mechanisms (chemoreflex
!! blunting, upper-airway tone) that this isolated model does not contain. The respiratory
!! axis therefore cannot be calibrated against human whole-body ventilation at all; it
!! needs isolated preBotC slice / en-bloc brainstem dose-response data.
!!
!! Until that anchor exists, this script is a correct instrument with no calibrated scale.

WHAT WAS WRONG BEFORE. `simulator.py` escalates dose by scaling the PAM gain linearly and
clipping at a hard-coded `gaba_a_efficacy_cap`. That cap is an invented number (2.5), it is
applied to both receptor pools although their headroom differs ~200-fold, and linear
scaling with a clip has no mechanistic saturation -- the "ceiling" is a fiat boundary, so
the resulting overdose index measures the boundary, not the pharmacology.

THE CORRECT MECHANISM. A positive allosteric modulator cannot exceed its own INTRINSIC
ALLOSTERIC EFFICACY: the maximum shift it can impose on the GABA dose-response when every
modulator site is occupied. Call it `S_max`. Raising the dose raises site occupancy, which
asymptotes at 1. So:

    occupancy(D) = D*c / (1 + D*c)          Langmuir; c fixes occupancy at the reference dose
    population   = occ * (modulated receptors) + (1 - occ) * (unmodulated receptors)

Because the population is a MIXTURE of shifted and unshifted receptors, each pool's
conductance gain is EXACTLY LINEAR in occupancy:

    gain(occ) = 1 + occ * (gain_at_full_occupancy - 1)

and therefore saturates at `gain_max`, which is set by `S_max` through the kinetic scheme.
That is the ceiling -- derived from the ligand, not asserted. (Linearity is exact for the
conductance gains. For the decay tau it is first-order only: a mixture of two exponentials
is not one exponential. Flagged, not hidden.)

THE TEST. Dose each arm to a MATCHED subjective effect (0.50), call that D = 1, then
escalate. Two possible outcomes, and they mean completely different things:

    survives at full occupancy  -> genuinely ceiling-protected. No dose kills it, because
                                   the drug runs out of mechanism before the network fails.
                                   Overdose index is INFINITE, not merely large.
    dies below full occupancy   -> a finite lethal multiple exists and the ceiling argument
                                   does NOT protect this compound.

This is the session's central question made testable: the project assumed overdose
protection is a property of the PAM CLASS. If it depends on `S_max`, it is a property of
the individual MOLECULE and has to be measured per candidate.

Run:  python scripts/overdose_kinetic.py
"""
import sys, os, itertools; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.cpg import Drug
from circuitpharm.subtypes import PROFILES
from circuitpharm.gabaa_kinetics import fit_scheme, derive, calibrate_pam

OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
          w=dict(ee_ampa=0.45, ee_nmda=0.2475))
PULSE = dict(peak_um=3000.0, clear_ms=1.00)
AMBIENT = 0.40
T, DT, WARM = 14000.0, 0.1, 4000.0
N_SEED = 2

# gaba_sens from the Session 7c re-calibration valley. It is NOT identified (0.06-0.30),
# so the mid value is used here and the extremes are reported as a sensitivity band.
SENS_TOTAL = 0.10
SUBJ_TARGET = 0.50
ARMS = ("alogabat", "nonselective_bz")
S_MAX = (2.5, 4.0, 6.0)          # ligand intrinsic allosteric efficacy (max EC50 shift)
DOSES = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, np.inf)   # inf = full occupancy

SCHEME = fit_scheme(verbose=False, pulse=PULSE)
MAXGAIN = {}
for _s in S_MAX:
    _aff = calibrate_pam(SCHEME, _s, "affinity", pulse=PULSE, ambient_um=AMBIENT)
    MAXGAIN[_s] = derive(SCHEME, affinity=_aff, ambient_um=AMBIENT, pulse=PULSE)


def arm_sens(key):
    """preBotC tonic/phasic sensitivities, partitioned by subtype localisation."""
    p = PROFILES[key]
    t, ph = p.regional_sens_split("prebotc")
    tot = t + ph
    # rescale so the TOTAL matches the calibrated SENS_TOTAL for a non-selective reference,
    # preserving each arm's selectivity ratio. Without this the arms are not on the same
    # calibration footing as the midazolam anchor.
    ref = sum(PROFILES["nonselective_bz"].regional_sens_split("prebotc"))
    k = SENS_TOTAL / ref
    return t * k, ph * k, tot * k


def occ_ref(key, s_max):
    """Occupancy that delivers the matched subjective target for this arm.

    delivered(occ) = occ * [subj_tonic*(gT-1) + subj_phasic*(gP-1)]  -- linear in occ.
    """
    g = MAXGAIN[s_max]
    st, sp = PROFILES[key].subjective_index_split()
    per_occ = st * (g["tonic_gain"] - 1.0) + sp * (g["phasic_gain"] - 1.0)
    if per_occ <= 1e-12:
        return None
    o = SUBJ_TARGET / per_occ
    return o if o <= 1.0 else None       # cannot reach the target even at full occupancy


def occupancy(D, o_ref):
    if not np.isfinite(D):
        return 1.0
    c = o_ref / max(1e-9, 1.0 - o_ref)
    return float(D * c / (1.0 + D * c))


def job(a):
    key, s_max, D, seed = a
    g = MAXGAIN[s_max]
    o_ref = occ_ref(key, s_max)
    if o_ref is None:
        return key, s_max, D, None, None
    occ = occupancy(D, o_ref)
    lin = lambda gmax: 1.0 + occ * (gmax - 1.0)
    st, sp, _ = arm_sens(key)
    d = Drug(gaba_a_gain=lin(g["phasic_gain"]),
             gaba_a_gain_tonic=lin(g["tonic_gain"]),
             gaba_a_tau=lin(g["tau_ratio"]),
             gaba_a_efficacy_cap=1e9,        # the MECHANISM caps it now, not a fiat number
             gaba_a_cap_tonic=1e9)
    b = PreBotC(drug=d, gaba_sens_tonic=st, gaba_sens_phasic=sp, seed=seed, **OP)
    for i in range(int(T / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > WARM
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return key, s_max, D, r["mean"], float(r["alive"])


def ctrl(seed):
    b = PreBotC(drug=Drug(), gaba_sens=SENS_TOTAL, seed=seed, **OP)
    for i in range(int(T / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > WARM
    return resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])["mean"]


if __name__ == "__main__":
    print("ligand intrinsic efficacy -> pool gains AT FULL OCCUPANCY (the real ceiling)")
    print(f"  {'S_max':>6}{'tonic gain':>12}{'phasic gain':>13}{'tau ratio':>11}")
    for s in S_MAX:
        g = MAXGAIN[s]
        print(f"  {s:6.1f}{g['tonic_gain']:12.2f}{g['phasic_gain']:13.3f}"
              f"{g['tau_ratio']:11.2f}")

    print(f"\nmatched-subjective reference occupancy (delivers subj {SUBJ_TARGET}):")
    for key, s in itertools.product(ARMS, S_MAX):
        o = occ_ref(key, s)
        txt = f"{o:.3f}" if o else "UNREACHABLE even at full occupancy"
        print(f"  {key:<20} S_max {s:<5} -> occupancy {txt}")
    print("\npreBotC sensitivities, partitioned by subtype localisation:")
    for key in ARMS:
        st, sp, tot = arm_sens(key)
        print(f"  {key:<20} tonic {st:.5f}  phasic {sp:.5f}  total {tot:.5f}")

    jobs = [(k, s, D, sd) for k, s, D in itertools.product(ARMS, S_MAX, DOSES)
            for sd in range(N_SEED)]
    print(f"\n{len(jobs)} sims + {N_SEED} controls on {os.cpu_count()} cores")
    with Pool(os.cpu_count()) as p:
        cv = np.mean(p.map(ctrl, range(N_SEED)))
        out = p.map(job, jobs)

    R = {}
    for key, s, D, mn, al in out:
        if mn is None:
            continue
        R.setdefault((key, s, D), []).append((mn, al))
    print(f"\ndrug-free control ventilation {cv:.2f}\n")

    print("VENTILATION % of control, by dose multiple of the matched-subjective dose")
    print("('X' = respiratory rhythm failed; 'inf' column = FULL modulator occupancy)")
    hdr = (f"{'arm':<18}{'S_max':>6}" +
           "".join(f"{('inf' if not np.isfinite(D) else f'{D:g}x'):>9}" for D in DOSES))
    print(hdr); print("-" * len(hdr))
    verdict = {}
    for key, s in itertools.product(ARMS, S_MAX):
        if (key, s, DOSES[0]) not in R:
            print(f"{key:<18}{s:6.1f}   subjective target unreachable")
            continue
        row = f"{key:<18}{s:6.1f}"
        died_at = None
        for D in DOSES:
            mn = np.mean([x[0] for x in R[(key, s, D)]])
            al = np.mean([x[1] for x in R[(key, s, D)]])
            cell = f"{100*mn/cv:.0f}%" if al >= 0.5 else "X"
            if al < 0.5 and died_at is None:
                died_at = D
            row += f"{cell:>9}"
        print(row)
        verdict[(key, s)] = died_at

    print(f"\n{'='*78}\nVERDICT — is the ceiling actually protective?\n{'='*78}")
    for (key, s), died in sorted(verdict.items()):
        if died is None:
            print(f"  {key:<20} S_max {s:<5} SURVIVES full occupancy -> "
                  f"overdose index INFINITE")
            print(f"  {'':<20} the drug exhausts its mechanism before the network fails.")
        elif not np.isfinite(died):
            print(f"  {key:<20} S_max {s:<5} fails ONLY at full occupancy -> "
                  f"protected at any finite dose")
        else:
            print(f"  {key:<20} S_max {s:<5} FAILS at {died:g}x the "
                  f"matched-subjective dose")
            print(f"  {'':<20} the ceiling does NOT protect this compound.")
    print("\nIf the verdict depends on S_max, then overdose protection is a property of")
    print("the INDIVIDUAL MOLECULE (its intrinsic allosteric efficacy, a measurable")
    print("quantity) and NOT of the PAM class. A novel candidate cannot inherit it.")
    print(f"\nCAVEAT: gaba_sens is fixed at {SENS_TOTAL} from a valley that spans")
    print("0.06-0.30 and is NOT identified (session 7c). These doses carry that")
    print("uncertainty; the QUALITATIVE verdict per S_max is the robust part.")
