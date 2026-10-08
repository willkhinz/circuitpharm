#!/usr/bin/env python3
"""Worked example: evaluate a compound, and see what the tool will and will not tell you.

Run:  python examples/quickstart.py
"""
from circuitpharm import calibration_report, Tier
from circuitpharm.evaluation import Compound, evaluate
from circuitpharm.results import VoidQuantityError


def main():
    # ------------------------------------------------------------------ 1. the ground rules
    print(__doc__)
    print("=" * 78)
    print("STEP 1 — ask the tool what it actually knows before trusting any number")
    print("=" * 78)
    print(calibration_report())

    # ------------------------------------------------- 2. describe a compound by MEASURED data
    print("\n" + "=" * 78)
    print("STEP 2 — describe a compound by its MEASURED receptor activity")
    print("=" * 78)
    print("""Efficacy is an INPUT, not a prediction. For an allosteric modulator, efficacy
is not affinity: it is the shift in the conformational equilibrium between shut, open and
desensitised states. Free-energy methods reproduce affinities and Markov state models have
captured a putative GABA-A open state, but there is no predictive efficacy calculation --
the field still publishes experimental screening designs for exactly this quantity, and
measured efficacy is probe-dependent. So you measure it and pass it in.

`s_max` is the ligand's intrinsic allosteric efficacy (max GABA EC50 shift at full
occupancy). It is what actually determines overdose protection, it is a property of the
individual molecule, and it also must be measured -- it cannot be inherited from the PAM
class.""")

    cand = Compound(
        name="hypothetical a5-selective PAM",
        a1=0.0,      # no a1 -> no sedation/respiratory arm
        a23=0.10,
        a5=1.00,     # a5 carries ethanol's discriminative stimulus
        s_max=2.5,       # classical BZ-site range
        occupancy=0.35,  # fraction of modulator sites bound -- NOT a mass dose
    )
    g = cand.pool_gains()
    print(f"\nderived from the Markov gating scheme (one calibrated parameter):")
    print(f"  tonic conductance gain   {g['tonic']:.2f}x")
    print(f"  phasic conductance gain  {g['phasic']:.3f}x")
    print(f"  decay tau ratio          {g['tau']:.2f}x")
    print(f"""
The pools differ ~{g['tonic']/g['phasic']:.0f}-fold and that is the central mechanistic result: an
affinity-type modulator barely moves a near-saturated synaptic peak but strongly
potentiates sub-saturating tonic current. Since a5 is predominantly EXTRASYNAPTIC, an
a5-selective drug acts almost entirely on the tonic pool.""")

    # ----------------------------------------------------------------- 3. evaluate
    print("\n" + "=" * 78)
    print("STEP 3 — evaluate, and read the tiers")
    print("=" * 78)
    rs = evaluate(cand, n_seed=2)
    print(rs)

    # ----------------------------------------------------------------- 4. the refusal
    print("\n" + "=" * 78)
    print("STEP 4 — what happens when you ask for something unjustifiable")
    print("=" * 78)
    try:
        rs.quantity("overdose_index").value
    except VoidQuantityError as e:
        print("Requesting the overdose index raises:\n")
        print("    " + "\n    ".join(str(e).split("\n")[:4]))
    print("""
That is the design. A caveat printed next to a number gets dropped when the number is
copied; this project's worst reporting error was exactly that. A tier that refuses does
not get dropped.""")

    # ----------------------------------------------------------------- 5. what to use
    print("\n" + "=" * 78)
    print("STEP 5 — the one quantity that is calibration-independent")
    print("=" * 78)
    for ref in ("nonselective_bz",):
        print(f"  selectivity ratio vs {ref}: "
              f"{cand.selectivity_ratio(reference=ref):.2f}x")
    print("""
Subjective drive per unit respiratory burden is a RATIO, so the unknown lumped sensitivity
factor cancels. It survives 20,000-draw propagation over every estimated parameter at the
5th percentile. Use it to RANK candidates.

It bounds nothing: no safety margin, no overdose multiple, no dose. And there is no
pharmacokinetics anywhere in this package -- "dose" means receptor occupancy.""")

    n_void = len(rs.by_tier(Tier.VOID))
    print(f"\n{len(rs.quotable())} quantities quotable, {n_void} withheld.")


if __name__ == "__main__":
    main()
