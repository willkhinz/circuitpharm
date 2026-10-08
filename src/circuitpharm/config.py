"""Single source of truth for operating points and calibrations.

WHY. These values were previously duplicated as literals across at least five scripts
(`simulator.py`, `subtype_sweep.py`, `calib_*.py`, `overdose_kinetic.py`,
`recalibrate_kinetic.py`). Duplicated constants are how a calibration gets carried across
a change of parameterisation without being refitted -- recurring error E12 in this project,
committed twice, the second time one worklog entry after writing down that the migration
would be needed. A single definition with its provenance attached makes that mistake
visible instead of silent.

Every calibration here carries: the value, the observable it was fitted to, and whether
that fit is still VALID. Read `CALIBRATIONS` before trusting any absolute number.
"""
from __future__ import annotations
from dataclasses import dataclass
from types import MappingProxyType

from .results import Tier


# ---------------------------------------------------------------- operating points
# preBotC operating point. `drive` MUST stay sub-rheobase (~200 pA at gL=10 nS with tonic
# GABA 2 nS): I_NaP carries the cell over threshold, so adaptation and mutual inhibition
# can switch a burst off. Supra-rheobase drive makes the network fire tonically and renders
# every weight sweep flat -- recurring error E1, committed twice.
RESP_OP = MappingProxyType(dict(
    drive=170.0,            # pA, sub-rheobase
    g_adapt=2.5,            # nS per spike, burst termination
    tau_adapt=400.0,        # ms
    w=MappingProxyType(dict(ee_ampa=0.45, ee_nmda=0.2475)),
))

# CONDUCTANCE-SUBSTRATE RESPIRATORY OPERATING POINT (roadmap link 5).
#
# None means NOT YET ANCHORED, and `PreBotC(substrate="cond")` raises rather than falling
# back to RESP_OP above. That refusal is the point: RESP_OP's `drive` is in pA and its
# weights in nS, both relative to the LIF cell (C=200 pF, g_L=10 nS). The Butera cell is
# C=21 pF, g_L=2.8 nS, so reusing those numbers is an order-of-magnitude error that produces
# a running network rather than an exception -- predicted as E14.
#
# WHY THIS IS ANCHORED RATHER THAN INHERITED. Butera-Rinzel-Smith part II supplies network
# coupling conductances for a population of these cells, which would have made these weights
# published rather than ours. They are not retrievable: the journal full text returns HTTP
# 403 and every accessible encoding (CellML, ModelDB 247647) is single-cell only. Inventing
# them was not an option, so the weights here are OURS, anchored to a matched operating
# point. The CELL is published; the COUPLING is not, and link 5 is partial for that reason.
#
# Matched operating point, not matched parameters: the same nS weight means different things
# on the two cells, so comparing substrates at equal weights compares two differently-broken
# networks. Each substrate is anchored independently to the SAME OBSERVABLE (control burst
# frequency inside EUPNOEA_BAND, with a comparable duty cycle), and only then is the drug
# applied and the fractional change from each substrate's own control compared.
# CONDUCTANCE-SUBSTRATE RESPIRATORY OPERATING POINT (roadmap link 5).
#
# None means NOT ANCHORED, and `PreBotC(substrate="cond")` raises rather than falling back
# to RESP_OP above -- whose `drive` is in pA and weights in nS relative to the LIF cell
# (C=200 pF, g_L=10 nS) rather than the Butera cell (C=21 pF, g_L=2.8 nS).
#
# WHY THIS IS None AGAIN, having briefly been a registered operating point (2026-10-07).
# Recorded in full because it is the clearest example of this project's signature failure
# that the project has produced, and I produced it in the same session that wrote the
# warning.
#
# An operating point was anchored (drive=20 pA, weights 0.9x the LIF table), verified at
# 30 s across 4 seeds with every seed required alive, and read 1.342 +/- 0.132 Hz against
# the LIF's 1.271 Hz. It was wrong, for two compounding reasons:
#
#   1. THE SEARCH NEVER VISITED THE BURSTING REGIME. The grid swept Exc drive over
#      {5, 10, 15, 20, 25, 30} pA. The isolated Butera cell bursts at -5..0 pA and is
#      TONIC from +5 pA upward. So every point searched held the cells depolarised out of
#      pacemaking, and "only 3 of 50 points were alive" -- which I recorded as the correct
#      biology of a narrow entrainment window -- was really the signature of searching
#      almost entirely outside the regime the cell can oscillate in.
#
#   2. THE MEASUREMENT WAS INSIDE A TRANSIENT. tau_h is 10 SECONDS, the slowest timescale
#      in the model by a factor of 25 against the LIF's tau_adapt = 400 ms. Every
#      conductance measurement used a 4 s warm-up, i.e. 0.4 time constants. Measured in
#      successive 10 s windows, the network's mean h falls 0.265 -> 0.080 over the first
#      20 s and the rhythm then wanders: frequency 1.33, 1.33, 2.86, 2.86, 3.67, 0.82,
#      0.31, 0.82, 0.61 Hz, alive in only 5 of 9 windows. The LIF, measured identically,
#      holds 1.327 Hz and mean 27-29 in every window.
#
# So the registered number was a plausible reading of a decaying transient, produced by a
# three-stage search with a seed-robustness check -- none of which could see the problem,
# because all three stages measured inside the same transient.
#
# WHAT A VALID ANCHOR MUST DO, enforced in scripts/anchor_cond_resp.py:
#   * search Exc drive in the range where the isolated cell BURSTS (verify with
#     neuron.run_isolated before choosing the grid, do not assume);
#   * warm up for at least 3 x tau_h = 30 s before measuring anything;
#   * measure two consecutive late windows and require the frequency to AGREE between
#     them, so a drifting network cannot pass as a stable one;
#   * require every seed alive, as before -- necessary but, as this showed, nowhere near
#     sufficient.
COND_RESP_OP = None

# Synaptic GABA transient seen by the PHASIC receptor pool. This particular (peak,
# clearance) pair is one of the 9/27 cells in which the Markov scheme reproduced every
# anchored benzodiazepine observable from a single calibrated parameter; the passing cells
# all sat at a saturating cleft, which the data selected rather than us imposing it.
SYNAPTIC_PULSE = MappingProxyType(dict(peak_um=3000.0, clear_ms=1.00))

# Ambient extrasynaptic GABA seen by the TONIC pool, uM. Literature range ~0.1-1.
# This value materially changes the derived efficacy ceiling: tonic headroom is ~3295x at
# 0.1 uM, ~211x at 0.4 uM and ~1.5x at 10 uM, so it is not a cosmetic parameter.
AMBIENT_GABA_UM = 0.40

# Ia -> Mn weight scale used in the stretch-reflex assay. Needed because the motoneuron
# pool's dynamic peak otherwise pins at the tref=8 ms firing ceiling (125 Hz), which makes
# every drug read ~100% of control -- recurring error E5.
IA_SCALE = 0.30

# GluN2B fraction of local NMDA receptors. The asymmetry between these two numbers IS the
# NMDA selectivity window: forebrain is GluN2B-rich, brainstem GluN2D-dominant.
FOREBRAIN_GLUN2B = 0.70
BRAINSTEM_GLUN2B = 0.15

# Physiologically admissible respiratory frequency, Hz (rat eupnoea 1-2 Hz). A hard gate on
# the `alive` flag: under heavy block the burst train fragments and the FFT correctly
# reports ~4 Hz, which is not tachypnoea but a disintegrated rhythm. Because `alive` is
# what the overdose scan keys on, omitting this gate made overdose indices optimistic.
EUPNOEA_BAND = (0.30, 2.50)


# ----------------------------------------------------------------- calibrations
@dataclass(frozen=True)
class Calibration:
    name: str
    value: float
    observable: str          # what it was fitted to
    tier: Tier
    status: str              # plain-language validity
    promote_by: str = ""

    def __str__(self) -> str:
        return (f"{self.name} = {self.value!r}  [{self.tier.label}]\n"
                f"    fitted to:  {self.observable}\n"
                f"    status:     {self.status}"
                + (f"\n    promote by: {self.promote_by}" if self.promote_by else ""))


CALIBRATIONS = {
    "prebotc_gaba_sens": Calibration(
        name="prebotc_gaba_sens",
        value=0.15,
        observable="human midazolam 2 mg IV, minute ventilation -14.3 +/- 5.9% "
                   "(second cohort -19 +/- 7%)",
        tier=Tier.VOID,
        status="INVALID. Whole-body ventilation is the wrong observable for an isolated "
               "preBotC. Much of benzodiazepine respiratory depression is chemoreflex "
               "blunting plus upper-airway and cortical drive, none of which exist in this "
               "circuit, so this parameter was a lumped factor absorbing mechanisms the "
               "model does not contain -- not the subunit fraction its name implies. It "
               "therefore cannot be extrapolated to a compound with a different subtype "
               "profile, which is exactly what it was being used for. Separately it is "
               "unidentified: 18 of 35 grid cells fit the anchor equally well, spanning "
               "gaba_sens 0.06-0.30.",
        promote_by="a structure-matched anchor: muscimol or GABA concentration-response on "
                   "inspiratory burst frequency AND amplitude in a rhythmic preBotC slice "
                   "or perfused preparation, in animals older than P12 (the preBotC "
                   "NKCC1/KCC2 crossover) so GABA-A is reliably inhibitory. A direct "
                   "orthosteric agonist pins this parameter at 1.0 by construction, "
                   "leaving ZERO free parameters -- see predict_muscimol.",
    ),
    "spinal_gaba_sens": Calibration(
        name="spinal_gaba_sens",
        value=1.00,
        observable="the value at which the four stretch-reflex phenotype checks passed "
                   "(strychnine hyperreflexia, benzodiazepine depression, GluN2B sparing)",
        tier=Tier.UNCALIBRATED,
        status="Never quantitatively fitted. Set where qualitative validations passed, "
               "which constrains the ordering but not the scale.",
        promote_by="a dose-response anchor for reflex depression, e.g. diazepam H-reflex "
                   "dose-response in a preparation matching the modelled arc",
    ),
    "forebrain_gaba_sens": Calibration(
        name="forebrain_gaba_sens",
        value=1.00,
        observable="nothing — reference region, defined as 1.0",
        tier=Tier.UNCALIBRATED,
        status="A unit convention, not a measurement. Consequently the subjective index "
               "has arbitrary units and a 'subjective target' of 0.50 has no physical "
               "meaning; only RATIOS between compounds are interpretable.",
        promote_by="not promotable by simulation. The subjective axis requires drug "
                   "discrimination in trained animals, which needs learning and is "
                   "therefore outside what this architecture can represent.",
    ),
    "gabaa_kinetics_baseline": Calibration(
        name="gabaa_kinetics_baseline",
        value=float("nan"),
        observable="GABA EC50 ~20 uM, max open probability ~0.75, IPSC deactivation tau "
                   "~15 ms (all three reproduced exactly by the fit)",
        tier=Tier.UNCALIBRATED,
        status="Structurally sound and its PHASIC predictions are externally corroborated "
               "(midazolam 1 uM prolongs IPSC decay ~1.5-2x without changing amplitude; "
               "the scheme predicts tau 1.66x, peak 1.06x from one parameter). BUT the "
               "three fit targets were taken from recalled literature RANGES, not "
               "digitised from specific published figures, so no quantitative claim here "
               "is citable.",
        promote_by="digitise the three target curves from named sources and refit; the "
                   "scheme itself does not change",
    ),
    "ranking": Calibration(
        name="selectivity_ranking",
        value=float("nan"),
        observable="nothing — it is calibration-INDEPENDENT by construction",
        tier=Tier.VALIDATED,
        status="Subjective drive per unit respiratory burden is a RATIO, so the unknown "
               "lumped sensitivity scale cancels. Survives 20,000 draws over every "
               "invented parameter (regional subunit fractions, extrasynaptic fractions "
               "with eps sampled uniform on [0,1] because it was flagged a guess, and the "
               "tonic:phasic gain ratio): all three a5 arms clear a non-selective "
               "benzodiazepine at the 5th percentile (>=2.49x). Ranks candidates; bounds "
               "nothing.",
    ),
}


def calibration_report() -> str:
    out = ["CALIBRATION REGISTRY", "=" * 72]
    for tier in (Tier.VALIDATED, Tier.UNCALIBRATED, Tier.VOID):
        group = [c for c in CALIBRATIONS.values() if c.tier is tier]
        if not group:
            continue
        out.append(f"\n-- {tier.label} " + "-" * (68 - len(tier.label)))
        out.extend(str(c) for c in group)
    out.append("\nNOTE: the respiratory axis currently has NO valid quantitative anchor. "
               "Absolute\nventilation numbers are not quotable; the selectivity ranking is.")
    return "\n".join(out)


if __name__ == "__main__":
    print(calibration_report())
