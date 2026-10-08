"""Public API: evaluate a compound's receptor profile and get back TIERED results.

This replaces the ad-hoc evaluation that lived in `simulator.py`, which had drifted out of
date in three ways that all flattered the model:

  1. it passed a SINGLE `gaba_sens` to each circuit, so one gain was applied to both the
     tonic and phasic receptor pools although an affinity-type PAM potentiates them by
     ~7.2x and ~1.06x respectively -- wrong by about sevenfold on whichever pool it had not
     been calibrated against;
  2. it clipped dose escalation at a hand-set `ceiling` of 2.5, a number that is only
     correct at an ambient GABA roughly ten times physiological, and which conflates a
     ligand property (intrinsic allosteric efficacy) with a synapse property (distance from
     agonist saturation);
  3. it printed an overdose index and a ventilation percentage as plain numbers with a
     warning block beneath them. The warning was correct and was ignored: a configuration
     the model itself flagged UNREACHABLE was summarised elsewhere as a clean result.

So evaluation now returns a `ResultSet` whose VOID entries cannot be read without an
explicit acknowledgement, and whose printed form withholds them entirely.

DOSE IS OCCUPANCY, NOT MASS. There is no pharmacokinetics anywhere in this package: no
mg/kg, no brain:plasma ratio, no time course. `occupancy` is the fraction of modulator
sites bound, and the pool gains are exactly linear in it because a partially occupied
population is a MIXTURE of modulated and unmodulated receptors. They therefore saturate at
full occupancy, and that saturation -- set by the ligand's own intrinsic allosteric
efficacy `s_max` -- is the only honest ceiling.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

from .cpg import Drug
from .subtypes import GabaProfile, PROFILES
from .results import Quantity, ResultSet, Tier
from .config import (RESP_OP, SYNAPTIC_PULSE, AMBIENT_GABA_UM,
                     BRAINSTEM_GLUN2B, FOREBRAIN_GLUN2B)


# ----------------------------------------------------------------------- compound
@dataclass
class Compound:
    """A compound described by MEASURED receptor activity.

    `s_max` is the ligand's intrinsic allosteric efficacy: the maximum factor by which it
    shifts the GABA dose-response at full site occupancy. ~2-3x for classical
    benzodiazepine-site ligands. This is the parameter that actually determines overdose
    protection, it is a property of the individual molecule, and it must be MEASURED --
    it cannot be inherited from the PAM class or predicted from structure.
    """
    name: str
    a1: float = 0.0
    a23: float = 0.0
    a5: float = 0.0
    d_a4: float = 0.0
    eps: float = 0.0
    s_max: float = 2.5             # intrinsic allosteric efficacy (max EC50 shift)
    occupancy: float = 0.5         # fraction of modulator sites bound
    nmda_block: float = 0.0
    glun2b_sel: float = 1.0
    glyr: float = 1.0
    _gains: dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_profile(cls, key: str, **kw) -> "Compound":
        p = PROFILES[key]
        return cls(name=p.name, a1=p.a1, a23=p.a23, a5=p.a5,
                   d_a4=p.d_a4, eps=p.eps, **kw)

    @property
    def gaba(self) -> GabaProfile:
        return GabaProfile(self.name, a1=self.a1, a23=self.a23, a5=self.a5,
                           d_a4=self.d_a4, eps=self.eps)

    # ---- kinetics: one mechanistic parameter -> both pool gains --------------------
    def pool_gains(self, ambient_um=AMBIENT_GABA_UM, pulse=None) -> dict:
        """Derive (tonic, phasic, tau) gains at this occupancy from the Markov scheme.

        Cached, because fitting the scheme costs a few seconds and the result depends only
        on (s_max, ambient, pulse) -- not on occupancy, which enters linearly afterwards.
        """
        pulse = pulse or dict(SYNAPTIC_PULSE)
        key = (self.s_max, ambient_um, tuple(sorted(pulse.items())))
        if key not in self._gains:
            from .gabaa_kinetics import fit_scheme, derive, calibrate_pam
            s = fit_scheme(verbose=False, pulse=pulse)
            aff = calibrate_pam(s, self.s_max, "affinity",
                                pulse=pulse, ambient_um=ambient_um)
            self._gains[key] = derive(s, affinity=aff, ambient_um=ambient_um,
                                      pulse=pulse)
        full = self._gains[key]
        lin = lambda g: 1.0 + self.occupancy * (g - 1.0)
        return dict(tonic=lin(full["tonic_gain"]),
                    phasic=lin(full["phasic_gain"]),
                    tau=lin(full["tau_ratio"]),
                    tonic_at_full=full["tonic_gain"],
                    phasic_at_full=full["phasic_gain"],
                    tonic_headroom=full["tonic_headroom"])

    def drug(self, region="prebotc", **kw) -> Drug:
        g = self.pool_gains(**kw)
        return Drug(gaba_a_gain=g["phasic"],
                    gaba_a_gain_tonic=g["tonic"],
                    gaba_a_tau=g["tau"],
                    # the MECHANISM caps the gain (it saturates at full occupancy), so no
                    # fiat cap is imposed here
                    gaba_a_efficacy_cap=1e9, gaba_a_cap_tonic=1e9,
                    nmda_block=self.nmda_block,
                    glun2b_selectivity=self.glun2b_sel,
                    glun2b_fraction=(BRAINSTEM_GLUN2B if region != "forebrain"
                                     else FOREBRAIN_GLUN2B),
                    glyr_gain=self.glyr)

    def sens(self, region: str) -> tuple:
        """(tonic, phasic) drug-modulatable fraction in this region."""
        return self.gaba.regional_sens_split(region)

    # ---- the calibration-independent quantity -------------------------------------
    def selectivity_ratio(self, reference="nonselective_bz",
                          tonic_phasic_ratio=6.0) -> float:
        """Subjective drive per unit preBotC respiratory burden, relative to `reference`.

        THE ONLY SAFETY-RELEVANT QUANTITY IN THIS PACKAGE THAT IS CALIBRATION-INDEPENDENT.
        Both numerator and denominator scale with the unknown lumped sensitivity factor, so
        it cancels in the ratio. Survives 20,000-draw propagation over every estimated
        parameter at the 5th percentile. It RANKS compounds; it bounds nothing.
        """
        def score(p: GabaProfile) -> float:
            st, sp = p.subjective_index_split()
            bt, bp = p.regional_sens_split("prebotc")
            subj = tonic_phasic_ratio * st + sp
            burden = tonic_phasic_ratio * bt + bp
            return subj / burden if burden > 1e-12 else float("nan")
        return score(self.gaba) / score(PROFILES[reference])


# ------------------------------------------------------------------- simulation
# CONTROL-RUN CACHE. Every evaluation needs drug-free controls to normalise against, and
# recomputing them per call made evaluate() six simulations deep; the test suite grew from
# ~4.5 to ~9.5 minutes as endpoints were added.
#
# Caching them is SAFE FOR A PROVABLE REASON, not an approximation. A drug-free `Drug()`
# has gaba_scale == gaba_scale_tonic == nmda_scale == glyr_gain == 1, and every sensitivity
# enters as eff = 1 + sens*(x - 1). With x == 1 that is 1 for ANY sens, and the decay taus
# scale by the same rule. So the drug-free control is completely independent of
# gaba_sens/tonic/phasic, and depends only on the seed and the circuit configuration.
#
# Keyed on seed alone for that reason. If a future change makes a control depend on
# sensitivity, this cache becomes wrong — tests/test_evaluate.py pins the invariance.
_CTRL_CACHE: dict = {}


def clear_control_cache() -> None:
    """Drop cached drug-free controls. Call after changing circuit configuration."""
    _CTRL_CACHE.clear()


def _simulate_resp(cand: Compound, seed=0, duration_ms=14000.0, warm_ms=4000.0):
    from .resp import PreBotC, resp_metrics
    st, sp = cand.sens("prebotc")
    b = PreBotC(drug=cand.drug("prebotc"), gaba_sens_tonic=st, gaba_sens_phasic=sp,
                seed=seed, **RESP_OP)
    for i in range(int(duration_ms / 0.1)):
        b.step(0.1)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > warm_ms
    return resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])


def evaluate(cand: Compound, n_seed=4, reference="nonselective_bz",
             include_motor=True) -> ResultSet:
    """Evaluate a compound and return tiered results.

    Ventilation is simulated and reported as UNCALIBRATED: the mechanism is sound but the
    respiratory axis currently has NO valid quantitative anchor (its only one, human
    whole-body minute ventilation, is the wrong observable for an isolated preBotC). The
    shape and the ordering between compounds may be used; the absolute scale may not.
    """
    rs = ResultSet(f"{cand.name}  (occupancy {cand.occupancy:.2f}, "
                   f"s_max {cand.s_max:.1f})")

    # --- VALIDATED: the calibration-independent ranking -----------------------------
    ratio = cand.selectivity_ratio(reference=reference)
    rs.add(Quantity(
        name="selectivity_ratio", _value=ratio, tier=Tier.VALIDATED, units="x",
        provenance=f"subjective drive per unit preBotC burden, relative to {reference}; "
                   "a ratio, so the unknown lumped sensitivity factor cancels. Survives "
                   "20,000-draw propagation over every estimated parameter at the 5th pct.",
        caveats=("ranks compounds; bounds nothing — no margin, no dose",)))

    # --- UNCALIBRATED: the receptor-level quantities --------------------------------
    g = cand.pool_gains()
    st, sp = cand.sens("prebotc")
    rs.add(Quantity("gain_tonic", g["tonic"], Tier.UNCALIBRATED, "x",
                    provenance="GABA-A Markov scheme, one parameter calibrated on the "
                               "EC50 left-shift; phasic predictions externally "
                               "corroborated (midazolam prolongs IPSC decay without "
                               "changing amplitude)",
                    promote_by="digitise the three baseline fit targets from named "
                               "sources rather than recalled ranges"))
    rs.add(Quantity("gain_phasic", g["phasic"], Tier.UNCALIBRATED, "x",
                    provenance="same scheme; near-saturating cleft transient leaves little "
                               "headroom, which is why this is ~1 while tonic is ~7"))
    rs.add(Quantity("prebotc_sens_tonic", st, Tier.UNCALIBRATED,
                    provenance="regional subunit composition x extrasynaptic localisation"))
    rs.add(Quantity("prebotc_sens_phasic", sp, Tier.UNCALIBRATED,
                    provenance="regional subunit composition x synaptic localisation"))

    subj_t, subj_p = cand.gaba.subjective_index_split()
    delivered = subj_t * (g["tonic"] - 1.0) + subj_p * (g["phasic"] - 1.0)
    rs.add(Quantity("subjective_index", delivered, Tier.UNCALIBRATED,
                    provenance="algebra over discrimination-literature weights, NOT a "
                               "simulation; the forebrain reference sensitivity is a unit "
                               "convention so the units are arbitrary",
                    promote_by="not promotable by simulation — drug discrimination "
                               "requires an animal that learns",
                    caveats=("only ratios between compounds are interpretable",)))

    vent = [_simulate_resp(cand, seed=s) for s in range(n_seed)]
    ctrl = []
    for s in range(n_seed):
        k = ("resp", s)
        if k not in _CTRL_CACHE:
            _CTRL_CACHE[k] = _simulate_resp(Compound("control", occupancy=0.0), seed=s)
        ctrl.append(_CTRL_CACHE[k])
    pct = 100.0 * np.mean([v["mean"] for v in vent]) / np.mean([c["mean"] for c in ctrl])
    rs.add(Quantity("ventilation", pct, Tier.UNCALIBRATED, "% of control",
                    provenance="preBotC minute-ventilation proxy, tonic/phasic resolved",
                    promote_by="muscimol or GABA concentration-response on burst frequency "
                               "AND amplitude in a rhythmic preBotC slice or perfused prep, "
                               "animals older than P12 (the preBotC NKCC1/KCC2 crossover). "
                               "A direct orthosteric agonist pins the sensitivity at 1.0, "
                               "leaving zero free parameters.",
                    caveats=("the only anchor this axis ever had was invalidated: human "
                             "whole-body ventilation is the wrong observable for an "
                             "isolated preBotC",)))

    alive = float(np.mean([v["alive"] for v in vent]))
    rs.add(Quantity("rhythm_alive_fraction", alive, Tier.UNCALIBRATED,
                    provenance="eupnoea-band gated; a fragmented fast train does not count "
                               "as a live rhythm"))

    # --- motor endpoints, if the body plant is installed ----------------------------
    # These are an optional extra (mujoco), so their absence must not fail an evaluation.
    if include_motor:
        try:
            from .assays import stretch_reflex, locomotion
        except Exception as e:                       # pragma: no cover
            rs.add(Quantity("motor_endpoints", None, Tier.UNCALIBRATED,
                            provenance=f"body plant unavailable ({type(e).__name__}); "
                                       "install the 'plant' extra for motor endpoints"))
            return rs
        stm, spm = cand.sens("spinal")
        d = cand.drug("spinal")
        try:
            refl = stretch_reflex(d, gaba_sens_tonic=stm, gaba_sens_phasic=spm)
            if ("reflex", 1) not in _CTRL_CACHE:
                _CTRL_CACHE[("reflex", 1)] = stretch_reflex(Drug())
            ctrl_refl = _CTRL_CACHE[("reflex", 1)]
            rs.add(Quantity(
                "reflex_gain", 100.0 * refl["gain"] / max(1e-9, ctrl_refl["gain"]),
                Tier.UNCALIBRATED, "% of control",
                provenance="ramp-and-hold stretch reflex, mirroring the servo-imposed "
                           "experiment; reproduces strychnine hyperreflexia and "
                           "benzodiazepine depression without being fitted to them",
                promote_by="a dose-response anchor, e.g. diazepam H-reflex depression in a "
                           "preparation matching the modelled arc. The spinal sensitivity "
                           "was never quantitatively fitted — it was set where the "
                           "qualitative validations passed."))

            loco = locomotion(d, gaba_sens_tonic=stm, gaba_sens_phasic=spm)
            if ("loco", 1) not in _CTRL_CACHE:
                _CTRL_CACHE[("loco", 1)] = locomotion(Drug())
            ctrl_loco = _CTRL_CACHE[("loco", 1)]
            rs.add(Quantity(
                "step_period", 100.0 * loco["step_period_ms"]
                / max(1e-9, ctrl_loco["step_period_ms"]),
                Tier.UNCALIBRATED, "% of control",
                provenance="free-running closed loop: receptor -> circuit -> muscle -> "
                           "joint -> spindle feedback, nothing imposed. Slowing is a valid "
                           "impairment signal.",
                promote_by="gait analysis under a sedative at matched receptor occupancy"))
            rs.add(Quantity(
                "coordination", loco["alternation"], Tier.UNCALIBRATED,
                provenance="flexor/extensor motoneuron correlation; more positive = worse "
                           "antagonist coordination. Degrades correctly under sedation "
                           "(-0.67 control -> -0.34 at a 4x PAM).",
                promote_by="gait analysis; this is the closest thing here to an ataxia "
                           "measure"))
            rs.add(Quantity(
                "walking", loco["walking"], Tier.UNCALIBRATED,
                provenance="whether rhythmic joint movement persists at all"))
            # the metric that looks right and is not
            rs.add(Quantity(
                "joint_excursion", float("nan"), Tier.VOID, "rad",
                provenance="reports the WRONG SIGN: a sedative INCREASES joint excursion "
                           "here (1.607 rad control -> 1.940 at a 4x PAM), because less "
                           "antagonist co-contraction leaves the joint less stiff. The "
                           "mechanism is real but the metric is an artefact of a "
                           "single-joint preparation with the body fixed, no gravitational "
                           "load and no ground contact — in an animal, lost co-contraction "
                           "presents as instability, and this model has nothing to "
                           "collapse against.",
                promote_by="a whole-body preparation with ground reaction forces; this "
                           "cannot be fixed by reweighting the metric"))
        except (ImportError, FileNotFoundError) as e:
            # ImportError MUST be caught HERE, not only around the import of the assay
            # functions above. `assays` imports mujoco lazily INSIDE each function, so the
            # import above succeeds on a lean install and the failure surfaces on the
            # first CALL -- which is outside the earlier guard. That made evaluate() crash
            # on exactly the install CI uses (`.[dev]`, no plant extra), while passing
            # locally where mujoco happens to be present.
            rs.add(Quantity(
                "motor_endpoints", None, Tier.UNCALIBRATED,
                provenance=f"body plant unavailable ({type(e).__name__}: {e}). The "
                           "respiratory and receptor-level results above are unaffected.",
                promote_by="pip install -e '.[plant]' for the reflex and locomotor "
                           "endpoints"))

    # --- VOID: quantities that rest on something invalid ----------------------------
    rs.add(Quantity(
        "overdose_index", float("nan"), Tier.VOID, "x",
        provenance="the efficacy ceiling it would clip against is pool-dependent "
                   "(synaptic headroom ~1.1x, extrasynaptic ~211x), and the project's two "
                   "core assumptions — a tight ceiling as the sole overdose protection AND "
                   "full subjective delivery by an a5-selective PAM — are the same "
                   "parameter and mutually exclusive. The dose-escalation formulation is "
                   "additionally incomplete: benzodiazepine potentiation is non-monotonic "
                   "(bell-shaped, antagonistic above ~10 uM) and this scheme is monotonic.",
        promote_by="measure the ligand's intrinsic allosteric efficacy by patch clamp, and "
                   "give the scheme a non-monotonic concentration-response. Do not refit."))

    if cand.nmda_block > 0:
        rs.add(Quantity(
            "nmda_respiratory_contribution", float("nan"), Tier.VOID, "% of control",
            provenance="established negative result: this preBotC ties burst maintenance "
                       "to the long (100 ms) NMDA conductance, so ANY block depresses it. "
                       "Clinically ketamine-class blockers PRESERVE or stimulate breathing. "
                       "Adding the missing disinhibition pathway (ei_nmda) did not help, "
                       "and removing 100% of synaptic inhibition recovers only ~10 of a "
                       "~55 point deficit.",
            promote_by="take the NMDA arm's respiratory contribution from clinical data; "
                       "this circuit cannot produce it"))
    return rs
