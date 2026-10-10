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
import math
from dataclasses import dataclass, field
import numpy as np

from .cpg import Drug
from .subtypes import GabaProfile, PROFILES
from .results import Quantity, ResultSet, Tier
from .config import (RESP_OP, SYNAPTIC_PULSE, AMBIENT_GABA_UM,
                     BRAINSTEM_GLUN2B, FOREBRAIN_GLUN2B)


# An s_max at or above this marks a DIRECT orthosteric agonist rather than an allosteric
# modulator: PROFILES uses 1e9 to mean "no efficacy ceiling", because an agonist does not
# require endogenous GABA and so never saturates the way a PAM does.
DIRECT_AGONIST_CEILING = 1e6


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
    # "auto" tries affinity (the BZ-site default) then gating (neurosteroid-like, which
    # raises maximal current and so reaches shifts affinity alone cannot). Force one with
    # "affinity" or "gating" when the mechanism is known.
    modality: str = "auto"
    nmda_block: float = 0.0
    glun2b_sel: float = 1.0
    glyr: float = 1.0
    _gains: dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_profile(cls, key: str, **kw) -> "Compound":
        """Build from a named profile, CARRYING ITS CEILING into `s_max`.

        `s_max` used to default to 2.5 for every profile, discarding each compound's own
        measured efficacy ceiling: imepitoin (1.25) was modelled at twice its real ceiling
        and the neurosteroid arm (6.0) at under half of its. Since `s_max` is what sets
        achievable effect and overdose protection, that silently flattened the distinction
        between a low-efficacy partial modulator and a high-efficacy one -- which is the
        distinction the project's safety argument turns on.
        """
        p = PROFILES[key]
        kw.setdefault("s_max", p.ceiling)
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
        # MODALITY IS IN THE KEY. Without it, setting `modality` after a first call
        # returned the earlier mechanism's gains from cache -- an affinity-type result
        # served for a gating-type request, silently (roadmap P0-12).
        key = (self.s_max, ambient_um, self.modality, tuple(sorted(pulse.items())))
        if key not in self._gains:
            from .gabaa_kinetics import fit_scheme, derive, calibrate_pam
            s = fit_scheme(verbose=False, pulse=pulse)
            # MODALITY. Not every compound in PROFILES is an affinity-type modulator, and
            # assuming so crashed two of them. Direct orthosteric agonists carry
            # ceiling=1e9 (they need no endogenous GABA, so they do not saturate the way a
            # PAM does) and neurosteroids act on GATING, raising maximal current rather
            # than apparent affinity.
            #
            # This was introduced by the interaction of two earlier fixes: carrying each
            # profile's ceiling into s_max (correct) met a hard raise on unreachable
            # affinity shifts (also correct), and together they made gaboxadol raise
            # instead of evaluate. Neither fix was wrong; the combination was unhandled.
            if self.s_max >= DIRECT_AGONIST_CEILING:
                raise ValueError(
                    f"{self.name!r} has s_max={self.s_max:g}, which marks a DIRECT "
                    f"orthosteric agonist, not an allosteric modulator. The Markov "
                    f"scheme's PAM machinery does not apply: an agonist opens receptors "
                    f"without endogenous GABA, so it has no affinity-shift ceiling and no "
                    f"tonic/phasic asymmetry of the kind modelled here. Model it as a "
                    f"standing conductance instead -- scripts/predict_muscimol.py does "
                    f"exactly that, and for an agonist the drug-modulatable fraction is "
                    f"1.0 by construction since every GABA-A receptor has the orthosteric "
                    f"site.")

            kind = self.modality
            if kind == "auto":
                # affinity first (the BZ-site default); fall back to gating, which raises
                # maximal current and so reaches shifts affinity alone cannot
                kind = "affinity"
                probe = calibrate_pam(s, self.s_max, "affinity",
                                      pulse=pulse, ambient_um=ambient_um)
                if not np.isfinite(probe):
                    kind = "gating"
            aff_kw = {kind: calibrate_pam(s, self.s_max, kind,
                                          pulse=pulse, ambient_um=ambient_um)}
            aff = aff_kw[kind]
            # calibrate_pam returns NaN when the requested EC50 shift is UNREACHABLE.
            # Unchecked, that NaN becomes koff=NaN, then NaN conductances, then NaN
            # voltages inside the LIF integrator -- surfacing as "SVD did not converge in
            # Linear Least Squares" and raw LAPACK complaints far from the cause.
            #
            # The check now lives in gabaa_kinetics.require_reachable so this module and
            # cpg.Drug.from_kinetics share one implementation; cpg had no guard at all,
            # which is recurring error E12 (roadmap P0-8). `kind` names whichever mechanism
            # was tried last, since "auto" has already fallen back by this point.
            from .gabaa_kinetics import require_reachable
            aff = require_reachable(aff, self.s_max, ambient_um, kind)
            aff_kw[kind] = aff
            self._gains[key] = derive(s, ambient_um=ambient_um, pulse=pulse, **aff_kw)
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



def _pct_of_control(value, control):
    """value as a percentage of control, NaN-safe.

    `100.0 * v / max(1e-9, control)` was wrong in a way that manufactures an enormous
    number from a missing one: NaN comparisons are False, so max(1e-9, nan) returns 1e-9
    and a perfectly ordinary value divided by it becomes ~1.5e11 percent. A control reflex
    gain IS legitimately NaN when the Ia response did not exceed noise, so this path is
    reachable, and an undefined ratio must stay undefined.
    """
    if value is None or control is None:
        return float("nan")
    if not (np.isfinite(value) and np.isfinite(control)) or abs(control) < 1e-9:
        return float("nan")
    return 100.0 * value / control


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
# Keyed on (endpoint, seed, SUBSTRATE). If a future change makes a control depend on
# sensitivity, this cache becomes wrong — tests/test_evaluate.py pins the invariance.
#
# THE SUBSTRATE IS IN THE KEY, added with roadmap link 4. The provable-safety argument above
# covers the DRUG parameters only; the neuron model is a different axis entirely, and it
# changes the drug-free control completely (the conductance network's control runs at a
# different frequency and mean output by construction, since it is anchored separately).
# Keyed on seed alone, a conductance arm would have been normalised against a LIF control,
# making every reported "percent of control" silently nonsense -- and the substrate
# comparison is precisely a comparison of those percentages, so the one number the upgrade
# exists to produce would have been the one corrupted.
_CTRL_CACHE: dict = {}


def clear_control_cache() -> None:
    """Drop cached drug-free controls. Call after changing circuit configuration."""
    _CTRL_CACHE.clear()


# THREE THINGS CHANGE TOGETHER WITH THE SUBSTRATE, and getting any one of them wrong is
# silent. Kept in one table so they cannot drift apart.
#
#   operating point  -- pA and nS are relative to the cell (C=200 pF/g_L=10 nS vs
#                       C=21 pF/g_L=2.8 nS), so the LIF's RESP_OP must not be reused.
#   validity band    -- EUPNOEA_BAND is an IN VIVO rat band; the Butera cell is neonatal
#                       rodent IN VITRO, where control burst frequency is ~0.11-0.24 Hz.
#                       Gating a cond run against the in vivo band reports every healthy
#                       rhythm as dead (frequency out of band).
#   settling time    -- tau_h is 10 s on the conductance cell, 25x the LIF's tau_adapt of
#                       400 ms. A 4 s warm-up is 0.4 time constants and measures a decaying
#                       transient: mean h falls 0.265 -> 0.080 over the first 20 s while the
#                       frequency reading looks perfectly stable. That produced a confident,
#                       seed-verified, WRONG operating point once already.
_SUBSTRATE = {
    "lif":  dict(band=None, warm_ms=4000.0,  duration_ms=14000.0),
    "cond": dict(band="invitro", warm_ms=30000.0, duration_ms=60000.0),
}


def _simulate_resp(cand: Compound, seed=0, duration_ms=None, warm_ms=None,
                   substrate="lif"):
    """Simulate the respiratory arm. `substrate` selects the neuron model (roadmap link 4).

    `duration_ms` and `warm_ms` default PER SUBSTRATE (see `_SUBSTRATE`): the conductance
    cell's slowest timescale is 25x the LIF's, so it needs a far longer warm-up. Explicit
    values still override, but the defaults must not be shared.
    """
    from .config import INVITRO_BAND
    from .resp import PreBotC, resp_metrics
    cfg = _SUBSTRATE[substrate]
    warm_ms = cfg["warm_ms"] if warm_ms is None else warm_ms
    duration_ms = cfg["duration_ms"] if duration_ms is None else duration_ms
    band = INVITRO_BAND if cfg["band"] == "invitro" else None

    st, sp = cand.sens("prebotc")
    kw = dict(**RESP_OP) if substrate == "lif" else dict(substrate="cond")
    b = PreBotC(drug=cand.drug("prebotc"), gaba_sens_tonic=st, gaba_sens_phasic=sp,
                seed=seed, **kw)
    for i in range(int(duration_ms / 0.1)):
        b.step(0.1)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > warm_ms
    return resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m], band=band)


def evaluate(cand: Compound, n_seed=4, reference="nonselective_bz",
             include_motor=True, substrate="lif") -> ResultSet:
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
    # A NaN or infinite ratio must NOT carry the VALIDATED tier. It arises when a compound
    # has essentially no respiratory burden, so the denominator vanishes -- an undefined
    # quantity, not a validated one. Printing "selectivity_ratio: nan x [VALIDATED]" would
    # be precisely the kind of unearned authority the tier system exists to prevent.
    ratio_ok = isinstance(ratio, float) and math.isfinite(ratio)
    rs.add(Quantity(
        name="selectivity_ratio", _value=ratio,
        tier=(Tier.VALIDATED if ratio_ok else Tier.VOID), units="x",
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
    gaba_term = subj_t * (g["tonic"] - 1.0) + subj_p * (g["phasic"] - 1.0)
    rs.add(Quantity("subjective_index_gaba", gaba_term, Tier.UNCALIBRATED,
                    provenance="algebra over discrimination-literature weights, NOT a "
                               "simulation; the forebrain reference sensitivity is a unit "
                               "convention so the units are arbitrary",
                    promote_by="not promotable by simulation — drug discrimination "
                               "requires an animal that learns",
                    caveats=("only ratios between compounds are interpretable",)))

    # NMDA ARM OF THE SUBJECTIVE EFFECT. This was dropped when evaluation moved out of the
    # old simulator, which silently under-reported any compound with an NMDA component --
    # and the two-arm design is the whole point, since ethanol's discriminative stimulus is
    # a COMPOUND stimulus (GABA-A positive modulator + NMDA antagonist) and mixtures of the
    # two generalise to ethanol in trained rats.
    #
    # Reported SEPARATELY and the total marked VOID, rather than quietly summed, because
    # the two arms are in INCOMMENSURABLE UNITS: the GABA term is a conductance-weighted
    # quantity and the NMDA term is a blocked-receptor fraction divided by an invented
    # scale. Summing them produces a number, and comparing them -- which is what the
    # project's "GABA salience >= NMDA salience" substitution constraint does -- compares
    # quantities with no common unit. That constraint is what selects the recommended
    # ratio, so the defect is load-bearing and must not be hidden behind a plausible total.
    if cand.nmda_block > 0:
        fb_nmda = cand.nmda_block * (FOREBRAIN_GLUN2B * cand.glun2b_sel
                                     + (1.0 - cand.glun2b_sel))
        rs.add(Quantity("subjective_index_nmda_raw", fb_nmda, Tier.UNCALIBRATED,
                        provenance="forebrain NMDA block fraction, GluN2B-weighted. A "
                                   "receptor occupancy, NOT in the same units as the GABA "
                                   "term above.",
                        promote_by="an ethanol dose-substitution curve that places both "
                                   "arms on one scale"))
        rs.add(Quantity(
            "subjective_index_total", float("nan"), Tier.VOID,
            provenance="the GABA and NMDA arms are in incommensurable units (a "
                       "conductance-weighted quantity vs a receptor-occupancy fraction "
                       "over an invented scale), so their sum is not a quantity. The "
                       "project's substitution constraint 'GABA salience >= NMDA salience' "
                       "compares these two directly and therefore rests on the same "
                       "defect -- and that constraint is what selects the recommended "
                       "ratio.",
            promote_by="anchor both arms to ethanol dose-substitution data so they share "
                       "a scale; do not sum them before that exists"))

    vent = [_simulate_resp(cand, seed=s, substrate=substrate) for s in range(n_seed)]
    ctrl = []
    for s in range(n_seed):
        k = ("resp", s, substrate)
        if k not in _CTRL_CACHE:
            _CTRL_CACHE[k] = _simulate_resp(Compound("control", occupancy=0.0), seed=s,
                                            substrate=substrate)
        ctrl.append(_CTRL_CACHE[k])
    pct = _pct_of_control(float(np.mean([v["mean"] for v in vent])),
                      float(np.mean([c["mean"] for c in ctrl])))
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
    #
    # NO EARLY RETURN. This branch used to `return rs` when `assays` would not import,
    # which skipped the VOID block below -- so the result SHAPE depended on the install,
    # in a package whose whole thesis is that the VOID entries always appear.
    # tests/test_evaluate.py asserts `overdose_index` is present, and on a lean install it
    # was present only because `assays` imports mujoco lazily (roadmap P0-10).
    motor_available = True
    if include_motor:
        try:
            from .assays import stretch_reflex, locomotion
        except Exception as e:                       # pragma: no cover
            motor_available = False
            rs.add(Quantity("motor_endpoints", None, Tier.UNCALIBRATED,
                            provenance=f"body plant unavailable ({type(e).__name__}); "
                                       "install the 'plant' extra for motor endpoints"))
    if include_motor and motor_available:
        stm, spm = cand.sens("spinal")
        d = cand.drug("spinal")
        try:
            # AVERAGED OVER SEEDS. These were single-seed point estimates (seed=1 hard
            # coded) while ventilation averaged over n_seed -- recurring error E6, which
            # this project catalogued and then committed in its own public API. Motor
            # readouts are noisier than ventilation, so a single seed is worse here.
            seeds = range(n_seed)
            refl_runs = [stretch_reflex(d, gaba_sens_tonic=stm, gaba_sens_phasic=spm,
                                        seed=s) for s in seeds]
            for s in seeds:
                if ("reflex", s) not in _CTRL_CACHE:
                    _CTRL_CACHE[("reflex", s)] = stretch_reflex(Drug(), seed=s)
            ctrl_refl_runs = [_CTRL_CACHE[("reflex", s)] for s in seeds]
            refl = {"gain": float(np.mean([r["gain"] for r in refl_runs])),
                    "mn_dyn": float(np.mean([r["mn_dyn"] for r in refl_runs]))}
            ctrl_refl = {"gain": float(np.mean([r["gain"] for r in ctrl_refl_runs]))}
            rs.add(Quantity(
                "reflex_gain", _pct_of_control(refl["gain"], ctrl_refl["gain"]),
                Tier.UNCALIBRATED, "% of control",
                provenance="ramp-and-hold stretch reflex, mirroring the servo-imposed "
                           "experiment; reproduces strychnine hyperreflexia and "
                           "benzodiazepine depression without being fitted to them",
                promote_by="a dose-response anchor, e.g. diazepam H-reflex depression in a "
                           "preparation matching the modelled arc. The spinal sensitivity "
                           "was never quantitatively fitted — it was set where the "
                           "qualitative validations passed."))

            loco_runs = [locomotion(d, gaba_sens_tonic=stm, gaba_sens_phasic=spm,
                                    seed=s) for s in seeds]
            for s in seeds:
                if ("loco", s) not in _CTRL_CACHE:
                    _CTRL_CACHE[("loco", s)] = locomotion(Drug(), seed=s)
            ctrl_loco_runs = [_CTRL_CACHE[("loco", s)] for s in seeds]
            def _m(runs, k):
                """Mean over seeds, NaN-safe without the warning.

                np.nanmean over an ALL-NaN list emits RuntimeWarning: Mean of empty slice
                and returns NaN. That happens legitimately -- a high-dose sedative can
                abolish locomotion in every seed, so step_period_ms and alternation are
                NaN throughout -- so the warning is noise on a correct result, and noise
                on correct results is how real warnings get ignored.
                """
                vals = [r[k] for r in runs
                        if isinstance(r[k], (int, float)) and np.isfinite(r[k])]
                return float(np.mean(vals)) if vals else float("nan")
            loco = {k: _m(loco_runs, k) for k in
                    ("step_period_ms", "alternation", "duty", "excursion_rad_INVALID")}
            loco["walking"] = bool(np.mean([r["walking"] for r in loco_runs]) >= 0.5)
            ctrl_loco = {"step_period_ms": _m(ctrl_loco_runs, "step_period_ms")}
            rs.add(Quantity(
                "step_period", _pct_of_control(loco["step_period_ms"], ctrl_loco["step_period_ms"]),
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


        # UNCONDITIONAL, and outside the try/except on purpose. This Quantity records that
        # `joint_excursion` is permanently invalid (wrong sign), which is a property of the
        # METRIC, not of whether the body plant happens to be installed. It used to sit
        # inside the try, after the assay calls, so on a lean install (`.[dev]`, no plant --
        # i.e. exactly what CI's lean jobs use) the ImportError jumped past it and the
        # quantity was simply absent: `rs.quantity("joint_excursion")` raised KeyError with
        # mujoco missing and returned a VOID quantity with it present.
        #
        # An environment-dependent result SCHEMA is the problem, not the missing number. A
        # caller cannot write `if rs.quantity("joint_excursion").tier is Tier.VOID` and have
        # it mean the same thing on two machines, and the VOID tier exists precisely so that
        # an invalid metric is visible rather than absent.
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
