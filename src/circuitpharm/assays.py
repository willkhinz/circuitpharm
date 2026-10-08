"""Behavioural assays: protocols that map a circuit onto a measurable animal readout.

WHY THESE ARE THE PACKAGE'S BEHAVIOURAL ENDPOINTS. The chain from a compound to arbitrary
behaviour is not closed -- controllers that produce naturalistic rodent behaviour are
artificial networks whose units have no conductances, so there is nothing to inject a
receptor current into. What IS closed is the chain to behaviours whose generating circuits
are biophysically modellable, and each assay here deliberately mirrors a real protocol so
the model's output and an animal's output are the same KIND of measurement:

    stretch_reflex   ramp-and-hold with imposed kinematics   <-> H-reflex / tendon jerk
    (respiratory metrics live in resp.py)                    <-> plethysmography

MODEL PATH RESOLUTION. `MODEL_XML` was previously the relative string
"models/rodent_muscle.xml", which silently broke whenever anything ran from outside the
repository root -- fine while this was a pile of scripts, not fine for an installed
package. It now resolves against the package location with a repo-root fallback.
"""
from __future__ import annotations
import os
from pathlib import Path

import numpy as np

from .cpg import Drug


def _find_model(name="rodent_muscle.xml") -> str:
    """Locate the MuJoCo body model.

    Searched in order: an explicit CIRCUITPHARM_MODELS override, a `models/` directory
    beside the installed package, then walking up from the package towards a repository
    root. Raising a clear error beats a MuJoCo parse failure on a path that does not exist.
    """
    env = os.environ.get("CIRCUITPHARM_MODELS")
    cands = []
    if env:
        cands.append(Path(env) / name)
    here = Path(__file__).resolve()
    cands.append(here.parent / "models" / name)
    for up in list(here.parents)[:5]:
        cands.append(up / "models" / name)
    for c in cands:
        if c.is_file():
            return str(c)
    raise FileNotFoundError(
        f"could not locate {name}. Looked in: "
        + ", ".join(str(c) for c in cands)
        + ". Set CIRCUITPHARM_MODELS to the directory holding it, or build it with "
          "scripts/port_muscle.py (the stock dm_control rodent has position servos, not "
          "muscle actuators, and is unusable for pharmacology -- see that script)."
    )


# resolved lazily so importing this module does not fail when the model is absent
# (the body plant is an optional extra)
def model_xml() -> str:
    return _find_model()


# ------------------------------------------------------------------ stretch reflex
# Protocol mirrors the standard experiment: a servo imposes a joint rotation while the
# reflex response is recorded. The rhythm generator is OFF (rg_gain=0) so this isolates
# the reflex arc from the locomotor rhythm.
DT, HOLD, RAMP, PRE, POST = 0.002, 0.30, 0.05, 0.40, 0.45     # seconds


def trajectory(q0: float, dq: float) -> list:
    """Pre-hold, fast ramp, hold, return. Returns a list of (t, q, v)."""
    out, t = [], 0.0

    def seg(T, f):
        """Sample a segment at s = DT, 2*DT, ... T, i.e. INCLUDING its endpoint.

        It previously sampled s = 0, DT, ... T-DT, so a ramp never reached its final
        value: with RAMP=0.05 and DT=0.002 the last ramp sample sat at 0.96*dq and the
        first HOLD sample jumped straight to dq. That is a 4% step in one 2 ms timestep --
        0.018 rad here -- which injects an impulse into joint velocity and therefore into
        the spindle's velocity term, which is the dominant input to the Ia afferent.
        """
        nonlocal t
        n = int(round(T / DT))
        for i in range(1, n + 1):
            q, v = f(i * DT)
            out.append((t, q, v))
            t += DT

    seg(PRE,  lambda s: (q0, 0.0))
    seg(RAMP, lambda s: (q0 + dq * s / RAMP, dq / RAMP))
    seg(HOLD, lambda s: (q0 + dq, 0.0))
    seg(POST, lambda s: (q0 + dq * (1 - s / POST), -dq / POST))
    return out


def stretch_reflex(drug: Drug | None = None, q0=0.1, dq=0.45, seed=1,
                   gaba_sens=1.0, gaba_sens_tonic=None, gaba_sens_phasic=None, glyr_sens=1.0,
                   ia_scale=None, xml=None) -> dict:
    """Run the ramp-and-hold stretch reflex and return phase-resolved measures.

    `ia_scale` rescales the Ia->Mn weight. It defaults to config.IA_SCALE because without
    it the motoneuron pool's dynamic peak pins at the tref=8 ms firing ceiling (125 Hz),
    which makes every drug read ~100% of control -- recurring error E5. Pass 1.0 only if
    you want to see that ceiling effect.
    """
    from .plant import JointPlant            # imports mujoco; optional extra
    from .circuit import SpinalCircuit
    from .config import IA_SCALE

    ia_scale = IA_SCALE if ia_scale is None else ia_scale
    # Passed as a per-instance override. This used to MUTATE SpinalCircuit.W (the class
    # attribute) and restore it in a `finally`, which corrupts every other instance in the
    # process if two assays run concurrently or an exception lands mid-way.
    w = {("Ia", "Mn", rec): SpinalCircuit.W[("Ia", "Mn", rec)] * ia_scale
         for rec in ("ampa", "nmda")}
    p = JointPlant(xml or model_xml(), drug=drug or Drug(), rg_gain=0.0, seed=seed,
                   gaba_sens=gaba_sens, gaba_sens_tonic=gaba_sens_tonic,
                   gaba_sens_phasic=gaba_sens_phasic, glyr_sens=glyr_sens, w=w)
    for (t, q, v) in trajectory(q0, dq):
        p.step(DT, impose=(q, v))
    A = p.arrays()

    tt = A["t"] / 1000.0
    base = (tt > 0.15) & (tt < PRE)                          # pre-stretch baseline
    dyn = (tt > PRE) & (tt < PRE + RAMP + 0.03)              # dynamic (ramp) phase
    sta = (tt > PRE + RAMP + 0.05) & (tt < PRE + RAMP + HOLD)  # static hold
    # the stretched muscle is the EXTENSOR (gear +1 lengthens as qpos rises)
    out = dict(
        ia_base=A["ia_E"][base].mean(), ia_dyn=A["ia_E"][dyn].max(),
        ia_sta=A["ia_E"][sta].mean(),
        mn_base=A["mn_E"][base].mean(), mn_dyn=A["mn_E"][dyn].max(),
        mn_sta=A["mn_E"][sta].mean(),
        mn_anta=A["mn_F"][sta].mean(),                        # reciprocal inhibition
        f_base=abs(A["f_ext"][base].mean()), f_dyn=abs(A["f_ext"][dyn]).max(),
        f_sta=abs(A["f_ext"][sta].mean()),
    )
    # REFLEX GAIN. The denominator is the Ia afferent response to the stretch. Clamping it
    # with max(1e-6, ...) was wrong in a way that inverts the result: if the Ia response is
    # NEGATIVE (noise, or severe depression dropping dynamic firing below baseline) the
    # clamp divides by 1e-6 and turns a small drop into a ~5,000,000 "gain". A reflex gain
    # is undefined when the afferent drive did not increase, so say so instead.
    ia_response = out["ia_dyn"] - out["ia_base"]
    if ia_response <= IA_RESPONSE_FLOOR:
        out["gain"] = float("nan")
        out["gain_undefined_reason"] = (
            f"Ia afferent response was {ia_response:+.3f} Hz, at or below the "
            f"{IA_RESPONSE_FLOOR} Hz floor: the stretch did not drive the afferent, so "
            "reflex gain has no denominator. Check the imposed trajectory and the spindle "
            "parameters rather than reading the ratio.")
    else:
        out["gain"] = (out["mn_dyn"] - out["mn_base"]) / ia_response
    return out


# reference conditions with known phenotypes, used by the reflex validation panel
REFLEX_CASES = [
    ("control",                       Drug()),
    ("GABA-A PAM 2x  (benzo-like)",   Drug(gaba_a_gain=2.0, gaba_a_tau=1.6)),
    ("GlyR PAM 1.6x  (ethanol-like)", Drug(glyr_gain=1.6)),
    ("GlyR block 0.4x (strychnine)",  Drug(glyr_gain=0.4)),
    ("NMDA 60% non-selective",        Drug(nmda_block=0.6, glun2b_selectivity=0.0)),
    ("NMDA 60% GluN2B-sel (spinal)",  Drug(nmda_block=0.6, glun2b_selectivity=1.0,
                                           glun2b_fraction=0.15)),
]


# ------------------------------------------------------------------- locomotion
# THE SECOND BEHAVIOURAL ENDPOINT, and the one the project lacked for a long time. Motor
# impairment was previously represented only by stretch-reflex gain, which is a thin proxy
# for what a sedative actually does to an animal's movement.
#
# This closes the chain receptor -> circuit -> muscle -> BODY MOVEMENT for locomotion:
# the group-pacemaker rhythm generator drives motoneuron pools, which drive Hill-type
# muscle actuators, which move the joint, whose spindles feed back. Nothing is imposed.
#
# Metrics are chosen to match what rodent gait analysis actually measures, so the model's
# output and an animal's output are the same kind of quantity:
#     step_period_ms  <-> cadence                      VALID as an impairment measure
#     alternation     <-> antagonist coordination      VALID
#     duty            <-> stance/swing ratio           VALID
#     excursion_rad   <-> stride length proxy          *** NOT VALID -- see below ***
#
# EXCURSION REPORTS THE WRONG SIGN AND MUST NOT BE USED AS AN IMPAIRMENT MEASURE.
# Measured here: control 1.607 rad, benzodiazepine-like PAM 2x -> 1.726, PAM 4x -> 1.940.
# A sedative INCREASES joint excursion in this preparation. The mechanism is real rather
# than a bug -- more inhibition means less co-contraction of the antagonist pair, so the
# joint is less stiff and swings further -- but it is an artefact of the preparation: a
# single joint with the rest of the body held fixed, no gravitational load and no ground
# contact. In an animal, reduced co-contraction shows up as instability and collapse, and
# this model has nothing to collapse against.
#
# A valid stride/ataxia measure needs the whole body and ground reaction forces. Until
# then use PERIOD and ALTERNATION, which degrade in the correct direction:
#     PAM 4x -> period 1250 -> 1667 ms (slowing) and alternation -0.69 -> -0.33
#     (coordination loss). Both are what a sedative does.
#
# A note on testing glycine here: `Drug(glyr_gain=0.4)` is a 60% block and leaves
# locomotion essentially unchanged (alternation -0.66 vs -0.69 control). That is recurring
# error E10, not a negative result -- this coupling is only informative under COMPLETE
# removal, which at circuit level means ie_gly = 0 (see tests/test_phenotypes.py).
# Minimum Ia afferent response (Hz) for a reflex gain to be defined. Below this the
# denominator is noise and the ratio is meaningless rather than large.
IA_RESPONSE_FLOOR = 1.0

LOCO_SETTLE_S = 1.0


def locomotion(drug: Drug | None = None, duration_s=6.0, seed=1, rg_gain=30.0,
               gaba_sens=1.0, gaba_sens_tonic=None, gaba_sens_phasic=None, glyr_sens=1.0,
               joint="knee_L", xml=None) -> dict:
    """Free-running closed-loop locomotion. No imposed kinematics.

    `rg_gain` converts the rhythm generator's PER-NEURON rate in Hz into an equivalent
    total presynaptic rate, so it is the presynaptic population size (30). Treating it as
    a free gain is recurring error E2.
    """
    from .plant import JointPlant

    p = JointPlant(xml or model_xml(), joint=joint, drug=drug or Drug(),
                   rg_gain=rg_gain, seed=seed, gaba_sens=gaba_sens,
                   gaba_sens_tonic=gaba_sens_tonic, gaba_sens_phasic=gaba_sens_phasic,
                   glyr_sens=glyr_sens)
    for _ in range(int(duration_s / DT)):
        p.step(DT)
    A = p.arrays()
    t = A["t"] / 1000.0
    m = t > LOCO_SETTLE_S
    q, mf, me = A["qpos"][m], A["mn_F"][m], A["mn_E"][m]

    out = dict(
        # NOT an impairment measure -- reports the wrong sign in this single-joint
        # preparation. Named to say so, because the plain name invites misuse.
        excursion_rad_INVALID=float(np.ptp(q)),
        mn_F_peak=float(mf.max()), mn_E_peak=float(me.max()),
        torque_range=float(np.ptp(A["torque"][m])))

    # frequency from the FFT of joint angle: amplitude-invariant, so a drug that weakens
    # the rhythm cannot masquerade as a faster one (recurring error E4)
    if q.std() > 1e-4:
        x = (q - q.mean()) * np.hanning(len(q))
        F = np.abs(np.fft.rfft(x))
        fr = np.fft.rfftfreq(len(q), DT)
        band = (fr > 0.2) & (fr < 8.0)
        f0 = float(fr[band][np.argmax(F[band])]) if band.any() else 0.0
        out["step_period_ms"] = 1000.0 / f0 if f0 > 0 else float("nan")
        out["freq_hz"] = f0
    else:
        out["step_period_ms"] = float("nan")
        out["freq_hz"] = 0.0

    out["alternation"] = (float(np.corrcoef(mf, me)[0, 1])
                          if mf.std() > 1e-9 and me.std() > 1e-9 else float("nan"))
    # duty: fraction of the cycle the extensor pool is above its mid-range
    thr = me.min() + 0.5 * np.ptp(me)
    out["duty"] = float((me > thr).mean()) if np.ptp(me) > 1e-9 else float("nan")
    # a drug can abolish locomotion outright; say so rather than returning a tiny number
    out["walking"] = bool(out["excursion_rad_INVALID"] > 0.1 and out["freq_hz"] > 0.2)
    return out
