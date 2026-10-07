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
        nonlocal t
        for i in range(int(T / DT)):
            q, v = f(i * DT)
            out.append((t, q, v))
            t += DT

    seg(PRE,  lambda s: (q0, 0.0))
    seg(RAMP, lambda s: (q0 + dq * s / RAMP, dq / RAMP))
    seg(HOLD, lambda s: (q0 + dq, 0.0))
    seg(POST, lambda s: (q0 + dq * (1 - s / POST), -dq / POST))
    return out


def stretch_reflex(drug: Drug | None = None, q0=0.1, dq=0.45, seed=1,
                   gaba_sens=1.0, gaba_sens_tonic=None, gaba_sens_phasic=None,
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
    orig = dict(SpinalCircuit.W)
    try:
        SpinalCircuit.W = dict(orig)
        for rec in ("ampa", "nmda"):
            SpinalCircuit.W[("Ia", "Mn", rec)] = orig[("Ia", "Mn", rec)] * ia_scale
        p = JointPlant(xml or model_xml(), drug=drug or Drug(), rg_gain=0.0, seed=seed,
                       gaba_sens=gaba_sens, gaba_sens_tonic=gaba_sens_tonic,
                       gaba_sens_phasic=gaba_sens_phasic)
        for (t, q, v) in trajectory(q0, dq):
            p.step(DT, impose=(q, v))
        A = p.arrays()
    finally:
        SpinalCircuit.W = orig

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
    out["gain"] = ((out["mn_dyn"] - out["mn_base"])
                   / max(1e-6, out["ia_dyn"] - out["ia_base"]))
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
