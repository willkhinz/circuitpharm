"""Substrate selection, shared by every circuit.

WHY THIS IS SHARED RATHER THAN COPIED. `resp.py` grew this logic first: validate a
conductance operating point, build the right population class, hand NMDA over unevaluated.
`rg2.py` and `circuit.py` need exactly the same thing, and copying it into them is how
recurring error E12 happened -- a fix applied everywhere except the one module written later,
where the same defect reappeared and was found only by a review pass two sessions on.

WHAT THE MIGRATION ESTABLISHED, and what this module enforces. A single multiplicative weight
scale cannot carry a weight table between the LIF cell (C=200 pF, g_L=10 nS) and the Butera
cell (C=21 pF, g_L=2.8 nS). Six anchoring attempts on the preBotC found four quantities each
needing a DIFFERENT transformation:

    AMPA / GABA / glycine   x GL_RATIO                  = 0.28
    NMDA                    x GL_RATIO x RELIEF_RATIO   = 0.025      (11x smaller)
    excitatory drive (pA)   from the cell's measured bursting window, + tonic compensation
    interneuron bias (pA)   SMALL, and bounded ABOVE by depolarisation block

The NMDA factor is the ratio of the LIF's measured mean Mg2+ relief (0.063,
scripts/diag_substrate_limits.py) to a depolarised conductance cell's (~0.70). The LIF's
weights were tuned on a substrate where NMDA delivered ~6% of its nominal conductance and
could not relieve; on a cell that does relieve, the same weight is ~11x too strong -- and
because relief RISES with depolarisation it is positive feedback into depolarisation block,
not a scale error. That one mistake produced three consecutive sweeps in which nothing
oscillated at any of 48, 72 and 81 operating points.

`derive_weights` therefore exists so no caller has to rediscover it.
"""
from __future__ import annotations

import numpy as np

#: Input-conductance ratio between the Butera cell (2.8 nS) and the LIF (10 nS). Every
#: conductance is "some fraction of the cell's leak", so this is how one carries across.
GL_RATIO = 2.8 / 10.0

#: The LIF's measured mean Mg2+ relief (0.063) over a depolarised conductance cell's (~0.70).
#: Applies to NMDA weights ONLY, on top of GL_RATIO.
RELIEF_RATIO = 0.063 / 0.70

#: Slowest timescale on the Butera cell: tau_h = 10 s, 25x the LIF's tau_adapt of 400 ms.
#: Anything measured with less than a few multiples of this is measuring a transient.
TAU_H_MS = 10000.0
MIN_SETTLE_MS = 3.0 * TAU_H_MS


def derive_weights(lif_w: dict, s_ampa=GL_RATIO, s_nmda=GL_RATIO * RELIEF_RATIO,
                   s_inh=GL_RATIO) -> dict:
    """Carry a LIF weight table onto the conductance cell, by receptor type.

    Defaults are the derived values. They are separated by receptor because the NMDA
    correction is ~11x the others and nothing in a single scale can express that.
    """
    out = {}
    for k, v in lif_w.items():
        if k.endswith("_nmda"):
            out[k] = v * s_nmda
        elif k.endswith("_ampa"):
            out[k] = v * s_ampa
        else:                                   # gaba, gly
            out[k] = v * s_inh
    return out


def resolve_op(substrate: str, op, registered, required_w, required_scalars, what: str):
    """Validate a conductance operating point, or refuse.

    Returns `None` for the LIF substrate (nothing to resolve), else the validated op.

    THE PARTIAL OPERATING POINT IS THE DANGEROUS CASE, not the missing one. Supply `drive`
    and the excitatory weights while the inhibitory ones quietly keep their LIF values and
    the network is internally inconsistent by ~an order of magnitude on the arm the drug acts
    through. It runs, produces a rhythm, and is wrong about pharmacology. So the full set is
    required rather than defaulted.
    """
    if substrate == "lif":
        return None
    if substrate != "cond":
        raise ValueError(f"substrate must be 'lif' or 'cond', got {substrate!r}")
    if op is None:
        op = registered
    if op is None:
        raise ValueError(
            f"substrate='cond' needs its own operating point for {what}. The LIF's is in pA "
            f"and nS relative to a C=200 pF / g_L=10 nS cell; the Butera cell is "
            f"C=21 pF / g_L=2.8 nS, so reusing it is an order-of-magnitude error that fails "
            f"silently. Use substrate.derive_weights() and anchor the drives, or pass "
            f"op={{...}} explicitly.")
    missing_w = [k for k in required_w if k not in (op.get("w") or {})]
    missing_s = [k for k in required_scalars if k not in op]
    if missing_w or missing_s:
        raise ValueError(
            f"conductance operating point for {what} is incomplete. missing weights: "
            f"{missing_w or 'none'}; missing scalars: {missing_s or 'none'}. Every pA and nS "
            f"quantity is relative to the cell, so a partial operating point leaves part of "
            f"the network LIF-scaled -- including the inhibitory arm the drug acts through.")
    return op


def make_pop(substrate: str, n: int, name: str, *, nap=False, g_adapt,
             tau_adapt=None, tref=5.0, cell=None):
    """Build the population class this substrate calls for.

    ADAPTATION IS NOT CARRIED ACROSS, and this is the whole reason population construction
    had to be routed through a factory rather than parameterised in place. Every circuit
    assigned `g_adapt` AFTER constructing a `Pop`, so pointing one at a conductance cell
    would switch spike-triggered adaptation back on at its LIF-tuned strength on top of a
    now-real I_NaP inactivation. Burst termination would be double-counted and the symptom
    would be bursts ending slightly early -- a plausible duty cycle, which is this project's
    signature failure shape. Predicted as E18 before the cell was written.

    `g_adapt` IS REQUIRED, AND `None` MEANS "KEEP THE CELL CLASS'S OWN DEFAULT". Both halves
    of that are load-bearing, and the reason is a regression this function caused. `resp.py`
    and `rg2.py` both set `g_adapt` explicitly before the migration, so routing them through
    a factory with a default was harmless. `circuit.py` did NOT -- it constructed bare
    `Pop`s and set only `tref`, inheriting the dataclass default of 0.55 nS -- so passing
    this function's own default of 0.0 deleted spike-triggered adaptation from every spinal
    population. Two phenotype tests caught it (the stretch reflex and the step cycle);
    nothing else did, and the LIF path is supposed to be byte-identical across this
    migration. Writing `g_adapt=0.55` here instead would have duplicated `Pop`'s default
    into a second file, which is the divergence this module exists to prevent -- so the
    default is deleted and callers that want the class's value ask for it by name.

    `tau_adapt=None` means the same thing. It happened to equal `Pop`'s default, so it hid
    behind the same mistake without contributing to it.
    """
    if substrate == "lif":
        from .cpg import Pop
        p = Pop(n=n, name=name, nap=nap)
        if g_adapt is not None:
            p.g_adapt = g_adapt
        if tau_adapt is not None:
            p.tau_adapt = tau_adapt
        p.tref = tref
        return p
    if substrate == "cond":
        from .neuron import BRS1999_MODEL1, CondPop
        return CondPop(n=n, name=name, params=cell or BRS1999_MODEL1, nap=nap,
                       g_adapt=0.0)
    raise ValueError(f"substrate must be 'lif' or 'cond', got {substrate!r}")


def raw_receptors(substrate: str) -> tuple:
    """Receptors whose conductance must be handed over UNEVALUATED on this substrate.

    `Syn.conductance(V)` applies the Mg2+ block at the PRE-STEP voltage. Harmless on the
    LIF, where V moves a few mV per step; on a cell that sub-steps through a ~70 mV spike it
    discards the entire voltage-dependent relief. The symptom would be "the substrate change
    did not move the NMDA result" -- a reassuring robustness check that is in fact the bug.
    Predicted as E15.
    """
    return ("nmda",) if substrate == "cond" else ()


def tonic_gaba(gaba_tonic: float, sens_tonic: float, scale_tonic: float, g_syn):
    """Standing extrasynaptic GABA-A, added to the synaptic conductance and clamped at zero.

    A NAM (scale < 1) meeting a sensitivity above 1 makes the multiplier negative, and a
    negative conductance inverts an inhibitory shunt into regenerative negative damping --
    the voltage diverges instead of the run failing visibly. Same class as the clamp in
    `Drug.nmda_scale`, which was reachable the same way.
    """
    eff = max(0.0, 1.0 + sens_tonic * (scale_tonic - 1.0))
    return np.maximum(0.0, g_syn + gaba_tonic * eff)
