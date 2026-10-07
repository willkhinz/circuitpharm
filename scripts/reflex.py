"""Monosynaptic stretch reflex: ramp-and-hold, imposed kinematics.

Protocol mirrors the standard experiment -- a servo imposes a joint rotation while the
reflex response is recorded. RG drive is off (rg_gain=0) so this isolates the reflex arc
from the locomotor rhythm.
"""
import sys; sys.path.insert(0,".")
import numpy as np
from circuitpharm.plant import JointPlant
from circuitpharm.cpg import Drug

XML = "models/rodent_muscle.xml"
DT, HOLD, RAMP, PRE, POST = 0.002, 0.30, 0.05, 0.40, 0.45   # s

def trajectory(q0, dq):
    """pre-hold, fast ramp, hold, return. Returns list of (t, q, v)."""
    out, t = [], 0.0
    def seg(T, f):
        nonlocal t
        n = int(T/DT)
        for i in range(n):
            tau = i*DT
            q, v = f(tau)
            out.append((t, q, v)); t += DT
    seg(PRE,  lambda s: (q0, 0.0))
    seg(RAMP, lambda s: (q0 + dq*s/RAMP, dq/RAMP))
    seg(HOLD, lambda s: (q0 + dq, 0.0))
    seg(POST, lambda s: (q0 + dq*(1 - s/POST), -dq/POST))
    return out

def run(drug, q0=0.1, dq=0.45, seed=1):
    p = JointPlant(XML, drug=drug, rg_gain=0.0, seed=seed)
    traj = trajectory(q0, dq)
    for (t, q, v) in traj:
        p.step(DT, impose=(q, v))
    A = p.arrays()
    tt = A["t"]/1000.0
    base = (tt > 0.15) & (tt < PRE)                 # pre-stretch baseline
    dyn  = (tt > PRE) & (tt < PRE+RAMP+0.03)        # dynamic (ramp) phase
    sta  = (tt > PRE+RAMP+0.05) & (tt < PRE+RAMP+HOLD)   # static hold
    # the stretched muscle here is the EXTENSOR (gear +1 lengthens as qpos rises)
    return dict(
        ia_base=A["ia_E"][base].mean(), ia_dyn=A["ia_E"][dyn].max(),
        ia_sta=A["ia_E"][sta].mean(),
        mn_base=A["mn_E"][base].mean(), mn_dyn=A["mn_E"][dyn].max(),
        mn_sta=A["mn_E"][sta].mean(),
        mn_anta=A["mn_F"][sta].mean(),               # reciprocal inhibition target
        f_base=abs(A["f_ext"][base].mean()), f_dyn=abs(A["f_ext"][dyn]).max(),
        f_sta=abs(A["f_ext"][sta].mean()),
    )

CASES = [
    ("control",                       Drug()),
    ("GABA-A PAM 2x  (benzo-like)",   Drug(gaba_a_gain=2.0, gaba_a_tau=1.6)),
    ("GlyR PAM 1.6x  (ethanol-like)", Drug(glyr_gain=1.6)),
    ("GlyR block 0.4x (strychnine)",  Drug(glyr_gain=0.4)),
    ("NMDA 60% non-selective",        Drug(nmda_block=0.6, glun2b_selectivity=0.0)),
    ("NMDA 60% GluN2B-sel (spinal)",  Drug(nmda_block=0.6, glun2b_selectivity=1.0,
                                           glun2b_fraction=0.15)),
]

if __name__ == "__main__":
    print(f"{'condition':<32} {'Ia base':>8} {'Ia dyn':>7} {'Mn base':>8} {'Mn dyn':>7} "
          f"{'Mn sta':>7} {'Mn ant':>7} {'F dyn':>7} {'refl gain':>10}")
    print("-"*104)
    ctrl = None
    for lab, d in CASES:
        r = run(d)
        gain = (r['mn_dyn'] - r['mn_base']) / max(1e-6, r['ia_dyn'] - r['ia_base'])
        if ctrl is None: ctrl = gain
        print(f"{lab:<32} {r['ia_base']:8.1f} {r['ia_dyn']:7.1f} {r['mn_base']:8.1f} "
              f"{r['mn_dyn']:7.1f} {r['mn_sta']:7.1f} {r['mn_anta']:7.1f} "
              f"{r['f_dyn']:7.3f} {gain:7.3f} ({100*gain/ctrl:3.0f}%)")
