"""Closed-loop coupling: MuJoCo muscle plant <-> spiking spinal circuit.

    Mn_F / Mn_E firing rate  --normalised-->  muscle activation (d.ctrl)
    actuator length, velocity --spindle-->    Ia afferent rate --> circuit

The CPG integrates at 0.1 ms and MuJoCo at 2 ms, so the circuit takes 20 substeps per
physics step. Muscle length is normalised against the actuator's own lengthrange, so the
spindle sees a dimensionless fascicle length in [0,1] with velocity in units of 1/s.

Joint convention: for a joint-transmission muscle, actuator_length = gear * qpos and
actuator_velocity = gear * qvel, so '_ext' (gear +1) lengthens as qpos rises and '_flx'
(gear -1) shortens. Each muscle's own spindle therefore reports ITS fascicle length,
which is what reciprocal Ia inhibition needs.
"""
import numpy as np
import mujoco
from .cpg import spindle_ia
from .circuit import SpinalCircuit


class JointPlant:
    """One antagonist joint pair, closed loop, with the rest of the body held fixed."""

    def __init__(self, xml, joint="knee_L", drug=None, rg_gain=900.0,
                 mn_sat=80.0, ia_kl=55.0, ia_kv=70.0, ia_bias=12.0, seed=0, gaba_sens=1.0, **kw):
        self.m = mujoco.MjModel.from_xml_path(xml)
        self.d = mujoco.MjData(self.m)
        self.joint = joint
        self.jid = mujoco.mj_name2id(self.m, mujoco.mjtObj.mjOBJ_JOINT, joint)
        self.qa = self.m.jnt_qposadr[self.jid]
        self.va = self.m.jnt_dofadr[self.jid]
        self.lo, self.hi = self.m.jnt_range[self.jid]
        aid = lambda n: mujoco.mj_name2id(self.m, mujoco.mjtObj.mjOBJ_ACTUATOR, n)
        self.a_ext, self.a_flx = aid(f"{joint}_ext"), aid(f"{joint}_flx")
        self.LR = {h: self.m.actuator_lengthrange[a].copy()
                   for h, a in (("E", self.a_ext), ("F", self.a_flx))}
        self.circuit = SpinalCircuit(drug=drug, rg_gain=rg_gain, seed=seed,
                                    gaba_sens=gaba_sens, **kw)
        self.mn_sat = mn_sat          # Hz of Mn pool rate mapped to activation 1.0
        self.ia = dict(kl=ia_kl, kv=ia_kv, bias=ia_bias)
        self.log = {k: [] for k in
                    ("t", "qpos", "qvel", "ia_F", "ia_E", "mn_F", "mn_E",
                     "act_F", "act_E", "f_ext", "f_flx", "torque")}

    # ---- spindle: normalise this muscle's fascicle length/velocity ----
    def _spindle(self, half):
        a = self.a_flx if half == "F" else self.a_ext
        L0, L1 = self.LR[half]
        span = max(1e-9, L1 - L0)
        Ln = (self.d.actuator_length[a] - L0) / span
        Vn = self.d.actuator_velocity[a] / span          # 1/s
        return spindle_ia(Ln, Vn, l0=0.5, **self.ia)

    def step(self, dt_phys=0.002, n_sub=20, impose=None):
        """impose: optional (qpos, qvel) to kinematically drive the joint, as a servo
        motor does in a stretch-reflex experiment."""
        if impose is not None:
            q, v = impose
            self.d.qpos[self.qa] = q
            self.d.qvel[self.va] = v
        mujoco.mj_forward(self.m, self.d)
        ia_F, ia_E = self._spindle("F"), self._spindle("E")
        dt_c = dt_phys * 1000.0 / n_sub                  # ms
        for _ in range(n_sub):
            self.circuit.step(dt_c, ia_F=ia_F, ia_E=ia_E)
        mnF = self.circuit.pops["Mn_F"].rate
        mnE = self.circuit.pops["Mn_E"].rate
        aF = float(np.clip(mnF / self.mn_sat, 0.0, 1.0))
        aE = float(np.clip(mnE / self.mn_sat, 0.0, 1.0))
        self.d.ctrl[self.a_flx] = aF
        self.d.ctrl[self.a_ext] = aE
        if impose is None:
            mujoco.mj_step(self.m, self.d)
        else:
            mujoco.mj_forward(self.m, self.d)
            self.circuit.t  # keep clock
        fe = float(self.d.actuator_force[self.a_ext])
        ff = float(self.d.actuator_force[self.a_flx])
        for k, v in (("t", self.circuit.t), ("qpos", float(self.d.qpos[self.qa])),
                     ("qvel", float(self.d.qvel[self.va])), ("ia_F", ia_F), ("ia_E", ia_E),
                     ("mn_F", mnF), ("mn_E", mnE), ("act_F", aF), ("act_E", aE),
                     ("f_ext", fe), ("f_flx", ff),
                     ("torque", fe * 1.0 + ff * -1.0)):
            self.log[k].append(v)
        return mnF, mnE

    def arrays(self):
        return {k: np.asarray(v) for k, v in self.log.items()}
