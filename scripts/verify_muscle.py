"""Verify the muscle-ported rodent: compiles, muscles are real, antagonists oppose."""
import numpy as np, mujoco, time, pathlib
# SIDE-EFFECT GUARD added 2026-10-07. This script previously did its work at MODULE level,
# so merely importing it ran it. For `port_muscle.py` that regenerated the MuJoCo body
# model, and for `build_kb.py`/`build_compounds.py` it rewrote the pharmacology database --
# both destructive, both triggered by any tool that imports or scans the package. The body
# is unchanged; it is now reached only when the file is run directly.


def main():
    M = pathlib.Path(__file__).resolve().parent.parent / "models" / "rodent_muscle.xml"
    m = mujoco.MjModel.from_xml_path(str(M)); d = mujoco.MjData(m)
    print(f"compiled OK: nq={m.nq} nv={m.nv} nu={m.nu} mass={m.body_subtreemass[0]:.4f}kg")

    gt = {}
    for i in range(m.nu): gt[int(m.actuator_gaintype[i])] = gt.get(int(m.actuator_gaintype[i]),0)+1
    print("gaintype histogram:", gt, " (2 = muscle)")
    nm = int((m.actuator_gaintype == mujoco.mjtGain.mjGAIN_MUSCLE).sum())
    print(f"muscle actuators: {nm}   dyntype=muscle: {int((m.actuator_dyntype==mujoco.mjtDyn.mjDYN_MUSCLE).sum())}")

    ids = {mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i): i for i in range(m.nu)}
    print("\nknee_L pair lengthrange / force:")
    for n in ("knee_L_ext","knee_L_flx"):
        i = ids[n]
        print(f"  {n:14s} gear={m.actuator_gear[i][0]:+.1f} lengthrange={m.actuator_lengthrange[i]}"
              f" force={m.actuator_gainprm[i][2]:.4g}")

    # antagonists must produce OPPOSING torque on the same joint
    jid = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "knee_L")
    dofadr = m.jnt_dofadr[jid]
    for n in ("knee_L_ext","knee_L_flx"):
        d2 = mujoco.MjData(m)
        mujoco.mj_forward(m, d2)
        d2.act[:] = 0; d2.ctrl[:] = 0
        d2.ctrl[ids[n]] = 1.0
        d2.act[m.actuator_actadr[ids[n]]] = 1.0      # fully activated
        mujoco.mj_forward(m, d2)
        print(f"  {n:14s} -> actuator_force={d2.actuator_force[ids[n]]:+.5f}"
              f"  qfrc_actuator[knee_L]={d2.qfrc_actuator[dofadr]:+.5f}")

    t0=time.time(); N=2000
    for _ in range(N):
        d.ctrl[:] = 0.0; mujoco.mj_step(m, d)
    el=time.time()-t0
    print(f"\nphysics: {N/el:,.0f} steps/s ({N/el*m.opt.timestep:.0f}x realtime)")
    print("OK")


if __name__ == "__main__":
    main()
