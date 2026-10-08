"""Verify the dm_control rodent (Aldarondo et al. 2024 body) loads and steps."""
import os, time
# macOS: leave MUJOCO_GL unset (defaults to glfw); physics runs headless regardless
import numpy as np
import mujoco
from dm_control import mjcf
import dm_control
# SIDE-EFFECT GUARD added 2026-10-07. This script previously did its work at MODULE level,
# so merely importing it ran it. For `port_muscle.py` that regenerated the MuJoCo body
# model, and for `build_kb.py`/`build_compounds.py` it rewrote the pharmacology database --
# both destructive, both triggered by any tool that imports or scans the package. The body
# is unchanged; it is now reached only when the file is run directly.


def main():

    xml = os.path.join(os.path.dirname(dm_control.__file__), "locomotion", "walkers", "assets", "rodent.xml")
    print("model xml:", xml, "exists:", os.path.exists(xml))

    m = mujoco.MjModel.from_xml_path(xml)
    d = mujoco.MjData(m)
    print(f"nq(joint coords)={m.nq}  nv(dof)={m.nv}  nu(actuators)={m.nu}  nbody={m.nbody}")
    print(f"total mass = {m.body_subtreemass[0]:.4f} kg")
    print(f"timestep = {m.opt.timestep} s")

    # actuator types: 0=motor(torque) 1=position 2=velocity 3=muscle ...
    types = {}
    for i in range(m.nu):
        g = m.actuator_gaintype[i]
        types[int(g)] = types.get(int(g), 0) + 1
    print("actuator gaintype histogram:", types, " (0=fixed/torque, 2=muscle)")
    print("has muscle actuators:", bool((m.actuator_gaintype == 2).any()))

    # step physics
    t0 = time.time()
    N = 2000
    for _ in range(N):
        d.ctrl[:] = 0.0
        mujoco.mj_step(m, d)
    el = time.time() - t0
    print(f"stepped {N} steps in {el:.2f}s -> {N/el:,.0f} steps/s  ({N*m.opt.timestep:.2f}s sim)")
    print("torso height after settle:", round(float(d.qpos[2]), 4))
    print("OK")


if __name__ == "__main__":
    main()
