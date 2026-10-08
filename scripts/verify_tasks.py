import numpy as np, time
from dm_control.locomotion.examples import basic_rodent_2020 as br
# SIDE-EFFECT GUARD added 2026-10-07. This script previously did its work at MODULE level,
# so merely importing it ran it. For `port_muscle.py` that regenerated the MuJoCo body
# model, and for `build_kb.py`/`build_compounds.py` it rewrote the pharmacology database --
# both destructive, both triggered by any tool that imports or scans the package. The body
# is unchanged; it is now reached only when the file is run directly.


def main():

    for name in ["rodent_escape_bowl", "rodent_run_gaps", "rodent_maze_forage", "rodent_two_touch"]:
        try:
            env = getattr(br, name)()
            ts = env.reset()
            spec = env.action_spec()
            obs = ts.observation
            t0 = time.time(); n = 200
            for _ in range(n):
                env.step(np.zeros(spec.shape))
            dt = time.time() - t0
            print(f"{name:22s} act={spec.shape[0]:3d}  n_obs={len(obs):2d}  {n/dt:7.0f} env-steps/s")
            if name == "rodent_escape_bowl":
                print("   observation keys:", sorted(obs.keys()))
        except Exception as e:
            print(f"{name:22s} FAILED: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
