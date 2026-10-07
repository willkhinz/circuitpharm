import numpy as np, time
from dm_control.locomotion.examples import basic_rodent_2020 as br

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
