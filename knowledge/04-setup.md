# Environment

macOS arm64 (Apple Silicon), 32 GB RAM. System Python is 3.13; the project venv pins
**3.11** because dm_control/mujoco wheels are reliable there.

```bash
cd /Users/whinz/Biochemistry
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python mujoco dm_control numpy
```

Installed: mujoco 3.15.0, dm_control, numpy 2.4.6, scipy 1.17.1.

## Gotchas

- **Do not set `MUJOCO_GL=osmesa` on macOS** — not a valid backend here; it raises at
  import. Leave `MUJOCO_GL` unset (defaults to glfw). Physics runs headless regardless.
- `dm_control.locomotion.walkers.assets` is a namespace package with no `__file__`.
  Resolve the model path from `dm_control.__file__` instead:
  `os.path.join(os.path.dirname(dm_control.__file__), "locomotion/walkers/assets/rodent.xml")`
- The egocentric camera dominates env-step cost (~41,700 physics steps/s bare vs.
  160-310 env-steps/s with vision). Disable vision for motor assays.

## Scripts

| Script | Purpose |
|---|---|
| `scripts/verify_rodent.py` | Loads rodent.xml, prints model dims, checks for muscle actuators, times physics |
| `scripts/verify_tasks.py` | Instantiates and steps the four rodent benchmark tasks |
| `scripts/build_kb.py` | Rebuilds `data/pharmacology.db` from inline source data (idempotent) |
| `scripts/kb.py` | Query helper for the knowledge DB |

## Upstream references

- `rodent.xml` ships with dm_control; it is the Aldarondo et al. 2024 body.
- `github.com/diegoaldarondo/virtual_rodent` — STAC skeletal registration, motion-mapper
  behavioural classification, inverse-dynamics inference.
- `github.com/emiwar/julia-virtual-rodent` — PPO reimplementation, same body and data.
- Trained policies and mocap data are NOT in dm_control and are large; not downloaded.
