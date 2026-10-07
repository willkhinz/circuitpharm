"""Stretch-reflex panel. The protocol itself now lives in `circuitpharm.assays`.

This file was a library masquerading as a script: nine other scripts did
`import scripts.reflex as R` to get `XML`, `trajectory` and `run`. That made a `scripts/`
module a dependency of the package's own analyses, and its `XML = "models/..."` relative
path broke anything run from outside the repository root.

The protocol moved to `circuitpharm.assays`, which resolves the model path properly. The
names below are re-exported so existing callers keep working; new code should import from
the package.
"""
from circuitpharm.assays import (                      # noqa: F401
    trajectory, stretch_reflex, model_xml,
    REFLEX_CASES as CASES,
    DT, HOLD, RAMP, PRE, POST,
)

# backward-compatible aliases
run = stretch_reflex


def XML():                                             # noqa: N802
    """Deprecated: `XML` was a module-level relative path string. Call `model_xml()`."""
    return model_xml()


if __name__ == "__main__":
    print(f"{'condition':<32} {'Ia base':>8} {'Ia dyn':>7} {'Mn base':>8} "
          f"{'Mn dyn':>7} {'Mn sta':>7} {'Mn ant':>7} {'F dyn':>7} {'refl gain':>10}")
    print("-" * 104)
    ctrl = None
    for lab, d in CASES:
        r = stretch_reflex(d)
        gain = r["gain"]
        if ctrl is None:
            ctrl = gain
        print(f"{lab:<32} {r['ia_base']:8.1f} {r['ia_dyn']:7.1f} {r['mn_base']:8.1f} "
              f"{r['mn_dyn']:7.1f} {r['mn_sta']:7.1f} {r['mn_anta']:7.1f} "
              f"{r['f_dyn']:7.3f} {gain:7.3f} ({100*gain/ctrl:3.0f}%)")
