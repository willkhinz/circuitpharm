"""Whether the Barberis deactivation trace is SAMPLED well enough to carry its own information.

WHY THIS SCRIPT EXISTS. `fitting.data.BARBERIS2007_DEACT` reproduces Barberis et al. 2007's
published triple exponential EXACTLY -- refitting the stored trace returns
tau = (2.800, 33.400, 221.350) ms and areas = (0.570, 0.230, 0.200) to the digit, and
sum(A_i tau_i) = 53.55 ms against the paper's 52.5 +/- 2.9 ms (0.36 SEM). By every check of
its VALUES the dataset is correct.

It was nevertheless built on a UNIFORM 10.1 ms grid over 0-1000 ms, and the components span
two decades (2.8 to 221 ms). The fast component carries the LARGEST area, 0.57, and a uniform
10.1 ms grid puts exactly ONE sample at or before one tau_fast. So the trace that looks like
it constrains the fast decay does not, and the dataset exists precisely to make an absolute
rate identifiable.

A noiseless fit hides this completely, because recovering three exponentials from an
exactly-determined noiseless inverse problem seeded at the truth succeeds on any grid. The
defect only appears at the project's own noise level. That is the whole point: a dataset can
reproduce its source perfectly and still be unable to do the job it was added for, and no
check of its values would notice.

Run:  python scripts/deact_sampling.py
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

# `python scripts/deact_sampling.py` puts scripts/ on sys.path, not the repo root, so the
# `scripts.` import below fails -- while `python -m scripts.deact_sampling` and pytest (via
# pyproject's pythonpath) both work. That discrepancy already cost this project a CI failure
# once (5643f74: tests/test_literature.py passed under `python -m pytest` and died under bare
# `pytest`), so it is handled here rather than left for the next person to rediscover.
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.literature_gamma2 import DEACTIVATIONS  # noqa: E402

_SRC = next(d for d in DEACTIVATIONS if d.label == "Barberis 2007, brief pulse")
#: The published components. Imported, never retyped -- see literature_gamma2 for provenance.
TAUS, AREAS = _SRC.components
TAU_W_PUBLISHED, TAU_W_SEM = _SRC.tau_w_ms, _SRC.err_ms
TAU_W_COMPONENTS = sum(a * t for a, t in zip(AREAS, TAUS))

#: The noise the project actually works at. knowledge/13-preregistration.md's P6 conclusion
#: is stated at sigma = 0.02 and asks whether sigma <= 0.005 is needed.
SIGMA_NOMINAL = 0.02
SIGMA_PREREG = 0.005

N_POINTS = 100
WINDOW_MS = 1000.0


def grids(n=N_POINTS, window=WINDOW_MS):
    """The two candidate time grids, at identical point counts."""
    return {
        "uniform": np.linspace(0.0, window, n),
        "log": np.r_[0.0, np.geomspace(0.1, window, n - 1)],
    }


def truth(t):
    return sum(a * np.exp(-t / tau) for a, tau in zip(AREAS, TAUS))


def _f3(t, a1, t1, a2, t2, a3, t3):
    return a1 * np.exp(-t / t1) + a2 * np.exp(-t / t2) + a3 * np.exp(-t / t3)


def recover(t, sigma, n_rep=400, seed=12345):
    """Refit the trace under noise, seeded AT THE TRUTH so this measures information
    content rather than optimiser luck. Returns per-replicate (tau_fast, A_fast, tau_w)."""
    y0 = truth(t)
    rng = np.random.default_rng(seed)
    p0 = [AREAS[0], TAUS[0], AREAS[1], TAUS[1], AREAS[2], TAUS[2]]
    out = []
    with np.errstate(over="ignore", invalid="ignore"):
        for _ in range(n_rep):
            yn = y0 + rng.normal(0.0, sigma, t.size)
            try:
                p, _ = curve_fit(_f3, t, yn, p0=p0, maxfev=100000)
            except Exception:
                continue
            if not np.all(np.isfinite(p)):
                continue
            out.append((p[1], p[0], p[0] * p[1] + p[2] * p[3] + p[4] * p[5]))
    return np.asarray(out)


def _robust(v):
    """Median and a normal-equivalent sigma from the IQR. Robust because the uniform grid
    produces runaway slow components, and a mean +/- SD over those says more about the
    blow-ups than about the typical fit."""
    med = float(np.median(v))
    q1, q3 = np.percentile(v, [25, 75])
    return med, float((q3 - q1) / 1.349)


def _within(v, true, rel=0.20):
    return 100.0 * float(np.mean(np.abs(np.asarray(v) - true) <= rel * abs(true)))


def sampling(n_rep=400):
    print("=" * 88)
    print("IS THE TRACE SAMPLED WELL ENOUGH TO CONSTRAIN ITS OWN COMPONENTS?")
    print("=" * 88)
    print(f"\npublished components  tau = {TAUS} ms, areas = {AREAS}")
    print(f"published tau_w       {TAU_W_PUBLISHED} +/- {TAU_W_SEM} ms (SEM, n = {_SRC.n})")
    print(f"sum(A_i tau_i)        {TAU_W_COMPONENTS:.2f} ms  "
          f"({abs(TAU_W_COMPONENTS - TAU_W_PUBLISHED) / TAU_W_SEM:.2f} SEM from published)")
    print(f"\nthe slow component carries {AREAS[2] * TAUS[2]:.1f} of {TAU_W_COMPONENTS:.1f} ms "
          f"({100 * AREAS[2] * TAUS[2] / TAU_W_COMPONENTS:.0f}%), so the {WINDOW_MS:.0f} ms "
          f"window is required;")
    print(f"the fast component carries the LARGEST area ({AREAS[0]}), so the EARLY sampling "
          f"is required too.")

    for name, t in grids().items():
        d = np.diff(t)
        n_fast = int((t <= TAUS[0]).sum())
        print(f"\n  {name:<8} {t.size} points, spacing {d.min():.3g}-{d.max():.3g} ms, "
              f"{n_fast} sample(s) at or before tau_fast = {TAUS[0]} ms")

    for sigma in (0.0, SIGMA_PREREG, SIGMA_NOMINAL):
        print(f"\n{'-' * 88}")
        print(f"noise sd = {sigma}" + ("  (noiseless -- hides the defect entirely)" if sigma == 0
              else f"  ({'P6 pre-registration target' if sigma == SIGMA_PREREG else 'the project nominal'})"))
        print(f"{'-' * 88}")
        print(f"  {'grid':<9} {'tau_fast (ms)':>20} {'A_fast':>16} {'tau_w (ms)':>22} "
              f"{'tau_w within 20%':>18}")
        for name, t in grids().items():
            r = recover(t, sigma, n_rep=1 if sigma == 0 else n_rep)
            tf, sf = _robust(r[:, 0])
            af, sa = _robust(r[:, 1])
            tw, sw = _robust(r[:, 2])
            print(f"  {name:<9} {f'{tf:7.2f} +/- {sf:<6.2f}':>20} "
                  f"{f'{af:5.3f} +/- {sa:<5.3f}':>16} {f'{tw:8.2f} +/- {sw:<8.2f}':>22} "
                  f"{_within(r[:, 2], TAU_W_COMPONENTS):>17.0f}%")
        print(f"  {'TRUTH':<9} {TAUS[0]:>13.2f} {AREAS[0]:>19.3f} {TAU_W_COMPONENTS:>22.2f}")

    print(f"\n{'=' * 88}")
    print("VERDICT")
    print("=" * 88)
    r_u = recover(grids()["uniform"], SIGMA_NOMINAL, n_rep=n_rep)
    r_l = recover(grids()["log"], SIGMA_NOMINAL, n_rep=n_rep)
    tfu, sfu = _robust(r_u[:, 0])
    tfl, sfl = _robust(r_l[:, 0])
    print(f"  At the project's own noise (sd = {SIGMA_NOMINAL}), with the SAME {N_POINTS} points:")
    print(f"    uniform  tau_fast = {tfu:.2f} +/- {sfu:.2f} ms  "
          f"({100 * sfu / tfu:.0f}% scatter), tau_w within 20%: {_within(r_u[:, 2], TAU_W_COMPONENTS):.0f}%")
    print(f"    log      tau_fast = {tfl:.2f} +/- {sfl:.2f} ms  "
          f"({100 * sfl / tfl:.0f}% scatter), tau_w within 20%: {_within(r_l[:, 2], TAU_W_COMPONENTS):.0f}%")
    print(f"\n  A uniform grid over a two-decade decay spends its points where nothing is")
    print(f"  happening and takes one sample where 57% of the amplitude lives. Log spacing")
    print(f"  costs nothing -- same point count, same window, same source values -- and is")
    print(f"  what the dataset uses.")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-rep", type=int, default=400)
    sampling(ap.parse_args().n_rep)


if __name__ == "__main__":
    main()
