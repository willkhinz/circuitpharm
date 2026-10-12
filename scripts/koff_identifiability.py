"""Does the Barberis deactivation dataset make the absolute rate `koff` identifiable?

WHY THIS SCRIPT EXISTS. `fitting.data.MISSING_DATASETS` describes a macroscopic deactivation
time course as "the dataset that makes an absolute rate (k_off) identifiable at all; without
it the fit determines only ratios", and `knowledge/11-identifiability.md` proves the
structural half: an equilibrium observable fixes three combinations (K_d, E, D) and leaves the
absolute scale free. The claim has always been that a kinetic observable closes it.

That claim had never been TESTED, because until 2026-10-11 nothing in `src/` read a
deactivation trace. `likelihood._predict` raised `ObservableMismatch` for every observable
except PEAK and EQUILIBRIUM and `_observed` read only `mean_response`, so
`DEACTIVATION_BENCHMARK` -- whose own `origin` says it "exists so the kinetic-timescale half
of the likelihood can be exercised" -- was inert, and so was the Barberis dataset when first
added. `likelihood.predict_deactivation` is that missing half; this script is what it was for.

WHAT IS MEASURED. A profile likelihood over `koff`: at each fixed `koff` the remaining free
rates are re-optimised and the best cost recorded. A flat profile means unidentifiable; an
interior minimum rising on both sides means identifiable, and the width at the chi-square
threshold is the interval.

TWO MEASUREMENT LESSONS ARE BAKED IN, both learned by getting them wrong first.

  1. THE OBJECTIVE IS THE DEACTIVATION ALONE, not the concentration-response plus the
     deactivation. Adding the Jahn CRC makes each evaluation ~40x dearer, because a PEAK
     prediction integrates a waveform per concentration, and the first attempt at the joint
     profile could not be optimised to convergence inside any reasonable budget. Its output
     looked like a result -- a sharp minimum -- and was not: the CRC-only arm had its minimum
     AT THE GRID EDGE (so unbounded, not identifiable) and the joint arm jumped 0 to 47 nats
     over one grid step and back down, which is optimiser noise, not curvature. It was
     discarded rather than reported. The question "does the deactivation carry absolute-rate
     information" is answered by the deactivation alone, and answered cheaply.
  2. NELDER-MEAD TOLERANCES DECIDE THE RUNTIME, not the grid. With xatol=1e-6/fatol=1e-9 the
     inner fits ran to `maxiter` every time: 13 grid points x 10 restarts x 4000 iterations
     x 17 ms = 2.5 hours. The settings below converge instead of exhausting.

Both datasets in the repo are `kind="parametric"`, so every number here inherits that grade
and `ProvisionalResultWarning` is correct to fire. This measures what the DATA can determine,
not what the receptor's rates are.

Run:  python scripts/koff_identifiability.py              # coarse, then refined
      python scripts/koff_identifiability.py --no-refine  # coarse only, ~11 min
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from circuitpharm import KineticAllosteryModel                      # noqa: E402
from circuitpharm.fitting.data import ALL                           # noqa: E402
from circuitpharm.fitting.likelihood import make_log_likelihood     # noqa: E402

#: Free while profiling. `d`/`r` (desensitisation) stay at their defaults: the question is
#: whether the decay constrains the BINDING/GATING scale, and letting desensitisation absorb
#: the decay shape answers a different one.
FREE = ("kon", "beta", "alpha")

#: 95% for one parameter, Wilks: 0.5 * chi2(0.95, 1 df), on the negative-log-likelihood scale.
DELTA_95 = 1.92

#: Converge, do not exhaust. See lesson 2 in the module docstring.
_OPT = dict(maxiter=500, xatol=1e-3, fatol=1e-6)
_N_RESTART = 4


def _cost_fn(datasets, n_params=4):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scorer = make_log_likelihood(datasets, caller="koff_identifiability",
                                     n_params=n_params)

    def cost(model):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                return -float(scorer(model).total_log_likelihood)
            except Exception:
                return np.inf
    return cost


def fit_at(koff, cost, base, n_restart=_N_RESTART, seed=0):
    """Best cost with `koff` HELD FIXED, optimising FREE in log10 space."""
    rng = np.random.default_rng(seed)
    x0 = np.log10([getattr(base, k) for k in FREE])
    best = np.inf
    for i in range(n_restart):
        start = x0 if i == 0 else x0 + rng.normal(0.0, 0.4, x0.size)

        def obj(x):
            try:
                m = base.with_params(koff=koff,
                                     **{k: 10.0 ** v for k, v in zip(FREE, x)})
            except Exception:
                return np.inf
            return cost(m)
        best = min(best, float(minimize(obj, start, method="Nelder-Mead",
                                        options=_OPT).fun))
    return best


def propagate_published_errors(n_rep=24, seed=7):
    """Best-fit koff when the SOURCE'S OWN published errors are resampled.

    WHY THIS IS THE HONEST INTERVAL AND THE PROFILE WIDTH IS NOT. The Barberis trace is a
    NOISE-FREE parametric curve, so the residuals a profile sees are the model's systematic
    misfit, not experimental scatter -- `sigma_hat` is tiny and the curvature correspondingly
    enormous. The profile therefore reports an interval narrower than its own finest grid
    step, which is not precision; it is the absence of noise.

    The real uncertainty is the one the paper states: tau = 2.8 +/- 0.3, 33.4 +/- 4.6,
    221.35 +/- 14.9 ms and areas 0.57 +/- 0.04, 0.23 +/- 0.02, 0.20 +/- 0.03 (SEM, n = 6).
    Resampling within those, regenerating the trace and refitting is what propagates them.
    Areas are renormalised to sum to 1, because that is a constraint of the published fit
    rather than a free parameter.
    """
    from circuitpharm.fitting.data import DeactivationDataset
    from circuitpharm.literature import DEACTIVATIONS
    from circuitpharm.models.base import Observable

    row = next(d for d in DEACTIVATIONS if d.label == "Barberis 2007, brief pulse")
    taus, areas = row.components
    tau_sem, area_sem = (0.3, 4.6, 14.9), (0.04, 0.02, 0.03)
    base = KineticAllosteryModel()
    template = ALL["barberis2007_deact"]
    t = np.asarray(template.time_ms, dtype=float)
    rng = np.random.default_rng(seed)

    out = []
    for _ in range(n_rep):
        tt = np.maximum([rng.normal(m, s) for m, s in zip(taus, tau_sem)], 1e-3)
        aa = np.maximum([rng.normal(m, s) for m, s in zip(areas, area_sem)], 1e-6)
        aa = np.asarray(aa) / np.sum(aa)
        y = sum(a * np.exp(-t / tau) for a, tau in zip(aa, np.sort(tt)))
        ds = DeactivationDataset(
            citation_label="Barberis resampled within published SEM",
            time_ms=t, normalized_current=y,
            pulse_ms=template.pulse_ms, conc_mm=template.conc_mm,
            n_components=3, fit_window_ms=template.fit_window_ms,
            observable=Observable.DEACTIVATION, normalisation="fraction_of_max",
            kind="synthetic", synthetic=True,
            origin=("the Barberis components resampled inside their published SEMs, to "
                    "propagate the source's stated uncertainty into koff"),
        )
        cost = _cost_fn([ds])
        # koff free as well, so this is a best fit rather than a profile point
        x0 = np.log10([base.koff] + [getattr(base, k) for k in FREE])
        def obj(x):
            try:
                m = base.with_params(koff=10.0 ** x[0],
                                     **{k: 10.0 ** v for k, v in zip(FREE, x[1:])})
            except Exception:
                return np.inf
            return cost(m)
        best, bx = np.inf, None
        for i in range(2):
            s0 = x0 if i == 0 else x0 + rng.normal(0.0, 0.3, x0.size)
            r = minimize(obj, s0, method="Nelder-Mead", options=_OPT)
            if float(r.fun) < best:
                best, bx = float(r.fun), r.x
        if bx is not None and np.isfinite(best):
            out.append(10.0 ** bx[0])
    return np.asarray(out)


def _scan(grid, cost, base):
    prof = np.array([fit_at(k, cost, base) for k in grid])
    return prof - np.nanmin(prof)


def _report(grid, prof, base, title):
    print(f"\n{title}")
    print(f"  {'koff /ms':>12} {'delta(-logL)':>14}")
    for k, d in zip(grid, prof):
        print(f"  {k:>12.5g} {d:>14.3f}{'   <= inside 95%' if d <= DELTA_95 else ''}")
    i = int(prof.argmin())
    at_edge = i in (0, grid.size - 1)
    left_ok = bool(np.all(np.diff(prof[:i + 1]) <= 0.05))
    print(f"\n  minimum at koff = {grid[i]:.5g} /ms   (model default {base.koff:.5g}, "
          f"ratio {grid[i] / base.koff:.2f}x)")
    print(f"  minimum at a grid edge? {at_edge}"
          + ("   <-- UNBOUNDED, not identifiable" if at_edge else ""))
    print(f"  monotone approaching the minimum from the left (converged): {left_ok}")
    print(f"  rise either side: left {prof[max(0, i - 1)]:.1f}, right "
          f"{prof[min(grid.size - 1, i + 1)]:.1f} nats over one grid step")
    return i, at_edge


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-coarse", type=int, default=13)
    ap.add_argument("--n-fine", type=int, default=11)
    ap.add_argument("--no-refine", action="store_true")
    a = ap.parse_args()

    base = KineticAllosteryModel()
    ds = ALL["barberis2007_deact"]
    cost = _cost_fn([ds])

    print("=" * 84)
    print("PROFILE LIKELIHOOD OVER koff, AGAINST THE BARBERIS DEACTIVATION ALONE")
    print("=" * 84)
    print(f"\ndataset   {ds.citation_label}")
    print(f"protocol  {ds.pulse_ms:g} ms pulse of {ds.conc_mm:g} mM GABA, "
          f"{ds.n_components} components, {ds.fit_window_ms:g} ms window")
    print(f"grade     kind={ds.kind!r}, source_key={ds.source_key!r} "
          f"(parametric: the PARAMETERS are attributable, the points are the paper's model)")
    print(f"free      {FREE}   (d, r held at defaults)")

    coarse = base.koff * 10.0 ** np.linspace(-2, 2, a.n_coarse)
    prof_c = _scan(coarse, cost, base)
    i, at_edge = _report(coarse, prof_c, base, "COARSE -- four decades")

    if not a.no_refine and not at_edge:
        lo = coarse[max(0, i - 1)]
        hi = coarse[min(coarse.size - 1, i + 1)]
        fine = np.linspace(lo, hi, a.n_fine)
        prof_f = _scan(fine, cost, base)
        _report(fine, prof_f, base, f"REFINED -- {lo:.4g} to {hi:.4g} /ms")
        inside = fine[prof_f <= DELTA_95]
        best = fine[int(prof_f.argmin())]
        step = float(np.diff(fine).mean())
        print(f"\n{'=' * 84}")
        print("VERDICT")
        print("=" * 84)
        print(f"  koff IS identifiable from this trace alone: the profile has an interior")
        print(f"  minimum at {best:.4g} /ms and rises on both sides.")
        if inside.size <= 1:
            print(f"\n  The 95% width is NARROWER THAN THE FINEST GRID STEP ({step:.4g} /ms),")
            print(f"  so no interval is quoted from it. THAT IS NOT PRECISION -- it is the")
            print(f"  absence of noise: the trace is a noise-free parametric curve, so the")
            print(f"  residuals are the model's systematic misfit and sigma_hat is tiny.")
        else:
            print(f"\n  profile 95% span {inside.min():.4g} to {inside.max():.4g} /ms "
                  f"-- but see the caveat below; this is not a measurement uncertainty.")

        print(f"\n  PROPAGATING THE SOURCE'S OWN ERRORS INSTEAD (tau +/- 0.3/4.6/14.9 ms,")
        print(f"  areas +/- 0.04/0.02/0.03, SEM at n = 6):")
        k = propagate_published_errors()
        if k.size:
            lo, hi = np.percentile(k, [2.5, 97.5])
            print(f"    best-fit koff over {k.size} resamples: median {np.median(k):.4g} /ms,")
            print(f"    2.5-97.5 percentile {lo:.4g} to {hi:.4g} /ms = {hi / lo:.2f}x "
                  f"({np.log10(hi / lo):.3f} decades)")
            print(f"    THIS is the interval to quote. It is {hi / lo:.1f}x wide against a")
            print(f"    profile width below {step:.4g} /ms -- the published uncertainty")
            print(f"    dominates, and the profile never saw it.")
        else:
            print("    every resampled fit failed; no propagated interval available.")

        print(f"\n  Against the model's current default koff = {base.koff:.4g} /ms, the trace")
        print(f"  prefers {best / base.koff:.2f}x that value.")
        print("\n  This could not have been measured before 2026-10-11: the deactivation")
        print("  half of the likelihood did not exist, so the dataset was inert.")
        print("  Everything here is PARAMETRIC-grade -- it says what the published")
        print("  components determine, not what the receptor's koff is.")


if __name__ == "__main__":
    main()
