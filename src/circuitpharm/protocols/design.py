"""The discriminating experiment, scored against parameter uncertainty (roadmap P6-3).

WHAT THIS REPLACES, AND WHY IT IS NOT A REFINEMENT. `oed.find_discriminating_protocol`
scores a protocol as

    |y_A - y_B| / sigma_measurement

at each model's point prediction. That quantity cannot see parameter uncertainty at all,
so it reports a protocol as maximally discriminating precisely where the two models'
PREDICTIONS are far apart -- which, for posteriors as wide as this project's, is routinely
where each model's own predictive interval is wider still. A protocol chosen that way
cannot refute anything, and the +/- 2*sigma_noise it offers as a falsification interval is
narrower than the real spread in exactly the direction that makes a pre-registration
overconfident. It stays callable and VOID; this is the version whose numbers mean something.

The score here is a posterior-predictive separation:

    score = E|y_A - y_B| / sqrt( var(y_A) + var(y_B) + sigma_meas^2 )

with `y_A`, `y_B` drawn from each model's own posterior. Three properties the old one did
not have:

  * IT FALLS WHEN A POSTERIOR WIDENS, because the denominator grows. That is the regression
    guard the roadmap names (C6), and it is the whole point: a protocol is only
    discriminating if the models disagree by more than either of them is unsure.
  * IT IS DIMENSIONLESS AND INTERPRETABLE. A score of 1 means the expected difference
    equals one standard deviation of the thing you would measure, so roughly 1.96 is the
    smallest value at which a single measurement could tell the models apart at all.
  * WHAT TO MEASURE IS PART OF THE DESIGN. The search ranges over the observable as well as
    over agonist and modulator, because the answer to "which experiment" is sometimes
    "a different quantity", and a grid that fixes the observable cannot return it.

AND THE MULTIPLE-COMPARISON COST IS CARRIED, not mentioned. The best of `n` searched
protocols is biased upward; `DesignResult.n_searched` and `score_penalised` report the
score after a Bonferroni-style correction on the implied z, so a pre-registration quotes
the corrected number.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from ..models.base import Observable, ObservableMismatch
from ..provisional import void
from ..results import Quantity, ResultSet, Tier

#: Measurement noise on a NORMALISED open-probability scale, per observation. A stated
#: convention, not a measurement: it is the one number in this module that does not come
#: from the posterior, and the score is reported against a sweep of it rather than at a
#: single value for exactly that reason (see `sigma_sensitivity`).
DEFAULT_SIGMA_MEAS = 0.02


@dataclass(frozen=True)
class ProtocolPoint:
    """One candidate experiment."""

    gaba_um: float
    pam_factor: float
    observable: Observable
    application_ms: float = 300.0

    def describe(self) -> str:
        return (f"{self.gaba_um:g} uM GABA, PAM {self.pam_factor:g}x (EQUILIBRIUM EC50 "
                f"shift), measured as {self.observable.value}"
                + (f" over a {self.application_ms:g} ms application"
                   if self.observable is Observable.PEAK else ""))


@dataclass(frozen=True)
class DesignResult:
    """The chosen protocol, its score, and what each model predicts there."""

    protocol: ProtocolPoint
    score: float
    #: The score after a multiple-comparison correction for having searched `n_searched`
    #: protocols. See `_penalise`.
    score_penalised: float
    n_searched: int
    #: Per model: the posterior-predictive median and 95% interval at the protocol. These
    #: are the falsification boundaries, and they are two-sided and include both parameter
    #: uncertainty and measurement noise.
    predictive_intervals: dict[str, tuple[float, float, float]]
    expected_separation: float
    predictive_sd: dict[str, float]
    sigma_meas: float
    n_draws: int
    #: The runner-up protocol and its score, so a reader can see whether the maximum is a
    #: peak or a plateau. A design chosen off a plateau is not the design, it is one of
    #: many equally good ones.
    runner_up: tuple[str, float] | None
    note: str

    @property
    def discriminating(self) -> bool:
        """Is the CORRECTED score above the 1.96 a single measurement would need?"""
        return self.score_penalised >= 1.96


def _penalise(score: float, n_searched: int) -> float:
    """Discount the winning score for the size of the search.

    The best of `n` protocols is biased upward, and the bias is not small: at `n = 200` the
    maximum of standard normals sits near 2.75 sigma, so a raw score of 3 is barely
    distinguishable from noise in the SEARCH. Treating the score as a z-statistic, the
    Bonferroni-corrected equivalent is the z whose one-sided tail probability is
    `n * P(Z > score)`. The result is a score that can be quoted in a pre-registration
    without the reader having to know how many protocols were tried.

    Returns `score` unchanged for `n_searched <= 1`, and 0.0 when the correction exhausts
    the tail -- which is itself the finding: the search found nothing the search cannot
    explain.
    """
    from scipy.stats import norm

    n = max(int(n_searched), 1)
    if n == 1 or not np.isfinite(score):
        return float(score)
    tail = float(norm.sf(score)) * n
    if tail >= 1.0:
        return 0.0
    return float(max(0.0, norm.isf(tail)))


def _predict_draws(model_factory, draws: np.ndarray, point: ProtocolPoint) -> np.ndarray:
    """Each posterior draw's prediction at one protocol, on the protocol's observable."""
    out = np.empty(draws.shape[0], dtype=float)
    conc = np.array([point.gaba_um], dtype=float)
    for i, row in enumerate(draws):
        try:
            m = model_factory(row)
            if point.observable is Observable.PEAK:
                val = float(np.atleast_1d(m.peak_dose_response(
                    conc, pam_factor=point.pam_factor,
                    application_ms=point.application_ms))[0])
            elif point.observable is Observable.EQUILIBRIUM:
                if m.native_observable is not Observable.EQUILIBRIUM:
                    raise ObservableMismatch(
                        f"{type(m).__name__} is {m.native_observable.value}-native and has "
                        f"no equilibrium response; design on PEAK instead.")
                val = float(m.steady_state(point.gaba_um, pam_factor=point.pam_factor))
            else:
                raise ObservableMismatch(
                    f"{point.observable.value} is not a single-number response; a CHARGE "
                    f"design needs a waveform and a window (protocols/waveforms.py).")
            out[i] = val if np.isfinite(val) else np.nan
        except ObservableMismatch:
            raise
        except Exception:
            out[i] = np.nan
    return out


def posterior_predictive_score(
    candidates: dict[str, tuple],
    point: ProtocolPoint,
    *,
    sigma_meas: float = DEFAULT_SIGMA_MEAS,
    rng_seed: int = 0,
) -> tuple[float, dict[str, tuple[float, float, float]], dict[str, float], float]:
    """Score one protocol. `candidates` maps a model name to `(factory, draws)`.

    Returns `(score, intervals, predictive_sds, expected_separation)`. Exactly two models:
    the score is a pairwise separation and averaging it over three pairs would hide which
    pair the protocol actually separates.
    """
    if len(candidates) != 2:
        raise ValueError(
            f"a discrimination score is pairwise; got {len(candidates)} models "
            f"{sorted(candidates)}. Score each pair and report them separately -- an "
            f"average over pairs cannot say which two the protocol tells apart.")
    if not (np.isfinite(sigma_meas) and sigma_meas > 0):
        raise ValueError(
            f"sigma_meas must be finite and positive, got {sigma_meas!r}: a zero "
            f"measurement noise makes every protocol infinitely discriminating.")

    (name_a, (fac_a, draws_a)), (name_b, (fac_b, draws_b)) = candidates.items()
    ya = _predict_draws(fac_a, np.atleast_2d(draws_a), point)
    yb = _predict_draws(fac_b, np.atleast_2d(draws_b), point)
    ya, yb = ya[np.isfinite(ya)], yb[np.isfinite(yb)]
    if ya.size < 2 or yb.size < 2:
        raise RuntimeError(
            f"fewer than two usable draws at {point.describe()} ({ya.size} and {yb.size}); "
            f"a predictive spread cannot be estimated and the protocol is not scored.")

    # the posteriors are independent samples of different sizes, so the expected absolute
    # difference is taken over matched pairs drawn from each -- not over the outer product,
    # which is the same expectation at n^2 cost.
    n = min(ya.size, yb.size)
    rng = np.random.default_rng(rng_seed)
    pairs_a = rng.permutation(ya)[:n]
    pairs_b = rng.permutation(yb)[:n]
    noise = rng.normal(0.0, sigma_meas, (2, n))
    expected_sep = float(np.mean(np.abs((pairs_a + noise[0]) - (pairs_b + noise[1]))))

    var_a = float(np.var(ya, ddof=1))
    var_b = float(np.var(yb, ddof=1))
    denom = float(np.sqrt(var_a + var_b + sigma_meas ** 2))
    score = expected_sep / denom if denom > 0 else float("inf")

    intervals = {}
    sds = {}
    for name, y in ((name_a, ya), (name_b, yb)):
        # the predictive interval includes measurement noise: it is what a single
        # experiment would return, which is what a pre-registration has to bound.
        draws_plus_noise = y + rng.normal(0.0, sigma_meas, y.size)
        intervals[name] = (float(np.percentile(draws_plus_noise, 2.5)),
                           float(np.median(draws_plus_noise)),
                           float(np.percentile(draws_plus_noise, 97.5)))
        sds[name] = float(np.std(draws_plus_noise, ddof=1))
    return float(score), intervals, sds, expected_sep


def find_discriminating_protocol_pp(
    candidates: dict[str, tuple],
    *,
    gaba_grid: Sequence[float] | None = None,
    pam_grid: Sequence[float] | None = None,
    observables: Sequence[Observable] = (Observable.PEAK, Observable.EQUILIBRIUM),
    application_grid: Sequence[float] = (300.0,),
    sigma_meas: float = DEFAULT_SIGMA_MEAS,
    rng_seed: int = 0,
) -> DesignResult:
    """Search protocol space for the maximum posterior-predictive separation.

    The search includes THE OBSERVABLE, because what to measure is part of the design: a
    PEAK protocol and an EQUILIBRIUM one at the same agonist and modulator are different
    experiments and can differ by more than any two agonist concentrations do.

    A protocol that raises `ObservableMismatch` -- a model with no such quantity -- is
    skipped with the observable recorded, not scored as zero: "this model cannot be
    measured that way" and "the two models agree" are different statements.
    """
    gaba = list(gaba_grid or (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0))
    pam = list(pam_grid or (1.0, 1.25, 1.5, 2.0, 2.5, 4.0))

    scored: list[tuple[float, ProtocolPoint, dict, dict, float]] = []
    skipped: dict[str, str] = {}
    for obs in observables:
        # WHETHER AN OBSERVABLE IS AVAILABLE AT ALL IS DECIDED ONCE, not rediscovered
        # inside the grid. A model that has no such quantity has none at every
        # concentration, and the earlier version broke out of the innermost loop and relied
        # on the enclosing ones happening to have a single iteration -- correct by accident.
        probe = ProtocolPoint(gaba_um=float(gaba[0]), pam_factor=float(pam[0]),
                              observable=obs,
                              application_ms=float(application_grid[0]))
        try:
            posterior_predictive_score(candidates, probe, sigma_meas=sigma_meas,
                                       rng_seed=rng_seed)
        except ObservableMismatch as exc:
            skipped[obs.value] = str(exc)
            continue
        except Exception:
            pass        # a bad probe point is not a reason to drop the observable

        apps = application_grid if obs is Observable.PEAK else (300.0,)
        for app in apps:
            for g in gaba:
                for p in pam:
                    pt = ProtocolPoint(gaba_um=float(g), pam_factor=float(p),
                                       observable=obs, application_ms=float(app))
                    try:
                        s, iv, sd, sep = posterior_predictive_score(
                            candidates, pt, sigma_meas=sigma_meas, rng_seed=rng_seed)
                    except Exception:
                        continue
                    scored.append((s, pt, iv, sd, sep))

    if not scored:
        raise RuntimeError(
            f"no protocol could be scored. Observables skipped because a model has no such "
            f"quantity: {skipped or 'none'}.")

    scored.sort(key=lambda t: -t[0])
    best_score, best_pt, intervals, sds, sep = scored[0]
    runner = ((scored[1][1].describe(), scored[1][0]) if len(scored) > 1 else None)
    penalised = _penalise(best_score, len(scored))
    n_draws = min(np.atleast_2d(d).shape[0] for _, d in candidates.values())

    note = (f"best of {len(scored)} protocols: {best_pt.describe()}. Score "
            f"{best_score:.3f}, {penalised:.3f} after correcting for the search. "
            + ("This clears the 1.96 a single measurement needs."
               if penalised >= 1.96 else
               "THIS DOES NOT CLEAR 1.96, so a single measurement at the best protocol in "
               "this space could not separate the models -- the honest conclusion is that "
               "the experiment does not exist yet, not that this is the experiment.")
            + (f" Observables skipped: {sorted(skipped)}." if skipped else ""))

    return DesignResult(
        protocol=best_pt, score=best_score, score_penalised=penalised,
        n_searched=len(scored), predictive_intervals=intervals,
        expected_separation=sep, predictive_sd=sds, sigma_meas=sigma_meas,
        n_draws=int(n_draws), runner_up=runner, note=note)


def sigma_sensitivity(candidates: dict[str, tuple], point: ProtocolPoint,
                      sigmas: Sequence[float] = (0.005, 0.01, 0.02, 0.05),
                      rng_seed: int = 0) -> dict[float, float]:
    """The score against measurement noise, because `sigma_meas` is the one assumed input.

    If the design's verdict flips across a plausible range of patch-clamp noise, the design
    is a statement about the amplifier and not about the receptor. Reported alongside the
    score so a reader can see which it is.
    """
    out = {}
    for s in sigmas:
        score, _, _, _ = posterior_predictive_score(candidates, point, sigma_meas=s,
                                                    rng_seed=rng_seed)
        out[float(s)] = float(score)
    return out


def design_report(result: DesignResult, *, posterior_tier: Tier = Tier.UNCALIBRATED,
                  posterior_note: str = "") -> ResultSet:
    """The design as tiered quantities, with the falsification bounds as the deliverable.

    `posterior_tier` is the tier of the POSTERIORS the draws came from and cannot be
    improved here: a falsification boundary from an unconverged chain is VOID however
    carefully the search was run.
    """
    rs = ResultSet(f"discriminating protocol ({result.note})")
    prov = (f"posterior-predictive separation E|y_A - y_B| / sqrt(var(y_A) + var(y_B) + "
            f"sigma^2) over {result.n_draws} draws per model at sigma_meas = "
            f"{result.sigma_meas:g}; best of {result.n_searched} protocols searched"
            + (f". {posterior_note}" if posterior_note else ""))
    caveats = (
        f"the search covered {result.n_searched} protocols, so the raw score "
        f"{result.score:.3f} is biased upward; {result.score_penalised:.3f} is the "
        f"corrected value and is the one to quote",
        f"sigma_meas = {result.sigma_meas:g} on a normalised open-probability scale is a "
        f"stated convention, not a measurement -- see `sigma_sensitivity`",
        "PAM strength is an EQUILIBRIUM EC50 shift here; gabaa_kinetics.calibrate_pam uses "
        "the PEAK shift and the two are not interchangeable",
    )

    def wrap(name, value, extra=()):
        if posterior_tier is Tier.VOID:
            return void(name, value,
                        defect=("drawn from a posterior the caller declared VOID: "
                                + (posterior_note or "its chain failed a diagnostic")),
                        register_item="P6-3",
                        promote_by="run the chain to convergence (fitting/posterior.py)",
                        caveats=caveats + tuple(extra))
        return Quantity(name=name, _value=value, tier=posterior_tier, provenance=prov,
                        promote_by=("fit both models to digitised data on their own "
                                    "observables, then re-run the design"),
                        caveats=caveats + tuple(extra))

    rs.add(wrap("discrimination_score", result.score_penalised,
                (f"uncorrected {result.score:.3f}",)))
    rs.add(wrap("protocol", result.protocol.describe()))
    for name, (lo, med, hi) in result.predictive_intervals.items():
        rs.add(wrap(f"falsification_interval_{name}", (lo, med, hi),
                    (f"two-sided 95% posterior-predictive interval including measurement "
                     f"noise; predictive SD {result.predictive_sd[name]:.4g}",
                     "an outcome outside BOTH models' intervals refutes both, which is the "
                     "outcome worth pre-registering")))
    return rs
