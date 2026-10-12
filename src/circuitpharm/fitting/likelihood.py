"""Per-observable likelihood, MLE with a reported spread, and a posterior that gates itself.

WHAT THIS FIXES relative to the objective in `identifiability.py` (roadmap P4-1):

  * EACH DATASET CONTRIBUTES ON ITS OWN OBSERVABLE. A PEAK concentration-response is
    compared against `peak_dose_response`, an EQUILIBRIUM one against `dose_response`. The
    previous objective compared a PEAK dataset against the equilibrium curve, which is the
    mismatch that drives D to its bound by deleting desensitisation (P3 §2.2).
  * SIGMA IS ESTIMATED, NOT ASSUMED. One scale per dataset, profiled out analytically:
    for a Gaussian with unknown scale the MLE is sigma^2 = RSS/n, and substituting it gives
    a concentrated log-likelihood that depends only on the parameters. So no AIC difference
    or Akaike weight depends on a hand-set noise level (roadmap C5, item 3). Each estimated
    scale still counts as a parameter for P5's k.
  * NORMALISATION CONSTANTS ARE INCLUDED, because P5 compares likelihoods ACROSS models
    and a dropped constant is only harmless within one.
  * THE HOLDOUT IS UNREACHABLE. `holdout_guard` raises rather than relying on a caller to
    remember, because P5's cross-validation and P6's falsification bound are both destroyed
    by a fit that has seen it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..models.base import Observable, ObservableMismatch
from .data import assert_real_data, holdout_guard
from .identifiability import IDENTIFIABLE_BOUNDS
from .reparam import EQUILIBRIUM_IDENTIFIABLE, IdentifiableParams


@dataclass(frozen=True)
class DatasetFit:
    """One dataset's contribution, with the noise scale that was estimated for it."""

    label: str
    observable: Observable
    n: int
    rss: float
    sigma_hat: float
    log_likelihood: float


@dataclass(frozen=True)
class LikelihoodResult:
    total_log_likelihood: float
    per_dataset: tuple[DatasetFit, ...]
    #: Estimated parameters plus one noise scale per dataset. This is the k an information
    #: criterion needs -- not the number of fields on a dataclass (roadmap C5, item 2).
    n_estimated: int


#: Simulation grid for a predicted deactivation. The pulse is resolved densely because it is
#: short (2 ms in the Barberis protocol) and sets the peak the decay is normalised by; the
#: tail is log-spaced because the components it has to resolve span two decades. These are
#: INTEGRATION settings, not data -- the prediction is interpolated onto the dataset's own
#: grid, so changing them must not move the score. `tests/test_deactivation_likelihood.py`
#: asserts that.
_DEACT_N_PULSE = 400
_DEACT_N_TAIL = 1600


def predict_deactivation(model, ds, *, n_pulse=_DEACT_N_PULSE, n_tail=_DEACT_N_TAIL):
    """Predicted deactivation decay for `ds`'s DECLARED protocol, normalised to its own peak.

    WHY THIS EXISTS. Until now nothing in `src/` read a deactivation trace's `time_ms` or
    `normalized_current` at all. `MISSING_DATASETS` describes the deactivation dataset as
    "the dataset that makes an absolute rate (k_off) identifiable at all", and
    `DEACTIVATION_BENCHMARK`'s own `origin` says it "exists so the kinetic-timescale half of
    the likelihood can be exercised" -- but that half of the likelihood did not exist, so a
    deactivation dataset was inert however good it was. `_predict` raised
    `ObservableMismatch` for anything that was not PEAK or EQUILIBRIUM.

    THE PROTOCOL COMES FROM THE DATASET, NEVER FROM A DEFAULT. `pulse_ms` and `conc_mm` are
    required fields on a DEACTIVATION dataset precisely so this function cannot invent them:
    Barberis 2007 measured the same receptor at 52.5 ms (2 ms pulse) and 364 ms (3 s pulse),
    6.9x apart, so a default pulse length here would silently pick which measurement the
    model is being compared against.

    The decay is re-zeroed at the model's PEAK and divided by it, because that is what the
    source's fit is: a triple exponential whose areas sum to 1 at the start of the decay.
    """
    # THE GUARD IS ON `q_matrix`, NOT ON `simulate_waveform`, and the difference matters.
    # Every model here HAS a `simulate_waveform`: `OperationalScalarModel` inherits one that
    # relaxes towards an instantaneous Hill target at an assumed time constant. So it would
    # happily return a "deactivation" that is a property of that invented relaxation rate
    # rather than of any binding or gating step -- scored against Barberis and relabelled on
    # the way, which is exactly the mismatch the P1 contract exists to stop. A deactivation
    # shape only means something for a model with states to move between, so the requirement
    # is a Markov generator. (It is also what stops Model A's relaxation integrator being
    # driven over this grid, where it is pathologically slow -- but slowness is the symptom,
    # not the reason.)
    if not hasattr(model, "q_matrix"):
        raise ObservableMismatch(
            f"{ds.citation_label!r} is DEACTIVATION data -- a decay SHAPE -- but "
            f"{type(model).__name__} has no `q_matrix`, so it has no states to relax "
            f"between. Its `simulate_waveform` relaxes towards a Hill target at an assumed "
            f"time constant, so any decay it produces measures that assumption and not a "
            f"rate. Compare it on PEAK, which every scheme here can produce, or drop the "
            f"model that has no such quantity.")
    pulse_ms = getattr(ds, "pulse_ms", None)
    conc_mm = getattr(ds, "conc_mm", None)
    window_ms = getattr(ds, "fit_window_ms", None)
    if not (pulse_ms and conc_mm and window_ms):
        raise ValueError(
            f"{ds.citation_label!r} is DEACTIVATION data but does not declare "
            f"pulse_ms={pulse_ms!r}, conc_mm={conc_mm!r}, fit_window_ms={window_ms!r}. "
            f"All three are part of the observable, not metadata about it.")

    t_sim = np.unique(np.concatenate([
        np.linspace(0.0, float(pulse_ms), n_pulse),
        float(pulse_ms) + np.geomspace(1e-4, float(window_ms), n_tail),
    ]))
    gaba_um = np.where(t_sim <= float(pulse_ms), float(conc_mm) * 1000.0, 0.0)
    p_open = np.asarray(model.simulate_waveform(t_sim, gaba_um).p_open, dtype=float)

    i_peak = int(np.argmax(p_open))
    peak = float(p_open[i_peak])
    if not (np.isfinite(peak) and peak > 0.0):
        raise ValueError(
            f"{ds.citation_label!r}: the model's peak open probability under its declared "
            f"protocol is {peak!r}, so there is no decay to normalise. A zero peak is a "
            f"failed prediction, not a very small one.")
    t_decay = t_sim[i_peak:] - t_sim[i_peak]
    y_decay = p_open[i_peak:] / peak
    return np.interp(np.asarray(ds.time_ms, dtype=float), t_decay, y_decay)


def _predict(model, ds):
    """Model prediction for one dataset, ON THAT DATASET'S OBSERVABLE AND SCALE."""
    obs = ds.observable
    if obs is Observable.PEAK:
        # every scheme here genuinely predicts a peak: B and C by integrating a square
        # application, A because its Hill curve IS its peak curve. So PEAK is an axis all
        # three can be compared on, and this branch needs no native-observable check.
        pred = np.asarray(model.peak_dose_response(ds.concs_um, pam_factor=1.0), float)
    elif obs is Observable.EQUILIBRIUM:
        # EQUILIBRIUM IS NOT. `dose_response` means "this model's native observable", and
        # `OperationalScalarModel.dose_response(concs, pam)` called positionally returns
        # its PEAK Hill curve without complaint -- its `observable` guard is a keyword
        # argument defaulting to PEAK, so the raise that exists for exactly this case never
        # fires through the base-class signature. A PEAK-native model would therefore have
        # its peak curve scored against equilibrium data and relabelled on the way, which
        # is the mismatch the P1 contract exists to stop. Checked here rather than in one
        # caller, because "the guard was applied in one of several places" is this
        # project's recurring error E12.
        if model.native_observable is not Observable.EQUILIBRIUM:
            raise ObservableMismatch(
                f"{ds.citation_label!r} is EQUILIBRIUM data but {type(model).__name__} is "
                f"{model.native_observable.value}-native, so it has no equilibrium curve "
                f"to be scored against -- `dose_response` would hand back its peak curve "
                f"under an equilibrium label. Compare on PEAK, which every scheme here can "
                f"produce, or drop the model that has no such quantity.")
        pred = np.asarray(model.dose_response(ds.concs_um, pam_factor=1.0), float)
    elif obs is Observable.DEACTIVATION:
        # A time course, not a curve. Already normalised to its own peak by construction, so
        # it returns here rather than falling through to the concentration-based
        # normalisation below -- a deactivation dataset has no `concs_um` to index.
        return predict_deactivation(model, ds)
    else:
        raise ObservableMismatch(
            f"{ds.citation_label!r} is tagged {obs.value}; this likelihood handles "
            f"concentration-response datasets (PEAK or EQUILIBRIUM) and DEACTIVATION time "
            f"courses. A CHARGE dataset is an integral over a window rather than a decay "
            f"shape, and still needs its own door -- see protocols/waveforms.py.")

    if getattr(ds, "normalisation", "absolute") == "fraction_of_max":
        # The data is I/I_max, so the prediction must be P/P(c_max) -- divided at the SAME
        # concentration the experimenters normalised by, which is the largest one in the
        # dataset. Comparing an absolute prediction to a normalised curve asks the model to
        # reach an asymptote of 1 that this scheme cannot have, and the optimiser obliges by
        # deleting desensitisation (measured: D -> 1e-11). See `data.Normalisation`.
        #
        # The consequence is intended: a normalised dataset then constrains only the SHAPE
        # of the curve, which is the only thing it carries information about. Absolute open
        # probability has to come from somewhere else.
        ref = float(pred[int(np.argmax(np.asarray(ds.concs_um, dtype=float)))])
        if not (np.isfinite(ref) and ref > 0.0):
            raise ValueError(
                f"{ds.citation_label!r} is normalised, but the prediction at its highest "
                f"concentration is {ref!r}, so there is nothing to normalise by. A zero "
                f"maximum is a failed prediction, not a very small one.")
        pred = pred / ref
    return pred


def _observed(ds):
    # A deactivation dataset carries a decay, not a response column. Dispatched on the
    # declared observable rather than on which attribute happens to exist, so a dataset that
    # grows a `mean_response` later cannot silently change which column is scored.
    obs = getattr(ds, "observable", None)
    if obs is Observable.DEACTIVATION:
        return np.asarray(ds.normalized_current, dtype=float)
    if not hasattr(ds, "mean_response"):
        # Reached before `_predict` can explain itself, because the observed column is read
        # first. Without this, a CHARGE dataset fails with
        # "AttributeError: 'DeactivationDataset' object has no attribute 'mean_response'",
        # which says nothing about the actual problem -- the observable has no door here.
        raise ObservableMismatch(
            f"{getattr(ds, 'citation_label', ds)!r} is tagged "
            f"{getattr(obs, 'value', obs)!r} and carries no response column this likelihood "
            f"can read. PEAK and EQUILIBRIUM datasets supply `mean_response`; DEACTIVATION "
            f"supplies `normalized_current`. A CHARGE dataset is an integral over a window "
            f"rather than either, and still needs its own door -- see protocols/waveforms.py.")
    return np.asarray(ds.mean_response, dtype=float)


def make_log_likelihood(datasets, *, caller: str, n_params: int):
    """Run the guards ONCE and return a scoring closure.

    THE GUARDS BELONG HERE, NOT IN THE INNER LOOP. `assert_real_data` formats a
    multi-paragraph `ProvisionalResultWarning` every time it fires and `provisional.py`
    installs `simplefilter("always")` for that category precisely so it cannot be shown
    once and then suppressed. An optimiser calling the scorer ten thousand times would
    therefore pay for ten thousand warnings and bury the one that matters -- a warning
    nobody can read is the same as no warning. So the taint check and the holdout check
    run once, against the dataset objects the whole fit will use, and the returned closure
    is pure arithmetic.

    Returns `score(model) -> LikelihoodResult`.
    """
    datasets = list(datasets)
    holdout_guard(datasets, caller=caller)
    assert_real_data(datasets, caller=caller)

    def score(model) -> LikelihoodResult:
        return _score(model, datasets, n_params=n_params)

    return score


def concentrated_log_likelihood(model, datasets, *, n_params: int) -> LikelihoodResult:
    """Gaussian log-likelihood with each dataset's noise scale profiled out.

    For `n` points with residual sum of squares `RSS` and unknown scale, the profile MLE is
    `sigma^2 = RSS/n` and the concentrated log-likelihood is

        -n/2 * (log(2*pi) + 1 + log(RSS/n))

    which is exact, carries its constants, and leaves nothing to tune.

    For a single evaluation. Inside a fit or a chain, call `make_log_likelihood` once and
    reuse the closure, so the synthetic-data warning fires once rather than per iteration.
    """
    return make_log_likelihood(datasets, caller="concentrated_log_likelihood",
                               n_params=n_params)(model)


def score_holdout(model, datasets, *, caller: str) -> LikelihoodResult:
    """Score an ALREADY-FITTED model against held-out data. The guard runs backwards here.

    `make_log_likelihood` refuses held-out datasets, because anything that can fit must not
    see them. Scoring a finished model against them is the one thing the holdout exists
    FOR, so it needs its own door -- and that door refuses TRAINING data, so neither
    function can be called with the wrong set and a typo cannot turn an out-of-sample score
    into an in-sample one.

    `n_params` is deliberately not a parameter: a held-out score is a predictive density,
    not an information criterion, and nothing is estimated here.
    """
    datasets = list(datasets)
    trained_on = [getattr(d, "citation_label", repr(d)) for d in datasets
                  if getattr(d, "role", "train") != "holdout"]
    if trained_on:
        raise ValueError(
            f"{caller} scored {trained_on} as holdout, but they are tagged role='train'. "
            f"An in-sample score reported as out-of-sample is the failure this whole split "
            f"exists to prevent; pass `fitting.data.HOLDOUT`.")
    assert_real_data(datasets, caller=f"{caller} (holdout score)")
    return _score(model, datasets, n_params=0)


def _score(model, datasets, *, n_params: int) -> LikelihoodResult:
    fits, total = [], 0.0
    for ds in datasets:
        y = _observed(ds)
        pred = _predict(model, ds)
        if pred.shape != y.shape:
            raise ValueError(
                f"{ds.citation_label!r}: prediction shape {pred.shape} does not match "
                f"observation shape {y.shape}")
        n = int(y.size)
        rss = float(np.sum((pred - y) ** 2))
        if not np.isfinite(rss):
            raise ValueError(
                f"{ds.citation_label!r}: residual sum of squares is {rss!r}. A non-finite "
                f"likelihood is a failure, not a very bad fit; it is not scored.")
        # a perfect fit sends sigma_hat to 0 and the log-likelihood to +inf, which is a
        # degenerate optimum rather than a good one. Floored at a level far below any real
        # measurement precision, and the floor is reported.
        sigma2 = max(rss / n, 1e-12)
        ll = -0.5 * n * (np.log(2.0 * np.pi) + 1.0 + np.log(sigma2))
        fits.append(DatasetFit(label=ds.citation_label, observable=ds.observable, n=n,
                               rss=rss, sigma_hat=float(np.sqrt(sigma2)),
                               log_likelihood=float(ll)))
        total += ll

    return LikelihoodResult(total_log_likelihood=float(total), per_dataset=tuple(fits),
                            n_estimated=int(n_params) + len(fits))


@dataclass(frozen=True)
class MLEResult:
    """An MLE, with the evidence that it IS one.

    `top_spread` is the range of the best five optima in each coordinate. If it exceeds the
    profile-likelihood intervals, the optimisation has not converged on a unique answer and
    reporting the single best one would be picking a winner from a tie (roadmap P4-2).
    """

    #: Fitted values by name, in whatever parameterisation was fitted.
    values: dict[str, float]
    #: The same thing as an `IdentifiableParams`, for the equilibrium-identifiable triple
    #: only. `None` for any other parameterisation -- a model with two desensitisation
    #: sinks has four combinations and does not fit in that dataclass.
    params: IdentifiableParams | None
    log_likelihood: float
    n_starts: int
    n_converged: int
    top_spread: dict[str, float]
    per_dataset: tuple[DatasetFit, ...]
    converged_uniquely: bool
    #: Per parameter: did the optimum land on its lower / upper bound? A parameter at its
    #: bound was not estimated -- the box stopped it -- and reporting it as a fitted value
    #: is reporting the box. On the normalised Jahn curve `log10_D` does exactly this.
    at_bound: dict[str, tuple[bool, bool]]
    bounds: dict[str, tuple[float, float]]
    note: str


def fit_mle(model_factory, datasets, *, n_starts: int = 32, seed: int = 0,
            bounds=None, spread_tol: float = 0.05) -> MLEResult:
    """Multi-start MLE over the equilibrium-identifiable combinations.

    `model_factory` takes an `IdentifiableParams` and returns a model. The caller supplies
    it, because turning three ratios into six rates needs conventions and this module will
    not choose them (see `reparam.IdentifiableParams.to_microscopic`).

    A thin wrapper over `fit_mle_vector`, which is the general form: P5 fits three models
    in three different parameterisations and only one of them is this triple.
    """
    names = list(EQUILIBRIUM_IDENTIFIABLE)
    return fit_mle_vector(
        lambda v: model_factory(IdentifiableParams.from_vector(v)), datasets,
        param_names=names, bounds=bounds or [IDENTIFIABLE_BOUNDS[k] for k in names],
        n_starts=n_starts, seed=seed, spread_tol=spread_tol, caller="fit_mle")


def fit_mle_vector(model_factory, datasets, *, param_names, bounds,
                   n_starts: int = 32, seed: int = 0, spread_tol: float = 0.05,
                   caller: str = "fit_mle_vector") -> MLEResult:
    """Multi-start MLE over an arbitrary parameter vector, reporting the top-5 spread.

    `model_factory` takes a plain `np.ndarray` in the order of `param_names` and returns a
    model. Every parameter is expected to be in log10 -- the spread and the bound checks
    are reported in decades and would be meaningless on a linear axis -- and `bounds` is
    required, because a search with no box is how an optimiser reports a desensitisation
    ratio of 1e-9 as a fit.
    """
    from scipy.optimize import minimize
    from scipy.stats import qmc

    datasets = list(datasets)
    names = list(param_names)
    bnds = list(bounds)
    if len(bnds) != len(names):
        raise ValueError(
            f"one bound per parameter: got {len(bnds)} for {len(names)} names {names}")
    score = make_log_likelihood(datasets, caller=caller, n_params=len(names))

    def neg_ll(v):
        try:
            return -score(model_factory(np.asarray(v, dtype=float))).total_log_likelihood
        except Exception:
            return 1e12

    # Latin hypercube rather than uniform random: 32 points over three decades-wide axes
    # cover far better, and the starts are reproducible from `seed`.
    sampler = qmc.LatinHypercube(d=len(bnds), seed=seed)
    lo = np.array([b[0] for b in bnds])
    hi = np.array([b[1] for b in bnds])
    starts = lo + sampler.random(n_starts) * (hi - lo)

    results = []
    for s in starts:
        r = minimize(neg_ll, s, method="L-BFGS-B", bounds=bnds,
                     options=dict(ftol=1e-14, gtol=1e-11, maxiter=4000))
        # THE POLISH RUNS WHETHER OR NOT L-BFGS-B SUCCEEDED. On a surface with a flat
        # direction it returns ABNORMAL_TERMINATION_IN_LNSRCH from a perfectly good point,
        # and treating that as a dead start threw away 9 of 12 starts on the Jahn curve --
        # which does not merely waste work, it computes the top-5 spread over three optima
        # and calls the result a tie between five.
        #
        # BOUNDS ON THE POLISH TOO. Without them Nelder-Mead walks straight out of the box:
        # measured, log10_D reached -9.03 against a stated lower bound of -3, so the
        # reported "fit" sat at a desensitisation ratio the search space excluded. An
        # optimiser silently leaving its feasible region is worse than one that fails,
        # because the number looks fitted.
        start2 = r.x if np.all(np.isfinite(r.x)) else s
        r2 = minimize(neg_ll, start2, method="Nelder-Mead", bounds=bnds,
                      options=dict(xatol=1e-9, fatol=1e-11, maxiter=2000))
        if r2.success and (not r.success or r2.fun <= r.fun):
            r = r2
        if r.success and np.isfinite(r.fun) and r.fun < 1e11:
            results.append((float(r.fun), np.asarray(r.x, dtype=float)))

    if not results:
        raise RuntimeError(
            f"none of {n_starts} starts converged. A failed optimisation is not a result "
            f"and is not reported as one; check the model factory and the bounds.")

    results.sort(key=lambda t: t[0])
    best_cost, best_x = results[0]
    top = np.vstack([x for _, x in results[:min(5, len(results))]])
    spread = {n: float(top[:, i].max() - top[:, i].min()) for i, n in enumerate(names)}
    unique = all(v <= spread_tol for v in spread.values())

    values = {n: float(best_x[i]) for i, n in enumerate(names)}
    ll = score(model_factory(best_x))
    params = (IdentifiableParams.from_vector(best_x)
              if names == list(EQUILIBRIUM_IDENTIFIABLE) else None)

    # "within 0.1% of the box width" rather than an exact comparison: L-BFGS-B returns a
    # point a hair inside the bound it is pressed against, so an equality test would never
    # fire and the check would be decorative.
    width = hi - lo
    at_bound = {n: (bool(best_x[i] - lo[i] <= 1e-3 * width[i]),
                    bool(hi[i] - best_x[i] <= 1e-3 * width[i]))
                for i, n in enumerate(names)}
    pinned = [n for n, p in at_bound.items() if any(p)]

    parts = []
    parts.append("the best five optima agree to within {:.3g} decades in every coordinate"
                 .format(max(spread.values())) if unique else
                 ("THE TOP FIVE OPTIMA DISAGREE by up to {:.3g} decades ({}), so this is "
                  "not a unique answer -- report the spread, not the winner"
                  .format(max(spread.values()), spread)))
    if pinned:
        parts.append("{} sit(s) ON a bound, so that value is the edge of the search space "
                     "and not an estimate -- the data does not constrain it"
                     .format(", ".join(pinned)))
    if len(results) < n_starts:
        parts.append(f"{n_starts - len(results)} of {n_starts} starts did not converge and "
                     f"are not scored")

    return MLEResult(values=values, params=params,
                     log_likelihood=ll.total_log_likelihood,
                     n_starts=n_starts, n_converged=len(results), top_spread=spread,
                     per_dataset=ll.per_dataset, converged_uniquely=unique,
                     at_bound=at_bound,
                     bounds={n: (float(lo[i]), float(hi[i])) for i, n in enumerate(names)},
                     note="; ".join(parts))
