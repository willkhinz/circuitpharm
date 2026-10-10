"""Model comparison that means something (roadmap P5).

WHAT THIS REPLACES. `protocols/oed.py::evaluate_model_fit` scores each model at whatever
parameters it was handed, with `k = len(param_names)` and a hand-set sigma. All three are
defects rather than approximations, and the module docstring there lists them. It stays
callable, loudly, returning VOID conclusions, because it is what the manuscript's numbers
came from and deleting it would make those numbers unreproducible rather than wrong.

WHAT THIS DOES INSTEAD, in order:

  1. FITS EACH MODEL INDEPENDENTLY, in its own parameterisation, to its own MLE. An
     information criterion is a statement about a model AT its maximum; evaluated anywhere
     else it ranks starting guesses.
  2. COUNTS k AS ESTIMATED IDENTIFIABLE PARAMETERS PLUS ESTIMATED NOISE SCALES. Not
     `len(param_names)`: that charges Model C for a `pam_desens_factor` the pam_factor = 1
     curve never touches, and charges Model B for `d` and `r` separately when only their
     ratio enters -- a penalty wrong in the direction that systematically disfavours the
     more detailed scheme.
  3. REPORTS AICc ALONGSIDE AIC AND BIC. n is 9-30 points here. At n/k below ~40 the
     second-order correction is not a refinement, it is the difference between AIC and a
     criterion that holds.
  4. CROSS-VALIDATES BY CONCENTRATION BLOCK. Random point-wise folds leak across a smooth
     curve -- a neighbouring point is almost the answer -- and flatter the flexible model.
     Contiguous blocks make each fold an extrapolation, which is what out-of-sample means.
  5. SCORES THE PROSPECTIVE HOLDOUT ONCE, separately, and never uses it to choose.
  6. REPORTS AKAIKE WEIGHTS ONLY ALONGSIDE THE CV RANKING, AND ONLY IF THE TWO AGREE. When
     they disagree the comparison is reported as UNRESOLVED, which is a result.

AND THE METHOD IS CHECKED AGAINST A KNOWN ANSWER. `recovery_check` generates data from one
model at known parameters and asks the protocol to find it. A protocol that cannot recover
a truth it was given cannot adjudicate one it was not. That is the one legitimate use of
generated data here and it is labelled as a method check, never as evidence about
receptors.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

from ..models.base import Observable, ObservableMismatch, ReceptorModel
from ..provisional import void
from ..results import Quantity, ResultSet, Tier
from .data import assert_real_data, holdout_guard
from .likelihood import _predict, fit_mle_vector, make_log_likelihood, score_holdout

#: Conventions held fixed while fitting, with what they rest on. EVERY ONE IS A CONVENTION
#: AND NONE IS A MEASUREMENT. They are here because no dataset in this repository carries a
#: timescale (`data.MISSING_DATASETS`), so the absolute rates are not identifiable and the
#: fit would wander along a flat direction if they were free. Fixing them is the honest
#: alternative to pretending they were estimated -- but it does mean every absolute rate
#: below is conditional on these, and only the RATIOS are fitted.
FIT_CONVENTIONS = {
    "koff": (0.180415, "gabaa_kinetics.fit_scheme's three-anchor result"),
    "alpha": (0.30, "gabaa_kinetics.FIT_FIXED_ALPHA, a declared mean open time of 3.33 ms"),
    "r": (0.002, "Model B's own default; no recorded origin (P0-13)"),
    "r_fast": (0.0050, "Model C's own default; no recorded origin (P0-13)"),
    "r_slow": (0.0004, "Model C's own default; no recorded origin (P0-13)"),
}


@dataclass(frozen=True)
class ModelSpec:
    """How to fit one model: its free parameters, their box, and what is held fixed."""

    name: str
    #: Free parameters, all in log10. The spread and bound checks are in decades.
    param_names: tuple[str, ...]
    bounds: tuple[tuple[float, float], ...]
    #: Takes a vector in `param_names` order and returns a model.
    factory: Callable[[np.ndarray], ReceptorModel]
    #: What is held fixed, and why. Keys index `FIT_CONVENTIONS`.
    conventions: tuple[str, ...]
    note: str

    @property
    def n_free(self) -> int:
        return len(self.param_names)


def _model_b(v: np.ndarray) -> ReceptorModel:
    from ..models.kinetic_jw95 import KineticAllosteryModel

    kd, e, d = (10.0 ** float(x) for x in v)
    koff = FIT_CONVENTIONS["koff"][0]
    alpha = FIT_CONVENTIONS["alpha"][0]
    r = FIT_CONVENTIONS["r"][0]
    return KineticAllosteryModel(kon=koff / kd, koff=koff, beta=alpha * e, alpha=alpha,
                                 d=r * d, r=r)


def _model_c(v: np.ndarray) -> ReceptorModel:
    from ..models.extended_desens import ExtendedDesensitizationModel

    kd, e, df, ds = (10.0 ** float(x) for x in v)
    koff = FIT_CONVENTIONS["koff"][0]
    alpha = FIT_CONVENTIONS["alpha"][0]
    rf = FIT_CONVENTIONS["r_fast"][0]
    rs = FIT_CONVENTIONS["r_slow"][0]
    return ExtendedDesensitizationModel(
        kon=koff / kd, koff=koff, beta=alpha * e, alpha=alpha,
        d_fast=rf * df, r_fast=rf, d_slow=rs * ds, r_slow=rs)


def _model_a(v: np.ndarray) -> ReceptorModel:
    from ..models.operational import OperationalScalarModel

    ec50, n_h, po_max = (10.0 ** float(x) for x in v)
    return OperationalScalarModel(ec50_um=ec50, hill_n=n_h, po_max=min(po_max, 1.0))


#: The three models, each in the parameterisation its own structure admits.
MODEL_SPECS = {
    "operational": ModelSpec(
        name="operational",
        param_names=("log10_ec50_um", "log10_hill_n", "log10_po_max"),
        bounds=((-1.0, 4.0), (-0.3, 0.7), (-2.0, 0.0)),
        factory=_model_a,
        conventions=(),
        note=("Model A. Three free Hill parameters and nothing held fixed -- it has no "
              "microscopic rates to be non-identifiable about, which is exactly why it is "
              "the model to beat rather than the model to prefer."),
    ),
    "kinetic_jw95": ModelSpec(
        name="kinetic_jw95",
        param_names=("log10_kd", "log10_E", "log10_D"),
        bounds=((-1.0, 4.0), (-2.0, 3.0), (-3.0, 3.0)),
        factory=_model_b,
        conventions=("koff", "alpha", "r"),
        note=("Model B. The three equilibrium-identifiable ratios (P3); the absolute rates "
              "follow from the conventions and are NOT fitted."),
    ),
    "extended_desens": ModelSpec(
        name="extended_desens",
        param_names=("log10_kd", "log10_E", "log10_D_fast", "log10_D_slow"),
        bounds=((-1.0, 4.0), (-2.0, 3.0), (-3.0, 3.0), (-4.0, 3.0)),
        factory=_model_c,
        conventions=("koff", "alpha", "r_fast", "r_slow"),
        note=("Model C. NESTS MODEL B: at D_slow -> 0 the sixth state empties and the "
              "schemes coincide, so C cannot fit worse than B at its own optimum and the "
              "whole question is whether the improvement pays for the extra parameter."),
    ),
}


@dataclass(frozen=True)
class ModelFit:
    """One model at its MLE, with everything an information criterion needs."""

    name: str
    values: dict[str, float]
    log_likelihood: float
    #: Free parameters. What the model was allowed to estimate.
    n_free: int
    #: Estimated noise scales, one per dataset, profiled out of the likelihood.
    n_sigma: int
    n_data: int
    aic: float
    aicc: float
    bic: float
    #: Parameters that landed on a bound. These were NOT estimated -- the box stopped them
    #: -- so `k_effective` discounts them, and the discount is reported rather than quiet.
    at_bound: tuple[str, ...]
    converged_uniquely: bool
    note: str

    @property
    def k(self) -> int:
        """Parameters an information criterion is charged for."""
        return self.n_free + self.n_sigma

    @property
    def k_effective(self) -> int:
        """`k` less the parameters that sat on a bound."""
        return self.k - len(self.at_bound)


def _information_criteria(log_lik: float, k: int, n: int) -> tuple[float, float, float]:
    aic = 2.0 * k - 2.0 * log_lik
    bic = k * np.log(max(n, 1)) - 2.0 * log_lik
    # AICc is undefined at n <= k + 1 and returns +inf rather than a negative correction
    # that would make an over-parameterised model look good.
    aicc = aic + (2.0 * k * (k + 1) / (n - k - 1) if n > k + 1 else np.inf)
    return float(aic), float(aicc), float(bic)


def fit_one(spec: ModelSpec, datasets, *, n_starts: int = 8, seed: int = 0) -> ModelFit:
    """Fit one model to its MLE and compute its information criteria."""
    datasets = list(datasets)
    fit = fit_mle_vector(spec.factory, datasets, param_names=spec.param_names,
                         bounds=spec.bounds, n_starts=n_starts, seed=seed,
                         caller=f"compare_models/{spec.name}")
    n_data = int(sum(f.n for f in fit.per_dataset))
    n_sigma = len(fit.per_dataset)
    k = spec.n_free + n_sigma
    aic, aicc, bic = _information_criteria(fit.log_likelihood, k, n_data)
    at_bound = tuple(n for n, p in fit.at_bound.items() if any(p))
    return ModelFit(name=spec.name, values=fit.values,
                    log_likelihood=fit.log_likelihood, n_free=spec.n_free,
                    n_sigma=n_sigma, n_data=n_data, aic=aic, aicc=aicc, bic=bic,
                    at_bound=at_bound, converged_uniquely=fit.converged_uniquely,
                    note=fit.note)


# ------------------------------------------------------------------ blocked CV
def _blocks(ds, n_folds: int) -> list[np.ndarray]:
    """Contiguous concentration blocks, with the normalising point never held out.

    POINT-WISE FOLDS LEAK. On a smooth concentration-response a point's neighbours are
    almost its answer, so a random fold scores interpolation and the most flexible model
    wins by construction. Contiguous blocks in log-concentration make each fold an
    extrapolation.

    AND THE NORMALISING CONCENTRATION STAYS IN TRAINING. For an `I/I_max` dataset the
    response column is DEFINED by the response at its highest concentration; a fold that
    held that point out would score the model against a scale nobody gave it, and the
    prediction would be normalised at a different concentration than the data was. So for
    a normalised dataset the top point is excluded from the held-out blocks and always
    trains. It is one of nine points on the Jahn curve and the alternative is a fold whose
    residual measures the fold, not the model.
    """
    concs = np.asarray(ds.concs_um, dtype=float)
    order = np.argsort(concs)
    if getattr(ds, "normalisation", "absolute") == "fraction_of_max":
        order = order[:-1]
    if n_folds < 2:
        raise ValueError(f"n_folds must be at least 2, got {n_folds}")
    if order.size < n_folds:
        raise ValueError(
            f"{ds.citation_label!r} has {order.size} foldable points for {n_folds} folds. "
            f"A fold with no points is not a fold; lower n_folds.")
    return [b for b in np.array_split(order, n_folds) if b.size]


def _subset(ds, idx: np.ndarray):
    """A dataset restricted to `idx`, keeping every provenance field."""
    idx = np.asarray(idx, dtype=int)
    return dataclasses.replace(
        ds, concs_um=np.asarray(ds.concs_um, dtype=float)[idx],
        mean_response=np.asarray(ds.mean_response, dtype=float)[idx],
        sem=np.asarray(ds.sem, dtype=float)[idx])


@dataclass(frozen=True)
class CVScore:
    """Blocked cross-validated held-out log predictive density, summed over folds."""

    name: str
    total_lpd: float
    per_fold: tuple[float, ...]
    n_folds: int
    n_failed_folds: int


def blocked_cv(spec: ModelSpec, datasets, *, n_folds: int = 4, n_starts: int = 4,
               seed: int = 0) -> CVScore:
    """Fit on all blocks but one, score the one. Higher `total_lpd` is better.

    The held-out score uses the SIGMA ESTIMATED ON THE TRAINING BLOCKS, not one re-estimated
    on the held-out points. Re-estimating it would let a model that predicts the held-out
    block badly declare the noise large and score as well as one that predicts it well --
    which is not cross-validation, it is grading your own exam.
    """
    datasets = list(datasets)
    if len(datasets) != 1:
        raise NotImplementedError(
            "blocked CV is implemented for a single concentration-response dataset. With "
            "several, the blocks have to be formed per dataset and the folds crossed, and "
            "doing that silently on this repository's one real curve would be untested "
            "machinery pretending to be a result.")
    ds = datasets[0]
    folds = _blocks(ds, n_folds)
    all_idx = np.arange(np.asarray(ds.concs_um).size)

    per_fold, failed = [], 0
    for held in folds:
        train_idx = np.setdiff1d(all_idx, held)
        train = _subset(ds, train_idx)
        test = _subset(ds, held)
        try:
            fit = fit_mle_vector(spec.factory, [train], param_names=spec.param_names,
                                 bounds=spec.bounds, n_starts=n_starts, seed=seed,
                                 caller=f"blocked_cv/{spec.name}")
            sigma = float(fit.per_dataset[0].sigma_hat)
            model = spec.factory(np.array([fit.values[n] for n in spec.param_names]))
            scored = make_log_likelihood([test], caller=f"blocked_cv/{spec.name}/score",
                                         n_params=spec.n_free)(model)
            # the held-out density at the TRAINING sigma, not the test one
            rss = float(scored.per_dataset[0].rss)
            n = int(scored.per_dataset[0].n)
            lpd = float(-0.5 * (rss / sigma ** 2 + n * np.log(2.0 * np.pi * sigma ** 2)))
        except Exception:
            failed += 1
            continue
        per_fold.append(lpd)

    if not per_fold:
        raise RuntimeError(
            f"every fold failed for {spec.name}. A cross-validation with no scored folds is "
            f"not a score of zero.")
    return CVScore(name=spec.name, total_lpd=float(np.sum(per_fold)),
                   per_fold=tuple(per_fold), n_folds=len(folds), n_failed_folds=failed)


# ------------------------------------------------------------------ the comparison
@dataclass(frozen=True)
class ComparisonResult:
    fits: dict[str, ModelFit]
    cv: dict[str, CVScore]
    holdout_log_likelihood: dict[str, float]
    akaike_weights: dict[str, Quantity]
    #: The model the CV ranking prefers, and the one AICc prefers.
    cv_winner: str
    cv_verdict: "CVVerdict"
    aicc_winner: str
    rankings_agree: bool
    verdict: Quantity
    note: str


def compare_models(models: Sequence[str], train, holdout=(), *, n_folds: int = 4,
                   n_starts: int = 8, cv_starts: int = 4, seed: int = 0,
                   min_nats: float = 2.0,
                   observable: Observable | None = None) -> ComparisonResult:
    """Rank models by blocked CV, with AICc as secondary evidence and a stated verdict.

    `models` are keys into `MODEL_SPECS`. `train` and `holdout` are datasets; the holdout is
    scored once at each model's training MLE and never used to choose.

    THE VERDICT IS VOID UNLESS THE TWO RANKINGS AGREE. Cross-validation and an information
    criterion answer slightly different questions, and when they disagree the honest report
    is that the data do not discriminate -- not whichever number favours the conclusion.
    """
    train = list(train)
    holdout = list(holdout)
    specs = [MODEL_SPECS[m] for m in models]
    if len(specs) < 2:
        raise ValueError("a comparison needs at least two models")

    holdout_guard(train, caller="compare_models(train=...)")
    data_obs = _require_one_observable(specs, train, observable)
    synthetic = assert_real_data(train + holdout, caller="compare_models")

    fits = {s.name: fit_one(s, train, n_starts=n_starts, seed=seed) for s in specs}
    cv = {s.name: blocked_cv(s, train, n_folds=n_folds, n_starts=cv_starts, seed=seed)
          for s in specs}

    held = {}
    for s in specs:
        if not holdout:
            continue
        model = s.factory(np.array([fits[s.name].values[n] for n in s.param_names]))
        # `score_holdout`, not `make_log_likelihood`: the holdout guard refuses held-out
        # data to anything that could fit on it, and scoring an ALREADY-FITTED model is the
        # one thing that is allowed. The inverted guard there refuses training data, so
        # neither call can be made with the wrong dataset.
        held[s.name] = score_holdout(
            model, holdout, caller=f"compare_models/{s.name}").total_log_likelihood

    cv_verdict = resolve_cv(cv, fits, min_nats=min_nats)
    cv_winner = cv_verdict.winner
    aicc_winner = min(fits, key=lambda n: fits[n].aicc)
    agree = cv_winner == aicc_winner

    best_aicc = min(f.aicc for f in fits.values())
    raw = {n: float(np.exp(-0.5 * (f.aicc - best_aicc))) for n, f in fits.items()}
    total = sum(raw.values())

    cv_margin = cv_verdict.margin_nats
    caveats = (_caveats(fits, synthetic, cv) + _cross_native_note(specs, data_obs)
               + (cv_verdict.note,))

    weights = {}
    for n in fits:
        value = raw[n] / total if total > 0 else float("nan")
        if agree:
            weights[n] = Quantity(
                name=f"akaike_weight_{n}", _value=value, tier=Tier.UNCALIBRATED,
                provenance=(f"AICc weight at each model's own MLE, k = estimated "
                            f"identifiable parameters plus {fits[n].n_sigma} profiled noise "
                            f"scale(s); reported because the blocked-CV ranking agrees"),
                promote_by=("fit against digitised data -- see fitting.data."
                            "MISSING_DATASETS"),
                caveats=caveats)
        else:
            weights[n] = void(
                f"akaike_weight_{n}", value,
                defect=(f"the blocked-CV ranking prefers {cv_winner} and AICc prefers "
                        f"{aicc_winner}. A weight quoted from the criterion that happens to "
                        f"agree with a conclusion is the conclusion, not evidence"),
                register_item="P5", promote_by=(
                    "acquire data that discriminates -- see protocols/oed.py for the "
                    "protocol that would -- or report the comparison as unresolved"),
                caveats=caveats)

    if not agree:
        verdict = void(
            "model_comparison_verdict",
            f"UNRESOLVED: CV prefers {cv_winner}, AICc prefers {aicc_winner}",
            defect=("the two rankings disagree, so the data do not discriminate these "
                    "models by this comparison"),
            register_item="P5", promote_by="as above", caveats=caveats)
        note = (f"UNRESOLVED. Blocked CV prefers {cv_winner} ({cv_verdict.note}); AICc "
                f"prefers {aicc_winner}. Both are reported; neither is the answer.")
    else:
        verdict = Quantity(
            name="model_comparison_verdict", _value=cv_winner, tier=Tier.UNCALIBRATED,
            provenance=(f"blocked {n_folds}-fold CV and AICc both prefer {cv_winner}; CV "
                        f"margin over the runner-up {cv_margin:.3f} nats, paired SE "
                        f"{cv_verdict.margin_se:.3f}"),
            promote_by="fit against digitised data, and score the prospective holdout",
            caveats=caveats)
        note = (f"Blocked CV and AICc agree on {cv_winner}. {cv_verdict.note}.")

    return ComparisonResult(fits=fits, cv=cv, holdout_log_likelihood=held,
                            akaike_weights=weights, cv_winner=cv_winner,
                            cv_verdict=cv_verdict,
                            aicc_winner=aicc_winner, rankings_agree=agree,
                            verdict=verdict, note=note)


@dataclass(frozen=True)
class CVVerdict:
    """Who wins the cross-validation, and whether the margin is big enough to say so."""

    winner: str
    #: Every model whose held-out score is not distinguishable from the leader's.
    tied: tuple[str, ...]
    #: Leader minus runner-up, in nats of held-out log predictive density.
    margin_nats: float
    #: Standard error of that difference, from the PAIRED per-fold differences.
    margin_se: float
    note: str


def resolve_cv(cv: dict[str, CVScore], fits: dict[str, ModelFit], *,
               min_nats: float = 2.0) -> CVVerdict:
    """Name a CV winner only if it beats the field, and prefer parsimony in a tie.

    TWO THRESHOLDS, BOTH NECESSARY:

      * THE PAIRED STANDARD ERROR. Each fold scores every model, so the difference is a
        paired quantity and its uncertainty is the fold-to-fold scatter of the DIFFERENCE,
        which is far smaller than the scatter of either score. Comparing totals without it
        is how a 0.8-nat gap becomes a ranking.
      * A FLOOR OF 2 NATS. The roadmap's own forbidden list says a 2-point AIC difference
        is not decisive, and a held-out log-density difference below ~2 nats is a likelihood
        ratio under 8 on 24 points. Measured on the P5 recovery check: Model C, which NESTS
        Model B, scored 0.84 nats better out of sample with per-fold differences of
        [0.0, 0.7, 0.0, 0.1]. Calling that a win would have selected the larger model over
        the one the data was generated from, on nothing.

    IN A TIE, FEWEST ESTIMATED PARAMETERS WINS. This is the one-standard-error rule, and it
    is the only tie-break that does not quietly reward flexibility: a model that nests
    another can never score much worse, so "best total" alone always drifts upward in
    complexity.
    """
    ordered = sorted(cv.values(), key=lambda c: -c.total_lpd)
    best = ordered[0]
    tied = [best.name]
    margin = margin_se = 0.0
    reasons = []
    for i, other in enumerate(ordered[1:]):
        a = np.asarray(best.per_fold, dtype=float)
        b = np.asarray(other.per_fold, dtype=float)
        if a.size == b.size and a.size > 1:
            d = a - b
            diff = float(d.sum())
            se = float(np.std(d, ddof=1) * np.sqrt(d.size))
        else:
            # folds failed for one of them, so the comparison is not paired and the only
            # honest standard error is none. Said, not assumed.
            diff = float(best.total_lpd - other.total_lpd)
            se = float("nan")
            reasons.append(f"{other.name} was scored on a different set of folds, so its "
                           f"comparison with {best.name} is unpaired and has no stated "
                           f"standard error")
        if i == 0:
            margin, margin_se = diff, se
        threshold = min_nats if not np.isfinite(se) else max(min_nats, se)
        if diff < threshold:
            tied.append(other.name)

    winner = min(tied, key=lambda n: (fits[n].k, -cv[n].total_lpd))
    if len(tied) == 1:
        note = (f"{winner} leads by {margin:.3f} nats (paired SE {margin_se:.3f}), which "
                f"clears both the standard error and the {min_nats:g}-nat floor")
    else:
        note = (f"{', '.join(tied)} are not distinguishable out of sample (leader by "
                f"{margin:.3f} nats, paired SE {margin_se:.3f}, floor {min_nats:g}); "
                f"{winner} is selected as the one with the fewest estimated parameters "
                f"(k = {fits[winner].k})")
    if reasons:
        note += ". " + "; ".join(reasons)
    return CVVerdict(winner=winner, tied=tuple(tied), margin_nats=margin,
                     margin_se=margin_se, note=note)


def _caveats(fits, synthetic, cv) -> tuple[str, ...]:
    out = [
        "the absolute rates are NOT fitted: koff, alpha and the resensitisation rates are "
        "held at the conventions in comparison.FIT_CONVENTIONS, two of which have no "
        "recorded origin, so only the ratios are compared",
    ]
    if synthetic:
        out.append("fitted against parametric or synthetic data (" +
                   "; ".join(synthetic) + "), so this ranks models against a formula")
    pinned = {n: f.at_bound for n, f in fits.items() if f.at_bound}
    if pinned:
        out.append(f"parameters on a bound, which were not estimated: {pinned}")
    ties = [n for n, f in fits.items() if not f.converged_uniquely]
    if ties:
        out.append(f"the MLE is not unique for {ties}: the top five optima disagree, so "
                   f"its likelihood is one of several equally good ones")
    failed = {n: c.n_failed_folds for n, c in cv.items() if c.n_failed_folds}
    if failed:
        out.append(f"CV folds that failed to fit and are not scored: {failed}")
    return tuple(out)


def _require_one_observable(specs, datasets, observable) -> Observable:
    """ONE comparison axis, and every model must be able to produce it (P1 contract).

    NOT "every model must be NATIVE to it". Models B and C are EQUILIBRIUM-native and
    Model A is PEAK-native, and PEAK is nevertheless the one axis all three can produce --
    for B and C by integrating a square application, for A because its Hill curve IS a peak
    curve. Comparing them there is legitimate and is what the roadmap asks for. What is not
    legitimate, and what raised the chi-squared of 9290 this contract exists to prevent, is
    comparing Model A's PEAK against Model B's EQUILIBRIUM on one axis -- so the datasets
    must span a single observable, and a model that cannot speak it is named rather than
    silently given its own.
    """
    obs = {getattr(d, "observable", None) for d in datasets}
    if len(obs) > 1:
        raise ObservableMismatch(
            f"the training datasets span observables {sorted(o.value for o in obs)}. A "
            f"single comparison axis is required; fit per observable and say so, rather "
            f"than summing likelihoods across quantities that are not the same quantity.")
    data_obs = next(iter(obs))
    if observable is not None and observable is not data_obs:
        raise ObservableMismatch(
            f"asked to compare on {observable.value} against data tagged "
            f"{data_obs.value}.")

    # Probed through `_predict`, the same path the likelihood uses, so the check cannot
    # disagree with what the fit will actually do.
    probe = dataclasses.replace(
        datasets[0], concs_um=np.array([1.0, 10.0, 100.0]),
        mean_response=np.zeros(3), sem=np.full(3, 1e-3))
    cannot = {}
    for spec in specs:
        model = spec.factory(np.array([0.5 * (lo + hi) for lo, hi in spec.bounds]))
        try:
            _predict(model, probe)
        except ObservableMismatch as exc:
            cannot[spec.name] = str(exc).split(".")[0]
    if cannot:
        raise ObservableMismatch(
            f"the data is {data_obs.value} and these models cannot produce it: {cannot}. "
            f"Compare on PEAK, which all three can produce, or drop the model that has no "
            f"such quantity -- do not score it on a different one and put the numbers in "
            f"the same table.")
    return data_obs


def _cross_native_note(specs, data_obs: Observable) -> tuple[str, ...]:
    """Say so when a model is being scored off its own native observable."""
    off = []
    for spec in specs:
        model = spec.factory(np.array([0.5 * (lo + hi) for lo, hi in spec.bounds]))
        if model.native_observable is not data_obs:
            off.append(f"{spec.name} ({model.native_observable.value}-native)")
    if not off:
        return ()
    return (f"compared on {data_obs.value}, which is not the native observable of "
            + ", ".join(off) + ". That axis is legitimate -- every model here predicts it "
            "-- but a model is neither credited nor penalised for a quantity it does not "
            "have, and a PEAK comparison cannot see a difference that lives only at "
            "equilibrium",)


def comparison_report(result: ComparisonResult) -> ResultSet:
    """The comparison as a tiered result set."""
    rs = ResultSet(f"model comparison ({result.note})")
    rs.add(result.verdict)
    for name in result.fits:
        rs.add(result.akaike_weights[name])
    return rs


# ------------------------------------------------------------------ the method check
@dataclass(frozen=True)
class RecoveryCheck:
    truth_model: str
    selected: str
    rankings_agree: bool
    cv: dict[str, float]
    aicc: dict[str, float]
    recovered: bool
    note: str


def recovery_check(truth: str, *, noise_sd: float = 0.005, n_points: int = 24,
                   seed: int = 0, models=("operational", "kinetic_jw95",
                                          "extended_desens"),
                   n_folds: int = 4, n_starts: int = 6) -> RecoveryCheck:
    """Generate data from `truth` at its own MLE-free defaults and see what gets selected.

    A METHOD CHECK AND NOTHING ELSE. If the protocol cannot recover a model it generated
    the data from, it cannot adjudicate a real comparison -- and if it can, that says
    nothing whatever about receptors. The data is labelled synthetic and taints everything
    it touches, as it should.

    `truth` is generated at the midpoint of its own fitting box, so the truth is inside
    every model's search space and a failure to recover it cannot be blamed on the bounds.
    """
    from .data import peak_crc_from_model

    spec = MODEL_SPECS[truth]
    centre = np.array([0.5 * (lo + hi) for lo, hi in spec.bounds])
    truth_model = spec.factory(centre)
    # PEAK, because Model A has no equilibrium and refuses to produce one. See
    # `data.peak_crc_from_model`.
    ds = peak_crc_from_model(
        truth_model, concs_um=np.logspace(-1.0, 3.5, n_points), noise_sd=noise_sd,
        seed=seed, label=f"recovery check, truth = {truth}")

    res = compare_models(list(models), [ds], n_folds=n_folds, n_starts=n_starts,
                         cv_starts=max(2, n_starts // 2), seed=seed)
    recovered = res.rankings_agree and res.cv_winner == truth
    return RecoveryCheck(
        truth_model=truth, selected=res.cv_winner, rankings_agree=res.rankings_agree,
        cv={n: c.total_lpd for n, c in res.cv.items()},
        aicc={n: f.aicc for n, f in res.fits.items()},
        recovered=recovered,
        note=(f"truth {truth}; CV chose {res.cv_winner}, AICc chose {res.aicc_winner}. "
              + ("recovered" if recovered else "NOT recovered -- " + res.note)))
