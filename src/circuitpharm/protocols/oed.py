"""Statistical model comparison and optimal experimental design.

EVERYTHING THIS MODULE RETURNS IS CURRENTLY `Tier.VOID` (roadmap P0-5 disposition, C5, C6).
The machinery stays callable and the raw predictions are real; the model-selection and
falsification CONCLUSIONS are not available, for four reasons that are defects rather than
uncertainties:

  1. NO FIT. `evaluate_model_fit` scores each model at whatever parameters it was handed,
     at pam_factor 1.0. AIC and BIC are statements about a model AT ITS MLE; computed at
     arbitrary parameters they rank starting guesses. The three models' defaults are
     uncalibrated and mutually inconsistent (P0-13), so the ranking is close to arbitrary.
  2. k COUNTS DECLARED, NOT ESTIMATED, PARAMETERS. `len(m.param_names)` charges Model C for
     9 parameters including `pam_desens_factor`, which the pam_factor=1.0 equilibrium curve
     does not use at all, and charges Model B for `d` and `r`, which enter only through
     their ratio. The complexity penalty is therefore wrong in a direction that
     systematically disfavours the more detailed schemes.
  3. SIGMA IS ASSUMED, NOT ESTIMATED. The Gaussian log-likelihood takes a fixed
     `measurement_noise_std`, so every AIC difference and every Akaike weight scales with
     an unjustified knob: at sigma 0.03 the weights are near-uniform and at 0.003 they are
     a delta. A ranking that moves that far on a free parameter is not a ranking.
  4. THE DISCRIMINATION SCORE IGNORES PARAMETER UNCERTAINTY. It is a difference of point
     predictions divided by a fixed measurement noise, so it reports a protocol as
     discriminating when the models' own posteriors overlap completely. The intervals it
     emits as "95% confidence tolerance" are +/- 2 sigma_noise, which is narrower than the
     real predictive spread and therefore overconfident in exactly the direction that
     matters for a pre-registration.

Additionally, until P1's contract was added these comparisons ran ACROSS OBSERVABLES:
Model A is PEAK-native and Models B and C are EQUILIBRIUM-native, and the three were
compared on one axis. That now raises instead.

P5 is the fix: an independent MLE per model, k counting estimated identifiable parameters,
sigma estimated or profiled out, blocked cross-validation, and a prospective holdout. P6
replaces the discrimination score with a posterior-predictive one.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence
import numpy as np

from ..models.base import Observable, ObservableMismatch, ReceptorModel
from ..provisional import void, warn_provisional
from ..results import Quantity, ResultSet, Tier

_DEFECT = ("AIC/BIC/weights computed at UNFITTED parameters, with k counting declared "
           "rather than estimated parameters and sigma assumed rather than estimated -- so "
           "this ranks starting guesses, not models")
_PROMOTE = ("P5: fit each model independently to its MLE, count only estimated identifiable "
            "parameters, estimate or profile out sigma, and score by blocked "
            "cross-validation plus a prospective holdout")
_OED_DEFECT = ("the score is a difference of POINT predictions over a fixed measurement "
               "noise, so it ignores parameter uncertainty entirely and its intervals are "
               "narrower than the real predictive spread")
_OED_PROMOTE = ("P6: score by posterior-predictive separation, "
                "E|y_A - y_B| / sqrt(var(y_A) + var(y_B) + sigma^2), using P4's posterior")


@dataclass(frozen=True)
class ModelComparisonResult:
    """Model-selection metrics. The SELECTION fields are VOID; see the module docstring.

    `rss` and `log_likelihood` are plain -- they are real computations against the data
    handed in. `aic`, `bic` and `akaike_weight` are VOID quantities, because an information
    criterion at unfitted parameters with a declared-parameter count and an assumed sigma
    is not a model-selection result.
    """
    model_name: str
    n_params: int
    n_data: int
    rss: float
    log_likelihood: float
    aic: Quantity
    bic: Quantity
    akaike_weight: Quantity
    observable: Observable
    defect: str = _DEFECT


@dataclass(frozen=True)
class DiscriminatingProtocol:
    """A candidate protocol. The SCORE and the BOUNDARIES are VOID; see the module docstring.

    `gaba_um`, `pam_factor` and `predicted_responses` are plain -- they are the protocol and
    the models' point predictions at it. `discrimination_score` and
    `falsification_boundaries` are VOID: the score divides a difference of point
    predictions by a fixed measurement noise, so it cannot see parameter uncertainty, and
    the boundaries are +/- 2 sigma_noise, which is narrower than the real predictive spread.
    Pre-registering against them would be overconfident in precisely the direction that
    matters.
    """
    gaba_um: float
    pam_factor: float
    pulse_duration_ms: float
    predicted_responses: dict[str, float]
    observable: Observable
    falsification_boundaries: Quantity
    discrimination_score: Quantity
    defect: str = _OED_DEFECT


def evaluate_model_fit(
    models: Sequence[ReceptorModel],
    concs_um: np.ndarray,
    observed_responses: np.ndarray,
    measurement_noise_std: float = 0.03,
    observable: Observable | None = None,
) -> list[ModelComparisonResult]:
    """Calculate log-likelihood, AIC, BIC, and Akaike weights across candidate models.

    Parameters
    ----------
    models : Sequence[ReceptorModel]
        List of competing models (e.g. Model A, Model B, Model C).
    concs_um : np.ndarray
        Agonist concentrations of the test data.
    observed_responses : np.ndarray
        Measured open probabilities or normalized current amplitudes.
    measurement_noise_std : float
        Estimated standard deviation of measurement noise.
    """
    y_obs = np.asarray(observed_responses, dtype=float)
    concs = np.asarray(concs_um, dtype=float)
    n_data = len(y_obs)
    sigma = max(measurement_noise_std, 1e-6)

    # ONE OBSERVABLE, OR IT RAISES. Model A is PEAK-native and Models B and C are
    # EQUILIBRIUM-native; comparing them on one axis was comparing a peak current against a
    # desensitised stationary level, which differ ~3x in EC50 at identical parameters
    # (roadmap §2.3). `observable` says which axis the comparison is on, and every model
    # must be able to speak it.
    natives = {m.name: m.native_observable for m in models}
    if observable is None:
        distinct = set(natives.values())
        if len(distinct) > 1:
            raise ObservableMismatch(
                f"these models are native to different observables ({natives}), so there "
                f"is no default axis to compare them on. Pass `observable=` explicitly -- "
                f"Observable.PEAK is the one all three can produce.")
        observable = next(iter(distinct))

    warn_provisional(what="model comparison (AIC/BIC/Akaike weights)", defect=_DEFECT,
                     register_item="P0-5 / C5", promote_by=_PROMOTE)

    raw_results = []
    aics = []

    for m in models:
        if observable is Observable.PEAK:
            y_pred = m.peak_dose_response(concs, pam_factor=1.0)
        elif observable is Observable.EQUILIBRIUM:
            y_pred = m.dose_response(concs, pam_factor=1.0)
        else:
            raise ObservableMismatch(
                f"a concentration-response comparison is on PEAK or EQUILIBRIUM; "
                f"{observable.value} is a charge and needs a window and a waveform.")
        residuals = y_obs - y_pred
        rss = float(np.sum(residuals ** 2))
        
        # Gaussian log-likelihood
        log_lik = float(-0.5 * np.sum((residuals / sigma) ** 2) - n_data * np.log(sigma * np.sqrt(2.0 * np.pi)))
        k = len(m.param_names)
        aic = 2.0 * k - 2.0 * log_lik
        bic = k * np.log(max(n_data, 1)) - 2.0 * log_lik
        
        aics.append(aic)
        raw_results.append((m.name, k, rss, log_lik, aic, bic))

    # Akaike weights: exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
    min_aic = min(aics)
    delta_aics = np.array([a - min_aic for a in aics])
    w = np.exp(-0.5 * delta_aics)
    weights = w / np.sum(w)

    final_results = []
    for (name, k, rss, log_lik, aic, bic), weight in zip(raw_results, weights):
        def _v(label, value, extra=()):
            return void(f"{name}_{label}", value, defect=_DEFECT,
                        register_item="P0-5 / C5", promote_by=_PROMOTE,
                        caveats=(f"k = {k} counts DECLARED parameters; the estimated "
                                 f"identifiable count is what an information criterion "
                                 f"needs",
                                 f"sigma = {sigma:g} was assumed, not estimated, and every "
                                 f"difference here scales with it") + tuple(extra))
        final_results.append(ModelComparisonResult(
            model_name=name,
            n_params=k,
            n_data=n_data,
            rss=rss,
            log_likelihood=log_lik,
            aic=_v("aic", aic),
            bic=_v("bic", bic),
            akaike_weight=_v("akaike_weight", float(weight),
                             ("a weight computed from VOID AICs is VOID; report it only "
                              "alongside a cross-validated ranking that agrees (P5)",)),
            observable=observable,
        ))

    return final_results


def find_discriminating_protocol(
    models: Sequence[ReceptorModel],
    gaba_candidates_um: Sequence[float] | None = None,
    pam_candidates: Sequence[float] | None = None,
    noise_sigma: float = 0.02,
    observable: Observable | None = None,
) -> DiscriminatingProtocol:
    """Find the GABA and PAM exposure combination that maximizes separation between models.

    Parameters
    ----------
    models : Sequence[ReceptorModel]
        The candidate models to discriminate between.
    gaba_candidates_um : Sequence[float], optional
        Candidate GABA concentrations to test (defaults to 0.1 to 1000 uM sweep).
    pam_candidates : Sequence[float], optional
        Candidate PAM multipliers (defaults to 1.5, 2.5, 4.0).
    noise_sigma : float
        Typical experimental patch-clamp noise level.
    """
    gaba_grid = gaba_candidates_um or [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 1000.0]
    pam_grid = pam_candidates or [1.5, 2.0, 2.5, 3.0, 4.0]

    # ONE OBSERVABLE, as in evaluate_model_fit. `steady_state` is EQUILIBRIUM for Models B
    # and C and PEAK for Model A, so the old code was maximising a difference partly
    # manufactured by the axis meaning different things for different models.
    natives = {m.name: m.native_observable for m in models}
    if observable is None:
        distinct = set(natives.values())
        if len(distinct) > 1:
            raise ObservableMismatch(
                f"these models are native to different observables ({natives}); pass "
                f"`observable=Observable.PEAK`, which all three can produce. Maximising a "
                f"difference across observables finds an axis mismatch, not a mechanism.")
        observable = next(iter(distinct))

    warn_provisional(what="the discriminating-protocol search", defect=_OED_DEFECT,
                     register_item="C6 / P6", promote_by=_OED_PROMOTE)

    def _predict(m, g, p):
        if observable is Observable.PEAK:
            return float(np.atleast_1d(m.peak_dose_response(np.array([g]), pam_factor=p))[0])
        return float(m.steady_state(g, pam_factor=p))

    best_score = -1.0
    best_gaba = 1.0
    best_pam = 2.5
    best_preds: dict[str, float] = {}

    for g in gaba_grid:
        for p in pam_grid:
            preds = {m.name: _predict(m, g, p) for m in models}
            vals = list(preds.values())
            # Spread = max pairwise distance normalized by noise
            max_dist = max(vals) - min(vals)
            score = max_dist / max(noise_sigma, 1e-4)

            if score > best_score:
                best_score = score
                best_gaba = g
                best_pam = p
                best_preds = preds

    # +/- 2 sigma_noise around each point prediction. NOT a 95% predictive interval: it
    # omits parameter uncertainty, which for a posterior as flat as P0-5 describes is the
    # dominant term. Labelled VOID rather than renamed, because the number is also
    # computed from an unfitted model.
    boundaries = {
        name: (max(0.0, val - 2.0 * noise_sigma), min(1.0, val + 2.0 * noise_sigma))
        for name, val in best_preds.items()
    }

    return DiscriminatingProtocol(
        gaba_um=best_gaba,
        pam_factor=best_pam,
        pulse_duration_ms=100.0,
        predicted_responses=best_preds,
        observable=observable,
        falsification_boundaries=void(
            "falsification_boundaries", boundaries, defect=_OED_DEFECT,
            register_item="C6 / P6", promote_by=_OED_PROMOTE,
            caveats=("+/- 2 sigma_noise only: parameter uncertainty is omitted, so these "
                     "are narrower than the real predictive spread",
                     "and they are computed at unfitted parameters (C5)")),
        discrimination_score=void(
            "discrimination_score", float(best_score), defect=_OED_DEFECT,
            register_item="C6 / P6", promote_by=_OED_PROMOTE,
            caveats=("the protocol space was searched, so even a posterior-predictive "
                     "version of this score would owe a multiple-comparison correction",)),
    )
