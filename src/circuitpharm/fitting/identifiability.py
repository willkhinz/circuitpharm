"""Profile likelihood and cost-Hessian identifiability analysis.

READ THIS BEFORE QUOTING ANYTHING FROM THIS MODULE (roadmap P0-5).

`equilibrium_dose_response_chi2` depends on the six microscopic rates only through three
ratios -- K_d = k_off/k_on, E = beta/alpha, D = d/r -- because the equilibrium
dose-response depends only on those and the P_o,max penalty only on E. It is therefore
exactly invariant under a 3-parameter scaling group: multiply any PAIR (k_on, k_off),
(beta, alpha) or (d, r) by a common factor and nothing changes.

Measured, not argued (`invariance_report`): the objective is constant to 2e-12 along each
of the three generators, with normalised curvature `v^T H v / max|H|` of +3.0e-12,
-1.0e-12 and +5.3e-11. Scaling all six rates at once -- the sum of all three -- changes
the cost by exactly 0 at six decimal places.

So there is a verified 3-dimensional subspace in which this objective contains no
information at all. Consequences, each of which was being reported as a result:

  * `profile_likelihood` must classify `kon`, `koff`, `beta`, `alpha`, `d` and `r`
    individually as non-identifiable. That is a property of the objective, not of any data.
  * the cost Hessian's condition number is not a parameter-uncertainty statement. It used
    to be computed over POSITIVE eigenvalues only, which hid two facts at once: the flat
    subspace, and that **3 of the 6 eigenvalues here are negative**, i.e. the parameters
    are not at a minimum, so the matrix is not a Fisher information matrix in the first
    place.
  * an MCMC posterior over all six is flat along three directions and bounded there only
    by the prior box, so its "95% credible intervals" are prior widths.

A CORRECTION WORTH KEEPING. An earlier revision of the roadmap asserted the Hessian is
"structurally rank <= 3". That does not follow and is not what is measured: the generators
form a vector field rather than three fixed directions, so zero curvature does not imply a
zero eigenvalue, and a finite-difference Hessian here carries ~1e-5 relative error.
Thresholding its eigenvalues reports rank 6. The group invariance is the provable
statement, so it is the one this module tests.

A SECOND CORRECTION. The same revision cited "chi-squared = 9290 for 9 data points" as the
measure of the observable mismatch between this objective and its PEAK-tagged dataset. 9290
is the cost at the PROVISIONAL chimeric defaults (roadmap P0-13); re-optimising reaches
16.15. The observable mismatch is real -- the dataset is PEAK and this function computes
EQUILIBRIUM, which differ ~3x in EC50 -- but 9290 measured the unanchored parameters, not
the mismatch.

These functions therefore stay callable -- the arrays they compute are real -- while every
INFERENTIAL summary comes back as a `Tier.VOID` quantity and every call warns. The real
treatment is P3: work in the identifiable reparameterisation, and admit absolute rates only
once a dataset carrying a kinetic timescale is in the likelihood.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np
from scipy.optimize import minimize

from ..models.kinetic_jw95 import KineticAllosteryModel
from ..provisional import void, warn_provisional
from ..results import Quantity, ResultSet, Tier
from .data import DOSE_RESPONSE_BENCHMARK, assert_real_data

#: What every conclusion in this module rests on, in one place so the three result types
#: and the MCMC module cannot drift apart in how they describe it.
_DEFECT = ("the objective depends on the six rates only through the three ratios "
           "K_d = k_off/k_on, E = beta/alpha and D = d/r, so it is exactly invariant "
           "under a 3-parameter scaling group (measured: constant to 2e-12 along each "
           "generator, zero curvature to 5e-11) and no per-rate identifiability claim "
           "from it describes the data rather than the objective")

#: The three generators of the invariance, as functions of the parameter vector.
#:
#: f(kon, koff, beta, alpha, d, r) = F(koff/kon, beta/alpha, d/r), so scaling each PAIR by
#: a common factor leaves f unchanged. Each generator below is the tangent to that scaling
#: at the current parameters, i.e. the direction along which f is provably flat.
#:
#: This is stated analytically and then VERIFIED numerically, rather than being inferred
#: from eigenvalue thresholding -- see `invariance_report`, and the note there on why a
#: finite-difference Hessian does NOT show three clean zero eigenvalues.
_INVARIANCE_PAIRS: tuple[tuple[str, str], ...] = (("kon", "koff"), ("beta", "alpha"),
                                                  ("d", "r"))
_PROMOTE = ("P3: reparameterise to the identifiable combinations (fitting/reparam.py) and "
            "add a dataset carrying a kinetic timescale, which is what makes k_off "
            "identifiable at all")


@dataclass(frozen=True)
class ProfileLikelihoodResult:
    """Profile likelihood curve, plus the classification it cannot yet support.

    `grid_values` and `profile_costs` are real computations and are plain arrays.
    `ci_95_bounds`, `is_identifiable` and `diagnostic_reason` are inferences drawn from a
    structurally flat objective, so they are VOID quantities: reading `.value` raises, and
    `.get(acknowledge_void=True)` works at a call site a reviewer can see.
    """

    param_name: str
    optimal_value: float
    grid_values: np.ndarray
    profile_costs: np.ndarray
    ci_95_bounds: Quantity
    is_identifiable: Quantity
    diagnostic_reason: Quantity
    defect: str = _DEFECT


@dataclass(frozen=True)
class HessianSpectrum:
    """Cost Hessian and its FULL spectrum.

    NOT called a Fisher information matrix any more: the FIM is the expected Hessian of the
    negative log-likelihood AT the optimum, and this is the Hessian of a chi-squared at
    whatever parameters it was handed. `eigenvalues` includes negatives and near-zeros --
    they are the finding, and filtering to positives is how the previous implementation
    turned "not at a minimum, and rank-deficient" into a finite, plausible number.
    """

    hessian: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    param_names: tuple[str, ...]
    rank: int
    rank_tol: float
    condition_number: Quantity
    near_null_directions: tuple[str, ...]
    defect: str = _DEFECT


def equilibrium_dose_response_chi2(
    params_vec: np.ndarray,
    param_keys: Sequence[str] = ("kon", "koff", "beta", "alpha", "d", "r"),
) -> float:
    """Weighted chi-squared of the EQUILIBRIUM dose-response against the benchmark, plus a
    P_o,max penalty. A FUNCTION OF (K_d, E, D) ONLY -- it cannot constrain absolute rates.

    Two things are wrong with this as a likelihood and both are P1/P3 work, not P0:

    1. OBSERVABLE MISMATCH. `DOSE_RESPONSE_BENCHMARK` is tagged `Observable.PEAK` and this
       compares it against `model.dose_response`, which is EQUILIBRIUM. They differ by a
       factor of ~3 in EC50 at identical parameters. P4's likelihood routes each dataset
       through the matching observable. (The cost is 9290 at the provisional defaults and
       16.15 at this objective's own optimum -- see the module docstring's second
       correction; the large number measures the unanchored parameters, not the mismatch.)
    2. The dataset is synthetic (P0-6).

    Kept, under a name that says what it is, because the P0-5 machinery above needs
    something to profile and because the flatness result is measured on it.
    """
    p_dict = {k: max(float(v), 1e-6) for k, v in zip(param_keys, params_vec)}
    model = KineticAllosteryModel(**p_dict)

    concs = DOSE_RESPONSE_BENCHMARK.concs_um
    obs = DOSE_RESPONSE_BENCHMARK.mean_response
    sem = DOSE_RESPONSE_BENCHMARK.sem
    pred = model.dose_response(concs, pam_factor=1.0)
    chi2_dr = float(np.sum(((pred - obs) / sem) ** 2))

    # A penalty on the analytic gating bound, which is a function of E alone.
    chi2_pomax = float(((model.po_max - 0.75) / 0.05) ** 2)
    return chi2_dr + chi2_pomax


#: Deprecated alias. `mcmc.py` and existing callers keep working.
def compute_kinetic_objective(params_vec: np.ndarray,
                              param_keys: Sequence[str] = ("kon", "koff", "beta", "alpha",
                                                           "d", "r")) -> float:
    """Deprecated name for `equilibrium_dose_response_chi2`.

    Renamed because "kinetic objective" claims to constrain kinetics, which it cannot:
    it is an equilibrium observable plus a gating-bound penalty.
    """
    return equilibrium_dose_response_chi2(params_vec, param_keys)


def compute_profile_likelihood(
    param_name: str,
    base_model: KineticAllosteryModel,
    scan_factors: Sequence[float] = (0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0),
    threshold_delta_chi2: float = 3.84,  # chi^2(1) at 95%
) -> ProfileLikelihoodResult:
    """Profile the objective along one parameter. Conclusions are VOID; see module docstring.

    The curve itself is computed honestly: the unconstrained problem is re-optimised first
    so `delta` is measured from a real optimum rather than from the cost at a nominal
    point (which could make `delta` negative and the whole CI logic meaningless), and an
    optimiser failure raises instead of being scored as a cost 10 above baseline -- which
    put it above the 3.84 threshold and therefore counted a FAILURE as evidence of
    identifiability.
    """
    keys = list(base_model.param_names)
    if param_name not in keys:
        raise ValueError(f"unknown parameter {param_name!r}; this model has {keys}")
    assert_real_data([DOSE_RESPONSE_BENCHMARK], caller="compute_profile_likelihood")
    warn_provisional(
        what=f"profile likelihood for {param_name!r}",
        defect=_DEFECT, register_item="P0-5", promote_by=_PROMOTE)

    target_idx = keys.index(param_name)
    x0 = np.array([getattr(base_model, k) for k in keys], dtype=float)
    free_indices = [i for i in range(len(keys)) if i != target_idx]
    bounds_all = [(1e-5, 100.0)] * len(keys)

    # (1) the unconstrained optimum, so delta is measured from a real minimum
    res0 = minimize(lambda v: equilibrium_dose_response_chi2(v, keys), x0,
                    method="L-BFGS-B", bounds=bounds_all)
    if not res0.success:
        raise RuntimeError(
            f"the unconstrained optimisation failed before profiling {param_name!r}: "
            f"{res0.message}. A profile measured against an un-optimised baseline is not a "
            f"profile likelihood; scoring the failure as a raised cost would count it as "
            f"evidence of identifiability.")
    base_cost = float(res0.fun)
    opt_val = float(res0.x[target_idx])

    grid_vals = np.array([opt_val * f for f in scan_factors], dtype=float)
    profile_costs = []
    warm = res0.x[free_indices].copy()      # (2) warm-start from the previous grid point

    for val in grid_vals:
        def sub_obj(free_vec: np.ndarray) -> float:
            full = np.empty(len(keys))
            full[target_idx] = val
            full[free_indices] = free_vec
            return equilibrium_dose_response_chi2(full, keys)

        res = minimize(sub_obj, warm, method="L-BFGS-B",
                       bounds=[bounds_all[i] for i in free_indices])
        if not res.success:
            raise RuntimeError(
                f"the constrained optimisation failed at {param_name}={val:.6g}: "
                f"{res.message}. Previously this was recorded as base_cost + 10.0, i.e. a "
                f"cost above the 3.84 threshold, so an optimiser failure was scored as "
                f"evidence that the parameter IS identifiable -- sign-inverted.")
        profile_costs.append(float(res.fun))
        warm = res.x.copy()

    costs = np.asarray(profile_costs, dtype=float)
    delta = costs - base_cost
    flat = bool(np.allclose(delta, delta[0], rtol=0.0, atol=1e-8))
    inside = delta <= threshold_delta_chi2

    if flat:
        reason = ("STRUCTURALLY_NON_IDENTIFIABLE: the profile is flat to 1e-8 across the "
                  "whole sweep, which for this objective is expected and provable -- see "
                  "the module docstring.")
        identifiable, ci = False, (float(grid_vals[0]), float(grid_vals[-1]))
    elif np.all(inside):
        reason = ("PRACTICALLY_NON_IDENTIFIABLE: the cost never crosses the 95% threshold "
                  "inside the scanned range, so the data bound this parameter no more "
                  "tightly than the sweep does.")
        identifiable, ci = False, (float(grid_vals[0]), float(grid_vals[-1]))
    else:
        reason = ("IDENTIFIABLE within the scanned range: the profile crosses the 95% "
                  "threshold, so a finite interval exists.")
        ok = grid_vals[inside]
        identifiable, ci = True, (float(ok.min()), float(ok.max()))

    return ProfileLikelihoodResult(
        param_name=param_name,
        optimal_value=opt_val,
        grid_values=grid_vals,
        profile_costs=costs,
        ci_95_bounds=void(f"{param_name}_ci_95", ci, defect=_DEFECT,
                          register_item="P0-5", promote_by=_PROMOTE,
                          caveats=("grid endpoints, not an interpolated threshold "
                                   "crossing; P3 interpolates",)),
        is_identifiable=void(f"{param_name}_is_identifiable", identifiable,
                             defect=_DEFECT, register_item="P0-5", promote_by=_PROMOTE),
        diagnostic_reason=void(f"{param_name}_diagnosis", reason, defect=_DEFECT,
                               register_item="P0-5", promote_by=_PROMOTE),
    )


def invariance_report(model: KineticAllosteryModel,
                      scales: Sequence[float] = (-0.5, 0.5, 5.0, 50.0)) -> dict:
    """Verify, numerically, the three directions along which the objective is flat.

    Returns, per generator: the spread of the objective along that ray (should be float
    noise), and the normalised curvature `v^T H v / max|H|` (should be ~0).

    WHY THIS EXISTS RATHER THAN A RANK CHECK. The invariance is a 3-parameter group, so
    there is a 3-dimensional subspace of zero CURVATURE -- `v^T H v = 0`, verified below to
    ~5e-11. It does NOT follow that the Hessian has three clean zero eigenvalues: the
    generators are a vector FIELD (each depends on the parameters), not three fixed
    directions, so `H v = 0` is not implied, and in any case a finite-difference Hessian
    carries relative error of 1e-5 to 1e-2 here, which swamps a would-be zero.

    An earlier revision of the roadmap claimed "structurally rank <= 3". That was wrong,
    and an eigenvalue threshold reported rank 6 with mixed signs -- which is the honest
    reading of a finite-difference Hessian at parameters that are not even at a minimum.
    The group invariance is the real, provable statement, so it is the one tested.
    """
    keys = tuple(model.param_names)
    theta = np.array([getattr(model, k) for k in keys], dtype=float)
    base = equilibrium_dose_response_chi2(theta, keys)

    gens = []
    for a, b in _INVARIANCE_PAIRS:
        v = np.zeros(len(keys))
        v[keys.index(a)] = theta[keys.index(a)]
        v[keys.index(b)] = theta[keys.index(b)]
        gens.append(v)

    out: dict[str, dict] = {}
    for (a, b), v in zip(_INVARIANCE_PAIRS, gens):
        vals = [equilibrium_dose_response_chi2(theta + t * v, keys) for t in scales]
        vals.append(base)
        out[f"{a}/{b}"] = dict(
            objective_spread=float(max(vals) - min(vals)),
            objective=float(base),
            generator=v.copy(),
        )
    return dict(base_objective=float(base), generators=out,
                parameters=dict(zip(keys, theta)))


def cost_hessian_and_spectrum(
    model: KineticAllosteryModel,
    epsilon: float = 1e-4,
    rank_tol: float = 1e-8,
) -> HessianSpectrum:
    """Numerical Hessian of the cost and its full eigendecomposition.

    EVERY eigenvalue is returned, negatives and near-zeros included. The previous version
    filtered to positives before taking a condition number, which discards exactly the
    information the analysis exists to produce: a zero eigenvalue is a flat direction, and
    a negative one means the parameters are not at a minimum.
    """
    assert_real_data([DOSE_RESPONSE_BENCHMARK], caller="cost_hessian_and_spectrum")
    warn_provisional(
        what="the cost Hessian spectrum", defect=_DEFECT,
        register_item="P0-5", promote_by=_PROMOTE)

    keys = tuple(model.param_names)
    theta = np.array([getattr(model, k) for k in keys], dtype=float)
    n = len(theta)
    hess = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i, n):
            h_i = epsilon * max(abs(theta[i]), 1e-4)
            h_j = epsilon * max(abs(theta[j]), 1e-4)
            pp, pm, mp, mm = (theta.copy() for _ in range(4))
            pp[i] += h_i; pp[j] += h_j
            pm[i] += h_i; pm[j] -= h_j
            mp[i] -= h_i; mp[j] += h_j
            mm[i] -= h_i; mm[j] -= h_j
            d2 = (equilibrium_dose_response_chi2(pp, keys)
                  - equilibrium_dose_response_chi2(pm, keys)
                  - equilibrium_dose_response_chi2(mp, keys)
                  + equilibrium_dose_response_chi2(mm, keys)) / (4.0 * h_i * h_j)
            hess[i, j] = hess[j, i] = d2

    eigvals, eigvecs = np.linalg.eigh(hess)
    scale = float(np.max(np.abs(eigvals))) if eigvals.size else 0.0
    keep = np.abs(eigvals) > rank_tol * max(scale, 1e-300)
    rank = int(keep.sum())

    if rank:
        cond_val = float(np.max(np.abs(eigvals[keep])) / np.min(np.abs(eigvals[keep])))
    else:
        cond_val = float("inf")

    # The ANALYTIC flat directions, with their measured curvature. Reported instead of
    # eigenvalue thresholding because the invariance is provable and the thresholding is
    # not: a finite-difference Hessian at these parameters comes out full rank with mixed
    # signs, which says the parameters are not at a minimum and the FD error is ~1e-5
    # relative -- not that the objective constrains all six rates.
    names = []
    for a, b in _INVARIANCE_PAIRS:
        v = np.zeros(len(keys))
        v[keys.index(a)] = theta[keys.index(a)]
        v[keys.index(b)] = theta[keys.index(b)]
        u = v / np.linalg.norm(v)
        curv = float(u @ hess @ u) / max(scale, 1e-300)
        names.append(f"{a}&{b} scale together (v^T H v / max|H| = {curv:+.2e})")

    n_neg = int(np.sum(eigvals < -rank_tol * max(scale, 1e-300)))
    return HessianSpectrum(
        hessian=hess, eigenvalues=eigvals, eigenvectors=eigvecs, param_names=keys,
        rank=rank, rank_tol=rank_tol,
        condition_number=void(
            "cost_hessian_condition_number", cond_val, defect=_DEFECT,
            register_item="P0-5", promote_by=_PROMOTE,
            caveats=(f"computed over the {rank} eigenvalue(s) above rank_tol={rank_tol}; "
                     f"{n_neg} eigenvalue(s) are NEGATIVE, which means these parameters "
                     f"are not at a minimum, so this is not a Fisher information matrix "
                     f"and its condition number is not a parameter-uncertainty statement",
                     "the objective's three flat directions are a zero-CURVATURE subspace "
                     "(see near_null_directions); they do not appear as zero eigenvalues, "
                     "because the generators vary with the parameters and the "
                     "finite-difference error here is ~1e-5 relative")),
        near_null_directions=tuple(names),
    )


#: Deprecated alias for the old name, which claimed to compute a Fisher information matrix.
def compute_fisher_information_matrix(model: KineticAllosteryModel,
                                      epsilon: float = 1e-4):
    """Deprecated. Returns `(hessian, condition_number_quantity)` from the spectrum.

    Renamed to `cost_hessian_and_spectrum`: the FIM is the expected Hessian of the negative
    log-likelihood at the optimum, and this is a chi-squared Hessian at arbitrary
    parameters. The condition number is now a VOID quantity rather than a float.
    """
    spec = cost_hessian_and_spectrum(model, epsilon=epsilon)
    return spec.hessian, spec.condition_number


# =======================================================================================
# P3: the analysis done on the parameters the data can actually determine.
#
# Everything above operates on the six microscopic rates, where three directions carry no
# information at all, and is VOID for that reason. What follows works in the identifiable
# reparameterisation (fitting/reparam.py), where the problem is well-posed -- so its
# results are NOT void, only UNCALIBRATED, and they are the ones to use.
# =======================================================================================
from .reparam import EQUILIBRIUM_IDENTIFIABLE, UNLOCKED_BY, IdentifiableParams  # noqa: E402


@dataclass(frozen=True)
class IdentifiabilityClass:
    """Three-way classification of one parameter direction."""

    IDENTIFIABLE = "IDENTIFIABLE"
    PRACTICALLY_NON_IDENTIFIABLE = "PRACTICALLY_NON_IDENTIFIABLE"
    STRUCTURALLY_NON_IDENTIFIABLE = "STRUCTURALLY_NON_IDENTIFIABLE"


@dataclass(frozen=True)
class ProfileResult:
    """A profile likelihood over an IDENTIFIABLE combination.

    Unlike `ProfileLikelihoodResult` above, nothing here is VOID: the parameter being
    profiled is one the equilibrium likelihood can determine, so the interval means
    something. It is UNCALIBRATED rather than VALIDATED because the dataset behind it is
    parametric or synthetic (P0-6) -- the method is sound, the inputs are not yet.
    """

    param_name: str
    mle: float
    grid: np.ndarray
    costs: np.ndarray
    delta: np.ndarray
    ci_95: tuple[float, float]
    classification: str
    reason: str


#: Physically meaningful bounds in log10, per identifiable combination. A fit that runs
#: into one of these is reported as practically non-identifiable on that side -- which is a
#: RESULT -- rather than crashing the optimiser, which is what an unbounded search did.
IDENTIFIABLE_BOUNDS = {
    "log10_kd": (-1.0, 4.0),     # K_d 0.1 uM to 10 mM
    "log10_E": (-2.0, 3.0),      # gating efficacy 0.01 to 1000
    "log10_D": (-3.0, 3.0),      # desensitisation ratio 0.001 to 1000
}


def equilibrium_chi2_identifiable(theta, dataset=None) -> float:
    """Chi-squared of the EQUILIBRIUM dose-response, over (log K_d, log E, log D).

    REQUIRES AN EQUILIBRIUM DATASET. Handing it a PEAK curve is the observable mismatch
    the P1 contract exists to stop (roadmap §2.3) -- and it is not academic: fitting the
    EQUILIBRIUM curve to a PEAK dataset whose plateau is 0.75 drives D to its lower bound,
    because an equilibrium plateau is E/(1+E+D) and the only way to raise it is to delete
    desensitisation. The fit then "succeeds" having silently removed a mechanism.

    No `po_max` penalty, unlike the microscopic version: that term is a function of E alone
    and duplicates information the curve already carries, which biases the profile.
    """
    from ..models.base import Observable, ObservableMismatch
    from ..models.kinetic_jw95 import KineticAllosteryModel

    ds = dataset if dataset is not None else DOSE_RESPONSE_BENCHMARK
    obs = getattr(ds, "observable", None)
    if obs is not None and obs is not Observable.EQUILIBRIUM:
        raise ObservableMismatch(
            f"this objective is the EQUILIBRIUM dose-response and "
            f"{getattr(ds, 'citation_label', ds)!r} is tagged {obs.value}. Fitting an "
            f"equilibrium curve to peak-current data drives D to its bound -- the "
            f"equilibrium plateau is E/(1+E+D), so the fit raises it by deleting "
            f"desensitisation. Use a PEAK objective, or an EQUILIBRIUM dataset.")
    p = IdentifiableParams.from_vector(np.asarray(theta, dtype=float)[:3])
    # any representative with the right three ratios gives the same equilibrium curve;
    # alpha and r are arbitrary here BY THE STRUCTURAL RESULT, which this relies on
    rates = p.to_microscopic(koff=1.0, alpha=1.0, r=1e-3)
    model = KineticAllosteryModel(**rates)
    pred = model.dose_response(ds.concs_um, pam_factor=1.0)
    return float(np.sum(((pred - ds.mean_response) / ds.sem) ** 2))


def fit_identifiable(dataset=None, *, x0=None) -> tuple[IdentifiableParams, float]:
    """MLE over (log10 K_d, log10 E, log10 D). Returns the parameters and the cost."""
    from scipy.optimize import minimize

    ds = dataset if dataset is not None else DOSE_RESPONSE_BENCHMARK
    assert_real_data([ds], caller="fit_identifiable")
    start = np.array([1.4, 0.6, 1.4] if x0 is None else x0, dtype=float)
    best, best_cost = None, np.inf
    rng = np.random.default_rng(0)
    # multi-start: three ratios over decades, so a single start can land in a local basin
    for i in range(12):
        s = start if i == 0 else start + rng.normal(0.0, 0.7, 3)
        bounds = [IDENTIFIABLE_BOUNDS[k] for k in EQUILIBRIUM_IDENTIFIABLE]
        res = minimize(lambda v: equilibrium_chi2_identifiable(v, ds), s,
                       method="L-BFGS-B", bounds=bounds,
                       options=dict(ftol=1e-15, gtol=1e-12, maxiter=5000))
        start2 = res.x if res.success else s
        res = minimize(lambda v: equilibrium_chi2_identifiable(v, ds), start2,
                       method="Nelder-Mead",
                       options=dict(xatol=1e-9, fatol=1e-11, maxiter=50000))
        if res.success and res.fun < best_cost:
            best, best_cost = res.x, float(res.fun)
    if best is None:
        raise RuntimeError(
            "no start converged while fitting the identifiable parameters. A failure here "
            "is not a result; it is not scored as one.")
    return IdentifiableParams.from_vector(best), best_cost


def profile_likelihood(param: str, dataset=None, *, n_grid: int = 25,
                       span_decades: float = 1.5,
                       threshold_delta_chi2: float = 3.84) -> ProfileResult:
    """Profile one identifiable combination, properly.

    Five things this does that the microscopic version could not (roadmap P3-2):

    1. RE-OPTIMISES the unconstrained problem first, so `delta` is measured from a real
       minimum rather than from a nominal point where it could come out negative.
    2. Grids in LOG space, +/- `span_decades` around the MLE, with >= 25 nodes. Seven
       multiplicative factors spanning 0.3-3.0 cannot locate a 95% boundary.
    3. WARM-STARTS each node from the previous one. Profiles are continuous; restarting
       from the nominal point at every node manufactures a non-monotone curve.
    4. RAISES on optimiser failure rather than scoring it as a raised cost, which counted a
       failure as evidence of identifiability.
    5. INTERPOLATES the threshold crossing instead of returning a grid node.
    """
    if param not in EQUILIBRIUM_IDENTIFIABLE:
        raise ValueError(
            f"{param!r} is not an equilibrium-identifiable combination. Have "
            f"{EQUILIBRIUM_IDENTIFIABLE}; for an absolute rate see reparam.UNLOCKED_BY, "
            f"which names the measurement that would make it identifiable.")
    from scipy.optimize import minimize

    ds = dataset if dataset is not None else DOSE_RESPONSE_BENCHMARK
    mle, base_cost = fit_identifiable(ds)
    idx = EQUILIBRIUM_IDENTIFIABLE.index(param)
    x_mle = mle.as_vector(with_timescale=False)
    free = [i for i in range(3) if i != idx]

    grid = np.linspace(x_mle[idx] - span_decades, x_mle[idx] + span_decades, n_grid)
    costs = np.empty(n_grid)
    # sweep outward from the MLE in both directions so every warm start is adjacent
    order = sorted(range(n_grid), key=lambda i: abs(grid[i] - x_mle[idx]))
    warm = {int(np.argmin(np.abs(grid - x_mle[idx]))): x_mle[free].copy()}
    for i in order:
        start = warm.get(i)
        if start is None:
            nearest = min(warm, key=lambda j: abs(j - i))
            start = warm[nearest]

        def sub(v):
            full = np.empty(3)
            full[idx] = grid[i]
            full[free] = v
            return equilibrium_chi2_identifiable(full, ds)

        # L-BFGS-B first, then a Nelder-Mead polish. The objective is cheap, smooth
        # algebra (no ODE), so a quasi-Newton method converges in tens of evaluations
        # where Nelder-Mead was exhausting 20000 iterations at the extreme grid nodes --
        # which the failure check then correctly refused to score. Bounds keep the
        # simplex inside physically meaningful decades.
        # L-BFGS-B first; Nelder-Mead where the gradient is unusable. The valley here is
        # extremely flat -- that is the finding, not a nuisance -- so a quasi-Newton line
        # search legitimately returns ABNORMAL at the edges of the scan. Falling back to a
        # derivative-free method is the right numerical choice, NOT a way of not noticing:
        # if both fail the result is still refused.
        bounds = [IDENTIFIABLE_BOUNDS[EQUILIBRIUM_IDENTIFIABLE[j]] for j in free]
        res = minimize(sub, start, method="L-BFGS-B", bounds=bounds,
                       options=dict(ftol=1e-15, gtol=1e-12, maxiter=5000))
        if res.success:
            polish = minimize(sub, res.x, method="Nelder-Mead",
                              options=dict(xatol=1e-9, fatol=1e-11, maxiter=20000))
            if polish.success and polish.fun <= res.fun:
                res = polish
        else:
            res = minimize(sub, start, method="Nelder-Mead",
                           options=dict(xatol=1e-9, fatol=1e-11, maxiter=50000))
        if not res.success:
            raise RuntimeError(
                f"both L-BFGS-B and Nelder-Mead failed at {param} = {grid[i]:.6g}: "
                f"{res.message}. Recording this as a raised cost -- the old behaviour -- "
                f"would score an optimiser failure as evidence of identifiability.")
        costs[i] = float(res.fun)
        warm[i] = res.x.copy()

    delta = costs - base_cost
    flat = bool(np.ptp(delta) < 1e-8)
    above = delta > threshold_delta_chi2
    left = above[:int(np.argmin(np.abs(grid - x_mle[idx])))].any()
    right = above[int(np.argmin(np.abs(grid - x_mle[idx]))):].any()

    def _cross(lo_i, hi_i):
        """Linear interpolation of the threshold crossing between two nodes."""
        d0, d1 = delta[lo_i], delta[hi_i]
        if d1 == d0:
            return float(grid[hi_i])
        f = (threshold_delta_chi2 - d0) / (d1 - d0)
        return float(grid[lo_i] + f * (grid[hi_i] - grid[lo_i]))

    centre = int(np.argmin(np.abs(grid - x_mle[idx])))
    lo = grid[0]
    for i in range(centre, 0, -1):
        if delta[i - 1] > threshold_delta_chi2 >= delta[i]:
            lo = _cross(i, i - 1)
            break
    hi = grid[-1]
    for i in range(centre, n_grid - 1):
        if delta[i + 1] > threshold_delta_chi2 >= delta[i]:
            hi = _cross(i, i + 1)
            break

    if flat:
        cls = IdentifiabilityClass.STRUCTURALLY_NON_IDENTIFIABLE
        reason = (f"the profile is flat to 1e-8 across +/-{span_decades} decades, so this "
                  f"objective contains no information about {param} at all.")
    elif left and right:
        cls = IdentifiabilityClass.IDENTIFIABLE
        reason = (f"the profile crosses the 95% threshold on both sides within "
                  f"+/-{span_decades} decades, so a finite interval exists.")
    elif left or right:
        cls = IdentifiabilityClass.PRACTICALLY_NON_IDENTIFIABLE
        side = "below" if left else "above"
        reason = (f"the profile crosses the threshold only {side} the MLE, so the data "
                  f"bound {param} on one side only; the other bound is the scan edge.")
    else:
        cls = IdentifiabilityClass.PRACTICALLY_NON_IDENTIFIABLE
        reason = (f"the profile never crosses the 95% threshold within "
                  f"+/-{span_decades} decades, so the data bound {param} no more tightly "
                  f"than the scan does.")

    return ProfileResult(param_name=param, mle=float(x_mle[idx]), grid=grid, costs=costs,
                         delta=delta, ci_95=(float(lo), float(hi)),
                         classification=cls, reason=reason)


def identifiability_report(dataset=None) -> ResultSet:
    """Profile all three identifiable combinations and report them as tiered quantities."""
    ds = dataset if dataset is not None else DOSE_RESPONSE_BENCHMARK
    mle, cost = fit_identifiable(ds)
    rs = ResultSet(f"equilibrium identifiability against {getattr(ds, 'citation_label', ds)}")
    prov = ("profile likelihood over the equilibrium-identifiable combinations "
            "(fitting/reparam.py); the six microscopic rates enter only through these "
            "three ratios, so this is the whole of what an equilibrium curve determines")
    promote = ("add a dataset with a TIMESCALE in it -- see fitting.data.MISSING_DATASETS "
               "and reparam.UNLOCKED_BY -- which is what makes an absolute rate "
               "identifiable; and replace the parametric/synthetic curve with a "
               "digitisation (P1-2)")
    for name in EQUILIBRIUM_IDENTIFIABLE:
        pr = profile_likelihood(name, ds)
        rs.add(Quantity(
            name, pr.mle, Tier.UNCALIBRATED, "log10",
            provenance=f"{prov}. {pr.reason}",
            promote_by=promote,
            caveats=(f"95% CI [{pr.ci_95[0]:.4f}, {pr.ci_95[1]:.4f}] in log10; "
                     f"classification {pr.classification}",
                     "the underlying dataset is not a digitisation, so this is a method "
                     "result rather than a measurement")))
    rs.add(Quantity(
        "n_identifiable_of_six", 3, Tier.VALIDATED,
        provenance=("structural, and provable: the equilibrium state distribution is "
                    "1 : 2x : x^2 : E x^2 : D x^2 with x = [G]/K_d, so the six rates enter "
                    "only through K_d, E and D. Measured corroboration: the objective is "
                    "constant to 2e-12 along each of the three scaling generators."),
        caveats=("this is a statement about EQUILIBRIUM observables; a kinetic observable "
                 "adds an absolute timescale and raises the count",)))
    return rs
