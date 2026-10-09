"""Parameter estimation, empirical datasets, identifiability, and MCMC sampling.

EVERY INFERENTIAL CONCLUSION IN THIS SUBPACKAGE IS CURRENTLY `Tier.VOID`, and every call
emits `provisional.ProvisionalResultWarning`. Two reasons, both recorded in the P0
register: the objective is a function of (K_d, E, D) only, so three of six parameter
directions are structurally flat (P0-5); and both datasets are synthetic (P0-6). The
computations are real and the arrays are plain; the conclusions are not available without
`.get(acknowledge_void=True)`. P3 and P4 are the fix.
"""
from .data import (
    ALL,
    equilibrium_crc_from_model,
    DEACTIVATION_BENCHMARK,
    DOSE_RESPONSE_BENCHMARK,
    HOLDOUT,
    JAHN1997_PEAK_CRC,
    MISSING_DATASETS,
    TRAIN,
    DataKind,
    DeactivationDataset,
    DoseResponseDataset,
    assert_real_data,
    holdout_guard,
)
from .identifiability import (
    EQUILIBRIUM_IDENTIFIABLE,
    IDENTIFIABLE_BOUNDS,
    HessianSpectrum,
    IdentifiabilityClass,
    ProfileResult,
    ProfileLikelihoodResult,
    compute_fisher_information_matrix,
    compute_kinetic_objective,
    compute_profile_likelihood,
    cost_hessian_and_spectrum,
    equilibrium_chi2_identifiable,
    equilibrium_dose_response_chi2,
    fit_identifiable,
    identifiability_report,
    invariance_report,
    profile_likelihood,
)
from .reparam import UNLOCKED_BY, IdentifiableParams
from .mcmc import (
    LOG10_BOUNDS,
    MCMCChainResult,
    integrated_autocorr_time,
    log_posterior_log10,
    log_prior_log10,
    run_ensemble_mcmc,
    split_rhat,
)

__all__ = [
    # data
    "ALL",
    "DOSE_RESPONSE_BENCHMARK",
    "DEACTIVATION_BENCHMARK",
    "TRAIN",
    "HOLDOUT",
    "JAHN1997_PEAK_CRC",
    "MISSING_DATASETS",
    "DataKind",
    "DoseResponseDataset",
    "DeactivationDataset",
    "assert_real_data",
    "holdout_guard",
    "equilibrium_crc_from_model",
    # identifiability
    "equilibrium_dose_response_chi2",
    "compute_kinetic_objective",
    "ProfileLikelihoodResult",
    "compute_profile_likelihood",
    "HessianSpectrum",
    "cost_hessian_and_spectrum",
    "compute_fisher_information_matrix",
    "invariance_report",
    # P3: the analysis on the parameters the data can determine
    "EQUILIBRIUM_IDENTIFIABLE",
    "IDENTIFIABLE_BOUNDS",
    "IdentifiabilityClass",
    "IdentifiableParams",
    "ProfileResult",
    "UNLOCKED_BY",
    "equilibrium_chi2_identifiable",
    "fit_identifiable",
    "profile_likelihood",
    "identifiability_report",
    # mcmc
    "LOG10_BOUNDS",
    "MCMCChainResult",
    "run_ensemble_mcmc",
    "log_prior_log10",
    "log_posterior_log10",
    "integrated_autocorr_time",
    "split_rhat",
]
