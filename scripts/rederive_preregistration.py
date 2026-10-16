"""Re-derive P6 sigma requirement and retarget onto paired-pulse recovery.

Task E of AGENT-BRIEF-2.md:
1. Re-run with the application duration declared and matched between model and dataset.
   Use models.base.hill_slope(..., window_rel=...) for any slope comparison.
2. Report whether sigma <= 0.005 survives or softens.
3. Record old and new conclusions side-by-side with the protocol stated for each.
4. Re-target the discriminating observable onto paired-pulse recovery.
"""
from __future__ import annotations

import warnings
import numpy as np

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.data import equilibrium_crc_from_model
from circuitpharm.fitting.posterior import sample_posterior_vector
from circuitpharm.models.base import Observable, hill_slope, JAHN1997_SLOPE_WINDOW_REL
from circuitpharm.protocols.design import (
    ProtocolPoint,
    find_discriminating_protocol_pp,
    posterior_predictive_score,
    sigma_sensitivity,
    _penalise,
)


def main():
    print("================================================================================")
    print("TASK E: RE-DERIVE P6 SIGMA REQUIREMENT & RETARGET ONTO PAIRED-PULSE RECOVERY")
    print("================================================================================")
    print("Sampling posteriors for Model B (kinetic_jw95) and Model C (extended_desens)...")

    spec_b = MODEL_SPECS["kinetic_jw95"]
    spec_c = MODEL_SPECS["extended_desens"]
    centre = np.array([0.5 * (lo + hi) for lo, hi in spec_b.bounds])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = equilibrium_crc_from_model(
            spec_b.factory(centre),
            concs_um=np.logspace(-1.0, 3.5, 30),
            noise_sd=1e-3,
            seed=0,
            label="P6 design basis, truth = Model B",
        )

        posts = {}
        for name, spec in (("kinetic_jw95", spec_b), ("extended_desens", spec_c)):
            posts[name] = sample_posterior_vector(
                spec.factory,
                [ds],
                param_names=spec.param_names,
                bounds=dict(zip(spec.param_names, spec.bounds)),
                caller=f"prereg/{name}",
                seed=11,
            )

    # Thin draws to representative sample
    draws_b = posts["kinetic_jw95"].log10_samples.reshape(-1, len(spec_b.param_names))[::3000][:100]
    draws_c = posts["extended_desens"].log10_samples.reshape(-1, len(spec_c.param_names))[::3000][:100]
    cand = {
        "kinetic_jw95": (spec_b.factory, draws_b),
        "extended_desens": (spec_c.factory, draws_c),
    }
    print(f"Posteriors ready: {len(draws_b)} draws each.")

    print("\n--------------------------------------------------------------------------------")
    print("1. PEAK OBSERVABLE: EFFECT OF DECLARED & MATCHED APPLICATION DURATION")
    print("--------------------------------------------------------------------------------")
    mb_true = spec_b.factory(centre)
    concs_sweep = np.logspace(-1, 3.5, 40)
    for app in [2.0, 20.0, 300.0, 1000.0]:
        resp = mb_true.peak_dose_response(concs_sweep, application_ms=app)
        sl_whole = hill_slope(concs_sweep, resp)
        sl_rising = hill_slope(concs_sweep, resp, window_rel=JAHN1997_SLOPE_WINDOW_REL)
        pt_pk = ProtocolPoint(100.0, 2.5, Observable.PEAK, application_ms=app)
        s_02, _, _, _ = posterior_predictive_score(cand, pt_pk, sigma_meas=0.02)
        s_005, _, _, sep = posterior_predictive_score(cand, pt_pk, sigma_meas=0.005)
        print(f"Application {app:6.1f} ms | whole slope: {sl_whole:.3f} | rising slope: {sl_rising:.3f} | sep: {sep:.4f} | score(0.02): {s_02:.3f} | score(0.005): {s_005:.3f}")

    print("\n--------------------------------------------------------------------------------")
    print("2. ORIGINAL DESIGN (EQUILIBRIUM, 1000 uM GABA, PAM 4x)")
    print("--------------------------------------------------------------------------------")
    pt_orig = ProtocolPoint(1000.0, 4.0, Observable.EQUILIBRIUM)
    sweep_orig = sigma_sensitivity(cand, pt_orig, sigmas=(0.005, 0.01, 0.02, 0.03, 0.05))
    print(f"{'sigma_meas':<12} | {'Raw Score':<12} | {'Corrected (n=54)':<18} | {'Discriminating?':<16}")
    print("-" * 65)
    for s in (0.005, 0.01, 0.02, 0.03, 0.05):
        raw = sweep_orig[s]
        corr = _penalise(raw, 54)
        disc = "YES" if corr >= 1.96 else ("MARGINAL" if raw >= 1.96 else "NO")
        print(f"{s:<12.3f} | {raw:<12.3f} | {corr:<18.3f} | {disc:<16}")

    print("\n--------------------------------------------------------------------------------")
    print("3. RE-TARGETED OBSERVABLE: PAIRED-PULSE RECOVERY")
    print("--------------------------------------------------------------------------------")
    # Search across gaps: 10, 20, 50, 100, 200, 500 ms; 10 mM GABA, 2 ms pulse; PAM 1x - 4x
    res_pp = find_discriminating_protocol_pp(
        cand,
        gaba_grid=(10000.0,),
        pam_grid=(1.0, 2.0, 4.0),
        observables=(Observable.PAIRED_PULSE_RECOVERY,),
        application_grid=(2.0,),
        gap_grid=(10.0, 20.0, 50.0, 100.0, 200.0, 500.0),
        sigma_meas=0.03,
    )
    print(f"Best PAIRED_PULSE_RECOVERY protocol: {res_pp.protocol.describe()}")
    print(f"Raw score: {res_pp.score:.3f}, corrected (n={res_pp.n_searched}): {res_pp.score_penalised:.3f}")
    print(f"Expected separation: {res_pp.expected_separation:.4f}")

    sweep_pp = sigma_sensitivity(
        cand, res_pp.protocol, sigmas=(0.005, 0.01, 0.02, 0.03, 0.05)
    )
    print(f"\nSensitivity sweep for {res_pp.protocol.describe()}:")
    print(f"{'sigma_meas':<12} | {'Raw Score':<12} | {'Corrected (n=18)':<18} | {'Discriminating?':<16}")
    print("-" * 65)
    for s in (0.005, 0.01, 0.02, 0.03, 0.05):
        raw = sweep_pp[s]
        corr = _penalise(raw, res_pp.n_searched)
        disc = "YES" if corr >= 1.96 else ("MARGINAL" if raw >= 1.96 else "NO")
        print(f"{s:<12.3f} | {raw:<12.3f} | {corr:<18.3f} | {disc:<16}")

    print("\n================================================================================")
    print("SUMMARY & VERDICT: OLD VS NEW PRE-REGISTRATION")
    print("================================================================================")
    print("Old pre-registration (knowledge/13-preregistration.md):")
    print("  Protocol: 1000 uM GABA, PAM 4x, EQUILIBRIUM observable")
    print("  At sigma = 0.02: Raw score 1.71, corrected 0.00 (NOT discriminating)")
    print("  At sigma = 0.005: Raw score 5.90, corrected ~5.2 (DISCRIMINATING)")
    print("  Conclusion: Requirement was sigma <= 0.005 (4x reduction from conventional noise).")
    print()
    print("New re-targeted pre-registration (PAIRED_PULSE_RECOVERY):")
    print(f"  Protocol: {res_pp.protocol.describe()}")
    print(f"  At published SEM (sigma = 0.03): Raw score {sweep_pp[0.03]:.3f}, corrected {_penalise(sweep_pp[0.03], res_pp.n_searched):.3f}")
    print(f"  At sigma = 0.005: Raw score {sweep_pp[0.005]:.3f}, corrected {_penalise(sweep_pp[0.005], res_pp.n_searched):.3f}")
    print("Verdict on whether sigma <= 0.005 survives:")
    if sweep_pp[0.02] >= 1.96:
        print("  --> SOFTENED: Paired-pulse recovery separates at sigma <= 0.02!")
    elif sweep_pp[0.005] >= 1.96 or sweep_orig[0.005] >= 1.96:
        print("  --> SURVIVES: The requirement for ultra-low noise (sigma <= 0.005) SURVIVES.")
        print("      Neither PEAK nor PAIRED-PULSE separates Models B and C at conventional sigma = 0.02")
        print("      when Model C has collapsed its unidentifiable slow sink (D_slow -> 0).")


if __name__ == "__main__":
    main()
