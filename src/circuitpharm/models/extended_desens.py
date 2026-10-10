"""Model C: Extended Kinetic Scheme with Dual Desensitisation Pathways and State-Dependent Modulation.

Incorporates fast and slow desensitisation branches (D_fast, D_slow) motivated by
structural and functional studies (e.g., Gielen et al., Nature Comms), allowing
modulators to differentially alter desensitisation entry/recovery alongside agonist binding.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp

from .base import (OPEN_STATE_INDEX, PEAK_APPLICATION_MS, DecayFit, Observable,
                   ReceptorModel, WaveformResult, fit_biexponential_decay,
                   peak_open_probability_constant, resolve_initial_state)


STATES_6 = ("R", "AR", "A2R", "A2O", "A2D_fast", "A2D_slow")
OPEN_IDX_6 = 3  # "A2O"
assert STATES_6[OPEN_STATE_INDEX] == "A2O", "open-state index desynchronised from base"


@dataclass(frozen=True)
class ExtendedDesensitizationModel(ReceptorModel):
    """6-state Markov scheme with dual-pathway desensitisation and state-dependent PAM actions.

    EVERY DEFAULT BELOW IS PROVISIONAL AND IS NOT A CALIBRATION. They are round numbers
    with no recorded origin and reproduce neither anchor: measured EC50 5.18 uM
    (EQUILIBRIUM) against the 20 uM anchor, P_o,max 0.8571 against 0.750. Obtain
    parameters from `circuitpharm.parameters.get(...)`; these exist for unit-testing the
    algebra, and the mismatch is pinned by a test so it cannot be mistaken for a fit
    (roadmap P0-13).
    """
    kon: float = 0.012          # uM^-1 ms^-1  PROVISIONAL - see class docstring
    koff: float = 0.35          # ms^-1        PROVISIONAL
    beta: float = 0.60          # ms^-1 (channel opening)  PROVISIONAL
    alpha: float = 0.10         # ms^-1 (channel closing)  PROVISIONAL
    d_fast: float = 0.080       # ms^-1 (rapid desensitisation entry)
    r_fast: float = 0.0050      # ms^-1 (rapid resensitisation)
    d_slow: float = 0.0080      # ms^-1 (slow deep desensitisation entry)
    r_slow: float = 0.0004      # ms^-1 (slow resensitisation recovery)
    # State-dependent modulation parameters
    pam_desens_factor: float = 1.0  # factor by which PAM alters d_fast/r_fast ratio

    @property
    def name(self) -> str:
        return "ExtendedDesensitization_DualD"

    @property
    def native_observable(self) -> Observable:
        """`dose_response` is the stationary distribution, both sinks included."""
        return Observable.EQUILIBRIUM

    @property
    def param_names(self) -> tuple[str, ...]:
        return (
            "kon", "koff", "beta", "alpha",
            "d_fast", "r_fast", "d_slow", "r_slow",
            "pam_desens_factor"
        )

    def get_params(self) -> dict[str, float]:
        return {
            "kon": self.kon,
            "koff": self.koff,
            "beta": self.beta,
            "alpha": self.alpha,
            "d_fast": self.d_fast,
            "r_fast": self.r_fast,
            "d_slow": self.d_slow,
            "r_slow": self.r_slow,
            "pam_desens_factor": self.pam_desens_factor,
        }

    def with_params(self, **kwargs: float) -> "ExtendedDesensitizationModel":
        return replace(self, **kwargs)

    @property
    def kd_um(self) -> float:
        return self.koff / self.kon

    @property
    def gating_efficacy(self) -> float:
        return self.beta / self.alpha

    @property
    def desens_fast_ratio(self) -> float:
        return self.d_fast / self.r_fast

    @property
    def desens_slow_ratio(self) -> float:
        return self.d_slow / self.r_slow

    @property
    def total_desens_ratio(self) -> float:
        return self.desens_fast_ratio + self.desens_slow_ratio

    @property
    def po_max(self) -> float:
        return self.beta / (self.beta + self.alpha)

    @property
    def po_inf(self) -> float:
        e = self.gating_efficacy
        d_tot = self.total_desens_ratio
        return e / (1.0 + e + d_tot)

    def apply_pam(self, pam_factor: float = 1.0) -> "ExtendedDesensitizationModel":
        """Apply state-dependent modulation.

        `pam_factor` acts as an affinity enhancer on koff (an EQUILIBRIUM EC50 fold-shift,
        as in Model B), while `pam_desens_factor` scales how much that shift also slows
        entry into the fast desensitised state.

        POSITIVE MODULATION ONLY, and the desensitisation denominator is GUARDED. Both
        matter. `pam_desens_factor` is in `param_names`, so a fitter or sampler can move
        it: with the old expression `d_fast / (1 + 0.1*(pf-1)*f)` a value of f = 20 at
        pf = 0.5 drives the denominator to exactly zero (infinite rate), and larger f with
        pf < 1 drives it negative -- an unphysical negative rate constant that propagates
        silently. The same class of defect as the clamps already in
        `cpg.Drug.nmda_scale` and `substrate.tonic_gaba`, which are the precedent for
        guarding it here rather than documenting it (roadmap P0-12).
        """
        pf = float(pam_factor)
        if not np.isfinite(pf):
            raise ValueError(f"pam_factor must be finite, got {pf!r}")
        if pf < 1.0:
            raise ValueError(
                f"pam_factor={pf!r} is below 1.0, i.e. negative allosteric modulation, "
                f"which this scheme's callers are not written for. Pass >= 1.0.")
        denom = 1.0 + 0.1 * (pf - 1.0) * self.pam_desens_factor
        if denom <= 1e-9:
            raise ValueError(
                f"pam_desens_factor={self.pam_desens_factor!r} with pam_factor={pf!r} "
                f"gives a desensitisation scaling denominator of {denom!r}, which would "
                f"make d_fast infinite or negative. A negative rate constant is not a "
                f"parameter choice; bound pam_desens_factor so "
                f"1 + 0.1*(pam_factor-1)*pam_desens_factor stays positive.")
        return replace(
            self,
            koff=self.koff / pf,
            # Dual-pathway state dependence: the PAM stabilises the open state relative to
            # fast desensitisation.
            d_fast=self.d_fast / denom,
        )

    def state_distribution(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        """Equilibrium probability vector P = [P_R, P_AR, P_A2R, P_A2O, P_A2D_fast, P_A2D_slow]."""
        mod = self.apply_pam(pam_factor=pam_factor)
        g = max(float(gaba_um), 0.0)
        x = g / mod.kd_um
        e = mod.gating_efficacy
        df = mod.desens_fast_ratio
        ds = mod.desens_slow_ratio

        denom = 1.0 + 2.0 * x + (x ** 2) * (1.0 + e + df + ds)
        if denom <= 0.0:
            return np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        p_r = 1.0 / denom
        p_ar = (2.0 * x) / denom
        p_a2r = (x ** 2) / denom
        p_a2o = (e * (x ** 2)) / denom
        p_df = (df * (x ** 2)) / denom
        p_ds = (ds * (x ** 2)) / denom

        return np.array([p_r, p_ar, p_a2r, p_a2o, p_df, p_ds], dtype=float)

    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float:
        dist = self.state_distribution(gaba_um, pam_factor)
        return float(dist[OPEN_IDX_6])

    def q_matrix(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        """Row-generator transition matrix Q: dP/dt = P Q.

        The same generator `simulate_waveform`'s `rhs` writes out term by term. It is
        written twice because the ODE path evaluates a time-varying concentration and
        building a 6x6 per step would cost more than the eight multiplications it saves,
        and `tests/test_p4_likelihood.py` asserts the two agree at several concentrations
        so the duplication cannot drift -- which is the only acceptable form of it.
        """
        mod = self.apply_pam(pam_factor=pam_factor)
        g = max(float(gaba_um), 0.0)
        kon, koff = mod.kon, mod.koff
        beta, alpha = mod.beta, mod.alpha
        df, rf = mod.d_fast, mod.r_fast
        ds, rs = mod.d_slow, mod.r_slow

        return np.array([
            [-2.0 * kon * g, 2.0 * kon * g, 0.0, 0.0, 0.0, 0.0],
            [koff, -(koff + kon * g), kon * g, 0.0, 0.0, 0.0],
            [0.0, 2.0 * koff, -(2.0 * koff + beta + df + ds), beta, df, ds],
            [0.0, 0.0, alpha, -alpha, 0.0, 0.0],
            [0.0, 0.0, rf, 0.0, -rf, 0.0],
            [0.0, 0.0, rs, 0.0, 0.0, -rs],
        ], dtype=float)

    def dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0) -> np.ndarray:
        """EQUILIBRIUM open probability, both desensitisation sinks included.

        Not the peak-current CRC that published data reports -- see Model B's note and
        roadmap §2.3. Use `peak_dose_response` for PEAK data.
        """
        concs = np.asarray(concs_um, dtype=float)
        mod = self.apply_pam(pam_factor=pam_factor)
        x = np.maximum(concs, 0.0) / mod.kd_um
        e = mod.gating_efficacy
        df = mod.desens_fast_ratio
        ds = mod.desens_slow_ratio
        denom = 1.0 + 2.0 * x + (x ** 2) * (1.0 + e + df + ds)
        return (e * (x ** 2)) / denom

    def peak_dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0,
                           application_ms: float = PEAK_APPLICATION_MS) -> np.ndarray:
        """PEAK open probability during a square application. See the protocol docstring."""
        # Solved exactly with one matrix exponential per concentration rather than an ODE
        # solve: at constant agonist the generator is constant. Agrees with the previous
        # solve_ivp implementation to 1.8e-7 and is ~300x faster -- see
        # `base.peak_open_probability_constant` for the measurement and for why `expm`
        # rather than an eigendecomposition.
        return peak_open_probability_constant(
            self.q_matrix, self.state_distribution(0.0), concs_um,
            pam_factor=pam_factor, application_ms=application_ms,
            open_index=OPEN_IDX_6)

    def simulate_waveform(
        self,
        t: np.ndarray,
        gaba_t: np.ndarray,
        pam_factor: float = 1.0,
        initial_state: np.ndarray | None = None,
    ) -> WaveformResult:
        """Integrate 6-state dual-desensitisation scheme across dynamic GABA waveform."""
        t_arr = np.asarray(t, dtype=float)
        gaba_arr = np.asarray(gaba_t, dtype=float)

        if len(t_arr) < 2:
            dist = (self.state_distribution(float(gaba_arr[0]), pam_factor)
                    if len(gaba_arr) > 0 else np.array([1.0, 0, 0, 0, 0, 0]))
            p0 = float(dist[OPEN_IDX_6])
            # decay is UNDEFINED on a single sample -> NaN, never a nominal tau (P0-4)
            return WaveformResult(
                t=t_arr, gaba=gaba_arr, p_open=np.full_like(t_arr, p0), peak_p_open=p0,
                charge_integral=0.0, decay_tau_ms=float("nan"), steady_state_tail=p0,
                states=dist[None, :], decay=DecayFit.failed())

        # P0-1: scalar honoured, vector validated, anything else raises. `len()` on a
        # float used to crash here before the comparison was even reached.
        p_init = resolve_initial_state(
            initial_state, 6, self.state_distribution(float(gaba_arr[0]), pam_factor))

        mod = self.apply_pam(pam_factor=pam_factor)
        kon = mod.kon
        koff = mod.koff
        beta = mod.beta
        alpha = mod.alpha
        df = mod.d_fast
        rf = mod.r_fast
        ds = mod.d_slow
        rs = mod.r_slow

        def rhs(time_val: float, y: np.ndarray) -> np.ndarray:
            g = float(np.interp(time_val, t_arr, gaba_arr))
            dp0 = -2.0 * kon * g * y[0] + koff * y[1]
            dp1 = 2.0 * kon * g * y[0] - (koff + kon * g) * y[1] + 2.0 * koff * y[2]
            dp2 = (
                kon * g * y[1]
                - (2.0 * koff + beta + df + ds) * y[2]
                + alpha * y[3]
                + rf * y[4]
                + rs * y[5]
            )
            dp3 = beta * y[2] - alpha * y[3]
            dp4 = df * y[2] - rf * y[4]
            dp5 = ds * y[2] - rs * y[5]
            return np.array([dp0, dp1, dp2, dp3, dp4, dp5])

        sol = solve_ivp(
            rhs,
            (float(t_arr[0]), float(t_arr[-1])),
            p_init,
            t_eval=t_arr,
            method="Radau",
            rtol=1e-5,
            atol=1e-7,
        )

        # P0-3: a failed integration RAISES rather than being replaced by a tiled
        # constant trace, which looked like a settled experiment that never ran.
        if not sol.success:
            raise RuntimeError(
                f"Radau failed to integrate the 6-state scheme over "
                f"[{t_arr[0]:.4g}, {t_arr[-1]:.4g}] ms: {sol.message}. The slowest "
                f"timescale here is 1/r_slow = {1.0 / self.r_slow:.0f} ms against a cleft "
                f"transient of ~1 ms, so this is a stiffer problem than the 5-state "
                f"scheme; densify `t` rather than loosening rtol.")

        states = sol.y.T
        row_sums = states.sum(axis=1)
        worst = float(np.max(np.abs(row_sums - 1.0))) if row_sums.size else 0.0
        if worst > 1e-3:
            raise RuntimeError(
                f"state probabilities drifted off the simplex by {worst:.3e} "
                f"(max |sum - 1|) while integrating the 6-state scheme. Renormalising "
                f"would hide an integration failure behind a plausible trace.")
        states = np.clip(states, 0.0, 1.0)
        states = states / np.maximum(states.sum(axis=1, keepdims=True), 1e-12)

        p_open = states[:, OPEN_IDX_6]
        peak_p = float(np.max(p_open))

        _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
        charge = float(_trapz(p_open, t_arr))

        # P0-4: real biexponential fit, NaN on failure -- never 15.0, the tau fit target.
        decay = fit_biexponential_decay(t_arr, p_open)

        return WaveformResult(
            t=t_arr,
            gaba=gaba_arr,
            p_open=p_open,
            peak_p_open=peak_p,
            charge_integral=charge,
            decay_tau_ms=decay.tau_weighted_ms,
            steady_state_tail=float(p_open[-1]),
            states=states,
            decay=decay,
        )
