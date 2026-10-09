"""Model C: Extended Kinetic Scheme with Dual Desensitisation Pathways and State-Dependent Modulation.

Incorporates fast and slow desensitisation branches (D_fast, D_slow) motivated by
structural and functional studies (e.g., Gielen et al., Nature Comms), allowing
modulators to differentially alter desensitisation entry/recovery alongside agonist binding.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp

from .base import ReceptorModel, WaveformResult


STATES_6 = ("R", "AR", "A2R", "A2O", "A2D_fast", "A2D_slow")
OPEN_IDX_6 = 3  # "A2O"


@dataclass(frozen=True)
class ExtendedDesensitizationModel(ReceptorModel):
    """6-state Markov scheme with dual-pathway desensitisation and state-dependent PAM actions."""
    kon: float = 0.012          # uM^-1 ms^-1
    koff: float = 0.35          # ms^-1
    beta: float = 0.60          # ms^-1 (channel opening)
    alpha: float = 0.10         # ms^-1 (channel closing)
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
        
        pam_factor acts as affinity enhancer on koff, while pam_desens_factor
        can modulate desensitisation trapping.
        """
        pf = max(float(pam_factor), 1e-6)
        return replace(
            self,
            koff=self.koff / pf,
            # Dual-pathway state dependence: PAM slightly stabilizes open state over fast desensitisation
            d_fast=self.d_fast / (1.0 + 0.1 * (pf - 1.0) * self.pam_desens_factor),
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

    def dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0) -> np.ndarray:
        concs = np.asarray(concs_um, dtype=float)
        mod = self.apply_pam(pam_factor=pam_factor)
        x = np.maximum(concs, 0.0) / mod.kd_um
        e = mod.gating_efficacy
        df = mod.desens_fast_ratio
        ds = mod.desens_slow_ratio
        denom = 1.0 + 2.0 * x + (x ** 2) * (1.0 + e + df + ds)
        return (e * (x ** 2)) / denom

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
            dist = self.state_distribution(float(gaba_arr[0]), pam_factor) if len(gaba_arr) > 0 else np.array([1.0, 0, 0, 0, 0, 0])
            p0 = float(dist[OPEN_IDX_6])
            return WaveformResult(t_arr, gaba_arr, np.full_like(t_arr, p0), p0, 0.0, 0.0, p0, dist[None, :])

        if initial_state is not None and len(initial_state) == 6:
            p_init = np.asarray(initial_state, dtype=float)
        else:
            p_init = self.state_distribution(float(gaba_arr[0]), pam_factor)

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

        states = sol.y.T if sol.success else np.tile(p_init, (len(t_arr), 1))
        states = np.clip(states, 0.0, 1.0)
        row_sums = states.sum(axis=1, keepdims=True)
        states = np.divide(states, np.maximum(row_sums, 1e-12))

        p_open = states[:, OPEN_IDX_6]
        peak_p = float(np.max(p_open))

        _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
        charge = float(_trapz(p_open, t_arr))

        peak_idx = int(np.argmax(p_open))
        if peak_idx < len(t_arr) - 2:
            t_tail = t_arr[peak_idx:] - t_arr[peak_idx]
            p_tail = p_open[peak_idx:] - p_open[-1]
            if np.max(p_tail) > 1e-4:
                half_mask = p_tail <= 0.5 * p_tail[0]
                decay_tau = float(t_tail[half_mask][0] / np.log(2)) if np.any(half_mask) else 15.0
            else:
                decay_tau = 15.0
        else:
            decay_tau = 15.0

        return WaveformResult(
            t=t_arr,
            gaba=gaba_arr,
            p_open=p_open,
            peak_p_open=peak_p,
            charge_integral=charge,
            decay_tau_ms=decay_tau,
            steady_state_tail=float(p_open[-1]),
            states=states,
        )
