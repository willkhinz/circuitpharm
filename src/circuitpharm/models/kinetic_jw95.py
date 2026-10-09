"""Model B: 5-State Del Castillo-Katz / Jones-Westbrook (1995) Kinetic Allostery Model.

Positive allosteric modulation acts by shifting microscopic rate constants
(e.g., slowing agonist unbinding koff, accelerating channel opening beta),
allowing open-probability ceilings and compartment differences to emerge naturally
from state redistribution.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp

from .base import ReceptorModel, WaveformResult


STATES_5 = ("R", "AR", "A2R", "A2O", "A2D")
OPEN_IDX_5 = 3  # "A2O"


@dataclass(frozen=True)
class KineticAllosteryModel(ReceptorModel):
    """5-state Markov gating scheme with microscopic allosteric modulation."""
    kon: float = 0.011244       # uM^-1 ms^-1 (association rate)
    koff: float = 0.333069      # ms^-1 (dissociation rate)
    beta: float = 0.559023      # ms^-1 (channel opening rate)
    alpha: float = 0.107891     # ms^-1 (channel closing rate)
    d: float = 0.050            # ms^-1 (desensitisation rate)
    r: float = 0.0020           # ms^-1 (resensitisation rate)

    @property
    def name(self) -> str:
        return "KineticAllostery_JW95"

    @property
    def param_names(self) -> tuple[str, ...]:
        return ("kon", "koff", "beta", "alpha", "d", "r")

    def get_params(self) -> dict[str, float]:
        return {
            "kon": self.kon,
            "koff": self.koff,
            "beta": self.beta,
            "alpha": self.alpha,
            "d": self.d,
            "r": self.r,
        }

    def with_params(self, **kwargs: float) -> "KineticAllosteryModel":
        return replace(self, **kwargs)

    @property
    def kd_um(self) -> float:
        """Microscopic dissociation constant Kd = koff / kon in uM."""
        return self.koff / self.kon

    @property
    def gating_efficacy(self) -> float:
        """E = beta / alpha."""
        return self.beta / self.alpha

    @property
    def desens_ratio(self) -> float:
        """D = d / r."""
        return self.d / self.r

    @property
    def po_max(self) -> float:
        """Theoretical maximum instantaneous open probability at saturating GABA before desensitisation."""
        return self.beta / (self.beta + self.alpha)

    @property
    def po_inf(self) -> float:
        """Steady-state open probability at infinite GABA after desensitisation equilibrium."""
        e = self.gating_efficacy
        d = self.desens_ratio
        return e / (1.0 + e + d)

    def apply_pam(self, affinity_factor: float = 1.0, gating_factor: float = 1.0) -> "KineticAllosteryModel":
        """Apply allosteric modulation.

        affinity_factor: fold-decrease in koff (slowing unbinding).
        gating_factor: fold-increase in beta (accelerating opening).
        """
        return replace(
            self,
            koff=self.koff / max(float(affinity_factor), 1e-6),
            beta=self.beta * float(gating_factor),
        )

    def q_matrix(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        """Row-generator transition matrix Q: dP/dt = P Q."""
        mod = self.apply_pam(affinity_factor=pam_factor)
        g = max(float(gaba_um), 0.0)
        kon = mod.kon
        koff = mod.koff
        beta = mod.beta
        alpha = mod.alpha
        d = mod.d
        r = mod.r

        q = np.array([
            [-2.0 * kon * g,  2.0 * kon * g,        0.0,              0.0,    0.0],
            [koff,            -(koff + kon * g),     kon * g,          0.0,    0.0],
            [0.0,             2.0 * koff,           -(2.0*koff+beta+d), beta,   d],
            [0.0,             0.0,                  alpha,            -alpha, 0.0],
            [0.0,             0.0,                  r,                0.0,    -r],
        ], dtype=float)
        return q

    def state_distribution(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        """Equilibrium probability vector P = [P_R, P_AR, P_A2R, P_A2O, P_A2D]."""
        mod = self.apply_pam(affinity_factor=pam_factor)
        g = max(float(gaba_um), 0.0)
        x = g / mod.kd_um
        e = mod.gating_efficacy
        d = mod.desens_ratio

        denom = 1.0 + 2.0 * x + (x ** 2) * (1.0 + e + d)
        if denom <= 0.0:
            return np.array([1.0, 0.0, 0.0, 0.0, 0.0])

        p_r = 1.0 / denom
        p_ar = (2.0 * x) / denom
        p_a2r = (x ** 2) / denom
        p_a2o = (e * (x ** 2)) / denom
        p_a2d = (d * (x ** 2)) / denom

        return np.array([p_r, p_ar, p_a2r, p_a2o, p_a2d], dtype=float)

    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float:
        """Stationary open probability P(A2O)."""
        dist = self.state_distribution(gaba_um, pam_factor)
        return float(dist[OPEN_IDX_5])

    def dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0) -> np.ndarray:
        concs = np.asarray(concs_um, dtype=float)
        mod = self.apply_pam(affinity_factor=pam_factor)
        x = np.maximum(concs, 0.0) / mod.kd_um
        e = mod.gating_efficacy
        d = mod.desens_ratio
        denom = 1.0 + 2.0 * x + (x ** 2) * (1.0 + e + d)
        return (e * (x ** 2)) / denom

    def simulate_waveform(
        self,
        t: np.ndarray,
        gaba_t: np.ndarray,
        pam_factor: float = 1.0,
        initial_state: np.ndarray | None = None,
    ) -> WaveformResult:
        """Integrate 5-state Markov scheme across time-varying GABA transient."""
        t_arr = np.asarray(t, dtype=float)
        gaba_arr = np.asarray(gaba_t, dtype=float)

        if len(t_arr) < 2:
            dist = self.state_distribution(float(gaba_arr[0]), pam_factor) if len(gaba_arr) > 0 else np.array([1.0, 0, 0, 0, 0])
            p0 = float(dist[OPEN_IDX_5])
            return WaveformResult(t_arr, gaba_arr, np.full_like(t_arr, p0), p0, 0.0, 0.0, p0, dist[None, :])

        if initial_state is not None:
            p_init = np.asarray(initial_state, dtype=float)
            if len(p_init) != 5:
                # If scalar open prob provided, distribute according to resting steady state
                p_init = self.state_distribution(float(gaba_arr[0]), pam_factor)
        else:
            p_init = self.state_distribution(float(gaba_arr[0]), pam_factor)

        # ODE: dP/dt = P * Q(gaba(t))
        # Note: solve_ivp works with column state vectors y, so dy/dt = Q^T y
        mod = self.apply_pam(affinity_factor=pam_factor)
        kon = mod.kon
        koff = mod.koff
        beta = mod.beta
        alpha = mod.alpha
        d_rate = mod.d
        r_rate = mod.r

        def rhs(time_val: float, y: np.ndarray) -> np.ndarray:
            g = float(np.interp(time_val, t_arr, gaba_arr))
            # Row vector P: dP_R = -2 kon g P_R + koff P_AR
            dp0 = -2.0 * kon * g * y[0] + koff * y[1]
            dp1 = 2.0 * kon * g * y[0] - (koff + kon * g) * y[1] + 2.0 * koff * y[2]
            dp2 = kon * g * y[1] - (2.0 * koff + beta + d_rate) * y[2] + alpha * y[3] + r_rate * y[4]
            dp3 = beta * y[2] - alpha * y[3]
            dp4 = d_rate * y[2] - r_rate * y[4]
            return np.array([dp0, dp1, dp2, dp3, dp4])

        sol = solve_ivp(
            rhs,
            (float(t_arr[0]), float(t_arr[-1])),
            p_init,
            t_eval=t_arr,
            method="Radau",  # Radau handles stiff multi-timescale transitions accurately
            rtol=1e-5,
            atol=1e-7,
        )

        states = sol.y.T if sol.success else np.tile(p_init, (len(t_arr), 1))
        # Enforce non-negativity and simplex constraint
        states = np.clip(states, 0.0, 1.0)
        row_sums = states.sum(axis=1, keepdims=True)
        states = np.divide(states, np.maximum(row_sums, 1e-12))

        p_open = states[:, OPEN_IDX_5]
        peak_p = float(np.max(p_open))

        _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
        charge = float(_trapz(p_open, t_arr))

        # Decay tau post-peak
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
