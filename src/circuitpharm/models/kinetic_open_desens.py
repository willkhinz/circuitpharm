"""Alternative Model: 5-State Kinetic Scheme with Open-State Desensitisation.

Unlike Model B (KineticAllosteryModel), where desensitisation connects to the doubly-bound
closed state A2R (A2R <-> A2D), here desensitisation connects strictly to the open channel
state A2O:

    R <--> AR <--> A2R <--> A2O (open) <--> A2D (desensitised)

This provides an explicit alternative topology for Brief 5 Task C to test whether the
conflict between fast deactivation and slow paired-pulse recovery is resolved by
open-channel desensitisation.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp

from .base import (
    OPEN_STATE_INDEX,
    PEAK_APPLICATION_MS,
    DecayFit,
    Observable,
    ObservableMismatch,
    ReceptorModel,
    WaveformResult,
    fit_biexponential_decay,
    peak_open_probability_constant,
    resolve_initial_state,
)


STATES_5 = ("R", "AR", "A2R", "A2O", "A2D")
OPEN_IDX_5 = 3
assert STATES_5[OPEN_STATE_INDEX] == "A2O", "open-state index desynchronised from base"


@dataclass(frozen=True)
class KineticOpenDesensModel(ReceptorModel):
    """5-state Markov gating scheme with desensitisation accessed from the OPEN state.

    Topology:
        R <--> AR <--> A2R <--> A2O <--> A2D
    """
    kon: float = 0.011244       # uM^-1 ms^-1
    koff: float = 0.333069      # ms^-1
    beta: float = 0.559023      # ms^-1
    alpha: float = 0.107891     # ms^-1
    d: float = 0.050            # ms^-1 (A2O -> A2D)
    r: float = 0.0020           # ms^-1 (A2D -> A2O)

    @property
    def name(self) -> str:
        return "KineticOpenDesens"

    @property
    def native_observable(self) -> Observable:
        return Observable.EQUILIBRIUM

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

    def with_params(self, **kwargs: float) -> "KineticOpenDesensModel":
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
        """Theoretical maximum instantaneous open probability before desensitisation."""
        return self.beta / (self.beta + self.alpha)

    @property
    def po_inf(self) -> float:
        """Steady-state open probability at infinite GABA after desensitisation equilibrium."""
        e = self.gating_efficacy
        d = self.desens_ratio
        return e / (1.0 + e + e * d)

    def apply_pam(
        self,
        affinity_factor: float = 1.0,
        gating_factor: float = 1.0,
        *,
        pam_factor: float | None = None,
    ) -> "KineticOpenDesensModel":
        """Apply allosteric modulation."""
        if pam_factor is not None:
            affinity_factor = pam_factor
        for nm, v in (("affinity_factor", affinity_factor), ("gating_factor", gating_factor)):
            v = float(v)
            if not np.isfinite(v):
                raise ValueError(f"{nm} must be finite, got {v!r}")
            if v < 1.0:
                raise ValueError(f"{nm}={v!r} is below 1.0. Pass >= 1.0.")
        return replace(
            self,
            koff=self.koff / float(affinity_factor),
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
            [-2.0 * kon * g,  2.0 * kon * g,        0.0,              0.0,            0.0],
            [koff,            -(koff + kon * g),     kon * g,          0.0,            0.0],
            [0.0,             2.0 * koff,           -(2.0*koff+beta), beta,           0.0],
            [0.0,             0.0,                  alpha,            -(alpha + d),   d],
            [0.0,             0.0,                  0.0,              r,              -r],
        ], dtype=float)
        return q

    def state_distribution(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        """Equilibrium probability vector P = [P_R, P_AR, P_A2R, P_A2O, P_A2D]."""
        mod = self.apply_pam(affinity_factor=pam_factor)
        g = max(float(gaba_um), 0.0)
        x = g / mod.kd_um
        e = mod.gating_efficacy
        d = mod.desens_ratio

        denom = 1.0 + 2.0 * x + (x ** 2) * (1.0 + e + e * d)
        if denom <= 0.0:
            return np.array([1.0, 0.0, 0.0, 0.0, 0.0])

        p_r = 1.0 / denom
        p_ar = (2.0 * x) / denom
        p_a2r = (x ** 2) / denom
        p_a2o = (e * (x ** 2)) / denom
        p_a2d = (e * d * (x ** 2)) / denom

        return np.array([p_r, p_ar, p_a2r, p_a2o, p_a2d], dtype=float)

    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float:
        """Stationary open probability P(A2O)."""
        dist = self.state_distribution(gaba_um, pam_factor)
        return float(dist[OPEN_IDX_5])

    def dose_response(
        self,
        concs_um: np.ndarray,
        pam_factor: float = 1.0,
        *,
        observable: Observable = Observable.EQUILIBRIUM,
    ) -> np.ndarray:
        """EQUILIBRIUM open probability."""
        if observable is not Observable.EQUILIBRIUM:
            raise ObservableMismatch(
                f"{self.name} is EQUILIBRIUM-native and has no {observable.value} curve in "
                f"dose_response. Use peak_dose_response for PEAK data."
            )
        concs = np.asarray(concs_um, dtype=float)
        mod = self.apply_pam(affinity_factor=pam_factor)
        x = np.maximum(concs, 0.0) / mod.kd_um
        e = mod.gating_efficacy
        d = mod.desens_ratio
        denom = 1.0 + 2.0 * x + (x ** 2) * (1.0 + e + e * d)
        return (e * (x ** 2)) / denom

    def peak_dose_response(
        self,
        concs_um: np.ndarray,
        pam_factor: float = 1.0,
        application_ms: float = PEAK_APPLICATION_MS,
    ) -> np.ndarray:
        """PEAK open probability during a square application."""
        return peak_open_probability_constant(
            self.q_matrix,
            self.state_distribution(0.0),
            concs_um,
            pam_factor=pam_factor,
            application_ms=application_ms,
            open_index=OPEN_IDX_5,
        )

    def simulate_waveform(
        self,
        t: np.ndarray,
        gaba_t: np.ndarray,
        pam_factor: float = 1.0,
        initial_state: np.ndarray | None = None,
    ) -> WaveformResult:
        """Integrate 5-state Markov scheme with open-state desensitisation."""
        t_arr = np.asarray(t, dtype=float)
        gaba_arr = np.asarray(gaba_t, dtype=float)

        if len(t_arr) < 2:
            dist = (
                self.state_distribution(float(gaba_arr[0]), pam_factor)
                if len(gaba_arr) > 0
                else np.array([1.0, 0, 0, 0, 0])
            )
            p0 = float(dist[OPEN_IDX_5])
            return WaveformResult(
                t=t_arr,
                gaba=gaba_arr,
                p_open=np.full_like(t_arr, p0),
                peak_p_open=p0,
                charge_integral=0.0,
                decay_tau_ms=float("nan"),
                steady_state_tail=p0,
                states=dist[None, :],
                decay=DecayFit.failed(),
            )

        p_init = resolve_initial_state(
            initial_state, 5, self.state_distribution(float(gaba_arr[0]), pam_factor)
        )

        mod = self.apply_pam(affinity_factor=pam_factor)
        kon = mod.kon
        koff = mod.koff
        beta = mod.beta
        alpha = mod.alpha
        d_rate = mod.d
        r_rate = mod.r

        def rhs(time_val: float, y: np.ndarray) -> np.ndarray:
            g = float(np.interp(time_val, t_arr, gaba_arr))
            dp0 = -2.0 * kon * g * y[0] + koff * y[1]
            dp1 = 2.0 * kon * g * y[0] - (koff + kon * g) * y[1] + 2.0 * koff * y[2]
            dp2 = kon * g * y[1] - (2.0 * koff + beta) * y[2] + alpha * y[3]
            dp3 = beta * y[2] - (alpha + d_rate) * y[3] + r_rate * y[4]
            dp4 = d_rate * y[3] - r_rate * y[4]
            return np.array([dp0, dp1, dp2, dp3, dp4])

        sol = solve_ivp(
            rhs,
            (float(t_arr[0]), float(t_arr[-1])),
            p_init,
            t_eval=t_arr,
            method="Radau",
            rtol=1e-5,
            atol=1e-7,
        )

        if not sol.success:
            raise RuntimeError(
                f"Radau failed to integrate KineticOpenDesensModel over "
                f"[{t_arr[0]:.4g}, {t_arr[-1]:.4g}] ms: {sol.message}."
            )

        states = sol.y.T
        row_sums = states.sum(axis=1)
        worst = float(np.max(np.abs(row_sums - 1.0))) if row_sums.size else 0.0
        if worst > 1e-3:
            raise RuntimeError(
                f"state probabilities drifted off the simplex by {worst:.3e} "
                f"while integrating KineticOpenDesensModel."
            )

        states = np.clip(states, 0.0, 1.0)
        states = states / np.maximum(states.sum(axis=1, keepdims=True), 1e-12)

        p_open = states[:, OPEN_IDX_5]
        peak_p = float(np.max(p_open))

        _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
        charge = float(_trapz(p_open, t_arr))

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
