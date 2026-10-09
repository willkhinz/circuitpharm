"""Model A: Phenomenological Operational Scalar Ceiling Model.

Represents the hypothesis that positive allosteric modulation acts via a simple
operational potency/efficacy shift subject to an extrinsic scalar ceiling, without
explicit microscopic Markov state redistribution.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp

from .base import ReceptorModel, WaveformResult


@dataclass(frozen=True)
class OperationalScalarModel(ReceptorModel):
    """Operational / Hill model with an imposed scalar ceiling on allosteric gain."""
    ec50_um: float = 25.0       # basal GABA EC50 (uM)
    hill_n: float = 1.4         # Hill coefficient
    po_max: float = 0.75        # maximum open probability ceiling
    s_max: float = 2.50         # hard ceiling on operational potency shift
    tau_act_ms: float = 0.5     # activation time constant (ms)
    tau_deact_ms: float = 15.0  # deactivation time constant (ms)

    @property
    def name(self) -> str:
        return "OperationalScalar"

    @property
    def param_names(self) -> tuple[str, ...]:
        return ("ec50_um", "hill_n", "po_max", "s_max", "tau_act_ms", "tau_deact_ms")

    def get_params(self) -> dict[str, float]:
        return {
            "ec50_um": self.ec50_um,
            "hill_n": self.hill_n,
            "po_max": self.po_max,
            "s_max": self.s_max,
            "tau_act_ms": self.tau_act_ms,
            "tau_deact_ms": self.tau_deact_ms,
        }

    def with_params(self, **kwargs: float) -> "OperationalScalarModel":
        return replace(self, **kwargs)

    def effective_ec50(self, pam_factor: float = 1.0) -> float:
        """Apply operational potency shift with ceiling s_max."""
        shift = min(float(pam_factor), self.s_max)
        return self.ec50_um / max(shift, 1e-6)

    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float:
        """Stationary Hill open probability."""
        if gaba_um <= 0.0:
            return 0.0
        eff_ec50 = self.effective_ec50(pam_factor)
        conc_ratio = (gaba_um / eff_ec50) ** self.hill_n
        return float(self.po_max * (conc_ratio / (1.0 + conc_ratio)))

    def dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0) -> np.ndarray:
        concs = np.asarray(concs_um, dtype=float)
        eff_ec50 = self.effective_ec50(pam_factor)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = np.where(concs > 0, (concs / eff_ec50) ** self.hill_n, 0.0)
            res = self.po_max * (ratios / (1.0 + ratios))
        return np.nan_to_num(res, nan=0.0)

    def simulate_waveform(
        self,
        t: np.ndarray,
        gaba_t: np.ndarray,
        pam_factor: float = 1.0,
        initial_state: np.ndarray | None = None,
    ) -> WaveformResult:
        """Simulate dynamic open probability relaxation towards instantaneous Hill target."""
        t_arr = np.asarray(t, dtype=float)
        gaba_arr = np.asarray(gaba_t, dtype=float)

        if len(t_arr) < 2:
            p0 = self.steady_state(float(gaba_arr[0]), pam_factor) if len(gaba_arr) > 0 else 0.0
            return WaveformResult(t_arr, gaba_arr, np.full_like(t_arr, p0), p0, 0.0, 0.0, p0)

        p_init = float(initial_state[0]) if initial_state is not None else self.steady_state(float(gaba_arr[0]), pam_factor)
        
        # Precompute target at all timepoints
        target_p = self.dose_response(gaba_arr, pam_factor)

        # Asymmetric relaxation ODE: fast activation, slower deactivation
        def ode(time_val: float, y: np.ndarray) -> np.ndarray:
            # Interpolate target
            p_targ = np.interp(time_val, t_arr, target_p)
            tau = self.tau_act_ms if p_targ >= y[0] else self.tau_deact_ms
            return np.array([(p_targ - y[0]) / max(tau, 1e-4)])

        dt_min = float(np.min(np.diff(t_arr))) if len(t_arr) > 1 else 1.0
        max_step = min(self.tau_act_ms, dt_min, 0.5)

        sol = solve_ivp(
            ode,
            (float(t_arr[0]), float(t_arr[-1])),
            [p_init],
            t_eval=t_arr,
            method="RK23",
            max_step=max_step,
            rtol=1e-5,
            atol=1e-7,
        )

        p_open = np.clip(sol.y[0] if sol.success else target_p, 0.0, 1.0)
        peak_p = float(np.max(p_open))
        
        # Integral via trapezoid
        _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
        charge = float(_trapz(p_open, t_arr))

        # Estimate decay tau post-peak
        peak_idx = int(np.argmax(p_open))
        if peak_idx < len(t_arr) - 2:
            t_tail = t_arr[peak_idx:] - t_arr[peak_idx]
            p_tail = p_open[peak_idx:] - p_open[-1]
            if np.max(p_tail) > 1e-4:
                half_mask = p_tail <= 0.5 * p_tail[0]
                decay_tau = float(t_tail[half_mask][0] / np.log(2)) if np.any(half_mask) else self.tau_deact_ms
            else:
                decay_tau = self.tau_deact_ms
        else:
            decay_tau = self.tau_deact_ms

        return WaveformResult(
            t=t_arr,
            gaba=gaba_arr,
            p_open=p_open,
            peak_p_open=peak_p,
            charge_integral=charge,
            decay_tau_ms=decay_tau,
            steady_state_tail=float(p_open[-1]),
        )
