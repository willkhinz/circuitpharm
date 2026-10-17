"""6-State Kinetic Scheme with Pre-Open ("Flip") Conformational State.

Agent Brief 8B:
Investigates whether adding a pre-open, shut, fully-liganded intermediate between
doubly-bound closed state A2R and open state A2O produces a fast, large-amplitude
relaxation component and so reproduces the published triexponential deactivation of
alpha1beta2gamma2 (Barberis et al. 2007).

Topology:
    R  <-->  AR  <-->  A2R  <-->  A2F  <-->  A2O
                         |
                       A2D

States and Order:
    0: "R"    -- unliganded closed receptor
    1: "AR"   -- singly-liganded closed receptor
    2: "A2R"  -- doubly-liganded resting shut receptor
    3: "A2F"  -- doubly-liganded shut pre-open intermediate ("flipped")
    4: "A2O"  -- doubly-liganded open conducting channel (OPEN STATE, open_index = 4)
    5: "A2D"  -- doubly-liganded desensitised shut receptor (connects strictly to A2R)

Microscopic Reversibility:
    The scheme's transition graph is an acyclic tree (loop-free).
    Because there are no closed cycles in the network, Kolmogorov's cycle condition
    is vacuously satisfied. Consequently, detailed balance (microscopic reversibility)
    holds unconditionally for any set of strictly positive rate constants; there are
    no algebraic loop constraints beyond positivity.
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


STATES_6 = ("R", "AR", "A2R", "A2F", "A2O", "A2D")
OPEN_IDX_6 = 4  # "A2O"


@dataclass(frozen=True)
class KineticFlipModel(ReceptorModel):
    """6-state Markov gating scheme with pre-open 'flip' state A2F.

    EVERY DEFAULT BELOW IS PROVISIONAL AND IS NOT A CALIBRATION.
    They exist so algebra and connectivity can be unit-tested.

    State order:
        0: R, 1: AR, 2: A2R, 3: A2F, 4: A2O, 5: A2D
    Open channel:
        State index 4 (A2O). Note this differs from 5-state JW95 where open_index is 3.
        Methods and analysis requiring open_index must use open_index = 4.
    """
    kon: float = 0.011244       # uM^-1 ms^-1  PROVISIONAL
    koff: float = 0.333069      # ms^-1        PROVISIONAL
    delta: float = 5.0          # ms^-1 (A2R -> A2F flipping forward) PROVISIONAL
    gamma: float = 1.0          # ms^-1 (A2F -> A2R flipping backward) PROVISIONAL
    beta: float = 0.559023      # ms^-1 (A2F -> A2O channel opening) PROVISIONAL
    alpha: float = 0.107891     # ms^-1 (A2O -> A2F channel closing) PROVISIONAL
    d: float = 0.050            # ms^-1 (A2R -> A2D desensitisation entry) PROVISIONAL
    r: float = 0.0020           # ms^-1 (A2D -> A2R resensitisation recovery) PROVISIONAL

    @property
    def name(self) -> str:
        return "KineticFlip"

    @property
    def open_index(self) -> int:
        return OPEN_IDX_6

    @property
    def native_observable(self) -> Observable:
        """`dose_response` is the stationary distribution, desensitisation included."""
        return Observable.EQUILIBRIUM

    @property
    def param_names(self) -> tuple[str, ...]:
        return ("kon", "koff", "delta", "gamma", "beta", "alpha", "d", "r")

    def get_params(self) -> dict[str, float]:
        return {
            "kon": self.kon,
            "koff": self.koff,
            "delta": self.delta,
            "gamma": self.gamma,
            "beta": self.beta,
            "alpha": self.alpha,
            "d": self.d,
            "r": self.r,
        }

    def with_params(self, **kwargs: float) -> "KineticFlipModel":
        return replace(self, **kwargs)

    @property
    def kd_um(self) -> float:
        """Microscopic dissociation constant Kd = koff / kon in uM."""
        return self.koff / self.kon

    @property
    def flipping_ratio(self) -> float:
        """Flipping equilibrium constant F = delta / gamma."""
        return self.delta / self.gamma

    @property
    def gating_efficacy(self) -> float:
        """Gating efficacy E = beta / alpha."""
        return self.beta / self.alpha

    @property
    def desens_ratio(self) -> float:
        """Desensitisation equilibrium ratio D = d / r."""
        return self.d / self.r

    @property
    def po_max(self) -> float:
        """Theoretical maximum instantaneous open probability at saturating GABA before desensitisation.

        In the flip scheme at infinite agonist before desensitisation, equilibrium between
        A2R, A2F, and A2O satisfies P(A2F) = F P(A2R) and P(A2O) = E P(A2F) = E F P(A2R).
        Normalising over {A2R, A2F, A2O}:
            po_max = E * F / (1 + F + E * F) = (F * E) / (1 + F * (1 + E)).
        As delta -> inf (or F -> inf), this collapses to E / (1 + E) = beta / (beta + alpha),
        recovering the JW95 po_max exactly.
        """
        f = self.flipping_ratio
        e = self.gating_efficacy
        return (f * e) / (1.0 + f * (1.0 + e))

    @property
    def po_inf(self) -> float:
        """Steady-state open probability at infinite GABA after desensitisation equilibrium."""
        f = self.flipping_ratio
        e = self.gating_efficacy
        d = self.desens_ratio
        return (f * e) / (1.0 + d + f * (1.0 + e))

    def apply_pam(
        self,
        affinity_factor: float = 1.0,
        gating_factor: float = 1.0,
        *,
        pam_factor: float | None = None,
    ) -> "KineticFlipModel":
        """Apply positive allosteric modulation (affinity shifts koff, gating shifts beta)."""
        if pam_factor is not None:
            affinity_factor = pam_factor
        for nm, v in (("affinity_factor", affinity_factor),
                      ("gating_factor", gating_factor)):
            v = float(v)
            if not np.isfinite(v):
                raise ValueError(f"{nm} must be finite, got {v!r}")
            if v < 1.0:
                raise ValueError(
                    f"{nm}={v!r} is below 1.0, i.e. NEGATIVE allosteric modulation, which "
                    f"this scheme's callers are not written for. Pass >= 1.0."
                )
        return replace(
            self,
            koff=self.koff / float(affinity_factor),
            beta=self.beta * float(gating_factor),
        )

    def q_matrix(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        """Row-generator transition matrix Q: dP/dt = P Q.

        States: 0: R, 1: AR, 2: A2R, 3: A2F, 4: A2O, 5: A2D.
        """
        mod = self.apply_pam(affinity_factor=pam_factor)
        g = max(float(gaba_um), 0.0)
        kon = mod.kon
        koff = mod.koff
        delta = mod.delta
        gamma = mod.gamma
        beta = mod.beta
        alpha = mod.alpha
        d = mod.d
        r = mod.r

        q = np.array([
            [-2.0 * kon * g,  2.0 * kon * g,        0.0,                    0.0,             0.0,    0.0],
            [koff,            -(koff + kon * g),     kon * g,                0.0,             0.0,    0.0],
            [0.0,             2.0 * koff,           -(2.0*koff+delta+d),     delta,           0.0,    d  ],
            [0.0,             0.0,                  gamma,                  -(gamma + beta), beta,   0.0],
            [0.0,             0.0,                  0.0,                    alpha,           -alpha, 0.0],
            [0.0,             0.0,                  r,                      0.0,             0.0,    -r ],
        ], dtype=float)
        return q

    def state_distribution(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        """Equilibrium probability vector P = [P_R, P_AR, P_A2R, P_A2F, P_A2O, P_A2D]."""
        mod = self.apply_pam(affinity_factor=pam_factor)
        g = max(float(gaba_um), 0.0)
        x = g / mod.kd_um
        f = mod.flipping_ratio
        e = mod.gating_efficacy
        d = mod.desens_ratio

        denom = 1.0 + 2.0 * x + (x ** 2) * (1.0 + d + f * (1.0 + e))
        if denom <= 0.0:
            return np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        p_r = 1.0 / denom
        p_ar = (2.0 * x) / denom
        p_a2r = (x ** 2) / denom
        p_a2f = (f * (x ** 2)) / denom
        p_a2o = (e * f * (x ** 2)) / denom
        p_a2d = (d * (x ** 2)) / denom

        return np.array([p_r, p_ar, p_a2r, p_a2f, p_a2o, p_a2d], dtype=float)

    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float:
        """Stationary open probability P(A2O)."""
        dist = self.state_distribution(gaba_um, pam_factor)
        return float(dist[OPEN_IDX_6])

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
        f = mod.flipping_ratio
        e = mod.gating_efficacy
        d = mod.desens_ratio
        denom = 1.0 + 2.0 * x + (x ** 2) * (1.0 + d + f * (1.0 + e))
        return (e * f * (x ** 2)) / denom

    def peak_dose_response(
        self,
        concs_um: np.ndarray,
        pam_factor: float = 1.0,
        application_ms: float = PEAK_APPLICATION_MS,
    ) -> np.ndarray:
        """PEAK open probability during a square agonist application."""
        return peak_open_probability_constant(
            self.q_matrix, self.state_distribution(0.0), concs_um,
            pam_factor=pam_factor, application_ms=application_ms,
            open_index=OPEN_IDX_6,
        )

    def simulate_waveform(
        self,
        t: np.ndarray,
        gaba_t: np.ndarray,
        pam_factor: float = 1.0,
        initial_state: np.ndarray | None = None,
    ) -> WaveformResult:
        """Integrate 6-state Markov scheme across time-varying GABA transient."""
        t_arr = np.asarray(t, dtype=float)
        gaba_arr = np.asarray(gaba_t, dtype=float)

        if len(t_arr) < 2:
            dist = (self.state_distribution(float(gaba_arr[0]), pam_factor)
                    if len(gaba_arr) > 0 else np.array([1.0, 0, 0, 0, 0, 0]))
            p0 = float(dist[OPEN_IDX_6])
            return WaveformResult(
                t=t_arr, gaba=gaba_arr, p_open=np.full_like(t_arr, p0), peak_p_open=p0,
                charge_integral=0.0, decay_tau_ms=float("nan"), steady_state_tail=p0,
                states=dist[None, :], decay=DecayFit.failed(),
            )

        p_init = resolve_initial_state(
            initial_state, 6, self.state_distribution(float(gaba_arr[0]), pam_factor),
            open_index=OPEN_IDX_6,
        )

        mod = self.apply_pam(affinity_factor=pam_factor)
        kon = mod.kon
        koff = mod.koff
        delta = mod.delta
        gamma = mod.gamma
        beta = mod.beta
        alpha = mod.alpha
        d_rate = mod.d
        r_rate = mod.r

        def rhs(time_val: float, y: np.ndarray) -> np.ndarray:
            g = float(np.interp(time_val, t_arr, gaba_arr))
            dp0 = -2.0 * kon * g * y[0] + koff * y[1]
            dp1 = 2.0 * kon * g * y[0] - (koff + kon * g) * y[1] + 2.0 * koff * y[2]
            dp2 = kon * g * y[1] - (2.0 * koff + delta + d_rate) * y[2] + gamma * y[3] + r_rate * y[5]
            dp3 = delta * y[2] - (gamma + beta) * y[3] + alpha * y[4]
            dp4 = beta * y[3] - alpha * y[4]
            dp5 = d_rate * y[2] - r_rate * y[5]
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

        if not sol.success:
            raise RuntimeError(
                f"Radau failed to integrate the 6-state scheme over "
                f"[{t_arr[0]:.4g}, {t_arr[-1]:.4g}] ms: {sol.message}. The stiffest "
                f"timescale here is 1/r = {1.0 / self.r:.0f} ms against a cleft transient "
                f"of ~1 ms, so this is a genuinely stiff problem; if the parameters are "
                f"legitimate, increase the density of `t` rather than loosening rtol."
            )

        states = sol.y.T
        row_sums = states.sum(axis=1)
        worst = float(np.max(np.abs(row_sums - 1.0))) if row_sums.size else 0.0
        if worst > 1e-3:
            raise RuntimeError(
                f"state probabilities drifted off the simplex by {worst:.3e} (max |sum - 1|) "
                f"while integrating the 6-state scheme. Renormalising would hide an "
                f"integration failure behind a plausible trace; tighten rtol/atol or "
                f"densify `t` instead."
            )
        states = np.clip(states, 0.0, 1.0)
        states = states / np.maximum(states.sum(axis=1, keepdims=True), 1e-12)

        p_open = states[:, OPEN_IDX_6]
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
