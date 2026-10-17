"""6-State Kinetic Scheme with Cooperative Binding and Two Sequential Open States.

Agent Brief 13 (Task 1):
Investigates whether adding binding cooperativity to the 6-state two-open scheme
(kinetic_two_open) relieves the Hill slope conflict (nH = 0.66 predicted vs 2.2 +/- 0.4
published, Jahn et al. 1997) while preserving the deactivation timescale, recovery,
and Po,max gating band.

Topology:
    R  <-->  AR  <-->  A2R  <-->  A2O (open)  <-->  A2O' (open)
                         |
                       A2D (desensitised)

States and Order:
    0: "R"     -- unliganded closed receptor
    1: "AR"    -- singly-liganded closed receptor
    2: "A2R"   -- doubly-liganded resting shut receptor
    3: "A2O"   -- primary doubly-liganded open conducting channel
    4: "A2O'"  -- secondary doubly-liganded open conducting channel
    5: "A2D"   -- doubly-liganded desensitised shut receptor (connects strictly to A2R)

Binding Cooperativity:
    The first binding step has microscopic association rate kon and dissociation rate koff.
    With statistical factor 2 for two identical unliganded sites, R -> AR is 2*kon*g and
    AR -> R is koff.
    The second binding step has association rate kon and its own microscopic dissociation rate koff2.
    AR -> A2R is kon*g and A2R -> AR is 2*koff2.
    The cooperativity factor is defined as:
        c = koff / koff2
    * c > 1: POSITIVE COOPERATIVITY -- the second ligand dissociates more slowly than the first
             (koff2 < koff), meaning the receptor exhibits higher affinity for the second binding step.
             This steepens the concentration-response curve (increases Hill slope).
    * c < 1: NEGATIVE COOPERATIVITY -- the second ligand dissociates more rapidly (koff2 > koff).
    * c = 1: NON-COOPERATIVE BINDING -- identical, independent sites (koff2 = koff). At c = 1,
             the scheme reduces EXACTLY to KineticTwoOpenModel.

Conductance & Observable:
    Both A2O (state 3) and A2O' (state 4) are open, conducting states with identical
    unit conductance. The physiological observable (current / open probability) is the
    SUM of both open-state occupancies:
        P_open(t) = P(A2O)(t) + P(A2O')(t)

Microscopic Reversibility:
    The scheme's transition graph is an acyclic tree (6 vertices, 5 edges, loop-free).
    Because the network contains no closed cycles (cycle space dimension E - V + 1 = 5 - 6 + 1 = 0),
    Kolmogorov's loop product criterion is trivially satisfied. Detailed balance (microscopic
    reversibility) holds unconditionally for any set of strictly positive rate constants;
    there are no algebraic cycle constraints beyond positivity.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Sequence

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

STATES_6 = ("R", "AR", "A2R", "A2O", "A2O'", "A2D")
OPEN_IDX_1 = 3  # "A2O"
OPEN_IDX_2 = 4  # "A2O'"
OPEN_INDICES_6 = (OPEN_IDX_1, OPEN_IDX_2)


@dataclass(frozen=True)
class KineticCoopModel(ReceptorModel):
    """6-state Markov gating scheme with binding cooperativity and two open conducting states.

    Topology:
        R <--> AR <--> A2R <--> A2O <--> A2O'
                         |
                       A2D

    Both A2O (state 3) and A2O' (state 4) are open conducting states. The observable
    is the sum of both occupancies: P_open = P(A2O) + P(A2O').
    The second dissociation rate koff2 defines cooperativity c = koff / koff2.

    EVERY DEFAULT BELOW IS PROVISIONAL AND IS NOT A CALIBRATION.
    They exist so algebra and connectivity can be unit-tested.
    """
    kon: float = 0.011244       # uM^-1 ms^-1  PROVISIONAL
    koff: float = 0.333069      # ms^-1 (first unbinding AR -> R) PROVISIONAL
    koff2: float = 0.333069     # ms^-1 (second unbinding A2R -> AR) PROVISIONAL (default c = 1.0)
    beta: float = 0.559023      # ms^-1 (A2R -> A2O opening) PROVISIONAL
    alpha: float = 0.107891     # ms^-1 (A2O -> A2R closing) PROVISIONAL
    delta: float = 1.0          # ms^-1 (A2O -> A2O' second opening) PROVISIONAL
    gamma: float = 1.0          # ms^-1 (A2O' -> A2O second closing) PROVISIONAL
    d: float = 0.050            # ms^-1 (A2R -> A2D desensitisation entry) PROVISIONAL
    r: float = 0.0020           # ms^-1 (A2D -> A2R resensitisation recovery) PROVISIONAL

    @property
    def name(self) -> str:
        return "KineticCoop"

    @property
    def open_indices(self) -> tuple[int, ...]:
        """Indices of all conducting open states: (3, 4) for (A2O, A2O')."""
        return OPEN_INDICES_6

    @property
    def open_index(self) -> tuple[int, ...]:
        """Open state indices for the two-open conducting scheme (A2O and A2O')."""
        return OPEN_INDICES_6

    @property
    def native_observable(self) -> Observable:
        """`dose_response` is the stationary distribution, desensitisation included."""
        return Observable.EQUILIBRIUM

    @property
    def param_names(self) -> tuple[str, ...]:
        return ("kon", "koff", "koff2", "beta", "alpha", "delta", "gamma", "d", "r")

    @property
    def beta2(self) -> float:
        """Alias for delta (forward rate A2O -> A2O')."""
        return self.delta

    @property
    def alpha2(self) -> float:
        """Alias for gamma (reverse rate A2O' -> A2O)."""
        return self.gamma

    def get_params(self) -> dict[str, float]:
        return {
            "kon": self.kon,
            "koff": self.koff,
            "koff2": self.koff2,
            "beta": self.beta,
            "alpha": self.alpha,
            "delta": self.delta,
            "gamma": self.gamma,
            "d": self.d,
            "r": self.r,
        }

    def with_params(self, **kwargs: float) -> "KineticCoopModel":
        if "c" in kwargs and "koff2" not in kwargs:
            c_val = float(kwargs.pop("c"))
            koff_val = float(kwargs.get("koff", self.koff))
            kwargs["koff2"] = koff_val / c_val
        return replace(self, **kwargs)

    @property
    def kd_um(self) -> float:
        """First microscopic dissociation constant Kd1 = koff / kon in uM."""
        return self.koff / self.kon

    @property
    def kd2_um(self) -> float:
        """Second microscopic dissociation constant Kd2 = koff2 / kon in uM."""
        return self.koff2 / self.kon

    @property
    def c(self) -> float:
        """Binding cooperativity factor c = koff / koff2.

        c > 1: positive cooperativity (koff2 < koff, tighter second binding).
        c < 1: negative cooperativity (koff2 > koff).
        c = 1: non-cooperative binding (reduces exactly to KineticTwoOpenModel).
        """
        return self.koff / self.koff2

    @property
    def gating_efficacy(self) -> float:
        """First gating efficacy E = beta / alpha."""
        return self.beta / self.alpha

    @property
    def two_open_ratio(self) -> float:
        """Equilibrium ratio between open states K_open2 = delta / gamma."""
        return self.delta / self.gamma

    @property
    def desens_ratio(self) -> float:
        """Desensitisation equilibrium ratio D = d / r."""
        return self.d / self.r

    @property
    def po_max(self) -> float:
        """Theoretical maximum instantaneous open probability at saturating GABA before desensitisation.

        In the two-open scheme at infinite agonist before desensitisation, equilibrium between
        A2R, A2O, and A2O' satisfies P(A2O) = E P(A2R) and P(A2O') = K_open2 P(A2O) = E K_open2 P(A2R).
        Both A2O and A2O' conduct. Normalising over {A2R, A2O, A2O'}:
            po_max = (P(A2O) + P(A2O')) / (P(A2R) + P(A2O) + P(A2O'))
                   = (E * (1 + K_open2)) / (1 + E * (1 + K_open2)).
        As delta -> 0 (or K_open2 -> 0), this collapses to E / (1 + E) = beta / (beta + alpha),
        recovering the JW95 po_max exactly. Unchanged by binding cooperativity because agonist is saturating.
        """
        e = self.gating_efficacy
        k2 = self.two_open_ratio
        eff_e = e * (1.0 + k2)
        return eff_e / (1.0 + eff_e)

    @property
    def po_inf(self) -> float:
        """Steady-state open probability at infinite GABA after desensitisation equilibrium."""
        e = self.gating_efficacy
        k2 = self.two_open_ratio
        d = self.desens_ratio
        eff_e = e * (1.0 + k2)
        return eff_e / (1.0 + d + eff_e)

    def apply_pam(
        self,
        affinity_factor: float = 1.0,
        gating_factor: float = 1.0,
        *,
        pam_factor: float | None = None,
    ) -> "KineticCoopModel":
        """Apply positive allosteric modulation (affinity shifts koff and koff2, gating shifts beta)."""
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
            koff2=self.koff2 / float(affinity_factor),
            beta=self.beta * float(gating_factor),
        )

    def q_matrix(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        """Row-generator transition matrix Q: dP/dt = P Q.

        States: 0: R, 1: AR, 2: A2R, 3: A2O, 4: A2O', 5: A2D.
        """
        mod = self.apply_pam(affinity_factor=pam_factor)
        g = max(float(gaba_um), 0.0)
        kon = mod.kon
        koff = mod.koff
        koff2 = mod.koff2
        beta = mod.beta
        alpha = mod.alpha
        delta = mod.delta
        gamma = mod.gamma
        d = mod.d
        r = mod.r

        q = np.array([
            [-2.0 * kon * g,  2.0 * kon * g,        0.0,                    0.0,             0.0,    0.0],
            [koff,            -(koff + kon * g),     kon * g,                0.0,             0.0,    0.0],
            [0.0,             2.0 * koff2,          -(2.0 * koff2 + beta + d), beta,         0.0,    d  ],
            [0.0,             0.0,                  alpha,                  -(alpha + delta), delta, 0.0],
            [0.0,             0.0,                  0.0,                    gamma,           -gamma, 0.0],
            [0.0,             0.0,                  r,                      0.0,             0.0,    -r ],
        ], dtype=float)
        return q

    def state_distribution(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        """Equilibrium probability vector P = [P_R, P_AR, P_A2R, P_A2O, P_A2O', P_A2D].

        Detailed balance on tree:
            P(AR)  = 2 * (kon * g / koff) * P(R) = 2 * x * P(R)
            P(A2R) = (kon * g / (2 * koff2)) * P(AR) = (kon * g)^2 / (koff * koff2) * P(R)
                   = (koff / koff2) * (kon * g / koff)^2 * P(R) = c * x^2 * P(R)
            P(A2O)  = E * P(A2R) = c * E * x^2 * P(R)
            P(A2O') = K_open2 * P(A2O) = c * K_open2 * E * x^2 * P(R)
            P(A2D)  = D * P(A2R) = c * D * x^2 * P(R)
        """
        mod = self.apply_pam(affinity_factor=pam_factor)
        g = max(float(gaba_um), 0.0)
        x = g / mod.kd_um
        c = mod.c
        e = mod.gating_efficacy
        k2 = mod.two_open_ratio
        d = mod.desens_ratio

        denom = 1.0 + 2.0 * x + c * (x ** 2) * (1.0 + d + e * (1.0 + k2))
        if denom <= 0.0:
            return np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        p_r = 1.0 / denom
        p_ar = (2.0 * x) / denom
        p_a2r = (c * (x ** 2)) / denom
        p_a2o = (c * e * (x ** 2)) / denom
        p_a2o_prime = (c * k2 * e * (x ** 2)) / denom
        p_a2d = (c * d * (x ** 2)) / denom

        return np.array([p_r, p_ar, p_a2r, p_a2o, p_a2o_prime, p_a2d], dtype=float)

    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float:
        """Stationary total open probability P(A2O) + P(A2O')."""
        dist = self.state_distribution(gaba_um, pam_factor)
        return float(dist[OPEN_IDX_1] + dist[OPEN_IDX_2])

    def dose_response(
        self,
        concs_um: np.ndarray,
        pam_factor: float = 1.0,
        *,
        observable: Observable = Observable.EQUILIBRIUM,
    ) -> np.ndarray:
        """EQUILIBRIUM total open probability P(A2O) + P(A2O')."""
        if observable is not Observable.EQUILIBRIUM:
            raise ObservableMismatch(
                f"{self.name} is EQUILIBRIUM-native and has no {observable.value} curve in "
                f"dose_response. Use peak_dose_response for PEAK data."
            )
        concs = np.asarray(concs_um, dtype=float)
        mod = self.apply_pam(affinity_factor=pam_factor)
        x = np.maximum(concs, 0.0) / mod.kd_um
        c = mod.c
        e = mod.gating_efficacy
        k2 = mod.two_open_ratio
        d = mod.desens_ratio
        denom = 1.0 + 2.0 * x + c * (x ** 2) * (1.0 + d + e * (1.0 + k2))
        return (c * e * (1.0 + k2) * (x ** 2)) / denom

    def peak_dose_response(
        self,
        concs_um: np.ndarray,
        pam_factor: float = 1.0,
        application_ms: float = PEAK_APPLICATION_MS,
    ) -> np.ndarray:
        """PEAK total open probability P(A2O) + P(A2O') during a square agonist application."""
        return peak_open_probability_constant(
            self.q_matrix,
            self.state_distribution(0.0),
            concs_um,
            pam_factor=pam_factor,
            application_ms=application_ms,
            open_index=OPEN_INDICES_6,
        )

    def simulate_waveform(
        self,
        t: np.ndarray,
        gaba_t: np.ndarray,
        pam_factor: float = 1.0,
        initial_state: np.ndarray | None = None,
    ) -> WaveformResult:
        """Integrate 6-state cooperative two-open scheme across time-varying GABA transient.

        The resulting p_open is the SUM of the two conducting states:
            p_open = states[:, 3] + states[:, 4]
        """
        t_arr = np.asarray(t, dtype=float)
        gaba_arr = np.asarray(gaba_t, dtype=float)

        if len(t_arr) < 2:
            dist = (self.state_distribution(float(gaba_arr[0]), pam_factor)
                    if len(gaba_arr) > 0 else np.array([1.0, 0, 0, 0, 0, 0]))
            p0 = float(dist[OPEN_IDX_1] + dist[OPEN_IDX_2])
            return WaveformResult(
                t=t_arr, gaba=gaba_arr, p_open=np.full_like(t_arr, p0), peak_p_open=p0,
                charge_integral=0.0, decay_tau_ms=float("nan"), steady_state_tail=p0,
                states=dist[None, :], decay=DecayFit.failed(),
            )

        p_init = resolve_initial_state(
            initial_state, 6, self.state_distribution(float(gaba_arr[0]), pam_factor),
            open_index=OPEN_INDICES_6,
        )

        mod = self.apply_pam(affinity_factor=pam_factor)
        kon = mod.kon
        koff = mod.koff
        koff2 = mod.koff2
        beta = mod.beta
        alpha = mod.alpha
        delta = mod.delta
        gamma = mod.gamma
        d_rate = mod.d
        r_rate = mod.r

        def rhs(time_val: float, y: np.ndarray) -> np.ndarray:
            g = float(np.interp(time_val, t_arr, gaba_arr))
            dp0 = -2.0 * kon * g * y[0] + koff * y[1]
            dp1 = 2.0 * kon * g * y[0] - (koff + kon * g) * y[1] + 2.0 * koff2 * y[2]
            dp2 = kon * g * y[1] - (2.0 * koff2 + beta + d_rate) * y[2] + alpha * y[3] + r_rate * y[5]
            dp3 = beta * y[2] - (alpha + delta) * y[3] + gamma * y[4]
            dp4 = delta * y[3] - gamma * y[4]
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
                f"Radau failed to integrate the 6-state cooperative scheme over "
                f"[{t_arr[0]:.4g}, {t_arr[-1]:.4g}] ms: {sol.message}."
            )

        states = sol.y.T
        row_sums = states.sum(axis=1)
        worst = float(np.max(np.abs(row_sums - 1.0))) if row_sums.size else 0.0
        if worst > 1e-3:
            raise RuntimeError(
                f"state probabilities drifted off the simplex by {worst:.3e} (max |sum - 1|) "
                f"while integrating the 6-state cooperative scheme. Tighten rtol/atol or densify `t`."
            )
        states = np.clip(states, 0.0, 1.0)
        states = states / np.maximum(states.sum(axis=1, keepdims=True), 1e-12)

        # Total open probability is the sum of both conducting states
        p_open = states[:, OPEN_IDX_1] + states[:, OPEN_IDX_2]
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
