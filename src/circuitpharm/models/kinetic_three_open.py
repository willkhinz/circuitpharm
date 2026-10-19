"""7-State Kinetic Schemes with Cooperative Binding and Three Open Conducting States.

Agent Brief 18 (Task 1):
Investigates whether adding a third conducting open state to the cooperative scheme
(kinetic_coop), as independently required by Keramidas & Harrison 2008 (PMC2213567)
for alpha1beta2gamma2S bursts (three shut and three open components in M-mode,
two shut and three open states in H-mode), relieves the residual Hill slope ceiling
(nH <= 1.25 against published 2.2 +/- 0.4, Jahn et al. 1997).

Topology Choices and Justification:
-----------------------------------
The literature establishes that bursts require three open components, but does not
uniquely determine the kinetic topology of the third open state. Two primary candidate
topologies exist:

1. SEQUENTIAL (SERIES) TOPOLOGY (KineticThreeOpenModel - Primary):
       R  <-->  AR  <-->  A2R  <-->  A2O  <-->  A2O'  <-->  A2O''
                            |
                          A2D
   - Justification:
     In single-channel burst kinetics (Colquhoun & Hawkes 1981, 1982; Jones & Westbrook 1995),
     sequential open states reflect progressive conformational relaxation of the channel pore
     along a single activation coordinate once unconstrained by shut gates. It is the direct
     natural extension of the project's own progression: kinetic_jw95 (1 open) ->
     kinetic_two_open (2 open) -> kinetic_three_open (3 open).
   - Gateway Structure:
     The open conducting manifold {A2O, A2O', A2O''} has a SINGLE gateway transition to shut
     states: A2R <--> A2O (rates beta, alpha).
   - Physical prediction:
     Channel closing to the resting/shut manifold can ONLY occur from A2O. Excursions into
     A2O' and A2O'' dwell within the conducting manifold and prolong open bursts without
     closing. Entry flux from shut states is 100% concentrated in A2O: pi_O = (1, 0, 0).
     The exact mean open dwell time before closing is:
         <tau_O> = (1 / alpha) * (1 + delta/gamma + (delta/gamma)*(delta2/gamma2)).

2. PARALLEL TOPOLOGY (KineticThreeOpenParallelModel - Alternative):
       R  <-->  AR  <-->  A2R  <-->  A2O  <-->  A2O'
                            |          ^
                          A2D          | (beta_p, alpha_p)
                                       v
                                     A2O''
   - Rationale:
     GABA-A receptors possess two non-equivalent agonist binding interfaces (alpha-beta and
     beta-alpha). Alternative conformations could open directly from the doubly-bound shut
     state A2R into distinct conducting states (e.g., subconductance or alternate gating modes).
   - Gateway Structure:
     TWO independent opening pathways from shut state A2R: A2R <--> A2O (beta, alpha) and
     A2R <--> A2O'' (beta_p, alpha_p).
   - What the alternative predicts differently:
     1. Openings into A2O'' can terminate directly back to A2R at rate alpha_p without passing
        through A2O.
     2. Entry flux from shut states is partitioned between A2O and A2O'':
        pi_O = (beta / (beta + beta_p), 0, beta_p / (beta + beta_p)).
     3. In deactivation, A2O'' adds an independent relaxation mode coupled directly to A2R,
        altering macroscopic deactivation amplitudes rather than prolonging intraburst dwells.

Microscopic Reversibility:
--------------------------
Both the sequential and parallel transition graphs are acyclic trees:
- Vertices V = 7 states
- Edges E = 6 reversible transitions
- Cycle space dimension = E - V + 1 = 6 - 7 + 1 = 0.
Because the transition network is strictly loop-free (no closed cycles), Kolmogorov's
loop product criterion is satisfied unconditionally. Detailed balance (microscopic
reversibility) holds for any strictly positive rate constants without algebraic constraints.

Conductance & Observable:
-------------------------
All three open states A2O (3), A2O' (4), and A2O'' (5) are conducting with identical
unit conductance. The physiological observable is the sum of all three open occupancies:
    P_open(t) = P(A2O)(t) + P(A2O')(t) + P(A2O'')(t).
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

STATES_7 = ("R", "AR", "A2R", "A2O", "A2O'", "A2O''", "A2D")
OPEN_IDX_1 = 3  # "A2O"
OPEN_IDX_2 = 4  # "A2O'"
OPEN_IDX_3 = 5  # "A2O''"
OPEN_INDICES_7 = (OPEN_IDX_1, OPEN_IDX_2, OPEN_IDX_3)


@dataclass(frozen=True)
class KineticThreeOpenModel(ReceptorModel):
    """7-state Markov gating scheme with binding cooperativity and three sequential open states.

    Topology:
        R <--> AR <--> A2R <--> A2O <--> A2O' <--> A2O''
                         |
                       A2D

    States:
        0: R     -- unliganded closed receptor
        1: AR    -- singly-liganded closed receptor
        2: A2R   -- doubly-liganded resting shut receptor
        3: A2O   -- primary doubly-liganded open conducting channel
        4: A2O'  -- secondary doubly-liganded open conducting channel
        5: A2O'' -- tertiary doubly-liganded open conducting channel
        6: A2D   -- doubly-liganded desensitised shut receptor (connects strictly to A2R)

    Conductance:
        All three states A2O, A2O', and A2O'' conduct.
        P_open = P(A2O) + P(A2O') + P(A2O'').

    Nesting:
        In the limit delta2 -> 0, transitions to A2O'' cease, and the scheme reduces
        EXACTLY to KineticCoopModel.

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
    delta2: float = 1.0         # ms^-1 (A2O' -> A2O'' third opening) PROVISIONAL
    gamma2: float = 1.0         # ms^-1 (A2O'' -> A2O' third closing) PROVISIONAL
    d: float = 0.050            # ms^-1 (A2R -> A2D desensitisation entry) PROVISIONAL
    r: float = 0.0020           # ms^-1 (A2D -> A2R resensitisation recovery) PROVISIONAL

    @property
    def name(self) -> str:
        return "KineticThreeOpen"

    @property
    def open_indices(self) -> tuple[int, ...]:
        """Indices of all conducting open states: (3, 4, 5) for (A2O, A2O', A2O'')."""
        return OPEN_INDICES_7

    @property
    def open_index(self) -> tuple[int, ...]:
        """Open state indices for the three-open conducting scheme."""
        return OPEN_INDICES_7

    @property
    def native_observable(self) -> Observable:
        """`dose_response` is the stationary distribution, desensitisation included."""
        return Observable.EQUILIBRIUM

    @property
    def param_names(self) -> tuple[str, ...]:
        return (
            "kon", "koff", "koff2", "beta", "alpha",
            "delta", "gamma", "delta2", "gamma2", "d", "r",
        )

    @property
    def beta2(self) -> float:
        """Alias for delta (forward rate A2O -> A2O')."""
        return self.delta

    @property
    def alpha2(self) -> float:
        """Alias for gamma (reverse rate A2O' -> A2O)."""
        return self.gamma

    @property
    def beta3(self) -> float:
        """Alias for delta2 (forward rate A2O' -> A2O'')."""
        return self.delta2

    @property
    def alpha3(self) -> float:
        """Alias for gamma2 (reverse rate A2O'' -> A2O')."""
        return self.gamma2

    def get_params(self) -> dict[str, float]:
        return {
            "kon": self.kon,
            "koff": self.koff,
            "koff2": self.koff2,
            "beta": self.beta,
            "alpha": self.alpha,
            "delta": self.delta,
            "gamma": self.gamma,
            "delta2": self.delta2,
            "gamma2": self.gamma2,
            "d": self.d,
            "r": self.r,
        }

    def with_params(self, **kwargs: float) -> "KineticThreeOpenModel":
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
        """Binding cooperativity factor c = koff / koff2."""
        return self.koff / self.koff2

    @property
    def gating_efficacy(self) -> float:
        """First gating efficacy E = beta / alpha."""
        return self.beta / self.alpha

    @property
    def two_open_ratio(self) -> float:
        """Equilibrium ratio between second and first open states K_open2 = delta / gamma."""
        return self.delta / self.gamma

    @property
    def three_open_ratio(self) -> float:
        """Equilibrium ratio between third and second open states K_open3 = delta2 / gamma2."""
        return (self.delta2 / self.gamma2) if self.gamma2 > 0.0 else 0.0

    @property
    def desens_ratio(self) -> float:
        """Desensitisation equilibrium ratio D = d / r."""
        return self.d / self.r

    @property
    def po_max(self) -> float:
        """Theoretical maximum instantaneous open probability at saturating GABA before desensitisation.

        In the sequential 3-open scheme at saturating agonist before desensitisation:
            P(A2O)   = E * P(A2R)
            P(A2O')  = K_open2 * P(A2O)   = E * K_open2 * P(A2R)
            P(A2O'') = K_open3 * P(A2O')  = E * K_open2 * K_open3 * P(A2R)
        Total open occupancy:
            P_open = E * (1 + K_open2 + K_open2 * K_open3) * P(A2R).
        Normalising over {A2R, A2O, A2O', A2O''}:
            eff_e = E * (1 + K_open2 + K_open2 * K_open3)
            po_max = eff_e / (1 + eff_e).
        As delta2 -> 0, this continuously reduces to E * (1 + K_open2) / (1 + E * (1 + K_open2)),
        recovering KineticCoopModel exactly.
        """
        e = self.gating_efficacy
        k2 = self.two_open_ratio
        k3 = self.three_open_ratio
        eff_e = e * (1.0 + k2 + k2 * k3)
        return eff_e / (1.0 + eff_e)

    @property
    def po_inf(self) -> float:
        """Steady-state open probability at infinite GABA after desensitisation equilibrium."""
        e = self.gating_efficacy
        k2 = self.two_open_ratio
        k3 = self.three_open_ratio
        d = self.desens_ratio
        eff_e = e * (1.0 + k2 + k2 * k3)
        return eff_e / (1.0 + d + eff_e)

    def apply_pam(
        self,
        affinity_factor: float = 1.0,
        gating_factor: float = 1.0,
        *,
        pam_factor: float | None = None,
    ) -> "KineticThreeOpenModel":
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

        States: 0: R, 1: AR, 2: A2R, 3: A2O, 4: A2O', 5: A2O'', 6: A2D.
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
        delta2 = mod.delta2
        gamma2 = mod.gamma2
        d = mod.d
        r = mod.r

        q = np.array([
            [-2.0 * kon * g,  2.0 * kon * g,        0.0,                        0.0,             0.0,               0.0,     0.0],
            [koff,            -(koff + kon * g),     kon * g,                    0.0,             0.0,               0.0,     0.0],
            [0.0,             2.0 * koff2,          -(2.0 * koff2 + beta + d),   beta,            0.0,               0.0,     d  ],
            [0.0,             0.0,                  alpha,                      -(alpha + delta), delta,             0.0,     0.0],
            [0.0,             0.0,                  0.0,                        gamma,           -(gamma + delta2),  delta2,  0.0],
            [0.0,             0.0,                  0.0,                        0.0,             gamma2,            -gamma2,  0.0],
            [0.0,             0.0,                  r,                          0.0,             0.0,               0.0,     -r ],
        ], dtype=float)
        return q

    def state_distribution(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        """Equilibrium probability vector P = [P_R, P_AR, P_A2R, P_A2O, P_A2O', P_A2O'', P_A2D].

        Detailed balance on tree:
            P(AR)   = 2 * x * P(R)
            P(A2R)  = c * x^2 * P(R)
            P(A2O)  = c * E * x^2 * P(R)
            P(A2O') = c * K_open2 * E * x^2 * P(R)
            P(A2O'')= c * K_open2 * K_open3 * E * x^2 * P(R)
            P(A2D)  = c * D * x^2 * P(R)
        """
        mod = self.apply_pam(affinity_factor=pam_factor)
        g = max(float(gaba_um), 0.0)
        x = g / mod.kd_um
        c = mod.c
        e = mod.gating_efficacy
        k2 = mod.two_open_ratio
        k3 = mod.three_open_ratio
        d = mod.desens_ratio

        eff_o = e * (1.0 + k2 + k2 * k3)
        denom = 1.0 + 2.0 * x + c * (x ** 2) * (1.0 + d + eff_o)
        if denom <= 0.0:
            return np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        p_r = 1.0 / denom
        p_ar = (2.0 * x) / denom
        p_a2r = (c * (x ** 2)) / denom
        p_a2o = (c * e * (x ** 2)) / denom
        p_a2o_prime = (c * k2 * e * (x ** 2)) / denom
        p_a2o_double_prime = (c * k2 * k3 * e * (x ** 2)) / denom
        p_a2d = (c * d * (x ** 2)) / denom

        return np.array([
            p_r, p_ar, p_a2r, p_a2o, p_a2o_prime, p_a2o_double_prime, p_a2d
        ], dtype=float)

    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float:
        """Stationary total open probability P(A2O) + P(A2O') + P(A2O'')."""
        dist = self.state_distribution(gaba_um, pam_factor)
        return float(dist[OPEN_IDX_1] + dist[OPEN_IDX_2] + dist[OPEN_IDX_3])

    def dose_response(
        self,
        concs_um: np.ndarray,
        pam_factor: float = 1.0,
        *,
        observable: Observable = Observable.EQUILIBRIUM,
    ) -> np.ndarray:
        """EQUILIBRIUM total open probability P(A2O) + P(A2O') + P(A2O'')."""
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
        k3 = mod.three_open_ratio
        d = mod.desens_ratio
        eff_o = e * (1.0 + k2 + k2 * k3)
        denom = 1.0 + 2.0 * x + c * (x ** 2) * (1.0 + d + eff_o)
        return (c * eff_o * (x ** 2)) / denom

    def peak_dose_response(
        self,
        concs_um: np.ndarray,
        pam_factor: float = 1.0,
        application_ms: float = PEAK_APPLICATION_MS,
    ) -> np.ndarray:
        """PEAK total open probability P(A2O) + P(A2O') + P(A2O'') during a square agonist application."""
        return peak_open_probability_constant(
            self.q_matrix,
            self.state_distribution(0.0),
            concs_um,
            pam_factor=pam_factor,
            application_ms=application_ms,
            open_index=OPEN_INDICES_7,
        )

    def simulate_waveform(
        self,
        t: np.ndarray,
        gaba_t: np.ndarray,
        pam_factor: float = 1.0,
        initial_state: np.ndarray | None = None,
    ) -> WaveformResult:
        """Integrate 7-state cooperative three-open scheme across time-varying GABA transient.

        The resulting p_open is the SUM of all three conducting states:
            p_open = states[:, 3] + states[:, 4] + states[:, 5]
        """
        t_arr = np.asarray(t, dtype=float)
        gaba_arr = np.asarray(gaba_t, dtype=float)

        if len(t_arr) < 2:
            dist = (self.state_distribution(float(gaba_arr[0]), pam_factor)
                    if len(gaba_arr) > 0 else np.array([1.0, 0, 0, 0, 0, 0, 0]))
            p0 = float(dist[OPEN_IDX_1] + dist[OPEN_IDX_2] + dist[OPEN_IDX_3])
            return WaveformResult(
                t=t_arr, gaba=gaba_arr, p_open=np.full_like(t_arr, p0), peak_p_open=p0,
                charge_integral=0.0, decay_tau_ms=float("nan"), steady_state_tail=p0,
                states=dist[None, :], decay=DecayFit.failed(),
            )

        p_init = resolve_initial_state(
            initial_state, 7, self.state_distribution(float(gaba_arr[0]), pam_factor),
            open_index=OPEN_INDICES_7,
        )

        mod = self.apply_pam(affinity_factor=pam_factor)
        kon = mod.kon
        koff = mod.koff
        koff2 = mod.koff2
        beta = mod.beta
        alpha = mod.alpha
        delta = mod.delta
        gamma = mod.gamma
        delta2 = mod.delta2
        gamma2 = mod.gamma2
        d_rate = mod.d
        r_rate = mod.r

        def rhs(time_val: float, y: np.ndarray) -> np.ndarray:
            g = float(np.interp(time_val, t_arr, gaba_arr))
            dp0 = -2.0 * kon * g * y[0] + koff * y[1]
            dp1 = 2.0 * kon * g * y[0] - (koff + kon * g) * y[1] + 2.0 * koff2 * y[2]
            dp2 = kon * g * y[1] - (2.0 * koff2 + beta + d_rate) * y[2] + alpha * y[3] + r_rate * y[6]
            dp3 = beta * y[2] - (alpha + delta) * y[3] + gamma * y[4]
            dp4 = delta * y[3] - (gamma + delta2) * y[4] + gamma2 * y[5]
            dp5 = delta2 * y[4] - gamma2 * y[5]
            dp6 = d_rate * y[2] - r_rate * y[6]
            return np.array([dp0, dp1, dp2, dp3, dp4, dp5, dp6])

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
                f"Radau failed to integrate the 7-state cooperative scheme over "
                f"[{t_arr[0]:.4g}, {t_arr[-1]:.4g}] ms: {sol.message}."
            )

        states = sol.y.T
        row_sums = states.sum(axis=1)
        worst = float(np.max(np.abs(row_sums - 1.0))) if row_sums.size else 0.0
        if worst > 1e-3:
            raise RuntimeError(
                f"state probabilities drifted off the simplex by {worst:.3e} (max |sum - 1|) "
                f"while integrating the 7-state scheme. Tighten rtol/atol or densify `t`."
            )
        states = np.clip(states, 0.0, 1.0)
        states = states / np.maximum(states.sum(axis=1, keepdims=True), 1e-12)

        # Total open probability is the sum of all three conducting states
        p_open = states[:, OPEN_IDX_1] + states[:, OPEN_IDX_2] + states[:, OPEN_IDX_3]
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


@dataclass(frozen=True)
class KineticThreeOpenParallelModel(ReceptorModel):
    """7-state Markov gating scheme with binding cooperativity and PARALLEL third open state.

    Topology:
        R <--> AR <--> A2R <--> A2O <--> A2O'
                         |        ^
                       A2D        | (beta_p, alpha_p)
                                  v
                                A2O''

    States:
        0: R     -- unliganded closed receptor
        1: AR    -- singly-liganded closed receptor
        2: A2R   -- doubly-liganded resting shut receptor
        3: A2O   -- primary doubly-liganded open conducting channel
        4: A2O'  -- secondary doubly-liganded open conducting channel
        5: A2O'' -- parallel doubly-liganded open conducting channel (connects to A2R)
        6: A2D   -- doubly-liganded desensitised shut receptor (connects strictly to A2R)

    Conductance:
        All three states A2O, A2O', and A2O'' conduct.
        P_open = P(A2O) + P(A2O') + P(A2O'').

    Nesting:
        In the limit beta_p -> 0, transitions to A2O'' cease, and the scheme reduces
        EXACTLY to KineticCoopModel.
    """
    kon: float = 0.011244       # uM^-1 ms^-1  PROVISIONAL
    koff: float = 0.333069      # ms^-1 (first unbinding AR -> R) PROVISIONAL
    koff2: float = 0.333069     # ms^-1 (second unbinding A2R -> AR) PROVISIONAL
    beta: float = 0.559023      # ms^-1 (A2R -> A2O opening) PROVISIONAL
    alpha: float = 0.107891     # ms^-1 (A2O -> A2R closing) PROVISIONAL
    delta: float = 1.0          # ms^-1 (A2O -> A2O' second opening) PROVISIONAL
    gamma: float = 1.0          # ms^-1 (A2O' -> A2O second closing) PROVISIONAL
    beta_p: float = 0.559023    # ms^-1 (A2R -> A2O'' parallel opening) PROVISIONAL
    alpha_p: float = 0.107891   # ms^-1 (A2O'' -> A2R parallel closing) PROVISIONAL
    d: float = 0.050            # ms^-1 (A2R -> A2D desensitisation entry) PROVISIONAL
    r: float = 0.0020           # ms^-1 (A2D -> A2R resensitisation recovery) PROVISIONAL

    @property
    def name(self) -> str:
        return "KineticThreeOpenParallel"

    @property
    def open_indices(self) -> tuple[int, ...]:
        return OPEN_INDICES_7

    @property
    def open_index(self) -> tuple[int, ...]:
        return OPEN_INDICES_7

    @property
    def native_observable(self) -> Observable:
        return Observable.EQUILIBRIUM

    @property
    def param_names(self) -> tuple[str, ...]:
        return (
            "kon", "koff", "koff2", "beta", "alpha",
            "delta", "gamma", "beta_p", "alpha_p", "d", "r",
        )

    def get_params(self) -> dict[str, float]:
        return {
            "kon": self.kon,
            "koff": self.koff,
            "koff2": self.koff2,
            "beta": self.beta,
            "alpha": self.alpha,
            "delta": self.delta,
            "gamma": self.gamma,
            "beta_p": self.beta_p,
            "alpha_p": self.alpha_p,
            "d": self.d,
            "r": self.r,
        }

    def with_params(self, **kwargs: float) -> "KineticThreeOpenParallelModel":
        if "c" in kwargs and "koff2" not in kwargs:
            c_val = float(kwargs.pop("c"))
            koff_val = float(kwargs.get("koff", self.koff))
            kwargs["koff2"] = koff_val / c_val
        return replace(self, **kwargs)

    @property
    def kd_um(self) -> float:
        return self.koff / self.kon

    @property
    def kd2_um(self) -> float:
        return self.koff2 / self.kon

    @property
    def c(self) -> float:
        return self.koff / self.koff2

    @property
    def gating_efficacy(self) -> float:
        return self.beta / self.alpha

    @property
    def parallel_gating_efficacy(self) -> float:
        return (self.beta_p / self.alpha_p) if self.alpha_p > 0.0 else 0.0

    @property
    def two_open_ratio(self) -> float:
        return self.delta / self.gamma

    @property
    def desens_ratio(self) -> float:
        return self.d / self.r

    @property
    def po_max(self) -> float:
        e = self.gating_efficacy
        k2 = self.two_open_ratio
        ep = self.parallel_gating_efficacy
        eff_e = e * (1.0 + k2) + ep
        return eff_e / (1.0 + eff_e)

    @property
    def po_inf(self) -> float:
        e = self.gating_efficacy
        k2 = self.two_open_ratio
        ep = self.parallel_gating_efficacy
        d = self.desens_ratio
        eff_e = e * (1.0 + k2) + ep
        return eff_e / (1.0 + d + eff_e)

    def apply_pam(
        self,
        affinity_factor: float = 1.0,
        gating_factor: float = 1.0,
        *,
        pam_factor: float | None = None,
    ) -> "KineticThreeOpenParallelModel":
        if pam_factor is not None:
            affinity_factor = pam_factor
        for nm, v in (("affinity_factor", affinity_factor),
                      ("gating_factor", gating_factor)):
            v = float(v)
            if not np.isfinite(v):
                raise ValueError(f"{nm} must be finite, got {v!r}")
            if v < 1.0:
                raise ValueError(f"{nm}={v!r} must be >= 1.0.")
        return replace(
            self,
            koff=self.koff / float(affinity_factor),
            koff2=self.koff2 / float(affinity_factor),
            beta=self.beta * float(gating_factor),
            beta_p=self.beta_p * float(gating_factor),
        )

    def q_matrix(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        mod = self.apply_pam(affinity_factor=pam_factor)
        g = max(float(gaba_um), 0.0)

        kon = mod.kon
        koff = mod.koff
        koff2 = mod.koff2
        beta = mod.beta
        alpha = mod.alpha
        delta = mod.delta
        gamma = mod.gamma
        beta_p = mod.beta_p
        alpha_p = mod.alpha_p
        d = mod.d
        r = mod.r

        q = np.array([
            [-2.0 * kon * g,  2.0 * kon * g,        0.0,                               0.0,             0.0,     0.0,      0.0],
            [koff,            -(koff + kon * g),     kon * g,                           0.0,             0.0,     0.0,      0.0],
            [0.0,             2.0 * koff2,          -(2.0 * koff2 + beta + beta_p + d), beta,            0.0,     beta_p,   d  ],
            [0.0,             0.0,                  alpha,                             -(alpha + delta), delta,   0.0,      0.0],
            [0.0,             0.0,                  0.0,                               gamma,           -gamma,  0.0,      0.0],
            [0.0,             0.0,                  alpha_p,                           0.0,             0.0,     -alpha_p, 0.0],
            [0.0,             0.0,                  r,                                 0.0,             0.0,     0.0,      -r ],
        ], dtype=float)
        return q

    def state_distribution(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        mod = self.apply_pam(affinity_factor=pam_factor)
        g = max(float(gaba_um), 0.0)
        x = g / mod.kd_um
        c = mod.c
        e = mod.gating_efficacy
        k2 = mod.two_open_ratio
        ep = mod.parallel_gating_efficacy
        d = mod.desens_ratio

        eff_o = e * (1.0 + k2) + ep
        denom = 1.0 + 2.0 * x + c * (x ** 2) * (1.0 + d + eff_o)
        if denom <= 0.0:
            return np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        p_r = 1.0 / denom
        p_ar = (2.0 * x) / denom
        p_a2r = (c * (x ** 2)) / denom
        p_a2o = (c * e * (x ** 2)) / denom
        p_a2o_prime = (c * k2 * e * (x ** 2)) / denom
        p_a2o_parallel = (c * ep * (x ** 2)) / denom
        p_a2d = (c * d * (x ** 2)) / denom

        return np.array([
            p_r, p_ar, p_a2r, p_a2o, p_a2o_prime, p_a2o_parallel, p_a2d
        ], dtype=float)

    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float:
        dist = self.state_distribution(gaba_um, pam_factor)
        return float(dist[OPEN_IDX_1] + dist[OPEN_IDX_2] + dist[OPEN_IDX_3])

    def dose_response(
        self,
        concs_um: np.ndarray,
        pam_factor: float = 1.0,
        *,
        observable: Observable = Observable.EQUILIBRIUM,
    ) -> np.ndarray:
        if observable is not Observable.EQUILIBRIUM:
            raise ObservableMismatch(
                f"{self.name} is EQUILIBRIUM-native. Use peak_dose_response for PEAK data."
            )
        concs = np.asarray(concs_um, dtype=float)
        mod = self.apply_pam(affinity_factor=pam_factor)
        x = np.maximum(concs, 0.0) / mod.kd_um
        c = mod.c
        e = mod.gating_efficacy
        k2 = mod.two_open_ratio
        ep = mod.parallel_gating_efficacy
        d = mod.desens_ratio
        eff_o = e * (1.0 + k2) + ep
        denom = 1.0 + 2.0 * x + c * (x ** 2) * (1.0 + d + eff_o)
        return (c * eff_o * (x ** 2)) / denom

    def peak_dose_response(
        self,
        concs_um: np.ndarray,
        pam_factor: float = 1.0,
        application_ms: float = PEAK_APPLICATION_MS,
    ) -> np.ndarray:
        return peak_open_probability_constant(
            self.q_matrix,
            self.state_distribution(0.0),
            concs_um,
            pam_factor=pam_factor,
            application_ms=application_ms,
            open_index=OPEN_INDICES_7,
        )

    def simulate_waveform(
        self,
        t: np.ndarray,
        gaba_t: np.ndarray,
        pam_factor: float = 1.0,
        initial_state: np.ndarray | None = None,
    ) -> WaveformResult:
        t_arr = np.asarray(t, dtype=float)
        gaba_arr = np.asarray(gaba_t, dtype=float)

        if len(t_arr) < 2:
            dist = (self.state_distribution(float(gaba_arr[0]), pam_factor)
                    if len(gaba_arr) > 0 else np.array([1.0, 0, 0, 0, 0, 0, 0]))
            p0 = float(dist[OPEN_IDX_1] + dist[OPEN_IDX_2] + dist[OPEN_IDX_3])
            return WaveformResult(
                t=t_arr, gaba=gaba_arr, p_open=np.full_like(t_arr, p0), peak_p_open=p0,
                charge_integral=0.0, decay_tau_ms=float("nan"), steady_state_tail=p0,
                states=dist[None, :], decay=DecayFit.failed(),
            )

        p_init = resolve_initial_state(
            initial_state, 7, self.state_distribution(float(gaba_arr[0]), pam_factor),
            open_index=OPEN_INDICES_7,
        )

        mod = self.apply_pam(affinity_factor=pam_factor)
        kon = mod.kon
        koff = mod.koff
        koff2 = mod.koff2
        beta = mod.beta
        alpha = mod.alpha
        delta = mod.delta
        gamma = mod.gamma
        beta_p = mod.beta_p
        alpha_p = mod.alpha_p
        d_rate = mod.d
        r_rate = mod.r

        def rhs(time_val: float, y: np.ndarray) -> np.ndarray:
            g = float(np.interp(time_val, t_arr, gaba_arr))
            dp0 = -2.0 * kon * g * y[0] + koff * y[1]
            dp1 = 2.0 * kon * g * y[0] - (koff + kon * g) * y[1] + 2.0 * koff2 * y[2]
            dp2 = kon * g * y[1] - (2.0 * koff2 + beta + beta_p + d_rate) * y[2] + alpha * y[3] + alpha_p * y[5] + r_rate * y[6]
            dp3 = beta * y[2] - (alpha + delta) * y[3] + gamma * y[4]
            dp4 = delta * y[3] - gamma * y[4]
            dp5 = beta_p * y[2] - alpha_p * y[5]
            dp6 = d_rate * y[2] - r_rate * y[6]
            return np.array([dp0, dp1, dp2, dp3, dp4, dp5, dp6])

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
                f"Radau failed to integrate the 7-state parallel scheme over "
                f"[{t_arr[0]:.4g}, {t_arr[-1]:.4g}] ms: {sol.message}."
            )

        states = sol.y.T
        row_sums = states.sum(axis=1)
        worst = float(np.max(np.abs(row_sums - 1.0))) if row_sums.size else 0.0
        if worst > 1e-3:
            raise RuntimeError(
                f"state probabilities drifted off the simplex by {worst:.3e} (max |sum - 1|) "
                f"while integrating the 7-state parallel scheme. Tighten rtol/atol or densify `t`."
            )
        states = np.clip(states, 0.0, 1.0)
        states = states / np.maximum(states.sum(axis=1, keepdims=True), 1e-12)

        p_open = states[:, OPEN_IDX_1] + states[:, OPEN_IDX_2] + states[:, OPEN_IDX_3]
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
