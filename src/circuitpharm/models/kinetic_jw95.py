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

from .base import (OPEN_STATE_INDEX, PEAK_APPLICATION_MS, DecayFit, Observable,
                   ReceptorModel, WaveformResult, fit_biexponential_decay,
                   resolve_initial_state)


STATES_5 = ("R", "AR", "A2R", "A2O", "A2D")
OPEN_IDX_5 = 3  # "A2O"
assert STATES_5[OPEN_STATE_INDEX] == "A2O", "open-state index desynchronised from base"


@dataclass(frozen=True)
class KineticAllosteryModel(ReceptorModel):
    """5-state Markov gating scheme with microscopic allosteric modulation.

    EVERY DEFAULT BELOW IS PROVISIONAL AND IS NOT A CALIBRATION. They reproduce neither
    of the project's two anchors -- measured EC50 6.34 uM (EQUILIBRIUM) against the 20 uM
    anchor, and P_o,max 0.8382 against 0.750 -- and they are a CHIMERA of two different
    fits: `beta`/`alpha` are the manuscript's current config-pulse fit, while `kon`/`koff`
    are approximately the superseded 1000 uM / 0.30 ms fit. The resulting
    K_d = 29.62 uM appears in no fit, no commit and no document.

    Obtain parameters from `circuitpharm.parameters.get(...)` instead. These exist so the
    algebra can be unit-tested, and `tests/test_nextgen_models.py` pins the mismatch above
    so they cannot be mistaken for a calibration (roadmap P0-13).
    """
    kon: float = 0.011244       # uM^-1 ms^-1  PROVISIONAL - see class docstring
    koff: float = 0.333069      # ms^-1        PROVISIONAL
    beta: float = 0.559023      # ms^-1        PROVISIONAL
    alpha: float = 0.107891     # ms^-1        PROVISIONAL
    d: float = 0.050            # ms^-1        PROVISIONAL (not fitted by any anchor)
    r: float = 0.0020           # ms^-1        PROVISIONAL (not fitted by any anchor)

    @property
    def name(self) -> str:
        return "KineticAllostery_JW95"

    @property
    def native_observable(self) -> Observable:
        """`dose_response` is the stationary distribution, desensitisation included."""
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

    def apply_pam(self, affinity_factor: float = 1.0,
                  gating_factor: float = 1.0) -> "KineticAllosteryModel":
        """Apply allosteric modulation.

        affinity_factor: fold-decrease in koff (slowing unbinding). At equilibrium this
            divides K_d exactly, so it IS the EQUILIBRIUM EC50 fold-shift -- a factor of
            2.5 here gives a 2.5x left-shift of `dose_response`. Note this is NOT the same
            modulator strength as `gabaa_kinetics.calibrate_pam`'s output, which calibrates
            against the PEAK EC50 shift and returns 2.945 for the same nominal 2.5x
            (roadmap P6-2). Both are internally consistent; the two must not be swapped.
        gating_factor: fold-increase in beta (accelerating opening). Raises P_o,max, so it
            is not bounded the way an affinity shift is.

        POSITIVE MODULATION ONLY. Values below 1.0 are rejected rather than clamped: the
        surrounding code (headroom, decomposition, dose-escalation) is written for
        potentiation, `calibrate_pam` documents NAMs as out of domain, and the old
        `max(factor, 1e-6)` admitted 0.0, which turned k_off into 1e6 x its value --
        a silent super-agonist. One decision, stated here, for all three models
        (roadmap P0-12).
        """
        for nm, v in (("affinity_factor", affinity_factor),
                      ("gating_factor", gating_factor)):
            v = float(v)
            if not np.isfinite(v):
                raise ValueError(f"{nm} must be finite, got {v!r}")
            if v < 1.0:
                raise ValueError(
                    f"{nm}={v!r} is below 1.0, i.e. NEGATIVE allosteric modulation, which "
                    f"this scheme's callers are not written for. `calibrate_pam` searches "
                    f"multipliers >= 1 and returns NaN below that by design; extending to "
                    f"NAMs is a feature, not a clamp. Pass >= 1.0.")
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
        """EQUILIBRIUM open probability. See `native_observable`.

        This is NOT the peak-current concentration-response that published CRCs report --
        it differs by a factor of ~3 in EC50 at identical parameters, because the
        stationary distribution has desensitisation competing throughout. Fitting a
        peak-current dataset with this curve is the observable mismatch that produced a
        chi-squared of 9290 for 9 data points (roadmap C2). Use `peak_dose_response` for
        PEAK data.
        """
        concs = np.asarray(concs_um, dtype=float)
        mod = self.apply_pam(affinity_factor=pam_factor)
        x = np.maximum(concs, 0.0) / mod.kd_um
        e = mod.gating_efficacy
        d = mod.desens_ratio
        denom = 1.0 + 2.0 * x + (x ** 2) * (1.0 + e + d)
        return (e * (x ** 2)) / denom

    def peak_dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0,
                           application_ms: float = PEAK_APPLICATION_MS) -> np.ndarray:
        """PEAK open probability during a square application. See the protocol docstring."""
        concs = np.atleast_1d(np.asarray(concs_um, dtype=float))
        out = np.empty(concs.shape, dtype=float)
        for i, c in enumerate(concs.ravel()):
            t = np.linspace(0.0, float(application_ms), 600)
            g = np.full_like(t, max(float(c), 0.0))
            res = self.simulate_waveform(t, g, pam_factor=pam_factor,
                                         initial_state=self.state_distribution(0.0))
            out.ravel()[i] = res.peak_p_open
        return out

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
            dist = (self.state_distribution(float(gaba_arr[0]), pam_factor)
                    if len(gaba_arr) > 0 else np.array([1.0, 0, 0, 0, 0]))
            p0 = float(dist[OPEN_IDX_5])
            # decay is UNDEFINED on a single sample, so it is NaN rather than 0.0 or a
            # nominal tau (P0-4).
            return WaveformResult(
                t=t_arr, gaba=gaba_arr, p_open=np.full_like(t_arr, p0), peak_p_open=p0,
                charge_integral=0.0, decay_tau_ms=float("nan"), steady_state_tail=p0,
                states=dist[None, :], decay=DecayFit.failed())

        # P0-1: a scalar open probability is honoured, a vector is validated, anything
        # else raises with what it accepts. `len()` on a 0-d array used to crash here.
        p_init = resolve_initial_state(
            initial_state, 5, self.state_distribution(float(gaba_arr[0]), pam_factor))

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

        # P0-3: a failed integration RAISES. It used to be replaced by a constant trace
        # tiled from p_init, which yields a flat p_open, a near-zero charge and (with the
        # old estimator) a decay tau of exactly the fit target -- three plausible numbers
        # from an experiment that never ran.
        if not sol.success:
            raise RuntimeError(
                f"Radau failed to integrate the 5-state scheme over "
                f"[{t_arr[0]:.4g}, {t_arr[-1]:.4g}] ms: {sol.message}. The stiffest "
                f"timescale here is 1/r = {1.0 / self.r:.0f} ms against a cleft transient "
                f"of ~1 ms, so this is a genuinely stiff problem; if the parameters are "
                f"legitimate, increase the density of `t` rather than loosening rtol.")

        states = sol.y.T
        # Conservation is a property of the generator (rows sum to zero), so a drift here
        # is integration error, not something to normalise away quietly. Check it, then
        # clip only the tiny negatives a stiff solver legitimately produces.
        row_sums = states.sum(axis=1)
        worst = float(np.max(np.abs(row_sums - 1.0))) if row_sums.size else 0.0
        if worst > 1e-3:
            raise RuntimeError(
                f"state probabilities drifted off the simplex by {worst:.3e} (max |sum - 1|) "
                f"while integrating the 5-state scheme. Renormalising would hide an "
                f"integration failure behind a plausible trace; tighten rtol/atol or "
                f"densify `t` instead.")
        states = np.clip(states, 0.0, 1.0)
        states = states / np.maximum(states.sum(axis=1, keepdims=True), 1e-12)

        p_open = states[:, OPEN_IDX_5]
        peak_p = float(np.max(p_open))

        _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
        charge = float(_trapz(p_open, t_arr))

        # P0-4: a real biexponential fit, NaN when it does not converge. No fallback to
        # 15.0, which is gabaa_kinetics.FIT_TARGETS["tau_ms"].
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
