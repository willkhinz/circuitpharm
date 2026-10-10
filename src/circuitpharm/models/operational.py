"""Model A: Phenomenological Operational Scalar Ceiling Model.

Represents the hypothesis that positive allosteric modulation acts via a simple
operational potency/efficacy shift subject to an extrinsic scalar ceiling, without
explicit microscopic Markov state redistribution.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp

from .base import (PEAK_APPLICATION_MS, DecayFit, Observable, ObservableMismatch,
                   ReceptorModel, WaveformResult, fit_biexponential_decay)


@dataclass(frozen=True)
class OperationalScalarModel(ReceptorModel):
    """Operational / Hill model with an imposed scalar ceiling on allosteric gain.

    DEFAULTS ARE PROVISIONAL. `ec50_um = 25.0` is a declared round number, not a fit, and
    `po_max`/`tau_deact_ms` happen to coincide with the project's fit TARGETS rather than
    having been fitted to anything. Obtain parameters from
    `circuitpharm.parameters.get(...)` (roadmap P0-13).
    """
    ec50_um: float = 25.0       # basal GABA EC50 (uM)   PROVISIONAL
    hill_n: float = 1.4         # Hill coefficient       PROVISIONAL
    po_max: float = 0.75        # max open probability   PROVISIONAL (= the fit target)
    s_max: float = 2.50         # ceiling on potency shift  PROVISIONAL
    tau_act_ms: float = 0.5     # activation time constant (ms)    PROVISIONAL
    tau_deact_ms: float = 15.0  # deactivation time constant (ms)  PROVISIONAL (= target)

    @property
    def name(self) -> str:
        return "OperationalScalar"

    @property
    def native_observable(self) -> Observable:
        """PEAK-native: this model has no desensitisation, so no separate equilibrium.

        That is a real structural difference from Models B and C, and P5's comparison must
        see it rather than average over it. Asking this model for an EQUILIBRIUM curve
        raises (see `dose_response`).
        """
        return Observable.PEAK

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
        """Apply operational potency shift with ceiling s_max.

        REJECTS pam_factor < 1.0, as models B and C already do. This model clamped with
        `max(shift, 1e-6)` instead, so `effective_ec50(0.5)` returned 2x the EC50 -- a
        silent NEGATIVE allosteric modulator -- while `KineticAllosteryModel.apply_pam` and
        `ExtendedDesensitizationModel.apply_pam` both raise ValueError on the same input.
        Verified before the fix: 25 -> 50 uM at pam_factor=0.5.

        An asymmetric interface across the three models is worse than any one of them being
        wrong, because model COMPARISON is what this layer exists for: a sweep that silently
        measures a NAM in model A and a rejection in B and C compares two different
        questions.
        """
        if float(pam_factor) < 1.0:
            raise ValueError(
                f"pam_factor={pam_factor!r} is below 1.0, which is negative allosteric "
                f"modulation. This operational model represents POSITIVE modulation only "
                f"(roadmap P0-12), as do the kinetic models -- a NAM needs its own "
                f"parameterisation, not a reciprocal potency shift.")
        shift = min(float(pam_factor), self.s_max)
        return self.ec50_um / max(shift, 1e-6)

    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float:
        """Stationary Hill open probability."""
        if gaba_um <= 0.0:
            return 0.0
        eff_ec50 = self.effective_ec50(pam_factor)
        conc_ratio = (gaba_um / eff_ec50) ** self.hill_n
        return float(self.po_max * (conc_ratio / (1.0 + conc_ratio)))

    def dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0, *,
                      observable: Observable = Observable.PEAK) -> np.ndarray:
        """The Hill curve, which for this model is the PEAK response.

        An EQUILIBRIUM request RAISES rather than returning this curve relabelled. The
        model has no desensitisation, so it has no stationary state distinct from its
        peak; handing the same numbers back under a different name is exactly the
        cross-observable comparison §2.3 forbids, and it is what let a PEAK-native model
        be ranked against two EQUILIBRIUM-native ones on one axis in `oed.py`.
        """
        if observable is not Observable.PEAK:
            raise ObservableMismatch(
                f"{self.name} is PEAK-native and has no {observable.value} curve: it "
                f"contains no desensitisation, so its stationary response IS its peak. "
                f"Returning the Hill curve under an EQUILIBRIUM label would assert a "
                f"structural property the model does not have. If you are comparing "
                f"models, compare them on PEAK, which all three can produce.")
        concs = np.asarray(concs_um, dtype=float)
        eff_ec50 = self.effective_ec50(pam_factor)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = np.where(concs > 0, (concs / eff_ec50) ** self.hill_n, 0.0)
            res = self.po_max * (ratios / (1.0 + ratios))
        return np.nan_to_num(res, nan=0.0)

    def peak_dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0,
                           application_ms: float = PEAK_APPLICATION_MS) -> np.ndarray:
        """Identical to `dose_response` here, and that identity is the model's content.

        For Models B and C this costs an integration and lands ~3x away from the
        equilibrium curve. For this model the two coincide by construction, because there
        is no desensitisation to separate them. `application_ms` is accepted and ignored,
        which is itself a prediction: this model says the CRC does not depend on
        application length, and Models B and C say it does.
        """
        return self.dose_response(concs_um, pam_factor, observable=Observable.PEAK)

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
            p0 = (self.steady_state(float(gaba_arr[0]), pam_factor)
                  if len(gaba_arr) > 0 else 0.0)
            return WaveformResult(
                t=t_arr, gaba=gaba_arr, p_open=np.full_like(t_arr, p0), peak_p_open=p0,
                charge_integral=0.0, decay_tau_ms=float("nan"), steady_state_tail=p0,
                decay=DecayFit.failed())

        # P0-1 (mirror form). `float(initial_state[0])` read element 0 of a state vector as
        # an open probability -- for the Markov models' resting distribution that is P_R,
        # so handing this model a resting state started it at p_open = 1.0. A scalar is
        # the only meaningful initial condition for a one-variable relaxation, so a vector
        # is rejected rather than silently indexed.
        if initial_state is None:
            p_init = self.steady_state(float(gaba_arr[0]), pam_factor)
        else:
            arr = np.asarray(initial_state, dtype=float)
            if arr.ndim != 0 and arr.size != 1:
                raise ValueError(
                    f"{self.name} has ONE state variable (open probability), so "
                    f"initial_state must be None or a scalar; got shape {arr.shape}. "
                    f"Passing a multi-state vector here previously read its first element "
                    f"as the open probability, which for a Markov resting distribution is "
                    f"P_R -- i.e. it started the trace fully open.")
            p_init = float(arr.reshape(()))
            if not (0.0 <= p_init <= 1.0):
                raise ValueError(
                    f"initial_state is an open probability and must lie in [0, 1]; got "
                    f"{p_init}")

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

        # P0-3, and this one was the worst of the three fallbacks: substituting `target_p`
        # replaces the relaxation with the instantaneous-Hill limit, i.e. it DELETES the
        # dynamics this model exists to express and biases peak and charge upward. A
        # failure here raises.
        if not sol.success:
            raise RuntimeError(
                f"RK23 failed to integrate the operational relaxation over "
                f"[{t_arr[0]:.4g}, {t_arr[-1]:.4g}] ms: {sol.message}. max_step was "
                f"{max_step:.4g} ms against tau_act {self.tau_act_ms:.4g} ms. Substituting "
                f"the instantaneous Hill target would remove the relaxation entirely and "
                f"overstate both peak and charge.")

        p_open = np.clip(sol.y[0], 0.0, 1.0)
        peak_p = float(np.max(p_open))

        # Integral via trapezoid
        _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
        charge = float(_trapz(p_open, t_arr))

        # P0-4: real biexponential fit, NaN on failure. The old fallback to
        # `self.tau_deact_ms` was the same defect as Models B and C falling back to 15.0,
        # only self-referential: the model returned its own input parameter as a
        # "measurement" of itself, so a failed fit looked like perfect agreement.
        decay = fit_biexponential_decay(t_arr, p_open)

        return WaveformResult(
            t=t_arr,
            gaba=gaba_arr,
            p_open=p_open,
            peak_p_open=peak_p,
            charge_integral=charge,
            decay_tau_ms=decay.tau_weighted_ms,
            steady_state_tail=float(p_open[-1]),
            decay=decay,
        )
