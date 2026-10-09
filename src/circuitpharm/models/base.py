"""Abstract protocol and base result structures for GABA-A receptor models.

CHLORIDE REVERSAL, STATED ONCE. `cpg.E_REV["gabaa"] = -75.0` mV is the circuit layer's
GABA-A reversal, and the receptor layer uses the same number so the two cannot drift (see
`DEFAULT_E_CL_MV`). A GABA-A current is INWARD-NEGATIVE at a holding potential below E_Cl
and the magnitude scales with |V_hold - E_Cl|, so a driving force must be stated; it has
no safe default. See `WaveformResult.current_pA`.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Protocol, runtime_checkable, Any
import numpy as np


#: GABA-A chloride reversal, mV. Kept equal to `cpg.E_REV["gabaa"]` deliberately -- a
#: second copy of a constant is how this project's most expensive error happened (see
#: `gabaa_kinetics.py`'s note on SYNAPTIC_PULSE). Asserted equal in tests.
DEFAULT_E_CL_MV = -75.0

#: Holding potential for a whole-cell GABA-A recording, mV. Typical voltage-clamp holding
#: potential; with `DEFAULT_E_CL_MV` this gives a 15 mV driving force. PROVISIONAL -- it is
#: a convention for turning P_open into a current, not a measurement.
DEFAULT_V_HOLD_MV = -60.0

#: Square-application duration for a PEAK dose-response, ms. The same 300 ms
#: `gabaa_kinetics.Scheme.po_peak` uses; see `ReceptorModel.peak_dose_response` for why it
#: must not be shortened.
PEAK_APPLICATION_MS = 300.0


class Observable(str, Enum):
    """What a number is a number OF.

    These are NOT interchangeable and must never be compared, fitted across, or ratioed
    (roadmap §2.3). The repository's signature failure is exactly this class of mistake:
    a preBotzinger circuit calibrated against whole-body minute ventilation, now VOID.

    PEAK         max P_open during an agonist application -- fast-perfusion
                 concentration-response curves report this.
    EQUILIBRIUM  stationary P_open with desensitisation included -- the tonic /
                 extrasynaptic current reports this. For the 5-state scheme it differs
                 from PEAK by a factor of ~3 at the same parameters.
    CHARGE       integral of P_open over a DECLARED window -- synaptic charge transfer,
                 and the bridge to the circuit layer's conductance.
    """

    PEAK = "PEAK"
    EQUILIBRIUM = "EQUILIBRIUM"
    CHARGE = "CHARGE"


class ObservableMismatch(TypeError):
    """Raised when two quantities on different observables are combined or compared."""


@dataclass(frozen=True)
class ObservedQuantity:
    """A number together with the observable it is a number of.

    Arithmetic between two of these raises unless the observables match. This is the
    enforcement half of roadmap §2.3; the declarations live on the model methods.
    """

    value: float
    observable: Observable
    window_ms: float | None = None

    def __post_init__(self):
        if self.observable is Observable.CHARGE:
            if self.window_ms is None or not (self.window_ms > 0):
                raise ValueError(
                    "a CHARGE quantity needs an explicit positive `window_ms`: the integral "
                    "of P_open depends on the window it was taken over, so a charge without "
                    "one is not a quantity. Pass the integration span in ms.")

    def _check(self, other: "ObservedQuantity", op: str) -> None:
        if not isinstance(other, ObservedQuantity):
            raise ObservableMismatch(
                f"cannot {op} an ObservedQuantity with {type(other).__name__}; wrap the "
                f"other operand and declare its observable.")
        if other.observable is not self.observable:
            raise ObservableMismatch(
                f"cannot {op} {self.observable.value} with {other.observable.value}. These "
                f"are different measurements, not different units: for the 5-state scheme "
                f"PEAK and EQUILIBRIUM EC50 differ by ~3x at identical parameters. If you "
                f"mean to compare them, say which one the claim is about and convert "
                f"explicitly.")
        if (self.observable is Observable.CHARGE
                and other.window_ms is not None
                and self.window_ms is not None
                and abs(self.window_ms - other.window_ms) > 1e-9):
            raise ObservableMismatch(
                f"cannot {op} two CHARGE quantities over different windows "
                f"({self.window_ms} ms vs {other.window_ms} ms).")

    def __truediv__(self, other):
        self._check(other, "divide")
        return self.value / other.value

    def __sub__(self, other):
        self._check(other, "subtract")
        return ObservedQuantity(self.value - other.value, self.observable, self.window_ms)

    def __add__(self, other):
        self._check(other, "add")
        return ObservedQuantity(self.value + other.value, self.observable, self.window_ms)

    def __lt__(self, other):
        self._check(other, "compare")
        return self.value < other.value

    def __gt__(self, other):
        self._check(other, "compare")
        return self.value > other.value


@dataclass(frozen=True)
class DecayFit:
    """Biexponential decay of a post-peak tail.

    EVERY FIELD IS NaN WHEN THE FIT DOES NOT CONVERGE, and `r_squared` is 0.0. It must
    never fall back to a nominal value: the previous implementation returned the literal
    15.0 on failure, which IS `gabaa_kinetics.FIT_TARGETS["tau_ms"]`, so a failed estimate
    was indistinguishable from a perfect reproduction of the anchor (roadmap P0-4).
    """

    tau_fast_ms: float
    tau_slow_ms: float
    weight_fast: float
    tau_weighted_ms: float
    r_squared: float
    n_points: int

    @classmethod
    def failed(cls) -> "DecayFit":
        nan = float("nan")
        return cls(nan, nan, nan, nan, 0.0, 0)

    @property
    def ok(self) -> bool:
        return bool(np.isfinite(self.tau_weighted_ms))


def resolve_initial_state(initial_state, n_states: int, resting: np.ndarray) -> np.ndarray:
    """Coerce a caller's initial condition into a full state vector on the simplex.

    Accepts:
      * None                  -> `resting`
      * a scalar open prob    -> that value in the open state, with the remaining mass
                                 distributed over `resting`'s NON-open states in their
                                 resting proportions
      * a length-n vector     -> validated (non-negative, sums to 1 within 1e-8)

    WHY A SCALAR IS HONOURED RATHER THAN REPLACED. `models/base.py` documents
    `initial_state` as "state vector OR initial open probability", and the previous code
    called `len()` on it -- `TypeError: len() of unsized object` for a 0-d array, and
    `object of type 'float' has no len()` for a float. The near-miss fix is worse than the
    crash: `kinetic_jw95.py` silently substituted the resting distribution, discarding the
    caller's instruction without a word (roadmap P0-1).

    `open_index` is taken as the index of `resting`'s own open state by convention: both
    schemes put A2O at index 3, which is asserted in the model modules.
    """
    resting = np.asarray(resting, dtype=float)
    if resting.shape != (n_states,):
        raise ValueError(
            f"resting distribution has shape {resting.shape}, expected ({n_states},)")

    if initial_state is None:
        return resting.copy()

    arr = np.asarray(initial_state, dtype=float)

    if arr.ndim == 0:
        p_open = float(arr)
        if not (0.0 <= p_open <= 1.0):
            raise ValueError(
                f"a scalar initial_state is an open probability and must lie in [0, 1]; "
                f"got {p_open}. Pass a length-{n_states} state vector to set the full "
                f"distribution.")
        out = np.zeros(n_states, dtype=float)
        out[OPEN_STATE_INDEX] = p_open
        others = [i for i in range(n_states) if i != OPEN_STATE_INDEX]
        mass = float(resting[others].sum())
        if mass > 1e-12:
            out[others] = resting[others] * ((1.0 - p_open) / mass)
        else:
            # resting is entirely in the open state (degenerate); spread the remainder
            # uniformly rather than silently returning a vector that does not sum to 1
            out[others] = (1.0 - p_open) / len(others)
        return out

    if arr.ndim != 1 or arr.shape[0] != n_states:
        raise ValueError(
            f"initial_state must be None, a scalar open probability, or a length-"
            f"{n_states} state vector; got shape {arr.shape}. This model has {n_states} "
            f"states.")
    if np.any(arr < -1e-12):
        raise ValueError(f"initial_state has negative entries: {arr}")
    total = float(arr.sum())
    if abs(total - 1.0) > 1e-8:
        raise ValueError(
            f"initial_state sums to {total!r}, not 1.0 (tolerance 1e-8). A state vector is "
            f"a probability distribution; rescale it or pass a scalar open probability "
            f"instead.")
    return np.clip(arr, 0.0, 1.0)


#: Index of the open state in both the 5- and 6-state schemes (A2O). Asserted in each
#: model module against its own STATES tuple, so a reordering cannot desynchronise them.
OPEN_STATE_INDEX = 3


def fit_biexponential_decay(t_ms, y, *, hi_frac: float = 1.00,
                            lo_frac: float = 0.01) -> DecayFit:
    """Fit A*exp(-t/tau_f) + B*exp(-t/tau_s) to a decaying tail.

    WINDOW: the whole post-peak decay by default, down to 1% of peak.

    This deliberately DIFFERS from `gabaa_kinetics.Scheme.ipsc_metrics`, which fits a
    SINGLE exponential over the 90%->10% window, and the difference is not a second
    convention for one quantity -- it is the right convention for a different quantity.
    A 90->10% window is the standard experimental window for a mono-exponential
    deactivation tau, and it necessarily truncates the slow component of a two-component
    decay, because the slow component only dominates below ~10% of peak. Measured on
    0.7*exp(-t/15) + 0.3*exp(-t/70) over 400 ms at 0.1 ms resolution, whose true weighted
    tau is 31.5 ms: the 90->10% window recovers 32.605 ms (+3.5%) while the full tail
    recovers 31.241 ms (-0.8%).

    Both numbers are defensible answers to different questions, so the window is a
    parameter and the default is the one that recovers the parameters. Pass
    `hi_frac=0.90, lo_frac=0.10` for the experimental deactivation window.

    Returns `DecayFit.failed()` -- every field NaN, `r_squared` 0.0 -- when the tail is
    flat, too short, or the optimiser does not converge. It NEVER returns a nominal value:
    the previous single-point estimator fell back to the literal 15.0, which is the project's
    IPSC tau fit target, so a failure read as a perfect fit (roadmap P0-4).

    `tau_weighted_ms` is the amplitude-weighted mean, (A*tf + B*ts)/(A + B), which is the
    quantity experimental papers report as "weighted tau".
    """
    from scipy.optimize import least_squares

    t_ms = np.asarray(t_ms, dtype=float)
    y = np.asarray(y, dtype=float)
    if t_ms.shape != y.shape or t_ms.ndim != 1 or t_ms.size < 8:
        return DecayFit.failed()
    if not (np.all(np.isfinite(t_ms)) and np.all(np.isfinite(y))):
        return DecayFit.failed()

    # work from the peak onward, with the tail's own asymptote removed
    i_pk = int(np.argmax(y))
    tt = t_ms[i_pk:] - t_ms[i_pk]
    yy = y[i_pk:] - float(y[-1])
    pk = float(yy[0]) if yy.size else 0.0
    if pk <= 1e-12 or yy.size < 8:
        return DecayFit.failed()

    mask = (yy <= hi_frac * pk) & (yy >= lo_frac * pk)
    if int(mask.sum()) < 6:
        return DecayFit.failed()
    tw, yw = tt[mask], yy[mask]
    tw = tw - tw[0]
    span = float(tw[-1]) if tw[-1] > 0 else 1.0

    def resid(x):
        a, log_tf, b, log_ts = x
        return (a * np.exp(-tw / 10.0 ** log_tf)
                + b * np.exp(-tw / 10.0 ** log_ts)) - yw

    y0 = float(yw[0])
    x0 = [0.7 * y0, np.log10(max(span / 4.0, 1e-3)),
          0.3 * y0, np.log10(max(span * 2.0, 1e-2))]
    try:
        out = least_squares(
            resid, x0,
            bounds=([0.0, -3.0, 0.0, -3.0], [10.0 * abs(y0) + 1e-9, 5.0,
                                             10.0 * abs(y0) + 1e-9, 5.0]),
            xtol=1e-12, ftol=1e-12, max_nfev=20000)
    except Exception:
        return DecayFit.failed()
    if not out.success:
        return DecayFit.failed()

    a, log_tf, b, log_ts = out.x
    tf, ts = 10.0 ** log_tf, 10.0 ** log_ts
    # order the components so "fast" is genuinely the faster one
    if tf > ts:
        a, tf, b, ts = b, ts, a, tf
    amp = a + b
    if amp <= 1e-15:
        return DecayFit.failed()

    ss_res = float(np.sum(out.fun ** 2))
    ss_tot = float(np.sum((yw - yw.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else 0.0
    if not np.isfinite(r2):
        return DecayFit.failed()

    return DecayFit(
        tau_fast_ms=float(tf), tau_slow_ms=float(ts),
        weight_fast=float(a / amp),
        tau_weighted_ms=float((a * tf + b * ts) / amp),
        r_squared=float(r2), n_points=int(mask.sum()))


@dataclass(frozen=True)
class WaveformResult:
    """Standardized electrophysiological output from a receptor model simulation."""
    t: np.ndarray             # time array in ms
    gaba: np.ndarray          # agonist concentration array in uM
    p_open: np.ndarray        # open probability over time
    peak_p_open: float        # maximum open probability attained
    charge_integral: float    # integral of open probability over time (ms)
    decay_tau_ms: float       # weighted decay constant after the peak (ms); NaN if unfit
    steady_state_tail: float  # final open probability in the trace
    states: np.ndarray | None = None  # optional full state trajectory (n_times x n_states)
    decay: DecayFit | None = None     # the full biexponential fit behind decay_tau_ms

    def current_pA(self, g_max_ns: float, *, v_hold_mv: float = DEFAULT_V_HOLD_MV,
                   e_cl_mv: float = DEFAULT_E_CL_MV) -> np.ndarray:
        """Macroscopic current trace in pA.

        `g_max_ns` IS REQUIRED and the driving force is keyword-only, because this method
        previously defaulted to `v_hold_mv = e_rev_mv = -70.0` and therefore returned an
        identically zero trace for every input -- a silent zero, which
        `protocols.waveforms.extract_electrophys_metrics` then worked around with its own
        separate `driving_force_mv` argument, leaving two conventions for one quantity
        (roadmap P0-11).

        Sign convention: I = g * (V_hold - E_Cl). At V_hold = -60 mV and E_Cl = -75 mV the
        driving force is +15 mV, so an open GABA-A channel passes a positive (outward,
        chloride-inward) current here. Flip `v_hold_mv` below `e_cl_mv` for the inward case.
        """
        if not np.isfinite(g_max_ns):
            raise ValueError(f"g_max_ns must be finite, got {g_max_ns!r}")
        driving = v_hold_mv - e_cl_mv
        if abs(driving) < 1e-12:
            raise ValueError(
                f"v_hold_mv ({v_hold_mv}) equals e_cl_mv ({e_cl_mv}), so the driving force "
                f"is zero and every current would be 0 pA regardless of P_open. That was "
                f"this method's old default and it silently produced empty traces. State a "
                f"holding potential away from the reversal.")
        # nS * mV = pA
        return g_max_ns * self.p_open * driving


@runtime_checkable
class ReceptorModel(Protocol):
    """Protocol satisfied by all competing GABA-A kinetic and operational models."""

    @property
    def native_observable(self) -> Observable:
        """Which observable this model's `dose_response` speaks natively.

        Model A (operational Hill, no desensitisation) is PEAK-native: its curve IS the
        peak response and it has no separate equilibrium. Models B and C are
        EQUILIBRIUM-native: `dose_response` is the stationary distribution with
        desensitisation included, and the PEAK curve has to be integrated for.

        That asymmetry is a real model difference, not an inconvenience, and P5's
        comparison must see it rather than average over it.
        """
        ...

    @property
    def name(self) -> str:
        """Human-readable identifier of the model."""
        ...

    @property
    def param_names(self) -> tuple[str, ...]:
        """Names of the model's tunable parameters."""
        ...

    def get_params(self) -> dict[str, float]:
        """Return a dictionary of current parameters."""
        ...

    def with_params(self, **kwargs: float) -> "ReceptorModel":
        """Return a copy of the model with updated parameter values."""
        ...

    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float:
        """Calculate equilibrium open probability at constant GABA concentration.

        Parameters
        ----------
        gaba_um : float
            Agonist concentration in micromolar.
        pam_factor : float
            Allosteric modulation multiplier (>= 1.0 for positive modulation).
        """
        ...

    def dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0) -> np.ndarray:
        """Open probability across concentrations, on this model's `native_observable`."""
        ...

    def peak_dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0,
                           application_ms: float = PEAK_APPLICATION_MS) -> np.ndarray:
        """PEAK open probability during a square agonist application.

        `application_ms` defaults to the same 300 ms that
        `gabaa_kinetics.Scheme.po_peak` uses, and for the reason recorded there: a 1 ms
        step is binding-rate-limited at low agonist, which inflates the fitted EC50,
        forces an implausibly high microscopic affinity, and compresses every derived PAM
        gain. Published concentration-response curves come from applications of hundreds
        of ms. Do not pick a new number here.
        """
        ...

    def simulate_waveform(
        self,
        t: np.ndarray,
        gaba_t: np.ndarray,
        pam_factor: float = 1.0,
        initial_state: np.ndarray | None = None,
    ) -> WaveformResult:
        """Simulate time-domain open probability in response to a time-varying GABA transient.

        Parameters
        ----------
        t : np.ndarray
            1D array of timepoints in ms.
        gaba_t : np.ndarray
            1D array of GABA concentrations in uM matching `t`.
        pam_factor : float
            Allosteric modulation factor.
        initial_state : np.ndarray, optional
            Initial state vector or initial open probability.
        """
        ...
