"""GABA-A receptor gating as an explicit Markov scheme — replacing three hand-set numbers
with one mechanism.

WHY THIS EXISTS. The model currently carries three INDEPENDENT free parameters for a
positive allosteric modulator:

    Drug.gaba_a_gain          multiplier on peak conductance     (hand-set, e.g. 2.0)
    Drug.gaba_a_tau           multiplier on decay tau            (hand-set as 1 + 0.6*(g-1))
    Drug.gaba_a_efficacy_cap  the ceiling                        (hand-set, 2.5 or 6.0)

In reality those are not three things. They are three CONSEQUENCES of one thing: how the
modulator shifts the rate constants of receptor gating. Tying them to a single kinetic
scheme removes two free parameters, lets the ceiling EMERGE instead of being asserted
(and the ceiling is what the uncertainty analysis identified as the binding constraint),
and makes the whole GABA-A arm fittable to patch-clamp data that already exists.

THE MECHANISTIC PAYOFF — the synaptic/extrasynaptic asymmetry comes out for free.
Classical benzodiazepine-site PAMs act predominantly by increasing APPARENT AGONIST
AFFINITY: they shift the GABA dose-response curve leftward with little change in the
maximum current. That single fact, pushed through the scheme, predicts all of:

  * SYNAPTIC receptors see ~1 mM GABA for ~1 ms -- essentially saturating. A leftward
    shift cannot raise an already-maximal peak, so peak amplitude barely changes. But
    slower agonist unbinding prolongs the DECAY. -> tau effect, little gain effect.
  * EXTRASYNAPTIC/TONIC receptors see ~0.1-1 uM ambient GABA -- far below EC50, on the
    steep part of the curve. A leftward shift there produces a LARGE fractional increase.
    -> big gain effect.
  * The CEILING is then not a free parameter at all -- but it is also NOT ONE NUMBER, and
    getting this wrong is easy. A pure affinity modulator can do no more than move the
    receptor up its own dose-response curve, so its maximum effect is set by HOW FAR BELOW
    SATURATION the agonist already sits. That differs enormously between the two pools:

      SYNAPTIC  pool sees a near-saturating transient -> headroom ~1.1x. The drug can do
                almost nothing. The "needs endogenous GABA / saturates" property holds.
      TONIC     pool sees ~0.4 uM against an EC50 of ~20 uM -> headroom ~210x (measured
                below; the asymptote is set by steady-state desensitisation, not by
                Po_max). The drug has a LOT of room.

    So the project's safety argument -- "the GABA arm's protection comes entirely from the
    efficacy ceiling, you cannot dose-shift out of it" -- is strong for synaptic receptors
    and WEAK for extrasynaptic ones. alpha5 is predominantly extrasynaptic. That is the
    uncomfortable consequence, and it is the opposite of what an earlier draft of this
    docstring claimed (it reported the saturating-agonist gain, ~1.0, as if it were the
    ceiling; that number answers a different question -- what the drug does to an ALREADY
    saturated receptor, not how far it can push an unsaturated one).

    The usable ceiling is therefore a LIGAND property (its maximum achievable EC50 shift,
    i.e. intrinsic allosteric efficacy, ~2-3x for classical BZ-site ligands -> tonic gain
    ~8) convolved with a SYNAPSE property (distance from saturation). The model's single
    `gaba_a_efficacy_cap = 2.5` conflates the two and matches neither pool.
  * A modulator acting on GATING instead (raising beta/alpha, as neurosteroids do) RAISES
    the maximum current and therefore has a much weaker ceiling. The project already
    encodes this difference by hand (nonselective_bz ceiling 2.5 vs neurosteroid 6.0).
    Here it is a prediction, which makes it a VALIDATION TARGET rather than an assumption.

SCHEME. Five states, two sequential agonist binding steps, one open state, one
desensitised state reached from the doubly-bound shut state:

    R  <--2*kon*A / koff-->  AR  <--kon*A / 2*koff-->  A2R  <--beta/alpha-->  A2O (open)
                                                        |
                                                   d /  r
                                                        v
                                                       A2D (desensitised)

The 2x / 2x factors on the binding steps are the statistical factors for two equivalent
sites. Desensitisation is placed on the doubly-bound SHUT state (Jones & Westbrook
topology) rather than the open state.

HONESTY ABOUT THE RATE CONSTANTS. The values in DEFAULT_RATES are provisional, of the
order used in published synaptic GABA-A schemes, and are then FITTED HERE to macroscopic
observables that are robust and well established rather than trusted from memory:

    GABA EC50 (peak current, alpha1beta2gamma2)   ~10-30 uM
    maximal open probability at saturating GABA   ~0.7-0.8
    synaptic IPSC deactivation tau                ~10-30 ms

Any quantitative claim from this module is only as good as that fit, and the fit targets
are ranges, not points. What the module is for is the STRUCTURE of the dependencies -- the
sign and rough magnitude of how gain, tau and ceiling move together -- which is far more
robust than any individual rate constant.

    python -m circuitpharm.gabaa_kinetics            # fit, validate, and derive the three numbers
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp

# np.trapezoid was added in NumPy 2.0 and np.trapz was REMOVED in it, while pyproject
# allows numpy>=1.26 -- so neither name alone works across the supported range. Resolved
# once at import rather than per call.
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
if _trapz is None:                                    # pragma: no cover
    raise ImportError("numpy exposes neither trapezoid nor trapz")

STATES = ("R", "AR", "A2R", "A2O", "A2D")
OPEN = STATES.index("A2O")

# provisional; fitted in fit_scheme(). Units: ms^-1, and uM^-1 ms^-1 for kon.
DEFAULT_RATES = dict(kon=0.010, koff=1.00, beta=1.00, alpha=0.30,
                     d=0.050, r=0.0020)

# Agonist concentrations the two receptor pools actually experience.
#
# ONE DEFINITION, IMPORTED. These were literals here (1000 uM / 0.30 ms) while
# `config.SYNAPTIC_PULSE` held 3000 uM / 1.00 ms and `cpg.Drug.from_kinetics` carried a THIRD
# copy as an inline literal matching config's. Every pharmacology consumer passed config's
# pulse explicitly, so the module defaults were reached only by callers that did not pass one
# -- and those callers were therefore fitting a DIFFERENT scheme from the one the whole
# pipeline uses, with no indication that they were:
#
#     module defaults (1000 / 0.30)   kon 0.0112  koff 0.3331  tonic headroom 184.6x
#     config pulse    (3000 / 1.00)   kon 0.0147  koff 0.4692  tonic headroom 210.7x
#
# A manuscript draft was accused of quoting unreproducible numbers on the strength of the
# first column. It was quoting the second -- the pipeline's own calibration -- and the audit
# that "refuted" it had grepped for `SYNAPTIC_PEAK_UM = 3000`, a symbol that never held the
# value. Three copies of one constant, and the disagreement surfaced only when a script
# compared a Drug's tonic gain against the manuscript's.
#
# config.SYNAPTIC_PULSE is authoritative because it is CALIBRATED: it is one of the 9 of 27
# (peak, clearance) cells in which the scheme reproduced every anchored benzodiazepine
# observable from a single free parameter. The former literals here were a provisional
# literature estimate. `config` imports only `results`, so this direction is safe.
from .config import AMBIENT_GABA_UM as AMBIENT_UM          # noqa: E402
from .config import SYNAPTIC_PULSE as _PULSE               # noqa: E402

SYNAPTIC_PEAK_UM = float(_PULSE["peak_um"])
SYNAPTIC_CLEAR_MS = float(_PULSE["clear_ms"])


@dataclass(frozen=True)
class Scheme:
    kon: float = DEFAULT_RATES["kon"]
    koff: float = DEFAULT_RATES["koff"]
    beta: float = DEFAULT_RATES["beta"]
    alpha: float = DEFAULT_RATES["alpha"]
    d: float = DEFAULT_RATES["d"]
    r: float = DEFAULT_RATES["r"]

    # ---- allosteric modulation -------------------------------------------------------
    def pam(self, affinity: float = 1.0, gating: float = 1.0) -> "Scheme":
        """Apply a positive allosteric modulator.

        affinity > 1 : increases apparent agonist affinity by SLOWING UNBINDING (koff/=f).
                       This is the dominant action of classical BZ-site ligands: it shifts
                       EC50 leftward and, because it cannot exceed the receptor's own
                       maximum, it SATURATES -- the origin of the efficacy ceiling.
        gating   > 1 : increases the opening/closing ratio (beta*=f). This RAISES maximal
                       open probability, so it is not bounded the same way -- the
                       neurosteroid-like action, with a correspondingly weaker ceiling.
        """
        return replace(self, koff=self.koff / affinity, beta=self.beta * gating)

    # ---- machinery -------------------------------------------------------------------
    def Q(self, conc_um: float) -> np.ndarray:
        """Transition-rate matrix. Q[i,j] = rate i->j for i!=j; rows sum to zero."""
        a = conc_um
        Q = np.zeros((5, 5))
        Q[0, 1] = 2.0 * self.kon * a          # R   -> AR
        Q[1, 0] = self.koff                   # AR  -> R
        Q[1, 2] = self.kon * a                # AR  -> A2R
        Q[2, 1] = 2.0 * self.koff             # A2R -> AR
        Q[2, 3] = self.beta                   # A2R -> A2O
        Q[3, 2] = self.alpha                  # A2O -> A2R
        Q[2, 4] = self.d                      # A2R -> A2D
        Q[4, 2] = self.r                      # A2D -> A2R
        np.fill_diagonal(Q, -Q.sum(axis=1))
        return Q

    def steady_state(self, conc_um: float) -> np.ndarray:
        """Equilibrium occupancy. Solved as the stationary distribution of Q."""
        Q = self.Q(conc_um)
        A = np.vstack([Q.T[:-1], np.ones(5)])
        b = np.zeros(5); b[-1] = 1.0
        p, *_ = np.linalg.lstsq(A, b, rcond=None)
        return np.clip(p, 0.0, 1.0)

    def po_tonic(self, conc_um: float = AMBIENT_UM) -> float:
        """Steady-state open probability — what a TONIC/extrasynaptic current reports."""
        return float(self.steady_state(conc_um)[OPEN])

    def po_peak(self, conc_um: float, dur_ms: float = 300.0) -> float:
        """Peak open probability during a sustained step of agonist.

        PROTOCOL MATTERS, and getting it wrong cost a debug cycle here. The first version
        used a 1 ms step. At low agonist concentration binding itself takes longer than
        that (kon*[A] sets the rate), so the measured peak was BINDING-RATE-LIMITED rather
        than equilibrium-limited. That inflated the fitted EC50, which forced the fit to an
        implausibly high microscopic affinity, which in turn made the baseline tonic open
        probability far too large and compressed every derived PAM gain. It was caught by a
        known-phenomenology check: the scheme failed to reproduce the 2-3x leftward EC50
        shift that defines a classical benzodiazepine, predicting only 1.03-1.08x.

        Published GABA EC50 values come from applications of hundreds of ms, long enough
        for the peak to be genuinely reached at every concentration. 300 ms matches that.
        The max over the window handles both regimes: at saturating agonist the peak occurs
        within ~1 ms and then desensitises, at low agonist it arrives late.
        """
        p0 = np.zeros(5); p0[0] = 1.0
        Q = self.Q(conc_um)
        sol = solve_ivp(lambda t, p: Q.T @ p, (0.0, dur_ms), p0,
                        method="LSODA", rtol=1e-9, atol=1e-12, dense_output=True)
        ts = np.linspace(0.0, dur_ms, 600)
        return float(sol.sol(ts)[OPEN].max())

    def synaptic_waveform(self, t_end_ms: float = 200.0,
                          peak_um: float = SYNAPTIC_PEAK_UM,
                          clear_ms: float = SYNAPTIC_CLEAR_MS):
        """Open probability after a synaptic GABA transient (exponentially cleared pulse).

        Returns (t, Po). This is the phasic IPSC shape: its peak gives the amplitude and
        its late decay gives the deactivation tau.
        """
        def conc(t):
            return peak_um * np.exp(-t / clear_ms)

        def rhs(t, p):
            return self.Q(conc(t)).T @ p

        p0 = np.zeros(5); p0[0] = 1.0
        sol = solve_ivp(rhs, (0.0, t_end_ms), p0, method="LSODA",
                        rtol=1e-8, atol=1e-11, dense_output=True, max_step=1.0)
        t = np.linspace(0.0, t_end_ms, 4000)
        return t, sol.sol(t)[OPEN]

    def ipsc_metrics(self, **kw) -> dict:
        """Peak amplitude, charge, and deactivation tau of the phasic response.

        Deactivation tau is fitted to the decay between 90% and 10% of peak, which is the
        standard experimental window and avoids both the rising phase and the noise floor.
        """
        t, po = self.synaptic_waveform(**kw)
        i_pk = int(np.argmax(po)); pk = float(po[i_pk])
        if pk <= 1e-12:
            return dict(peak=0.0, tau=np.nan, charge=0.0)
        tail_t, tail_p = t[i_pk:], po[i_pk:]
        hi, lo = 0.90 * pk, 0.10 * pk
        m = (tail_p <= hi) & (tail_p >= lo)
        tau = np.nan
        if m.sum() > 5:
            # single-exponential fit in log space over the 90->10% window
            c = np.polyfit(tail_t[m] - tail_t[m][0], np.log(tail_p[m]), 1)
            if c[0] < 0:
                tau = float(-1.0 / c[0])
        return dict(peak=pk, tau=tau, charge=float(_trapz(po, t)))

    def dose_response(self, concs_um) -> np.ndarray:
        return np.array([self.po_peak(c) for c in concs_um])

    def ec50_um(self) -> float:
        """GABA EC50 for PEAK current, by interpolation on a log concentration grid."""
        c = np.logspace(-1, 4, 30)
        y = self.dose_response(c)
        ymax = y.max()
        if ymax <= 1e-12:
            return np.nan
        half = 0.5 * ymax
        i = int(np.argmax(y >= half))
        if i == 0:
            return float(c[0])
        x0, x1, y0, y1 = c[i - 1], c[i], y[i - 1], y[i]
        f = (half - y0) / max(1e-18, y1 - y0)
        return float(np.exp(np.log(x0) + f * (np.log(x1) - np.log(x0))))

    def po_max(self) -> float:
        return float(self.po_peak(1e4))


# ---------------------------------------------------------------------------- fitting
FIT_TARGETS = dict(ec50_um=20.0, po_max=0.75, tau_ms=15.0)
FIT_RANGES = dict(ec50_um=(10.0, 30.0), po_max=(0.70, 0.80), tau_ms=(10.0, 30.0))


def fit_scheme(verbose=True, pulse=None) -> Scheme:
    """Fit kon, koff, beta, alpha to the three macroscopic anchors.

    d and r (desensitisation) are NOT fitted here -- they shape the response to prolonged
    agonist, which none of the three anchors constrains. Leaving them at provisional values
    and saying so is more honest than fitting them to data that cannot identify them.
    """
    from scipy.optimize import least_squares

    pulse = pulse or {}

    def resid(x):
        s = Scheme(kon=10 ** x[0], koff=10 ** x[1], beta=10 ** x[2], alpha=10 ** x[3])
        m = s.ipsc_metrics(**pulse)
        tau = m["tau"] if np.isfinite(m["tau"]) else 1e3
        return [np.log(s.ec50_um() / FIT_TARGETS["ec50_um"]),
                (s.po_max() - FIT_TARGETS["po_max"]) / 0.05,
                np.log(tau / FIT_TARGETS["tau_ms"])]

    x0 = [np.log10(DEFAULT_RATES[k]) for k in ("kon", "koff", "beta", "alpha")]
    out = least_squares(resid, x0, bounds=([-4, -2, -2, -2.5], [-1, 1.2, 1.2, 1.0]),
                        xtol=1e-10, ftol=1e-10)
    s = Scheme(kon=10 ** out.x[0], koff=10 ** out.x[1],
               beta=10 ** out.x[2], alpha=10 ** out.x[3])
    if verbose:
        m = s.ipsc_metrics(**pulse)
        print("FIT to macroscopic anchors")
        print(f"  {'observable':<22}{'fitted':>10}{'target':>10}{'accepted range':>18}")
        for lab, got, key in (("GABA EC50 (uM)", s.ec50_um(), "ec50_um"),
                              ("max open prob", s.po_max(), "po_max"),
                              ("IPSC deact tau (ms)", m["tau"], "tau_ms")):
            lo, hi = FIT_RANGES[key]
            ok = "OK" if lo <= got <= hi else "OUT OF RANGE"
            print(f"  {lab:<22}{got:10.3f}{FIT_TARGETS[key]:10.3f}"
                  f"{f'{lo}-{hi}':>14}  {ok}")
        print(f"\n  rates (ms^-1; kon in uM^-1 ms^-1): kon={s.kon:.5f} koff={s.koff:.4f} "
              f"beta={s.beta:.4f} alpha={s.alpha:.4f}")
        print(f"  d={s.d} r={s.r}  (NOT fitted -- unconstrained by these anchors)")
    return s


# ------------------------------------------------------------- derived drug parameters
def derive(scheme: Scheme, affinity=1.0, gating=1.0,
           ambient_um=AMBIENT_UM, pulse=None) -> dict:
    """Derive the model's three GABA-A drug numbers from the kinetic scheme.

    Returns the quantities circuitpharm/cpg.py currently takes as independent free parameters:
      tonic_gain  -> what Drug.gaba_a_gain means for the EXTRASYNAPTIC/tonic conductance
      phasic_gain -> ... and for the SYNAPTIC peak. These differ, which is the point.
      tau_ratio   -> Drug.gaba_a_tau
      ceiling     -> gain at SATURATING agonist. NOTE: this is NOT the analogue of
                     Drug.gaba_a_efficacy_cap. It answers "what does the drug do to an
                     already-saturated receptor?" (answer, for affinity-type: nothing).
      tonic_headroom -> THIS is the quantity the model's efficacy cap is trying to be:
                     the maximum tonic gain available as the affinity shift grows without
                     bound, at the given ambient agonist. Bounded by steady-state
                     desensitisation rather than by Po_max. See the module docstring --
                     conflating these two was a real error, caught numerically.
    """
    pulse = pulse or {}
    base, mod = scheme, scheme.pam(affinity=affinity, gating=gating)
    b_t, m_t = base.po_tonic(ambient_um), mod.po_tonic(ambient_um)
    b_i, m_i = base.ipsc_metrics(**pulse), mod.ipsc_metrics(**pulse)
    b_sat, m_sat = base.po_max(), mod.po_max()
    f = lambda new, old: float(new / old) if old > 1e-15 else float("nan")
    sat = base.pam(affinity=1e5)          # the mechanism pushed to its limit
    return dict(tonic_gain=f(m_t, b_t),
                tonic_headroom=f(sat.po_tonic(ambient_um), b_t),
                phasic_headroom=f(sat.ipsc_metrics(**pulse)["peak"], b_i["peak"]),
                phasic_gain=f(m_i["peak"], b_i["peak"]),
                tau_ratio=f(m_i["tau"], b_i["tau"]),
                charge_ratio=f(m_i["charge"], b_i["charge"]),
                ceiling=f(m_sat, b_sat),
                ec50_shift=f(base.ec50_um(), mod.ec50_um()))


def calibrate_pam(scheme: Scheme, target_shift: float = 2.5,
                  kind: str = "affinity", pulse=None, ambient_um=AMBIENT_UM) -> float:
    """Find the ONE modulator parameter that reproduces the measured EC50 left-shift.

    This is the structural point of the whole module. The baseline scheme is fitted to
    BASELINE observables (EC50, max Po, IPSC tau). The drug then gets exactly ONE free
    parameter, calibrated against exactly ONE drug observable -- the leftward shift of the
    GABA dose-response, which is the defining and most robustly measured action of a
    classical benzodiazepine-site ligand (~2-3x).

    Everything else -- tonic gain, phasic gain, tau ratio, ceiling -- is then a PREDICTION,
    not a fit. That replaces three independently hand-set numbers with one calibrated
    number and four falsifiable outputs.
    """
    from scipy.optimize import brentq
    f = lambda x: (derive(scheme, pulse=pulse, ambient_um=ambient_um,
                          **{kind: x})["ec50_shift"] - target_shift)
    lo, hi = 1.0001, 1e4
    if f(hi) < 0:
        return float("nan")      # cannot reach the target shift by this mechanism
    # THE LOWER END NEEDS THE SAME CHECK AS THE UPPER END (review pass 7).
    #
    # Only `f(hi) < 0` was tested, so the upper, saturating end of the bracket was handled
    # and the lower end was not. The search domain is x >= 1.0001, i.e. POTENTIATION only,
    # so the smallest reachable shift is ~1.0. Ask for `target_shift = 1.0` (a neutral
    # ligand) or anything below it (a NEGATIVE allosteric modulator, which right-shifts
    # EC50 and therefore needs x < 1) and both bracket endpoints come out positive.
    # `brentq` requires opposite signs, so it raised a bare
    # `ValueError: f(a) and f(b) must have different signs` -- an scipy internal, not a
    # statement about pharmacology, from a function whose entire other failure path is a
    # documented NaN.
    #
    # NAMs are not merely unhandled here, they are OUT OF DOMAIN: this function searches
    # affinity/gating multipliers above 1. Returning NaN uses the convention the caller
    # already understands ("unreachable by this mechanism"; see evaluation.pool_gains,
    # which tests np.isfinite and falls back). Extending the bracket below 1 to cover NAMs
    # would be a feature, not a bug fix, and is deliberately not done here.
    if f(lo) > 0:
        return float("nan")
    return float(brentq(f, lo, hi, xtol=1e-6))


if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    s = fit_scheme()
    print(f"\n  SANITY: baseline open probability at ambient {AMBIENT_UM} uM "
          f"(tonic pool) = {s.po_tonic():.5f}")
    print(f"          baseline open probability at saturating GABA = {s.po_max():.3f}")
    print("          the tonic pool must sit FAR below saturation, else there is no")
    print("          headroom for a PAM and the derived ceiling will be spuriously tight.")

    print("\n" + "=" * 78)
    print("PREDICTION 1 — the synaptic/extrasynaptic asymmetry, from the kinetics alone")
    print("=" * 78)
    print(f"ambient (tonic) GABA {AMBIENT_UM} uM vs synaptic peak "
          f"{SYNAPTIC_PEAK_UM:.0f} uM; EC50 {s.ec50_um():.1f} uM")
    print(f"\n{'affinity x':>10}{'EC50 shift':>12}{'TONIC gain':>12}"
          f"{'PHASIC gain':>13}{'tau ratio':>11}{'charge':>9}{'ceiling':>9}")
    print("-" * 78)
    for aff in (1.0, 1.5, 2.0, 3.0, 5.0, 10.0):
        d = derive(s, affinity=aff)
        print(f"{aff:10.1f}{d['ec50_shift']:12.2f}{d['tonic_gain']:12.2f}"
              f"{d['phasic_gain']:13.2f}{d['tau_ratio']:11.2f}"
              f"{d['charge_ratio']:9.2f}{d['ceiling']:9.2f}")
    print("\nREAD: an affinity-type PAM potentiates the TONIC current strongly, barely")
    print("moves the already-saturated SYNAPTIC peak, and prolongs decay -- the measured")
    print("benzodiazepine signature. 'ceiling' is its gain at saturating agonist: near")
    print("1.0, i.e. it can do NOTHING to a saturated receptor. That is the overdose-")
    print("protection property, derived rather than assumed.")

    print("\n" + "=" * 78)
    print("PREDICTION 2 — gating-type modulation (neurosteroid-like) should differ")
    print("=" * 78)
    print(f"{'gating x':>10}{'EC50 shift':>12}{'TONIC gain':>12}"
          f"{'PHASIC gain':>13}{'tau ratio':>11}{'charge':>9}{'ceiling':>9}")
    print("-" * 78)
    for g in (1.0, 1.5, 2.0, 3.0, 5.0):
        d = derive(s, gating=g)
        print(f"{g:10.1f}{d['ec50_shift']:12.2f}{d['tonic_gain']:12.2f}"
              f"{d['phasic_gain']:13.2f}{d['tau_ratio']:11.2f}"
              f"{d['charge_ratio']:9.2f}{d['ceiling']:9.2f}")
    print("\nVALIDATION TARGET: the project hand-sets nonselective_bz ceiling=2.5 and")
    print("neurosteroid ceiling=6.0. If gating-type modulation shows a materially higher")
    print("ceiling than affinity-type here, that ordering is PREDICTED, not assumed.")

    print("\n" + "=" * 78)
    print("CALIBRATION — one drug parameter, one drug observable, four predictions")
    print("=" * 78)
    aff = calibrate_pam(s, target_shift=2.5, kind="affinity")
    print(f"affinity multiplier reproducing a 2.5x EC50 left-shift: {aff:.3f}")
    if not np.isfinite(aff):
        print("  UNREACHABLE by an affinity-only mechanism in this scheme.")
    else:
        d = derive(s, affinity=aff)
        print(f"\n{'prediction':<26}{'value':>8}   known benzodiazepine phenomenology")
        print("-" * 78)
        checks = [
            ("phasic peak gain", d["phasic_gain"], (0.95, 1.25),
             "little change in mIPSC amplitude"),
            ("phasic decay tau ratio", d["tau_ratio"], (1.3, 2.5),
             "mIPSC decay prolonged ~1.5-2x"),
            ("tonic current gain", d["tonic_gain"], (1.5, 6.0),
             "tonic current strongly potentiated"),
            ("ceiling (gain at sat.)", d["ceiling"], (1.0, 1.3),
             "no effect on a saturated receptor"),
        ]
        n_ok = 0
        for lab, val, (lo, hi), note in checks:
            ok = lo <= val <= hi
            n_ok += ok
            print(f"{lab:<26}{val:8.2f}   {'OK ' if ok else 'OUT'}  "
                  f"[{lo}-{hi}]  {note}")
        print(f"\n  {n_ok}/4 predictions inside the expected range, from ONE "
              f"calibrated parameter.")

    print("\n" + "=" * 78)
    print("COMPARISON with the hand-set parameters now in the model")
    print("=" * 78)
    print("the model currently uses, for a 'clinically sedative BZ dose':")
    print("  gaba_a_gain 2.0 applied to BOTH tonic and synaptic pools")
    print("  gaba_a_tau  1 + 0.6*(2-1) = 1.60")
    print("  cap         2.5")
    if np.isfinite(aff):
        d = derive(s, affinity=aff)
        print(f"\nthe kinetic scheme, calibrated only on the EC50 shift, instead says:")
        print(f"  tonic gain   {d['tonic_gain']:.2f}  vs synaptic gain "
              f"{d['phasic_gain']:.2f}  -- these MUST differ; the model uses one number")
        print(f"  tau ratio    {d['tau_ratio']:.2f}  (model hand-sets 1.60)")
        print(f"  ceiling      {d['ceiling']:.2f}  (model hand-sets 2.5)")
        print("\nThe ceiling is the consequential disagreement. The model's 2.5 means a PAM")
        print("can multiply the GABA-A effect 2.5-fold at saturation; the kinetics say an")
        print("affinity-type modulator can do essentially NOTHING once the receptor is")
        print("saturated. If that holds, the model has been OVERSTATING how much subjective")
        print("effect a PAM can deliver -- which bears directly on the reachability")
        print("finding, where the a5 arm needed gain 2.54 against a ceiling of 2.5.")
