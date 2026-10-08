"""Conductance-based single-compartment neuron, and the parameter registry that owns it.

ROADMAP LINK 4. The LIF population in `cpg.Pop` cannot carry a voltage-gated mechanism:
V is reset at threshold, so a slow voltage-gated gate equilibrates at the subthreshold mean
and never moves (recurring error E7, measured h-gate span 0.01 against a needed 0.7). That
forced spike-triggered adaptation as a lumped stand-in for I_NaP inactivation PLUS
Ca-dependent K current PLUS synaptic depression -- documented, defensible, and not the
biophysics.

WHAT THE LIF ACTUALLY COSTS, measured rather than argued (`scripts/diag_substrate_limits.py`;
the full case is knowledge/05-design-conductance-substrate.md section 1). On the shipped
preBotC the LIF holds V in [-71, -44] mV, and therefore:

  * NMDA conductance sits in near-permanent Mg2+ block -- mean relief 0.063, dynamic range
    4.5x against the 15.5x a spiking cell has. `cpg.py` says "voltage dependence is why NMDA
    block is state-dependent"; on that substrate it is not.
  * The PHASIC GABA-A pool's driving force is truncated exactly when it matters. Mean GABA-A
    driving force is 10.4 mV; during a real burst it would be up to 95 mV. Both pools share
    one V, so at rest both are right -- the asymmetry is in TIMING. Tonic is always on at the
    resting driving force; phasic arrives correlated with the burst.
  * Excitation never self-limits: the glutamate driving force never approaches zero.

The objection that the hand-tuned weights absorbed all this is correct about MAGNITUDE and
wrong about DEPENDENCE. A scalar weight rescales a mean; it cannot turn a 4.5x relief range
into 15.5x, nor make one pool's driving force swing 9x while the other's stays flat, because
both pools share a single weight-independent E_rev and a single V.

So this module is not "more biophysics for its own sake". It exists so the question can be
asked: is the selectivity ranking -- this project's one VALIDATED result -- substrate
independent as well as calibration independent? Nobody has checked.

WHAT THIS MODULE DOES NOT DO. It does not create a respiratory calibration anchor: that is
blocked on an OBSERVABLE (whole-body ventilation is wrong for an isolated preBotC), not on
what the neurons are made of. It makes the developmental confound WORSE, not better -- see
`PROVENANCE AND AGE` on CellParams.
"""
from __future__ import annotations

import math
import zlib
from dataclasses import dataclass, field, replace

import numpy as np


# ======================================================================== parameter sets
@dataclass(frozen=True)
class CellParams:
    """One published cell, with the provenance that makes it citable.

    PROVENANCE AND AGE IS LOAD-BEARING, not bookkeeping. Adopting a parameter set adopts its
    preparation. Butera-Rinzel-Smith is neonatal rodent in vitro (~P0-P4), and this project
    has already been bitten once by the developmental chloride switch (NKCC1/KCC2 crossover
    at ~P11-P12, over-corrected on and then corrected back). Two consequences that must
    travel with these numbers:

      1. E_GABA should be the NEONATAL value, which is more depolarised than the -75 mV in
         `cpg.E_REV`, because the chloride switch has not happened yet. So GABA-A here is
         LESS inhibitory, possibly shunting-only. The choice of parameter set is therefore
         NOT pharmacologically neutral.
      2. This cell is at the wrong age for the >P12 muscimol concentration-response anchor
         the respiratory axis is blocked on. Link 4 moves that anchor FURTHER away.

    `scale_ref` is the input conductance the synaptic weight table must have been tuned
    against. It exists to make a mixed cell/weight pair unconstructible -- see `ParamSet`.
    """
    name: str
    provenance: str                 # publication, species, age, preparation
    C: float                        # pF
    g_na: float                     # nS, fast transient Na
    g_k: float                      # nS, delayed rectifier K
    g_nap: float                    # nS, persistent Na
    g_l: float                      # nS, leak
    e_na: float                     # mV
    e_k: float                      # mV
    e_l: float                      # mV
    # (theta, sigma) for each gate, plus tau_bar for the two that are not instantaneous.
    # x_inf(V) = 1 / (1 + exp((V - theta) / sigma));  tau(V) = tau_bar / cosh((V-theta)/(2 sigma))
    m: tuple                        # fast Na activation, INSTANTANEOUS
    n: tuple                        # K activation (and Na inactivation via 1-n)
    mp: tuple                       # persistent Na activation, INSTANTANEOUS
    h: tuple                        # persistent Na inactivation -- the slow variable
    tau_n: float                    # ms, tau_bar for n
    tau_h: float                    # ms, tau_bar for h
    dt_max: float                   # ms, largest integration substep that passes E13
    v_detect: float = -20.0         # mV, spike detection crossing (see E16)
    v_reset_lock: float = -30.0     # mV, must fall below this to re-arm the detector
    sigma_i: float = 10.0           # pA/sqrt(ms), membrane current noise

    @property
    def g_input(self) -> float:
        """Resting input conductance, the scale every nS synaptic weight is relative to."""
        return self.g_l


# BUTERA-RINZEL-SMITH 1999 MODEL 1.
#
# EVERY NUMBER BELOW WAS READ OFF A MACHINE-READABLE IMPLEMENTATION, NOT RECALLED. That
# distinction matters here: this project has already scored an invented range as if it were
# data, and the design document flagged these values as UNVERIFIED precisely so they would be
# checked before use. Source: the curated CellML encoding of the model
# (models.cellml.org/exposure/293a909eeeca07d0ca8ad583839eb0bc).
#
# The check caught one error in the recalled set: E_L is -57.5 mV, not -65 mV. That is not
# cosmetic -- E_L is the bifurcation parameter of this model, the knob that moves the cell
# quiescent -> bursting -> tonic spiking. Starting from -65 would have put the cell in the
# wrong regime and the rhythm would then have been retuned around a wrong resting drive,
# which is the exact shape of this project's recurring failures.
BRS1999_MODEL1 = CellParams(
    name="butera_rinzel_smith_1999_model1",
    provenance=("Butera, Rinzel & Smith 1999, J Neurophysiol 82:382-397, model 1 "
                "(I_NaP-h); neonatal rodent preBotC in vitro, ~P0-P4. Parameters read "
                "from the curated CellML encoding, not from recall."),
    C=21.0,
    g_na=28.0, g_k=11.2, g_nap=2.8, g_l=2.8,
    e_na=50.0, e_k=-85.0, e_l=-57.5,
    m=(-34.0, -5.0),
    n=(-29.0, -4.0),
    mp=(-40.0, -6.0),
    h=(-48.0, 6.0),
    tau_n=10.0,
    tau_h=10000.0,
    dt_max=0.05,
)

CELLS = {BRS1999_MODEL1.name: BRS1999_MODEL1}


class ProvenanceMismatch(ValueError):
    """Raised when a cell is paired with a weight table tuned against a different cell."""


@dataclass(frozen=True)
class ParamSet:
    """A cell and the synaptic weights tuned against it, as ONE unconstructible-if-mixed unit.

    WHY THIS GUARD EXISTS. The LIF cell is C=200 pF, g_L=10 nS. This cell is C=21 pF,
    g_L=2.8 nS -- a ~3.6x difference in input conductance and ~10x in capacitance. Every
    synaptic weight in cpg/resp/rg2/circuit is in nS, hand-tuned against the LIF cell. Swap
    the cell alone and every weight is wrong by roughly that factor; the network will not
    oscillate, or will oscillate for the wrong reason, and NOTHING RAISES.

    That is the single most likely defect in this migration (predicted as E14 before any code
    was written), so it is made structurally impossible rather than guarded by a comment.
    """
    cell: CellParams
    weights: dict = field(default_factory=dict)
    weights_provenance: str = ""

    def __post_init__(self):
        if self.weights and not self.weights_provenance:
            raise ProvenanceMismatch(
                f"weights supplied for cell {self.cell.name!r} with no provenance. A weight "
                f"table is only meaningful relative to the cell it was tuned against; an "
                f"untagged one cannot be checked and must not be used.")
        if self.weights and self.cell.name not in self.weights_provenance:
            raise ProvenanceMismatch(
                f"weight table provenance {self.weights_provenance!r} does not name cell "
                f"{self.cell.name!r}. nS weights are relative to the cell's input "
                f"conductance ({self.cell.g_input} nS here); pairing them across cells is a "
                f"silent order-of-magnitude error, not a configuration choice.")


# ============================================================================ the cell
def _x_inf(V, gate):
    theta, sigma = gate
    return 1.0 / (1.0 + np.exp((V - theta) / sigma))


def _tau(V, gate, tau_bar):
    theta, sigma = gate
    return tau_bar / np.cosh((V - theta) / (2.0 * sigma))


@dataclass
class CondPop:
    """Conductance-based population. Duck-type compatible with `cpg.Pop`.

    Same call signature as the LIF: `step(dt, g, E, Idrive, rng)`, exposing `V`, `spk`,
    `rate` and `a`. A circuit can therefore swap substrates without touching its stepping
    loop -- which is what makes the LIF-vs-conductance comparison practical, and that
    comparison is the actual deliverable (see module docstring).

    THREE STATE VARIABLES, not a full Hodgkin-Huxley model: V, n and h. Fast-Na activation
    `m` and persistent-Na activation `mp` are instantaneous, and Na inactivation is carried
    by `(1-n)`. That economy is why this costs roughly 20x the LIF rather than 60x.

    ADAPTATION DEFAULTS OFF (`g_adapt=0.0`), deliberately, and that is a behavioural change
    from the LIF rather than an oversight. Spike-triggered adaptation was a STAND-IN for
    I_NaP inactivation; now that `h` is real, leaving adaptation at its LIF-tuned strength
    double-counts burst termination. The symptom would be bursts terminating slightly early,
    i.e. a plausible duty cycle -- predicted as E18.
    """
    n: int
    name: str
    params: CellParams = BRS1999_MODEL1
    # pA, standing applied current. POSITIVE DEPOLARISES.
    #
    # THIS SIGN DIFFERS FROM THE SOURCE ENCODING and the difference is silent. The CellML
    # writes dV/dt = -(i_NaP + i_Na + i_K + i_L + i_tonic_e + i_app)/C, putting i_app INSIDE
    # the negated sum, so there a positive i_app HYPERPOLARISES. Here it is outside, matching
    # `cpg.Pop`, where `Idrive` is added to an inward-positive current sum and a positive
    # value depolarises.
    #
    # Matching the LIF is the right choice, because the whole point of this class is that a
    # circuit can swap substrates without changing its stepping loop -- and every circuit in
    # this package already passes a positive `Idrive` meaning excitation. But anyone
    # comparing a number here against the published figures must flip the sign, and a
    # sign error on the bifurcation parameter would move the cell between regimes while
    # looking like a plausible result.
    i_app: float = 0.0
    g_adapt: float = 0.0            # nS per spike; see the docstring -- off by design
    tau_adapt: float = 450.0        # ms
    nap: bool = True                # kept for Pop compatibility; False zeroes g_NaP
    sigma_i: float | None = None    # pA/sqrt(ms); None -> params.sigma_i

    def __post_init__(self):
        p = self.params
        # Stable seed, as in cpg.Pop: `hash(str)` is randomised per process via
        # PYTHONHASHSEED, which made initial voltages differ between runs even under an
        # explicit seed. crc32 is stable across processes, versions and platforms.
        rng = np.random.default_rng(zlib.crc32(self.name.encode()) % (2**32))
        self.V = p.e_l + rng.normal(0.0, 3.0, self.n)
        self.nn = _x_inf(self.V, p.n)
        self.h = _x_inf(self.V, p.h)
        self.a = np.zeros(self.n)
        self.spk = np.zeros(self.n, bool)
        self.rate = 0.0
        self._armed = np.ones(self.n, bool)     # spike detector re-armed below v_reset_lock
        self.EK = p.e_k                         # Pop compatibility
        self.EL = p.e_l
        self.C = p.C

    # ---- compatibility shims ---------------------------------------------------------
    @property
    def Vth(self):
        """There is no threshold in a conductance model. Exposed as the DETECTION voltage so
        code that reports `p.Vth` keeps working, but it is not a reset threshold: V genuinely
        exceeds it, which is the whole point (see test_E7 inverted)."""
        return self.params.v_detect

    # ---- dynamics --------------------------------------------------------------------
    def _currents(self, V, g, E, raw):
        p = self.params
        m3 = _x_inf(V, p.m) ** 3
        i = (p.g_na * m3 * (1.0 - self.nn) * (V - p.e_na)
             + p.g_k * self.nn ** 4 * (V - p.e_k)
             + p.g_l * (V - p.e_l))
        if self.nap:
            i = i + p.g_nap * _x_inf(V, p.mp) * self.h * (V - p.e_na)
        i = i + self.a * (V - p.e_k)                      # adaptation, outward
        for k, gk in g.items():
            gg = gk
            if k in raw and k == "nmda":
                # MG2+ BLOCK APPLIED AT THIS CELL'S OWN CURRENT V, EVERY SUBSTEP.
                #
                # This is the difference between the upgrade working and silently doing
                # nothing for NMDA. The LIF interface hands over a conductance already
                # evaluated at the PRE-STEP voltage, which is harmless when V moves a few mV
                # per step and destroys the entire voltage-dependent relief when V sweeps
                # 85 mV through a spike. The failure mode would have been "the upgrade did
                # not change the NMDA result" -- which reads as a reassuring robustness
                # check and is in fact the bug (predicted as E15).
                gg = gk / (1.0 + 0.28 * np.exp(-0.062 * V))
            i = i + gg * (V - E[k])
        return i

    def step(self, dt, g, E, Idrive, rng, raw=()):
        """Advance by `dt` ms, substepping internally so callers need not change their loop.

        `g`, `E`: dicts of conductance (nS) and reversal (mV) per receptor, as in `cpg.Pop`.
        `raw`: receptor keys in `g` whose conductance is UNEVALUATED, so this cell applies
        the voltage-dependent factor itself at its own V. Pass `raw=("nmda",)` to get the
        Mg2+ relief the LIF substrate cannot express.
        """
        p = self.params
        n_sub = max(1, int(math.ceil(dt / p.dt_max)))
        hs = dt / n_sub
        sig = p.sigma_i if self.sigma_i is None else self.sigma_i
        fired = np.zeros(self.n, bool)

        for _ in range(n_sub):
            i_ion = self._currents(self.V, g, E, raw)
            # Noise enters as a CURRENT (pA/sqrt(ms)), not as a voltage kick. The LIF adds
            # sigma*sqrt(dt) directly in mV, whose physical meaning depends on C, so the
            # same parameter means different things on cells of different capacitance.
            # Applied only here; the LIF path is left byte-identical on purpose (A7).
            #
            # The Wiener increment is sigma*xi*sqrt(hs), divided by C to reach mV. Writing
            # it as a current amplitude sigma/sqrt(hs) multiplied back by hs/C is the same
            # number by two cancelling factors, and I had those factors cancel wrongly on
            # the first pass -- leaving noise that did not scale with the substep at all,
            # so halving dt would have quietly changed the noise amplitude and contaminated
            # the E13 dt-convergence test with the very thing it exists to detect.
            dv_det = hs * (-i_ion + self.i_app + Idrive) / p.C
            dv_noise = sig * math.sqrt(hs) * rng.standard_normal(self.n) / p.C
            self.V = self.V + dv_det + dv_noise

            # Rush-Larsen (exponential Euler) on the gates: exact for a linear relaxation
            # with tau and x_inf frozen over the substep, so the gates are unconditionally
            # stable and dt is set by the accuracy of the spike upstroke, not by gate
            # stiffness. tau_h is 10 s, so forward Euler would be wasteful here, not unstable.
            for attr, gate, tau_bar in (("nn", p.n, p.tau_n), ("h", p.h, p.tau_h)):
                x = getattr(self, attr)
                xinf = _x_inf(self.V, gate)
                tau = _tau(self.V, gate, tau_bar)
                setattr(self, attr, xinf + (x - xinf) * np.exp(-hs / tau))

            # SPIKE DETECTION WITH A RE-ARM LOCKOUT (predicted as E16). A bare upward
            # crossing of v_detect counts one broad spike several times as V jitters across
            # the line, inflating every firing rate -- and an inflated rate looks like a
            # plausible number. A detection is only possible once V has fallen back below
            # v_reset_lock.
            crossed = self._armed & (self.V >= p.v_detect)
            fired |= crossed
            self._armed = np.where(crossed, False, self._armed)
            self._armed = np.where(self.V < p.v_reset_lock, True, self._armed)

            if self.g_adapt:
                self.a = self.a + self.g_adapt * crossed
            self.a = self.a - hs * self.a / self.tau_adapt

        self.spk = fired
        inst = fired.sum() / self.n / (dt * 1e-3)          # Hz
        self.rate += dt * (inst - self.rate) / 20.0        # 20 ms low-pass, as in cpg.Pop
        return self.spk

    # ---- diagnostics -----------------------------------------------------------------
    def gate_state(self) -> dict:
        """Current gating variables. `h` is the one E7 said a LIF cannot move."""
        return dict(n=np.array(self.nn), h=np.array(self.h), V=np.array(self.V))


def isolated(params: CellParams = BRS1999_MODEL1, **kw) -> CondPop:
    """A single cell with no noise -- the deterministic preparation the acceptance tests use."""
    return CondPop(n=1, name=kw.pop("name", "isolated"), params=params, sigma_i=0.0, **kw)


# ======================================================================= characterisation
def spike_regime(spike_times_ms, settle_ms=2000.0, min_gap_ms=50.0,
                 gap_factor=5.0) -> dict:
    """Classify a spike train as quiescent / tonic / bursting.

    A SETTLING PERIOD IS DISCARDED, and that is not cosmetic. The first version of this
    analysis kept the whole trace, and the startup transient -- one long interval while `h`
    relaxes from its initial value -- registered as a burst boundary. A cell firing tonically
    at 100 Hz was reported as "2 bursts". Since the regime classification is what the
    excitability-sequence test keys on, that would have scored a tonic cell as bursting at
    high drive and inverted the published sequence.

    A burst boundary is a gap exceeding both `min_gap_ms` and `gap_factor` x the median
    interval, so neither an absolute nor a relative criterion can be fooled alone: a purely
    relative threshold finds "bursts" in a regular train whose intervals jitter, and a purely
    absolute one misses bursting at a frequency it was not chosen for.
    """
    s = np.asarray(spike_times_ms, float)
    s = s[s >= settle_ms]
    if s.size == 0:
        return dict(regime="quiescent", n_bursts=0, n_spikes=0, rate_hz=0.0,
                    spikes_per_burst=0.0, burst_period_ms=float("nan"))
    span_s = (s[-1] - s[0]) / 1000.0
    rate = s.size / span_s if span_s > 0 else 0.0
    if s.size < 3:
        return dict(regime="quiescent", n_bursts=0, n_spikes=int(s.size), rate_hz=rate,
                    spikes_per_burst=0.0, burst_period_ms=float("nan"))

    isi = np.diff(s)
    thr = max(min_gap_ms, gap_factor * float(np.median(isi)))
    boundary = isi > thr
    # group spikes into bursts
    starts = np.r_[0, np.flatnonzero(boundary) + 1]
    ends = np.r_[np.flatnonzero(boundary), s.size - 1]
    sizes = ends - starts + 1
    multi = sizes >= 2
    n_bursts = int(multi.sum())
    spb = float(sizes[multi].mean()) if n_bursts else 0.0
    onsets = s[starts[multi]]
    period = float(np.mean(np.diff(onsets))) if onsets.size > 2 else float("nan")

    # Bursting needs several MULTI-spike groups. A tonic train is one group; a train of
    # singletons is not bursting however regularly spaced, which is why spikes-per-burst
    # must exceed 2 rather than merely being defined.
    regime = "bursting" if (n_bursts >= 3 and spb > 2.0) else "tonic"
    return dict(regime=regime, n_bursts=n_bursts, n_spikes=int(s.size), rate_hz=rate,
                spikes_per_burst=spb, burst_period_ms=period)


def run_isolated(seconds=20.0, dt=0.1, params: CellParams = BRS1999_MODEL1,
                 freeze_h=None, seed=0, record_every=5, **kw) -> dict:
    """Drive one deterministic cell open-loop and return its spike train and gate trace.

    `freeze_h` clamps the slow I_NaP inactivation gate at a constant value, which is the
    manipulation that tests whether `h` is what terminates the burst. It has no biological
    counterpart -- it is a model dissection, and the strongest available check that the
    rhythm comes from the mechanism claimed rather than from something incidental.
    """
    c = isolated(params=params, **kw)
    rng = np.random.default_rng(seed)
    spk, H, V, T = [], [], [], []
    for i in range(int(round(seconds * 1000.0 / dt))):
        if freeze_h is not None:
            c.h[:] = freeze_h
        if c.step(dt, {}, {}, 0.0, rng)[0]:
            spk.append(i * dt)
        if record_every and i % record_every == 0:
            H.append(float(c.h[0])); V.append(float(c.V[0])); T.append(i * dt)
    H = np.asarray(H); V = np.asarray(V)
    out = dict(spikes_ms=np.asarray(spk), h=H, V=V, t_ms=np.asarray(T),
               h_span=float(H.max() - H.min()) if H.size else 0.0,
               v_max=float(V.max()) if V.size else float("nan"))
    out.update(spike_regime(out["spikes_ms"]))
    return out
