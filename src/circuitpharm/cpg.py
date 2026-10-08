"""Spinal half-centre CPG with receptor-resolved pharmacology.

This is where the drug acts. Deliberately NOT in the policy network: conductances here
have real units (nS, mV, ms), so a GABA-A PAM or an NMDA antagonist is expressed as the
physical quantity it actually changes, not as a dimensionless activation-function gain.

Architecture (two-layer, after Rybak/McCrea-style rat hindlimb CPG models), per
antagonist joint pair:

    RG-F <--| InRG-F        RG-E <--| InRG-E       rhythm generator, I_NaP-driven
      |  reciprocal inhibition (glycinergic)  |
      v                                        v
    PF-F                                     PF-E   pattern formation
      |                                        |
      v                                        v
    Mn-F <--| IaIn-F  ...  IaIn-E |--> Mn-E         motoneuron pools
      |____> Rc-F |                | Rc-E <____|     Renshaw recurrent inhibition
      ^                                        ^
    Ia afferent (spindle, length+velocity) ----+     monosynaptic, AMPA+NMDA

Receptor map -- the pharmacological interface:
  GABA-A   : tonic + phasic inhibition onto RG / PF / Mn        <- GABA-A PAM acts here
  Glycine  : IaIn reciprocal inhibition, Renshaw recurrent      <- ethanol potentiates
  NMDA     : Ia->Mn monosynaptic, RG excitatory drive           <- NMDA antagonist acts here
  AMPA     : fast excitation (not a drug target, needed for function)

SIMPLIFICATIONS, stated honestly:
  - LIF neurons with conductance synapses; only the RG carries I_NaP (the rhythmogenic
    current). No dendritic compartments, no calcium dynamics, no short-term plasticity.
  - Population sizes are small (10-30) -- enough for rate coding, not for finite-size
    network effects.
  - Ia afferent model is a Prochazka-style static rate law, not a Mileusnic/Loeb spindle.
  - No gamma motoneurons, so no fusimotor gain control.
  - Connectivity weights are hand-tuned to produce alternating rhythm, NOT fitted to
    rat data. They are a working substrate, not a validated model. Validation against
    the four motor anchors (see knowledge/03-model-spec.md) has NOT been done.
"""
import zlib
from dataclasses import dataclass, field
import numpy as np


# ----------------------------------------------------------------------------- drug
@dataclass
class Drug:
    """Receptor-level perturbation. All fields default to drug-free (identity)."""
    # GABA-A positive allosteric modulation. A PAM raises peak conductance AND
    # prolongs decay; both matter and they are not interchangeable.
    #
    # THE TWO POOLS ARE NOT EQUIVALENT (circuitpharm/gabaa_kinetics.py). A PAM acting by
    # increasing apparent agonist affinity barely touches the SYNAPTIC peak, because the
    # cleft transient is already near-saturating (derived phasic gain ~1.06), while it
    # potentiates the TONIC/extrasynaptic current strongly, because ambient GABA (~0.4 uM)
    # sits far below EC50 (~20 uM) (derived tonic gain ~7.2 at clinical BZ potency). The
    # model previously applied ONE multiplier to both, which is wrong by a factor of ~7 on
    # whichever pool it was not calibrated against.
    #
    # `gaba_a_gain_tonic = None` means "follow gaba_a_gain", reproducing the old behaviour
    # exactly so prior results stay reproducible. Use Drug.from_kinetics() for the derived
    # two-pool version.
    gaba_a_gain: float = 1.0        # multiplier on PHASIC (synaptic) peak g_GABA-A
    gaba_a_gain_tonic: float | None = None   # multiplier on TONIC g_GABA-A
    gaba_a_tau: float = 1.0         # multiplier on GABA-A decay tau
    gaba_a_efficacy_cap: float = 2.5  # ceiling on the phasic gain
    gaba_a_cap_tonic: float | None = None    # ceiling on the tonic gain; None -> follow
    # Glycine receptor potentiation (ethanol does this; most tool compounds do not)
    glyr_gain: float = 1.0
    # NMDA antagonism, resolved by subunit so selectivity can be expressed
    nmda_block: float = 0.0         # fraction of NMDA conductance blocked (0-1)
    glun2b_fraction: float = 0.7    # fraction of local NMDA that is GluN2B
    glun2b_selectivity: float = 1.0 # 1.0 = GluN2B-selective, 0.0 = non-selective

    def gaba_scale(self) -> float:
        """PHASIC (synaptic) GABA-A conductance multiplier, capped."""
        return min(self.gaba_a_gain, self.gaba_a_efficacy_cap)

    def gaba_scale_tonic(self) -> float:
        """TONIC (extrasynaptic) GABA-A conductance multiplier, capped.

        Falls back to the phasic value when unset, which is the legacy one-pool
        behaviour. See the field comments above for why the two differ.
        """
        g = self.gaba_a_gain if self.gaba_a_gain_tonic is None else self.gaba_a_gain_tonic
        cap = (self.gaba_a_efficacy_cap if self.gaba_a_cap_tonic is None
               else self.gaba_a_cap_tonic)
        return min(g, cap)

    @classmethod
    def from_kinetics(cls, affinity=None, ec50_shift=2.5, ambient_um=0.4,
                      pulse=None, scheme=None, **kw) -> "Drug":
        """Build a Drug whose GABA-A numbers are DERIVED from the Markov gating scheme
        rather than hand-set.

        One mechanistic input (the ligand's EC50 left-shift, or its affinity multiplier
        directly) replaces three independent parameters. The efficacy caps are set to the
        pool-specific HEADROOM -- the most that mechanism can achieve in that pool at the
        given ambient agonist -- instead of a single invented 2.5.

        Imported lazily: scipy is needed only on this path, so the circuits stay importable
        without it.
        """
        from .gabaa_kinetics import fit_scheme, derive, calibrate_pam
        pulse = pulse or dict(peak_um=3000.0, clear_ms=1.00)   # a validated passing cell
        s = scheme or fit_scheme(verbose=False, pulse=pulse)
        if affinity is None:
            affinity = calibrate_pam(s, ec50_shift, "affinity",
                                     pulse=pulse, ambient_um=ambient_um)
        d = derive(s, affinity=affinity, ambient_um=ambient_um, pulse=pulse)
        return cls(gaba_a_gain=d["phasic_gain"],
                   gaba_a_gain_tonic=d["tonic_gain"],
                   gaba_a_tau=d["tau_ratio"],
                   gaba_a_efficacy_cap=d["phasic_headroom"],
                   gaba_a_cap_tonic=d["tonic_headroom"], **kw)

    def nmda_scale(self) -> float:
        """Surviving fraction of NMDA conductance.

        A GluN2B-selective antagonist can only block the GluN2B pool, so its maximum
        effect is bounded by glun2b_fraction. That bound IS the selectivity window --
        in forebrain (GluN2B-rich) it blocks a lot, in brainstem (GluN2D-dominant,
        low glun2b_fraction) it blocks little. This is the mechanism the whole
        margin hypothesis rests on.
        """
        s = np.clip(self.glun2b_selectivity, 0.0, 1.0)
        accessible = s * self.glun2b_fraction + (1.0 - s) * 1.0
        return float(1.0 - self.nmda_block * accessible)


# ------------------------------------------------------------------------ population
@dataclass
class Pop:
    """LIF population with conductance synapses. Voltages in mV, conductances in nS."""
    n: int
    name: str
    nap: bool = False               # persistent sodium (rhythmogenic)
    C: float = 200.0                # pF
    gL: float = 10.0                # nS
    EL: float = -65.0               # mV
    Vth: float = -50.0
    Vreset: float = -65.0
    tref: float = 2.0               # ms
    gnap: float = 14.0              # nS, only if nap (regenerative excitability)
    g_adapt: float = 0.55           # nS added per spike (burst termination)
    tau_adapt: float = 450.0        # ms, adaptation decay
    EK: float = -85.0               # mV, adaptation reversal
    sigma: float = 0.20             # noise, mV/sqrt(ms)

    def __post_init__(self):
        # STABLE hash, not Python's. `hash(str)` is randomised per interpreter process via
        # PYTHONHASHSEED, so `abs(hash(self.name))` gave DIFFERENT initial membrane
        # voltages on every run even when the caller passed an explicit seed -- measured
        # V0[0] = -71.05, -62.73, -69.90 across three processes for the same population.
        # That makes every result irreproducible between runs, which for a tool whose
        # output is meant to be checked by others is fatal rather than cosmetic.
        # crc32 is stable across processes, versions and platforms.
        rng = np.random.default_rng(zlib.crc32(self.name.encode()) % (2**32))
        self.V = self.EL + rng.normal(0, 3.0, self.n)
        self.ref = np.zeros(self.n)
        self.a = np.zeros(self.n)            # spike-triggered adaptation (nS)
        self.spk = np.zeros(self.n, bool)
        self.rate = 0.0                      # low-pass firing rate, Hz

    # Burst mechanism. I_NaP (m-gated, instantaneous) supplies the regenerative
    # subthreshold depolarisation that lets the cell fire on SUB-rheobase drive. Burst
    # TERMINATION is a spike-triggered adaptation conductance, not voltage-gated I_NaP
    # inactivation.
    #
    # Why: in a LIF neuron V is reset at threshold and never reaches spike voltages, so
    # a voltage-gated h gate equilibrates at the subthreshold mean (h stays ~0.44, span
    # ~0.01) and never deeply inactivates. Voltage-gated inactivation cannot work on a
    # truncated LIF trajectory. The spike-triggered form is the standard substitute and
    # is a defensible lumped stand-in for I_NaP inactivation plus Ca-dependent K current
    # plus synaptic depression -- all activity-dependent. It is a SIMPLIFICATION, not
    # the biophysics.
    def _inap_m(self):
        return 1.0 / (1.0 + np.exp(-(self.V + 47.1) / 3.1))

    def step(self, dt, g, E, Idrive, rng):
        """g, E: dicts of conductance (nS) and reversal (mV) per receptor."""
        I = self.gL * (self.EL - self.V) + Idrive
        for k, gk in g.items():
            I = I + gk * (E[k] - self.V)
        if self.nap:
            I = I + self.gnap * self._inap_m() * (55.0 - self.V)
        I = I + self.a * (self.EK - self.V)          # adaptation (outward)
        self.V += dt * I / self.C
        self.V += self.sigma * np.sqrt(dt) * rng.standard_normal(self.n)
        live = self.ref <= 0
        self.spk = live & (self.V >= self.Vth)
        self.V[self.spk] = self.Vreset
        self.ref[self.spk] = self.tref
        self.ref[~self.spk] = np.maximum(0.0, self.ref[~self.spk] - dt)
        self.a += self.g_adapt * self.spk            # spike-triggered increment
        self.a -= dt * self.a / self.tau_adapt
        inst = self.spk.sum() / self.n / (dt * 1e-3)      # Hz
        self.rate += dt * (inst - self.rate) / 20.0       # 20 ms low-pass
        return self.spk


# --------------------------------------------------------------------------- synapse
E_REV = {"ampa": 0.0, "nmda": 0.0, "gabaa": -75.0, "gly": -80.0}
TAU   = {"ampa": 5.0, "nmda": 100.0, "gabaa": 12.0, "gly": 8.0}


class Syn:
    """Exponential conductance synapse. One instance per (target, receptor)."""
    def __init__(self, n, kind, drug: Drug, tau_scale=1.0, sens=1.0):
        """sens = fraction of this receptor population that is drug-modulatable.

        It must scale BOTH the peak-conductance change and the decay prolongation. Scaling
        only the conductance leaves tau fully prolonged, which silently floors the
        achievable sensitivity (observed: ventilation stuck at -25% however low sens went).
        """
        self.kind, self.n = kind, n
        self.g = np.zeros(n)
        self.tau = TAU[kind] * tau_scale
        eff = lambda x: 1.0 + sens * (x - 1.0)
        if kind == "gabaa":
            self.tau *= eff(drug.gaba_a_tau)
        self.w_scale = {"gabaa": eff(drug.gaba_scale()),
                        "gly":   eff(drug.glyr_gain),
                        "nmda":  1.0 - sens * (1.0 - drug.nmda_scale()),
                        "ampa":  1.0}[kind]

    def inject(self, pre_spikes, w):
        """w is nS of peak conductance per presynaptic spike (a lumped population
        weight, not a unitary synapse)."""
        k = int(pre_spikes.sum())
        if k:
            self.g += w * self.w_scale * k

    def decay(self, dt):
        self.g *= np.exp(-dt / self.tau)

    def conductance(self, V):
        if self.kind == "nmda":
            # Mg2+ block: voltage dependence is why NMDA block is state-dependent
            return self.g / (1.0 + 0.28 * np.exp(-0.062 * V))
        return self.g
def spindle_ia(length, velocity, l0=0.0, kl=55.0, kv=70.0, bias=10.0):
    """Prochazka-style static Ia rate law (Hz). Velocity exponent 0.6 captures the
    strong, saturating velocity sensitivity of primary spindle afferents.
    NOT a Mileusnic/Loeb spindle: no fusimotor drive, no sub-fibre dynamics."""
    v = np.sign(velocity) * (abs(velocity) ** 0.6)
    return float(max(0.0, bias + kl * (length - l0) + kv * v))


def burst_metrics(t, rate, thresh_frac=0.35):
    """Period, duty cycle and burst onsets from a low-passed rate trace."""
    r = np.asarray(rate); t = np.asarray(t)
    if r.max() <= 0: return dict(period=np.nan, duty=np.nan, n_bursts=0, onsets=[])
    hi = r > thresh_frac * r.max()
    edges = np.diff(hi.astype(int))
    on = t[1:][edges == 1]; off = t[1:][edges == -1]
    period = float(np.mean(np.diff(on))) if len(on) > 2 else np.nan
    duty = np.nan
    if len(on) > 1 and len(off) > 1:
        n = min(len(on), len(off))
        d = [off[i] - on[i] for i in range(n) if off[i] > on[i]]
        if d and period == period: duty = float(np.mean(d) / period)
    return dict(period=period, duty=duty, n_bursts=int(len(on)), onsets=on.tolist())
