"""preBotzinger complex: respiratory rhythm generator and the respiratory-depression axis.

ARCHITECTURE CHOICE. Unlike the locomotor half-centre (which needs plateau potentials
that a LIF cannot produce -- see circuitpharm/rg.py), the preBotC rhythm is modelled as a
"group pacemaker": a RECURRENTLY EXCITATORY glutamatergic population whose synchronous
burst is terminated by spike-triggered adaptation. Recurrent excitation gives
regenerative recruitment and synchrony, so this mechanism DOES work in a LIF network --
no plateau required. That is why this module is all-spiking while the locomotor RG is not.

    Exc (glutamatergic, recurrent AMPA+NMDA, adapting)  -->  Out (inspiratory motor)
      |                                   ^
      +--> Inh (GABA-A / glycinergic) -----+

Pharmacology, and why the margin lives here:
  NMDA block    -- recurrent excitation is the rhythmogenic core, so NMDA antagonism
                   attacks rhythm generation DIRECTLY. The literature is explicit that
                   glutamatergic neurons are the critical component for both
                   rhythmogenesis and opioid-induced respiratory depression.
  GABA-A PAM    -- tonic + phasic inhibition. The respiratory network expresses alpha1,
                   alpha2/3, alpha4/delta and epsilon subunits, so there is NO subtype
                   window here; safety must come from the PAM efficacy ceiling.
  Glycine       -- inhibitory cotransmission in the respiratory core.

THE KEY ASYMMETRY the whole project rests on: brainstem is GluN2D-dominant (low GluN2B
fraction), so a GluN2B-selective antagonist should barely touch this network while
producing a substantial forebrain effect. `glun2b_fraction` is the knob that encodes it.
"""
import numpy as np
from . import substrate as sub
from .cpg import Pop, Syn, Drug, E_REV

# Physiologically admissible respiratory frequency, Hz. Rat eupnoea is 1-2 Hz (60-120
# breaths/min). The band is deliberately generous at both ends: a respiratory depressant
# SLOWS breathing, so the lower bound only has to exclude near-apnoea, and the upper bound
# only has to exclude spectral peaks that cannot be breathing at all. Used as a hard gate
# on `alive` -- see E4 SECOND FORM in resp_metrics.
EUPNOEA_BAND = (0.30, 2.50)

# Every weight the network needs, so a conductance operating point can be checked for
# completeness rather than silently half-inherited from the LIF (see __init__).
REQUIRED_W = ("ee_ampa", "ee_nmda", "ei_ampa", "ei_nmda",
              "ie_gaba", "ie_gly", "eo_ampa", "eo_nmda")


class PreBotC:
    def __init__(self, drug: Drug | None = None, n_exc=50, n_inh=20, n_out=20,
                 drive=None, gaba_tonic=1.5, g_adapt=None, tau_adapt=None,
                 gaba_sens=1.0, gaba_sens_tonic=None, gaba_sens_phasic=None,
                 glyr_sens=1.0,
                 resp_drive_boost=0.0, w=None, seed=0,
                 substrate="lif", cell=None, op=None):
        self.drug = drug or Drug()
        # SUBSTRATE SELECTION IS EXPLICIT AND DEFAULTS TO THE VALIDATED PATH.
        #
        # 'lif'  -- cpg.Pop, byte-identical to every result this project has published.
        # 'cond' -- neuron.CondPop, a real I_NaP (roadmap link 4).
        #
        # A conductance operating point is NOT the LIF one. The LIF cell is C=200 pF,
        # g_L=10 nS; the Butera cell is C=21 pF, g_L=2.8 nS. `drive` is in pA and every
        # synaptic weight is in nS, so running the LIF's RESP_OP on the conductance cell is
        # wrong by roughly an order of magnitude -- and it would not raise, it would just
        # produce a network that does not oscillate or oscillates for the wrong reason
        # (predicted as E14). So `substrate='cond'` REFUSES to run without an operating
        # point, rather than silently inheriting numbers tuned against a different cell.
        from .config import COND_RESP_OP
        _op = sub.resolve_op(substrate, op, COND_RESP_OP, required_w=REQUIRED_W,
                             required_scalars=("drive", "drive_other", "gaba_tonic"),
                             what="the respiratory rhythm generator")
        if _op is not None:
            drive = _op["drive"]
            gaba_tonic = _op["gaba_tonic"]
            self.drive_other = _op["drive_other"]
            w = dict(_op["w"], **(w or {}))
        # DEFAULTS COME FROM config.RESP_OP, not from literals here. They used to be
        # drive=190, g_adapt=1.6, tau_adapt=450, ee_ampa=0.16, ee_nmda=0.09 -- an obsolete
        # untuned set with ~3x weaker recurrent excitation than the calibrated operating
        # point. Every analysis passes **RESP_OP, so this was invisible in project results,
        # but any caller writing a plain PreBotC() silently got the wrong network.
        from .config import RESP_OP as _OP
        if drive is None:
            drive = _OP["drive"]
        if g_adapt is None:
            g_adapt = _OP["g_adapt"]
        if tau_adapt is None:
            tau_adapt = _OP["tau_adapt"]
        self.W = dict(_OP["w"])                       # recurrent excitation (rhythmogenic)
        self.W.update(ei_ampa=0.55, ei_nmda=0.0,      # exc -> inh (see DISINHIBITION below)
                      ie_gaba=0.45, ie_gly=0.35,      # inh -> exc
                      eo_ampa=0.70, eo_nmda=0.30)     # exc -> output (phrenic analogue)
        # DISINHIBITION PATHWAY. ei_nmda defaults to 0.0, which reproduces the model as it
        # stood when the "NMDA block always depresses ventilation" limitation was declared.
        # That declaration was reached by sweeping the share of recurrent excitation
        # (ee_nmda) only -- and with ei_nmda absent, NMDA block CANNOT reduce inhibition by
        # any amount, so the disinhibition arm of ketamine's action was structurally
        # unrepresentable rather than merely mis-parameterised. preBotC inhibitory
        # interneurons do carry NMDA receptors, so a nonzero share is the biologically
        # correct default; it is left at 0 only so earlier runs stay reproducible, and
        # scripts/calib_ei_nmda.py sweeps it.
        if w: self.W.update(w)
        self.drive, self.gaba_tonic = drive, gaba_tonic
        # gaba_sens = FRACTION OF preBotC GABA-A CONDUCTANCE THAT IS DRUG-MODULATABLE.
        # Not a fudge factor: a benzodiazepine-site PAM requires a gamma2 subunit, and the
        # respiratory network also expresses delta, alpha4 and epsilon subunits -- and
        # epsilon-containing receptors (enriched on NK1R+ rhythm-generating neurons) are
        # BZ-insensitive. So only part of the local GABA-A conductance can be modulated.
        # Calibrated against human midazolam 2 mg IV (a clearly sedative dose), which
        # reduces minute ventilation by only ~14-19%. Uncalibrated (gaba_sens=1.0) the
        # model gave ~50% at PAM 1.5x, i.e. 3-4x too sensitive.
        self.gaba_sens = gaba_sens
        # The two pools need SEPARATE sensitivities, because subtype localisation differs
        # (circuitpharm/subtypes.py EXTRASYN) and the two pools are driven by gains that differ
        # ~200-fold in headroom. Both default to `gaba_sens`, which reproduces the old
        # single-pool behaviour exactly.
        self.gaba_sens_tonic = gaba_sens if gaba_sens_tonic is None else gaba_sens_tonic
        self.gaba_sens_phasic = gaba_sens if gaba_sens_phasic is None else gaba_sens_phasic
        # GLYCINE GETS ITS OWN SENSITIVITY. `gly` used to be lumped with `gabaa` and
        # scaled by the GABA-A-derived fraction, which is a category error: glycine
        # receptors contain no GABA-A subunits, so a fraction computed from alpha1/
        # alpha2-3/alpha5 expression says nothing about how much of the glycinergic
        # conductance a drug can reach. The effect was severe and silent -- at alogabat's
        # preBotC phasic sensitivity a 1.6x glycine potentiation became 1.0019x, i.e.
        # ethanol's glycine mechanism and strychnine were ~99.8% deleted from the circuit.
        #
        # Default 1.0: a glycine-site ligand reaches the glycinergic conductance. Narrow it
        # only with a glycine-specific measurement, never with a GABA-A number.
        self.glyr_sens = glyr_sens

        # EMPIRICAL OVERRIDE, not derived from this circuit.
        # Ketamine-class NMDA channel blockers do NOT depress breathing -- clinically they
        # PRESERVE respiratory drive and airway reflexes, and s-ketamine actively
        # STIMULATES breathing and ATTENUATES propofol- and opioid-induced hypoventilation.
        # This model contains only the direct glutamatergic rhythmogenic drive to the
        # preBotC, so it predicts the opposite sign. The real stimulation arises through
        # mechanisms outside the modelled core: block of NMDA on inhibitory interneurons
        # (disinhibition), raised sympathetic tone, and upper-airway dilator activation.
        # resp_drive_boost adds that known effect as a fractional increase in excitatory
        # drive. It is imposed from clinical data, NOT predicted, and is flagged as such
        # wherever it is used.
        self.resp_drive_boost = resp_drive_boost
        # LIF value, unchanged; the cond path overwrote this from its operating point above.
        self.drive_other = getattr(self, "drive_other", 60.0)
        self.rng = np.random.default_rng(seed)
        self.substrate = substrate
        self.pops = {}
        for nm, n, nap in (("Exc", n_exc, True), ("Inh", n_inh, False), ("Out", n_out, False)):
            p = sub.make_pop(substrate, n, f"{nm}{seed}", nap=nap,
                             g_adapt=(g_adapt if nm == "Exc" else 0.0),
                             tau_adapt=tau_adapt, tref=5.0, cell=cell)
            self.pops[nm] = p
        self.syn = {(nm, rec): Syn(p.n, rec, self.drug,
                                   sens=(self.gaba_sens_phasic if rec == "gabaa"
                                         else self.glyr_sens if rec == "gly"
                                         else 1.0))
                    for nm, p in self.pops.items()
                    for rec in ("ampa", "nmda", "gabaa", "gly")}
        self.t = 0.0
        self.trace = {k: [] for k in ("t", "Exc", "Inh", "Out")}

    def step(self, dt, drive_scale=1.0):
        se, si = self.pops["Exc"].spk, self.pops["Inh"].spk
        self.syn[("Exc","ampa")].inject(se, self.W["ee_ampa"])
        self.syn[("Exc","nmda")].inject(se, self.W["ee_nmda"])
        self.syn[("Inh","ampa")].inject(se, self.W["ei_ampa"])
        if self.W["ei_nmda"]:
            self.syn[("Inh","nmda")].inject(se, self.W["ei_nmda"])
        self.syn[("Out","ampa")].inject(se, self.W["eo_ampa"])
        self.syn[("Out","nmda")].inject(se, self.W["eo_nmda"])
        self.syn[("Exc","gabaa")].inject(si, self.W["ie_gaba"])
        self.syn[("Exc","gly")].inject(si, self.W["ie_gly"])
        # ON THE CONDUCTANCE SUBSTRATE, NMDA IS HANDED OVER UNEVALUATED.
        #
        # `Syn.conductance(V)` applies the Mg2+ block at the PRE-STEP voltage. That is
        # harmless on the LIF, where V moves a few mV per step -- and it would discard the
        # entire voltage-dependent relief on a cell whose V sweeps ~70 mV through a spike,
        # because the cell sub-steps internally at 0.05 ms while the block would stay frozen
        # at the value it had before the upstroke.
        #
        # That is the single most likely way for this upgrade to do nothing for NMDA while
        # appearing to work: the symptom is "the substrate change did not move the NMDA
        # result", which reads as a reassuring robustness check and is in fact the bug.
        # Predicted as E15 before the cell was written. `raw=("nmda",)` makes CondPop apply
        # the block at its own V, every substep.
        raw = sub.raw_receptors(self.substrate)
        for nm, p in self.pops.items():
            g, E = {}, {}
            for rec in ("ampa", "nmda", "gabaa", "gly"):
                s = self.syn[(nm, rec)]
                g[rec] = s.g if rec in raw else s.conductance(p.V)
                E[rec] = E_REV[rec]
            g["gabaa"] = sub.tonic_gaba(self.gaba_tonic, self.gaba_sens_tonic,
                                        self.drug.gaba_scale_tonic(), g["gabaa"])
            # `drive_other` was the literal 60.0, which is pA against the LIF cell -- a
            # third LIF-scaled quantity hiding in the step loop after `drive` and the weight
            # table had both been moved into the operating point. On the conductance cell
            # (C=21 pF against 200 pF) it is the Inh and Out populations' entire excitatory
            # drive, so leaving it behind would have left two of three populations driven at
            # roughly the wrong order of magnitude while the rhythm still ran.
            Id = (self.drive * drive_scale * (1.0 + self.resp_drive_boost)
                  if nm == "Exc" else self.drive_other)
            if raw:
                p.step(dt, g, E, Id, self.rng, raw=raw)
            else:
                p.step(dt, g, E, Id, self.rng)
        for s in self.syn.values():
            s.decay(dt)
        self.t += dt

    def record(self):
        self.trace["t"].append(self.t)
        for k in ("Exc", "Inh", "Out"):
            self.trace[k].append(self.pops[k].rate)

    def arrays(self):
        return {k: np.asarray(v) for k, v in self.trace.items()}


def resp_metrics(t, out, exc, ctrl_mean=None, band=None):
    """Robust respiratory metrics. Rat eupnoea is ~1-2 Hz (60-120 breaths/min).

    A RELATIVE threshold (frac * signal.max()) fails catastrophically here: when a drug
    degrades the rhythm, the lowered max puts the threshold into the noise and spurious
    crossings read as a HIGHER frequency -- so respiratory depression reports as
    tachypnoea (observed: 13 Hz = 780 breaths/min). Frequency is therefore taken from the
    FFT peak, which is amplitude-invariant, and burst structure from a modulation index.
    Apnoea is declared on collapse of modulation or of mean output, never on crossings.

    E4, SECOND FORM (found 2026-10-07). The FFT fix above cures *spurious* crossings but
    not genuine high-frequency power. Under heavy NMDA block the burst train FRAGMENTS, and
    the FFT correctly reports a 4.06 Hz peak -- ~240 breaths/min, with mean output at 45% of
    control. That is not tachypnoea, it is a disintegrated rhythm, and `mod > 0.8 and n >= 3`
    passed it as alive. Because `alive` is what the overdose scan keys on, this made the
    overdose index OPTIMISTIC: a pathological rhythm counted as survival.

    Fix: a physiological plausibility band. Rat eupnoea is 1-2 Hz; a respiratory DEPRESSANT
    lowers frequency, so a high peak combined with reduced mean output is incoherent as
    breathing whatever its spectral power. EUPNOEA_BAND is therefore a hard gate on `alive`,
    and `frag` is reported so the caller can see which failure mode occurred. This makes the
    overdose index strictly more conservative, which is the honest direction.
    """
    # THE BAND IS A PROPERTY OF THE PREPARATION, not of the metric. EUPNOEA_BAND is an
    # in vivo rat band and the default, so every existing caller is unchanged. The
    # conductance substrate passes config.INVITRO_BAND, because the Butera cell is neonatal
    # rodent in vitro and its network's physiological frequency is ~0.1-0.3 Hz; see the long
    # note on INVITRO_BAND in config.py for why that is a preparation difference and not a
    # relaxed gate.
    band_lo, band_hi = EUPNOEA_BAND if band is None else band
    t = np.asarray(t, float); o = np.asarray(out, float)
    dt_ms = float(t[1] - t[0]) if len(t) > 1 else 1.0
    w = max(1, int(100.0 / max(1e-9, dt_ms)))          # ~100 ms smoothing
    sm = np.convolve(o, np.ones(w) / w, mode="same")
    if len(sm) > 3 * w:                                 # drop convolution edges
        sm, t = sm[w:-w], t[w:-w]
    mean = float(sm.mean())
    if mean < 1.0:
        return dict(freq=0.0, amp=0.0, mod=0.0, mean=mean, alive=False, n=0,
                    frag=False, reason="mean output collapsed")
    p5, p95 = np.percentile(sm, 5), np.percentile(sm, 95)
    mod = float((p95 - p5) / max(1e-9, mean))
    dtr = (t[1] - t[0]) / 1000.0 if len(t) > 1 else dt_ms / 1000.0
    x = (sm - sm.mean()) * np.hanning(len(sm))
    F = np.abs(np.fft.rfft(x)); fr = np.fft.rfftfreq(len(sm), dtr)
    # The FFT SEARCH range must reach below the validity band's floor, or a rhythm
    # inside the in vitro band (down to 0.05 Hz) cannot be found at all -- the peak would
    # be picked from whatever sits above 0.25 Hz, which is noise or a harmonic.
    # 0.25 Hz for the in vivo default, so the LIF path is untouched; lowered only when the
    # validity band itself reaches below it. My first version wrote min(0.25, 0.5*band_lo),
    # which LOWERS the floor to 0.15 Hz for the in vivo band too -- changing the search
    # range on the substrate whose results must stay byte-identical, for no reason.
    fft_lo = 0.25 if band_lo >= 0.25 else 0.5 * band_lo
    band = (fr > fft_lo) & (fr < 5.0)
    freq = float(fr[band][np.argmax(F[band])]) if band.any() else 0.0
    thr = p5 + 0.5 * (p95 - p5)                         # midpoint of the ACTUAL range
    n = int((np.diff((sm > thr).astype(int)) == 1).sum())
    # physiological plausibility gate -- see E4 SECOND FORM in the docstring
    frag = bool(freq > band_hi)
    ok_mod = mod > 0.8
    ok_n = n >= 3
    ok_mean = ctrl_mean is None or mean > 0.2 * ctrl_mean
    ok_band = band_lo <= freq <= band_hi
    alive = bool(ok_mod and ok_n and ok_mean and ok_band)
    reason = ("alive" if alive else
              "mean below 20% of control" if not ok_mean else
              f"frequency {freq:.2f} Hz outside band "
              f"{band_lo}-{band_hi} Hz "
              f"({'fragmented rhythm' if frag else 'too slow'})" if not ok_band else
              "modulation collapsed" if not ok_mod else "too few bursts")
    return dict(freq=freq, amp=float(p95), mod=mod,
                mean=mean, alive=alive, n=n, frag=frag, reason=reason)
