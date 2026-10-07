"""Locomotor rhythm generator, architecture v2 — group pacemaker with phase coupling.

SUPERSEDES circuitpharm/rg.py (Matsuoka). Reason for the rewrite, from literature:

    In lamprey fictive locomotion, strychnine (glycine receptor BLOCK) ELIMINATES
    left-right alternation while robust rhythmic activity PERSISTS -- same burst
    proportion, same rostro-caudal coordination, just co-activation instead of
    alternation. Strychnine also lengthened cycle period.

So in the real cord, glycinergic reciprocal inhibition sets the PHASE RELATIONSHIP; it
does NOT generate the rhythm. The Matsuoka oscillator conflates the two: there, mutual
inhibition IS the oscillator, so removing or over-strengthening it stops the rhythm. That
made the earlier prediction (glycine potentiation slows the cycle) rest on a pathway that
does not work this way in biology. Retracted.

v2 architecture -- each half-centre is an INDEPENDENT group pacemaker:

    RG_F: recurrent excitatory population + spike-triggered adaptation  -> bursts alone
    RG_E: same
    InRG_F/E: inhibitory interneurons, glycinergic, cross-projecting    -> sets anti-phase

Rhythm comes from recurrent excitation + adaptation (the preBotC mechanism already
validated in circuitpharm/resp.py, and it works in LIF because it needs no plateau potential).
Inhibition only enforces alternation. Therefore blocking glycine must leave the rhythm
intact and only destroy the phase relationship -- which is the published, falsifiable
validation target for this module.
"""
import numpy as np
from .cpg import Pop, Syn, Drug, E_REV


class GroupPacemakerRG:
    def __init__(self, drug: Drug | None = None, n_exc=30, n_inh=12,
                 drive=150.0, gaba_tonic=1.0, g_adapt=2.2, tau_adapt=500.0,
                 gaba_sens=1.0, asym=0.08, w=None, seed=0):
        self.drug = drug or Drug()
        self.W = dict(ee_ampa=0.42, ee_nmda=0.23,   # recurrent excitation (rhythmogenic)
                      ei_ampa=0.80,                 # half-centre -> own interneurons
                      ie_gly=1.20)                  # interneurons -> OTHER half-centre
        if w:
            self.W.update(w)
        # Deliberate intrinsic-frequency MISMATCH between the half-centres.
        # Without it the coupling test is meaningless: two independent oscillators with
        # identical intrinsic period hold whatever arbitrary phase offset they start in,
        # so the flexor/extensor correlation is the same coupled or uncoupled (observed:
        # -0.55 both ways). With a mismatch, coupled oscillators ENTRAIN to a common
        # frequency in anti-phase while uncoupled ones drift, so correlation -> 0.
        # That is what makes 'blocking glycine abolishes alternation' testable at all.
        self.drive_half = {"F": drive * (1.0 + asym), "E": drive * (1.0 - asym)}
        self.drive, self.gaba_tonic = drive, gaba_tonic
        self.gaba_sens = gaba_sens
        self.rng = np.random.default_rng(seed)
        self.pops, self.syn = {}, {}
        for half in ("F", "E"):
            for base, n, adapt in (("RG", n_exc, g_adapt), ("InRG", n_inh, 0.0)):
                nm = f"{base}_{half}"
                p = Pop(n=n, name=nm + str(seed))
                p.g_adapt, p.tau_adapt, p.tref = adapt, tau_adapt, 5.0
                self.pops[nm] = p
                for rec in ("ampa", "nmda", "gabaa", "gly"):
                    self.syn[(nm, rec)] = Syn(
                        n, rec, self.drug,
                        sens=(gaba_sens if rec in ("gabaa", "gly") else 1.0))
        self.t = 0.0
        self.trace = {k: [] for k in list(self.pops) + ["t"]}

    def step(self, dt, drive_scale=1.0):
        other = {"F": "E", "E": "F"}
        for h in ("F", "E"):
            srg = self.pops[f"RG_{h}"].spk
            sin = self.pops[f"InRG_{h}"].spk
            self.syn[(f"RG_{h}", "ampa")].inject(srg, self.W["ee_ampa"])
            self.syn[(f"RG_{h}", "nmda")].inject(srg, self.W["ee_nmda"])
            self.syn[(f"InRG_{h}", "ampa")].inject(srg, self.W["ei_ampa"])
            # glycinergic projection to the CONTRALATERAL half-centre: phase only
            self.syn[(f"RG_{other[h]}", "gly")].inject(sin, self.W["ie_gly"])
        for nm, p in self.pops.items():
            g, E = {}, {}
            for rec in ("ampa", "nmda", "gabaa", "gly"):
                s = self.syn[(nm, rec)]
                g[rec] = s.conductance(p.V); E[rec] = E_REV[rec]
            eff = 1.0 + self.gaba_sens * (self.drug.gaba_scale() - 1.0)
            g["gabaa"] = g["gabaa"] + self.gaba_tonic * eff
            Id = (self.drive_half[nm[-1]] * drive_scale if nm.startswith("RG") else 70.0)
            p.step(dt, g, E, Id, self.rng)
        for s in self.syn.values():
            s.decay(dt)
        self.t += dt
        return self.pops["RG_F"].rate, self.pops["RG_E"].rate

    def record(self):
        self.trace["t"].append(self.t)
        for nm, p in self.pops.items():
            self.trace[nm].append(p.rate)

    def arrays(self):
        return {k: np.asarray(v) for k, v in self.trace.items()}
