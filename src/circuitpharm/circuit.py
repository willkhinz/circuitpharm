"""Spinal reflex arc and pattern-formation stages, driven by a spiking rhythm generator.

    GroupPacemakerRG  ->  PF_F/PF_E (spiking)  ->  Mn_F/Mn_E (spiking)
     (circuitpharm/rg2)    |  InPF (gly)            |  Rc (gly, recurrent)
                           +--------------+         |  IaIn (gly, reciprocal)
                                          Ia afferent -> Mn (AMPA + NMDA)

HISTORY, because it explains a retraction. This circuit originally used a rate-based
Matsuoka oscillator, chosen after an exhaustive 288-point search showed an all-spiking LIF
half-centre could not reach a physiological locomotor rhythm (best score 0.83; near-target
period only without alternation, good alternation only at ~90 ms, duty never above 0.24).
That reasoning was sound about LIF half-centres -- a LIF has no plateau potential, so the
half-centre is never bistable on the fast timescale -- but the conclusion was wrong,
because it assumed the rhythm must COME FROM mutual inhibition.

A group pacemaker (recurrent excitation + spike-triggered adaptation) needs no plateau and
so works in LIF, which is why circuitpharm/rg2.py replaced the oscillator outright. The
difference is not cosmetic: in Matsuoka, inhibition IS the oscillator, so glycine block
stops the rhythm. The lamprey strychnine phenotype shows the opposite -- block abolishes
alternation while the rhythm persists. A prediction of this project rested on the Matsuoka
pathway and was retracted. The replacement reproduces the phenotype
(tests/test_phenotypes.py).

All pharmacologically load-bearing nodes are conductance-based:
  tonic extrasynaptic GABA-A on PF and Mn      <- GABA-A PAM, via gaba_scale_tonic
  phasic synaptic GABA-A                       <- GABA-A PAM, via gaba_scale (much weaker;
                                                  the cleft transient is near-saturating)
  glycinergic InPF / IaIn / Renshaw inhibition <- glycine potentiation
  NMDA on RG->PF, PF->Mn and Ia->Mn            <- NMDA antagonist (subunit-resolved)
"""
import numpy as np
from .cpg import Pop, Syn, Drug, E_REV, spindle_ia
from .rg2 import GroupPacemakerRG

# RG OUTPUT UNITS — the thing that makes swapping rhythm generators dangerous.
#
# This circuit converts the rhythm generator's output into an equivalent TOTAL presynaptic
# firing rate before turning it into conductance. How to do that depends on what the RG
# returns, and getting it wrong is recurring error E2 (omitting the presynaptic population
# size, which leaves the downstream pool silent and every weight sweep flat):
#
#   GroupPacemakerRG  returns a PER-NEURON population rate in Hz, so the conversion factor
#                     is the presynaptic POPULATION SIZE (n_exc = 30). Physical units.
#   MatsuokaRG        (removed) returned a dimensionless activation in ~[0,1], so its
#                     factor (170) was a lumped fudge absorbing both population size and an
#                     unknown rate scale.
#
# The Matsuoka oscillator was removed outright rather than kept as an option: it conflates
# rhythm generation with phase setting, so in it glycinergic inhibition IS the oscillator.
# That contradicts the lamprey strychnine phenotype (glycine block abolishes alternation
# while the rhythm PERSISTS) and it invalidated an earlier prediction of this project. Its
# replacement reproduces that phenotype and is covered by tests/test_phenotypes.py.
RG_GAIN_PER_NEURON = 30.0      # = GroupPacemakerRG n_exc


class SpinalCircuit:
    POPS = dict(PF=20, InPF=10, Mn=30, IaIn=10, Rc=10)

    W = {  # nS per presynaptic spike
        ("RG", "PF", "ampa"): 2.2, ("RG", "PF", "nmda"): 1.1,
        ("PF", "InPF", "ampa"): 1.7, ("InPF", "PF*", "gly"): 3.5,
        ("PF", "Mn", "ampa"): 2.2, ("PF", "Mn", "nmda"): 1.1,
        ("Mn", "Rc", "ampa"): 1.3, ("Rc", "Mn", "gly"): 0.45,
        ("IaIn", "Mn*", "gly"): 0.6,
        ("Ia", "Mn", "ampa"): 0.55, ("Ia", "Mn", "nmda"): 0.30,
        ("Ia", "IaIn", "ampa"): 0.45,
    }

    def __init__(self, drug: Drug | None = None, rg_gain=RG_GAIN_PER_NEURON,
                 gaba_tonic=2.0, ia_n=25,
                 gaba_sens=1.0, gaba_sens_tonic=None, gaba_sens_phasic=None,
                 glyr_sens=1.0,
                 drive_pf=110.0, drive_mn=100.0, drive_in=75.0,
                 rg_kw=None, w=None, seed=0):
        self.drug = drug or Drug()
        # separate tonic/phasic sensitivities; both default to gaba_sens, which reproduces
        # the legacy single-pool behaviour exactly
        self.gaba_sens_tonic = gaba_sens if gaba_sens_tonic is None else gaba_sens_tonic
        self.gaba_sens_phasic = gaba_sens if gaba_sens_phasic is None else gaba_sens_phasic
        # Glycine gets its own sensitivity -- it contains no GABA-A subunits, so a
        # GABA-A-derived fraction says nothing about it. See resp.py for the full note.
        self.glyr_sens = glyr_sens
        # FORWARD THE SENSITIVITIES. These were not passed at all, so the rhythm
        # generator always ran at gaba_sens=1.0 / glyr_sens=1.0 -- fully drug-sensitive --
        # while the pattern-formation and motoneuron layers used the calibrated spinal
        # values. Under sedation the RG therefore saw ~7.2x tonic conductance where the
        # rest of the circuit saw ~1.5x, so the locomotor rhythm slowed or arrested far
        # earlier than the tissue pharmacology implies, and every drug effect on step
        # period and coordination was overstated.
        #
        # `rg_kw` still wins, so a caller can override deliberately.
        rg_defaults = dict(gaba_sens_tonic=self.gaba_sens_tonic,
                           gaba_sens_phasic=self.gaba_sens_phasic,
                           glyr_sens=self.glyr_sens)
        rg_defaults.update(rg_kw or {})
        self.rg = GroupPacemakerRG(drug=self.drug, seed=seed, **rg_defaults)
        # rg_gain maps RG output (arb units) to an equivalent TOTAL presynaptic rate,
        # i.e. it already absorbs the presynaptic population size. ia_n is kept explicit
        # because a muscle has ~25 Ia afferents and omitting that factor silently makes
        # the stretch reflex ~25x too weak to reach motoneuron threshold.
        # PER-INSTANCE weights. `W` is still declared at class level as the documented
        # default, but each instance now gets its own copy and an optional override --
        # matching PreBotC and GroupPacemakerRG, which always worked this way.
        #
        # Why this changed: the stretch-reflex assay needs to rescale the Ia->Mn weight
        # (recurring error E5 -- at full weight the motoneuron peak pins at the tref
        # ceiling and every drug reads ~100% of control). With no instance copy, the only
        # way to do that was to MUTATE THE CLASS ATTRIBUTE and restore it in a `finally`.
        # That corrupts weights for every other instance in the process if two assays run
        # concurrently, or if an exception lands between the assignment and the restore.
        self.W = dict(self.W)
        if w:
            self.W.update(w)
        self.rg_gain = rg_gain
        self.ia_n = ia_n
        self.gaba_sens = gaba_sens
        self.gaba_tonic = gaba_tonic
        self.drive = dict(PF=drive_pf, Mn=drive_mn, InPF=drive_in,
                          IaIn=drive_in, Rc=drive_in)
        self.rng = np.random.default_rng(seed)
        self.pops, self.syn = {}, {}
        for half in ("F", "E"):
            for base, n in self.POPS.items():
                nm = f"{base}_{half}"
                self.pops[nm] = Pop(n=n, name=nm + str(seed))
                # physiological firing ceilings: rat motoneurons saturate near
                # 50-100 Hz, so tref=2 ms (=> 500 Hz cap) is far too permissive
                self.pops[nm].tref = 8.0 if base == "Mn" else 4.0
                for rec in ("ampa", "nmda", "gabaa", "gly"):
                    self.syn[(nm, rec)] = Syn(
                        n, rec, self.drug,
                        sens=(self.gaba_sens_phasic if rec == "gabaa"
                              else self.glyr_sens if rec == "gly" else 1.0))
        self.t = 0.0
        self.trace = {k: [] for k in self.pops}
        self.trace.update(t=[], RG_F=[], RG_E=[])

    def _tgt(self, post, half):
        other = "E" if half == "F" else "F"
        return f"{post[:-1]}_{other}" if post.endswith("*") else f"{post}_{half}"

    def _rate_inject(self, target, rec, w, rate_hz, dt):
        k = max(0.0, rate_hz) * dt * 1e-3
        if k > 0:
            s = self.syn[(target, rec)]
            s.g += w * k * (s.w_scale if rec in ("nmda", "gly", "gabaa") else 1.0)

    def step(self, dt, ia_F=0.0, ia_E=0.0):
        y = self.rg.step(dt)                       # (flexor, extensor), arb units
        self.rg_out = y
        # RG -> PF as an equivalent presynaptic rate
        for h, yi in (("F", y[0]), ("E", y[1])):
            for rec in ("ampa", "nmda"):
                self._rate_inject(f"PF_{h}", rec, self.W[("RG","PF",rec)],
                                  yi * self.rg_gain, dt)
        # spike-driven connections
        for (pre, post, rec), w in self.W.items():
            if pre in ("RG", "Ia"):
                continue
            for h in ("F", "E"):
                self.syn[(self._tgt(post, h), rec)].inject(self.pops[f"{pre}_{h}"].spk, w)
        # Ia afferents -- rate is per-afferent, so multiply by the afferent count
        for h, ia in (("F", ia_F), ("E", ia_E)):
            total = np.clip(ia, 0, 200) * self.ia_n
            for post, rec in (("Mn","ampa"), ("Mn","nmda"), ("IaIn","ampa")):
                self._rate_inject(f"{post}_{h}", rec, self.W[("Ia",post,rec)], total, dt)
        # integrate
        for nm, p in self.pops.items():
            base = nm.split("_")[0]
            g, E = {}, {}
            for rec in ("ampa", "nmda", "gabaa", "gly"):
                s = self.syn[(nm, rec)]
                g[rec] = s.conductance(p.V); E[rec] = E_REV[rec]
            if base in ("PF", "Mn"):               # tonic extrasynaptic GABA-A
                # Clamped at zero; see the note in resp.py. A negative tonic term can drive
                # the total gabaa conductance negative, which inverts the inhibitory shunt
                # into regenerative negative damping rather than failing visibly.
                eff = max(0.0, 1.0 + self.gaba_sens_tonic
                          * (self.drug.gaba_scale_tonic() - 1.0))
                g["gabaa"] = np.maximum(0.0, g["gabaa"] + self.gaba_tonic * eff)
            p.step(dt, g, E, self.drive.get(base, 0.0), self.rng)
        for s in self.syn.values():
            s.decay(dt)
        self.t += dt
        return {nm: p.rate for nm, p in self.pops.items()}

    def record(self):
        self.trace["t"].append(self.t)
        self.trace["RG_F"].append(float(self.rg_out[0]))
        self.trace["RG_E"].append(float(self.rg_out[1]))
        for nm, p in self.pops.items():
            self.trace[nm].append(p.rate)

    def arrays(self):
        return {k: np.asarray(v) for k, v in self.trace.items()}
