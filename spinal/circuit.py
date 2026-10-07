"""Hybrid spinal circuit: rate-based rhythm generator + spiking output stages.

    MatsuokaRG (clock)  ->  PF_F/PF_E (spiking)  ->  Mn_F/Mn_E (spiking)
                             |  InPF (gly)          |  Rc (gly, recurrent)
                             +--------------+       |  IaIn (gly, reciprocal)
                                            Ia afferent -> Mn (AMPA + NMDA)

Rationale for the split is in spinal/rg.py: an exhaustive 288-point search showed the
all-spiking LIF rhythm generator cannot reach a physiological locomotor rhythm (best
score 0.83; near-target period only with no alternation, good alternation only at
~90 ms, duty never above 0.24). The cause is structural, not parametric -- a LIF has no
plateau potential, so the half-centre is never bistable on the fast timescale and
adaptation cannot set the period.

All pharmacologically load-bearing nodes remain conductance-based:
  tonic extrasynaptic GABA-A on PF and Mn      <- GABA-A PAM (with efficacy ceiling)
  glycinergic InPF / IaIn / Renshaw inhibition <- glycine potentiation
  NMDA on RG->PF, PF->Mn and Ia->Mn            <- NMDA antagonist (subunit-resolved)
  RG reciprocal inhibition + tonic bias        <- glycine / GABA-A (inside MatsuokaRG)
"""
import numpy as np
from .cpg import Pop, Syn, Drug, E_REV, spindle_ia
from .rg import MatsuokaRG


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

    def __init__(self, drug: Drug | None = None, rg_gain=170.0, gaba_tonic=2.0, ia_n=25,
                 gaba_sens=1.0,
                 drive_pf=110.0, drive_mn=100.0, drive_in=75.0,
                 rg_kw=None, seed=0):
        self.drug = drug or Drug()
        self.rg = MatsuokaRG(drug=self.drug, seed=seed, **(rg_kw or {}))
        # rg_gain maps RG output (arb units) to an equivalent TOTAL presynaptic rate,
        # i.e. it already absorbs the presynaptic population size. ia_n is kept explicit
        # because a muscle has ~25 Ia afferents and omitting that factor silently makes
        # the stretch reflex ~25x too weak to reach motoneuron threshold.
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
                        sens=(gaba_sens if rec in ("gabaa", "gly") else 1.0))
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
                eff = 1.0 + self.gaba_sens * (self.drug.gaba_scale_tonic() - 1.0)
                g["gabaa"] = g["gabaa"] + self.gaba_tonic * eff
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
