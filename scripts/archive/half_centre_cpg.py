"""The REJECTED all-spiking locomotor CPG. Archived; not part of the package.

WHY THIS IS HERE AND NOT IN circuitpharm. `HalfCentreCPG` generates its rhythm from mutual
inhibition between flexor and extensor half-centres. An exhaustive 288-point search showed
it cannot reach a physiological locomotor rhythm in LIF neurons (best score 0.83:
near-target period only without alternation, good alternation only at ~90 ms, duty never
above 0.24). The cause is structural -- a LIF resets V at every spike, so no plateau
potential forms, the half-centre is never bistable on the fast timescale, and adaptation
cannot set the period.

That reasoning was right about LIF HALF-CENTRES and wrong about the conclusion, which
assumed the rhythm must come from inhibition. A GROUP PACEMAKER (recurrent excitation plus
spike-triggered adaptation) needs no plateau and works fine in LIF -- see
circuitpharm/rg2.py, which also reproduces the lamprey strychnine phenotype that this
whole class of architecture gets backwards.

It was moved out of the package in session 8 because it was 103 lines of dead code (28% of
cpg.py) reachable only from the three archived scripts beside this file.
"""
import numpy as np

from circuitpharm.cpg import Pop, Syn, Drug, E_REV




# --------------------------------------------------------------------------- network
class HalfCentreCPG:
    """One antagonist joint pair. Flexor/extensor half-centres with glycinergic
    reciprocal inhibition, Renshaw recurrent inhibition, Ia afferent feedback, and a
    tonic GABA-A background. Step with dt in ms."""

    POPS = dict(RG=(20, True), InRG=(10, False), PF=(20, False), InPF=(10, False),
                Mn=(30, False), IaIn=(10, False), Rc=(10, False))

    # (pre, post, receptor, weight) -- 'F'/'E' suffixes added per limb side
    WIRING = [
        ("RG",   "RG",   "ampa",  0.25),     # recurrent excitation within half-centre
        ("RG",   "InRG", "ampa",  1.70),
        ("InRG", "RG*",  "gly",   6.00),     # * = contralateral half-centre
        ("RG",   "PF",   "ampa",  2.20),
        ("RG",   "PF",   "nmda",  1.10),
        ("PF",   "InPF", "ampa",  1.70),
        ("InPF", "PF*",  "gly",   3.50),
        ("PF",   "Mn",   "ampa",  2.20),
        ("PF",   "Mn",   "nmda",  1.10),
        ("Mn",   "Rc",   "ampa",  1.30),
        ("Rc",   "Mn",   "gly",   0.45),     # recurrent
        ("IaIn", "Mn*",  "gly",   0.60),     # reciprocal Ia inhibition
    ]

    def __init__(self, drug: Drug | None = None, drive=135.0, gaba_tonic=2.0,
                 drive_pf=110.0, drive_mn=100.0, drive_in=75.0,
                 tau_adapt=800.0, seed=0):
        self.drug = drug or Drug()
        self.drive = drive            # pA tonic drive to RG. MUST be sub-rheobase (~200 pA):
        # I_NaP carries the cell over threshold, so mutual inhibition + adaptation can
        # switch a half-centre off. Supra-rheobase drive makes it un-switchable.
        self.gaba_tonic = gaba_tonic  # nS standing extrasynaptic GABA-A conductance
        # rheobase with tonic GABA: I = (gL+gt)*Vth - (gL*EL + gt*Egaba) ~= 200 pA
        # PF and Mn carry their own sub-rheobase tonic drive (descending +
        # propriospinal), so CPG input only has to push them the rest of the way.
        # Without it the RG->PF->Mn chain needs implausibly large synaptic weights.
        self.drive_pf, self.drive_mn = drive_pf, drive_mn
        # Interneurons need tonic drive too. Without it InRG sits at ~-52 mV against a
        # -50 mV threshold, barely fires, and reciprocal inhibition never engages --
        # which makes the inhibitory weight irrelevant and silently breaks the
        # half-centre mechanism.
        self.drive_in = drive_in
        self.rng = np.random.default_rng(seed)
        self.pops, self.syn = {}, {}
        for half in ("F", "E"):
            for base, (n, nap) in self.POPS.items():
                nm = f"{base}_{half}"
                self.pops[nm] = Pop(n=n, name=nm + str(seed), nap=nap)
                self.pops[nm].tau_adapt = tau_adapt
                for rec in ("ampa", "nmda", "gabaa", "gly"):
                    self.syn[(nm, rec)] = Syn(n, rec, self.drug)
        self.t = 0.0
        self.trace = {k: [] for k in self.pops}
        self.trace["t"] = []

    def _post(self, base, half):
        other = "E" if half == "F" else "F"
        return f"{base[:-1]}_{other}" if base.endswith("*") else f"{base}_{half}"

    def step(self, dt, ia_F=0.0, ia_E=0.0, drive_scale=1.0):
        """ia_*: Ia afferent firing rate (Hz) from the spindle model."""
        # 1. route spikes from the previous step
        for pre, post, rec, w in self.WIRING:
            for half in ("F", "E"):
                p_pre = self.pops[f"{pre}_{half}"]
                self.syn[(self._post(post, half), rec)].inject(p_pre.spk, w)
        # 2. Ia afferents: monosynaptic AMPA+NMDA to Mn, AMPA to IaIn
        for half, ia in (("F", ia_F), ("E", ia_E)):
            k = np.clip(ia, 0, 200) * dt * 1e-3            # expected spikes
            if k > 0:
                self.syn[(f"Mn_{half}",   "ampa")].g += 0.55 * k
                s = self.syn[(f"Mn_{half}", "nmda")]
                s.g += 0.30 * k * s.w_scale
                self.syn[(f"IaIn_{half}", "ampa")].g += 0.45 * k
        # 3. integrate (tonic GABA-A is added as a constant conductance below)
        for nm, p in self.pops.items():
            g, E = {}, {}
            for rec in ("ampa", "nmda", "gabaa", "gly"):
                s = self.syn[(nm, rec)]
                g[rec] = s.conductance(p.V); E[rec] = E_REV[rec]
            # TONIC extrasynaptic GABA-A conductance on RG/PF/Mn. This is a constant,
            # not a Poisson train -- tonic inhibition is literally a standing
            # conductance. It is also the single most important drug node: ethanol's
            # putative GABA-A action is on extrasynaptic delta-containing receptors
            # mediating exactly this current, and a PAM scales it directly.
            if nm.split("_")[0] in ("RG", "PF", "Mn"):
                g["gabaa"] = g["gabaa"] + self.gaba_tonic * self.drug.gaba_scale()
            base = nm.split("_")[0]
            Id = {"RG": self.drive * drive_scale, "PF": self.drive_pf,
                  "Mn": self.drive_mn}.get(base, self.drive_in)
            p.step(dt, g, E, Id, self.rng)
        for s in self.syn.values():
            s.decay(dt)
        self.t += dt
        return {nm: p.rate for nm, p in self.pops.items()}

    def record(self):
        self.trace["t"].append(self.t)
        for nm, p in self.pops.items():
            self.trace[nm].append(p.rate)

    def arrays(self):
        return {k: np.asarray(v) for k, v in self.trace.items()}
