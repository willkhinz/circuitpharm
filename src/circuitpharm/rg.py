"""Rate-based rhythm generator (Matsuoka half-centre) as a robust locomotor clock.

WHY THIS EXISTS. The spiking LIF rhythm generator in cpg.py oscillates at the
mutual-inhibition LOOP timescale (~60-120 ms), not the adaptation timescale, because a
LIF neuron has no plateau potential: V resets to Vreset on every spike, so the
regenerative I_NaP "up state" that latches a real half-centre ON cannot form. Without
that hysteresis the network cannot be bistable on the fast timescale, so adaptation never
gets to set the period, and no amount of parameter search fixes it -- it is structural.

This module replaces the RG with a Matsuoka oscillator, which oscillates robustly over a
wide parameter range with a period that is analytically controllable:
    period ~ 2*pi*sqrt(tau * tau_a)
so tau=50 ms, tau_a=400 ms gives ~890 ms -- squarely in the rat locomotor range
(300-1500 ms).

This is a deliberate abstraction and it follows the project's own principle: keep
biophysical detail WHERE THE DRUG ACTS, abstract elsewhere. The RG is a clock. The
pharmacologically load-bearing nodes are the tonic extrasynaptic GABA-A conductance on
PF/Mn, the glycinergic reciprocal (IaIn) and recurrent (Renshaw) inhibition, and NMDA on
the Ia->Mn monosynaptic pathway -- all of which stay conductance-based and spiking.

The drug still reaches the RG in a meaningful way: reciprocal inhibition is glycinergic
(scaled by glyr_gain) and the tonic term is GABA-A (scaled by the PAM, with its ceiling).
"""
import numpy as np
from .cpg import Drug


class MatsuokaRG:
    def __init__(self, drug: Drug | None = None, drive=1.0, tau=50.0, tau_a=400.0,
                 beta=2.5, w=2.2, gaba_tonic=0.25, seed=0):
        self.drug = drug or Drug()
        self.tau, self.tau_a, self.beta = tau, tau_a, beta
        self.w = w * self.drug.glyr_gain        # reciprocal inhibition is glycinergic
        self.drive = drive
        # tonic GABA-A appears as a subtractive bias; a PAM deepens it (with its ceiling)
        self.tonic = gaba_tonic * self.drug.gaba_scale()
        rng = np.random.default_rng(seed)
        self.x = rng.normal(0, 0.1, 2)          # [flexor, extensor]
        self.v = np.zeros(2)

    @property
    def y(self):
        return np.maximum(0.0, self.x)

    def step(self, dt):
        y = self.y
        cross = y[::-1]                          # contralateral output
        dx = (-self.x - self.beta * self.v - self.w * cross
              + self.drive - self.tonic) / self.tau
        dv = (-self.v + y) / self.tau_a
        self.x = self.x + dt * dx
        self.v = self.v + dt * dv
        return self.y                            # (flexor, extensor) in [0, inf)
