# Design: conductance-based neuron substrate (roadmap links 4 + 5)

**Status: DESIGN ONLY. Do not implement yet.** The gate in `KNOWLEDGE.md` still stands —
one review pass must find nothing before a new surface this large goes onto this base. This
document exists so the design is not lost the way the roadmap was (recorded only in
conversation) and the way the tuned RG parameters were (recorded only in a past run's
stdout). It is the third time this project has nearly lost work by leaving it unwritten.

Written 2026-10-07. Every diagnostic measurement in §1 is reproducible with
`scripts/diag_substrate_limits.py`, which changes no model state and fits nothing.

---

## 1. The problem, measured rather than asserted

The roadmap's stated reason for links 4–5 was "LIF is not the biophysics." That is true but
weak — it argues from principle, and this project has been wrong from principle before. So I
measured what the LIF substrate actually does to the voltage-dependent mechanisms the
pharmacology depends on. Instrumented `PreBotC(**RESP_OP)`, 15 s, 250k samples of the `Exc`
population after settling:

```
V range              -71.11 .. -44.17 mV        (Vth = -50, so nothing above it but noise)
V  p1/p50/p99        -69.72 / -65.77 / -50.75

NMDA Mg2+ relief     mean 0.0632   range 0.0417 .. 0.1877
  dynamic range      4.50x         (a spiking cell spanning -65..+20 mV: 15.5x)
  relief at +20 mV   0.9250        unreachable on this substrate

GABA-A driving force mean 10.4 mV  max 30.8 mV
  at V = +20 mV      95.0 mV       = 9.11x the LIF mean
glutamate driving f. mean 64.6 mV  min 44.2 mV  (never collapses toward 0)
```

Three findings, in increasing order of how much they matter.

**(a) NMDA conductance is held in near-permanent Mg²⁺ block.** Mean relief 0.063 — the NMDA
conductance delivers ~6% of its unblocked value, always. `cpg.py:243` carries the comment
"voltage dependence is why NMDA block is state-dependent." On this substrate it is not:
available relief range is 4.5× against the 15.5× a spiking cell has, and the upper end
(0.925 at +20 mV) is structurally unreachable.

**(b) Excitation never self-limits.** Glutamate driving force stays in [44, 65] mV and never
approaches zero, because V never approaches E_glut = 0. In a real cell the excitatory current
collapses at the burst peak; that collapse is part of what shapes a burst. Here it is absent.

**(c) The phasic GABA-A pool's driving force is truncated exactly when it matters.** This is
the one that threatens a headline result.

Both GABA-A pools see the same V, so at rest both get the same ~10 mV driving force. The
asymmetry is in *when* each arrives. The tonic conductance is a standing term — always on,
always at the resting driving force, which the LIF represents correctly. The phasic
conductance arrives in events correlated with the burst (`Inh` fires driven by `Exc`), i.e.
precisely when a real cell would be depolarized to spike voltages and the driving force would
be up to 95 mV. On this substrate it arrives at a cell clamped below -44 mV.

**Why this is not absorbable into the hand-tuned weights.** The standing objection to all
three findings is that the weights were tuned on this substrate, so absolute magnitudes were
fitted to compensate. That is correct and it is why (a) and (b) are listed first and
lightly — a scalar weight can rescale a mean. What a scalar weight *cannot* restore is a
voltage-**dependence**: it cannot turn a 4.5× relief range into 15.5×, and it cannot make one
pool's driving force swing 9× while the other's stays flat, because the two pools share a
single weight-independent `E_rev` and a single V.

So the specific thing the LIF gets wrong is **state-dependence**, and the results at risk are
exactly the ones that rest on state-dependence:

| result | exposure |
|---|---|
| tonic/phasic gain ratio (~7×; ~200× headroom) — **the central finding** | conductance ratio is from `gabaa_kinetics` and stands; but its *behavioural* consequence runs through current = g·(E−V), and the substrate suppresses the phasic pool's driving force by up to 9× at the moment it arrives. The ratio may be right and its circuit-level consequence still wrong. |
| NMDA respiratory arm (already **VOID**) | was tested on a substrate that cannot express the mechanism its premise names. VOID for a different reason than recorded — worth amending the provenance string either way. |
| efficacy ceiling (**VOID**) | saturation of a PAM depends partly on driving force collapsing as V → E_Cl. Never happens here. |
| selectivity ranking (**VALIDATED**, 20k draws) | ordinal, derived from occupancy, least exposed. Testing whether it survives a substrate change is the single best available check on it. |

That last row is the real prize and is argued in §6.

---

## 2. Scope, and why links 4 and 5 are one commit

**Link 4** — replace the LIF `Pop` with a conductance-based single-compartment cell.
**Link 5** — replace the hand-tuned preBötC with a published parameter set.

These cannot ship separately. The current cell is C = 200 pF, g_L = 10 nS (τ_m = 20 ms). The
Butera–Rinzel–Smith preBötC cell is ~C = 21 pF, g_L = 2.8 nS — roughly a 10× scale change in
input conductance. Every synaptic weight in `cpg.py`, `resp.py`, `rg2.py` and `circuit.py` is
in nS tuned against the 200 pF cell. Swap the cell alone and every weight is wrong by about
an order of magnitude; the network will not oscillate, or will oscillate for the wrong
reason, and the failure will be slow to diagnose because nothing raises.

So the cell and its weights are **one atomic unit**, which is also the whole point of link 5:
we are not rescaling hand-tuned weights, we are *replacing* them with weights that came from
the same publication as the cell. That is the credibility gain — "these are published and
independently validated" instead of "these produce a rhythm on my machine."

Mechanically this means the parameter registry must make a mixed pair **unconstructible**
(§5, E14).

---

## 3. The cell

**Choice: Butera–Rinzel–Smith (1999) Model 1** — J Neurophysiol 82:382–397 (part I, single
cell) and 82:398–415 (part II, network). Reasons:

- **I_NaP is the mechanism we need.** It is the current the LIF provably cannot carry (E7).
- **Only three state variables** — V, n (delayed-rectifier K activation), h (I_NaP
  inactivation). Na activation `m∞(V)` and I_NaP activation `mp∞(V)` are instantaneous, and
  Na inactivation uses the `(1 − n)` economy. This is far cheaper than full Hodgkin–Huxley
  and is the difference between a ~20× and a ~60× cost multiplier (§7).
- **Part II supplies network synaptic conductances** for a coupled heterogeneous population,
  so link 5's weights come from the same source as link 4's cell.
- **Same modelling lineage `cpg.py` already claims to follow.** Its docstring says
  "after Rybak/McCrea-style rat hindlimb CPG models"; Rybak-lineage locomotor neurons are
  HH cells with I_NaP, so one cell class serves both the respiratory and locomotor circuits.
  Currently that docstring claim is aspirational — the cells are LIF. This makes it true.

**Structure** (all parameter values below are from recall and **must be checked against the
papers before any code is written** — this project has already scored an invented range as if
it were data; see WORKLOG):

```
C dV/dt = -[ I_Na + I_K + I_NaP + I_L + I_syn + I_adapt? + I_app ]

I_Na  = g_Na  · m∞(V)³ · (1-n) · (V - E_Na)
I_K   = g_K   · n⁴            · (V - E_K)
I_NaP = g_NaP · mp∞(V) · h    · (V - E_Na)
I_L   = g_L                   · (V - E_L)

x∞(V) = 1 / (1 + exp((V - V_x)/σ_x))        for x in {m, n, mp, h}
τ_x(V) = τ̄_x / cosh((V - V_x) / (2 σ_x))    for x in {n, h}
dx/dt = (x∞(V) - x) / τ_x(V)
```

Recalled values, **UNVERIFIED**: C 21 pF; g_Na 28, g_K 11.2, g_NaP 2.8, g_L 2.8 nS;
E_Na +50, E_K −85, E_L −65 mV (some printings −57.5 — this one matters and must be read off
the paper, not recalled); (V_m, σ_m) = (−34, −5); (V_n, σ_n) = (−29, −4), τ̄_n 10 ms;
(V_mp, σ_mp) = (−40, −6); (V_h, σ_h) = (−48, +6), τ̄_h 10000 ms.

**Spike-triggered adaptation is retained but demoted.** It stays available as an optional
conductance so existing phenotype tests have a continuous path, but the acceptance criterion
in §4 is that the cell bursts with `g_adapt = 0`. Adaptation was a *stand-in* for I_NaP
inactivation plus Ca-dependent K plus synaptic depression; once I_NaP inactivation is real,
leaving adaptation on at its tuned strength double-counts burst termination. Default off for
the conductance backend.

**Integration: Rush–Larsen** (exponential Euler) on n and h — exact for a linear relaxation
with frozen τ and x∞, so the gates are unconditionally stable and dt is set by the accuracy
of the spike upstroke, not by gate stiffness. Minimum membrane time constant at full
activation is 21/(28 + 11.2 + 2.8 + 2.8) ≈ 0.47 ms, so forward Euler on V is *stable* up to
dt ≈ 0.9 ms but nowhere near *accurate*. Expect dt ≈ 0.05 ms (vs the current 0.1 ms) and
prove it with the dt-halving test (E13), never assume it.

**Noise becomes a current, not a voltage.** The LIF adds `sigma · √dt · randn` in mV, which
is a voltage kick whose physical meaning depends on C. For the conductance cell it should be
an injected current in pA so that the same noise parameter means the same thing across cells
with different capacitance. This is a small change with a sharp edge: it alters results on
the LIF backend too if applied there, so **apply it only to the new backend** and leave the
LIF path byte-identical.

---

## 4. Acceptance tests — written before the code

The point of listing these first is that this project's errors are overwhelmingly of the form
"plausible and wrong," and the only defence that has worked is a test that fails loudly. Each
of these is an observable the cell was **not fitted to**, so passing them is tier-promoting
evidence rather than a tautology.

**A1 — E7 inverted. The headline test.** E7 currently asserts that LIF voltage never exceeds
Vth, documenting *why* voltage-gated mechanisms are unavailable, with the note "If this ever
fails, voltage-gated mechanisms become available and the modelling choice should be
revisited." The conductance cell must do what E7 proves the LIF cannot: h-gate span > 0.5
over a burst cycle (measured LIF span: 0.01; needed: ~0.7), and V reaching > 0 mV on a spike.
The failure that documented the limit becomes the test that the replacement cleared it.
`test_E7_*` stays, parameterised by backend, asserting the LIF still behaves as recorded.

**A2 — intrinsic bursting with no adaptation.** A single isolated cell, `g_adapt = 0`,
I_app = 0, must burst. The LIF cannot (it requires spike-triggered adaptation). This is the
test that separates "rhythm from I_NaP" from "rhythm from a lumped stand-in," and it is the
reason the upgrade is not cosmetic.

**A3 — the published excitability sequence.** Raising tonic drive must move the isolated cell
quiescent → bursting → tonic spiking, with burst frequency increasing monotonically within
the bursting regime. Published, and not fitted by us.

**A4 — the riluzole/CdCl₂ dissociation. A new falsifiable prediction.** In vitro, blocking
I_NaP does **not** abolish the network rhythm, although it abolishes isolated-cell pacemaking
— the pacemaker-versus-network-rhythm controversy. So: with `g_NaP → 0`, A2 must fail
(isolated cell stops bursting) while the *network* rhythm persists. The LIF model cannot even
pose this test, because it has no I_NaP to block. If the network rhythm dies with I_NaP, the
network is a pacemaker-driven model and does not match the preparation — a real, informative
negative.

**A5 — f–I curve** of the isolated cell against the published curve.

**A6 — existing phenotypes survive, or fail for a stated reason.** Decide *in advance* what a
legitimate failure looks like. The strychnine phenotype in `rg2.py` (glycine removal destroys
alternation, rhythm persists) is a structural claim about group pacemakers and should survive
a substrate change; if it does not, that is a bug, not a discovery. The respiratory
`EUPNOEA_BAND` gate is a frequency band and will need re-anchoring to the new operating
point — a legitimate re-anchor. Pre-committing this distinction is what stops the
post-hoc rationalising that produced the retracted Matsuoka prediction.

**A7 — the LIF path is byte-identical.** Every current result must reproduce exactly on the
LIF backend after the refactor. Checked by hash of the trace arrays, not by eye.

---

## 5. Architecture

### 5.1 Two backends, one interface

`Pop` becomes a thin interface with `LIFPop` (the current code, unchanged) and `CondPop`.
Not for politeness — for four concrete reasons:

1. **Nothing existing breaks**, and A7 makes that checkable.
2. **The 20k-draw Monte Carlo stays affordable** (§7). Its output is an *ordering*, derived
   from occupancy; it does not need the spike waveform.
3. **Mechanism work runs on the conductance backend** where state-dependence matters.
4. **Both can be run and compared** — which converts the upgrade from "trust me, biophysics
   is better" into a measurable claim about which conclusions are substrate-dependent. This
   is the real payoff and it is argued in §6.

Selection is explicit per construction site, never a global default and never an env var. A
`substrate=` argument that must be passed, with the current LIF as the default, so an
un-updated caller keeps the behaviour it was validated against.

### 5.2 Parameter sets are atomic and carry provenance

A registry entry binds, as one unconstructible-if-mixed unit: cell parameters, synaptic
weights, E_rev values, dt, and a provenance string naming the publication and the preparation
(species, age, in vitro/in vivo). `CondPop` must refuse to be built with a weight table whose
provenance tag does not match its cell's. This directly prevents E14 (§8), which is the
failure mode §2 makes likely.

This also slots into the existing `config.CALIBRATIONS` and `results.Tier` machinery rather
than inventing a parallel mechanism — a published parameter set is a calibration with a
provenance and a tier like any other.

### 5.3 The synaptic interface must pass conductances *unevaluated*

Current contract: `Pop.step(dt, g, E, Idrive, rng)` where the caller has already called
`Syn.conductance(p.V)` — i.e. the NMDA Mg²⁺ block is evaluated at the **pre-step** voltage.

On the LIF that is harmless: V moves a few mV per step. On a conductance cell integrating a
spike in sub-steps of 0.05 ms, evaluating the Mg block once at the pre-spike voltage
**discards the entire voltage-dependent relief** — finding (a) in §1, the main NMDA benefit of
the upgrade, silently thrown away. The symptom would be "the upgrade didn't change the NMDA
result," which reads as a reassuring robustness check and is in fact the bug. That shape —
reassuring and wrong — is this project's signature failure.

Fix: `step` takes the raw conductance plus its receptor kind and applies the voltage factor at
its own current V, each sub-step. Keep the pre-evaluated dict form as a shim on the LIF path
for A7.

### 5.4 Chloride becomes available as a state variable

Not in scope for the first commit, but the design should not foreclose it, because two of this
project's open problems live here:

- **The developmental chloride switch** (NKCC1/KCC2, ~P11–P12) is currently a prose note that
  this project already over-corrected on once. A conductance cell with intracellular Cl⁻
  dynamics can represent it mechanistically.
- **Cl⁻ accumulation during sustained GABA-A activation self-limits a strong PAM**, which is
  exactly the saturation the currently-VOID efficacy ceiling is about.

Reserve the interface (E_rev as a per-population attribute rather than a module-level
constant `E_REV`), implement later.

---

## 6. What this is actually for: substrate-independence of the one validated result

The project's single VALIDATED headline is the calibration-independent selectivity ranking —
α5 arms beat non-selective BZ in 99.9% of 20,000 draws, 5th-percentile floor ≥ 2.49×. It is
calibration-independent. It is **not** substrate-independent, and nobody has checked.

Running the ranking on both backends is the only way to find out, and either outcome is worth
having:

- **Holds on both** → the result is substrate-independent as well as calibration-independent,
  which is a materially stronger claim than the one currently in the paper.
- **Flips or narrows** → the LIF result was partly an artifact of a substrate that suppresses
  the phasic pool's driving force by up to 9×. Better to find that ourselves.

**The comparison must be at matched operating points, not matched parameters.** §2 establishes
that the same nS weight means different things on the two cells, so comparing at equal weights
compares two differently-broken networks. Each substrate must first be independently anchored
to the *same observable* — control burst frequency and duty cycle — and only then is the drug
applied and the **fractional change from its own control** compared. Also: same random stream,
or the substrate gets credit for seed noise. Both of these are easy to get wrong and would
produce a confident, wrong answer.

Reduced draw count (~2,000) is sufficient for the conductance arm if the 5th-percentile floor
is the statistic of interest; the full 20k stays on LIF.

---

## 7. Compute budget

The constraint here is wall-clock and heat, not RAM — state grows from 2 arrays of n floats
per population to 3. Memory is a non-issue; duration is not.

**MEASURED: 4.1×** (`CondPop` vs `cpg.Pop`, 2 s simulated, same conductance dict, n = 1 and
n = 50 — 3.67× and 4.12×). Rule 1 below said to measure before committing, and the estimate
it was checking against was **~20×**, i.e. 5× too pessimistic.

The estimate reasoned "~10 transcendental evaluations per cell against the LIF's ~1, times 2
substeps." Both factors were real; the error was treating them as the whole cost. Per-step
overhead the two substrates share — the loop over receptor conductances, RNG draws, numpy
dispatch on small arrays — dominates at these population sizes, so the extra arithmetic is
marginal rather than multiplicative. The Butera three-state economy (V, n, h, with `m` and
`mp` instantaneous and Na inactivation carried by `1−n`) is what keeps it this cheap.

This materially loosens the budget. At 4×, running the conductance arm at the **full** 20,000
draws is plausible rather than out of reach, so §6's reduced-draw compromise should be
re-examined against a timed run before being accepted.

Rules, to be enforced rather than intended:

1. **Measure before committing.** Done, above. Re-measure for the network: the 4.1× figure is
   for one population stepped in isolation, and link 5 changes the shape of the work.
2. **The Monte Carlo stays on LIF** unless a timed run says otherwise — which, at 4.1×, it
   now might. The conductance arm runs at reduced draws (§6) only if the full count proves
   too slow in practice.
3. **The plant tolerates this.** `plant.py` runs 2 ms physics with `n_sub=20` neural substeps
   at 0.1 ms; dt = 0.05 means `n_sub=40`. A 2× cost on the motor path is fine.
4. **Hard rule:** if a script on the conductance backend exceeds its LIF wall-clock by more
   than the measured multiplier, it is a bug (runaway, non-convergence, or a dt that the
   integrator is silently fighting) — not an expected cost. Fail, do not wait it out.
5. Check memory before any parallel sweep, per the standing instruction, and keep the
   `__main__` guard in a real file — E8 bit this project twice, the second time through a
   heredoc-fed script that produced 13,334 tracebacks and 144 MB of output.

---

## 8. Predicted new failure modes, to be written as tests first

Based on the four review passes (38 defects at a flat rate of 9 / 11 / 8 / 10), the right
assumption is that this change will introduce defects of the project's characteristic shape:
silent, plausible, and reassuring. Naming them in advance is cheaper than finding them.

**E13 — silent dt-dependence.** An HH cell integrated with too large a dt does not blow up; it
produces a plausible but wrong firing rate. *Test:* run at dt and dt/4, require < 2%
difference in burst frequency and duty. This is the single most likely new defect.

**E14 — weight-scale carryover.** A nS weight tuned on the 200 pF cell survives into a 21 pF
network, and the result is a quiet 10× error. *Test:* constructing a cell and a weight table
with mismatched provenance tags must raise (§5.2).

**E15 — Mg block evaluated at the stale voltage.** §5.3. *Test:* assert the NMDA Mg relief
factor measured *during* a spike exceeds the resting value by > 5× on the conductance backend.
Pins the actual mechanism, not the plumbing.

**E16 — spike double-counting at the detector.** Upward crossing of a detection voltage
without a refractory lockout counts one spike several times on a broad spike, inflating every
rate. *Test:* an isolated cell driven to a known firing rate must report that rate; and spike
count must be invariant to the detection threshold across a range of plausible values.

**E17 — substrate difference attributed from an unmatched comparison.** §6. *Test:* the
comparison harness must assert identical random streams and matched control operating points
before it reports any difference.

**E18 — double-counted burst termination.** Adaptation left on at its LIF-tuned strength while
I_NaV inactivation is now real. Symptom: bursts terminate too early, which looks like a
plausible duty cycle. *Test:* `g_adapt = 0` must be the conductance backend's default, and A2
must pass at that default.

---

## 9. What this does **not** fix

Stating this because the adjacent error — thinking a substrate upgrade unblocks a calibration
— is easy and expensive.

- **It does not create a respiratory calibration anchor.** `prebotc_gaba_sens` is VOID because
  whole-body ventilation is the wrong *observable* for an isolated preBötC. That is a problem
  about what was measured, not about what the neurons are made of. Still wet-lab blocked on
  the > P12 muscimol concentration-response curve.
- **It makes the age confound structural, not incidental.** Butera parameters are neonatal
  rodent in vitro (~P0–P4). Adopting them imports that age, which means E_GABA must be the
  *neonatal* value — more depolarised than the −75 mV currently in `E_REV`, because the
  NKCC1/KCC2 switch has not yet happened. Consequences: (i) GABA-A becomes less inhibitory,
  possibly shunting-only, so the choice of parameter set is **not pharmacologically neutral**;
  (ii) the cell is at the wrong age for the > P12 muscimol anchor we are blocked on, so link 4
  moves that anchor *further* away, not closer. Both must be stated wherever the new backend
  is used.
- **It does not predict efficacy from structure.** Link 2, declined by design.
- **It does not give a subjective endpoint.** That needs an animal that learns.
- **It does not validate the locomotor weights** unless Rybak-lineage weights are adopted
  too — and those are for a different preparation, so that is a separate decision with its own
  provenance problem, not a freebie.

---

## 10. Sequencing

1. **Gate: one review pass finds nothing.** **Not met.** Passes 7 and 10 (2026-10-07) found
   eight and three, taking the total to **49 across eight passes**: 9 / 11 / 8 / 10 / 8 / 3.
   The rate is falling at last, but both passes found defects in classes already raised,
   and each restated the argument for this gate using this project's own code:

   * Pass 7 found a *second* instance of a defect class already hit and already "fixed" —
     a destructive database rebuild — where the guard added the first time was counting the
     wrong tables and so could not see the collateral it destroyed.
   * Pass 10 found a NaN latching the grid search in all three calibration sweeps (the
     scripts this project's calibration constants came from), and in one of them **the
     pass-7 fix is what made the path reachable**, by converting a loud `ZeroDivisionError`
     into a silent NaN.

   Three is not zero, and a fix that creates the next finding is the exact failure this
   gate exists to keep away from a conductance-based rewrite.
2. Verify the Butera parameters against the papers. Read them off; do not recall them. §3.
3. Write §4 acceptance tests and §8 failure-mode tests. They fail. That is correct.
4. `CondPop` as an isolated cell only. A1, A2, A3, A5, E13, E16 go green.
5. Parameter registry with provenance locking (§5.2). E14 green.
6. Synaptic interface change (§5.3). E15 green. A7 — LIF byte-identical — green.
7. Network: preBötC on published weights. A4 (riluzole dissociation) and A6.
8. The substrate-independence comparison (§6). This is the deliverable, not step 7.
9. Only then consider the locomotor circuit and chloride dynamics (§5.4).

Steps 4–6 are the bulk and have nothing to do with neurons; they are the reason this is
described as a large new surface. Step 8 is the only step that produces a result worth
publishing, which is an argument for not stopping at step 7 and calling the model improved.
