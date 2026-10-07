# Model specification

## Three modules — only two are simulated

| Module | Method | Rationale |
|---|---|---|
| Subjective match | **Empirical constraint** | Drug discrimination requires learning; a frozen policy cannot model it. The literature already defines the admissible region. |
| Motor impairment | Simulated | What embodied simulation is genuinely good at. |
| Respiratory depression | Simulated | preBotC models already exist with exactly this perturbation structure. |

Principle: simulation only where it has grounding; empirical data everywhere else.

## Motor module

**Plant.** The dm_control rodent (`rodent.xml`, the Aldarondo et al. 2024 body) is
installed and verified: 67 joints, 67 dof, 38 actuators, 0.3378 kg, 2 ms timestep,
~41,700 physics steps/s single-core (~83x realtime) on Apple Silicon.

**VERIFIED LIMITATION:** all 38 actuators are gaintype 0 — **torque motors, no muscles**.
No Hill-type force-velocity, no spindle or Golgi afferents, no monosynaptic reflex loop.
Reflex phenomena (hyperreflexia, clonus, tremor) therefore cannot emerge for the right
reasons in this plant. Proprioceptive-ish observations that DO exist: joints_pos/vel,
tendons_pos/vel, sensors_force/torque/touch, accelerometer, gyro, velocimeter.

**Controller.** Do NOT apply drug effects to a trained ANN policy (see reasoning trail).
Use a biophysical spinal CPG with identified cell types and transmitters. Candidate:
the two-layer rat hindlimb CPG model (rhythm generator + pattern formation layers,
validated against resetting and non-resetting deletions, with afferented vs. deafferented
perturbation recovery). The 38-muscle SimTK rat hindlimb model is the muscle-level asset.

**Porting decision (OPEN).** Existing rat motor models are OpenSim/Simbody, not MuJoCo.
Port to MuJoCo (costs ~a week; unifies with the rodent, gives fast parallel rollouts for
the ratio sweep) vs. work in OpenSim (no porting; slower, loses the rodent's trained
descending controller). Recommendation: port, because the ratio optimization needs many
rollouts.

**Readouts.** Gait variability, limb coordination, perturbation-recovery dynamics.
Use **probe perturbations** (pushes, slips, treadmill speed steps, uneven terrain) rather
than more baseline metrics — recovery dynamics discriminate sedation from ataxia from
weakness, which collapse together at baseline. Grade difficulty to keep controls at ~80%;
ceiling and floor effects are how behavioural pharmacology studies fail.

**Benchmark tasks verified working:** `rodent_escape_bowl`, `rodent_run_gaps`,
`rodent_maze_forage`, `rodent_two_touch` (160-310 env-steps/s; the egocentric camera
dominates cost — disable vision for motor assays). Map sim tasks to standardized assays
(rotarod -> narrowing beam, open field -> arena locomotion, gait -> treadmill footfall
timing) so there is comparable ground truth. Do NOT invent novel tasks.

## Respiratory module

**Framework.** preBotzinger complex models from the opioid-induced-respiratory-depression
literature are already structured as "apply a receptor-level perturbation, observe rhythm
degradation." Swap mu-opioid for GABA-A / NMDA.

**Inherited methodological finding:** opioid sensitivity was best predicted by **network
topology**, not cellular properties. Therefore run across topology ensembles, not one
network. (Degeneracy-as-ensemble, now empirically motivated rather than assumed.)

**Readouts.** Rhythm frequency, amplitude, failure threshold, and critically **whether
failure is graded or catastrophic** — i.e. does the dose-response plateau.

## Validation anchors

**Respiratory — the one that matters: the benzodiazepine/barbiturate ceiling dissociation.**
Benzodiazepines have a respiratory ceiling; barbiturates do not; the mechanism (PAM
requiring endogenous GABA and saturating, vs. direct agonist) is textbook. If the model
reproduces *why* one plateaus and the other does not, it can be trusted on a novel
compound's ceiling — the single property the product's safety depends on.

**Motor — four arms, two positive and two negative:** diazepam and ethanol (both
impairing), CPP (ataxia and reduced muscle tone at higher doses), Ro25-6981 (NR2B-selective,
did NOT potentiate motor impairment in combination). The last is the key one: it is the
positive control for the GluN2B-sparing hypothesis.

## Kill criteria

- No ratio achieves >10x margin on **both** axes -> this compound class cannot deliver
  the product; different mechanisms needed.
- Combination effects come out **supra-additive on respiratory depression** -> the mixture
  approach is more dangerous than single agents and the whole premise inverts.
- Model cannot reproduce the benzodiazepine/barbiturate ceiling dissociation -> no
  trustworthy forward model; stop before interpreting any novel prediction.

## Why a model at all

Combination effects are non-additive and endpoint-specific (triazolam + pregnanolone was
supra-additive for sedation yet showed *attenuated* ataxia; NR2B antagonist + morphine
showed no rotarod potentiation; CPP alone causes ataxia). Single-agent dose-response does
not predict combination margins, and the ratio space is too large to test empirically.
Triaging it is a legitimate job for simulation.

**The structural worry the model exists to quantify:** both the respiratory rhythm
generator and the locomotor CPG depend on glutamatergic excitation and GABA/glycine
reciprocal inhibition — the same two mechanisms being modulated for the wanted effect.
The margin may be intrinsically narrow, and regional/subunit selectivity is the only
thing that could open it.
