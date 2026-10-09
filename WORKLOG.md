# Work log — append-only

Running record of what was tried, what happened, and what was concluded. Purpose: survive
context compaction, prevent repeated steps, and prevent re-deriving dead ends.

**Conventions.** Append at the bottom. Never delete an entry — if a conclusion is later
overturned, add a new entry marked `CORRECTION` and leave the original. Every claim says
whether it is *literature*, *model-derived*, or *verified locally*.

Companion files: `KNOWLEDGE.md` (state + decisions), `knowledge/02-reasoning-trail.md`
(abandoned approaches and why), `data/pharmacology.db` (structured facts + citations).

---

## RECURRING ERRORS — read this before debugging anything

I have now made each of these more than once. Check them first.

| # | Error | Symptom | Fix |
|---|---|---|---|
| E1 | **Suprathreshold tonic drive** in a network that must switch or burst | population fires tonically; mutual inhibition or adaptation has no effect; `w` sweeps come out flat | drive must be **sub-rheobase** (~200 pA with gL=10 nS, tonic GABA 2 nS). Regenerative/switchable dynamics need a silent phase. **Made twice: locomotor RG, then preBötC.** |
| E2 | **Rate→conductance scaling omits presynaptic population size** | downstream population silent; conductance ~100× too small | a rate input represents a *population*. Multiply by afferent/unit count (`ia_n=25`) or absorb it in a documented lumped gain (`rg_gain`). **Made twice: RG→PF, then Ia→Mn.** |
| E3 | **Interneurons left without tonic drive** | they sit just below threshold (~−52 mV vs −50), barely fire, so the pathway they mediate silently does nothing; the weight sweep looks flat | give interneurons their own sub-rheobase drive (`drive_in`). Diagnostic: *if a weight sweep changes nothing, the presynaptic population is probably silent.* |
| E4 | **Relative threshold on a degrading signal** | drug that depresses a rhythm reports *higher* frequency (observed 13 Hz = 780 breaths/min) | never use `frac * signal.max()`. Use FFT peak (amplitude-invariant) + modulation index; declare failure on modulation/mean collapse, not crossing count. |
| E5 | **Ceiling effect masks drug effects** | all conditions ≈100% of control; readout sits exactly at a hard limit | check whether the readout equals a structural cap (Mn dynamic peak was 125.1 Hz with `tref=8 ms` → cap 125 Hz). Re-scale so controls sit mid-range. |
| E6 | **Over-reading a single seed / single operating point** | a dramatic effect that doesn't replicate | ≥4–6 seeds, paired per-seed ratios, and check the effect across ≥2 operating points before claiming it. |
| E8 | **Pool worker defined inside `__main__`** | `AttributeError: Can't get attribute '<fn>'` in every worker | macOS uses spawn; all multiprocessing worker functions must be at MODULE level |
| E9 | **Coupling test without a frequency mismatch** | correlation identical coupled vs uncoupled; `w` sweeps flat | give the oscillators different intrinsic frequencies, else they hold their arbitrary starting phase and the test is blind |
| E10 | **Testing a coupling by PARTIAL block** | 95% conductance reduction left entrainment intact (corr -0.49) | entrainment is cheap; remove the coupling COMPLETELY to test it |
| E11 | **Scaling a drug effect on conductance but not on kinetics** | sensitivity floors out (ventilation stuck at -25% at any `gaba_sens`) | a PAM changes peak conductance AND decay tau; scale both |
| E12 | **A calibration constant carried across a change of parameterisation** (regions, or pool decompositions) | applying the preBotC constant to the spinal cord abolished every reflex drug effect (all arms ~100% of control); later, applying the one-pool `gaba_sens` to the tonic/phasic split collapsed the drug effect to -2% and made every arm "survive" every overdose dose | re-calibrate against the anchor whenever the parameterisation changes. **Made twice: regions (session 5), then the tonic/phasic split (session 7d) — the second time ONE ENTRY after writing down that this migration would be needed.** A suspiciously reassuring result is the symptom to watch for. |
| E7 | **Voltage-gated gates on LIF neurons** | slow gate equilibrates at subthreshold mean and never moves (h span 0.01 vs needed 0.7) | LIF resets V at threshold, so V never reaches spike voltages. Use spike-triggered adaptation instead of voltage-gated inactivation. |

**General diagnostic heuristic learned here:** *a parameter sweep that produces identical
results across a wide range means the mechanism is not engaged at all* — look for a silent
population upstream, not for a subtler parameter.

---

## 2026-10-06 — Session 1

### Setup
- `/Users/whinz/Biochemistry`, uv venv pinned **Python 3.11** (3.13 wheels unreliable for
  dm_control). mujoco 3.15.0, dm_control, numpy 2.4.6, scipy 1.17.1. 330 MB total.
- Hardware: Apple M5, 10 cores, 32 GB. Verified during a 10-way sweep: 938% CPU,
  **0.52 GB RSS total**. Hardware is not the constraint; single-threaded Python is.
- Gotcha: **do not set `MUJOCO_GL=osmesa` on macOS** (invalid backend, raises at import).
- Gotcha: `dm_control.locomotion.walkers.assets` is a namespace package with no
  `__file__`; resolve paths via `dm_control.__file__`.

### Body model — VERIFIED
`rodent.xml` in dm_control *is* the Aldarondo et al. 2024 body. Measured: nq=67, nv=67,
nu=38, nbody=66, mass **0.3378 kg**, timestep 2 ms, **~41,700 physics steps/s (~83×
realtime)**. All four benchmark tasks run (`escape_bowl`, `run_gaps`, `maze_forage`,
`two_touch`) at 160–310 env-steps/s — the egocentric camera dominates that cost, so
disable vision for motor assays.

**Confirmed empirically (was previously only asserted):** `actuator_gaintype` histogram
is `{0: 38}` — all torque/servo, **zero muscle actuators**.

### Muscle port — DONE
`scripts/port_muscle.py` → `models/rodent_muscle.xml`. 24 antagonist muscle actuators
replace 12 hindlimb actuators (nu 38→50).

Three reasons the stock model is wrong for pharmacology (found by XML inspection — they
are **position servos**, `dyntype=filter`/`biastype=affine`, i.e. proportional feedback to
a setpoint):
1. flexor+extensor collapse into one signed command, so glycinergic reciprocal
   inhibition — our drug target — is unrepresentable;
2. the servo *compensates* for drug-induced weakness, masking the effect;
3. no force-length/velocity, so no substrate for spindle afferents.

**Port gotcha (cost a debug cycle):** the stock `<default><general gainprm="0.01"/>` leaks
into `<muscle>` (a `general` shortcut) and overwrites muscle `range[0]`, leaving a dead
zone at neutral posture where neither antagonist can pull. Must set `range` explicitly —
used `0.7 1.3`, centred so the force-length peak sits at the joint midpoint.

Verified by qpos sweep on knee_L: both antagonists pull across the whole range, each
strongest when lengthened (correct Hill force-length), equal and balanced at the midpoint,
co-contraction gives the expected stiffness profile.

### Locomotor rhythm generator — all-spiking version ABANDONED
Sequence of failures, each fixed and each revealing the next: no activity (E1) → silent
interneurons (E3) → weights too small (E2) → h gate frozen (E7).

Then an **exhaustive 288-point parallel grid search** (670 s, 10 cores): 288/288 produced
*a* rhythm, **0 reached the target box**. Best score 0.83 (0 = on target). Hard tradeoff:
near-target period (304, 262 ms) only with no alternation (corr −0.17, −0.20); good
alternation (−0.5) only at 81–111 ms; duty never above 0.24 vs physiological 0.3–0.7.
`w_gly` 8→40 changed nothing.

**Conclusion (structural, not parametric):** a LIF resets V every spike, so no plateau
potential forms, the half-centre is never bistable on the fast timescale, and adaptation
can never set the period. The oscillation is a *delayed-inhibition loop* oscillation at
~2×(5 ms AMPA + 8 ms gly + 20 ms membrane). **Do not retry tuning a spiking RG.**

### Locomotor RG — Matsuoka replacement, WORKS
`circuitpharm/rg.py`. Period law `≈ 2π√(τ·τₐ)` confirmed: predicted 688/889/1192 vs observed
**704/920/1299 ms**. Duty 0.48, alternation −0.90.

Free finding: at `w=4.0` the rhythm **stops entirely**. Reciprocal inhibition is
glycinergic, so a GlyR PAM (`glyr_gain≈1.8`) pushes w from 2.2→4.0 and abolishes the
rhythm. Ethanol potentiates GlyR → mechanistically grounded prediction.

### Hybrid circuit — WORKS
`circuitpharm/circuit.py`: Matsuoka RG → spiking PF/Mn/InPF/IaIn/Rc. Operating point
`rg_gain=900`, `gaba_tonic=2.0 nS`, Mn `tref=8 ms`. Output: **period 920 ms, Mn 35.8 Hz
(peak 123), duty 0.33, flex/ext corr −0.53** — inside the rat physiological box that 288
all-spiking combos could not reach.

Principle applied: biophysical detail **where the drug acts** (tonic extrasynaptic GABA-A
on PF/Mn, glycinergic IaIn/Renshaw, NMDA on RG→PF / PF→Mn / Ia→Mn), abstraction elsewhere.

### Drug interface
`spinal.cpg.Drug` — receptor- and subunit-resolved. `nmda_scale()` bounds a selective
antagonist's effect by `glun2b_fraction`, which *is* the selectivity window. Verified
numerically: 60% GluN2B-selective block removes **42% of NMDA conductance in forebrain
(2B=0.7) but only 9% in brainstem (2B=0.15)**; non-selective removes 60% everywhere.

### Open-loop drug panel — first results (model-derived)
`scripts/drug_panel.py`. Control: Mn 35.8 Hz, period 920 ms, duty 0.33.
- GABA-A PAM 1.5×/2.0×/3.0× → 77% / 50% / 17% of control. Efficacy ceiling (capped 2.5×)
  **limits but does not prevent** severe depression.
- NMDA 60% non-selective → 57%. NMDA 60% GluN2B-sel **spinal** → **93%**. Same block,
  forebrain profile → 66%.
- **Only** the ethanol-like condition (includes `glyr_gain=1.3`) changed the **period**
  (920→1336 ms). GABA-A PAM and NMDA antagonism changed amplitude/duty only.

### CORRECTION (supersedes the open-loop supra-additivity claim)
Open-loop single-seed suggested PAM+non-selective was strongly supra-additive (predicted
28.5%, observed 17%). **Multi-seed closed-loop shows this is operating-point dependent**:
strongly supra-additive near threshold (scale 0.22: predicted 64%, observed 36%) but
absent at higher gain (scale 0.40: predicted 61%, observed 63%). It is *not* a general
property. This was error E6.

### Closed loop — DONE
`circuitpharm/plant.py`. Mn pool rate → muscle activation (`d.ctrl`); actuator length/velocity
normalised against `actuator_lengthrange` → `spindle_ia()` → Ia afferents. Circuit takes
20 substeps of 0.1 ms per 2 ms physics step.

Spindle verified physiological: Ia resting **10.9 Hz**, dynamic peak **161.8 Hz**.

Hit E2 again (missing `ia_n=25` afferent count → reflex 25× too weak) and E5 (Mn dynamic
peak pinned at the 125 Hz refractory cap, so every drug read ~100% of control).

### Stretch reflex — 4 independent validations PASS (model-derived)
`scripts/reflex_panel.py`, 6 seeds, paired per-seed ratios, 3 operating points. None of
these directions were tuned for:

| condition | scale 0.22 | 0.30 | 0.40 | expected |
|---|---|---|---|---|
| Strychnine (GlyR block) | 155±18% | 122±11% | 116±3% | hyperreflexia ✓ |
| GABA-A PAM 2× (benzo-like) | 75±5% | 67±7% | 79±3% | muscle relaxant ✓ |
| NMDA 60% non-selective | 85±8% | 71±9% | 77±5% | depresses ✓ |
| NMDA 60% GluN2B-sel (spinal) | 99±5% | 98±4% | 93±4% | spared ✓ |

**The selectivity window reproduces on a second independent readout.**

Unexplained minor anomaly, logged not chased: GlyR PAM at scale 0.40 raised Mn *static*
firing 0.62±0.12 → 1.76±0.26 Hz (~4 sd) while depressing the dynamic response. Possible
Renshaw negative-feedback self-limitation. Small absolute quantities; not pursued.

### preBötC respiratory module — built, metric fixed, now bursting
`circuitpharm/resp.py`. Deliberately **all-spiking**, unlike the locomotor RG: a *group
pacemaker* (recurrently excitatory glutamatergic population + spike-triggered adaptation)
needs no plateau, so the mechanism works in LIF.

First tuning run reported a clean 1.58 Hz rhythm — **it was an artifact of E4**. Direct
trace inspection showed no bursts at all: fast low-amplitude ripple (FFT peak 3.5–3.67 Hz,
unchanged by any drug) around a tonic mean, because drive=240 pA was suprathreshold (E1
again). Drugs reported 400–660% of control frequency.

After rewriting `resp_metrics` (FFT peak + modulation index, failure on modulation/mean
collapse) and moving to sub-rheobase drive: **15/36 combos burst in the eupnoea range**.
Operating point `drive=170, g_adapt=2.5, tau_adapt=400, w_ee=0.45` →
**1.33 Hz = 80 breaths/min, modulation index 4.92**, deep bursting with silent
interburst phases.

### Current state
| Component | Status |
|---|---|
| rodent body | verified |
| muscle port | done, verified |
| locomotor circuit (hybrid) | working, physiological |
| drug interface | working, subunit-resolved |
| closed loop + spindle | done, verified |
| stretch reflex validation | 4/4 pass |
| preBötC | bursting at 80 breaths/min |
| respiratory drug panel | **in progress** |
| margin optimisation | not started |
| external validation vs gait data | not started |

### Next
1. Respiratory dose-response → apnoea thresholds (safety axis). **running now**
2. Compute the actual margin: forebrain effect (subjective proxy) vs spinal reflex
   depression (motor) vs respiratory depression, over the (PAM, NMDA block, selectivity)
   space — subject to the GABA≥NMDA salience constraint from drug discrimination.
3. External falsification test: model predicts ethanol lengthens step-cycle duration (via
   GlyR) while a PAM+NMDA mixture reduces amplitude *without* slowing the cycle. Checkable
   against existing rodent gait data, no new experiment.

---

## 2026-10-06 — Session 1 (continued)

### Respiratory drug panel — metric limitation found and bounded
`scripts/resp_panel.py`, operating point `drive=170, g_adapt=2.5, tau_adapt=400,
w_ee=0.45`. Control **1.30 Hz = 78 breaths/min**, mean output 28.0 Hz, modulation 4.91.

**Frequency is NOT a usable readout in the depressed regime** — a second, narrower form of
E4. Even after switching to the FFT peak (amplitude-invariant), frequency rises to
314–347% of control under depression, because once the burst pattern degrades the FFT peak
follows the faster residual fluctuations. Mean output and modulation index both fall
correctly and monotonically.

**Resolution (and it is the physiologically better choice anyway):** primary respiratory
readout is **mean output = minute-ventilation proxy** (rate x amplitude, which is what
matters clinically). Modulation index is the burst-integrity measure. Frequency is
reported nowhere.

Dose-response on ventilation (% of control):
- PAM 1.5x/2.0x/2.5x → 50 / 45 / 29. **PAM 3.0x and 4.0x also 29** — the ceiling bites.
- NMDA non-sel 30/50/70/90% → 67 / 50 / 43 / 46
- NMDA GluN2B-sel 30/50/70/90% → **93 / 90 / 86 / 82**

### Ceiling dissociation — VALIDATION ANCHOR PASSES
`scripts/ceiling_test.py`. The anchor specified in `knowledge/03-model-spec.md` for the
respiratory module: does the model reproduce *why* a PAM plateaus and a direct agonist
does not? Paired comparison, identical in every respect except
`gaba_a_efficacy_cap` (2.5 vs effectively infinite).

| dose | PAM ventilation | direct agonist |
|---|---|---|
| 1.0x | 100% | 99% |
| 2.5x | 29% | 29% |
| 3.0x | 29% | 25% |
| 4.0x | 28% | 18% |
| **6.0x** | **27%** | **APNOEA** |
| 9.0x | 24% | APNOEA |
| **14.0x** | **20%, still breathing** | APNOEA |

Apnoea threshold: PAM **not reached at 14x**; direct agonist **6x**. The arms are
identical below the cap and diverge above it, exactly as the mechanism predicts.

**This is the benzodiazepine-vs-barbiturate overdose safety difference reproduced from
first principles.** It was the stated criterion for trusting the model on a novel
compound's ceiling, which is the single property the product's overdose safety depends on.
Strongest validation result in the project so far.

Note the ceiling bounds but does not abolish harm: the PAM plateaus at ~20-29% of control
ventilation. Severe, but not apnoeic at any dose. A direct agonist (gaboxadol, muscimol,
barbiturate-like) is therefore excluded on overdose-safety grounds independently of its
drug-discrimination failure.

### Margin calculation — running
`scripts/margin.py`. Compares **cost at matched subjective effect**, not matched dose,
because a selective and a non-selective antagonist at equal occupancy produce different
forebrain effects.

Design: fix a forebrain NMDA conductance reduction (15/25/35/45%, GluN2B fraction 0.7),
solve for the `nmda_block` each arm needs, then measure simulated respiratory and
stretch-reflex cost with peripheral GluN2B fraction 0.15. GABA-A component fixed at
PAM 2.0x in every arm so the comparison isolates the NMDA selectivity choice.

Arithmetic of the window, before simulation: to reach 30% forebrain reduction a selective
antagonist needs `block = 0.30/0.70 = 0.43`, which removes `0.43 x 0.15 = 6.4%` of
peripheral NMDA; a non-selective antagonist needs `block = 0.30`, which removes 30%
peripherally. **~4.7x less off-target exposure at matched on-target effect.**

**Stated weakness:** the subjective axis is a linear receptor-level index, not a circuit
model, because drug discrimination requires learning and cannot be simulated with a frozen
policy. All margin numbers are therefore "cost at a *nominal* matched subjective effect".

---

## 2026-10-06 — Session 1, final block

### Bug: multiprocessing worker must be module-level
`scripts/margin.py` first run died with `AttributeError: Can't get attribute 'dispatch'`
in all 38 workers. macOS multiprocessing uses **spawn**, so a function defined inside
`if __name__ == "__main__"` is not picklable. Fixed by hoisting `dispatch` to module level.
Added to the recurring-error class: **all Pool worker functions at module level.**

### FINDING 15 IS FALSIFIED — and the cause is an architectural error

Finding 15 predicted ethanol lengthens step-cycle duration via GlyR potentiation, while a
PAM+NMDA mixture reduces amplitude without slowing the cycle. Two literature checks:

1. **Gait (confounded, points the wrong way).** Ethanol in increasing doses *increased*
   stride frequency and *decreased* stride length in D2 mice (not B6) on ventral-plane
   videography. Increased frequency = SHORTER cycle, opposite to prediction. BUT this was
   a motorised treadmill at fixed belt speed, where `speed = stride_length x frequency`
   is a kinematic identity — if ethanol shortens stride length, frequency MUST rise to
   keep up with the belt. So this cannot measure the intrinsic CPG period and is not
   decisive. Also strain-dependent (D2 but not B6).

2. **Fictive locomotion (decisive, and it refutes the mechanism).** In lamprey fictive
   swimming, strychnine (glycine receptor BLOCK) **increased cycle period**, and
   critically: strychnine **eliminates left-right alternation while robust rhythmic
   activity persists** — same burst proportion, same rostro-caudal coordination, just
   left-right co-activation instead of alternation.

**The architectural error.** In the real cord, glycinergic reciprocal inhibition sets the
PHASE RELATIONSHIP; it does not generate the rhythm. Remove it and the rhythm continues.
In my Matsuoka RG, mutual inhibition IS the oscillator — remove it (or over-strengthen it)
and oscillation stops. My model conflates rhythm generation with phase setting, which is
exactly what the two-layer Rybak architecture exists to separate.

Consequence: the `w = 4.0 -> rhythm stops` result and finding 15 both rest on a pathway
that does not work that way in biology. **Retract finding 15.**

**The fix is already built.** The preBotC module uses the correct architecture: a group
pacemaker (recurrent excitation + spike-triggered adaptation) generates the rhythm with no
inhibition required, so blocking inhibition would leave the rhythm intact and only disturb
phase. The locomotor RG should be rebuilt on that same mechanism, with glycinergic
reciprocal inhibition demoted to a phase-setting layer that does not control the period.
Validation target for the rebuild: **blocking glycine must abolish alternation while
preserving rhythm** — a direct, published, falsifiable check.

**Scope of the damage — limited, and this matters.** The margin results are NOT affected:
the reflex arm runs with `rg_gain=0` (locomotor RG entirely off) and the respiratory arm
uses the group-pacemaker architecture. Affected: the open-loop locomotor drug panel
(period/duty/alternation rows) and finding 15.

### MARGIN RESULT — negative signal for the candidate, but confounded
`scripts/margin.py`. Controls: ventilation 27.8 Hz, reflex gain 0.306. GABA-A fixed at
PAM 2.0x in every arm; cost measured at matched *forebrain* NMDA reduction.

| forebrain target | arm | nmda block | peripheral NMDA cut | ventilation % | reflex % |
|---|---|---|---|---|---|
| 15% | GluN2B-sel | 0.21 | 3.2% | 43 | 66 |
| 15% | non-selective | 0.15 | 15.0% | 35 | 65 |
| 25% | GluN2B-sel | 0.36 | 5.4% | 41 | 66 |
| 25% | non-selective | 0.25 | 25.0% | 32 | 63 |
| 35% | GluN2B-sel | 0.50 | 7.5% | 39 | 64 |
| 35% | non-selective | 0.35 | 35.0% | 32 | 58 |
| 45% | GluN2B-sel | 0.64 | 9.6% | 38 | 64 |
| 45% | non-selective | 0.45 | 45.0% | 31 | 58 |

Selective advantage: only **+7 to +9 points ventilation, +1 to +6 points reflex.**

**Why so small.** The GABA-A arm dominates and has no selectivity window. PAM 2.0x alone
already puts ventilation at ~45% and reflex at ~67%; the NMDA component then moves things
only a little, so the (large, real) GluN2B selectivity advantage is swamped. This is the
literature-derived tension, now quantified: the drug-discrimination constraint
(GABA salience >= NMDA salience) forces mass onto the arm with no selectivity escape.

**And the absolute level is alarming:** ventilation 31-43% of control at a merely nominal
subjective effect, with no dose headroom anywhere in the tested space.

### THE PROBLEM I CANNOT RESOLVE FROM THE MODEL
The preBotC **GABA-A sensitivity is uncalibrated.** PAM 1.5x -> 50% ventilation; 2.0x ->
45%. Real benzodiazepines at anxiolytic or recreational doses do not halve ventilation in
humans or rodents. So either:
  (a) the candidate genuinely cannot deliver a subjective effect without severe
      respiratory depression (real negative result), or
  (b) my preBotC is several-fold too GABA-A-sensitive (model calibration artifact).

**I cannot distinguish these from inside the model.** The ceiling dissociation validated
the *shape* of the dose-response (plateau vs runaway) but says nothing about absolute
sensitivity. Resolving it requires external calibration: measured benzodiazepine (or
propofol) dose vs minute ventilation in rodents, fitted to set the preBotC GABA-A gain.
Until that is done, **no margin number here should be believed in absolute terms** — only
the selective-vs-non-selective *contrast*, which is a paired comparison and survives a
common calibration error.

---

## 2026-10-06 — Session 2: calibration, architecture fix, first viable combination

### CALIBRATION DONE — the (a)/(b) question from Session 1 is answered: it was (b)

Literature anchors found:
- **Human midazolam 2 mg IV** (clearly sedative; in the range that substitutes for ethanol
  in discrimination): minute ventilation **-14.3 +/- 5.9%** (second cohort -19 +/- 7%;
  younger subjects only -6 to -9%), tidal volume -16 to -22%, and a **COMPENSATORY
  respiratory rate INCREASE of ~10%**.
- **Propofol** at a sedative dose: tidal volume ~**-60%**.

Uncalibrated model gave -50% ventilation at PAM 1.5x => **3-4x too GABA-A-sensitive.**

**New parameter `gaba_sens`** = fraction of preBotC GABA-A conductance that is
drug-modulatable. Biologically motivated, not a fudge: BZ-site PAMs require gamma2, and
the respiratory network also expresses delta, alpha4 and **epsilon** subunits, with
epsilon (enriched on NK1R+ rhythm-generating neurons) BZ-insensitive.

**Sub-bug found during fitting:** scaling only the conductance left the GABA-A decay
prolongation (`gaba_a_tau`) unscaled, which floored achievable sensitivity -- ventilation
stuck at -25% however low `gaba_sens` went. Sensitivity must scale BOTH peak conductance
AND tau. Centralised in `Syn(..., sens=)`. Added to recurring errors.

**Fitted: gaba_sens = 0.15.** Validation at two anchors:

| | model | measured |
|---|---|---|
| capped PAM 2.0x, ventilation | **-17.3%** | -14 to -19% |
| capped PAM 2.0x, tidal volume | **-20.8%** | -16 to -22% |
| capped PAM 2.0x, rate | +25% (direction right, magnitude ~2x high) | +10% |
| uncapped 4.0x, amplitude | **-70.5%** | propofol ~-60% |

Arms identical at 2.0x (below the cap) and divergent above -- the ceiling dissociation
survives calibration. Propofol must NOT be modelled as a capped BZ-site PAM (it also hits
glycine receptors, HCN and Na channels and has direct agonist activity).

### FIRST VIABLE COMBINATION (calibrated)
`scripts/optimise.py`, 270 sims. Subjective index = gaba_term + nmda_term with
gaba_term=(gaba_scale-1)/1.5, nmda_term=forebrain_reduction/0.5, constraint
gaba_term >= nmda_term (the salience requirement).

| subjective threshold | best feasible combination | ventilation | reflex |
|---|---|---|---|
| >=0.6 | PAM 1.8x + 10% forebrain NMDA, GluN2B-sel | 83% | 79% |
| >=0.8 | PAM 1.8x + 20% forebrain NMDA, GluN2B-sel | 82% | 75% |
| >=1.0 | PAM 2.2x + 10% forebrain NMDA, non-sel | 71% | 67% |
| **>=1.2** | **PAM 2.2x + 20% forebrain NMDA, GluN2B-sel** | **75%** | **66%** |
| >=1.4 | PAM 2.2x + 40% forebrain NMDA, GluN2B-sel | 73% | 66% |

**The central tension, fully quantified.** The GluN2B selectivity advantage GROWS as weight
shifts to the NMDA arm, but the salience constraint forbids exactly that shift:

| configuration | subj | 2B-sel vent | non-sel vent | advantage | salience |
|---|---|---|---|---|---|
| PAM 2.2 + 20% fb | 1.20 | 75% | 62% | +13 pts | OK |
| PAM 1.8 + 40% fb | 1.33 | 75% | 49% | +26 pts | VIOLATES |
| PAM 1.0 + 40% fb | 0.80 | **91%** (reflex **103%**) | 56% | +35 pts | VIOLATES |

So a GluN2B-selective antagonist ALONE at 40% forebrain block has essentially **no
physiological cost** (ventilation 91%, reflex 103% = no motor impairment at all) -- but it
is NMDA-only and per the discrimination literature would not feel like alcohol. That is the
"different pleasant state" option flagged early in the project, now with numbers.

### LOCOMOTOR RG v2 — ARCHITECTURE FIXED AND VALIDATED
`circuitpharm/rg2.py`. Each half-centre is an independent group pacemaker (recurrent excitation +
spike-triggered adaptation, the preBotC mechanism); glycinergic interneurons only enforce
anti-phase. Operating point: drive=260, g_adapt=1.2, tau_adapt=280, ie_gly=3.0,
ee_ampa=0.55, asym=0.08 -> period 1245 ms, duty ~0.3, corr -0.51.

| condition | corr | period F | period E |
|---|---|---|---|
| coupled (ie_gly=3.0) | -0.51 | 1245 | 1245 |
| glycine 5% residual | -0.49 | 1240 | 1241 |
| **glycine fully OFF** | **-0.09** | **1156** | **1245** |
| uncoupled, identical oscillators | **+0.99** | 1192 | 1192 |

**Validation target MET**: with coupling both half-centres entrain to a common period in
anti-phase; remove it and periods diverge, alternation is lost, and the rhythm PERSISTS.
Identical uncoupled oscillators co-activate (+0.99) -- the strychnine phenotype.

Two testing errors found and fixed along the way, both now in the recurring-error table:
- **E9** A coupling test needs an intrinsic FREQUENCY MISMATCH between the oscillators.
  Without one, two identical oscillators hold their arbitrary starting phase offset and the
  correlation is the same coupled or uncoupled (observed -0.55 both ways) -- the test cannot
  detect the thing it exists to test.
- **E10** **Entrainment is cheap**: it survived a 95% reduction in coupling conductance
  (corr -0.49 at glyr_gain=0.05). Testing a coupling by partially reducing it fails; it must
  be removed completely. This is why several `w_gly` sweeps came out flat.

**Revised prediction replacing retracted finding 15:** in v2 the period is set by
adaptation, so glycine potentiation changes PHASE/coordination, not cycle period -- and
since entrainment saturates at weak coupling, a GlyR PAM should have little effect on
locomotor timing at all. Testable against gait data as inter-limb COORDINATION indices
rather than cycle duration.

---

## 2026-10-06 — Session 2 conclusion: RECOMMENDED COMBINATION

### Final confirmation (8 seeds, tight error bars)
`scripts/confirm.py`. Controls n=8: ventilation 27.74, amplitude 137.0, reflex gain 0.302.

| condition | ventilation | amplitude | reflex |
|---|---|---|---|
| control | 100±1% | 100±1% | 100±8% |
| *reference:* sedative benzodiazepine alone (PAM 2.0x) | 83±1% | 79±0% | 68±3% |
| **CANDIDATE A** BZ-PAM 2.2x + 20% forebrain GluN2B-sel | **75±1%** | 71±0% | **63±3%** |
| same but NON-selective NMDA | 62±1% | 58±1% | 60±5% |
| **CANDIDATE B** neurosteroid 1.8x + 20% forebrain GluN2B-sel | **82±0%** | 78±0% | **75±5%** |
| **OPTION C** GluN2B-selective alone, 40% forebrain | **90±1%** | 88±0% | **99±7%** |

Readings:
- Candidate B is **respiratorily equivalent to a benzodiazepine alone** (82% vs 83%) while
  delivering the full alcohol-like subjective index, and is BETTER than a benzodiazepine on
  motor impairment (75% vs 68%).
- Selectivity is worth **13 points of ventilation** (A 75% vs non-selective 62%).
- Option C is nearly harm-free (ventilation 90%, reflex 99% = no motor impairment) but is
  NMDA-only and so violates the salience constraint -- it would not feel like alcohol.

### Overdose test — the decisive safety result
`scripts/overdose.py`, dose multiplier applied to the WHOLE product.

| dose | BZ-PAM | neurosteroid (cap 6) | uncapped (barbiturate-like) |
|---|---|---|---|
| 5x | 69% | 61% | 54% |
| 12x | 66% | 59% | 31% |
| 26x | **61%** | **57%** | **8%** |

Both capped arms PLATEAU and remain survivable to 26x; the uncapped reference keeps
degrading toward death. Alcohol's own intoxicating:lethal ratio is roughly 5x (20 mM strong
intoxication, ~100 mM lethal), so **both candidates eliminate the overdose-death harm** --
one of the four physiological harms this project set out to remove.

### GABA-arm choice
`scripts/final_combo.py`, 648 sims, swept the assumed subjective efficiency of the
neurosteroid arm because that assumption is only qualitatively supported:

| assumed neurosteroid subj. efficiency | best neurosteroid combination | vent | reflex |
|---|---|---|---|
| 1.0x (no advantage) | PAM 2.2 + 20% fb | 76% | 66% |
| 1.3x | PAM 1.8 + 30% fb | 80% | 74% |
| 1.5x | PAM 1.8 + 20% fb | 82% | 75% |
| 1.8x | PAM 1.8 + 20% fb | 82% | 75% |

Advantage is +6 ventilation / +9 reflex and is **robust for any efficiency >= 1.3x**. Basis:
ethanol is more likely to substitute for neuroactive steroids than for benzodiazepines, and
ganaxolone's preclinical sedation is "comparable to midazolam" so there is no extra
respiratory penalty per unit effect.

Counterweight: brexanolone carries a boxed warning for excessive sedation and SUDDEN LOSS
OF CONSCIOUSNESS (4% vs 0% placebo) with mandatory continuous pulse oximetry -- consistent
with high-efficacy neurosteroids having direct agonist activity and therefore LOSING the
ceiling, which finding 22 established as the single most important safety property.

### RECOMMENDATION
**Primary: Candidate B, built with a LOW-EFFICACY neurosteroid** -- one without direct
agonist activity, so the overdose ceiling is preserved while the subjective efficiency
advantage is kept. Target ~1.8x GABA-A potentiation equivalent.
**Fallback: Candidate A** (benzodiazepine-site PARTIAL agonist, ceiling empirically proven
for the class) if a ceiling cannot be assured for the neurosteroid.
**NMDA arm either way: GluN2B-selective antagonist** at ~20% forebrain NMDA conductance
reduction (= ~29% receptor block given forebrain GluN2B fraction 0.7). Existing tool
compounds: ifenprodil, traxoprodil (CP-101,606), Ro25-6981, radiprodil, rislenemdaz.
**Separate product option: Option C**, a GluN2B-selective antagonist alone -- not
alcohol-like, but near-zero motor and respiratory cost.

### What is NOT established
- The subjective axis is never simulated. It is a receptor-level index constrained by the
  drug-discrimination literature. No model here can tell you what the compound feels like.
- Hepatotoxicity and carcinogenicity are out of scope for every model in this repo and
  need a conventional toolchain. They are where the top-four harms actually get removed.
- Dependence liability is explicitly out of scope by user decision and is NOT addressed.
  Any GABA-A PAM will produce tolerance and a withdrawal syndrome with chronic use.
- The spinal GABA-A sensitivity was never calibrated (only the preBotC was). Reflex numbers
  are qualitatively validated (4/4 anchors) but not quantitatively anchored.
- Formulation, dosing and PK matching are regulated drug development, deliberately not done.

---

## 2026-10-07 — Side-effect review: the recommendation changes

### SERIOUS: the BZ-site partial agonist fallback is in a class with a hepatotoxicity record
- **Ocinaplon**: anxiolytic without typical BZ side effects in a 4-week GAD trial, then
  development HALTED for hepatic toxicity (terminated early after liver problems in one
  patient).
- **Alpidem**: reached market as an anxiolytic, then WITHDRAWN for hepatotoxicity.
- Bretazenil / abecarnil / pazinaclone: failed clinically for other reasons.

Probably idiosyncratic DILI rather than mechanism-based (zolpidem, from alpidem's own
chemical family, is fine), but liver injury is THE harm this project exists to remove, so a
class prior like this is disqualifying without compound-level screening.

**EFFECT ON THE RECOMMENDATION: inverts the ranking.** The neurosteroid arm is now primary
for a SECOND independent reason -- brexanolone, zuranolone and ganaxolone are all approved
with no hepatotoxicity signal. The BZ-site partial agonist is no longer a clean fallback;
it is a per-compound DILI screening problem.

### GluN2B arm: the model was partly wrong, plus a liability it cannot see
| compound | outcome |
|---|---|
| Traxoprodil (CP-101,606) | discontinued for **QT prolongation**; dose-related dissociation and amnesia; depersonalisation/confusion at high dose |
| Rislenemdaz (CERC-301) | safety pharmacology and neurotoxicity clean; failed on EFFICACY only |
| Ro 25-6981 | safe for further development in infantile rats |
| Ifenprodil | alpha1-adrenergic off-target -> hypotension; poor selectivity |

**CORRECTION to the model's prediction.** Finding 26's mechanism (GluN2B-selective block
spares GluN2D-enriched PV interneurons -> no disinhibition -> non-psychotomimetic) is
directionally right but QUANTITATIVELY TOO OPTIMISTIC: traxoprodil produced ketamine-like
dissociation at higher doses, described as a *favourable* profile relative to unselective
antagonists, not an absent one. Expect reduced, not zero, dissociation.

**QTc is invisible to every model in this repo** -- a hERG off-target property. It killed
traxoprodil and is now a mandatory screening gate.

Encouraging: rislenemdaz cleared safety pharmacology and failed only on efficacy, which is
irrelevant here (we want a subjective effect, not an antidepressant one). The class can be
clean.

### Neurosteroid arm side effects
Somnolence, dizziness, sedation, dry mouth, flushing. Zuranolone: boxed warning for
IMPAIRED DRIVING. Brexanolone: boxed warning for excessive sedation and sudden LOC.

Mechanistically specific withdrawal (relevant even though dependence is out of scope):
neurosteroid withdrawal causes a **threefold increase in GABA-A alpha4-subunit expression**
via an Egr3 pathway, with increased anxiety and increased seizure susceptibility, explicitly
likened to alcohol/benzodiazepine/barbiturate withdrawal. Critically it also causes a
**significant reduction in the antiseizure efficacy of diazepam** -- so if a withdrawing
user seizes, the standard rescue drug works less well. That is a treatment problem, not just
a user problem, and should be stated regardless of scoping.

Also: neurosteroids interconvert with progesterone metabolites -> plausible endocrine
effects on chronic use, which no model here touches.

### One FAVOURABLE interaction
GABA-A agonists are protective against NMDA-antagonist neurotoxicity (Olney's lesions,
rodent retrosplenial/posterior cingulate, seen with channel blockers). GluN2B-selective
antagonists are not associated with those lesions anyway, so selectivity and GABA
co-administration push the same direction. The COMBINATION should be safer on this axis
than the NMDA component alone -- a rare case where mixing helps.

### Combination risks outside the models
- Additive sedation and driving impairment (both components; both approved neurosteroids
  carry driving warnings).
- **Users will mix it with actual alcohol.** Both components are additive with ethanol at
  the same targets. This is the most likely real-world serious-harm route and ceiling
  engineering does not fix it.
- Carcinogenicity: genuinely favourable -- the harm came from acetaldehyde and neither
  component is metabolised to it.

### Three screening gates, none evaluable by any model here
1. hERG/QTc on the NMDA candidate (killed traxoprodil).
2. DILI screening on the GABA candidate, whatever its class.
3. Human dissociation dose-response: find the dose where dissociation sits below the
   subjective-effect dose.

---

## 2026-10-07 — Session 3: compound survey, the alpha5 finding, and the simulator

### THE KEY FINDING: alpha5 carries ethanol's subjective effect; alpha1 carries the harm

This corrects a conclusion I got wrong twice (findings 20, 29) and it reframes the project.

| subtype | role | regional distribution |
|---|---|---|
| **a5** | **mediates ethanol's DISCRIMINATIVE STIMULUS** -- a5 agonists (QH-ii-066, panadiplon) mimic it; a5 inverse agonist L-655,708 blocks it | **forebrain-restricted**: hippocampus (>25% of CA1/CA3 neurons) and olfactory bulb highest; **pons and medulla much lower** |
| a2/a3 | also substitutes -- HZ-166 and YT-III-31 substituted FULLY for ethanol in rhesus **without altering response rates** | widespread |
| **a1** | **sedation, RESPIRATORY DEPRESSION, sleep, ataxia, motor impairment**; a1 antagonists do NOT attenuate ethanol's stimulus | **predominant adult preBotC subtype**; medulla |
| delta | NOT necessary for ethanol's stimulus | respiratory network expresses it |

I had previously concluded "the GABA arm has NO selectivity window" because the preBotC
expresses a1, a2/3, delta, a4 and epsilon. That was the wrong question. The right question
is which subtype carries the SUBJECTIVE effect -- and it is a5, which the medulla barely
expresses. The window was there all along, pointing the other way.

Confirming compound: **MP-III-022** is 300% efficacious at a5b3g2, **nonmodulatory at a1**,
and at 1-10 mg/kg is "devoid of ataxia, sedation or an influence on the anxiety level".
Contrast QH-ii-066 (older, less clean a5 agonist) which DOES cause sedation and ataxia --
consistent with residual a1 activity being the culprit.

### MODEL UPGRADE: subunit-resolved GABA-A (`circuitpharm/subtypes.py`)
Regional subunit fractions x per-compound subtype efficacies, with PER-REGION calibration
constants anchored so a non-selective benzodiazepine reproduces each region's empirical
value (preBotC 0.15 from midazolam ventilation; spinal 1.00 where the reflex validations
passed).

**Bug caught by this:** my first version used ONE global constant, which asserted a spinal
sensitivity of 0.166 and silently abolished every reflex drug effect (all arms came out at
~100% of control, contradicting the earlier validated benzodiazepine result of 67-79%).
Added as recurring error E12: *one calibration constant cannot serve two regions whose
empirical anchors differ by an order of magnitude.*

Resulting selectivity ratios (subjective effect per unit respiratory burden):

| compound | subj/resp |
|---|---|
| ideal a5-only PAM | 77.0 |
| **MP-III-022 (a5-selective)** | **40.8** |
| L-838417 (a2/a3/a5, a1 antagonist) | 16.6 |
| HZ-166 / KRM-II-81 (a2/a3) | 8.6 |
| non-selective benzodiazepine | 3.7 |
| **neurosteroid** | **2.9** |
| zolpidem (a1-preferring) | 0.6 |

### NEW RECOMMENDATION (8 seeds, % of drug-free control)

| condition | ventilation | stretch reflex |
|---|---|---|
| drug-free control | 100±1% | 100±8% |
| ref: non-selective benzodiazepine, subj 0.50 | 83±1% | 72±5% |
| ref: neurosteroid, subj 0.50 | 80±1% | 67±5% |
| **a5-selective PAM ALONE, subj 0.50** | **100±0%** | **98±7%** |
| **RECOMMENDED: a5-sel PAM 0.50 + GluN2B 20%** | **92±1%** | **96±6%** |
| a5-sel PAM 0.50 + GluN2B 35% | 89±1% | 94±7% |
| same but NON-selective NMDA | 73±0% | 89±7% |
| *previous recommendation* (neurosteroid + GluN2B) | 76±0% | 68±5% |
| a2/a3 PAM 0.35 + GluN2B 20% | 90±1% | 90±7% |

An a5-selective PAM alone at a substantial subjective effect is **indistinguishable from
drug-free on both axes**. The new combination improves on the previous recommendation by
16 points of ventilation and 28 points of reflex.

**Constraint:** a ceiling-limited a5-selective compound saturates at subjective ~0.51,
because it carries the whole GABA load on one subtype. Stronger effects need the NMDA arm
or a higher ceiling (which forfeits overdose protection).

### NEUROSTEROIDS DISQUALIFIED as primary (reverses the Session 2 recommendation)
Two independent reasons:
1. They **differ in potency but NOT in intrinsic efficacy**, so there is no "low-efficacy
   neurosteroid" to select -- which was exactly what my Session 2 recommendation specified.
2. They are **subtype non-selective** and act at a distinct site that also engages delta
   and epsilon, so they cannot spare a1. Their subj/resp ratio (2.9) is WORSE than a plain
   benzodiazepine (3.7).

### HEPATOTOXICITY: four GABA-A anxiolytics, four unrelated chemotypes
| compound | outcome |
|---|---|
| Alpidem | reached market, WITHDRAWN for hepatotoxicity |
| Ocinaplon | HALTED; liver injury in one patient |
| Panadiplon | HALTED Phase 1; **human hepatic toxicity NOT predicted by rat, dog OR monkey**; mitochondrial/Reye's-like via a carboxylic acid metabolite |
| Kava / kavalactones | banned/withdrawn in several EU markets; transient enzyme rise to fulminant failure and death; **excessive alcohol intake is a listed risk factor** |

Probably idiosyncratic rather than mechanism-based (zolpidem, alpidem's own chemical
family, is fine), but the base rate among novel non-benzodiazepine GABA-A PAMs reaching
the clinic is alarming, and panadiplon shows standard preclinical species can miss it.
Kava is the directly analogous warning: a GABA-A-active hepatotoxin consumed by drinkers.

### NMDA ARM: better options than GluN2B found
- **Esmethadone (REL-1017)** -- Phase 3. Low-affinity low-potency uncompetitive blocker,
  the **opioid-INACTIVE d-isomer of methadone**. NO meaningful abuse potential in
  recreational drug users (tested against oxycodone AND ketamine), NOT neurotoxic in rats,
  no reinforcement/dependence/withdrawal. A clean enantiomer separation -- its mirror image
  is an opioid, which closes the loop on this project's original chirality thread.
- **Lanicemine (AZD6765)** -- low-trapping channel blocker (54% vs ketamine's 86%);
  150 mg IV in 22 patients with no psychosis or dissociation. Non-psychotomimetic by a
  mechanism OTHER than subunit selectivity.
- **Rislenemdaz (CERC-301)** -- GluN2B-selective, safety pharmacology and neurotoxicity
  CLEAN, failed only on antidepressant efficacy (irrelevant here).
- REJECTED: memantine and dextromethorphan have **PCP-like**, not ethanol-like,
  discriminative stimulus effects. Glycine-site compounds (AV-101, rapastinel) avoid
  dissociation by being essentially non-psychoactive -- disqualifying here.

### EXISTENCE PROOF
**Nitrous oxide** is an NMDA antagonist AND a GABA enhancer in one molecule, is
intoxicating, and causes **no respiratory depression or desaturation at conscious-sedation
levels**. It proves the GABA+NMDA combination can intoxicate without respiratory cost.
Disqualified as a product by irreversible B12/cobalamin inactivation causing neuropathy
(mechanism specific to N2O, not shared by the approach), hypoxia risk, and delivery.
Isoflurane is a second such proof.

### OTHER ETHANOL TARGETS SURVEYED
- **ML297** (first selective GIRK1/2 activator): GIRK IS an ethanol target; anxiolytic and
  antiepileptic **without sedation or addiction-related behaviour**. Promising adjunct arm.
- **5-HT3**: ethanol is a PAM of 5-HT3A but NOT 5-HT3AB; 5-chloroindole is an available
  PAM. Subunit-selective in principle, but 5-HT3 potentiation causes nausea/emesis.
- **Glycine receptor**: still the main undeveloped ethanol mechanism. No clean selective
  PAM exists; ivermectin is far too promiscuous; tropeines potentiate a1 GlyR but INHIBIT a3.
- **GHB / sodium oxybate**: genuinely alcohol-mimetic and clinically used for alcohol
  withdrawal in Italy and Austria, but a very narrow therapeutic index -- the overdose
  problem reproduced.
- **Orexin antagonists**: the safest profile in the whole survey (no respiratory
  abnormality to 1000 mg/kg; abuse potential below Schedule IV hypnotics) but they produce
  SLEEP, not intoxication. Instructive: safety and the wanted subjective effect are in
  genuine tension.

### REGULATORY: the question from Session 1 has an answer
GABA Labs / Alcarelle are pursuing **novel food** registration (UK/EU/Canada) and **GRAS**
(US), not a drug pathway, with **SENTIA Plus announced for a 2026 US launch**. Composition
undisclosed; claims not independently verified here.

### DELIVERABLE: `simulator.py`
CLI and library for novel substances. Takes a receptor-activity profile
(GABA-A a1/a2-3/a5/delta/epsilon efficacies, efficacy ceiling, NMDA block + GluN2B
selectivity, glycine gain) and returns predicted subjective index with the salience
constraint checked, ventilation, stretch reflex, and overdose index -- plus an explicit
list of what it CANNOT predict.

    python simulator.py subtypes
    python simulator.py profile --a5 1.0 --a23 0.15 --nmda-block 0.29 --subj 0.50
    python simulator.py compound mp_iii_022 --nmda-block 0.29
    python simulator.py refs
    python simulator.py optimise --subj 0.50

Known limitation: `optimise` is slow because the overdose scan runs serially per profile.

---

## 2026-10-07 — Session 4: the non-hepatotoxic compound

### SCAFFOLD HYPOTHESIS: hepatotoxicity tracks CHEMOTYPE, not the GABA-A mechanism
All four hepatotoxic GABA-A anxiolytics are NON-benzodiazepine chemotypes:
alpidem (imidazopyridine), ocinaplon (pyrazolopyrimidine), panadiplon (quinoxalinone),
kava (kavalactones). The benzodiazepine and IMIDAZOBENZODIAZEPINE scaffolds -- diazepam,
midazolam, zolpidem, and the whole Cook-lab subtype-selective series (MP-III-022,
QH-ii-066, HZ-166, KRM-II-81, GL-II-73) -- carry no hepatotoxicity signal. And MRK-016's
a5 programme failed on TOLERABILITY in the elderly, not on liver injury.

A hypothesis from a small sample, not established. But it means the a5-selective candidates
sit on scaffolds with clean records rather than on the chemotypes that failed.

### THE ANSWER: Alogabat (RG7816 / RO7017773)
Roche, **Phase 2** (autism spectrum disorder; Angelman syndrome). **Highly selective
GABA-A a5 positive allosteric modulator.**
- sustained **>=80% a5 receptor occupancy in humans**
- **safe and well tolerated** so far
- being dosed in **CHILDREN aged 5-17 for 12 weeks** (ALDEBARAN trial) -- chronic pediatric
  dosing is about as strong a hepatic safety signal as exists short of approval
- no hepatotoxicity signal

Simulated (subjective 0.45): **ventilation 101%, reflex 103%, overdose >26x.**
With a GluN2B-selective antagonist added (subjective 0.86): **ventilation 93%, reflex 100%.**

### Other non-hepatotoxic candidates found
| compound | status | hepatic | simulated |
|---|---|---|---|
| **Alogabat** | Phase 2, a5-selective PAM | no signal; pediatric chronic dosing | vent 101%, reflex 103% |
| **GL-II-73** | preclinical, a5-preferring, NO a1 affinity | no signal | vent 100%, reflex 101% |
| **Basmisanil** | Phase 2 a5 NAM (wrong direction) | AEs at PLACEBO level, max occupancy by PET | de-risks the target |
| **KRM-II-81 analogs** | IND-enabling (RespireRx) | human hepatocyte metabolism mapped | minimal respiratory depression (explicit) |
| **Imepitoin** | APPROVED (vet) | **no rise in ALT/AP/AST/GGT across doses** | UNREACHABLE (subj 0.028) |

### GL-II-73 is NEUROTROPHIC on chronic dosing
Increased dendritic length and spine density of pyramidal neurons in prefrontal cortex and
hippocampus in stress and aging models -- the OPPOSITE of chronic alcohol's dendritic
atrophy and white-matter loss. If it holds for the class, the a5 route offers a positive
claim rather than merely an absence of harm.

### ASP-8062 addresses the one risk nothing else did
GABA-B PAM, Phase 2 for alcohol use disorder. Over 6 weeks in AUD patients: adverse events
no different from placebo, little sedation. And critically, a dedicated co-administration
study found **safety findings with alcohol alone were NOT AUGMENTED when ASP8062 was given
with alcohol.** Users drinking on top of the product was the dominant unmitigated risk;
this is the first candidate with direct evidence against it. Candidate adjunct arm.

### IMEPITOIN: instructive boundary case
Genuinely non-hepatotoxic with explicit dose-ranging liver enzyme data, approved, no
tolerance or withdrawal. But its maximal efficacy is only **20% of diazepam**, giving a
simulated subjective index of 0.028 -- the simulator correctly flags it UNREACHABLE.
It is non-hepatotoxic PRECISELY BECAUSE it is too weak to intoxicate. Reinforces the
theme: the hard constraint in this project is the subjective requirement, not safety.

### a5-PAM EEG signature converges with ethanol
a5-PAMs decrease peak THETA power in freely moving rats, in contrast to diazepam. Ethanol
slows hippocampal theta. So the a5 route shares an EEG signature with ethanol that diazepam
lacks -- independent support arriving from the axis I had discounted (finding 1), after the
discrimination and regional-expression evidence.

### Simulator bug fixed
`reachable()` returned True whenever the dose stayed under the ceiling, even when clamping
meant the requested subjective effect was never DELIVERED. The optimiser consequently
reported a "best" profile that only reached subjective 0.281 against a 0.50 target. Now
checks delivery, and the optimiser filters on it. Caught by imepitoin, whose low ceiling
makes the failure obvious.

---

## 2026-10-07 — Session 5: chirality

### THE CENTRAL CHIRAL FINDING: one stereocentre switches GABA-A subtype selectivity
| enantiomer | profile |
|---|---|
| **SH-053-2'F-R-CH3** | **a5-SELECTIVE** -- high a5 affinity and efficacy, separation window vs a1, LOW affinity at a1/a2/a3 |
| **SH-053-2'F-S-CH3** | **a2/a3/a5** partial PAM; "completely avoids the memory impairment commonly caused by benzodiazepine site agonists" |

The two viable GABA arms of this project are the two enantiomers of the same molecule.
**MP-III-022, the lead a5 candidate, was derived from the R enantiomer** -- so the entire
a5-selective line traces back to a single carbon.

Simulated at matched subjective effect (0.35), paired contrast R minus S:
ventilation **+1.5 points**, reflex **+6.6 points**, driven by 3.6x lower preBotC
sensitivity and 3.5x lower spinal sensitivity for R.

Interpretation: in this pair the stereocentre switches between two GOOD options rather than
between good and bad, because BOTH enantiomers avoid a1. The outcome gap is therefore
modest even though the receptor-level gap is large. That is a different shape from the
classic eutomer/distomer story.

### CHIRALITY SOLVES THE NMDA ARM'S DISSOCIATION PROBLEM
In both of the best NMDA candidates, the clean option IS the enantiomer:
- **Esmethadone (d-methadone)** -- NMDA blocker, OPIOID-INACTIVE. Its mirror image
  (levomethadone) is the clinical opioid. Complete target separation.
- **Arketamine (R-ketamine)** -- lower psychotomimetic incidence than esketamine; esketamine
  produces depersonalisation and hallucinations at the same dose where arketamine does not;
  and esketamine but NOT arketamine enhances striatal dopamine release, implicated in both
  dissociation and abuse potential. **It inverts the naive prediction: the WEAKER-binding
  enantiomer is both more potent and safer.**

That gives three independent routes to non-psychotomimetic NMDA antagonism -- subunit
selectivity (GluN2B), low trapping (lanicemine), and **chirality** (arketamine,
esmethadone) -- which was previously the arm's uncorrectable liability (finding 32).

Also noted: dextromethorphan/levomethorphan is a second case where one enantiomer is an
NMDA blocker and its mirror is an opioid. Both rejected here (PCP-like), but the pattern
is striking.

### CHIRAL PAIR DATABASE (12 pairs, `kb.py chiral`)
Classified by KIND, which turns out to matter more than potency ratio:
- **SWITCH (4)** -- enantiomers hit DIFFERENT targets: SH-053 (a5 vs a2/3/a5), methadone
  (NMDA vs opioid), methorphan (NMDA vs opioid), salsolinol (R active, S not)
- **SAFETY (1)** -- same target, one enantiomer cleaner: ketamine (arketamine vs esketamine)
- **OPPOSITE (2)** -- opposing actions: allopregnanolone 3a-OH PAM vs 3b-OH antagonist;
  certain N-alkyl barbiturates, one depressant and one CONVULSANT
- **POTENCY (5)** -- same target, different potency: etomidate (10-20x; given clinically as
  the SINGLE R isomer, unlike most chiral anaesthetics), baclofen (~100x),
  levetiracetam (~1000x), zopiclone, amphetamine (3-5x)

The GABA-A modulatory site is strongly chirally discriminating -- "the chiral binding site
provides differential stabilization of the R-enantiomer" (etomidate).

### SIMULATOR: `pair` subcommand
    python simulator.py pair sh053_R sh053_S --subj 0.35

Evaluates an enantiomer pair and reports a PAIRED CONTRAST. Rationale, which is the
statistically strongest thing this model can do: both arms run through the same pipeline
with the same wrong parameters, so shared systematic error (wrong conductances, wrong
connectivity, wrong plant, uncalibrated spinal sensitivity) cancels to first order in the
contrast. Enantiomers are also a near-perfect matched control -- identical mass, logP, pKa
and polar surface area, differing only in 3D fit. Absolute columns inherit full calibration
uncertainty; the contrast column does not.

Named profiles added: `sh053_R`, `sh053_S`, `alogabat`, `gl_ii_73`, `imepitoin`.

### DESIGN CONCEPT: the self-protecting racemate
If a stereocentre can switch subtype selectivity (SH-053 proves it can), then a molecule
whose **(R) enantiomer is an a5 PAM and whose (S) enantiomer is an a1 ANTAGONIST** would be
self-protecting as a racemate: the distomer would block the sedation/respiratory subtype
while the eutomer delivers the subjective effect. Precedent that one molecule can do both
exists -- L-838417 is a partial agonist at a2/a3/a5 AND an antagonist at a1. Speculative as
a chiral design, but not implausible, and it would convert the distomer from dead weight
into an active safety feature.

### OPEN
- **Alogabat's stereochemistry is not publicly documented** in anything I could find.
  If it is racemic there may be a chiral-switch opportunity; if a single enantiomer, the
  distomer is an unexplored control compound. Needs PubChem/patent lookup.
- The **parity-invariance unit test** from the project's first mirror-compound discussion
  still applies to any structure-based chiral work here: Mirror(ligand + receptor) is an
  exact identity, so a docking engine MUST score mirror-image complexes identically.
  Build the D-receptor by reflecting coordinates and dock the L-ligand into it; any score
  difference means the engine is broken. A free, exact validation harness.

---

## 2026-10-07 — Session 6: closing the open chiral items

### OPEN ITEM 1 RESOLVED — alogabat is ACHIRAL
Structure found: C21H23N5O4, CAS 2230009-48-8, PubChem CID 134588268, MW 409.45.
SMILES `CC1=NC=C(C=C1)C2=NOC(=C2COC3=NN=C(C=C3)C(=O)NC4CCOCC4)C`.
RDKit: **zero tetrahedral stereocentres.**

So there is no chiral switch to make and no distomer control compound. More useful than a
null result: **a5 selectivity does NOT require a stereocentre.** SH-053 achieves it with
chirality (R = a5, S = a2/3/a5); alogabat achieves it with an achiral scaffold
(isoxazole + pyridazine + tetrahydropyran amide). Two independent routes to the same
subtype selectivity, and the achiral route avoids all chiral-switch manufacturing and
regulatory complexity.

### OPEN ITEM 2 BUILT — the parity-invariance harness (`chirality.py`)
`python chirality.py parity` — EXACT pass on ketamine, levetiracetam and allopregnanolone:

| test | result |
|---|---|
| internal geometry, max abs(delta d) | **0.00e+00** |
| MMFF energy in vacuum, abs(delta E) | **0.00e+00 kcal/mol** |
| CIP labels | all flip (allopregnanolone: all EIGHT simultaneously) |

The MMFF energy being EXACTLY identical confirms the force field is parity-even, as it must
be. Usage as a docking harness: reflect BOTH receptor and ligand coordinates, re-score, and
the score must match to numerical precision. Any difference means the engine, the
force-field parameters, or the setup is broken. Costs nothing, and it is an identity rather
than a benchmark.

### Also demonstrated numerically: why 2D QSAR cannot rank enantiomers
`python chirality.py descriptors`. For levetiracetam/R-etiracetam (~1000x activity
difference) and R-/S-baclofen (~100x): MolWt, MolLogP, TPSA, HBD, HBA, RotB, RingCount and
HeavyAtomCount are **identical to four decimal places**, and the default Morgan fingerprint
is **bit-identical**. Chirality-aware fingerprints do differ but encode only the LABEL,
carrying no information about which label a chiral site prefers. This was asserted early in
the project; now proven.

### SELF-PROTECTING RACEMATE — works mechanically, REJECTED
Modelled R = a5 PAM, S = a1 NAM (negative efficacy):

| variant | preBotC sens | ventilation | subjective |
|---|---|---|---|
| eutomer alone (pure a5 PAM) | +0.0039 | 98% | 0.300 |
| racemate 1:1 with a1-NAM distomer | **-0.0273** | **109%** | 0.225 |
| stronger a1 NAM | **-0.0565** | **116%** | 0.225 |
| imperfect eutomer (residual a1 +0.25) | +0.0331 | 96% | 0.300 |
| same, RESCUED by the a1-NAM distomer | +0.0019 | 99% | 0.225 |

The distomer drives preBotC sensitivity NEGATIVE, pushing ventilation ABOVE control, and
the real use case is the last row: rescuing an imperfectly selective eutomer.

**Rejected anyway**, for three reasons: it halves subjective potency (the distomer is
subjectively inert); the stronger version made the overdose index FINITE for the first time
(26x); and decisively, **a1 inverse agonism is PROCONVULSANT** -- beta-carbolines such as
DMCM convulse through exactly this mechanism, and alcohol withdrawal already carries seizure
risk. Not worth it when alogabat and MP-III-022 are already cleanly a1-null.

### ARKETAMINE (R-ketamine / PCN-101) IS THE BEST NMDA ARM
- **Ketamine-class compounds PRESERVE and even STIMULATE breathing**: central respiratory
  drive depression "scant", upper airway/pharyngeal/laryngeal reflexes preserved,
  bronchodilator effect; and s-ketamine **stimulates breathing and ATTENUATES propofol- and
  opioid-induced hypoventilation**. The NMDA arm is respiratory-PROTECTIVE, not additive
  harm.
- **Ketamine produces dose-related ETHANOL-LIKE effects IN HUMANS** (recently detoxified
  alcoholics) -- direct human evidence for the NMDA arm's subjective contribution, not just
  rodent discrimination. PCP and ketamine dose-dependently substitute for ethanol.
- **"The simultaneous action of GABA-A agonist and NMDA antagonist mechanisms produce a
  GREATER ETHANOL-SPECIFIC discriminative stimulus than activation of either component
  individually"** -- direct support for the two-arm design being more ethanol-specific than
  either arm alone.
- NMDA's contribution to ethanol's subjective effect **grows with ethanol dose**, so the
  optimal GABA:NMDA ratio should shift toward NMDA for stronger intoxication targets.
- Clinical: atai/Perception PCN-101. Phase 1 SAD in 58 healthy volunteers, safe and well
  tolerated; IV-to-SC bridging at 60/90/120 mg. **Phase 2a missed its primary endpoint** --
  irrelevant here, as with rislenemdaz, because we want a subjective effect not an
  antidepressant one. FDA cleared a DDI study, so development continues.
- Caveat: ketamine enantiomers "share overlapping but NOT isomorphic discriminative stimulus
  effects", so arketamine's ethanol-likeness specifically is UNTESTED.

### MODEL LIMITATION DECLARED (could not be fixed)
The preBotC model **cannot represent the clinical fact that NMDA antagonists spare
ventilation.** It predicted 93% for a GluN2B arm and 70% for a non-selective one, where
clinically the answer is ~100% or better.

Two fix attempts, both failed:
1. **Fudge** -- a `resp_drive_boost` term. Too small to matter (+3% against a -30%
   prediction), and tuning it until it matched a known answer would be circular.
2. **Structural** -- reducing the NMDA share of preBotC recurrent excitation. Swept
   0.55 -> 0.06; at 6% share a 60% block STILL gave 56% ventilation.

The cause is architectural: the long-tau (100 ms) NMDA conductance is what sustains the
burst plateau, so removing it shortens bursts regardless of its weight share. No parameter
fixes this.

**Resolution (honest):** the simulator now prints an explicit MODEL LIMITATION block
whenever `nmda_block > 0`, the fudge term is removed, and the NMDA arm's respiratory
contribution is to be taken as ~100% FROM CLINICAL DATA. The model's reliable output is the
GABA arm, run with `--nmda-block 0`.

**Net effect on the recommendation: it IMPROVES.** The a5 arm alone gives
**ventilation 100%, reflex 100%, overdose >26x** (6 seeds), and the NMDA arm contributes
~100% on clinical grounds rather than the 93% the model wrongly predicted.

---

## 2026-10-07 — Session 7: the NMDA/ventilation limitation, properly tested

Goal restated by the user: **make the model as honest as possible; "there is no answer"
is an acceptable outcome.** This session produced exactly that, and in the process found
that the previous session's headline limitation had been declared on insufficient grounds.

### CORRECTION — the Session 6 limitation was UNSUPPORTED as argued

Session 6 declared: "the preBotC cannot represent NMDA antagonists sparing ventilation,
and no parameter fixes it." The evidence offered was a sweep of `ee_nmda` (the NMDA share
of recurrent excitation) from 0.55 down to 0.06.

That evidence did not support the claim. Inspection of `PreBotC.W` found the weight set to
be `ee_ampa, ee_nmda, ei_ampa, ie_gaba, ie_gly, eo_ampa, eo_nmda` — **there was no
`ei_nmda`**. The projection from the excitatory population onto the inhibitory population
was pure AMPA, so NMDA block could not reduce inhibition *by any amount*. The
disinhibition arm of ketamine's action — the mechanism most often invoked to explain why
NMDA blockers spare or stimulate breathing — was **structurally absent from the circuit**,
not mis-parameterised. Sweeping `ee_nmda` could never have found it.

The conclusion may still be right, but as argued it was an artefact of a missing pathway.

### Pathway added (verified locally)
`circuitpharm/resp.py`: `ei_nmda` added to `self.W`, injected onto `("Inh","nmda")`. Default
**0.0**, which exactly reproduces prior behaviour so earlier runs stay reproducible.

### Test 1 — does disinhibition change the sign? NO (model-derived)
`scripts/calib_ei_nmda.py`. `ei_nmda/ei_ampa` share swept 0 → 2.0 with the TOTAL exc→inh
weight held constant (so the drug-free rhythm is unchanged and each column is comparable
to its own control). Non-selective block, the harshest case. 3 seeds.

| ei_nmda share | ctrl mod | ctrl Hz | blk 30% | blk 60% | blk 90% |
|---|---|---|---|---|---|
| 0.00 (as declared) | 4.96 | 1.33 | 68% | 44% | 47% |
| 0.20 | 4.84 | 1.22 | 64% | 46% | 40% |
| 0.50 | 4.80 | 1.22 | 63% | 43% | 39% |
| 1.00 | 4.73 | 1.22 | 62% | 43% | 39% |
| 2.00 | 4.68 | 1.22 | 63% | 44% | 39% |

Adding the missing mechanism changes nothing, and at large shares is marginally **worse**.

### Test 2 — the limit case, and a REFUTED hypothesis (model-derived)
`scripts/diag_inhib_load.py`. Hypothesis: inhibition is not load-bearing for burst
amplitude, so disinhibition has nothing to release. Tested per **E10** — remove the
coupling completely, never partially. 4 seeds.

| condition | vent % | mod | freq Hz | alive |
|---|---|---|---|---|
| intact (reference) | 100% | 4.92 | 1.33 | 1.00 |
| inhibition REMOVED | **167%** | 4.12 | 0.71 | 1.00 |
| adaptation REMOVED | 689% | 0.00 | 0.00 | 0.00 |
| NMDA 90% block | 45% | 2.85 | 4.06 | 1.00 |
| NMDA 90% + inhib REMOVED | 55% | 7.74 | 1.63 | 1.00 |

**The hypothesis was wrong.** Inhibition IS load-bearing — removing it is worth +67 points
drug-free, holding ventilation to ~60% of its uninhibited value. So there genuinely is
inhibition available for NMDA block to release.

**The conclusion survives for a stronger reason.** Removing *100% of synaptic inhibition*
— vastly more disinhibition than any drug could produce — recovers only **45% → 55%**,
~10 points of a ~55 point deficit. The same removal is worth 67 points drug-free. The
asymmetry is the finding: once the long-tau NMDA conductance is gone, the rhythm's
**character** is degraded, and disinhibition can only amplify what remains, not restore
synchronous bursting. `adaptation REMOVED` (689%, mod 0.00, not alive) confirms adaptation
terminates bursts, so inhibition and adaptation are doing different jobs and only the
latter is in the NMDA pathway's way.

### Status of the limitation: UPGRADED from asserted to ESTABLISHED, and BOUNDED
It now rests on the mechanism being present and tested rather than absent, and it carries a
quantitative headroom: **disinhibition buys at most ~10 percentage points.** The model
cannot reproduce NMDA-sparing of ventilation, and the gap is far too large for the
mechanism that was missing to close. This is a genuine negative result.

### NEW DEFECT FOUND — frequency readout under severe rhythm degradation
`NMDA 90%` reports **freq 4.06 Hz** = ~240 breaths/min, not physiological for a depressed
rat. The ventilation conclusions use `mean`, not `freq`, so nothing above is affected — but
the frequency readout is untrustworthy once the rhythm fragments. This is **E4 in a form
the FFT fix does not solve**: the FFT is amplitude-invariant, which cured spurious
threshold crossings, but a fragmented burst train has genuine spectral power at high
frequency. `resp_metrics` declares `alive=True` here on `mod > 0.8` and `n >= 3`, neither
of which catches it. Do not quote preBotC frequency under heavy block until fixed.

### Metric fix — E4 second form, and it was making the overdose index OPTIMISTIC
`resp_metrics` now gates `alive` on a physiological frequency band
(`EUPNOEA_BAND = (0.30, 2.50)` Hz) and returns `frag` plus a `reason` string.

Why this matters beyond tidiness: `alive` is what the **overdose scan keys on**. A
fragmented 4.08 Hz train at 46% of control mean output was passing `mod > 0.8 and n >= 3`
and counting as **survival**, so every overdose index computed with an NMDA component was
biased optimistic. After the fix the same condition correctly reads `alive=0`, while the
control (1.33 Hz), the inhibition-removed case (0.71 Hz) and the matched NMDA+disinhibition
case (1.63 Hz) are unaffected. `scripts/tune_resp.py` had independently been applying
`0.6 < freq < 3.0` by hand, which corroborates the band.

**Calibration anchor re-verified after the change** (`scripts/calibrate_resp.py`):
`gaba_sens=0.15` still gives PAM2 ventilation **-18.4%** against the clinical -16% target
and amplitude -21.1% inside the measured -16 to -22% range. Known residual mismatches are
unchanged: PAM2 rate +25% vs +10% measured, PAM5 amplitude -38% vs -60%.

---

## 2026-10-07 — Session 7b: uncertainty propagation, and a knife-edge in the recommendation

`scripts/uncertainty.py`. The project had never asked whether the spread induced by its own
weak parameters is wider than the gaps it was reporting between compounds. It is cheap to
ask and should have been asked first. 120 parameter draws x 2 arms x 2 seeds = 480 sims.

Sampled: regional subunit fractions (Dirichlet, kappa=15, preserving nominal means and
ordering) and the preBotC anchor (lognormal, median 0.15, 90% CI ~[0.09, 0.25], matching
the midazolam clinical spread). Ceiling and circuit parameters held FIXED, so every
interval below is a **lower bound** on true spread. Arms evaluated under the SAME draw and
seeds so the paired difference cancels shared systematic error.

### FINDING 1 (the important one) — the a5 arm is UNREACHABLE at nominal parameters
Analytic, no simulation required:

| arm | gain needed | ceiling | margin | |
|---|---|---|---|---|
| a5-selective (alogabat-like, a23=0.10, a5=1.00) | 2.54 | 2.5 | **-1.5%** | UNREACHABLE |
| non-selective BZ | 1.91 | 2.5 | +23.6% | reachable |

At `subj_target=0.50` the alogabat-like profile **cannot deliver the target subjective
effect within its efficacy ceiling.** The project's own `refs` entry uses `a23=0.15`, which
needs gain 2.48 — reachable by a **+0.8% margin**. The headline recommendation therefore
sits on a knife-edge decided by the third significant figure of a literature-estimated
subunit fraction. Under parameter uncertainty the target is unreachable in **61% of draws**
for the a5 arm versus 7% for the non-selective one.

This was foreseeable and was in fact foreseen — `scripts/subtype_sweep.py` says in its
docstring that "a highly selective compound with an efficacy ceiling has a LOWER MAXIMUM
achievable subjective effect, because it carries the whole subjective load on fewer
subtypes." That consequence was never propagated into the recommendation.

**The binding constraint on this project is NOT respiratory safety. It is efficacy ceiling
versus subjective load.** Selectivity buys safety and spends reachability, and at
subj_target 0.50 the a5 arm has already spent more than it has.

### FINDING 2 — conditional on reachability, the ranking IS robust
Restricting to the 47/120 draws where BOTH arms reach the target (matched subjective
effect — comparing a full dose of one against a clamped, under-delivering dose of the other
is not a comparison, and the first version of this script made exactly that error before it
was caught):

    a5 minus non-selective, per draw:  median +11.4 pp, 90% CI [+4.9, +20.3] pp
    a5 arm better in 100% of draws  ->  interval excludes zero

So the model **can** rank these two arms on respiratory burden, and the ordering survives
its own parameter uncertainty. That vindicates the `circuitpharm/subtypes.py` claim that the
fractions' *ordering* is solid even though their values are not.

### FINDING 3 — marginal intervals, and why they mislead on their own
| arm | median vent | 5th | 95th | width | preBotC sens (median) |
|---|---|---|---|---|---|
| a5-selective | 100% | 97% | 102% | 5 pp | 0.005 |
| non-selective BZ | 84% | 74% | 93% | 19 pp | 0.150 |

The a5 interval is narrow mostly because its preBotC sensitivity is ~30x lower — it is
barely touching the network — and partly because 61% of those draws are **under-dosed**.
A good marginal number can mean "safe" or merely "not delivering the effect"; the script
now says so in its output.

### Consequence for the recommendation
The previously reported "a5 arm: ventilation 100%, reflex 100%, overdose >26x" is not
wrong, but it is **incomplete in a way that flatters it**: at the subjective target used,
that arm is at or past its ceiling. The honest statement is that the a5 arm is safer *per
unit of subjective effect delivered* and may be unable to deliver enough of it. Either the
subjective target must be lowered, or the ceiling raised (which spends the overdose
protection that motivated the PAM choice), or the load shared with a second mechanism.
**That trade is now the live question, and it is a harder one than the safety comparison.**

---

## 2026-10-07 — Session 7c: GABA-A as an explicit Markov scheme

New module `circuitpharm/gabaa_kinetics.py`. Motivation: the model carries THREE independent
hand-set numbers for a PAM (`gaba_a_gain`, `gaba_a_tau`, `gaba_a_efficacy_cap`) which are
in reality three consequences of ONE thing — how the modulator shifts gating rate
constants. Tying them to a single kinetic scheme removes two free parameters and lets the
ceiling emerge, and the ceiling is what Session 7b identified as the binding constraint.

**Scheme.** Five states, two sequential agonist binding steps with statistical factors,
one open state, desensitisation from the doubly-bound shut state (Jones & Westbrook
topology): `R <-> AR <-> A2R <-> A2O`, with `A2R <-> A2D`.

**Structure of the calibration — this is the point.** Baseline scheme fitted to three
BASELINE observables (GABA EC50 20 uM, max open probability 0.75, IPSC deactivation tau
15 ms; all three hit exactly). The drug then gets exactly ONE free parameter — an
affinity multiplier slowing agonist unbinding — calibrated against exactly ONE drug
observable, the ~2.5x leftward EC50 shift that defines a classical BZ-site ligand.
Everything else is a PREDICTION.

### Protocol error, caught by a known-phenomenology check
First version used a 1 ms agonist step for the dose-response. At low concentration the
peak is then BINDING-RATE limited, not equilibrium limited, which inflated the fitted
EC50, forced an implausibly high microscopic affinity (Kd 1.9 uM, below the EC50 — the
wrong ordering, since gating should make EC50 < Kd), and gave a baseline tonic open
probability ~0.38, which is absurd for a tonic current. Caught because the scheme then
failed to reproduce the BZ EC50 shift at all (predicted 1.03-1.08x vs measured 2-3x).
Fixed to a 300 ms application, matching the experimental protocol. After the fix:
Kd 33 uM > EC50 20 uM (correct ordering), baseline tonic Po **0.00076** (correct order of
magnitude), EC50 shift physiological.

### Scoring error I made, and corrected
The first sweep scored four predictions alike and announced a NEGATIVE RESULT. But only
three are anchored to robust measurements (mIPSC amplitude barely changes; mIPSC decay
prolonged ~1.5-2x; no effect at saturating agonist). The fourth — the *magnitude* of tonic
current potentiation — was a range I invented from the verbal description "strongly
potentiated," and it was the only one failing. Declaring a model wrong against a range
invented from a phrase is pretending to a measurement that does not exist. Rescored with
unanchored quantities REPORTED but not scored.

### Result (`scripts/kinetics_sensitivity.py`, 27 cells, baseline refit per cell)
**9/27 physiological settings reproduce every anchored BZ observable from one calibrated
parameter.** All passing cells sit at a SATURATING cleft (peak 1000-3000 uM with slow
clearance) — the physiologically expected regime, selected by the data rather than imposed.

    cleft peak 3000 uM, clearance 1.00 ms, affinity x2.79
      phasic peak gain 1.06   tau ratio 1.66   ceiling-at-saturation 1.00

The predicted tau ratio of **1.66** independently reproduces the model's hand-set 1.60 —
a number that was previously pure invention. Standing unvalidated prediction: tonic-current
potentiation of 6-9x at clinical BZ potency.

### THE CONSEQUENTIAL FINDING — the efficacy-ceiling safety argument is mis-specified
I first wrote this backwards in the module docstring and caught it numerically. The gain at
*saturating* agonist (~1.00) answers "what does the drug do to an already-saturated
receptor?" It is NOT the analogue of `gaba_a_efficacy_cap`, which is trying to bound the
gain at *physiological* agonist. Those are different quantities. Measured properly:

| ambient GABA (uM) | tonic headroom | tonic gain at BZ potency | phasic gain |
|---|---|---|---|
| 0.1 | 3295 | 7.66 | 1.06 |
| 0.4 (physiological) | **211** | 7.21 | 1.06 |
| 1.0 | 35.8 | 6.01 | 1.06 |
| 3.0 | 5.3 | 3.10 | 1.06 |
| 10.0 | 1.5 | 1.37 | 1.06 |

**The model's cap of 2.5 is only correct at an ambient GABA of ~4-5 uM, an order of
magnitude above physiological.** At realistic ambient the headroom is ~211x.

KNOWLEDGE.md states the project's central safety claim: *"the GABA arm's safety must come
entirely from the ceiling property — a low-efficacy PAM (requires endogenous GABA,
saturates), never a direct agonist. You cannot dose-shift out of it."* The kinetics say
that argument is strong for SYNAPTIC receptors (near-saturating transient, headroom ~1.1x)
and **weak for EXTRASYNAPTIC ones (headroom ~211x). alpha5 is predominantly
extrasynaptic**, so the protection is weakest exactly where the a5 strategy places the drug.

**What this does NOT overturn.** The clinical fact that BZ overdose is rarely fatal alone
is an observation and stands. A PAM still requires endogenous GABA, so at zero ambient it
does nothing regardless of headroom — it is genuinely unlike a direct agonist. What changes
is the MECHANISM the model attributes the protection to. If protection comes not from
agonist saturation but from the ligand's bounded intrinsic allosteric efficacy (maximum
achievable EC50 shift, ~2-3x for classical BZ-site ligands), then it is **a property of the
specific molecule that must be measured, not a property of the PAM class that can be
assumed.** For a project designing a novel compound that distinction is the whole point:
the ceiling cannot be inherited from the class, it has to be demonstrated per candidate.

### Status
Module and sweep built and self-consistent on anchored data; NOT yet wired into
`circuitpharm/cpg.py`. Wiring it requires splitting GABA-A into separate tonic and phasic
conductances (they now have different gains, 7.2x vs 1.06x, where the model applies one
number to both) and then RE-CALIBRATING `gaba_sens` against the midazolam anchor, since
the decomposition changes even though the current calibration is internally consistent.

### Two-pool drug wired in (verified locally)
`Drug` now carries `gaba_a_gain` (phasic) and `gaba_a_gain_tonic`, with matching caps, plus
`Drug.from_kinetics()` which derives all four from the scheme. `gaba_a_gain_tonic=None`
means "follow the phasic value", reproducing the old one-pool behaviour EXACTLY — verified:
drug-free gives 1.000/1.000 and legacy PAM2 gives 2.000/2.000. The tonic sites in
`resp.py` and `circuit.py` now call `gaba_scale_tonic()`.

Derived BZ at a 2.5x EC50 shift, ambient 0.4 uM: phasic **1.062x**, tonic **7.214x**,
tau **1.656x**, tonic cap 210.7. The model previously applied ONE number (2.0) to both
pools, i.e. it was wrong by ~7x on whichever pool it had not been calibrated against.

### RE-CALIBRATION — and a target that should never have been fitted
`scripts/recalibrate_kinetic.py`. The old `gaba_sens=0.15` is attached to a
parameterisation that no longer exists, so a re-fit was mandatory.

**First methodological finding: the respiratory-RATE target is unreachable in principle.**
The clinical +10% is a COMPENSATORY increase — midazolam cuts tidal volume, PaCO2 rises,
chemoreceptors drive rate up. This model is an isolated preBotC with fixed drive, no gas
exchange, no CO2 compartment, no chemoreceptor. It cannot produce compensatory tachypnoea
at any parameter setting. The old one-pool fit reported rate **+25%** vs target +10% and
logged it as a "known residual mismatch"; the two-pool fit gives **-7.4%**, the opposite
sign. Neither number is evidence about the drug. Including an unreachable target in the
objective pushes the fitted parameters to absorb an error they cannot remove, so it was
biasing the calibration. Rate is now reported as a diagnostic and excluded from the fit.

**Second finding: the calibration is NOT IDENTIFIED.** With rate excluded there are two
unknowns (the ligand's effective EC50 shift at the sedative dose, and `gaba_sens`) and
effectively ONE usable constraint, because the amplitude/rate SPLIT is set in vivo by the
same chemoreflex the model lacks — their product is meaningful, the components are not.
Result: **18 of 35 grid cells put ventilation inside the clinical spread (-8% to -26%)**
and are therefore equally acceptable:

    acceptable ec50_shift  1.5 - 4.0
    acceptable gaba_sens   0.060 - 0.300   (a FIVE-FOLD range)

They trade off along a diagonal valley — a weaker ligand in a more sensitive network is
indistinguishable from the reverse. The nominal best cell (shift 3.0, `gaba_sens` 0.100,
ventilation -19.9%, amplitude -12.5%) is one point in that valley, not a determination.

**This retrospectively weakens the project's single real calibration.** `gaba_sens=0.15`
has been carried everywhere as "the fitted value anchored on human midazolam data". It is
one point in a 5-fold-wide valley, and that uncertainty was never propagated into anything.
Session 7b sampled the anchor over a 90% CI of roughly [0.09, 0.25]; the true identified
range is wider still.

**What would fix it:** an anchor from a preparation that ALSO lacks chemoreflex feedback —
benzodiazepine or GABA dose-response on burst frequency and amplitude in an isolated
preBotC slice or en-bloc brainstem. Structure-matched, so both components become
comparable and the two parameters separate. This is the highest-value single measurement
identified so far, and it is in-vitro electrophysiology, not simulation.

### Session 7 summary — what improved and what got worse
IMPROVED: the NMDA limitation is established and bounded rather than asserted; the overdose
index is no longer inflated by fragmented rhythms counting as survival; the GABA-A arm runs
on one mechanistic parameter with four falsifiable outputs instead of three invented
numbers, and independently reproduces the hand-set tau ratio (1.66 vs 1.60).

GOT WORSE, in the sense of being honestly worse than believed: the a5 arm cannot reach its
subjective target at nominal parameters; the efficacy-ceiling safety argument is weak
exactly where a5 acts; the sole calibration anchor is 5-fold underdetermined; and one of
its three fit targets was never reachable.

Net: the model is more accurate and considerably less reassuring. Those are the same thing.

### Compartment split implemented (verified locally) — STOPPED MID-TASK HERE
`circuitpharm/subtypes.py`: added `EXTRASYN` (per-subtype extrasynaptic fraction) plus
`regional_sens_split()` and `subjective_index_split()`. Split identities verified exactly —
tonic + phasic reproduces the old total for every profile, so the per-region calibration is
untouched. `circuitpharm/resp.py`: `PreBotC` now accepts `gaba_sens_tonic` / `gaba_sens_phasic`,
both defaulting to `gaba_sens`; legacy identity verified exactly (mean 22.7744 both ways).

**The finding this enables — the project's two core assumptions are the same parameter and
cannot both hold.** With compartment-correct gains (tonic 7.21x, phasic 1.062x at one
clinical-potency ligand):

| arm | subj tonic | subj phasic | delivered | vs target 0.50 |
|---|---|---|---|---|
| alogabat (a5) | 0.245 | 0.080 | **1.527** | reached |
| mp_iii_022 | 0.248 | 0.090 | 1.544 | reached |
| ideal_a5 | 0.240 | 0.060 | 1.495 | reached |
| non-selective BZ | 0.290 | 0.260 | 1.818 | reached |

So the Session 7b reachability wall was an **artefact of applying a synaptic-flavoured
ceiling (2.5) to an extrasynaptically-acting drug.** a5 is ~80% extrasynaptic; its
subjective effect is a TONIC-current phenomenon with ~211x headroom, not 2.5x.

Corollary, and it is not optional: if the subjective target is that easy to reach, **the
efficacy ceiling is not providing overdose protection.** The two assumptions — (a) tight
ceiling as the sole GABA-arm safety mechanism, (b) full subjective delivery by an
a5-selective PAM — are the same number. Make it tight enough for (a) and (b) fails; loose
enough for (b) and (a) fails.

Good news that survives: the a5 selectivity advantage is INTACT and larger in the split
model — preBotC tonic sens 0.0037 (a5) vs 0.0265 (non-selective), a 7x advantage, and
0.0031 vs 0.1235 on the phasic pool, a 40x advantage.

### OPEN — next step, designed but NOT RUN
Re-test the overdose axis on the correct mechanism. The current `simulator.py` scales dose
linearly and clips at a hard cap, which has no saturation. The honest formulation:
a ligand is characterised by its MAXIMUM EC50 shift `S_max` (its intrinsic allosteric
efficacy); dose sets occupancy; effective shift = `1 + (S_max - 1) * occupancy`, which
SATURATES at `S_max` as dose grows. Then:
  * if the rhythm survives at `S_max`, the compound is genuinely ceiling-protected
    (overdose index infinite);
  * if it fails at some shift below `S_max`, there is a finite lethal multiple.
This makes overdose protection a measurable LIGAND property (`S_max`) rather than an
assumed class property — which is the session's main conclusion. Sweep `S_max` over
2.5 (classical BZ) / 4.0 / 6.0 and both arms. Needs `gaba_sens_tonic/phasic` fed from
`regional_sens_split()`, which is now wired and ready.

Also still open: `simulator.py` and `scripts/uncertainty.py` have NOT been migrated to the
split sensitivities — they still pass a single `gaba_sens`, so they run the legacy
one-pool drug. Their numbers remain internally consistent but do not yet reflect the
compartment finding. Migrate before quoting any new endpoint number.

---

## 2026-10-07 — Session 7d: the overdose mechanism, and what `gaba_sens` really is

### Overdose rebuilt on the correct mechanism (`scripts/overdose_kinetic.py`)
The old overdose scan scaled gain linearly and clipped at the invented
`gaba_a_efficacy_cap`, so it measured the clip, not the pharmacology. Correct formulation:
a modulator cannot exceed its own **intrinsic allosteric efficacy** `S_max` (maximum EC50
shift at full site occupancy). Dose sets occupancy (Langmuir, asymptote 1). At partial
occupancy the receptor population is a **mixture** of modulated and unmodulated receptors,
so each pool's conductance gain is **exactly linear in occupancy** and saturates at its
full-occupancy value. That saturation IS the ceiling — derived from the ligand rather than
asserted. (Exact for conductance gains; first-order only for decay tau, since a mixture of
two exponentials is not one exponential. Flagged in the script.)

Arms dosed to MATCHED subjective effect (0.50) = dose 1x, then escalated to full occupancy.

### First run was UNCALIBRATED — recurring error E12, made again, one entry after warning
I scaled the split sensitivities so their TOTAL equalled 0.10, the value fitted for the
ONE-POOL drug. Carrying a calibration constant across a change of parameterisation is
exactly E12, and the previous worklog entry explicitly flagged this migration as required.

Damage: the split gives the non-selective arm tonic 0.0265 / phasic 0.1235, so tonic
sensitivity fell ~5.6x below what it was fitted at. The kinetic drug acts almost entirely
on the tonic pool (7.21x vs 1.062x), so the drug effect collapsed — the non-selective arm
read **-2% ventilation** where the anchor demands -16 to -19%, and consequently EVERY arm
"survived" EVERY dose up to full occupancy. A maximally reassuring result manufactured by
an uncalibrated parameter. Discarded.

### Calibrating the split hits a hard wall (`scripts/calib_split.py`)
One scalar `k` fitted against the midazolam anchor, applied identically to all arms so
selectivity ratios are preserved:

| k | tonic sens | phasic sens | ventilation |
|---|---|---|---|
| 3 | 0.0795 | 0.3705 | -5.4% |
| 6 | 0.1590 | 0.7410 | -11.8% |
| **10** | **0.2649** | **1.0000 (CLIPPED)** | **-17.4%** |
| 20 | 0.5299 | 1.0000 | -26.6% |

The fit lands at k=10, but **phasic sensitivity is clipped at 1.0** — i.e. reproducing the
anchor requires that essentially 100% of preBotC phasic GABA-A conductance be
benzodiazepine-modulatable. That directly contradicts the biological rationale for
`gaba_sens` existing at all: it is below 1 precisely BECAUSE the respiratory network
expresses delta, alpha4 and BZ-insensitive epsilon subunits. The calibration is riding a
boundary that the model's own justification forbids.

### Mechanism mix does NOT rescue it — and my predicted fix was wrong
Hypothesis: affinity-only modulation understates the phasic effect; adding a gating
component (which raises maximal current) would rebalance the pools. I predicted gating-only
would give phasic ~2.2x, matching the one-pool 2.0. **Wrong.** Calibrated to the SAME 2.5x
EC50 shift throughout:

| mix (0=affinity, 1=gating) | tonic | phasic | tonic:phasic |
|---|---|---|---|
| 0.00 | 7.21 | 1.06 | 6.8 |
| 0.50 | 6.41 | 1.25 | 5.1 |
| 1.00 | 5.10 | **1.40** | 3.6 |

The earlier figure of ~2.2x came from reading a gating=5.0 row whose EC50 shift was 2.46 —
but at matched shift the phasic gain is only 1.40x. The imbalance improves from 6.8:1 to
3.6:1 and no further. **Under NO mechanism mix does the kinetic scheme reproduce a uniform
2.0x on both pools.** The one-pool model's central drug parameterisation is not derivable
from receptor kinetics at all.

### WHAT `gaba_sens` ACTUALLY IS — the structural finding
Its docstring claims it is "the fraction of preBotC GABA-A conductance that is
drug-modulatable", justified by subunit composition. That cannot be right. Taken together:

* the kinetics say a BZ potentiates tonic current several-fold and synaptic current by
  1.06-1.40x, under any mechanism;
* the preBotC tonic pool is a minority (~18%) of its GABA-A conductance under the EXTRASYN
  partition;
* so the achievable preBotC rhythm depression is far smaller than the measured -16% whole-
  body ventilation drop, unless sensitivity is pushed past its own physiological bound.

The most likely resolution is that **benzodiazepine respiratory depression is substantially
NOT preBotC rhythm depression.** A large part is blunting of the hypercapnic ventilatory
response (chemoreflex), plus upper-airway dilator tone and cortical/behavioural drive —
none of which exist in this model. `gaba_sens=0.15` was therefore never a subunit fraction;
it is a **lumped fudge factor absorbing several extra-preBotC mechanisms into one in-model
parameter.** It fits, but it does not mean what it says, and it cannot be extrapolated to a
novel compound with a different subtype profile — which is exactly what the project has
been using it for.

This is the same conclusion as session 7c's identifiability finding but for a stronger
reason: it is not merely that whole-body ventilation under-determines two parameters, it is
that whole-body ventilation **is the wrong observable for this model entirely.**

### Consequences
1. The respiratory axis cannot be honestly calibrated against human whole-body ventilation.
   It needs an observable generated by the structure the model actually contains: BZ or GABA
   dose-response on burst frequency/amplitude in an **isolated preBotC slice or en-bloc
   brainstem**. Now the single highest-value measurement in the project, for two independent
   reasons.
2. Every absolute ventilation percentage in the project inherits this. The PAIRED CONTRASTS
   between arms are less affected, because the lumped factor is shared — which is again why
   contrasts were always the defensible output.
3. The a5 selectivity ADVANTAGE survives all of this and is larger in the split model:
   preBotC tonic 0.0370 vs 0.2649 (7.2x) and phasic 0.0312 vs 1.0000 (32x) at k=10.
   Selectivity is a ratio of sensitivities, so the shared lumped factor cancels in it.
4. The overdose verdict is NOT yet established. The calibrated re-run has not been done,
   and doing it at k=10 would rest on a clipped boundary. Overdose protection remains
   untested on the correct mechanism.

### Status
`overdose_kinetic.py` is correct in FORM but its numbers are void pending a defensible
split calibration, which the anchor cannot currently supply. Do not quote them.

### WHAT SURVIVES — the calibration-independent selectivity ratio
Since `sens(arm) = k * raw(arm)` with a shared unknown lumped `k`, any RATIO between arms
is independent of `k` and therefore immune to the session-7d calibration failure.
Respiratory burden weighted 6:1 tonic:phasic (the drug acts 3.6-6.8x more on tonic under
every mechanism mix tested; the weighting is itself a ratio, so `k` still cancels):

| arm | tonic burden | phasic burden | subj (tonic) | subj per unit burden |
|---|---|---|---|---|
| non-selective BZ | 1.000 | 1.000 | 1.000 | 1.0 |
| neurosteroid | 2.338 | 1.050 | 1.000 | **0.6** |
| HZ-166 (a2/a3) | 0.221 | 0.189 | 0.172 | 1.2 |
| MP-III-022 (a5) | 0.151 | 0.035 | 0.853 | 7.9 |
| alogabat (a5) | 0.140 | 0.025 | 0.845 | **8.6** |
| IDEAL a5-only | 0.118 | 0.006 | 0.828 | 10.9 |

All relative to a non-selective benzodiazepine. The a5 arms retain roughly an
order-of-magnitude advantage in subjective drive per unit of respiratory burden, and the
neurosteroid arm is WORSE than a non-selective BZ (0.6x) — consistent with the earlier
retraction of the neurosteroid recommendation, now on independent grounds.

**This is the only safety-relevant quantity in the project that session 7d does not
invalidate.** It ranks compounds; it does not bound anything. It yields no absolute
ventilation percentage, no overdose multiple, and no claim that any dose is safe. The
project's defensible output has narrowed to a RANKING.

---

## 2026-10-07 — Session 7e: the in-vitro anchor exists, and it contradicts the model

Literature search for the structure-matched respiratory anchor identified in 7c/7d (a
GABA-A modulator concentration-response in a preparation with NO chemoreflex). Three
findings, one of which overturns a conclusion from 7c.

### 1. CORRECTION to session 7c — "rate increase unreachable in principle" was TOO STRONG
7c concluded that the clinical +10% respiratory-rate rise under midazolam is necessarily
chemoreflex-mediated (reduced tidal volume -> PaCO2 up -> chemoreceptor drive up), and that
since this model has no CO2 loop the target is unreachable in PRINCIPLE and must be
excluded from the objective.

**Literature (corroborated across two independent search summaries):** acute **diazepam
1 uM applied to an isolated medulla-spinal cord preparation INCREASED respiratory-like
frequency** in newborn rats. That preparation has no chemoreflex, no gas exchange and no
higher CNS. So an INTRINSIC network mechanism for a frequency increase exists, and the
clinical rate rise need not be purely chemoreflex compensation.

Revised position: the model's frequency prediction (two-pool drug gave **-7.4%**) has the
**WRONG SIGN** against the only structure-matched measurement found. That is a model
DEFECT, not an unreachable target. Excluding rate from the objective is still defensible on
other grounds (the clinical number mixes mechanisms the model lacks), but the stated
justification was wrong and the exclusion must not be used to dismiss the sign error.

**CAVEAT that keeps this from being usable as an adult anchor yet.** The preparation is
NEONATAL. In neonates GABA-A signalling can be depolarising/excitatory because of high
NKCC1 / low KCC2 and the resulting high intracellular chloride. A GABA-A PAM speeding up a
neonatal rhythm may therefore be a developmental chloride-gradient effect that does NOT
generalise to adult, where the project's clinical anchors sit. Supporting hint: the same
searches surfaced neonatal-rat KCC2/NKCC1 manipulation papers (CLP290, bumetanide) in the
context of sedative action. **Resolve the developmental chloride question before adopting
this as the calibration anchor.** Do not simply refit to the neonatal sign.

### 2. INDEPENDENT VALIDATION of the kinetic scheme's phasic predictions (literature)
Bath-applied **midazolam 1 uM prolonged the decay phase of evoked and miniature IPSCs
(GABA-A mediated) WITHOUT A CHANGE IN AMPLITUDE**; GABA/muscimol currents were enhanced at
0.1 uM, flumazenil-sensitive.

This is exactly what `circuitpharm/gabaa_kinetics.py` predicted from ONE calibrated parameter:
phasic peak gain **1.06x** (no amplitude change) with decay tau ratio **1.66x** (decay
prolonged). The two anchored checks in `scripts/kinetics_sensitivity.py` were scored against
ranges; they now have a specific matching measurement behind them. The kinetic scheme's
phasic arm is corroborated.

### 3. NEW FIDELITY LIMITATION — benzodiazepine potentiation is NON-MONOTONIC
Midazolam potentiation of GABA-activated Cl- current is **bell-shaped**: potentiation from
1 nM with a **maximum near 0.1 uM**, declining through 10 uM, and **antagonistic above
10 uM**.

The kinetic scheme models modulation as a monotonic affinity (or gating) shift, so it cannot
produce a bell. More seriously, this breaks an assumption built into the whole 7d overdose
formulation: that raising occupancy raises effect monotonically to a saturating ceiling. For
midazolam the real concentration-response TURNS OVER. In reality that is protective; in the
model it is simply absent. Any overdose statement derived from monotonic occupancy scaling
is therefore not merely uncalibrated (7d) but structurally incomplete.

### Net effect on the project's position
* the kinetic scheme's PHASIC predictions are now externally corroborated;
* its FREQUENCY prediction has the wrong sign against the one matched measurement found;
* the anchor needed to fix the calibration exists in the literature but the accessible
  instance is neonatal and confounded by developmental chloride gradients;
* the overdose formulation needs a non-monotonic concentration-response to be honest.

The defensible output remains the calibration-independent SELECTIVITY RANKING from 7d.
Nothing found here affects it, because it is a ratio.

### CORRECTION to the session-7e correction — the model's frequency sign is probably RIGHT
I over-corrected on a single neonatal datapoint. The developmental confound I flagged as
needing resolution is now resolved, and it goes against the 7e reading.

**Developmental chloride switch, measured in the preBotC specifically** (Liu & Wong-Riley,
postnatal NKCC1/KCC2 immunoreactivity across eight brainstem respiratory nuclei):
NKCC1 stays high through P0-P11 then **falls precipitously at P12**; KCC2 rises from P0 and
peaks ~P12; the two curves **intersect at or close to P11**, with KCC2 achieving clear
dominance after P12 in all eight nuclei. The authors call this "the definitive starting
point when the action of KCC2 dominates over that of NKCC1."

The en-bloc medulla-spinal cord preparation used in the diazepam study is NEWBORN, i.e.
**pre-switch**, where high intracellular chloride makes GABA-A weakly inhibitory or
depolarising. A GABA-A PAM speeding up that rhythm is therefore the expected developmental
artefact, not evidence about mature respiratory pharmacology.

**Mature-direction evidence agrees with the model.** In rhythmic preBotC preparations the
GABA-A agonist muscimol SUPPRESSES rhythm — bath 5 uM fully or partially suppressed
spontaneous rhythmic activity in all cultures tested (N=5; ANOVA F=12.76, p=0.017), and
bicuculline 10 uM completely reactivated it (N=7, p=0.0004). Bilateral muscimol into
preBotC halts the rhythm; unilateral does not. In preterm rabbit, muscimol decreased
respiratory frequency and minute volume with increased expiratory and respiratory time.

So GABA-A activation slows/stops the mature rhythm, which is the **same sign the model
produces** (two-pool drug: frequency -7.4%). The 7e claim that the model's sign is a defect
is **withdrawn**. What stands from 7e: excluding rate from the objective remains right, but
for the narrower reason that the CLINICAL +10% mixes in mechanisms the model lacks
(chemoreflex, airway, cortical drive) — not because a rate change is unreachable in
principle, and not because the model's sign is wrong.

Two corrections in two sessions on this point. Recording the pattern: **I twice drew a
sign conclusion from a single preparation without first checking its developmental stage.**
For any future respiratory anchor, establish the animal age relative to P12 BEFORE using it.

### STILL NO USABLE QUANTITATIVE ANCHOR — and now a precise specification of what is needed
What the literature offers is single-concentration results, not graded curves: muscimol
5 uM (organotypic, cultures from P2.5-P10.5 donors, 14-35 days in vitro — age relative to
the functional switch AMBIGUOUS), muscimol 10 uM (arterially perfused), muscimol 100 uM
bilateral microinjection in vivo. The eNeuro organotypic paper states explicitly that no
dose-response relationship was established; only single concentrations were tested.

**METHODOLOGICAL ADVANCE — use a DIRECT AGONIST, not a benzodiazepine.** `gaba_sens` is
defined as the network's sensitivity to a change in GABA-A conductance. Calibrating it with
a PAM requires simultaneously knowing what the PAM does to the conductance (tonic/phasic
split, intrinsic efficacy, and the bell-shaped concentration-response from 7e) — which is
exactly the entanglement that made the 7d calibration unidentifiable. A direct agonist
bypasses all of it and separates the two unknowns:

    muscimol/GABA concentration-response  ->  how sensitive is the NETWORK to GABA-A conductance
    circuitpharm/gabaa_kinetics.py              ->  what does THIS PAM do to GABA-A conductance

**Specification of the needed experiment:** muscimol (or GABA) concentration-response on
inspiratory burst frequency AND amplitude in a rhythmic preBotC slice or perfused
preparation, spanning partial to complete suppression (roughly 0.1-10 uM given the 5 uM
partial/full suppression datapoint), in animals **older than P12** so KCC2 dominance is
established. That single curve would resolve the 7d identifiability failure.

### Independent corroboration of the gaba_sens < 1 rationale (literature)
"Increased GABA-A receptor epsilon-subunit expression on ventral respiratory column
neurons protects breathing during pregnancy" — epsilon-containing receptors in the ventral
respiratory column are protective of breathing. The model's `eps=0.08` in `prebotc` and the
`subtypes.py` rationale that BZ-insensitive epsilon/delta/alpha4 subunits are WHY
`gaba_sens` is below 1 are supported by an independent line of evidence. The structural
argument is sound even though its fitted magnitude is not identified.

### ZERO-FREE-PARAMETER PREDICTION generated (`scripts/predict_muscimol.py`)
The direct-agonist insight is stronger than first stated: for an orthosteric agonist
`gaba_sens` is not merely better constrained, it is **fixed at 1.0**. `gaba_sens` is below 1
only because BZ-site ligands need gamma2 and are inactive at delta/alpha4/epsilon
receptors; GABA and muscimol bind the orthosteric site that every GABA-A receptor has. So a
bath agonist reaches the whole population and the model has **no free parameter left**.

Bath muscimol modelled as a standing tonic conductance (not a phasic PAM — a bath agonist
does not prolong synaptic decay), baseline 1.5 nS, sens 1.0, 4 seeds:

| g_tonic x | ventilation | amplitude | frequency | rhythm |
|---|---|---|---|---|
| 1.0 | 100% | 100% | 100% | intact |
| 2.0 | 87% | 95% | 90% | intact |
| 3.0 | 66% | 88% | 77% | intact |
| 4.0 | 53% | 82% | 69% | intact |
| **6.0** | **1%** | 0% | 6% | **FAILED** |
| 8.0+ | 0% | 0% | 0% | FAILED |

**PREDICTION: the curve is extremely steep.** Half-suppression at x4.0, total failure by
x6.0 — the whole graded range spans about **1.5x in conductance**, and the model falls from
87% to 1% over a 3x range. That is a very high effective Hill slope, far steeper than the
decade-wide curve a Hill coefficient near 1 would give.

**This is the falsifiable test, and the comparison must be NORMALISED.** The experiment
yields micromolar muscimol, the model yields nS, and the link between them is a receptor
density the model does not contain. So compare the parameter-free ratio
(concentration causing failure) / (concentration causing 50% suppression). The model says
**~1.5**. A real curve spanning a decade would falsify the circuit as over-sensitive to
tonic inhibition — and would mean `gaba_sens` has been absorbing that error all along, which
would explain why it had to sit so low (0.15) to match clinical data.

**Weak corroboration of the steepness, noted honestly as weak.** The one literature
datapoint — muscimol 5 uM gave "fully OR PARTIALLY suppressed" activity across N=5 cultures
— is a MIXTURE of outcomes at a single concentration. Cultures straddling a sharp
bifurcation would produce exactly that mixture, whereas a shallow curve would give uniform
partial suppression in all five. Suggestive of steepness, nowhere near sufficient, and from
cultures of ambiguous functional age.

**If the model is too steep, the identified culprit is structural, not parametric:** a
single lumped tonic conductance. Real tonic inhibition is distributed across subtypes with
different GABA affinities (delta-containing receptors have notably high agonist affinity),
and that heterogeneity would broaden the response. Fixing it means splitting the tonic pool
by subtype affinity, which the `EXTRASYN`/`REGIONS` machinery can already express.

### State of the respiratory axis after session 7
Calibration against human ventilation: **abandoned as structurally invalid** (wrong
observable for an isolated preBotC). Calibration against a mature muscimol curve:
**specified, zero free parameters, prediction published above, data not yet in hand.**
Defensible output meanwhile: the calibration-independent selectivity ranking (7d).

---

## 2026-10-07 — Session 7f: ranking robustness, and the project verdict

### The decisive test (`scripts/ranking_robustness.py`)
After 7a-7e the only surviving output was the calibration-independent selectivity ranking.
But it rests on invented numbers — `REGIONS` fractions, the `EXTRASYN` localisation
fractions added in 7c (with `eps` flagged in the source as an outright guess), and the
tonic:phasic gain ratio (3.6-6.8 depending on an unknown mechanism mix). If the ranking
does not survive those, nothing does.

20,000 analytic draws. `REGIONS` Dirichlet(kappa=15); `EXTRASYN` Beta(conc=8) per subtype
except `d_a4` fixed at 1.0 (exclusively extrasynaptic, solid) and **`eps` sampled UNIFORM
on [0,1]** — if it was declared a guess it should be treated as one; tonic:phasic ratio
uniform [3.0, 7.5]. Score = subjective drive per unit preBotC burden, relative to a
non-selective benzodiazepine, so the lumped calibration scale cancels in every draw.

| arm | median | 5th | 95th | P(better than ref) |
|---|---|---|---|---|
| IDEAL a5-only | 40.88 | **2.58** | 1218.65 | 99.9% |
| alogabat (a5) | 16.05 | **2.51** | 84.85 | 99.9% |
| MP-III-022 (a5) | 12.59 | **2.49** | 60.98 | 99.9% |
| HZ-166 (a2/a3) | 1.34 | 0.34 | 6.44 | 63.8% |
| neurosteroid | 0.60 | 0.33 | 0.84 | **0.0%** |

a5 vs a2/a3 pairwise: median +14.56, 90% CI [-0.09, +81.19], a5 better in **94.8%** of
draws. Neurosteroid worse than a plain benzodiazepine in **100%** of draws, independently
reproducing the earlier neurosteroid retraction on new grounds.

**The ranking is ROBUST.** All three a5 arms clear the reference at the 5th percentile
(>=2.49x), i.e. the ordering holds under adverse draws of every invented parameter
simultaneously. This is a real result and it is the project's surviving deliverable.

### VERDICT — none exists, and the reason is structural
Recorded in full at the top of `KNOWLEDGE.md`. Summary of the three established limits:

1. **The subjective endpoint is not reachable by this method at all.** Drug discrimination
   needs an animal that learns. The subjective index is algebra over literature weights. A
   forebrain circuit would not fix it, because its readout would still need a learning
   animal to calibrate. Category limit, not precision limit — which is also the complete
   answer to whether atomic-resolution simulation would help: it would not.
2. **The safety endpoints cannot be bounded in silico.** The respiratory calibration is
   against the wrong observable (7d); the right anchor does not exist in usable form (7e);
   the replacement prediction is untestable without wet-lab work (7e).
3. **Overdose protection is a per-molecule measured property, not a class property** (7c),
   and the project's two central assumptions are the same parameter and mutually exclusive
   (7d).

**What the method can do, proven: rank. What it cannot do: bound.** Ranking orders
candidates; it yields no margin, no overdose multiple, no dose.

### STOP CONDITION — recorded deliberately
Further elaboration of the simulator is NOT the next step and should be declined. The
dominant uncertainty now lives entirely in measurements the simulator cannot generate, so
added detail buys apparent authority without reducing it. Four things are needed and all
four are wet-lab: the >P12 muscimol curve (zero free parameters, prediction on record),
per-candidate intrinsic allosteric efficacy by patch clamp, drug discrimination in trained
animals, and the external gates (DILI, hERG/QTc, human dose-response).

Session 7 net: the model is substantially more accurate and substantially less reassuring.
Those turned out to be the same thing, every time.

---

## 2026-10-07 — Session 8: packaging, test suite, and the locomotor behavioural endpoint

Goal: make the model publishable as a generalised tool. Five commits; the repo is now
under version control (it was not), installable, and covered by 92 tests.

### Packaging (verified)
- renamed `spinal` -> **`circuitpharm`** (the package spans receptor kinetics, respiratory
  and spinal circuits, so `spinal` was a misnomer), src/ layout, `pip install -e .`
- importable from ANYWHERE; previously `import spinal.cpg` failed outside the repo root
- removed all **33** `sys.path.insert` hacks
- pyproject with optional extras (`plant`, `chem`, `dev`), MIT LICENSE, README
- `config.py`: single source for operating points that were duplicated as literals across
  five scripts. That duplication is HOW E12 happened twice.
- archived 18 superseded scripts with a README explaining each retirement; 23 remain
- `.github/workflows/tests.yml`: 3.11 and 3.12, fast and slow suites separately
- `scripts/reflex.py` was a library nine other scripts imported, with a relative model
  path that broke outside the repo root -> moved to `circuitpharm/assays.py` with proper
  path resolution (verified resolving from /tmp); the old module is a re-exporting shim

### Reliability tiers, enforced in code (`results.py`) — the distinguishing feature
`Quantity` carries tier + provenance + promotion path. Reading a **VOID** quantity RAISES,
and VOID values never appear in printed output. This exists because the project's worst
reporting error was not a modelling mistake: a configuration the model itself flagged
UNREACHABLE was summarised as a clean result, because the flag was printed NEXT TO the
number instead of blocking it. Caveats get dropped when numbers are copied.

`CALIBRATIONS` registry now records in code that `prebotc_gaba_sens` is VOID.

### TWO REAL BUGS FOUND
**1. The "validated" locomotor RG shipped with UNTUNED defaults.** The validated
parameters existed only in the stdout of a past `tune_rg2.py` run and were never written
back. `RG_E` was silent entirely — winner-take-all, not alternation — so every importer
inherited a broken rhythm while this worklog described the module as validated. Defaults
are now the validated set (drive 260, g_adapt 1.2, tau_adapt 280, ie_gly 3.0, ee_ampa
0.55) and a test guards it.

**2. `tune_rg2.py` tests the strychnine phenotype with a 95% block — error E10**, the
project's own catalogued mistake, in the script that was supposed to validate the module.
Retested by COMPLETE removal:

| condition | corr | per_F | per_E |
|---|---|---|---|
| coupled | -0.51 | 1245 | 1245 |
| 95% block | -0.50 | 1242 | 1243 (**no effect**) |
| COMPLETE removal | **-0.13** | **1154** | **1245** |

Phenotype HOLDS (rhythm intact, 20 bursts each). The original claim was sound; the tuning
script's test was not. E10 is now an executable demonstration in the suite.

Also: `rg.py` (MatsuokaRG) deleted; `RG_GAIN_PER_NEURON` documents the units change (rg2
returns Hz per neuron where Matsuoka returned dimensionless — getting this wrong is E2).

### Public API with tiered results (`evaluate.py`)
`simulator.py` had drifted in three ways that ALL flattered the model: a single
`gaba_sens` applied to both receptor pools (wrong by ~7x on one), dose escalation clipped
at a hand-set ceiling of 2.5, and an overdose index printed as a plain number with an
ignored warning beneath. Now dose is modulator OCCUPANCY, pool gains come from the Markov
scheme and are exactly linear in occupancy (a partly occupied population is a MIXTURE), so
they saturate at full occupancy — and that saturation, set by the ligand's own `s_max`, is
the only honest ceiling. `simulator.py` is a thin CLI.

The ranking now separates the central chiral pair directly: **SH-053-R 9.35 vs
SH-053-S 3.63** from one stereocentre. Sanity checks pass too: zolpidem 0.10 (a1-preferring
= all burden, no subjective drive), gaboxadol 0.00.

### NEW: the locomotor behavioural endpoint (link 6, closed for locomotion)
`circuitpharm.assays.locomotion` — free-running closed loop, NOTHING imposed:
receptor -> circuit -> motoneurons -> Hill-type muscle -> joint -> spindle feedback.

    control    excursion 1.61 rad, 0.80 Hz, alternation -0.67, duty 0.58
               both motoneuron pools 26-146 Hz (physiological)

Motor impairment now has a body-level readout instead of reflex gain alone:

| condition | step period | alternation | duty |
|---|---|---|---|
| control | 1250 ms | -0.67 | 0.58 |
| PAM 2x | 1250 ms | -0.64 | 0.49 |
| PAM 4x | **1667 ms** | **-0.34** | 0.35 |

### ...and the assay immediately caught its own INVALID metric
Joint excursion goes the WRONG WAY under sedation: control 1.607 rad -> 1.726 at PAM 2x ->
1.940 at PAM 4x. A sedative increases it. The mechanism is real, not a coding bug — more
inhibition means less antagonist co-contraction, so the joint is less stiff and swings
further — but it is an artefact of the preparation: ONE joint, body fixed, no gravitational
load, no ground contact. In an animal, lost co-contraction presents as instability and
collapse; this model has nothing to collapse against.

The key is therefore named **`excursion_rad_INVALID`** so the number cannot be used
casually, and a test pins both the name and the wrong-sign behaviour so nobody renames it
without fixing the preparation. Valid impairment metrics are PERIOD and ALTERNATION.

A valid stride/ataxia measure needs the whole body and ground reaction forces.

### Test suite: 92 tests
- `test_identities.py` (41) exact identities: the tonic/phasic split is a partition, legacy
  single-pool behaviour reproduced exactly, the GluN2B selectivity bound, per-region anchors
- `test_failure_modes.py` **E1-E12 as executable regressions**, two-sided where possible.
  Caught an error in its own E11 test (the default efficacy cap binds a gain of 3.0).
- `test_tiers.py` (12) the refusal machinery, incl. that VOID never leaks into output
- `test_evaluate.py` the public API; occupancy linearity; scale-invariance of the ranking
- `test_phenotypes.py` (11) **the only tests that check biology rather than code:**
  strychnine hyperreflexia (gain 138% of control), benzodiazepine reflex depression (77%),
  the **GluN2B selectivity window** (non-selective 84% vs selective 102% at the SAME 60%
  block), adaptation-not-inhibition rhythmogenesis, closed-loop walking, sedative
  coordination loss, and the E5 guard (Mn peak 40.6 Hz, not pinned at the 125 Hz ceiling)

### E10 fixed in `tune_rg2.py` (verified) — the sweep is now informative
The tuning script tested the phenotype with `Drug(glyr_gain=0.05)`, a 95% block, in the
very script meant to validate the module. Replaced with complete removal (`ie_gly = 0`).
Before, every row reported control and "blocked" as identical. After:

| drive | g_adapt | tau_a | ctrl per / corr | removed per / corr / bursts |
|---|---|---|---|---|
| 260 | 1.2 | 280 | 1243 / -0.51 | 1154 / **-0.13** / 20 |
| 200 | 1.2 | 280 | 1507 / -0.35 | 1321 / **-0.01** / 18 |
| 260 | 1.2 | 380 | 1907 / -0.51 | 1785 / -0.06 / 13 |

Every row now detects the coupling. The top-scoring cell is 260 / 1.2 / 280, which is
independently what was written in as the module defaults. Note `ie_gly` 3.0 and 7.0 give
near-identical results — the coupling saturates by 3.0.

### Session 8 close — the tool
Nine commits. 114 tests, parallel, 2m35s wall (was 9m32s serial before the control cache
and xdist). The repo is installable, importable from anywhere, has CI for 3.11 and 3.12,
a worked example, and 18 superseded scripts archived with reasons.

**New this session beyond packaging:**
- `results.py` — reliability tiers enforced in code; VOID quantities RAISE on access and
  never print. The distinguishing feature of the tool.
- `config.py` — single source for operating points (their duplication is how E12 happened
  twice) plus a calibration registry recording that the preBotC anchor is VOID.
- `evaluate.py` — public API; dose is modulator OCCUPANCY, pool gains derived from the
  Markov scheme and exactly linear in occupancy, saturating at full occupancy, which is
  the only honest ceiling.
- `assays.py` — behavioural protocols, with model-path resolution that works from anywhere.
- **the locomotor endpoint**: free-running closed loop, nothing imposed. Motor impairment
  finally has a body-level readout (period 1250 -> 1667 ms and coordination -0.67 -> -0.34
  under a 4x PAM) rather than reflex gain alone.
- `tests/test_chirality.py` — the parity-invariance identity pinned exactly.

**Three bugs found by doing the engineering:**
1. the "validated" locomotor RG shipped with UNTUNED defaults; one half-centre was silent
2. its validation script tested the phenotype with a 95% block (E10)
3. the new locomotor assay's joint-excursion metric reports the WRONG SIGN under sedation

Each of those looked fine from the outside, which is the recurring lesson of this project
and the reason the failure modes are now executable tests rather than prose.

### Session 8b — hygiene found three more real problems

**1. Scripts did destructive work at IMPORT time.** Eight scripts had no `__main__` guard,
so merely importing them ran them: `port_muscle.py` regenerated the MuJoCo body model and
`build_kb.py` / `build_compounds.py` rewrote the pharmacology database. Any tool that
imports or scans the package would have triggered that. All eight guarded; verified that
importing is now silent and that running them directly still works, and that
`port_muscle.py` is deterministic (re-running it produces no diff in `models/`).

**2. A stale knowledge-base fact actively recommending the deleted architecture.**
`build_kb.py` carried "Do not retry tuning the spiking RG -- use the Matsuoka RG in
circuitpharm/rg.py", pointing at a file deleted this session for being wrong. A future
agent reading the KB would have been sent straight back to the retracted approach. Fixed,
and two `model_facts` rows (the Matsuoka operating point and its period law) marked
SUPERSEDED rather than deleted, per the convention that conclusions are corrected, never
removed. Four facts added that now matter: the VOID preBotC calibration, the
pool-dependent ceiling, the tonic-vs-phasic gain asymmetry, and the locomotor endpoint
including its invalid excursion metric. `model_facts` 23 -> 27 rows.

**3. 108 lines of dead code shipping in the package.** `HalfCentreCPG` — the REJECTED
all-spiking architecture, 28% of `cpg.py` — was reachable only from three archived
scripts. Moved to `scripts/archive/half_centre_cpg.py` with a header recording why the
288-point search that rejected it was right about LIF half-centres and wrong in its
conclusion (it assumed the rhythm must come from inhibition; a group pacemaker needs no
plateau and works fine in LIF).

**Coverage is now measured and gated:** 77% -> 90% (omitting the optional-extra chirality
module) -> **95%** after the dead-code extraction, which alone took `cpg.py` from 67% to
94%. `fail_under = 90` in pyproject, and CI reports it. `cpg.py` is 253 lines, was 362.

All 114 tests still pass.

### I DESTROYED DATA, AND RECOVERED IT — and the lesson is the entry, not the recovery

While verifying the import-side-effect fix (above), I ran `scripts/build_kb.py` directly.
It contained an unconditional `if DB.exists(): DB.unlink()` and rebuilds the database from
the literal tables written into it: **42 sources and 17 findings**. The live database held
**90 sources and 58 findings**, accumulated across later sessions by other means.

So I silently discarded **48 sources and 41 findings**, then committed the damaged file in
`c86bfdb`. I noticed only because a later count came back 42/17 when I expected 90/58.

Recovered from `e1fd284` (the pre-refactor snapshot), and the model_facts work plus the new
verification columns re-applied on top. Final state: sources 90, findings 58,
model_facts 32, compounds 53. **This was recoverable ONLY because the database is tracked
in git** — which it is only because I initialised version control at the start of this
session for an unrelated reason.

**The irony is the point:** this happened during the step where I was fixing the exact class
of bug it belongs to. Finding that a script has destructive side effects on import is not
the same as making it safe to run.

**Lasting fix — a destruction guard in both builders.** Each now compares the LIVE row count
against what it is about to write and refuses unless `--force`, taking a timestamped backup
either way:

    build_kb.py       REFUSING: sources 90 live vs 42 written; findings 58 vs 17
    build_compounds.py REFUSING: compounds 53 live vs 46 written

The second guard immediately earned itself: `build_compounds.py` had the same
`DROP TABLE IF EXISTS compounds` pattern and would have destroyed 7 compounds.

`tests/test_data_safety.py` (6 tests): a canary on the accumulated row counts, both
builders must contain a guard and a `--force` override, both must have `__main__` guards,
and the citation verification column must exist with UNVERIFIED rows in it.

### Citations: verification status is now explicit rather than implied
Added a `verification` column to `sources` and `findings`, defaulting to **UNVERIFIED**.
The 90 sources carry resolvable URLs to real papers and were read at the time, but nobody
has gone back to each figure or table to confirm the attributed claim is what the paper
supports *at the stated precision*. Recording that honestly is the same discipline as the
VOID tier on quantities: omitting the column lets a reader assume a check that never
happened. Three findings that rest on the VOID calibration, the pool-dependent ceiling or
the overdose index are marked SUPERSEDED in the row itself, not only in this log.

Before publication every load-bearing quantitative citation needs checking against the
actual figure. Session 7c found my recall wrong by a factor of 1.6 on one number
(predicted gating-only phasic gain ~2.2x, actual 1.40x), so the error rate is not zero.

### Single source for the operating point (verified)
Nine live scripts still carried the preBotC operating point as literals — the exact
mechanism behind E12, committed twice. All nine now import from `circuitpharm.config`;
each verified to produce identical values, and `predict_muscimol.py` reproduces its
published curve exactly. `tests/test_config_single_source.py` pins it: no LIVE script may
contain the literal (archived ones are exempt as frozen records), the shared constants are
immutable so an in-process caller cannot poison later runs, and the VOID calibration cannot
be quietly promoted by an edit. My grep found 4 offenders; the test found 9 — five used
formatting the pattern missed.

### Pushed
Remote `willkhinz/circuitpharm` (**private**) already existed with my history through
`c86bfdb`. Fast-forward, no force, nothing overwritten. CI runs on 3.11 and 3.12.

### A bug that passed locally and would have failed CI
Simulating a lean install (hiding `mujoco`) showed `evaluate()` **crashes** without the
body plant. The guard wrapped only the IMPORT of the assay functions — but `assays`
imports mujoco lazily INSIDE each function, so on a lean install that import succeeds and
the `ImportError` surfaces on the first CALL, outside the guard, where only
`FileNotFoundError` was caught.

This is the install CI uses (`.[dev]`, no plant extra), and it passed here only because
mujoco happens to be present on this machine. Fixed by catching `ImportError` at the call
site; a lean install now returns 10 quantities including a `motor_endpoints` notice saying
what is missing and how to get it, with the respiratory and receptor-level results intact.

`tests/conftest.py` auto-skips plant-dependent tests by inspecting what each test calls,
rather than by hand-marking — so a newly added plant test cannot forget its marker and
turn a healthy lean install into a red build. The regression test simulates the absence
even when mujoco IS installed, otherwise it would pass for the wrong reason.

### Clean-install verification found 5 more failures the local run hides
The package had only ever been installed in-place in the development venv, which masks
missing dependencies and layout assumptions. Built a fresh venv in /tmp, installed
`"/Users/whinz/Biochemistry[dev]"` (LEAN: no plant, no chem) and ran the suite. Five
failures, both classes real:

1. **`test_evaluate_includes_all_three_endpoints`** asserted `reflex_gain` and `step_period`
   exist, which on a lean install they legitimately do not. The conftest auto-skip missed
   it because the heuristic looks for calls to `stretch_reflex`/`locomotion` and this test
   never names them — it calls `evaluate()` and inspects result KEYS. So the heuristic
   alone is insufficient; plant-dependent tests now carry an explicit
   `@pytest.mark.needs_plant`, which the conftest honours, with the heuristic kept as a
   safety net.
2. **Four `test_data_safety` tests** read the SOURCE of the builder scripts via
   `parents[1]/scripts`, so they assume a repository checkout. They errored rather than
   skipped for anyone who pip-installs the package and runs the tests from elsewhere. Now
   skipped behind `needs_repo`.

Result after fixing: **full install 147 passed / 0 skipped; lean install 91 passed /
18 skipped.** Both correct.

### CI split into LEAN and FULL jobs
A single job could not catch this, because whichever configuration it tested would hide
bugs in the other. Now:
- **lean** (3.11 and 3.12): core only, and it ASSERTS mujoco and rdkit are absent before
  running, so the job cannot silently become a full install and stop testing what it exists
  to test.
- **full** (3.11): installs `.[all]` with `MUJOCO_GL=osmesa` for the headless runner, runs
  coverage, and FAILS if anything is skipped — a skip on the full install means an optional
  dependency silently went missing.

Also `-n auto` instead of `-n 8` (eight workers each running a MuJoCo simulation thrash on
a 2-core runner, which is likely why early CI runs sat at 17+ minutes against 2m39s
locally), and `timeout-minutes` so a stuck runner fails instead of hanging.

`CITATION.cff` added.

---

## 2026-10-07 — Session 8c: two external code reviews, 20 defects, all real

Two reviews arrived. **All 20 findings were legitimate and all 20 were mine.** Several were
scientifically consequential rather than cosmetic. Fixed, with regressions in
`tests/test_review_regressions.py`.

### Review 1 (9 findings) — worst first
1. **Reproducibility broken.** `Pop.__post_init__` seeded its RNG with
   `abs(hash(self.name))`, and `hash(str)` is randomised per interpreter via PYTHONHASHSEED.
   Measured V0[0] = **-71.05, -62.73, -69.90** across three processes for the same
   population with the same explicit seed. Every published number was irreproducible
   between runs. Fixed with `zlib.crc32`; verified identical across processes. The
   regression test runs in SUBPROCESSES, because within one process the hash is stable and
   an in-process check would have passed while the bug was live.
2. **Global class-state mutation.** The stretch-reflex assay rescaled Ia→Mn by MUTATING
   `SpinalCircuit.W` (a class attribute) and restoring it in a `finally` — corrupting every
   other instance in the process, permanently if an exception landed first. `SpinalCircuit`
   now takes per-instance weights like `PreBotC` and `GroupPacemakerRG` always did.
   Verified results identical (0.269, 0.370) and the class attribute untouched.
3. **Motor endpoints were single-seed** (`seed=1` hard-coded) while ventilation averaged
   over `n_seed` — recurring error **E6, committed in the public API**. Now averaged.
4. **Booleans printed as 1/0** (`bool` subclasses `int`, so `f"{True:.4g}"` → "1").
5. **The NMDA subjective arm was dropped** when evaluation moved out of the old simulator,
   silently under-reporting every compound with an NMDA component — and the two-arm design
   is the entire point. Restored, but reported SEPARATELY with the total marked **VOID**,
   because the arms are in incommensurable units and the project's "GABA salience ≥ NMDA
   salience" constraint compares them directly. That constraint selects the recommended
   ratio, so the defect is load-bearing and must not hide behind a plausible total.
6. **`max(1e-6, ia_dyn - ia_base)`** divided by 1e-6 when the Ia response was NEGATIVE,
   turning a 0.5 Hz drop into a gain of **5,000,000**. Now NaN with a stated reason.
7. **NaN propagation**: an unreachable `s_max` made `calibrate_pam` return NaN → NaN
   conductances → NaN voltages, surfacing as "SVD did not converge" and raw LAPACK errors.
   Now raises at the cause.
8. **Over-broad skip heuristic**: any test with "reflex" or "walk" in its NAME was skipped
   on a lean install. A test skipped when it should run is worse than a failing one.
9. **`fail_under` in pyproject** fired on any coverage run, so `pytest -m "not slow" --cov`
   failed at 65% with every test passing. Moved to CI on the full suite.

### Review 2 (11 findings)
1. **`np.trapezoid` is NumPy 2.0+** while pyproject allows `numpy>=1.26`; `np.trapz` was
   removed in 2.0, so neither name alone spans the supported range. Resolved once at import.
2. **GLYCINE SENSITIVITY WAS BOUND TO GABA-A SUBUNIT FRACTIONS** — the most serious finding
   in either review. `sens=(gaba_sens_phasic if rec in ("gabaa", "gly") else 1.0)` scaled
   glycinergic drug action by a fraction computed from α1/α2-3/α5 expression. Glycine
   receptors contain no GABA-A subunits, so this is a category error. Measured: a 1.6x
   glycine potentiation became **1.0019x** at alogabat's preBötC sensitivity — ethanol's
   glycine mechanism and strychnine were ~99.8% deleted from the circuit, silently. Now a
   separate `glyr_sens`, defaulting to 1.0, in `resp.py`, `rg2.py`, `circuit.py` and
   `plant.py`.
3. **`from_profile` discarded each compound's ceiling**, defaulting `s_max=2.5` for all:
   imepitoin (1.25) was modelled at twice its real ceiling and the neurosteroid arm (6.0)
   at 42% of its. Since `s_max` sets achievable effect AND overdose protection, this
   flattened exactly the distinction the safety argument turns on.
4. **Trajectory fencepost**: segments sampled s = 0 … T-DT, so a ramp never reached its
   endpoint and the next segment jumped to it — a 0.018 rad step in one 2 ms timestep,
   injecting an impulse into spindle velocity, which is the dominant Ia input. Now 0.000.
5. **`mj_name2id` returns -1** for an unknown name and Python indexes the LAST actuator, so
   a typo'd joint silently drove the wrong muscle. Now raises.
6. **`reflex_panel.py` double-scaled Ia** (its own `scale` AND the assay's 0.30 default),
   so every point on that sweep was mislabelled by >3x — and it still mutated the class
   attribute.
7. **`kb.py` opened the database read-write** while accepting arbitrary SQL. Now `mode=ro`.
8. **`build_compounds.py` had no backup** and used `>` so `live == WILL_WRITE` wiped edits.
9. **A NaN selectivity ratio was tagged VALIDATED**, printing
   `selectivity_ratio: nan x [VALIDATED]` — unearned authority of exactly the kind the tier
   system exists to prevent. Non-finite now VOID.
10. **`np.int64` is not `isinstance` of `int`**, so numpy scalars lost units and precision.
11. **CLI `--s-max` defaulted to 2.5**, overwriting every profile's ceiling.

### And fixing them corrected a previously reported result
Removing the trajectory impulse (R2 #4) and averaging over seeds changed the GluN2B
selectivity window: **75% non-selective vs 100% selective, a 24-point window**, against the
18 points previously recorded from single-seed runs (84% vs 102%). The window is **wider
and cleaner** than reported. The single-seed value after the fix was 92%, which tripped my
own threshold — reading it as the answer would have been E6 again, inside the test suite
that exists to catch E6.

**160 tests passing.** The lesson I take: every one of these 20 produced plausible output.
Twelve sessions of scrutiny by the same eyes did not surface a single one of them.

### Review 3 (8 findings) — and two were caused by my review-2 fixes
All 8 legitimate. **28 defects across three reviews, all mine.** The new lesson is in the
two self-inflicted ones.

1. **Duty cycle silently NaN whenever a trace started mid-burst.** `burst_metrics` paired
   onsets and offsets BY INDEX, which is only correct when the trace starts below
   threshold. Analysis windows open after a settling period at an arbitrary phase, so
   starting high is routine: every pair shifted by one, `off[i] > on[i]` failed for every
   burst, and duty returned NaN for a strongly, regularly bursting circuit. Verified
   phase-dependence directly. **The silence was the damage** — `scripts/tune_rg2.py`
   scores a NaN duty as 0.0 and penalises it, so locomotor parameter sets were ranked
   partly on where their analysis window happened to open. Fixed by pairing each onset
   with the first offset that follows it; duty is now identical across five phases and a
   50% square wave reads 0.500.

2. **Unreachability undercounted.** `scripts/uncertainty.py` ran `if np.isnan(mn):
   continue` BEFORE `if not reached:`, so a draw whose forebrain drive collapsed never
   incremented the counter and never entered `bad_draw`. That deflated the headline
   "target unreachable in X% of draws" AND let failed draws into `reached_all`,
   contaminating the matched-subjective paired contrast — the one comparison in that
   script that is meant to be apples-to-apples.

3. **`resp_panel.py` ran at the UNCALIBRATED default** `gaba_sens=1.0`, making the network
   3-4x too drug-sensitive: it reported severe respiratory depression for doses that
   clinically produce mild sedation. Every other analysis passes 0.15 or the split value.

4. **Two forebrain scales that both sound like the same thing.**
   `regional_sens("forebrain")` includes `K_REGION` (1.00 for a non-selective BZ);
   `subjective_index()` omits it (0.55) — a 1.8x gap. The omission is CORRECT, because the
   index is a weighted SUBSET of subtypes and a whole-conductance normaliser would be
   meaningless on it. But undocumented it invites exactly the comparison it cannot support.
   Now stated at the source, including that the subjective index has no absolute scale at
   all and a "target of 0.50" is a number in invented units.

5. **SELF-INFLICTED: direct agonists crashed.** Carrying each profile's ceiling into
   `s_max` (review-2 fix #3, correct) met a hard raise on unreachable affinity shifts
   (review-1 fix #7, also correct), and together they made gaboxadol — a direct orthosteric
   agonist with `ceiling=1e9` — raise instead of evaluate. Neither fix was wrong; the
   interaction was unhandled. Now: a modality field with affinity→gating auto-fallback, and
   agonists raise with the reason and a pointer to `predict_muscimol.py`, which models them
   correctly as a standing conductance. Side effect worth noting: the neurosteroid arm now
   evaluates at its real ceiling of 6.0 and gives tonic gain **25.76** rather than being
   silently capped at 2.5.

6. **SELF-INFLICTED, incomplete: `glyr_sens` never reached the assays.** Review 2's most
   serious finding was glycine sensitivity being scaled by a GABA-A-derived fraction. I
   added `glyr_sens` to the circuits and the plant but not to `stretch_reflex` or
   `locomotion`, so the public entry points could not vary it — the fix was invisible from
   outside, which looks identical to no fix.

7. Unguarded division by a possibly-zero control mean in `calib_ei_nmda.py`.
8. `np.nanmean` over an all-NaN seed list warned "Mean of empty slice" on a LEGITIMATE
   case (a heavy sedative abolishing locomotion in every seed). Noise on correct results is
   how real warnings get ignored.

**175 tests passing.**

### What three reviews establish
Twenty-eight defects. Every one produced plausible output; none announced itself. Twelve
sessions of my own scrutiny surfaced none of them, and the third review found two that my
own second-review fixes had created or left half-done. The specific failure mode to
remember: **a correct fix can interact badly with another correct fix, and a partial fix is
indistinguishable from a complete one from the outside.** Neither is caught by re-reading
your own diff.

This is now the strongest argument in the repo for external review before publication, and
it belongs in the paper's methods rather than being quietly fixed.

### Review 4 (10 findings) — 38 defects total across four reviews
All 10 legitimate. The rate of genuine findings has not dropped across four passes.

1. **`SpinalCircuit` never forwarded ANY sensitivity to its rhythm generator.** It accepted
   `gaba_sens`, `gaba_sens_tonic/phasic` and `glyr_sens` and passed none of them to
   `GroupPacemakerRG`, which therefore always ran at **1.0** — fully drug-sensitive — while
   the pattern-formation and motoneuron layers used the calibrated spinal values. Under
   sedation the RG saw ~7.2x tonic conductance where the rest of the circuit saw ~1.5x, so
   the locomotor rhythm slowed or arrested far earlier than the tissue pharmacology implies
   and **every drug effect on step period and coordination was overstated.** Verified:
   circuit tonic 0.08 → RG 1.0.
2. **`GroupPacemakerRG` had no tonic/phasic split** — one lumped `gaba_sens` applied to both
   pools, which differ ~200x in PAM headroom. **E12 again**, reintroduced in the locomotor
   RG after being fixed everywhere else.
3. **`PreBotC` defaults diverged from `config.RESP_OP`**: drive 190 vs 170, g_adapt 1.6 vs
   2.5, tau_adapt 450 vs 400, ee_ampa **0.16 vs 0.45** — ~3x weaker recurrent excitation.
   Invisible in project results (everything passes `**RESP_OP`) but any caller writing a
   plain `PreBotC()` silently got an obsolete untuned network. Defaults now come from config.
4. **`max(1e-9, nan)` returns 1e-9**, so a NaN control turned an ordinary value into
   **150,000,000,000%**. A control reflex gain IS legitimately NaN when the Ia response does
   not exceed noise, so this path is reachable. Now a helper that propagates NaN.
5. **`nmda_scale()` returned NEGATIVE** for `nmda_block > 1` (verified: -0.05 at 1.5),
   reachable from an optimiser sweep or a sampling tail. A negative `w_scale` makes an
   excitatory conductance negative, which drives positive-feedback voltage divergence
   rather than failing visibly. Now clipped.
6. **The muscle force-sign convention was documented backwards.** MuJoCo muscle actuators
   are PULL-ONLY (force <= 0), so `gear=+1` ('_ext') yields NEGATIVE joint torque, not
   positive as the docstrings claimed. The physics was always right — nothing depended on
   the comment — but anyone reasoning about torque direction from the docs, or adding a
   joint by analogy, would have had the sign inverted.
7. **`burst_metrics` NaN guard never fired**: `np.nan <= 0` is False, so a NaN-containing
   trace slid past, `r > thresh_frac * nan` gave an all-False mask, and the function
   returned plausible empty results instead of declaring the input unusable.
8. **Locomotion ran with UNSCALED Ia weights** while `stretch_reflex` scaled them by
   `IA_SCALE=0.30` — the same circuit received **3.33x stronger afferent feedback** in one
   assay than the other. IA_SCALE is a measurement fix (keeping the motoneuron readout off
   its tref ceiling, E5), so it belongs wherever that readout is used. **This changed the
   locomotion numbers: step period 1250 → 1001 ms, alternation -0.67 → -0.79.**
9. **`subtypes.py` CLI printed un-normalised ratios** mixing the two forebrain scales — the
   exact comparison the `subjective_index` docstring warns against, printed by the module
   that contains the warning. Now normalised to a non-selective benzodiazepine.
10. **`Compound` and `evaluate` were not exported**, so `from circuitpharm import
    Compound, evaluate` raised ImportError.

### My fix for #10 was wrong TWICE before it was right
Worth recording in full, because both wrong versions looked fine.

* **Attempt 1, lazy via `from . import evaluate`** → infinite recursion. The submodule
  `circuitpharm.evaluate` and the function `evaluate` share a name, so the `from` form
  re-entered `__getattr__` on the same attribute until the stack died.
* **Attempt 2, lazy via `importlib`** → stopped the recursion but produced something
  **worse than the original ImportError**: `from circuitpharm import evaluate` bound to the
  **MODULE**, which is not callable, and which of the two you got depended on import order.
  Verified: `callable(evaluate)` was False.
* **Attempt 3** — renamed the submodule to `circuitpharm.evaluation` and imported eagerly,
  removing the collision instead of working around it. `callable(evaluate)` is now True.

A fix that replaces a loud failure with a quiet wrong answer is worse than no fix, and I
shipped that state briefly.

**175 tests passing.**

### Four reviews: what the pattern says
38 defects. Every one produced plausible output. My own scrutiny across twelve sessions
found none of them, and each review found defects created or left incomplete by the
previous round's fixes. The findings-per-review rate has not fallen: 9, 11, 8, 10.

That last number is the important one. It does not support "the code is nearly clean now";
it supports "this code has a defect density my own review cannot measure." Any claim that
this is publishable needs at least one review pass that finds nothing, and that has not
happened yet.

---

## Session 13 — roadmap recorded, substrate diagnosed, review pass 7

### The roadmap existed only in conversation
Asked what the next step was, I searched WORKLOG and KNOWLEDGE for the six-link plan and
found **nothing**. The compound-to-behaviour roadmap had never been written down. That is
the third time this project has nearly lost work to the same failure mode: tuned RG
parameters that lived only in a past run's stdout, destroyed KB rows that survived only by
accident of version control, and now the plan itself. Written into `KNOWLEDGE.md` as a
per-link status table (`cc6d3f9`).

### Link 4/5 design, written before any code
`knowledge/05-design-conductance-substrate.md`. Design only — the gate below still stands.

The roadmap's stated reason for the conductance upgrade was "LIF is not the biophysics,"
which argues from principle, and this project has been wrong from principle before. So I
**measured** what the LIF substrate does to the voltage-dependent mechanisms the
pharmacology actually depends on. `PreBotC(**RESP_OP)`, 15 s, 250k samples of `Exc`:

| quantity | LIF substrate | a spiking cell |
|---|---|---|
| V range | −71.1 .. −44.2 mV | −65 .. +20 mV |
| NMDA Mg²⁺ relief | mean **0.063** | up to 0.925 |
| Mg relief dynamic range | **4.50×** | 15.5× |
| GABA-A driving force | mean **10.4** mV | up to 95 mV |
| glutamate driving force | 44–65 mV, never collapses | collapses at burst peak |

Three consequences, in increasing order of how much they matter:

1. **NMDA conductance is held in near-permanent Mg²⁺ block** — ~6% of unblocked, always.
   `cpg.py:243` says "voltage dependence is why NMDA block is state-dependent." On this
   substrate it is not.
2. **Excitation never self-limits**: driving force never approaches zero.
3. **The phasic GABA-A pool's driving force is truncated exactly when it matters.** Both
   pools share one V, so at rest both are correct. The asymmetry is in *timing*: tonic is
   always on at the resting driving force, which the LIF gets right; phasic arrives
   correlated with the burst, i.e. when a real cell would be at spike voltages and the
   driving force would be up to **9.1× larger**.

**Why this is not absorbable into the hand-tuned weights.** The obvious objection is that
the weights were fitted on this substrate, so magnitudes were compensated. True, and that
is why (1) and (2) are listed first and lightly — a scalar weight rescales a mean. What a
scalar weight cannot restore is a voltage-**dependence**: it cannot turn 4.5× into 15.5×,
and it cannot make one pool's driving force swing 9× while the other's stays flat, because
both share a single weight-independent `E_rev` and a single V.

So the defect is **state-dependence**, and the exposed result is the project's central one:
the tonic/phasic ratio is a conductance claim and stands, but its *behavioural* consequence
runs through g·(E−V), and the substrate suppresses the phasic pool's driving force at the
moment it arrives. The design's deliverable is therefore not "a better neuron" — it is
**whether the selectivity ranking is substrate-independent as well as
calibration-independent**, which nobody has checked. Either answer is worth having.

Also recorded there: links 4 and 5 cannot ship separately (a 21 pF Butera cell against the
current 200 pF one makes every nS weight wrong by ~10×, silently); six acceptance tests,
including **E7 inverted** — the test that documented the LIF's limit becomes the test that
the replacement cleared it; six predicted new failure modes (E13–E18); and an explicit list
of what the upgrade does **not** fix, chiefly that it moves the >P12 muscimol anchor
*further* away by importing Butera's neonatal age.

### Review pass 7 — eight findings
All eight legitimate and reproduced before fixing. Two were latent rather than live and are
labelled as such in the code.

1. **`calibrate_pam` raised an scipy internal below its search domain.** Only `f(hi) < 0`
   was checked. The bracket is x ≥ 1.0001 (potentiation only), so `target_shift = 1.0` or a
   NAM made both endpoints positive and `brentq` raised a bare *"f(a) and f(b) must have
   different signs"* — from a function whose other failure path is a documented NaN. Now
   returns NaN, the convention callers already test.
2. **`build_kb.py --force` destroyed tables it did not own — and worse than the review
   found.** It called `DB.unlink()`, deleting the whole *file*, then recreated only its own
   8 tables. Collateral: `compounds` (53 rows, owned by `build_compounds.py`) and
   `chiral_pairs` (12 rows). I checked for a builder for `chiral_pairs`: **there is none
   anywhere in the repo**, so those 12 rows exist only in the live `.db` and would have
   been **unrecoverable**. The row-count guard I added the *last* time this script destroyed
   data did not help, because it only counts the tables the script owns — it was watching
   the wrong thing. A guard on `sources` and `findings` says nothing about `compounds`. Now
   drops only its own 8 tables; verified on a copy that both foreign tables survive.
3. **`np.percentile([], 5)` raises `IndexError`, it does not return NaN.** An arm whose
   `s_max` is unreachable yields all-NaN draws — a legitimate pharmacological outcome — and
   crashed the entire robustness report.
4. **Negative total GABA-A conductance** in `resp.py`, `rg2.py`, `circuit.py`. A NAM meeting
   sensitivity > 1 makes `eff` negative; since current is g·(E−V), a negative conductance
   inverts an inhibitory shunt into **regenerative negative damping** and the voltage
   diverges rather than failing visibly. Same class as the clamp already in
   `Drug.nmda_scale`. *Latent*: measured max sensitivity across all `PROFILES` × regions is
   0.815, so no current path reaches it — but `gaba_sens_tonic` and
   `Drug(gaba_a_gain_tonic=...)` are both public.
5. **`overdose_kinetic.py` fed a NaN affinity into `derive()`**, making `koff = NaN` and
   poisoning every rate matrix — surfacing much later as an opaque ODE failure. Now exits
   with the cause named.
6. **Locomotor FFT spacing was the module-level `DT`, not the recorded time vector.** Same
   number today, so not wrong — *unpinned*. Decimate the recording and `rfftfreq` would be
   told the undecimated spacing and report a frequency N× too fast, silently. `resp_metrics`
   already derives it correctly. A drug that slows the step cycle reading as faster is E4,
   which this project has now hit in two forms. *Latent*: needs a decimation path that does
   not exist yet.
7. **Division by zero on quiescent controls** in `calib_split.py` and `reflex.py`.
   `stretch_reflex` legitimately returns NaN gain below `IA_RESPONSE_FLOOR`, so reachable
   without a bug.
8. Dead `src` local in `conftest.py` — residue of the name-based skip heuristic review 1
   removed.

Three new regression tests cover the findings with a library surface (R7.1–R7.3).

### The gate is not met
**46 defects across seven passes: 9, 11, 8, 10, —, —, 8.** The rate is still flat. Pass 7
also found a *second* instance of the exact class of bug (destructive DB rebuild) that I had
already been bitten by and had already "fixed" — and my fix was watching the wrong tables.

That is the strongest available evidence for the sequencing argument in the design doc:
a conductance-based rewrite plus a circuit swap is a large new surface on a base whose
defect density my own review cannot measure. **Design is written; implementation stays
gated on a review pass that finds nothing.**

### Review pass 10 — three findings, and one of them was mine from pass 7
All three reproduced before fixing.

1. **A NaN latched the grid search in all three calibration sweeps**
   (`calib_split.py`, `calibrate_resp.py`, `recalibrate_kinetic.py`). The idiom was
   `if best is None or err < best[0]`. `err < nan` is **False**, so one non-finite error in
   the *first* grid cell latched `best` permanently and every later finite, better
   candidate was silently discarded. The script then printed that poisoned cell as `BEST`
   with `NaN%` beside it — which reads as a converged answer. These scripts are what the
   project's calibration constants came from, so this is a defect in the provenance of
   numbers already in use, not just in a tool.

   **And `calib_split.py` is where I made the path reachable.** My pass-7 fix turned a
   `ZeroDivisionError` on a quiescent control into a NaN `dv`, which lands directly in
   `e = abs(dv - TARGET_VENT)`. Before that fix the script crashed loudly; after it, it
   would have reported a wrong best silently. **A fix that converts a loud failure into a
   quiet wrong answer is worse than no fix, and this is the second time I have done it**
   (the first was the `evaluate` lazy-import attempt that bound a non-callable module).
   Both sweeps now test `np.isfinite` before comparing, and refuse with a stated reason
   when no cell is finite rather than reporting a winner.

2. **`locomotion()` crashed on any duration at or below the settle time.** `m = t >
   LOCO_SETTLE_S` is all-False for `duration_s <= 1.0`, so every analysis slice was shape
   `(0,)` and the first reduction raised *"zero-size array to reduction operation maximum
   which has no identity"* — an obscure numpy message from deep inside the function naming
   nothing the caller controls. `duration_s=0.5` is exactly what someone writes for a fast
   unit test or a latency sweep, so this is a likely call. Now refused at the function head
   with the settle time and the default named.

3. **Remaining unguarded baseline divisions** in `predict_muscimol.py` and
   `recalibrate_kinetic.py` (`100*mn/c`, `100*amp/ca`, `100*(v-cv)/cv`,
   `100*(am-ca)/ca`). Same class as pass 7 finding 7 and the same fix;
   `evaluation._pct_of_control` has done this correctly for several sessions while the
   scripts had not. In `recalibrate_kinetic.py` the sibling `df` was *already* guarded with
   `max(1e-9, cf)` on the line below — so the guard was present, applied to one of three
   divisions, and nobody noticed the other two.

Ten new regression tests (R10.1–R10.3). The NaN-latching one pins the *pattern* in
isolation plus a source check on each sweep, because running the three sweeps takes minutes
each.

### Running total: 49 defects across eight passes
**9, 11, 8, 10, 8, 3.** The rate is finally falling, and pass 10 was reported as the last
of this review round. But two of its three findings were in the same classes as pass 7's
(unguarded baseline division; a silent-NaN path), and one was *created by my own pass-7
fix* — so the fall is partly the reviewer running out of new surface, not the code becoming
clean. The gate in `knowledge/05-design-conductance-substrate.md` asks for a pass that
finds **nothing**; 3 is not 0.

---

## Session 13 (cont.) — LINK 4 BUILT: conductance-based cell

`src/circuitpharm/neuron.py`, `tests/test_conductance_cell.py` (21 tests). The design in
`knowledge/05-design-conductance-substrate.md` was written first and the acceptance tests
(A1–A7) and predicted failure modes (E13–E18) were specified before any code existed.

### Reading the parameters off the source caught a real error
The design doc flagged my recalled Butera–Rinzel–Smith values as UNVERIFIED and said "read
them off the paper, do not recall them." Doing that found one wrong number:

**E_L is −57.5 mV, not −65 mV.**

That is not cosmetic. E_L is this model's bifurcation parameter — the knob that moves the
cell quiescent → bursting → tonic. Starting from −65 mV would have put the cell in the
wrong regime, and the rhythm would then have been retuned around a wrong resting drive,
which is the exact shape of this project's recurring failures.

How narrow that window is, measured:

| E_L (mV) | −65.0 | −60.0 | **−57.5** | −55.0 | −52.5 |
|---|---|---|---|---|---|
| regime | **silent, 0 spikes** | tonic | **bursting** | tonic | tonic |

−65 mV is not merely the wrong regime — it is a dead cell. And the bursting band is only a
few mV wide, with tonic spiking on both sides, which is what makes this the single parameter
most expensive to get wrong.
 Everything else in the
recalled set was right (C 21 pF; g_Na 28, g_K 11.2, g_NaP 2.8, g_L 2.8 nS; E_Na +50,
E_K −85; gate θ/σ; τ̄_n 10 ms, τ̄_h 10 s), as were the current equations and exponents
(`g_Na m∞³(1−n)`, `g_K n⁴`, `g_NaP mp∞ h`). Source: the curated CellML encoding, which is
machine-readable; the journal full text returned HTTP 403. `test_E14_the_one_parameter_that_
was_misremembered_is_right` pins it, and also asserts that the wrong value does **not**
burst — so if the pin ever stops testing what it claims, that shows up too.

### The acceptance tests, and one I had to replace honestly
**A1, E7 inverted.** The E7 test asserts the LIF's V never exceeds threshold and says "if
this ever fails, voltage-gated mechanisms become available and the modelling choice should
be revisited." This is that revisit from the other side. Measured: **V peaks at +6.49 mV**,
h-gate span **0.4544** against the LIF's 0.01.

But the design doc pre-committed a span of **>0.5**, and the cell gives 0.454. The 0.7
figure in the original E7 note was **never sourced** — it was my own estimate of what
inactivation would need. So the pre-committed number would have failed a cell that is
behaving correctly.

Rather than quietly relax the bar I replaced it with a functional one, which is what the
0.7 was a proxy for anyway: **freezing h abolishes bursting in every direction.**

| h | regime | rate |
|---|---|---|
| free | **bursting**, 12 bursts | 4.9 Hz |
| frozen 0.46 | quiescent | 0 Hz |
| frozen 0.60 | tonic | 68.8 Hz |
| frozen 0.90 | tonic | 118.9 Hz |

So the rhythm comes from the mechanism claimed — slow voltage-dependent inactivation of a
persistent sodium current — not from something incidental. The LIF can neither pass nor fail
this test; it has no h gate to clamp. That is the acceptance criterion for link 4. The span
assertion is kept only as a coarse floor at 30× the LIF's value, with the measured number
recorded so drift is visible.

**A2.** Bursts with `g_adapt = 0`, which is now the default (E18). The LIF cannot burst
without spike-triggered adaptation. And `g_NaP → 0` makes the isolated cell **quiescent**,
so the current is load-bearing. The other half of that prediction — that the NETWORK rhythm
should *persist* without I_NaP (the pacemaker-vs-network controversy) — needs the coupled
population and is deliberately **not** asserted yet; asserting half a prediction and calling
it validated is how the retracted Matsuoka result happened.

**A3, the published excitability sequence**, reproduced and not fitted by us:

| i_app (pA) | −30 | −15 | −5 | 0 | +10 | +25 | +50 |
|---|---|---|---|---|---|---|---|
| regime | quiescent | quiescent | **bursting** | **bursting** | tonic | tonic | tonic |

**E13.** Burst period moves **0.02%** between dt_max = 0.05 and 0.0125 ms. Converged. Note
dt_max 0.2 and 0.1 give *identical* output at an external dt of 0.1 — both reduce to one
substep — so that agreement is not evidence of dt-independence, and a second test pins the
substep count to dt_max to stop that being misread.

**E15.** Mg²⁺ relief now spans **11.80×** (0.071 → 0.842) against the **4.50×** measured on
the LIF preBotC, with peak/rest 6.52×. The design doc projected 15.5× for a cell spanning
−65..+20 mV; this cell spans −61.9..+6.5, so the realised gain is 2.6× rather than 3.4×.
Reporting the measured value, not the projection.

### Two traps found while building, both silent
1. **The noise term.** I wrote the current-noise amplitude as `sigma/sqrt(hs)` and then
   multiplied back by `sqrt(hs)`, so the two factors cancelled and the noise did **not scale
   with the substep at all**. Halving dt would then have quietly changed the noise amplitude
   and contaminated the E13 dt-convergence test with the very thing that test exists to
   detect. Now written as the Wiener increment directly.
2. **The i_app sign differs from the source encoding.** The CellML writes
   `dV/dt = −(i_NaP + i_Na + i_K + i_L + i_tonic + i_app)/C`, putting i_app *inside* the
   negated sum, so there a positive i_app **hyperpolarises**. Here it is outside, matching
   `cpg.Pop` where a positive `Idrive` depolarises. Matching the LIF is right — the whole
   point is that a circuit can swap substrates without changing its stepping loop — but
   anyone comparing against published figures must flip the sign, and a sign error on the
   bifurcation parameter would move the cell between regimes while looking plausible.

### Compute: measured 4.1×, not the estimated 20×
`CondPop` vs `cpg.Pop`, 2 s simulated: **3.67× at n=1, 4.12× at n=50.** The estimate reasoned
"~10 transcendentals against ~1, times 2 substeps"; both factors were real, and the error was
treating them as the whole cost. Shared per-step overhead (the receptor loop, RNG draws,
numpy dispatch on small arrays) dominates at these population sizes. **This loosens the
budget enough that running the conductance arm at the full 20,000 draws is plausible**, so
§6's reduced-draw compromise should be re-examined against a timed network run rather than
assumed.

### E14: mixed cell/weight pairs are now unconstructible
The predicted most-likely defect of this migration. The LIF cell is C=200 pF, g_L=10 nS;
this one is C=21 pF, g_L=2.8 nS. Every synaptic weight in the package is in nS, hand-tuned
against the LIF. `ParamSet` refuses a weight table whose provenance does not name its cell,
so the silent order-of-magnitude error is structurally impossible rather than guarded by a
comment.

### What is NOT done
Link 5. The cell is built and characterised **in isolation**; no circuit uses it yet.
`resp.py`, `rg2.py` and `circuit.py` still construct `cpg.Pop`, and they assign `g_adapt`
*after* construction — which would switch adaptation back on and re-introduce E18 the moment
a circuit is pointed at `CondPop`. That is link 5's first problem, not a defect in link 4.
The substrate-independence comparison (§6), which is the actual deliverable, needs the
network.

---

## Session 13 (cont.) — LINK 5: conductance cell in the network

### I had the deliverable of links 4–5 wrong
Recorded first because it is the most important thing in this session.

The design doc, `KNOWLEDGE.md` and my report to the user all said the point of links 4–5 was
testing whether **the selectivity ranking** — the project's one VALIDATED result — is
substrate-independent as well as calibration-independent.

**That is false.** `scripts/ranking_robustness.py` never imports a circuit module. Its
`score()` is pure algebra over subunit expression fractions, measured efficacies,
extrasynaptic fractions and a gain ratio. Its own docstring says so on line 21: *"This is
analytic — no circuit simulation — so it runs at thousands of draws."* I had read that file
earlier in the same session and quoted a different line from it.

So the ranking is substrate-independent **by construction** — it cannot change when the
neuron model changes, because it never runs a neuron. **The VALIDATED headline was never at
risk from the substrate, and links 4–5 cannot strengthen it.**

**Corrected deliverable.** What does run through the circuit is everything
`evaluation.evaluate()` simulates — the UNCALIBRATED respiratory and motor endpoints, where
shape and *ordering between compounds* are usable and absolute scale is not. Those go
through g·(E−V), which is what the LIF truncates. So the question is **whether the ordering
of compounds by simulated respiratory burden survives the substrate change**
(`scripts/compare_substrates.py`). Genuinely at risk; a smaller prize. Link 5's reach is the
UNCALIBRATED tier only. Corrected in the design doc §6 and KNOWLEDGE.md, both marked as
corrections rather than silently rewritten.

### The published coupling weights do not exist in reachable form
Butera–Rinzel–Smith **part II** supplies network coupling conductances for a population of
exactly these cells, which would have made link 5's weights *published* rather than ours —
the entire point of link 5. They are not retrievable: the journal full text returns HTTP
403, and every accessible encoding (the curated CellML, ModelDB 247647 via its GitHub
mirror) is **single-cell only** — `g_tonic_e = 0.0`, no coupling term.

Inventing them was not an option. So link 5 delivers a **published CELL in a network whose
COUPLING is ours**, anchored to a matched operating point. Recorded as a partial completion,
not a tick, and said plainly in `config.COND_RESP_OP` and `scripts/anchor_cond_resp.py`
because a reader would otherwise reasonably assume the whole network came from the paper.

### Plumbing: `PreBotC(substrate="cond")`
Building the switch forced out **three LIF-scaled quantities** that would each have been
silent:

1. `drive` — pA against a C=200 pF / g_L=10 nS cell.
2. `gaba_tonic` — nS, likewise.
3. **A hardcoded `60.0`** supplying the entire excitatory drive to the Inh and Out
   populations. Found *after* the first two had been moved into the operating point, which
   is the point: two of three populations would have run at roughly the wrong order of
   magnitude while the rhythm still looked fine.

So the conductance substrate **refuses to start** without a complete operating point, and
rejects a *partial* one too. The partial case is the dangerous one: supply `drive` and the
excitatory weights while `ie_gaba`, `ie_gly` and `gaba_tonic` keep their LIF values, and the
network is inconsistent by ~an order of magnitude on the **inhibitory arm — the arm the drug
acts through**. It would run, produce a rhythm, and be wrong about pharmacology.

**E18 handled:** the LIF assigned `g_adapt` after construction, so pointing it at a
conductance cell would have switched spike-triggered adaptation back on at its LIF-tuned
2.5 nS/spike on top of a now-real I_NaP inactivation — double-counted burst termination,
symptom a plausible duty cycle.

**E15 handled:** on the cond substrate NMDA is handed over **unevaluated** (`raw=("nmda",)`)
so the cell applies the Mg²⁺ block at its own V each 0.05 ms substep. Pre-evaluating it at
the step-entry voltage would have discarded the entire voltage-dependent relief, and the
symptom would have been "the substrate change did not move the NMDA result" — a reassuring
robustness check that is in fact the bug.

**A7 verified, not assumed:** the LIF `PreBotC` trace hashes identically to the committed
pre-refactor module (`cb808accc11db8462b7f` both sides).

### A bug I was about to introduce
`evaluation._CTRL_CACHE` is keyed on **seed alone**, argued safe because a drug-free run
gives `eff = 1` for any sensitivity. That argument covers the *drug* parameters only;
substrate is a different axis and changes the drug-free control completely. Keyed on seed
alone, a conductance arm would have been normalised against a **LIF control**, making every
percent-of-control silently nonsense — and the substrate comparison is precisely a
comparison of those percentages, so **the one number the upgrade exists to produce would
have been the one corrupted.** Substrate is now in the key, and a test pins it.

### Anchoring: the viable window is narrow
`scripts/anchor_cond_resp.py`, three stages. Target is the LIF's *measured* control
(1.271 Hz, mod 4.74) rather than a literature value — the question is whether the same drug
gives the same *fractional* change, which needs the two controls to agree on the observable
and says nothing about either being the right absolute frequency.

* **Stage 1** (50 points, 8 s): **3 alive.**
* **Stage 2** (81 points, 12 s): **16 alive.** Best 1.356 Hz.
* **Stage 3** (finalists, 30 s, 4 seeds, every seed required alive): in progress.

**3/50 is a result, not a search artifact.** A population of intrinsically bursting Butera
cells has a far smaller region of synchronised in-band behaviour than the LIF network did,
because each cell is already an oscillator and the coupling has to *entrain* rather than
*create* the rhythm. That is the correct biology for a coupled-pacemaker preBötC, and it
means this operating point is more fragile than the LIF's.

**Short sweeps are not trustworthy for frequency:** the stage-1 winner read 1.026 Hz at 8 s
and 1.102 Hz at 12 s — a 7% move from duration alone. Hence stage 3, and hence the rule that
a sweep winner is a *candidate*, never the answer.

**Also noted:** the cond control's mean inspiratory output is ~5 against the LIF's ~29. The
absolute level is UNCALIBRATED on both, so the comparison is unaffected, but the compressed
range leaves less room between a drug effect and the `mean < 1.0` collapse floor. Worth
watching when reading the comparison.

---

## Session 13 (cont.) — the anchoring was wrong, and the provenance hole is deeper than recorded

### I produced this project's signature failure, in the session that wrote the warning
`COND_RESP_OP` was anchored, verified across 4 seeds at 30 s with every seed required alive,
read 1.342 ± 0.132 Hz against the LIF's 1.271, and registered. **It was wrong**, for two
compounding reasons, and it is now `None` again with the reasoning written into `config.py`.

**1. The search never visited the bursting regime.** The grid swept Exc drive over
{5, 10, 15, 20, 25, 30} pA. Measured with `neuron.run_isolated`:

| i_app (pA) | −15 | −10 | **−5** | **0** | +5 | +10 | +20 | +30 |
|---|---|---|---|---|---|---|---|---|
| regime | quiescent | quiescent | **bursting** | **bursting** | tonic | tonic | tonic | tonic |

Every point searched held the cells depolarised out of pacemaking. And "only 3 of 50 points
were alive", which I recorded as *the correct biology of a narrow entrainment window*, was
really the signature of searching almost entirely outside the regime the cell can oscillate
in. I wrote a confident mechanistic explanation for an artifact of my own grid.

**2. Every measurement was inside a transient.** τ_h is **10 seconds** — the slowest
timescale in the Butera cell, 25× the LIF's τ_adapt of 400 ms. Every conductance measurement
used a 4 s warm-up, i.e. **0.4 time constants**. Measured in successive 10 s windows:

| t₀ (s) | 0 | 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 |
|---|---|---|---|---|---|---|---|---|---|
| **cond** freq | 1.33 | 1.33 | 2.86 | 2.86 | 3.67 | 0.82 | 0.31 | 0.82 | 0.61 |
| **cond** mean | 11.7 | 2.4 | 2.7 | 3.2 | 2.2 | 2.4 | 3.5 | 2.6 | 3.1 |
| **cond** mean h | 0.265 | 0.080 | 0.080 | 0.086 | 0.082 | 0.081 | 0.087 | 0.078 | 0.090 |
| **LIF** freq | 1.327 | 1.326 | 1.327 | 1.327 | 1.224 | 1.327 | 1.327 | 1.326 | 1.326 |
| **LIF** mean | 28.7 | 29.1 | 29.3 | 29.4 | 28.2 | 27.0 | 27.3 | 28.2 | 29.3 |

The conductance network's mean h collapses 0.265 → 0.080 over the first 20 s and the rhythm
then wanders, alive in 5 of 9 windows. The LIF holds 1.327 Hz in every window.

So the registered number was **a plausible reading of a decaying transient**, produced by a
three-stage search with a seed-robustness check — none of which could see it, because all
three stages measured inside the same transient. A longer verification would not have helped
either: stage 3 used 30 s, which is still only 3 τ_h *total*, with 4 s of it as warm-up.

**Consequences:** the A4 I_NaP-dissociation negative result is **withdrawn** — it was measured
on the invalid operating point with a 4 s warm-up, so it tests nothing. The substrate
comparison is blocked until re-anchored.

**The anchoring script now:** searches drive where the isolated cell bursts (verified with
`run_isolated`, not assumed); warms up 3 τ_h = 30 s; measures two consecutive 15 s windows and
**gates hard on frequency agreement between them**, so a drifting network cannot pass.

### A real inefficiency in `compare_substrates.py`
`burden()` recomputed the drug-free controls **once per arm** — six times the necessary work
on the default arm list, on the substrate that costs 4× per simulation. Controls do not depend
on the arm. This is what made the first runs look hung rather than merely slow. Fixed, and the
script now prints per-arm as it goes and flushes, because a twelve-minute run that prints only
at the end is indistinguishable from a hung one.

### The provenance hole: 6 of 28
The recorded gap was "90 sources marked UNVERIFIED". That understates it by a layer.

**The numbers driving the one VALIDATED result are not unverified — they are unsourced.**
`scripts/ranking_robustness.py` scores on exactly `PROFILES[...].eff()`, `REGIONS`,
`SUBJECTIVE_WEIGHT` and `EXTRASYN`. `subtypes.py` *appears* to cite five sources; all five —
`alogabat`, `gl_ii_73`, `imepitoin`, `sh053`, `tpa023` — are **compound names that collide
with source keys**. The true citation count in the module holding the project's central
numbers is **zero**. The database could not have supplied them either: `subunit_expression`
holds 8 rows and every one is qualitative ("predominant", "present", "enriched").

`src/circuitpharm/provenance.py` now records a basis for each: **6 of 28 name a source (21%);
18 UNSOURCED, 8 FROM_QUALITATIVE, 1 FITTED, 1 GUESS.** `tests/test_provenance.py` (12 tests)
pins the count so it cannot move quietly in either direction, and fails if a new parameter is
added without a record or if a record names a source absent from the database.

**The worst single entry:** `REGIONS["prebotc"]["a5"] = 0.02`. An α5-selective compound's
modelled respiratory burden is roughly proportional to it, so it sets the entire safety
margin — and no source in the knowledge base states it.

**`a5_dist` is the trap.** It is the obvious source for the regional fractions, it resolves
perfectly by DOI (title similarity 1.00), and it **does not support them**: a 1988 study using
a single probe for *"the alpha subunit"*, reporting TOTAL α mRNA by region (cerebellum >
thalamus = cortex = hippocampus ≫ pons = striatum = medulla). It neither resolves subtypes nor
measures composition — `REGIONS` rows sum to 1 by construction, so they encode composition,
while the paper measures level. A metadata-only pipeline would have marked it VERIFIED. That
is why `scripts/verify_sources.py` reports `METADATA_*` and `CLAIM_SUPPORT_UNASSESSED` as
separate fields and never collapses them.

### Is the headline an artifact of its prior? No — but half of it is
The Dirichlet prior puts **38% of draws below a tenth of the nominal preBötC α5**, and the
forebrain:preBötC α5 ratio has a p99 of ~3.5×10⁷, so the α5 arms' scores should diverge.
Flooring the drawn fraction at 0.005 / 0.01 / 0.02 (= nominal):

| floor | alogabat p5 | median | p95 |
|---|---|---|---|
| none | **2.51** | 16.05 | 84.9 |
| 0.005 | **2.51** | 13.56 | 36.0 |
| 0.01 | **2.51** | 10.83 | 23.5 |
| 0.02 | **2.51** | 7.63 | 14.5 |

`corr(score, 1/α5) = 0.02`. **The 5th-percentile floor and the ordering are robust; my
suspicion was wrong.** But the median and p95 are prior artifacts — the previously-quoted
*"median 33.76×, p95 9264.80×"* for the ideal α5 arm are properties of the prior's tail, not
pharmacological claims, and are now marked NOT QUOTABLE in the script's own output.

Also found: for a perfectly α5-selective compound the burden can reach zero, so the score is
undefined, and the `isfinite` filter silently discarded **exactly the draws most favourable to
that arm**. Conservative, which is the honest direction — but it was invisible. Now counted
and reported.

**What remains unrescued:** `KAPPA = 15` and `EXTRASYN_CONC = 8` are unsourced choices of
mine. "Robust across 20,000 draws" means robust across a prior whose centre *and* width I
invented.

### Source verification, mechanised
`scripts/verify_sources.py` resolves each source by DOI (Crossref) or PMID/PMCID/title
(Europe PMC) and compares the resolved title against the recorded citation. Two bugs in my own
script, both found because sources I had *just cited* came back UNRESOLVED:

1. **URL-encoded DOIs were invisible.** PLOS writes `?id=10.1371%2Fjournal.pone.0030608`; the
   DOI regex needs a literal slash. `pbc_eps` went UNRESOLVED → RESOLVED_MATCH after
   URL-decoding.
2. **Nature URLs carry the DOI suffix without the prefix.** `/articles/s41598-017-17379-x` is
   `10.1038/s41598-017-17379-x`. `pbc_delta` and `vrodent` recovered.

A verification tool that under-reports resolution is not a safe failure: it makes the
provenance gap look worse than it is, which misleads in the opposite direction.

Also added `NON_LITERATURE` for trial registries, encyclopedia pages and software resources —
reporting a clinicaltrials.gov NCT number as UNRESOLVED conflates "could not check" with "not
the kind of thing checked this way".

### Why the conductance network collapsed: two causes at once, pulling opposite ways
The corrected anchoring (bursting-regime drive, 30 s warm-up, drift gate) found **0 of 48
points alive**, every one "mean output collapsed". Diagnosed rather than widened, by probing
one point and stripping components:

| configuration | Exc | Inh | Out |
|---|---|---|---|
| as searched (gaba_tonic 1.5, d_oth 25) | 0 Hz (−61.9 mV) | 0 (−57.8) | 0 (−57.8) |
| **gaba_tonic = 0** | 34.5 Hz (−28.4) | 57.1 (−29.3) | **0 (−17.8)** |
| gaba_tonic = 0, no synaptic inhibition | 11.3 (−22.9) | 111.3 (−32.9) | **0 (−22.1)** |

**Cause 1 — a conductance cannot be carried between cells by a weight scale.** The silenced
voltage is exactly predictable: with `g_L = 2.8 nS` to −57.5 mV, `gaba_tonic = 1.5 nS` to
−75 mV and +25 pA, steady state is `2.8(V+57.5) + 1.5(V+75) = 25` → **V = −57.8 mV**, which
is what was measured. The tonic GABA was cancelling the drive exactly.

A conductance means "this fraction of the cell's leak". The LIF's `gaba_tonic = 1.5 nS` is
**15% of its 10 nS leak**; on the Butera cell it is **54% of 2.8 nS**. The correct scale
factor is the g_L ratio, **0.28** — so 1.5 nS becomes 0.42 nS. My grids searched 0.75–2.25 nS,
i.e. **1.8× to 5.4× too strong**. The "weight-block scale" of 0.5–2.0 I had been sweeping
cannot fix this, because what matters is the ratio to g_L and that differs by 3.6× between
the cells.

**Cause 2 — this cell has an UPPER bound on drive, and the LIF has none.** With
`gaba_tonic = 0`, `Out` sits at **−17.8 mV and fires at 0 Hz**: depolarisation block. At that
potential the Na inactivation term `(1−n)` is nearly zero, so no spike can be produced however
much current is injected. **The LIF cannot do this at all** — more drive simply means a higher
rate, capped only by `tref`.

So `mean output collapsed` had two distinct causes operating simultaneously and in **opposite
directions**: tonic GABA silencing Exc from below, and `drive_other` blocking Out from above.
A one-dimensional reading of either leads to the wrong correction — and my first *derived*
grid did exactly that, raising `drive_other` to 40–70 pA on a threshold calculation, which
would have driven Out further into block. Corrected to 5–20 pA, with Inh and Out recruited by
`eo_ampa`/`ei_ampa` from Exc rather than by injected current.

This is the clearest single lesson of the migration: **the LIF→conductance change is not a
rescaling.** Three distinct quantities (conductances, the Exc drive, the Inh/Out bias) each
need a different transformation, and two of them are bounded above as well as below.

### The cause: NMDA, and it is the §1 finding biting back
Two derived grids failed (0/48, 0/72), both "mean output collapsed". Diagnosed by
instrumenting per-population state across the range instead of sweeping further:

| scale | Exc | Inh | Out |
|---|---|---|---|
| 0.056 | 112–119 Hz, V −30 | 0.6–56 Hz | 124–127 Hz, V −27 |
| 0.280 | 85–87 Hz, V −27 | 125 Hz | **0 Hz, V −19.9** |

Two regimes, neither usable: at 0.28 `Out` is in **depolarisation block**; at 0.056 everything
fires **tonically at 110–127 Hz**. `Exc` never bursts at any scale — V sits at −27 to −31 mV,
far above its bursting range.

**The cause is NMDA, and the arithmetic is exact.** `tau_nmda` is **100 ms, 20× the AMPA
tau**, so with lumped population weights (`inject` adds `w × total_spike_count` to every
postsynaptic cell) the NMDA conductance accumulates 20× more per unit firing rate. At scale
0.28 and 112 Hz/cell:

```
g_ampa  =  0.126 nS × 5600 spikes/s × 0.005 s =  3.53 nS
g_nmda  =  0.069 nS × 5600 spikes/s × 0.100 s = 38.8 nS  (unblocked)
```

against `g_L = 2.8 nS`. Whether that silences the cell depends entirely on the Mg²⁺ relief:

| Mg relief | total g | V_ss |
|---|---|---|
| **0.063** — the LIF's measured mean | 5.97 nS | −18.4 mV |
| **0.70** — a depolarised conductance cell | 30.7 nS | **−4.8 mV** |

**So the weights were tuned on a substrate where NMDA was effectively 6% of nominal.** §1 of
the design document identified precisely this as the LIF's central defect: *"NMDA conductance
is held in near-permanent Mg²⁺ block — mean relief 0.063, dynamic range 4.5× against the
15.5× a spiking cell has."* Fixing that defect — correctly, by evaluating the block at the
cell's own voltage each substep (E15) — made the inherited NMDA weights **11× too strong.**

And it is positive feedback, not merely a scale error: relief *rises* with depolarisation, so
more NMDA current depolarises further, relieving more block. The network runs away into
depolarisation block, which the LIF cannot do at all.

**Derived correction — the receptor types need different scales:**

```
ee_ampa_cond = ee_ampa_LIF × (g_L ratio)                     = LIF × 0.280
ee_nmda_cond = ee_nmda_LIF × (0.063/0.70) × (g_L ratio)      = LIF × 0.025
```

**an 11× difference between AMPA and NMDA.** This is the third and sharpest instance of the
same lesson: a single multiplicative weight scale cannot carry a weight table between these
two cells. The quantities needing independent transformation are now four — conductances by
the g_L ratio, NMDA additionally by the relief ratio, the Exc drive from the bursting window,
and the Inh/Out bias bounded above by depolarisation block.

Worth stating plainly: **the LIF's respiratory NMDA results were produced with NMDA at ~6% of
its nominal conductance and unable to relieve.** That does not overturn them — the weights
were fitted in that regime, so the *net* excitation was right — but it means the NMDA arm's
state-dependence was absent, which is exactly what the already-VOID NMDA respiratory
contribution was declared VOID for. The substrate change does not rescue that result; it
explains why it was unrescuable.

### The conductance preBötC is anchorable. Six causes, five silent, each needing a different fix
Six anchoring attempts. Each step came from a diagnosis, not a wider grid — and the one time
I widened derivedly instead of diagnosing, I made it worse.

| attempt | alive | cause found |
|---|---|---|
| 1 | 3/50 | searched outside the bursting regime **and** measured inside a 10 s transient |
| 2 | 0/48 | tonic GABA 1.8–5.4× too strong — conductances scale by the **g_L ratio**, not a free weight scale |
| 3 | 0/72 | `Out` in **depolarisation block** — this cell has an *upper* bound on drive; the LIF has none |
| 4 | 9/81 | **NMDA 11× too strong** — τ=100 ms, and the LIF's weights were tuned against 6% Mg relief |
| 5 | 16/36 | anchored to an **in vivo** frequency; the Butera cell is neonatal **in vitro** |
| 6 | **41/48** | optimising modulation alone selected a rhythm no drug could be measured against |

**What the attempts 2–4 share:** a single multiplicative weight scale cannot carry a weight
table between these two cells. Four quantities need four different transformations:

```
AMPA / GABA / glycine  ->  x (g_L ratio)              = 0.28
NMDA                   ->  x (g_L ratio) x (0.063/0.70) = 0.025      <- 11x smaller
Exc drive (pA)         ->  set from the measured bursting window
Inh/Out bias (pA)      ->  small; bounded ABOVE by depolarisation block
```

The NMDA factor is the interesting one: it is the ratio of the LIF's *measured* mean Mg²⁺
relief (0.063) to a depolarised conductance cell's (~0.70). **The design document's §1 named
that 0.063 as the LIF's central defect; fixing it correctly made the inherited weights 11×
too strong.** And because relief rises with depolarisation it is positive feedback, not a
scale error — the network runs away into block, which the LIF cannot do.

**Attempt 5 was a wrong target, not a wrong network.** `EUPNOEA_BAND` (0.30–2.50 Hz) is an
*in vivo* rat band and correct for the LIF, tuned to 1.27 Hz. Neonatal rat preBötC slices
run at 6.6 ± 3.1 to 14.6 ± 2.0 bursts/min ≈ 0.11–0.24 Hz (`pbc_invitro_freq`: Revill et al.
2021, Front Physiol 12:626470). Forcing an in vitro preparation 5–12× above its own physiology
is what destroyed modulation (0.9 vs the LIF's 4.7) and stationarity (drift 40–164%). The
source was read and entered in the knowledge base *before* the band was introduced, the new
`INVITRO_BAND` is a separate constant, `EUPNOEA_BAND` is untouched, and the new band still
excludes the fragmented 3.8–4.9 Hz rhythms the gate exists to catch. A test
(`test_the_in_vitro_band_is_sourced_not_convenient`) fails unless that source carries a
hand-read verdict — because a band added to rescue a failing point looks identical in code to
one added because the preparation differs, and only the record distinguishes them.

**The consequence that travels with it:** the two substrates no longer share a frequency, so
absolute frequency comparisons between them are meaningless. Only fractional change from each
substrate's own control is comparable — which design doc §6 required all along. Matching the
frequencies was my addition, and invalid.

**Attempt 6 and a threshold I had to correct against myself.** `resp_metrics` calls an arm
dead below `max(1.0, 0.2 × control)`, so a graded measurement range exists only when
`ctrl_mean > 5.0`. That is derived. I set the gate at **10.0** "for margin", and at
verification it disqualified the best rhythm in the finalist set — modulation 4.44 (closest to
the LIF's 4.736) and the tightest frequency SD of the four (0.015 Hz) — for coming in at mean
9.5, i.e. **0.5 below a number I had chosen**, while satisfying the derived criterion with 7.6
of graded range. Selection then fell to a rhythm with modulation 2.43.

The gate is now the derived 5.0 with 10.0 as a scored preference. Recorded in the code and
here because adjusting a threshold after seeing which candidate it excludes is precisely the
move that needs its reason on the record, and the reason must be that 5.0 is derived and 10.0
was not — not that I preferred the outcome. I could see which point the change favoured before
making it. `MIN_CTRL_MEAN` is now marked as the one threshold that must not move.

### Also fixed: the three substrate-coupled settings
`_simulate_resp(substrate="cond")` was still using the in vivo band and the LIF's 4 s warm-up,
so every conductance evaluation would have mis-gated healthy rhythms as dead *and* measured
inside the transient. Operating point, validity band and settling time now live in one
`_SUBSTRATE` table with four tests pinning the coupling, including that the cond warm-up is at
least 3 τ_h.

### A4: the dissociation is NOT reproduced, and that is an inherited limitation, not a bug
Re-run on the anchored operating point, 60 s with a 30 s settle, gated against the in vitro
band:

| | alive | freq | mean | mod |
|---|---|---|---|---|
| isolated cell, g_NaP 2.8 nS | bursting, 37 bursts | — | — | h-span 0.4544 |
| isolated cell, g_NaP 0 | **quiescent** | — | — | h-span 0.0820 |
| network, g_NaP 2.8 nS | **3/3** | 0.280 Hz | 9.73 | 4.24 |
| network, g_NaP 0 | **0/3** | 0.000 | collapsed | 0.00 |

In vitro, riluzole abolishes isolated-cell pacemaking while the **network rhythm persists**.
Here it abolishes both. The dissociation is not reproduced.

**This was predicted before the run**, in `anchor_cond_resp.py`'s stage-2 comment: *"if the
viable scale turns out to be this low, the rhythm in this network IS pacemaker-driven, and A4
will fail STRUCTURALLY rather than numerically."* The viable coupling is NMDA × 0.025 and
AMPA × 0.40 — weak recurrent excitation — so the rhythm rests entirely on the cells' intrinsic
pacemaking. Remove I_NaP and the cells cannot fire at all, so there is nothing for the network
to synchronise.

**And it is an inherited limitation rather than an implementation defect.** Butera–Rinzel–Smith
model 1 *is* the pacemaker hypothesis: bursting arises from fast activation and slow
inactivation of I_NaP in individual cells. The riluzole experiments are the principal published
argument *against* pacemaker-driven rhythmogenesis. A faithful implementation of a pacemaker
model therefore *must* fail this test — passing it would mean we had implemented something
other than the model we adopted.

So A4 is evidence that the implementation is faithful, and simultaneously a real limitation of
what link 5 bought: the network reproduces a published *cell* correctly, including that cell
model's known inability to account for the riluzole result. Reaching rhythm persistence under
I_NaP block would need either much stronger recurrent excitation (which depolarisation block
forbids on this cell) or a burst-terminating mechanism independent of I_NaP — Ca-dependent K
current, or synaptic depression. Both are additions beyond link 5.

**Not retuned.** Fitting the coupling until the rhythm survives and then citing the survival as
validation is how the Matsuoka locomotor prediction had to be retracted.

### THE DELIVERABLE: the subtype ordering is substrate-independent; the non-selective BZ's SIGN is not
`scripts/compare_substrates.py`, 3 seeds, occupancy 1.0, each arm as a fractional reduction
in mean inspiratory output from **its own** substrate's control.

| arm | LIF | COND | Δ |
|---|---|---|---|
| neurosteroid | +0.610 | +0.994 *(alive 0/3)* | +0.384 |
| mp_iii_022 | +0.020 | +0.071 | +0.051 |
| ideal_a5 | +0.005 | +0.051 | +0.046 |
| alogabat | +0.002 | +0.043 | +0.041 |
| hz_166 | +0.001 | −0.100 | −0.101 |
| **nonselective_bz** | **+0.079** | **−0.377** | **−0.456** |

The script's verdict ("ORDERING MOVED", Spearman +0.43) buries the actual result:

* **Spearman over the five SUBTYPE-SELECTIVE arms = +1.0000.** Identical order on both
  substrates: neurosteroid > mp_iii_022 > ideal_a5 > alogabat > hz_166. Two neuron models
  that differ in every mechanism — real I_NaP, voltage-gated inactivation, 11.8× Mg²⁺ relief
  span, depolarisation block, versus none of those — give the same ranking.
* **Exactly one arm moves, and only because its sign flips.** The non-selective BZ goes from
  2nd-most-burdensome (+0.079) to least (−0.377, a 38% *increase* in output). It is the arm
  with the largest α1 efficacy, acting on a preBötC the model treats as α1-predominant
  (`REGIONS` α1 = 0.60) — so the flip is largest exactly where the drug effect is largest.

So the ordering result is a genuine (if narrow) strengthening of the UNCALIBRATED tier for
subtype-selective compounds, and a clear warning for the non-selective reference arm.

### A mechanism I claimed, then had to retract on resolution grounds
I attributed the flip to GABA-A shunting shortening bursts → less I_NaP inactivation → faster
recovery → higher rate, and presented frequency readings as confirming it:

| | difference | FFT bins |
|---|---|---|
| COND 0.302 → 0.336 Hz | 0.034 Hz | **1.02 bins** |
| LIF 1.327 → 1.327 Hz | 0.000 Hz | below 1 bin (0.102 Hz = 7.7%) |

**Both at the resolution limit.** The conductance "+11.1%" is one bin, indistinguishable from
quantisation; the LIF's "+0.0%" means only "under 7.7%", not "frequency-insensitive". I had
stated both as findings. Retracted.

### Direct burst counting: what survives, and a fourth inherited-protocol failure
`burst_metrics` over 180 s (cond) / 90 s (LIF) instead of the FFT:

| | COND | LIF |
|---|---|---|
| mean output | **+41.7%** (comparison said +35.3%) | **−8.6%** (said −7.4%) |
| duty cycle | **+65.5%** | −6.7% |
| peak | −7.3% | −1.4% |

**Robust, by two independent methods:** the sign flip is real, it is a change in *how much of
the time* the network is active, and peak burst height is essentially unchanged on both.

**Not established:** whether that is more bursts or longer bursts. `burst_metrics` reported
5.63 Hz on cond against the FFT's 0.302 Hz — a 19× discrepancy — because the 20 ms rate
low-pass leaves intra-burst fluctuation and a 0.35×max threshold fires repeatedly *within* one
respiratory burst. So its "rate" is an event rate, not a burst rate, and the +25.3% / −17.3%
cannot be read as burst frequency. A burst-level detector is needed.

**That is the fourth measurement protocol inherited from the LIF that needed re-deriving for
this substrate**, after the 4 s warm-up against a 10 s τ_h, the in vivo band on an in vitro
preparation, and the FFT resolution. The pattern is consistent enough to state as a rule: on
the conductance substrate, assume every inherited protocol is wrong until re-derived. None of
the four failed loudly; each produced a plausible number.

---

# Session 8: the E20 regression, and reviewing the manuscript against the model

## E20 — the substrate factory's default erased the LIF cell's adaptation

The `substrate.py` migration left the suite at **2 failed, 258 passed**:
`test_benzodiazepine_depresses_the_reflex` and
`test_sedative_degrades_coordination_and_slows_the_step_cycle`. Both LIF-path motor
phenotypes, both green at 251/251 before the refactor, so a regression I introduced.

**Cause, exactly as predicted before the break.** `sub.make_pop` carried `g_adapt=0.0` as a
default. `resp.py` and `rg2.py` had always set adaptation explicitly, so routing them through
the factory was faithful. `circuit.py` had **not** — it built bare `Pop`s and set only `tref`,
inheriting the dataclass's `g_adapt = 0.55` — so the factory's default silently deleted
spike-triggered adaptation from every spinal population.

**The fix is structural, not a number.** Writing `g_adapt=0.55` in `circuit.py` would have
duplicated `Pop`'s default into a second file, which is the divergence `substrate.py` exists to
prevent. So `g_adapt` now has **no default at all** — every caller states its intent — and
`None` means "keep the cell class's own value". `tau_adapt` likewise; it happened to equal
`Pop`'s default, so it hid behind the same mistake without contributing to it.

**Then I ran the check I should have run before claiming the LIF path was preserved.** I had
verified `resp.py` byte-identical and asserted the same for `circuit.py`/`rg2.py` without
testing it. A git worktree at HEAD, 28 trace hashes over `rg2` and `circuit` (control and BZ
arms, every population, 2 s at dt = 0.25 ms): **all 28 identical**. That is the check that
would have caught E20 in the first place, and it costs one worktree.

Pinned as two tests: `make_pop` must have no `g_adapt` default, and every spinal population
must carry `Pop`'s value. Suite: **262 passed, 0 failed**.

## Reviewing `initial_paper.md` against the code it claims to describe

Full review in `knowledge/07-paper-review.md`. Regenerator in `scripts/paper_numbers.py`.

**The algebraic half reproduces exactly.** Every nominal R (10.8750, 8.6442, 7.8750, 1.2083,
1.0000, 0.5633, 0.0000), the stereoisomer pair (9.3487 / 3.6250), every Table 6A input vector,
all three Monte Carlo floor scenarios to the digit (5th percentiles 2.51/2.49/2.57; alogabat
median 16.05→10.79→7.57, 95th 84.85→23.42→14.31), the 38% tail, and the declared
Python/NumPy/SciPy versions and seed. §4.2–4.5 is publishable as it stands, and its Dirichlet
dispersion derivation is better than the repo's own documentation — it found the 175% CV on
preBötC α5 that explains the 38% tail `ranking_robustness.py` reports but never explained.

**The kinetic half quotes numbers no commit of this repository produces.** The draft calibrated
to `Po_max = 0.84` against a 3.0 mM / τ_clear = 1.0 ms transient; the repo says **0.75** and
**1.0 mM / 0.30 ms**. `git log -S` finds neither 0.84 nor 3000 anywhere in history,
`FIT_RANGES` declares po_max acceptable only on **[0.70, 0.80]** so 0.84 would have been
**rejected by the module's own validation**, and the commit hash the draft cites twice
(`c8f17a9`) **exists in no ref**.

Not a rounding complaint. The thesis is the contrast between a saturated synapse and an open
extrasynaptic space, and the synaptic half — the part the abstract leads with — moved most:

| | draft | repo |
|---|---|---|
| phasic peak gain | 1.062× | **1.319×** |
| phasic max headroom | 1.124× | **1.668×** |
| charge ratio | 1.655× | **2.419×** |
| tonic max headroom | 210.7× | **184.6×** |
| the headline contrast | 187× | **111×** |

The qualitative conclusion survives. Every quoted figure does not.

**The worst finding is in §5.2, the falsification test.** The draft pre-registers four intervals
for a wet lab to run. Two are violated by the model itself:

| `[GABA]_bath` | draft | repo (s_max 2.40–2.50) | asymptote |
|---|---|---|---|
| 0.1 µM | **> 15×** | **7.66–8.47×** | 2881× |
| 0.4 µM | [7, 15] | 7.18–7.88× | 184.6× |
| 3.0 µM | **< 3.0×** | **2.94–3.03×** | 4.8× |
| 10.0 µM | < 1.5× | 1.34–1.35× | 1.5× |

A lab measuring 8× at 0.1 µM would have reported the model falsified when the model predicts
8×. The intervals read like asymptote values, and the cause is structural: **as ambient GABA
falls, asymptotic headroom grows without limit while the gain a finite s_max can reach barely
moves.** From 0.4 to 0.1 µM the asymptote rises 15.6× and the reachable gain 1.08×. Headroom
and reachable gain are different quantities and the draft's own Table 4 keeps them apart —
then §5.2 pre-registers against the wrong one. Pinned as a test.

**8 of 20 references are defective**, four carrying load-bearing parameters: Walters 2000 (the
diazepam s_max = 2.50 that is the index's denominator — actually Nat. Neurosci. 3:1274–1281, not
Br. J. Pharmacol. 131:1307–1314), Haas & Macdonald 1999 (the kinetic topology — J. Physiol.
514:27–45, not J. Neurosci. 19:2435–2445), Kasugai 2010 (`EXTRASYN['a5'] = 0.80`, the whole
headroom argument — Eur. J. Neurosci. 32:1868–1888, not J. Neurosci. 30:14024–14035), and three
**unresolvable**: Nutt 2007, Fischer 2010 for MP-III-022 (anachronistic — the earliest primary
MP-III-022 literature Crossref indexes is 2024–2026), and Saba 2017 (which carries
`w_subj,δ = 0.0`, the assignment producing gaboxadol's R = 0.00). Crestani is cited to the wrong
year *and* pages; Otis & Mody to a non-existent coordinate.

The Kasugai error is already on record: `provenance.report()` says *"a5_dist resolves perfectly
by DOI and does not support the numbers it is the obvious candidate for."* The draft reached the
same source independently and cited it to the wrong journal.

**The draft also quotes the two statistics the model refuses to stand behind.**
`ranking_robustness.py` prints "NOT QUOTABLE: the median and 95th percentile"; the draft's
Table 6 makes both of them columns. §4.4.3 then explains the prior-sensitivity correctly — so
the draft knows — but a reader who stops at Table 6 has taken numbers the model disowns.

**And it overstates provenance.** `provenance.report()` says 6 of 28 parameters name a source
(21%). The draft attributes Table 5 to Pirker 2000 and Kasugai 2010; `prebotc.a5`, four of five
forebrain fractions, `EXTRASYN.a1`, `EXTRASYN.a23`, four of five `w_subj` weights and `KAPPA`
are all **UNSOURCED** in the audit.

**§4.4 is titled "Multiscale Circuit Selectivity Pipeline" and is not multiscale.**
`ranking_robustness.py` imports `circuitpharm.subtypes` and nothing else — no neuron, no
circuit, no simulation; its own docstring says so. The draft's §4.4.1 subheading already says
"Algebraic Pipeline", so the section heading contradicts its own first line.

**What the draft is missing is this project's best result.** It contains no circuit simulation
at all, and so omits link 5: **Spearman +1.0000 over the five subtype-selective arms** across
an integrate-and-fire cell and a Butera–Rinzel–Smith conductance cell that share almost no
mechanism. That perturbs the model's *structure*, not its parameters, and is a stronger
robustness argument than the Dirichlet sweep. It also omits the finding that cuts against the
draft's own reference arm: the non-selective BZ's sign on respiratory output **flips**
(+0.079 → −0.377), and diazepam is the denominator of every R in Table 6.

## The pattern, stated once

E20 and the manuscript are the same failure at two scales. A number that lives in two places
diverges, and the symptom is plausible rather than loud: a spinal circuit that still runs
without adaptation, a paper whose numbers are internally consistent and match no commit. The
answer in both cases was to delete the second copy — `g_adapt` has no default, and the
manuscript's tables are printed by `scripts/paper_numbers.py` rather than typed.

## Regenerating the manuscript on the 0.75 anchor

Decision: **keep `FIT_TARGETS['po_max'] = 0.75`** and bring the manuscript to the code, rather
than moving the code to the draft's 0.84. 0.84 sits outside the module's own declared
`FIT_RANGES` of [0.70, 0.80], so adopting it would have meant widening a validation range to
accommodate a number that arrived from outside the repository.

Output: `knowledge/08-manuscript.md`, with every figure emitted by `scripts/paper_numbers.py`
at tag **`manuscript-v2`** (commit `1c0faa5`) — a tag created for the purpose, since the draft's
`c8f17a9` resolves to nothing.

### A conflation the review had not caught

The draft used one number for two quantities, and so did my review of it. `po_max()` in the code
is `po_peak(1e4)` — the **simulated** peak during a 300 ms saturating step, with desensitisation
competing throughout the rise — which is **0.7500**. The **analytic** gating bound β/(α+β) is
**0.8282**. The draft wrote "P_o,max = β/(α+β) ≈ 0.84" and used 0.8382 for both. They differ by
10%, the analytic bound is never attained, and Table 1 of the regenerated manuscript now lists
them as separate rows with the distinction spelled out. Related: the analytic equilibrium
midpoint [G]₁/₂ = 18.43 µM against the simulated peak EC₅₀ = 20.00 µM, an 8% gap of the same kind.

### What changed, beyond the numbers

* **§5.2 falsification rules rewritten against the reachable column**, with the asymptote
  printed beside it so the two cannot be confused again. Added two rules the draft lacked, both
  stronger than the original because they test the mechanism rather than the fit: R_PAM must
  decline monotonically with ambient GABA (follows from the topology, unrescuable by refitting),
  and R_PAM at 0.1 µM must not exceed its 0.4 µM value by more than ~1.2× or a gating component
  should be suspected.
* **Table 6 drops the median and 95th-percentile columns.** They now appear only in §4.4.3 as a
  prior-sensitivity diagnostic, with the five-fold swing (alogabat 95th: 84.85 → 14.31 on one
  prior assumption) stated as the reason.
* **§4.4 retitled "Algebraic selectivity index"**, with a scope note that the script imports
  `circuitpharm.subtypes` and nothing else.
* **Table 5 carries a per-cell provenance column** generated from `provenance.py`: 7 UNSOURCED,
  1 GUESS, 7 FROM_QUALITATIVE, **0 QUANTITATIVE** over 15 cells.
* **§4.6 added** — substrate independence (Spearman +1.0000) and the BZ sign flip, including the
  retracted frequency mechanism and the burst-detector limitation.
* **§5.3 expanded** to four categories, with the A4 riluzole failure, the unpublished coupling
  weights, the in vitro band, and the reference-arm instability all declared.
* **References: 5 corrected, 3 withdrawn.** The withdrawn ones are listed in a dedicated section
  with the parameters they carried now labelled UNSOURCED, rather than quietly substituting
  plausible replacements. MP-III-022's real primary source turned out to be Stamenić et al.
  (2016) *Eur. J. Pharmacol.* 791:433–443, ten years later than the draft's citation.
* ρ = 6.00 is now reported honestly as sitting *at* the peak-based ratio (5.97) rather than
  between peak and charge, which is what it did under the draft's parameterisation (6.79/4.36).

### My own error, caught by my own test

`tests/test_manuscript_consistency.py` greps the manuscript for the decimal representation of
each figure, computed fresh. Ten checks. Two failed on first run and both were mine:

1. An over-naive assertion that `"R_max > 15"` never appears — but §5.2 legitimately **quotes**
   the discredited criterion while retracting it. Rewritten to require every occurrence to sit
   within 400 characters of the retraction, and to be absent from the pre-registered rules block.
2. A provenance assertion written from the manuscript's own prose ("nine of fifteen UNSOURCED")
   rather than from the audit. **The prose was wrong** — it is 7 UNSOURCED, 1 GUESS,
   7 FROM_QUALITATIVE. A test that trusts the document it checks checks nothing, so the test now
   derives every expected count from `REGIONS_PROV`/`EXTRASYN_PROV` and asserts the table
   reprints them exactly. Fixed three places in the manuscript where I had written "nine".

That second one is the same failure I spent this session documenting, committed inside the test
written to prevent it. Worth recording for that reason alone.

Suite: **276 passed, 0 failed.**

### Still outstanding

* Replace the UNSOURCED `REGIONS`/`EXTRASYN`/`SUBJECTIVE_WEIGHT` point values with sourced
  ranges and re-run the robustness analysis. This is the highest-value remaining work, and §5.3
  says so.
* `KAPPA = 15.0` governs the whole uncertainty model and has no source.
* Read `pbc_alpha` (403) and `a5_disc` properly; resolve the 12 UNRESOLVED sources in the
  registry; `ganaxolone`'s mis-citation.
* Burst-level detector for the conductance substrate, to settle whether the BZ sign flip is
  more bursts or longer bursts.
* Anchor `COND_LOCO_OP` and `COND_SPINAL_OP`.
