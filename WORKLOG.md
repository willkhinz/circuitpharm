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
