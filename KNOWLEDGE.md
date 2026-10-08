# Alcohol-Substitute Project — Knowledge Base

**Entry point for any agent picking this up.** Read this file first, then
`knowledge/02-reasoning-trail.md` (it will stop you re-deriving several dead ends).

---

## ROADMAP — the chain from compound to behaviour, and where it stands

Recorded here because it previously existed only in conversation. This project has already
been bitten twice by results living outside the repo: the tuned locomotor RG parameters sat
in a past run's stdout while the module shipped untuned defaults, and the destroyed KB rows
were recoverable only by accident of version control. A plan in chat is the same failure.

| # | link | status |
|---|---|---|
| 1 | structure → binding affinity | **not attempted.** Achievable with ~1 log unit error; worse for subtype selectivity, which is the hardest regime |
| 2 | affinity → **functional efficacy** | **DECLINED BY DESIGN.** Not predictable from structure; taken as a measured INPUT (`Compound(a1=…, a5=…, s_max=…)`) |
| 3 | efficacy → conductance change | **BUILT** (`gabaa_kinetics.py`). Calibration BLOCKED on wet-lab data |
| 4 | conductance → neuron excitability | **BUILT** (`neuron.py`). Butera–Rinzel–Smith 1999 model 1, conductance-based, real I_NaP. Characterised in isolation; no circuit uses it yet |
| 5 | neurons → circuit dynamics | **NOT DONE.** Weights are "hand-tuned to produce alternating rhythm, NOT fitted to rat data". Unblocked by link 4; `ParamSet` now refuses a mixed cell/weight pair |
| 6 | circuit → behaviour | **BUILT for motor** — stretch reflex and closed-loop locomotion |

### Link 4 is built; link 5 is the remaining half
Link 4 gated link 5. The LIF neurons structurally cannot carry voltage-gated mechanisms
(**E7**: V is reset at threshold, so a slow voltage-gated gate equilibrates at the
subthreshold mean and never moves — measured h-gate span 0.01). That forced spike-triggered
adaptation as a lumped stand-in for I_NaP inactivation *plus* calcium-dependent potassium
current *plus* synaptic depression — defensible, documented, and not the biophysics. It also
blocked link 5, because the published models worth adopting (Butera–Rinzel–Smith,
Rybak-style locomotor CPGs) are conductance-based with a real I_NaP and cannot drop into a
LIF substrate.

**`src/circuitpharm/neuron.py` now holds a conductance-based cell** — Butera–Rinzel–Smith
1999 model 1, three state variables (V, n, h; `m` and `mp` instantaneous, Na inactivation
via `1−n`), Rush–Larsen on the gates, 21 acceptance tests. Design and full rationale:
`knowledge/05-design-conductance-substrate.md`.

What it demonstrably does that the LIF cannot:

* **V reaches +6.49 mV** and the h gate spans **0.4544** against the LIF's 0.01 — E7's test,
  inverted.
* **Freezing h abolishes bursting in every direction**: quiescent at h=0.46, tonic at 68.8
  and 118.9 Hz at h=0.60 and 0.90. So the rhythm comes from slow voltage-dependent I_NaP
  inactivation, not from something incidental. The LIF can neither pass nor fail this — it
  has no h gate to clamp. *This replaced a pre-committed h-span threshold of 0.5, which the
  cell misses at 0.454; that 0.7/0.5 figure was never sourced, and the functional test is
  what it was a proxy for. Recorded in WORKLOG rather than quietly adjusted.*
* **Bursts with adaptation entirely off**, and `g_NaP → 0` makes the isolated cell quiescent.
* **Reproduces the published excitability sequence** quiescent → bursting → tonic in drive,
  which was not fitted by us.
* **Mg²⁺ relief spans 11.80×** against the LIF's 4.50× — the state-dependence the whole NMDA
  arm rests on.

Reading the parameters off a machine-readable source rather than reciting them caught a real
error: **E_L is −57.5 mV, not −65 mV.** E_L is this model's bifurcation parameter, and the
recalled value gives a **completely silent cell** — zero spikes. The bursting band is only a
few mV wide (−60 and −55 are both tonic; only −57.5 bursts), so the rhythm would then have
been retuned around a wrong resting drive to compensate.

**Measured cost: 4.1× the LIF**, not the estimated ~20× — loose enough that the conductance
arm of the substrate-independence comparison may not need reduced draws after all.

**Still to do (link 5):** no circuit uses the cell yet. `resp.py`, `rg2.py` and `circuit.py`
construct `cpg.Pop`, and they assign `g_adapt` *after* construction, which would switch
adaptation back on and re-introduce E18 the moment a circuit is pointed at `CondPop`. The
substrate-independence comparison — the actual deliverable — needs the network.

### Why the upgrade is not cosmetic — measured, not argued
`scripts/diag_substrate_limits.py`. The LIF preBötC holds V in [−71, −44] mV, and:

* **NMDA conductance sits in near-permanent Mg²⁺ block** — mean relief **0.063**, dynamic
  range **4.5x** against the 15.5x a spiking cell has. `cpg.py` says "voltage dependence is
  why NMDA block is state-dependent"; on this substrate it is not.
* **The phasic GABA-A pool's driving force is truncated exactly when it matters.** Mean
  GABA-A driving force is **10.4 mV**; during a real burst it would be up to **95 mV**. Both
  pools share one V, so at rest both are right — the asymmetry is in TIMING. Tonic is always
  on at the resting driving force; phasic arrives correlated with the burst.
* **Excitation never self-limits**: glutamate driving force never approaches zero.

The standing objection is that the weights were tuned on this substrate so magnitudes were
compensated. True — and that is why the claim is about voltage-**dependence**, not
magnitude. A scalar weight rescales a mean; it cannot turn 4.5x into 15.5x, nor make one
pool's driving force swing 9x while the other's stays flat, because both share a single
weight-independent `E_rev` and a single V.

So the exposed result is the project's central one. The tonic/phasic ratio is a conductance
claim and stands; its *behavioural* consequence runs through g·(E−V), which this substrate
truncates for the phasic pool. **The deliverable of links 4–5 is therefore not "a better
neuron" — it is whether the selectivity ranking is substrate-independent as well as
calibration-independent.** Nobody has checked. Either answer is worth having.

### But sequencing matters more than the next feature right now
Eight review passes have found **49 defects**, at a rate of 9, 11, 8, 10, 8, 3. Each review
found defects created or left incomplete by the previous round's fixes. My own scrutiny
across thirteen sessions found none of them.

Two findings make the point better than the count does:

* **Pass 7** found a *second instance* of a defect class this project had already been
  bitten by and already "fixed" — a destructive database rebuild — where the guard added
  the first time was counting the wrong tables and so could not see the collateral it was
  destroying (including 12 rows with no builder anywhere in the repo, which would have been
  unrecoverable).
* **Pass 10** found a NaN latching the grid search in all three calibration sweeps — the
  scripts the project's calibration constants came from — and in one of them **my own
  pass-7 fix is what made the path reachable**, by converting a loud `ZeroDivisionError`
  into a silent NaN. That is the second time in this project I have replaced a loud failure
  with a quiet wrong answer.

The rate is falling, but pass 10 was the last of this review round, and two of its three
findings were in classes pass 7 had already raised — so the fall is partly the reviewer
running out of new surface rather than the code becoming clean. **The gate asks for a pass
that finds nothing. Three is not zero.**

A conductance-based rewrite of the neuron model plus a circuit swap is a large new surface
on a base whose defect density I cannot measure. **Get one review pass that finds nothing
first.** Until that happens, more model is more unvalidated model, and the flat findings
rate is better evidence about this codebase than any amount of added mechanism.

---

## THE TOOL (session 8) — what exists now

The code is packaged as **`circuitpharm`** and is the durable output of this project,
independent of whether the alcohol-substitute question is ever answered.

```bash
pip install -e ".[dev]"      # core + tests;  ".[all]" adds the body plant and RDKit
pytest                       # 114 tests, parallel, ~2.5 min
python examples/quickstart.py
```

**What it does.** Takes a compound's MEASURED receptor activity and propagates it to three
behavioural endpoints, attaching a reliability tier to every number:

| endpoint | readout | maps to |
|---|---|---|
| respiratory | preBötC rhythm, ventilation proxy | plethysmography |
| spinal reflex | ramp-and-hold stretch reflex gain | H-reflex / tendon jerk |
| locomotor | free-running closed loop: step period, coordination | gait analysis |

**What makes it unusual.** It refuses. Reading a `VOID` quantity raises rather than
returning a caveated number, because caveats get dropped when numbers are copied — which
is exactly how this project's worst reporting error happened. `python simulator.py
calibration` prints what is anchored, what is not, and what would promote each one.

**Phenotypes it reproduces without having been fitted to them** (the only route to
VALIDATED): strychnine hyperreflexia (138% of control), benzodiazepine reflex depression
(77%), the GluN2B selectivity window (84% non-selective vs 102% selective at the *same*
60% block), the lamprey strychnine phenotype in the locomotor RG (alternation lost, rhythm
intact, only under COMPLETE glycine removal), adaptation-not-inhibition rhythmogenesis, and
the parity-invariance identity exactly (0.00e+00).

**Failure modes E1–E12 are executable regression tests**, not prose. Every one of them was
made during development, several twice, and every one produced a result that looked fine.
Doing session 8's engineering found three more of the same kind — see the worklog.

---

## VERDICT (2026-10-07, after sessions 7a-7f) — READ BEFORE DOING ANY MORE MODELLING

**Asked: is there a perfect answer, or does none exist? Answer: none exists by this
method, and the reason is structural rather than a shortage of effort or compute.**

Three limits, each ESTABLISHED this session rather than assumed:

**1. The target endpoint is not a physical observable this method can reach.**
Alcohol's subjective effect is measured by drug discrimination, which requires an animal
that LEARNS. The "subjective index" in this repo is algebra over literature weights, not a
simulation — there is no forebrain circuit, and building one would not fix it, because its
readout would still need calibrating against a learning animal. No increase in resolution
(atomic, molecular-dynamics or otherwise) changes this; it is a category limit, not a
precision limit.

**2. The safety endpoints cannot be bounded in silico.** Calibrating the respiratory axis
against human whole-body ventilation is structurally invalid: `gaba_sens` was absorbing
chemoreflex, upper-airway and cortical mechanisms that an isolated preBotC does not
contain, so it never meant what its docstring claimed and cannot be extrapolated to a novel
subtype profile. The structure-matched anchor that would fix this does not exist in usable
form in the literature (single concentrations, developmentally ambiguous preparations).
A zero-free-parameter replacement prediction was generated instead
(`scripts/predict_muscimol.py`) but it is untestable without new wet-lab work.

**3. The mechanism the entire safety case rested on is a per-molecule MEASURED property.**
The efficacy ceiling is pool-dependent (synaptic headroom ~1.1x, extrasynaptic ~211x); the
hard-coded 2.5 is correct only at an ambient GABA ~10x above physiological; and the
project's two core assumptions — a tight ceiling as the sole source of overdose protection,
AND full subjective delivery by an a5-selective PAM — **are the same parameter and are
mutually exclusive.** Overdose protection must therefore be demonstrated per candidate. It
cannot be inherited from the PAM class or derived from structure.

### What the method CAN do, and this is now proven rather than hoped
**It ranks candidates robustly.** 20,000-draw propagation over every invented parameter —
regional subunit fractions, extrasynaptic fractions (with `eps` sampled uniform on [0,1]
because it was flagged a guess), and the tonic:phasic gain ratio — gives, for subjective
drive per unit preBotC respiratory burden relative to a non-selective benzodiazepine:

| arm | median | 5th pct | P(better than non-selective) |
|---|---|---|---|
| IDEAL a5-only PAM | 40.9 | 2.58 | 99.9% |
| alogabat (a5) | 16.1 | 2.51 | 99.9% |
| MP-III-022 (a5) | 12.6 | 2.49 | 99.9% |
| HZ-166 (a2/a3) | 1.34 | 0.34 | 63.8% |
| neurosteroid | 0.60 | 0.33 | **0.0%** |

a5 beats a2/a3 in 94.8% of draws. Neurosteroid is worse than a plain benzodiazepine in
100% of draws — independently reproducing an earlier retraction.

**The ordering is the deliverable.** It survives adverse draws of every parameter because it
is a RATIO, so the unknown lumped calibration scale cancels. It orders candidates; it does
NOT bound risk — no absolute margin, no overdose multiple, no dose.

### What would be needed to go further (all wet-lab, none of it compute)
1. Muscimol or GABA concentration-response on inspiratory burst frequency AND amplitude in
   a rhythmic preBotC slice or perfused preparation, **in animals older than P12** so KCC2
   dominance is established. Zero free parameters; falsifies or calibrates the respiratory
   axis in one experiment. Prediction already on record: the ratio of failure concentration
   to 50%-suppression concentration should be ~1.5.
2. Per-candidate intrinsic allosteric efficacy (maximum achievable GABA EC50 shift) by
   patch clamp — this IS the overdose-protection parameter.
3. Drug discrimination in trained animals for any subjective claim. Not substitutable.
4. The mandatory external gates that no simulation touches: DILI/hepatotoxicity screening
   (four GABA-A anxiolytics died of it, one with human toxicity invisible in rat, dog and
   monkey), hERG/QTc, and human dose-response.

**Do not spend further effort making the simulator more elaborate.** The dominant
uncertainty is in measurements it cannot generate, and added detail buys apparent authority
without reducing it. The repo's remaining value is the ranking above, the error catalogue
(E1-E12), the parity-invariance harness, and the specified experiments.

## The goal

Find a compound or mixture that reproduces alcohol's acute subjective effect while
eliminating its **physiological** harms: hepatotoxicity, carcinogenicity, thiamine/
nutritional damage, and overdose death.

**Explicitly out of scope** (user decision, not an oversight): dependence, withdrawal,
accidents, violence. These are accepted as consequences of any intoxicant in a society
that drinks. Do not re-open them as objections.

## Project state (2026-10-06)

The mechanism question is **answered** and did not require novel chemistry.
Ethanol's discriminative stimulus is a *compound stimulus* = GABA-A positive modulator +
NMDA antagonist. Mixtures of the two generalize to ethanol in trained rats.

The project is now a **constrained ratio optimization**, not a discovery problem:

```
maximize   M_motor = D_impairment / D_subjective
           M_resp  = D_respiratory / D_subjective
over       (GABA-A occupancy, NMDA occupancy, subunit selectivity)
subject to GABA-A salience >= NMDA salience      [subjective-match constraint]
```

### The central tension — this IS the project

| | Subjective match requires | Safety margin favours |
|---|---|---|
| Ratio | **GABA-A >= NMDA** salience | **NMDA-heavy** |
| Why | Ethanol only substituted for mixtures where the GABA component had equal or greater salience | NMDA has a real regional selectivity window (GluN2B forebrain vs GluN2D brainstem); GABA-A has **none** |

The two constraints pull opposite ways. Resolving this is the whole optimization.
Because GABA-A has no respiratory selectivity window, the GABA arm's safety must come
**entirely from the ceiling property** — a low-efficacy positive allosteric modulator
(requires endogenous GABA, saturates), never a direct agonist. You cannot dose-shift
out of it.

> **CORRECTION 2026-10-07 (session 7c) — read this before relying on the ceiling claim
> above.** An explicit Markov gating model (`circuitpharm/gabaa_kinetics.py`) shows the ceiling
> is not one number and is much weaker than assumed in the pool that matters here:
>
> | pool | agonist seen | headroom for an affinity-type PAM |
> |---|---|---|
> | synaptic | near-saturating cleft transient | **~1.1x** — argument holds |
> | extrasynaptic / tonic | ~0.4 uM vs EC50 ~20 uM | **~211x** — argument fails |
>
> The model's hard-coded `gaba_a_efficacy_cap = 2.5` is only correct at an ambient GABA of
> ~4-5 uM, an order of magnitude above physiological. **α5 is predominantly extrasynaptic,
> so the protection is weakest precisely where the α5 strategy places the drug.**
>
> What survives: a PAM still needs endogenous GABA and does nothing without it, so it
> remains genuinely unlike a direct agonist, and the clinical observation that BZ overdose
> is rarely fatal alone is an observation and stands. What changes is the *mechanism*. If
> the protection comes from the ligand's bounded intrinsic allosteric efficacy (max
> achievable EC50 shift) rather than from agonist saturation, then it is **a property of
> the specific molecule, to be measured per candidate, not a property of the PAM class that
> a novel compound inherits.** For a design project that distinction is the whole point.
>
> Also note (session 7b): at the subjective target used throughout, the α5 arm needs a PAM
> gain of 2.54 against a ceiling of 2.5 — **it cannot deliver the target effect**, and is
> unreachable in 61% of parameter draws. The binding constraint on this project is efficacy
> ceiling versus subjective load, not respiratory safety.

### Decisions already taken — do not silently revisit

1. **EEG-signature matching is abandoned.** The spectral signature is anti-correlated
   with subjective substitution across the three best-characterized classes. See finding 1.
2. **Drug discrimination is the subjective proxy**, taken as an *empirical constraint*
   rather than simulated (it requires learning; a frozen policy cannot model it).
3. **Accept divergence on the reward/VTA axis.** Matching it would require mu-opioid
   agonism, reintroducing respiratory depression and dependence. See findings 5.
4. **Only two things get simulated**: motor impairment and respiratory depression.
5. **Direction chosen**: retarget the simulation (vs. going purely empirical, or pivoting
   to publish the dissociation).

## Primary hypothesis to test

A **GluN2B-selective antagonist** (ifenprodil, traxoprodil/CP-101,606, Ro25-6981,
radiprodil, rislenemdaz) plus a **low-efficacy GABA-A PAM**, at the most NMDA-weighted
ratio that still satisfies the substitution constraint, gives a wider margin than any
alternative. Rationale: GluN2B is ethanol's own NMDA mechanism, is forebrain-restricted
(overlapping the NAc core / CA1 sites where intracranial MK-801 fully substitutes), and
should spare GluN2D-dominant brainstem respiratory circuits. It may also be
non-psychotomimetic for a mechanistic reason worth testing: forebrain PV+/SST+
interneurons are GluN2D-enriched, so GluN2B-selective block spares them, avoiding the
disinhibition that is thought to drive ketamine's gamma surge and dissociation.

## RECOMMENDED COMBINATION (revised 2026-10-07)

**ALOGABAT (RG7816 / RO7017773) — a highly selective GABA-A α5 PAM already in Roche Phase 2
— plus a GluN2B-selective or low-trapping NMDA antagonist.**

Alogabat: ≥80% sustained α5 occupancy in humans, safe and well tolerated, dosed in children
5–17 for 12 weeks, **no hepatotoxicity signal**. Simulated at subjective 0.45: ventilation
101%, reflex 103%, overdose >26×. With the NMDA arm (subjective 0.86): ventilation 93%,
reflex 100%.

Backups on the same target: **GL-II-73** (α5-preferring, no α1 affinity, and NEUROTROPHIC on
chronic dosing — the opposite of alcohol), **MP-III-022**, **KRM-II-81 analogs** (IND-enabling).
Adjunct worth adding: **ASP-8062** (GABA-B PAM, Phase 2 for AUD) — the only candidate with
direct evidence of **no augmentation when co-administered with alcohol**.

| condition (8 seeds, % of drug-free control) | ventilation | stretch reflex |
|---|---|---|
| α5-selective PAM ALONE, subjective 0.50 | **100±0%** | **98±7%** |
| **RECOMMENDED: α5 PAM + GluN2B 20%** | **92±1%** | **96±6%** |
| ref: non-selective benzodiazepine, same subjective effect | 83±1% | 72±5% |
| ref: neurosteroid, same subjective effect | 80±1% | 67±5% |
| superseded recommendation (neurosteroid + GluN2B) | 76±0% | 68±5% |

**Why**: α5 mediates ethanol's discriminative stimulus and is forebrain-restricted; α1
mediates sedation, respiratory depression and ataxia and dominates the medulla. Sparing α1
separates the subjective effect from the harm. Overdose index >26× (alcohol's is ~5×).

**Constraint**: a ceiling-limited α5-selective compound saturates at subjective ≈0.51, so
stronger effects need the NMDA arm or a higher ceiling (forfeiting overdose protection).

**NMDA arm options, best human safety evidence first**: esmethadone (REL-1017; Phase 3, no
abuse potential in recreational users, opioid-inactive d-isomer of methadone), lanicemine
(low-trapping, no dissociation at 150 mg IV), rislenemdaz (GluN2B-selective, safety clean).
Reject memantine and dextromethorphan — PCP-like, not ethanol-like.

**DOMINANT RISK — hepatotoxicity, now partly resolved.** Four GABA-A anxiolytics caused
liver injury: alpidem (withdrawn post-market), ocinaplon (halted), panadiplon (halted; human
toxicity invisible to rat, dog AND monkey), kava (EU bans; alcohol co-use a listed risk
factor). **But all four are NON-benzodiazepine chemotypes** (imidazopyridine,
pyrazolopyrimidine, quinoxalinone, kavalactone), while the benzodiazepine and
imidazobenzodiazepine scaffolds — including every α5-selective candidate above — carry no
signal. Hepatotoxicity appears to track chemotype, not the GABA-A mechanism. A hypothesis
from a small sample: DILI screening remains mandatory, and standard preclinical species
missed panadiplon.

Other mandatory gates no model here can evaluate: **hERG/QTc** (killed traxoprodil) and a
**human dissociation dose-response** (GluN2B block is less dissociative than ketamine, not
non-dissociative).

## Simulator

`python simulator.py --help` — predicts subjective index (with the salience constraint),
ventilation, stretch reflex and overdose index from a receptor-activity profile, and lists
what it cannot predict. Subcommands: `subtypes`, `profile`, `compound`, `refs`, `optimise`.

## Compound survey

46 researched compounds in the `compounds` table: 7 candidates, 7 backups, 8 tools,
22 rejects, 2 gaps. Query with `.venv/bin/python scripts/kb.py compounds` or `cand`.

## Build state (what runs today)

| Component | Status |
|---|---|
| dm_control rodent body | verified: 67 joints, 38 actuators, 0.3378 kg, ~41,700 steps/s (~83x realtime) |
| **Muscle port** | **done** — `models/rodent_muscle.xml`, 24 antagonist muscle actuators replace 12 hindlimb position servos (nu 38→50), Hill force-length verified |
| **Spinal circuit** | **running** — `spinal/` : Matsuoka RG + spiking PF/Mn/interneurons, period 920 ms, Mn 36 Hz, duty 0.33, inside the rat physiological box |
| Drug interface | working, receptor-resolved and subunit-resolved (`spinal.cpg.Drug`) |
| MuJoCo ↔ circuit coupling | **NOT DONE** — all drug results so far are open-loop |
| preBötC respiratory module | not started |
| Validation vs the four motor anchors | **not done** — nothing here is validated yet |

First in-silico results are in the `sim_results` table and findings 13-17. Headline: a
60% GluN2B-selective NMDA block costs only 7% motoneuron output in GluN2D-dominant spinal
tissue vs 43% for a non-selective block at the same occupancy — the margin hypothesis,
quantified in-model. And the PAM + non-selective combination is supra-additive on motor
depression while the PAM + GluN2B-selective combination is ~multiplicative.

**RETRACTED: finding 15** (ethanol lengthens step-cycle duration via GlyR). The
mechanism is architecturally wrong — in the real cord, strychnine abolishes left-right
alternation while the rhythm *persists*, so glycinergic reciprocal inhibition sets PHASE,
not period. The Matsuoka RG conflates the two. Rebuild the locomotor RG on the
group-pacemaker architecture already working in `circuitpharm/resp.py`. Margin results are
unaffected (reflex runs with the RG off; respiratory uses the correct architecture).

**BLOCKING UNKNOWN: preBötC GABA-A sensitivity is uncalibrated** (finding 21). Absolute
margin numbers are not believable until it is fitted against measured benzodiazepine or
propofol dose-vs-ventilation data in rodents. The selective-vs-non-selective *contrast* is
a paired comparison and does survive a common calibration error.

**Do not retry tuning an all-spiking rhythm generator.** A 288-point exhaustive search
proved it cannot work, for a structural reason (finding 17).

## How to query the knowledge base

```bash
.venv/bin/python scripts/kb.py findings          # claims + project bearing
.venv/bin/python scripts/kb.py targets           # ethanol mechanism panel w/ concentrations
.venv/bin/python scripts/kb.py disc              # drug-discrimination substitution data
.venv/bin/python scripts/kb.py subunits          # the selectivity-window evidence
.venv/bin/python scripts/kb.py eeg               # EEG study comparability audit
.venv/bin/python scripts/kb.py model             # locally VERIFIED facts about rodent.xml
.venv/bin/python scripts/kb.py sql "SELECT ..."  # arbitrary query
```

Rebuild from source: `.venv/bin/python scripts/build_kb.py` (idempotent).
Every row carries a `source_key` joining to `sources` (39 citations with URLs).

## Open items, highest value first

1. **Assemble the full admissible-region constraint set.** Only 3 mixture data points are
   in the DB so far (see `discrimination` table). Need the complete set of (GABA, NMDA)
   dose combinations with substitution percentages — that is the box the optimizer searches.
2. **The epsilon-subunit idea is unvalidated and is mine, not the literature's.** It is the
   only potential selectivity lever on the GABA arm. Read the epsilon pharmacology properly
   before building around it.
3. **Aperiodic (1/f) reanalysis of existing ethanol EEG recordings.** No new experiment
   needed; nobody has done it. Cheapest publishable item. (Note: now peripheral to the
   main line, since EEG matching was abandoned.)
4. **Close the MuJoCo ↔ circuit loop.** Mn pool rate → muscle activation; joint/tendon
   length and velocity → `spindle_ia()` → Ia afferents. Until this is done IaIn is silent
   and reciprocal inhibition is untested.
5. **Test finding 15 against existing rodent gait data.** The model predicts ethanol
   lengthens step-cycle duration (via GlyR) while a GABA-A PAM + NMDA antagonist mixture
   reduces burst amplitude *without* slowing the cycle. Cheapest external validation
   available, and it is a genuine falsification test.
6. **Build the preBötC respiratory module** — the safety-critical half, and the
   benzodiazepine/barbiturate ceiling dissociation is its validation anchor.

## Hard constraint on scope

Mechanism- and class-level work only. Specific formulations, ratios-as-product, or doses
intended for human use are regulated drug development and need institutional context.
Earlier in this project's history the phrase "mirror compounds" was briefly ambiguous with
mirror-image *biology* (D-amino-acid proteins, mirror life); that was ruled out by the user
and is not part of this work. Do not reintroduce it.
