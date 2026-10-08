# circuitpharm

Receptor-resolved pharmacology coupled to spiking neural circuit models, with calibration
and identifiability analysis — and with reliability tiers enforced in code.

## What it is for

Given a compound's **measured** receptor activity, propagate it to circuit-level and
motor/autonomic output, and report every number with an explicit statement of how much it
can be trusted.

```
measured subtype efficacies        <- INPUT (an experiment, not a prediction)
  |  gabaa_kinetics   Markov gating scheme: one calibrated parameter, four predictions
  v
tonic / phasic conductance changes
  |  subtypes         regional subunit composition and synaptic localisation
  v
per-region drug sensitivity
  |  cpg, resp, circuit   spiking networks (LIF, conductance synapses)
  v
circuit dynamics
  |  plant            MuJoCo muscle body (optional extra)
  v
motor / autonomic output
```

## The distinguishing feature: it refuses

Most in-silico pharmacology tools emit a point estimate for anything you ask. This one
attaches a tier to every quantity and **raises** rather than returning a number it cannot
justify:

```python
from circuitpharm import calibration_report
print(calibration_report())        # what is anchored, what is not, and what would fix it
```

```python
q = result.quantity("overdose_index")
q.value          # -> VoidQuantityError, with the reason and the promotion path
q.get(acknowledge_void=True)       # debugging access, explicit at the call site
```

Three tiers, defined by evidence rather than precision:

| tier | meaning |
|---|---|
| `VALIDATED` | reproduces an observable it was not fitted to, or is an exact identity, or survives propagation of every parameter it depends on |
| `UNCALIBRATED` | mechanism sound, no valid quantitative anchor — trust the ordering, not the scale |
| `VOID` | rests on something known to be invalid; not quotable even with a caveat |

This exists because the project's worst reporting error was not a modelling mistake: a
configuration the model itself had flagged `UNREACHABLE` was summarised as a clean result,
because the flag was *printed next to* the number instead of blocking it. Caveats get
dropped when numbers are copied.

## Current state, honestly

**The respiratory axis has no valid quantitative anchor.** Its single calibration (fitted
to human midazolam minute ventilation) was invalidated: whole-body ventilation is the wrong
observable for an isolated preBötzinger complex, because much of benzodiazepine respiratory
depression is chemoreflex blunting plus upper-airway and cortical drive — none of which
exist in this circuit. The parameter was a lumped factor absorbing mechanisms the model does
not contain, so it cannot be extrapolated to a compound with a different subtype profile.
It is additionally unidentified: 18 of 35 grid cells fit the anchor equally well.

**What is validated** is the calibration-independent selectivity *ranking*: subjective drive
per unit respiratory burden is a ratio, so the unknown lumped scale cancels. It survives
20,000 draws over every estimated parameter. It **ranks** candidates; it **bounds** nothing
— no absolute margin, no overdose multiple, no dose.

## What it deliberately does not do

Two links from chemical structure to behaviour are not closed, and the package declines to
fake either.

**structure → efficacy.** For an allosteric modulator, efficacy is not affinity: it is the
shift in the conformational equilibrium between shut, open and desensitised states.
Free-energy methods reproduce affinities and Markov state models have captured a putative
GABA-A open state, but there is no predictive efficacy calculation — the field still
publishes experimental *screening designs* for exactly this quantity, and measured efficacy
is probe-dependent (it varies with the agonist concentration assayed). Efficacy is an input
here, by design.

**circuit → arbitrary behaviour.** Controllers that produce naturalistic rodent behaviour
are artificial networks whose units have no conductances, so there is nothing to inject a
receptor current into. Fully biophysical closed loops exist (BAAIWorm) but in a 302-neuron
animal. This package therefore targets behaviours whose generating circuits *are*
biophysically modellable — breathing, spinal reflexes, locomotor rhythm — which map onto
standard assays (plethysmography, H-reflex, gait analysis).

Subjective or affective state is outside the architecture entirely: drug discrimination
requires an animal that learns, and no increase in resolution crosses that boundary. Also
not modelled: hepatotoxicity, hERG/QTc, dependence, and **any pharmacokinetics** — "dose"
means receptor-occupancy multiple, not mg/kg.

## Worked example

```bash
python examples/quickstart.py
```

It walks the five steps the tool is built around: ask what is actually anchored, describe a
compound by **measured** receptor activity, evaluate, see the refusal fire when you ask for
something unjustifiable, and then use the one quantity that is calibration-independent.

## Install

```bash
pip install -e ".[dev]"            # core + tests
pip install -e ".[all]"            # + MuJoCo body plant and RDKit chemistry tools
```

Requires Python 3.11 or 3.12.

## Tests

```bash
pytest                             # full suite, parallel by default (~2.5 min)
pytest -m "not slow"               # skip circuit integrations (~15 s)
pytest -n 0                        # serial, for debugging (xdist hides stdout)
```

97 tests. The suite is dominated by circuit integration and is embarrassingly parallel —
every test builds its own network with its own seed and shares no state — so it runs on
8 workers by default (2m35s versus 7m16s serial).

It is in five parts, and the second and fifth are the point:

- `test_identities.py` — exact identities (the tonic/phasic split is a partition; legacy
  single-pool behaviour is reproduced exactly; the GluN2B selectivity bound; each region's
  calibration anchor)
- `test_failure_modes.py` — **the catalogued failure modes E1–E12 as executable
  regressions.** Every one of these errors was made during development, several twice, and
  every one produced a result that *looked fine* — a flat sweep, a plausible percentage, a
  reassuring safety margin. Documenting them was not enough. Tests are two-sided where
  possible (the mechanism is engaged, *and* it would disengage if the error returned).
- `test_tiers.py` — the refusal machinery, including that VOID values never leak into
  printed output
- `test_evaluate.py` — the public API: occupancy linearity, scale-invariance of the
  ranking, and the invariance the control cache rests on
- `test_phenotypes.py` — **the only tests that check biology rather than code.** Published
  results the circuits were not fitted to: strychnine hyperreflexia (138% of control),
  benzodiazepine reflex depression (77%), the GluN2B selectivity window (84% non-selective
  versus 102% selective at the *same* 60% block), adaptation-not-inhibition
  rhythmogenesis, closed-loop walking, and sedative-induced coordination loss

## Repository layout

```
src/circuitpharm/
  results.py          reliability tiers (Quantity, ResultSet, VoidQuantityError)
  config.py           single source for operating points + the calibration registry
  cpg.py              LIF populations, conductance synapses, Drug, spindle afferents
  gabaa_kinetics.py   GABA-A Markov gating scheme (5 states)
  subtypes.py         regional subunit composition, synaptic/extrasynaptic localisation
  resp.py             preBötzinger complex, respiratory metrics
  rg2.py              locomotor rhythm generator (group pacemaker)
  circuit.py          spinal reflex arc
  plant.py            MuJoCo closed loop
  chirality.py        stereochemistry tools, parity-invariance harness
tests/                see above
scripts/              analysis and calibration runs (see WORKLOG.md for which are current)
WORKLOG.md            append-only development record, including the E1-E12 catalogue
KNOWLEDGE.md          project state, decisions, and the verdict on what is achievable
```

## Reading the development record

`WORKLOG.md` is append-only and conclusions are never deleted — where one was later
overturned, a `CORRECTION` entry follows and the original stands. Several headline results
in this project were corrected, twice in one case, and the trail is the useful part.

## Licence

MIT.
