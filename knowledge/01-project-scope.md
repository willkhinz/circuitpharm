# Scope and how the project got here

## In scope

Reproduce alcohol's **acute subjective effect** while eliminating its **physiological**
harms: hepatotoxicity, carcinogenicity, nutritional/thiamine damage, overdose death.

Those four are genuinely separable from the receptor pharmacology because they live in
**metabolism and PK**, not at the receptor. Ethanol -> acetaldehyde (ADH) -> acetate
(ALDH), CYP2E1 induction, NAD+/NADH redox shift causing steatosis, and acetaldehyde as an
IARC Group 1 carcinogen. A non-ethanol compound metabolized differently avoids all of it.
Overdose is addressable via a **ceiling effect**.

## Out of scope — by user decision, not oversight

Dependence, withdrawal, accidents, violence. Accepted as consequences of any intoxicant
in a society that drinks. Do not raise these as objections to the project; they were
considered and deliberately excluded.

Note the honest framing this forces: the claim is **not** "alcohol without harm." It is
"alcohol without the hepatotoxicity, carcinogenicity, and overdose death, with dependence
liability that is hopefully lower and must be measured." Any GABA-A positive modulator
will produce receptor adaptation and a withdrawal syndrome with chronic use — that is the
observed behaviour of every drug class in this niche (barbiturates, benzodiazepines,
Z-drugs, GHB, all cross-tolerant with alcohol). Benzodiazepines were explicitly marketed
as the safer substitute for alcohol and barbiturates.

## Also out of scope — different toolchain

Hepatotoxicity and carcinogenicity screening are **not** circuit problems and cannot be
addressed by any model in this repo. They need CYP profiling, reactive-metabolite
screening, hepatocyte assays, DILI models. Route them elsewhere. This is where the
top-four harms actually get eliminated, and it is where the real win lives — the
simulation work here is about not *introducing* new harms while doing it.

## Prior art — read before anything else

David Nutt (Imperial, ex-chair UK ACMD) has pursued GABA-A subtype-selective partial
agonists as alcohol substitutes for roughly a decade, via Alcarelle / Alcohol Research
Ltd, with a botanical precursor product already sold. A decade of accumulated obstacles
is documented there.

## Regulatory reality

A novel psychoactive intended for recreational use fits no existing category: not a
medicine (no disease indication), not a food additive (psychoactive), not alcohol. The
grey-market "alcosynth" episode made regulators actively hostile. If the binding
constraint on this project turns out to be regulatory rather than scientific, simulation
is not the first step — it is the step taken *after* knowing what evidence a regulator
would need.

## One terminology note

Early in this project's history the phrase "mirror compounds" was briefly ambiguous
between enantiomer pairs (in scope, and the thread that led here) and mirror-image
*biology* — D-amino-acid proteins, mirror life. The latter was explicitly ruled out by the
user. It is not part of this work and should not be reintroduced.
