# Reasoning trail — corrections, dead ends, and why

The most valuable file here. Each entry is something that was believed, then revised.
Reading this prevents re-deriving dead ends that cost real effort.

## Abandoned: the virtual rodent's ANN as a pharmacological substrate

**Believed:** apply drug-induced f-I curve shifts to the policy network's activation
functions (gain/bias), then run behaviour.

**Why it failed:** the policy has no cell types, no E/I partition, no Dale's law, no
laminar structure, no region labels. There is no referent for "apply the D2-MSN shift to
the striatal subpopulation." The Nature 2024 result is *representational* (network
activations predict recorded activity via fitted decoders), which does not license
"units correspond to neurons, so perturbing units ~ perturbing neurons." Also: gain/bias
have no biophysical units, and changing them moves the policy off its training manifold,
so any perturbation produces generic degradation that looks vaguely drug-like.

**Replaced by:** apply pharmacology only where biology is grounded (spinal CPG, preBotC).
The ANN stays untouched as a descending-command generator — a claim the Nature result
does support.

## Abandoned: EEG / neural-signature matching as the subjective proxy

**Believed:** match ethanol's network signature (spectra, firing rates) as a proxy for
"feels like", since no animal can report subjective state.

**Why it failed:** the spectral signature is *anti-correlated* with subjective
substitution. Gaboxadol reproduces ethanol's distinguishing delta increase but does NOT
substitute for ethanol in trained rats. Benzodiazepines diverge on EEG (beta, not delta)
but DO substitute. NMDA antagonists diverge (gamma) but DO substitute. Opposite rank
ordering on the three best-characterized classes.

**Replaced by:** drug discrimination as the subjective proxy, taken as an empirical
constraint rather than a simulated output.

## Over-read: the 2025 subtype-comparison negative result

**Believed:** "none of five GABA-A subtype-selective compounds matched ethanol" falsifies
the GABA-A route (and by extension the Nutt premise).

**Correction:** four of the five (alprazolam, zolpidem, KRM-II-81, MP-III-022) act at the
**benzodiazepine site**, which requires a gamma subunit — and gamma-containing GABA-A
receptors are not ethanol targets below ~100 mM, while 20 mM is strong intoxication. The
study tested compounds against a receptor ethanol barely engages. The negative result is
mundane, not deep.

**Also:** the fifth compound, gaboxadol, was reported as having *no* spectral effect at
1-10 mg/kg. But THIP at 4 and 6 mg/kg produces a dramatic sub-6-Hz power surge that
**dissipates within 60-120 min**. Effective doses sit inside their range, so a long
averaging window plausibly diluted a transient effect to nothing. Do not trust that arm.

## Wrong: "NMDA antagonism is the missing ingredient" (mechanism, not ingredient)

**Believed:** the ethanol-vs-benzodiazepine gap is explained by ethanol's NMDA antagonism.

**Correction:** quantitatively, recombinant GluN2B inhibition is only **6.5 +/- 0.8% at
20 mM**, and even 80 mM is not maximal. A 6% receptor effect cannot flip a cortical
signature.

**But the ingredient was right:** drug discrimination confirms NMDA antagonism IS a real
component of ethanol's subjective stimulus. The reconciliation is an amplification
gradient: ~6.5% at the recombinant receptor -> ~50% of NMDA-mediated synaptic
transmission in BNST (GluN2B-dependent) -> a full behavioural component. **That gradient
is exactly what a circuit model exists to explain** and is an argument for keeping the
simulation.

## Wrong: "weight the mixture NMDA-heavy"

**Believed:** since NMDA has the selectivity window (GluN2B forebrain / GluN2D brainstem)
and GABA-A has none, shift the subjective burden to the NMDA component.

**Correction:** the discrimination data forbids it. Ethanol substituted only for mixtures
where the GABA-A component had **equal or greater salience** than the NMDA component.
The documented failure case was NMDA-weighted (pentobarbitone 8.0 + dizocilpine 0.08 ->
33%), and the fix was to *raise* GABA and *lower* NMDA (12 + 0.04 -> 75%).

**Consequence:** this is now the project's binding constraint and its central tension.
See KNOWLEDGE.md.

## Resolved, favourably: the acetaldehyde threat

**Worry:** acetaldehyde is both the carcinogenic/hepatotoxic metabolite to be eliminated
AND a driver of VTA dopamine activity — so removing the harm might remove the reward.

**Resolution:** every acetaldehyde-blocking result is on a **reinforcement** measure —
preference, acquisition of self-administration (not maintenance), relapse/ADE,
conditioned place preference, sensitization. Brain acetaldehyde during ordinary
intoxication is low micromolar and critics argue too low for direct effect. Reinforcement
is already out of scope, so the premise survives — and a substitute being *less*
reinforcing is a feature. **Caveat:** D-penicillamine also blocked acute behavioural
stimulation, so some acute component is acetaldehyde-dependent.

**Side effect:** this demoted VTA dopamine firing, which had been promoted to primary
observable one step earlier. Ethanol's pVTA excitation requires local catalase-derived
acetaldehyde condensing with dopamine to salsolinol, acting via **mu-opioid receptors**.
A non-metabolized substitute cannot reproduce it without mu-opioid agonism, which
reintroduces the exact harms being removed. Decision: accept divergence there.

## Parked side-quest (real, not on the critical path)

(R)-salsolinol **stereospecifically** induces behavioural sensitization; salsolinol forms
as racemic R/S by non-enzymatic Pictet-Spengler condensation of dopamine with
acetaldehyde. A chiral pair with activity confined to one enantiomer, sitting at the
centre of alcohol's own reward mechanism. This is a genuine open research question.

## Methodological principle that survived everything

**Target dimensionless contrasts, never absolutes.** Running both arms of a comparison
through the same pipeline with the same wrong parameters cancels shared systematic error
to first order. The pipeline can be badly miscalibrated and still get the *difference*
approximately right. Applies to enantiomer pairs, to the chosen validation anchors, and
to the ratio optimization. Corollary: prefer ordinal/categorical predictions, which
survive multiplicative monotone error.
