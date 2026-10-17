# Finding — WITHDRAWN

**Status: WITHDRAWN by Claude Opus 5, 2026-10-17.** The characterisation below — that the
deactivation failure is an *amplitude* obstruction while the time constants fit exactly — does
not survive its own canonical test. The sections after this one are preserved as written, so
that what was claimed stays legible. Do not cite them.

## Why it is withdrawn

**1. The χ² = 452.99 figure is an arithmetic identity, not a property of any model.**

The "amplitudes only" objective scores fractional amplitudes, which sum to 1 by construction.
Any model that can drive its fastest selected amplitude to zero is then left with amplitudes
`(0, x, 1-x)`, and minimising the published χ² over that one free number gives:

```
chi^2 = 452.9856  at x = 0.4054   ->  amplitudes (0, 0.4054, 0.5946)
   of which the zeroed fast term alone contributes (0.57/0.04)^2 = 203.06
```

The measured value was 452.99, in **every** scheme — `kinetic_jw95_free_rates`,
`kinetic_flip` and `kinetic_open_desens_free_rates` alike, to two decimal places. **It matches
the identity to 0.0044.** The "achievable amplitude triple" reported below as
(0, 0.405, 0.595) is simply that arithmetic optimum.

I cited the invariance across topologies as evidence that the constraint was structural. It is
the opposite: the invariance is what an arithmetic floor looks like. Every scheme hits the same
number because every scheme can zero its fast amplitude, and after that no model structure
enters at all.

**2. Under the canonical fit route the time constants do NOT fit exactly.**
`DECISION-deactivation-objective.md` ruled the triexponential-fit route canonical, because it
is the pipeline the published numbers came out of. Measured by agent 9 at 50 LHS + 2 polish:

| scheme | route | all six | τ only | amplitudes only |
|---|---|---|---|---|
| `kinetic_jw95_free_rates` | eigen | 527.86 | **0.00** | 452.99 |
| `kinetic_jw95_free_rates` | **fit (canonical)** | 411.20 | **123.19** | **71.06** |
| `kinetic_open_desens_free_rates` | eigen | 530.64 | **0.00** | 452.99 |
| `kinetic_open_desens_free_rates` | **fit (canonical)** | 1144.72 | **106.77** | 1031.15 |

The τ-only χ² of 0.00 is a feature of the **eigen route only**. Canonically it is 123.19 and
106.77 — nowhere near a fit.

**3. For the reference scheme the canonical route reverses the finding's direction.** For
`kinetic_jw95_free_rates`, amplitudes-only (**71.06**) fits *better* than τ-only (**123.19**).
The obstruction is not "in the amplitudes" there at all. For `kinetic_open_desens_free_rates`
the amplitudes are indeed worse (1031.15 vs 106.77), so the pattern is **scheme-dependent** —
precisely what a structural claim cannot be.

Agent 9 reported the finding as standing, summarising these as "driven by amplitudes
(amplitudes-only χ² = 71.06 to 1031.15; taus-only χ² = 106.77 to 123.19)". Quoting ranges
across schemes conceals that for the reference scheme the amplitude number is the lower one.
That verdict is rejected, as 7A's was, and for the same reason: a conclusion read off a
favourable presentation of the objective rather than off the objective.

## Correction to "What survives" — it survives less than I first said

Written an hour after the withdrawal, having followed the arithmetic one step further.

**The amplitude floor dominates the headline χ² as well.** Of the all-six eigen-route figure of
527.86 for `kinetic_jw95_free_rates`, **452.99 is the arithmetic floor** and only **74.87** is
attributable to anything else. So "χ² in the hundreds against ~6 for a genuine fit" is mostly a
restatement of the floor, not evidence about the model.

**And the amplitude term is misspecified at a more basic level.** The published fractional
amplitudes are 0.57, 0.23, 0.20 — they **sum to exactly 1.0**. They therefore carry **two**
independent degrees of freedom, not three, and the objective scores all three as independent
measurements with independent errors. Their published errors (0.04, 0.02, 0.03) cannot all be
independent either, for the same reason. Scoring a constrained composition as if unconstrained
double-counts and inflates χ², which is precisely the effect seen.

So the claim that **no scheme reproduces the published three-component deactivation** is no
longer established either. It may well be true — 74.87 on a handful of effective degrees of
freedom is still poor — but the number that was supposed to demonstrate it is contaminated by
the same flaw. Brief 10 re-derives it with the composition handled properly.

What is left standing, and it is the only thing in this episode that has survived every test:
the methodological result below about pseudo-data leverage.

## What survives (as originally written — now qualified by the correction above)

* **No scheme tested reproduces the published three-component deactivation of α1β2γ2**, under
  either route: χ² in the hundreds against ~6 for a genuine fit. That is a real model-adequacy
  result and it is unaffected by the above.
* `kinetic_flip`'s failure to improve on `kinetic_jw95_free_rates` (527.41 against 527.86, two
  extra free parameters, unchanged at 6× budget) stands **as a measurement**. Its
  interpretation is now unclear, since part of the objective it was scored on has an arithmetic
  floor.
* The methodological point stands untouched and is the most durable thing here: scoring a
  101-point reconstruction as data, when it is deterministic given six published numbers,
  manufactured a 10× leverage imbalance and drove the brief 5 / brief 6 dispute.

## What a sound version would need

1. An amplitude objective **not** built on a sum-to-one normalisation, so that a floor set by
   the published numbers alone cannot masquerade as a model limit. Score raw amplitudes against
   `P_o(0)`, or score the waveform's shape directly.
2. A fit-route budget adequate to the dimension. 50 Latin-hypercube samples with 2 polished, in
   a six-parameter box, is not enough to separate a result from a failed search — and the
   canonical route costs ~0.9 s per evaluation even warm-started, so this needs planning rather
   than optimism.
3. Any decomposition claim stated **per scheme**, never as a range across schemes.

---

# (preserved below as originally written)

# Finding — the deactivation failure is an amplitude obstruction, not a rate one

**Status: PROVISIONAL.** Measured by Claude Opus 5, 2026-10-17, verifying agents 8A and 8B.
One re-derivation is required before this is stated as settled — see "What would overturn
this".

## The result

Scored against the published Barberis 2007 deactivation parameters at their published errors
(τ = 2.8 ± 0.3, 33.4 ± 4.6, 221.35 ± 14.9 ms; fractional amplitudes 0.57 ± 0.04, 0.23 ± 0.02,
0.20 ± 0.03; n = 6), no scheme in this project reproduces α1β2γ2 deactivation:

| scheme | free params | best χ² over the six numbers |
|---|---|---|
| `kinetic_jw95_free_rates` (5-state) | 6 | **527.86** |
| `kinetic_flip` (6-state, pre-open) | 8 | **527.41** |
| `extended_desens` (6-state, two sinks) | 4 | 1107.73 |
| `kinetic_open_desens` (5-state) | 3 | 2181.92 |

A genuine fit would be of order 6.

**Decomposing the objective localises the obstruction completely.** Fitting the time constants
and the amplitudes separately, same code path, same budget:

| fitted | `kinetic_jw95_free_rates` | `kinetic_flip` |
|---|---|---|
| all six | χ² = 527.86 | χ² = 527.86 |
| **time constants only** | **χ² = 0.00** | **χ² = 0.00** |
| amplitudes only | χ² = 452.99 | χ² = 452.99 |

**The time constants are not the problem at all.** Both schemes place the three relaxation time
constants at exactly (2.8, 33.4, 221.35) ms. The obstruction is entirely in the **amplitudes**,
and it is not relieved by freeing the time constants: fitting amplitudes alone, with τ free to
go anywhere — the optimiser sent one to 8.5e5 ms — still cannot do better than χ² = 452.99.

**The achievable amplitude triple is ≈ (0, 0.405, 0.595); the measurement requires
(0.57, 0.23, 0.20).** The fastest relaxation mode cannot carry appreciable amplitude in these
schemes, while the published fit gives the fastest component the **largest** amplitude. When τ
is forced onto the published values, amplitudes diverge into large cancelling pairs —
(−0.001, −20.5, +21.5) for JW95 — which is not a decay waveform.

**Adding states does not help, and the invariance is the striking part.** `kinetic_flip`
inserts a shut, fully-liganded pre-open state between `A2R` and `A2O` — the flip mechanism
(Lape, Krashia, Colquhoun & Sivilotti 2008), the standard remedy for exactly this kind of
shortfall — and differs from `kinetic_jw95` in that one state and nothing else. With **two
extra free parameters** it improves χ² by **0.45 in 527.86**, i.e. 0.09%. In the decomposed
runs the two schemes agree to **two decimal places on all three numbers**. Re-run at **6× the
search budget (9000 Latin-hypercube starts, 40 polished) and a different seed**, the figures
were unchanged: 527.86 and 527.41. This is not an under-search.

## Why — WITHDRAWN, I was wrong twice

This section originally argued that open probability starts at its maximum, that a
0.57-amplitude component at τ = 2.8 ms implies an initial slope near 0.2/ms, and that the
liganded shut states buffer the open state and forbid such a decay. **I tested that and it is
false.** Both attempts to explain the obstruction failed:

**Attempt 1 — an initial-slope limit.** The published components fix the initial normalised
decay rate at Σ aᵢ/τᵢ = **0.2114 ± 0.0261 /ms** (96.3% of it from the fast component). The same
quantity for a model is `-(p_end @ Q₀)[open] / p_end[open]`, computable from the generator with
no decomposition, no fit and no selection rule. Maximised over each box:

| scheme | max achievable Σaᵢ/τᵢ | vs required | `po_max` there |
|---|---|---|---|
| `kinetic_jw95_free_rates` | 3.0925 /ms | **14.63×** | 0.0099 |
| `kinetic_flip` | 3.3538 /ms | 15.87× | 0.0000 |
| `kinetic_open_desens` | 1.9392 /ms | 9.17× | 0.9990 |
| `extended_desens` | 0.2917 /ms | 1.38× | 0.0099 |

Every scheme can decay far faster than required, so there is no slope limit. The test also
measured the wrong proxy: the model reaches a large Σaᵢ/τᵢ with a *tiny* amplitude at a *tinier*
τ (0.003 at 0.0051 ms), where the published value comes from a moderate amplitude at a moderate
τ.

**Attempt 2 — a fast mode that cannot carry amplitude.** Maximising the fractional amplitude
carried by modes with τ inside the published τ₁ ± 1σ window (2.5–3.1 ms), with `po_max`
constrained to [0.70, 0.80], `kinetic_jw95_free_rates` reaches **1.2404** at `po_max` = 0.7000,
against the 0.57 required. So that is not the obstruction either.

**A methodological hazard surfaced instead, and it bears on the numbers above.** The same
search returned "amplitudes" of **194** and **9013**. Fractional amplitudes are `aᵢ / Σaⱼ`, and
where `Σaⱼ → 0` that normalisation is ill-conditioned and explodes. An optimiser will find such
regions. **Whether the χ² = 452.99 amplitudes-only figure is partly exploiting the same
ill-conditioning has not been checked**, and brief 9 must check it: report the raw `Σaⱼ` at
every reported optimum and reject any cell where it is near zero.

**So the measurements in the table above stand, and the explanation does not.** The obstruction
is some joint constraint over the components, `po_max` and the waveform that I have not
identified. I am recording the two failed explanations rather than deleting them, because a
plausible mechanism written into a findings document and left unchallenged is exactly the
failure this project exists to prevent — and I wrote it.

## Why this matters

It reframes what briefs 5 and 6 were arguing about. Those briefs asked whether the Barberis
deactivation and paired-pulse recovery measurements are mutually incompatible, and got
opposite answers under different weightings of a 101-point reconstructed curve against a
1-point scalar. **The curve is deterministic given the six published numbers and the grid**, so
scoring it as 101 observations manufactured a 10× leverage imbalance with no principle
selecting its direction (measured: 4764, 11170 and 110.6 units of `n/σ²` against the scalar's
1111). Scored instead against the six real measurements, the apparent dataset conflict is a
symptom: the schemes get the deactivation **shape** wrong, and "incompatibility" is what that
failure looks like when an objective is free to trade one dataset against another.

## What would overturn this

1. **The selection rule.** These numbers come from the eigen route, which yields four non-zero
   modes for a five-state scheme and then selects three by largest |amplitude|. That rule has
   no experimental counterpart, and `DECISION-deactivation-objective.md` has already ruled that
   the canonical path is the triexponential-fit route, which replicates the measurement
   procedure instead. **This finding must be re-derived through the fit route before it is
   stated as settled.** The identical χ² across schemes is suggestive of a real constraint, but
   it is also exactly what a binding selection rule would produce.
2. **Bounded search.** No vector was *found* within the declared bounds at the stated budget.
   That is not a proof that none exists.
3. **Single open state.** Every scheme tested has one. The amplitude argument above predicts
   that this is the binding restriction.
4. **The amplitude convention.** Published "areas" are treated as fractional amplitudes, which
   the repo supports by recomputing Σ A_i τ_i = 53.55 ms against the published τ_w of
   52.5 ± 2.9. If the paper meant fractional *integrated areas* instead, the comparison target
   changes and so does every number here.

Numbers above were produced by `scripts/deact_structural_search.py` and
`scripts/flip_structural_search.py` (committed, on `claude/agent-brief-8a` and
`claude/agent-brief-8b`), and by the rate/amplitude decomposition that brief 9 is tasked with
committing as `scripts/deact_obstruction_decomposition.py`.
