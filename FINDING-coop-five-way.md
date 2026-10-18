# Finding — `kinetic_coop` satisfies five observables, on a convention that decides the result

**Status: PROVISIONAL.** Recorded by Claude Opus 5, 2026-10-17, verifying agent 15 against its
own code and arithmetic. This is the project's first fit with **positive residual degrees of
freedom**, which makes it the first one that is a test at all. It also turns on a scoring
convention, and that is stated before the result rather than after.

## The result

Five observables, each scored against its own published values and errors — deactivation
components, CRC (EC50 and nH), paired-pulse recovery, `po_max` band, and now **mean open time**:

| scheme | k | n_data | **residual dof** | five-way minimax | satisfies all five? |
|---|---|---|---|---|---|
| `kinetic_jw95_free_rates` | 6 | 10 | +4 | **14.25σ** | no |
| `kinetic_two_open` | 8 | 10 | +2 | **3.23σ** | no |
| `kinetic_coop` | 9 | 10 | **+1** | **2.58σ** | **yes** |

The nesting assertion held on both the joint and minimax objectives
(`coop ≤ two_open + 1e-4`).

## The convention that decides it, stated plainly

The deciding observable is **mean open time**, and the pass depends on how an unstated error
kind is handled.

Keramidas & Harrison 2008 Table III, α1β2γ2S — the **same subunit combination as Barberis** —
gives M-mode components τ = (0.49 ± 0.02, 2.58 ± 0.09, 5.17 ± 0.38) ms with areas
(0.26, 0.49, 0.25) and `err_kind = UNSTATED`. The area-weighted mean open time is
**2.684 ms** (verified independently: 0.26·0.49 + 0.49·2.58 + 0.25·5.17 = 2.684), and
propagating the published component errors gives the band **[2.54, 2.83] ms**. **The band is
derived from published errors, not invented** — I checked, because an invented σ here would have
manufactured the whole result.

`kinetic_coop` predicts **τ_O = 3.20 ms**, which is **outside** that band:

| convention | residual | verdict |
|---|---|---|
| **band** — zero penalty inside, excess beyond the edge over the half-width | (3.20 − 2.83)/0.145 = **2.55σ** | **passes** |
| **point** — deviation from the midpoint over the half-width | (3.20 − 2.684)/0.145 = **3.56σ** | **fails** |

Both are in agent 15's own table, which is the honest way to present it. But it means the
headline is convention-dependent: **the model misses the measured mean open time by 0.52 ms, a
19% relative error, and passes only because the source does not state whether its errors are SD
or SEM.** The band treatment is defensible precisely because of that silence — committing to a
Gaussian σ the paper never claimed would be worse — but the result should be quoted as
"2.58σ under a band treatment of an unstated error", never as "2.58σ".

## One genuine corroboration

`kinetic_coop` needs `koff` = 85.67 ms⁻¹ and `koff2` = 0.0857 ms⁻¹. Brief 15 was forbidden from
gating `koff` against the interval `koff_identifiability.py` propagates from the Jahn CRC
(0.0848–0.1421 ms⁻¹) because that propagation assumes the JW95 topology — see
`AGENT-BRIEF-15.md` §2, where I withdrew the gate before it was implemented.

Agent 15 then noticed something I had missed: **`koff2` = 0.0857 ms⁻¹ falls inside that window.**
Under cooperativity the quantity the CRC actually constrains is the *doubly-bound* unbinding
rate, which is `koff2`, not `koff`. So the JW95-derived interval and the cooperative fit agree
on the rate they have in common. That is a real consistency check and it was not engineered —
the interval was computed under a different topology entirely.

Caveat: 0.0857 sits at the **bottom edge** of 0.0848–0.1421, not centrally, so this is
corroboration and not confirmation.

## What is still wrong with it

1. **Cooperativity has grown, not shrunk.** `c = koff/koff2` is now **1000** (was 815.94 at four
   observables), against published site-to-site affinity ratios of roughly 23× (Olsen et al.
   1981) and 5× (Agey & Dunn 1989) — so **43–200× beyond anything published**, while positive
   cooperativity itself is confirmed (Mozrzymas et al. 2003). All three are abstracts, now
   registered with `where="abstract"` and explicitly barred from use as fitting constraints.
2. **Two of the ten data are bands** (`po_max` and mean open time), scored with half-width as σ.
   That convention is pre-existing in the project for `po_max`, but a band carries less
   information than a measurement with a stated error, so `res_dof = +1` overstates how
   over-determined the fit really is.
3. **A per-component amplitude residual exceeds 3σ.** The five-way row has `z_a2 = +3.08σ`
   while the deactivation family aggregates to √χ²_a = 2.58σ. That is not a contradiction —
   the amplitude χ² is the rank-2 quadratic form on the 2-simplex from brief 10, since three
   fractional amplitudes summing to one carry two degrees of freedom — but a naive reading of
   the per-component column would look like a violation, and it should be labelled as such
   wherever it is quoted.
4. **Li et al. 2008 would be more favourable and was correctly not used preferentially**:
   γ2L, sub-saturating 50 µM, weighted mean 2.96 ms — closer to coop's 3.20 ms. Choosing the
   source that fits better would be exactly the wrong move; Keramidas is the right source
   because it matches the preparation.

## Status

The first result in this project with positive residual degrees of freedom and a scheme that
satisfies every observable. **Not settled**, because it rests on a band convention that, if
replaced by the standard point treatment, fails at 3.56σ; and because it requires cooperativity
1.5–2 orders of magnitude beyond published estimates. The honest headline is that **a scheme
with two conducting open states and strong binding cooperativity is the first to come within
reach of all five observables at once, and the cooperativity it demands is its main liability.**

---

## Addendum, 2026-10-17 — the cooperativity liability evaporates; the Hill slope is the real tension

Agent 16 profiled `c` across 23 log-spaced points, 1–10,000, with four-pass bidirectional
warm-starting (down from 1000, up from 1000, up from 1, down from 10,000) to remove optimiser
noise and hysteresis. 18.4 min, ~48 s per point. I read the profile JSON directly; the numbers
below are from it, not from the report.

### c = 1000 was an artefact of an under-constrained plateau

| c | minimax | | c | minimax |
|---|---|---|---|---|
| 1.0 | 3.232σ | | 100.0 | **2.751σ** |
| 5.0 | 3.078σ | | 158.5 | **2.709σ** |
| **23.0** | **3.013σ** | | 1000.0 | **2.585σ** |
| **35.0** | **2.909σ** | | **1584.9** | **2.373σ** |
| 50.0 | **2.841σ** | | 2511.9 | **2.868σ** |
| 75.0 | **2.745σ** | | 3981.1 | 4.314σ |

An asymmetric trough: a broad **plateau over c ∈ [35, 2500]** where all five observables pass,
a shallow global minimum at c = 1585 (2.373σ), and a steep wall above 2500 bound by mean open
time and deactivation. **So c = 1000 was never required** — it is simply where a minimax
optimiser settled inside a flat basin.

**The liability largely evaporates.** The lowest passing cooperativity is **c = 35**, and at
Olsen's published ratio of **23×** the minimax is **3.013σ** — missing the 3σ criterion by
**0.013σ**, which is inside any honest optimiser tolerance. The model needs roughly **1.5×** the
published cooperativity, not 43–200×.

### But the Hill slope never reaches the published value, anywhere

Reading the Hill slope along the same profile:

| c | nH | z_nH | mean open time | P(AR) peak |
|---|---|---|---|---|
| 1.0 | 0.91 | −3.23σ | 2.22 ms | 8.81% |
| 23.0 | 0.99 | −3.01σ | 3.27 ms | 2.62% |
| 35.0 | 1.04 | −2.91σ | 2.83 ms | 2.34% |
| 100.0 | 1.10 | −2.75σ | 2.52 ms | 1.61% |
| 1000.0 | 1.18 | −2.55σ | 3.20 ms | 0.55% |
| 10000.0 | 1.63 | −1.43σ | 1.37 ms | 0.39% |

* **Maximum Hill slope anywhere on the profile: 1.63** — and only at c = 10,000, where the
  minimax is 8.10σ and mean open time collapses to 1.37 ms.
* **Maximum Hill slope among passing points: 1.25.** Published: **2.2 ± 0.4**.
* At **every** passing point the Hill slope residual is between **−2.55σ and −2.91σ**.

**So "satisfies all five observables" means, precisely: its worst residual is just under 3σ, and
that worst residual is always the Hill slope.** Cooperativity steepens nH monotonically, exactly
as brief 13 predicted, but cannot carry it past ~1.25 without destroying mean open time and
deactivation. **The concentration–response steepness is the irreducible tension** — not
deactivation, not recovery, not `po_max`.

That is a sharper and more useful conclusion than the pass/fail. It is also a prediction: any
scheme in this class should fall ~2.5σ short on Hill slope.

### A falsifiable prediction that the fitted data cannot settle

Because the plateau is flat, the five fitted observables **cannot determine `c`**. Agent 16
identified something that can: the **10–90% current risetime** under rapid application of 10 µM
GABA on outside-out macropatches — **3.54 ms at c = 23 against 9.19 ms at c = 1000**, a 2.6×
difference. The mechanism is explicit: at high `c` the first unbinding is ultra-fast
(`koff` = 85.7 ms⁻¹, 11.7 µs dwell, Kd₁ = 401 µM), so many abortive unbindings occur before a
second agonist binds, and the rise is slow.

This is an experimental-design result: one measurement discriminates across a basin that the
entire fitted dataset leaves flat. **No literature value was sought or assigned for it** — it
stands as a prediction.

### Status

`kinetic_coop` passes five observables with residual dof +1 at any c ∈ [35, 2500], needing only
~1.5× the published cooperativity. Two caveats remain, and the second is now the headline:
the mean-open-time pass still depends on the band treatment of an unstated error (see above),
and **the Hill slope is 2.5–2.9σ short at every passing point and cannot be fixed within this
model class.**
