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
