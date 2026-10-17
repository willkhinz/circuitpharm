# Agent brief 15 — a fifth observable, so the coop result becomes testable

**Read first:** `FINDING-two-open-states.md` **including the 2026-10-17 addendum**,
`RESULTS-BRIEF-14.md`, `DECISION-deactivation-objective.md` with its amendment,
`AGENT-BRIEF-6.md` §0. `FINDING-amplitude-obstruction.md` and
`FINDING-minimax-conflict.md` are withdrawn — do not cite their claims.
`PYTHONPATH=<your worktree>/src`.

**Your worktree:** branch `claude/agent-brief-15`, based on `claude/agent-brief-14`. One agent.

Brief 14 is accepted in full — the four-way objective, the nesting assertion, the manuscript
update and the derived falsification threshold. It produced the project's first positive result:
`kinetic_coop` satisfies all four observables at ≤2.39σ with every parameter interior.

**It is not yet a test, and that is the whole of this brief.**

## 1. The problem: residual dof = 0

`kinetic_coop` fits **9 free parameters** to **9 independent data**. Exactly saturated. An exact
or near-exact fit is then guaranteed by parameter count and carries no evidential weight — it is
an existence statement. The remedy is **another observable**, not fewer parameters.

## 2. A correction before you start: the koff interval does NOT transfer

I had intended this brief to gate `koff` against the interval `scripts/koff_identifiability.py`
propagates from the Jahn CRC. **I regenerated it and then talked myself out of it; do not
implement it.**

The 24 replicates give `koff` ∈ **0.0848–0.1421 ms⁻¹**, against `kinetic_coop`'s fitted
**75.55 ms⁻¹** — about 530× above the top. That looks decisive and is not, because **the
propagation is model-dependent**: it derives `koff` from the CRC *under the JW95 structure*.
`kinetic_coop` has cooperative binding and two open states, so the CRC-to-`koff` relationship is
different by construction — and the proof is that `kinetic_coop` **already satisfies the CRC**
at `z_EC50` = +2.39σ *while* holding `koff` = 75.55. A constraint derived under one structure
cannot be imported to exclude a parameter in another.

**Record this reasoning in `RESULTS-BRIEF-15.md`** with both numbers. It is the correct handling
of a cross-model inference, and this project has been bitten repeatedly by the opposite.

## 3. Task 1 — add mean open time, honestly

`circuitpharm.literature.OPEN_TIME_FITS` holds three records:

| label | receptor | errors |
|---|---|---|
| Keramidas 2008 M-Mode | α1β2γ2S | **UNSTATED** |
| Keramidas 2008 H-Mode | α1β2γ2S | **UNSTATED** |
| Li 2008 / Li 2007b control | α1β2γ2**L** | SD |

Keramidas is on **α1β2γ2S — the same subunit combination as Barberis**, which is why it is the
right choice. But:

1. **Its errors are UNSTATED, and you may not invent one.** A Gaussian term needs a σ this source
   does not provide. Score it as an **interval/band constraint** instead — the project already
   does this for `po_max` via `FIT_RANGES` — or state clearly why some other honest treatment is
   better. **Do not assign a plausible-looking SD.**
2. **Do not assume two reported modes means two open states.** Keramidas reports an M-mode and an
   H-mode, which in the single-channel literature usually means **modal gating** — the channel
   switching between kinetic modes over seconds — not two open states within a burst. These are
   different claims and conflating them would manufacture exactly the agreement we are testing
   for. Read what the record actually says (`where` field), state your interpretation
   explicitly, and if the record is insufficient to decide, say so and treat the two modes as
   **two separate band constraints** rather than one bi-exponential.
3. Li 2008 is γ2**L**, a different splice variant. The project has substituted γ2S for γ2L
   before and **recorded it**; if you use Li, record the substitution the same way.

Implement the predictor: for a scheme with one open state the mean open time is 1/α; for
`kinetic_two_open` and `kinetic_coop` with two conducting states it is not, and you must derive
it properly from the generator — the mean dwell in the open set before leaving it. State the
derivation in the docstring.

## 4. Task 2 — re-run with five observables

Same code path and budget as brief 14's `obj_fourway_minimax`, extended to five. Report per
scheme — `kinetic_jw95_free_rates`, `kinetic_two_open`, `kinetic_coop` — with every residual:

* The five-way minimax, and **`n_data`, `k_fitted`, residual dof**. For `kinetic_coop` this
  should now be **+1 or better**, which is the point of the brief.
* **Does `kinetic_coop` still satisfy everything?** If it does, it is for the first time a
  *tested* fit rather than a saturated one, and that is a real result. If adding mean open time
  breaks it, say so plainly — that is equally a result, and a sharper one.
* Keep the nesting assertion armed: `coop` must not score worse than `two_open` on the shared
  objective.

## 5. Task 3 — register the cooperativity literature

The addendum to `FINDING-two-open-states.md` records that the fit needs a site-to-site affinity
ratio of **815.94**, against published ratios of roughly **23×** (Olsen et al. 1981,
*Mol Pharmacol* — K_D 13 ± 6 nM and 300 ± 150 nM) and **5×** (Agey et al. 1989, *Biochemistry* —
monoliganded forms at 16 and 82 nM), with positive cooperativity itself confirmed by Mozrzymas
et al. 2003, *J Neurophysiol*.

**I read abstracts only.** Register these three in the `sources` table and `provenance.py` with
that recorded honestly — `where="abstract"`, and the basis reflecting what an abstract supports.
**Check what `provenance.py` actually uses; do not assume the Basis enum member.** Resolve each
DOI/PMID and record it. If a DOI does not resolve, **say so and do not register it** — a
resolving identifier proves existence, nothing more.

Do **not** use these numbers as a fitting constraint. They are a plausibility comparison on a
different preparation, and Olsen reads its two affinity classes as distinct receptor populations
rather than intramolecular cooperativity.

## 6. Verify

Both tiers, **final pytest line pasted verbatim**. Baselines measured by me on `main`:
574 fast / 646 full; brief 14 reported 577 / 651 on its branch. Do not touch or gate the Apple
Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] Mean-open-time predictor derived properly for multi-open-state schemes, derivation in the docstring.
- [ ] Keramidas scored **without an invented σ**; treatment justified; modal-gating vs two-open-states distinction stated explicitly rather than assumed.
- [ ] Any γ2L→γ2S substitution recorded as the project has done before.
- [ ] Five-way fit per scheme with `n_data`, `k_fitted`, residual dof; `kinetic_coop` at residual dof ≥ +1.
- [ ] A plain statement of whether `kinetic_coop` survives as a *tested* fit.
- [ ] Nesting assertion armed and its outcome reported.
- [ ] The koff non-transfer reasoning recorded with both numbers (0.0848–0.1421 vs 75.55 ms⁻¹).
- [ ] Three sources registered with `where="abstract"`, identifiers resolved, basis checked against the code.
- [ ] Both tiers pass; pytest lines verbatim.
- [ ] Commit per task; `kind=blocked` with what landed rather than running to the timeout.

## Not yours to decide

`kind=question`: adoption of any scheme into the canonical set; promotion of any finding to
settled; any treatment of Keramidas that would require inventing an error bar.
