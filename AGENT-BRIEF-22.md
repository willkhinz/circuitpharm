# Agent brief 22 — adversarial audit: try to disprove every claim in the paper

**Read first:** `MANUSCRIPT-DRAFT.md` in full, `STATE-OF-RESULTS.md`,
`knowledge/15-amplitude-constraint-theorem.md`, `knowledge/16-positive-direction-theorem.md`,
`AGENT-BRIEF-6.md` §0. `PYTHONPATH=<your worktree>/src`.

**Your worktree:** branch `claude/agent-brief-22`, based on `main`. One agent.

Brief 21 is accepted in full, and reporting `blocked` with zero fabricated points was exactly
right. Your Task 3 finding also corrected a mistake of mine — I had wrongly treated the
provisional-data gate as applying to the live results — and I verified it independently before
accepting it.

## 0. The standing instruction has changed

The requirement is now: **nothing goes in the paper that can be disproven.** Your job this round
is **not** to add results. It is to **attack** the ones we have.

**Breaking a claim is the best possible outcome of this brief.** If you falsify something, you
have saved the paper. Report it immediately and prominently. Do not soften it, do not look for a
reading under which it survives, and do not balance it with reassurance. I would rather lose
three claims now than defend one bad one later.

Treat every number and every sentence in `MANUSCRIPT-DRAFT.md` as a target.

## 1. Task 1 — attack the two theorems at their corners

Both theorems assume strictly positive rates and rely on Jacobi-matrix structure. Attack the
assumptions, not the algebra (the algebra has been checked independently).

Specifically probe:

1. **Degenerate eigenvalues.** The Sturm/Gantmacher–Krein alternation argument needs simple
   eigenvalues. For a Jacobi matrix with non-zero off-diagonals eigenvalues are simple — but the
   off-diagonals are √(dr) and √(αβ). **Search for parameter sets where an off-diagonal is
   numerically negligible** (d → 0, r → 0, α → 0, β → 0, or at the declared bounds) and check
   whether the conclusion still holds or merely becomes unverifiable.
2. **Boundary behaviour.** Evaluate the constraint as each rate approaches its bound and as
   ratios diverge. Does `a₁ ≤ 0` survive? Is the supremum ever attained at exactly 0 and could
   numerical error put it above?
3. **The necessity theorem's dimensional step.** It argues 3 modes ⇒ at most 2 positive. Check:
   can a five-state single-open scheme produce three positive *fitted* triexponential amplitudes
   even though it has only two positive *eigenmode* amplitudes? The published amplitudes come
   from a triexponential **fit**, not from an eigendecomposition. **If a three-exponential fit to
   a two-positive-mode waveform can return three positive amplitudes, the necessity theorem does
   not apply to the published data and must be withdrawn.** This is the single most dangerous
   gap in the paper; test it directly and numerically.
4. Whether either theorem's conclusion depends on the open state being a **leaf**, and what
   happens for topologies where it is not.

## 2. Task 2 — attack the fifth observable

The five-observable result depends on mapping Keramidas's **burst-level dwell-time components**
onto **scheme states**. `FINDING-keramidas-fulltext.md` records this as an assumption. Attack it.

* Re-run the five-way fit scoring mean open time against (a) the **M-mode** band [2.54, 2.83]
  (current choice), (b) the **H-mode** band [5.86, 8.63], and (c) the **full modal span**
  [2.54, 8.63].
* **Does `kinetic_coop` still pass under all three?** If it passes only under the M-mode band —
  the narrowest and the one that happens to be the current choice — then the result is
  convention-dependent in a second, undisclosed way, and the manuscript must say so.
* Modal gating means the receptor switches between modes over seconds. A macroscopic
  deactivation measurement averages over whatever modes were present. **State whether scoring
  against a single mode is defensible at all**, and if not, which band the manuscript should
  use.

## 3. Task 3 — attack the Hill-slope bracket and its falsification boundary

The paper claims any two-site scheme is confined to n_H ∈ (4 − 2√2, 2.000), and offers
n_H > 2.000 as a class-level falsification boundary. The published value is 2.2 ± 0.4, whose
central value is **above the ceiling**.

* **Check the upper bound is really 2.000** and not an artefact of the Hill-slope *definition*
  used. The project extracts n_H over the measured 1–10 µM window, while the analytic limit was
  derived at EC50. Those are different quantities. **Does the 2.000 ceiling hold for the
  windowed definition the paper actually reports?** If the windowed slope can exceed 2.000, the
  falsification boundary is wrong and must be restated.
* Does the bracket depend on the all-or-none readout assumption (only doubly-liganded states
  conduct)? **Test a scheme where the mono-liganded state also conducts** — if that lifts the
  ceiling above 2.000, the boundary is a property of the assumption, not of two-site binding.
* State the bracket's exact hypotheses, and whether α1β2γ2 is known to satisfy them.

## 4. Task 4 — attack everything else

Go through `MANUSCRIPT-DRAFT.md` claim by claim. For **each** numbered or bolded claim, record in
`RESULTS-BRIEF-22.md`:

| claim | what would falsify it | attempted? | survived? |

Pay particular attention to:
* the **leverage figures** (4764 / 11170 / 110.6 vs 1111) — are they definition-dependent?
* the **203.0625 floor** — it is an arithmetic consequence of a sum-to-one normalisation; is it
  still meaningful under the corrected 2-simplex objective, or is it a residue of the withdrawn
  one?
* the claim that **open states cannot change the normalised CRC** — does it hold if a conducting
  state is reachable from a *shut* state other than A2R?
* every σ figure quoted in §3 and §5.

## 5. Verify

Both tiers, **final pytest line pasted verbatim, run twice** (xdist gives inconsistent totals).
Baseline on `main`: 647 fast / 721 full after brief 21 — confirm. Do not touch or gate the Apple
Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] Task 1's four probes attempted, with the triexponential-fit-vs-eigenmode test (1.3) answered **numerically and explicitly** — it can withdraw a theorem.
- [ ] Five-way fit re-run under all three mean-open-time bands; pass/fail stated for each.
- [ ] Hill-slope ceiling checked for the **windowed** definition the paper reports, and under a conducting mono-liganded state.
- [ ] A falsification table covering every claim in the manuscript.
- [ ] Every claim that survives only under a convention **named as such**.
- [ ] Both tiers pass; pytest lines verbatim, twice.
- [ ] Commit per task.

## Not yours to decide

`kind=question`: any rewrite of the manuscript — report what broke and I will cut or rescope it.
**Do not repair a broken claim by weakening the test.** Report the break.
