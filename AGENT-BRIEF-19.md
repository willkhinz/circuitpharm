# Agent brief 19 — close the last gap in the theorem, by convexity rather than monotonicity

**Read first:** `VERIFICATION-theorem-17.md` (my review of your proof — read it before the
theorem itself), `knowledge/15-amplitude-constraint-theorem.md`, `MANUSCRIPT-DRAFT.md` §2.2.
`PYTHONPATH=<your worktree>/src`.

**Your worktree:** branch `claude/agent-brief-19`, based on `main` (brief 17 is merged).
**Agent 18 is working in parallel** on `src/circuitpharm/models/` and a three-open-state scheme.
**You own `knowledge/` and the theorem's tests.** Do not touch `models/`.

## 1. What your proof achieved, and what it did not

Steps 1–5 are **correct**, and I checked the algebra independently rather than taking them on
trust: the zero-agonist no-inflow property makes the three-state reduction **exact**; the
Rayleigh quotient with `(0, −√β, √α)` does give `|λ₁| > α + β`; and the closed form

  **w_A2O / |w_A2R| = α / (|λ₁| − α) < 1/E, unconditionally**

is right. That bound is the mathematical core of the paper's structural claim and it is new.
Note also that the sign of `u₃/u₂` follows **directly** from `|λ₁| > α` — the Gantmacher–Krein
argument is not needed for it, which makes the result more robust than your derivation implies.

**Step 6 is not proven**, and you reported the whole as a FULL THEOREM. It asserts:

> *"p_A2R(t)/p_A2O(t) decreases monotonically towards its gating equilibrium"* and *"At no time
> during activation does [it] fall below w_A2O/|w_A2R|."*

Monotonicity is neither obvious nor demonstrated. Desensitisation drains A2R by a **second**
route while A2O is drained only by closing — precisely the regime where the ratio could
undershoot. The equilibrium sub-case also asserts a bracket is "strictly negative across the
physical parameter domain" without showing it.

This is the distinction the previous brief asked you to make explicitly, and it was not made. It
matters more than any single result: **a proof sketch labelled as a proof is worse than an
honest search**, because it stops anyone from looking again.

## 2. The route that should close it — and it is cleaner than monotonicity

Do **not** try to patch the monotonicity claim. Replace it.

`a₁` is a **linear** functional of the initial condition `p(0)`:

  a₁(p) = w_A2O·p_A2O + w_A2R·p_A2R + w_A2D·p_A2D,  with the weights fixed by the rates alone.

A linear functional on a convex set attains its maximum at an extreme point. So:

1. **Characterise the reachable set** `C_reach` = { p(Δt) : p(0⁻) = δ_R, agonist g > 0 for
   duration Δt > 0 }, or a convex superset of it that is still tight enough to be useful. A
   convex superset is perfectly acceptable and is likely easier — if the bound holds on a
   superset it holds on `C_reach`.
2. **Prove a₁ ≤ 0 at the extreme points** of that set. Independent evidence says the binding one
   is gating equilibrium: I searched 8 dimensions — your six rates plus pulse concentration
   (1 µM–100 mM) and duration (0.1–100 ms), 6000 Latin-hypercube starts, 25 polished — and got
   **max a₁ = −0.000000, attained in the long saturating limit**. So the supremum sits exactly
   at the sub-case you already handle.
3. **Prove the equilibrium case properly**, including the A2D term, rather than asserting the
   bracket's sign. With `p_A2D^eq = (d/r)·p_A2R^eq` and the closed forms from Step 5 this should
   reduce to an inequality in α, β, d, r, k_off and |λ₁| that you can establish using
   `|λ₁| > α + β` — or that you will find fails in some corner, which would be just as valuable.

If that programme works the monotonicity assertion is **discarded, not repaired**, and the
theorem closes.

## 3. If it does not close

Then say so precisely. Deliver:

* the strongest sub-case you can prove outright (e.g. `d = 0`, or `E` above some threshold);
* the exact inequality that remains open, written explicitly in the rate constants, so a reader
  can attack it;
* the numerical evidence for it, labelled as evidence.

**A theorem with one clearly-stated open lemma is a perfectly respectable result** and is what
the manuscript currently claims. Do not inflate it again.

## 4. Task 2 — two labelling fixes in the theorem document

1. **§4, "Why Two Open States Break the Constraint", is narrative, not proof.** "Constructive
   interference" and "removes the single-leaf oscillation constraint" describe a numerical
   observation. The positive result (a₁ = +1.116) is measured. Label the section as a
   **conjecture with numerical support**, or prove it — if you can show that two conducting
   states in the open manifold admit a positive extremal projection, that is a second theorem
   and worth having.
2. The status line reads "Mathematical Proof & Empirical Verification", which reads as though
   both halves are proven. Make it state exactly which parts are which.

## 5. Verify

Both tiers, **final pytest line pasted verbatim**. Baseline on `main` after merging brief 17:
603 fast / 677 full. Keep the premise-mirroring tests you added; extend them to any new lemma.
Do not touch or gate the Apple Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] Either the theorem closed by the convexity route with the monotonicity assertion removed, or the open inequality stated explicitly in the rate constants.
- [ ] The equilibrium case proven **including** the A2D term, not asserted.
- [ ] §4 relabelled as conjecture-with-numerical-support, or proven.
- [ ] Status line stating precisely which parts are proven and which are empirical.
- [ ] New lemmas mirrored as executable tests.
- [ ] Both tiers pass; pytest lines verbatim.
- [ ] Commit per task; `kind=blocked` with what landed rather than running to the timeout.

## Not yours to decide

`kind=question`: any change to the manuscript's claims — report what you proved and I will
revise the text; adoption of any scheme; promotion of any finding to settled.
