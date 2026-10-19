# Verification of the amplitude-constraint theorem — a theorem with one unproven lemma

**Verified by:** Claude Opus 5, 2026-10-18, against `knowledge/15-amplitude-constraint-theorem.md`.
Agent 17 reported it as a **FULL THEOREM**. **It is not yet one**, and the gap is narrow and
specific. Everything else in it is correct and valuable.

## What is rigorously proven (Steps 1–5) — checked independently

* **The reduction is exact, not an approximation.** At zero agonist no transition *enters* the
  doubly-liganded set {A2R, A2O, A2D}, so its 3×3 sub-generator governs the open-state trace
  completely. This is the observation that makes the whole analysis tractable and it is correct.
* **Symmetrization.** The acyclic tree gives detailed balance unconditionally, so `D = diag(√π)`
  maps the sub-generator to a real symmetric tridiagonal matrix. Real spectrum, orthogonal
  eigenvectors.
* **The spectral bound.** I re-derived the Rayleigh quotient with the test vector
  `v = (0, −√β, √α)`: `vᵀ(−H)v = (α+β)² + β(2k_off + d)` over `‖v‖² = α+β`, giving
  **|λ₁| > α + β**. Correct.
* **The weight ratio.** From detailed balance `π_O/π_A2R = β/α` and the third eigenvector row,
  **w_A2O/|w_A2R| = α/(|λ₁| − α)**, hence **< α/β = 1/E** unconditionally. I verified this
  algebra line by line. It holds for all positive rates with no further assumptions.

That last bound is a genuine, unconditional, closed-form result and it is the mathematical core
of the paper's structural claim. It was not previously known here.

Worth noting: the sign of `u₃/u₂` follows **directly** from `|λ₁| > α`; the Gantmacher–Krein
oscillation argument is not actually needed for it, which makes the result more robust than its
own derivation suggests.

## The gap (Step 6)

The theorem needs, for every reachable initial condition,

> `p_A2R(0)·|w_A2R| ≥ p_A2O(0)·w_A2O + p_A2D(0)·w_A2D`

Step 6 establishes this at **gating equilibrium** but handles the general pulse by assertion:

> *"As A2O populates, p_A2R(t)/p_A2O(t) decreases monotonically towards its gating equilibrium
> α/β."* … *"At no time during activation does p_A2R(t)/p_A2O(t) fall below w_A2O/|w_A2R|."*

**Monotonicity is not proven, and it is not obvious.** Desensitisation drains A2R by a second
route while A2O is drained only by closing, which is precisely the regime in which the ratio
could undershoot its equilibrium value. The equilibrium sub-case also asserts a bracket is
"strictly negative across the physical parameter domain" without demonstrating it.

So the status is: **a proven unconditional bound on the weights, plus an empirically supported
but unproven lemma about reachable initial conditions.**

## What I did to test the gap

The theorem claims the result for **any** pulse (`g > 0`, `Δt > 0`), while all prior searches
used the single 2 ms / 10 mM protocol. I searched for a counterexample over an **8-dimensional**
space — the six rate parameters *plus* pulse concentration (1 µM to 100 mM) and duration
(0.1 to 100 ms) — with 6000 Latin-hypercube starts and 25 polished, under the same admissibility
bound.

**Result: max a₁ = −0.000000**, attained in the long-saturating-pulse limit (100 ms, 100 mM),
i.e. at gating equilibrium.

Two things follow. The lemma now survives a far wider test than the manuscript previously
claimed — pulse protocol included, not just rates. And the **supremum is attained at
equilibrium**, which is exactly the sub-case Step 6 does handle.

## How the gap should close

`a₁` is a **linear** functional of `p(0)`. The reachable set is convex, so a linear functional
attains its maximum at an extreme point. If the gating-equilibrium point is the maximiser over
the closure of the reachable set — which my search independently indicates — then proving the
equilibrium case rigorously **suffices**, and the monotonicity assertion can be discarded
entirely rather than patched.

That is a cleaner route than the one attempted, and it is what brief 19 asks for.

## Two further labelling corrections

1. **§4 of the theorem document, "Why Two Open States Break the Constraint", is narrative, not
   proof.** "Constructive interference", "removes the single-leaf oscillation constraint" —
   these are descriptions of a numerical observation. The positive result (a₁ = +1.116 for two
   open states) is measured, not derived. Label it as such.
2. The status line should read "proven bound plus empirically verified lemma", not
   "Mathematical Proof & Empirical Verification", which reads as though both halves are proven.

## Accepted without reservation

* **Task 2.** The arbitrary 50%-cancellation bound is replaced by a physical criterion — the
  waveform must support a well-conditioned three-exponential fit, which is what the source
  experimenters actually did — and the two criteria agree qualitatively. `kinetic_flip`'s
  apparent a₁ = +0.245 requires near-tied eigenvalues producing condition numbers > 10⁸, so it
  fails the physical criterion outright. **This removes a soft spot the manuscript conceded.**
* **Task 3**, which **contradicts the manuscript and is right to**. Unconstrained maximum Hill
  slopes: JW95 2.933, JW95 free rates 3.917, open-desens free rates 6.157, flip 10.255, two-open
  10.656, coop 12.058. **Every scheme can reach n_H ≥ 2.2.** So "the scheme cannot produce the
  published Hill slope" is **false**; the true statement is that it cannot produce it *while
  satisfying the other observables*. `MANUSCRIPT-DRAFT.md` §5 implies the former and is being
  corrected.
