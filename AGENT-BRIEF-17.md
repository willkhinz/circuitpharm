# Agent brief 17 — prove the amplitude constraint instead of searching for it

**Read first:** `MANUSCRIPT-DRAFT.md` §2 and §8, `STATE-OF-RESULTS.md`,
`FINDING-two-open-states.md`, `AGENT-BRIEF-6.md` §0. `PYTHONPATH=<your worktree>/src`.

**Your worktree:** branch `claude/agent-brief-17`, based on `main` (brief 16 is merged). One agent.

Brief 16 is accepted in full. The bidirectional profile, the plateau finding that dissolved the
c = 1000 liability, the risetime prediction, and `STATE-OF-RESULTS.md` are all good work.

The manuscript is drafted. **The standard has now been raised to airtight**, and this brief
closes the largest gap: the paper's central structural claim is currently a bounded numerical
search.

## 1. The claim to prove

Measured, over the declared parameter boxes, with at most 50% inter-mode cancellation:

> In five-state schemes with a **single** conducting open state, the fastest relaxation mode's
> fractional amplitude in the open-state projection satisfies **a₁ ≤ 0**, with zero cancellation
> at the optimum. The measured value is +0.57.

This is reported as "no admissible vector was found", which is not a proof of non-existence. **A
proof would make it the strongest result in the paper.** Attempt it.

## 2. The route that should work

The pieces are standard and they fit this problem well:

1. **The scheme is detailed-balanced.** At zero agonist the five-state generator is an acyclic
   tree (`R – AR – A2R – A2O`, with `A2D` branching off `A2R`), so Kolmogorov's criterion holds
   unconditionally — you verified exactly this reasoning for `kinetic_coop` in brief 13.
2. **Therefore Q is symmetrizable.** With π the stationary distribution and
   `S = diag(√π)`, the matrix `S Q S⁻¹` is **symmetric**. Hence the spectrum is real and
   non-positive and the eigenvectors are orthogonal under the π-weighted inner product. You
   already observe real eigenvalues numerically; this is why.
3. **Sign structure.** For a symmetric tridiagonal (Jacobi) matrix, Sturm oscillation theory /
   Gantmacher–Krein gives that the eigenvector of the k-th eigenvalue has exactly k−1 sign
   changes along the chain; the extremal (fastest) mode alternates. The discrete nodal domain
   theorem extends this to trees. **The open state is a leaf at the end of the path**, which is
   what makes its component's sign determined rather than free.
4. **The amplitude is a product.** `a₁ = c₁ · v₁[O]`, where `c₁` is the projection of the
   end-of-pulse occupancy onto the corresponding eigendirection. So the proof needs the sign of
   **both** factors, and the second depends on the initial condition — which, after a saturating
   2 ms pulse of 10 mM GABA, is concentrated in the doubly-liganded states. That restriction on
   reachable initial conditions is part of the theorem, not an aside.

**Do not assume any of the above is correct because I wrote it.** Check each step; if the nodal
domain argument does not transfer to this tree, say so and say why.

## 3. What to deliver, in descending order of value

1. **A full proof**, stated as a theorem with explicit hypotheses (topology, positivity of rates,
   the class of reachable initial conditions) and a complete argument. Put it in
   `knowledge/`, in a file named for what it proves, and **mirror the key steps as executable
   assertions** in `tests/` — e.g. that the symmetrized generator is symmetric to machine
   precision, that eigenvalues are real and non-positive, that the sign pattern is as claimed —
   so the proof's premises are checked against the code on every run.
2. **A partial proof**, if that is as far as it goes: prove whatever sub-case is reachable (for
   instance the fully saturating-pulse limit, or the chain without the desensitised branch) and
   state precisely what remains conjectural.
3. **A counterexample.** If a₁ > 0 is in fact attainable with an admissible decomposition, that
   overturns §2 of the manuscript. **Report it immediately as `kind=done` with the vector** —
   finding one would be a better outcome than a proof, because it would mean the paper is wrong
   in a way we can fix before anyone else finds it.
4. Only if 1–3 all fail: a far larger search with explicit convergence evidence, reported
   honestly as a search.

**Be clear which of these you delivered.** A proof sketch presented as a proof is worse than an
honest search.

## 4. Task 2 — replace my arbitrary admissibility bound

Manuscript §8 concedes that the `Σ|aᵢ|/|Σaᵢ| ≤ 1.5` criterion is a convention I chose, and that
the two six-state results press against it and depend on it. That is not airtight.

**Replace it with a physical criterion.** The natural one is already in the repo: a
decomposition is admissible if the model's waveform **supports a well-conditioned
three-exponential fit** — because that is precisely what the source experimenters did, and a
waveform that cannot support such a fit could not have produced the published numbers. Brief 10
already pinned the degenerate case as a test.

Re-derive the §2 table under that criterion, report it alongside the L1 version, and state
whether the qualitative conclusion — 5-state single-open forced ≤ 0, two-open exceeding 0.57 —
holds under both. **If the two criteria disagree, say so prominently**; the manuscript then needs
rewriting and I would rather know now.

## 5. Task 3 — the Hill slope ceiling

Manuscript §5 claims the maximum achievable Hill slope among admissible fits is 1.25, against a
published 2.2 ± 0.4. That too is a bounded search.

Maximise n_H directly over each scheme's parameter box, with no other observable constrained,
and report the unconstrained maximum per scheme. This separates two very different statements:
"the scheme cannot produce n_H = 2.2 at all" from "it cannot produce it while satisfying the
other observables". **The manuscript currently implies the first and has only measured the
second.** Say which is true.

## 6. Verify

Both tiers, **final pytest line pasted verbatim**. Baseline on `main`: brief 16 measured 588
fast / 662 full. Do not touch or gate the Apple Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] A theorem, a partial theorem, or a counterexample for the a₁ constraint, with its status stated unambiguously.
- [ ] Proof premises mirrored as executable assertions in `tests/`.
- [ ] §2's table re-derived under the fit-conditioning criterion; agreement or disagreement with the L1 version stated prominently.
- [ ] Unconstrained maximum n_H per scheme, with the two Hill-slope statements distinguished.
- [ ] Both tiers pass; pytest lines verbatim.
- [ ] Commit per task; `kind=blocked` with what landed rather than running to the timeout.

## Not yours to decide

`kind=question`: any rewrite of the manuscript's claims — report what you measured or proved and
I will revise the text; adoption of any scheme; promotion of any finding to settled.
