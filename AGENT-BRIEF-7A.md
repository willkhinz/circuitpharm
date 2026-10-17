# Agent brief 7A — predicted deactivation components from the Q matrix eigenvalues

**Read `/Users/whinz/Biochemistry/BRIEF-7-SHARED.md` first, in full.** It holds the discipline
rules, the published measurement table, why weighting (iv) exists, and the shared acceptance
criteria. This file only adds what is specific to your method.

**Your worktree:** you are already in it, on branch `claude/agent-brief-7a`.
**Another agent (7B) is solving the same problem by a different method, in a different
worktree.** Do not look for its work, do not coordinate, do not try to match it.

---

## Your method: analytic, from the generator matrix

For a linear Markov scheme the relaxation after a step change in agonist is a sum of
exponentials whose **time constants are −1/λ_i for the non-zero eigenvalues λ_i of the Q
matrix** at the post-step agonist concentration, and whose **amplitudes follow from projecting
the initial state occupancy onto the eigenvectors** and then onto the observable.

For the Barberis protocol that means:

1. **Initial condition.** Occupancy at the end of a **2 ms pulse of 10 mM GABA**, starting
   from equilibrium at zero agonist. Integrate or exponentiate the bound-state generator for
   2 ms. State the method you used. The pulse is brief and the concentration is saturating, so
   do not assume the receptor reaches steady state during it — check, and report whether it
   does.
2. **Relaxation.** At zero agonist, decompose the Q matrix: `Q = V diag(λ) V⁻¹`. The zero
   eigenvalue is the stationary direction. The remaining four give τ_i = −1/λ_i.
3. **Observable.** The measured current is proportional to open probability, i.e. occupancy of
   `A2O` (index 3 in `[R, AR, A2R, A2O, A2D]`). Project the decomposition onto that
   coordinate to get the four amplitudes, then normalise to peak as the dataset's
   `normalized_current` is normalised. Confirm how `BARBERIS2007_DEACT` normalises before you
   choose; do not guess.
4. **The 4-vs-3 problem.** You will have four components; the paper reports three. Options
   include reporting the three largest by area, merging components whose τ differ by less than
   their published errors, or comparing only derived quantities such as
   τ_w = Σ aᵢτᵢ / Σ aᵢ. **Pick one, justify it, and state what the choice discards.** Report
   all four components in your results file regardless of which three you score.

Relevant existing code: `src/circuitpharm/models/kinetic_jw95.py` has `q_matrix`;
`src/circuitpharm/models/kinetic_open_desens.py` is the alternative topology;
`scripts/deact_sampling.py` already measures how well components are recovered from a trace
and is worth reading for the conventions it uses.

**Numerical care, not optional.** `D` can be large, which makes Q badly conditioned: measured
in brief 6, an enlarged fit pushed `log10_D` to its bound of 3.0. Use `scipy.linalg.eig` and
**check your decomposition**: assert the eigenvalues are real and non-positive to tolerance
(this scheme is a reversible chain, so complex pairs indicate a conditioning problem or a bug,
not oscillation), assert `Σ amplitudes` reproduces the peak open probability, and reconstruct
the waveform from your components and compare it against `simulate_waveform` on the dataset's
own grid. **Report the max absolute deviation of that reconstruction.** If it exceeds 1e-6,
say so and investigate rather than proceeding.

## Specific deliverables

1. `src/circuitpharm/fitting/deact_components.py` (new) with
   `predict_components_eigen(model, *, pulse_ms=2.0, conc_mm=10.0) -> tuple[taus, areas]`.
   Pure function, no global state, documented with the three steps above.
2. A likelihood term — put it where it belongs, your call, and say where — scoring predicted
   components against `circuitpharm.literature.DEACTIVATIONS[1]`'s published values and
   errors. Read the errors from that record and from its `note` field, which carries the
   component and area errors; **do not retype them as literals.** If they are only available
   in the `note` as prose, parse them there or promote them to structured fields — and if you
   promote them, that is a change to a literature record, so state it prominently in your
   report and keep the prose intact.
3. The self-consistency test from the shared criteria, plus the reconstruction check above as
   its own test.
4. `RESULTS-BRIEF-7A.md`, per the shared criteria, additionally reporting: your 4 components
   at each of the three models' brief-6 optima, the max reconstruction deviation, whether the
   receptor reaches steady state during the 2 ms pulse, and your 4-vs-3 choice with what it
   discards.

## What is yours to decide and what is not

Yours: the 4-vs-3 handling, where the likelihood term lives, the eigen-decomposition details.

**Not yours:** whether weighting (iv) becomes the project default. Report `kind=question`
rather than adopting it.
