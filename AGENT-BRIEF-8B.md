# Agent brief 8B — does a pre-open conformational step reproduce the measured deactivation?

**Read `/Users/whinz/Biochemistry/BRIEF-7-SHARED.md` and `AGENT-BRIEF-6.md` §0 first.** The
discipline rules apply unchanged.

**Your worktree:** already created, branch `claude/agent-brief-8b`, based on
`claude/agent-brief-7a`, so `src/circuitpharm/fitting/deact_components.py` and
`predict_components_eigen_all` are already there and **verified** — an independent
reimplementation agrees to 6.66e-16 and its components reconstruct the true relaxation to
6.72e-15. Use it; do not rewrite it.

**Agent 8A is working in parallel** and owns `src/circuitpharm/fitting/deact_components.py` and
`src/circuitpharm/literature.py`. **Do not modify either.** You own
`src/circuitpharm/models/` and your own script and tests.

Run tests with `PYTHONPATH=<your worktree>/src`. Interpreter
`/Users/whinz/Biochemistry/.venv/bin/python`.

---

## The question

Measured by the planner, using the verified predictor: **no parameter vector was found, in any
of the project's three schemes, that reproduces the published three-component deactivation of
α1β2γ2.** Best χ² over the six published numbers — τ = 2.8 ± 0.3, 33.4 ± 4.6, 221.35 ± 14.9 ms
and fractional amplitudes 0.57 ± 0.04, 0.23 ± 0.02, 0.20 ± 0.03 — was **527.9** for
`kinetic_jw95_free_rates` (6 free rates), 1107.7 for `extended_desens`, 2181.9 for
`kinetic_open_desens`, against ~6 for a genuine fit.

The failure has a consistent signature. **In every scheme the fastest mode's amplitude
collapses towards zero** (σ ≈ −14.3 on area₁), while the data require the **fastest** component
(2.8 ms) to carry the **largest** amplitude (0.57). At the best JW95 vector one single mode
carries **98.2%** of the amplitude, so the predicted deactivation is effectively
mono-exponential where the measurement is clearly triexponential.

**Your question:** does adding a **pre-open conformational step** — a shut, fully-liganded
intermediate between `A2R` and `A2O` — produce a fast, large-amplitude component and so
reproduce the measurement?

This is the "flip" / pre-open mechanism (Lape, Krashia, Colquhoun & Sivilotti 2008 and the
subsequent literature on glycine and nicotinic receptors), and it is the obvious candidate
precisely because it puts a fast conformational equilibrium *between* binding and opening,
where the current schemes have nothing.

**Both answers are publishable and neither is expected.** If it does not work, that is a
sharper negative result than the one we have, because it rules out the standard remedy.

## Your tasks

### Task 1 — the scheme

Add `src/circuitpharm/models/kinetic_flip.py` with a 6-state model:

```
R  <->  AR  <->  A2R  <->  A2F  <->  A2O
                   |
                 A2D
```

`A2F` is shut and fully liganded; `A2O` is open; desensitisation stays on `A2R` as in JW95
(`A2R <-> A2D`), so that **the only structural difference from `kinetic_jw95` is the inserted
pre-open state.** That isolates the variable. Follow `kinetic_jw95.py`'s interface exactly:
`q_matrix(conc_um)`, `po_max`, `native_observable`, `state_distribution`, `simulate_waveform`.
State the state order in the docstring; the open index is what `deact_components` needs.

Rate constants: `kon`, `koff` (binding), `delta`/`gamma` (A2R→A2F / A2F→A2R, the flipping
equilibrium `F = delta/gamma`), `beta`/`alpha` (A2F→A2O / A2O→A2F), `d`/`r` (desensitisation).
**Microscopic reversibility:** there is no cycle in this scheme, so no constraint beyond
positivity — confirm that for yourself and say so, rather than taking my word.

Add a connectivity test in the style of
`tests/test_nextgen_models.py::test_kinetic_jw95_a2d_connectivity_is_strictly_to_a2r`, with
explicit negative assertions: `A2F` connects only to `A2R` and `A2O`; `A2D` connects only to
`A2R`; `A2O` connects only to `A2F`. Assert the zeros, not just the non-zeros — that is how the
earlier topology error was caught.

### Task 2 — a spec, and the search

Add a `MODEL_SPECS` entry `kinetic_flip` in `src/circuitpharm/fitting/comparison.py`.
**Important:** 8A owns `deact_components.py` and `literature.py` but *not* `comparison.py`, so
this is yours — keep the edit to the new spec and its factory, additive, nothing else touched.

Declare its nesting honestly. At `delta → ∞` with `gamma` fixed, or whatever limit actually
collapses `A2F` onto `A2R`, the scheme should reduce to `kinetic_jw95` — **verify whether it
does** using the `nests` / `nesting_embedding` machinery already in `comparison.py`, and if the
limit sits on a bound or does not hold exactly, **record the measured discrepancy and exclude
it**, as brief 6 did for `extended_desens` (which turned out not to nest JW95 at all, by 16.66
nats, because `r_fast` = 0.005 against `r` = 0.002).

Then write `scripts/flip_structural_search.py`, modelled on the search the planner ran:
Latin-hypercube over the spec's bounds, best 12 polished by Nelder–Mead, minimising χ² over all
six published numbers read from `circuitpharm.literature.DEACTIVATIONS[1]` — **no literals** —
using `predict_components_eigen_all` and selecting the three modes with the largest
|fractional amplitude|.

Report in `RESULTS-BRIEF-8B.md`: the command, runtime, search budget, best χ², the component
table with per-quantity σ, how many modes carry |amplitude| ≥ 0.10, whether three amplitudes
are positive, and the fitted rate constants with any that sit on a bound flagged. **A parameter
on a bound is the edge of the search box, not an estimate** — `at_bound` exists for this.

For comparison in the same table, re-run the same search on `kinetic_jw95_free_rates` and
report its χ² alongside; a cross-scheme claim needs both numbers from the same code path and
the same budget.

### Task 3 — if it works, test it honestly

If `kinetic_flip` reproduces the six numbers at χ² of order tens rather than hundreds, then
**check it against the rest of the evidence before anyone gets excited**:

1. Paired-pulse recovery at 100 ms — `BARBERIS2007_PAIRED_PULSE`, published 0.33 ± 0.03. Report
   the σ miss at the deactivation optimum.
2. The Jahn 1997 peak CRC — `JAHN1997_PEAK_CRC`. Report RMSE against its baseline of 0.004917.
   Brief 6's reconciling vectors missed this by 97.6× and 124.4×; this is where previous
   "reconciliations" died.
3. `po_max` against `FIT_RANGES` [0.70, 0.80]. **0.8004 is outside the band**, not "aligned
   with" it — an earlier report made exactly that slip.

A scheme that fits deactivation and then misses recovery or the CRC has not solved the problem,
and the results file must say so plainly.

### Task 4 — verify

Both tiers with `PYTHONPATH` set, counts stated separately and labelled. Baselines: 524 fast /
596 full on `main`, ~535 fast on the 7A branch you are based on. Do not touch or gate the known
Apple Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] `kinetic_flip.py` with the stated topology; connectivity test with explicit zero assertions.
- [ ] `kinetic_flip` spec in `MODEL_SPECS`, additive; nesting against `kinetic_jw95` either verified or excluded with the measured discrepancy recorded.
- [ ] `scripts/flip_structural_search.py` committed; no published value or error hardcoded.
- [ ] `RESULTS-BRIEF-8B.md` reports χ² for `kinetic_flip` **and** `kinetic_jw95_free_rates` from the same code path and budget, with per-quantity σ, mode counts, and bound flags.
- [ ] If deactivation is fitted, Task 3's three cross-checks are reported with no favourable subset selected — report all of them whichever way they fall.
- [ ] Both tiers pass, counts labelled.

Report `kind=question` rather than deciding: whether `kinetic_flip` is adopted into the
project's model set. Building and measuring it is yours; adopting it is not.
