# Agent brief 13 — cooperative binding, and an unsourced target that needs resolving

**Read first:** `FINDING-minimax-conflict.md`, `FINDING-two-open-states.md`,
`DECISION-deactivation-objective.md` **with its amendment**, `AGENT-BRIEF-6.md` §0.
`FINDING-amplitude-obstruction.md` is withdrawn — do not cite it.
`PYTHONPATH=<your worktree>/src`. Interpreter `/Users/whinz/Biochemistry/.venv/bin/python`.

**Your worktree:** branch `claude/agent-brief-13`, based on `main` (the whole brief 6–12 chain is
now merged there). One agent.

Your brief-12 work is accepted, including the first accurate test-count report in this project —
574 fast, verified independently. One framing correction is in `FINDING-minimax-conflict.md`:
the conflict is a **trade-off**, not "structural and topological", because your own sensitivity
table reaches nH = 2.09 (−0.28σ) on the CRC alone. Each target is individually attainable; no
single vector gets all four.

## Task 1 — the mechanism the broken observable points at

The minimax is now **4.34σ** for `kinetic_two_open`, and what breaks is the **Hill slope**:
predicted **0.66** against published **2.2 ± 0.4**. Everything else is met — deactivation 0.50σ,
EC50 −0.02σ, `po_max` +0.10σ.

**A Hill slope of 2.2 on a two-site receptor requires strong positive cooperativity in binding,
and none of the project's schemes has any.** Every one of them uses independent, identical sites
with the statistical factors (`2·kon·g` for the first binding, `2·koff` for the first
unbinding). That is by construction non-cooperative, so a steep CRC has to be manufactured by
gating instead — which is exactly the trade-off that collapses τ₃ to 64.9 s in your CRC-alone
row.

So: **add binding cooperativity as a single isolated structural change**, the same way
`kinetic_flip` and `kinetic_two_open` isolated theirs.

Add `src/circuitpharm/models/kinetic_coop.py`: take `kinetic_two_open` exactly as it is and give
the **second** binding step its own dissociation rate, `koff2`, with the first unchanged. Define
the cooperativity factor `c = koff / koff2` and say in the docstring which direction `c > 1`
means. Change **nothing else** — same states, same topology, same two conducting open states, so
that the only difference from `kinetic_two_open` is the one parameter.

**Microscopic reversibility:** check whether this introduces any cycle requiring a constraint.
`kinetic_two_open` is a branched chain, so adding `koff2` probably does not — **verify that
yourself and state it**, do not take my word.

Connectivity test with explicit zero assertions, as in `tests/test_nextgen_models.py`. Declare
nesting: at `c = 1` (`koff2 = koff`) the scheme must reduce **exactly** to `kinetic_two_open`.
Use `nests` / `nesting_embedding` and **assert the embedding reproduces the sub-model's
likelihood to numerical tolerance** — that test exists and has caught real defects.

Then add a `kinetic_coop` spec and re-run `scripts/joint_fourway_fit.py` unchanged, same code
path and budget, for `kinetic_coop`, `kinetic_two_open` and `kinetic_jw95_free_rates`.

Answer in `RESULTS-BRIEF-13.md`, with every residual in σ, per scheme:

* **What is the minimax residual for `kinetic_coop`?** Compare against 4.34σ and 14.25σ.
* **Does it close the 3σ gap?** Yes or no.
* What `c` does it need, and is that value physically sensible — report `koff` and `koff2`
  explicitly, not just the ratio.
* Report `n_data`, `k_fitted`, residual dof. **`kinetic_coop` has 9 parameters against `n_data`
  = 9, so residual dof is 0 — the fit is exactly saturated and is not a test at all.** Say so
  prominently. If it fits, that is an existence result and nothing more, and the honest
  comparison is against `kinetic_two_open`'s 8 parameters at the same budget.
* Every bound flag, and `k_effective`.

**Do not invent a weighting.** Every quantity has a published error. If you need a relative
weight, stop and report `kind=question`.

## Task 2 — an unsourced number is being used as a project target

`FIT_TARGETS['ec50_um'] = 20.0` µM. The sourced value is Jahn 1997's **11.6 ± 0.9 µM**. Your
own discrepancy register records 20.0 as *"not stated / recalled literature middle from earlier
legacy tuning files"* — i.e. **UNSOURCED**, and it disagrees with the measurement by **+9.33σ**.

This project exists to separate what is measured from what is assumed, and an unsourced number
is currently serving as a fitting target. **Resolve it in favour of the sourced value**, and do
it carefully:

1. Find **every** use of `FIT_TARGETS['ec50_um']` across `src/`, `scripts/` and `tests/`.
2. Change the target to the sourced 11.6 µM with its 0.9 µM SEM, attributed to
   `jahn_a1b2g2_kinetics`, recorded with `Basis.QUANTITATIVE` if that is what `provenance.py`
   actually uses for a sourced measurement — **check the code, do not assume; a previous brief
   of mine stated `Basis.CONVENTION` where the code records `FITTED`.**
3. **Expect this to move other numbers.** Run the full tier and report **every** test whose
   expected value changes, with the old and new figures. Do not update a pinned number without
   listing it.
4. If any test's expectation cannot be reconciled with the sourced value, **report
   `kind=question` rather than editing the test.** A target change that quietly rewrites
   assertions is how a project loses its anchors.

Keep the old value visible in a comment with its history and why it was replaced.

## Task 3 — verify

Both tiers, **final pytest line pasted verbatim**. Baseline on `main` after the brief-12 merge:
measure it yourself and state it. Do not touch or gate the Apple Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] `kinetic_coop.py`: `kinetic_two_open` plus `koff2` only; connectivity test with explicit zeros; microscopic reversibility checked and stated.
- [ ] Nesting declared and the embedding asserted to reproduce `kinetic_two_open`'s likelihood at `c = 1`.
- [ ] Four-way fit re-run for three schemes, same code path and budget, residuals in σ, no invented weighting.
- [ ] `RESULTS-BRIEF-13.md` states the minimax for `kinetic_coop` and whether it closes 3σ, with `koff`/`koff2` values and all bound flags.
- [ ] Residual dof = 0 for `kinetic_coop` stated prominently, with the saturation caveat.
- [ ] `FIT_TARGETS['ec50_um']` resolved to the sourced value with provenance verified against the code; every affected test listed with old and new numbers; old value preserved in a comment.
- [ ] Both tiers pass; pytest lines verbatim.
- [ ] Commit per task; `kind=blocked` with what landed rather than running to the timeout.

## Not yours to decide

`kind=question`: adoption of any scheme into the canonical set; promotion of any finding to
settled; any test whose expectation cannot be reconciled with the sourced EC50.
