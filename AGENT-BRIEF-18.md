# Agent brief 18 — three open states, because the source says three

**Read first:** `FINDING-keramidas-fulltext.md` (new, and the reason for this brief),
`MANUSCRIPT-DRAFT.md` §3 and §5, `FINDING-coop-five-way.md`, `AGENT-BRIEF-6.md` §0.
`PYTHONPATH=<your worktree>/src`.

**Your worktree:** branch `claude/agent-brief-18`, based on `main`. One agent.

**Agent 17 is working in parallel** on the analytic proof of the a₁ constraint. It owns
`knowledge/`, the proof tests, and the admissibility criterion inside `deact_components.py`.
**You own `src/circuitpharm/models/`, your own script, and an additive spec entry.** Keep your
edit to `comparison.py` to the new spec and its factory and nothing else, so the merge is clean.

## 1. Why this brief exists

I read Keramidas & Harrison 2008 in full text (PMC2213567). The paper our mean-open-time
observable comes from reports, for **α1β2γ2S** — the same subunit combination as our
deactivation source — that:

> "GABA-gated bursts required three shut and **three open components** to adequately describe
> M-Mode activation (… = 10 patches). This is the same number of components reported for
> αβγGABARs."

and that the high-P_o mode has "two shut and **three open** states".

Our best scheme, `kinetic_coop`, has **two** conducting open states. An independent
single-channel experiment on the same receptor requires three, and says other studies agree.

This is simultaneously the first external corroboration of the manuscript's structural claim —
derived here from macroscopic amplitudes, corroborated there by dwell-time analysis — and a
strong hint that two is not enough. **The manuscript's one unresolved tension is the Hill slope
ceiling of 1.25 against a published 2.2 ± 0.4, and adding open states raised the achievable
slope in every test so far.**

## 2. Task 1 — build it

Add `src/circuitpharm/models/kinetic_three_open.py`: take `kinetic_coop` exactly as it is and
add a **third** conducting open state, keeping everything else — the cooperative binding step,
the desensitisation branch, the connectivity of the existing states — unchanged, so that the
only structural difference is the extra open state.

**You must decide and justify the connectivity of the third open state**, because the source
does not dictate it. Candidates: in series after the second open state, or in parallel from the
same shut state. **State which you chose, why, and what the alternative would predict
differently.** If both are defensible, implement both and report both — this is exactly the kind
of choice that should not be made silently.

Requirements as before: connectivity test with explicit **zero** assertions; microscopic
reversibility checked and stated (a cycle here would need a constraint, unlike the acyclic trees
so far — **check, do not assume**); the observable is the summed occupancy of all three open
states; `open_indices` exposed so the component predictors receive it without a hard-coded index.

Declare nesting: the scheme must reduce exactly to `kinetic_coop` in the limit that collapses
the third open state, verified through `nesting_embedding` with the likelihood asserted to
match. Warm-start from `kinetic_coop`'s optimum and **refuse to report a worse value** — that
assertion has caught real failures twice in this project.

## 3. Task 2 — the question

Run the five-way fit (`scripts/joint_fiveway_fit.py`), same code path and budget, for
`kinetic_three_open` alongside `kinetic_coop` and `kinetic_two_open`. Report per scheme, every
residual in σ, plus `n_data`, `k_fitted`, **residual dof**.

Then answer plainly in `RESULTS-BRIEF-18.md`:

1. **What is the maximum achievable Hill slope for the three-open scheme?** Measure it both
   unconstrained over the box and at the five-way optimum, as brief 17 is doing for the others.
   Published target: **2.2 ± 0.4**. Current ceiling among admissible two-open fits: **1.25**.
2. Does the five-way minimax improve on `kinetic_coop`'s 2.58σ?
3. **Watch the degrees of freedom.** A third open state adds parameters, and `kinetic_coop` is
   already at residual dof +1 with ten data. If the third open state pushes residual dof to ≤ 0,
   **the fit stops being a test** and you must say so prominently — a better minimax at zero
   dof is not an improvement, it is saturation. Report dof before reporting the minimax.
4. Does the cooperativity requirement change? `kinetic_coop` needs c ∈ [35, 2500]; if three open
   states reduce it toward published ratios of 5–23×, that is a substantive result.

## 4. Task 3 — do not overclaim the corroboration

`FINDING-keramidas-fulltext.md` §2 records a caveat you must respect and carry into your report:
Keramidas's three open **components** are dwell-time components fitted **within bursts** at
saturating GABA. Our open **states** are states in a macroscopic gating scheme. **These are not
automatically the same objects.** State the mapping you are assuming between a burst-level
dwell-time decomposition and a whole-receptor scheme, and flag it as an assumption.

Note also that the two gating modes differ in their **shut** state count — three for M-mode, two
for H-mode. That is modal gating, and **nothing in this project represents it**. Do not attempt
to model it here; just record that it is unrepresented, because a reader will ask.

Update the `circuitpharm.literature` record for Keramidas to note that the full text has now
been read and that the source reports three open components — the record currently carries only
the parameter values and `where` should reflect the full-text read.

## 5. Verify

Both tiers, **final pytest line pasted verbatim**. Baseline on `main`: 588 fast / 662 full. Do
not touch or gate the Apple Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] `kinetic_three_open.py` differing from `kinetic_coop` only in the added open state; connectivity choice justified, alternative stated.
- [ ] Connectivity test with explicit zero assertions; microscopic reversibility checked, not assumed.
- [ ] Nesting to `kinetic_coop` declared, embedding verified, warm-start armed with refusal to report a worse nested optimum.
- [ ] Five-way fit for three schemes, same code path and budget, all residuals in σ.
- [ ] **Residual dof reported before the minimax**, with explicit statement if it reaches ≤ 0.
- [ ] Maximum achievable Hill slope reported, unconstrained and at the optimum.
- [ ] Cooperativity requirement compared against `kinetic_coop`'s [35, 2500].
- [ ] The dwell-time-component vs scheme-state mapping stated as an assumption; modal gating recorded as unrepresented.
- [ ] Keramidas literature record updated to full-text read.
- [ ] Both tiers pass; pytest lines verbatim.
- [ ] Commit per task; `kind=blocked` with what landed rather than running to the timeout.

## Not yours to decide

`kind=question`: adoption of any scheme; promotion of any finding to settled; any rewrite of the
manuscript's claims — report what you measured and I will revise the text.
