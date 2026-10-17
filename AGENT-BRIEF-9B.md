# Agent brief 9B — continue brief 9, and stop polling your own background jobs

**You are continuing `AGENT-BRIEF-9.md`, same worktree, same branch
(`claude/agent-brief-9`).** Read that brief again; everything in it still applies. Read also
the updated `FINDING-amplitude-obstruction.md`, whose "Why" section I have **withdrawn** — see
§3 below.

## 0. Why your last run produced nothing, and what to do differently

Your previous run hit its 4-hour timeout without reporting. Its final output was a sequence of:

> *"I have launched X in the background and will await its completion."*
> *"X is progressing. I will await its completion."*

**You spent the run polling your own background jobs.** That is the single thing to change.

* **Run things in the foreground.** A command that takes ten minutes and blocks is cheaper than
  one you launch and then check on ten times.
* If something genuinely must run detached, **launch it, then do other useful work, and read
  its output file once.** Never poll a job whose output you are not yet ready to use.
* **Commit as you finish each task**, not at the end. Your merge (934ab2d) survived because it
  was committed; tasks 2 and 3 were uncommitted when the run died and I had to preserve them
  myself as `cd7af80`, unverified.
* If you are running out of room, **report `kind=blocked` with what is done** rather than
  continuing until the timeout. A partial result reported is worth more than a complete one
  that is never delivered.

## 1. What is already committed on your branch

* `934ab2d` — task 1, the 8A/8B merge. Done.
* `cd7af80` — **WIP, yours, unverified by me**: `scripts/deact_obstruction_decomposition.py`,
  `src/circuitpharm/models/kinetic_two_open.py`, and edits to `comparison.py`,
  `deact_components.py`, `models/base.py`, `models/__init__.py`, `tests/test_nesting.py`,
  `tests/test_nextgen_models.py`. Pick up from there. No test tier has been run against it.

## 2. The fit route is expensive — budget it explicitly

I smoke-tested your decomposition script at a deliberately tiny budget
(`--n-samples-eigen 150 --n-polish-eigen 3 --n-samples-fit 12 --n-polish-fit 2`) and it had not
finished in two minutes. Each fit-route objective evaluation simulates a waveform and then
multi-start fits a triexponential to it, so it is orders of magnitude slower than the eigen
route.

**Therefore:** do fewer schemes properly rather than four badly. Priority order:

1. `kinetic_jw95_free_rates` — the reference.
2. `kinetic_open_desens_free_rates` — **the new spec from task 3.1, and the most interesting
   cell in the table**, because that topology has never had its rates freed.
3. `kinetic_flip`.
4. `extended_desens`.

Get 1 and 2 done and reported before touching 3 and 4. State the budget used per cell and the
runtime. A cell you did not run is reported as not run — **never as a number carried over from
the eigen route.**

## 3. A hazard that must now be checked, because it may undermine the finding

`FINDING-amplitude-obstruction.md`'s mechanism is withdrawn: I proposed one, tested it, and it
was false twice. Details are in that file. **Do not resurrect either explanation.** The
measurements stand; the "why" is open.

What came out of those tests is a methodological problem you must handle:

**Fractional amplitudes `aᵢ / Σaⱼ` are ill-conditioned where `Σaⱼ → 0`.** My searches returned
"amplitudes" of **194** and **9013** from exactly that. An optimiser will find those regions,
and **it is unchecked whether the χ² = 452.99 amplitudes-only figure is partly exploiting the
same ill-conditioning.**

So, for **every** optimum your decomposition reports:

* Report the **raw** `Σaⱼ` alongside the fractional amplitudes, and the true `P_o` at the start
  of the relaxation. They should agree closely — `Σaⱼ` *is* `P_o(0)` for the eigen route.
* **Reject any cell where `|Σaⱼ|` is small relative to `P_o(0)`** — pick and state a threshold —
  and report it as ill-conditioned rather than as a χ².
* Add a test pinning this: a vector whose raw amplitude sum is near zero must be rejected by
  the objective, not scored.

Then state explicitly in `RESULTS-BRIEF-9.md` whether the 452.99 figure survives this check.
**If it does not, the finding is withdrawn** and that file must say so.

## 4. Acceptance criteria (replacing brief 9's, which still otherwise stand)

- [ ] Each task committed as it completes, not batched.
- [ ] Fit-route vs eigen-route decomposition for at least `kinetic_jw95_free_rates` and `kinetic_open_desens_free_rates`, all three objectives (all six / τ only / amplitudes only), with the budget and runtime stated per cell.
- [ ] Raw `Σaⱼ` and `P_o(0)` reported at every optimum; ill-conditioned cells rejected and labelled, with the threshold stated.
- [ ] A test rejecting near-zero-amplitude-sum vectors.
- [ ] `RESULTS-BRIEF-9.md` states plainly whether the obstruction is still in the amplitudes under the canonical route, and whether χ² = 452.99 survives the conditioning check.
- [ ] `FINDING-amplitude-obstruction.md` updated in place if the finding falls, preserving what it claimed.
- [ ] Degenerate fit-route cells reported as degenerate, never backfilled from the eigen route.
- [ ] Both test tiers run with `PYTHONPATH` set, counts stated, labelled and accounted for against 546 fast / 618 full on `claude/agent-brief-8a`.
- [ ] Task 4 of brief 9 (recovery, Jahn CRC, `po_max`, bound flags) reported in full for any scheme that reaches χ² of order tens.

## 5. Not yours to decide

`kind=question`, do not decide: adoption of any new scheme into the model set; promotion of
`FINDING-amplitude-obstruction.md` from provisional to settled.

If you hit the wall again, **report `kind=blocked` with what landed.** That is the correct
outcome, not a failure.
