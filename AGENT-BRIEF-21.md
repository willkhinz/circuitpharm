# Agent brief 21 — build the digitisation pipeline; do NOT invent a single data point

**Read first:** `STATE-OF-RESULTS.md`, `MANUSCRIPT-DRAFT.md` §8, `AGENT-BRIEF-6.md` §0, and
`src/circuitpharm/fitting/data.py` around line 1106 (the `promote_by` text that defines this
task). `PYTHONPATH=<your worktree>/src`.

**Your worktree:** branch `claude/agent-brief-21`, based on `main`. One agent.

Brief 20 is accepted in full, including both withdrawals — the risetime prediction and the
over-determined-test claim. Reporting those cost you nothing and improved the paper.

## 0. THE RULE FOR THIS BRIEF, ABOVE ALL OTHERS

The repository's own gate currently says, of the two principal datasets:

> *"IS PROVISIONAL — its conclusions are not usable… the following datasets are generated from
> formulae, not measured… promote by: digitise the concentration-response and deactivation
> curves from named published figures, register them in provenance.py, and set
> synthetic=False"*

`JAHN1997_PEAK_CRC` and `BARBERIS2007_DEACT` are both `kind="parametric"`, `synthetic=True`.
Only `BARBERIS2007_PAIRED_PULSE` is digitised.

**This is the task where fabrication is easiest and most damaging.** A plausible-looking set of
digitised points is trivial to generate and nearly impossible to detect later. So:

**If you cannot obtain the actual published figure as an image file on disk, you must NOT produce
any digitised data point. Report `kind=blocked` and build the machinery in §2 instead.** That is
the correct, expected, and fully acceptable outcome of this brief.

**You may not:** estimate points "from the published parameters", reconstruct a curve and call
it digitised, infer points from the text, or produce any number not read off an image you
actually have. Writing `synthetic=False` on anything not read from a real figure would
invalidate every result in this project at once.

I have checked the retrieval route myself so you do not waste the run: the PubMed Central
full-text tool returns **text only** for PMC1950087 — no figure URLs, no figure captions, zero
`Fig`/`Figure` references. The text is available; the images are not, through that route.

## 1. Task 1 — attempt retrieval, and record exactly what happened

Barberis et al. 2007, *Eur J Neurosci* 25(9):2726–40, **PMC1950087**, DOI
10.1111/j.1460-9568.2007.05530.x. Jahn et al. 1997, *NeuroReport* 8(16):3443–6 (no PMC record
known; likely not open access).

Attempt, in order, and **record the outcome of each attempt verbatim** in
`RESULTS-BRIEF-21.md`, including failures:

1. Whether any locally available tool returns figure images for PMC1950087.
2. Whether the article's figures exist as files anywhere already on this machine (search the
   repo and any `data/` or download directories — **do not** reach outside the project for
   arbitrary files).

**Do not scrape, and do not circumvent any access control or paywall.** If the images are not
available through a legitimate route you already have, that is the answer: report it.

## 2. Task 2 — build the pipeline so the data can be dropped in later

This is the deliverable that does not depend on obtaining anything, and it is the real value of
this brief. Build it so that a digitised dataset **cannot be registered without its evidence**.

1. **A digitised-dataset loader** in `src/circuitpharm/fitting/` that constructs a dataset of
   `kind="digitised"`, `synthetic=False` **only** when all of the following are present, and
   raises a clear error otherwise:
   * a source image file that exists on disk, with its **SHA-256 recorded** in the dataset
     record;
   * the figure identifier (journal, DOI, figure number, panel);
   * the axis calibration actually used (the two reference points per axis, in both pixel and
     data coordinates), so the transformation is reproducible and auditable;
   * the digitisation method and operator recorded as a string.
   **Make the missing-image case a hard failure, not a warning.** The point is that a future
   agent or human physically cannot set `synthetic=False` without the artefact.
2. **A self-validation check.** A correctly digitised Barberis deactivation trace, re-fitted
   with a triexponential by the project's own `predict_components_fit` machinery, must reproduce
   the published τ = (2.8, 33.4, 221.35) ms and areas (0.57, 0.23, 0.20) **within their
   published errors** (0.3, 4.6, 14.9 and 0.04, 0.02, 0.03). Implement that as the acceptance
   test for any future digitisation. A fabricated or badly-calibrated curve fails it.
   Likewise a digitised Jahn CRC must reproduce EC50 = 11.6 ± 0.9 µM and n_H = 2.2 ± 0.4.
   **Read all published values from the literature records; no literals.**
3. **Tests** covering: the loader refuses without an image; refuses on checksum mismatch;
   accepts a correctly-formed record; and the self-validation catches a deliberately perturbed
   curve (perturb a synthetic one in the test itself — that is a fixture, not data, and label it
   as such).

## 3. Task 3 — quantify what the gate actually blocks

Separately and usefully: establish **which live claims depend on the synthetic datasets**. Not
all of them do.

For each of the project's current findings — the pseudo-data leverage result, the amplitude
constraint theorem, the necessity theorem, the Hill-slope bracket, and `kinetic_coop`'s
five-observable consistency — state whether it depends on the **synthetic curves**, on the
**published parameters** (which are real measurements), or on **neither** (the theorems depend
on neither; they are statements about model structure).

This matters because weighting (iv) scores against published parameters rather than the
generated curves, so some findings may be unaffected by the synthetic flag entirely.
**Report it per claim**, and say plainly which conclusions would change if real digitised curves
replaced the reconstructions and which would not.

## 4. Verify

Both tiers, **final pytest line pasted verbatim**; `pytest-xdist` is active and has given
inconsistent totals, so **run twice and report both**. Baseline on `main`: I measure 631 fast /
705 full after brief 20 — confirm it. Do not touch or gate the Apple Silicon split-R̂ artefact.

## Acceptance criteria

- [ ] Retrieval attempts recorded verbatim, including every failure. **No digitised data point produced unless a real image was obtained.**
- [ ] Loader built; registration **impossible** without an on-disk image plus recorded SHA-256, figure identity, axis calibration and method.
- [ ] Self-validation implemented against the published components and CRC parameters, read from the literature records with no literals.
- [ ] Tests covering refusal-without-image, checksum mismatch, acceptance, and detection of a perturbed curve.
- [ ] Per-claim dependency analysis: synthetic curves / published parameters / neither.
- [ ] Both tiers pass; pytest lines verbatim, twice.
- [ ] Commit per task.

## Not yours to decide

`kind=question` or `kind=blocked`: anything requiring the figures themselves. **A blocked report
here is a success, not a failure** — it tells the user exactly what to supply.
