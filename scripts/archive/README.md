# Archived scripts

These are **superseded one-off runs**, kept for provenance only. Several of them produced
results that were later overturned, and the trail matters — but none of them should be run
to produce a current number, and several will not even reproduce their original output
because the modules they call have since changed.

`WORKLOG.md` records what each one established and where it was corrected.

## Why each was retired

**Old all-spiking locomotor CPG** — `run_cpg.py`, `sweep_cpg.py`, `tune_cpg.py`
Explored `HalfCentreCPG`, whose rhythm depends on mutual inhibition. An exhaustive
288-point search showed a LIF half-centre cannot reach a physiological locomotor rhythm
(best score 0.83). The reasoning was right about LIF half-centres and wrong about the
conclusion: a *group pacemaker* needs no plateau potential and works fine in LIF. Replaced
by `circuitpharm/rg2.py`. The Matsuoka oscillator that followed was also removed, because
in it inhibition IS the oscillator, which contradicts the lamprey strychnine phenotype and
invalidated a prediction this project had to retract.

**Pre-split pharmacology panels** — `drug_panel.py`, `subtype_run.py`, `subtype_sweep.py`
All pass a single `gaba_sens` to the circuits, so one gain is applied to both the tonic and
phasic receptor pools. The GABA-A Markov scheme shows those gains differ by ~7x, so these
are wrong by about sevenfold on whichever pool they were not calibrated against. Superseded
by `circuitpharm.evaluate`.

**One-off result runs** — `final2.py`, `final_combo.py`, `confirm.py`, `margin.py`,
`optimise.py`, `overdose.py`
Produced the headline numbers of sessions 3-6, including overdose indices and absolute
ventilation percentages. Those quantities are now **VOID**: the efficacy ceiling they clip
against is pool-dependent, and the respiratory axis's only calibration anchor was
invalidated (human whole-body ventilation is the wrong observable for an isolated preBötC).
`overdose.py` additionally assumes dose escalation is monotonic, which benzodiazepine
potentiation is not — it is bell-shaped and antagonistic above ~10 µM.

**Calibration and tuning one-offs** — `calibrate_resp.py` is NOT here (it documents the
now-VOID anchor and is still worth reading); these are:
`ceiling_test.py`, `sweep_gain.py`, `tune_reflex.py`, `tune_resp.py`,
`calib_nmda_resp.py`, `propofol_check.py`
Single-purpose parameter searches whose outcomes are recorded in `WORKLOG.md`.
`calib_nmda_resp.py` is the sweep that reduced the NMDA share of recurrent excitation from
0.55 to 0.06 and still gave 56% ventilation at 60% block — the search that led to declaring
the NMDA respiratory limitation. Note that declaration was initially made on insufficient
grounds (the disinhibition pathway `ei_nmda` was absent from the circuit entirely, so this
sweep could not have found it); see `scripts/calib_ei_nmda.py` and
`scripts/diag_inhib_load.py` for the proper test, which upheld the conclusion for a
stronger reason.

## A recurring lesson these encode

Several of these scripts produced results that *looked fine* — a flat sweep, a plausible
percentage, a reassuring safety margin — and were wrong. That is why the failure modes
E1–E12 are now executable regression tests in `tests/test_failure_modes.py` rather than
prose, and why every quantity the package reports carries a reliability tier that refuses
to be read when it is void.
