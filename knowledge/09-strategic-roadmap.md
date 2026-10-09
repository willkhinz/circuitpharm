# Strategic Roadmap & Implementation Directive
### From a numerical phenomenon to rigorous computational neuropharmacology

**Document status:** PRIMARY ENTRY DIRECTIVE. Revision 2 (2026-10-09).
**Audience:** every agent working in this repository, human or model.
**Supersedes:** Revision 1 (2026-10-08), whose six advancements are preserved but
**re-ordered, re-scoped and made executable**. §1 records exactly what changed and why.

---

## 0. How to use this document

Read in this order, then stop and work:

1. §1 — what changed in this revision (5 min; tells you why the order below is not the
   order you may remember).
2. §2 — the seven invariants. **These are not style preferences. A change that violates
   one is reverted, not reviewed.**
3. §3 — the P0 defect register. **No new feature work starts until P0 is closed.**
4. §4 — the phase you have been assigned. Phases are dependency-ordered. Do not start a
   phase whose predecessor is not marked CLOSED in §8.
5. §5 — forbidden patterns. Every one is drawn from a real defect in this repository.
6. §7 — the workflow protocol: branch, test, commit, report.

**If you are an implementation agent and you have been given only a task number, you must
still read §2, §5 and §7 in full before writing code.** They are short. They are the
difference between code that is merged and code that is reverted.

### 0.1 Environment facts you need before the first command

```bash
# this project requires Python 3.11 or 3.12 (pyproject: >=3.11,<3.13)
python3.11 -m venv .venv && . .venv/bin/activate
pip install -e '.[dev]'          # core + pytest; this is what CI's lean job uses
pip install -e '.[all]'          # adds mujoco (body plant) + rdkit (chirality)

pytest -q                        # whole suite, parallel by default (-n auto in pyproject)
pytest -q -n 0 tests/test_x.py   # single file, serial, stdout visible — use when debugging
pytest -m "not slow" -q          # fast subset; the circuit-integration tests are excluded
```

The full suite takes **~15 minutes** on a 4-core machine and is dominated by circuit
integration. Expect that. Do not reduce `-n` to "speed it up"; do not add `-x` to a final
verification run (you need the whole failure list, not the first one).

### 0.2 The two systems that define this codebase

You must understand both before you add a number or a result anywhere.

**`src/circuitpharm/results.py` — reliability tiers.** Every model *output* carries a
`Tier`: `VALIDATED` (reproduces an observable it was not fitted to, or is an exact
identity, or survives full parameter propagation), `UNCALIBRATED` (mechanism sound, no
valid quantitative anchor — ordering yes, scale no), `VOID` (rests on something known to
be invalid; reading `.value` **raises**). A tier is a claim about *evidence*, never about
precision. A number can be exact to ten decimals and still be VOID.

**`src/circuitpharm/provenance.py` — parameter bases.** Every load-bearing *input* carries
a `Basis`: `QUANTITATIVE` (read as a number from a named source), `FROM_QUALITATIVE` (a
source says "predominant", we turned it into a number — the conversion is ours),
`FITTED`, `CONVENTION`, `UNSOURCED`, `GUESS`. `tests/test_provenance.py` fails if a
registered parameter has no record.

These two systems are the package's entire scientific identity. The reason they exist is
written in `results.py`'s docstring and is worth quoting, because it is the failure this
roadmap exists to prevent at a larger scale:

> a stale ventilation percentage was quoted in project summaries for several sessions
> while sitting directly beneath an "UNREACHABLE" flag that the model itself had printed.

**The new `models/`, `protocols/`, `fitting/` and `dynamic_range.py` code added on
2026-10-08 bypasses both systems completely.** 3,177 lines of new science code, and not
one `Quantity`, not one `Basis` record. Closing that gap is phase P2 and it is not
optional.

---

## 1. What changed in this revision, and why

Revision 1 was a good research programme in the wrong order, and it under-specified
exactly the things this repository has historically got wrong. The implementation sprint
that followed it (commits `1328902`, `e2f2d16`) is the evidence: it built all three
models, the decomposition engine, the OED layer and the MCMC sampler **on top of
parameter sets that reproduce none of the project's calibration anchors**, and shipped a
test suite that checks the code runs rather than that the numbers are right.

That is not a criticism of the implementation agent. It did what Revision 1 said. The
document was the defect. Eleven changes follow.

**C1. The fit moves to the front, and becomes the only source of parameters.**
Revision 1's pyramid put "Empirical Data Fitting" at level [1] and
"Mechanistic Decomposition" at [2], but listed Advancement 1 (decomposition) *before*
Advancement 2 (fitting) and marked both "Highest Priority". With no stated dependency, the
sprint built downstream first. Measured consequence, from the shipped defaults:

| | EC₅₀ (equilibrium) | P_o,max | project anchor |
|---|---|---|---|
| `KineticAllosteryModel()` | 6.34 µM | 0.8382 | 20 µM / 0.750 |
| `ExtendedDesensitizationModel()` | 5.18 µM | 0.8571 | 20 µM / 0.750 |
| `OperationalScalarModel()` | 25.0 µM (declared) | 0.75 | 20 µM / 0.750 |

`evaluate_dynamic_range()` and `compare_phasic_tonic_mechanisms()` both default to
`KineticAllosteryModel()`. So the three-tiered headroom hierarchy — the whole point of
Advancement 5 — is currently computed at an EC₅₀ three times too low, which puts ambient
0.4 µM only 16× below EC₅₀ instead of ~50×. Tonic headroom is *precisely* a function of
that distance. The headline quantity is being produced at parameters nothing anchored.

Worse, those defaults are a **chimera of two different fits**. `KineticAllosteryModel`'s
`beta=0.559023, alpha=0.107891` are character-for-character the manuscript's current
(config-pulse) fit; its `kon=0.011244, koff=0.333069` are approximately the *superseded*
1000 µM / 0.30 ms fit. The resulting K_d is 29.62 µM, which is neither the manuscript's
31.95 nor the old draft's 29.69 — a parameter set that exists in no fit, no commit and no
document. That is why its EC₅₀ and P_o,max match nothing.

**Fix: P4 produces a parameter registry; after P4, instantiating a model with literal
defaults outside a unit test is a violation (§2.5).**

**C2. An explicit observable contract.** Revision 1 never said which *observable* each
quantity is defined on. This repository's signature failure is exactly that: a preBötzinger
circuit calibrated against whole-body minute ventilation (now `Tier.VOID`, see
`config.CALIBRATIONS["prebotc_gaba_sens"]`), and a `subjective_index` that cannot be
compared with `regional_sens` although both sound like "forebrain drug engagement".

The new code repeated it twice:
* `fitting/data.py` describes `DOSE_RESPONSE_BENCHMARK` as "macroscopic **peak** current"
  (EC₅₀ 20–30 µM), and `fitting/identifiability.py` fits it with
  `model.dose_response()`, which is the **desensitisation-inclusive equilibrium** curve.
  Different observable. Measured misfit: **χ² = 9290 for 9 data points**, and
  `tests/test_nextgen_models.py::test_fitting_objective_and_identifiability` asserts only
  `cost >= 0.0` and `np.isfinite(cost)`, so it passes.
* `protocols/oed.py::evaluate_model_fit` compares Model A's Hill *peak* response against
  Models B/C's *equilibrium* open probability on one axis and calls the difference
  evidence for a mechanism.

**Fix: P1 defines an `Observable` enum and every model method, every dataset and every
comparison declares which one it speaks. A comparison across two observables raises.**

**C3. Identifiability becomes a question before it is a procedure.** Revision 1 said
"evaluate profile likelihood paths for all 6 microscopic parameters". That guarantees a
useless answer, because the equilibrium dose-response depends on the 6 rates **only
through 3 combinations** — K_d = k_off/k_on, E = β/α, and D = d/r. Verified against the
shipped objective:

```
cost at nominal              = 9290.379443
cost with ALL rates x2       = 9290.379443
cost with ALL rates x10      = 9290.379443
cost with ALL rates x100     = 9290.379443
```

Three of six directions are *exactly* flat — not "practically non-identifiable", but
structurally invisible to that objective. Any profile-likelihood or FIM result computed
there is a property of the objective, not of the data. **Fix: P3 requires a structural
identifiability proof first, performs the analysis in the identifiable
reparameterisation, and only admits absolute rates once a dataset with a kinetic
timescale (deactivation, paired-pulse, single-channel) is in the likelihood.**

**C4. Datasets must be real and registered.** Revision 1 said "digitized patch-clamp
data" with no sourcing requirement. `fitting/data.py` now contains:

```python
_I_DEACT = 0.70 * np.exp(-_T_DEACT / 15.0) + 0.30 * np.exp(-_T_DEACT / 70.0)
DEACTIVATION_BENCHMARK = DeactivationDataset(
    citation="Haas & Macdonald 1999 (J. Physiol. 514:27-45) / Jones & Westbrook 1995", ...)
```

That is a formula, not a digitisation, and the dominant 15 ms component is attributed to
the one paper `knowledge/07-paper-review.md` records as **measuring 76.1 ms** for that
quantity. The audit that found that is three commits old. `DOSE_RESPONSE_BENCHMARK` has
the same problem in a subtler form: its plateau is exactly `0.750`, the model's own fit
target, so "fitting the model to the data" is partly circular. **Fix: P1 requires every
dataset to carry a `provenance.Record`, a figure/table reference, and a digitisation
method; a trace synthesised from the model family under test is forbidden (§5.4).**

**C5. AIC/BIC only at each model's MLE.** Revision 1 said "benchmark using BIC/AIC".
`oed.py` computes them at `pam_factor=1.0` on each model's *current, unfitted*
parameters, penalises `len(m.param_names)` (every *declared* parameter, including
Model C's `pam_desens_factor`, which the pam_factor=1.0 equilibrium curve does not use),
and holds σ fixed at a default `0.03` that nothing estimated — so the entire ranking and
all Akaike weights scale with an unjustified knob. **Fix: P5 defines the protocol —
independent MLE per model, k = number of *estimated identifiable* parameters, σ either
estimated jointly or profiled out, plus held-out cross-validation, which is what the
claim actually needs.**

**C6. Falsification intervals come from the posterior, not from the noise.**
`oed.py::find_discriminating_protocol` emits `(val - 2σ, val + 2σ)` and labels it a "95%
confidence tolerance interval". That ignores parameter uncertainty entirely — the very
quantity P3/P4 exist to produce. An interval narrower than the real predictive spread
makes a model look falsifiable when it is not. **Fix: P6 requires posterior predictive
intervals and a discrimination criterion that accounts for both models' spreads.**

**C7. The headline number becomes a formula with an interval.** Revision 1 defended a
single asymptotic fold-change (Rev 1 lines 13, 85). That number is a point on the flat
manifold of C3, and it **moves 13% with the SciPy version alone**. Measured in this
container at `e2f2d16`, on the *current* config pulse, against the manuscript's own
published rate table:

```
                                        anchors reproduced          derived
manuscript rate table (§2.1 of 08-)   EC50 20.0000  Po_max 0.7500  tau 15.000  |resid| 2.4e-06
                                      -> Kd 31.95 uM   E 5.1814   HEADROOM 210.7x
refit here (scipy 1.17.1/numpy 2.4.6) EC50 20.0000  Po_max 0.7500  tau 15.000  |resid| 2.2e-09
                                      -> Kd 29.48 uM   E 4.7837   HEADROOM 182.2x
```

**Both parameter sets reproduce all three calibration anchors to the precision the
manuscript quotes.** They differ by 7.7% in K_d and gating efficacy, and by **13.5% in
the headline asymptotic dynamic range**. `fit_scheme` optimises **4 parameters against 3
residuals**, so the reported 7 significant figures are the optimiser's path, not the
model's content — and `P_open,∞/P_0` is not constant along the family, which is why a
library upgrade moves it this far.

Six tests are currently red on a fresh install for this reason alone: the five
`test_manuscript_consistency` numeric assertions, and
`test_review_regressions.py::test_the_synaptic_pulse_has_exactly_one_definition`, which
asserts `abs(Kd - 31.95) < 0.05` — a ±0.16% tolerance on a quantity that moves 7.7%
between environments. That test was written three commits ago and already fails. The code
is not wrong; the claim "this number is a property of the model" is.

**Fix: Tier 1 is reported as `P_open,∞/P_0` with its inputs and a credible interval (P6-2);
the environment is pinned so the tests mean something (P0-7); P3/P4 replace the point fit
with the identifiable reparameterisation and a posterior; §6 governs what may be claimed.**

**C8. No silent fallbacks, ever.** Four places in the new code substitute a
plausible-looking result for a failure. The worst returns the literal `15.0` — which *is*
`FIT_TARGETS["tau_ms"]` — when the decay fit fails, so a failed estimate is
indistinguishable from a perfect one. **Fix: §2.4 and P0-4.**

**C9. Tests must be able to fail for the right reason.** The ten new tests assert shapes,
finiteness, monotonicity and ordering. Not one would fail if a number were wrong by 3×.
Contrast `tests/test_manuscript_consistency.py`, which greps for the decimal
representation of every claim and is currently — correctly — red. **Fix: §2.7 anchor-test
requirement; every phase's Definition of Done names the test that would catch the error
the phase exists to prevent.**

**C10. Numerical reproducibility is a requirement, not an aspiration.** No dependency
pin, no tolerance policy, no seed policy anywhere in Rev 1. The manuscript's own
reproduction instruction (`git checkout manuscript-v3`) does not reproduce its numbers.
**Fix: P0-7.**

**C11. The circuit boundary is stated.** Rev 1 is entirely receptor-level, in a package
whose distinguishing feature is receptor → circuit → behaviour. Leaving that implicit
invites an agent to wire `charge_integral` into `resp.py` without an operating point —
the exact class of error `substrate.resolve_op` exists to refuse. **Fix: §4.P7 declares
circuit coupling OUT OF SCOPE for this programme and specifies the hand-off contract for
when it is in scope.**

---

## 2. The seven invariants

### 2.1 Every reported quantity carries a Tier

Any function whose output a human might quote returns `results.Quantity` or
`results.ResultSet`, not a bare float or a plain dataclass of floats. Internal
computation uses floats freely; the **boundary** carries tiers.

```python
# WRONG — a bare float that a reader will copy into a slide
def tonic_headroom(model) -> float: ...

# RIGHT
def tonic_headroom(model) -> Quantity:
    return Quantity(
        name="tonic_headroom", _value=v, tier=Tier.UNCALIBRATED, units="x",
        provenance="P_open,inf / P_open(ambient) at the P4 posterior median; "
                   "ambient 0.40 uM from config.AMBIENT_GABA_UM (UNSOURCED, 0.1-1 range)",
        promote_by="an ambient-GABA measurement in the preparation being modelled",
        caveats=("asymptotic: assumes k_off -> 0+, which no ligand reaches",))
```

Tier assignment rules for this programme, which you apply without asking:
* Exact algebraic identity, or reproduces a held-out observable → `VALIDATED`.
* Depends on a fitted parameter, or on an `UNSOURCED`/`GUESS` input → at best
  `UNCALIBRATED`, however tight its interval.
* Compares or combines two different observables, or rests on an input recorded invalid
  → `VOID`. Do not "fix" it by widening the error bar.

### 2.2 Every parameter and dataset carries a Basis

New numbers get a `provenance.Record` in the same commit that introduces them. Add new
tables to `provenance.ALL` so `tests/test_provenance.py` covers them. A number with no
source is `Basis.UNSOURCED` with a note saying so — that is an acceptable, honest state.
Attributing it to a paper you have not read against the specific claim is not.

> The pattern that has bitten this project twice, and is named in
> `knowledge/06-source-provenance.md`: *a source whose title matches the claim, in the
> right journal, by the right group, measuring a neighbouring quantity.* Metadata
> verification cannot catch it. Only reading the paper against the number can.

### 2.3 Every model method declares its observable

No function may compare, fit, or ratio two quantities defined on different observables.
Enforced in code from P1 onward (`models/base.Observable`). Three observables exist in
this programme and they are not interchangeable:

| Observable | Meaning | Where it is right |
|---|---|---|
| `PEAK` | max P_open during an agonist application | concentration-response curves from fast perfusion |
| `EQUILIBRIUM` | stationary P_open, desensitisation included | tonic/extrasynaptic current |
| `CHARGE` | ∫P_open dt over a defined window | synaptic charge transfer, the circuit bridge |

The legacy `gabaa_kinetics.Scheme.ec50_um()` is `PEAK` (it maxes `po_peak` over a 300 ms
window). The new `models.*.dose_response()` is `EQUILIBRIUM`. **They differ by a factor of
~3 at the same parameters.** That is not a bug in either; it is two observables, and the
bug is comparing them.

### 2.4 Failure raises; it never substitutes

A numerical failure (solver non-convergence, a fit that will not converge, an empty
analysis window, a NaN from an unreachable calibration) **raises with a message naming
the cause**, or returns a `Quantity` with `Tier.VOID` and the reason. It never returns a
default, a target value, a clamped value, or the input unchanged.

```python
# WRONG — the project has shipped all four of these
states = sol.y.T if sol.success else np.tile(p_init, (len(t_arr), 1))
p_open = sol.y[0] if sol.success else target_p
decay_tau = ... if np.any(half_mask) else 15.0        # 15.0 IS the fit target
cost = res.fun if res.success else base_cost + 10.0   # failure scored as "identifiable"

# RIGHT
if not sol.success:
    raise RuntimeError(
        f"Radau failed to integrate the 5-state scheme over "
        f"[{t_arr[0]}, {t_arr[-1]}] ms: {sol.message}. The stiffest timescale here is "
        f"1/r = {1/self.r:.0f} ms against a cleft transient of ~1 ms; if this is a "
        f"legitimate parameter set, raise t_eval density rather than loosening rtol.")
```

**Never choose a sentinel that is also a plausible correct answer.** If you need one, use
`float("nan")` and let it propagate to a `VOID` quantity.

### 2.5 Parameters come from a registry, never from literal defaults

After P4 closes, every dataclass default in `models/` is marked
`# PROVISIONAL - not a calibration; see parameters.py` and production code obtains
parameters from `circuitpharm.parameters.get(...)`. Literal defaults remain legal in unit
tests that are testing algebra, and nowhere else. The reason is C1: three uncalibrated
default sets are already in the tree and are already producing headline numbers.

### 2.6 Backwards compatibility of the validated path

`gabaa_kinetics.Scheme`, `cpg.Pop`, `resp.PreBotC(substrate="lif")` and everything
downstream of them are **byte-identical-output** surfaces. The project's one
`VALIDATED` result (the calibration-independent selectivity ranking) and the manuscript
both depend on them. You may add to them; you may not change what they return. If a
refactor would change a number on the LIF path, stop and report instead.

### 2.7 Every phase ships at least one anchor test

An **anchor test** asserts a *value*, against something the code was not fitted to, with
a tolerance you justify in a comment. A test that asserts `isfinite`, `> 0`, `shape ==`,
or a monotone ordering is a **mechanics test**. Mechanics tests are welcome and
insufficient. Each phase's Definition of Done names its anchor tests explicitly; a phase
without one does not close.

```python
# mechanics test — necessary, not sufficient
assert np.isfinite(cost) and cost >= 0.0

# anchor test — this is what the phase is judged on
def test_equilibrium_ec50_matches_closed_form():
    """[G]_1/2 = Kd(1+sqrt(2+E))/(1+E) is exact for the 5-state scheme with no
    desensitisation. Tolerance 1e-6 relative: both sides are analytic, so anything
    looser would hide an algebra error in the Q matrix."""
    m = KineticAllosteryModel(d=0.0, r=1.0)      # desensitisation off
    E, Kd = m.gating_efficacy, m.kd_um
    assert np.isclose(ec50_numeric(m), Kd*(1+np.sqrt(2+E))/(1+E), rtol=1e-6)
```

---

## 3. P0 — Defect register. Close this before any phase work.

Thirteen defects, all verified in this container at commit `e2f2d16` unless marked
STATIC (read, not run). Each has a required fix and the test that must accompany it.
Fix them in one branch, one commit per numbered item, in this order.

### P0-1 — `len()` on a scalar initial state crashes both Markov models
**Files:** `models/kinetic_jw95.py:161-165`, `models/extended_desens.py:158`
**Severity:** high (documented API, hard crash)

`models/base.py:83` documents `initial_state` as "Initial state vector **or initial open
probability**". Reproduced:

```
B: TypeError: len() of unsized object
C: TypeError: object of type 'float' has no len()
```

In B, `np.asarray(0.1)` is 0-d and `len()` raises before the `!= 5` test; in C, `len()` is
called on the raw float.

**Required fix.** One shared helper, because two copies of this logic is how C3-class
divergence starts:

```python
# models/base.py
def resolve_initial_state(initial_state, n_states: int, resting: np.ndarray) -> np.ndarray:
    """Coerce a user initial condition to a full state vector on the simplex.

    Accepts None (use `resting`), a scalar open probability (place it in the open state
    and distribute the remainder across `resting`'s non-open mass, renormalised), or a
    full length-`n_states` vector (validated: non-negative, sums to 1 within 1e-8).
    Anything else raises ValueError naming what was passed and what is accepted.
    """
```

A scalar must *not* be silently replaced by the resting distribution, which is what
`kinetic_jw95.py:165` does today — that discards the caller's instruction without a word.

**Tests:** `tests/test_nextgen_models.py` — scalar accepted and honoured
(`p_open[0] == pytest.approx(0.1)`); 5-vector accepted; 4-vector raises `ValueError`;
a vector summing to 1.3 raises `ValueError`.

### P0-2 — `rise_ms == clear_ms` silently produces a zero GABA transient
**File:** `protocols/waveforms.py:39-47`
**Severity:** high (silent wrong input to every downstream simulation)

Reproduced:
```
rise=1.0   clear=1.0 : peak = 0 uM      (asked for 1000)
rise=0.999 clear=1.0 : peak = 1000 uM
rise=0.9   clear=1.0 : peak = 998.3 uM
```
`raw = exp(-t/c) - exp(-t/r)` is identically zero when `c == r`; `norm_factor` is then
clamped to `1e-6` and every sample is `0`. A receptor model handed an all-zero agonist
trace returns a flat resting trace, a `charge_integral` near zero, and a `decay_tau_ms` of
`15.0` (P0-4) — three plausible numbers from an empty experiment.

**Required fix.** Implement the analytic `c → r` limit, which is the alpha function:

```
raw(t) = (t / r) * exp(-t / r)        normalised to 1 at its peak t = r
```
Switch on `abs(c - r) < 1e-9 * max(c, r)`, not `c == r`. Add a comment recording that the
near-equal case (0.999 vs 1.0) is numerically fine, so a reader does not "simplify" the
guard into an equality test again.

**Tests:** peak within 0.1% of `peak_um` for `clear_ms/rise_ms` in
`{1.0, 1.0+1e-12, 1.001, 1.1, 10.0}`; continuity — the peak and the charge integral of
the `c = r` waveform are within 0.5% of those at `c = r*(1+1e-6)`.

### P0-3 — solver failure substitutes a plausible trace
**Files:** `models/kinetic_jw95.py:199`, `models/extended_desens.py:199`,
`models/operational.py:109`
**Severity:** high (STATIC — the failure path is not reachable at default parameters, which
is why it has not bitten yet)

Models B and C tile `p_init` into a constant trace; Model A substitutes `target_p`, the
instantaneous-Hill limit, which **deletes the relaxation dynamics the model exists to
express** and biases peak and charge upward.

**Required fix.** `raise RuntimeError` per §2.4, with `sol.message`, the integration span
and the stiffest timescale in the message. No `try`/`except` around the raise anywhere
upstream.

**Tests:** monkeypatch `solve_ivp` to return `success=False` and assert
`pytest.raises(RuntimeError)` for all three models; assert the message contains
`sol.message`.

### P0-4 — `decay_tau_ms` returns the fit target on failure
**Files:** `models/kinetic_jw95.py:216-222`, `models/extended_desens.py:216-220`
(literal `15.0`), `models/operational.py:123-127` (`self.tau_deact_ms`)
**Severity:** high (a failed estimate is indistinguishable from a perfect fit)

`15.0` is `gabaa_kinetics.FIT_TARGETS["tau_ms"]`. Three separate branches return it. Any
check of the form "does the model reproduce the 15 ms IPSC decay?" passes on failure.

**Required fix.** Return `float("nan")` from every failure branch. Additionally, replace
the single-point half-decay estimate with a weighted biexponential fit, because §4.P6 and
the technical spec both require `tau_fast`, `tau_slow`, `tau_weighted`:

```python
@dataclass(frozen=True)
class DecayFit:
    tau_fast_ms: float; tau_slow_ms: float; weight_fast: float
    tau_weighted_ms: float; r_squared: float; n_points: int
    # every field NaN, r_squared 0.0, when the fit does not converge or the tail is flat
```
Fit on the 90%→10% window of the post-peak tail (the standard experimental window, and
what `gabaa_kinetics.ipsc_metrics` already uses — reuse that convention, do not invent a
second one). Report `r_squared` so a caller can see a bad fit instead of inferring it.

**Tests:** anchor — a synthetic `0.7*exp(-t/15) + 0.3*exp(-t/70)` trace recovers
`tau_weighted` within 2% and both components within 5%; a flat trace yields all-NaN with
`r_squared == 0.0`; `np.isnan` on the failure path of all three models.

### P0-5 — the fitting objective cannot see the absolute rates
**File:** `fitting/identifiability.py:25-46`
**Severity:** critical (invalidates every identifiability and MCMC result produced so far)

Reproduced: scaling all six rates by 2, 10 and 100 changes the cost by **zero** at 6
decimal places. The objective's only terms are the equilibrium dose-response (a function
of K_d, E, D alone) and a `po_max` penalty (a function of E alone). Consequences, all of
which are currently being reported as results:
* `compute_profile_likelihood` must return "non-identifiable" for `kon`, `koff`, `beta`,
  `alpha`, `d`, `r` individually — a property of the objective, not the data.
* `compute_fisher_information_matrix` is structurally rank ≤ 3, and its condition number
  is computed over *positive eigenvalues only* (`identifiability.py:166-167`), silently
  dropping the zeros that are the actual finding.
* `run_ensemble_mcmc` samples a posterior that is flat along three directions, bounded
  only by `log_prior`'s box. The reported "95% credible intervals" are prior widths.

**Required fix (P0 part — the full treatment is P3).** Do not attempt the real fit here.
In P0: (a) raise `NotImplementedError` from `compute_profile_likelihood`,
`compute_fisher_information_matrix` and `run_ensemble_mcmc` with a message pointing at
this register entry and at P3; (b) keep `compute_kinetic_objective` but rename it
`equilibrium_dose_response_chi2` and document in one line that it constrains only
`(Kd, E, D)`; (c) delete or `@pytest.mark.skip(reason="blocked on P3")` the three tests
that assert those functions return plausible values. Shipping a disabled analysis is
honest; shipping a flat one as a result is not.

**Tests:** anchor — `test_objective_is_invariant_to_uniform_rate_scaling` asserts the
cost is unchanged to `rtol=1e-12` under ×100, with a docstring explaining that this is the
*structural* reason the analysis is blocked. That test is the permanent record of why P3
exists; it must not be deleted when P3 lands, only amended.

### P0-6 — the deactivation "benchmark" is model-generated and mis-cited
**File:** `fitting/data.py:43-56`
**Severity:** critical (scientific integrity; a repeat of a defect retracted 3 commits ago)

`_I_DEACT` is `0.70*exp(-t/15) + 0.30*exp(-t/70)`, cited to Haas & Macdonald 1999 —
recorded in `knowledge/07-paper-review.md` as measuring **76.1 ms** for this quantity, one
of the six citations that audit found did not support the number it carried. The
`DOSE_RESPONSE_BENCHMARK` plateau is exactly `0.750`, the model's own fit target, which
makes the fit partly circular; its `sem` values have no stated origin.

**Required fix.** (a) Rename both to `*_SYNTHETIC` and add a module-level docstring
stating in the first line that these are **synthetic traces for mechanics testing, not
data**; (b) register both in `provenance.ALL` as `Basis.UNSOURCED` with the note
explaining the construction; (c) add the `Observable` tag (`PEAK` for the
concentration-response, `CHARGE`/current-decay for the trace) so P1's contract applies;
(d) remove the `citation=` strings or move them to a `motivated_by=` field whose name
cannot be mistaken for provenance. Real digitisation is P1-2 and must not be faked here.

**Tests:** `tests/test_provenance.py` must fail if either dataset loses its record;
a test asserting the word "synthetic" appears in the module docstring (blunt on purpose —
the same technique `test_manuscript_consistency.py` uses, for the same reason).

### P0-7 — the environment is not pinned, so the manuscript does not reproduce
**Files:** `pyproject.toml`, `.github/workflows/tests.yml`,
`knowledge/08-manuscript.md`
**Severity:** high (6 currently-failing tests, a reproducibility claim that does not hold,
and CI that has not completed in six pushes)

Measured above (C7): identical anchors, K_d 31.95 vs 29.48 µM, headline range 210.7× vs
182.2×, from a SciPy upgrade alone. `pyproject` allows `numpy>=1.26`, `scipy>=1.11`; this
container resolves to numpy 2.4.6 / scipy 1.17.1. Full current state of the suite on a
lean install here: **6 failed, 270 passed, 19 skipped** — the five
`test_manuscript_consistency` numeric tests plus
`test_the_synaptic_pulse_has_exactly_one_definition`.

**Required fix.**
1. Add a `[project.optional-dependencies] repro` extra pinning the exact
   numpy/scipy/python used to generate the manuscript, and state in `08-manuscript.md`'s
   reproduction section that `pip install -e '.[repro]'` is required — the tag alone is
   insufficient.
2. In CI, add a `repro` job that installs that extra and runs
   `pytest tests/test_manuscript_consistency.py`. The lean and full jobs continue to run
   unpinned, so drift is *visible* rather than silent.
3. Fix `actions/checkout@v4` to fetch tags (`fetch-tags: true` or `fetch-depth: 0`).
   `test_the_cited_commit_and_tag_exist` greps for `manuscript-v3` and runs
   `git cat-file -t` on it; the default checkout fetches no tags, so this test fails in CI
   for an environment reason. Verified: the tag exists on the remote and points at
   `02e684a`; `git fetch --tags` makes the test pass locally.
4. Separately, P3 must decide the real answer to C7 — a non-identifiable fit should not be
   quoted to 7 significant figures at all. Pinning makes the number stable; it does not
   make it meaningful. Both are needed.

**And a separate CI defect found while checking this one: CI has not completed on this
repository for at least the last six pushes.** Every run from `b7f78a5` to `02e684a`
reports `conclusion: cancelled`, which is what GitHub reports for a job that hits its
`timeout-minutes`. On run 19 (`02e684a`, the commit before the sprint):

```
lean (3.11)  "Fast tests"              04:02:00 -> 04:32:01   30m01s  CANCELLED (limit 30)
lean (3.12)  "Fast tests"              04:02:03 -> 04:31:58   29m55s  CANCELLED (limit 30)
full         "Full suite with coverage" 04:02:13 -> 04:46:59   44m46s  CANCELLED (limit 45)
```

The `Circuit integration tests` step never ran in either lean job, and
`Assert nothing was skipped` never ran in `full`. So **nobody has seen this suite's real
CI verdict since 2026-10-08**, and the `--cov-fail-under=90` floor and the
no-skips assertion have not been enforced since then either. Note what that means for the
sprint: the 3,177 new lines were merged without a completed CI run.

The `-m "not slow"` subset alone exceeding 30 minutes is the signal — it is supposed to be
the fast one. **Required fix, in this order:** (a) measure it —
`pytest -m "not slow" -q --durations=25` — and report the table; (b) mark whatever
dominates it `@pytest.mark.slow` if it is genuinely circuit-scale, or make it cheap if it
is not (the new `fitting/` tests run optimisers and samplers inside the fast subset and are
unmarked); (c) only then adjust `timeout-minutes`, and raise it because the measured
runtime justifies the number, not to make the red go away. Do not raise the timeout first:
a 30-minute "fast" subset is a defect in the subset, and the timeout is the only thing
currently reporting it.

### P0-8 — `Drug.from_kinetics` is missing the NaN guard that `pool_gains` has
**File:** `cpg.py:88-113`
**Severity:** medium (hard crash with a misleading message, on the legacy path)

Reproduced:
```
Drug.from_kinetics(ec50_shift=1.0)
 ** On entry to DLASCL parameter number  4 had an illegal value
 LinAlgError: SVD did not converge in Linear Least Squares
```
Same for `ec50_shift=1e6`. `calibrate_pam` returns NaN for any unreachable shift (≤1.0 is
a neutral ligand or a NAM — out of domain by design); `evaluation.py:146-154` catches
exactly this and raises a readable error whose comment predicts the LAPACK message above.
`from_kinetics` is the second public entry to the same machinery and never got the guard.
This is recurring error E12 (a fix applied everywhere but one module).

**Required fix.** Extract the check into `gabaa_kinetics.require_reachable(aff, s_max,
ambient_um)` and call it from both sites. Do not copy the error text.

**Tests:** `pytest.raises(ValueError, match="not reachable")` for `ec50_shift` in
`{0.5, 1.0, 1e6}`; `Drug.from_kinetics()` with defaults still returns
`gaba_a_gain ≈ 1.0495`, `gaba_a_gain_tonic ≈ 7.9377` (anchor — these are the current
values on the unified `config.SYNAPTIC_PULSE`).

### P0-9 — `SpinalCircuit` accepts a conductance operating point with no weights
**File:** `circuit.py:83-91`
**Severity:** medium (latent; `COND_SPINAL_OP is None` today, but `op=` is public and
anchoring the spinal circuit is the documented next step)

It passes `required_w=()` to `substrate.resolve_op`. Reproduced with a scalars-only `op`:

```
('RG','PF','ampa') 2.2    ('PF','Mn','ampa') 2.2    ('InPF','PF*','gly') 3.5
cell g_L = 2.8 nS; LIF g_L = 10 nS
```

The circuit runs the LIF-tuned nS table on the Butera cell — precisely the partial
operating point `resolve_op`'s docstring calls "the dangerous case, not the missing one".
`resp.py` requires all 8 of its weights; `rg2.py` all 4 of its; `circuit.py` requires none.

**Required fix.** Pass the full tuple of weight keys the circuit actually uses (the 12
`(pre, post, rec)` keys in `SpinalCircuit.W`). `resolve_op`'s existing `k not in op["w"]`
check handles tuple keys unchanged.

**Tests:** `pytest.raises(ValueError, match="incomplete")` for a scalars-only cond `op`;
the LIF path unaffected (§2.6).

### P0-10 — `evaluate()` returns a different result shape on a lean install
**File:** `evaluation.py:416-420`
**Severity:** low (STATIC; the branch needs `assays` itself to fail to import)

The `except Exception` branch `return rs` early, before the `VOID overdose_index` and
`nmda_respiratory_contribution` entries are added. `tests/test_evaluate.py:110` asserts
`overdose_index` is always present. A package whose thesis is that the VOID entries always
appear must not have an install-dependent result shape.

**Required fix.** Replace the early `return` with a flag so control reaches the VOID
block; or hoist the VOID additions above the motor section. Add
`test_void_quantities_present_without_plant`.

### P0-11 — `WaveformResult.current_pA()` returns all zeros with its documented defaults
**File:** `models/base.py:21-24`
**Severity:** low
`v_hold_mv=-70.0, e_rev_mv=-70.0` → zero driving force → identically zero current.
Reproduced: `max |I| = 0 pA`. `protocols/waveforms.extract_electrophys_metrics` works
around it with a separate `driving_force_mv=30.0` argument, so there are now two
conventions for one quantity.

**Required fix.** Make `driving_force_mv` a required argument of one function and delete
the other path, or default `v_hold_mv=-60.0` with a comment giving the physiological
basis. State the chloride reversal convention once, in `models/base.py`, and reference
`cpg.E_REV["gabaa"] = -75.0` so the receptor layer and the circuit layer cannot drift.

### P0-12 — smaller items, fix in one sweep commit
* `models/extended_desens.py:103` — `d_fast / (1 + 0.1*(pf-1)*pam_desens_factor)` has an
  unguarded denominator. With `pam_desens_factor` free to a fitter (it is in
  `param_names`), a modest negative or large positive value drives the denominator to
  zero or negative, giving an infinite or **negative** rate constant. Clamp as
  `cpg.Drug.nmda_scale` and `substrate.tonic_gaba` already do, and say in the comment that
  those two clamps are the precedent.
* `models/*.apply_pam` — `max(float(pam_factor), 1e-6)` admits `pam_factor < 1`, i.e. a
  NAM, which the surrounding code is not written for (`decomposition.py:52` separately
  clamps to `>= 1.0`, so the two disagree). Decide once: either support NAMs end to end or
  raise below 1.0. Document the decision in `models/base.py`.
* `fitting/identifiability.py:84` — `cost = res.fun if res.success else base_cost + 10.0`
  scores an optimiser **failure** as evidence of identifiability (cost above the 3.84
  threshold). Sign-inverted. Raise instead (§2.4).
* `fitting/identifiability.py:65,88` — `base_cost` is the cost at the *unoptimised*
  `base_model`, so `delta_costs` can be negative and the CI logic is not a profile
  likelihood. Re-optimise the unconstrained problem first in P3.
* `fitting/mcmc.py:26-41` — `log_prior` is documented "log-uniform" and returns `0.0`
  inside the box, which is uniform on the rates. For rate constants spanning decades the
  difference is material. Either return `-sum(log(theta))` or fix the docstring; prefer
  sampling in `log10` space (P3 does).
* `chirality.py:64` — `hasattr(Chem, "EnumerateStereoisomers")` is always False (the
  submodule is not imported), so `iso` is always `None` and is never read. Delete.
* `resp.py:37` — re-declares `EUPNOEA_BAND`, duplicating `config.EUPNOEA_BAND`, against
  the single-source rule `config.py` exists to enforce. Import it.
* `resp.py:270` — `band` (the validity tuple parameter) is rebound to the FFT boolean
  mask. Rename the mask.
* `config.py:33-94` — two copies of the `COND_RESP_OP` header comment; the first still
  says `None means NOT YET ANCHORED` while the value below it is anchored.
* `evaluation.py:103` — the `pool_gains` cache key omits `modality`, so changing a
  compound's modality after a first call silently returns the earlier mechanism's gains.

### P0-13 — the model defaults are a chimera of two fits; label them before P4
**Files:** `models/kinetic_jw95.py:24-29`, `models/extended_desens.py:23-30`,
`models/operational.py:19-24`
**Severity:** medium now, critical the moment anyone quotes a number from them

Verified (C1): `KineticAllosteryModel`'s `beta`/`alpha` are the manuscript's current
config-pulse fit; its `kon`/`koff` are approximately the superseded 1000 µM / 0.30 ms fit.
K_d = 29.62 µM, a value from no fit, no commit and no document. `ExtendedDesensitizationModel`'s
defaults are round numbers with no stated origin at all (`kon=0.012, koff=0.35,
beta=0.60, alpha=0.10`).

**Required fix (P0 part only — the real fix is P4).** Add to each dataclass a
module-level comment and a `# PROVISIONAL` marker on every default, stating: these are not
a calibration, they reproduce neither the EC₅₀ nor the P_o,max anchor, the measured values
are EC₅₀ 6.34 µM / P_o,max 0.8382 (B) and 5.18 µM / 0.8571 (C), and P4 replaces them.
Then make `evaluate_dynamic_range`, `compare_phasic_tonic_mechanisms` and
`decompose_modulation_gain` require an explicit `model` argument — delete the
`model or KineticAllosteryModel()` fallbacks (`dynamic_range.py:51`,
`decomposition.py:117`) so no caller can reach the defaults by accident.

**Tests:** **ANCHOR** `test_model_defaults_are_not_the_anchors` asserts the measured
EC₅₀/P_o,max of each default set and documents in its docstring that this test exists to
prevent the defaults being mistaken for a calibration — it should be *amended*, not
deleted, when P4 lands. Plus `pytest.raises(TypeError)` for each entry point called with
no model.

**P0 Definition of Done.** All thirteen closed; `pytest -q` shows **zero failures** other
than any explicitly `skip`-marked in P0-5; the three new anchor tests
(P0-2 continuity, P0-4 biexponential recovery, P0-5 scaling invariance) present and
passing; `knowledge/09-strategic-roadmap.md` §8 updated to `P0: CLOSED <sha>`.

---

## 4. Phases P1 – P7

Each phase states: **Goal**, **Why here** (the dependency), **Deliverables** (files and
signatures), **Algorithm**, **Acceptance** (named tests, with the anchor tests marked),
**Definition of Done**, and **Forbidden** (shortcuts that have been tried).

Phases are sequential. P2 may run in parallel with P1. Nothing else may.

---

### P1 — The observable contract and real data

**Goal.** Make it impossible to fit or compare across observables, and replace the
synthetic benchmarks with digitised data carrying provenance.

**Why here.** C2 and C4. Every downstream phase fits, compares or ratios these
quantities. A χ² of 9290 from an observable mismatch is not a fitting problem; it is a
contract problem, and no amount of optimiser work fixes it.

#### P1-1 The `Observable` enum and model declarations

```python
# models/base.py
class Observable(str, Enum):
    PEAK = "PEAK"                # max P_open during an application; fast-perfusion CRCs
    EQUILIBRIUM = "EQUILIBRIUM"  # stationary P_open incl. desensitisation; tonic current
    CHARGE = "CHARGE"            # integral P_open dt over a declared window

@dataclass(frozen=True)
class ObservedQuantity:
    """A number plus the observable it is a number OF. Arithmetic between two of these
    with different `observable` raises TypeError."""
    value: float
    observable: Observable
    window_ms: float | None = None   # required and > 0 when observable is CHARGE
```

Every `ReceptorModel` gains:

```python
@property
def native_observable(self) -> Observable: ...
def dose_response(self, concs_um, pam_factor=1.0, *,
                  observable: Observable = Observable.EQUILIBRIUM) -> np.ndarray: ...
def peak_dose_response(self, concs_um, pam_factor=1.0,
                       application_ms: float = 300.0) -> np.ndarray: ...
```

`peak_dose_response` integrates a square application and takes the max — the same protocol
as `gabaa_kinetics.Scheme.po_peak`, whose docstring records why 300 ms and not 1 ms (at low
agonist a 1 ms step is binding-rate limited, which inflates the fitted EC₅₀ and
compresses every derived PAM gain). **Reuse that constant; do not pick a new one.**

Model A is `PEAK`-native with no desensitisation; its `EQUILIBRIUM` request must raise
`NotImplementedError` with a one-line explanation, not silently return the Hill curve.
That asymmetry *is* a model difference and P5 must see it.

#### P1-2 Digitised datasets with provenance

```python
# fitting/data.py
@dataclass(frozen=True)
class Dataset:
    key: str
    observable: Observable
    source_key: str              # key into provenance; must resolve
    figure: str                  # e.g. "Fig 2B, filled circles"
    digitisation: str            # tool + method + date + who
    preparation: str             # subunit composition, expression system, temperature
    x: np.ndarray                # concentration (uM) or time (ms)
    y: np.ndarray
    y_err: np.ndarray            # SEM or SD — say which in `digitisation`
    n_cells: int | None
    window_ms: float | None = None
```

Required datasets, minimum three, each a *separate* published source:
1. **Concentration–response, `PEAK`.** α1β2γ2, recombinant, fast perfusion. Must report
   the figure digitised and the number of cells.
2. **Deactivation time course, `CHARGE`.** A real current trace after a brief saturating
   pulse. This is the dataset that makes the absolute rates identifiable (P3) — without a
   kinetic timescale in the likelihood, P0-5's flat manifold persists no matter what else
   you add.
3. **A held-out set, never used in any fit.** Candidates: paired-pulse recovery at a
   stated interval, a desensitisation onset time course, single-channel mean open time, or
   a PAM concentration-response at a *different* ambient GABA. P5's cross-validation and
   P6's falsification both consume this. Reserve it in code —
   `Dataset.role: Literal["train", "holdout"]` — and add a test that no fitting function
   can reach a `holdout` row.

**If you cannot obtain a real digitisation for one of these, stop and report.** Do not
synthesise. An honest P1 with two real datasets and a written gap is worth more than three
fabricated ones, and §5.4 makes fabrication a revert.

#### P1-3 Acceptance

* **ANCHOR** `test_peak_and_equilibrium_ec50_differ_as_expected`: for the P4 parameter
  set, `peak_dose_response` EC₅₀ ≈ 20 µM and `dose_response` (EQUILIBRIUM) EC₅₀ ≈ 6 µM,
  both within 5%. The test's docstring states that these are two observables and that
  their ratio is a prediction of the scheme, not an error.
* **ANCHOR** `test_closed_form_equilibrium_ec50`: as in §2.7, `rtol=1e-6`, desensitisation
  off.
* `test_cross_observable_arithmetic_raises`: `ObservedQuantity(PEAK) / ObservedQuantity(EQUILIBRIUM)`
  raises `TypeError`.
* `test_every_dataset_has_a_resolvable_provenance_record` — iterate `fitting.data.ALL`.
* `test_holdout_unreachable_from_fitting` — introspect the fitting entry points.
* `test_model_a_equilibrium_raises`.

**Done when** no fitting or comparison function in the tree can be called with two
different observables, and `fitting/data.py` contains no expression that generates `y`
from a formula.

**Forbidden.** Adding an `observable=` argument that defaults to something and is never
checked. Tagging the synthetic traces `PEAK` and calling P1 closed.

---

### P2 — Tier and provenance integration for the new layer

**Goal.** Every quantity the new layer reports comes back as a `Quantity` with a tier,
and every new parameter has a `Basis`. May run in parallel with P1.

**Why here.** §2.1 and §2.2. This is the package's identity and 3,177 lines currently
bypass it. It is also cheap, and doing it now means P3–P6 are written tier-aware instead
of being retrofitted.

**Deliverables.**
* `dynamic_range.evaluate_dynamic_range() -> ResultSet` with three entries:
  `theoretical_asymptotic_headroom` (**UNCALIBRATED** until P4, with the caveat that
  `k_off → 0⁺` is unreachable by any ligand), `reachable_gain_at_smax_2.5`
  (**UNCALIBRATED**), `physiological_charge_ratio` (**UNCALIBRATED**, and
  `window_ms` stated in its provenance because the ratio depends on it).
* `protocols/decomposition.py` → `ResultSet`. The Shapley percentages are
  **UNCALIBRATED** at best; see P4-3 for why the current normalisation is not a Shapley
  value and must be renamed until it is one.
* `protocols/oed.py` → every AIC/BIC/weight is **VOID** until P5 closes, with
  `provenance` naming C5 and `promote_by` naming P5. This is the correct tier today: they
  are computed at unfitted parameters, so they rest on something known to be invalid.
* `fitting/` → `NotImplementedError` per P0-5 until P3.
* `provenance.ALL` gains `MODEL_PARAM_PROV` and `DATASET_PROV` tables.

**Acceptance.**
* **ANCHOR** `test_oed_metrics_are_void_until_fitted`: reading `.value` on an AIC quantity
  raises `VoidQuantityError`, and the message names the unfitted-parameter reason.
* `test_every_new_public_result_is_a_quantity` — walk the public API of `models`,
  `protocols`, `fitting`, `dynamic_range`; any function returning a bare `float` whose
  name is not prefixed `_` fails the test, with an explicit allowlist for pure algebra
  helpers (`kd_um`, `po_max`, …).
* `test_provenance_covers_model_params`.
* `test_str_of_resultset_withholds_void` — the printed form must not contain the number.

**Done when** `provenance_report()` includes the new tables and no new public function
returns an untiered quoteable number.

**Forbidden.** Marking anything `VALIDATED` in this phase. Nothing in the new layer has
reproduced a held-out observable yet; that is P5's job at the earliest.

---

### P3 — Identifiability, done as a question

**Goal.** Determine *what the data can determine*, before fitting anything. Produce the
identifiable reparameterisation and a written classification of every direction.

**Why here.** C3/P0-5. Fitting 6 rates to data that constrain 3 combinations produces a
number with a false interval, which is worse than no number. The project has already
published one of those (`prebotc_gaba_sens`, 18 of 35 grid cells fitting the anchor
equally well, now `VOID`).

#### P3-1 Structural identifiability, analytically, first

Write the equilibrium and peak observables in terms of `(K_d, E, D)` and prove on paper —
in `knowledge/11-identifiability.md` — which groupings each observable determines. The
result is already half-known: `dose_response` depends only on `(K_d, E, D)`; a
deactivation time course additionally determines an absolute timescale (`k_off`, hence
`k_on` via `K_d`); `β` and `α` separate only through a kinetic observable that resolves
gating (rise time, single-channel open time), not through `E` alone.

Deliverable: `fitting/reparam.py`

```python
@dataclass(frozen=True)
class IdentifiableParams:
    """The combinations the data can determine, with the map to microscopic rates.

    log10_kd, log10_E, log10_D are determined by an equilibrium or peak CRC.
    log10_koff is determined only when a kinetic timescale is in the likelihood.
    beta and alpha separate only with a gating-resolving observable; absent one,
    `from_microscopic` is many-to-one and `to_microscopic` requires `koff_fixed`.
    """
    log10_kd: float; log10_E: float; log10_D: float
    log10_koff: float | None = None
```
with `to_microscopic()` / `from_microscopic()` and a round-trip test.

#### P3-2 Profile likelihood, correctly

```python
def profile_likelihood(param: str, data: Sequence[Dataset], *,
                       n_grid: int = 25, span_decades: float = 1.5,
                       optimiser: str = "L-BFGS-B") -> ProfileLikelihoodResult
```
Rules, each of which fixes a specific defect in the current implementation:
1. **Re-optimise the unconstrained problem first.** `base_cost` is the global MLE cost,
   not the cost at a nominal point (`identifiability.py:65`).
2. **Grid in log space**, ±`span_decades` around the MLE, `n_grid ≥ 25`. Seven
   multiplicative factors spanning 0.3–3.0 cannot locate a 95% boundary; the CI it
   reports is a grid artefact.
3. **Warm-start each grid point from the previous point's solution.** Profile likelihoods
   are continuous; restarting from the nominal point at every node is how you get a
   non-monotone profile that looks like noise.
4. **Optimiser failure raises** (P0-12), never scores as identifiable.
5. **Classify in three ways, not two**: `IDENTIFIABLE` (profile exceeds Δχ² = 3.84 on
   both sides within the span), `PRACTICALLY_NON_IDENTIFIABLE` (exceeds on one side —
   report which, and the bound that exists), `STRUCTURALLY_NON_IDENTIFIABLE` (flat to
   `rtol=1e-8` across the span — assert this analytically from P3-1 as well; a numerical
   flat is evidence, the algebra is proof).
6. **Interpolate the threshold crossing** between grid nodes; do not return a node.

#### P3-3 The Fisher information matrix, honestly

Rename to `cost_hessian_and_spectrum`. Return **all** eigenvalues and eigenvectors, the
rank at a stated tolerance, and the condition number over the retained spectrum with the
tolerance recorded. Never filter to positive eigenvalues (`identifiability.py:166`) — away
from a minimum, negative eigenvalues are the finding. Report each near-null eigenvector in
the `(log K_d, log E, log D, log k_off)` basis, so a flat direction is *named*
("k_on and k_off scale together") rather than left as a number.

**Acceptance.**
* **ANCHOR** `test_structural_flatness_is_detected`: on the equilibrium-only likelihood,
  the uniform rate-scaling direction is classified `STRUCTURALLY_NON_IDENTIFIABLE`, and
  the Hessian rank is exactly 3 at `tol=1e-8`.
* **ANCHOR** `test_koff_becomes_identifiable_with_a_kinetic_dataset`: adding the P1-2
  deactivation dataset moves `log10_koff` to `IDENTIFIABLE` with a finite CI. This is the
  phase's whole claim; if it fails, the deactivation dataset is not constraining what you
  think it is.
* `test_reparam_round_trip` (`rtol=1e-12`).
* `test_profile_is_monotone_away_from_the_mle` (warm-start regression).
* `test_optimiser_failure_raises`.

**Done when** `knowledge/11-identifiability.md` states, per parameter, which class it is
in and under which dataset combination, and the code's classification agrees with the
algebra for every case the algebra covers.

**Forbidden.** Reporting a CI for a `STRUCTURALLY_NON_IDENTIFIABLE` direction. Fitting all
six rates because the sampler happens to run.

---

### P4 — The fit, the posterior, and the parameter registry

**Goal.** One MLE and one posterior per model, over the identifiable parameters only,
from the P1 training data. Publish them as the single source of model parameters.

**Why here.** Everything downstream quotes numbers; they must come from here (C1, §2.5).

#### P4-1 Likelihood

```python
# fitting/likelihood.py
def log_likelihood(theta: IdentifiableParams, model_cls, data: Sequence[Dataset],
                   *, sigma: float | None = None) -> float
```
* Each dataset contributes on **its own** observable, through the matching model method.
* Weight by the dataset's own `y_err`. If `sigma` is `None`, estimate a per-dataset scale
  jointly (one extra parameter per dataset, counted in `k` for P5) — do not hard-code
  `0.03` (C5).
* Normalisation constants included, because P5 compares likelihoods across models.
* `holdout` datasets raise if passed (P1-2).

#### P4-2 MLE and posterior

* MLE: multi-start (≥ 32 Latin-hypercube starts in log space), `L-BFGS-B` then
  Nelder–Mead polish, and **report the spread of the top 5 optima**. If they disagree
  beyond the P3 CIs, the problem is not solved — say so rather than taking the best.
* Posterior: keep the affine-invariant sampler, but fix it:
  - sample in `log10` space (removes the prior mismatch of P0-12 and the positivity
    constraints at once);
  - priors stated explicitly per parameter in a table, each with a `Basis` record;
  - `n_walkers >= 4 * dim`, `n_steps` set by convergence, not by a literal — current
    defaults are 16 walkers × 150 steps, which for 6 dimensions is far short;
  - convergence diagnostics are **mandatory outputs**: integrated autocorrelation time
    per parameter, chain length in multiples of it (require ≥ 50), split-R̂ (require
    < 1.01), acceptance fraction in 0.2–0.5. A run that fails any of these returns a
    `VOID` result naming the failing diagnostic — it does not return a credible interval.
* Save chains to `data/posterior_<model>_<sha>.npz` with the git sha, the dataset keys and
  the dependency versions in the archive.

#### P4-3 The parameter registry

```python
# parameters.py  (new, top-level in the package)
def get(model: str, which: Literal["mle", "posterior_median", "nominal"] = "mle"
        ) -> tuple[dict[str, float], Quantity]
```
Returns the parameters and a `Quantity` describing their provenance. After this lands:
* every `models/*` dataclass default gets the `# PROVISIONAL` comment of §2.5;
* `dynamic_range`, `decomposition` and `oed` take a required `params` argument — no
  `model or KineticAllosteryModel()` fallbacks (`dynamic_range.py:51`,
  `decomposition.py:117`);
* and the `shapley_pct_*` fields in `decomposition.py` are renamed
  `log_share_pct_*` until someone implements an actual Shapley value. The present
  computation is `|Δ_i| / Σ|Δ_j|` (`decomposition.py:85-89`), which is a normalised
  absolute-log share: it is not permutation-averaged, it discards sign, and it does not
  sum to the total log gain. Three of the four terms also do not add up to `log_gain`
  by construction — the identity in the technical spec,
  `ln G = Δln(Occ) + Δln(Kinet) + Δln(Headroom)`, is **not** satisfied by the current
  definitions, because `headroom_factor = (Po_max - Po_base)/(Po_max - Po_pam)` is not a
  factor of `Po_pam/Po_base`. Either derive a decomposition that is exact and prove it in
  a test, or rename the fields to say they are a heuristic attribution. **Do not ship an
  identity in a docstring that the code does not satisfy.**

**Acceptance.**
* **ANCHOR** `test_decomposition_identity_is_exact`: `log_gain` equals the sum of the
  three terms to `rtol=1e-10`, for 20 parameter draws and both regimes. If you cannot make
  this pass, the decomposition is wrong and P4-3's renaming is mandatory.
* **ANCHOR** `test_mle_reproduces_training_peak_ec50`: the fitted Model B reproduces the
  P1-2 concentration-response EC₅₀ within its digitisation error.
* **ANCHOR** `test_posterior_median_in_profile_ci`: for each identifiable parameter, the
  posterior median lies inside the P3 profile-likelihood 95% CI. Two independent
  uncertainty machineries agreeing is the strongest check available here; disagreement
  means one of them is wrong.
* `test_mcmc_diagnostics_gate_the_result` — a deliberately under-run chain (10 steps)
  returns `VOID`, not an interval.
* `test_no_model_instantiated_without_registry_params` — AST-walk `src/` for
  `KineticAllosteryModel(` / `ExtendedDesensitizationModel(` / `OperationalScalarModel(`
  with no arguments, outside `tests/` and `parameters.py`.

**Done when** `parameters.get("kinetic_jw95")` returns a fitted set, the registry is the
only source of production parameters, and every credible interval shipped carries its
diagnostics.

**Forbidden.** Quoting a posterior from a chain that failed a diagnostic. Fitting to the
holdout. Keeping the 6-parameter fit if P3 says only 4 directions are identifiable —
sample the identifiable ones and *derive* the rest, stating the constraint.

---

### P5 — Model comparison that means something

**Goal.** Rank Models A, B, C by held-out predictive performance, with AIC/BIC as
secondary evidence, each at its own MLE.

**Why here.** C5. Requires P1 (observables + holdout), P3 (what k is) and P4 (the MLEs).

**Deliverables.** Rewrite `protocols/oed.py::evaluate_model_fit` as:

```python
def compare_models(models: Sequence[type[ReceptorModel]], train: Sequence[Dataset],
                   holdout: Sequence[Dataset], *, n_folds: int = 5) -> ResultSet
```
Protocol, in order:
1. Fit each model independently on `train` via P4-1/P4-2. Record each MLE and its
   convergence evidence.
2. `k` = number of **estimated identifiable** parameters (P3) **plus** estimated noise
   scales — not `len(param_names)`. Model C's `pam_desens_factor` counts only if the
   likelihood actually constrains it; prove that with a profile, do not assume it.
3. AIC, BIC and AICc (n is small — AICc is the honest choice and the current code omits
   it).
4. K-fold cross-validated held-out log predictive density, by *concentration block* for
   CRC data and by *time segment* for traces. Random point-wise folds leak information
   across a smooth curve and will flatter the flexible model.
5. The prospective holdout (P1-2 item 3), scored once, reported separately, never used to
   choose.
6. Akaike weights only alongside the CV ranking, and only if the two agree. If they
   disagree, report both and say the comparison is unresolved — that is a finding.

**Acceptance.**
* **ANCHOR** `test_nested_model_recovers_the_truth`: generate data from Model B at known
  parameters with realistic noise; the protocol must select Model B over Model A, and must
  *not* select Model C (which nests B) once the complexity penalty applies. If your
  protocol cannot recover a known answer on synthetic data, it cannot adjudicate a real
  one. This is the one legitimate use of synthetic data in this programme, and it must be
  labelled as a method check, never as evidence about receptors.
* **ANCHOR** `test_k_counts_identifiable_params_only`.
* `test_cv_folds_are_blocked_not_pointwise`.
* `test_sigma_is_estimated_not_assumed` — the ranking must be invariant to the removal of
  the old `measurement_noise_std` default.
* `test_comparison_requires_same_observable` (P1 contract).

**Done when** the ResultSet names a winner with a CV margin and a stated uncertainty, or
states that the data do not discriminate. Both are acceptable results. Neither may be a
number computed at unfitted parameters.

**Forbidden.** Reporting Akaike weights alone. Comparing a `PEAK`-native model to an
`EQUILIBRIUM`-native one on one axis. Treating a 2-point AIC difference as decisive.

---

### P6 — Waveforms, charge, and the discriminating experiment

**Goal.** Drive each fitted model with realistic waveforms, report charge transfer with
intervals, and compute the protocol that maximally separates the models under *both*
measurement noise and parameter uncertainty.

**Why here.** The design needs P4's posteriors to make a falsifiable interval, and P5's
fitted models to be worth separating.

#### P6-1 Waveforms (extends `protocols/waveforms.py`, post-P0-2)
* Biexponential synaptic transient with the `τ_rise → τ_clear` limit handled (P0-2).
* **Biexponential clearance**, which the technical spec requires (`τ_fast ≈ 1 ms`,
  `τ_slow ≈ 10–30 ms`) and the current code does not implement — `synaptic_transient` has
  a single `clear_ms`. Add `clear_fast_ms`, `clear_slow_ms`, `weight_fast`, and keep the
  single-tau call as a documented special case so existing results are reproducible
  (§2.6).
* Pulse trains at 10/50/100 Hz, ≥ 20 pulses, with paired-pulse ratio `I₂/I₁` and a
  steady-state-to-first-pulse ratio.
* Ambient + spillover, with the spillover amplitude and τ carrying `Basis` records.
* `extract_electrophys_metrics` gains the P0-4 `DecayFit` and declares `Observable.CHARGE`
  with its window.

#### P6-2 Three-tiered dynamic range, with intervals
Rewrite `dynamic_range.evaluate_dynamic_range` to take `n_posterior_draws: int = 1000`
and return each tier as a `Quantity` whose value is the posterior median with a
`[2.5, 97.5]` percentile interval in its provenance:
1. **Theoretical asymptotic headroom** — `P_open,∞ / P_open(ambient)`. State in the
   caveat that `k_off → 0⁺` is reachable by no ligand, and give the *formula*, not just
   the number (C7).
2. **Pharmacologically reachable gain** at `s_max ∈ {1.25, 1.50, 2.50}`. Note in the
   provenance that `s_max` here is an **equilibrium EC₅₀ shift** and that the legacy
   `gabaa_kinetics.calibrate_pam` calibrates against the **peak** EC₅₀ shift, so the same
   nominal `s_max` is a different `k_off` divisor in the two layers (2.5 vs 2.945 at the
   manuscript parameters). Both are internally consistent; numbers must not be moved
   between them without restating the observable.
3. **Physiological charge ratio** `ΔQ_tonic / ΔQ_phasic`, with the integration window in
   the provenance, plus a boundary scan over ambient GABA and `s_max` locating where the
   ratio reaches 1.0 (the collapse surface). Report the collapse boundary as a *curve*,
   not a boolean — `headroom_collapses: bool` throws away the information the scan
   produced.

#### P6-3 The discriminating protocol
Replace `find_discriminating_protocol`. The score must be a **posterior-predictive
separation**, not a difference of point predictions over a fixed noise:

```
score(protocol) = E_posterior | y_A - y_B |  /  sqrt( var_posterior(y_A)
                                                    + var_posterior(y_B)
                                                    + sigma_meas^2 )
```
Search over `[G]`, `[PAM]` (expressed as `s_max`, with its observable stated), pulse
duration, pulse frequency and the measured observable itself — the choice of *what to
measure* is part of the design and the current grid omits it. Use the posterior draws from
P4; 2–3 dimensions can be gridded, more needs a coarse-to-fine search.

Falsification boundaries come from the posterior predictive, two-sided, with the
multiple-comparison cost of having searched the protocol space stated. A pre-registration
block goes in `knowledge/12-preregistration.md`: the protocol, the predicted intervals per
model, the decision rule, and what outcome would refute the kinetic hypothesis **before**
any such experiment is run.

**Acceptance.**
* **ANCHOR** `test_charge_ratio_recovers_steady_state_limit`: for a long square
  application the `CHARGE` gain converges to the `EQUILIBRIUM` gain within 1% as the
  window grows past `5/r`. Ties the new machinery to the analytic result.
* **ANCHOR** `test_paired_pulse_depression_sign`: at 100 Hz the second pulse peak is lower
  than the first for the fitted model, and the depression deepens with frequency —
  monotonically across 10/50/100 Hz. A desensitisation scheme that does not do this is
  mis-wired.
* **ANCHOR** `test_discrimination_score_falls_when_posterior_widens`: inflate the
  posterior spread 10× and the best score must drop. The current implementation is
  insensitive to parameter uncertainty by construction, so this test is the regression
  guard for C6.
* `test_biexponential_clearance_reduces_to_single_tau` (§2.6 reproducibility).
* `test_collapse_boundary_is_monotone_in_ambient`.

**Done when** each tier reports median + CI + its observable + its window, and the
pre-registration document exists with numbers in it.

**Forbidden.** `± 2σ_noise` as a falsification interval. A boolean collapse flag with no
boundary. Quoting a single headline fold-change without its interval and its tier.

---

### P7 — Manuscript integration, and the circuit boundary

**Goal.** Fold P1–P6 into `knowledge/08-manuscript.md` with the consistency tests green,
and state explicitly what this programme does *not* claim.

**Deliverables.**
* Extend `tests/test_manuscript_consistency.py` to the new quantities, using the same
  blunt string-grep technique and the same `main_text` / Supplementary split. Every number
  in the manuscript must be recomputable by a test.
* Quote intervals, not point estimates, for everything that depends on a fitted parameter.
  Where a point value is quoted (an algebraic identity), say so in the sentence.
* A scope box, in the pattern the current manuscript already uses, stating: the fitted
  parameters are from N recombinant datasets in a stated preparation; the circuit-level
  consequences remain conditional on the `UNSOURCED` subunit fractions
  (`provenance.REGIONS_PROV` — 19 of 28 parameters are `UNSOURCED`); nothing here is
  evidence that α5-selective compounds preserve respiration.
* Update `knowledge/10-model-technical-spec.md` to match what was built. It currently
  describes several things that do not exist (biexponential clearance in `waveforms.py`,
  single-channel datasets, state-dependent PAM binding "via a thermodynamic cycle" in
  Model C — the implementation scales `d_fast` by a linear factor, which is not a
  thermodynamic cycle and does not satisfy detailed balance around the R/O/D loop).
  **A spec that overstates the implementation is a defect of the same kind as a
  mis-cited number.**

**The circuit boundary — OUT OF SCOPE for P1–P7, and here is the contract for later.**
Coupling `CHARGE` into `resp.py` requires: a conductance scale (nS per unit charge) with
its own anchor; an operating point for the target substrate (`substrate.resolve_op`, and
see P0-9); a validity band appropriate to the preparation (`config.INVITRO_BAND` vs
`EUPNOEA_BAND` — the two are a preparation difference, not a tolerance); and settling at
≥ 3·τ_h = 30 s on the conductance substrate. Anyone who wires a receptor-level charge
into a circuit without all four is repeating the error documented at length in
`config.py:59-94`, where a confident, seed-verified operating point turned out to be a
reading of a decaying transient. **Do not start this without a new roadmap revision.**

---

## 5. Forbidden patterns

Each is a real defect from this repository. The file references are so you can read the
original.

**5.1 The silent fallback.** Returning a default, a target, or the input when a
computation fails. `kinetic_jw95.py:199`, `operational.py:109`,
`identifiability.py:84`. See §2.4.

**5.2 The sentinel that is also a plausible answer.** `decay_tau = 15.0` where 15.0 is
the fit target (`kinetic_jw95.py:218`). Use NaN.

**5.3 The second copy of a constant.** `SYNAPTIC_PULSE` existed in three places with two
different values for weeks; the disagreement surfaced only when a script compared a
`Drug`'s tonic gain against the manuscript's, and an audit "refuted" the manuscript on the
strength of the wrong copy (see the note now in `gabaa_kinetics.py:99-121`). Import from
`config.py`. Always. `resp.py:37` still violates this.

**5.4 The fabricated dataset.** Generating `y` from a formula and attaching a citation.
`fitting/data.py:43-56`. If you cannot get the real data, say so and stop.

**5.5 The neighbouring-quantity citation.** Attaching a source that resolves perfectly by
DOI and measures something adjacent. Twice in this project: `a5_dist` (regional *level*,
cited for regional *composition*) and Kasugai 2010 (measured α1/α2/β3, cited for α5).
Read the paper against the specific number.

**5.6 The cross-observable comparison.** `subjective_index` vs `regional_sens` (1.8×
apart, both sound like "forebrain engagement"); whole-body ventilation vs isolated preBötC
(now `VOID`); peak-current data fitted with an equilibrium curve
(`identifiability.py:37`). See §2.3.

**5.7 The mechanics test that cannot fail.** `assert cost >= 0.0` on a χ² of 9290
(`test_nextgen_models.py:241`). See §2.7.

**5.8 The unanchored parameter set that produces a headline.** Three uncalibrated model
defaults feeding `evaluate_dynamic_range()`. See §2.5.

**5.9 The identity asserted in a docstring and not satisfied by the code.**
`ln G = Δln(Occ) + Δln(Kinet) + Δln(Headroom)` (technical spec §3,
`decomposition.py:1-7`). Prove it in a test or rename the fields.

**5.10 The precision that the inputs cannot carry.** Seven significant figures from a
fit with more parameters than residuals (C7). Round to what the inputs support and say
why — `test_manuscript_consistency.py` already enforces both halves of this for Table 6.

**5.11 Measuring inside a transient.** A 4 s warm-up on a cell whose slowest timescale is
10 s, which produced a confident wrong operating point verified across four seeds
(`config.py:70-85`). On any new dynamic measurement, state the slowest timescale and warm
up for ≥ 3× it.

**5.12 The partial operating point.** Supplying some conductance-substrate parameters and
letting the rest stay LIF-scaled. `circuit.py:83` (P0-9). The full set is required, not
defaulted.

---

## 6. What may be claimed, and when

| Claim | Permitted after | Tier |
|---|---|---|
| "The 5-state scheme's equilibrium EC₅₀ is `Kd(1+√(2+E))/(1+E)`" | now (algebra) | VALIDATED |
| "Phasic and tonic pools sit in different agonist regimes" | now (definitional) | VALIDATED |
| "Tonic headroom exceeds phasic headroom by ~X×" | P4 | UNCALIBRATED + CI |
| "The divergence holds across the data-consistent ensemble" | P4 (posterior) | UNCALIBRATED + CI |
| "Occupancy/gating/headroom contribute X/Y/Z%" | P4-3 (exact identity) | UNCALIBRATED |
| "Model B is preferred over A and C" | P5 (CV + holdout) | UNCALIBRATED |
| "The scheme predicts a held-out observable it was not fitted to" | P5 holdout | VALIDATED |
| "Protocol P discriminates the models; outcome O refutes kinetic allostery" | P6 (pre-registered) | UNCALIBRATED |
| "This bears on respiratory safety of α5-selective PAMs" | **never, on this programme** | — |

The last row is not pessimism. The circuit-level argument rests on
`REGIONS["prebotc"]["a5"] = 0.02` and `EXTRASYN["a5"] = 0.80`, both `UNSOURCED`
(`provenance.py:82`, `:134`), and the respiratory axis has no valid quantitative anchor at
all. Receptor-kinetic results are untouched by that — which is exactly why the manuscript's
scope box separates them, and why you must not quietly re-merge them.

---

## 7. Agent workflow protocol

**Branch.** Work on the branch you were assigned. If none, `claude/<topic>-<short-id>`.
Never push to `main`.

**One concern per commit.** A P0 item, or one numbered task. Commit message: what changed,
why, and the measured evidence. The repository's own history is the style guide — read
`git log --oneline -20`; messages name the defect and its consequence, not the files
touched.

**Before every commit:**
```bash
pytest -q -n 0 tests/test_<the_file_you_changed>.py   # fast, stdout visible
pytest -q                                             # the whole suite, ~15 min
```
A commit that leaves a test red is acceptable **only** if the message says which test,
why, and which task will fix it. Silence is not.

**Report format.** When you finish a task, state: what you changed; the test names that
now pass and their measured values; what you did **not** do and why; anything you found
that is outside your task (do not fix it — record it here, and it becomes a P0 entry).

**When you are blocked, stop and report.** Specifically, stop if: a required dataset is
not obtainable (do not synthesise — §5.4); a test demands a number the model cannot
produce (do not loosen the tolerance — report the discrepancy); the fit will not converge
(do not take the best of 32 bad optima — report the spread); or an instruction here
conflicts with the code (the code is evidence, this document may be wrong — say so and
quote both).

**Never** relax a validity band, a tolerance, or a threshold to make a result qualify. The
one time this project moved a band (`INVITRO_BAND`) it was legitimate, and the note
justifying it runs 30 lines and cites the measurements and the source that were collected
*before* the constant was introduced (`config.py:183-212`). That is the standard.

---

## 8. Phase status — update this table in the same commit that closes a phase

| Phase | Status | Closing commit | Notes |
|---|---|---|---|
| P0 Defect register | **OPEN** | — | 13 items; blocks everything |
| P1 Observable contract + data | BLOCKED on P0 | — | may start with P2 |
| P2 Tier/provenance integration | BLOCKED on P0 | — | parallel with P1 |
| P3 Identifiability | BLOCKED on P1 | — | |
| P4 Fit + posterior + registry | BLOCKED on P3 | — | |
| P5 Model comparison | BLOCKED on P4 | — | |
| P6 Waveforms + OED | BLOCKED on P4 | — | |
| P7 Manuscript + spec | BLOCKED on P6 | — | |

---

## 9. Architectural guardrail (unchanged from Revision 1)

> **Never add biological complexity solely for sophistication.**
> Do not introduce more receptor states, unanchored conductance parameters, or larger
> circuit models unless:
> 1. it explains empirical data that the simpler model structurally fails to explain, **or**
> 2. it generates a distinct, testable, falsifiable prediction.

Revision 2 adds a second guardrail, because the 2026-10-08 sprint showed the first is not
sufficient on its own:

> **Never add analytical sophistication on top of an unanchored parameter set.**
> An MCMC posterior, a Shapley attribution, an AIC ranking and an optimal experimental
> design are all machinery for turning data into conclusions. Run on uncalibrated
> parameters, a fabricated dataset, or a flat likelihood, each one produces a confident,
> well-formatted, publishable-looking number that means nothing — and is far harder to
> detect as wrong than the simple error underneath it.
