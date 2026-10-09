# The likelihood, the MLE and the posterior

**Roadmap phase:** P4. **Status:** machinery built and method-validated; one structural
result about the scheme that does not depend on any fit.

Read `knowledge/11-identifiability.md` first. It establishes what the data can determine;
this document is what was built on top of that and what it says.

---

## 1. What the likelihood does differently

`fitting/likelihood.py`, replacing the chi-squared in `identifiability.py`:

| | before | now |
|---|---|---|
| observable | every dataset scored against the equilibrium curve | each dataset scored on its own `Observable` |
| scale | normalised peak data compared to absolute predictions | normalised data compared to normalised predictions |
| σ | hand-set | profiled out analytically, one per dataset |
| constants | dropped | kept, because P5 compares across models |
| k for an information criterion | parameter count | parameters **plus** each estimated σ |
| holdout | remembered | `holdout_guard` raises |

### 1.1 σ profiled out is σ marginalised out, here

Substituting the MLE `σ² = RSS/n` gives `const − (n/2) log RSS`. Profiling a nuisance
parameter is **not** in general the same as marginalising it, so the posterior in §3 would
be using the wrong object if that were all. It happens to be the same here: marginalising σ
under the scale-invariant prior `1/σ` gives `Γ(n/2) 2^(n/2−1) RSS^(−n/2)`, whose logarithm
is `const − (n/2) log RSS` — the same parameter dependence, and the constant it drops
depends only on `n`, which is fixed. So §3 is the marginal posterior over the three ratios
with σ integrated out under a Jeffreys prior. Stated explicitly because assuming it would
have been a real defect.

### 1.2 The second observable mismatch: normalisation

P1 caught PEAK-against-EQUILIBRIUM. It did not catch **PEAK against normalised PEAK**.

A published concentration-response reports `I/I_max`, with an asymptote of 1 by
construction. This scheme's absolute peak open probability saturates near 0.75 at the
three-anchor parameters. Compared directly, the residual at saturation is ~0.25 and the
only way to close it is to raise the plateau — which the optimiser did, by **deleting
desensitisation**: `D → 1e-11`. Same failure mode as the first mismatch, a different cause.

`fitting.data.Normalisation` now makes a dataset declare which it is, and
`likelihood._predict` divides the prediction at the same concentration the experimenters
normalised by (the largest in the dataset). The consequence is intended: a normalised
dataset then constrains only the **shape** of the curve, which is the only thing it carries
information about. Absolute open probability has to come from somewhere else.

---

## 2. The structural conflict: the slope and the plateau

**This is the main result of P4, and it is not about a fit.**

Jahn 1997 reports a peak Hill slope of **nH = 2.2 ± 0.4** for α1β2γ2L. This project anchors
the absolute peak open probability at **P_o,max = 0.750**. Scanning `E ∈ [0.1, 1000]` and
`D ∈ [1e-3, 1000]` — four and six decades, 33 × 31 points, at fixed `K_d`:

| constraint imposed | best the scheme can do on the other |
|---|---|
| absolute peak `P_o,max ∈ [0.73, 0.77]` | **nH = 1.33** (at `E = 3.16`, `D = 1e-3`) |
| `nH ≥ 1.8` | **`P_o,max ≥ 0.99`** (lowest found 0.9898, at `E = 100`, `D = 3.98`) |
| neither | nH up to 2.46 is reachable, but only with `P_o,max` near 1 |

**The two cannot both hold.** The measured slope is more than two of its own standard
deviations away from the steepest curve the scheme can produce at the assumed plateau. A
single open state reached only from the doubly-liganded closed state cannot be both that
steep and that leaky.

Two things this is not:

* **Not an artefact of the equilibrium algebra.** The equilibrium curve's shape depends
  only on `F = 1 + E + D` while its plateau is `E/F`, so there those two *are* independent
  and no conflict exists. The conflict is specific to the PEAK observable, where
  desensitisation competes during the application.
* **Not settled in favour of either side.** Note which is sourced: the slope is published;
  the 0.750 is one of the project's own fit **targets** (`parameters._NOMINAL_DEFECT`,
  roadmap P0-13), carried for years as if it were a measurement. The honest statement is
  *the scheme cannot satisfy both*, and the next step is to find out which is wrong —
  which is a P6 preregistration question, not something to resolve by fitting harder.

Pinned in `test_the_slope_and_the_plateau_cannot_both_hold`.

### 2.1 What the fit to the Jahn curve gives

Fitted to the normalised Jahn PEAK curve with `k_off`, `α` and `r` held at declared
conventions (`fit_mle`, 8 starts, 11 of 12 converging and the best five agreeing to
1.2e-6 decades):

| | fitted to Jahn shape | three-anchor fit | Jahn 1997 | project anchor |
|---|---|---|---|---|
| `K_d` (µM) | 104.4 | 23.81 | — | — |
| `E` | 163.9 | 3.925 | — | — |
| `D` | **0.001, on its bound** | 25.0 | — | — |
| PEAK EC₅₀ (µM) | 11.663 | 23.16 | 11.6 ± 0.9 | 20 |
| PEAK nH | 2.030 | 1.314 | 2.2 ± 0.4 | — |
| absolute peak `P_o,max` | 0.9939 | 0.7483 | — | 0.750 |

The fit reproduces both published numbers to well inside their errors — and pays for it
with `D` on its lower bound and a plateau of 0.99. `MLEResult.at_bound` reports the pinning
and `note` says in words that the value is the edge of the search space and not an estimate.
That is the same conflict as §2, arrived at from the other direction.

---

## 3. The posterior, and why it needed a sampler fix

`fitting/posterior.py` samples `(log₁₀ K_d, log₁₀ E, log₁₀ D)` — the combinations §1 of
`11-identifiability.md` proves an equilibrium observable determines — against the
per-observable likelihood, with a flat prior in log₁₀ (log-uniform in the ratio).

Walkers start **from the prior**, not from the MLE. A chain started in a small ball at the
optimum reports a split-R̂ near 1 whether or not it mixed, because the statistic compares
between-walker to within-walker variance and the walkers began identical.

### 3.1 The frozen-walker trap

Over-dispersed starts exposed a real defect in the sampler. The Goodman & Weare stretch
move proposes `y = x_j + z(x_k − x_j)` with `z ∈ [1/a, a]`, so at the canonical `a = 2` it
**cannot propose a contraction tighter than one half**. Measured on this posterior:

* 2 of 24 walkers accepted **94 and 405** moves against a median of **3350**;
* they froze 2.5 decades from the mode while the other 22 sat on the right answer;
* split-R̂ came out **62**, and the reported 97.5th percentile was a stuck walker's
  position rather than a posterior tail;
* along the line from one of them to the mode the density fell by **16 log units** before
  rising, so the move it needed was `z = 0.14` and the move it could make was `z = 0.5`.

| `a` | acceptance | τ | worst split-R̂ | frozen |
|---|---|---|---|---|
| 2 | 0.49 | 57 | 62 | 2 of 24 |
| 50 | 0.10 | 78 | 1.026 | 0 |
| 200 | 0.045 | 175 | 1.034 | 0 |
| **(2, 20, 200) mixed** | **0.28** | **59** | **1.006** | **0** |

`a` may now be several scales, sampled uniformly per proposal: each fixed-`a` stretch move
satisfies detailed balance with respect to the target, and choosing between them
independently of the current state leaves that intact, so the mixture is exact rather than
a tuning trick. It beats every single scale on every diagnostic.

The freezing is now **detected**: per-walker acceptance is reported, and a walker that
moved less than a tenth as often as the median is a diagnostic failure, because an ensemble
with a frozen member has not sampled anything. `test_a_frozen_walker_is_a_diagnostic_failure`
reproduces it on purpose; `test_mixing_stretch_scales_is_what_fixes_it` pins the fix.

### 3.2 Two machineries agreeing

On a method-check dataset (30 concentrations, 0.1% noise, `sem` equal to the noise that was
added), truth `K_d = 25 µM, E = 4, D = 25`:

| | posterior median | posterior 95% | profile 95% | truth | covered |
|---|---|---|---|---|---|
| log₁₀ K_d | 1.3703 | [1.3001, 1.4448] | [1.3087, 1.4455] | 1.3979 | yes |
| log₁₀ E | 0.5420 | [0.4259, 0.6681] | [0.4008, 0.6961] | 0.6021 | yes |
| log₁₀ D | 1.3347 | [1.2115, 1.4664] | [1.1858, 1.4953] | 1.3979 | yes |

Profile likelihood re-optimises the other parameters at each node of a grid; MCMC
integrates over them. They share no optimiser, no grid and no convergence criterion, so
agreement is the strongest check available here — and they agree on the intervals, not just
the point. `agreement_with_profiles` computes it and
`test_posterior_median_in_profile_ci` asserts on it.

**`sem` has to be the noise.** At a declared `sem` of 1e-3 against 2e-4 of actual noise the
profile intervals came out **five times too wide** and the two machineries disagreed for a
reason that was about neither of them. `equilibrium_crc_from_model` no longer floors it.

### 3.3 The gate

A chain that fails τ, split-R̂, acceptance or the frozen-walker check returns intervals that
are **VOID** — reading one raises — rather than a number with a footnote. A deliberately
under-run chain (10 steps) is tested for exactly this. An unconverged chain's percentiles
are not an interval; they are where the walkers happened to be.

Tier, in order of precedence: VOID if a diagnostic failed; UNCALIBRATED if any dataset is
parametric or synthetic; VALIDATED only if the chain converged against measured data. No
posterior in this repository is VALIDATED, because no dataset here is a digitisation.

---

## 4. The parameter registry

`circuitpharm/parameters.py` is the single source. `get(model, which)` returns
`(rates, Quantity)`, so there is no way to take the numbers without being handed their
provenance — and for `which="nominal"` the quantity is VOID, so reading it raises.

| model | fitted | declared | nominal |
|---|---|---|---|
| `kinetic_jw95` | UNCALIBRATED (3 anchors, α a convention) | — | VOID |
| `operational` | — | UNCALIBRATED (two values are the project's own targets) | VOID |
| `extended_desens` | **raises** | — | VOID |

Model C refuses a `fitted` set because it has never been fitted and its defaults have no
recorded origin. Returning them with a tier on them would let P5 rank a fitted model
against an unfitted one and call the difference evidence about mechanism.

`test_no_model_is_instantiated_from_its_defaults_in_production_code` AST-walks `src/` and
fails on any bare `KineticAllosteryModel()`, which is how the defect it guards (P0-13 — the
project's headline numbers produced at a chimeric parameter set that reproduced no anchor)
gets to stay fixed.

---

## 5. A performance change with a correctness test attached

`peak_dose_response` solved an ODE per concentration. At constant agonist the generator is
constant, so `P(t) = P(0) exp(Qt)` and there is nothing to integrate. Replaced with one
matrix exponential plus repeated doubling:

* **agreement with the `solve_ivp` path: 1.8e-7 absolute** over nine concentrations
  spanning four decades — the difference is the integrator's own tolerance, not this
  method's error;
* **~300× faster**: 0.30 s → 1.0 ms for nine concentrations;
* `scipy.linalg.expm` turned out to cost 8 ms on a 5×5, essentially all Python-level
  overhead, so the exponential is done by scaling and squaring in a dozen lines — checked
  against `scipy.linalg.expm` to 2.6e-10 on the generators this project builds.

This is why a PEAK likelihood evaluation costs under 2 ms instead of 240 ms, which is the
difference between the Jahn MLE being a 20-minute job and a 50-second one, and between P5's
per-model fits being feasible and not. It also cut the full test suite from 13m36s to
6m24s. Both paths are kept and `test_exact_peak_agrees_with_the_ode_it_replaced` compares
them, because a 300× speedup with no test is how every PEAK number in the project moves for
a reason that is not physics.

---

## 6. Reproducing this

```python
import numpy as np
from circuitpharm.models import KineticAllosteryModel
from circuitpharm.fitting.data import JAHN1997_PEAK_CRC, equilibrium_crc_from_model
from circuitpharm.fitting.likelihood import fit_mle
from circuitpharm.fitting.posterior import (
    agreement_with_profiles, sample_identifiable_posterior)

TRUTH = KineticAllosteryModel(kon=0.01, koff=0.25, beta=0.8, alpha=0.2, d=0.05, r=0.002)
conv = dict(koff=TRUTH.koff, alpha=TRUTH.alpha, r=TRUTH.r)
factory = lambda p: KineticAllosteryModel(**p.to_microscopic(**conv))

ds = equilibrium_crc_from_model(TRUTH, concs_um=np.logspace(-1, 3.5, 30), noise_sd=1e-3)
post = sample_identifiable_posterior(factory, [ds])     # §3.2, ~20 s
print(post.note)
print(agreement_with_profiles(post, ds))

print(fit_mle(factory, [JAHN1997_PEAK_CRC], n_starts=8).note)   # §2.1, ~35 s
```

All of it is pinned in `tests/test_p4_likelihood.py`. The sampling tests are marked `slow`.
