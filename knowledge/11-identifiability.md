# Identifiability of the 5-state GABA-A scheme

**Roadmap phase:** P3. **Status:** structural result proved and verified; practical result
measured on a method-check dataset, because no real EQUILIBRIUM dataset exists yet.

This answers the question the roadmap insists comes before any fitting: *what can the data
determine?* Doing it the other way round is how the 2026-10-08 sprint produced profile
likelihoods, a Fisher information matrix and an MCMC posterior over six rates that an
objective could only constrain three combinations of.

---

## 1. The structural result, derived

The equilibrium state distribution of the scheme

```
R  <--2k_on[G]/k_off-->  AR  <--k_on[G]/2k_off-->  A2R  <--beta/alpha-->  A2O
                                                    |  d / r
                                                    v
                                                   A2D
```

is, writing `x = [G]/K_d`,

```
P_R : P_AR : P_A2R : P_A2O : P_A2D   =   1 : 2x : x² : E x² : D x²
```

with

```
K_d = k_off / k_on        E = beta / alpha        D = d / r
```

**Every equilibrium observable is therefore a function of exactly those three quantities.**
The six microscopic rates enter only through the three ratios, so the equilibrium
likelihood is invariant under a three-parameter scaling group: multiply any of the pairs
`(k_on, k_off)`, `(beta, alpha)` or `(d, r)` by a common factor and nothing changes.

Numerical corroboration (`fitting.identifiability.invariance_report`): the objective is
constant to **2e-12** along each generator, with normalised curvature `vᵀHv / max|H|` of
**+3.0e-12, −1.0e-12 and +5.3e-11**. Scaling all six rates at once — the sum of all three
generators — changes the cost by exactly 0 at six decimal places.

### 1.1 What this is *not*

An earlier revision of the roadmap claimed the cost Hessian is "structurally rank ≤ 3".
**That does not follow and is not what is measured.** The generators are a vector *field* —
each depends on the current parameters — so zero curvature along them does not imply a zero
eigenvalue, and a finite-difference Hessian here carries ~1e-5 relative error anyway.
Thresholding its eigenvalues reports rank 6, with three of them negative (the parameters
are not at a minimum, so the matrix is not a Fisher information matrix either). The group
invariance is the provable statement and is the one the code tests.

### 1.2 The ladder: what each observable unlocks

A **kinetic** observable does not share the invariance, because scaling a pair changes how
fast the system equilibrates while leaving where it equilibrates alone.

| likelihood contains | determines | of 6 |
|---|---|---|
| equilibrium CRC only | `K_d`, `E`, `D` | 3 |
| + a deactivation time course | + `k_off`, hence `k_on` | 4 combinations |
| + a mean open time (1/alpha) | + `alpha`, hence `beta` | 5 |
| + a desensitisation time course | + `d`, hence `r` | 6 |

`fitting.reparam.UNLOCKED_BY` holds this table in code. **This repository has none of those
three kinetic datasets** (`fitting.data.MISSING_DATASETS`), which is why
`gabaa_kinetics.FIT_FIXED_ALPHA` is a `CONVENTION` and why acquiring a single-channel mean
open time is the highest-value data acquisition in the project.

### 1.3 Scope

This is about the **equilibrium** likelihood in `fitting/identifiability.py`.
`gabaa_kinetics.fit_scheme` fits PEAK and IPSC-decay observables, which *do* carry
timescales, so it is not degenerate in this way — it was degenerate for the simpler reason
of having four unknowns against three residuals (roadmap P0-7, fixed by holding `alpha`).

---

## 2. The practical result, measured

Structural identifiability is necessary and nowhere near sufficient. The method was
verified against a dataset generated from the model at known parameters
(`fitting.data.equilibrium_crc_from_model` — the one legitimate use of generated data, a
method check), truth `K_d = 25 µM, E = 4, D = 25`:

| | log₁₀ K_d | log₁₀ E | log₁₀ D |
|---|---|---|---|
| **truth** | 1.3979 | 0.6021 | 1.3979 |
| **MLE, noiseless** | 1.3979 | 0.6021 | 1.3979 |
| **MLE, 0.2% noise** | 1.7339 | 1.2168 | 2.0336 |

**Noiseless recovery is exact** (cost 5.6e-18), so the estimator is correct and all three
combinations are structurally identifiable, as §1 says.

**At 0.2% noise on a normalised curve, it falls apart:**

| | 95% CI width (decades), noiseless | at 0.2% noise | classification at 0.2% |
|---|---|---|---|
| log₁₀ K_d | 0.195 | 1.335 | IDENTIFIABLE |
| log₁₀ E | 0.474 | 2.115 | PRACTICALLY NON-IDENTIFIABLE |
| log₁₀ D | 0.491 | 2.137 | PRACTICALLY NON-IDENTIFIABLE |

`E` and `D` span **more than two decades** — a factor of ~130 — and the true values sit at
the very edge of each interval. The fitted set is 2–4× away from the truth in every
coordinate and fits the noise realisation *better* than the truth does (cost 4.96 against
11.41), while the two predicted curves differ by at most 2.4e-3, which is the noise.

### 2.1 The consequence for the programme

**More precision on an equilibrium concentration-response will not fix this.** Two
parameter sets a hundredfold apart produce curves separated by less than realistic
measurement noise, so the information is not there to be extracted. The fix is a
*different observable* — one with a timescale in it — exactly as §1.2 says. That is the
quantitative argument for roadmap P1-2 being the blocking task, and it is stronger than
the qualitative one the roadmap originally gave.

### 2.2 And a trap found while doing it

Fitting the **equilibrium** curve to a **PEAK** dataset drives `D` to its lower bound. The
equilibrium plateau is `E/(1 + E + D)`; a peak-current curve plateaus near 0.75, and the
only way an equilibrium curve reaches that is by deleting desensitisation. Measured:
`log10_D → −15`, i.e. `D → 0`. The fit "succeeds" having silently removed a mechanism.

`equilibrium_chi2_identifiable` now refuses a dataset whose `observable` is not
`EQUILIBRIUM`. Both concentration-response curves in `fitting/data.py` are `PEAK`, so the
equilibrium analysis has **no real dataset to run on** — which is itself the finding, and
is why §2 is a method check rather than a result about receptors.

### 2.3 And a second one, found in P4

The `Observable` tag was not enough. A published concentration-response reports `I/I_max`,
asymptote 1 by construction; this scheme's **absolute** peak open probability saturates near
0.75. Comparing them directly leaves a ~0.25 residual at saturation that can only be closed
by raising the plateau, and the optimiser did it the same way as above: `D → 1e-11`. So
PEAK-against-EQUILIBRIUM was caught here and **PEAK-against-normalised-PEAK was not**.
`fitting.data.Normalisation` now makes a dataset declare its scale and the likelihood
normalises the prediction to match. See `knowledge/12-inference.md` §1.2, and §2 there for
the apparent conflict this exposed between the published Hill slope and the assumed
plateau — **since WITHDRAWN**: the published slope is a local one over 1–10 µM and was
being compared against a whole-curve regression. §2 there also records a THIRD mismatch
class found while correcting it: the PEAK observable is protocol-dependent, moving its
EC₅₀ sevenfold with the application duration.

---

## 3. What is VOID, what is UNCALIBRATED, and why

| quantity | tier | reason |
|---|---|---|
| per-rate profile likelihood, FIM condition number, 6-rate MCMC intervals | **VOID** | computed over directions the objective cannot see (§1); the intervals are prior widths |
| `n_identifiable_of_six = 3` | **VALIDATED** | structural, proved in §1 and corroborated to 2e-12 |
| profiles over `K_d`, `E`, `D` | **UNCALIBRATED** | the method is sound and verified (§2), the dataset is not a digitisation |
| anything about a real receptor from §2 | — | not claimable at all; the dataset is generated |

The VOID machinery stays callable behind `provisional.ProvisionalResultWarning` rather than
being deleted — see roadmap §2.4 and forbidden pattern 5.4a.

---

## 4. Reproducing this

```python
from circuitpharm.models import KineticAllosteryModel
from circuitpharm.fitting.data import equilibrium_crc_from_model
from circuitpharm.fitting.identifiability import (
    fit_identifiable, invariance_report, profile_likelihood)

truth = KineticAllosteryModel(kon=0.01, koff=0.25, beta=0.8, alpha=0.2, d=0.05, r=0.002)
print(invariance_report(truth))                      # §1: the three flat directions

for sd in (0.0, 0.002):                              # §2: structural vs practical
    ds = equilibrium_crc_from_model(truth, noise_sd=sd)
    print(fit_identifiable(ds))
    for p in ("log10_kd", "log10_E", "log10_D"):
        print(profile_likelihood(p, ds))
```

Pinned in `tests/test_identifiability_p3.py`. The noiseless-recovery test is the method
validation; the 0.2%-noise test is the result.

One number here moved in P4: `equilibrium_crc_from_model` used to floor the declared `sem`
at 1e-3, so a dataset generated with less noise than that misstated its own uncertainty.
That widened the profile intervals by the same factor — five times, at `noise_sd = 2e-4` —
and made profile likelihood and MCMC disagree for a reason that was about neither. The
0.2%-noise numbers in §2 are unaffected (0.002 > 1e-3), and `knowledge/12-inference.md`
§3.2 has the corrected comparison.
