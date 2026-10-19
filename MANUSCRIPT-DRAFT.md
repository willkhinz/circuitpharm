# DRAFT — not submitted, not peer reviewed

**Authorship, affiliation and acknowledgements are deliberately left blank.** They are the
user's to assign; nothing here should be read as claiming authorship on anyone's behalf.

**Provenance of this draft.** Every number below was produced by a script committed in this
repository and was checked independently against that script before being written here. Where a
source was read only as an abstract, the text says so. Nothing is cited that was not read, and
no number appears that a committed script does not print. Three claims made earlier in this work
were retracted; they are described in §6 rather than omitted.

---

## Working title

Parametric reconstructions inflate leverage in multi-observable fits of GABA-A receptor
kinetics: a second conducting open state is necessary but not sufficient

## Abstract (draft)

Kinetic schemes for the GABA-A receptor are commonly fitted to published current waveforms that
are themselves reconstructions — curves regenerated from a small number of fitted parameters in
the source paper. We show that scoring such a reconstruction as if it were raw data inflates its
statistical leverage by roughly an order of magnitude relative to a co-fitted scalar
measurement, and that this inflation, not the data, produced six successive conclusions that two
α1β2γ2 datasets were mutually incompatible. Every one of those conclusions dissolved when the
objective was reformulated to score each dataset against the published parameters and errors
from which it was generated.

With a corrected objective we then ask what the deactivation waveform actually requires. In
five-state schemes with a single conducting open state, the fastest relaxation mode's fractional
amplitude is constrained to be non-positive, whereas the measured fast component carries 0.57 of
the amplitude; this imposes an exact and irreducible penalty of (0.57/0.04)² = 203.06 on any such
scheme. Adding a shut pre-open ("flip") state at matched parameter count does not relieve it. A
second conducting open state does.

A scheme with two conducting open states and positive binding cooperativity satisfies five
independent observables — deactivation components, concentration–response EC50 and Hill slope,
paired-pulse recovery, peak open probability, and mean open time — with positive residual
degrees of freedom. Profiling the cooperativity shows a broad plateau: the fit is insensitive to
the site-to-site ratio over c ∈ [35, 2500], so published estimates are not excluded. One
quantity resists: the Hill slope reaches at most 1.25 among admissible fits against a published
2.2 ± 0.4, a deficit of 2.5–2.9σ that no parameter choice in this model class removes. We
identify a single measurement — the 10–90% current risetime under rapid 10 µM GABA application —
that discriminates across the cooperativity plateau that the entire fitted dataset leaves flat.

## 1. The leverage problem

### 1.1 What the datasets are

Two α1β2γ2 datasets are in routine use here.

* **Deactivation** (Barberis et al. 2007, *Eur J Neurosci* 25(9):2726–2740; read in full text).
  The paper reports a three-exponential fit to the current decay following a 2 ms pulse of 10 mM
  GABA at 22–24 °C: τ = 2.8 ± 0.3, 33.4 ± 4.6, 221.35 ± 14.9 ms with fractional amplitudes
  0.57 ± 0.04, 0.23 ± 0.02, 0.20 ± 0.03, n = 6, and a weighted τ_w of 52.5 ± 2.9 ms. **Six
  numbers.** The "trace" used in fitting is that triexponential evaluated on a 101-point grid.
* **Concentration–response** (Jahn et al. 1997, *NeuroReport* 8(16):3443–6). EC50 = 11.6 ± 0.9
  µM and n_H = 2.2 ± 0.4. **Two numbers.** The "curve" is a nine-point Hill function generated
  from them.
* **Paired-pulse recovery** (Barberis et al. 2007). A single digitised scalar, 0.33 ± 0.03, n = 8.
  This one is a measurement, not a reconstruction.

### 1.2 The inflation, measured

Both reconstructions were scored with the residual variance profiled out, so the likelihood
rewards curve-matching without bound: there is no measurement noise in a curve generated from
parameters. Leverage per dataset, as n/σ², printed by `scripts/pseudodata_leverage.py`:

| weighting of the deactivation reconstruction | deactivation | recovery scalar | imbalance |
|---|---|---|---|
| σ profiled out (as used) | 4764 | 1111 | 4.3× |
| σ fixed at the single-dataset baseline | 11170 | 1111 | **10.1×** |
| each dataset divided by its own n | 110.6 | 1111 | **0.10×** |

The first two over-weight the reconstruction; the third over-weights the scalar by the same
order. No principle in the data selects among them. **The relative weight of two measurements
was therefore being set by a convention nobody had declared.**

### 1.3 Six verdicts that did not survive the correction

Scoring each parametric dataset against the published parameters and errors it was generated
from removes the convention entirely: the relative leverage is then fixed by the errors the
original authors reported. Applied successively, this dissolved every incompatibility
conclusion reached in this work:

| conclusion as first reached | after the objective was corrected |
|---|---|
| deactivation vs recovery irreconcilable | artefact of profiling σ on a reconstruction |
| conflict survives freeing rate conventions | conflict shrinks ~3× when the CRC is scored symmetrically |
| "no vector satisfies all four", 4.34σ | 28.6% of the gap was recovery being absent from the objective |
| same, 17.69σ for a cooperative scheme | **86.5%** of the gap was the same omission |

The recurring failure is specific and worth naming: **a conclusion about whether two observables
are compatible was repeatedly drawn from an objective that did not contain both of them.**

## 2. What the deactivation waveform requires

### 2.1 An amplitude constraint in single-open schemes

Relaxation components were obtained analytically from the generator matrix (eigenvalues of Qᵀ for
row-vector propagation; the decomposition reconstructs the true relaxation to 6.7×10⁻¹⁵, and an
independent reimplementation agrees on fractional amplitudes to 6.7×10⁻¹⁶).

Maximising the fastest component's fractional amplitude a₁ over each scheme's full parameter
box, admitting only decompositions with Σ|aᵢ|/|Σaᵢ| ≤ 1.5 (i.e. at most 50% cancellation between
modes), gives:

| scheme | states | conducting open states | free params | max admissible a₁ | cancellation at optimum |
|---|---|---|---|---|---|
| JW95, all rates free | 5 | 1 | 6 | **−0.00000** | 1.000 (none) |
| open-state desensitisation | 5 | 1 | 6 | **−0.00000** | 1.000 (none) |
| pre-open ("flip") | 6 | 1 | 8 | +0.245 | 1.500 (at the limit) |
| **two open states** | 6 | **2** | 8 | **+1.116** | 1.500 (at the limit) |

Measured requirement: **a₁ = 0.57**.

In both five-state single-open schemes a₁ is constrained to be non-positive, with no cancellation
at all. Because the fractional amplitudes sum to one, a scheme forced to a₁ = 0 incurs an exact
and irreducible penalty of (0.57/0.04)² = **203.0625** on the amplitude term, which is observed
to 2 d.p. in every single-open scheme tested.

**The negative control matters more than the positive one.** The pre-open scheme has the *same*
eight free parameters as the two-open scheme and does not relieve the constraint: it reaches
a₁ = 0.245, still short of 0.57, and only at the cancellation limit. So this is a property of
the topology, not of parameter count — adding a shut state upstream of the gate does not help;
adding a second conducting state does.

### 2.2 The constraint is a theorem, not an observation

At zero agonist **no transition enters** the doubly-liganded set {A2R, A2O, A2D}, so its 3×3
sub-generator governs the open-state relaxation exactly. That subset is an acyclic tree, so
detailed balance holds unconditionally and `D = diag(√π)` symmetrizes the sub-generator to a real
symmetric tridiagonal matrix. A Rayleigh quotient with test vector `(0, −√β, √α)` gives
`|λ₁| > α + β`, and the open-state eigenvector row then yields in closed form

  w_A2O / |w_A2R| = α / (|λ₁| − α) < α/β = 1/E,

**unconditionally, for all positive rate constants.** Since the fastest mode's amplitude is
a₁ = p_A2O(0)·w_A2O − p_A2R(0)·|w_A2R| + p_A2D(0)·w_A2D, non-positivity follows whenever the
initial condition satisfies p_A2R(0)/p_A2O(0) ≥ w_A2O/|w_A2R|, which holds with strict inequality
at gating equilibrium.

**The reachable-initial-condition step is also proven**, by an argument that needs no
monotonicity assumption. At zero agonist probability enters the doubly-liganded set **only**
through the gateway shut state A2R. By Duhamel's formula every reachable p(0) is therefore a
non-negative combination of impulse responses ρ(u) = e_A2Rᵀ e^{Bu}, u ≥ 0. The weight vector is
proportional to B's right eigenvector for λ₁, so

  a₁(ρ(u)) = e^{λ₁u} · w_A2R < 0 for every u ≥ 0,

because the gateway weight w_A2R is strictly negative. Since a₁ is linear, it is negative on the
entire reachable cone. The equilibrium case follows from an exact secular identity obtained from
the Schur complement of H,

  β/(|λ₁| − α) + d/(|λ₁| − r) − 1 = −2k_off/|λ₁| < 0,

which gives a₁^eq = −(2k_off/|λ₁|)·|w_A2R|·p_A2R^eq < 0 for all positive rates, and shows that
the quantity controlling the sign is precisely the leak out of the doubly-liganded set by
unbinding. The same construction yields the result for the open-state-desensitisation topology.

Both identities were verified numerically across 400 random parameter sets with no violations
(2000 checks of the Duhamel identity, 400 of the secular identity), and every premise of the
proof is mirrored as an executable assertion in the test suite.

### 2.3 Two conducting open states are *necessary*, and that is also a theorem

The converse direction is now proven as well, and by a sharper route than the amplitude bound.

For a single-open **five-state** scheme the doubly-liganded block has dimension 3, so the
relaxation has at most **three** non-zero modes. The extremal mode's amplitude is strictly
negative (§2.2). Therefore **at most two modes can carry positive amplitude.**

The measurement reports **three** positive fractional amplitudes — 0.57 ± 0.04, 0.23 ± 0.02,
0.20 ± 0.03. A five-state single-open scheme cannot produce three positive decay components at
any parameter values whatsoever. **Observing them necessitates at least two conducting open
states.**

A sufficient condition is also available: a closed-form gateway projection gives a spectral sign
partition, and a root in the required interval is guaranteed whenever γ+δ > r and
2k_off/(γ+δ) + d/(γ+δ−r) > 1 — conditions the GABA-A rate ranges satisfy comfortably, so the
published a₁ = 0.57 is reachable rather than merely not excluded.

**Limits.** The necessity argument is dimensional and therefore scoped to **five-state**
single-open schemes. A six-state single-open scheme such as the pre-open variant has a
four-dimensional block and could in principle supply three positive components; it is excluded
by measurement (a₁ ≤ 0.245, and only with condition numbers above 10⁸), not by this theorem. The admissibility
criterion is no longer arbitrary: a decomposition is admitted only if the model's waveform
supports a well-conditioned three-exponential fit, which is the procedure the source
experimenters applied. Under that criterion the pre-open scheme's apparent a₁ = +0.245 disappears
— reaching it requires near-degenerate eigenvalues with condition numbers above 10⁸ — and the
qualitative table is unchanged.

## 3. A scheme that satisfies five observables

Scoring five independent observables, each against its own published values and errors —
deactivation components, EC50 and n_H, paired-pulse recovery, peak open probability band, and
mean open time — with no relative weighting introduced:

| scheme | k | n_data | residual dof | five-way minimax | all five ≤ 3σ |
|---|---|---|---|---|---|
| JW95, all rates free | 6 | 10 | +4 | 14.25σ | no |
| two open states | 8 | 10 | +2 | 3.23σ | no |
| **two open states + cooperativity** | 9 | 10 | **+1** | **2.58σ** | **yes** |

Positive residual degrees of freedom matter here: an earlier four-observable version of this fit
had nine parameters against nine data and was exactly saturated, so its success carried no
evidential weight. Mean open time was added as a fifth observable specifically to recover a
degree of freedom.

**That recovery does not survive an honest accounting of the band constraints, and the
manuscript states so rather than relying on it.** Two of the ten data — `po_max` and mean open
time — are bands, not measurements with errors. Scoring a band by treating its half-width as a
σ is what yields `res_dof = +1`. Under the alternative reading, a uniform prior over the band
with effective σ = half-width/√3, a flat prior carries **zero Fisher information in the
interior**: `po_max` = 0.767 sits interior and so constrains nothing, while mean open time is
*active* on its upper bound. Counting only the eight genuine point measurements gives
**res_dof = 0 or −1**. The fit still passes numerically at 2.58σ, and `kinetic_two_open` still
fails at 3.28σ, so the *ranking* is robust — but **`kinetic_coop` should not be described as an
over-determined statistical test.** It is a scheme that is consistent with five observables, not
one that has survived a test with degrees of freedom to spare.

### 3.1 Mean open time, and a convention that affects the verdict

Mean open time was taken from Keramidas & Harrison 2008 on α1β2γ2S — the same subunit
combination as the deactivation source. Its M-mode components are
τ = (0.49 ± 0.02, 2.58 ± 0.09, 5.17 ± 0.38) ms with areas (0.26, 0.49, 0.25), giving an
area-weighted mean open time of **2.684 ms**; the source does not state whether its errors are
SD or SEM.

**The absence of an error convention was verified against the full text** (PubMed Central
PMC2213567): searching the complete article for `s.e.m`, "standard error", "standard deviation",
`SD` and `mean ±` returns no genuine occurrence. The band treatment is therefore not a
conservative choice under ambiguity but the only treatment the source supports; assigning a
Gaussian σ would require inventing a convention the authors never published.

The observable was scored as a band [2.54, 2.83] ms propagated from the published component
errors, rather than as a Gaussian with an assumed σ.
The fitted scheme predicts **3.20 ms**, outside that band. **The verdict depends on this
choice:**

* band treatment (excess beyond the edge, over the half-width): 2.55σ — **passes**
* point treatment (deviation from the midpoint, over the half-width): 3.56σ — **fails**

The model misses the measured mean open time by 0.52 ms, a 19% relative error. The band
treatment is defensible precisely because the source is silent on its error kind, but the result
should be reported as "2.58σ under a band treatment of an unstated error", never as "2.58σ".

Two further caveats. Keramidas reports an M-mode and an H-mode, which in the single-channel
literature generally denotes modal gating rather than two open states within a burst; the two
modes were therefore treated as separate band constraints and **not** as evidence for the
two-open-state topology, which would have been circular. A γ2L source (Li et al. 2008, weighted
mean 2.96 ms) sits closer to the model's prediction and was deliberately **not** preferred, being
a different splice variant at sub-saturating agonist.

## 4. The cooperativity is not tightly determined

The cooperative scheme gives the second binding step its own dissociation rate; c = koff/koff2 is
the site-to-site ratio, and because kon is shared it is directly a ratio of dissociation
constants, comparable with equilibrium binding data.

Profiling c over 23 log-spaced points from 1 to 10,000 with four-pass bidirectional
warm-starting (`scripts/coop_profile.py`, 18.4 min):

| c | minimax | | c | minimax |
|---|---|---|---|---|
| 1 | 3.232σ | | 100 | 2.751σ |
| 5 | 3.078σ | | 158 | 2.709σ |
| **23** | **3.013σ** | | 1000 | 2.585σ |
| **35** | **2.909σ** ✓ | | **1585** | **2.373σ** |
| 50 | 2.841σ ✓ | | 2512 | 2.868σ ✓ |
| 75 | 2.745σ ✓ | | 3981 | 4.314σ |

An asymmetric trough with a **plateau over c ∈ [35, 2500]**, a shallow minimum at c = 1585, and a
steep wall above 2500 bound by mean open time and deactivation.

This matters for interpretation. A single minimax optimum had returned c ≈ 1000, which would
require cooperativity one to two orders of magnitude beyond published site-to-site affinity
ratios. **The profile shows that value was not required** — it is where an optimiser settled
inside a flat basin. The lowest admissible value is c = 35, and at an published ratio of ~23 the
minimax is 3.013σ, missing a 3σ criterion by 0.013σ, which is within optimiser tolerance.

For comparison, and read **as abstracts only**: positive cooperativity between GABA-A agonist
binding sites is reported (Mozrzymas et al. 2003, *J Neurophysiol*), with equilibrium binding
studies describing two affinity classes differing by ~23× (Olsen et al. 1981, *Mol Pharmacol*)
and ~5× (Agey & Dunn 1989, *Biochemistry*). Neither binding study is on α1β2γ2S or at the
relevant temperature, and Olsen interprets the two classes as distinct receptor populations
rather than intramolecular cooperativity. These are an order-of-magnitude orientation, not a
constraint, and were not used as one.

## 5. The Hill slope is the residual tension

Along the same profile:

| c | n_H | z(n_H) | mean open time | peak singly-bound occupancy |
|---|---|---|---|---|
| 1 | 0.91 | −3.23σ | 2.22 ms | 8.81% |
| 23 | 0.99 | −3.01σ | 3.27 ms | 2.62% |
| 35 | 1.04 | −2.91σ | 2.83 ms | 2.34% |
| 100 | 1.10 | −2.75σ | 2.52 ms | 1.61% |
| 1000 | 1.18 | −2.55σ | 3.20 ms | 0.55% |
| 10000 | 1.63 | −1.43σ | 1.37 ms | 0.39% |

The maximum Hill slope **anywhere** on the profile is 1.63, reached only at c = 10,000 where the
minimax is 8.10σ and mean open time has collapsed to 1.37 ms. Among admissible fits the maximum
is **1.25**, against a published 2.2 ± 0.4.

**At every passing point the Hill slope residual lies between −2.55σ and −2.91σ.** So the claim
"satisfies five observables" means precisely that the worst residual is just below a 3σ
criterion, and that worst residual is always the Hill slope. Cooperativity steepens n_H
monotonically, as expected, but cannot carry it past ~1.25 without destroying mean open time and
deactivation.

**The ceiling has an analytic explanation, and open states cannot lift it.**

The normalised concentration–response of these schemes depends **only on the binding steps**.
If every conducting state is reached from the doubly-liganded shut state A2R, the equilibrium
open probability is a fixed multiple of p_A2R, and that multiple is independent of agonist
concentration; normalising it away leaves a curve determined entirely by binding. **Adding open
states — a second, a third, in series or in parallel — therefore cannot change the Hill slope at
all.** This is a derivation, not an observation, and it is confirmed numerically: a third
conducting state moves n_H from 1.178 to 1.184.

The achievable range is then bracketed analytically. For two **independent** identical sites with
an all-or-none readout, the Hill slope at EC50 is

  n_H = 2/(1 + 1/√2) = **4 − 2√2 ≈ 1.1716**,

verified numerically to 1.17157. In the **fully cooperative** limit the two sites bind
concertedly and n_H → **2.000** exactly. So the entire reachable interval for any two-site scheme
of this family is (1.1716, 2.000), with cooperativity the only parameter that moves within it.

Two consequences follow, and they sharpen the tension rather than resolving it.

* The published 2.2 ± 0.4 has a 1σ interval of [1.8, 2.6]. The fully cooperative asymptote of
  2.000 lies **inside** it, at −0.5σ. The measurement is therefore **not** inconsistent with a
  two-site scheme in principle. **This is a falsification boundary:** a Hill slope reliably
  above 2.000 would exclude the entire two-site class, not a parameter choice. Making that test
  decisive at 3σ requires measuring n_H to σ ≤ 0.067 — about a sixfold improvement on the
  published precision, and the single most informative experiment this analysis points to.
* But the jointly-fitted value is **1.18** — within 0.01 of the *non-cooperative* value 1.1716,
  despite c ≈ 1000. **The cooperativity is not being spent on the concentration–response.** The
  deactivation and mean-open-time constraints consume it, leaving the CRC essentially where an
  entirely non-cooperative scheme would put it.

An earlier version of this section reported unconstrained maxima of 2.93 to 12.06 and concluded
that every scheme can reach n_H = 2.2 in isolation. That is arithmetically true and
scientifically misleading: those values are attained at EC50 values far from the measured
11.6 ± 0.9 µM. Constraining EC50 to a physiological window collapses the maximum to ≤ 1.18.

**The concentration–response steepness is therefore the irreducible tension in this model class**
— not deactivation, not recovery, not peak open probability — and it is a conflict *between*
observables rather than a limit of any one scheme. Stated as a prediction: schemes of this
family should fall roughly 2.5σ short on Hill slope when fitted jointly.

### 5.1 A discriminating measurement — SOUGHT AND NOT FOUND

An earlier version of this section proposed the 10–90% current risetime under rapid 10 µM GABA
application as a measurement that would resolve the cooperativity plateau, quoting 3.54 ms at
c = 23 against 9.20 ms at c = 1000 — an apparent 2.6× separation.

**That claim is withdrawn.** Those were two single optima drawn from a flat basin. Propagating
across the *admissible region* at each cooperativity gives

| c | predicted 10–90% risetime |
|---|---|
| 23 | [3.05, 3.60] ms |
| 1000 | [1.71, 15.24] ms |

The c = 23 interval is **entirely contained** within the c = 1000 interval. The minimum
resolvable difference is **0.00 ms**: no experimental precision can discriminate, because the
model's own predictions overlap completely. Ultrafast application (60–100 µs exchange) is far
more than fast enough; the limitation is the model, not the method.

This is recorded rather than deleted because the error is instructive and is the same one this
paper documents elsewhere: **a point estimate quoted from a flat likelihood basin, mistaken for
a determination.** The cooperativity plateau means the fitted data cannot determine c; it also
means predictions derived from a single point on that plateau cannot be trusted. No
discriminating measurement has been identified.

## 6. Claims retracted during this work

Reported here because the retractions are part of the result, and because two of them are
instances of the same failure the paper is about.

1. **An "amplitude obstruction" at χ² = 452.99.** Withdrawn: this is an arithmetic identity, the
   minimum of the published χ² over amplitudes (0, x, 1−x) under a sum-to-one normalisation, at
   452.9856. Its invariance across topologies — which had been offered as evidence that the
   constraint was structural — is exactly what an arithmetic floor produces.
2. **"The time constants fit exactly; the obstruction is in the amplitudes."** Withdrawn: true
   only on the analytic route, and the canonical route (replicating the source's own
   three-exponential fitting procedure) gives a time-constant residual far from zero. For the
   reference scheme the direction even reverses.
3. **"No scheme satisfies all four observables", minimax 4.34σ.** Withdrawn: paired-pulse
   recovery was not in the objective. It was evaluated after optimisation, so the reported figure
   was an upper bound on the achievable worst-case residual, not a minimum.

A fourth error was caught before it propagated: an interval for the first dissociation rate,
derived by inverting a non-cooperative five-state scheme against the concentration–response
curve, was nearly used to exclude the cooperative model. That inference is invalid across
topologies. Independently, the cooperative fit's **second** dissociation rate falls inside that
interval — which is the quantity the concentration–response curve constrains once cooperativity
is present, and is a genuine consistency check between two independently derived calculations.

## 7. Methods

Model schemes, objectives and search procedures are in `src/circuitpharm/`; every figure quoted
above is printed by a script in `scripts/` and recorded in the corresponding `RESULTS-BRIEF-*.md`.
`STATE-OF-RESULTS.md` indexes each claim with its status, its number, the script that prints it,
and its principal caveat.

Three practices are worth stating as method rather than detail, because each caught a real error:

* **Nested models are asserted, not assumed.** Where one scheme reduces to another in a limit,
  the embedding is checked to reproduce the sub-model's likelihood numerically, the enlarged fit
  is warm-started from the sub-model's optimum, and a nested fit scoring worse is treated as a
  failed optimisation rather than a finding.
* **Degrees of freedom are reported as n_data, k_fitted and residual dof,** with p-values given
  only where residual dof > 0. An exact fit at k ≥ n is reported as an existence statement.
* **Constrained compositions are scored as such.** Three fractional amplitudes summing to one
  carry two degrees of freedom, not three; scoring them as independent inflates χ² and was the
  source of retraction (1).

## 8. What this does not establish

* No scheme examined here reproduces the published Hill slope while satisfying the other
  observables; the best jointly-admissible value is ~1.18 against 2.2 ± 0.4. Open states cannot
  help — that is proven. Cooperativity can, in principle, reach 2.000, which is inside the
  published interval; it does not, because the other observables consume it.
* The amplitude constraint of §2 is **proven**, for both the weights and the reachable initial
  conditions, for the two single-open topologies stated. What is **not** proven is the
  complementary positive claim: that a second conducting open state *must* lift the constraint.
  That direction rests on measurement (a₁ = +1.116) and a conjecture about mode phases.
* The five-observable result rests on a band treatment of mean open time. The absence of a stated
  error kind is verified from the full text, so the band is source-mandated; but a reader who
  insists on a point treatment would see a 3.56σ failure, and that disagreement is about
  inference, not about the measurement.
* **`kinetic_coop` is consistent with five observables; it has not survived an over-determined
  test.** Under a defensible reading of the two band constraints its residual degrees of freedom
  are 0 or −1 (§3.1).
* **No discriminating measurement for the cooperativity plateau has been found** (§5.1). The
  previously proposed risetime experiment does not discriminate.
* All searches are bounded. Statements of the form "no vector was found" are not proofs of
  non-existence.
* Two of the ten data are band constraints scored with half-width as σ, so residual dof = +1
  overstates how over-determined the fit is.
* The binding-cooperativity comparisons are drawn from abstracts of studies on different
  preparations.
