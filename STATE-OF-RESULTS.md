# State of Results — Consolidated Index of Findings

**Date:** 2026-10-18 (updated by Claude Opus 5 after briefs 17–19 merged)  
**Branch:** `main`  
**Purpose:** Single-page reference index defining the live, provisional, resolved, and withdrawn status of all major project findings, with their supporting numbers, scripts, and primary caveats.

---

| Item | Claim / Topic | Status | Supporting Number(s) | Producing Script | Main Caveat |
| :--- | :--- | :---: | :--- | :--- | :--- |
| **1** | **Pseudo-data leverage imbalance** | **LIVE** | 101 trace points generated from 6 published numbers created 4.3× to 10× artificial leverage over recovery scalar in weightings (i)–(iii); eliminated in (iv). | `scripts/pseudodata_leverage.py` | Model must be evaluated through the same three-exponential fitting pipeline as the experiment (`predict_components_fit` or eigen-decomposition matching published components), never by scoring raw waveforms directly. |
| **2** | **Forced $a_1 \le 0$ & $(0.57/0.04)^2$ amplitude floor** | **LIVE** | 5-state single-open schemes structurally forced to $a_1 \le 0$ under admissible decomposition ($L_1 \le 1.5$), enforcing an exact arithmetic $\chi^2$ floor of $(0.57/0.04)^2 = \mathbf{203.0625}$. Two open states (`kinetic_two_open`) break it ($a_1 = +1.12$). | `scripts/amplitude_admissibility.py`<br>`scripts/deact_obstruction_decomposition.py` | **Now PROVEN** (`knowledge/15-amplitude-constraint-theorem.md`): at zero agonist nothing enters the doubly-liganded set, so its sub-generator is exact and symmetrizable; $\|\lambda_1\| > \alpha+\beta$ gives $w_{A2O}/\|w_{A2R}\| < 1/E$ unconditionally, and Duhamel confines reachable states to a cone on which $a_1 < 0$. The $L_1$ convention is retired in favour of fit-conditioning. **Unproven:** the converse, that two open states *must* lift it. |
| **3** | **`kinetic_coop` five-way fit** | **PROVISIONAL**<br>*(Convention-dependent)* | **$2.58\sigma$** under band evaluation of Keramidas 2008 mean open time ($[2.54, 2.83]$ ms, passes $\le 3\sigma$); **$3.56\sigma$** under midpoint evaluation ($2.684 \pm 0.145$ ms, fails $\le 3\sigma$). $\text{res\_dof} = +1$. | `scripts/joint_fiveway_fit.py` | Pass hinges on the band treatment — but the **absence of a stated error kind is now verified against the full text** (PMC2213567: no occurrence of SEM, SD, standard error/deviation or `mean ±`), so the band is source-mandated rather than chosen. Model still misses the midpoint by $0.52$ ms (19%). Two of ten data are bands scored half-width-as-$\sigma$, which inflates `res_dof`. |
| **4** | **`FINDING-amplitude-obstruction.md`** | **WITHDRAWN** | $\chi^2 = 452.99$ was identical across schemes. | *None (withdrawn)* | Pure arithmetic identity arising from sum-to-one amplitude normalization; carried zero biophysical or kinetic information. |
| **5** | **`FINDING-minimax-conflict.md`** | **WITHDRAWN** | Claimed $17.69\sigma$ irreconcilable four-way conflict. | *None (withdrawn)* | Paired-pulse recovery was omitted from the optimizer objective and evaluated only downstream; incorporating recovery directly dropped minimax to $2.39\sigma$. |
| **6** | **`FIT_TARGETS['ec50_um']` correction** | **RESOLVED** | $20.0\ \mu\text{M} \to \mathbf{11.6 \pm 0.9\ \mu\text{M}}$ (Jahn et al. 1997 peak CRC, $9.33\sigma$ discrepancy eliminated). | `src/circuitpharm/fitting/crc_components.py`<br>(`surface_ec50_discrepancy()`) | Earlier project fits against $20.0\ \mu$M targeted an unverified literature midpoint; all current analyses enforce the sourced $11.6 \pm 0.9\ \mu$M and $n_H = 2.2 \pm 0.4$. |
| **7** | **Cooperativity $c$ profile** | **LIVE** | $c = 1000$ is **an artifact of an under-constrained plateau**, not a sharp minimum. $c = 50$ achieves **$2.841\sigma$**; $c = 35$ achieves **$2.909\sigma$**; $c = 23$ achieves **$3.013\sigma$** (misses by $0.013\sigma$). Plateau spans $c \in [35, 2500]$. | `scripts/coop_profile.py`<br>(`scripts/coop_profile_results.json`) | While the extreme requirement $c = 1000$ largely evaporates, values $c \le 23$ fail to satisfy the five-way minimax at $\le 3\sigma$ due to low CRC Hill slope ($n_H \le 0.99$ vs $2.2 \pm 0.4$). |

---

## Notes on Taxonomy and Conventions
- **Settled vs Provisional:** A finding is *settled* only if its conclusion is invariant to permissible evaluation conventions and holds positive residual degrees of freedom. `kinetic_coop` remains *provisional* because its pass depends on the MOT band convention.
- **Negative Literature Rules:** Binding cooperativity literature sources (Olsen 1981, Agey & Dunn 1989, Mozrzymas 2003) remain registered as `where="abstract"` in `src/circuitpharm/provenance.py` and are strictly barred from use as numerical fitting constraints.

---

## Added after briefs 17–19 (Claude Opus 5, 2026-10-18)

| Item | Claim / Topic | Status | Supporting Number(s) | Producing Script | Main Caveat |
| :--- | :--- | :---: | :--- | :--- | :--- |
| **8** | **Amplitude constraint is a theorem** | **LIVE (proven)** | $w_{A2O}/\|w_{A2R}\| = \alpha/(\|\lambda_1\|-\alpha) < 1/E$ unconditionally; Duhamel cone gives $a_1(\rho(u)) = e^{\lambda_1 u} w_{A2R} < 0$; equilibrium secular identity $\beta/(\|\lambda_1\|-\alpha) + d/(\|\lambda_1\|-r) - 1 = -2k_{off}/\|\lambda_1\|$. Verified: 0 violations in 2000 + 400 numerical checks. | `knowledge/15-amplitude-constraint-theorem.md`<br>`tests/test_amplitude_constraint_theorem.py` (13 tests) | Covers the two single-open topologies stated. The positive direction (two open states *can*) remains a labelled conjecture. |
| **9** | **Open states cannot change the Hill slope** | **LIVE (derived)** | Normalised CRC depends only on binding when all conducting states are reached from A2R. Third open state moves $n_H$ only $1.178 \to 1.184$. | `scripts/three_open_hill_max.py` | Applies to schemes where every open state hangs off A2R; a conducting state reachable from a *different* shut state is not covered. |
| **10** | **Analytic Hill-slope bracket** | **LIVE (derived)** | Independent two sites: $n_H(EC_{50}) = 4-2\sqrt{2} = 1.17157$. Fully cooperative limit: exactly $2.000$. Any two-site scheme confined to $(1.172, 2.000)$. Jointly fitted value $1.18$ — within $0.01$ of the **non-cooperative** value despite $c \approx 1000$. | verified analytically and numerically; see `RESULTS-BRIEF-18.md` §5 | Published $2.2 \pm 0.4$ straddles the upper end: $2.000$ is inside the 1σ interval at $-0.5\sigma$. **Falsification boundary:** $n_H$ reliably $> 2.0$ would exclude the whole model class. |
| **11** | **Three open states** | **LIVE (negative)** | Series and parallel both reach `res_dof = -1` with a parameter pinned at bound; minimax $2.54\sigma$ vs coop's $2.58\sigma$. | `src/circuitpharm/models/kinetic_three_open.py` | Saturated, so not a test. Adds nothing and costs testability. Keramidas's three open *components* are burst-level dwell-time fits, not necessarily scheme states. |
| **12** | **Keramidas 2008 full text read** | **LIVE** | No error convention stated anywhere in the article. Source reports **three** open components for α1β2γ2S in both gating modes. | `FINDING-keramidas-fulltext.md` | First **external** corroboration that one open state is insufficient — from dwell-time analysis, independent of the macroscopic amplitudes fitted here. Modal gating (differing shut-state counts between modes) is unrepresented by any scheme in this project. |

**Test baseline on `main`, measured 2026-10-18: 614 fast passed, 0 skipped.**
Note `pytest-xdist` is active and has produced inconsistent totals across runs; confirm counts
rather than transcribing them.
