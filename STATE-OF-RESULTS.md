# State of Results — Consolidated Index of Findings

**Date:** 2026-10-17  
**Branch:** `claude/agent-brief-16`  
**Purpose:** Single-page reference index defining the live, provisional, resolved, and withdrawn status of all major project findings, with their supporting numbers, scripts, and primary caveats.

---

| Item | Claim / Topic | Status | Supporting Number(s) | Producing Script | Main Caveat |
| :--- | :--- | :---: | :--- | :--- | :--- |
| **1** | **Pseudo-data leverage imbalance** | **LIVE** | 101 trace points generated from 6 published numbers created 4.3× to 10× artificial leverage over recovery scalar in weightings (i)–(iii); eliminated in (iv). | `scripts/pseudodata_leverage.py` | Model must be evaluated through the same three-exponential fitting pipeline as the experiment (`predict_components_fit` or eigen-decomposition matching published components), never by scoring raw waveforms directly. |
| **2** | **Forced $a_1 \le 0$ & $(0.57/0.04)^2$ amplitude floor** | **LIVE** | 5-state single-open schemes structurally forced to $a_1 \le 0$ under admissible decomposition ($L_1 \le 1.5$), enforcing an exact arithmetic $\chi^2$ floor of $(0.57/0.04)^2 = \mathbf{203.0625}$. Two open states (`kinetic_two_open`) break it ($a_1 = +1.12$). | `scripts/amplitude_admissibility.py`<br>`scripts/deact_obstruction_decomposition.py` | Established by bounded Latin hypercube box optimization under the $L_1 \le 1.5$ mode-cancellation convention; analytical proof from infinitesimal generator structure is not yet formulated. |
| **3** | **`kinetic_coop` five-way fit** | **PROVISIONAL**<br>*(Convention-dependent)* | **$2.58\sigma$** under band evaluation of Keramidas 2008 mean open time ($[2.54, 2.83]$ ms, passes $\le 3\sigma$); **$3.56\sigma$** under midpoint evaluation ($2.684 \pm 0.145$ ms, fails $\le 3\sigma$). $\text{res\_dof} = +1$. | `scripts/joint_fiveway_fit.py` | Pass hinges entirely on treating Keramidas 2008 unstated error kind as an empirical flat band rather than an invented Gaussian SEM. Model misses midpoint by $0.52$ ms (19% relative error). |
| **4** | **`FINDING-amplitude-obstruction.md`** | **WITHDRAWN** | $\chi^2 = 452.99$ was identical across schemes. | *None (withdrawn)* | Pure arithmetic identity arising from sum-to-one amplitude normalization; carried zero biophysical or kinetic information. |
| **5** | **`FINDING-minimax-conflict.md`** | **WITHDRAWN** | Claimed $17.69\sigma$ irreconcilable four-way conflict. | *None (withdrawn)* | Paired-pulse recovery was omitted from the optimizer objective and evaluated only downstream; incorporating recovery directly dropped minimax to $2.39\sigma$. |
| **6** | **`FIT_TARGETS['ec50_um']` correction** | **RESOLVED** | $20.0\ \mu\text{M} \to \mathbf{11.6 \pm 0.9\ \mu\text{M}}$ (Jahn et al. 1997 peak CRC, $9.33\sigma$ discrepancy eliminated). | `src/circuitpharm/fitting/crc_components.py`<br>(`surface_ec50_discrepancy()`) | Earlier project fits against $20.0\ \mu$M targeted an unverified literature midpoint; all current analyses enforce the sourced $11.6 \pm 0.9\ \mu$M and $n_H = 2.2 \pm 0.4$. |
| **7** | **Cooperativity $c$ profile** | **LIVE** | $c = 1000$ is **an artifact of an under-constrained plateau**, not a sharp minimum. $c = 50$ achieves **$2.841\sigma$**; $c = 35$ achieves **$2.909\sigma$**; $c = 23$ achieves **$3.013\sigma$** (misses by $0.013\sigma$). Plateau spans $c \in [35, 2500]$. | `scripts/coop_profile.py`<br>(`scripts/coop_profile_results.json`) | While the extreme requirement $c = 1000$ largely evaporates, values $c \le 23$ fail to satisfy the five-way minimax at $\le 3\sigma$ due to low CRC Hill slope ($n_H \le 0.99$ vs $2.2 \pm 0.4$). |

---

## Notes on Taxonomy and Conventions
- **Settled vs Provisional:** A finding is *settled* only if its conclusion is invariant to permissible evaluation conventions and holds positive residual degrees of freedom. `kinetic_coop` remains *provisional* because its pass depends on the MOT band convention.
- **Negative Literature Rules:** Binding cooperativity literature sources (Olsen 1981, Agey & Dunn 1989, Mozrzymas 2003) remain registered as `where="abstract"` in `src/circuitpharm/provenance.py` and are strictly barred from use as numerical fitting constraints.
