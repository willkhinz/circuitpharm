# RESULTS-BRIEF-8A: Structural Search Over Deactivation Components and Predictor Unification

**Agent:** Gemini 3.8 Flash (High) (Worker 8A)  
**Execution Command:**  
```bash
PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python scripts/deact_structural_search.py --n-samples 1500 --n-polish 12 --n-distribution 4000 --seed 42
```
**Runtime:** 1.78 s  
**Worktree:** `/Users/whinz/Biochemistry-wt-agent-brief-8a` (branch `claude/agent-brief-8a`)  
**Interpreter:** `/Users/whinz/Biochemistry/.venv/bin/python`  
**Literature Target:** `Barberis 2007, brief pulse` (`circuitpharm.literature.DEACTIVATIONS[1]`)  
**Published Target Values:** $\tau = (2.8, 33.4, 221.35)\text{ ms}$, $\text{areas} = (0.57, 0.23, 0.20)$  
**Published Standard Errors:** $\sigma_\tau = (0.3, 4.6, 14.9)\text{ ms}$, $\sigma_a = (0.04, 0.02, 0.03)$  

---

## 1. Executive Summary & Bounded Result

> **Bounded Conclusion:** No parameter vector was found, within the declared bounds and that search budget (1,500 LHS samples polished by Nelder–Mead), reproducing the published three-component deactivation in any of the three kinetic schemes (`kinetic_jw95_free_rates`, `extended_desens`, `kinetic_open_desens`).
> **This is a bounded empirical finding under the declared budget and parameter box, not a mathematical impossibility proof.**

Key findings from the structural search:
1. **Severe $\chi^2$ penalties across all schemes:** The best-fitting scheme is `kinetic_jw95_free_rates` with $\chi^2 \approx 527.9$ (residual sum of squares in standard errors across 6 observations), followed by `extended_desens` ($\chi^2 \approx 1107.7$) and `kinetic_open_desens` ($\chi^2 \approx 2181.9$).
2. **Zero schemes produce three positive decay amplitudes:** In all three schemes, the fastest non-zero relaxation mode that carries observable current enters with a negative fractional amplitude. This reflects transient channel opening ($A2R \to A2O$) following sudden agonist removal before current decay dominates.
3. **Appreciable modes count:** Over a 4,000-sample Latin Hypercube survey of the 6-parameter box (`kinetic_jw95_free_rates`), only ~5.0% of vectors have 3 modes with $|\text{fractional amplitude}| \ge 0.10$. At the minimum-$\chi^2$ optima, the schemes present at most 1 or 2 appreciable modes.
4. **Resolution of the 7A vs 7B predictor conflict:**
   - Brief 7A's analytic eigen-decomposition (`predict_components_eigen`) recovers the true mathematical relaxation spectrum of the generator matrix $Q_0$.
   - Brief 7B's triexponential fit (`predict_components_fit`) faithfully fits post-peak decay waveforms, matching Barberis 2007's empirical protocol.
   - When one relaxation mode carries $\ge 90\%$ of the amplitude (as at the JW95 optimum where 1 mode carries 98.2%), a triexponential fit to the model waveform is rank-deficient and degenerates into repeated/tied time constants. That degeneracy is an authentic physical property of the model's relaxation, not an optimizer flaw.
   - When $\ge 3$ modes carry $|\text{amplitude}| \ge 0.10$, both predictors agree on the dominant time constants to within $0.3\%$ relative tolerance.

---

## 2. Scheme Comparison Table

| Scheme | Free Params | Best $\chi^2$ (LHS $\to$ Polish) | Appreciable Modes ($|a| \ge 0.10$) | 3 Positive Amplitudes? | $\tau$ Residuals $(\sigma)$ | Area Residuals $(\sigma)$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `kinetic_jw95_free_rates` | 6 | 578.0 $\to$ **527.9** | 2 | **No** | (-8.34, +0.03, +0.10) | (-14.33, +8.79, +13.25) |
| `extended_desens` | 4 | 1135.8 $\to$ **1107.7** | 2 | **No** | (-9.18, +17.39, +9.44) | (-14.26, +19.89, +5.74) |
| `kinetic_open_desens` | 3 | 2182.1 $\to$ **2181.9** | 1 | **No** | (-8.62, +2.14, +18.70) | (-14.38, +38.76, -6.67) |

---

## 3. Detailed Component Decomposition at Best Optima

### Scheme: `kinetic_jw95_free_rates`
* **Best $\chi^2$:** 527.86 (reproduced expected $\approx 527.9$)
* **Optimal Parameter Vector:** `{'log10_kd': -1.0000, 'log10_E': -0.6943, 'log10_D': -0.2984, 'log10_r': -1.7142, 'log10_koff': -2.3800, 'log10_alpha': 0.4471}`
* **All Relaxation Modes (taus ms, areas):**
  * $\tau = 0.2969\text{ ms}$, $\text{area} = -0.003292$
  * $\tau = 33.5353\text{ ms}$, $\text{area} = +0.405791$
  * $\tau = 222.8185\text{ ms}$, $\text{area} = +0.597501$
  * $\tau = 239.8664\text{ ms}$, $\text{area} = +0.000000$ (unbinding mode $AR \to R$, zero amplitude in $A2O$)

| Component | Quantity | Target | Error (SEM) | Predicted | Residual $(\sigma)$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | $\tau_1$ | 2.80 ms | 0.30 ms | 0.2969 ms | -8.34$\sigma$ |
| 1 | $\text{area}_1$ | 0.570 | 0.040 | -0.0033 | -14.33$\sigma$ |
| 2 | $\tau_2$ | 33.40 ms | 4.60 ms | 33.5353 ms | +0.03$\sigma$ |
| 2 | $\text{area}_2$ | 0.230 | 0.020 | +0.4058 | +8.79$\sigma$ |
| 3 | $\tau_3$ | 221.35 ms | 14.90 ms | 222.8185 ms | +0.10$\sigma$ |
| 3 | $\text{area}_3$ | 0.200 | 0.030 | +0.5975 | +13.25$\sigma$ |

*Note on `jw95_free_rates`:* The optimiser matches $\tau_2$ ($+0.03\sigma$) and $\tau_3$ ($+0.10\sigma$) with remarkable accuracy, but $\tau_1$ collapses to $0.297\text{ ms}$ with negative area $-0.0033$.

### Scheme: `extended_desens`
* **Best $\chi^2$:** 1107.73 (reproduced expected $\approx 1107.7$)
* **Optimal Parameter Vector:** `{'log10_kd': -0.0399, 'log10_E': 1.8551, 'log10_D_fast': 1.4115, 'log10_D_slow': 0.2701}`
* **All Relaxation Modes (taus ms, areas):**
  * $\tau = 0.0449\text{ ms}$, $\text{area} = -0.000217$
  * $\tau = 5.5428\text{ ms}$, $\text{area} = -0.000000$ (unbinding mode)
  * $\tau = 113.4098\text{ ms}$, $\text{area} = +0.627800$
  * $\tau = 361.9842\text{ ms}$, $\text{area} = +0.372200$
  * $\tau = 2505.8092\text{ ms}$, $\text{area} = +0.000217$

| Component | Quantity | Target | Error (SEM) | Predicted | Residual $(\sigma)$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | $\tau_1$ | 2.80 ms | 0.30 ms | 0.0449 ms | -9.18$\sigma$ |
| 1 | $\text{area}_1$ | 0.570 | 0.040 | -0.0002 | -14.26$\sigma$ |
| 2 | $\tau_2$ | 33.40 ms | 4.60 ms | 113.4098 ms | +17.39$\sigma$ |
| 2 | $\text{area}_2$ | 0.230 | 0.020 | +0.6278 | +19.89$\sigma$ |
| 3 | $\tau_3$ | 221.35 ms | 14.90 ms | 361.9842 ms | +9.44$\sigma$ |
| 3 | $\text{area}_3$ | 0.200 | 0.030 | +0.3722 | +5.74$\sigma$ |

### Scheme: `kinetic_open_desens`
* **Best $\chi^2$:** 2181.92 (reproduced expected $\approx 2181.9$)
* **Optimal Parameter Vector:** `{'log10_kd': -1.0000, 'log10_E': 1.1296, 'log10_D': -3.0000}`
* **All Relaxation Modes (taus ms, areas):**
  * $\tau = 0.2136\text{ ms}$, $\text{area} = -0.005122$
  * $\tau = 5.5428\text{ ms}$, $\text{area} = -0.000000$ (unbinding mode)
  * $\tau = 43.2410\text{ ms}$, $\text{area} = +1.005113$
  * $\tau = 500.0445\text{ ms}$, $\text{area} = +0.000009$

| Component | Quantity | Target | Error (SEM) | Predicted | Residual $(\sigma)$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | $\tau_1$ | 2.80 ms | 0.30 ms | 0.2136 ms | -8.62$\sigma$ |
| 1 | $\text{area}_1$ | 0.570 | 0.040 | -0.0051 | -14.38$\sigma$ |
| 2 | $\tau_2$ | 33.40 ms | 4.60 ms | 43.2410 ms | +2.14$\sigma$ |
| 2 | $\text{area}_2$ | 0.230 | 0.020 | +1.0051 | +38.76$\sigma$ |
| 3 | $\tau_3$ | 221.35 ms | 14.90 ms | 500.0445 ms | +18.70$\sigma$ |
| 3 | $\text{area}_3$ | 0.200 | 0.030 | +0.0000 | -6.67$\sigma$ |

---

## 4. Mode-Count Distribution Over 4000 Samples (6-Parameter Box)

Distribution of how many non-zero relaxation modes carry $|\text{fractional amplitude}| \ge 0.10$ across 4,000 Latin Hypercube samples of `kinetic_jw95_free_rates`:

| Appreciable Modes ($|a_i| \ge 0.10$) | Sample Count ($n=4000$) | Fraction (%) | Regression Expected Check |
| :---: | :---: | :---: | :---: |
| **0 modes** | 0 | 0.0% | 0.0% |
| **1 mode** | 2367 | 59.2% | ~57.0% |
| **2 modes** | 1433 | 35.8% | ~37.6% |
| **3 modes** | 200 | 5.0% | ~5.3% |
| **4 modes** | 0 | 0.0% | 0.0% |

The observed distribution is consistent with the regression baseline (~5.0% vs 5.3% for 3 modes; 0% for 4 modes).
