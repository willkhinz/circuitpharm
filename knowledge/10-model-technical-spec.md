# Technical Specification: Next-Generation GABA-A Kinetic & Neuropharmacology Model

**Scope:** Architectural and algorithmic design for the computational model in `src/circuitpharm/`.  
**Consensus Date:** 2026-10-08 (via structured design tree alignment).

---

## 1. Architectural Overview & Modular Layout

To move from an isolated 5-state Markov calculation to a rigorously tested, data-constrained, multi-model neuropharmacology platform without breaking existing tests (277 passing), the codebase is organized into clean modular subpackages under `src/circuitpharm/`:

```
src/circuitpharm/
├── models/
│   ├── base.py                 # Abstract ReceptorModel protocol & base interface
│   ├── operational.py          # Model A: Phenomenological Operational Scalar Ceiling
│   ├── kinetic_jw95.py         # Model B: 5-state Jones & Westbrook Kinetic Allostery
│   └── extended_desens.py      # Model C: Dual-pathway desensitization & state-dependent PAM
├── fitting/
│   ├── data.py                 # Digitized patch-clamp & concentration-response datasets
│   ├── identifiability.py      # Profile likelihood & Fisher Information Matrix
│   └── mcmc.py                 # Pure NumPy/SciPy affine-invariant ensemble MCMC
├── protocols/
│   ├── waveforms.py            # Transient, pulse-train & ambient spillover generators
│   ├── decomposition.py        # Log-additive & Shapley mechanistic attribution
│   └── oed.py                  # Optimal experimental design & falsification protocol
└── gabaa_kinetics.py           # Preserved backwards-compatible façade wrapping Model B
```

---

## 2. Model Hierarchy (Models A, B, C)

All models implement the polymorphic `ReceptorModel` protocol (`models/base.py`):
```python
class ReceptorModel(Protocol):
    def steady_state(self, gaba_um: float, pam_factor: float = 1.0) -> float: ...
    def simulate_waveform(self, t: np.ndarray, gaba_t: np.ndarray, pam_factor: float = 1.0) -> WaveformResult: ...
    def dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0) -> np.ndarray: ...
```

### Model A: Operational Scalar Ceiling (`models/operational.py`)
* Phenomenological operational model of allosterism (Black & Leff / modern operational Hill).
* Agonist response: $E([G]) = E_{\max} \frac{[G]^n}{[G]^n + EC_{50}^n}$.
* Modulator acts via empirical potency-shift factor $s \le s_{\max}$ and scaling factor $\alpha_{\max}$, imposing an extrinsic efficacy ceiling blind to microscopic state flux.

### Model B: Kinetic Allostery Scheme (`models/kinetic_jw95.py`)
* 5-state Del Castillo–Katz / Jones–Westbrook (1995) Markov scheme:
  $$R \underset{k_{\text{off}}}{\overset{2 k_{\text{on}} [G]}{\rightleftharpoons}} AR \underset{2 k_{\text{off}}}{\overset{k_{\text{on}} [G]}{\rightleftharpoons}} A_2R \underset{\alpha}{\overset{\beta}{\rightleftharpoons}} A_2O$$
  with desensitization from the shut doubly-bound state: $A_2R \underset{r}{\overset{d}{\rightleftharpoons}} A_2D$.
* PAM acts via accelerated binding / slowed unbinding ($k_{\text{off}} \leftarrow k_{\text{off}} / s$) and/or gating acceleration ($\beta \leftarrow \beta \cdot g$), with zero ad-hoc gain caps.

### Model C: Extended Desensitization & State-Dependent Affinity (`models/extended_desens.py`)
* Motivated by the *Nature Communications* GABAA desensitization pathway.
* Dual desensitization branches: fast ($D_{\text{fast}}$) and slow ($D_{\text{slow}}$) conformational sinks.
* Modulator binds with state-dependent affinity across shut ($R$), open ($O$), and desensitized ($D$) states via a thermodynamic cycle, testing whether differential stabilization of desensitized states alters the phasic-tonic headroom divergence.

---

## 3. Mechanistic Decomposition Engine (`protocols/decomposition.py`)

Isolates the three physical contributors to the phasic vs. tonic modulation difference:
1. **Agonist Occupancy ($\Delta \ln \text{Occ}$):** Fraction of receptors bound by GABA.
2. **Kinetic Redistribution ($\Delta \ln \text{Kinet}$):** Conformational equilibrium between bound-shut, open, and desensitized states ($\beta/\alpha$, $d/r$).
3. **Operating Point Headroom ($\Delta \ln \text{Headroom}$):** Proximity to the maximal open probability ceiling $P_{\text{o,max}} = \beta/(\alpha+\beta)$.

**Algorithm:**
* Exact factorial grid decomposing total logarithmic gain:
  $$\ln G = \Delta \ln(\text{Occupancy}) + \Delta \ln(\text{Kinetic Redistribution}) + \Delta \ln(\text{Operating Point Headroom})$$
* Computes exact Shapley attributions (% of total log-gain) for both synaptic transients and steady-state ambient GABA.

---

## 4. Parameter Estimation, Identifiability & Uncertainty (`fitting/`)

* **Datasets (`fitting/data.py`):**
  * Digitized concentration-response curves (EC₅₀, Hill slope).
  * Rapid perfusion macroscopic deactivation & desensitization current traces (decay $\tau_1, \tau_2$).
  * Single-channel open probability and mean burst durations.
* **Profile Likelihood (`fitting/identifiability.py`):**
  * Evaluates profile likelihood paths for all 6 microscopic parameters.
  * Formally classifies parameters into structurally identifiable, practically identifiable, or non-identifiable manifold directions.
* **Affine-Invariant MCMC Ensemble (`fitting/mcmc.py`):**
  * Lightweight pure NumPy/SciPy Goodman & Weare ensemble sampler (no heavy PPL dependencies).
  * Produces joint posterior parameter chains $\mathcal{P}(\mathbf{\theta} \mid \mathcal{D})$.
  * Validates that the predicted phasic-tonic headroom divergence holds across the full 95% posterior credible interval.

---

## 5. Realistic Waveform & Electrophysiology Engine (`protocols/waveforms.py`)

Provides standardized physiological stimulation protocols:
1. **Biexponential Synaptic Pulse:** Variable rise ($\tau_{\text{rise}} = 0.1\text{–}0.5\text{ ms}$) and biexponential clearance ($\tau_{\text{fast}} = 1\text{ ms}, \tau_{\text{slow}} = 10\text{–}30\text{ ms}$).
2. **High-Frequency Pulse Trains:** 10 Hz, 50 Hz, 100 Hz bursts to quantify frequency-dependent desensitization accumulation ($P_D$) and paired-pulse depression.
3. **Activity-Dependent Spillover:** Baseline ambient GABA with transient spillover elevations.

**Standard Output Metrics:**
* Peak open probability and peak current: $P_{\text{peak}}$, $I_{\text{peak}} = g_{\max} P_{\text{peak}} (V - E_{\text{Cl}})$.
* Integrated inhibitory charge: $Q = \int I(t)\,dt$.
* Biexponential decay time constants: $\tau_{\text{fast}}, \tau_{\text{slow}}, \tau_{\text{weighted}}$.
* Paired-pulse recovery ratio: $I_2 / I_1$.

---

## 6. Three-Tiered Dynamic Range Metric (`DynamicRangeEvaluation`)

Replaces single headline point estimates with a structured metric dataclass:
1. **Theoretical Asymptotic Headroom:** $P_{\text{open},\infty} / P_0$ as $k_{\text{off}} \to 0^+$ at stationary ambient $[G]$.
2. **Pharmacologically Reachable Gain:** Gain evaluated at ligand-specific operational efficacy ceilings ($s_{\max} = 1.25, 1.50, 2.50$).
3. **Physiological Charge Transfer Ratio:** Realized ratio $\Delta Q_{\text{tonic}} / \Delta Q_{\text{phasic}}$ under realistic driving forces and time-varying waveforms.
* All quantities report median and [2.5%, 97.5%] posterior credible intervals, alongside parameter boundary scans showing where the ratio collapses ($\le 1.0$).

---

## 7. Model Discrimination & Optimal Experimental Design (`protocols/oed.py`)

* **Model Selection:** Computes AIC, BIC, and cross-validated out-of-sample log-likelihood across training and held-out validation datasets for Models A, B, and C.
* **Optimal Experimental Design (OED):**
  * Optimizes agonist concentration $[G]$, PAM concentration $[PAM]$, and pulse duration $\Delta t$ to maximize predictive divergence (e.g. Jensen-Shannon divergence relative to experimental noise $\sigma_{\text{noise}}$).
  * Outputs a concrete experimental testing protocol with quantitative, non-overlapping falsification boundaries for each model hypothesis.
