# Results Brief 21 — Digitisation Pipeline & Data Provenance Gate

**Date:** 2026-10-19  
**Agent:** Antigravity Worker (brief21)  
**Worktree:** `/Users/whinz/Biochemistry-wt-agent-brief-21` (`claude/agent-brief-21`)  

---

## 1. Task 1 — Attempt Retrieval (Verbatim Records)

Target sources:
- Barberis et al. 2007, *Eur J Neurosci* 25(9):2726–40, **PMC1950087**, PMID 17561840, DOI 10.1111/j.1460-9568.2007.05530.x.
- Jahn et al. 1997, *NeuroReport* 8(16):3443–6, PMID 9427304 (no PMC record known; paywalled / closed access).

### 1.1 Locally Available Tools for PMC1950087

#### Attempt 1: EuropePMC Tool (`europepmc_api.py download_pdf`)
Command executed:
```bash
uv run ~/.gemini/config/plugins/science/skills/literature_search_europepmc/scripts/europepmc_api.py download_pdf PMC1950087 --output /tmp/pmc1950087.pdf
```
Outcome (verbatim error):
```
Downloading PDF for PMC1950087...
Traceback (most recent call last):
  ...
urllib.error.HTTPError: HTTP Error 403: Forbidden
polite_http.http_client.HttpError: HTTP Error 403 while fetching https://europepmc.org/articles/PMC1950087?pdf=render
Server body: <!DOCTYPE html><html lang="en-US"><head><title>Just a moment...</title><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta http-equiv="X-UA-Compatible" content="IE=Edge"><meta name="robots" content="noindex,nofollow"><meta name="viewport" content="width=device-width,initial-scale=1">...
```
Result: **FAILED (HTTP 403 Forbidden / Cloudflare challenge)**.

#### Attempt 2: PubMed Central BioC Full-Text Tool (`pubmed_api.py get_full_text_pmc`)
Command executed:
```bash
uv run ~/.gemini/config/plugins/science/skills/pubmed_database/scripts/pubmed_api.py /tmp/barberis_pmc.json get_full_text_pmc "17561840"
```
Inspection:
```python
import json, re
with open('/tmp/barberis_pmc.json') as f:
    data = json.load(f)
text = data.get('full_text', '')
urls = re.findall(r'https?://[^\s<>"\'\)]+', text)
images = re.findall(r'[^\s]+\.(?:jpg|jpeg|png|tif|tiff|gif)', text, re.IGNORECASE)
```
Output:
- Text length: 68,189 characters.
- URLs count: 0.
- Image file extensions count: 0.
Result: **TEXT ONLY**. As predicted by Opus in the brief, the BioC API tool provides plain text passages only; no figure URLs, no figure captions, and no image files.

#### Attempt 3: NCBI PMC Open Access Web Service (`oa.fcgi`)
Command executed:
```python
import urllib.request
url = 'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC1950087'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
resp = urllib.request.urlopen(req)
```
Outcome (verbatim error):
```
urllib.error.HTTPError: HTTP Error 404: Not Found
```
Result: **FAILED (HTTP 404 Not Found)**. PMC1950087 is an author manuscript / publisher-deposited paper in PMC, but is NOT part of the PMC Open Access subset that distributes image archives.

#### Attempt 4: Jahn et al. 1997 (*NeuroReport* 8(16):3443–6)
- PubMed search confirms PMID 9427304.
- No PMC copy exists.
- The article is closed access / paywalled behind Wolters Kluwer / Ovid.
Result: **NOT RETRIEVED**. No open-access retrieval route exists.

---

### 1.2 Search for Figures on Local Machine

Command executed across repo, `data/`, and download directories:
```bash
find . -iname "*barberis*" -o -iname "*jahn*" -o -iname "*1950087*" -o -iname "*.png" -o -iname "*.jpg" -o -iname "*.tif*" -o -iname "*.pdf"
find /Users/whinz/Biochemistry/data -type f
```
Outcome:
- Files found in worktree:
  - `./tests/test_config_single_source.py`
  - `./scripts/jahn_options.py`
  - `./src/circuitpharm/config.py`
- Files found in `data/`:
  - `pharmacology.db` (and backups)
- Image / PDF files count: **0**.
Result: **ZERO figure images or PDFs exist on disk**.

### 1.3 Task 1 Summary

In accordance with Section 0 of Agent Brief 21 ("If you cannot obtain the actual published figure as an image file on disk, you must NOT produce any digitised data point. Report kind=blocked and build the machinery in §2 instead"):
- **Zero digitised data points have been produced or fabricated.**
- Neither Barberis 2007 nor Jahn 1997 figure images could be retrieved legitimately via existing tools.
- The gate remains closed, protecting the integrity of all empirical analyses in this project.

---

## 2. Task 2 — Digitisation Pipeline & Evidence Gate

The pipeline machinery has been built in `src/circuitpharm/fitting/digitisation.py` and exported through `src/circuitpharm/fitting/__init__.py`.

### 2.1 Implemented Components

1. **Evidence Verification Data Structures:**
   - `FigureIdentifier`: Enforces non-empty `journal`, `doi`, `figure_number`, and `panel`.
   - `AxisReference` & `AxisCalibration`: Records reference coordinates (2 reference points per axis, in both pixel and data coordinates) supporting linear and log10 scales. Provides forward (`pixel_to_data_x`, `pixel_to_data_y`) and inverse transformations with strict non-degeneracy validation.
   - `DigitisationEvidence`: Encapsulates on-disk `image_path`, recorded `image_sha256`, `figure_id`, `calibration`, `method`, and `operator`. Hard failure (`MissingImageError`, `ChecksumMismatchError`) is triggered upon construction or inspection if the file does not exist on disk or the SHA-256 digest fails to match.
   - `compute_file_sha256(path)` & `verify_image_evidence(path, sha256)`.

2. **Digitised Dataset Loaders:**
   - `load_digitised_deactivation(...)`: Constructs `DeactivationDataset` with `kind="digitised"`, `synthetic=False`, and verified `evidence`. Accepts data directly or raw pixel coordinates transformed via `AxisCalibration`.
   - `load_digitised_dose_response(...)`: Constructs `DoseResponseDataset` with `kind="digitised"`, `synthetic=False`, and verified `evidence`. Accepts data directly or raw pixel coordinates transformed via `AxisCalibration`.

3. **Registration Gate:**
   - `assert_digitised_evidence(dataset)` & `register_digitised_dataset(name, dataset, registry)`: Enforce that any dataset claiming `kind="digitised"` and `synthetic=False` physically cannot be registered into `ALL` without attached, verified `DigitisationEvidence`.
   - `_validate_dataset` in `src/circuitpharm/fitting/data.py` calls `assert_digitised_evidence` whenever `evidence` is present.

4. **Self-Validation Check Against Literature Records (Zero Hardcoded Literals):**
   - `validate_barberis_deactivation(dataset)`: Re-fits deactivation waveform using `predict_components_fit` and verifies that fitted time constants $(\tau_1, \tau_2, \tau_3)$ and fractional areas $(A_1, A_2, A_3)$ match published values from `circuitpharm.literature.DEACTIVATIONS[1]` within published errors:
     - $\tau = (2.8 \pm 0.3, 33.4 \pm 4.6, 221.35 \pm 14.9)$ ms
     - areas $= (0.57 \pm 0.04, 0.23 \pm 0.02, 0.20 \pm 0.03)$
   - `validate_jahn_crc(dataset)`: Fits the concentration-response curve and verifies that $EC_{50}$ and $n_H$ match published values from `get_jahn_published_records()`:
     - $EC_{50} = 11.6 \pm 0.9\ \mu$M
     - $n_H = 2.2 \pm 0.4$
   - `assert_self_validation(dataset)`: Raises `SelfValidationError` if validation fails.

### 2.2 Test Suite (`tests/test_digitisation.py`)

16 dedicated tests pass covering:
- Refusal without image (`MissingImageError`).
- Refusal on checksum mismatch (`ChecksumMismatchError`).
- Refusal on empty operator/method or degenerate calibration.
- Acceptance of valid deactivation and CRC records.
- Bidirectional pixel-to-data calibration accuracy.
- Registration refusal when evidence is missing.
- Self-validation acceptance of true literature parameters.
- Self-validation rejection of deliberately perturbed curves (labeled as test fixtures, NOT real data):
  - Perturbed $\tau_1 = 4.5$ ms ($+5.67\sigma$) caught.
  - Perturbed $A_1 = 0.35$ ($-5.50\sigma$) caught.
  - Perturbed $EC_{50} = 20.0\ \mu$M ($+9.33\sigma$) caught.
  - Perturbed $n_H = 1.2$ ($-2.50\sigma$) caught.

---

## 3. Task 3 — Per-Claim Dependency Analysis (Quantifying What the Gate Blocks)

We evaluate each of the project's current live findings against the three dependency categories:
1. **Synthetic curves** (e.g. 101-point deactivation reconstruction `BARBERIS2007_DEACT`, 9-point CRC `JAHN1997_PEAK_CRC`).
2. **Published parameters** (e.g. Barberis Table 1 components, Jahn abstract $EC_{50}$ and $n_H$, Barberis paired-pulse scalar $R=0.33$, Keramidas MOT).
3. **Neither** (pure mathematical theorems about model topology/generator structure).

| Finding / Claim | Category | Depends on Synthetic Curves? | Depends on Published Parameters? | Would conclusion change if real digitised curves replaced reconstructions? |
| :--- | :---: | :---: | :---: | :--- |
| **1. Pseudo-data leverage imbalance** (`scripts/pseudodata_leverage.py`) | **Synthetic Curves (diagnostic) & Published Parameters (solution)** | **Yes (as diagnostic artifact)** | **Yes (as resolution)** | **NO.** The diagnosis demonstrates that fitting against 101 synthetic curve points created 4.3× to 10× artificial leverage over scalar recovery; the resolution (Weighting iv) scores directly against the 6 published parameters. If real digitised curves replaced the reconstructions, naive point-wise fitting would still suffer from sample-density leverage imbalance (hundreds of points vs 1 scalar). The finding that waveforms must be evaluated by exponential decomposition / eigen-matching rather than raw point-wise SSE is invariant. |
| **2. Amplitude constraint theorem** (`knowledge/15-amplitude-constraint-theorem.md`) | **Neither** | **NO** | **NO** | **NO.** The theorem proves analytically that for 5-state single-open topologies at zero agonist, the fastest relaxation mode amplitude satisfies $a_1 \le 0$ unconditionally across all positive rate vectors and all reachable end-of-pulse initial occupancies. This is an exact mathematical deduction from continuous-time Markov generator structure and Cauchy interlacing of Jacobi matrices. |
| **3. Necessity theorem** (Two open states required to lift $a_1 \le 0$) | **Published Parameters** | **NO** | **YES** | **NO.** The structural non-positivity $a_1 \le 0$ is a mathematical theorem (depends on *neither*). The empirical obstruction ($\chi^2 \ge (0.57/0.04)^2 = 203.0625$) arises because Barberis Table 1 reports a positive fast amplitude $a_1 = +0.57 \pm 0.04$. It does not depend on the synthetic curve points. Any faithful digitisation of Figure 1 must reproduce $a_1 \approx 0.57 > 0$ (enforced by the self-validation gate), so the necessity of two open states is unaffected. |
| **4. Analytic Hill-slope bracket** (`scripts/three_open_hill_max.py`) | **Neither (bracket) / Published Parameters (tension)** | **NO** | **YES (for residual)** | **NO.** The analytical bracket $(1.1716, 2.000)$ is a mathematical theorem derived from binding polynomial factorization when all open states branch off $A_2R$: $[O_\text{total}] = E_\text{eff}[A_2R]$, where $E_\text{eff}$ cancels identically in normalized response $Y(x)$, bounding $n_H(EC_{50})$ between $4 - 2\sqrt{2} \approx 1.1716$ and $2.000$. The empirical tension (model achieves $1.18$ vs published $2.2 \pm 0.4$) targets Jahn's published abstract parameters, not the 9 synthetic curve points. |
| **5. `kinetic_coop` five-observable consistency** (`scripts/joint_fiveway_fit.py`) | **Published Parameters** | **NO** | **YES** | **NO.** Under Weighting (iv) extended, all five observables in the minimax objective are evaluated against published parameters: (1) Barberis Table 1 components ($\tau_i, a_i$); (2) Jahn 1997 parameters ($EC_{50}, n_H$); (3) Barberis paired-pulse scalar ($R=0.33$); (4) $P_{o,\text{max}}$ band; (5) Keramidas MOT band. Neither the 101-point deactivation curve nor the 9-point CRC enters the joint likelihood. Replacing the synthetic reconstructions with real digitised curves leaves the five-way minimax score ($2.58\sigma$) and optimal parameters completely unchanged. |

### 3.1 Plain Summary of Gate Impact
The provenance gate blocks the use of synthetic point arrays as empirical measurements of receptors. However, because the project's primary inferential framework (Weighting iv) was previously restructured to evaluate models directly against published exponential components and Hill parameters rather than raw synthetic waveforms, **none of the project's live theoretical or empirical conclusions depend on the synthetic curves**.

If and when real digitised figure images are provided and ingested:
1. The synthetic flag on `BARBERIS2007_DEACT` and `JAHN1997_PEAK_CRC` will switch from `True` to `False`, eliminating `ProvisionalResultWarning` in point-wise checks.
2. The self-validation gate will guarantee that digitised curves reproduce published parameters within published errors.
3. The five-way minimax fit ($2.58\sigma$), the amplitude constraint theorem ($a_1 \le 0$), the necessity theorem, and the Hill slope ceiling ($1.18$) will remain completely intact.

---

## 4. Task 4 — Test Suite Verification

Both test tiers were executed twice in `/Users/whinz/Biochemistry-wt-agent-brief-21` using `/Users/whinz/Biochemistry/.venv/bin/pytest` with `PYTHONPATH=src`. All tests passed on every run.

### Fast Tier: `pytest -q -m "not slow"`
Baseline on main: 631 passed.
With Brief 21 (16 new digitisation tests in `tests/test_digitisation.py`): 647 passed.

- **Run 1:**
```
647 passed, 32 warnings in 137.13s (0:02:17)
```
- **Run 2:**
```
647 passed, 42 warnings in 136.61s (0:02:16)
```

### Full Tier: `pytest -q`
Baseline on main: 705 passed.
With Brief 21: 721 passed.

- **Run 1:**
```
721 passed, 120 warnings in 904.99s (0:15:04)
```
- **Run 2:**
```
721 passed, 137 warnings in 743.08s (0:12:23)
```

