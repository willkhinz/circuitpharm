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
