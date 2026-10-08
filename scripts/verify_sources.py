#!/usr/bin/env python3
"""Verify the 90 knowledge-base sources, mechanically where that is possible.

THE HOLE THIS FILLS. Every source in `data/pharmacology.db` was marked UNVERIFIED -- nobody
had checked that the cited work exists, that the recorded citation matches it, or that it
supports the claim attributed to it. The project's whole factual base rests on them,
including the subtype assignments that drive its one VALIDATED result.

VERIFICATION HAS TWO HALVES AND ONLY ONE IS MECHANICAL. This script does the mechanical
half, honestly labelled, and refuses to pretend to the other:

  METADATA (mechanical, done here)
    Resolve the work by DOI through Crossref or by PMID/title through Europe PMC, then
    compare the resolved title against the recorded citation. Outcomes:
      RESOLVED_MATCH     the work exists and the recorded citation matches it
      RESOLVED_MISMATCH  a work resolves but its title does not match what we recorded
                         -- a possible mis-citation, and the dangerous case
      UNRESOLVED         nothing resolves; the citation may still be correct (old papers,
                         book chapters, industry pages) but it cannot be checked this way
      INTERNAL           this project's own output, not literature; nothing to resolve

  CLAIM SUPPORT (judgement, NOT done here)
    Whether the abstract or full text actually supports the claim attributed to it is a
    reading task. This script RETRIEVES abstracts so that reading is possible, and records
    their availability, but it does not score support. Any source whose claim-support has
    been assessed by reading carries that verdict written by hand, with the reasoning.

WHY THE DISTINCTION IS LOAD-BEARING. A script that resolved a DOI and then wrote "VERIFIED"
would be worse than the UNVERIFIED it replaced: it would convert "nobody checked" into
"checked", while having checked only that the paper exists. The first source audited by hand
in this project resolved perfectly and did NOT support the number attributed to it (see
`a5_dist` in knowledge/06-source-provenance.md). Existence and support are different
questions and this script answers only the first.

Run:  python scripts/verify_sources.py [--write] [--limit N] [--key KEY]
      --write updates the `verification` column; without it, the run is read-only.
"""
import argparse
import datetime
import difflib
import json
import pathlib
import re
import shutil
import sqlite3
import sys
import time
import urllib.parse
import urllib.request

DB = pathlib.Path(__file__).resolve().parent.parent / "data" / "pharmacology.db"
UA = "circuitpharm-source-verifier/1.0 (academic use; contact via repository)"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF = "https://api.crossref.org/works/"


def _get(url, timeout=20):
    req = urllib.request.Request(url, headers={"User-Agent": UA,
                                               "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8", "replace"))


def _norm(t):
    """Normalise a title for comparison: strip markup, greek-letter spellings, punctuation.

    `alpha1` and the literal Greek alpha-1 must compare equal, or every GABA-A subunit paper
    reads as a mismatch -- which would bury the real mis-citations in false positives.
    """
    t = re.sub(r"<[^>]+>", "", t or "").lower()
    for g, a in (("α", "alpha"), ("β", "beta"), ("γ", "gamma"),
                 ("δ", "delta"), ("ε", "epsilon"), ("ö", "o")):
        t = t.replace(g, a)
    t = t.replace("gabaa", "gaba-a").replace("gaba a", "gaba-a")
    return re.sub(r"[^a-z0-9]+", " ", t).strip()


def _similar(a, b):
    return difflib.SequenceMatcher(None, _norm(a), _norm(b)).ratio()


def _classify(resolved_title, citation, floor=0.72):
    """RESOLVED_MATCH / RESOLVED_TRUNCATED / RESOLVED_MISMATCH.

    TRUNCATION IS NOT MIS-CITATION and separating them is the difference between a useful
    report and a noisy one. The first run flagged five mismatches; four were a recorded
    citation that is a PREFIX of the real title (`bk`), an abbreviation of it (`girk`), or a
    summary of the finding used in place of a title (`gl_ii_73`, similarity 0.35). Those are
    citation-hygiene problems -- they cannot be checked mechanically by anyone else -- but
    they are not wrong sources, and lumping them in with a genuine mis-citation hides the
    one case that matters.

    A containment test catches them: if either normalised title contains the other, the
    recorded citation is a fragment of the resolved work rather than a different work.
    """
    nr, nc = _norm(resolved_title), _norm(citation)
    sim = _similar(resolved_title, citation)
    if sim >= floor:
        return "RESOLVED_MATCH", sim
    if nr and nc and (nc[:60] in nr or nr[:60] in nc):
        return "RESOLVED_TRUNCATED", sim
    return "RESOLVED_MISMATCH", sim


def _ids_from_url(url):
    """Pull a DOI and/or PMID out of whatever shape the recorded URL has."""
    # URL-DECODE FIRST. Publisher links routinely carry the DOI percent-encoded, e.g. PLOS
    # writes `?id=10.1371%2Fjournal.pone.0030608`, where the slash is %2F. The DOI regex
    # below needs a literal slash, so without this the DOI is invisible and the source falls
    # through to an unreliable title search and is reported UNRESOLVED.
    #
    # Caught because two sources this audit had just cited as provenance for REGIONS
    # entries -- pbc_delta and pbc_eps -- came back UNRESOLVED. A verification tool that
    # under-reports resolution is not a safe failure: it makes real citations look missing
    # and so makes the provenance gap look worse than it is, which is just as misleading as
    # the opposite.
    url = urllib.parse.unquote(url or "")
    doi = None
    m = re.search(r"(10\.\d{4,9}/[^\s\"'>?#]+)", url)
    if m:
        doi = m.group(1).rstrip(".,;)")
    else:
        # Publishers whose article URLs embed the DOI SUFFIX without the prefix. Nature's
        # `/articles/s41598-017-17379-x` is `10.1038/s41598-017-17379-x`; this covers the
        # Nature portfolio (Sci Rep, Nature, Nat Neurosci, Nat Commun), which supplies
        # several sources here. Without it they fall through to an unreliable title search
        # and are reported UNRESOLVED, overstating the provenance gap.
        m = re.search(r"nature\.com/articles/([a-z0-9\-]+)", url)
        if m:
            doi = "10.1038/" + m.group(1)
    pmid = None
    m = re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", url)
    if m:
        pmid = m.group(1)
    pmcid = None
    m = re.search(r"(PMC\d+)", url)
    if m:
        pmcid = m.group(1)
    return doi, pmid, pmcid


def resolve(citation, url, kind):
    """Try to find the work. Returns (status, resolved_title, abstract_len, how)."""
    if kind == "internal":
        return "INTERNAL", "", 0, "this project's own output"
    # NOT EVERYTHING CITED IS A PAPER. A clinical-trial registry entry, an encyclopedia
    # page, a software repository or a dataset has no DOI and no abstract, and reporting
    # them as UNRESOLVED alongside a paper that genuinely failed to resolve conflates "we
    # could not check this" with "this is not the kind of thing that is checked this way".
    # They still need assessing -- a registry entry can be misdescribed like anything else
    # -- but by following the link, not by DOI resolution.
    if kind in ("trial", "reference", "resource", "industry"):
        return "NON_LITERATURE", "", 0, f"kind={kind}; verify by following the link"
    doi, pmid, pmcid = _ids_from_url(url)

    # 1. DOI through Crossref -- the most authoritative identifier when present
    if doi:
        try:
            m = _get(CROSSREF + urllib.parse.quote(doi, safe="/"))["message"]
            title = (m.get("title") or [""])[0]
            if title:
                status, sim = _classify(title, citation)
                return status, title, 0, f"crossref doi={doi} sim={sim:.2f}"
        except Exception:
            pass

    # 2. PMID / PMCID through Europe PMC, which also carries abstracts
    for q in ([f"EXT_ID:{pmid}"] if pmid else []) + ([f"PMCID:{pmcid}"] if pmcid else []):
        try:
            res = _get(f"{EPMC}?query={urllib.parse.quote(q)}&resultType=core"
                       f"&format=json&pageSize=1").get("resultList", {}).get("result", [])
            if res:
                r = res[0]
                title = r.get("title") or ""
                ab = re.sub(r"<[^>]+>", "", r.get("abstractText") or "")
                status, sim = _classify(title, citation)
                return status, title, len(ab), f"europepmc {q} sim={sim:.2f}"
        except Exception:
            pass

    # 3. Title search, last resort -- a match here is weaker evidence, so the threshold is
    #    higher: we are asking "does a paper with this title exist" rather than confirming
    #    an identifier we already had.
    try:
        q = f'TITLE:"{citation[:120]}"'
        res = _get(f"{EPMC}?query={urllib.parse.quote(q)}&resultType=core"
                   f"&format=json&pageSize=1").get("resultList", {}).get("result", [])
        if res:
            r = res[0]
            title = r.get("title") or ""
            ab = re.sub(r"<[^>]+>", "", r.get("abstractText") or "")
            sim = _similar(title, citation)
            # Title search is weaker evidence than an identifier, so the bar is higher:
            # here we are asking "does a paper with this title exist" rather than
            # confirming an id we already had.
            status, _ = _classify(title, citation, floor=0.85)
            return status, title, len(ab), f"europepmc title sim={sim:.2f}"
    except Exception:
        pass
    return "UNRESOLVED", "", 0, "no doi/pmid resolved and no title match"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true",
                    help="update the verification column (a timestamped backup is kept)")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--key", type=str, default=None)
    ap.add_argument("--force-overwrite", action="store_true",
                    help="overwrite hand-written (HAND:) verdicts too. Almost never right; "
                         "a hand assessment is the only thing here that is not reproducible "
                         "by re-running.")
    ap.add_argument("--sleep", type=float, default=0.34,
                    help="pause between lookups; these are free public APIs")
    a = ap.parse_args()

    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row
    rows = con.execute(
        "select key, citation, url, year, kind, verification from sources "
        + ("where key = ? " if a.key else "") + "order by key",
        (a.key,) if a.key else ()).fetchall()
    if a.limit:
        rows = rows[:a.limit]

    if a.write:
        # Never touch this database without a backup. It has been destroyed once already,
        # and once more in a near miss that only version control caught.
        bk = DB.with_suffix(f".backup-{datetime.datetime.now():%Y%m%d-%H%M%S}.db")
        shutil.copy2(DB, bk)
        print(f"backed up database to {bk.name}\n")

    out, tally = [], {}
    for i, r in enumerate(rows, 1):
        status, title, ablen, how = resolve(r["citation"], r["url"], r["kind"])
        tally[status] = tally.get(status, 0) + 1
        out.append((r["key"], status, title, ablen, how))
        flag = ("  <-- MIS-CITED?" if status == "RESOLVED_MISMATCH"
                else "  <-- citation is a fragment, not the title"
                if status == "RESOLVED_TRUNCATED" else "")
        print(f"[{i:3d}/{len(rows)}] {r['key']:<18} {status:<18} "
              f"abs={ablen:>5}  {how}{flag}", flush=True)
        if status in ("RESOLVED_MISMATCH", "RESOLVED_TRUNCATED"):
            print(f"            recorded: {r['citation'][:92]}")
            print(f"            resolved: {title[:92]}")
        if status != "INTERNAL":
            time.sleep(a.sleep)

    print("\n" + "=" * 78)
    print("METADATA RESOLUTION SUMMARY")
    print("=" * 78)
    for k in ("RESOLVED_MATCH", "RESOLVED_TRUNCATED", "RESOLVED_MISMATCH", "UNRESOLVED",
              "NON_LITERATURE", "INTERNAL"):
        if k in tally:
            print(f"  {k:<20} {tally[k]:3d}")
    print(f"\n  abstracts retrieved: {sum(1 for _,_,_,n,_ in out if n > 0)}")
    print("\nThis says NOTHING about whether each source supports the claim attributed to")
    print("it. That is a reading task; see knowledge/06-source-provenance.md for the ones")
    print("assessed by hand, and note that the first one assessed resolved perfectly and")
    print("did NOT support its attributed number.")

    if a.write:
        now = f"{datetime.datetime.now():%Y-%m-%d}"
        # NEVER CLOBBER A HAND-WRITTEN VERDICT.
        #
        # Found by running this script twice. The first run's output was read by hand, one
        # genuine mis-citation was found (`ganaxolone`: its recorded title and its recorded
        # PMCID are two different papers by the same author) and written into the column
        # with the reasoning. The second run overwrote it with the automatic
        # METADATA_RESOLVED_MISMATCH, discarding the assessment and leaving only the
        # mechanical status.
        #
        # That is the same failure as build_kb.py destroying tables it did not own: an
        # automated writer overwriting curated data, where the loss is silent because the
        # replacement looks like a valid value. A verdict is preserved if it carries the
        # HAND: sentinel, and the count of skipped rows is reported so the preservation is
        # visible rather than assumed.
        kept = []
        for key, status, title, ablen, how in out:
            prior = con.execute("select verification from sources where key = ?",
                                (key,)).fetchone()[0] or ""
            if "HAND:" in prior and not a.force_overwrite:
                kept.append(key)
                continue
            con.execute(
                "update sources set verification = ? where key = ?",
                (f"METADATA_{status} ({now}; {how}); CLAIM_SUPPORT_UNASSESSED", key))
        con.commit()
        print(f"\nwrote metadata verification status for {len(out) - len(kept)} sources")
        if kept:
            print(f"PRESERVED {len(kept)} hand-written verdict(s), not overwritten: "
                  f"{', '.join(kept)}")
        print("CLAIM_SUPPORT_UNASSESSED is recorded explicitly, so a resolved DOI cannot be")
        print("mistaken for a verified claim.")
    else:
        print("\n(read-only run; pass --write to record these in the database)")


# __main__ guard: E8 bit this project twice, once through a script that did destructive work
# at import time -- and this one writes to the database.
if __name__ == "__main__":
    main()
