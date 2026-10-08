#!/usr/bin/env python3
"""Read-only query helper for the pharmacology database.

OPENED READ-ONLY (`mode=ro`). This used to connect read-write while accepting arbitrary
SQL from the command line, so `kb.py sql "DROP TABLE sources"` would have worked. Given
that this project has already destroyed 48 sources and 41 findings once by running a
rebuild script, a query tool that can write is not a theoretical risk.
Query helper for data/pharmacology.db. Usage: kb.py <view> | kb.py sql "<query>" """
import sqlite3, sys, pathlib, textwrap

DB = pathlib.Path(__file__).resolve().parent.parent / "data" / "pharmacology.db"
VIEWS = {
 "findings": """SELECT f.id, f.confidence AS conf, f.claim, f.bearing, s.url
                FROM findings f LEFT JOIN sources s ON s.key=f.source_key ORDER BY f.id""",
 "targets":  """SELECT target, subunit, direction,
                       conc_low_mM||'-'||COALESCE(conc_high_mM,'inf') AS mM,
                       magnitude, relevant_at_social_dose AS rel_5_20mM,
                       evidence_quality AS qual, note
                FROM ethanol_targets ORDER BY relevant_at_social_dose DESC, target""",
 "disc":     """SELECT training_drug, test_drug, pct_substitution AS pct,
                       substitutes, note FROM discrimination""",
 "subunits": """SELECT receptor, subunit, region, level, bearing FROM subunit_expression
                ORDER BY receptor, subunit""",
 "eeg":      """SELECT study_key, dose_route, bec_measured AS bec, strain,
                       state_scored AS scored, bands, aperiodic AS aper, key_result
                FROM eeg_studies""",
 "model":    """SELECT item, value, verified_how FROM model_facts""",
 "compounds": """SELECT name, arm, target, status, eth_like, sedation, resp,
                        hepatotox, ceiling, verdict, note FROM compounds
                 ORDER BY CASE verdict WHEN 'CANDIDATE' THEN 1 WHEN 'CANDIDATE (adjunct)' THEN 2
                 WHEN 'BACKUP' THEN 3 WHEN 'TOOL' THEN 4 ELSE 5 END, name""",
 "cand":      """SELECT name, arm, target, status, eth_like, sedation, verdict
                 FROM compounds WHERE verdict LIKE 'CANDIDATE%' OR verdict='BACKUP'""",
 "chiral":    """SELECT pair, eutomer, eutomer_action, distomer, distomer_action,
                        kind, ratio, relevance FROM chiral_pairs
                 ORDER BY CASE kind WHEN 'SWITCH' THEN 1 WHEN 'SAFETY' THEN 2
                 WHEN 'OPPOSITE' THEN 3 ELSE 4 END""",
 "sim":      """SELECT condition, mn_hz, pct_control AS pct, period_ms AS per,
                       duty, flex_ext_corr AS corr, note FROM sim_results""",
 "sources":  """SELECT key, year, kind, citation, url FROM sources ORDER BY year""",
}

def show(sql):
    c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True); c.row_factory = sqlite3.Row
    rows = c.execute(sql).fetchall()
    if not rows: print("(no rows)"); return
    for i, r in enumerate(rows, 1):
        print(f"\n\033[1m[{i}]\033[0m")
        for k in r.keys():
            v = r[k]
            if v is None or v == "": continue
            body = textwrap.fill(str(v), 94, subsequent_indent=" " * 18)
            print(f"  {k:<15} {body}")
    print(f"\n{len(rows)} rows")

if __name__ == "__main__":
    a = sys.argv[1:]
    if not a or a[0] in ("-h", "--help"):
        print("views:", ", ".join(VIEWS), "\n   or: kb.py sql \"SELECT ...\""); sys.exit(0)
    if a[0] == "sql": show(" ".join(a[1:]))
    elif a[0] in VIEWS: show(VIEWS[a[0]])
    else: print(f"unknown view '{a[0]}'. views: {', '.join(VIEWS)}"); sys.exit(1)
