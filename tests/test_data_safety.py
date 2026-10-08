"""The knowledge-base builders must not be able to destroy accumulated data.

WHY THIS FILE EXISTS. `scripts/build_kb.py` rebuilds the pharmacology database from the
literal tables written into it: 42 sources and 17 findings. The LIVE database had
accumulated 90 sources and 58 findings across later sessions, added by other means. The
script's `if DB.exists(): DB.unlink()` therefore silently discarded 48 sources and 41
findings -- and before the `__main__` guard was added, merely IMPORTING the script did it.

That actually happened during session 8, while checking for exactly this class of bug. It
was recoverable only because the database is tracked in git.

These tests are cheap and they protect real, hard-to-reproduce data.
"""
import pathlib
import sqlite3

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
DB = REPO / "data" / "pharmacology.db"
BUILDERS = ["build_kb.py", "build_compounds.py"]

# These tests inspect the SOURCE of the builder scripts, so they need the repository
# checkout rather than an installed package. Found by running the suite against a clean
# `pip install` in /tmp, where tests/ was present and scripts/ was not -- they errored
# instead of skipping, which would fail CI for anyone testing an installed copy.
HAVE_SCRIPTS = (REPO / "scripts").is_dir()
needs_repo = pytest.mark.skipif(
    not HAVE_SCRIPTS,
    reason="needs the repository checkout (inspects scripts/ source), not an installed copy")


@pytest.mark.skipif(not DB.exists(), reason="no database in this checkout")
def test_database_still_holds_the_accumulated_rows():
    """A canary. If these counts drop, something rebuilt over the live data."""
    c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    assert c.execute("select count(*) from sources").fetchone()[0] >= 90
    assert c.execute("select count(*) from findings").fetchone()[0] >= 58
    assert c.execute("select count(*) from model_facts").fetchone()[0] >= 28


@needs_repo
@pytest.mark.parametrize("name", BUILDERS)
def test_builders_have_a_destruction_guard(name):
    """Any script that drops or unlinks must first compare against the live row count."""
    src = (REPO / "scripts" / name).read_text()
    destructive = ("unlink()" in src or "DROP TABLE" in src)
    if not destructive:
        pytest.skip(f"{name} does not drop or unlink")
    assert "REFUSING" in src, (
        f"{name} destroys data without a guard. It must compare the live row count "
        "against what it will write and refuse unless --force.")
    assert "--force" in src, f"{name} offers no explicit override"


@needs_repo
@pytest.mark.parametrize("name", BUILDERS)
def test_builders_do_nothing_on_import(name):
    """Before session 8 these ran at module level, so importing them rebuilt the database.
    Any tool that imports or scans the package would have triggered it."""
    src = (REPO / "scripts" / name).read_text()
    assert '__name__ == "__main__"' in src, f"{name} lacks a __main__ guard"


@pytest.mark.skipif(not DB.exists(), reason="no database in this checkout")
def test_verification_status_is_explicit_not_implied():
    """No source may imply its CLAIM has been verified when only its metadata has.

    THIS TEST EARNED ITS PLACE AND THEN NEEDED REWRITING, which is worth recording. It used
    to assert `count(verification = 'UNVERIFIED') > 0` -- a tripwire against quietly marking
    everything verified. `scripts/verify_sources.py` then resolved all 90 sources by
    DOI/PMID and replaced that literal with `METADATA_<status>; CLAIM_SUPPORT_UNASSESSED`,
    so the count went to zero and the test fired.

    It was right to fire: the string had changed and nothing else would have noticed. But
    the invariant it was defending is not "the word UNVERIFIED appears somewhere" -- it is
    **that resolving a citation must not be mistakable for verifying a claim**. Those are
    different questions, and the first source audited by hand in this project proves it:
    `a5_dist` resolves at similarity 1.00 and does NOT support the numbers attributed to it
    (a 1988 generic alpha-subunit probe measuring regional LEVEL, where the model needs
    per-subtype COMPOSITION). A pipeline that wrote "VERIFIED" on a resolved DOI would have
    marked that one verified.

    So the assertion is now on the invariant itself: every source must carry an explicit
    claim-support state, and a claim may only be called assessed where a hand verdict with
    reasoning exists (the `HAND:` sentinel).
    """
    c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    for tbl in ("sources", "findings"):
        cols = [r[1] for r in c.execute(f"PRAGMA table_info({tbl})")]
        assert "verification" in cols, f"{tbl} has no verification column"

    rows = c.execute("select key, verification from sources").fetchall()
    silent = [k for k, v in rows if not (v or "").strip()]
    assert not silent, f"sources with an empty verification state: {silent}"

    # Every source states its claim-support position, one way or the other.
    unstated = [k for k, v in rows
                if "CLAIM_SUPPORT" not in (v or "") and "UNVERIFIED" not in (v or "")]
    assert not unstated, (
        f"sources whose verification says nothing about CLAIM support: {unstated}. "
        f"Metadata resolution is not claim verification; see a5_dist.")

    # A claim may only be marked assessed where a hand verdict is recorded. Anything else
    # claiming assessment is the exact conflation this test exists to prevent.
    assessed = [k for k, v in rows if "CLAIM_SUPPORT_UNASSESSED" not in (v or "")
                and "UNVERIFIED" not in (v or "")]
    for k in assessed:
        v = dict(rows)[k]
        assert "HAND:" in v, (
            f"source {k!r} implies its claim was assessed without a HAND: verdict "
            f"recording who assessed it and how: {v[:120]}")
