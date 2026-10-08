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
    """Citations carry a verification column whose DEFAULT is UNVERIFIED.

    The sources have resolvable URLs to real papers, but nobody has gone back to each
    figure or table to confirm the attributed claim is what the paper supports at the
    stated precision. Recording that as UNVERIFIED is honest; omitting the column would
    let a reader assume otherwise -- the same discipline as the VOID tier on quantities.
    """
    c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    for tbl in ("sources", "findings"):
        cols = [r[1] for r in c.execute(f"PRAGMA table_info({tbl})")]
        assert "verification" in cols, f"{tbl} has no verification column"
    n_unver = c.execute(
        "select count(*) from sources where verification='UNVERIFIED'").fetchone()[0]
    assert n_unver > 0, (
        "every source is marked verified; if that is now true, the claim-by-claim check "
        "should be recorded in the worklog with how it was done")
