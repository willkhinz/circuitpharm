"""Verify the manuscript's numbers against a CLEAN CHECKOUT of the tag it cites.

WHY A SEPARATE SCRIPT FROM THE TEST. `tests/test_manuscript_consistency.py` checks the
manuscript against the WORKING TREE. That catches retyped tables, but it cannot catch the
claim the manuscript actually makes, which is that a reader starting from a named tag
reproduces the figures. Those are different assertions, and only the second is the
reproducibility statement. A reviewer put it exactly right: the reproducibility statement is a
testable claim, not a substitute for testing it.

So this script checks out the cited ref into a throwaway git worktree, runs the generators
from THAT tree with its own `src` on the path, and asserts every figure the manuscript quotes
appears in the output of code it did not write.

It also checks the ref contains the manuscript and its tests. The first tag created for this
purpose did not: `manuscript-v2` pointed at the commit BEFORE the manuscript was written, so a
reader following the citation would have found the generator but not the document. The numbers
reproduced; the citation was still wrong, and this check exists because nothing else caught it.

Run:  python scripts/verify_manuscript.py                  # the ref the manuscript cites
      python scripts/verify_manuscript.py --ref HEAD        # or any other
      python scripts/verify_manuscript.py --keep            # leave the worktree for inspection
"""
from __future__ import annotations

import argparse
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

MANUSCRIPT = pathlib.Path("knowledge/08-manuscript.md")

#: Files the cited ref MUST contain for the citation to mean what the manuscript says.
REQUIRED_AT_REF = (
    "knowledge/08-manuscript.md",
    "scripts/paper_numbers.py",
    "scripts/ranking_robustness.py",
    "tests/test_manuscript_consistency.py",
    "src/circuitpharm/gabaa_kinetics.py",
    "src/circuitpharm/subtypes.py",
)

#: Figures quoted in the manuscript, grouped by the reviewer's four audit categories, each as
#: (label, regex to find in the manuscript, how to find the same value in generator output).
#: Written as literal strings rather than recomputed, because recomputing here would just
#: reproduce the generator and prove nothing. The point is a THIRD party: clean-tree output.
CATEGORIES = {
    "calibration targets and fitted rates": [
        "20.0000", "0.7500", "15.000",
        "0.0112168", "0.333062", "0.646493", "0.134142",
        "29.69", "4.8195", "0.20749",
    ],
    "phasic trace: peak, deactivation, integrated open probability": [
        "0.419609", "0.553593", "1.319", "1.895", "2.419", "1.668", "28.426",
    ],
    "steady-state asymptote and finite-modulator tonic gains": [
        "0.156377", "184.6", "7.876", "2.9321",
        "2.21", "3.84", "7.18", "20.25", "56.46", "158.08",
        "2881.1", "725.8", "48.1", "15.0", "4.8",
    ],
    "selectivity indices and Monte Carlo lower tail": [
        "10.87", "8.64", "7.87", "1.20", "0.56", "9.34", "3.62",
        "2.51", "2.49", "2.57", "2.48", "2.55", "63.8", "99.9",
    ],
}


def run(cmd, cwd=None, env=None):
    r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        raise SystemExit(f"FAILED: {' '.join(cmd)}\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    return r.stdout


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", default=None,
                    help="git ref to verify; default is the tag the manuscript cites")
    ap.add_argument("--draws", type=int, default=20000)
    ap.add_argument("--keep", action="store_true", help="do not delete the worktree")
    a = ap.parse_args()

    text = MANUSCRIPT.read_text()

    ref = a.ref
    if ref is None:
        tags = re.findall(r"tag \*\*`(manuscript-v[\d.]+)`\*\*", text)
        if not tags:
            raise SystemExit("the manuscript cites no tag; pass --ref explicitly")
        ref = tags[0]
    print(f"verifying the manuscript against a clean checkout of {ref!r}\n")

    sha = run(["git", "rev-parse", f"{ref}^{{commit}}"]).strip()
    print(f"  {ref} -> {sha[:12]}")

    # ---- the citation must be self-contained -------------------------------------------
    at_ref = set(run(["git", "ls-tree", "-r", "--name-only", ref]).split())
    missing = [f for f in REQUIRED_AT_REF if f not in at_ref]
    print(f"\n[1/4] does {ref} contain what the citation implies?")
    for f in REQUIRED_AT_REF:
        print(f"      {'ok  ' if f in at_ref else 'MISS'} {f}")
    if missing:
        print(f"\n  FAIL: {ref} is missing {len(missing)} file(s). A reader following the "
              f"citation\n        cannot reproduce from it alone.")
        return 1

    # ---- run the generators from the clean tree ----------------------------------------
    tmp = pathlib.Path(tempfile.mkdtemp(prefix="verify-manuscript-"))
    tree = tmp / "tree"
    try:
        run(["git", "worktree", "add", "--detach", str(tree), sha])
        import os
        env = dict(os.environ)
        env["PYTHONPATH"] = str(tree / "src")
        env.pop("PYTHONHOME", None)

        print(f"\n[2/4] environment")
        ver = run([sys.executable, "-c",
                   "import sys,numpy,scipy,circuitpharm;"
                   "print(sys.version.split()[0],numpy.__version__,scipy.__version__);"
                   "print(circuitpharm.__file__)"], cwd=tree, env=env).split("\n")
        py, np_v, sp_v = ver[0].split()
        print(f"      python {py}  numpy {np_v}  scipy {sp_v}")
        print(f"      circuitpharm loaded from {ver[1]}")
        assert str(tree) in ver[1], (
            f"circuitpharm resolved to {ver[1]}, not the clean tree -- the verification would "
            f"be testing the working tree instead")
        for claimed in (py, np_v, sp_v):
            if claimed not in text:
                print(f"      WARNING: the manuscript does not state version {claimed}")

        print(f"\n[3/4] running the generators from the clean tree")
        out = run([sys.executable, "scripts/paper_numbers.py"], cwd=tree, env=env)
        print(f"      paper_numbers.py            {len(out.splitlines())} lines")
        for floor in (0.0, 0.01, 0.02):
            out += run([sys.executable, "scripts/ranking_robustness.py",
                        "--draws", str(a.draws), "--a5-floor", str(floor)], cwd=tree, env=env)
            print(f"      ranking_robustness.py       --a5-floor {floor}")

        # ---- every quoted figure must appear in that output ----------------------------
        print(f"\n[4/4] checking the manuscript's figures against that output")
        bad = 0
        for cat, vals in CATEGORIES.items():
            misses = [v for v in vals if v not in out]
            absent = [v for v in vals if v not in text]
            status = "ok" if not misses and not absent else "FAIL"
            print(f"      {status:4s} {cat}  ({len(vals)} figures)")
            for v in misses:
                print(f"           {v!r} is in the manuscript but NOT in clean-tree output")
            for v in absent:
                print(f"           {v!r} is expected in the manuscript but absent from it")
            bad += len(misses) + len(absent)

        print()
        if bad:
            print(f"VERDICT: {bad} figure(s) do not reproduce from {ref}.")
            return 1
        total = sum(len(v) for v in CATEGORIES.values())
        print(f"VERDICT: all {total} audited figures reproduce from a clean checkout of {ref}")
        print(f"         across {len(CATEGORIES)} categories, {a.draws} Monte Carlo draws.")
        return 0
    finally:
        if a.keep:
            print(f"\nworktree kept at {tree}")
        else:
            run(["git", "worktree", "remove", "--force", str(tree)])
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
