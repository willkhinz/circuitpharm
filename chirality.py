#!/usr/bin/env python3
"""Chirality toolkit: stereocentre analysis, the parity-invariance validation harness, and
a demonstration of why standard QSAR descriptors cannot distinguish enantiomers.

    python chirality.py stereo          # stereocentres in the project's compounds
    python chirality.py parity          # the parity-invariance unit test
    python chirality.py descriptors     # why achiral descriptors fail on enantiomers

THE PARITY-INVARIANCE HARNESS. The molecular Hamiltonian is parity-invariant (weak-force
parity violation between enantiomers is ~1e-17 in relative energy, far below thermal
noise). Therefore, EXACTLY:

    Mirror(ligand + receptor)  ==  (mirrored ligand) + (mirrored receptor)

with identical energy. This is an identity, not an approximation, so it is a free and exact
unit test for any structure-based scoring pipeline: build the mirrored receptor by
reflecting coordinates, dock the mirrored ligand into it, and the score MUST match the
original to numerical precision. Any difference means the engine, the force field
parameters, or the setup is broken. Equivalently, since the receptor is built from
L-amino acids:

    (D-ligand + L-receptor)  ==  Mirror of  (L-ligand + D-receptor)

which is why a chiral binding site discriminates enantiomers at all, and why the
discrimination is purely a property of the complex rather than of either partner alone.

WHY THE DESCRIPTOR TEST MATTERS. Standard QSAR descriptors — molecular weight, logP, TPSA,
H-bond donors/acceptors, ring counts — are computed from the molecular graph and are
therefore IDENTICAL for enantiomers. A fingerprint model is structurally incapable of
distinguishing a eutomer from a distomer, which is why chiral pharmacology needs
mechanistic or 3D methods rather than 2D ML. This script proves it numerically.
"""
import argparse, sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import rdMolTransforms  # noqa: F401
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

# project-relevant molecules
MOLS = {
    "alogabat (RG7816) — lead a5 PAM":
        "CC1=NC=C(C=C1)C2=NOC(=C2COC3=NN=C(C=C3)C(=O)NC4CCOCC4)C",
    "ketamine (racemic)":            "CNC1(c2ccccc2Cl)CCCCC1=O",
    "methadone (racemic)":           "CCC(=O)C(CC(C)N(C)C)(c1ccccc1)c1ccccc1",
    "baclofen (racemic)":            "NCC(CC(=O)O)c1ccc(Cl)cc1",
    "etomidate":                     "CCOC(=O)c1cncn1C(C)c1ccccc1",
    "allopregnanolone":              "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4[C@@H]3CC[C@H](O)C4)[C@@H]1CC[C@@H]2C(C)=O",
    "levetiracetam (S)":             "CC[C@H](N1CCCC1=O)C(N)=O",
    "salsolinol":                    "CC1NCCc2cc(O)c(O)cc21",
}


def stereo_report(name, smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return f"{name:<42} INVALID SMILES"
    Chem.AssignStereochemistry(m, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    centres = Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
    si = Chem.FindPotentialStereo(m)
    n_at = sum(1 for e in si if str(e.type) == "Atom_Tetrahedral")
    n_bd = sum(1 for e in si if str(e.type) == "Bond_Double")
    iso = len(tuple(Chem.EnumerateStereoisomers.EnumerateStereoisomers(m))) \
        if hasattr(Chem, "EnumerateStereoisomers") else None
    tag = "ACHIRAL — no stereocentres" if n_at == 0 else f"{n_at} tetrahedral stereocentre(s)"
    return (f"{name:<42} {rdMolDescriptors.CalcMolFormula(m):<17}"
            f"MW {Descriptors.MolWt(m):7.2f}  {tag}"
            + (f", {n_bd} stereo double bond(s)" if n_bd else "")
            + (f"\n{'':<42} centres: {centres}" if centres else ""))


def cmd_stereo():
    print("STEREOCENTRE ANALYSIS\n" + "=" * 100)
    for name, smi in MOLS.items():
        print(stereo_report(name, smi))
    print("\nThe lead compound is the question that was open: see whether alogabat has any.")


def reflect(conf_pos):
    """Improper rotation: negate x. Any reflection works; this is the simplest."""
    p = conf_pos.copy()
    p[:, 0] *= -1.0
    return p


def internal_geometry(pos):
    """All pairwise distances — invariant under ANY proper or improper isometry."""
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    return d[np.triu_indices_from(d, k=1)]


def cmd_parity():
    print("PARITY-INVARIANCE UNIT TEST\n" + "=" * 78)
    print("Claim: reflection is an exact isometry, so a mirrored molecule has identical")
    print("internal geometry and identical energy in an ACHIRAL environment. Only the")
    print("CIP label and any interaction with a CHIRAL partner change.\n")
    failures = 0
    for name, smi in (("ketamine", MOLS["ketamine (racemic)"]),
                      ("levetiracetam (S)", MOLS["levetiracetam (S)"]),
                      ("allopregnanolone", MOLS["allopregnanolone"])):
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        if AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE) != 0:
            print(f"  {name}: embedding failed"); continue
        AllChem.MMFFOptimizeMolecule(m)
        pos = m.GetConformer().GetPositions()

        mirror = Chem.Mol(m)
        c = mirror.GetConformer()
        mp = reflect(pos)
        for i in range(mirror.GetNumAtoms()):
            c.SetAtomPosition(i, mp[i].tolist())

        # 1. internal geometry must be EXACTLY preserved
        d0, d1 = internal_geometry(pos), internal_geometry(mp)
        geom_err = float(np.abs(d0 - d1).max())

        # 2. MMFF energy in vacuum (an achiral environment) must be identical
        def mmff(mol):
            props = AllChem.MMFFGetMoleculeProperties(mol)
            ff = AllChem.MMFFGetMoleculeForceField(mol, props)
            return ff.CalcEnergy()
        e0, e1 = mmff(m), mmff(mirror)
        e_err = abs(e0 - e1)

        # 3. the CIP label MUST flip
        def labels(mol):
            mm = Chem.Mol(mol)
            Chem.AssignStereochemistryFrom3D(mm)
            return tuple(a.GetPropsAsDict().get("_CIPCode", "")
                         for a in mm.GetAtoms() if a.HasProp("_CIPCode"))
        l0, l1 = labels(m), labels(mirror)
        flipped = bool(l0) and all(x != y for x, y in zip(l0, l1))

        ok_geom = geom_err < 1e-9
        ok_e = e_err < 1e-6
        failures += 0 if (ok_geom and ok_e and flipped) else 1
        print(f"  {name}")
        print(f"    internal geometry   max |Δd| = {geom_err:.2e}   "
              f"{'PASS' if ok_geom else 'FAIL'}")
        print(f"    MMFF energy vacuum  |ΔE| = {e_err:.2e} kcal/mol  "
              f"{'PASS' if ok_e else 'FAIL'}")
        print(f"    CIP labels          {l0} -> {l1}   "
              f"{'PASS (flipped)' if flipped else 'FAIL (did not flip)'}")
    print(f"\n  {'ALL PASS' if failures == 0 else f'{failures} FAILURE(S)'}")
    print("""
HOW TO USE THIS AS A DOCKING HARNESS
  1. take the receptor structure and the ligand pose
  2. reflect BOTH sets of coordinates through any plane
  3. re-score
  4. the score must be identical to numerical precision
A difference means the engine, the force-field parameters, or the setup is broken --
most commonly chirality-naive parameters, or a scoring term that is not parity-even.
This test costs nothing and is exact, so run it before trusting any chiral docking result.""")


def cmd_descriptors():
    print("WHY STANDARD QSAR DESCRIPTORS CANNOT DISTINGUISH ENANTIOMERS\n" + "=" * 78)
    pairs = [("levetiracetam (S)", "CC[C@H](N1CCCC1=O)C(N)=O",
              "R-etiracetam (~1000x less active)", "CC[C@@H](N1CCCC1=O)C(N)=O"),
             ("R-baclofen (arbaclofen, active)", "NC[C@@H](CC(=O)O)c1ccc(Cl)cc1",
              "S-baclofen (~100x weaker)", "NC[C@H](CC(=O)O)c1ccc(Cl)cc1")]
    desc = [("MolWt", Descriptors.MolWt), ("MolLogP", Descriptors.MolLogP),
            ("TPSA", Descriptors.TPSA), ("HBD", Descriptors.NumHDonors),
            ("HBA", Descriptors.NumHAcceptors),
            ("RotB", Descriptors.NumRotatableBonds),
            ("RingCount", Descriptors.RingCount),
            ("HeavyAtoms", Descriptors.HeavyAtomCount)]
    for n1, s1, n2, s2 in pairs:
        a, b = Chem.MolFromSmiles(s1), Chem.MolFromSmiles(s2)
        print(f"\n  {n1}\n  vs {n2}")
        print(f"    {'descriptor':<14}{'eutomer':>12}{'distomer':>12}{'Δ':>10}")
        for nm, fn in desc:
            x, y = fn(a), fn(b)
            print(f"    {nm:<14}{x:12.4f}{y:12.4f}{x-y:10.4f}")
        fp1 = AllChem.GetMorganFingerprintAsBitVect(a, 2, 2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(b, 2, 2048)
        same = list(fp1.GetOnBits()) == list(fp2.GetOnBits())
        print(f"    Morgan fingerprint (radius 2, achiral): "
              f"{'IDENTICAL' if same else 'different'}")
        fp1c = AllChem.GetMorganFingerprintAsBitVect(a, 2, 2048, useChirality=True)
        fp2c = AllChem.GetMorganFingerprintAsBitVect(b, 2, 2048, useChirality=True)
        samec = list(fp1c.GetOnBits()) == list(fp2c.GetOnBits())
        print(f"    Morgan fingerprint with useChirality=True: "
              f"{'IDENTICAL' if samec else 'different'}")
    print("""
CONSEQUENCE. Every scalar descriptor above is identical, and the default fingerprint is
identical, for compounds differing up to 1000-fold in activity. A 2D QSAR or fingerprint
model is therefore STRUCTURALLY INCAPABLE of ranking enantiomers -- not merely untrained
for it. Chirality-aware fingerprints exist (useChirality=True) and do differ, but they
encode only the LABEL, carrying no information about which label a chiral site prefers.
That is why chiral pharmacology needs mechanistic or 3D structure-based methods, and why
this project's receptor-profile simulator takes measured subtype activities as INPUT rather
than trying to predict them from structure.""")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["stereo", "parity", "descriptors"])
    a = ap.parse_args()
    {"stereo": cmd_stereo, "parity": cmd_parity, "descriptors": cmd_descriptors}[a.cmd]()
