"""Stereochemistry: the parity-invariance identity, and why 2D descriptors cannot rank
enantiomers.

These are the package's only EXACT results. The molecular Hamiltonian is parity-invariant
(weak-force parity violation between enantiomers is ~1e-17 in relative energy, far below
thermal noise), so reflecting a molecule must leave its internal geometry and its energy in
an achiral environment UNCHANGED -- not approximately, exactly. That makes it a free and
unarguable unit test for any structure-based scoring pipeline: reflect both receptor and
ligand, re-score, and the score must match to numerical precision. A difference means the
engine, the force-field parameters, or the setup is broken.

RDKit is an optional extra (`pip install -e ".[chem]"`), so these skip when it is absent.
"""
import pytest

rdkit = pytest.importorskip("rdkit", reason="RDKit is the optional 'chem' extra")

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

from circuitpharm.chirality import reflect, internal_geometry, MOLS


def _embed(smi, seed=0xC0FFEE):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    assert AllChem.EmbedMolecule(m, randomSeed=seed) == 0, "embedding failed"
    AllChem.MMFFOptimizeMolecule(m)
    return m


def _mirror(m):
    mm = Chem.Mol(m)
    c = mm.GetConformer()
    mp = reflect(m.GetConformer().GetPositions())
    for i in range(mm.GetNumAtoms()):
        c.SetAtomPosition(i, mp[i].tolist())
    return mm, mp


def _mmff(mol):
    props = AllChem.MMFFGetMoleculeProperties(mol)
    return AllChem.MMFFGetMoleculeForceField(mol, props).CalcEnergy()


CHIRAL = ["ketamine (racemic)", "levetiracetam (S)", "allopregnanolone"]


@pytest.mark.parametrize("name", CHIRAL)
def test_reflection_preserves_internal_geometry_exactly(name):
    """Every pairwise distance is invariant under ANY isometry, proper or improper.
    Tolerance is float round-off, not 'close enough'."""
    m = _embed(MOLS[name])
    pos = m.GetConformer().GetPositions()
    _, mp = _mirror(m)
    err = float(np.abs(internal_geometry(pos) - internal_geometry(mp)).max())
    assert err < 1e-9, f"reflection changed internal geometry by {err:.2e} A"


@pytest.mark.parametrize("name", CHIRAL)
def test_reflection_preserves_mmff_energy_exactly(name):
    """A mirrored molecule must have identical energy in vacuum, which is achiral.
    If this fails, the force field is not parity-even and no chiral docking result from
    it can be trusted."""
    m = _embed(MOLS[name])
    mm, _ = _mirror(m)
    err = abs(_mmff(m) - _mmff(mm))
    assert err < 1e-6, f"reflection changed MMFF energy by {err:.2e} kcal/mol"


@pytest.mark.parametrize("name", CHIRAL)
def test_reflection_flips_every_cip_label(name):
    """The geometry and energy are invariant, but the LABELS must all invert -- that is
    what makes the mirrored molecule a different compound."""
    m = _embed(MOLS[name])
    mm, _ = _mirror(m)

    def labels(mol):
        x = Chem.Mol(mol)
        Chem.AssignStereochemistryFrom3D(x)
        return tuple(a.GetPropsAsDict().get("_CIPCode", "")
                     for a in x.GetAtoms() if a.HasProp("_CIPCode"))

    l0, l1 = labels(m), labels(mm)
    assert l0, f"{name} has no assigned stereocentres to flip"
    assert all(a != b for a, b in zip(l0, l1)), f"{l0} -> {l1} did not fully invert"


def test_allopregnanolone_has_eight_stereocentres_all_flipping():
    """A strong case: eight centres must invert together, not individually."""
    m = _embed(MOLS["allopregnanolone"])
    mm, _ = _mirror(m)

    def labels(mol):
        x = Chem.Mol(mol)
        Chem.AssignStereochemistryFrom3D(x)
        return tuple(a.GetPropsAsDict().get("_CIPCode", "")
                     for a in x.GetAtoms() if a.HasProp("_CIPCode"))

    l0, l1 = labels(m), labels(mm)
    assert len(l0) == 8, f"expected 8 stereocentres, found {len(l0)}"
    assert all(a != b for a, b in zip(l0, l1))


# ------------------------------------------------- why 2D descriptors cannot rank them
PAIRS = [("levetiracetam (S)", "CC[C@H](N1CCCC1=O)C(N)=O",
          "R-etiracetam", "CC[C@@H](N1CCCC1=O)C(N)=O"),
         ("R-baclofen (arbaclofen)", "NC[C@@H](CC(=O)O)c1ccc(Cl)cc1",
          "S-baclofen", "NC[C@H](CC(=O)O)c1ccc(Cl)cc1")]


@pytest.mark.parametrize("n1,s1,n2,s2", PAIRS)
def test_scalar_descriptors_are_identical_for_enantiomers(n1, s1, n2, s2):
    """These pairs differ up to 1000-fold in activity. Every scalar descriptor is
    identical, so a 2D QSAR model is STRUCTURALLY incapable of ranking them -- not merely
    untrained for it. This is why the package takes measured subtype activities as INPUT
    instead of predicting them from structure."""
    a, b = Chem.MolFromSmiles(s1), Chem.MolFromSmiles(s2)
    for fn in (Descriptors.MolWt, Descriptors.MolLogP, Descriptors.TPSA,
               Descriptors.NumHDonors, Descriptors.NumHAcceptors,
               Descriptors.NumRotatableBonds, Descriptors.RingCount,
               Descriptors.HeavyAtomCount):
        assert fn(a) == pytest.approx(fn(b), abs=1e-9), fn.__name__


@pytest.mark.parametrize("n1,s1,n2,s2", PAIRS)
def test_default_morgan_fingerprint_is_bit_identical(n1, s1, n2, s2):
    a, b = Chem.MolFromSmiles(s1), Chem.MolFromSmiles(s2)
    fa = AllChem.GetMorganFingerprintAsBitVect(a, 2, 2048)
    fb = AllChem.GetMorganFingerprintAsBitVect(b, 2, 2048)
    assert list(fa.GetOnBits()) == list(fb.GetOnBits())


@pytest.mark.parametrize("n1,s1,n2,s2", PAIRS)
def test_chirality_aware_fingerprint_differs_but_only_encodes_the_label(n1, s1, n2, s2):
    """useChirality=True does distinguish them -- but it encodes only WHICH label, carrying
    no information about which label a chiral binding site prefers. Distinguishable is not
    rankable."""
    a, b = Chem.MolFromSmiles(s1), Chem.MolFromSmiles(s2)
    fa = AllChem.GetMorganFingerprintAsBitVect(a, 2, 2048, useChirality=True)
    fb = AllChem.GetMorganFingerprintAsBitVect(b, 2, 2048, useChirality=True)
    assert list(fa.GetOnBits()) != list(fb.GetOnBits())


def test_alogabat_is_achiral():
    """Recorded finding: the lead a5-selective compound has NO stereocentres, so a5
    selectivity does not require a stereocentre (SH-053 achieves it chirally, alogabat
    achirally). Pinned because the ranking separates the SH-053 enantiomers and it would
    be easy to assume chirality is load-bearing for the mechanism."""
    m = Chem.MolFromSmiles(MOLS["alogabat (RG7816) — lead a5 PAM"])
    assert m is not None
    si = Chem.FindPotentialStereo(m)
    n_tet = sum(1 for e in si if str(e.type) == "Atom_Tetrahedral")
    assert n_tet == 0, f"alogabat reported {n_tet} tetrahedral stereocentres"
