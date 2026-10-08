"""Shared fixtures and optional-dependency guards.

WHY THIS MATTERS FOR CI. The MuJoCo body plant is an optional extra (`.[plant]`), so a lean
install -- which is what CI uses, and what most users will have -- cannot run the reflex or
locomotion assays. Without a guard those tests ERROR rather than skip, which turns a
perfectly healthy lean install into a red build and trains people to ignore CI.

`requires_plant` skips them cleanly instead. The same applies to RDKit for the chirality
tests, which handle it with importorskip at module level.
"""
import pytest


def _have(mod: str) -> bool:
    try:
        __import__(mod)
        return True
    except Exception:
        return False


HAVE_MUJOCO = _have("mujoco")
HAVE_RDKIT = _have("rdkit")

requires_plant = pytest.mark.skipif(
    not HAVE_MUJOCO,
    reason="needs the optional body plant: pip install -e '.[plant]'")

requires_chem = pytest.mark.skipif(
    not HAVE_RDKIT,
    reason="needs the optional chemistry tools: pip install -e '.[chem]'")


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "needs_plant: requires the optional MuJoCo body plant")


def pytest_collection_modifyitems(config, items):
    """Auto-skip anything that touches the body plant when MuJoCo is absent.

    Done by inspection rather than by hand-marking every test, so a newly added
    plant-dependent test cannot forget its marker and break the lean build.
    """
    if HAVE_MUJOCO:
        return
    skip = pytest.mark.skip(reason="MuJoCo not installed (optional 'plant' extra)")
    for item in items:
        src = ""
        try:
            src = item.function.__code__.co_consts and str(item.function.__code__.co_names)
        except Exception:
            pass
        names = set(getattr(item.function, "__code__", None).co_names) \
            if getattr(item.function, "__code__", None) else set()
        # Two signals only, both specific:
        #   1. an explicit @pytest.mark.needs_plant
        #   2. the test actually NAMES a plant-dependent callable
        #
        # The NAME of the test is deliberately NOT used any more. It was -- any test whose
        # name contained "reflex" or "walk" got skipped -- which silently skips a test of,
        # say, reflex arithmetic that needs no MuJoCo at all. A skipped test that should
        # have run is worse than a failing one, because nothing reports it.
        #
        # The call-name check is kept as a safety net for a new test that forgets the
        # marker, but it cannot catch everything: a test can depend on the plant without
        # naming an assay (one that calls evaluate() and asserts the motor keys exist).
        # That case was found by running the suite against a clean install. Hence the
        # explicit marker is the primary mechanism.
        PLANT_CALLS = {"stretch_reflex", "locomotion", "JointPlant", "model_xml",
                       "REFLEX_CASES"}
        if item.get_closest_marker("needs_plant") or (names & PLANT_CALLS):
            item.add_marker(skip)
