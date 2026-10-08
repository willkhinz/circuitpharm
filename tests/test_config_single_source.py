"""The operating point must have exactly ONE definition.

Recurring error E12 is a calibration constant carried across a change of parameterisation,
and it happened twice in this project. Both times the mechanism was the same: the operating
point and the sensitivities were duplicated as literals in several scripts, so updating one
left the others silently stale.

These tests pin the single source. They are cheap and they fail loudly the moment someone
pastes a literal back in.
"""
import pathlib
import re

import pytest

from circuitpharm.config import (RESP_OP, SYNAPTIC_PULSE, AMBIENT_GABA_UM,
                                 EUPNOEA_BAND, IA_SCALE, CALIBRATIONS)

REPO = pathlib.Path(__file__).resolve().parents[1]
LIVE_SCRIPTS = sorted(p for p in (REPO / "scripts").glob("*.py"))


def test_config_values_are_what_the_project_calibrated():
    """Guards against a silent edit to the shared constants themselves."""
    assert RESP_OP["drive"] == 170.0
    assert RESP_OP["g_adapt"] == 2.5
    assert RESP_OP["tau_adapt"] == 400.0
    assert RESP_OP["w"]["ee_ampa"] == 0.45
    assert RESP_OP["w"]["ee_nmda"] == 0.2475
    assert SYNAPTIC_PULSE == {"peak_um": 3000.0, "clear_ms": 1.00}
    assert AMBIENT_GABA_UM == 0.40
    assert IA_SCALE == 0.30
    assert EUPNOEA_BAND == (0.30, 2.50)


def test_resp_op_is_immutable():
    """It is shared, so a caller mutating it would poison every later run in the process."""
    with pytest.raises(TypeError):
        RESP_OP["drive"] = 999.0
    with pytest.raises(TypeError):
        RESP_OP["w"]["ee_ampa"] = 999.0


@pytest.mark.parametrize("path", LIVE_SCRIPTS, ids=lambda p: p.name)
def test_no_live_script_redefines_the_operating_point(path):
    """No LIVE script may contain the operating point as a literal. Archived scripts are
    exempt: they are frozen records of superseded runs and must keep the values they ran
    with."""
    src = path.read_text()
    assert "drive=170.0" not in src, (
        f"{path.name} hard-codes the preBotC drive. Import RESP_OP from "
        "circuitpharm.config instead -- duplicated constants are how E12 happened twice.")
    assert not re.search(r"peak_um\s*=\s*3000", src), (
        f"{path.name} hard-codes the cleft pulse; import SYNAPTIC_PULSE.")


def test_the_void_calibration_cannot_be_quietly_promoted():
    """`prebotc_gaba_sens` is VOID because whole-body ventilation is the wrong observable
    for an isolated preBotC. Promoting it requires a structure-matched anchor, not an edit,
    so this test exists to make a quiet promotion fail."""
    c = CALIBRATIONS["prebotc_gaba_sens"]
    from circuitpharm.results import Tier
    assert c.tier is Tier.VOID
    assert "wrong observable" in c.status.lower()
