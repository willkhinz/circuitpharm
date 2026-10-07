"""The reliability-tier machinery must actually refuse, not merely annotate.

This exists because the project's worst reporting error was not a modelling mistake: a
configuration the model itself had flagged UNREACHABLE was summarised as a clean result,
because the flag was printed next to the number instead of blocking it. Caveats get dropped
when numbers are copied. These tests pin the refusal behaviour.
"""
import pytest

from circuitpharm.results import Tier, Quantity, ResultSet, VoidQuantityError


def _q(tier, value=42.0, name="x"):
    return Quantity(name=name, _value=value, tier=tier, units="%",
                    provenance="test", promote_by="do the experiment")


def test_validated_and_uncalibrated_values_are_readable():
    assert _q(Tier.VALIDATED).value == 42.0
    assert _q(Tier.UNCALIBRATED).value == 42.0
    assert _q(Tier.VALIDATED).get() == 42.0


def test_void_value_raises_on_plain_access():
    q = _q(Tier.VOID)
    with pytest.raises(VoidQuantityError):
        q.value
    with pytest.raises(VoidQuantityError):
        q.get()


def test_void_raise_message_carries_the_reason_and_the_fix():
    q = Quantity(name="overdose_index", _value=26.0, tier=Tier.VOID,
                 provenance="clips against a pool-dependent ceiling",
                 promote_by="fix the structural defect, do not refit")
    with pytest.raises(VoidQuantityError) as e:
        q.value
    msg = str(e.value)
    assert "overdose_index" in msg
    assert "pool-dependent" in msg          # why it is void
    assert "do not refit" in msg            # how to promote it


def test_void_is_readable_only_with_explicit_acknowledgement():
    """Debugging access exists, but the acknowledgement must be written at the call site
    where a reviewer will see it."""
    assert _q(Tier.VOID).get(acknowledge_void=True) == 42.0


def test_is_quotable_matches_tier():
    assert _q(Tier.VALIDATED).is_quotable()
    assert _q(Tier.UNCALIBRATED).is_quotable()
    assert not _q(Tier.VOID).is_quotable()


def test_void_str_never_leaks_the_number():
    """Printing must not emit the value, or the refusal is cosmetic."""
    s = str(_q(Tier.VOID, value=123.456))
    assert "123" not in s
    assert "VOID" in s


def test_resultset_withholds_void_values_when_printed():
    rs = (ResultSet("demo")
          .add(_q(Tier.VALIDATED, 1.0, "ranking"))
          .add(_q(Tier.UNCALIBRATED, 2.0, "ventilation"))
          .add(_q(Tier.VOID, 999.0, "overdose")))
    out = str(rs)
    assert "999" not in out, "a VOID value leaked into the printed report"
    assert "withheld" in out
    assert "ranking" in out and "ventilation" in out


def test_resultset_quotable_excludes_void():
    rs = (ResultSet("demo")
          .add(_q(Tier.VALIDATED, 1.0, "a"))
          .add(_q(Tier.VOID, 2.0, "b")))
    assert [q.name for q in rs.quotable()] == ["a"]


def test_resultset_lookup_and_missing_key():
    rs = ResultSet("demo").add(_q(Tier.VALIDATED, 1.0, "a"))
    assert rs.quantity("a").value == 1.0
    with pytest.raises(KeyError):
        rs.quantity("nope")


def test_tier_ordering_is_meaningful():
    assert Tier.VALIDATED < Tier.UNCALIBRATED < Tier.VOID


def test_quantity_is_immutable():
    """A tier must not be downgradeable in passing by downstream code."""
    q = _q(Tier.VOID)
    with pytest.raises(Exception):
        q.tier = Tier.VALIDATED


def test_calibration_report_renders_and_states_the_gap():
    from circuitpharm.config import calibration_report
    txt = calibration_report()
    assert "VOID" in txt and "VALIDATED" in txt
    assert "NO valid quantitative anchor" in txt
