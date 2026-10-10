"""Machinery that runs, and whose conclusions rest on a recorded defect.

WHY THIS IS NOT AN EXCEPTION. A FAILURE is a computation that did not happen; it raises
(roadmap §2.4, and see the three `RuntimeError`s added to `models/` in P0-3). This module
is for the other case: code that computes exactly what it says it computes, on an input or
an objective that cannot support the inference drawn from it. The arrays are real and
someone may legitimately want them; the conclusion is not available.

Two ways of handling that were considered and rejected:

  * `NotImplementedError`. Looks rigorous, destroys a real computation, and hides the
    defect behind an exception nobody reads. Recorded as forbidden pattern §5.4a.
  * A docstring note. Does not survive a copy-paste into a slide, which is the exact
    failure this project has already made: a stale ventilation percentage was quoted for
    several sessions while sitting directly beneath an UNREACHABLE flag the model itself
    had printed (see `results.py`).

So: the function stays callable, it warns on EVERY call, and its inferential summaries come
back as `Tier.VOID` quantities whose `.value` raises. The raw computed arrays stay plain.
"""
from __future__ import annotations

import warnings

__all__ = ["ProvisionalResultWarning", "warn_provisional", "void"]


class ProvisionalResultWarning(UserWarning):
    """Emitted by machinery whose conclusions rest on a defect in the P0 register."""


# "always", NOT the default "once per location".
#
# Python's default filter for UserWarning shows a given warning once per unique code
# location per process. That means the person who makes the first call sees the banner and
# every later caller in the same session does not -- including the loop that generates the
# table that goes into the write-up. That is precisely how a caveat gets dropped between
# the run and the paper, so the filter is widened for this category only.
#
# `tests/test_p0_register.py` and `tests/test_provisional.py` pin this: a test calls twice
# inside catch_warnings and asserts TWO records. Without it, a future tidy-up silently
# restores once-per-location and nothing fails.
warnings.simplefilter("always", ProvisionalResultWarning)


def warn_provisional(what: str, defect: str, register_item: str, promote_by: str,
                     stacklevel: int = 3) -> None:
    """Emit the standard provisional-result banner.

    Every argument is mandatory and must be specific: a banner that says "results may be
    unreliable" is worth nothing. `defect` says what is wrong, `register_item` says where
    it is written down, and `promote_by` says what would make the conclusion usable.
    """
    for name, value in (("what", what), ("defect", defect),
                        ("register_item", register_item), ("promote_by", promote_by)):
        if not isinstance(value, str) or not value.strip():
            raise ValueError(
                f"warn_provisional() needs a specific {name!r}; a vague banner is worse "
                f"than none because it trains people to ignore the category.")

    warnings.warn(
        "\n" + "=" * 78 + "\n"
        f"{what} IS PROVISIONAL -- its conclusions are not usable.\n"
        + "=" * 78 + "\n"
        f"  defect:     {defect}\n"
        f"  recorded:   roadmap {register_item}\n"
        f"  promote by: {promote_by}\n"
        "  The computation below is real; the inference from it is not. Reading a VOID\n"
        "  field of the result raises unless you pass acknowledge_void=True.\n"
        + "=" * 78,
        ProvisionalResultWarning, stacklevel=stacklevel)


def void(name: str, value, *, defect: str, register_item: str, promote_by: str,
         units: str = "", caveats: tuple = ()):
    """Wrap an unsupported inference as a `Tier.VOID` `Quantity`.

    Imported lazily so this module stays dependency-free for `models/`.
    """
    from .results import Quantity, Tier

    return Quantity(
        name=name, _value=value, tier=Tier.VOID, units=units,
        provenance=f"{defect} (roadmap {register_item})",
        promote_by=promote_by, caveats=caveats)
