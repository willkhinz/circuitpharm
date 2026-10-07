"""Reliability tiers, enforced in code rather than in a comment.

WHY THIS MODULE IS THE POINT OF THE PACKAGE. Most in-silico pharmacology tools emit a
point estimate for anything you ask them, and the user cannot tell from the output whether
the number rests on a calibration, on an uncalibrated-but-sound mechanism, or on an anchor
that has since been shown invalid. This project has produced all three, and keeping the
distinction in docstrings and a worklog was not enough -- a stale ventilation percentage
was quoted in project summaries for several sessions while sitting directly beneath an
"UNREACHABLE" flag that the model itself had printed.

So the distinction is a type here. A `Quantity` carries its tier, its provenance, and what
would be needed to promote it. Reading the value of a VOID quantity RAISES unless the
caller explicitly acknowledges it. That is deliberate friction at exactly the point where
this project repeatedly went wrong.

TIER DEFINITIONS, and what it takes to move between them:

  VALIDATED     Reproduces an observable it was NOT fitted to, or is an exact identity, or
                survives propagation of every parameter it depends on. Examples: the
                parity-invariance identity; the selectivity ranking (20k draws); the
                kinetic scheme's phasic predictions (externally corroborated by midazolam
                prolonging IPSC decay without changing amplitude).

  UNCALIBRATED  The mechanism is structurally sound but no valid quantitative anchor
                exists, so the SHAPE and ORDERING may be trusted while the SCALE may not.
                Examples: absolute ventilation percentages; stretch-reflex gains; the
                subjective index (which is algebra over literature weights in arbitrary
                units, not a simulation).
                -> promote by supplying a structure-matched anchor.

  VOID          Known to rest on something invalid. Must not be quoted even with a caveat,
                because caveats get dropped when numbers are copied. Examples: overdose
                multiples (the efficacy ceiling they clip against is pool-dependent and
                the single-pool value is wrong); the NMDA respiratory contribution (an
                established negative result -- this circuit ties burst maintenance to the
                long NMDA conductance, so any block depresses it, where clinically
                ketamine-class blockers preserve or stimulate breathing).
                -> promote by fixing the structural defect, not by refitting.

A tier is a claim about EVIDENCE, never about precision. A number can be reported to four
decimals and still be VOID.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from enum import IntEnum
from typing import Any


class Tier(IntEnum):
    VALIDATED = 1
    UNCALIBRATED = 2
    VOID = 3

    @property
    def label(self) -> str:
        return {1: "VALIDATED", 2: "UNCALIBRATED", 3: "VOID"}[int(self)]


class VoidQuantityError(RuntimeError):
    """Raised when a VOID quantity's value is read without explicit acknowledgement."""


@dataclass(frozen=True)
class Quantity:
    """A model output that knows how much it can be trusted.

    `value` is deliberately NOT a plain attribute read for VOID quantities: use `.get()`
    and, for VOID, pass `acknowledge_void=True`, which forces the caller to write the
    acknowledgement at the call site where a reviewer will see it.
    """
    name: str
    _value: Any
    tier: Tier
    units: str = ""
    provenance: str = ""            # what this rests on
    promote_by: str = ""            # what would raise its tier
    caveats: tuple = field(default_factory=tuple)

    # ---- access ----------------------------------------------------------------------
    def get(self, acknowledge_void: bool = False) -> Any:
        if self.tier is Tier.VOID and not acknowledge_void:
            raise VoidQuantityError(
                f"'{self.name}' is VOID and must not be quoted.\n"
                f"  why:        {self.provenance}\n"
                f"  promote by: {self.promote_by}\n"
                f"  If you genuinely need the raw number for debugging, call "
                f".get(acknowledge_void=True)."
            )
        return self._value

    @property
    def value(self):
        """Convenience for non-VOID quantities; raises for VOID ones. See `get`."""
        return self.get()

    def is_quotable(self) -> bool:
        return self.tier is not Tier.VOID

    # ---- presentation ----------------------------------------------------------------
    def __str__(self) -> str:
        if self.tier is Tier.VOID:
            body = "VOID — not quotable"
        else:
            v = self._value
            body = f"{v:.4g}{(' ' + self.units) if self.units else ''}" \
                if isinstance(v, (int, float)) else str(v)
        return f"{self.name}: {body}  [{self.tier.label}]"

    def report(self) -> str:
        lines = [str(self)]
        if self.provenance:
            lines.append(f"    rests on:   {self.provenance}")
        if self.tier is not Tier.VALIDATED and self.promote_by:
            lines.append(f"    promote by: {self.promote_by}")
        for c in self.caveats:
            lines.append(f"    caveat:     {c}")
        return "\n".join(lines)


@dataclass
class ResultSet:
    """A collection of quantities from one model evaluation.

    Printing a ResultSet never prints a VOID value -- it prints the refusal and the reason.
    That is the behaviour that would have prevented this project's own worst reporting
    error, where a flagged-unreachable configuration was summarised as a clean result.
    """
    title: str
    items: list = field(default_factory=list)

    def add(self, q: Quantity) -> "ResultSet":
        self.items.append(q)
        return self

    def quantity(self, name: str) -> Quantity:
        for q in self.items:
            if q.name == name:
                return q
        raise KeyError(f"no quantity named {name!r}; have "
                       f"{[q.name for q in self.items]}")

    def quotable(self) -> list:
        return [q for q in self.items if q.is_quotable()]

    def by_tier(self, tier: Tier) -> list:
        return [q for q in self.items if q.tier is tier]

    def __str__(self) -> str:
        w = max(72, len(self.title) + 4)
        out = ["=" * w, self.title, "=" * w]
        for tier in (Tier.VALIDATED, Tier.UNCALIBRATED, Tier.VOID):
            group = self.by_tier(tier)
            if not group:
                continue
            out.append(f"\n-- {tier.label} " + "-" * max(0, w - len(tier.label) - 4))
            for q in group:
                out.append(q.report())
        n_void = len(self.by_tier(Tier.VOID))
        if n_void:
            out.append(f"\n{n_void} quantity(ies) withheld as VOID. Quoting them would "
                       f"misrepresent\nthe model; see 'promote by' for what each needs.")
        return "\n".join(out)
