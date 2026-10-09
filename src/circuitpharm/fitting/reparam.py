"""The combinations the data can determine, and the map back to microscopic rates.

THE STRUCTURAL RESULT, derived rather than discovered numerically.

The 5-state scheme's EQUILIBRIUM state distribution is, with x = [G]/K_d,

    P_R : P_AR : P_A2R : P_A2O : P_A2D  =  1 : 2x : x^2 : E x^2 : D x^2

so every equilibrium observable is a function of exactly three quantities:

    K_d = k_off / k_on        the microscopic dissociation constant
    E   = beta / alpha        the gating efficacy
    D   = d / r               the desensitisation ratio

The six microscopic rates enter ONLY through those three ratios. Equivalently, the
equilibrium likelihood is invariant under a three-parameter scaling group -- multiply any
of the pairs (k_on, k_off), (beta, alpha) or (d, r) by a common factor and nothing changes.
`fitting.identifiability.invariance_report` measures this: the objective is constant to
2e-12 along each generator, with normalised curvature of order 1e-12 to 1e-11.

SO AN EQUILIBRIUM-ONLY LIKELIHOOD DETERMINES THREE NUMBERS, NOT SIX, and any per-rate
confidence interval from one is a statement about the prior box.

WHAT BREAKS THE DEGENERACY. A KINETIC observable -- anything with a timescale in it --
does not share the invariance, because scaling a pair changes how FAST the system
equilibrates while leaving where it equilibrates alone. One deactivation time course adds
one absolute timescale, which pins one more combination:

    log10 k_off        from a deactivation or relaxation time constant

With K_d known, k_on follows. beta and alpha separate only with an observable that resolves
GATING specifically -- a single-channel mean open time (1/alpha), a burst analysis, or a
rise time at saturating agonist -- and d and r separate only with a desensitisation onset
or recovery time course. So the ladder is:

    equilibrium CRC alone          -> K_d, E, D                      (3 of 6)
    + a deactivation time course   -> + k_off, hence k_on            (5 of 6, as 4 free)
    + a mean open time             -> + alpha, hence beta            (6 of 6)
    + a desensitisation time course-> + d, hence r                   (fully determined)

`fitting.data.MISSING_DATASETS` records which of those this repository has: none of them.
That is why `gabaa_kinetics.FIT_FIXED_ALPHA` is a CONVENTION, and why fixing it is the
single most valuable dataset acquisition in the project (roadmap P1-2).

NOTE THE SCOPE. This is about the EQUILIBRIUM likelihood in `identifiability.py`.
`gabaa_kinetics.fit_scheme` fits PEAK and IPSC-decay observables, which do carry
timescales, so it is not degenerate in the same way -- it was degenerate for the simpler
reason of having four unknowns and three residuals (roadmap P0-7).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: The three combinations an equilibrium observable can determine, in log10 space.
EQUILIBRIUM_IDENTIFIABLE = ("log10_kd", "log10_E", "log10_D")

#: What each additional observable adds. Keys are the combination that becomes
#: identifiable; values say which measurement supplies it.
UNLOCKED_BY = {
    "log10_koff": ("a deactivation or relaxation TIME COURSE. Adds one absolute "
                   "timescale; with K_d already known this gives k_on too."),
    "log10_alpha": ("a GATING-resolving observable: single-channel mean open time "
                    "(1/alpha), a burst analysis, or a rise time at saturating agonist. "
                    "With E known this separates beta."),
    "log10_d": ("a desensitisation onset or recovery time course. With D known this "
                "separates r."),
}


@dataclass(frozen=True)
class IdentifiableParams:
    """The equilibrium-identifiable combinations, with an optional absolute timescale.

    `log10_koff` is `None` when no kinetic observable is in the likelihood, and in that
    case `to_microscopic` cannot produce rates without being told one -- which is the
    point: the map is many-to-one and the code refuses to pretend otherwise.
    """

    log10_kd: float
    log10_E: float
    log10_D: float
    log10_koff: float | None = None

    @property
    def kd_um(self) -> float:
        return float(10.0 ** self.log10_kd)

    @property
    def gating_efficacy(self) -> float:
        return float(10.0 ** self.log10_E)

    @property
    def desens_ratio(self) -> float:
        return float(10.0 ** self.log10_D)

    # ---- the map -----------------------------------------------------------------
    @classmethod
    def from_microscopic(cls, kon: float, koff: float, beta: float, alpha: float,
                         d: float, r: float, *, with_timescale: bool = True
                         ) -> "IdentifiableParams":
        """Project microscopic rates onto the identifiable combinations.

        MANY-TO-ONE BY CONSTRUCTION. Any (kon, koff) with the same ratio maps here
        identically, and so for the other two pairs -- that is the whole content of the
        structural result. `with_timescale=False` drops `log10_koff` as well, giving the
        strictly equilibrium-identifiable triple.
        """
        for name, v in (("kon", kon), ("koff", koff), ("beta", beta), ("alpha", alpha),
                        ("d", d), ("r", r)):
            if not (np.isfinite(v) and v > 0):
                raise ValueError(f"{name} must be finite and positive, got {v!r}")
        return cls(
            log10_kd=float(np.log10(koff / kon)),
            log10_E=float(np.log10(beta / alpha)),
            log10_D=float(np.log10(d / r)),
            log10_koff=float(np.log10(koff)) if with_timescale else None,
        )

    def to_microscopic(self, *, koff: float | None = None, alpha: float | None = None,
                       r: float | None = None) -> dict[str, float]:
        """Recover microscopic rates, given whatever the data could not determine.

        `koff` may be omitted if this object carries `log10_koff`. `alpha` and `r` must
        ALWAYS be supplied, because no equilibrium or single-timescale observable separates
        them from `beta` and `d` -- see `UNLOCKED_BY`. Supplying them is a declared
        convention, exactly as `gabaa_kinetics.FIT_FIXED_ALPHA` is, and the caller has to
        write it down rather than have a default quietly chosen here.
        """
        if koff is None:
            if self.log10_koff is None:
                raise ValueError(
                    "this IdentifiableParams carries no absolute timescale (log10_koff is "
                    "None), so microscopic rates cannot be recovered. Either fit against a "
                    "dataset with a timescale in it -- see fitting.data.MISSING_DATASETS "
                    "-- or pass koff= explicitly as a declared convention.")
            koff = float(10.0 ** self.log10_koff)
        if alpha is None or r is None:
            raise ValueError(
                "alpha and r must be supplied: no equilibrium observable, and no single "
                "deactivation timescale, separates alpha from beta or r from d. See "
                "reparam.UNLOCKED_BY for which measurement would. Passing them is a "
                "convention and must be recorded as one (provenance.Basis.CONVENTION).")
        for name, v in (("koff", koff), ("alpha", alpha), ("r", r)):
            if not (np.isfinite(v) and v > 0):
                raise ValueError(f"{name} must be finite and positive, got {v!r}")

        return dict(
            kon=float(koff / self.kd_um),
            koff=float(koff),
            beta=float(alpha * self.gating_efficacy),
            alpha=float(alpha),
            d=float(r * self.desens_ratio),
            r=float(r),
        )

    def as_vector(self, *, with_timescale: bool | None = None) -> np.ndarray:
        use = self.log10_koff is not None if with_timescale is None else with_timescale
        if use and self.log10_koff is None:
            raise ValueError("no timescale to include")
        base = [self.log10_kd, self.log10_E, self.log10_D]
        return np.array(base + ([self.log10_koff] if use else []), dtype=float)

    @classmethod
    def from_vector(cls, v) -> "IdentifiableParams":
        v = np.asarray(v, dtype=float)
        if v.size == 3:
            return cls(*map(float, v))
        if v.size == 4:
            return cls(float(v[0]), float(v[1]), float(v[2]), float(v[3]))
        raise ValueError(f"expected 3 or 4 entries, got {v.size}")

    @property
    def names(self) -> tuple[str, ...]:
        base = EQUILIBRIUM_IDENTIFIABLE
        return base + ("log10_koff",) if self.log10_koff is not None else base
