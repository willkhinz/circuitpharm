"""Mechanistic decomposition of allosteric modulation into physical drivers.

THE IDENTITY IS NOW EXACT, which it was not. The technical spec (§3) and this module's own
docstring both asserted

    ln G = D_ln(Occupancy) + D_ln(Kinetic redistribution) + D_ln(Operating-point headroom)

and the code did not satisfy it. The third term was `(Po_max - Po_base)/(Po_max - Po_pam)`,
a proximity-to-ceiling ratio, which is not a factor of `Po_pam/Po_base` -- so the three
terms did not sum to the log gain, and could not at any parameters. The percentages were
then formed as `|D_i| / sum|D_j|`, which discards sign and is not permutation-averaged, so
it was not a Shapley value either (roadmap §5.9, P4-3).

WHAT REPLACES IT. For the 5-state scheme, write the open probability as a product of three
conditional factors that are exact by construction:

    Po  =  Occ  x  (P_A2 / Occ)  x  (Po / P_A2)

      Occ   = P_AR + P_A2R + P_A2O + P_A2D      any ligand bound
      P_A2  = P_A2R + P_A2O + P_A2D             doubly liganded
      Po    = P_A2O                             doubly liganded AND open

so, between modulated and baseline,

    ln G = D_ln(occupancy) + D_ln(double-occupancy share) + D_ln(gating share)

which is an identity rather than an approximation -- `tests/test_decomposition.py` holds it
to rtol 1e-10 over random parameter draws in both regimes. Each term reads cleanly:

  occupancy           how much ligand is bound at all. At 0.4 uM against a ~24 uM K_d
                      almost nothing is, which is where a tonic PAM has room.
  double-occupancy    what fraction of bound receptors carry TWO agonists. Slowing
                      unbinding shifts the population up the binding ladder; this term
                      reports that shift.
  gating share        what fraction of doubly-liganded receptors are open rather than shut
                      or desensitised. Bounded above by E/(1 + E + D) -- the quantity
                      "headroom" was reaching for. At a saturating synaptic transient the
                      baseline already sits near that bound, which is why the phasic pool
                      has so little room.

BECAUSE THE DECOMPOSITION IS EXACTLY ADDITIVE IN LOGS, the Shapley value of each factor
equals its own term: for an additive game the marginal contribution is order-independent,
so permutation averaging changes nothing. The shares below are therefore genuine Shapley
shares, and they are signed. The fields are still renamed `log_share_pct_*`, because
"Shapley" invites a reader to assume a cooperative-game computation that is not happening
here; this docstring states the equivalence instead.

WHAT THE EXACT DECOMPOSITION IMMEDIATELY SHOWS, and it is not what the roadmap expected.
Measured at the provisional parameters, PAM 2.5x:

    regime   fold gain     occupancy   double-occupancy   gating
    tonic      5.8513x       +61.3%          +38.7%        0.0%
    phasic     1.0012x        +2.0%          +98.0%        0.0%

THE GATING TERM IS IDENTICALLY ZERO IN BOTH REGIMES, and the gating share sits exactly at
its ceiling (0.166169 = E/(1+E+D)) in both. That is not a numerical accident: an
affinity-type modulator divides k_off and touches neither beta/alpha nor d/r, so the
equilibrium distribution WITHIN the doubly-liganded states is untouched, and
P(open | doubly bound) is a constant independent of agonist concentration and of the PAM.

So at equilibrium there is no "operating point headroom" to spend: the gating share is
pinned at its ceiling always. Roadmap Advancement 1 expected "operating point saturation
and desensitisation trapping govern the phasic ceiling"; the decomposition says the phasic
ceiling is entirely a BINDING-side phenomenon -- the baseline already sits high on the
binding ladder, so there is nowhere further to move it. The whole phasic/tonic divergence
lives in the first two terms.

The third factor only becomes informative where P(open | A2) can actually move: a
GATING-type modulator (beta scaled, as neurosteroids do), or a time-varying waveform where
the distribution is away from stationarity. Both are in scope for P6, and the decomposition
is now the right instrument to measure them with.

AND IT SEPARATES THE TWO MECHANISMS CLEANLY, which is what Advancement 1 asked for --
"determine which factor is necessary and which is sufficient". Measured at 0.4 uM:

    condition         occupancy     double-occ share   gating share   (its ceiling)
    baseline          3.1658e-02        0.173915          0.166169      0.166169
    affinity 2.5x     9.3426e-02        0.344830          0.166169      0.166169
    gating   2.5x     3.2985e-02        0.208237          0.332536      0.332536

An affinity modulator moves terms 1 and 2 and leaves term 3 at its ceiling. A gating
modulator roughly doubles term 3 (tracking its own raised ceiling) and barely touches term
1. So the decomposition is a MECHANISM ASSAY, not just an attribution: given a measured
gain and the three shares, the two modulator types are distinguishable.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..models.kinetic_jw95 import KineticAllosteryModel
from ..results import Quantity, ResultSet, Tier

#: Indices into the 5-state vector [R, AR, A2R, A2O, A2D].
_AR, _A2R, _A2O, _A2D = 1, 2, 3, 4


@dataclass(frozen=True)
class FactorialDecompositionResult:
    """Exact three-factor attribution of a modulation gain. See the module docstring."""

    condition: str
    gaba_conc_um: float
    pam_factor: float
    total_fold_gain: float
    log_gain: float
    delta_log_occupancy: float
    delta_log_double_occupancy: float
    delta_log_gating: float
    #: Signed shares of the log gain. Genuine Shapley shares, because the decomposition is
    #: exactly additive in logs -- see the module docstring on why the name changed anyway.
    log_share_pct_occupancy: float
    log_share_pct_double_occupancy: float
    log_share_pct_gating: float
    #: The gating share's analytic ceiling E/(1 + E + D), and the baseline's level against
    #: it. This is what the old "headroom" term was reaching for, reported as a LEVEL
    #: rather than smuggled into the product where it does not belong.
    gating_share_baseline: float
    gating_share_ceiling: float
    residual: float             # ln G minus the three terms; must be ~0

    def check(self, rtol: float = 1e-10) -> None:
        """Raise if the identity does not hold. Called by `decompose_modulation_gain`."""
        if abs(self.residual) > rtol * max(abs(self.log_gain), 1.0):
            raise AssertionError(
                f"the log-gain identity failed: ln G = {self.log_gain!r} but the three "
                f"terms sum to {self.log_gain - self.residual!r} (residual "
                f"{self.residual!r}). This is an identity by construction, so a failure "
                f"means the state partition is wrong -- not that a tolerance needs widening.")


def _factors(dist: np.ndarray) -> tuple[float, float, float]:
    """(occupancy, double-occupancy share, gating share) from a 5-state vector."""
    occ = float(dist[_AR] + dist[_A2R] + dist[_A2O] + dist[_A2D])
    a2 = float(dist[_A2R] + dist[_A2O] + dist[_A2D])
    po = float(dist[_A2O])
    if occ <= 0.0 or a2 <= 0.0:
        return occ, 0.0, 0.0
    return occ, a2 / occ, po / a2


def decompose_modulation_gain(
    model: KineticAllosteryModel,
    gaba_um: float,
    pam_factor: float = 2.50,
    condition_label: str = "custom",
) -> FactorialDecompositionResult:
    """Decompose a modulation gain into occupancy, double-occupancy and gating terms.

    `model` is REQUIRED (roadmap P0-13). `pam_factor` is an EQUILIBRIUM EC50 fold-shift,
    which for this scheme is exactly the factor dividing K_d -- NOT the same quantity as
    `gabaa_kinetics.calibrate_pam`'s PEAK-calibrated multiplier, which is 2.945 for a
    nominal 2.5 (roadmap P6-2).
    """
    if model is None:
        raise TypeError(
            "decompose_modulation_gain() requires an explicit `model`; see roadmap P0-13.")
    g = max(float(gaba_um), 1e-9)
    pf = float(pam_factor)
    if pf < 1.0:
        raise ValueError(
            f"pam_factor={pf!r} is below 1.0. Negative allosteric modulation is out of "
            f"domain for this scheme's callers; see models/base.py.")

    base = model.state_distribution(g, pam_factor=1.0)
    pam = model.state_distribution(g, pam_factor=pf)
    occ_b, dbl_b, gate_b = _factors(base)
    occ_p, dbl_p, gate_p = _factors(pam)

    po_b, po_p = float(base[_A2O]), float(pam[_A2O])
    if po_b <= 0.0:
        raise ValueError(
            f"baseline open probability is {po_b!r} at {g} uM, so a fold-gain is "
            f"undefined. Raise the agonist concentration or check the parameters.")

    gain = po_p / po_b
    log_gain = float(np.log(gain))
    d_occ = float(np.log(occ_p / occ_b))
    d_dbl = float(np.log(dbl_p / dbl_b))
    d_gate = float(np.log(gate_p / gate_b))
    residual = log_gain - (d_occ + d_dbl + d_gate)

    # NORMALISED BY THE SIGNED TOTAL, not by the sum of absolute values. The docstring and
    # `test_the_shares_sum_to_one_hundred_percent` both assert that the three shares sum to
    # 100%, and that identity holds only against the signed sum: with any negative term,
    # sum(|d_i|) > |sum(d_i)| and the shares silently sum to less than 100%.
    #
    # HARDENING, NOT A LIVE BUG FIX, and the distinction is worth recording. Probing 175
    # (gaba_um, pam_factor) combinations across 10^-3 to 10^5 uM found NO negative term: for
    # an affinity-type PAM the gating term is identically 0.0 and both occupancy terms are
    # non-negative, so abs-sum equals signed-sum and the shares do sum to 100% today. The
    # premise becomes reachable the moment a modulator with a negative component exists --
    # a NAM, or a gating modulator that trades double occupancy for gating -- and NAMs are
    # currently rejected by every model's `apply_pam`. So this is correct-by-construction
    # rather than a repair, and the identity no longer depends on an invariant enforced
    # somewhere else.
    total = d_occ + d_dbl + d_gate
    # abs(), because the total is now SIGNED. `total > 1e-12` was correct for a sum of
    # magnitudes and wrong the moment the sum could be negative: a modulator with a net
    # NEGATIVE log gain would fall through to the zero branch and report 0%/0%/0% shares for
    # a real effect. My own regression, introduced by the line above.
    if abs(total) > 1e-12:
        pct = (100.0 * d_occ / total, 100.0 * d_dbl / total, 100.0 * d_gate / total)
    else:
        # All three terms vanish: the modulator did nothing. Reporting thirds (the old
        # behaviour) would invent structure, so the shares are zero and the unit gain says
        # why.
        pct = (0.0, 0.0, 0.0)

    e = model.gating_efficacy
    d_ratio = model.desens_ratio
    res = FactorialDecompositionResult(
        condition=condition_label,
        gaba_conc_um=g,
        pam_factor=pf,
        total_fold_gain=gain,
        log_gain=log_gain,
        delta_log_occupancy=d_occ,
        delta_log_double_occupancy=d_dbl,
        delta_log_gating=d_gate,
        log_share_pct_occupancy=pct[0],
        log_share_pct_double_occupancy=pct[1],
        log_share_pct_gating=pct[2],
        gating_share_baseline=gate_b,
        gating_share_ceiling=float(e / (1.0 + e + d_ratio)),
        residual=float(residual),
    )
    res.check()
    return res


def compare_phasic_tonic_mechanisms(
    model: KineticAllosteryModel,
    tonic_gaba_um: float = 0.40,
    phasic_gaba_um: float = 1000.0,
    pam_factor: float = 2.50,
) -> dict[str, FactorialDecompositionResult]:
    """Run the decomposition in both regimes.

    `model` is REQUIRED: the former `model or KineticAllosteryModel()` fallback ran the
    decomposition at provisional chimeric defaults (roadmap P0-13).
    """
    if model is None:
        raise TypeError(
            "compare_phasic_tonic_mechanisms() requires an explicit `model`; see roadmap "
            "P0-13 for why the default was removed.")
    return {
        "tonic": decompose_modulation_gain(model, tonic_gaba_um, pam_factor=pam_factor,
                                           condition_label="tonic"),
        "phasic": decompose_modulation_gain(model, phasic_gaba_um, pam_factor=pam_factor,
                                            condition_label="phasic"),
    }


def decomposition_report(model: KineticAllosteryModel, **kw) -> ResultSet:
    """The decomposition as tiered quantities (roadmap P2).

    UNCALIBRATED, not VALIDATED. The identity is exact -- that part is algebra -- but the
    NUMBERS depend on the fitted rates, which rest on three anchors taken from recalled
    literature ranges plus one declared convention (`gabaa_kinetics.FIT_FIXED_ALPHA`), and
    on `d`/`r`, which no anchor constrains at all. The ordering and the mechanism are
    usable; the magnitudes are not quotable as measurements.
    """
    both = compare_phasic_tonic_mechanisms(model, **kw)
    rs = ResultSet(f"mechanistic decomposition at PAM "
                   f"{both['tonic'].pam_factor:.2f}x (EQUILIBRIUM EC50 shift)")
    prov = ("exact three-factor identity Po = Occ x (P_A2/Occ) x (Po/P_A2); the "
            "decomposition is algebra, the magnitudes come from the fitted rates")
    promote = ("fit the rates to sourced data on their own observables (P4) -- the "
               "identity itself needs nothing")
    for name, r in both.items():
        rs.add(Quantity(f"{name}_fold_gain", r.total_fold_gain, Tier.UNCALIBRATED, "x",
                        provenance=prov, promote_by=promote))
        for lab, val in (("occupancy", r.log_share_pct_occupancy),
                         ("double_occupancy", r.log_share_pct_double_occupancy),
                         ("gating", r.log_share_pct_gating)):
            rs.add(Quantity(f"{name}_share_{lab}", val, Tier.UNCALIBRATED, "%",
                            provenance=prov, promote_by=promote,
                            caveats=("a signed share of the log gain; the three sum to "
                                     "100% because the decomposition is exactly additive",)))
        rs.add(Quantity(
            f"{name}_gating_share_baseline", r.gating_share_baseline, Tier.UNCALIBRATED,
            provenance=f"{prov}; analytic ceiling E/(1+E+D) = {r.gating_share_ceiling:.6f}",
            promote_by=promote,
            caveats=("distance from that ceiling is what 'headroom' meant; it is reported "
                     "as a level because it is not a factor of the gain",)))
    return rs
