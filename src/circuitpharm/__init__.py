"""circuitpharm — receptor-resolved pharmacology coupled to spiking neural circuits.

WHAT THIS PACKAGE IS. A multiscale instrument that takes a compound's MEASURED receptor
activity and propagates it to circuit-level and motor/autonomic output, carrying an
explicit reliability tier on every number it produces.

    measured subtype efficacies          <- INPUT (an experiment, not a prediction)
      |  gabaa_kinetics : Markov gating scheme
      v
    tonic / phasic conductance changes
      |  subtypes : regional subunit composition and localisation
      v
    per-region drug sensitivity
      |  cpg, resp, circuit : spiking networks
      v
    circuit dynamics
      |  plant : MuJoCo muscle body
      v
    motor / autonomic output

WHAT IT DELIBERATELY DOES NOT DO. Two links in the chain from structure to behaviour are
not closed, and the package declines to fake either:

  structure -> EFFICACY. For an allosteric modulator, efficacy is not affinity: it is the
    shift in the conformational equilibrium between shut, open and desensitised states.
    Free-energy methods reproduce affinities, and Markov state models have captured a
    putative GABA-A open state, but there is no predictive efficacy calculation -- the
    field still publishes experimental screening designs for exactly this quantity, and
    measured efficacy is additionally PROBE-DEPENDENT (it varies with the agonist
    concentration assayed). So efficacy is an INPUT here, by design.

  circuit -> ARBITRARY BEHAVIOUR. Behaviour-generating controllers that produce naturalistic
    rodent behaviour are artificial networks whose units have no conductances, so there is
    nothing to inject a receptor current into. Full biophysical closed loops exist
    (BAAIWorm, in a 302-neuron animal) but not for mammals. This package therefore targets
    behaviours whose generating circuits ARE biophysically modellable -- breathing, spinal
    reflexes, locomotor rhythm and gait -- which happen to be the physiological endpoints
    of interest, and which map onto standard assays (plethysmography, H-reflex, gait
    analysis, righting reflex, rotarod).

  Subjective or affective state is outside the architecture entirely. Drug discrimination
  requires an animal that learns; no increase in resolution crosses that boundary.

RELIABILITY TIERS ARE ENFORCED IN CODE (see `results`). Reading a VOID quantity raises
rather than returning a number with a caveat attached, because caveats get dropped when
numbers are copied -- which happened in this project's own reporting.

Quick start:
    from circuitpharm import calibration_report, Tier
    print(calibration_report())          # what is and is not anchored

    from circuitpharm.subtypes import PROFILES
    from circuitpharm.gabaa_kinetics import fit_scheme, derive, calibrate_pam
"""
from .results import Tier, Quantity, ResultSet, VoidQuantityError
from .config import (
    RESP_OP, SYNAPTIC_PULSE, AMBIENT_GABA_UM, IA_SCALE,
    FOREBRAIN_GLUN2B, BRAINSTEM_GLUN2B, EUPNOEA_BAND,
    CALIBRATIONS, Calibration, calibration_report,
)

__version__ = "0.1.0"

__all__ = [
    "Tier", "Quantity", "ResultSet", "VoidQuantityError",
    "RESP_OP", "SYNAPTIC_PULSE", "AMBIENT_GABA_UM", "IA_SCALE",
    "FOREBRAIN_GLUN2B", "BRAINSTEM_GLUN2B", "EUPNOEA_BAND",
    "CALIBRATIONS", "Calibration", "calibration_report",
    "__version__",
]
