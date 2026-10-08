#!/usr/bin/env python3
"""Port the dm_control rodent hindlimb from position-servo actuators to antagonist
MuJoCo muscle pairs.

Why: the stock model drives each hindlimb DOF with ONE signed position-servo
(dyntype=filter, biastype=affine => proportional feedback to a setpoint). That is wrong
for pharmacology in three ways:
  1. flexor and extensor collapse into one signed number, so the spinal CPG's reciprocal
     inhibition (glycine / GABA-A) -- our drug target -- is unrepresentable;
  2. the servo's feedback COMPENSATES for drug-induced weakness, masking the effect;
  3. no force-length / force-velocity dependence, so no substrate for spindle afferents.

Each hindlimb DOF becomes two muscles: <joint>_ext (gear +1) and <joint>_flx (gear -1).
MuJoCo muscles pull only (force in [-F,0]), so the opposing gears span both directions.
NOTE: _ext/_flx is a LENGTH convention: ext LENGTHENS as qpos rises (gear +1). NOTE muscles are pull-only (force <= 0), so gear +1 yields NEGATIVE joint torque -- see circuitpharm/plant.py; an earlier comment here had this sign backwards; the anatomical
flexor/extensor assignment must be validated against joint kinematics before mapping
motoneuron pools onto it. Non-hindlimb actuators are left untouched.
"""
import os, pathlib, dm_control
from lxml import etree
# SIDE-EFFECT GUARD added 2026-10-07. This script previously did its work at MODULE level,
# so merely importing it ran it. For `port_muscle.py` that regenerated the MuJoCo body
# model, and for `build_kb.py`/`build_compounds.py` it rewrote the pharmacology database --
# both destructive, both triggered by any tool that imports or scans the package. The body
# is unchanged; it is now reached only when the file is run directly.


def main():

    SRC = pathlib.Path(dm_control.__file__).parent / "locomotion/walkers/assets/rodent.xml"
    OUT = pathlib.Path(__file__).resolve().parent.parent / "models" / "rodent_muscle.xml"
    OUT.parent.mkdir(exist_ok=True)

    HINDLIMB = [f"{j}_{s}" for s in ("L", "R")
                for j in ("hip_{}_supinate", "hip_{}_abduct", "hip_{}_extend",
                          "knee_{}", "ankle_{}", "toe_{}")]
    HINDLIMB = [f"hip_{s}_supinate" for s in "LR"] + [f"hip_{s}_abduct" for s in "LR"] \
             + [f"hip_{s}_extend"   for s in "LR"] + [f"knee_{s}"       for s in "LR"] \
             + [f"ankle_{s}"        for s in "LR"] + [f"toe_{s}"        for s in "LR"]

    tree = etree.parse(str(SRC)); root = tree.getroot()

    # muscle actuators need a lengthrange; let the compiler compute it
    comp = root.find("compiler")
    comp.set("autolimits", "true")
    lr = etree.SubElement(comp, "lengthrange")
    lr.set("mode", "muscle"); lr.set("useexisting", "true"); lr.set("accel", "20")
    lr.set("timeconst", "1"); lr.set("tolrange", "0.05")

    # asset paths: the stock XML uses relative meshdir + a render-only skin. We run headless,
    # so point meshdir at the upstream asset dir and drop the skin entirely.
    for c in root.iter("compiler"):
        c.set("meshdir", str(SRC.parent))
    for a in root.iter("asset"):
        for sk in a.findall("skin"):
            a.remove(sk)

    act = root.find("actuator")
    replaced, pairs = [], 0
    for el in list(act):
        name = el.get("name")
        if name not in HINDLIMB:
            continue
        joint = el.get("joint")
        # peak isometric force from the servo's forcerange, with headroom for the
        # force-length curve (a muscle is rarely at optimum length)
        fr = el.get("forcerange", "-0.2 0.2").split()
        F = abs(float(fr[1])) * 2.0
        idx = list(act).index(el)
        act.remove(el)
        for k, (suf, gear) in enumerate((("ext", "1"), ("flx", "-1"))):
            m = etree.Element("muscle")
            m.set("name", f"{joint}_{suf}")
            m.set("joint", joint)
            m.set("gear", gear)
            m.set("force", f"{F:.6g}")
            m.set("ctrllimited", "true")
            m.set("ctrlrange", "0 1")
            m.set("timeconst", "0.01 0.04")      # tau_act, tau_deact (s)
            # MUST set range explicitly: the stock <default><general gainprm="0.01"/> (for its
            # position servos) otherwise leaks into gainprm[0] and corrupts the muscle
            # operating range, leaving a dead zone at neutral posture. Centred on 1.0 so the
            # force-length peak sits at the joint midpoint and both antagonists can pull
            # across the whole range.
            m.set("range", "0.7 1.3")
            m.set("lmin", "0.5"); m.set("lmax", "1.6")
            m.set("vmax", "1.5"); m.set("fpmax", "1.3"); m.set("fvmax", "1.2")
            act.insert(idx + k, m)
        replaced.append(name); pairs += 1

    root.set("model", "rat_muscle_hindlimb")
    tree.write(str(OUT), pretty_print=True, encoding="utf-8", xml_declaration=False)
    print(f"replaced {pairs} servo actuators with {pairs*2} muscles")
    print("replaced:", ", ".join(replaced))
    print("->", OUT)


if __name__ == "__main__":
    main()
