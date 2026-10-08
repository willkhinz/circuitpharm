#!/usr/bin/env python3
"""Build the compound survey table: every researched compound evaluated as a candidate
arm for an alcohol substitute, with why it does or does not qualify.

Columns of judgement:
  eth_like     does it substitute for ethanol's discriminative stimulus? (the subjective proxy)
  sedation     sedation / ataxia / motor impairment
  resp         respiratory depression
  hepatotox    liver injury signal
  ceiling      does the dose-response plateau (overdose protection)?
  verdict      CANDIDATE / BACKUP / REJECT / TOOL / GAP
"""
import sqlite3, pathlib
# SIDE-EFFECT GUARD added 2026-10-07. This script previously did its work at MODULE level,
# so merely importing it ran it. For `port_muscle.py` that regenerated the MuJoCo body
# model, and for `build_kb.py`/`build_compounds.py` it rewrote the pharmacology database --
# both destructive, both triggered by any tool that imports or scans the package. The body
# is unchanged; it is now reached only when the file is run directly.


def main():

    DB = pathlib.Path(__file__).resolve().parent.parent / "data" / "pharmacology.db"
    c = sqlite3.connect(DB)
    c.execute("DROP TABLE IF EXISTS compounds")
    c.execute("""CREATE TABLE compounds(
      name TEXT, arm TEXT, target TEXT, status TEXT,
      eth_like TEXT, sedation TEXT, resp TEXT, hepatotox TEXT, ceiling TEXT,
      verdict TEXT, note TEXT)""")

    # arm: GABA | NMDA | BOTH | OTHER
    X = [
    # ---------------------------------------------------------------- GABA-A, subtype-selective
    ("MP-III-022","GABA","alpha5 PAM; 300% efficacy at a5b3g2; NONMODULATORY at a1","preclinical",
     "a5 route: a5 agonists mimic ethanol","NONE at 1-10 mg/kg (mild myorelaxation at 10)","untested",
     "untested","PAM (ceiling expected)","CANDIDATE",
     "Best GABA arm found. a5 mediates ethanol's discriminative stimulus; a1 (sedation/respiratory) untouched. a5 is forebrain-restricted, low in pons/medulla."),
    ("QH-ii-066","GABA","alpha5-selective agonist","research tool",
     "YES - mimics ethanol in squirrel monkeys, blocked by L-655,708","YES - sedation and ataxia","untested",
     "untested","agonist","TOOL",
     "Proves the a5 route works subjectively, but retains sedation/ataxia - likely residual a1 activity. Use as the positive control, not the candidate."),
    ("TPA023 / MK-0777","GABA","a2/a3 partial agonist; SILENT ANTAGONIST at a1 and a5","Phase 2 (GAD)",
     "a2/a3 route plausible (see HZ-166)","healthy volunteers: no effect on alertness/memory/body sway; PATIENTS: dizziness, drowsiness, incoordination",
     "untested","untested","partial agonist = ceiling","BACKUP",
     "Clinical precedent that a1-sparing works in humans. Mixed sedation picture between volunteers and patients. Note it is a SILENT ANTAGONIST at a5, so it forfeits the a5 subjective route."),
    ("TPA023B","GABA","a2/a3 partial agonist","Phase 1",
     "untested vs ethanol","well tolerated at >50% occupancy","untested","untested","partial agonist","BACKUP",
     "Improved TPA023. Sedation in the related MRK-409 was attributed specifically to a1 efficacy - direct evidence that a1 is the sedation subtype."),
    ("HZ-166","GABA","a2/a3 PAM","preclinical (ester, poor bioavailability)",
     "YES - substituted FULLY for ethanol in rhesus monkeys WITHOUT altering response rates","none at effective dose",
     "untested","untested","PAM","CANDIDATE",
     "Second clean route: full ethanol substitution with no rate suppression (= no sedation) in primates. Metabolically labile ester limits it; KRM-II-81 is the fix."),
    ("KRM-II-81","GABA","a2/a3 PAM (oxazole bioisostere of HZ-166)","DISCONTINUED",
     "parent HZ-166 substitutes fully","devoid of sedation/myorelaxation/cognitive impairment in rat pain models; sedation at HIGH dose",
     "untested","untested","PAM","BACKUP",
     "Orally bioavailable. Explicitly avoids a1 AND a5 - so it forfeits the a5 subjective route and must rely on a2/a3. Failed ethanol EEG match in the 2025 study."),
    ("YT-III-31","GABA","a3 PAM","preclinical",
     "YES - substituted fully for ethanol in rhesus, no rate change","none reported","untested","untested","PAM","CANDIDATE",
     "a3-only route also fully substitutes. Narrowest selectivity that still works subjectively."),
    ("L-838417","GABA","partial agonist a2/a3/a5; ANTAGONIST at a1","research tool",
     "plausible (hits a5)","reduces anxiety with few sedative effects","untested","untested","partial","TOOL",
     "Parent of TPA023. Hits a5 AND spares a1 - the theoretically ideal profile. Poor pharmacokinetics prevented development."),
    ("L-655,708","GABA","a5 inverse agonist / antagonist","research tool",
     "BLOCKS ethanol's stimulus","n/a","n/a","untested","n/a","TOOL",
     "The pharmacological proof that a5 carries ethanol's subjective effect. Also reduced alcohol drinking in rhesus monkeys."),
    ("Basmisanil","GABA","a5 NAM, highly selective","Phase 2, target engagement shown in man",
     "would OPPOSE ethanol","n/a","n/a","no signal reported","n/a","TOOL",
     "Demonstrates that highly a5-selective drugs can be developed and dosed in humans. Wrong direction (negative modulator) but de-risks the target."),
    ("alpha5IA","GABA","a5 inverse agonist","Phase 1","would oppose ethanol","n/a","n/a","untested","n/a","TOOL",
     "Cognition enhancer. Same target, opposite direction."),
    # ---------------------------------------------------------------- GABA-A, non-selective / rejected
    ("Zolpidem","GABA","a1-preferring PAM","approved",
     "NO - a1 is not the ethanol subtype","STRONG (a1)","yes","no class signal","PAM","REJECT",
     "Wrong subtype entirely. Increased BETA power, not ethanol's delta, in the 2025 EEG comparison."),
    ("Alprazolam / diazepam","GABA","non-selective BZ-site PAM","approved",
     "ethanol generalises TO them, but NOT the reverse (one-way asymmetry)","yes","modest (-14 to -19% MV at sedative dose)",
     "no class signal","PAM - plateaus, no apnoea","BACKUP",
     "The calibration reference. Non-selective means a1 is engaged, so sedation and respiratory cost come along."),
    ("Ocinaplon","GABA","low-potency non-selective BZ-site PAM","HALTED - hepatotoxicity",
     "untested","anxiolytic WITHOUT typical BZ side effects","untested","YES - halted, liver injury in one patient","PAM","REJECT",
     "The profile we want, killed by the harm we are trying to remove."),
    ("Alpidem","GABA","BZ-site PAM (imidazopyridine)","WITHDRAWN - hepatotoxicity",
     "untested","anxiolytic, non-sedating","untested","YES - market withdrawal","PAM","REJECT",
     "Reached market then withdrawn for liver injury. Zolpidem, same chemical family, is fine - so idiosyncratic, not mechanism-based."),
    ("Panadiplon","GABA","quinoxalinone high-affinity partial agonist","HALTED Phase 1 - hepatotoxicity",
     "YES - a5 agonist, mimics ethanol in monkeys","little sedative or amnestic effect","untested",
     "YES - human hepatic toxicity NOT predicted by rat, dog OR monkey; mitochondrial/Reye's-like via a carboxylic acid metabolite","partial","REJECT",
     "Most alarming entry in the table: exactly the desired profile (a5, ethanol-mimetic, non-sedating) and it injured human livers invisibly to three preclinical species."),
    ("Kava / kavalactones","GABA","GABA-A potentiation + MAO-B inhibition + NA/DA reuptake","marketed supplement; banned/withdrawn in several EU markets",
     "anxiolytic, alcohol-adjacent use","mild","untested",
     "YES - transient enzyme rise to fulminant failure and death; excessive ALCOHOL INTAKE is a listed risk factor","unknown","REJECT",
     "Directly analogous warning: a GABA-A-active hepatotoxin consumed by drinkers. 98% of reported cases had liver injury; 83% recovered on withdrawal."),
    ("Bretazenil / abecarnil / pazinaclone","GABA","BZ-site partial agonists","failed clinically (various)",
     "untested","higher doses needed for BZ-like side effects","untested","class concern","partial = ceiling","REJECT",
     "The partial-agonist generation that did not succeed; useful as the cautionary base rate."),
    ("Gaboxadol / THIP","GABA","a4/delta extrasynaptic DIRECT AGONIST","failed (insomnia)",
     "NO - ethanol-trained rats respond on the VEHICLE lever","yes","yes","untested",
     "NO CEILING (direct agonist)","REJECT",
     "Double disqualification: fails ethanol substitution AND forfeits the overdose ceiling. Also delta is not required for ethanol's stimulus."),
    ("Muscimol / barbiturates","GABA","direct agonists","n/a (tool / legacy)",
     "barbiturates: ethanol generalises to them","strong","strong, NO ceiling","barbiturates no","NO CEILING","REJECT",
     "No ceiling = the alcohol overdose problem reproduced. The modelled uncapped arm fell to 8% ventilation at 26x."),
    # ---------------------------------------------------------------- neurosteroids
    ("Allopregnanolone / brexanolone","GABA","neurosteroid PAM, distinct site, non-subtype-selective","approved (IV)",
     "ethanol substitutes MORE readily for neurosteroids than benzodiazepines","yes","BOXED WARNING: excessive sedation, SUDDEN LOSS OF CONSCIOUSNESS (4%)",
     "no signal","HIGH efficacy - ceiling likely LOST","REJECT as primary",
     "Best subjective efficiency, worst ceiling risk. Mandatory continuous pulse oximetry in use."),
    ("Zuranolone (SAGE-217)","GABA","oral neurosteroid PAM","approved",
     "neurosteroid route","yes","boxed warning for impaired driving","no signal","high efficacy","BACKUP",
     "Oral and approved. Produces tonic-current potentiation lasting BEYOND application (unlike ganaxolone)."),
    ("Ganaxolone","GABA","3b-methyl allopregnanolone analogue","approved (CDKL5)",
     "neurosteroid route","sedation/hypoactivity COMPARABLE TO MIDAZOLAM","no extra penalty per unit effect","no signal",
     "high efficacy","BACKUP",
     "Key datum: no extra respiratory penalty versus a benzodiazepine. Lacks the persistent tonic effect the others have."),
    ("neurosteroids as a class","GABA","GABA-A PAM, distinct transmembrane site","various",
     "best subjective match","yes - non-selective across subtypes so a1 IS engaged","per-unit comparable to midazolam",
     "no class signal","differ in POTENCY but NOT intrinsic efficacy","REJECT as primary",
     "DISQUALIFYING for my earlier recommendation: there is no 'low-efficacy neurosteroid' to choose, because they do not differ in intrinsic efficacy. And being subtype-non-selective they cannot spare a1."),
    # ---------------------------------------------------------------- NMDA arm
    ("Esmethadone (REL-1017)","NMDA","low-affinity low-potency uncompetitive channel blocker; opioid-INACTIVE d-isomer of methadone","Phase 3",
     "untested vs ethanol","mild","no signal reported",
     "no signal","low affinity = self-limiting","CANDIDATE",
     "Strongest NMDA arm on safety evidence: NO meaningful abuse potential in recreational drug users (tested against oxycodone AND ketamine), NOT neurotoxic in rats, no reinforcement/dependence/withdrawal. Its mirror image is an opioid - a clean enantiomer separation."),
    ("Lanicemine (AZD6765)","NMDA","LOW-TRAPPING channel blocker (54% vs ketamine 86%)","Phase 2b",
     "untested vs ethanol","minimal","no signal","no signal","low trapping","CANDIDATE",
     "Non-psychotomimetic by a mechanism other than subunit selectivity: 150 mg IV in 22 patients with no psychosis or dissociation."),
    ("Traxoprodil (CP-101,606)","NMDA","GluN2B-selective","DISCONTINUED - QT prolongation",
     "GluN2B is ethanol's own NMDA mechanism","dose-related dissociation and amnesia; favourable RELATIVE to unselective",
     "untested","untested","n/a","REJECT (compound), target OK",
     "Killed by hERG/QTc, which no model here can predict. Also corrects my prediction: GluN2B-selective is LESS dissociative, not non-dissociative."),
    ("Rislenemdaz (CERC-301)","NMDA","oral GluN2B-selective","Phase 2 (failed on efficacy)",
     "GluN2B route","not reported as problematic","no signal","no signal","n/a","CANDIDATE",
     "Safety pharmacology and neurotoxicity studies CLEAN; failed only on antidepressant efficacy, which is irrelevant here since a subjective effect is wanted."),
    ("Radiprodil","NMDA","GluN2B-selective","Phase 2 (epilepsy/TSC)",
     "GluN2B route","somnolence, hypotonia","untested","untested","n/a","BACKUP","Active development in neurodevelopmental indications."),
    ("Ro 25-6981","NMDA","GluN2B-selective","research tool",
     "GluN2B route","did NOT potentiate motor impairment in combination with morphine","untested","untested","n/a","TOOL",
     "Explicitly described as safe for further development in infantile rats. The motor-sparing positive control."),
    ("Ifenprodil","NMDA","GluN2B-selective but dirty (a1-adrenergic, sigma)","marketed (vasodilator, Japan)",
     "GluN2B route","dizziness","untested","untested","n/a","REJECT",
     "a1-adrenergic off-target gives hypotension; selectivity too poor to attribute effects."),
    ("Memantine","NMDA","low-affinity uncompetitive channel blocker","approved",
     "NO - PCP-LIKE discriminative stimulus; did not modulate ethanol discrimination","mild","mild","no","voltage-dependent","REJECT",
     "Subjectively PCP-like rather than ethanol-like. Preferentially blocks extrasynaptic NMDA receptors."),
    ("Dextromethorphan","NMDA","channel blocker + sigma","approved",
     "NO - PCP-like substitution","yes","yes","no","none","REJECT","PCP-like subjective profile; also widely misused."),
    ("Ketamine / MK-801 / PCP","NMDA","high-affinity high-trapping channel blockers","various",
     "MK-801 gives COMPLETE generalisation to ethanol","strong","yes","no","none","REJECT as product, TOOL",
     "MK-801 is the pharmacological proof that NMDA antagonism is part of ethanol's stimulus; all are far too psychotomimetic to be products."),
    ("CGP-40116 / CPP / selfotel","NMDA","competitive antagonists","failed (stroke)",
     "CGP-40116 substitutes 88% for ethanol","CPP: ataxia and reduced muscle tone at higher doses","untested","untested","none","REJECT",
     "Competitive antagonism also substitutes for ethanol, but the class failed on psychotomimetic effects and poor therapeutic index."),
    ("AV-101 / rapastinel / apimostinel","NMDA","glycine-site (GlyB) modulators","FAILED Phase 3 (efficacy)",
     "unlikely - essentially non-psychoactive","none","none","no","n/a","REJECT",
     "Designed to avoid dissociation and they succeed - but they are non-psychoactive, which for this project is the disqualifier rather than the feature."),
    # ---------------------------------------------------------------- both arms in one molecule
    ("Nitrous oxide","BOTH","NMDA antagonist + GABA enhancement + opioid interaction","approved / widely misused",
     "intoxicating, euphoric, alcohol-adjacent","yes","MINIMAL - conscious sedation causes no respiratory depression or desaturation",
     "no","n/a (gas, self-limiting exposure)","TOOL / existence proof",
     "EXISTENCE PROOF that GABA+NMDA together can intoxicate with minimal respiratory depression. Disqualified as a product by irreversible B12/cobalamin inactivation causing neuropathy, hypoxia risk, and delivery."),
    ("Isoflurane / volatile anaesthetics","BOTH","GABA-A PAM + NMDA antagonist","approved",
     "shows BOTH GABA-A-PAM-like and NMDA-antagonist-like discriminative effects","yes","yes","halothane: yes","none","REJECT",
     "Second existence proof that one molecule can carry both components."),
    ("GHB / sodium oxybate","OTHER","GABA-B agonist + extrasynaptic GABA-A","approved (narcolepsy); used for alcohol withdrawal in IT/AT",
     "YES - substitutes for ethanol, cross-tolerant; but symmetrical generalisation only within NARROW dose ranges",
     "yes","yes - narrow therapeutic index","no","NO - very narrow TI","REJECT",
     "Genuinely alcohol-mimetic and clinically used as substitution therapy, but the narrow therapeutic index is precisely the overdose problem this project exists to remove."),
    ("Baclofen","OTHER","GABA-B agonist","approved",
     "substitutes for GHB but NOT equivalent to ethanol","yes","yes","no","n/a","REJECT",
     "R-baclofen is the active enantiomer (~100x). Used in alcohol use disorder, but not an ethanol-like subjective substitute."),
    # ---------------------------------------------------------------- other ethanol targets
    ("ML297 (VU0456810)","OTHER","first selective GIRK1/2 activator, IC50 160 nM","preclinical",
     "GIRK IS an ethanol target (ethanol activates GIRK G-protein-independently)",
     "NO - anxiolytic WITHOUT sedation","untested","untested","n/a","CANDIDATE (adjunct)",
     "Anxiolytic and antiepileptic with no sedation and no addiction-related behaviour, acting on a genuine ethanol target. Promising third arm to offload burden from GABA-A."),
    ("5-chloroindole","OTHER","5-HT3 PAM","research tool",
     "ethanol is a PAM of 5-HT3A (but NOT of 5-HT3AB, even at 200 mM)","n/a","n/a","untested","n/a","GAP",
     "A real ethanol mechanism with an available PAM, but 5-HT3 potentiation causes nausea/emesis - an unwanted effect. Subunit-specific (A vs AB) so selectivity is in principle possible."),
    ("Ivermectin / tropeines","OTHER","glycine receptor PAMs (non-selective; tropeines potentiate a1 GlyR but INHIBIT a3)","approved for other uses",
     "GlyR is an ethanol target","n/a","n/a","n/a","n/a","GAP",
     "No clean selective GlyR PAM exists. Ivermectin is far too promiscuous (GlyR + GABA-A + a7 nACh + GluCl). This remains the main undeveloped ethanol mechanism."),
    ("Dexmedetomidine","OTHER","a2-adrenergic agonist","approved",
     "NO - not an ethanol target","sedation (the distinct 'arousable' kind)","minimal respiratory depression","no","n/a","REJECT",
     "Notable for sedation with minimal respiratory depression, but not ethanol-like and causes bradycardia/hypotension."),
    ("Suvorexant / lemborexant / daridorexant","OTHER","dual orexin receptor antagonists","approved",
     "NO - orexin is not an ethanol target","sleepiness, not intoxication",
     "NONE - no respiratory abnormality to 1000 mg/kg; unlikely to cause overdose death","no","n/a","REJECT (instructive)",
     "The safest profile in the entire table - no respiratory depression, abuse potential below Schedule IV hypnotics - but it produces sleep, not intoxication. Shows safety and the wanted subjective effect are genuinely in tension."),
    ("Alcarelle / 'alcosynth' (GABA Labs)","GABA","undisclosed GABA-targeting","novel food / GRAS route; SENTIA Plus announced for 2026 US launch",
     "designed for it","claimed low","claimed low","claimed low","claimed",
     "COMPETITOR / precedent",
     "Answers the regulatory question raised earlier: they are pursuing NOVEL FOOD (UK/EU/Canada) and GRAS (US), not a drug pathway. Composition undisclosed, claims not independently verified here."),
    ]
    c.executemany("INSERT INTO compounds VALUES(?,?,?,?,?,?,?,?,?,?,?)", X)
    c.commit()
    n = c.execute("SELECT COUNT(*) FROM compounds").fetchone()[0]
    print(f"compounds: {n}")
    for v in ("CANDIDATE", "BACKUP", "TOOL", "REJECT", "GAP"):
        k = c.execute("SELECT COUNT(*) FROM compounds WHERE verdict LIKE ?", (f"%{v}%",)).fetchone()[0]
        print(f"  {v:<10} {k}")
    print("\nhepatotoxicity signals:")
    for r in c.execute("SELECT name, hepatotox FROM compounds WHERE hepatotox LIKE 'YES%'"):
        print(f"  {r[0]:<34} {r[1][:70]}")
    c.close()


if __name__ == "__main__":
    main()
