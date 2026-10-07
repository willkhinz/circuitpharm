#!/usr/bin/env python3
"""Build data/pharmacology.db — structured knowledge base for the alcohol-substitute project.
Idempotent: drops and rebuilds. Every row carries a source_key into `sources`."""
import sqlite3, os, pathlib

DB = pathlib.Path(__file__).resolve().parent.parent / "data" / "pharmacology.db"
DB.parent.mkdir(exist_ok=True)
if DB.exists(): DB.unlink()
c = sqlite3.connect(DB); x = c.executescript

x("""
CREATE TABLE sources(
  key TEXT PRIMARY KEY, citation TEXT, url TEXT, year INT, kind TEXT);

CREATE TABLE ethanol_targets(          -- the mechanism panel
  target TEXT, subunit TEXT, direction TEXT,
  conc_low_mM REAL, conc_high_mM REAL, magnitude TEXT,
  relevant_at_social_dose INT,         -- 1 = active in 5-20 mM window
  evidence_quality TEXT, note TEXT, source_key TEXT);

CREATE TABLE discrimination(           -- ethanol-substitution: the subjective proxy
  training_drug TEXT, test_drug TEXT, doses TEXT,
  pct_substitution REAL, substitutes TEXT, note TEXT, source_key TEXT);

CREATE TABLE eeg_studies(              -- comparability table / admission criteria audit
  study_key TEXT, dose_route TEXT, bec_measured INT, strain TEXT,
  state_scored INT, region_ref TEXT, bands TEXT, limb TEXT,
  aperiodic INT, key_result TEXT);

CREATE TABLE subunit_expression(       -- the selectivity-window question
  receptor TEXT, subunit TEXT, region TEXT, level TEXT,
  bearing TEXT, source_key TEXT);

CREATE TABLE findings(                 -- claims with project bearing
  id INT PRIMARY KEY, claim TEXT, endpoint TEXT, confidence TEXT,
  bearing TEXT, source_key TEXT);

CREATE TABLE sim_results(      -- in-silico, NOT literature. Open-loop spinal circuit.
  condition TEXT, mn_hz REAL, pct_control REAL, period_ms REAL,
  duty REAL, flex_ext_corr REAL, rhythm TEXT, note TEXT);

CREATE TABLE model_facts(              -- verified locally, not from literature
  item TEXT, value TEXT, verified_how TEXT);
""")

S = [
 ("gaba_conc","Ethanol acts directly on extrasynaptic GABA-A subtypes to increase tonic inhibition","https://pmc.ncbi.nlm.nih.gov/articles/PMC2040048/",2007,"primary"),
 ("gaba_decade","The role of GABA-A receptors in acute and chronic effects of ethanol: a decade of progress","https://pmc.ncbi.nlm.nih.gov/articles/PMC2814770/",2010,"review"),
 ("nmda_drivelimit","Inhibition of rat recombinant GluN1/GluN2A and GluN1/GluN2B NMDA receptors by ethanol at drink-drive-limit concentrations","https://pubmed.ncbi.nlm.nih.gov/19394328",2009,"primary"),
 ("nmda_bnst","GluN2B subunit deletion reveals key role in acute and chronic ethanol sensitivity of glutamate synapses in BNST","https://www.pnas.org/doi/10.1073/pnas.1113820109",2012,"primary"),
 ("glyr_ki","Changes in ethanol effects in knock-in mice expressing ethanol-insensitive alpha1 and alpha2 glycine receptor subunits","https://pubmed.ncbi.nlm.nih.gov/38679193/",2024,"primary"),
 ("bk","Ethanol modulation of mammalian BK channels in excitable tissues","https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4256990/",2014,"review"),
 ("girk","Potassium channels as targets for ethanol: GIRK2 null mutant mice","https://pubmed.ncbi.nlm.nih.gov/11454913/",2001,"primary"),
 ("aden","Adenosine A1 receptor blockade mimics caffeine's attenuation of ethanol-induced motor incoordination","https://onlinelibrary.wiley.com/doi/10.1111/j.1742-7843.2004.pto950509.x",2004,"primary"),
 ("subtype2025","Ethanol-induced REM sleep suppression in rats: comparison with subtype-selective GABA-A compounds","https://pmc.ncbi.nlm.nih.gov/articles/PMC12935073/",2025,"primary"),
 ("bec1982","Relationships between blood ethanol levels and ethanol-induced changes in cortical EEG power spectra in the rat","https://pubmed.ncbi.nlm.nih.gov/7121742/",1982,"primary"),
 ("thip_rat","The GABA-A agonist THIP (gaboxadol) increases non-REM sleep and enhances delta activity in the rat","https://pubmed.ncbi.nlm.nih.gov/8930997/",1996,"primary"),
 ("theta1977","Changes in frequency of EEG rhythms of the rat caused by single IP injections of ethanol","https://pubmed.ncbi.nlm.nih.gov/931462",1977,"primary"),
 ("ethbenzo1992","Effects of ethanol on human sleep EEG power spectra differ from benzodiazepine receptor agonists","https://pubmed.ncbi.nlm.nih.gov/1326982/",1992,"primary"),
 ("mixture","Generalisation of ethanol with drug mixtures containing a positive modulator of the GABA-A receptor and an NMDA antagonist","https://pubmed.ncbi.nlm.nih.gov/11077078/",2000,"primary"),
 ("mk801","The NMDA receptor antagonist MK-801 produces ethanol-like discrimination in the rat","https://pubmed.ncbi.nlm.nih.gov/8507387",1993,"primary"),
 ("cgp","Competitive NMDA receptor antagonist CGP 40116 substitutes for the discriminative stimulus effects of ethanol","https://www.sciencedirect.com/science/article/abs/pii/S0014299996006589",1996,"primary"),
 ("disc_review","Discriminative stimulus effects of ethanol: neuropharmacological characterization","https://www.sciencedirect.com/science/article/abs/pii/S0741832998000354",1998,"review"),
 ("disc_regions","Discriminative stimulus effects of ethanol are mediated by NMDA and GABA-A receptors in specific limbic brain regions","https://link.springer.com/article/10.1007/s002130050694",1998,"primary"),
 ("gabox_disc","Comparing discriminative stimulus effects of modulators of GABA-A receptors containing alpha4-delta subunits with those of gaboxadol","https://pmc.ncbi.nlm.nih.gov/articles/PMC5054722/",2016,"primary"),
 ("glun2","Regulation of NMDA glutamate receptor functions by the GluN2 subunits","https://onlinelibrary.wiley.com/doi/10.1111/jnc.14970",2020,"review"),
 ("glun2d_stn","NMDA receptors containing the GluN2D subunit control neuronal function in the subthalamic nucleus","https://pmc.ncbi.nlm.nih.gov/articles/PMC4666920/",2015,"primary"),
 ("glun2d_hipp","GluN2D-containing NMDA receptors mediate synaptic currents in hippocampal interneurons","https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4373385/",2015,"primary"),
 ("pbc_alpha","Developmental changes in expression of GABA-A receptor subunits alpha1, alpha2, alpha3 in the rat pre-Botzinger complex","https://journals.physiology.org/doi/full/10.1152/japplphysiol.01264.2003",2004,"primary"),
 ("pbc_delta","delta-Subunit containing GABA-A receptors modulate respiratory networks","https://www.nature.com/articles/s41598-017-17379-x",2017,"primary"),
 ("pbc_eps","Increased GABA-A receptor epsilon-subunit expression on ventral respiratory column neurons protects breathing during pregnancy","https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0030608",2012,"primary"),
 ("pbc_a4","Respiratory and behavioral dysfunction following loss of the GABA-A receptor alpha4 subunit","https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3607152/",2013,"primary"),
 ("pbc_model","Modeling effects of variable preBotzinger complex network topology and cellular properties on opioid-induced respiratory depression","https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10921262/",2024,"primary"),
 ("cpg_rat","Neuromechanical model of rat hindlimb walking with two-layer CPGs","https://pmc.ncbi.nlm.nih.gov/articles/PMC6477610/",2019,"primary"),
 ("simtk_rat","SimTK musculoskeletal model of the rat hindlimb (38 muscles, hip/knee/ankle, 3D)","https://simtk.org/projects/rat_hlimb_model",None,"resource"),
 ("dpen","Pre-clinical studies with D-penicillamine as a novel pharmacological strategy to treat alcoholism","https://pubmed.ncbi.nlm.nih.gov/28326026/",2017,"review"),
 ("salsolinol","Key role of salsolinol in ethanol actions on dopamine neuronal activity of the posterior VTA","https://onlinelibrary.wiley.com/doi/10.1111/adb.12097",2015,"primary"),
 ("salsolinol_mu","Salsolinol stimulates dopamine neurons in pVTA indirectly by activating mu-opioid receptors","https://www.sciencedirect.com/science/article/abs/pii/S0022356524181542",None,"primary"),
 ("salsolinol_R","(R)-Salsolinol stereospecifically induces behavioral sensitization and leads to excessive alcohol intake","https://onlinelibrary.wiley.com/doi/10.1111/adb.12268",2016,"primary"),
 ("acet_crit","Acetaldehyde: deja vu du jour","https://pubmed.ncbi.nlm.nih.gov/15536764",2004,"commentary"),
 ("acet_rev","Role of acetaldehyde in the neurobehavioral effects of ethanol: comprehensive review of animal studies","https://www.sciencedirect.com/science/article/abs/pii/S0301008205000249",2005,"review"),
 ("vta_invivo","Ethanol stimulates the firing rate of nigral dopaminergic neurons in unanesthetized rats","https://www.sciencedirect.com/science/article/abs/pii/0006899384908904",1984,"primary"),
 ("memantine","Memantine preferentially blocks extrasynaptic over synaptic NMDA receptor currents","https://www.jneurosci.org/content/30/33/11246",2010,"primary"),
 ("vrodent","A virtual rodent predicts the structure of neural activity across behaviours (Aldarondo et al.)","https://www.nature.com/articles/s41586-024-07633-4",2024,"primary"),
 ("aperiodic","Reliable and dynamic monitoring of cortical E/I balance using aperiodic 1/f slope","https://www.sciencedirect.com/science/article/abs/pii/S0306452225012175",2025,"review"),
]

S += [
 ("matsuoka","Matsuoka half-centre oscillator (rate-based reciprocal inhibition + adaptation)",None,1985,"method"),
 ("sim_port","THIS PROJECT: muscle port of dm_control rodent hindlimb (scripts/port_muscle.py)",None,2026,"internal"),
 ("sim_cpg","THIS PROJECT: hybrid spinal circuit, Matsuoka RG + spiking PF/Mn (circuitpharm/circuit.py)",None,2026,"internal"),
]
c.executemany("INSERT INTO sources VALUES(?,?,?,?,?)", S)

T = [
 ("BK/Slo1 channel",None,"modulated",10,50,"EC50 ~22 mM",1,"good","Best concentration match of any target. Behavioural relevance partly unestablished (K361 residue study negative).","bk"),
 ("GABA-A extrasynaptic","delta / alpha4-alpha6","potentiated",0,30,"tonic current increase",1,"CONTESTED","Borghese & Harris could not replicate low-conc sensitivity incl. human a4b3d clones. In vivo support exists (Ro15-4513 binding, NAc shell intake).","gaba_conc"),
 ("GABA-A synaptic (benzodiazepine site)","gamma1-gamma3","potentiated",100,None,"only >100 mM, 'if at all'",0,"good","NOT an ethanol target at relevant concentrations. Explains why benzodiazepine-site compounds fail to mimic ethanol on EEG.","gaba_conc"),
 ("Glycine receptor","alpha1, alpha2 (+beta)","potentiated",5,50,"potentiation at 43.2 mM tested",1,"good","KI mice with ethanol-insensitive a1/a2 show 30% shorter LORR. a2b potentiated; a2 alone insensitive. Gbg-dependent.","glyr_ki"),
 ("NMDA receptor (recombinant)","GluN2B","inhibited",5,50,"6.5+/-0.8% at 20 mM",1,"good","Quantitatively weak at social doses. Concentration-dependent 5-50 mM.","nmda_drivelimit"),
 ("NMDA receptor (synaptic, BNST)","GluN2B","inhibited",None,None,"~50% of NMDA transmission",1,"good","Abolished in GluN2B KO, retained in GluN2A KO. Amplification gradient vs recombinant 6.5%.","nmda_bnst"),
 ("5-HT3 receptor","5-HT3A","potentiated",25,100,"potentiation",0,"moderate","Marginal at 20 mM. Ethanol POTENTIATES, so antagonists (ondansetron) are the wrong direction for mimicry.","gaba_decade"),
 ("GIRK2 channel","GIRK2","enhanced",10,50,"enhanced at intoxicating conc",1,"moderate","Recombinant homomeric and heteromeric channels.","girk"),
 ("Adenosine A1 (indirect via ENT1)","A1","increased signalling",None,None,"motor incoordination",1,"good","Indirect: ethanol inhibits ENT1 -> adenosine accumulates. A1 blockade / caffeine attenuates incoordination. A2A handles reward.","aden"),
 ("nAChR",None,"potentiated",100,None,"only >=100 mM",0,"good","NOT relevant at social doses. Drop from panel.","gaba_conc"),
]
c.executemany("INSERT INTO ethanol_targets VALUES(?,?,?,?,?,?,?,?,?,?)", T)

D = [
 ("ethanol","MK-801 (dizocilpine)","various",None,"yes","Complete generalization. Noncompetitive NMDA antagonist.","mk801"),
 ("ethanol","CGP 40116","various",88.0,"yes","Competitive NMDA antagonist, active d-stereoisomer of CGP 37849.","cgp"),
 ("ethanol","diazepam / chlordiazepoxide / pentobarbital","various",None,"yes","GABA-A positive modulators substitute.","disc_review"),
 ("ethanol","neuroactive steroids","various",None,"yes","Ethanol MORE likely to substitute for neurosteroids than for benzodiazepines.","disc_review"),
 ("ethanol","GHB","various",None,"yes","Complete substitution reported.","disc_review"),
 ("ethanol","gaboxadol (THIP)","various",None,"NO","Ethanol-trained rats respond on VEHICLE lever. delta-selective direct agonist does NOT substitute.","gabox_disc"),
 ("midazolam","gaboxadol","various",None,"NO","Also fails in midazolam- and pregnanolone-trained rats.","gabox_disc"),
 ("benzodiazepine / barbiturate","ethanol","various",None,"NO","ASYMMETRY: ethanol generalizes to benzos/barbs, but not the reverse. Ethanol stimulus CONTAINS theirs plus more.","disc_review"),
 ("chlordiazepoxide 5.0 + dizocilpine 0.08 mg/kg","ethanol 3.0 g/kg","CDP 5.0 + DIZ 0.08 mg/kg",76.0,"yes","MIXTURE: GABA-salient ratio works.","mixture"),
 ("pentobarbitone 8.0 + dizocilpine 0.08 mg/kg","ethanol 3.0 g/kg","PB 8.0 + DIZ 0.08 mg/kg",33.0,"no","MIXTURE FAILS: NMDA-weighted ratio.","mixture"),
 ("pentobarbitone 12 + dizocilpine 0.04 mg/kg","ethanol 3.0 g/kg","PB 12 + DIZ 0.04 mg/kg",75.0,"yes","MIXTURE RESCUED by raising GABA and lowering NMDA.","mixture"),
 ("ethanol","MK-801 intra-NAc core or CA1","intracranial",None,"yes","FULL substitution. Localizes the NMDA component.","disc_regions"),
 ("ethanol","MK-801 intra-prelimbic cortex or amygdala","intracranial",None,"NO","Failed to substitute. Region-specific.","disc_regions"),
]
c.executemany("INSERT INTO discrimination VALUES(?,?,?,?,?,?,?)", D)

E = [
 ("bec1982","acute, Sprague-Dawley females",1,"Sprague-Dawley",0,"cortical","0-4, 8-13 Hz","both (time course)",0,
  "Delta 0-4 Hz increase NOT correlated with BEC and persisted after BEC=0. 8-13 Hz dropped then recovered, significant NEGATIVE linear correlation with BEC. => 8-13 Hz is the only BEC-locked spectral observable."),
 ("subtype2025","EtOH 1-3 g/kg IP; alprazolam/zolpidem/MP-III-022/gaboxadol 1-10, KRM-II-81 3-30 mg/kg",0,"Sprague-Dawley",1,
  "centro-frontal, occipital ref","d1-5.5 t5.5-8.5 a8.5-12.5 b12.5-30 g30-50","no",0,
  "Ethanol uniquely increased delta. Alprazolam+zolpidem increased BETA. Gaboxadol: no spectral effect at all (likely time-averaging artifact, see thip_rat). 'none of the compounds tested induced overall similar effects compared to ethanol'. NB delta band to 5.5 Hz overlaps others' theta."),
 ("thip_rat","THIP 2/4/6 mg/kg IP",0,"rat",1,"frontal+parietal","SWA 0.75-4 Hz","no",0,
  "Dramatic SWA increase at 4 and 6 mg/kg but NOT 2 mg/kg; massive power increase <6 Hz in waking AND NREM; dissipates within 60-120 min. CONTRADICTS subtype2025 null -> that study's long averaging window likely diluted a transient effect."),
 ("theta1977","2-4 g/kg IP (0.5-1 g/kg no effect)",0,"rat",0,"hippocampus","theta peak","no",0,
  "Theta slowing only at 2-4 g/kg; peak shifts to 4.50-5.75 Hz. BLIND in the social-drinking range. Burst-resolved analyses give opposite sign to continuous spectra."),
 ("ethbenzo1992","human",0,"human",1,"scalp","power spectra","no",0,
  "Ethanol's sleep EEG spectral effects DIFFER from benzodiazepine receptor agonists. Converges with subtype2025 across species."),
]
c.executemany("INSERT INTO eeg_studies VALUES(?,?,?,?,?,?,?,?,?,?)", E)

X = [
 ("NMDA","GluN2B","hippocampus, cortex, striatum","high","FAVOURABLE: overlaps the regions where intra-cranial MK-801 fully substitutes for ethanol (NAc core, CA1). Also ethanol's own NMDA mechanism.","glun2"),
 ("NMDA","GluN2D","brainstem, substantia nigra, diencephalon, cerebellum, spinal cord","high","FAVOURABLE: brainstem respiratory circuits are GluN2D-dominant, so GluN2B-selective antagonism should spare them.","glun2d_stn"),
 ("NMDA","GluN2D","forebrain PV+/SST+ interneurons","restricted","Predicts GluN2B-selective antagonists SPARE forebrain PV interneurons -> no disinhibition -> non-psychotomimetic. TESTABLE in circuit model.","glun2d_hipp"),
 ("GABA-A","alpha1","preBotzinger complex (adult)","predominant","UNFAVOURABLE: no respiratory-sparing window at alpha1.","pbc_alpha"),
 ("GABA-A","alpha2/alpha3","preBotzinger complex (specific neurons)","present","UNFAVOURABLE.","pbc_alpha"),
 ("GABA-A","delta (extrasynaptic)","respiratory network","present","UNFAVOURABLE: delta-containing receptors modulate respiratory networks.","pbc_delta"),
 ("GABA-A","alpha4","respiratory network","present","UNFAVOURABLE: alpha4 loss produces respiratory dysfunction.","pbc_a4"),
 ("GABA-A","epsilon","ventral respiratory column, NK1R+ preBotC rhythm neurons","enriched","ONLY POTENTIAL GABA-ARM LEVER: a PAM inactive at epsilon-containing receptors would be respiratory-sparing. UNVALIDATED IDEA - needs the epsilon pharmacology read properly.","pbc_eps"),
]
c.executemany("INSERT INTO subunit_expression VALUES(?,?,?,?,?,?)", X)

F = [
 (1,"EEG spectral signature is ANTI-CORRELATED with subjective substitution across the three best-characterized classes: gaboxadol matches ethanol's delta signature but does not substitute; benzodiazepines diverge on EEG (beta) but do substitute; NMDA antagonists diverge (gamma) but do substitute.","EEG vs drug discrimination","high",
  "KILLS the neural-signature-matching strategy. Drug discrimination is the better subjective proxy but cannot be simulated (requires learning) - so take it as an EMPIRICAL CONSTRAINT instead.","gabox_disc"),
 (2,"Ethanol's discriminative stimulus is a COMPOUND stimulus with a GABA-A-positive-modulator component and an NMDA-antagonist component.","drug discrimination","high",
  "Answers 'do we need new chemicals' = NO. The mechanism pair is known and mixtures have been tested.","disc_review"),
 (3,"Ethanol substituted only for mixtures in which the GABA-A component had EQUAL OR GREATER salience than the NMDA-antagonist component.","drug discrimination","high",
  "BINDING CONSTRAINT. Conflicts with the safety optimum: NMDA has the regional selectivity window (GluN2B/GluN2D) while GABA-A has none, which would favour NMDA-heavy. The subjective constraint forbids that. THIS TENSION IS THE PROJECT.","mixture"),
 (4,"Acetaldehyde-blocking manipulations affect preference, acquisition of self-administration, relapse-like drinking, conditioned place preference and sensitization - all REINFORCEMENT measures - while brain acetaldehyde during ordinary intoxication is low micromolar and critics argue too low for direct effect.","reinforcement vs acute state","moderate",
  "RESOLVES the acetaldehyde threat to the project premise: its role is 'wanting', not 'feeling', and 'wanting' is already out of scope. A substitute is predicted to be LESS reinforcing - a feature. Caveat: D-penicillamine also blocked acute behavioural stimulation.","dpen"),
 (5,"Ethanol's pVTA dopamine excitation requires local catalase-mediated acetaldehyde production, condensing with dopamine to salsolinol, which acts via mu-opioid receptors.","VTA dopamine firing","moderate",
  "DEMOTES VTA dopamine as a matching target. A non-metabolized substitute will not reproduce it without mu-opioid agonism, which reintroduces respiratory depression and dependence. DECISION TAKEN: accept divergence on the reward axis.","salsolinol"),
 (6,"(R)-salsolinol stereospecifically induces behavioral sensitization; salsolinol forms as racemic R/S by non-enzymatic Pictet-Spengler condensation.","sensitization","moderate",
  "Parked side-quest: a chiral pair with activity confined to one enantiomer at the centre of alcohol's reward mechanism. Not on the critical path.","salsolinol_R"),
 (7,"Benzodiazepine-site GABA-A receptors (gamma-containing) are not ethanol targets below ~100 mM, while 20 mM is strong intoxication and 100 mM is lethal.","receptor pharmacology","high",
  "Explains the subtype2025 negative result mundanely: it tested four benzodiazepine-site compounds against a receptor ethanol barely engages. Do NOT over-read that study.","gaba_conc"),
 (8,"Combination effects are non-additive and endpoint-specific (triazolam+pregnanolone supra-additive for sedation yet ATTENUATED ataxia; NR2B antagonist + morphine no rotarod potentiation; CPP alone causes ataxia).","motor / sedation","moderate",
  "JUSTIFIES the simulation: single-agent dose-response does not predict combination margins, and the ratio space is too large to test empirically.","disc_review"),
 (9,"Both the preBotzinger respiratory rhythm generator and the spinal locomotor CPG depend on glutamatergic excitation and GABA/glycine reciprocal inhibition - the same two mechanisms being modulated for the wanted effect.","architecture","high",
  "The wanted effect and both unwanted effects run through the same machinery. Margin may be intrinsically narrow; regional/subunit selectivity is the only thing that can open it.","pbc_model"),
 (10,"Isoflurane shows BOTH GABA-A-positive-modulator-like and NMDA-antagonist-like discriminative stimulus effects.","drug discrimination","moderate",
  "Proof of concept that a single molecule can carry both components - a mixture is not strictly required.","disc_review"),
 (11,"Reference concentrations: 1 drink ~5 mM (already alters mood/coordination); 20 mM strong intoxication; 100 mM lethal for most.","concentration anchors","high",
  "Design window is 5-20 mM. Any target requiring >50 mM is irrelevant.","gaba_conc"),
 (12,"No published ethanol EEG study uses periodic/aperiodic (1/f) spectral decomposition, though the aperiodic exponent is validated as a dose-dependent E/I-balance proxy.","EEG methods","high",
  "Cheapest high-value open item: reanalysis of existing recordings, no new experiment. Would also resolve relative-power artifacts that generate several literature contradictions.","aperiodic"),
]

F += [
 (13,"In the spinal circuit model, a 60% GluN2B-selective NMDA block depresses motoneuron output by only 7% in GluN2D-dominant spinal tissue (2B fraction 0.15) but 34% with a forebrain-like subunit profile (2B 0.7); a non-selective block at the same 60% depresses spinal output 43%.","simulated Mn firing rate","model-derived",
  "QUANTIFIES the central margin hypothesis in-model. GluN2B selectivity converts a 43% motor cost into a 7% one at matched receptor block. NOT validated against the four motor anchors yet.","sim_cpg"),
 (14,"The GABA-A PAM + non-selective NMDA antagonist combination is SUPRA-ADDITIVE on motor depression (predicted 28.5% of control from single agents, observed 17%), whereas the GABA-A PAM + GluN2B-selective combination is approximately multiplicative (predicted 46.5%, observed 42%).","simulated Mn firing rate","model-derived",
  "Directly engages the kill criterion about supra-additivity. The non-selective route carries a combination hazard the selective route escapes. Also confirms the model's reason for existing: single-agent data does not predict combination margins.","sim_cpg"),
 (15,"Glycine-receptor potentiation is the ONLY one of the three mechanisms that changes locomotor PERIOD (920 -> 1336 ms); GABA-A PAM and NMDA antagonism change burst amplitude and duty cycle while leaving the cycle period untouched.","simulated locomotor period","model-derived",
  "FALSIFIABLE PREDICTION, testable against existing rodent gait data: ethanol (which potentiates GlyR) should lengthen step-cycle duration, whereas a GABA-A PAM + NMDA antagonist mixture should reduce burst amplitude WITHOUT slowing the cycle. This is the cheapest external validation available.","sim_cpg"),
 (16,"A GABA-A efficacy ceiling limits but does not prevent severe motor depression: capping the PAM at 2.5x still leaves motoneuron output at 17% of control.","simulated Mn firing rate","model-derived",
  "Consistent with the earlier conclusion that the GABA arm has NO selectivity window and only a ceiling. The ceiling is necessary for overdose safety but is not sufficient for a motor margin.","sim_cpg"),
 (17,"An exhaustive 288-point parameter search showed an all-spiking LIF rhythm generator CANNOT produce a physiological locomotor rhythm: near-target period only with no alternation (corr -0.17), good alternation only at 81-111 ms, duty cycle never above 0.24. Reciprocal inhibition weight from 8 to 40 changed nothing.","locomotor rhythm","model-derived",
  "STRUCTURAL, not parametric: a LIF resets V every spike so no plateau potential forms, the half-centre is never bistable on the fast timescale, and adaptation can never set the period. Do not retry tuning the spiking RG -- use the Matsuoka RG in circuitpharm/rg.py.","sim_cpg"),
]
c.executemany("INSERT INTO findings VALUES(?,?,?,?,?,?)", F)
R = [
 ("control", 35.8, 100., 920., 0.33, -0.53, "yes", "operating point: rg_gain=900, Matsuoka tau=50/tau_a=400"),
 ("GABA-A PAM 1.5x", 27.4, 77., 920., 0.26, -0.39, "yes", "dose-dependent amplitude depression, period unchanged"),
 ("GABA-A PAM 2.0x", 17.7, 50., 920., 0.22, -0.32, "yes", ""),
 ("GABA-A PAM 3.0x (capped at 2.5)", 6.0, 17., 921., 0.16, -0.21, "yes", "efficacy ceiling LIMITS but does not PREVENT severe depression"),
 ("NMDA 30% non-selective", 27.8, 78., 920., 0.27, -0.41, "yes", ""),
 ("NMDA 60% non-selective", 20.4, 57., 920., 0.21, -0.31, "yes", "43% motor depression"),
 ("NMDA 60% GluN2B-sel, spinal (2B frac 0.15)", 33.4, 93., 920., 0.31, -0.49, "yes", "SELECTIVITY WINDOW: only 7% depression at the same 60% block"),
 ("NMDA 60% GluN2B-sel, forebrain-like (2B 0.7)", 23.7, 66., 920., 0.23, -0.34, "yes", "same drug, forebrain subunit profile -> 34% effect"),
 ("ethanol-like (GABA+NMDA+GlyR)", 17.5, 49., 1336., 0.17, -0.22, "yes", "ONLY condition that changed PERIOD (920->1336 ms), via GlyR potentiation"),
 ("candidate: PAM 2x + GluN2B-sel 60%", 15.1, 42., 920., 0.19, -0.28, "yes", "~multiplicative (predicted 46.5%, observed 42%)"),
 ("candidate: PAM 2x + non-selective 60%", 6.2, 17., 920., 0.14, -0.19, "yes", "SUPRA-ADDITIVE (predicted 28.5%, observed 17%) -- hazard signal"),
]
c.executemany("INSERT INTO sim_results VALUES(?,?,?,?,?,?,?,?)", R)


M = [
 ("rodent.xml path",".venv/lib/python3.11/site-packages/dm_control/locomotion/walkers/assets/rodent.xml","file exists"),
 ("nq (joint coords)","67","mujoco.MjModel introspection"),
 ("nv (dof)","67","mujoco.MjModel introspection"),
 ("nu (actuators)","38","mujoco.MjModel introspection"),
 ("nbody","66","mujoco.MjModel introspection"),
 ("total mass","0.3378 kg","body_subtreemass[0]"),
 ("timestep","0.002 s","m.opt.timestep"),
 ("muscle actuators","NONE - all 38 are gaintype 0 (torque motors)","actuator_gaintype histogram {0:38}"),
 ("physics speed","~41,700 steps/s single core = ~83x realtime on Apple Silicon","2000-step timing loop"),
 ("tasks verified","rodent_escape_bowl, rodent_run_gaps, rodent_maze_forage, rodent_two_touch","basic_rodent_2020, 200 steps each"),
 ("task speed with vision","160-310 env-steps/s (egocentric camera dominates cost)","timing loop"),
 ("proprioceptive obs available","joints_pos/vel, tendons_pos/vel, sensors_force/torque/touch, accelerometer, gyro, velocimeter, end_effectors_pos","observation spec"),
 ("CRITICAL GAP","no muscles, no spindle/Golgi afferents, no reflex loop -> reflex phenomena (hyperreflexia, clonus, tremor) cannot emerge for the right reasons","verified above"),
]

M += [
 ("PORTED MODEL","models/rodent_muscle.xml -- 24 antagonist muscle actuators replace 12 hindlimb position servos; nu 38->50","mujoco compile + force sweep"),
 ("why ported","stock hindlimb actuators are POSITION SERVOS (dyntype=filter, biastype=affine): they (a) collapse flexor+extensor into one signed command so glycinergic reciprocal inhibition is unrepresentable, (b) compensate for drug-induced weakness, masking it, (c) have no force-length/velocity for spindles","XML inspection"),
 ("muscle verification","both antagonists pull across the whole joint range, each strongest when lengthened, equal and balanced at joint midpoint, co-contraction gives expected stiffness","qpos sweep, knee_L"),
 ("PORT GOTCHA","stock <default><general gainprm=0.01/> LEAKS into <muscle> (a general shortcut) and overwrites muscle range[0], leaving a dead zone at neutral posture. MUST set range explicitly (0.7 1.3).","force sweep debugging"),
 ("CPG operating point","Matsuoka RG tau=50 tau_a=400 beta=2.5 w=2.2; rg_gain=900; gaba_tonic=2.0 nS; Mn tref=8 ms","scripts/sweep_gain.py"),
 ("CPG verified output","period 920 ms, Mn 35.8 Hz, peak 123 Hz, duty 0.33, flex/ext corr -0.53 -- inside the rat physiological box","scripts/drug_panel.py control row"),
 ("Matsuoka period law","period ~ 2*pi*sqrt(tau*tau_a); predicted 688/889/1192 vs observed 704/920/1299 ms","parameter sweep"),
 ("CPG speed","~0.7x realtime single-threaded; MuJoCo physics is ~41,700 steps/s, so the brain is ~140x slower than the body. Batch across drug conditions (NumPy dim) + multiprocess for the ratio sweep.","timing"),
 ("hardware","Apple M5, 10 cores, 32 GB. 10-way parallel sweep = 938% CPU, only 0.52 GB RSS total. Hardware is NOT the constraint; single-threaded Python is.","ps/memory_pressure during sweep"),
 ("OPEN-LOOP CAVEAT","all drug results so far are OPEN-LOOP: Ia afferents are zero, so IaIn reciprocal inhibition is silent and flex/ext alternation (-0.53) is weaker than the RG's (-0.90). Closing the MuJoCo loop should sharpen it.","drug_panel.py"),
]
c.executemany("INSERT INTO model_facts VALUES(?,?,?)", M)

c.commit()
for t in ["sources","ethanol_targets","discrimination","eeg_studies","subunit_expression","findings","sim_results","model_facts"]:
    print(f"{t:20s} {c.execute(f'SELECT COUNT(*) FROM {t}').fetchone()[0]:3d} rows")
c.close(); print(f"\n-> {DB}")
