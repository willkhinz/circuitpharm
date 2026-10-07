"""Run the subtype-resolved GABA profiles through the actual circuits at MATCHED
subjective effect, combined with the best NMDA arm.

Each GABA compound is dosed so its forebrain subjective index equals a common target.
Its regional sensitivities (preBotC, spinal) then follow from its subunit selectivity,
with the non-selective benzodiazepine anchored at the calibrated value of 0.15.
"""
import sys, os; sys.path.insert(0,".")
import numpy as np
from multiprocessing import Pool
from spinal.resp import PreBotC, resp_metrics
from spinal.circuit import SpinalCircuit
from spinal.cpg import Drug
from spinal.subtypes import PROFILES
import scripts.reflex as R

RESP_OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0, w=dict(ee_ampa=0.45, ee_nmda=0.2475))
IA_SCALE, N_SEED = 0.30, 4
SUBJ_TARGET = 0.30          # common forebrain subjective index
NMDA_FB = 0.20              # 20% forebrain NMDA reduction, GluN2B-selective
PERIPH_2B, FB_2B = 0.15, 0.70

ARMS = ["nonselective_bz","neurosteroid","mp_iii_022","hz_166","l_838417","tpa023","ideal_a5","zolpidem"]

def dose_for(key, target):
    """PAM gain needed for this compound's subjective index to reach `target`."""
    p = PROFILES[key]; si = p.subjective_index()
    if si <= 1e-9: return None
    # subjective index scales with (gain-1); unit gain increment gives si
    gain = 1.0 + target / si
    return min(gain, p.ceiling)

def job(a):
    key, kind, seed, with_nmda = a
    p = PROFILES[key]
    gain = dose_for(key, SUBJ_TARGET)
    if gain is None: return None
    kw = dict(gaba_a_gain=gain, gaba_a_tau=1.0+0.6*(gain-1.0), gaba_a_efficacy_cap=p.ceiling)
    if with_nmda:
        kw.update(nmda_block=min(NMDA_FB/FB_2B,1.0), glun2b_selectivity=1.0,
                  glun2b_fraction=PERIPH_2B)
    d = Drug(**kw)
    if kind=="resp":
        b=PreBotC(drug=d, gaba_sens=p.regional_sens("prebotc"), seed=seed, **RESP_OP)
        for i in range(int(14000/0.1)):
            b.step(0.1)
            if i%10==0: b.record()
        A=b.arrays(); m=A["t"]>4000
        r=resp_metrics(A["t"][m],A["Out"][m],A["Exc"][m])
        return key,kind,with_nmda,r["mean"],float(r["alive"])
    orig=dict(SpinalCircuit.W)
    try:
        SpinalCircuit.W=dict(orig)
        for rec in ("ampa","nmda"):
            SpinalCircuit.W[("Ia","Mn",rec)]=orig[("Ia","Mn",rec)]*IA_SCALE
        # patch the spinal sensitivity for this compound
        import scripts.reflex as RR
        from spinal.plant import JointPlant
        pl=JointPlant(RR.XML, drug=d, rg_gain=0.0, seed=seed, gaba_sens=p.regional_sens("spinal"))
        traj=RR.trajectory(0.1,0.45)
        for (t,q,v) in traj: pl.step(RR.DT, impose=(q,v))
        A=pl.arrays(); tt=A["t"]/1000.
        base=(tt>0.15)&(tt<RR.PRE); dyn=(tt>RR.PRE)&(tt<RR.PRE+RR.RAMP+0.03)
        g=(A["mn_E"][dyn].max()-A["mn_E"][base].mean())/max(1e-6,A["ia_E"][dyn].max()-A["ia_E"][base].mean())
        return key,kind,with_nmda,g,1.0
    finally:
        SpinalCircuit.W=orig

if __name__=="__main__":
    jobs=[(k,kind,s,wn) for k in ARMS for kind in ("resp","reflex")
          for s in range(N_SEED) for wn in (False,True)]
    jobs=[("nonselective_bz","resp",s,False) for s in range(N_SEED)]+jobs  # control base
    print(f"{len(jobs)} sims on {os.cpu_count()} cores")
    with Pool(os.cpu_count()) as p: out=[r for r in p.map(job,jobs) if r]
    A={}
    for key,kind,wn,v,al in out: A.setdefault((key,kind,wn),[]).append((v,al))
    # controls: drug-free
    ctrl_jobs=[("ideal_a5","resp",s,False) for s in range(1)]
    print()
    print(f"{'GABA arm':<46} {'PAM':>5} {'resp sens':>10} | {'GABA only':>20} | {'+ GluN2B 20%':>20}")
    print(f"{'':<46} {'dose':>5} {'(preBotC)':>10} | {'vent%':>9} {'reflex%':>10} | {'vent%':>9} {'reflex%':>10}")
    print("-"*118)
    # reference = non-selective BZ without NMDA, normalised to itself for readability
    base_v=np.mean([x[0] for x in A[("nonselective_bz","resp",False)]])
    base_r=np.mean([x[0] for x in A[("nonselective_bz","reflex",False)]])
    for k in ARMS:
        p=PROFILES[k]; g=dose_for(k,SUBJ_TARGET)
        if g is None:
            print(f"{p.name:<46} {'n/a':>5}  (no subjective route)"); continue
        row=f"{p.name:<46} {g:5.2f} {p.regional_sens('prebotc'):10.4f} |"
        for wn in (False,True):
            v=np.mean([x[0] for x in A[(k,'resp',wn)]]); r=np.mean([x[0] for x in A[(k,'reflex',wn)]])
            row+=f" {100*v/base_v:9.0f} {100*r/base_r:10.0f} |"
        print(row)
    print("\n(all values relative to the NON-SELECTIVE BENZODIAZEPINE arm at the same")
    print(" subjective target = 100%; higher is better, i.e. more function preserved)")
