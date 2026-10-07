"""Final confirmation: alpha5-selective GABA arm + NMDA arm, 8 seeds, vs references."""
import sys, os; sys.path.insert(0,".")
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.circuit import SpinalCircuit
from circuitpharm.plant import JointPlant
from circuitpharm.cpg import Drug
from circuitpharm.subtypes import PROFILES
import scripts.reflex as R

RESP_OP=dict(drive=170.0,g_adapt=2.5,tau_adapt=400.0,w=dict(ee_ampa=0.45,ee_nmda=0.2475))
IA_SCALE,N_SEED=0.30,8; FB_2B,PERIPH_2B=0.70,0.15

def mk(gkey, subj, nmda_fb, sel):
    if gkey is None: g,p=1.0,PROFILES["nonselective_bz"]
    else:
        p=PROFILES[gkey]; g=1.0+subj/p.subjective_index(); g=min(g,p.ceiling)
    kw=dict(gaba_a_gain=g,gaba_a_tau=1.0+0.6*(g-1.0),gaba_a_efficacy_cap=p.ceiling)
    if nmda_fb>0:
        kw.update(nmda_block=min(nmda_fb/(FB_2B if sel else 1.0),1.0),
                  glun2b_selectivity=1.0 if sel else 0.0, glun2b_fraction=PERIPH_2B)
    return Drug(**kw), (p.regional_sens("prebotc") if gkey else 0.15), \
           (p.regional_sens("spinal") if gkey else 1.0)

CASES=[("drug-free control",None,0.0,0.0,True),
       ("ref: non-selective BZ alone, subj 0.50","nonselective_bz",0.50,0.0,True),
       ("ref: neurosteroid alone, subj 0.50","neurosteroid",0.50,0.0,True),
       ("A5 ALONE: a5-sel PAM, subj 0.50","mp_iii_022",0.50,0.0,True),
       ("CANDIDATE: a5-sel PAM 0.50 + GluN2B 20%","mp_iii_022",0.50,0.20,True),
       ("CANDIDATE: a5-sel PAM 0.50 + GluN2B 35%","mp_iii_022",0.50,0.35,True),
       ("  same but NON-selective NMDA 20%","mp_iii_022",0.50,0.20,False),
       ("old rec: neurosteroid 0.50 + GluN2B 20%","neurosteroid",0.50,0.20,True),
       ("a2/a3 PAM 0.35 + GluN2B 20%","hz_166",0.35,0.20,True)]

def job(a):
    i,kind,seed=a
    lab,gk,subj,nf,sel=CASES[i]
    d,sr,ss=mk(gk,subj,nf,sel)
    if kind=="resp":
        b=PreBotC(drug=d,gaba_sens=sr,seed=seed,**RESP_OP)
        for k in range(int(14000/0.1)):
            b.step(0.1)
            if k%10==0: b.record()
        A=b.arrays(); m=A["t"]>4000
        r=resp_metrics(A["t"][m],A["Out"][m],A["Exc"][m])
        return i,kind,r["mean"],float(r["alive"])
    orig=dict(SpinalCircuit.W)
    try:
        SpinalCircuit.W=dict(orig)
        for rec in ("ampa","nmda"):
            SpinalCircuit.W[("Ia","Mn",rec)]=orig[("Ia","Mn",rec)]*IA_SCALE
        pl=JointPlant(R.XML,drug=d,rg_gain=0.0,seed=seed,gaba_sens=ss)
        for (t,q,v) in R.trajectory(0.1,0.45): pl.step(R.DT,impose=(q,v))
        A=pl.arrays(); tt=A["t"]/1000.
        base=(tt>0.15)&(tt<R.PRE); dyn=(tt>R.PRE)&(tt<R.PRE+R.RAMP+0.03)
        g=(A["mn_E"][dyn].max()-A["mn_E"][base].mean())/max(1e-6,A["ia_E"][dyn].max()-A["ia_E"][base].mean())
        return i,kind,g,1.0
    finally: SpinalCircuit.W=orig

if __name__=="__main__":
    jobs=[(i,k,s) for i in range(len(CASES)) for k in ("resp","reflex") for s in range(N_SEED)]
    with Pool(os.cpu_count()) as p: out=p.map(job,jobs)
    A={}
    for i,k,v,al in out: A.setdefault((i,k),[]).append((v,al))
    cv=np.mean([x[0] for x in A[(0,"resp")]]); cr=np.mean([x[0] for x in A[(0,"reflex")]])
    print(f"drug-free control (n={N_SEED}): ventilation {cv:.2f}, reflex gain {cr:.3f}\n")
    print(f"{'condition':<44} {'ventilation':>16} {'stretch reflex':>18}")
    print("-"*80)
    for i,(lab,*_ ) in enumerate(CASES):
        rv=100*np.array([x[0] for x in A[(i,"resp")]])/cv
        rf=100*np.array([x[0] for x in A[(i,"reflex")]])/cr
        print(f"{lab:<44} {rv.mean():10.0f}±{rv.std():3.0f}% {rf.mean():12.0f}±{rf.std():3.0f}%")
