import os
import numpy as np
from multiprocessing import Pool
from circuitpharm.circuit import SpinalCircuit
from circuitpharm.cpg import burst_metrics

def run(args):
    gain, gt = args
    c = SpinalCircuit(rg_gain=gain, gaba_tonic=gt, seed=1); dt=0.1
    for i in range(int(9000/dt)):
        c.step(dt)
        if i%10==0: c.record()
    A=c.arrays(); m=A["t"]>3000
    # smooth before burst detection: the 20 ms population low-pass leaves sub-peaks
    # that get miscounted as bursts on a ~900 ms locomotor cycle
    k=np.ones(60)/60.0
    sm=lambda v: np.convolve(v,k,mode="same")
    mf,me = sm(A["Mn_F"][m]), sm(A["Mn_E"][m])
    if mf.std()<1e-9 or me.std()<1e-9: return None
    bm=burst_metrics(A["t"][m],mf)
    return dict(gain=gain, gt=gt, pf=float(A["PF_F"][m].mean()),
                mnf=float(mf.mean()), mne=float(me.mean()),
                mnmax=float(mf.max()), rc=float(A["Rc_F"][m].mean()),
                per=bm["period"], duty=bm["duty"], nb=bm["n_bursts"],
                cm=float(np.corrcoef(mf,me)[0,1]))

GRID=[(g,gt) for g in (250.,400.,550.,700.,900.) for gt in (2.0,)]
if __name__=="__main__":
    with Pool(os.cpu_count()) as p: res=[r for r in p.map(run,GRID) if r]
    print(f"{'rg_gain':>8} {'PF Hz':>6} {'Mnf':>6} {'Mne':>6} {'Mnmax':>6} {'Rc':>6} "
          f"{'per':>6} {'duty':>5} {'nb':>3} {'corr':>6}")
    for r in sorted(res,key=lambda x:x['gain']):
        d=r['duty'] if r['duty']==r['duty'] else float('nan')
        pe=r['per'] if r['per']==r['per'] else float('nan')
        print(f"{r['gain']:8.0f} {r['pf']:6.1f} {r['mnf']:6.1f} {r['mne']:6.1f} "
              f"{r['mnmax']:6.1f} {r['rc']:6.1f} {pe:6.0f} {d:5.2f} {r['nb']:3d} {r['cm']:+6.2f}")
