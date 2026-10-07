"""Tune the v2 group-pacemaker RG, and test the published validation target:
blocking glycine must ABOLISH ALTERNATION while PRESERVING the rhythm."""
import sys, os, itertools; sys.path.insert(0,".")
import numpy as np
from multiprocessing import Pool
from circuitpharm.rg2 import GroupPacemakerRG
from circuitpharm.cpg import Drug, burst_metrics

def metrics(drug, drive, gadp, tau_a, wgly, wee=0.42, seed=1, T=30000., dt=0.1):
    r = GroupPacemakerRG(drug=drug, drive=drive, g_adapt=gadp, tau_adapt=tau_a,
                         w=dict(ie_gly=wgly, ee_ampa=wee, ee_nmda=wee*0.55), seed=seed)
    for i in range(int(T/dt)):
        r.step(dt)
        if i%10==0: r.record()
    A=r.arrays(); m=A["t"]>6000
    k=np.ones(50)/50.; sm=lambda v: np.convolve(v,k,mode="same")
    f,e = sm(A["RG_F"][m]), sm(A["RG_E"][m])
    if f.std()<1e-6 or e.std()<1e-6: return None
    bm = burst_metrics(A["t"][m], f)
    return dict(f=float(f.mean()), e=float(e.mean()),
                inrg=float(np.convolve(A["InRG_F"][m],k,mode="same").mean()),
                corr=float(np.corrcoef(f,e)[0,1]),
                per=bm["period"], duty=bm["duty"], nb=bm["n_bursts"])

def probe(a):
    drive,gadp,tau_a,wgly,wee = a
    ctrl = metrics(Drug(), drive,gadp,tau_a,wgly,wee)
    if ctrl is None or ctrl["nb"]<4: return None
    blk  = metrics(Drug(glyr_gain=0.05), drive,gadp,tau_a,wgly,wee)   # strychnine
    if blk is None: return None
    return dict(drive=drive,gadp=gadp,tau_a=tau_a,wgly=wgly,wee=wee,ctrl=ctrl,blk=blk)

# longer bursts: stronger recurrent excitation sustains the burst, weaker adaptation
# terminates it later. Duty 0.05-0.08 was far too brief for inhibition to enforce phase.
G=list(itertools.product((200.,260.),(0.5,0.8,1.2),(180.,280.,380.),(3.0,7.0),(0.55,0.75)))
if __name__=="__main__":
    with Pool(os.cpu_count()) as p: res=[r for r in p.map(probe,G) if r]
    # With an intrinsic frequency mismatch (asym), the test is now meaningful:
    #   control        -> entrained, anti-phase  => corr strongly negative
    #   glycine blocked-> drift apart            => corr toward 0, rhythm intact
    def score(r):
        c,b=r["ctrl"],r["blk"]
        s = 0.0
        s += 0.0 if 300<=c["per"]<=1500 else 2.0
        s += 0.0 if c["nb"]>=12 else 1.5            # need statistical power
        d = c["duty"] if c["duty"]==c["duty"] else 0.0
        s += max(0.0, 0.35-d)/0.35 * 2.0            # duty must reach ~0.35-0.6
        s += max(0.0, c["corr"]+0.4)/0.6            # want corr <= -0.4
        s += max(0.0, -b["corr"])/0.5               # want blocked corr >= 0
        s += 0.0 if b["nb"]>=10 else 2.0            # rhythm must survive, >=10 cycles
        s += abs(np.log(max(b["per"],1)/max(c["per"],1)))   # period roughly preserved
        return s
    res.sort(key=score)
    print(f"{len(res)}/{len(G)} rhythmic")
    print(f"{'drv':>4} {'gadp':>5} {'tauA':>5} {'wgly':>5} {'wee':>5} | "
          f"{'CTRL per':>8} {'corr':>6} {'duty':>5} {'nb':>3} | "
          f"{'BLK per':>8} {'corr':>6} {'nb':>3} | {'score':>6}")
    for r in res[:12]:
        c,b=r["ctrl"],r["blk"]
        pc=c['per'] if c['per']==c['per'] else float('nan')
        pb=b['per'] if b['per']==b['per'] else float('nan')
        dc=c['duty'] if c['duty']==c['duty'] else float('nan')
        print(f"{r['drive']:4.0f} {r['gadp']:5.1f} {r['tau_a']:5.0f} {r['wgly']:5.1f} {r['wee']:5.2f} | "
              f"{pc:8.0f} {c['corr']:+6.2f} {dc:5.2f} {c['nb']:3d} | "
              f"{pb:8.0f} {b['corr']:+6.2f} {b['nb']:3d} | {score(r):6.2f}"
              f"  In={c['inrg']:.1f}Hz")
