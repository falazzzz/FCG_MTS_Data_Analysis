# 实验参数估计函数库

import numpy as np
import mts_analysis


def DeltaKPredict(w, thickness, start, final, load, r, step=1):
    a = np.arange(start, final+step, step)
    pmin = load * r
    dk = mts_analysis.DeltaKCalculating(b=thickness, w=w, a=a, pmax=load, pmin=pmin)
    return dk


def CycleIntegrateBySimpson(c, m, dk, cracklength):
    cyclelist = []
    cycle = 0
    def parisintegrate(c, m, dki):
        f = 1/(c*(dki**m))
        return f

    for seq, value in enumerate(cracklength[1:]):
        fa = parisintegrate(c=c, m=m, dki=dk[seq])
        fm = parisintegrate(c=c, m=m, dki=0.5*(dk[seq]+dk[seq+1]))
        fb = parisintegrate(c=c, m=m, dki=dk[seq+1])
        a = cracklength[seq]
        b = value
        cycle_increased = (b - a)/6 * (fa + 4*fm + fb)
        cycle = cycle + cycle_increased
        cyclelist.append(cycle)
    cyclelist.insert(0, 0)
    return np.array(cyclelist)


def PredictPrint(load, a, dk, cycles, dadn, valid):
    print('************* Max Load =', load, 'N. *****************')
    seq = 0
    while seq < len(a):
        print('cracklength = ', str(a[seq]),
              'mm,Valid =', str(valid[seq]),
              ',DeltaK =', str(dk[seq]),
              'Mpa.m0.5, cycles =', str(cycles[seq]),
              ',dadn=', str(dadn[seq]), 'mm/cycle.')
        if valid[seq] == False:
            print('**************************************************')
            break
        seq = seq + 1
    return seq




