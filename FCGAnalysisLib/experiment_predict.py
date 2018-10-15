# 实验参数估计函数库
# Last Update: 2018/10/15 Version 2.1.5
# Lu Yunchao

import numpy as np

from FCGAnalysisLib import mts_analysis


def DeltaKPredict(w, thickness, start, final, load, r, step=1):
    """
    对于CT试件，根据试件几何尺寸和载荷参数，输入的需要计算的起止裂纹长度，按照一定的步长，输出计算得到的dk
    :param w: CT试件几何参数width
    :param thickness: CT试件几何参数
    :param start: 输入的需要计算的裂纹长度起点，单位mm，该长度含notch的长度
    :param final: 输入的需要计算的裂纹长度终点，单位mm，该长度含notch的长度
    :param load: 载荷峰值
    :param r: 应力比
    :param step: 计算的步长，默认为1mm
    :return: dk: 裂纹长度范围内对应的应力强度因子变幅
    """
    a = np.arange(start, final+step, step)
    pmin = load * r
    dk = mts_analysis.DeltaKCalculating(b=thickness, w=w, a=a, pmax=load, pmin=pmin)
    return dk


def CycleIntegrateBySimpson(c, m, dk, cracklength):
    """
    利用辛普森公式由dadn积分得到a-N，首先输入dk，根据Paris公式计算得到dadn，再将dadn按步长积分
    :param c: Paris公式参数c
    :param m: Paris公式参数m
    :param dk: 需要计算的应力强度因子变幅
    :param cracklength: 裂纹长度
    :return: cyclelist: 各dk对应的循环次数列表
    """
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


def PredictPrint(load, a, dk, cycles, dadn, valid, alpha=-1):
    # 输出预测的结果
    if alpha != -1:
        print('********* Max Load =', load, 'N,alpha=', alpha, '. **********')
    else:
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




