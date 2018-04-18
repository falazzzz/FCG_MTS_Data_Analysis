# 高载迟滞分析函数库，包括塑性区计算及各高载模型的计算
# Last Update: 2018/4/11 Version 2.2.1
# Lu Yunchao

import numpy as np
from FCGAnalysisLib import mts_analysis
import matplotlib.pyplot as plt


def PlasticZoneWithFactor(kmax, ys, factor):
    # usage: 计算Irwin及其改进型的塑性区尺寸，表达式r = factor*(kmax/ys)^2，其中factor为塑性区系数
    # input parameter:
    # kmax：应力强度因子/MPa.m0.5
    # ys：屈服强度/GPa
    # factor：塑性区尺寸系数
    # return parameter:
    # r：塑形区尺寸/mm
    r = []
    ys = ys * 1e3
    for kmaxs in kmax:
        r.append((factor * (kmaxs/ys)**2)*1e3)
    return np.array(r)


def WheelerFittingBaseParis(a, dadn, dk, rm, aol, rol, c, m):
    # usage: 由高载为起始点的数据，以Paris公式为基础裂纹扩展速率模型进行Wheeler模型拟合
    # NOTE：数据的起始点必须是高载后！
    # input parameter:
    # a, dadn, dk, rm：FCG基础数据及对应时刻的塑性区尺寸rm
    # aol，rol：高载时刻的裂纹长度aol和塑形区尺寸rol
    # c，m：基础Paris公式参数
    # return parameter:
    # m1：Wheeler迟滞指数
    # b：拟合时图形的截距
    left = a + rm
    right = aol + rol
    if left[-1] < right:
        dadn_fitting = dadn
        dk_fitting = dk
        a_fitting = a
        rm_fitting = rm
        print("Wheeler Region Check: Fitting Region:dk_start=" + str(dk_fitting[0]) + ",to the END.")
    elif left[0] >= right:
        print("No Wheeler Region, Check Input Parameter!")
    else:
        for seq, _ in enumerate(left):
            if left[seq] >= right:
                dadn_fitting = dadn[0:seq]
                dk_fitting = dk[0:seq]
                a_fitting = a[0:seq]
                rm_fitting = rm[0:seq]
                print("Wheeler Region Check: Fitting Region:dk_start=" + str(dk_fitting[0]) + ",dk_final=" + str(
                    dk_fitting[-1]))
                break
    # 根据Wheeler公式筛选迟滞参数适用区域（确定拟合采用的数据）

    X = []
    Y = []
    for seq, _ in enumerate(dadn_fitting):
        y = np.log(dadn_fitting[seq]) - m * np.log(dk_fitting[seq]) - np.log(c)
        Y.append(y)
        x = np.log(rm_fitting[seq]) - np.log(aol + rol - a_fitting[seq])
        X.append(x)
    X = np.array(X).reshape(-1, )
    Y = np.array(Y)
    p = np.polyfit(X, Y, deg=1)
    m1, b = p
    print("Wheeler Model Fitting Result: m1=", m1, ",intercept=", b)
    # Wheeler模型 m1系数拟合
    return m1, b


def DelayFittingBaseParis(a, dadn, dk, cp, rd, aol, rdol, c, m):
    # usage: 由高载为起始点的数据，以Paris公式为基础裂纹扩展速率模型进行Wheeler模型拟合
    # NOTE：数据的起始点必须是高载后！
    # input parameter:
    # a, dadn, dk, rd，cp：FCG基础数据及对应时刻的延迟区尺寸rd，延迟区系数cp
    # aol，rol：高载时刻的裂纹长度aol和延迟区尺寸rdol
    # c，m：基础Paris公式参数
    # return parameter:
    # m_mod：Salvati延迟模型指数
    # b：拟合时图形的截距
    left = a + rd
    right = aol + rdol
    if left[-1] < right:
        dadn_fitting = dadn
        dk_fitting = dk
        a_fitting = a
        rd_fitting = rd
        cp_fitting = cp
        print("Delay Region Check: Fitting Region:dk_start=" + str(dk_fitting[0]) + ",to the END.")
    elif left[0] >= right:
        print("No Delay Region, Check Input Parameter!")
    else:
        for seq, _ in enumerate(left):
            if left[seq] >= right:
                dadn_fitting = dadn[0:seq]
                dk_fitting = dk[0:seq]
                a_fitting = a[0:seq]
                rd_fitting = rd[0:seq]
                cp_fitting = cp[0:seq]
                print("Delay Region Check: Fitting Region:dk_start=" + str(dk_fitting[0]) + ",dk_final=" + str(
                    dk_fitting[-1]))
                break
    # 根据Salvati延迟模型筛选迟滞参数适用区域（确定拟合采用的数据）

    X = []
    Y = []
    for seq, _ in enumerate(dadn_fitting):
        y = np.log(dadn_fitting[seq]) - m * np.log(dk_fitting[seq]) - np.log(c) - np.log(cp_fitting[seq])
        Y.append(y)
        x = np.log(aol + rdol - a_fitting[seq]) - np.log(rd_fitting[seq])
        X.append(x)
    X = np.array(X).reshape(-1, )
    Y = np.array(Y)
    p = np.polyfit(X, Y, deg=1)
    m_mod, b = p
    print("Delay Model Fitting Result: m_mod=", m_mod, ",intercept=", b)
    # Salvati延迟模型 m1系数拟合
    return m_mod, b


def WheelerCalculatingBaseParis(astart, afinal, b, w, ys, pmax, r, aol, rol, plasticzonefactor, m1, c, m):
    # usage: 以裂纹长度a指定计算范围，计算基于Paris公式的Wheeler模型
    # NOTE：数据的起始点必须是高载后！
    # input parameter:
    # astart, afinal: 计算的裂纹长度范围的起点和终点，起点推荐为高载时裂纹长度，终点建议为有效的最大值
    # b, w, ys：试件基本参数thickness，width和yield strength
    # pmax，r：加载基本参数循环载荷峰值pmax和应力比r
    # aol，rol：高载时裂纹长度aol及其对应的塑性区尺寸rol
    # plasticzonefactor：塑性区尺寸计算时的塑性区尺寸系数
    # m1：由WheelerFitting函数拟合得到的Wheeler迟滞参数的指数系数，表征迟滞恢复的快慢
    # c，m：基础Paris公式参数
    # return parameter:
    # dadn_wheeler, dk_wheeler, a_wheeler，cp_wheeler模型计算结果对应的dadn，SIF变幅dk和裂纹长度a，及迟滞系数cp
    a_wheeler = np.arange(astart, afinal, 0.02)
    dk_wheeler = mts_analysis.DeltaKCalculating(b=b, w=w, a=a_wheeler, pmax=pmax, pmin=pmax * r)
    kmax_wheeler = np.array([dk / (1 - r) for dk in dk_wheeler])
    rm_wheeler = PlasticZoneWithFactor(kmax=kmax_wheeler, ys=ys, factor=plasticzonefactor)
    left_wheeler = a_wheeler + rm_wheeler
    right_wheeler = aol + rol
    # 迟滞计算参数和数据准备，由裂纹长度a出发

    cp = []
    for seq, _ in enumerate(a_wheeler):
        if left_wheeler[seq] < right_wheeler:
            cpi = (rm_wheeler[seq] / (aol + rol - a_wheeler[seq])) ** m1
        else:
            cpi = 1
        cp.append(cpi)
    cp = np.array(cp)
    # 迟滞系数cp计算

    dadn_wheeler = []
    for seq, _ in enumerate(dk_wheeler):
        dadn = cp[seq] * (c * dk_wheeler[seq] ** m)
        dadn_wheeler.append(dadn)
    dadn_wheeler = np.array(dadn_wheeler)
    # Wheeler模型 拟合函数计算
    return dadn_wheeler, dk_wheeler, a_wheeler, cp


def CpCalculatingBaseParis(a, b, w, ys, pmax, r, aol, rol, plasticzonefactor, m1):
    # usage: 计算裂纹长度a对应的Wheeler迟滞系数cp
    # NOTE：数据的起始点必须是高载后！
    # input parameter:
    # a: 裂纹长度
    # b, w, ys：试件基本参数thickness，width和yield strength
    # pmax，r：加载基本参数循环载荷峰值pmax和应力比r
    # aol，rol：高载时裂纹长度aol及其对应的塑性区尺寸rol
    # plasticzonefactor：塑性区尺寸计算时的塑性区尺寸系数
    # m1：由WheelerFitting函数拟合得到的Wheeler迟滞参数的指数系数，表征迟滞恢复的快慢
    # return parameter:
    # dadn_wheeler, dk_wheeler, a_wheeler，cp_wheeler模型计算结果对应的dadn，SIF变幅dk和裂纹长度a，及迟滞系数cp
    a_wheeler = a
    dk_wheeler = mts_analysis.DeltaKCalculating(b=b, w=w, a=a_wheeler, pmax=pmax, pmin=pmax * r)
    kmax_wheeler = np.array([dk / (1 - r) for dk in dk_wheeler])
    rm_wheeler = PlasticZoneWithFactor(kmax=kmax_wheeler, ys=ys, factor=plasticzonefactor)
    left_wheeler = a_wheeler + rm_wheeler
    right_wheeler = aol + rol
    # 迟滞计算参数和数据准备，由裂纹长度a出发

    cp = []
    for seq, _ in enumerate(a_wheeler):
        if left_wheeler[seq] < right_wheeler:
            cpi = (rm_wheeler[seq] / (aol + rol - a_wheeler[seq])) ** m1
        else:
            cpi = 1
        cp.append(cpi)
    cp = np.array(cp)
    # 迟滞系数cp计算

    return cp


def SalvatiWheelerCalculatingBaseParis(astart, afinal, b, w, ys, pmax, r, aol, rol, plasticzonefactor, m1, c, m,
                                       apre, rdol, delayfactor, m_mod):
    # usage: 以裂纹长度a指定计算范围，计算基于Paris公式的Wheeler模型
    # NOTE：数据的起始点必须是高载后！
    # input parameter:
    # astart, afinal: 计算的裂纹长度范围的起点和终点，起点推荐为高载时裂纹长度，终点建议为有效的最大值
    # b, w, ys：试件基本参数thickness，width和yield strength
    # pmax，r：加载基本参数循环载荷峰值pmax和应力比r
    # aol，rol：高载时裂纹长度aol及其对应的塑性区尺寸rol
    # plasticzonefactor：塑性区尺寸计算时的塑性区尺寸系数
    # m1：由WheelerFitting函数拟合得到的Wheeler迟滞参数的指数系数，表征迟滞恢复的快慢
    # c，m：基础Paris公式参数
    # apre：高载起点（dadN最高值）
    # rdol：采用延迟系数计算的延迟区rd,OL
    # delayfactor：Yuen等提出的延迟系数
    # m_mod：由Salvati延迟系数拟合得到的延迟指数
    # return parameter:
    # dadn_wheeler, dk_wheeler, a_wheeler，cp, cd:模型计算结果对应的dadn，SIF变幅dk和裂纹长度a，及迟滞和延迟系数cp,cd
    a_wheeler = np.arange(astart, afinal, 0.02)
    dk_wheeler = mts_analysis.DeltaKCalculating(b=b, w=w, a=a_wheeler, pmax=pmax, pmin=pmax * r)
    kmax_wheeler = np.array([dk / (1 - r) for dk in dk_wheeler])
    rm_wheeler = PlasticZoneWithFactor(kmax=kmax_wheeler, ys=ys, factor=plasticzonefactor)
    left_wheeler = a_wheeler + rm_wheeler
    right_wheeler = aol + rol
    rd_wheeler = PlasticZoneWithFactor(kmax=kmax_wheeler, ys=ys, factor=delayfactor)
    left_delay = a_wheeler + rd_wheeler
    right_delay = apre + rdol
    # 迟滞计算参数和数据准备，由裂纹长度a出发

    cp = []
    for seq, _ in enumerate(a_wheeler):
        if left_wheeler[seq] < right_wheeler:
            cpi = (rm_wheeler[seq] / (aol + rol - a_wheeler[seq])) ** m1
        else:
            cpi = 1
        cp.append(cpi)
    cp = np.array(cp)
    # 迟滞系数cp计算

    cd = []
    for seq, _ in enumerate(a_wheeler):
        if left_delay[seq] < right_delay:
            cdi = ((apre + rdol - a_wheeler[seq])/rd_wheeler[seq]) ** m_mod
        else:
            cdi = 1
        cd.append(cdi)
    cd = np.array(cd)

    dadn_salvati = []
    for seq, _ in enumerate(dk_wheeler):
        dadn = cd[seq] * cp[seq] * (c * dk_wheeler[seq] ** m)
        dadn_salvati.append(dadn)
    dadn_salvati = np.array(dadn_salvati)
    # Wheeler模型 拟合函数计算
    return dadn_salvati, dk_wheeler, a_wheeler, cd, cp