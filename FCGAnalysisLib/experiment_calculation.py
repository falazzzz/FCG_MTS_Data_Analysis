# 实验结果原始值计算函数
from FCGAnalysisLib import mts_analysis
import numpy as np


def FCGRandDKbyOriginalData(b, w, a, n, pmax, pmin, ys, r, kclosure=[], threshold=0):
    # usage: 依据标准E647-11由载荷峰谷值、裂纹长度、循环次数计算FCGR并进行有效性检查
    # input parameter:
    # b：试件厚度thickness，mm
    # w：试件宽度Width，mm
    # a：裂纹长度数组，mm
    # n：循环次数数组
    # pmax：载荷峰值数组，N
    # pmin：载荷谷值数组，N
    # ys：材料屈服强度，GPa
    # r：应力比
    # kclosure：裂纹闭合时应力强度因子，默认为空数组（不输出）
    # threshold：DeltaK阈值筛选，将删去小于阈值的数据，默认为0（不筛选）
    # return parameter：
    # dadn_ligament_valid：裂纹扩展速率
    # n_ligament_valid：循环次数
    # dk_ligament_valid：SIF变幅
    # a_ligament_valid：裂纹长度
    dk = mts_analysis.DeltaKCalculating(b=b, w=w, a=a, pmax=pmax, pmin=pmin)
    # 由裂纹长度依标准计算SIF变幅
    dadn_poly, a_poly, n_poly, dk_poly = mts_analysis.FCGRateByPolynomial(a=a, n=n, dk=dk, num=4)
    # 由裂纹长度和循环次数依多项式法计算FCGR
    n_ligament_valid = mts_analysis.DataSelectByLigament(w=w, a=a_poly, dk=dk_poly, ys=ys, data=n_poly, r=r)
    a_ligament_valid = mts_analysis.DataSelectByLigament(w=w, a=a_poly, dk=dk_poly, ys=ys, data=a_poly, r=r)
    dadn_ligament_valid = mts_analysis.DataSelectByLigament(w=w, a=a_poly, dk=dk_poly, ys=ys, data=dadn_poly, r=r)
    dk_ligament_valid = mts_analysis.DataSelectByLigament(w=w, a=a_poly, dk=dk_poly, ys=ys, data=dk_poly, r=r)
    if len(kclosure) != 0:
        kc_ligament_valid = mts_analysis.DataSelectByLigament(w=w, a=a_poly, dk=dk_poly, ys=ys, data=kclosure, r=r)
    print("Data Deleted by Ligament Check:", str(len(a_poly) - len(a_ligament_valid)))
    # 依据标准进行韧带长度有效性筛选
    n_valid = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk_ligament_valid, data=n_ligament_valid)
    a_valid = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk_ligament_valid, data=a_ligament_valid)
    dadn_valid = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk_ligament_valid, data=dadn_ligament_valid)
    dk_valid = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk_ligament_valid, data=dk_ligament_valid)
    if len(kclosure) != 0:
        kc_valid = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk_ligament_valid, data=kc_ligament_valid)
    # 依据DeltaK阈值进行筛选（可用于筛去DeltaK小于Paris区域的值）
    if len(kclosure) == 0:
        return np.array(dadn_valid), np.array(n_valid), np.array(dk_valid), np.array(a_valid)
    else:
        return np.array(dadn_valid), np.array(n_valid), np.array(dk_valid), np.array(a_valid), np.array(kc_valid)

