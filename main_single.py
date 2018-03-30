# MTS单组实验数据处理函数
# 适用采用MTS内置Fatigue Crack Growth程序实验得到的数据格式
# 读取的文件包括原始数据（CrackLength等），MTS处理得到FCGR-DK数据，以及Test Summary
# 主要内容：
# (1)读取文件中的相应数据
# (2)将MTS的数据进行Paris拟合，得到参数c_MTS, m_MTS
# (3)读取原始数据中的裂纹长度，依据E647-11计算SIF变幅，采用割线法计算FCGR
# (4)根据E647-11中对韧带长度有效性要求筛去不符合的数据，根据指定的阈值（目前为无）筛去不符合的数据，得到满足要求的FCGR和DK
# (5)将由原始数据手动计算的数据进行Paris拟合，得到参数c_manual，m_Manual
# (6)根据指定的裂纹长度和SIF变幅，分别依据原始数据和MTS数据推算相应的循环数
# (7)绘图
# 全局单位设定：
# length: mm
# dadn: mm/cycle
# SIF: MPa.m0.5
# stress, modulus: GPa
# Load：N


import read_data
import mts_analysis
import experiment_calculation
import matplotlib.pyplot as plt
import numpy as np


# 分析和绘图参数
show = 1
save = 0
way0 = 1                                   # 采用da dN, Delta K by Cycles.csv文件提供的结果计算
way1 = 1                                   # 采用Crack Length,Min Max K, Load by Cycles.csv文件提供的结果计算
way2 = 1                                   # 采用Channels by Cycles.csv文件提供的结果进行计算
sequence = "yang-baoban_Lu-420-05"         # Graph Saving Sequence
cracklength_for_overload = 12              # For Finding Part, Find the cycles when length = 10mm
dk_for_overload = 17                       # For Finding Part, Find the cycles when SIF = 30 MPa.m0.5
stress_ratio = 0.1                         # Stress Ratio
threshold = 0                           # Threshold for paris region selecting


# 实验基本参数读取
specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence)
print('specimen name:', str(specimen))
# mm for length, GPa for strength and modulus


# MTS数据读取和处理（Way0）
if way0:
    dadn_MTS, cycles_MTS, dk_MTS = read_data.ReadMtsResult(sequence=sequence)
    c_MTS, m_MTS = mts_analysis.ParisFitting(dadn=dadn_MTS, dk=dk_MTS)
    c_MTS = 7.9713e-10  # Temp
    m_MTS = 3.6797
    dadn_paris_by_MTS = mts_analysis.ParisCalculating(c=c_MTS, m=m_MTS, dk=dk_MTS)
    print("MTS Results Fitting Result by Paris:c=", str(c_MTS), ",m=", str(m_MTS))


# 由裂纹长度、载荷、循环次数计算（Way1）
if way1:
    cycles, cracklength, kmax, kmin, pmax, pmin, codmax, codmin, = \
        read_data.ReadOriginResult(sequence=sequence)

    dadn_Manual, n_Manual, dk_Manual, a_Manual = \
        experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin,
                                                       a=cracklength,
                                                       ys=yield_strength, r=stress_ratio, threshold=0)

    c_Manual, m_Manual = mts_analysis.ParisFitting(dadn=dadn_Manual, dk=dk_Manual)
    dadn_paris_by_Manual = mts_analysis.ParisCalculating(c=c_Manual, m=m_Manual, dk=dk_Manual)
    print("Manual Results Fitting by Paris: c=", str(c_Manual), ",m=", str(m_Manual))


# 由循环数、轴力值、COD值详表计算（Way2）：
if way2:
    cycles_r, load_r, cod_r = read_data.ReadCodData(sequence)

    cycles_way2, cod_way2, cod_way2d, loadmax_way2, loadmin_way2 = \
        mts_analysis.OutputCodandLoad(cycles_r, load_r, cod_r, rate1=1, rate2=0)

    cracklength_r = mts_analysis.Compliance(e=elastic_modulus, b=thickness, w=width, p=(loadmax_way2 - loadmin_way2),
                                            v=cod_way2 - cod_way2d)

    dadn_way2, n_way2, dk_way2, a_way2 = \
        experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, a=cracklength_r,
                                                       n=cycles_way2, pmax=loadmax_way2, pmin=loadmin_way2,
                                                       ys=yield_strength, r=stress_ratio, threshold=threshold)

    c_way2, m_way2 = mts_analysis.ParisFitting(dadn=dadn_way2, dk=dk_way2)
    dadn_paris_by_way2 = mts_analysis.ParisCalculating(c=c_way2, m=m_way2, dk=dk_way2)
    print("Way2 Results Fitting by Paris: c=", str(c_way2), ",m=", str(m_way2))


#  Finding Information for Overload
'''
cycle_for_overload1 = int(mts_analysis.FindAscentDataBySeq(value=cracklength_for_overload, item=cracklength, target=cycles))
print("Predicting Cycle when Crack Length =", str(cracklength_for_overload), "mm:", str(cycle_for_overload1))
# Finding Cycles by CrackLength(by Origin Data)
cycle_for_overload2 = int(mts_analysis.FindAscentDataBySeq(value=dk_for_overload, item=dk_MTS, target=cycles_MTS))
print("Predicting Cycle when DeltaK =", str(dk_for_overload), "Mpa.m0.5:", str(cycle_for_overload2))
'''
# Finding Cycles by CrackLength(by MTS Data)


# Plotting
# Plot: (1)Crack Length - Cycles (Way1, 2)
plt.figure(num=1, figsize=(10, 8))
if way1:
    plt.scatter(n_Manual, a_Manual, s=1, label='$Way1$', color='green', lw=1)
if way2:
    plt.scatter(n_way2, a_way2, s=1, label='$Way2$', color='blue', lw=1)
plt.xlabel("Cycles")
plt.ylabel("Crack Length (mm)")
plt.title('Crack Length - Cycles ' + sequence)
plt.legend()
plt.grid()
if save:
    plt.savefig('a-n_' + sequence + '_main_single.png', dpi=320)
if show:
    plt.show()

# Plot: (2)da/dN - dk plot (Only Way0)
if way0:
    plt.figure(num=2, figsize=(10, 8))
    plt.scatter(dk_MTS, dadn_MTS, s=1, label='$ExperimentData$', color='red', lw=1)
    plt.plot(dk_MTS, dadn_paris_by_MTS, label='$Fitting By Paris$', color='blue', linewidth=2)
    plt.axis([min(dk_MTS) * 0.95, max(dk_MTS) * 1.05, min(dadn_MTS) * 0.95, max(dadn_MTS) * 1.05])
    textplacex = (max(dk_MTS) - min(dk_MTS)) * 0.6 + min(dk_MTS)
    plt.text(textplacex, min(dadn_MTS), 'c=%.4e' % c_MTS + ',m=%.4f' % m_MTS)
    plt.xlabel("DeltaK Applied (MPa*m^0.5)")
    plt.xscale('log')
    plt.ylabel("da/dN (mm/cycle)")
    plt.yscale('log')
    plt.yticks(np.linspace(min(dadn_MTS), max(dadn_MTS), 6))
    plt.xticks(np.linspace(min(dk_MTS), max(dk_MTS), 6))
    plt.title('da/dN - dK ' + sequence + '(MTS Result)')
    plt.legend()
    plt.grid()
    if save:
        plt.savefig('dadn-dk_MTS_' + sequence + '_main_single.png', dpi=320)
    if show:
        plt.show()



# Plot: da/dN - dk plot (Way0, 1, 2)
plt.figure(num=3, figsize=(10, 8))
if way0:
    plt.scatter(dk_MTS, dadn_MTS, s=1, label='$MTS Results$', color='red', lw=1)
    plt.plot(dk_MTS, dadn_paris_by_MTS, label='$MTS Results Fitting$', color='red', linewidth=2)
if way1:
    plt.scatter(dk_Manual, dadn_Manual, s=1, label='$way1 Results$', color='green', lw=1)
    plt.plot(dk_Manual, dadn_paris_by_Manual, label='$Way1 Results Fitting$', color='green', linewidth=2)
if way2:
    plt.scatter(dk_way2, dadn_way2, s=1, label='$way2 Results$', color='blue', lw=1)
    plt.plot(dk_way2, dadn_paris_by_way2, label='$way2 Results Fitting$', color='blue', linewidth=2)
plt.axis([min(dk_MTS)*0.95, max(dk_MTS)*1.05, min(dadn_MTS)*0.95, max(dadn_MTS)*1.05])
plt.text(textplacex, min(dadn_MTS), 'MTS:c=%.4e' % c_MTS+',m=%.4f' % m_MTS+'\nWay2:c=%.4e' % c_way2+',m=%.4f' % m_way2)
plt.xlabel("DeltaK Applied (MPa*m^0.5)")
plt.xscale('log')
plt.ylabel("da/dN (mm/cycle)")
plt.yscale('log')
plt.yticks(np.linspace(min(dadn_MTS), max(dadn_MTS), 6))
plt.xticks(np.linspace(min(dk_MTS), max(dk_MTS), 6))
plt.title('da/dN - dK '+sequence+'(MTS compare with Manual)')
plt.legend()
plt.grid()
if save:
    plt.savefig('dadn-dk_MTS_and_Manual_'+sequence+'_main_single.png', dpi=320)
if show:
    plt.show()
