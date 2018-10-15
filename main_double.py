from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import read_data
from FCGAnalysisLib import mts_analysis
import matplotlib.pyplot as plt
import numpy as np

# main1: for origin data dispose, drawing figures
# 用于CA数据处理.
# Lu Yunchao

# global unit
# length: mm
# dadn: mm/cycle
# SIF: MPa.m0.5
# stress, modulus: GPa
# Load：N

#Switch
show = 1
save = 0
sequence = ["yang-baoban_Lu-420-11", "yang-baoban_Lu-420-13"]         # Graph Saving Sequence
stress_ratio = 0.5     # Stress Ratio
threshold = [0, 0]
print('Specimens:', sequence)
# Specimen1
dadn_MTS1, cycles_MTS1, dk_MTS1 = read_data.ReadMtsResult(sequence=sequence[0])
cycles1, cracklength1, kmax1, kmin1, pmax1, pmin1, codmax1, codmin1, = \
    read_data.ReadOriginResult(sequence=sequence[0])
specimen1, width1, notch_length1, thickness1, elastic_modulus1, yield_strength1, precrack_length1 = \
    read_data.ReadTestInf(sequence=sequence[0])
# 读取三个文件
dadn_Manual1, n_Manual1, dk_Manual1, a_Manual1 = \
    experiment_calculation.FCGRandDKbyOriginalData(b=thickness1, w=width1, a=cracklength1, n=cycles1,
                                                   pmax=pmax1, pmin=pmin1, ys=yield_strength1, r=stress_ratio,
                                                   threshold=threshold[0])
# 输出FCGR、循环次数、SIF变幅、裂纹长度（Manual）

# Specimen2
dadn_MTS2, cycles_MTS2, dk_MTS2 = read_data.ReadMtsResult(sequence=sequence[1])
cycles2, cracklength2, kmax2, kmin2, pmax2, pmin2, codmax2, codmin2, = \
    read_data.ReadOriginResult(sequence=sequence[1])
specimen2, width2, notch_length2, thickness2, elastic_modulus2, yield_strength2, precrack_length2 = \
    read_data.ReadTestInf(sequence=sequence[1])
# 读取三个文件
dadn_Manual2, n_Manual2, dk_Manual2, a_Manual2 = \
    experiment_calculation.FCGRandDKbyOriginalData(b=thickness2, w=width2, a=cracklength2, n=cycles2,
                                                   pmax=pmax2, pmin=pmin2, ys=yield_strength2, r=stress_ratio,
                                                   threshold=threshold[1])

# 合并两组Manual数据进行Paris拟合
dadn_Manual = np.concatenate((dadn_Manual1, dadn_Manual2))
dk_Manual = np.concatenate((dk_Manual1, dk_Manual2))
c_Manual, m_Manual = mts_analysis.ParisFitting(dadn=dadn_Manual, dk=dk_Manual)
print("Manual Results Fitting Result by Paris:c=", str(c_Manual), ",m=", str(m_Manual))
# Paris拟合得到c，m
dk_paris_sorted_manual = sorted(dk_Manual)
dadn_paris_manual = mts_analysis.ParisCalculating(c=c_Manual, m=m_Manual, dk=dk_paris_sorted_manual)
# 得到Paris曲线

# 合并两组MTS数据进行Paris拟合
dadn_MTS = np.concatenate((dadn_MTS1, dadn_MTS2))
dk_MTS = np.concatenate((dk_MTS1, dk_MTS2))
c_MTS, m_MTS = mts_analysis.ParisFitting(dadn=dadn_MTS, dk=dk_MTS)
print("MTS Results Fitting Result by Paris:c=", str(c_MTS), ",m=", str(m_MTS))
# Paris拟合得到c，m
dk_paris_sorted_MTS = sorted(dk_MTS)
dadn_paris_MTS = mts_analysis.ParisCalculating(c=c_MTS, m=m_MTS, dk=dk_paris_sorted_MTS)
# 得到Paris曲线

# 绘图部分
name = 'R=' + str(stress_ratio) + '_QSTE420TM'
textplacex = (max(dk_paris_sorted_manual) - min(dk_paris_sorted_manual)) * 0.5 + min(dk_paris_sorted_manual)
# Plotting: (1)da/dN - dk(Manual) plot(2次实验结果对比)
plt.figure(num=1, figsize=(7, 5))
plt.scatter(dk_Manual1, dadn_Manual1, s=1, label='$Pmax=3.6kN$', color='red', lw=1)
plt.scatter(dk_Manual2, dadn_Manual2, s=1, label='$Pmax=2.8kN$', color='blue', lw=1)
plt.plot(dk_paris_sorted_manual, dadn_paris_manual, label='$Fitting By Paris$', color='black', linewidth=2)
plt.axis([min(dk_paris_sorted_manual)*0.95, max(dk_paris_sorted_manual)*1.05,
          min(dadn_paris_manual)*0.95, max(dadn_paris_manual)*1.05])
#plt.text(textplacex, min(dadn_Manual), 'Manual:c=%.4e' % c_Manual+',m=%.4f' % m_Manual)
plt.xlabel("DeltaK Applied (MPa*m^0.5)")
plt.xscale('log')
plt.xticks(np.linspace(min(dk_paris_sorted_manual), max(dk_paris_sorted_manual), 5))
plt.ylabel("da/dN (mm/cycle)")
plt.yscale('log')
plt.title('da/dN - dK ' + name + '(Manual Result)')
plt.legend()
plt.grid()
if save:
    plt.savefig('CA_dadn_dk_R07_Manual_main_double.png', dpi=320)
if show:
    plt.show()

# Plotting: (2)da/dN - dk(MTS) plot(2次实验结果对比)
textplacex = (max(dk_paris_sorted_MTS) - min(dk_paris_sorted_MTS)) * 0.5 + min(dk_paris_sorted_MTS)
plt.figure(num=2, figsize=(7, 5))
plt.scatter(dk_MTS1, dadn_MTS1, s=1, label='$Pmax=3.6kN$', color='red', lw=1)
plt.scatter(dk_MTS2, dadn_MTS2, s=1, label='$Pmax=2.8kN$', color='blue', lw=1)
plt.plot(dk_paris_sorted_MTS, dadn_paris_MTS, label='$Fitting By Paris$', color='black', linewidth=2)
plt.axis([min(dk_paris_sorted_MTS)*0.95, max(dk_paris_sorted_MTS)*1.05,
          min(dadn_paris_MTS)*0.95, max(dadn_paris_MTS)*1.05])
#plt.text(textplacex, min(dadn_MTS), 'MTS:c=%.4e' % c_MTS+',m=%.4f' % m_MTS)
plt.xlabel("DeltaK Applied (MPa*m^0.5)")
plt.xscale('log')
plt.ylabel("da/dN (mm/cycle)")
plt.yscale('log')
#plt.xticks(np.linspace(min(dk_paris_sorted_MTS), max(dk_paris_sorted_MTS), 6))
plt.title('da/dN - dK ' + name + '(MTS Result)')
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()
plt.grid()
if save:
    plt.savefig('CA_dadn_dk_R07_MTS_main_double.png', dpi=320)
if show:
    plt.show()
