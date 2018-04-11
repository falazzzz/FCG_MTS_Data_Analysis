import math
import matplotlib.pyplot as plt
import numpy as np
from FCGAnalysisLib import read_data
from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import overload_analysis


sequence = ["yang-baoban_Lu-420-09"]
save = 1
stress_ratio = 0.7
threshold = 0
c_ca = 5.6665e-09
m_ca = 3.0965
aol = 13
pol = 8000
ppeak = 4000
# 参数设置

specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence[0])
cycles, cracklength, kmax, kmin, pmax, pmin = read_data.ReadOriginResult(sequence[0], closure=False, cod=False)
dadn_Manual, n_Manual, dk_Manual, a_Manual = \
    experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin,
                                                   a=cracklength,
                                                   ys=yield_strength, r=stress_ratio, threshold=threshold)
# 实验数据读取和初步处理

kol = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=aol, pmax=pol, pmin=0)
rol = overload_analysis.PlasticZoneWithFactor(kmax=np.array([kol]), ys=yield_strength, factor=1/math.pi)
# 高载参数计算

dadn_ola, dk_ola, n_ola, a_ola = \
    mts_analysis.FCGDataSelectByThreshold(dadn=dadn_Manual, dk=dk_Manual, n=n_Manual, a=a_Manual,
                                          threshold=aol, target='a', keepbigger=1)
# 以高载时裂纹长度aol筛选出高载后的数据

kmax_ola = np.array([dk/(1-stress_ratio) for dk in dk_ola])
rm_ola = overload_analysis.PlasticZoneWithFactor(kmax=kmax_ola, ys=yield_strength, factor=1/math.pi)
# 高载段参数计算

m1, _ = overload_analysis.WheelerFittingBaseParis(a=a_ola, dadn=dadn_ola, dk=dk_ola, rm=rm_ola, aol=aol, rol=rol, c=c_ca, m=m_ca)
# Wheeler模型指数系数m1拟合

dadn_wheeler, dk_wheeler, a_wheeler = \
    overload_analysis.WheelerCalculatingBaseParis(astart=aol, afinal=max(a_Manual),
                                                  b=thickness, w=width, ys=yield_strength, pmax=ppeak, r=stress_ratio,
                                                  aol=aol, rol=rol, plasticzonefactor=1/math.pi,
                                                  m1=m1, c=c_ca, m=m_ca)
# Wheeler模型计算

# Plotting
plt.figure(num=1, figsize=(10, 8))
plt.scatter(dk_Manual, dadn_Manual, lw=1, marker='+', label='Experiment')
plt.plot(dk_wheeler, dadn_wheeler, label='Wheeler Model Fitting', color='black', linewidth=2)
plt.title("FCG Rates - deltaK(Wheeler Model),OLR="+str(pol/ppeak))
plt.ylabel("FCG Rates/mm per cycle")
plt.xlabel("DeltaSIF/MPa.m0.5")
plt.xscale('log')
plt.yscale('log')
plt.axis([min(10, min(dk_Manual)), max(20, max(dk_Manual)), min(dadn_Manual)*0.9, max(dadn_Manual)*1.1])
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()
if save:
    plt.savefig('WheelerModel_' + sequence[0] + '.png', dpi=320)
plt.show()
