# LU YUNCHAO
# 2018/10/15 V1.0.0
# Li等提出的基于Osgood模型的Wheeler模型优化的拟合模型实现程序，效果一般

from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import read_data
from FCGAnalysisLib import tensiontest_analysis
from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import overload_analysis
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import optimize

ys = 0.446
E = 190
sequence = ["b-1#"]
extension, load, stress, strain = read_data.ReadTensionResult(sequence=sequence[0])
n, alpha = experiment_calculation.RambergOsgoodFitFromOriginal(s=stress, e=strain, ys=ys, E=E, uplimit=0.25, lowlimit=0.01)
truestress = tensiontest_analysis.TrueStress(engineeringstress=stress, engineeringstrain=strain)
truestrain = tensiontest_analysis.TrueStrain(engineeringstrain=strain)
stress_osgood = np.arange(0, 700, 1)
strain_osgood = tensiontest_analysis.RambergOsgoodCalculation(n=n, alpha=alpha, s=stress_osgood, ys=ys, E=E)

plt.figure(num=1, figsize=(10, 8))
plt.scatter(strain, stress, label='$Engineering$', marker='1')
plt.scatter(truestrain, truestress, label='$True$', marker='2')
plt.plot(strain_osgood, stress_osgood, label='$OsgoodFit$', color='blue', linewidth=2)
plt.title("Stress-Strain Curve")
plt.ylabel("Stress/MPa")
plt.xlabel("Strain")
plt.legend()
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')

n = 6.7278
sequence = ["yang-baoban_Lu-420-06"]
save = 0
stress_ratio = 0.1
threshold = 0
if stress_ratio == 0.1:
    c_ca = 7.9713e-10
    m_ca = 3.6797
elif stress_ratio == 0.7:
    c_ca = 5.6665e-09
    m_ca = 3.0965
aol = 12.3495   # 高载时dadN最低值对应的裂纹长度
apre = 12.10  # 高载时dadN最高值对应的裂纹长度
ad = aol - apre  # retardation delay长度
pol = 4000
ppeak = 2000
factor = (n - 1)/(n + 1) * 1/math.pi
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
rol = overload_analysis.PlasticZoneWithFactor(kmax=np.array([kol]), ys=yield_strength, factor=factor)
# 高载参数计算

dadn_ola, dk_ola, n_ola, a_ola = \
    mts_analysis.FCGDataSelectByThreshold(dadn=dadn_Manual, dk=dk_Manual, n=n_Manual, a=a_Manual,
                                          threshold=aol, target='a', keepbigger=1)
# 以高载时裂纹长度aol筛选出高载后的数据

kmax_ola = np.array([dk/(1-stress_ratio) for dk in dk_ola])
rm_ola = overload_analysis.PlasticZoneWithFactor(kmax=kmax_ola, ys=yield_strength, factor=factor)
# 高载段参数计算

m1, _ = overload_analysis.WheelerFittingBaseParis(a=a_ola, dadn=dadn_ola, dk=dk_ola, rm=rm_ola, aol=aol, rol=rol, c=c_ca, m=m_ca)
# Wheeler模型指数系数m1拟合

dadn_wheeler, dk_wheeler, a_wheeler, cp_wheeler = \
    overload_analysis.WheelerCalculatingBaseParis(astart=apre, afinal=max(a_Manual),
                                                  b=thickness, w=width, ys=yield_strength, pmax=ppeak, r=stress_ratio,
                                                  aol=aol, rol=rol, plasticzonefactor=factor,
                                                  m1=m1, c=c_ca, m=m_ca)
# Wheeler模型计算

# Plotting
plt.figure(num=2, figsize=(10, 8))
plt.scatter(dk_Manual, dadn_Manual, lw=1, marker='+', label='Experiment')
plt.plot(dk_wheeler, dadn_wheeler, label='Wheeler Fitting', color='black', linewidth=2)
plt.title("FCG Rates - deltaK(Wheeler Model) with Li's factor,OLR="+str(pol/ppeak))
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
