import math
import matplotlib.pyplot as plt
import numpy as np
from FCGAnalysisLib import read_data
from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import overload_analysis


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
apre = 11.9671  # 高载时dadN最高值对应的裂纹长度
ad = 12.4516 - apre  # retardation delay长度
pol = 4000
ppeak = 2000
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

dadn_wheeler, dk_wheeler, a_wheeler, cp_wheeler = \
    overload_analysis.WheelerCalculatingBaseParis(astart=apre, afinal=max(a_Manual),
                                                  b=thickness, w=width, ys=yield_strength, pmax=ppeak, r=stress_ratio,
                                                  aol=aol, rol=rol, plasticzonefactor=1/math.pi,
                                                  m1=m1, c=c_ca, m=m_ca)
# Wheeler模型计算
'''
# 采用Mehrzadi, 2013提出的延迟系数计算
kd = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=aol + ad, pmax=ppeak, pmin=ppeak*stress_ratio)
alpha_d = ad*1e-3 / ((kol/yield_strength*1e-3)**2 - (kd/yield_strength*1e-3)**2)
rd_ol = overload_analysis.PlasticZoneWithFactor(kmax=[kol], ys=yield_strength, factor=alpha_d)
rd_i =overload_analysis.PlasticZoneWithFactor(kmax=dk_ola/(1-stress_ratio), ys=yield_strength, factor=alpha_d)
cpmin = np.min(cp_wheeler)
left = a_ola + rd_i
right = aol + rd_ol
if left[-1] < right:
    dadn_old = dadn_ola
    dk_old = dk_ola
    a_old = a_ola
    rm_old = rm_ola
    print("Delay Region Check: Fitting Region:dk_start=" + str(dk_old[0]) + ",to the END.")
elif left[0] >= right:
    print("No Delay Region, Check Input Parameter!")
else:
    for seq, _ in enumerate(left):
        if left[seq] >= right:
            dadn_old = dadn_ola[0:seq]
            dk_old = dk_ola[0:seq]
            a_old = a_ola[0:seq]
            rm_old = rm_ola[0:seq]
            print("Delay Region Check: Fitting Region:dk_start=" + str(dk_old[0]) + ",dk_final=" + str(dk_old[-1]))
            break
cd_wheeler = np.exp((a_old - aol)/ad * np.log(cpmin))
cd_wheeler = np.concatenate((cd_wheeler, np.full(len(cp_wheeler)-len(cd_wheeler), 1)))
# 延迟系数D的计算

c_wheeler = cd_wheeler*cp_wheeler
dadn_wheelermodified = c_wheeler * (c_ca * dk_wheeler**m_ca)
# 用带延迟系数D的模型计算dadn
'''
# 采用Salvati, 2016提出的延迟系数计算
dadn_old, dk_old, n_old, a_old = \
    mts_analysis.FCGDataSelectByThreshold(dadn=dadn_Manual, dk=dk_Manual, n=n_Manual, a=a_Manual,
                                          threshold=apre, target='a', keepbigger=1)
kd = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=apre+ad, pmax=ppeak, pmin=ppeak*stress_ratio)
alpha_d = ad*1e-3 / ((kol/yield_strength*1e-3)**2 - (kd/yield_strength*1e-3)**2)
rd_ol = overload_analysis.PlasticZoneWithFactor(kmax=[kol], ys=yield_strength, factor=alpha_d)
rd_i = overload_analysis.PlasticZoneWithFactor(kmax=dk_old/(1-stress_ratio), ys=yield_strength, factor=alpha_d)
cp_old = overload_analysis.CpCalculatingBaseParis(a=a_old, b=thickness, w=width, ys=yield_strength, pmax=ppeak,
                                                  r=stress_ratio, aol=aol, rol=rol, plasticzonefactor=1/math.pi, m1=m1)
m_mod, _ = overload_analysis.DelayFittingBaseParis(a=a_old, dadn=dadn_old, dk=dk_old, cp=cp_old,
                                                   rd=rd_i, aol=apre, rdol=rd_ol, c=c_ca, m=m_ca)

dadn_salvati, dk_salvati, a_salvati, cd_salvati, cp_salvati = \
    overload_analysis.SalvatiWheelerCalculatingBaseParis(astart=apre, afinal=max(a_Manual),
                                                         b=thickness, w=width, ys=yield_strength,
                                                         pmax=ppeak, r=stress_ratio, aol=aol, rol=rol,
                                                         plasticzonefactor=1/math.pi, m1=m1, c=c_ca, m=m_ca,
                                                         apre=apre, rdol=rd_ol, delayfactor=alpha_d, m_mod=m_mod)
# Salvati修正的Wheeler模型拟合

# Plotting
plt.figure(num=1, figsize=(10, 8))
plt.scatter(dk_Manual, dadn_Manual, lw=1, marker='+', label='Experiment')
plt.plot(dk_wheeler, dadn_wheeler, label='Wheeler Fitting', color='black', linewidth=2)
plt.plot(dk_salvati, dadn_salvati, label='Modified Wheeler Fitting', linewidth=2)
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
