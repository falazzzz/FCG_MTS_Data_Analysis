# 采用Salvati提出的带延迟的Wheeler模型计算
# 2018/7/2 Version 2.0
# Lu Yunchao

import math
import matplotlib.pyplot as plt
import numpy as np
from FCGAnalysisLib import read_data
from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import overload_analysis
from FCGAnalysisLib import paris_and_walker


# 实验参数
sequence = ["yang-baoban_Lu-420-06"]
stress_ratio = 0.1
threshold = 0
pol = 4000
ppeak = 2000

# 程序参数
figsave = 0         # 图像保存开关
show = 1            # 图像显示开关

# paris公式参数读取
c_ca, m_ca = paris_and_walker.ParisParameter(r=stress_ratio)

# Wheeler模型拟合设定参数
a_ol = 12.0057
a_ol_applied = 12.0057
a_fcgr_min = 12.4516   # 高载时dadN最低值对应的裂纹长度
a_fcgr_max = 12.0057  # 高载时dadN最高值对应的裂纹长度

ad = a_fcgr_min - a_fcgr_max  # retardation delay长度

# 实验数据读取和初步处理(本函数不能用BasicDataDispose封装API，因函数要调用细节数据）
specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence[0])
cycles, cracklength, kmax, kmin, pmax, pmin = read_data.ReadOriginResult(sequence[0], closure=False, cod=False)
dadn_Manual, n_Manual, dk_Manual, a_Manual = \
    experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin,
                                                   a=cracklength,
                                                   ys=yield_strength, r=stress_ratio, threshold=threshold)

# 高载参数计算
kol = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=a_ol, pmax=pol, pmin=0)
rol = overload_analysis.PlasticZoneWithFactor(kmax=np.array([kol]), ys=yield_strength, factor=1/math.pi)

# 以高载时裂纹长度aol筛选出高载后的数据
dadn_ola, dk_ola, n_ola, a_ola = \
    mts_analysis.FCGDataSelectByThreshold(dadn=dadn_Manual, dk=dk_Manual, n=n_Manual, a=a_Manual,
                                          threshold=a_fcgr_min, target='a', keepbigger=1)

# 高载段参数计算
kmax_ola = np.array([dk/(1-stress_ratio) for dk in dk_ola])
rm_ola = overload_analysis.PlasticZoneWithFactor(kmax=kmax_ola, ys=yield_strength, factor=1/math.pi)

# Wheeler模型指数系数m1拟合
m1, _ = overload_analysis.WheelerFittingBaseParis(a=a_ola, dadn=dadn_ola, dk=dk_ola, rm=rm_ola, aol=a_ol, rol=rol,
                                                  c=c_ca, m=m_ca)

# 原始Wheeler模型计算
alist = np.arange(a_fcgr_max, max(a_Manual), 0.02)
dadn_wheeler, dk_wheeler, a_wheeler, cp_wheeler = \
    overload_analysis.WheelerCalculatingBaseParis(a_wheeler=alist,
                                                  b=thickness, w=width, ys=yield_strength, pmax=ppeak, r=stress_ratio,
                                                  aol=a_ol, rol=rol, plasticzonefactor=1/math.pi,
                                                  m1=m1, c=c_ca, m=m_ca)

# 筛选出自dadN最高点往后的数据
dadn_old, dk_old, n_old, a_old = \
    mts_analysis.FCGDataSelectByThreshold(dadn=dadn_Manual, dk=dk_Manual, n=n_Manual, a=a_Manual,
                                          threshold=a_fcgr_max, target='a', keepbigger=1)

# 计算延迟结束位置（即dadN最低值对应的裂纹长度）对应的SIF
kd = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=a_fcgr_max+ad, pmax=ppeak, pmin=ppeak*stress_ratio)

# 以alhpa_d为系数的延迟区域尺寸计算
alpha_d = ad*1e-3 / ((kol/yield_strength*1e-3)**2 - (kd/yield_strength*1e-3)**2)
rd_ol = overload_analysis.PlasticZoneWithFactor(kmax=[kol], ys=yield_strength, factor=alpha_d)
rd_i = overload_analysis.PlasticZoneWithFactor(kmax=dk_old/(1-stress_ratio), ys=yield_strength, factor=alpha_d)

# 按照已拟合得到的Wheeler模型计算实验点对应的Retardation Parameter
cp_old = overload_analysis.CpCalculatingBaseParis(a=a_old, b=thickness, w=width, ys=yield_strength, pmax=ppeak,
                                                  r=stress_ratio, aol=a_ol, rol=rol, plasticzonefactor=1/math.pi, m1=m1)

# 拟合延迟系数m_mod
m_mod, _ = overload_analysis.DelayFittingBaseParis(a=a_old, dadn=dadn_old, dk=dk_old, cp=cp_old,
                                                   rd=rd_i, aol=a_ol_applied, rdol=rd_ol, c=c_ca, m=m_ca)

# Salvati修正的Wheeler模型拟合
dadn_salvati, dk_salvati, a_salvati, cd_salvati, cp_salvati = \
    overload_analysis.SalvatiWheelerCalculatingBaseParis(astart=a_fcgr_max, afinal=max(a_Manual),
                                                         b=thickness, w=width, ys=yield_strength,
                                                         pmax=ppeak, r=stress_ratio, aol=a_ol, rol=rol,
                                                         plasticzonefactor=1/math.pi, m1=m1, c=c_ca, m=m_ca,
                                                         apre=a_ol_applied, rdol=rd_ol, delayfactor=alpha_d, m_mod=m_mod)

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
if figsave:
    plt.savefig('WheelerModel_' + sequence[0] + '.png', dpi=320)
if show:
    plt.show()
