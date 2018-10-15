# Willenborg模型的实现程序，效果一般
# 2018/10/15 Version 1.0
# Lu Yunchao

from FCGAnalysisLib import read_data
from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import overload_analysis
from FCGAnalysisLib import paris_and_walker
import numpy as np
import matplotlib.pyplot as plt

# 实验参数
sequence = ["yang-baoban_Lu-420-06"]         # Graph Saving Sequence
pol = 4000
r = 0.1
threshold = 0
aol_hold = 12.0057  # 高载长度：按Hold时记录的裂纹长度
aol_maxrate = 12.0057   # 高载长度：对应FCGRate最大时的裂纹长度
aol_minrate = 12.7927   # 高载长度：对应FCGRate最小时的裂纹长度

# 程序参数
show = 1            # 绘图开关

# 裂纹扩展公式数据库
c_w, m_w, gamma = paris_and_walker.WalkerParameter()
c_p, m_p = paris_and_walker.ParisParameter(r=r)

# 读取数据并完成基本的FCG处理
specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence[0])
cycles, cracklength, kmax, kmin, pmax, pmin, kclosure, closureload = \
    read_data.ReadOriginResult(sequence=sequence[0], closure=True, cod=False)
dadn_Manual, n_Manual, dk_Manual, a_Manual, kc_Manual = \
    experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin,
                                                   a=cracklength,
                                                   ys=yield_strength, r=r, kclosure=kclosure,
                                                   threshold=threshold)

# Kmax at Overload计算
k_max_ol = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=aol_hold, pmax=pol, pmin=0)

# 塑性区系数及尺寸计算
factor = overload_analysis.FactorXiaoping(kmax=k_max_ol, t=thickness, ys=yield_strength)
r_ol = overload_analysis.PlasticZoneWithFactor(kmax=[k_max_ol], ys=yield_strength, factor=factor)

# 计算dK_eff
dkeff = []
dk_willenborg = []
for seq, value in enumerate(a_Manual):
    # 只计算裂纹长度扩展到aol_hold以后的数据的dk_eff(未施加高载以前的不进行计算)
    if value > aol_minrate:
        k_red = k_max_ol*np.sqrt(1 - (a_Manual[seq]-aol_minrate)*1e-3/r_ol) - dk_Manual[seq]/(1-r)
        kmax_eff = dk_Manual[seq]/(1-r) - k_red
        if dk_Manual[seq]*r/(1-r) > k_red:
            kmin_eff = dk_Manual[seq]*r/(1-r) - k_red
        else:
            kmin_eff = 0
        dkeffseq = kmax_eff - kmin_eff
        dkeff.append(dkeffseq)
        dk_willenborg.append(dk_Manual[seq])
        print(k_red, kmax_eff, kmin_eff, dkeffseq, dk_Manual[seq])
dkeff = np.array(dkeff)
dk_willenborg = np.array(dk_willenborg)

# 采用Walker模型计算
dadn_willenborg = mts_analysis.WalkerCalculating(c0=c_w, m0=m_w, gamma=gamma, dk=dkeff, r=r)

plt.figure(num=1, figsize=(10, 8))
plt.scatter(dk_Manual, dadn_Manual, lw=1, marker='+', label='Experiment')
plt.scatter(dk_willenborg, dadn_willenborg, lw=1, marker='*', label='Willenborg Fitting')
plt.title("FCG Rates - deltaK(Willenborg Model)")
plt.ylabel("FCG Rates/mm per cycle")
plt.xlabel("DeltaSIF/MPa.m0.5")
plt.xscale('log')
plt.yscale('log')
plt.axis([min(dk_Manual), max(dk_Manual)*1.2, min(dadn_Manual)*0.5, max(dadn_Manual)*1.5])
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()
if show:
    plt.show()




