# 2018/7/2 Version 2.0
# 本函数用于多组OL实验绘图，不含计算过程
# Lu Yunchao

from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import paris_and_walker
from FCGAnalysisLib import closureanalysis
import numpy as np
import matplotlib.pyplot as plt

# 实验参数
sequence = ["yang-baoban_Lu-420-13", "yang-baoban_Lu-420-15", "yang-baoban_Lu-420-16",
            "yang-baoban_Lu-420-07", "yang-baoban_Lu-420-10"]
stress_ratio = [0.5, 0.5, 0.5, 0.5, 0.5]

# 程序参数
show = 1            # 绘图开关

# paris参数读取，以stress_ratio中第一组的应力比为准
c_CA, m_CA = paris_and_walker.ParisParameter(r=stress_ratio[0])

# seqence中数据的读取和处理，未对SIFclosure进行读取(记为-)
dadn_Manual1, n_Manual1, dk_Manual1, a_Manual1, _ = \
    closureanalysis.BasicDataDispose(sequence=sequence[0], r=stress_ratio[0])

dadn_Manual2, n_Manual2, dk_Manual2, a_Manual2, _ = \
    closureanalysis.BasicDataDispose(sequence=sequence[1], r=stress_ratio[1])

dadn_Manual3, n_Manual3, dk_Manual3, a_Manual3, _ = \
    closureanalysis.BasicDataDispose(sequence=sequence[2], r=stress_ratio[2])

dadn_Manual4, n_Manual4, dk_Manual4, a_Manual4, _ = \
    closureanalysis.BasicDataDispose(sequence=sequence[3], r=stress_ratio[3])

dadn_Manual5, n_Manual5, dk_Manual5, a_Manual5, _ = \
    closureanalysis.BasicDataDispose(sequence=sequence[4], r=stress_ratio[4])

# CA曲线计算
dk_paris = np.linspace(15, 25, 9)
dadn_paris = mts_analysis.ParisCalculating(c=c_CA, m=m_CA, dk=dk_paris)

plt.figure(num=1, figsize=(10, 8))
plt.scatter(dk_Manual1, dadn_Manual1, lw=1, marker='+', label='CA')
plt.scatter(dk_Manual2, dadn_Manual2, lw=1, marker='*', label='OLR=1.5')
plt.scatter(dk_Manual3, dadn_Manual3, lw=1, marker='2', label='OLR=2.0')
#plt.scatter(dk_Manual4, dadn_Manual4, lw=1, marker='3', label='OLR=2.5')
#plt.scatter(dk_Manual5, dadn_Manual5, lw=1, marker='1', label='OLR=4.0')
plt.plot(dk_paris, dadn_paris, label='CA_Paris', color='black', linewidth=2)
plt.axis([min(dk_Manual1), max(dk_Manual1), min(dadn_Manual1), max(dadn_Manual1)*1.2])
plt.yticks(np.linspace(min(dadn_Manual3)*0.5, max(dadn_Manual1), 6))
plt.xticks(np.linspace(10, 30, 6))
plt.title("FCG Rates - Delta SIF")
#plt.title("CrackLength - Cycles")
plt.ylabel("FCG Rates/ mm per cycle")
#plt.ylabel("CrackLength / mm")
plt.xlabel("DeltaK/ MPa.m0.5")
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
if show:
    plt.show()


plt.figure(num=2, figsize=(10, 8))
plt.scatter(n_Manual1, a_Manual1, lw=1, marker='+', label='CA')
plt.scatter(n_Manual2, a_Manual2, lw=1, marker='*', label='OLR=1.5')
plt.scatter(n_Manual3, a_Manual3, lw=1, marker='2', label='OLR=2.0')
#plt.scatter(n_Manual4, a_Manual4, lw=1, marker='3', label='OLR=2.5')
#plt.scatter(n_Manual5, a_Manual5, lw=1, marker='1', label='OLR=4.0')
plt.axis([0, max(n_Manual3)*1.1, min(a_Manual1), max(a_Manual1)*1.1])
plt.yticks(np.linspace(min(a_Manual3), max(a_Manual1)*1.1, 6))
plt.xticks(np.linspace(0, max(n_Manual3)*1.1, 6))
plt.title("CrackLength - Cycles")
plt.ylabel("CrackLength / mm")
plt.xlabel("Cycle")
#plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
if show:
    plt.show()

