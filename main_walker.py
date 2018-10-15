# CT件实验裂纹扩展实验参数预测函数
# 主要内容：
# (1)根据E647-11标准对SIF变幅进行估算；
# (2)根据Paris公式对裂纹扩展速率（FCGR）进行预测；
# (3)采用Simpson插值积分公式由FCGR估算裂纹扩展速率；
# (4)根据E647-11标准对韧带长度有效性要求确定有效的裂纹扩展长度
# mm for length, MPa.m0.5 for SIFs

# 禁止改动！论文使用.

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from FCGAnalysisLib import read_data
from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import write_data

show = 0
filesave = 0
sequence = np.array(["yang-baoban_Lu-420-01", "yang-baoban_Lu-420-02", "yang-baoban_Lu-420-11", "yang-baoban_Lu-420-13",
                     "yang-baoban_Lu-420-20"])

dadn1, cycles1, dk1 = read_data.ReadMtsResult(sequence=sequence[0])
c1, m1 = mts_analysis.ParisFitting(dadn=dadn1, dk=dk1)
print("Specimen:", sequence[0], ",Paris:c=", str(c1), ",m=", str(m1))
# 读#01数据并拟合Paris

dadn2, cycles2, dk2 = read_data.ReadMtsResult(sequence=sequence[1])
c2, m2 = mts_analysis.ParisFitting(dadn=dadn2, dk=dk2)
print("Specimen:", sequence[1], ",Paris:c=", str(c2), ",m=", str(m2))
# 读#02数据并拟合Paris

dadn3, cycles3, dk3 = read_data.ReadMtsResult(sequence=sequence[2])
c3, m3 = mts_analysis.ParisFitting(dadn=dadn3, dk=dk3)
print("Specimen:", sequence[2], ",Paris:c=", str(c3), ",m=", str(m3))
# 读#11数据，筛去DK<15的部分，拟合Paris

threshold = 15.0
dadn4_temp, cycles4_temp, dk4_temp = read_data.ReadMtsResult(sequence=sequence[3])
dadn4 = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk4_temp, data=dadn4_temp)
cycles4 = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk4_temp, data=cycles4_temp)
dk4 = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk4_temp, data=dk4_temp)
c4, m4 = mts_analysis.ParisFitting(dadn=dadn4, dk=dk4)
print("Specimen:", sequence[3], ",Paris:c=", str(c4), ",m=", str(m4))
# 读#13数据，筛去DK<15的部分，拟合Paris

dadn5, cycles5, dk5 = read_data.ReadMtsResult(sequence=sequence[4])
c5, m5 = mts_analysis.ParisFitting(dadn=dadn5, dk=dk5)
print("Specimen:", sequence[4], ",Paris:c=", str(c5), ",m=", str(m5))
# 读#20数据并拟合Paris

dadn_r01 = np.array(dadn1)
dk_r01 = np.array(dk1)
c_r01, m_r01 = mts_analysis.ParisFitting(dadn=dadn_r01, dk=dk_r01)
print("r = 0.1:", ",Paris:c=", str(c_r01), ",m=", str(m_r01))
# 合并#01和#02数据(r=0.1)并拟合Paris


dadn_r05 = np.array(dadn4)
dk_r05 = np.array(dk4)
# dk_r05 = np.concatenate((dk3, dk4))
c_r05, m_r05 = mts_analysis.ParisFitting(dadn=dadn_r05, dk=dk_r05)
print("r = 0.5:", ",Paris:c=", str(c_r05), ",m=", str(m_r05))
# 合并#11和#13数据(r=0.5)并拟合Paris


dadn_r03 = np.array(dadn5)
dk_r03 = np.array(dk5)
c_r03, m_r03 = mts_analysis.ParisFitting(dadn=dadn_r03, dk=dk_r03)
print("r = 0.3:", ",Paris:c=", str(c_r03), ",m=", str(m_r03))
# 合并#20数据(r=0.3)并拟合Paris


c0b, m0b, gammab = mts_analysis.WalkerFittingByRegression(dadn1=dadn_r01, dk1=dk_r01, r1=0.1,
                                                          dadn2=dadn_r05, dk2=dk_r05, r2=0.5,
                                                          dadn3=dadn_r03, dk3=dk_r03, r3=0.3)
print("Walker's model by Regression: c0 =", str(c0b), ",m0 =", str(m0b), ",gamma=", str(gammab))
# 直接由dadn和dk回归得Walker模型参数
dadn_walker = []
dk_walker, r_walker = np.mgrid[13:46:3, 0.1:0.8:0.2]
dadn_walker = mts_analysis.WalkerCalculating(c0=c0b, m0=m0b, gamma=gammab, dk=dk_walker, r=r_walker)
dk3d = []
r3d = []
dadn3d = []
for dkt in dk_walker:
    for rt in r_walker:
        for dki in dkt:
            for ri in rt:
                dadni = mts_analysis.WalkerCalculating(c0=c0b, m0=m0b, gamma=gammab, dk=dki, r=ri)
                dk3d.append(dki)
                r3d.append(ri)
                dadn3d.append(dadni)
dk3d = np.array(dk3d)
r3d = np.array(r3d)
dadn3d = np.array(dadn3d)
# Walker模型计算平面

dk_walker_r01 = np.arange(np.min(dk_r01), np.max(dk_r01), 1)
dadn_walker_r01 = mts_analysis.WalkerCalculating(c0=c0b, m0=m0b, gamma=gammab, dk=dk_walker_r01, r=0.1)
dk_walker_r05 = np.arange(np.min(dk_r05), np.max(dk_r05), 1)
dadn_walker_r05 = mts_analysis.WalkerCalculating(c0=c0b, m0=m0b, gamma=gammab, dk=dk_walker_r05, r=0.5)
dk_walker_r03 = np.arange(np.min(dk_r03), np.max(dk_r03), 1)
dadn_walker_r03 = mts_analysis.WalkerCalculating(c0=c0b, m0=m0b, gamma=gammab, dk=dk_walker_r03, r=0.3)
# 用Walker模型分别拟合R=0.1,0.5,0.7的情况

dk_eq_r01 = dk_r01*(1 - 0.1)**(gammab - 1)
dk_eq_r05 = dk_r05*(1 - 0.5)**(gammab - 1)
dk_eq_r03 = dk_r03*(1 - 0.3)**(gammab - 1)
dk_walker_eq = np.arange(min(dk_eq_r05), max(dk_eq_r01), 1)
dadn_walker_eq = c0b*dk_walker_eq**m0b
# 用Walker模型的等效SIF表示结果


# Plotting: (1)da/dN - dk - r 3D plot
r01 = np.full(len(dk_r01), 0.1)
r05 = np.full(len(dk_r05), 0.5)
r03 = np.full(len(dk_r03), 0.3)
fig = plt.figure(figsize=(7, 5), dpi=320, num=1)
ax = plt.subplot(111, projection='3d')
ax.scatter(np.log(dk_r01), np.log(1-r01), np.log(dadn_r01), s=1, label='$r=0.1$', color='red', lw=1)
ax.scatter(np.log(dk_r05), np.log(1-r05), np.log(dadn_r05), s=1, label='$r=0.5$', color='orange', lw=1)
ax.scatter(np.log(dk_r03), np.log(1-r03), np.log(dadn_r03), s=1, label='$r=0.3$', color='blue', lw=1)
Axes3D.plot_surface(self=ax, X=np.log(dk_walker), Y=np.log(1-r_walker), Z=np.log(dadn_walker), rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
ax.set_xlabel("log(DeltaK Applied)/MPa*m^0.5")
ax.set_ylabel("log(r)")
ax.set_zlabel("log(da/dN)/mm per cycle")
if show:
    plt.show()

# file save:
if filesave:
    data2 = np.array([1-r3d, dk3d, dadn3d])
    name2 = ["walker_1-r", "walker_dk", "walker_dadn"]
    _ = write_data.SaveData(dataset=data2, name=name2)

# Plotting: (2)dadN - dk 2D plot in dk
if show:
    plt.figure(num=1, figsize=(10, 8))
    plt.scatter(dk_r01, dadn_r01, lw=1, marker='+', label='R=0.1')
    plt.scatter(dk_r05, dadn_r05, lw=1, marker='2', label='R=0.5')
    # plt.scatter(dk_r07, dadn_r07, lw=1, marker='*', label='R=0.7')
    plt.plot(dk_walker_r01, dadn_walker_r01, label='WalkerFit R=0.1', color='black', linewidth=2)
    plt.plot(dk_walker_r05, dadn_walker_r05, label='WalkerCalculation R=0.5', color='blue', linewidth=2)
    # plt.plot(dk_walker_r07, dadn_walker_r07, label='WalkerFit R=0.7', color='red', linewidth=2)
    plt.axis([min(dk_r05), max(dk_r01), min(dadn_r05), max(dadn_r01) * 1.2])
    plt.yticks(np.linspace(min(dadn_r05) * 0.9, max(dadn_r01), 6))
    plt.xticks(np.linspace(min(dk_r05), max(dk_r01), 6))
    plt.title("FCG Rates - Delta SIF(Real)")
    plt.ylabel("FCG Rates/ mm per cycle")
    plt.xlabel("DeltaK/ MPa.m0.5")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.grid(which='minor', linestyle='--')
    plt.grid(which='major', linestyle='--')
    plt.show()

# Plotting: (2)dadN - dk 2D plot
if show:
    plt.figure(num=3, figsize=(10, 8))
    plt.scatter(dk_eq_r01, dadn_r01, lw=1, marker='+', label='R=0.1')
    plt.scatter(dk_eq_r05, dadn_r05, lw=1, marker='2', label='R=0.5')
    # plt.scatter(dk_eq_r07, dadn_r07, lw=1, marker='*', label='R=0.7')
    plt.plot(dk_walker_eq, dadn_walker_eq, label='WalkerFit', color='black', linewidth=2)
    plt.axis([min(dk_eq_r01), max(dk_eq_r01), min(dadn_r01), max(dadn_r01) * 1.2])
    plt.yticks(np.linspace(min(dadn_r05) * 0.9, max(dadn_r01), 6))
    plt.xticks(np.linspace(min(dk_eq_r05), max(dk_eq_r01), 6))
    plt.title("FCG Rates - Delta SIF(Equivalent)")
    plt.ylabel("FCG Rates/ mm per cycle")
    plt.xlabel("DeltaKeq/ MPa.m0.5")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.grid(which='minor', linestyle='--')
    plt.grid(which='major', linestyle='--')
    plt.show()
