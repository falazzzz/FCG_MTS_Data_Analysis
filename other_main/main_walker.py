# CT件实验裂纹扩展实验参数预测函数
# 主要内容：
# (1)根据E647-11标准对SIF变幅进行估算；
# (2)根据Paris公式对裂纹扩展速率（FCGR）进行预测；
# (3)采用Simpson插值积分公式由FCGR估算裂纹扩展速率；
# (4)根据E647-11标准对韧带长度有效性要求确定有效的裂纹扩展长度
# mm for length, MPa.m0.5 for SIFs

import read_data
import mts_analysis
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

show = 1
save = 0
sequence = np.array(["yang-baoban_Lu-420-01","yang-baoban_Lu-420-02","yang-baoban_Lu-420-03","yang-baoban_Lu-420-04"])

dadn1, cycles1, dk1 = read_data.ReadMtsResult(sequence=sequence[0])
c1, m1 = mts_analysis.ParisFitting(dadn=dadn1, dk=dk1)
print("Specimen:", sequence[0], ",Paris:c=", str(c1), ",m=", str(m1))
# 读#01数据并拟合Paris

dadn2, cycles2, dk2 = read_data.ReadMtsResult(sequence=sequence[1])
c2, m2 = mts_analysis.ParisFitting(dadn=dadn2, dk=dk2)
print("Specimen:", sequence[1], ",Paris:c=", str(c2), ",m=", str(m2))
# 读#02数据并拟合Paris

threshold = 15
dadn3_temp, cycles3_temp, dk3_temp = read_data.ReadMtsResult(sequence=sequence[2])
dadn3 = dadn3_temp
cycles3 =cycles3_temp
dk3 = dk3_temp
'''
dadn3 = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk3_temp, data=dadn3_temp)
cycles3 = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk3_temp, data=cycles3_temp)
dk3 = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk3_temp, data=dk3_temp)
'''
c3, m3 = mts_analysis.ParisFitting(dadn=dadn3, dk=dk3)
print("Specimen:", sequence[2], ",Paris:c=", str(c3), ",m=", str(m3))
# 读#03数据，筛去DK<15的部分，拟合Paris

dadn4, cycles4, dk4 = read_data.ReadMtsResult(sequence=sequence[3])
c4, m4 = mts_analysis.ParisFitting(dadn=dadn4, dk=dk4)
print("Specimen:", sequence[3], ",Paris:c=", str(c4), ",m=", str(m4))
# 读#04数据并拟合Paris

dadn_r01 = np.concatenate((dadn1, dadn2))
dk_r01 = np.concatenate((dk1, dk2))
c_r01, m_r01 = mts_analysis.ParisFitting(dadn=dadn_r01, dk=dk_r01)
print("r = 0.1:", ",Paris:c=", str(c_r01), ",m=", str(m_r01))
# 合并#01和#02数据(r=0.1)并拟合Paris


dadn_r07 = np.concatenate((dadn3, dadn4))
dk_r07 = np.concatenate((dk3, dk4))
c_r07, m_r07 = mts_analysis.ParisFitting(dadn=dadn_r07, dk=dk_r07)
print("r = 0.7:", ",Paris:c=", str(c_r07), ",m=", str(m_r07))
# 合并#03和#04数据(r=0.7)并拟合Paris


c0a, m0a, gammaa = mts_analysis.WalkerFittingFromParis(c=[c_r01,c_r07], m=[m_r01, m_r07], r=[0.1, 0.7])
print("Walker's model From Paris: c0 =", str(c0a), ",m0 =",str(m0a), ",gamma=", str(gammaa))
# 由Paris参数拟合Walker模型参数
c0b, m0b, gammab = mts_analysis.WalkerFittingByRegression(dadn1=dadn_r01, dk1=dk_r01, r1=0.1, dadn2=dadn_r07, dk2=dadn_r07, r2=0.7)
print("Walker's model by Regression: c0 =", str(c0b), ",m0 =", str(m0b), ",gamma=", str(gammab))
# 直接由dadn和dk回归得Walker模型参数

dadn_walker = []
dk_walker, r_walker = np.mgrid[10:41:1, 0:0.8:0.1]
dadn_walker = mts_analysis.WalkerCalculating(c0=c0a, m0=m0a, gamma=gammaa, dk=dk_walker, r=r_walker)

# Walker模型计算



# Plotting: (1)da/dN - dk - r (MTS) plot
r01 = np.full(len(dk_r01), 0.1)
r07 = np.full(len(dk_r07), 0.7)
fig = plt.figure(figsize=(7, 5), dpi=320, num=1)
ax = plt.subplot(111, projection='3d')
ax.scatter(dk_r01, r01, dadn_r01, s=1, label='$r=0.1$', color='red', lw=1)
ax.scatter(dk_r07, r07, dadn_r07, s=1, label='$r=0.7$', color='blue', lw=1)
Axes3D.plot_surface(self=ax, X=dk_walker, Y=r_walker, Z=dadn_walker, rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
ax.set_xlabel("DeltaK Applied (MPa*m^0.5)")
ax.set_ylabel("r")
ax.set_zlabel("da/dN (mm/cycle)")
ax.set_zlim(1e-6, 1e-3)
if save:
    plt.savefig('Walker.png', dpi=320)
if show:
    plt.show()

