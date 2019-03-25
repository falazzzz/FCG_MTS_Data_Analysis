# 2018/7/2 Version 3.0
# Crcak Closure模型计算（包括basic、2/PI法和alpha法，程序源于main_crackclosure）
# 本函数用于alpha模型输出。通过给定衰减系数和幅值（亦可给0使得通过拟合得到），并输出曲线图像及具体值文件（csv格式）
# Lu Yunchao

import matplotlib.pyplot as plt
import numpy as np

from FCGAnalysisLib import closureanalysis
from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import paris_and_walker
from FCGAnalysisLib import write_data

'''
参数输入
'''
# 实验基本参数
sequence = "yang-baoban_Lu-420-19"         # Graph Saving Sequence

# 计算参数
threshold = 0                              # 最小保留数据阈值
numofaverage = 7                           # 移动平均的个数
ratio = 0.2                                # 误差带半宽度
fade_rate = 7.2                           # alpha模型参数（改为0时将通过拟合得到）
amplitude = 2/np.pi                        # alpha模型参数（改为0时将通过拟合得到）

# 显示参数
show = 1                                   # 图像显示开关
datasave = 1                               # 文件保存开关，保存时文件将以temp+保存时间命名。

'''
数据预处理
'''

# alpha模型及实验参数
r, overload_ratio, a_range = paris_and_walker.testinf(sequence)     # 由试件名从数据库调应力比、高载比及拟合范围
aol = a_range[0]-0.4                # alpha模型拟合起点，取k_closure的最低值
aol_end = a_range[1]            # alpha模型拟合终点，可通过该值筛除后段不稳定的数据

# 材料Paris参数和Walker参数
c_w, m_w, gamma = paris_and_walker.WalkerParameter()
c_p, m_p = paris_and_walker.ParisParameter(r=r)


# 基本数据读取和计算
dadn_Manual, n_Manual, dk_Manual, a_Manual, kc_Manual = closureanalysis.BasicDataDispose(sequence=sequence, r=r)

'''
2pi法
'''
# dkeff计算
dkeff_Manual_2pi = closureanalysis.dKeffCalculation(kmax=dk_Manual/(1-r), kclose=kc_Manual, kmin=dk_Manual*(r/(1-r)),
                                                    method='2/PI')
# 移动平均
dkeff_Manual_averaged0, dk_Manual_averaged0 = \
    closureanalysis.DoubleDataSelectMovingAverage(data=dkeff_Manual_2pi, reference=dk_Manual,
                                                  ratio=ratio, numofaverage=numofaverage)
_, a_Manual_averaged0 = \
    closureanalysis.DoubleDataSelectMovingAverage(data=dkeff_Manual_2pi, reference=a_Manual,
                                                  ratio=ratio, numofaverage=numofaverage)
# 计算基于裂纹闭合效应的dadn(2/PI法)
dadn_paris_dkeff_2pi = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=dkeff_Manual_averaged0)

'''
传统dkeff法(basic)
'''
# dkeff计算
dkeff_Manual_bac = closureanalysis.dKeffCalculation(kmax=dk_Manual/(1-r), kclose=kc_Manual, kmin=dk_Manual*(r/(1-r)),
                                                    method='basic')
# 移动平均
dkeff_Manual_averaged1, dk_Manual_averaged1 = \
    closureanalysis.DoubleDataSelectMovingAverage(data=dkeff_Manual_bac, reference=dk_Manual,
                                                  ratio=ratio, numofaverage=numofaverage)
_, a_Manual_averaged1 = \
    closureanalysis.DoubleDataSelectMovingAverage(data=dkeff_Manual_bac, reference=a_Manual,
                                                  ratio=ratio, numofaverage=numofaverage)
# 计算基于裂纹闭合效应的dadn(basic法)
dadn_paris_dkeff_bac = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=dkeff_Manual_averaged1)
# 计算不进行移动平均的dadn(basic法)
dadn_paris_dkeff_noaverage = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=dkeff_Manual_bac)

'''
alpha Method
'''
# Alpha法API
D, n, dk_selected_averaged2, dkeff_alpha_e, dadn_alpha_e, alpha_e = \
    closureanalysis.Alpha2(sequence=sequence, stress_ratio=r, error_ratio=ratio,
                           numofaverage=numofaverage, a_range=[aol, aol_end],
                           fade_rate=fade_rate, amplitude=amplitude, model='walker')
U = dkeff_alpha_e/dk_selected_averaged2

'''
数据筛选
'''
dadn_paris_dkeff_2pi_selected = mts_analysis.DataSelectByThreshold(threshold=0, parameter=a_Manual_averaged0,
                                                                   data=dadn_paris_dkeff_2pi)
dadn_paris_dkeff_2pi_selected = np.array(dadn_paris_dkeff_2pi_selected)
dk_Manual_averaged0_selected = mts_analysis.DataSelectByThreshold(threshold=0, parameter=a_Manual_averaged0,
                                                                  data=dk_Manual_averaged0)
dk_Manual_averaged0_selected = np.array(dk_Manual_averaged0_selected)
dadn_paris_dkeff_bac_selected = mts_analysis.DataSelectByThreshold(threshold=0, parameter=a_Manual_averaged1,
                                                                   data=dadn_paris_dkeff_bac)
dadn_paris_dkeff_bac_selected = np.array(dadn_paris_dkeff_bac_selected)
dk_Manual_averaged1_selected = mts_analysis.DataSelectByThreshold(threshold=0, parameter=a_Manual_averaged1,
                                                                  data=dk_Manual_averaged1)
dk_Manual_averaged1_selected = np.array(dk_Manual_averaged1_selected)

'''
绘图命令
'''
# 常幅循环(CA)下Parsi的dadn
dadn_paris_Manual = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=dk_Manual)

# 评估移动平均对dkeff的修正效果和Basic、2PI法的效果
plt.figure(num=2, figsize=(6, 6))
plt.scatter(dk_Manual, dadn_Manual, label='$Experiment$', marker='.')
plt.plot(dk_Manual, dadn_paris_Manual, label='$Paris for CA$', color='black', linewidth=2)
#plt.scatter(dk_Manual, dadn_paris_dkeff_noaverage, label='$CrackClosure Without Average$', marker='1')
plt.scatter(dk_Manual_averaged0_selected, dadn_paris_dkeff_2pi_selected, label='$CrackClosure 2/PI Method$', marker='*')
plt.scatter(dk_Manual_averaged1_selected, dadn_paris_dkeff_bac_selected, label='$CrackClosure Basic Method$', marker='*')
plt.scatter(dk_selected_averaged2, dadn_alpha_e, label='$CrackClosure Alpha Method$', marker='2')
plt.axis([min(dk_Manual)*0.9, max(dk_Manual),
          min(min(dadn_paris_Manual),min(dadn_paris_dkeff_bac))*0.9, max(dadn_Manual)*1.1])
plt.yticks(np.linspace(min(min(dadn_paris_Manual), min(dadn_paris_dkeff_bac))*0.9, max(dadn_Manual), 6))
plt.xticks(np.linspace(min(dk_Manual)*0.9, max(dk_Manual), 6))
plt.title("FCG Rates - Delta SIF_"+sequence)
plt.ylabel("FCG Rates/ mm per cycle")
plt.xlabel("DeltaK/ MPa.m0.5")
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')

# 显示图像
if show:
    plt.show()

'''
文件输出
'''
if datasave:
    data = np.array([dk_Manual, dadn_Manual,  dk_selected_averaged2, dadn_alpha_e,
                     dk_Manual_averaged1_selected, dadn_paris_dkeff_bac_selected, a_Manual_averaged1,
                     dk_Manual_averaged0_selected, dadn_paris_dkeff_2pi_selected, a_Manual_averaged0,
                     dk_Manual, dadn_paris_Manual])
    name = ['dK_exp', 'dadn_exp', 'dk_alpha', 'dadn_alpha',
            'dk_bac', 'dadn_bac', 'a_bac',
            'dk_2pi', 'dadn_2pi', 'a_2pi',
            'dk_CA', 'dadn_CA']
    _ = write_data.SaveData(dataset=data, name=name)