# MTS单组实验数据处理函数
# 适用采用MTS内置Fatigue Crack Growth程序实验得到的数据格式
# 读取的文件包括原始数据（CrackLength等），MTS处理得到FCGR-DK数据，以及Test Summary
# 主要内容：
# (1)读取文件中的相应数据
# (2)将MTS的数据进行Paris拟合，得到参数c_MTS, m_MTS
# (3)读取原始数据中的裂纹长度，依据E647-11计算SIF变幅，采用割线法计算FCGR
# (4)根据E647-11中对韧带长度有效性要求筛去不符合的数据，根据指定的阈值（目前为无）筛去不符合的数据，得到满足要求的FCGR和DK
# (5)将由原始数据手动计算的数据进行Paris拟合，得到参数c_manual，m_Manual
# (6)根据指定的裂纹长度和SIF变幅，分别依据原始数据和MTS数据推算相应的循环数
# (7)绘图
# 全局单位设定：
# length: mm
# dadn: mm/cycle
# SIF: MPa.m0.5
# stress, modulus: GPa
# Load：N

# Updatad on 2018/5/30
# 代码完全重构，添加了数据输出功能

from FCGAnalysisLib import mts_analysis
import matplotlib.pyplot as plt
from FCGAnalysisLib import read_data
from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import write_data
from FCGAnalysisLib import paris_and_walker
import numpy as np

'''
参数输入
'''
# 实验基本参数
sequence = "yang-baoban_Lu-420-19"         # Graph Saving Sequence
r = 0.3

# 计算参数
threshold = 0                              # 最小保留数据阈值

# 控制参数
fitting = 0
show = 1
filesave = 0
figsave = 0

# 材料Paris参数和Walker参数
c_w, m_w, gamma = paris_and_walker.WalkerParameter()
c_p, m_p = paris_and_walker.ParisParameter(r=r)

'''
数据读取处理
'''
# 基本数据读取和计算
specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence)
dadn_MTS, n_MTS, dk_MTS = read_data.ReadMtsResult(sequence=sequence)
cycles, cracklength, kmax, kmin, pmax, pmin, kclosure, closureload = \
    read_data.ReadOriginResult(sequence=sequence, closure=True, cod=False)
dadn_Manual, n_Manual, dk_Manual, a_Manual, kc_Manual = \
    experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin,
                                                   a=cracklength,
                                                   ys=yield_strength, r=r, kclosure=kclosure,
                                                   threshold=threshold)

'''
Paris拟合和Walker模型效果对比
'''
if fitting:
    c_p, m_p = mts_analysis.ParisFitting(dadn=dadn_Manual, dk=dk_Manual)
    print("c=", c_p, "m=", m_p)
dadn_paris_Manual = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=dk_Manual)
dadn_walker_Manual = mts_analysis.WalkerCalculating(c0=c_w, m0=m_w, gamma=gamma, r=r, dk=dk_Manual)

'''
数据保存
'''
data = np.array([n_Manual, a_Manual, dk_Manual, dadn_Manual,
                 kc_Manual,
                 ])
name = ["Original_Cycles", "Original_CrackLength/mm", "Original_SIFRange/MPam0.5", "Original_FCGRate/mm per cycle",
        "Original_ClosureSIF/MPam0.5",
        ]
if filesave:
    _ = write_data.SaveData(dataset=data, name=name, filename=sequence)

'''
绘图命令
'''
# 评估移动平均对dadN的修正效果
plt.figure(num=1, figsize=(10, 8))
plt.scatter(dk_Manual, dadn_Manual, label='$ExperimentData$', color='red', marker='.')
# plt.scatter(dk_MTS, dadn_MTS, label='$ExperimentDataByMTS$', color='orange', marker='.')
plt.plot(dk_Manual, dadn_paris_Manual, label='$Paris for CA$', color='black', linewidth=2)
# plt.plot(dk_Manual, dadn_walker_Manual, label='$Walker MOdel$', color='blue', linewidth=2)
plt.axis([min(dk_Manual) * 0.95, max(dk_Manual) * 1.05, min(dadn_Manual) * 0.95, max(dadn_Manual) * 1.05])
plt.xlabel("DeltaK Applied (MPa*m^0.5)")
plt.xscale('log')
plt.ylabel("da/dN (mm/cycle)")
plt.yscale('log')
# plt.xticks(np.linspace(min(dk_Manual), max(dk_Manual), 6))
plt.title('da/dN - dK ' + sequence)
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()
plt.grid()
if figsave:
    plt.savefig(filename='dadNEvaluate'+sequence)

plt.figure(num=2, figsize=(10, 8))
plt.scatter(n_Manual, a_Manual, label='$ExperimentData$', color='red', marker='.')
plt.xlabel("Cycles")
plt.ylabel("Crack Length/mm")
plt.title('a - N ' + sequence)
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
if figsave:
    plt.savefig(filename='a_N_'+sequence)

# 显示图像
if show:
    plt.show()
