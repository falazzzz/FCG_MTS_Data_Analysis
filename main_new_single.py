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
# Updated on 2019/1/13 with new class method

from FCGAnalysisLib import mts_class
from FCGAnalysisLib import mts_analysis
import matplotlib.pyplot as plt
from FCGAnalysisLib import write_data
from FCGAnalysisLib import paris_and_walker
import numpy as np

'''
参数输入
'''
# 实例基本参数
sequences = ["yang-baoban_Lu-420-05"]         # Graph Saving Sequence
r = 0.1
threshold = [0]                              # 最小保留数据阈值(SIF Range)

# 控制参数
fitting = 0
show = 1
filesave = 0
figsave = 0

# 材料Paris参数
c_p, m_p = paris_and_walker.ParisParameter(r=r)

'''
数据读取处理
'''
specimen_data = []
for sequence in sequences:
    specimen_data.append(mts_class.SpecimenBasic(name=sequence, stress_ratio=r, threshold=0))   # 创建实例

for specimen in specimen_data:
    specimen.read_experimental_data()     # 读取和计算

'''
Paris拟合和Walker模型效果对比
'''
dadn_paris_Manual = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=np.array(specimen.dk_Manual))

'''
数据保存
'''
if filesave:
    for specimen in specimen_data:
        data = np.array([specimen.n_Manual, specimen.a_Manual, specimen.dk_Manual, specimen.dadn_Manual,
                         specimen.kc_Manual])
        name = ["Original_Cycles", "Original_CrackLength/mm", "Original_SIFRange/MPam0.5",
                "Original_FCGRate/mm per cycle", "Original_ClosureSIF/MPam0.5"]
    _ = write_data.SaveData(dataset=data, name=name, filename=specimen.name)

'''
绘图命令
'''
# 评估移动平均对dadN的修正效果
plt.figure(num=1, figsize=(10, 8))
for specimen in specimen_data:
    plt.scatter(specimen.dk_Manual, specimen.dadn_Manual, label=specimen.name, marker='.')
    plt.plot(specimen.dk_Manual, specimen.dadn_Manual)
    #plt.scatter(specimen.dk_MTS, specimen.dadn_MTS, label='$ExperimentDataByMTS$', color='orange', marker='.')
plt.plot(specimen.dk_Manual, dadn_paris_Manual, label='$Paris for CA$', color='black', linewidth=2)
# plt.plot(dk_Manual, dadn_walker_Manual, label='$Walker Model$', color='blue', linewidth=2)
plt.axis([min(specimen.dk_Manual) * 0.95, max(specimen.dk_Manual) * 1.05,
          min(specimen.dadn_Manual) * 0.95, max(specimen.dadn_Manual) * 1.05])
plt.xlabel("DeltaK Applied (MPa*m^0.5)")
plt.xscale('log')
plt.ylabel("da/dN (mm/cycle)")
plt.yscale('log')
plt.title('da/dN - dK ')
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
if figsave:
    plt.savefig(filename='dadNEvaluate'+sequence)

plt.figure(num=2, figsize=(10, 8))
for specimen in specimen_data:
    plt.scatter(specimen.n_Manual, specimen.a_Manual, label=specimen.name, marker='.')
plt.xlabel("Cycles")
plt.ylabel("Crack Length/mm")
plt.title('a - N ')
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
if figsave:
    plt.savefig(filename='a_N_'+sequence)

# 显示图像
if show:
    plt.legend()
    plt.show()
