# 全局单位设定：
# length: mm
# dadn: mm/cycle
# SIF: MPa.m0.5
# stress, modulus: GPa
# Load：N
# Updated on 2019/1/13

from FCGAnalysisLib import mts_class
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
# 实例基本参数
sequences = ["yang-baoban_Lu-420-01", "yang-baoban_Lu-420-05", "yang-baoban_Lu-340-3-7", "yang-baoban_Lu-340-3-6", "yang-baoban_Lu-340-3-5", "yang-baoban_Lu-340-3-4", "yang-baoban_Lu-340-3-3", "yang-baoban_Lu-340-3-2", "yang-baoban_Lu-340-3-1"]
ol_1_cycle = [60041, 64006, 57172, 52371, 52014, 51758, 50234, 58414, 66728]
r = 0.1
threshold = np.full(len(sequences), 0)                              # 最小保留数据阈值(SIF Range)

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
for seq, sequence in enumerate(sequences):
    specimen_data.append(mts_class.DoubleOLSpecimenBasic(name=sequence, stress_ratio=r,
                                                         threshold=threshold[seq], ol_1_cycle=ol_1_cycle[seq]))  # 实例

for seq, specimen in enumerate(specimen_data):
    if seq == 0 or seq == 1:
        specimen.read_experimental_data(ys=365)
    else:
        specimen.read_experimental_data()     # 读取和计算
    specimen.ol_1_fix()
    max_n = np.max(specimen.n_Fixed_by_ol_1)
    print(specimen.name, max_n, (max_n - 124025)/124025)

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
    #plt.scatter(specimen.n_Manual, specimen.a_Manual, label='$ExperimentData$', marker='.')
    plt.scatter(specimen.n_Fixed_by_ol_1, specimen.a_Manual, label=specimen.name, marker='.')
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
