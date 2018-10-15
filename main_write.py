# MTS单组实验数据筛选保存函数
# 适用采用MTS内置Fatigue Crack Growth程序实验得到的数据格式
# 读取的文件包括原始数据（CrackLength等），MTS处理得到FCGR-DK数据，以及Test Summary
# 主要内容：
# (1)读取文件中的相应数据
# (2)将MTS的数据进行Paris拟合，得到参数c_MTS, m_MTS
# (3)读取原始数据中的裂纹长度，依据E647-11计算SIF变幅，采用割线法计算FCGR
# (4)根据E647-11中对韧带长度有效性要求筛去不符合的数据，根据指定的阈值（目前为无）筛去不符合的数据，得到满足要求的FCGR和DK
# (5)将由原始数据手动计算的数据进行Paris拟合，得到参数c_manual，m_Manual
# (6)按照指定的保存步长step，将Manual数据和MTS数据保存为csv文件

# 全局单位设定：
# length: mm
# dadn: mm/cycle
# SIF: MPa.m0.5
# stress, modulus: GPa
# Load：N


import numpy as np
from FCGAnalysisLib import closureanalysis
from FCGAnalysisLib import write_data

# 文件参数
sequence = "yang-baoban_Lu-420-20"         # Graph Saving Sequence
stress_ratio = 0.3                         # Stress Ratio
filename = sequence + "_basic_inf"


# 实验基本参数读取
dadn_Manual, n_Manual, dk_Manual, a_Manual, kop_Manual = closureanalysis.BasicDataDispose(sequence=sequence, r=stress_ratio)
print('specimen name:', sequence)

# 结果文件输出
data = np.array([a_Manual, n_Manual, dk_Manual, dadn_Manual, kop_Manual])
name = ['crack length', 'cycles', 'dK', 'dadn', 'K_open']
_ = write_data.SaveData(dataset=data, name=name, filename=filename)

