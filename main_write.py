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


import read_data
import mts_analysis
import experiment_calculation
import write_data


# 文件参数
sequence = "yang-baoban_Lu-420-01"         # Graph Saving Sequence
stress_ratio = 0.1                         # Stress Ratio
threshold = 0                              # Threshold for paris region selecting

# 保存参数
save = True                                # Save File or not
dkstep = 1                                 # Recording Step for MTSResult
astep = 0.5                                # Recording Step for ManualResult

# 实验基本参数读取
specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence)
print('specimen name:', str(specimen))
# mm for length, GPa for strength and modulus


# MTS数据读取和处理
dadn_MTS, cycles_MTS, dk_MTS = read_data.ReadMtsResult(sequence=sequence)
c_MTS, m_MTS = mts_analysis.ParisFitting(dadn=dadn_MTS, dk=dk_MTS)
print("MTS Results Fitting Result by Paris:c=", str(c_MTS), ",m=", str(m_MTS))


# 由裂纹长度、载荷、循环次数计算
cycles, cracklength, kmax, kmin, pmax, pmin, codmax, codmin, = \
    read_data.ReadOriginResult(sequence=sequence)

dadn_Manual, n_Manual, dk_Manual, a_Manual = \
    experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin, a=cracklength,
                                                   ys=yield_strength, r=stress_ratio, threshold=threshold)

c_Manual, m_Manual = mts_analysis.ParisFitting(dadn=dadn_Manual, dk=dk_Manual)
print("Manual Restlt Fitting by Paris: c=", str(c_Manual), ",m=", str(m_Manual))


# 结果文件输出
MTSResult = write_data.ArrangeData(dadn=dadn_MTS, cycles=cycles_MTS, dk=dk_MTS,
                                   option="dk", step=dkstep, save=save, name="MTSResult_"+sequence)
ManualResult = write_data.ArrangeData(dadn=dadn_Manual, cycles=n_Manual, dk=dk_Manual, a=a_Manual,
                                      option="a", step=astep, save=save, name="ManualResult_"+sequence)
