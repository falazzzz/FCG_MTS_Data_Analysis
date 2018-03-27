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


import read_data
import mts_analysis
import plot
import matplotlib.pyplot as plt


# 分析和绘图参数
show = 0
save = 0
sequence = "yang-baoban_Lu-420-02"         # Graph Saving Sequence
cracklength_for_overload = 14              # For Finding Part, Find the cycles when length = 10mm
dk_for_overload = 25                       # For Finding Part, Find the cycles when SIF = 30 MPa.m0.5
stress_ratio = 0.1                         # Stress Ratio

dadn_MTS, cycles_MTS, dk_MTS = read_data.ReadMtsResult(sequence=sequence)
# mm for length, MPa.m0.5 for SIFs

cycles, cracklength, kmax, kmin, pmax, pmin, codmax, codmin, = \
    read_data.ReadOriginResult(sequence=sequence)
# mm for length, MPa.m0.5 for SIFs

specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence)
# mm for length, GPa for strength and modulus

print('specimen name:',str(specimen))

c_MTS, m_MTS = mts_analysis.ParisFitting(dadn=dadn_MTS, dk=dk_MTS)
# Fitting Paris by MTS Result, returning Paris Parameters C and m
print("MTS Results Fitting Result by Paris:c=", str(c_MTS), ",m=", str(m_MTS))

dadn_paris_by_MTS = mts_analysis.ParisCalculating(c=c_MTS, m=m_MTS, dk=dk_MTS)
dadn_reference = mts_analysis.ParisCalculating(c=4.4668e-8, m=2.3671, dk=dk_MTS)
# Calculate MTS Result by Paris(Fitting), Walker(Fitting), Paris(Reference)


# Original Data Analysis

cracklength_complience_max = mts_analysis.Compliance(e=elastic_modulus, b=thickness, w=width, p=pmax, v=codmax)
cracklength_complience_min = mts_analysis.Compliance(e=elastic_modulus, b=thickness, w=width, p=pmin, v=codmin)
cracklength_compliance = cracklength_complience_min     # 若取裂纹平均值改为(mina + maxa)/2
# Calculate crack length by complience, return the max as the cracklength
# 内部处理时采用柔度法计算裂纹长度时，考虑了90%/10%的Upper和Lower区间，此处处理无法实现，以最大最小进行计算误差较大.

dk_compliance = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=cracklength, pmax=pmax, pmin=pmin)
# Calculate delta K

dadn_secant, cracklength_secant, cycles_secant, dk_secant = \
    mts_analysis.FCGRateBySecant(a=cracklength, n=cycles, dk=dk_compliance)
# Calculate FCG Rate by Secant, data selected at the same time, negative results are discarded.

cycles_ligament_valid = \
    mts_analysis.DataSelectByLigament(w=width, a=cracklength_secant, dk=dk_secant, ys=yield_strength, data=cycles_secant, r=stress_ratio)
cracklength_ligament_valid = \
    mts_analysis.DataSelectByLigament(w=width, a=cracklength_secant, dk=dk_secant, ys=yield_strength, data=cracklength_secant, r=stress_ratio)
dadn_ligament_valid = \
    mts_analysis.DataSelectByLigament(w=width, a=cracklength_secant, dk=dk_secant, ys=yield_strength, data=dadn_secant, r=stress_ratio)
dk_ligament_valid = \
    mts_analysis.DataSelectByLigament(w=width, a=cracklength_secant, dk=dk_secant, ys=yield_strength, data=dk_secant, r=stress_ratio)
print("Ligament Vaild Check,", str(len(cycles_secant) - len(cycles_ligament_valid)), "records deleted.")
# Check by Ligament Validation

threshold = int(cycles.tail(1))
cycles_threshold_valid = \
    mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=cycles_ligament_valid, data=cycles_ligament_valid)
cracklength_threshold_valid = \
    mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=cycles_ligament_valid, data=cracklength_ligament_valid)
dadn_threshold_valid = \
    mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=cycles_ligament_valid, data=dadn_ligament_valid)
dk_threshold_valid = \
    mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=cycles_ligament_valid, data=dk_ligament_valid)
print("Threshold Vaild Check,", str(len(cycles_ligament_valid) - len(cycles_threshold_valid)), "records deleted.")
# Selected by Threshold，不做筛选（选择本身的最大值进行筛选，相当于不做筛选）

c_Manual, m_Manual = mts_analysis.ParisFitting(dadn=dadn_threshold_valid, dk=dk_threshold_valid)
print("Manual Restlt Fitting by Paris: c=", str(c_Manual), ",m=", str(m_Manual))
dadn_paris_by_Manual = mts_analysis.ParisCalculating(c=c_Manual, m=m_Manual, dk=dk_threshold_valid)
# Fitting Manual Result by Paris


# Finding Information for Overload
cycle_for_overload1 = int(mts_analysis.FindAscentDataBySeq(value=cracklength_for_overload, item=cracklength, target=cycles))
print("Predicting Cycle when Crack Length =", str(cracklength_for_overload), "mm:", str(cycle_for_overload1))
# Finding Cycles by CrackLength(by Origin Data)
cycle_for_overload2 = int(mts_analysis.FindAscentDataBySeq(value=dk_for_overload, item=dk_MTS, target=cycles_MTS))
print("Predicting Cycle when DeltaK =", str(dk_for_overload), "Mpa.m0.5:", str(cycle_for_overload2))
# Finding Cycles by CrackLength(by MTS Data)



# Ploting
plot.CrackLength_Cycles(sequence=sequence,
                        n_MTS=cycles,
                        a_MTS=cracklength,
                        save=save,
                        show=show)

plot.CrackLengthError_Cycles(sequence=sequence,
                             a_Manual=cracklength_compliance,
                             a_MTS=cracklength,
                             n=cycles,
                             save=save,
                             show=show)

plot.FCGR_DeltaK_MTS(sequence=sequence,
                     dk=dk_MTS,
                     dadn=dadn_MTS,
                     dadn_paris=dadn_paris_by_MTS,
                     save=save,
                     show=show)

plot.FCGR_DeltaK_Comparation(sequence=sequence,
                             dk_MTS=dk_MTS,
                             dadn_MTS=dadn_MTS,
                             dk_Manual=dk_threshold_valid,
                             dadn_manual=dadn_threshold_valid,
                             dadn_paris_MTS=dadn_paris_by_MTS,
                             dadn_paris_Manual=dadn_paris_by_Manual,
                             save=save,
                             show=show)

if show == 0:
    plt.close()