import read_data
import mts_analysis
import plot
import matplotlib.pyplot as plt

# global unit
# length: mm
# dadn: mm/cycle
# SIF: MPa.m0.5
# stress, modulus: GPa
# Load：N

#Switch
show = 0
save = 1
sequence = "yang-baoban_Lu-420-02"         # Graph Saving Sequence

# Temp Parameter
stress_ratio = 0.1     # Stress Ratio

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

c0_MTS, m0_MTS, gamma_MTS = mts_analysis.WalkerFitting(dadn=dadn_MTS, dk=dk_MTS, r=stress_ratio)
# Fitting Walker's Model by MTS Result, returning Walker Parameters C0, m0 adn gamma
print("MTS Results Fitting by Walker's Model:c0 = ", str(c0_MTS), ',gamma=', str(gamma_MTS), ',m=', str(m0_MTS))

dadn_paris_by_MTS = mts_analysis.ParisCalculating(c=c_MTS, m=m_MTS, dk=dk_MTS)
dadn_walker_by_MTS = mts_analysis.WalkerCalculating(c0=c0_MTS, m0=m0_MTS, gamma=gamma_MTS, dk=dk_MTS, r=stress_ratio)
dadn_reference = mts_analysis.ParisCalculating(c=4.4668e-8, m=2.3671, dk=dk_MTS)
# Calculate MTS Result by Paris(Fitting), Walker(Fitting), Paris(Reference)

# Original Data Analysis

cracklength_complience_max = mts_analysis.Compliance(e=elastic_modulus, b=thickness, w=width, p=pmax, v=codmax)
cracklength_complience_min = mts_analysis.Compliance(e=elastic_modulus, b=thickness, w=width, p=pmin, v=codmin)
cracklength_compliance = cracklength_complience_max     # 若取裂纹平均值改为(mina + maxa)/2
# Calculate crack length by complience, return the max as the cracklength


dk_compliance = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=cracklength_compliance, pmax=pmax, pmin=pmin)
# Calculate delta K

dadn_secant, cracklength_secant, cycles_secant, dk_secant = \
    mts_analysis.FCGRateBySecant(a=cracklength_compliance, n=cycles, dk=dk_compliance)
# Calculate FCG Rate by Secant, data selected at the same time, negative results are discarded.

cycles_ligament_valid = \
    mts_analysis.DataSelectByLigament(w=width, a=cracklength_secant, dk=dk_secant, ys=yield_strength, data=cycles_secant, r=stress_ratio)
cracklength_ligament_valid = \
    mts_analysis.DataSelectByLigament(w=width, a=cracklength_secant, dk=dk_secant, ys=yield_strength, data=cracklength_secant, r=stress_ratio)
dadn_ligament_valid = \
    mts_analysis.DataSelectByLigament(w=width, a=cracklength_secant, dk=dk_secant, ys=yield_strength, data=dadn_secant, r=stress_ratio)
dk_ligament_valid = \
    mts_analysis.DataSelectByLigament(w=width, a=cracklength_secant, dk=dk_secant, ys=yield_strength, data=dk_secant, r=stress_ratio)
# Check by Ligament Validation, 函数需要检查！从上一版留下的毛病！

threshold = int(cycles_MTS.tail(1))
cycles_threshold_valid = \
    mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=cycles_ligament_valid, data=cycles_ligament_valid)
cracklength_threshold_valid = \
    mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=cycles_ligament_valid, data=cracklength_ligament_valid)
dadn_threshold_valid = \
    mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=cycles_ligament_valid, data=dadn_ligament_valid)
dk_threshold_valid = \
    mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=cycles_ligament_valid, data=dk_ligament_valid)
# Selected by Threshold，选取了MTS认定有效的循环数以内的数据

c_Manual, m_Manual = mts_analysis.ParisFitting(dadn=dadn_threshold_valid, dk=dk_threshold_valid)
print("Manual Restlt Fitting by Paris: c=", str(c_Manual), ",m=", str(m_Manual))
dadn_paris_by_Manual = mts_analysis.ParisCalculating(c=c_Manual, m=m_Manual, dk=dk_threshold_valid)

# 绘图部分
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