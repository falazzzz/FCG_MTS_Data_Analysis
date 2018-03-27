import read_data
import mts_analysis
import matplotlib.pyplot as plt
import numpy as np

# main1: for origin data dispose, drawing figures

# global unit
# length: mm
# dadn: mm/cycle
# SIF: MPa.m0.5
# stress, modulus: GPa
# Load：N

#Switch
show = 1
save = 0
sequence = "yang-baoban_Lu-420-01"         # Graph Saving Sequence
cracklength_for_overload = 14              # For Finding Part, Find the cycles when length = 10mm
dk_for_overload = 25                       # For Finding Part, Find the cycles when SIF = 30 MPa.m0.5

# Temp Parameter
stress_ratio = 0.1     # Stress Ratio

dadn_MTS1, cycles_MTS1, dk_MTS1 = read_data.ReadMtsResult(sequence=sequence)
# mm for length, MPa.m0.5 for SIFs

cycles, cracklength, kmax, kmin, pmax, pmin, codmax, codmin, = \
    read_data.ReadOriginResult(sequence=sequence)
# mm for length, MPa.m0.5 for SIFs

specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence)
# mm for length, GPa for strength and modulus

sequence = "yang-baoban_Lu-420-02"         # Graph Saving Sequence

# Temp Parameter
stress_ratio = 0.1     # Stress Ratio

dadn_MTS2, cycles_MTS2, dk_MTS2 = read_data.ReadMtsResult(sequence=sequence)
# mm for length, MPa.m0.5 for SIFs

cycles2, cracklength2, kmax2, kmin2, pmax2, pmin2, codmax2, codmin2, = \
    read_data.ReadOriginResult(sequence=sequence)
# mm for length, MPa.m0.5 for SIFs

specimen2, width2, notch_length2, thickness2, elastic_modulus2, yield_strength2, precrack_length2 = \
    read_data.ReadTestInf(sequence=sequence)
# mm for length, GPa for strength and modulus

dadn_MTS = np.concatenate((dadn_MTS1, dadn_MTS2))
dk_MTS = np.concatenate((dk_MTS1, dk_MTS2))

c_MTS, m_MTS = mts_analysis.ParisFitting(dadn=dadn_MTS, dk=dk_MTS)
print("MTS Results Fitting Result by Paris:c=", str(c_MTS), ",m=", str(m_MTS))
dk_paris_sorted = sorted(dk_MTS)
dadn_paris = mts_analysis.ParisCalculating(c=c_MTS, m=m_MTS, dk=dk_paris_sorted)

sequence = 'R=0.1_QSTE420TM'

# Plotting: da/dN - dk plot(2次实验结果对比)
plt.figure(num=3, figsize=(7, 5))
plt.scatter(dk_MTS1, dadn_MTS1, s=1, label='$Pmax=2.0kN$', color='red', lw=1)
plt.scatter(dk_MTS2, dadn_MTS2, s=1, label='$Pmax=1.6kN$', color='blue', lw=1)
plt.plot(dk_paris_sorted, dadn_paris, label='$Fitting By Paris$', color='black', linewidth=2)
plt.xlabel("DeltaK Applied (MPa*m^0.5)")
plt.xscale('log')
plt.ylabel("da/dN (mm/cycle)")
plt.yscale('log')
plt.title('da/dN - dK ' + sequence + '(MTS Result)')
plt.legend()
plt.grid()
if save:
    plt.savefig('dadn_dk_MTS' + sequence + '.png', dpi=320)
if show:
    plt.show()
