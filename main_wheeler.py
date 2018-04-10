import math
import matplotlib.pyplot as plt
import numpy as np
from FCGAnalysisLib import read_data
from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import experiment_calculation


def PlasticZoneWithFactor(kmax, ys, factor):
    # usage: 计算Irwin及其改进型的塑性区尺寸，表达式r = factor*(kmax/ys)^2，其中factor为塑性区系数
    # input parameter:
    # kmax：应力强度因子/MPa.m0.5
    # ys：屈服强度/GPa
    # factor：塑性区尺寸系数
    # return parameter:
    # r：塑形区尺寸/mm
    r = []
    ys = ys * 1e3
    for kmaxs in kmax:
        r.append((factor * (kmaxs/ys)**2)*1e3)
    return np.array(r)


def WheelerRetardationParameter(rm, aol, rol, a, m):
    # usage:
    # input parameter:
    # sequence:
    # return parameter:
    #
    return 0

sequence = ["yang-baoban_Lu-420-06"]
stress_ratio = 0.1
threshold = 0
c_ca = 7.9713e-10
m_ca = 3.6797
aol = 12.0057
pol = 4000
ppeak = 2000
# 参数设置

specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence[0])
cycles, cracklength, kmax, kmin, pmax, pmin = read_data.ReadOriginResult(sequence[0], closure=False, cod=False)
dadn_Manual, n_Manual, dk_Manual, a_Manual = \
    experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin,
                                                   a=cracklength,
                                                   ys=yield_strength, r=stress_ratio, threshold=threshold)

kol = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=aol, pmax=pol, pmin=0)
rol = PlasticZoneWithFactor(kmax=np.array([kol]), ys=yield_strength, factor=1/math.pi)
# 高载值计算

dadn_ola = mts_analysis.DataSelectByThreshold(threshold=aol, parameter=a_Manual, data=dadn_Manual, keepbigger=1)
dk_ola = mts_analysis.DataSelectByThreshold(threshold=aol, parameter=a_Manual, data=dk_Manual, keepbigger=1)
n_ola = mts_analysis.DataSelectByThreshold(threshold=aol, parameter=a_Manual, data=n_Manual, keepbigger=1)
a_ola = mts_analysis.DataSelectByThreshold(threshold=aol, parameter=a_Manual, data=a_Manual, keepbigger=1)
# 以DeltaK筛选出高载段的属性

kmax_ol = np.array([dk/(1-stress_ratio) for dk in dk_ola])
rm = PlasticZoneWithFactor(kmax=kmax_ol, ys=yield_strength, factor=1/math.pi)
# 高载段参数计算

left = a_ola + rm
right = aol + rol
for seq, _ in enumerate(left):
    if left[seq] >= right:
        dadn_ol = dadn_ola[0:seq]
        dk_ol = dk_ola[0:seq]
        n_ol = n_ola[0:seq]
        a_ol = a_ola[0:seq]
        rm = rm[0:seq]
        print("Wheeler Region Check: Invalid Results Detected and Deleted.")
        break
# Wheeler迟滞参数适用区域计算

X = []
Y = []
for seq, _ in enumerate(dadn_ol):
    y = np.log(dadn_ol[seq]) - m_ca*np.log(dk_ol[seq]) - np.log(c_ca)
    Y.append(y)
    x = np.log(rm[seq]) - np.log(aol + rol - a_ol[seq])
    X.append(x)
X = np.array(X).reshape(-1, )
Y = np.array(Y)
p = np.polyfit(X, Y, deg=1)
m1, b = p
print("Wheeler Model Fitting Result: m1=", m1, ",intercept=",b)
# Wheeler模型 m1系数拟合


a_wheeler = np.arange(aol, np.max(a_Manual), 0.1)
dk_wheeler = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=a_wheeler, pmax=ppeak, pmin=ppeak * stress_ratio)
kmax_wheeler = np.array([dk/(1 - stress_ratio) for dk in dk_wheeler])
rm_wheeler = PlasticZoneWithFactor(kmax=kmax_wheeler, ys=yield_strength, factor=1/math.pi)
left_wheeler = a_wheeler + rm_wheeler
right_wheeler = aol + rol
# 迟滞计算参数准备，由裂纹长度a出发

cp = []
for seq, _ in enumerate(a_wheeler):
    if left[seq] < right:
        cpi = (rm_wheeler[seq]/(aol + rol - a_wheeler[seq]))**m1
    else:
        cpi = 1
    cp.append(cpi)
cp = np.array(cp)
# 迟滞系数cp计算
dadn_wheeler = []
for seq, _ in enumerate(dk_wheeler):
    dadn = cp[seq]*(c_ca*dk_wheeler[seq]**m_ca)
    dadn_wheeler.append(dadn)
dadn_wheeler = np.array(dadn_wheeler)
# Wheeler模型 拟合函数计算

# Plotting

plt.figure(num=3, figsize=(10, 8))
plt.scatter(dk_Manual, dadn_Manual, lw=1, marker='+', label='Experiment')
plt.plot(dk_wheeler, dadn_wheeler, label='Wheeler Model Fitting', color='black', linewidth=2)
plt.title("FCG Rates - deltaSIF(Wheeler Model)")
plt.ylabel("FCG Rates/mm per cycle")
plt.xlabel("DeltaSIF/MPa.m0.5")
plt.xscale('log')
plt.yscale('log')
plt.axis([min(dk_Manual), max(dk_Manual), min(dadn_Manual)*0.9, max(dadn_Manual)*1.1])
#plt.yticks(np.linspace(min(dadn_Manual)*0.8, max(dadn_Manual), 6))
#plt.xticks(np.linspace(min(dk_Manual), max(dk_Manual), 6))
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()
plt.show()
