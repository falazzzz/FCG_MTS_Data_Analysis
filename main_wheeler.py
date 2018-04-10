import pandas as pd
import numpy as np
import mts_analysis
import math
import read_data
import experiment_calculation
import matplotlib.pyplot as plt


def ReadTensionResult(sequence):
    # usage:读取简单拉伸实验结果数据
    # input parameter:
    # sequence: 试件序号
    # return parameter:
    # extension：位移/mm
    # load：载荷/N
    # stress：工程应力/MPa
    # strain：工程应变/MPa
    tensionresult = pd.read_csv(u"QSTE420TM_TensionTest_" + sequence + u".csv")  # Data File Reading
    extension = tensionresult["Extension (mm)"]
    load = tensionresult["Load (N)"]
    stress = tensionresult["Stress (MPa)"]
    strain = tensionresult["Strain (%)"] * 0.01
    return extension, load, stress, strain


def TrueStress(engineeringstress, engineeringstrain):
    # usage：由工程应力和工程应变计算真实应力
    truestress = engineeringstress * (1 + engineeringstrain)
    return truestress


def TrueStrain(engineeringstrain):
    # usage：由工程应变计算真实应变
    truestrain = np.log(1 + engineeringstrain)
    return truestrain


def RambergOsgoodFitting(s, e, ys, E):
    # usage: 由真实应力和真实应变拟合Ramberg-Osgood本构模型
    # input parameter:
    # s: 真实应力/MPa
    # e：真实应变
    # ys：屈服强度/GPa
    # E：弹性模量/GPa
    # return parameter:
    # n: Ramberg-Osgood模型参数（硬化指数）
    # alpha: Ramberg-Osgood模型参数
    ys = ys * 1e3
    E = E * 1e3
    Y0 = np.array([(E/ys)*strain for strain in e])
    Y1 = np.array([stress/ys for stress in s])
    Y = np.log(Y0 + Y1)
    del Y0, Y1
    X = np.log([stress/ys for stress in s])
    p = np.polyfit(X, Y, deg=1)
    n, B = p
    alpha = np.exp(B)
    return n, alpha


def RambergOsgoodCalculation(n, alpha, s, ys, E):
    # usage: 由已经得到的Ramberg-Osgood参数和需要拟合的真实应力值，以该本构关系计算真实应变
    # input parameter:
    # n，alpha：Ramberg-Osgood模型参数
    # s：真实应力/MPa
    # ys：屈服强度/GPa
    # E：弹性模量/GPa
    # return parameter:
    # e：由Ramberg-Osgood本构关系计算得到的真实应变
    e = []
    ys = ys * 1e3
    E = E * 1e3
    for stress in s:
        left = stress/ys + alpha*(stress/ys)**n
        e.append(left*ys/E)
    return np.array(e)


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

sequence = ["b-1#", "yang-baoban_Lu-420-06"]
stress_ratio = 0.1
threshold = 18
c_ca = 7.9713e-10
m_ca = 3.6797
# 参数设置

l, f, s_engineering, e_engineering = ReadTensionResult(sequence=sequence[0])
specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence[1])
cycles, cracklength, kmax, kmin, pmax, pmin = read_data.ReadOriginResult(sequence[1], closure=False, cod=False)
dadn_Manual, n_Manual, dk_Manual, a_Manual = \
    experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin,
                                                   a=cracklength,
                                                   ys=yield_strength, r=stress_ratio, threshold=threshold)
e_true1 = TrueStrain(e_engineering)
s_true1 = TrueStress(s_engineering, e_engineering)
e_true = mts_analysis.DataSelectByThreshold(threshold=0.25, parameter=e_true1, data=e_true1, keepbigger=0)
s_true = mts_analysis.DataSelectByThreshold(threshold=0.25, parameter=e_true1, data=s_true1, keepbigger=0)
e_true0 = mts_analysis.DataSelectByThreshold(threshold=446, parameter=s_true, data=e_true, keepbigger=1)
s_true0 = mts_analysis.DataSelectByThreshold(threshold=446, parameter=s_true, data=s_true, keepbigger=1)
n, alpha = RambergOsgoodFitting(s=s_true0, e=e_true0, ys=yield_strength, E=elastic_modulus)
# 筛选简单拉伸数据的拟合部分，仅采用塑性段进行拟合
print("Ramberg-Osgood Fitting Result:n=%.4f" % n , ", alpha=%.4f" % alpha)
slist = np.arange(0, 725, 25)
elist = RambergOsgoodCalculation(n=n, alpha=alpha, s=slist, ys=yield_strength, E=elastic_modulus)
print(yield_strength)
# 拟合Ramberg-Osgood模型，得到n

factor_plasticzone = (n - 1)/(n + 1)
# 塑形区尺寸因子

aol = 12.0057
pol = 4000
kol = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=aol, pmax=pol, pmin=0)
rol = PlasticZoneWithFactor(kmax=np.array([kol]), ys=yield_strength, factor=factor_plasticzone)
# 高载值计算

olrange = [20.55, 27.15]
dadn_temp = mts_analysis.DataSelectByThreshold(threshold=olrange[0], parameter=dk_Manual, data=dadn_Manual, keepbigger=1)
dk_temp = mts_analysis.DataSelectByThreshold(threshold=olrange[0], parameter=dk_Manual, data=dk_Manual, keepbigger=1)
n_temp = mts_analysis.DataSelectByThreshold(threshold=olrange[0], parameter=dk_Manual, data=n_Manual, keepbigger=1)
a_temp = mts_analysis.DataSelectByThreshold(threshold=olrange[0], parameter=dk_Manual, data=a_Manual, keepbigger=1)
dadn_ol = mts_analysis.DataSelectByThreshold(threshold=olrange[1], parameter=dk_temp, data=dadn_temp, keepbigger=0)
dk_ol = mts_analysis.DataSelectByThreshold(threshold=olrange[1], parameter=dk_temp, data=dk_temp, keepbigger=0)
n_ol = mts_analysis.DataSelectByThreshold(threshold=olrange[1], parameter=dk_temp, data=n_temp, keepbigger=0)
a_ol = mts_analysis.DataSelectByThreshold(threshold=olrange[1], parameter=dk_temp, data=a_temp, keepbigger=0)
# 以DeltaK筛选出高载段的属性


kmax_ol = np.array([dk/(1-stress_ratio) for dk in dk_ol])
rm = PlasticZoneWithFactor(kmax=kmax_ol, ys=yield_strength, factor=1/math.pi)
# 高载段参数计算

left = a_ol + rm
right = a_ol + rol
for seq, _ in enumerate(left):
    if left[seq] >= right[seq]:
        dadn_ol = dadn_ol[0:seq]
        dk_ol = dk_ol[0:seq]
        n_ol = n_ol[0:seq]
        a_ol = a_ol[0:seq]
        print("Invalid Results Detected and Deleted.")
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
print(m1, b)
# Wheeler模型 m1系数拟合

cp = []
for seq, _ in enumerate(rm):
    cpi = (rm[seq]/(aol + rol - a_ol[seq]))**m1
    cp.append(cpi)
cp = np.array(cp)
dadn_wheeler = []
for seq, _ in enumerate(dk_ol):
    dadn = cp[seq]*(c_ca*dk_ol[seq]**m_ca)
    dadn_wheeler.append(dadn)
dadn_wheeler = np.array(dadn_wheeler)
# Wheeler模型 拟合函数计算

# Plotting
plt.figure(num=2, figsize=(10, 8))
plt.scatter(e_engineering, s_engineering, lw=1, marker='+', label='Engineering')
plt.scatter(e_true, s_true, lw=1, marker='*', label='True')
plt.plot(elist, slist, label='Ramberg-Osgood Fitting', color='black', linewidth=2)
plt.title("Stress - Strain Curve")
plt.ylabel("Stress / MPa")
plt.xlabel("Strain")
plt.legend()
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.show()

plt.figure(num=3, figsize=(10, 8))
plt.scatter(dk_ol, dadn_ol, lw=1, marker='+', label='Experiment')
plt.plot(dk_ol, dadn_wheeler, label='Wheeler Model Fitting', color='black', linewidth=2)
plt.title("FCG Rates - deltaSIF(Wheeler Model)")
plt.ylabel("FCG Rates/mm per cycle")
plt.xlabel("DeltaSIF/MPa.m0.5")

plt.xscale('log')
plt.yscale('log')
plt.axis([min(dk_ol), max(dk_ol), min(dadn_wheeler), max(dadn_wheeler)*1.2])
plt.yticks(np.linspace(min(dadn_wheeler)*0.8, max(dadn_wheeler), 6))
plt.xticks(np.linspace(min(dk_ol), max(dk_ol), 6))
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()
plt.show()
