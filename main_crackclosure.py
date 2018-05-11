from FCGAnalysisLib import mts_analysis
import matplotlib.pyplot as plt
from FCGAnalysisLib import read_data
from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import overload_analysis
from FCGAnalysisLib import closureanalysis
import numpy as np
import math
from scipy import optimize

sequence = "yang-baoban_Lu-420-06"         # Graph Saving Sequence
a_ol = 12.0116
pol = 4000
pi = math.pi
r = 0.1
threshold = 0
specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence)
cycles, cracklength, kmax, kmin, pmax, pmin, kclosure, closureload = \
    read_data.ReadOriginResult(sequence=sequence, closure=True, cod=False)
dadn_Manual, n_Manual, dk_Manual, a_Manual, kc_Manual= \
    experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin,
                                                   a=cracklength,
                                                   ys=yield_strength, r=r, kclosure=kclosure,
                                                   threshold=threshold)
dkeff_Manual = closureanalysis.dKeffCalculation(kmax=dk_Manual/(1-r), kclose=kc_Manual, kmin=dk_Manual*(r/(1-r)),
                                                method='2/PI')
# 读取结果，计算有效DK，采用2/PI方法

c_w = 8.0448e-10
m_w = 3.6655
gamma = 0.9055
if r == 0.1:
    c_p = 7.9713e-10
    m_p = 3.6797
    load = 2000 * (1 - r)
elif r == 0.7:
    c_p = 5.6665e-09
    m_p = 3.0965
    load = 4000 * (1 - r)
elif r == 0.5:
    c_p = 0.9627e-09
    m_p = 3.7073
    load = 2800 * (1 - r)
dadn_paris_dkeff = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=dkeff_Manual)
dadn_walker_dkeff = mts_analysis.WalkerCalculating(c0=c_w, m0=m_w, gamma=gamma, dk=dkeff_Manual, r=r)
# 计算裂纹闭合理论得到的dadn结果

beta = []
beta_w = []
for seq, dk in enumerate(a_Manual):
    dkeqeff0 = (dadn_paris_dkeff[seq] / c_w) ** (1 / m_w)
    dkeqeff1 = (dadn_walker_dkeff[seq] / c_w) ** (1 / m_w)
    dkeq0 = dk_Manual[seq]*(1 - r)**(gamma - 1)
    beta.append((dkeq0 - dkeqeff0)/dkeqeff0)
    beta_w.append((dkeq0 - dkeqeff1)/dkeqeff1)
# 计算beta值

beta_afterol_, dkeff_afterol_, dk_afterol_, a_afterol_ = \
    mts_analysis.FCGDataSelectByThreshold(dadn=beta, dk=dkeff_Manual, n=dk_Manual, a= a_Manual,
                                          threshold=a_ol, target='a')
beta_afterol, dkeff_afterol, dk_afterol , a_afterol = \
    mts_analysis.FCGDataSelectByThreshold(dadn=beta_afterol_, dk=dkeff_afterol_, n=dk_afterol_, a=a_afterol_,
                                          threshold=a_ol + 8, target='a', keepbigger=0)
# 筛去高载以前和末尾段不稳定的数据

X = a_afterol - a_ol
Y = beta_afterol

def residuals(p):
    k, miu, sigma = p
    return Y - (k/(sigma*X*math.sqrt(2*pi)))*np.exp(-(np.log(X) - miu)**2/(2*sigma**2))

result = optimize.leastsq(residuals, [1, 1, 1])
k, miu, sigma = result[0]
print("Regression Result: ", k, miu, sigma)
# 用最小二乘法进行对数正态分布拟合

aplot = np.arange(0.001, 8, 0.001)
betafit = np.array([(k/(sigma*a*math.sqrt(2*pi)))*np.exp(-(np.log(a) - miu)**2/(2*sigma**2)) for a in aplot])
# 拟合曲线计算

kol = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=a_ol, pmax=pol, pmin=0)
xiaoping = 0.35 - 0.29/(1+(1.08*kol**2/(thickness*yield_strength**2))**2.15)
plasticzone = overload_analysis.PlasticZoneWithFactor(kmax=[kol], ys=yield_strength, factor=xiaoping)
print("Plastic Zone Size with Xiaoping Factor: ", str(plasticzone[0]), "mm.")
# 计算高载产生的塑性区尺寸

a_predict = np.arange(a_ol+0.02, 18.02, 0.02)
load = np.full(len(a_predict), load)
# 利用beta拟合dadn曲线部分，给出a的数组，其余参数均在程序头部给出了

dk_predict = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=a_predict, pmax=load, pmin=0)
# 计算对应的dk

dkeq_predict = dk_predict*(1 - r)**(gamma - 1)
# 由公式得到dkeq

aafterol_predict = a_predict - a_ol
beta_predict = [(k/(sigma*a*math.sqrt(2*pi)))*np.exp(-(np.log(a) - miu)**2/(2*sigma**2)) for a in aafterol_predict]
# 由对数正态分布拟合式计算beta

dkeqeff_predict = np.array([dkeq_predict[seq]/(1+beta_predict[seq]) for seq, _ in enumerate(dkeq_predict)])
# 由beta的表达式计算dkeq_eff

dadn_predict = mts_analysis.WalkerCalculating(c0=c_w, m0=m_w, gamma=gamma, dk=dkeqeff_predict, r=r)
# 由Walker模型得到dadn

############
sequence = "yang-baoban_Lu-420-06"         # Graph Saving Sequence
a_ol2 = 12.3495
pol = 4000
pi = math.pi
r = 0.1
threshold = 0
specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence)
cycles, cracklength, kmax, kmin, pmax, pmin, kclosure, closureload = \
    read_data.ReadOriginResult(sequence=sequence, closure=True, cod=False)
dadn_Manual2, n_Manual2, dk_Manual2, a_Manual2, kc_Manual2= \
    experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin,
                                                   a=cracklength,
                                                   ys=yield_strength, r=r, kclosure=kclosure,
                                                   threshold=threshold)
dkeff_Manual2 = dk_Manual2/(1 - r) - kc_Manual2
# 读取结果，计算有效DK

dadn_paris_dkeff2 = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=dkeff_Manual2)
dadn_walker_dkeff2 = mts_analysis.WalkerCalculating(c0=c_w, m0=m_w, gamma=gamma, dk=dkeff_Manual2, r=r)
# 计算裂纹闭合理论得到的dadn结果

beta2 = []
beta_w2 = []
for seq, dk in enumerate(a_Manual2):
    dkeqeff02 = (dadn_paris_dkeff2[seq] / c_w) ** (1 / m_w)
    dkeqeff12 = (dadn_walker_dkeff2[seq] / c_w) ** (1 / m_w)
    dkeq02 = dk_Manual2[seq]*(1 - r)**(gamma - 1)
    beta2.append((dkeq02 - dkeqeff02)/dkeqeff02)
    beta_w2.append((dkeq02 - dkeqeff12)/dkeqeff12)
# 计算beta值

beta_afterol_, dkeff_afterol_, dk_afterol_, a_afterol_ = \
    mts_analysis.FCGDataSelectByThreshold(dadn=beta2, dk=dkeff_Manual2, n=dk_Manual2, a=a_Manual2,
                                          threshold=a_ol2 + 0, target='a')
beta_afterol2, dkeff_afterol2, dk_afterol2 , a_afterol2 = \
    mts_analysis.FCGDataSelectByThreshold(dadn=beta_afterol_, dk=dkeff_afterol_, n=dk_afterol_, a=a_afterol_,
                                          threshold=a_ol2 + 8, target='a', keepbigger=0)
# 筛去高载以前和末尾段不稳定的数据

X = a_afterol2 - a_ol2
Y = beta_afterol2

def residuals(p):
    k, miu, sigma = p
    return Y - (k/(sigma*X*math.sqrt(2*pi)))*np.exp(-(np.log(X) - miu)**2/(2*sigma**2))

result = optimize.leastsq(residuals, [1, 1, 1])
k, miu, sigma = result[0]
print("Regression Result: ", k, miu, sigma)
# 用最小二乘法进行对数正态分布拟合

aplot2 = np.arange(0.001, 8, 0.001)
betafit2 = np.array([(k/(sigma*a*math.sqrt(2*pi)))*np.exp(-(np.log(a) - miu)**2/(2*sigma**2)) for a in aplot2])
# 拟合曲线计算

kol2 = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=a_ol2, pmax=pol, pmin=0)
xiaoping2 = 0.35 - 0.29/(1+(1.08*kol2**2/(thickness*yield_strength**2))**2.15)
plasticzone2 = overload_analysis.PlasticZoneWithFactor(kmax=[kol2], ys=yield_strength, factor=2)
print("Plastic Zone Size with Xiaoping Factor: ", str(plasticzone2[0]), "mm.")
# 计算高载产生的塑性区尺寸

a_predict2 = np.arange(a_ol2+0.02, 18.02, 0.02)
load2 = np.full(len(a_predict2), min(load))
# 利用beta拟合dadn曲线部分，给出a的数组，其余参数均在程序头部给出了

dk_predict2 = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=a_predict2, pmax=load2, pmin=0)
# 计算对应的dk

dkeq_predict2 = dk_predict2*(1 - r)**(gamma - 1)
# 由公式得到dkeq

aafterol_predict2 = a_predict2 - a_ol2
beta_predict2 = [(k/(sigma*a*math.sqrt(2*pi)))*np.exp(-(np.log(a) - miu)**2/(2*sigma**2)) for a in aafterol_predict2]
# 由对数正态分布拟合式计算beta

dkeqeff_predict2 = np.array([dkeq_predict2[seq]/(1+beta_predict2[seq]) for seq, _ in enumerate(dkeq_predict2)])
# 由beta的表达式计算dkeq_eff

dadn_predict2 = mts_analysis.WalkerCalculating(c0=c_w, m0=m_w, gamma=gamma, dk=dkeqeff_predict2, r=r)
# 由Walker模型得到dadn
############

plt.figure(num=1, figsize=(10, 8))
plt.scatter(a_Manual-a_ol, beta_w, label='$Beta R=0.1$', marker='2')
plt.plot(aplot, betafit, label='$beta fit R=0.1$', color='blue', linewidth=2)
plt.plot([plasticzone, plasticzone], [0, 1.5], label="$Plastic Zone Size R=0.1$", color="orange", linestyle="--")
plt.axis([0, 8, 0, 4])
plt.title("SIF Reduction Ratio(OLR=1.5)")
plt.ylabel("beta")
plt.xlabel("cracklength after overload/mm")
plt.legend()
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')

plt.figure(num=2, figsize=(10, 8))
plt.scatter(dk_Manual, dadn_Manual, label='$Experiment OLR=1.5$', marker='.')
#plt.scatter(dk_Manual2, dadn_Manual2, label='$Experiment OLR=2.0$', marker='.')
plt.scatter(dk_Manual, dadn_paris_dkeff, label='$CrackClosure OLR=1.5$', marker='*')
#plt.scatter(dk_Manual2, dadn_paris_dkeff2, label='$CrackClosure OLR=2.0$', marker='*')
plt.plot(dk_predict, dadn_predict, label='$FitbyBeta OLR=1.5$', color='black')
#plt.plot(dk_predict2, dadn_predict2, label='$FitbyBeta OLR=2.0$', color='black')
plt.axis([min(dk_Manual), max(dk_Manual), min(dadn_predict), max(dadn_Manual)*1.2])
plt.yticks(np.linspace(min(dadn_predict2)*0.9, max(dadn_Manual), 6))
plt.xticks(np.linspace(min(dk_Manual)*0.9, max(dk_Manual), 6))
plt.title("FCG Rates - Delta SIF")
#plt.title("CrackLength - Cycles")
plt.ylabel("FCG Rates/ mm per cycle")
#plt.ylabel("CrackLength / mm")
plt.xlabel("DeltaK/ MPa.m0.5")
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.show()
