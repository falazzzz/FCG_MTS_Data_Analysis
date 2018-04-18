from FCGAnalysisLib import mts_analysis
import matplotlib.pyplot as plt
from FCGAnalysisLib import read_data
from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import overload_analysis
import numpy as np
import math
from scipy import optimize

sequence = "yang-baoban_Lu-420-05"         # Graph Saving Sequence
a_ol = 12.3669
pol = 3000
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
dkeff_Manual = dk_Manual/(1 - r) - kc_Manual
# 读取结果，计算有效DK

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
betafit = [(k/(sigma*a*math.sqrt(2*pi)))*np.exp(-(np.log(a) - miu)**2/(2*sigma**2)) for a in aplot]
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

plt.figure(num=1, figsize=(10,8))
plt.scatter(a_Manual-a_ol, beta, label='$BetaByParis$', marker='1')
plt.scatter(a_Manual-a_ol, beta_w, label='$BetaByWalker$', marker='2')
plt.plot(aplot, betafit, label='$beta fit$', color='blue', linewidth=2)
plt.plot([plasticzone, plasticzone], [0, 1.5], label="$Plastic Zone Size$", color="orange", linestyle="--")
plt.axis([0, 8, 0, 4])
plt.title("SIF Reduction Ratio")
plt.ylabel("beta")
plt.xlabel("cracklength after overload/mm")
plt.legend()
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')

plt.figure(num=2, figsize=(10, 8))
plt.scatter(dk_Manual, dadn_Manual, label='$ExperimentData$', marker='+')
plt.scatter(dk_Manual, dadn_paris_dkeff, label='$dkeffResult$', marker='*')
plt.plot(dk_predict, dadn_predict, label='$Fitting$')
plt.axis([min(dk_Manual), max(dk_Manual), min(dadn_predict), max(dadn_Manual)*1.2])
plt.yticks(np.linspace(min(dadn_Manual)*0.9, max(dadn_Manual), 6))
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
