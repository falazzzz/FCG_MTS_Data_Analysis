# Residual stress - crack closure model
# 修正的Wheeler模型和原始的Alpha模型实现及与其它模型（Wheeler和2/pi法）对比
# V1.1    2018/10/15

from FCGAnalysisLib import read_data
from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import paris_and_walker
from FCGAnalysisLib import overload_analysis
from FCGAnalysisLib import closureanalysis

import numpy as np
import matplotlib.pyplot as plt


# 类定义
class OverloadSpecimen:
    # 由实验结果导数基本参数及实验数据
    # 需使用库read_data, experiment_calculation, mts_analysis

    def __init__(self, name, stress_ratio=1):
        self.name = name
        # Loading condition
        self.stress_ratio = stress_ratio
        self.maxload = -1
        # Geometry
        self.width = -1
        self.thickness = -1
        self.notch_length = -1
        # Material properties
        self.elastic_modulus = -1  # Gpa
        self.yield_stress = 0.365  # GPa
        self.possion_ratio = 0.3
        # Precrack information
        self.precrack_length = -1
        # Overload status
        self.overload = -1
        self.a_ol_applied = -1
        self.a_dadn_max = -1
        self.a_dadn_min = -1
        self.a_kc_max = -1
        # Alpha Fitting Region
        self.a_alpha_begin = -1
        self.a_alpha_end = -1
        # Experiment data
        self.dadn = []
        self.a = []
        self.n = []
        self.dk = []
        self.kc = []

    def basic_result_calculate(self):
        _, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
            read_data.ReadTestInf(sequence=self.name)
        cycles, cracklength, kmax, kmin, pmax, pmin, kc_Manual, pc_Manual = \
            read_data.ReadOriginResult(self.name, closure=True, cod=False)
        dadn_Manual, n_Manual, dk_Manual, a_Manual, kop_Manual = \
            experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin,
                                                           a=cracklength,
                                                           ys=yield_strength, r=self.stress_ratio, kclosure=kc_Manual,
                                                           threshold=0)
        self.width = width
        self.thickness = thickness
        self.notch_length = notch_length
        self.elastic_modulus = elastic_modulus
        self.precrack_length = precrack_length
        self.dadn = dadn_Manual
        self.n = n_Manual
        self.dk = dk_Manual
        self.a = a_Manual
        self.kc = kop_Manual

    def load_condition_input(self, maxload, stress_ratio=1):
        self.maxload = maxload
        if stress_ratio != 1:
            self.stress_ratio = stress_ratio
            print('Stress ratio changed to:' + str(self.stress_ratio))

    def overload_status_input(self, overload, a_ol_applid, a_dadn_min, a_dadn_max, a_kc_max):
        self.overload = overload
        self.a_ol_applied = a_ol_applid
        self.a_dadn_min = a_dadn_min
        self.a_dadn_max = a_dadn_max
        self.a_kc_max = a_kc_max

    def status_input_from_database(self):
        database = paris_and_walker.test_database2(specimen=self.name)
        self.maxload = database["maxload"]
        self.overload = database["overload"]
        self.stress_ratio = database["stress_ratio"]
        self.a_ol_applied = database["a_ol_applied"]
        self.a_dadn_max = database["a_dadn_max"]
        self.a_dadn_min = database["a_dadn_min"]
        self.a_kc_max = database["a_kc_max"]
        self.a_alpha_end = database["a_alpha_end"]

    def specimen_print(self):
        print('***********************************************')
        print('Specimen name:' + self.name)
        print('Geometry width:' + str(self.width) + 'mm, thickness:' + str(self.thickness) + 'mm')
        print('Material, E=' + str(self.elastic_modulus) + 'GPa,YS=' + str(self.yield_stress) + 'GPa')
        print('Loading condition: stress ratio=' + str(self.stress_ratio) + ',maxload=' + str(self.maxload) + 'N')
        print('Overload status:,overload=' + str(self.overload)
              + 'N, crack length WHEN OL applied:' + str(self.a_ol_applied)
              + 'mm / FCGRate reach maximum:' + str(self.a_dadn_max)
              + 'mm / FCGRate reach minimum:' + str(self.a_dadn_min)
              + 'mm / SIF at closure reach max:' + str(self.a_kc_max) + 'mm')
        print('************************************************')

    def data_select_by_cracklength(self, lowerlimit=-1, upperlimit=-1):
        a = self.a
        n = self.n
        dk = self.dk
        dadn = self.dadn
        kc = self.kc
        if lowerlimit != -1:
            dadn = mts_analysis.DataSelectByThreshold(threshold=lowerlimit, parameter=a, data=dadn, keepbigger=1)
            dk = mts_analysis.DataSelectByThreshold(threshold=lowerlimit, parameter=a, data=dk, keepbigger=1)
            n = mts_analysis.DataSelectByThreshold(threshold=lowerlimit, parameter=a, data=n, keepbigger=1)
            kc = mts_analysis.DataSelectByThreshold(threshold=lowerlimit, parameter=a, data=kc, keepbigger=1)
            a = mts_analysis.DataSelectByThreshold(threshold=lowerlimit, parameter=a, data=a, keepbigger=1)
        if upperlimit != -1:
            dadn = mts_analysis.DataSelectByThreshold(threshold=upperlimit, parameter=a, data=dadn, keepbigger=0)
            dk = mts_analysis.DataSelectByThreshold(threshold=upperlimit, parameter=a, data=dk, keepbigger=0)
            n = mts_analysis.DataSelectByThreshold(threshold=upperlimit, parameter=a, data=n, keepbigger=0)
            kc = mts_analysis.DataSelectByThreshold(threshold=upperlimit, parameter=a, data=kc, keepbigger=0)
            a = mts_analysis.DataSelectByThreshold(threshold=upperlimit, parameter=a, data=a, keepbigger=0)
        return np.array(dadn), np.array(dk), np.array(kc), np.array(n), np.array(a)


# 参数输入
sequence = "yang-baoban_Lu-420-19"
show = 1  # 图像显示开关（没用）
moving_average = 0  # 移动平均开关

# 模型参数
bias = 0  # R = 0.1时调整wheeler calculation的开关
eta = 2.5  # Modified-wheeler模型修正
kappa = 100
plz_factor = 'Xiaoping'

# 类建立
specimen = OverloadSpecimen(name=sequence)
specimen.status_input_from_database()
specimen.basic_result_calculate()
# specimen.a_ol_applied -= 1.5
# specimen.a_dadn_min = 13.5055
# specimen.a_kc_max = specimen.a_alpha_end
specimen.a_alpha_end += 2
specimen.specimen_print()

# paris参数
c_ca, m_ca = paris_and_walker.ParisParameter(r=specimen.stress_ratio)

# Overload PLZ相关参数，塑性区系数为Irwin
kol = mts_analysis.DeltaKCalculating(b=specimen.thickness, w=specimen.width,
                                     a=specimen.a_ol_applied, pmax=specimen.overload, pmin=0)
plz_ol = overload_analysis.PlasticZoneWithFactor(kmax=np.array([kol]), ys=specimen.yield_stress,
                                                 factor=plz_factor, t=specimen.thickness)

# 数据筛选（筛选出dadn_min后的数据，标记为1）
dadn_1, dk_1, kc_1, n_1, a_1 = specimen.data_select_by_cracklength(lowerlimit=specimen.a_dadn_min,
                                                                   upperlimit=specimen.a_kc_max)
# 移动平均
if moving_average:
    print('Notice: moving average working!')
    dadn_ave, dk_ave = closureanalysis.DoubleDataSelectMovingAverage(data=specimen.dadn, reference=specimen.dk,
                                                                     ratio=0.2, numofaverage=7)
    _, a_ave = closureanalysis.DoubleDataSelectMovingAverage(data=specimen.dadn, reference=specimen.a,
                                                             ratio=0.2, numofaverage=7)
    dadn_temp, dk_temp = closureanalysis.DoubleDataSelectMovingAverage(data=dadn_1, reference=dk_1,
                                                                       ratio=0.2, numofaverage=7)
    _, a_temp = closureanalysis.DoubleDataSelectMovingAverage(data=dadn_1, reference=a_1,
                                                              ratio=0.2, numofaverage=7)
    dadn_1 = dadn_temp[:]
    dk_1 = dk_temp[:]
    a_1 = a_temp[:]

# 高载段参数计算，塑性区系数为Irwin
kmax_1 = np.array([dk / (1 - specimen.stress_ratio) for dk in dk_1])
plz_reverse_1 = overload_analysis.PlasticZoneWithFactor(kmax=kmax_1, ys=specimen.yield_stress,
                                                        factor=plz_factor, t=specimen.thickness)

# Wheeler模型指数系数m1拟合
m1 = overload_analysis.WheelerFittingBaseParis2(a=a_1, dadn=dadn_1, dk=dk_1, rm=plz_reverse_1,
                                                aol=specimen.a_ol_applied, rol=plz_ol, c=c_ca, m=m_ca, eta=eta)

# Alpha模型部分计算
dk_alpha_e, dkeff_alpha_e, dadn_alpha_e, alpha_e, betas, a_alpha_e = \
    closureanalysis.AlphaCalculate(fade_rate=7.2, amplitude=2 / np.pi, dadn=specimen.dadn, dk=specimen.dk,
                                   kc=specimen.kc, a=specimen.a,
                                   a_range=[specimen.a_dadn_min - bias, specimen.a_alpha_end],
                                   a_amplitude=specimen.a_kc_max, stress_ratio=specimen.stress_ratio)

# Modified Wheeler模型计算,塑性区系数为Irwin
dadn_wheeler, dk_wheeler, a_wheeler, cp_wheeler = \
    overload_analysis.WheelerCalculatingBaseParis(a_wheeler=a_alpha_e,
                                                  b=specimen.thickness, w=specimen.width, ys=specimen.yield_stress,
                                                  pmax=specimen.maxload, r=specimen.stress_ratio,
                                                  aol=specimen.a_ol_applied, rol=plz_ol, plasticzonefactor=plz_factor,
                                                  m1=m1, c=c_ca, m=m_ca, eta=eta)
# 原始模型计算:
# Wheeler:
dadn_2, dk_2, kc_2, n_2, a_2 = specimen.data_select_by_cracklength(lowerlimit=specimen.a_dadn_min,
                                                                   upperlimit=specimen.a_alpha_end)
kmax_2 = np.array([dk / (1 - specimen.stress_ratio) for dk in dk_2])
plz_reverse_2 = overload_analysis.PlasticZoneWithFactor(kmax=kmax_2, ys=specimen.yield_stress,
                                                        factor=plz_factor, t=specimen.thickness)
m1_2 = overload_analysis.WheelerFittingBaseParis2(a=a_2, dadn=dadn_2, dk=dk_2, rm=plz_reverse_2,
                                                  aol=specimen.a_ol_applied, rol=plz_ol, c=c_ca, m=m_ca, eta=0)
dadn_wheeler2, dk_wheeler2, a_wheeler2, cp_wheeler2 = \
    overload_analysis.WheelerCalculatingBaseParis(a_wheeler=np.arange(specimen.a_dadn_min, specimen.a_alpha_end, 0.02),
                                                  b=specimen.thickness, w=specimen.width, ys=specimen.yield_stress,
                                                  pmax=specimen.maxload, r=specimen.stress_ratio,
                                                  aol=specimen.a_ol_applied, rol=plz_ol, plasticzonefactor=plz_factor,
                                                  m1=m1_2, c=c_ca, m=m_ca, eta=0)
# 2/PI法:
# dkeff计算
dkeff_Manual_2pi = closureanalysis.dKeffCalculation(kmax=specimen.dk/(1-specimen.stress_ratio),
                                                    kclose=specimen.kc,
                                                    kmin=specimen.dk*(specimen.stress_ratio/(1-specimen.stress_ratio)),
                                                    method='2/PI')
dkeff_2pi, dk_2pi = \
    closureanalysis.DoubleDataSelectMovingAverage(data=dkeff_Manual_2pi, reference=specimen.dk,
                                                  ratio=0.2, numofaverage=7)
_, a_2pi = \
    closureanalysis.DoubleDataSelectMovingAverage(data=dkeff_Manual_2pi, reference=specimen.a,
                                                  ratio=0.2, numofaverage=7)
dadn_2pi = mts_analysis.ParisCalculating(c=c_ca, m=m_ca, dk=dkeff_2pi)

# Try: mix try 1
# 采用logistic sigmoid函数分配两种
# logistic sigmoid单调递增，表征crack closure
sigmas = 1 / (1 + np.exp(-kappa * betas))
dadn_mixed_1 = dadn_wheeler * (1 - sigmas) + dadn_alpha_e * sigmas

# 绘图部分
# Plotting
plt.figure(num=1, figsize=(8, 6))
plt.scatter(specimen.dk, specimen.dadn, lw=1, marker='+', label='Experiment')
if moving_average:
    plt.scatter(dk_ave, dadn_ave, lw=1, marker='*', label='Moving average')
plt.plot(dk_wheeler2, dadn_wheeler2, label='Original Wheeler', color='black', linewidth=2)
plt.scatter(dk_2pi, dadn_2pi, label='$CrackClosure 2/PI Method$', marker='2')
plt.plot(dk_alpha_e, dadn_mixed_1, label='Mix Try 1', color='red', linewidth=2)
plt.title("FCG Rates - deltaK(Wheeler Model),OLR=" + str(specimen.overload / specimen.maxload))
plt.ylabel("FCG Rates/mm per cycle")
plt.xlabel("DeltaSIF/MPa.m0.5")
plt.xscale('log')
plt.yscale('log')
plt.axis([max(10, min(specimen.dk)), max(20, max(specimen.dk)), min(specimen.dadn) * 0.9, max(specimen.dadn) * 1.1])
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()
if show:
    plt.show()
