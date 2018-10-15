# Residual stress - crack closure model
# 带与原始模型对比，修正尝试不再需要详细的闭合数据，模型与Try3相同，但由figure读入
# 改为采用dk位置来标记分割
# V1.1.0    2018/8/6
# partly from main_wheeler
# 6082T6

from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import overload_analysis
from FCGAnalysisLib import closureanalysis
from FCGAnalysisLib import mixed_model
from FCGAnalysisLib import write_data

import numpy as np
import matplotlib.pyplot as plt

'''
基本参数
'''
sequence = "6082T6_data.csv"
show = 1  # 图像显示开关（没用）
moving_average = 0  # 移动平均开关
savedata = 1        # 数据保存开关

'''
模型参数
'''
bias = 0  # R = 0.1时调整wheeler calculation的开关
eta = 0  # Modified-wheeler模型修正
kappa = 100
plz_factor = 'Xiaoping'


'''
初始化实验数据
'''
# 建立实例
specimens = mixed_model.FigureDataReader(filename=sequence)
specimen = specimens[0]

# 材料参数
specimen.width = 60
specimen.thickness = 3
specimen.yield_stress = 0.323
specimen.elastic_modulus = 70

# 载荷参数
specimen.maxload = 2000
specimen.stress_ratio = 0.125
specimen.overload = 4000
specimen.a_ol_applied = 17.1
specimen.a_dadn_min = specimen.a_ol_applied + 0.5376
specimen.a_retard_end = specimen.a_dadn_min + 8
specimen.specimen_print()

# paris参数拟合
c_ca = 6.3307e-8
m_ca = 3.3236

'''
高载塑性区尺寸计算
'''
# Overload PLZ相关参数，塑性区系数由参数指定
kol = mts_analysis.DeltaKCalculating(b=specimen.thickness, w=specimen.width,
                                     a=specimen.a_ol_applied, pmax=specimen.overload, pmin=0)
plz_ol = overload_analysis.PlasticZoneWithFactor(kmax=np.array([kol]), ys=specimen.yield_stress,
                                                 factor=plz_factor, t=specimen.thickness)


'''
Try: 由塑性区比例特点计算Kop，代替实验实测值，进行Alpha模型计算
'''
# step1：数据筛选（拟合用4，计算用3，筛选出的数据上下限见参数lowerlimit和upperlimit）
dk_dadn_min = mts_analysis.DeltaKCalculating(b=specimen.thickness, w=specimen.width, a=specimen.a_dadn_min,
                                             pmax=specimen.maxload, pmin=specimen.maxload*specimen.stress_ratio)
dadn_3, dk_3 = specimen.data_select_by_sif(lowerlimit=dk_dadn_min, upperlimit=-1)

# step2：由比例采用二分法计算出主导位置交换的对应的裂纹长度a_kc_max的值
if specimen.overload != -1 and specimen.maxload != -1:
    ol_ratio = specimen.overload / specimen.maxload
    kmax_evaluate = mts_analysis.DeltaKCalculating(b=specimen.thickness, w=specimen.width,
                                                   a=plz_ol / 2 + specimen.a_ol_applied,
                                                   pmax=specimen.maxload, pmin=0)
    plz_evaluate = overload_analysis.PlasticZoneWithFactor(factor=plz_factor,
                                                           kmax=kmax_evaluate,
                                                           ys=specimen.yield_stress,
                                                           t=specimen.thickness)
    if ol_ratio == 1.5:
        olz_precentage = 0.80
    elif ol_ratio == 2:
        olz_precentage = 0.44
    else:
        print("ERROR! OLR is incorrect. Check again.")
        olz_precentage = -1
    a_kc_max_calculate, _ = mixed_model.BiSelection(a=0, b=plz_ol*olz_precentage, const=plz_ol*olz_precentage,
                                                    t=specimen.thickness, w=specimen.width, ys=specimen.yield_stress,
                                                    a_ol_applied=specimen.a_ol_applied, pmax=specimen.maxload,
                                                    plz_factor=plz_factor, accu=1e-4)
    a_kc_max_calculate += specimen.a_ol_applied
    print("OLR="+str(ol_ratio)+",exp.a_kc_max="+str(specimen.a_kc_max)+",cal.a_kc_max="+str(a_kc_max_calculate))

# step3：从a_1中提取计算得到的裂纹主导交换位置对应的扩展速率和SIF变幅
dk_a_kc_max_cal = mts_analysis.DeltaKCalculating(b=specimen.thickness, w=specimen.width, a=a_kc_max_calculate,
                                                 pmax=specimen.maxload, pmin=specimen.maxload*specimen.stress_ratio)
dadn_a_kc_max_cal = mts_analysis.FindAscentDataBySeq(value=dk_a_kc_max_cal, item=dk_3, target=dadn_3)

# step4: 由dadn_a_kc_max_cal采用Paris逆推此时的dkeff
dkeff_a_kc_max_cal = mts_analysis.ParisInverseCalculation(c=c_ca, m=m_ca, dadn=dadn_a_kc_max_cal)
print("The Effective SIF when beta=0:"+str(dkeff_a_kc_max_cal)+",The dSIF when beta=0:"+str(dk_a_kc_max_cal))

# step5: 由dkeff_a_kc_max_cal依据Alpha模型计算Kop
pi = np.pi
kop_a_kc_max_cal = pi/2*(dk_a_kc_max_cal/(1-specimen.stress_ratio) - dkeff_a_kc_max_cal - (1-2/pi)*dk_a_kc_max_cal*specimen.stress_ratio/(1-specimen.stress_ratio))
print("The Open SIF when beta=0 based on Alpha Method:"+str(kop_a_kc_max_cal))

# step6: 拟合实验Kop得到斜率和截距，并得到计算版本的Kop值
dk_a_retard_end = mts_analysis.DeltaKCalculating(b=specimen.thickness, w=specimen.width, a=specimen.a_retard_end,
                                                 pmax=specimen.maxload, pmin=specimen.maxload*specimen.stress_ratio)
kop_a_retard_end_cal = dk_a_retard_end * (specimen.stress_ratio/(1-specimen.stress_ratio))
dk_4 = np.array([dk_a_kc_max_cal[0], dk_a_retard_end])
kc_4 = np.array([kop_a_kc_max_cal[0], kop_a_retard_end_cal])
k_kopfit_cal, b_kopfit_cal = mixed_model.LogLinearFit(y=kc_4, x=dk_4)
kc_cal = mixed_model.LogLinearFittedCalculate(k=k_kopfit_cal, b=b_kopfit_cal, x=dk_3)

# 校对: 用拟合结果计算2/Pi法(检验用，与模型无关)，与实验给出的Kop对应的计算结果进行对比
dkeff_Manual_2pi_cal = closureanalysis.dKeffCalculation(kmax=dk_3/(1-specimen.stress_ratio),
                                                        kclose=kc_cal,
                                                        kmin=dk_3*(specimen.stress_ratio/(1-specimen.stress_ratio)),
                                                        method='2/PI')
dkeff_2pi_cal, dk_2pi_cal = \
    closureanalysis.DoubleDataSelectMovingAverage(data=dkeff_Manual_2pi_cal, reference=dk_3,
                                                  ratio=0.2, numofaverage=7)
dadn_2pi_cal = mts_analysis.ParisCalculating(c=c_ca, m=m_ca, dk=dkeff_2pi_cal)

# step7: 用计算版的Kop和主导位置进行Alpha模型计算
a_alpha_e = np.arange(specimen.a_dadn_min - bias, specimen.a_retard_end, 0.02)
dk_alpha_e, dkeff_alpha_e, dadn_alpha_e, alpha_e, betas =\
    closureanalysis.AlphaModifiedCalculate(alpha_fade_rate=7.2, alpha_amplitude=2/np.pi,
                                           pmax=specimen.maxload, stress_ratio=specimen.stress_ratio,
                                           thickness=specimen.thickness, width=specimen.width,
                                           a_amplitude=a_kc_max_calculate, a_list=a_alpha_e,
                                           c_ca=c_ca, m_ca=m_ca, a_start=a_kc_max_calculate,
                                           a_end=specimen.a_retard_end, dadn_start=dadn_a_kc_max_cal)

'''
残余应力部分：modified wheeler模型参数拟合和计算
'''
# 数据筛选（筛选出dadn_min-a_kc_max范围的数据，标记为1）
dadn_1, dk_1 = specimen.data_select_by_sif(lowerlimit=dk_dadn_min, upperlimit=dk_a_kc_max_cal)
# 移动平均（开关控制）
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

# 迟滞段塑性区计算
kmax_1 = np.array([dk / (1 - specimen.stress_ratio) for dk in dk_1])
plz_reverse_1 = overload_analysis.PlasticZoneWithFactor(kmax=kmax_1, ys=specimen.yield_stress,
                                                        factor=plz_factor, t=specimen.thickness)

# Wheeler模型指数系数m1拟合
a_1 = []
for sif in dk_1:
    a_1_temp, _ = mixed_model.BiCalculateCracklength(a=specimen.a_dadn_min, b=specimen.a_retard_end,
                                                     dk=sif, t=specimen.thickness, w=specimen.width,
                                                     pmax=specimen.maxload, r=specimen.stress_ratio)
    a_1.append(a_1_temp)
a_1 = np.array(a_1)
m1 = overload_analysis.WheelerFittingBaseParis2(a=a_1, dadn=dadn_1, dk=dk_1, rm=plz_reverse_1,
                                                aol=specimen.a_ol_applied, rol=plz_ol, c=c_ca, m=m_ca, eta=eta)

# Modified Wheeler模型（引入参数eta）计算,塑性区系数由参数指定
dadn_wheeler, dk_wheeler, a_wheeler, cp_wheeler = \
    overload_analysis.WheelerCalculatingBaseParis(a_wheeler=a_alpha_e,
                                                  b=specimen.thickness, w=specimen.width, ys=specimen.yield_stress,
                                                  pmax=specimen.maxload, r=specimen.stress_ratio,
                                                  aol=specimen.a_ol_applied, rol=plz_ol, plasticzonefactor=plz_factor,
                                                  m1=m1, c=c_ca, m=m_ca, eta=eta)

'''
采用logistic sigmoid函数分配两种主导因素
'''
# logistic sigmoid单调递增，表征crack closure
sigmas = 1 / (1 + np.exp(-kappa * betas))
dadn_mixed_1 = dadn_wheeler * (1 - sigmas) + dadn_alpha_e * sigmas

'''
绘图部分
'''

# 图1：综合效果
plt.figure(num=1, figsize=(8, 6))
plt.scatter(specimen.dk, specimen.dadn, lw=1, marker='+', label='Experiment')
plt.scatter(dk_2pi_cal, dadn_2pi_cal, label='Calculated Keff Value')
plt.scatter(dk_a_kc_max_cal, dadn_a_kc_max_cal, label='ResidualStress-ClosureEffect Convert')
plt.plot(dk_alpha_e, dadn_mixed_1, label='Modified Wheeler-Alpha Model', color='red', linewidth=2)
plt.title("FCG Rates - deltaK(Modified Wheeler-Alpha Model),Data:" + specimen.name)
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

if savedata:
    data = np.array([specimen.dk, specimen.dadn, dk_alpha_e, dadn_mixed_1])
    name = ['dk_experiment', 'dadn_experiment', 'dk_ModifiedAlpha', 'dadn_ModifiedAlpha']
    _ = write_data.SaveData(dataset=data, name=name, filename=sequence + "_modified_wheeler_alpha")