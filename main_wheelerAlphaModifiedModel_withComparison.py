# Residual stress - crack closure model
# 采用修正的Wheeler和修正的Alpha模型的程序实现，基于归纳的Kop峰值位置和归纳的dKeff线性变化规律，详见英文论文
# 带与原始Wheeler模型和2/pi模型的对比
# V1.1    2018/10/15

from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import paris_and_walker
from FCGAnalysisLib import overload_analysis
from FCGAnalysisLib import closureanalysis
from FCGAnalysisLib import mixed_model
from FCGAnalysisLib import write_data

import numpy as np
import matplotlib.pyplot as plt

'''
基本参数
'''
sequence = "yang-baoban_Lu-420-18"
show = 1  # 图像显示开关（没用）
moving_average = 0  # 移动平均开关
savedata = 0    # 保存数据

'''
模型参数
'''
bias = 0  # R = 0.1时调整alpha calculation的开关
eta = -0.5  # Modified-wheeler模型修正
kappa = 100
plz_factor = 'Xiaoping'


'''
初始化实验数据
'''
# 建立实例
specimen = mixed_model.OverloadSpecimen(name=sequence)
specimen.status_input_from_database()
specimen.basic_result_calculate()
# specimen.a_ol_applied -= 1.5
# specimen.a_dadn_min = 12.3506
# specimen.a_kc_max = specimen.a_alpha_end
specimen.a_alpha_end += 2
specimen.specimen_print()
# paris参数读取
c_ca, m_ca = paris_and_walker.ParisParameter(r=specimen.stress_ratio)

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
dadn_3, dk_3, kc_3, n_3, a_3 = specimen.data_select_by_cracklength(lowerlimit=specimen.a_dadn_min-0.2,
                                                                   upperlimit=specimen.a_retard_end+1.5)
dadn_4, dk_4, kc_4, n_4, a_4 = specimen.data_select_by_cracklength(lowerlimit=specimen.a_kc_max,
                                                                   upperlimit=specimen.a_retard_end)

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
dadn_a_kc_max_cal = mts_analysis.FindAscentDataBySeq(value=a_kc_max_calculate, item=a_3, target=dadn_3)
dk_a_kc_max_cal = mts_analysis.FindAscentDataBySeq(value=a_kc_max_calculate, item=a_3, target=dk_3)

# step4: 由dadn_a_kc_max_cal采用Paris逆推此时的dkeff
dkeff_a_kc_max_cal = mts_analysis.ParisInverseCalculation(c=c_ca, m=m_ca, dadn=dadn_a_kc_max_cal)
print("The Effective SIF when beta=0:"+str(dkeff_a_kc_max_cal))

# step5: 由dkeff_a_kc_max_cal依据Alpha模型计算Kop
pi = np.pi
kop_a_kc_max_cal = pi/2*(dk_a_kc_max_cal/(1-specimen.stress_ratio) - dkeff_a_kc_max_cal - (1-2/pi)*dk_a_kc_max_cal*specimen.stress_ratio/(1-specimen.stress_ratio))
print("The Open SIF when beta=0 based on Alpha Method:"+str(kop_a_kc_max_cal))

# step6: 拟合实验Kop得到斜率和截距，并得到计算版本的Kop值
k_kopfit_cal, b_kopfit_cal = mixed_model.LogLinearFit(y=kc_4, x=dk_4)
kc_cal = mixed_model.LogLinearFittedCalculate(k=k_kopfit_cal, b=b_kopfit_cal, x=dk_3)

# 校对: 用拟合结果计算2/Pi法(检验用，与模型无关)，与实验给出的Kop对应的计算结果进行对比
dkeff_2pi_cal = closureanalysis.dKeffCalculation(kmax=dk_3/(1-specimen.stress_ratio),
                                                 kclose=kc_cal,
                                                 kmin=dk_3*(specimen.stress_ratio/(1-specimen.stress_ratio)),
                                                 method='2/PI')
dadn_2pi_cal = mts_analysis.ParisCalculating(c=c_ca, m=m_ca, dk=dkeff_2pi_cal)

# step7: 用计算版的Kop和主导位置进行Alpha模型计算
dk_alpha_e, dkeff_alpha_e, dadn_alpha_e, alpha_e, betas, a_alpha_e = \
    closureanalysis.AlphaCalculate(fade_rate=7.2, amplitude=2 / np.pi, dadn=dadn_3, dk=dk_3,
                                   kc=kc_cal, a=a_3,
                                   a_range=[specimen.a_dadn_min - bias, specimen.a_alpha_end],
                                   a_amplitude=a_kc_max_calculate, stress_ratio=specimen.stress_ratio)


'''
残余应力部分：modified wheeler模型参数拟合和计算
'''
# 数据筛选（筛选出dadn_min-a_kc_max范围的数据，标记为1）
dadn_1, dk_1, kc_1, n_1, a_1 = specimen.data_select_by_cracklength(lowerlimit=specimen.a_dadn_min,
                                                                   upperlimit=specimen.a_kc_max)
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
m1 = overload_analysis.WheelerFittingBaseParis2(a=a_1, dadn=dadn_1, dk=dk_1, rm=plz_reverse_1,
                                                aol=specimen.a_ol_applied, rol=plz_ol, c=c_ca, m=m_ca, eta=eta,
                                                savedata=0)

# Modified Wheeler模型（引入参数eta）计算,塑性区系数由参数指定
dadn_wheeler, dk_wheeler, a_wheeler, cp_wheeler = \
    overload_analysis.WheelerCalculatingBaseParis(a_wheeler=a_alpha_e,
                                                  b=specimen.thickness, w=specimen.width, ys=specimen.yield_stress,
                                                  pmax=specimen.maxload, r=specimen.stress_ratio,
                                                  aol=specimen.a_ol_applied, rol=plz_ol, plasticzonefactor=plz_factor,
                                                  m1=m1, c=c_ca, m=m_ca, eta=eta)

'''
原始模型计算
'''
# Wheeler:
dadn_2, dk_2, kc_2, n_2, a_2 = specimen.data_select_by_cracklength(lowerlimit=specimen.a_dadn_min,
                                                                   upperlimit=specimen.a_alpha_end)
kmax_2 = np.array([dk / (1 - specimen.stress_ratio) for dk in dk_2])
plz_reverse_2 = overload_analysis.PlasticZoneWithFactor(kmax=kmax_2, ys=specimen.yield_stress,
                                                        factor=plz_factor, t=specimen.thickness)
m1_2 = overload_analysis.WheelerFittingBaseParis2(a=a_2, dadn=dadn_2, dk=dk_2, rm=plz_reverse_2,
                                                  aol=specimen.a_ol_applied, rol=plz_ol, c=c_ca, m=m_ca, eta=0,
                                                  savedata=0)
dadn_wheeler2, dk_wheeler2, a_wheeler2, cp_wheeler2 = \
    overload_analysis.WheelerCalculatingBaseParis(a_wheeler=np.arange(specimen.a_dadn_min, specimen.a_alpha_end, 0.02),
                                                  b=specimen.thickness, w=specimen.width, ys=specimen.yield_stress,
                                                  pmax=specimen.maxload, r=specimen.stress_ratio,
                                                  aol=specimen.a_ol_applied, rol=plz_ol, plasticzonefactor=plz_factor,
                                                  m1=m1_2, c=c_ca, m=m_ca, eta=0)
# 2/PI法:
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

'''
采用logistic sigmoid函数分配两种主导因素
'''
# logistic sigmoid单调递增，表征crack closure
sigmas = 1 / (1 + np.exp(-kappa * betas))
dadn_mixed_1 = dadn_wheeler * (1 - sigmas) + dadn_alpha_e * sigmas

'''
绘图部分
'''
# 图2: Kop拟合效果
plt.figure(num=2, figsize=(8, 6))
plt.scatter(dk_2, kc_2, label='$KcExperiment$', marker='2')
plt.scatter(dk_3, kc_cal, label='$KcCalculation$', marker='2')
plt.title("dKeff - deltaK(Wheeler Model),OLR=" + str(specimen.overload / specimen.maxload))
plt.ylabel("dKeff/MPa.m0.5")
plt.xlabel("DeltaSIF/MPa.m0.5")
plt.xscale('log')
plt.yscale('log')
#plt.axis([max(10, min(specimen.dk)), max(20, max(specimen.dk)), min(dkeff_2pi), max(20, max(specimen.dk))])
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()

# 图1：综合效果
plt.figure(num=1, figsize=(8, 6))
plt.scatter(specimen.dk, specimen.dadn, lw=1, marker='+', label='Experiment')
if moving_average:
    plt.scatter(dk_ave, dadn_ave, lw=1, marker='*', label='Moving average')
plt.plot(dk_wheeler2, dadn_wheeler2, label='Original Wheeler', color='black', linewidth=2)
plt.scatter(dk_2pi, dadn_2pi, label='$CrackClosure 2/PI Method$', marker='2')
#plt.scatter(dk_alpha_e, dadn_alpha_e, label='Modified Alpha Model')
#plt.scatter(dk_3, dadn_2pi_cal, label='Kop Fitting Result')
plt.scatter(dk_a_kc_max_cal, dadn_a_kc_max_cal, label='Point where Kop reach maximum')
plt.plot(dk_alpha_e, dadn_mixed_1, label='Modified Wheeler-Alpha Model', color='red', linewidth=2)
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
print(dk_alpha_e)
if savedata:
    data = np.array([specimen.dk, specimen.dadn, dk_alpha_e, dadn_mixed_1,
                     dk_wheeler2, dadn_wheeler2, dk_2pi, dadn_2pi])
    name = ['dk_experiment', 'dadn_experiment', 'dk_ModifiedAlpha', 'dadn_ModifiedAlpha',
            'dk_OriginalWheeler', 'dadn_OriginalWheeler', 'dk_2/PI', 'dadn_2/PI']
    _ = write_data.SaveData(dataset=data, name=name, filename=specimen.name+"_modified_wheeler_alpha")
