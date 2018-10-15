# 裂纹闭合效应和柔度法计算相关函数库
# Last Update: 2018/10/15 Version 2.0.0
# Lu Yunchao

import math

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import paris_and_walker
from FCGAnalysisLib import read_data
from FCGAnalysisLib import mixed_model


def OutputOffset(cycles, load, cod, OpenCrackRange=[0.75, 1.00], ComplianceRange=0.1, ComplianceStep=0.5, plot=False):
    # usage: 根据详细的循环内轴力-COD（引伸计）位移数组拟合出指定位置的COD值
    # 方法：取出每个循环内[upperLSF, lowerLSF]之间的载荷和COD位移值，对其进行线性拟合
    # 然后输出比例rate的值在拟合线上的值返回
    # input parameter:
    # cycles: 循环数组
    # load：载荷数组，单位kN
    # cod：COD位移值数组，单位 mm，要求以上三个数组长度相同
    # upperLSF：根据标准E647-11（GB/T6398-2017），对于CT试件，柔度法上下存在非线性段。应舍去，upperLSF指上界的比例
    # lowerLSF：柔度法下界比例
    # rate1, 2：最终输出的载荷值的位置，eg.若输出rate=0.9，则输出的是Pmin+0.9*(Pmax-Pmin)的载荷及其对应的COD值
    # plot：是否绘制线性拟合图，用于检查upperLSF和lowerLSF是否合适，以及是否有其它逻辑错误
    # return parameter:
    # cycles_result: 输出循环次数
    # cod_result: 输出COD位移，单位 mm
    # loadrate1_result: 输出rate1比例对应的载荷值,单位 N
    # loadrate2_result: 输出rate2比例对应的载荷值，单位 N
    #
    # 放弃
    #
    n_list, n_seq = np.unique(cycles, return_index=True)
    n_couple_list = n_list.reshape(-1, 2)
    n_couple_seq = n_seq.reshape(-1, 2)

    cycles_result = []
    loadrate1_result = []
    loadrate2_result = []
    codrate1_result = []
    codrate2_result = []
    for seq in range(len(n_list)):
        start = n_seq[seq]
        if seq == len(n_list) - 1:
            end = len(cycles)
        else:
            end = n_seq[seq + 1]
        n_temp = cycles[start:end]
        load_temp = load[start:end]
        cod_temp = cod[start:end]
        # 提取出该循环内的所有数据
        load_lowerlimit = np.min(load_temp) + lowerLSF * np.ptp(load_temp)
        load_upperlimit = np.min(load_temp) + upperLSF * np.ptp(load_temp)
        load_target1 = np.min(load_temp) + rate1 * np.ptp(load_temp)
        load_target2 = np.min(load_temp) + rate2 * np.ptp(load_temp)
        # 设定对循环内数据的筛选标准及载荷峰值谷值
        n_selected = []
        load_selected = []
        cod_selected = []
        for seq2, _ in enumerate(load_temp):
            if (load_temp[seq2] > load_lowerlimit) and (load_temp[seq2] < load_upperlimit):
                n_selected.append(n_temp[seq2, 0])
                load_selected.append(load_temp[seq2, 0])
                cod_selected.append(cod_temp[seq2, 0])
        n_selected = np.array(n_selected)
        load_selected = np.array(load_selected)
        cod_selected = np.array(cod_selected)
        # 筛选留下载荷介于[load_lowerlimit, load_upperlimit]之间的所有数据
        p = np.polyfit(load_selected, cod_selected, deg=1)
        k, b = p
        cod_target1 = k * load_target1 + b
        cod_target2 = k * load_target2 + b
        # 将数据进行拟合并代入load_target输出目标COD值
        if plot:
            plt.figure(num=9, figsize=(10, 8))
            plt.scatter(cod_selected, load_selected, lw=1, marker='+', label='Data:' + str(n_list[seq]))
            load_fit = np.sort(load_selected)
            cod_fit = k * load_fit + b
            plt.plot(cod_fit, load_fit, linewidth=2, label='Fit' + str(n_list[seq]))
            plt.title("Cycles:" + str(n_selected[0]))
            plt.xlabel("COD/mm")
            plt.ylabel("LOAD/kN")
            plt.legend()
            plt.grid()

        if n_list[seq] != n_selected[0]:
            print("Something went WRONG in loop selecting! Check!")
            return 0, 0, 0
        # 程序逻辑检查，再次确认处理的循环数正确
        else:
            cycles_result.append(n_list[seq])
            loadrate1_result.append(load_target1 * 1e3)
            loadrate2_result.append(load_target2 * 1e3)
            codrate1_result.append(cod_target1)
            codrate2_result.append(cod_target2)
            # 将每个循环的计算结果添加到最终数组内
    return np.array(cycles_result), np.array(codrate1_result), np.array(codrate2_result), np.array(
        loadrate1_result), np.array(loadrate2_result)


def dKeffCalculation(kmax, kclose, kmin, method='basic'):
    # Usage: 采用不同的方法计算有效应力强度因子变幅Dk_eff.
    # 包括的方法有basic, 2/PI和2/PI0
    # basic: dkeff = kmax - kclosure
    # 2/PI: dkeff = kmax - 2/PI*kclosure - (1-2/PI)*kmin
    # 2/PI0: dkeff = kmax - 2/PI*kclosure
    # input parameter:
    # kmax: 应力强度因子峰值，单位MPa.m0.5
    # kclose: 裂纹闭合时应力强度因子值，单位MPa.m0.5
    # kmin: 应力强度因子谷值，单位MPa.m0.5
    # method: 计算的方法，字符，可填basic, 2/PI, 2/PI0
    # return parameter:
    # dkeff: 计算得到的有效应力强度因子变幅
    pi = math.pi
    if method == 'basic':
        dkeff = kmax - kclose
    elif method == '2/PI':
        dkeff = kmax - (2 / pi) * kclose - (1 - 2 / pi) * kmin
    elif method == '2/PI0':
        dkeff = kmax - (2 / pi) * kclose
    else:
        print("Error Method Parameter, Calculated by basic method.")
        dkeff = kmax - kclose
    return dkeff


def WeightedMovingAverage(data, numofaverage=7, weight=[]):
    # usage: 对数组data的数据进行带权移动平均，降低由于随机噪声带来的离散性
    # input parameter：
    # data：输入的待平均的数据
    # numofaverage：向前平均的数据个数，默认为7个数据
    # weight：计算移动平均的权，默认为等权
    # return parameter:
    # averageddata：完成移动平均的数据，数量较远数组减少(numofaverage-1)个
    if len(weight) == 0:
        weight = np.full(numofaverage, 1 / numofaverage)
    averageddata = np.zeros(len(data) - numofaverage + 1)
    for seq, _ in enumerate(data[(numofaverage - 1):]):
        data_to_average = data[seq: seq + numofaverage]
        averageddata[seq] = np.average(data_to_average, weights=weight, axis=0)

    return averageddata


def DoubleDataSelectMovingAverage(data, reference, ratio, numofaverage=7, weight=[]):
    # usage: 首先计算一次数据的移动平均，然后根据移动平均产生的误差带对数据进行筛选，保留误差带内的数据重新进行移动平均
    # Note：前(numofaverage-1)个数据由于不存在移动平均的值而无法进行筛选
    # input parameter：
    # data：进行移动平均和筛选的数据
    # reference：对数据进行标定位置的另一轴数据，该数据会根据列序号同步进行筛选，故需要保证data和reference一一对应
    # ratio:误差带半宽度，eg.0.2（保留移动平均下80%至120%的数据），该参数需要根据该数据的波动率进行调整，若ratio=0不筛选
    # numofaverage：移动平均采用的数据个数，默认为7
    # weight：移动平均的权，默认为直接平均（不带权）
    # return parameter:
    # data_f：输出的二次移动平均后的数据
    # raference_f：输出的二次移动平均后的参考数组

    # 第一次移动平均
    data_0 = WeightedMovingAverage(data=data, numofaverage=numofaverage, weight=weight)
    # reference_0 = np.array(reference[(numofaverage-1):])
    reference_0 = WeightedMovingAverage(data=reference, numofaverage=numofaverage, weight=weight)

    # 数据筛选
    if ratio > 0:
        reference_f = mts_analysis.DataSelectByRatio(ratio=ratio, points=data[(numofaverage - 1):],
                                                     average=data_0, data=reference_0)
        reference_f = np.array(reference_f)
        data_1 = mts_analysis.DataSelectByRatio(ratio=ratio, points=data[(numofaverage - 1):],
                                                average=data_0, data=data[(numofaverage - 1):])
    else:
        reference_f = reference_0
        data_1 = data_0

    # 恢复前(numofaverage-1)个数据
    data_2 = np.concatenate((data[0:numofaverage - 1], data_1))

    # 第二次移动平均
    data_f = WeightedMovingAverage(data=data_2, numofaverage=numofaverage, weight=weight)

    return data_f, reference_f


def BasicDataDispose(sequence, r, threshold=0):
    # usage: 常用的基础数据处理
    # 由基本的TestInf文件和原始实验数据OriginResult读取并计算，得到主要的裂纹扩展速率实验参数
    # input parameters:
    # sequence: 试件名（与MTS内设置一致）
    # r：实验应力比
    # threshold：数据筛选最小值，默认为0（即不筛选），可用于前段数据质量差时将其筛去
    # return paramater：
    # dadn、n、dk、a、kop：常见参数
    specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
        read_data.ReadTestInf(sequence=sequence)
    cycles, cracklength, kmax, kmin, pmax, pmin, kclosure, closureload = \
        read_data.ReadOriginResult(sequence=sequence, closure=True, cod=False)
    dadn_Manual, n_Manual, dk_Manual, a_Manual, kop_Manual = \
        experiment_calculation.FCGRandDKbyOriginalData(b=thickness, w=width, n=cycles, pmax=pmax, pmin=pmin,
                                                       a=cracklength,
                                                       ys=yield_strength, r=r, kclosure=kclosure,
                                                       threshold=threshold)
    return dadn_Manual, n_Manual, dk_Manual, a_Manual, kop_Manual


def Alpha2(sequence, stress_ratio, error_ratio, numofaverage, a_range, amplitude=0, fade_rate=0, fitting=0,
           model='paris'):
    # usage: 采用提出的Alpha方法先拟合出幅度系数M，再根据其计算dadn结果
    # input parameter:
    # sequence: 试件名（与MTS中一致，将直接读取文件）
    # stress_ratio：实验应力比
    # error_ratio：进行数据筛选的误差带半宽度，详细注释见函数DoubleDataSelectMovingAverage
    # numofaverage：进行移动平均的数据个数，详细注释见函数DoubleDataSelectMovingAverage
    # a_range：拟合的数据范围，格式为数组[a_begin, a_end]，a_begin推荐为dkeff_basic的谷值或dadn的谷值，
    # a_end要注意截止在数据开始离散之前
    # fade_rate：幂函数中的指数系数，直接指定，常见[5,8]
    # amplitude: 幂函数中的幅度系数，可定2/pi
    # 注：函数自动将对输入为0的项先拟合参数再进行计算，若均不为零则直接计算不进行拟合
    # fitting: 若fitting为1，则输出拟合时对标的dkeff_experiment值，以用于计算lost
    # model: 设定计算过程中采用的model，默认采用paris，可选walker
    # return parameter：
    # Ce1: 计算得到的幅度系数
    # dk_selected_averaged2, dkeff_alpha_e, dadn_alpha_e：对应的裂纹扩展参数
    # alpha_e：对应的参数alpha

    # 数据读入
    dadn_Manual, n_Manual, dk_Manual, a_Manual, kc_Manual = \
        BasicDataDispose(sequence=sequence, r=stress_ratio)
    if model == 'walker':
        c_w, m_w, gamma = paris_and_walker.WalkerParameter(data=1)
    else:
        c_p, m_p = paris_and_walker.ParisParameter(r=stress_ratio)

    # 移动平均
    kc_Manual_averaged2, dk_Manual_averaged2 = \
        DoubleDataSelectMovingAverage(data=kc_Manual, reference=dk_Manual,
                                      ratio=error_ratio, numofaverage=numofaverage)
    _, a_Manual_averaged2 = \
        DoubleDataSelectMovingAverage(data=kc_Manual, reference=a_Manual,
                                      ratio=error_ratio, numofaverage=numofaverage)
    _, dadn_Manual_averaged2 = \
        DoubleDataSelectMovingAverage(data=kc_Manual, reference=dadn_Manual,
                                      ratio=error_ratio, numofaverage=numofaverage)

    # 数据筛选a_range
    dadn_selecting_average2, dk_selecting_averaged2, kc_selecting_averaged2, a_selecting_averaged2 = \
        mts_analysis.FCGDataSelectByThreshold(dadn=dadn_Manual_averaged2, dk=dk_Manual_averaged2,
                                              n=kc_Manual_averaged2, a=a_Manual_averaged2,
                                              threshold=a_range[0], target='a', keepbigger=1)
    dadn_selected_averaged2, dk_selected_averaged2, kc_selected_averaged2, a_selected_averaged2 = \
        mts_analysis.FCGDataSelectByThreshold(dadn=dadn_selecting_average2, dk=dk_selecting_averaged2,
                                              n=kc_selecting_averaged2, a=a_selecting_averaged2,
                                              threshold=a_range[1], target='a', keepbigger=0)

    # 拟合数据准备
    kmin_selected_averaged2 = dk_selected_averaged2 * (stress_ratio / (1 - stress_ratio))  # Kmin计算
    kmax_selected_averaged2 = dk_selected_averaged2 * (1 / (1 - stress_ratio))  # Kmax计算
    if model == 'walker':
        dkeff_experiment_averaged2 = ((dadn_selected_averaged2 / c_w) ** (1 / m_w)) / (
            (1 - stress_ratio) ** (gamma - 1))
    else:
        dkeff_experiment_averaged2 = (dadn_selected_averaged2 / c_p) ** (1 / m_p)  # 拟合时的真值dkeff，由Paris给出
    beta = (a_selected_averaged2 - a_range[0]) / a_range[0]  # 裂纹长度无量纲量

    # 拟合，得到幅度或衰减系数
    if fade_rate != 0 and amplitude == 0:
        def residual_alpha_e(p):
            Ce = p
            return dkeff_experiment_averaged2 - (
                kmax_selected_averaged2 - (Ce * np.exp(-fade_rate * beta)) * kc_selected_averaged2 - (
                    1 - (Ce * np.exp(-fade_rate * beta))) * kmin_selected_averaged2)

        k = opt.leastsq(residual_alpha_e, np.array([1]))
        Ce1 = k[0]
    elif amplitude != 0 and fade_rate == 0:
        def residual_alpha_e(p):
            n = p
            return dkeff_experiment_averaged2 - (
                kmax_selected_averaged2 - (amplitude * np.exp(-n * beta)) * kc_selected_averaged2 - (
                    1 - (amplitude * np.exp(-n * beta))) * kmin_selected_averaged2)

        k = opt.leastsq(residual_alpha_e, np.array([1]))
        fade_rate = k[0]
        Ce1 = amplitude

    elif amplitude != 0 and fade_rate != 0:
        Ce1 = amplitude
    elif amplitude == 0 and fade_rate == 0:
        def residual_alpha_e(p):
            Ce, n = p
            return dkeff_experiment_averaged2 - (
                kmax_selected_averaged2 - (Ce * np.exp(-n * beta)) * kc_selected_averaged2 - (
                    1 - (Ce * np.exp(-n * beta))) * kmin_selected_averaged2)

        k = opt.leastsq(residual_alpha_e, np.array([1, 1]))
        Ce1, fade_rate = k[0]

    # 计算参数alpha
    alpha_e = Ce1 * np.exp(-fade_rate * beta)
    print("alpha幂函数拟合:试件名:", sequence, "幅度:", Ce1, "衰减速率:", fade_rate)

    # 计算dkeff_alpha和对应的dadn
    dkeff_alpha_e = kmax_selected_averaged2 - alpha_e * kc_selected_averaged2 - (1 - alpha_e) * kmin_selected_averaged2
    if model == 'walker':
        dadn_alpha_e = mts_analysis.WalkerCalculating(c0=c_w, gamma=gamma, m0=m_w, r=stress_ratio, dk=dkeff_alpha_e)
    else:
        dadn_alpha_e = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=dkeff_alpha_e)
    if fitting:
        return Ce1, fade_rate, dk_selected_averaged2, dkeff_alpha_e, dadn_alpha_e, alpha_e, dkeff_experiment_averaged2
    else:
        return Ce1, fade_rate, dk_selected_averaged2, dkeff_alpha_e, dadn_alpha_e, alpha_e


def AlphaCalculate(fade_rate, amplitude, dadn, dk, kc, a, a_range,
                   a_amplitude, stress_ratio, numofaverage=7, error_ratio=0.2):
    # 放弃
    # 移动平均
    kc_Manual_averaged2, dk_Manual_averaged2 = \
        DoubleDataSelectMovingAverage(data=kc, reference=dk,
                                      ratio=error_ratio, numofaverage=numofaverage)
    _, a_Manual_averaged2 = \
        DoubleDataSelectMovingAverage(data=kc, reference=a,
                                      ratio=error_ratio, numofaverage=numofaverage)
    _, dadn_Manual_averaged2 = \
        DoubleDataSelectMovingAverage(data=kc, reference=dadn,
                                      ratio=error_ratio, numofaverage=numofaverage)
    # 数据筛选a_range
    dadn_selecting_average2, dk_selecting_averaged2, kc_selecting_averaged2, a_selecting_averaged2 = \
        mts_analysis.FCGDataSelectByThreshold(dadn=dadn_Manual_averaged2, dk=dk_Manual_averaged2,
                                              n=kc_Manual_averaged2, a=a_Manual_averaged2,
                                              threshold=a_range[0], target='a', keepbigger=1)
    dadn_selected_averaged2, dk_selected_averaged2, kc_selected_averaged2, a_selected_averaged2 = \
        mts_analysis.FCGDataSelectByThreshold(dadn=dadn_selecting_average2, dk=dk_selecting_averaged2,
                                              n=kc_selecting_averaged2, a=a_selecting_averaged2,
                                              threshold=a_range[1], target='a', keepbigger=0)
    # 拟合数据准备
    kmin_selected_averaged2 = dk_selected_averaged2 * (stress_ratio / (1 - stress_ratio))  # Kmin计算
    kmax_selected_averaged2 = dk_selected_averaged2 * (1 / (1 - stress_ratio))  # Kmax计算
    c_p, m_p = paris_and_walker.ParisParameter(r=stress_ratio)                  # Paris参数读取
    betas = (a_selected_averaged2 - a_amplitude) / a_amplitude  # 裂纹长度无量纲量
    # 计算参数alpha
    alpha_e = []
    for beta in betas:
        if beta < 0:
            alpha_e.append(2/np.pi)
        else:
            alpha_e.append(amplitude * np.exp(-fade_rate * beta))
    alpha_e = np.array(alpha_e)
    # 基于参数计算有效应力强度因子
    dkeff_alpha_e = kmax_selected_averaged2 - alpha_e * kc_selected_averaged2 - (1 - alpha_e) * kmin_selected_averaged2
    # 计算dadn
    dadn_alpha_e = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=dkeff_alpha_e)
    return dk_selected_averaged2, dkeff_alpha_e, dadn_alpha_e, alpha_e, betas, a_selected_averaged2


def AlphaModifiedCalculate(alpha_fade_rate, alpha_amplitude, pmax, stress_ratio, thickness, width, a_amplitude, a_list,
                           c_ca, m_ca, a_start, a_end, dadn_start):
    # Modified Alpha模型计算函数
    # 该函数通过由已知起点(通常为a_kc_max_cal）和终点（通常为a_retard_end）的裂纹长度值和裂纹扩展速率
    # 采用直线拟合得到Kop的值，代入已知的Alpha模型参数值计算
    # 得到该段裂纹闭合效应主导的裂纹扩展速率模型
    # input parameters:
    # alpha_fade_rate: 模型衰减率参数
    # alpha_amplitude: 模型峰值参数，常取2/pi
    # pmax, stress_ratio: 载荷峰值，应力比
    # thickness, width: CT试件几何参数厚度和宽度
    # a_amplitude: 模型起点a0
    # a_list：要计算的裂纹长度值的数组
    # c_ca, m_ca：该应力比下的Paris公式参数
    # a_start, a_end: 要输出的裂纹长度起止点
    # dadn_start: 输入a0对应的dadn速率，通过反转Paris公式计算其dkeff
    # return parameters:
    # dk: 输入的a_list计算对应的dk
    # dkeff_alpha_e: 输入的a_list对应的alpha模型下的dkeff
    # dadn_alpha_e: 输入的a_list由alpha模型计算得到的dadn
    # alpha_e: a_list各值对应的alpha值
    # betas：a_list各值对应的无量纲裂纹长度参数值
    kmax = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=a_list, pmax=pmax, pmin=0)
    kmin = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=a_list, pmax=pmax*stress_ratio, pmin=0)
    dk = kmax - kmin
    # 计算无量纲量
    betas = (a_list - a_amplitude)/a_amplitude
    # 计算系数alpha
    alpha_e = []
    for beta in betas:
        if beta < 0:
            alpha_e.append(2/np.pi)
        else:
            alpha_e.append(alpha_amplitude * np.exp(-alpha_fade_rate * beta))
    alpha_e = np.array(alpha_e)
    # 计算Kop
    dk_start = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=a_start, pmax=pmax, pmin=pmax*stress_ratio)
    dk_end = mts_analysis.DeltaKCalculating(b=thickness, w=width, a=a_end, pmax=pmax, pmin=pmax*stress_ratio)
    dkeff_start = mts_analysis.ParisInverseCalculation(c=c_ca, m=m_ca, dadn=dadn_start)
    print("The Effective SIF when beta=0:" + str(dkeff_start))
    pi = np.pi
    kop_start = pi / 2 * (dk_start / (1 - stress_ratio) - dkeff_start - (1 - 2/pi) * dk_start * stress_ratio / (1 - stress_ratio))
    print("The Open SIF when beta=0 based on Alpha Method:" + str(kop_start))
    kop_end = dk_end * (stress_ratio / (1 - stress_ratio))
    dk_fit = [dk_start[0], dk_end]
    kop_fit = [kop_start[0], kop_end]
    k_kopfit_cal, b_kopfit_cal = mixed_model.LogLinearFit(y=kop_fit, x=dk_fit)
    kop_cal = mixed_model.LogLinearFittedCalculate(k=k_kopfit_cal, b=b_kopfit_cal, x=dk)

    # 计算dKeff
    dkeff_alpha_e = kmax - alpha_e * kop_cal - (1 - alpha_e) * kmin
    # 计算dadn
    dadn_alpha_e = mts_analysis.ParisCalculating(c=c_ca, m=m_ca, dk=dkeff_alpha_e)

    return dk, dkeff_alpha_e, dadn_alpha_e, alpha_e, betas
