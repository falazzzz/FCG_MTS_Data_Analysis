import numpy as np
import matplotlib.pyplot as plt
import math


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
    # ！！！未完工！！！
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
    return np.array(cycles_result), np.array(codrate1_result), np.array(codrate2_result), np.array(loadrate1_result), np.array(loadrate2_result)


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
        dkeff = kmax - (2/pi)*kclose - (1-2/pi)*kmin
    elif method == '2/PI0':
        dkeff = kmax - (2/pi)*kclose
    else:
        print("Error Method Parameter, Calculated by basic method.")
        dkeff = kmax - kclose
    return dkeff
