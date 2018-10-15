# 处理数据结果函数库
# Last Update: 2018/4/11 Version 2.2.0
# Lu Yunchao

import numpy as np
import scipy.optimize as opt
import math
import matplotlib.pyplot as plt

factor_for_k = 1e3 * math.sqrt(1e-3)  # convert kN.mm**-2 to Mpa.m**-0.5


def ParisFitting(dadn, dk):
    # usage: 通过离散的dadn, dk数据点拟合Paris公式，将Paris公式对数化以进行线性拟合，输出Paris公式参数C，m
    # input parameter:
    # dadn：裂纹扩展速率，单位 mm/cycles
    # dk：SIF变幅，单位 MPa.m0.5
    # return parameter:
    # Paris公式参数C，单位 mm/cycles
    # Paris公式参数m，无量纲
    lndadn = np.log(dadn)
    lndk = np.log(dk)
    p = np.polyfit(lndk, lndadn, deg=1)
    m, b = p
    c = np.exp(b)
    return c, m


def ParisCalculating(c, m, dk):
    # usage: 给定SIF变幅dk，依据Paris公式计算裂纹扩展速率dadn
    # input parameter:
    # c: Paris公式参数C，单位推荐 mm/cycles
    # m: Paris公式参数m，无量纲
    # dk: 需要计算的对应SIF变幅dk，单位推荐MPa.m0.5
    # return parameter:
    # dadn: 裂纹扩展速率，按推荐输入对应单位为 mm/cycles
    dadn = c * (dk ** m)
    return dadn


def ParisInverseCalculation(c, m, dadn):
    # 完成依据Paris公式的逆运算
    # input arguments:
    # c, m: Paris参数
    # dadn：要计算的裂纹扩展速率
    # return argument:
    # dk：由Paris逆运算的计算结果
    dk = (dadn/c)**(1/m)
    return dk


def WalkerFittingFromParis(dadn1, dk1, r1, dadn2, dk2, r2, m):
    # usage: 通过不同应力比下的Paris参数拟合Walker模型，输出Walker模型参数C0，m0, gamma
    # NOTE: 该方法参考论文 Modelling the fatigue crack growth in friction stir welded joint of 2024-T351 Al alloy
    # m0取不同应力比下paris参数m的平均值
    # c0和gamma则由不同应力比R下的Paris参数C关系（应变量R-自变量C）的线性拟合得到
    # input parameter:
    # c: 数组，不同应力比下对应的Paris参数C
    # m：数组，不同应力比下对应的Paris参数m
    # r：应力比
    # return parameter:
    # Walker公式参数C0，单位 mm/cycles
    # Walker公式参数m0，无量纲
    # Walker公式参数gamma，无量纲
    m = np.array(m)
    m0 = np.average(m)
    Y = np.concatenate((np.log(dadn1), np.log(dadn2))) - m0 * np.concatenate((np.log(dk1), np.log(dk2)))
    r_1 = np.full(len(dadn1), r1)
    r_2 = np.full(len(dadn2), r2)
    X = np.log(1 - np.concatenate((r_1, r_2)))

    def residual_walker(p):
        K, B = p
        return Y - (K * X + B)

    r = opt.leastsq(residual_walker, np.array([1, 1]))
    k, b = r[0]
    c0 = np.exp(b)
    gamma = k/m0 + 1
    return c0, m0, gamma


def WalkerFittingByRegression(dadn1, dk1, r1, dadn2, dk2, r2, dadn3=[], dk3=[], r3=0):
    # usage: 通过不同应力比实验得到的dadn和dk数据拟合Walker模型，输出Walker模型参数C0，m0, gamma
    # NOTE: 目前只能拟合两个应力比的数据，值对数化后直接进行三元平面拟合
    #      log(da/dN) = m*log(dk) + m(gamma-1)*log(1-r) + logc0
    # =>        Z     = A*   X    +      B    *    Y    +   C
    # input parameter:
    # dadn1,dadn2：应力比r1和r2得到的裂纹扩展速率，单位 mm/cycles
    # dk1,dadn2：应力比r1和r2得到的SIF变幅，单位 MPa.m0.5
    # r1,r2：应力比，无量纲，默认r = 0.1
    # return parameter:
    # Walker公式参数C0，单位 mm/cycles
    # Walker公式参数m0，无量纲
    # Walker公式参数gamma，无量纲
    #
    # 函数通过检验
    #
    if len(dadn3) == 0:
        Z = np.concatenate((np.log(dadn1), np.log(dadn2)))
        X = np.concatenate((np.log(dk1), np.log(dk2)))
        r_1 = np.full(len(dadn1), r1)
        r_2 = np.full(len(dadn2), r2)
        Y = np.log(1 - np.concatenate((r_1, r_2)))
    else:
        Z = np.concatenate((np.log(dadn1), np.log(dadn2), np.log(dadn3)))
        X = np.concatenate((np.log(dk1), np.log(dk2), np.log(dk3)))
        r_1 = np.full(len(dadn1), r1)
        r_2 = np.full(len(dadn2), r2)
        r_3 = np.full(len(dadn3), r3)
        Y = np.log(1 - np.concatenate((r_1, r_2, r_3)))

    def residual_walker(p):
        A, B, C = p
        return np.abs(A*X + B*Y - Z + C)/math.sqrt(A**2 + B**2 + 1)

    r = opt.leastsq(residual_walker, np.array([4., 0., -20.]))
    a, b, c = r[0]
    m0 = a
    gamma = (b + m0) / m0
    c0 = np.exp(c)
    return c0, m0, gamma


def WalkerCalculating(c0, m0, gamma, dk, r):
    # usage: 给定SIF变幅dk和应力比r，依据Walker模型计算裂纹扩展速率dadn
    # input parameter:
    # c0: Walker模型参数C0，单位推荐 mm/cycles
    # m0: Walker模型参数m0，无量纲
    # gamma：Walker模型参数gamma，无量纲
    # dk: 需要计算的对应SIF变幅dk，单位推荐MPa.m0.5
    # r：需要计算的dk数据对应的应力比r，无量纲
    # return parameter:
    # dadn: 裂纹扩展速率，按推荐输入对应单位为 mm/cycles
    dadn = c0 * ((dk * ((1 - r) ** (gamma - 1))) ** m0)
    return dadn


def OutputCodandLoad(cycles, load, cod, upperLSF=0.9, lowerLSF=0.1, rate1=0.9, rate2=0.1, plot=False):
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
    n_list, n_seq = np.unique(cycles, return_index=True)
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


def Compliance(e, b, w, p, v):
    # usage: 依据E647-11采用柔度法计算裂纹长度，有效范围0.2 <= a/W <= 0.975
    # input parameter:
    # e: 弹性模量，单位 GPa
    # b: 试件厚度thickness，单位 mm
    # w: 试件宽度Width，单位 mm
    # p：载荷p，单位 N
    # v: cod读数，单位 mm
    # return parameter:
    # a = alpha*w: 裂纹长度，单位 mm
    c0 = 1.0010
    c1 = -4.6695
    c2 = 18.460
    c3 = -236.82
    c4 = 1214.9
    c5 = -2143.6
    p = p * 1e-3  # convert N to kN
    ux = 1 / (np.sqrt(e * v * b / p) + 1)
    alpha = ((((c5 * ux + c4) * ux + c3) * ux + c2) * ux + c1) * ux + c0
    return alpha * w


def DeltaKCalculating(b, w, a, pmax, pmin):
    # usage: 依据E647-11标准计算CT试件的SIF变幅，适用范围a/W >= 0.2且韧带长度有效
    # input parameter:
    # b: 试件厚度thickness，单位 mm
    # w: 试件宽度Width，单位 mm
    # a：对应裂纹长度，单位 mm
    # pmax：对应载荷最大值，单位 N
    # pmin：对应载荷最小值，单位 N
    # return parameter:
    # deltaK：SIF变幅，单位 MPa.m0.5
    kc0 = 0.886
    kc1 = 4.64
    kc2 = -13.32
    kc3 = 14.72
    kc4 = -5.6
    pmax = pmax * 1e-3  # convert N to kN
    pmin = pmin * 1e-3  # convert N to kN
    m1 = (pmax - pmin) / (b * np.sqrt(w))
    alpha = a / w
    m2 = (2 + alpha) / (1 - alpha) ** 1.5
    m3 = kc0 + kc1 * alpha + kc2 * alpha ** 2 + kc3 * alpha ** 3 + kc4 * alpha ** 4
    deltaK = m1 * m2 * m3 * factor_for_k
    return deltaK


def DeltaKCalculatingByMurakami(b, w, a, pmax, pmin):
    # usage: 依据Murakami在1987年给出的公式计算CT试件的SIF变幅，适用范围a/W >= 0.2且韧带长度有效
    # input parameter:
    # b: 试件厚度thickness，单位 mm
    # w: 试件宽度Width，单位 mm
    # a：对应裂纹长度，单位 mm
    # pmax：对应载荷最大值，单位 N
    # pmin：对应载荷最小值，单位 N
    # return parameter:
    # deltaK：SIF变幅，单位 MPa.m0.5
    alphas = a/w
    p = pmax - pmin
    s1 = 16.7
    s2 = -104.7
    s3 = 369.9
    s4 = -573.8
    s5 = 360.5
    k1 = []
    for seq, alpha in enumerate(alphas):
        k1l = (p[seq]/b)*math.sqrt(math.pi/w)
        k1r = (s1*alpha**0.5 + s2*alpha**1.5 + s3*alpha**2.5 + s4*alpha**3.5 + s5*alpha**4.5)
        k1.append(k1l * k1r / math.sqrt(1e3))
    return np.array(k1)


def FCGRateBySecant(a, n, dk, select=True):
    # usage: 依据E647-11标准采用割线法计算裂纹扩展速率（FCGR）dadn
    # Note: 由于Secant对数据进行了合并（两组计算一次）和筛选（略去FCGR<0的数据），
    # 因此输出将对应更新FCGR，循环次数，SIF变幅，裂纹长度数组
    # input parameter:
    # a: 裂纹长度（数组），单位 mm
    # n: 循环次数
    # dk：SIF变幅，单位 MPa.m0.5
    # return parameter:
    # 更新的dadn, cracklength, cycles, dk
    dadn_secant = []
    cracklength_secant = []
    cycles_secant = []
    dk_secant = []
    i = 0

    def secant(a_, n_, dk_, i_):
        # 割线法函数
        dadn_ = (a_[i_ + 1] - a_[i_]) / (n_[i_ + 1] - n_)
        return dadn_[i_], a_[i_], n_[i_], dk_[i_]

    while i < len(n) - 1:
        result = secant(a, n, dk, i)
        i = i + 1
        if select:
            if result[0] < 0:
                continue
        dadn_secant.append(result[0])
        cracklength_secant.append(result[1])
        cycles_secant.append(result[2])
        dk_secant.append(result[3])
    return dadn_secant, cracklength_secant, cycles_secant, dk_secant


def FCGRateByPolynomial(a, n, dk, select=True, num=3):
    # usage: 依据E647-11标准采用7点多项式拟合法计算裂纹扩展速率（FCGR）dadn
    # Note: 由于Polynomial对数据进行了合并（7组计算一次）和筛选（略去FCGR<0的数据），
    # 因此输出将对应更新FCGR，循环次数，SIF变幅，裂纹长度数组（缩短前后共2*num个数据点）
    # input parameter:
    # a: 裂纹长度（数组），单位 mm
    # n: 循环次数
    # dk：SIF变幅，单位 MPa.m0.5
    # select：是否筛除计算得到dadn < 0的值，默认将筛除
    # num：多项式拟合时选取的数据点数为2num+1，若取7点，则num=3，若取9点，则num=4，默认3
    # return parameter:
    # 更新的dadn, cracklength, cycles, dk
    seq = num
    dadn_poly = []
    a_poly = []
    n_poly = []
    dk_poly = []
    while seq < len(a) - num:
        c1 = 0.5 * (n[seq - num] + n[seq + num])
        c2 = 0.5 * (n[seq + num] - n[seq - num])
        a_temp = a[seq - num:seq + num]
        n_temp = n[seq - num:seq + num]
        x = (n_temp - c1) / c2
        p = np.polyfit(x, a_temp, deg=2)
        b2, b1, _ = p
        dadn_temp = b1 / c2 + 2 * b2 * (n[seq] - c1) / c2 ** 2
        dadn_poly.append(dadn_temp)
        seq = seq + 1
        if select:
            if dadn_temp < 0:
                continue
        a_poly.append(a[seq])
        n_poly.append(n[seq])
        dk_poly.append(dk[seq])
    return dadn_poly, a_poly, n_poly, dk_poly


def LigamentValidCheck(w, a, dk, ys, r=0.1):
    # usage: 依据E647-11标准建议韧带长度是否符合要求
    # input parameter:
    # w: 试件宽度Width，单位 mm
    # a：对应裂纹长度，单位 mm
    # dk：SIF变幅，单位 MPa.m0.5
    # ys：屈服强度yield_strength，单位 GPa
    # r：应力比
    # return parameter:
    # 布尔值，若符合韧带条件则返回True，反正返回False
    pi = 3.1415926
    kmax = dk / (1 - r)
    if (w - a) * 1e-3 >= (4. / pi) * (kmax / (ys * 1000)) ** 2:
        return True
    else:
        return False


def DataSelectByLigament(w, a, dk, ys, data, r=0.1):
    # usage: 根据LigamentValidCheck函数返回的布尔值结果对数据进行筛选，True则保留数据，False则丢弃
    # input parameter:
    # w: 试件宽度Width，单位 mm
    # a：对应裂纹长度，单位 mm
    # dk：SIF变幅，单位 MPa.m0.5
    # ys：屈服强度yield_strength，单位 GPa
    # data：需要被筛选的数据（数组）
    # r：应力比
    # return parameter:
    # newdata：经过筛选的data数组
    seq = 0
    newdata = []
    while seq < len(a):
        judge = LigamentValidCheck(w=w, a=a[seq], dk=dk[seq], ys=ys, r=r)
        if judge:
            newdata.append(data[seq])
        seq = seq + 1
    return newdata


def DataSelectByThreshold(threshold, parameter, data, keepbigger=1):
    # usage: 根据门槛值threshold进行筛选，
    # 若参考数组parameter的值大于等于阈值，则将排序对应的data数组中的数据保留，反之则丢弃
    # input parameter:
    # threshold: 阈值
    # parameter：参考数组，将该数组的值与阈值对比
    # data：被筛选的数组，通过数组序号进行筛选
    # return parameter：
    # newdata：筛选后的data数组
    newdata = []
    for seq, value in enumerate(parameter):
        if keepbigger == 0:
            if value < threshold:
                newdata.append(data[seq])
        elif value >= threshold:
            newdata.append(data[seq])
    return newdata


def DataSelectByRatio(ratio, points, average, data):
    # usage: 根据比例ratio进行筛选，
    # 若参考数组parameter的值在半宽度为ratio误差带内，则将排序对应的data数组中的数据保留，反之则丢弃
    # input parameter:
    # ratio: 误差带半宽度
    # points：被筛选的数据（波动性大的点）
    # average：被筛选的数据的判断标准（即误差带的中轴）
    # data：被筛选的数组，通过数组序号进行筛选
    # return parameter：
    # newdata：筛选后的data数组
    newdata = []
    for seq, value in enumerate(points):
        if (1 - ratio)*value < average[seq] < (1 + ratio)*value:
            newdata.append(data[seq])
    return newdata


def FCGDataSelectByThreshold(dadn, dk, threshold, a=[], n=[], target='dk', keepbigger=1):
    # usage: 一组筛选FCG的四个主要参数dadn、dk、n和a，保留阈值以上/下的值
    # input parameter:
    # dadn, dk, n, a: 筛选前四个主要参数的原始数据，其中a、n可不输入，但不输入a、n时只能用dk为筛选参量
    # threshold: 筛选的阈值
    # target：筛选的参数，可选:'dk', 'a', 'n'
    # keepbigger：1：保留大于阈值的部分，0：保留小于阈值的部分
    # return parameter:
    # dadn_new, dk_new, n_new, a_new：筛选保留下的新数据
    if target == 'dk':
        parameter = dk
    elif target == 'a':
        parameter = a
        if len(a) == 0:
            print('FCGDataSelectByThreshold Function Error: No cracklength data input while selecting by cracklength.')
    elif target == 'n':
        parameter = n
        if len(n) == 0:
            print('FCGDataSelectByThreshold Function Error: No cycles data input while selecting by cycles.')
    else:
        print('Error parameter: target, set dk as default.')
        parameter = 'dk'
    dadn_new = DataSelectByThreshold(threshold=threshold, parameter=parameter,
                                     data=dadn, keepbigger=keepbigger)
    dk_new = DataSelectByThreshold(threshold=threshold, parameter=parameter,
                                   data=dk, keepbigger=keepbigger)
    if len(n) != 0:
        n_new = DataSelectByThreshold(threshold=threshold, parameter=parameter,
                                      data=n, keepbigger=keepbigger)
    if len(a) != 0:
        a_new = DataSelectByThreshold(threshold=threshold, parameter=parameter,
                                      data=a, keepbigger=keepbigger)
    if len(n) == 0 and len(a) == 0:
        return np.array(dadn_new), np.array(dk_new)
    elif len(n) != 0 and len(a) == 0:
        return np.array(dadn_new), np.array(dk_new), np.array(n_new)
    elif len(n) == 0 and len(a) != 0:
        return np.array(dadn_new), np.array(dk_new), np.array(a_new)
    else:
        return np.array(dadn_new), np.array(dk_new), np.array(n_new), np.array(a_new)


def FindAscentDataBySeq(value, item, target):
    # usage: 在item数组中找出与value值最接近的值，在target数组中找到对应的值返回
    # 注意：要求item数组和target数组的排序是相同且是升序的，即行行数据对应
    # 找出最近的两个值，进行线性插值返回
    # input parameter:
    # value: 要寻找的值
    # item: 数组，函数将在item中寻找与value相近的值
    # target: 数组，函数将返回item中与value相近值的序号相同的值并进行插值
    # return parameter:
    # result: 搜寻并拟合的值
    item = np.array(item)
    seq = 0
    if value > item[-1] or value < item[0]:
        result = 0
        print('FunctionFindAscentDataBySeqERROR: Value cannot be found in the list.')
        return result
    while seq < len(item):
        if item[seq] >= value:
            ratio = (value - item[seq - 1]) / (item[seq] - item[seq - 1])
            result = target[seq - 1] + ratio * (target[seq] - target[seq - 1])
            break
        seq = seq + 1
    return result


def CTS_Richard_K1_calculate(f, r, a, w, b, alpha):
    # 根据Richard提出的公式，计算对应裂纹长度下CTS试件的K1值
    # Testing 20180705
    # 输入单位：力N，长度mm，角度为弧度制
    # 适用范围: 0.5<=a/W<=0.7
    df = f * (1 - r)
    x1 = (df / w / b) * np.sqrt(np.pi * a)
    x2 = np.cos(alpha) / (1 - a / w)
    factor = a / (w - a)
    x3 = np.sqrt((0.26 + 2.65 * factor)) / (1 + 0.55*factor - 0.08*factor**2)
    ctsK1 = x1*x2*x3/factor_for_k
    return ctsK1


def CTS_Richard_K2_calculate(f, r, a, w, b, alpha):
    # 根据Richard提出的公式，计算对应裂纹长度下CTS试件的K2值
    # Testing 20180705
    # 输入单位：力N，长度mm，角度为弧度制
    # 适用范围: 0.5<=a/W<=0.7
    df = f * (1 - r)
    x1 = (df / w / b) * np.sqrt(np.pi * a)
    x2 = np.sin(alpha) / (1 - a / w)
    factor = a / (w - a)
    x3 = np.sqrt((-0.23 + 1.40 * factor)) / (1 - 0.67*factor + 2.08*factor**2)
    ctsK2 = x1*x2*x3/factor_for_k
    return ctsK2


def CTS_Soon_K1_calculate(f, r, a, w, b, alpha):
    # 根据Soon提出的公式，计算对应裂纹长度下CTS试件的K1值
    # Testing 20180706
    # 输入单位：力N，长度mm，角度为弧度制
    # 适用范围: 0.3<=a/W<=0.5
    df = f * (1 - r)
    x1 = df * np.sqrt(np.pi*a) / (w*b)
    x2 = np.cos(alpha) * np.sqrt(np.cos(alpha/3))
    f1 = 2.32158 - 14.36777*a/w + 66.85752*(a/w)**2 - 117.66921*(a/w)**3 + 89.72502*(a/w)**4
    ctsK1 = x1*x2*f1/factor_for_k
    return ctsK1


def CTS_Soon_K2_calculate(f, r, a, w, b, alpha):
    # 根据Soon提出的公式，计算对应裂纹长度下CTS试件的K2值
    # Testing 20180706
    # 输入单位：力N，长度mm，角度为弧度制
    # 适用范围: 0.3<=a/W<=0.5
    df = f * (1 - r)
    x1 = df * np.sqrt(np.pi * a) / (w * b)
    x2 = np.sin(alpha)
    f2 = -0.05741 + 4.36076*(a/w) - 4.46168*(a/w)**2 + 2.48807*(a/w)**3
    ctsK2 = x1 * x2 * f2 / factor_for_k
    return ctsK2

