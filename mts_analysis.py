# 处理数据结果函数

import numpy as np
import scipy.optimize as opt
import math

factor_for_k = 1e3*math.sqrt(1e-3)          # convert kN.mm**-2 to Mpa.m**-0.5

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


def WalkerFitting(dadn, dk, r=0.1):
    # usage: 通过离散的dadn, dk，r数据点拟合Walker模型，输出Paris公式参数C0，m0, gamma
    # input parameter:
    # dadn：裂纹扩展速率，单位 mm/cycles
    # dk：SIF变幅，单位 MPa.m0.5
    # r：应力比，无量纲，默认r = 0.1
    # return parameter:
    # Walker公式参数C0，单位 mm/cycles
    # Walker公式参数m0，无量纲
    # Walker公式参数gamma，无量纲
    #
    # 该函数尚未完工！！！
    #
    m0 = 3.39829
    lndadn = np.log(dadn)
    lndk = np.log(dk)
    def residual_walker(p):
        gamma_, logc0_ = p
        return (lndadn - m0 * lndk) - (logc0_ - (1 - gamma_) * m0 * np.log(1 - r))  # 应力比R = 0.1

    r = opt.leastsq(residual_walker, np.array([1e-7, 1.]))
    gamma, logc0 = r[0]
    c0 = np.exp(logc0)
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
    dadn = c0 * (dk * (1 - r) ** (gamma - 1)) ** m0
    return dadn


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
    p = p * 1e-3        # convert N to kN
    ux = 1/(np.sqrt(e*v*b/p)+1)
    alpha = ((((c5*ux + c4)*ux + c3)*ux + c2)*ux + c1)*ux + c0
    return alpha * w


def DeltaKCalculating(b, w, a, pmax, pmin):
    # usage: 依据E647-11标准计算CT试件的SIF变幅，适用范围a/W >= 0.2
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
    pmax = pmax * 1e-3      # convert N to kN
    pmin = pmin * 1e-3      # convert N to kN
    m1 = (pmax - pmin)/(b*np.sqrt(w))
    alpha = a/w
    m2 = (2 + alpha)/(1 - alpha)**1.5
    m3 = kc0 + kc1*alpha + kc2*alpha**2 + kc3*alpha**3 + kc4*alpha**4
    deltaK = m1*m2*m3*factor_for_k
    return deltaK


def FCGRateBySecant(a, n, dk, select=True):
    # usage: 依据E647-11标准采用割线法计算裂纹扩展速率（FCGR）dadn
    # Note: 由于Secant对数据进行了合并（两组计算一次）和筛选（略去FCGR<0的数据），
    # 因此输出将对应更新FCGR，循环次数，SIF变幅，裂纹长度数组
    # input parameter:
    # a: 裂纹长度（数组），单位 mm
    # n: 循环次数
    # dk：SIF变幅，单位 MPa.m0.5
    # return parameter:
    # deltaK：SIF变幅，单位 MPa.m0.5
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
    if (w - a) >= (4 / pi) * (kmax / (ys * 1000)**2):
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


def DataSelectByThreshold(threshold, parameter, data):
    # usage: 根据门槛值threshold进行筛选，
    # 若参考数组parameter的值小于等于阈值，则将排序对应的data数组中的数据保留，反之则丢弃
    # input parameter:
    # threshold: 阈值
    # parameter：参考数组，将该数组的值与阈值对比
    # data：被筛选的数组，通过数组序号进行筛选
    # return parameter：
    # newdata：筛选后的data数组
    newdata = []
    for seq, value in enumerate(parameter):
        if value <= threshold:
            newdata.append(data[seq])
    return newdata