# read_data部分，负责读取文件并整理格式，返回数组
# Last Update: 2018/4/13 Version 2.2.2
# Lu Yunchao

import math
import pandas as pd
import numpy as np
import sys

factor_for_k = 1e3 * math.sqrt(1e-3)  # convert kN.mm**-2 to Mpa.m**-0.5
inputfile = "\\inputfile\\"             # 读取文件放置的子文件夹名
path = sys.path[0] + inputfile  # 设定读取的绝对路径


def ReadMtsResult(sequence, dataselect=1, dkeffread=0):
    # usage: 读取MTS输出的拟合计算结果，共有裂纹扩展速率da/dN，循环数cycles，SIF变幅dk三项
    # 若dataselect=1，则会根据数据Validity一栏对结果的有效性进行筛选，删去扩展速率为负的结果
    # input parameter:
    # sequence: 文件名前缀，通常由batch名和specimen名组成
    # dataselect: 是否根据文件内的标签对数据的有效性进行筛选
    # dkeffread :是否读取等效SIF变幅数据，默认不读取
    # return parameter:
    # 裂纹扩展速率dadn，循环数cycles，SIF变幅dk，等效SIF变幅dkeff
    mtsresult = pd.read_csv(path + sequence + u"_da dN, Delta K by Cycles.csv")  # Data File Reading
    dadn = mtsresult["da/dN (mm/cycle)"]
    cycles = mtsresult["Cycles"]
    dk = mtsresult["Delta K Applied (kN/mm^1.5)"] * factor_for_k
    dkeff = mtsresult["Delta K Effective (kN/mm^1.5)"] * factor_for_k
    if dataselect:
        dadn_s = []
        cycles_s = []
        dk_s = []
        if dkeffread == 1:
            dkeff_s = []
        validity = mtsresult["Validity"]
        for seq, value in enumerate(validity):
            if value == "Valid":
                dadn_s.append(dadn[seq])
                cycles_s.append(cycles[seq])
                dk_s.append(dk[seq])
                if dkeffread == 1:
                    dkeff_s.append(dkeff[seq])
        dadn_s = np.array(dadn_s)
        cycles_s = np.array(cycles_s)
        dk_s = np.array(dk_s)
        if dkeffread == 1:
            dkeff_s = np.array(dkeff_s)
    else:
        dadn_s = dadn
        cycles_s = cycles
        dk_s = dk
        if dkeffread == 1:
            dkeff_s = dkeff
    if dkeffread == 1:
        return dadn_s, cycles_s, dk_s, dkeff_s
    else:
        return dadn_s, cycles_s, dk_s


def ReadOriginResult(sequence, closure=False, cod=True):
    # usage: 读取MTS输出的原始数据，共有循环次数，裂纹长度，SIF最大最小值，K Closure， Closure Load，载荷最大最小值，
    # COD最大最小值共10项，其中cod参数默认输出，closure默认不输出，需要注意只有在MTS实验机设置为store timed data时才
    # 会有closure部分的结果
    # input parameter:
    # sequence: 文件名前缀，通常由batch名和specimen名组成
    # closure: 控制输出项K closure和Closure Load，默认为不输出
    # cod：控制引伸计输出项codmin和codmax，默认输出
    # return parameter:
    # 循环次数cycle，裂纹长度cracklength，SIF最大值kmax，SIF最小值kmin，载荷最大值pmax，载荷最小值pmin
    # (若closure打开）裂纹闭合效应参数kclosure，closureload，（若cod打开）张开位移最大值codmax，张开位移最小值codmin
    mtsdata = pd.read_csv(path + sequence + u"_Crack Length, Min Max K, Load by Cycles.csv")   # Data File Reading
    cycles = mtsdata["Cycles"]      # put Data into array
    cracklength = mtsdata["Crack Length (mm)"]
    kmax = mtsdata["Maximum K (kN/mm^1.5)"] * factor_for_k  # MPa.m^0.5
    kmin = mtsdata["Minimum K (kN/mm^1.5)"] * factor_for_k  # MPa.m^0.5
    pmax = mtsdata["Maximum Axial Force (kN)"] * 1e3  # N
    pmin = mtsdata["Minimum Axial Force (kN)"] * 1e3  # N
    if closure:
        kclosure = mtsdata["K Closure (kN/mm^1.5)"] * factor_for_k  # MPa.m^0.5
        closureload = mtsdata["Closure Load (kN/mm^1.5)"]    # kN, 推测MTS软件输出显示的单位有误
    if cod:
        codmax = mtsdata["Maximum Axial COD (mm)"]  # mm
        codmin = mtsdata["Minimum Axial COD (mm)"]  # mm
    if closure and cod:
        return cycles, cracklength, kmax, kmin, pmax, pmin, kclosure, closureload, codmax, codmin
    elif closure and (not cod):
        return cycles, cracklength, kmax, kmin, pmax, pmin, kclosure, closureload
    elif (not closure) and cod:
        return cycles, cracklength, kmax, kmin, pmax, pmin, codmax, codmin
    else:
        return cycles, cracklength, kmax, kmin, pmax, pmin


def ReadTestInf(sequence):
    # usage: 读取MTS输出的实验参数，包括试件名，尺寸参数（Width，notch_length，thickness），材料参数（弹性模量E，
    # 屈服强度YS），预制裂纹长度
    # input parameter:
    # sequence: 文件名前缀，通常由batch名和specimen名组成
    # return parameter:
    # 试件名specimen(str), 试件宽width，缺口长度notch_length，厚度thickness，弹性模量elastic_modulus，
    # 屈服强度yield_strength，预制裂纹长度precrack_length(均为float，单位见函数内部)
    mtsinf = pd.read_csv(path + sequence + u"_Test Summary.csv", index_col=0)  # Data File Reading
    specimen = mtsinf.loc["Specimen"]   # Material and Specimen Data From Summary File
    width = mtsinf.loc["Width (W)"]
    notch_length = mtsinf.loc["Notch Length (a0)"]
    thickness = mtsinf.loc["Thickness (B)"]
    elastic_modulus = mtsinf.loc["Elastic Modulus (E)"]
    yield_strength = mtsinf.loc["Yield Strength"]
    precrack_length = mtsinf.loc["Precrack crack length"]

    def datatransfer(data):
        s = data.str.split(" ")
        s1 = s['Unnamed: 1']
        return s1[0]

    specimen = datatransfer(specimen)
    width = float(datatransfer(width))  # mm
    notch_length = float(datatransfer(notch_length))  # mm
    thickness = float(datatransfer(thickness))  # mm
    elastic_modulus = float(datatransfer(elastic_modulus))  # GPa
    yield_strength = float(datatransfer(yield_strength))  # GPa
    precrack_length = float(datatransfer(precrack_length))  # mm

    return specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length


def ReadCodData(sequence):
    # usage: 读取COD和轴力的数据，文件名以Channels by Cycles结尾，包括指定循环内250组力与COD的数据
    # input parameter:
    # sequence: 试件Batch-specimen名
    # return parameter:
    # 循环次数数组、载荷数组、COD位移值数组
    cycles = []
    load = []       # kN
    cod = []        # mm
    for df in pd.read_csv(path + sequence + u"_Channels by Cycles.csv",
                          chunksize=250):       # Read 250 records at one time
        cycles.append(df["Nominal Cycle"])
        load.append(df["Axial Force (kN)"])
        cod.append(df["Axial COD (mm)"])
    cycles = np.array(cycles).reshape(-1, 1)
    load = np.array(load).reshape(-1, 1)
    cod = np.array(cod).reshape(-1, 1)
    return cycles, load, cod


def ReadTensionResult(sequence):
    # usage:读取简单拉伸实验结果数据
    # input parameter:
    # sequence: 试件序号
    # return parameter:
    # extension：位移/mm
    # load：载荷/N
    # stress：工程应力/MPa
    # strain：工程应变/MPa
    tensionresult = pd.read_csv(path + u"QSTE420TM_TensionTest_" + sequence + u".csv")  # Data File Reading
    extension = tensionresult["Extension (mm)"]
    load = tensionresult["Load (N)"]
    stress = tensionresult["Stress (MPa)"]
    strain = tensionresult["Strain (%)"] * 0.01
    return extension, load, stress, strain
