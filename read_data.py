# read_data部分，负责读取文件并整理格式，返回数组

import math
import pandas as pd

factor_for_k = 1e3 * math.sqrt(1e-3)  # convert kN.mm**-2 to Mpa.m**-0.5


def ReadMtsResult(sequence):
    # usage: 读取MTS输出的拟合计算结果，共有裂纹扩展速率da/dN，循环数cycles，SIF变幅dk三项
    # input parameter:
    # sequence: 文件名前缀，通常由batch名和specimen名组成
    # return parameter:
    # 裂纹扩展速率dadn，循环数cycles，SIF变幅dk
    mtsresult = pd.read_csv(sequence + u"_da dN, Delta K by Cycles.csv")  # Data File Reading
    dadn = mtsresult["da/dN (mm/cycle)"]
    cycles = mtsresult["Cycles"]
    dk = mtsresult["Delta K Applied (kN/mm^1.5)"] * factor_for_k
    return dadn, cycles, dk


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
    mtsdata = pd.read_csv(sequence + u"_Crack Length, Min Max K, Load by Cycles.csv")   # Data File Reading
    cycles = mtsdata["Cycles"]      # put Data into array
    cracklength = mtsdata["Crack Length (mm)"]
    kmax = mtsdata["Maximum K (kN/mm^1.5)"] * factor_for_k  # MPa.m^0.5
    kmin = mtsdata["Minimum K (kN/mm^1.5)"] * factor_for_k  # MPa.m^0.5
    pmax = mtsdata["Maximum Axial Force (kN)"] * 1e3  # N
    pmin = mtsdata["Minimum Axial Force (kN)"] * 1e3  # N
    if closure:
        kclosure = mtsdata["K Closure (kN/mm^1.5)"] * factor_for_k  # MPa.m^0.5
        closureload = mtsdata["Closure Load (kN/mm^1.5)"] * factor_for_k    # MPa.m^0.5
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
    mtsinf = pd.read_csv(sequence + u"_Test Summary.csv", index_col=0)  # Data File Reading
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