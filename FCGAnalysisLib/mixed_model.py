# mixed residual stress and crack closure model相关函数
# 2018.8.2 Version 1.0.0
# Lu Yunchao

from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import read_data
from FCGAnalysisLib import experiment_calculation
from FCGAnalysisLib import paris_and_walker
from FCGAnalysisLib import overload_analysis

import numpy as np
import pandas as pd
import sys


# 类定义
class OverloadSpecimen:
    # 由实验结果导数基本参数及实验数据
    # 需使用库read_data, experiment_calculation, mts_analysis, paris_and_walker
    # 包含函数：
    # 1.basic_result_calculate: 读取试件名为self.name的相应文件，读取试件材料和几何参数，并初始化计算得到实验数据
    # 2.load_condition_input：单独输入载荷峰值和应力比
    # 3.overload_status_input：单独输入高载值及各裂纹长度节点
    # 4.status_input_from_database：由paris_and_walker数据库读取主要参数
    # 5.specimen_print：输出类的实例已存储的数据，可用于确认数据是否正确
    # 6.data_select_by_cracklength：数据筛选函数。可设定筛选的裂纹长度上限 and/or 下限，输出相应的5类FCG实验数据

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
        self.a_retard_end = -1
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
        self.a_retard_end = database["a_retard_end"]

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


class FCGDataFromFigure:
    # 由论文中采点得到的数据
    # 参数order表示该数据在数据csv文件中的列号，由于数据通常包括不止一列，因此可选择读入的数据
    # 类包括一个函数
    # data_reader: 由csv读入的dataFrame格式中读取试件名、dk和dadn参数
    def __init__(self, order):
        self.order = order
        # FCG basic data
        self.name = ""
        self.dk = ""
        self.dadn = ""
        # Loading condition
        self.stress_ratio = 1
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
        self.a_retard_end = -1
        # Alpha Fitting Region
        self.a_alpha_begin = -1
        self.a_alpha_end = -1

    def data_reader(self, database):
        self.name = database.columns[self.order]
        self.dk = database[database.columns[0]]
        self.dadn = database[self.name]

    def data_select_by_sif(self, lowerlimit=-1, upperlimit=-1):
        dk = self.dk
        dadn = self.dadn
        if lowerlimit != -1:
            dadn = mts_analysis.DataSelectByThreshold(threshold=lowerlimit, parameter=dk, data=dadn, keepbigger=1)
            dk = mts_analysis.DataSelectByThreshold(threshold=lowerlimit, parameter=dk, data=dk, keepbigger=1)
        if upperlimit != -1:
            dadn = mts_analysis.DataSelectByThreshold(threshold=upperlimit, parameter=dk, data=dadn, keepbigger=0)
            dk = mts_analysis.DataSelectByThreshold(threshold=upperlimit, parameter=dk, data=dk, keepbigger=0)
        return np.array(dadn), np.array(dk)

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


def FigureDataReader(filename):
    # 利用定义的数据类，包括SIF和FCGR两项，类中定义函数读取从pandas Dataframe中读取数据
    # 函数主体为循环，循环内定义类，讲类存储在数组中.
    # 函数读取csv文件存入database的pandas dataframe中，再利用循环读取所有数据
    # 最后返回由类构成的数组
    inputfile = "\\inputfile\\figure_data\\"  # 读取文件放置的子文件夹名
    path = sys.path[0] + inputfile  # 设定读取的绝对路径
    database = pd.read_csv(path + filename)
    length = len(database.columns)
    data = []
    for seq in range(1, length):
        data.append(FCGDataFromFigure(order=seq))
        data[-1].data_reader(database=database)
    return data


def BiSelection(a, b, const, t, w, ys, a_ol_applied, pmax, plz_factor, accu=0.001):
    # 利用二分法计算主导位置变化
    # 裂纹自高载后的扩展长度+循环塑性区尺寸是已知值，但循环区尺寸是扩展长度的函数且不易解，故采用二分法求解
    # 求解的方程可表示为求解长度da使方程：
    # da + plz,reversing(da) - plz,ol * precentage(44% for OLR = 2, 80% for OLR = 1.5）= 0
    # input arguments:
    # a, b: 二分法的左右两端
    # const：高载塑性区尺寸乘比例（44% or 80%）
    # t, w, ys: 试件的thickness、width和yielding stress，可由类读取
    # a_ol_applied：试验高载时对应裂纹长度，可由类读取
    # pmax：试验中循环载荷峰值，可由类读取
    # plz_factor：塑性区尺寸，可输入常数或Xiaoping，若输入Xiaoping则会自动根据试件尺寸计算
    # accu：二分法精度，默认为0.001
    # return arguments:
    # result: 数组，第0位为解的值x，第1位为此时函数的值
    i = 0   # 计数器
    # 计算二分法左端的函数值
    ka = mts_analysis.DeltaKCalculating(b=t, w=w, a=a_ol_applied + a,pmax=pmax, pmin=0)
    plza = overload_analysis.PlasticZoneWithFactor(kmax=[ka], factor=plz_factor, ys=ys, t=t)
    ya = a + plza - const
    # 计算二分法右端的函数值
    kb = mts_analysis.DeltaKCalculating(b=t, w=w, a=a_ol_applied + b, pmax=pmax, pmin=0)
    plzb = overload_analysis.PlasticZoneWithFactor(kmax=[kb], factor=plz_factor, ys=ys, t=t)
    yb = b + plzb - const
    result = []
    # 若两端为0的情况：
    if ya == 0:
        result.append(a)
        result.append(ya)
    elif yb == 0:
        result.append(b)
        result.append(yb)
    else:
        # 若两端不为零->进入二分法
        while i >= 0:
            c = (a+b)/2
            # 计算中点的函数值
            kc = mts_analysis.DeltaKCalculating(b=t, w=w, a=a_ol_applied + c, pmax=pmax, pmin=0)
            plzc = overload_analysis.PlasticZoneWithFactor(kmax=[kc], factor=plz_factor, ys=ys, t=t)
            yc = c + plzc - const
            # 判断解的范围精度是否满足accu的设定
            if np.abs(b-a) < accu:
                result.append(c)
                result.append(yc)
                print('Find Zero_point due to tolrance at c='+str(c))
                break
            if yc == 0:
                result.append(c)
                result.append(yc)
                print('Find Zero_point due to Midpoint_value=0 at c='+str(c))
                break
            elif yc*ya < 0:
                b = c
                i = i + 1
                continue
            elif yc*yb < 0:
                i = i + 1
                a = c
                continue
            else:
                print('ERROR: Unexpected situation come across.Check the process.')
                break
    return result


def LogLinearFit(x, y):
    # 输入x数组和y数组，对数化后线性拟合，返回斜率和截距
    # 目前用于mixed_model迟滞段的Kop拟合
    # 对数化后进行线性拟合
    # 返回斜率和截距
    if len(x) > 2:
        pp = np.polyfit(x=np.log(x), y=np.log(y), deg=1)
        k, b = pp
    elif len(x) == 2:
        x1, x2 = np.log(x)
        y1, y2 = np.log(y)
        k = (y2 - y1)/(x2 - x1)
        b = (y1*x2 - y2*x1)/(x2 - x1)
    print("NoticeFromFunction'LogLinearFit':Kop linear fit:k="+str(k)+",b="+str(b))
    return k, b


def LogLinearFittedCalculate(k, b, x):
    # 输入在对数下拟合的斜率和截距参数，输入x数组，计算输出y
    # 目前用于mixed_model迟滞段的Kop拟合
    # 计算x对应的y
    y = np.exp(b)*x**k
    return y


def BiCalculateCracklength(a, b, dk, t, w, pmax, r, accu=0.001):
    # 对于CT试件，依据E647标准，已知SIF时利用二分法计算裂纹长度
    # 裂纹自高载后的扩展长度+循环塑性区尺寸是已知值，但循环区尺寸是扩展长度的函数且不易解，故采用二分法求解
    # 求解的方程可表示为求解长度a使方程：
    # dk(a) - dk(输入值) = 0
    # input arguments:
    # a, b: 二分法的左右两端(裂纹长度可能存在的范围)
    # dk：已知的SIF值
    # t, w: 试件的thickness、width，可由类读取
    # pmax, r：试验中循环载荷峰值和应力比，可由类读取
    # accu：二分法精度，默认为0.001
    # return arguments:
    # result: 数组，第0位为解的值x，第1位为此时函数的值
    i = 0   # 计数器
    # 计算二分法左端的函数值
    ya = mts_analysis.DeltaKCalculating(b=t, w=w, a=a, pmax=pmax, pmin=pmax * r) - dk
    # 计算二分法右端的函数值
    yb = mts_analysis.DeltaKCalculating(b=t, w=w, a=b, pmax=pmax, pmin=pmax * r) - dk
    result = []
    # 若两端为0的情况：
    if ya == 0:
        result.append(a)
        result.append(ya)
    elif yb == 0:
        result.append(b)
        result.append(yb)
    else:
        # 若两端不为零->进入二分法
        while i >= 0:
            c = (a+b)/2
            # 计算中点的函数值
            yc = mts_analysis.DeltaKCalculating(b=t, w=w, a=c, pmax=pmax, pmin=pmax * r) - dk
            # 判断解的范围精度是否满足accu的设定
            if np.abs(b-a) < accu:
                result.append(c)
                result.append(yc)
                print('Find Zero_point due to tolrance at c='+str(c))
                break
            if yc == 0:
                result.append(c)
                result.append(yc)
                print('Find Zero_point due to Midpoint_value=0 at c='+str(c))
                break
            elif yc*ya < 0:
                b = c
                i = i + 1
                continue
            elif yc*yb < 0:
                i = i + 1
                a = c
                continue
            else:
                print('ERROR: Unexpected situation come across.Check the process.')
                break
    return result
