# 2019/1/3
# 为数据处理建立类，完成标准部分流程
# Lu Yunchao

from FCGAnalysisLib import read_data
from FCGAnalysisLib import experiment_calculation


class SpecimenBasic:
    def __init__(self, name, stress_ratio, threshold=0):
        # 基本参数，创建实例时输入
        self.name = name                    # MTS内部试件名, example: "yang-baoban_Lu-340-3-4"
        self.stress_ratio = stress_ratio    # 应力比
        self.threshold = threshold          # 数据处理可将SIF range低于阈值的数据删除
        # Test information，由ReadTestInf函数读入
        self.specimen = 0                   # 试件名
        self.width = 0                      # 试件宽度
        self.notch_length = 0               # 缺口长度
        self.thickness = 0                  # 试件厚度
        self.elastic_modulus = 0            # 弹性模量
        self.yield_strength = 0             # 屈服强度
        self.precrack_length = 0            # 预制裂纹长度
        # MTS计算数据，由ReadMtsResult函数读入
        self.dadn_MTS = []                   # 裂纹扩展速率
        self.n_MTS = []                      # 循环次数
        self.dk_MTS = []                     # SIF Range
        # MTS原始数据，由ReadOriginResult函数读入
        self.cycles = []
        self.cracklength = []
        self.kmax = []
        self.kmin = []
        self.pmax = []
        self.pmin = []
        self.kclosure = []
        self.closureload = []
        # 由MTS原始数据计算得到的结果，由FCGRandDKbyOriginalData函数计算
        self.dadn_Manual = []
        self.n_Manual = []
        self.dk_Manual = []
        self.a_Manual = []
        self.kc_Manual = []

    def read_experimental_data(self, read_test_inf=True, ys=0):
        # 读取实验结果，默认读取实验参数文件，也可选择不读，若不读需在执行该函数前先存入相关参数
        # 若传入ys参数不为零，则存数据时存为传入的yield strength值，输入单位MPa
        if read_test_inf:
            specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
                read_data.ReadTestInf(sequence=self.name)
            self.specimen = specimen
            self.width = width
            self.notch_length = notch_length
            self.thickness = thickness
            self.elastic_modulus = elastic_modulus
            if ys == 0:
                self.yield_strength = yield_strength
            else:
                self.yield_strength = ys * 1e-3     # TestInf记录单位为GPa
            self.precrack_length = precrack_length
        dadn_MTS, n_MTS, dk_MTS = read_data.ReadMtsResult(sequence=self.name)
        cycles, cracklength, kmax, kmin, pmax, pmin, kclosure, closureload = \
            read_data.ReadOriginResult(sequence=self.name, closure=True, cod=False)
        dadn_Manual, n_Manual, dk_Manual, a_Manual, kc_Manual = \
            experiment_calculation.FCGRandDKbyOriginalData(b=self.thickness, w=self.width,
                                                           n=cycles, pmax=pmax, pmin=pmin, a=cracklength,
                                                           ys=self.yield_strength, r=self.stress_ratio,
                                                           kclosure=kclosure, threshold=self.threshold)
        self.dadn_MTS = dadn_MTS
        self.n_MTS = n_MTS
        self.dk_MTS = dk_MTS

        self.cycles = cycles
        self.cracklength = cracklength
        self.kmax = kmax
        self.kmin = kmin
        self.pmax = pmax
        self.pmin = pmin
        self.kclosure = kclosure
        self.closureload = closureload

        self.dadn_Manual = dadn_Manual
        self.n_Manual = n_Manual
        self.dk_Manual = dk_Manual
        self.a_Manual = a_Manual
        self.kc_Manual = kc_Manual


class DoubleOLSpecimenBasic(SpecimenBasic):
    def __init__(self, name, stress_ratio, threshold, ol_1_cycle):
        super().__init__(name, stress_ratio, threshold)
        self.ol_1_cycle = ol_1_cycle        # 施加第一次高载时的循环数
        self.n_Fixed_by_ol_1 = []

    def ol_1_fix(self):
        self.n_Fixed_by_ol_1 = self.n_Manual - self.ol_1_cycle

