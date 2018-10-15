# 保留
# rename: main_predict2
# 预测CTS模型试件的裂纹扩展，采用class模块

import numpy as np
import matplotlib.pyplot as plt

from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import paris_and_walker
from FCGAnalysisLib import experiment_predict


class PredictCTS:
    # Predict CTS试件的裂纹扩展模型, 适用QSTE340TM材料
    # 输入角度为角度制，内部转化为弧度制计算
    # 力单位为N，长度mm，应力GPa
    def __init__(self, load, stressratio, width, thickness, alpha, ys=0.365):
        self.load = load
        self.stressratio = stressratio
        self.width = width
        self.thickness = thickness
        self.alpha = alpha / 180 * np.pi
        self.yieldstrength = ys
        self.cracklength = []
        self.sif_range_k1 = []
        self.sif_range_k2 = []
        self.sif_range_eff = []
        self.fcg_rate = []
        self.cycle = []
        self.valid = []

    def set_crack_length(self, cracklength):
        self.cracklength = np.array(cracklength)

    def evaluate_sif_range_k1(self):
        self.sif_range_k1 = \
            mts_analysis.CTS_Richard_K1_calculate(f=self.load, r=self.stressratio, a=self.cracklength,
                                                  w=self.width, b=self.thickness, alpha=self.alpha)
        return self.sif_range_k1

    def evaluate_sif_range_k2(self):
        self.sif_range_k2 = \
            mts_analysis.CTS_Richard_K2_calculate(f=self.load, r=self.stressratio, a=self.cracklength,
                                                  w=self.width, b=self.thickness, alpha=self.alpha)
        return self.sif_range_k2

    def calculate_sif_range_eff(self, method='tanaka'):
        # 默认采用tanaka方法计算I-II型混合应力强度因子，未来可继续添加其它方法.
        if method == 'tanaka':
            self.sif_range_eff = (self.sif_range_k1**4 + 8*self.sif_range_k2**4)**0.25
        return self.sif_range_eff

    def evaluate_fcg_rate_by_paris(self, c_p, m_p):
        self.fcg_rate = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=self.sif_range_eff)
        return self.fcg_rate

    def evaluate_cycle(self, c_p, m_p):
        self.cycle = \
            experiment_predict.CycleIntegrateBySimpson(c=c_p, m=m_p, dk=self.sif_range_eff,
                                                       cracklength=self.cracklength)
        return self.cycle

    def predict(self, c_p, m_p, method='tanaka' ):
        # usage：完成预测过程，输出SIF变幅、FCGR、循环次数、有效性判定
        # NOTE：参数全部引用自函数外部，不适合从其它文件调用
        # input parameter：
        # load：实验施加的载荷峰值，应当为1个数，可从外部的数组中取一个调用
        # return parameter：
        # cracklength：给定的裂纹长度，mm
        # dk_predict：计算的SIF变幅，Mpa.m0.5
        # dadn：估计的FCGR，mm/cycles
        # valid：有效性判断结果
        _ = self.evaluate_sif_range_k1()
        _ = self.evaluate_sif_range_k2()
        _ = self.calculate_sif_range_eff(method=method)
        self.fcg_rate = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=self.sif_range_eff)
        self.cycle = experiment_predict.CycleIntegrateBySimpson(c=c_p, m=m_p, dk=self.sif_range_eff,
                                                                cracklength=self.cracklength)
        for seq, value in enumerate(self.cracklength):
            self.valid.append(
                mts_analysis.LigamentValidCheck(w=self.width, a=value, dk=self.sif_range_eff[seq],
                                                ys=self.yieldstrength, r=self.stressratio))
        return self.cracklength, self.sif_range_eff, self.fcg_rate, self.cycle, self.valid


# 基础参数
load = 2000
alphas = [0, 30, 60, 90]
cracklist = np.arange(18, 31, 0.25)
c_p, m_p = paris_and_walker.ParisParameter(r=0.1)

# 程序参数
show = 1      # 绘图开关

# 各组计算
predicts = []
fail = []
for alpha in alphas:
    predicts.append(PredictCTS(load=load, stressratio=0.1, width=40, thickness=2.5, alpha=alpha))
    predicts[-1].set_crack_length(cracklength=cracklist)
    a, dk, dadn, n, valid = predicts[-1].predict(c_p=c_p, m_p=m_p)
    fail.append(experiment_predict.PredictPrint(load=2000, a=a, dk=dk, cycles=n, dadn=dadn, valid=valid, alpha=alpha))

# 绘图部分
# DeltaK - CrackLength Plotting
plt.figure(num=1, figsize=(7, 5))
for seq, predict in enumerate(predicts):
    plt.plot(predict.cracklength, predict.sif_range_eff, label='alpha='+str(predict.alpha), linewidth=2)
plt.xlabel("a/mm")
plt.xlim(np.min(cracklist), np.max(cracklist))
plt.ylabel("SIF/MPa.m0.5")
plt.title('DeltaK - CrackLength,load='+str(load))
plt.legend()
plt.grid()

# Cycle - CrackLength Plotting
plt.figure(num=2, figsize=(7, 5))
for seq, predict in enumerate(predicts):
    plt.plot(predict.cycle, predict.cracklength, label='alpha='+str(predict.alpha), linewidth=2)
plt.ylabel("Cracklength/mm")
plt.xlabel("Cycle")
plt.title('Cycle - CrackLength,load='+str(load))
plt.legend()
plt.grid()

if show:
    plt.show()
