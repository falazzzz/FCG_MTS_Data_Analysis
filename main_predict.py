# CT件实验裂纹扩展实验参数预测函数
# 主要内容：
# (1)根据E647-11标准对SIF变幅进行估算；
# (2)根据Paris公式对裂纹扩展速率（FCGR）进行预测；
# (3)采用Simpson插值积分公式由FCGR估算裂纹扩展速率；
# (4)根据E647-11标准对韧带长度有效性要求确定有效的裂纹扩展长度


import matplotlib.pyplot as plt
from FCGAnalysisLib import mts_analysis
import numpy as np

from FCGAnalysisLib import experiment_predict

# 估计参数
c = 8.0448e-10       # Walker模型参数，单位 mm/cycles
m = 3.6655           # Walker模型参数，无量纲
gamma = 0.9055
width = 39.99                  # CT试件尺寸Width，单位 mm
thickness = 2.54             # CT试件尺寸厚度，单位 mm
start = 10                  # Precrack结束长度a0，单位 mm
final = 21                  # 估算的最大裂纹长度，单位 mm
load = [3200, 3300, 3400, 4000]         # 计算的载荷峰值，增加或删减需要对下方的主程序部分进行对应增减，单位 N
stress_ratio = 0.7          # 应力比
step = 0.25                    # 预估的裂纹扩展步长，单位 mm
yield_strength = 0.446      # 材料的屈服强度，单位 Gpa


# 绘图参数
sequence = 'R=0.7_QSTE420TM'            # 保存文件名和绘图标题名
save = 0            # 保存开关
show = 1            # 显示开关


# 函数部分
def predict(load):
    # usage：完成预测过程，输出SIF变幅、FCGR、循环次数、有效性判定
    # NOTE：参数全部引用自函数外部，不适合从其它文件调用
    # input parameter：
    # load：实验施加的载荷峰值，应当为1个数，可从外部的数组中取一个调用
    # return parameter：
    # cracklength：给定的裂纹长度，mm
    # dk_predict：计算的SIF变幅，Mpa.m0.5
    # dadn：估计的FCGR，mm/cycles
    # valid：有效性判断结果
    cracklength = np.arange(start, final+step, step)
    dk_predict = experiment_predict.DeltaKPredict(w=width, thickness=thickness, start=start, final=final, load=load, r=stress_ratio, step=step)
    dadn = mts_analysis.WalkerCalculating(c0=c, m0=m, gamma=gamma, dk=dk_predict, r=stress_ratio)
    cycle = experiment_predict.CycleIntegrateBySimpson(c=c, m=m, dk=dk_predict, cracklength=cracklength)
    valid = []
    for seq, value in enumerate(cracklength):
        valid.append(
            mts_analysis.LigamentValidCheck(w=width, a=value, dk=dk_predict[seq], ys=yield_strength, r=stress_ratio))
    return cracklength, dk_predict, dadn, cycle, valid


# 主程序
a1, dk1, dadn1, cycle1, valid1 = predict(load[0])
invalidseq1 = experiment_predict.PredictPrint(load[0], a=a1, dk=dk1, cycles=cycle1, dadn=dadn1, valid=valid1)

a2, dk2, dadn2, cycle2, valid2 = predict(load[1])
invalidseq2 = experiment_predict.PredictPrint(load[1], a=a2, dk=dk2, cycles=cycle2, dadn=dadn2, valid=valid2)

a3, dk3, dadn3, cycle3, valid3 = predict(load[2])
invalidseq3 = experiment_predict.PredictPrint(load[2], a=a3, dk=dk3, cycles=cycle3, dadn=dadn3, valid=valid3)

a4, dk4, dadn4, cycle4, valid4 = predict(load[3])
invalidseq4 = experiment_predict.PredictPrint(load[3], a=a4, dk=dk4, cycles=cycle4, dadn=dadn4, valid=valid4)


# 绘图部分
# DeltaK - CrackLength Plotting
plt.figure(num=1, figsize=(7, 5))
plt.plot(a1, dk1, label='$P_m = $'+str(load[0])+'$N$', color='red', linewidth=2)
plt.plot(a2, dk2, label='$P_m = $'+str(load[1])+'$N$', color='blue', linewidth=2)
plt.plot(a3, dk3, label='$P_m = $'+str(load[2])+'$N$', color='purple', linewidth=2)
plt.plot(a4, dk4, label='$P_m = $'+str(load[3])+'$N$', color='orange', linewidth=2)
plt.scatter([a1[invalidseq1], a2[invalidseq2], a3[invalidseq3], a4[invalidseq4]],
            [dk1[invalidseq1], dk2[invalidseq2], dk3[invalidseq3], dk4[invalidseq4]],
            marker='*', label='Invalid Points')
plt.xlabel("a/mm")
plt.xlim(start, final)
plt.ylabel("SIF/MPa.m0.5")
plt.title('DeltaK - CrackLength_'+sequence)
plt.legend()
plt.grid()
if save:
    plt.savefig('dk_a'+sequence+'.png', dpi=320)
if show:
    plt.show()


# FCGR - CrackLength Plotting
plt.figure(num=2, figsize=(7, 5))
plt.plot(dk1, dadn1, label='$P_m = $'+str(load[0])+'$N$', color='red', linewidth=2)
plt.plot(dk2, dadn2, label='$P_m = $'+str(load[1])+'$N$', color='blue', linewidth=2)
plt.plot(dk3, dadn3, label='$P_m = $'+str(load[2])+'$N$', color='purple', linewidth=2)
plt.plot(dk4, dadn4, label='$P_m = $'+str(load[3])+'$N$', color='orange', linewidth=2)
plt.scatter([dk1[invalidseq1], dk2[invalidseq2], dk3[invalidseq3], dk4[invalidseq4]],
            [dadn1[invalidseq1], dadn2[invalidseq2], dadn3[invalidseq3], dadn4[invalidseq4]],
            marker='*', label='Invalid Points')
plt.xlabel("SIF Range/MPa.m0.5")
plt.axis([min(dk1), dk4[invalidseq4]*1.1, min(dadn1), dadn4[invalidseq4]*1.1])
plt.ylabel("FCG rate/mm per cycle")
plt.xscale('log')
plt.yscale('log')
plt.title('FCGR - SIF Range_'+sequence)
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()
plt.grid()
if save:
    plt.savefig('dadn_dk'+sequence+'.png', dpi=320)
if show:
    plt.show()
