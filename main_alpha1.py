# 2018/7/2 Version 1.0
# 用于alpha模型中寻找最佳衰减速率的命令
# 幅度指定为2/pi，也可通过改为0，改为0时将通过拟合确定适合各实验组的最佳幅值
# Lu Yunchao

from FCGAnalysisLib import closureanalysis
from FCGAnalysisLib import paris_and_walker
import numpy as np
import matplotlib.pyplot as plt


sequence = ["yang-baoban_Lu-420-05",
            "yang-baoban_Lu-420-06",
            "yang-baoban_Lu-420-18",
            "yang-baoban_Lu-420-19",
            "yang-baoban_Lu-420-15",
            "yang-baoban_Lu-420-16"
            ]
error_ratio = 0.2
numofaverage = 7
fade_rates = np.arange(7, 12, 0.1)
amplitude = 2/np.pi
show = 1

lost_list = []
for fade_rate in fade_rates:
    lost = []
    for seq, value in enumerate(sequence):
        stress_ratio, overload_ratio, a_range = paris_and_walker.testinf(sequence[seq])
        D, n, dk, dkeff_alpha, dadn, alpha, dkeff_exp = closureanalysis.Alpha2(sequence=sequence[seq],
                                                                               stress_ratio=stress_ratio,
                                                                               error_ratio=error_ratio,
                                                                               numofaverage=numofaverage,
                                                                               a_range=a_range,
                                                                               fade_rate=fade_rate,
                                                                               amplitude=amplitude,
                                                                               fitting=1, model='paris')
        lost_seq = np.sum((dkeff_alpha - dkeff_exp) ** 2) / len(dkeff_alpha)
        lost.append(lost_seq)
    lost_list.append(np.sum(lost))
    print('指定的衰减速率:', fade_rate, '总损失值:', np.sum(lost), '损失值：', lost)

plt.figure(num=1, figsize=(10, 6))
plt.plot(fade_rates, lost_list, label='lost in different fade rates', color='black', linewidth=2)
plt.ylabel("lost")
plt.xlabel("fade_rate")
plt.legend()
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
if show:
    plt.show()
