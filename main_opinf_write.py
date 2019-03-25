# 2018/11/30
# Lu Yunchao
# 脚本完成以下工作：
# 1. 由文件中读取文件，读取其Closure Load和Kop，读取裂纹长度a，计算对应的dk
# 2. 绘图
# 3. 将以上4项存储至csv文件中

from FCGAnalysisLib import mixed_model
from FCGAnalysisLib import write_data
from FCGAnalysisLib import closureanalysis

import numpy as np
import matplotlib.pyplot as plt

'''
*****基本参数设定*****
'''
sequence = ["yang-baoban_Lu-420-08", "yang-baoban_Lu-420-09"]
moving_average = [1, 1]     # 移动平均开关，针对每个试件
numofaverage = 7        # 移动平均参数：移动平均区间
ratio = 0.5             # 移动平均参数：误差带宽度
figshow = 1        # 图像显示开关
datasave = 0    # 保存数据


'''
*****初始化*****
'''
# 创建实例
specimens = []
for name in sequence:
    specimens.append(mixed_model.OPInfOverloadSpecimen(name=name))
# 读取数据
for specimen in specimens:
    specimen.status_input_from_database()
    specimen.basic_result_calculate()
    specimen.specimen_print()

'''
*****数据移动平均*****
'''
for i, specimen in enumerate(specimens):
    if moving_average[i]:
        specimen.ma_fop, specimen.ma_dk = \
            closureanalysis.DoubleDataSelectMovingAverage(data=specimen.fop, reference=specimen.dk,
                                                          ratio=ratio, numofaverage=numofaverage)
        specimen.ma_a = specimen.a[numofaverage-2:]
        specimen.ma_kop = specimen.kc[numofaverage-2:]
        specimen.ma_dadn = specimen.dadn[numofaverage-2:]
    else:
        specimen.ma_fop = specimen.fop
        specimen.ma_dk = specimen.dk
        specimen.ma_a = specimen.a
        specimen.ma_kop = specimen.kc
        specimen.ma_dadn = specimen.dadn


'''
*****数据存储*****
'''
if datasave:
    for specimen in specimens:
        filename = specimen.name + "_closure_load_inf"
        data = np.array([specimen.ma_a, specimen.ma_dk, specimen.ma_kop, specimen.ma_fop, specimen.ma_dadn])
        name = ['crack length', 'dK', 'Kop', 'closure load', 'dadn']
        _ = write_data.SaveData(dataset=data, name=name, filename=filename)


'''
*****绘图*****
'''
'''
plt.figure(num=1, figsize=(8, 6))   # kop-dk曲线
for specimen in specimens:
    plt.scatter(specimen.dk, specimen.kc, marker='.', label=specimen.name)
plt.title("Kop-dK Diagram")
plt.ylabel("Kop/MPa.m0.5")
plt.xlabel("dK/MPa.m0.5")
plt.xscale('log')
plt.yscale('log')
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()

plt.figure(num=2, figsize=(8, 6))   # kop-a曲线
for specimen in specimens:
    plt.scatter(specimen.a, specimen.kc, marker='.', label=specimen.name)
plt.title("Kop-a Diagram")
plt.ylabel("Kop/MPa.m0.5")
plt.xlabel("crack length/mm")
#plt.xscale('log')
#plt.yscale('log')
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()

plt.figure(num=3, figsize=(8, 6))   # fop-a曲线
for specimen in specimens:
    plt.scatter(specimen.a, specimen.fop, marker='.', label=specimen.name)
plt.title("Fop-a Diagram")
plt.ylabel("Fop/kN")
plt.xlabel("crack length/mm")
#plt.xscale('log')
#plt.yscale('log')
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()
'''
plt.figure(num=4, figsize=(8, 6))   # fop-dk曲线
for specimen in specimens:
    plt.scatter(specimen.ma_dk, specimen.ma_fop, marker='.', label=specimen.name)
plt.title("Fop-dK Diagram")
plt.ylabel("Fop/kN")
plt.xlabel("dK/MPa.m0.5")
plt.xscale('log')
#plt.yscale('log')
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()

if figshow:
    plt.show()

