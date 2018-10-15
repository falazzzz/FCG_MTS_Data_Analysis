# 显微镜输出csv文件数据处理
# V1.0.0, 2018/7/18, LU YUNCHAO
from FCGAnalysisLib import paris_and_walker

import pandas as pd
import sys
import matplotlib.pyplot as plt


def MicroscopeDataReader(filename):
    # 数据读取
    inputfile = "\\inputfile\\microscope_data\\"  # 读取文件放置的子文件夹名
    path = sys.path[0] + inputfile  # 设定读取的绝对路径
    csvdata = pd.read_csv(path + filename, skiprows=3)
    xdata = csvdata["X"]
    ydata = csvdata["Y"]
    zdata = csvdata["Height(mm)"]
    return [xdata, ydata, zdata]


def CrossSectionDataReader(specimenname, location):
    # 为本次实验订制的输出函数
    filenames = [specimenname+"A_"+location+"_data.csv", specimenname+"B_"+location+"_data.csv"]
    x_coordinate = []
    height = []
    for filename in filenames:
        data = MicroscopeDataReader(filename=filename)
        x_coordinate.append(data[0])
        height.append(data[2])
    return x_coordinate, height


# 控制参数
show = 1

# 输入参数
specimen = "19"
location = "edge"

# 数据读取
x_coor, z_coor = CrossSectionDataReader(specimenname=specimen, location=location)
a_inf = paris_and_walker.test_database2(specimen="yang-baoban_Lu-420-"+specimen)
a_ol_applied = a_inf["a_ol_applied"]
a_dadn_max = a_inf["a_dadn_max"]
a_dadn_min = a_inf["a_dadn_min"]
a_kc_max = a_inf["a_kc_max"]

# 数据处理
# B侧翻面及上修
z_coor[1] = -z_coor[1]
z_coor[1] = z_coor[1] - min(z_coor[1])
z_coor[1] += 0.03
# A侧下修
z_coor[0] = z_coor[0] - min(z_coor[0])

# 绘图关键位置标识
ol_site = max(x_coor[0])/2 - 1.2
dadn_max_site = (a_dadn_max - a_ol_applied) + ol_site
dadn_min_site = (a_dadn_min - a_ol_applied) + ol_site
kc_max_site = (a_kc_max - a_ol_applied) + ol_site

# 绘图部分
plt.figure(num=1, figsize=(10, 8))
plt.plot(x_coor[0], z_coor[0], label="A side", linewidth=2)
plt.plot(x_coor[1], z_coor[1], label="B side", linewidth=2)
plt.plot([ol_site, ol_site], [min(z_coor[0]), max(z_coor[1])], label="Overload Site", linestyle='--')
plt.plot([dadn_max_site, dadn_max_site], [min(z_coor[0]), max(z_coor[1])], label="FCGR Max Site", linestyle='--')
plt.plot([dadn_min_site, dadn_min_site], [min(z_coor[0]), max(z_coor[1])], label="FCGR Min Site", linestyle='--')
if a_kc_max != -1:
    plt.plot([kc_max_site, kc_max_site], [min(z_coor[0]), max(z_coor[1])], label="Kop Max Site", linestyle='--')
plt.xlabel("Location/mm")
plt.ylabel("Height/mm")
plt.title("Height of cross section of specimen:"+specimen+",at location:"+location)
plt.legend()
plt.grid()
if show:
    plt.show
