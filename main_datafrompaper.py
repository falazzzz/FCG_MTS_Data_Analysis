# LU YUNCHAO
# 2018/7/30 V1.0.0
# 由 Engauge Digitizer 软件读取的论文曲线点的csv文件读取数据作图

from FCGAnalysisLib import mixed_model
from FCGAnalysisLib import mts_analysis
import matplotlib.pyplot as plt

# 主程序
filename = "6082T6_CA.csv"
data = mixed_model.FigureDataReader(filename)
specimen = data[0]
c_ca, m_ca = mts_analysis.ParisFitting(dadn=specimen.dadn, dk=specimen.dk)
print("Paris Fitting Result: c_ca="+str(c_ca)+"mm/cycle, m_ca="+str(m_ca))

# 绘图
plt.figure(num=1, figsize=(10, 8))
for datum in data:
    plt.scatter(datum.dk, datum.dadn, label=datum.name, marker='.')
plt.axis([5, 30, 1e-5, 1e-2])
plt.xlabel("DeltaK Applied (MPa*m^0.5)")
plt.xscale('log')
plt.ylabel("da/dN (mm/cycle)")
plt.yscale('log')
plt.title('da/dN - dK ' + filename)
plt.grid(which='minor', linestyle='--')
plt.grid(which='major', linestyle='--')
plt.legend()
plt.grid()
plt.show()
