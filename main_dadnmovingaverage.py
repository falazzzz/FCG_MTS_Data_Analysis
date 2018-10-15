# 2018/7/2 Version 2.0
# 对结果进行移动平均，基于封装API进行了部分优化
# Lu Yunchao

from FCGAnalysisLib import closureanalysis
from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import paris_and_walker
import matplotlib.pyplot as plt
import numpy as np

# parameter
sequence = ["yang-baoban_Lu-420-09"]
r = 0.7
threshold = 0

# 材料Paris参数和Walker参数
c_w, m_w, gamma = paris_and_walker.WalkerParameter()
c_p, m_p = paris_and_walker.ParisParameter(r=r)

# 数据读取
dadn_Manual, n_Manual, dk_Manual, a_Manual, kc_Manual = closureanalysis.BasicDataDispose(sequence=sequence[0], r=r)
dkeff_Manual = closureanalysis.dKeffCalculation(kmax=dk_Manual/(1-r), kclose=kc_Manual, kmin=dk_Manual*(r/(1-r)),
                                                method='basic')

# moving average
numofaverage = 7
ratio = 0.5
dadn_Manual_averaged, dk_Manual_averaged = \
    closureanalysis.DoubleDataSelectMovingAverage(data=dadn_Manual, reference=dk_Manual,
                                                  ratio=ratio, numofaverage=numofaverage)

# 无高载Paris虚线
dadn_paris_Manual = mts_analysis.ParisCalculating(c=c_p, m=m_p, dk=dk_Manual)

# ploting
plt.figure(num=1, figsize=(10, 8))
plt.scatter(dk_Manual, dadn_Manual, label='$ExperimentData$', color='red', marker='.')
plt.scatter(dk_Manual_averaged, dadn_Manual_averaged, label='$MovingAverage$', color='blue', marker='.')
plt.plot(dk_Manual, dadn_paris_Manual, label='$Paris for CA$', color='black', linewidth=2)
plt.axis([min(dk_Manual) * 0.95, max(dk_Manual) * 1.05, min(dadn_Manual) * 0.95, max(dadn_Manual) * 1.05])
plt.xlabel("DeltaK Applied (MPa*m^0.5)")
plt.xscale('log')
plt.ylabel("da/dN (mm/cycle)")
plt.yscale('log')
plt.xticks(np.linspace(min(dk_Manual), max(dk_Manual), 6))
plt.title('da/dN - dK ' + sequence[0] + '(MTS Result)')
plt.legend()
plt.grid()
plt.show()
