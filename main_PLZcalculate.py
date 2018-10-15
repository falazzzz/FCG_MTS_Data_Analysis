# 计算塑性区大小
# 2018/10/15 Version 1.0
# Lu Yunchao

from FCGAnalysisLib import overload_analysis
from FCGAnalysisLib import mts_analysis
from FCGAnalysisLib import paris_and_walker
import numpy as np

sequence = ["yang-baoban_Lu-420-05"]
specimen_data = paris_and_walker.test_database2(specimen=sequence[0])
ys = 0.365          # GPa
a_ol = np.array([specimen_data["a_ol_applied"]])      # mm
p_ol = specimen_data["overload"]         # N
p_max = specimen_data["maxload"]
a_kc_max = specimen_data["a_kc_max"]
print('specimen:'+sequence[0]+',a_ol_applied:'+str(a_ol)+'mm,overload:'+str(p_ol)+'N, maxload:'+str(p_max))
b = 2.5             # mm
w = 40              # mm

k_ol = np.array([mts_analysis.DeltaKCalculating(b=b, w=w, a=a_ol, pmax=p_ol, pmin=0)])
k_kc_max = np.array([mts_analysis.DeltaKCalculating(b=b, w=w, a=a_kc_max, pmax=p_max, pmin=0)])
factor = overload_analysis.FactorXiaoping(kmax=k_ol, t=b, ys=ys)
print(factor)
print(str(1/np.pi))
r_ol = overload_analysis.PlasticZoneWithFactor(kmax=k_ol, ys=ys, factor='Xiaoping', t=b)
r_kc_max = overload_analysis.PlasticZoneWithFactor(kmax=k_kc_max, ys=ys, factor='Xiaoping', t=b)
print("r_ol:", r_ol)
print("r_kc_max:", r_kc_max)