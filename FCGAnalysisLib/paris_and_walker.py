# 2018/6/6
# Version 1.0
# Lu Yunchao
# QSTE340TM, Paris公式和Walker模型的参数库
# 该文件不公开
import numpy as np


def ParisParameter(r):
    # 输入应力比r，返回Paris参数c和m
    # 输入库中没有的应力比将返回0
    paris_c = {
        0.1: 7.9713e-10,
        0.3: 2.3232e-09,
        0.5: 0.9627e-09,
        0.7: 5.6665e-09
    }
    paris_m = {
        0.1: 3.6797,
        0.3: 3.3784,
        0.5: 3.7073,
        0.7: 3.0965
    }
    try:
        c_p = paris_c[r]
        m_p = paris_m[r]
    except KeyError:
        print("Error Stress Ratio Input: No match result in database, returning 0.")
        c_p = 0
        m_p = 0
    return c_p, m_p


def WalkerParameter(data=0):
    # 调用返回Walker模型的参数c、m和gamma
    # data=0返回由r=0.1和r=0.7拟合的数据（默认）
    # data=1返回由r=0.1,r=0.3,r=0.5拟合的数据
    if data == 0:
        c_w = 8.0448e-10
        m_w = 3.6655
        gamma = 0.9055
    elif data == 1:
        c_w = 8.2791e-10
        m_w =3.6531
        gamma = 0.8743

    return c_w, m_w, gamma


# alpha模型拟合函数参数
def testinf(testname):
    # 输入实验试件名，返回应力比、alpha模型拟合范围
    test_information = {
        "yang-baoban_Lu-420-05": {
            "stress_ratio": 0.1,
            "overload_ratio": 1.5,
            "a_range": np.array([13.1121, 18.3370])     # [12.8016, 18.3370]
        },
        "yang-baoban_Lu-420-06": {
            "stress_ratio": 0.1,
            "overload_ratio": 2.0,
            "a_range": np.array([13.1361, 18.3370])     # [12.9516, 18.3370]
        },
        "yang-baoban_Lu-420-15": {
            "stress_ratio": 0.5,
            "overload_ratio": 1.5,
            "a_range": np.array([13.6319, 16.3370])     # [13.3016, 16.3370]
        },
        "yang-baoban_Lu-420-16": {
            "stress_ratio": 0.5,
            "overload_ratio": 2.0,
            "a_range": np.array([13.8580, 18.3370])     # [13.8016, 18.3370]
        },
        "yang-baoban_Lu-420-18": {
            "stress_ratio": 0.3,
            "overload_ratio": 1.5,
            "a_range": np.array([13.1596, 17.3370])     # [12.8016, 17.3370]
        },
        "yang-baoban_Lu-420-19": {
            "stress_ratio": 0.3,
            "overload_ratio": 2.0,
            "a_range": np.array([13.4774, 18.3370])     #
        }
    }
    try:
        test = test_information[testname]
    except KeyError:
        print("Test Information was not collected in the database. Returning -1.")
        stress_ratio = -1
        overload_ratio = -1
        a_range = np.array([-1, -1])
    stress_ratio = test["stress_ratio"]
    overload_ratio = test["overload_ratio"]
    a_range = test["a_range"]
    return stress_ratio, overload_ratio, a_range


# alpha_modified wheeler模型拟合参数库
def test_database2(specimen):
    # 输入试件名，返回相关数据的字典，由类函数直接读取字典
    test_information = {
        "yang-baoban_Lu-420-05": {"stress_ratio": 0.1, "maxload": 2000, "overload": 3000,
                                  "a_ol_applied": 11.9961, "a_dadn_max": 12.0537, "a_dadn_min": 12.4761,
                                  "a_kc_max": 13.1121, "a_alpha_end": 18.3370, "a_retard_end": 16.1509},
        "yang-baoban_Lu-420-06": {"stress_ratio": 0.1, "maxload": 2000, "overload": 4000,
                                  "a_ol_applied": 12.0057, "a_dadn_max": 12.0190, "a_dadn_min": 12.4516,
                                  "a_kc_max": 13.1361, "a_alpha_end": 18.3370, "a_retard_end": 16.9345},
        "yang-baoban_Lu-420-18": {"stress_ratio": 0.3, "maxload": 2400, "overload": 3600,
                                  "a_ol_applied": 12.0225, "a_dadn_max": 12.0660, "a_dadn_min": 12.5081,
                                  "a_kc_max": 13.1596, "a_alpha_end": 17.3370, "a_retard_end": 14.7394},
        "yang-baoban_Lu-420-19": {"stress_ratio": 0.3, "maxload": 2400, "overload": 4800,
                                  "a_ol_applied": 12.0225, "a_dadn_max": 12.1313, "a_dadn_min": 12.6006,
                                  "a_kc_max": 13.4774, "a_alpha_end": 18.3370, "a_retard_end": 16.3330},
        "yang-baoban_Lu-420-15": {"stress_ratio": 0.5, "maxload": 2800, "overload": 4200,
                                  "a_ol_applied": 12.0116, "a_dadn_max": 12.0770, "a_dadn_min": 12.6605,
                                  "a_kc_max": 13.6319, "a_alpha_end": 16.3370, "a_retard_end": 16.6605},
        "yang-baoban_Lu-420-16": {"stress_ratio": 0.5, "maxload": 2800, "overload": 5600,
                                  "a_ol_applied": 12.0133, "a_dadn_max": 12.2137, "a_dadn_min": 13.3306,
                                  "a_kc_max": 13.8580, "a_alpha_end": 18.3370, "a_retard_end": 17.6820},
        "yang-baoban_Lu-420-08": {"stress_ratio": 0.7, "maxload": 4000, "overload": 6000,
                                  "a_ol_applied": 12.0053, "a_dadn_max": 12.1116, "a_dadn_min": 12.4773,
                                  "a_kc_max": -1, "a_alpha_end": 18.3370, "a_retard_end": -1},
        "yang-baoban_Lu-420-17": {"stress_ratio": 0.7, "maxload": 3400, "overload": 5100,
                                  "a_ol_applied": 12.0038, "a_dadn_max": 12.0947, "a_dadn_min": 12.4599,
                                  "a_kc_max": -1, "a_alpha_end": 18.3370, "a_retard_end": 15.8853},
        "yang-baoban_Lu-420-09": {"stress_ratio": 0.7, "maxload": 4000, "overload": 8000,
                                  "a_ol_applied": 12.0159, "a_dadn_max": 12.3526, "a_dadn_min": 13.1055,
                                  "a_kc_max": -1, "a_alpha_end": 18.3370, "a_retard_end": -1}
    }
    try:
        test = test_information[specimen]
    except KeyError:
        print("Test Information was not collected in the database. Returning empty dict.")
        test = {}
    return test
