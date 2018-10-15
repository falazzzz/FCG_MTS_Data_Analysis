# 结果输出函数库
# Last Update: 2018/5/25 Version 2.3.0
# Lu Yunchao

import pandas as pd
import numpy as np
import time


def timestring():
    # usage：输出当前时间，以YYMMDD_HHMMSS的格式
    time_string = time.strftime('%Y%m%d_%H%M%S', time.localtime(time.time()))
    return time_string


def ArrangeData(dadn, cycles, dk, a=[], option="dk", step=1, save=True, name="temp"):
    if option == "dk":
        start = dk[0]
    elif option == "a":
        start = a[0]
    else:
        print("WRONG Option Input: Must be dk or a.")
        return False
    dadn_arranged = []
    cycles_arranged = []
    dk_arranged = []
    if len(a):
        a_arranged = [a[0]]
    dadn_arranged.append(dadn[0])
    cycles_arranged.append(cycles[0])
    dk_arranged.append(dk[0])
    seq = 0
    while seq < len(dk):
        if seq == len(dk)-1:
            dadn_arranged.append(dadn[seq])
            cycles_arranged.append(cycles[seq])
            dk_arranged.append(dk[seq])
            if len(a):
                a_arranged.append(a[seq])
            break
        if option == "dk" and dk[seq] >= start + step:
            dadn_arranged.append(dadn[seq])
            cycles_arranged.append(cycles[seq])
            dk_arranged.append(dk[seq])
            if len(a):
                a_arranged.append(a[seq])
            start = start + step
        elif option == "a" and a[seq] >= start + step:
            dadn_arranged.append(dadn[seq])
            cycles_arranged.append(cycles[seq])
            dk_arranged.append(dk[seq])
            if len(a):
                a_arranged.append(a[seq])
            start = start + step
        seq = seq + 1
    dadn_arranged = pd.Series(dadn_arranged, name="dadn(mm/cycle)")
    cycles_arranged = pd.Series(cycles_arranged, name="cycle")
    dk_arranged = pd.Series(dk_arranged, name="dK(Mpa.m0.5)")
    if len(a):
        a_arranged = pd.Series(a_arranged, name="cracklength(mm)")
        if save:
            con = pd.concat([dadn_arranged, cycles_arranged, dk_arranged, a_arranged], axis=1, )
            con.to_csv("FCG_ManualData_"+name+".csv", index=False, sep=",")
        return pd.DataFrame({"dadN": dadn_arranged, "cycles": cycles_arranged,
                             "dk": dk_arranged, "cracklength": a_arranged})
    else:
        if save:
            con = pd.concat([dadn_arranged, cycles_arranged, dk_arranged], axis=1)
            con.to_csv("FCG_MTSData_"+name+".csv", index=False, sep=",")
        return pd.DataFrame({"dadN": dadn_arranged, "cycles": cycles_arranged, "dk": dk_arranged})


def SaveData(dataset, name, filename="_"):
    # dataset和name需为numpy数组格式，dataset的数据和name位置一一对应
    # Note:对于长度不一的数据，最后部分会被填零

    # 寻找数组长度最大值
    maxlength = 0
    for data in dataset:
        if maxlength < data.shape[0]:
            maxlength = data.shape[0]

    # 数组方阵化
    datamatrix = []
    for data in dataset:
        fixed_data = []
        if data.shape[0] < maxlength:
            fixlength = maxlength - data.shape[0]
            fixlist = np.zeros(fixlength)
            fixed_data = np.concatenate((data, fixlist))
        else:
            fixed_data = data
        datamatrix.append(fixed_data)

    # 数据格式转化
    datamatrix = np.array(datamatrix)
    datamatrix = pd.DataFrame(data=datamatrix.T, columns=name)
    if filename == "_":
        time_string = timestring()
        filename = "temp_"+time_string


    datamatrix.to_csv(filename+".csv", index=False, sep=',')
    print("Data Saved to File:"+filename+".csv")
    return 1
