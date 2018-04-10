# 结果输出函数库

import pandas as pd


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
