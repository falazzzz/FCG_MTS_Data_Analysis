import numpy as np
import matplotlib.pyplot as plt
from FCGAnalysisLib import read_data


def ClosureloadFromOffset(offset, n_list, compliance_load, criteria0):
    for seq4, value in enumerate(offset):
        if value > criteria0:
            # 若是第2个~第n-1个offset，连续两个超过判断标准criteria0的，则认定其之前为closure load
            if (1 < seq4 < len(offset)-1) and (offset[seq4+1] > criteria0):
                # 采用插值法确定criteria0对应的载荷
                closeload = ((criteria0 - offset[seq4])/(offset[seq4-1] - offset[seq4]))*(compliance_load[seq4-1] - compliance_load[seq4]) + compliance_load[seq4]
                close_cycle = n_list[seq, 1]
                break
            # 若是第n个offset，其本身超过标准criteria0的，则认定其之前为closure load
            elif seq4 == len(offset)-1:
                # 采用插值法确定criteria0对应的载荷
                closeload = ((criteria0 - offset[seq4]) / (offset[seq4 - 1] - offset[seq4])) * (compliance_load[seq4 - 1] - compliance_load[seq4]) + compliance_load[seq4]
                close_cycle = n_list[seq, 1]
                break
    return closeload, close_cycle


plot = 1

sequence = ["yang-baoban_Lu-420-06"]
stress_ratio = 0.1
specimen, width, notch_length, thickness, elastic_modulus, yield_strength, precrack_length = \
    read_data.ReadTestInf(sequence=sequence[0])
print('specimen name:', str(specimen))

unloadingrange = [0.70, 0.95]
compliance_range = 0.1
compliance_step = 0.05
criteria0 = 0.02    # 判断标准，ASTM方法
criteria1 = 0.08    # 判断标准，nASTM方法
cycles_r, load_r, cod_r = read_data.ReadCodData(sequence[0])
n_list, n_seq = np.unique(cycles_r, return_index=True)
n_couple_list = n_list.reshape(-1, 2)
n_couple_seq = n_seq.reshape(-1, 2)
length = n_couple_list.shape[0]
closureload = []
closureload_cycle = []
closureload_r = []
closureload_cycle_r = []
for seq in range(length):
    start = n_couple_seq[seq, 0]
    mid = n_couple_seq[seq, 1] - 1     # 两个循环为一组，start为前循环第一位，mid为前循环最后一位，end为后循环最后一位
    if seq == length - 1:
        end = len(cycles_r) - 1
    else:
        end = n_couple_seq[seq+1, 0] - 1
    load_cycle0 = load_r[start:mid+1]
    cod_cycle0 = cod_r[start:mid+1]
    load_cycle1 = load_r[mid+1:end+1]
    cod_cycle1 = cod_r[mid+1:end+1]
    max_load_seq_cycle1 = np.argmax(load_cycle1)
    max_load_seq_cycle0 = np.argmax(load_cycle0)
    min_load_seq_cycle1 = np.argmin(load_cycle1)
    min_load_seq_cycle0 = np.argmin(load_cycle0)
    if min_load_seq_cycle0 < max_load_seq_cycle0:
        # 即连续的两个循环中第一个循环先出现谷值的情况,数据选择两个循环之间的卸载部分
        unloadingpart_load = np.concatenate((load_cycle0[max_load_seq_cycle0:], load_cycle1[:min_load_seq_cycle1]))
        unloadingpart_cod = np.concatenate((cod_cycle0[max_load_seq_cycle0:], cod_cycle1[:min_load_seq_cycle1]))
        loadingpart_load = load_cycle0[min_load_seq_cycle0: max_load_seq_cycle0]
        loadingpart_cod = cod_cycle0[min_load_seq_cycle0: max_load_seq_cycle0]
    elif max_load_seq_cycle0 < min_load_seq_cycle0:
        # 即连续的两个循环中第一个循环先出现峰值的情况，数据选择第一个循环的卸载部分
        unloadingpart_load = load_cycle0[max_load_seq_cycle0:min_load_seq_cycle0+1]
        unloadingpart_cod = cod_cycle0[max_load_seq_cycle0:min_load_seq_cycle0+1]
        loadingpart_load = np.concatenate((load_cycle0[min_load_seq_cycle0:], load_cycle1[:max_load_seq_cycle1]))
        loadingpart_cod = np.concatenate((cod_cycle0[min_load_seq_cycle0:], cod_cycle1[:max_load_seq_cycle1]))
    # 完成卸载段数据和加载段数据的筛选

    unload_load_valid = []
    unload_cod_valid = []
    upperlimit = unloadingrange[1] * np.max(unloadingpart_load)
    lowerlimit = unloadingrange[0] * np.max(unloadingpart_load)
    for seq1, value in enumerate(unloadingpart_load):
        if lowerlimit < value < upperlimit:
            unload_load_valid.append(unloadingpart_load[seq1])
            unload_cod_valid.append(unloadingpart_cod[seq1])
    unload_load_valid = np.array(unload_load_valid).reshape(1, -1)
    unload_load_valid = unload_load_valid[0]
    unload_cod_valid = np.array(unload_cod_valid).reshape(1, -1)
    unload_cod_valid = unload_cod_valid[0]
    # 筛选卸载曲线内拟合区间内的数据

    p = np.polyfit(unload_load_valid, unload_cod_valid, deg=1)
    compliance_op, _ = p
    # 拟合得到open crack compliance

    begin = unloadingrange[1]       # 指定compliance计算区间S1的开头值为计算unloadingrange的最高位置
    compliance_section = []
    while begin - compliance_range > 0:
        right = begin * np.max(loadingpart_load)
        left = (begin - compliance_range) * np.max(loadingpart_load)
        begin = begin - compliance_step
        compliance_section.append([left, right])
    compliance_section = np.array(compliance_section)
    # 计算出每一个要计算compliance的载荷区间

    compliance = []
    compliance_load = []
    for seq2 in range(compliance_section.shape[0]):
        fitting_load = []
        fitting_cod = []
        upperlimit = compliance_section[seq2, 1]
        lowerlimit = compliance_section[seq2, 0]
        for seq3, value in enumerate(loadingpart_load):
            if lowerlimit < value < upperlimit:
                fitting_load.append(loadingpart_load[seq3])
                fitting_cod.append(loadingpart_cod[seq3])
        fitting_load = np.array(fitting_load).reshape(1, -1)
        fitting_load = fitting_load[0]
        fitting_cod = np.array(fitting_cod).reshape((1, -1))
        fitting_cod = fitting_cod[0]
        p = np.polyfit(fitting_load, fitting_cod, deg=1)
        compliance1, _ = p
        load1 = np.mean(fitting_load)
        compliance.append(compliance1)
        compliance_load.append(load1)
    compliance = np.array(compliance)
    compliance_load = np.array(compliance_load)
    # 完成各段的compliance计算

    offset = (compliance_op - compliance) / compliance_op
    relativeoffset = offset / np.max(offset)
    # 计算偏移量offset

    closeload, close_cycle = ClosureloadFromOffset(offset=offset, n_list=n_couple_list, compliance_load=compliance_load, criteria0=criteria0)
    closureload.append(closeload)
    closureload_cycle.append(close_cycle)

    closeload, close_cycle = ClosureloadFromOffset(offset=relativeoffset, n_list=n_couple_list, compliance_load=compliance_load, criteria0=criteria1)
    closureload_r.append(closeload)
    closureload_cycle_r.append(close_cycle)


    if plot:
        if seq % 5 == 0 and 120000 < n_couple_list[seq, 0] < 511200:
            plt.figure(num=seq, figsize=(10, 8))
            plt.scatter(offset, compliance_load, lw=1, marker='*', label='Offset')
            plt.plot(offset, compliance_load)
            plt.scatter(relativeoffset, compliance_load, lw=1, marker='*', label='RelativeOffset')
            plt.plot(relativeoffset, compliance_load)
            plt.plot(np.full(2, 0.02), np.array([max(compliance_load), min(compliance_load)]))
            plt.plot(np.full(2, 0.08), np.array([max(compliance_load), min(compliance_load)]))
            plt.title("Cycles:" + str(n_couple_list[seq, 1]))
            plt.xlabel("OFFSET")
            plt.ylabel("LOAD/kN")
            plt.legend()
            plt.grid()

closureload = np.array(closureload)
closureload_cycle = np.array(closureload_cycle)

cycles_MTS, cracklength_MTS, _, _, _, _, kclosure, closureload_MTS = \
    read_data.ReadOriginResult(sequence=sequence[0], closure=True, cod=False)

if plot:
    plt.figure(num=seq+1, figsize=(10, 8))
    plt.scatter(closureload_cycle, closureload, lw=1, marker='*', label='Closureload_ASTM')
    plt.scatter(cycles_MTS, closureload_MTS, lw=1, marker='*', label='Closureload_MTS')
    plt.scatter(closureload_cycle_r, closureload_r, lw=1, marker='*', label='Closureload_nASTM')
    plt.title("Closure: " + sequence[0])
    plt.xlabel("Cycle")
    plt.ylabel("Closure Load/kN")
    plt.legend()
    plt.grid()

if plot:
    plt.show()

