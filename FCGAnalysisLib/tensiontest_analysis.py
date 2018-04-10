# 简单拉伸实验拟合Ramberg-Osgood模型函数库
# Last Update: 2018/4/10 Version 2.1.4
# Lu Yunchao
import numpy as np


def TrueStress(engineeringstress, engineeringstrain):
    # usage：由工程应力和工程应变计算真实应力
    truestress = engineeringstress * (1 + engineeringstrain)
    return truestress


def TrueStrain(engineeringstrain):
    # usage：由工程应变计算真实应变
    truestrain = np.log(1 + engineeringstrain)
    return truestrain


def RambergOsgoodFitting(s, e, ys, E):
    # usage: 由真实应力和真实应变拟合Ramberg-Osgood本构模型
    # input parameter:
    # s: 真实应力/MPa
    # e：真实应变
    # ys：屈服强度/GPa
    # E：弹性模量/GPa
    # return parameter:
    # n: Ramberg-Osgood模型参数（硬化指数）
    # alpha: Ramberg-Osgood模型参数
    ys = ys * 1e3
    E = E * 1e3
    Y0 = np.array([(E/ys)*strain for strain in e])
    Y1 = np.array([stress/ys for stress in s])
    Y = np.log(Y0 + Y1)
    del Y0, Y1
    X = np.log([stress/ys for stress in s])
    p = np.polyfit(X, Y, deg=1)
    n, B = p
    alpha = np.exp(B)
    return n, alpha


def RambergOsgoodCalculation(n, alpha, s, ys, E):
    # usage: 由已经得到的Ramberg-Osgood参数和需要拟合的真实应力值，以该本构关系计算真实应变
    # input parameter:
    # n，alpha：Ramberg-Osgood模型参数
    # s：真实应力/MPa
    # ys：屈服强度/GPa
    # E：弹性模量/GPa
    # return parameter:
    # e：由Ramberg-Osgood本构关系计算得到的真实应变
    e = []
    ys = ys * 1e3
    E = E * 1e3
    for stress in s:
        left = stress/ys + alpha*(stress/ys)**n
        e.append(left*ys/E)
    return np.array(e)

