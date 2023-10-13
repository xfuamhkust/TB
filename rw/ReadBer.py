import numpy as np


def ReadBerryPara(FileName):
    f = [line.strip().split() for line in open(FileName,"r").readlines()]
    nk0 = int(f[0][0])
    Tem = np.array([float(f[1][i]) for i in range(len(f[1]))])
    EF  = np.array([float(f[2][i]) for i in range(len(f[2]))])
    ParaBerIn = {"nk0":    nk0,
                 "Tem":    Tem,
                 "EF":     EF,
                 }
    return ParaBerIn