import numpy as np

def ReadQcPara(FileName):
    f = [line.strip().split() for line in open(FileName,"r").readlines()]
    Method = f[0][0]
    Emin, Emax, nE, eta = float(f[1][0]), float(f[1][1]), int(f[1][2]), float(f[1][3]),
    Cell = np.array([int(f[2][0]), int(f[2][1])])
    Dir  = int(f[2][2])
    Wid, Len = int(f[3][0]), int(f[3][1])
    ParaQcIn  = {"Method": Method,
                 "nE":     nE,
                 "Emin":   Emin,
                 "Emax":   Emax,
                 "eta":    eta,
                 "Cell":   Cell,
                 "Dir":    Dir,
                 "Wid":    Wid,
                 "Len":    Len,
                 }

    return ParaQcIn