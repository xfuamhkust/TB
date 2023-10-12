import numpy as np


def ReadKPoint(FileName):
    f = [line.strip() for line in open(FileName,"r").readlines()]
    Method = f[0] # "Volume" or "Line"
    if Method == "Volume":
        fs = f[1].split()
        nk = [int(fs[i]) for i in range(len(fs))]
        if len(fs) == 2:
            nk.append(1)
        nk = np.array(nk)
        ParaKptIn = {"Method":       Method,
                     "nk":           nk,
                     }
    elif Method == "Line":
        nk = int(f[1])
        Lv = np.array([f[i].split() for i in range(2,5)],float)
        KptLvEnd = []
        KptLvEndName = []
        for fi in f[5:]:
            fs = fi.split()
            KptLvEnd.append([float(fs[i]) for i in range(3)])
            if len(fs) == 4:
                KptLvEndName.append(fs[3])
            else:
                KptLvEndName.append("")
        KptLvEnd      = np.array(KptLvEnd)
        KptLvEndName  = np.array(KptLvEndName)
        ParaKptIn = {"Method":       Method,
                     "nk":           nk, #number of points
                     "Lv":           Lv, # basis for primitive cell
                     "KptLvEnd":     KptLvEnd, # coordinates of special points in BZ under reciprocal basis
                     "KptLvEndName": KptLvEndName, #name of special points along the path in BZ
                     }
    elif Method == "Point":
        KptLv = []
        for fi in f[1:]:
            fs = fi.split()
            KptLv.append([float(fs[i]) for i in range(3)])
        KptLv = np.array(KptLv)
        ParaKptIn = {"Method":       Method,
                     "KptLv":        KptLv,
                     }
        
    return ParaKptIn