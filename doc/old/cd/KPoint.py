import numpy as np
pi2 = 2*np.pi

def GetKPoint(ParaKptIn):
    
    Method           = ParaKptIn["Method"]
    if Method == "Volume":
        nk           = ParaKptIn["nk"]
        KptLv = GetKptUni(nk)
        ParaKpt = {"Method":      Method,
                   "KptLv":       KptLv,
                   }
    elif Method == "Line":
        nk           = ParaKptIn["nk"]
        Lv           = ParaKptIn["Lv"]
        KptLvEnd     = ParaKptIn["KptLvEnd"]
        KptLvEndName = ParaKptIn["KptLvEndName"]
        Kpt0R,KptR,KptLv,KptXyz = GetKptPath(Lv,KptLvEnd,nk)
        ParaKpt = {"Method":      Method,
                   "Kpt0Name":    KptLvEndName,
                   "Kpt0R":       Kpt0R,
                   "KptR":        KptR,
                   "KptLv":       KptLv,
                   "KptXyz":      KptXyz,
                   }
    elif Method == "Point":
        KptLv        = ParaKptIn["KptLv"]
        ParaKpt = {"Method":      Method,
                   "KptLv":       KptLv,
                   }
    return ParaKpt    
        
def GetKptUni(nk=[10,10,1]):
    nk1, nk2, nk3 = nk
    Nk = nk1*nk2*nk3
    c1 = np.linspace(1/nk1-1,1-1/nk1,nk1)*0.5
    c2 = np.linspace(1/nk2-1,1-1/nk2,nk2)*0.5
    c3 = np.linspace(1/nk3-1,1-1/nk3,nk3)*0.5
    k = np.zeros((Nk,3))
    count = 0
    for ik1 in range(nk1):
        for ik2 in range(nk2):
            for ik3 in range(nk3):
                k[count] = [c1[ik1], c2[ik2], c3[ik3]]
                count += 1
    return k

def GetKptHiSym():
    pass

def GetRecLv(Lv):
    a1,a2,a3 = Lv
    Vol = np.dot(a1,np.cross(a2,a3))
    b1 = pi2 * np.cross(a2,a3) / Vol
    b2 = pi2 * np.cross(a3,a1) / Vol
    b3 = pi2 * np.cross(a1,a2) / Vol
    RecLv = np.array([b1,b2,b3])
    return RecLv

def GetKptPath(Lv,KptPath0Lv,NumKpt=100):
    RecLv = GetRecLv(Lv)
    KptPath0Xyz = KptPath0Lv @ RecLv
    KptDis0 = np.linalg.norm(KptPath0Xyz[1:] - KptPath0Xyz[:-1],axis=1)
    NumPath = len(KptDis0)
    KptPath0R = np.zeros(NumPath+1)
    for iPath in range(NumPath):
        KptPath0R[iPath+1] = KptPath0R[iPath] + KptDis0[iPath]
    NumKptPath = np.round(NumKpt * KptDis0 / np.sum(KptDis0)).astype(int)
    NumKpt = np.sum(NumKptPath)
    KptPathR = np.zeros(NumKpt+1)
    KptPathLv = np.zeros((NumKpt+1,3))
    KptPathXyz = np.zeros((NumKpt+1,3))
    count = 0
    XyzTemp = KptPathXyz[0]
    PathTemp = 0
    for iPath in range(NumPath):
        KptRate = np.linspace(0,1-1/NumKptPath[iPath],NumKptPath[iPath])
        for iKpt in range(NumKptPath[iPath]):
            KptPathLv[count]  = (1-KptRate[iKpt]) * KptPath0Lv[iPath]  + KptRate[iKpt] * KptPath0Lv[iPath+1]
            KptPathXyz[count] = (1-KptRate[iKpt]) * KptPath0Xyz[iPath] + KptRate[iKpt] * KptPath0Xyz[iPath+1]
            PathTemp += np.linalg.norm(KptPathXyz[count] - XyzTemp)
            KptPathR[count] = PathTemp
            XyzTemp = KptPathXyz[count]
            count += 1
    KptPathLv[-1]  = KptPath0Lv[-1]
    KptPathXyz[-1] = KptPath0Xyz[-1]
    PathTemp += np.linalg.norm(KptPathXyz[count] - XyzTemp)
    KptPathR[count] = PathTemp
    return KptPath0R,KptPathR,KptPathLv,KptPathXyz