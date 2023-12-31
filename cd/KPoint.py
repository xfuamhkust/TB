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
        nk           = ParaKptIn["nk"] # number of total points after interpolation along the path in BZ
        Lv           = ParaKptIn["Lv"] # basis for primitive cell
        KptLvEnd     = ParaKptIn["KptLvEnd"] # special points along the path in  BZ under basis RecLv
        KptLvEndName = ParaKptIn["KptLvEndName"]# names of the special points along the path in  BZ
        Kpt0R,KptR,KptLv,KptXyz = GetKptPath(Lv,KptLvEnd,nk)
        ParaKpt = {"Method":      Method,
                   "Kpt0Name":    KptLvEndName,
                   "Kpt0R":       Kpt0R, #cumsum of distances  between special points  starting from KptPath0Xyz[0]
                   "KptR":        KptR, # cumulative distances along the path in BZ
                   "KptLv":       KptLv, #interpolated points along the path under RecLv basis in BZ
                   "KptXyz":      KptXyz, #  interpolated points along the path under Cartesian basis in BZ
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
    """

    :param Lv: Basis of primitive cell
    :param KptPath0Lv: # special points along the path in  BZ under reciprocal basis
    :param NumKpt: total number of points along the path in BZ after interpolation
    :return:
    """

    RecLv = GetRecLv(Lv) # basis in BZ  (TODO: check if it is primitive in reciprocal space)
    KptPath0Xyz = KptPath0Lv @ RecLv # Cartesian coordinates for special points in BZ
    KptDis0 = np.linalg.norm(KptPath0Xyz[1:] - KptPath0Xyz[:-1],axis=1)# distance between neighboring special points in KptPath0Xyz
    NumPath = len(KptDis0) # number of distances along the path
    KptPath0R = np.zeros(NumPath+1)#cumsum of distances starting between special points  starting from KptPath0Xyz[0]
    for iPath in range(NumPath):
        KptPath0R[iPath+1] = KptPath0R[iPath] + KptDis0[iPath]
    NumKptPath = np.round(NumKpt * KptDis0 / np.sum(KptDis0)).astype(int)# number of points along each segment, according to length/total distance
    NumKpt = np.sum(NumKptPath) # number of total points
    KptPathR = np.zeros(NumKpt+1)# cumulative distances along the path in BZ
    KptPathLv = np.zeros((NumKpt+1,3))# interpolated points along the path under RecLv basis in BZ
    KptPathXyz = np.zeros((NumKpt+1,3))  # interpolated points along the path under Cartesian basis in BZ
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
    KptPathLv[-1]  = KptPath0Lv[-1] # make a closed path
    KptPathXyz[-1] = KptPath0Xyz[-1]# make a closed path
    PathTemp += np.linalg.norm(KptPathXyz[count] - XyzTemp)
    KptPathR[count] = PathTemp
    return KptPath0R,KptPathR,KptPathLv,KptPathXyz
    # KptPath0R: cumsum of distances starting between special points  starting from KptPath0Xyz[0]
    # KptPathR: cumulative distances along the path in BZ
    # KptPathLv: interpolated points along the path under RecLv basis in BZ
    # KptPathXyz: interpolated points along the path under Cartesian basis in BZ