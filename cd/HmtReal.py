import numpy as np

def GetHamiltonianReal(ParaRel,HopValClas):
    
    # Parameters
    HopIndClas = ParaRel["LvAtAtOO"]
    NumClas = len(HopIndClas)
    
    # Combine Hopping terms in different groups
    HopInd = []
    for iClas in range(NumClas):
        HopInd += HopIndClas[iClas].tolist()
    HopInd = np.array(HopInd)
    
    # Combine Hopping values in different groups
    HopVal = []
    for iClas in range(NumClas):
        HopVal += HopValClas[iClas].tolist()
    HopVal = np.array(HopVal)
    
    # Region of n1 n2 n3 At Orb
    n1_min = np.min(HopInd[:,0]); n1_max = np.max(HopInd[:,0])
    n2_min = np.min(HopInd[:,1]); n2_max = np.max(HopInd[:,1])
    n3_min = np.min(HopInd[:,2]); n3_max = np.max(HopInd[:,2])
    NumLv = (n1_max-n1_min+1) * (n2_max-n2_min+1) * (n3_max-n3_min+1)
    # NumAt = np.max(HopInd[:,3:5])
    NumAt=len(set(HopInd[:,3].tolist()+HopInd[:,4].tolist()))
    AtOrbInd = []
    for iAt in range(0,NumAt):
        Orbi = np.unique(HopInd[np.where(HopInd[:,3]==iAt)[0],5])
        for Orbii in Orbi:
            AtOrbInd.append([iAt,Orbii])
    AtOrbInd = np.array(AtOrbInd)# AtOrbInd is sorted first by the value of iAt, then by the value of Orbii
    NumStt = len(AtOrbInd)
    # print(AtOrbInd)
    
    # Get Hamiltonian in real space (HR)
    HmtRealSpace = np.zeros((NumLv*NumStt*NumStt,7))
    count = 0
    for n1 in range(n1_min,n1_max+1):
        for n2 in range(n2_min,n2_max+1):
            for n3 in range(n3_min,n3_max+1):
                for iStt in range(0,NumStt):
                    for jStt in range(0,NumStt):
                        iAt, iOrb = AtOrbInd[iStt];
                        jAt, jOrb = AtOrbInd[jStt];
                        HopIndi = np.array([n1, n2, n3, iAt, jAt, iOrb, jOrb])
                        IndHopi = FindIndex(HopIndi, HopInd)
                        if IndHopi == -1 :
                            HopVali = 0
                        else:
                            HopVali = HopVal[IndHopi]
                        HmtRealSpace[count] = [n1, n2, n3, iStt, jStt, HopVali, 0] # No Spin & Real Hopping Terms!!!
                        count += 1
    
    # Get Hamiltonian in real space (HmtMatLv & HmtLv)
    HmtMatLv = np.zeros((NumLv,NumStt,NumStt))
    HmtLv    = np.zeros((NumLv,3),int)
    count = 0
    for n1 in range(n1_min,n1_max+1):
        for n2 in range(n2_min,n2_max+1):
            for n3 in range(n3_min,n3_max+1):
                HmtLv[count] = [n1, n2, n3]
                for iStt in range(0,NumStt):
                    for jStt in range(0,NumStt):
                        iAt, iOrb = AtOrbInd[iStt];
                        jAt, jOrb = AtOrbInd[jStt];
                        HopIndi = np.array([n1, n2, n3, iAt, jAt, iOrb, jOrb])
                        IndHopi = FindIndex(HopIndi, HopInd)
                        if IndHopi == -1 : continue
                        HmtMatLv[count,iStt,jStt] = HopVal[IndHopi]
                count += 1
        
    # Output
    ParaHmtR = {"AtOrbInd":     AtOrbInd,
                "HmtRealSpace": HmtRealSpace,
                "HmtMatLv":     HmtMatLv,
                "HmtLv":        HmtLv,
                }
    
    return ParaHmtR
    
def FindIndex(x,X,tol = 1e-3):
    dX = np.sum(abs(X - x), axis=1)
    min_dX = min(dX)
    if min_dX < tol:
        ind = np.where(dX == min_dX)[0][0]
        return ind
    else:
        return -1
    