import numpy as np

def FindNeighbor(ParaIn):
    
    # Parameters
    Nbr  = ParaIn["NeighborNumber"]
    Lv   = ParaIn["LatticeVector"]#primitive cell vector
    AtLv = ParaIn["AtomSite"]
    originBilbao = ParaIn["origin Bilbao"]
    AtLv = [vec - originBilbao for vec in AtLv]  # shift origin
    N1, N2, N3 = Nbr
    NumLv = (2*N1+1)*(2*N2+1)*(2*N3+1)#cell number
    NumAt = len(AtLv)
    
    # Find all atoms in [[-N1, N1],[-N2, N2],[-N3, N3]]
    LvAt = np.zeros((NumLv*NumAt,4),int)#[unitCelln1,unitCelln2,unitCelln3,atom]
    count = 0
    for n1 in range(-N1,N1+1):
        for n2 in range(-N2,N2+1):
            for n3 in range(-N3,N3+1):
                for jAt in range(1,NumAt+1):
                    LvAt[count] = [n1,n2,n3,jAt]
                    count += 1
    
    # Find all neighboring pairs between atoms in [0,0,0] and atoms in [[-N1, N1],[-N2, N2],[-N3, N3]]
    LvAtAt = np.zeros((NumLv*NumAt*NumAt,5),int)
    LvAtAtInd = np.zeros((NumLv*NumAt*NumAt,2),int)
    count = 0
    for iLvAt in range(NumLv*NumAt):
        n1,n2,n3,jAt = LvAt[iLvAt]
        for iAt in range(1,NumAt+1):
            LvAtAt[count] = [n1,n2,n3,iAt,jAt]#iAt: 0,0,0 cell
            LvAtAtInd[count] = [iAt,iLvAt]
            count += 1
    
    # Calculate the position of atoms in [0,0,0] cell
    AtXyz = AtLv @ Lv

    # Calculate the position of atoms in [[-N1, N1],[-N2, N2],[-N3, N3]]
    LvAtXyz = np.zeros((NumLv*NumAt,3))
    for iLvAt in range(NumLv*NumAt):
        n1,n2,n3,jAt = LvAt[iLvAt]
        LvAtXyz[iLvAt] = n1*Lv[0] + n2*Lv[1] + n3*Lv[2] + AtLv[jAt-1] @ Lv
    
    # Calculate the distance of the neighboring pairs
    LvAtAtDis = np.zeros((NumLv*NumAt*NumAt))
    for iLvAtAt in range(NumLv*NumAt*NumAt):
        iAt,iLvAt = LvAtAtInd[iLvAtAt]
        LvAtAtDis[iLvAtAt] = np.linalg.norm(LvAtXyz[iLvAt] - AtXyz[iAt-1])
    
    # Find all unique distances
    error = 1e-6
    DisUnq = [LvAtAtDis[0]]
    for iLvAtAt in range(NumLv*NumAt*NumAt):
        dr = abs(np.array(DisUnq) - LvAtAtDis[iLvAtAt])
        DisUnq.append(LvAtAtDis[iLvAtAt]) if np.min(dr) > error else None
    DisUnq = np.sort(np.array(DisUnq))
            
    # Classify neighboring pairs by distance
    LvAtAtDisInd = np.zeros((NumLv*NumAt*NumAt),int)
    for iLvAtAt in range(NumLv*NumAt*NumAt):
        dr = abs(DisUnq - LvAtAtDis[iLvAtAt])
        LvAtAtDisInd[iLvAtAt] = np.where(dr == np.min(dr))[0][0]
    
    # Output
    ParaNbr = {"LvAt":         LvAt,
               "AtXyz":        AtXyz,
               "LvAtXyz":      LvAtXyz,
               "LvAtAt":       LvAtAt,
               "LvAtAtInd":    LvAtAtInd,
               "LvAtAtDis":    LvAtAtDis,
               "LvAtAtDisInd": LvAtAtDisInd,
               "DisUnq":       DisUnq,
               }
    
    return ParaNbr