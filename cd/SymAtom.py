import numpy as np

def FindAtomSymmetry(ParaIn,ParaSym,ParaNbr):
    
    # Parameters
    AtLv         = ParaIn["AtomSite"]
    originBilbao = ParaSym["origin Bilbao"]
    AtLv = [vec - originBilbao for vec in AtLv]  # shift origin
    SymLv        = ParaSym["SymLv"]# space group operators under primitive cell basis
    LvAt         = ParaNbr["LvAt"] #all atoms in [[-N1, N1],[-N2, N2],[-N3, N3]]: [unitCelln1,unitCelln2,unitCelln3,atom]
    LvAtAt       = ParaNbr["LvAtAt"]# 1-direction hopping relation table: one entry of LvAtAt is [n1,n2,n3,iAt,jAt], iAt: 0,0,0 cell; jAt: n1,n2,n3 cell.
    LvAtAtInd    = ParaNbr["LvAtAtInd"]# one entry of LvAtAtInd is [iAt,iLvAt], iAt: 0,0,0 cell; iLvAt: row number of LvAt, taking values from 0 to NumLv*NumAt-1
    # LvAtAtDis    = ParaNbr["LvAtAtDis"]
    LvAtAtDisInd = ParaNbr["LvAtAtDisInd"]#index of entry in LvAtAtDis in DisUnq
    DisUnq       = ParaNbr["DisUnq"]#unique distances
    NumAt        = len(AtLv)
    NumSym       = len(SymLv)
    NumLvAt      = len(LvAt)
    NumLvAtAt    = len(LvAtAt)
    
    ''' Influence of Symmetry Operations on Atoms '''
    
    # Symmetry operations acting on atoms in 000
    error = 1e-6
    SymAt0 = np.zeros((NumSym,NumAt,4),int) # atom iAt in 0,0,0 cell --> atom jAt in dn1, dn2, dn3 cell under matrix iSym
    for iSym in range(NumSym):
        for iAt in range(0,NumAt):
            SymAt0ii = SymLv[iSym,0:3,0:3] @ AtLv[iAt] + SymLv[iSym,0:3,3]
            flag = 1
            for jAt in range(0,NumAt):
                dLv = SymAt0ii - AtLv[jAt]
                if np.linalg.norm(dLv - np.round(dLv)) < error:
                    dn1, dn2, dn3 = np.round(dLv).astype(int)
                    SymAt0[iSym,iAt] = [dn1, dn2, dn3, jAt]
                    flag = 0; break
            if flag:
                print("Wrong Symmetry: %d"%(iSym+1) + "!")

    # Symmetry operations acting on atoms in n1 n2 n3
    #atom iAt in primitive cell n1,n2,n3 --> atom OiAt in primitive cell n1+dn1, n2+dn2, n3+dn3
    #iLvAt is the index for atom i in primitive cell n1,n2,n3
    SymLvAt = np.zeros((NumSym,NumLvAt,4),int)
    for iSym in range(NumSym):
        for iLvAt in range(NumLvAt):
            n1,   n2,  n3,  iAt = LvAt[iLvAt]
            dn1, dn2, dn3, OiAt = SymAt0[iSym,iAt]
            On1, On2, On3 = SymLv[iSym,0:3,0:3] @ np.array([n1, n2, n3])
            SymLvAt[iSym,iLvAt] = [On1+dn1, On2+dn2, On3+dn3, OiAt]
    
    # Symmetry operations acting on neihboring pairs between atoms in 000 and atoms in n1 n2 n3
    # NumSym: number of matrices
    # NumLvAtAt: total 1-direction hopping numbers
    # under the action of the space group matrices, i atom in 0,0,0 cell --> Oi atom in dn1,dn2,dn3 cell, i atom in n1,n2,n3 cell --> Oj atom in On1+dn1, On2+dn2,On3+dn3 cell
    # SymLvAtAt is the rotation of the 1-direction hopping table LvAtAt
    SymLvAtAt = np.zeros((NumSym,NumLvAtAt,5),int)
    SymLvAtAtInd = np.zeros((NumSym,NumLvAtAt),int)#index of hopping relation Oi atom in dn1,dn2,dn3 cell --> Oj atom in On1+dn1,On2+dn2,On3+dn3 cell in hopping relation table LvAtAt
    for iSym in range(NumSym):
        for iLvAtAt in range(NumLvAtAt):
            iAt,iLvAt = LvAtAtInd[iLvAtAt]
            dn1,       dn2,    dn3, OiAt = SymAt0[iSym,iAt]
            On1dn1, On2dn2, On3dn3, OjAt = SymLvAt[iSym,iLvAt]
            SymLvAtAt[iSym,iLvAtAt] = [On1dn1-dn1, On2dn2-dn2, On3dn3-dn3, OiAt, OjAt]
            SymLvAtAtInd[iSym,iLvAtAt] = FindIndex(SymLvAtAt[iSym,iLvAtAt],LvAtAt)
    
    ''' Exclusion of Redundant Terms '''
    
    # Exclude neihboring pairs outside the [[-N1, N1],[-N2, N2],[-N3, N3]] after symmetry operations
    InclTF = np.min(SymLvAtAtInd,axis=0)>=0
    NumLvAtAtIncl = np.sum(InclTF)
    SymLvAtAtInd = np.copy(SymLvAtAtInd[:,InclTF])
    # InclInd = np.where(InclTF)[0] = SymLvAtAtInd[0]
    # print(SymLvAtAtInd)
    ''' Classify Neighboring Pairs by Symmetry '''
    
    # Group neighboring pairs by all symmetries, and sort them by distance
    LvAtAtClasTemp = [np.unique(SymLvAtAtInd[:,iLvAtAt]) for iLvAtAt in range(NumLvAtAtIncl)]#unique indices in each column
    LvAtAtClasTemp0 = np.array([LvAtAtClasTemp[iLvAtAt][0] for iLvAtAt in range(NumLvAtAtIncl)])#0th index for each entry in LvAtAtClasTemp
    _, ClasIndTemp = np.unique(LvAtAtClasTemp0, return_index=1)# position of the 0th ocurrence of the same indices in LvAtAtClasTemp0
    DisTemp = DisUnq[LvAtAtDisInd[LvAtAtClasTemp0[ClasIndTemp]]]# the distances corresponding to hopping relation indices in LvAtAtClasTemp0
    ClasInd = ClasIndTemp[np.argsort(DisTemp)]# entries in ClasIndTemp sorted by distances
    LvAtAtClas = [LvAtAtClasTemp[ClasInd[i]] for i in range(len(ClasInd))]# each class sorted by ascending order of distance
    # NumClas = len(LvAtAtClas)
    # LvAtAtClas0 = np.array([LvAtAtClas[iLvAtAt][0] for iLvAtAt in range(NumClas)])
    # Dis = DisUnq[LvAtAtDisInd[LvAtAtClas0]]
    
    # Consider the effect of Hermitian Conjugate Operator (<iAt|H|jAt> & <jAt|H|iAt>)
    NumClas = len(LvAtAtClas) 
    ComInd = []
    for iClas in range(NumClas):
        n1, n2, n3, iAt, jAt = LvAtAt[LvAtAtClas[iClas][0]]
        iLvAtAtConj = FindIndex(np.array([-n1, -n2, -n3, jAt, iAt]),LvAtAt)
        for jClas in range(iClas+1,NumClas):
            if iLvAtAtConj in LvAtAtClas[jClas]:
                ComInd.append([iClas,jClas]); break
    ComInd = np.array(ComInd)
    if len(ComInd):
        LvAtAtClasNew = []
        for iClas in range(NumClas):
            if iClas in ComInd[:,0]:
                Ind1,Ind2 = ComInd[np.where(ComInd[:,0]==iClas)[0][0]]
                LvAtAtClasNew.append(np.r_[LvAtAtClas[Ind1],LvAtAtClas[Ind2]])
            elif 1 - (iClas in ComInd[:,1]):
                LvAtAtClasNew.append(LvAtAtClas[iClas])
        LvAtAtClas = LvAtAtClasNew
    
    # Output
    ParaSymAt = {"SymLvAtAtInd": SymLvAtAtInd,
                 "LvAtAtClas":   LvAtAtClas,
                 "SymAt0":       SymAt0, # Used in CheckH
                 }
    
    return ParaSymAt

def FindIndex(x,X,tol = 1e-3):
    dX = np.sum(abs(X - x), axis=1)
    min_dX = min(dX)
    if min_dX < tol:
        ind = np.where(dX == min_dX)[0][0]
        return ind
    else:
        return -1