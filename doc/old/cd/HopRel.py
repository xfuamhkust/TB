import numpy as np
import sympy as sp

def FindRelation(ParaIn,ParaSym,ParaNbr,ParaSymAt):
    
    # Parameters
    AtNum        = ParaIn["AtomNumber"]
    AtOrb        = ParaIn["AtomOrbital"]
    AtTypeInd    = ParaIn["AtomTypeIndex"]
    LvAtAt       = ParaNbr["LvAtAt"]
    SymOrb       = ParaSym["SymOrb"]
    SymLvAtAtInd = ParaSymAt["SymLvAtAtInd"]
    LvAtAtClas   = ParaSymAt["LvAtAtClas"]
    NumSym       = len(SymOrb[0])
    NumAtType    = len(AtOrb)
    NumClas      = len(LvAtAtClas)
    
    ''' Get the Representations of Input Orbitals '''
    
    # Write Symmetries acting on SPDF together
    SymSPDF = np.zeros((NumSym,16,16))
    for iOrb in range(4):
        SymSPDF[:,iOrb**2:(iOrb+1)**2,iOrb**2:(iOrb+1)**2] = SymOrb[iOrb]
    
    # Write Symmetries acting on input orbitals
    Ind  = [np.ix_(np.where(AtOrb[iAt])[0],np.where(AtOrb[iAt])[0]) for iAt in range(NumAtType)]
    SymAtOrb = []
    for iAt in range(NumAtType):
        SymAtOrbi = []
        for iSym in range(NumSym):
            SymAtOrbi.append(SymSPDF[iSym][Ind[iAt]])
        for i in range(AtNum[iAt]):
            SymAtOrb.append(np.array(SymAtOrbi))
    
    ''' Find Relations of Hopping Terms in Groups '''
    
    # Find all hopping terms including atoms and orbitals (no spin)
    LvAtAtOO = []
    for iClas in range(NumClas):
        LvAtAtClasi = LvAtAtClas[iClas]
        LvAtAtOOi = []
        NumHopi = len(LvAtAtClasi)
        for iHop in range(NumHopi):
            n1, n2, n3, iAt, jAt = LvAtAt[LvAtAtClasi[iHop]]
            Orbi = np.where(AtOrb[AtTypeInd[iAt-1]])[0]
            Orbj = np.where(AtOrb[AtTypeInd[jAt-1]])[0]
            for iOrb in Orbi:
                for jOrb in Orbj:
                    LvAtAtOOi.append([n1, n2, n3, iAt, jAt, iOrb, jOrb])
        LvAtAtOO.append(np.array(LvAtAtOOi,"int"))
    
    # Find relations between hopping terms
    NumHopClas = np.array([len(LvAtAtOO[iClas]) for iClas in range(NumClas)])
    HopRelClas = []
    for iClas in range(NumClas):
        NumHopi = NumHopClas[iClas]
        HopRelClasi = np.zeros((NumHopi*(NumSym+1),NumHopi))
        count = 0
        for iHop in range(NumHopi):
            n1, n2, n3, iAt, jAt, iOrb, jOrb = LvAtAtOO[iClas][iHop]
            LvAtAtii = np.array([n1, n2, n3, iAt, jAt])
            LvAtAtIndii = FindIndex(LvAtAtii,LvAtAt)
            Orbi = np.where(AtOrb[AtTypeInd[iAt-1]])[0]
            Orbj = np.where(AtOrb[AtTypeInd[jAt-1]])[0]
            iOrbInd = np.where(Orbi == iOrb)[0][0]
            jOrbInd = np.where(Orbj == jOrb)[0][0]
            # Effect of symmetry
            SymLvAtAtIndii = SymLvAtAtInd[:,np.where(SymLvAtAtInd[0] == LvAtAtIndii)[0][0]]
            for iSym in range(NumSym):
                SymLvAtAtIndiii = SymLvAtAtIndii[iSym]
                SymLvAtAtiii = LvAtAt[SymLvAtAtIndiii]
                n1_, n2_, n3_, iAt_, jAt_ = SymLvAtAtiii
                SymiOrb = SymAtOrb[iAt-1][iSym,:,iOrbInd]
                SymjOrb = SymAtOrb[jAt-1][iSym,:,jOrbInd]
                Orbi_ = np.where(AtOrb[AtTypeInd[iAt_-1]])[0]
                Orbj_ = np.where(AtOrb[AtTypeInd[jAt_-1]])[0]
                HopRelClasiii = np.zeros(NumHopi)
                HopRelClasiii[iHop] = -1
                for io in range(len(SymiOrb)):
                    for jo in range(len(SymjOrb)):
                        SymLvAtAtOOiii = np.array([n1_, n2_, n3_, iAt_, jAt_, Orbi_[io], Orbj_[jo]])
                        SymLvAtAtOOIndiii = FindIndex(SymLvAtAtOOiii,LvAtAtOO[iClas])
                        HopRelClasiii[SymLvAtAtOOIndiii] += SymiOrb[io] * SymjOrb[jo] # No Spin & Real Hopping Terms!!!
                HopRelClasi[count] += HopRelClasiii
                count += 1
            # Effect of h.c.
            n1_, n2_, n3_, iAt_, jAt_, iOrb_, jOrb_ = -n1, -n2, -n3, jAt, iAt, jOrb, iOrb
            ConjLvAtAtOOiii = np.array([n1_, n2_, n3_, iAt_, jAt_, iOrb_, jOrb_])
            ConjLvAtAtOOIndiii = FindIndex(ConjLvAtAtOOiii,LvAtAtOO[iClas])
            HopRelClasiii = np.zeros(NumHopi)
            HopRelClasiii[iHop] = -1
            HopRelClasiii[ConjLvAtAtOOIndiii] += 1 # No Spin & Real Hopping Terms!!!
            HopRelClasi[count] += HopRelClasiii
            count += 1
        HopRelClas.append(HopRelClasi)
        
    # Simplification of hopping relations
    # 1. Delete zeros lines
    for iClas in range(NumClas):
        HopRelClasi = HopRelClas[iClas]
        SumAbsHRi = np.sum(abs(HopRelClasi),axis=1)
        IndNonZeroi = np.where(SumAbsHRi)[0]
        HopRelClasNewi = HopRelClasi[IndNonZeroi]
        HopRelClas[iClas] = HopRelClasNewi
    # 2. Delete directly repeated lines
    for iClas in range(NumClas):
        HopRelClasi = HopRelClas[iClas]
        HopRelClasNewi = np.unique(HopRelClasi,axis=0)
        HopRelClas[iClas] = HopRelClasNewi
    # 3. Reduced row echelon form
    error = 1e-6
    for iClas in range(NumClas):
        HopRelClasi = HopRelClas[iClas]
        HopRelClasNewi = Rref(HopRelClasi,error)
        HopRelClas[iClas] = HopRelClasNewi
    
    # Alternative form of hopping relations
    HopRelAltClas = []
    for iClas in range(NumClas):
        NumHopi = NumHopClas[iClas]
        HopRelClasi = HopRelClas[iClas]
        FreIndi = [] # Indices of free terms
        UnfIndi = [] # Indices of unfree terms
        UnfVali = [] # Values of unfree terms. [[Ind1,Val1],...]
        for HopRelClasij in HopRelClasi:
            IndNonZero = np.where(HopRelClasij)[0]
            UnfIndi.append(IndNonZero[0])
            UnfVali.append([])
            for INZi in IndNonZero[1:]:
                UnfVali[-1].append([INZi,-HopRelClasij[INZi]])
        for iHop in range(NumHopi):
            if 1-(iHop in UnfIndi): FreIndi.append(iHop)
        HopRelAltClas.append([FreIndi,UnfIndi,UnfVali])
    
    # Output
    ParaRel = {"LvAtAtOO":      LvAtAtOO,
               "HopRelAltClas": HopRelAltClas,
               }
    
    return ParaRel

def FindIndex(x,X,tol = 1e-3):
    dX = np.sum(abs(X - x), axis=1)
    min_dX = min(dX)
    if min_dX < tol:
        ind = np.where(dX == min_dX)[0][0]
        return ind
    else:
        return -1

def Rref(M,error):
    m,n = M.shape
    ErrInt = round(-np.log10(error)); error = 10**-ErrInt
    # Rref
    M1 = np.array(sp.Matrix(M).rref()[0].tolist(),float)
    # Eliminate error
    for i in range(m):
        for j in range(n):
            Mij = M1[i,j]; rMij = np.round(Mij,ErrInt)
            if abs(Mij - rMij) < error*1e-2:
                M1[i,j] = rMij
    # Delete zero lines
    if len(M1):
        In0 = np.where(np.sum(abs(M1),axis=1))[0]
        return M1[In0]
        
    return M1