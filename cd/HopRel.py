import numpy as np
import sympy as sp

def FindRelation(ParaIn,ParaSym,ParaNbr,ParaSymAt):
    
    # Parameters
    # AtNum        = ParaIn["AtomNumber"]
    AtOrb        = ParaIn["AtomOrbital"]
    AtTypeInd    = ParaIn["AtomTypeIndex"]
    LvAtAt       = ParaNbr["LvAtAt"]
    SymOrb       = ParaSym["SymOrb"]
    SymAtOrb = ParaSym["SymAtOrb"]
    SymLvAtAtInd = ParaSymAt["SymLvAtAtInd"]
    LvAtAtClas   = ParaSymAt["LvAtAtClas"]
    NumSym       = len(SymOrb[0])
    # NumAtType    = len(AtOrb)
    NumClas      = len(LvAtAtClas)
    
    ''' Find Relations of Hopping Terms in Groups '''
    
    # Find all hopping terms including atoms and orbitals (no spin)
    LvAtAtOO = []
    for iClas in range(NumClas):
        LvAtAtClasi = LvAtAtClas[iClas]
        LvAtAtOOi = []
        NumHopi = len(LvAtAtClasi)
        for iHop in range(NumHopi):
            n1, n2, n3, iAt, jAt = LvAtAt[LvAtAtClasi[iHop]]
            Orbi = np.where(AtOrb[AtTypeInd[iAt]])[0]
            Orbj = np.where(AtOrb[AtTypeInd[jAt]])[0]
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
            Orbi = np.where(AtOrb[AtTypeInd[iAt]])[0]
            Orbj = np.where(AtOrb[AtTypeInd[jAt]])[0]
            iOrbInd = np.where(Orbi == iOrb)[0][0]
            jOrbInd = np.where(Orbj == jOrb)[0][0]
            # Effect of symmetry
            SymLvAtAtIndii = SymLvAtAtInd[:,np.where(SymLvAtAtInd[0] == LvAtAtIndii)[0][0]]
            for iSym in range(NumSym):
                SymLvAtAtIndiii = SymLvAtAtIndii[iSym]
                SymLvAtAtiii = LvAtAt[SymLvAtAtIndiii]
                n1_, n2_, n3_, iAt_, jAt_ = SymLvAtAtiii
                SymiOrb = SymAtOrb[iAt][iSym,:,iOrbInd]
                SymjOrb = SymAtOrb[jAt][iSym,:,jOrbInd]
                Orbi_ = np.where(AtOrb[AtTypeInd[iAt_]])[0]
                Orbj_ = np.where(AtOrb[AtTypeInd[jAt_]])[0]
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
        
    # Reduced row echelon form
    error = 1e-6
    for iClas in range(NumClas):
        HopRelClasi = HopRelClas[iClas]
        HopRelClasNewi = Rref(HopRelClasi,error)
        In0 = np.where(np.sum(abs(HopRelClasNewi),axis=1)>error)[0]
        HopRelClas[iClas] = HopRelClasNewi[In0]
    
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

def Rref(A,tol):
    # This code is changed from that on the following website: 
    # https://stackoverflow.com/questions/7664246/python-built-in-function-to-do-matrix-reduction
    # The original code is created by stackoverflow[user:1887559]
    m,n = A.shape
    pcol = -1 #pivot colum
    for i in range(m):
        pcol += 1
        if pcol >= n : break
        #pivot index
        pid = np.argmax( abs(A[i:,pcol]) )
        #Row exchange
        A[i,:],A[pid+i,:] = A[pid+i,:].copy(),A[i,:].copy()
        #pivot with given precision
        while pcol < n and abs(A[i,pcol]) < tol:
            pcol += 1
            if pcol >= n : break
            #pivot index
            pid = np.argmax( abs(A[i:,pcol]) )
            #Row exchange
            A[i,:],A[pid+i,:] = A[pid+i,:].copy(),A[i,:].copy()
        if pcol >= n : break
        pivot = float(A[i,pcol])
        for j in range(m):
            if j == i: continue
            mul = float(A[j,pcol])/pivot
            A[j,:] = A[j,:] - A[i,:]*mul
        A[i,:] /= pivot
    # Transfer all integer-like numbers to integers
    for i in range(m):
        for j in range(n):
            if np.abs(A[i,j]-np.round(A[i,j])) < tol:
                A[i,j] = np.round(A[i,j])
    return A