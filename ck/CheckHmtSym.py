import numpy as np
from cd.HmtK import HamiltonianK
pi2 = np.pi * 2

def CheckHamiltonianSymmetry(ParaIn,ParaSym,ParaSymAt,ParaHmtR):
    
    '''
    This code is to check the symmetry of Hamiltonian in k-space.
    Ref: PRM 2, 103805 (2018)
        H(k) = D^k(g^−1) H(gk) D^k(g) 
    If we choose symmetric gauge, D^k(g) = e^i(gk)·t D(g)
        H(k) = D(g^−1) H(gk) D(g)
    Hamiltonian from periodic gauge to symmetric gauge:
        PrdGag: Hij(k) = Hij[R] e^ik·(R)
        SymGag: Hij(k) = Hij[R] e^ik·(R+tj−ti)
        H^S(k) = Fct(k)^+ H^P(k) Fct(k)
        Fct(k)ij = e^ik·ti if i == j else 0
    Symmetry operation matrix D(g)
        Only atom and orbital freedoms are considered (no spin).
        Atom part: D^at_ij = delta_i,gj
        Orbital part: D^orb_ab = <a|gb>
    '''
    
    '''################   Parameters   ################'''
    # Parameters
    Lv            = ParaIn["LatticeVector"]
    SymXyz        = ParaSym["SymXyz"]
    HmtMatLv      = ParaHmtR["HmtMatLv"]
    HmtLv         = ParaHmtR["HmtLv"]
    AtLv          = ParaIn["AtomSite"]
    SymAt0        = ParaSymAt["SymAt0"]
    AtNum         = ParaIn["AtomNumber"]
    AtOrb         = ParaIn["AtomOrbital"]
    AtTypeInd     = ParaIn["AtomTypeIndex"]
    SymOrb        = ParaSym["SymOrb"]
    NumSym        = len(SymXyz)
    NumStt        = len(HmtMatLv[0])
    NumAt         = len(SymAt0[0])
    NumAtType     = len(AtOrb)
    AtOrbAll = np.zeros((NumAt,16),int)
    for iAt in range(1,NumAt+1):
        AtOrbAll[iAt-1] = AtOrb[AtTypeInd[iAt-1]]
    # Trans uStt to iAt & iOrb
    Stt2AtOrb = np.zeros((NumStt,2),int)
    count = 0
    for iAt in range(1,NumAt+1):
        Ind = np.where(AtOrbAll[iAt-1]==1)[0]
        for Indi in Ind:
            Stt2AtOrb[count] = [iAt,Indi]
            count += 1
    # Trans iAt to uStt
    At2Stt = np.zeros((NumAt,2),int)
    for iAt in range(1,NumAt+1):
        Ind = np.where(Stt2AtOrb[:,0]==iAt)[0]
        i1 = np.min(Ind)
        i2 = np.max(Ind) + 1
        At2Stt[iAt-1] = [i1,i2]
    '''################################################'''
    
    '''###################   H(gk)   ##################'''
    # k, Random k-points
    NumKpt = 10
    KptLvRand = 2*np.random.random((NumKpt,3)) - 1
    RecLv = GetRecLv(Lv)
    KptRand = KptLvRand @ RecLv
    
    # gk, k-points acted by symmetry operations
    RecLvInv = np.linalg.inv(RecLv)
    SymKptRand = np.zeros((NumSym,NumKpt,3))
    SymKptLvRand = np.zeros((NumSym,NumKpt,3))
    for iSym in range(NumSym):
        SymKptRand[iSym] = KptRand @ SymXyz[iSym].T
        SymKptLvRand[iSym] = SymKptRand[iSym] @ RecLvInv
    
    # H(gk) - periodic gauge
    HmtKptFunc  = HamiltonianK(HmtLv, HmtMatLv)
    HmtSymKpt0   = np.zeros((NumSym,NumKpt,NumStt,NumStt),complex)
    for iSym in range(NumSym):
        for iKpt in range(NumKpt):
            HmtSymKpt0[iSym,iKpt] = HmtKptFunc.Hk(SymKptLvRand[iSym,iKpt])
    
    # H(gk) - symmetric gauge
    FctSymKpt = np.zeros((NumSym,NumKpt,NumStt,NumStt),complex)
    for iSym in range(NumSym):
        for iKpt in range(NumKpt):
            for iAt in range(1,NumAt+1):
                i1, i2 = At2Stt[iAt-1]
                FctSymKpt[iSym,iKpt,i1:i2,i1:i2] = np.eye(i2-i1) * np.exp(1j * pi2 * (SymKptLvRand[iSym,iKpt] @ AtLv[iAt-1]))
    HmtSymKpt = np.zeros((NumSym,NumKpt,NumStt,NumStt),complex)
    for iSym in range(NumSym):
        for iKpt in range(NumKpt):
            HmtSymKpt[iSym,iKpt] = FctSymKpt[iSym,iKpt].conjugate() @ HmtSymKpt0[iSym,iKpt] @ FctSymKpt[iSym,iKpt]
    '''################################################'''
    
    '''###################   D(g)   ##################'''
    # D^orb_ab = <a|gb> = <a|g|b> --- From HopRel
    # Write Symmetries acting on SPDF together
    SymSPDF = np.zeros((NumSym,16,16))
    for iOrb in range(4):
        SymSPDF[:,iOrb**2:(iOrb+1)**2,iOrb**2:(iOrb+1)**2] = SymOrb[iOrb]
    
    # Write Symmetries acting on input orbitals
    Ind  = [np.ix_(np.where(AtOrb[iAt])[0],np.where(AtOrb[iAt])[0]) for iAt in range(NumAtType)]
    SymAtOrb = [] # DgOrb
    for iAt in range(NumAtType):
        SymAtOrbi = []
        for iSym in range(NumSym):
            SymAtOrbi.append(SymSPDF[iSym][Ind[iAt]])
        for i in range(AtNum[iAt]):
            SymAtOrb.append(np.array(SymAtOrbi))
    
    # D(g)
    Dg = np.zeros((NumSym,NumStt,NumStt))
    for iSym in range(NumSym):
        for iAt in range(1,NumAt+1):
            OiAt = SymAt0[iSym,iAt-1,-1]
            i1, i2 = At2Stt[ iAt-1]
            j1, j2 = At2Stt[OiAt-1]
            # Dg[iSym,i1:i2,j1:j2] = SymAtOrb[iAt-1][iSym]
            Dg[iSym,j1:j2,i1:i2] = SymAtOrb[OiAt-1][iSym]
    '''################################################'''
    
    '''#############  D(g^−1) H(gk) D(g)  #############'''
    # D(g^−1) H(gk) D(g)
    DHD = np.zeros((NumSym,NumKpt,NumStt,NumStt),complex)
    for iSym in range(NumSym):
        for iKpt in range(NumKpt):
            DHD[iSym,iKpt] = Dg[iSym].T @ HmtSymKpt[iSym,iKpt] @ Dg[iSym]
    
    # D(g^−1) H(gk) D(g) == H(k)
    dH = np.zeros(NumSym)
    for iSym in range(NumSym):#[4]:#
        for iKpt in range(NumKpt):#[0]:#
            dH[iSym] += np.linalg.norm(DHD[iSym,iKpt] - HmtSymKpt[0,iKpt])
    print(np.linalg.norm(dH))
    '''################################################'''
    
def GetRecLv(Lv):
    a1,a2,a3 = Lv
    Vol = np.dot(a1,np.cross(a2,a3))
    b1 = pi2 * np.cross(a2,a3) / Vol
    b2 = pi2 * np.cross(a3,a1) / Vol
    b3 = pi2 * np.cross(a1,a2) / Vol
    RecLv = np.array([b1,b2,b3])
    return RecLv