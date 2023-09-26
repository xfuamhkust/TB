import sys
import numpy as np
import matplotlib.pyplot as plt
from cd.HmtK import HamiltonianK
pi2 = np.pi * 2

def CheckEnergySymmetry(ParaIn,ParaSym,ParaSymAt,ParaHmtR,Visual = True):
    
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
    np.random.seed(10)
    # Parameters
    Lv            = ParaIn["LatticeVector"]
    SymXyz        = ParaSym["SymXyz"]
    HmtMatLv      = ParaHmtR["HmtMatLv"]
    HmtLv         = ParaHmtR["HmtLv"]
    AtLv          = ParaIn["AtomSite"]
    SymAt0        = ParaSymAt["SymAt0"]
    AtOrb         = ParaIn["AtomOrbital"]
    AtTypeInd     = ParaIn["AtomTypeIndex"]
    NumSym        = len(SymXyz)
    NumStt        = len(HmtMatLv[0])
    NumAt         = len(SymAt0[0])
    AtomOrbNum    = np.sum(AtOrb,axis=1)
    AtOrbAll = np.zeros((NumAt,94),int)
    for iAt in range(1,NumAt+1):
        AtOrbAll[iAt-1] = AtOrb[AtTypeInd[iAt-1]]
    AtomOrbNumAll = np.sum(AtOrbAll,axis=1)
    AtomOrbNum = AtomOrbNumAll
    
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
    FctKpt = np.zeros((NumKpt,NumStt,NumStt),complex)
    for iKpt in range(NumKpt):
        for iAt in range(1,NumAt+1):
            i1 = np.sum(AtomOrbNum[: iAt-1])
            i2 = np.sum(AtomOrbNum[: iAt-0])
            FctKpt[iKpt,i1:i2,i1:i2] = np.eye(i2-i1) * np.exp(1j * pi2 * (KptLvRand[iKpt] @ AtLv[iAt-1]))
    HmtSymKpt = np.zeros((NumSym,NumKpt,NumStt,NumStt),complex)
    for iSym in range(NumSym):
        for iKpt in range(NumKpt):
            HmtSymKpt[iSym,iKpt] = FctKpt[iKpt].conjugate() @ HmtSymKpt0[iSym,iKpt] @ FctKpt[iKpt]
    
    # E(gk) - periodic gauge & symmetric gauge
    EigVal0 = np.zeros((NumSym,NumKpt,NumStt))
    EigVal  = np.zeros((NumSym,NumKpt,NumStt))
    for iSym in range(NumSym):
        for iKpt in range(NumKpt):
            EigVal0[iSym,iKpt],_ = np.linalg.eigh(HmtSymKpt0[iSym,iKpt])
            EigVal[iSym,iKpt],_  = np.linalg.eigh(HmtSymKpt[iSym,iKpt])
            
    # E(gk) == E(k)
    dE0 = np.linalg.norm(EigVal0 - EigVal0[0])
    dE  = np.linalg.norm(EigVal - EigVal[0])
    print(dE0)
    # print(dE)
    # print(np.linalg.norm(EigVal0-EigVal)) # E^P(gk) == E^S(gk)
    
    # Plot
    if Visual:
        plt.figure()
        marker0 = ['.',',','o','v','^','<','>','1','2','3','4','s','p','*','h','H','+','x','D','d','|','_','.',',']
        xs = np.arange(NumKpt)
        LW = np.linspace(5,1,NumSym)
        MS = np.linspace(10,1,NumSym)
        for iSym in range(NumSym):
            for iStt in range(NumStt):
                plt.plot(xs,EigVal0[iSym,:,iStt],linewidth=LW[iSym],markersize = MS[iSym], marker = marker0[iStt%len(marker0)])
        
def GetRecLv(Lv):
    a1,a2,a3 = Lv
    Vol = np.dot(a1,np.cross(a2,a3))
    b1 = pi2 * np.cross(a2,a3) / Vol
    b2 = pi2 * np.cross(a3,a1) / Vol
    b3 = pi2 * np.cross(a1,a2) / Vol
    RecLv = np.array([b1,b2,b3])
    return RecLv
