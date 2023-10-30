import numpy as np
from cd.HmtK import HamiltonianK
pi = np.pi

def GetDos(ParaIn,ParaHmtR,ParaDosIn):
    
    # Parameters
    Lv         = ParaIn["LatticeVector"]
    Dim        = ParaIn["Dimension"]
    HmtMatLv   = ParaHmtR["HmtMatLv"]
    HmtLv      = ParaHmtR["HmtLv"]
    Method     = ParaDosIn["Method"]
    nE         = ParaDosIn["nE"]
    Emin       = ParaDosIn["Emin"]
    Emax       = ParaDosIn["Emax"]
    nk0        = ParaDosIn["nk0"]
    
    # Calculate DOS
    if Method == "Diag":
        eta        = ParaDosIn["eta"]
        EigVal, DosE = GetDosDiag(Dim,Lv,HmtLv,HmtMatLv,nE,Emin,Emax,eta,nk0)
    elif Method == "Inv":
        eta        = ParaDosIn["eta"]
        EigVal, DosE = GetDosInv( Dim,Lv,HmtLv,HmtMatLv,nE,Emin,Emax,eta,nk0)
    elif Method == "KPM":
        nsr        = ParaDosIn["nsr"]
        EigVal, DosE = GetDosKPM( Dim,Lv,HmtLv,HmtMatLv,nE,Emin,Emax,nsr,nk0)
    
    # Output
    ParaDos = {"EigVal": EigVal,
               "DosE": DosE,
               }
    
    return ParaDos

def GetDosDiag(Dim,Lv,HmtLv,HmtMatLv,nE,Emin,Emax,eta,nk0):
    # Parameters
    Vr  = Lv[0] @ np.cross(Lv[1],Lv[2])
    Hmt = HamiltonianK(HmtLv, HmtMatLv, Lv)
    Ein = np.linspace(Emin,Emax,nE)
    # k points
    k0 = np.linspace(1/nk0-1,1-1/nk0,nk0)*0.5
    Kpt = np.array(np.meshgrid(k0,k0,k0 if Dim == 3 else np.array([0]))).reshape(3,nk0**Dim).T
    # Energies
    NumKpt, NumStt = len(Kpt), len(Hmt.Hk([0,0,0]))
    EigVal = np.zeros((NumKpt,NumStt))
    for iKpt in range(NumKpt):
        EigVal[iKpt] = np.linalg.eigvalsh(Hmt.Hk(Kpt[iKpt]))
    # DOS
    DosE = np.zeros(nE)
    for iE in range(nE):
        DosE[iE] = np.sum(1 / ((EigVal - Ein[iE])**2 + eta**2))
    DosE *= (eta/pi/NumKpt/Vr)
    return Ein, DosE

def GetDosInv(Dim,Lv,HmtLv,HmtMatLv,nE,Emin,Emax,eta,nk0):
    # Parameters
    Vr  = Lv[0] @ np.cross(Lv[1],Lv[2])
    Hmt = HamiltonianK(HmtLv, HmtMatLv, Lv)
    Ein = np.linspace(Emin,Emax,nE)
    # k points
    k0 = np.linspace(1/nk0-1,1-1/nk0,nk0)*0.5
    Kpt = np.array(np.meshgrid(k0,k0,k0 if Dim == 3 else np.array([0]))).reshape(3,nk0**Dim).T
    # Traces: Tr[1/(E-H+i*eta)]
    NumKpt, NumStt = len(Kpt), len(Hmt.Hk([0,0,0]))
    Ein_ = np.zeros((nE,NumStt,NumStt),complex)
    for iE in range(nE):
        Ein_[iE] = (Ein[iE]+1j*eta)*np.eye(NumStt) 
    TrG = np.zeros((nE,NumKpt))
    for iKpt in range(NumKpt):
        HmtKpt = Hmt.Hk(Kpt[iKpt])
        for iE in range(nE):
            TrG[iE,iKpt] = np.imag(np.trace(np.linalg.inv(Ein_[iE]-HmtKpt)))
    # DOS
    DosE = np.sum(TrG,axis=1) * (-1/pi/NumKpt/Vr)
    return Ein, DosE

def GetDosKPM(Dim,Lv,HmtLv,HmtMatLv,nE,Emin,Emax,N,nk0):
    '''
        This is a function to calculate DOS by KPM.
        Refrences:
            https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.78.275
            https://docs.pybinding.site/en/stable/tutorial/kpm.html
            Peter's thesis
    '''
    # Parameters
    Vr  = Lv[0] @ np.cross(Lv[1],Lv[2])
    Hmt = HamiltonianK(HmtLv, HmtMatLv, Lv)
    Ein = np.linspace(Emin,Emax,nE)
    # k points
    k0 = np.linspace(1/nk0-1,1-1/nk0,nk0)*0.5
    Kpt = np.array(np.meshgrid(k0,k0,k0 if Dim == 3 else np.array([0]))).reshape(3,nk0**Dim).T
    # Scale energy
    E_ = np.linspace(-1,1,nE) 
    dE = (Emax-Emin)/2
    kH = 1/dE
    bH = -(Emax+Emin)/dE/2
    # g_m
    gm = np.zeros(N+1)
    for m in range(N+1):
        gm[m] = 1/(N+1)*((N-m+1) * np.cos(pi*m/(N+1)) + np.sin(pi*m/(N+1)) / np.tan(pi/(N+1)))
    # T_m(E)
    TmE = np.zeros((N+1,nE))
    TmE[0] = 1; TmE[1] = E_
    for m in range(2,N+1):
        TmE[m] = 2*E_*TmE[m-1] - TmE[m-2]
    # Tr[T_m(H_k)]
    NumKpt, NumStt = len(Kpt), len(Hmt.Hk([0,0,0]))
    TrTmHk = np.zeros((N+1,NumKpt))
    for iKpt in range(NumKpt):
        HmtKpt = Hmt.Hk(Kpt[iKpt]) 
        HmtKpt_ = kH * HmtKpt + bH # Scale Hamiltonian
        TmH1 = np.eye(NumStt)     # T_0(H_k)
        TmH2 = np.copy(HmtKpt_)   # T_1(H_k)
        TrTmHk[0,iKpt] = np.real(np.trace(TmH1))
        TrTmHk[1,iKpt] = np.real(np.trace(TmH2))
        for m in range(2,N):
            TmH3 = 2 * HmtKpt_ @ TmH2 - TmH1
            TrTmHk[m,iKpt] =  np.real(np.trace(TmH3))
            TmH1 = np.copy(TmH2)
            TmH2 = np.copy(TmH3)
    # amV
    amV = np.sum(TrTmHk,axis=1) * (1/NumKpt/Vr)
    # DOS
    DosE = np.zeros(nE)
    for iE in range(nE):
        DosE[iE] = (1/dE/pi/np.sqrt(1-E_[iE]**2)) *(gm[0]*amV[0] + 2 * np.sum(gm[1:]*amV[1:]*TmE[1:,iE]))
    return Ein, DosE

