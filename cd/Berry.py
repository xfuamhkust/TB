import numpy as np
from cd.HmtK import HamiltonianK
e, kB, h_  = 1.6021766208e-19, 1.3806505e-23, 1.05457266e-34
np.seterr(divide='ignore', invalid='ignore')

def GetBerryK(ParaIn,ParaHmtR,ParaKpt):
    
    # Parameters
    Lv         = ParaIn["LatticeVector"]
    Dim        = ParaIn["Dimension"]
    HmtMatLv   = ParaHmtR["HmtMatLv"]
    HmtLv      = ParaHmtR["HmtLv"]
    KptLv      = ParaKpt["KptLv"]
    
    # Calculate Berry curvatures
    BerKptFunc = BerryK(Dim,Lv,HmtLv,HmtMatLv)
    BerKpt     = BerKptFunc.OmegaKpt(KptLv)
    
    # Output
    ParaBerK = {"BerKpt": BerKpt,
                }
    
    return ParaBerK

def GetAHNC(ParaIn,ParaHmtR,ParaBerIn):
    
    # Parameters
    Lv         = ParaIn["LatticeVector"]
    Dim        = ParaIn["Dimension"]
    HmtMatLv   = ParaHmtR["HmtMatLv"]
    HmtLv      = ParaHmtR["HmtLv"]
    nk0        = ParaBerIn["nk0"]
    Tem        = ParaBerIn["Tem"]
    EF         = ParaBerIn["EF"]
    
    
    # Calculate Berry curvatures
    BerKptFunc = BerryK(Dim,Lv,HmtLv,HmtMatLv)
    _, _, AHC, ANC = BerKptFunc.AHNC(nk0,EF,Tem)
    
    # Output
    ParaBerK = {"AHC": AHC,
                "ANC": ANC,
                }
    
    return ParaBerK

class BerryK():

    def __init__(self,Dim,Lv,HmtLv,HmtMatLv):
        self.Dim = Dim
        self.Vr  = Lv[0] @ np.cross(Lv[1],Lv[2])
        self.Hmt = HamiltonianK(HmtLv, HmtMatLv, Lv)
        
    def AHNC(self,nk0,EF,T):
        k0 = np.linspace(1/nk0-1,1-1/nk0,nk0)*0.5
        Kpt = np.array(np.meshgrid(k0,k0,k0 if self.Dim == 3 else np.array([0]))).reshape(3,nk0**self.Dim).T
        NumKpt, NumPara = len(Kpt), len(EF)
        Omegaf = np.zeros((NumKpt, NumPara))
        Omegag = np.zeros((NumKpt, NumPara))
        for iKpt in range(NumKpt):
            Omegaf[iKpt],Omegag[iKpt] = self.GetOmega(Kpt[iKpt],EF,T,Berry=0,Berryfg=1)
        AHC = -(e**2/h_/self.Vr) * np.mean(Omegaf,axis=0)
        ANC =  (e*kB/h_/self.Vr) * np.mean(Omegag,axis=0)
        return Omegaf, Omegag, AHC, ANC
    
    def OmegaKpt(self,Kpt):
        NumKpt, NumStt = len(Kpt), len(self.Hmt.Hk([0,0,0]))
        Omega = np.zeros((NumKpt, NumStt, NumStt))
        for iKpt in range(NumKpt):
            Omega[iKpt] = self.GetOmega(Kpt[iKpt],Berry=1,Berryfg=0)
        return Omega
    
    def GetOmega(self,Kpt,EF=0,T=0,Berry=1,Berryfg=0):
        H_k, [dH_dkx, dH_dky] = self.Hmt.Hk(Kpt), self.Hmt.dHk(Kpt)
        EigVal, EigVct = np.linalg.eigh(H_k)
        EigVct_H = EigVct.conjugate().T
        dH_dkx_psi, dH_dky_psi = EigVct_H @ dH_dkx @ EigVct, EigVct_H @ dH_dky @ EigVct
        dE = np.triu(1./((EigVal[:,None] - EigVal[None,:])**2 + np.eye(len(EigVal))), k=1) #1/(Em-En)^2
        IndInf = np.isinf(dE); IndNan = np.isnan(dE)
        dE[IndInf] = 0; dE[IndNan] = 0
        Omega = -2*np.imag(dH_dkx_psi * dH_dky_psi.T) * dE # Omega = i*[<psi_m|dH/dkx|psi_n><psi_n|dH/dky|psi_m> - (x<->y)]/(Em-En)^2
        if Berry == 1 and Berryfg == 0:
            return Omega
        Omegaf = np.zeros(len(T))
        Omegag = np.zeros(len(T))
        for iPara in range(len(T)): 
            E = EigVal-EF[iPara]
            f = np.zeros(len(E))
            g = np.zeros(len(E))
            f[E < 0], f[E == 0] = 1., 0.5 
            g[E == 0] = 0.6931471805599453
            if T[iPara]:
                E_ = (e/kB)*(E/T[iPara])
                Ind = np.where(abs(E_)< 100)
                f[Ind] = 1./(1.+np.exp(E_[Ind]))
                g[Ind] = E_[Ind]*f[Ind]+np.log(1+np.exp(-E_[Ind]))
                # g = E_*f+np.logaddexp(0.0,-E_)
            df = f[:,None] - f[None,:]
            dg = g[:,None] - g[None,:]
            Omegaf[iPara] = np.sum(np.triu(Omega*df, k=1))
            Omegag[iPara] = np.sum(np.triu(Omega*dg, k=1))
        if Berry == 0 and Berryfg == 1:
            return Omegaf, Omegag
        if Berry == 1 and Berryfg == 1:
            return Omega, Omegaf, Omegag