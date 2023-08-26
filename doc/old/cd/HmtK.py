import numpy as np

def GetHamiltonianK(ParaHmtR,ParaKpt):
    
    # Parameters
    HmtMatLv      = ParaHmtR["HmtMatLv"]
    HmtLv         = ParaHmtR["HmtLv"]
    KptLv         = ParaKpt["KptLv"]
    NumKpt = len(KptLv)
    NumStt = len(HmtMatLv[0])
    
    # Calculate Hamiltonian in k space
    HmtKptFunc = HamiltonianK(HmtLv, HmtMatLv)
    HmtKpt   = np.zeros((NumKpt,NumStt,NumStt),complex)
    for iKpt in range(NumKpt):
        HmtKpt[iKpt] = HmtKptFunc.Hk(KptLv[iKpt])
    
    # Output
    ParaHmtK = {"HmtKpt": HmtKpt,
                }
    
    return ParaHmtK

class HamiltonianK():
    
    def __init__(self, HmtLv, HmtMatLv, Lv=np.eye(3)):
        self.HmtLv = HmtLv
        self.NumLv, self.NumStt = HmtMatLv[:,:,0].shape
        self.HmtMatLv = HmtMatLv.reshape(self.NumLv, self.NumStt**2)
        self.Lvx_HmtMatLv = 1j * (self.HmtLv @ Lv[:,0])[:,None] * self.HmtMatLv # i * Rn * hk
        self.Lvy_HmtMatLv = 1j * (self.HmtLv @ Lv[:,1])[:,None] * self.HmtMatLv
        
    def Hk(self,Kpt):#Kpt is m = [m1,m2,m3], HmtLv is n = [n1,n2,n3]. k·Rn = 2 pi (n1m1+n2m2+n3m3)
        exp_kRn = np.exp(1j * 2 * np.pi * np.sum(Kpt * self.HmtLv,axis=1))
        Hk      = (exp_kRn @ self.HmtMatLv).reshape(self.NumStt, self.NumStt) # np.tensordot(a, b, axes=([0],[0]))
        return Hk
    
    def dHk(self,Kpt):
        exp_kRn = np.exp(1j * 2 * np.pi * np.sum(Kpt * self.HmtLv,axis=1))
        dH_dkx  = (exp_kRn @ self.Lvx_HmtMatLv).reshape(self.NumStt, self.NumStt)
        dH_dky  = (exp_kRn @ self.Lvy_HmtMatLv).reshape(self.NumStt, self.NumStt)
        return dH_dkx, dH_dky
    
    