import numpy as np
import sympy as smp
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


def HkSp2Np(HkMatSp,k1,k2,k3,tValsSpAll,k1Val,k2Val,k3Val,HopValIn):
    """

    :param HkMatSp: sympy matrix of Hk
    :param k1: momentum symbol
    :param k2: momentum symbol
    :param k3: momentum symbol
    :param tValsSpAll: a list of symbols of free hopping coefficients
    :param k1Val: momentum value
    :param k2Val: momentum value
    :param k3Val: momentum value
    :param HopValIn: a list of free hopping coefficients
    :return: numpy matrix of Hk
    """
    inSymbols=[k1,k2,k3]
    for ti in tValsSpAll:
        inSymbols.append(ti)
    HkNp=smp.lambdify(inSymbols,HkMatSp,"numpy")


    HopValInList=list(HopValIn)
    numToAdd=len(inSymbols)-len(HopValInList)-3
    for i in range(0,numToAdd):
        HopValInList.append(0)

    return HkNp(k1Val,k2Val,k3Val,*HopValInList)


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
    
    