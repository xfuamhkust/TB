import matplotlib.pyplot as plt
import sympy as smp
import numpy as np
from cd.HmtK import HkSp2Np


plt.rcParams['figure.dpi']= 300

def PlotEnergyBand(ParaKpt,ParaEig,folderName=""):
    
    # Parameters
    Method = ParaKpt["Method"]
    EigVal = ParaEig["EigVal"]
    NumKpt, NumStt = EigVal.shape
    
    # Plot energy bands
    if Method == "Volume":
        pass
    elif Method == "Line":
        Kpt0Name = ParaKpt["Kpt0Name"]# names of the special points along the path in  BZ
        Kpt0R    = ParaKpt["Kpt0R"]#cumulative distances between special points  starting from KptPath0Xyz[0]
        KptR     = ParaKpt["KptR"]# cumulative distances along the path in BZ
        NumKpt0  = len(Kpt0R)# number of special points along the path
        fig = plt.figure()
        for iKpt0 in range(NumKpt0):
            plt.axvline(Kpt0R[iKpt0],c='grey',ls='--')
        for iStt in range(NumStt):
            plt.plot(KptR,EigVal[:,iStt],"k-")# TODO:  magnitude of eigenvalues
        plt.xticks(Kpt0R,Kpt0Name)
        plt.ylabel("Energy")
        plt.xlim([KptR[0],KptR[-1]])
        plt.savefig(folderName+ "/EnergyBand.svg")
        
    # Output
    ParaEigPlt = {"FigureEB": fig,
                  }
    
    return ParaEigPlt



def plotEigsFromHkSp(HkMatSp,k1,k2,k3,tValsSpAll,ParaKpt,HopValIn,ParaRel,folderName="./"):
    """

    :param HkMatSp: sympy matrix of Hk
    :param k1: momentum symbol
    :param k2: momentum symbol
    :param k3: momentum symbol
    :param tValsSpAll: a list of symbols of free hopping coefficients
    :param ParaKpt: contains momentum values
    :param HopValIn:  a list of free hopping coefficients
    :param ParaRel:
    :return: eigenvalues for each k
    """

    Method = ParaKpt["Method"]
    # Plot energy bands
    if Method == "Volume":
        pass
    elif Method == "Line":
        Kpt0Name = ParaKpt["Kpt0Name"]  # names of the special points along the path in  BZ
        Kpt0R = ParaKpt["Kpt0R"]  # cumsum of distances between special points  starting from KptPath0Xyz[0]
        KptR = ParaKpt["KptR"]  # cumulative distances along the path in BZ
        NumKpt0 = len(Kpt0R)  # number of special points along the path
        KptXyz = ParaKpt["KptXyz"]#interpolated points along the path under Cartesian basis in BZ



        NumKpt=len(KptXyz)
        NumStt,_=HkMatSp.shape
        # Calculate eigenvalues
        EigVal = np.zeros((NumKpt, NumStt))

        for iKpt in range(NumKpt):
            k1Val,k2Val,k3Val=KptXyz[iKpt]
            HkTmp = HkSp2Np(HkMatSp, k1, k2, k3, tValsSpAll, k1Val, k2Val, k3Val, HopValIn)
            EigVal[iKpt],_=np.linalg.eigh(HkTmp)
        fig = plt.figure()
        for iKpt0 in range(NumKpt0):
            plt.axvline(Kpt0R[iKpt0], c='grey', ls='--')

        for iStt in range(NumStt):
            plt.plot(KptR, EigVal[:, iStt], "b-")  # TODO:  magnitude of eigenvalues

        plt.xticks(Kpt0R, Kpt0Name)
        plt.ylabel("Energy")
        plt.xlim([KptR[0], KptR[-1]])
        plt.savefig(folderName + "/Sympy2EnergyBand.svg")

        plt.close()
        return fig







