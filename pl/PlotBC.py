import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['figure.dpi']= 300

def PlotBerryCurvature(ParaKpt,ParaBerK,folderName=""):
    
    # Parameters
    Method = ParaKpt["Method"]
    BerKpt = ParaBerK["BerKpt"]
    NumKpt, NumStt, _ = BerKpt.shape
    
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
            for jStt in range(iStt+1,NumStt):
                plt.plot(KptR,BerKpt[:,iStt,jStt],label = "%d"%(iStt+1) + "_" + "%d"%(jStt+1))
        plt.xticks(Kpt0R,Kpt0Name)
        plt.ylabel("Berry Curvature")
        plt.xlim([KptR[0],KptR[-1]])
        plt.legend()
        plt.savefig(folderName+ "/BerryCurvature.svg")
        
    # Output
    ParaBerPlt = {"FigureBC": fig,
                  }
    
    return ParaBerPlt
