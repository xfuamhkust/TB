import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']= 300

def PlotEnergyBand(ParaKpt,ParaEig,Name=""):
    
    # Parameters
    Method = ParaKpt["Method"]
    EigVal = ParaEig["EigVal"]
    NumKpt, NumStt = EigVal.shape
    
    # Plot energy bands
    if Method == "Volume":
        pass
    elif Method == "Line":
        Kpt0Name = ParaKpt["Kpt0Name"]
        Kpt0R    = ParaKpt["Kpt0R"]
        KptR     = ParaKpt["KptR"]
        NumKpt0  = len(Kpt0R)
        fig = plt.figure()
        for iKpt0 in range(NumKpt0):
            plt.axvline(Kpt0R[iKpt0],c='grey',ls='--')
        for iStt in range(NumStt):
            plt.plot(KptR,EigVal[:,iStt],"k-")
        plt.xticks(Kpt0R,Kpt0Name)
        plt.ylabel("Energy")
        plt.xlim([KptR[0],KptR[-1]])
        plt.savefig("data/" + Name + "/EnergyBand.svg")
        
    # Output
    ParaEigPlt = {"FigureEB": fig,
                  }
    
    return ParaEigPlt