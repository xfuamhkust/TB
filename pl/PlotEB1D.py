import matplotlib.pyplot as plt
import numpy as np


plt.rcParams['figure.dpi']= 300

def PlotEnergyBand1D(ParaEig1D,folderName=""):
    
    # Parameters
    Eig = ParaEig1D["Eig"]
    nk, NumStt = Eig.shape
    Kpt = np.linspace(-1,1,nk)
    
    # Plot energy bands of a strip
    fig = plt.figure()
    plt.plot(Kpt,Eig,"k")
    plt.xlabel("k (" + chr(960) + "/a)")
    plt.ylabel("Energy")
    plt.savefig(folderName+ "/Qc.svg")
        
    # Output
    ParaEig1dPlt = {"FigureEB1D": fig,
                  }
    
    return ParaEig1dPlt