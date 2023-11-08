import matplotlib.pyplot as plt

plt.rcParams['figure.dpi']= 300

def PlotQuanCond(ParaQc,folderName=""):
    
    # Parameters
    Ener = ParaQc["Ener"]
    TLR  = ParaQc["TLR"]
    
    # Plot energy bands
    fig = plt.figure()
    plt.plot(Ener,TLR)
    plt.xlabel("Energy")
    plt.ylabel("Quantum Conductance $(e^2/h)$")
    plt.savefig(folderName+ "/Qc.svg")
        
    # Output
    ParaQcPlt = {"FigureQc": fig,
                  }
    
    return ParaQcPlt
