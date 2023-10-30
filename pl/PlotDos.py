import matplotlib.pyplot as plt

plt.rcParams['figure.dpi']= 300

def PlotDOS(ParaDos,folderName=""):
    
    # Parameters
    EigVal = ParaDos["EigVal"]
    DosE   = ParaDos["DosE"]
    
    # Plot energy bands
    fig = plt.figure()
    plt.plot(EigVal,DosE)
    plt.xlabel("E")
    plt.ylabel("DOS")
    plt.savefig(folderName+ "/DOS.svg")
        
    # Output
    ParaDosPlt = {"FigureDOS": fig,
                  }
    
    return ParaDosPlt
