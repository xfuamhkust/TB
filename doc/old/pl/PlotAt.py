import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']= 300

def PlotAtoms(ParaNbr,Name=""):
    
    # Parameters
    LvAt    = ParaNbr["LvAt"]
    LvAtXyz = ParaNbr["LvAtXyz"]
    
    # Classify different atoms
    LvAtType = LvAt[:,-1]
    NumType = np.max(LvAtType)
    IndType = np.array([np.where(LvAtType==i)[0] for i in range(1,NumType+1)])
    
    # Plot atoms
    fig = plt.figure()
    x = LvAtXyz[:,0]; y = LvAtXyz[:,1]; z = LvAtXyz[:,2]
    if np.linalg.norm(z):
        # 3D
        ax = fig.add_subplot(111, projection='3d')
        for i in range(NumType):
            ax.scatter(x[IndType[i]],y[IndType[i]],z[IndType[i]], label='122')
        xyzmax = np.max(abs(LvAtXyz))*1.1
        ax.set_xlim([-xyzmax,xyzmax])
        ax.set_ylim([-xyzmax,xyzmax])
        ax.set_zlim([-xyzmax,xyzmax])
    else:
        # 2D
        for i in range(NumType):
            plt.plot(x[IndType[i]],y[IndType[i]],".")
        plt.axis("equal")
    plt.show()
    
    # Output
    ParaEigPlt = {"FigureAt": fig,
                  }
    
    return ParaEigPlt